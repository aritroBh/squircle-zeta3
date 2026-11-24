import os
cwd = os.getcwd()
lib = r'''# w5verify.py -- expand a weight-5 LLL relation to explicit terms + verify
import sys, json, math, itertools
from fractions import Fraction as Ff
from collections import defaultdict
import mpmath as mp
sys.path.insert(0,'.')
import w5blib as W
from w5blib import QD

def expand_a(D, cands, K, comps, keepc, atom_names, a, names):
    """a over [T(D^5), z5, K_j(j in keepc), product atoms...] ->
    explicit: a_T*T = a_z5*z5 + sum li5coef*Li5 + sum prodcoef*products.
    products named 'lX*lY*...' and 'pi2*...', 'pi4*...', 'z3*...', 'T3*...'.
    Returns dict with li5={i:coef}, prod={name:coef}, a_T, a_z5, and aux atoms (z3/T3/pi coefs)."""
    d=len(names); nK=len(keepc)
    a_T=int(a[0]); a_z5=int(a[1])
    a_combos=[int(x) for x in a[2:2+nK]]
    a_atoms=[int(x) for x in a[2+nK:]]
    atom_prod_names = atom_names[2+nK:]
    # li5 coefficients: sum_j a_combos_j * K[keepc_j][i]
    li5c=defaultdict(Ff)
    for aj,kj in zip(a_combos, keepc):
        if aj:
            row=K[kj]
            for i in range(len(cands)):
                if row[i]: li5c[i]+=aj*row[i]
    # products: companions (from combos) with MINUS sign, plus explicit product atoms
    prodc=defaultdict(Ff)  # key: tuple of factors like ('lE13','lE13',...,'pi2')
    for aj,kj in zip(a_combos, keepc):
        if not aj: continue
        for ms_s, cf_s in comps[kj].items():
            ms=eval(ms_s); cf=Ff(cf_s)
            key=tuple(sorted(names[ii] for ii in ms))
            prodc[key] += -aj*cf
    # explicit product atoms
    aux=defaultdict(Ff)
    for aval, pname in zip(a_atoms, atom_prod_names):
        if not aval: continue
        parts=tuple(pname.split('*'))
        prodc[parts] += aval
    # Full relation: a_T*T + a_z5*z5 + sum li5c*Li5 + sum prodc*prod = 0  =>  a_T*T = -a_z5*z5 - li5 - prodc
    out_li5={i:-li5c[i] for i in li5c if li5c[i]}
    out_prod={k:-prodc[k] for k in prodc if prodc[k]}
    return dict(a_T=a_T, a_z5=Ff(-a_z5), li5=out_li5, prod=out_prod)

def prod_value(parts, lg_by_name, pi2, z3, T3):
    """parts: tuple of factor names. Returns mpf value."""
    v=mp.mpf(1)
    for p in parts:
        if p=='pi2': v*=pi2
        elif p=='pi4': v*=pi2**2
        elif p=='z3': v*=z3
        elif p=='T3': v*=T3
        else: v*=lg_by_name[p]
    return v

def explicit_resid(D, F, cands, e, P):
    """Compute a_T*T - a_z5*z5 - sum li5 - sum prod at precision P. e in RHS-form."""
    mp.mp.dps=P+40
    sD=mp.sqrt(D)
    names=F['names']
    lg={n: mp.log(abs(g.emb(P+40))) for n,g in zip(names, F['letters'])}
    pi2=mp.pi**2; z3=mp.zeta(3)
    L5=W.Lchi_hurwitz(D,5,P); T3=W.Lchi_hurwitz(D,3,P)/sD
    tot = Ff(e['a_T'])*(L5/sD) - e['a_z5']*mp.zeta(5)
    for i,cf in e['li5'].items():
        tot -= cf*W.li5(cands[i]['x'].emb(P+40), P+20)
    for parts,cf in e['prod'].items():
        tot -= cf*prod_value(parts, lg, pi2, z3, T3)
    return tot

def symbol_check5(cands, e, names):
    """Exact full symbol of RHS li5+log^5 products (over S-unit letters). pi2/z3/T3 excluded.
    Returns nonzero dict (should be empty)."""
    d=len(names)
    nidx={n:i for i,n in enumerate(names)}
    T=defaultdict(Ff)
    for i,cf in e['li5'].items():
        s=W.symbol_liw_dense(cands[i]['vx'],cands[i]['v1m'],d,5)
        for idx,v in enumerate(s):
            if v: T[idx]+=cf*v
    for parts,cf in e['prod'].items():
        if any(p in ('pi2','pi4','z3','T3') for p in parts): continue
        ms=[nidx[p] for p in parts]
        v=W.prod_symbol_dense(ms,d,5)
        for i2,vv in enumerate(v):
            if vv: T[i2]+=cf*vv
    return {k:v for k,v in T.items() if v}

def cobracket_check5(cands, e, names):
    """(3,2) cobracket of the full symbol dict; returns (full_count, cut_count)."""
    d=len(names)
    nidx={n:i for i,n in enumerate(names)}
    T=defaultdict(Ff)
    for i,cf in e['li5'].items():
        s=W.symbol_liw_dense(cands[i]['vx'],cands[i]['v1m'],d,5)
        for idx,v in enumerate(s):
            if v: T[idx]+=cf*v
    for parts,cf in e['prod'].items():
        if any(p in ('pi2','pi4','z3','T3') for p in parts): continue
        ms=[nidx[p] for p in parts]
        v=W.prod_symbol_dense(ms,d,5)
        for i2,vv in enumerate(v):
            if vv: T[i2]+=cf*vv
    dic={W.idx_to_tup(k,d,5):v for k,v in T.items()}
    rf=W.cobracket32(dic,d,'full')
    rc=W.cobracket32(dic,d,'cut')
    return len(rf), len(rc)

def factorization_check(cands, e, F):
    """Verify for every used argument: x (up to sign) and 1-x are S-units with the stored coords."""
    bad=0; nchk=0
    for i in e['li5']:
        c=cands[i]; nchk+=1
        sgn=1 if c['x'].emb(60)>0 else -1
        xabs=c['x'] if sgn==1 else -c['x']
        x=QD(1,0,F['D'])
        for g,ee in zip(F['letters'], c['vx']):
            if ee: x=x*(g**int(ee))
        okx=(x==xabs)
        f1=W.factor_sunit(1-c['x'], F)
        ok1=(f1==c['v1m'])
        if not (okx and ok1): bad+=1
    return nchk, bad
'''
open('w5verify.py','w').write(lib)
print("w5verify.py written")
import time
time.sleep(60)
print("=== D13 vals log after wait ===")
print(open('step_19_bg.log').read())
