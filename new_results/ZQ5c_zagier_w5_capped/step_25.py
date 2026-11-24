import os
cwd = os.getcwd()
# Rewrite w3lib_c.py to save raw relation vectors for later expansion/verification
script = r'''import sys, json, time, itertools, math
from fractions import Fraction as Ff
from collections import defaultdict
import mpmath as mp
import flint
sys.path.insert(0,'.')
import w5blib as W
from w5blib import QD

def w3_run(D, F, cands, P):
    """Returns dict with ALL LLL relations (a-vectors), combos, kernel data."""
    d=len(F['letters']); N=len(cands); names=F['names']
    t0=time.time()
    sym_dicts=[W.dense_to_dict(W.symbol_liw_dense(c['vx'],c['v1m'],d,3)) for c in cands]
    K, orbits = W.symbol_kernel(sym_dicts, d, 3, N)
    mp.mp.dps=P+40
    sD=mp.sqrt(D)
    lg=[mp.log(abs(g.emb(P+40))) for g in F['letters']]
    li3v=[W.li3(c['x'].emb(P+40),P+20) for c in cands]
    z3=mp.zeta(3); pi2=mp.pi**2
    L3=W.Lchi_hurwitz(D,3,P)
    combos=[]
    comps=[]
    for row in K:
        v=sum(mp.mpf(c)*li3v[i] for i,c in enumerate(row) if c)
        stot=W.combo_symbol(row, sym_dicts, N)
        comp=W.companion(stot, d, 3, orbits)
        comps.append({str(ms): str(cf) for ms,cf in comp.items()})
        for ms, cf in comp.items():
            pv=mp.mpf(1)
            for ii in ms: pv*=lg[ii]
            v = v - mp.mpf(cf.numerator)/cf.denominator*pv
        combos.append(v)
    keep=[j for j,v in enumerate(combos) if abs(v)>mp.mpf(10)**(-(P-60))]
    tv = L3/sD
    vals_list=[tv, z3]+[combos[j] for j in keep]+[pi2*l for l in lg]
    rows,dt=W.lll_rows(vals_list, P)
    rels=[]
    for row in rows:
        a=row[:-1]
        if not any(a) or a[0]==0: continue
        mp.mp.dps=P+40
        res=sum(mp.mpf(x)*v for x,v in zip(a,vals_list))
        if abs(res)<mp.mpf(10)**(-(P-40)):
            rels.append([int(x) for x in a])
    print(f"  D={D} P={P}: kernel={len(K)} combos_kept={len(keep)} LLLdim={len(vals_list)} ({dt:.0f}s) rels={len(rels)} ({time.time()-t0:.0f}s)", flush=True)
    return dict(rels=rels, K=K, keep=keep, comps=comps, N=N, d=d, names=names,
                kernel=K)

def expand_a(D, F, cands, run, a):
    """Expand one relation a-vector to explicit exact terms. a_T*T = a_z3*z3 + li3 + prod + pi2."""
    K=run['K']; keep=run['keep']; comps=run['comps']; N=run['N']; d=run['d']
    a_T = int(a[0]); a_z3 = int(a[1])
    a_combos = [int(x) for x in a[2:2+len(keep)]]
    a_pi2 = [int(x) for x in a[2+len(keep):2+len(keep)+d]]
    li3c = defaultdict(Ff)
    for aj, kj in zip(a_combos, keep):
        if aj:
            row=K[kj]
            for i in range(N):
                if row[i]: li3c[i] += aj*row[i]
    prodc = defaultdict(Ff)
    for aj, kj in zip(a_combos, keep):
        if not aj: continue
        for ms_s, cf_s in comps[kj].items():
            prodc[eval(ms_s)] += -aj*Ff(cf_s)
    # a_T*T + a_z3*z3 + sum li3c*li3 + sum prodc*prod + sum a_pi2*pi2log = 0
    out_li3 = {i: -li3c[i] for i in li3c if li3c[i]}
    out_prod = {ms: -prodc[ms] for ms in prodc if prodc[ms]}
    out_pi2 = {i: -a_pi2[i] for i in range(d) if a_pi2[i]}
    return dict(a_T=a_T, a_z3=Ff(-a_z3), li3=out_li3, prod=out_prod, pi2=out_pi2)

def explicit_resid(D, F, cands, e, P):
    mp.mp.dps=P+40
    sD=mp.sqrt(D)
    lg=[mp.log(abs(g.emb(P+40))) for g in F['letters']]
    L3=W.Lchi_hurwitz(D,3,P)
    tot = Ff(e['a_T'])*(L3/sD) - e['a_z3']*mp.zeta(3)
    for i,cf in e['li3'].items():
        tot -= cf*W.li3(cands[i]['x'].emb(P+40), P+20)
    for ms,cf in e['prod'].items():
        pv=mp.mpf(1)
        for ii in ms: pv*=lg[ii]
        tot -= cf*pv
    for ii,cf in e['pi2'].items():
        tot -= cf*mp.pi**2*lg[ii]
    return tot

def symbol_check(cands, e, d):
    T=defaultdict(Ff)
    for i,cf in e['li3'].items():
        s=W.symbol_liw_dense(cands[i]['vx'],cands[i]['v1m'],d,3)
        for idx,v in enumerate(s):
            if v: T[idx]+=cf*v
    for ms,cf in e['prod'].items():
        v=W.prod_symbol_dense(list(ms),d,3)
        for i2,vv in enumerate(v):
            if vv: T[i2]+=cf*vv
    return {k:v for k,v in T.items() if v}
'''
open('w3lib_c.py','w').write(script)
print("rewritten")
