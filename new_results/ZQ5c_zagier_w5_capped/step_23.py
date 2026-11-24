import os
cwd = os.getcwd()
script = r'''import sys, json, time, itertools, math
from fractions import Fraction as Ff
from collections import defaultdict
import mpmath as mp
import flint
sys.path.insert(0,'.')
import w5blib as W
from w5blib import QD

# sanity: li3 vs mpmath
mp.mp.dps=120
for xs in ('0.3','0.6','0.9','-0.7'):
    x=mp.mpf(xs)
    assert abs(W.li3(x,100)-mp.polylog(3,x))<mp.mpf('1e-95'), xs
print("li3 sanity ok", flush=True)

def w3_rel(D, F, cands, P):
    """Run LLL for T3 = L(chi_D,3)/sqrt D over [T, z3, combos K_j-P_j, pi2*log_i].
    Returns (a, namelist, K, keep, sym_dicts, orbits) with a the best (smallest maxcoef) relation."""
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
    for row in K:
        v=sum(mp.mpf(c)*li3v[i] for i,c in enumerate(row) if c)
        stot=W.combo_symbol(row, sym_dicts, N)
        comp=W.companion(stot, d, 3, orbits)
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
            rels.append((max(abs(x) for x in a), a, res))
    rels.sort(key=lambda t:t[0])
    print(f"  D={D} P={P}: N={N} d={d} kernel={len(K)} combos_kept={len(keep)} LLLdim={len(vals_list)} ({dt:.0f}s) rels={len(rels)} ({time.time()-t0:.0f}s)", flush=True)
    for mx,a,res in rels[:2]:
        print(f"    maxcoef {mx}, T-coef {a[0]}, z3-coef {a[1]}, resid {mp.nstr(res,3)}, nterms {sum(1 for x in a if x)}", flush=True)
    best = rels[0] if rels else None
    return dict(best=best, K=K, keep=keep, sym_dicts=sym_dicts, orbits=orbits,
                vals=dict(L3=L3, tv=tv), N=N, d=d, names=names)

def expand(D, F, cands, run):
    """Expand best relation a to explicit li3 coefficients and product coefficients (exact rationals)."""
    mx, a, res = run['best']
    K=run['K']; keep=run['keep']; N=run['N']; d=run['d']; names=run['names']
    sym_dicts=run['sym_dicts']; orbits=run['orbits']
    a_T = int(a[0]); a_z3 = int(a[1])
    a_combos = [int(x) for x in a[2:2+len(keep)]]
    a_pi2 = [int(x) for x in a[2+len(keep):2+len(keep)+d]]
    # li3 coefficient of x_i = sum_j a_combos_j * K[keep_j][i]
    li3c = [Ff(0)]*N
    for aj, kj in zip(a_combos, keep):
        if aj:
            row=K[kj]
            for i in range(N):
                if row[i]: li3c[i] += aj*row[i]
    # products: from companions sum_j a_combos_j * P_j  (with MINUS sign: combo = K - P, so a_j*(K-P) -> product coef = -a_j*comp)
    prodc = defaultdict(Ff)
    for aj, kj in zip(a_combos, keep):
        if not aj: continue
        stot=W.combo_symbol(K[kj], sym_dicts, N)
        comp=W.companion(stot, d, 3, orbits)
        for ms, cf in comp.items():
            prodc[ms] += -aj*cf
    # pi2*log_i: coefficient a_pi2[i] on pi2*log_i
    # Full relation: a_T*T + a_z3*z3 + sum a_j*(K_j - P_j) + sum a_pi2_i*pi2*log_i = 0
    # => a_T*T = -a_z3*z3 - sum_i li3c_i*li3(x_i) - sum_ms prodc_ms*prod log + ... wait sign
    # sum_j a_j*(K_j-P_j) = sum_i li3c_i li3(x_i) + sum_ms prodc_ms prodms(log)
    # So a_T*T = -a_z3*z3 - li3c - prodc - a_pi2*pi2log
    out_li3 = {i: -li3c[i] for i in range(N) if li3c[i]}
    out_prod = {ms: -prodc[ms] for ms in prodc if prodc[ms]}
    out_pi2 = {i: -a_pi2[i] for i in range(d) if a_pi2[i]}
    return dict(a_T=a_T, a_z3=-a_z3, li3=out_li3, prod=out_prod, pi2=out_pi2,
                maxcoef=mx, resid=res)
'''
open('w3lib_c.py','w').write(script)
print("w3lib_c.py written", len(script))
