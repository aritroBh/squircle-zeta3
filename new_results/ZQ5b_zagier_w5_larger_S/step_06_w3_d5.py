
import sys, json, time, itertools
from fractions import Fraction as Ff
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
print("li3 sanity ok")

def w3_pipeline(D, F, cands, P=600, targets=('T3over','T3times','L')):
    d=len(F['letters']); N=len(cands); names=F['names']
    t0=time.time()
    sym_dicts=[W.dense_to_dict(W.symbol_liw_dense(c['vx'],c['v1m'],d,3)) for c in cands]
    K, orbits = W.symbol_kernel(sym_dicts, d, 3, N)
    print(f"  D={D}: N={N} d={d} kernel dim={len(K)} ({time.time()-t0:.0f}s)", flush=True)
    # companions: for each K row, symbol combo, companion in log^3 products
    mp.mp.dps=P+40
    sD=mp.sqrt(D)
    lg=[mp.log(abs(g.emb(P+40))) for g in F['letters']]
    li3v=[W.li3(c['x'].emb(P+40),P+20) for c in cands]
    z3=mp.zeta(3); pi2=mp.pi**2
    L3=W.Lchi_hurwitz(D,3,P)
    # build combo values K_j - P_j
    combos=[]
    for row in K:
        v=sum(mp.mpf(c)*li3v[i] for i,c in enumerate(row) if c)
        stot=W.combo_symbol(row, sym_dicts, N)
        comp=W.companion(stot, d, 3, orbits)  # ms -> coef
        for ms, cf in comp.items():
            pv=mp.mpf(1)
            for ii in ms: pv*=lg[ii]
            v = v - mp.mpf(cf.numerator)/cf.denominator*pv
        combos.append(v)
    # drop exact-zero combos
    keep=[j for j,v in enumerate(combos) if abs(v)>mp.mpf(10)**(-(P-60))]
    dropped=[j for j,v in enumerate(combos) if abs(v)<=mp.mpf(10)**(-(P-60))]
    print(f"  D={D}: {len(keep)} nonzero combos, {len(dropped)} exact-zero ({time.time()-t0:.0f}s)")
    tmap={'T3over': L3/sD, 'T3times': L3*sD, 'L': L3}
    found={}
    for tname in targets:
        tv = tmap[tname]
        vals_list=[tv, z3]+[combos[j] for j in keep]+[pi2*l for l in lg]
        nm_list=[tname,'z3']+[f'K{keep[j]}' for j in range(len(keep))]+[f'pi2*{names[i]}' for i in range(d)]
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
        found[tname]=rels
        print(f"  D={D} target {tname}: LLL dim={len(vals_list)} {dt:.0f}s, relations={len(rels)}", flush=True)
        for mx,a,res in rels[:2]:
            print(f"    maxcoef {mx}, T-coef {a[0]}, residual {mp.nstr(res,3)}, nterms {sum(1 for x in a if x)}")
    return dict(K=K, keep=keep, found=found, L3=L3)

# D=5 control
F5 = W.make_field(5, split=[], inert=[2])
print("D=5 letters:", F5['names'])
c5 = W.enum_exceptional_fast(F5, B=8)
print(f"D=5: {len(c5)} exceptional")
r5 = w3_pipeline(5, F5, c5, P=500)
