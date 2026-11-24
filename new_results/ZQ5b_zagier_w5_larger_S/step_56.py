import os, json
from fractions import Fraction as Ff
os.chdir('<project>/research/08_theory/ZQ5b_zagier_w5_larger_S')
step = r'''
import sys, json, time
from fractions import Fraction as Ff
import mpmath as mp
sys.path.insert(0,'.')
import w5blib as W
from w5blib import QD
exec(open('step_06_w3_d5.py').read().split('# D=5 control')[0].split('def w3_pipeline')[0])  # li3 sanity only
def w3_pipeline(D, F, cands, P=600, targets=('T3over','L')):
    d=len(F['letters']); N=len(cands); names=F['names']
    t0=time.time()
    sym_dicts=[W.dense_to_dict(W.symbol_liw_dense(c['vx'],c['v1m'],d,3)) for c in cands]
    K, orbits = W.symbol_kernel(sym_dicts, d, 3, N)
    print(f"  D={D}: N={N} d={d} kernel dim={len(K)} ({time.time()-t0:.0f}s)", flush=True)
    mp.mp.dps=P+40
    sD=mp.sqrt(D)
    lg=[mp.log(abs(g.emb(P+40))) for g in F['letters']]
    li3v=[W.li3(c['x'].emb(P+40),P+20) for c in cands]
    z3=mp.zeta(3); pi2=mp.pi**2
    L3=W.Lchi_hurwitz(D,3,P)
    print(f"  L(chi{D},3) = {mp.nstr(L3,30)}")
    print(f"  L/sqrt{D} = {mp.nstr(L3/sD,30)}")
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
    print(f"  {len(keep)} nonzero combos, {len(combos)-len(keep)} exact-zero ({time.time()-t0:.0f}s)", flush=True)
    tmap={'T3over': L3/sD, 'L': L3}
    for tname in targets:
        tv=tmap[tname]
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
        print(f"  D={D} target {tname}: LLL dim={len(vals_list)} {dt:.0f}s, relations={len(rels)}", flush=True)
        for mx,a,res in rels[:2]:
            print(f"    maxcoef {mx}, T-coef {a[0]}, residual {mp.nstr(res,3)}, nterms {sum(1 for x in a if x)}")
F13 = W.make_field(13, split=[3], inert=[2])
print("D=13 letters:", F13['names'])
c13 = W.enum_exceptional_fast(F13, B=6)
c13 = [c for c in c13 if c['h']<=3]
print(f"D=13: using {len(c13)} exceptional (h<=3, B=6)")
r13 = w3_pipeline(13, F13, c13, P=500)
'''
open('step_07_w3_d13.py','w').write(step)
import subprocess
r = subprocess.run(['nice','-n','10','<project>/.venv/bin/python','step_07_w3_d13.py'], capture_output=True, text=True, cwd='.', timeout=1700)
print(r.stdout); print("ERR:", r.stderr[-1500:])
