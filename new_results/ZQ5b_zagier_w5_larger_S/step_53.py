import os, json
from fractions import Fraction as Ff
os.chdir('<project>/research/08_theory/ZQ5b_zagier_w5_larger_S')
step = r'''
import sys, json, time
from fractions import Fraction as Ff
from collections import defaultdict
import mpmath as mp
sys.path.insert(0,'.')
import w5blib as W
from w5blib import QD
ZQ='../ZQ5_zagier_weight5'
J = json.load(open(f'{ZQ}/out/d17/candidates.json'))
cands=[dict(x=QD(Ff(c['x'][0]),Ff(c['x'][1]),17), vx=c['vx'], v1m=c['v1m']) for c in J['cands']]
names=J['names']; d=len(names); N=len(cands)
S = json.load(open('d17_min_setup.json'))
M=S['M']; prod_list=[(p[0],tuple(p[1]),p[2]) for p in S['prod_list']]
idx_of={n:k for k,n in enumerate(names)}
def load_c(fn):
    return [Ff(x) for x in json.load(open(fn))['c']]
forms = {
  'planner (94 Li5)': [Ff(x) for x in S['c_rel']],
  'maxnum-min (103 Li5)': load_c('d17_min_maxnum.json'),
  'term-min (74 Li5)': load_c('d17_min_cand.json'),
}
def stats(c):
    li5n=sum(1 for i in range(N) if c[1+i])
    mx=max(abs(x.numerator) for x in c)
    nt=sum(1 for x in c if x)
    import math
    lc=1
    for x in c:
        if x: lc=lc*x.denominator//math.gcd(lc,x.denominator)
    return li5n, mx, nt, lc
print("=== three forms (Li5 terms / max numerator / total terms / lcm den) ===")
for nm,c in forms.items():
    li5n,mx,nt,lc=stats(c)
    print(f"  {nm:24s}: Li5 {li5n:3d}  maxnum {mx:7d}  terms {nt:3d}  lcm {lc}")
# finalise: maxnum form
c = forms['maxnum-min (103 Li5)']
# exact symbol
T=defaultdict(Ff)
for i in range(N):
    if c[1+i]:
        s=W.symbol_liw_dense(cands[i]['vx'],cands[i]['v1m'],d,5)
        for idx,v in enumerate(s):
            if v: T[idx]+=c[1+i]*v
for pi in range(M):
    kind,ms,pnm=prod_list[pi]
    if kind!="log5" or not c[1+N+pi]: continue
    v=W.prod_symbol_dense(list(ms),d,5)
    for i2,vv in enumerate(v):
        if vv: T[i2]+=c[1+N+pi]*vv
Tz={k:v for k,v in T.items() if v}
print("exact full symbol of maxnum form:", "ZERO" if not Tz else f"NONZERO {len(Tz)}")
# cobracket
dic={W.idx_to_tup(k,d,5):v for k,v in T.items()}
for mode in ("full","cut"):
    r=W.cobracket32(dic,d,mode)
    print(f"delta_{{3,2}} [{mode}]: {len(r)} nonzero")
# write d17_final.json
li5out={}
for i in range(N):
    if c[1+i]:
        x=cands[i]['x']
        li5out[f"{x.a}+{x.b}*sqrt17"]=str(c[1+i])
prodout={prod_list[pi][2]: str(c[1+N+pi]) for pi in range(M) if c[1+N+pi]}
json.dump(dict(D=17, T_coef=-1419857, z5_coef=str(c[0]), li5=li5out, products=prodout,
               letters=J['disp'],
               note="minimum max-numerator form over the 43-dim exact-zero lattice; NUMERICAL (LLL), not a proof"),
          open('d17_final.json','w'), indent=1)
print("wrote d17_final.json: Li5 terms", len(li5out), " products", len(prodout), " z5", c[0])
# fresh residual at 1300 digits
P=1300
mp.mp.dps=P+30
s=mp.sqrt(17)
lg={n: mp.log(abs(g.emb(P+30))) for n,g in zip(names, W.make_field(17,split=[2],inert=[])['letters'])}
Lchi=W.Lchi_hurwitz(17,5,P)
tot = Ff(-1419857)*Lchi/s + c[0]*mp.zeta(5)
t0=time.time()
for i in range(N):
    if c[1+i]:
        tot += c[1+i]*W.li5(cands[i]['x'].emb(P+30), P+10)
for pi in range(M):
    if c[1+N+pi]:
        kind,ms,pnm=prod_list[pi]
        v=mp.mpf(c[1+N+pi].numerator)/c[1+N+pi].denominator
        if kind=="log5":
            for ii in ms: v*=lg[names[ii]]
        elif kind=="pi2log3":
            v*=mp.pi**2
            for ii in ms: v*=lg[names[ii]]
        else:
            v*=mp.pi**4*lg[names[ms[0]]]
        tot += v
print(f"fresh residual of d17_final @1300d: {mp.nstr(tot,4)}  (L chi17,5 = {mp.nstr(Lchi,26)})  ({time.time()-t0:.0f}s)")
'''
open('step_05_d17_final.py','w').write(step)
import subprocess
r = subprocess.run(['nice','-n','10','<project>/.venv/bin/python','step_05_d17_final.py'], capture_output=True, text=True, cwd='.', timeout=1700)
print(r.stdout); print("ERR:", r.stderr[-1200:])
