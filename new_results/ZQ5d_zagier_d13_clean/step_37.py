script = r'''# verify_relation.py -- independent verification of a d13 relation JSON (d17_final format):
# residuals at 1300 & 2200 digits (planner_li5 + Hurwitz L), EXACT symbol, exact S-unit factorisations, stats.
import os, sys, json, math, time, itertools
os.nice(10)
from fractions import Fraction as Ff
from collections import defaultdict
import mpmath as mp
sys.path.insert(0,'.')
from planner_li5 import li5 as pli5
path = sys.argv[1]
rel = json.load(open(path))
J = json.load(open('out/d13/candidates_cap.json'))
GEN = {'lE13': (Ff(3,2),Ff(1,2)), 'ls13': (Ff(0),Ff(1)), 'l2': (Ff(2),Ff(0)), 'lp3': (Ff(5,2),Ff(1,2)), 'lp3c': (Ff(5,2),Ff(-1,2))}
assert rel['letters'] == {n:f'{a}+{b}*sqrt13' for n,(a,b) in GEN.items()}, 'letters mismatch'
T_coef=Ff(rel['T_coef']); z5_coef=Ff(rel['z5_coef'])
def parse(s):
    s=s.replace('+-','-'); i=s.rfind('*sqrt13'); head=s[:i]
    k=max(head.rfind('+'), head.rfind('-'))
    return Ff(head[:k]), Ff(head[k:])
li5=[(parse(k), Ff(v)) for k,v in rel['li5'].items()]
prods=[(tuple(k.split('*')), Ff(v)) for k,v in rel['products'].items()]
LET=list(GEN)
# match args to pool (exact) for coords
pool={(Ff(c['x'][0]),Ff(c['x'][1])):(c['vx'],c['v1m']) for c in J['cands']}
# stats
coefs=[c for _,c in li5]+[c for _,c in prods]+[z5_coef]
lm=1
for c in coefs: lm=math.lcm(lm,c.denominator)
print(f'{path}: T_coef={T_coef} (={T_coef/13**5}*13^5), z5={z5_coef}, li5 {len(li5)}, prods {len(prods)}, '
      f'max|num| {max(abs(c.numerator) for c in coefs)}, lcm {lm}', flush=True)
def chi13(a):
    a%=13; return 0 if a==0 else (1 if pow(a,6,13)==1 else -1)
for dps in (1300, 2200):
    t0=time.time(); mp.mp.dps=dps+30
    s=mp.sqrt(13)
    val=lambda ab: mp.mpf(ab[0].numerator)/ab[0].denominator + mp.mpf(ab[1].numerator)/ab[1].denominator*s
    lv={n: mp.log(abs(val(g))) for n,g in GEN.items()}; lv['pi2']=mp.pi**2; lv['pi4']=mp.pi**4
    L5=mp.fsum(chi13(a)*mp.zeta(5, mp.mpf(a)/13) for a in range(1,13))/mp.mpf(13)**5
    tot = T_coef.numerator*L5/s/T_coef.denominator + z5_coef.numerator*mp.zeta(5)/z5_coef.denominator
    for ab,c in li5:
        tot += c.numerator*pli5(val(ab), dps+10)/c.denominator
    for parts,c in prods:
        tot += c.numerator*mp.fprod(lv[x] for x in parts)/c.denominator
    print(f'  residual @{dps}: {mp.nstr(tot,5)}   ({time.time()-t0:.0f}s)', flush=True)
# EXACT symbol over the 5 letters (integer arithmetic)
sym=defaultdict(Ff); miss=0
for ab,c in li5:
    if ab not in pool: miss+=1; continue
    ex,e1=pool[ab]
    for i0 in range(5):
        c0=-e1[i0]
        if not c0: continue
        for t in itertools.product(range(5), repeat=4):
            w=c0
            for j in t: w*=ex[j]
            if w: sym[(i0,)+t]+=c*w
for parts,c in prods:
    if 'pi2' in parts or 'pi4' in parts: continue
    ms=[LET.index(p) for p in parts]
    m=1; from collections import Counter
    for k2 in Counter(ms).values(): m*=math.factorial(k2)
    for p2 in set(itertools.permutations(ms)): sym[p2]+=c*m
nz={k:v for k,v in sym.items() if v}
print(f'  args matched to pool: {len(li5)-miss}/{len(li5)}; EXACT symbol nonzero entries: {len(nz)}', flush=True)
# EXACT factorisations: for each arg, x and 1-x are S-units with the pool coords
sys.path.insert(0,'.')
import w5blib as W
from w5blib import QD
F=W.make_field(13, split=[3], inert=[2])
bad=0
for ab,c in li5:
    x=QD(ab[0],ab[1],13); vx,v1m=pool[ab]
    sgn=1 if x.emb(60)>0 else -1
    xabs=x if sgn==1 else -x
    xr=QD(1,0,13)
    for g,ee in zip(F['letters'], vx):
        if ee: xr=xr*(g**int(ee))
    okx=(xr==xabs); ok1=(W.factor_sunit(1-x, F)==v1m)
    if not (okx and ok1): bad+=1
print(f'  S-unit factorisations: {len(li5)} checked, {bad} bad', flush=True)
print('VERIFY DONE', path, flush=True)
'''
open('verify_relation.py','w').write(script)
import subprocess, os
log = open('verify_gcd.log','w')
p = subprocess.Popen(['nice','-n','10','<project>/.venv/bin/python',
                      '-c','import subprocess,sys\nfor f in ["d13_gcd_minnum.json","d13_gcd.json"]:\n subprocess.run([sys.executable,"verify_relation.py",f])'],
                     stdout=log, stderr=subprocess.STDOUT, cwd='.')
open('verify_gcd_pid.txt','w').write(str(p.pid))
print('launched pid', p.pid)
