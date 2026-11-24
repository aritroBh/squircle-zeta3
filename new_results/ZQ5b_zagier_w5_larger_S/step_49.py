import os, json
from fractions import Fraction as Ff
os.chdir('<project>/research/08_theory/ZQ5b_zagier_w5_larger_S')
step = r'''
import sys, json, time, random
from fractions import Fraction as Ff
sys.path.insert(0,'.')
S = json.load(open('d17_min_setup.json'))
N=S['N']; M=S['M']; L=1+N+M
c0=[Ff(x) for x in S['c_rel']]
Z=[[Ff(x) for x in row] for row in S['Z']]
nZ=len(Z)
Zsp=[[(k,v) for k,v in enumerate(row) if v] for row in Z]
rows_with=[[] for _ in range(L)]
for zi,zs in enumerate(Zsp):
    for k,v in zs: rows_with[k].append(zi)
def maxnum(c): return max(abs(x.numerator) for x in c)
def nli5(c): return sum(1 for i in range(N) if c[1+i])
def nterms(c): return sum(1 for x in c if x)
def apply(c,zi,t):
    cc=list(c)
    for k,v in Zsp[zi]: cc[k]=cc[k]+t*v
    return cc
def reduce_terms(c, cap, t0, budget, rng_li5=True):
    c=list(c); rng=range(1,1+N) if rng_li5 else range(L)
    while True:
        if time.time()-t0>budget: break
        best=None
        for i in rng:
            if c[i]==0: continue
            for zi in rows_with[i]:
                zv=dict(Zsp[zi])
                t=-c[i]/zv[i]
                cand=apply(c,zi,t)
                co=maxnum(cand)
                if co<=cap and (best is None or co<best[0]): best=(co,cand,i)
        if best is None: break
        c=best[1]
    return c
t0=time.time()
bestA=[Ff(x) for x in json.load(open('d17_min_maxnum.json'))['c']]
print("start maxnum sol:", maxnum(bestA), "li5", nli5(bestA), "terms", nterms(bestA))
# term reduction from c0 (94 li5) at several caps
results={}
for cap in [90000, 120000, 180000, 232217]:
    cb=reduce_terms(c0, cap, t0, 1400)
    results[cap]=(maxnum(cb), nterms(cb), nli5(cb))
    print(f"cap {cap}: maxnum={maxnum(cb)} terms={nterms(cb)} li5={nli5(cb)} ({time.time()-t0:.0f}s)", flush=True)
    json.dump(dict(c=[str(x) for x in cb]), open(f'd17_min_cap{cap}.json','w'))
# also reduce terms from the maxnum solution at cap 90000
cb=reduce_terms(bestA, 90000, t0, 1400)
print(f"from maxnum, cap 90000: maxnum={maxnum(cb)} terms={nterms(cb)} li5={nli5(cb)}")
json.dump(dict(c=[str(x) for x in cb]), open('d17_min_cap90000b.json','w'))
print("done")
'''
open('step_04f_d17_minB.py','w').write(step)
import subprocess
r = subprocess.run(['nice','-n','10','<project>/.venv/bin/python','step_04f_d17_minB.py'], capture_output=True, text=True, cwd='.', timeout=1700)
print(r.stdout); print("ERR:", r.stderr[-1500:])
