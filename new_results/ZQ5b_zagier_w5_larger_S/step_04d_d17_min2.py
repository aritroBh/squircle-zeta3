
import sys, json, time, random
from fractions import Fraction as Ff
sys.path.insert(0,'.')
S = json.load(open('d17_min_setup.json'))
N=S['N']; M=S['M']; L=1+N+M
c0=[Ff(x) for x in S['c_rel']]
Z=[[Ff(x) for x in row] for row in S['Z']]
nZ=len(Z)
def maxnum(c): return max(abs(x.numerator) for x in c)
def nli5(c): return sum(1 for i in range(N) if c[1+i])
def nterms(c): return sum(1 for x in c if x)
def reduce_terms(c, cap, alldim=True):
    c=list(c); rng=range(L) if alldim else range(1,1+N)
    while True:
        best=None
        for i in rng:
            if c[i]==0: continue
            for z in Z:
                if z[i]==0: continue
                t=-c[i]/z[i]
                cand=[c[k]+t*z[k] for k in range(L)]
                co=maxnum(cand)
                if co<=cap and (best is None or co<best[0]):
                    best=(co,cand)
        if best is None: break
        c=best[1]
    return c
def greedy_maxnum(c, maxit=400):
    c=list(c)
    for it in range(maxit):
        cur=maxnum(c); best=None
        for i in range(L):
            if c[i]==0: continue
            for z in Z:
                if z[i]==0: continue
                t=-c[i]/z[i]
                cand=[c[k]+t*z[k] for k in range(L)]
                co=maxnum(cand)
                if best is None or co<best[0]: best=(co,cand)
        if best and best[0]<cur: c=best[1]
        else: break
    return c
random.seed(11)
t0=time.time()
# stage A: minimise maxnum
cA=greedy_maxnum(c0)
best=list(cA); bestc=maxnum(best)
print(f"stage A greedy: maxnum={bestc} terms={nterms(best)} li5={nli5(best)}")
for it in range(30000):
    cand=list(best)
    for _ in range(random.randint(1,5)):
        z=Z[random.randrange(nZ)]
        t=Ff(random.randint(-4,4), random.choice([1,1,1,2,3]))
        cand=[cand[k]+t*z[k] for k in range(L)]
    co=maxnum(cand)
    if co<=bestc:
        if co<bestc: bestc=co
        best=cand
best=greedy_maxnum(best)
print(f"stage A +random: maxnum={maxnum(best)} terms={nterms(best)} li5={nli5(best)} ({time.time()-t0:.0f}s)")
# stage B: term reduction at caps
frontier=[]
for cap in [maxnum(best), int(maxnum(best)*1.5), int(maxnum(best)*3), 232217]:
    cb=reduce_terms(best if cap<232217 else c0, cap)
    frontier.append((cap, maxnum(cb), nterms(cb), nli5(cb), cb))
    print(f"cap {cap}: maxnum={maxnum(cb)} total terms={nterms(cb)} li5={nli5(cb)}")
# choose: best maxnum solution and a balanced one
json.dump(dict(c=[str(x) for x in best]), open('d17_min_maxnum.json','w'))
# balanced: cap*3 row
cbal = frontier[2][4]
json.dump(dict(c=[str(x) for x in cbal]), open('d17_min_balanced.json','w'))
print("saved d17_min_maxnum.json and d17_min_balanced.json")
