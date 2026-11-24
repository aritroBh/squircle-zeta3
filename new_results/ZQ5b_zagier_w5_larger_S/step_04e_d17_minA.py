
import sys, json, time, random
from fractions import Fraction as Ff
sys.path.insert(0,'.'
)
S = json.load(open('d17_min_setup.json'))
N=S['N']; M=S['M']; L=1+N+M
c0=[Ff(x) for x in S['c_rel']]
Z=[[Ff(x) for x in row] for row in S['Z']]
nZ=len(Z)
# sparse representation of Z rows
Zsp=[]
for row in Z:
    Zsp.append([(k,v) for k,v in enumerate(row) if v])
rows_with = [[] for _ in range(L)]
for zi,zs in enumerate(Zsp):
    for k,v in zs: rows_with[k].append(zi)
def maxnum(c): return max(abs(x.numerator) for x in c)
def nli5(c): return sum(1 for i in range(N) if c[1+i])
def nterms(c): return sum(1 for x in c if x)
def apply(c, zi, t):
    cc=list(c)
    for k,v in Zsp[zi]:
        cc[k]=cc[k]+t*v
    return cc
def greedy_maxnum(c, maxit=80, t0=None, budget=1200):
    c=list(c)
    for it in range(maxit):
        if t0 and time.time()-t0>budget: break
        cur=maxnum(c); best=None
        for i in range(L):
            if c[i]==0: continue
            for zi in rows_with[i]:
                zv = dict(Zsp[zi])
                t=-c[i]/zv[i]
                cand=apply(c,zi,t)
                co=maxnum(cand)
                if best is None or co<best[0]: best=(co,cand)
        if best and best[0]<cur: c=best[1]
        else: break
    return c
t0=time.time()
cA=greedy_maxnum(c0, maxit=60, t0=t0, budget=600)
print(f"stage A greedy: maxnum={maxnum(cA)} terms={nterms(cA)} li5={nli5(cA)} ({time.time()-t0:.0f}s)", flush=True)
random.seed(11)
best=list(cA); bestc=maxnum(best)
for it in range(20000):
    if it%2000==0: print(f"  rand {it} bestc={bestc} ({time.time()-t0:.0f}s)", flush=True)
    if time.time()-t0>1100: break
    cand=list(best)
    for _ in range(random.randint(1,4)):
        zi=random.randrange(nZ)
        t=Ff(random.randint(-4,4), random.choice([1,1,1,2]))
        cand=apply(cand,zi,t)
    co=maxnum(cand)
    if co<=bestc:
        if co<bestc: bestc=co
        best=cand
best=greedy_maxnum(best, maxit=40, t0=t0, budget=1400)
print(f"stage A final: maxnum={maxnum(best)} terms={nterms(best)} li5={nli5(best)} ({time.time()-t0:.0f}s)", flush=True)
json.dump(dict(c=[str(x) for x in best]), open('d17_min_maxnum.json','w'))
print("saved d17_min_maxnum.json")
