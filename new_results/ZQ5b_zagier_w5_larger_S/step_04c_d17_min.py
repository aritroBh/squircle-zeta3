
import sys, json, time, random
from fractions import Fraction as Ff
sys.path.insert(0,'.')
S = json.load(open('d17_min_setup.json'))
N=S['N']; M=S['M']
L=1+N+M
c0=[Ff(x) for x in S['c_rel']]
Z=[[Ff(x) for x in row] for row in S['Z']]
def costA(c):   # max |numerator|
    return max(abs(x.numerator) for x in c)
def costB(c):   # (nonzero Li5, max num)
    return (sum(1 for i in range(N) if c[1+i]), max(abs(x.numerator) for x in c))
print("start: costA =", costA(c0), " Li5 terms:", sum(1 for i in range(N) if c0[1+i]))
def greedy(c, costf, maxit=600, li5only=True):
    c=list(c)
    rng = range(1,1+N) if li5only else range(L)
    for it in range(maxit):
        cur = costf(c)
        best=None
        for i in rng:
            if c[i]==0: continue
            for z in Z:
                if z[i]==0: continue
                t = -c[i]/z[i]
                cand = [c[k]+t*z[k] for k in range(L)]
                co = costf(cand)
                if best is None or co<best[0]:
                    best=(co,i,t,cand)
        if best and best[0]<cur:
            c=best[3]
        else:
            break
    return c
t0=time.time()
c1=greedy(c0, costA)
print(f"greedy costA: costA={costA(c1)} Li5 terms={sum(1 for i in range(N) if c1[1+i])} ({time.time()-t0:.0f}s)")
c2=greedy(c1, costB)
print(f"greedy costB: costB={costB(c2)} costA={costA(c2)} ({time.time()-t0:.0f}s)")
# random refinement on costA then costB
random.seed(7)
best=list(c2); bestc=costA(best)
nZ=len(Z)
for it in range(4000):
    cand=list(best)
    for _ in range(random.randint(1,4)):
        z=Z[random.randrange(nZ)]
        t=Ff(random.randint(-3,3), random.choice([1,1,1,2]))
        cand=[cand[k]+t*z[k] for k in range(L)]
    co=costA(cand)
    if co<bestc:
        bestc=co; best=cand
print(f"random refine: costA={bestc} Li5 terms={sum(1 for i in range(N) if best[1+i])}")
c3=greedy(best, costB)
print(f"final greedy costB: costB={costB(c3)} costA={costA(c3)}")
# save best of c3 and c2
fin = c3 if (costB(c3)<=costB(c2)) else c2
json.dump(dict(c=[str(x) for x in fin]), open('d17_min_cand.json','w'))
used=[i for i in range(N) if fin[1+i]]
print("saved d17_min_cand.json;  Li5 terms:", len(used), " max num:", costA(fin))
