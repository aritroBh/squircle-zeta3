import sys, math, time
sys.path.insert(0,'<project>/research/08_theory/DCYA1_corner_decay_proof/')
from arb_ei import E_I_sc, iv
def thin_cell(N, dw=0.1):
    s=math.sqrt(N); s2=s+dw
    n1=int(math.floor(s*s)); n2=int(math.ceil(s2*s2))
    if n2<=n1: n2=n1+1
    return iv(str(n1),str(n2))
t=time.time()
best=(0,0)
for N in range(200, 2601, 100):
    v=float(E_I_sc(thin_cell(N)))
    if v>best[1]: best=(N,v)
    print(f"N={N}: {v:.2f}")
print("PEAK:", best, "time:", round(time.time()-t,1),"s")
