import sys, math
sys.path.insert(0,'<project>/research/08_theory/DCYA1_corner_decay_proof/')
from arb_ei import E_I_sc, iv
def thin_cell(N, dw=0.1):
    s=math.sqrt(N); s2=s+dw
    n1=int(math.floor(s*s)); n2=int(math.ceil(s2*s2))
    if n2<=n1: n2=n1+1
    return iv(str(n1),str(n2))
best=(0,0)
for N in range(1380, 1701, 5):
    v=float(E_I_sc(thin_cell(N)))
    if v>best[1]: best=(N,v)
print("PEAK (1380-1700):", best)
# print the values around the peak
for N in range(1415,1441,5):
    print(f"N={N}: {float(E_I_sc(thin_cell(N))):.3f}")
