import os, sys, time
base='<project>/research/08_theory'
W=base+'/CHN5_p5_w12_exact'
os.chdir(W); sys.path.insert(0,W)
import point_engine as PE, mzv_exact as M
from fractions import Fraction as F
def four_coords(pt):
    c=PE.reduced_at(pt,12)[12]
    return (c[7], c[8], c[9]+F(14,9)*c[10], c[11])
tests = {
 'a2=0':  [(1,0,2,3,4),(2,0,1,5,3),(3,0,4,1,2),(5,0,2,1,3)],
 'a3=0':  [(1,2,0,3,4),(2,1,0,5,3),(3,4,0,1,2),(5,2,0,1,3)],
 'a4=0':  [(1,3,2,0,4),(2,5,1,0,3),(3,1,4,0,2),(5,3,2,0,1)],
 'Sa=0':  [(1,2,3,4,-10),(2,1,4,3,-10),(1,3,2,5,-11),(3,1,5,2,-11)],
}
t0=time.time()
allok=True
for name,pts in tests.items():
    for pt in pts:
        fc=four_coords(pt)
        ok=all(v==0 for v in fc)
        allok &= ok
        print(name, pt, "all four zero:" , ok, ("" if ok else fc))
print(f"ALL HYPERPLANE CHECKS PASS: {allok}  [{time.time()-t0:.0f}s]")
