import os, sys, time
base='<project>/research/08_theory'
W=base+'/CHN5_p5_w12_exact'
os.chdir(W); sys.path.insert(0,W)
import point_engine as PE, mzv_exact as M
from fractions import Fraction as F
t0=time.time()
T12=PE.table(12)
print("W12 table: nullity", T12.nullity, "n words", T12.n, f"{time.time()-t0:.1f}s")
t0=time.time()
c=PE.reduced_at((1,2,3,4,5),12)
print("reduced_at(1..5) W12:", f"{time.time()-t0:.1f}s")
print("A =", c[12][7]); print("B =", c[12][8])
print("Cz9 =", c[12][9]); print("Cz57 =", c[12][10])
print("CLASS =", c[12][9]+F(14,9)*c[12][10]); print("Z =", c[12][11])
