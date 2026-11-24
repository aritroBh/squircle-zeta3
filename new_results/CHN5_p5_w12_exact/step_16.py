import os, sys, time
base='<project>/research/08_theory'
W=base+'/CHN5_p5_w12_exact'
os.chdir(W); sys.path.insert(0,W)
import point_engine as PE, mzv_exact as M
from fractions import Fraction as F
import pickle
# reversal invariance check at two generic points
for pt in [(1,2,3,4,5),(2,3,5,7,11)]:
    c1=PE.reduced_at(pt,12)[12]
    c2=PE.reduced_at(pt[::-1],12)[12]
    inv=all(c1[i]==c2[i] for i in range(12))
    print(pt, "full-coord reversal invariant:", inv)
    if not inv:
        for i in range(12):
            if c1[i]!=c2[i]: print("  idx",i,"diff", c1[i]-c2[i])
print("A,B,CLASS,Z reversal:", )
