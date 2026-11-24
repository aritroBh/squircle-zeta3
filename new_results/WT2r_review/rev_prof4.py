
import sys, time
sys.path.insert(0,'.')
from math import comb, factorial
from fractions import Fraction as Fr
from genintegrand import S_poly
import buildintegrand as bi
from fastred import *
IQ = QSI(0,0,1,0)
t0=time.time()
for m in range(0, 11):
    tm=time.time()
    S = S_poly(10,m)
    St, D = bi.Stilde(S)
    print(f"m={m}: S_poly+Stilde {time.time()-tm:.2f}s (D={D}, len(St)={len(St)})", flush=True)
print("total", time.time()-t0)
