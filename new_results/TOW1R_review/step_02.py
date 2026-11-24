
import sys, os, time
sys.dont_write_bytecode = True
WTS = "<project>/research/08_theory/WTS_weight_k_structure/theory/"
sys.path.insert(0, WTS)
from fractions import Fraction as Fr
from fastred import *
from buildintegrand import build_integrand_pf

def claimed_R20(k):
    d = {}
    d[ZERO]  = QSI(Fr(1,4), Fr(-1,4))
    d[AQ]    = QSI(Fr(1,4))
    d[-AQ]   = QSI(Fr(2*k-7,4))
    d[-BQ]   = QSI(Fr(-(k-2),2))
    return d
def claimed_R11(k):
    d = {}
    d[ZERO]  = QSI(Fr(2*k-5,4), Fr(1,4))
    d[-AQ]   = QSI(Fr(-(2*k-5),2))
    return d

for k in [16,18,20]:
    t0=time.time()
    ipf = build_integrand_pf(k)
    fin = reduce_k_fast(k, ipf)
    P20 = fin[(2,0)][0]; P11 = fin[(1,1)][0]
    e20 = claimed_R20(k); e11 = claimed_R11(k)
    ok20 = all(P20.get(r,ZERO)==c for r,c in e20.items()) and all(c.is_zero() or r in e20 for r,c in P20.items())
    ok11 = all(P11.get(r,ZERO)==c for r,c in e11.items()) and all(c.is_zero() or r in e11 for r,c in P11.items())
    # residual (max abs numerator over mismatches)
    def resid(P,E):
        rs=[]
        for r in set(list(P)+list(E)):
            diff = P.get(r,ZERO)-E.get(r,ZERO)
            rs.append((r.sympy(), diff.sympy()))
        return [x for x in rs if x[1]!=0]
    print(f"k={k}: R20 match={ok20}  R11 match={ok11}   ({time.time()-t0:.1f}s)")
    if not ok20: print("   R20 mismatches:", resid(P20,e20))
    if not ok11: print("   R11 mismatches:", resid(P11,e11))
    print("   higher(2,0):", fin[(2,0)][1], fin[(2,0)][2], "higher(1,1):", fin[(1,1)][1], fin[(1,1)][2])
