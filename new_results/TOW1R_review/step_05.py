
import sys, os
sys.dont_write_bytecode = True
WTS = "<project>/research/08_theory/WTS_weight_k_structure/theory/"
sys.path.insert(0, WTS)
from fractions import Fraction as Fr
from math import comb, factorial
from fastred import *
from buildintegrand import build_integrand_pf

def doublefact(n):  # (2d-3)!!
    r=1
    for x in range(n,0,-2): r*=x
    return r
def Kd(d):
    return QSI(0,0,0,0)+QSI(Fr((-1)**d * doublefact(2*d-3), 2*factorial(d))) * (IQ**(d+1)) / SQ2Q
# Wait: K_d = (-1)^d (2d-3)!! i^{d+1}/(2 sqrt2 d!). 1/sqrt2 = sqrt2/2.
def Kd2(d):
    return (QSI(Fr((-1)**d * doublefact(2*d-3), 2*factorial(d))) * (IQ**(d+1))) * QSI(0,Fr(1,2))  # *(sqrt2/2)
def r_p(d,p):
    return QSI(2*(d-2*p),0,-2*d,0) - QSI(Fr(4*(d-1)*(d-2*p), 2*d-3))
def rho_closed(p,q,delta):
    d=p+q; mult = 2 if p>q else 1
    r0 = Kd2(d) * QSI(Fr((-1)**p * comb(d,p) * mult))
    if delta==0: return r0
    if delta==1: return r0 * QSI(0,0,-(d-2),0)
    if delta==2:
        val = QSI(Fr(3*d-3,2)) - r_p(d,p)*(IQ.scale(Fr(1,2))) - QSI((d-2)**2)
        return (r0*val).scale(Fr(1,2))

bad=0; tested=0
for d in range(2,13):
    ipf = build_integrand_pf(d)
    for p in range(d, (d-1)//2, -1):
        q=d-p
        P = ipf[(p,q)]
        ipos = P[1]
        for delta in (0,1,2):
            r = 2*d-2-delta
            if r<2: continue
            tested+=1
            data = ipos[r-1] if r-1 < len(ipos) else ZERO
            cl = rho_closed(p,q,delta)
            if data != cl:
                bad+=1
                print(f"MISMATCH d={d} key=({p},{q}) delta={delta}")
                print("   data:", data.sympy())
                print("   closed:", cl.sympy())
print(f"initial-data check d=2..12: {tested} cells, {bad} bad")
