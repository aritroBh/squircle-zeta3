
import sys, time
sys.path.insert(0,'.')
from math import comb, factorial
from fractions import Fraction as Fr
from fastred import *

IQ = QSI(0,0,1,0); SQ = QSI(0,1,0,0)
tau1 = QSI(-1,0,-1,0)/QSI(2); tau2 = QSI(-1,0,1,0)/QSI(2)
def Kf(d):
    odd=1
    for m in range(1,2*d-2,2): odd*=m
    return (IQ**(d+1))*QSI(Fr((-1)**d*odd,2*factorial(d)))/SQ
def rho_cf(p,q):
    d=p+q
    return Kf(d)*((-ONE)**p)*QSI(comb(d,p)*(1+(1 if p>q else 0)))

t0=time.time()
okI1 = all( Kf(d+1)*QSI(d+1) == Kf(d)*QSI(0,0,-(2*d-1),0) for d in range(2,401) )
print(f"(I1) K_{{d+1}}(d+1) = -i(2d-1)K_d  exact d=2..400: {okI1}  [{time.time()-t0:.0f}s]")

badI2=[]
for d in range(2,301):
    for q in range(0,d):
        p=d-q
        if p < q+2: continue
        L = rho_cf(p,q)*QSI(2*d-1)
        R = rho_cf(p+1,q)*QSI(p+1)*tau1 + rho_cf(p,q+1)*QSI(q+1)*tau2
        if L!=R: badI2.append((p,q))
print(f"(I2) non-adjacent recursion, d=2..300 all q: {'ALL HOLD' if not badI2 else badI2[:5]}  [{time.time()-t0:.0f}s]")

badI3=[]
for p in range(2,201):
    L = rho_cf(p,p-1)*QSI(4*p-3)
    rpp = rho_cf(p,p)
    R = rho_cf(p+1,p-1)*QSI(p+1)*tau1 + (rpp - rpp.conj_i())*QSI(p)*tau2
    if L!=R: badI3.append(p)
print(f"(I3) adjacent recursion, p=2..200: {'ALL HOLD' if not badI3 else badI3[:5]}  [{time.time()-t0:.0f}s]")

badI4=[]
for p in range(1,201):
    hat = rho_cf(p+1,p)*QSI(p+1)*tau1/QSI(4*p-1)
    if (hat-hat.conj_i()) != (rho_cf(p,p)-rho_cf(p,p).conj_i()): badI4.append(p)
print(f"(I4) Im(hat_p)=Im(rho_pp), p=1..200: {'ALL HOLD' if not badI4 else badI4[:5]}  [{time.time()-t0:.0f}s]")
