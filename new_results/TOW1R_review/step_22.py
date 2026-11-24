
import sys, os
sys.dont_write_bytecode=True
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from myengine import *
# I3: (4p-3) rho_{p,p-1} = (p+1) tau1 rho_{p+1,p-1} + p tau2 (rho_pp - conj rho_pp)   [depth-0 adjacent]
def tau2(n):
    s=n+1; return ((-1)**s)*((IQ-AQ)**(-(s+1))+(IQ+BQ)**(-(s+1))-IQ**(-(s+1)))
bad=0
for p in range(2,20):
    d=2*p-1
    lhs = rho(p,p-1,0).scale(4*p-3)
    rhs = rho(p+1,p-1,0)*tau1(0).scale(p+1) + (rho(p,p,0)-rho(p,p,0).conj_i())*tau2(0).scale(p)
    if lhs!=rhs: bad+=1; print("I3 fail",p)
print("I3 (adjacent stabilization), p=2..19:", "0 bad" if bad==0 else f"{bad} bad")
# I4: hhat_p = (p+1) tau1 rho_{p+1,p}/(4p-1) = (1-i) rho_pp ; and hhat-conj(hhat) = rho_pp-conj(rho_pp)
bad=0
for p in range(1,20):
    hhat = rho(p+1,p,0)*tau1(0).scale(Fr(p+1,4*p-1))
    target = rho(p,p,0)*(QSI(1)-IQ)
    if hhat!=target: bad+=1; print("I4a fail",p, hhat.sympy(), target.sympy())
    if (hhat-hhat.conj_i()) != (rho(p,p,0)-rho(p,p,0).conj_i()): bad+=1; print("I4b fail",p)
print("I4 (diagonal stabilization), p=1..19:", "0 bad" if bad==0 else f"{bad} bad")
# I2: non-adjacent depth-0: (2d-1) rho_{pq} = (p+1) tau1 rho_{p+1,q} + (q+1) tau2 rho_{p,q+1}
bad=0
for d in range(4,30):
    for p in range(d,(d)//2,-1):
        q=d-p
        if p<q+2: continue
        lhs=rho(p,q,0).scale(2*d-1)
        rhs=rho(p+1,q,0)*tau1(0).scale(p+1)+rho(p,q+1,0)*tau2(0).scale(q+1)
        if lhs!=rhs: bad+=1
print("I2 (non-adjacent depth-0), d=4..29:", "0 bad" if bad==0 else f"{bad} bad")
