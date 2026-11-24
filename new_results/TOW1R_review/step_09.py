
import sys, os, time
sys.dont_write_bytecode = True
WTS = "<project>/research/08_theory/WTS_weight_k_structure/theory/"
sys.path.insert(0, WTS)
from fractions import Fraction as Fr
from math import comb, factorial
from fastred import *

def doublefact(n):
    r=1
    for x in range(n,0,-2): r*=x
    return r
def Kd(d): return (QSI(Fr((-1)**d*doublefact(2*d-3),2*factorial(d)))*(IQ**(d+1)))*QSI(0,Fr(1,2))
def rp(d,p): return QSI(2*(d-2*p),0,-2*d,0)-QSI(Fr(4*(d-1)*(d-2*p),2*d-3))
def rho(p,q,delta):
    d=p+q; mult=2 if p>q else 1
    r0=Kd(d)*QSI(Fr((-1)**p*comb(d,p)*mult))
    if delta==0: return r0
    if delta==1: return r0*QSI(0,0,-(d-2),0)
    if delta==2:
        val=QSI(Fr(3*d-3,2))-rp(d,p)*(IQ.scale(Fr(1,2)))-QSI((d-2)**2)
        return (r0*val).scale(Fr(1,2))
def tau1(n):
    s=n+1; return ((-1)**s)*((IQ-BQ)**(-(s+1))+(IQ+AQ)**(-(s+1))-IQ**(-(s+1)))
def tau2(n):
    s=n+1; return ((-1)**s)*((IQ-AQ)**(-(s+1))+(IQ+BQ)**(-(s+1))-IQ**(-(s+1)))
def e_off(p,q,delta,s):
    if (delta,s)==(0,0): return rho(p,q,0)
    if (delta,s)==(1,0): return rho(p,q,1)
    if (delta,s)==(2,0): return rho(p,q,2)
    if (delta,s)==(2,1): return rho(p,q,0).scale(Fr(-1,2*(p+q)-3))
    return ZERO
def ehat(p,delta,s):
    P=p+1; par=(p+1,p)
    if (delta,s)==(0,0): return e_off(*par,0,0)*tau1(0).scale(Fr(P,4*p-1))
    if (delta,s)==(1,0): return (e_off(*par,1,0)*tau1(0).scale(Fr(P,4*p-2))+e_off(*par,0,0)*tau1(1).scale(Fr(P,4*p-1)))
    if (delta,s)==(2,0): return (e_off(*par,2,0)*tau1(0).scale(Fr(P,4*p-3))+e_off(*par,1,0)*tau1(1).scale(Fr(P,4*p-2))+e_off(*par,0,0)*tau1(2).scale(Fr(P,4*p-1)))
    if (delta,s)==(2,1): return e_off(*par,2,1)*tau1(0).scale(Fr(P,4*p-3))
    return ZERO
def e_any(p,q,delta,s):
    return ehat(p,delta,s) if p==q else e_off(p,q,delta,s)
def Ediag_conv(p,q,delta,sigma):
    if p!=q: return e_off(p,q,delta,sigma)
    out=ZERO; t=0
    while sigma+t<=delta//2:
        out = out + ehat(p,delta,sigma+t).scale((-1)**t); t+=1
    return out
def LambdaE(p,q,delta,sigma):
    d=p+q; tot=ZERO
    for dp in range(0,delta+1):
        tot = tot + e_any(p+1,q,dp,sigma)*tau1(delta-dp).scale(Fr(p+1,2*d-dp-1))
        if p>q:
            tot = tot + Ediag_conv(p,q+1,dp,sigma)*tau2(delta-dp).scale(Fr(q+1,2*d-dp-1))
    if q==p-1:
        m=4*p-4-delta; ft=ZERO
        for mu in range(m,4*p-3):
            g=ZERO
            for dp in range(0,4*p-4-mu+1):
                g = g - Ediag_conv(p,p,dp,sigma)*tau1(4*p-4-mu-dp).scale(Fr(1,4*p-3-dp))
            ft = ft + g.conj_i()*(IQ**(-mu))*Fr(comb(mu-2,mu-m))*(IQ**(m-2))
        tot = tot + ft.scale(-p)
    return tot

t0=time.time(); bad=0; tested=0
for d in range(2,26):
    for p in range(d,(d-1)//2,-1):
        q=d-p
        if p==q: continue
        for delta in (0,1,2):
            if delta>2*d-4: continue
            for sigma in range(0,delta//2+1):
                tested+=1
                lhs = e_off(p,q,delta,sigma)+e_off(p,q,delta,sigma+1)
                rhs = LambdaE(p,q,delta,sigma)
                if lhs!=rhs:
                    bad+=1
                    if bad<=6: print("  MISMATCH",(p,q),delta,sigma)
print(f"Newton identities E1,E2a,E2b (my encoding), d=2..25: {tested} cells, {bad} bad  ({time.time()-t0:.0f}s)")
