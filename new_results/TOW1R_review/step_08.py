
import sys, os, time
sys.dont_write_bytecode = True
WTS = "<project>/research/08_theory/WTS_weight_k_structure/theory/"
sys.path.insert(0, WTS)
from fractions import Fraction as Fr
from math import comb, factorial
from fastred import *
from buildintegrand import build_integrand_pf

# ---- my own tower builder (fastred + record higher part) ----
def build_tower(k):
    ipf = build_integrand_pf(k)
    work = {}
    for key,R in ipf.items():
        kk,RR = fold_key(key,R); work[kk]=pf_add(work.get(kk,pf_zero()),RR)
    Hrec={}
    for deg in range(k,-1,-1):
        for key in sorted([kk for kk in work if sum(kk)==deg], reverse=True):
            R=work.pop(key)
            if deg==0: continue
            simple,higher,const = pf_split(R)
            Hrec[key]={r:e for r,e in enumerate(higher[1],start=1) if not e.is_zero()}
            if not higher[1] and not higher[2]: continue
            G=pf_integrate_higher(higher); C=pf_eval_at(G,BQ); Gm=[G[0],G[1],G[2],-C]
            i,j=key
            if i>=1:
                kk,RR=fold_key((i-1,j),pf_scale(pf_mul_DL(Gm,DL1),Fr(-i))); work[kk]=pf_add(work.get(kk,pf_zero()),RR)
            if j>=1:
                kk,RR=fold_key((i,j-1),pf_scale(pf_mul_DL(Gm,DL2),Fr(-j))); work[kk]=pf_add(work.get(kk,pf_zero()),RR)
    return Hrec
def u_of(H,key,delta):
    d=sum(key); r=2*d-2-delta
    return H.get(key,{}).get(r,ZERO)

# ---- independent closed-form e / U_e ----
def doublefact(n):
    r=1
    for x in range(n,0,-2): r*=x
    return r
def Kd(d):
    return (QSI(Fr((-1)**d*doublefact(2*d-3),2*factorial(d)))*(IQ**(d+1)))*QSI(0,Fr(1,2))
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
def U_e(p,q,k,delta):
    d=p+q
    if p==q:
        if k==2*p: return rho(p,p,delta)
        return sum((ehat(p,delta,s).scale(comb(k-2*p-1,s)) for s in range(0,delta//2+1)),ZERO)
    return sum((e_off(p,q,delta,s).scale(comb(k-d,s)) for s in range(0,delta//2+1)),ZERO)

for k in [15,16,17]:
    t0=time.time(); H=build_tower(k)
    bad=0; tested=0
    for d in range(2,k):
        for p in range(d,(d-1)//2,-1):
            q=d-p
            for delta in (0,1,2):
                if delta>2*d-4: continue
                onlaw = (p>q) or (k>=2*p+1)
                if not onlaw: continue
                tested+=1
                if U_e(p,q,k,delta)!=u_of(H,(p,q),delta):
                    bad+=1
                    if bad<=5: print("  MISMATCH",k,(p,q),delta)
    print(f"k={k}: tower law U_e==data, {tested} cells, {bad} bad  ({time.time()-t0:.0f}s)")
