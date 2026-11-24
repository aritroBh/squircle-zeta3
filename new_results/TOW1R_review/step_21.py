
import sys, os
sys.dont_write_bytecode=True
WTS="<project>/research/08_theory/WTS_weight_k_structure/theory/"
sys.path.insert(0,WTS)
from fractions import Fraction as Fr
from math import comb
import fastred
from fastred import *
from buildintegrand import build_integrand_pf
def build_tower(k):
    ipf=build_integrand_pf(k); work={}
    for key,R in ipf.items():
        kk,RR=fastred.fold_key(key,R); work[kk]=pf_add(work.get(kk,pf_zero()),RR)
    Hrec={}
    for deg in range(k,-1,-1):
        for key in sorted([kk for kk in work if sum(kk)==deg],reverse=True):
            R=work.pop(key)
            if deg==0: continue
            simple,higher,const=pf_split(R)
            Hrec[key]={r:e for r,e in enumerate(higher[1],start=1) if not e.is_zero()}
            if not higher[1] and not higher[2]: continue
            G=fastred.pf_integrate_higher(higher); C=pf_eval_at(G,BQ); Gm=[G[0],G[1],G[2],-C]
            i,j=key
            if i>=1:
                kk,RR=fastred.fold_key((i-1,j),pf_scale(pf_mul_DL(Gm,DL1),Fr(-i))); work[kk]=pf_add(work.get(kk,pf_zero()),RR)
            if j>=1:
                kk,RR=fold_key((i,j-1),pf_scale(pf_mul_DL(Gm,DL2),Fr(-j))); work[kk]=pf_add(work.get(kk,pf_zero()),RR)
    return Hrec
def u(Hrec,key,delta):
    d=sum(key); return Hrec.get(key,{}).get(2*d-2-delta,ZERO)
# closed forms
import sys as _s
def doublefact(n):
    r=1
    for x in range(n,0,-2): r*=x
    return r
from math import factorial
def Kd(dd): return (QSI(Fr((-1)**dd*doublefact(2*dd-3),2*factorial(dd)))*(IQ**(dd+1)))*QSI(0,Fr(1,2))
def rp(dd,pp): return QSI(2*(dd-2*pp),0,-2*dd,0)-QSI(Fr(4*(dd-1)*(dd-2*pp),2*dd-3))
def rho(pp,qq,delta):
    dd=pp+qq; mult=2 if pp>qq else 1
    r0=Kd(dd)*QSI(Fr((-1)**pp*comb(dd,pp)*mult))
    if delta==0: return r0
    if delta==1: return r0*QSI(0,0,-(dd-2),0)
    val=QSI(Fr(3*dd-3,2))-rp(dd,pp)*(IQ.scale(Fr(1,2)))-QSI((dd-2)**2)
    return (r0*val).scale(Fr(1,2))
def tau1(n):
    s=n+1; return ((-1)**s)*((IQ-BQ)**(-(s+1))+(IQ+AQ)**(-(s+1))-IQ**(-(s+1)))
def ehat(pp,delta,s):
    P=pp+1
    if (delta,s)==(0,0): return rho(pp+1,pp,0)*tau1(0).scale(Fr(P,4*pp-1))
    if (delta,s)==(1,0): return (rho(pp+1,pp,1)*tau1(0).scale(Fr(P,4*pp-2))+rho(pp+1,pp,0)*tau1(1).scale(Fr(P,4*pp-1)))
    if (delta,s)==(2,0): return (rho(pp+1,pp,2)*tau1(0).scale(Fr(P,4*pp-3))+rho(pp+1,pp,1)*tau1(1).scale(Fr(P,4*pp-2))+rho(pp+1,pp,0)*tau1(2).scale(Fr(P,4*pp-1)))
    if (delta,s)==(2,1): return rho(pp+1,pp,0).scale(Fr(-1,2*(2*pp+1)-3))*tau1(0).scale(Fr(P,4*pp-3))
    return ZERO

# towers
Tw = {k: build_tower(k) for k in range(2,10)}
# 1) diagonal base case k=2p+1: u_pp^(2p+1)(delta) == ehat_pp(delta,0)
print("--- diagonal base case u_pp^(2p+1)==ehat(delta,0) ---")
bad=0
for p in range(1,5):
    k=2*p+1
    for delta in (0,1,2):
        if delta>2*(2*p)-4: continue
        got=u(Tw[k],(p,p),delta); exp=ehat(p,delta,0)
        if got!=exp: bad+=1; print(f"  p={p} delta={delta}: got {got.sympy()} exp {exp.sympy()}")
print("  ", "0 bad" if bad==0 else f"{bad} bad")
# 2) diagonal tower law k>=2p+1: u=sum_s binom(k-2p-1,s) ehat(delta,s)
print("--- diagonal tower law k=2p+1..2p+4 ---")
bad=0
for p in range(1,5):
    for k in range(2*p+1, min(2*p+5,10)):
        for delta in (0,1,2):
            if delta>2*(2*p)-4: continue
            got=u(Tw[k],(p,p),delta)
            exp=sum((ehat(p,delta,s).scale(comb(k-2*p-1,s)) for s in range(0,delta//2+1)),ZERO)
            if got!=exp: bad+=1; print(f"  p={p} k={k} delta={delta}: got {got.sympy()} exp {exp.sympy()}")
print("  ", "0 bad" if bad==0 else f"{bad} bad")
# 3) blindness: adjacent child (p,p-1) at k=2p equals tower-law value rho_child(delta)+...
print("--- blindness: adjacent child at k=2p ---")
bad=0
for p in range(2,6):
    k=2*p
    if k>9: continue
    for delta in (0,1,2):
        if delta>2*(2*p-1)-4: continue
        got=u(Tw[k],(p,p-1),delta)
        dd=2*p-1
        # tower-law off-diagonal: rho(delta) for delta<=1; rho(2)+(k-dd)*e(2,1)
        if delta<2: exp=rho(p,p-1,delta)
        else: exp=rho(p,p-1,2)+rho(p,p-1,0).scale(Fr(-(k-dd),2*dd-3))
        if got!=exp: bad+=1; print(f"  p={p} delta={delta}: got {got.sympy()} exp {exp.sympy()}")
print("  ", "0 bad" if bad==0 else f"{bad} bad")
