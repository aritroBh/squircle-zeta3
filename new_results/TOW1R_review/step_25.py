
import sys, os, time
sys.dont_write_bytecode=True
WTS="<project>/research/08_theory/WTS_weight_k_structure/theory/"
sys.path.insert(0,WTS)
from fractions import Fraction as Fr
from math import comb, factorial
import fastred
from fastred import *
from buildintegrand import build_integrand_pf
def build(k):
    ipf=build_integrand_pf(k); work={}
    for key,R in ipf.items():
        kk,RR=fastred.fold_key(key,R); work[kk]=pf_add(work.get(kk,pf_zero()),RR)
    Hrec={}; Gm={}
    for deg in range(k,-1,-1):
        for key in sorted([kk for kk in work if sum(kk)==deg],reverse=True):
            R=work.pop(key)
            if deg==0: continue
            simple,higher,const=pf_split(R)
            Hrec[key]={r:e for r,e in enumerate(higher[1],start=1) if not e.is_zero()}
            if not higher[1] and not higher[2]: continue
            G=fastred.pf_integrate_higher(higher); C=pf_eval_at(G,BQ); Gm[key]=[G[0],G[1],G[2],-C]
            i,j=key
            if i>=1:
                kk,RR=fastred.fold_key((i-1,j),pf_scale(pf_mul_DL(Gm[key],DL1),Fr(-i))); work[kk]=pf_add(work.get(kk,pf_zero()),RR)
            if j>=1:
                kk,RR=fold_key((i,j-1),pf_scale(pf_mul_DL(Gm[key],DL2),Fr(-j))); work[kk]=pf_add(work.get(kk,pf_zero()),RR)
    return Hrec,Gm
def u(Hrec,key,delta):
    d=sum(key); return Hrec.get(key,{}).get(2*d-2-delta,ZERO)
def doublefact(n):
    r=1
    for x in range(n,0,-2): r*=x
    return r
def Kd(dd): return (QSI(Fr((-1)**dd*doublefact(2*dd-3),2*factorial(dd)))*(IQ**(dd+1)))*QSI(0,Fr(1,2))
def rp(dd,pp): return QSI(2*(dd-2*pp),0,-2*dd,0)-QSI(Fr(4*(dd-1)*(dd-2*pp),2*dd-3))
def rho(pp,qq,delta):
    dd=pp+qq; mult=2 if pp>qq else 1
    r0=Kd(dd)*QSI(Fr((-1)**pp*comb(dd,pp)*mult))
    if delta==0: return r0
    if delta==1: return r0*QSI(0,0,-(dd-2),0)
    val=QSI(Fr(3*dd-3,2))-rp(dd,pp)*(IQ.scale(Fr(1,2)))-QSI((dd-2)**2)
    return (r0*val).scale(Fr(1,2))
def U_e(pp,qq,k,delta):
    dd=pp+qq
    if delta==0: return rho(pp,qq,0)
    if delta==1: return rho(pp,qq,1)
    if delta==2: return rho(pp,qq,2)+rho(pp,qq,0).scale(Fr(-(k-dd),2*dd-3))
s30=[{},[QSI(0,Fr(1,12))],[QSI(0,Fr(1,12))],QSI(Fr(-1,12))]
s21=[{},[QSI(0,Fr(-1,4))],[QSI(0,Fr(-1,4))],QSI(Fr(1,4))]
for k in [18,20]:
    t0=time.time(); Hrec,Gm=build(k)
    # degree-3 tower law
    bad=0
    for key in [(3,0),(2,1)]:
        for delta in (0,1,2):
            if u(Hrec,key,delta)!=U_e(*key,k,delta): bad+=1; print("  tower mismatch",k,key,delta)
    # affine Gm law: Gm^(k) = Gm^(3)+(k-3)s  -- need Gm^(3); build once at k=3
    print(f"k={k}: degree-3 tower law δ≤2: {bad} bad  ({time.time()-t0:.0f}s)")
    # stash Gm30,Gm21 for affine check across k via sympy compare
    import sympy as sp
    for key,s in (((3,0),s30),((2,1),s21)):
        # compare Gm[k]-Gm[3] computed now vs (k-3)*s : but Gm^(3) computed at k=3 below
        pass
    pickle_gm = {kk: pf_to_sympy(Gm[kk]) for kk in [(3,0),(2,1)]}
    print("   Gm30(k) =", sp.simplify(list(pickle_gm.values())[0]))
