import sys, io, contextlib
sys.path.insert(0,'.')
import t1_direct as H
from fractions import Fraction as Fr
def run(label, dmax=12, kextra=5):
    bad=0; n=0
    for dd in range(4,dmax+1):
        for p in range(dd,(dd-1)//2,-1):
            q=dd-p
            ks = range(dd+1,dd+kextra+1) if p!=q else range(2*p+1,2*p+kextra+1)
            for k in ks:
                for delta in range(0,min(4,2*dd-4)+1):
                    n+=1
                    if not (H.U(p,q,k,delta)==H.T1rhs(p,q,k,delta)): bad+=1
    print(f'{label}: {n} cells, {bad} bad')
run('baseline')
# mutation 1: e(4,1) with wrong constant
oldF=H.F
def F1(delta,s,dd,mm):
    v=oldF(delta,s,dd,mm)
    if (delta,s)==(4,1): v=v+H.G(Fr(1,1000))
    return v
H.F=F1; run('mutate e(4,1) +1/1000'); H.F=oldF
# mutation 2: e(4,2) constant 3 instead of 2
def F2(delta,s,dd,mm):
    if (delta,s)==(4,2): return H.G(Fr(3,(2*dd-3)*(2*dd-5)))
    return oldF(delta,s,dd,mm)
H.F=F2; run('mutate e(4,2) const 3'); H.F=oldF
# mutation 3: drop E3: set rho_pp = polynomial extension at k=2p
oldU=H.U
def U3(p,q,k,delta):
    if p==q and k==2*p:
        tot=H.G(0)
        for s in range(delta//2+1):
            tot=tot+H.ehat(p,delta,s)*((-1)**s)
        return tot
    return oldU(p,q,k,delta)
H.U=U3; run('mutate: diag k=2p uses polynomial ext (E3 swapped)'); H.U=oldU
# mutation 4: perturb rho(4) at one place
oldA=H.arho
def A4(delta,dd,mm):
    v=oldA(delta,dd,mm)
    if delta==4 and dd==9 and mm==5: v=v+H.G(Fr(1,10**6))
    return v
H.arho=A4; H._arho_cache.clear(); run('mutate rho(4) at (9,5)'); H.arho=oldA; H._arho_cache.clear()
