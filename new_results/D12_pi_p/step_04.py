# step_04.py -- TASK 3: exact series arithmetic for [h^n] Psi (rational coeffs in s, L=log(1+e^-s), log2)
from fractions import Fraction
import math, pickle
import sympy as sp
KMAX = 44
def addC(A,B):
    C=dict(A)
    for k,v in B.items():
        C[k]=C.get(k,Fraction(0))+v
        if C[k]==0: del C[k]
    return C
def scaleC(A,c): return {k:v*c for k,v in A.items()}
def mulC(A,B):
    C={}
    for (i1,j1,k1),v1 in A.items():
        for (i2,j2,k2),v2 in B.items():
            key=(i1+i2,j1+j2,k1+k2); C[key]=C.get(key,Fraction(0))+v1*v2
    return {k:v for k,v in C.items() if v!=0}
def serMul(A,B,K=KMAX):
    C=[{} for _ in range(K+1)]
    for i in range(len(A)):
        if not A[i]: continue
        for j in range(len(B)):
            if not B[j]: continue
            if i+j>K: break
            C[i+j]=addC(C[i+j],mulC(A[i],B[j]))
    return C
def serExp(A,K=KMAX):
    R=[{} for _ in range(K+1)]; R[0]={(0,0,0):Fraction(1)}
    for n in range(1,K+1):
        acc={}
        for k in range(1,n+1):
            if A[k]: acc=addC(acc,mulC(scaleC(A[k],Fraction(k)),R[n-k]))
        R[n]=scaleC(acc,Fraction(1,n))
    return R
y=sp.symbols('y')
lc=sp.log(sp.cosh(y)).series(y,0,KMAX+2).removeO()
c=[None]*(KMAX+1); c[1]=Fraction(-1,2)
for r in range(2,KMAX+1): c[r]=Fraction(sp.expand(lc).coeff(y,r))/Fraction(2)**r
T=[{} for _ in range(KMAX+1)]; T[0]={(0,0,1):Fraction(1)}   # log(1+e^{-sa}) = log2 + sum c_r s^r a^r
for r in range(1,KMAX+1):
    if c[r]==0: continue
    for t in range(0,KMAX+1-r):
        T[r+t]=addC(T[r+t],{(r,0,0):Fraction(math.comb(r+t-1,t))*c[r]})
Q=[{} for _ in range(KMAX+1)]
for n in range(KMAX+1):
    Q[n]=addC(Q[n],scaleC(T[n],Fraction(-1)))
    if n+1<=KMAX: Q[n+1]=addC(Q[n+1],scaleC(T[n],Fraction(-1)))
psi2=[scaleC(co,Fraction(1,2)) for co in serExp([{}]+Q[1:])]   # (1+e^{-sa})^{-1-h}, exp(Q0)=1/2
B=[{} for _ in range(KMAX+1)]
for m in range(2,KMAX+1): B[m]={(1,0,0):Fraction(-1)}
psi1=serExp(B)                                                  # exp(-s h^2/(1-h))
psi3=[{} for _ in range(KMAX+1)]
for m in range(1,KMAX+1): psi3[m]={(0,m,0):Fraction(1,math.factorial(m))}  # exp(hL)-1
PSI=serMul(serMul(psi1,psi2),psi3)
pickle.dump({'PSI':PSI,'KMAX':KMAX},open('/tmp/psi.pkl','wb'))
print("low orders:", PSI[1], PSI[2])
print("saved PSI to /tmp/psi.pkl")
