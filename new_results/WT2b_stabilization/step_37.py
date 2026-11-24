import sys, pickle
sys.path.insert(0,'.')
from fastred import *
from fractions import Fraction as Fr
record = pickle.load(open('record_W2b.pkl','rb'))
from math import comb, factorial
IQ = QSI(0,0,1,0); SQ = QSI(0,1,0,0)
def K(d):
    odd = 1
    for m in range(1,2*d-2,2): odd*=m   # (2d-3)!!
    return (IQ**(d+1)) * QSI(Fr((-1)**d * odd, 2*factorial(d))) / SQ
def rho_cf(p,q):
    d=p+q
    return K(d) * ((-ONE)**p) * QSI(comb(d,p)*(1+(1 if p>q else 0)))

bad=[]
for d in range(2,13):
    for p in range(d,(d-1)//2,-1):
        q=d-p
        rec = record[(d,(p,q))]
        H = rec['H'][1]; top=2*d-2
        rv = H[top-1] if len(H)>=top else ZERO
        if rv!=rho_cf(p,q): bad.append((p,q,rv.sympy(),rho_cf(p,q).sympy()))
print("closed form rho_pq vs data (all canonical keys, d=2..12):", "ALL MATCH" if not bad else bad)

tau1 = QSI(-1,0,-1,0)/QSI(2)   # (-1-i)/2
tau2 = QSI(-1,0,1,0)/QSI(2)    # (-1+i)/2
den = IQ*(IQ-BQ)*(IQ+AQ)
t1 = QSI(0,0,2,0)/den
print("tau1_0 =", t1.sympy(), " == (-1-i)/2:", t1==tau1)

okI1 = all( K(d+1)*QSI(d+1) == K(d)*QSI(0,0,-(2*d-1),0) for d in range(2,14))
print("(I1) K_{d+1}(d+1) = -i(2d-1) K_d  for d=2..13:", okI1)
badI2=[]
for d in range(2,14):
    for q in range(0,d):
        p=d-q
        if p < q+2: continue
        L = rho_cf(p,q)*QSI(2*d-1)
        R = rho_cf(p+1,q)*QSI(p+1)*tau1 + rho_cf(p,q+1)*QSI(q+1)*tau2
        if L!=R: badI2.append((p,q))
print("(I2) ID1 non-adjacent off-diagonal, d<=14:", "ALL HOLD" if not badI2 else badI2)
badI3=[]
for p in range(2,8):
    L = rho_cf(p,p-1)*QSI(4*p-3)
    rpp = rho_cf(p,p)
    R = rho_cf(p+1,p-1)*QSI(p+1)*tau1 + (rpp - rpp.conj_i())*QSI(p)*tau2
    if L!=R: badI3.append((p,L.sympy(),R.sympy()))
print("(I3) ID-adjacent, p=2..7:", "ALL HOLD" if not badI3 else badI3)
badI4=[]
for p in range(1,8):
    hat = rho_cf(p+1,p)*QSI(p+1)*tau1/QSI(4*p-1)
    if (hat-hat.conj_i()) != (rho_cf(p,p)-rho_cf(p,p).conj_i()): badI4.append(p)
print("(I4) Im hat = Im rho_pp, p=1..7:", "ALL HOLD" if not badI4 else badI4)
