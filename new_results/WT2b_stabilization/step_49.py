step=r'''import sys, pickle
sys.path.insert(0,'.')
import sympy as sp
from fastred import *
from fastred import v
from fractions import Fraction as Fr
from math import comb, factorial
record = pickle.load(open('record_W2b.pkl','rb'))
S2=sp.sqrt(2); I=sp.I
tau1 = QSI(-1,0,-1,0)/QSI(2); tau2 = QSI(-1,0,1,0)/QSI(2)
IQ=QSI(0,0,1,0); SQ=QSI(0,1,0,0)
def Kf(d):
    odd=1
    for m in range(1,2*d-2,2): odd*=m
    return (IQ**(d+1))*QSI(Fr((-1)**d*odd,2*factorial(d)))/SQ
def rho_cf(p,q):
    d=p+q
    return Kf(d)*((-ONE)**p)*QSI(comb(d,p)*(1+(1 if p>q else 0)))

ok = True
# 1. closed form vs all initial data
for d in range(2,13):
    for p in range(d,(d-1)//2,-1):
        q=d-p; H=record[(d,(p,q))]['H'][1]; top=2*d-2
        rv = H[top-1] if len(H)>=top else ZERO
        ok &= (rv==rho_cf(p,q))
print("1. rho_pq closed form == initial data (all keys d=2..12):", ok)
# 2. identities
okI1 = all( Kf(d+1)*QSI(d+1) == Kf(d)*QSI(0,0,-(2*d-1),0) for d in range(2,16))
okI2 = all( rho_cf(p,q)*QSI(2*(p+q)-1) == rho_cf(p+1,q)*QSI(p+1)*tau1 + rho_cf(p,q+1)*QSI(q+1)*tau2
            for d in range(2,16) for q in range(0,d) for p in [d-q] if p>=q+2 )
okI3 = all( rho_cf(p,p-1)*QSI(4*p-3) == rho_cf(p+1,p-1)*QSI(p+1)*tau1 + (rho_cf(p,p)-rho_cf(p,p).conj_i())*QSI(p)*tau2
            for p in range(2,10) )
okI4 = all( (lambda hat: (hat-hat.conj_i())==(rho_cf(p,p)-rho_cf(p,p).conj_i()))(rho_cf(p+1,p)*QSI(p+1)*tau1/QSI(4*p-1))
            for p in range(1,10) )
print("2. (I1):", okI1, " (I2):", okI2, " (I3):", okI3, " (I4):", okI4)
# 3. top-order stabilization across all data
def h(k,key):
    i,j=key; d=i+j
    rec=record.get((k,key))
    if rec is None: return None
    H=rec['H'][1]; top=2*d-2
    return H[top-1] if len(H)>=top else ZERO
def h_hat(p): return rho_cf(p+1,p)*QSI(p+1)*tau1/QSI(4*p-1)
okD=True
for k in range(2,13):
    for i in range(1,k+1):
        for j in range(0,i+1):
            d=i+j
            if d<2 or d>k: continue
            hv=h(k,(i,j))
            if hv is None: continue
            pred = rho_cf(i,j) if i>j else (rho_cf(i,i) if k==2*i else h_hat(i))
            okD &= (hv==pred)
print("3. Lemma D stabilization vs all data k<=12:", okD)
# 4. stabilized Gm's
g20   = (1+S2)/4 - S2/(2*(v**2+1))
g11   = -S2*(v**2+2*v-1)/(4*(v**2+1))
g20_pf = pf_from_sympy(g20); g11_pf = pf_from_sympy(g11); ng20_pf = pf_neg(g20_pf)
def pf_nz(P): 
    return bool(P[0]) or any(not e.is_zero() for e in P[1]) or any(not e.is_zero() for e in P[2]) or not P[3].is_zero()
okG = all(not pf_nz(pf_sub(record[(k,(2,0))]['Gm'], g20_pf)) for k in range(2,13))
okG11 = all(not pf_nz(pf_sub(record[(k,(1,1))]['Gm'], ng20_pf if k==2 else g11_pf)) for k in range(2,13))
print("4. Gm_20==g20 (k=2..12):", okG, "  Gm_11==-g20 (k=2), ==g11 (k>=3):", okG11)
# 5. W_10^{(k)} == Rw2 for k=2..12 (direct, from reducer final kernel at (1,0))
Rw2   = 1/(2*v) + 1/(2*(v-(1+S2))) - 1/(v+(1+S2)) - 1/(2*(v+(S2-1)))
Rw2_pf = pf_from_sympy(Rw2)
okW=True
for k in range(2,13):
    W10 = record[(k,(1,0))]['W']
    # SimplePart = W minus HigherPart
    s_,h_,c_ = pf_split(W10)
    okW &= (not pf_nz(pf_sub(s_, Rw2_pf)))
print("5. SimplePart(W_10^(k)) == Rw2 for k=2..12:", okW)
# 6. eps-series identity
eps = sp.Symbol('eps')
a=1+S2; b=S2-1
R10   = (1-S2/2)/v - 2/(v+a)
DL1s = 1/(v-b) + 1/(v+a) - 1/v
DL2s = 1/(v-a) + 1/(v+b) - 1/v
def F(R): return sp.together(R.subs(v, 1/v)/v**2)
Gm20 = eps**2/(1-eps)*g20
Gm11 = -eps**2*g20 + eps**3/(1-eps)*g11
Wcas = -2*Gm20*DL1s - Gm11*DL2s - F(sp.together(Gm11*DL1s))
num = sp.expand(sp.together((1-eps)*Wcas))
c2 = sp.together(num.coeff(eps,2)); c3 = sp.together(num.coeff(eps,3))
rem = sp.together(num - c2*eps**2 - c3*eps**3)
ok6 = (rem==0) and sp.simplify(sp.radsimp(sp.together(c2-Rw2)))==0 and sp.simplify(sp.radsimp(sp.together(c3)))==0
print("6. (1-eps)(W_10 - eps R_10) = eps^2 Rw2 exactly:", ok6)
print()
print("ALL CHECKS:", all([ok, okI1, okI2, okI3, okI4, okD, okG, okG11, okW, ok6]))
'''
open('step_42.py','w').write(step)
import subprocess
r=subprocess.run(['<project>/.venv/bin/python','step_42.py'],capture_output=True,text=True)
print(r.stdout); print(r.stderr[-500:] if r.stderr else '')
