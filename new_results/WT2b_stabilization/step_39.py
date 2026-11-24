import sys, pickle
sys.path.insert(0,'.')
from fastred import *
record = pickle.load(open('record_W2b.pkl','rb'))
tau1 = QSI(-1,0,-1,0)/QSI(2); tau2 = QSI(-1,0,1,0)/QSI(2)
def h(k,key):
    i,j=key; d=i+j
    rec=record.get((k,key))
    if rec is None: return None
    H=rec['H'][1]; top=2*d-2
    return H[top-1] if len(H)>=top else ZERO
bad=[]
for k in range(3,13):
    for i in range(1,k+1):
        for j in range(0,i+1):
            d=i+j
            if d<2 or k==d: continue
            hv = h(k,(i,j))
            if hv is None: continue
            s = ZERO
            hA = h(k,(i+1,j))
            if hA is not None: s = s + hA*QSI(i+1)*tau1
            if j==i-1:                      # adjacent: B and F share key (i,i)
                hd = h(k,(i,i))
                if hd is not None: s = s + (hd - hd.conj_i())*QSI(i)*tau2
            elif i>j:
                hB = h(k,(i,j+1))
                if hB is not None: s = s + hB*QSI(j+1)*tau2
            s = s/QSI(2*d-1)
            if hv != s: bad.append((k,(i,j),hv.sympy(),s.sympy()))
print("(V3a fixed) top-order recursion vs data (k=3..12):", "ALL MATCH" if not bad else bad[:6])

# (V3b) full top-order STABILIZATION as predicted by the theorem:
# off-diagonal: h_pq(k)=rho_pq for all k>=d ; diagonal: h_pp(k)=rho_pp (k=2p), hat(k>=2p+1)
from math import comb, factorial
IQ=QSI(0,0,1,0); SQ=QSI(0,1,0,0)
def K(d):
    odd=1
    for m in range(1,2*d-2,2): odd*=m
    return (IQ**(d+1))*QSI(Fr((-1)**d*odd,2*factorial(d)))/SQ
from fractions import Fraction as Fr
def rho_cf(p,q):
    d=p+q
    return K(d)*((-ONE)**p)*QSI(comb(d,p)*(1+(1 if p>q else 0)))
def h_hat(p):
    return rho_cf(p+1,p)*QSI(p+1)*tau1/QSI(4*p-1)
viol=[]
for k in range(2,13):
    for i in range(1,k+1):
        for j in range(0,i+1):
            d=i+j
            if d<2 or d>k: continue
            hv=h(k,(i,j))
            if hv is None: continue
            if i>j:
                pred = rho_cf(i,j)
            else:
                pred = rho_cf(i,i) if k==2*i else h_hat(i)
            if hv!=pred: viol.append((k,(i,j),hv.sympy(),pred.sympy()))
print("(V3b) top-order stabilization (all keys, k=2..12):", "ALL AS PREDICTED" if not viol else viol[:8])
