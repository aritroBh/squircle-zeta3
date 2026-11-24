
import sys, time
sys.path.insert(0,'.')
from math import comb, factorial
from fractions import Fraction as Fr
from genintegrand import S_poly
import buildintegrand as bi
from fastred import *
IQ = QSI(0,0,1,0); SQ = QSI(0,1,0,0)
def sqrt_num_den(k, S):
    St, D = bi.Stilde(S)
    if not St: return None, None
    e = 2*k - 1 - D
    num = bi.qsi_pscale(St, bi.SQ2Q * QSI(Fr(-(8**k), 2*4**D)))
    dev = 0
    if e >= 0:
        for _ in range(e): num = bi.qsi_pmul(num, bi.VPOLY)
    else: dev = -e
    return num, dev
def qsi_peval(P, x):
    tot = ZERO
    for c in reversed(P): tot = tot*x + c
    return tot
def top_plus_i(k, S):
    num, dev = sqrt_num_den(k, S)
    if num is None: return ZERO
    ni = qsi_peval(num, IQ)
    d0 = (IQ**2 - AQ**2)*(IQ**2 - BQ**2) * (IQ**dev) * ((2*IQ)**(2*k-2))
    return ni / d0
def Kf(d):
    odd=1
    for m in range(1,2*d-2,2): odd*=m
    return (IQ**(d+1))*QSI(Fr((-1)**d*odd,2*factorial(d)))/SQ
def rho_cf(p,q):
    d=p+q
    return Kf(d)*((-ONE)**p)*QSI(comb(d,p)*(1+(1 if p>q else 0)))
bad6=[]; bad4=[]
t0=time.time()
for d in range(2,27):
    Ss = [S_poly(d,m) for m in range(d+1)]
    tops = [top_plus_i(d,S) for S in Ss]
    for i in range(d, (d-1)//2, -1):
        j = d-i
        tot = tops[i]
        if i!=j:
            tot = tops[i] + (-tops[j].conj_i())
        if tot != rho_cf(i,j): bad6.append((d,(i,j)))
    for q in range(0,(d-1)//2+1):
        p=d-q
        if p==q: continue
        if tops[p] != -tops[q].conj_i(): bad4.append((d,q))
    print(f"d={d} done [{time.time()-t0:.0f}s]", flush=True)
print("(C6) rho_pq closed form == top of ACTUAL level-d integrand, d=2..26:", "ALL MATCH" if not bad6 else bad6[:5])
print("(C4) folded top == unfolded top (via -conj), d=2..26:", "ALL MATCH" if not bad4 else bad4[:5])
