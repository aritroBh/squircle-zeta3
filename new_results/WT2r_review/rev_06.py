
import sys, time
sys.path.insert(0,'.')
from math import comb, factorial
from fractions import Fraction as Fr
from genintegrand import S_poly
from fastred import *

# Replicate sqrt_part_PF's num/den construction (from buildintegrand source) WITHOUT pf_from_poly
import buildintegrand as bi
def sqrt_num_den(k, S):
    St, D = bi.Stilde(S)
    if not St: return None, None
    e = 2*k - 1 - D
    num = bi.qsi_pscale(St, bi.SQ2Q * QSI(Fr(-(8**k), 2*4**D)))
    den_extra_v = 0
    if e >= 0:
        for _ in range(e): num = bi.qsi_pmul(num, bi.VPOLY)
    else:
        den_extra_v = -e
    return num, den_extra_v

def qsi_peval(P, x):
    tot = ZERO
    for c in reversed(P): tot = tot*x + c
    return tot

IQ = QSI(0,0,1,0)
SQ = QSI(0,1,0,0)
def top_plus_i(k, S):
    """top coefficient at +i of order 2k-2 of the sqrt part for coefficient S."""
    num, dev = sqrt_num_den(k, S)
    if num is None: return ZERO
    # den = (v^2-a^2)(v^2-b^2) v^{dev} (v^2+1)^{2k-2}
    # top at +i = num(i) / [ (2i)^{2k-2} * (i^2-a^2)(i^2-b^2) * i^{dev} ]
    ni = qsi_peval(num, IQ)
    d0 = (IQ**2 - AQ**2)*(IQ**2 - BQ**2) * (IQ**dev) * ((2*IQ)**(2*k-2))
    return ni / d0

# sanity: agree with sqrt_part_PF top for small d
ok = True
for d in range(2,9):
    for m in range(d+1):
        P = bi.sqrt_part_PF(d, S_poly(d,m))
        t1 = P[1][2*d-3] if len(P[1])>=2*d-2 else ZERO
        t2 = top_plus_i(d, S_poly(d,m))
        ok &= (t1==t2)
print("direct top-coefficient evaluator == pf_from_poly top, d=2..8:", ok)

# (C6) rho closed form vs actual integrand, now d<=40
def Kf(d):
    odd=1
    for m in range(1,2*d-2,2): odd*=m
    return (IQ**(d+1))*QSI(Fr((-1)**d*odd,2*factorial(d)))/SQ
def rho_cf(p,q):
    d=p+q
    return Kf(d)*((-ONE)**p)*QSI(comb(d,p)*(1+(1 if p>q else 0)))
bad=[]; t0=time.time()
for d in range(2,41):
    Ss = [S_poly(d,m) for m in range(d+1)]
    for i in range(d, (d-1)//2, -1):
        j = d-i
        Tp = top_plus_i(d, Ss[i])
        tot = Tp
        if i!=j:
            Tfold = -top_plus_i(d, Ss[j]).conj_i()   # folded piece top = -conj of (j,i) top (reality)
            tot = Tp + Tfold
        # AB parts: no +-i poles -> no contribution
        if tot != rho_cf(i,j): bad.append((d,(i,j)))
print("(C6) rho_pq closed form == top of ACTUAL level-d integrand, d=2..40:", "ALL MATCH" if not bad else bad[:5], f"[{time.time()-t0:.0f}s]")

# (C4) fold doubling == -conj relation, d<=40 (already used above; state separately)
bad=[]
for d in range(2,41):
    Ss = [S_poly(d,m) for m in range(d+1)]
    for q in range(0,(d-1)//2+1):
        p=d-q
        if p==q: continue
        if top_plus_i(d,Ss[p]) != -top_plus_i(d,Ss[q]).conj_i(): bad.append((d,q))
print("(C4) top of folded piece == top of unfolded piece (=-conj relation), d=2..40:", "ALL MATCH" if not bad else bad[:5])
