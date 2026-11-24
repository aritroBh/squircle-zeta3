code = '''
import sys, time
sys.path.insert(0,'.')
import sympy as sp
from math import comb, factorial
import buildintegrand
from genintegrand import S_poly, binom_half
from buildintegrand import sqrt_part_PF
from fastred import *
I = sp.I
t0=time.time()
# ---- (C4) fold doubling, d<=30 ----
def top_of(P, top):
    return P[1][top-1] if len(P[1])>=top else ZERO
bad=[]
for d in range(2,31):
    for q in range(0, (d-1)//2 + 1):
        p = d-q
        if p==q: continue
        Tp = top_of(sqrt_part_PF(d,S_poly(d,p)), 2*d-2)
        Tfold = top_of(pf_fold(sqrt_part_PF(d,S_poly(d,q))), 2*d-2)
        if Tp != Tfold: bad.append((d,q))
print("(C4) fold doubles the top coefficient, d=2..30:", "ALL MATCH" if not bad else bad[:5], f"[{time.time()-t0:.0f}s]")
# ---- (C5) binom(1/2,d) closed form ----
dd = sp.Symbol('d', positive=True, integer=True)
bh = sp.binomial(sp.Rational(1,2), dd)
closed = (-1)**(dd-1) * sp.factorial(2*dd-2) / (2**(dd-1)*sp.factorial(dd-1)) / (2**dd * sp.factorial(dd))
print("(C5) sympy symbolic: binom(1/2,d) / [(-1)^{d-1}(2d-3)!!/(2^d d!)] =", sp.simplify(sp.expand_func(bh)/closed))
ok5 = True
for d in range(1,101):
    odd=1
    for m in range(1,2*d-2,2): odd*=m
    ok5 &= (sp.Rational(binom_half(d)) == sp.Rational((-1)**(d-1)*odd, 2**d*factorial(d)))
print("(C5) exact numeric d=1..100:", ok5)
# ---- (C6) rho vs actual level-d integrand, d<=20 ----
from fractions import Fraction as Fr
IQ=QSI(0,0,1,0); SQ=QSI(0,1,0,0)
def Kf(d):
    odd=1
    for m in range(1,2*d-2,2): odd*=m
    return (IQ**(d+1))*QSI(Fr((-1)**d*odd,2*factorial(d)))/SQ
def rho_cf(p,q):
    d=p+q
    return Kf(d)*((-ONE)**p)*QSI(comb(d,p)*(1+(1 if p>q else 0)))
bad=[]
t0=time.time()
for d in range(2,21):
    integ = buildintegrand.build_integrand_pf(d)
    for (i,j),P in integ.items():
        tv = top_of(P, 2*d-2)
        if tv != rho_cf(i,j): bad.append((d,(i,j)))
print("(C6) rho_pq closed form == top of ACTUAL level-d integrand, d=2..20:", "ALL MATCH" if not bad else bad[:5], f"[{time.time()-t0:.0f}s]")
'''
open('rev_05b.py','w').write(code)
import sys, subprocess
work = '<project>/research/08_theory/WT2r_review'
r = subprocess.run([sys.executable, f'{work}/rev_05b.py'], capture_output=True, text=True, cwd=work, timeout=1700)
print(r.stdout); print(r.stderr[-800:] if r.returncode else "")
