
import sys, time
sys.path.insert(0,'.')
import sympy as sp
d = sp.Symbol('d', positive=True, integer=True)
z = sp.Symbol('z')   # fresh
I = sp.I; S2 = sp.sqrt(2)

# (C5) ratio recurrence proof, gamma reduced manually
b = sp.binomial(sp.Rational(1,2), d)
b1 = sp.binomial(sp.Rational(1,2), d+1)
ratio = sp.gammasimp(b1/b)   # gamma(3/2-d)/((d+1) gamma(1/2-d))
ratio2 = sp.simplify(ratio.subs(sp.gamma(sp.Rational(3,2)-d), (sp.Rational(1,2)-d)*sp.gamma(sp.Rational(1,2)-d)))
def df(n): return sp.factorial(2*n-2)/(2**(n-1)*sp.factorial(n-1))
f = lambda n: (-1)**(n-1)*df(n)/(2**n*sp.factorial(n))
ratio_rhs = sp.simplify(f(d+1)/f(d))
print("(C5) binom(1/2,d+1)/binom(1/2,d) =", ratio2, "; f(d+1)/f(d) =", ratio_rhs,
      "; equal:", sp.simplify(ratio2-ratio_rhs)==0,
      "; base d=1: binom(1/2,1)=", sp.binomial(sp.Rational(1,2),1), "= f(1)=", f(1))
print("     => binom(1/2,d) = (-1)^{d-1} (2d-3)!!/(2^d d!) for all d>=1 (induction)")

# (C7) only j=d contributes to the top of R~, exact check d<=15 via fast evaluator
from math import comb, factorial
from fractions import Fraction as Fr
from genintegrand import S_poly, w, binom_half
import buildintegrand as bi
from fastred import *
IQ = QSI(0,0,1,0)
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
def Ud_at_tau(dd, m):
    # U_d^{(m)}(tau) via w's, tau=(1+i)/2: tau^2=i/2,(1-tau)^2=-i/2
    out = QSI(0)
    t2 = QSI(0,0,Fr(1,2),0); o2 = QSI(0,0,Fr(-1,2),0)   # i/2, -i/2
    for s in range(dd+1):
        w1=w(s,m); w2=w(dd-s,dd-m)
        if w1==0 or w2==0: continue
        out = out + QSI(comb(dd,s)*w1*w2)*(o2**s)*(t2**(dd-s))
    return out
bad=[]
for dd in range(2,16):
    for m in range(dd+1):
        tfull = top_plus_i(dd, S_poly(dd,m))
        Ud = Ud_at_tau(dd,m)
        pred = QSI(Fr(binom_half(dd))*Fr(2**(dd-1))) * QSI(0,0,-2,0) * Ud / QSI(0,2,0,0)  # /(2 sqrt2)
        if pred != tfull: bad.append((dd,m))
print("(C7) top of sqrt part == (1/(2sqrt2)) binom(1/2,d) 2^{d-1} (-2i) U_d(tau) [j=d term only], d=2..15:",
      "ALL MATCH" if not bad else bad[:5])
