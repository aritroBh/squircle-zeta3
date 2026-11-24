
import sys
sys.path.insert(0,'.')
import sympy as sp
from math import comb, factorial
import buildintegrand
from genintegrand import S_poly, w, binom_half, p_comp_1mt
from buildintegrand import sqrt_part_PF
from fastred import *
from fastred import v
I = sp.I; S2 = sp.sqrt(2)

# ---- (C3) analytic R~ == buildintegrand sqrt_part, d<=10 (extend step_41's d<=6) ----
t_ = (v**2+2*v-1)/(4*v)
z_ = (v**2+1)/(2*S2*v)
def Uj_poly(j,m,d):
    tt = sp.Symbol('T'); out=0
    for s in range(0,j+1):
        w1=w(s,m); w2=w(j-s,d-m)
        if w1==0 or w2==0: continue
        out += comb(j,s)*sp.Rational(w1.numerator,w1.denominator)*sp.Rational(w2.numerator,w2.denominator)*(1-tt)**(2*s)*tt**(2*(j-s))
    return sp.expand(out)
def Rtilde(d,m):
    tt=sp.Symbol('T'); ssum=0
    for j in range(1,d+1):
        bh=binom_half(j)
        U=Uj_poly(j,m,d)
        ssum += sp.Rational(bh.numerator,bh.denominator)*sp.expand(U)*z_**(2-2*j)
    ssum = sp.expand(ssum).subs(tt,t_)
    return sp.together(ssum/(2*S2*v*t_*(1-t_)))
bad=[]
for d in range(2,11):
    for m in range(0,d+1):
        a = sp.together(pf_to_sympy(sqrt_part_PF(d,S_poly(d,m))))
        if sp.simplify(sp.together(a-Rtilde(d,m)))!=0: bad.append((d,m))
print("(C3) R~(v) analytic formula == sqrt_part_PF, d=2..10 all m:", "ALL MATCH" if not bad else bad)

# ---- (C4) fold doubling: top of pf_fold(sqrt_part(S_q)) == top of sqrt_part(S_{d-q}), d<=30 ----
# (both at +i order 2d-2)
def top_of(P, top):
    return P[1][top-1] if len(P[1])>=top else ZERO
bad=[]
for d in range(2,31):
    for q in range(0, d//2 + (0 if d%2 else 1)):
        p = d-q
        if p==q: continue
        Tp = top_of(sqrt_part_PF(d,S_poly(d,p)), 2*d-2)
        Tfold = top_of(pf_fold(sqrt_part_PF(d,S_poly(d,q))), 2*d-2)
        if Tp != Tfold: bad.append((d,q))
print("(C4) fold doubles the top coefficient (fold top == unfolded top), d=2..30:", "ALL MATCH" if not bad else bad[:5])

# ---- (C5) binom(1/2,d) = (-1)^{d-1} (2d-3)!!/(2^d d!) : symbolic + exact ----
dd = sp.Symbol('d', positive=True, integer=True)
bh = sp.binomial(sp.Rational(1,2), dd)
# use (2d-3)!! = (2d-2)!/(2^{d-1}(d-1)!)
closed = (-1)**(dd-1) * sp.factorial(2*dd-2) / (2**(dd-1)*sp.factorial(dd-1)) / (2**dd * sp.factorial(dd))
print("(C5) sympy: binom(1/2,d)/closed_form simplifies to 1:", sp.simplify(sp.expand_func(bh)/closed))
ok5 = True
for d in range(1,101):
    odd=1
    for m in range(1,2*d-2,2): odd*=m
    ok5 &= (sp.Rational(binom_half(d)) == sp.Rational((-1)**(d-1)*odd, 2**d*factorial(d)))
print("(C5) exact numeric d=1..100:", ok5)

# ---- (C6) rho closed form vs ACTUAL LEVEL-d INTEGRAND top coefficients, d=2..20 ----
IQ=QSI(0,0,1,0); SQ=QSI(0,1,0,0)
def Kf(d):
    odd=1
    for m in range(1,2*d-2,2): odd*=m
    return (IQ**(d+1))*QSI(Fr((-1)**d*odd,2*factorial(d)))/SQ
def rho_cf(p,q):
    d=p+q
    return Kf(d)*((-ONE)**p)*QSI(comb(d,p)*(1+(1 if p>q else 0)))
bad=[]
for d in range(2,21):
    integ = buildintegrand.build_integrand_pf(d)
    for (i,j),P in integ.items():
        tv = top_of(P, 2*d-2)
        if tv != rho_cf(i,j): bad.append((d,(i,j)))
print("(C6) rho_pq closed form == top coeff of ACTUAL level-d integrand, d=2..20 all keys:",
      "ALL MATCH" if not bad else bad[:5])

# ---- (C7) only j=d contributes to the top of R~ (d<=15) ----
bad=[]
for d in range(2,16):
    for m in range(0,d+1):
        # full top:
        full = sqrt_part_PF(d,S_poly(d,m))
        tfull = top_of(full,2*d-2)
        # j=d-only evaluation: (1/(2 sqrt2)) binom(1/2,d) 2^{d-1} (-2i) U_d(tau)
        Ud = Uj_poly(d,m,d).subs(sp.Symbol('T'), (1+I)/2)
        pred = sp.simplify(sp.Rational(binom_half(d)) * 2**(d-1) * (-2*I) * Ud / (2*S2))
        # compare as QSI
        if qsi_from_sympy_real(sp.expand(sp.radsimp(pred.rewrite(sp.exp)))) if False else None: pass
        ps = sp.expand(pred)
        val = QSI(0)
        # pred in Q(i): extract
        pr = sp.simplify(ps.subs(I,0)); pi_ = sp.simplify((ps-pr)/I)
        val = QSI(Fr(int(pr.numerator),int(pr.denominator)),0,Fr(int(pi_.numerator),int(pi_.denominator)),0)
        if val != tfull: bad.append((d,m, val.sympy(), tfull.sympy()))
print("(C7) top of R~ equals j=d term alone, d=2..15 all m:", "ALL MATCH" if not bad else bad[:4])
