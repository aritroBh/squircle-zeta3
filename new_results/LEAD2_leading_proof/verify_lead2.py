"""LEAD2 verification script. Reruns every identity used in proof_lead.md:
  A. D-table extraction (d=1..9) == binomial x zigzag pattern (exact rationals, Newton differences)
  B. u-measures (m_odd, m_e, m_diag, m_oddlast), t = (1+u)/2, m_odd = 1/(1-u) - u m_e  (symbolic)
  C. even-d chain (d=2,4,6,8) and odd-d chain (d=1,3,5,7,9): key sum == collapsed forms == (d!/2)e_{d+1}
  D. closed forms T_d = d!(1-2^{-d})zeta(d+1), piece2 = T_d - (d!/2)e_{d+1}  (d=1..11 odd)
  E. Master Lemma: WW = E+E-, FF = (E-E-)/s, GG = 2pi/sin - (E+E-)/s; Beta sub; IBP; F-series; recursion
  F. full pipeline: exact kernels -> D-table -> quadrature -> c(w), w=2..10, vs the law (dps 100)
  G. E(y) Gamma identities, duplication, log-E expansion, zigzag signs, Hadamard GF form
Prints PASS/FAIL per item. Expected runtime: a few minutes.
"""
import os, sys, json, time, math
os.nice(10)
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import mpmath
mpmath.mp.dps = 70
from mpmath import mpf, mpc, sqrt, log, pi, zeta, gamma, sin, quad, nsum, inf, exp as mexp, bernoulli
from fractions import Fraction as Fr
from math import comb
t0 = time.time()
ok_all = True
def report(name, ok, resid=None):
    global ok_all
    ok_all &= bool(ok)
    print(f"[{'PASS' if ok else 'FAIL'}] {name}" + (f"   residual={mpmath.nstr(resid,5)}" if resid is not None else ""), flush=True)

from lead2_lib import (E_coeffs, Efun, zigzag, eps, load_kernels, extract_Dtop, delta_frac,
                       predicted_D, m_e, m_odd, m_oddlast, L1, L2)
E = E_coeffs(40)

# ---------- A. D-table ----------
kern = load_kernels(18)
Dtop = extract_Dtop(kern, 9)
bad = []
for (d, key), rec in Dtop.items():
    i, j = key
    de = delta_frac(d)
    for pole in ['-a', '0', '-b', 'a']:
        got = rec[pole]['Dtop'] if pole in rec else Fr(0)
        if got != predicted_D(d, i, j, pole, de): bad.append((d, key, pole))
poly_ok = all(pr['poly_ok'] for rec in Dtop.values() for pr in rec.values())
report("A: D-table (d=1..9, all keys/poles) == binomial x zigzag; kernel degrees == d (d=9: determined fit)", not bad and poly_ok)

# ---------- B. u-measures (symbolic) ----------
import sympy as sp
u = sp.symbols('u'); S2 = sp.sqrt(2); S = sp.sqrt(1+u**2)
av = 1+S2; bv = S2-1; v = u+S; dvdu = sp.diff(v,u)
def sym0(expr): return sp.radsimp(sp.simplify(expr)) == 0
m_e_s    = (1-S2/S)/(1-u**2)
m_odd_s  = (1+S2*u/S)/(1-u**2)
checks = [
  ('odd',  (1/(v+av)-1/(v-av))*dvdu, m_odd_s),
  ('even', (1/(v+av)-1/(v+bv))*dvdu, m_e_s),
  ('diag', (1/(v+av)-1/(2*v))*dvdu, (1-u)/2*m_e_s),
  ('oddlast', (1/(v+av)-1/(2*(v-av))+1/(2*(v+bv))-1/(2*v))*dvdu, m_odd_s-1/(2*(1-u))),
]
okB = all(sym0(meas-tgt) for _, meas, tgt in checks)
okB &= sym0(m_odd_s - (1/(1-u) - u*m_e_s))
okB &= sym0((v-bv)*(v+av)/(4*v) - (1+u)/2)
okB &= sym0((av-v)*(v+bv)/(4*v) - (1-u)/2)
okB &= sp.simplify(v*(u-S)) == -1
report("B: u-measures, t=(1+u)/2, 1-t=(1-u)/2, m_odd=1/(1-u)-u m_e, deck v+*v-=-1 (symbolic)", okB)

# ---------- helpers for C/D/E ----------
def hfun(w): return ((1+w)/sqrt(1+w**2)-1)/w
def qfun(w): return (1+w)/(1+w**2)**mpf('1.5')
def g1fun(w): return (1-w)/(1+w)*hfun(w)
def Mint(i,j,m): return quad(lambda uu: L1(uu)**i*L2(uu)**j*m(uu), [-1,0,1])
def Aint(p,q):
    return (-1)**p*quad(lambda r: r**p*log(1-mexp(-r))**q/(mexp(r)-1), [0,1,mpmath.inf])
def T_odd(d):
    mid=(d-1)//2
    return sum((-1)**(j+1)*comb(d,j)*Aint(d-j,j) for j in range(mid+1)) \
           - mpf('0.5')*(-1)**(mid+1)*comb(d,mid)*Aint(mid+1,mid)
def piece2(d):
    return mpf('0.5')*quad(lambda w: (-log(w))**d*(1-w)/(1+w)*hfun(w), [0,1])

# ---------- C. even-d and odd-d chains ----------
for d in [2,4,6,8,10]:
    tot = mpf(0)
    for j in range(d//2+1):
        m = m_e if j < d//2 else (lambda uu: (1-uu)/2*m_e(uu))
        tot += (-1)**(j+1)*comb(d,j)*Mint(d-j,j,m)
    coll = -mpf('0.5')*quad(lambda uu: (L1(uu)-L2(uu))**d*m_e(uu), [-1,0,1])
    tgt = mpmath.factorial(d)/2*E[d+1]
    report(f"C: even d={d}: key sum == -(1/2)int log^d m_e == (d!/2)e_{{d+1}}",
           abs(tot-coll) < mpf(10)**(-48) and abs(coll-tgt) < mpf(10)**(-48), max(abs(tot-coll),abs(coll-tgt)))
for d in [1,3,5,7,9,11]:
    mid=(d-1)//2; tot = mpf(0)
    for j in range(mid+1):
        m = m_odd if j < mid else m_oddlast
        tot += (-1)**(j+1)*comb(d,j)*Mint(d-j,j,m)
    T = T_odd(d); p2 = piece2(d); tgt = mpmath.factorial(d)/2*E[d+1]
    report(f"C: odd d={d}: key sum == T_d - piece2 == (d!/2)e_{{d+1}}",
           abs(T-p2-tot) < mpf(10)**(-48) and abs(tot-tgt) < mpf(10)**(-48), max(abs(T-p2-tot),abs(tot-tgt)))

# ---------- D. closed forms ----------
for d in [1,3,5,7,9,11]:
    T = T_odd(d); p2 = piece2(d)
    f1 = mpmath.factorial(d)*(1-mpf(2)**(-d))*zeta(d+1)
    f2 = f1 - mpmath.factorial(d)/2*E[d+1]
    report(f"D: odd d={d}: T_d == d!(1-2^-d)zeta(d+1); piece2 == T_d - (d!/2)e_{{d+1}}",
           abs(T-f1) < mpf(10)**(-45) and abs(p2-f2) < mpf(10)**(-45), max(abs(T-f1),abs(p2-f2)))

# ---------- E. Master Lemma + Beta + IBP + F-series + recursion ----------
mpmath.mp.dps = 60
pts = [mpf('0.3'), mpf('-0.4'), mpf('0.73'), mpc(0,'0.2')]
for s in pts:
    W = quad(lambda w: w**(-s)*qfun(w), [0,1,mpmath.inf]); rhsW = Efun(s)+Efun(-s)
    FF = quad(lambda w: w**(-s)*hfun(w), [0,1,mpmath.inf]); rhsF = (Efun(s)-Efun(-s))/s
    G = quad(lambda w: w**(-s)*g1fun(w), [0,1,mpmath.inf]); rhsG = 2*pi/sin(pi*s)-(Efun(s)+Efun(-s))/s
    tol = mpf(10)**(-10) if abs(s) > mpf('0.7') else mpf(10)**(-35)
    report(f"E: Master Lemma at s={mpmath.nstr(s,6)}",
           abs(W-rhsW)<tol and abs(FF-rhsF)<tol and abs(G-rhsG)<tol,
           max(abs(W-rhsW),abs(FF-rhsF),abs(G-rhsG)))
for s in pts:
    bsub = quad(lambda w: (w**(1+s)+w**(-s))/(1+w**2)**mpf('1.5'), [0,1])
    F1 = quad(lambda w: w**(-s)*hfun(w), [0,1]); F2 = quad(lambda w: w**(s)*hfun(w), [0,1])
    tol = mpf(10)**(-10) if abs(s) > mpf('0.7') else mpf(10)**(-35)
    report(f"E: Beta sub + IBP at s={mpmath.nstr(s,6)}",
           abs(bsub-Efun(s))<tol and abs((Efun(s)-Efun(-s))-s*(F1+F2))<tol,
           max(abs(bsub-Efun(s)), abs((Efun(s)-Efun(-s))-s*(F1+F2))))
def cn(n): return mpf(comb(2*int(n),int(n)))/4**int(n)
def F_series(y):
    return nsum(lambda n: (-1)**int(n)*cn(n)/(2*n-y), [1, inf]) + \
           nsum(lambda n: (-1)**int(n)*cn(n)/(2*n+1-y), [0, inf])
for y in [mpf('0.3'), mpf('-0.45')]:
    fq = quad(lambda w: w**(-y)*hfun(w), [0,1])
    report(f"E: F partial-fraction series (n>=0 second sum) at y={y}", abs(F_series(y)-fq) < mpf(10)**(-40), abs(F_series(y)-fq))
worst = mpf(0)
for p in range(1,6):
    for q in range(0,5):
        worst = max(worst, abs(Aint(p,q)/(mpmath.factorial(p)*mpmath.factorial(q))
                                - Aint(q+1,p-1)/(mpmath.factorial(q+1)*mpmath.factorial(p-1))))
report("E: recursion aa_{p,q} = aa_{q+1,p-1}, p=1..5 q=0..4", worst < mpf(10)**(-40), worst)
# B(1+z,-z) = -pi/sin(pi z) and [z^d] pi/sin = 2(1-2^{-d})zeta(d+1) for odd d
z0 = mpf('0.37')
lhs = mpmath.beta(1+z0, -z0); rhs = -pi/sin(pi*z0)
report("E: B(1+z,-z) == -pi/sin(pi z) at z=0.37", abs(lhs-rhs) < mpf(10)**(-40), abs(lhs-rhs))
for d in [1,3,5]:
    n=(d+1)//2
    ser = (-1)**(n+1)*(2**(2*n)-2)*bernoulli(2*n)*pi**(2*n)/mpmath.factorial(2*n)
    report(f"E: [z^{d}] pi/sin(pi z) == 2(1-2^-{d}) zeta({d+1})", abs(ser-2*(1-mpf(2)**(-d))*zeta(d+1)) < mpf(10)**(-50), ser-2*(1-mpf(2)**(-d))*zeta(d+1))

# ---------- F. full pipeline kernels -> c(w) ----------
mpmath.mp.dps = 100
E = E_coeffs(40)   # recompute at dps 100
av2 = 1+sqrt(2); bv2 = sqrt(2)-1
POLES = {'-a':-av2,'0':mpf(0),'-b':-bv2,'a':av2}
def tfun(vv): return (vv-bv2)*(vv+av2)/(4*vv)
def ufun(vv): return (av2-vv)*(vv+bv2)/(4*vv)
jobs = set()
for (d,key),rec in Dtop.items():
    i,j=key
    for pole,pr in rec.items():
        if pr['Dtop']!=0: jobs.add((i,j,pole))
N = {}
for (i,j,pole) in sorted(jobs):
    r = POLES[pole]
    N[(i,j,pole)] = quad(lambda vv: log(tfun(vv))**i*log(ufun(vv))**j/(vv-r), [bv2, av2])
cw = {}
for d in range(1,10):
    tot = mpf(0)
    for (dd,key),rec in Dtop.items():
        if dd!=d: continue
        i,j=key
        for pole,pr in rec.items():
            cf = pr['Dtop']
            if cf: tot += mpf(cf.numerator)/cf.denominator*N[(i,j,pole)]
    cw[d+1] = tot/mpmath.factorial(d-1)
A = zigzag(11)
worst = mpf(0)
for w in range(2,11):
    m = w-2
    law = eps(m)*A[m+1]/(2*mpmath.factorial(m))*E[w]
    worst = max(worst, abs(cw[w]-law))
report("F: kernels -> D-table -> c(w) == law, w=2..10 (dps 100)", worst < mpf(10)**(-85), worst)

# ---------- G. Gamma identities + zigzag + Hadamard ----------
mpmath.mp.dps = 50
x = sp.symbols('x')
ser = sp.series((1+sp.sinh(x))/sp.cosh(x), x, 0, 22).removeO()
okG = all(sp.Rational(ser.coeff(x,d)) == sp.Rational(([1,-1,-1,1][(d-1)%4])*zigzag(d+1)[d], sp.factorial(d)) for d in range(1,21))
ser2 = sp.series((1-sp.sinh(x))/sp.cosh(x)**2, x, 0, 22).removeO()
A = zigzag(22)
okG &= all(sp.Rational(ser2.coeff(x,m)) == sp.Rational(([1,-1,-1,1][m%4])*A[m+1], sp.factorial(m)) for m in range(0,20))
report("G: delta_d = eps_{d-1}A_d/d! and eps_mA_{m+1}/m! = [z^m](1-sinh)/cosh^2, d,m <= 20 (symbolic)", okG)
y = mpf('0.31')
lhs = Efun(y); rhs = mpf(2)**y*gamma(1-y)*gamma(1+y/2)/gamma(1-y/2)
report("G: duplication E(y) = 2^y Gamma(1-y)Gamma(1+y/2)/Gamma(1-y/2) at y=0.31", abs(lhs-rhs) < mpf(10)**(-45), lhs-rhs)
NN=40
y = mpf('0.1')
logE = [mpf(0)]*(NN+1); logE[1] = log(2)
for n in range(2,NN+1): logE[n] = (zeta(n) if n%2==0 else (1-mpf(2)**(1-n))*zeta(n))/n
sersum = sum(logE[n]*y**n for n in range(NN+1))
report("G: log E = y log2 + sum s_n y^n/n at y=0.1 (truncation O(1e-41))", abs(log(Efun(y))-sersum) < mpf(10)**(-30), log(Efun(y))-sersum)
y = mpf('0.31')
# Hadamard GF form vs law coefficients: c(w) u^w summed = 1/2 [u^2(1-sinh u)/cosh^2 u ⊙ E(u)]
A = zigzag(14)
uu = mpf('0.21')
lhs = sum(eps(w-2)*A[w-1]/(2*mpmath.factorial(w-2))*E[w]*uu**w for w in range(2,14))
co1 = sp.series((1-sp.sinh(x))/sp.cosh(x)**2, x, 0, 14).removeO()
rhs = mpf('0.5')*sum(mpf(str(sp.N(co1.coeff(x,m),40)))*mpf(str(sp.N(E[m+2],40)))*uu**(m+2) for m in range(0,12))
report("G: sum c(w)u^w == (1/2)[u^2(1-sinh)/cosh^2 ⊙ E](u) at u=0.21 (law coefficients)", abs(lhs-rhs) < mpf(10)**(-35), lhs-rhs)

print(f"\nALL: {'PASS' if ok_all else 'FAIL'}   ({time.time()-t0:.0f}s)")
sys.exit(0 if ok_all else 1)
