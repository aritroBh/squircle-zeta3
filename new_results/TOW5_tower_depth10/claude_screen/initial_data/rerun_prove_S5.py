"""prove_S5.py -- PROOF for all (d,m) of
   R_5 := S^{(5)}(t*)/S(t*) = N_5(d,m) / (3 * prod_{j=1..5}(2d-2j-1)).
Extends prove_S34.py's method to l=5: the normalized derivative is a finite sum of explicit
rational functions of (d,m) over Q(i) (TOW2 proof_all_depths.md Sec.1, proved method); the new
ingredient is the Stirling tail P_5(x) = S(x,x-5), PROVED below by the same recurrence + anchors.
FIT-FREE: N_5 is read off the computed total (gcd-reduced), not fitted. The identity
total(5) == N_5/(3 Pi_5) is proved symbolically (cancel'd difference == 0), then CONFIRMED by
direct exact differentiation of S(t) at 24 integer keys d=6..10 (independent of the derivation).
"""
import sympy as sp, time, pickle, os
d, m = sp.symbols('d m')
i = sp.I
_HERE = os.path.dirname(os.path.abspath(__file__))

def ff(x, c):
    r = sp.Integer(1)
    for u in range(c): r *= (x-u)
    return r

def tf(A):
    return (i/2)**(A//2) if A % 2 == 0 else (1+i)/2 * (i/2)**((A-1)//2)
def tf1(A):
    return (-i/2)**(A//2) if A % 2 == 0 else (1-i)/2 * (-i/2)**((A-1)//2)

x = sp.Symbol('x')
P = {0: sp.Integer(1), 1: x*(x-1)/2, 2: x*(x-1)*(x-2)*(3*x-5)/24,
     3: x*(x-1)*(x-2)**2*(x-3)**2/48,
     4: x*(x-1)*(x-2)*(x-3)*(x-4)*(15*x**3-150*x**2+485*x-502)/5760,
     5: x*(x-1)*(x-2)*(x-3)*(x-4)**2*(x-5)**2*(3*x**2-23*x+38)/11520}
# PROOF of the Stirling tails for all r <= 5:
#   (i) recurrence P_r(x+1) = (x+1-r) P_{r-1}(x) + P_r(x) as exact polynomial identities;
#   (ii) anchors P_r(a) = 0 for a = 0..r;
# => P_r(n) = S(n, n-r) for all integers n >= 0: both sides satisfy
#    a_{r}(n+1) = (n+1-r) a_{r-1}(n) + a_r(n)  (Stirling recurrence S(n+1,k) = k S(n,k) + S(n,k-1)
#    at k = n+1-r) and agree at n = 0 (induction on n, using P_{r-1}(n) = S(n,n-r+1)).
for r in range(1,6):
    assert sp.expand(P[r].subs(x,x+1) - ((x+1-r)*P[r-1] + P[r])) == 0
for r in range(1,6):
    for a in range(0, r+1):
        assert P[r].subs(x,a) == 0
print('P_5 (and P_r, r<=5) PROVED: recurrence + anchors hold', flush=True)

def normalized_summand(l, n, r, a, c):
    b = l - a
    jj = d - n; s = m - r
    t = sp.Integer(1)
    for rr in range(0, n):
        R = d - n + 1 + rr
        t *= -2*R/(2*R-3)
    t *= ff(m, r) * ff(d-m, n-r) / ff(d, n)
    t *= P[r].subs(x, m) / ff(m, r)
    t *= P[n-r].subs(x, d-m) / ff(d-m, n-r)
    t *= sp.binomial(l, a) * sp.Rational(2)**n * sp.binomial(n, a-n) * i**(2*n-a) * sp.factorial(a)
    A1 = 2*(r-n) - c; A2 = -2*r - b + c
    t *= sp.binomial(b, c) * ff(2*(jj-s), c) * (-1)**(b-c) * ff(2*s, b-c) * tf(A1) * tf1(A2)
    return sp.expand(t)

def Pi(l):
    out = sp.Integer(1)
    for jj in range(1, l+1): out *= (2*d-2*jj-1)
    return out

def total(l):
    tot = sp.Integer(0)
    for n in range(0, l+1):
        for r in range(0, n+1):
            for a in range(n, l+1):
                for c in range(0, l-a+1):
                    tot += normalized_summand(l, n, r, a, c)
    return sp.expand(tot)

t0 = time.time()
T = total(5)
print(f'l=5: total built ({round(time.time()-t0,1)}s)', flush=True)
R5 = sp.cancel(T)
num5, den5 = sp.fraction(R5)
Pi5 = Pi(5)
ratio = sp.cancel(den5 / Pi5)
assert ratio == 3, f'den5/Pi_5 = {ratio} (expected 3)'
N5 = sp.expand(num5)
print(f'R_5 = N_5/(3 Pi_5) PROVED symbolically: N_5 degree {sp.degree(N5,d)} in d, '
      f'{len(sp.Poly(N5,d,m).terms())} terms, lc {sp.LC(sp.Poly(N5,d))} d^10  ({round(time.time()-t0,1)}s)', flush=True)
# symbolic identity assertion (cancel'd difference is exactly 0):
assert sp.cancel(T - N5/(3*Pi5)) == 0
pickle.dump({'N5': str(N5), 'R5': str(R5)}, open(os.path.join(_HERE, 'N5_proved.pkl'), 'wb'))

# ---- independent exact check: direct differentiation of S(t) at integer (d,m), d=6..10 ----
from math import comb, factorial
def stir2(n, k):
    if n == k == 0: return 1
    if n <= 0 or k <= 0 or k > n: return 0
    return sum((-1)**(k-j) * comb(k, j) * j**n for j in range(k+1)) // factorial(k)
def S_poly(dd, pp):
    qq = dd - pp
    t = sp.Symbol('t')
    W = t**2 + (1-t)**2
    out = sp.Integer(0)
    for j in range(1, dd+1):
        Cj = sp.Rational(sp.binomial(sp.Rational(1,2), j) * factorial(j)) * 2**dd
        slo = max(0, j-qq); shi = min(j, pp)
        for s in range(slo, shi+1):
            Sps = stir2(pp, s); Sqj = stir2(qq, j-s)
            if Sps == 0 or Sqj == 0: continue
            out += Cj * Sps * Sqj * t**(2*(j-s)) * (1-t)**(2*s) * W**(dd-j)
    return sp.expand(out) / (factorial(pp)*factorial(qq))
tstar = (1+i)/2
ok = True; tested = 0
for dd in range(6, 11):
    for pp in range((dd+1)//2, dd+1):
        Sp = S_poly(dd, pp)
        direct = sp.simplify(sp.diff(Sp, 't', 5).subs('t', tstar) / Sp.subs('t', tstar))
        pred = sp.simplify(R5.subs({d: dd, m: pp}))
        tested += 1
        if sp.simplify(direct - pred) != 0:
            ok = False; print('MISMATCH', dd, pp)
print(f'independent direct-differentiation check, {tested} integer keys d=6..10: '
      f'{"ALL MATCH" if ok else "FAILED"} ({round(time.time()-t0,1)}s)', flush=True)
assert ok
print('VERDICT prove_S5: R_5 = N_5/(3 Pi_5) PROVED for all (d,m); P_5 proved; direct check passed.', flush=True)
