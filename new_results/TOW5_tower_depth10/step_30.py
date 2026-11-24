prove_S5 = r'''"""prove_S5.py -- PROOF for all (d,m) of
   S^{(5)}(t*)/S(t*) = N_5(d,m) / prod_{j=1..5}(2d-2j-1).
Extends prove_S34.py's method to l=5: normalized derivative = finite sum of explicit
rational functions of (d,m) over Q(i); Stirling tail P_5(x) = S(x,x-5) PROVED here by the
same recurrence + anchors. FIT-FREE: N_5 is read off the computed total (no fitted input).
Also: independent exact check vs direct symbolic differentiation of S(t) at integer (d,m).
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
# recurrence proof for all r <= 5 (symbolic identities)
for r in range(1,6):
    assert sp.expand(P[r].subs(x,x+1) - ((x+1-r)*P[r-1] + P[r])) == 0
# anchors: P_r vanishes at x = 0..r
for r in range(1,6):
    for a in range(0, r+1):
        assert P[r].subs(x,a) == 0

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
cand = sp.together(T*Pi(5))
num, den = sp.fraction(cand)
assert den == 1, f'denominator not 1: {den}'
N5 = sp.expand(num)
print(f'S^5/S identity PROVED: S^(5)(t*)/S(t*) = N5/Pi5 with N5 polynomial '
      f'(deg in d: {sp.degree(sp.Poly(N5, d))}, terms: {len(sp.Poly(N5, d, m).terms())})  ({round(time.time()-t0,1)}s)', flush=True)
R5 = sp.cancel(T)
# save
pickle.dump({'N5': str(N5)}, open(os.path.join(_HERE, 'N5_proved.pkl'), 'wb'))

# ---- independent exact check: direct differentiation of S(t) at integer (d,m) ----
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
ok = True
tested = 0
for dd in range(5, 13):
    for pp in range((dd+1)//2, dd+1):
        qq = dd - pp
        Sp = S_poly(dd, pp)
        S5 = sp.diff(Sp, 't', 5).subs('t', tstar)
        S0 = Sp.subs('t', tstar)
        direct = sp.simplify(S5/S0)
        pred = sp.simplify(R5.subs({d: dd, m: pp}))
        tested += 1
        if sp.simplify(direct - pred) != 0:
            ok = False; print('MISMATCH', dd, pp)
print(f'independent direct-differentiation check at {tested} integer keys d=5..12: {"ALL MATCH" if ok else "FAILED"} ({round(time.time()-t0,1)}s)', flush=True)
assert ok
print('VERDICT prove_S5: R_5 = N_5/Pi_5 PROVED for all (d,m); P_5 proved; direct check passed.', flush=True)
'''
with open('<project>/research/08_theory/TOW5_tower_depth10/prove_S5.py','w') as f:
    f.write(prove_S5)
print("written", len(prove_S5))
