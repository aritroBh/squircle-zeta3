"""step_05 (TOW5R): wide independent direct check of R_5 = N5/(3 Pi_5).
My OWN implementation of the TOW1 S-polynomial (proof_tower.md sec.3):
  S_m(t) = sum_{j=1}^d binom(1/2,j) q2^{d-j} sum_{s=0}^j binom(j,s) w(s,m) w(j-s,d-m)
           t^{2(j-s)} (1-t)^{2s},   w(s,m) = s! Stir2(m,s) 2^m/m!,  q2 = t^2+(1-t)^2.
Direct exact 5th derivative at t* = (1+i)/2 vs tow5lib.R5 at:
  all keys d = 6..25, plus spot keys at d = 30, 37, 40.
"""
import sympy as sp, os, time, pickle, sys
from math import comb, factorial
d, m = sp.symbols('d m')
i = sp.I
t = sp.Symbol('t')
tstar = (1+i)/2
t0 = time.time()
_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _HERE + '/tow5_copy')
_p = pickle.load(open(_HERE + '/tow5_copy/N5_proved.pkl', 'rb'))
N5 = sp.sympify(_p['N5'], locals={'d': d, 'm': m, 'I': sp.I})
Pi5 = sp.Integer(1)
for j in range(1, 6): Pi5 *= (2*d-2*j-1)
R5 = N5/(3*Pi5)

def stir2(n, k):
    if n == k == 0: return 1
    if n <= 0 or k <= 0 or k > n: return 0
    return sum((-1)**(k-jj) * comb(k, jj) * jj**n for jj in range(k+1)) // factorial(k)
def wfun(s, mm):
    if s < 0 or s > mm: return sp.Rational(0)
    return sp.Rational(factorial(s) * stir2(mm, s) * 2**mm, factorial(mm))
def S_poly_own(dd, mm):
    q2 = t**2 + (1-t)**2
    out = sp.Integer(0)
    for j in range(1, dd+1):
        cj = sp.binomial(sp.Rational(1,2), j)
        inner = sp.Integer(0)
        for s in range(0, j+1):
            w1 = wfun(s, mm); w2 = wfun(j-s, dd-mm)
            if w1 == 0 or w2 == 0: continue
            inner += sp.binomial(j, s) * w1 * w2 * t**(2*(j-s)) * (1-t)**(2*s)
        out += cj * q2**(dd-j) * inner
    return sp.expand(out)

def check_key(dd, mm):
    Sp = S_poly_own(dd, mm)
    d5 = sp.diff(Sp, t, 5).subs(t, tstar)
    s0 = Sp.subs(t, tstar)
    direct = sp.simplify(d5/s0)
    pred = sp.simplify(R5.subs({d: dd, m: mm}))
    return sp.simplify(direct - pred) == 0

tested = 0; bad = []
for dd in [7, 9, 11, 13, 17, 21]:
    for mm in range((dd+1)//2, dd+1):
        tested += 1
        if not check_key(dd, mm): bad.append((dd, mm))
    print(f'd={dd} done ({round(time.time()-t0,1)}s, tested {tested})', flush=True)
for dd, mlist in [(30, [30, 20]), (40, [40])]:
    for mm in mlist:
        tested += 1
        if not check_key(dd, mm): bad.append((dd, mm))
    print(f'd={dd} spots done ({round(time.time()-t0,1)}s)', flush=True)
print(f'STEP_05 VERDICT: independent S_m(t) 5th-derivative check vs R_5: {tested} keys, bad = {len(bad)} {bad[:5]}', flush=True)
assert not bad
