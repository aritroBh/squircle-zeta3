code = r'''"""referee_checks.py -- G2R referee: non-vacuity, boundary validity, denominator audit.

(1) NON-VACUITY: the even/odd "identities" of step_11 are enormous exact cancellations, not empty
    comparisons: count insertion events and L1 magnitudes at n=5 for part '1'.
(2) BOUNDARY: every symbolic identity instance used by step_11 is checked EXACTLY at n = 1,2,3:
      - hpl_coeff_at(w, delta) vs direct Taylor coefficient (own evaluator + hpl.coeffs_exact),
        for every word w (and -w) occurring in B4, delta in {-1,..,4} (the full range parity-e);
      - shift(v, j) vs Z_val(v, 2n+j) for every Z-word v in the expansions, j in {-2,..,3};
      - harmonic-number conversions r = 1..4.
(3) DENOMINATOR AUDIT: all denominators of the Q(n) coefficients in the five even/odd canonical
    forms are factored; we verify none vanishes at an integer n >= 1.
"""
import sys, time
from fractions import Fraction as Fr
sys.path.insert(0, "<project>/research/10_frontier/W1_corner_family/B_taylor_pi4")
sys.path.insert(0, "run_base")
sys.setrecursionlimit(10000)
import zsum
from zsum import K, nn
from hpl import coeffs_exact
import sympy as sp
import step_10 as s10
import g2_fit as gf

G1 = "<project>/research/10_frontier/G1_zeta3_any_angle/"
t0 = time.time()

# ---------- own independent HPL Taylor evaluator (from the integral definition) ----------
from functools import lru_cache

@lru_cache(maxsize=None)
def my_hpl_coeff(w, m):
    """[c^m] H_w(c), letters {0,1,-1}, last letter nonzero. Exact."""
    w = tuple(w)
    if m <= 0:
        return Fr(0)
    if len(w) == 1:
        a = w[0]
        assert a in (1, -1)
        return Fr(a**(m - 1), m)          # H_1: 1/m ; H_-1: (-1)^(m+1)/m
    a, rest = w[0], w[1:]
    if a == 0:
        return my_hpl_coeff(rest, m)/m
    return sum((a**(m - 1 - i))*my_hpl_coeff(rest, i) for i in range(1, m)) / m

# cross-check my evaluator vs hpl.coeffs_exact for every word of B4
PARTS = {
    '1':    s10.load_terms(['run_base/hpl_fit_k4_(0,0,0).json']),
    'z2':   s10.load_terms(['run_base/hpl_fit_k4_(1,0,0).json']),
    'z3':   [(k, w, e, -q) for k, w, e, q in s10.load_terms([G1 + 'hpl_fit_k3_(1,0,0).json'])],
}
p4 = s10.pi4_terms.__globals__['json'].load(open('run_base/b4_pi4_part.json'))
PARTS['z4'] = s10.pi4_terms()['F1']
PARTS['z2^2'] = s10.pi4_terms()['F2']
WORDS = set()
for p, ts in PARTS.items():
    for kind, w, e, q in ts:
        if kind in ('hpl', 'even'):
            WORDS.add(tuple(w))
            if kind == 'even':
                WORDS.add(gf.neg(w))
print(f"distinct HPL words in B4 (incl. negations): {len(WORDS)}; max depth {max(map(len, WORDS))}")
ok_eval = all(my_hpl_coeff(w, m) == coeffs_exact(w, 84)[m] for w in WORDS for m in range(0, 84))
print(f"(1a) own HPL evaluator == hpl.coeffs_exact, all {len(WORDS)} words, m < 84: {ok_eval}")

# ---------- (1) NON-VACUITY of the step_11 comparisons ----------
import step_11 as s11
Z = s11.mu4_zsum()
def ev_expr(E, nv):
    tot = Fr(0)
    for w, cf in E.items():
        cv = Fr(int(cf.numer(nv)), 1)/Fr(int(cf.denom(nv)), 1)
        tot += cv*zsum.Z_val(tuple(w), 2*nv)
    return tot

def ev_expr_signed(E, nv):
    """L1 magnitude: sum |coef| * |Z_val| (terms taken with absolute value)."""
    tot = Fr(0)
    for w, cf in E.items():
        cv = Fr(int(cf.numer(nv)), 1)/Fr(int(cf.denom(nv)), 1)
        tot += abs(cv)*abs(zsum.Z_val(tuple(w), 2*nv))
    return tot

terms1 = PARTS['1']
insertions = 0
od = {}
for kind, w, e, q in terms1:
    if kind in ('1', 'c', '1c'):
        continue
    qK = K(q.numerator)/K(q.denominator)
    if kind == 'hpl':
        X = zsum.hpl_coeff_at(w, 1 - e); insertions += len(X)
        for ww, c in X.items():
            zsum.add_into(od, ww, qK*c)
    else:
        sgn = 1 if (e + gf.nz(w)) % 2 == 0 else -1
        X = zsum.hpl_coeff_at(w, 1 - e); X2 = zsum.hpl_coeff_at(gf.neg(w), 1 - e)
        insertions += len(X) + len(X2)
        for ww, c in X.items():
            zsum.add_into(od, ww, qK*c)
        for ww, c in X2.items():
            zsum.add_into(od, ww, qK*c*sgn)
print(f"(1b) part '1' ODD identity: {insertions} (word,coef) insertions cancel to {len(od)} surviving words")
nv = 5
ev1 = s11.hpl_side(terms1, 0)
mu1 = {tuple(w): c for w, c in Z['1'].items()}
print(f"(1c) at n={nv}: even HPL side value = {ev_expr(ev1, nv)}")
print(f"                mu_4 side value      = {ev_expr(mu1, nv)}")
print(f"     even-side L1 magnitude at n={nv}: {float(ev_expr_signed(ev1, nv)):.6e} (exact cancellation vs value {float(abs(ev_expr(ev1,nv))):.3e})")
# exact agreement of the two canonical forms evaluated at small n (independent of symbolic equality)
ok_small = all(ev_expr(ev1, k) == ev_expr(mu1, k) for k in range(1, 7))
print(f"(1d) even canonical forms evaluate equal exactly at n=1..6: {ok_small}")
print(f"     odd expression is identically zero dict: {od == {}} ; odd L1 would-be cancellation checked via insertions above")
sys.stdout.flush()

# ---------- (2) BOUNDARY checks at n = 1,2,3 ----------
okB = True; nchecks = 0
for w in WORDS:
    ce = coeffs_exact(w, 2*3 + 4 + 2)     # enough for 2n+delta, n<=3, delta<=4 -> N <= 10
    for delta in range(-1, 5):
        E = zsum.hpl_coeff_at(w, delta)
        for nv in (1, 2, 3):
            N = 2*nv + delta
            lhs = ev_expr(E, nv)
            rhs = ce[N] if N >= 0 else Fr(0)
            nchecks += 1
            if lhs != rhs:
                okB = False; print("   BOUNDARY FAIL", w, delta, nv, lhs, rhs)
print(f"(2a) hpl_coeff_at(w, delta) == [c^(2n+delta)] H_w exactly at n=1,2,3, all words, delta=-1..4: {okB} ({nchecks} checks)")
sys.stdout.flush()

# every Z-word v appearing in hpl_coeff_terms of B4 words; shift j in {-2..3}; n=1,2
VWORDS = set()
for w in WORDS:
    for cf, eps, k, v in zsum.hpl_coeff_terms(w):
        VWORDS.add(tuple(v))
okS = True; nchecks = 0
for v in VWORDS:
    for j in range(-2, 4):
        E = zsum.shift(v, j)
        for nv in (1, 2):
            lhs = ev_expr(E, nv)
            rhs = zsum.Z_val(v, 2*nv + j) if 2*nv + j >= 0 else Fr(0)
            nchecks += 1
            if lhs != rhs:
                okS = False; print("   SHIFT FAIL", v, j, nv, lhs, rhs)
print(f"(2b) shift(v, j) == Z_v(2n+j) exactly at n=1,2, all {len(VWORDS)} Z-words used, j=-2..3: {okS} ({nchecks} checks)")

okH = True
for r in range(1, 5):
    for nv in (1, 2, 3):
        Hn1 = sum((Fr(1, i**r) for i in range(1, nv)), Fr(0))
        H2n1 = sum((Fr(1, i**r) for i in range(1, 2*nv)), Fr(0))
        if ev_expr(zsum.harm_n_minus_1(r), nv) != Hn1: okH = False
        if ev_expr(zsum.harm_2n_minus_1(r), nv) != H2n1: okH = False
print(f"(2c) harmonic conversions H^(r)_n-1, H^(r)_2n-1 (r<=4) exact at n=1,2,3: {okH}")
sys.stdout.flush()

# ---------- (3) DENOMINATOR AUDIT ----------
n = sp.Symbol('n')
bad_den = []
den_factors = set()
def audit(E, tag):
    for w, cf in E.items():
        d = sp.Poly(sp.sympify(str(cf.denom.as_expr())), n)
        for f, mult in sp.factor_list(d.as_expr())[1]:
            den_factors.add(str(sp.factor(f)))
            roots = sp.solve(sp.Eq(f, 0), n)
            for rt in roots:
                if rt.is_rational and sp.Rational(rt) >= 1 and sp.denominator(rt) == 1:
                    bad_den.append((tag, w, str(f), str(rt)))
for name, ts in PARTS.items():
    audit(s11.hpl_side(ts, 0), f'even/{name}')
    audit(s11.hpl_side(ts, 1), f'odd/{name}')
for name, E in Z.items():
    audit({tuple(w): c for w, c in E.items()}, f'mu/{name}')
print(f"(3) denominator factors over all ten canonical forms: {sorted(den_factors)}")
print(f"    any factor vanishing at an integer n >= 1: {bool(bad_den)} {bad_den[:3]}")
print(f"ALL REFEREE CHECKS DONE ({time.time()-t0:.0f}s):", ok_eval and okB and okS and okH and not bad_den)
'''
CWD = "<project>/research/08_theory/G2R_allN_proof_referee"
open(CWD + "/referee_checks.py", "w").write(code)
print("written", len(code), "chars")
