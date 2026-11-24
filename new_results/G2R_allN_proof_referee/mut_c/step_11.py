"""step_11.py -- PROVED for all n: B4 as an identity of analytic functions on |c| < 1.

Method of b3_allN_proof.py, extended to k=4. mu_{4,n} is an explicit Q(n)-polynomial in
H^(r)_{n-1}, H^(r)_{2n-1} (r <= 4) (step_05), rewritten exactly as Q(n)-combinations of alternating
Z-sums Z_w(2n). On the HPL side [c^N] c^e H_w = h_w(N - e) with h_w(N) = sum coef eps^N N^{-k} Z_v(N-1),
and shifts Z_v(2n+j) reduce exactly to Z-sums at 2n. For N = 2n (even) and N = 2n+1 (odd) we check,
as exact identities of Q(n)-linear combinations of Z-sum words, for each of the five parts:
  (even)  sum_{terms} q h(2n - e)     == mu_{4,n} part       for all n >= 1,
  (odd)   sum_{terms} q h(2n+1 - e)  == 0                    for all n >= 1,
and verify the remaining monomial coefficients c^{-3}..c^3 directly in exact arithmetic.
The even-sector terms c^e (H_w + s H_{-w}) contribute h_w(N-e) + s h_{-w}(N-e).
Equality of canonical forms is sufficient (no independence statement needed).
Writes mu4_zsum.pkl; prints the PROVED verdict per part.
"""
import sys, json, re, pickle
from fractions import Fraction as Fr
import sympy as sp
sys.path.insert(0, "<project>/research/10_frontier/W1_corner_family/B_taylor_pi4")
sys.setrecursionlimit(10000)
import zsum
from zsum import K, nn
from hpl import coeffs_exact
import g2_fit as gf

G1 = "<project>/research/10_frontier/G1_zeta3_any_angle/"
sys.path.insert(0, '.')
from step_10 import load_terms, pi4_terms


def toK(s):
    e = sp.sympify(s)
    num, den = sp.fraction(sp.together(e))
    def polyK(p):
        P = sp.Poly(sp.expand(p), sp.Symbol('n'))
        out = K(0)
        for (k,), c in P.terms():
            c = sp.Rational(c)
            out += K(int(c.p))/K(int(c.q))*nn**k
        return out
    return polyK(num)/polyK(den)


def mu4_zsum():
    terms = pickle.load(open('mu4_formula.pkl', 'rb'))
    sub = {f'h{r}': zsum.harm_n_minus_1(r) for r in range(1, 5)}
    sub.update({f'g{r}': zsum.harm_2n_minus_1(r) for r in range(1, 5)})
    NAMES = ['h1', 'h2', 'h3', 'h4', 'g1', 'g2', 'g3', 'g4']
    parts = {'1': {}, 'z2': {}, 'z3': {}, 'z4': {}, 'z2^2': {}}
    for mon, cf in terms:
        E = {(): K(1)}
        for X, e in zip(NAMES, mon[:8]):
            for _ in range(e):
                E = zsum.mul(E, sub[X])
        e2, e3, e4 = mon[8], mon[9], mon[10]
        key = 'z2^2' if e2 == 2 else ('z2' if e2 else ('z3' if e3 else ('z4' if e4 else '1')))
        cfK = toK(cf)
        for w, c in E.items():
            zsum.add_into(parts[key], w, cfK*c)
    return parts


def hpl_side(terms, parity):
    E = {}
    for kind, w, e, q in terms:
        if kind in ('1', 'c', '1c'):
            continue
        qK = K(q.numerator)/K(q.denominator)
        if kind == 'hpl':
            X = zsum.hpl_coeff_at(w, parity - e)
            for ww, c in X.items():
                zsum.add_into(E, ww, qK*c)
        else:  # even-sector: h_w + sgn * h_{-w}
            sgn = 1 if (e + gf.nz(w)) % 2 == 0 else -1
            X = zsum.hpl_coeff_at(w, parity - e)
            X2 = zsum.hpl_coeff_at(gf.neg(w), parity - e)
            for ww, c in X.items():
                zsum.add_into(E, ww, qK*c)
            for ww, c in X2.items():
                zsum.add_into(E, ww, qK*c*sgn)
    return E


def low_order_check(terms, part, tab):
    """exact Taylor coefficients of the combination at c^N, N = -3..3, vs mu (mu_{4,0} := 0, odd := 0)."""
    ok = True
    for N in range(-3, 4):
        tot = Fr(0)
        for kind, w, e, q in terms:
            if kind in ('1', 'c'):
                if kind == '1' and N == 0: tot += q
                if kind == 'c' and N == 1: tot += q
                continue
            if kind == '1c':
                if N == e: tot += q
                continue
            idx = N - e
            if idx >= 1:
                if kind == 'hpl':
                    tot += q*coeffs_exact(w, idx + 2)[idx]
                else:
                    sgn = 1 if (e + gf.nz(w)) % 2 == 0 else -1
                    tot += q*(coeffs_exact(w, idx + 2)[idx] + sgn*coeffs_exact(gf.neg(w), idx + 2)[idx])
        target = Fr(0)
        if N >= 2 and N % 2 == 0:
            target = Fr(tab[str(N//2)].get(part, '0'))
        if tot != target:
            ok = False
            print(f'   low-order mismatch part {part} at c^{N}: {tot} vs {target}')
    return ok


if __name__ == '__main__':
    import time
    t0 = time.time()
    Z = mu4_zsum()
    print(f'mu_4 Z-sum canonical forms built ({time.time()-t0:.0f}s):',
          {p: len(E) for p, E in Z.items()}, flush=True)
    pickle.dump({p: [(w, str(c.as_expr())) for w, c in E.items()] for p, E in Z.items()},
                open('mu4_zsum.pkl', 'wb'))
    tab = json.load(open('flat_exact_N640_K4_formula.json'))['4']
    p4 = pi4_terms()
    PARTS = [
        ('1', '(0, 0, 0)', load_terms(['hpl_fit_k4_(0,0,0).json'])),
        ('z2', '(1, 0, 0)', load_terms(['hpl_fit_k4_(1,0,0).json'])),
        ('z3', '(0, 1, 0)', [(k, w, e, -q) for k, w, e, q in load_terms([G1 + 'hpl_fit_k3_(1,0,0).json'])]),
        ('z4', '(0, 0, 1)', p4['F1']),
        ('z2^2', '(2, 0, 0)', p4['F2']),
    ]
    allok = True
    for name, key, terms in PARTS:
        muE = {tuple(w): c for w, c in Z[name].items()}   # already K-elements
        ev = hpl_side(terms, 0)
        diff = zsum.lin((1, ev), (-1, muE))
        od = hpl_side(terms, 1)
        lo = low_order_check(terms, key, tab)
        print(f'part {name}: {len(terms)} terms; even words {len(ev)}, mu words {len(muE)}')
        print(f'   even identity [c^(2n)] (B4 - mu_4)|_{name} == 0 identically in n : {len(diff) == 0}')
        if diff:
            print('   residual:', list(diff.items())[:3])
        print(f'   odd identity  [c^(2n+1)] == 0 identically in n                : {len(od) == 0}')
        if od:
            print('   odd residual:', list(od.items())[:3])
        print(f'   low orders c^-3..c^3 exact                                   : {lo}')
        allok = allok and not diff and not od and lo
        sys.stdout.flush()
    print('B4 PROVED for all n (all five parts):', allok, f'({time.time()-t0:.0f}s)')
