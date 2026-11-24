"""b3_allN_proof.py -- PROVED (exact symbolic-n identity of alternating nested sums).

Claim (G1's B3, stated there as EXACT through c^640 + NUMERICAL):
   M_3(alpha) = R_3(c) + z2 Psi(c) + z3 (-Phi(c)),     c = cos(alpha/2), |c| < 1,
with R_3, Psi, -Phi the HPL combinations stored in G1's hpl_fit_k3_{(0,0,0),(1,0,0),(0,1,0)}.json.

Proof.  By Theorem A (mu3_symbolic.py) M_3 = sum_{n>=1} mu_{3,n} c^{2n} with mu_{3,n} = r_n + q_n z2 + p_n z3,
each part an explicit Q(n)-combination of products of H^(r)_{n-1}, H^(r)_{2n-1}, rewritten exactly as a Q(n)-
combination of alternating Z-sums Z_w(2n) (mu_zsum.pkl).  On the other side [c^N] c^e H_w(c) = h_w(N - e), and
h_w(N) = sum coef eps^N N^{-k} Z_v(N-1) exactly (zsum.hpl_coeff_terms).  For N = 2n (resp. 2n+1) the signs
eps^N are constants and the shifts Z_v(2n + j) reduce exactly to Z-sums at 2n.  We check, as identities of
Q(n)-linear combinations of Z-sum words (an exact, finite computation), that
   (even)  sum_{e,w} q_{e,w} h_w(2n - e)   ==  part of mu_{3,n}      for all n >= 1,
   (odd)   sum_{e,w} q_{e,w} h_w(2n+1 - e) ==  0                      for all n >= 1,
for each of the three parts separately, and verify the finitely many remaining coefficients
c^{-2}, c^{-1}, c^0, c^1 directly in exact arithmetic.  Equality of the canonical forms is sufficient
(no independence statement is needed).  This proves B3 for all n, hence as an identity of analytic functions
on |c| < 1 (both sides converge there), and by continuity on [0,1).
Also prints the Taylor coefficients of the weight-4 component W_4(c) in the same canonical form.
"""
import sys, json, re, pickle
from fractions import Fraction as Fr
import sympy as sp
sys.path.insert(0, '.')
sys.setrecursionlimit(10000)
import zsum
from zsum import K, nn
from hpl import coeffs_exact

G1 = 'research/10_frontier/G1_zeta3_any_angle/'
FILES = {'1': 'hpl_fit_k3_(0,0,0).json', 'z2': 'hpl_fit_k3_(1,0,0).json', 'z3': 'hpl_fit_k3_(0,1,0).json'}


def load(name):
    out = []
    for nm, v in json.load(open(G1 + name)):
        if nm == '1':
            out.append((0, (), Fr(v))); continue
        if nm == 'c':
            out.append((1, (), Fr(v))); continue
        m = re.match(r'c\^(-?\d+) H\[(.*)\]', nm)
        out.append((int(m.group(1)), tuple(int(t) for t in m.group(2).split(',')), Fr(v)))
    return out


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


def hpl_side(terms, parity, weight=None):
    """sum q h_w(2n + parity - e) as canonical Z-sum expression (constants excluded: they only hit N=0,1)."""
    E = {}
    for e, w, q in terms:
        if w == ():
            continue
        if weight is not None and len(w) != weight:
            continue
        X = zsum.hpl_coeff_at(w, parity - e)
        qK = K(q.numerator)/K(q.denominator)
        for ww, c in X.items():
            zsum.add_into(E, ww, qK*c)
    return E


def low_order_check(terms, part, mu_exact_fn):
    """exact Taylor coefficients of the HPL side at c^N for N = -2..3 vs mu (mu_{3,0} := 0, odd := 0)."""
    ok = True
    for N in range(-2, 4):
        tot = Fr(0)
        for e, w, q in terms:
            if w == ():
                if N == e:
                    tot += q
                continue
            idx = N - e
            if idx >= 1:
                tot += q*coeffs_exact(w, idx + 2)[idx]
        target = Fr(0)
        if N >= 2 and N % 2 == 0:
            target = mu_exact_fn(N//2)[part]
        if tot != target:
            ok = False
            print(f'   low-order mismatch part {part} at c^{N}: {tot} vs {target}')
    return ok


if __name__ == '__main__':
    import mu3_symbolic as M
    Zmu = pickle.load(open('mu_zsum.pkl', 'rb'))
    allok = True
    for part, fname in FILES.items():
        terms = load(fname)
        muE = {tuple(w): toK(c) for w, c in Zmu[part]}
        ev = hpl_side(terms, 0)
        diff = zsum.lin((1, ev), (-1, muE))
        od = hpl_side(terms, 1)
        lo = low_order_check(terms, part, lambda nv: M.mu_exact(3, nv))
        print(f'part {part}: {len(terms)} HPL terms; even canonical form: {len(ev)} words; mu form: {len(muE)} words')
        print(f'   even identity  [c^(2n)] B3 - mu_(3,n) == 0 identically in n : {len(diff) == 0}')
        print(f'   odd identity   [c^(2n+1)] B3 == 0 identically in n          : {len(od) == 0}')
        print(f'   low orders c^-2..c^3 exact                                  : {lo}')
        if diff:
            print('   residual words:', list(diff.items())[:5])
        if od:
            print('   odd residual words:', list(od.items())[:5])
        allok = allok and not diff and not od and lo
        sys.stdout.flush()
    print('B3 PROVED for all n (all three parts):', allok)

    # ---- the weight-4 component: W_4(c) = R_3|_{wt 4} + z2 Psi|_{wt 2} + z3 (-Phi)|_{wt 1}
    out = {}
    for part, fname, wt in (('1', FILES['1'], 4), ('z2', FILES['z2'], 2), ('z3', FILES['z3'], 1)):
        terms = load(fname)
        E = hpl_side(terms, 0, weight=wt)
        out[part] = E
        print(f'\n[c^(2n)] W_4(c), {part}-part: {len(E)} Z-sum words (canonical form, coefficients in Q(n))')
        for w, c in sorted(E.items(), key=lambda x: (len(x[0]), x[0])):
            print(f'   Z{list(w)}(2n) : {sp.factor(c.as_expr())}')
    pickle.dump({p: [(w, str(c.as_expr())) for w, c in E.items()] for p, E in out.items()}, open('w4_taylor_zsum.pkl', 'wb'))
