"""WT6S step D2: depth-two polynomial N(W6(k)) = sum_m C(k-6,m-1) N(E_{5,m}) in closed form.

Reads n_e5m_solution.json; builds each atom's coefficient as a polynomial in k (power basis);
verifies (exactly, via pairing vectors) that N(W6(k)) has the W6(k) pairing polynomial at
k=6,7,10,21,52; prints the k=5 extrapolation vs actual-N6 observation."""
import json
from fractions import Fraction as F
import sympy as sp
import qsqrt2
from w4_symbol import Q, R, B, ONE as ONEQ, ban_float, self_check
import wt6s_core as C6
from wt6s_core import depth2_symbol, pairing_vector6, add, scale, ZEROQ, log

k = sp.Symbol('k')
ATOMS = [
    ("Li_{3,3}(b,1)", 3, 3, B, ONEQ),
    ("Li_{3,3}(1/sqrt2,1)", 3, 3, R / 2, ONEQ),
    ("Li_{3,3}(b,2)", 3, 3, B, Q(2)),
    ("Li_{4,2}(b,1/sqrt2)", 4, 2, B, R / 2),
]
BINOM2 = (k - 6) * (k - 7) / 2

def binom(m1):   # C(k-6, m-1)
    from math import comb
    n = m1
    return sp.expand(sp.prod([k - 6 - j for j in range(n)]) / sp.factorial(n))

def main():
    ban_float(); self_check()
    sol = json.load(open('n_e5m_solution.json'))['solutions']
    # coeff[atom] = polynomial in k (QS2 coefficients)
    table = {}
    for name, *_ in ATOMS:
        pr = sq = sp.Integer(0)
        for m in range(1, 6):
            c = sol[f"E5,{m}"]['coeffs'].get(name)
            if c is None:
                continue
            pr += sp.Rational(c[0]) * binom(m - 1)
            sq += sp.Rational(c[1]) * binom(m - 1)
        table[name] = (sp.expand(pr), sp.expand(sq))
    print("N(W6(k)) = sum over atoms of coeff(k) * atom:")
    for name, (pr, sq) in table.items():
        print(f"  ({pr} + ({sq})*sqrt2) * {name}")
    # verify pairing polynomial at integer k via depth2 symbols
    ePv = json.load(open('e5m_pairings.json'))['pairings_pq']
    epv = {m: [Q(F(p), F(q)) for p, q in ePv[f"E5,{m}"]] for m in range(1, 6)}
    ok = True
    for kv in [6, 7, 10, 21, 52]:
        S = {}
        for name, m, n, x, y in ATOMS:
            pr, sq = table[name]
            c = Q(F(sp.Rational(pr.subs(k, kv))), F(sp.Rational(sq.subs(k, kv))))
            if c == ZEROQ:
                continue
            S = add(S, scale(depth2_symbol(m, n, x, y), c))
        got = pairing_vector6(S) if S else [ZEROQ] * 4
        want = [sum((epv[m][t] * int(sp.binomial(kv - 6, m - 1)) for m in range(1, 6)), ZEROQ)
                for t in range(4)]
        ok &= all(got[t] == want[t] for t in range(4))
    print(f"N(W6(k)) pairings == W6(k) pairing polynomial at k=6,7,10,21,52: {ok}")
    assert ok
    # k=5 observation
    S5 = {}
    for name, m, n, x, y in ATOMS:
        pr, sq = table[name]
        c = Q(F(sp.Rational(pr.subs(k, 5))), F(sp.Rational(sq.subs(k, 5))))
        if c == ZEROQ:
            continue
        S5 = add(S5, scale(depth2_symbol(m, n, x, y), c))
    pv5 = pairing_vector6(S5)
    print("N(W6(k)) extrapolated to k=5 pairings:", [str(v) for v in pv5])
    print("  (actual N6: 449/4*sqrt2, -30, 15*sqrt2, -9; difference only in (3,3))")
    json.dump({name: {'rat': str(table[name][0]), 'sqrt2': str(table[name][1])} for name, *_ in ATOMS},
              open('n_w6k_polynomial.json', 'w'), indent=1)
    with open('progress.log', 'a') as f:
        f.write('step_D2: N(W6(k)) atom-coefficient polynomials; pairing match k=6,7,10,21,52; k=5 observation\n')
    print("saved n_w6k_polynomial.json")

if __name__ == '__main__':
    main()
