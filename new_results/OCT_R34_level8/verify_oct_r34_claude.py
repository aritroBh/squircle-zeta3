"""Second, independent exact proof of the octagon depth-2 reductions R3/R4 (Claude, 2026-10-22).

Targets (OCT6_octagon_closed/reductions.pkl, keys 'Lam', 'swap'), with Li_{2,1}(x,y) = sum_{n>m>=1} x^n y^m/(n^2 m):
  Lam  = Li_{2,1}(1-sqrt2, -sqrt2/2),   swap = Li_{2,1}(1-sqrt2/2, -sqrt2).
Method: Li_{2,1}(x,y) = G(0, 1/x, 1/(xy); 1) (series expansion of the iterated integral). For both targets the letters are
{0, -A, 2+sqrt2} (A = 1+sqrt2). The Moebius frame (0,1,inf) -> (0, B, A) (B = sqrt2-1) sends -A -> -A and 2+sqrt2 -> 1, so
G pulls back to level-8 H-words by cmzv.forms_path with forms dlog(g(a), A), path 0 -> B (no configuration point inside
(0, B), and A is not between 0 and B). The reductions are built from cmzv atoms; Li_n(-sqrt2/2) uses duplication
Li_n(-z) = 2^(1-n) Li_n(z^2) - Li_n(z). DIFF = vector(target) - vector(reduction) is then tested for membership in the
proved level-8 base and given an exact certificate, as for the chi3 bridge.
Run from research/08_theory:  ../../.venv/bin/python OCT_R34_level8/verify_oct_r34_claude.py
"""
import sys, pickle
from fractions import Fraction as Fr
sys.path.insert(0, 'research/08_theory/A3f4_glanois_basis/Q')
import cmzv as C
import mpmath as mp
import sympy as sp
import flint
from sympy.ntheory import sqrt_mod

mp.mp.dps = 150
Q, A, B, S2 = C.Q, C.A, C.B, C.SQ2
V = pickle.load(open('A3f6_level16/vals8_150digit.pickle', 'rb'))
val = lambda vec: mp.fsum(mp.mpf(c.numerator) / c.denominator * mp.mpc(*V[w]) for w, c in vec.items())
scale = lambda vec, c: {w: a * Fr(c) for w, a in vec.items()}


def li(n, z):
    try:
        return C.atom(('li', n, z))
    except ValueError:
        # duplication for z = -t: Li_n(-t) = 2^(1-n) Li_n(t^2) - Li_n(t)
        t = -z
        return C.add((C.atom(('li', n, t * t)), Fr(1, 2 ** (n - 1))), (C.atom(('li', n, t)), -1))


pts = {"x1": Q(1) - S2, "x2": Q(1) - S2 / 2, "k": B, "d": Q(2) - S2, "r2h": S2 / 2, "mr2h": -S2 / 2}
fpts = {k: mp.mpf(v.p.numerator) / v.p.denominator + mp.mpf(v.q.numerator) / v.q.denominator * mp.sqrt(2) for k, v in pts.items()}
atoms = {"l2": C.atom(('log', Q(2))), "lu": C.atom(('log', A)), "z2": scale(C.atom(('pi2',)), Fr(1, 6)), "z3": C.atom(('z3',))}
for nm in pts:
    for n in (2, 3):
        atoms[f"L{n}{nm}"] = li(n, pts[nm])
ref = {"l2": mp.log(2), "lu": mp.log(1 + mp.sqrt(2)), "z2": mp.zeta(2), "z3": mp.zeta(3)}
ref.update({f"L{n}{nm}": mp.polylog(n, fpts[nm]) for nm in pts for n in (2, 3)})
w3 = [k for k in atoms if all(len(w) == 3 for w in atoms[k])]   # the value table holds weight-3 words only
worst = max(abs(val(atoms[k]) - ref[k]) for k in w3)
print(f"weight-3 atoms {w3}: worst |vector - mpmath| = {mp.nstr(worst, 3)} (lower-weight atoms are checked via the products below)")


def depth2(x, y):
    rho = max(abs(x), abs(x * y)); N = int(mp.ceil(mp.log(mp.mpf(10) ** (-mp.mp.dps - 5) * (1 - rho)) / mp.log(rho)))
    h, xyp, tot = mp.mpf(0), x * y, mp.mpf(0)
    for n in range(1, N + 1):
        tot += h / n ** 2; h = x * h + x * xyp / n; xyp *= x * y
    return tot


def expr_vec(e):
    """Evaluate a sympy polynomial in the atom symbols in the shuffle algebra."""
    out = {}
    for term in sp.Add.make_args(sp.expand(e)):
        coeff, rest = term.as_coeff_Mul()
        vec = {(): Fr(1)}
        for base, ex in rest.as_powers_dict().items():
            for _ in range(int(ex)):
                vec = C.mul(vec, atoms[str(base)])
        out = C.add((out, 1), (vec, Fr(sp.Rational(coeff).p, sp.Rational(coeff).q)))
    return {w: c for w, c in out.items() if c != 0}


red = pickle.load(open('OCT6_octagon_closed/reductions.pkl', 'rb'))
P = C.POINTS
targets = {
    "Lam": ([C.dlog(Q(0), A), C.dlog(-A, A), C.dlog(Q(1), A)], (1 - mp.sqrt(2), -mp.sqrt(2) / 2)),
    "swap": ([C.dlog(Q(0), A), C.dlog(Q(1), A), C.dlog(-A, A)], (1 - mp.sqrt(2) / 2, -mp.sqrt(2))),
}
W = pickle.load(open('A3f3_level8_compute/weight3_with_lifts.pickle', 'rb'))
K = pickle.load(open('A3f3_level8_compute/k3_newHrows.pickle', 'rb'))['krows']
cols = list(W['cols']); idx = {c: i for i, c in enumerate(cols)}; n = len(cols)


def modp(x, p):
    x = Fr(x); return x.numerator % p * pow(x.denominator % p, -1, p) % p


def accum(r):
    out = {}
    for w, a in r: out[w] = out.get(w, 0) + Fr(a)
    return out


Wrows = [accum(r) for r in W['rows']]
for key, (forms, (x, y)) in targets.items():
    G = C.forms_path(forms, P[Q(0)], P[B])
    print(f"\n{key}: vector {len(G)} words; vector - depth2 series = {mp.nstr(abs(val(G) - depth2(x, y)), 3)}")
    R = expr_vec(red[key])
    print(f"{key}: reduction vector - depth2 series = {mp.nstr(abs(val(R) - depth2(x, y)), 3)}")
    DIFF = {w: c for w, c in C.add((G, 1), (R, -1)).items() if c != 0}
    assert all(w in idx for w in DIFF)
    print(f"{key}: DIFF {len(DIFF)} words, DIFF . values = {mp.nstr(abs(val(DIFF)), 3)}")
    inspan = True
    for p in (1000000033, 1000000241):
        M = [[0] * n for _ in Wrows]
        for i, r in enumerate(Wrows):
            for w, a in r.items(): M[i][idx[w]] = (M[i][idx[w]] + modp(a, p)) % p
        d = [modp(DIFF.get(c, 0), p) for c in cols]
        r0 = flint.nmod_mat(len(M), n, [v for r in M for v in r], p).rank()
        r1 = flint.nmod_mat(len(M) + 1, n, [v for r in M for v in r] + d, p).rank()
        print(f"{key}: p={p} rational base rank {r0}, +DIFF {r1}, in span: {r0 == r1}")
        inspan &= (r0 == r1)
    if not inspan:
        print(f"{key}: not in the rational base span; (Kummer / R1-R5 test not implemented here)"); continue
    # exact certificate: independent rows mod p, pivot columns, exact solve, exact check on all words
    p = 1000000033
    M = [[0] * n for _ in Wrows]
    for i, r in enumerate(Wrows):
        for w, a in r.items(): M[i][idx[w]] = (M[i][idx[w]] + modp(a, p)) % p
    E, rk = flint.nmod_mat(n, len(M), [M[i][j] for j in range(n) for i in range(len(M))], p).rref()
    piv = [next(j for j in range(len(M)) if int(E[i, j]) != 0) for i in range(rk)]
    E2, _ = flint.nmod_mat(rk, n, [v for i in piv for v in M[i]], p).rref()
    pc = [next(j for j in range(n) if int(E2[i, j]) != 0) for i in range(rk)]
    fq = lambda f: flint.fmpq(Fr(f).numerator, Fr(f).denominator)
    Am = flint.fmpq_mat(rk, rk, [fq(Wrows[piv[i]].get(cols[j], 0)) for j in pc for i in range(rk)])
    bv = flint.fmpq_mat(rk, 1, [fq(DIFF.get(cols[j], 0)) for j in pc])
    sol = Am.solve(bv)
    c = [Fr(int(sol[i, 0].p), int(sol[i, 0].q)) for i in range(rk)]
    comb = {}
    for ci, i in zip(c, piv):
        if ci:
            for w, a in Wrows[i].items(): comb[w] = comb.get(w, 0) + ci * a
    comb = {w: a for w, a in comb.items() if a != 0}
    used = [i for ci, i in zip(c, piv) if ci]
    print(f"{key}: EXACT CERTIFICATE {'holds on all 576 words' if comb == DIFF else 'FAILED'}; "
          f"{len(used)} rows used, max row index {max(used) if used else None} (Rogers lifts are >= 1298)")
    pickle.dump({'piv_rows': piv, 'multipliers': c, 'DIFF': DIFF, 'ok': comb == DIFF},
                open(f'OCT_R34_level8/certificate_{key}.pickle', 'wb'))
