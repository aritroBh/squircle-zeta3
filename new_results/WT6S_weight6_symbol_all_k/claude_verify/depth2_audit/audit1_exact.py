"""Depth-two audit, part 1: exact recomputation in Q(sqrt2) with my own arithmetic.
Reads (read-only): M5a atom_pairings.json, M5a/M5b n6_solution.json, WT6S result.md numbers
(transcribed below by hand from result.md sections B, C, D), WT6S data files for cross-compare.
"""
import json, hashlib, itertools
from fractions import Fraction as F
from math import comb
import sympy as sp

TH = "<project>/research/08_theory/"
M5A = TH + "M5a_next_coefficient/"
M5B = TH + "M5b_w6_closed_form/"
WT = TH + "WT6S_weight6_symbol_all_k/"


class QS:
    """p + q*sqrt2 with Fraction p, q (own implementation)."""
    __slots__ = ("p", "q")
    def __init__(s, p=0, q=0):
        s.p = F(p); s.q = F(q)
    def __add__(s, o):
        o = qs(o); return QS(s.p + o.p, s.q + o.q)
    __radd__ = __add__
    def __sub__(s, o):
        o = qs(o); return QS(s.p - o.p, s.q - o.q)
    def __rsub__(s, o):
        return qs(o) - s
    def __neg__(s):
        return QS(-s.p, -s.q)
    def __mul__(s, o):
        o = qs(o); return QS(s.p * o.p + 2 * s.q * o.q, s.p * o.q + s.q * o.p)
    __rmul__ = __mul__
    def inv(s):
        n = s.p * s.p - 2 * s.q * s.q
        return QS(s.p / n, -s.q / n)
    def __truediv__(s, o):
        return s * qs(o).inv()
    def __eq__(s, o):
        o = qs(o); return s.p == o.p and s.q == o.q
    def __hash__(s):
        return hash((s.p, s.q))
    def iszero(s):
        return s.p == 0 and s.q == 0
    def __repr__(s):
        if s.q == 0: return str(s.p)
        if s.p == 0: return f"({s.q})*sqrt2"
        return f"({s.p} + {s.q}*sqrt2)"


def qs(x):
    return x if isinstance(x, QS) else QS(x, 0)

R = QS(0, 1)
ZERO = QS(0, 0)


def parse(sx):
    """parse a string such as '-21/2*sqrt2', '15*sqrt2', '-75/2', '(-1+1*sqrt2)' into QS via sympy"""
    r = sp.Symbol("r")
    e = sp.sympify(sx.replace("sqrt2", "r"), locals={"r": r})
    e = sp.expand(e)
    pol = sp.Poly(e, r)
    c = pol.all_coeffs()[::-1]
    assert len(c) <= 2, sx
    c += [0] * (2 - len(c))
    return QS(F(str(sp.Rational(c[0]))), F(str(sp.Rational(c[1]))))


def vec(lst):
    return [parse(x) if isinstance(x, str) else qs(x) for x in lst]


def vadd(u, v): return [a + b for a, b in zip(u, v)]
def vscale(c, u): return [qs(c) * a for a in u]
def vzero(): return [ZERO] * 4
def vfmt(u): return "[" + ", ".join(repr(a) for a in u) + "]"


def lin(coeffs, vecs):
    out = vzero()
    for c, v in zip(coeffs, vecs):
        out = vadd(out, vscale(c, v))
    return out


def rat_det(M):
    """exact determinant of a Fraction matrix by Bareiss-free Gaussian elimination"""
    n = len(M); A = [row[:] for row in M]; det = F(1)
    for i in range(n):
        piv = next((r for r in range(i, n) if A[r][i] != 0), None)
        if piv is None: return F(0)
        if piv != i:
            A[i], A[piv] = A[piv], A[i]; det = -det
        det *= A[i][i]
        for r in range(i + 1, n):
            f = A[r][i] / A[i][i]
            for c in range(i, n):
                A[r][c] -= f * A[i][c]
    return det


results = {}
def report(name, ok, extra=""):
    results[name] = ok
    print(f"[{'PASS' if ok else 'FAIL'}] {name} {extra}")


# ---------------------------------------------------------------- (1) atoms
print("=" * 70); print("(1) ATOM PAIRING VECTORS from M5a atom_pairings.json (order 33, 42a, 42b, 42c)")
ap = json.load(open(M5A + "atom_pairings.json"))
ATOM_KEYS = ["Li_33(b,1)", "Li_33(1/sqrt2,1)", "Li_33(b,2)", "Li_42(b,1/sqrt2)"]
ATOM_LBL = ["Li33(b,1)", "Li33(1/sqrt2,1)", "Li33(b,2)", "Li42(b,1/sqrt2)"]
A = []
for key in ATOM_KEYS:
    v = vec(ap[key]); A.append(v); print(f"  {key:20s} {vfmt(v)}")
# sanity: M5a's quoted controls
report("M5a table: Li42(b,1) = (0,0,0,-1/4)", vec(ap["Li_42(b,1)"]) == vec(["0", "0", "0", "-1/4"]))
report("M5a table: Li33(b,1) = (0,0,0,3/4)", A[0] == vec(["0", "0", "0", "3/4"]))

# does M5a's n6_solution.json atom list correspond to these four?
sa = json.load(open(M5A + "n6_solution.json")); sb = json.load(open(M5B + "n6_solution.json"))
report("M5a n6_solution.json == M5b n6_solution.json (parsed)", sa == sb,
       "sha256 equal bytes: %s" % (hashlib.sha256(open(M5A + 'n6_solution.json', 'rb').read()).hexdigest() ==
                                 hashlib.sha256(open(M5B + 'n6_solution.json', 'rb').read()).hexdigest()))
names = [(a["m"], a["n"], parse(a["x"]), parse(a["y"])) for a in sa["atoms"]]
expect = [(3, 3, QS(-1, 1), QS(1)), (3, 3, QS(0, F(1, 2)), QS(1)), (3, 3, QS(-1, 1), QS(2)), (4, 2, QS(-1, 1), QS(0, F(1, 2)))]
report("n6_solution.json atoms == (Li33(b,1), Li33(1/sqrt2,1), Li33(b,2), Li42(b,1/sqrt2)) with b=sqrt2-1",
       names == expect)
n6c = [QS(F(a["coeff_p"]), F(a["coeff_q"])) for a in sa["atoms"]]
print("  N6 coefficients (M5a/M5b):", [repr(c) for c in n6c])
n6pv = lin(n6c, A)
tgt6 = vec(sa["target33"].split()[:1]) + vec(sa["target42"])
tgt6 = [parse(sa["target33"])] + [parse(x) for x in sa["target42"]]
print("  pairings of N6 :", vfmt(n6pv)); print("  M5a W6(k=5)    :", vfmt(tgt6))
report("N6 pairings == W6(k=5) pairings (449sqrt2/4, -30, 15sqrt2, -9)", n6pv == tgt6)

M = [[a.p for a in row] for row in A]   # all atom pairings are rational
report("atom pairing entries all rational (sqrt2 part 0)", all(a.q == 0 for row in A for a in row))
det = rat_det(M)
print("  4x4 determinant (rows = atoms in WT6S order, cols = 33,42a,42b,42c):", det)
report("|det| == 9/512", abs(det) == F(9, 512), f"(det = {det})")
dsp = sp.Matrix([[sp.Rational(x.numerator, x.denominator) for x in row] for row in M]).det()
report("sympy det agrees", sp.Rational(det.numerator, det.denominator) == dsp, f"(sympy {dsp})")

# solver: c^T M = t, c in Q(sqrt2)^4
Minv = sp.Matrix([[sp.Rational(x.numerator, x.denominator) for x in row] for row in M]).inv()
def solve(t):
    t0 = sp.Matrix([[sp.Rational(a.p.numerator, a.p.denominator) for a in t]])
    t1 = sp.Matrix([[sp.Rational(a.q.numerator, a.q.denominator) for a in t]])
    c0 = t0 * Minv; c1 = t1 * Minv
    return [QS(F(str(c0[i])), F(str(c1[i]))) for i in range(4)]
c = solve(tgt6)
print("  independent solve for M5a target:", [repr(x) for x in c])
report("solve(M5a target) == M5a N6 coefficients", c == n6c)

# ---------------------------------------------------------------- (1b) WT6S N(E_{5,m}) from result.md section D
print("=" * 70); print("(1b) WT6S N(E_{5,m}) (result.md section D, hand-transcribed) -> pairings, vs section B table")
r2 = QS(0, 1)
ND = {  # coefficients on (Li33(b,1), Li33(1/sqrt2,1), Li33(b,2), Li42(b,1/sqrt2))
    1: [-56 * r2, QS(-200, 80), 16 * r2, 120 * r2],
    2: [QS(28, -42), QS(F(-280, 3)), QS(0, F(-21, 2)), ZERO],
    3: [QS(16, -24), QS(F(-160, 3)), QS(0, -6), ZERO],
    4: [ZERO] * 4, 5: [ZERO] * 4,
}
TABLE_B = {  # result.md section B, (3,3),(4,2)a,(4,2)b,(4,2)c
    1: [QS(0, 46), QS(F(-75, 2)), QS(0, 15), ZERO],
    2: [QS(0, F(-21, 2)), QS(F(-35, 2)), ZERO, QS(21)],
    3: [QS(0, -6), QS(-10), ZERO, QS(12)],
    4: vzero(), 5: vzero(),
}
for m in range(1, 6):
    got = lin(ND[m], A)
    print(f"  m={m}: pairings of stated N(E_5,{m}) = {vfmt(got)}   table B = {vfmt(TABLE_B[m])}")
    report(f"pairings(N(E_5,{m})) == table B row", got == TABLE_B[m])
    sol = solve(TABLE_B[m])
    ref = ND[m]
    report(f"independent solve of table-B row {m} reproduces stated N(E_5,{m}) coefficients", sol == ref,
           f"(solve = {[repr(x) for x in sol]})")

# cross-compare with WT6S data files (not independent; transcription guard)
ne = json.load(open(WT + "n_e5m_solution.json"))
atom_json = ["Li_{3,3}(b,1)", "Li_{3,3}(1/sqrt2,1)", "Li_{3,3}(b,2)", "Li_{4,2}(b,1/sqrt2)"]
okj = True
for m in range(1, 6):
    cj = ne["solutions"][f"E5,{m}"]["coeffs"]
    vv = [QS(F(cj[a][0]), F(cj[a][1])) if a in cj else ZERO for a in atom_json]
    okj &= (vv == ND[m])
report("transcription of section D == n_e5m_solution.json", okj)
report("n_e5m_solution.json det string is 9/512", ne["det"] == "9/512")
ep = json.load(open(WT + "e5m_pairings.json"))["pairings"]
report("transcription of section B == e5m_pairings.json", all(vec(ep[f"E5,{m}"]) == TABLE_B[m] for m in range(1, 6)))

# ---------------------------------------------------------------- (2) P(k)
print("=" * 70); print("(2) P(k) = sum_m C(k-6,m-1) pv(E_5,m)  vs stated power-basis polynomials")
k = sp.Symbol("k"); rs = sp.Symbol("r")  # rs = sqrt2
def tosym(x): return sp.Rational(x.p.numerator, x.p.denominator) + sp.Rational(x.q.numerator, x.q.denominator) * rs
def binomsym(m1):
    return sp.expand(sp.prod([k - 6 - j for j in range(m1)]) / sp.factorial(m1))
Psym = []
for t in range(4):
    expr = sum(binomsym(m - 1) * tosym(TABLE_B[m][t]) for m in range(1, 6))
    Psym.append(sp.expand(expr))
stated = [
    (-3 * k**2 + sp.Rational(57, 2) * k - 17) * rs,
    -5 * k**2 + sp.Rational(95, 2) * k - sp.Rational(285, 2),
    15 * rs,
    6 * k**2 - 57 * k + 126,
]
names4 = ["P33", "P42a", "P42b", "P42c"]
for t in range(4):
    diff = sp.expand(Psym[t] - stated[t])
    print(f"  {names4[t]}(k) from table = {Psym[t]};  stated = {sp.expand(stated[t])};  diff = {diff}")
    report(f"{names4[t]}(k) polynomial identity (symbolic in k)", diff == 0)
# pointwise at every integer k in 6..52 with math.comb (independent of sympy binomial polynomials)
okpt = True
for kv in range(6, 53):
    for t in range(4):
        s = ZERO
        for m in range(1, 6):
            s = s + comb(kv - 6, m - 1) * TABLE_B[m][t]
        st = stated[t].subs(k, kv)
        st = sp.expand(st)
        sq = QS(F(str(st.coeff(rs, 0))), F(str(st.coeff(rs, 1))))
        okpt &= (s == sq)
report("P(k) == stated polynomial, exact, at every integer k = 6..52 (math.comb)", okpt)

# depth-two coefficient polynomials
print("-" * 70); print("depth-two part of W6(k): sum_m C(k-6,m-1) N(E_5,m) coefficient polynomials")
D2_stated = {
    0: sp.expand((8 * k**2 - 76 * k + 168) + (-12 * k**2 + 114 * k - 308) * rs),
    1: sp.expand(-sp.Rational(80, 3) * k**2 + sp.Rational(760, 3) * k - 760 + 80 * rs),
    2: sp.expand((-3 * k**2 + sp.Rational(57, 2) * k - 47) * rs),
    3: sp.expand(120 * rs),
}
for a in range(4):
    poly = sp.expand(sum(binomsym(m - 1) * tosym(ND[m][a]) for m in range(1, 6)))
    diff = sp.expand(poly - D2_stated[a])
    print(f"  {ATOM_LBL[a]:18s}: {poly}   diff vs stated: {diff}")
    report(f"depth-two coefficient polynomial of {ATOM_LBL[a]}", diff == 0)
nw = json.load(open(WT + "n_w6k_polynomial.json"))
for a, jk in enumerate(atom_json):
    pj = sp.expand(sp.sympify(nw[jk]["rat"], locals={"k": k}) + sp.sympify(nw[jk]["sqrt2"], locals={"k": k}) * rs)
    report(f"n_w6k_polynomial.json {ATOM_LBL[a]} == stated", sp.expand(pj - D2_stated[a]) == 0)
# pairing of the depth-two polynomial combination equals P(k), symbolically
for t in range(4):
    expr = sp.expand(sum(D2_stated[a] * sp.Rational(A[a][t].p.numerator, A[a][t].p.denominator) for a in range(4)))
    # reduce rs^2 -> 2 (only products of D2 coefficient (linear in rs) with rational atoms: no rs^2 appears)
    report(f"pairing of depth-two polynomial combination == {names4[t]}(k) (symbolic)", sp.expand(expr - stated[t]) == 0)
okk = True
for kv in range(6, 53):
    co = [QS(F(str(sp.expand(D2_stated[a].subs(k, kv)).coeff(rs, 0))), F(str(sp.expand(D2_stated[a].subs(k, kv)).coeff(rs, 1)))) for a in range(4)]
    got = lin(co, A)
    want = []
    for t in range(4):
        st = sp.expand(stated[t].subs(k, kv))
        want.append(QS(F(str(st.coeff(rs, 0))), F(str(st.coeff(rs, 1)))))
    okk &= (got == want)
report("depth-two combination pairing == P(k) exactly at every k = 6..52", okk)

# ---------------------------------------------------------------- (3) k = 5
print("=" * 70); print("(3) k = 5 observation")
P5 = [QS(F(str(sp.expand(stated[t].subs(k, 5)).coeff(rs, 0))), F(str(sp.expand(stated[t].subs(k, 5)).coeff(rs, 1)))) for t in range(4)]
print("  stated polynomials at k=5:", vfmt(P5))
claimed = [QS(0, F(101, 2)), QS(-30), QS(0, 15), QS(-9)]
report("polynomials at k=5 == claimed (101sqrt2/2, -30, 15sqrt2, -9)", P5 == claimed)
Palt = [ZERO] * 4
for m in range(1, 6):
    Palt = vadd(Palt, vscale(comb(-1 + 0, 0) if False else (-1) ** (m - 1), TABLE_B[m]))
print("  via C(-1,m-1) = (-1)^(m-1) alternating sum of table B:", vfmt(Palt))
report("alternating-sum extrapolation == polynomial at k=5", Palt == P5)
true5 = [QS(0, F(449, 4)), QS(-30), QS(0, 15), QS(-9)]
report("M5a true k=5 target == (449sqrt2/4, -30, 15sqrt2, -9)", true5 == tgt6)
diff5 = vadd(true5, vscale(-1, P5))
print("  true(k=5) - extrapolated =", vfmt(diff5))
report("difference is purely in (3,3), equal to +247sqrt2/4 (extrapolated - true = -247sqrt2/4)",
       diff5 == [QS(0, F(247, 4)), ZERO, ZERO, ZERO])
# depth-two coefficient at k=5 vs N6
co5 = [QS(F(str(sp.expand(D2_stated[a].subs(k, 5)).coeff(rs, 0))), F(str(sp.expand(D2_stated[a].subs(k, 5)).coeff(rs, 1)))) for a in range(4)]
print("  depth-two coefficients extrapolated to k=5:", [repr(x) for x in co5])
print("  M5a N6 coefficients                        :", [repr(x) for x in n6c])
dco = [n6c[i] - co5[i] for i in range(4)]
print("  N6 - extrapolated:", [repr(x) for x in dco])
report("N6 - extrapolated(k=5) = 247sqrt2*Li33(b,1) + (247sqrt2/4)*Li33(b,2), others 0",
       dco == [QS(0, 247), ZERO, QS(0, F(247, 4)), ZERO])
print("  pairing of that difference:", vfmt(lin(dco, A)))

print("=" * 70)
print("SUMMARY: %d checks, %d FAIL" % (len(results), sum(1 for v in results.values() if not v)))
for kname, v in results.items():
    if not v: print("  FAILED:", kname)
