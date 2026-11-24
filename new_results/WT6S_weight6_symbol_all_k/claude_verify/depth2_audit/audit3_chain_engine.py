"""Depth-two audit, part 3 (extra skepticism): regenerate the 12 pairing vectors pv(N_{(i,j),r}) with the W4-session
engine + my own pairing code; redo the binomial-basis rewrite of the TOWN1 deg-5 residue polynomials myself; assemble
pv(E_{5,m}) and P(k); compare to result.md table B / C; and use M5a's omega5.json to reproduce the k=5 pairings.
G-expression of N_{(i,j),r}: G_1((r,w);v) with (w,kk) in log_product(LT^i, MT^j), times log(2)^kk; endpoints e=a minus e=b
(same construction as M5a gexpr5.py / planner symbol_of, re-implemented here).
"""
import sys, json, time
from fractions import Fraction as F
from math import comb
sys.dont_write_bytecode = True
TH = "<project>/research/08_theory/"
WT = TH + "WT6S_weight6_symbol_all_k/"
sys.path.insert(0, TH + "A3_W4_proof")
from w4_symbol import Q, R, A, B, ZERO, ONE, add, scale, shuffle, log_symbol, concat, g_symbol, log_product, LT, MT, clean
import sympy as sp

X = {(0,): ONE}; Y = {(1,): ONE}
def br(a, b): return add(concat(a, b), scale(concat(b, a), -1))
P = br(X, Y); XP = br(X, P); YP = br(Y, P)
C42 = [br(X, XP), br(Y, XP), br(Y, YP)]
def pair(S, U, V):
    tot = ZERO
    for u, a in U.items():
        for v, b in V.items():
            tot = tot + a * b * (S.get(u + v, ZERO) - S.get(v + u, ZERO))
    return tot
def pv(S): return [pair(S, XP, YP)] + [pair(S, C, P) for C in C42]

POLE = {"0": ZERO, "a": A, "-a": -A, "b": B, "-b": -B}
LOG2 = log_symbol(Q(2))

def N_symbol(i, j, r):
    z = POLE[r]
    factors = (LT,) * i + (MT,) * j
    out = {}
    for (w, kk), d in log_product(*factors).items():
        word = (z,) + w
        assert len(word) + kk == 6
        s = add(g_symbol(ONE, word, A), scale(g_symbol(ONE, word, B), -1))
        for _ in range(kk):
            s = shuffle(s, LOG2)
        out = add(out, scale(s, d))
    return out

def fq(v):  # Q -> string
    return str(v.p) if v.q == 0 else (f"{v.q}*sqrt2" if v.p == 0 else f"{v.p}+{v.q}*sqrt2")

t0 = time.time()
KEYS = [(5, 0), (4, 1), (3, 2)]
POLES = ["0", "a", "-a", "-b"]
PVN = {}
for (i, j) in KEYS:
    for r in POLES:
        t1 = time.time()
        S = N_symbol(i, j, r)
        v = pv(S)
        PVN[(i, j, r)] = v
        print(f"pv(N_{i}{j},{r:>2s}) words={len(S):5d}  {[fq(x) for x in v]}   [{time.time() - t1:.1f}s]", flush=True)
print("total engine time %.1fs" % (time.time() - t0), flush=True)
json.dump({f"{i}{j}|{r}": [fq(x) for x in v] for (i, j, r), v in PVN.items()}, open("audit3_pvN_regenerated.json", "w"), indent=1)

# compare with WT6S n6_symbols_pairings.json (these were rational)
wt = json.load(open(WT + "n6_symbols_pairings.json"))["pairing_vectors"]
okN = True
for (i, j, r), v in PVN.items():
    ref = wt[f"{i}{j}|{r}"]
    mine = [fq(x) for x in v]
    okN &= (mine == ref)
print("pv(N_{ij,r}) regenerated == WT6S n6_symbols_pairings.json for all 12:", okN, flush=True)

# --- my own binomial rewrite of deg5 residue polynomials
k = sp.Symbol("k")
dp = json.load(open(WT + "deg5_residue_polynomials.json"))
def cellpoly(key, r):
    cell = dp[key][r]
    rat = sp.sympify(cell.get("rat", "0"), locals={"k": k})
    sq = sp.sympify(cell.get("sqrt2", "0"), locals={"k": k})
    return rat, sq
def newton_coeffs(poly):
    vals = [sp.Rational(poly.subs(k, 6 + t)) for t in range(5)]
    # forward differences at k=6
    out = []; row = vals[:]
    for m in range(5):
        out.append(row[0]); row = [row[t + 1] - row[t] for t in range(len(row) - 1)]
    # check
    chk = sum(out[m] * sp.binomial(k - 6, m) for m in range(5))
    assert sp.expand(sp.expand_func(chk) - sp.expand(poly)) == 0, "newton rewrite failed"
    return out
E = {}   # (i,j,r) -> list over m=1..5 of QS (p,q)
for (i, j) in KEYS:
    for r in POLES:
        rat, sq = cellpoly(f"{i}{j}", r)
        cr = newton_coeffs(rat); cq = newton_coeffs(sq)
        E[(i, j, r)] = [(cr[m], cq[m]) for m in range(5)]
eb = json.load(open(WT + "e5m_binomial_coeffs.json"))
okE = True
for (i, j, r), lst in E.items():
    ref = eb[f"{i}{j}|{r}"]
    for m in range(5):
        okE &= (sp.Rational(ref["rat"][m]) == lst[m][0] and sp.Rational(ref["sqrt2"][m]) == lst[m][1])
print("my binomial-basis rewrite of residue polynomials == e5m_binomial_coeffs.json (all 12 cells x 5 m):", okE)

def QQ(p, q): return Q(F(int(p.p), int(p.q)), F(int(q.p), int(q.q))) if False else Q(F(p.p, p.q), F(q.p, q.q))
def vscale_q(c, v): return [c * x for x in v]
def vadd(u, v): return [a + b for a, b in zip(u, v)]
def vzero(): return [ZERO] * 4
def fmt(v): return "[" + ", ".join(fq(x) for x in v) + "]"

PVE = {}
for m in range(1, 6):
    tot = vzero()
    for key, lst in E.items():
        c = QQ(*[sp.Rational(x) for x in lst[m - 1]]) if False else Q(F(lst[m - 1][0].p, lst[m - 1][0].q), F(lst[m - 1][1].p, lst[m - 1][1].q))
        tot = vadd(tot, vscale_q(c, PVN[key]))
    PVE[m] = tot
    print(f"pv(E_5,{m}) from engine-regenerated pv(N) and my e's: {fmt(tot)}")
TABLE_B = {
    1: [Q(0, 46), Q(F(-75, 2)), Q(0, 15), ZERO],
    2: [Q(0, F(-21, 2)), Q(F(-35, 2)), ZERO, Q(21)],
    3: [Q(0, -6), Q(-10), ZERO, Q(12)],
    4: [ZERO] * 4, 5: [ZERO] * 4,
}
okB = all(PVE[m] == TABLE_B[m] for m in range(1, 6))
print("REGENERATED pv(E_5,m) == result.md table B (all m=1..5):", okB, flush=True)

# --- direct P(k) = sum res_{ij,r}(k) pv(N_{ij,r}) at k = 6..52 vs stated polynomials, independent of binomial basis
rs = sp.Symbol("r")
stated = [(-3 * k**2 + sp.Rational(57, 2) * k - 17) * rs, -5 * k**2 + sp.Rational(95, 2) * k - sp.Rational(285, 2), 15 * rs, 6 * k**2 - 57 * k + 126]
def Qsp(x):
    x = sp.expand(x)
    return Q(F(str(x.coeff(rs, 0))), F(str(x.coeff(rs, 1))))
okP = True
for kv in range(6, 53):
    tot = vzero()
    for (i, j, r) in PVN:
        rat, sq = cellpoly(f"{i}{j}", r)
        c = Q(F(str(sp.Rational(rat.subs(k, kv)))), F(str(sp.Rational(sq.subs(k, kv)))))
        tot = vadd(tot, vscale_q(c, PVN[(i, j, r)]))
    want = [Qsp(stated[t].subs(k, kv)) for t in range(4)]
    okP &= (tot == want)
print("DIRECT sum res(k)*pv(N) == stated power-basis P(k) at every k=6..52:", okP, flush=True)

# --- k=5 from M5a's omega5.json (independent derivation of k=5 kernel), and comparison of residues with the quartic family at k=5
om = json.load(open(TH + "M5a_next_coefficient/omega5.json"))
tot5 = vzero()
resid5 = {}
for (i, j) in KEYS:
    for name, (p, q) in om[f"{i},{j}"]:
        c = Q(F(p), F(q))
        resid5[(i, j, name)] = resid5.get((i, j, name), ZERO) + c   # pole '0' is listed twice (rational part, sqrt2 part): accumulate
        tot5 = vadd(tot5, vscale_q(c, PVN[(i, j, name)]))
print("k=5 pairings from M5a omega5.json residues with regenerated pv(N):", fmt(tot5), "  (M5a: [449/4*sqrt2, -30, 15*sqrt2, -9])")
okk5 = tot5 == [Q(0, F(449, 4)), Q(-30), Q(0, 15), Q(-9)]
print("matches M5a k=5 target:", okk5)
print("k=5 residue cells: omega5.json residue vs quartic-family polynomial evaluated at k=5, and difference:")
tot_diff = vzero()
for (i, j) in KEYS:
    for r in POLES:
        rat, sq = cellpoly(f"{i}{j}", r)
        poly5 = Q(F(str(sp.Rational(rat.subs(k, 5)))), F(str(sp.Rational(sq.subs(k, 5)))))
        true5 = resid5.get((i, j, r), ZERO)
        d = true5 - poly5
        if d != ZERO:
            print(f"   ({i},{j}) pole {r:>2s}: omega5 = {fq(true5):>14s}   quartic(k=5) = {fq(poly5):>14s}   diff = {fq(d)}")
        tot_diff = vadd(tot_diff, vscale_q(d, PVN[(i, j, r)]))
print("pairing of (true k=5 kernel - quartic-family k=5 kernel):", fmt(tot_diff), "  expected [247/4*sqrt2,0,0,0]")
print("expected match:", tot_diff == [Q(0, F(247, 4)), ZERO, ZERO, ZERO])
print("poles in omega5 absent from the 12-cell set (should be none):", [kk for kk in resid5 if (kk[0], kk[1], kk[2]) not in PVN])
print("ALL DONE")
