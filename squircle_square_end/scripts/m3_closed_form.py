"""Identify W4 (hence M3) in a natural Li4 basis, check stability, and write the experimental closed form.

1. Greedy selection of independent Li4's in a 'natural' order (1/2, b, -b, 1/sqrt2, b^2, -b^2, 2-sqrt2, ...)
   modulo the 22 weight-4 products (LLL at 1000 digits).
2. LLL on W4 + products + selected Li4's at D digits, for D in a list, with M3 from method A and from method B:
   the relation must be the same (up to sign) in every run.
3. Exact rational form of W4 and of M3 = W4 + W3' + T, printed in plain text and LaTeX; numeric residuals.
"""
import os, sys, pickle, itertools, time
import mpmath as mp
import flint
from fractions import Fraction as Fr
from qsqrt2 import QS2, SQ2, A, B, sunit_log
import polylog_exact as PE
import sympy as sp

HERE = os.path.dirname(os.path.abspath(__file__))
OUTD = os.path.join(HERE, "..", "outputs")
mp.mp.dps = 1220
T0 = time.time()
lines = []


def say(*a):
    s = " ".join(str(x) for x in a)
    print(s, flush=True)
    lines.append(s)


tails = pickle.load(open(os.path.join(OUTD, "m3_tails_exact.pkl"), "rb"))
mp.mp.dps = 1220
W3p, T = tails["W3p"].num(), tails["T"].num()
M3 = {m: mp.mpf(open(os.path.join(OUTD, f"M3_method{m}_1200.txt")).read().strip()) for m in "AB"}
W4 = {m: M3[m] - W3p - T for m in "AB"}

r2 = mp.sqrt(2); lam = mp.log(2); L = mp.log(1 + r2); b = r2 - 1
consts = {"lam": lam, "L": L, "pi2": mp.pi ** 2, "Li2b": mp.polylog(2, b), "z3": mp.zeta(3),
          "Li3b": mp.polylog(3, b), "Li3mb": mp.polylog(3, -b), "Li3h": mp.polylog(3, 1 / r2)}
w1, w2, w3 = ["lam", "L"], ["pi2", "Li2b"], ["z3", "Li3b", "Li3mb", "Li3h"]
prods = [tuple(c) for c in itertools.combinations_with_replacement(w1, 4)]
prods += [(x,) + c for x in w2 for c in itertools.combinations_with_replacement(w1, 2)]
prods += [tuple(c) for c in itertools.combinations_with_replacement(w2, 2)]
prods += [(x, y) for x in w3 for y in w1]
pval = {p: mp.fprod(consts[x] for x in p) for p in prods}

natural = [QS2(Fr(1, 2)), B, -B, QS2(0, Fr(1, 2)), B * B, -B * B, QS2(2) - SQ2, QS2(1) - QS2(0, Fr(1, 2)),
           (QS2(2) - SQ2) / 4, (QS2(2) + SQ2) / 4, 2 * B, -B / 2, QS2(0, Fr(-1, 2)), QS2(-1), -SQ2 / 2 * B,
           QS2(17, -12), QS2(-16, 12)]
li4 = {z: mp.polylog(4, PE.DE.float_exact(z)) for z in natural}


def lll(vals, digits):
    C = mp.mpf(10) ** digits
    n = len(vals)
    rows = [[1 if i == j else 0 for j in range(n)] + [int(mp.nint(C * vals[i]))] for i in range(n)]
    Bm = flint.fmpz_mat(rows).lll()
    best = [int(Bm[0, j]) for j in range(n)]
    with mp.workdps(digits + 30):
        res = mp.fsum(c * v for c, v in zip(best, vals))
    return best, max(abs(c) for c in best), res


def qs2vals(vals):
    out = []
    for v in vals:
        out += [v, r2 * v]
    return out


base = qs2vals([pval[p] for p in prods])
kept = []
for z in natural:
    trial = base + qs2vals([li4[k] for k in kept] + [li4[z]])
    vec, mx, res = lll(trial, 1000)
    dep = mx < 10 ** (1000 / (2 * len(trial))) and (vec[-1] or vec[-2])
    say(f"  Li4({z}): {'dependent' if dep else 'kept'} (max coef {mx:.2e}) ({time.time()-T0:.0f}s)")
    if not dep:
        kept.append(z)
say("kept:", [str(z) for z in kept])

names = ["W4"] + [f"{s}{'*'.join(p)}" for p in prods for s in ("", "sqrt2*")] + \
        [f"{s}Li4({z})" for z in kept for s in ("", "sqrt2*")]
rels = {}
for m in "AB":
    for D in (700, 1150):
        vals = [W4[m]] + base + qs2vals([li4[z] for z in kept])
        vec, mx, res = lll(vals, D)
        if vec[0] < 0:
            vec = [-c for c in vec]
        rels[(m, D)] = vec
        say(f"method {m}, {D} digits: W4 coef {vec[0]}, max|coef| {mx:.3e}, residual {mp.nstr(res, 3)}, "
            f"random size ~1e{D/len(vals):.1f}")
same = len({tuple(v) for v in rels.values()}) == 1
say("relation identical in all four runs:", same)
vec = rels[("A", 1150)]

# exact form: W4 = -(1/vec0) * sum c_i x_i
c0 = vec[0]
S2 = sp.sqrt(2)
symb = {k: sp.Symbol(k) for k in ["lam", "L", "pi2", "Li2b", "z3", "Li3b", "Li3mb", "Li3h"]}
expr = 0
idx = 1
for p in prods:
    mono = sp.Mul(*[symb[x] for x in p])
    expr += sp.Rational(-vec[idx], c0) * mono + sp.Rational(-vec[idx + 1], c0) * S2 * mono
    idx += 2
li4sym = {}
for z in kept:
    s_ = sp.Symbol(f"Li4[{z}]")
    li4sym[z] = s_
    expr += sp.Rational(-vec[idx], c0) * s_ + sp.Rational(-vec[idx + 1], c0) * S2 * s_
    idx += 2
expr = sp.expand(expr)
say("\nW4 =", expr)
# numeric residual of the exact form at 1200 digits
subs = {symb[k]: consts[k] for k in symb}
f = sp.lambdify(list(symb.values()) + list(li4sym.values()), expr, modules="mpmath")
val = f(*[consts[k] for k in symb], *[li4[z] for z in kept])
for m in "AB":
    say(f"exact form - W4(method {m}) = {mp.nstr(val - W4[m], 3)}")
# M3 = W4 + W3' + T  (exact tails as sympy)
pickle.dump({"kept": kept, "vec": vec, "prods": prods}, open(os.path.join(OUTD, "m3_relation.pkl"), "wb"))
with open(os.path.join(OUTD, "m3_closed_form.out"), "w") as fh:
    fh.write("\n".join(lines) + "\n")
