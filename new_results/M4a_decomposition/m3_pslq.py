"""Integer-relation search for the weight-4 part of M3.

M3 = W4 + W3' + (pi^2 + 6 log^2 2)/24 with W3' and the tail EXACT (m3_tails_exact.py); hence
W4 = M3 - W3' - tail is known to the precision of M3 (method A, mk_numerics.py) and is a pure weight-4
iterated integral on P^1 minus {0, +-(1+sqrt2), +-(sqrt2-1), infinity}.

Basis over Q(sqrt2) (every element x enters as x and sqrt2*x):
  * the 22 weight-4 products of  {log2, L}, {pi^2, Li2(b)}, {z3, Li3(b), Li3(-b), Li3(1/sqrt2)}
  * Li4(z) for real S-units z < 1 of Q(sqrt2) with 1-z also an S-unit (bounded exponents), after removing
    candidates that are dependent (LLL, small coefficients) on the products and previously kept Li4's.
Lattice reduction: python-flint LLL on the standard integer-relation lattice.

Usage: python m3_pslq.py <M3 value file> <digits to use>
"""
import os, sys, pickle, time, itertools
import mpmath as mp
import flint
from fractions import Fraction as Fr
from qsqrt2 import QS2, SQ2, A, B, sunit_log
import polylog_exact as PE

HERE = os.path.dirname(os.path.abspath(__file__))
OUTD = os.path.join(HERE, "..", "outputs")
M3FILE = sys.argv[1]
D = int(sys.argv[2])
mp.mp.dps = D + 40
T0 = time.time()
lines = []


def say(*a):
    s = " ".join(str(x) for x in a)
    print(s, flush=True)
    lines.append(s)


M3 = mp.mpf(open(M3FILE).read().strip())
tails = pickle.load(open(os.path.join(OUTD, "m3_tails_exact.pkl"), "rb"))
W3p, T = tails["W3p"], tails["T"]
W4 = M3 - W3p.num() - T.num()
say(f"M3 from {os.path.basename(M3FILE)}; using {D} digits")
say("W4 = M3 - W3' - T =", mp.nstr(W4, 50))

r2 = mp.sqrt(2); lam = mp.log(2); L = mp.log(1 + r2); b = r2 - 1
pi2 = mp.pi ** 2; Li2b = mp.polylog(2, b)
z3 = mp.zeta(3); Li3b = mp.polylog(3, b); Li3mb = mp.polylog(3, -b); Li3h = mp.polylog(3, 1 / r2)
w1 = [("log2", lam), ("L", L)]
w2 = [("pi^2", pi2), ("Li2(b)", Li2b)]
w3 = [("z3", z3), ("Li3(b)", Li3b), ("Li3(-b)", Li3mb), ("Li3(1/sqrt2)", Li3h)]
prods = []
for c in itertools.combinations_with_replacement(range(2), 4):
    nm = "*".join(w1[i][0] for i in c); val = mp.fprod(w1[i][1] for i in c); prods.append((nm, val))
for (n2, v2) in w2:
    for c in itertools.combinations_with_replacement(range(2), 2):
        prods.append((n2 + "*" + "*".join(w1[i][0] for i in c), v2 * mp.fprod(w1[i][1] for i in c)))
for c in itertools.combinations_with_replacement(range(2), 2):
    prods.append(("*".join(w2[i][0] for i in c), mp.fprod(w2[i][1] for i in c)))
for (n3, v3) in w3:
    for (n1, v1) in w1:
        prods.append((n3 + "*" + n1, v3 * v1))
say(f"{len(prods)} weight-4 products")

# candidate Li4 arguments
cands = set()
for j in range(-6, 7):
    for m in range(-4, 5):
        for sg in (1, -1):
            z = QS2(sg) * SQ2 ** j * A ** m
            if not (z < QS2(1)) or z == QS2(0):
                continue
            try:
                sunit_log(QS2(1) - z)
            except ValueError:
                continue
            if abs(float(z)) > 1:      # use inversion: keep |z| <= 1 representatives
                continue
            cands.add(z)
cands = sorted(cands, key=lambda z: (abs(float(z)) > 0.9, -abs(float(z))))
say(f"{len(cands)} candidate Li4 arguments:", [str(z) for z in cands])


def qs2(lst):
    out = []
    for nm, v in lst:
        out += [(nm, v), ("sqrt2*" + nm, r2 * v)]
    return out


def lll(vals, digits):
    C = mp.mpf(10) ** digits
    n = len(vals)
    rows = [[1 if i == j else 0 for j in range(n)] + [int(mp.nint(C * vals[i]))] for i in range(n)]
    Bm = flint.fmpz_mat(rows).lll()
    best = [int(Bm[0, j]) for j in range(n)]
    with mp.workdps(digits + 20):
        resid = mp.fsum(c * v for c, v in zip(best, vals))
    return best, max(abs(c) for c in best), resid


# greedy selection of independent Li4's modulo products
base = qs2(prods)
kept = []
for z in cands:
    val = mp.polylog(4, PE.DE.float_exact(z))
    trial = base + qs2(kept) + qs2([(f"Li4({z})", val)])
    DS = min(D, 1000)
    vec, mx, res = lll([v for _, v in trial], DS)
    thresh = 10 ** (DS / (2 * len(trial)))         # random vectors have max coef ~ 10^(DS/n)
    dependent = mx < thresh and (vec[-1] != 0 or vec[-2] != 0)
    say(f"  Li4({z}): {'dependent' if dependent else 'kept'} (max coef {mx:.2e}) ({time.time()-T0:.0f}s)")
    if not dependent:
        kept.append((f"Li4({z})", val))
say(f"kept {len(kept)} independent Li4's:", [n for n, _ in kept])

# the search
full = [("W4", W4)] + qs2(prods) + qs2(kept)
vec, mx, res = lll([v for _, v in full], D)
thresh = 10 ** (D / (2 * len(full)))
say(f"\nLLL on W4 + {len(full)-1} basis numbers at {D} digits: max|coef| = {mx:.3e}, residual = {mp.nstr(res, 3)}"
    f"  (random size ~ 1e{D/len(full):.1f}; genuine-relation threshold {thresh:.1e})")
if vec[0] != 0 and mx < thresh:
    say("RELATION FOUND (W4 coefficient", vec[0], "):")
    say("   " + " ".join(f"{c:+d}*{nm}" for c, nm in zip(vec, [n for n, _ in full]) if c))
else:
    say("no relation with small coefficients involving W4 (W4 coefficient %d)" % vec[0])

with open(os.path.join(OUTD, f"m3_pslq_{D}.out"), "w") as fh:
    fh.write("\n".join(lines) + "\n")
