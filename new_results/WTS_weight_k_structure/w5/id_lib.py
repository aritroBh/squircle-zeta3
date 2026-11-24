# Helper library for the identification step (python-flint arb constants + exact-Gram LLL relation search).
# No mpmath here; precision is flint ctx.prec (bits) and every verification prints it.
import sys, os, json, time
sys.dont_write_bytecode = True
from flint import arb, ctx, fmpz, fmpz_mat

def setprec(D):
    ctx.prec = int(D * 3.3219281) + 256

def S(x): return arb(x)

class Consts:
    """weight-graded constants at the current precision (call after setprec)."""
    def __init__(self):
        one = arb(1)
        self.sq2 = arb(2).sqrt()
        self.a = 1 + self.sq2; self.b = self.sq2 - 1
        self.pi = arb.pi()
        self.lam = arb(2).log()
        self.L = self.a.log()
        self.z3 = arb(3).zeta(); self.z5 = arb(5).zeta()
        half = one / 2
        self.Li4h = half.polylog(4); self.Li5h = half.polylog(5)
    def Li(self, s, z): return arb(z).polylog(s)

def to_int(x, D):
    return (x * arb(10) ** D).floor().unique_fmpz()

def lll_relation(target_and_basis, D, extra_scale=0):
    """integer relation search among vals (list of arb) using LLL (gram='exact') on [I | round(10^D x)].
    returns (coeffs list of python ints, last entry (residual integer), max|coef|)."""
    vals = target_and_basis
    n = len(vals)
    M = fmpz_mat(n, n + 1)
    for i in range(n):
        M[i, i] = 1
        M[i, n] = to_int(vals[i], D)
    R = M.lll(gram="exact")
    best = None
    for r in range(R.nrows()):
        row = [int(R[r, c]) for c in range(n + 1)]
        mx = max(abs(x) for x in row[:n])
        if best is None or mx < best[2]:
            best = (row[:n], row[n], mx)
    return best, [[int(R[r, c]) for c in range(n + 1)] for r in range(min(3, R.nrows()))]

def residual(coeffs, vals):
    s = arb(0)
    for c, v in zip(coeffs, vals): s += c * v
    return s

def load_cm(fn="results/cm_D1500_m10.json"):
    R = json.load(open(fn))
    return R

# ---------------------------------------------------------------------------------------------
# weight-5 atoms (analogue of the weight-4 atom list of the paper / sibling w4 step)
def atoms5(C, with_sqrt2=True, li4=True, li5=True, li3=True, li2=True, chi8=False):
    lam, L, pi2 = C.lam, C.L, C.pi ** 2
    b = C.b; s2 = C.sq2
    Li2b = C.Li(2, b)
    Li3 = {"b": C.Li(3, b), "-b": C.Li(3, -b), "1/s2": C.Li(3, 1 / s2)}
    pts = {"1/2": arb(1) / 2, "b": b, "-b": -b, "1/s2": 1 / s2, "-b2": -b * b, "2-s2": 2 - s2, "1-1/s2": 1 - 1 / s2, "1/s2-1": 1 / s2 - 1}
    out = []
    # (a) lam^a L^(5-a)
    for a_ in range(6): out.append((f"lam^{a_} L^{5-a_}", lam ** a_ * L ** (5 - a_)))
    for a_ in range(4): out.append((f"pi2 lam^{a_} L^{3-a_}", pi2 * lam ** a_ * L ** (3 - a_)))
    out.append(("pi4 lam", pi2 ** 2 * lam)); out.append(("pi4 L", pi2 ** 2 * L))
    if li2:
        for a_ in range(4): out.append((f"Li2b lam^{a_} L^{3-a_}", Li2b * lam ** a_ * L ** (3 - a_)))
        out.append(("Li2b pi2 lam", Li2b * pi2 * lam)); out.append(("Li2b pi2 L", Li2b * pi2 * L))
        out.append(("Li2b^2 lam", Li2b ** 2 * lam)); out.append(("Li2b^2 L", Li2b ** 2 * L))
    out.append(("z3 lam^2", C.z3 * lam ** 2)); out.append(("z3 lam L", C.z3 * lam * L)); out.append(("z3 L^2", C.z3 * L ** 2))
    out.append(("z3 pi2", C.z3 * pi2))
    if li2: out.append(("z3 Li2b", C.z3 * Li2b))
    if li3:
        for nm, v in Li3.items():
            out.append((f"Li3({nm}) lam^2", v * lam ** 2)); out.append((f"Li3({nm}) lam L", v * lam * L)); out.append((f"Li3({nm}) L^2", v * L ** 2))
            out.append((f"Li3({nm}) pi2", v * pi2))
            if li2: out.append((f"Li3({nm}) Li2b", v * Li2b))
    if li4:
        for nm, z in pts.items():
            if li4 == "half" and nm != "1/2": continue
            v = z.polylog(4)
            out.append((f"Li4({nm}) lam", v * lam)); out.append((f"Li4({nm}) L", v * L))
    if li5:
        out.append(("z5", C.z5))
        for nm, z in pts.items():
            if li5 == "half" and nm != "1/2": continue
            out.append((f"Li5({nm})", z.polylog(5)))
    if chi8:
        tot = arb(0)
        for r, ch in ((1, 1), (3, -1), (5, -1), (7, 1)):
            tot += ch * arb(5).zeta(arb(r) / 8)
        out.append(("L(chi8,5)", tot / arb(8) ** 5))
    if with_sqrt2:
        out = out + [("sq2*" + n, v * s2) for n, v in out]
    return out
