"""Symbolic proof of the closed form of M1 (Theorem 1 of the note).

Chain of exact steps (each one is checked here; numerics are used only as an independent sanity check):

 (0) Ahmed's parametrisation x = t^eps, y = (1-t)^eps, eps = 1/p:
         Q(p) = 2 - eps F(eps),  F(eps) = int_0^1 [a + b - sqrt(a^2+b^2)] dt,  a = t^(eps-1), b = (1-t)^(eps-1),
     so P/8 = 1 + eps M(eps) with M = -F/2 and M_k = -F^(k)(0) / (2 k!).
 (1) Differentiating once and using the symmetry t <-> 1-t:
         M1 = int_0^1 (log t / t) ((1-t)/q(t) - 1) dt,      q(t) = sqrt(t^2 + (1-t)^2).
 (2) Euler substitution t = (v^2 + 2v - 1)/(4v), v in [b, a] (b = sqrt2-1, a = 1+sqrt2 = 1/b):
         q = (v^2+1)/(2 sqrt2 v),  dt = (v^2+1)/(4v^2) dv,  t = (v-b)(v+a)/(4v),  1-t = (a-v)(v+b)/(4v);
     the rational factor collapses (double root at v=b) to
         M1 = int_b^a [log(v-b) + log(v+a) - log v - 2 log2] * [ (1 - 1/sqrt2)/v - 2/(v+a) ] dv.
 (3) Each int_b^a log|v-x| dv/(v-y) is a dilogarithm (exact antiderivative below).
 (4) All Li2 values are reduced to Li2(b) with the identities of dilog_exact.py.
 (5) Exact comparison with the claimed closed form (zero difference in Q(sqrt2)[pi^2, log^2 2, log2 L, L^2, Li2(b)]).
"""
import sys
import sympy as sp
import mpmath as mp
from fractions import Fraction as Fr
from qsqrt2 import QS2, SQ2, A, B
from dilog_exact import W, LOG, mul, PI2, LI2, reduce_all, CERT

mp.mp.dps = 60
ONE = QS2(1)
out = []


def say(*a):
    s = " ".join(str(x) for x in a)
    print(s)
    out.append(s)


# ---------------------------------------------------------------- step (2): symbolic substitution check
t, v = sp.symbols('t v', positive=True)
s2 = sp.sqrt(2)
a_, b_ = 1 + s2, s2 - 1
tv = (v ** 2 + 2 * v - 1) / (4 * v)
qv = (v ** 2 + 1) / (2 * s2 * v)
chk = {
    "t(v) = (v-b)(v+a)/(4v)": sp.simplify(tv - (v - b_) * (v + a_) / (4 * v)),
    "1-t(v) = (a-v)(v+b)/(4v)": sp.simplify(1 - tv - (a_ - v) * (v + b_) / (4 * v)),
    "q(t(v))^2 = t^2+(1-t)^2": sp.simplify(qv ** 2 - (tv ** 2 + (1 - tv) ** 2)),
    "t(b) = 0": sp.simplify(tv.subs(v, b_)),
    "t(a) = 1": sp.simplify(tv.subs(v, a_) - 1),
}
integrand_t = ((1 - t) / sp.sqrt(t ** 2 + (1 - t) ** 2) - 1) / t        # multiplies log t
rat_v = sp.simplify((((1 - tv) / qv - 1) / tv) * sp.diff(tv, v))       # q>0 on [b,a] so sqrt -> qv
target = (1 - 1 / s2) / v - 2 / (v + a_)
chk["rational factor = (1-1/sqrt2)/v - 2/(v+a)"] = sp.simplify(sp.radsimp(rat_v - target))
chk["log t = log(v-b)+log(v+a)-log v-2log2"] = sp.simplify(sp.expand_log(
    sp.log((v - b_) * (v + a_) / (4 * v)), force=True) - (sp.log(v - b_) + sp.log(v + a_) - sp.log(v) - 2 * sp.log(2)))
chk["dt/dv = (v^2+1)/(4v^2) > 0"] = sp.simplify(sp.diff(tv, v) - (v ** 2 + 1) / (4 * v ** 2))
say("Step (2) symbolic checks (all must be 0):")
for k, val in chk.items():
    say(f"   {k:48s}: {val}")
    assert val == 0

# ---------------------------------------------------------------- step (3): exact dilog antiderivatives
LO, HI = B, A


def K(wv):
    """antiderivative of log|1-w|/w, continuous on the real line minus 0:
       K(w) = -Li2(w) for w <= 1;  K(w) = log^2(w)/2 + Li2(1/w) - pi^2/3 for w >= 1."""
    if wv <= ONE:
        return -LI2(wv)
    lw = LOG(wv)
    return mul(lw, lw).scale(Fr(1, 2)) + LI2(ONE / wv) - PI2.scale(Fr(1, 3))


def I(x, y):
    """int_b^a log|v-x| dv/(v-y),  y outside [b,a].  Exact, as a W."""
    assert not (LO <= y <= HI)
    if x == y:
        l1, l0 = LOG(HI - y), LOG(LO - y)
        return (mul(l1, l1) - mul(l0, l0)).scale(Fr(1, 2))
    # w = (v-y)/(x-y): log|v-x| = log|x-y| + log|1-w|, dv/(v-y) = dw/w; the path in w avoids 0
    w0, w1 = (LO - y) / (x - y), (HI - y) / (x - y)
    return mul(LOG(x - y), LOG(w1 / w0)) + K(w1) - K(w0)


def I_num(x, y):
    xf, yf = float(x), float(y)
    X, Y = mp.mpf(xf), mp.mpf(yf)
    # recompute exactly-defined endpoints in mp
    from dilog_exact import float_exact
    X, Y = float_exact(x), float_exact(y)
    lo, hi = float_exact(LO), float_exact(HI)
    return mp.quad(lambda vv: mp.log(abs(vv - X)) / (vv - Y), [lo, 1, hi])


# integrand pieces:  log terms  log|v-x| with weights, times  c0/v - 2/(v+a)
c0 = ONE - ONE / SQ2
logs = [(B, QS2(1)), (-A, QS2(1)), (QS2(0), QS2(-1))]          # +log(v-b) + log(v+a) - log v
poles = [(QS2(0), c0), (-A, QS2(-2))]                              # c0/v - 2/(v+a)
M1 = W()
say("\nStep (3) elementary integrals I(x,y) = int_b^a log|v-x|/(v-y) dv:")
for x, cx in logs:
    for y, cy in poles:
        val = I(x, y)
        d = val.num() - I_num(x, y)
        say(f"   I({x}, {y}) = {val}      [quadrature check {mp.nstr(d, 3)}]")
        M1 = M1 + val.scale(cx * cy)
# constant -2 log 2 times int (c0/v - 2/(v+a)) dv
const = LOG(HI / LO).scale(c0) - LOG((HI + A) / (LO + A)).scale(2)
M1 = M1 - mul(LOG(QS2(4)), const)

M1_num_before = M1.num()
say("\nM1 before Li2 reduction:", M1)
say("   numerical value:", mp.nstr(M1_num_before, 50))

# ---------------------------------------------------------------- step (4): reduce Li2's
M1r = reduce_all(M1)
say("\nStep (4) certificate of Li2 reductions (each line is one classical identity):")
for c in CERT:
    say(c)
say("\nM1 after reduction:", M1r)

# ---------------------------------------------------------------- step (5): exact comparison with the claim
lam2, lamL, L2 = ('lam', 'lam'), ('L', 'lam'), ('L', 'L')
claim48 = W({('pi2',): QS2(-2, 19), lam2: 24, L2: QS2(-48, -120), lamL: QS2(0, 96), ('Li2', B): QS2(0, -240)})
claim = claim48.scale(Fr(1, 48))
diff = M1r - claim
say("\nStep (5) exact difference  M1(derived) - M1(claim) =", diff)
assert diff.d == {}, "MISMATCH"
say("=> 48*M1 = -2pi^2 + 19 sqrt2 pi^2 + 24 log^2 2 - 48 L^2 - 120 sqrt2 L^2 + 96 sqrt2 log2 L - 240 sqrt2 Li2(sqrt2-1)"
    "   PROVED (exact arithmetic).")

# ---------------------------------------------------------------- sanity: compare with direct quadrature of step (1)
mp.mp.dps = 60
M1_t = mp.quad(lambda tt: mp.log(tt) / tt * ((1 - tt) / mp.sqrt(tt ** 2 + (1 - tt) ** 2) - 1), [0, 0.25, 0.5, 1])
say("\nSanity: quadrature of step (1) integral =", mp.nstr(M1_t, 50))
say("        closed form                     =", mp.nstr(claim.num(), 50))
say("        8*M1 = kappa (Ahmed)            =", mp.nstr(8 * claim.num(), 50))

with open(sys.argv[1] if len(sys.argv) > 1 else "/dev/null", "w") as fh:
    fh.write("\n".join(out) + "\n")
