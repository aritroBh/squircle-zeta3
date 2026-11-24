"""Attempted symbolic evaluation of the pure weight-3 part W3 of M2.

From m2_decomposition_and_tail.py (exact):
    M2 = W3 + (pi^2 + 6 log^2 2)/24,
    W3 = int_b^a [ log^2 t * ( -1/(v+a) + (4-3 sqrt2)/(8v) ) + (sqrt2/8) log t log(1-t) / v ] dv,
    log t    = log(v-b) + log(v+a) - log v - 2 log2,
    log(1-t) = log(a-v) + log(v+b) - log v - 2 log2.
Every product log|v-x| log|v-y| dv/(v-z) is integrated exactly on [b,1] and [1,a] (polylog_exact.py),
then Li2's are reduced with dilog_exact.py and Li3's by Gaussian elimination over classical Li3 relations.
"""
import sys
import mpmath as mp
from fractions import Fraction as Fr
from qsqrt2 import QS2, SQ2, A, B
import polylog_exact as PE
from polylog_exact import P, LOG, LAM, LL, PI2, Z3, LI3

mp.mp.dps = 60
ONE, ZERO = QS2(1), QS2(0)
out = []


def say(*a):
    s = " ".join(str(x) for x in a)
    print(s, flush=True)
    out.append(s)


say("Antiderivative checks (sympy, must be 0):", PE.verify_antiderivatives())

# logs as {letter: coeff}, const P
LT = ({B: ONE, -A: ONE, ZERO: -ONE}, LAM * (-2))
L1T = ({A: ONE, -B: ONE, ZERO: -ONE}, LAM * (-2))
SUB = [(B, ONE), (ONE, A)]


def int_prod(l1, l2, z, coef):
    """coef * int_b^a l1(v) l2(v) dv/(v-z) as P."""
    (c1, k1), (c2, k2) = l1, l2
    tot = P()
    for v0, v1 in SUB:
        for x, cx in c1.items():
            for y, cy in c2.items():
                tot = tot + PE.int_loglog(x, y, z, v0, v1) * (cx * cy)
        for x, cx in c1.items():
            tot = tot + PE.int_log(x, z, v0, v1) * cx * k2
        for y, cy in c2.items():
            tot = tot + PE.int_log(y, z, v0, v1) * cy * k1
        tot = tot + PE.int_one(z, v0, v1) * k1 * k2
    return tot * coef


# --- numerical spot checks of the basic integrals
def numint(f, v0, v1):
    return mp.quad(f, [PE.DE.float_exact(v0), PE.DE.float_exact(v1)])


fe = PE.DE.float_exact
for (x, y, z) in [(B, -A, ZERO), (B, A, ZERO), (ZERO, -B, -A), (A, A, ZERO), (B, B, -A), (-A, ZERO, -A), (A, -B, ZERO)]:
    for v0, v1 in SUB:
        ex = PE.int_loglog(x, y, z, v0, v1).num()
        nu = numint(lambda vv: mp.log(abs(vv - fe(x))) * mp.log(abs(vv - fe(y))) / (vv - fe(z)), v0, v1)
        assert abs(ex - nu) < mp.mpf(10) ** -40, (x, y, z, v0, v1, ex, nu)
say("basic weight-3 integrals: spot checks against quadrature OK")

c_a = QS2(-1)                       # coefficient of dv/(v+a) on log^2 t
c_0 = (QS2(4) - SQ2 * 3) / 8        # coefficient of dv/v   on log^2 t
c_x = SQ2 / 8                       # coefficient of dv/v   on log t log(1-t)
W3 = int_prod(LT, LT, -A, c_a) + int_prod(LT, LT, ZERO, c_0) + int_prod(LT, L1T, ZERO, c_x)

mp.mp.dps = 60
r2 = mp.sqrt(2); aN, bN = 1 + r2, r2 - 1
lt = lambda x: mp.log((x - bN) * (x + aN) / (4 * x))
l1t = lambda x: mp.log((aN - x) * (x + bN) / (4 * x))
W3_num = mp.quad(lambda x: lt(x) ** 2 * (-1 / (x + aN) + (4 - 3 * r2) / (8 * x)) + r2 / 8 * lt(x) * l1t(x) / x,
                 [bN, 1, aN])
say("W3 exact-expression numeric =", mp.nstr(W3.num(), 45))
say("W3 quadrature               =", mp.nstr(W3_num, 45))
assert abs(W3.num() - W3_num) < mp.mpf(10) ** -40
say("Li2 atoms:", sorted(W3.atoms('Li2'), key=float))
say("Li3 atoms:", sorted(W3.atoms('Li3'), key=float))

W3b = PE.reduce_li2_all(W3)
say("after Li2 reduction, Li2 atoms:", sorted(W3b.atoms('Li2'), key=float))
basis = [B, -B, QS2(0, Fr(1, 2))]          # Li3(b), Li3(-b), Li3(1/sqrt2)
import kummer
KREL = []
for (x, y) in kummer.search():
    KREL.append((f"Kummer (K) at (x,y)=({x}, {y})", kummer.kummer_relation(x, y)))
say(f"{len(KREL)} Spence-Kummer instances with S-unit arguments added")
try:
    W3c, used, rk = PE.solve_li3_linear(W3b, basis, extra_relations=KREL)
    say(f"Li3 reduction: rank {rk}; relations used with multipliers:")
    for mu, nm in used:
        say(f"    {mu} x [{nm}]")
    say("W3 reduced =", W3c)
    with mp.workdps(60):
        say("   numeric check:", mp.nstr(W3c.num() - W3_num, 5))
except RuntimeError as e:
    say("Li3 reduction FAILED:", e)
    W3c = None

with open(sys.argv[1] if len(sys.argv) > 1 else "/dev/null", "w") as fh:
    fh.write("\n".join(out) + "\n")

# ---------------------------------------------------------------------------------------------
# Final exact comparison:  M2 = W3 + T  with  T = (pi^2 + 6 log^2 2)/24  (m2_decomposition_and_tail.py)
if W3c is not None:
    T = (PI2 + LAM * LAM * 6) * Fr(1, 24)
    M2 = W3c + T
    say("\nM2 = W3 + T =", M2)
    from polylog_exact import LI2
    r = SQ2
    Q = lambda p, q=0: QS2(Fr(p), Fr(q))
    # claimed closed form, real-Li3 form ("form D" of research/05_final_verify, with Li3(1/2) -> zeta(3)):
    claim = (PI2 * Fr(1, 24) + LAM * LAM * Fr(1, 4)
             + Z3 * QS2(Fr(-19, 16), Fr(189, 64))
             + LI3(B) * QS2(11, 15) + LI3(-B) * QS2(11, -15) + LI3(QS2(0, Fr(1, 2))) * QS2(0, -27)
             + LAM * PI2 * QS2(Fr(-1, 24), Fr(-65, 48)) + LL * PI2 * QS2(Fr(-7, 12), Fr(15, 8))
             + LAM * LI2(B) * QS2(0, Fr(-7, 2)) + LL * LI2(B) * 12
             + LAM * LAM * LAM * QS2(Fr(1, 6), Fr(9, 16)) + LAM * LAM * LL * QS2(0, Fr(11, 8))
             + LAM * LL * LL * QS2(-1, Fr(-7, 4)) + LL * LL * LL * QS2(Fr(25, 6), Fr(-7, 3)))
    diff = M2 - claim
    say("M2(derived) - M2(claimed real-Li3 form) =", diff)
    assert not diff.d
    say("=> THEOREM: M2 = (pi^2+6log^2 2)/24 + (-19/16+189 sqrt2/64) z3 + (11+15sqrt2) Li3(b) + (11-15sqrt2) Li3(-b)"
        " - 27 sqrt2 Li3(1/sqrt2) + products   PROVED (exact arithmetic, classical + Kummer relations).")
    with mp.workdps(60):
        say("   numeric value of the proved expression:", mp.nstr(M2.num(), 50))

    # per-atom table of the Li3 reductions (for the paper)
    say("\nPer-atom Li3 reductions (each in the basis Li3(b), Li3(-b), Li3(1/sqrt2), z3, products):")
    for z in sorted(W3b.atoms('Li3'), key=float):
        if z in basis:
            continue
        try:
            red, used_z, _ = PE.solve_li3_linear(LI3(z), basis, extra_relations=KREL)
            kinds = sorted({nm.split(' (')[0] for _, nm in used_z})
            say(f"  Li3({z}) = {red}    [uses: {', '.join(kinds)}]")
        except RuntimeError as e:
            say(f"  Li3({z}): not individually reducible ({str(e)[-60:]})")

with open(sys.argv[1] if len(sys.argv) > 1 else "/dev/null", "w") as fh:
    fh.write("\n".join(out) + "\n")
