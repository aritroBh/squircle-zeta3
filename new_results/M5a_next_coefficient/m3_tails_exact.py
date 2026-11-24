"""M3 = W4 + W3' + T3 (ibp_reduce.py): exact evaluation of the weight-3 piece W3' and the weight-2 piece T3.
The weight-4 piece W4 is then M3 - W3' - T3 (numerically, from the 1200-digit M3) and goes to PSLQ.

Writes ../outputs/m3_tails_exact.out and a pickle with the exact tails for the PSLQ script.
"""
import sys, pickle, os
import sympy as sp
import mpmath as mp
from fractions import Fraction as Fr
from qsqrt2 import QS2, SQ2, A, B
import polylog_exact as PE
from polylog_exact import P, LOG, LAM, LL, PI2, Z3, LI3, LI2
from w2_integrals import integrate as w2_integrate, partial_fractions, LOG_T, LOG_1MT, to_qs2
from ibp_reduce import reduce_k
from mk_integrands import v
import kummer

HERE = os.path.dirname(os.path.abspath(__file__))
mp.mp.dps = 60
ONE, ZERO = QS2(1), QS2(0)
out = []


def say(*a):
    s = " ".join(str(x) for x in a)
    print(s, flush=True)
    out.append(s)


K = int(sys.argv[1]) if len(sys.argv) > 1 else 3
pieces = reduce_k(K, verbose=False)
say(f"k = {K}: pieces", {k: str(v_) for k, v_ in pieces.items()})

LT = ({B: ONE, -A: ONE, ZERO: -ONE}, LAM * (-2))
L1T = ({A: ONE, -B: ONE, ZERO: -ONE}, LAM * (-2))
SUB = [(B, ONE), (ONE, A)]


def residues(R):
    """sum of terms c/(v-p) (output of ibp_reduce) -> {pole (QS2): residue (QS2)}; exact."""
    res = {}
    for term in sp.Add.make_args(sp.expand(R)):
        num, den = sp.fraction(sp.together(term))
        dp = sp.Poly(sp.expand(den), v)
        assert dp.degree() == 1 and not sp.expand(num).has(v), term
        lc = dp.LC()
        p = to_qs2(sp.radsimp(-dp.nth(0) / lc))
        c = to_qs2(sp.radsimp(num / lc))
        res[p] = res.get(p, QS2(0)) + c
    # exactness check
    back = sum(c.sympy() / (v - p.sympy()) for p, c in res.items())
    assert sp.simplify(sp.radsimp(sp.together(back - R))) == 0
    return res


def as_clean(R):
    return sum(c.sympy() / (v - p.sympy()) for p, c in residues(R).items())


pieces = {k: (as_clean(R) if k != (0, 0) else R) for k, R in pieces.items()}


def std_prod(l1, l2, z, v0, v1):
    """int_{v0}^{v1} l1 l2 dv/(v-z), z outside [v0, v1]."""
    (c1, k1), (c2, k2) = l1, l2
    tot = P()
    for x, cx in c1.items():
        for y, cy in c2.items():
            tot = tot + PE.int_loglog(x, y, z, v0, v1) * (cx * cy)
    for x, cx in c1.items():
        tot = tot + PE.int_log(x, z, v0, v1) * cx * k2
    for y, cy in c2.items():
        tot = tot + PE.int_log(y, z, v0, v1) * cy * k1
    tot = tot + PE.int_one(z, v0, v1) * k1 * k2
    return tot


def lam_at(lam, e):
    co, k = lam
    val = k
    for x, c in co.items():
        val = val + LOG(e - x) * c
    return val


def reg_prod(l1, l2, e, v0, v1):
    """int l1 l2 dv/(v-e) with e an endpoint where l1 and l2 both vanish."""
    assert not lam_at(l1, e).d and not lam_at(l2, e).d
    tot = P()
    for x, cx in l1[0].items():
        for y, cy in l2[0].items():
            tot = tot + PE.int_reg2(x, y, e, v0, v1) * (cx * cy)
    return tot


def weight3_piece(i, j, R):
    l_a = LT
    l_b = LT if i == 2 else L1T
    assert (i, j) in ((2, 0), (1, 1))
    tot = P()
    for z, r in residues(R).items():
        for v0, v1 in SUB:
            if not (v0 <= z <= v1):
                tot = tot + std_prod(l_a, l_b, z, v0, v1) * r
            else:
                assert z in (v0, v1), z
                tot = tot + reg_prod(l_a, l_b, z, v0, v1) * r
    return tot


# ---------------------------------------------------------------- numerics for checks
r2 = mp.sqrt(2); aN, bN = 1 + r2, r2 - 1
lt = lambda x: mp.log((x - bN) * (x + aN) / (4 * x))
l1t = lambda x: mp.log((aN - x) * (x + bN) / (4 * x))
fe = PE.DE.float_exact


def quad_piece(keys):
    fs = [(i, j, sp.lambdify(v, pieces[(i, j)], modules="mpmath")) for (i, j) in keys]
    return mp.quad(lambda x: mp.fsum(f(x) * lt(x) ** i * l1t(x) ** j for i, j, f in fs), [bN, 1, aN])


# weight 2
T3 = PE.from_W(w2_integrate(LOG_T, pieces[(1, 0)]))
say("\nweight-2 piece before reduction:", T3)
T3r = PE.reduce_li2_all(T3)
say("weight-2 piece T =", T3r, "   [quadrature check:", mp.nstr(T3r.num() - quad_piece([(1, 0)]), 3), "]")
tail2 = (PI2 + LAM * LAM * 6) * Fr(1, 24)
say("T - (pi^2+6log^2 2)/24 =", T3r - tail2)

# weight 3
W3p = weight3_piece(2, 0, pieces[(2, 0)]) + weight3_piece(1, 1, pieces[(1, 1)])
q3 = quad_piece([(2, 0), (1, 1)])
say("\nweight-3 piece: exact-expression numeric =", mp.nstr(W3p.num(), 40), "  quadrature =", mp.nstr(q3, 40))
assert abs(W3p.num() - q3) < mp.mpf(10) ** -40
W3pb = PE.reduce_li2_all(W3p)
KREL = [(f"Kummer (K) at (x,y)=({x}, {y})", kummer.kummer_relation(x, y)) for (x, y) in kummer.search()]
basis = [B, -B, QS2(0, Fr(1, 2))]
W3pc, used, rk = PE.solve_li3_linear(W3pb, basis, extra_relations=KREL)
say(f"Li3 reduction: rank {rk}; {len(used)} relations used:")
for mu, nm in used:
    say(f"    {mu} x [{nm}]")
say("weight-3 piece W3' =", W3pc)
with mp.workdps(60):
    say("   check vs quadrature:", mp.nstr(W3pc.num() - q3, 3))

with open(os.path.join(HERE, f"m{K}_tails_exact.pkl"), "wb") as fh:
    pickle.dump({"T": T3r, "W3p": W3pc}, fh)
with open(sys.argv[2] if len(sys.argv) > 2 else os.path.join(HERE, f"m{K}_tails_exact.out"), "w") as fh:
    fh.write("\n".join(out) + "\n")
