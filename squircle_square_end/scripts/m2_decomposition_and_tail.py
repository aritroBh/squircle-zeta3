"""M2: exact decomposition into a pure weight-3 integral over the real alphabet {0, +-a, +-b}
plus a weight-2 integral T, and exact evaluation  T = (pi^2 + 6 log^2 2)/24   (Theorem 2(a) of the note).

Steps (all symbolic; numerics only as sanity checks):
 (0) M2 = -F''(0)/4 and, using t <-> 1-t,
       F''(0) = 2 int_0^1 log^2 t (1 - 2c + c^3) dt/t + 2 int_0^1 log t log(1-t) t(1-t)/q^3 dt,
       c = (1-t)/q,  q = sqrt(t^2+(1-t)^2).
 (1) Euler substitution t = (v^2+2v-1)/(4v):  F''(0) = 2 int_b^a [log^2 t R1(v) + log t log(1-t) R2(v)] dv with
       R1 = 2/(v+a) + (3 sqrt2 - 4)/(4v) + D(v),   R2 = -sqrt2/(4v) - D(v),   D(v) = -2 sqrt2 v/(v^2+1)^2.
     NO simple poles at v = +-i; the double poles of R1 and R2 cancel in pairs.
 (2) Hence  M2 = W3 + T with
       W3 = -1/2 int_b^a [ log^2 t (2/(v+a) + (3sqrt2-4)/(4v)) - (sqrt2/4) log t log(1-t)/v ] dv   (pure weight 3)
       T  = -1/2 int_b^a log t * log(t/(1-t)) D(v) dv.
     D dv = dg with g = sqrt2/(v^2+1) - a/2 (g(b) = 0, g(a) = -1); integrating by parts (boundary terms vanish)
       T = 1/2 int_b^a g(v) [ (2 log t - log(1-t)) dt/t + log t dt/(1-t) ].
 (3) The rational functions g dt/t, g dt/(1-t) have only real simple poles; T is evaluated exactly
     with w2_integrals.py and the Li2's are reduced with dilog_exact.py.
"""
import sys
import sympy as sp
import mpmath as mp
from fractions import Fraction as Fr
from qsqrt2 import QS2, SQ2, A, B
from dilog_exact import W, LOG, mul, PI2, reduce_all, CERT
from w2_integrals import integrate, partial_fractions, LOG_T, LOG_1MT, v

out = []


def say(*a):
    s = " ".join(str(x) for x in a)
    print(s)
    out.append(s)


s2 = sp.sqrt(2)
a_, b_ = 1 + s2, s2 - 1
tv = (v ** 2 + 2 * v - 1) / (4 * v)
qv = (v ** 2 + 1) / (2 * s2 * v)
dt = sp.diff(tv, v)
c = (1 - tv) / qv
R1 = sp.together((1 - 2 * c + c ** 3) / tv * dt)
R2 = sp.together(tv * (1 - tv) / qv ** 3 * dt)
D = -2 * s2 * v / (v ** 2 + 1) ** 2
R1_claim = 2 / (v + a_) + (3 * s2 - 4) / (4 * v) + D
R2_claim = -s2 / (4 * v) - D
say("Step (1) symbolic checks (must be 0):")
for nm, e in (("R1 - [2/(v+a) + (3sqrt2-4)/(4v) + D]", R1 - R1_claim), ("R2 - [-sqrt2/(4v) - D]", R2 - R2_claim)):
    z = sp.simplify(sp.radsimp(sp.together(e)))
    say(f"   {nm:42s}: {z}")
    assert z == 0
g = s2 / (v ** 2 + 1) - a_ / 2
chk = {"dg/dv - D": sp.simplify(sp.diff(g, v) - D), "g(b)": sp.simplify(sp.radsimp(g.subs(v, b_))),
       "g(a) + 1": sp.simplify(sp.radsimp(g.subs(v, a_) + 1))}
for k, z in chk.items():
    say(f"   {k:42s}: {z}")
    assert z == 0

# ---------------------------------------------------------------- numerical sanity of the decomposition
mp.mp.dps = 40
r2 = mp.sqrt(2); aN, bN = 1 + r2, r2 - 1; lam = mp.log(2); LN = mp.log(aN)
def lt(x): return mp.log((x - bN) * (x + aN) / (4 * x))
def l1t(x): return mp.log((aN - x) * (x + bN) / (4 * x))
W3_num = -mp.quad(lambda x: lt(x) ** 2 * (2 / (x + aN) + (3 * r2 - 4) / (4 * x)) - r2 / 4 * lt(x) * l1t(x) / x,
                  [bN, 1, aN]) / 2
T_num = -mp.quad(lambda x: lt(x) * (lt(x) - l1t(x)) * (-2 * r2 * x / (x ** 2 + 1) ** 2), [bN, 1, aN]) / 2
M2_ref = mp.mpf("-0.022161911087269399423754937739791514977955046328743456207867634148706991765943994309262616296552742025630667691293284821742851205290465367552390676970896303323914341510399819277290948522489025009550006567579318412923373908941089144005488472431894895856637895216705883246453272686942878734765615672745362707573251174787438333759515199438421146347275664477270882574501108359203680625805755453910796094470060004588321440880533568814598619340343848014343853180251535074928278627991165064660975431526453403166772142422885348586980752848224717754393592376755320733766119744721911020740710518875060960812042588017769053777847325281507975172567795008212595617995266518757504473579523178008146096578787402")
say("\nSanity (numerics): W3 =", mp.nstr(W3_num, 35), "  T =", mp.nstr(T_num, 35))
say("   W3 + T - M2(reference, 690 digits) =", mp.nstr(W3_num + T_num - M2_ref, 5))

# ---------------------------------------------------------------- step (3): exact evaluation of T
# T = 1/2 int g [ (2 log t - log(1-t)) dt/t + log t dt/(1-t) ]
Gt = sp.together(g * dt / tv)            # g dt/t
G1t = sp.together(g * dt / (1 - tv))     # g dt/(1-t)
for nm, G in (("g dt/t", Gt), ("g dt/(1-t)", G1t)):
    poly, res = partial_fractions(G)
    say(f"   partial fractions of {nm}: poly={poly}, residues={res}")


def lin(l1, c1, l2=None, c2=None):
    """c1*Lambda1 + c2*Lambda2 as a Lambda."""
    co = {k: v_ * c1 for k, v_ in l1[0].items()}
    const = l1[1].scale(c1)
    if l2 is not None:
        for k, v_ in l2[0].items():
            co[k] = co.get(k, QS2(0)) + v_ * c2
        const = const + l2[1].scale(c2)
    co = {k: v_ for k, v_ in co.items() if v_ != QS2(0)}
    return co, const


T = (integrate(lin(LOG_T, QS2(2), LOG_1MT, QS2(-1)), Gt) + integrate(LOG_T, G1t)).scale(Fr(1, 2))
say("\nT before Li2 reduction:", T)
say("   numerical:", mp.nstr(T.num(), 35))
Tr = reduce_all(T)
say("\nLi2 reduction certificate:")
for line in CERT:
    say(line)
say("\nT after reduction:", Tr)
tail = (PI2 + mul(LOG(QS2(2)), LOG(QS2(2))).scale(6)).scale(Fr(1, 24))
d = Tr - tail
say("T - (pi^2 + 6 log^2 2)/24 =", d)
assert d.d == {}
say("=> T = (pi^2 + 6 log^2 2)/24 = (zeta(2) - Li2(1/2))/2   PROVED (exact arithmetic).")
say("   Hence M2 = W3 + (pi^2+6 log^2 2)/24 with W3 the pure weight-3 integral above.")

with open(sys.argv[1] if len(sys.argv) > 1 else "/dev/null", "w") as fh:
    fh.write("\n".join(out) + "\n")
