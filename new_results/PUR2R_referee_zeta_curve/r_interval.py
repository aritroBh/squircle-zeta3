"""Referee PUR2R: rigorous ball-arithmetic (python-flint arb) certificates for kappa^2 = 9/7.
(1) admissibility: Phi(pi/(2k)) > 1, so V < pi/(2k) and phi > 0 on [0, V];
(2) an enclosure [Vlo, Vhi] of V;
(3) the E2(i) condition  kappa^2 <= 3 e^v cos^3(kappa v)  for ALL v in [0, Vhi], by subdivision
    (no unimodality needed); also at the endpoint ball.
(4) kappa_max bracket: q(kappa) = kappa e^{pi/(2kappa)} - 2 - kappa^2 changes sign."""
from flint import arb, ctx
ctx.prec = 200

k2 = arb(9) / 7
k = k2.sqrt()
pi = arb.pi()

def Phi(v):
    return (v.exp() * ((k * v).cos() + k * (k * v).sin()) - 1) / (1 + k2)

vtop = pi / (2 * k)
val = Phi(vtop)
print("(1) Phi(pi/(2k)) =", val.str(20), " > 1 certified:", bool(val > 1))
print("    (Phi' = e^v cos(kv) > 0 on [0, pi/(2k)), so Phi is a strictly increasing bijection [0,V] -> [0,1].)")

# bisection with exact-ish rationals for V
lo, hi = arb(0), arb("1.5")
for _ in range(190):
    mid = (lo + hi) / 2
    mid = arb(mid.mid())
    if Phi(mid) < 1:
        lo = mid
    elif Phi(mid) > 1:
        hi = mid
    else:
        break
print("(2) Phi(Vlo) < 1:", bool(Phi(lo) < 1), " Phi(Vhi) > 1:", bool(Phi(hi) > 1), " Vhi < pi/(2k):", bool(hi < vtop))
print("    V in [%s, %s]" % (lo.str(40, radius=False), hi.str(40, radius=False)))
Vball = arb((lo + hi) / 2, ((hi - lo) / 2).upper() if hasattr((hi - lo), 'upper') else (hi - lo))
Vball = lo.union(hi)
g = lambda v: 3 * v.exp() * (k * v).cos() ** 3
gV = g(Vball)
print("    3 e^V cos^3(kV) in", gV.str(25), "; > 9/7 certified:", bool(gV > k2))
print("    phi(V) = e^V cos(kV) in", (Vball.exp() * (k * Vball).cos()).str(20))

# (3) subdivision of [0, Vhi]
n = 4000
ok = True
worst = None
for i in range(n):
    a = hi * i / n
    b = hi * (i + 1) / n
    I = a.union(b)
    gi = g(I)
    if not (gi > k2):
        ok = False
        print("    FAIL on", I.str(10))
        break
    lowb = gi.lower()
    if worst is None or lowb < worst:
        worst = lowb
print("(3) kappa^2 < 3 e^v cos^3(kappa v) for all v in [0, Vhi]:", ok, "; certified lower bound of min:", worst.str(15))

# (4) kappa_max
def q(kk):
    return kk * (pi / (2 * kk)).exp() - 2 - kk * kk
a, b = arb("1.50781081090595335"), arb("1.50781081090595336")
print("(4) q(1.50781081090595335) =", q(a).str(5), " >0:", bool(q(a) > 0), "; q(...36) =", q(b).str(5), " <0:", bool(q(b) < 0))
print("    q(sqrt(9/7)) =", q(k).str(15), " > 0:", bool(q(k) > 0))

# E2(i) as an inequality chain needs only (3); also report the k'(p) bound margin at V
