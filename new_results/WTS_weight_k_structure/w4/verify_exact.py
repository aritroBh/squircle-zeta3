"""verify_exact.py <dir_with_red_11_to_15>
EXACT (Fraction arithmetic) re-derivation, from the pickled final pieces, of the k-dependence of the weight-4 pieces:
  W_4(k) = sum_{poles r in {-a,0,-b,a,b}}  [ cJ_r(k) * J_r + cK_r(k) * K_r ],   J_r = int dv/(v-r) log^2 t log(1-t) (key (2,1)),
                                                                                 K_r = int dv/(v-r) log^3 t        (key (3,0)).
Coefficients cJ_r(k), cK_r(k) in Q(sqrt2) are read for k=3..15 (k=3..10 from the author's red_k.pkl, k=11..15 from a fresh
re-run of the exact reducer done by me), fitted exactly by a quadratic through k=4,5,6 and tested EXACTLY at every other k,
and compared with the table claimed in identify_w4.md.  Also tests that no piece of weight 4 has an i-pole or constant."""
import sys, os, pickle
sys.dont_write_bytecode = True
from fractions import Fraction as Fr
W4 = "<project>/research/08_theory/WTS_weight_k_structure/w4"
sys.path.insert(0, W4)
import fastred

extra = sys.argv[1] if len(sys.argv) > 1 else None
POLE = {(-1, -1): "-a", (0, 0): "0", (1, -1): "-b", (1, 1): "a", (-1, 1): "b"}


def load(k):
    fn = os.path.join(W4, f"red_{k}.pkl")
    if not os.path.exists(fn):
        fn = os.path.join(extra, f"red_{k}.pkl")
    return pickle.load(open(fn, "rb"))


def table(k):
    red = load(k)
    out = {}
    for (i, j), P in red.items():
        if i + j != 3:
            continue
        real, ip, im, const = P
        assert not ip and not im and const.is_zero(), ("non-final weight-4 piece", k, (i, j))
        assert (i, j) in ((3, 0), (2, 1))
        for r, c in real.items():
            assert r.p2 == 0 and r.p3 == 0 and c.p2 == 0 and c.p3 == 0
            name = ("K" if (i, j) == (3, 0) else "J") + "_" + POLE[(r.p0, r.p1)]
            out[name] = (Fr(c.p0), Fr(c.p1))          # coefficient p0 + p1*sqrt2
    return out


KS = list(range(3, 16))
T = {k: table(k) for k in KS}
names = sorted(set().union(*[set(T[k]) for k in KS]))
print("constants occurring for k=3..15:", names)
print("k : which constants occur")
for k in KS:
    print(f"  k={k:2d}: {sorted(T[k])}")


def quad_fit(ys):               # exact quadratic through k=4,5,6
    (k1, k2, k3), (y1, y2, y3) = (4, 5, 6), ys
    c2 = ((y3 - y2) / (k3 - k2) - (y2 - y1) / (k2 - k1)) / (k3 - k1)
    c1 = (y2 - y1) / (k2 - k1) - c2 * (k1 + k2)
    c0 = y1 - c1 * k1 - c2 * k1 * k1
    return c0, c1, c2


ALL_OK = True
FIT = {}
for nm in names:
    fit = []
    for comp in (0, 1):
        ys = [T[k].get(nm, (Fr(0), Fr(0)))[comp] for k in (4, 5, 6)]
        fit.append(quad_fit(ys))
    ok_k = []
    for k in KS:
        if k < 4:
            continue
        y = T[k].get(nm, (Fr(0), Fr(0)))
        pred = tuple(fit[c][0] + fit[c][1] * k + fit[c][2] * k * k for c in (0, 1))
        ok_k.append(pred == y)
    ALL_OK &= all(ok_k)
    FIT[nm] = fit
    print(f"{nm:5s}: quadratic through k=4,5,6 reproduces EXACTLY k=4..15 : {all(ok_k)}   "
          f"c0={fit[0][0]}+{fit[1][0]}*sqrt2   c1={fit[0][1]}+{fit[1][1]}*sqrt2   c2={fit[0][2]}+{fit[1][2]}*sqrt2")
print("ALL constants exactly quadratic for k=4..15:", ALL_OK)

# the table claimed in identify_w4.md  (c0, c1, c2) each as (rational, rational*sqrt2)
F = Fr
CLAIM = {
    "J_-a": ((F(-23, 8), 0), (F(5, 2), 0), (F(-1, 2), 0)),
    "J_0": ((F(23, 16), F(1, 8)), (F(-5, 4), 0), (F(1, 4), 0)),
    "J_-b": ((F(-37, 16), 0), (F(3, 2), 0), (F(-1, 4), 0)),
    "J_a": ((F(37, 16), 0), (F(-3, 2), 0), (F(1, 4), 0)),
    "K_-a": ((F(5, 24), 0), (F(-2, 3), 0), (F(1, 6), 0)),
    "K_0": ((F(3, 8), F(-1, 8)), (F(-1, 12), 0), (0, 0)),
    "K_-b": ((F(1, 2), 0), (F(-1, 6), 0), (0, 0)),
    "K_a": ((F(-35, 24), 0), (F(1), 0), (F(-1, 6), 0)),
}
ok = True
for nm in names:
    mine = tuple((FIT[nm][0][j], FIT[nm][1][j]) for j in range(3))
    claim = CLAIM.get(nm)
    same = claim is not None and all((Fr(mine[j][0]), Fr(mine[j][1])) == (Fr(claim[j][0]), Fr(claim[j][1])) for j in range(3))
    ok &= same
    print(f"  claimed table row {nm}: matches my exact fit: {same}")
print("claimed table == my exact table:", ok and set(names) == set(CLAIM))

# k = 3 row
print("k=3 actual coefficients (p0 + p1 sqrt2):", {nm: (str(v[0]), str(v[1])) for nm, v in T[3].items()})
for nm in names:
    pred3 = tuple(FIT[nm][c][0] + 3 * FIT[nm][c][1] + 9 * FIT[nm][c][2] for c in (0, 1))
    act = T[3].get(nm, (Fr(0), Fr(0)))
    print(f"   {nm:5s}: poly(3) = {pred3[0]} + {pred3[1]} sqrt2 ; actual k=3 = {act[0]} + {act[1]} sqrt2 ; equal: {pred3 == act}")
