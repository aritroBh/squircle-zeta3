"""EXACT k-dependence of the coefficients of the 10 fixed constants I_{key,pole} = int_b^a dv/(v-r) log^i t log^j(1-t), (i,j) in {(3,0),(2,1)},
r in {0,a,-a,b,-b}.  W_4(k) = sum_n coef_n(k) I_n with coef_n(k) in Q(sqrt2).  Fit quadratic in k through k=4,5,6 (exact Fractions),
verify at k=7..10, and report k=3 separately.  No floating point."""
import sys, os, pickle
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from fractions import Fraction as Fr
import fastred
def tab(k):
    red = pickle.load(open(os.path.join(HERE, f"red_{k}.pkl"), "rb"))
    out = {}
    for (i, j), P in red.items():
        if i + j != 3: continue
        real, ip, im, const = P
        assert not ip and not im and const.is_zero()
        for r, c in real.items():
            assert r.p2 == 0 and r.p3 == 0 and c.p2 == 0 and c.p3 == 0
            out[((i, j), (r.p0, r.p1))] = (c.p0, c.p1)    # coef = p0 + p1 sqrt2
    return out
T = {k: tab(k) for k in range(3, 11)}
keys = sorted(set().union(*[set(T[k]) for k in T]))
print("n constants:", len(keys))
def lag(ks, ys, x):   # Lagrange through points, exact
    tot = Fr(0)
    for a, ya in zip(ks, ys):
        t = Fr(ya)
        for b in ks:
            if b != a: t *= Fr(x - b, a - b)
        tot += t
    return tot
def coeffs(ks, ys):   # monomial coeffs of quadratic through 3 points
    (k1, k2, k3), (y1, y2, y3) = ks, ys
    c2 = ((y3 - y2) / (k3 - k2) - (y2 - y1) / (k2 - k1)) / (k3 - k1)
    c1 = (y2 - y1) / (k2 - k1) - c2 * (k1 + k2)
    c0 = y1 - c1 * k1 - c2 * k1 * k1
    return c0, c1, c2
import json
out = {}
for key in keys:
    row = {}
    ok = True
    for comp in (0, 1):
        ys = [T[k].get(key, (Fr(0), Fr(0)))[comp] for k in (4, 5, 6)]
        c0, c1, c2 = coeffs((4, 5, 6), ys)
        for k in (7, 8, 9, 10):
            y = T[k].get(key, (Fr(0), Fr(0)))[comp]
            if c0 + c1 * k + c2 * k * k != y: ok = False
        row[comp] = (c0, c1, c2)
    y3 = T[3].get(key, (Fr(0), Fr(0)))
    p3 = [row[c][0] + 3 * row[c][1] + 9 * row[c][2] for c in (0, 1)]
    print(f"{key[0]} pole {key[1]}: quad_ok(k=4..10)={ok}")
    print(f"   c0 = {row[0][0]} + {row[1][0]} sqrt2")
    print(f"   c1 = {row[0][1]} + {row[1][1]} sqrt2")
    print(f"   c2 = {row[0][2]} + {row[1][2]} sqrt2")
    print(f"   k=3 actual {y3[0]} + {y3[1]} sqrt2 ; poly(3) {p3[0]} + {p3[1]} sqrt2")
    out[str(key)] = {"c0": [str(row[0][0]), str(row[1][0])], "c1": [str(row[0][1]), str(row[1][1])],
                     "c2": [str(row[0][2]), str(row[1][2])], "k3": [str(y3[0]), str(y3[1])], "quad_ok": ok}
json.dump(out, open(os.path.join(HERE, "results", "id_exact_coeffs.json"), "w"), indent=1)
