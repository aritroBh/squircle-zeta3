"""w6 identification step 25: the three combinations X_key = N[key|0] - N[key|a] - N[key|-b] (= -int_b^a log^i t log^j(1-t) dlog(1-t),
because dlog(1-t) = w_a + w_{-b} - w_0 with w_r = dv/(v-r)) are pure MZV constants: solve the three LLL relations of s18 exactly,
and compare with a direct quadrature in the variable t (int_0^1 log^i(1-u) log^j(u)/u du).  All printed numbers from this script."""
import sys, os, json
sys.dont_write_bytecode = True
W = os.path.dirname(os.path.abspath(__file__))
import sympy as sp
import mpmath as mp
mp.mp.dps = 400
N = {k: mp.mpf(v) for k, v in json.load(open(os.path.join(W, "Nints_400.json")))["N"].items()}
X = {key: N[f"{key}|0"] - N[f"{key}|a"] - N[f"{key}|-b"] for key in ("5,0", "4,1", "3,2")}
x50, x41, x32, z, p = sp.symbols("x50 x41 x32 z p")      # z = zeta(3)^2, p = pi^6
sol = sp.solve([8 * x50 + 3 * x41 - 4 * x32 + 12 * z + p, -16 * x50 + 17 * x41 - 16 * x32 - 12 * z - 2 * p,
                23 * x50 + 8 * x41 - 4 * x32 - 48 * z + 3 * p], [x50, x41, x32])
print("exact solution of the three relations:", sol)
z3sq = mp.zeta(3) ** 2; pi6 = mp.pi ** 6
for key, sym in (("5,0", x50), ("4,1", x41), ("3,2", x32)):
    val = sol[sym].subs({z: 1, p: 0}) * z3sq + sol[sym].subs({z: 0, p: 1}) * pi6
    pred = mp.mpf(sp.Rational(sol[sym].coeff(z)).p) / mp.mpf(sp.Rational(sol[sym].coeff(z)).q) * z3sq + mp.mpf(sp.Rational(sol[sym].coeff(p)).p) / mp.mpf(sp.Rational(sol[sym].coeff(p)).q) * pi6
    print(f"X[{key}] - (rational*zeta3^2 + rational*pi^6) = {mp.nstr(X[key] - pred, 4)}   (dps {mp.mp.dps})")
# independent check: X_ij = int_0^1 log^i(1-u) log^j(u) / u du
mp.mp.dps = 60
for (i, j) in ((5, 0), (4, 1), (3, 2)):
    q = mp.quad(lambda u: mp.log1p(-u) ** i * mp.log(u) ** j / u, [0, 0.5, 1])
    print(f"int_0^1 log^{i}(1-u) log^{j}(u)/u du = {mp.nstr(q, 30)}   vs X[{i},{j}] from N: {mp.nstr(X[f'{i},{j}'], 30)}   diff {mp.nstr(q - X[f'{i},{j}'], 3)} (dps {mp.mp.dps})")
print("X50 vs -120 zeta(6):", mp.nstr(X["5,0"] + 120 * mp.zeta(6), 4))
