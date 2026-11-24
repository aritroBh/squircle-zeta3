"""w6 identification step 24: final verification pass.  All printed numbers from this script.
 A. c_4 closed form vs the fitted c_4 of the four independent high-precision runs (dps 220, 250, 400, 900), residual printed with mp.mp.dps;
    LLL relation of c_4 against the 12-element basis {8 Gamma-world monomials, Li4(1/2) pi^2, Li4(1/2) lam^2, Li5(1/2) lam, Li6(1/2)}
    at D=850 (dps-900 run) and D=380 (dps-400 run) and D=235 (dps-250 run).
 B. weight-4 next-to-leading coefficient c_1^(4): closed form vs fitted (dps-400 data, 395 digits; dps-150 data, 140 digits).
 C. leading-coefficient law  g_m = eps_m A_{m+1} / (2 m!) [y^(m+2)] E(y)  (A_n Euler zigzag numbers, eps_m = +,-,-,+ periodic),
    checked against the fitted leading coefficients for weights 2..8 (residuals printed)."""
import sys, os, json
sys.dont_write_bytecode = True
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fractions import Fraction as Fr
from w6lib import *
from w6lll import find_rel2
def cf4(lam, pi2, z3, z5):
    return (168 * lam ** 6 + 420 * pi2 * lam ** 4 + 378 * pi2 ** 2 * lam ** 2 + 61 * pi2 ** 3 + 5040 * z3 * lam ** 3
            + 2520 * pi2 * z3 * lam + 3780 * z3 ** 2 + 22680 * z5 * lam) / 362880
print("=== A. c_4")
for run, dps in ((220, 220), (250, 250), (400, 400), (900, 900)):
    mp.mp.dps = dps - 10
    cf = cf4(mp.log(2), mp.pi ** 2, mp.zeta(3), mp.zeta(5))
    c4 = mp.mpf(coef(run, 4))
    print(f"run {run}: closed form - fitted c_4 = {mp.nstr(cf - c4, 5)}   (mp.mp.dps = {mp.mp.dps})")
for run, D in ((900, 850), (400, 380), (250, 235)):
    mp.mp.dps = D + 25
    lam = mp.log(2); z3 = mp.zeta(3); z5 = mp.zeta(5); pi2 = mp.pi ** 2; h = mp.mpf(1) / 2
    Li4 = mp.polylog(4, h); Li5 = mp.polylog(5, h); Li6 = mp.polylog(6, h)
    B = [("lam^6", lam ** 6), ("pi^2 lam^4", pi2 * lam ** 4), ("pi^4 lam^2", pi2 ** 2 * lam ** 2), ("pi^6", pi2 ** 3), ("z3 lam^3", z3 * lam ** 3),
         ("pi^2 z3 lam", pi2 * z3 * lam), ("z3^2", z3 ** 2), ("z5 lam", z5 * lam), ("Li4(1/2) pi^2", Li4 * pi2), ("Li4(1/2) lam^2", Li4 * lam ** 2),
         ("Li5(1/2) lam", Li5 * lam), ("Li6(1/2)", Li6)]
    c4 = mp.mpf(coef(run, 4)); mp.mp.dps = D + 25
    good, internal = find_rel2(c4, B, D)
    r = good[0]
    print(f"LLL c_4 vs 12-element basis, run {run}, D={D}, mp.mp.dps={mp.mp.dps}: max|coeff|={r['maxcoeff']}, residual {mp.nstr(r['res'],4)}, relation {r['terms']}")
print("=== B. c_1^(4)")
def c1_4(S):
    lam, L, pi2, z3 = S["lam"], S["L"], S["pi2"], S["z3"]
    B1, B2, B3 = S["B1"], S["B2"], S["B3"]
    return (11 * pi2 ** 2 - 210 * L * z3 - 100 * L ** 2 * pi2 + 160 * L ** 4 + 360 * lam * z3 + 100 * lam * L * pi2 + 80 * lam * L ** 3
            + 55 * lam ** 2 * pi2 - 120 * lam ** 2 * L ** 2 - 40 * lam ** 3 * L + 40 * lam ** 4
            - 480 * lam * (B1 + B2) - 960 * L * (B1 - B2) + 1920 * L * B3) / 480
def setupS(dps):
    mp.mp.dps = dps
    r2 = mp.sqrt(2); b = r2 - 1
    return dict(lam=mp.log(2), L=mp.log(1 + r2), pi2=mp.pi ** 2, z3=mp.zeta(3), B1=mp.polylog(3, b), B2=mp.polylog(3, -b), B3=mp.polylog(3, 1 / r2))
S = setupS(395)
c = mp.mpf(json.load(open(os.path.join(W, "lower_coef_400.json")))["W4"][1])
print("c_1^(4) closed form - fitted (dps-400 data):", mp.nstr(c1_4(S) - c, 5), f"(mp.mp.dps = {mp.mp.dps})")
S = setupS(140)
c = mp.mpf(json.load(open(os.path.join(W, "lower_coef.json")))["W4"][1])
print("c_1^(4) closed form - fitted (dps-150 data):", mp.nstr(c1_4(S) - c, 5), f"(mp.mp.dps = {mp.mp.dps})")
print("c_1^(4) =", mp.nstr(c1_4(setupS(60)), 50))
print("=== C. leading-coefficient law")
mp.mp.dps = 150
A = [1, 1, 1, 2, 5, 16, 61, 272, 1385]            # Euler zigzag numbers A_0..A_8 (sec x + tan x = sum A_n x^n / n!)
import sympy
x = sympy.symbols("x")
ser = sympy.series(sympy.sec(x) + sympy.tan(x), x, 0, 10).removeO()
assert [int(ser.coeff(x, n) * sympy.factorial(n)) for n in range(9)] == A, "zigzag check"
eps = [1, -1, -1, 1]
E = lambda y: mp.power(2, y) * mp.gamma(1 - y) * mp.gamma(1 + y / 2) / mp.gamma(1 - y / 2)
Ec = mp.taylor(E, 0, 9)
# fitted leading coefficients g_m, weights 2..8 : reuse s20's definitions
lam = mp.log(2); z3 = mp.zeta(3); pi2 = mp.pi ** 2
def fitlead(d, ks):
    def poly_through(pts):
        n = len(pts); co = [Fr(0)] * n
        for a in range(n):
            ka, ya = pts[a]; num = [Fr(1)]; den = Fr(1)
            for cc in range(n):
                if cc == a: continue
                kc = pts[cc][0]; new = [Fr(0)] * (len(num) + 1)
                for t, v in enumerate(num): new[t + 1] += v; new[t] -= kc * v
                num = new; den *= (ka - kc)
            for t, v in enumerate(num): co[t] += ya * v / den
        return co
    Y = {k: Fr(json.load(open(os.path.join(W, "vals_150", f"k{k}_d{d}.json")))["value"]) for k in ks}
    cc = poly_through([(Fr(k), Y[k]) for k in ks[:d]])
    return mp.mpf(cc[d - 1].numerator) / mp.mpf(cc[d - 1].denominator)
G = {0: (pi2 + 6 * lam ** 2) / 24, 1: -(z3 / 8 + pi2 * lam / 24 + lam ** 3 / 12),
     2: fitlead(3, list(range(4, 16))), 3: fitlead(4, list(range(5, 16))),
     4: mp.mpf(coef(900, 4)), 5: fitlead(6, list(range(7, 16))), 6: fitlead(7, list(range(8, 16)))}
for m in range(7):
    pred = mp.mpf(eps[m % 4]) * A[m + 1] / (2 * mp.factorial(m)) * Ec[m + 2]
    print(f"m={m} (weight {m+2}): g_m - eps_m A_{m+1}/(2 m!) E_{m+2} = {mp.nstr(G[m] - pred, 4)}   (|g_m| = {mp.nstr(abs(G[m]), 8)}, dps {mp.mp.dps})")
