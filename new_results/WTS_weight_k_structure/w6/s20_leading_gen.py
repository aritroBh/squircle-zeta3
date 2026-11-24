"""w6 identification step 13 (bonus, structure of the leading coefficients): fit W_w(k) for w=7 (d=6, k=7..15) and w=8 (d=7, k=8..15)
at dps 150 and test the conjecture  g_m := leading k-coefficient of W_{m+2}(k) = lambda_m * [y^(m+2)] E(y),
E(y) = 2^y Gamma(1-y) Gamma(1+y/2) / Gamma(1-y/2)  (= exp(sum_n s_n y^n/n), s_n = eta(n) odd n, zeta(n) even n), lambda_m rational.
Uses vals_150 (weights 4,5 from s12; weights 7,8 new) and the 400-dps c_4 for weight 6.  All printed numbers from this script."""
import sys, os, json, glob
sys.dont_write_bytecode = True
W = os.path.dirname(os.path.abspath(__file__))
from fractions import Fraction as Fr
import mpmath as mp
DPS = 150
mp.mp.dps = DPS + 10
print("mp.mp.dps =", mp.mp.dps)
def tomp(f): return mp.mpf(f.numerator) / mp.mpf(f.denominator)
def poly_through(pts):
    n = len(pts); co = [Fr(0)] * n
    for a in range(n):
        ka, ya = pts[a]; num = [Fr(1)]; den = Fr(1)
        for c in range(n):
            if c == a: continue
            kc = pts[c][0]; new = [Fr(0)] * (len(num) + 1)
            for t, v in enumerate(num): new[t + 1] += v; new[t] -= kc * v
            num = new; den *= (ka - kc)
        for t, v in enumerate(num): co[t] += ya * v / den
    return co
def load(d, ks):
    return {k: Fr(json.load(open(os.path.join(W, "vals_150", f"k{k}_d{d}.json")))["value"]) for k in ks}
def fit_weight(d, ks):
    Y = load(d, ks); deg = d - 1            # degree w-2 = d-1
    c = poly_through([(Fr(k), Y[k]) for k in ks[:deg + 1]])
    chk = [tomp(Y[k] - sum(c[j] * Fr(k) ** j for j in range(deg + 1))) for k in ks[deg + 1:]]
    return c, chk
G = {}
for d, k0 in ((3, 4), (4, 5), (6, 7), (7, 8)):
    ks = list(range(k0, 16))
    c, chk = fit_weight(d, ks)
    print(f"weight {d+1}: fit degree {d-1} through k={ks[:d]}, check residuals at k={ks[d:]}: {[mp.nstr(x, 3) for x in chk]}")
    G[d - 1] = tomp(c[d - 1])            # g_m, m = d-1-... weight w=d+1: m = w-2 = d-1
# weight 6 (m=4) from the 400-dps coefficient (identified closed form)
c4 = mp.mpf(json.load(open(os.path.join(W, "coef_400.json")))["coefficients"][4])
G[4] = c4
# weight 3 (m=1) and 2 (m=0), exact
lam = mp.log(2); z3 = mp.zeta(3); pi2 = mp.pi ** 2
G[0] = (pi2 + 6 * lam ** 2) / 24
G[1] = -(z3 / 8 + pi2 * lam / 24 + lam ** 3 / 12)
# E(y) Taylor coefficients
E = lambda y: mp.power(2, y) * mp.gamma(1 - y) * mp.gamma(1 + y / 2) / mp.gamma(1 - y / 2)
Ec = mp.taylor(E, 0, 9)
print("E coefficients E_n (n=2..9):", [mp.nstr(Ec[n], 12) for n in range(2, 10)])
res = {}
for m in sorted(G):
    n = m + 2
    ratio = G[m] / Ec[n]
    # nearest rational with small denominator
    fr = Fr(str(mp.nstr(ratio, 40))).limit_denominator(10 ** 6)
    resid = G[m] - tomp(fr) * Ec[n]
    print(f"m={m} (weight {n}): g_m/E_{n} = {mp.nstr(ratio, 30)}   ~ {fr}   residual g_m - lambda*E_n = {mp.nstr(resid, 4)}  (dps {mp.mp.dps})")
    res[m] = str(fr)
json.dump({"lambda": res}, open(os.path.join(W, "lambda_m.json"), "w"))
