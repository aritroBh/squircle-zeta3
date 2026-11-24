# Planner check of DEC5E's FINAL closed form: rebuild P5(1)/L = Wt1+...+Wt5 from the stored fits with own mpmath atom values,
# compare with DEC5D's independent 190-digit gold value; then substitute ZQ5's golden ladder for L(chi_5,5) (fully classical form).
import pickle, re, mpmath as mp
from fractions import Fraction as Fr
mp.mp.dps = 210
W = "<project>/research/08_theory/DEC5E_decagon_p5_closed/"
s5 = mp.sqrt(5); phi = (1 + s5)/2; rho = 1/phi; lf = mp.log(phi); l2 = mp.log(2)
chi5 = lambda a: [0, 1, -1, -1, 1][a % 5]
Lchi5 = mp.mpf(5)**-5 * mp.fsum(chi5(a) * mp.zeta(5, mp.mpf(a)/5) for a in range(1, 5))
base = {"L3m": mp.polylog(3, -rho), "L4m": mp.polylog(4, -rho), "L4q": mp.polylog(4, rho**2), "L5m": mp.polylog(5, -rho),
        "L5q": mp.polylog(5, rho**2), "Lchi5": Lchi5, "l2": l2, "lf": lf, "z2": mp.zeta(2), "z3": mp.zeta(3), "z5": mp.zeta(5)}
tok = re.compile(r"(Lchi5|L3m|L4m|L4q|L5m|L5q|l2|lf|z2|z3|z5)(\^(\d+))?")
def atom(name):
    v, pos = mp.mpf(1), 0
    for m in tok.finditer(name):
        assert m.start() == pos, name; pos = m.end()
        v *= base[m.group(1)] ** int(m.group(3) or 1)
    assert pos == len(name), name
    return v
def coef(pq): p, q = pq; return mp.mpf(Fr(p).numerator)/Fr(p).denominator + mp.mpf(Fr(q).numerator)/Fr(q).denominator * s5
Wt = {1: 1200*l2 + (3000 - 600*s5)*lf,
      2: (-mp.mpf(1425)/16 - mp.mpf(5415)/16*s5)*mp.zeta(2) - 2400*l2**2 + (-12000 + 2400*s5)*l2*lf + (-mp.mpf(86925)/4 + mp.mpf(17325)/2*s5)*lf**2}
for k, f in [(3, "fit_wt3.pkl"), (4, "fit_wt4.pkl"), (5, "fit_wt5_raw.pkl")]:
    rel = pickle.load(open(W + f, "rb"))["rel"]
    Wt[k] = mp.fsum(coef(pq) * atom(a) for a, pq in rel.items())
    if k == 5: c_l = rel["Lchi5"]; c_z5 = rel["z5"]
gold = mp.mpf(pickle.load(open("<project>/research/08_theory/DEC5D_decagon_order5_depth/p5_gold.pkl", "rb"))["P5overL"])
tot = mp.fsum(Wt.values())
print("P5(1)/L from the closed form :", mp.nstr(tot, 50))
print("DEC5D gold (independent)     :", mp.nstr(gold, 50))
print("|closed form - gold| =", mp.nstr(abs(tot - gold), 3), "  (gold has about 190 digits)")
print("L(chi_5,5) coefficient (1, sqrt5 parts):", c_l, "| zeta(5) coefficient:", c_z5)
# ZQ5 ladder: L(chi5,5) = sqrt5*{(1728/6875)[Li5(rho)-Li5(-rho)] - (512/6875)Li5(rho^3) + (16/6875)Li5(rho^6) + (648/34375)lf^5 - (48/6875)pi^2 lf^3 + (64/20625)pi^4 lf}
lad = s5*(mp.mpf(1728)/6875*(mp.polylog(5, rho) - mp.polylog(5, -rho)) - mp.mpf(512)/6875*mp.polylog(5, rho**3) + mp.mpf(16)/6875*mp.polylog(5, rho**6)
          + mp.mpf(648)/34375*lf**5 - mp.mpf(48)/6875*mp.pi**2*lf**3 + mp.mpf(64)/20625*mp.pi**4*lf)
print("ZQ5 ladder - L(chi_5,5) =", mp.nstr(lad - Lchi5, 3), "-> the closed form is classical with Li5 at -rho, rho^2, rho^3, rho^6 (rho = 1/phi)")
