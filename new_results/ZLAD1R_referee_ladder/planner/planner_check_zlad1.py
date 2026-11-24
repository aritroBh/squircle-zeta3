# Planner check of ZLAD1/ZLAD1R: (1) the Gamma identity behind Theorem ZLAD1' (numerically, own code);
# (2) the single-zeta rung formula and two closed forms on the exact table (all n <= 150, k <= 8).
import json, mpmath as mp
from fractions import Fraction as Fr
from math import comb, factorial
d = json.load(open("<project>/research/08_theory/ZLAD1_zeta_ladder/mu_exact_N150_K8.json"))
print("monomial keys at k=6, n=3:", sorted(d["6"]["3"].keys()))
def c(k, n, mono):
    if k < 0: return Fr(0)
    return Fr(d[str(k)][str(n)].get(mono, "0"))
def z2pow(e): return "1" if e == 0 else ("z2" if e == 1 else f"z2^{e}")
# rung formula: [zeta(w)] mu_{k,n} = (-1)^w sum_{e=1}^{w//2} (e!/(w-e)) C(w-e,e) [zeta(2)^e] mu_{k-w+2e,n}
bad = 0; tested = 0
for w in range(3, 9):
    for k in range(w, 9):
        for n in range(1, 151):
            lhs = c(k, n, f"z{w}")
            rhs = (-1)**w * sum(Fr(factorial(e), w - e) * comb(w - e, e) * c(k - w + 2*e, n, z2pow(e)) for e in range(1, w // 2 + 1))
            tested += 1; bad += (lhs != rhs)
print(f"single-zeta rungs w=3..8, k<=8, n<=150: {tested} tested, {bad} failures")
cf5 = all(c(5, n, "z5") == Fr(-3, (2*n-1)*(2*n+1)*(2*n+3)) for n in range(1, 151))
cf6 = all(c(6, n, "z2^3") == Fr(n*(5*n-1), 3*(2*n-1)*(2*n+1)*(2*n+3)*(2*n+5)) for n in range(1, 151))
print("closed forms [z5]mu_5 and [z2^3]mu_6 for n<=150:", cf5, cf6)
mp.mp.dps = 50
for (xi, eta) in [(mp.mpf('0.13'), mp.mpf('0.21')), (mp.mpf('-0.3'), mp.mpf('0.17')), (mp.mpc('0.1', '0.2'), mp.mpc('0.25', '-0.05'))]:
    u = xi + eta; s = xi * eta
    # beta_w = u^w - p_w(u, s) with p_w = xi^w + eta^w (Lucas), i.e. beta_w = (xi+eta)^w - xi^w - eta^w
    lhs = mp.nsum(lambda w: (-1)**(int(w)+1) * mp.zeta(int(w)) * ((xi+eta)**int(w) - xi**int(w) - eta**int(w)) / int(w), [2, mp.inf])
    rhs = mp.log(mp.gamma(1+xi) * mp.gamma(1+eta) / mp.gamma(1+u))
    print("Gamma identity at", mp.nstr(xi, 4), mp.nstr(eta, 4), ":", mp.nstr(abs(lhs - rhs), 3))
