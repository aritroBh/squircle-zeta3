"""WT4 step G: assemble closed_forms.json; exact rational check of the remark c2 = W2^2 - lam*S - (7/16)zeta(4);
W_4(k) table k=4..14 from the proved formula; k=3 comparison."""
import os, json
HERE = os.path.dirname(os.path.abspath(__file__))
from fractions import Fraction as Fr
import mpmath as mp
mp.mp.dps = 950

lll = json.load(open(os.path.join(HERE, "lll_results.json")))
D1cf = {nm: Fr(v) for nm, v in lll["D1_closed_form"].items()}
D2cf = {nm: Fr(v) for nm, v in lll["D2_closed_form"].items()}

# ---- exact remark check: c2 = D2/2 == W2^2 - lam*S - (7/16)*zeta(4) ----
# basis monomials: pi^4, pi^2*lam^2, lam^4, z3*lam.  zeta(4) = pi^4/90.
W2 = {"pi^4": Fr(1, 576), "pi^2*lam^2": Fr(12, 576), "lam^4": Fr(36, 576)}  # ((pi^2+6 lam^2)/24)^2
lamS = {"z3*lam": Fr(1, 8), "pi^2*lam^2": Fr(1, 24), "lam^4": Fr(1, 12)}
z4 = {"pi^4": Fr(7, 16*90)}
c2 = {}
for src, sgn in ((W2, 1), (lamS, -1), (z4, -1)):
    for m, v in src.items():
        c2[m] = c2.get(m, Fr(0)) + sgn*v
# D2/2 in the same monomial names
alias = {"pi^4": "pi^4", "pi^2*log2^2": "pi^2*lam^2", "log2^4": "lam^4", "z3*log2": "z3*lam"}
half = {alias[nm]: v/2 for nm, v in D2cf.items()}
print("remark c2 - D2/2 =", {m: c2.get(m, Fr(0)) - half.get(m, Fr(0)) for m in set(c2)|set(half)})
remark_ok = all(c2.get(m, Fr(0)) == half.get(m, Fr(0)) for m in set(c2)|set(half))
print("remark c2 = W2^2 - lam*S - (7/16)zeta(4) exact:", remark_ok)
assert remark_ok

# ---- values ----
D1V = mp.mpf(open(os.path.join(HERE, "d1_860.txt")).read().strip())
D2V = mp.mpf(open(os.path.join(HERE, "d2_860.txt")).read().strip())
W4pp = mp.mpf(open("<project>/research/08_theory/M4a_decomposition/W4pp_1060.txt").read().strip())
from math import comb
tab = {}
for k in range(4, 15):
    tab[k] = mp.nstr(W4pp + (k-4)*D1V + comb(k-4, 2)*D2V, 50)
san = json.load(open(os.path.join(HERE, "sanity_k3_table.json")))
tab["3_actual"] = san["W43_thm_950"][:52]
tab["3_quadratic_prediction"] = san["P3"][:52]
tab["3_offset_delta"] = san["delta_P3_minus_W43"][:52]

cf_out = {
 "proved": True,
 "formula": "W_4(k) = W4pp + (k-4)*D1 + binom(k-4,2)*D2  (exact for 4 <= k <= 14; kernel polynomial law EXACT, KSPAN + wt4_A_kernels.py)",
 "W4pp": "PROVED in M4b2_w4pp_nice (see nice_display.md); value W4pp_1060.txt",
 "D2_closed_form": lll["D2_closed_form"],
 "D2_pretty": "D2 = -(1/24)*log2^4 - (1/24)*pi^2*log2^2 - (1/160)*pi^4 - (1/4)*zeta(3)*log2",
 "D2_level": "level 2 (alternating MZV weight 4); Li4(1/2) coefficient is 0",
 "D1_closed_form": lll["D1_closed_form"],
 "D1_pretty": "D1 = -(5/48)lam^4 - (1/12)lam^3 L - (1/4)lam^2 L^2 + (1/6)lam L^3 + (1/3)L^4 - (7/96)pi^2 lam^2 + (5/24)pi^2 lam L - (5/24)pi^2 L^2 - (1/192)pi^4 - (3/8)z3 lam - (7/16)z3 L - lam Li3(b) - 2L Li3(b) - lam Li3(-b) + 2L Li3(-b) + 4L Li3(1/sqrt2)",
 "D1_level": "level 8, but NO sqrt2 and no Li4/Li2 atoms",
 "remark": "c2 = D2/2 = W2^2 - lam*S - (7/16)*zeta(4) with W2 = (pi^2+6 lam^2)/24, S = z3/8 + pi^2 lam/24 + lam^3/12 (exact rational identity, checked)",
 "certificates": ["D1_certificate.json", "D2_certificate.json"],
 "verifier": "verify_wt4_cert.py <D1|D2>",
 "W4_table_50d": tab,
}
json.dump(cf_out, open(os.path.join(HERE, "closed_forms.json"), "w"), indent=1)
print("closed_forms.json written")
for k in list(range(4, 15)) + ["3_actual", "3_quadratic_prediction", "3_offset_delta"]:
    print(f"  W_4({k}) = {tab[k]}")
with open(os.path.join(HERE, "progress.log"), "a") as f:
    f.write("step_G: closed_forms.json written; remark identity exact; W_4(k) table k=4..14 assembled\n")
