"""w6 step 6: assemble values_w6.json and result_w6.md from the stored script outputs (no new numerics)."""
import sys, os, json, re
sys.dont_write_bytecode = True
W = os.path.dirname(os.path.abspath(__file__))
import mpmath as mp
mp.mp.dps = 110

A70 = json.load(open(os.path.join(W, "analysis_vals.json")))
A110 = json.load(open(os.path.join(W, "analysis_vals_110.json")))

def val(dirn, k, d=5):
    return json.load(open(os.path.join(W, dirn, f"k{k}_d{d}.json")))

ks = list(range(5, 16))
v70 = {k: val("vals", k) for k in ks}
v110 = {k: val("vals_110", k) for k in ks}
VDIFF = {k: abs(mp.mpf(v70[k]["value"]) - mp.mpf(v110[k]["value"])) for k in ks}
VAL66 = {k: mp.nstr(mp.mpf(v110[k]["value"]), 100) for k in ks}
VAL80 = {k: mp.nstr(mp.mpf(v110[k]["value"]), 80) for k in ks}
VDMAX = max(VDIFF.values())

# coefficients: 110-run, rounded to the digits that agree with the 70-run
RELIABLE = 100
c110 = [mp.mpf(x) for x in A110["fit"]["coefficients"]]
coef_str = [mp.nstr(x, RELIABLE) for x in c110]
coef_full = [mp.nstr(x, 100) for x in c110]

def maxres(prefix):
    rs = [abs(mp.mpf(c["residual"])) for c in A70["controls"] if c["name"].startswith(prefix)]
    return mp.nstr(max(rs), 4), len(rs), all(c["pass"] for c in A70["controls"] if c["name"].startswith(prefix))

w2max, w2n, w2ok = maxres("W_2(k=")
w3max, w3n, w3ok = maxres("W_3(k=")
by_name = {c["name"]: c for c in A70["controls"]}
ctl_rows = []
ctl_rows.append({"name": f"weight 2: W_2(k) = (pi^2+6 log^2 2)/24, all {w2n} k = 3..15 (max |residual|)", "residual": w2max, "pass": w2ok})
ctl_rows.append({"name": f"weight 3: W_3(k) = W3' - (k-3) S (W3' from thm:M3 in main.tex), all {w3n} k = 3..15 (max |residual|)", "residual": w3max, "pass": w3ok})
for nm in ["W_4(k=3) - W_4 printed in thm:M3 (48 digits)", "W_4+W_3+W_2 at k=3 - M3 printed in thm:M3 (51 digits)",
           "W_6(k=5) - stored M5a quad5_100.json (100 digits)", "W_5(k=5) - stored M5a quad5_100.json (100 digits)",
           "W_4(k=5) - stored M5a quad5_100.json (100 digits)",
           "sum_{w=2..6} W_w(5) - M5 (M5a/M5_400.txt, direct integrand quadrature)"]:
    c = by_name[nm]
    ctl_rows.append({"name": nm, "residual": c["residual"], "pass": c["pass"]})
by110 = {c["name"]: c for c in A110["controls"]}
ctl110 = [{"name": nm + " [110-dps rerun]", "residual": by110[nm]["residual"], "pass": by110[nm]["pass"]}
          for nm in ["W_6(k=5) - stored M5a quad5_100.json (100 digits)",
                     "sum_{w=2..6} W_w(5) - M5 (M5a/M5_400.txt, direct integrand quadrature)"]]
ctl_rows += ctl110
ds70 = A70["degree_search"]; ds110 = A110["degree_search"]
fit70 = A70["fit"]; fit110 = A110["fit"]

out = {
    "weight": 6,
    "dps_used_primary": 70, "dps_used_rerun": 110,
    "k_requested": [5, 6, 7, 8, 9, 10], "k_extension": [11, 12, 13, 14, 15],
    "values_70dps": {str(k): v70[k]["value"] for k in ks},
    "values_100digits_110dps_run": {str(k): VAL66[k] for k in ks},
    "max_abs_diff_values_70run_vs_110run": mp.nstr(VDMAX, 4),
    "values_110dps": {str(k): v110[k]["value"] for k in ks},
    "quadrature_err_est_70dps": {str(k): v70[k]["err_est"] for k in ks},
    "quadrature_err_est_110dps": {str(k): v110[k]["err_est"] for k in ks},
    "pieces_70dps": {str(k): v70[k]["pieces"] for k in ks},
    "controls": ctl_rows,
    "degree_search_70dps": ds70, "degree_search_110dps": ds110,
    "polynomial_degree": fit110["degree"], "regime": fit110["regime"], "regime_k": fit110["regime_k"],
    "k5_status": "OUTSIDE the polynomial regime (top-weight case d=k)",
    "k5_defect_W6_minus_poly": "-0.219366913442151079860842560263653086600349859821953115922874",
    "coefficients_100digits": coef_str,
    
    "coefficients_70dps_run_65digits": fit70["coefficients"],
    "reliable_digits": RELIABLE,
    "fit_residuals_70dps": fit70["residuals"], "fit_residuals_110dps": fit110["residuals"],
    "leading_coefficient_closed_form": "c_4 = (168 lam^6 + 420 pi^2 lam^4 + 378 pi^4 lam^2 + 61 pi^6 + 5040 zeta3 lam^3 + 2520 pi^2 zeta3 lam + 3780 zeta3^2 + 22680 zeta5 lam)/362880, lam = log 2 (PSLQ, verified to 3.1e-105 at dps 110)",
}
json.dump(out, open(os.path.join(W, "values_w6.json"), "w"), indent=1)
print("wrote values_w6.json")

# ------------------------------------------------------------------ result_w6.md
VDIFF = {k: abs(mp.mpf(v70[k]["value"]) - mp.mpf(v110[k]["value"])) for k in ks}
VAL66 = {k: mp.nstr(mp.mpf(v110[k]["value"]), 100) for k in ks}
VAL80 = {k: mp.nstr(mp.mpf(v110[k]["value"]), 80) for k in ks}
VDMAX = max(VDIFF.values())

def row_vals():
    lines = ["| k | W_6(k) (110-dps run, 80 digits shown; 100 in values_w6.json) | mp.quad err_est (70 / 110 dps run) | actual |70-run - 110-run| |", "|---|---|---|---|"]
    for k in ks:
        tag = "" if k <= 10 else " (extension)"
        lines.append(f"| {k}{tag} | {VAL80[k]} | {v70[k]['err_est']} / {v110[k]['err_est']} | {mp.nstr(VDIFF[k], 3)} |")
    return "\n".join(lines)

def row_ctl():
    lines = ["| control | residual | pass |", "|---|---|---|"]
    for c in ctl_rows:
        lines.append(f"| {c['name']} | {c['residual']} | {'PASS' if c['pass'] else 'FAIL'} |")
    return "\n".join(lines)

def row_fd(ds):
    lines = ["| set of k | order n | max abs Delta^n |", "|---|---|---|"]
    for name in ("requested_k5_10", "k6_10", "k6_to_15"):
        if name not in ds: continue
        kk = ds[name]["k"]
        for n, m in ds[name]["max_abs_diff_by_order"].items():
            lines.append(f"| {kk[0]}..{kk[-1]} | {n} | {m} |")
    return "\n".join(lines)

res70 = fit70["residuals"]; res110 = fit110["residuals"]
fitrows = ["| k | W_6(k) - p(k), 70-dps run | W_6(k) - p(k), 110-dps run |", "|---|---|---|"]
for k in ks:
    role = "outside regime" if k == 5 else ("fitted" if k <= 10 else "CHECK")
    fitrows.append(f"| {k} ({role}) | {res70[str(k)]} | {res110[str(k)]} |")

md = f"""# W_6(k): the weight-six part of M_k, k = 5..15

## VERDICT

**The polynomial law holds for k >= 6 with degree 4, and it FAILS at k = 5.**
(Numerics to about 100 digits plus an exact rational check of the kernel residues for k <= 15; this is evidence for k <= 15, not an all-k proof.)

* On the requested window k = 5..10 alone the law is NOT confirmed: the 5th finite difference is
  Delta^5 = {ds70['requested_k5_10']['max_abs_diff_by_order']['5']} (first window: 0.2193669134421510798608426), not zero, so no degree d is found
  there. The cause is k = 5 only.
* k = 5 is the top-weight case: reading the reducer, the weight-six keys (i+j = 5) at k = 5 are just the simple-pole parts of the closed-form integrand,
  whereas for k >= 6 they also receive the integration-by-parts contributions of the degree 6..k keys. It is the analogue of k = 1 for weight two and
  k = 2 for weight three, which the WT2/WT3 notes already record as off their laws (not recomputed here). Its defect is
  W_6(5) - p(5) = -0.219366913442151079860842560263653086600349859821953115922874 (110-dps run), where p is the degree-4
  polynomial fixed by k = 6..10.
* k = 6..10 contains only 5 points, which fix a quartic with NO check left over. The check therefore
  comes from the extension k = 11..15 (computed beyond the requested range, same code path). On k = 6..15 (10 points)
  the 5th differences vanish: max |Delta^5| = {ds70['k6_to_15']['max_abs_diff_by_order']['5']} at 70 dps and
  {ds110['k6_to_15']['max_abs_diff_by_order']['5']} at 110 dps (5 windows each), and 4th differences are nonzero
  (|Delta^4| ~ 8.0565), so the smallest degree is d = 4, as expected from weight 2 (degree 0) and weight 3 (degree 1).
* The quartic fitted through k = 6..10 predicts k = 11..15 with absolute error <= {mp.nstr(max(abs(mp.mpf(res110[str(k)])) for k in range(11,16)), 4)} (110-dps run;
  the values themselves are of size up to 2.7e3, so this is rounding level).

* EXACT kernel-level confirmation (rational arithmetic, `s07_exact_kernel_poly.py`, s07_output.log): every residue of the three degree-5 kernels
  K_50, K_41, K_32 (15 nonzero Q(sqrt2)-components over the letters 0, a, -a, -b; none at b) is exactly a polynomial in k of degree <= 4 on k = 6..15
  (10 points, all 5th and higher exact finite differences are 0), degree exactly 4 for several (e.g. the residues at a and -a), and NONE of the 15 is
  polynomial once k = 5 is included. Since W_6(k) is linear in these residues, W_6(k) is exactly a polynomial of degree <= 4 in k for 6 <= k <= 15
  (the high-precision numerics only supply its coefficients). This is a finite check for k <= 15, not an all-k proof.

## Result

W_6(k) = c_0 + c_1 k + c_2 k^2 + c_3 k^3 + c_4 k^4  for k >= 6, from the 110-dps run, printed to {RELIABLE} digits.
Reliability: in the 110-dps run the fit through k = 6..10 and the fit through k = 11..15 agree to {A110['fit']['worst_relative_difference_first_vs_last_fit']} relative
(worst coefficient), and c_4 agrees with its PSLQ closed form below to 3.1e-105. The requested 70-dps run gives the same coefficients but only to
64 digits (relative difference 70-run vs 110-run: 1.9e-65 worst, s04_output.log), because its quadrature loses about 3 digits and
the fit amplifies the rounding at k = 0.

```
c_0 = {coef_str[0]}
c_1 = {coef_str[1]}
c_2 = {coef_str[2]}
c_3 = {coef_str[3]}
c_4 = {coef_str[4]}
```

(Same strings in values_w6.json under coefficients_100digits.)

**Leading coefficient identified (exploratory PSLQ, s04/s05; relation printed by the script, closed form checked
against the fitted c_4 to 3.1e-105 at dps 110):**

c_4 = ( 168 L^6 + 420 pi^2 L^4 + 378 pi^4 L^2 + 61 pi^6 + 5040 zeta(3) L^3 + 2520 pi^2 zeta(3) L + 3780 zeta(3)^2 + 22680 zeta(5) L ) / 362880,   L = log 2.

This is pure (no sqrt2, no Li_n(1/2) needed; the PSLQ basis was {{L^6, pi^2 L^4, pi^4 L^2, pi^6, zeta3 L^3, pi^2 zeta3 L, zeta3^2, zeta5 L}}).
Same family as the lower weights: weight 2 constant (pi^2+6L^2)/24, weight 3 slope -(zeta3/8 + pi^2 L/24 + L^3/12). No identification of c_3..c_0
or of the k = 5 defect was found (PSLQ over the 8-element basis: no relation; over a 12-element basis only a spurious relation with coefficients up to 7 digits
(5-7 digits each), discarded).

## Values (requested k = 5..10, extension k = 11..15)

{row_vals()}

Note on precision: mp.quad's own error estimates (1e-79..1e-84 at 70 dps) are over-optimistic. Comparing the 70-dps values with the 110-dps rerun,
the true absolute error of the 70-dps values is between 2.8e-67 (k = 5) and {mp.nstr(VDMAX, 3)} (k = 15), about 66 significant digits (last column).
The 110-dps values are good to about 1e-104 absolute: W_6(5) and the sum W_2..W_6 at k = 5 agree with the stored M5a data to 4.3e-101 (limited by the 100 stored digits) and 7.7e-105
(400 stored digits), the 5th differences on k = 6..15 are <= 4.0e-106, and the two fits agree to 5e-105 relative.

Pipeline per k: `s01_reduce.py` (exact fast reducer; the copied pickle binteg_k.pkl was verified equal to the closed-form rebuild
`build_integrand_pf(k)` in exact arithmetic for every k) -> finals_k.json (exact rational data; at degree 5 only real simple poles occur,
poles at b never, so every key converges) -> `s02_eval.py` (tanh-sinh, maxdegree 12, split at v = 1, integration in the exact distance to the
endpoint so there is no cancellation near v = b, a) -> `s03_analyze.py`.

## Controls (same code path)

{row_ctl()}

(`mp.mp.dps` printed as 70 for the main run, 110 for the rerun; set after all imports. Tolerances: 1e-58 at 70 dps, 1e-98 at 110 dps
except where the comparison datum has fewer digits.)

## Finite differences of W_6 (70-dps run)

{row_fd(ds70)}

## Fit residuals W_6(k) - p(k), p through k = 6..10 (exact rational interpolation of the decimal strings)

{chr(10).join(fitrows)}

Residuals at k = 11..15 are rounding level of the stored 70 (resp. 110) significant digits times the size of W_6(k); k = 5 is the defect above.

## Cross-check at k = 5

M5 is stored in research/08_theory/M5a_next_coefficient (M5_400.txt, direct integrand quadrature, independent of the reducer). The sum
W_2(5)+W_3(5)+W_4(5)+W_5(5)+W_6(5) from the fast reducer equals it to 2.6e-67 (70 dps) and 7.7e-105 (110 dps); W_6(5) also
equals the stored slow-reducer value (quad5_100.json) to 2.8e-67 (70 dps).

## Bonus: lower-weight degrees from the same data (70-dps run)

Same finite-difference test, restricted to k >= weight (the top-weight case k = weight - 1 excluded), on the data computed here:
weight 5: degree {A70['bonus_lower_weights']['5']['degree']} (k = {A70['bonus_lower_weights']['5']['regime_k'][0]}..{A70['bonus_lower_weights']['5']['regime_k'][-1]});
weight 4: degree {A70['bonus_lower_weights']['4']['degree']} (k = {A70['bonus_lower_weights']['4']['regime_k'][0]}..{A70['bonus_lower_weights']['4']['regime_k'][-1]});
weight 3: degree {A70['bonus_lower_weights']['3']['degree']} (k = {A70['bonus_lower_weights']['3']['regime_k'][0]}..{A70['bonus_lower_weights']['3']['regime_k'][-1]}).
Consistent with degree = weight - 2. These windows are short (weight 5: 6 points, 2 windows of Delta^4; weight 4: 7 points; see analysis_vals70.log), so this is weaker evidence than for weight 6.

## Files

s01_reduce.py, s02_eval.py, s03_analyze.py, s04_compare_pslq.py, s05_verify_c4_defect.py, s06_write_outputs.py, s07_exact_kernel_poly.py, run_reduce_all.sh,
run_eval_all.sh, run_eval_110.sh; finals_k.json, vals/ (70 dps), vals_110/, analysis_vals.json, analysis_vals_110.json, logs
(reduce_log.txt, eval_log.txt, eval_log_110.txt, analysis_vals70.log, analysis_vals110.log, s04_output.log, s05_output.log, s07_output.log); exact_kernel_poly.json; values_w6.json.
Machine: load average ~200-300; all runs nice -n 10, one process at a time; every run took seconds.
"""
open(os.path.join(W, "result_w6.md"), "w").write(md)
print("wrote result_w6.md")
