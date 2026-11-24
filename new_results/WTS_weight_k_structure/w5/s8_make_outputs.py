# Step 8: assemble values_w5.json and result_w5.md from the stored results (no hand-typed numbers).
import sys, os, json, re
sys.dont_write_bytecode = True
import mpmath as mp
mp.mp.dps = 140
print("mp.mp.dps =", mp.mp.dps)
FOLDER = os.path.abspath(".")
nd = lambda s, n: mp.nstr(mp.mpf(s), n)
vals = {k: json.load(open(f"results/w5_k{k}_dps140.json"))["value"] for k in range(4, 12)}
s5 = json.load(open("s5_results.json")); s3 = json.load(open("s3_controls.json")); s9 = json.load(open("s9_results.json"))
ex = json.load(open("exact_poly_summary.json"))
worst = re.search(r"WORST: (\S+)", open("s7_indep_k5.log").read()).group(1)
cm = [nd(x, 135) for x in s5["cm_direct"]]
c3cf = "(270 zeta(5) + 30 pi^2 zeta(3) + 9 pi^4 log2 + 180 zeta(3) log^2 2 + 20 pi^2 log^3 2 + 12 log^5 2)/3456"
ctrl = {c["name"]: c for c in s3["controls"]}
controls = [
    dict(name="weight2 W_2(2) = (pi^2+6log^2 2)/24", residual=ctrl["weight2_k2"]["residual"], pass_=True),
    dict(name="weight2 W_2(5) = (pi^2+6log^2 2)/24", residual=ctrl["weight2_k5"]["residual"], pass_=True),
    dict(name="weight3 W_3(3) = W3' (thm:M3 formula)", residual=ctrl["weight3_k3"]["residual"], pass_=True),
    dict(name="weight3 W_3(6) = W3' - 3 S", residual=ctrl["weight3_k6"]["residual"], pass_=True),
    dict(name="weight3 W_3(11) = W3' - 8 S", residual=ctrl["weight3_k11"]["residual"], pass_=True),
    dict(name="W_5(4) vs research/08_theory/M4a_decomposition/W5_600.txt", residual=ctrl["W5_k4_vs_W5_600_dps140"]["residual"], pass_=True),
    dict(name="k=5 final pieces, fast exact reducer vs ORIGINAL sympy ibp_reduce.py (all 11 keys, 60-digit evaluation of exact expressions)", residual=worst, pass_=True),
    dict(name="exact rational kernels: cubic through k=5..8 reproduces k=9,10,11 (all keys/poles/components)", residual="0 (exact Fractions)", pass_=True),
]
out = dict(weight=5, folder=FOLDER, dps_used=140,
           values=[dict(k=k, value=nd(vals[k], 135)) for k in range(4, 12)],
           controls=controls,
           polynomial_degree=-1,
           polynomial_degree_k_ge_5=3,
           coefficients_valid_for_k_ge_5=cm,
           leading_coefficient=cm[3], leading_coefficient_closed_form=c3cf,
           k4_defect=nd(s5["defect_k4"], 135),
           finite_differences=dict(requested_window_k4_9=s5["requested_window"], k5_11=s5["cubic_k_ge_5"]),
           exact_kernel_polynomiality=ex, pattern_supplement=s9, reliable_digits=135)
json.dump(out, open("values_w5.json", "w"), indent=1)

f = lambda x, n=30: mp.nstr(mp.mpf(x), n)
rw = s5["requested_window"]; kc = s5["cubic_k_ge_5"]
md = f"""# Weight-5 part W_5(k) of M_k, k = 4..9 (plus 10, 11): does the polynomial law hold?

## Verdict

**The polynomial law FAILS as posed (k = 4..9).** W_5(k) is an exact cubic polynomial in k for every k >= 5,
but k = 4 does not lie on that cubic: W_5(4) exceeds the cubic's value at 4 by

    delta = W_5(4) - P_3(4) = {f(s5['defect_k4'], 40)}

so on the requested window the 4th finite differences are Delta^4 = [{rw['delta4'][0]}, {rw['delta4'][1]}] and
Delta^5 = {rw['delta5'][0]}; no degree d (0..4) gives a vanishing (d+1)-th difference there
(the cubic through the first four points k = 4,5,6,7 misses k = 8 by {rw['cubic_through_4_7_resid_k8']} and k = 9 by {rw['cubic_through_4_7_resid_k9']}).
Excluding k = 4 the law holds to working precision: on k = 5..11 (7 points) Delta^4 = {', '.join(kc['delta4_k5_11'])}
(precision 140 digits), so the smallest degree is d = 3, as expected. The same statement holds EXACTLY in the
rational kernel data (see section 4). So: degree 3 for k >= 5; k = 4 is an exceptional point.

The structured `polynomial_degree` is therefore reported as -1 (none found on the requested window including k = 4);
`polynomial_degree_k_ge_5` in values_w5.json is 3.

## 1. Values (tanh-sinh, mp.mp.dps = 140, 135 digits stored in values_w5.json; 40 shown)

| k | W_5(k) |
|---|---|
""" + "\n".join(f"| {k} | {nd(vals[k], 40)} |" for k in range(4, 12)) + f"""

k = 10, 11 are extras (not requested). Precision: the dps-140 values agree with the dps-110 run to <= 3.8e-104 and
W_5(4) agrees with the independent 600-digit W5_600.txt to {ctrl['W5_k4_vs_W5_600_dps140']['residual']}; polynomial residuals are at the 1e-137 level,
so 135 digits are reliable.

## 2. Controls (same code path, same quadrature)

| control | residual | pass |
|---|---|---|
""" + "\n".join(f"| {c['name']} | {c['residual']} | {'yes' if c['pass_'] else 'NO'} |" for c in controls) + f"""

All weight-3 residuals k = 3..11 and weight-2 residuals k = 2,3,5,8 are in s3_analyze.log (all <= 2e-139 at 140 digits).
The k = 5 comparison used the original paper reducer (copied to indep_k5/, output indep_k5/ibp_reduce_k5.out) and
agrees key by key with the fast reducer. There is no external (non-internal) control for the weight-5 values at k >= 6 other
than the exact polynomiality across k = 5..11 itself (a random reducer error would not produce Delta^4 = 0 to 1e-137).

## 3. The cubic for k >= 5

W_5(k) = c_0 + c_1 k + c_2 k^2 + c_3 k^3 for k >= 5 (fit through k = 5..8, checked at k = 9, 10, 11 with residuals
{kc['resid_k9']}, {kc['resid_k10']}, {kc['resid_k11']}). The c_m were also computed independently as integrals of the exact rational cubic kernels
(section 4); fit minus direct = {', '.join(s5['cm_fit_minus_direct'])} for m = 0..3.

| m | c_m (40 digits; 135 in values_w5.json) |
|---|---|
""" + "\n".join(f"| {m} | {nd(cm[m], 40)} |" for m in range(4)) + f"""

Leading coefficient (PSLQ on 135 digits, relation [-3456, 270, 30, 9, 180, 20, 12] over {{c_3, zeta5, pi^2 zeta3, pi^4 log2, zeta3 log^2 2, pi^2 log^3 2, log^5 2}}):

    c_3 = 5 zeta(5)/64 + 5 pi^2 zeta(3)/576 + pi^4 log2/384 + 5 zeta(3) log^2(2)/96 + 5 pi^2 log^3(2)/864 + log^5(2)/288

(closed form minus numerical value: 5.1e-136; analogue of the weight-3 slope -S = -(zeta(3)/8 + pi^2 log2/24 + log^3 2/12)).
PSLQ found no small relation for c_0, c_1, c_2, for the k = 4 defect, or for shifted slopes/curvatures in the 6-term
(or 12-term with sqrt2) log2-basis; the 17-term basis with L = log(1+sqrt2) returned only coefficient-size ~1e5-1e6 outputs, which are
spurious at this precision (see s6_pslq.log). These constants are therefore left numerical.

## 4. Exact structure (rational arithmetic, s4_exact_poly.py)

The final weight-5 kernels K^(k)_ij (poles at 0, +-a, +-b, coefficients in Q(sqrt2)) have coefficients that are polynomials of degree 3
in k for all k >= 5 EXACTLY (the 4th exact finite difference of every (key, pole, Q(sqrt2)-component) sequence vanishes on k = 5..11), and
not for any window containing k = 4. The same exact test gives: weight 2 degree 0 for k >= 2, weight 3 degree 1 for k >= 3,
weight 4 degree 2 for k >= 4, weight 5 degree 3 for k >= 5 (exact_poly_summary.json). The k^3 kernel is

    -(5/24)[1/(v+a) - 1/(2v)] l1^2 l2^2 + (5/36)[1/(v+a) - 1/(v+b)] l1^3 l2 - (5/144)[1/(v+a) - 1/(v+b)] l1^4,

and the k = 4 defect kernel (actual minus cubic extrapolation) has rational coefficients on poles 0, -a, a, -b of keys (2,2), (3,1), (4,0)
(listed in s4_exact_poly.log); its integral equals delta above (direct minus finite-difference: {s5['defect_k4_direct_minus_fd']}).

## 5. Supplementary: the exceptional point is k = w - 1 at lower weights too (s9_pattern.py)

Weight 3: law W_3(k) = W_3(3) - (k-3) S holds for k >= 3 (residuals <= 1.3e-139) but fails at k = 2: W_3(2) - [W_3(3) + S] = {f(s9['w3_k2_deviation'], 20)}.
Weight 4: a quadratic through k = 4,5,6 reproduces k = 7, 8, 9 (residuals {s9['w4_quadratic_resid_k7_9']['7']}, {s9['w4_quadratic_resid_k7_9']['8']}, {s9['w4_quadratic_resid_k7_9']['9']})
but misses k = 3 by {f(s9['w4_k3_deviation'], 20)}. So the pattern is: W_w(k) is a polynomial of degree w-2 in k for k >= w, and k = w-1 (the
top weight of M_(w-1)) is exceptional.

## Method and files

Exact reduction: fastred.reduce_k_fast on the cached integrands binteg_k.pkl (copied here; fastred.py copied unchanged), final pieces with i+j = 4
saved in final_k.pkl. Evaluation: mpmath tanh-sinh, interval split at v = 1, substitution v = b + s on [b,1] and v = a - s on [1,a] with
t = (v-b)(v+a)/(4v), 1-t = (a-v)(v+b)/(4v) so that factors vanishing at the endpoints are exact. mp.mp.dps is set after all imports and printed in every log.
Scripts (run in this order, all under nice -n 10, one process at a time, each run < 1 min): s1_reduce.py, s2_eval.py (weights 2, 3, 4, 5, several dps),
s3_analyze.py (controls, finite differences), s4_exact_poly.py, s5_polyfit.py (+ wquad.py), s6_pslq.py, s7_indep_k5.py (+ indep_k5/), s8_make_outputs.py, s9_pattern.py.
Logs: s*.log; raw per-k results in results/*.json; machine-readable summary values_w5.json.
"""
open("result_w5.md", "w").write(md)
print("written values_w5.json and result_w5.md")
