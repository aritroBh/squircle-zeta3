# Referee report: remainder O(log n / n^3) in thm:P2asym

Object: `new_results/P2R_remainder_bound/proof.md` (with check.py, check_output.txt).
Paper: `paper/ieee/main.tex`, thm:P2asym (l. 1399) and its sketch. Own numerics: `referee_check.py`, `referee_check_output.txt`.

## Verdict: PROVED (two small repairs to make it self-contained; neither changes the result)

E_n := L_theta(n) - [1+cos th + 2z(mu/n)^eps + c_2/n^2] satisfies, for n >= 9 and every th in (0,pi),
|E_n| <= C z log n / n^3 with C absolute (z = sin^2(th/2)), and E_n = -2 z a_2 log n / n^3 + O(z/n^3).
a_2 > 0 for every th in (0,pi), so the log is genuine at every angle. No angle makes it vanish.

## Line-by-line checks (all passed)

- Lemma A. |(-eps)_k| <= a (k-1)! e^{a H_{k-1}}, |(1+eps)_k| <= k! e^{a H_k}, so |t_k| <= a e^{2a} k^{2a}/(k(k+1)) <= (e^{1/2}/4) k^{-3/2} for |eps| <= 1/4: correct. Uniform convergence on {|eps|<=1/4} x [0,1] gives holomorphy in eps and joint continuity, including z=1. Cauchy circles of radius 1/8 around [0,1/8] stay inside |eps|<=1/4: correct.
- Positivity: k_c(s) = sin^2 th/(s+c+sqrt(1+2cs+s^2)) > 0 since sqrt(...) > |s+c| when c^2<1; z=1 via Gauss, g = 1/(1+eps); z=0 trivially. Compactness gives g_min. Correct.
- Lemma B. Lagrange remainder on the real C^3 function log g; |x| <= (A + K_3/48) eps^2; exp remainder. Correct.
- a_1 in (-1,0): uses -log(1-z) < z/(1-z), so u in (-1,0). Correct. Hence q = (mu/n)^eps <= 1 and lambda_n in (log n, log n + 1).
- Decomposition T_1..T_4 is exact (c_2 = 2z a_2 - 2z(1-z); eps/(1-eps) = 1/(n-2); Lambda_1 = sin^2 th/2 = 2z(1-z)). Each numerical constant checked: (2n-1)n/(n-1)^2 <= 2.39 <= 3; eps^3 <= 1.42/n^3; |T_3| = 4z(1-z)/(n^2(n-2)) <= 6z/n^3; |T_4| <= 2.53 z/n^3.
- (F2) |Lambda_m| <= 2z holds because Lambda_m = int_c^1 P_m and |P_m| <= 1 (the stated reason is right, though the bound is 1-c, not the 2/(2m+1) one gets from the polynomial form).
- Sharpness: (1-q)eps^2 = eps^3 lambda_n + O(eps^4 lambda_n^2), and eps^3 lambda_n - log n/n^3 = O(1/n^3). Only T_1 carries a log. Correct.
- n^{-3} constant d = 2z[a_2(2+a_1)+a_3] - z(3-z-2z^2): re-derived independently (4z a_2 + 2z a_1 a_2 + 2z a_3 - 4z(1-z) + Lambda_2/2, Lambda_2 = 2z(1-z)(1-2z)). Correct, given a_3.
- Uniformity: on all of (0,pi), with the factor z; continuous extension to th = pi. Correct.

## Repairs (minor)

R1. a_2's closed form is cited from the notes, yet a_2 > 0 (the sharpness) rests on it. Add the two-line derivation:
"From t_k = -eps/(k(k+1)) (1 + eps/k) + O(eps^3) we get F = 1 - (1+u) eps + (1+u-Li_2(z)) eps^2 + O(eps^3), because sum z^k/(k^2(k+1)) = Li_2(z) - u - 1. With log(pi eps/sin pi eps) = zeta(2) eps^2 + O(eps^4), a_2 = zeta(2) + (1+u-Li_2(z)) - (1+u)^2/2 = zeta(2) + 1/2 - Li_2(z) - u^2/2."
(I checked this; it also matches a numerical Taylor coefficient to 1e-81 at th = 1 and 2.5.)

R2. Status table row "n^{-3} coefficient d: PROVED given a_3" should say that a_3's closed form is from the notes; d is not needed for the theorem.

## Independent numerics (own code, 80 digits)

Exact length from thm:P2 (hyp2f1 + Legendre series), checked against direct quadrature of the arc length at n = 50, 400 (diff <= 4e-81). a_1, a_2 taken from numerical Taylor coefficients of log g, not from closed forms.

| th | -2 z a_2 | E n^3/log n at n=1e6 | at 1e12 | slope of E n^3 vs log n, 1e11 -> 1e12 |
|---|---|---|---|---|
| 1.0 | -0.6975586640 | -0.64739 | -0.67247 | -0.6975586647 |
| 2.5 | -1.4612325598 | -1.38546 | -1.42335 | -1.4612325613 |

E n^3/log n creeps to the limit only at rate 1/log n (the constant d). The slope of E n^3 in log n matches -2 z a_2 to about 1e-9 at n = 1e12. The log term is real and nonzero at both angles.

## The paper's wording

- Statement of thm:P2asym: correct as written. It could add "uniformly in th, the error is at most C sin^2(th/2) log n/n^3, and its log n/n^3 coefficient is -2 sin^2(th/2) a_2 != 0".
- Sketch: "Taylor's theorem gives ... (1+a_2 eps^2+O(eps^3)) uniformly" is fine for fixed th (uniform in n). It does not say uniform in th; if the paper ever wants that, it needs Lemma A.
- Sketch: "It is analytic and positive for 0 <= eps <= 1/2" is true; positivity comes from the Mellin integral of a positive function (first identity of thm:P2), not from the 2F1 series. Optional to say so.
- Paragraph after: "L(n)=1+(2/(en))^{1/(n-1)}+pi^2/(12n^2)+O(log n/n^3)" -- correct, and the log is genuine there too (at z = 1/2, -2z a_2 = -(pi^2/12 + 1/2) = -1.32247).
- "The terms of order n^{-3} and n^{-4} ... are explicit too." This proof covers the n^{-3} terms (log n and constant, given a_3) only. It does not prove anything about n^{-4} or about the size of the remainder after it. That sentence is not supported by the refereed proof; it needs its own source or should be softened.
