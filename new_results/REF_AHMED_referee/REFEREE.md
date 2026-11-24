# Referee report: four Ahmed results (thm:P2, thm:fn, thm:xnyn, thm:S), paper draft 36

I am an independent referee. I wrote my own code (all scripts are in this folder). I did not import any project
script. I read the paper `paper/ieee/main.tex` and the notes named in the brief. I did not modify anything outside
this folder. I do not have Ahmed's preprint, so every claim about what Ahmed wrote or asked is flagged, not checked.

## Summary

| Result | Verdict | My statement check (digits) | Remainder order |
|---|---|---|---|
| thm:P2 (Thm 15) exact length at every angle | **PROVED** | length: 50-61 digits (7 angles x 8 values of n, n = 1.3 to 1e8); Mellin identity: 48 digits (7 angles x 5 eps) | exact formula, no remainder |
| thm:fn (Thm 27) graph of f^n | **PROVED** (theorem). Two text repairs around it (below) | D_n = 2 - L_n to 60+ digits (self-check 2e-64), 5 functions, n = 1e3..1e10 | O(log^2 n/n^2) confirmed; ratio tends to -(a^2-f''(1))/(2a^3) |
| thm:xnyn (Thm 29) Ahmed's Open Problem 4 | **PROVED** (proof = notes + `verify_proof.md` closures) | A(n)-1 to 45 digits (self-check 7e-47), n = 2..1e10, direct from the defining integral | n^4 * error / (1+log n)^4 stays in [-0.28, -0.088] for 2 <= n <= 1e10; tends to -1/8 |
| thm:S (Thm 30) superellipsoid | **PROVED** (Q0 closed form is computer-assisted; I checked it to 35 digits, not the certificate) | exact formula vs direct surface integral: 30-31 digits at p = 1, 2, 2.5, 3, 4, 10; Q0: 35; Q1: 37; s1, s2: 24-26 (limited by printing), s3: 24 | p^4 (S - 24 - s1/p - s2/p^2 - s3/p^3) = -77.17, -77.76, -77.98, -78.04 at p = 100..3000 (tends to s4 = -78.07) |

---

## 1. thm:P2 (Theorem 15, exact length at every angle)

### (a) Statement check (`s1_P2.py`, `s1b_mellin.py`)
* **Length.** I computed L_theta(n) = int_0^1 |u + n x^{n-1} v| dx directly (variable tau = -log x, breakpoints
  around the slope-one point), at 60 digits. I compared it with the printed formula (2F1 for Lambda_eps,
  Legendre recurrence for Lambda_m, the full series). Grid: theta = 0.3, pi/3, pi/2, 2pi/3, 3pi/4, 2.9, 3.1 and
  n = 1.3, 1.7, 2.5, 3, 7.5, 100, 1e4, 1e8. Worst |difference| 1.75e-51 (at pi/2, n = 1e8, the quadrature limit);
  typical 1e-61. So the formula holds also for 1 < n < 2 (eps > 1), as claimed.
* **Removable singularities** n = 2 and n = 3/2: formula at n0 +- 1e-30 differs from the direct value by +-3e-31,
  antisymmetric, so the limit exists and equals L(n0).
* **Mellin identity** (first display of the theorem): 48 digits at 35 points (7 angles, eps = 0.05..0.97).
  Note: the first block of `s1_P2.out` ("Mellin identity") is invalid. It used r - s - c directly and lost all
  digits to cancellation at large s. `s1b_mellin.py` redoes it with k_c = sin^2(theta)/(s + c + r) and with the
  endpoint singularities integrated exactly; that is the valid check.
* Lambda_eps = int_c^1 P_eps = 2 sin^2(theta/2) 2F1(...): 62 digits. Lambda_1 = sin^2(theta)/2: exact.
  |Lambda_m| <= 1 - c: holds for m <= 58 at two angles.
* At theta = pi/2 the formula equals 1 + R(eps) n^{-eps} - sum_k binom(1/2,k)/(n^{2k-1}((2k-1)(n-1)-1)) to 1e-61.
  That is the form the notes attribute to Ahmed. I cannot confirm it is Ahmed's printed formula.

### (b) Proof check (notes `P2_general_angle/claude/result.md`, Thm 1, Thm 2)
* Lemma 0 (forms of Lambda): steps 1-4 re-derived by hand. Correct.
* Thm 1, step 1 (two integrations by parts): re-derived. Boundary terms vanish for 0 < eps < 1 as stated.
  phi'' = -int s^{eps+1} r^{-3} also at c = 1 (convergent for eps < 1). ODE (1 - c^2) Phi'' + eps(eps+1) Phi = 0 correct.
* Step 3 data at c = 1 and step 4 (Psi solves the same problem) correct.
* Step 5 (Gronwall at the singular end c = 1): correct, because |D| <= (1 - c) m and 1 - c^2 = (1 - c)(1 + c) cancel.
* Thm 2: Gegenbauer expansion k_c(s) = sum_{m>=2} C_m^{(-1/2)}(-c) s^{1-m} and C_m^{(-1/2)}(-c) = (-1)^m Lambda_{m-1}
  re-derived (uses int_{-1}^1 P_{m-1} = 0 for m >= 2). Term-by-term integral and re-indexing re-derived; the sign
  of the series is right. Extension to 1 < n <= 2 by the identity theorem: correct (radicand zeros at y = -e^{+-i theta}
  stay off the path for n near the real axis when theta < pi; the right side is meromorphic in Re n > 1).
* Everything is by hand or exact; nothing rests on numerics. **No gap.**

### (c) Wording
* "At theta = pi/2 this is Ahmed's exact formula." and "He found the exact length L_{pi/2}(n) through a Gamma
  function", "expecting Legendre functions": these depend on Ahmed's text. My check only shows agreement with the
  notes' transcription. Keep, but someone with the preprint should confirm once.
* The citation "[GRtab, sect. 3.25]" for the convergent Mellin transforms: not checked (no access).

---

## 2. thm:fn (Theorem 27, length of the graph of f^n)

### (a) Statement check (`s2_fn.py`, `s2b_slope.py`)
* I computed D_n = 2 - L_n = int_0^1 [1 - 1/(y' + sqrt(1 + y'^2))] dx (no cancellation), y = f^n, directly in
  tau = 1 - x, at 60 digits (80-digit rerun agrees to 2e-64). Functions: x, (e^x - 1)/(e - 1), 2x/(1 + x),
  sin(pi x/4)/sin(pi/4), (x + x^3)/2. n = 1e3, 1e4, 1e5, 1e6, 1e8, 1e10.
* (D_n - printed main term)/(log^2 n/n^2) is bounded and drifts to -(a^2 - f''(1))/(2a^3), the predicted
  log^2 coefficient. Example 2x/(1+x): -2.24, -2.41, -2.51, -2.58, -2.68, -2.74 -> -3. So the constant
  log f'(1) + 1 - log 2 and the remainder order O(log^2 n/n^2) are confirmed.
* Next orders (text claim "we also prove the next order"): with J_1, J_2, J_3 of Sec. sec:fn at z = 1/2, c = 0,
  (D_n - eps J1 - eps^2 J2 - eps^3 J3)/(eps^4 log^4 n) stays bounded for all five functions (e.g. 2x/(1+x):
  -1.96 ... -4.37, slowly approaching the predicted top coefficient -6.25 that I computed from theta_40 of the notes).
* Cumulants of log(1/S): k1, k2, k3, k4 confirmed to 32, 30, 28, 26 digits; geometric mean of S = 2/e to 32 digits.
* Slope-2/e law (notes T3; only named in the paper): (L_n - 1 - x*) n^2 -> c_2(f) for three functions
  (e.g. f = x: 0.822467 = pi^2/12), and the residual times n^3/log n is bounded.
* Contact law at theta = pi/2, f = sin(pi x/2): D/(log n/(m beta n))^{1/m} = 1.18, 1.14, 1.12, 1.09, 1.07 at
  n = 1e4..1e16, and D divided by the 1/Lambda-corrected form of the notes = 0.986 ... 0.9993. Consistent.
* Hypothesis probe: a C^1 (not C^2) kink at x = 1/2, C^2 near 1: the ratio stays bounded (-0.45, -0.43, -0.42), as
  the proof (only the germ at 1 matters) predicts.

### (b) Proof check (notes `P3_fn_graph/claude/result.md`, Prop 1, Prop 2, Thm 1, Thm 2)
* Prop 1 (layer cake): by hand. Correct for y in C^1 (the proof uses int y' = 1).
* Prop 2: nondecreasing f with f(1) = 1 and f'(1) > 0 forces f < 1 on [0,1) (else f = 1 on an interval ending at 1).
  So theta = f(1 - eta) < 1 and the bound y' <= n theta^{n-1} ||f'|| on [0, 1 - eta] holds. (ln y')' > 0 on
  [1 - eta, 1] uses f'' there: covered by "C^2 near 1". (iii) correct.
* Thm 1: the bound tau = mu(lambda + ln a)/a + O(mu^2 (1 + |lambda|)^2) is uniform also when lambda is of order n
  (s near sigma_n), because tau <= eta and A(tau) = a tau + O(tau^2) holds on all of [0, eta]. Integration against g':
  the four bullets re-derived; tail bound mu/(2 a^2 n) correct (g' <= 1/(2 s^2)). E ln S = ln 2 - 1 checked.
* Hypotheses of the printed theorem (C^1 on [0,1], nondecreasing, f(0) = 0, f(1) = 1, f'(1) > 0, C^2 near 1)
  are exactly what the proof uses (f(0) = 0 is not even needed). **No gap.**
* Uniformity: constants depend on f (eta, n_0, ||f'||, ||f''|| near 1). The statement claims nothing uniform in f.
  Fine.
* Thm 2 (all orders; quoted in the text after the theorem): Theta is C^{K+1} if f is C^{K+2}; the Taylor bound is
  valid on the whole compact range of W_s (not only near 0); the T-part bookkeeping (i + j <= k - 1) re-derived.
  Dependence on f'(1)..f^{(k)}(1) correct.

### (c) Wording
1. Line 1764: "the laws are stated in Table~\ref{tab:verify}". **False.** The table only names them
   ("expansion, slope law, contact law"). Suggested text: state them, e.g.
   "There is a unique x_n^* with (f^n)'(x_n^*) = 2/e, and L_n = 1 + x_n^* + c_2(f)/n^2 + O(log n/n^3) with
   c_2(f) = ((a^2 - f''(1)) pi^2/12 + (a^2 - a - f''(1))/2)/a^3, a = f'(1)." and the contact law as in item 2.
2. Line 1802 (contact law): "When f touches 1 with higher contact, -log f(1-tau) ~ beta tau^m, the deficit changes
   regime ..." The notes (Thm 4) prove it only under Phi(tau) = beta tau^m (1 + O(tau^delta)) **and**
   Phi'(tau) = m beta tau^{m-1}(1 + O(tau^delta)), m > 0, m != 1. A bare "~" does not control the slope (f^n)'.
   Suggested: "When f touches 1 with contact of order m != 1, that is -log f(1-tau) = beta tau^m (1 + O(tau^delta))
   with the same relation differentiated once (true for f in C^{m+1} when m >= 2 is an integer), ...".
   This sentence belongs to the every-angle discussion (thm:fnangle), but the pi/2 case is cited for thm:fn in
   the table.
3. Line 1742: "Ahmed conjectured yes: for any smooth increasing f ... 2 - L_n ~ log n/(f'(1) n); he asked for a proof
   and the next term". Depends on Ahmed's text. What is answered, as printed here: the conjectured asymptotic, under
   weaker hypotheses (C^1, C^2 near 1), plus the O(1/n) constant (log f'(1) + 1 - log 2)/f'(1) and an
   O(log^2 n/n^2) error; the full expansion is in the text after the theorem. If Ahmed's "next term" means the
   n^{-2} log^2 n order, then the theorem itself does not give it; the text does (J_2). Suggest: "This proves
   Ahmed's conjecture and gives the next term" only if "next term" means the constant.
4. Optional: "For any nondecreasing y with y(1) - y(0) = 1" -> "For any absolutely continuous nondecreasing y ...".

---

## 3. thm:xnyn (Theorem 29, Ahmed's Open Problem 4)

### (a) Statement check (`s3_xnyn.py`)
* A(n) - 1 = iint (sqrt(1 + G^2) - 1) dx dy from the definition, in x = e^{-s/n}, y = e^{-t/n}, sigma = s + t,
  delta = s - t, written as G^2/(1 + sqrt(1 + G^2)). 45 digits (60-digit rerun agrees to 7e-47 relative).
  No decomposition, no Mellin transforms.
* Error r(n) = A(n) - [printed expansion through R_3]: n^4 r/(1 + log n)^4 = -0.279, -0.200, -0.149, -0.117, -0.096,
  -0.089, -0.088, -0.091, -0.094, -0.097, -0.102, -0.105 at n = 2, 3, 5, 10, 30, 100, 1e3, 1e4, 1e5, 1e6, 1e8, 1e10.
  Bounded for all real n >= 2, as claimed. n^4 r/log^4 n -> -1/8 (the leading -log^4 n/8 of R_4).
  This pins every printed constant: an error of size delta in d3 would make n^4 r grow like n delta.
* D(1e5) = -76.7411466, R + R_3/n = -76.7411464: matches the printed "-76.74" and "-76.7411".
* Printed decimals of beta1, beta0, d3 agree with the closed forms (truncated, not rounded; fine).

### (b) Proof check
* I read `result.md` and `verify_proof.md`. `result.md` alone has the two real gaps that `verify_proof.md` found
  (G1: no uniform bound on the Delta_1 Taylor remainder; G3: e^2 coefficient not derived). `verify_proof.md`
  closes them with Lemmas 4-6. I re-checked the key inequalities by hand:
  kappa = e sigma/2 + (1/2) log cosh(e delta) <= e sigma (so e^{-e sigma} e^kappa <= 1);
  psi(w) = w f'(w) <= min(w, 1/(2w)); psi(u) <= e^kappa min(W0, 1/(2 W0)) for u in [W0, W0 e^kappa].
  The Mellin transform of phi(w)/w equals Ahmed's Gamma factor (as transcribed in the notes) at eps = z - 1,
  divided by 1 - z: checked by hand.
* The paper's proof sketch matches this route. **No gap, provided "in our notes" points to the closed proof in
  `verify_proof.md` §2 (Lemmas 3-6), not to `result.md` §2 alone.** `verify_proof.md` also flags a false relation
  in `result.md` §3 (G5); it is not used by the paper.

### (c) Wording
* Line 1811: "Yes, and one more order as well. This solves Ahmed's Open Problem 4." What is proved: beta_1, beta_0 in
  closed form, the n^{-3} term R_3 with d_3 (containing -zeta(3)/2), and an O((1 + log n)^4/n^4) error valid for all
  real n >= 2. Whether this is all of "Problem 4" depends on Ahmed's wording (the paper quotes it as "how does the
  area approach its limit" and then as "find beta_1 and beta_0"). Suggested: "This answers Ahmed's Open Problem 4
  (the constants beta_1 and beta_0), and gives one more order."
* Line 1809: "-(1/2 log^2 n + O(log n log log n))/n^2" is true but odd (the error is O(log n)). Harmless.

---

## 4. thm:S (Theorem 30, superellipsoid)

### (a) Statement check (`s4_S.py`, `s4b_constants.py`)
* **Exact formula.** Direct surface area from the definition (graph over the part of an octant where z is the
  largest coordinate and y <= x, times 48; integrand sqrt(1 + (x/z)^{2p-2} + (y/z)^{2p-2})) against
  6 G(P - 4) - 8 eps^2 Q(eps), with P from the arc-length integral and Q from the printed double integral
  (h_eps by quadrature, then int h^3). Relative differences: p = 1: 2.3e-31 (vs 4 sqrt3), p = 2: 1.3e-31 (vs 4 pi),
  p = 2.5: < 1e-30, p = 3: 2.0e-31, p = 4: 3.6e-31, p = 10: 1.5e-31. So 30-31 digits.
* **Q0.** Q(0) from the printed integral = -2.83293568349484776046205725151577; the printed closed form agrees to
  1.6e-35 (35 digits).
* **Q1** = Q'(0), from differentiating the printed integral under the integral sign:
  Q1 = -(3/sqrt(pi)) int h_0^2 dh + 2 gamma Q0 = 6.50432557374100975876610779051. Agrees with the notes' 400-digit value
  to 7.6e-38. A crude central difference (delta = 1e-4) agrees to 7e-8, as expected.
* **M0, M1, M2** from the perimeter integral (derivatives under the integral): M0 and M1 agree with their closed forms
  to 31-32 digits, M2 with the printed -0.02216191108726939942 to all 20 printed digits.
* **s1, s2, s3** from these numbers by the printed formulas: s1 = -26.558558386584754433798786,
  s2 = -4.6892155184605022499600570, s3 = 48.287432486014447759215859. All printed decimals agree.
  The table formula (T6) for s2 equals 48 M1 - 4 pi^2 - 8 Q0 (closed forms) to 7e-50.
  The algebra s_k from G = 1 - zeta(2) eps^2 + 2 zeta(3) eps^3 + ... and P = 8(1 + M0 eps + ...) re-derived by hand.
* **Expansion and remainder order.** R4(p) = p^4 (S(p) - 24 - s1/p - s2/p^2 - s3/p^3), with S(p) from the exact
  formula at 35 digits: -77.168, -77.765, -77.976, -78.036 at p = 100, 300, 1000, 3000. Bounded and converging to
  s4 = -78.066 like s5/p (s5 about 90). So O(p^{-4}) after s3 holds, and no log p term appears through order p^{-3}.

### (b) Proof check (notes `S3_superellipsoid/claude/result.md`)
* Area element: re-derived (dA = eps^2 (uvw)^{eps-1} ||X|| du dv with x = u^eps). Correct.
* Gaussian linearisation ||X|| = (1/(2 sqrt pi)) int (1 - e^{-t||X||^2}) t^{-3/2} dt: correct. Inclusion-exclusion and
  Tonelli for each term (each term is finite separately): correct.
* Laplace-Dirichlet lemma: I re-derived the scaling phi_i(tau x; t) = tau^{eps-1} phi_i(x; t tau^{2-2eps}), the factor
  tau^{3eps-1}, and the factor Gamma(1 + 2 eps) after the t-integration. H_1 = 1 re-derived. Correct.
* No-log claim: H_2, H_3 are C^infinity at eps = 0 by dominated differentiation; the final formula has no 1/eps.
  Correct (and confirmed numerically above).
* Q0 closed form: one-dimensional reduction (Step 1) is by hand; the reduction of 54 dilogarithm values (Step 3) is an
  exact computer-algebra certificate (`13_proof_Q0.py`), which I did not re-run. My 35-digit check supports it.
* **No gap** in the theorem as printed. Q1 (hence the value of s3) is numerical, and the paper says so.

### (c) Wording
* Line 1877: "The closed form of the vertex constant Q0 is proved by an exact reduction ..." Suggest "is proved, with
  computer assistance, by an exact reduction ...", and "Proved (computer-assisted)$^\dagger$" in the table row.
* Theorem statement: "a, b, L, lambda as in Section~\ref{sec:square}" - a is not used in the theorem; harmless.

---

## Files
* `s1_P2.py` / `.out` (length, removable points, Ahmed form at pi/2; ignore its first, invalid Mellin block),
  `s1b_mellin.py` / `.out` (valid Mellin identity check).
* `s2_fn.py` / `.out` (Thm 27, J_2, J_3, cumulants, contact law), `s2b_slope.py` / `.out` (slope-2/e law).
* `s3_xnyn.py` / `.out` (Thm 29).
* `s4_S.py` / `.out` (Thm 30; its closed-form comparisons used 15-digit constants and are superseded by)
  `s4b_constants.py` / `.out` (closed forms at 50 digits against the numbers of `s4_S.out`).
* `*.log` are the stdout copies.
