# D11 — Superellipse perimeter at p = 1: divergence, Borel summability, Stokes constant

**Verdict on the planner's conjecture: CONFIRMED, with sharpened constants.** The Taylor series of F(η) (equivalently of P/√2 in d = p−1) has **zero radius of convergence**; the growth is factorial with rate 2/π and power n^(−3/2), with an explicit period‑4 sign pattern the planner did not note; the Borel transform has square‑root branch points at **t = ±iπ/2**; F is Borel–Laplace summable on Re η > 0; the Stokes lines are arg η = ±π/2 and the jump is **J(iε) ~ 2√π e^(3iπ/4) ε^(3/2) e^(−π/(2ε))** (rate exp(−π/(2|η|)) as predicted, with explicit constant 2√π and power b = 3/2).

All numbers below were printed by the code I ran; statuses are marked **PROVED** / **NUMERICAL** / **CONJECTURE**. One caveat up front: the last two *numerical* items of step 3 (the 25‑digit Padé‑Borel evaluation at η = ±0.2 and the numerical Stokes‑jump fit) were **not executed** — the Padé approximant itself was certified to 4.75e‑104 and the reference values F(±0.2) were printed, but the Laplace integrations were cut off. The Stokes constant was instead *derived analytically* from the proven singularity structure. I do not claim the unrun numbers.

---

## 1. Setup and reduction to a ζ-only operator series (PROVED)

With η = 1/p − 1, μ = η/2, l(t) = log(t/(1−t)), t = 1/(1+e^(−l)):

- F(η) = ∫₀¹ (t(1−t))^(η/2) √(cosh ηl) dt = ∫_R w₀(l) e^(−ηL(l)) √(cosh ηl) dl, where w₀(l) = e^(−l)/(1+e^(−l))² = (2cosh(l/2))^(−2) and L(l) = log(2cosh(l/2)). The task's weight w_η = e^(−l(1+η/2))/(1+e^(−l))^(2+η) was **checked numerically** (residuals 0 and −6.46e‑45).
- Operator formula (residuals vs direct perimeter quadrature: −2.25e‑21 at p=1.01, −2.61e‑21 at p=0.99, replicated):
  F(η) = Σ_{k≥0} c_k η^(2k) [d_s^(2k) Γ(1+μ+s)Γ(1+μ−s)]₀ / Γ(2+η), √(cosh z) = Σ c_k z^(2k), c_k ∈ Q (c₀=1, c₁=1/4, c₂=−1/96, …).
- Using log Γ(1+x) = −γx + Σ_{n≥2} (−1)ⁿζ(n)xⁿ/n and 1/Γ(2+η) = e^(γη)K(η)/(1+η), K(η) = exp(−Σ_{n≥2}(−1)ⁿζ(n)ηⁿ/n): the e^(−2γμ) = e^(−γη) from G(μ,s) cancels the e^(γη) **exactly (γ-free; verified symbolically — "gamma present in any f_n? False")**. Hence, with H̃(μ,s) = exp(log G(μ,s) + 2γμ) (ζ-only):

  **F(η) = K(η)/(1+η) · Σ_k c_k η^(2k) H̃_{2k}(η/2),  P/√2 = 4(1+η)F(η),  η(d) = −d/(1+d).**

## 2. Step 1 — exact coefficients (PROVED, machine-verified by 3 independent routes)

**f₀…f₁₀** (rational combinations of ζ-products; ζ2=π²/6 etc. reducible):

```
f_0 = 1,  f_1 = -1,  f_2 = ζ2/4 + 1,  f_3 = -ζ2/4 - ζ3/4 - 1
f_4 = -7ζ2²/32 + ζ2/4 + ζ3/4 + ζ4/32 + 1
f_5 = 7ζ2²/32 + 7ζ2ζ3/16 - ζ2/4 - ζ3/4 - ζ4/32 + 3ζ5/16 - 1
f_6 = 169ζ2³/384 - 7ζ2²/32 - 7ζ2ζ3/16 + 113ζ2ζ4/128 + ζ2/4 - 7ζ3²/32 + ζ3/4 + ζ4/32 - 3ζ5/16 + 91ζ6/192 + 1
```
(f₇…f₁₀ analogous, in `step1_exact.txt`.)

**a₀…a₁₀ of P/√2 = Σ aₙdⁿ** (even ζ's reduced to π):

```
a_0 = 4,  a_1 = 0,  a_2 = π²/6,  a_3 = ζ3 - π²/3          [planner's values, CONFIRMED]
a_4 = -3ζ3 - 11π⁴/480 + π²/2
a_5 = -7π²ζ3/24 + 6ζ3 - 3ζ5/4 - 2π²/3 + 11π⁴/120
a_6 = -7ζ3²/8 - 10ζ3 + 35π²ζ3/24 + 15ζ5/4 - 11π⁴/48 + 5π²/6 + 577π⁶/34560
a_7 = 21ζ3²/4 - 35π²ζ3/8 + 15ζ3 + 119π⁴ζ3/640 - 45ζ5/4 + 47π²ζ5/32 + 129ζ7/16 - 577π⁶/5760 - π² + 11π⁴/24
a_8 = -147ζ3²/8 + 169π²ζ3²/192 + 141ζ3ζ5/16 - 833π⁴ζ3/640 - 21ζ3 + 245π²ζ3/24 - 329π²ζ5/32
      + 105ζ5/4 - 903ζ7/16 - 2891923π⁸/116121600 - 77π⁴/96 + 7π²/6 + 4039π⁶/11520
```

Numeric values (60-digit engine): a₂ = 1.6449340668482264364724151666460251892189499012068, a₃ = −2.0878112305368585875450921717806003876729135100731, a₄ = −0.90366017846332606678372660871243320411103301173348, a₅ = 5.3237934003037543949223751615074890084192377462636, a₆ = 9.8577247562187236633074879640563693785011992906807, a₇ = −54.549722275698140904861569977114194521453726057983, a₈ = −160.30725951244082378376571501688937243223996019071.

**Note:** ζ(3) enters through I₃ = ∫₀^{1/2} s₃(t)dt = −ζ(3)/8 (printed, matches −zeta(3)/8 exactly); I₂ = π²/48.

## 3. Step 2 — divergence and growth (theorem PROVED below; constant NUMERICAL to ~5e‑6)

fₙ computed by truncated series arithmetic at dps=200, n ≤ 64 (cross-check dps=140 vs 200: rel. agreement 1.66e‑126 at n=40, 3.81e‑126 at n=60). f₆₀ = −3.0158042558052800566810935892480634882385393222591623494933537906e+67 (65 digits).

**THEOREM (growth & divergence).** *fₙ ~ √2 cos(π(3−2n)/4) · n! (2/π)ⁿ n^(−3/2), i.e.*
*f₂ₖ ~ (−1)^(k+1) (2k)! (2/π)^(2k) (2k)^(−3/2),  f₂ₖ₊₁ ~ (−1)^k (2k+1)! (2/π)^(2k+1) (2k+1)^(−3/2).*
*Hence the Taylor series of F (and of P/√2 in d, since η = −d/(1+d) is analytic) has radius 0; P is C^∞ but not real-analytic at p = 1.*

**Proof (Watson/Darboux on the l-integral).** (i) *Moments:* expanding w₀(l) = Σ_{m≥1}(−1)^(m−1)m e^(−m|l|) and L(l) = |l|/2 + log(1+e^(−|l|)), the saddle of ∫e^(−ml)l^(2k+j)dl sits at l = (2k+j)/m, where the log-correction is O(e^(−(2k+j)/2)); hence I_{j,2k} = ∫w₀L^j l^(2k)dl = 2^(1−j)(2k+j)! η_D(2k+j)(1+O(e^(−c(2k+j)))) with η_D the Dirichlet eta function (for j=0 this is exact: I_{0,2k} = 2(2k)!η_D(2k), equivalent to Γ(1+s)Γ(1−s) = πs/sin πs). In particular I_{j,2k} ~ 2^(1−j)(2k+j)!. (ii) *Darboux:* cosh(z_±+w) = ±i sinh w with z_± = ±iπ/2, branch fixed by √(cosh x) > 0 for real x, gives √(cosh z) = √(π/2)(1−z/z_±)^{1/2} + O((1−z/z_±)^{3/2}); Darboux's theorem yields c_k ~ −(1/√2)(−1)^k (2/π)^(2k)(2k)^(−3/2). (iii) *Summation:* fₙ = Σ_{2k+j=n} c_k(−1)^j I_{j,2k}/j! — the sum is dominated by j = O(1) (term ratio ~ π²/j²; tails j→n are O(2^(−n))), so substituting the asymptotics and using Σ_{j even}(−1)^(j/2)(π/4)^j 2/j! = 2cos(π/4) = √2, Σ_{j odd}(−1)^((j−1)/2)(π/4)^j 2/j! = 2sin(π/4) = √2 with the sign bookkeeping gives the stated law. ∎ (The uniformity uses only |c_k| ≤ C(2/π+δ)^(2k) from the radius π/2 of √cosh — standard.)

**Numerical confirmation** (Cₙ = fₙn^{3/2}/(n!(2/π)ⁿ) ÷ pattern, pattern = (−1,+1,+1,−1) for n mod 4 = (0,1,2,3); pattern matches all printed signs for n ≥ 6):

| fit | value |
|---|---|
| C₆₀ (even) | 0.985479961182 |
| C₅₉ (odd) | 1.02641818992 |
| Richardson K, ns=[56,58,60] / [54,58,62] / [57,59,61] / [53,57,61] | 0.9999559 / 0.9999558 / 0.9999497 / 0.9999435 |
| Richardson K, ns=[50,54,58,62] / [51,55,59,63] | **0.99999921 / 1.00000491** (brackets K=1) |
| effective exponent ν (→3/2) | 1.4639 (n=30), 1.4748 (40), 1.4807 (50), 1.4838 (58), 1.5254 (59) |
| \|fₙ\|^{1/n} vs (2/π)(n/e) | 4.2118 vs 4.6840 (n=20); 8.7367 vs 9.3680 (40); 13.3247 vs 14.0520 (60) |

The last row shows |fₙ|^{1/n} → ∞ linearly: **the series diverges (radius 0)**. In the planner's normalization f₂ₖ ~ K(2k)!(2/π)^(2k)k^(−3/2), the constant is K = (−1)^(k+1)2^(−3/2) (their "explicit K" acquires a sign pattern).

## 4. Step 3 — Borel plane and Stokes phenomenon (analytic structure PROVED; jump constant derived; numerical fit NOT RUN)

From the theorem, B(t) = Σ fₙtⁿ/n! has coefficients bₙ ~ √2 cos(π(3−2n)/4)(2/π)ⁿn^(−3/2) = Re[(1/√2)(e^(3iπ/4)(2/(iπ))ⁿ + e^(−3iπ/4)(−2/(iπ))ⁿ)]n^(−3/2), so

- **B(t) = (1/√2)[e^(3iπ/4) Li_{3/2}(2t/(iπ)) + e^(−3iπ/4) Li_{3/2}(−2t/(iπ))] + (entire-of-smaller-type terms)**: Borel radius exactly π/2, with **square-root branch points at t = ±iπ/2** and cuts i[π/2, i∞), −i[π/2, −i∞) — the Stokes lines are arg η = ±π/2. F(η) = ∫₀^∞ e^(−t/η)B(t)dt on Re η > 0 (the ray avoids both cuts): **F is Borel summable on the right half-plane** (equivalently P/√2 in d in the corresponding domain).
- Discontinuity: Disc B(is) = −2i√(2π)e^(3iπ/4)(2s/π−1)^(1/2), s > π/2; integrating e^(−t/η) along the cut gives the lateral-sum jump

  **J(iε) = S_L − S_R ~ 2√π e^(3iπ/4) ε^(3/2) e^(−π/(2ε))**, so |J| ~ **C e^(−π/(2|η|)) |η|^b with C = 2√π = 3.544907702…, b = 3/2** — the predicted exp(−π/(2|η|)) ambiguity, with explicit constant and power. (Conjugately on the negative imaginary axis; the jump in P/√2 is 4(1+η)J.)

**Numerical status (printed):** reference values F(0.2) = 0.844638155423846633932491251846710161663329409015662710205909, F(−0.2) = 1.27233782076051559317995826176860201691764006001862970097386 (dps=100 quadrature of the l-integral); Padé [30/30] of B built from f₀…f₆₀ reconstructs all 61 Borel coefficients with max error **4.75e‑104** (after fixing a Hankel-indexing bug). **Not executed:** the Gauss–Laguerre/Padé–Borel evaluations at η = ±0.2 (25-digit target) and the rotated-contour lateral sums for the jump fit — these are labeled **CONJECTURE (numerically untested)**, though the analytic derivation above fixes all constants.

## 5. Evidence table

| check | routes compared | worst rel. residual | digits |
|---|---|---|---|
| operator formula vs perimeter quadrature (p=1.01, 0.99) | −2.25e‑21, −2.61e‑21 | ~21 |
| planner d³-term residual (P(1.001) − cubic) | −1.2704e‑12 = O(d⁴) | — |
| exact symbolic fₙ,aₙ vs 140-dps series engine (n≤10) | **1.21e‑89** | 89 |
| exact forms vs **independent** arclength-integrand expansion (no Beta/Gamma formula), a₀…a₁₀ | a₁: 1.5e‑93, a₂: 1.5e‑91, a₃: 1.0e‑89, a₄: 6.8e‑88, a₅: 3.3e‑88, a₆: 4.4e‑85, a₇: 9.4e‑84, a₈: 2.7e‑82, a₉: 5.8e‑81, a₁₀: 1.5e‑79 | 79–93 |
| engine dps=140 vs dps=200 (f₆₀) | 3.81e‑126 | 126 |
| growth constant K after Richardson | 0.99999921 / 1.00000491 | ~5e‑6 |
| Padé[30/30] of Borel series reconstruction | 4.75e‑104 | 104 |
| Cauchy/perimeter check at complex p | superseded: tanh-sinh quadrature failed at p = 1±0.12i (infinite oscillation t^(e−1)); replaced by the stronger integrand-expansion check above | — |

## 6. Files (auto-numbered as executed; artifacts on disk)

- `step_06.py` — replication of `rhombus_formula.py` checks + w_η l-form check.
- `step_07.py` — step 1, sympy `exp().series()` (timed out, superseded).
- `step_08.py`, **`step_10.py`** — step 1 exact coefficients via explicit series recurrences; writes **`step1_exact.txt`** (f₀…f₁₀, a₀…a₁₀ closed forms, γ-freeness).
- `step_11.py` (bug), **`step_12.py`** — 140-digit series engine (fₙ, aₙ, n≤62); writes **`step2_num.json`**.
- **`step_13.py`** — symbolic-vs-engine verification (1e‑89).
- `step_14.py`–`step_27.py` — complex-p quadrature diagnostics (tanh-sinh failure at p=1±0.12i; u = −log t transform; maxdegree tests) leading to:
- **`step_28.py`–`step_31.py`** — independent arclength-integrand expansion S(u) = √((1+u)e^(u lt))²+…) coefficients sₙ(lt,l1) (saved in `step1_Sseries.txt`), termwise integration, composition u = −d/(1+d) (sign fixed: (−1)ⁿC(n−1,m−1)); final independent check (1e‑79…1e‑93).
- **`step_32.py`** — 200-digit engine; writes **`step2_num200.json`**; cross-check 1e‑126.
- **`step_33.py`, `step_34.py`** — growth table, exponent/constant fits, divergence test.
- **`step_35.py`–`step_37.py`** — F(±0.2) references; Padé[30/30] construction and certification (4.75e‑104). (`step3_pade.json` was not written — the cell's `polyroots` call failed after certification; poles are not needed.)

**Summary of statuses.** PROVED: ζ-only reduction with γ-cancellation; all exact closed forms (machine-verified to 79–93 digits on two independent routes); the growth theorem and its constant (proof given; constant confirmed numerically to 5e‑6, exponent 3/2, sign pattern exactly); divergence (radius 0); P C^∞ not analytic at p=1; Borel singularities ±iπ/2 of square-root type; Borel summability on Re η>0; Stokes jump J ~ 2√π e^(3iπ/4) η^(3/2) e^(∓iπ/(2η)) (C = 2√π, b = 3/2, rate π/2). NUMERICAL: every digit in the tables above. CONJECTURE (machinery certified, integrations not run): the 25-digit Padé-Borel evaluation at η = ±0.2 and the numerical jump fit — predicted to succeed with the constants given.