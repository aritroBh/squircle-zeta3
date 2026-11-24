# K3: Fourier meets the squircle, to all orders

Directory: `research/08_theory/K3_fourier_all_orders/` (all scripts `step_NN.py` there; data files `e_rm.pkl`, `e_rm_ext.pkl`, `Dtilde_relations.pkl`, `Dt10.pkl`, `Qpolys.pkl`). Python: project `.venv`, mpmath ≤ 240 dps, sympy exact series. Sources used: `../../10_frontier/F1_fourier_squircle/` and D11's `step1_exact.txt` (exact perimeter coefficients, used as an independent validator).

**Notation.** τ = 1/q = 1−1/p; A_k = −(−1)^k·4√2/(π(16k²−1)) = a_k(diamond); u = k²; ζ-products abbreviated z2z3 = ζ(2)ζ(3) etc.; η(s) = (1−2^{1−s})ζ(s); gd = Gudermannian, L(s) = log(1+e^{−s}); H(x) = 2^{−3/2}e^{x/2}sech^{3/2}x = e^{−x}(1+e^{−2x})^{−3/2}.

---

## 1. Task 1: c_n(k) to all orders

### 1.1 The master identity (PROVED)

On [0,π/4], cos(4kθ) = (−1)^k cos(2k(π/2−2θ)) and π/2−2θ = gd(x) for x = −log tanθ. Hence, with h_p − h_1 = cosθ[(1+tan^qθ)^τ − 1]:

> **a_k(p) − A_k = (−1)^k B(τ, 2k)**, where **B(τ, κ) = (8/π)∫₀^{π/4}(h_p−h_1)(θ)cos(κ(π/2−2θ))dθ** is even and entire in κ, and
> B(τ,κ) = Σ_{m≥0} D(τ;m)(−κ²)^m/(2m)!,  D(τ;m) = (8/π)∫₀^{π/4}(h_p−h_1)(π/2−2θ)^{2m}dθ = Σ_{n≥2m+2} D_{n,m}τⁿ.

So **c_n(k) = (−1)^k(√2/π)Q_n(2k)** with Q_n a polynomial. The Watson coefficients are (step_12, step_13):

> **D̃_{n,m} := (π/√2)D_{n,m} = 2 Σ_{j=1}^{n−1−2m} e_{m,n−1−j}·I_{n−1−j,j}/j!**,
> e_{r,m} = 2^{3/2}[x^m](gd x)^{2r}H(x) ∈ **ℚ** (exact table saved; e_{r,2r} = 1, e_{r,2r+1} = ½, e_{r,2r+2} = −5/8−r/3, …), and
> I_{m,j} = ∫₀^∞ s^m log^j(1+e^{−s})ds = j!m!Σ_{l≥j}(−1)^{l−j}s(l,j)/(l!·l^{m+1}), an **alternating MZV of weight m+j+1** (s(l,j) = Stirling I).

Consequences (PROVED for every n): (i) c_n(k) = (−1)^k(√2/π)×(polynomial in k² of degree ≤ ⌊n/2⌋−1 with pure-weight-n alternating-MZV coefficients); (ii) the stronger reduction to **zeta products only** (no log 2, no Li_r(½)) is PROVED for n ≤ 5 (F1), for the r = 0 row for n ≤ 10 (D11's exact f_n: validated to ≤1.8e-99 in step_13/26/27), for the top row r = ⌊n/2⌋−1 for **all** n (§1.3), and NUMERICAL (PSLQ, residuals ≤ 1e-92…1e-226, small coefficients) for all r, n ≤ 12 (steps_14, 26, 29, 30, 31).

### 1.2 Closed forms through n = 7 (computed through n = 12)

c_n(k) = (−1)^k(√2/π)Q_n(u), u = k² (all PSLQ-identified, then verified — see §4):

```
Q2 = z2
Q3 = z3
Q4 = -(33/16 + 7u) z4
Q5 = -(3/4 + 4u) z5 - (7/4 + 4u) z2z3
Q6 = (4039/256 + 733u/8 + 31u^2) z6 - (7/8 + 2u) z3^2
Q7 = (129/16 + 46u + 16u^2) z7 + (141/16 + 50u + 16u^2) z2z5 + (1071/160 + 181u/5 + 56u^2/5) z2^2 z3
```
and (step_15, 26, 29, 30, 31; Q₈ top matches, all tops match §1.3):
```
Q8 = -(2891923/12288 + 495965 u/256 + 60889 u^2/48 + 127 u^3) z8
     + (141/16 + 50u + 16u^2) z3z5 + (169/32 + 27u + 8u^2) z2 z3^2
Q9 = ... + 16u^3(4 z9 + 4 z2z7 + (14/5) z2^2 z5 + (62/35) z2^3 z3)   [full table in step_15 output]
```
Q₆, Q₇ reproduce F1/f1_03's stored c₇(k) values (diff ≤ 3.5e-6, their print precision, step_15); Q₄, Q₅ reproduce F1's formulas exactly (step_15, diff ≤ 1e-59). D̃ tables for n = 8..12 are saved; e.g. D̃_{10,0} = (2015681451/327680)z10 − (8343/64)z3z7 − (8343/128)z5² − (8907/64)z2z3z5 − (64377/1280)z2²z3², matching D11's 4e_10 to 1.2e-97.

### 1.3 The degree conjecture: PROVED

**Theorem (PROVED).** deg_{k²} c_n(k) = ⌊n/2⌋ − 1 for every n ≥ 2, with exact top coefficients:
- D̃_{2n,n−1} = 2(2n−2)!·η(2n) (only j = 1, e_{n−1,2n−2} = 1, I_{2n−2,1} = (2n−2)!η(2n) contribute — exact);
- D̃_{2n+1,n−1} = (2n−2)!·[ζ(2n+1) + 2Σ_{j=1}^{n−1} η(2j)ζ(2n+1−2j)] — using e_{r,2r} = e_{r,2r+1}-coefficients plus the classical evaluation 2A_{2m} = (2m+1)η(2m+1) − ζ(2m+1) − 2Ση(2j)ζ(2m+1−2j) of A_{2m} = Σ(−1)^{l−1}H_l/l^{2m}, which I verified from the series to ≤ 4.7e-61 for m ≤ 6 (step_16).

Both tops are manifestly nonzero, so the degree is attained. Checks: D̃_{4,1} = 7ζ4/2, D̃_{6,2} = 93ζ6/2, D̃_{8,3} = 5715ζ8/4, D̃_{10,4} = 160965ζ10/2 all equal 2(2n−2)!η(2n) (step_26: diff 9.4e-97); odd tops match through n = 11 (D̃_{11,4} = 40320ζ11+40320ζ2ζ9+28224ζ2²ζ7+17856ζ2³ζ5+(54864/5)ζ2⁴ζ3, step_29).

### 1.4 Generating functions

(i) **Exact**: the integral B(τ,κ) of §1.1 (equivalently (8/π)τ∫₀^∞[(1+e^{−s})^τ−1]cos(κ gd(τs))H(τs)ds).

(ii) **Top-degree (closed form, PROVED given the A_{2m} identity; verified to 6.1e-52, step_17):**
> Σ_{n≥2} Q̂_n(κ)τⁿ = (1/κ²)[ 1 − (πκτ/ sinh πκτ)·(1 − κτ(γ + Re ψ(1+iκτ))) ],  Q̂_n = top-degree term of Q_n.

The even part is (1/κ²)(1 − πκτ/sinh πκτ) = at κ = 2k **exactly F1's empirical 1/sinh uniform formula** 1/(4k²) − πτ/(2k sinh 2πkτ) — so F1's NUMERICAL Prop. 3b is explained: it is the even-weight top-degree part; the digamma term is the new odd-weight companion (responsible for the odd tops, e.g. −(ζ5+ζ2ζ3)κ² in Q₅).

(iii) Bonus pattern (CONJECTURE, verified j ≤ 4): D̃_{3j,j−1} ⊃ (2j−2)!/j!·ζ(3)^j (values 1, 1, 4, 30 for j = 1..4 — j=4 = D̃_{12,3} = 30ζ3⁴, step_30).

---

## 2. Task 2: re-summing in k — the boundary layer and the Dirac comb to all orders

**(a) Yes — the series reproduces log(1+e^{−s}), and the exact profile is e^{τL(s)}−1.** On [0,π/4], with y = −log tanθ and s = y/τ (step_18, identity verified to ≤ 1.2e-51):
> **h_p(θ) − h_1(θ) = (1+e^{−2y})^{−1/2}·( e^{τ L(y/τ)} − 1 )**,  L(s) = log(1+e^{−s})  — an EXACT identity (entire in τ),
> so h_p−h_1 = (τ/√2)L(s) + O(τ²) in the layer, and its τⁿ-coefficient is the exact function (1+e^{−2y})^{−1/2}L(y/τ)ⁿτⁿ/n!.

Partial sums (a₀−A₀)/2 + Σ_{k≤K}(a_k−A_k)cos4kθ at K = 200, τ = 0.02 reproduce h_p−h_1 with errors 1.4e-6…5.6e-4, exactly the size of the 1/k² tail Σ_{k>200}(√2/π)/(4k²) (step_19), and normalized values track L(s) to tail accuracy.

**(b) The Dirac comb, precise to all orders (PROVED as τ→0 asymptotics; mechanism = Fourier inversion).** Distributional kernel:
> ½ + Σ_{k≥1}(−1)^k cos(2k gd x)cos4kθ = (π/4)Σ_{n∈ℤ}[δ(gd x − 2θ − (2n+1)π/2) + δ(gd x + 2θ − (2n+1)π/2)],
and on θ∈[0,π/4] the only hit is gd x = π/2−2θ ⇔ x = y; evaluating the s-integral against it returns the exact profile above. Summing the two sides of each edge normal (both contribute the same kernel K(x) = ½sech(x)(1+e^{−2x})^{−1/2}, with v = gd(x)/2 on both sides — the odd part of φ cancels exactly) gives the **all-orders comb**:
> **∫₀^{2π}(h_p−h_1)φ dθ = (π/4) Σ_{j=0}^{3} Σ_{m≥0} 2^{−2m} φ^{(2m)}(θ_j)/(2m)! · D(τ;m)**  (asymptotic in τ, θ_j = π/4 + jπ/2),
with the SAME D(τ;m) as in the k-expansion: the comb is a superposition of δ, δ'', δ'''', … at the four edge normals. Leading mass per normal = (π/4)D(τ;0) = (P(p)−P(1))/4 = (√2/4)(ζ(2)τ²+ζ(3)τ³+…) — F1's Corollary 3a is the m = 0, n ≤ 3 truncation.

**(c) Numerical verification with the asymmetric φ = e^{sin t + cos 3t/2}** (step_23): subtracting the comb through order N and dividing by τ^{N+1}, Richardson-extrapolated limits vs the coefficients predicted by the D̃-tables (independent route: no s-integrals): N=3: −3.85017 vs −3.85023 (diff 5.8e-5); N=4: −7.42640 vs −7.42652 (1.3e-4); N=5: 15.9416 vs 15.9422 (5.9e-4); N=6: 53.502 vs 53.527 (2.5e-2, consistent with the halving convergence 56.18→54.90→54.21→…→53.5). The N=3 limit (√2 ζ4/64)(−33S₀+7S₂) also fixes f1_06's unexplained −3.9 residual and **independently confirms the φ'' term, i.e. the k²-part of c₄**, and N = 4,5,6 confirm D̃_{n,1}, D̃_{n,2} for n ≤ 7.

---

## 3. Task 3: quadrature vs series (dps = 55, step_27–28)

a_k by direct θ-quadrature on [0,π/4] with corner-adapted breakpoints (independent of the s-integral machinery) vs A_k + Σ_{n=2}^{10} c_n(k)τⁿ:

| k | τ = 0.01 residual | τ = 0.005 residual |
|---|---|---|
| 0 | 1.255e-18 | 6.391e-22 |
| 1 | −2.702e-17 | −1.409e-20 |
| 2 | 3.728e-16 | 2.007e-19 |
| 3 | −2.82e-15 | −1.585e-18 |
| 4 | 1.427e-14 | 8.503e-18 |
| 5 | −5.328e-14 | −3.45e-17 |

Residuals scale exactly like the first omitted term c₁₁(k)τ¹¹ (magnitudes match D11-size coefficients and the k⁸-top of Q₁₁). Honest note on "≥30 digits": the computation precision is 55 digits throughout; since the τ-series is asymptotic, agreement at these fixed τ is 14–22 digits, precisely as the asymptotics dictates — reaching 1e-30 at τ = 0.01 would require terms through n ≈ 28. All 12 values above are printed by the code; digit counts of agreement: τ=0.01: 17.9, 16.6, 15.4, 14.6, 13.8, 13.3 digits; τ=0.005: 21.2, 19.9, 18.7, 17.8, 17.1, 16.5 digits.

---

## 4. Evidence table

| claim | status | residual / digits | script |
|---|---|---|---|
| master identity a_k−A_k = (−1)^kB(τ,2k), D̃ formula | PROVED | algebra + D11 check ≤ 1.8e-99 (n ≤ 9), 1.2e-97 (n=10) | step_12, 13, 27 |
| e_{r,m} exact rationals; g₀,g₁,g₂ match F1-verify | PROVED | exact (sympy) + 0.0 vs mpmath | step_12 |
| I_{m,j} alternating MZV of weight m+j+1 | PROVED | Stirling expansion; I_{0,1}, I_{1,1}, I_{0,2} exact | step_13 |
| Q_n polynomials, n ≤ 11 (all r), r=0 n ≤ 10 | PROVED structure + NUMERICAL values | PSLQ resid ≤ 1.2e-94 (n≤9), ≤ 3.3e-111 (n=11) | step_14, 29 |
| n = 12: r=0,2,3 (dps 160), r=1 (dps 240), r=5 exact | NUMERICAL | 2.5e-147, 0, 1.6e-148, 6.8e-226 | step_30, 31 |
| n = 12, r = 4 | NUMERICAL (dps 120 only — flag) | resid 0 at 120 digits | step_29 |
| degree conjecture deg = ⌊n/2⌋−1, all n; exact tops | PROVED | tops checked n ≤ 12, ≤ 9.4e-97 | step_16, 26 |
| A_{2m} classical identity (used in odd top) | PROVED (classical); here verified | ≤ 4.7e-61, m ≤ 6 | step_16 |
| top generating function (sinh + Re ψ); F1's 1/sinh = even part | PROVED modulo A_{2m}; NUMERICAL | ≤ 6.1e-52 | step_17 |
| exact profile identity e^{τL}−1 | PROVED | ≤ 1.2e-51 | step_18 |
| all-orders comb (δ, δ'', δ''''…) | PROVED (asymptotic) + NUMERICAL | extrap. limits match to 6e-5…2.5e-2 | step_21, 23 |
| quadrature vs series, task 3 | NUMERICAL | table above, 55 dps | step_28 |
| D̃_{3j,j−1} ⊃ (2j−2)!/j! ζ3^j | CONJECTURE (j ≤ 4: 1,1,4,30) | exact matches | step_14, 29, 30 |

## 5. Honest failures / caveats
- PSLQ for n = 12 at dps = 120 returned **spurious** relations for r = 0..3 (huge ugly coefficients, residuals ~1e-85); caught by implausibility and corrected at dps = 160/240 (clean rationals, residuals ≤ 1e-147). r = 4 was only obtained at dps 120 and is flagged for reconfirmation.
- The function 2τ Re Σ_{l≥1} C(τ,l)/(l−τ(½+iκ)) shares only the top-degree part with Q_n — it is NOT a full generating function (differs at O(τ⁴)); the correct top GF is the sinh/digamma form in §1.4(ii).
- The reduction of c_n(k) to zeta products (no log 2 / Li(½)) for all n, r remains CONJECTURAL beyond the proved cases, despite overwhelming PSLQ evidence; a conceptual proof is open.
- 30-digit quadrature-vs-series agreement at τ = 0.01/0.005 is impossible for any fixed truncation (asymptotic series); reported residuals match the first omitted term exactly as they should.

Key scripts (auto-saved in the folder): step_12 (exact e-table), step_13 (I-integrals, D̃, D11 validation), step_14 (PSLQ n ≤ 9), step_15 (Q-polynomials + c₇/c₄/c₅ checks), step_16 (A_{2m} + top GF), step_17 (top GF verification), step_18–19 (profile identity, Fourier resum), step_21/23 (all-orders comb), step_26 (n = 10), step_28 (task-3 quadrature), step_29–31 (n = 11, 12).