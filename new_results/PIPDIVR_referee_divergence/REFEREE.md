# PIPDIVR — referee report on `PIPDIV_divergence_proof/result.md`

Object: h = 1/p, κ = h/(1−h), π_p/4 = 2^{−h} + Ê(h), π_p/4 ~ 1 + Σ_{k≥0} N_k h^{k+1} as h → 0+.
I wrote all code myself. I read the author's scripts only to see which methods they used; I did not import or run any of them.

## Verdicts

| Item | Verdict |
|---|---|
| Theorem A (i) Laplace form | **PROVED** |
| Theorem A (ii) N_k = (−log2)^{k+1}/(k+1)! + T_k + U_k, bound on U_k | **PROVED** |
| Lemma 3 (local form at ±iπ) | **PROVED AFTER A MINOR REPAIR** (the domain for the identity theorem is not stated; R1 below) |
| Lemma 4 (Hankel bound), Proposition 5 | **PROVED** |
| Theorem B (all-orders large-k expansion) | **PROVED** (one parenthetical sentence is misleading; W2) |
| Theorem C, two-term law, both parities, sign period | **PROVED** (the even-k error terms are true but not sharp; W4) |
| Theorem C, three-term refinement | **PROVED** |
| Corollary D, p → ∞ (radius 0) | **PROVED** |
| Corollary D, p → 1 | **PROVED AFTER REPAIR R2**: it rests on the duality π_p = π_q. The cited source (D12 §1.2/§1.5) labels that NUMERICAL, so a proof must be cited. |

I found no mathematical error that changes a statement. Every step of the proof was checked by hand (section 2 below). Every
numerical claim I could test agrees with my independent computation to all printed digits.

---

## 1. Numerical tests

### 1.1 Independent N_k (k ≤ 150)
`rf1_Nk.py` evaluates the D12 moment expansion in the s-variable. G(y,h) = e^{−hy}(1+e^{−y})^{−1−h} = Σ g_j(h) y^j, and
(1+e^{−s})^h − 1 = Σ_r h^r L^r/r! with L = log(1+e^{−s}). This gives
Ê ~ Σ_j h^{j+2}(1−h)^{−j−1} g_j(h) Σ_r (h^r/r!) M_{j,r}, with M_{j,r} = ∫_0^∞ s^j L^r ds.
Everything is done in **arb ball arithmetic** (python-flint). Each M_{j,r} comes from the rigorous integrator `acb.integral` on
[0,S], plus the proven tail bound 2S^j e^{−rS}/r (valid for S ≥ 2j). So every N_k is a certified ball, given the moment
expansion, which D12 proved and which §1.3 below tests independently.
This route does not use Theorem A or the Laplace variable.

| run | k range | prec | S | max relative ball radius |
|---|---|---|---|---|
| `rf1_N_K81_p400_S300.txt` | 0–80 | 400 bits | 300 | — |
| `rf1_N_K101_p520_S420.txt` | 0–100 | 520 bits | 420 | 5.0e−52 |
| `rf1_N_K151_p700_S400.txt` | 0–150 | 700 bits | 400 | 7.4e−46 |

The runs agree pairwise to 3.0e−52 (k ≤ 100) and 1.2e−41 (k ≤ 80), inside the balls.
N_0..N_5 agree with the D12 closed forms inside the ball radius (about 1e−61).
N_2..N_30 agree with REF_PIP `r5_log.txt` to 4.6e−13 relative, which is all 13 of its printed digits.
Every ν_k, 2-term, 3-term, r₂, r₃ and L-column value in the author's table is reproduced to all printed digits.

### 1.2 Theorem C (`rf2_log_K151.txt`)
ν_k = N_k π^{k−1}/(2(k−2)!). ν₂ and ν₃ are the two- and three-term laws, which I implemented from the statement.
r₂ = (ν−ν₂)k²/log²k. r₃ = (ν−ν₃)k³/log³k.

| k | ν_k | two-term | three-term | r₂ | r₃ |
|---|---|---|---|---|---|
| 21 | 1.283311526 | 1.259726432 | 1.285540308 | 1.122 | −0.731 |
| 22 | 0.5086477974 | 0.3774638922 | 0.4733539929 | 6.645 | 12.72 |
| 39 | −1.123929765 | −1.133373032 | −1.126229627 | 1.070 | 2.775 |
| 40 | −0.2903389951 | −0.2527209846 | −0.2861291375 | −4.423 | −5.367 |
| 59 | −1.078163107 | −1.086575477 | −1.079153993 | 1.761 | 3.002 |
| 60 | −0.2059596364 | −0.1887275674 | −0.2049844765 | −3.701 | −3.069 |
| 61 | 1.075444187 | 1.083640715 | 1.076365003 | −1.805 | −3.009 |
| 62 | 0.2004248641 | 0.1842268605 | 0.1995614995 | 3.656 | 2.927 |
| 99 | −1.045892322 | −1.050874249 | −1.046194275 | 2.312 | 3.020 |
| 100 | −0.135279012 | −0.1286239372 | −0.1351507703 | −3.138 | −1.313 |
| 119 | −1.038208083 | −1.042177797 | −1.038402362 | 2.461 | 2.999 |
| 120 | −0.1165508475 | −0.1117908358 | −0.1164963293 | −2.991 | −0.859 |
| 139 | −1.032776981 | −1.036020454 | −1.032910163 | 2.574 | 2.977 |
| 140 | −0.1027566283 | −0.09916760249 | −0.1027336791 | −2.881 | −0.522 |
| 149 | 1.030614375 | 1.033570083 | 1.030726711 | −2.621 | −2.966 |
| 150 | 0.09712089391 | 0.09395730689 | 0.09710661392 | 2.835 | 0.383 |

For 21 ≤ k ≤ 150: max|r₂| is 2.62 (odd k) and 6.65 (even k); max|r₃| is 3.03 (odd k) and 12.7 (even k).
All scaled residuals stay bounded.

- **Odd k:** r₂ creeps up towards its predicted limit −σ_k π²/2 (4.93 in absolute value). The i = 2 term gives
  ν − ν₂ ≈ −σ_k π²(X² − 5π²/12 − 3)/(2k²), which is 2.44 at k = 149; the observed value is 2.62. So the odd-k error
  O(log²k/k²) is **sharp**.
- **Even k, bracket** B_k = τ_k N_k π^{k−2}/(2(k−3)!) against X_{k−2} = ψ(k−2) + γ − log π:

  | k | e = B_k − X_{k−2} | e·k/log²k (claimed bounded) | e·k/log k | e(k−3)/X_{k−3} |
  |---|---|---|---|---|
  | 20 | 0.9361 | 2.086 | 6.249 | 7.117 |
  | 60 | 0.3181 | 1.139 | 4.662 | 5.231 |
  | 100 | 0.2076 | 0.979 | 4.508 | 5.032 |
  | 150 | 0.1490 | 0.890 | 4.462 | 4.957 |

  The claimed O(log²k/k) holds. The true error is π²X_{k−3}/(2(k−3)) + … = O(log k/k): the last column tends to π²/2.
- **Even k, ν_k:** ν_k/(τ_kπ(log k+γ−log π)/k) − 1 times k is 9.57, 7.37, 6.65, 6.56 at k = 20, 40, 100, 140. So the
  relative error is O(1/k), not merely the claimed O(log k/k).
- **Signs** (`rf5_log.txt`): the pattern + + − − for k ≡ 1, 2, 3, 0 (mod 4) holds for k = 7 and for every 9 ≤ k ≤ 150. It
  fails at k = 3, 4, 5, 6, 8. For odd k, |N_k| ≥ (k−2)!/π^{k−1} holds for every odd k from 3 to 149.
- **Growth:** |N_k|^{1/k}/(k/(eπ)) = 0.893, 0.914, 0.940, 0.969 at k = 20, 50, 100, 149, tending to 1 slowly, as Stirling predicts.

### 1.3 The N_k really are the asymptotic coefficients (`rf3_log_K151.txt`)
I computed π_p/4 directly from the definition π_p = (2/p)∫_0^1(t^{1−p}+(1−t)^{1−p})^{1/p}dt, substituting t = v^p.
This direct value agrees with the D12 split to ≤ 8e−62. I then compared it with the optimally truncated series:

| p | best K | min_K |π_p/4 − S_K| | e^{−πp} | ratio |
|---|---|---|---|---|
| 6 | 14 | 3.06e−11 | 6.51e−9 | 4.7e−3 |
| 10 | 32 | 7.86e−17 | 2.27e−14 | 3.5e−3 |
| 20 | 64 | 3.53e−31 | 5.16e−28 | 6.8e−4 |
| 30 | 96 | 3.03e−45 | 1.17e−41 | 2.6e−4 |
| 40 | 126 | 3.43e−59 | 2.66e−55 | 1.3e−4 |

This end-to-end check is independent of both the author's formula and mine. The ratio falls roughly like p^{−2.3}, so the
optimal-truncation error is e^{−πp} times a negative power of p (W5).

### 1.4 Theorem A(ii): the split T_k + U_k (`rf2_log_K151.txt`, Table 1)
I computed T_k exactly from my own bivariate series of A = e^{2y}(1+e^y)^{−1−h}, and U_k = N_k − T_k − (−log2)^{k+1}/(k+1)!.
Since N_k comes from a different formula (§1.1), this is a genuine test of A(ii).

| k | T_k/N_k − 1 | U_k(2π)^{k−1}/(k−2)! | U_k/(k(k−2)!6^{−k} + 2^k) |
|---|---|---|---|
| 20 | −1.97e−6 | −1.162 | −4.9e−6 |
| 60 | −1.38e−18 | −0.3266 | −2.2e−3 |
| 100 | −1.25e−30 | −0.2149 | −1.3e−4 |
| 140 | −1.15e−42 | −0.1645 | −1.2e−5 |

At k = 150, U_150 is below the ball radius of N_150 (about 7e−46 relative), so that row carries no information.
The proven bound on U_k holds with a lot of room. The sharper behaviour |U_k| ≈ c(k−2)!/(2π)^{k−1}, with c slowly
decreasing, matches the n = 2 singularities at ±2πi.

### 1.5 Theorem B (Table 4 of `rf2_log_K151.txt`)
I built T_k^{(L)} from scratch. φ_l(h) come from my series of Φ. The factors (1+h−l)_m are expanded as exact polynomials
in h, and e^{−h Log y₀} as a series.

Column L = 1..4 gives |T_k − T_k^{(L)}|/(2(k−2)!/π^{k−1}); parentheses give that value times k^{L+1}/log k, which Theorem B says is bounded.

| k | L = 1 | L = 2 | L = 3 | L = 4 |
|---|---|---|---|---|
| 21 | 5.39e−2 (7.81) | 2.31e−3 (7.03) | 1.62e−3 (103) | 4.73e−5 (63.4) |
| 22 | 2.46e−2 (3.85) | 6.81e−3 (23.5) | 1.15e−4 (8.69) | 2.55e−4 (425) |
| 81 | 2.53e−3 (3.78) | 9.22e−6 (1.12) | 2.89e−6 (28.3) | 5.74e−9 (4.56) |
| 82 | 4.19e−4 (0.639) | 7.23e−5 (9.05) | 2.32e−7 (2.38) | 1.13e−7 (94.8) |
| 141 | 7.96e−4 (3.20) | 1.08e−6 (0.613) | 2.83e−7 (22.6) | 2.63e−10 (2.96) |
| 142 | 8.51e−5 (0.346) | 1.29e−5 (7.47) | 1.67e−8 (1.37) | 6.29e−9 (73.3) |

Every scaled column is bounded and decreasing.

### 1.6 Lemmas (`rf4_log.txt`)
- **(a) Theorem A(i).** I used the closed form I_n = ₂F₁(1+h, a; a+1; −1)/a with a = n(1−h)/h + h, and summed over n
  with nsum. The result matches the s-integral to 2e−43 or better at p = 3, 4, 10.
- **(b) Lemma 3.** I continued A analytically from y = 0 by integrating e^w/(1+e^w) along [0,y], rather than using a
  principal log. Principal Log(1+e^y) jumps on Im y = π, Re y > 0, which lies inside V_+.
  - Test points: 8 points in V_+, including both sides of the slit and both sides of the line Im y = π.
  - Values of h: 0.2+0.15i, −0.25, 0.1−0.3i, 0.7.
  - Result: max |A/RHS − 1| = 7.9e−41.
- **(c) Lemma 4.** I took h = ¼e^{iθ} for θ = 0, π/4, π/2, π, 3π/2 and m ≤ 400 (1200-bit series).
  - Quantity tested: sup over windows of m of |[y^m]R_L|π^m m^{L−|h|}.
  - L = 1, 2, 3 all stay bounded. For example, θ = 0, L = 1 gives 0.529, 0.528, 0.527 on the windows 101–200, 201–300, 301–400.
  - For Re h < 0 these values decay. That is expected, since the true exponent is Re h and the lemma uses |h| ≥ Re h.
- **(d) Gamma bound.** max |Γ(M+h)/Γ(M)|/M^{1/4} = 1.2254 on |h| = 1/4, against the claimed 1.4.

---

## 2. Proof, step by step

**Lemma 1 / Theorem A(i).** Correct. For 0 < h < 1 the coefficients C(h,n), n ≥ 1, alternate as (−1)^{n−1}, and
Σ|C(h,n)| = 1. Domination by hκe^{−hκs} justifies Fubini. The substitution y = κs gives
e^{−ns} = e^{−ny/h}e^{ny} (because 1/κ = 1/h − 1), and e^{−hy}(1+e^{−y})^{−1−h}e^{ny} = e^{(n+1)y}(1+e^y)^{−1−h}.
Each Laplace integral converges because n − h − n/h < 0.

**Lemma 2 / Theorem A(ii).** Correct.
- **Cauchy radii.** Radius 3 < π for n ≥ 1, using |t/n| ≤ 3 and e^{(n+1)|t|/n} ≤ e^6. Radius 6 for n ≥ 2, giving e^9.
  In both cases the branch ℓ is defined on |w| < π.
- **Binomial bound.** |C(h,n)| ≤ |h|e^{|h|}n^{|h|−1}, using H_{n−1} ≤ 1 + log n.
- **Watson remainder, uniform in n.** On [0,2] the geometric tail is at most 3M₀(t/3)^M. On [2,∞),
  H_n(t/n,h) ≤ e^t and 1/h − 1 ≥ 1/(2h) for h ≤ 1/2, so the tail is at most 2he^{−1/h}. The polynomial tails are at
  most e^{−1/h}m!(2h)^{m+1}.
- **The D_m.** Each D_m is analytic on |h| ≤ 1/2 with D_m(0) = 0 and |D_m| ≤ C3^{−m}|h|.
- **Uniqueness of asymptotic coefficients.** Used correctly. Truncating each D_m at order R and M at the target order gives
  an O(h^{J+1}) expansion, and asymptotic coefficients are unique. Index bookkeeping for (2.3): m + 2 + r = k + 1 with r ≥ 1.
- **The n = 1 part.** C(h,1) = h and b_{1,m} = [y^m]A. This gives T_k with [h^{k−2−m}].
- **U_k.** Cauchy on |h| = 1/2 gives the 2^r factor. The split of Σ m!12^{−m} is right: terms ≤ 1 for m ≤ 12, and
  increasing from m = 11 on. Since 2^{k−1}·144·12^{−k} = 72·6^{−k}, the stated bound follows.
- **The log 2 term.** Its sign and index, (−log2)^{k+1}/(k+1)!, are right.

**Lemma 3.** The algebra is correct.
- Since e^{y₀} = −1, we have 1 + e^y = (y₀ − y)φ(u) and e^{2y} = e^{2u}.
- Base-point check: log π + iπ/2 + log(2/π) − iπ/2 = log 2. Also arg φ(−iθ') = −θ'/2.
- **Gap (R1).** The proof compares the two logarithms "near y = 0", but y = 0 ∉ V_±, and no domain is named that contains
  both 0 and V_± on which both expressions are analytic.
- **Repair (R1).** Use W := D ∩ {|y − y₀| < 2π} ∖ σ_±. On W, ℓ, Log y₀ + Log(1 − y/y₀) and log φ(y − y₀) are all
  analytic, and W contains 0 and V_±. W is connected: a convex set minus one segment running to its boundary, since
  σ_∓ lies outside the disc. Two continuous logarithms of 1 + e^y on W differ by a constant in 2πiZ, and the constant is
  0 at y = 0.
- The numerics (1.6b) confirm the branch bookkeeping, including on both sides of the slit.

**(2.4).** Correct. For integer l the definition of (y₀−y)^a is additive in the exponent. Expanding
[x^m](1−x)^{−β} = (β)_m/m! gives the stated m![y^m]S_l^±. The −iπ coefficients are the complex conjugates of the +iπ ones
for real h, because φ_l has real coefficients. That justifies "2 Re".

**Lemma 4.** Correct.
- **h = 0.** The residues are e^{2y₀}/e^{y₀} = −1 and −1. Both cancel, and the next poles are at ±3iπ.
- **General h, local bound.** Cauchy on |u| = π (Φ is analytic for |u| < 2π) gives the tail bound 2M_Φ(|u|/π)^L on
  |u| ≤ π/2. The power satisfies |(y₀−y)^{−1−h}| ≤ C|u|^{−1−Re h} with a bounded e^{π|Im h|} factor, so α ∈ [L − 5/4, L − 3/4].
- **Contour.** It lies in D, and R_L has continuous boundary values on the slits.
- **The pieces of the contour integral:**
  - Outer arcs: (5π/4)^{−m}.
  - g_+, moved across the slit: (5π/4)^{−m}.
  - Small circle: (1 − 1/m)^{−m−1} ≤ e^{5/3} for m ≥ 4.
  - Slit sides: the Beta integral Γ(α+1)Γ(m−α)/Γ(m+1) ~ m^{−α−1}, uniform for α in a compact subset of (−1, ∞).
- Small m is handled by Cauchy on |y| = 1. Joint analyticity of R_L in (y, h) gives analyticity of [y^m]R_L in h.
- Verified numerically in 1.6c.

**Proposition 5 / Theorem B.** Correct.
- The i = 0 term uses the h = 0 bound, giving (1/2)^k relative.
- For i ≥ 1, Cauchy on |h| = 1/log k (≤ 1/4 for k ≥ e⁴) uses m^{1/log k} ≤ e.
- The ratio test is fine: (m/(m−1))^L ≤ 2^L for m ≥ 2, with threshold m₁ = 2^{L+1}π log k. This gives the geometric sum 2t₁.
- t₁ is (log k)π²(k−3)^{−L}/(k−2) relative, which is O(k^{−L−1}log k).
- The m < m₁ block: m₁·m₁!(log k)^k = exp(O(k log log k)) = o((k−3)!π^{−k}k^{−L−1}).
- The log 2 and U_k terms are O(k(π/6)^k) relative.

**Theorem C.** Correct.
- **Setup.** With L = 1: (1+h)_m = Γ(m+1+h)/Γ(1+h), M = m + 1 = k − 1 − i, and y₀^{−1−m} = y₀^{−(k−1−i)}.
- **Explicit terms.** [h⁰]F_M = Γ(M). [h¹]F_M = Γ(M)(ψ(M) + γ − Log y₀), using (log Γ(1+h))' = −γ at h = 0.
- **Gamma bound.** Checked case by case: x > 0 modulus, log-convexity both ways, and 1 ≤ M < 2 with Γ(3/4)/Γ_min = 1.384.
- **Constants.** |Log iπ| = 1.944, so |e^{−hLog y₀}| ≤ e^{0.49} on |h| ≤ 1/4.
- **Error.** The i ≥ 2 tail is O(log²k/k²) relative; adding Theorem B with L = 1 gives the stated two-term law.

**Parities.** I re-derived these.
- Odd k: i^{1−k} = σ_k and i^{2−k} = iσ_k. Re[iσ(X − iπ/2)] = σπ/2, which gives π²/(2(k−2)).
- Even k: Re i^{1−k} = 0 and i^{2−k} = τ_k.
- The even-k error is genuinely smaller than the main term. Main term X ≍ log k; error O(log²k/k), since
  (k−2)/(2π)·log²k/k² = O(log²k/k). So the main term dominates for both parities.
- Signs: + + − − is right for large k (numerically from k = 9 on).

**Three-term refinement.** I re-derived every coefficient.
- l = 0, i = 2: [h²]F_M = Γ(M)·½[Y² + ψ'(M) − ζ(2)].
- l = 1, i = 1: −(3/2)(k−4)!y₀^{3−k}. Its real part vanishes for even k and gives +3σπ²/(2(k−2)(k−3)) for odd k.
- Together these give exactly −π²(X² − 5π²/12 − 3 + ψ')/(2(k−2)(k−3)) for odd k, and π³X_{k−3}/(2(k−2)(k−3)) for even k.
- The error bounds hold (Cauchy with |(h)_m| ≤ C|h|Γ(m)). For odd k they are not sharp: the true next term is O(log²k/k³).

**Corollary D.**
- **p → ∞.** Correct. Stirling gives |N_k|^{1/k} ~ k/(eπ) for both parities.
- **p → 1.** The algebra is correct: h = 1 − 1/q = −η, so π_q/4 = 1 + Σ_k (−1)^k N_{k−1} η^k, an asymptotic expansion as
  η → 0−. But the step uses π_p = π_q, which D12 §1.2 and §1.5 label NUMERICAL.
- **Repair (R2).** Cite a proof of the duality: Keller and Vakil, "π_p, the value of π in ℓ_p", Amer. Math. Monthly 116 (2009),
  where π_p = π_q is proved. The earlier referee (REF_PIP) treated it as a literature fact. Optionally add: radius 0 is
  invariant under any analytic change of local variable (for example η ↦ q − 1 = −η/(1+η)), so the series in powers of
  q − 1 also has radius 0.

---

## 3. Sentences that claim more than is established (W-list), with suggested wording

1. **W1, §4 "Statement for the paper": "…with signs of period four."** False for small k: N_3, N_4, N_5, N_6 and N_8 break the
   pattern. Suggested: "…with signs following the period-four pattern + + − − (k ≡ 1, 2, 3, 0 mod 4) for all sufficiently
   large k (numerically for every k ≥ 9)." Also "smaller by a factor π(log k + γ − log π)/k" →
   "smaller by a factor asymptotic to π(log k + γ − log π)/k".
2. **W2, Theorem B parenthetical: "Each T_k^{(L)} is an explicit finite sum of rational combinations of log π, π, γ and values
   ψ^{(j)}(integer)."** This is misleading. The coefficients of (1+h−l)_m and φ_l(h) are rational, and γ and the ψ^{(j)}
   appear only after rewriting harmonic-type sums (ψ(n) + γ = H_{n−1} ∈ Q). Suggested: "Each T_k^{(L)} is an explicit
   finite sum. It is a polynomial in log π and π^{±1} with rational coefficients."
3. **W3, Corollary D, p → 1, "(duality, D12 §1.5)".** That source has the duality as NUMERICAL. Replace with a reference
   to a proof (R2).
4. **W4, Theorem C even-k error terms.** These are correct but not sharp. That is not an overclaim. The author's own
   three-term refinement already implies the bracket error is O(log k/k), and that ν_k = τ_kπ(log k+γ−log π)/k·(1 + O(1/k)).
   Suggested: either keep the statement as is, or state the sharper forms with "(by the three-term refinement)".
   Numerics: e·(k−3)/X_{k−3} → π²/2, and k·(relative error) ≈ 6.6.
5. **W5, Remark (Borel picture): "Optimal truncation error is then of size e^{−π/h} = e^{−πp}".** It is labelled
   heuristic, which is fine. Numerically the error is e^{−πp} times a decreasing power of p (ratio 4.7e−3 at p = 6,
   1.3e−4 at p = 40). Suggested: "of size e^{−πp} up to a power of p".
6. **W6, preamble: "All constants hidden in O(·) below exist but are not computed".** This is accurate, but the sign-period
   and lower-bound statements in Theorem C and Corollary D are therefore ineffective ("for large k", with no explicit
   k₀). That is harmless for radius 0. Do not let the paper state an explicit range unless it says the range was checked
   numerically, for example: "the pattern holds for 9 ≤ k ≤ 150 by direct computation".
7. **W7, Lemma 3 proof: "Near y = 0 both are exp(…)".** The point y = 0 is not in V_±. Use the repair R1 wording:
   "on the connected open set W = D ∩ {|y − y₀| < 2π} ∖ σ_±, which contains 0 and V_±, both exponents are continuous
   logarithms of 1 + e^y, so they differ by a constant in 2πiZ; at y = 0 the constant is 0."
8. **W8, §2.2 parenthetical: "Numerically U_k ≈ ±(k−2)!/(2π)^{k−1}".** The constant is not 1. It decreases slowly: 1.16, 0.33,
   0.21, 0.16 at k = 20, 60, 100, 140. Suggested: "U_k is of order (k−2)!/(2π)^{k−1} (n = 2 singularities at ±2πi)".

## 4. Repairs

- **R1 (Lemma 3).** Name the domain W as in W7. This is a one-sentence fix and does not change any statement.
- **R2 (Corollary D, p → 1).** Cite a proof of π_p = π_q, for example Keller–Vakil (2009). Without it, the p → 1 half is
  conditional on a numerically verified identity.

## 5. What I consider proved (exact statement)

With h = 1/p, π_p/4 ~ 1 + Σ_{k≥0} N_k h^{k+1} as p → ∞, y₀ = iπ, Log principal, X_j = ψ(j) + γ − log π:

1. N_k = (−log2)^{k+1}/(k+1)! + T_k + U_k, where T_k = Σ_{m=0}^{k−2} m![y^m h^{k−2−m}] e^{2y}(1+e^y)^{−1−h} and
   |U_k| ≤ C(2^k + k(k−2)!6^{−k}).
2. For every L ≥ 1: N_k = T_k^{(L)} + O_L((k−2)!π^{−(k−1)}k^{−L−1}log k), with T_k^{(L)} as in Theorem B.
3. N_k = 2Re[(k−2)! y₀^{1−k} + (k−3)! y₀^{2−k}(ψ(k−2) + γ − Log y₀)] + O((k−2)!π^{−(k−1)} log²k/k²). Equivalently:
   - odd k: N_k = (−1)^{(k−1)/2}·(2(k−2)!/π^{k−1})·[1 + π²/(2(k−2)) + O(log²k/k²)], and this error order is attained;
   - even k: N_k = (−1)^{(k−2)/2}·(2(k−3)!/π^{k−2})·[X_{k−2} + O(log²k/k)]. The three-term refinement improves this to
     O(log k/k).
   - The three-term refinements of Theorem C hold as stated.
   - For all sufficiently large k (no explicit threshold), sign N_k is + + − − for k ≡ 1, 2, 3, 0 (mod 4).
4. |N_k|^{1/k} ~ k/(eπ), so Σ N_k h^{k+1} has radius of convergence 0.
5. Given the duality π_p = π_q (which must be cited, R2), the p → 1 series Σ(−1)^k N_{k−1}η^k, with η = 1/q − 1, also has radius 0.

## 6. What I did not check

- I did not re-check D12's moment expansion (the formula behind my N_k) beyond the end-to-end test of §1.3; the brief said not to re-referee it.
- I did not re-check the D12 closed forms of N_0..N_5; I only compared numbers.
- I did not re-check the literature proof of the duality π_p = π_q.
- My N_k balls are rigorous only given the moment expansion. Nothing beyond k = 150 was computed.
- I did not check the existence or size of the uncomputed O-constants. In particular I did not determine an explicit k₀
  for the sign pattern.
- I did not check the author's numerical files `pd*` or their quadrature settings, except by reproducing every number in
  the result.md table.
- I did not try to prove the sharper error orders noted in W4 and W5 beyond what follows from the author's own three-term
  refinement.

## 7. Files (all in this folder)

- `rf_series.py`: truncated power-series helpers (arb/acb).
- `rf1_Nk.py K PREC S OUT`: certified N_k from the moment expansion. Outputs `rf1_N_K81_p400_S300.txt`,
  `rf1_N_K101_p520_S420.txt`, `rf1_N_K151_p700_S400.txt`; logs `rf1_log_K*.txt`.
- `rf2_laws.py NFILE KMAX`: closed forms, r5 comparison, T_k/U_k, Theorem C laws and residuals, Theorem B approximants.
  Logs `rf2_log_K101.txt`, `rf2_log_K151.txt`.
- `rf3_direct.py NFILE plist dps`: direct π_p from the definition vs optimal truncation. Logs `rf3_log_K101.txt`, `rf3_log_K151.txt`.
- `rf4_lemmas.py`: Theorem A(i), Lemma 3 at complex h, Lemma 4 at complex h, Gamma bound. Log `rf4_log.txt`.
- `rf5_signs.py NFILE`: sign pattern, odd-k lower bound, growth rate. Log `rf5_log.txt`.

Run with `../../../.venv/bin/python`. No processes were left running.
