# RSTSR: referee report on RSTS (sharpness of the rest exponent √(2π))

Refereed: `new_results/RSTS_rest_sharpness/result.md` (+ scripts). Inputs treated as refereed: CX6 (8), (10)–(18).
All computations in this folder use the referee's own code (Arb integrals, own mpmath integrators, own sympy derivation).
N = k+1, θ₂ = √(2πN), λ = √(N/π), Λ = √(πN), ω = e^{iπ/4}, a = π/(2√2).

## Verdicts

| Claim | Verdict |
|---|---|
| Decomposition rest_k = S_k + U_k (§2.1), incl. A^{[2]} closed form and contour move | CONFIRMED |
| (B3) bound on A^{[≥3]} | CONFIRMED |
| (ii) \|U_k\| ≤ 1.0584·10⁷ e^{−(3/2)√(π(N−3))}, N ≥ 100 | CONFIRMED (C_U = 10 583 366.83 reproduced with Arb) |
| (i) S_k = (√π/N)e^{−θ₂}[cos(θ₂−π/4) + O(N^{−1/2})] | CONFIRMED (O-constant not effective, as stated) |
| Corollary: limsup = √π, liminf = 0; exponent √(2π) and power N^{−1} sharp | CONFIRMED |
| ϖ₂ closed form, status "DERIVED" | CONFIRMED as derived; status label correct (not proved) |
| RST1 phase "cos√(2πN)" wrong, correct phase −π/4 | CONFIRMED |
| §0: "+O(e^{−√(2πN)})" in eq:Mkfull only meaningful for the exact M_k^main | CONFIRMED |
| §6 proposed paper text | CONFIRMED AFTER REPAIR (one sentence misstates what is bounded; repair below) |

## 1. Decomposition (§2.1)

**A^{[2]}.** Re-derived by hand. In CX6 (10) expand (sL−β)^k binomially. i=0: (sL)^k(1/(1+w)−1) = (sL)^k(−w + w²/(1+w)). i=1: −k(sL)^{k−1}β/(1+w) = −k(sL)^{k−1}(w + r₁). The "−w" parts give
(1/k!)∫₀^∞ e^{−(s+2)L}[−(sL)^k − k(sL)^{k−1}]dL = −s^k/(s+2)^{k+1} − s^{k−1}/(s+2)^k = −2(s+1)s^{k−1}/(s+2)^{k+1}. Correct.
Checks: sympy exact integral at k=7 (difference 0); numerically A_k − A^{[≥3]} − A^{[2]} = O(10^{−60}) relative at 6 points s on Γ, k = 20, 40, with A_k from the full CX6 kernel and A^{[≥3]} from the remainder kernel (`r4_A2_and_U_direct.log`).

**Contour move.** HEA^{[2]} is O(s^{k−2}) at 0, analytic in Re s > 0, pole −2 outside; the strip between Γ and the π/4-ray meets ℝ in [0,½]; connecting segments at infinity vanish (|E| = e^{−π Im s/2}). Orientation gives ∫_Γ = ∫_ray − ∫₀^{1/2}, consistent with CX6's own −∫₀^{1/2}HEQ_k. Ray vs positive real axis: arc contribution O(R^{−3/2}) by Stirling for H and |A^{[2]}| ≤ 2/|s|. Correct.

**Direct numerical test of the whole identity** at k = 100: U_k computed from its definition, −J_k + (2/π)Im{I₁+I₂+I₃−∫₀^{1/2}HEQ_k−∫₀^{1/2}HEA^{[2]}}, with I₂ taken from the remainder kernel (no closed form), gives −4.97511773·10^{−17}. Independently, M_k(data) − M_main(Arb) − S_k(Arb) = −4.9751177·10^{−17}. Agreement to all printed digits.

## 2. (B3)

|β−w| ≤ |w|² from β−w = −w²∫₀¹t/(1+tw)dt and |1+tw| ≥ ½ (CX6 (11)); r₁ = (β−w−w²)/(1+w), so |r₁| ≤ 4|w|²; |w²/(1+w)| ≤ 2|w|². Extracting e^{−2L} analytically and rotating to arg L = −arg(s+3) ∈ (−π/4, 0] is the CX6 §5 argument verbatim (same wedge, no zeros of 1+w, Re((s+3)L) > 0 during rotation). Factorials cancel: 2ρ^k/|s+3|^{k+1} + 4ρ^{k−1}/|s+3|^k, and the i ≥ 2 terms are those of CX6 (12). Correct. Numerically |A^{[≥3]}|/B3 ∈ [0.002, 0.27] at the test points.

## 3. U_k bound (ii)

Each term re-derived:
- I₁: |H||K−E||A^{[2]}| ≤ 2·2e^{−3ar}·2ρ^{k−1}/|s+2|^k (using |s+1| ≤ |s+2|) → 8F(k−1,2;3a); exponent √(3·3a·2(k−1)/√2) = √(9π(N−2)/2). Correct.
- I₂: |HK| ≤ 4e^{−ar}; terms 2F(k,3;a), 4F(k−1,3;a) (exponents √(9π(N−1)/4), √(9π(N−2)/4)); i-terms F(k−i,i+1;a) with prefactor (4/a)e^{a/2}e^{π(i+1)/4} since a/√2 = π/4, so 2^i e^{πi/4} = x^i. Identity (k−i)(i+1) − 3(k−2) = (i−2)(k−3−i) verified symbolically; it is ≥ 0 for 2 ≤ i ≤ k−3, which contains 2 ≤ i ≤ k/2 once k ≥ 6. i > k/2: CX6 (16). Correct.
- I₃: CX6 (17) with √(9πk/4) ≥ (3/2)√(π(N−3)). Correct.
- Elementary: |J_k| ≤ e^{−(N−1)/2}; segment bounds (2/(πk))3^{−k} (CX6 (18)) and (6/π)4^{−k}/(k−1) (re-derived: 2(s+1)/(s+2) ≤ 3/2 on [0,½], s^{k−2}/(s+2)^k ≤ s^{k−2}2^{−k}). Both ≤ e^{−(N−1)/2}, total coefficient 2. Correct.
- Range: (N−1)/2 ≥ (3/2)√(π(N−3)) for N ≥ 100 (49.5 vs 26.18; derivative positive once √(N−3) > 1.5√π, i.e. N ≥ 11). Also checked every exponent comparison on N ∈ [100, 5000] ∪ {10⁴,…,10¹²}: no failures.

Arb (`r1_CU_arb.log`): I₁ 3601.126, I₂ (i=0,1) 1011.623, I₂ (2 ≤ i ≤ k/2) 5633.368, I₂ (i > k/2) 10 572 948.05, I₃ 170.658, elementary 2; **C_U = [10583366.8267599533 ± 2.7e−11]**, matches RSTS. Actual |U_k|/bound ≤ 1.8·10^{−12} for 100 ≤ k ≤ 400 (`r2_rest_table.log`): the bound is valid and extremely loose, as RSTS says.

## 4. Saddle (i)

Recomputed independently (by hand, then sympy in `r5_varpi2.py`):
- A^{[2]} = −(2/s)(1+x/2)(1+x)^{−N}, x = 2/s; Nx = Λω̄/r; Nx²/2 = −iπ/(2r²). Correct.
- H = √(π/2)s^{−3/2}(1 − 3/(4s) + …) (Γ(z)/Γ(z+3/2) = z^{−3/2}(1 − 3/(8z) + …), z = s/2); uniform in |arg s| ≤ π/4 by DLMF 5.11.13 plus compactness for 1 ≤ |s| ≤ R₀. Correct.
- ψ(r) = iωr − ω̄/r, Re ψ = −(r+1/r)/√2; ψ(1) = −√2 + i√2; ψ′(1) = 0; ψ″(1) = −2ω̄ (sympy: f″(1) = −ψ″(1) = 2e^{−iπ/4}). p(1) = −(√π/2)e^{−3iπ/8}e^{−iπ/2} with principal s^{−3/2}. √(2π/(−Λψ″(1))) = √(π/Λ)e^{iπ/8} (principal root, Re(−ψ″) = √2 > 0). Product −(√π/2)e^{−3iπ/4}; λ^{−3/2}Λ^{−1/2} = √π/N. Hence Z_k ~ −(π^{3/2}/(2N))e^{−θ₂}e^{i(θ₂−3π/4)} and (2/π)Im Z_k ~ −(√π/N)e^{−θ₂}sin(θ₂−3π/4) = (√π/N)e^{−θ₂}cos(θ₂−π/4). **Amplitude √π and phase −π/4 confirmed.**
- (a): |x| ≤ 4/λ ≤ ½, N|τ| ≤ (2/3)π/(λr³) ≤ 134/λ. Note: at λ = 8 this is 16.8, so "|ε_N| ≤ C/λ" holds with a huge C (≈ 134e^{16.8}); harmless for an O-statement, but it confirms the O-constant is not usable numerically. Perturbation O(λ^{−1}·λ^{−3/2}Λ^{−1/2}e^{−θ₂}) = O(N^{−3/2}e^{−θ₂}). Correct.
- (b) Laplace lemma: window |r−1| ≤ Λ^{−2/5} (then Λ|t|³ ≤ Λ^{−1/5}), tails ≤ e^{−√2Λ−Λ^{1/5}/(2√2)} and, for r ≥ 2, ≤ e^{−√2Λ−Λ(r−1)/(2√2)}; odd terms vanish on the symmetric window; remainder O(Λ^{−3/2}). Standard and correct.
- (c): |1+2/s|² = 1 + 2√2/|s| + 4/|s|² on the ray; gives |s/(s+2)|^N ≤ e^{−2√2Λ+8π}; |s| ≤ 1 piece 4√2/(2^{k+1}(k−1)). Both re-derived; correct and o(N^{−3/2}e^{−√2Λ}).

Numerics: Arb enclosures of S_k agree with RSTS at k = 100, 400 (ρ_k identical to 13 digits). mpmath log-form integrator, k up to 6 553 600 (`r3b_Sk_large_mp.log`): √N·[N e^{θ₂}S_k/√π − cos(θ₂−π/4)] stays in [−2.48, 2.16], consistent with O(N^{−1/2}) and with |ϖ₂| = 2.549.

## 5. Corollary

N e^{θ₂}|U_k| ≤ C_U N e^{θ₂−(3/2)√(π(N−3))} → 0 since (3/2)√π − √(2π) = 0.152 > 0. So N e^{θ₂}rest_k = √π cos(θ₂−π/4) + o(1). θ₂(N) is increasing, unbounded, with gaps √(2π)(√(N+1)−√N) ≤ √(2π)/(2√N) → 0. For every large target t (a maximum of |cos(·−π/4)|, or a zero) the first N with θ₂(N) ≥ t is within one gap of t. Along these N, |cos| → 1 or |cos| ≤ gap → 0. **No equidistribution is needed**; monotonicity plus vanishing gaps suffice, and the O(N^{−1/2}) correction vanishes in both limits. limsup = √π, liminf = 0: correct. Sharpness of the power: limsup > 0 rules out C N^b e^{−√(2πN)} for b < −1, and finiteness of limsup gives the bound with b = −1. Correct.

## 6. ϖ₂

Independent derivation (`r5_varpi2.py`): classical derivative-form second-order Laplace formula in the real variable r (not RSTS's Gaussian-moment code), plus the O(1/λ) amplitude terms −3/(8rω) (Stirling), +1/(2rω) (1+x/2), −π/(3r³ω³) (e^{−Nx³/3}). Result equals RSTS's ϖ₂ to 160 digits. The same code reproduces the paper's ϖ₁ (difference 0 at 170 digits). Numerical: Richardson fit of (ρ_k−1)√N over k = 1.0·10⁵ … 6.55·10⁶ gives ϖ₂ ≈ 2.2508949528675 − 1.1957189247196i, i.e. agreement to 11 digits; next coefficient ≈ 3.193 − 2.731i (bounded).
Status: "DERIVED" is the right label. A proof would need the Laplace expansion to the next order with a remainder bound and ε_N to O(λ^{−2}); this is routine but not written. Do not call it proved.

## 7. §0 and the paper text (§6 of result.md)

**§0 is right.** In eq:Mkfull the first line, with O(N^{−1}) inside the bracket, carries an error of size N^{−2}e^{−√(πN)}, which is much larger than e^{−√(2πN)}. The same holds for any finite truncation. So "+O(e^{−√(2πN)})" only makes sense if the first line *means* the exact M_k^main. The paper already says this in words (line 1273 of paper/ieee/main.tex: "the part that gives the first line"). RSTS's added sentence is correct and useful.

**One sentence in the proposed text is not what is proved.** "The pole s = −2 gives the term [display of M_k−M_k^main] … and everything else is at most 1.06·10⁷e^{−(3/2)√(π(N−3))}". Read literally, "everything else" is M_k − M_k^main minus the displayed leading term (√π/N)e^{−θ₂}cos(θ₂−π/4). That difference is of order N^{−3/2}e^{−√(2πN)}, which is *not* ≤ 1.06·10⁷e^{−(3/2)√(π(N−3))} for large N. The bound holds only for U_k = rest − S_k, where S_k is the **exact** integral of the pole s = −2. The displayed sentence never defines S_k. (1.06·10⁷ ≥ C_U is fine.)

**Repair (replacement text):**

> Let $M_k^{\rm main}=\frac2\pi\operatorname{Im}\int_0^\infty H(s)e^{i\pi s/2}s^{-1}(\frac{s}{1+s})^N ds$ be the exact contribution of $\pm i\pi/2$, and let $S_k=\frac2\pi\operatorname{Im}\int_0^\infty H(s)e^{i\pi s/2}A_k^{[2]}(s)\,ds$, $A_k^{[2]}(s)=-2(s+1)s^{k-1}/(s+2)^{k+1}$, be the exact contribution of the pole $s=-2$. Then $|M_k-M_k^{\rm main}-S_k|\le1.06\cdot10^{7}e^{-\frac32\sqrt{\pi(N-3)}}$ for $N\ge100$, and a saddle point, as for the main term, gives $S_k=\frac{\sqrt\pi}{N}e^{-\sqrt{2\pi N}}[\cos(\sqrt{2\pi N}-\frac\pi4)+O(N^{-1/2})]$. Hence
> $$M_k-M_k^{\rm main}=\frac{\sqrt\pi}{N}e^{-\sqrt{2\pi N}}\Bigl[\cos\bigl(\sqrt{2\pi N}-\tfrac\pi4\bigr)+O(N^{-1/2})\Bigr],$$
> and $\limsup N e^{\sqrt{2\pi N}}|M_k-M_k^{\rm main}|=\sqrt\pi$: the exponent $\sqrt{2\pi}$ and the power $N^{-1}$ are sharp. The constant in $O(N^{-1/2})$ is not computed. The correction inside the bracket is $\operatorname{Re}[e^{i(\sqrt{2\pi N}-\pi/4)}\varpi_2]/\sqrt N$, with $\varpi_2=\sqrt2[(12+33\pi+2\pi^2)+i(12-33\pi+2\pi^2)]/(48\sqrt\pi)$ (derived; checked numerically on $S_k$ up to $k=6.5\cdot10^6$ and on $M_k-M_k^{\rm main}$ for $k\le400$). The $O(e^{-\sqrt{2\pi N}})$ in \eqref{eq:Mkfull} refers to this exact $M_k^{\rm main}$, not to a truncation of the first line; it can now be written $O(N^{-1}e^{-\sqrt{2\pi N}})$.

Smaller points: (a) the original wording "checked numerically to k = 6.5·10⁶" was a check of S_k, not of M_k; the repair says so. (b) Making the non-effective O-constant explicit is honest and costs one clause. (c) The consequential edits (Status row → Proved; Future Work item reduced) are justified.

## 8. Numerical summary (own code)

| k | rest_k (Arb, data M_k − own M_main) | rest/S_k − 1 | N e^{θ₂}rest/√π | ϖ₂ law | \|U_k\|/(C_U bound) |
|---|---|---|---|---|---|
| 100 | 1.692524470668146e−13 | −2.939e−4 | 0.840863 | 0.835588 | 1.3e−12 |
| 200 | −3.811407938780106e−18 | 3.48e−5 | −1.173506 | −1.153193 | 2.2e−13 |
| 300 | 1.927871694236297e−22 | −1.09e−5 | 0.252250 | 0.258636 | 1.7e−14 |
| 400 | 4.799125365553638e−25 | −1.17e−5 | 0.684268 | 0.684319 | 5.8e−14 |

All 13 rows (k = 100, 125, …, 400) agree with the RSTS table to every printed digit (`r2_rest_table.log`). rest radii ≤ 2·10^{−31}.

## Files

- `rcommon.py`: Arb integrators for M_k^main and S_k on the π/4 ray with hand head/tail bounds.
- `r1_CU_arb.py/.log`: C_U, exponent comparisons, index identity.
- `r2_rest_table.py/.log`: k = 100…400 table (rigorous balls).
- `r3_Sk_large.py/.log`: Arb S_k, k = 100, 400 (larger k too slow in Arb; aborted).
- `r3b_Sk_large_mp.py/.log`: mpmath S_k to k = 6.55·10⁶, Richardson ϖ₂.
- `r4_A2_and_U_direct.py/.log` (+ `r4_time.log`): A^{[2]} symbolic and numeric, (B3) test, direct U_100.
- `r5_varpi2.py/.log`: independent ϖ₂ and control ϖ₁.
