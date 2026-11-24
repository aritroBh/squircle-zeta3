# FN2R referee report — FN2 "Ahmed's Problem 3 at every angle"

**Verdicts: Theorem 1 — PROVED. Theorem 2 — PROVED AFTER REPAIRS (one minor justification gap, repair written out below; conclusions unchanged). Lemmas A–E, the Corollary (displayed 𝒥₁,𝒥₂,𝒥₃), and the Proposition — CORRECT.** Every claim was re-checked with independent code (my own derivations and engines, in `new_results/FN2R_referee_every_angle/`); nothing is WRONG.

## 1. Verdict table

| item | verdict | key evidence (step files) |
|---|---|---|
| Lemma A (layer cake, mass 1−cos θ) | CORRECT (PROVED) | Tonelli, hand-checked; engine cross-check bump vs direct length 1.2e-39 (step_36–41) |
| Lemma B(i) (φ_θ forms, monotone, mass) | CORRECT | sympy exact (step_14) |
| Lemma B(ii) (φ_θ ≤ (1+\|c\|)/(2(1−\|c\|)) s⁻²) | CORRECT | ratio to bound < 1, 4 c-values (step_15) |
| Lemma B(iii) (tail series, d_l = (−1)ˡ(1−c²)P′_{l−1}(c)/l) | CORRECT | exact l≤12; Legendre identity exact; c=0 → −C(−1/2,m) exact; remainder ~ next term (step_14, 15) |
| Lemma C(i) (E[S^w]<∞ iff −1<w<1; all log-moments) | CORRECT | cutoff growth at w=±1 (step_28) |
| Lemma C(ii) (Mellin=Legendre; k₁,k₂,k₃ closed forms) | CORRECT | Mellin residuals ≤7.2e-59 (w=±0.1), 2.5e-60 (w=±0.9 smoothed); cumulants ≤5.5e-46, 7 angles incl. π/12, 11π/12; π/2 values exact (step_28, 32–34) |
| Lemma D (window, θ-independent) | CORRECT | contains only f, n; the one θ-line (cutoff) correctly carries φ_θ ≤ 1−c |
| Lemma E (bivariate inverse Θ, θ-independent) | CORRECT | my own series reversion reproduces θ_{ij} exactly through order 4 (step_22–24, 50) |
| **Theorem 1** (expansion at every angle) | **CORRECT (PROVED)** | master→displayed exact at general c (step_26); thm:P2 series (f=x) through ε³ at general c = displayed 𝒥's, diff 0 (step_27); numerics (step_36–41) |
| Corollary (displayed 𝒥₁,𝒥₂,𝒥₃) | CORRECT | T₂=(1−c²)/(2a²), T₃=−(4A+ac)(1−c²)/(4a⁴) re-derived; π/2 → [P3]'s Q's diff 0 |
| **Theorem 2**, steps 1–4, 6 | CORRECT | numerics (step_43–46, 49, 51) |
| **Theorem 2**, step 5 (inner-root bound) | **GAP (minor)** — justification incomplete for 1<m<2; conclusion correct | the line "density O(s⁻²) ⇒ V-term O(ΛN^{−min(1,1/(m−1))})" silently needs E[S^{1/(m−1)}]<∞, false for 1<m≤2. Repair below; regime tested at m=3/2 (step_49, 51) |
| Theorem 2 (a),(b),(c) as stated | CORRECT after repair | (a): exact-Φ models match D to ~1e-33; (b): D/L₃→1, τ⁺(e^{−k₁}) form inside O(logΛ/Λ³); (c): refined-law scaled residuals bounded |
| Proposition (u* limit; Klayer=(1−c)k₁/2) | CORRECT | IBP identity to 1e-29; AH3's guess off by exactly (1−c²)/4 (step_53) |

## 2. Overall verdicts and the repair

**Theorem 1: PROVED.** The θ-audit is complete: the only θ-sensitive inputs (layer-cake density, log-moments/cumulants, density tail series) are proved in Lemmas A–C; the window {y′>s} and the inverse Θ are genuinely θ-free. All interchanges are justified: Tonelli (layer cake); finite Taylor sum; remainder integral dominated by log-moments only (all finite by Lemma C(i)) giving exactly O(ε^{K+1}log^{K+1}n); tail series integrated termwise on s ≥ na ≥ 2 via the truncated form φ = Σ_{l≤L}d_ls^{−l}+O_L(s^{−L−1}); cutoffs exponentially small. No higher moment E[S^w] (w≥1) is ever used. Constants blow up only as θ→0,π — the statement is for fixed θ, as claimed.

**Theorem 2: PROVED AFTER REPAIRS.** The repair (step 5 of [P3]'s Theorem 4, inherited verbatim by FN2's audit, row 5):

> **Repair (inner-root bound for 1<m≤2).** Replace the moment argument by the two-sided density bound for ξ=ξ_θ: f_ξ(t)=φ_θ(e^{−t})e^{−t}/(1−c) ≤ Ce^{−|t|} for all t (t≥0: φ_θ≤1−c; t<0: Lemma B(ii) gives φ_θ(e^{−t})≤C₂e^{2t}). With V^{1/m}(y)≤Ce^{−y/(m−1)} and ℓ=0 for y<y_m,
> E[V^{1/m}(Λ+ξ)] ≤ C∫_{y_m}^∞ V^{1/m}(y)e^{−|y−Λ|}dy = e^{−Λ}∫_{y_m}^{Λ}V^{1/m}(y)e^y dy + O(e^{−Λ/(m−1)}).
> For 1<m<2, ∫V^{1/m}e^y dy converges (V^{1/m}e^y ~ Ce^{−y(2−m)/(m−1)}), so the first term is O(e^{−Λ})=O(N⁻¹); for m≥2 it is O(e^{−Λ/(m−1)})=O(N^{−1/(m−1)}). Hence the V-term is O(N^{−min(1,1/(m−1))}) — exactly the bound claimed — power-small against F(Λ)~Λ^{1/m}, so (b) and (c) stand as stated.

All other Thm-2 focus points check out: the Lambert-W forms (W₋₁/W₀ for m>1; W₀ for m<1) re-derived; the hypothesis with derivative is used exactly in step 1's |ε_n(σ)|≤Cσ^δ(1+σ^m)N^{−δ}; the m<1 case (single root, V:=0, y_m=−∞, f absolutely continuous) works; the refined law's O((loglog n/log n)²) follows via Λ=(log n+log(m^mβ))/m and log(m^mβ)−(m−1)log m=log(mβ) (re-derived), and is numerically genuine (scaled residuals bounded away from 0 and ∞); approach from above for m>1 and from below for m<1 — both confirmed.

## 3. Evidence table (all numbers printed by my code; 60-digit arithmetic unless stated)

| check | residual / value | digits | step |
|---|---|---|---|
| d_l(c) Legendre form, l≤12; c=0→−C(−1/2,m); Legendre identity | 0 (sympy exact) | exact | step_14 |
| θ_{ij} (own reversion) vs FN2/[P3], orders 3 and 4 | 0 (exact) | exact | step_22–24, 50 |
| master formula → displayed 𝒥₁,𝒥₂,𝒥₃ at general c | 0 (exact) | exact | step_26 |
| π/2 reduction to [P3]'s Q₁,Q₂,Q₃ | 0 (exact) | exact | step_26 |
| thm:P2 exact series (f=x) through ε³ = displayed 𝒥's, general c | 0 (exact) | exact | step_27 |
| Mellin E[S^w] = (πw/sinπw)·2z₂F₁/(1−c), w=±0.1 … ±0.9, 5 angles | ≤7.2e-59; 2.5e-60 | 59–60 | step_28, 32–34 |
| k₁,k₂,k₃(θ) vs closed forms, 7 angles (incl. π/12, 11π/12) | ≤5.5e-46 | 46 | step_28 |
| k_j = (−1)ʲj!a_j (paper's a_j) | 1.1e-50 | 50 | step_28 |
| k₄(π/2) = 6+π⁴/15 | 6.4e-50 | 50 | step_51 |
| engine: bump ∫h_c(y') vs direct 2−L_n−f(0)ⁿ | 1.2e-39 | 39 | step_36–41 |
| resid after 𝒥₃ = ε⁴𝒥₄(1+o(1)) (own θ_{4j}, T₄(c), k₄) | ratios 0.99990→1.00000 | — | step_54, 56 |

**Residuals after 𝒥₃, scaled by ε⁴log⁴n (bounded as required; n=10³…10⁶):**

| f | θ | 10³ | 10⁴ | 10⁵ | 10⁶ |
|---|---|---|---|---|---|
| (eˣ−1)/(e−1) | π/3 | −0.000378 | 0.000848 | 0.001344 | 0.001585 |
| (eˣ−1)/(e−1) | 2π/3 | −0.001109 | 0.002987 | 0.004536 | 0.005245 |
| (1+x)/2·e^{x−1} (a=3/2, A=2, B=5/2) | π/3 | −0.001278 | −0.000548 | −0.000277 | −0.000158 |
| (1+x)/2·e^{x−1} | 2π/3 | −0.004567 | −0.001868 | −0.000918 | −0.000515 |

(exp-f values equal FN2's step_39 numbers to all their printed digits — independent engines cross-validate.)

**Contact law, f=sin(πx/2), m=2, θ=π/3, 2π/3 (control π/2), n=10⁴…10⁷:** D/L₀ at n=10⁵: 1.143549, 1.181202 (control 1.159762) — exactly AH3's printed values; D/full-model: 0.9999141→0.99999985 with error ≈ 0.092·log n/n (the sharper O(δ_N/Λ), δ=m); D/L₃: 0.98998→0.99935; refined-law scaled residuals (D/L₀−pred_Λ)(Λ/logΛ)² ∈ [0.70,0.85], (D/L₀−pred_logn)(log n/loglog n)² ∈ [−0.44,−0.19] — both bounded; D/L₀>1 throughout (approach from above).

**Gap regime and m<1 (exact-Φ models, θ=π/3):** |1−D/Fm| = 9.7e-35 (m=3/2), 4.3e-33 (m=1/2), 6.9e-35 (m=2), 5.2e-35 (m=3); m=3/2 refined-law scaled residual −0.132, −0.083 (bounded); m=1/2: D/L₀ = 0.684, 0.752 < 1 (approach from below, as predicted).

**Proposition:** Klayer−(1−c)k₁/2 = −9.8e-29 (π/3), 5.7e-30 (π/2); (I−(1−c)u*)u* → (1−c)k₁/2 (slow, from below as F″<0); AH3's guess off by exactly (1−c²)/4 (π/3: −0.1875 = −3/16).

## 4. Files produced

`result.md` (full report), `verify_fn2r.py` + `verify_fn2r.out` (**86 s < 5 min; R1–R9 all pass**), `progress.log`, and the step files: step_14 (Lemma B symbolic), step_15 (bounds), step_22–24 (θ_{ij} reversion), step_26 (master→displayed, π/2 reduction), step_27 (thm:P2 series end-to-end), step_28+32–34 (Lemma C numerics), step_36–41 (defect engine, residuals, validation), step_43–46 (contact sin), step_49+51 (m=3/2, 1/2, 2, 3 exact-Φ), step_53 (Proposition), step_50+54+56 (order-4, resid=ε⁴𝒥₄), step_51 (k₄, τ⁺ form); data: fn2r_residuals*.json, fn2r_contact.json.

Transparency notes (all three were my own bugs, caught and fixed; none affect FN2): (i) my first m-sum sign error in step_06; (ii) wrong (a,A,B) for (1+x)/2·e^{x−1} (correct: a=3/2, A=2, B=5/2) producing fake residuals — the expansion is extremely sensitive to derivative values; (iii) two quadrature settings (unresolved √-singularity of ℓ at y_m; v-domain leaks past t=1) — after fixes, everything is consistent with FN2's numbers.