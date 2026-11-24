report = r"""# FN2R referee report — FN2 "Ahmed's Problem 3 at every angle" (proof.md)

**Under review:** `research/08_theory/FN2_every_angle_proofs/proof.md` (Theorem 1: every-angle
expansion with displayed 𝒥₁,𝒥₂,𝒥₃; Theorem 2: contact law with refined correction;
Lemmas A–E; Proposition), its `result.md`, `paper_proposal.tex`, `verify_fn2.py`; the
θ=π/2 template `research/07_push/P3_fn_graph/claude/result.md` ([P3]); the paper's
`sec:fn`, `thm:P2`, `thm:P2asym`.
**This folder:** independent code `step_04..step_41`, consolidated `verify_fn2r.py` (< 5 min).
Every number below was printed by my own code (60-digit mpmath unless stated).

## 1. Verdict table

| item | verdict | evidence (my steps) |
|---|---|---|
| Lemma A (layer cake, mass 1−c) | CORRECT (PROVED) | hand-checked; engine cross-check bump vs direct length: 1.2e-39 |
| Lemma B(i) (φ_θ forms, monotone, mass) | CORRECT | step_04: sympy exact (rationalization, φ′<0, φ(0)=1−c) |
| Lemma B(ii) (s⁻² bound) | CORRECT | step_04b: ratio to bound < 1 at 4 c-values, s=1..100 |
| Lemma B(iii) (tail series, d_l) | CORRECT | step_04: d_l=(−1)ˡ(1−c²)P′_{l−1}(c)/l exact l≤12; Legendre identity exact k≤9; c=0→−C(−1/2,m) exact m≤6; remainder ~ next term (step_04b) |
| Lemma C(i) (moment range −1<w<1) | CORRECT | step_07: cutoff growth at w=±1; all log-moments finite — this is the input used everywhere |
| Lemma C(ii) (Mellin=Legendre; k₁,k₂,k₃) | CORRECT | step_07: ≤7.2e-59 (w=±0.1, 5 angles), 2.5e-60 (w=±0.9 smoothed); cumulants vs closed forms ≤5.5e-46 (7 angles incl. π/12, 11π/12); k_j=(−1)ʲj!a_j consistent 1.1e-50; π/2 reduction exact |
| Lemma D (window, θ-independent) | CORRECT | restated from [P3] Prop. 2; the only θ-line (σ_n cutoff) correctly carries φ_θ≤1−c |
| Lemma E (bivariate inverse Θ, θ-independent) | CORRECT | step_05a: my own series reversion gives FN2's θ_{ij} exactly; θ_{1j}=a⁻¹(−A/a²)ʲ; θ_{k0}=Acal-reversion; extended to order 4 (step_12a) matching [P3]'s θ_{40},θ_{31},θ_{22},θ_{13} exactly |
| Theorem 1 (expansion at every angle) | CORRECT (PROVED) | line-by-line §3 below; master→displayed exact at general c (step_05b); π/2→[P3] exact; thm:P2 exact series (f=x) through ε³ at general c = displayed 𝒥's, difference 0 (step_06); numerics step_08 |
| Corollary (displayed 𝒥₁,𝒥₂,𝒥₃) | CORRECT | same checks; T₂=(1−c²)/(2a²), T₃=−(4A+ac)(1−c²)/(4a⁴) re-derived |
| Theorem 2 (contact), steps 1–4, 6 | CORRECT | step_09/09b/10/13 numerics; hand audit §3 |
| Theorem 2, step 5 (inner-root bound) | **GAP (minor) — conclusion correct, justification incomplete for 1<m<2** | the line "density O(s⁻²) ⇒ V-term O(ΛN^{−min(1,1/(m−1))})" silently uses E[S^{1/(m−1)}]<∞, which fails for 1<m≤2 (Lemma C(i)). Repair below; conclusions unchanged (verified numerically at m=3/2: step_10) |
| Theorem 2 (a),(b),(c) as stated | CORRECT after repair | (a): exact-Φ models match D to ~1e-33 (step_10b); (b): D/L3→1, τ⁺(e^{−k₁}) form within O(logΛ/Λ³) (step_13); (c): refined laws, scaled residuals bounded (step_09) |
| Proposition (u* limit, Klayer) | CORRECT | step_11: Klayer=(1−c)k₁/2 to 1e-29 (IBP by hand); limit →(1−c)k₁/2 (slow, as proved); AH3 guess off by exactly (1−c²)/4 |

## 2. Overall verdicts

* **Theorem 1: PROVED.** The θ-audit is complete: the three θ-sensitive inputs (layer-cake
  density, log-moments/cumulants, density tail series) are proved in Lemmas A–C; the window
  and the bivariate inverse are genuinely θ-free (Lemma D's equation contains only f and n;
  Θ solves Acal(Θ)−μB̃(Θ)=W with Acal,B̃ from f alone). No repair needed.

* **Theorem 2: PROVED AFTER REPAIRS.** One justification gap (step 5, 1<m<2), repair:

  **Repair of step 5 (inner root, 1<m≤2).** Replace the moment argument by the two-sided
  density bound for ξ=ξ_θ: f_ξ(t)=φ_θ(e^{−t})e^{−t}/(1−c) satisfies f_ξ(t)≤Ce^{−|t|} for
  all t (for t≥0 use φ_θ≤1−c; for t<0 use Lemma B(ii), φ_θ(e^{−t})≤C₂e^{2t}). With
  V^{1/m}(y)≤C e^{−y/(m−1)} and ℓ(y)=0 for y<y_m,
    E[V^{1/m}(Λ+ξ)] ≤ C∫_{y_m}^∞ V^{1/m}(y)e^{−|y−Λ|}dy
      = e^{−Λ}∫_{y_m}^{Λ}V^{1/m}(y)e^{y}dy + O(e^{−Λ/(m−1)}).
  For 1<m<2 the integral ∫V^{1/m}e^y dy converges (V^{1/m}e^y ~ Ce^{−y(2−m)/(m−1)}), so the
  first term is O(e^{−Λ})=O(N⁻¹); for m≥2 it is O(e^{Λ(1−1/(m−1))}), giving
  O(e^{−Λ/(m−1)})=O(N^{−1/(m−1)}). In both cases the V-term is
  O(N^{−min(1,1/(m−1))}) — exactly the bound [P3] and FN2's audit claim — hence power-small
  against F(Λ)~Λ^{1/m}, and (b) stands as stated. (For m<1 there is no inner root.)
  *Numerically the conclusion in the gap regime is confirmed: m=3/2, exact-Φ model:
  |1−D/Fm|≤9.7e-35; refined law scaled residual −0.13,−0.08 (bounded); D/L3→1.*

  No other moment is ever needed beyond Lemma C(i): every other use is a log-moment
  (all finite) or the two tails P(S>s)=O(s⁻¹), P(S≤s)≤s.

## 3. Line-by-line referee notes (focus points)

**Interchanges of limit/sum/integral.** (i) Layer cake: Tonelli, integrand ≥0 — fine.
(ii) Taylor polynomial (finite sum) integrated against φ_θ — trivial. (iii) Taylor remainder
integral: O(μ^{K+1}(1−c)E(1+|log(na/S_θ)|)^{K+1}) = O(μ^{K+1}log^{K+1}n); only log-moments —
finite by Lemma C(i). Uniformity in s∈[σ_n,na] holds because Θ is C^{K+1} on the compact
segment [0,W_s] (W_s≤log(1/θ_f)+o(1); strict monotonicity of Acal−μB̃ on [0,η] gives the
extension of Θ beyond the IFT neighbourhood — [P3]'s remark is adequate). (iv) Tail series
termwise integration on [na,∞): the power series in 1/s has radius 1 and converges uniformly
on s≥na≥2; with the logʲ factor use the truncated-L form φ=Σ_{l≤L}d_ls^{−l}+O_L(s^{−L−1})
and integrate the finite sum; remainder O((na)^{−L}logʲ) — fine; ∫_{na}^∞s^{−l}logʲ(na/s)ds
=(na)^{1−l}(−1)ʲj!/(l−1)^{j+1} checked. (v) Extension ∫_{σ_n}→∫₀: ≤(1−c)σ_n(C+|logσ_n|)^K,
exponentially small — checked ( ∫₀^{σ_n}logʲ(na/s)ds = σ_nΣ C(j,r)log^{j−r}(na/σ_n)r! ).
(vi) Thm 2 step 4 sandwich: F monotone, F′=O(Λ^{1/m−1}); the √-singularity of ℓ at y_m meets
the density of ξ at y_m−Λ, which is O(e^{−Λ}) — checked via f_ξ(t)~Ce^t as t→−∞.
(vii) Thm 2 step 6: third-order Taylor of F on {ξ>−Λ/2}, remainder CΛ^{1/m−3}|ξ−k₁|³ with
E|ξ−k₁|³<∞ (log-moment); complement O(Λ^{1/m}e^{−Λ/2}); Cauchy–Schwarz for the polynomial
terms — fine.

**Moment conditions.** E[S^w]<∞ exactly for −1<w<1 (Lemma C(i), verified numerically).
The only higher-moment use in the whole note is the step-5 line flagged above. Everything
else uses log-moments or the two tails — no other gap.

**Error terms.** Thm 1: genuinely O(ε^{K+1}log^{K+1}n): the Taylor-remainder integral gives
exactly that; the tail-series truncation gives O(μ^{K+1}) (no logs) and is absorbed; cutoffs
exponentially small. Implied constants depend on θ and the germ of f — they blow up as
θ→0,π (Lemma B(ii) constant (1+|c|)/(2(1−|c|))); the theorems are for fixed θ, as stated.
Thm 2(c): the O((loglog n/log n)²) follows from the Λ-form by
Λ=(log n+log(m^mβ))/m, logΛ=loglog n−log m+O(1/log n), and
log(m^mβ)−(m−1)log m=log(mβ) — re-derived by hand, matches; scaled residuals bounded
(see tables).

**θ-independence of window and inverse.** Confirmed (Lemma D equation and Θ contain only
f, n). Also the factor (na)^{1−l}=(μ/(a(1+μ)))^{l−1} and the sign (−1)^{i+1} in T — checked.

**Thm 2 specifics.** Lambert-W forms of U (branches W₋₁/W₀ for m>1; W₀ for m<1) re-derived
from u−c_m log u=y and used in code; the derivative hypothesis is used exactly in step 1's
|ε_n(σ)|≤Cσ^δ(1+σ^m)N^{−δ} (re-derived); m<1 case (single root, V:=0, y_m=−∞, f' (1)=∞,
y only absolutely continuous — layer cake still applies) tested numerically (m=1/2 below);
approach from above for m>1 (D/L0>1 at all 12 (θ,n) points) and from below for m<1
(D/L0=0.68,0.75) — confirmed.

## 4. Independent re-derivations (sympy, exact unless stated)

| check | result | step |
|---|---|---|
| φ_θ rationalization, φ′<0, φ(0)=1−c | exact | step_04 |
| d_l(c)=(−1)ˡ(1−c²)P′_{l−1}(c)/l, l≤12; d₂=2z(1−z), d₃=−4cz(1−z), d₄=3(1−c²)(5c²−1)/8; c=0→−C(−1/2,m), m≤6; Legendre identity k(P_k−xP_{k−1})=(x²−1)P′_{k−1}, k≤9 | exact | step_04 |
| θ_{ij}, i+j≤3 (my own series reversion) = FN2's = [P3]'s; θ_{1j}=a⁻¹(−A/a²)ʲ | exact | step_05a |
| θ_{40},θ_{31},θ_{22},θ_{13} (my reversion) = [P3]'s quoted values | exact | step_12a |
| T₂=(1−c²)/(2a²), T₃=−(4A+ac)(1−c²)/(4a⁴), my T₄(c); T₄(0)=P3's T₄ | exact | step_05b, step_12a |
| master formula → paper's displayed 𝒥₁,𝒥₂,𝒥₃ at general c (z=(1−c)/2): difference 0 | exact | step_05b |
| π/2 reduction to [P3]'s Q₁,Q₂,Q₃: difference 0 | exact | step_05b |
| **thm:P2 exact series (f=x) expanded through ε³ at general c = displayed 𝒥's: difference 0** (strongest end-to-end check: cumulant conversion k_j=(−1)ʲj!a_j, M_i structure, tail T₂,T₃, ε=1/(n−1) bookkeeping all at once) | exact | step_06 |
| Mellin E[S^w]=(πw/sinπw)·2z₂F₁/(1−c): residuals ≤7.2e-59 (w=±0.1), 2.5e-60 (w=±0.9, smoothed), 5 angles | 59–60 digits | step_07 |
| k₁,k₂,k₃(θ) vs closed forms: worst 5.5e-46 over 7 angles (incl. π/12, 11π/12); k_j vs (−1)ʲj!a_j: 1.1e-50; π/2 values exact; k₄(π/2)=6+π⁴/15: 6.4e-50 | 46–50 digits | step_07, step_13 |

## 5. Independent numerics (own engines)

**Engine validation:** bump ∫h_c(y')dx vs direct length 2−L_n−f(0)ⁿ agree to 1.2e-39
(step_08c); my exp-f residuals reproduce FN2's step_39 numbers to all printed digits
(cross-validation of independent engines).

**Residual after 𝒥₃, scaled by ε⁴log⁴n** (task: bounded). n=10³…10⁶:

| f | θ | n=10³ | 10⁴ | 10⁵ | 10⁶ |
|---|---|---|---|---|---|
| (eˣ−1)/(e−1) | π/3 | −0.000378 | 0.000848 | 0.001344 | 0.001585 |
| (eˣ−1)/(e−1) | 2π/3 | −0.001109 | 0.002987 | 0.004536 | 0.005245 |
| (eˣ−1)/(e−1) | π/2 | −0.000786 | 0.001807 | 0.002823 | 0.003304 |
| (1+x)/2·e^{x−1} (a=3/2,A=2,B=5/2) | π/3 | −0.001278 | −0.000548 | −0.000277 | −0.000158 |
| (1+x)/2·e^{x−1} | 2π/3 | −0.004567 | −0.001868 | −0.000918 | −0.000515 |
| (1+x)/2·e^{x−1} | π/2 | −0.002814 | −0.001176 | −0.000585 | −0.000332 |

All bounded (and ≈𝒥₄(0)/log⁴n + o(1)); the exp-f values equal FN2's to all digits.
**Refined:** resid/(ε⁴𝒥₄) with my own θ_{4j}, T₄(c), k₄: 0.99990, 0.999996, 1.0000001
(exp, π/3, n=10⁴,10⁵,10⁶), 0.99990→0.999998 (linexp) — resid=ε⁴𝒥₄(1+o(1)) confirmed,
so the error term O(ε⁴log⁴n) has the right shape.

**Contact law, f=sin(πx/2) (m=2, β=π²/8), θ=π/3, 2π/3 (control π/2), n=10⁴…10⁷:**
D/L₀ at n=10⁵: 1.143549, 1.181202 (control 1.159762) — exactly AH3's printed values.
D/full-model (a): 0.9999141 → 0.99999985 (π/3), error/(log n/n)=0.093→0.090 (sharper
O(δ_N/Λ) with δ=m confirmed). D/L₃ (law (b)): 0.98998→0.99935. Refined-law scaled
residuals: (D/L₀−pred_Λ)(Λ/logΛ)² ∈ [0.70,0.85] (π/3: 0.704,0.746,0.768,0.781);
(D/L₀−pred_logn)(log n/loglog n)² ∈ [−0.44,−0.19] — both bounded, so the
O((loglog n/log n)²) remainder is genuine (even Θ of it). D/L₀>1 throughout: approach
from above for m=2>1.

**Gap regime m=3/2 and infinite-slope m=1/2 (exact-Φ models, θ=π/3, π/2):**
(a) |1−D/Fm| ≤ 4.3e-33 (m=1.5: ≤9.7e-35; m=0.5: ≤4.3e-33; m=2: ≤6.9e-35; m=3: ≤5.2e-35)
— the Lambert-W reduction is essentially exact for exact Φ. m=3/2: D/L₀=1.117,1.090,
scaled residual −0.132,−0.083, D/L₃→1 — the theorem holds where the step-5 justification
needed the repair. m=1/2: D/L₀=0.684,0.752 < 1 (approach from below, as predicted for
m<1), scaled residual 1.05,0.90 (bounded), D/L₃=0.999993,0.9999993.

**τ⁺(e^{−k₁(θ)}) form of (b):** D/((1−c)τ⁺(e^{−k₁})) vs 1−(m−1)k₂/(2m²Λ²): differences
−0.0033,−0.0020,−0.0006 (n=10⁵) shrinking to −0.0011,−0.0005,+0.0001 (n=10⁶), inside the
stated O(logΛ/Λ³) (0.0067, 0.0045).

**Proposition:** Klayer=(1−c)k₁/2: diffs −9.8e-29 (π/3), 5.7e-30 (π/2), 2.0e-10 (2π/3,
quadrature-limited). (I(V,c)−(1−c)u*)u* → (1−c)k₁/2: at V=10⁸ residuals −0.0096..−0.0138
and shrinking (log-slow, sign from F″<0 as in proof); AH3's guess Klayer−(1−c²)/4 is off
by exactly (1−c²)/4 (π/3: −0.1875=−3/16) — FN2's correction is right.

## 6. Minor remarks (not affecting verdicts)

1. proof.md's Theorem 2 statement could flag that for 1<m<2 the step-5 bound needs the
   two-sided density argument (repair in §2). One-line fix in the audit table, row 5.
2. The B(ii) constant and all implied constants blow up at θ→0,π — consistent with
   fixed-θ statements (already noted in FN2's result.md).
3. My own first linexp run used wrong derivatives (a,A,B)=(2,5/2,3) instead of
   (3/2,2,5/2) and produced fake residuals ~ε·log n — a reminder that the expansion is
   extremely sensitive to the derivative values; with correct values everything bounds.
   (My bug, not FN2's; recorded for transparency.)
4. FN2's verify S6 "pred" is the Λ-form of the refined law; at n=10⁵ the logn-form and
   Λ-form differ by O((loglog/log)²) ≈ 6% — both are valid truncations; both bounded
   after scaling (my step_09).
5. Klayer quadrature at θ=2π/3 limited to 2e-10 (tail ∫₁^∞ρ/s) — the IBP identity is
   proved in proof.md anyway and holds to 1e-29 at the other angles.

## 7. Files

step_04 (Lemma B symbolic), step_04b (bounds), step_05a (θ_{ij} reversion), step_05b
(master→displayed, π/2 reduction), step_06 (thm:P2 series f=x end-to-end), step_07/07b
(Lemma C numerics), step_08/08b/08c (defect engine + residuals + validation),
step_09/09b (contact sin), step_10/10b (m=3/2, 1/2, 2, 3 exact-Φ), step_11 (Proposition),
step_12a/12b (order-4 coefficients, resid=ε⁴𝒥₄), step_13 (k₄ sanity, τ⁺ form),
fn2r_residuals*.json, fn2r_contact.json, verify_fn2r.py (consolidated, < 5 min).
"""
open("result.md","w").write(report)
import time
with open("progress.log","a") as f:
    f.write(f"{time.strftime('%H:%M')} result.md written: verdicts Theorem 1 PROVED, Theorem 2 PROVED AFTER REPAIRS (step-5 inner-root justification for 1<m<2; repair written out) [referee]\n")
print("result.md written", len(report))
