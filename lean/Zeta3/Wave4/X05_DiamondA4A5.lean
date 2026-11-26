import Mathlib
/-!
# X05 — `a_4` and `a_5` at the diamond, as genuine derivatives of the perimeter

Target id: X05.
Paper labels: `prop:an` (main.tex lines 1163–1174), `thm:pw` (main.tex lines 1183–1195),
the derivation "Steps 1–4" of Section `sec:diamond` (main.tex lines 1099–1131), the
perimeter integral `eq:Pint` (main.tex lines 825–827), and the detailed proof of `thm:pw` in
paper/ieee/proofs/sec_pureweight.tex (Lemmas `lem:pw:logA`, `lem:pw:remainder`,
`lem:pw:asymp`, `prop:pw:rhombus` with `q = 0`, and the table of Section `sec:pw:series`).

Quoted verbatim.

main.tex 825–827:

  parametrize the quarter arc in the first quadrant by $x=t^{1/p}$, $y=(1-t)^{1/p}$,
  $0\le t\le1$. Then $x^p+y^p=1$, and the arc-length formula gives
  P(p)=\frac4p\int_0^1\sqrt{t^{2/p-2}+(1-t)^{2/p-2}}\;dt .

main.tex 1102–1105 (Step 1):

  Put $\eta=\frac1p-1$ ... Let's factor out $[t(1-t)]^{\eta/2}$ and use the logit coordinate
  $l=\log\frac{t}{1-t}$. ...
  P(p)=4\sqrt2\,(1+\eta)\int_0^1\bigl(t(1-t)\bigr)^{\eta/2}\sqrt{\cosh(\eta l)}\;dt
  exactly, for every $p>0$. At $\eta=0$ the integrand is $1$ and $P(1)=4\sqrt2$.

main.tex 1108–1112 (Step 2):

  A(\eta,s)&=(1+\eta)\int_0^1\bigl(t(1-t)\bigr)^{\eta/2}e^{sl}\dd t\\
  &=\frac{\Gamma\bigl(1+\frac\eta2+s\bigr)\,\Gamma\bigl(1+\frac\eta2-s\bigr)}{\Gamma(1+\eta)} .

main.tex 1119–1124 (Step 3):

  \log\Gamma(1+x)=-\gamma x+\sum_{n\ge2}\frac{(-1)^n\zeta(n)}{n}x^n\qquad(|x|<1).
  ...
  \log A(\eta,s)=\sum_{n\ge2}\frac{(-1)^n\zeta(n)}{n}\Bigl[\Bigl(\frac\eta2+s\Bigr)^n
     +\Bigl(\frac\eta2-s\Bigr)^n-\eta^n\Bigr]

main.tex 1164–1168 (`prop:an`):

  Let $a_n=P^{(n)}(1)/(\sqrt2\,n!)$, the Taylor coefficients of $P(1+d)/\sqrt2$ at $d=0$
  (the series diverges, by Theorem~\ref{thm:zero}). Then $a_0=4$, $a_1=0$, $a_2=\pi^2/6$,
  $a_3=\zeta(3)-\pi^2/3$,
  \begin{gather*}
  a_4=-3\zeta(3)-\frac{11\pi^4}{480}+\frac{\pi^2}{2},\\
  a_5=-\frac{7\pi^2\zeta(3)}{24}+6\zeta(3)-\frac{3\zeta(5)}{4}-\frac{2\pi^2}{3}+\frac{11\pi^4}{120},
  \end{gather*}

main.tex 1186–1189 (`thm:pw`, the displayed series; `\tau=1-\frac1p` from line 1181):

  \frac{P}{\sqrt2}=4+\zeta(2)\tau^2+\zeta(3)\tau^3-\frac{33}{16}\zeta(4)\tau^4\\
  -\Bigl(\frac34\zeta(5)+\frac74\zeta(2)\zeta(3)\Bigr)\tau^5+O(\tau^6).

sec_pureweight.tex (Lemma `lem:pw:remainder`, `lem:pw:asymp`, table of `sec:pw:series`):

  \Bigl|g(z)-\sum_{m<M}g_mz^m\Bigr|\le C_M|z|^M\eu^{|z|/2}.
  \frac{P_{A,B}(p)}{4r}=\sum_{m<M,\ m\ \mathrm{even}}g_m(q)\,\eta^m\,\partial_s^m\Psi(\eta,0)+O(|\eta|^M).
  For $q=0$: $g_0=1$, $g_2=\frac14$, $g_4=-\frac1{96}$.
  (table) weight 2: -\tfrac14Z_2, 2Z_2; weight 3: \tfrac14Z_3, -2Z_3;
  weight 4: \tfrac1{32}Z_2^2-\tfrac7{32}Z_4, -\tfrac12Z_2^2+\tfrac32Z_4, 12Z_2^2+12Z_4;
  weight 5: -\tfrac1{16}Z_2Z_3+\tfrac3{16}Z_5, Z_2Z_3-Z_5, -24Z_2Z_3-24Z_5

Formalization choices.
* `ζ(n)` is `(riemannZeta n).re`, the real part of Mathlib's complex `riemannZeta` (for
  `n ≥ 2` it is real and equals `∑ 1/k^n`, by `zeta_nat_eq_tsum_of_gt_one`). As the job allows,
  `ζ(2)` and `ζ(4)` are written `π²/6` and `π⁴/90` (Mathlib `riemannZeta_two`,
  `riemannZeta_four`); `ζ(3)`, `ζ(5)` are `z 3`, `z 5`; in stmt_3 the general `ζ(n)` is `z n`.
* `P` is the paper's `eq:Pint` with `Real.rpow`. For `p > 0` the integrand is
  `≤ t^{1/p-1} + (1-t)^{1/p-1}`, integrable, so the `intervalIntegral` is not a junk value;
  stmt_1 records this integrability. The conventions `0^x` and `log 0 = 0` affect only the
  endpoints `t = 0, 1` (measure zero).
* `A η s` is the paper's `A(η,s)` (supplement: `Ψ(η,s)`), defined by the integral; stmt_2 is
  the Beta evaluation, for `η > -1`, `|s| < 1 + η/2` (where the integral converges).
* Derivatives: `a n = iteratedDeriv n P 1 / (√2 · n!)`. Because `iteratedDeriv` is `0` at
  non-differentiable points, stmt_8 states that `P` is `C^n` at `1` for every `n`; this is the
  paper's "Differentiating under the integral is justified by an integrable bound"
  (proof of `thm:rhombus`), and it makes stmt_9–stmt_11 statements about genuine derivatives.
* All `O(·)` claims are explicit: `∃ C (and a radius δ)`, `|error| ≤ C · |var|^k` on the
  two-sided neighbourhood (the supplement's Lemma `lem:pw:asymp` is two-sided, `|η| ≤ 1/4`).
  In stmt_5 the point `p = 1/(1+η)` is the inverse of `η = 1/p - 1`.
* `∂_s^m A(η,0)` is `iteratedDeriv m (fun s => A η s) 0`.

Statements.
* stmt_1: Step 1, the exact representation, for every `p > 0` (plus integrability).
* stmt_2: Step 2, the Beta moments `A(η,s) = Γ(1+η/2+s)Γ(1+η/2-s)/Γ(1+η)`.
* stmt_3: Step 3 / Lemma `lem:pw:logA`, the series for `log A(η,s)`, `|η|,|s| < 1/2`.
* stmt_4: Lemma `lem:pw:remainder` at `q = 0`, `M = 6`:
  `|√cosh z - (1 + z²/4 - z⁴/96)| ≤ C |z|^6 e^{|z|/2}` for real `z`.
* stmt_5: Lemma `lem:pw:asymp` at `q = 0`, `M = 6`.
* stmt_6: the moment table (q = 0, weights ≤ 5) as explicit expansions of `A(η,0)`,
  `∂_s²A(η,0)`, `∂_s⁴A(η,0)` in `η`.
* stmt_7: the displayed `τ`-series of `thm:pw` through `τ^5` with an explicit `O(τ^6)`.
* stmt_8: `P` is `C^n` at `p = 1` for every `n`.
* stmt_9: `a_0 = 4`, `a_1 = 0`, `a_2 = π²/6`, `a_3 = ζ(3) - π²/3`.
* stmt_10: `a_4 = -3ζ(3) - 11π⁴/480 + π²/2`.
* stmt_11: `a_5 = -7π²ζ(3)/24 + 6ζ(3) - 3ζ(5)/4 - 2π²/3 + 11π⁴/120` (the final claim).

NOT COVERED: the parenthetical divergence of the Taylor series (`thm:zero`); `a_6`, `a_7`
(Table `tab:closed`); the general polynomial/weight claims and the `ω_n` formula of
`prop:an`; the general-order `log(P/4√2) = Σ w_k τ^k` purity claim of `thm:pw`, its rhombus
(`q ≠ 0`) and `N`-dimensional parts; the logarithmic form `w_2, …, w_5` (only the expansion of
`P/√2` itself is stated).
-/
-- BEGIN STATEMENTS
namespace X05

open Real

/-- `ζ(n)` as a real number: the real part of Mathlib's `riemannZeta n`. -/
noncomputable def z (n : ℕ) : ℝ := (riemannZeta n).re

/-- The perimeter integral `eq:Pint`:
`P(p) = (4/p) ∫_0^1 sqrt(t^{2/p-2} + (1-t)^{2/p-2}) dt`. -/
noncomputable def P (p : ℝ) : ℝ :=
  4 / p * ∫ t in (0 : ℝ)..1, Real.sqrt (t ^ (2 / p - 2) + (1 - t) ^ (2 / p - 2))

/-- The logit coordinate `l = log (t/(1-t))`. -/
noncomputable def logit (t : ℝ) : ℝ := Real.log (t / (1 - t))

/-- Step 2: `A(η,s) = (1+η) ∫_0^1 (t(1-t))^{η/2} e^{s l} dt`. -/
noncomputable def A (η s : ℝ) : ℝ :=
  (1 + η) * ∫ t in (0 : ℝ)..1, (t * (1 - t)) ^ (η / 2) * Real.exp (s * logit t)

/-- `a_n = P^{(n)}(1) / (√2 n!)` (`prop:an`). -/
noncomputable def a (n : ℕ) : ℝ := iteratedDeriv n P 1 / (Real.sqrt 2 * (n.factorial : ℝ))

/-- Step 1 (main.tex 1102–1105): for every `p > 0`, with `η = 1/p - 1`,
`P(p) = 4√2 (1+η) ∫_0^1 (t(1-t))^{η/2} sqrt(cosh(η l)) dt`; the integrand of `eq:Pint` is
integrable on `[0,1]`. -/
def stmt_1 : Prop :=
  ∀ p : ℝ, 0 < p →
    IntervalIntegrable (fun t : ℝ => Real.sqrt (t ^ (2 / p - 2) + (1 - t) ^ (2 / p - 2)))
      MeasureTheory.volume 0 1 ∧
    P p = 4 * Real.sqrt 2 * (1 + (1 / p - 1)) *
      ∫ t in (0 : ℝ)..1,
        (t * (1 - t)) ^ ((1 / p - 1) / 2) * Real.sqrt (Real.cosh ((1 / p - 1) * logit t))

/-- Step 2 (main.tex 1108–1112): the Beta moments
`A(η,s) = Γ(1+η/2+s) Γ(1+η/2-s) / Γ(1+η)` for `η > -1`, `|s| < 1 + η/2`. -/
def stmt_2 : Prop :=
  ∀ η s : ℝ, -1 < η → |s| < 1 + η / 2 →
    A η s = Real.Gamma (1 + η / 2 + s) * Real.Gamma (1 + η / 2 - s) / Real.Gamma (1 + η)

/-- Step 3 (main.tex 1121–1124; supplement Lemma `lem:pw:logA`): for `|η| < 1/2`, `|s| < 1/2`,
`log A(η,s) = Σ_{n≥2} (-1)^n ζ(n)/n [(η/2+s)^n + (η/2-s)^n - η^n]` (Euler's `γ` cancels). -/
def stmt_3 : Prop :=
  ∀ η s : ℝ, |η| < 1 / 2 → |s| < 1 / 2 →
    0 < A η s ∧
    HasSum (fun k : ℕ => (-1 : ℝ) ^ (k + 2) * z (k + 2) / ((k + 2 : ℕ) : ℝ) *
        ((η / 2 + s) ^ (k + 2) + (η / 2 - s) ^ (k + 2) - η ^ (k + 2)))
      (Real.log (A η s))

/-- Supplement Lemma `lem:pw:remainder` with `q = 0`, `M = 6` (`g_0 = 1`, `g_2 = 1/4`,
`g_4 = -1/96`, odd `g_m = 0`): `|√cosh z - (1 + z²/4 - z⁴/96)| ≤ C |z|^6 e^{|z|/2}`. -/
def stmt_4 : Prop :=
  ∃ C : ℝ, ∀ x : ℝ,
    |Real.sqrt (Real.cosh x) - (1 + x ^ 2 / 4 - x ^ 4 / 96)| ≤ C * |x| ^ 6 * Real.exp (|x| / 2)

/-- Supplement Lemma `lem:pw:asymp` with `q = 0`, `M = 6`: for `|η| ≤ 1/4` and `p = 1/(1+η)`,
`P(p)/(4√2) = A(η,0) + (1/4) η² ∂_s²A(η,0) - (1/96) η⁴ ∂_s⁴A(η,0) + O(|η|^6)`
(main.tex 1114–1116: "as an identity at every finite order"). -/
def stmt_5 : Prop :=
  ∃ C : ℝ, ∀ η : ℝ, |η| ≤ 1 / 4 →
    |P (1 / (1 + η)) / (4 * Real.sqrt 2) -
        (A η 0 + 1 / 4 * η ^ 2 * iteratedDeriv 2 (fun s => A η s) 0
          - 1 / 96 * η ^ 4 * iteratedDeriv 4 (fun s => A η s) 0)| ≤ C * |η| ^ 6

/-- The moment table of sec_pureweight.tex (Section `sec:pw:series`, from Lemma
`lem:pw:logA`), with `ζ(2) = π²/6`, `ζ(4) = π⁴/90`:
`A(η,0) = 1 - ζ2/4 η² + ζ3/4 η³ + (ζ2²/32 - 7ζ4/32) η⁴ + (-ζ2ζ3/16 + 3ζ5/16) η⁵ + O(η⁶)`,
`∂_s²A(η,0) = 2ζ2 - 2ζ3 η + (-ζ2²/2 + 3ζ4/2) η² + (ζ2ζ3 - ζ5) η³ + O(η⁴)`,
`∂_s⁴A(η,0) = 12ζ2² + 12ζ4 + O(η)` (main.tex 1128 gives the first two to lower order). -/
def stmt_6 : Prop :=
  ∃ C δ : ℝ, 0 < δ ∧ ∀ η : ℝ, |η| ≤ δ →
    |A η 0 - (1 - (π ^ 2 / 6) / 4 * η ^ 2 + z 3 / 4 * η ^ 3
        + ((π ^ 2 / 6) ^ 2 / 32 - 7 * (π ^ 4 / 90) / 32) * η ^ 4
        + (-(π ^ 2 / 6) * z 3 / 16 + 3 * z 5 / 16) * η ^ 5)| ≤ C * |η| ^ 6 ∧
    |iteratedDeriv 2 (fun s => A η s) 0 - (2 * (π ^ 2 / 6) - 2 * z 3 * η
        + (-(π ^ 2 / 6) ^ 2 / 2 + 3 * (π ^ 4 / 90) / 2) * η ^ 2
        + ((π ^ 2 / 6) * z 3 - z 5) * η ^ 3)| ≤ C * |η| ^ 4 ∧
    |iteratedDeriv 4 (fun s => A η s) 0 - (12 * (π ^ 2 / 6) ^ 2 + 12 * (π ^ 4 / 90))|
      ≤ C * |η|

/-- `thm:pw`, displayed series (main.tex 1186–1189), with `τ = 1 - 1/p` on both sides of
`p = 1`: `P/√2 = 4 + ζ(2)τ² + ζ(3)τ³ - (33/16)ζ(4)τ⁴ - (3/4 ζ(5) + 7/4 ζ(2)ζ(3))τ⁵ + O(τ⁶)`. -/
def stmt_7 : Prop :=
  ∃ C δ : ℝ, 0 < δ ∧ ∀ p : ℝ, 0 < p → |1 - 1 / p| ≤ δ →
    |P p / Real.sqrt 2 - (4 + (π ^ 2 / 6) * (1 - 1 / p) ^ 2 + z 3 * (1 - 1 / p) ^ 3
        - 33 / 16 * (π ^ 4 / 90) * (1 - 1 / p) ^ 4
        - (3 / 4 * z 5 + 7 / 4 * (π ^ 2 / 6) * z 3) * (1 - 1 / p) ^ 5)|
      ≤ C * |1 - 1 / p| ^ 6

/-- `P` is `C^n` at `p = 1` for every `n` (proof of `thm:rhombus`: "Differentiating under the
integral is justified by an integrable bound"); so `P^{(n)}(1)` are genuine derivatives. -/
def stmt_8 : Prop :=
  ∀ n : ℕ, ContDiffAt ℝ n P 1

/-- `prop:an`: `a_0 = 4`, `a_1 = 0`, `a_2 = π²/6`, `a_3 = ζ(3) - π²/3`. -/
def stmt_9 : Prop :=
  a 0 = 4 ∧ a 1 = 0 ∧ a 2 = π ^ 2 / 6 ∧ a 3 = z 3 - π ^ 2 / 3

/-- `prop:an`: `a_4 = -3ζ(3) - 11π⁴/480 + π²/2`. -/
def stmt_10 : Prop :=
  a 4 = -3 * z 3 - 11 * π ^ 4 / 480 + π ^ 2 / 2

/-- `prop:an`: `a_5 = -7π²ζ(3)/24 + 6ζ(3) - 3ζ(5)/4 - 2π²/3 + 11π⁴/120`. -/
def stmt_11 : Prop :=
  a 5 = -(7 * π ^ 2 * z 3) / 24 + 6 * z 3 - 3 * z 5 / 4 - 2 * π ^ 2 / 3 + 11 * π ^ 4 / 120

end X05
-- END STATEMENTS


open MeasureTheory Real Set

namespace X05

lemma pointwise1 {η t : ℝ} (ht0 : 0 < t) (ht1 : t < 1) :
    Real.sqrt (t ^ (2 * η) + (1 - t) ^ (2 * η)) =
      Real.sqrt 2 * ((t * (1 - t)) ^ (η / 2) * Real.sqrt (Real.cosh (η * logit t))) := by
  have h1 : 0 < 1 - t := by linarith
  have hc : 0 ≤ Real.cosh (η * logit t) := (Real.cosh_pos _).le
  have hpow : 0 ≤ (t * (1 - t)) ^ (η / 2) := by positivity
  rw [← Real.sqrt_sq hpow, ← Real.sqrt_mul (by positivity), ← Real.sqrt_mul (by norm_num)]
  congr 1
  rw [Real.cosh_eq, logit, ← Real.rpow_natCast, ← Real.rpow_mul (by positivity)]
  have e1 : Real.exp (η * Real.log (t / (1 - t))) = t ^ η * (1 - t) ^ (-η) := by
    rw [Real.log_div ht0.ne' h1.ne', mul_sub, Real.exp_sub, Real.rpow_neg h1.le,
      Real.rpow_def_of_pos ht0, Real.rpow_def_of_pos h1, mul_comm η, mul_comm η]
    field_simp
  have e2 : Real.exp (-(η * Real.log (t / (1 - t)))) = t ^ (-η) * (1 - t) ^ η := by
    rw [Real.exp_neg, e1, Real.rpow_neg ht0.le, Real.rpow_neg h1.le]
    have : 0 < t ^ η := Real.rpow_pos_of_pos ht0 _
    have : 0 < (1 - t) ^ η := Real.rpow_pos_of_pos h1 _
    field_simp
  rw [e1, e2, Real.mul_rpow ht0.le h1.le]
  push_cast
  have q1 : t ^ (η / 2 * 2) = t ^ η := by ring_nf
  have q2 : (1 - t) ^ (η / 2 * 2) = (1 - t) ^ η := by ring_nf
  rw [q1, q2]
  have a1 : t ^ (2 * η) = t ^ η * t ^ η := by rw [← Real.rpow_add ht0]; ring_nf
  have a2 : (1 - t) ^ (2 * η) = (1 - t) ^ η * (1 - t) ^ η := by rw [← Real.rpow_add h1]; ring_nf
  rw [a1, a2]
  have hA : 0 < t ^ η := Real.rpow_pos_of_pos ht0 _
  have hB : 0 < (1 - t) ^ η := Real.rpow_pos_of_pos h1 _
  rw [Real.rpow_neg ht0.le, Real.rpow_neg h1.le]
  field_simp

lemma integrand_bound {c t : ℝ} (ht0 : 0 < t) (ht1 : t < 1) :
    Real.sqrt (t ^ (2 * c) + (1 - t) ^ (2 * c)) ≤ t ^ c + (1 - t) ^ c := by
  have h1 : 0 < 1 - t := by linarith
  have ha : 0 < t ^ c := Real.rpow_pos_of_pos ht0 _
  have hb : 0 < (1 - t) ^ c := Real.rpow_pos_of_pos h1 _
  have a1 : t ^ (2 * c) = t ^ c * t ^ c := by rw [← Real.rpow_add ht0]; ring_nf
  have a2 : (1 - t) ^ (2 * c) = (1 - t) ^ c * (1 - t) ^ c := by rw [← Real.rpow_add h1]; ring_nf
  rw [a1, a2, Real.sqrt_le_iff]
  constructor
  · positivity
  · nlinarith

end X05

theorem X05.proof_1 : X05.stmt_1 := by
  intro p hp
  set η := 1 / p - 1 with hη
  have hex : 2 / p - 2 = 2 * η := by rw [hη]; ring
  have hη1 : 0 < η + 1 := by rw [hη]; simp [hp]
  have hint : IntervalIntegrable (fun t : ℝ => t ^ η + (1 - t) ^ η) volume 0 1 := by
    have ha : IntervalIntegrable (fun t : ℝ => t ^ η) volume 0 1 :=
      intervalIntegral.intervalIntegrable_rpow' (by linarith)
    have hb : IntervalIntegrable (fun t : ℝ => (1 - t) ^ η) volume 0 1 := by
      have := ha.comp_sub_left 1
      simpa using this.symm
    exact ha.add hb
  refine ⟨?_, ?_⟩
  · refine IntervalIntegrable.mono_fun' hint (by fun_prop) ?_
    rw [uIoc_of_le zero_le_one, ← restrict_Ioo_eq_restrict_Ioc]
    refine ae_restrict_of_forall_mem measurableSet_Ioo fun t ht => ?_
    dsimp only
    rw [Real.norm_eq_abs, abs_of_nonneg (Real.sqrt_nonneg _), hex]
    exact X05.integrand_bound ht.1 ht.2
  · unfold X05.P
    have h4 : (4 : ℝ) / p = 4 * (1 + η) := by rw [hη]; ring
    rw [h4]
    have : ∫ t in (0 : ℝ)..1, Real.sqrt (t ^ (2 / p - 2) + (1 - t) ^ (2 / p - 2)) =
        ∫ t in (0 : ℝ)..1, Real.sqrt 2 *
          ((t * (1 - t)) ^ (η / 2) * Real.sqrt (Real.cosh (η * X05.logit t))) := by
      refine intervalIntegral.integral_congr_uIoo fun t ht => ?_
      rw [uIoo_of_le zero_le_one] at ht
      simp only [hex]
      exact X05.pointwise1 ht.1 ht.2
    rw [this, intervalIntegral.integral_const_mul]
    ring

open MeasureTheory Real Set

namespace X05

lemma pointwise2 {η s t : ℝ} (ht0 : 0 < t) (ht1 : t < 1) :
    (t * (1 - t)) ^ (η / 2) * Real.exp (s * logit t) =
      t ^ (η / 2 + s) * (1 - t) ^ (η / 2 - s) := by
  have h1 : 0 < 1 - t := by linarith
  rw [logit, Real.log_div ht0.ne' h1.ne', Real.mul_rpow ht0.le h1.le,
    Real.rpow_add ht0, Real.rpow_sub h1, Real.rpow_def_of_pos ht0 s, Real.rpow_def_of_pos h1 s,
    mul_sub, Real.exp_sub, mul_comm s, mul_comm s]
  field_simp

/-- the real Beta integral equals the Gamma quotient -/
lemma beta_real {a b : ℝ} (ha : 0 < a) (hb : 0 < b) :
    (∫ t in (0 : ℝ)..1, t ^ (a - 1) * (1 - t) ^ (b - 1)) * Real.Gamma (a + b) =
      Real.Gamma a * Real.Gamma b := by
  have key := Complex.Gamma_mul_Gamma_eq_betaIntegral (s := (a : ℂ)) (t := (b : ℂ))
    (by simpa using ha) (by simpa using hb)
  have hB : Complex.betaIntegral a b =
      ((∫ t in (0 : ℝ)..1, t ^ (a - 1) * (1 - t) ^ (b - 1) : ℝ) : ℂ) := by
    rw [Complex.betaIntegral, ← intervalIntegral.integral_ofReal]
    refine intervalIntegral.integral_congr_uIoo fun t ht => ?_
    rw [uIoo_of_le zero_le_one] at ht
    push_cast
    rw [Complex.ofReal_cpow ht.1.le, Complex.ofReal_cpow (by linarith [ht.2] : (0:ℝ) ≤ 1 - t)]
    push_cast
    ring_nf
  rw [hB, ← Complex.ofReal_add, Complex.Gamma_ofReal, Complex.Gamma_ofReal,
    Complex.Gamma_ofReal] at key
  have k2 := key.symm
  rw [mul_comm] at k2
  exact_mod_cast k2

end X05

theorem X05.proof_2 : X05.stmt_2 := by
  intro η s hη hs
  have hu : 0 < 1 + η / 2 + s := by cases abs_lt.mp hs; linarith
  have hv : 0 < 1 + η / 2 - s := by cases abs_lt.mp hs; linarith
  have key := X05.beta_real hu hv
  have hsum : 1 + η / 2 + s + (1 + η / 2 - s) = (1 + η) + 1 := by ring
  rw [hsum, Real.Gamma_add_one (by linarith)] at key
  have hG : 0 < Real.Gamma (1 + η) := Real.Gamma_pos_of_pos (by linarith)
  unfold X05.A
  have heq : ∫ t in (0 : ℝ)..1, (t * (1 - t)) ^ (η / 2) * Real.exp (s * X05.logit t) =
      ∫ t in (0 : ℝ)..1, t ^ (1 + η / 2 + s - 1) * (1 - t) ^ (1 + η / 2 - s - 1) := by
    refine intervalIntegral.integral_congr_uIoo fun t ht => ?_
    rw [uIoo_of_le zero_le_one] at ht
    rw [X05.pointwise2 ht.1 ht.2]
    ring_nf
  rw [heq, eq_div_iff hG.ne', ← key]
  ring


open MeasureTheory Real Set

namespace X05

lemma cosh_taylor {x : ℝ} (hx : |x| ≤ 1) :
    |Real.cosh x - (1 + x ^ 2 / 2 + x ^ 4 / 24)| ≤ (x ^ 2) ^ 3 / 500 := by
  have h1 := Real.exp_bound hx (n := 6) (by norm_num)
  have hx' : |-x| ≤ 1 := by rwa [abs_neg]
  have h2 := Real.exp_bound hx' (n := 6) (by norm_num)
  simp only [Finset.sum_range_succ, Finset.sum_range_zero, Nat.factorial, abs_neg] at h1 h2
  norm_num at h1 h2
  rw [Real.cosh_eq]
  have e6 : |x| ^ 6 = (x ^ 2) ^ 3 := by rw [show |x| ^ 6 = (|x| ^ 2) ^ 3 by ring, sq_abs]
  rw [e6] at h1 h2
  have hu : 0 ≤ (x ^ 2) ^ 3 := by positivity
  rw [abs_le] at h1 h2 ⊢
  constructor <;> nlinarith [h1.1, h1.2, h2.1, h2.2]

lemma sqrt_one_add_bounds {y : ℝ} (hy0 : 0 ≤ y) (hy1 : y ≤ 1) :
    1 + y / 2 - y ^ 2 / 8 ≤ Real.sqrt (1 + y) ∧
      Real.sqrt (1 + y) ≤ 1 + y / 2 - y ^ 2 / 8 + y ^ 3 / 16 := by
  constructor
  · rw [Real.le_sqrt (by nlinarith) (by linarith)]
    nlinarith [pow_nonneg hy0 3, pow_nonneg hy0 4, mul_nonneg (pow_nonneg hy0 3) (sub_nonneg.2 hy1)]
  · rw [Real.sqrt_le_left (by nlinarith [pow_nonneg hy0 3])]
    nlinarith [pow_nonneg hy0 3, pow_nonneg hy0 4, pow_nonneg hy0 5, pow_nonneg hy0 6,
      mul_nonneg (pow_nonneg hy0 4) (sub_nonneg.2 hy1)]

lemma sqrt_le_exp_half (x : ℝ) : Real.sqrt (Real.cosh x) ≤ Real.exp (|x| / 2) := by
  rw [Real.exp_half]
  apply Real.sqrt_le_sqrt
  rw [Real.cosh_eq]
  rcases le_total 0 x with h | h
  · rw [abs_of_nonneg h]
    have : Real.exp (-x) ≤ Real.exp x := Real.exp_le_exp.2 (by linarith)
    linarith
  · rw [abs_of_nonpos h]
    have : Real.exp x ≤ Real.exp (-x) := Real.exp_le_exp.2 (by linarith)
    linarith

end X05

lemma X05.stmt4_small (x : ℝ) (hx : |x| ≤ 1) :
    |Real.sqrt (Real.cosh x) - (1 + x ^ 2 / 4 - x ^ 4 / 96)| ≤ 3 * |x| ^ 6 * Real.exp (|x| / 2) := by
  have hE : 1 ≤ Real.exp (|x| / 2) := Real.one_le_exp (by positivity)
  have e6 : |x| ^ 6 = (x ^ 2) ^ 3 := by rw [show |x| ^ 6 = (|x| ^ 2) ^ 3 by ring, sq_abs]
  have hu0 : 0 ≤ x ^ 2 := sq_nonneg x
  have hu1 : x ^ 2 ≤ 1 := by
    have := sq_abs x; nlinarith [abs_nonneg x]
  have hc := X05.cosh_taylor hx
  rw [abs_le] at hc
  set u := x ^ 2 with hu
  have x4 : x ^ 4 = u ^ 2 := by rw [hu]; ring
  rw [x4] at hc ⊢
  set y := Real.cosh x - 1 with hy
  have hcy : Real.cosh x = 1 + y := by rw [hy]; ring
  rw [hcy] at hc ⊢
  have hu3 : 0 ≤ u ^ 3 := by positivity
  have hy0 : 0 ≤ y := by nlinarith [hc.1, pow_le_one₀ hu0 hu1 (n := 2), pow_le_one₀ hu0 hu1 (n := 3)]
  have hy1 : y ≤ u := by nlinarith [hc.2, pow_le_one₀ hu0 hu1 (n := 2), pow_le_one₀ hu0 hu1 (n := 3),
    mul_le_mul_of_nonneg_left hu1 hu0, mul_le_mul_of_nonneg_left hu1 (sq_nonneg u)]
  have hyle : y ≤ 1 := le_trans hy1 hu1
  obtain ⟨b1, b2⟩ := X05.sqrt_one_add_bounds hy0 hyle
  set r := y - (u / 2 + u ^ 2 / 24) with hr
  have hyr : y = u / 2 + u ^ 2 / 24 + r := by rw [hr]; ring
  have hr1 : |r| ≤ u ^ 3 / 500 := by rw [abs_le]; constructor <;> linarith [hc.1, hc.2]
  have hy3 : y ^ 3 ≤ u ^ 3 := pow_le_pow_left₀ hy0 hy1 3
  have hD : |(1 + y / 2 - y ^ 2 / 8) - (1 + u / 4 - u ^ 2 / 96)| ≤ u ^ 3 := by
    rw [hyr]
    have hra : r ^ 2 ≤ u ^ 3 / 500 := by
      have : r ^ 2 = |r| * |r| := by rw [← sq_abs]; ring
      rw [this]
      have hr2 : |r| ≤ 1 / 500 := by
        have : u ^ 3 ≤ 1 := pow_le_one₀ hu0 hu1
        linarith
      nlinarith [abs_nonneg r]
    have hrb : |r * (u + u ^ 2 / 12)| ≤ u ^ 3 / 500 * (13 / 12) := by
      rw [abs_mul, abs_of_nonneg (by positivity : 0 ≤ u + u ^ 2 / 12)]
      have : u + u ^ 2 / 12 ≤ 13 / 12 := by nlinarith
      calc |r| * (u + u ^ 2 / 12) ≤ u ^ 3 / 500 * (u + u ^ 2 / 12) :=
            mul_le_mul_of_nonneg_right hr1 (by positivity)
        _ ≤ u ^ 3 / 500 * (13 / 12) := mul_le_mul_of_nonneg_left this (by positivity)
    have hu4 : u ^ 4 ≤ u ^ 3 := by
      have := mul_le_mul_of_nonneg_left hu1 hu3; nlinarith
    have hpoly : (1 + (u / 2 + u ^ 2 / 24 + r) / 2 - (u / 2 + u ^ 2 / 24 + r) ^ 2 / 8)
        - (1 + u / 4 - u ^ 2 / 96) =
        r / 2 - (u ^ 3 / 24 + u ^ 4 / 576 + r * (u + u ^ 2 / 12) + r ^ 2) / 8 := by ring
    rw [hpoly]
    rw [abs_le] at hr1 hrb ⊢
    constructor <;> nlinarith [sq_nonneg r]
  have hfin : |Real.sqrt (1 + y) - (1 + u / 4 - u ^ 2 / 96)| ≤ 2 * u ^ 3 := by
    rw [abs_le] at hD ⊢
    constructor <;> nlinarith
  rw [e6]
  have a2 := mul_le_mul_of_nonneg_left hE (by positivity : (0:ℝ) ≤ 3 * u ^ 3)
  linarith

lemma X05.stmt4_large (x : ℝ) (hx : 1 < |x|) :
    |Real.sqrt (Real.cosh x) - (1 + x ^ 2 / 4 - x ^ 4 / 96)| ≤ 3 * |x| ^ 6 * Real.exp (|x| / 2) := by
  have hE : 1 ≤ Real.exp (|x| / 2) := Real.one_le_exp (by positivity)
  have e6 : |x| ^ 6 = (x ^ 2) ^ 3 := by rw [show |x| ^ 6 = (|x| ^ 2) ^ 3 by ring, sq_abs]
  have hu0 : 0 ≤ x ^ 2 := sq_nonneg x
  have hs := X05.sqrt_le_exp_half x
  have hs0 := Real.sqrt_nonneg (Real.cosh x)
  have hx6 : 1 ≤ |x| ^ 6 := one_le_pow₀ hx.le
  have hx2 : x ^ 2 ≤ |x| ^ 6 := by
    rw [← sq_abs]; exact pow_le_pow_right₀ hx.le (by norm_num)
  have hx4 : x ^ 4 ≤ |x| ^ 6 := by
    have h44 : |x| ^ 4 = x ^ 4 := by
      rw [show |x| ^ 4 = (|x| ^ 2) ^ 2 by ring, sq_abs]; ring
    rw [← h44]
    exact pow_le_pow_right₀ hx.le (by norm_num)
  have h4 : 0 ≤ x ^ 4 := by positivity
  set Q := 1 + x ^ 2 / 4 - x ^ 4 / 96 with hQ
  have hP : |Q| ≤ 2 * |x| ^ 6 := by
    rw [abs_le, hQ]; constructor <;> linarith only [h4, hx2, hx4, hx6, hu0]
  have t1 : |Real.sqrt (Real.cosh x) - Q| ≤ Real.sqrt (Real.cosh x) + |Q| := by
    have := abs_sub (Real.sqrt (Real.cosh x)) Q
    rwa [abs_of_nonneg hs0] at this
  clear_value Q
  generalize Real.sqrt (Real.cosh x) = S at t1 hs hs0 ⊢
  generalize Real.exp (|x| / 2) = E at hs hE ⊢
  generalize |x| ^ 6 = W at hP hx6 ⊢
  have a1 := mul_le_mul_of_nonneg_right hx6 (le_trans zero_le_one hE)
  have a2 := mul_le_mul_of_nonneg_left hE (le_trans zero_le_one hx6)
  rw [one_mul] at a1
  rw [mul_one] at a2
  have : 3 * W * E = W * E + 2 * (W * E) := by ring
  rw [this]
  linarith only [a1, a2, t1, hs, hP]

theorem X05.proof_4 : X05.stmt_4 := by
  refine ⟨3, fun x => ?_⟩
  rcases le_or_gt |x| 1 with hx | hx
  · exact X05.stmt4_small x hx
  · exact X05.stmt4_large x hx


open Filter Topology Finset

namespace X05

lemma z_hasSum {k : ℕ} (hk : 1 < k) :
    HasSum (fun n : ℕ => 1 / ((n : ℝ) + 1) ^ k) (z k) := by
  have hs : Summable (fun n : ℕ => 1 / (n : ℝ) ^ k) := Real.summable_one_div_nat_pow.mpr hk
  have hz : z k = ∑' n : ℕ, 1 / (n : ℝ) ^ k := by
    unfold z
    rw [zeta_nat_eq_tsum_of_gt_one hk]
    rw [show (fun n : ℕ => 1 / (n : ℂ) ^ k) = fun n : ℕ => ((1 / (n : ℝ) ^ k : ℝ) : ℂ) by
      ext n; push_cast; rfl]
    rw [← Complex.ofReal_tsum, Complex.ofReal_re]
  rw [hz]
  have := hs.hasSum
  rw [← hasSum_nat_add_iff' 1] at this
  simpa [zero_pow (show k ≠ 0 by omega)] using this

lemma prod_shift (x : ℝ) (n : ℕ) :
    ∏ j ∈ range (n + 1), (1 + x + (j : ℝ)) =
      ((n + 1).factorial : ℝ) * ∏ k ∈ range (n + 1), (1 + x / ((k : ℝ) + 1)) := by
  rw [show (((n + 1).factorial : ℕ) : ℝ) = ∏ k ∈ range (n + 1), ((k : ℝ) + 1) by
    rw [← Finset.prod_range_add_one_eq_factorial]; push_cast; rfl, ← Finset.prod_mul_distrib]
  refine Finset.prod_congr rfl fun k _ => ?_
  have : (0 : ℝ) < (k : ℝ) + 1 := by positivity
  field_simp
  ring

lemma fac_pos {x : ℝ} (hx : |x| < 1) (k : ℕ) : 0 < 1 + x / ((k : ℝ) + 1) := by
  have h1 : (1 : ℝ) ≤ (k : ℝ) + 1 := by have := k.cast_nonneg (α := ℝ); linarith
  have : |x / ((k : ℝ) + 1)| < 1 := by
    rw [abs_div, abs_of_pos (by linarith : (0:ℝ) < (k : ℝ) + 1), div_lt_one (by linarith)]
    linarith
  linarith [neg_abs_le (x / ((k : ℝ) + 1))]

/-- the partial products `R_m` -/
noncomputable def R (η u v : ℝ) (m : ℕ) : ℝ :=
  ∏ k ∈ range m, ((1 + η / ((k : ℝ) + 1)) / ((1 + u / ((k : ℝ) + 1)) * (1 + v / ((k : ℝ) + 1))))

lemma gammaSeq_ratio {η u v : ℝ} (huv : u + v = η) (hu : |u| < 1) (hv : |v| < 1)
    (hη : |η| < 1) {n : ℕ} (hn : 1 ≤ n) :
    Real.GammaSeq (1 + u) n * Real.GammaSeq (1 + v) n / Real.GammaSeq (1 + η) n =
      (n : ℝ) / ((n : ℝ) + 1) * R η u v (n + 1) := by
  unfold Real.GammaSeq R
  rw [prod_shift, prod_shift, prod_shift, Finset.prod_div_distrib, Finset.prod_mul_distrib]
  have hn0 : (0 : ℝ) < n := by exact_mod_cast hn
  have pu : 0 < ∏ k ∈ range (n + 1), (1 + u / ((k : ℝ) + 1)) :=
    Finset.prod_pos fun k _ => fac_pos hu k
  have pv : 0 < ∏ k ∈ range (n + 1), (1 + v / ((k : ℝ) + 1)) :=
    Finset.prod_pos fun k _ => fac_pos hv k
  have pη : 0 < ∏ k ∈ range (n + 1), (1 + η / ((k : ℝ) + 1)) :=
    Finset.prod_pos fun k _ => fac_pos hη k
  have e : ∀ x : ℝ, (n : ℝ) ^ (1 + x) = n * (n : ℝ) ^ x := fun x => by
    rw [Real.rpow_add hn0, Real.rpow_one]
  have euv : (n : ℝ) ^ u * (n : ℝ) ^ v = (n : ℝ) ^ η := by rw [← Real.rpow_add hn0, huv]
  rw [e, e, e]
  have hpη : 0 < (n : ℝ) ^ η := Real.rpow_pos_of_pos hn0 _
  have hf : (0 : ℝ) < (n.factorial : ℝ) := by exact_mod_cast n.factorial_pos
  have hf1 : (0 : ℝ) < ((n + 1).factorial : ℝ) := by exact_mod_cast (n + 1).factorial_pos
  rw [Nat.factorial_succ] at hf1 ⊢
  push_cast at hf1 ⊢
  field_simp
  rw [← euv]

lemma log_series_k {η u v : ℝ} (huv : u + v = η) (hu : |u| < 1) (hv : |v| < 1)
    (hη : |η| < 1) (k : ℕ) :
    HasSum (fun j : ℕ => (-1 : ℝ) ^ (j + 2) * (1 / ((k : ℝ) + 1) ^ (j + 2)) / ((j + 2 : ℕ) : ℝ) *
        (u ^ (j + 2) + v ^ (j + 2) - η ^ (j + 2)))
      (Real.log (1 + η / ((k : ℝ) + 1)) - Real.log (1 + u / ((k : ℝ) + 1))
        - Real.log (1 + v / ((k : ℝ) + 1))) := by
  have hk : (0 : ℝ) < (k : ℝ) + 1 := by positivity
  have sm : ∀ x : ℝ, |x| < 1 → |-(x / ((k : ℝ) + 1))| < 1 := by
    intro x hx
    have h1 : (1 : ℝ) ≤ (k : ℝ) + 1 := by have := k.cast_nonneg (α := ℝ); linarith
    rw [abs_neg, abs_div, abs_of_pos hk, div_lt_one hk]; linarith
  have Lu := Real.hasSum_pow_div_log_of_abs_lt_one (sm u hu)
  have Lv := Real.hasSum_pow_div_log_of_abs_lt_one (sm v hv)
  have Lη := Real.hasSum_pow_div_log_of_abs_lt_one (sm η hη)
  have L := (Lu.add Lv).sub Lη
  simp only [sub_neg_eq_add] at L
  rw [← hasSum_nat_add_iff' 1] at L
  simp only [Finset.range_one, Finset.sum_singleton] at L
  convert L using 1
  · ext j
    simp only [neg_pow (_ / _), div_pow]
    push_cast
    ring
  · rw [← huv]
    ring_nf

lemma log_R {η u v : ℝ} (hu : |u| < 1) (hv : |v| < 1) (hη : |η| < 1) (m : ℕ) :
    Real.log (R η u v m) = ∑ k ∈ range m, (Real.log (1 + η / ((k : ℝ) + 1))
      - Real.log (1 + u / ((k : ℝ) + 1)) - Real.log (1 + v / ((k : ℝ) + 1))) := by
  unfold R
  rw [Real.log_prod]
  · refine Finset.sum_congr rfl fun k _ => ?_
    have a := fac_pos hu k; have b := fac_pos hv k; have c := fac_pos hη k
    rw [Real.log_div c.ne' (by positivity), Real.log_mul a.ne' b.ne']
    ring
  · intro k _
    have a := fac_pos hu k; have b := fac_pos hv k; have c := fac_pos hη k
    positivity

lemma R_tendsto {η u v : ℝ} (huv : u + v = η) (hu : |u| < 1) (hv : |v| < 1)
    (hη : |η| < 1) (hη1 : -1 < η) :
    Tendsto (fun n : ℕ => R η u v (n + 1)) atTop
      (𝓝 (Real.Gamma (1 + u) * Real.Gamma (1 + v) / Real.Gamma (1 + η))) := by
  have hG : Real.Gamma (1 + η) ≠ 0 := (Real.Gamma_pos_of_pos (by linarith)).ne'
  have h1 := ((Real.GammaSeq_tendsto_Gamma (1 + u)).mul
    (Real.GammaSeq_tendsto_Gamma (1 + v))).div (Real.GammaSeq_tendsto_Gamma (1 + η)) hG
  have h2 : Tendsto (fun n : ℕ => ((n : ℝ) + 1) / n) atTop (𝓝 1) := by
    have := (tendsto_one_div_atTop_nhds_zero_nat).const_add (1 : ℝ)
    rw [add_zero] at this
    refine this.congr' ?_
    filter_upwards [eventually_ge_atTop 1] with n hn
    have : (0 : ℝ) < n := by exact_mod_cast hn
    field_simp
  have h3 := h1.mul h2
  rw [mul_one] at h3
  refine h3.congr' ?_
  filter_upwards [eventually_ge_atTop 1] with n hn
  have : (0 : ℝ) < n := by exact_mod_cast hn
  simp only [Pi.div_apply]
  rw [gammaSeq_ratio huv hu hv hη hn]
  field_simp

end X05

theorem X05.proof_3 : X05.stmt_3 := by
  intro η s hη hs
  set u := η / 2 + s with hudef
  set v := η / 2 - s with hvdef
  have huv : u + v = η := by rw [hudef, hvdef]; ring
  obtain ⟨hη1, hη2⟩ := abs_lt.mp hη
  obtain ⟨hs1, hs2⟩ := abs_lt.mp hs
  have hu : |u| < 3 / 4 := by rw [abs_lt]; constructor <;> linarith
  have hv : |v| < 3 / 4 := by rw [abs_lt]; constructor <;> linarith
  have hu' : |u| < 1 := by linarith
  have hv' : |v| < 1 := by linarith
  have hη' : |η| < 1 := by linarith
  have hA : X05.A η s = Real.Gamma (1 + u) * Real.Gamma (1 + v) / Real.Gamma (1 + η) := by
    rw [X05.proof_2 η s (by linarith) (by rw [abs_lt]; constructor <;> linarith)]
    rw [hudef, hvdef, ← add_assoc, ← add_sub_assoc]
  have hApos : 0 < X05.A η s := by
    rw [hA]
    have := Real.Gamma_pos_of_pos (show 0 < 1 + u by cases abs_lt.mp hu; linarith)
    have := Real.Gamma_pos_of_pos (show 0 < 1 + v by cases abs_lt.mp hv; linarith)
    have := Real.Gamma_pos_of_pos (show 0 < 1 + η by linarith)
    positivity
  refine ⟨hApos, ?_⟩
  -- the double series
  set f : ℕ × ℕ → ℝ := fun p => (-1 : ℝ) ^ (p.2 + 2) * (1 / ((p.1 : ℝ) + 1) ^ (p.2 + 2)) /
      ((p.2 + 2 : ℕ) : ℝ) * (u ^ (p.2 + 2) + v ^ (p.2 + 2) - η ^ (p.2 + 2)) with hfdef
  have hfs : Summable f := by
    have hg1 : Summable (fun k : ℕ => 1 / ((k : ℝ) + 1) ^ 2) := (X05.z_hasSum (by norm_num)).summable
    have hg2 : Summable (fun j : ℕ => 3 * (3 / 4 : ℝ) ^ j) :=
      (summable_geometric_of_lt_one (by norm_num) (by norm_num)).mul_left 3
    refine Summable.of_norm_bounded (hg1.mul_of_nonneg hg2 (fun k => by positivity)
      (fun j => by positivity)) ?_
    rintro ⟨k, j⟩
    simp only [hfdef, Real.norm_eq_abs]
    have hk1 : (1 : ℝ) ≤ (k : ℝ) + 1 := by have := k.cast_nonneg (α := ℝ); linarith
    have hj : (1 : ℝ) ≤ ((j + 2 : ℕ) : ℝ) := by push_cast; have := j.cast_nonneg (α := ℝ); linarith
    rw [abs_mul, abs_div, abs_mul, abs_pow, abs_neg, abs_one, one_pow, one_mul,
      abs_of_pos (by positivity : (0 : ℝ) < 1 / ((k : ℝ) + 1) ^ (j + 2)),
      abs_of_pos (by positivity : (0 : ℝ) < ((j + 2 : ℕ) : ℝ))]
    have b1 : 1 / ((k : ℝ) + 1) ^ (j + 2) ≤ 1 / ((k : ℝ) + 1) ^ 2 := by
      apply one_div_le_one_div_of_le (by positivity)
      exact pow_le_pow_right₀ hk1 (by omega)
    have b2 : |u ^ (j + 2) + v ^ (j + 2) - η ^ (j + 2)| ≤ 3 * (3 / 4 : ℝ) ^ j := by
      have pu : |u| ^ (j + 2) ≤ (3 / 4 : ℝ) ^ j :=
        (pow_le_pow_left₀ (abs_nonneg u) hu.le _).trans
          (pow_le_pow_of_le_one (by norm_num) (by norm_num) (by omega))
      have pv : |v| ^ (j + 2) ≤ (3 / 4 : ℝ) ^ j :=
        (pow_le_pow_left₀ (abs_nonneg v) hv.le _).trans
          (pow_le_pow_of_le_one (by norm_num) (by norm_num) (by omega))
      have pη : |η| ^ (j + 2) ≤ (3 / 4 : ℝ) ^ j :=
        (pow_le_pow_left₀ (abs_nonneg η) (show |η| ≤ 3 / 4 by linarith) _).trans
          (pow_le_pow_of_le_one (by norm_num) (by norm_num) (by omega))
      calc |u ^ (j + 2) + v ^ (j + 2) - η ^ (j + 2)|
          ≤ |u ^ (j + 2)| + |v ^ (j + 2)| + |η ^ (j + 2)| := by
            have := abs_sub (u ^ (j + 2) + v ^ (j + 2)) (η ^ (j + 2))
            have := abs_add_le (u ^ (j + 2)) (v ^ (j + 2))
            linarith
        _ ≤ 3 * (3 / 4 : ℝ) ^ j := by rw [abs_pow, abs_pow, abs_pow]; linarith
    calc 1 / ((k : ℝ) + 1) ^ (j + 2) / ((j + 2 : ℕ) : ℝ) * |u ^ (j + 2) + v ^ (j + 2) - η ^ (j + 2)|
        ≤ 1 / ((k : ℝ) + 1) ^ 2 / 1 * (3 * (3 / 4 : ℝ) ^ j) := by
          gcongr
      _ = _ := by simp
  have hf := hfs.hasSum
  -- fibers in k
  have hT := hf.prod_fiberwise (fun k => X05.log_series_k huv hu' hv' hη' k)
  -- fibers in j
  have hf' : HasSum (f ∘ Equiv.prodComm ℕ ℕ) (∑' p, f p) := (Equiv.prodComm ℕ ℕ).hasSum_iff.mpr hf
  have hZ := hf'.prod_fiberwise (fun j => by
    have := ((X05.z_hasSum (show 1 < j + 2 by omega)).mul_left ((-1 : ℝ) ^ (j + 2))).div_const
      ((j + 2 : ℕ) : ℝ) |>.mul_right (u ^ (j + 2) + v ^ (j + 2) - η ^ (j + 2))
    exact this)
  -- identify the sum
  have hlim1 := (tendsto_add_atTop_iff_nat 1).mpr hT.tendsto_sum_nat
  have hlim2 : Tendsto (fun n : ℕ => ∑ k ∈ range (n + 1), (Real.log (1 + η / ((k : ℝ) + 1))
      - Real.log (1 + u / ((k : ℝ) + 1)) - Real.log (1 + v / ((k : ℝ) + 1)))) atTop
      (𝓝 (Real.log (X05.A η s))) := by
    have := ((Real.continuousAt_log hApos.ne').tendsto).comp
      (hA ▸ X05.R_tendsto huv hu' hv' hη' (by linarith))
    refine this.congr fun n => ?_
    simp only [Function.comp]
    exact X05.log_R hu' hv' hη' (n + 1)
  have heq := tendsto_nhds_unique hlim1 hlim2
  rw [heq] at hZ
  convert hZ using 1


namespace X05

lemma ae_ne_one : ∀ᵐ t ∂(volume : Measure ℝ), t ≠ 1 := by
  rw [ae_iff]; simp

lemma log_bd {x : ℝ} (hx0 : 0 < x) (hx1 : x ≤ 1) :
    |Real.log x| ≤ 48 * x ^ (-(1 / 48 : ℝ)) := by
  have h := Real.neg_rpow_div_le_log hx0.le (show (0 : ℝ) < 1 / 48 by norm_num)
  have hl : Real.log x ≤ 0 := Real.log_nonpos hx0.le hx1
  rw [abs_of_nonpos hl]
  have : -x ^ (-(1 / 48 : ℝ)) / (1 / 48) = -(48 * x ^ (-(1 / 48 : ℝ))) := by ring
  linarith

lemma logit_bd {t : ℝ} (ht0 : 0 < t) (ht1 : t < 1) :
    |logit t| ≤ 96 * (t * (1 - t)) ^ (-(1 / 48 : ℝ)) := by
  have h1 : 0 < 1 - t := by linarith
  have hy : 0 < t * (1 - t) := mul_pos ht0 h1
  rw [logit, Real.log_div ht0.ne' h1.ne']
  have a := log_bd ht0 ht1.le
  have b := log_bd h1 (by linarith)
  have a' : t ^ (-(1 / 48 : ℝ)) ≤ (t * (1 - t)) ^ (-(1 / 48 : ℝ)) :=
    Real.rpow_le_rpow_of_nonpos hy (by nlinarith) (by norm_num)
  have b' : (1 - t) ^ (-(1 / 48 : ℝ)) ≤ (t * (1 - t)) ^ (-(1 / 48 : ℝ)) :=
    Real.rpow_le_rpow_of_nonpos hy (by nlinarith) (by norm_num)
  calc |Real.log t - Real.log (1 - t)| ≤ |Real.log t| + |Real.log (1 - t)| := abs_sub _ _
    _ ≤ _ := by linarith

lemma logit_pow_bd {t : ℝ} (ht0 : 0 < t) (ht1 : t < 1) {m : ℕ} (hm : m ≤ 6) :
    |logit t| ^ m ≤ 96 ^ 6 * (t * (1 - t)) ^ (-(1 / 8 : ℝ)) := by
  have h1 : 0 < 1 - t := by linarith
  have hy : 0 < t * (1 - t) := mul_pos ht0 h1
  have hy1 : t * (1 - t) ≤ 1 := by nlinarith
  set X := (t * (1 - t)) ^ (-(1 / 48 : ℝ)) with hX
  have hX1 : 1 ≤ X := Real.one_le_rpow_of_pos_of_le_one_of_nonpos hy hy1 (by norm_num)
  have hb := logit_bd ht0 ht1
  have h96 : 1 ≤ 96 * X := by linarith
  calc |logit t| ^ m ≤ (96 * X) ^ m := pow_le_pow_left₀ (abs_nonneg _) hb m
    _ ≤ (96 * X) ^ 6 := pow_le_pow_right₀ h96 hm
    _ = 96 ^ 6 * (t * (1 - t)) ^ (-(1 / 8 : ℝ)) := by
      rw [mul_pow, hX, ← Real.rpow_mul_natCast hy.le]
      norm_num

lemma pow_pair_bd {t a b : ℝ} (ht0 : 0 < t) (ht1 : t < 1) (ha : |a| ≤ 1 / 4) (hb : |b| ≤ 1 / 4) :
    t ^ a * (1 - t) ^ b ≤ (t * (1 - t)) ^ (-(1 / 2 : ℝ)) := by
  have h1 : 0 < 1 - t := by linarith
  have hy : 0 < t * (1 - t) := mul_pos ht0 h1
  have ta : t ^ a ≤ (t * (1 - t)) ^ (-(1 / 4 : ℝ)) :=
    (Real.rpow_le_rpow_of_exponent_ge ht0 ht1.le (by linarith [neg_abs_le a])).trans
      (Real.rpow_le_rpow_of_nonpos hy (by nlinarith) (by norm_num))
  have tb : (1 - t) ^ b ≤ (t * (1 - t)) ^ (-(1 / 4 : ℝ)) :=
    (Real.rpow_le_rpow_of_exponent_ge h1 (by linarith) (by linarith [neg_abs_le b])).trans
      (Real.rpow_le_rpow_of_nonpos hy (by nlinarith) (by norm_num))
  calc t ^ a * (1 - t) ^ b ≤ (t * (1 - t)) ^ (-(1 / 4 : ℝ)) * (t * (1 - t)) ^ (-(1 / 4 : ℝ)) :=
        mul_le_mul ta tb (by positivity) (by positivity)
    _ = _ := by rw [← Real.rpow_add hy]; norm_num

/-- the dominating function `(t(1-t))^{-5/8}` -/
noncomputable def q (t : ℝ) : ℝ := (t * (1 - t)) ^ (-(5 / 8 : ℝ))

lemma main_bd {t a b : ℝ} (ht0 : 0 < t) (ht1 : t < 1) (ha : |a| ≤ 1 / 4) (hb : |b| ≤ 1 / 4)
    {m : ℕ} (hm : m ≤ 6) :
    |logit t| ^ m * (t ^ a * (1 - t) ^ b) ≤ 96 ^ 6 * q t := by
  have h1 : 0 < 1 - t := by linarith
  have hy : 0 < t * (1 - t) := mul_pos ht0 h1
  have e1 := logit_pow_bd ht0 ht1 hm
  have e2 := pow_pair_bd ht0 ht1 ha hb
  calc |logit t| ^ m * (t ^ a * (1 - t) ^ b)
      ≤ (96 ^ 6 * (t * (1 - t)) ^ (-(1 / 8 : ℝ))) * (t * (1 - t)) ^ (-(1 / 2 : ℝ)) :=
        mul_le_mul e1 e2 (by positivity) (by positivity)
    _ = 96 ^ 6 * q t := by
      rw [q, mul_assoc, ← Real.rpow_add hy]; norm_num

lemma q_int : IntervalIntegrable q volume 0 1 := by
  have ha : IntervalIntegrable (fun t : ℝ => t ^ (-(5 / 8 : ℝ))) volume 0 1 :=
    intervalIntegral.intervalIntegrable_rpow' (by norm_num)
  have hb : IntervalIntegrable (fun t : ℝ => (1 - t) ^ (-(5 / 8 : ℝ))) volume 0 1 := by
    have := ha.comp_sub_left 1
    simpa using this.symm
  have hab := (ha.add hb).const_mul 2
  refine IntervalIntegrable.mono_fun' hab
    (by unfold q; exact (by fun_prop : Measurable fun t : ℝ => (t * (1 - t)) ^ (-(5 / 8 : ℝ))).aestronglyMeasurable) ?_
  rw [uIoc_of_le zero_le_one, ← restrict_Ioo_eq_restrict_Ioc]
  refine ae_restrict_of_forall_mem measurableSet_Ioo fun t ht => ?_
  have ht0 := ht.1; have ht1 := ht.2
  have h1 : 0 < 1 - t := by linarith
  have hq : 0 ≤ q t := by unfold q; positivity
  dsimp only
  rw [Real.norm_eq_abs, abs_of_nonneg hq, q, Real.mul_rpow ht0.le h1.le]
  have hA : 0 < t ^ (-(5 / 8 : ℝ)) := Real.rpow_pos_of_pos ht0 _
  have hB : 0 < (1 - t) ^ (-(5 / 8 : ℝ)) := Real.rpow_pos_of_pos h1 _
  have h2 : (1 / 2 : ℝ) ^ (-(5 / 8 : ℝ)) ≤ 2 := by
    have := Real.rpow_le_rpow_of_exponent_ge (show (0:ℝ) < 1 / 2 by norm_num)
      (show (1 / 2 : ℝ) ≤ 1 by norm_num) (show (-1 : ℝ) ≤ -(5 / 8) by norm_num)
    rw [Real.rpow_neg_one] at this; norm_num at this ⊢; linarith
  rcases le_total t (1 / 2) with h | h
  · have : (1 - t) ^ (-(5 / 8 : ℝ)) ≤ 2 :=
      (Real.rpow_le_rpow_of_nonpos (by norm_num) (by linarith) (by norm_num)).trans h2
    nlinarith
  · have : t ^ (-(5 / 8 : ℝ)) ≤ 2 :=
      (Real.rpow_le_rpow_of_nonpos (by norm_num) h (by norm_num)).trans h2
    nlinarith

/-- `G_k(s) = ∫ (t(1-t))^{η/2} l^k e^{s l} dt` -/
noncomputable def G (η : ℝ) (k : ℕ) (s : ℝ) : ℝ :=
  ∫ t in (0 : ℝ)..1, (t * (1 - t)) ^ (η / 2) * logit t ^ k * Real.exp (s * logit t)

lemma G_integrand_bd {η s t : ℝ} (hη : |η| ≤ 1 / 4) (hs : |s| ≤ 1 / 8) (ht0 : 0 < t) (ht1 : t < 1)
    {k : ℕ} (hk : k ≤ 6) :
    ‖(t * (1 - t)) ^ (η / 2) * logit t ^ k * Real.exp (s * logit t)‖ ≤ 96 ^ 6 * q t := by
  have e : (t * (1 - t)) ^ (η / 2) * logit t ^ k * Real.exp (s * logit t) =
      logit t ^ k * ((t * (1 - t)) ^ (η / 2) * Real.exp (s * logit t)) := by ring
  rw [e, pointwise2 ht0 ht1, Real.norm_eq_abs, abs_mul, abs_pow,
    abs_of_pos (by have : 0 < 1 - t := by linarith
                   positivity : (0:ℝ) < t ^ (η / 2 + s) * (1 - t) ^ (η / 2 - s))]
  have h1 : |η / 2 + s| ≤ 1 / 4 := by
    rw [abs_le] at hη hs ⊢; constructor <;> linarith
  have h2 : |η / 2 - s| ≤ 1 / 4 := by
    rw [abs_le] at hη hs ⊢; constructor <;> linarith
  exact main_bd ht0 ht1 h1 h2 hk

lemma G_int {η s : ℝ} (hη : |η| ≤ 1 / 4) (hs : |s| ≤ 1 / 8) {k : ℕ} (hk : k ≤ 6) :
    IntervalIntegrable (fun t => (t * (1 - t)) ^ (η / 2) * logit t ^ k * Real.exp (s * logit t))
      volume 0 1 := by
  refine IntervalIntegrable.mono_fun' (q_int.const_mul (96 ^ 6)) (by unfold logit; fun_prop) ?_
  rw [uIoc_of_le zero_le_one, ← restrict_Ioo_eq_restrict_Ioc]
  exact ae_restrict_of_forall_mem measurableSet_Ioo fun t ht =>
    G_integrand_bd hη hs ht.1 ht.2 hk

lemma G_deriv {η s : ℝ} (hη : |η| ≤ 1 / 4) (hs : |s| < 1 / 8) {k : ℕ} (hk : k ≤ 5) :
    HasDerivAt (G η k) (G η (k + 1) s) s := by
  have hS : Set.Ioo (-(1 / 8) : ℝ) (1 / 8) ∈ 𝓝 s :=
    Ioo_mem_nhds (by linarith [neg_abs_le s]) (by linarith [le_abs_self s])
  have := intervalIntegral.hasDerivAt_integral_of_dominated_loc_of_deriv_le
    (μ := volume) (a := 0) (b := 1)
    (F := fun x t => (t * (1 - t)) ^ (η / 2) * logit t ^ k * Real.exp (x * logit t))
    (F' := fun x t => (t * (1 - t)) ^ (η / 2) * logit t ^ (k + 1) * Real.exp (x * logit t))
    (bound := fun t => 96 ^ 6 * q t) hS
    (Filter.Eventually.of_forall fun x => by unfold logit; fun_prop)
    (G_int hη hs.le (by omega)) (by unfold logit; fun_prop) ?_ (q_int.const_mul _) ?_
  · exact this.2
  · filter_upwards [ae_ne_one] with t ht1 htm x hx
    rw [uIoc_of_le zero_le_one] at htm
    have hx' : |x| ≤ 1 / 8 := by rw [abs_le]; exact ⟨hx.1.le, hx.2.le⟩
    exact G_integrand_bd hη hx' htm.1 (lt_of_le_of_ne htm.2 ht1) (by omega)
  · refine Filter.Eventually.of_forall fun t _ x _ => ?_
    have h := ((hasDerivAt_id' x).mul_const (logit t)).exp.const_mul
      ((t * (1 - t)) ^ (η / 2) * logit t ^ k)
    convert h using 1
    ring

lemma A_eq_G (η : ℝ) : (fun s => A η s) = fun s => (1 + η) * G η 0 s := by
  ext s; simp [A, G]

lemma iter_A {η : ℝ} (hη : |η| ≤ 1 / 4) :
    ∀ k ≤ 5, ∀ s, |s| < 1 / 8 → iteratedDeriv k (fun s => A η s) s = (1 + η) * G η k s := by
  intro k
  induction k with
  | zero => intro _ s _; simp [A_eq_G]
  | succ k ih =>
    intro hk s hs
    rw [iteratedDeriv_succ]
    have hev : iteratedDeriv k (fun s => A η s) =ᶠ[𝓝 s] fun s => (1 + η) * G η k s := by
      have hS : Set.Ioo (-(1 / 8) : ℝ) (1 / 8) ∈ 𝓝 s :=
        Ioo_mem_nhds (by linarith [neg_abs_le s]) (by linarith [le_abs_self s])
      filter_upwards [hS] with x hx
      exact ih (by omega) x (by rw [abs_lt]; exact hx)
    rw [hev.deriv_eq]
    exact ((G_deriv hη hs (by omega)).const_mul (1 + η)).deriv

lemma sqrtcosh_bd {η t : ℝ} (hη : |η| ≤ 1 / 4) (ht0 : 0 < t) (ht1 : t < 1) {k : ℕ} (hk : k ≤ 6) :
    |logit t| ^ k * ((t * (1 - t)) ^ (η / 2) * Real.exp (|η * logit t| / 2)) ≤
      2 * (96 ^ 6 * q t) := by
  have hη2 : |η / 2| ≤ 1 / 8 := by rw [abs_div]; norm_num; linarith
  have hη2' : |-(η / 2)| ≤ 1 / 8 := by rwa [abs_neg]
  have a := G_integrand_bd hη hη2 ht0 ht1 hk
  have b := G_integrand_bd hη hη2' ht0 ht1 hk
  have hw : 0 ≤ (t * (1 - t)) ^ (η / 2) := by
    have : 0 < 1 - t := by linarith
    positivity
  rw [Real.norm_eq_abs, abs_mul, abs_mul, abs_pow, abs_of_nonneg hw,
    abs_of_pos (Real.exp_pos _)] at a b
  have he : Real.exp (|η * logit t| / 2) ≤
      Real.exp (η / 2 * logit t) + Real.exp (-(η / 2) * logit t) := by
    rcases le_total 0 (η * logit t) with h | h
    · rw [abs_of_nonneg h]
      have : Real.exp (η * logit t / 2) = Real.exp (η / 2 * logit t) := by ring_nf
      rw [this]; linarith [Real.exp_pos (-(η / 2) * logit t)]
    · rw [abs_of_nonpos h]
      have : Real.exp (-(η * logit t) / 2) = Real.exp (-(η / 2) * logit t) := by ring_nf
      rw [this]; linarith [Real.exp_pos (η / 2 * logit t)]
  have hl : 0 ≤ |logit t| ^ k := by positivity
  calc |logit t| ^ k * ((t * (1 - t)) ^ (η / 2) * Real.exp (|η * logit t| / 2))
      ≤ |logit t| ^ k * ((t * (1 - t)) ^ (η / 2) *
          (Real.exp (η / 2 * logit t) + Real.exp (-(η / 2) * logit t))) :=
        mul_le_mul_of_nonneg_left (mul_le_mul_of_nonneg_left he hw) hl
    _ = (t * (1 - t)) ^ (η / 2) * |logit t| ^ k * Real.exp (η / 2 * logit t) +
        (t * (1 - t)) ^ (η / 2) * |logit t| ^ k * Real.exp (-(η / 2) * logit t) := by ring
    _ ≤ _ := by linarith

end X05

theorem X05.proof_5 : X05.stmt_5 := by
  obtain ⟨C4, hC4⟩ := X05.proof_4
  have hC4' : 0 ≤ C4 := by
    have h := hC4 1
    have h0 := abs_nonneg (Real.sqrt (Real.cosh 1) - (1 + 1 ^ 2 / 4 - 1 ^ 4 / 96))
    have he := Real.exp_pos (|(1:ℝ)| / 2)
    simp only [abs_one, one_pow, mul_one] at h he
    by_contra hneg
    push Not at hneg
    nlinarith
  set I := ∫ t in (0 : ℝ)..1, X05.q t with hI
  refine ⟨5 / 4 * (C4 * (2 * 96 ^ 6) * I), fun η hη => ?_⟩
  have hη0 : 0 < 1 + η := by rw [abs_le] at hη; linarith
  have hη5 : 1 + η ≤ 5 / 4 := by rw [abs_le] at hη; linarith
  have hp : 0 < 1 / (1 + η) := by positivity
  have h1 := (X05.proof_1 (1 / (1 + η)) hp).2
  have e : 1 / (1 / (1 + η)) - 1 = η := by field_simp; ring
  rw [e] at h1
  have hA0 : X05.A η 0 = (1 + η) * X05.G η 0 0 := congrFun (X05.A_eq_G η) 0
  have h2 := X05.iter_A hη 2 (by norm_num) 0 (by norm_num)
  have h4 := X05.iter_A hη 4 (by norm_num) 0 (by norm_num)
  rw [h1, hA0, h2, h4]
  -- integrability
  have i0 := X05.G_int hη (s := 0) (by norm_num) (k := 0) (by norm_num)
  have i2 := X05.G_int hη (s := 0) (by norm_num) (k := 2) (by norm_num)
  have i4 := X05.G_int hη (s := 0) (by norm_num) (k := 4) (by norm_num)
  have iJ : IntervalIntegrable (fun t => (t * (1 - t)) ^ (η / 2) *
      Real.sqrt (Real.cosh (η * X05.logit t))) volume 0 1 := by
    refine IntervalIntegrable.mono_fun' ((X05.q_int.const_mul (96 ^ 6)).const_mul 2)
      (by unfold X05.logit; fun_prop) ?_
    rw [uIoc_of_le zero_le_one, ← restrict_Ioo_eq_restrict_Ioc]
    refine ae_restrict_of_forall_mem measurableSet_Ioo fun t ht => ?_
    have hb := X05.sqrtcosh_bd hη ht.1 ht.2 (k := 0) (by norm_num)
    have hs := X05.sqrt_le_exp_half (η * X05.logit t)
    have hw : 0 ≤ (t * (1 - t)) ^ (η / 2) := by
      have : 0 < 1 - t := by linarith [ht.2]
      have := ht.1
      positivity
    rw [pow_zero, one_mul] at hb
    dsimp only
    rw [Real.norm_eq_abs, abs_of_nonneg (mul_nonneg hw (Real.sqrt_nonneg _))]
    exact (mul_le_mul_of_nonneg_left hs hw).trans hb
  set J := ∫ t in (0 : ℝ)..1, (t * (1 - t)) ^ (η / 2) * Real.sqrt (Real.cosh (η * X05.logit t))
  have hlin : J - X05.G η 0 0 - 1 / 4 * η ^ 2 * X05.G η 2 0 + 1 / 96 * η ^ 4 * X05.G η 4 0 =
      ∫ t in (0 : ℝ)..1, ((t * (1 - t)) ^ (η / 2) * Real.sqrt (Real.cosh (η * X05.logit t))
        - (t * (1 - t)) ^ (η / 2) * X05.logit t ^ 0 * Real.exp (0 * X05.logit t)
        - 1 / 4 * η ^ 2 * ((t * (1 - t)) ^ (η / 2) * X05.logit t ^ 2 * Real.exp (0 * X05.logit t))
        + 1 / 96 * η ^ 4 *
          ((t * (1 - t)) ^ (η / 2) * X05.logit t ^ 4 * Real.exp (0 * X05.logit t))) := by
    rw [intervalIntegral.integral_add ((iJ.sub i0).sub (i2.const_mul _)) (i4.const_mul _),
      intervalIntegral.integral_sub (iJ.sub i0) (i2.const_mul _),
      intervalIntegral.integral_sub iJ i0, intervalIntegral.integral_const_mul,
      intervalIntegral.integral_const_mul]
    rfl
  have hbound : |J - X05.G η 0 0 - 1 / 4 * η ^ 2 * X05.G η 2 0 + 1 / 96 * η ^ 4 * X05.G η 4 0|
      ≤ C4 * (2 * 96 ^ 6) * I * |η| ^ 6 := by
    rw [hlin, ← Real.norm_eq_abs]
    refine le_trans (intervalIntegral.norm_integral_le_of_norm_le
      (g := fun t => C4 * |η| ^ 6 * (2 * (96 ^ 6 * X05.q t))) zero_le_one ?_
      ((X05.q_int.const_mul (96 ^ 6)).const_mul 2 |>.const_mul (C4 * |η| ^ 6))) (le_of_eq ?_)
    rotate_left
    · rw [intervalIntegral.integral_const_mul, intervalIntegral.integral_const_mul,
        intervalIntegral.integral_const_mul]
      ring
    · filter_upwards [X05.ae_ne_one] with t ht1 htm
      have ht0 := htm.1
      have ht1' : t < 1 := lt_of_le_of_ne htm.2 ht1
      have hw : 0 ≤ (t * (1 - t)) ^ (η / 2) := by
        have : 0 < 1 - t := by linarith
        positivity
      have hb := X05.sqrtcosh_bd hη ht0 ht1' (k := 6) (by norm_num)
      have hc := hC4 (η * X05.logit t)
      have hr : (t * (1 - t)) ^ (η / 2) * Real.sqrt (Real.cosh (η * X05.logit t))
          - (t * (1 - t)) ^ (η / 2) * X05.logit t ^ 0 * Real.exp (0 * X05.logit t)
          - 1 / 4 * η ^ 2 * ((t * (1 - t)) ^ (η / 2) * X05.logit t ^ 2 * Real.exp (0 * X05.logit t))
          + 1 / 96 * η ^ 4 *
            ((t * (1 - t)) ^ (η / 2) * X05.logit t ^ 4 * Real.exp (0 * X05.logit t)) =
          (t * (1 - t)) ^ (η / 2) * (Real.sqrt (Real.cosh (η * X05.logit t)) -
            (1 + (η * X05.logit t) ^ 2 / 4 - (η * X05.logit t) ^ 4 / 96)) := by
        simp only [zero_mul, Real.exp_zero, pow_zero, mul_one]; ring
      rw [hr, Real.norm_eq_abs, abs_mul, abs_of_nonneg hw]
      calc (t * (1 - t)) ^ (η / 2) * |Real.sqrt (Real.cosh (η * X05.logit t)) -
            (1 + (η * X05.logit t) ^ 2 / 4 - (η * X05.logit t) ^ 4 / 96)|
          ≤ (t * (1 - t)) ^ (η / 2) * (C4 * |η * X05.logit t| ^ 6 *
              Real.exp (|η * X05.logit t| / 2)) := mul_le_mul_of_nonneg_left hc hw
        _ = C4 * |η| ^ 6 * (|X05.logit t| ^ 6 *
              ((t * (1 - t)) ^ (η / 2) * Real.exp (|η * X05.logit t| / 2))) := by
            rw [abs_mul, mul_pow]; ring
        _ ≤ _ := mul_le_mul_of_nonneg_left hb (by positivity)
  have e2 : 4 * Real.sqrt 2 * (1 + η) * J / (4 * Real.sqrt 2) -
      ((1 + η) * X05.G η 0 0 + 1 / 4 * η ^ 2 * ((1 + η) * X05.G η 2 0) -
        1 / 96 * η ^ 4 * ((1 + η) * X05.G η 4 0)) =
      (1 + η) * (J - X05.G η 0 0 - 1 / 4 * η ^ 2 * X05.G η 2 0 + 1 / 96 * η ^ 4 * X05.G η 4 0) := by
    have : Real.sqrt 2 ≠ 0 := by positivity
    field_simp
    ring
  rw [e2, abs_mul, abs_of_pos hη0]
  have hX := abs_nonneg (J - X05.G η 0 0 - 1 / 4 * η ^ 2 * X05.G η 2 0 + 1 / 96 * η ^ 4 * X05.G η 4 0)
  calc (1 + η) * |J - X05.G η 0 0 - 1 / 4 * η ^ 2 * X05.G η 2 0 + 1 / 96 * η ^ 4 * X05.G η 4 0|
      ≤ 5 / 4 * (C4 * (2 * 96 ^ 6) * I * |η| ^ 6) :=
        mul_le_mul hη5 hbound hX (by norm_num)
    _ = _ := by ring


namespace X05

lemma abs_mono_bd {a b x B : ℝ} (ha : |a| ≤ B) (hb : |b| ≤ B) (hx : |x| ≤ 1) (m : ℕ) :
    |a * b * x ^ m| ≤ B ^ 2 := by
  rw [abs_mul, abs_mul, abs_pow]
  have h1 : |x| ^ m ≤ 1 := pow_le_one₀ (abs_nonneg _) hx
  have h2 : |a| * |b| ≤ B * B := mul_le_mul ha hb (abs_nonneg _) ((abs_nonneg _).trans ha)
  calc |a| * |b| * |x| ^ m ≤ (B * B) * 1 :=
        mul_le_mul h2 h1 (by positivity) (by nlinarith [abs_nonneg a])
    _ = B ^ 2 := by ring

lemma abs_coef_bd {a x B : ℝ} (ha : |a| ≤ B) (hx : |x| ≤ 1) (m : ℕ) :
    |a * x ^ (m + 2)| ≤ B * x ^ 2 := by
  rw [abs_mul, abs_pow]
  have : |x| ^ (m + 2) ≤ x ^ 2 := by
    rw [pow_add, ← sq_abs x, mul_comm]
    exact mul_le_of_le_one_right (by positivity) (pow_le_one₀ (abs_nonneg _) hx)
  exact mul_le_mul ha this (by positivity) ((abs_nonneg _).trans ha)

lemma exp_trunc {L a0 a1 a2 a3 x B K : ℝ} (hB : 0 ≤ B) (hK : 0 ≤ K)
    (ha0 : |a0| ≤ B) (ha1 : |a1| ≤ B) (ha2 : |a2| ≤ B) (ha3 : |a3| ≤ B)
    (hx : |x| ≤ 1) (hsmall : (4 * B + K) * x ^ 2 ≤ 1)
    (hL : |L - (a0 * x ^ 2 + a1 * x ^ 3 + a2 * x ^ 4 + a3 * x ^ 5)| ≤ K * x ^ 6) :
    |Real.exp L - (1 + a0 * x ^ 2 + a1 * x ^ 3 + (a2 + a0 ^ 2 / 2) * x ^ 4
        + (a3 + a0 * a1) * x ^ 5)| ≤ (K + K * (8 * B + K) + 7 * B ^ 2 + (4 * B + K) ^ 3) * x ^ 6 := by
  obtain ⟨P, hPdef⟩ : ∃ P, P = a0 * x ^ 2 + a1 * x ^ 3 + a2 * x ^ 4 + a3 * x ^ 5 := ⟨_, rfl⟩
  rw [← hPdef] at hL
  have hx2 : 0 ≤ x ^ 2 := sq_nonneg x
  have hx21 : x ^ 2 ≤ 1 := by rw [← sq_abs]; exact pow_le_one₀ (abs_nonneg _) hx
  have hx6 : x ^ 6 ≤ x ^ 2 := by
    rw [show x ^ 6 = x ^ 2 * (x ^ 2) ^ 2 by ring]
    exact mul_le_of_le_one_right hx2 (pow_le_one₀ hx2 hx21)
  have hP : |P| ≤ 4 * B * x ^ 2 := by
    have e0 := abs_coef_bd ha0 hx 0
    have e1 := abs_coef_bd ha1 hx 1
    have e2 := abs_coef_bd ha2 hx 2
    have e3 := abs_coef_bd ha3 hx 3
    simp only [zero_add] at e0 e1 e2 e3
    have t1 := abs_add_le (a0 * x ^ 2 + a1 * x ^ 3 + a2 * x ^ 4) (a3 * x ^ 5)
    have t2 := abs_add_le (a0 * x ^ 2 + a1 * x ^ 3) (a2 * x ^ 4)
    have t3 := abs_add_le (a0 * x ^ 2) (a1 * x ^ 3)
    rw [hPdef]
    linarith
  have hLb : |L| ≤ (4 * B + K) * x ^ 2 := by
    have : |L| ≤ |P| + |L - P| := by
      have := abs_add_le P (L - P); rwa [add_sub_cancel] at this
    nlinarith
  have hL1 : |L| ≤ 1 := hLb.trans hsmall
  have he := Real.exp_bound hL1 (n := 3) (by norm_num)
  simp only [Finset.sum_range_succ, Finset.sum_range_zero, Nat.factorial] at he
  norm_num at he
  have hL3 : |L| ^ 3 ≤ (4 * B + K) ^ 3 * x ^ 6 := by
    calc |L| ^ 3 ≤ ((4 * B + K) * x ^ 2) ^ 3 := pow_le_pow_left₀ (abs_nonneg _) hLb 3
      _ = (4 * B + K) ^ 3 * x ^ 6 := by ring
  have hE1 : |Real.exp L - (1 + L + L ^ 2 / 2)| ≤ (4 * B + K) ^ 3 * x ^ 6 := by
    have : |L| ^ 3 * (2 / 9) ≤ |L| ^ 3 := by nlinarith [pow_nonneg (abs_nonneg L) 3]
    linarith
  have hE2 : |(L ^ 2 - P ^ 2) / 2| ≤ K * (8 * B + K) * x ^ 6 := by
    have hLP : |L + P| ≤ (8 * B + K) * x ^ 2 := by
      calc |L + P| ≤ |L| + |P| := abs_add_le _ _
        _ ≤ (8 * B + K) * x ^ 2 := by linarith
    rw [show (L ^ 2 - P ^ 2) / 2 = (L - P) * (L + P) / 2 by ring, abs_div, abs_mul,
      abs_of_pos (by norm_num : (0:ℝ) < 2)]
    have := mul_le_mul hL hLP (abs_nonneg _) (by positivity)
    have h8 : 0 ≤ K * (8 * B + K) * x ^ 6 := by positivity
    have : K * x ^ 6 * ((8 * B + K) * x ^ 2) ≤ K * (8 * B + K) * x ^ 6 := by
      have := mul_le_mul_of_nonneg_left hx21 h8
      nlinarith
    linarith
  obtain ⟨W, hW⟩ : ∃ W : ℝ, W = a1 ^ 2 / 2 + a0 * a2 + (a0 * a3 + a1 * a2) * x + (a2 ^ 2 / 2 + a1 * a3) * x ^ 2
    + a2 * a3 * x ^ 3 + a3 ^ 2 / 2 * x ^ 4 := ⟨_, rfl⟩
  have hWb : |W| ≤ 7 * B ^ 2 := by
    have b1 := abs_mono_bd ha1 ha1 hx 0
    have b2 := abs_mono_bd ha0 ha2 hx 0
    have b3 := abs_mono_bd ha0 ha3 hx 1
    have b4 := abs_mono_bd ha1 ha2 hx 1
    have b5 := abs_mono_bd ha2 ha2 hx 2
    have b6 := abs_mono_bd ha1 ha3 hx 2
    have b7 := abs_mono_bd ha2 ha3 hx 3
    have b8 := abs_mono_bd ha3 ha3 hx 4
    have eW : W = (a1 * a1 * x ^ 0) / 2 + a0 * a2 * x ^ 0 + a0 * a3 * x ^ 1 + a1 * a2 * x ^ 1
        + (a2 * a2 * x ^ 2) / 2 + a1 * a3 * x ^ 2 + a2 * a3 * x ^ 3 + (a3 * a3 * x ^ 4) / 2 := by
      rw [hW]; ring
    rw [eW]
    rw [abs_le] at b1 b2 b3 b4 b5 b6 b7 b8 ⊢
    constructor <;> linarith [b1.1, b1.2, b2.1, b2.2, b3.1, b3.2, b4.1, b4.2, b5.1, b5.2,
      b6.1, b6.2, b7.1, b7.2, b8.1, b8.2]
  have hE3 : |x ^ 6 * W| ≤ 7 * B ^ 2 * x ^ 6 := by
    rw [abs_mul, abs_of_nonneg (by positivity : (0:ℝ) ≤ x ^ 6), mul_comm]
    exact mul_le_mul_of_nonneg_right hWb (by positivity)
  have key : Real.exp L - (1 + a0 * x ^ 2 + a1 * x ^ 3 + (a2 + a0 ^ 2 / 2) * x ^ 4
        + (a3 + a0 * a1) * x ^ 5) =
      (Real.exp L - (1 + L + L ^ 2 / 2)) + (L - P) + (L ^ 2 - P ^ 2) / 2 + x ^ 6 * W := by
    rw [hW, hPdef]; ring
  rw [key]
  have t1 := abs_add_le ((Real.exp L - (1 + L + L ^ 2 / 2)) + (L - P) + (L ^ 2 - P ^ 2) / 2)
    (x ^ 6 * W)
  have t2 := abs_add_le ((Real.exp L - (1 + L + L ^ 2 / 2)) + (L - P)) ((L ^ 2 - P ^ 2) / 2)
  have t3 := abs_add_le (Real.exp L - (1 + L + L ^ 2 / 2)) (L - P)
  have hk : 0 ≤ K * (8 * B + K) * x ^ 6 := by positivity
  linarith

end X05

namespace X05

lemma z_two : z 2 = π ^ 2 / 6 := by
  unfold z
  rw [show ((2 : ℕ) : ℂ) = 2 by norm_num, riemannZeta_two]
  rw [show (π : ℂ) ^ 2 / 6 = ((π ^ 2 / 6 : ℝ) : ℂ) by push_cast; ring, Complex.ofReal_re]

lemma z_four : z 4 = π ^ 4 / 90 := by
  unfold z
  rw [show ((4 : ℕ) : ℂ) = 4 by norm_num, riemannZeta_four]
  rw [show (π : ℂ) ^ 4 / 90 = ((π ^ 4 / 90 : ℝ) : ℂ) by push_cast; ring, Complex.ofReal_re]

lemma z_bd {n : ℕ} (hn : 2 ≤ n) : |z n| ≤ 2 := by
  have h := z_hasSum (k := n) (by omega)
  have h2 := z_hasSum (k := 2) (by omega)
  have hnn : 0 ≤ z n := h.nonneg fun m => by positivity
  have hle : z n ≤ z 2 := hasSum_le (fun m => by
    have h1 : (1 : ℝ) ≤ (m : ℝ) + 1 := by have := m.cast_nonneg (α := ℝ); linarith
    exact one_div_le_one_div_of_le (by positivity) (pow_le_pow_right₀ h1 hn)) h h2
  rw [z_two] at hle
  have := Real.pi_lt_d2
  have : π ^ 2 / 6 ≤ 2 := by nlinarith [Real.pi_pos]
  rw [abs_of_nonneg hnn]; linarith

/-- coefficient `a_k` of the one-variable series -/
noncomputable def cf (α β γ : ℝ) (k : ℕ) : ℝ :=
  (-1 : ℝ) ^ (k + 2) * z (k + 2) / ((k + 2 : ℕ) : ℝ) * (α ^ (k + 2) + β ^ (k + 2) - γ ^ (k + 2))

lemma cf_bd {α β γ : ℝ} (hα : |α| ≤ 3 / 2) (hβ : |β| ≤ 3 / 2) (hγ : |γ| ≤ 3 / 2) (k : ℕ) :
    |cf α β γ k| ≤ 3 * (3 / 2) ^ (k + 2) := by
  unfold cf
  have hz := z_bd (n := k + 2) (by omega)
  have hn : (2 : ℝ) ≤ ((k + 2 : ℕ) : ℝ) := by push_cast; have := k.cast_nonneg (α := ℝ); linarith
  have hb : |α ^ (k + 2) + β ^ (k + 2) - γ ^ (k + 2)| ≤ 3 * (3 / 2) ^ (k + 2) := by
    have a := pow_le_pow_left₀ (abs_nonneg α) hα (k + 2)
    have b := pow_le_pow_left₀ (abs_nonneg β) hβ (k + 2)
    have c := pow_le_pow_left₀ (abs_nonneg γ) hγ (k + 2)
    rw [← abs_pow] at a b c
    have t1 := abs_sub (α ^ (k + 2) + β ^ (k + 2)) (γ ^ (k + 2))
    have t2 := abs_add_le (α ^ (k + 2)) (β ^ (k + 2))
    linarith
  rw [abs_mul, abs_div, abs_mul, abs_pow, abs_neg, abs_one, one_pow, one_mul,
    abs_of_pos (by linarith : (0:ℝ) < ((k + 2 : ℕ) : ℝ))]
  have : |z (k + 2)| / ((k + 2 : ℕ) : ℝ) ≤ 1 := by
    rw [div_le_one (by linarith)]; linarith
  calc |z (k + 2)| / ((k + 2 : ℕ) : ℝ) * |α ^ (k + 2) + β ^ (k + 2) - γ ^ (k + 2)|
      ≤ 1 * (3 * (3 / 2) ^ (k + 2)) := mul_le_mul this hb (abs_nonneg _) (by norm_num)
    _ = _ := one_mul _

lemma series_trunc {α β γ x L : ℝ} (hα : |α| ≤ 3 / 2) (hβ : |β| ≤ 3 / 2) (hγ : |γ| ≤ 3 / 2)
    (hx : |x| ≤ 1 / 3)
    (hL : HasSum (fun k : ℕ => (-1 : ℝ) ^ (k + 2) * z (k + 2) / ((k + 2 : ℕ) : ℝ) *
      ((α * x) ^ (k + 2) + (β * x) ^ (k + 2) - (γ * x) ^ (k + 2))) L) :
    |L - (cf α β γ 0 * x ^ 2 + cf α β γ 1 * x ^ 3 + cf α β γ 2 * x ^ 4 + cf α β γ 3 * x ^ 5)|
      ≤ 69 * x ^ 6 := by
  have hf : (fun k : ℕ => (-1 : ℝ) ^ (k + 2) * z (k + 2) / ((k + 2 : ℕ) : ℝ) *
      ((α * x) ^ (k + 2) + (β * x) ^ (k + 2) - (γ * x) ^ (k + 2))) =
      fun k => cf α β γ k * x ^ (k + 2) := by
    ext k; unfold cf; rw [mul_pow, mul_pow, mul_pow]; ring
  rw [hf] at hL
  set r := 3 / 2 * |x| with hr
  have hr0 : 0 ≤ r := by positivity
  have hr1 : r ≤ 1 / 2 := by rw [hr]; linarith
  have hbd : ∀ k, |cf α β γ k * x ^ (k + 2)| ≤ 3 * r ^ (k + 2) := by
    intro k
    rw [abs_mul, abs_pow, hr, mul_pow]
    have := cf_bd hα hβ hγ k
    calc |cf α β γ k| * |x| ^ (k + 2) ≤ 3 * (3 / 2) ^ (k + 2) * |x| ^ (k + 2) :=
          mul_le_mul_of_nonneg_right this (by positivity)
      _ = _ := by ring
  have hg := (hasSum_nat_add_iff' 4).mpr hL
  have hgeo : HasSum (fun k : ℕ => 3 * r ^ 6 * r ^ k) (3 * r ^ 6 * (1 - r)⁻¹) :=
    (hasSum_geometric_of_lt_one hr0 (by linarith)).mul_left _
  have hb1 : ∀ k : ℕ, cf α β γ (k + 4) * x ^ (k + 4 + 2) ≤ 3 * r ^ 6 * r ^ k := by
    intro k
    have := (abs_le.mp (hbd (k + 4))).2
    calc _ ≤ 3 * r ^ (k + 4 + 2) := this
      _ = _ := by ring
  have hb2 : ∀ k : ℕ, -(cf α β γ (k + 4) * x ^ (k + 4 + 2)) ≤ 3 * r ^ 6 * r ^ k := by
    intro k
    have := (abs_le.mp (hbd (k + 4))).1
    calc _ ≤ 3 * r ^ (k + 4 + 2) := by linarith
      _ = _ := by ring
  have u1 := hasSum_le hb1 hg hgeo
  have u2 := hasSum_le hb2 hg.neg hgeo
  have hinv : 3 * r ^ 6 * (1 - r)⁻¹ ≤ 6 * r ^ 6 := by
    rw [mul_inv_le_iff₀ (by linarith)]
    nlinarith [pow_nonneg hr0 6]
  have hr6 : 6 * r ^ 6 ≤ 69 * x ^ 6 := by
    rw [hr, mul_pow, pow_abs, abs_of_nonneg (by positivity : (0:ℝ) ≤ x ^ 6)]
    nlinarith [pow_nonneg (sq_nonneg x) 3]
  simp only [Finset.sum_range_succ, Finset.sum_range_zero, zero_add] at u1 u2
  rw [abs_le]
  constructor <;> norm_num at u1 u2 ⊢ <;> linarith

lemma A_expand {α β γ x L : ℝ} (hα : |α| ≤ 3 / 2) (hβ : |β| ≤ 3 / 2) (hγ : |γ| ≤ 3 / 2)
    (hx : |x| ≤ 1 / 20)
    (hL : HasSum (fun k : ℕ => (-1 : ℝ) ^ (k + 2) * z (k + 2) / ((k + 2 : ℕ) : ℝ) *
      ((α * x) ^ (k + 2) + (β * x) ^ (k + 2) - (γ * x) ^ (k + 2))) L) :
    |Real.exp L - (1 + cf α β γ 0 * x ^ 2 + cf α β γ 1 * x ^ 3
        + (cf α β γ 2 + cf α β γ 0 ^ 2 / 2) * x ^ 4
        + (cf α β γ 3 + cf α β γ 0 * cf α β γ 1) * x ^ 5)|
      ≤ (69 + 69 * (8 * 23 + 69) + 7 * 23 ^ 2 + (4 * 23 + 69) ^ 3) * x ^ 6 := by
  have h := series_trunc hα hβ hγ (by linarith) hL
  have hB : ∀ k ≤ 3, |cf α β γ k| ≤ 23 := by
    intro k hk
    have := cf_bd hα hβ hγ k
    have : (3 / 2 : ℝ) ^ (k + 2) ≤ (3 / 2) ^ 5 := pow_le_pow_right₀ (by norm_num) (by omega)
    norm_num at this ⊢; linarith
  have hx2 : x ^ 2 ≤ 1 / 400 := by
    have := pow_le_pow_left₀ (abs_nonneg x) hx 2
    rw [sq_abs] at this; norm_num at this ⊢; linarith
  exact exp_trunc (by norm_num) (by norm_num) (hB 0 (by norm_num)) (hB 1 (by norm_num))
    (hB 2 (by norm_num)) (hB 3 (by norm_num)) (by linarith) (by nlinarith) h

end X05

namespace X05

lemma real_exp_rem (x : ℝ) :
    |Real.exp x - ∑ m ∈ Finset.range 6, x ^ m / (m.factorial : ℝ)| ≤ |x| ^ 6 * Real.exp |x| := by
  have h := Complex.norm_exp_sub_sum_le_norm_mul_exp (x : ℂ) 6
  have e : Complex.exp (x : ℂ) - ∑ m ∈ Finset.range 6, (x : ℂ) ^ m / (m.factorial : ℂ) =
      ((Real.exp x - ∑ m ∈ Finset.range 6, x ^ m / (m.factorial : ℝ) : ℝ) : ℂ) := by
    push_cast; rfl
  rw [e, Complex.norm_real, Real.norm_eq_abs, Complex.norm_real, Real.norm_eq_abs] at h
  exact h

lemma s_taylor {η s : ℝ} (hη : |η| ≤ 1 / 4) (hs : |s| ≤ 1 / 8) :
    |A η s - ∑ m ∈ Finset.range 6, s ^ m / (m.factorial : ℝ) *
        iteratedDeriv m (fun s => A η s) 0|
      ≤ 5 / 4 * (2 * 96 ^ 6) * (∫ t in (0 : ℝ)..1, q t) * |s| ^ 6 := by
  have hη0 : 0 < 1 + η := by rw [abs_le] at hη; linarith
  have hη5 : 1 + η ≤ 5 / 4 := by rw [abs_le] at hη; linarith
  have hA : A η s = (1 + η) * G η 0 s := congrFun (A_eq_G η) s
  have hD : ∀ m ∈ Finset.range 6, s ^ m / (m.factorial : ℝ) *
      iteratedDeriv m (fun s => A η s) 0 = (1 + η) * (s ^ m / (m.factorial : ℝ) * G η m 0) := by
    intro m hm
    rw [iter_A hη m (by simp at hm; omega) 0 (by norm_num)]; ring
  rw [Finset.sum_congr rfl hD, ← Finset.mul_sum, hA, ← mul_sub, abs_mul, abs_of_pos hη0]
  have hint : ∀ m ∈ Finset.range 6, IntervalIntegrable (fun t => s ^ m / (m.factorial : ℝ) *
      ((t * (1 - t)) ^ (η / 2) * logit t ^ m * Real.exp (0 * logit t))) volume 0 1 :=
    fun m hm => (G_int hη (s := 0) (by norm_num) (k := m) (by simp at hm; omega)).const_mul _
  have hsum : ∑ m ∈ Finset.range 6, s ^ m / (m.factorial : ℝ) * G η m 0 =
      ∫ t in (0 : ℝ)..1, ∑ m ∈ Finset.range 6, s ^ m / (m.factorial : ℝ) *
        ((t * (1 - t)) ^ (η / 2) * logit t ^ m * Real.exp (0 * logit t)) := by
    rw [intervalIntegral.integral_finsetSum hint]
    refine Finset.sum_congr rfl fun m _ => ?_
    rw [G, intervalIntegral.integral_const_mul]
  have i0 := G_int hη hs (k := 0) (by norm_num)
  have isum : IntervalIntegrable (fun t => ∑ m ∈ Finset.range 6, s ^ m / (m.factorial : ℝ) *
      ((t * (1 - t)) ^ (η / 2) * logit t ^ m * Real.exp (0 * logit t))) volume 0 1 := by
    have := IntervalIntegrable.sum (Finset.range 6) hint
    convert this using 1
    ext t; simp [Finset.sum_apply]
  rw [hsum, G, ← intervalIntegral.integral_sub i0 isum]
  have hbound : ‖∫ t in (0 : ℝ)..1, ((t * (1 - t)) ^ (η / 2) * logit t ^ 0 * Real.exp (s * logit t)
      - ∑ m ∈ Finset.range 6, s ^ m / (m.factorial : ℝ) *
        ((t * (1 - t)) ^ (η / 2) * logit t ^ m * Real.exp (0 * logit t)))‖ ≤
      ∫ t in (0 : ℝ)..1, |s| ^ 6 * (2 * (96 ^ 6 * q t)) := by
    refine intervalIntegral.norm_integral_le_of_norm_le zero_le_one ?_
      ((q_int.const_mul (96 ^ 6)).const_mul 2 |>.const_mul _)
    filter_upwards [ae_ne_one] with t ht1 htm
    have ht0 := htm.1
    have ht1' : t < 1 := lt_of_le_of_ne htm.2 ht1
    have hw : 0 ≤ (t * (1 - t)) ^ (η / 2) := by
      have : 0 < 1 - t := by linarith
      positivity
    have e : (t * (1 - t)) ^ (η / 2) * logit t ^ 0 * Real.exp (s * logit t)
        - ∑ m ∈ Finset.range 6, s ^ m / (m.factorial : ℝ) *
          ((t * (1 - t)) ^ (η / 2) * logit t ^ m * Real.exp (0 * logit t)) =
        (t * (1 - t)) ^ (η / 2) * (Real.exp (s * logit t) -
          ∑ m ∈ Finset.range 6, (s * logit t) ^ m / (m.factorial : ℝ)) := by
      simp only [Finset.sum_range_succ, Finset.sum_range_zero, zero_mul, Real.exp_zero]
      ring
    rw [e, Real.norm_eq_abs, abs_mul, abs_of_nonneg hw]
    have r := real_exp_rem (s * logit t)
    have hs' : |-s| ≤ 1 / 8 := by rwa [abs_neg]
    have a := G_integrand_bd hη hs ht0 ht1' (k := 6) (by norm_num)
    have b := G_integrand_bd hη hs' ht0 ht1' (k := 6) (by norm_num)
    rw [Real.norm_eq_abs, abs_mul, abs_mul, abs_pow, abs_of_nonneg hw,
      abs_of_pos (Real.exp_pos _)] at a b
    have he : Real.exp |s * logit t| ≤ Real.exp (s * logit t) + Real.exp (-s * logit t) := by
      rcases le_total 0 (s * logit t) with h | h
      · rw [abs_of_nonneg h]; linarith [Real.exp_pos (-s * logit t)]
      · rw [abs_of_nonpos h, show -(s * logit t) = -s * logit t by ring]
        linarith [Real.exp_pos (s * logit t)]
    calc (t * (1 - t)) ^ (η / 2) * |Real.exp (s * logit t) -
          ∑ m ∈ Finset.range 6, (s * logit t) ^ m / (m.factorial : ℝ)|
        ≤ (t * (1 - t)) ^ (η / 2) * (|s * logit t| ^ 6 * Real.exp |s * logit t|) :=
          mul_le_mul_of_nonneg_left r hw
      _ ≤ (t * (1 - t)) ^ (η / 2) * (|s * logit t| ^ 6 *
            (Real.exp (s * logit t) + Real.exp (-s * logit t))) :=
          mul_le_mul_of_nonneg_left (mul_le_mul_of_nonneg_left he (by positivity)) hw
      _ = |s| ^ 6 * ((t * (1 - t)) ^ (η / 2) * |logit t| ^ 6 * Real.exp (s * logit t) +
            (t * (1 - t)) ^ (η / 2) * |logit t| ^ 6 * Real.exp (-s * logit t)) := by
          rw [abs_mul, mul_pow]; ring
      _ ≤ |s| ^ 6 * (2 * (96 ^ 6 * q t)) :=
          mul_le_mul_of_nonneg_left (by linarith) (by positivity)
  rw [Real.norm_eq_abs] at hbound
  rw [intervalIntegral.integral_const_mul, intervalIntegral.integral_const_mul,
    intervalIntegral.integral_const_mul] at hbound
  have hX := abs_nonneg (∫ t in (0 : ℝ)..1, ((t * (1 - t)) ^ (η / 2) * logit t ^ 0 *
    Real.exp (s * logit t) - ∑ m ∈ Finset.range 6, s ^ m / (m.factorial : ℝ) *
      ((t * (1 - t)) ^ (η / 2) * logit t ^ m * Real.exp (0 * logit t))))
  calc _ ≤ 5 / 4 * (|s| ^ 6 * (2 * (96 ^ 6 * ∫ t in (0 : ℝ)..1, q t))) :=
        mul_le_mul hη5 hbound hX (by norm_num)
    _ = _ := by ring

lemma A_neg {η s : ℝ} (hη : -1 < η) (hs : |s| < 1 + η / 2) : A η (-s) = A η s := by
  rw [proof_2 η s hη hs, proof_2 η (-s) hη (by rwa [abs_neg]), sub_neg_eq_add,
    ← sub_eq_add_neg, mul_comm (Real.Gamma (1 + η / 2 - s))]

lemma s_even {η s : ℝ} (hη : |η| ≤ 1 / 4) (hs : |s| ≤ 1 / 8) :
    |A η s - (iteratedDeriv 0 (fun s => A η s) 0 + s ^ 2 / 2 * iteratedDeriv 2 (fun s => A η s) 0
      + s ^ 4 / 24 * iteratedDeriv 4 (fun s => A η s) 0)|
      ≤ 5 / 4 * (2 * 96 ^ 6) * (∫ t in (0 : ℝ)..1, q t) * |s| ^ 6 := by
  have h1 := s_taylor hη hs
  have hs' : |-s| ≤ 1 / 8 := by rwa [abs_neg]
  have h2 := s_taylor hη hs'
  rw [A_neg (by rw [abs_le] at hη; linarith)
    (by rw [abs_le] at hη; linarith), abs_neg] at h2
  simp only [Finset.sum_range_succ, Finset.sum_range_zero, Nat.factorial] at h1 h2
  rw [abs_le] at h1 h2 ⊢
  push_cast at h1 h2 ⊢
  constructor <;> nlinarith [h1.1, h1.2, h2.1, h2.2]

lemma eq_zero_of_le {a K d : ℝ} (hd : 0 < d) (hK : 0 ≤ K)
    (h : ∀ s : ℝ, 0 < s → s ≤ d → |a| ≤ K * s ^ 2) : a = 0 := by
  by_contra ha
  have ha' : 0 < |a| := abs_pos.mpr ha
  set s := min d (min 1 (|a| / (K + 1))) with hs
  have hs0 : 0 < s := lt_min hd (lt_min one_pos (by positivity))
  have hs1 : s ≤ 1 := (min_le_right _ _).trans (min_le_left _ _)
  have hs2 : s ≤ |a| / (K + 1) := (min_le_right _ _).trans (min_le_right _ _)
  have := h s hs0 (min_le_left _ _)
  have hsq : s ^ 2 ≤ s := by nlinarith
  have : K * s ^ 2 ≤ K * (|a| / (K + 1)) := mul_le_mul_of_nonneg_left (hsq.trans hs2) hK
  have : K * (|a| / (K + 1)) < |a| := by
    rw [mul_div_assoc', div_lt_iff₀ (by linarith)]; nlinarith
  linarith

end X05

namespace X05

/-- the constant of `A_expand` -/
noncomputable def C1 : ℝ := 69 + 69 * (8 * 23 + 69) + 7 * 23 ^ 2 + (4 * 23 + 69) ^ 3

/-- the constant of `s_even` -/
noncomputable def K6 : ℝ := 5 / 4 * (2 * 96 ^ 6) * (∫ t in (0 : ℝ)..1, q t)

lemma K6_nonneg : 0 ≤ K6 := by
  unfold K6
  have : 0 ≤ ∫ t in (0 : ℝ)..1, q t :=
    intervalIntegral.integral_nonneg zero_le_one fun t ht =>
      Real.rpow_nonneg (mul_nonneg ht.1 (by linarith [ht.2])) _
  exact mul_nonneg (by norm_num) this

lemma C1_pos : 0 < C1 := by unfold C1; norm_num

lemma expand0 {η : ℝ} (hη : |η| ≤ 1 / 20) :
    |A η 0 - (1 - (π ^ 2 / 6) / 4 * η ^ 2 + z 3 / 4 * η ^ 3
        + ((π ^ 2 / 6) ^ 2 / 32 - 7 * (π ^ 4 / 90) / 32) * η ^ 4
        + (-(π ^ 2 / 6) * z 3 / 16 + 3 * z 5 / 16) * η ^ 5)| ≤ C1 * η ^ 6 := by
  obtain ⟨hpos, h⟩ := proof_3 η 0 (by linarith) (by norm_num)
  rw [show η / 2 + 0 = 1 / 2 * η by ring, show η / 2 - 0 = 1 / 2 * η by ring] at h
  have := A_expand (α := 1 / 2) (β := 1 / 2) (γ := 1) (by norm_num) (by norm_num) (by norm_num) hη
    (by simpa only [one_mul] using h)
  rw [Real.exp_log hpos] at this
  rw [show C1 = 69 + 69 * (8 * 23 + 69) + 7 * 23 ^ 2 + (4 * 23 + 69) ^ 3 from rfl]
  convert this using 3
  norm_num [cf, z_two, z_four]
  ring

lemma expand1 {η : ℝ} (hη : |η| ≤ 1 / 20) :
    |A η η - (1 + cf (3 / 2) (-1 / 2) 1 0 * η ^ 2 + cf (3 / 2) (-1 / 2) 1 1 * η ^ 3
        + (cf (3 / 2) (-1 / 2) 1 2 + cf (3 / 2) (-1 / 2) 1 0 ^ 2 / 2) * η ^ 4
        + (cf (3 / 2) (-1 / 2) 1 3 + cf (3 / 2) (-1 / 2) 1 0 * cf (3 / 2) (-1 / 2) 1 1) * η ^ 5)|
      ≤ C1 * η ^ 6 := by
  obtain ⟨hpos, h⟩ := proof_3 η η (by linarith [abs_nonneg η]) (by linarith [abs_nonneg η])
  rw [show η / 2 + η = 3 / 2 * η by ring, show η / 2 - η = -1 / 2 * η by ring] at h
  have := A_expand (α := 3 / 2) (β := -1 / 2) (γ := 1) (by norm_num) (by norm_num) (by norm_num) hη
    (by simpa only [one_mul] using h)
  rw [Real.exp_log hpos] at this
  exact this

lemma expand_half {η : ℝ} (hη : |η| ≤ 1 / 4) : A η (η / 2) = 1 := by
  rw [abs_le] at hη
  rw [proof_2 η (η / 2) (by linarith) (by rw [abs_lt]; constructor <;> linarith)]
  rw [show 1 + η / 2 + η / 2 = 1 + η by ring, show 1 + η / 2 - η / 2 = 1 by ring,
    Real.Gamma_one, mul_one, div_self (Real.Gamma_pos_of_pos (by linarith)).ne']

lemma M_zero :
    iteratedDeriv 2 (fun s => A 0 s) 0 = 2 * (π ^ 2 / 6) ∧
    iteratedDeriv 4 (fun s => A 0 s) 0 = 12 * (π ^ 2 / 6) ^ 2 + 12 * (π ^ 4 / 90) := by
  set M2 := iteratedDeriv 2 (fun s => A 0 s) 0
  set M4 := iteratedDeriv 4 (fun s => A 0 s) 0
  have hA00 : A 0 0 = 1 := by
    rw [proof_2 0 0 (by norm_num) (by norm_num)]; simp
  have key : ∀ s : ℝ, 0 < s → s ≤ 1 / 20 →
      |(M2 / 2 - π ^ 2 / 6) * s ^ 2 + (M4 / 24 - (π ^ 4 / 90) / 2 - (π ^ 2 / 6) ^ 2 / 2) * s ^ 4|
        ≤ (C1 + K6) * s ^ 6 := by
    intro s hs0 hs1
    have hsa : |s| ≤ 1 / 20 := by rw [abs_of_pos hs0]; exact hs1
    obtain ⟨hpos, h⟩ := proof_3 0 s (by norm_num) (by rw [abs_of_pos hs0]; linarith)
    simp only [zero_div, zero_add, zero_sub] at h
    have e1 := A_expand (α := 1) (β := -1) (γ := 0) (by norm_num) (by norm_num) (by norm_num) hsa
      (by simpa only [one_mul, neg_one_mul, zero_mul] using h)
    rw [Real.exp_log hpos] at e1
    have e2 := s_even (η := 0) (s := s) (by norm_num) (by linarith)
    rw [iteratedDeriv_zero, hA00] at e2
    have hc0 : cf 1 (-1) 0 0 = π ^ 2 / 6 := by norm_num [cf, z_two]
    have hc1 : cf 1 (-1) 0 1 = 0 := by norm_num [cf]
    have hc2 : cf 1 (-1) 0 2 = (π ^ 4 / 90) / 2 := by norm_num [cf, z_four]; ring
    have hc3 : cf 1 (-1) 0 3 = 0 := by norm_num [cf]
    rw [hc0, hc1, hc2, hc3,
      show (69 + 69 * (8 * 23 + 69) + 7 * 23 ^ 2 + (4 * 23 + 69) ^ 3 : ℝ) = C1 from rfl] at e1
    rw [abs_of_pos hs0] at e2
    change |A 0 s - (1 + s ^ 2 / 2 * M2 + s ^ 4 / 24 * M4)| ≤ K6 * s ^ 6 at e2
    rw [abs_le] at e1 e2 ⊢
    constructor <;> linarith [e1.1, e1.2, e2.1, e2.2]
  have hCK : 0 ≤ C1 + K6 := by linarith [C1_pos, K6_nonneg]
  set α := M2 / 2 - π ^ 2 / 6
  set β := M4 / 24 - (π ^ 4 / 90) / 2 - (π ^ 2 / 6) ^ 2 / 2
  have hα : α = 0 := by
    refine eq_zero_of_le (K := |β| + (C1 + K6)) (d := 1 / 20) (by norm_num) (by positivity)
      fun s hs0 hs1 => ?_
    have h := key s hs0 hs1
    have hs2 : 0 < s ^ 2 := by positivity
    have hs21 : s ^ 2 ≤ 1 := by nlinarith
    have e : α * s ^ 2 + β * s ^ 4 = s ^ 2 * (α + β * s ^ 2) := by ring
    rw [e, abs_mul, abs_of_pos hs2] at h
    have h' : |α + β * s ^ 2| ≤ (C1 + K6) * s ^ 4 := by
      have : (C1 + K6) * s ^ 6 = s ^ 2 * ((C1 + K6) * s ^ 4) := by ring
      rw [this] at h
      exact le_of_mul_le_mul_left h hs2
    have t := abs_sub (α + β * s ^ 2) (β * s ^ 2)
    rw [add_sub_cancel_right, abs_mul, abs_of_pos hs2] at t
    have : (C1 + K6) * s ^ 4 ≤ (C1 + K6) * s ^ 2 := by
      apply mul_le_mul_of_nonneg_left _ hCK
      nlinarith
    nlinarith [abs_nonneg β]
  have hβ : β = 0 := by
    refine eq_zero_of_le (K := C1 + K6) (d := 1 / 20) (by norm_num) hCK fun s hs0 hs1 => ?_
    have h := key s hs0 hs1
    rw [hα, zero_mul, zero_add, abs_mul, abs_of_pos (by positivity : 0 < s ^ 4)] at h
    have hs4 : 0 < s ^ 4 := by positivity
    have : (C1 + K6) * s ^ 6 = s ^ 4 * ((C1 + K6) * s ^ 2) := by ring
    rw [this, mul_comm] at h
    exact le_of_mul_le_mul_left h hs4
  constructor
  · have : M2 / 2 - π ^ 2 / 6 = 0 := hα
    linarith
  · have : M4 / 24 - (π ^ 4 / 90) / 2 - (π ^ 2 / 6) ^ 2 / 2 = 0 := hβ
    linarith

end X05

namespace X05

lemma div_bound {X e W : ℝ} (he : 0 < e) (h : |e * X| ≤ e * W) : |X| ≤ W := by
  rw [abs_mul, abs_of_pos he] at h
  exact le_of_mul_le_mul_left h he

lemma alg_b {a0 a1 M2 M4 T0 T1 pb η C K : ℝ} (hη : η ≠ 0)
    (hid : 16 - 15 * T0 - T1 = 3 / 2 * η ^ 2 * pb)
    (E0 : |a0 - T0| ≤ C * η ^ 6) (E1 : |a1 - T1| ≤ C * η ^ 6)
    (D1 : |a1 - (a0 + η ^ 2 / 2 * M2 + η ^ 4 / 24 * M4)| ≤ K * η ^ 6)
    (Dh : |1 - (a0 + (η / 2) ^ 2 / 2 * M2 + (η / 2) ^ 4 / 24 * M4)| ≤ K * (η ^ 6 / 64)) :
    |M2 - pb| ≤ 2 / 3 * (5 * K / 4 + 16 * C) * |η| ^ 4 := by
  have hpos : 0 < 3 / 2 * η ^ 2 := by positivity
  refine div_bound hpos ?_
  have e : 3 / 2 * η ^ 2 * (M2 - pb) = (a1 - (a0 + η ^ 2 / 2 * M2 + η ^ 4 / 24 * M4))
      - 16 * (1 - (a0 + (η / 2) ^ 2 / 2 * M2 + (η / 2) ^ 4 / 24 * M4)) - (a1 - T1)
      - 15 * (a0 - T0) := by linear_combination hid
  have e4 : |η| ^ 4 = η ^ 4 := by rw [show |η| ^ 4 = (|η| ^ 2) ^ 2 by ring, sq_abs]; ring
  rw [e, e4]
  generalize a1 - (a0 + η ^ 2 / 2 * M2 + η ^ 4 / 24 * M4) = d1 at D1
  generalize 1 - (a0 + (η / 2) ^ 2 / 2 * M2 + (η / 2) ^ 4 / 24 * M4) = dh at Dh
  generalize a1 - T1 = f1 at E1
  generalize a0 - T0 = f0 at E0
  have : 3 / 2 * η ^ 2 * (2 / 3 * (5 * K / 4 + 16 * C) * η ^ 4) =
      K * η ^ 6 + 16 * (K * (η ^ 6 / 64)) + C * η ^ 6 + 15 * (C * η ^ 6) := by ring
  rw [this]
  rw [abs_le] at D1 Dh E1 E0 ⊢
  constructor <;> linarith [D1.1, D1.2, Dh.1, Dh.2, E1.1, E1.2, E0.1, E0.2]

lemma alg_c {a0 a1 M2 M4 T0 T1 pc w η C K : ℝ} (hη : η ≠ 0) (hη1 : |η| ≤ 1) (hC : 0 ≤ C)
    (hK : 0 ≤ K) (hw : |w| ≤ 6)
    (hid : 4 * T1 - 16 + 12 * T0 = η ^ 4 / 8 * pc - 3 * w * η ^ 5)
    (E0 : |a0 - T0| ≤ C * η ^ 6) (E1 : |a1 - T1| ≤ C * η ^ 6)
    (D1 : |a1 - (a0 + η ^ 2 / 2 * M2 + η ^ 4 / 24 * M4)| ≤ K * η ^ 6)
    (Dh : |1 - (a0 + (η / 2) ^ 2 / 2 * M2 + (η / 2) ^ 4 / 24 * M4)| ≤ K * (η ^ 6 / 64)) :
    |M4 - pc| ≤ 8 * (18 + 16 * C + 17 * K / 4) * |η| := by
  have hpos : 0 < η ^ 4 / 8 := by positivity
  refine div_bound hpos ?_
  have e : η ^ 4 / 8 * (M4 - pc) = 4 * (a1 - T1) + 12 * (a0 - T0) - 3 * w * η ^ 5
      - 4 * (a1 - (a0 + η ^ 2 / 2 * M2 + η ^ 4 / 24 * M4))
      + 16 * (1 - (a0 + (η / 2) ^ 2 / 2 * M2 + (η / 2) ^ 4 / 24 * M4)) := by
    linear_combination hid
  rw [e]
  generalize a1 - (a0 + η ^ 2 / 2 * M2 + η ^ 4 / 24 * M4) = d1 at D1
  generalize 1 - (a0 + (η / 2) ^ 2 / 2 * M2 + (η / 2) ^ 4 / 24 * M4) = dh at Dh
  generalize a1 - T1 = f1 at E1
  generalize a0 - T0 = f0 at E0
  have e4 : η ^ 4 = |η| ^ 4 := by rw [show |η| ^ 4 = (|η| ^ 2) ^ 2 by ring, sq_abs]; ring
  have e6 : η ^ 6 = |η| ^ 6 := by rw [show |η| ^ 6 = (|η| ^ 2) ^ 3 by ring, sq_abs]; ring
  have h65 : |η| ^ 6 ≤ |η| ^ 5 := pow_le_pow_of_le_one (abs_nonneg _) hη1 (by norm_num)
  have hw5 : |3 * w * η ^ 5| ≤ 18 * |η| ^ 5 := by
    rw [abs_mul, abs_mul, abs_pow]
    have := mul_le_mul_of_nonneg_right hw (pow_nonneg (abs_nonneg η) 5)
    norm_num; linarith
  rw [e6] at D1 Dh E1 E0
  have : η ^ 4 / 8 * (8 * (18 + 16 * C + 17 * K / 4) * |η|) =
      (18 + 16 * C + 17 * K / 4) * |η| ^ 5 := by rw [e4]; ring
  rw [this]
  have hC6 := mul_le_mul_of_nonneg_left h65 hC
  have hK6 := mul_le_mul_of_nonneg_left h65 hK
  rw [abs_le] at D1 Dh E1 E0 hw5 ⊢
  constructor <;> linarith [D1.1, D1.2, Dh.1, Dh.2, E1.1, E1.2, E0.1, E0.2, hw5.1, hw5.2]

end X05

theorem X05.proof_6 : X05.stmt_6 := by
  have hC1 := X05.C1_pos
  have hK := X05.K6_nonneg
  have hpi := Real.pi_lt_d2
  have hpi0 := Real.pi_pos
  have hz2 : π ^ 2 / 6 ≤ 2 := by nlinarith
  have hz3 := X05.z_bd (n := 3) (by norm_num)
  have hz5 := X05.z_bd (n := 5) (by norm_num)
  refine ⟨X05.C1 + 2 / 3 * (5 * X05.K6 / 4 + 16 * X05.C1) + 8 * (18 + 16 * X05.C1 + 17 * X05.K6 / 4),
    1 / 20, by norm_num, fun η hη => ?_⟩
  have hη4 : |η| ≤ 1 / 4 := by linarith
  have hη8 : |η| ≤ 1 / 8 := by linarith
  have hη1 : |η| ≤ 1 := by linarith
  have hη2' : |η / 2| ≤ 1 / 8 := by rw [abs_div]; norm_num; linarith
  have E0 := X05.expand0 hη
  have E1 := X05.expand1 hη
  have Eh := X05.expand_half hη4
  have D1 := X05.s_even hη4 hη8
  have Dh := X05.s_even hη4 hη2'
  change _ ≤ X05.K6 * _ at D1 Dh
  rw [iteratedDeriv_zero] at D1 Dh
  rw [Eh] at Dh
  have e6 : |η| ^ 6 = η ^ 6 := by rw [show |η| ^ 6 = (|η| ^ 2) ^ 3 by ring, sq_abs]; ring
  have eh6 : |η / 2| ^ 6 = η ^ 6 / 64 := by
    rw [show |η / 2| ^ 6 = (|η / 2| ^ 2) ^ 3 by ring, sq_abs]; ring
  rw [e6] at D1
  rw [eh6] at Dh
  have hη6 : 0 ≤ η ^ 6 := by positivity
  have hCpos : 0 ≤ 2 / 3 * (5 * X05.K6 / 4 + 16 * X05.C1) := by positivity
  have hCpos' : 0 ≤ 8 * (18 + 16 * X05.C1 + 17 * X05.K6 / 4) := by positivity
  -- identities
  have id2 : 16 - 15 * (1 - (π ^ 2 / 6) / 4 * η ^ 2 + X05.z 3 / 4 * η ^ 3
        + ((π ^ 2 / 6) ^ 2 / 32 - 7 * (π ^ 4 / 90) / 32) * η ^ 4
        + (-(π ^ 2 / 6) * X05.z 3 / 16 + 3 * X05.z 5 / 16) * η ^ 5)
      - (1 + X05.cf (3 / 2) (-1 / 2) 1 0 * η ^ 2 + X05.cf (3 / 2) (-1 / 2) 1 1 * η ^ 3
        + (X05.cf (3 / 2) (-1 / 2) 1 2 + X05.cf (3 / 2) (-1 / 2) 1 0 ^ 2 / 2) * η ^ 4
        + (X05.cf (3 / 2) (-1 / 2) 1 3 + X05.cf (3 / 2) (-1 / 2) 1 0 * X05.cf (3 / 2) (-1 / 2) 1 1)
          * η ^ 5)
      = 3 / 2 * η ^ 2 * (2 * (π ^ 2 / 6) - 2 * X05.z 3 * η
        + (-(π ^ 2 / 6) ^ 2 / 2 + 3 * (π ^ 4 / 90) / 2) * η ^ 2
        + ((π ^ 2 / 6) * X05.z 3 - X05.z 5) * η ^ 3) := by
    norm_num [X05.cf, X05.z_two, X05.z_four]; ring
  have id4 : 4 * (1 + X05.cf (3 / 2) (-1 / 2) 1 0 * η ^ 2 + X05.cf (3 / 2) (-1 / 2) 1 1 * η ^ 3
        + (X05.cf (3 / 2) (-1 / 2) 1 2 + X05.cf (3 / 2) (-1 / 2) 1 0 ^ 2 / 2) * η ^ 4
        + (X05.cf (3 / 2) (-1 / 2) 1 3 + X05.cf (3 / 2) (-1 / 2) 1 0 * X05.cf (3 / 2) (-1 / 2) 1 1)
          * η ^ 5) - 16
      + 12 * (1 - (π ^ 2 / 6) / 4 * η ^ 2 + X05.z 3 / 4 * η ^ 3
        + ((π ^ 2 / 6) ^ 2 / 32 - 7 * (π ^ 4 / 90) / 32) * η ^ 4
        + (-(π ^ 2 / 6) * X05.z 3 / 16 + 3 * X05.z 5 / 16) * η ^ 5)
      = η ^ 4 / 8 * (12 * (π ^ 2 / 6) ^ 2 + 12 * (π ^ 4 / 90))
        - 3 * (X05.z 5 + (π ^ 2 / 6) * X05.z 3) * η ^ 5 := by
    norm_num [X05.cf, X05.z_two, X05.z_four]; ring
  have hwb : |X05.z 5 + (π ^ 2 / 6) * X05.z 3| ≤ 6 := by
    have : |(π ^ 2 / 6) * X05.z 3| ≤ 4 := by
      rw [abs_mul, abs_of_nonneg (by positivity)]
      nlinarith [abs_nonneg (X05.z 3)]
    have := abs_add_le (X05.z 5) ((π ^ 2 / 6) * X05.z 3)
    linarith
  refine ⟨E0.trans (by rw [e6]; exact mul_le_mul_of_nonneg_right (by linarith) hη6), ?_, ?_⟩
  · rcases eq_or_ne η 0 with h0 | h0
    · subst h0
      rw [(X05.M_zero).1]; simp
    exact (X05.alg_b h0 id2 E0 E1 D1 Dh).trans
      (mul_le_mul_of_nonneg_right (by linarith) (pow_nonneg (abs_nonneg η) 4))
  · rcases eq_or_ne η 0 with h0 | h0
    · subst h0
      rw [(X05.M_zero).2]; simp
    exact (X05.alg_c h0 hη1 hC1.le hK hwb id4 E0 E1 D1 Dh).trans
      (mul_le_mul_of_nonneg_right (by linarith) (abs_nonneg η))


namespace X05

lemma alg_c2 {a0 a1 M2 M4 T0 T1 pc w η C K : ℝ} (hη : η ≠ 0) (_hη1 : |η| ≤ 1)
    (hid : 4 * T1 - 16 + 12 * T0 = η ^ 4 / 8 * pc - 3 * w * η ^ 5)
    (E0 : |a0 - T0| ≤ C * η ^ 6) (E1 : |a1 - T1| ≤ C * η ^ 6)
    (D1 : |a1 - (a0 + η ^ 2 / 2 * M2 + η ^ 4 / 24 * M4)| ≤ K * η ^ 6)
    (Dh : |1 - (a0 + (η / 2) ^ 2 / 2 * M2 + (η / 2) ^ 4 / 24 * M4)| ≤ K * (η ^ 6 / 64)) :
    |M4 - (pc - 24 * w * η)| ≤ 8 * (16 * C + 17 * K / 4) * η ^ 2 := by
  have hpos : 0 < η ^ 4 / 8 := by positivity
  refine div_bound hpos ?_
  have e : η ^ 4 / 8 * (M4 - (pc - 24 * w * η)) = 4 * (a1 - T1) + 12 * (a0 - T0)
      - 4 * (a1 - (a0 + η ^ 2 / 2 * M2 + η ^ 4 / 24 * M4))
      + 16 * (1 - (a0 + (η / 2) ^ 2 / 2 * M2 + (η / 2) ^ 4 / 24 * M4)) := by
    linear_combination hid
  rw [e]
  generalize a1 - (a0 + η ^ 2 / 2 * M2 + η ^ 4 / 24 * M4) = d1 at D1
  generalize 1 - (a0 + (η / 2) ^ 2 / 2 * M2 + (η / 2) ^ 4 / 24 * M4) = dh at Dh
  generalize a1 - T1 = f1 at E1
  generalize a0 - T0 = f0 at E0
  have : η ^ 4 / 8 * (8 * (16 * C + 17 * K / 4) * η ^ 2) =
      4 * (C * η ^ 6) + 12 * (C * η ^ 6) + 4 * (K * η ^ 6) + 16 * (K * (η ^ 6 / 64)) := by ring
  rw [this]
  rw [abs_le] at D1 Dh E1 E0 ⊢
  constructor <;> linarith [D1.1, D1.2, Dh.1, Dh.2, E1.1, E1.2, E0.1, E0.2]

/-- `∂_s⁴A(η,0)` to second order (needed for the `τ⁵` coefficient). -/
lemma M4_refined {η : ℝ} (hη : |η| ≤ 1 / 20) :
    |iteratedDeriv 4 (fun s => A η s) 0 - (12 * (π ^ 2 / 6) ^ 2 + 12 * (π ^ 4 / 90)
      - 24 * (z 5 + (π ^ 2 / 6) * z 3) * η)| ≤ 8 * (16 * C1 + 17 * K6 / 4) * η ^ 2 := by
  rcases eq_or_ne η 0 with h0 | h0
  · subst h0; rw [M_zero.2]; simp
  have hη4 : |η| ≤ 1 / 4 := by linarith
  have hη8 : |η| ≤ 1 / 8 := by linarith
  have hη1 : |η| ≤ 1 := by linarith
  have hη2' : |η / 2| ≤ 1 / 8 := by rw [abs_div]; norm_num; linarith
  have E0 := expand0 hη
  have E1 := expand1 hη
  have Eh := expand_half hη4
  have D1 := s_even hη4 hη8
  have Dh := s_even hη4 hη2'
  change _ ≤ K6 * _ at D1 Dh
  rw [iteratedDeriv_zero] at D1 Dh
  rw [Eh] at Dh
  have e6 : |η| ^ 6 = η ^ 6 := by rw [show |η| ^ 6 = (|η| ^ 2) ^ 3 by ring, sq_abs]; ring
  have eh6 : |η / 2| ^ 6 = η ^ 6 / 64 := by
    rw [show |η / 2| ^ 6 = (|η / 2| ^ 2) ^ 3 by ring, sq_abs]; ring
  rw [e6] at D1
  rw [eh6] at Dh
  have id4 : 4 * (1 + cf (3 / 2) (-1 / 2) 1 0 * η ^ 2 + cf (3 / 2) (-1 / 2) 1 1 * η ^ 3
        + (cf (3 / 2) (-1 / 2) 1 2 + cf (3 / 2) (-1 / 2) 1 0 ^ 2 / 2) * η ^ 4
        + (cf (3 / 2) (-1 / 2) 1 3 + cf (3 / 2) (-1 / 2) 1 0 * cf (3 / 2) (-1 / 2) 1 1)
          * η ^ 5) - 16
      + 12 * (1 - (π ^ 2 / 6) / 4 * η ^ 2 + z 3 / 4 * η ^ 3
        + ((π ^ 2 / 6) ^ 2 / 32 - 7 * (π ^ 4 / 90) / 32) * η ^ 4
        + (-(π ^ 2 / 6) * z 3 / 16 + 3 * z 5 / 16) * η ^ 5)
      = η ^ 4 / 8 * (12 * (π ^ 2 / 6) ^ 2 + 12 * (π ^ 4 / 90))
        - 3 * (z 5 + (π ^ 2 / 6) * z 3) * η ^ 5 := by
    norm_num [cf, z_two, z_four]; ring
  exact alg_c2 h0 hη1 id4 E0 E1 D1 Dh

end X05

namespace X05

lemma alg7 {X a M2 M4 η C5 C6 C7 z3 z5 c : ℝ}
    (h5 : |X - (a + 1 / 4 * η ^ 2 * M2 - 1 / 96 * η ^ 4 * M4)| ≤ C5 * η ^ 6)
    (ha : |a - (1 - c ^ 2 / 6 / 4 * η ^ 2 + z3 / 4 * η ^ 3 +
          ((c ^ 2 / 6) ^ 2 / 32 - 7 * (c ^ 4 / 90) / 32) * η ^ 4 +
          (-(c ^ 2 / 6) * z3 / 16 + 3 * z5 / 16) * η ^ 5)| ≤ C6 * η ^ 6)
    (hb : |M2 - (2 * (c ^ 2 / 6) - 2 * z3 * η +
          (-(c ^ 2 / 6) ^ 2 / 2 + 3 * (c ^ 4 / 90) / 2) * η ^ 2 +
          (c ^ 2 / 6 * z3 - z5) * η ^ 3)| ≤ C6 * η ^ 4)
    (hc : |M4 - (12 * (c ^ 2 / 6) ^ 2 + 12 * (c ^ 4 / 90) -
          24 * (z5 + c ^ 2 / 6 * z3) * η)| ≤ C7 * η ^ 2) :
    |4 * X - (4 + c ^ 2 / 6 * (-η) ^ 2 + z3 * (-η) ^ 3 - 33 / 16 * (c ^ 4 / 90) * (-η) ^ 4
      - (3 / 4 * z5 + 7 / 4 * (c ^ 2 / 6) * z3) * (-η) ^ 5)| ≤
      (4 * C5 + 4 * C6 + C6 + C7 / 24) * η ^ 6 := by
  have key : 4 * X - (4 + c ^ 2 / 6 * (-η) ^ 2 + z3 * (-η) ^ 3 - 33 / 16 * (c ^ 4 / 90) * (-η) ^ 4
      - (3 / 4 * z5 + 7 / 4 * (c ^ 2 / 6) * z3) * (-η) ^ 5) =
      4 * (X - (a + 1 / 4 * η ^ 2 * M2 - 1 / 96 * η ^ 4 * M4))
      + 4 * (a - (1 - c ^ 2 / 6 / 4 * η ^ 2 + z3 / 4 * η ^ 3 +
          ((c ^ 2 / 6) ^ 2 / 32 - 7 * (c ^ 4 / 90) / 32) * η ^ 4 +
          (-(c ^ 2 / 6) * z3 / 16 + 3 * z5 / 16) * η ^ 5))
      + η ^ 2 * (M2 - (2 * (c ^ 2 / 6) - 2 * z3 * η +
          (-(c ^ 2 / 6) ^ 2 / 2 + 3 * (c ^ 4 / 90) / 2) * η ^ 2 +
          (c ^ 2 / 6 * z3 - z5) * η ^ 3))
      - η ^ 4 * (M4 - (12 * (c ^ 2 / 6) ^ 2 + 12 * (c ^ 4 / 90) -
          24 * (z5 + c ^ 2 / 6 * z3) * η)) / 24 := by ring
  rw [key]
  have hη2 : 0 ≤ η ^ 2 := sq_nonneg η
  have hη4 : 0 ≤ η ^ 4 := by positivity
  generalize X - (a + 1 / 4 * η ^ 2 * M2 - 1 / 96 * η ^ 4 * M4) = u1 at h5
  generalize a - _ = u2 at ha
  generalize M2 - _ = u3 at hb
  generalize M4 - _ = u4 at hc
  have hb' : |η ^ 2 * u3| ≤ C6 * η ^ 6 := by
    rw [abs_mul, abs_of_nonneg hη2]
    calc η ^ 2 * |u3| ≤ η ^ 2 * (C6 * η ^ 4) := mul_le_mul_of_nonneg_left hb hη2
      _ = _ := by ring
  have hc' : |η ^ 4 * u4| ≤ C7 * η ^ 6 := by
    rw [abs_mul, abs_of_nonneg hη4]
    calc η ^ 4 * |u4| ≤ η ^ 4 * (C7 * η ^ 2) := mul_le_mul_of_nonneg_left hc hη4
      _ = _ := by ring
  generalize η ^ 2 * u3 = v3 at hb'
  generalize η ^ 4 * u4 = v4 at hc'
  rw [abs_le] at h5 ha hb' hc' ⊢
  constructor <;> linarith [h5.1, h5.2, ha.1, ha.2, hb'.1, hb'.2, hc'.1, hc'.2]

end X05

theorem X05.proof_7 : X05.stmt_7 := by
  obtain ⟨C5, hC5⟩ := X05.proof_5
  obtain ⟨C6, δ6, hδ6, hC6⟩ := X05.proof_6
  refine ⟨4 * C5 + 4 * C6 + C6 + 8 * (16 * X05.C1 + 17 * X05.K6 / 4) / 24, min (1 / 20) δ6,
    lt_min (by norm_num) hδ6, fun p hp hτ => ?_⟩
  have hτη : 1 - 1 / p = -(1 / p - 1) := by ring
  rw [hτη, abs_neg] at hτ ⊢
  have hη20 : |1 / p - 1| ≤ 1 / 20 := hτ.trans (min_le_left _ _)
  have hηδ : |1 / p - 1| ≤ δ6 := hτ.trans (min_le_right _ _)
  have hpη : 1 / (1 + (1 / p - 1)) = p := by field_simp; ring
  have h5 := hC5 (1 / p - 1) (by linarith)
  rw [hpη] at h5
  obtain ⟨ha, hb, -⟩ := hC6 (1 / p - 1) hηδ
  have hc := X05.M4_refined hη20
  have e4 : |1 / p - 1| ^ 4 = (1 / p - 1) ^ 4 := by
    rw [show |1 / p - 1| ^ 4 = (|1 / p - 1| ^ 2) ^ 2 by ring, sq_abs]; ring
  have e6 : |1 / p - 1| ^ 6 = (1 / p - 1) ^ 6 := by
    rw [show |1 / p - 1| ^ 6 = (|1 / p - 1| ^ 2) ^ 3 by ring, sq_abs]; ring
  rw [e6] at h5 ha ⊢
  rw [e4] at hb
  have hs2 : Real.sqrt 2 ≠ 0 := by positivity
  have hP : X05.P p / Real.sqrt 2 = 4 * (X05.P p / (4 * Real.sqrt 2)) := by field_simp
  rw [hP]
  exact X05.alg7 h5 ha hb hc


namespace X05

/-- `√cosh` -/
noncomputable def hc (x : ℝ) : ℝ := Real.sqrt (Real.cosh x)

lemma hc_contDiff (n : ℕ) : ContDiff ℝ n hc :=
  Real.contDiff_cosh.sqrt fun x => (Real.cosh_pos x).ne'

lemma hc_pos (x : ℝ) : 0 < hc x := Real.sqrt_pos.2 (Real.cosh_pos x)

lemma tanh_hasDerivAt (x : ℝ) : HasDerivAt Real.tanh (1 - Real.tanh x ^ 2) x := by
  have hc0 := (Real.cosh_pos x).ne'
  have h := (Real.hasDerivAt_sinh x).div (Real.hasDerivAt_cosh x) hc0
  have e : Real.sinh / Real.cosh = Real.tanh := by
    ext y; rw [Pi.div_apply, Real.tanh_eq_sinh_div_cosh]
  rw [e] at h
  convert h using 1
  rw [Real.tanh_eq_sinh_div_cosh]
  have := Real.cosh_sq_sub_sinh_sq x
  field_simp

lemma hc_hasDerivAt (x : ℝ) : HasDerivAt hc (Real.tanh x / 2 * hc x) x := by
  have hc0 := (Real.cosh_pos x).ne'
  have h := (Real.hasDerivAt_cosh x).sqrt hc0
  unfold hc
  convert h using 1
  rw [Real.tanh_eq_sinh_div_cosh]
  have hs : Real.sqrt (Real.cosh x) ^ 2 = Real.cosh x := Real.sq_sqrt (Real.cosh_pos x).le
  have hs0 : Real.sqrt (Real.cosh x) ≠ 0 := (Real.sqrt_pos.2 (Real.cosh_pos x)).ne'
  field_simp
  rw [hs]

lemma hc_struct (k : ℕ) : ∃ Q : Polynomial ℝ, ∀ x, iteratedDeriv k hc x =
    Q.eval (Real.tanh x) * hc x := by
  induction k with
  | zero => exact ⟨1, fun x => by simp⟩
  | succ k ih =>
    obtain ⟨Q, hQ⟩ := ih
    refine ⟨Q.derivative * (1 - Polynomial.X ^ 2) + Q * (Polynomial.C (1 / 2) * Polynomial.X),
      fun x => ?_⟩
    rw [iteratedDeriv_succ, show iteratedDeriv k hc = fun x => Q.eval (Real.tanh x) * hc x from
      funext hQ]
    have h1 := ((Q.hasDerivAt (Real.tanh x)).comp x (tanh_hasDerivAt x)).mul (hc_hasDerivAt x)
    have h2 : HasDerivAt (fun y => Q.eval (Real.tanh y) * hc y) _ x := h1
    rw [h2.deriv]
    simp only [Polynomial.eval_add, Polynomial.eval_mul, Polynomial.eval_sub, Polynomial.eval_pow,
      Polynomial.eval_X, Polynomial.eval_C, Polynomial.eval_one, Function.comp]
    ring

lemma hc_deriv_bd (k : ℕ) : ∃ C : ℝ, 0 ≤ C ∧ ∀ x, |iteratedDeriv k hc x| ≤ C * hc x := by
  obtain ⟨Q, hQ⟩ := hc_struct k
  obtain ⟨C, hC⟩ := (isCompact_Icc (a := (-1 : ℝ)) (b := 1)).exists_bound_of_continuousOn
    (Q.continuous.continuousOn)
  refine ⟨max C 0, le_max_right _ _, fun x => ?_⟩
  rw [hQ, abs_mul, abs_of_pos (hc_pos x)]
  have hx : Real.tanh x ∈ Set.Icc (-1 : ℝ) 1 :=
    ⟨(Real.neg_one_lt_tanh x).le, (Real.tanh_lt_one x).le⟩
  have := hC _ hx
  rw [Real.norm_eq_abs] at this
  exact mul_le_mul_of_nonneg_right (this.trans (le_max_left _ _)) (hc_pos x).le

/-- a monotone choice of the constants -/
lemma hc_deriv_bd' (N : ℕ) : ∃ C : ℝ, 0 ≤ C ∧ ∀ k ≤ N, ∀ x, |iteratedDeriv k hc x| ≤ C * hc x := by
  induction N with
  | zero =>
    obtain ⟨C, hC0, hC⟩ := hc_deriv_bd 0
    exact ⟨C, hC0, fun k hk x => by rw [Nat.le_zero.1 hk]; exact hC x⟩
  | succ N ih =>
    obtain ⟨C, hC0, hC⟩ := ih
    obtain ⟨D, hD0, hD⟩ := hc_deriv_bd (N + 1)
    refine ⟨C + D, by positivity, fun k hk x => ?_⟩
    have hp := hc_pos x
    rcases Nat.lt_or_ge k (N + 1) with h | h
    · have := hC k (by omega) x; nlinarith
    · have := hD x; rw [show k = N + 1 by omega]; nlinarith

end X05

namespace X05

/-- `log(t(1-t))/2` -/
noncomputable def la (t : ℝ) : ℝ := Real.log (t * (1 - t)) / 2

/-- the building blocks `a^j l^k e^{η a} (√cosh)^{(k)}(η l)` of the `η`-derivatives -/
noncomputable def φ (j k : ℕ) (η t : ℝ) : ℝ :=
  la t ^ j * logit t ^ k * Real.exp (η * la t) * iteratedDeriv k hc (η * logit t)

lemma φ_hasDerivAt (j k : ℕ) (η t : ℝ) :
    HasDerivAt (fun η => φ j k η t) (φ (j + 1) k η t + φ j (k + 1) η t) η := by
  have hD : ∀ y, HasDerivAt (iteratedDeriv k hc) (iteratedDeriv (k + 1) hc y) y := by
    intro y
    have := ((hc_contDiff (k + 1)).differentiable_iteratedDeriv' k y).hasDerivAt
    rwa [← iteratedDeriv_succ] at this
  have h1 := (hD (η * logit t)).comp η ((hasDerivAt_id' η).mul_const (logit t))
  have h2 := ((hasDerivAt_id' η).mul_const (la t)).exp
  have h3 := (h2.mul h1).const_mul (la t ^ j * logit t ^ k)
  unfold φ
  convert h3 using 1
  · ext y; simp only [Pi.mul_apply, Function.comp_apply]; ring
  · simp only [Function.comp_apply]; ring

lemma log_bd_eps {x ε : ℝ} (hx0 : 0 < x) (hx1 : x ≤ 1) (hε : 0 < ε) :
    |Real.log x| ≤ x ^ (-ε) / ε := by
  have h := Real.neg_rpow_div_le_log hx0.le hε
  have hl : Real.log x ≤ 0 := Real.log_nonpos hx0.le hx1
  rw [abs_of_nonpos hl]
  have : -x ^ (-ε) / ε = -(x ^ (-ε) / ε) := by ring
  linarith

lemma al_bd (N : ℕ) : ∃ K : ℝ, 0 ≤ K ∧ ∀ t : ℝ, 0 < t → t < 1 →
    (|la t| + |logit t|) ^ N ≤ K * (t * (1 - t)) ^ (-(1 / 8 : ℝ)) := by
  set ε : ℝ := 1 / (8 * ((N : ℝ) + 1)) with hε
  have hε0 : 0 < ε := by positivity
  refine ⟨(3 / ε) ^ N, by positivity, fun t ht0 ht1 => ?_⟩
  have h1 : 0 < 1 - t := by linarith
  have hy : 0 < t * (1 - t) := mul_pos ht0 h1
  have hy1 : t * (1 - t) ≤ 1 := by nlinarith
  have a := log_bd_eps ht0 ht1.le hε0
  have b := log_bd_eps h1 (by linarith) hε0
  have a' : t ^ (-ε) ≤ (t * (1 - t)) ^ (-ε) :=
    Real.rpow_le_rpow_of_nonpos hy (by nlinarith) (by linarith)
  have b' : (1 - t) ^ (-ε) ≤ (t * (1 - t)) ^ (-ε) :=
    Real.rpow_le_rpow_of_nonpos hy (by nlinarith) (by linarith)
  have hla : |la t| ≤ (|Real.log t| + |Real.log (1 - t)|) / 2 := by
    unfold la
    rw [Real.log_mul ht0.ne' h1.ne', abs_div, abs_of_pos (by norm_num : (0:ℝ) < 2)]
    gcongr
    exact abs_add_le _ _
  have hl : |logit t| ≤ |Real.log t| + |Real.log (1 - t)| := by
    unfold logit; rw [Real.log_div ht0.ne' h1.ne']; exact abs_sub _ _
  have hsum : |la t| + |logit t| ≤ 3 / ε * (t * (1 - t)) ^ (-ε) := by
    have hdiv : (t * (1 - t)) ^ (-ε) / ε + (t * (1 - t)) ^ (-ε) / ε =
        2 / ε * (t * (1 - t)) ^ (-ε) := by ring
    have : 0 ≤ (t * (1 - t)) ^ (-ε) := by positivity
    have : 3 / ε * (t * (1 - t)) ^ (-ε) = 3 * ((t * (1 - t)) ^ (-ε) / ε) := by ring
    have ha'' : |Real.log t| ≤ (t * (1 - t)) ^ (-ε) / ε :=
      a.trans (div_le_div_of_nonneg_right a' hε0.le)
    have hb'' : |Real.log (1 - t)| ≤ (t * (1 - t)) ^ (-ε) / ε :=
      b.trans (div_le_div_of_nonneg_right b' hε0.le)
    linarith
  have hpos : 0 ≤ |la t| + |logit t| := by positivity
  calc (|la t| + |logit t|) ^ N ≤ (3 / ε * (t * (1 - t)) ^ (-ε)) ^ N :=
        pow_le_pow_left₀ hpos hsum N
    _ = (3 / ε) ^ N * (t * (1 - t)) ^ (-ε * N) := by
        rw [mul_pow, ← Real.rpow_mul_natCast hy.le]
    _ ≤ (3 / ε) ^ N * (t * (1 - t)) ^ (-(1 / 8 : ℝ)) := by
        refine mul_le_mul_of_nonneg_left (Real.rpow_le_rpow_of_exponent_ge hy hy1 ?_) (by positivity)
        rw [hε]
        have hN : (0 : ℝ) ≤ N := N.cast_nonneg
        rw [neg_mul, neg_le_neg_iff, div_mul_eq_mul_div, one_mul, div_le_div_iff₀ (by positivity)
          (by norm_num)]
        nlinarith

end X05

namespace X05

lemma exp_hc_bd {η t : ℝ} (hη : |η| ≤ 1 / 4) (ht0 : 0 < t) (ht1 : t < 1) :
    Real.exp (η * la t) * hc (η * logit t) ≤ 2 * (t * (1 - t)) ^ (-(1 / 2 : ℝ)) := by
  have h1 : 0 < 1 - t := by linarith
  have hy : 0 < t * (1 - t) := mul_pos ht0 h1
  have hw : Real.exp (η * la t) = (t * (1 - t)) ^ (η / 2) := by
    rw [Real.rpow_def_of_pos hy, la]; ring_nf
  rw [hw]
  have hs := sqrt_le_exp_half (η * logit t)
  have he : Real.exp (|η * logit t| / 2) ≤
      Real.exp (η / 2 * logit t) + Real.exp (-(η / 2) * logit t) := by
    rcases le_total 0 (η * logit t) with h | h
    · rw [abs_of_nonneg h]
      have : Real.exp (η * logit t / 2) = Real.exp (η / 2 * logit t) := by ring_nf
      rw [this]; linarith [Real.exp_pos (-(η / 2) * logit t)]
    · rw [abs_of_nonpos h]
      have : Real.exp (-(η * logit t) / 2) = Real.exp (-(η / 2) * logit t) := by ring_nf
      rw [this]; linarith [Real.exp_pos (η / 2 * logit t)]
  have p1 := pointwise2 (η := η) (s := η / 2) ht0 ht1
  have p2 := pointwise2 (η := η) (s := -(η / 2)) ht0 ht1
  have hb1 : |η / 2 + η / 2| ≤ 1 / 4 := by rw [show η / 2 + η / 2 = η by ring]; exact hη
  have hb2 : |η / 2 - η / 2| ≤ 1 / 4 := by rw [sub_self, abs_zero]; norm_num
  have hb3 : |η / 2 + -(η / 2)| ≤ 1 / 4 := by rw [add_neg_cancel, abs_zero]; norm_num
  have hb4 : |η / 2 - -(η / 2)| ≤ 1 / 4 := by rw [sub_neg_eq_add, show η / 2 + η / 2 = η by ring]; exact hη
  have q1 := pow_pair_bd ht0 ht1 hb1 hb2
  have q2 := pow_pair_bd ht0 ht1 hb3 hb4
  have hw0 : 0 ≤ (t * (1 - t)) ^ (η / 2) := by positivity
  calc (t * (1 - t)) ^ (η / 2) * hc (η * logit t)
      ≤ (t * (1 - t)) ^ (η / 2) * (Real.exp (η / 2 * logit t) + Real.exp (-(η / 2) * logit t)) :=
        mul_le_mul_of_nonneg_left (hs.trans he) hw0
    _ = (t * (1 - t)) ^ (η / 2) * Real.exp (η / 2 * logit t) +
        (t * (1 - t)) ^ (η / 2) * Real.exp (-(η / 2) * logit t) := by ring
    _ ≤ _ := by rw [p1, p2]; linarith

lemma φ_bd (j k : ℕ) : ∃ B : ℝ, 0 ≤ B ∧ ∀ η t : ℝ, |η| ≤ 1 / 4 → 0 < t → t < 1 →
    |φ j k η t| ≤ B * q t := by
  obtain ⟨C, hC0, hC⟩ := hc_deriv_bd k
  obtain ⟨K, hK0, hK⟩ := al_bd (j + k)
  refine ⟨C * K * 2, by positivity, fun η t hη ht0 ht1 => ?_⟩
  have h1 : 0 < 1 - t := by linarith
  have hy : 0 < t * (1 - t) := mul_pos ht0 h1
  unfold φ
  rw [abs_mul, abs_mul, abs_mul, abs_pow, abs_pow, abs_of_pos (Real.exp_pos _)]
  have e1 := hC (η * logit t)
  have e2 := exp_hc_bd hη ht0 ht1
  have e3 := hK t ht0 ht1
  have hal : |la t| ^ j * |logit t| ^ k ≤ (|la t| + |logit t|) ^ (j + k) := by
    rw [pow_add]
    exact mul_le_mul (pow_le_pow_left₀ (abs_nonneg _) (by linarith [abs_nonneg (logit t)]) _)
      (pow_le_pow_left₀ (abs_nonneg _) (by linarith [abs_nonneg (la t)]) _) (by positivity)
      (by positivity)
  have hexp : Real.exp (η * la t) * |iteratedDeriv k hc (η * logit t)| ≤
      C * (2 * (t * (1 - t)) ^ (-(1 / 2 : ℝ))) := by
    calc Real.exp (η * la t) * |iteratedDeriv k hc (η * logit t)|
        ≤ Real.exp (η * la t) * (C * hc (η * logit t)) :=
          mul_le_mul_of_nonneg_left e1 (Real.exp_pos _).le
      _ = C * (Real.exp (η * la t) * hc (η * logit t)) := by ring
      _ ≤ _ := mul_le_mul_of_nonneg_left e2 hC0
  calc |la t| ^ j * |logit t| ^ k * Real.exp (η * la t) * |iteratedDeriv k hc (η * logit t)|
      = (|la t| ^ j * |logit t| ^ k) * (Real.exp (η * la t) * |iteratedDeriv k hc (η * logit t)|) := by
        ring
    _ ≤ (K * (t * (1 - t)) ^ (-(1 / 8 : ℝ))) * (C * (2 * (t * (1 - t)) ^ (-(1 / 2 : ℝ)))) :=
        mul_le_mul (hal.trans e3) hexp (by positivity) (by positivity)
    _ = C * K * 2 * q t := by
        rw [q, show C * K * 2 * (t * (1 - t)) ^ (-(5 / 8 : ℝ)) =
          C * K * 2 * ((t * (1 - t)) ^ (-(1 / 8 : ℝ)) * (t * (1 - t)) ^ (-(1 / 2 : ℝ))) by
          rw [← Real.rpow_add hy]; norm_num]
        ring

lemma φ_meas (j k : ℕ) (η : ℝ) : Measurable (fun t => φ j k η t) := by
  have hc' : Measurable (iteratedDeriv k hc) :=
    ((hc_contDiff k).continuous_iteratedDeriv' k).measurable
  have hl : Measurable logit := by unfold logit; fun_prop
  have ha : Measurable la := by unfold la; fun_prop
  unfold φ
  exact ((ha.pow_const j).mul (hl.pow_const k)).mul ((measurable_const.mul ha).exp) |>.mul
    (hc'.comp (measurable_const.mul hl))

lemma φ_int (j k : ℕ) {η : ℝ} (hη : |η| ≤ 1 / 4) :
    IntervalIntegrable (fun t => φ j k η t) volume 0 1 := by
  obtain ⟨B, hB0, hB⟩ := φ_bd j k
  refine IntervalIntegrable.mono_fun' (q_int.const_mul B) (φ_meas j k η).aestronglyMeasurable ?_
  rw [uIoc_of_le zero_le_one, ← restrict_Ioo_eq_restrict_Ioc]
  exact ae_restrict_of_forall_mem measurableSet_Ioo fun t ht => by
    dsimp only; rw [Real.norm_eq_abs]; exact hB η t hη ht.1 ht.2

/-- `Φ_{j,k}(η) = ∫ φ_{j,k}(η,t) dt` -/
noncomputable def Φ (j k : ℕ) (η : ℝ) : ℝ := ∫ t in (0 : ℝ)..1, φ j k η t

lemma Φ_hasDerivAt (j k : ℕ) {η₀ : ℝ} (h : |η₀| < 1 / 4) :
    HasDerivAt (Φ j k) (Φ (j + 1) k η₀ + Φ j (k + 1) η₀) η₀ := by
  obtain ⟨B1, hB1, hb1⟩ := φ_bd (j + 1) k
  obtain ⟨B2, hB2, hb2⟩ := φ_bd j (k + 1)
  have hS : Set.Ioo (-(1 / 4) : ℝ) (1 / 4) ∈ 𝓝 η₀ :=
    Ioo_mem_nhds (by linarith [neg_abs_le η₀]) (by linarith [le_abs_self η₀])
  have := intervalIntegral.hasDerivAt_integral_of_dominated_loc_of_deriv_le
    (μ := volume) (a := 0) (b := 1)
    (F := fun η t => φ j k η t) (F' := fun η t => φ (j + 1) k η t + φ j (k + 1) η t)
    (bound := fun t => (B1 + B2) * q t) hS
    (Filter.Eventually.of_forall fun x => (φ_meas j k x).aestronglyMeasurable)
    (φ_int j k h.le) (((φ_meas _ _ _).add (φ_meas _ _ _)).aestronglyMeasurable) ?_
    (q_int.const_mul _) (Filter.Eventually.of_forall fun t _ x _ => φ_hasDerivAt j k x t)
  · rw [intervalIntegral.integral_add (φ_int _ _ h.le) (φ_int _ _ h.le)] at this
    exact this.2
  · filter_upwards [ae_ne_one] with t ht1 htm x hx
    rw [uIoc_of_le zero_le_one] at htm
    have hx' : |x| ≤ 1 / 4 := by rw [abs_le]; exact ⟨hx.1.le, hx.2.le⟩
    have ht1' : t < 1 := lt_of_le_of_ne htm.2 ht1
    have a := hb1 x t hx' htm.1 ht1'
    have b := hb2 x t hx' htm.1 ht1'
    rw [Real.norm_eq_abs]
    calc |φ (j + 1) k x t + φ j (k + 1) x t| ≤ |φ (j + 1) k x t| + |φ j (k + 1) x t| :=
          abs_add_le _ _
      _ ≤ _ := by linarith

lemma Φ_contDiffOn (m : ℕ) : ∀ j k : ℕ, ContDiffOn ℝ m (Φ j k) (Set.Ioo (-(1 / 4)) (1 / 4)) := by
  induction m with
  | zero =>
    intro j k
    rw [Nat.cast_zero, contDiffOn_zero]
    intro x hx
    have hx' : |x| < 1 / 4 := abs_lt.mpr hx
    exact (Φ_hasDerivAt j k hx').continuousAt.continuousWithinAt
  | succ m ih =>
    intro j k
    rw [show ((m + 1 : ℕ) : WithTop ℕ∞) = (m : WithTop ℕ∞) + 1 by push_cast; rfl,
      contDiffOn_succ_iff_deriv_of_isOpen isOpen_Ioo]
    refine ⟨fun x hx => (Φ_hasDerivAt j k (abs_lt.mpr hx)).differentiableAt.differentiableWithinAt,
      fun h => absurd h (by simp), ?_⟩
    refine ((ih (j + 1) k).add (ih j (k + 1))).congr fun x hx => ?_
    exact (Φ_hasDerivAt j k (abs_lt.mpr hx)).deriv

end X05

lemma X05.P_eq_near {p : ℝ} (hp : 0 < p) :
    X05.P p = 4 * Real.sqrt 2 * p⁻¹ * X05.Φ 0 0 (p⁻¹ - 1) := by
  rw [(X05.proof_1 p hp).2]
  have hJ : ∫ t in (0 : ℝ)..1, (t * (1 - t)) ^ ((1 / p - 1) / 2) *
      Real.sqrt (Real.cosh ((1 / p - 1) * X05.logit t)) = X05.Φ 0 0 (p⁻¹ - 1) := by
    unfold X05.Φ
    refine intervalIntegral.integral_congr_uIoo fun t ht => ?_
    rw [uIoo_of_le zero_le_one] at ht
    have hy : 0 < t * (1 - t) := mul_pos ht.1 (by linarith [ht.2])
    simp only [X05.φ, pow_zero, one_mul, iteratedDeriv_zero, X05.hc, X05.la]
    rw [Real.rpow_def_of_pos hy, one_div]
    ring_nf
  rw [hJ, one_div]
  ring

lemma X05.g_contDiffAt (n : ℕ) :
    ContDiffAt ℝ n (fun p : ℝ => 4 * Real.sqrt 2 * p⁻¹ * X05.Φ 0 0 (p⁻¹ - 1)) 1 := by
  have hΦ : ContDiffAt ℝ n (X05.Φ 0 0) 0 :=
    (X05.Φ_contDiffOn n 0 0).contDiffAt (Ioo_mem_nhds (by norm_num) (by norm_num))
  have hinner : ContDiffAt ℝ n (fun p : ℝ => p⁻¹ - 1) 1 :=
    (contDiffAt_inv ℝ one_ne_zero).sub contDiffAt_const
  have hΦ' : ContDiffAt ℝ n (X05.Φ 0 0) ((fun p : ℝ => p⁻¹ - 1) 1) := by
    convert hΦ using 2; norm_num
  have h2 : ContDiffAt ℝ n (fun p : ℝ => X05.Φ 0 0 (p⁻¹ - 1)) 1 :=
    ContDiffAt.comp (g := X05.Φ 0 0) (f := fun p : ℝ => p⁻¹ - 1) 1 hΦ' hinner
  have h3 : ContDiffAt ℝ n (fun p : ℝ => 4 * Real.sqrt 2 * p⁻¹) 1 :=
    contDiffAt_const.mul (contDiffAt_inv ℝ one_ne_zero)
  exact h3.mul h2

theorem X05.proof_8 : X05.stmt_8 := by
  intro n
  refine (X05.g_contDiffAt n).congr_of_eventuallyEq ?_
  filter_upwards [Ioi_mem_nhds (show (0:ℝ) < 1 by norm_num)] with p hp
  exact X05.P_eq_near hp


-- Copied (renamed) from Zeta3/Wave1/T02_SeriesDiamond.lean: Taylor coefficients from the tau-series.
namespace X05T

open Asymptotics Filter

/-- `ζ(2)`, the actual Mathlib value (real part of the complex `riemannZeta 2`). -/
noncomputable def ζ2 : ℝ := (riemannZeta 2).re

/-- `ζ(4)`, the actual Mathlib value (real part of the complex `riemannZeta 4`). -/
noncomputable def ζ4 : ℝ := (riemannZeta 4).re

/-- Line 1181: `τ = d/(1+d)`. -/
noncomputable def tau (d : ℝ) : ℝ := d / (1 + d)

/-- The degree-5 polynomial of `thm:pw`:
`4 + ζ(2)τ² + ζ(3)τ³ - (33/16)ζ(4)τ⁴ - (3/4 ζ(5) + 7/4 ζ(2)ζ(3))τ⁵`,
with `ζ(3) = z3`, `ζ(5) = z5` formal. -/
noncomputable def pw (z3 z5 t : ℝ) : ℝ :=
  4 + ζ2 * t ^ 2 + z3 * t ^ 3 - 33 / 16 * ζ4 * t ^ 4
    - (3 / 4 * z5 + 7 / 4 * ζ2 * z3) * t ^ 5

/-- The coefficients `a_0, …, a_5` of `prop:an`, with `ζ(3) = z3`, `ζ(5) = z5` formal
(the value `0` for `n ≥ 6` is a placeholder, never used). -/
noncomputable def a (z3 z5 : ℝ) : ℕ → ℝ
  | 0 => 4
  | 1 => 0
  | 2 => Real.pi ^ 2 / 6
  | 3 => z3 - Real.pi ^ 2 / 3
  | 4 => -3 * z3 - 11 * Real.pi ^ 4 / 480 + Real.pi ^ 2 / 2
  | 5 => -(7 * Real.pi ^ 2 * z3) / 24 + 6 * z3 - 3 * z5 / 4 - 2 * Real.pi ^ 2 / 3
          + 11 * Real.pi ^ 4 / 120
  | _ => 0

/-- The degree-5 polynomial `Σ_{n=0}^{5} a_n d^n` of `prop:an`. -/
noncomputable def anPoly (z3 z5 d : ℝ) : ℝ :=
  ∑ n ∈ Finset.range 6, a z3 z5 n * d ^ n

/-- Re-expanding the `thm:pw` polynomial in `τ = d/(1+d)` in powers of `d` gives exactly
`a_0, …, a_5` of `prop:an`: the difference is `O(d^6)` as `d → 0`. -/
def st_1 : Prop :=
  ∀ z3 z5 : ℝ,
    (fun d : ℝ => pw z3 z5 (tau d) - anPoly z3 z5 d) =O[nhds 0] (fun d : ℝ => d ^ 6)

/-- Consistency of the two results: if `P(p)/√2 = (thm:pw polynomial in τ) + O(τ^6)` with
`p = 1/(1-τ)`, then `P(1+d)/√2 = Σ_{n≤5} a_n d^n + O(d^6)` with the `a_n` of `prop:an`. -/
def st_2 : Prop :=
  ∀ (z3 z5 : ℝ) (F : ℝ → ℝ),
    (fun t : ℝ => F (1 / (1 - t)) - pw z3 z5 t) =O[nhds 0] (fun t : ℝ => t ^ 6) →
    (fun d : ℝ => F (1 + d) - anPoly z3 z5 d) =O[nhds 0] (fun d : ℝ => d ^ 6)

/-- `prop:an` as literally defined (`a_n = F^{(n)}(1)/n!` with `F = P/√2`) follows from
`thm:pw`: if `F` is `C^5` at `1` and satisfies the `thm:pw` expansion, then
`F^{(n)}(1)/n! = a_n` for `n = 0, …, 5`. -/
def st_3 : Prop :=
  ∀ (z3 z5 : ℝ) (F : ℝ → ℝ), ContDiffAt ℝ 5 F 1 →
    (fun t : ℝ => F (1 / (1 - t)) - pw z3 z5 t) =O[nhds 0] (fun t : ℝ => t ^ 6) →
    ∀ n ≤ 5, iteratedDeriv n F 1 / (n.factorial : ℝ) = a z3 z5 n

end X05T

section X05THelpers
open Asymptotics Filter

noncomputable def X05T.Q (z3 z5 d : ℝ) : ℝ :=
  ((5/6)*Real.pi^2 + (-10)*z3 + (-11/48)*Real.pi^4 + (35/24)*z3*Real.pi^2 + (15/4)*z5)
  + ((19/6)*Real.pi^2 + (-35)*z3 + (-11/16)*Real.pi^4 + (35/12)*z3*Real.pi^2 + (15/2)*z5) * d
  + ((-46)*z3 + (9/2)*Real.pi^2 + (-77/96)*Real.pi^4 + (35/12)*z3*Real.pi^2 + (15/2)*z5) * d^2
  + ((-27)*z3 + (-209/480)*Real.pi^4 + (17/6)*Real.pi^2 + (35/24)*z3*Real.pi^2 + (15/4)*z5) * d^3
  + ((7/24)*z3*Real.pi^2 + (-6)*z3 + (3/4)*z5 + (2/3)*Real.pi^2 + (-11/120)*Real.pi^4) * d^4

lemma X05T.ζ2_eq : X05T.ζ2 = Real.pi ^ 2 / 6 := by
  unfold X05T.ζ2
  rw [riemannZeta_two]
  have : ((Real.pi : ℂ) ^ 2 / 6) = ((Real.pi ^ 2 / 6 : ℝ) : ℂ) := by push_cast; ring
  rw [this, Complex.ofReal_re]

lemma X05T.ζ4_eq : X05T.ζ4 = Real.pi ^ 4 / 90 := by
  unfold X05T.ζ4
  rw [riemannZeta_four]
  have : ((Real.pi : ℂ) ^ 4 / 90) = ((Real.pi ^ 4 / 90 : ℝ) : ℂ) := by push_cast; ring
  rw [this, Complex.ofReal_re]

lemma X05T.key (z3 z5 d : ℝ) (hd : 1 + d ≠ 0) :
    X05T.pw z3 z5 (X05T.tau d) - X05T.anPoly z3 z5 d = d ^ 6 * (X05T.Q z3 z5 d / (1 + d) ^ 5) := by
  unfold X05T.pw X05T.tau X05T.anPoly X05T.Q
  rw [X05T.ζ2_eq, X05T.ζ4_eq]
  simp only [Finset.sum_range_succ, Finset.sum_range_zero, X05T.a]
  field_simp
  ring

lemma X05T.tauO : Filter.Tendsto X05T.tau (nhds 0) (nhds 0) := by
  have : ContinuousAt X05T.tau 0 := by
    unfold X05T.tau
    have hc : ContinuousAt (fun d : ℝ => 1 + d) 0 := by fun_prop
    exact continuousAt_id.div hc (by norm_num : (1:ℝ) + 0 ≠ 0)
  simpa [X05T.tau] using this.tendsto

lemma X05T.ev_ne : ∀ᶠ d in nhds (0:ℝ), 1 + d ≠ 0 :=
  have hc : ContinuousAt (fun d : ℝ => 1 + d) 0 := by fun_prop
  hc.eventually_ne (by norm_num : (1:ℝ) + 0 ≠ 0)

lemma X05T.h1 (z3 z5 : ℝ) :
    (fun d : ℝ => X05T.pw z3 z5 (X05T.tau d) - X05T.anPoly z3 z5 d) =O[nhds 0] (fun d : ℝ => d ^ 6) := by
  have hg : ContinuousAt (fun d : ℝ => X05T.Q z3 z5 d / (1 + d) ^ 5) 0 := by
    have hQ : ContinuousAt (fun d : ℝ => X05T.Q z3 z5 d) 0 :=
      (by unfold X05T.Q; fun_prop : Continuous (fun d : ℝ => X05T.Q z3 z5 d)).continuousAt
    have hden : ContinuousAt (fun d : ℝ => (1 + d) ^ 5) 0 := by fun_prop
    exact hQ.div hden (by simp)
  have hb := hg.isBigO (F := ℝ)
  have h2 : (fun d : ℝ => d ^ 6 * (X05T.Q z3 z5 d / (1 + d) ^ 5)) =O[nhds 0] (fun d : ℝ => d ^ 6) := by
    simpa using (Asymptotics.isBigO_refl (fun d : ℝ => d ^ 6) (nhds 0)).mul hb
  refine h2.congr' ?_ Filter.EventuallyEq.rfl
  filter_upwards [X05T.ev_ne] with d hd
  exact (X05T.key z3 z5 d hd).symm

lemma X05T.poly_zero : ∀ (m : ℕ) (b : ℕ → ℝ),
    (fun d : ℝ => ∑ k ∈ Finset.range (m + 1), b k * d ^ k) =o[nhdsWithin 0 {0}ᶜ] (fun d : ℝ => d ^ m) →
    ∀ k ≤ m, b k = 0 := by
  intro m
  induction m with
  | zero =>
    intro b h k hk
    have hk0 : k = 0 := by omega
    subst hk0
    simp only [zero_add, Finset.range_one, Finset.sum_singleton, pow_zero, mul_one] at h
    simpa using h
  | succ m ih =>
    intro b h k hk
    have hcont : Continuous (fun d : ℝ => ∑ k ∈ Finset.range (m + 1 + 1), b k * d ^ k) := by
      fun_prop
    have hf0 : (∑ k ∈ Finset.range (m + 1 + 1), b k * (0:ℝ) ^ k) = b 0 := by
      rw [Finset.sum_range_succ']; simp
    have ht1 : Filter.Tendsto (fun d : ℝ => ∑ k ∈ Finset.range (m + 1 + 1), b k * d ^ k)
        (nhdsWithin 0 {0}ᶜ) (nhds (b 0)) := by
      have := (hcont.tendsto 0).mono_left (nhdsWithin_le_nhds (s := ({0}ᶜ : Set ℝ)))
      rwa [hf0] at this
    have hpow : Filter.Tendsto (fun d : ℝ => d ^ (m + 1)) (nhdsWithin 0 {0}ᶜ) (nhds 0) := by
      have := ((continuous_pow (m + 1)).tendsto (0:ℝ)).mono_left
        (nhdsWithin_le_nhds (s := ({0}ᶜ : Set ℝ)))
      simpa using this
    have ht2 := h.trans_tendsto hpow
    have hb0 : b 0 = 0 := tendsto_nhds_unique ht1 ht2
    have hsum : ∀ d : ℝ, ∑ k ∈ Finset.range (m + 1 + 1), b k * d ^ k
        = d * ∑ k ∈ Finset.range (m + 1), b (k + 1) * d ^ k := by
      intro d
      rw [Finset.sum_range_succ', hb0, Finset.mul_sum]
      simp only [zero_mul, add_zero]
      exact Finset.sum_congr rfl (fun k _ => by ring)
    have h' : (fun d : ℝ => ∑ k ∈ Finset.range (m + 1), b (k + 1) * d ^ k)
        =o[nhdsWithin 0 {0}ᶜ] (fun d : ℝ => d ^ m) := by
      rw [Asymptotics.isLittleO_iff]
      intro c hc
      have h1 := h.def hc
      have hne : ∀ᶠ d in nhdsWithin (0:ℝ) {0}ᶜ, d ≠ 0 := self_mem_nhdsWithin
      filter_upwards [h1, hne] with d hd hd0
      rw [hsum d, norm_mul, pow_succ, norm_mul] at hd
      have hpos : 0 < ‖d‖ := norm_pos_iff.mpr hd0
      have : ‖d‖ * ‖∑ k ∈ Finset.range (m + 1), b (k + 1) * d ^ k‖ ≤ ‖d‖ * (c * ‖d ^ m‖) := by
        calc _ ≤ c * (‖d ^ m‖ * ‖d‖) := hd
          _ = _ := by ring
      exact le_of_mul_le_mul_left this hpos
    have := ih (fun k => b (k + 1)) h'
    rcases Nat.eq_zero_or_pos k with rfl | hk0
    · exact hb0
    · obtain ⟨j, rfl⟩ := Nat.exists_eq_succ_of_ne_zero hk0.ne'
      exact this j (by omega)

lemma X05T.pf_1 : X05T.st_1 := by
  intro z3 z5
  exact X05T.h1 z3 z5

lemma X05T.ev_eq (F : ℝ → ℝ) (d : ℝ) (hd : 1 + d ≠ 0) : F (1 / (1 - X05T.tau d)) = F (1 + d) := by
  congr 1
  unfold X05T.tau
  field_simp
  ring

lemma X05T.tau6 : (fun d : ℝ => X05T.tau d ^ 6) =O[nhds 0] (fun d : ℝ => d ^ 6) := by
  have hg : ContinuousAt (fun d : ℝ => ((1 + d)⁻¹) ^ 6) 0 := by
    have hc : ContinuousAt (fun d : ℝ => 1 + d) 0 := by fun_prop
    have hi : ContinuousAt (fun d : ℝ => (1 + d)⁻¹) 0 := hc.inv₀ (by simp)
    exact hi.pow 6
  have hb := hg.isBigO (F := ℝ)
  have h2 : (fun d : ℝ => d ^ 6 * ((1 + d)⁻¹) ^ 6) =O[nhds 0] (fun d : ℝ => d ^ 6) := by
    simpa using (Asymptotics.isBigO_refl (fun d : ℝ => d ^ 6) (nhds 0)).mul hb
  refine h2.congr' ?_ Filter.EventuallyEq.rfl
  filter_upwards with d
  unfold X05T.tau
  rw [div_eq_mul_inv, mul_pow]

lemma X05T.pf_2 : X05T.st_2 := by
  intro z3 z5 F hF
  have h1 : (fun d : ℝ => F (1 / (1 - X05T.tau d)) - X05T.pw z3 z5 (X05T.tau d))
      =O[nhds 0] (fun d : ℝ => X05T.tau d ^ 6) := hF.comp_tendsto X05T.tauO
  have h4 := (h1.trans X05T.tau6).add (X05T.h1 z3 z5)
  refine h4.congr' ?_ Filter.EventuallyEq.rfl
  filter_upwards [X05T.ev_ne] with d hd
  rw [X05T.ev_eq F d hd]
  ring

lemma X05T.pf_3 : X05T.st_3 := by
  intro z3 z5 F hF hO n hn
  have h6 := X05T.pf_2 z3 z5 F hO
  have hev : ∀ᶠ x in nhds (1:ℝ), ContDiffAt ℝ 5 F x := hF.eventually (by simp)
  obtain ⟨r, hr, hball⟩ := Metric.eventually_nhds_iff_ball.mp hev
  have hcd : ContDiffOn ℝ 5 F (Metric.ball 1 r) := fun x hx => (hball x hx).contDiffWithinAt
  have hT := taylor_isLittleO (n := 5) (convex_ball (1:ℝ) r) (Metric.mem_ball_self hr)
    (by exact_mod_cast hcd)
  rw [nhdsWithin_eq_nhds.mpr (Metric.isOpen_ball.mem_nhds (Metric.mem_ball_self hr))] at hT
  have htend : Tendsto (fun d : ℝ => 1 + d) (nhds 0) (nhds 1) := by
    have : ContinuousAt (fun d : ℝ => 1 + d) 0 := by fun_prop
    simpa using this.tendsto
  have hT2 : (fun d : ℝ => F (1 + d) -
      ∑ k ∈ Finset.range 6, (iteratedDeriv k F 1 / (k.factorial : ℝ)) * d ^ k)
      =o[nhds 0] (fun d : ℝ => d ^ 5) := by
    have := hT.comp_tendsto htend
    refine this.congr (fun d => ?_) (fun d => ?_)
    · show F (1 + d) - taylorWithinEval F 5 (Metric.ball 1 r) 1 (1 + d) = _
      rw [taylor_within_apply]
      congr 1
      apply Finset.sum_congr rfl
      intro k _
      rw [iteratedDerivWithin_of_isOpen Metric.isOpen_ball (Metric.mem_ball_self hr)]
      simp only [smul_eq_mul, add_sub_cancel_left]
      ring
    · simp
  have h6' : (fun d : ℝ => F (1 + d) - X05T.anPoly z3 z5 d) =o[nhds 0] (fun d : ℝ => d ^ 5) :=
    h6.trans_isLittleO (Asymptotics.isLittleO_pow_pow (m := 5) (n := 6) (by norm_num))
  have hdiff : (fun d : ℝ => ∑ k ∈ Finset.range 6,
      (iteratedDeriv k F 1 / (k.factorial : ℝ) - X05T.a z3 z5 k) * d ^ k)
      =o[nhds 0] (fun d : ℝ => d ^ 5) := by
    refine (h6'.sub hT2).congr (fun d => ?_) (fun _ => rfl)
    simp only [X05T.anPoly, sub_mul, Finset.sum_sub_distrib]
    ring
  have := X05T.poly_zero 5 _ (hdiff.mono nhdsWithin_le_nhds) n hn
  exact sub_eq_zero.mp this

end X05THelpers

namespace X05

lemma tau_bigO : (fun t : ℝ => (fun p => P p / Real.sqrt 2) (1 / (1 - t)) -
    X05T.pw (z 3) (z 5) t) =O[𝓝 0] (fun t : ℝ => t ^ 6) := by
  obtain ⟨C, δ, hδ, hC⟩ := proof_7
  refine Asymptotics.IsBigO.of_bound C ?_
  have hev : ∀ᶠ t in 𝓝 (0 : ℝ), |t| < min δ (1 / 2) :=
    (continuous_abs.tendsto 0).eventually (gt_mem_nhds (by simp [hδ]))
  filter_upwards [hev] with t ht
  have ht1 : t < 1 := by
    have := (abs_lt.mp (ht.trans_le (min_le_right _ _))).2; linarith
  have hp : 0 < 1 / (1 - t) := by apply div_pos one_pos; linarith
  have e : 1 - 1 / (1 / (1 - t)) = t := by rw [one_div_one_div]; ring
  have h := hC (1 / (1 - t)) hp (by rw [e]; exact (ht.trans_le (min_le_left _ _)).le)
  rw [e] at h
  simp only [Real.norm_eq_abs]
  have e6 : |t| ^ 6 = t ^ 6 := by rw [show |t| ^ 6 = (|t| ^ 2) ^ 3 by ring, sq_abs]; ring
  rw [e6] at h
  rw [abs_of_nonneg (by positivity : (0:ℝ) ≤ t ^ 6)]
  unfold X05T.pw
  rw [X05T.ζ2_eq, X05T.ζ4_eq]
  exact h

lemma a_eq (n : ℕ) (hn : n ≤ 5) : a n = X05T.a (z 3) (z 5) n := by
  have hF : ContDiffAt ℝ 5 (fun p => P p / Real.sqrt 2) 1 := (proof_8 5).div_const _
  have := X05T.pf_3 (z 3) (z 5) (fun p => P p / Real.sqrt 2) hF tau_bigO n hn
  rw [iteratedDeriv_div_const] at this
  rw [a, ← this]
  field_simp

end X05

theorem X05.proof_9 : X05.stmt_9 := by
  unfold X05.stmt_9
  refine ⟨?_, ?_, ?_, ?_⟩
  · rw [X05.a_eq 0 (by norm_num)]; rfl
  · rw [X05.a_eq 1 (by norm_num)]; rfl
  · rw [X05.a_eq 2 (by norm_num)]; rfl
  · rw [X05.a_eq 3 (by norm_num)]; rfl

theorem X05.proof_10 : X05.stmt_10 := by
  unfold X05.stmt_10
  rw [X05.a_eq 4 (by norm_num)]; rfl

theorem X05.proof_11 : X05.stmt_11 := by
  unfold X05.stmt_11
  rw [X05.a_eq 5 (by norm_num)]; rfl
