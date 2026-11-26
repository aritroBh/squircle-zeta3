import Mathlib
/-!
# T02 — consistency of Proposition `prop:an` with Theorem `thm:pw`

Target id: T02.
Paper labels: `prop:an` (main.tex lines 1163–1174) and `thm:pw` (main.tex lines 1183–1191),
with the change of variable stated at main.tex line 1181 and the remark at line 1197.

Claims quoted verbatim.

`prop:an` (lines 1164–1168):

  Let $a_n=P^{(n)}(1)/(\sqrt2\,n!)$, the Taylor coefficients of $P(1+d)/\sqrt2$ at $d=0$
  (the series diverges, by Theorem~\ref{thm:zero}). Then $a_0=4$, $a_1=0$, $a_2=\pi^2/6$,
  $a_3=\zeta(3)-\pi^2/3$,
  \begin{gather*}
  a_4=-3\zeta(3)-\frac{11\pi^4}{480}+\frac{\pi^2}{2},\\
  a_5=-\frac{7\pi^2\zeta(3)}{24}+6\zeta(3)-\frac{3\zeta(5)}{4}-\frac{2\pi^2}{3}+\frac{11\pi^4}{120},
  \end{gather*}

Line 1181: Let's use $\tau=1-\frac1p=\frac{d}{1+d}=-\eta$

`thm:pw` (lines 1186–1189):

  \frac{P}{\sqrt2}=4+\zeta(2)\tau^2+\zeta(3)\tau^3-\frac{33}{16}\zeta(4)\tau^4\\
  -\Bigl(\frac34\zeta(5)+\frac74\zeta(2)\zeta(3)\Bigr)\tau^5+O(\tau^6).

Line 1197: going back to $d$ means $\tau=d-d^2+d^3-\cdots$, so
$\zeta(2)\tau^2=\zeta(2)(d^2-2d^3+\cdots)$, and $-2\zeta(2)d^3=-\frac{\pi^2}{3}d^3$ is the
whole difference between $a_3$ and $\zeta(3)$.

Formalization choices:
* `ζ(2)` and `ζ(4)` are the genuine Mathlib values `(riemannZeta 2).re`, `(riemannZeta 4).re`
  (Mathlib's `riemannZeta` is complex-valued; these values are real). The `a_n` are written,
  as in the paper, with powers of `π`; matching them against `ζ(2)`, `ζ(4)` therefore uses
  `riemannZeta_two`, `riemannZeta_four` and is part of the content.
* `ζ(3)` and `ζ(5)` are universally quantified real parameters `z3`, `z5`: the claim is an
  identity of formal coefficients, valid for every value of these symbols (in particular for
  `(riemannZeta 3).re`, `(riemannZeta 5).re`).
* "`O(τ^6)`" and "Taylor coefficients through order 5" are expressed as `IsBigO` at the
  two-sided neighbourhood filter `𝓝 0`. The paper uses both sides of `p = 1` (Fig. `fig:rhombus`
  plots both sides), so the two-sided filter is the faithful reading.
* `F` plays the role of `p ↦ P(p)/√2`; it is an arbitrary real function, so no property of the
  actual perimeter is assumed beyond the displayed hypothesis. `p = 1/(1-τ)` is the inverse of
  `τ = 1 - 1/p`; only values near `τ = 0` (resp. `d = 0`) matter, so the conventions
  `x/0 = 0` play no role.
* stmt_3 adds the hypothesis that `F` is `C^5` near `p = 1`, which is needed for
  `P^{(n)}(1)` to be meaningful (the paper establishes smoothness elsewhere); without it
  `iteratedDeriv` would be the junk value `0`.

NOT COVERED:
* The statements of `thm:pw` and `prop:an` themselves about the actual perimeter `P(p)` of
  `|x|^p+|y|^p=1` (not defined here); this target only checks that the two expansions agree.
* `a_6`, `a_7` (Table `tab:closed`), the weight/polynomiality claims, the `ω_n` part of
  `prop:an` (that is target T03), the rhombus and higher-dimensional extensions of `thm:pw`,
  and the `log(P/(4√2)) = Σ w_k τ^k` form.
-/
-- BEGIN STATEMENTS
namespace T02

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
def stmt_1 : Prop :=
  ∀ z3 z5 : ℝ,
    (fun d : ℝ => pw z3 z5 (tau d) - anPoly z3 z5 d) =O[nhds 0] (fun d : ℝ => d ^ 6)

/-- Consistency of the two results: if `P(p)/√2 = (thm:pw polynomial in τ) + O(τ^6)` with
`p = 1/(1-τ)`, then `P(1+d)/√2 = Σ_{n≤5} a_n d^n + O(d^6)` with the `a_n` of `prop:an`. -/
def stmt_2 : Prop :=
  ∀ (z3 z5 : ℝ) (F : ℝ → ℝ),
    (fun t : ℝ => F (1 / (1 - t)) - pw z3 z5 t) =O[nhds 0] (fun t : ℝ => t ^ 6) →
    (fun d : ℝ => F (1 + d) - anPoly z3 z5 d) =O[nhds 0] (fun d : ℝ => d ^ 6)

/-- `prop:an` as literally defined (`a_n = F^{(n)}(1)/n!` with `F = P/√2`) follows from
`thm:pw`: if `F` is `C^5` at `1` and satisfies the `thm:pw` expansion, then
`F^{(n)}(1)/n! = a_n` for `n = 0, …, 5`. -/
def stmt_3 : Prop :=
  ∀ (z3 z5 : ℝ) (F : ℝ → ℝ), ContDiffAt ℝ 5 F 1 →
    (fun t : ℝ => F (1 / (1 - t)) - pw z3 z5 t) =O[nhds 0] (fun t : ℝ => t ^ 6) →
    ∀ n ≤ 5, iteratedDeriv n F 1 / (n.factorial : ℝ) = a z3 z5 n

end T02
-- END STATEMENTS

section Helpers
open Asymptotics Filter

noncomputable def T02.Q (z3 z5 d : ℝ) : ℝ :=
  ((5/6)*Real.pi^2 + (-10)*z3 + (-11/48)*Real.pi^4 + (35/24)*z3*Real.pi^2 + (15/4)*z5)
  + ((19/6)*Real.pi^2 + (-35)*z3 + (-11/16)*Real.pi^4 + (35/12)*z3*Real.pi^2 + (15/2)*z5) * d
  + ((-46)*z3 + (9/2)*Real.pi^2 + (-77/96)*Real.pi^4 + (35/12)*z3*Real.pi^2 + (15/2)*z5) * d^2
  + ((-27)*z3 + (-209/480)*Real.pi^4 + (17/6)*Real.pi^2 + (35/24)*z3*Real.pi^2 + (15/4)*z5) * d^3
  + ((7/24)*z3*Real.pi^2 + (-6)*z3 + (3/4)*z5 + (2/3)*Real.pi^2 + (-11/120)*Real.pi^4) * d^4

lemma T02.ζ2_eq : T02.ζ2 = Real.pi ^ 2 / 6 := by
  unfold T02.ζ2
  rw [riemannZeta_two]
  have : ((Real.pi : ℂ) ^ 2 / 6) = ((Real.pi ^ 2 / 6 : ℝ) : ℂ) := by push_cast; ring
  rw [this, Complex.ofReal_re]

lemma T02.ζ4_eq : T02.ζ4 = Real.pi ^ 4 / 90 := by
  unfold T02.ζ4
  rw [riemannZeta_four]
  have : ((Real.pi : ℂ) ^ 4 / 90) = ((Real.pi ^ 4 / 90 : ℝ) : ℂ) := by push_cast; ring
  rw [this, Complex.ofReal_re]

lemma T02.key (z3 z5 d : ℝ) (hd : 1 + d ≠ 0) :
    T02.pw z3 z5 (T02.tau d) - T02.anPoly z3 z5 d = d ^ 6 * (T02.Q z3 z5 d / (1 + d) ^ 5) := by
  unfold T02.pw T02.tau T02.anPoly T02.Q
  rw [T02.ζ2_eq, T02.ζ4_eq]
  simp only [Finset.sum_range_succ, Finset.sum_range_zero, T02.a]
  field_simp
  ring

lemma T02.tauO : Filter.Tendsto T02.tau (nhds 0) (nhds 0) := by
  have : ContinuousAt T02.tau 0 := by
    unfold T02.tau
    have hc : ContinuousAt (fun d : ℝ => 1 + d) 0 := by fun_prop
    exact continuousAt_id.div hc (by norm_num : (1:ℝ) + 0 ≠ 0)
  simpa [T02.tau] using this.tendsto

lemma T02.ev_ne : ∀ᶠ d in nhds (0:ℝ), 1 + d ≠ 0 :=
  have hc : ContinuousAt (fun d : ℝ => 1 + d) 0 := by fun_prop
  hc.eventually_ne (by norm_num : (1:ℝ) + 0 ≠ 0)

lemma T02.h1 (z3 z5 : ℝ) :
    (fun d : ℝ => T02.pw z3 z5 (T02.tau d) - T02.anPoly z3 z5 d) =O[nhds 0] (fun d : ℝ => d ^ 6) := by
  have hg : ContinuousAt (fun d : ℝ => T02.Q z3 z5 d / (1 + d) ^ 5) 0 := by
    have hQ : ContinuousAt (fun d : ℝ => T02.Q z3 z5 d) 0 :=
      (by unfold T02.Q; fun_prop : Continuous (fun d : ℝ => T02.Q z3 z5 d)).continuousAt
    have hden : ContinuousAt (fun d : ℝ => (1 + d) ^ 5) 0 := by fun_prop
    exact hQ.div hden (by simp)
  have hb := hg.isBigO (F := ℝ)
  have h2 : (fun d : ℝ => d ^ 6 * (T02.Q z3 z5 d / (1 + d) ^ 5)) =O[nhds 0] (fun d : ℝ => d ^ 6) := by
    simpa using (Asymptotics.isBigO_refl (fun d : ℝ => d ^ 6) (nhds 0)).mul hb
  refine h2.congr' ?_ Filter.EventuallyEq.rfl
  filter_upwards [T02.ev_ne] with d hd
  exact (T02.key z3 z5 d hd).symm

lemma T02.poly_zero : ∀ (m : ℕ) (b : ℕ → ℝ),
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

theorem T02.proof_1 : T02.stmt_1 := by
  intro z3 z5
  exact T02.h1 z3 z5

lemma T02.ev_eq (F : ℝ → ℝ) (d : ℝ) (hd : 1 + d ≠ 0) : F (1 / (1 - T02.tau d)) = F (1 + d) := by
  congr 1
  unfold T02.tau
  field_simp
  ring

lemma T02.tau6 : (fun d : ℝ => T02.tau d ^ 6) =O[nhds 0] (fun d : ℝ => d ^ 6) := by
  have hg : ContinuousAt (fun d : ℝ => ((1 + d)⁻¹) ^ 6) 0 := by
    have hc : ContinuousAt (fun d : ℝ => 1 + d) 0 := by fun_prop
    have hi : ContinuousAt (fun d : ℝ => (1 + d)⁻¹) 0 := hc.inv₀ (by simp)
    exact hi.pow 6
  have hb := hg.isBigO (F := ℝ)
  have h2 : (fun d : ℝ => d ^ 6 * ((1 + d)⁻¹) ^ 6) =O[nhds 0] (fun d : ℝ => d ^ 6) := by
    simpa using (Asymptotics.isBigO_refl (fun d : ℝ => d ^ 6) (nhds 0)).mul hb
  refine h2.congr' ?_ Filter.EventuallyEq.rfl
  filter_upwards with d
  unfold T02.tau
  rw [div_eq_mul_inv, mul_pow]

theorem T02.proof_2 : T02.stmt_2 := by
  intro z3 z5 F hF
  have h1 : (fun d : ℝ => F (1 / (1 - T02.tau d)) - T02.pw z3 z5 (T02.tau d))
      =O[nhds 0] (fun d : ℝ => T02.tau d ^ 6) := hF.comp_tendsto T02.tauO
  have h4 := (h1.trans T02.tau6).add (T02.h1 z3 z5)
  refine h4.congr' ?_ Filter.EventuallyEq.rfl
  filter_upwards [T02.ev_ne] with d hd
  rw [T02.ev_eq F d hd]
  ring

theorem T02.proof_3 : T02.stmt_3 := by
  intro z3 z5 F hF hO n hn
  have h6 := T02.proof_2 z3 z5 F hO
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
  have h6' : (fun d : ℝ => F (1 + d) - T02.anPoly z3 z5 d) =o[nhds 0] (fun d : ℝ => d ^ 5) :=
    h6.trans_isLittleO (Asymptotics.isLittleO_pow_pow (m := 5) (n := 6) (by norm_num))
  have hdiff : (fun d : ℝ => ∑ k ∈ Finset.range 6,
      (iteratedDeriv k F 1 / (k.factorial : ℝ) - T02.a z3 z5 k) * d ^ k)
      =o[nhds 0] (fun d : ℝ => d ^ 5) := by
    refine (h6'.sub hT2).congr (fun d => ?_) (fun _ => rfl)
    simp only [T02.anPoly, sub_mul, Finset.sum_sub_distrib]
    ring
  have := T02.poly_zero 5 _ (hdiff.mono nhdsWithin_le_nhds) n hn
  exact sub_eq_zero.mp this

end Helpers
