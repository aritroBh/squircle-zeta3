import Mathlib
import Zeta3.Wave3.W02_Holo
import Zeta3.Wave8.Q05_Decay
/-!
# Q07 -- the rest of the decay law of `M_k` (continuation of Q05, `thm:decay`)

Target id: Q07. Continuation of `Q05` (paper/ieee/main.tex `thm:decay`, supplement
paper/ieee/proofs/sec_g5_decay.tex and sec_g5_varpi1.tex). The definitions (`H`, `K`, `Lpt`,
`lineInt`, `P`, `W`, `D`, `M`, `Mmain`, `B1`, `B2`, `ξ`, `Cstar`, `X`, ...) are those of `Q05`;
the statements of Q05 that are proved stay in `Q05_Decay.lean`. The statements below were
audited as Q05.stmt_9, 10, 12, 14, 15, 16, 18 and are stated here, unproved (planner split,
2026-11-26); the original docstrings are kept.

## Statements

* stmt_1 (was Q05.stmt_9): g5:lem:kernel.
* stmt_2 (was Q05.stmt_10): g5:lem:Mk.
* stmt_3 (was Q05.stmt_12): g5:lem:W.
* stmt_4 (was Q05.stmt_14): g5:prop:rest, the displayed bound.
* stmt_5 (was Q05.stmt_15): g5:prop:rest (a) as used by the paper: `|M_k - M_k^main| ≤ 10,600,000 e^{-2√N}`.
* stmt_6 (was Q05.stmt_16): g5:prop:main: `|M_k^main + X| ≤ 71,000 N^{-2} e^{-θ}`.
* stmt_7 (was Q05.stmt_18): thm:decay in full: the inequality with `1.21·10^9`, and eq:Mk in the
  explicit form of the supplement (`k ≥ 99`, constant `1.21·10^9 + 1.9√π`).

## NOT COVERED

* The paper's interval-arithmetic bound with the constant `254.23` (main.tex lines 1281-1284,
  85,210 cells; audited as old Q05.stmt_17) has no written proof in the supplement and is not
  stated.
* g5:prop:rest (b) (`c_rest < 120,000`) and the internals of g5:prop:main (Regions I and II of
  sec_g5_varpi1.tex) are not stated.
-/
-- BEGIN STATEMENTS
namespace Q07

open Real MeasureTheory Set Filter Q05

/-- g5:lem:kernel: "For every real $z\ge0$,
$C(z)=1-\frac2\pi\operatorname{Im}\int_{\Lline}H(s)K(s)\,\eu^{-sz}\,ds$, and the integral
converges absolutely." (was Q05.stmt_9) -/
def stmt_1 : Prop :=
  ∀ z : ℝ, 0 ≤ z →
    IntegrableOn (fun r : ℝ => H (Lpt r) * K (Lpt r) * Complex.exp (-(Lpt r) * z)) (Ioi 0) ∧
      Cf z = 1 - 2 / π * (lineInt fun s => H s * K s * Complex.exp (-s * z)).im

/-- g5:lem:Mk: "For every $k\ge0$, $M_k=-D_k+\frac2\pi\operatorname{Im}\int_{\Lline}H(s)K(s)P_k(s)\,ds$,
and the integral converges absolutely." Stated for the complex coefficient `W02.Mk k`, so it
also says that `M_k` is real. (was Q05.stmt_10) -/
def stmt_2 : Prop :=
  ∀ k : ℕ,
    IntegrableOn (fun r : ℝ => H (Lpt r) * K (Lpt r) * P k (Lpt r)) (Ioi 0) ∧
      W02.Mk k = ((-D k + 2 / π * (lineInt fun s => H s * K s * P k s).im : ℝ) : ℂ)

/-- g5:lem:W: for `k ≥ 1`, `Re s > 0`, `0 ≤ arg s ≤ π/4`, `ρ = |s|`,
`|W_k(s)| ≤ 2ρ^k/|s+2|^{k+1} + 2 Σ_{i=1}^k (2^i/i!) ρ^{k-i}/|s+i+1|^{k-i+1}`. (was Q05.stmt_12) -/
def stmt_3 : Prop :=
  ∀ k : ℕ, 1 ≤ k → ∀ s : ℂ, 0 < s.re → 0 ≤ s.arg → s.arg ≤ π / 4 →
    ‖W k s‖ ≤ 2 * ‖s‖ ^ k / ‖s + 2‖ ^ (k + 1) +
      2 * ∑ i ∈ Finset.Icc 1 k,
        (2 : ℝ) ^ i / (Nat.factorial i : ℝ) * (‖s‖ ^ (k - i) / ‖s + (i + 1 : ℂ)‖ ^ (k - i + 1))

/-- g5:prop:rest, the displayed inequality, for every `k` with `N = k+1 ≥ 100`. (was Q05.stmt_14) -/
def stmt_4 : Prop :=
  ∀ k : ℕ, 100 ≤ k + 1 →
    |M k - Mmain k| ≤
      B1 * (Real.exp (π / 4) * Real.exp (-√(3 * π * k / 2)) +
          ξ * Real.exp (-√(3 * π * ((k : ℝ) - 1) / 2)) +
          (Real.exp ξ - 1 - ξ) * Real.exp (-√(9 * π * ((k : ℝ) - 2) / 4)) +
          Real.exp (Real.exp 1 * ξ) * Real.exp (-(k : ℝ) / 2)) +
        B2 * Real.exp (-√(9 * π * k / 4)) + 2 * Real.exp (-(((k + 1 : ℕ) : ℝ) - 1) / 2)

/-- g5:prop:rest (a) with the constant of the paper: "$|M_k-M_k^{\mathrm{main}}|\le10{,}600{,}000\,
e^{-2\sqrt N}$ for $N\ge100$" (main.tex proof idea), together with the supplement's
`|M_k - M_k^main| ≤ C_⋆ e^{-2√N}`. (was Q05.stmt_15) -/
def stmt_5 : Prop :=
  ∀ k : ℕ, 100 ≤ k + 1 →
    |M k - Mmain k| ≤ Cstar * Real.exp (-2 * √((k + 1 : ℕ) : ℝ)) ∧
      |M k - Mmain k| ≤ 10600000 * Real.exp (-2 * √((k + 1 : ℕ) : ℝ))

/-- g5:prop:main (sec_g5_varpi1.tex): "For every integer $N\ge100$,
$|M_k^{\mathrm{main}}+X|\le71{,}000\,N^{-2}\eu^{-\theta}$", `θ = √(πN)`. (was Q05.stmt_16) -/
def stmt_6 : Prop :=
  ∀ k : ℕ, 100 ≤ k + 1 →
    |Mmain k + X ((k + 1 : ℕ) : ℝ)| ≤
      71000 / (((k + 1 : ℕ) : ℝ)) ^ 2 * Real.exp (-√(π * ((k + 1 : ℕ) : ℝ)))

/-- thm:decay in full: "For every $k$ with $N=k+1\ge100$,
$|M_k+\frac{\sqrt\pi}{N}e^{-\sqrt{\pi N}}\cos\sqrt{\pi N}|\le1.21\cdot10^{9}N^{-3/2}
e^{-\sqrt{\pi N}}$. In particular \eqref{eq:Mk} holds", with eq:Mk in the explicit form of the
supplement: for `k ≥ 99`, `|M_k + (√π/k) e^{-√(πk)} cos √(πk)| ≤ (1.21·10^9 + 1.9√π)
k^{-3/2} e^{-√(πk)}`, i.e. `M_k = -(√π/k)e^{-√(πk)}[cos√(πk) + O(k^{-1/2})]` with an explicit
constant. (was Q05.stmt_18) -/
def stmt_7 : Prop :=
  (∀ k : ℕ, 100 ≤ k + 1 →
      |M k + √π / ((k + 1 : ℕ) : ℝ) * Real.exp (-√(π * ((k + 1 : ℕ) : ℝ))) *
          Real.cos (√(π * ((k + 1 : ℕ) : ℝ)))| ≤
        1.21e9 * (((k + 1 : ℕ) : ℝ)) ^ (-(3 / 2) : ℝ) * Real.exp (-√(π * ((k + 1 : ℕ) : ℝ)))) ∧
    ∀ k : ℕ, 99 ≤ k →
      |M k + √π / (k : ℝ) * Real.exp (-√(π * k)) * Real.cos (√(π * k))| ≤
        (1.21e9 + √π * 1.9) * (k : ℝ) ^ (-(3 / 2) : ℝ) * Real.exp (-√(π * k))

end Q07
-- END STATEMENTS

/-! ### Merged from work/Q07/PartA.lean -/
namespace Q07.PA
open Real MeasureTheory Set Filter Q05 Topology

/-- Removing finitely many simple poles with known residues. -/
lemma glue (S : Finset ℂ) : ∀ (U : Set ℂ), IsOpen U → (↑S ⊆ U) → ∀ (f c : ℂ → ℂ),
    DifferentiableOn ℂ f (U \ ↑S) →
    (∀ p ∈ S, Tendsto (fun s => (s - p) * f s) (𝓝[≠] p) (𝓝 (c p))) →
    ∃ g : ℂ → ℂ, DifferentiableOn ℂ g U ∧ ∀ s ∈ U \ ↑S, g s = f s - ∑ p ∈ S, c p / (s - p) := by
  classical
  induction S using Finset.induction_on with
  | empty =>
    intro U _ _ f c hf _
    exact ⟨f, by simpa using hf, fun s _ => by simp⟩
  | insert p S' hp ih =>
    intro U hU hS f c hf hres
    have hpU : p ∈ U := hS (Finset.mem_insert_self p S')
    have hU' : IsOpen (U \ {p}) := hU.sdiff isClosed_singleton
    have hS' : (↑S' : Set ℂ) ⊆ U \ {p} := by
      intro q hq
      refine ⟨hS (Finset.mem_insert_of_mem hq), ?_⟩
      rintro rfl; exact hp hq
    have hf' : DifferentiableOn ℂ f ((U \ {p}) \ (↑S' : Set ℂ)) := by
      refine hf.mono ?_
      intro s hs; simp only [Set.mem_sdiff, mem_singleton_iff, Finset.coe_insert, mem_insert_iff,
        Finset.mem_coe] at hs ⊢
      tauto
    obtain ⟨g', hg'd, hg'eq⟩ := ih (U \ {p}) hU' hS' f c hf'
      (fun q hq => hres q (Finset.mem_insert_of_mem hq))
    set h : ℂ → ℂ := fun s => g' s - c p / (s - p) with hh
    have hhd : DifferentiableOn ℂ h (U \ {p}) := by
      refine hg'd.sub ?_
      intro s hs
      have : s - p ≠ 0 := sub_ne_zero.2 hs.2
      exact ((differentiableAt_const _).div (differentiableAt_id.sub_const p) this).differentiableWithinAt
    -- (s - p) h s → 0
    have hev : ∀ᶠ s in 𝓝[≠] p, (s - p) * h s =
        (s - p) * f s - ∑ q ∈ S', c q * ((s - p) / (s - q)) - c p := by
      have h1 : ∀ᶠ s in 𝓝 p, s ∈ U := hU.mem_nhds hpU
      have h2 : ∀ᶠ s in 𝓝 p, s ∉ (↑S' : Set ℂ) :=
        (S'.finite_toSet.isClosed.isOpen_compl).mem_nhds (by simpa using hp)
      have h1' : ∀ᶠ s in 𝓝[≠] p, s ∈ U := h1.filter_mono nhdsWithin_le_nhds
      have h2' : ∀ᶠ s in 𝓝[≠] p, s ∉ (↑S' : Set ℂ) := h2.filter_mono nhdsWithin_le_nhds
      filter_upwards [h1', h2', self_mem_nhdsWithin] with s hs1 hs2 hs3
      have hsp : s - p ≠ 0 := sub_ne_zero.2 hs3
      rw [hh]; simp only
      rw [hg'eq s ⟨⟨hs1, hs3⟩, hs2⟩, mul_sub, mul_sub, Finset.mul_sum]
      have e1 : (s - p) * (c p / (s - p)) = c p := by field_simp
      rw [e1]
      congr 2
      exact Finset.sum_congr rfl (fun q _ => by ring)
    have hlim0 : Tendsto (fun s => (s - p) * h s) (𝓝[≠] p) (𝓝 0) := by
      have hsum : Tendsto (fun s => ∑ q ∈ S', c q * ((s - p) / (s - q))) (𝓝[≠] p) (𝓝 0) := by
        have : (0 : ℂ) = ∑ q ∈ S', c q * ((p - p) / (p - q)) := by simp
        rw [this]
        refine tendsto_nhdsWithin_of_tendsto_nhds ?_
        refine tendsto_finsetSum _ (fun q hq => ?_)
        have hqp : p - q ≠ 0 := sub_ne_zero.2 (fun h => hp (h ▸ hq))
        exact ((continuousAt_const.mul ((continuousAt_id.sub continuousAt_const).div
          (continuousAt_id.sub continuousAt_const) hqp))).tendsto
      have := ((hres p (Finset.mem_insert_self p S')).sub hsum).sub_const (c p)
      rw [sub_zero, sub_self] at this
      exact this.congr' (hev.mono fun s hs => hs.symm)
    have ho : (fun z => h z - h p) =o[𝓝[≠] p] fun z => (z - p)⁻¹ := by
      rw [Asymptotics.isLittleO_iff_tendsto']
      · have : Tendsto (fun z => (z - p) * h z - (z - p) * h p) (𝓝[≠] p) (𝓝 0) := by
          have h2 : Tendsto (fun z => (z - p) * h p) (𝓝[≠] p) (𝓝 0) := by
            have : Tendsto (fun z : ℂ => (z - p) * h p) (𝓝 p) (𝓝 ((p - p) * h p)) :=
              ((continuousAt_id.sub continuousAt_const).mul continuousAt_const).tendsto
            rw [sub_self, zero_mul] at this
            exact tendsto_nhdsWithin_of_tendsto_nhds this
          simpa using hlim0.sub h2
        refine this.congr' ?_
        filter_upwards [self_mem_nhdsWithin] with z hz
        rw [div_inv_eq_mul]; ring
      · filter_upwards [self_mem_nhdsWithin] with z hz
        intro h0
        exact absurd h0 (inv_ne_zero (sub_ne_zero.2 hz))
    have hfin := Complex.differentiableOn_update_limUnder_of_isLittleO (hU.mem_nhds hpU) hhd ho
    refine ⟨_, hfin, fun s hs => ?_⟩
    have hsp : s ≠ p := fun h => hs.2 (by simp [h])
    rw [Function.update_of_ne hsp]
    have hs' : s ∈ (U \ {p}) \ (↑S' : Set ℂ) := by
      refine ⟨⟨hs.1, hsp⟩, fun h => hs.2 (by simp [Finset.mem_coe.1 h])⟩
    simp only [hh]
    rw [hg'eq s hs', Finset.sum_insert hp]
    ring


lemma alg (A w s c : ℝ) (hw : 0 < w) (hc : 0 < c) (hc2 : c ^ 2 = 1 - s ^ 2) (hA2 : A ^ 2 = 1 + w ^ 2)
    (hD : 0 < 1 + (A - 1) * s ^ 2) :
    -(((A - 1) * (1 / c) + A * (1 / (1 + ((A - 1) * s * c / (1 + (A - 1) * s ^ 2)) ^ 2) *
      ((((A - 1) * 1 * c + (A - 1) * s * (-(2 * s) / (2 * c))) * (1 + (A - 1) * s ^ 2) -
        (A - 1) * s * c * ((A - 1) * (2 * s))) / (1 + (A - 1) * s ^ 2) ^ 2))) / w) =
      -(c * w) / (1 + w ^ 2 * s ^ 2) := by
  have e1 : (1 + (A - 1) * s ^ 2) ^ 2 + ((A - 1) * s * c) ^ 2 = 1 + w ^ 2 * s ^ 2 := by
    linear_combination (A - 1) ^ 2 * s ^ 2 * hc2 + s ^ 2 * hA2
  have hE : 0 < 1 + w ^ 2 * s ^ 2 := by positivity
  have e2 : 1 / (1 + ((A - 1) * s * c / (1 + (A - 1) * s ^ 2)) ^ 2) =
      (1 + (A - 1) * s ^ 2) ^ 2 / (1 + w ^ 2 * s ^ 2) := by
    rw [← e1]; field_simp
  rw [e2]
  field_simp
  linear_combination (A ^ 3 * s ^ 2 - 2 * A ^ 2 * s ^ 2 - A ^ 2 + A * s ^ 2 + A + w ^ 2) * hc2 + (A * s ^ 2 - 1) * hA2

/-- The primitive used for the real integral, as a function of `s = e^{-u}`. -/
noncomputable def Fc (w s : ℝ) : ℝ :=
  -(((√(1 + w ^ 2) - 1) * Real.arcsin s + √(1 + w ^ 2) *
    Real.arctan ((√(1 + w ^ 2) - 1) * s * √(1 - s ^ 2) / (1 + (√(1 + w ^ 2) - 1) * s ^ 2))) / w)

lemma Fc_deriv (w s : ℝ) (hw : 0 < w) (hs0 : 0 < s) (hs1 : s < 1) :
    HasDerivAt (Fc w) (-(√(1 - s ^ 2) * w) / (1 + w ^ 2 * s ^ 2)) s := by
  set A := √(1 + w ^ 2) with hAdef
  have hA2 : A ^ 2 = 1 + w ^ 2 := Real.sq_sqrt (by positivity)
  have hA1 : 1 ≤ A := by rw [hAdef]; exact Real.one_le_sqrt.2 (by nlinarith)
  set c := √(1 - s ^ 2) with hcdef
  have hc2 : c ^ 2 = 1 - s ^ 2 := Real.sq_sqrt (by nlinarith)
  have hcpos : 0 < c := Real.sqrt_pos.2 (by nlinarith)
  have hD : 0 < 1 + (A - 1) * s ^ 2 := by nlinarith
  have hsq : HasDerivAt (fun y : ℝ => √(1 - y ^ 2)) (-(2 * s) / (2 * c)) s := by
    have h1 : HasDerivAt (fun y : ℝ => 1 - y ^ 2) (-(2 * s)) s := by
      simpa using ((hasDerivAt_id s).pow 2).const_sub 1
    exact h1.sqrt (by nlinarith)
  have hN : HasDerivAt (fun y : ℝ => (A - 1) * y * √(1 - y ^ 2))
      ((A - 1) * 1 * c + (A - 1) * s * (-(2 * s) / (2 * c))) s :=
    (((hasDerivAt_id s).const_mul (A - 1)).mul hsq)
  have hDd : HasDerivAt (fun y : ℝ => 1 + (A - 1) * y ^ 2) ((A - 1) * (2 * s)) s := by
    simpa using (((hasDerivAt_id s).pow 2).const_mul (A - 1)).const_add 1
  have hq := (hN.div hDd hD.ne').arctan
  have has := Real.hasDerivAt_arcsin (x := s) (by linarith) hs1.ne
  have htot := ((has.const_mul (A - 1)).add (hq.const_mul A)).div_const w |>.neg
  convert htot using 1
  · ext y; simp [Fc, hAdef]
  · simp only [Pi.div_apply]
    exact (alg A w s c hw hcpos hc2 hA2 hD).symm

lemma Fc_cont (w : ℝ) : Continuous (Fc w) := by
  unfold Fc
  have hD : ∀ s : ℝ, 1 + (√(1 + w ^ 2) - 1) * s ^ 2 ≠ 0 := by
    intro s
    have : 1 ≤ √(1 + w ^ 2) := Real.one_le_sqrt.2 (by nlinarith)
    have : 0 ≤ (√(1 + w ^ 2) - 1) * s ^ 2 := by positivity
    linarith
  refine Continuous.neg (Continuous.div_const ?_ _)
  refine Continuous.add (continuous_const.mul Real.continuous_arcsin)
    (continuous_const.mul (Real.continuous_arctan.comp ?_))
  have h1 : Continuous (fun s : ℝ => (√(1 + w ^ 2) - 1) * s * √(1 - s ^ 2)) := by fun_prop
  have h2 : Continuous (fun s : ℝ => 1 + (√(1 + w ^ 2) - 1) * s ^ 2) := by fun_prop
  exact h1.div h2 hD

lemma Fc_zero (w : ℝ) : Fc w 0 = 0 := by simp [Fc]

lemma Fc_one (w : ℝ) : Fc w 1 = -((√(1 + w ^ 2) - 1) * (π / 2) / w) := by simp [Fc]

/-- The integrand of the real integral, `x = w e^{-u}`. -/
noncomputable def Gint (w u : ℝ) : ℝ :=
  √(1 - Real.exp (-u) ^ 2) * (w * Real.exp (-u)) / (1 + w ^ 2 * Real.exp (-u) ^ 2)

lemma Gint_cont (w : ℝ) : Continuous (Gint w) := by
  unfold Gint
  have h1 : Continuous (fun u : ℝ => √(1 - Real.exp (-u) ^ 2) * (w * Real.exp (-u))) := by
    fun_prop
  have h2 : Continuous (fun u : ℝ => 1 + w ^ 2 * Real.exp (-u) ^ 2) := by fun_prop
  exact h1.div h2 (fun u => by positivity)

lemma Gint_le (w u : ℝ) (hw : 0 < w) : |Gint w u| ≤ w * Real.exp (-u) := by
  unfold Gint
  have h1 : √(1 - Real.exp (-u) ^ 2) ≤ 1 := by
    rw [Real.sqrt_le_one]; nlinarith [sq_nonneg (Real.exp (-u))]
  have h0 : 0 ≤ √(1 - Real.exp (-u) ^ 2) := Real.sqrt_nonneg _
  have hD : 1 ≤ 1 + w ^ 2 * Real.exp (-u) ^ 2 := by nlinarith [sq_nonneg (w * Real.exp (-u))]
  have hx : 0 < w * Real.exp (-u) := by positivity
  rw [abs_of_nonneg (by positivity), div_le_iff₀ (by linarith)]
  nlinarith

lemma real_int (w : ℝ) (hw : 0 < w) :
    ∫ u in Ioi (0 : ℝ), Gint w u = π / 2 * (w / (1 + √(1 + w ^ 2))) := by
  have hder : ∀ u ∈ Ioi (0 : ℝ), HasDerivAt (fun u => Fc w (Real.exp (-u))) (Gint w u) u := by
    intro u hu
    have hs0 : 0 < Real.exp (-u) := Real.exp_pos _
    have hs1 : Real.exp (-u) < 1 := Real.exp_lt_one_iff.2 (by simp at hu; linarith)
    have he : HasDerivAt (fun u : ℝ => Real.exp (-u)) (-Real.exp (-u)) u := by
      simpa using (hasDerivAt_neg u).exp
    have := (Fc_deriv w _ hw hs0 hs1).comp u he
    exact this.congr_deriv (by unfold Gint; ring)
  have hint : IntegrableOn (Gint w) (Ioi 0) := by
    refine Integrable.mono' ((exp_neg_integrableOn_Ioi 0 one_pos).const_mul w)
      (Gint_cont w).aestronglyMeasurable ?_
    refine Filter.Eventually.of_forall fun u => ?_
    rw [Real.norm_eq_abs]
    simpa using Gint_le w u hw
  have ht : Tendsto (fun u => Fc w (Real.exp (-u))) atTop (𝓝 (Fc w 0)) :=
    ((Fc_cont w).tendsto 0).comp Real.tendsto_exp_neg_atTop_nhds_zero
  rw [integral_Ioi_of_hasDerivAt_of_tendsto
    ((Fc_cont w).comp (by fun_prop : Continuous fun u : ℝ => Real.exp (-u))).continuousWithinAt
    hder hint ht]
  simp only [Function.comp, Fc_zero, neg_zero, Real.exp_zero, Fc_one]
  have hA2 : √(1 + w ^ 2) ^ 2 = 1 + w ^ 2 := Real.sq_sqrt (by positivity)
  have hA : 0 < 1 + √(1 + w ^ 2) := by positivity
  field_simp
  linear_combination (π * 2 - π) * hA2

/-- The `n`-th term `√(1-e^{-2u}) x (-x^2)^n`, `x = w e^{-u}`. -/
noncomputable def tn (w u : ℝ) (n : ℕ) : ℝ :=
  √(1 - Real.exp (-u) ^ 2) * (w * Real.exp (-u) * (-(w * Real.exp (-u)) ^ 2) ^ n)

lemma tn_cont (w : ℝ) (n : ℕ) : Continuous (fun u => tn w u n) := by
  unfold tn; fun_prop

lemma tn_abs (w u : ℝ) (n : ℕ) (hw0 : 0 < w) (hw1 : w ≤ 1) (hu : 0 ≤ u) :
    |tn w u n| ≤ Real.exp (-u) := by
  unfold tn
  have h1 : √(1 - Real.exp (-u) ^ 2) ≤ 1 := by
    rw [Real.sqrt_le_one]; nlinarith [sq_nonneg (Real.exp (-u))]
  have he1 : Real.exp (-u) ≤ 1 := Real.exp_le_one_iff.2 (by linarith)
  have he0 : 0 < Real.exp (-u) := Real.exp_pos _
  have hx1 : w * Real.exp (-u) ≤ 1 := by nlinarith
  have hy : |(-(w * Real.exp (-u)) ^ 2) ^ n| ≤ 1 := by
    rw [abs_pow]; apply pow_le_one₀ (abs_nonneg _)
    rw [abs_neg, abs_of_nonneg (by positivity)]; exact pow_le_one₀ (by positivity) hx1
  rw [abs_mul, abs_mul, abs_of_nonneg (Real.sqrt_nonneg _), abs_of_pos (by positivity)]
  have : |(-(w * Real.exp (-u)) ^ 2) ^ n| * (w * Real.exp (-u)) ≤ 1 * Real.exp (-u) := by
    have := mul_le_mul hy (by nlinarith : w * Real.exp (-u) ≤ Real.exp (-u)) (by positivity)
      zero_le_one
    linarith
  have h0 := Real.sqrt_nonneg (1 - Real.exp (-u) ^ 2)
  calc √(1 - Real.exp (-u) ^ 2) * (w * Real.exp (-u) * |(-(w * Real.exp (-u)) ^ 2) ^ n|)
      ≤ 1 * (w * Real.exp (-u) * |(-(w * Real.exp (-u)) ^ 2) ^ n|) := by
        apply mul_le_mul_of_nonneg_right h1; positivity
    _ ≤ Real.exp (-u) := by linarith

lemma tn_int (w : ℝ) (n : ℕ) (hw0 : 0 < w) (hw1 : w ≤ 1) :
    IntegrableOn (fun u => tn w u n) (Ioi 0) := by
  refine Integrable.mono' (exp_neg_integrableOn_Ioi 0 one_pos)
    (tn_cont w n).aestronglyMeasurable ?_
  refine (ae_restrict_iff' measurableSet_Ioi).2 (Filter.Eventually.of_forall fun u hu => ?_)
  rw [Real.norm_eq_abs]
  simpa using tn_abs w u n hw0 hw1 (le_of_lt hu)

lemma term_eq (z : ℝ) (n : ℕ) :
    (-1 : ℂ) ^ n * (H (2 * (n : ℂ) + 1) * Complex.exp (-(2 * (n : ℂ) + 1) * z)) =
      ((∫ u in Ioi (0 : ℝ), tn (Real.exp (-z)) u n : ℝ) : ℂ) := by
  have hre : 0 < (2 * (n : ℂ) + 1).re := by simp; positivity
  obtain ⟨-, hH⟩ := H_rep _ hre
  rw [hH, ← integral_complex_ofReal, ← integral_mul_const, ← integral_const_mul]
  refine integral_congr_ae (Filter.Eventually.of_forall fun u => ?_)
  simp only [Hint, tn]
  have e1 : Complex.exp (-(2 * (n : ℂ) + 1) * u) = ((Real.exp (-u) ^ (2 * n + 1) : ℝ) : ℂ) := by
    rw [← Real.exp_nat_mul, Complex.ofReal_exp]; push_cast; ring_nf
  have e2 : Complex.exp (-(2 * (n : ℂ) + 1) * z) = ((Real.exp (-z) ^ (2 * n + 1) : ℝ) : ℂ) := by
    rw [← Real.exp_nat_mul, Complex.ofReal_exp]; push_cast; ring_nf
  have e3 : Real.exp (-2 * u) = Real.exp (-u) ^ 2 := by
    rw [← Real.exp_nat_mul]; ring_nf
  rw [e1, e2, e3]
  push_cast
  ring

lemma sum_eq (z : ℝ) (hz : 0 ≤ z) (J : ℕ) :
    ∑ n ∈ Finset.range J,
        (-1 : ℂ) ^ n * (H (2 * (n : ℂ) + 1) * Complex.exp (-(2 * (n : ℂ) + 1) * z)) =
      ((∫ u in Ioi (0 : ℝ), ∑ n ∈ Finset.range J, tn (Real.exp (-z)) u n : ℝ) : ℂ) := by
  have hw0 : 0 < Real.exp (-z) := Real.exp_pos _
  have hw1 : Real.exp (-z) ≤ 1 := Real.exp_le_one_iff.2 (by linarith)
  rw [integral_finsetSum _ (fun n _ => tn_int _ n hw0 hw1)]
  push_cast
  exact Finset.sum_congr rfl (fun n _ => term_eq z n)

lemma sum_tendsto (z : ℝ) (hz : 0 ≤ z) :
    Tendsto (fun J : ℕ => ∫ u in Ioi (0 : ℝ), ∑ n ∈ Finset.range J, tn (Real.exp (-z)) u n)
      atTop (𝓝 (∫ u in Ioi (0 : ℝ), Gint (Real.exp (-z)) u)) := by
  set w := Real.exp (-z) with hwdef
  have hw0 : 0 < w := Real.exp_pos _
  have hw1 : w ≤ 1 := Real.exp_le_one_iff.2 (by linarith)
  refine tendsto_integral_of_dominated_convergence (fun u => 2 * Real.exp (-u)) ?_ ?_ ?_ ?_
  · intro J
    exact (continuous_finsetSum _ (fun n _ => tn_cont w n)).aestronglyMeasurable
  · have : IntegrableOn (fun u : ℝ => Real.exp (-u)) (Ioi 0) := by
      simpa using exp_neg_integrableOn_Ioi 0 one_pos
    exact this.const_mul 2
  · intro J
    refine (ae_restrict_iff' measurableSet_Ioi).2 (Filter.Eventually.of_forall fun u hu => ?_)
    have hu0 : 0 < u := hu
    set x := w * Real.exp (-u)
    have he0 : 0 < Real.exp (-u) := Real.exp_pos _
    have he1 : Real.exp (-u) ≤ 1 := Real.exp_le_one_iff.2 (by linarith)
    have hx0 : 0 < x := by positivity
    have hx1 : x ≤ Real.exp (-u) := by nlinarith
    have hx1' : x ≤ 1 := by linarith
    have hs : ∑ n ∈ Finset.range J, tn w u n =
        √(1 - Real.exp (-u) ^ 2) * x * ∑ n ∈ Finset.range J, (-x ^ 2) ^ n := by
      rw [Finset.mul_sum]; exact Finset.sum_congr rfl (fun n _ => by simp only [tn, x]; ring)
    have hne : -x ^ 2 ≠ 1 := by nlinarith
    have hg := geom_sum_eq hne J
    have hb : |∑ n ∈ Finset.range J, (-x ^ 2) ^ n| ≤ 2 := by
      rw [hg, abs_div]
      have hyJ : |(-x ^ 2) ^ J| ≤ 1 := by
        rw [abs_pow]; apply pow_le_one₀ (abs_nonneg _)
        rw [abs_neg, abs_of_nonneg (by positivity)]; nlinarith
      have h1 : |(-x ^ 2) ^ J - 1| ≤ 2 := by
        calc |(-x ^ 2) ^ J - 1| ≤ |(-x ^ 2) ^ J| + |(1 : ℝ)| := abs_sub _ _
          _ ≤ 2 := by rw [abs_one]; linarith
      have h2 : 1 ≤ |-x ^ 2 - 1| := by
        rw [abs_of_neg (by nlinarith)]; nlinarith
      rw [div_le_iff₀ (by linarith)]; nlinarith
    have hsq : √(1 - Real.exp (-u) ^ 2) ≤ 1 := by
      rw [Real.sqrt_le_one]; nlinarith [sq_nonneg (Real.exp (-u))]
    have hsq0 := Real.sqrt_nonneg (1 - Real.exp (-u) ^ 2)
    rw [Real.norm_eq_abs, hs, abs_mul, abs_mul, abs_of_nonneg hsq0, abs_of_pos hx0]
    calc √(1 - Real.exp (-u) ^ 2) * x * |∑ n ∈ Finset.range J, (-x ^ 2) ^ n|
        ≤ 1 * Real.exp (-u) * 2 := by gcongr
      _ = 2 * Real.exp (-u) := by ring
  · refine (ae_restrict_iff' measurableSet_Ioi).2 (Filter.Eventually.of_forall fun u hu => ?_)
    have hu0 : 0 < u := hu
    set x := w * Real.exp (-u)
    have he0 : 0 < Real.exp (-u) := Real.exp_pos _
    have he1 : Real.exp (-u) < 1 := Real.exp_lt_one_iff.2 (by linarith)
    have hx0 : 0 < x := by positivity
    have hx1 : x < 1 := by nlinarith
    have hs : ∀ J, ∑ n ∈ Finset.range J, tn w u n =
        √(1 - Real.exp (-u) ^ 2) * x * ∑ n ∈ Finset.range J, (-x ^ 2) ^ n := by
      intro J
      rw [Finset.mul_sum]; exact Finset.sum_congr rfl (fun n _ => by simp only [tn, x]; ring)
    simp_rw [hs]
    have habs : |-x ^ 2| < 1 := by rw [abs_neg, abs_of_nonneg (by positivity)]; nlinarith
    have hgeo := (hasSum_geometric_of_abs_lt_one habs).tendsto_sum_nat
    have := hgeo.const_mul (√(1 - Real.exp (-u) ^ 2) * x)
    convert this using 2
    simp only [Gint, x]
    have h1 : 1 + w ^ 2 * Real.exp (-u) ^ 2 ≠ 0 := by positivity
    have h2 : (1 - -x ^ 2) = 1 + w ^ 2 * Real.exp (-u) ^ 2 := by simp only [x]; ring
    rw [h2, div_eq_mul_inv]


lemma memR (R x : ℝ) (hR : 1 ≤ R) (hx : x ∈ Icc (1 / 2 : ℝ) R) :
    (x : ℂ) ∈ Metric.ball ((2 * R ^ 2 : ℝ) : ℂ) (2 * R ^ 2 - 1 / 4) := by
  have := mem_ball_aux x 0 (2 * R ^ 2) (by nlinarith) (by nlinarith [hx.1, hx.2])
  simpa using this

lemma memV (R y : ℝ) (hR : 1 ≤ R) (hy : y ∈ Icc (0 : ℝ) (R - 1 / 2)) :
    ((R : ℂ) + (y : ℂ) * Complex.I) ∈ Metric.ball ((2 * R ^ 2 : ℝ) : ℂ) (2 * R ^ 2 - 1 / 4) :=
  mem_ball_aux R y (2 * R ^ 2) (by nlinarith) (by nlinarith [hy.1, hy.2])

lemma memL (R r : ℝ) (hR : 1 ≤ R) (hr : r ∈ Icc (0 : ℝ) ((R - 1 / 2) * √2)) :
    Lpt r ∈ Metric.ball ((2 * R ^ 2 : ℝ) : ℂ) (2 * R ^ 2 - 1 / 4) := by
  have h22 := sqrt2_sq
  have hs2 : (0 : ℝ) < √2 := by positivity
  rw [Lpt_eq]
  set u := r * (√2 / 2)
  have hu0 : 0 ≤ u := by have := hr.1; positivity
  have hu1 : u ≤ R - 1 / 2 := by
    have := hr.2
    have : r * (√2 / 2) ≤ (R - 1 / 2) * √2 * (√2 / 2) :=
      mul_le_mul_of_nonneg_right this (by positivity)
    rw [sq2half] at this
    exact this
  apply mem_ball_aux _ _ (2 * R ^ 2) (by nlinarith)
  have h1 : u * u ≤ u * R := mul_le_mul_of_nonneg_left (by linarith) hu0
  have h2 : u * R ≤ u * R ^ 2 := mul_le_mul_of_nonneg_left (by nlinarith) hu0
  nlinarith

lemma ftc_L (F f : ℂ → ℂ) (T : ℝ) (hT : 0 ≤ T)
    (hd : ∀ r ∈ Icc (0 : ℝ) T, HasDerivAt F (f (Lpt r)) (Lpt r))
    (hc : ContinuousOn (fun r : ℝ => f (Lpt r)) (Icc 0 T)) :
    ∫ r in (0 : ℝ)..T, f (Lpt r) * ω = F (Lpt T) - F (Lpt 0) := by
  apply intervalIntegral.integral_eq_sub_of_hasDerivAt (f := fun r : ℝ => F (Lpt r))
  · intro r hr
    rw [uIcc_of_le hT] at hr
    have h1 : HasDerivAt (fun z : ℂ => (1 / 2 : ℂ) + z * ω) ω (r : ℂ) := by
      simpa using ((hasDerivAt_id (r : ℂ)).mul_const ω).const_add (1 / 2 : ℂ)
    have hFr : HasDerivAt F (f (Lpt r)) ((1 / 2 : ℂ) + (r : ℂ) * ω) := hd r hr
    exact (hFr.comp (r : ℂ) h1).comp_ofReal
  · apply ContinuousOn.intervalIntegrable
    rw [uIcc_of_le hT]
    exact hc.mul continuousOn_const

lemma ftc_V (F f : ℂ → ℂ) (R Y : ℝ) (hY : 0 ≤ Y)
    (hd : ∀ y ∈ Icc (0 : ℝ) Y, HasDerivAt F (f (R + y * Complex.I)) (R + y * Complex.I))
    (hc : ContinuousOn (fun y : ℝ => f (R + y * Complex.I)) (Icc 0 Y)) :
    ∫ y in (0 : ℝ)..Y, f (R + y * Complex.I) * Complex.I = F (R + Y * Complex.I) - F R := by
  have := intervalIntegral.integral_eq_sub_of_hasDerivAt
    (f := fun y : ℝ => F (R + y * Complex.I)) (a := 0) (b := Y)
    (f' := fun y : ℝ => f (R + y * Complex.I) * Complex.I) ?_ ?_
  · simpa using this
  · intro y hy
    rw [uIcc_of_le hY] at hy
    have h1 : HasDerivAt (fun z : ℂ => (R : ℂ) + z * Complex.I) Complex.I (y : ℂ) := by
      simpa using ((hasDerivAt_id (y : ℂ)).mul_const Complex.I).const_add (R : ℂ)
    exact ((hd y hy).comp (y : ℂ) h1).comp_ofReal
  · apply ContinuousOn.intervalIntegrable
    rw [uIcc_of_le hY]
    exact hc.mul continuousOn_const

lemma ftc_R (F f : ℂ → ℂ) (a b : ℝ) (hab : a ≤ b)
    (hd : ∀ x ∈ Icc a b, HasDerivAt F (f x) x)
    (hc : ContinuousOn (fun x : ℝ => f x) (Icc a b)) :
    ∫ x in a..b, f x = F b - F a := by
  apply intervalIntegral.integral_eq_sub_of_hasDerivAt (f := fun x : ℝ => F x)
  · intro x hx
    rw [uIcc_of_le hab] at hx
    exact (hd x hx).comp_ofReal
  · apply ContinuousOn.intervalIntegrable
    rw [uIcc_of_le hab]
    exact hc


/-- The integrand `H K e^{-sz}`. -/
noncomputable def fK (z : ℝ) (s : ℂ) : ℂ := H s * K s * Complex.exp (-s * z)

/-- The residue of `fK` at an odd integer `p`. -/
noncomputable def cres (z : ℝ) (p : ℂ) : ℂ :=
  -(H p * Complex.exp (-p * z)) / (π * Complex.sin (π * p / 2))

lemma fK_diff (z : ℝ) (s : ℂ) (hs : 0 < s.re) (hodd : ∀ n : ℤ, s ≠ 2 * n + 1) :
    DifferentiableAt ℂ (fK z) s := by
  have hop : IsOpen {s : ℂ | 0 < s.re} := isOpen_lt continuous_const Complex.continuous_re
  have hH : DifferentiableAt ℂ H s := (H_diff s hs).differentiableAt (hop.mem_nhds hs)
  have hK : DifferentiableAt ℂ K s := Q05.proof_6.1 s hodd
  have hE : DifferentiableAt ℂ (fun s : ℂ => Complex.exp (-s * z)) s := by fun_prop
  exact (hH.mul hK).mul hE

lemma fK_res (z : ℝ) (n : ℕ) :
    Tendsto (fun s => (s - (2 * (n : ℂ) + 1)) * fK z s) (𝓝[≠] (2 * (n : ℂ) + 1))
      (𝓝 (cres z (2 * (n : ℂ) + 1))) := by
  have hK := K_residue (n : ℤ)
  push_cast at hK
  set p : ℂ := 2 * (n : ℂ) + 1
  have hre : 0 < p.re := by simp [p]; positivity
  have hop : IsOpen {s : ℂ | 0 < s.re} := isOpen_lt continuous_const Complex.continuous_re
  have hφ : ContinuousAt (fun s => H s * Complex.exp (-s * z)) p := by
    have hH : DifferentiableAt ℂ H p := (H_diff p hre).differentiableAt (hop.mem_nhds hre)
    have hE : DifferentiableAt ℂ (fun s : ℂ => Complex.exp (-s * z)) p := by fun_prop
    exact (hH.mul hE).continuousAt
  have ht := (hφ.tendsto.mono_left nhdsWithin_le_nhds).mul hK
  have hsin : Complex.sin (π * p / 2) = (-1 : ℂ) ^ n := by
    have := sin_odd (n : ℤ); push_cast at this; simpa [p] using this
  have hval : H p * Complex.exp (-p * z) * ((-1 : ℂ) ^ ((n : ℤ) + 1) / π) = cres z p := by
    rw [cres, hsin, zpow_add_one₀ (by norm_num), zpow_natCast]
    have hp : (π : ℂ) ≠ 0 := Complex.ofReal_ne_zero.2 Real.pi_ne_zero
    have h1 : ((-1 : ℂ) ^ n) ^ 2 = 1 := by rw [← pow_mul, mul_comm, pow_mul]; norm_num
    have h2 : ((-1 : ℂ) ^ n) ≠ 0 := pow_ne_zero _ (by norm_num)
    field_simp
    linear_combination (-(H p)) * h1
  rw [← hval]
  refine ht.congr' (Filter.Eventually.of_forall fun s => ?_)
  simp only [fK]; ring

lemma ofReal_H (x : ℝ) (hx : 0 < x) : ∃ r : ℝ, H x = r := by
  refine ⟨1 / 2 * (Real.Gamma (x / 2) * Real.Gamma (3 / 2) / Real.Gamma (x / 2 + 3 / 2)), ?_⟩
  rw [H_eq_Gamma _ (by simpa using hx)]
  have e1 : (x : ℂ) / 2 = ((x / 2 : ℝ) : ℂ) := by push_cast; ring
  have e2 : (x : ℂ) / 2 + 3 / 2 = ((x / 2 + 3 / 2 : ℝ) : ℂ) := by push_cast; ring
  have e3 : (3 / 2 : ℂ) = ((3 / 2 : ℝ) : ℂ) := by push_cast; ring
  rw [e2, e1, e3, Complex.Gamma_ofReal, Complex.Gamma_ofReal, Complex.Gamma_ofReal]
  push_cast; ring

lemma fK_real (z x : ℝ) (hx : 0 < x) : (fK z x).im = 0 := by
  obtain ⟨r, hr⟩ := ofReal_H x hx
  have : fK z x = ((r * (1 / (2 * Real.cos (π * x / 2))) * Real.exp (-x * z) : ℝ) : ℂ) := by
    simp only [fK, K, hr]; push_cast; ring_nf
  rw [this, Complex.ofReal_im]

lemma K_vert (J : ℕ) (y : ℝ) : ‖K ((2 * J : ℝ) + y * Complex.I)‖ ≤ Real.exp (-(π / 2 * y)) := by
  have e : (π : ℂ) * (((2 * J : ℝ) : ℂ) + y * Complex.I) / 2 =
      ((J * π : ℝ) : ℂ) + ((π / 2 * y : ℝ) : ℂ) * Complex.I := by push_cast; ring
  rw [K, e, Complex.cos_add_mul_I, ← Complex.ofReal_cos, ← Complex.ofReal_sin,
    Real.sin_nat_mul_pi, Real.cos_nat_mul_pi, ← Complex.ofReal_cosh]
  simp only [Complex.ofReal_zero, zero_mul, sub_zero]
  set b := π / 2 * y
  have hc : Real.exp b / 2 ≤ Real.cosh b := by
    rw [Real.cosh_eq]; linarith [Real.exp_pos (-b)]
  have hcp : 0 < Real.cosh b := Real.cosh_pos b
  rw [norm_div, norm_one, norm_mul, norm_mul, Complex.norm_two, Complex.norm_real,
    Complex.norm_real, Real.norm_eq_abs, Real.norm_eq_abs, abs_pow, abs_neg, abs_one, one_pow,
    one_mul, abs_of_pos hcp, div_le_iff₀ (by positivity)]
  have : Real.exp (-b) * Real.exp b = 1 := by rw [← Real.exp_add]; simp
  nlinarith [Real.exp_pos (-b)]


lemma not_odd_of_im (s : ℂ) (h : s.im ≠ 0) : ∀ n : ℤ, s ≠ 2 * n + 1 := by
  intro n hn; apply h; rw [hn]; simp

lemma exp_le_one' (z : ℝ) (hz : 0 ≤ z) (s : ℂ) (hs : 0 ≤ s.re) : ‖Complex.exp (-s * z)‖ ≤ 1 := by
  rw [Complex.norm_exp]
  apply Real.exp_le_one_iff.2
  simp; nlinarith

lemma vert_bd (z : ℝ) (hz : 0 ≤ z) (J : ℕ) (hJ : 1 ≤ J) :
    ‖∫ y in (0 : ℝ)..((2 * J : ℝ) - 1 / 2), fK z ((2 * J : ℝ) + y * Complex.I) * Complex.I‖ ≤
      1 / (2 * J) := by
  have hJ' : (1 : ℝ) ≤ J := by exact_mod_cast hJ
  set R : ℝ := 2 * J with hRdef
  have hR : 2 ≤ R := by linarith
  have hY : 0 ≤ R - 1 / 2 := by linarith
  have hb : ∀ y : ℝ, 0 ≤ y →
      ‖fK z (R + y * Complex.I) * Complex.I‖ ≤ 1 / R * Real.exp (-(π / 2 * y)) := by
    intro y hy
    set s : ℂ := (R : ℂ) + (y : ℂ) * Complex.I
    have hre : s.re = R := by simp [s]
    have hH := H_le_one_div s (by rw [hre]; linarith)
    rw [hre] at hH
    have hK := K_vert J y
    have hE := exp_le_one' z hz s (by rw [hre]; linarith)
    rw [norm_mul, Complex.norm_I, mul_one, fK, norm_mul, norm_mul]
    calc ‖H s‖ * ‖K s‖ * ‖Complex.exp (-s * z)‖ ≤ (1 / R) * Real.exp (-(π / 2 * y)) * 1 := by
          gcongr
      _ = _ := mul_one _
  have h1 := intervalIntegral.norm_integral_le_of_norm_le (μ := volume) hY
    (Filter.Eventually.of_forall fun y hy => hb y (le_of_lt hy.1))
    (g := fun y => 1 / R * Real.exp (-(π / 2 * y))) (by
      apply Continuous.intervalIntegrable; fun_prop)
  refine h1.trans ?_
  have hpi := Real.pi_gt_three
  have hint : ∫ y in (0 : ℝ)..(R - 1 / 2), Real.exp (-(π / 2 * y)) =
      (1 - Real.exp (-(π / 2 * (R - 1 / 2)))) / (π / 2) := by
    rw [intervalIntegral.integral_eq_sub_of_hasDerivAt
      (f := fun y => -Real.exp (-(π / 2 * y)) / (π / 2))]
    · simp; field_simp; ring
    · intro y _
      have h0 : HasDerivAt (fun y : ℝ => -(π / 2 * y)) (-(π / 2)) y := by
        simpa using ((hasDerivAt_id y).const_mul (π / 2)).fun_neg
      have := (h0.exp.fun_neg).div_const (π / 2)
      convert this using 1
      field_simp
    · apply Continuous.intervalIntegrable; fun_prop
  rw [intervalIntegral.integral_const_mul, hint]
  have he : 0 < Real.exp (-(π / 2 * (R - 1 / 2))) := Real.exp_pos _
  have hR0 : 0 < R := by linarith
  rw [div_mul_eq_mul_div, one_mul, div_le_div_iff₀ (by positivity) hR0, one_mul, div_mul_eq_mul_div,
    div_le_iff₀ (by positivity)]
  nlinarith

lemma int_part (z : ℝ) (hz : 0 ≤ z) :
    IntegrableOn (fun r : ℝ => H (Lpt r) * K (Lpt r) * Complex.exp (-(Lpt r) * z)) (Ioi 0) := by
  have hLc : Continuous Lpt := by unfold Lpt; fun_prop
  have hc : ContinuousOn (fun r : ℝ => fK z (Lpt r)) (Ioi 0) := by
    intro r hr
    have hr0 : 0 < r := hr
    have him : (Lpt r).im ≠ 0 := by rw [Lpt_im]; positivity
    exact ((fK_diff z (Lpt r) (Lpt_re_pos r hr0.le) (not_odd_of_im _ him)).continuousAt.comp
      hLc.continuousAt).continuousWithinAt
  have hg : IntegrableOn (fun r : ℝ => 4 * Real.exp (-a * r)) (Ioi 0) :=
    (exp_neg_integrableOn_Ioi 0 a_pos).const_mul 4
  refine Integrable.mono' hg (hc.aestronglyMeasurable measurableSet_Ioi) ?_
  refine (ae_restrict_iff' measurableSet_Ioi).2 (Filter.Eventually.of_forall fun r hr => ?_)
  have hr0 : 0 ≤ r := le_of_lt hr
  obtain ⟨hH, -, hK, -⟩ := stmt7_at r hr0
  have hE := exp_le_one' z hz (Lpt r) (Lpt_re_pos r hr0).le
  simp only [norm_mul]
  calc ‖H (Lpt r)‖ * ‖K (Lpt r)‖ * ‖Complex.exp (-Lpt r * z)‖ ≤
        2 * (2 * Real.exp (-a * r)) * 1 := by gcongr
    _ = 4 * Real.exp (-a * r) := by ring

lemma cres_real (z : ℝ) (n : ℕ) :
    cres z (2 * (n : ℂ) + 1) = ((-(∫ u in Ioi (0 : ℝ), tn (Real.exp (-z)) u n) / π : ℝ) : ℂ) := by
  have ht := term_eq z n
  set p : ℂ := 2 * (n : ℂ) + 1
  have hsin : Complex.sin (π * p / 2) = (-1 : ℂ) ^ n := by
    have := sin_odd (n : ℤ); push_cast at this; simpa [p] using this
  have h1 : ((-1 : ℂ) ^ n) ^ 2 = 1 := by rw [← pow_mul, mul_comm, pow_mul]; norm_num
  have h2 : ((-1 : ℂ) ^ n) ≠ 0 := pow_ne_zero _ (by norm_num)
  have hp : (π : ℂ) ≠ 0 := Complex.ofReal_ne_zero.2 Real.pi_ne_zero
  have hφ : H p * Complex.exp (-p * z) =
      (-1 : ℂ) ^ n * ((∫ u in Ioi (0 : ℝ), tn (Real.exp (-z)) u n : ℝ) : ℂ) := by
    rw [← ht, ← mul_assoc, ← sq, h1, one_mul]
  rw [cres, hsin, hφ]
  push_cast
  field_simp

lemma arg_half (n : ℕ) :
    Complex.arg (-Complex.I * (((1 / 2 : ℝ) : ℂ) - (2 * (n : ℂ) + 1))) = π / 2 := by
  rw [Complex.arg_eq_pi_div_two_iff]
  constructor
  · simp
  · simp
    have : (0 : ℝ) ≤ n := Nat.cast_nonneg n
    linarith

lemma arg_R (n J : ℕ) :
    Complex.arg (-Complex.I * (((2 * J : ℝ) : ℂ) - (2 * (n : ℂ) + 1))) =
      if n < J then -(π / 2) else π / 2 := by
  split_ifs with h
  · rw [Complex.arg_eq_neg_pi_div_two_iff]
    have : (n : ℝ) + 1 ≤ J := by exact_mod_cast h
    constructor
    · simp
    · simp; linarith
  · rw [Complex.arg_eq_pi_div_two_iff]
    have : (J : ℝ) ≤ n := by exact_mod_cast (not_lt.1 h)
    constructor
    · simp
    · simp; linarith

lemma pole_im (z : ℝ) (n J : ℕ) :
    (cres z (2 * (n : ℂ) + 1) *
      (Complex.log (-Complex.I * (((2 * J : ℝ) : ℂ) - (2 * (n : ℂ) + 1))) -
        Complex.log (-Complex.I * (((1 / 2 : ℝ) : ℂ) - (2 * (n : ℂ) + 1))))).im =
      if n < J then ∫ u in Ioi (0 : ℝ), tn (Real.exp (-z)) u n else 0 := by
  rw [cres_real, Complex.im_ofReal_mul, Complex.sub_im, Complex.log_im, Complex.log_im, arg_half,
    arg_R]
  have hp : π ≠ 0 := Real.pi_ne_zero
  split_ifs
  · field_simp; ring
  · ring

lemma contour_J (z : ℝ) (hz : 0 ≤ z) (J : ℕ) (hJ : 1 ≤ J) :
    (∫ r in (0 : ℝ)..(((2 * J : ℝ) - 1 / 2) * √2), fK z (Lpt r) * ω).im =
      (∫ y in (0 : ℝ)..((2 * J : ℝ) - 1 / 2),
          fK z ((2 * J : ℝ) + y * Complex.I) * Complex.I).im +
        ∫ u in Ioi (0 : ℝ), ∑ n ∈ Finset.range J, tn (Real.exp (-z)) u n := by
  classical
  have hJ' : (1 : ℝ) ≤ J := by exact_mod_cast hJ
  set R : ℝ := 2 * J with hRdef
  have hR : 1 ≤ R := by linarith
  have hY : 0 ≤ R - 1 / 2 := by linarith
  set M : ℕ := 8 * J ^ 2 with hMdef
  have hMX : (M : ℝ) = 2 * R ^ 2 := by simp only [hMdef, hRdef]; push_cast; ring
  set B := Metric.ball ((2 * R ^ 2 : ℝ) : ℂ) (2 * R ^ 2 - 1 / 4) with hBdef
  set pt : ℕ → ℂ := fun n => 2 * (n : ℂ) + 1 with hpt
  set S : Finset ℂ := (Finset.range M).image pt with hSdef
  have hptinj : Function.Injective pt := by
    intro a b h
    simp only [hpt] at h
    have : (a : ℂ) = b := by linear_combination h / 2
    exact_mod_cast this
  have hpt_re : ∀ n, (pt n).re = 2 * n + 1 := fun n => by simp [hpt]
  have hpt_im : ∀ n, (pt n).im = 0 := fun n => by simp [hpt]
  -- the poles in the ball
  have hSB : (↑S : Set ℂ) ⊆ B := by
    intro p hp
    rw [hSdef, Finset.coe_image, Set.mem_image] at hp
    obtain ⟨n, hn, rfl⟩ := hp
    have hn' : (n : ℝ) + 1 ≤ M := by exact_mod_cast (Finset.mem_range.1 hn)
    have hn0 : (0 : ℝ) ≤ n := Nat.cast_nonneg n
    have h1 : 0 < 2 * (2 * R ^ 2) - 1 / 4 - (2 * n + 1) := by linarith
    have h2 : 0 < (2 * (n : ℝ) + 1) - 1 / 4 := by linarith
    have := mem_ball_aux (2 * n + 1) 0 (2 * R ^ 2) (by nlinarith) (by nlinarith [mul_pos h1 h2])
    convert this using 1
    simp [hpt]
  have hodd : ∀ s ∈ B, s ∉ S → ∀ m : ℤ, s ≠ 2 * m + 1 := by
    intro s hs hsS m hm
    apply hsS
    have hsB := hs
    rw [hm, hBdef, Metric.mem_ball, dist_eq_norm] at hsB
    have e : (2 * (m : ℂ) + 1) - ((2 * R ^ 2 : ℝ) : ℂ) = (((2 * m + 1 - 2 * R ^ 2 : ℝ)) : ℂ) := by
      push_cast; ring
    rw [e, Complex.norm_real, Real.norm_eq_abs] at hsB
    have h1 := abs_lt.1 hsB
    have hm1 : (-1 : ℝ) < m := by linarith [h1.1]
    have hm1' : -1 < m := by exact_mod_cast hm1
    have hm0 : 0 ≤ m := by omega
    have hmM : (m : ℝ) < ((M : ℤ) : ℝ) := by push_cast; linarith [h1.2]
    have hmM' : m < (M : ℤ) := by exact_mod_cast hmM
    rw [hSdef, Finset.mem_image]
    refine ⟨m.toNat, Finset.mem_range.2 (by omega), ?_⟩
    have : ((m.toNat : ℕ) : ℂ) = (m : ℂ) := by
      rw [← Int.cast_natCast, Int.toNat_of_nonneg hm0]
    simp only [hpt, this, hm]
  have hdiff : DifferentiableOn ℂ (fK z) (B \ ↑S) := fun s hs =>
    (fK_diff z s (ball_sub _ hs.1) (hodd s hs.1 hs.2)).differentiableWithinAt
  have hres : ∀ p ∈ S, Tendsto (fun s => (s - p) * fK z s) (𝓝[≠] p) (𝓝 (cres z p)) := by
    intro p hp
    obtain ⟨n, -, rfl⟩ := Finset.mem_image.1 hp
    exact fK_res z n
  obtain ⟨g, hgd, hgeq⟩ := glue S B Metric.isOpen_ball hSB (fK z) (cres z) hdiff hres
  obtain ⟨G, hG⟩ := hgd.isExactOn_ball
  set P : ℂ → ℂ := fun s => ∑ p ∈ S, cres z p * Complex.log (-Complex.I * (s - p)) with hPdef
  have hP : ∀ s, (∀ p ∈ S, -Complex.I * (s - p) ∈ Complex.slitPlane) →
      HasDerivAt P (∑ p ∈ S, cres z p / (s - p)) s := by
    intro s hs
    refine HasDerivAt.fun_sum (fun p hp => ?_)
    have h1 : HasDerivAt (fun s => -Complex.I * (s - p)) (-Complex.I * 1) s :=
      ((hasDerivAt_id s).sub_const p).const_mul (-Complex.I)
    have h2 := (h1.clog (hs p hp)).const_mul (cres z p)
    convert h2 using 1
    have hsp : s - p ≠ 0 := by
      intro h0; have := Complex.slitPlane_ne_zero (hs p hp); rw [h0, mul_zero] at this
      exact this rfl
    field_simp
  have hGP : ∀ s ∈ B \ ↑S, (∀ p ∈ S, -Complex.I * (s - p) ∈ Complex.slitPlane) →
      HasDerivAt (fun s => G s + P s) (fK z s) s := by
    intro s hs hsl
    have := (hG s hs.1).add (hP s hsl)
    convert this using 1
    rw [hgeq s hs]; ring
  have hs2 : (0 : ℝ) < √2 := by positivity
  set T : ℝ := (R - 1 / 2) * √2 with hTdef
  have hT : 0 ≤ T := by positivity
  have hLc : Continuous Lpt := by unfold Lpt; fun_prop
  have hLS : ∀ r ∈ Icc (0 : ℝ) T, Lpt r ∈ B \ ↑S ∧
      ∀ p ∈ S, -Complex.I * (Lpt r - p) ∈ Complex.slitPlane := by
    intro r hr
    have hr0 : 0 ≤ r := hr.1
    refine ⟨⟨memL R r hR hr, ?_⟩, ?_⟩
    · intro hmem
      obtain ⟨n, -, hn⟩ := Finset.mem_image.1 (Finset.mem_coe.1 hmem)
      have him := congrArg Complex.im hn
      have hre := congrArg Complex.re hn
      rw [hpt_im, Lpt_im] at him
      rw [hpt_re, Lpt_re] at hre
      have hr0' : r = 0 := by
        rcases (mul_eq_zero.1 him.symm) with h | h
        · exact h
        · linarith
      rw [hr0'] at hre
      have : (0 : ℝ) ≤ n := Nat.cast_nonneg n
      linarith
    · intro p hp
      obtain ⟨n, -, rfl⟩ := Finset.mem_image.1 hp
      rw [Complex.mem_slitPlane_iff]
      have hre : (-Complex.I * (Lpt r - pt n)).re = r * (√2 / 2) := by
        simp [Complex.mul_re, Lpt_im, hpt_im]
      have him : (-Complex.I * (Lpt r - pt n)).im = (2 * n + 1) - (1 / 2 + r * (√2 / 2)) := by
        simp [Complex.mul_im, Lpt_re, hpt_re]
      rcases hr0.eq_or_lt with h | h
      · right; rw [him, ← h]
        have : (0 : ℝ) ≤ n := Nat.cast_nonneg n
        intro h'; linarith
      · left; rw [hre]; positivity
  have hVS : ∀ y ∈ Icc (0 : ℝ) (R - 1 / 2), ((R : ℂ) + y * Complex.I) ∈ B \ ↑S ∧
      ∀ p ∈ S, -Complex.I * (((R : ℂ) + y * Complex.I) - p) ∈ Complex.slitPlane := by
    intro y hy
    have hpar : ∀ n : ℕ, R ≠ 2 * n + 1 := by
      intro n h
      have : ((2 * J : ℤ) : ℝ) = ((2 * n + 1 : ℤ) : ℝ) := by push_cast; rw [← h]
      have := Int.cast_injective this
      omega
    refine ⟨⟨memV R y hR hy, ?_⟩, ?_⟩
    · intro hmem
      obtain ⟨n, -, hn⟩ := Finset.mem_image.1 (Finset.mem_coe.1 hmem)
      have hre := congrArg Complex.re hn
      rw [hpt_re] at hre
      simp at hre
      exact hpar n hre.symm
    · intro p hp
      obtain ⟨n, -, rfl⟩ := Finset.mem_image.1 hp
      rw [Complex.mem_slitPlane_iff]
      right
      have him : (-Complex.I * (((R : ℂ) + y * Complex.I) - pt n)).im = (2 * n + 1) - R := by
        simp [Complex.mul_im, hpt_re]
      rw [him]
      intro h; exact hpar n (by linarith)
  have hcL : ContinuousOn (fun r : ℝ => fK z (Lpt r)) (Icc 0 T) := by
    intro r hr
    have h := hLS r hr
    exact ((fK_diff z _ (ball_sub _ h.1.1) (hodd _ h.1.1 h.1.2)).continuousAt.comp
      hLc.continuousAt).continuousWithinAt
  have hVc : Continuous (fun y : ℝ => (R : ℂ) + y * Complex.I) := by fun_prop
  have hcV : ContinuousOn (fun y : ℝ => fK z ((R : ℂ) + y * Complex.I)) (Icc 0 (R - 1 / 2)) := by
    intro y hy
    have h := hVS y hy
    exact (ContinuousAt.comp (f := fun y : ℝ => (R : ℂ) + y * Complex.I)
      (fK_diff z _ (ball_sub _ h.1.1) (hodd _ h.1.1 h.1.2)).continuousAt
      hVc.continuousAt).continuousWithinAt
  have eL := ftc_L (fun s => G s + P s) (fK z) T hT
    (fun r hr => hGP _ (hLS r hr).1 (hLS r hr).2) hcL
  have eV := ftc_V (fun s => G s + P s) (fK z) R (R - 1 / 2) hY
    (fun y hy => hGP _ (hVS y hy).1 (hVS y hy).2) hcV
  have hgc : ContinuousOn (fun x : ℝ => g x) (Icc (1 / 2) R) :=
    hgd.continuousOn.comp Complex.continuous_ofReal.continuousOn (fun x hx => memR R x hR hx)
  have eR := ftc_R G g (1 / 2) R (by linarith) (fun x hx => hG _ (memR R x hR hx)) hgc
  rw [hTdef, Lpt_T] at eL
  have h0 : Lpt 0 = ((1 / 2 : ℝ) : ℂ) := by rw [Lpt]; push_cast; ring
  rw [h0] at eL
  have key : ∫ r in (0 : ℝ)..((R - 1 / 2) * √2), fK z (Lpt r) * ω =
      (∫ y in (0 : ℝ)..(R - 1 / 2), fK z ((R : ℂ) + y * Complex.I) * Complex.I) +
        (∫ x in (1 / 2 : ℝ)..R, g x) + (P R - P ((1 / 2 : ℝ) : ℂ)) := by
    rw [eL, eV, eR]; ring
  -- the real-axis integral is real
  have him_g : (∫ x in (1 / 2 : ℝ)..R, g x).im = 0 := by
    have hE : {x : ℝ | (x : ℂ) ∈ S}.Finite :=
      S.finite_toSet.preimage Complex.ofReal_injective.injOn
    have hae : ∀ᵐ x ∂volume, x ∈ Set.uIoc (1 / 2 : ℝ) R → g x = ((g x).re : ℂ) := by
      filter_upwards [measure_eq_zero_iff_ae_notMem.1 (hE.measure_zero volume)] with x hx hxI
      rw [Set.uIoc_of_le (by linarith)] at hxI
      have hxI' : x ∈ Icc (1 / 2 : ℝ) R := Ioc_subset_Icc_self hxI
      have hxS : (x : ℂ) ∉ S := hx
      have hgx := hgeq x ⟨memR R x hR hxI', hxS⟩
      apply Complex.ext
      · simp
      · rw [Complex.ofReal_im, hgx, Complex.sub_im, fK_real z x (by linarith [hxI'.1]),
          Complex.im_sum]
        rw [Finset.sum_eq_zero, sub_zero]
        intro p hp
        obtain ⟨n, -, rfl⟩ := Finset.mem_image.1 hp
        rw [cres_real]
        have : (x : ℂ) - pt n = ((x - (2 * n + 1) : ℝ) : ℂ) := by simp [hpt]
        rw [this, ← Complex.ofReal_div, Complex.ofReal_im]
    rw [intervalIntegral.integral_congr_ae hae, intervalIntegral.integral_ofReal,
      Complex.ofReal_im]
  -- the principal parts
  have hJM : J ≤ M := by simp only [hMdef]; nlinarith
  obtain ⟨k, hk⟩ : ∃ k, M = J + k := ⟨M - J, by omega⟩
  have hP_im : (P R - P ((1 / 2 : ℝ) : ℂ)).im =
      ∫ u in Ioi (0 : ℝ), ∑ n ∈ Finset.range J, tn (Real.exp (-z)) u n := by
    have hw0 : 0 < Real.exp (-z) := Real.exp_pos _
    have hw1 : Real.exp (-z) ≤ 1 := Real.exp_le_one_iff.2 (by linarith)
    rw [integral_finsetSum _ (fun n _ => tn_int _ n hw0 hw1)]
    simp only [hPdef]
    rw [← Finset.sum_sub_distrib, hSdef, Finset.sum_image (fun a _ b _ h => hptinj h),
      Complex.im_sum]
    have : ∀ n ∈ Finset.range M, (cres z (pt n) * Complex.log (-Complex.I * ((R : ℂ) - pt n)) -
        cres z (pt n) * Complex.log (-Complex.I * (((1 / 2 : ℝ) : ℂ) - pt n))).im =
        if n < J then ∫ u in Ioi (0 : ℝ), tn (Real.exp (-z)) u n else 0 := by
      intro n _
      rw [← mul_sub]
      exact pole_im z n J
    rw [Finset.sum_congr rfl this, hk, Finset.sum_range_add]
    have h1 : ∑ n ∈ Finset.range J,
        (if n < J then ∫ u in Ioi (0 : ℝ), tn (Real.exp (-z)) u n else 0) =
        ∑ n ∈ Finset.range J, ∫ u in Ioi (0 : ℝ), tn (Real.exp (-z)) u n :=
      Finset.sum_congr rfl (fun n hn => by simp [Finset.mem_range.1 hn])
    have h2 : ∑ n ∈ Finset.range k,
        (if J + n < J then ∫ u in Ioi (0 : ℝ), tn (Real.exp (-z)) u (J + n) else 0) = 0 :=
      Finset.sum_eq_zero (fun n _ => by simp)
    rw [h1, h2, add_zero]
  rw [key, Complex.add_im, Complex.add_im, him_g, hP_im, add_zero]


lemma lim_im (z : ℝ) (hz : 0 ≤ z) :
    (lineInt (fK z)).im = π / 2 * (Real.exp (-z) / (1 + √(1 + Real.exp (-z) ^ 2))) := by
  have hI : IntegrableOn (fun r : ℝ => fK z (Lpt r) * ω) (Ioi 0) := (int_part z hz).mul_const ω
  have hT : Tendsto (fun J : ℕ => ((2 * J : ℝ) - 1 / 2) * √2) atTop atTop := by
    have h1 : Tendsto (fun J : ℕ => (2 * J : ℝ) - 1 / 2) atTop atTop :=
      tendsto_atTop_add_const_right _ _
        (tendsto_natCast_atTop_atTop.const_mul_atTop (by norm_num : (0 : ℝ) < 2))
    exact h1.atTop_mul_const (by positivity)
  have tL : Tendsto (fun J : ℕ => (∫ r in (0 : ℝ)..(((2 * J : ℝ) - 1 / 2) * √2),
      fK z (Lpt r) * ω).im) atTop (𝓝 (lineInt (fK z)).im) :=
    (Complex.continuous_im.tendsto _).comp (intervalIntegral_tendsto_integral_Ioi 0 hI hT)
  have tV : Tendsto (fun J : ℕ => (∫ y in (0 : ℝ)..((2 * J : ℝ) - 1 / 2),
      fK z ((2 * J : ℝ) + y * Complex.I) * Complex.I).im) atTop (𝓝 0) := by
    have h0 : Tendsto (fun J : ℕ => 1 / (2 * (J : ℝ))) atTop (𝓝 0) := by
      have := (tendsto_natCast_atTop_atTop (R := ℝ)).const_mul_atTop (by norm_num : (0 : ℝ) < 2)
      refine this.inv_tendsto_atTop.congr (fun J => ?_)
      simp
    refine squeeze_zero_norm' ?_ h0
    filter_upwards [eventually_ge_atTop 1] with J hJ
    rw [Real.norm_eq_abs]
    exact (Complex.abs_im_le_norm _).trans (vert_bd z hz J hJ)
  have tS := sum_tendsto z hz
  have tR := tV.add tS
  rw [zero_add] at tR
  have heq : (fun J : ℕ => (∫ r in (0 : ℝ)..(((2 * J : ℝ) - 1 / 2) * √2),
      fK z (Lpt r) * ω).im) =ᶠ[atTop] fun J : ℕ =>
        (∫ y in (0 : ℝ)..((2 * J : ℝ) - 1 / 2),
          fK z ((2 * J : ℝ) + y * Complex.I) * Complex.I).im +
        ∫ u in Ioi (0 : ℝ), ∑ n ∈ Finset.range J, tn (Real.exp (-z)) u n := by
    filter_upwards [eventually_ge_atTop 1] with J hJ
    exact contour_J z hz J hJ
  rw [tendsto_nhds_unique tL (tR.congr' heq.symm), real_int _ (Real.exp_pos _)]

theorem main : Q07.stmt_1 := by
  intro z hz
  refine ⟨int_part z hz, ?_⟩
  have h := lim_im z hz
  have e : (lineInt fun s => H s * K s * Complex.exp (-s * z)) = lineInt (fK z) := rfl
  rw [e, h, Cf]
  have e2 : Real.exp (-2 * z) = Real.exp (-z) ^ 2 := by rw [← Real.exp_nat_mul]; ring_nf
  rw [e2]
  have hp : π ≠ 0 := Real.pi_ne_zero
  field_simp

end Q07.PA

/-! ### Merged from work/Q07/PartB.lean -/
/-!
# Q07, part B: `Q07.stmt_1 → Q07.stmt_2` (g5:lem:Mk from the kernel identity)

Plan. For real `ε` with `|ε| ≤ 1/4`, put `w_ε(x) = e^{-x}(1+e^{-x})^{-1-ε}`,
`ℓ(s,x) = log(1+e^{-x}) - s x`, `base(s,x) = e^{-x}/(1+e^{-x}) · e^{-sx}` and
`Φ(ε,s) = ∫_0^∞ base(s,x) e^{-εℓ(s,x)} dx = ∫_0^∞ w_ε(x) e^{-s(1-ε)x} dx`.
* (A) logit form + kernel identity + Fubini:
  `-(F ε)/2 = -Re Φ(ε,0) + (2/π) Im ∫_Λ H K Φ(ε,·)`.
* (B) expansion of `e^{-εℓ}`: `Φ(ε,s) = Σ ε^n P_n(s)` (change of variables `t = e^{-x}`).
* (C) dominated convergence on `Λ`: `∫_Λ HKΦ(ε,·) = Σ ε^n ∫_Λ HK P_n`.
* `P_n(0) = D_n`.
* (E) Taylor series of `FC` at `0`: `-(FC ε)/2 = Σ M_n ε^n`.
* (G) uniqueness of the coefficients of a real power series.
-/

namespace Q07.PB

open Real MeasureTheory Set Filter Q05

/-- `ℓ(s,x) = log(1+e^{-x}) - s x`. -/
noncomputable def ell (s : ℂ) (x : ℝ) : ℂ := (Real.log (1 + Real.exp (-x)) : ℂ) - s * x

/-- `e^{-x}/(1+e^{-x}) · e^{-sx}`. -/
noncomputable def base (s : ℂ) (x : ℝ) : ℂ :=
  ((Real.exp (-x) / (1 + Real.exp (-x)) : ℝ) : ℂ) * Complex.exp (-s * x)

/-- The terms of the expansion. -/
noncomputable def fn (e : ℝ) (s : ℂ) (n : ℕ) (x : ℝ) : ℂ :=
  base s x * ((-(e : ℂ) * ell s x) ^ n / (n.factorial : ℂ))

/-- `Φ(ε,s)`. -/
noncomputable def Phi (e : ℝ) (s : ℂ) : ℂ :=
  ∫ x in Ioi (0 : ℝ), base s x * Complex.exp (-(e : ℂ) * ell s x)

/-- `∫_0^∞ base · ℓ^n`. -/
noncomputable def Pa (n : ℕ) (s : ℂ) : ℂ := ∫ x in Ioi (0 : ℝ), base s x * ell s x ^ n

lemma image_exp_neg : (fun x : ℝ => Real.exp (-x)) '' Ioi 0 = Ioo 0 1 := by
  ext t
  simp only [mem_image, mem_Ioi, mem_Ioo]
  constructor
  · rintro ⟨x, hx, rfl⟩
    exact ⟨Real.exp_pos _, Real.exp_lt_one_iff.2 (by linarith)⟩
  · rintro ⟨h0, h1⟩
    refine ⟨-Real.log t, ?_, ?_⟩
    · have := Real.log_neg h0 h1; linarith
    · rw [neg_neg, Real.exp_log h0]

lemma cov (g : ℝ → ℂ) :
    ∫ t in (0 : ℝ)..1, g t = ∫ x in Ioi (0 : ℝ), (Real.exp (-x) : ℂ) * g (Real.exp (-x)) := by
  have hder : ∀ x ∈ Ioi (0 : ℝ), HasDerivWithinAt (fun x : ℝ => Real.exp (-x))
      (-Real.exp (-x)) (Ioi 0) x := by
    intro x _
    have := (hasDerivAt_neg x).exp
    exact (by simpa using this : HasDerivAt (fun x : ℝ => Real.exp (-x)) (-Real.exp (-x)) x).hasDerivWithinAt
  have hinj : InjOn (fun x : ℝ => Real.exp (-x)) (Ioi 0) := by
    intro x _ y _ h
    have := Real.exp_injective h
    linarith
  have hEq := integral_image_eq_integral_abs_deriv_smul measurableSet_Ioi hder hinj g
  rw [image_exp_neg] at hEq
  rw [intervalIntegral.integral_of_le zero_le_one, integral_Ioc_eq_integral_Ioo, hEq]
  refine setIntegral_congr_fun measurableSet_Ioi (fun x _ => ?_)
  rw [abs_neg, abs_of_pos (Real.exp_pos _), Complex.real_smul]

lemma P_integrand (n : ℕ) (s : ℂ) (x : ℝ) :
    (Real.exp (-x) : ℂ) * ((Real.exp (-x) : ℂ) ^ s / (1 + (Real.exp (-x) : ℂ)) *
      (s * (Real.log (Real.exp (-x)) : ℂ) + (Real.log (1 + Real.exp (-x)) : ℂ)) ^ n) =
    base s x * ell s x ^ n := by
  have hE : (0 : ℝ) < Real.exp (-x) := Real.exp_pos _
  have h1 : (Real.exp (-x) : ℂ) ^ s = Complex.exp (-s * x) := by
    rw [Complex.cpow_def_of_ne_zero (by exact_mod_cast hE.ne'), ← Complex.ofReal_log hE.le,
      Real.log_exp]
    push_cast; ring_nf
  have h2 : (1 + (Real.exp (-x) : ℂ)) ≠ 0 := by exact_mod_cast (by positivity : (1 + Real.exp (-x)) ≠ 0)
  rw [h1, Real.log_exp, base, ell]
  push_cast
  field_simp
  ring

lemma P_eq (n : ℕ) (s : ℂ) : P n s = (-1) ^ n / (n.factorial : ℂ) * Pa n s := by
  rw [P, cov, Pa]
  congr 1
  refine setIntegral_congr_fun measurableSet_Ioi (fun x _ => ?_)
  exact P_integrand n s x

lemma norm_base_le (s : ℂ) (x : ℝ) : ‖base s x‖ ≤ Real.exp (-x) * Real.exp (-(s.re * x)) := by
  rw [base, norm_mul, Complex.norm_real, Complex.norm_exp, Real.norm_eq_abs,
    abs_of_pos (by positivity)]
  have h1 : Real.exp (-x) / (1 + Real.exp (-x)) ≤ Real.exp (-x) := by
    rw [div_le_iff₀ (by positivity)]; nlinarith [Real.exp_pos (-x)]
  have h2 : (-s * (x : ℂ)).re = -(s.re * x) := by simp
  rw [h2]
  exact mul_le_mul_of_nonneg_right h1 (Real.exp_pos _).le

lemma log_le (x : ℝ) (hx : 0 ≤ x) :
    0 ≤ Real.log (1 + Real.exp (-x)) ∧ Real.log (1 + Real.exp (-x)) ≤ Real.log 2 := by
  have h1 : Real.exp (-x) ≤ 1 := Real.exp_le_one_iff.2 (by linarith)
  have h0 := Real.exp_pos (-x)
  exact ⟨Real.log_nonneg (by linarith), Real.log_le_log (by linarith) (by linarith)⟩

lemma norm_ell_le (s : ℂ) (x : ℝ) (hx : 0 ≤ x) : ‖ell s x‖ ≤ Real.log 2 + ‖s‖ * x := by
  obtain ⟨h0, h1⟩ := log_le x hx
  calc ‖ell s x‖ ≤ ‖(Real.log (1 + Real.exp (-x)) : ℂ)‖ + ‖s * (x : ℂ)‖ := norm_sub_le _ _
    _ ≤ Real.log 2 + ‖s‖ * x := by
      rw [Complex.norm_real, Real.norm_eq_abs, abs_of_nonneg h0, norm_mul, Complex.norm_real,
        Real.norm_eq_abs, abs_of_nonneg hx]
      linarith

lemma bnd (s : ℂ) (hs : ‖s‖ ≤ 2 * s.re) (x : ℝ) (hx : 0 ≤ x) (n : ℕ) :
    ‖base s x‖ * ‖ell s x‖ ^ n / (n.factorial : ℝ) ≤ 2 ^ (n + 1) * Real.exp (-x) := by
  have hl := norm_ell_le s x hx
  have hl0 : 0 ≤ ‖ell s x‖ := norm_nonneg _
  have hf : (0 : ℝ) < n.factorial := by exact_mod_cast Nat.factorial_pos n
  have h1 : ‖ell s x‖ ^ n / (n.factorial : ℝ) ≤ 2 ^ n * Real.exp (‖ell s x‖ / 2) := by
    have := Real.pow_div_factorial_le_exp _ (by positivity : 0 ≤ ‖ell s x‖ / 2) n
    rw [div_pow] at this
    have e : ‖ell s x‖ ^ n / (n.factorial : ℝ) = 2 ^ n * (‖ell s x‖ ^ n / 2 ^ n / n.factorial) := by
      field_simp
    rw [e]
    exact mul_le_mul_of_nonneg_left this (by positivity)
  have h2 : Real.exp (‖ell s x‖ / 2) ≤ 2 * Real.exp (s.re * x) := by
    have hl2 : Real.log 2 ≥ 0 := Real.log_nonneg (by norm_num)
    have : ‖ell s x‖ / 2 ≤ Real.log 2 + s.re * x := by nlinarith
    calc Real.exp (‖ell s x‖ / 2) ≤ Real.exp (Real.log 2 + s.re * x) := Real.exp_le_exp.2 this
      _ = 2 * Real.exp (s.re * x) := by rw [Real.exp_add, Real.exp_log (by norm_num)]
  have hb := norm_base_le s x
  have hE : Real.exp (-(s.re * x)) * Real.exp (s.re * x) = 1 := by
    rw [← Real.exp_add]; simp
  calc ‖base s x‖ * ‖ell s x‖ ^ n / (n.factorial : ℝ)
      = ‖base s x‖ * (‖ell s x‖ ^ n / (n.factorial : ℝ)) := by ring
    _ ≤ (Real.exp (-x) * Real.exp (-(s.re * x))) * (2 ^ n * (2 * Real.exp (s.re * x))) := by
      gcongr
      exact h1.trans (mul_le_mul_of_nonneg_left h2 (by positivity))
    _ = 2 ^ (n + 1) * Real.exp (-x) := by
      rw [pow_succ]
      linear_combination (2 ^ n * 2 * Real.exp (-x)) * hE

lemma fn_bnd (e : ℝ) (he : |e| ≤ 1 / 4) (s : ℂ) (hs : ‖s‖ ≤ 2 * s.re) (x : ℝ) (hx : 0 ≤ x)
    (n : ℕ) : ‖fn e s n x‖ ≤ 2 * (1 / 2) ^ n * Real.exp (-x) := by
  have hb := bnd s hs x hx n
  have hf : (0 : ℝ) < n.factorial := by exact_mod_cast Nat.factorial_pos n
  have e1 : ‖fn e s n x‖ = |e| ^ n * (‖base s x‖ * ‖ell s x‖ ^ n / (n.factorial : ℝ)) := by
    rw [fn, norm_mul, norm_div, norm_pow, norm_mul, norm_neg, Complex.norm_real,
      Real.norm_eq_abs, Complex.norm_natCast]
    ring
  rw [e1]
  calc |e| ^ n * (‖base s x‖ * ‖ell s x‖ ^ n / (n.factorial : ℝ))
      ≤ (1 / 4) ^ n * (2 ^ (n + 1) * Real.exp (-x)) := by
        gcongr
    _ = 2 * (1 / 2) ^ n * Real.exp (-x) := by
        rw [pow_succ, show (1 / 4 : ℝ) = 1 / 2 * (1 / 2) by norm_num, mul_pow]
        field_simp
        rw [← mul_pow]; norm_num

lemma fn_meas (e : ℝ) (s : ℂ) (n : ℕ) : Measurable (fn e s n) := by
  unfold fn base ell
  fun_prop

lemma fn_meas2 (e : ℝ) (n : ℕ) : Measurable (fun p : ℝ × ℝ => fn e (Lpt p.1) n p.2) := by
  unfold fn base ell Lpt
  fun_prop

lemma exp_series (e : ℝ) (s : ℂ) (x : ℝ) :
    HasSum (fun n => fn e s n x) (base s x * Complex.exp (-(e : ℂ) * ell s x)) := by
  have h := NormedSpace.expSeries_div_hasSum_exp (-(e : ℂ) * ell s x)
  rw [← Complex.exp_eq_exp_ℂ] at h
  exact h.mul_left (base s x)

lemma tsum_bound (x : ℝ) : ∑' n : ℕ, 2 * (1 / 2 : ℝ) ^ n * Real.exp (-x) = 4 * Real.exp (-x) := by
  rw [tsum_mul_right, tsum_mul_left, tsum_geometric_two]; ring

lemma phi_hasSum (e : ℝ) (he : |e| ≤ 1 / 4) (s : ℂ) (hs : ‖s‖ ≤ 2 * s.re) :
    HasSum (fun n => ∫ x in Ioi (0 : ℝ), fn e s n x) (Phi e s) := by
  refine hasSum_integral_of_dominated_convergence (fun n x => 2 * (1 / 2 : ℝ) ^ n * Real.exp (-x))
    (fun n => (fn_meas e s n).aestronglyMeasurable) (fun n => ?_) ?_ ?_ ?_
  · exact (ae_restrict_iff' measurableSet_Ioi).2 (Eventually.of_forall fun x hx =>
      fn_bnd e he s hs x (le_of_lt hx) n)
  · exact Eventually.of_forall fun x =>
      (summable_geometric_two.mul_left 2).mul_right _
  · simp only [tsum_bound]
    have := (exp_neg_integrableOn_Ioi 0 one_pos).const_mul 4
    exact IntegrableOn.congr_fun this (fun x _ => by simp) measurableSet_Ioi
  · exact Eventually.of_forall fun x => exp_series e s x

lemma int_fn (e : ℝ) (s : ℂ) (n : ℕ) :
    ∫ x in Ioi (0 : ℝ), fn e s n x = (e : ℂ) ^ n * P n s := by
  have h : ∀ x, fn e s n x = (-(e : ℂ)) ^ n / (n.factorial : ℂ) * (base s x * ell s x ^ n) := by
    intro x; rw [fn, mul_pow]; ring
  simp_rw [h]
  rw [integral_const_mul, P_eq, ← Pa, neg_pow]
  ring

lemma int_log (n : ℕ) :
    ∫ t in (0 : ℝ)..1, Real.log (1 + t) ^ n / (1 + t) = Real.log 2 ^ (n + 1) / (n + 1) := by
  have hpos : ∀ t ∈ uIcc (0 : ℝ) 1, 1 + t ≠ 0 := by
    intro t ht; rw [uIcc_of_le zero_le_one] at ht; linarith [ht.1]
  have hd : ∀ t ∈ uIcc (0 : ℝ) 1, HasDerivAt (fun t : ℝ => Real.log (1 + t) ^ (n + 1) / (n + 1))
      (Real.log (1 + t) ^ n / (1 + t)) t := by
    intro t ht
    have h1 : HasDerivAt (fun t : ℝ => 1 + t) 1 t := (hasDerivAt_id t).const_add 1
    have h2 := ((h1.log (hpos t ht)).pow (n + 1)).div_const ((n : ℝ) + 1)
    convert h2 using 1
    have : (n : ℝ) + 1 ≠ 0 := by positivity
    field_simp
    push_cast
    ring
  have hc : ContinuousOn (fun t : ℝ => Real.log (1 + t) ^ n / (1 + t)) (uIcc 0 1) := by
    apply ContinuousOn.div _ (by fun_prop) hpos
    apply ContinuousOn.pow
    exact ContinuousOn.log (by fun_prop) hpos
  rw [intervalIntegral.integral_eq_sub_of_hasDerivAt hd (hc.intervalIntegrable)]
  norm_num

lemma P_zero (n : ℕ) : P n 0 = ((D n : ℝ) : ℂ) := by
  have h : (fun t : ℝ => (t : ℂ) ^ (0 : ℂ) / (1 + (t : ℂ)) *
      ((0 : ℂ) * (Real.log t : ℂ) + (Real.log (1 + t) : ℂ)) ^ n) =
      fun t : ℝ => ((Real.log (1 + t) ^ n / (1 + t) : ℝ) : ℂ) := by
    funext t; rw [Complex.cpow_zero]; push_cast; ring
  rw [P, h, intervalIntegral.integral_ofReal, int_log, D, Nat.factorial_succ]
  push_cast
  have : (n.factorial : ℂ) ≠ 0 := by exact_mod_cast (Nat.factorial_pos n).ne'
  have : (n : ℂ) + 1 ≠ 0 := by exact_mod_cast (Nat.succ_ne_zero n)
  field_simp

/-- `w_ε(x) = e^{-x}(1+e^{-x})^{-1-ε}`. -/
noncomputable def w (e x : ℝ) : ℝ := Real.exp (-x) * (1 + Real.exp (-x)) ^ (-1 - e)

lemma phi_integrand (e : ℝ) (s : ℂ) (x : ℝ) :
    base s x * Complex.exp (-(e : ℂ) * ell s x) =
      ((w e x : ℝ) : ℂ) * Complex.exp (-s * (((1 - e) * x : ℝ) : ℂ)) := by
  set E := Real.exp (-x)
  have hL0 : 0 < 1 + E := by positivity
  set A : ℂ := ((Real.log (1 + E) : ℝ) : ℂ)
  have hL : (((1 + E : ℝ)) : ℂ) = Complex.exp A := by
    rw [← Complex.ofReal_exp, Real.exp_log hL0]
  have hw : ((w e x : ℝ) : ℂ) = (E : ℂ) * Complex.exp (A * (-1 - e)) := by
    rw [w, Real.rpow_def_of_pos hL0, Complex.ofReal_mul, Complex.ofReal_exp]
    push_cast; rfl
  have key : Complex.exp (-A) * Complex.exp (-s * x) * Complex.exp (-(e : ℂ) * (A - s * x)) =
      Complex.exp (A * (-1 - e)) * Complex.exp (-s * (((1 - e) * x : ℝ) : ℂ)) := by
    rw [← Complex.exp_add, ← Complex.exp_add, ← Complex.exp_add]
    congr 1; push_cast; ring
  rw [hw, base, ell, Complex.ofReal_div, hL, div_eq_mul_inv, ← Complex.exp_neg]
  linear_combination (E : ℂ) * key

lemma w_nonneg (e x : ℝ) : 0 ≤ w e x := by
  unfold w; positivity

lemma w_le (e x : ℝ) (he : -1 ≤ e) : w e x ≤ Real.exp (-x) := by
  unfold w
  have h1 : (1 + Real.exp (-x)) ^ (-1 - e) ≤ 1 :=
    Real.rpow_le_one_of_one_le_of_nonpos (by linarith [Real.exp_pos (-x)]) (by linarith)
  have := Real.exp_pos (-x)
  nlinarith

lemma w_meas (e : ℝ) : Measurable (w e) := by
  unfold w; fun_prop

lemma w_int (e : ℝ) (he : -1 ≤ e) : IntegrableOn (w e) (Ioi 0) := by
  have hg := exp_neg_integrableOn_Ioi 0 (one_pos : (0 : ℝ) < 1)
  refine Integrable.mono' hg (w_meas e).aestronglyMeasurable ?_
  refine (ae_restrict_iff' measurableSet_Ioi).2 (Eventually.of_forall fun x hx => ?_)
  rw [Real.norm_eq_abs, abs_of_nonneg (w_nonneg e x), neg_one_mul]
  exact w_le e x he

lemma Phi_zero (e : ℝ) : Phi e 0 = (((∫ x in Ioi (0 : ℝ), w e x) : ℝ) : ℂ) := by
  rw [Phi, ← integral_complex_ofReal]
  refine setIntegral_congr_fun measurableSet_Ioi (fun x _ => ?_)
  rw [phi_integrand]; simp

lemma Lpt_cont : Continuous Lpt := by unfold Lpt; fun_prop

lemma HK_cont : ContinuousOn (fun r : ℝ => H (Lpt r) * K (Lpt r)) (Ioi 0) := by
  have hH : ContinuousOn (fun r : ℝ => H (Lpt r)) (Ioi 0) :=
    H_diff.continuousOn.comp Lpt_cont.continuousOn (fun r hr => Lpt_re_pos r (le_of_lt hr))
  have hK : ContinuousOn (fun r : ℝ => K (Lpt r)) (Ioi 0) := by
    intro r hr
    have hnot : ∀ n : ℤ, Lpt r ≠ 2 * n + 1 := by
      intro n h
      have := congrArg Complex.im h
      rw [Lpt_im] at this
      simp at this
      have : 0 < r := hr
      have : (0 : ℝ) < √2 := by positivity
      nlinarith
    exact ((Q05.proof_6.1 _ hnot).continuousAt.comp Lpt_cont.continuousAt).continuousWithinAt
  exact hH.mul hK

lemma HK_meas : AEStronglyMeasurable (fun r : ℝ => H (Lpt r) * K (Lpt r)) (volume.restrict (Ioi 0)) :=
  HK_cont.aestronglyMeasurable measurableSet_Ioi

lemma HK_le (r : ℝ) (hr : 0 ≤ r) : ‖H (Lpt r) * K (Lpt r)‖ ≤ 4 * Real.exp (-a * r) := by
  obtain ⟨hH, -, hK, -⟩ := stmt7_at r hr
  rw [norm_mul]
  calc ‖H (Lpt r)‖ * ‖K (Lpt r)‖ ≤ 2 * (2 * Real.exp (-a * r)) :=
        mul_le_mul hH hK (norm_nonneg _) (by norm_num)
    _ = 4 * Real.exp (-a * r) := by ring

lemma norm_ω : ‖ω‖ = 1 := by rw [ω, Complex.norm_exp]; simp

lemma Lpt_good (r : ℝ) (hr : 0 ≤ r) : ‖Lpt r‖ ≤ 2 * (Lpt r).re := by
  have h1 := norm_Lpt_le r hr
  rw [Lpt_re]
  have : (1 : ℝ) < √2 := by rw [Real.lt_sqrt (by norm_num)]; norm_num
  nlinarith

/-- The integrand of the double integral over `x ∈ (0,∞)` and `s ∈ Λ`. -/
noncomputable def G (e : ℝ) (x r : ℝ) : ℂ :=
  ((w e x : ℝ) : ℂ) * (H (Lpt r) * K (Lpt r) * Complex.exp (-(Lpt r) * (((1 - e) * x : ℝ) : ℂ)) * ω)

lemma norm_exp_le (e x r : ℝ) (he : e ≤ 1) (hx : 0 ≤ x) (hr : 0 ≤ r) :
    ‖Complex.exp (-(Lpt r) * (((1 - e) * x : ℝ) : ℂ))‖ ≤ 1 := by
  rw [Complex.norm_exp]
  apply Real.exp_le_one_iff.2
  have h1 := Lpt_re_pos r hr
  have h2 : 0 ≤ (1 - e) * x := mul_nonneg (by linarith) hx
  have : (-(Lpt r) * (((1 - e) * x : ℝ) : ℂ)).re = -((Lpt r).re * ((1 - e) * x)) := by
    rw [Complex.mul_re]; simp
  rw [this]
  nlinarith

lemma G_int (e : ℝ) (he : |e| ≤ 1 / 4) :
    Integrable (Function.uncurry (G e)) ((volume.restrict (Ioi 0)).prod (volume.restrict (Ioi 0))) := by
  have he' := abs_le.1 he
  have hg1 : Integrable (fun x : ℝ => Real.exp (-1 * x)) (volume.restrict (Ioi 0)) :=
    exp_neg_integrableOn_Ioi 0 one_pos
  have hg2 : Integrable (fun r : ℝ => 4 * Real.exp (-a * r)) (volume.restrict (Ioi 0)) :=
    (exp_neg_integrableOn_Ioi 0 a_pos).const_mul 4
  have hg := hg1.mul_prod hg2
  have hmeas : AEStronglyMeasurable (Function.uncurry (G e))
      ((volume.restrict (Ioi 0)).prod (volume.restrict (Ioi 0))) := by
    have hw : Measurable (fun p : ℝ × ℝ => ((w e p.1 : ℝ) : ℂ)) := by
      unfold w; fun_prop
    have hHK := HK_meas.comp_quasiMeasurePreserving
      (Measure.quasiMeasurePreserving_snd (μ := volume.restrict (Ioi (0 : ℝ)))
        (ν := volume.restrict (Ioi (0 : ℝ))))
    have hex : Measurable (fun p : ℝ × ℝ =>
        Complex.exp (-(Lpt p.2) * (((1 - e) * p.1 : ℝ) : ℂ))) := by
      unfold Lpt; fun_prop
    exact hw.aestronglyMeasurable.mul ((hHK.mul hex.aestronglyMeasurable).mul_const ω)
  refine Integrable.mono' hg hmeas ?_
  have hae : ∀ᵐ p ∂((volume.restrict (Ioi (0 : ℝ))).prod (volume.restrict (Ioi (0 : ℝ)))),
      p ∈ Ioi (0 : ℝ) ×ˢ Ioi (0 : ℝ) := by
    rw [Measure.prod_restrict]
    exact ae_restrict_mem (measurableSet_Ioi.prod measurableSet_Ioi)
  filter_upwards [hae] with p hp
  obtain ⟨hx, hr⟩ := hp
  have hx0 : 0 ≤ p.1 := le_of_lt hx
  have hr0 : 0 ≤ p.2 := le_of_lt hr
  show ‖G e p.1 p.2‖ ≤ _
  unfold G
  rw [norm_mul, norm_mul, norm_mul, Complex.norm_real, Real.norm_eq_abs,
    abs_of_nonneg (w_nonneg _ _), norm_ω, mul_one, neg_one_mul]
  have h1 := w_le e p.1 (by linarith)
  have h2 := HK_le p.2 hr0
  have h3 := norm_exp_le e p.1 p.2 (by linarith) hx0 hr0
  calc w e p.1 * (‖H (Lpt p.2) * K (Lpt p.2)‖ *
        ‖Complex.exp (-(Lpt p.2) * (((1 - e) * p.1 : ℝ) : ℂ))‖)
      ≤ Real.exp (-p.1) * ((4 * Real.exp (-a * p.2)) * 1) := by
        gcongr
    _ = Real.exp (-p.1) * (4 * Real.exp (-a * p.2)) := by ring

lemma G_inner (e r : ℝ) :
    ∫ x in Ioi (0 : ℝ), G e x r = H (Lpt r) * K (Lpt r) * Phi e (Lpt r) * ω := by
  have h : ∀ x, G e x r = (H (Lpt r) * K (Lpt r) * ω) *
      (base (Lpt r) x * Complex.exp (-(e : ℂ) * ell (Lpt r) x)) := by
    intro x; rw [phi_integrand, G]; ring
  simp_rw [h]
  rw [integral_const_mul, Phi]
  ring

lemma stepA (h1 : Q07.stmt_1) (e : ℝ) (he : |e| ≤ 1 / 4) :
    -(W02.FR e) / 2 = -(Phi e 0).re +
      2 / π * (∫ r in Ioi (0 : ℝ), H (Lpt r) * K (Lpt r) * Phi e (Lpt r) * ω).im := by
  have he' := abs_le.1 he
  have h8 := (Q05.proof_8 e (by linarith)).2
  have hpt : ∀ x ∈ Ioi (0 : ℝ), Real.exp (-x) * (1 + Real.exp (-x)) ^ (-1 - e) * Cf ((1 - e) * x) =
      w e x - 2 / π * (∫ r in Ioi (0 : ℝ), G e x r).im := by
    intro x hx
    have hx0 : (0 : ℝ) ≤ x := le_of_lt hx
    show w e x * Cf ((1 - e) * x) = _
    rw [(h1 ((1 - e) * x) (mul_nonneg (by linarith) hx0)).2, lineInt]
    simp only [G]
    rw [integral_const_mul, Complex.im_ofReal_mul]
    ring
  have hGi := G_int e he
  have hxint : Integrable (fun x => ∫ r in Ioi (0 : ℝ), G e x r) (volume.restrict (Ioi 0)) :=
    hGi.integral_prod_left
  have hxim : Integrable (fun x => 2 / π * (∫ r in Ioi (0 : ℝ), G e x r).im)
      (volume.restrict (Ioi 0)) := hxint.im.const_mul _
  have him := integral_im hxint
  simp only [RCLike.im_to_complex] at him
  rw [h8, setIntegral_congr_fun measurableSet_Ioi hpt, integral_sub (w_int e (by linarith)) hxim,
    integral_const_mul, him, integral_integral_swap hGi]
  simp_rw [G_inner]
  rw [Phi_zero, Complex.ofReal_re]
  ring

/-- The `n`-th term on `Λ`. -/
noncomputable def Fr (e : ℝ) (n : ℕ) (r : ℝ) : ℂ :=
  H (Lpt r) * K (Lpt r) * (∫ x in Ioi (0 : ℝ), fn e (Lpt r) n x) * ω

lemma Fr_eq (e : ℝ) (n : ℕ) (r : ℝ) :
    Fr e n r = (e : ℂ) ^ n * (H (Lpt r) * K (Lpt r) * P n (Lpt r) * ω) := by
  rw [Fr, int_fn]; ring

lemma Fr_meas (e : ℝ) (n : ℕ) : AEStronglyMeasurable (Fr e n) (volume.restrict (Ioi 0)) := by
  have hi : StronglyMeasurable (fun r : ℝ => ∫ x in Ioi (0 : ℝ), fn e (Lpt r) n x) := by
    apply StronglyMeasurable.integral_prod_right
    exact (fn_meas2 e n).stronglyMeasurable
  exact (HK_meas.mul hi.aestronglyMeasurable).mul_const ω

lemma int_fn_le (e : ℝ) (he : |e| ≤ 1 / 4) (s : ℂ) (hs : ‖s‖ ≤ 2 * s.re) (n : ℕ) :
    ‖∫ x in Ioi (0 : ℝ), fn e s n x‖ ≤ 2 * (1 / 2) ^ n := by
  have hg : Integrable (fun x : ℝ => 2 * (1 / 2 : ℝ) ^ n * Real.exp (-x)) (volume.restrict (Ioi 0)) := by
    have := (exp_neg_integrableOn_Ioi 0 (one_pos : (0 : ℝ) < 1)).const_mul (2 * (1 / 2 : ℝ) ^ n)
    exact IntegrableOn.congr_fun this (fun x _ => by simp) measurableSet_Ioi
  have h := norm_integral_le_of_norm_le hg ((ae_restrict_iff' measurableSet_Ioi).2
    (Eventually.of_forall fun x hx => fn_bnd e he s hs x (le_of_lt hx) n))
  rw [integral_const_mul, integral_exp_neg_Ioi_zero, mul_one] at h
  exact h

lemma Fr_le (e : ℝ) (he : |e| ≤ 1 / 4) (n : ℕ) (r : ℝ) (hr : 0 ≤ r) :
    ‖Fr e n r‖ ≤ 4 * Real.exp (-a * r) * (2 * (1 / 2) ^ n) := by
  rw [Fr, norm_mul, norm_mul, norm_ω, mul_one]
  exact mul_le_mul (HK_le r hr) (int_fn_le e he _ (Lpt_good r hr) n) (norm_nonneg _) (by positivity)

lemma Fr_int (e : ℝ) (he : |e| ≤ 1 / 4) (n : ℕ) : IntegrableOn (Fr e n) (Ioi 0) := by
  have hg : IntegrableOn (fun r : ℝ => 4 * Real.exp (-a * r) * (2 * (1 / 2 : ℝ) ^ n)) (Ioi 0) :=
    ((exp_neg_integrableOn_Ioi 0 a_pos).const_mul 4).mul_const _
  refine Integrable.mono' hg (Fr_meas e n) ?_
  exact (ae_restrict_iff' measurableSet_Ioi).2 (Eventually.of_forall fun r hr =>
    Fr_le e he n r (le_of_lt hr))

lemma stepC (e : ℝ) (he : |e| ≤ 1 / 4) :
    HasSum (fun n => (e : ℂ) ^ n * lineInt (fun s => H s * K s * P n s))
      (∫ r in Ioi (0 : ℝ), H (Lpt r) * K (Lpt r) * Phi e (Lpt r) * ω) := by
  have h := hasSum_integral_of_dominated_convergence (μ := volume.restrict (Ioi (0 : ℝ)))
    (F := Fr e) (f := fun r => H (Lpt r) * K (Lpt r) * Phi e (Lpt r) * ω)
    (fun n r => 4 * Real.exp (-a * r) * (2 * (1 / 2 : ℝ) ^ n)) (Fr_meas e)
    (fun n => (ae_restrict_iff' measurableSet_Ioi).2 (Eventually.of_forall fun r hr =>
      Fr_le e he n r (le_of_lt hr)))
    (Eventually.of_forall fun r => (summable_geometric_two.mul_left 2).mul_left _) ?_ ?_
  · convert h using 1
    funext n
    simp_rw [Fr_eq]
    rw [integral_const_mul, lineInt]
  · have ht : ∀ r : ℝ, ∑' n : ℕ, 4 * Real.exp (-a * r) * (2 * (1 / 2 : ℝ) ^ n) =
        16 * Real.exp (-a * r) := by
      intro r; rw [tsum_mul_left, tsum_mul_left, tsum_geometric_two]; ring
    simp_rw [ht]
    exact (exp_neg_integrableOn_Ioi 0 a_pos).const_mul 16
  · refine (ae_restrict_iff' measurableSet_Ioi).2 (Eventually.of_forall fun r hr => ?_)
    have := ((phi_hasSum e he (Lpt r) (Lpt_good r (le_of_lt hr))).mul_left
      (H (Lpt r) * K (Lpt r))).mul_right ω
    exact this

lemma HKP_int (n : ℕ) : IntegrableOn (fun r : ℝ => H (Lpt r) * K (Lpt r) * P n (Lpt r)) (Ioi 0) := by
  have h := (Fr_int (1 / 4) (by rw [abs_of_pos (by norm_num)]) n).const_mul
    ((4 : ℂ) ^ n * ω⁻¹)
  have hω : ω ≠ 0 := by rw [ω]; exact Complex.exp_ne_zero _
  refine IntegrableOn.congr_fun h (fun r _ => ?_) measurableSet_Ioi
  rw [Fr_eq]
  push_cast
  field_simp
  rw [← mul_pow]; norm_num

/-- The claimed value of `M_n`. -/
noncomputable def acoef (n : ℕ) : ℝ := -D n + 2 / π * (lineInt fun s => H s * K s * P n s).im

lemma real_series (h1 : Q07.stmt_1) (e : ℝ) (he : |e| ≤ 1 / 4) :
    HasSum (fun n => e ^ n * acoef n) (-(W02.FR e) / 2) := by
  have hD : HasSum (fun n => e ^ n * D n) (Phi e 0).re := by
    have := Complex.hasSum_re (phi_hasSum e he 0 (by simp))
    convert this using 1
    funext n
    rw [int_fn, P_zero, ← Complex.ofReal_pow, ← Complex.ofReal_mul, Complex.ofReal_re]
  have hI : HasSum (fun n => e ^ n * (lineInt fun s => H s * K s * P n s).im)
      (∫ r in Ioi (0 : ℝ), H (Lpt r) * K (Lpt r) * Phi e (Lpt r) * ω).im := by
    have := Complex.hasSum_im (stepC e he)
    convert this using 1
    funext n
    rw [← Complex.ofReal_pow, Complex.im_ofReal_mul]
  rw [stepA h1 e he]
  convert hD.neg.add (hI.mul_left (2 / π)) using 1
  funext n
  unfold acoef
  ring

lemma Mk_series (e : ℝ) (he : |e| < 1) :
    HasSum (fun n => W02.Mk n * (e : ℂ) ^ n) (-(W02.FC e) / 2) := by
  have hf : DifferentiableOn ℂ W02.FC (Metric.ball 0 1) := by
    refine W02.proof_5.2.mono (fun z hz => ?_)
    rw [mem_ball_zero_iff] at hz
    show z.re < 1
    exact lt_of_le_of_lt (Complex.re_le_norm z) hz
  have hz : (e : ℂ) ∈ Metric.ball (0 : ℂ) 1 := by
    rw [mem_ball_zero_iff, Complex.norm_real, Real.norm_eq_abs]; exact he
  have := (Complex.hasSum_taylorSeries_on_ball hf hz).mul_left (-1 / 2 : ℂ)
  convert this using 1
  · funext n
    rw [W02.Mk, sub_zero, smul_eq_mul, smul_eq_mul]
    have : (n.factorial : ℂ) ≠ 0 := by exact_mod_cast (Nat.factorial_pos n).ne'
    field_simp
  · ring

lemma diff_series (h1 : Q07.stmt_1) (e : ℝ) (he : |e| ≤ 1 / 4) :
    HasSum (fun n => (W02.Mk n - (acoef n : ℂ)) * (e : ℂ) ^ n) 0 := by
  have hM := Mk_series e (by linarith)
  have hR := (Complex.hasSum_ofReal).2 (real_series h1 e he)
  have h4 := (W02.proof_4 e (by linarith [(abs_le.1 he).2])).2
  rw [h4] at hM
  convert hM.sub hR using 1
  · funext n; push_cast; ring
  · push_cast; ring

lemma coeff_zero (b : ℕ → ℂ) (h : ∀ e : ℝ, |e| ≤ 1 / 4 → HasSum (fun n => b n * (e : ℂ) ^ n) 0) :
    ∀ n, b n = 0 := by
  let p : FormalMultilinearSeries ℝ ℝ ℂ := fun n => ContinuousMultilinearMap.mkPiRing ℝ (Fin n) (b n)
  have happ : ∀ (n : ℕ) (y : ℝ), p n (fun _ => y) = b n * (y : ℂ) ^ n := by
    intro n y
    simp only [p, ContinuousMultilinearMap.mkPiRing_apply, Finset.prod_const, Finset.card_univ,
      Fintype.card_fin, Complex.real_smul]
    push_cast; ring
  set r0 : NNReal := ⟨1 / 4, by norm_num⟩
  have hr0 : ((r0 : ℝ)) = 1 / 4 := rfl
  have hrad : (r0 : ENNReal) ≤ p.radius := by
    apply p.le_radius_of_tendsto (l := 0)
    have h0 := ((h (1 / 4) (by rw [abs_of_pos (by norm_num)])).summable.tendsto_atTop_zero).norm
    rw [norm_zero] at h0
    convert h0 using 1
    funext n
    rw [ContinuousMultilinearMap.norm_mkPiRing, norm_mul, norm_pow, Complex.norm_real, hr0,
      Real.norm_eq_abs, abs_of_pos (by norm_num)]
  have hps : HasFPowerSeriesOnBall (0 : ℝ → ℂ) p 0 r0 := by
    refine ⟨hrad, by rw [ENNReal.coe_pos, ← NNReal.coe_pos, hr0]; norm_num, fun {y} hy => ?_⟩
    rw [Metric.eball_coe, mem_ball_zero_iff, hr0, Real.norm_eq_abs] at hy
    simp only [happ, Pi.zero_apply]
    exact h y hy.le
  have hpz := hps.hasFPowerSeriesAt.eq_zero
  intro n
  have := congrArg (fun q : FormalMultilinearSeries ℝ ℝ ℂ => q n (fun _ => (1 : ℝ))) hpz
  simp only [happ] at this
  simpa using this

theorem main : Q07.stmt_1 → Q07.stmt_2 := by
  intro h1 k
  refine ⟨HKP_int k, ?_⟩
  have := coeff_zero (fun n => W02.Mk n - (acoef n : ℂ)) (fun e he => diff_series h1 e he) k
  exact sub_eq_zero.1 this

end Q07.PB

/-! ### Merged from work/Q07/PartC.lean -/
namespace Q07.PC

open Real MeasureTheory Set Filter

/-! ## Facts on the sector `|Im L| ≤ Re L` -/

lemma re_one_add_ge (L : ℂ) (hL : |L.im| ≤ L.re) (v : ℝ) (hv0 : 0 ≤ v) (hv1 : v ≤ 1) :
    1 / 2 ≤ (1 + (v : ℂ) * Complex.exp (-L)).re := by
  have hre : (1 + (v : ℂ) * Complex.exp (-L)).re =
      1 + v * (Real.exp (-L.re) * Real.cos L.im) := by
    simp [Complex.exp_re, Real.cos_neg]
  rw [hre]
  have he : 0 < Real.exp (-L.re) := Real.exp_pos _
  by_cases hc : 0 ≤ Real.cos L.im
  · have : 0 ≤ v * (Real.exp (-L.re) * Real.cos L.im) := by positivity
    linarith
  · push Not at hc
    have hw : π / 2 < |L.im| := by
      by_contra h
      push Not at h
      have := Real.cos_nonneg_of_mem_Icc ⟨by linarith [neg_abs_le L.im], le_trans (le_abs_self _) h⟩
      linarith
    have hu : 1 ≤ L.re := by linarith [Real.pi_gt_three]
    have hE : Real.exp (-L.re) ≤ 1 / 2 := by
      rw [Real.exp_neg, inv_eq_one_div, div_le_div_iff₀ (Real.exp_pos _) (by norm_num)]
      have := Real.add_one_le_exp L.re
      linarith
    have hc1 : -1 ≤ Real.cos L.im := Real.neg_one_le_cos _
    have h1 : -(Real.exp (-L.re)) ≤ Real.exp (-L.re) * Real.cos L.im := by nlinarith
    have h2 : -(Real.exp (-L.re)) ≤ v * (Real.exp (-L.re) * Real.cos L.im) := by
      have h3 : Real.exp (-L.re) * Real.cos L.im ≤ 0 := by nlinarith
      nlinarith
    linarith

lemma norm_one_add_ge (L : ℂ) (hL : |L.im| ≤ L.re) (v : ℝ) (hv0 : 0 ≤ v) (hv1 : v ≤ 1) :
    1 / 2 ≤ ‖1 + (v : ℂ) * Complex.exp (-L)‖ :=
  (re_one_add_ge L hL v hv0 hv1).trans (Complex.re_le_norm _)

lemma slit_of_cone (L : ℂ) (hL : |L.im| ≤ L.re) : 1 + Complex.exp (-L) ∈ Complex.slitPlane := by
  rw [Complex.mem_slitPlane_iff]
  left
  have := re_one_add_ge L hL 1 zero_le_one le_rfl
  simp only [Complex.ofReal_one, one_mul] at this
  linarith

lemma inv_one_add_le (L : ℂ) (hL : |L.im| ≤ L.re) : ‖(1 + Complex.exp (-L))⁻¹‖ ≤ 2 := by
  have h := norm_one_add_ge L hL 1 zero_le_one le_rfl
  simp only [Complex.ofReal_one, one_mul] at h
  rw [norm_inv]
  rw [inv_le_comm₀ (by linarith) (by norm_num)]
  linarith

/-- `ℓ(L) = Log(1 + e^{-L})`. -/
noncomputable def ell (L : ℂ) : ℂ := Complex.log (1 + Complex.exp (-L))

lemma ell_le (L : ℂ) (hL : |L.im| ≤ L.re) : ‖ell L‖ ≤ 2 * Real.exp (-L.re) := by
  rw [ell, Complex.log_eq_integral (slit_of_cone L hL), norm_mul]
  have hz : ‖Complex.exp (-L)‖ = Real.exp (-L.re) := by
    rw [Complex.norm_exp]; simp
  rw [hz]
  have hI : ‖∫ t in (0 : ℝ)..1, (1 + t • Complex.exp (-L))⁻¹‖ ≤ 2 * |(1 : ℝ) - 0| := by
    apply intervalIntegral.norm_integral_le_of_norm_le_const
    intro t ht
    rw [uIoc_of_le zero_le_one] at ht
    have h := norm_one_add_ge L hL t (le_of_lt ht.1) ht.2
    rw [Complex.real_smul, norm_inv, inv_le_comm₀ (by linarith) (by norm_num)]
    linarith
  rw [sub_zero, abs_one, mul_one] at hI
  have := Real.exp_pos (-L.re)
  nlinarith [norm_nonneg (∫ t in (0 : ℝ)..1, (1 + t • Complex.exp (-L))⁻¹)]

/-! ## The domain of holomorphy and the triangle -/

/-- A set on which `1 + e^{-L}` lies in the slit plane. -/
def Dom : Set ℂ := {L | 0 < L.re ∨ |L.im| < π}

lemma dom_slit (L : ℂ) (hL : L ∈ Dom) : 1 + Complex.exp (-L) ∈ Complex.slitPlane := by
  rw [Complex.mem_slitPlane_iff]
  have hre : (1 + Complex.exp (-L)).re = 1 + Real.exp (-L.re) * Real.cos L.im := by
    simp [Complex.exp_re, Real.cos_neg]
  have him : (1 + Complex.exp (-L)).im = -(Real.exp (-L.re) * Real.sin L.im) := by
    simp [Complex.exp_im, Real.sin_neg]
  have he : 0 < Real.exp (-L.re) := Real.exp_pos _
  rcases hL with h | h
  · left
    rw [hre]
    have hE : Real.exp (-L.re) < 1 := Real.exp_lt_one_iff.2 (by linarith)
    have hc1 : -1 ≤ Real.cos L.im := Real.neg_one_le_cos _
    nlinarith
  · by_cases hs : Real.sin L.im = 0
    · left
      rw [hre]
      have h0 : L.im = 0 := (Real.sin_eq_zero_iff_of_lt_of_lt (by linarith [neg_abs_le L.im])
        (lt_of_le_of_lt (le_abs_self _) h)).1 hs
      rw [h0, Real.cos_zero]
      linarith
    · right
      rw [him]
      intro h'
      apply hs
      have : Real.exp (-L.re) * Real.sin L.im = 0 := by linarith
      rcases mul_eq_zero.1 this with h1 | h1
      · linarith
      · exact h1

lemma ball_sub_dom (R : ℝ) (hR : 1 ≤ R) : Metric.ball (R : ℂ) (R + 1 / R) ⊆ Dom := by
  intro L hL
  rw [Metric.mem_ball, dist_eq_norm, Complex.norm_eq_sqrt_sq_add_sq] at hL
  have hR0 : 0 < R := by linarith
  by_cases h : 0 < L.re
  · exact Or.inl h
  · right
    push Not at h
    have h1 : (L - R).re ^ 2 + (L - R).im ^ 2 < (R + 1 / R) ^ 2 := by
      rwa [Real.sqrt_lt' (by positivity)] at hL
    simp only [Complex.sub_re, Complex.ofReal_re, Complex.sub_im, Complex.ofReal_im,
      sub_zero] at h1
    have hRR : R * (1 / R) = 1 := by field_simp
    have hR1 : 1 / R ≤ 1 := by rw [div_le_one hR0]; exact hR
    have hR2 : 0 < 1 / R := by positivity
    have h2 : L.im ^ 2 < π ^ 2 := by nlinarith [Real.pi_gt_three]
    have := sq_lt_sq.1 h2
    rwa [abs_of_pos Real.pi_pos] at this

lemma cone_mem_ball (R : ℝ) (hR : 1 ≤ R) (L : ℂ) (hc : |L.im| ≤ L.re) (hn : ‖L‖ ≤ R) :
    L ∈ Metric.ball (R : ℂ) (R + 1 / R) := by
  have hR0 : 0 < R := by linarith
  rw [Metric.mem_ball, dist_eq_norm]
  have hre : L.re ≤ R := (Complex.re_le_norm L).trans hn
  have hre0 : 0 ≤ L.re := (abs_nonneg _).trans hc
  have hsq : L.im ^ 2 ≤ L.re ^ 2 := by
    have := sq_le_sq' (neg_le_of_abs_le hc) (le_of_abs_le hc)
    exact this
  have h1 : ‖L - R‖ ≤ R := by
    rw [Complex.norm_eq_sqrt_sq_add_sq]
    calc √((L - R).re ^ 2 + (L - R).im ^ 2) ≤ √(R ^ 2) := by
          apply Real.sqrt_le_sqrt
          simp only [Complex.sub_re, Complex.ofReal_re, Complex.sub_im, Complex.ofReal_im,
            sub_zero]
          nlinarith
      _ = R := Real.sqrt_sq hR0.le
  have : 0 < 1 / R := by positivity
  linarith

/-! ## The Gamma integral -/

lemma gamma_int (n : ℕ) (b : ℝ) (hb : 0 < b) :
    IntegrableOn (fun l : ℝ => l ^ n * Real.exp (-(b * l))) (Ioi 0) ∧
      ∫ l in Ioi (0 : ℝ), l ^ n * Real.exp (-(b * l)) = n.factorial / b ^ (n + 1) := by
  constructor
  · have h := integrableOn_rpow_mul_exp_neg_mul_rpow (s := (n : ℝ)) (p := 1) (b := b)
      (by have := Nat.cast_nonneg (α := ℝ) n; linarith) one_pos hb
    refine h.congr_fun (fun x _ => ?_) measurableSet_Ioi
    simp only [Real.rpow_natCast, Real.rpow_one, neg_mul]
  · have h := Real.integral_rpow_mul_exp_neg_mul_Ioi (a := (n : ℝ) + 1) (r := b) (by positivity) hb
    have h1 : ∫ l in Ioi (0 : ℝ), l ^ n * Real.exp (-(b * l)) =
        ∫ t in Ioi (0 : ℝ), t ^ ((n : ℝ) + 1 - 1) * Real.exp (-(b * t)) := by
      refine setIntegral_congr_fun measurableSet_Ioi (fun x _ => ?_)
      simp only [add_sub_cancel_right, Real.rpow_natCast]
    rw [h1, h, Real.Gamma_nat_eq_factorial]
    rw [show ((n : ℝ) + 1) = ((n + 1 : ℕ) : ℝ) by push_cast; ring, Real.rpow_natCast, div_pow,
      one_pow]
    ring

/-! ## Rotation of the half line `(0,∞)` to the ray through `conj σ / |σ|` -/

/-- The direction `conj σ / |σ|`. -/
noncomputable def dir (σ : ℂ) : ℂ := (starRingEnd ℂ σ) / (‖σ‖ : ℂ)

lemma dir_re (σ : ℂ) : (dir σ).re = σ.re / ‖σ‖ := by
  simp [dir, Complex.div_ofReal_re]

lemma dir_im (σ : ℂ) : (dir σ).im = -σ.im / ‖σ‖ := by
  simp [dir, Complex.div_ofReal_im, neg_div]

lemma norm_dir (σ : ℂ) (hσ : σ ≠ 0) : ‖dir σ‖ = 1 := by
  rw [dir, norm_div, Complex.norm_conj, Complex.norm_real, Real.norm_eq_abs, abs_norm,
    div_self (norm_ne_zero_iff.2 hσ)]

lemma mul_dir (σ : ℂ) (hσ : σ ≠ 0) : σ * dir σ = (‖σ‖ : ℂ) := by
  have h0 : (‖σ‖ : ℂ) ≠ 0 := by exact_mod_cast norm_ne_zero_iff.2 hσ
  rw [dir, mul_div_assoc', Complex.mul_conj, Complex.normSq_eq_norm_sq]
  push_cast
  field_simp

lemma dir_cone (σ : ℂ) (hσ : 0 < σ.re) (hσi : |σ.im| ≤ σ.re) : |(dir σ).im| ≤ (dir σ).re := by
  have hn : 0 < ‖σ‖ := lt_of_lt_of_le hσ (Complex.re_le_norm σ)
  rw [dir_re, dir_im, abs_div, abs_neg, abs_of_pos hn]
  exact div_le_div_of_nonneg_right hσi hn.le

lemma comb_cone (z : ℂ) (hzc : |z.im| ≤ z.re) (a b : ℝ) (ha : 0 ≤ a) (hb : 0 ≤ b) :
    |((a : ℂ) + (b : ℂ) * z).im| ≤ ((a : ℂ) + (b : ℂ) * z).re := by
  simp only [Complex.add_im, Complex.ofReal_im, zero_add, Complex.add_re, Complex.ofReal_re,
    Complex.re_ofReal_mul, Complex.im_ofReal_mul]
  rw [abs_mul, abs_of_nonneg hb]
  have := mul_le_mul_of_nonneg_left hzc hb
  linarith

lemma comb_norm (z : ℂ) (hz : ‖z‖ = 1) (a b : ℝ) (ha : 0 ≤ a) (hb : 0 ≤ b) :
    ‖(a : ℂ) + (b : ℂ) * z‖ ≤ a + b := by
  refine (norm_add_le _ _).trans ?_
  rw [norm_mul, hz, mul_one, Complex.norm_real, Complex.norm_real, Real.norm_eq_abs,
    Real.norm_eq_abs, abs_of_nonneg ha, abs_of_nonneg hb]

lemma rot_FTC (f : ℂ → ℂ) (z : ℂ) (hz : ‖z‖ = 1) (hzc : |z.im| ≤ z.re)
    (hd : ∀ L ∈ Dom, DifferentiableAt ℂ f L) (R : ℝ) (hR : 1 ≤ R) :
    ∫ x in (0 : ℝ)..R, f x =
      (∫ l in (0 : ℝ)..R, f (l * z) * z) -
        ∫ t in (0 : ℝ)..1, f (R + t * (R * z - R)) * (R * z - R) := by
  have hR0 : 0 ≤ R := by linarith
  set B := Metric.ball (R : ℂ) (R + 1 / R) with hB
  have hdiff : DifferentiableOn ℂ f B :=
    fun L hL => (hd L (ball_sub_dom R hR hL)).differentiableWithinAt
  obtain ⟨F, hF⟩ := hdiff.isExactOn_ball
  have hcont : ContinuousOn f B := hdiff.continuousOn
  have hmem : ∀ a b : ℝ, 0 ≤ a → 0 ≤ b → a + b ≤ R → ((a : ℂ) + (b : ℂ) * z) ∈ B :=
    fun a b ha hb hab => cone_mem_ball R hR _ (comb_cone z hzc a b ha hb)
      ((comb_norm z hz a b ha hb).trans hab)
  have hmR : ∀ x ∈ uIcc (0 : ℝ) R, (x : ℂ) ∈ B := by
    intro x hx
    rw [uIcc_of_le hR0] at hx
    have := hmem x 0 hx.1 le_rfl (by linarith [hx.2])
    simpa using this
  have hmL : ∀ l ∈ uIcc (0 : ℝ) R, ((l : ℂ) * z) ∈ B := by
    intro l hl
    rw [uIcc_of_le hR0] at hl
    have := hmem 0 l le_rfl hl.1 (by linarith [hl.2])
    simpa using this
  have hmC : ∀ t ∈ uIcc (0 : ℝ) 1, ((R : ℂ) + (t : ℂ) * (R * z - R)) ∈ B := by
    intro t ht
    rw [uIcc_of_le zero_le_one] at ht
    have h1 : 0 ≤ (1 - t) * R := mul_nonneg (by linarith [ht.2]) hR0
    have h2 : 0 ≤ t * R := mul_nonneg ht.1 hR0
    have := hmem ((1 - t) * R) (t * R) h1 h2 (by nlinarith)
    convert this using 1
    push_cast; ring
  have hA : ∫ x in (0 : ℝ)..R, f x = F R - F 0 := by
    have := intervalIntegral.integral_eq_sub_of_hasDerivAt (f := fun x : ℝ => F x)
      (f' := fun x : ℝ => f x) (a := 0) (b := R)
      (fun x hx => (hF _ (hmR x hx)).comp_ofReal)
      (ContinuousOn.intervalIntegrable
        (hcont.comp Complex.continuous_ofReal.continuousOn (fun x hx => hmR x hx)))
    simpa using this
  have hL : ∫ l in (0 : ℝ)..R, f (l * z) * z = F (R * z) - F 0 := by
    have := intervalIntegral.integral_eq_sub_of_hasDerivAt (f := fun l : ℝ => F (l * z))
      (f' := fun l : ℝ => f (l * z) * z) (a := 0) (b := R)
      (fun l hl => by
        have h1 : HasDerivAt (fun w : ℂ => w * z) z (l : ℂ) := by
          simpa using (hasDerivAt_id (l : ℂ)).mul_const z
        exact ((hF _ (hmL l hl)).comp (l : ℂ) h1).comp_ofReal)
      (ContinuousOn.intervalIntegrable (ContinuousOn.mul
        (hcont.comp (by fun_prop) (fun l hl => hmL l hl)) continuousOn_const))
    simpa using this
  have hC : ∫ t in (0 : ℝ)..1, f (R + t * (R * z - R)) * (R * z - R) = F (R * z) - F R := by
    have := intervalIntegral.integral_eq_sub_of_hasDerivAt
      (f := fun t : ℝ => F (R + t * (R * z - R)))
      (f' := fun t : ℝ => f (R + t * (R * z - R)) * (R * z - R)) (a := 0) (b := 1)
      (fun t ht => by
        have h1 : HasDerivAt (fun w : ℂ => (R : ℂ) + w * (R * z - R)) (R * z - R) (t : ℂ) := by
          simpa using ((hasDerivAt_id (t : ℂ)).mul_const ((R : ℂ) * z - R)).const_add (R : ℂ)
        exact ((hF _ (hmC t ht)).comp (t : ℂ) h1).comp_ofReal)
      (ContinuousOn.intervalIntegrable (ContinuousOn.mul
        (hcont.comp (by fun_prop) (fun t ht => hmC t ht)) continuousOn_const))
    simpa using this
  rw [hA, hL, hC]
  ring

lemma chord_le (f : ℂ → ℂ) (σ : ℂ) (n : ℕ) (c : ℝ) (hc : 0 ≤ c) (hσ : 0 < σ.re)
    (hσi : |σ.im| ≤ σ.re)
    (hb : ∀ L : ℂ, |L.im| ≤ L.re → ‖f L‖ ≤ c * ‖L‖ ^ n * Real.exp (-(σ * L).re))
    (R : ℝ) (hR : 1 ≤ R) :
    ‖∫ t in (0 : ℝ)..1, f (R + t * (R * dir σ - R)) * (R * dir σ - R)‖ ≤
      c * R ^ n * Real.exp (-(σ.re * R)) * (2 * R) := by
  have hR0 : 0 ≤ R := by linarith
  have hσ0 : σ ≠ 0 := fun h => by rw [h] at hσ; simp at hσ
  have hn : σ.re ≤ ‖σ‖ := Complex.re_le_norm σ
  have hz := norm_dir σ hσ0
  have hzc := dir_cone σ hσ hσi
  have h := intervalIntegral.norm_integral_le_of_norm_le_const
    (f := fun t : ℝ => f (R + t * (R * dir σ - R)) * (R * dir σ - R)) (a := 0) (b := 1)
    (C := c * R ^ n * Real.exp (-(σ.re * R)) * (2 * R)) (by
      intro t ht
      rw [uIoc_of_le zero_le_one] at ht
      have ht0 : 0 ≤ t := le_of_lt ht.1
      have h1 : 0 ≤ (1 - t) * R := mul_nonneg (by linarith [ht.2]) hR0
      have h2 : 0 ≤ t * R := mul_nonneg ht0 hR0
      have hpt : ((R : ℂ) + (t : ℂ) * (R * dir σ - R)) =
          (((1 - t) * R : ℝ) : ℂ) + ((t * R : ℝ) : ℂ) * dir σ := by push_cast; ring
      have hcone := comb_cone (dir σ) hzc _ _ h1 h2
      have hnorm := comb_norm (dir σ) hz _ _ h1 h2
      have hre : (σ * ((((1 - t) * R : ℝ) : ℂ) + ((t * R : ℝ) : ℂ) * dir σ)).re =
          (1 - t) * R * σ.re + t * R * ‖σ‖ := by
        rw [mul_add, mul_left_comm, mul_dir σ hσ0]
        simp [Complex.mul_re]
        ring
      have hfb := hb _ hcone
      rw [hre] at hfb
      have hN : ‖(((1 - t) * R : ℝ) : ℂ) + ((t * R : ℝ) : ℂ) * dir σ‖ ≤ R := by
        refine hnorm.trans (le_of_eq (by ring))
      have hE : Real.exp (-((1 - t) * R * σ.re + t * R * ‖σ‖)) ≤ Real.exp (-(σ.re * R)) := by
        apply Real.exp_le_exp.2
        have : t * R * σ.re ≤ t * R * ‖σ‖ := mul_le_mul_of_nonneg_left hn h2
        nlinarith
      have hP : ‖(((1 - t) * R : ℝ) : ℂ) + ((t * R : ℝ) : ℂ) * dir σ‖ ^ n ≤ R ^ n :=
        pow_le_pow_left₀ (norm_nonneg _) hN n
      have hD : ‖(R : ℂ) * dir σ - R‖ ≤ 2 * R := by
        refine (norm_sub_le _ _).trans ?_
        rw [norm_mul, hz, mul_one, Complex.norm_real, Real.norm_eq_abs, abs_of_nonneg hR0]
        linarith
      rw [norm_mul, hpt]
      have hf0 : ‖f ((((1 - t) * R : ℝ) : ℂ) + ((t * R : ℝ) : ℂ) * dir σ)‖ ≤
          c * R ^ n * Real.exp (-(σ.re * R)) := by
        refine hfb.trans ?_
        have := mul_le_mul hP hE (Real.exp_pos _).le (by positivity)
        calc c * ‖(((1 - t) * R : ℝ) : ℂ) + ((t * R : ℝ) : ℂ) * dir σ‖ ^ n *
              Real.exp (-((1 - t) * R * σ.re + t * R * ‖σ‖))
            = c * (‖(((1 - t) * R : ℝ) : ℂ) + ((t * R : ℝ) : ℂ) * dir σ‖ ^ n *
              Real.exp (-((1 - t) * R * σ.re + t * R * ‖σ‖))) := by ring
          _ ≤ c * (R ^ n * Real.exp (-(σ.re * R))) := mul_le_mul_of_nonneg_left this hc
          _ = c * R ^ n * Real.exp (-(σ.re * R)) := by ring
      exact mul_le_mul hf0 hD (norm_nonneg _) (by positivity))
  simpa using h

lemma ray_le (f : ℂ → ℂ) (σ : ℂ) (n : ℕ) (c : ℝ) (hc : 0 ≤ c) (hσ : 0 < σ.re)
    (hσi : |σ.im| ≤ σ.re)
    (hb : ∀ L : ℂ, |L.im| ≤ L.re → ‖f L‖ ≤ c * ‖L‖ ^ n * Real.exp (-(σ * L).re))
    (R : ℝ) (hR : 0 ≤ R) :
    ‖∫ l in (0 : ℝ)..R, f (l * dir σ) * dir σ‖ ≤ c * n.factorial / ‖σ‖ ^ (n + 1) := by
  have hσ0 : σ ≠ 0 := fun h => by rw [h] at hσ; simp at hσ
  have hnp : 0 < ‖σ‖ := lt_of_lt_of_le hσ (Complex.re_le_norm σ)
  have hz := norm_dir σ hσ0
  have hzc := dir_cone σ hσ hσi
  obtain ⟨hgi, hgv⟩ := gamma_int n ‖σ‖ hnp
  have hgc : Continuous (fun l : ℝ => c * (l ^ n * Real.exp (-(‖σ‖ * l)))) := by fun_prop
  have h1 := intervalIntegral.norm_integral_le_of_norm_le (a := 0) (b := R)
    (f := fun l : ℝ => f (l * dir σ) * dir σ)
    (g := fun l : ℝ => c * (l ^ n * Real.exp (-(‖σ‖ * l)))) hR
    (Filter.Eventually.of_forall fun l hl => by
      have hl0 : 0 ≤ l := le_of_lt hl.1
      have hcone := comb_cone (dir σ) hzc 0 l le_rfl hl0
      simp only [Complex.ofReal_zero, zero_add] at hcone
      have hfb := hb _ hcone
      have hre : (σ * ((l : ℂ) * dir σ)).re = ‖σ‖ * l := by
        rw [mul_left_comm, mul_dir σ hσ0]; simp [Complex.mul_re]; ring
      have hN : ‖(l : ℂ) * dir σ‖ = l := by
        rw [norm_mul, hz, mul_one, Complex.norm_real, Real.norm_eq_abs, abs_of_nonneg hl0]
      rw [hre, hN] at hfb
      rw [norm_mul, hz, mul_one]
      calc ‖f (l * dir σ)‖ ≤ c * l ^ n * Real.exp (-(‖σ‖ * l)) := hfb
        _ = c * (l ^ n * Real.exp (-(‖σ‖ * l))) := by ring)
    (hgc.intervalIntegrable (μ := volume) 0 R)
  refine h1.trans ?_
  rw [intervalIntegral.integral_of_le hR]
  have h2 : ∫ l in Ioc 0 R, c * (l ^ n * Real.exp (-(‖σ‖ * l))) ≤
      ∫ l in Ioi 0, c * (l ^ n * Real.exp (-(‖σ‖ * l))) := by
    apply setIntegral_mono_set (hgi.const_mul c)
    · refine (ae_restrict_iff' measurableSet_Ioi).2 (Filter.Eventually.of_forall fun l hl => ?_)
      have : 0 ≤ l := le_of_lt hl
      positivity
    · exact Filter.Eventually.of_forall Ioc_subset_Ioi_self
  refine h2.trans (le_of_eq ?_)
  rw [integral_const_mul, hgv]
  ring

lemma real_int (f : ℂ → ℂ) (σ : ℂ) (n : ℕ) (c : ℝ) (hσ : 0 < σ.re)
    (hd : ∀ L ∈ Dom, DifferentiableAt ℂ f L)
    (hb : ∀ L : ℂ, |L.im| ≤ L.re → ‖f L‖ ≤ c * ‖L‖ ^ n * Real.exp (-(σ * L).re)) :
    IntegrableOn (fun x : ℝ => f x) (Ioi 0) := by
  obtain ⟨hgi, -⟩ := gamma_int n σ.re hσ
  have hcont : ContinuousOn (fun x : ℝ => f x) (Ioi 0) := by
    intro x hx
    have hx0 : (0 : ℝ) < x := hx
    have : (x : ℂ) ∈ Dom := Or.inl (by simpa using hx0)
    exact ((hd _ this).continuousAt.comp Complex.continuous_ofReal.continuousAt).continuousWithinAt
  refine Integrable.mono' (hgi.const_mul c) (hcont.aestronglyMeasurable measurableSet_Ioi) ?_
  refine (ae_restrict_iff' measurableSet_Ioi).2 (Filter.Eventually.of_forall fun x hx => ?_)
  have hx0 : (0 : ℝ) ≤ x := le_of_lt hx
  have hcone : |(x : ℂ).im| ≤ (x : ℂ).re := by simpa using hx0
  have h := hb _ hcone
  have hre : (σ * (x : ℂ)).re = σ.re * x := by simp [Complex.mul_re]
  rw [hre, Complex.norm_real, Real.norm_eq_abs, abs_of_nonneg hx0] at h
  calc ‖f x‖ ≤ c * x ^ n * Real.exp (-(σ.re * x)) := h
    _ = c * (x ^ n * Real.exp (-(σ.re * x))) := by ring

lemma chord_tendsto (f : ℂ → ℂ) (σ : ℂ) (n : ℕ) (c : ℝ) (hc : 0 ≤ c) (hσ : 0 < σ.re)
    (hσi : |σ.im| ≤ σ.re)
    (hb : ∀ L : ℂ, |L.im| ≤ L.re → ‖f L‖ ≤ c * ‖L‖ ^ n * Real.exp (-(σ * L).re)) :
    Tendsto (fun R : ℝ => ∫ t in (0 : ℝ)..1, f (R + t * (R * dir σ - R)) * (R * dir σ - R))
      atTop (nhds 0) := by
  have hg : Tendsto (fun R : ℝ => c * R ^ n * Real.exp (-(σ.re * R)) * (2 * R)) atTop (nhds 0) := by
    have h := (tendsto_rpow_mul_exp_neg_mul_atTop_nhds_zero ((n + 1 : ℕ) : ℝ) σ.re hσ).const_mul
      (2 * c)
    rw [mul_zero] at h
    refine h.congr' ?_
    filter_upwards with R
    rw [Real.rpow_natCast, pow_succ]
    ring_nf
  refine squeeze_zero_norm' ?_ hg
  filter_upwards [eventually_ge_atTop (1 : ℝ)] with R hR
  exact chord_le f σ n c hc hσ hσi hb R hR

lemma rot_tendsto (f : ℂ → ℂ) (σ : ℂ) (n : ℕ) (c : ℝ) (hc : 0 ≤ c) (hσ : 0 < σ.re)
    (hσi : |σ.im| ≤ σ.re) (hd : ∀ L ∈ Dom, DifferentiableAt ℂ f L)
    (hb : ∀ L : ℂ, |L.im| ≤ L.re → ‖f L‖ ≤ c * ‖L‖ ^ n * Real.exp (-(σ * L).re)) :
    IntegrableOn (fun x : ℝ => f x) (Ioi 0) ∧
      Tendsto (fun R : ℝ => ∫ l in (0 : ℝ)..R, f (l * dir σ) * dir σ) atTop
        (nhds (∫ x in Ioi (0 : ℝ), f x)) := by
  have hσ0 : σ ≠ 0 := fun h => by rw [h] at hσ; simp at hσ
  have hI := real_int f σ n c hσ hd hb
  refine ⟨hI, ?_⟩
  have tI := intervalIntegral_tendsto_integral_Ioi (0 : ℝ) hI tendsto_id
  have tC := chord_tendsto f σ n c hc hσ hσi hb
  have h := tI.add tC
  rw [add_zero] at h
  refine h.congr' ?_
  filter_upwards [eventually_ge_atTop (1 : ℝ)] with R hR
  simp only [id]
  rw [rot_FTC f (dir σ) (norm_dir σ hσ0) (dir_cone σ hσ hσi) hd R hR]
  ring

lemma rot_bound (f : ℂ → ℂ) (σ : ℂ) (n : ℕ) (c : ℝ) (hc : 0 ≤ c) (hσ : 0 < σ.re)
    (hσi : |σ.im| ≤ σ.re) (hd : ∀ L ∈ Dom, DifferentiableAt ℂ f L)
    (hb : ∀ L : ℂ, |L.im| ≤ L.re → ‖f L‖ ≤ c * ‖L‖ ^ n * Real.exp (-(σ * L).re)) :
    IntegrableOn (fun x : ℝ => f x) (Ioi 0) ∧
      ‖∫ x in Ioi (0 : ℝ), f x‖ ≤ c * n.factorial / ‖σ‖ ^ (n + 1) := by
  obtain ⟨hI, hT⟩ := rot_tendsto f σ n c hc hσ hσi hd hb
  refine ⟨hI, le_of_tendsto hT.norm ?_⟩
  filter_upwards [eventually_ge_atTop (0 : ℝ)] with R hR
  exact ray_le f σ n c hc hσ hσi hb R hR

/-! ## `∫_0^∞ x^k e^{-σx} dx = k!/σ^{k+1}` for complex `σ` in the sector -/

lemma cgamma (k : ℕ) (σ : ℂ) (hσ : 0 < σ.re) (hσi : |σ.im| ≤ σ.re) :
    IntegrableOn (fun x : ℝ => (x : ℂ) ^ k * Complex.exp (-(σ * x))) (Ioi 0) ∧
      ∫ x in Ioi (0 : ℝ), (x : ℂ) ^ k * Complex.exp (-(σ * x)) =
        (k.factorial : ℂ) / σ ^ (k + 1) := by
  have hσ0 : σ ≠ 0 := fun h => by rw [h] at hσ; simp at hσ
  have hnp : 0 < ‖σ‖ := lt_of_lt_of_le hσ (Complex.re_le_norm σ)
  set f : ℂ → ℂ := fun L => L ^ k * Complex.exp (-(σ * L)) with hf
  have hd : ∀ L ∈ Dom, DifferentiableAt ℂ f L := fun L _ => by rw [hf]; fun_prop
  have hb : ∀ L : ℂ, |L.im| ≤ L.re → ‖f L‖ ≤ 1 * ‖L‖ ^ k * Real.exp (-(σ * L).re) := by
    intro L _
    rw [hf]; simp only
    rw [norm_mul, norm_pow, Complex.norm_exp, Complex.neg_re, one_mul]
  obtain ⟨hI, hT⟩ := rot_tendsto f σ k 1 zero_le_one hσ hσi hd hb
  refine ⟨hI, ?_⟩
  set z := dir σ with hz
  have hray : ∀ R : ℝ, ∫ l in (0 : ℝ)..R, f (l * z) * z =
      z ^ (k + 1) * ((∫ l in (0 : ℝ)..R, l ^ k * Real.exp (-(‖σ‖ * l)) : ℝ) : ℂ) := by
    intro R
    rw [← intervalIntegral.integral_ofReal, ← intervalIntegral.integral_const_mul]
    refine intervalIntegral.integral_congr (fun l _ => ?_)
    simp only [hf]
    have h1 : σ * ((l : ℂ) * z) = (l : ℂ) * (‖σ‖ : ℂ) := by
      rw [mul_left_comm, hz, mul_dir σ hσ0]
    rw [h1]
    push_cast
    ring_nf
  obtain ⟨hgi, hgv⟩ := gamma_int k ‖σ‖ hnp
  have tR := intervalIntegral_tendsto_integral_Ioi (0 : ℝ) hgi tendsto_id
  rw [hgv] at tR
  have tR2 : Tendsto (fun R : ℝ => z ^ (k + 1) *
      ((∫ l in (0 : ℝ)..R, l ^ k * Real.exp (-(‖σ‖ * l)) : ℝ) : ℂ)) atTop
      (nhds (z ^ (k + 1) * (((k.factorial / ‖σ‖ ^ (k + 1) : ℝ)) : ℂ))) :=
    ((Complex.continuous_ofReal.tendsto _).comp tR).const_mul _
  have hT' : Tendsto (fun R : ℝ => z ^ (k + 1) *
      ((∫ l in (0 : ℝ)..R, l ^ k * Real.exp (-(‖σ‖ * l)) : ℝ) : ℂ)) atTop
      (nhds (∫ x in Ioi (0 : ℝ), f x)) := hT.congr (fun R => hray R)
  have heq := tendsto_nhds_unique hT' tR2
  simp only [hf] at heq
  rw [heq]
  have hzs : z = (‖σ‖ : ℂ) / σ := by
    rw [eq_div_iff hσ0, mul_comm, hz, mul_dir σ hσ0]
  have h0 : (‖σ‖ : ℂ) ≠ 0 := by exact_mod_cast hnp.ne'
  rw [hzs, div_pow]
  push_cast
  field_simp

/-! ## `P_k` as an integral over `(0,∞)` (substitution `t = e^{-L}`) -/

/-- The integrand of `P_k` after `t = e^{-L}`. -/
noncomputable def Pint (k : ℕ) (s : ℂ) (L : ℂ) : ℂ :=
  Complex.exp (-((s + 1) * L)) * (s * L - ell L) ^ k / (1 + Complex.exp (-L))

lemma image_exp_neg : (fun u : ℝ => Real.exp (-u)) '' Ioi 0 = Ioo 0 1 := by
  ext x
  simp only [mem_image, mem_Ioi, mem_Ioo]
  constructor
  · rintro ⟨u, hu, rfl⟩
    exact ⟨Real.exp_pos _, Real.exp_lt_one_iff.2 (by linarith)⟩
  · rintro ⟨h0, h1⟩
    refine ⟨-Real.log x, ?_, ?_⟩
    · have := Real.log_neg h0 h1; linarith
    · rw [neg_neg, Real.exp_log h0]

lemma P_integrand (k : ℕ) (s : ℂ) (x : ℝ) :
    |(-Real.exp (-x))| • (((Real.exp (-x) : ℝ) : ℂ) ^ s / (1 + ((Real.exp (-x) : ℝ) : ℂ)) *
      (s * (Real.log (Real.exp (-x)) : ℂ) + (Real.log (1 + Real.exp (-x)) : ℂ)) ^ k) =
    (-1) ^ k * Pint k s x := by
  have he : 0 < Real.exp (-x) := Real.exp_pos _
  rw [abs_neg, abs_of_pos he, Complex.real_smul, Real.log_exp]
  have h1 : ((Real.exp (-x) : ℝ) : ℂ) ^ s = Complex.exp (-x * s) := by
    rw [Complex.cpow_def_of_ne_zero (by exact_mod_cast he.ne'), ← Complex.ofReal_log he.le,
      Real.log_exp]
    push_cast; ring_nf
  have h2 : ((Real.log (1 + Real.exp (-x)) : ℝ) : ℂ) = ell x := by
    rw [Complex.ofReal_log (by positivity), ell]
    push_cast; rfl
  rw [h1, h2, Pint]
  have h3 : Complex.exp (-((s + 1) * x)) = Complex.exp (-x) * Complex.exp (-x * s) := by
    rw [← Complex.exp_add]; ring_nf
  have h4 : s * ((-x : ℝ) : ℂ) + ell x = (-1) * (s * x - ell x) := by push_cast; ring
  rw [h3, h4, mul_pow]
  push_cast
  ring

lemma P_eq (k : ℕ) (s : ℂ) :
    Q05.P k s = 1 / (k.factorial : ℂ) * ∫ x in Ioi (0 : ℝ), Pint k s x := by
  set G : ℝ → ℂ := fun t => (t : ℂ) ^ s / (1 + (t : ℂ)) *
    (s * (Real.log t : ℂ) + (Real.log (1 + t) : ℂ)) ^ k with hG
  have hder : ∀ x ∈ Ioi (0 : ℝ), HasDerivWithinAt (fun u : ℝ => Real.exp (-u))
      (-Real.exp (-x)) (Ioi 0) x := by
    intro x _
    have := (hasDerivAt_neg x).exp
    exact (by simpa [mul_comm] using this : HasDerivAt (fun u : ℝ => Real.exp (-u))
      (-Real.exp (-x)) x).hasDerivWithinAt
  have hinj : InjOn (fun u : ℝ => Real.exp (-u)) (Ioi 0) := by
    intro x _ y _ h
    have := Real.exp_injective h
    linarith
  have hEq := integral_image_eq_integral_abs_deriv_smul measurableSet_Ioi hder hinj G
  rw [image_exp_neg] at hEq
  have hP : Q05.P k s = (-1) ^ k / (k.factorial : ℂ) * ∫ t in Ioo (0 : ℝ) 1, G t := by
    rw [Q05.P, intervalIntegral.integral_of_le zero_le_one, integral_Ioc_eq_integral_Ioo]
  rw [hP, hEq]
  have hc : EqOn (fun x : ℝ => |(-Real.exp (-x))| • G (Real.exp (-x)))
      (fun x : ℝ => (-1) ^ k * Pint k s x) (Ioi 0) := by
    intro x _
    exact P_integrand k s x
  rw [setIntegral_congr_fun measurableSet_Ioi hc, integral_const_mul]
  have : ((-1 : ℂ) ^ k) * (-1) ^ k = 1 := by rw [← mul_pow]; norm_num
  rw [div_eq_mul_inv, one_div]
  linear_combination ((k.factorial : ℂ)⁻¹ * ∫ x in Ioi (0 : ℝ), Pint k s x) * this

/-! ## `Q_k` as an integral, the decomposition of the bracket, and the term bounds -/

lemma sig_cone (s c : ℂ) (hcr : 0 ≤ c.re) (hci : c.im = 0) (hs : 0 < s.re) (h0 : 0 ≤ s.im)
    (h1 : s.im ≤ s.re) : 0 < (s + c).re ∧ |(s + c).im| ≤ (s + c).re := by
  simp only [Complex.add_re, Complex.add_im, hci, add_zero]
  refine ⟨by linarith, ?_⟩
  rw [abs_of_nonneg h0]; linarith

/-- The integrand of `Q_k`. -/
noncomputable def Qint (k : ℕ) (s L : ℂ) : ℂ := Complex.exp (-((s + 1) * L)) * (s * L) ^ k

lemma Q_eq (k : ℕ) (s : ℂ) (hs : 0 < s.re) (h0 : 0 ≤ s.im) (h1 : s.im ≤ s.re) :
    IntegrableOn (fun x : ℝ => Qint k s x) (Ioi 0) ∧
      Q05.Q k s = 1 / (k.factorial : ℂ) * ∫ x in Ioi (0 : ℝ), Qint k s x := by
  obtain ⟨hr, hi⟩ := sig_cone s 1 (by simp) (by simp) hs h0 h1
  obtain ⟨hI, hV⟩ := cgamma k (s + 1) hr hi
  have hfun : (fun x : ℝ => Qint k s x) =
      fun x : ℝ => s ^ k * ((x : ℂ) ^ k * Complex.exp (-((s + 1) * x))) := by
    funext x; simp only [Qint]; ring
  rw [hfun]
  refine ⟨hI.const_mul _, ?_⟩
  rw [integral_const_mul, hV, Q05.Q]
  have hne : s + 1 ≠ 0 := fun h => by
    have := congrArg Complex.re h; simp at this; linarith
  have hf : (k.factorial : ℂ) ≠ 0 := by exact_mod_cast (Nat.factorial_pos k).ne'
  rw [add_comm (1 : ℂ) s]
  field_simp

/-- The term `i = 0` of the bracket. -/
noncomputable def f0 (k : ℕ) (s L : ℂ) : ℂ :=
  -(Complex.exp (-((s + 2) * L)) * (s * L) ^ k / (1 + Complex.exp (-L)))

/-- The term `i ≥ 1` of the bracket. -/
noncomputable def fi (k : ℕ) (s : ℂ) (i : ℕ) (L : ℂ) : ℂ :=
  (k.choose i : ℂ) * Complex.exp (-((s + 1) * L)) * (s * L) ^ (k - i) * (-ell L) ^ i /
    (1 + Complex.exp (-L))

lemma decomp (k : ℕ) (s L : ℂ) (hL : 1 + Complex.exp (-L) ≠ 0) :
    Pint k s L - Qint k s L = f0 k s L + ∑ i ∈ Finset.Icc 1 k, fi k s i L := by
  have hadd : (s * L - ell L) ^ k = ∑ m ∈ Finset.range (k + 1),
      (-ell L) ^ m * (s * L) ^ (k - m) * (k.choose m : ℂ) := by
    rw [← add_pow]; ring
  have hsplit : ∑ m ∈ Finset.range (k + 1), (-ell L) ^ m * (s * L) ^ (k - m) * (k.choose m : ℂ) =
      (s * L) ^ k + ∑ m ∈ Finset.Icc 1 k, (-ell L) ^ m * (s * L) ^ (k - m) * (k.choose m : ℂ) := by
    rw [Finset.range_eq_Ico, Finset.sum_eq_sum_Ico_succ_bot (by omega)]
    simp only [pow_zero, one_mul, Nat.sub_zero, Nat.choose_zero_right, Nat.cast_one, mul_one]
    rfl
  have h3 : ∑ i ∈ Finset.Icc 1 k, fi k s i L = Complex.exp (-((s + 1) * L)) /
      (1 + Complex.exp (-L)) *
        ∑ m ∈ Finset.Icc 1 k, (-ell L) ^ m * (s * L) ^ (k - m) * (k.choose m : ℂ) := by
    rw [Finset.mul_sum]
    refine Finset.sum_congr rfl (fun i _ => ?_)
    unfold fi; ring
  have h2 : Complex.exp (-((s + 2) * L)) = Complex.exp (-((s + 1) * L)) * Complex.exp (-L) := by
    rw [← Complex.exp_add]; ring_nf
  rw [h3]
  unfold Pint Qint f0
  rw [hadd, hsplit, h2]
  field_simp
  ring

lemma f0_diff (k : ℕ) (s : ℂ) : ∀ L ∈ Dom, DifferentiableAt ℂ (f0 k s) L := by
  intro L hL
  have hne := Complex.slitPlane_ne_zero (dom_slit L hL)
  unfold f0
  fun_prop (disch := exact hne)

lemma fi_diff (k : ℕ) (s : ℂ) (i : ℕ) : ∀ L ∈ Dom, DifferentiableAt ℂ (fi k s i) L := by
  intro L hL
  have hsl := dom_slit L hL
  have hne := Complex.slitPlane_ne_zero hsl
  have hl : DifferentiableAt ℂ ell L := by
    unfold ell
    exact (by fun_prop : DifferentiableAt ℂ (fun L : ℂ => 1 + Complex.exp (-L)) L).clog hsl
  unfold fi
  fun_prop (disch := exact hne)

lemma f0_bound (k : ℕ) (s : ℂ) : ∀ L : ℂ, |L.im| ≤ L.re →
    ‖f0 k s L‖ ≤ (2 * ‖s‖ ^ k) * ‖L‖ ^ k * Real.exp (-((s + 2) * L).re) := by
  intro L hL
  have hinv := inv_one_add_le L hL
  unfold f0
  rw [norm_neg, div_eq_mul_inv, norm_mul, norm_mul, norm_pow, norm_mul, Complex.norm_exp,
    Complex.neg_re, mul_pow]
  calc Real.exp (-((s + 2) * L).re) * (‖s‖ ^ k * ‖L‖ ^ k) * ‖(1 + Complex.exp (-L))⁻¹‖
      ≤ Real.exp (-((s + 2) * L).re) * (‖s‖ ^ k * ‖L‖ ^ k) * 2 :=
        mul_le_mul_of_nonneg_left hinv (by positivity)
    _ = _ := by ring

lemma fi_bound (k : ℕ) (s : ℂ) (i : ℕ) : ∀ L : ℂ, |L.im| ≤ L.re →
    ‖fi k s i L‖ ≤ (2 * (k.choose i : ℝ) * 2 ^ i * ‖s‖ ^ (k - i)) * ‖L‖ ^ (k - i) *
      Real.exp (-((s + ((i : ℂ) + 1)) * L).re) := by
  intro L hL
  have hl := ell_le L hL
  have hinv := inv_one_add_le L hL
  have hli : ‖ell L‖ ^ i ≤ (2 * Real.exp (-L.re)) ^ i := pow_le_pow_left₀ (norm_nonneg _) hl i
  have hre : ((s + ((i : ℂ) + 1)) * L).re = ((s + 1) * L).re + i * L.re := by
    simp [Complex.mul_re]; ring
  have hexp : Real.exp (-((s + 1) * L).re) * (2 * Real.exp (-L.re)) ^ i =
      2 ^ i * Real.exp (-((s + ((i : ℂ) + 1)) * L).re) := by
    rw [mul_pow, ← Real.exp_nat_mul, hre,
      show -(((s + 1) * L).re + i * L.re) = -((s + 1) * L).re + i * (-L.re) by ring, Real.exp_add]
    ring
  unfold fi
  rw [div_eq_mul_inv]
  simp only [norm_mul, norm_pow, norm_neg, Complex.norm_natCast, Complex.norm_exp, Complex.neg_re]
  calc (k.choose i : ℝ) * Real.exp (-((s + 1) * L).re) * (‖s‖ * ‖L‖) ^ (k - i) * ‖ell L‖ ^ i *
        ‖(1 + Complex.exp (-L))⁻¹‖
      ≤ (k.choose i : ℝ) * Real.exp (-((s + 1) * L).re) * (‖s‖ * ‖L‖) ^ (k - i) *
        (2 * Real.exp (-L.re)) ^ i * 2 := by gcongr
    _ = _ := by
      rw [mul_pow]
      linear_combination ((k.choose i : ℝ) * (‖s‖ ^ (k - i) * ‖L‖ ^ (k - i)) * 2) * hexp

/-! ## Assembly -/

lemma W_eq (k : ℕ) (s : ℂ) (hs : 0 < s.re) (h0 : 0 ≤ s.im) (h1 : s.im ≤ s.re) :
    Q05.W k s = 1 / (k.factorial : ℂ) * ((∫ x in Ioi (0 : ℝ), f0 k s x) +
      ∑ i ∈ Finset.Icc 1 k, ∫ x in Ioi (0 : ℝ), fi k s i x) := by
  obtain ⟨c0r, c0i⟩ := sig_cone s 2 (by norm_num) (by simp) hs h0 h1
  have hR0 := rot_bound (f0 k s) (s + 2) k (2 * ‖s‖ ^ k) (by positivity) c0r c0i (f0_diff k s)
    (f0_bound k s)
  have hRi : ∀ i ∈ Finset.Icc 1 k, IntegrableOn (fun x : ℝ => fi k s i x) (Ioi 0) := by
    intro i _
    obtain ⟨cr, ci⟩ := sig_cone s ((i : ℂ) + 1) (by simp; positivity) (by simp) hs h0 h1
    exact (rot_bound (fi k s i) (s + ((i : ℂ) + 1)) (k - i)
      (2 * (k.choose i : ℝ) * 2 ^ i * ‖s‖ ^ (k - i)) (by positivity) cr ci (fi_diff k s i)
      (fi_bound k s i)).1
  obtain ⟨hQI, hQ⟩ := Q_eq k s hs h0 h1
  have hsumI : IntegrableOn (fun x : ℝ => ∑ i ∈ Finset.Icc 1 k, fi k s i x) (Ioi 0) :=
    integrable_finsetSum _ (fun i hi => hRi i hi)
  have hDI : IntegrableOn (fun x : ℝ => f0 k s x + ∑ i ∈ Finset.Icc 1 k, fi k s i x) (Ioi 0) :=
    hR0.1.add hsumI
  have hPeq : ∫ x in Ioi (0 : ℝ), Pint k s x =
      ∫ x in Ioi (0 : ℝ), ((f0 k s x + ∑ i ∈ Finset.Icc 1 k, fi k s i x) + Qint k s x) := by
    refine setIntegral_congr_fun measurableSet_Ioi (fun x hx => ?_)
    have hx0 : (0 : ℝ) ≤ x := le_of_lt hx
    have hne := Complex.slitPlane_ne_zero (slit_of_cone (x : ℂ) (by simpa using hx0))
    rw [← decomp k s x hne]
    ring
  rw [integral_add hDI hQI, integral_add hR0.1 hsumI, integral_finsetSum _ hRi] at hPeq
  rw [Q05.W, P_eq, hQ, hPeq]
  ring

lemma W_bound (k : ℕ) (s : ℂ) (hs : 0 < s.re) (h0 : 0 ≤ s.im) (h1 : s.im ≤ s.re) :
    ‖Q05.W k s‖ ≤ 2 * ‖s‖ ^ k / ‖s + 2‖ ^ (k + 1) +
      2 * ∑ i ∈ Finset.Icc 1 k,
        (2 : ℝ) ^ i / (Nat.factorial i : ℝ) * (‖s‖ ^ (k - i) / ‖s + (i + 1 : ℂ)‖ ^ (k - i + 1)) := by
  obtain ⟨c0r, c0i⟩ := sig_cone s 2 (by norm_num) (by simp) hs h0 h1
  have hR0 := (rot_bound (f0 k s) (s + 2) k (2 * ‖s‖ ^ k) (by positivity) c0r c0i (f0_diff k s)
    (f0_bound k s)).2
  have hRi : ∀ i ∈ Finset.Icc 1 k, ‖∫ x in Ioi (0 : ℝ), fi k s i x‖ ≤
      (2 * (k.choose i : ℝ) * 2 ^ i * ‖s‖ ^ (k - i)) * ((k - i).factorial : ℝ) /
        ‖s + ((i : ℂ) + 1)‖ ^ (k - i + 1) := by
    intro i _
    obtain ⟨cr, ci⟩ := sig_cone s ((i : ℂ) + 1) (by simp; positivity) (by simp) hs h0 h1
    exact (rot_bound (fi k s i) (s + ((i : ℂ) + 1)) (k - i)
      (2 * (k.choose i : ℝ) * 2 ^ i * ‖s‖ ^ (k - i)) (by positivity) cr ci (fi_diff k s i)
      (fi_bound k s i)).2
  have hkf : (0 : ℝ) < k.factorial := by exact_mod_cast Nat.factorial_pos k
  rw [W_eq k s hs h0 h1, norm_mul]
  have hn1 : ‖(1 : ℂ) / (k.factorial : ℂ)‖ = 1 / (k.factorial : ℝ) := by
    rw [norm_div, norm_one, Complex.norm_natCast]
  rw [hn1]
  have hsum : ‖(∫ x in Ioi (0 : ℝ), f0 k s x) + ∑ i ∈ Finset.Icc 1 k, ∫ x in Ioi (0 : ℝ), fi k s i x‖
      ≤ (2 * ‖s‖ ^ k) * (k.factorial : ℝ) / ‖s + 2‖ ^ (k + 1) +
        ∑ i ∈ Finset.Icc 1 k, (2 * (k.choose i : ℝ) * 2 ^ i * ‖s‖ ^ (k - i)) *
          ((k - i).factorial : ℝ) / ‖s + ((i : ℂ) + 1)‖ ^ (k - i + 1) := by
    refine (norm_add_le _ _).trans (add_le_add hR0 ?_)
    exact (norm_sum_le _ _).trans (Finset.sum_le_sum hRi)
  refine (mul_le_mul_of_nonneg_left hsum (by positivity)).trans (le_of_eq ?_)
  rw [mul_add, Finset.mul_sum, Finset.mul_sum]
  congr 1
  · field_simp
  · refine Finset.sum_congr rfl (fun i hi => ?_)
    have hik : i ≤ k := (Finset.mem_Icc.1 hi).2
    have hc : (k.choose i : ℝ) * (i.factorial : ℝ) * ((k - i).factorial : ℝ) = k.factorial := by
      exact_mod_cast Nat.choose_mul_factorial_mul_factorial hik
    have hif : (0 : ℝ) < i.factorial := by exact_mod_cast Nat.factorial_pos i
    have hc' : (k.choose i : ℝ) * ((k - i).factorial : ℝ) = k.factorial / i.factorial := by
      rw [eq_div_iff hif.ne']; linear_combination hc
    calc 1 / (k.factorial : ℝ) * (2 * (k.choose i : ℝ) * 2 ^ i * ‖s‖ ^ (k - i) *
          ((k - i).factorial : ℝ) / ‖s + ((i : ℂ) + 1)‖ ^ (k - i + 1))
        = 1 / (k.factorial : ℝ) * (2 * 2 ^ i * ‖s‖ ^ (k - i) / ‖s + ((i : ℂ) + 1)‖ ^ (k - i + 1)) *
          ((k.choose i : ℝ) * ((k - i).factorial : ℝ)) := by ring
      _ = _ := by
        rw [hc']
        field_simp

lemma cone_of_arg (s : ℂ) (hs : 0 < s.re) (ha0 : 0 ≤ s.arg) (ha1 : s.arg ≤ π / 4) :
    0 ≤ s.im ∧ s.im ≤ s.re := by
  have h0 : 0 ≤ s.im := Complex.arg_nonneg_iff.1 ha0
  refine ⟨h0, ?_⟩
  have hn : 0 < ‖s‖ := lt_of_lt_of_le hs (Complex.re_le_norm s)
  have hsin : Real.sin s.arg ≤ Real.sin (π / 4) :=
    Real.sin_le_sin_of_le_of_le_pi_div_two (by linarith [Real.pi_pos]) (by linarith [Real.pi_pos])
      ha1
  rw [Complex.sin_arg, Real.sin_pi_div_four, div_le_iff₀ hn] at hsin
  have hsq : ‖s‖ ^ 2 = s.re ^ 2 + s.im ^ 2 := by
    rw [Complex.sq_norm, Complex.normSq_apply]; ring
  have h2 : √2 * √2 = 2 := Real.mul_self_sqrt (by norm_num)
  have h3 : s.im ^ 2 ≤ (√2 / 2 * ‖s‖) ^ 2 := pow_le_pow_left₀ h0 hsin 2
  have h4 : s.im ^ 2 ≤ s.re ^ 2 := by nlinarith
  nlinarith

theorem main : Q07.stmt_3 := by
  intro k _ s hs ha0 ha1
  obtain ⟨h0, h1⟩ := cone_of_arg s hs ha0 ha1
  exact W_bound k s hs h0 h1

end Q07.PC

/-! ### Merged from work/Q07/PartD.lean -/
/-!
# Q07 Part D: `stmt_2 → stmt_3 → stmt_4` (g5:prop:rest, the displayed inequality)
-/

namespace Q07.PD

open Real MeasureTheory Set Filter Q05

/-- `2^i / i!`. -/
noncomputable def cc (i : ℕ) : ℝ := (2 : ℝ) ^ i / (Nat.factorial i : ℝ)

/-- The majorant of `|HKP - HEQ|` on the half line. -/
noncomputable def G (k : ℕ) (r : ℝ) : ℝ :=
  8 * FjqIntegrand k 2 a r +
    ∑ i ∈ Finset.Icc 1 k, (8 * cc i) * FjqIntegrand (k - i) ((i : ℝ) + 1) a r +
    4 * FjqIntegrand k 1 (3 * a) r

lemma arg_Lpt (r : ℝ) (hr : 0 ≤ r) : 0 ≤ (Lpt r).arg ∧ (Lpt r).arg ≤ π / 4 := by
  have hre := Lpt_re r
  have him := Lpt_im r
  have hs : 0 ≤ r * (√2 / 2) := by positivity
  have hrep : 0 < (Lpt r).re := by rw [hre]; linarith
  constructor
  · rw [Complex.arg_nonneg_iff, him]; exact hs
  · have h1 := Complex.abs_arg_lt_pi_div_two_iff.2 (Or.inl hrep)
    rw [abs_lt] at h1
    have ht := Complex.tan_arg (Lpt r)
    have he : (Lpt r).arg = Real.arctan ((Lpt r).im / (Lpt r).re) := by
      rw [← ht, Real.arctan_tan h1.1 h1.2]
    rw [he, ← Real.arctan_one]
    apply Real.arctan_strictMono.monotone
    rw [div_le_one hrep, hre, him]; linarith

lemma norm_ω : ‖ω‖ = 1 := by
  rw [ω, Complex.norm_exp]; simp

lemma F2_eq (k : ℕ) (r : ℝ) : FjqIntegrand k 2 a r =
    Real.exp (-a * r) * (‖Lpt r‖ ^ k / ‖Lpt r + 2‖ ^ (k + 1)) := by
  unfold FjqIntegrand; push_cast; ring

lemma Fi_eq (k i : ℕ) (r : ℝ) : FjqIntegrand (k - i) ((i : ℝ) + 1) a r =
    Real.exp (-a * r) * (‖Lpt r‖ ^ (k - i) / ‖Lpt r + ((i : ℂ) + 1)‖ ^ (k - i + 1)) := by
  unfold FjqIntegrand; push_cast; ring

lemma F1_eq (k : ℕ) (r : ℝ) : FjqIntegrand k 1 (3 * a) r =
    Real.exp (-3 * a * r) * ‖Q k (Lpt r)‖ := by
  unfold FjqIntegrand
  rw [Q_le', add_comm (1 : ℂ) (Lpt r)]
  push_cast
  ring_nf

lemma pt (h3 : stmt_3) (k : ℕ) (hk : 1 ≤ k) (r : ℝ) (hr : 0 ≤ r) :
    ‖H (Lpt r) * K (Lpt r) * P k (Lpt r) * ω - fM k (Lpt r) * ω‖ ≤ G k r := by
  obtain ⟨hH, -, hK, hKE⟩ := stmt7_at r hr
  have hre := Lpt_re_pos r hr
  obtain ⟨ha0, ha1⟩ := arg_Lpt r hr
  have hW := h3 k hk (Lpt r) hre ha0 ha1
  set s := Lpt r
  have heq : H s * K s * P k s * ω - fM k s * ω =
      (H s * K s * W k s + H s * (K s - E s) * Q k s) * ω := by
    unfold fM W; ring
  rw [heq, norm_mul, norm_ω, mul_one]
  refine (norm_add_le _ _).trans ?_
  rw [norm_mul, norm_mul, norm_mul, norm_mul]
  set Wb := 2 * ‖s‖ ^ k / ‖s + 2‖ ^ (k + 1) +
      2 * ∑ i ∈ Finset.Icc 1 k,
        (2 : ℝ) ^ i / (Nat.factorial i : ℝ) * (‖s‖ ^ (k - i) / ‖s + (i + 1 : ℂ)‖ ^ (k - i + 1))
  have h1 : ‖H s‖ * ‖K s‖ * ‖W k s‖ ≤ 2 * (2 * Real.exp (-a * r)) * Wb := by
    gcongr
  have h2 : ‖H s‖ * ‖K s - E s‖ * ‖Q k s‖ ≤ 2 * (2 * Real.exp (-3 * a * r)) * ‖Q k s‖ := by
    gcongr
  refine (add_le_add h1 h2).trans (le_of_eq ?_)
  unfold G
  rw [F2_eq, F1_eq, Finset.sum_congr rfl (fun i _ => by rw [Fi_eq])]
  simp only [Wb, cc]
  rw [Finset.mul_sum, mul_add, Finset.mul_sum]
  have : ∀ i ∈ Finset.Icc 1 k, 2 * (2 * Real.exp (-a * r)) *
      (2 * ((2 : ℝ) ^ i / (Nat.factorial i : ℝ) *
        (‖s‖ ^ (k - i) / ‖s + (i + 1 : ℂ)‖ ^ (k - i + 1)))) =
      8 * ((2 : ℝ) ^ i / (Nat.factorial i : ℝ)) *
        (Real.exp (-a * r) * (‖s‖ ^ (k - i) / ‖s + ((i : ℂ) + 1)‖ ^ (k - i + 1))) := by
    intro i _; ring
  rw [Finset.sum_congr rfl this]
  ring

lemma Fint (j : ℕ) (q l : ℝ) (hq : 1 ≤ q) (hl : 0 < l) :
    IntegrableOn (FjqIntegrand j q l) (Ioi 0) :=
  (Q05.proof_13 j q l hq hl).1

lemma Fi_int (k i : ℕ) : IntegrableOn
    (fun r => (8 * cc i) * FjqIntegrand (k - i) ((i : ℝ) + 1) a r) (Ioi 0) :=
  (Fint _ _ _ (by have : (0:ℝ) ≤ i := Nat.cast_nonneg i; linarith) a_pos).const_mul _

lemma Fsum_int (k : ℕ) : IntegrableOn
    (fun r => ∑ i ∈ Finset.Icc 1 k, (8 * cc i) * FjqIntegrand (k - i) ((i : ℝ) + 1) a r)
    (Ioi 0) :=
  integrable_finsetSum _ (fun i _ => Fi_int k i)

lemma F2_int (k : ℕ) : IntegrableOn (fun r => 8 * FjqIntegrand k 2 a r) (Ioi 0) :=
  (Fint k 2 a (by norm_num) a_pos).const_mul 8

lemma F1_int (k : ℕ) : IntegrableOn (fun r => 4 * FjqIntegrand k 1 (3 * a) r) (Ioi 0) :=
  (Fint k 1 (3 * a) le_rfl (by linarith [a_pos])).const_mul 4

lemma G_int (k : ℕ) : IntegrableOn (G k) (Ioi 0) :=
  Integrable.add (Integrable.add (F2_int k) (Fsum_int k)) (F1_int k)

lemma G_integral (k : ℕ) : ∫ r in Ioi (0:ℝ), G k r =
    8 * Fjq k 2 a + ∑ i ∈ Finset.Icc 1 k, (8 * cc i) * Fjq (k - i) ((i : ℝ) + 1) a +
      4 * Fjq k 1 (3 * a) := by
  unfold G
  set f1 : ℝ → ℝ := fun r => 8 * FjqIntegrand k 2 a r
  set f2 : ℝ → ℝ := fun r => ∑ i ∈ Finset.Icc 1 k, (8 * cc i) * FjqIntegrand (k - i) ((i : ℝ) + 1) a r
  set f3 : ℝ → ℝ := fun r => 4 * FjqIntegrand k 1 (3 * a) r
  have e1 : ∫ r in Ioi (0:ℝ), (f1 r + f2 r + f3 r) =
      (∫ r in Ioi (0:ℝ), (f1 r + f2 r)) + ∫ r in Ioi (0:ℝ), f3 r :=
    integral_add (Integrable.add (F2_int k) (Fsum_int k)) (F1_int k)
  have e2 : ∫ r in Ioi (0:ℝ), (f1 r + f2 r) =
      (∫ r in Ioi (0:ℝ), f1 r) + ∫ r in Ioi (0:ℝ), f2 r :=
    integral_add (F2_int k) (Fsum_int k)
  have e3 : ∫ r in Ioi (0:ℝ), f2 r = ∑ i ∈ Finset.Icc 1 k,
      ∫ r in Ioi (0:ℝ), (8 * cc i) * FjqIntegrand (k - i) ((i : ℝ) + 1) a r :=
    integral_finsetSum _ (fun i _ => Fi_int k i)
  show ∫ r in Ioi (0:ℝ), (f1 r + f2 r + f3 r) = _
  rw [e1, e2, e3]
  simp only [f1, f3]
  rw [integral_const_mul, integral_const_mul]
  unfold Fjq
  congr 2
  exact Finset.sum_congr rfl (fun i _ => integral_const_mul _ _)

lemma line_bound (h2 : stmt_2) (h3 : stmt_3) (k : ℕ) (hk : 1 ≤ k) :
    ‖lineInt (fun s => H s * K s * P k s) - lineInt (fM k)‖ ≤ ∫ r in Ioi (0:ℝ), G k r := by
  have hA := ((h2 k).1).mul_const ω
  have hB := (fM_line_int k).mul_const ω
  unfold lineInt
  rw [← integral_sub hA hB]
  apply norm_integral_le_of_norm_le (G_int k)
  refine (ae_restrict_iff' measurableSet_Ioi).2 (Eventually.of_forall fun r hr => ?_)
  exact pt h3 k hk r (le_of_lt hr)

lemma a_div : a / √2 = π / 4 := by
  rw [a, div_div, mul_assoc, sqrt2_sq]; ring

lemma Fa (j : ℕ) (q : ℝ) (hq : 1 ≤ q) :
    16 / π * Fjq j q a ≤
      B1 * Real.exp (π / 4 * (q - 1)) * Real.exp (-√(3 * π * j * q / 4)) := by
  have h := (Q05.proof_13 j q a hq a_pos).2
  have e1 : a * q / √2 = π / 4 * q := by
    rw [show a * q / √2 = q * (a / √2) by ring, a_div]; ring
  have e2 : 3 * a * j * q / √2 = 3 * π * j * q / 4 := by
    rw [show 3 * a * j * q / √2 = 3 * j * q * (a / √2) by ring, a_div]; ring
  rw [e1, e2] at h
  have hpi : 0 < 16 / π := by have := Real.pi_pos; positivity
  calc 16 / π * Fjq j q a ≤ 16 / π * (4 / a * Real.exp (a / 2 + π / 4 * q) *
        Real.exp (-√(3 * π * j * q / 4))) := mul_le_mul_of_nonneg_left h hpi.le
    _ = _ := by
        unfold B1
        have : Real.exp (a / 2 + π / 4 * q) =
            Real.exp (a / 2 + π / 4) * Real.exp (π / 4 * (q - 1)) := by
          rw [← Real.exp_add]; ring_nf
        rw [this]
        have := a_pos; have := Real.pi_pos
        field_simp
        ring

lemma F3a (j : ℕ) : 8 / π * Fjq j 1 (3 * a) ≤ B2 * Real.exp (-√(9 * π * j / 4)) := by
  have h := (Q05.proof_13 j 1 (3 * a) le_rfl (by linarith [a_pos])).2
  have e1 : 3 * a / 2 + 3 * a * 1 / √2 = 3 * a / 2 + 3 * π / 4 := by
    rw [show 3 * a * 1 / √2 = 3 * (a / √2) by ring, a_div]; ring
  have e2 : 3 * (3 * a) * j * 1 / √2 = 9 * π * j / 4 := by
    rw [show 3 * (3 * a) * j * 1 / √2 = 9 * j * (a / √2) by ring, a_div]; ring
  rw [e1, e2] at h
  have hpi : 0 < 8 / π := by have := Real.pi_pos; positivity
  calc 8 / π * Fjq j 1 (3 * a) ≤ 8 / π * (4 / (3 * a) * Real.exp (3 * a / 2 + 3 * π / 4) *
        Real.exp (-√(9 * π * j / 4))) := mul_le_mul_of_nonneg_left h hpi.le
    _ = _ := by
        unfold B2
        have := a_pos; have := Real.pi_pos
        field_simp
        ring

lemma xi_pow (i : ℕ) : cc i * Real.exp (π / 4 * (((i : ℝ) + 1) - 1)) =
    ξ ^ i / (Nat.factorial i : ℝ) := by
  unfold cc ξ
  rw [show π / 4 * (((i : ℝ) + 1) - 1) = (i : ℝ) * (π / 4) by ring, Real.exp_nat_mul, mul_pow]
  ring

lemma U_bound (h2 : stmt_2) (h3 : stmt_3) (k : ℕ) (hk : 1 ≤ k) :
    2 / π * ‖lineInt (fun s => H s * K s * P k s) - lineInt (fM k)‖ ≤
      B1 * (Real.exp (π / 4) * Real.exp (-√(3 * π * k * 2 / 4))) +
      B1 * ∑ i ∈ Finset.Icc 1 k, ξ ^ i / (Nat.factorial i : ℝ) *
        Real.exp (-√(3 * π * ((k - i : ℕ) : ℝ) * ((i : ℝ) + 1) / 4)) +
      B2 * Real.exp (-√(9 * π * k / 4)) := by
  have hpi : 0 < 2 / π := by have := Real.pi_pos; positivity
  have hL := mul_le_mul_of_nonneg_left (line_bound h2 h3 k hk) hpi.le
  rw [G_integral] at hL
  refine hL.trans ?_
  have eq : 2 / π * (8 * Fjq k 2 a + ∑ i ∈ Finset.Icc 1 k, (8 * cc i) * Fjq (k - i) ((i : ℝ) + 1) a +
      4 * Fjq k 1 (3 * a)) = 16 / π * Fjq k 2 a +
      ∑ i ∈ Finset.Icc 1 k, cc i * (16 / π * Fjq (k - i) ((i : ℝ) + 1) a) +
      8 / π * Fjq k 1 (3 * a) := by
    rw [mul_add, mul_add, Finset.mul_sum]
    have : ∀ i ∈ Finset.Icc 1 k, 2 / π * (8 * cc i * Fjq (k - i) ((i : ℝ) + 1) a) =
        cc i * (16 / π * Fjq (k - i) ((i : ℝ) + 1) a) := fun i _ => by ring
    rw [Finset.sum_congr rfl this]; ring
  rw [eq]
  have b1 := Fa k 2 (by norm_num)
  rw [show π / 4 * ((2 : ℝ) - 1) = π / 4 by ring,
    show 3 * π * (k : ℝ) * 2 / 4 = 3 * π * k * 2 / 4 by ring] at b1
  have b2 := F3a k
  have b3 : ∑ i ∈ Finset.Icc 1 k, cc i * (16 / π * Fjq (k - i) ((i : ℝ) + 1) a) ≤
      B1 * ∑ i ∈ Finset.Icc 1 k, ξ ^ i / (Nat.factorial i : ℝ) *
        Real.exp (-√(3 * π * ((k - i : ℕ) : ℝ) * ((i : ℝ) + 1) / 4)) := by
    rw [Finset.mul_sum]
    apply Finset.sum_le_sum
    intro i _
    have hc : 0 ≤ cc i := by unfold cc; positivity
    have hb := mul_le_mul_of_nonneg_left
      (Fa (k - i) ((i : ℝ) + 1) (by have : (0:ℝ) ≤ i := Nat.cast_nonneg i; linarith)) hc
    refine hb.trans (le_of_eq ?_)
    rw [← xi_pow]; ring
  linarith

lemma range_split (k : ℕ) (c : ℕ → ℝ) :
    ∑ i ∈ Finset.range (k + 1), c i = c 0 + ∑ i ∈ Finset.Icc 1 k, c i := by
  rw [Finset.range_eq_Ico, Finset.sum_eq_sum_Ico_succ_bot (by omega)]
  rfl

lemma xi_pos : 0 < ξ := by unfold ξ; positivity

lemma sum_Icc_le_exp (k : ℕ) (x : ℝ) (hx : 0 ≤ x) :
    ∑ i ∈ Finset.Icc 1 k, x ^ i / (Nat.factorial i : ℝ) ≤ Real.exp x - 1 := by
  have h := Real.sum_le_exp_of_nonneg hx (k + 1)
  rw [range_split] at h
  simp at h
  linarith

lemma sum_ge2 (k : ℕ) (hk : 1 ≤ k) (x : ℝ) (hx : 0 ≤ x) :
    ∑ i ∈ Finset.Icc 1 k, (if 2 ≤ i then x ^ i / (Nat.factorial i : ℝ) else 0) ≤
      Real.exp x - 1 - x := by
  have h1 := sum_Icc_le_exp k x hx
  have h2 : ∑ i ∈ Finset.Icc 1 k, x ^ i / (Nat.factorial i : ℝ) -
      ∑ i ∈ Finset.Icc 1 k, (if 2 ≤ i then x ^ i / (Nat.factorial i : ℝ) else 0) = x := by
    rw [← Finset.sum_sub_distrib]
    have : ∀ i ∈ Finset.Icc 1 k, x ^ i / (Nat.factorial i : ℝ) -
        (if 2 ≤ i then x ^ i / (Nat.factorial i : ℝ) else 0) =
        if 1 = i then x ^ i / (Nat.factorial i : ℝ) else 0 := by
      intro i hi
      rw [Finset.mem_Icc] at hi
      by_cases h : i = 1
      · subst h; simp
      · have : 2 ≤ i := by omega
        rw [if_pos this, if_neg (Ne.symm h)]; ring
    rw [Finset.sum_congr rfl this, Finset.sum_ite_eq]
    rw [if_pos (by rw [Finset.mem_Icc]; omega)]
    simp
  linarith

/-- The split of the sum over `i`. -/
lemma sum_split (k : ℕ) (hk : 99 ≤ k) :
    ∑ i ∈ Finset.Icc 1 k, ξ ^ i / (Nat.factorial i : ℝ) *
        Real.exp (-√(3 * π * ((k - i : ℕ) : ℝ) * ((i : ℝ) + 1) / 4)) ≤
      ξ * Real.exp (-√(3 * π * ((k : ℝ) - 1) / 2)) +
        (Real.exp ξ - 1 - ξ) * Real.exp (-√(9 * π * ((k : ℝ) - 2) / 4)) +
        Real.exp (Real.exp 1 * ξ) * Real.exp (-(k : ℝ) / 2) := by
  have hxi := xi_pos
  set E1 := Real.exp (-√(3 * π * ((k : ℝ) - 1) / 2))
  set E2 := Real.exp (-√(9 * π * ((k : ℝ) - 2) / 4))
  set E3 := Real.exp (-(k : ℝ) / 2)
  have hE2 : 0 < E2 := Real.exp_pos _
  have hE3 : 0 < E3 := Real.exp_pos _
  have hpt : ∀ i ∈ Finset.Icc 1 k, ξ ^ i / (Nat.factorial i : ℝ) *
      Real.exp (-√(3 * π * ((k - i : ℕ) : ℝ) * ((i : ℝ) + 1) / 4)) ≤
      (if 1 = i then ξ * E1 else 0) +
        (if 2 ≤ i then ξ ^ i / (Nat.factorial i : ℝ) else 0) * E2 +
        E3 * ((Real.exp 1 * ξ) ^ i / (Nat.factorial i : ℝ)) := by
    intro i hi
    rw [Finset.mem_Icc] at hi
    have hc : 0 ≤ ξ ^ i / (Nat.factorial i : ℝ) := by positivity
    have hf3 : 0 ≤ E3 * ((Real.exp 1 * ξ) ^ i / (Nat.factorial i : ℝ)) := by positivity
    have hcast : ((k - i : ℕ) : ℝ) = (k : ℝ) - i := by rw [Nat.cast_sub hi.2]
    rw [hcast]
    by_cases h1 : i = 1
    · subst h1
      rw [if_pos rfl, if_neg (by norm_num)]
      have : 3 * π * ((k : ℝ) - ((1 : ℕ) : ℝ)) * (((1 : ℕ) : ℝ) + 1) / 4 =
          3 * π * ((k : ℝ) - 1) / 2 := by push_cast; ring
      rw [this]
      simp only [pow_one, Nat.factorial_one, Nat.cast_one, div_one]
      linarith
    · have h2 : 2 ≤ i := by omega
      rw [if_neg (Ne.symm h1), if_pos h2]
      by_cases h3 : 2 * i ≤ k
      · have hle : Real.exp (-√(3 * π * ((k : ℝ) - i) * ((i : ℝ) + 1) / 4)) ≤ E2 := by
          apply Real.exp_le_exp.2
          apply neg_le_neg
          apply Real.sqrt_le_sqrt
          have hi2 : (2 : ℝ) ≤ i := by exact_mod_cast h2
          have hik : 2 * (i : ℝ) ≤ k := by exact_mod_cast h3
          have hk' : (99 : ℝ) ≤ k := by exact_mod_cast hk
          have hp := Real.pi_pos
          have : 3 * ((k : ℝ) - 2) ≤ ((k : ℝ) - i) * ((i : ℝ) + 1) := by
            nlinarith [mul_nonneg (by linarith : (0:ℝ) ≤ i - 2) (by linarith : (0:ℝ) ≤ k - i - 3)]
          nlinarith
        have := mul_le_mul_of_nonneg_left hle hc
        linarith
      · have hle : Real.exp (-√(3 * π * ((k : ℝ) - i) * ((i : ℝ) + 1) / 4)) ≤ 1 := by
          rw [Real.exp_le_one_iff]; exact neg_nonpos.2 (Real.sqrt_nonneg _)
        have h4 : ξ ^ i / (Nat.factorial i : ℝ) ≤
            E3 * ((Real.exp 1 * ξ) ^ i / (Nat.factorial i : ℝ)) := by
          rw [mul_pow, Real.exp_one_pow]
          have hik : (k : ℝ) < 2 * i := by exact_mod_cast (not_le.1 h3)
          have he : 1 ≤ E3 * Real.exp i := by
            rw [← Real.exp_add]; apply Real.one_le_exp; linarith
          calc ξ ^ i / (Nat.factorial i : ℝ) ≤ (E3 * Real.exp i) * (ξ ^ i / (Nat.factorial i : ℝ)) :=
                le_mul_of_one_le_left hc he
            _ = _ := by ring
        have := mul_le_mul_of_nonneg_left hle hc
        nlinarith
  refine (Finset.sum_le_sum hpt).trans ?_
  rw [Finset.sum_add_distrib, Finset.sum_add_distrib, Finset.sum_ite_eq,
    if_pos (by rw [Finset.mem_Icc]; omega), ← Finset.sum_mul, ← Finset.mul_sum]
  have s2 := sum_ge2 k (by omega) ξ hxi.le
  have s3 := sum_Icc_le_exp k (Real.exp 1 * ξ) (by positivity)
  have t2 := mul_le_mul_of_nonneg_right s2 hE2.le
  have t3 := mul_le_mul_of_nonneg_left s3 hE3.le
  have : 0 ≤ E3 := hE3.le
  nlinarith

lemma two_pow_le_fact (k : ℕ) : 2 ^ k ≤ Nat.factorial (k + 1) := by
  induction k with
  | zero => simp
  | succ n ih =>
    rw [Nat.factorial_succ, pow_succ]
    have : 2 ≤ n + 1 + 1 := by omega
    calc 2 ^ n * 2 ≤ Nat.factorial (n + 1) * (n + 1 + 1) := Nat.mul_le_mul ih this
      _ = (n + 1 + 1) * Nat.factorial (n + 1) := by ring

lemma half_pow_le (k : ℕ) : (1 / 2 : ℝ) ^ k ≤ Real.exp (-(k : ℝ) / 2) := by
  rw [show -(k : ℝ) / 2 = k * (-1 / 2) by ring, Real.exp_nat_mul]
  apply pow_le_pow_left₀ (by norm_num)
  have := Real.add_one_le_exp (-1 / 2 : ℝ)
  linarith

lemma D_bound (k : ℕ) : |D k| ≤ Real.exp (-(k : ℝ) / 2) := by
  refine le_trans ?_ (half_pow_le k)
  unfold D
  have hl0 : 0 < Real.log 2 := Real.log_pos (by norm_num)
  have hl1 : Real.log 2 ≤ 1 := by
    have := Real.log_two_lt_d9; linarith
  have hf : (2 : ℝ) ^ k ≤ (Nat.factorial (k + 1) : ℝ) := by exact_mod_cast two_pow_le_fact k
  rw [abs_div, abs_mul, abs_pow, abs_neg, abs_one, one_pow, one_mul, abs_pow,
    abs_of_pos hl0, Nat.abs_cast]
  have hp : Real.log 2 ^ (k + 1) ≤ 1 := pow_le_one₀ hl0.le hl1
  rw [div_pow, one_pow]
  calc Real.log 2 ^ (k + 1) / (Nat.factorial (k + 1) : ℝ) ≤ 1 / (Nat.factorial (k + 1) : ℝ) :=
        div_le_div_of_nonneg_right hp (by positivity)
    _ ≤ 1 / 2 ^ k := one_div_le_one_div_of_le (by positivity) hf

lemma U3_bound (k : ℕ) (hk : 1 ≤ k) :
    ‖∫ x in (0 : ℝ)..(1 / 2), fM k x‖ ≤ (1 / 2 : ℝ) ^ k := by
  have hb : ∀ x ∈ Set.uIoc (0 : ℝ) (1 / 2), ‖fM k x‖ ≤ (1 / 2 : ℝ) ^ (k - 1) := by
    intro x hx
    rw [uIoc_of_le (by norm_num)] at hx
    have hx0 : 0 < x := hx.1
    have hH := H_le_one_div (x : ℂ) (by simpa using hx0)
    simp only [Complex.ofReal_re] at hH
    have hE : ‖E (x : ℂ)‖ = 1 := by rw [E_norm]; simp
    have hQ : ‖Q k (x : ℂ)‖ ≤ x ^ k := by
      rw [Q_le', Complex.norm_real, Real.norm_eq_abs, abs_of_pos hx0]
      have h1 : 1 ≤ ‖1 + (x : ℂ)‖ ^ (k + 1) := by
        apply one_le_pow₀; have := one_le_norm_one_add (x : ℂ); simp at this; linarith
      calc x ^ k / ‖1 + (x : ℂ)‖ ^ (k + 1) ≤ x ^ k / 1 :=
            div_le_div_of_nonneg_left (by positivity) one_pos h1
        _ = x ^ k := div_one _
    rw [fM, norm_mul, norm_mul, hE, mul_one]
    have hxk : x ^ k = x * x ^ (k - 1) := by
      rw [← pow_succ']; congr 1; omega
    calc ‖H x‖ * ‖Q k x‖ ≤ 1 / x * x ^ k := mul_le_mul hH hQ (norm_nonneg _) (by positivity)
      _ = x ^ (k - 1) := by rw [hxk]; field_simp
      _ ≤ (1 / 2 : ℝ) ^ (k - 1) := pow_le_pow_left₀ hx0.le hx.2 _
  have := intervalIntegral.norm_integral_le_of_norm_le_const hb
  refine this.trans (le_of_eq ?_)
  rw [sub_zero, abs_of_pos (by norm_num), ← pow_succ]
  congr 1; omega

lemma M_eq (h2 : stmt_2) (k : ℕ) :
    M k = -D k + 2 / π * (lineInt fun s => H s * K s * P k s).im := by
  unfold M; rw [(h2 k).2, Complex.ofReal_re]

lemma Mmain_eq (k : ℕ) (hk : 1 ≤ k) :
    Mmain k = 2 / π * ((∫ x in (0 : ℝ)..(1 / 2), fM k x) + lineInt (fM k)).im := by
  obtain ⟨-, -, hid, hsh⟩ := Q05.proof_11 k hk
  unfold Mmain
  congr 2
  have hc : ∫ s in Ioi (0 : ℝ), H s * E s * (s : ℂ)⁻¹ * ((s : ℂ) / (1 + s)) ^ (k + 1) =
      ∫ s in Ioi (0 : ℝ), H s * E s * Q k s := by
    refine setIntegral_congr_fun measurableSet_Ioi (fun s hs => ?_)
    rw [mul_assoc (H s * E s), hid s hs]
  rw [hc, hsh]
  rfl

lemma rest_core (h2 : stmt_2) (k : ℕ) (hk : 1 ≤ k) :
    |M k - Mmain k| ≤ |D k| +
      2 / π * ‖lineInt (fun s => H s * K s * P k s) - lineInt (fM k)‖ +
      2 / π * ‖∫ x in (0 : ℝ)..(1 / 2), fM k x‖ := by
  rw [M_eq h2, Mmain_eq k hk]
  set A := lineInt fun s => H s * K s * P k s
  set B := lineInt (fM k)
  set C := ∫ x in (0 : ℝ)..(1 / 2), fM k x
  have e : -D k + 2 / π * A.im - 2 / π * (C + B).im = -D k + 2 / π * ((A - B).im - C.im) := by
    simp only [Complex.add_im, Complex.sub_im]; ring
  rw [e]
  have hpi : 0 < 2 / π := by have := Real.pi_pos; positivity
  have h1 := Complex.abs_im_le_norm (A - B)
  have h2' := Complex.abs_im_le_norm C
  calc |-D k + 2 / π * ((A - B).im - C.im)| ≤ |D k| + 2 / π * (|(A - B).im| + |C.im|) := by
        refine (abs_add_le _ _).trans ?_
        rw [abs_neg, abs_mul, abs_of_pos hpi]
        gcongr
        exact abs_sub _ _
    _ ≤ _ := by nlinarith

theorem main : Q07.stmt_2 → Q07.stmt_3 → Q07.stmt_4 := by
  intro h2 h3 k hk
  have hk99 : 99 ≤ k := by omega
  have hk1 : 1 ≤ k := by omega
  refine (rest_core h2 k hk1).trans ?_
  have hU := U_bound h2 h3 k hk1
  have hS := sum_split k hk99
  have hD := D_bound k
  have hU3 := U3_bound k hk1
  have hpi : 0 < 2 / π := by have := Real.pi_pos; positivity
  have hpi1 : 2 / π ≤ 1 := by
    rw [div_le_one Real.pi_pos]; linarith [Real.pi_gt_three]
  have h3' : 2 / π * ‖∫ x in (0 : ℝ)..(1 / 2), fM k x‖ ≤ Real.exp (-(k : ℝ) / 2) := by
    calc 2 / π * ‖∫ x in (0 : ℝ)..(1 / 2), fM k x‖ ≤ 1 * (1 / 2 : ℝ) ^ k :=
          mul_le_mul hpi1 hU3 (norm_nonneg _) zero_le_one
      _ ≤ _ := by rw [one_mul]; exact half_pow_le k
  have hB1 : 0 ≤ B1 := by unfold B1; have := a_pos; have := Real.pi_pos; positivity
  have hS' := mul_le_mul_of_nonneg_left hS hB1
  have e1 : 3 * π * (k : ℝ) * 2 / 4 = 3 * π * k / 2 := by ring
  have e2 : -((((k + 1 : ℕ) : ℝ)) - 1) / 2 = -(k : ℝ) / 2 := by push_cast; ring
  rw [e1] at hU
  rw [e2]
  nlinarith

end Q07.PD

/-! ### Merged from work/Q07/PartE.lean -/
namespace Q07.PE

open Real Q05

/-- If `4N ≤ x` then `e^{-√x} ≤ e^{-2√N}`. -/
lemma exp_neg_sqrt_le (x N : ℝ) (hN : 0 ≤ N) (h : 4 * N ≤ x) :
    Real.exp (-√x) ≤ Real.exp (-2 * √N) := by
  apply Real.exp_le_exp.2
  have h1 : √(4 * N) ≤ √x := Real.sqrt_le_sqrt h
  have h2 : √(4 * N) = 2 * √N := by
    rw [Real.sqrt_mul (by norm_num), show (4 : ℝ) = 2 ^ 2 by norm_num,
      Real.sqrt_sq (by norm_num)]
  linarith

/-- If `0 ≤ y` and `4N ≤ y²` then `e^{-y} ≤ e^{-2√N}`. -/
lemma exp_neg_le (y N : ℝ) (hN : 0 ≤ N) (hy : 0 ≤ y) (h : 4 * N ≤ y ^ 2) :
    Real.exp (-y) ≤ Real.exp (-2 * √N) := by
  have := exp_neg_sqrt_le (y ^ 2) N hN h
  rwa [Real.sqrt_sq hy] at this

lemma B1_pos : 0 < B1 := by
  unfold B1 Q05.a; positivity

lemma B2_pos : 0 < B2 := by
  unfold B2 Q05.a; positivity

lemma xi_nonneg : 0 ≤ ξ := by
  unfold ξ; positivity

lemma expxi_nonneg : 0 ≤ Real.exp ξ - 1 - ξ := by
  have := Real.add_one_le_exp ξ; linarith

lemma pi_gt3 : (3 : ℝ) < π := Real.pi_gt_three

/-- The displayed right-hand side is at most `C_⋆ e^{-2√N}` for `k ≥ 99`. -/
lemma rhs_le (k : ℕ) (hk : 100 ≤ k + 1) :
    B1 * (Real.exp (π / 4) * Real.exp (-√(3 * π * k / 2)) +
          ξ * Real.exp (-√(3 * π * ((k : ℝ) - 1) / 2)) +
          (Real.exp ξ - 1 - ξ) * Real.exp (-√(9 * π * ((k : ℝ) - 2) / 4)) +
          Real.exp (Real.exp 1 * ξ) * Real.exp (-(k : ℝ) / 2)) +
        B2 * Real.exp (-√(9 * π * k / 4)) + 2 * Real.exp (-(((k + 1 : ℕ) : ℝ) - 1) / 2) ≤
      Cstar * Real.exp (-2 * √((k + 1 : ℕ) : ℝ)) := by
  have hk' : (99 : ℝ) ≤ k := by exact_mod_cast (by omega : 99 ≤ k)
  have hp := pi_gt3
  set N : ℝ := ((k + 1 : ℕ) : ℝ) with hNdef
  have hN : N = (k : ℝ) + 1 := by rw [hNdef]; push_cast; ring
  have hN0 : 0 ≤ N := by rw [hN]; linarith
  set e := Real.exp (-2 * √N) with he
  have t1 : Real.exp (-√(3 * π * k / 2)) ≤ e :=
    exp_neg_sqrt_le _ _ hN0 (by rw [hN]; nlinarith)
  have t2 : Real.exp (-√(3 * π * ((k : ℝ) - 1) / 2)) ≤ e :=
    exp_neg_sqrt_le _ _ hN0 (by rw [hN]; nlinarith)
  have t3 : Real.exp (-√(9 * π * ((k : ℝ) - 2) / 4)) ≤ e :=
    exp_neg_sqrt_le _ _ hN0 (by rw [hN]; nlinarith)
  have t4 : Real.exp (-(k : ℝ) / 2) ≤ e := by
    rw [neg_div]
    exact exp_neg_le _ _ hN0 (by linarith) (by rw [hN]; nlinarith)
  have t5 : Real.exp (-√(9 * π * k / 4)) ≤ e :=
    exp_neg_sqrt_le _ _ hN0 (by rw [hN]; nlinarith)
  have t6 : Real.exp (-(N - 1) / 2) ≤ e := by
    rw [hN, show (k : ℝ) + 1 - 1 = k by ring]; exact t4
  have hB1 := B1_pos
  have hB2 := B2_pos
  have hxi := xi_nonneg
  have hex := expxi_nonneg
  have h41 : 0 ≤ Real.exp (π / 4) := (Real.exp_pos _).le
  have h42 : 0 ≤ Real.exp (Real.exp 1 * ξ) := (Real.exp_pos _).le
  have inner : Real.exp (π / 4) * Real.exp (-√(3 * π * k / 2)) +
          ξ * Real.exp (-√(3 * π * ((k : ℝ) - 1) / 2)) +
          (Real.exp ξ - 1 - ξ) * Real.exp (-√(9 * π * ((k : ℝ) - 2) / 4)) +
          Real.exp (Real.exp 1 * ξ) * Real.exp (-(k : ℝ) / 2) ≤
        (Real.exp (π / 4) + Real.exp ξ + Real.exp (Real.exp 1 * ξ)) * e := by
    have := mul_le_mul_of_nonneg_left t1 h41
    have := mul_le_mul_of_nonneg_left t2 hxi
    have := mul_le_mul_of_nonneg_left t3 hex
    have := mul_le_mul_of_nonneg_left t4 h42
    have : 0 ≤ e := (Real.exp_pos _).le
    nlinarith
  have := mul_le_mul_of_nonneg_left inner hB1.le
  have := mul_le_mul_of_nonneg_left t5 hB2.le
  unfold Cstar
  nlinarith

theorem main : Q07.stmt_4 → Q07.stmt_5 := by
  intro h4 k hk
  have h1 := (h4 k hk).trans (rhs_le k hk)
  refine ⟨h1, h1.trans ?_⟩
  exact mul_le_mul_of_nonneg_right Q05.proof_4.le (Real.exp_pos _).le

end Q07.PE

/-! ### Merged from work/Q07/PartF.lean -/
/-!
# Part F of Q07: g5:prop:main (`Q07.stmt_6`)

Plan (sec_g5_varpi1.tex, with simplifications that keep the constants below 71,000):
1. rotate `∫_0^∞ H E Q_k` to the ray `ωt` and substitute `t = ρ e^σ`;
2. pointwise and tail bounds for `f(σ) = H(s) E(s) (s/(1+s))^N`, `s = ωρe^σ`;
3. region I (`N ≤ N₁`): crude bound;
4. region II (`N ≥ N₁`): saddle point expansion with explicit remainder.
-/

namespace Q07.PF

open Real MeasureTheory Set Filter Q05

/-! ## 1. The rotation of the contour to the ray `ωt` -/

lemma ω_re : ω.re = √2 / 2 := by rw [ω_eq]; simp

lemma ω_im : ω.im = √2 / 2 := by rw [ω_eq]; simp

lemma norm_ω : ‖ω‖ = 1 := by
  rw [ω, Complex.norm_exp]; simp

lemma ωt_eq (t : ℝ) : ω * (t : ℂ) = ((t * (√2 / 2) : ℝ) : ℂ) + ((t * (√2 / 2) : ℝ) : ℂ) * Complex.I := by
  rw [ω_eq]; push_cast; ring

lemma ωt_re (t : ℝ) : (ω * (t : ℂ)).re = t * (√2 / 2) := by rw [ωt_eq]; simp

lemma ωt_im (t : ℝ) : (ω * (t : ℂ)).im = t * (√2 / 2) := by rw [ωt_eq]; simp

lemma norm_ωt (t : ℝ) : ‖ω * (t : ℂ)‖ = |t| := by
  rw [norm_mul, norm_ω, one_mul, Complex.norm_real, Real.norm_eq_abs]

lemma mem_ball_gen (u v X r : ℝ) (h : (u - X) ^ 2 + v ^ 2 < (X - r / 2) ^ 2) (hX : r / 2 < X) :
    ((u : ℂ) + (v : ℂ) * Complex.I) ∈ Metric.ball (X : ℂ) (X - r / 2) := by
  rw [Metric.mem_ball, dist_eq_norm, Complex.norm_eq_sqrt_sq_add_sq]
  rw [Real.sqrt_lt' (by linarith)]
  simpa using h

lemma ball_sub_gen (X r : ℝ) (hr : 0 < r) :
    Metric.ball (X : ℂ) (X - r / 2) ⊆ {s : ℂ | 0 < s.re} := by
  intro z hz
  rw [Metric.mem_ball, dist_eq_norm] at hz
  have := Complex.abs_re_le_norm (z - X)
  simp only [Complex.sub_re, Complex.ofReal_re] at this
  show 0 < z.re
  have := neg_abs_le (z.re - X)
  linarith

lemma sqrt2_pos : (0 : ℝ) < √2 := by positivity

lemma contour_sector (k : ℕ) (R r : ℝ) (hr : 0 < r) (hrR : r ≤ R) :
    (∫ x in r..R, fM k x) + ∫ v in (0 : ℝ)..R, fM k (R + v * Complex.I) * Complex.I =
      (∫ t in (r * √2)..(R * √2), fM k (ω * t) * ω) +
        ∫ v in (0 : ℝ)..r, fM k (r + v * Complex.I) * Complex.I := by
  set X : ℝ := 3 * R with hX
  have hR0 : 0 < R := lt_of_lt_of_le hr hrR
  have hXr : r / 2 < X := by rw [hX]; linarith
  have hdiff : DifferentiableOn ℂ (fM k) (Metric.ball (X : ℂ) (X - r / 2)) :=
    (fM_diff k).mono (ball_sub_gen X r hr)
  obtain ⟨F, hF⟩ := hdiff.isExactOn_ball
  have hcont : ContinuousOn (fM k) (Metric.ball (X : ℂ) (X - r / 2)) := hdiff.continuousOn
  have hmR : ∀ x ∈ uIcc r R, (x : ℂ) ∈ Metric.ball (X : ℂ) (X - r / 2) := by
    intro x hx
    rw [uIcc_of_le hrR] at hx
    have := mem_ball_gen x 0 X r (by nlinarith [hx.1, hx.2]) hXr
    simpa using this
  have hmV : ∀ v ∈ uIcc (0 : ℝ) R,
      ((R : ℂ) + (v : ℂ) * Complex.I) ∈ Metric.ball (X : ℂ) (X - r / 2) := by
    intro v hv
    rw [uIcc_of_le hR0.le] at hv
    exact mem_ball_gen R v X r (by nlinarith [hv.1, hv.2]) hXr
  have hmv : ∀ v ∈ uIcc (0 : ℝ) r,
      ((r : ℂ) + (v : ℂ) * Complex.I) ∈ Metric.ball (X : ℂ) (X - r / 2) := by
    intro v hv
    rw [uIcc_of_le hr.le] at hv
    exact mem_ball_gen r v X r (by nlinarith [hv.1, hv.2]) hXr
  have h22 := sqrt2_sq
  have hmL : ∀ t ∈ uIcc (r * √2) (R * √2), ω * (t : ℂ) ∈ Metric.ball (X : ℂ) (X - r / 2) := by
    intro t ht
    rw [uIcc_of_le (mul_le_mul_of_nonneg_right hrR sqrt2_pos.le)] at ht
    rw [ωt_eq]
    set u := t * (√2 / 2)
    have hu0 : r ≤ u := by
      have := mul_le_mul_of_nonneg_right ht.1 (by positivity : (0 : ℝ) ≤ √2 / 2)
      rwa [sq2half] at this
    have hu1 : u ≤ R := by
      have := mul_le_mul_of_nonneg_right ht.2 (by positivity : (0 : ℝ) ≤ √2 / 2)
      rwa [sq2half] at this
    apply mem_ball_gen _ _ X r _ hXr
    nlinarith
  have hA : ∫ x in r..R, fM k x = F R - F r := by
    apply intervalIntegral.integral_eq_sub_of_hasDerivAt
    · intro x hx; exact (hF _ (hmR x hx)).comp_ofReal
    · apply ContinuousOn.intervalIntegrable
      exact hcont.comp Complex.continuous_ofReal.continuousOn (fun x hx => hmR x hx)
  have hVgen : ∀ c : ℝ, 0 < c → (∀ v ∈ uIcc (0 : ℝ) c,
      ((c : ℂ) + (v : ℂ) * Complex.I) ∈ Metric.ball (X : ℂ) (X - r / 2)) →
      ∫ v in (0 : ℝ)..c, fM k (c + v * Complex.I) * Complex.I =
        F (c + ((c : ℝ) : ℂ) * Complex.I) - F (c + ((0 : ℝ) : ℂ) * Complex.I) := by
    intro c hc hm
    apply intervalIntegral.integral_eq_sub_of_hasDerivAt (f := fun v : ℝ => F (c + v * Complex.I))
    · intro v hv
      have h1 : HasDerivAt (fun z : ℂ => (c : ℂ) + z * Complex.I) Complex.I (v : ℂ) := by
        simpa using ((hasDerivAt_id (v : ℂ)).mul_const Complex.I).const_add (c : ℂ)
      exact ((hF _ (hm v hv)).comp (v : ℂ) h1).comp_ofReal
    · apply ContinuousOn.intervalIntegrable
      refine ContinuousOn.mul ?_ continuousOn_const
      exact hcont.comp (by fun_prop) (fun v hv => hm v hv)
  have hB : ∫ t in (r * √2)..(R * √2), fM k (ω * t) * ω =
      F (ω * ((R * √2 : ℝ) : ℂ)) - F (ω * ((r * √2 : ℝ) : ℂ)) := by
    apply intervalIntegral.integral_eq_sub_of_hasDerivAt (f := fun t : ℝ => F (ω * t))
    · intro t ht
      have h1 : HasDerivAt (fun z : ℂ => ω * z) ω (t : ℂ) := by
        simpa using (hasDerivAt_id (t : ℂ)).const_mul ω
      exact ((hF _ (hmL t ht)).comp (t : ℂ) h1).comp_ofReal
    · apply ContinuousOn.intervalIntegrable
      refine ContinuousOn.mul ?_ continuousOn_const
      exact hcont.comp (by fun_prop) (fun t ht => hmL t ht)
  have hpt : ∀ c : ℝ, ω * ((c * √2 : ℝ) : ℂ) = (c : ℂ) + ((c : ℝ) : ℂ) * Complex.I := by
    intro c; rw [ωt_eq, sq2half]
  rw [hA, hVgen R hR0 hmV, hVgen r hr hmv, hB, hpt, hpt]
  push_cast
  ring_nf

lemma HQ_ray (k : ℕ) (hk : 1 ≤ k) (t : ℝ) (ht : 0 < t) :
    ‖H (ω * (t : ℂ))‖ * ‖Q k (ω * (t : ℂ))‖ ≤ √2 := by
  have hre : (ω * (t : ℂ)).re = t * (√2 / 2) := ωt_re t
  have hre0 : 0 < (ω * (t : ℂ)).re := by rw [hre]; positivity
  have hH := H_le_one_div _ hre0
  rw [hre] at hH
  have hH' : ‖H (ω * (t : ℂ))‖ ≤ √2 / t := by
    refine hH.trans (le_of_eq ?_)
    field_simp
    rw [Real.sq_sqrt (by norm_num)]
  have h1 : 1 ≤ ‖1 + ω * (t : ℂ)‖ := by
    have := one_le_norm_one_add (ω * (t : ℂ)); linarith
  have hQ1 : ‖Q k (ω * (t : ℂ))‖ ≤ 1 := by
    refine (Q_le k _ hre0.le).trans ?_
    rw [div_le_one (by linarith)]; exact h1
  have hst : 0 < √2 / t := by positivity
  rcases le_or_gt t 1 with h | h
  · have hQt : ‖Q k (ω * (t : ℂ))‖ ≤ t := by
      rw [Q_le', norm_ωt, abs_of_pos ht]
      have h2 : 1 ≤ ‖1 + ω * (t : ℂ)‖ ^ (k + 1) := one_le_pow₀ h1
      have h3 : t ^ k ≤ t := by
        calc t ^ k ≤ t ^ 1 := pow_le_pow_of_le_one ht.le h hk
          _ = t := pow_one t
      calc t ^ k / ‖1 + ω * (t : ℂ)‖ ^ (k + 1) ≤ t ^ k / 1 :=
            div_le_div_of_nonneg_left (by positivity) one_pos h2
        _ ≤ t := by rw [div_one]; exact h3
    calc ‖H (ω * (t : ℂ))‖ * ‖Q k (ω * (t : ℂ))‖ ≤ (√2 / t) * t :=
          mul_le_mul hH' hQt (norm_nonneg _) hst.le
      _ = √2 := by field_simp
  · calc ‖H (ω * (t : ℂ))‖ * ‖Q k (ω * (t : ℂ))‖ ≤ (√2 / t) * 1 :=
          mul_le_mul hH' hQ1 (norm_nonneg _) hst.le
      _ ≤ √2 := by
          rw [mul_one, div_le_iff₀ ht]; nlinarith [sqrt2_pos]

/-- `c_E = π√2/4`, the decay rate of `E` on the ray. -/
noncomputable def cE : ℝ := π * √2 / 4

lemma cE_pos : 0 < cE := by unfold cE; positivity

lemma E_ray (t : ℝ) : ‖E (ω * (t : ℂ))‖ = Real.exp (-cE * t) := by
  rw [E_norm, ωt_im, cE]; congr 1; ring

lemma fM_ray_le (k : ℕ) (hk : 1 ≤ k) (t : ℝ) (ht : 0 < t) :
    ‖fM k (ω * (t : ℂ)) * ω‖ ≤ √2 * Real.exp (-cE * t) := by
  rw [norm_mul, norm_ω, mul_one, fM, norm_mul, norm_mul, E_ray]
  have := HQ_ray k hk t ht
  have he := Real.exp_pos (-cE * t)
  calc ‖H (ω * (t : ℂ))‖ * Real.exp (-cE * t) * ‖Q k (ω * (t : ℂ))‖
      = (‖H (ω * (t : ℂ))‖ * ‖Q k (ω * (t : ℂ))‖) * Real.exp (-cE * t) := by ring
    _ ≤ √2 * Real.exp (-cE * t) := mul_le_mul_of_nonneg_right this he.le

lemma fM_ray_cont (k : ℕ) : ContinuousOn (fun t : ℝ => fM k (ω * (t : ℂ)) * ω) (Ioi 0) := by
  refine ContinuousOn.mul ?_ continuousOn_const
  refine (fM_cont k).comp (by fun_prop) (fun t ht => ?_)
  show 0 < (ω * (t : ℂ)).re
  rw [ωt_re]; have : 0 < t := ht; positivity

lemma fM_ray_int (k : ℕ) (hk : 1 ≤ k) :
    IntegrableOn (fun t : ℝ => fM k (ω * (t : ℂ)) * ω) (Ioi 0) := by
  have hg : IntegrableOn (fun t : ℝ => √2 * Real.exp (-cE * t)) (Ioi 0) :=
    (exp_neg_integrableOn_Ioi 0 cE_pos).const_mul _
  refine Integrable.mono' hg ((fM_ray_cont k).aestronglyMeasurable measurableSet_Ioi) ?_
  refine (ae_restrict_iff' measurableSet_Ioi).2 (Filter.Eventually.of_forall fun t ht => ?_)
  exact fM_ray_le k hk t ht

lemma fM_real_le_one (k : ℕ) (hk : 1 ≤ k) (x : ℝ) (hx0 : 0 < x) (hx1 : x ≤ 1) :
    ‖fM k (x : ℂ)‖ ≤ 1 := by
  have hH := H_le_one_div (x : ℂ) (by simpa using hx0)
  simp only [Complex.ofReal_re] at hH
  have hE : ‖E (x : ℂ)‖ = 1 := by rw [E_norm]; simp
  have hQ : ‖Q k (x : ℂ)‖ ≤ x := by
    rw [Q_le', Complex.norm_real, Real.norm_eq_abs, abs_of_pos hx0]
    have h1 : 1 ≤ ‖1 + (x : ℂ)‖ ^ (k + 1) := by
      apply one_le_pow₀; have := one_le_norm_one_add (x : ℂ); simp at this; linarith
    have h2 : x ^ k ≤ x := by
      calc x ^ k ≤ x ^ 1 := pow_le_pow_of_le_one hx0.le hx1 hk
        _ = x := pow_one x
    calc x ^ k / ‖1 + (x : ℂ)‖ ^ (k + 1) ≤ x ^ k / 1 :=
          div_le_div_of_nonneg_left (by positivity) one_pos h1
      _ ≤ x := by rw [div_one]; exact h2
  rw [fM, norm_mul, norm_mul, hE, mul_one]
  calc ‖H x‖ * ‖Q k x‖ ≤ 1 / x * x := mul_le_mul hH hQ (norm_nonneg _) (by positivity)
    _ = 1 := by field_simp

lemma vert_big (k : ℕ) (R : ℝ) (hR : 1 ≤ R) :
    ‖∫ v in (0 : ℝ)..R, fM k (R + v * Complex.I) * Complex.I‖ ≤ R⁻¹ := by
  have hb : ∀ v ∈ Set.uIoc (0 : ℝ) R, ‖fM k (R + v * Complex.I) * Complex.I‖ ≤ 1 / R ^ 2 := by
    intro v hv
    rw [uIoc_of_le (by linarith)] at hv
    have hv0 : 0 ≤ v := le_of_lt hv.1
    set s : ℂ := (R : ℂ) + (v : ℂ) * Complex.I
    have hre : s.re = R := by simp [s]
    have him : s.im = v := by simp [s]
    have hH := H_le_one_div s (by rw [hre]; linarith)
    rw [hre] at hH
    have hE : ‖E s‖ ≤ 1 := by
      rw [E_norm, him]; apply Real.exp_le_one_iff.2
      have := Real.pi_pos; nlinarith
    have hQ : ‖Q k s‖ ≤ 1 / R := by
      refine (Q_le k s (by rw [hre]; linarith)).trans ?_
      have := one_le_norm_one_add s
      rw [hre] at this
      exact div_le_div_of_nonneg_left zero_le_one (by linarith) (by linarith)
    rw [norm_mul, Complex.norm_I, mul_one, fM, norm_mul, norm_mul]
    calc ‖H s‖ * ‖E s‖ * ‖Q k s‖ ≤ (1 / R) * 1 * (1 / R) := by gcongr
      _ = 1 / R ^ 2 := by ring
  have := intervalIntegral.norm_integral_le_of_norm_le_const hb
  refine this.trans ?_
  rw [sub_zero, abs_of_nonneg (by linarith)]
  have hR0 : 0 < R := by linarith
  rw [div_mul_eq_mul_div, one_mul, div_le_iff₀ (by positivity)]
  have : R⁻¹ * R ^ 2 = R := by field_simp
  rw [this]

lemma vert_small (k : ℕ) (hk : 1 ≤ k) (r : ℝ) (hr : 0 < r) (hr1 : r ≤ 1 / 2) :
    ‖∫ v in (0 : ℝ)..r, fM k (r + v * Complex.I) * Complex.I‖ ≤ √2 * r := by
  have hb : ∀ v ∈ Set.uIoc (0 : ℝ) r, ‖fM k (r + v * Complex.I) * Complex.I‖ ≤ √2 := by
    intro v hv
    rw [uIoc_of_le hr.le] at hv
    have hv0 : 0 ≤ v := le_of_lt hv.1
    set s : ℂ := (r : ℂ) + (v : ℂ) * Complex.I
    have hre : s.re = r := by simp [s]
    have him : s.im = v := by simp [s]
    have hH := H_le_one_div s (by rw [hre]; linarith)
    rw [hre] at hH
    have hE : ‖E s‖ ≤ 1 := by
      rw [E_norm, him]; apply Real.exp_le_one_iff.2
      have := Real.pi_pos; nlinarith
    have hns : ‖s‖ ≤ √2 * r := by
      rw [Complex.norm_eq_sqrt_sq_add_sq, hre, him]
      rw [Real.sqrt_le_left (by positivity)]
      rw [mul_pow, Real.sq_sqrt (by norm_num)]
      nlinarith [hv.2]
    have hs1 : ‖s‖ ≤ 1 := by
      refine hns.trans ?_
      have : √2 ≤ 2 := by rw [Real.sqrt_le_left (by norm_num)]; norm_num
      nlinarith
    have hQ : ‖Q k s‖ ≤ √2 * r := by
      rw [Q_le']
      have h1 : 1 ≤ ‖1 + s‖ ^ (k + 1) := by
        apply one_le_pow₀; have := one_le_norm_one_add s; rw [hre] at this; linarith
      have h2 : ‖s‖ ^ k ≤ ‖s‖ := by
        calc ‖s‖ ^ k ≤ ‖s‖ ^ 1 := pow_le_pow_of_le_one (norm_nonneg _) hs1 hk
          _ = ‖s‖ := pow_one _
      calc ‖s‖ ^ k / ‖1 + s‖ ^ (k + 1) ≤ ‖s‖ ^ k / 1 :=
            div_le_div_of_nonneg_left (by positivity) one_pos h1
        _ ≤ √2 * r := by rw [div_one]; exact h2.trans hns
    rw [norm_mul, Complex.norm_I, mul_one, fM, norm_mul, norm_mul]
    calc ‖H s‖ * ‖E s‖ * ‖Q k s‖ ≤ (1 / r) * 1 * (√2 * r) := by gcongr
      _ = √2 := by field_simp
  have := intervalIntegral.norm_integral_le_of_norm_le_const hb
  refine this.trans ?_
  rw [sub_zero, abs_of_nonneg hr.le, mul_comm]

/-- The rotation: `∫_0^∞ H E Q_k = ∫_0^∞ H(ωt) E(ωt) Q_k(ωt) ω dt`. -/
lemma rot (k : ℕ) (hk : 1 ≤ k) :
    ∫ x in Ioi (0 : ℝ), fM k x = ∫ t in Ioi (0 : ℝ), fM k (ω * (t : ℂ)) * ω := by
  have hI1 : IntegrableOn (fun x : ℝ => fM k x) (Ioi 0) := fM_real_int k hk
  have hI2 := fM_ray_int k hk
  -- the four pieces as functions of R, with r = 1/(2R)
  have tA := intervalIntegral_tendsto_integral_Ioi (0 : ℝ) hI1 tendsto_id
  have hT : Tendsto (fun R : ℝ => R * √2) atTop atTop :=
    tendsto_id.atTop_mul_const sqrt2_pos
  have tB := intervalIntegral_tendsto_integral_Ioi (0 : ℝ) hI2 hT
  have hr0 : Tendsto (fun R : ℝ => 1 / (2 * R)) atTop (nhds 0) := by
    exact (tendsto_inv_atTop_zero.comp (tendsto_id.const_mul_atTop two_pos)).congr
      (fun R => by simp [Function.comp, one_div])
  -- small pieces tend to 0
  have tA0 : Tendsto (fun R : ℝ => ∫ x in (0 : ℝ)..(1 / (2 * R)), fM k x) atTop (nhds 0) := by
    refine squeeze_zero_norm' ?_ hr0
    filter_upwards [eventually_ge_atTop (1 : ℝ)] with R hR
    have hp : 0 < 1 / (2 * R) := by positivity
    have hb : ∀ x ∈ Set.uIoc (0 : ℝ) (1 / (2 * R)), ‖fM k x‖ ≤ 1 := by
      intro x hx
      rw [uIoc_of_le hp.le] at hx
      refine fM_real_le_one k hk x hx.1 (hx.2.trans ?_)
      rw [div_le_one (by positivity)]; linarith
    have := intervalIntegral.norm_integral_le_of_norm_le_const hb
    rw [sub_zero, abs_of_pos hp, one_mul] at this
    exact this
  have tB0 : Tendsto (fun R : ℝ => ∫ t in (0 : ℝ)..(1 / (2 * R) * √2), fM k (ω * (t : ℂ)) * ω)
      atTop (nhds 0) := by
    have hr0' : Tendsto (fun R : ℝ => √2 * (1 / (2 * R) * √2)) atTop (nhds 0) := by
      simpa using (hr0.mul_const √2).const_mul √2
    refine squeeze_zero_norm' ?_ hr0'
    filter_upwards [eventually_ge_atTop (1 : ℝ)] with R hR
    have hp : 0 < 1 / (2 * R) * √2 := by positivity
    have hb : ∀ t ∈ Set.uIoc (0 : ℝ) (1 / (2 * R) * √2), ‖fM k (ω * (t : ℂ)) * ω‖ ≤ √2 := by
      intro t ht
      rw [uIoc_of_le hp.le] at ht
      refine (fM_ray_le k hk t ht.1).trans ?_
      have : Real.exp (-cE * t) ≤ 1 := by
        rw [Real.exp_le_one_iff]; have := cE_pos; have := ht.1; nlinarith
      nlinarith [sqrt2_pos]
    have := intervalIntegral.norm_integral_le_of_norm_le_const hb
    rw [sub_zero, abs_of_pos hp] at this
    exact this
  have tV : Tendsto (fun R : ℝ => ∫ v in (0 : ℝ)..R, fM k (R + v * Complex.I) * Complex.I)
      atTop (nhds 0) := by
    refine squeeze_zero_norm' ?_ tendsto_inv_atTop_zero
    filter_upwards [eventually_ge_atTop (1 : ℝ)] with R hR
    exact vert_big k R hR
  have tv : Tendsto (fun R : ℝ => ∫ v in (0 : ℝ)..(1 / (2 * R)),
      fM k ((1 / (2 * R) : ℝ) + v * Complex.I) * Complex.I) atTop (nhds 0) := by
    have hr0' : Tendsto (fun R : ℝ => √2 * (1 / (2 * R))) atTop (nhds 0) := by
      simpa using hr0.const_mul √2
    refine squeeze_zero_norm' ?_ hr0'
    filter_upwards [eventually_ge_atTop (1 : ℝ)] with R hR
    refine vert_small k hk _ (by positivity) ?_
    rw [div_le_div_iff₀ (by positivity) (by norm_num)]; linarith
  -- assemble
  have hL : Tendsto (fun R : ℝ => ((∫ x in (0 : ℝ)..R, fM k x) - ∫ x in (0 : ℝ)..(1 / (2 * R)), fM k x) +
      ∫ v in (0 : ℝ)..R, fM k (R + v * Complex.I) * Complex.I) atTop
      (nhds ((∫ x in Ioi (0 : ℝ), fM k x) - 0 + 0)) := (tA.sub tA0).add tV
  have hRt : Tendsto (fun R : ℝ => ((∫ t in (0 : ℝ)..(R * √2), fM k (ω * (t : ℂ)) * ω) -
      ∫ t in (0 : ℝ)..(1 / (2 * R) * √2), fM k (ω * (t : ℂ)) * ω) +
      ∫ v in (0 : ℝ)..(1 / (2 * R)), fM k ((1 / (2 * R) : ℝ) + v * Complex.I) * Complex.I) atTop
      (nhds ((∫ t in Ioi (0 : ℝ), fM k (ω * (t : ℂ)) * ω) - 0 + 0)) := (tB.sub tB0).add tv
  rw [sub_zero, add_zero] at hL hRt
  refine tendsto_nhds_unique hL (hRt.congr' ?_)
  filter_upwards [eventually_ge_atTop (1 : ℝ)] with R hR
  have hp : 0 < 1 / (2 * R) := by positivity
  have hpR : 1 / (2 * R) ≤ R := by
    rw [div_le_iff₀ (by positivity)]; nlinarith
  have hc := contour_sector k R (1 / (2 * R)) hp hpR
  have ii : ∀ {g : ℝ → ℂ}, IntegrableOn g (Ioi 0) → ∀ b : ℝ, 0 ≤ b →
      IntervalIntegrable g volume 0 b := fun hg b hb =>
    (intervalIntegrable_iff_integrableOn_Ioc_of_le hb).2 (hg.mono_set Ioc_subset_Ioi_self)
  have e1 := intervalIntegral.integral_interval_sub_left (ii hI1 R (by linarith)) (ii hI1 _ hp.le)
  have e2 := intervalIntegral.integral_interval_sub_left (ii hI2 (R * √2) (by positivity))
    (ii hI2 (1 / (2 * R) * √2) (by positivity))
  rw [e1, e2, hc]

/-! ## 2. The substitution `t = ρ e^σ` -/

lemma cov (g : ℝ → ℂ) (hg : IntegrableOn g (Ioi 0)) (ρ : ℝ) (hρ : 0 < ρ) :
    Integrable (fun σ : ℝ => ((ρ * Real.exp σ : ℝ) : ℂ) * g (ρ * Real.exp σ)) ∧
      ∫ σ : ℝ, ((ρ * Real.exp σ : ℝ) : ℂ) * g (ρ * Real.exp σ) = ∫ t in Ioi (0 : ℝ), g t := by
  set G : ℝ → ℂ := fun u => (ρ : ℂ) * g (ρ * u)
  have hG : IntegrableOn G (Ioi 0) := by
    have h1 : IntegrableOn (fun u : ℝ => g (ρ * u)) (Ioi 0) :=
      (integrableOn_Ioi_comp_mul_left_iff g 0 hρ).2 (by simpa using hg)
    exact h1.const_mul (ρ : ℂ)
  have heq : (fun σ : ℝ => ((ρ * Real.exp σ : ℝ) : ℂ) * g (ρ * Real.exp σ)) =
      fun σ : ℝ => Real.exp σ • G (Real.exp σ) := by
    funext σ; simp only [G, Complex.real_smul]; push_cast; ring
  rw [heq]
  refine ⟨(integrable_comp_exp G).2 hG, ?_⟩
  rw [integral_comp_exp G]
  simp only [G]
  rw [integral_const_mul, integral_comp_mul_left_Ioi g 0 hρ, mul_zero, Complex.real_smul]
  have : (ρ : ℂ) ≠ 0 := by exact_mod_cast hρ.ne'
  push_cast
  field_simp

/-- `ρ = √(2N/π)`. -/
noncomputable def ρv (N : ℝ) : ℝ := √(2 * N / π)

/-- The point `s = ωρe^σ`. -/
noncomputable def sv (N σ : ℝ) : ℂ := ω * ((ρv N * Real.exp σ : ℝ) : ℂ)

/-- `f(σ) = H(s) E(s) (s/(1+s))^N`, `s = ωρe^σ`. -/
noncomputable def fσ (n : ℕ) (σ : ℝ) : ℂ :=
  H (sv n σ) * E (sv n σ) * (sv n σ / (1 + sv n σ)) ^ n

lemma ρv_pos (N : ℝ) (hN : 0 < N) : 0 < ρv N := by unfold ρv; positivity

lemma fM_mul (k : ℕ) (s : ℂ) (hs : 1 + s ≠ 0) : fM k s * s = H s * E s * (s / (1 + s)) ^ (k + 1) := by
  rw [fM, Q, div_pow, pow_succ s k]
  field_simp

/-- g5:lem:sigma (a): `∫_0^∞ H E Q_k = ∫_ℝ f(σ) dσ`, absolutely convergent. -/
lemma J_eq (k : ℕ) (hk : 1 ≤ k) :
    Integrable (fσ (k + 1)) ∧ ∫ x in Ioi (0 : ℝ), fM k x = ∫ σ : ℝ, fσ (k + 1) σ := by
  have hρ := ρv_pos ((k + 1 : ℕ) : ℝ) (by positivity)
  obtain ⟨h1, h2⟩ := cov (fun t : ℝ => fM k (ω * (t : ℂ)) * ω) (fM_ray_int k hk) _ hρ
  have heq : (fun σ : ℝ => (((ρv ((k + 1 : ℕ) : ℝ)) * Real.exp σ : ℝ) : ℂ) *
      (fM k (ω * (((ρv ((k + 1 : ℕ) : ℝ)) * Real.exp σ : ℝ) : ℂ)) * ω)) = fσ (k + 1) := by
    funext σ
    have hs : 1 + sv (k + 1 : ℕ) σ ≠ 0 := by
      apply one_add_ne
      rw [sv, ωt_re]; positivity
    rw [fσ, ← fM_mul k _ hs, sv]; push_cast; ring
  rw [heq] at h1 h2
  exact ⟨h1, by rw [rot k hk, h2]⟩

/-! ## 3. Pointwise bounds for `f` -/

/-- `y = √(π/(2N))`. -/
noncomputable def yv (N : ℝ) : ℝ := √(π / (2 * N))

/-- `θ = √(πN)`. -/
noncomputable def θv (N : ℝ) : ℝ := √(π * N)

/-- `γ₀ = √2 π/12`. -/
noncomputable def γ0 : ℝ := √2 * π / 12

lemma yv_pos (N : ℝ) (hN : 0 < N) : 0 < yv N := by unfold yv; positivity

lemma θv_pos (N : ℝ) (hN : 0 < N) : 0 < θv N := by unfold θv; positivity

lemma ρy (N : ℝ) (hN : 0 < N) : ρv N * yv N = 1 := by
  unfold ρv yv
  rw [← Real.sqrt_mul (by positivity)]
  rw [show 2 * N / π * (π / (2 * N)) = 1 by field_simp]
  simp

lemma θ_sq (N : ℝ) (hN : 0 < N) : θv N ^ 2 = π * N := by
  unfold θv; rw [Real.sq_sqrt (by positivity)]

lemma y_sq (N : ℝ) (hN : 0 < N) : yv N ^ 2 = π / (2 * N) := by
  unfold yv; rw [Real.sq_sqrt (by positivity)]

lemma ρ_sq (N : ℝ) (hN : 0 < N) : ρv N ^ 2 = 2 * N / π := by
  unfold ρv; rw [Real.sq_sqrt (by positivity)]

lemma Ny2 (N : ℝ) (hN : 0 < N) : N * yv N ^ 2 = π / 2 := by
  rw [y_sq N hN]; field_simp

lemma Nyθ (N : ℝ) (hN : 0 < N) : N * yv N * √2 = θv N := by
  have h1 : 0 ≤ N * yv N * √2 := by have := yv_pos N hN; positivity
  have h2 : 0 ≤ θv N := (θv_pos N hN).le
  have : (N * yv N * √2) ^ 2 = θv N ^ 2 := by
    rw [mul_pow, mul_pow, y_sq N hN, θ_sq N hN, Real.sq_sqrt (by norm_num)]
    field_simp
  exact (pow_left_inj₀ h1 h2 two_ne_zero).1 this

lemma πρθ (N : ℝ) (hN : 0 < N) : π * ρv N * √2 / 4 = θv N / 2 := by
  have h1 : 0 ≤ π * ρv N * √2 := by have := ρv_pos N hN; positivity
  have h2 : 0 ≤ θv N := (θv_pos N hN).le
  have : (π * ρv N * √2) ^ 2 = (2 * θv N) ^ 2 := by
    rw [mul_pow, mul_pow, mul_pow, ρ_sq N hN, θ_sq N hN, Real.sq_sqrt (by norm_num)]
    field_simp
  have := (pow_left_inj₀ h1 (by positivity) two_ne_zero).1 this
  linarith

lemma sv_re (N σ : ℝ) : (sv N σ).re = ρv N * Real.exp σ * (√2 / 2) := by rw [sv, ωt_re]

lemma sv_im (N σ : ℝ) : (sv N σ).im = ρv N * Real.exp σ * (√2 / 2) := by rw [sv, ωt_im]

lemma sv_re_pos (N σ : ℝ) (hN : 0 < N) : 0 < (sv N σ).re := by
  rw [sv_re]; have := ρv_pos N hN; positivity

lemma norm_sv (N σ : ℝ) (hN : 0 < N) : ‖sv N σ‖ = ρv N * Real.exp σ := by
  rw [sv, norm_ωt, abs_of_pos (by have := ρv_pos N hN; positivity)]

lemma E_sv (N σ : ℝ) (hN : 0 < N) : ‖E (sv N σ)‖ = Real.exp (-(θv N * Real.exp σ / 2)) := by
  rw [E_norm, sv_im]
  congr 1
  have := πρθ N hN
  have : π * (ρv N * Real.exp σ * (√2 / 2)) / 2 = (π * ρv N * √2 / 4) * Real.exp σ := by ring
  rw [this, πρθ N hN]; ring

lemma H_sv (N σ : ℝ) (hN : 0 < N) : ‖H (sv N σ)‖ ≤ √2 * yv N * Real.exp (-σ) := by
  refine (H_le_one_div _ (sv_re_pos N σ hN)).trans (le_of_eq ?_)
  rw [sv_re]
  have h1 := ρy N hN
  have hρ := ρv_pos N hN
  have he : Real.exp σ * Real.exp (-σ) = 1 := by rw [← Real.exp_add]; simp
  have hy : yv N = 1 / ρv N := by field_simp; linarith
  rw [hy]
  field_simp
  rw [Real.sq_sqrt (by norm_num)]
  nlinarith [he]

/-- `u = √2 y e^{-σ}`. -/
noncomputable def uv (N σ : ℝ) : ℝ := √2 * yv N * Real.exp (-σ)

/-- `p = 1 + u + u²/2 = |1 + 1/s|²`. -/
noncomputable def pv (N σ : ℝ) : ℝ := 1 + uv N σ + uv N σ ^ 2 / 2

lemma uv_nonneg (N σ : ℝ) (hN : 0 < N) : 0 ≤ uv N σ := by
  unfold uv; have := yv_pos N hN; positivity

lemma pv_pos (N σ : ℝ) (hN : 0 < N) : 0 < pv N σ := by
  unfold pv; have := uv_nonneg N σ hN; positivity

lemma ratio_sq (N σ : ℝ) (hN : 0 < N) :
    ‖sv N σ / (1 + sv N σ)‖ ^ 2 * pv N σ = 1 := by
  have hs : 1 + sv N σ ≠ 0 := one_add_ne _ (sv_re_pos N σ hN).le
  have hn : 0 < ‖1 + sv N σ‖ := norm_pos_iff.2 hs
  rw [norm_div, div_pow, norm_sv N σ hN]
  have h1 : ‖1 + sv N σ‖ ^ 2 = (1 + ρv N * Real.exp σ * (√2 / 2)) ^ 2 +
      (ρv N * Real.exp σ * (√2 / 2)) ^ 2 := by
    rw [Complex.sq_norm, Complex.normSq_apply]
    simp only [Complex.add_re, Complex.one_re, Complex.add_im, Complex.one_im, zero_add, sv_re, sv_im]
    ring
  rw [h1, div_mul_eq_mul_div, div_eq_one_iff_eq (by rw [← h1]; positivity)]
  unfold pv uv
  have hρ := ρy N hN
  have he : Real.exp σ * Real.exp (-σ) = 1 := by rw [← Real.exp_add]; simp
  have h22 : √2 ^ 2 = 2 := Real.sq_sqrt (by norm_num)
  set A := ρv N * Real.exp σ
  set B := yv N * Real.exp (-σ)
  have hAB : A * B = 1 := by
    simp only [A, B]
    calc ρv N * Real.exp σ * (yv N * Real.exp (-σ)) = (ρv N * yv N) * (Real.exp σ * Real.exp (-σ)) := by ring
      _ = 1 := by rw [hρ, he, one_mul]
  have e1 : √2 * yv N * Real.exp (-σ) = √2 * B := by simp only [B]; ring
  rw [e1]
  have : A ^ 2 * (1 + √2 * B + (√2 * B) ^ 2 / 2) = A ^ 2 + √2 * A * (A * B) + (A * B) ^ 2 * (√2 ^ 2 / 2) := by ring
  rw [this, hAB, h22]
  have : (1 + A * (√2 / 2)) ^ 2 + (A * (√2 / 2)) ^ 2 = 1 + √2 * A + A ^ 2 * (√2 ^ 2 / 2) := by ring
  rw [this, h22]; ring

lemma norm_ratio_pow (n : ℕ) (N σ : ℝ) (hN : 0 < N) :
    ‖(sv N σ / (1 + sv N σ)) ^ n‖ = Real.exp (-(n / 2 * Real.log (pv N σ))) := by
  have hp := pv_pos N σ hN
  have h := ratio_sq N σ hN
  have hw : ‖sv N σ / (1 + sv N σ)‖ = Real.exp (-(Real.log (pv N σ) / 2)) := by
    have h0 : 0 ≤ ‖sv N σ / (1 + sv N σ)‖ := norm_nonneg _
    have h1 : 0 ≤ Real.exp (-(Real.log (pv N σ) / 2)) := (Real.exp_pos _).le
    have : ‖sv N σ / (1 + sv N σ)‖ ^ 2 = Real.exp (-(Real.log (pv N σ) / 2)) ^ 2 := by
      rw [← Real.exp_nat_mul]
      push_cast
      rw [show (2 : ℝ) * -(Real.log (pv N σ) / 2) = -Real.log (pv N σ) by ring, Real.exp_neg,
        Real.exp_log hp]
      field_simp; linarith
    exact (pow_left_inj₀ h0 h1 two_ne_zero).1 this
  rw [norm_pow, hw, ← Real.exp_nat_mul]
  congr 1; ring

lemma norm_fσ (n : ℕ) (hn : 0 < n) (σ : ℝ) :
    ‖fσ n σ‖ = ‖H (sv n σ)‖ *
      Real.exp (-(θv n * Real.exp σ / 2 + n / 2 * Real.log (pv n σ))) := by
  have hN : (0 : ℝ) < n := by exact_mod_cast hn
  rw [fσ, norm_mul, norm_mul, E_sv n σ hN, norm_ratio_pow n n σ hN, mul_assoc, ← Real.exp_add]
  congr 2; ring

lemma log_cubic (u : ℝ) (hu : 0 ≤ u) : u - u ^ 3 / 6 ≤ Real.log (1 + u + u ^ 2 / 2) := by
  let D : ℝ → ℝ := fun u => Real.log (1 + u + u ^ 2 / 2) - u + u ^ 3 / 6
  have hpos : ∀ x : ℝ, 0 < 1 + x + x ^ 2 / 2 := fun x => by nlinarith [sq_nonneg (x + 1)]
  have hd : ∀ x : ℝ, HasDerivAt D ((1 + x) / (1 + x + x ^ 2 / 2) - 1 + x ^ 2 / 2) x := by
    intro x
    have h1 : HasDerivAt (fun x : ℝ => 1 + x + x ^ 2 / 2) (1 + x) x := by
      have h := ((hasDerivAt_id' x).const_add 1).add ((hasDerivAt_pow 2 x).div_const 2)
      convert h using 1; norm_num
    have h2 := h1.log (hpos x).ne'
    have h3 := (h2.sub (hasDerivAt_id' x)).add ((hasDerivAt_pow 3 x).div_const 6)
    convert h3 using 1
    norm_num; ring
  have hmono : MonotoneOn D (Ici 0) := by
    apply monotoneOn_of_deriv_nonneg (convex_Ici 0)
    · exact fun x _ => (hd x).continuousAt.continuousWithinAt
    · intro x _; exact (hd x).differentiableAt.differentiableWithinAt
    · intro x hx
      rw [interior_Ici] at hx
      rw [(hd x).deriv]
      have hx0 : 0 < x := hx
      have hp := hpos x
      have : (1 + x) / (1 + x + x ^ 2 / 2) - 1 + x ^ 2 / 2 =
          x ^ 2 / 2 * (x + x ^ 2 / 2) / (1 + x + x ^ 2 / 2) := by
        field_simp; ring
      rw [this]; positivity
  have := hmono (by simp : (0 : ℝ) ∈ Ici 0) hu hu
  simp only [D] at this
  norm_num at this
  linarith

lemma exp_neg_pow3 (σ : ℝ) : Real.exp (-σ) ^ 3 = Real.exp (-(3 * σ)) := by
  rw [← Real.exp_nat_mul]; push_cast; ring_nf

lemma half_n_u (n : ℕ) (hn : 0 < n) (σ : ℝ) :
    (n : ℝ) / 2 * (uv n σ - uv n σ ^ 3 / 6) =
      θv n * Real.exp (-σ) / 2 - γ0 * yv n * Real.exp (-(3 * σ)) := by
  have hN : (0 : ℝ) < n := by exact_mod_cast hn
  have h1 := Nyθ n hN
  have h2 := Ny2 n hN
  unfold uv γ0
  rw [mul_pow, mul_pow, exp_neg_pow3]
  have h3 : √2 ^ 3 = 2 * √2 := by
    rw [pow_succ, Real.sq_sqrt (by norm_num)]
  rw [h3, ← h1]
  have : (n : ℝ) * yv n ^ 3 = π / 2 * yv n := by rw [pow_succ, ← mul_assoc, h2]
  linear_combination (-(√2 * Real.exp (-(3 * σ)) / 6)) * this

lemma norm_f_le (n : ℕ) (hn : 0 < n) (σ : ℝ) :
    ‖fσ n σ‖ ≤ √2 * yv n * Real.exp (-σ) *
      Real.exp (-(θv n * (Real.exp σ + Real.exp (-σ)) / 2) + γ0 * yv n * Real.exp (-(3 * σ))) := by
  have hN : (0 : ℝ) < n := by exact_mod_cast hn
  rw [norm_fσ n hn σ]
  have hlog := log_cubic (uv n σ) (uv_nonneg n σ hN)
  have hh := half_n_u n hn σ
  have hexp : Real.exp (-(θv n * Real.exp σ / 2 + n / 2 * Real.log (pv n σ))) ≤
      Real.exp (-(θv n * (Real.exp σ + Real.exp (-σ)) / 2) + γ0 * yv n * Real.exp (-(3 * σ))) := by
    apply Real.exp_le_exp.2
    have : (n : ℝ) / 2 * (uv n σ - uv n σ ^ 3 / 6) ≤ n / 2 * Real.log (pv n σ) :=
      mul_le_mul_of_nonneg_left hlog (by positivity)
    unfold pv at *
    linarith
  exact mul_le_mul (H_sv n σ hN) hexp (Real.exp_pos _).le (by have := yv_pos n hN; positivity)

/-- The slope `S` of the left tail at `-δ`. -/
noncomputable def Sl (n : ℕ) (δ : ℝ) : ℝ :=
  n / 2 * ((uv n (-δ) + uv n (-δ) ^ 2) / pv n (-δ)) - θv n * Real.exp (-δ) / 2

lemma norm_f_left (n : ℕ) (hn : 0 < n) (δ : ℝ) (hu : uv n (-δ) ^ 2 ≤ 2) (σ : ℝ)
    (hσ : σ ≤ -δ) :
    ‖fσ n σ‖ ≤ √2 * yv n * Real.exp (-σ) *
      Real.exp (-(θv n * (Real.exp δ + Real.exp (-δ)) / 2) + γ0 * yv n * Real.exp (3 * δ) -
        Sl n δ * (-δ - σ)) := by
  have hN : (0 : ℝ) < n := by exact_mod_cast hn
  rw [norm_fσ n hn σ]
  set u1 := uv n (-δ)
  set P1 := pv n (-δ)
  set τ := -δ - σ
  have hτ : 0 ≤ τ := by simp only [τ]; linarith
  set t := Real.exp τ
  have ht : 1 ≤ t := Real.one_le_exp hτ
  have hu1 : 0 ≤ u1 := uv_nonneg n (-δ) hN
  have hP1 : 0 < P1 := pv_pos n (-δ) hN
  set β := (u1 + u1 ^ 2) / P1
  have hβ0 : 0 ≤ β := by positivity
  have hβ1 : β ≤ 1 := by
    rw [div_le_one hP1]; simp only [P1, pv]; nlinarith
  have huσ : uv n σ = u1 * t := by
    simp only [u1, t, τ, uv]
    have : Real.exp (- -δ) * Real.exp (-δ - σ) = Real.exp (-σ) := by rw [← Real.exp_add]; ring_nf
    rw [mul_assoc (√2 * yv n), this]
  have hp : P1 * (1 + β * (t - 1)) ≤ pv n σ := by
    have : P1 * (1 + β * (t - 1)) = P1 + (u1 + u1 ^ 2) * (t - 1) := by
      simp only [β]; field_simp
    rw [this]; unfold pv; rw [huσ]; simp only [P1, pv]
    nlinarith [sq_nonneg (u1 * (t - 1))]
  have hrp : t ^ β ≤ 1 + β * (t - 1) := by
    have := rpow_one_add_le_one_add_mul_self (by linarith : (-1 : ℝ) ≤ t - 1) hβ0 hβ1
    simpa using this
  have ht0 : 0 < t := Real.exp_pos _
  have hlog : Real.log P1 + β * τ ≤ Real.log (pv n σ) := by
    have h1 : P1 * t ^ β ≤ pv n σ := (mul_le_mul_of_nonneg_left hrp hP1.le).trans hp
    have h2 := Real.log_le_log (by positivity) h1
    rwa [Real.log_mul hP1.ne' (Real.rpow_pos_of_pos ht0 β).ne', Real.log_rpow ht0, Real.log_exp] at h2
  have hlog1 := log_cubic u1 hu1
  have hh := half_n_u n hn (-δ)
  rw [neg_neg, show -(3 * -δ) = 3 * δ by ring] at hh
  have hE : Real.exp (-δ) * (1 - τ) ≤ Real.exp σ := by
    have h1 : 1 - τ ≤ Real.exp (-τ) := by linarith [Real.add_one_le_exp (-τ)]
    have h2 : Real.exp σ = Real.exp (-δ) * Real.exp (-τ) := by
      rw [← Real.exp_add]; congr 1; simp only [τ]; ring
    rw [h2]; exact mul_le_mul_of_nonneg_left h1 (Real.exp_pos _).le
  have hθ := θv_pos n hN
  have hexp : Real.exp (-(θv n * Real.exp σ / 2 + n / 2 * Real.log (pv n σ))) ≤
      Real.exp (-(θv n * (Real.exp δ + Real.exp (-δ)) / 2) + γ0 * yv n * Real.exp (3 * δ) -
        Sl n δ * τ) := by
    apply Real.exp_le_exp.2
    have a1 : (n : ℝ) / 2 * (Real.log P1 + β * τ) ≤ n / 2 * Real.log (pv n σ) :=
      mul_le_mul_of_nonneg_left hlog (by positivity)
    have a2 : (n : ℝ) / 2 * (u1 - u1 ^ 3 / 6) ≤ n / 2 * Real.log P1 :=
      mul_le_mul_of_nonneg_left hlog1 (by positivity)
    have a3 : θv n * (Real.exp (-δ) * (1 - τ)) ≤ θv n * Real.exp σ :=
      mul_le_mul_of_nonneg_left hE hθ.le
    have hS : Sl n δ = n / 2 * β - θv n * Real.exp (-δ) / 2 := rfl
    rw [hS]
    simp only [u1] at a2
    nlinarith
  exact mul_le_mul (H_sv n σ hN) hexp (Real.exp_pos _).le (by have := yv_pos n hN; positivity)

/-! ## 4. Tail integrals -/

lemma left_tail (n : ℕ) (hn : 0 < n) (δ : ℝ) (hu : uv n (-δ) ^ 2 ≤ 2) (hS : 1 < Sl n δ) :
    ‖∫ σ in Iic (-δ), fσ n σ‖ ≤ √2 * yv n * Real.exp δ *
      Real.exp (-(θv n * (Real.exp δ + Real.exp (-δ)) / 2) + γ0 * yv n * Real.exp (3 * δ)) /
        (Sl n δ - 1) := by
  have hN : (0 : ℝ) < n := by exact_mod_cast hn
  set A := -(θv n * (Real.exp δ + Real.exp (-δ)) / 2) + γ0 * yv n * Real.exp (3 * δ)
  set S := Sl n δ
  set C := √2 * yv n * Real.exp (A + S * δ)
  have hg : IntegrableOn (fun σ : ℝ => C * Real.exp ((S - 1) * σ)) (Iic (-δ)) :=
    (integrableOn_exp_mul_Iic (by linarith) (-δ)).const_mul C
  have hb : ∀ σ ∈ Iic (-δ), ‖fσ n σ‖ ≤ C * Real.exp ((S - 1) * σ) := by
    intro σ hσ
    refine (norm_f_left n hn δ hu σ hσ).trans (le_of_eq ?_)
    simp only [C]
    rw [mul_assoc (√2 * yv n), mul_assoc (√2 * yv n), ← Real.exp_add, ← Real.exp_add]
    congr 2; ring
  have h1 := norm_integral_le_of_norm_le hg
    ((ae_restrict_iff' measurableSet_Iic).2 (Filter.Eventually.of_forall hb))
  refine h1.trans (le_of_eq ?_)
  rw [integral_const_mul, integral_exp_mul_Iic (by linarith) (-δ)]
  simp only [C]
  rw [mul_div_assoc', mul_assoc (√2 * yv n), ← Real.exp_add]
  have : A + S * δ + (S - 1) * -δ = δ + A := by ring
  rw [this, Real.exp_add]; ring

lemma cosh_lin (δ σ : ℝ) :
    (Real.exp δ + Real.exp (-δ)) / 2 + (Real.exp δ - Real.exp (-δ)) / 2 * (σ - δ) ≤
      (Real.exp σ + Real.exp (-σ)) / 2 := by
  have h1 : Real.exp σ = Real.exp δ * Real.exp (σ - δ) := by rw [← Real.exp_add]; ring_nf
  have h2 : Real.exp (-σ) = Real.exp (-δ) * Real.exp (-(σ - δ)) := by rw [← Real.exp_add]; ring_nf
  have h3 := Real.add_one_le_exp (σ - δ)
  have h4 := Real.add_one_le_exp (-(σ - δ))
  rw [h1, h2]
  have := Real.exp_pos δ; have := Real.exp_pos (-δ)
  nlinarith

lemma right_tail (n : ℕ) (hn : 0 < n) (δ : ℝ) (hδ : 0 ≤ δ) :
    ‖∫ σ in Ioi δ, fσ n σ‖ ≤ √2 * yv n * Real.exp (-δ) *
      Real.exp (-(θv n * (Real.exp δ + Real.exp (-δ)) / 2) + γ0 * yv n * Real.exp (-(3 * δ))) /
        (1 + θv n * ((Real.exp δ - Real.exp (-δ)) / 2)) := by
  have hN : (0 : ℝ) < n := by exact_mod_cast hn
  have hθ := θv_pos n hN
  have hy := yv_pos n hN
  set c := (Real.exp δ + Real.exp (-δ)) / 2
  set sh := (Real.exp δ - Real.exp (-δ)) / 2
  have hsh : 0 ≤ sh := by
    simp only [sh]; have := Real.exp_le_exp.2 (by linarith : -δ ≤ δ); linarith
  set A := -(θv n * c) + γ0 * yv n * Real.exp (-(3 * δ))
  set L := 1 + θv n * sh
  have hL : 0 < L := by simp only [L]; positivity
  set C := √2 * yv n * Real.exp (A + θv n * sh * δ)
  have hg : IntegrableOn (fun σ : ℝ => C * Real.exp (-L * σ)) (Ioi δ) :=
    (integrableOn_exp_mul_Ioi (by linarith) δ).const_mul C
  have hb : ∀ σ ∈ Ioi δ, ‖fσ n σ‖ ≤ C * Real.exp (-L * σ) := by
    intro σ hσ
    have hσ' : δ < σ := hσ
    refine (norm_f_le n hn σ).trans ?_
    simp only [C]
    rw [mul_assoc (√2 * yv n), mul_assoc (√2 * yv n), ← Real.exp_add, ← Real.exp_add]
    apply mul_le_mul_of_nonneg_left _ (by positivity)
    apply Real.exp_le_exp.2
    have h1 := cosh_lin δ σ
    have h2 : Real.exp (-(3 * σ)) ≤ Real.exp (-(3 * δ)) := Real.exp_le_exp.2 (by linarith)
    have h3 : 0 ≤ γ0 := by unfold γ0; positivity
    simp only [A, L, c, sh] at *
    have h4 : γ0 * yv n * Real.exp (-(3 * σ)) ≤ γ0 * yv n * Real.exp (-(3 * δ)) :=
      mul_le_mul_of_nonneg_left h2 (by positivity)
    have h5 := mul_le_mul_of_nonneg_left h1 hθ.le
    nlinarith
  have h1 := norm_integral_le_of_norm_le hg
    ((ae_restrict_iff' measurableSet_Ioi).2 (Filter.Eventually.of_forall hb))
  refine h1.trans (le_of_eq ?_)
  rw [integral_const_mul, integral_exp_mul_Ioi (by linarith) δ]
  simp only [C]
  have e : Real.exp (A + θv n * sh * δ) * Real.exp (-L * δ) = Real.exp (-δ) * Real.exp A := by
    rw [← Real.exp_add, ← Real.exp_add]; congr 1; simp only [L]; ring
  rw [mul_assoc (√2 * yv n), neg_div_neg_eq, mul_div_assoc', e]; simp only [A, c]; ring_nf

lemma sinh_sq_ge (t : ℝ) : t ^ 2 ≤ Real.sinh t ^ 2 := by
  rcases le_total 0 t with h | h
  · have := Real.self_le_sinh_iff.2 h
    nlinarith
  · have h' : 0 ≤ -t := by linarith
    have := Real.self_le_sinh_iff.2 h'
    rw [Real.sinh_neg] at this
    nlinarith

lemma cosh_quad (σ : ℝ) : 1 + σ ^ 2 / 2 ≤ (Real.exp σ + Real.exp (-σ)) / 2 := by
  rw [← Real.cosh_eq]
  have h1 : Real.cosh σ = 1 + 2 * Real.sinh (σ / 2) ^ 2 := by
    have := Real.cosh_two_mul (σ / 2)
    rw [show 2 * (σ / 2) = σ by ring] at this
    rw [this, Real.cosh_sq]; ring
  have h2 := sinh_sq_ge (σ / 2)
  rw [h1]; nlinarith

lemma gauss_lin (b : ℝ) (hb : 0 < b) :
    Integrable (fun σ : ℝ => Real.exp (-σ - b * σ ^ 2)) ∧
      ∫ σ : ℝ, Real.exp (-σ - b * σ ^ 2) = √(π / b) * Real.exp (1 / (4 * b)) := by
  have heq : (fun σ : ℝ => Real.exp (-σ - b * σ ^ 2)) =
      fun σ : ℝ => Real.exp (1 / (4 * b)) * Real.exp (-b * (σ + 1 / (2 * b)) ^ 2) := by
    funext σ
    rw [← Real.exp_add]; congr 1; field_simp; ring
  rw [heq]
  have hi := integrable_exp_neg_mul_sq hb
  refine ⟨(hi.comp_add_right (1 / (2 * b))).const_mul _, ?_⟩
  rw [integral_const_mul, integral_add_right_eq_self (fun σ : ℝ => Real.exp (-b * σ ^ 2)),
    integral_gaussian]
  ring

lemma split_Iic_Ioi (g : ℝ → ℂ) (hg : Integrable g) (c : ℝ) :
    ∫ σ, g σ = (∫ σ in Iic c, g σ) + ∫ σ in Ioi c, g σ := by
  rw [← setIntegral_univ, ← Iic_union_Ioi (a := c)]
  exact setIntegral_union (Iic_disjoint_Ioi le_rfl) measurableSet_Ioi hg.integrableOn hg.integrableOn

lemma mid_bound (n : ℕ) (hn : 0 < n) (δ : ℝ) :
    ‖∫ σ in Ioi (-δ), fσ n σ‖ ≤ √2 * yv n * Real.exp (γ0 * yv n * Real.exp (3 * δ)) *
      Real.exp (-θv n) * (√(π / (θv n / 2)) * Real.exp (1 / (4 * (θv n / 2)))) := by
  have hN : (0 : ℝ) < n := by exact_mod_cast hn
  have hθ := θv_pos n hN
  have hy := yv_pos n hN
  set C := √2 * yv n * Real.exp (γ0 * yv n * Real.exp (3 * δ)) * Real.exp (-θv n)
  obtain ⟨hgi, hgv⟩ := gauss_lin (θv n / 2) (by positivity)
  have hg : Integrable (fun σ : ℝ => C * Real.exp (-σ - θv n / 2 * σ ^ 2)) := hgi.const_mul C
  have hb : ∀ σ ∈ Ioi (-δ), ‖fσ n σ‖ ≤ C * Real.exp (-σ - θv n / 2 * σ ^ 2) := by
    intro σ hσ
    have hσ' : -δ < σ := hσ
    refine (norm_f_le n hn σ).trans (le_of_eq_of_le rfl ?_)
    simp only [C]
    have e1 : √2 * yv n * Real.exp (γ0 * yv n * Real.exp (3 * δ)) * Real.exp (-θv n) *
        Real.exp (-σ - θv n / 2 * σ ^ 2) = √2 * yv n * Real.exp (-σ) *
        Real.exp (γ0 * yv n * Real.exp (3 * δ) - θv n - θv n / 2 * σ ^ 2) := by
      rw [mul_assoc (√2 * yv n), mul_assoc (√2 * yv n), mul_assoc (√2 * yv n), ← Real.exp_add,
        ← Real.exp_add, ← Real.exp_add]
      congr 2; ring
    rw [e1]
    apply mul_le_mul_of_nonneg_left _ (by positivity)
    apply Real.exp_le_exp.2
    have h1 := cosh_quad σ
    have h2 : Real.exp (-(3 * σ)) ≤ Real.exp (3 * δ) := Real.exp_le_exp.2 (by linarith)
    have h3 : 0 ≤ γ0 := by unfold γ0; positivity
    have h4 : γ0 * yv n * Real.exp (-(3 * σ)) ≤ γ0 * yv n * Real.exp (3 * δ) :=
      mul_le_mul_of_nonneg_left h2 (by positivity)
    have h5 := mul_le_mul_of_nonneg_left h1 hθ.le
    nlinarith
  have h1 := norm_integral_le_of_norm_le hg.integrableOn
    ((ae_restrict_iff' measurableSet_Ioi).2 (Filter.Eventually.of_forall hb))
  refine h1.trans ?_
  have hC : 0 ≤ C := by simp only [C]; positivity
  refine (setIntegral_le_integral hg (Filter.Eventually.of_forall fun σ => by
    simp only [Pi.zero_apply]; positivity)).trans (le_of_eq ?_)
  rw [integral_const_mul, hgv]

/-! ## 5. Numerical helpers -/

lemma exp_le_quad (x : ℝ) (hx : |x| ≤ 1) : Real.exp x ≤ 1 + x + x ^ 2 * (3 / 4) := by
  have := Real.exp_bound hx (n := 2) (by norm_num)
  simp [Finset.sum_range_succ] at this
  have := (abs_le.1 this).2
  norm_num at this
  linarith

lemma exp_quarter_lt : Real.exp (1 / 4) < 1.2841 := by
  have := Real.exp_bound (x := 1 / 4) (by norm_num [abs_of_pos]) (n := 5) (by norm_num)
  simp [Finset.sum_range_succ, Nat.factorial] at this
  rw [abs_le] at this
  norm_num at this
  linarith [this.2]

lemma exp_quarter_gt : 1.284 < Real.exp (1 / 4) := by
  have := Real.exp_bound (x := 1 / 4) (by norm_num [abs_of_pos]) (n := 5) (by norm_num)
  simp [Finset.sum_range_succ, Nat.factorial] at this
  rw [abs_le] at this
  norm_num at this
  linarith [this.1]

lemma exp_neg_quarter_lt : Real.exp (-(1 / 4)) < 0.77883 := by
  rw [Real.exp_neg, inv_lt_comm₀ (Real.exp_pos _) (by norm_num)]
  have := exp_quarter_gt; norm_num at this ⊢; linarith

lemma exp_neg_quarter_gt : 0.7787 < Real.exp (-(1 / 4)) := by
  rw [Real.exp_neg, lt_inv_comm₀ (by norm_num) (Real.exp_pos _)]
  have := exp_quarter_lt; norm_num at this ⊢; linarith

lemma exp_three_quarter_lt : Real.exp (3 * (1 / 4)) < 2.1174 := by
  have h : Real.exp (3 * (1 / 4)) = Real.exp (1 / 4) ^ 3 := by
    rw [← Real.exp_nat_mul]; norm_num
  rw [h]
  have := exp_quarter_lt
  have h0 := Real.exp_pos (1 / 4)
  calc Real.exp (1 / 4) ^ 3 < 1.2841 ^ 3 := by gcongr
    _ < 2.1174 := by norm_num

lemma γ0_lt : γ0 < 0.37025 := by
  unfold γ0
  have h2 : √2 < 1.41422 := by
    rw [Real.sqrt_lt' (by norm_num)]; norm_num
  have := Real.pi_lt_d6
  have h0 : 0 < √2 := sqrt2_pos
  nlinarith

lemma γ0_pos : 0 < γ0 := by unfold γ0; positivity

lemma sqrt2_lt : √2 < 1.41422 := by rw [Real.sqrt_lt' (by norm_num)]; norm_num

lemma sqrt2_gt : 1.41421 < √2 := by rw [Real.lt_sqrt (by norm_num)]; norm_num

lemma uv_neg (n : ℕ) (δ : ℝ) : uv n (-δ) = √2 * yv n * Real.exp δ := by
  unfold uv; rw [neg_neg]

lemma Sl_ge (n : ℕ) (hn : 0 < n) (δ : ℝ) :
    θv n * (Real.exp δ / 2 * (1 - uv n (-δ) ^ 2 / 2) - Real.exp (-δ) / 2) ≤ Sl n δ := by
  have hN : (0 : ℝ) < n := by exact_mod_cast hn
  have hu := uv_nonneg n (-δ) hN
  set u := uv n (-δ)
  have hP := pv_pos n (-δ) hN
  have hnu : (n : ℝ) / 2 * u = θv n * Real.exp δ / 2 := by
    simp only [u]; rw [uv_neg, ← Nyθ n hN]; ring
  have key : u * (1 - u ^ 2 / 2) ≤ (u + u ^ 2) / pv n (-δ) := by
    rw [le_div_iff₀ hP]
    simp only [pv]
    nlinarith [sq_nonneg u, mul_nonneg hu (sq_nonneg u), pow_nonneg hu 3, pow_nonneg hu 4,
      pow_nonneg hu 5]
  unfold Sl
  have := mul_le_mul_of_nonneg_left key (by positivity : (0 : ℝ) ≤ n / 2)
  have e : (n : ℝ) / 2 * (u * (1 - u ^ 2 / 2)) = θv n * (Real.exp δ / 2 * (1 - u ^ 2 / 2)) := by
    rw [← mul_assoc, hnu]; ring
  rw [e] at this
  linarith

/-! ## 6. Region I: `100 ≤ N ≤ 2000` -/

lemma Mmain_eq (k : ℕ) (hk : 1 ≤ k) : Mmain k = 2 / π * (∫ σ : ℝ, fσ (k + 1) σ).im := by
  unfold Mmain
  congr 2
  rw [← (J_eq k hk).2]
  refine setIntegral_congr_fun measurableSet_Ioi (fun s hs => ?_)
  have := (Q05.proof_11 k hk).2.2.1 s hs
  simp only [fM]
  rw [mul_assoc (H s * E s), this]

lemma varpi1_le : ‖varpi1‖ ≤ 0.843 := by
  have h := Q05.proof_3
  unfold Q05.stmt_3 at h
  have := Q05.sqrt_pi_gt
  have h0 := norm_nonneg varpi1
  nlinarith

lemma X_le (N : ℝ) (hN : 0 < N) :
    |X N| ≤ √π / N * Real.exp (-√(π * N)) * (1 + ‖varpi1‖ / √N) := by
  unfold X
  rw [abs_mul]
  rw [abs_of_nonneg (by positivity)]
  apply mul_le_mul_of_nonneg_left _ (by positivity)
  refine (Complex.abs_re_le_norm _).trans ?_
  rw [norm_mul, mul_comm Complex.I, Complex.norm_exp_ofReal_mul_I, one_mul]
  refine (norm_add_le _ _).trans ?_
  rw [norm_one, norm_div, Complex.norm_real, Real.norm_eq_abs, abs_of_nonneg (Real.sqrt_nonneg _)]

lemma regI_left_num (θ y n S g c : ℝ) (hθlo : 17.72 ≤ θ) (hn : 0 < n)
    (hyθ : √2 * y = θ / n) (hS : 0.2358 * θ ≤ S) (hg : Real.exp g ≤ 1.1057) (hc : 1 ≤ c) :
    2 / π * (√2 * y * Real.exp (1 / 4) * Real.exp (-(θ * c) + g) / (S - 1)) ≤
      5.1 / n * Real.exp (-θ) := by
  have hpi := Real.pi_gt_d6
  have e14 := exp_quarter_lt
  have hθ : 0 < θ := by linarith
  have hexp : Real.exp (-(θ * c) + g) ≤ Real.exp (-θ) * Real.exp g := by
    rw [← Real.exp_add]; apply Real.exp_le_exp.2; nlinarith
  have hS1' : 0 < S - 1 := by nlinarith
  have hθS : θ / (S - 1) ≤ 5.58 := by
    rw [div_le_iff₀ hS1']; nlinarith
  rw [hyθ]
  have hE := Real.exp_pos (-θ)
  calc 2 / π * (θ / n * Real.exp (1 / 4) * Real.exp (-(θ * c) + g) / (S - 1))
      ≤ 2 / π * (θ / n * Real.exp (1 / 4) * (Real.exp (-θ) * Real.exp g) / (S - 1)) := by
        gcongr
    _ = (2 / π * Real.exp (1 / 4) * Real.exp g * (θ / (S - 1))) / n * Real.exp (-θ) := by
        field_simp
    _ ≤ (2 / 3.14159 * 1.2841 * 1.1057 * 5.58) / n * Real.exp (-θ) := by
        gcongr; linarith
    _ ≤ 5.1 / n * Real.exp (-θ) := by
        gcongr; norm_num

lemma regI_mid_num (θ y n g : ℝ) (hθlo : 17.72 ≤ θ) (hθhi : θ ≤ 80) (hn : 0 < n)
    (hyθ : √2 * y = θ / n) (hg0 : 0 ≤ g) (hg : g ≤ 0.0984) :
    2 / π * (√2 * y * Real.exp g * Real.exp (-θ) *
      (√(π / (θ / 2)) * Real.exp (1 / (4 * (θ / 2))))) ≤ 16.3 / n * Real.exp (-θ) := by
  have hpi := Real.pi_gt_d6
  have hpi' := Real.pi_lt_d6
  have hθ : 0 < θ := by linarith
  rw [hyθ]
  have hsq : θ * √(π / (θ / 2)) ≤ 22.43 := by
    have h0 : 0 ≤ θ * √(π / (θ / 2)) := by positivity
    have : (θ * √(π / (θ / 2))) ^ 2 = 2 * π * θ := by
      rw [mul_pow, Real.sq_sqrt (by positivity)]; field_simp
    nlinarith
  have h4 : 1 / (4 * (θ / 2)) ≤ 0.0283 := by
    rw [div_le_iff₀ (by positivity)]; nlinarith
  have he2 : Real.exp g * Real.exp (1 / (4 * (θ / 2))) ≤ 1.139 := by
    rw [← Real.exp_add]
    have h5 : 0 ≤ 1 / (4 * (θ / 2)) := by positivity
    have := exp_le_quad (g + 1 / (4 * (θ / 2))) (by rw [abs_le]; constructor <;> linarith)
    have h6 : (g + 1 / (4 * (θ / 2))) ^ 2 ≤ 0.1267 ^ 2 := by
      apply pow_le_pow_left₀ (by positivity); linarith
    nlinarith
  have hE := Real.exp_pos (-θ)
  calc 2 / π * (θ / n * Real.exp g * Real.exp (-θ) *
        (√(π / (θ / 2)) * Real.exp (1 / (4 * (θ / 2)))))
      = (2 / π * (θ * √(π / (θ / 2))) * (Real.exp g * Real.exp (1 / (4 * (θ / 2))))) / n *
        Real.exp (-θ) := by field_simp
    _ ≤ (2 / 3.14159 * 22.43 * 1.139) / n * Real.exp (-θ) := by gcongr; linarith
    _ ≤ 16.3 / n * Real.exp (-θ) := by gcongr; norm_num

lemma regI_num (n : ℕ) (h1 : 100 ≤ n) (h2 : n ≤ 2000) :
    2 / π * ‖∫ σ : ℝ, fσ n σ‖ ≤ 21.4 / n * Real.exp (-θv n) := by
  have hn : 0 < n := by omega
  have hN : (0 : ℝ) < n := by exact_mod_cast hn
  have hN1 : (100 : ℝ) ≤ n := by exact_mod_cast h1
  have hN2 : (n : ℝ) ≤ 2000 := by exact_mod_cast h2
  have hθ := θv_pos n hN
  have hy := yv_pos n hN
  have hpi := Real.pi_gt_d6
  have hpi' := Real.pi_lt_d6
  have hθlo : 17.72 ≤ θv n := by
    have := θ_sq n hN; nlinarith
  have hθhi : θv n ≤ 80 := by
    have := θ_sq n hN; nlinarith
  have hylo : yv n ≤ 0.1254 := by
    have h := y_sq n hN
    have : π / (2 * n) ≤ 0.015709 := by
      rw [div_le_iff₀ (by positivity)]; nlinarith
    nlinarith
  have hu1 : uv n (-(1 / 4)) ≤ 0.2278 := by
    rw [uv_neg]
    have : √2 * yv n ≤ 1.41422 * 0.1254 := by
      have := sqrt2_lt; gcongr
    have := Real.exp_pos (1 / 4)
    have := exp_quarter_lt
    nlinarith
  have hu0 := uv_nonneg n (-(1 / 4)) hN
  have hu2 : uv n (-(1 / 4)) ^ 2 ≤ 2 := by nlinarith
  have hS := Sl_ge n hn (1 / 4)
  have hSθ : 0.2358 * θv n ≤ Sl n (1 / 4) := by
    refine le_trans ?_ hS
    rw [mul_comm (0.2358 : ℝ)]
    apply mul_le_mul_of_nonneg_left _ hθ.le
    have := exp_quarter_gt
    have := exp_neg_quarter_lt
    nlinarith
  have hS1 : 1 < Sl n (1 / 4) := by nlinarith
  have hL := left_tail n hn (1 / 4) hu2 hS1
  have hM := mid_bound n hn (1 / 4)
  have hint : Integrable (fσ n) := by
    obtain ⟨k, rfl⟩ : ∃ k, n = k + 1 := ⟨n - 1, by omega⟩
    exact (J_eq k (by omega)).1
  rw [split_Iic_Ioi (fσ n) hint (-(1 / 4))]
  have htri := norm_add_le (∫ σ in Iic (-(1 / 4)), fσ n σ) (∫ σ in Ioi (-(1 / 4)), fσ n σ)
  have hnyθ := Nyθ n hN
  have hyθ : √2 * yv n = θv n / n := by field_simp; linarith
  have hg : γ0 * yv n * Real.exp (3 * (1 / 4)) ≤ 0.0984 := by
    have : γ0 * yv n ≤ 0.37025 * 0.1254 := by
      have := γ0_lt; have := γ0_pos; gcongr
    have := Real.exp_pos (3 * (1 / 4))
    have := exp_three_quarter_lt
    nlinarith
  have hg0 : 0 ≤ γ0 * yv n * Real.exp (3 * (1 / 4)) := by have := γ0_pos; positivity
  have heg : Real.exp (γ0 * yv n * Real.exp (3 * (1 / 4))) ≤ 1.1057 := by
    have := exp_le_quad (γ0 * yv n * Real.exp (3 * (1 / 4)))
      (by rw [abs_le]; constructor <;> linarith)
    nlinarith
  have hcosh : 1 ≤ (Real.exp (1 / 4) + Real.exp (-(1 / 4))) / 2 := by
    have := cosh_quad (1 / 4); nlinarith
  have hLt := regI_left_num (θv n) (yv n) n (Sl n (1 / 4)) (γ0 * yv n * Real.exp (3 * (1 / 4)))
    ((Real.exp (1 / 4) + Real.exp (-(1 / 4))) / 2) hθlo hN hyθ hSθ heg hcosh
  have hMt := regI_mid_num (θv n) (yv n) n (γ0 * yv n * Real.exp (3 * (1 / 4))) hθlo hθhi hN
    hyθ hg0 hg
  have e1 : -(θv n * (Real.exp (1 / 4) + Real.exp (-(1 / 4))) / 2) =
      -(θv n * ((Real.exp (1 / 4) + Real.exp (-(1 / 4))) / 2)) := by ring
  rw [e1] at hL
  have hpos : (0 : ℝ) ≤ 2 / π := by positivity
  calc 2 / π * ‖(∫ σ in Iic (-(1 / 4)), fσ n σ) + ∫ σ in Ioi (-(1 / 4)), fσ n σ‖
      ≤ 2 / π * (√2 * yv n * Real.exp (1 / 4) *
          Real.exp (-(θv n * ((Real.exp (1 / 4) + Real.exp (-(1 / 4))) / 2)) +
            γ0 * yv n * Real.exp (3 * (1 / 4))) / (Sl n (1 / 4) - 1)) +
        2 / π * (√2 * yv n * Real.exp (γ0 * yv n * Real.exp (3 * (1 / 4))) * Real.exp (-θv n) *
          (√(π / (θv n / 2)) * Real.exp (1 / (4 * (θv n / 2))))) := by
        rw [← mul_add]
        exact mul_le_mul_of_nonneg_left (htri.trans (add_le_add hL hM)) hpos
    _ ≤ 5.1 / n * Real.exp (-θv n) + 16.3 / n * Real.exp (-θv n) := add_le_add hLt hMt
    _ ≤ 21.4 / n * Real.exp (-θv n) := by
        rw [← add_mul, ← add_div]; gcongr; norm_num

lemma X_le_num (n : ℕ) (h1 : 100 ≤ n) : |X n| ≤ 1.923 / n * Real.exp (-θv n) := by
  have hN : (0 : ℝ) < n := by have : (100 : ℝ) ≤ n := by exact_mod_cast h1
                              linarith
  have hN1 : (100 : ℝ) ≤ n := by exact_mod_cast h1
  refine (X_le n hN).trans ?_
  have hs : 10 ≤ √(n : ℝ) := by rw [Real.le_sqrt (by norm_num) hN.le]; linarith
  have hv := varpi1_le
  have hv0 := norm_nonneg varpi1
  have h3 : ‖varpi1‖ / √(n : ℝ) ≤ 0.0843 := by
    rw [div_le_iff₀ (by linarith)]; nlinarith
  have hsp := Q05.sqrt_pi_lt
  have hE := Real.exp_pos (-θv n)
  unfold θv
  calc √π / n * Real.exp (-√(π * n)) * (1 + ‖varpi1‖ / √(n : ℝ))
      ≤ 1.7724539 / n * Real.exp (-√(π * n)) * (1 + 0.0843) := by gcongr
    _ ≤ 1.923 / n * Real.exp (-√(π * n)) := by
        rw [mul_assoc, mul_comm (Real.exp _), ← mul_assoc]
        have : 1.7724539 / (n : ℝ) * (1 + 0.0843) ≤ 1.923 / n := by
          rw [div_mul_eq_mul_div]; gcongr; norm_num
        exact mul_le_mul_of_nonneg_right this (Real.exp_pos _).le

lemma regI (k : ℕ) (h1 : 100 ≤ k + 1) (h2 : k + 1 ≤ 2000) :
    |Mmain k + X ((k + 1 : ℕ) : ℝ)| ≤
      71000 / (((k + 1 : ℕ) : ℝ)) ^ 2 * Real.exp (-√(π * ((k + 1 : ℕ) : ℝ))) := by
  have hk : 1 ≤ k := by omega
  set n := k + 1 with hn
  have hN1 : (100 : ℝ) ≤ n := by exact_mod_cast h1
  have hN2 : (n : ℝ) ≤ 2000 := by exact_mod_cast h2
  have hN : (0 : ℝ) < n := by linarith
  have hJ := regI_num n h1 h2
  have hXb := X_le_num n h1
  have hM : |Mmain k| ≤ 21.4 / n * Real.exp (-θv n) := by
    rw [Mmain_eq k hk, abs_mul, abs_of_pos (by positivity : (0 : ℝ) < 2 / π)]
    exact (mul_le_mul_of_nonneg_left (Complex.abs_im_le_norm _) (by positivity)).trans hJ
  have hE := Real.exp_pos (-θv n)
  calc |Mmain k + X n| ≤ |Mmain k| + |X n| := abs_add_le _ _
    _ ≤ 21.4 / n * Real.exp (-θv n) + 1.923 / n * Real.exp (-θv n) := add_le_add hM hXb
    _ = (23.323 * n) / n ^ 2 * Real.exp (-θv n) := by field_simp; ring
    _ ≤ 71000 / n ^ 2 * Real.exp (-θv n) := by gcongr; linarith

/-! ## 7. Gaussian moments -/

lemma integrable_pow_gauss (m : ℕ) (c : ℝ) (hc : 0 < c) :
    Integrable (fun x : ℝ => x ^ m * Real.exp (-c * x ^ 2)) := by
  have := integrable_rpow_mul_exp_neg_mul_sq hc (s := (m : ℝ)) (by linarith [m.cast_nonneg (α := ℝ)])
  simpa [Real.rpow_natCast] using this

lemma integrable_cpow_gauss (m : ℕ) (b : ℂ) (hb : 0 < b.re) :
    Integrable (fun x : ℝ => (x : ℂ) ^ m * Complex.exp (-b * (x : ℂ) ^ 2)) := by
  refine (integrable_pow_gauss m b.re hb).norm.mono' (by fun_prop) ?_
  refine Filter.Eventually.of_forall fun x => le_of_eq ?_
  simp only [norm_mul, norm_pow, Complex.norm_real, Complex.norm_exp, Real.norm_eq_abs,
    Real.abs_exp]
  congr 2
  simp [Complex.mul_re, sq]

lemma gauss_rec (m : ℕ) (b : ℂ) (hb : 0 < b.re) :
    ∫ x : ℝ, (x : ℂ) ^ (m + 2) * Complex.exp (-b * (x : ℂ) ^ 2) =
      (m + 1) / (2 * b) * ∫ x : ℝ, (x : ℂ) ^ m * Complex.exp (-b * (x : ℂ) ^ 2) := by
  have hb0 : b ≠ 0 := by intro h; rw [h] at hb; simp at hb
  set g : ℂ → ℂ := fun w => w ^ (m + 1) * Complex.exp (-b * w ^ 2)
  set g' : ℂ → ℂ := fun w => (m + 1) * (w ^ m * Complex.exp (-b * w ^ 2)) -
    2 * b * (w ^ (m + 2) * Complex.exp (-b * w ^ 2))
  have hd : ∀ w : ℂ, HasDerivAt g (g' w) w := by
    intro w
    have h1 := hasDerivAt_pow (m + 1) w
    have h2 : HasDerivAt (fun w : ℂ => -b * w ^ 2) (-b * (2 * w)) w := by
      have := (hasDerivAt_pow 2 w).const_mul (-b)
      simpa using this
    have h3 := (h1.mul h2.cexp)
    convert h3 using 1
    simp only [g']
    push_cast
    ring
  have hdr : ∀ x : ℝ, HasDerivAt (fun x : ℝ => g x) (g' x) x := fun x => (hd x).comp_ofReal
  have hI1 := integrable_cpow_gauss m b hb
  have hI2 := integrable_cpow_gauss (m + 2) b hb
  have hI3 := integrable_cpow_gauss (m + 1) b hb
  have hg'i : Integrable (fun x : ℝ => g' x) := by
    simp only [g']
    exact (hI1.const_mul _).sub (hI2.const_mul _)
  have h0 := integral_eq_zero_of_hasDerivAt_of_integrable hdr hg'i hI3
  simp only [g'] at h0
  rw [integral_sub (hI1.const_mul _) (hI2.const_mul _), integral_const_mul, integral_const_mul,
    sub_eq_zero] at h0
  have h2b : (2 : ℂ) * b ≠ 0 := mul_ne_zero two_ne_zero hb0
  rw [div_mul_eq_mul_div, h0, mul_div_cancel_left₀ _ h2b]

lemma gauss_mom0 (b : ℂ) (hb : 0 < b.re) :
    ∫ x : ℝ, (x : ℂ) ^ 0 * Complex.exp (-b * (x : ℂ) ^ 2) = (π / b) ^ (1 / 2 : ℂ) := by
  simp only [pow_zero, one_mul]; exact integral_gaussian_complex hb

lemma gauss_mom2 (b : ℂ) (hb : 0 < b.re) :
    ∫ x : ℝ, (x : ℂ) ^ 2 * Complex.exp (-b * (x : ℂ) ^ 2) = (π / b) ^ (1 / 2 : ℂ) / (2 * b) := by
  rw [gauss_rec 0 b hb, gauss_mom0 b hb]; push_cast; ring

lemma gauss_mom4 (b : ℂ) (hb : 0 < b.re) :
    ∫ x : ℝ, (x : ℂ) ^ 4 * Complex.exp (-b * (x : ℂ) ^ 2) =
      3 * (π / b) ^ (1 / 2 : ℂ) / (4 * b ^ 2) := by
  rw [gauss_rec 2 b hb, gauss_mom2 b hb]
  have hb0 : b ≠ 0 := by intro h; rw [h] at hb; simp at hb
  push_cast; field_simp; ring

lemma gauss_rec_real (m : ℕ) (c : ℝ) (hc : 0 < c) :
    ∫ x : ℝ, x ^ (m + 2) * Real.exp (-c * x ^ 2) =
      (m + 1) / (2 * c) * ∫ x : ℝ, x ^ m * Real.exp (-c * x ^ 2) := by
  have h := gauss_rec m (c : ℂ) (by simpa using hc)
  have e : ∀ k : ℕ, ∫ x : ℝ, (x : ℂ) ^ k * Complex.exp (-(c : ℂ) * (x : ℂ) ^ 2) =
      ((∫ x : ℝ, x ^ k * Real.exp (-c * x ^ 2) : ℝ) : ℂ) := by
    intro k
    rw [show ((∫ x : ℝ, x ^ k * Real.exp (-c * x ^ 2) : ℝ) : ℂ) =
      ∫ x : ℝ, ((x ^ k * Real.exp (-c * x ^ 2) : ℝ) : ℂ) from integral_ofReal.symm]
    congr 1; funext x; push_cast; rfl
  rw [e, e] at h
  exact_mod_cast h

/-- `∫_ℝ σ^{2j} e^{-cσ²} = √(π/c) (2j-1)!!/(2c)^j`, as the list of the needed values. -/
lemma gauss_real_vals (c : ℝ) (hc : 0 < c) :
    (∫ x : ℝ, x ^ 0 * Real.exp (-c * x ^ 2)) = √(π / c) ∧
    (∫ x : ℝ, x ^ 2 * Real.exp (-c * x ^ 2)) = √(π / c) / (2 * c) ∧
    (∫ x : ℝ, x ^ 4 * Real.exp (-c * x ^ 2)) = 3 * √(π / c) / (2 * c) ^ 2 ∧
    (∫ x : ℝ, x ^ 6 * Real.exp (-c * x ^ 2)) = 15 * √(π / c) / (2 * c) ^ 3 ∧
    (∫ x : ℝ, x ^ 8 * Real.exp (-c * x ^ 2)) = 105 * √(π / c) / (2 * c) ^ 4 ∧
    (∫ x : ℝ, x ^ 10 * Real.exp (-c * x ^ 2)) = 945 * √(π / c) / (2 * c) ^ 5 := by
  have h0 : (∫ x : ℝ, x ^ 0 * Real.exp (-c * x ^ 2)) = √(π / c) := by
    simp only [pow_zero, one_mul]; exact integral_gaussian c
  have h2 := gauss_rec_real 0 c hc
  have h4 := gauss_rec_real 2 c hc
  have h6 := gauss_rec_real 4 c hc
  have h8 := gauss_rec_real 6 c hc
  have h10 := gauss_rec_real 8 c hc
  simp only [zero_add] at h2
  rw [h0] at h2
  rw [h2] at h4; rw [h4] at h6; rw [h6] at h8; rw [h8] at h10
  norm_num at h2 h4 h6 h8 h10
  refine ⟨h0, ?_, ?_, ?_, ?_, ?_⟩
  · simp only [neg_mul]; rw [h2]; field_simp
  · simp only [neg_mul]; rw [h4]; field_simp
  · simp only [neg_mul]; rw [h6]; field_simp; ring_nf
  · simp only [neg_mul]; rw [h8]; field_simp; ring_nf
  · simp only [neg_mul]; rw [h10]; field_simp; ring_nf

/-! ## 8. `H` near infinity (g5:lem:H) -/

lemma exp_neg_le2 (w : ℝ) (hw : 0 ≤ w) : Real.exp (-w) ≤ 1 - w + w ^ 2 / 2 := by
  let D : ℝ → ℝ := fun w => 1 - w + w ^ 2 / 2 - Real.exp (-w)
  have hd : ∀ x : ℝ, HasDerivAt D (-1 + x + Real.exp (-x)) x := by
    intro x
    have h1 := ((hasDerivAt_id' x).neg.const_add 1).add ((hasDerivAt_pow 2 x).div_const 2)
    have h2 := (hasDerivAt_id' x).neg.exp
    convert h1.sub h2 using 1
    norm_num
  have hmono : MonotoneOn D (Ici 0) := by
    apply monotoneOn_of_deriv_nonneg (convex_Ici 0)
    · exact fun x _ => (hd x).continuousAt.continuousWithinAt
    · intro x _; exact (hd x).differentiableAt.differentiableWithinAt
    · intro x _
      rw [(hd x).deriv]
      linarith [Real.add_one_le_exp (-x)]
  have := hmono (by simp : (0 : ℝ) ∈ Ici 0) hw hw
  simp only [D] at this
  norm_num at this
  linarith

lemma exp_neg_ge3 (w : ℝ) (hw : 0 ≤ w) : 1 - w + w ^ 2 / 2 - w ^ 3 / 6 ≤ Real.exp (-w) := by
  let D : ℝ → ℝ := fun w => Real.exp (-w) - (1 - w + w ^ 2 / 2 - w ^ 3 / 6)
  have hd : ∀ x : ℝ, HasDerivAt D (-Real.exp (-x) + 1 - x + x ^ 2 / 2) x := by
    intro x
    have h1 := ((((hasDerivAt_id' x).neg.const_add 1).add ((hasDerivAt_pow 2 x).div_const 2)).sub
      ((hasDerivAt_pow 3 x).div_const 6))
    have h2 := (hasDerivAt_id' x).neg.exp
    convert h2.sub h1 using 1
    norm_num; ring
  have hmono : MonotoneOn D (Ici 0) := by
    apply monotoneOn_of_deriv_nonneg (convex_Ici 0)
    · exact fun x _ => (hd x).continuousAt.continuousWithinAt
    · intro x _; exact (hd x).differentiableAt.differentiableWithinAt
    · intro x hx
      rw [interior_Ici] at hx
      rw [(hd x).deriv]
      linarith [exp_neg_le2 x (le_of_lt hx)]
  have := hmono (by simp : (0 : ℝ) ∈ Ici 0) hw hw
  simp only [D] at this
  norm_num at this
  linarith

lemma Hint_model_pt (u : ℝ) (hu : 0 < u) :
    |√(1 - Real.exp (-2 * u)) - √2 * √u * (1 - u / 2)| ≤ 3 / 2 * √2 * (u ^ 2 * √u) := by
  have hsu : 0 < √u := Real.sqrt_pos.2 hu
  have hsu2 : √u ^ 2 = u := Real.sq_sqrt hu.le
  have h22 : √2 ^ 2 = 2 := Real.sq_sqrt (by norm_num)
  have hs2 := sqrt2_pos
  set A := √(1 - Real.exp (-2 * u))
  have hA0 : 0 ≤ A := Real.sqrt_nonneg _
  have he1 : Real.exp (-2 * u) ≤ 1 := Real.exp_le_one_iff.2 (by linarith)
  have hA2 : A ^ 2 = 1 - Real.exp (-2 * u) := Real.sq_sqrt (by linarith)
  have hA1 : A ≤ 1 := by
    rw [Real.sqrt_le_one]; linarith [Real.exp_pos (-2 * u)]
  rcases le_or_gt u 1 with h | h
  · set B := √2 * √u * (1 - u / 2)
    have hB : √2 * √u / 2 ≤ B := by
      simp only [B]
      have : (1 : ℝ) / 2 ≤ 1 - u / 2 := by linarith
      have := mul_le_mul_of_nonneg_left this (by positivity : 0 ≤ √2 * √u)
      linarith
    have hB2 : B ^ 2 = 2 * u * (1 - u / 2) ^ 2 := by
      simp only [B]; rw [mul_pow, mul_pow, h22, hsu2]
    have e1 := exp_neg_le2 (2 * u) (by linarith)
    have e2 := exp_neg_ge3 (2 * u) (by linarith)
    rw [show -(2 * u) = -2 * u by ring] at e1 e2
    have hsq : |A ^ 2 - B ^ 2| ≤ 5 / 6 * u ^ 3 := by
      rw [hA2, hB2, abs_le]; constructor <;> nlinarith
    have hprod : |A - B| * (A + B) ≤ 5 / 6 * u ^ 3 := by
      have : |A - B| * (A + B) = |A ^ 2 - B ^ 2| := by
        rw [show A ^ 2 - B ^ 2 = (A - B) * (A + B) by ring, abs_mul,
          abs_of_nonneg (by have : 0 ≤ √2 * √u / 2 := by positivity
                            linarith : 0 ≤ A + B)]
      rw [this]; exact hsq
    have hAB : √2 * √u / 2 ≤ A + B := by linarith
    have h3 : |A - B| * (√2 * √u / 2) ≤ 5 / 6 * u ^ 3 :=
      le_trans (mul_le_mul_of_nonneg_left hAB (abs_nonneg _)) hprod
    have hu3 : u ^ 3 = (u ^ 2 * √u) * √u := by
      rw [mul_assoc, ← sq, hsu2]; ring
    rw [hu3] at h3
    have h4 : |A - B| * √2 ≤ 5 / 3 * (u ^ 2 * √u) := by
      have : (|A - B| * √2 - 5 / 3 * (u ^ 2 * √u)) * √u ≤ 0 := by nlinarith
      have := (mul_nonpos_iff.1 this)
      rcases this with ⟨_, h2⟩ | ⟨h1, _⟩
      · linarith
      · linarith
    have h5 : |A - B| ≤ 5 / 6 * √2 * (u ^ 2 * √u) := by
      have : |A - B| * √2 * √2 ≤ 5 / 3 * (u ^ 2 * √u) * √2 :=
        mul_le_mul_of_nonneg_right h4 hs2.le
      rw [mul_assoc, ← sq, h22] at this
      nlinarith
    have : 0 ≤ u ^ 2 * √u := by positivity
    nlinarith
  · have h1u : 1 ≤ √u := by rw [Real.le_sqrt (by norm_num) hu.le]; linarith
    have hB : |√2 * √u * (1 - u / 2)| ≤ √2 * √u * (u / 2) := by
      rw [abs_mul, abs_of_nonneg (by positivity : 0 ≤ √2 * √u)]
      apply mul_le_mul_of_nonneg_left _ (by positivity)
      rw [abs_le]; constructor <;> linarith
    have h2 : |A - √2 * √u * (1 - u / 2)| ≤ A + √2 * √u * (u / 2) := by
      refine (abs_sub _ _).trans ?_
      rw [abs_of_nonneg hA0]; linarith
    have h3 : √2 * √u * (u / 2) ≤ 1 / 2 * √2 * (u ^ 2 * √u) := by
      have : u ≤ u ^ 2 := by nlinarith
      have := mul_le_mul_of_nonneg_left this (by positivity : 0 ≤ √2 * √u)
      nlinarith
    have h4 : 1 ≤ √2 * (u ^ 2 * √u) := by
      have hs21 : 1 ≤ √2 := by rw [Real.le_sqrt (by norm_num) (by norm_num)]; norm_num
      have : 1 ≤ u ^ 2 * √u := by nlinarith
      nlinarith
    nlinarith

lemma rpow_half_pow (j : ℕ) (u : ℝ) (hu : 0 < u) : u ^ ((j : ℝ) + 1 / 2) = u ^ j * √u := by
  rw [Real.rpow_add hu, Real.rpow_natCast, ← Real.sqrt_eq_rpow]

lemma integrableOn_pow_sqrt (j : ℕ) (x : ℝ) (hx : 0 < x) :
    IntegrableOn (fun u : ℝ => u ^ j * √u * Real.exp (-(x * u))) (Ioi 0) := by
  have h := integrableOn_rpow_mul_exp_neg_mul_rpow (p := 1) (s := (j : ℝ) + 1 / 2) (b := x)
    (by have := j.cast_nonneg (α := ℝ); linarith) one_pos hx
  refine h.congr_fun (fun u hu => ?_) measurableSet_Ioi
  simp only [Real.rpow_one]
  rw [rpow_half_pow j u hu]; ring_nf

lemma integral_pow_sqrt (j : ℕ) (x : ℝ) (hx : 0 < x) :
    ∫ u in Ioi (0 : ℝ), u ^ j * √u * Real.exp (-(x * u)) =
      (1 / x) ^ ((j : ℝ) + 3 / 2) * Real.Gamma ((j : ℝ) + 3 / 2) := by
  rw [← Real.integral_rpow_mul_exp_neg_mul_Ioi (by positivity) hx]
  refine setIntegral_congr_fun measurableSet_Ioi (fun u hu => ?_)
  rw [show (j : ℝ) + 3 / 2 - 1 = j + 1 / 2 by ring, rpow_half_pow j u hu]

lemma even_integral (g : ℝ → ℂ) (hg : Integrable g) (he : ∀ x, g (-x) = g x) :
    ∫ x, g x = 2 * ∫ x in Ioi (0 : ℝ), g x := by
  rw [split_Iic_Ioi g hg 0]
  have : ∫ x in Iic (0 : ℝ), g x = ∫ x in Ioi (0 : ℝ), g x := by
    have h := integral_comp_neg_Ioi 0 g
    rw [neg_zero] at h
    rw [← h]
    exact setIntegral_congr_fun measurableSet_Ioi (fun x _ => he x)
  rw [this]; ring

/-- `G_j(u) = u^j √u e^{-su}` as a complex function. -/
noncomputable def Gj (j : ℕ) (s : ℂ) (u : ℝ) : ℂ := ((u ^ j * √u : ℝ) : ℂ) * Complex.exp (-s * u)

lemma Gj_int (j : ℕ) (s : ℂ) (hs : 0 < s.re) : IntegrableOn (Gj j s) (Ioi 0) := by
  refine (integrableOn_pow_sqrt j s.re hs).norm.mono' ?_ ?_
  · refine (Continuous.aestronglyMeasurable ?_)
    unfold Gj; fun_prop
  · refine (ae_restrict_iff' measurableSet_Ioi).2 (Filter.Eventually.of_forall fun u hu => ?_)
    have hu' : 0 < u := hu
    simp only [Gj, norm_mul, Complex.norm_real, Complex.norm_exp, Real.norm_eq_abs, Real.abs_exp]
    apply le_of_eq; congr 2; simp

lemma Gj_eq (j : ℕ) (s : ℂ) (hs : 0 < s.re) :
    ∫ u in Ioi (0 : ℝ), Gj j s u =
      ∫ x : ℝ, (x : ℂ) ^ (2 * j + 2) * Complex.exp (-s * (x : ℂ) ^ 2) := by
  rw [even_integral _ (integrable_cpow_gauss (2 * j + 2) s hs) (fun x => by
    simp only [Complex.ofReal_neg]
    rw [Even.neg_pow (by exact ⟨j + 1, by ring⟩), neg_sq])]
  rw [← integral_comp_rpow_Ioi_of_pos (g := Gj j s) (p := 2) two_pos, ← integral_const_mul]
  refine setIntegral_congr_fun measurableSet_Ioi (fun x hx => ?_)
  have hx' : 0 < x := hx
  simp only [Gj]
  rw [show (2 : ℝ) - 1 = 1 by norm_num, Real.rpow_one, Real.rpow_two, Real.sqrt_sq hx'.le,
    Complex.real_smul]
  push_cast
  ring

lemma G0_val (s : ℂ) (hs : 0 < s.re) :
    ∫ u in Ioi (0 : ℝ), Gj 0 s u = (π / s) ^ (1 / 2 : ℂ) / (2 * s) := by
  rw [Gj_eq 0 s hs]; simp only [mul_zero, zero_add]; rw [gauss_mom2 s hs]

lemma G1_val (s : ℂ) (hs : 0 < s.re) :
    ∫ u in Ioi (0 : ℝ), Gj 1 s u = 3 * (π / s) ^ (1 / 2 : ℂ) / (4 * s ^ 2) := by
  rw [Gj_eq 1 s hs, show 2 * 1 + 2 = 4 by norm_num, gauss_mom4 s hs]

lemma Gamma_seven_halves : Real.Gamma ((2 : ℕ) + 3 / 2) = 15 * √π / 8 := by
  rw [show ((2 : ℕ) : ℝ) + 3 / 2 = (1 / 2 + 1 + 1) + 1 by norm_num,
    Real.Gamma_add_one (by norm_num), Real.Gamma_add_one (by norm_num),
    Real.Gamma_add_one (by norm_num), Real.Gamma_one_half_eq]
  ring

/-- g5:lem:H: `|H(s) - (π/s)^{1/2}/(√2 s) (1 - 3/(4s))| ≤ (45/16)√(2π) (Re s)^{-7/2}`. -/
lemma H_approx (s : ℂ) (hs : 0 < s.re) :
    ‖H s - (π / s) ^ (1 / 2 : ℂ) / (√2 * s) * (1 - 3 / (4 * s))‖ ≤
      45 / 16 * √(2 * π) * (1 / s.re) ^ ((2 : ℕ) + 3 / 2 : ℝ) := by
  have hs0 : s ≠ 0 := by intro h; rw [h] at hs; simp at hs
  obtain ⟨hI, hH⟩ := H_rep s hs
  have hG0 := Gj_int 0 s hs
  have hG1 := Gj_int 1 s hs
  set Mint : ℝ → ℂ := fun u => (√2 : ℂ) * Gj 0 s u - ((√2 / 2 : ℝ) : ℂ) * Gj 1 s u
  have hMi : IntegrableOn Mint (Ioi 0) := (hG0.const_mul _).sub (hG1.const_mul _)
  have hMv : ∫ u in Ioi (0 : ℝ), Mint u = (π / s) ^ (1 / 2 : ℂ) / (√2 * s) * (1 - 3 / (4 * s)) := by
    simp only [Mint]
    rw [integral_sub (hG0.const_mul _) (hG1.const_mul _), integral_const_mul, integral_const_mul,
      G0_val s hs, G1_val s hs]
    have h2 : (√2 : ℂ) ≠ 0 := by exact_mod_cast sqrt2_pos.ne'
    have h22 : (√2 : ℂ) * √2 = 2 := by
      rw [← Complex.ofReal_mul, Real.mul_self_sqrt (by norm_num)]; norm_num
    push_cast
    field_simp
    linear_combination (π / s) ^ (1 / 2 : ℂ) * (4 * s - 3) * h22
  have hb : ∀ u ∈ Ioi (0 : ℝ), ‖Hint s u - Mint u‖ ≤
      3 / 2 * √2 * (u ^ 2 * √u * Real.exp (-(s.re * u))) := by
    intro u hu
    have hu' : 0 < u := hu
    have e : Hint s u - Mint u =
        ((√(1 - Real.exp (-2 * u)) - √2 * √u * (1 - u / 2) : ℝ) : ℂ) * Complex.exp (-s * u) := by
      simp only [Hint, Mint, Gj]; push_cast; ring
    rw [e, norm_mul, Complex.norm_real, Complex.norm_exp, Real.norm_eq_abs]
    have : (-s * (u : ℂ)).re = -(s.re * u) := by simp
    rw [this]
    have := Hint_model_pt u hu'
    have he := Real.exp_pos (-(s.re * u))
    nlinarith
  have hbi : IntegrableOn (fun u : ℝ => 3 / 2 * √2 * (u ^ 2 * √u * Real.exp (-(s.re * u))))
      (Ioi 0) := (integrableOn_pow_sqrt 2 s.re hs).const_mul _
  rw [hH, ← hMv, ← integral_sub hI hMi]
  refine (norm_integral_le_of_norm_le hbi
    ((ae_restrict_iff' measurableSet_Ioi).2 (Filter.Eventually.of_forall hb))).trans (le_of_eq ?_)
  rw [integral_const_mul, integral_pow_sqrt 2 s.re hs, Gamma_seven_halves]
  rw [Real.sqrt_mul (by norm_num)]
  ring

/-! ## 9. Taylor expansion of an entire function with Cauchy bounds -/

/-- The Taylor coefficient `f^{(k)}(0)/k!`. -/
noncomputable def tc (f : ℂ → ℂ) (k : ℕ) : ℂ := ((k.factorial : ℂ))⁻¹ * iteratedDeriv k f 0

lemma tc_le (f : ℂ → ℂ) (hf : Differentiable ℂ f) (M : ℝ)
    (hM : ∀ w ∈ Metric.sphere (0 : ℂ) (1 / 2), ‖f w‖ ≤ M) (k : ℕ) :
    ‖tc f k‖ ≤ M * 2 ^ k := by
  have h := Complex.norm_iteratedDeriv_le_of_forall_mem_sphere_norm_le k (by norm_num : (0 : ℝ) < 1 / 2)
    hf.diffContOnCl hM
  unfold tc
  rw [norm_mul, norm_inv, Complex.norm_natCast]
  have hk : (0 : ℝ) < k.factorial := by exact_mod_cast Nat.factorial_pos k
  rw [inv_mul_le_iff₀ hk]
  refine h.trans (le_of_eq ?_)
  rw [one_div, inv_pow, div_inv_eq_mul]; ring

lemma taylor4 (f : ℂ → ℂ) (hf : Differentiable ℂ f) (M : ℝ)
    (hM : ∀ w ∈ Metric.sphere (0 : ℂ) (1 / 2), ‖f w‖ ≤ M) (w : ℂ) (hw : ‖w‖ ≤ 1 / 4) :
    ‖f w - (tc f 0 + tc f 1 * w + tc f 2 * w ^ 2 + tc f 3 * w ^ 3)‖ ≤ 32 * M * ‖w‖ ^ 4 := by
  have hS := Complex.hasSum_taylorSeries_of_entire hf 0 w
  simp only [sub_zero] at hS
  have h4 := (hasSum_nat_add_iff' 4).2 hS
  have hsum : ∑ i ∈ Finset.range 4, ((i.factorial : ℂ))⁻¹ • w ^ i • iteratedDeriv i f 0 =
      tc f 0 + tc f 1 * w + tc f 2 * w ^ 2 + tc f 3 * w ^ 3 := by
    simp [Finset.sum_range_succ, tc, smul_eq_mul]; ring
  rw [hsum] at h4
  have hM0 : 0 ≤ M := by
    have := hM ((1 / 2 : ℂ)) (by simp); exact (norm_nonneg _).trans this
  set r := 2 * ‖w‖
  have hr0 : 0 ≤ r := by positivity
  have hr1 : r < 1 := by simp only [r]; linarith
  have hg : HasSum (fun n : ℕ => M * r ^ 4 * r ^ n) (M * r ^ 4 * (1 - r)⁻¹) :=
    (hasSum_geometric_of_lt_one hr0 hr1).mul_left _
  have hb : ∀ n : ℕ, ‖((n + 4).factorial : ℂ)⁻¹ • w ^ (n + 4) • iteratedDeriv (n + 4) f 0‖ ≤
      M * r ^ 4 * r ^ n := by
    intro n
    have := tc_le f hf M hM (n + 4)
    unfold tc at this
    rw [smul_eq_mul, smul_eq_mul, ← mul_assoc, mul_comm (((n + 4).factorial : ℂ)⁻¹), mul_assoc,
      norm_mul, norm_pow]
    calc ‖w‖ ^ (n + 4) * ‖((n + 4).factorial : ℂ)⁻¹ * iteratedDeriv (n + 4) f 0‖
        ≤ ‖w‖ ^ (n + 4) * (M * 2 ^ (n + 4)) := by gcongr
      _ = M * r ^ 4 * r ^ n := by simp only [r]; ring
  have := h4.norm_le_of_bounded hg hb
  refine this.trans ?_
  have h1 : (1 - r)⁻¹ ≤ 2 := by
    rw [inv_le_comm₀ (by linarith) (by norm_num)]; simp only [r]; linarith
  calc M * r ^ 4 * (1 - r)⁻¹ ≤ M * r ^ 4 * 2 := by gcongr
    _ = 32 * M * ‖w‖ ^ 4 := by simp only [r]; ring

/-! ## 10. The model amplitude `P = φ + y χ` -/

/-- `ω̄ = e^{-iπ/4}`. -/
noncomputable def ωb : ℂ := Complex.exp (-(π / 4) * Complex.I)

/-- `ψ(w) = -3w/2 - (iπ/4) e^{-2w}`. -/
noncomputable def ψf (w : ℂ) : ℂ := -(3 / 2) * w - Complex.I * π / 4 * Complex.exp (-2 * w)

/-- `φ = e^ψ`. -/
noncomputable def φf (w : ℂ) : ℂ := Complex.exp (ψf w)

/-- `χ(w) = -((3/4) ω̄ e^{-w} + (π/6) ω̄³ e^{-3w}) φ(w)`. -/
noncomputable def χf (w : ℂ) : ℂ :=
  -(3 / 4 * ωb * Complex.exp (-w) + π / 6 * ωb ^ 3 * Complex.exp (-3 * w)) * φf w

/-- `ν₀ = ((-3/2 + iπ/2)² - iπ)/2`. -/
noncomputable def ν0 : ℂ := ((-(3 / 2) + Complex.I * π / 2) ^ 2 - Complex.I * π) / 2

lemma φf_diff : Differentiable ℂ φf := by unfold φf ψf; fun_prop

lemma χf_diff : Differentiable ℂ χf := by
  unfold χf; have := φf_diff; fun_prop

lemma norm_ωb : ‖ωb‖ = 1 := by
  rw [ωb, Complex.norm_exp]; simp

/-- `M_φ = e^{3/4 + πe/4}`. -/
noncomputable def Mφ : ℝ := Real.exp (3 / 4 + π * Real.exp 1 / 4)

lemma sphere_re (w : ℂ) (hw : w ∈ Metric.sphere (0 : ℂ) (1 / 2)) : |w.re| ≤ 1 / 2 := by
  have h : ‖w‖ = 1 / 2 := by simpa using hw
  rw [← h]; exact Complex.abs_re_le_norm w

lemma φf_sphere (w : ℂ) (hw : w ∈ Metric.sphere (0 : ℂ) (1 / 2)) : ‖φf w‖ ≤ Mφ := by
  have hre := abs_le.1 (sphere_re w hw)
  unfold φf Mφ
  rw [Complex.norm_exp]
  apply Real.exp_le_exp.2
  have h1 : (ψf w).re = -(3 / 2) * w.re + π / 4 * (Complex.exp (-2 * w)).im := by
    simp [ψf, Complex.mul_re, Complex.mul_im]
  have h2 : (Complex.exp (-2 * w)).im ≤ Real.exp 1 := by
    refine (Complex.im_le_norm _).trans ?_
    rw [Complex.norm_exp]; apply Real.exp_le_exp.2; simp; linarith
  rw [h1]
  have := Real.pi_pos
  nlinarith

lemma norm_cexp_neg_le (c : ℝ) (hc : 0 ≤ c) (w : ℂ) (hw : w ∈ Metric.sphere (0 : ℂ) (1 / 2)) :
    ‖Complex.exp (-c * w)‖ ≤ Real.exp (c / 2) := by
  have hre := abs_le.1 (sphere_re w hw)
  rw [Complex.norm_exp]; apply Real.exp_le_exp.2
  simp only [Complex.mul_re, Complex.neg_re, Complex.ofReal_re, Complex.neg_im, Complex.ofReal_im]
  nlinarith

/-- `M_χ = ((3/4)e^{1/2} + (π/6)e^{3/2}) M_φ`. -/
noncomputable def Mχ : ℝ := (3 / 4 * Real.exp (1 / 2) + π / 6 * Real.exp (3 / 2)) * Mφ

lemma χf_sphere (w : ℂ) (hw : w ∈ Metric.sphere (0 : ℂ) (1 / 2)) : ‖χf w‖ ≤ Mχ := by
  unfold χf Mχ
  rw [norm_mul, norm_neg]
  have h1 : ‖Complex.exp (-w)‖ ≤ Real.exp (1 / 2) := by
    simpa using norm_cexp_neg_le 1 zero_le_one w hw
  have h3 : ‖Complex.exp (-(3 * w))‖ ≤ Real.exp (3 / 2) := by
    simpa using norm_cexp_neg_le 3 (by norm_num) w hw
  have hφ := φf_sphere w hw
  have hb : ‖3 / 4 * ωb * Complex.exp (-w) + π / 6 * ωb ^ 3 * Complex.exp (-(3 * w))‖ ≤
      3 / 4 * Real.exp (1 / 2) + π / 6 * Real.exp (3 / 2) := by
    refine (norm_add_le _ _).trans (add_le_add ?_ ?_)
    · rw [norm_mul, norm_mul, norm_ωb, mul_one]
      have : ‖(3 / 4 : ℂ)‖ = 3 / 4 := by norm_num
      rw [this]; gcongr
    · rw [norm_mul, norm_mul, norm_pow, norm_ωb, one_pow, mul_one]
      have : ‖((π : ℂ) / 6)‖ = π / 6 := by
        rw [norm_div, Complex.norm_real, Real.norm_eq_abs, abs_of_pos Real.pi_pos]; norm_num
      rw [this]; gcongr
  rw [show -(3 : ℂ) * w = -(3 * w) by ring]
  exact mul_le_mul hb hφ (norm_nonneg _) (by positivity)

lemma tc_zero (f : ℂ → ℂ) : tc f 0 = f 0 := by simp [tc]

lemma φf_zero : φf 0 = Complex.exp (-(Complex.I * π / 4)) := by
  simp [φf, ψf]

lemma χf_zero : χf 0 = -(3 / 4 * ωb + π / 6 * ωb ^ 3) * Complex.exp (-(Complex.I * π / 4)) := by
  simp [χf, φf_zero]

/-- `ψ'`. -/
noncomputable def ψd (w : ℂ) : ℂ := -(3 / 2) + Complex.I * π / 2 * Complex.exp (-2 * w)

lemma ψf_hasDeriv (w : ℂ) : HasDerivAt ψf (ψd w) w := by
  have h1 : HasDerivAt (fun w : ℂ => -2 * w) (-2) w := by
    simpa using (hasDerivAt_id w).const_mul (-2 : ℂ)
  have h2 := (h1.cexp).const_mul (Complex.I * π / 4)
  have h3 : HasDerivAt (fun w : ℂ => -(3 / 2 : ℂ) * w - Complex.I * π / 4 * Complex.exp (-2 * w))
      _ w := ((hasDerivAt_id' w).const_mul (-(3 / 2) : ℂ)).sub h2
  unfold ψf ψd
  exact h3.congr_deriv (by ring)

lemma ψd_hasDeriv (w : ℂ) : HasDerivAt ψd (-(Complex.I * π) * Complex.exp (-2 * w)) w := by
  have h1 : HasDerivAt (fun w : ℂ => -2 * w) (-2) w := by
    simpa using (hasDerivAt_id w).const_mul (-2 : ℂ)
  have h2 := ((h1.cexp).const_mul (Complex.I * π / 2)).const_add (-(3 / 2) : ℂ)
  unfold ψd
  convert h2 using 1
  ring

lemma φf_hasDeriv (w : ℂ) : HasDerivAt φf (φf w * ψd w) w := by
  have := (ψf_hasDeriv w).cexp
  unfold φf; convert this using 1

lemma deriv_φf : deriv φf = fun w => φf w * ψd w := by
  funext w; exact (φf_hasDeriv w).deriv

lemma tc_φf_two : tc φf 2 = ν0 * Complex.exp (-(Complex.I * π / 4)) := by
  have hd : HasDerivAt (fun w => φf w * ψd w)
      (φf 0 * ψd 0 * ψd 0 + φf 0 * (-(Complex.I * π) * Complex.exp (-2 * 0))) 0 :=
    (φf_hasDeriv 0).mul (ψd_hasDeriv 0)
  have h2 : iteratedDeriv 2 φf 0 = φf 0 * ψd 0 * ψd 0 + φf 0 * (-(Complex.I * π) *
      Complex.exp (-2 * 0)) := by
    rw [iteratedDeriv_succ, iteratedDeriv_one, deriv_φf]
    exact hd.deriv
  unfold tc
  rw [h2, φf_zero]
  simp only [ψd, mul_zero, Complex.exp_zero, mul_one, ν0]
  norm_num
  ring

/-! ## 11. The cosh factor in the window -/

/-- `q(σ) = cosh σ - 1 - σ²/2`. -/
noncomputable def qf (σ : ℝ) : ℝ := (Real.exp σ + Real.exp (-σ)) / 2 - 1 - σ ^ 2 / 2

lemma q_six (σ : ℝ) (hσ : |σ| ≤ 1 / 4) : |qf σ - σ ^ 4 / 24| ≤ 7 / 4320 * σ ^ 6 := by
  have h1 := Real.exp_bound (x := σ) (by linarith) (n := 6) (by norm_num)
  have h2 := Real.exp_bound (x := -σ) (by rw [abs_neg]; linarith) (n := 6) (by norm_num)
  simp [Finset.sum_range_succ, Nat.factorial] at h1 h2
  have e6 : |σ| ^ 6 = σ ^ 6 := by rw [show (6 : ℕ) = 2 * 3 by rfl, pow_mul, sq_abs, ← pow_mul]
  rw [e6] at h1 h2
  rw [abs_le] at h1 h2 ⊢
  unfold qf
  constructor <;> nlinarith [h1.1, h1.2, h2.1, h2.2]

lemma q_nonneg (σ : ℝ) (hσ : |σ| ≤ 1 / 4) : 0 ≤ qf σ := by
  have h := abs_le.1 (q_six σ hσ)
  have hs : σ ^ 2 ≤ 1 / 16 := by
    have := sq_abs σ; nlinarith [abs_nonneg σ]
  have : 7 / 4320 * σ ^ 6 ≤ σ ^ 4 / 24 := by
    have h4 : 0 ≤ σ ^ 4 := by positivity
    have : σ ^ 6 = σ ^ 4 * σ ^ 2 := by ring
    rw [this]; nlinarith
  linarith [h.1]

lemma q_four (σ : ℝ) (hσ : |σ| ≤ 1 / 4) : qf σ ≤ 0.0418 * σ ^ 4 := by
  have h := abs_le.1 (q_six σ hσ)
  have hs : σ ^ 2 ≤ 1 / 16 := by
    have := sq_abs σ; nlinarith [abs_nonneg σ]
  have h4 : 0 ≤ σ ^ 4 := by positivity
  have : σ ^ 6 = σ ^ 4 * σ ^ 2 := by ring
  rw [this] at h
  nlinarith [h.2]

/-- `z = θ(1 - i)`. -/
noncomputable def zv (θ : ℝ) : ℂ := (θ : ℂ) * (1 - Complex.I)

lemma norm_zv (θ : ℝ) (hθ : 0 ≤ θ) : ‖zv θ‖ = √2 * θ := by
  unfold zv
  rw [norm_mul, Complex.norm_real, Real.norm_eq_abs, abs_of_nonneg hθ, Complex.norm_eq_sqrt_sq_add_sq]
  simp; ring_nf

lemma zv_re (θ : ℝ) : (zv θ).re = θ := by simp [zv]

/-- `K(σ) = e^{-z(cosh σ - 1)}`. -/
noncomputable def Kf (θ σ : ℝ) : ℂ :=
  Complex.exp (-zv θ * (((Real.exp σ + Real.exp (-σ)) / 2 - 1 : ℝ) : ℂ))

/-- `G(σ) = e^{-zσ²/2}`. -/
noncomputable def Gf (θ σ : ℝ) : ℂ := Complex.exp (-zv θ * ((σ ^ 2 / 2 : ℝ) : ℂ))

lemma norm_cexp_mz (θ t : ℝ) : ‖Complex.exp (-zv θ * (t : ℂ))‖ = Real.exp (-(θ * t)) := by
  rw [Complex.norm_exp]; congr 1; simp [zv]

lemma norm_Kf_le (θ σ : ℝ) (hθ : 0 ≤ θ) : ‖Kf θ σ‖ ≤ Real.exp (-(θ * (σ ^ 2 / 2))) := by
  unfold Kf; rw [norm_cexp_mz]; apply Real.exp_le_exp.2
  have := cosh_quad σ; nlinarith

lemma norm_Gf (θ σ : ℝ) : ‖Gf θ σ‖ = Real.exp (-(θ * (σ ^ 2 / 2))) := by
  unfold Gf; rw [norm_cexp_mz]

lemma Kf_eq (θ σ : ℝ) : Kf θ σ = Gf θ σ * Complex.exp (-zv θ * (qf σ : ℂ)) := by
  unfold Kf Gf qf; rw [← Complex.exp_add]; congr 1; push_cast; ring

lemma er_bound (θ σ : ℝ) (hθ : 0 ≤ θ) (hσ : |σ| ≤ 1 / 4) :
    ‖Kf θ σ - Gf θ σ * (1 - zv θ * (σ : ℂ) ^ 4 / 24)‖ ≤
      Real.exp (-(0.4963 * θ * σ ^ 2)) * (0.0023 * θ * σ ^ 6 + 0.0035 * θ ^ 2 * σ ^ 8) := by
  set w := -zv θ * (qf σ : ℂ)
  have hq0 := q_nonneg σ hσ
  have hq4 := q_four σ hσ
  have hq6 := q_six σ hσ
  have hs : σ ^ 2 ≤ 1 / 16 := by
    have := sq_abs σ; nlinarith [abs_nonneg σ]
  have hs2 := sqrt2_lt
  have hs2' := sqrt2_gt
  have hw : ‖w‖ = √2 * θ * qf σ := by
    simp only [w]
    rw [norm_mul, norm_neg, norm_zv θ hθ, Complex.norm_real, Real.norm_eq_abs, abs_of_nonneg hq0]
  have hwle : ‖w‖ ≤ 1.41422 * 0.0418 * θ * σ ^ 4 := by
    rw [hw]
    have := mul_le_mul_of_nonneg_left hq4 (by positivity : 0 ≤ √2 * θ)
    nlinarith [mul_le_mul_of_nonneg_right hs2.le (by positivity : 0 ≤ θ * (0.0418 * σ ^ 4))]
  have e1 : Kf θ σ - Gf θ σ * (1 - zv θ * (σ : ℂ) ^ 4 / 24) =
      Gf θ σ * ((Complex.exp w - (1 + w)) + (w + zv θ * (σ : ℂ) ^ 4 / 24)) := by
    rw [Kf_eq]; ring
  rw [e1, norm_mul, norm_Gf]
  have h2 : ‖Complex.exp w - (1 + w)‖ ≤ ‖w‖ ^ 2 * Real.exp ‖w‖ := by
    have := Complex.norm_exp_sub_sum_le_norm_mul_exp w 2
    simpa [Finset.sum_range_succ] using this
  have h3 : ‖w + zv θ * (σ : ℂ) ^ 4 / 24‖ ≤ √2 * θ * (7 / 4320 * σ ^ 6) := by
    have : w + zv θ * (σ : ℂ) ^ 4 / 24 = -zv θ * ((qf σ - σ ^ 4 / 24 : ℝ) : ℂ) := by
      simp only [w]; push_cast; ring
    rw [this, norm_mul, norm_neg, norm_zv θ hθ, Complex.norm_real, Real.norm_eq_abs]
    exact mul_le_mul_of_nonneg_left hq6 (by positivity)
  have hw0 : 0 ≤ ‖w‖ := norm_nonneg _
  have hws : ‖w‖ ≤ 0.003695 * θ * σ ^ 2 := by
    have h4 : σ ^ 4 ≤ σ ^ 2 / 16 := by
      have : σ ^ 4 = σ ^ 2 * σ ^ 2 := by ring
      rw [this]; nlinarith [sq_nonneg σ]
    nlinarith
  have hexp : Real.exp (-(θ * (σ ^ 2 / 2))) * Real.exp ‖w‖ ≤ Real.exp (-(0.4963 * θ * σ ^ 2)) := by
    rw [← Real.exp_add]; apply Real.exp_le_exp.2; nlinarith [sq_nonneg σ]
  have hexp2 : Real.exp (-(θ * (σ ^ 2 / 2))) ≤ Real.exp (-(0.4963 * θ * σ ^ 2)) := by
    apply Real.exp_le_exp.2; nlinarith [sq_nonneg σ]
  have hw2 : ‖w‖ ^ 2 ≤ 0.0035 * θ ^ 2 * σ ^ 8 := by
    have : ‖w‖ ^ 2 ≤ (1.41422 * 0.0418 * θ * σ ^ 4) ^ 2 := pow_le_pow_left₀ hw0 hwle 2
    nlinarith [sq_nonneg (θ * σ ^ 4)]
  have hE1 := Real.exp_pos (-(θ * (σ ^ 2 / 2)))
  have hE2 := Real.exp_pos ‖w‖
  calc Real.exp (-(θ * (σ ^ 2 / 2))) * ‖Complex.exp w - (1 + w) + (w + zv θ * (σ : ℂ) ^ 4 / 24)‖
      ≤ Real.exp (-(θ * (σ ^ 2 / 2))) * (‖w‖ ^ 2 * Real.exp ‖w‖ + √2 * θ * (7 / 4320 * σ ^ 6)) := by
        gcongr; exact (norm_add_le _ _).trans (add_le_add h2 h3)
    _ = (Real.exp (-(θ * (σ ^ 2 / 2))) * Real.exp ‖w‖) * ‖w‖ ^ 2 +
        Real.exp (-(θ * (σ ^ 2 / 2))) * (√2 * θ * (7 / 4320 * σ ^ 6)) := by ring
    _ ≤ Real.exp (-(0.4963 * θ * σ ^ 2)) * (0.0035 * θ ^ 2 * σ ^ 8) +
        Real.exp (-(0.4963 * θ * σ ^ 2)) * (0.0023 * θ * σ ^ 6) := by
        have h5 : √2 * θ * (7 / 4320 * σ ^ 6) ≤ 0.0023 * θ * σ ^ 6 := by
          have h6 : 0 ≤ θ * σ ^ 6 := by positivity
          nlinarith [mul_le_mul_of_nonneg_right hs2.le h6]
        apply add_le_add
        · exact mul_le_mul hexp hw2 (by positivity) (by positivity)
        · exact mul_le_mul hexp2 h5 (by positivity) (by positivity)
    _ = _ := by ring

/-! ## 12. `f` on the window: the factorisation -/

lemma ωb_eq : ωb = ((√2 / 2 : ℝ) : ℂ) - ((√2 / 2 : ℝ) : ℂ) * Complex.I := by
  rw [ωb, neg_mul, Complex.exp_neg, ← Complex.exp_neg]
  rw [show -((π : ℂ) / 4 * Complex.I) = ((-(π / 4) : ℝ) : ℂ) * Complex.I by push_cast; ring,
    Complex.exp_mul_I]
  rw [← Complex.ofReal_cos, ← Complex.ofReal_sin, Real.cos_neg, Real.sin_neg, Real.cos_pi_div_four,
    Real.sin_pi_div_four]
  push_cast; ring

lemma ω_ωb : ω * ωb = 1 := by
  rw [ω, ωb, ← Complex.exp_add]; simp

lemma ωb_sq : ωb ^ 2 = -Complex.I := by
  rw [ωb_eq]
  have h : ((√2 / 2 : ℝ) : ℂ) ^ 2 = 1 / 2 := by
    rw [← Complex.ofReal_pow, div_pow, Real.sq_sqrt (by norm_num)]; norm_num
  linear_combination (1 - 2 * Complex.I + Complex.I ^ 2) * h + (1 / 2 : ℂ) * Complex.I_sq

/-- `ε = 1/s = ω̄ y e^{-σ}`. -/
noncomputable def εv (N σ : ℝ) : ℂ := ωb * ((yv N * Real.exp (-σ) : ℝ) : ℂ)

lemma s_ε (N σ : ℝ) (hN : 0 < N) : sv N σ * εv N σ = 1 := by
  unfold sv εv
  have h1 : ρv N * Real.exp σ * (yv N * Real.exp (-σ)) = 1 := by
    have he : Real.exp σ * Real.exp (-σ) = 1 := by rw [← Real.exp_add]; simp
    calc ρv N * Real.exp σ * (yv N * Real.exp (-σ)) = (ρv N * yv N) * (Real.exp σ * Real.exp (-σ)) := by ring
      _ = 1 := by rw [ρy N hN, he, one_mul]
  calc ω * ((ρv N * Real.exp σ : ℝ) : ℂ) * (ωb * ((yv N * Real.exp (-σ) : ℝ) : ℂ))
      = (ω * ωb) * ((ρv N * Real.exp σ * (yv N * Real.exp (-σ)) : ℝ) : ℂ) := by push_cast; ring
    _ = 1 := by rw [ω_ωb, h1]; simp

lemma norm_εv (N σ : ℝ) (hN : 0 < N) : ‖εv N σ‖ = yv N * Real.exp (-σ) := by
  unfold εv
  rw [norm_mul, norm_ωb, one_mul, Complex.norm_real, Real.norm_eq_abs,
    abs_of_pos (by have := yv_pos N hN; positivity)]

lemma sv_ne (N σ : ℝ) (hN : 0 < N) : sv N σ ≠ 0 := by
  intro h; have := s_ε N σ hN; rw [h, zero_mul] at this; exact zero_ne_one this

lemma ratio_eq (N σ : ℝ) (hN : 0 < N) : sv N σ / (1 + sv N σ) = (1 + εv N σ)⁻¹ := by
  have hs := sv_ne N σ hN
  have h1 := s_ε N σ hN
  have h2 : 1 + sv N σ ≠ 0 := one_add_ne _ (sv_re_pos N σ hN).le
  have h3 : 1 + εv N σ ≠ 0 := by
    intro h; apply h2
    have : 1 + sv N σ = sv N σ * (1 + εv N σ) := by rw [mul_add, h1]; ring
    rw [this, h, mul_zero]
  field_simp
  linear_combination h1

lemma one_add_ε_ne (N σ : ℝ) (hN : 0 < N) : 1 + εv N σ ≠ 0 := by
  have h2 : 1 + sv N σ ≠ 0 := one_add_ne _ (sv_re_pos N σ hN).le
  have h1 := s_ε N σ hN
  intro h; apply h2
  have : 1 + sv N σ = sv N σ * (1 + εv N σ) := by rw [mul_add, h1]; ring
  rw [this, h, mul_zero]

lemma inv_pow_eq (n : ℕ) (w : ℂ) (hw : w ≠ 0) : (w⁻¹) ^ n = Complex.exp (-(n * Complex.log w)) := by
  rw [Complex.exp_neg, Complex.exp_nat_mul, Complex.exp_log hw, inv_pow]

lemma log_expand (w : ℂ) (hw : ‖w‖ < 1) :
    ‖Complex.log (1 + w) - (w - w ^ 2 / 2 + w ^ 3 / 3)‖ ≤ ‖w‖ ^ 4 * (1 - ‖w‖)⁻¹ / 4 := by
  have h := Complex.norm_log_sub_logTaylor_le 3 hw
  have e : Complex.logTaylor (3 + 1) w = w - w ^ 2 / 2 + w ^ 3 / 3 := by
    simp [Complex.logTaylor, Finset.sum_range_succ]; ring
  rw [e] at h
  norm_num at h ⊢
  exact h

lemma E_sv' (N σ : ℝ) (hN : 0 < N) :
    E (sv N σ) = Complex.exp (-(zv (θv N) / 2 * ((Real.exp σ : ℝ) : ℂ))) := by
  unfold E sv
  congr 1
  rw [ω_eq, zv]
  have h := πρθ N hN
  have h2 : (θv N : ℂ) = ((π * ρv N * √2 / 2 : ℝ) : ℂ) := by
    congr 1; linarith
  rw [h2]
  push_cast
  linear_combination ((π : ℂ) * (√2 : ℂ) * (ρv N : ℂ) * Complex.exp (σ : ℂ) / 4) * Complex.I_sq

lemma nε1 (n : ℕ) (hn : 0 < n) (σ : ℝ) :
    (n : ℂ) * εv n σ = zv (θv n) / 2 * ((Real.exp (-σ) : ℝ) : ℂ) := by
  have hN : (0 : ℝ) < n := by exact_mod_cast hn
  have h := Nyθ n hN
  unfold εv zv
  rw [ωb_eq]
  have h2 : (θv n : ℂ) = ((n * yv n * √2 : ℝ) : ℂ) := by rw [h]
  rw [h2]
  have h22 : (√2 : ℂ) * √2 = 2 := by
    rw [← Complex.ofReal_mul, Real.mul_self_sqrt (by norm_num)]; norm_num
  push_cast
  linear_combination (n * yv n * Real.exp (-σ) * (1 - Complex.I) / 2 : ℂ) * h22 -
    (n * yv n * Real.exp (-σ) * (1 - Complex.I) / 2 : ℂ) * h22

lemma nε2 (n : ℕ) (hn : 0 < n) (σ : ℝ) :
    (n : ℂ) * εv n σ ^ 2 = -(Complex.I * π / 2) * ((Real.exp (-(2 * σ)) : ℝ) : ℂ) := by
  have hN : (0 : ℝ) < n := by exact_mod_cast hn
  have h := Ny2 n hN
  unfold εv
  rw [mul_pow, ωb_sq]
  have he : Real.exp (-σ) ^ 2 = Real.exp (-(2 * σ)) := by
    rw [← Real.exp_nat_mul]; push_cast; ring_nf
  have : ((yv n * Real.exp (-σ) : ℝ) : ℂ) ^ 2 = ((yv n ^ 2 : ℝ) : ℂ) * ((Real.exp (-(2 * σ)) : ℝ) : ℂ) := by
    rw [← he]; push_cast; ring
  rw [this]
  have h2 : (n : ℂ) * ((yv n ^ 2 : ℝ) : ℂ) = ((π / 2 : ℝ) : ℂ) := by
    rw [← h]; push_cast; ring
  push_cast at h2 ⊢
  linear_combination (-Complex.I * Complex.exp (-(2 * (σ : ℂ)))) * h2

lemma nε3 (n : ℕ) (hn : 0 < n) (σ : ℝ) :
    (n : ℂ) * εv n σ ^ 3 = (π / 2 * yv n : ℝ) * ωb ^ 3 * ((Real.exp (-(3 * σ)) : ℝ) : ℂ) := by
  have hN : (0 : ℝ) < n := by exact_mod_cast hn
  have h := Ny2 n hN
  unfold εv
  rw [mul_pow]
  have he : Real.exp (-σ) ^ 3 = Real.exp (-(3 * σ)) := exp_neg_pow3 σ
  have : ((yv n * Real.exp (-σ) : ℝ) : ℂ) ^ 3 = ((yv n ^ 3 : ℝ) : ℂ) * ((Real.exp (-(3 * σ)) : ℝ) : ℂ) := by
    rw [← he]; push_cast; ring
  rw [this]
  have h2 : (n : ℂ) * ((yv n ^ 3 : ℝ) : ℂ) = ((π / 2 * yv n : ℝ) : ℂ) := by
    have : (n : ℝ) * yv n ^ 3 = π / 2 * yv n := by rw [pow_succ, ← mul_assoc, h]
    rw [← this]; push_cast; ring
  linear_combination (ωb ^ 3 * ((Real.exp (-(3 * σ)) : ℝ) : ℂ)) * h2

lemma cpow_half_polar (r φ : ℝ) (hr : 0 < r) (h1 : -π < φ) (h2 : φ ≤ π) :
    ((r : ℂ) * Complex.exp (φ * Complex.I)) ^ (1 / 2 : ℂ) = (√r : ℂ) * Complex.exp (φ / 2 * Complex.I) := by
  have hx : (r : ℂ) * Complex.exp (φ * Complex.I) = Complex.exp ((Real.log r : ℂ) + φ * Complex.I) := by
    rw [Complex.exp_add, ← Complex.ofReal_exp, Real.exp_log hr]
  have hne : (r : ℂ) * Complex.exp (φ * Complex.I) ≠ 0 := by
    rw [hx]; exact Complex.exp_ne_zero _
  rw [Complex.cpow_def_of_ne_zero hne, hx, Complex.log_exp (by simpa using h1) (by simpa using h2)]
  have hs : (√r : ℂ) = Complex.exp ((Real.log r / 2 : ℝ) : ℂ) := by
    rw [← Complex.ofReal_exp]; congr 1
    rw [Real.sqrt_eq_rpow, Real.rpow_def_of_pos hr]; ring_nf
  rw [hs, ← Complex.exp_add]; congr 1; push_cast; ring

lemma sqrt_exp (x : ℝ) : √(Real.exp x) = Real.exp (x / 2) := by
  rw [show Real.exp x = Real.exp (x / 2) ^ 2 by rw [← Real.exp_nat_mul]; ring_nf]
  exact Real.sqrt_sq (Real.exp_pos _).le

/-- `C₀ = √(π/2) y^{3/2} e^{-3iπ/8}`. -/
noncomputable def C0v (N : ℝ) : ℂ :=
  ((√(π / 2) * yv N * √(yv N) : ℝ) : ℂ) * Complex.exp (-(3 * π / 8) * Complex.I)

lemma M0_eq (N σ : ℝ) (hN : 0 < N) :
    (π / sv N σ) ^ (1 / 2 : ℂ) / (√2 * sv N σ) = C0v N * ((Real.exp (-(3 / 2 * σ)) : ℝ) : ℂ) := by
  have hy := yv_pos N hN
  have hs := sv_ne N σ hN
  have hse := s_ε N σ hN
  have hinv : 1 / sv N σ = εv N σ := by
    rw [div_eq_iff hs]; linear_combination -1 * hse
  have h1 : (π : ℂ) / sv N σ = ((π * yv N * Real.exp (-σ) : ℝ) : ℂ) * Complex.exp ((-(π / 4) : ℝ) * Complex.I) := by
    rw [div_eq_mul_one_div, hinv, εv, ωb]; push_cast; ring
  rw [h1, cpow_half_polar _ _ (by positivity) (by linarith [Real.pi_pos]) (by linarith [Real.pi_pos])]
  rw [div_eq_mul_one_div (_ * _), show (1 : ℂ) / (√2 * sv N σ) = (1 / √2) * (1 / sv N σ) by
    field_simp, hinv, εv, ωb, C0v]
  have hsq : √(π * yv N * Real.exp (-σ)) = √π * √(yv N) * Real.exp (-σ / 2) := by
    rw [Real.sqrt_mul (by positivity), Real.sqrt_mul (by positivity), sqrt_exp]
  have hsq2 : √(π / 2) = √π / √2 := Real.sqrt_div' π (by norm_num)
  rw [hsq, hsq2]
  have he : Real.exp (-(3 / 2 * σ)) = Real.exp (-σ / 2) * Real.exp (-σ) := by
    rw [← Real.exp_add]; ring_nf
  rw [he]
  have h2 : (√2 : ℂ) ≠ 0 := by exact_mod_cast sqrt2_pos.ne'
  have hexp : Complex.exp ((-(π / 4) : ℝ) / 2 * Complex.I) * Complex.exp (-(π / 4) * Complex.I) =
      Complex.exp (-(3 * π / 8) * Complex.I) := by
    rw [← Complex.exp_add]; congr 1; push_cast; ring
  push_cast
  field_simp
  have hexp' : Complex.exp (-(↑π * Complex.I / (2 * 4))) * Complex.exp (-(↑π * Complex.I / 4)) =
      Complex.exp (-(↑π * Complex.I * 3 / 8)) := by
    rw [← Complex.exp_add]; congr 1; ring
  linear_combination ((√(yv N) : ℝ) * (yv N : ℝ) : ℂ) * hexp'

/-- `Φ(σ) = e^{-(iπ/4) e^{-2σ}}`. -/
noncomputable def Φv (σ : ℝ) : ℂ := Complex.exp (-(Complex.I * π / 4) * ((Real.exp (-(2 * σ)) : ℝ) : ℂ))

/-- `w₃ = -(πy/6) ω̄³ e^{-3σ}`. -/
noncomputable def w3v (N σ : ℝ) : ℂ := -((π / 6 * yv N : ℝ) : ℂ) * ωb ^ 3 * ((Real.exp (-(3 * σ)) : ℝ) : ℂ)

/-- the remainder of the log series. -/
noncomputable def rv (N σ : ℝ) : ℂ :=
  Complex.log (1 + εv N σ) - (εv N σ - εv N σ ^ 2 / 2 + εv N σ ^ 3 / 3)

/-- `b = w₃ - N r`. -/
noncomputable def bv (n : ℕ) (σ : ℝ) : ℂ := w3v n σ - (n : ℂ) * rv n σ

/-- `cosh σ`. -/
noncomputable def cv (σ : ℝ) : ℝ := (Real.exp σ + Real.exp (-σ)) / 2

lemma fσ_fact (n : ℕ) (hn : 0 < n) (σ : ℝ) :
    fσ n σ = Complex.exp (-zv (θv n) * (cv σ : ℂ)) * Φv σ * Complex.exp (bv n σ) * H (sv n σ) := by
  have hN : (0 : ℝ) < n := by exact_mod_cast hn
  unfold fσ
  rw [ratio_eq n σ hN, inv_pow_eq n _ (one_add_ε_ne n σ hN), E_sv' n σ hN]
  have hlog : Complex.log (1 + εv n σ) = εv n σ - εv n σ ^ 2 / 2 + εv n σ ^ 3 / 3 + rv n σ := by
    unfold rv; ring
  have e1 : (n : ℂ) * Complex.log (1 + εv n σ) = zv (θv n) / 2 * ((Real.exp (-σ) : ℝ) : ℂ) -
      (-(Complex.I * π / 2) * ((Real.exp (-(2 * σ)) : ℝ) : ℂ)) / 2 +
      ((π / 2 * yv n : ℝ) : ℂ) * ωb ^ 3 * ((Real.exp (-(3 * σ)) : ℝ) : ℂ) / 3 + n * rv n σ := by
    rw [hlog, ← nε1 n hn σ, ← nε2 n hn σ, ← nε3 n hn σ]; ring
  rw [e1]
  have key : Complex.exp (-(zv (θv n) / 2 * ((Real.exp σ : ℝ) : ℂ))) *
      Complex.exp (-(zv (θv n) / 2 * ((Real.exp (-σ) : ℝ) : ℂ) -
        (-(Complex.I * π / 2) * ((Real.exp (-(2 * σ)) : ℝ) : ℂ)) / 2 +
        ((π / 2 * yv n : ℝ) : ℂ) * ωb ^ 3 * ((Real.exp (-(3 * σ)) : ℝ) : ℂ) / 3 + n * rv n σ)) =
      Complex.exp (-zv (θv n) * (cv σ : ℂ)) * Φv σ * Complex.exp (bv n σ) := by
    simp only [Φv, bv, w3v, cv]
    rw [← Complex.exp_add, ← Complex.exp_add, ← Complex.exp_add]
    congr 1
    push_cast
    ring
  rw [mul_assoc, key]
  ring

/-- the explicit model amplitude at real `σ`. -/
noncomputable def Pv (N σ : ℝ) : ℂ := φf σ + (yv N : ℂ) * χf σ

lemma φf_real (σ : ℝ) : φf σ = ((Real.exp (-(3 / 2 * σ)) : ℝ) : ℂ) * Φv σ := by
  unfold φf ψf Φv
  rw [Complex.ofReal_exp, ← Complex.exp_add]
  congr 1
  push_cast; ring_nf

lemma C0P (N σ : ℝ) (hN : 0 < N) :
    C0v N * Pv N σ = (π / sv N σ) ^ (1 / 2 : ℂ) / (√2 * sv N σ) * Φv σ *
      (1 - 3 / (4 * sv N σ) + w3v N σ) := by
  rw [M0_eq N σ hN]
  have hs := sv_ne N σ hN
  have hse := s_ε N σ hN
  have hinv0 : 1 / sv N σ = εv N σ := by
    rw [div_eq_iff hs]; linear_combination -1 * hse
  have hinv : 3 / (4 * sv N σ) = 3 / 4 * εv N σ := by
    rw [← hinv0]; field_simp
  rw [hinv]
  unfold Pv χf w3v εv
  rw [φf_real]
  have e1 : Complex.exp (-(σ : ℂ)) = ((Real.exp (-σ) : ℝ) : ℂ) := by push_cast; rfl
  have e3 : Complex.exp (-3 * (σ : ℂ)) = ((Real.exp (-(3 * σ)) : ℝ) : ℂ) := by push_cast; ring_nf
  rw [e1, e3]
  push_cast
  ring

/-- the window error `err(σ)`. -/
noncomputable def errv (n : ℕ) (σ : ℝ) : ℂ :=
  Φv σ * ((π / sv n σ) ^ (1 / 2 : ℂ) / (√2 * sv n σ) *
      ((1 - 3 / (4 * sv n σ)) * Complex.exp (bv n σ) - 1 + 3 / (4 * sv n σ) - w3v n σ) +
    (H (sv n σ) - (π / sv n σ) ^ (1 / 2 : ℂ) / (√2 * sv n σ) * (1 - 3 / (4 * sv n σ))) *
      Complex.exp (bv n σ))

lemma fσ_split (n : ℕ) (hn : 0 < n) (σ : ℝ) :
    fσ n σ = Complex.exp (-zv (θv n) * (cv σ : ℂ)) * (C0v n * Pv n σ + errv n σ) := by
  have hN : (0 : ℝ) < n := by exact_mod_cast hn
  rw [fσ_fact n hn σ, C0P n σ hN, errv]
  ring

lemma norm_Φv (σ : ℝ) : ‖Φv σ‖ = 1 := by
  unfold Φv
  rw [show -(Complex.I * π / 4) * ((Real.exp (-(2 * σ)) : ℝ) : ℂ) =
    ((-(π / 4 * Real.exp (-(2 * σ))) : ℝ) : ℂ) * Complex.I by push_cast; ring,
    Complex.norm_exp_ofReal_mul_I]

lemma norm_C0v (N : ℝ) (hN : 0 < N) : ‖C0v N‖ = √(π / 2) * yv N * √(yv N) := by
  have hy := yv_pos N hN
  unfold C0v
  rw [norm_mul, Complex.norm_real, Real.norm_eq_abs, abs_of_nonneg (by positivity)]
  rw [Complex.norm_exp]; simp

lemma norm_cexp_sub_one (b : ℂ) : ‖Complex.exp b - 1‖ ≤ ‖b‖ * Real.exp ‖b‖ := by
  have := Complex.norm_exp_sub_sum_le_norm_mul_exp b 1
  simpa using this

lemma norm_cexp_sub_one_sub (b : ℂ) : ‖Complex.exp b - (1 + b)‖ ≤ ‖b‖ ^ 2 * Real.exp ‖b‖ := by
  have := Complex.norm_exp_sub_sum_le_norm_mul_exp b 2
  simpa [Finset.sum_range_succ] using this

/-- `e₁ = y e^{1/4}`, the bound for `|ε|` on the window. -/
noncomputable def e1v (N : ℝ) : ℝ := yv N * Real.exp (1 / 4)

/-- `R₄ = N e₁⁴/(4(1-e₁))`. -/
noncomputable def R4v (n : ℕ) : ℝ := n * (e1v n ^ 4 * (1 - e1v n)⁻¹ / 4)

/-- `B = (π/6) y e^{3/4} + R₄`. -/
noncomputable def Bbv (n : ℕ) : ℝ := π / 6 * yv n * Real.exp (3 / 4) + R4v n

/-- the bound for `|err|` on the window. -/
noncomputable def Errb (n : ℕ) : ℝ :=
  √(π / 2) * yv n * √(yv n) * Real.exp (3 / 8) *
      (Bbv n ^ 2 * Real.exp (Bbv n) + R4v n + 3 / 4 * e1v n * (Bbv n * Real.exp (Bbv n))) +
    45 / 16 * √(2 * π) * (√2 * e1v n) ^ ((2 : ℕ) + 3 / 2 : ℝ) * Real.exp (Bbv n)

lemma errv_bound (n : ℕ) (hn : 0 < n) (he1 : e1v n < 1) (σ : ℝ) (hσ : |σ| ≤ 1 / 4) :
    ‖errv n σ‖ ≤ Errb n := by
  have hN : (0 : ℝ) < n := by exact_mod_cast hn
  have hy := yv_pos n hN
  have hσ' := abs_le.1 hσ
  have hs := sv_ne n σ hN
  have hse := s_ε n σ hN
  have hinv0 : 1 / sv n σ = εv n σ := by
    rw [div_eq_iff hs]; linear_combination -1 * hse
  -- basic norms
  have hεn := norm_εv n σ hN
  have hε : ‖εv n σ‖ ≤ e1v n := by
    rw [hεn]; unfold e1v
    exact mul_le_mul_of_nonneg_left (Real.exp_le_exp.2 (by linarith)) hy.le
  have he10 : 0 ≤ e1v n := by unfold e1v; positivity
  have ha : ‖3 / (4 * sv n σ)‖ ≤ 3 / 4 * e1v n := by
    rw [show (3 : ℂ) / (4 * sv n σ) = 3 / 4 * (1 / sv n σ) by field_simp, hinv0, norm_mul]
    have : ‖(3 / 4 : ℂ)‖ = 3 / 4 := by norm_num
    rw [this]; gcongr
  have hw3 : ‖w3v n σ‖ ≤ π / 6 * yv n * Real.exp (3 / 4) := by
    unfold w3v
    rw [norm_mul, norm_mul, norm_neg, norm_pow, norm_ωb, one_pow, mul_one, Complex.norm_real,
      Complex.norm_real, Real.norm_eq_abs, Real.norm_eq_abs, abs_of_pos (by positivity),
      abs_of_pos (Real.exp_pos _)]
    exact mul_le_mul_of_nonneg_left (Real.exp_le_exp.2 (by linarith)) (by positivity)
  have hr : ‖rv n σ‖ ≤ e1v n ^ 4 * (1 - e1v n)⁻¹ / 4 := by
    have hε1 : ‖εv n σ‖ < 1 := lt_of_le_of_lt hε he1
    refine (log_expand _ hε1).trans ?_
    have h1 : (1 - ‖εv n σ‖)⁻¹ ≤ (1 - e1v n)⁻¹ :=
      inv_anti₀ (by linarith) (by linarith)
    have h2 : ‖εv n σ‖ ^ 4 ≤ e1v n ^ 4 := pow_le_pow_left₀ (norm_nonneg _) hε 4
    have h3 : 0 ≤ (1 - ‖εv n σ‖)⁻¹ := inv_nonneg.2 (by linarith)
    have : ‖εv n σ‖ ^ 4 * (1 - ‖εv n σ‖)⁻¹ ≤ e1v n ^ 4 * (1 - e1v n)⁻¹ :=
      mul_le_mul h2 h1 h3 (by positivity)
    linarith
  have hnr : ‖(n : ℂ) * rv n σ‖ ≤ R4v n := by
    rw [norm_mul, Complex.norm_natCast]; unfold R4v
    exact mul_le_mul_of_nonneg_left hr hN.le
  have hb : ‖bv n σ‖ ≤ Bbv n := by
    unfold bv Bbv
    exact (norm_sub_le _ _).trans (add_le_add hw3 hnr)
  have hB0 : 0 ≤ Bbv n := (norm_nonneg _).trans hb
  have heb : Real.exp ‖bv n σ‖ ≤ Real.exp (Bbv n) := Real.exp_le_exp.2 hb
  -- the M0 part
  have hM0 : ‖(π / sv n σ) ^ (1 / 2 : ℂ) / (√2 * sv n σ)‖ ≤
      √(π / 2) * yv n * √(yv n) * Real.exp (3 / 8) := by
    rw [M0_eq n σ hN, norm_mul, norm_C0v n hN, Complex.norm_real, Real.norm_eq_abs,
      abs_of_pos (Real.exp_pos _)]
    exact mul_le_mul_of_nonneg_left (Real.exp_le_exp.2 (by linarith)) (by positivity)
  have hinner : ‖(1 - 3 / (4 * sv n σ)) * Complex.exp (bv n σ) - 1 + 3 / (4 * sv n σ) - w3v n σ‖ ≤
      Bbv n ^ 2 * Real.exp (Bbv n) + R4v n + 3 / 4 * e1v n * (Bbv n * Real.exp (Bbv n)) := by
    have e : (1 - 3 / (4 * sv n σ)) * Complex.exp (bv n σ) - 1 + 3 / (4 * sv n σ) - w3v n σ =
        (Complex.exp (bv n σ) - (1 + bv n σ)) - (n : ℂ) * rv n σ -
          3 / (4 * sv n σ) * (Complex.exp (bv n σ) - 1) := by
      unfold bv; ring
    rw [e]
    have h1 := norm_cexp_sub_one_sub (bv n σ)
    have h2 := norm_cexp_sub_one (bv n σ)
    have h1' : ‖Complex.exp (bv n σ) - (1 + bv n σ)‖ ≤ Bbv n ^ 2 * Real.exp (Bbv n) :=
      h1.trans (mul_le_mul (pow_le_pow_left₀ (norm_nonneg _) hb 2) heb (by positivity)
        (by positivity))
    have h2' : ‖Complex.exp (bv n σ) - 1‖ ≤ Bbv n * Real.exp (Bbv n) :=
      h2.trans (mul_le_mul hb heb (by positivity) hB0)
    calc _ ≤ ‖Complex.exp (bv n σ) - (1 + bv n σ)‖ + ‖(n : ℂ) * rv n σ‖ +
          ‖3 / (4 * sv n σ)‖ * ‖Complex.exp (bv n σ) - 1‖ := by
          refine (norm_sub_le _ _).trans ?_
          rw [norm_mul]
          exact add_le_add (norm_sub_le _ _) le_rfl
      _ ≤ _ := by
          gcongr
  -- the rH part
  have hrH : ‖H (sv n σ) - (π / sv n σ) ^ (1 / 2 : ℂ) / (√2 * sv n σ) * (1 - 3 / (4 * sv n σ))‖ ≤
      45 / 16 * √(2 * π) * (√2 * e1v n) ^ ((2 : ℕ) + 3 / 2 : ℝ) := by
    refine (H_approx _ (sv_re_pos n σ hN)).trans ?_
    apply mul_le_mul_of_nonneg_left _ (by positivity)
    apply Real.rpow_le_rpow (by have := sv_re_pos n σ hN; positivity) _ (by positivity)
    rw [sv_re]
    have hρ := ρv_pos n hN
    have hρy := ρy n hN
    have he : Real.exp σ * Real.exp (-σ) = 1 := by rw [← Real.exp_add]; simp
    have hyρ : yv n = 1 / ρv n := by field_simp; linarith
    have : 1 / (ρv n * Real.exp σ * (√2 / 2)) = √2 * yv n * Real.exp (-σ) := by
      rw [hyρ]; field_simp; rw [Real.sq_sqrt (by norm_num)]; nlinarith
    rw [this]; unfold e1v
    have : yv n * Real.exp (-σ) ≤ yv n * Real.exp (1 / 4) :=
      mul_le_mul_of_nonneg_left (Real.exp_le_exp.2 (by linarith)) hy.le
    have h2 := sqrt2_pos
    nlinarith
  unfold errv Errb
  rw [norm_mul, norm_Φv, one_mul]
  refine (norm_add_le _ _).trans (add_le_add ?_ ?_)
  · rw [norm_mul]; exact mul_le_mul hM0 hinner (norm_nonneg _) (by positivity)
  · rw [norm_mul]
    refine mul_le_mul hrH ((Complex.norm_exp_le_exp_norm _).trans heb) (norm_nonneg _) (by positivity)

/-! ## 13. The window integral -/

lemma ii_le_real (g : ℝ → ℝ) (hg : Integrable g) (h0 : ∀ x, 0 ≤ g x) (δ : ℝ) (hδ : 0 ≤ δ) :
    ∫ x in -δ..δ, g x ≤ ∫ x, g x := by
  rw [intervalIntegral.integral_of_le (by linarith)]
  exact setIntegral_le_integral hg (Filter.Eventually.of_forall h0)

lemma norm_ii_le (f : ℝ → ℂ) (g : ℝ → ℝ) (hg : Integrable g) (h0 : ∀ x, 0 ≤ g x) (δ : ℝ)
    (hδ : 0 ≤ δ) (h : ∀ x, |x| ≤ δ → ‖f x‖ ≤ g x) : ‖∫ x in -δ..δ, f x‖ ≤ ∫ x, g x := by
  refine (intervalIntegral.norm_integral_le_of_norm_le (by linarith)
    (Filter.Eventually.of_forall fun x hx => h x ?_) hg.intervalIntegrable).trans
    (ii_le_real g hg h0 δ hδ)
  rw [abs_le]; exact ⟨by linarith [hx.1], hx.2⟩

lemma J_split3 (g : ℝ → ℂ) (hg : Integrable g) (δ : ℝ) (hδ : 0 ≤ δ) :
    ∫ σ, g σ = (∫ σ in Iic (-δ), g σ) + (∫ σ in -δ..δ, g σ) + ∫ σ in Ioi δ, g σ := by
  rw [split_Iic_Ioi g hg (-δ), add_assoc]
  congr 1
  rw [← Ioc_union_Ioi_eq_Ioi (by linarith : -δ ≤ δ),
    setIntegral_union (Ioc_disjoint_Ioi le_rfl) measurableSet_Ioi hg.integrableOn hg.integrableOn,
    intervalIntegral.integral_of_le (by linarith)]

lemma cexp_zc (θ σ : ℝ) : Complex.exp (-zv θ * (cv σ : ℂ)) = Complex.exp (-zv θ) * Kf θ σ := by
  unfold Kf cv; rw [← Complex.exp_add]; congr 1; push_cast; ring

lemma Kf_cont (θ : ℝ) : Continuous (Kf θ) := by unfold Kf; fun_prop

lemma Pv_cont (N : ℝ) : Continuous (Pv N) := by
  unfold Pv
  have h1 := φf_diff.continuous; have h2 := χf_diff.continuous
  fun_prop

lemma window_split (n : ℕ) (hf : Integrable (fσ n)) (δ : ℝ) :
    ∫ σ in -δ..δ, fσ n σ = C0v n * Complex.exp (-zv (θv n)) * (∫ σ in -δ..δ, Kf (θv n) σ * Pv n σ) +
      ∫ σ in -δ..δ, (fσ n σ - C0v n * Complex.exp (-zv (θv n)) * (Kf (θv n) σ * Pv n σ)) := by
  have hc : IntervalIntegrable (fun σ => C0v n * Complex.exp (-zv (θv n)) * (Kf (θv n) σ * Pv n σ))
      volume (-δ) δ :=
    ((Kf_cont _).mul (Pv_cont _)).intervalIntegrable _ _ |>.const_mul _
  have e1 : C0v n * Complex.exp (-zv (θv n)) * (∫ σ in -δ..δ, Kf (θv n) σ * Pv n σ) =
      ∫ σ in -δ..δ, C0v n * Complex.exp (-zv (θv n)) * (Kf (θv n) σ * Pv n σ) :=
    (intervalIntegral.integral_const_mul _ _).symm
  rw [e1, ← intervalIntegral.integral_add hc ((hf.intervalIntegrable).sub hc)]
  congr 1; funext σ; ring

lemma D_bound (n : ℕ) (hn : 0 < n) (he1 : e1v n < 1) (σ : ℝ) (hσ : |σ| ≤ 1 / 4) :
    ‖fσ n σ - C0v n * Complex.exp (-zv (θv n)) * (Kf (θv n) σ * Pv n σ)‖ ≤
      Errb n * Real.exp (-θv n) * Real.exp (-(θv n / 2) * σ ^ 2) := by
  have hN : (0 : ℝ) < n := by exact_mod_cast hn
  have hθ := θv_pos n hN
  have e : fσ n σ - C0v n * Complex.exp (-zv (θv n)) * (Kf (θv n) σ * Pv n σ) =
      Complex.exp (-zv (θv n) * (cv σ : ℂ)) * errv n σ := by
    rw [fσ_split n hn σ, cexp_zc]; ring
  rw [e, norm_mul, norm_cexp_mz]
  have h1 := errv_bound n hn he1 σ hσ
  have h2 : Real.exp (-(θv n * cv σ)) ≤ Real.exp (-θv n) * Real.exp (-(θv n / 2) * σ ^ 2) := by
    rw [← Real.exp_add]; apply Real.exp_le_exp.2
    have := cosh_quad σ; unfold cv; nlinarith
  calc Real.exp (-(θv n * cv σ)) * ‖errv n σ‖ ≤
      (Real.exp (-θv n) * Real.exp (-(θv n / 2) * σ ^ 2)) * Errb n :=
        mul_le_mul h2 h1 (norm_nonneg _) (by positivity)
    _ = _ := by ring

lemma W1 (n : ℕ) (hn : 0 < n) (he1 : e1v n < 1) :
    ‖∫ σ in -(1 / 4 : ℝ)..(1 / 4), (fσ n σ - C0v n * Complex.exp (-zv (θv n)) *
      (Kf (θv n) σ * Pv n σ))‖ ≤ Errb n * Real.exp (-θv n) * √(π / (θv n / 2)) := by
  have hN : (0 : ℝ) < n := by exact_mod_cast hn
  have hθ := θv_pos n hN
  have hEb : 0 ≤ Errb n := (norm_nonneg _).trans (errv_bound n hn he1 0 (by norm_num))
  have hg : Integrable (fun σ : ℝ => Errb n * Real.exp (-θv n) * Real.exp (-(θv n / 2) * σ ^ 2)) :=
    (integrable_exp_neg_mul_sq (by positivity)).const_mul _
  refine (norm_ii_le _ _ hg (fun x => by positivity) (1 / 4) (by norm_num)
    (fun x hx => D_bound n hn he1 x hx)).trans (le_of_eq ?_)
  rw [integral_const_mul, integral_gaussian]

lemma odd_zero (g : ℝ → ℂ) (hodd : ∀ x, g (-x) = -g x) (δ : ℝ) : ∫ x in -δ..δ, g x = 0 := by
  have h := intervalIntegral.integral_comp_neg (a := -δ) (b := δ) g
  rw [neg_neg] at h
  simp only [hodd, intervalIntegral.integral_neg] at h
  have : (2 : ℂ) * ∫ x in -δ..δ, g x = 0 := by linear_combination -h
  simpa using this

lemma Kf_even (θ σ : ℝ) : Kf θ (-σ) = Kf θ σ := by
  unfold Kf; congr 3; rw [neg_neg]; ring

/-- the coefficients `p_j = φ_j + y χ_j`. -/
noncomputable def pc (N : ℝ) (j : ℕ) : ℂ := tc φf j + (yv N : ℂ) * tc χf j

lemma norm_sigma (σ : ℝ) (hσ : |σ| ≤ 1 / 4) : ‖(σ : ℂ)‖ ≤ 1 / 4 := by
  rw [Complex.norm_real, Real.norm_eq_abs]; exact hσ

lemma W23 (n : ℕ) (hn : 0 < n) :
    ‖(∫ σ in -(1 / 4 : ℝ)..(1 / 4), Kf (θv n) σ * Pv n σ) -
      ∫ σ in -(1 / 4 : ℝ)..(1 / 4), Kf (θv n) σ * (pc n 0 + pc n 2 * (σ : ℂ) ^ 2)‖ ≤
      32 * (Mφ + yv n * Mχ) * ∫ x : ℝ, x ^ 4 * Real.exp (-(θv n / 2) * x ^ 2) := by
  have hN : (0 : ℝ) < n := by exact_mod_cast hn
  have hθ := θv_pos n hN
  have hy := yv_pos n hN
  set θ := θv n
  set y := yv n
  -- the odd part
  set g1 : ℝ → ℂ := fun σ => Kf θ σ * (pc n 1 * σ + pc n 3 * (σ : ℂ) ^ 3)
  have hodd : ∫ σ in -(1 / 4 : ℝ)..(1 / 4), g1 σ = 0 := by
    refine odd_zero g1 (fun x => ?_) _
    simp only [g1, Kf_even]; push_cast; ring
  -- the remainder
  set g2 : ℝ → ℂ := fun σ => Kf θ σ * ((φf σ - (tc φf 0 + tc φf 1 * σ + tc φf 2 * (σ : ℂ) ^ 2 +
    tc φf 3 * (σ : ℂ) ^ 3)) + (y : ℂ) * (χf σ - (tc χf 0 + tc χf 1 * σ + tc χf 2 * (σ : ℂ) ^ 2 +
    tc χf 3 * (σ : ℂ) ^ 3)))
  have hc : ∀ g : ℝ → ℂ, Continuous g → IntervalIntegrable g volume (-(1 / 4 : ℝ)) (1 / 4) :=
    fun g hg => hg.intervalIntegrable _ _
  have hK := Kf_cont θ
  have hφ := φf_diff.continuous
  have hχ := χf_diff.continuous
  have e : (∫ σ in -(1 / 4 : ℝ)..(1 / 4), Kf θ σ * Pv n σ) -
      ∫ σ in -(1 / 4 : ℝ)..(1 / 4), Kf θ σ * (pc n 0 + pc n 2 * (σ : ℂ) ^ 2) =
      (∫ σ in -(1 / 4 : ℝ)..(1 / 4), g1 σ) + ∫ σ in -(1 / 4 : ℝ)..(1 / 4), g2 σ := by
    rw [← intervalIntegral.integral_sub (hc _ (by unfold Pv; fun_prop)) (hc _ (by fun_prop)),
      ← intervalIntegral.integral_add (hc _ (by simp only [g1]; fun_prop))
        (hc _ (by simp only [g2]; fun_prop))]
    congr 1; funext σ
    simp only [g1, g2, Pv, pc]; ring
  rw [e, hodd, zero_add]
  have hint : Integrable (fun x : ℝ => 32 * (Mφ + y * Mχ) * (x ^ 4 * Real.exp (-(θ / 2) * x ^ 2))) :=
    (integrable_pow_gauss 4 (θ / 2) (by positivity)).const_mul _
  have hMφ : 0 ≤ Mφ := by unfold Mφ; positivity
  have hMχ : 0 ≤ Mχ := by unfold Mχ; positivity
  rw [← integral_const_mul]
  refine norm_ii_le _ _ hint (fun x => by positivity) (1 / 4) (by norm_num) (fun σ hσ => ?_)
  simp only [g2]
  rw [norm_mul]
  have h1 := taylor4 φf φf_diff Mφ φf_sphere σ (norm_sigma σ hσ)
  have h2 := taylor4 χf χf_diff Mχ χf_sphere σ (norm_sigma σ hσ)
  have hs4 : ‖(σ : ℂ)‖ ^ 4 = σ ^ 4 := by
    rw [Complex.norm_real, Real.norm_eq_abs, show (4 : ℕ) = 2 * 2 by rfl, pow_mul, sq_abs, ← pow_mul]
  rw [hs4] at h1 h2
  have hK' := norm_Kf_le θ σ hθ.le
  have hsum : ‖(φf σ - (tc φf 0 + tc φf 1 * σ + tc φf 2 * (σ : ℂ) ^ 2 + tc φf 3 * (σ : ℂ) ^ 3)) +
      (y : ℂ) * (χf σ - (tc χf 0 + tc χf 1 * σ + tc χf 2 * (σ : ℂ) ^ 2 + tc χf 3 * (σ : ℂ) ^ 3))‖ ≤
      32 * (Mφ + y * Mχ) * σ ^ 4 := by
    refine (norm_add_le _ _).trans ?_
    rw [norm_mul, Complex.norm_real, Real.norm_eq_abs, abs_of_pos hy]
    have := mul_le_mul_of_nonneg_left h2 hy.le
    nlinarith
  have hE : Real.exp (-(θ * (σ ^ 2 / 2))) = Real.exp (-(θ / 2) * σ ^ 2) := by ring_nf
  rw [hE] at hK'
  calc ‖Kf θ σ‖ * _ ≤ Real.exp (-(θ / 2) * σ ^ 2) * (32 * (Mφ + y * Mχ) * σ ^ 4) :=
        mul_le_mul hK' hsum (norm_nonneg _) (by positivity)
    _ = _ := by ring

/-- `Q₄(σ) = p₀ + p₂σ² - (z/24) p₀ σ⁴`. -/
noncomputable def Q4 (N : ℝ) (σ : ℝ) : ℂ :=
  pc N 0 + pc N 2 * (σ : ℂ) ^ 2 - zv (θv N) / 24 * pc N 0 * (σ : ℂ) ^ 4

lemma gm_int (j : ℕ) (c : ℝ) (hc : 0 < c) (a : ℝ) :
    Integrable (fun x : ℝ => a * (x ^ j * Real.exp (-c * x ^ 2))) :=
  (integrable_pow_gauss j c hc).const_mul a

lemma W4 (n : ℕ) (hn : 0 < n) :
    ‖∫ σ in -(1 / 4 : ℝ)..(1 / 4), (Kf (θv n) σ * (pc n 0 + pc n 2 * (σ : ℂ) ^ 2) -
      Gf (θv n) σ * Q4 n σ)‖ ≤
    0.0023 * θv n * ‖pc n 0‖ * (∫ x : ℝ, x ^ 6 * Real.exp (-(0.4963 * θv n) * x ^ 2)) +
    (0.0023 * θv n * ‖pc n 2‖ + 0.0035 * θv n ^ 2 * ‖pc n 0‖) *
      (∫ x : ℝ, x ^ 8 * Real.exp (-(0.4963 * θv n) * x ^ 2)) +
    0.0035 * θv n ^ 2 * ‖pc n 2‖ * (∫ x : ℝ, x ^ 10 * Real.exp (-(0.4963 * θv n) * x ^ 2)) +
    √2 * θv n / 24 * ‖pc n 2‖ * (∫ x : ℝ, x ^ 6 * Real.exp (-(θv n / 2) * x ^ 2)) := by
  have hN : (0 : ℝ) < n := by exact_mod_cast hn
  have hθ := θv_pos n hN
  set θ := θv n
  set A := ‖pc n 0‖
  set B := ‖pc n 2‖
  have hA : 0 ≤ A := norm_nonneg _
  have hB : 0 ≤ B := norm_nonneg _
  set c1 := 0.4963 * θ
  have hc1 : 0 < c1 := by simp only [c1]; positivity
  have hc2 : 0 < θ / 2 := by positivity
  set g : ℝ → ℝ := fun x => 0.0023 * θ * A * (x ^ 6 * Real.exp (-c1 * x ^ 2)) +
    (0.0023 * θ * B + 0.0035 * θ ^ 2 * A) * (x ^ 8 * Real.exp (-c1 * x ^ 2)) +
    0.0035 * θ ^ 2 * B * (x ^ 10 * Real.exp (-c1 * x ^ 2)) +
    √2 * θ / 24 * B * (x ^ 6 * Real.exp (-(θ / 2) * x ^ 2))
  have hgi : Integrable g :=
    (((gm_int 6 c1 hc1 _).add (gm_int 8 c1 hc1 _)).add (gm_int 10 c1 hc1 _)).add (gm_int 6 _ hc2 _)
  have hg0 : ∀ x, 0 ≤ g x := by
    intro x; simp only [g]
    have h6 : 0 ≤ x ^ 6 := by positivity
    have h8 : 0 ≤ x ^ 8 := by positivity
    have h10 : 0 ≤ x ^ 10 := by positivity
    positivity
  have hval : ∫ x, g x = 0.0023 * θ * A * (∫ x : ℝ, x ^ 6 * Real.exp (-c1 * x ^ 2)) +
      (0.0023 * θ * B + 0.0035 * θ ^ 2 * A) * (∫ x : ℝ, x ^ 8 * Real.exp (-c1 * x ^ 2)) +
      0.0035 * θ ^ 2 * B * (∫ x : ℝ, x ^ 10 * Real.exp (-c1 * x ^ 2)) +
      √2 * θ / 24 * B * (∫ x : ℝ, x ^ 6 * Real.exp (-(θ / 2) * x ^ 2)) := by
    simp only [g]
    have i1 := gm_int 6 c1 hc1 (0.0023 * θ * A)
    have i2 := gm_int 8 c1 hc1 (0.0023 * θ * B + 0.0035 * θ ^ 2 * A)
    have i3 := gm_int 10 c1 hc1 (0.0035 * θ ^ 2 * B)
    have i4 := gm_int 6 (θ / 2) hc2 (√2 * θ / 24 * B)
    have i12 : Integrable (fun x : ℝ => 0.0023 * θ * A * (x ^ 6 * Real.exp (-c1 * x ^ 2)) +
        (0.0023 * θ * B + 0.0035 * θ ^ 2 * A) * (x ^ 8 * Real.exp (-c1 * x ^ 2))) := i1.add i2
    have i123 : Integrable (fun x : ℝ => 0.0023 * θ * A * (x ^ 6 * Real.exp (-c1 * x ^ 2)) +
        (0.0023 * θ * B + 0.0035 * θ ^ 2 * A) * (x ^ 8 * Real.exp (-c1 * x ^ 2)) +
        0.0035 * θ ^ 2 * B * (x ^ 10 * Real.exp (-c1 * x ^ 2))) := i12.add i3
    rw [integral_add i123 i4, integral_add i12 i3, integral_add i1 i2]
    simp only [integral_const_mul]
  rw [← hval]
  refine norm_ii_le _ g hgi hg0 (1 / 4) (by norm_num) (fun σ hσ => ?_)
  have e : Kf θ σ * (pc n 0 + pc n 2 * (σ : ℂ) ^ 2) - Gf θ σ * Q4 n σ =
      (Kf θ σ - Gf θ σ * (1 - zv θ * (σ : ℂ) ^ 4 / 24)) * (pc n 0 + pc n 2 * (σ : ℂ) ^ 2) -
        Gf θ σ * (zv θ / 24 * pc n 2 * (σ : ℂ) ^ 6) := by
    simp only [Q4]; ring
  rw [e]
  have h1 := er_bound θ σ hθ.le hσ
  have hp : ‖pc n 0 + pc n 2 * (σ : ℂ) ^ 2‖ ≤ A + B * σ ^ 2 := by
    refine (norm_add_le _ _).trans ?_
    rw [norm_mul, norm_pow, Complex.norm_real, Real.norm_eq_abs, sq_abs]
  have h3 : ‖Gf θ σ * (zv θ / 24 * pc n 2 * (σ : ℂ) ^ 6)‖ ≤
      Real.exp (-(θ / 2) * σ ^ 2) * (√2 * θ / 24 * B * σ ^ 6) := by
    rw [norm_mul, norm_Gf, norm_mul, norm_mul, norm_div, norm_zv θ hθ.le, norm_pow,
      Complex.norm_real, Real.norm_eq_abs]
    have : |σ| ^ 6 = σ ^ 6 := by
      rw [show (6 : ℕ) = 2 * 3 by rfl, pow_mul, sq_abs, ← pow_mul]
    rw [this]
    have : ‖(24 : ℂ)‖ = 24 := by norm_num
    rw [this]
    apply le_of_eq; simp only [B]; ring_nf
  refine (norm_sub_le _ _).trans ?_
  rw [norm_mul]
  have h4 : ‖Kf θ σ - Gf θ σ * (1 - zv θ * (σ : ℂ) ^ 4 / 24)‖ * ‖pc n 0 + pc n 2 * (σ : ℂ) ^ 2‖ ≤
      Real.exp (-(0.4963 * θ * σ ^ 2)) * (0.0023 * θ * σ ^ 6 + 0.0035 * θ ^ 2 * σ ^ 8) *
        (A + B * σ ^ 2) := mul_le_mul h1 hp (norm_nonneg _) (by positivity)
  have e2 : Real.exp (-(0.4963 * θ * σ ^ 2)) = Real.exp (-c1 * σ ^ 2) := by
    simp only [c1]; ring_nf
  rw [e2] at h4
  have : g σ = Real.exp (-c1 * σ ^ 2) * (0.0023 * θ * σ ^ 6 + 0.0035 * θ ^ 2 * σ ^ 8) *
      (A + B * σ ^ 2) + Real.exp (-(θ / 2) * σ ^ 2) * (√2 * θ / 24 * B * σ ^ 6) := by
    simp only [g]; ring
  rw [this]
  exact add_le_add h4 h3

lemma Gf_eq (θ σ : ℝ) : Gf θ σ = Complex.exp (-(zv θ / 2) * (σ : ℂ) ^ 2) := by
  unfold Gf; congr 1; push_cast; ring

lemma zv_half_re (θ : ℝ) : (zv θ / 2).re = θ / 2 := by simp [zv]

lemma GQ_cont (n : ℕ) : Continuous (fun σ : ℝ => Gf (θv n) σ * Q4 n σ) := by
  unfold Gf Q4; fun_prop

lemma GQ_int (n : ℕ) (hn : 0 < n) : Integrable (fun σ : ℝ => Gf (θv n) σ * Q4 n σ) := by
  have hN : (0 : ℝ) < n := by exact_mod_cast hn
  have hθ := θv_pos n hN
  have hb : 0 < (zv (θv n) / 2).re := by rw [zv_half_re]; positivity
  have h0 := integrable_cpow_gauss 0 _ hb
  have h2 := integrable_cpow_gauss 2 _ hb
  have h4 := integrable_cpow_gauss 4 _ hb
  have : (fun σ : ℝ => Gf (θv n) σ * Q4 n σ) = fun σ : ℝ =>
      pc n 0 * ((σ : ℂ) ^ 0 * Complex.exp (-(zv (θv n) / 2) * (σ : ℂ) ^ 2)) +
      pc n 2 * ((σ : ℂ) ^ 2 * Complex.exp (-(zv (θv n) / 2) * (σ : ℂ) ^ 2)) -
      zv (θv n) / 24 * pc n 0 * ((σ : ℂ) ^ 4 * Complex.exp (-(zv (θv n) / 2) * (σ : ℂ) ^ 2)) := by
    funext σ; rw [Gf_eq]; unfold Q4; ring
  rw [this]
  exact ((h0.const_mul _).add (h2.const_mul _)).sub (h4.const_mul _)

/-- `m₀ = (π/(z/2))^{1/2} = ∫ e^{-zσ²/2}`. -/
noncomputable def m0c (θ : ℝ) : ℂ := (π / (zv θ / 2)) ^ (1 / 2 : ℂ)

lemma W6 (n : ℕ) (hn : 0 < n) :
    ∫ σ : ℝ, Gf (θv n) σ * Q4 n σ =
      m0c (θv n) * (pc n 0 + pc n 2 / zv (θv n) - pc n 0 / (8 * zv (θv n))) := by
  have hN : (0 : ℝ) < n := by exact_mod_cast hn
  have hθ := θv_pos n hN
  have hb : 0 < (zv (θv n) / 2).re := by rw [zv_half_re]; positivity
  have hz : zv (θv n) ≠ 0 := by
    intro h; have := congrArg Complex.re h; rw [zv_re] at this; simp at this; linarith
  have h0 := integrable_cpow_gauss 0 _ hb
  have h2 := integrable_cpow_gauss 2 _ hb
  have h4 := integrable_cpow_gauss 4 _ hb
  have : (fun σ : ℝ => Gf (θv n) σ * Q4 n σ) = fun σ : ℝ =>
      pc n 0 * ((σ : ℂ) ^ 0 * Complex.exp (-(zv (θv n) / 2) * (σ : ℂ) ^ 2)) +
      pc n 2 * ((σ : ℂ) ^ 2 * Complex.exp (-(zv (θv n) / 2) * (σ : ℂ) ^ 2)) -
      zv (θv n) / 24 * pc n 0 * ((σ : ℂ) ^ 4 * Complex.exp (-(zv (θv n) / 2) * (σ : ℂ) ^ 2)) := by
    funext σ; rw [Gf_eq]; unfold Q4; ring
  have j0 := h0.const_mul (pc n 0)
  have j2 := h2.const_mul (pc n 2)
  have j4 := h4.const_mul (zv (θv n) / 24 * pc n 0)
  have j02 : Integrable (fun σ : ℝ =>
      pc n 0 * ((σ : ℂ) ^ 0 * Complex.exp (-(zv (θv n) / 2) * (σ : ℂ) ^ 2)) +
      pc n 2 * ((σ : ℂ) ^ 2 * Complex.exp (-(zv (θv n) / 2) * (σ : ℂ) ^ 2))) := j0.add j2
  rw [this, integral_sub j02 j4, integral_add j0 j2, integral_const_mul, integral_const_mul,
    integral_const_mul, gauss_mom0 _ hb, gauss_mom2 _ hb, gauss_mom4 _ hb]
  unfold m0c
  field_simp
  ring

lemma W5 (n : ℕ) (hn : 0 < n) :
    ‖(∫ σ in -(1 / 4 : ℝ)..(1 / 4), Gf (θv n) σ * Q4 n σ) - ∫ σ : ℝ, Gf (θv n) σ * Q4 n σ‖ ≤
      Real.exp (-(θv n / 64)) * (‖pc n 0‖ * (∫ x : ℝ, x ^ 0 * Real.exp (-(θv n / 4) * x ^ 2)) +
        ‖pc n 2‖ * (∫ x : ℝ, x ^ 2 * Real.exp (-(θv n / 4) * x ^ 2)) +
        √2 * θv n / 24 * ‖pc n 0‖ * (∫ x : ℝ, x ^ 4 * Real.exp (-(θv n / 4) * x ^ 2))) := by
  have hN : (0 : ℝ) < n := by exact_mod_cast hn
  have hθ := θv_pos n hN
  set θ := θv n
  set A := ‖pc n 0‖
  set B := ‖pc n 2‖
  have hA : 0 ≤ A := norm_nonneg _
  have hB : 0 ≤ B := norm_nonneg _
  have hc : 0 < θ / 4 := by positivity
  set b : ℝ → ℝ := fun x => Real.exp (-(θ / 64)) * (A * (x ^ 0 * Real.exp (-(θ / 4) * x ^ 2)) +
    B * (x ^ 2 * Real.exp (-(θ / 4) * x ^ 2)) + √2 * θ / 24 * A * (x ^ 4 * Real.exp (-(θ / 4) * x ^ 2)))
  have i0 := gm_int 0 _ hc A
  have i2 := gm_int 2 _ hc B
  have i4 := gm_int 4 _ hc (√2 * θ / 24 * A)
  have i02 : Integrable (fun x : ℝ => A * (x ^ 0 * Real.exp (-(θ / 4) * x ^ 2)) +
      B * (x ^ 2 * Real.exp (-(θ / 4) * x ^ 2))) := i0.add i2
  have i024 : Integrable (fun x : ℝ => A * (x ^ 0 * Real.exp (-(θ / 4) * x ^ 2)) +
      B * (x ^ 2 * Real.exp (-(θ / 4) * x ^ 2)) +
      √2 * θ / 24 * A * (x ^ 4 * Real.exp (-(θ / 4) * x ^ 2))) := i02.add i4
  have hbi : Integrable b := i024.const_mul _
  have hb0 : ∀ x, 0 ≤ b x := by
    intro x; simp only [b]
    have : 0 ≤ x ^ 2 := by positivity
    have : 0 ≤ x ^ 4 := by positivity
    positivity
  have hbv : ∫ x, b x = Real.exp (-(θ / 64)) * (A * (∫ x : ℝ, x ^ 0 * Real.exp (-(θ / 4) * x ^ 2)) +
      B * (∫ x : ℝ, x ^ 2 * Real.exp (-(θ / 4) * x ^ 2)) +
      √2 * θ / 24 * A * (∫ x : ℝ, x ^ 4 * Real.exp (-(θ / 4) * x ^ 2))) := by
    simp only [b]
    rw [integral_const_mul, integral_add i02 i4, integral_add i0 i2]
    simp only [integral_const_mul]
  rw [← hbv]
  -- pointwise outside the window
  have hpt : ∀ x : ℝ, 1 / 4 ≤ |x| → ‖Gf θ x * Q4 n x‖ ≤ b x := by
    intro x hx
    have hx2 : 1 / 16 ≤ x ^ 2 := by
      have := sq_abs x; nlinarith [abs_nonneg x]
    rw [norm_mul, norm_Gf]
    have hQ : ‖Q4 n x‖ ≤ A + B * x ^ 2 + √2 * θ / 24 * A * x ^ 4 := by
      unfold Q4
      refine (norm_sub_le _ _).trans (add_le_add ((norm_add_le _ _).trans (le_of_eq ?_)) (le_of_eq ?_))
      · rw [norm_mul, norm_pow, Complex.norm_real, Real.norm_eq_abs, sq_abs]
      · rw [norm_mul, norm_mul, norm_div, norm_zv θ hθ.le, norm_pow, Complex.norm_real,
          Real.norm_eq_abs]
        have : |x| ^ 4 = x ^ 4 := by
          rw [show (4 : ℕ) = 2 * 2 by rfl, pow_mul, sq_abs, ← pow_mul]
        rw [this]
        have : ‖(24 : ℂ)‖ = 24 := by norm_num
        rw [this]
    have hE : Real.exp (-(θ * (x ^ 2 / 2))) ≤ Real.exp (-(θ / 64)) * Real.exp (-(θ / 4) * x ^ 2) := by
      rw [← Real.exp_add]; apply Real.exp_le_exp.2; nlinarith
    have : b x = Real.exp (-(θ / 64)) * Real.exp (-(θ / 4) * x ^ 2) *
        (A + B * x ^ 2 + √2 * θ / 24 * A * x ^ 4) := by simp only [b]; ring
    rw [this]
    exact mul_le_mul hE hQ (norm_nonneg _) (by positivity)
  have hgi := GQ_int n hn
  rw [J_split3 _ hgi (1 / 4) (by norm_num)]
  have e : (∫ σ in -(1 / 4 : ℝ)..(1 / 4), Gf θ σ * Q4 n σ) -
      ((∫ σ in Iic (-(1 / 4 : ℝ)), Gf θ σ * Q4 n σ) + (∫ σ in -(1 / 4 : ℝ)..(1 / 4), Gf θ σ * Q4 n σ) +
        ∫ σ in Ioi (1 / 4 : ℝ), Gf θ σ * Q4 n σ) =
      -((∫ σ in Iic (-(1 / 4 : ℝ)), Gf θ σ * Q4 n σ) + ∫ σ in Ioi (1 / 4 : ℝ), Gf θ σ * Q4 n σ) := by
    ring
  rw [e, norm_neg]
  have hL : ‖∫ σ in Iic (-(1 / 4 : ℝ)), Gf θ σ * Q4 n σ‖ ≤ ∫ σ in Iic (-(1 / 4 : ℝ)), b σ := by
    refine norm_integral_le_of_norm_le hbi.integrableOn
      ((ae_restrict_iff' measurableSet_Iic).2 (Filter.Eventually.of_forall fun x hx => hpt x ?_))
    have : x ≤ -(1 / 4) := hx
    rw [abs_of_neg (by linarith)]; linarith
  have hR : ‖∫ σ in Ioi (1 / 4 : ℝ), Gf θ σ * Q4 n σ‖ ≤ ∫ σ in Ioi (1 / 4 : ℝ), b σ := by
    refine norm_integral_le_of_norm_le hbi.integrableOn
      ((ae_restrict_iff' measurableSet_Ioi).2 (Filter.Eventually.of_forall fun x hx => hpt x ?_))
    have : 1 / 4 < x := hx
    rw [abs_of_pos (by linarith)]; linarith
  have hsplit : ∫ x, b x = (∫ x in Iic (-(1 / 4 : ℝ)), b x) + ∫ x in Ioi (-(1 / 4 : ℝ)), b x := by
    rw [← setIntegral_univ, ← Iic_union_Ioi (a := -(1 / 4 : ℝ))]
    exact setIntegral_union (Iic_disjoint_Ioi le_rfl) measurableSet_Ioi hbi.integrableOn
      hbi.integrableOn
  have hmono : ∫ x in Ioi (1 / 4 : ℝ), b x ≤ ∫ x in Ioi (-(1 / 4 : ℝ)), b x :=
    setIntegral_mono_set hbi.integrableOn (Filter.Eventually.of_forall hb0)
      (Filter.Eventually.of_forall (Ioi_subset_Ioi (by norm_num)))
  rw [hsplit]
  exact (norm_add_le _ _).trans (add_le_add hL (hR.trans hmono))

/-! ## 14. The model equals `T` up to `O(y²)` -/

/-- `η₁ = -(3/4) ω̄ - (π/6) ω̄³`. -/
noncomputable def η1 : ℂ := -(3 / 4 * ωb + π / 6 * ωb ^ 3)

lemma ωb_cube : ωb ^ 3 = ((√2 / 2 : ℝ) : ℂ) * (-1 - Complex.I) := by
  rw [pow_succ, ωb_sq, ωb_eq]; ring_nf; rw [Complex.I_sq]; ring

lemma η1_eq : η1 = ((√2 / 2 : ℝ) : ℂ) * ((-(3 / 4) + π / 6 : ℝ) + (3 / 4 + π / 6 : ℝ) * Complex.I) := by
  unfold η1; rw [ωb_cube, ωb_eq]; push_cast; ring

lemma θ_eq (N : ℝ) : θv N = √π * √N := by
  unfold θv; rw [Real.sqrt_mul Real.pi_pos.le]

lemma y_sqrt2 (N : ℝ) (hN : 0 < N) : yv N * √2 = √π / √N := by
  unfold yv
  rw [← Real.sqrt_mul (by positivity), show π / (2 * N) * 2 = π / N by field_simp,
    Real.sqrt_div' π hN.le]

lemma inv_zv (θ : ℝ) (hθ : 0 < θ) : 1 / zv θ = (1 + Complex.I) / (2 * θ) := by
  have h1 : (1 - Complex.I) * (1 + Complex.I) = 2 := by ring_nf; rw [Complex.I_sq]; ring
  have hθ' : (θ : ℂ) ≠ 0 := by exact_mod_cast hθ.ne'
  have h2 : (1 - Complex.I) ≠ 0 := by
    intro h; have := congrArg Complex.im h; simp at this
  unfold zv
  field_simp
  linear_combination -h1

lemma ν0_eq : ν0 = (((9 - π ^ 2) / 8 : ℝ) : ℂ) + ((-(5 * π / 4)) : ℝ) * Complex.I := by
  unfold ν0; push_cast; ring_nf; rw [Complex.I_sq]; ring

lemma K2_alg (a b : ℝ) (ha : 0 < a) (hb : 0 < b) :
    ((a / b / 2 : ℝ) : ℂ) * ((-(3 / 4) + a ^ 2 / 6 : ℝ) + (3 / 4 + a ^ 2 / 6 : ℝ) * Complex.I) +
      ((((9 - (a ^ 2) ^ 2) / 8 : ℝ) : ℂ) + ((-(5 * a ^ 2 / 4)) : ℝ) * Complex.I - 1 / 8) *
        ((1 + Complex.I) / (2 * ((a * b : ℝ) : ℂ))) =
      (⟨(24 + 12 * a ^ 2 + (a ^ 2) ^ 2) / (48 * a), (24 - 12 * a ^ 2 + (a ^ 2) ^ 2) / (48 * a)⟩ : ℂ) /
        (b : ℂ) := by
  have ha' : (a : ℂ) ≠ 0 := by exact_mod_cast ha.ne'
  have hb' : (b : ℂ) ≠ 0 := by exact_mod_cast hb.ne'
  apply Complex.ext
  · simp [Complex.div_re, Complex.div_im, Complex.mul_re, Complex.mul_im]
    simp only [← Complex.ofReal_pow, Complex.ofReal_re, Complex.ofReal_im]
    field_simp
    ring
  · simp [Complex.div_re, Complex.div_im, Complex.mul_re, Complex.mul_im]
    simp only [← Complex.ofReal_pow, Complex.ofReal_re, Complex.ofReal_im]
    field_simp
    ring

lemma K2 (N : ℝ) (hN : 0 < N) :
    (yv N : ℂ) * η1 + (ν0 - 1 / 8) / zv (θv N) = varpi1 / (√N : ℂ) := by
  have ha : √π ^ 2 = π := Real.sq_sqrt Real.pi_pos.le
  have hpa : 0 < √π := Real.sqrt_pos.2 Real.pi_pos
  have hpb : 0 < √N := Real.sqrt_pos.2 hN
  have h := K2_alg (√π) (√N) hpa hpb
  rw [ha] at h
  have hθ : 0 < θv N := θv_pos N hN
  rw [div_eq_mul_one_div (ν0 - 1 / 8), inv_zv _ hθ, θ_eq N, η1_eq, ν0_eq]
  have hy : (yv N : ℂ) * ((√2 / 2 : ℝ) : ℂ) = ((√π / √N / 2 : ℝ) : ℂ) := by
    rw [← Complex.ofReal_mul]; congr 1
    rw [← y_sqrt2 N hN]; ring
  rw [← mul_assoc, hy]
  unfold varpi1
  convert h using 3

lemma m0c_eq (θ : ℝ) (hθ : 0 < θ) :
    m0c θ = ((√(√2 * π / θ)) : ℂ) * Complex.exp ((π / 8 : ℝ) * Complex.I) := by
  unfold m0c
  have h1 : (π : ℂ) / (zv θ / 2) = ((√2 * π / θ : ℝ) : ℂ) * Complex.exp ((π / 4 : ℝ) * Complex.I) := by
    have e : (π : ℂ) / (zv θ / 2) = 2 * π * (1 / zv θ) := by
      have hz : zv θ ≠ 0 := by
        intro h; have := congrArg Complex.re h; rw [zv_re] at this; simp at this; linarith
      field_simp
    rw [e, inv_zv θ hθ]
    have hw : Complex.exp ((π / 4 : ℝ) * Complex.I) = ((√2 / 2 : ℝ) : ℂ) * (1 + Complex.I) := by
      have := ω_eq; unfold ω at this
      rw [show ((π / 4 : ℝ) : ℂ) * Complex.I = (π : ℂ) / 4 * Complex.I by push_cast; ring, this]
      ring
    rw [hw]
    have hθ' : (θ : ℂ) ≠ 0 := by exact_mod_cast hθ.ne'
    have h22 : (√2 : ℂ) * √2 = 2 := by
      rw [← Complex.ofReal_mul, Real.mul_self_sqrt (by norm_num)]; norm_num
    push_cast
    field_simp
    linear_combination (-(1 + Complex.I)) * h22
  rw [h1, cpow_half_polar _ _ (by positivity) (by linarith [Real.pi_pos]) (by linarith [Real.pi_pos])]
  congr 2; push_cast; ring

lemma K1_real (N : ℝ) (hN : 0 < N) :
    2 / π * (√(π / 2) * yv N * √(yv N)) * √(√2 * π / θv N) = √π / N := by
  have hy := yv_pos N hN
  have hθ := θv_pos N hN
  have h1 : 0 ≤ 2 / π * (√(π / 2) * yv N * √(yv N)) * √(√2 * π / θv N) := by positivity
  have h2 : 0 ≤ √π / N := by positivity
  refine (pow_left_inj₀ h1 h2 two_ne_zero).1 ?_
  have hy2 := y_sq N hN
  have hNy := Nyθ N hN
  rw [mul_pow, mul_pow, mul_pow, mul_pow, div_pow, Real.sq_sqrt (by positivity),
    Real.sq_sqrt hy.le, Real.sq_sqrt (by positivity), div_pow (√π) N 2, Real.sq_sqrt Real.pi_pos.le]
  have e : 2 ^ 2 / π ^ 2 * (π / 2 * yv N ^ 2 * yv N) * (√2 * π / θv N) = 2 * yv N ^ 2 / N := by
    rw [← hNy]; field_simp
  rw [e, hy2]; field_simp

lemma K1 (N : ℝ) (hN : 0 < N) :
    (2 / π : ℂ) * C0v N * m0c (θv N) * Complex.exp (-(Complex.I * π / 4)) =
      -Complex.I * ((√π / N : ℝ) : ℂ) := by
  rw [m0c_eq _ (θv_pos N hN), C0v, ← K1_real N hN]
  have hph : Complex.exp (-(3 * π / 8) * Complex.I) * Complex.exp ((π / 8 : ℝ) * Complex.I) *
      Complex.exp (-(Complex.I * π / 4)) = -Complex.I := by
    rw [← Complex.exp_add, ← Complex.exp_add]
    have : -(3 * (π : ℂ) / 8) * Complex.I + ((π / 8 : ℝ) : ℂ) * Complex.I + -(Complex.I * π / 4) =
        ((-(π / 2) : ℝ) : ℂ) * Complex.I := by push_cast; ring
    rw [this, Complex.exp_mul_I, ← Complex.ofReal_cos, ← Complex.ofReal_sin, Real.cos_neg,
      Real.sin_neg, Real.cos_pi_div_two, Real.sin_pi_div_two]
    simp
  push_cast at hph ⊢
  linear_combination (2 / π * (√(π / 2) * yv N * √(yv N)) * √(√2 * π / θv N) : ℂ) * hph

/-- `T = -i (√π/N) e^{-z} (1 + ϖ₁/√N)`. -/
noncomputable def Tv (N : ℝ) : ℂ :=
  -Complex.I * ((√π / N : ℝ) : ℂ) * Complex.exp (-zv (θv N)) * (1 + varpi1 / (√N : ℂ))

lemma T_im (N : ℝ) : (Tv N).im = -X N := by
  unfold Tv X
  have e : Complex.exp (-zv (θv N)) = ((Real.exp (-θv N) : ℝ) : ℂ) *
      Complex.exp (Complex.I * (θv N : ℂ)) := by
    rw [Complex.ofReal_exp, ← Complex.exp_add]; congr 1; unfold zv; push_cast; ring
  rw [e]
  unfold θv
  set w := Complex.exp (Complex.I * ((√(π * N) : ℝ) : ℂ)) * (1 + varpi1 / ((√N : ℝ) : ℂ))
  have : -Complex.I * ((√π / N : ℝ) : ℂ) * (((Real.exp (-√(π * N)) : ℝ) : ℂ) *
      Complex.exp (Complex.I * ((√(π * N) : ℝ) : ℂ))) * (1 + varpi1 / ((√N : ℝ) : ℂ)) =
      -Complex.I * ((√π / N * Real.exp (-√(π * N)) : ℝ) : ℂ) * w := by
    simp only [w]; push_cast; ring
  have h2 : ∀ (r : ℝ) (v : ℂ), (-Complex.I * (r : ℂ) * v).im = -(r * v.re) := by
    intro r v; simp [Complex.mul_im, Complex.mul_re]
  rw [this, h2]

lemma pc0_eq (N : ℝ) : pc N 0 = Complex.exp (-(Complex.I * π / 4)) * (1 + (yv N : ℂ) * η1) := by
  unfold pc; rw [tc_zero, tc_zero, φf_zero, χf_zero]; unfold η1; ring

lemma pc2_eq (N : ℝ) : pc N 2 = ν0 * Complex.exp (-(Complex.I * π / 4)) + (yv N : ℂ) * tc χf 2 := by
  unfold pc; rw [tc_φf_two]

lemma W7 (n : ℕ) (hn : 0 < n) :
    (2 / π : ℂ) * (C0v n * Complex.exp (-zv (θv n)) *
      (m0c (θv n) * (pc n 0 + pc n 2 / zv (θv n) - pc n 0 / (8 * zv (θv n))))) =
    Tv n + (2 / π : ℂ) * C0v n * m0c (θv n) * Complex.exp (-zv (θv n)) *
      ((yv n : ℂ) * (tc χf 2 - tc χf 0 / 8) / zv (θv n)) := by
  have hN : (0 : ℝ) < n := by exact_mod_cast hn
  have k1 := K1 n hN
  have k2 := K2 n hN
  have hχ0 : tc χf 0 = η1 * Complex.exp (-(Complex.I * π / 4)) := by
    rw [tc_zero, χf_zero]; unfold η1; ring
  rw [pc0_eq, pc2_eq, hχ0]
  unfold Tv
  set u := Complex.exp (-(Complex.I * π / 4))
  set E := Complex.exp (-zv (θv n))
  set z := zv (θv n)
  set y : ℂ := (yv n : ℂ)
  linear_combination (E * (1 + y * η1 + (ν0 - 1 / 8) / z)) * k1 +
    (E * (-Complex.I * ((√π / n : ℝ) : ℂ))) * k2

/-! ## 15. Region II: the assembly -/

lemma regII_core (n : ℕ) (hn : 0 < n) (he1 : e1v n < 1) (hu : uv n (-(1 / 4)) ^ 2 ≤ 2)
    (hS : 1 < Sl n (1 / 4)) (hf : Integrable (fσ n)) :
    ‖(2 / π : ℂ) * (∫ σ : ℝ, fσ n σ) - Tv n‖ ≤
      2 / π * (√2 * yv n * Real.exp (1 / 4) *
          Real.exp (-(θv n * (Real.exp (1 / 4) + Real.exp (-(1 / 4))) / 2) +
            γ0 * yv n * Real.exp (3 * (1 / 4))) / (Sl n (1 / 4) - 1) +
        √2 * yv n * Real.exp (-(1 / 4)) *
          Real.exp (-(θv n * (Real.exp (1 / 4) + Real.exp (-(1 / 4))) / 2) +
            γ0 * yv n * Real.exp (-(3 * (1 / 4)))) /
          (1 + θv n * ((Real.exp (1 / 4) - Real.exp (-(1 / 4))) / 2)) +
        Errb n * Real.exp (-θv n) * √(π / (θv n / 2))) +
      2 / π * (√(π / 2) * yv n * √(yv n)) * Real.exp (-θv n) *
        (32 * (Mφ + yv n * Mχ) * (∫ x : ℝ, x ^ 4 * Real.exp (-(θv n / 2) * x ^ 2)) +
          (0.0023 * θv n * ‖pc n 0‖ * (∫ x : ℝ, x ^ 6 * Real.exp (-(0.4963 * θv n) * x ^ 2)) +
            (0.0023 * θv n * ‖pc n 2‖ + 0.0035 * θv n ^ 2 * ‖pc n 0‖) *
              (∫ x : ℝ, x ^ 8 * Real.exp (-(0.4963 * θv n) * x ^ 2)) +
            0.0035 * θv n ^ 2 * ‖pc n 2‖ * (∫ x : ℝ, x ^ 10 * Real.exp (-(0.4963 * θv n) * x ^ 2)) +
            √2 * θv n / 24 * ‖pc n 2‖ * (∫ x : ℝ, x ^ 6 * Real.exp (-(θv n / 2) * x ^ 2))) +
          Real.exp (-(θv n / 64)) * (‖pc n 0‖ * (∫ x : ℝ, x ^ 0 * Real.exp (-(θv n / 4) * x ^ 2)) +
            ‖pc n 2‖ * (∫ x : ℝ, x ^ 2 * Real.exp (-(θv n / 4) * x ^ 2)) +
            √2 * θv n / 24 * ‖pc n 0‖ * (∫ x : ℝ, x ^ 4 * Real.exp (-(θv n / 4) * x ^ 2)))) +
      √π / n * Real.exp (-θv n) * (yv n * (‖tc χf 2‖ + ‖tc χf 0‖ / 8) / (√2 * θv n)) := by
  have hN : (0 : ℝ) < n := by exact_mod_cast hn
  have hθ := θv_pos n hN
  have hy := yv_pos n hN
  set θ := θv n
  set z := zv θ
  set E := Complex.exp (-z)
  set C0 := C0v n
  -- the pieces
  set L := ∫ σ in Iic (-(1 / 4 : ℝ)), fσ n σ
  set R := ∫ σ in Ioi (1 / 4 : ℝ), fσ n σ
  set ID := ∫ σ in -(1 / 4 : ℝ)..(1 / 4), (fσ n σ - C0 * E * (Kf θ σ * Pv n σ))
  set IKP := ∫ σ in -(1 / 4 : ℝ)..(1 / 4), Kf θ σ * Pv n σ
  set IK2 := ∫ σ in -(1 / 4 : ℝ)..(1 / 4), Kf θ σ * (pc n 0 + pc n 2 * (σ : ℂ) ^ 2)
  set IGQ := ∫ σ in -(1 / 4 : ℝ)..(1 / 4), Gf θ σ * Q4 n σ
  set IR := ∫ σ : ℝ, Gf θ σ * Q4 n σ
  have hJ : ∫ σ : ℝ, fσ n σ = L + (C0 * E * IKP + ID) + R := by
    rw [J_split3 _ hf (1 / 4) (by norm_num), window_split n hf (1 / 4)]
  have h4e : IK2 - IGQ = ∫ σ in -(1 / 4 : ℝ)..(1 / 4),
      (Kf θ σ * (pc n 0 + pc n 2 * (σ : ℂ) ^ 2) - Gf θ σ * Q4 n σ) := by
    rw [intervalIntegral.integral_sub]
    · exact ((Kf_cont θ).mul (by fun_prop)).intervalIntegrable _ _
    · exact (GQ_cont n).intervalIntegrable _ _
  have hIR : IR = m0c θ * (pc n 0 + pc n 2 / z - pc n 0 / (8 * z)) := W6 n hn
  have h7 := W7 n hn
  have key : (2 / π : ℂ) * (∫ σ : ℝ, fσ n σ) - Tv n =
      (2 / π : ℂ) * (L + R + ID) + (2 / π : ℂ) * C0 * E * ((IKP - IK2) + (IK2 - IGQ) + (IGQ - IR)) +
        (2 / π : ℂ) * C0 * m0c θ * E * ((yv n : ℂ) * (tc χf 2 - tc χf 0 / 8) / z) := by
    rw [hJ]
    have : (2 / π : ℂ) * (C0 * E * IR) = Tv n + (2 / π : ℂ) * C0 * m0c θ * E *
        ((yv n : ℂ) * (tc χf 2 - tc χf 0 / 8) / z) := by
      rw [hIR]; exact h7
    linear_combination this
  rw [key]
  -- the bounds
  have bL := left_tail n hn (1 / 4) hu hS
  have bR := right_tail n hn (1 / 4) (by norm_num)
  have bD := W1 n hn he1
  have b23 := W23 n hn
  have b4 := W4 n hn
  rw [← h4e] at b4
  have b5 := W5 n hn
  have hpi : ‖(2 / π : ℂ)‖ = 2 / π := by
    rw [norm_div, Complex.norm_real, Real.norm_eq_abs, abs_of_pos Real.pi_pos]; norm_num
  have hE : ‖E‖ = Real.exp (-θ) := by
    simp only [E, z]; rw [Complex.norm_exp]; congr 1; simp [zv]
  have hC0 : ‖C0‖ = √(π / 2) * yv n * √(yv n) := norm_C0v n hN
  have hΛ : ‖(2 / π : ℂ) * C0 * m0c θ‖ = √π / n := by
    have := K1 n hN
    have hu1 : ‖Complex.exp (-(Complex.I * π / 4))‖ = 1 := by
      rw [show -(Complex.I * π / 4) = ((-(π / 4) : ℝ) : ℂ) * Complex.I by push_cast; ring,
        Complex.norm_exp_ofReal_mul_I]
    have h2 := congrArg norm this
    simp only [norm_mul, hu1, mul_one, norm_neg, Complex.norm_I, one_mul, Complex.norm_real,
      Real.norm_eq_abs] at h2
    rw [norm_mul, norm_mul, h2, abs_of_pos (by positivity)]
  have hz : ‖z‖ = √2 * θ := norm_zv θ hθ.le
  have b7 : ‖(yv n : ℂ) * (tc χf 2 - tc χf 0 / 8) / z‖ ≤
      yv n * (‖tc χf 2‖ + ‖tc χf 0‖ / 8) / (√2 * θ) := by
    rw [norm_div, hz, norm_mul, Complex.norm_real, Real.norm_eq_abs, abs_of_pos hy]
    gcongr
    refine (norm_sub_le _ _).trans (le_of_eq ?_)
    rw [norm_div]; norm_num
  refine (norm_add_le _ _).trans (add_le_add ((norm_add_le _ _).trans (add_le_add ?_ ?_)) ?_)
  · rw [norm_mul, hpi]
    gcongr
    exact (norm_add_le _ _).trans (add_le_add ((norm_add_le _ _).trans (add_le_add bL bR)) bD)
  · rw [norm_mul, norm_mul, norm_mul, hpi, hC0, hE]
    gcongr
    refine (norm_add_le _ _).trans (add_le_add ((norm_add_le _ _).trans (add_le_add b23 b4)) ?_)
    exact b5
  · rw [norm_mul, norm_mul, hΛ, hE]
    gcongr

/-! ## 16. Region II: numerics -/

lemma facts2 (n : ℕ) (h : 2000 ≤ n) :
    79.26 ≤ θv n ∧ yv n ≤ 0.028025 ∧ 0 < yv n ∧ (n : ℝ) * yv n * √2 = θv n ∧
      (n : ℝ) * yv n ^ 2 = π / 2 ∧ θv n ^ 2 = π * n ∧ yv n * θv n * √2 = π := by
  have hN : (0 : ℝ) < n := by have : (2000 : ℝ) ≤ n := by exact_mod_cast h
                              linarith
  have hN1 : (2000 : ℝ) ≤ n := by exact_mod_cast h
  have hpi := Real.pi_gt_d6
  have hpi' := Real.pi_lt_d6
  have hθ := θv_pos n hN
  have hy := yv_pos n hN
  have h1 := Nyθ n hN
  have h2 := Ny2 n hN
  have h3 := θ_sq n hN
  refine ⟨?_, ?_, hy, h1, h2, h3, ?_⟩
  · nlinarith
  · have h4 := y_sq n hN
    have : π / (2 * n) ≤ 0.0007854 := by rw [div_le_iff₀ (by positivity)]; nlinarith
    nlinarith
  · have : yv n * θv n * √2 = yv n * (n * yv n * √2) * √2 := by rw [h1]
    rw [this]
    have h22 : √2 * √2 = 2 := Real.mul_self_sqrt (by norm_num)
    calc yv n * (n * yv n * √2) * √2 = (n * yv n ^ 2) * (√2 * √2) := by ring
      _ = π := by rw [h2, h22]; ring

lemma pow_exp_le (p : ℕ) (hp : 0 < p) (c x : ℝ) (hc : 0 < c) (hx : 0 < x) :
    x ^ p * Real.exp (-(c * x)) ≤ (p / c) ^ p * Real.exp (-(p : ℝ)) := by
  have hp' : (0 : ℝ) < p := by exact_mod_cast hp
  have h1 : c * x / p ≤ Real.exp (c * x / p - 1) := by linarith [Real.add_one_le_exp (c * x / p - 1)]
  have h2 : (c * x / p) ^ p ≤ Real.exp (c * x / p - 1) ^ p :=
    pow_le_pow_left₀ (by positivity) h1 p
  rw [← Real.exp_nat_mul] at h2
  have e1 : (p : ℝ) * (c * x / p - 1) = c * x - p := by field_simp
  rw [e1] at h2
  have e2 : x ^ p = (p / c) ^ p * (c * x / p) ^ p := by
    rw [← mul_pow]; congr 1; field_simp
  rw [e2, mul_assoc]
  apply mul_le_mul_of_nonneg_left _ (by positivity)
  calc (c * x / p) ^ p * Real.exp (-(c * x)) ≤ Real.exp (c * x - p) * Real.exp (-(c * x)) :=
        mul_le_mul_of_nonneg_right h2 (Real.exp_pos _).le
    _ = Real.exp (-(p : ℝ)) := by rw [← Real.exp_add]; ring_nf

lemma exp_one_lt' : Real.exp 1 < 2.7182818286 := Real.exp_one_lt_d9

lemma exp_one_gt' : 2.7182818283 < Real.exp 1 := Real.exp_one_gt_d9

lemma Mφ_le : Mφ ≤ 20.1 := by
  unfold Mφ
  have he := exp_one_lt'
  have hpi := Real.pi_lt_d6
  have h1 : 3 / 4 + π * Real.exp 1 / 4 ≤ 3 := by nlinarith [Real.exp_pos 1]
  have h2 : Real.exp 3 ≤ 20.1 := by
    have : Real.exp 3 = Real.exp 1 ^ 3 := by rw [← Real.exp_nat_mul]; norm_num
    rw [this]
    have h0 := Real.exp_pos 1
    calc Real.exp 1 ^ 3 ≤ 2.7182818286 ^ 3 := by gcongr
      _ ≤ 20.1 := by norm_num
  exact (Real.exp_le_exp.2 h1).trans h2

lemma Mχ_le : Mχ ≤ 72.1 := by
  unfold Mχ
  have he := exp_one_lt'
  have hpi := Real.pi_lt_d6
  have hM := Mφ_le
  have hM0 : 0 ≤ Mφ := by unfold Mφ; positivity
  have h12 : Real.exp (1 / 2) ≤ 1.6488 := by
    have h : Real.exp (1 / 2) ^ 2 = Real.exp 1 := by rw [← Real.exp_nat_mul]; norm_num
    nlinarith [Real.exp_pos (1 / 2)]
  have h32 : Real.exp (3 / 2) ≤ 4.482 := by
    have h : Real.exp (3 / 2) = Real.exp 1 * Real.exp (1 / 2) := by rw [← Real.exp_add]; norm_num
    rw [h]
    have h0 := Real.exp_pos (1 / 2)
    have h1 := Real.exp_pos 1
    calc Real.exp 1 * Real.exp (1 / 2) ≤ 2.7182818286 * 1.6488 := by gcongr
      _ ≤ 4.482 := by norm_num
  have h5 : π / 6 * Real.exp (3 / 2) ≤ 3.141593 / 6 * 4.482 := by gcongr
  have : 3 / 4 * Real.exp (1 / 2) + π / 6 * Real.exp (3 / 2) ≤ 3.587 := by linarith
  nlinarith [Real.exp_pos (1 / 2), Real.exp_pos (3 / 2)]

lemma exp_neg_two_le : Real.exp (-2) ≤ 0.13534 := by
  rw [Real.exp_neg, inv_le_comm₀ (Real.exp_pos _) (by norm_num)]
  have h : Real.exp 2 = Real.exp 1 ^ 2 := by rw [← Real.exp_nat_mul]; norm_num
  rw [h]; have := exp_one_gt'; nlinarith

lemma exp_neg_one_le : Real.exp (-1) ≤ 0.36788 := by
  rw [Real.exp_neg, inv_le_comm₀ (Real.exp_pos _) (by norm_num)]
  have := exp_one_gt'; linarith

lemma θ2_sup (θ : ℝ) (hθ : 0 < θ) (c : ℝ) (hc : 0 < c) :
    θ ^ 2 * Real.exp (-(c * θ)) ≤ (2 / c) ^ 2 * 0.13534 := by
  have := pow_exp_le 2 (by norm_num) c θ hc hθ
  push_cast at this
  exact this.trans (mul_le_mul_of_nonneg_left exp_neg_two_le (by positivity))

lemma θ1_sup (θ : ℝ) (hθ : 0 < θ) (c : ℝ) (hc : 0 < c) :
    θ * Real.exp (-(c * θ)) ≤ (1 / c) * 0.36788 := by
  have := pow_exp_le 1 (by norm_num) c θ hc hθ
  push_cast at this
  simp only [pow_one] at this
  exact this.trans (mul_le_mul_of_nonneg_left exp_neg_one_le (by positivity))

lemma termLR (n : ℕ) (h : 2000 ≤ n) :
    2 / π * (√2 * yv n * Real.exp (1 / 4) *
          Real.exp (-(θv n * (Real.exp (1 / 4) + Real.exp (-(1 / 4))) / 2) +
            γ0 * yv n * Real.exp (3 * (1 / 4))) / (Sl n (1 / 4) - 1) +
        √2 * yv n * Real.exp (-(1 / 4)) *
          Real.exp (-(θv n * (Real.exp (1 / 4) + Real.exp (-(1 / 4))) / 2) +
            γ0 * yv n * Real.exp (-(3 * (1 / 4)))) /
          (1 + θv n * ((Real.exp (1 / 4) - Real.exp (-(1 / 4))) / 2))) ≤
      1000 / (n : ℝ) ^ 2 * Real.exp (-θv n) := by
  obtain ⟨hθ1, hy1, hy0, hnyθ, hny2, hθ2, hyθ⟩ := facts2 n h
  have hN : (0 : ℝ) < n := by have : (2000 : ℝ) ≤ n := by exact_mod_cast h
                              linarith
  have hn : 0 < n := by omega
  set θ := θv n
  set y := yv n
  have hθ : 0 < θ := by linarith
  have hpi := Real.pi_gt_d6
  have hpi' := Real.pi_lt_d6
  have e14 := exp_quarter_lt
  have e14' := exp_quarter_gt
  have em14 := exp_neg_quarter_lt
  have em14' := exp_neg_quarter_gt
  have e34 := exp_three_quarter_lt
  have hg := γ0_lt
  have hg' := γ0_pos
  have hs2 := sqrt2_lt
  -- cosh(1/4) ≥ 1 + 1/32
  have hc : 1 + 1 / 32 ≤ (Real.exp (1 / 4) + Real.exp (-(1 / 4))) / 2 := by
    have := cosh_quad (1 / 4); norm_num at this ⊢; linarith
  have hsh : 0.2525 ≤ (Real.exp (1 / 4) - Real.exp (-(1 / 4))) / 2 := by linarith
  -- the γ₀ corrections
  have hgy : γ0 * y * Real.exp (3 * (1 / 4)) ≤ 0.022 := by
    have h1 : γ0 * y ≤ 0.37025 * 0.028025 := by gcongr
    have h2 : γ0 * y * Real.exp (3 * (1 / 4)) ≤ (0.37025 * 0.028025) * 2.1174 := by
      gcongr
    linarith
  have hgy0 : 0 ≤ γ0 * y * Real.exp (3 * (1 / 4)) := by positivity
  have hgy' : γ0 * y * Real.exp (-(3 * (1 / 4))) ≤ 0.022 := by
    have : Real.exp (-(3 * (1 / 4))) ≤ Real.exp (3 * (1 / 4)) := Real.exp_le_exp.2 (by norm_num)
    have := mul_le_mul_of_nonneg_left this (by positivity : 0 ≤ γ0 * y)
    linarith
  have heg : Real.exp 0.022 ≤ 1.0224 := by
    have := exp_le_quad 0.022 (by norm_num [abs_of_pos]); norm_num at this ⊢; linarith
  -- the slope
  have hu1 : uv n (-(1 / 4)) ≤ 0.051 := by
    rw [uv_neg]
    have h1 : √2 * yv n ≤ 1.41422 * 0.028025 := by gcongr
    have h2 : √2 * yv n * Real.exp (1 / 4) ≤ (1.41422 * 0.028025) * 1.2841 := by gcongr
    linarith
  have hu0 := uv_nonneg n (-(1 / 4)) hN
  have hS := Sl_ge n hn (1 / 4)
  have hSθ : 0.2516 * θ ≤ Sl n (1 / 4) := by
    refine le_trans ?_ hS
    rw [mul_comm (0.2516 : ℝ)]
    apply mul_le_mul_of_nonneg_left _ hθ.le
    have hu2 : uv n (-(1 / 4)) ^ 2 ≤ 0.051 ^ 2 := pow_le_pow_left₀ hu0 hu1 2
    have h3 : Real.exp (1 / 4) / 2 * (uv n (-(1 / 4)) ^ 2 / 2) ≤ 1.2841 / 2 * (0.051 ^ 2 / 2) := by
      gcongr
    have h4 : Real.exp (1 / 4) / 2 * (1 - uv n (-(1 / 4)) ^ 2 / 2) =
        Real.exp (1 / 4) / 2 - Real.exp (1 / 4) / 2 * (uv n (-(1 / 4)) ^ 2 / 2) := by ring
    rw [h4]; norm_num at h3 ⊢; linarith
  have hS1 : 0.2389 * θ ≤ Sl n (1 / 4) - 1 := by linarith
  -- the exponentials
  have hE1 : Real.exp (-(θ * (Real.exp (1 / 4) + Real.exp (-(1 / 4))) / 2) +
      γ0 * y * Real.exp (3 * (1 / 4))) ≤ Real.exp (-θ) * (Real.exp (-(1 / 32 * θ)) * 1.0224) := by
    calc Real.exp (-(θ * (Real.exp (1 / 4) + Real.exp (-(1 / 4))) / 2) +
          γ0 * y * Real.exp (3 * (1 / 4))) ≤ Real.exp (-θ + -(1 / 32 * θ) + 0.022) := by
          apply Real.exp_le_exp.2; nlinarith
      _ = Real.exp (-θ) * (Real.exp (-(1 / 32 * θ)) * Real.exp 0.022) := by
          rw [← Real.exp_add, ← Real.exp_add]; ring_nf
      _ ≤ Real.exp (-θ) * (Real.exp (-(1 / 32 * θ)) * 1.0224) := by gcongr
  have hE2 : Real.exp (-(θ * (Real.exp (1 / 4) + Real.exp (-(1 / 4))) / 2) +
      γ0 * y * Real.exp (-(3 * (1 / 4)))) ≤ Real.exp (-θ) * (Real.exp (-(1 / 32 * θ)) * 1.0224) := by
    calc Real.exp (-(θ * (Real.exp (1 / 4) + Real.exp (-(1 / 4))) / 2) +
          γ0 * y * Real.exp (-(3 * (1 / 4)))) ≤ Real.exp (-θ + -(1 / 32 * θ) + 0.022) := by
          apply Real.exp_le_exp.2; nlinarith
      _ = Real.exp (-θ) * (Real.exp (-(1 / 32 * θ)) * Real.exp 0.022) := by
          rw [← Real.exp_add, ← Real.exp_add]; ring_nf
      _ ≤ Real.exp (-θ) * (Real.exp (-(1 / 32 * θ)) * 1.0224) := by gcongr
  have hsy : √2 * y = θ / n := by field_simp; linarith
  have hsup := θ2_sup θ hθ (1 / 32) (by norm_num)
  have hnθ : (n : ℝ) = θ ^ 2 / π := by field_simp; linarith
  set X := Real.exp (-(1 / 32 * θ))
  have hX : 0 < X := Real.exp_pos _
  have hEθ := Real.exp_pos (-θ)
  -- the two terms
  have t1 : √2 * y * Real.exp (1 / 4) *
      Real.exp (-(θ * (Real.exp (1 / 4) + Real.exp (-(1 / 4))) / 2) + γ0 * y * Real.exp (3 * (1 / 4))) /
        (Sl n (1 / 4) - 1) ≤ (1.2841 * 1.0224 / 0.2389) / n * (Real.exp (-θ) * X) := by
    rw [hsy, div_le_iff₀ (by linarith)]
    calc θ / n * Real.exp (1 / 4) * Real.exp (-(θ * (Real.exp (1 / 4) + Real.exp (-(1 / 4))) / 2) +
          γ0 * y * Real.exp (3 * (1 / 4))) ≤ θ / n * 1.2841 * (Real.exp (-θ) * (X * 1.0224)) := by
          gcongr
      _ ≤ (1.2841 * 1.0224 / 0.2389) / n * (Real.exp (-θ) * X) * (Sl n (1 / 4) - 1) := by
          have : θ / n * 1.2841 * (Real.exp (-θ) * (X * 1.0224)) =
              (1.2841 * 1.0224 / 0.2389) / n * (Real.exp (-θ) * X) * (0.2389 * θ) := by
            field_simp
          rw [this]; gcongr
  have t2 : √2 * y * Real.exp (-(1 / 4)) *
      Real.exp (-(θ * (Real.exp (1 / 4) + Real.exp (-(1 / 4))) / 2) + γ0 * y * Real.exp (-(3 * (1 / 4)))) /
        (1 + θ * ((Real.exp (1 / 4) - Real.exp (-(1 / 4))) / 2)) ≤
        (0.77883 * 1.0224 / 0.2525) / n * (Real.exp (-θ) * X) := by
    have hpos : 0 < 1 + θ * ((Real.exp (1 / 4) - Real.exp (-(1 / 4))) / 2) := by
      have := mul_nonneg hθ.le (le_trans (by norm_num) hsh); linarith
    rw [hsy, div_le_iff₀ hpos]
    calc θ / n * Real.exp (-(1 / 4)) * Real.exp (-(θ * (Real.exp (1 / 4) + Real.exp (-(1 / 4))) / 2) +
          γ0 * y * Real.exp (-(3 * (1 / 4)))) ≤ θ / n * 0.77883 * (Real.exp (-θ) * (X * 1.0224)) := by
          gcongr
      _ ≤ (0.77883 * 1.0224 / 0.2525) / n * (Real.exp (-θ) * X) *
          (1 + θ * ((Real.exp (1 / 4) - Real.exp (-(1 / 4))) / 2)) := by
          have : θ / n * 0.77883 * (Real.exp (-θ) * (X * 1.0224)) =
              (0.77883 * 1.0224 / 0.2525) / n * (Real.exp (-θ) * X) * (0.2525 * θ) := by
            field_simp
          rw [this]
          apply mul_le_mul_of_nonneg_left _ (by positivity)
          have := mul_le_mul_of_nonneg_left hsh hθ.le
          linarith
  have hfin : 2 / π * ((1.2841 * 1.0224 / 0.2389) / n * (Real.exp (-θ) * X) +
      (0.77883 * 1.0224 / 0.2525) / n * (Real.exp (-θ) * X)) ≤ 1000 / (n : ℝ) ^ 2 * Real.exp (-θ) := by
    have e : 2 / π * ((1.2841 * 1.0224 / 0.2389) / n * (Real.exp (-θ) * X) +
        (0.77883 * 1.0224 / 0.2525) / n * (Real.exp (-θ) * X)) =
        (2 / π * (1.2841 * 1.0224 / 0.2389 + 0.77883 * 1.0224 / 0.2525) * (n * X)) / (n : ℝ) ^ 2 *
          Real.exp (-θ) := by field_simp
    rw [e]
    gcongr
    have hnX : (n : ℝ) * X ≤ (2 / (1 / 32)) ^ 2 * 0.13534 / π := by
      rw [hnθ, div_mul_eq_mul_div, div_le_div_iff_of_pos_right Real.pi_pos]; exact hsup
    calc 2 / π * (1.2841 * 1.0224 / 0.2389 + 0.77883 * 1.0224 / 0.2525) * (n * X)
        ≤ 2 / π * (1.2841 * 1.0224 / 0.2389 + 0.77883 * 1.0224 / 0.2525) *
          ((2 / (1 / 32)) ^ 2 * 0.13534 / π) := by gcongr
      _ ≤ 2 / 3.14159 * (1.2841 * 1.0224 / 0.2389 + 0.77883 * 1.0224 / 0.2525) *
          ((2 / (1 / 32)) ^ 2 * 0.13534 / 3.14159) := by gcongr <;> linarith
      _ ≤ 1000 := by norm_num
  calc _ ≤ 2 / π * ((1.2841 * 1.0224 / 0.2389) / n * (Real.exp (-θ) * X) +
      (0.77883 * 1.0224 / 0.2525) / n * (Real.exp (-θ) * X)) := by gcongr
    _ ≤ _ := hfin

lemma keyK (n : ℕ) (h : 2000 ≤ n) (a : ℝ) (ha : 0 < a) :
    (n : ℝ) ^ 2 * (yv n * √(yv n)) * √(π / (a * θv n)) = θv n ^ 2 / (√2 * √(√2 * a)) := by
  obtain ⟨hθ1, hy1, hy0, hnyθ, hny2, hθ2, hyθ⟩ := facts2 n h
  have hθ : 0 < θv n := by linarith
  have h1 : √(yv n) * √(π / (a * θv n)) = π / (θv n * √(√2 * a)) := by
    rw [← Real.sqrt_mul hy0.le]
    rw [show π / (θv n * √(√2 * a)) = √((π / (θv n * √(√2 * a))) ^ 2) by
      rw [Real.sqrt_sq (by positivity)]]
    congr 1
    rw [div_pow, mul_pow, Real.sq_sqrt (by positivity)]
    have hy : yv n = π / (√2 * θv n) := by
      field_simp; linarith
    rw [hy]; field_simp
  have hn : (n : ℝ) ^ 2 * π ^ 2 = θv n ^ 4 := by
    rw [show θv n ^ 4 = (θv n ^ 2) ^ 2 by ring, hθ2]; ring
  have hy : yv n = π / (√2 * θv n) := by field_simp; linarith
  calc (n : ℝ) ^ 2 * (yv n * √(yv n)) * √(π / (a * θv n))
      = (n : ℝ) ^ 2 * yv n * (√(yv n) * √(π / (a * θv n))) := by ring
    _ = θv n ^ 2 / (√2 * √(√2 * a)) := by
      rw [h1, hy]
      have h2 : (0 : ℝ) < √2 := sqrt2_pos
      have h3 : (0 : ℝ) < √(√2 * a) := by positivity
      field_simp
      linear_combination hn

lemma rpow72 (b y : ℝ) (hb : 0 < b) (hy : 0 < y) :
    (b * y) ^ ((2 : ℕ) + 3 / 2 : ℝ) = b ^ 3 * √b * y ^ 2 * (y * √y) := by
  have hby : 0 < b * y := by positivity
  rw [show ((2 : ℕ) : ℝ) + 3 / 2 = ((3 : ℕ) : ℝ) + 1 / 2 by norm_num, Real.rpow_add hby,
    Real.rpow_natCast, ← Real.sqrt_eq_rpow, Real.sqrt_mul hb.le]
  ring

lemma Errb_le (n : ℕ) (h : 2000 ≤ n) : Errb n ≤ 65.5 * yv n ^ 2 * (yv n * √(yv n)) := by
  obtain ⟨hθ1, hy1, hy0, hnyθ, hny2, hθ2, hyθ⟩ := facts2 n h
  have hN : (0 : ℝ) < n := by have : (2000 : ℝ) ≤ n := by exact_mod_cast h
                              linarith
  have hθ : 0 < θv n := by linarith
  have hpi := Real.pi_gt_d6
  have hpi' := Real.pi_lt_d6
  have e14 := exp_quarter_lt
  have e14' := exp_quarter_gt
  have e34 := exp_three_quarter_lt
  set y := yv n
  -- e₁, R₄, B
  have he1 : e1v n ≤ 0.035988 := by
    unfold e1v
    have : y * Real.exp (1 / 4) ≤ 0.028025 * 1.2841 := by gcongr
    linarith
  have he10 : 0 ≤ e1v n := by unfold e1v; positivity
  have hR4 : R4v n ≤ 1.108 * y ^ 2 := by
    unfold R4v e1v
    have h4 : (y * Real.exp (1 / 4)) ^ 4 = y ^ 4 * Real.exp 1 := by
      rw [mul_pow, ← Real.exp_nat_mul]; norm_num
    rw [h4]
    have hye : y * Real.exp (1 / 4) ≤ 0.035988 := he1
    have hinv : (1 - y * Real.exp (1 / 4))⁻¹ ≤ 1.0374 := by
      rw [inv_le_comm₀ (by linarith) (by norm_num)]; norm_num; linarith
    have hpos1 : 0 ≤ (1 - y * Real.exp (1 / 4))⁻¹ := inv_nonneg.2 (by linarith)
    have hny4 : (n : ℝ) * y ^ 4 = π / 2 * y ^ 2 := by
      rw [show y ^ 4 = y ^ 2 * y ^ 2 by ring, ← mul_assoc, hny2]
    have he := exp_one_lt'
    calc (n : ℝ) * (y ^ 4 * Real.exp 1 * (1 - y * Real.exp (1 / 4))⁻¹ / 4)
        = (n * y ^ 4) * Real.exp 1 * (1 - y * Real.exp (1 / 4))⁻¹ / 4 := by ring
      _ = π / 2 * Real.exp 1 * (1 - y * Real.exp (1 / 4))⁻¹ / 4 * y ^ 2 := by rw [hny4]; ring
      _ ≤ 3.141593 / 2 * 2.7182818286 * 1.0374 / 4 * y ^ 2 := by gcongr
      _ ≤ 1.108 * y ^ 2 := by gcongr; norm_num
  have hR40 : 0 ≤ R4v n := by
    unfold R4v; have : 0 ≤ (1 - e1v n)⁻¹ := inv_nonneg.2 (by linarith); positivity
  have hB : Bbv n ≤ 1.14 * y := by
    unfold Bbv
    have h1 : π / 6 * y * Real.exp (3 / 4) ≤ 3.141593 / 6 * 2.1174 * y := by
      have : Real.exp (3 / 4) ≤ 2.1174 := by
        have := exp_three_quarter_lt; norm_num at this ⊢; linarith
      calc π / 6 * y * Real.exp (3 / 4) = π / 6 * Real.exp (3 / 4) * y := by ring
        _ ≤ 3.141593 / 6 * 2.1174 * y := by gcongr
    have h2 : 1.108 * y ^ 2 ≤ 1.108 * 0.028025 * y := by
      rw [sq, ← mul_assoc]; gcongr
    have : (3.141593 / 6 * 2.1174 + 1.108 * 0.028025 : ℝ) ≤ 1.14 := by norm_num
    nlinarith
  have hB0 : 0 ≤ Bbv n := by unfold Bbv; positivity
  have hBs : Bbv n ≤ 0.032 := by nlinarith
  have heB : Real.exp (Bbv n) ≤ 1.033 := by
    have := exp_le_quad (Bbv n) (by rw [abs_le]; constructor <;> linarith)
    nlinarith
  have he38 : Real.exp (3 / 8) ≤ 1.481 := by
    have := exp_le_quad (3 / 8) (by norm_num [abs_of_pos]); norm_num at this ⊢; linarith
  -- the bracket
  have hbr : Bbv n ^ 2 * Real.exp (Bbv n) + R4v n + 3 / 4 * e1v n * (Bbv n * Real.exp (Bbv n)) ≤
      3.586 * y ^ 2 := by
    have hB2 : Bbv n ^ 2 ≤ (1.14 * y) ^ 2 := pow_le_pow_left₀ hB0 hB 2
    have h1 : Bbv n ^ 2 * Real.exp (Bbv n) ≤ (1.14 * y) ^ 2 * 1.033 :=
      mul_le_mul hB2 heB (Real.exp_pos _).le (by positivity)
    have he1y : e1v n ≤ 1.2841 * y := by
      unfold e1v; rw [mul_comm]; gcongr
    have h3 : 3 / 4 * e1v n * (Bbv n * Real.exp (Bbv n)) ≤ 3 / 4 * (1.2841 * y) * (1.14 * y * 1.033) := by
      gcongr
    nlinarith
  -- the rH part
  have hW : (√2 * e1v n) ^ ((2 : ℕ) + 3 / 2 : ℝ) ≤ 8.072 * y ^ 2 * (y * √y) := by
    have e : √2 * e1v n = (√2 * Real.exp (1 / 4)) * y := by unfold e1v; ring
    rw [e, rpow72 _ _ (by positivity) hy0]
    have hb : √2 * Real.exp (1 / 4) ≤ 1.81601 := by
      have := sqrt2_lt
      calc √2 * Real.exp (1 / 4) ≤ 1.41422 * 1.2841 := by gcongr
        _ ≤ 1.81601 := by norm_num
    have hb0 : 0 < √2 * Real.exp (1 / 4) := by positivity
    have hsb : √(√2 * Real.exp (1 / 4)) ≤ 1.3476 := by
      rw [Real.sqrt_le_left (by norm_num)]; nlinarith
    have : (√2 * Real.exp (1 / 4)) ^ 3 * √(√2 * Real.exp (1 / 4)) ≤ 1.81601 ^ 3 * 1.3476 := by
      gcongr
    have h2 : (1.81601 : ℝ) ^ 3 * 1.3476 ≤ 8.072 := by norm_num
    have hpos : 0 ≤ y ^ 2 * (y * √y) := by positivity
    nlinarith
  have hs2π : √(2 * π) ≤ 2.5067 := by
    rw [Real.sqrt_le_left (by norm_num)]; nlinarith
  have hspi2 : √(π / 2) ≤ 1.2534 := by
    rw [Real.sqrt_le_left (by norm_num)]; nlinarith
  show Errb n ≤ 65.5 * y ^ 2 * (y * √y)
  · unfold Errb
    have hpos : 0 ≤ y ^ 2 * (y * √y) := by positivity
    have h1 : √(π / 2) * y * √y * Real.exp (3 / 8) *
        (Bbv n ^ 2 * Real.exp (Bbv n) + R4v n + 3 / 4 * e1v n * (Bbv n * Real.exp (Bbv n))) ≤
        1.2534 * 1.481 * 3.586 * (y ^ 2 * (y * √y)) := by
      calc √(π / 2) * y * √y * Real.exp (3 / 8) *
          (Bbv n ^ 2 * Real.exp (Bbv n) + R4v n + 3 / 4 * e1v n * (Bbv n * Real.exp (Bbv n)))
          ≤ 1.2534 * y * √y * 1.481 * (3.586 * y ^ 2) := by gcongr
        _ = 1.2534 * 1.481 * 3.586 * (y ^ 2 * (y * √y)) := by ring
    have h2 : 45 / 16 * √(2 * π) * (√2 * e1v n) ^ ((2 : ℕ) + 3 / 2 : ℝ) * Real.exp (Bbv n) ≤
        45 / 16 * 2.5067 * (8.072 * y ^ 2 * (y * √y)) * 1.033 := by gcongr
    linarith
lemma sqrt_aux : (1.189 : ℝ) ≤ √2 * √(√2 * (1 / 2)) := by
  rw [← Real.sqrt_mul (by norm_num), Real.le_sqrt (by norm_num) (by positivity)]
  have := sqrt2_gt; norm_num; linarith

lemma termD (n : ℕ) (h : 2000 ≤ n) :
    2 / π * (Errb n * Real.exp (-θv n) * √(π / (θv n / 2))) ≤ 200 / (n : ℝ) ^ 2 * Real.exp (-θv n) := by
  obtain ⟨hθ1, hy1, hy0, hnyθ, hny2, hθ2, hyθ⟩ := facts2 n h
  have hN : (0 : ℝ) < n := by have : (2000 : ℝ) ≤ n := by exact_mod_cast h
                              linarith
  have hθ : 0 < θv n := by linarith
  have hpi := Real.pi_lt_d6
  have hpi' := Real.pi_gt_d6
  set y := yv n
  have hErr := Errb_le n h
  -- assemble
  have hK := keyK n h (1 / 2) (by norm_num)
  rw [show (1 / 2 : ℝ) * θv n = θv n / 2 by ring] at hK
  have h2r := sqrt_aux
  have hyθ2 : y ^ 2 * θv n ^ 2 = π ^ 2 / 2 := by
    have h22 : √2 ^ 2 = 2 := Real.sq_sqrt (by norm_num)
    have : (y * θv n * √2) ^ 2 = π ^ 2 := by rw [hyθ]
    rw [mul_pow, mul_pow, h22] at this; linarith
  have hE := Real.exp_pos (-θv n)
  have hsq : 0 < √(π / (θv n / 2)) := by positivity
  have hfin : 2 / π * ((n : ℝ) ^ 2 * Errb n * √(π / (θv n / 2))) ≤ 200 := by
    have h3 : (n : ℝ) ^ 2 * Errb n * √(π / (θv n / 2)) ≤
        65.5 * y ^ 2 * ((n : ℝ) ^ 2 * (y * √y) * √(π / (θv n / 2))) := by
      have := mul_le_mul_of_nonneg_right hErr hsq.le
      have h4 := mul_le_mul_of_nonneg_left this (by positivity : (0 : ℝ) ≤ (n : ℝ) ^ 2)
      calc (n : ℝ) ^ 2 * Errb n * √(π / (θv n / 2)) = (n : ℝ) ^ 2 * (Errb n * √(π / (θv n / 2))) := by ring
        _ ≤ (n : ℝ) ^ 2 * (65.5 * y ^ 2 * (y * √y) * √(π / (θv n / 2))) := h4
        _ = _ := by ring
    rw [hK] at h3
    have h5 : 65.5 * y ^ 2 * (θv n ^ 2 / (√2 * √(√2 * (1 / 2)))) ≤ 65.5 * (π ^ 2 / 2) / 1.189 := by
      calc 65.5 * y ^ 2 * (θv n ^ 2 / (√2 * √(√2 * (1 / 2))))
          = 65.5 * (y ^ 2 * θv n ^ 2) / (√2 * √(√2 * (1 / 2))) := by ring
        _ = 65.5 * (π ^ 2 / 2) / (√2 * √(√2 * (1 / 2))) := by rw [hyθ2]
        _ ≤ 65.5 * (π ^ 2 / 2) / 1.189 := by gcongr
    calc 2 / π * ((n : ℝ) ^ 2 * Errb n * √(π / (θv n / 2))) ≤ 2 / π * (65.5 * (π ^ 2 / 2) / 1.189) := by
          gcongr; exact h3.trans h5
      _ = 65.5 * π / 1.189 := by field_simp
      _ ≤ 200 := by rw [div_le_iff₀ (by norm_num)]; linarith [Real.pi_lt_d6]
  calc 2 / π * (Errb n * Real.exp (-θv n) * √(π / (θv n / 2)))
      = (2 / π * ((n : ℝ) ^ 2 * Errb n * √(π / (θv n / 2)))) / (n : ℝ) ^ 2 * Real.exp (-θv n) := by
        field_simp
    _ ≤ 200 / (n : ℝ) ^ 2 * Real.exp (-θv n) := by gcongr

lemma norm_η1 : ‖η1‖ ≤ 1.2736 := by
  unfold η1
  rw [norm_neg]
  refine (norm_add_le _ _).trans ?_
  rw [norm_mul, norm_mul, norm_pow, norm_ωb, one_pow, mul_one, mul_one]
  have h1 : ‖(3 / 4 : ℂ)‖ = 3 / 4 := by norm_num
  have h2 : ‖((π : ℂ) / 6)‖ = π / 6 := by
    rw [norm_div, Complex.norm_real, Real.norm_eq_abs, abs_of_pos Real.pi_pos]; norm_num
  rw [h1, h2]; linarith [Real.pi_lt_d6]

lemma norm_ν0 : ‖ν0‖ ≤ 4.1 := by
  rw [ν0_eq]
  refine (norm_add_le _ _).trans ?_
  rw [norm_mul, Complex.norm_I, mul_one, Complex.norm_real, Complex.norm_real, Real.norm_eq_abs,
    Real.norm_eq_abs]
  have hpi := Real.pi_gt_d6
  have hpi' := Real.pi_lt_d6
  have h1 : |(9 - π ^ 2) / 8| ≤ 0.11 := by rw [abs_le]; constructor <;> nlinarith
  have h2 : |-(5 * π / 4)| ≤ 3.93 := by rw [abs_le]; constructor <;> nlinarith
  linarith

lemma norm_u : ‖Complex.exp (-(Complex.I * π / 4))‖ = 1 := by
  rw [show -(Complex.I * π / 4) = ((-(π / 4) : ℝ) : ℂ) * Complex.I by push_cast; ring,
    Complex.norm_exp_ofReal_mul_I]

lemma norm_pc0 (n : ℕ) (h : 2000 ≤ n) : ‖pc n 0‖ ≤ 1.0357 := by
  obtain ⟨-, hy1, hy0, -⟩ := facts2 n h
  rw [pc0_eq, norm_mul, norm_u, one_mul]
  refine (norm_add_le _ _).trans ?_
  rw [norm_one, norm_mul, Complex.norm_real, Real.norm_eq_abs, abs_of_pos hy0]
  have := norm_η1
  have : yv n * ‖η1‖ ≤ 0.028025 * 1.2736 := by gcongr
  linarith

lemma norm_tcχ (k : ℕ) : ‖tc χf k‖ ≤ 72.1 * 2 ^ k :=
  (tc_le χf χf_diff Mχ χf_sphere k).trans (by gcongr; exact Mχ_le)

lemma norm_pc2 (n : ℕ) (h : 2000 ≤ n) : ‖pc n 2‖ ≤ 12.19 := by
  obtain ⟨-, hy1, hy0, -⟩ := facts2 n h
  rw [pc2_eq]
  refine (norm_add_le _ _).trans ?_
  rw [norm_mul, norm_u, mul_one, norm_mul, Complex.norm_real, Real.norm_eq_abs, abs_of_pos hy0]
  have h1 := norm_ν0
  have h2 := norm_tcχ 2
  have : yv n * ‖tc χf 2‖ ≤ 0.028025 * (72.1 * 2 ^ 2) := by gcongr
  norm_num at this ⊢; linarith

lemma term7 (n : ℕ) (h : 2000 ≤ n) :
    √π / n * Real.exp (-θv n) * (yv n * (‖tc χf 2‖ + ‖tc χf 0‖ / 8) / (√2 * θv n)) ≤
      300 / (n : ℝ) ^ 2 * Real.exp (-θv n) := by
  obtain ⟨hθ1, hy1, hy0, hnyθ, hny2, hθ2, hyθ⟩ := facts2 n h
  have hN : (0 : ℝ) < n := by have : (2000 : ℝ) ≤ n := by exact_mod_cast h
                              linarith
  have hθ : 0 < θv n := by linarith
  have e : yv n / (√2 * θv n) = 1 / (2 * n) := by
    rw [← hnyθ]
    have h22 : √2 * √2 = 2 := Real.mul_self_sqrt (by norm_num)
    have h2 : (0 : ℝ) < √2 := sqrt2_pos
    field_simp
    linear_combination -1 * h22
  have h2 := norm_tcχ 2
  have h0 := norm_tcχ 0
  have hsp := Q05.sqrt_pi_lt
  have hE := Real.exp_pos (-θv n)
  calc √π / n * Real.exp (-θv n) * (yv n * (‖tc χf 2‖ + ‖tc χf 0‖ / 8) / (√2 * θv n))
      = √π * (‖tc χf 2‖ + ‖tc χf 0‖ / 8) / 2 / (n : ℝ) ^ 2 * Real.exp (-θv n) := by
        rw [mul_div_right_comm, e]; field_simp
    _ ≤ 1.7724539 * (72.1 * 2 ^ 2 + 72.1 * 2 ^ 0 / 8) / 2 / (n : ℝ) ^ 2 * Real.exp (-θv n) := by
        gcongr
    _ ≤ 300 / (n : ℝ) ^ 2 * Real.exp (-θv n) := by gcongr; norm_num

lemma nG0_mom (n : ℕ) (h : 2000 ≤ n) (a : ℝ) (ha : 0 < a) (d : ℝ) (j : ℕ) :
    (n : ℝ) ^ 2 * (√(π / 2) * yv n * √(yv n)) * (d * √(π / (a * θv n)) / (2 * (a * θv n)) ^ j) =
      √(π / 2) * d / (√2 * √(√2 * a)) / (2 * a) ^ j * θv n ^ 2 / θv n ^ j := by
  have hK := keyK n h a ha
  obtain ⟨hθ1, -⟩ := facts2 n h
  have hθ : 0 < θv n := by linarith
  have e : (n : ℝ) ^ 2 * (√(π / 2) * yv n * √(yv n)) * (d * √(π / (a * θv n)) / (2 * (a * θv n)) ^ j) =
      √(π / 2) * d * ((n : ℝ) ^ 2 * (yv n * √(yv n)) * √(π / (a * θv n))) / (2 * (a * θv n)) ^ j := by
    ring
  rw [e, hK, mul_pow]
  have : (0 : ℝ) < √2 * √(√2 * a) := by positivity
  field_simp
  ring

lemma D_half : (1.189 : ℝ) ≤ √2 * √(√2 * (1 / 2)) := sqrt_aux

lemma D_quarter : (0.8408 : ℝ) ≤ √2 * √(√2 * (1 / 4)) := by
  rw [← Real.sqrt_mul (by norm_num), Real.le_sqrt (by norm_num) (by positivity)]
  have := sqrt2_gt; norm_num; linarith

lemma D_c1 : (1.1847 : ℝ) ≤ √2 * √(√2 * 0.4963) := by
  rw [← Real.sqrt_mul (by norm_num), Real.le_sqrt (by norm_num) (by positivity)]
  have := sqrt2_gt; norm_num; nlinarith

lemma termB23 (n : ℕ) (h : 2000 ≤ n) :
    2 / π * (√(π / 2) * yv n * √(yv n)) * Real.exp (-θv n) *
      (32 * (Mφ + yv n * Mχ) * (∫ x : ℝ, x ^ 4 * Real.exp (-(θv n / 2) * x ^ 2))) ≤
      1430 / (n : ℝ) ^ 2 * Real.exp (-θv n) := by
  obtain ⟨hθ1, hy1, hy0, -⟩ := facts2 n h
  have hN : (0 : ℝ) < n := by have : (2000 : ℝ) ≤ n := by exact_mod_cast h
                              linarith
  have hθ : 0 < θv n := by linarith
  obtain ⟨-, -, h4, -⟩ := gauss_real_vals (θv n / 2) (by positivity)
  rw [h4, show θv n / 2 = 1 / 2 * θv n by ring]
  have hm := nG0_mom n h (1 / 2) (by norm_num) 3 2
  have hM : Mφ + yv n * Mχ ≤ 22.13 := by
    have := Mφ_le; have := Mχ_le
    have hMχ0 : 0 ≤ Mχ := by unfold Mχ Mφ; positivity
    have : yv n * Mχ ≤ 0.028025 * 72.1 := by gcongr
    linarith
  have hM0 : 0 ≤ Mφ + yv n * Mχ := by
    have : 0 ≤ Mφ := by unfold Mφ; positivity
    have : 0 ≤ Mχ := by unfold Mχ Mφ; positivity
    positivity
  have hsp : √(π / 2) ≤ 1.2534 := by
    rw [Real.sqrt_le_left (by norm_num)]; nlinarith [Real.pi_lt_d6]
  have hD := D_half
  have hE := Real.exp_pos (-θv n)
  set G := √(π / 2) * yv n * √(yv n)
  have key : 2 / π * G * Real.exp (-θv n) * (32 * (Mφ + yv n * Mχ) * (3 * √(π / (1 / 2 * θv n)) /
      (2 * (1 / 2 * θv n)) ^ 2)) = (2 / π * 32 * (Mφ + yv n * Mχ) *
      ((n : ℝ) ^ 2 * G * (3 * √(π / (1 / 2 * θv n)) / (2 * (1 / 2 * θv n)) ^ 2))) / (n : ℝ) ^ 2 *
      Real.exp (-θv n) := by field_simp
  rw [key, hm]
  gcongr
  have e2 : √(π / 2) * 3 / (√2 * √(√2 * (1 / 2))) / (2 * (1 / 2)) ^ 2 * θv n ^ 2 / θv n ^ 2 =
      √(π / 2) * 3 / (√2 * √(√2 * (1 / 2))) := by field_simp
  rw [e2]
  calc 2 / π * 32 * (Mφ + yv n * Mχ) * (√(π / 2) * 3 / (√2 * √(√2 * (1 / 2))))
      ≤ 2 / 3.14159 * 32 * 22.13 * (1.2534 * 3 / 1.189) := by
        gcongr
        · linarith [Real.pi_gt_d6]
    _ ≤ 1430 := by norm_num

lemma reduce (n : ℕ) (hn : 0 < n) (G X C θ : ℝ) (h : 2 / π * ((n : ℝ) ^ 2 * G * X) ≤ C) :
    2 / π * G * Real.exp (-θ) * X ≤ C / (n : ℝ) ^ 2 * Real.exp (-θ) := by
  have hN : (0 : ℝ) < n := by exact_mod_cast hn
  have e : 2 / π * G * Real.exp (-θ) * X = (2 / π * ((n : ℝ) ^ 2 * G * X)) / (n : ℝ) ^ 2 *
      Real.exp (-θ) := by field_simp
  rw [e]; gcongr

lemma termB5 (n : ℕ) (h : 2000 ≤ n) :
    2 / π * (√(π / 2) * yv n * √(yv n)) * Real.exp (-θv n) *
      (Real.exp (-(θv n / 64)) * (‖pc n 0‖ * (∫ x : ℝ, x ^ 0 * Real.exp (-(θv n / 4) * x ^ 2)) +
        ‖pc n 2‖ * (∫ x : ℝ, x ^ 2 * Real.exp (-(θv n / 4) * x ^ 2)) +
        √2 * θv n / 24 * ‖pc n 0‖ * (∫ x : ℝ, x ^ 4 * Real.exp (-(θv n / 4) * x ^ 2)))) ≤
      2760 / (n : ℝ) ^ 2 * Real.exp (-θv n) := by
  obtain ⟨hθ1, hy1, hy0, -⟩ := facts2 n h
  have hn : 0 < n := by omega
  have hθ : 0 < θv n := by linarith
  apply reduce n hn
  obtain ⟨h0, h2, h4, -⟩ := gauss_real_vals (θv n / 4) (by positivity)
  rw [h0, h2, h4, show θv n / 4 = 1 / 4 * θv n by ring]
  set θ := θv n
  set G := √(π / 2) * yv n * √(yv n)
  set A := ‖pc n 0‖
  set B := ‖pc n 2‖
  have hA := norm_pc0 n h
  have hB := norm_pc2 n h
  have hA0 : 0 ≤ A := norm_nonneg _
  have hB0 : 0 ≤ B := norm_nonneg _
  have k0 := nG0_mom n h (1 / 4) (by norm_num) 1 0
  have k2 := nG0_mom n h (1 / 4) (by norm_num) 1 1
  have k4 := nG0_mom n h (1 / 4) (by norm_num) 3 2
  simp only [pow_zero, div_one, one_mul, pow_one] at k0 k2
  have e : (n : ℝ) ^ 2 * G * (Real.exp (-(θ / 64)) * (A * √(π / (1 / 4 * θ)) +
      B * (√(π / (1 / 4 * θ)) / (2 * (1 / 4 * θ))) +
      √2 * θ / 24 * A * (3 * √(π / (1 / 4 * θ)) / (2 * (1 / 4 * θ)) ^ 2))) =
      Real.exp (-(θ / 64)) * (A * ((n : ℝ) ^ 2 * G * √(π / (1 / 4 * θ))) +
        B * ((n : ℝ) ^ 2 * G * (√(π / (1 / 4 * θ)) / (2 * (1 / 4 * θ)))) +
        √2 * θ / 24 * A * ((n : ℝ) ^ 2 * G * (3 * √(π / (1 / 4 * θ)) / (2 * (1 / 4 * θ)) ^ 2))) := by
    ring
  rw [e, k0, k2, k4]
  set D := √2 * √(√2 * (1 / 4))
  have hD := D_quarter
  have hsp : √(π / 2) ≤ 1.2534 := by
    rw [Real.sqrt_le_left (by norm_num)]; nlinarith [Real.pi_lt_d6]
  have hs2 := sqrt2_lt
  set X := Real.exp (-(θ / 64))
  have hX0 : 0 < X := Real.exp_pos _
  have hs1 : θ ^ 2 * X ≤ 2217.5 := by
    have := θ2_sup θ hθ (1 / 64) (by norm_num)
    rw [show 1 / 64 * θ = θ / 64 by ring] at this; norm_num at this ⊢; linarith
  have hs0 : θ * X ≤ 23.55 := by
    have := θ1_sup θ hθ (1 / 64) (by norm_num)
    rw [show 1 / 64 * θ = θ / 64 by ring] at this; norm_num at this ⊢; linarith
  have e2 : X * (A * (√(π / 2) * 1 / D * θ ^ 2) + B * (√(π / 2) * 1 / D / (2 * (1 / 4)) * θ ^ 2 / θ) +
      √2 * θ / 24 * A * (√(π / 2) * 3 / D / (2 * (1 / 4)) ^ 2 * θ ^ 2 / θ ^ 2)) =
      √(π / 2) / D * (A * (θ ^ 2 * X) + (2 * B + √2 / 2 * A) * (θ * X)) := by
    field_simp; ring
  rw [e2]
  have hD0 : 0 < D := by positivity
  have h1 : A * (θ ^ 2 * X) + (2 * B + √2 / 2 * A) * (θ * X) ≤
      1.0357 * 2217.5 + (2 * 12.19 + 1.41422 / 2 * 1.0357) * 23.55 := by
    have hθX : 0 ≤ θ * X := by positivity
    have hθ2X : 0 ≤ θ ^ 2 * X := by positivity
    have a1 : A * (θ ^ 2 * X) ≤ 1.0357 * 2217.5 := mul_le_mul hA hs1 hθ2X (by norm_num)
    have a2 : 2 * B + √2 / 2 * A ≤ 2 * 12.19 + 1.41422 / 2 * 1.0357 := by gcongr
    have a3 : (2 * B + √2 / 2 * A) * (θ * X) ≤ (2 * 12.19 + 1.41422 / 2 * 1.0357) * 23.55 :=
      mul_le_mul a2 hs0 hθX (by norm_num)
    linarith
  have hpos : 0 ≤ A * (θ ^ 2 * X) + (2 * B + √2 / 2 * A) * (θ * X) := by positivity
  calc 2 / π * (√(π / 2) / D * (A * (θ ^ 2 * X) + (2 * B + √2 / 2 * A) * (θ * X)))
      ≤ 2 / 3.14159 * (1.2534 / 0.8408 * (1.0357 * 2217.5 +
          (2 * 12.19 + 1.41422 / 2 * 1.0357) * 23.55)) := by
        gcongr
        · linarith [Real.pi_gt_d6]
    _ ≤ 2760 := by norm_num

lemma termB4 (n : ℕ) (h : 2000 ≤ n) :
    2 / π * (√(π / 2) * yv n * √(yv n)) * Real.exp (-θv n) *
      (0.0023 * θv n * ‖pc n 0‖ * (∫ x : ℝ, x ^ 6 * Real.exp (-(0.4963 * θv n) * x ^ 2)) +
        (0.0023 * θv n * ‖pc n 2‖ + 0.0035 * θv n ^ 2 * ‖pc n 0‖) *
          (∫ x : ℝ, x ^ 8 * Real.exp (-(0.4963 * θv n) * x ^ 2)) +
        0.0035 * θv n ^ 2 * ‖pc n 2‖ * (∫ x : ℝ, x ^ 10 * Real.exp (-(0.4963 * θv n) * x ^ 2)) +
        √2 * θv n / 24 * ‖pc n 2‖ * (∫ x : ℝ, x ^ 6 * Real.exp (-(θv n / 2) * x ^ 2))) ≤
      20 / (n : ℝ) ^ 2 * Real.exp (-θv n) := by
  obtain ⟨hθ1, hy1, hy0, -⟩ := facts2 n h
  have hn : 0 < n := by omega
  have hθ : 0 < θv n := by linarith
  apply reduce n hn
  obtain ⟨-, -, -, g6, g8, g10⟩ := gauss_real_vals (0.4963 * θv n) (by positivity)
  obtain ⟨-, -, -, g6', -, -⟩ := gauss_real_vals (θv n / 2) (by positivity)
  rw [g6, g8, g10, g6', show θv n / 2 = 1 / 2 * θv n by ring]
  set θ := θv n
  set G := √(π / 2) * yv n * √(yv n)
  set A := ‖pc n 0‖
  set B := ‖pc n 2‖
  have hA := norm_pc0 n h
  have hB := norm_pc2 n h
  have hA0 : 0 ≤ A := norm_nonneg _
  have hB0 : 0 ≤ B := norm_nonneg _
  have k6 := nG0_mom n h 0.4963 (by norm_num) 15 3
  have k8 := nG0_mom n h 0.4963 (by norm_num) 105 4
  have k10 := nG0_mom n h 0.4963 (by norm_num) 945 5
  have k6' := nG0_mom n h (1 / 2) (by norm_num) 15 3
  set M6 := 15 * √(π / (0.4963 * θ)) / (2 * (0.4963 * θ)) ^ 3
  set M8 := 105 * √(π / (0.4963 * θ)) / (2 * (0.4963 * θ)) ^ 4
  set M10 := 945 * √(π / (0.4963 * θ)) / (2 * (0.4963 * θ)) ^ 5
  set M6' := 15 * √(π / (1 / 2 * θ)) / (2 * (1 / 2 * θ)) ^ 3
  have e : (n : ℝ) ^ 2 * G * (0.0023 * θ * A * M6 + (0.0023 * θ * B + 0.0035 * θ ^ 2 * A) * M8 +
      0.0035 * θ ^ 2 * B * M10 + √2 * θ / 24 * B * M6') =
      0.0023 * θ * A * ((n : ℝ) ^ 2 * G * M6) + (0.0023 * θ * B + 0.0035 * θ ^ 2 * A) *
        ((n : ℝ) ^ 2 * G * M8) + 0.0035 * θ ^ 2 * B * ((n : ℝ) ^ 2 * G * M10) +
        √2 * θ / 24 * B * ((n : ℝ) ^ 2 * G * M6') := by ring
  rw [e, k6, k8, k10, k6']
  set D5 := √2 * √(√2 * 0.4963)
  set D2 := √2 * √(√2 * (1 / 2))
  have hD5 := D_c1
  have hD2 := D_half
  have hsp : √(π / 2) ≤ 1.2534 := by
    rw [Real.sqrt_le_left (by norm_num)]; nlinarith [Real.pi_lt_d6]
  have hsp0 : 0 ≤ √(π / 2) := Real.sqrt_nonneg _
  have hs2 := sqrt2_lt
  have ht : 1 / θ ≤ 0.01262 := by rw [div_le_iff₀ hθ]; linarith
  have ht0 : 0 < 1 / θ := by positivity
  have hD50 : 0 < D5 := by positivity
  have hD20 : 0 < D2 := by positivity
  have e2 : 0.0023 * θ * A * (√(π / 2) * 15 / D5 / (2 * 0.4963) ^ 3 * θ ^ 2 / θ ^ 3) +
      (0.0023 * θ * B + 0.0035 * θ ^ 2 * A) * (√(π / 2) * 105 / D5 / (2 * 0.4963) ^ 4 * θ ^ 2 / θ ^ 4) +
      0.0035 * θ ^ 2 * B * (√(π / 2) * 945 / D5 / (2 * 0.4963) ^ 5 * θ ^ 2 / θ ^ 5) +
      √2 * θ / 24 * B * (√(π / 2) * 15 / D2 / (2 * (1 / 2)) ^ 3 * θ ^ 2 / θ ^ 3) =
      √(π / 2) / D5 * (0.0023 * 15 / (2 * 0.4963) ^ 3 * A +
        105 / (2 * 0.4963) ^ 4 * (0.0023 * B * (1 / θ) + 0.0035 * A) +
        0.0035 * 945 / (2 * 0.4963) ^ 5 * B * (1 / θ)) + √(π / 2) / D2 * (√2 * 15 / 24 * B) := by
    field_simp
  rw [e2]
  have b1 : 0.0023 * 15 / (2 * 0.4963) ^ 3 * A + 105 / (2 * 0.4963) ^ 4 * (0.0023 * B * (1 / θ) +
      0.0035 * A) + 0.0035 * 945 / (2 * 0.4963) ^ 5 * B * (1 / θ) ≤
      0.0023 * 15 / (2 * 0.4963) ^ 3 * 1.0357 + 105 / (2 * 0.4963) ^ 4 * (0.0023 * 12.19 * 0.01262 +
      0.0035 * 1.0357) + 0.0035 * 945 / (2 * 0.4963) ^ 5 * 12.19 * 0.01262 := by
    gcongr
  have b2 : √2 * 15 / 24 * B ≤ 1.41422 * 15 / 24 * 12.19 := by gcongr
  have hp1 : 0 ≤ 0.0023 * 15 / (2 * 0.4963) ^ 3 * A + 105 / (2 * 0.4963) ^ 4 * (0.0023 * B * (1 / θ) +
      0.0035 * A) + 0.0035 * 945 / (2 * 0.4963) ^ 5 * B * (1 / θ) := by positivity
  have hp2 : 0 ≤ √2 * 15 / 24 * B := by positivity
  calc 2 / π * (√(π / 2) / D5 * (0.0023 * 15 / (2 * 0.4963) ^ 3 * A +
        105 / (2 * 0.4963) ^ 4 * (0.0023 * B * (1 / θ) + 0.0035 * A) +
        0.0035 * 945 / (2 * 0.4963) ^ 5 * B * (1 / θ)) + √(π / 2) / D2 * (√2 * 15 / 24 * B))
      ≤ 2 / 3.14159 * (1.2534 / 1.1847 * (0.0023 * 15 / (2 * 0.4963) ^ 3 * 1.0357 +
        105 / (2 * 0.4963) ^ 4 * (0.0023 * 12.19 * 0.01262 + 0.0035 * 1.0357) +
        0.0035 * 945 / (2 * 0.4963) ^ 5 * 12.19 * 0.01262) +
        1.2534 / 1.189 * (1.41422 * 15 / 24 * 12.19)) := by
        gcongr
        · linarith [Real.pi_gt_d6]
    _ ≤ 20 := by norm_num

lemma regII_hyp (n : ℕ) (h : 2000 ≤ n) :
    e1v n < 1 ∧ uv n (-(1 / 4)) ^ 2 ≤ 2 ∧ 1 < Sl n (1 / 4) := by
  obtain ⟨hθ1, hy1, hy0, -⟩ := facts2 n h
  have hn : 0 < n := by omega
  have hN : (0 : ℝ) < n := by exact_mod_cast hn
  have e14 := exp_quarter_lt
  have e14' := exp_quarter_gt
  have em14 := exp_neg_quarter_lt
  have hu1 : uv n (-(1 / 4)) ≤ 0.051 := by
    rw [uv_neg]
    have h1 : √2 * yv n ≤ 1.41422 * 0.028025 := by have := sqrt2_lt; gcongr
    have h2 : √2 * yv n * Real.exp (1 / 4) ≤ (1.41422 * 0.028025) * 1.2841 := by gcongr
    linarith
  have hu0 := uv_nonneg n (-(1 / 4)) hN
  have hu2 : uv n (-(1 / 4)) ^ 2 ≤ 0.051 ^ 2 := pow_le_pow_left₀ hu0 hu1 2
  refine ⟨?_, by linarith, ?_⟩
  · unfold e1v
    have : yv n * Real.exp (1 / 4) ≤ 0.028025 * 1.2841 := by gcongr
    linarith
  · have hS := Sl_ge n hn (1 / 4)
    have h3 : Real.exp (1 / 4) / 2 * (uv n (-(1 / 4)) ^ 2 / 2) ≤ 1.2841 / 2 * (0.051 ^ 2 / 2) := by
      gcongr
    have h4 : 0.25 ≤ Real.exp (1 / 4) / 2 * (1 - uv n (-(1 / 4)) ^ 2 / 2) - Real.exp (-(1 / 4)) / 2 := by
      have : Real.exp (1 / 4) / 2 * (1 - uv n (-(1 / 4)) ^ 2 / 2) =
          Real.exp (1 / 4) / 2 - Real.exp (1 / 4) / 2 * (uv n (-(1 / 4)) ^ 2 / 2) := by ring
      rw [this]; norm_num at h3 ⊢; linarith
    have := mul_le_mul_of_nonneg_left h4 (by linarith : (0 : ℝ) ≤ θv n)
    linarith

lemma regII (k : ℕ) (h1 : 2000 ≤ k + 1) :
    |Mmain k + X ((k + 1 : ℕ) : ℝ)| ≤
      71000 / (((k + 1 : ℕ) : ℝ)) ^ 2 * Real.exp (-√(π * ((k + 1 : ℕ) : ℝ))) := by
  have hk : 1 ≤ k := by omega
  set n := k + 1 with hn
  have hn0 : 0 < n := by omega
  obtain ⟨he1, hu, hS⟩ := regII_hyp n h1
  have hf := (J_eq k hk).1
  have hcore := regII_core n hn0 he1 hu hS hf
  have hIm : Mmain k + X n = ((2 / π : ℂ) * (∫ σ : ℝ, fσ n σ) - Tv n).im := by
    rw [Complex.sub_im, T_im, Mmain_eq k hk]
    have : ((2 / π : ℂ) * ∫ σ : ℝ, fσ n σ).im = 2 / π * (∫ σ : ℝ, fσ n σ).im := by
      rw [show (2 / π : ℂ) = ((2 / π : ℝ) : ℂ) by push_cast; rfl, Complex.im_ofReal_mul]
    rw [this]; ring
  rw [hIm]
  refine (Complex.abs_im_le_norm _).trans (hcore.trans ?_)
  have tLR := termLR n h1
  have tD := termD n h1
  have t23 := termB23 n h1
  have t4 := termB4 n h1
  have t5 := termB5 n h1
  have t7 := term7 n h1
  have hE := Real.exp_pos (-θv n)
  have hN : (0 : ℝ) < n := by exact_mod_cast hn0
  show _ ≤ 71000 / (n : ℝ) ^ 2 * Real.exp (-θv n)
  rw [mul_add, mul_add (2 / π * _ * _), mul_add (2 / π * _ * _)]
  calc _ ≤ 1000 / (n : ℝ) ^ 2 * Real.exp (-θv n) + 200 / (n : ℝ) ^ 2 * Real.exp (-θv n) +
        (1430 / (n : ℝ) ^ 2 * Real.exp (-θv n) + 20 / (n : ℝ) ^ 2 * Real.exp (-θv n) +
          2760 / (n : ℝ) ^ 2 * Real.exp (-θv n)) + 300 / (n : ℝ) ^ 2 * Real.exp (-θv n) := by
        gcongr
    _ = 5710 / (n : ℝ) ^ 2 * Real.exp (-θv n) := by ring
    _ ≤ 71000 / (n : ℝ) ^ 2 * Real.exp (-θv n) := by gcongr; norm_num

end Q07.PF

/-- g5:prop:main. -/
theorem Q07.PF.main : Q07.stmt_6 := by
  intro k hk
  rcases le_or_gt (k + 1) 2000 with h | h
  · exact Q07.PF.regI k hk h
  · exact Q07.PF.regII k (by omega)

/-! ### Merged from work/Q07/PartG.lean -/
namespace Q07.PG

open Real Q05

lemma rpow_three_halves (n : ℝ) (hn : 0 < n) : n ^ (-(3 / 2) : ℝ) = 1 / (n * √n) := by
  rw [Real.rpow_neg hn.le, show (3 / 2 : ℝ) = 1 + 1 / 2 by norm_num, Real.rpow_add hn,
    Real.rpow_one, ← Real.sqrt_eq_rpow, one_div]

/-- `X n = (√π/n) e^{-θ} (cos θ + Re(e^{iθ} ϖ₁)/√n)`. -/
lemma X_expand (n : ℝ) :
    X n = √π / n * Real.exp (-√(π * n)) *
      (Real.cos (√(π * n)) + (Complex.exp (Complex.I * (√(π * n) : ℝ)) * varpi1).re / √n) := by
  unfold X
  congr 1
  rw [mul_add, mul_one, Complex.add_re, ← mul_div_assoc, Complex.div_ofReal_re, mul_comm Complex.I,
    Complex.exp_ofReal_mul_I_re]

lemma X_diff (n : ℝ) (hn : 0 < n) :
    |√π / n * Real.exp (-√(π * n)) * Real.cos (√(π * n)) - X n| ≤
      1.494 * (1 / (n * √n)) * Real.exp (-√(π * n)) := by
  rw [X_expand]
  set z := (Complex.exp (Complex.I * (√(π * n) : ℝ)) * varpi1).re
  have hz : |z| ≤ ‖varpi1‖ := by
    refine (Complex.abs_re_le_norm _).trans ?_
    rw [norm_mul, mul_comm Complex.I, Complex.norm_exp_ofReal_mul_I, one_mul]
  have hs : 0 < √n := Real.sqrt_pos.2 hn
  have e : √π / n * Real.exp (-√(π * n)) * Real.cos (√(π * n)) -
      √π / n * Real.exp (-√(π * n)) * (Real.cos (√(π * n)) + z / √n) =
      -(√π * z * (1 / (n * √n)) * Real.exp (-√(π * n))) := by
    field_simp; ring
  rw [e, abs_neg, abs_mul, abs_mul, abs_mul, abs_of_pos (Real.sqrt_pos.2 Real.pi_pos),
    abs_of_pos (by positivity : (0 : ℝ) < 1 / (n * √n)), abs_of_pos (Real.exp_pos _)]
  have h3 := Q05.proof_3
  unfold Q05.stmt_3 at h3
  have : √π * |z| ≤ 1.494 :=
    (mul_le_mul_of_nonneg_left hz (Real.sqrt_nonneg _)).trans h3
  have hp : 0 < 1 / (n * √n) * Real.exp (-√(π * n)) := by positivity
  nlinarith

lemma exp_two_sqrt (n : ℝ) (hn : 0 < n) :
    Real.exp (-2 * √n) ≤ 114.1 * (1 / (n * √n)) * Real.exp (-√(π * n)) := by
  obtain ⟨h1, -, h3⟩ := Q05.proof_1
  set u := √n
  have hu : 0 < u := Real.sqrt_pos.2 hn
  have hb := (h1 u hu).trans h3.le
  have hnu : n * u = u ^ 3 := by
    have : u * u = n := Real.mul_self_sqrt hn.le
    rw [← this]; ring
  have hpi : √(π * n) = √π * u := Real.sqrt_mul Real.pi_pos.le n
  rw [hnu, hpi]
  have he : Real.exp (-2 * u) = Real.exp (-(2 - √π) * u) * Real.exp (-(√π * u)) := by
    rw [← Real.exp_add]; ring_nf
  rw [he]
  have hu3 : 0 < u ^ 3 := by positivity
  rw [show 114.1 * (1 / u ^ 3) * Real.exp (-(√π * u)) =
    (114.1 / u ^ 3) * Real.exp (-(√π * u)) by ring]
  apply mul_le_mul_of_nonneg_right _ (Real.exp_pos _).le
  rw [le_div_iff₀ hu3]; linarith

lemma main_le (n : ℝ) (hn : 100 ≤ n) :
    71000 / n ^ 2 * Real.exp (-√(π * n)) ≤ 7100 * (1 / (n * √n)) * Real.exp (-√(π * n)) := by
  apply mul_le_mul_of_nonneg_right _ (Real.exp_pos _).le
  have hn0 : 0 < n := by linarith
  have hs : 10 ≤ √n := by
    rw [show (10 : ℝ) = √100 by rw [show (100 : ℝ) = 10 ^ 2 by norm_num, Real.sqrt_sq (by norm_num)]]
    exact Real.sqrt_le_sqrt hn
  have hsq : √n * √n = n := Real.mul_self_sqrt hn0.le
  rw [div_le_iff₀ (by positivity)]
  have : 7100 * (1 / (n * √n)) * n ^ 2 = 7100 * √n := by
    field_simp; nlinarith
  rw [this]; linarith

/-- The first conjunct, for real `n ≥ 100`. -/
lemma decay_real (n m mm : ℝ) (hn : 100 ≤ n)
    (h5 : |m - mm| ≤ 10600000 * Real.exp (-2 * √n))
    (h6 : |mm + X n| ≤ 71000 / n ^ 2 * Real.exp (-√(π * n))) :
    |m + √π / n * Real.exp (-√(π * n)) * Real.cos (√(π * n))| ≤
      1.21e9 * (1 / (n * √n)) * Real.exp (-√(π * n)) := by
  have hn0 : 0 < n := by linarith
  have hA := exp_two_sqrt n hn0
  have hB := main_le n hn
  have hC := X_diff n hn0
  have e : m + √π / n * Real.exp (-√(π * n)) * Real.cos (√(π * n)) =
      (m - mm) + (mm + X n) + (√π / n * Real.exp (-√(π * n)) * Real.cos (√(π * n)) - X n) := by
    ring
  rw [e]
  have hP : 0 < 1 / (n * √n) * Real.exp (-√(π * n)) := by
    have := Real.sqrt_pos.2 hn0; positivity
  have t := abs_add_three (m - mm) (mm + X n)
    (√π / n * Real.exp (-√(π * n)) * Real.cos (√(π * n)) - X n)
  nlinarith

/-- `δ = √(π(x+1)) - √(πx) ≤ √π/(2√x)`. -/
lemma delta_le (x : ℝ) (hx : 0 < x) :
    0 ≤ √(π * (x + 1)) - √(π * x) ∧ √(π * (x + 1)) - √(π * x) ≤ √π / (2 * √x) := by
  have hp := Real.pi_pos
  have h0 : 0 < √(π * x) := Real.sqrt_pos.2 (by positivity)
  have hle : √(π * x) ≤ √(π * (x + 1)) := Real.sqrt_le_sqrt (by nlinarith)
  refine ⟨by linarith, ?_⟩
  have s1 : √(π * (x + 1)) ^ 2 = π * (x + 1) := Real.sq_sqrt (by positivity)
  have s0 : √(π * x) ^ 2 = π * x := Real.sq_sqrt (by positivity)
  have hprod : (√(π * (x + 1)) - √(π * x)) * (√(π * (x + 1)) + √(π * x)) = π := by
    nlinarith
  have hpx : √(π * x) = √π * √x := Real.sqrt_mul hp.le x
  have hsx : 0 < √x := Real.sqrt_pos.2 hx
  have hsp : 0 < √π := Real.sqrt_pos.2 hp
  have hspp : √π * √π = π := Real.mul_self_sqrt hp.le
  rw [le_div_iff₀ (by positivity)]
  generalize hd : √(π * (x + 1)) - √(π * x) = d at hprod ⊢
  have hd0 : 0 ≤ d := by rw [← hd]; linarith
  have key : d * (2 * √(π * x)) ≤ d * (√(π * (x + 1)) + √(π * x)) :=
    mul_le_mul_of_nonneg_left (by linarith) hd0
  rw [hprod, hpx] at key
  have k2 : (d * (2 * √x)) * √π ≤ √π * √π := by rw [hspp]; linarith
  exact le_of_mul_le_mul_right k2 hsp

/-- `|e^{-t₁} cos t₁ - e^{-t₀} cos t₀| ≤ 2 (t₁ - t₀) e^{-t₀}` for `t₀ ≤ t₁`. -/
lemma g_diff (t0 t1 : ℝ) (h : t0 ≤ t1) :
    |Real.exp (-t1) * Real.cos t1 - Real.exp (-t0) * Real.cos t0| ≤
      2 * (t1 - t0) * Real.exp (-t0) := by
  have e : Real.exp (-t1) * Real.cos t1 - Real.exp (-t0) * Real.cos t0 =
      Real.exp (-t1) * (Real.cos t1 - Real.cos t0) + (Real.exp (-t1) - Real.exp (-t0)) * Real.cos t0 := by
    ring
  rw [e]
  have hc := Real.abs_cos_sub_cos_le t1 t0
  rw [abs_of_nonneg (by linarith : (0 : ℝ) ≤ t1 - t0)] at hc
  have he1 : Real.exp (-t1) ≤ Real.exp (-t0) := Real.exp_le_exp.2 (by linarith)
  have hexp : Real.exp (-t1) = Real.exp (-t0) * Real.exp (-(t1 - t0)) := by
    rw [← Real.exp_add]; ring_nf
  have hd := Real.add_one_le_exp (-(t1 - t0))
  have he2 : |Real.exp (-t1) - Real.exp (-t0)| ≤ (t1 - t0) * Real.exp (-t0) := by
    rw [abs_of_nonpos (by linarith)]
    rw [hexp]; nlinarith [Real.exp_pos (-t0)]
  have hcos : |Real.cos t0| ≤ 1 := Real.abs_cos_le_one t0
  calc |Real.exp (-t1) * (Real.cos t1 - Real.cos t0) + (Real.exp (-t1) - Real.exp (-t0)) * Real.cos t0|
      ≤ |Real.exp (-t1) * (Real.cos t1 - Real.cos t0)| +
          |(Real.exp (-t1) - Real.exp (-t0)) * Real.cos t0| := abs_add_le _ _
    _ = Real.exp (-t1) * |Real.cos t1 - Real.cos t0| +
          |Real.exp (-t1) - Real.exp (-t0)| * |Real.cos t0| := by
        rw [abs_mul, abs_mul, abs_of_pos (Real.exp_pos _)]
    _ ≤ Real.exp (-t0) * (t1 - t0) + (t1 - t0) * Real.exp (-t0) * 1 := by
        gcongr
    _ = 2 * (t1 - t0) * Real.exp (-t0) := by ring

/-- The shift between `N = x+1` and `x`. -/
lemma shift_diff (x : ℝ) (hx : 99 ≤ x) :
    |√π / (x + 1) * Real.exp (-√(π * (x + 1))) * Real.cos (√(π * (x + 1))) -
        √π / x * Real.exp (-√(π * x)) * Real.cos (√(π * x))| ≤
      √π * 1.9 * (1 / (x * √x)) * Real.exp (-√(π * x)) := by
  have hx0 : 0 < x := by linarith
  set t0 := √(π * x)
  set t1 := √(π * (x + 1))
  obtain ⟨hd0, hd1⟩ := delta_le x hx0
  have hg := g_diff t0 t1 (by linarith)
  have hsx : 0 < √x := Real.sqrt_pos.2 hx0
  have hsp : 0 < √π := Real.sqrt_pos.2 Real.pi_pos
  have hsq : √x * √x = x := Real.mul_self_sqrt hx0.le
  have h99 : 9.9 ≤ √x := by
    rw [show (9.9 : ℝ) = √(9.9 ^ 2) by rw [Real.sqrt_sq (by norm_num)]]
    exact Real.sqrt_le_sqrt (by norm_num; linarith)
  have e : √π / (x + 1) * Real.exp (-t1) * Real.cos t1 - √π / x * Real.exp (-t0) * Real.cos t0 =
      √π * ((1 / (x + 1) - 1 / x) * (Real.exp (-t1) * Real.cos t1) +
        1 / x * (Real.exp (-t1) * Real.cos t1 - Real.exp (-t0) * Real.cos t0)) := by
    field_simp; ring
  rw [e, abs_mul, abs_of_pos hsp, mul_assoc (√π) 1.9, mul_assoc (√π)]
  apply mul_le_mul_of_nonneg_left _ hsp.le
  have h1 : |1 / (x + 1) - 1 / x| ≤ 1 / (x * x) := by
    rw [abs_sub_comm, abs_of_nonneg (by
      rw [sub_nonneg]; exact one_div_le_one_div_of_le hx0 (by linarith))]
    rw [div_sub_div _ _ (by linarith) (by linarith)]
    apply div_le_div₀ (by norm_num) (by linarith) (by positivity) (by nlinarith)
  have h2 : |Real.exp (-t1) * Real.cos t1| ≤ Real.exp (-t0) := by
    rw [abs_mul, abs_of_pos (Real.exp_pos _)]
    calc Real.exp (-t1) * |Real.cos t1| ≤ Real.exp (-t1) * 1 :=
          mul_le_mul_of_nonneg_left (Real.abs_cos_le_one _) (Real.exp_pos _).le
      _ ≤ Real.exp (-t0) := by rw [mul_one]; exact Real.exp_le_exp.2 (by linarith)
  have hE := Real.exp_pos (-t0)
  have hgd : 2 * (t1 - t0) ≤ √π / √x := by
    have : √π / (2 * √x) * 2 = √π / √x := by field_simp
    linarith
  calc |(1 / (x + 1) - 1 / x) * (Real.exp (-t1) * Real.cos t1) +
          1 / x * (Real.exp (-t1) * Real.cos t1 - Real.exp (-t0) * Real.cos t0)|
      ≤ |(1 / (x + 1) - 1 / x) * (Real.exp (-t1) * Real.cos t1)| +
          |1 / x * (Real.exp (-t1) * Real.cos t1 - Real.exp (-t0) * Real.cos t0)| :=
        abs_add_le _ _
    _ = |1 / (x + 1) - 1 / x| * |Real.exp (-t1) * Real.cos t1| +
          1 / x * |Real.exp (-t1) * Real.cos t1 - Real.exp (-t0) * Real.cos t0| := by
        rw [abs_mul (1 / (x + 1) - 1 / x), abs_mul (1 / x),
          abs_of_pos (by positivity : (0 : ℝ) < 1 / x)]
    _ ≤ 1 / (x * x) * Real.exp (-t0) + 1 / x * ((√π / √x) * Real.exp (-t0)) := by
        gcongr
        exact hg.trans (mul_le_mul_of_nonneg_right hgd hE.le)
    _ = (1 / √x + √π) * (1 / (x * √x)) * Real.exp (-t0) := by
        field_simp; ring_nf; rw [show √x ^ 2 = x by rw [sq]; exact hsq]
    _ ≤ 1.9 * (1 / (x * √x)) * Real.exp (-t0) := by
        apply mul_le_mul_of_nonneg_right _ (by positivity)
        apply mul_le_mul_of_nonneg_right _ (by positivity)
        have : 1 / √x ≤ 1 / 9.9 := one_div_le_one_div_of_le (by norm_num) h99
        linarith [Q05.sqrt_pi_lt]

theorem main : Q07.stmt_5 → Q07.stmt_6 → Q07.stmt_7 := by
  intro h5 h6
  have first : ∀ k : ℕ, 100 ≤ k + 1 →
      |M k + √π / ((k + 1 : ℕ) : ℝ) * Real.exp (-√(π * ((k + 1 : ℕ) : ℝ))) *
          Real.cos (√(π * ((k + 1 : ℕ) : ℝ)))| ≤
        1.21e9 * (((k + 1 : ℕ) : ℝ)) ^ (-(3 / 2) : ℝ) * Real.exp (-√(π * ((k + 1 : ℕ) : ℝ))) := by
    intro k hk
    have hn : (100 : ℝ) ≤ ((k + 1 : ℕ) : ℝ) := by exact_mod_cast hk
    rw [rpow_three_halves _ (by linarith)]
    exact decay_real _ _ _ hn (h5 k hk).2 (h6 k hk)
  refine ⟨first, fun k hk => ?_⟩
  have hk1 : 100 ≤ k + 1 := by omega
  have h := first k hk1
  have hx : (99 : ℝ) ≤ k := by exact_mod_cast hk
  have hx0 : (0 : ℝ) < k := by linarith
  push_cast at h
  have hs := shift_diff (k : ℝ) hx
  rw [rpow_three_halves _ (by linarith)] at h
  rw [rpow_three_halves _ hx0]
  -- monotonicity of the bound in N
  have hsx : 0 < √(k : ℝ) := Real.sqrt_pos.2 hx0
  have hmono1 : 1 / (((k : ℝ) + 1) * √((k : ℝ) + 1)) ≤ 1 / ((k : ℝ) * √(k : ℝ)) := by
    apply one_div_le_one_div_of_le (by positivity)
    apply mul_le_mul (by linarith) (Real.sqrt_le_sqrt (by linarith)) hsx.le (by linarith)
  have hmono2 : Real.exp (-√(π * ((k : ℝ) + 1))) ≤ Real.exp (-√(π * (k : ℝ))) :=
    Real.exp_le_exp.2 (neg_le_neg (Real.sqrt_le_sqrt (by nlinarith [Real.pi_pos])))
  have hb : 1.21e9 * (1 / (((k : ℝ) + 1) * √((k : ℝ) + 1))) * Real.exp (-√(π * ((k : ℝ) + 1))) ≤
      1.21e9 * (1 / ((k : ℝ) * √(k : ℝ))) * Real.exp (-√(π * (k : ℝ))) := by
    gcongr
  have e : M k + √π / (k : ℝ) * Real.exp (-√(π * k)) * Real.cos (√(π * k)) =
      (M k + √π / ((k : ℝ) + 1) * Real.exp (-√(π * ((k : ℝ) + 1))) *
          Real.cos (√(π * ((k : ℝ) + 1)))) -
        (√π / ((k : ℝ) + 1) * Real.exp (-√(π * ((k : ℝ) + 1))) * Real.cos (√(π * ((k : ℝ) + 1))) -
          √π / (k : ℝ) * Real.exp (-√(π * (k : ℝ))) * Real.cos (√(π * (k : ℝ)))) := by ring
  rw [e]
  refine (abs_sub _ _).trans ?_
  have : (1.21e9 + √π * 1.9) * (1 / ((k : ℝ) * √(k : ℝ))) * Real.exp (-√(π * (k : ℝ))) =
      1.21e9 * (1 / ((k : ℝ) * √(k : ℝ))) * Real.exp (-√(π * (k : ℝ))) +
        √π * 1.9 * (1 / ((k : ℝ) * √(k : ℝ))) * Real.exp (-√(π * (k : ℝ))) := by ring
  rw [this]
  linarith

end Q07.PG

theorem Q07.proof_1 : Q07.stmt_1 := Q07.PA.main
theorem Q07.proof_2 : Q07.stmt_2 := Q07.PB.main Q07.proof_1
theorem Q07.proof_3 : Q07.stmt_3 := Q07.PC.main
theorem Q07.proof_4 : Q07.stmt_4 := Q07.PD.main Q07.proof_2 Q07.proof_3
theorem Q07.proof_5 : Q07.stmt_5 := Q07.PE.main Q07.proof_4
theorem Q07.proof_6 : Q07.stmt_6 := Q07.PF.main
theorem Q07.proof_7 : Q07.stmt_7 := Q07.PG.main Q07.proof_5 Q07.proof_6
