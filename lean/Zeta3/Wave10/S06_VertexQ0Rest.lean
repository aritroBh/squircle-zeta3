import Mathlib
import Zeta3.Wave10.S02_VertexQ0
/-!
# S06 -- continuation of S02: the vertex constant `Q₀`, the simplex and `J` integrals, `P₁, P₂ₐ, P₂c`, closed form and decimals of `s₂`

Target id: S06, the continuation of `Zeta3.Wave10.S02_VertexQ0` (planner split, 2026-11-26).
Paper label: `thm:S` (paper/ieee/main.tex lines 1921–1943), the closed form of `Q₀`; the sentence
"Then $s_2$ is \eqref{eq:s2} in Table~\ref{tab:closed}" (main.tex 1944) with `eq:s2`
(main.tex 949–950) and the decimals of `s₂` (main.tex 1944–1947).
Statements `stmt_1, 2, 4, 6, 7, 8, 9, 10, 11` of S02 keep their numbers; `stmt_3` and `stmt_5` stay in S02,
whose definitions (`S02.J`, `S02.P1`, …) are used here.
-/
-- BEGIN STATEMENTS
namespace S06

open Real MeasureTheory Set

/-- Lemma `g3:lem:Q0simplex`: `D ≥ 0` (on `[0,∞)³`), the integrand `D(σ)/(σ₁σ₂σ₃)` is
integrable on `Δ`, and `Q₀ = -∫_Δ D(σ)/(σ₁σ₂σ₃) dσ`. -/
def stmt_1 : Prop :=
  (∀ x₁ x₂ x₃ : ℝ, 0 ≤ x₁ → 0 ≤ x₂ → 0 ≤ x₃ → 0 ≤ S02.Dfun x₁ x₂ x₃) ∧
    IntegrableOn S02.simplexIntegrand Y05.Δ ∧
    Z04.Q 0 = -∫ t in Y05.Δ, S02.simplexIntegrand t

/-- Lemma `g3:lem:step1`: for `q ∈ [1/√2, 1]` the integral `J(q)` converges, `s ↦ J(q(s))`
is integrable on `[0,1]`, and `Q₀ = 6 ∫₀¹ J(q(s)) ds`. -/
def stmt_2 : Prop :=
  (∀ q : ℝ, 1 / Real.sqrt 2 ≤ q → q ≤ 1 → IntegrableOn (S02.Jintegrand q) (Ioi 0)) ∧
    IntervalIntegrable (fun s => S02.J (S02.qs s)) volume 0 1 ∧
    Z04.Q 0 = 6 * ∫ s in (0 : ℝ)..1, S02.J (S02.qs s)

/-- sec_g3_S.tex 168, 208–209: the three integrands are integrable on `[0,1]` and
`Q₀ = 6(P₁ - P_{2a} - P_{2c})`. -/
def stmt_4 : Prop :=
  IntervalIntegrable (fun s => S02.f1 (S02.qs s)) volume 0 1 ∧
    IntervalIntegrable (fun s => Real.arsinh (S02.qs s) / Real.sqrt (1 + S02.qs s ^ 2)) volume 0 1 ∧
    IntervalIntegrable (fun s => Real.arsinh (1 / S02.qs s) / Real.sqrt (1 + S02.qs s ^ 2)) volume 0 1 ∧
    Z04.Q 0 = 6 * (S02.P1 - S02.P2a - S02.P2c)

/-- Lemma `g3:lem:P`, first line: the integrand of `P₁` is integrable and
`P₁ = 3√2Li₂(b) - λ²/2 - (√2/2)λL + (1+3√2/2)L² - (1/48+3√2/16)π²`. -/
def stmt_6 : Prop :=
  IntervalIntegrable (fun s => S02.f1 (S02.qs s)) volume 0 1 ∧
    S02.P1 = 3 * Real.sqrt 2 * D00.Li2 D02.b - D02.lam ^ 2 / 2
      - Real.sqrt 2 / 2 * D02.lam * D02.L + (1 + 3 * Real.sqrt 2 / 2) * D02.L ^ 2
      - (1 / 48 + 3 * Real.sqrt 2 / 16) * π ^ 2

/-- Lemma `g3:lem:P`, second line: the integrand of `P_{2a}` is integrable and
`P_{2a} = √2[3Li₂(b) - ¼Li₂(-b⁴) - λL/2 + L log 3 + L²/2 - 5π²/24]`. -/
def stmt_7 : Prop :=
  IntervalIntegrable (fun s => Real.arsinh (S02.qs s) / Real.sqrt (1 + S02.qs s ^ 2)) volume 0 1 ∧
    S02.P2a = Real.sqrt 2 * (3 * D00.Li2 D02.b - (1 / 4) * D00.Li2 (-D02.b ^ 4)
      - D02.lam * D02.L / 2 + D02.L * Real.log 3 + D02.L ^ 2 / 2 - 5 * π ^ 2 / 24)

/-- Lemma `g3:lem:P`, third line: the integrand of `P_{2c}` is integrable and
`P_{2c} = √2[Li₂(b) - ⅛Li₂(-b⁴) + (L log 3)/2 + π²/96 - πα_c/4]`. -/
def stmt_8 : Prop :=
  IntervalIntegrable (fun s => Real.arsinh (1 / S02.qs s) / Real.sqrt (1 + S02.qs s ^ 2)) volume 0 1 ∧
    S02.P2c = Real.sqrt 2 * (D00.Li2 D02.b - (1 / 8) * D00.Li2 (-D02.b ^ 4)
      + D02.L * Real.log 3 / 2 + π ^ 2 / 96 - π * S02.alphac / 4)

/-- `thm:S`, the closed form of the vertex constant: `Q₀ = Q(0)` equals
`-π²/8 + (√2/16)π² + (3√2/2)πα_c - 3λ² - 9√2 L log 3 + (6+6√2)L² - 6√2Li₂(b)
+ (9√2/4)Li₂(-b⁴)`. -/
def stmt_9 : Prop :=
  Z04.Q 0 = S02.Q0closed

/-- main.tex 1944, "Then $s_2$ is \eqref{eq:s2}": the coefficient `s₂` of `p^{-2}` in the
expansion of `S(p)` (R04) equals the closed form (T6). -/
def stmt_10 : Prop :=
  R04.s 2 = S02.s2closed

/-- main.tex 1946, "$s_2=-4.6892155184\ldots$" (truncated decimals):
`-4.6892155185 < s₂ ≤ -4.6892155184`. -/
def stmt_11 : Prop :=
  -4.6892155185 < R04.s 2 ∧ R04.s 2 ≤ -4.6892155184

end S06
-- END STATEMENTS

/-! ### Merged from work/S06/PartA.lean -/
/-!
# S06 Part A: `S06.stmt_1` (Lemma `g3:lem:Q0simplex`)

Route: `h₀(y)³ = ∫_Q Π hIntegrand`, Tonelli in `(g, y)`, the `y`-integral
`∫₀^∞ Π(1 - e^{-xᵢ²/y²}) dy = √π D(x)` (inclusion–exclusion and `Z04.gy`), and the Dirichlet
scaling `g = r σ` (lintegral version of `Y05.dirichlet`). Integrability of `h₀³` is `R04.proof_6`.
-/

namespace S06.PA

open Real MeasureTheory Set Filter Topology
open scoped ENNReal

/-- `E_c(y) = 1 - e^{-c/y²}`. -/
noncomputable def E (c y : ℝ) : ℝ := 1 - Real.exp (-(c / y ^ 2))

lemma E_int {c : ℝ} (hc : 0 ≤ c) :
    IntegrableOn (E c) (Ioi 0) ∧ ∫ y in Ioi (0 : ℝ), E c y = Real.sqrt π * Real.sqrt c :=
  Z04.gy hc

lemma prod_expand (c₁ c₂ c₃ y : ℝ) :
    E c₁ y * E c₂ y * E c₃ y =
      E c₁ y + E c₂ y + E c₃ y - E (c₁ + c₂) y - E (c₂ + c₃) y - E (c₁ + c₃) y
        + E (c₁ + c₂ + c₃) y := by
  unfold E
  simp only [add_div, neg_add, Real.exp_add]
  ring

lemma prod_int {c₁ c₂ c₃ : ℝ} (h1 : 0 ≤ c₁) (h2 : 0 ≤ c₂) (h3 : 0 ≤ c₃) :
    IntegrableOn (fun y => E c₁ y * E c₂ y * E c₃ y) (Ioi 0) ∧
      ∫ y in Ioi (0 : ℝ), E c₁ y * E c₂ y * E c₃ y =
        Real.sqrt π * (Real.sqrt c₁ + Real.sqrt c₂ + Real.sqrt c₃ - Real.sqrt (c₁ + c₂)
          - Real.sqrt (c₂ + c₃) - Real.sqrt (c₁ + c₃) + Real.sqrt (c₁ + c₂ + c₃)) := by
  obtain ⟨i1, v1⟩ := E_int h1
  obtain ⟨i2, v2⟩ := E_int h2
  obtain ⟨i3, v3⟩ := E_int h3
  obtain ⟨i4, v4⟩ := E_int (add_nonneg h1 h2)
  obtain ⟨i5, v5⟩ := E_int (add_nonneg h2 h3)
  obtain ⟨i6, v6⟩ := E_int (add_nonneg h1 h3)
  obtain ⟨i7, v7⟩ := E_int (add_nonneg (add_nonneg h1 h2) h3)
  have e : (fun y => E c₁ y * E c₂ y * E c₃ y) = fun y =>
      E c₁ y + E c₂ y + E c₃ y - E (c₁ + c₂) y - E (c₂ + c₃) y - E (c₁ + c₃) y
        + E (c₁ + c₂ + c₃) y := funext fun y => prod_expand c₁ c₂ c₃ y
  rw [e]
  have j2 : IntegrableOn (fun y => E c₁ y + E c₂ y) (Ioi 0) := i1.add i2
  have j3 : IntegrableOn (fun y => E c₁ y + E c₂ y + E c₃ y) (Ioi 0) := j2.add i3
  have j4 : IntegrableOn (fun y => E c₁ y + E c₂ y + E c₃ y - E (c₁ + c₂) y) (Ioi 0) :=
    j3.sub i4
  have j5 : IntegrableOn (fun y => E c₁ y + E c₂ y + E c₃ y - E (c₁ + c₂) y - E (c₂ + c₃) y)
    (Ioi 0) := j4.sub i5
  have j6 : IntegrableOn (fun y => E c₁ y + E c₂ y + E c₃ y - E (c₁ + c₂) y - E (c₂ + c₃) y
    - E (c₁ + c₃) y) (Ioi 0) := j5.sub i6
  have j7 : IntegrableOn (fun y => E c₁ y + E c₂ y + E c₃ y - E (c₁ + c₂) y - E (c₂ + c₃) y
    - E (c₁ + c₃) y + E (c₁ + c₂ + c₃) y) (Ioi 0) := j6.add i7
  refine ⟨j7, ?_⟩
  rw [integral_add j6 i7, integral_sub j5 i6, integral_sub j4 i5, integral_sub j3 i4,
    integral_add j2 i3, integral_add i1 i2, v1, v2, v3, v4, v5, v6, v7]
  ring

/-- `hIntegrand 0 y x = x⁻¹ e^{-x} E_{x²}(y)` for `x > 0`. -/
lemma hI0 {x : ℝ} (_hx : 0 < x) (y : ℝ) :
    Z04.hIntegrand 0 y x = x⁻¹ * Real.exp (-x) * E (x ^ 2) y := by
  unfold Z04.hIntegrand E
  rw [show (2 : ℝ) - 2 * 0 = (2 : ℕ) by norm_num, Real.rpow_natCast,
    show (0 : ℝ) - 1 = -1 by norm_num, Real.rpow_neg_one]
  ring

/-- The Laplace-side integrand `e^{-(g₁+g₂+g₃)} D(g)/(g₁g₂g₃)`. -/
noncomputable def K (g : ℝ × ℝ × ℝ) : ℝ :=
  Real.exp (-(g.1 + g.2.1 + g.2.2)) * (S02.Dfun g.1 g.2.1 g.2.2 / (g.1 * g.2.1 * g.2.2))

/-- The `(g, y)` integrand `Π hIntegrand 0 y gᵢ`. -/
noncomputable def Φ (z : (ℝ × ℝ × ℝ) × ℝ) : ℝ :=
  Z04.hIntegrand 0 z.2 z.1.1 * Z04.hIntegrand 0 z.2 z.1.2.1 * Z04.hIntegrand 0 z.2 z.1.2.2

lemma Φ_meas : Measurable Φ := by
  unfold Φ Z04.hIntegrand; fun_prop

lemma K_meas : Measurable K := by
  unfold K S02.Dfun; fun_prop

lemma Φ_nn {g : ℝ × ℝ × ℝ} (hg : g ∈ Y05.Q) (y : ℝ) : 0 ≤ Φ (g, y) := by
  obtain ⟨ha, hb, hc⟩ := hg
  simp only [mem_Ioi] at ha hb hc
  unfold Φ
  have := Z04.hI_nn 0 y ha
  have := Z04.hI_nn 0 y hb
  have := Z04.hI_nn 0 y hc
  simp only
  positivity

lemma Φ_y {g : ℝ × ℝ × ℝ} (hg : g ∈ Y05.Q) :
    IntegrableOn (fun y => Φ (g, y)) (Ioi 0) ∧
      ∫ y in Ioi (0 : ℝ), Φ (g, y) = Real.sqrt π * K g := by
  obtain ⟨ha, hb, hc⟩ := hg
  simp only [mem_Ioi] at ha hb hc
  obtain ⟨i, v⟩ := prod_int (sq_nonneg g.1) (sq_nonneg g.2.1) (sq_nonneg g.2.2)
  have e : (fun y => Φ (g, y)) = fun y =>
      (g.1⁻¹ * Real.exp (-g.1) * (g.2.1⁻¹ * Real.exp (-g.2.1)) * (g.2.2⁻¹ * Real.exp (-g.2.2)))
        * (E (g.1 ^ 2) y * E (g.2.1 ^ 2) y * E (g.2.2 ^ 2) y) := by
    funext y
    unfold Φ
    simp only
    rw [hI0 ha, hI0 hb, hI0 hc]
    ring
  rw [e]
  refine ⟨i.const_mul _, ?_⟩
  rw [integral_const_mul, v, Real.sqrt_sq ha.le, Real.sqrt_sq hb.le, Real.sqrt_sq hc.le]
  unfold K S02.Dfun
  rw [neg_add, neg_add, Real.exp_add, Real.exp_add]
  field_simp

lemma h3_eq {y : ℝ} (hy : 0 < y) :
    ENNReal.ofReal (Z04.h 0 y ^ 3) = ∫⁻ g in Y05.Q, ENNReal.ofReal (Φ (g, y)) := by
  have hi := Z04.hI_int (ε := 0) (by norm_num) hy
  obtain ⟨i, v⟩ := Z04.prod3 hi hi hi
  have e : Z04.h 0 y ^ 3 = ∫ g in Y05.Q, Φ (g, y) := by
    unfold Φ; simp only; rw [v]; unfold Z04.h; ring
  have i' : IntegrableOn (fun g => Φ (g, y)) Y05.Q := i
  rw [e, ofReal_integral_eq_lintegral_ofReal i']
  · exact (ae_restrict_iff' Y05.Q_meas).mpr (Eventually.of_forall fun g hg => Φ_nn hg y)

lemma lint_h3 : ∫⁻ y in Ioi (0 : ℝ), ENNReal.ofReal (Z04.h 0 y ^ 3) =
    ENNReal.ofReal (Real.sqrt π) * ∫⁻ g in Y05.Q, ENNReal.ofReal (K g) := by
  rw [setLIntegral_congr_fun measurableSet_Ioi (fun y hy => h3_eq hy)]
  rw [lintegral_lintegral_swap ((Φ_meas.comp measurable_swap).ennreal_ofReal.aemeasurable)]
  rw [← lintegral_const_mul _ K_meas.ennreal_ofReal]
  refine setLIntegral_congr_fun Y05.Q_meas (fun g hg => ?_)
  obtain ⟨i, v⟩ := Φ_y hg
  rw [← ofReal_integral_eq_lintegral_ofReal i, v, ENNReal.ofReal_mul (Real.sqrt_nonneg _)]
  exact (ae_restrict_iff' measurableSet_Ioi).mpr (Eventually.of_forall fun y _ => Φ_nn hg y)

lemma sqrt_scale {r : ℝ} (hr : 0 ≤ r) (s : ℝ) :
    Real.sqrt (r ^ 2 * s) = r * Real.sqrt s := by
  rw [Real.sqrt_mul (sq_nonneg r), Real.sqrt_sq hr]

lemma Dfun_scale {r : ℝ} (hr : 0 ≤ r) (a b c : ℝ) :
    S02.Dfun (r * a) (r * b) (r * c) = r * S02.Dfun a b c := by
  unfold S02.Dfun
  rw [show (r * a) ^ 2 + (r * b) ^ 2 + (r * c) ^ 2 = r ^ 2 * (a ^ 2 + b ^ 2 + c ^ 2) by ring,
    show (r * a) ^ 2 + (r * b) ^ 2 = r ^ 2 * (a ^ 2 + b ^ 2) by ring,
    show (r * b) ^ 2 + (r * c) ^ 2 = r ^ 2 * (b ^ 2 + c ^ 2) by ring,
    show (r * a) ^ 2 + (r * c) ^ 2 = r ^ 2 * (a ^ 2 + c ^ 2) by ring,
    sqrt_scale hr, sqrt_scale hr, sqrt_scale hr, sqrt_scale hr]
  ring

lemma K_scale {r : ℝ} (hr : 0 < r) (t : ℝ × ℝ) :
    K (r * t.1, r * t.2, r - r * t.1 - r * t.2) =
      (r ^ 2)⁻¹ * Real.exp (-r) * S02.simplexIntegrand t := by
  unfold K S02.simplexIntegrand
  simp only
  rw [show r - r * t.1 - r * t.2 = r * (1 - t.1 - t.2) by ring,
    show r * t.1 + r * t.2 + r * (1 - t.1 - t.2) = r by ring, Dfun_scale hr.le,
    show r * t.1 * (r * t.2) * (r * (1 - t.1 - t.2)) = r ^ 3 * (t.1 * t.2 * (1 - t.1 - t.2)) by
      ring, ← div_mul_div_comm]
  have : r / r ^ 3 = (r ^ 2)⁻¹ := by field_simp
  rw [this]; ring

lemma mem_Q_iff {r : ℝ} (hr : 0 < r) (t : ℝ × ℝ) :
    (r * t.1, r * t.2, r - r * t.1 - r * t.2) ∈ Y05.Q ↔ t ∈ Y05.Δo := by
  constructor
  · rintro ⟨ha, hb, hc⟩
    simp only [mem_Ioi] at ha hb hc
    refine ⟨pos_of_mul_pos_right ha hr.le, pos_of_mul_pos_right hb hr.le, ?_⟩
    have : 0 < r * (1 - t.1 - t.2) := by linarith
    exact pos_of_mul_pos_right this hr.le
  · rintro ⟨ha, hb, hc⟩
    refine ⟨?_, ?_, ?_⟩ <;> simp only [mem_Ioi]
    · positivity
    · positivity
    · nlinarith

lemma sI_meas : Measurable S02.simplexIntegrand := by
  unfold S02.simplexIntegrand S02.Dfun; fun_prop

/-- `F = 1_Q · ofReal ∘ K`. -/
noncomputable def F : ℝ × ℝ × ℝ → ℝ≥0∞ := Y05.Q.indicator (fun g => ENNReal.ofReal (K g))

lemma F_meas : Measurable F := K_meas.ennreal_ofReal.indicator Y05.Q_meas

lemma inner {r : ℝ} (hr : 0 < r) :
    ∫⁻ u : ℝ × ℝ, F (u.1, u.2, r - u.1 - u.2) =
      ENNReal.ofReal (Real.exp (-r)) * ∫⁻ t in Y05.Δ, ENNReal.ofReal (S02.simplexIntegrand t) := by
  set f : ℝ × ℝ → ℝ≥0∞ := fun u => F (u.1, u.2, r - u.1 - u.2) with hf
  have hfm : Measurable f := F_meas.comp (by fun_prop)
  have hmap := Measure.map_addHaar_smul (volume : Measure (ℝ × ℝ)) hr.ne'
  have hfr : Module.finrank ℝ (ℝ × ℝ) = 2 := by simp
  rw [hfr] at hmap
  have h1 : ∫⁻ t, f (r • t) = ENNReal.ofReal |(r ^ 2)⁻¹| * ∫⁻ u, f u := by
    rw [← lintegral_map hfm (measurable_const_smul r), hmap, lintegral_smul_measure]
    rfl
  have hpt : (fun t : ℝ × ℝ => f (r • t)) = Y05.Δo.indicator
      (fun t => ENNReal.ofReal ((r ^ 2)⁻¹ * Real.exp (-r)) *
        ENNReal.ofReal (S02.simplexIntegrand t)) := by
    funext t
    simp only [hf, F, Prod.smul_fst, Prod.smul_snd, smul_eq_mul]
    by_cases ht : t ∈ Y05.Δo
    · rw [indicator_of_mem ((mem_Q_iff hr t).mpr ht), indicator_of_mem ht, K_scale hr t,
        ENNReal.ofReal_mul (by positivity)]
    · rw [indicator_of_notMem (fun h => ht ((mem_Q_iff hr t).mp h)), indicator_of_notMem ht]
  rw [hpt, lintegral_indicator Y05.Δo_meas,
    lintegral_const_mul _ sI_meas.ennreal_ofReal, setLIntegral_congr Y05.Δo_ae] at h1
  have hr2 : (0 : ℝ) < r ^ 2 := by positivity
  rw [abs_of_pos (inv_pos.mpr hr2)] at h1
  have key : ∫⁻ u, f u = ENNReal.ofReal (r ^ 2) * (ENNReal.ofReal |(r ^ 2)⁻¹| * ∫⁻ u, f u) := by
    rw [abs_of_pos (inv_pos.mpr hr2), ← mul_assoc, ← ENNReal.ofReal_mul hr2.le,
      mul_inv_cancel₀ hr2.ne', ENNReal.ofReal_one, one_mul]
  rw [abs_of_pos (inv_pos.mpr hr2), ← h1, ← mul_assoc, ← ENNReal.ofReal_mul hr2.le,
    ← mul_assoc, mul_inv_cancel₀ hr2.ne', one_mul] at key
  exact key

lemma F_zero {r : ℝ} (hr : r ≤ 0) (u : ℝ × ℝ) : F (u.1, u.2, r - u.1 - u.2) = 0 := by
  unfold F
  rw [indicator_of_notMem]
  rintro ⟨ha, hb, hc⟩
  simp only [mem_Ioi] at ha hb hc
  linarith

lemma lint_exp : ∫⁻ r in Ioi (0 : ℝ), ENNReal.ofReal (Real.exp (-r)) = 1 := by
  rw [← ofReal_integral_eq_lintegral_ofReal (integrableOn_exp_neg_Ioi 0)
    (Eventually.of_forall fun r => (Real.exp_pos _).le), integral_exp_neg_Ioi_zero,
    ENNReal.ofReal_one]

/-- The Dirichlet scaling `g = r σ`, Lebesgue-integral form. -/
lemma dirichlet : ∫⁻ g in Y05.Q, ENNReal.ofReal (K g) =
    ∫⁻ t in Y05.Δ, ENNReal.ofReal (S02.simplexIntegrand t) := by
  rw [← lintegral_indicator Y05.Q_meas]
  change ∫⁻ g, F g = _
  have hA := measurePreserving_prodAssoc (volume : Measure ℝ) (volume : Measure ℝ)
    (volume : Measure ℝ)
  have e1 : ∫⁻ g, F g = ∫⁻ v : (ℝ × ℝ) × ℝ, F (v.1.1, v.1.2, v.2) ∂(volume.prod volume) := by
    rw [show (volume : Measure (ℝ × ℝ × ℝ)) = volume.prod (volume.prod volume) from rfl,
      ← hA.lintegral_comp F_meas]
    rfl
  have hS : MeasurePreserving (fun w : (ℝ × ℝ) × ℝ => (w.1, w.2 - w.1.1 - w.1.2))
      (volume.prod volume) (volume.prod volume) := by
    refine MeasurePreserving.skew_product (MeasurePreserving.id _)
      (g := fun (u : ℝ × ℝ) (r : ℝ) => r - u.1 - u.2) (by fun_prop) (Eventually.of_forall fun u => ?_)
    have := map_sub_right_eq_self (volume : Measure ℝ) (u.1 + u.2)
    simpa [sub_sub] using this
  have hGm : Measurable (fun v : (ℝ × ℝ) × ℝ => F (v.1.1, v.1.2, v.2)) :=
    F_meas.comp (by fun_prop)
  have e2 : ∫⁻ v : (ℝ × ℝ) × ℝ, F (v.1.1, v.1.2, v.2) ∂(volume.prod volume) =
      ∫⁻ a : (ℝ × ℝ) × ℝ, F (a.1.1, a.1.2, a.2 - a.1.1 - a.1.2) ∂(volume.prod volume) :=
    (hS.lintegral_comp hGm).symm
  have hHm : Measurable (fun a : (ℝ × ℝ) × ℝ => F (a.1.1, a.1.2, a.2 - a.1.1 - a.1.2)) :=
    F_meas.comp (by fun_prop)
  rw [e1, e2, lintegral_prod_symm' _ hHm]
  have e3 : ∀ r : ℝ, ∫⁻ u : ℝ × ℝ, F (u.1, u.2, r - u.1 - u.2) =
      (Ioi (0 : ℝ)).indicator (fun r => ENNReal.ofReal (Real.exp (-r)) *
        ∫⁻ t in Y05.Δ, ENNReal.ofReal (S02.simplexIntegrand t)) r := by
    intro r
    by_cases hr : 0 < r
    · rw [indicator_of_mem (show r ∈ Ioi 0 from hr), inner hr]
    · rw [indicator_of_notMem (show r ∉ Ioi 0 from hr)]
      simp only [F_zero (not_lt.mp hr), lintegral_zero]
  simp_rw [e3]
  rw [lintegral_indicator measurableSet_Ioi, lintegral_mul_const _ (by fun_prop), lint_exp, one_mul]

lemma sI_nn {t : ℝ × ℝ} (ht : t ∈ Y05.Δ) : 0 ≤ S02.simplexIntegrand t := by
  obtain ⟨ha, hb, hab⟩ := ht
  unfold S02.simplexIntegrand
  exact div_nonneg (S02.hlawka3 _ _ _ ha hb (by linarith))
    (mul_nonneg (mul_nonneg ha hb) (by linarith))

lemma sI_ae_nn : 0 ≤ᵐ[volume.restrict Y05.Δ] S02.simplexIntegrand :=
  (ae_restrict_iff' Y05.Δ_meas).mpr (Eventually.of_forall fun _ ht => sI_nn ht)

lemma h3_int : IntegrableOn (fun y => Z04.h 0 y ^ 3) (Ioi 0) :=
  (R04.proof_6 3 (by norm_num) le_rfl).1 0 ⟨le_rfl, by norm_num⟩

lemma key : ENNReal.ofReal (∫ y in Ioi (0 : ℝ), Z04.h 0 y ^ 3) =
    ENNReal.ofReal (Real.sqrt π) * ∫⁻ t in Y05.Δ, ENNReal.ofReal (S02.simplexIntegrand t) := by
  rw [ofReal_integral_eq_lintegral_ofReal h3_int
    (Eventually.of_forall fun y => pow_nonneg (Z04.h_nn 0 y) 3), lint_h3, dirichlet]

lemma lint_ne_top : ∫⁻ t in Y05.Δ, ENNReal.ofReal (S02.simplexIntegrand t) ≠ ⊤ := by
  intro hX
  have k := key
  rw [hX, ENNReal.mul_top (by simpa using Real.sqrt_pos.mpr Real.pi_pos)] at k
  exact ENNReal.ofReal_ne_top k

lemma sI_int : IntegrableOn S02.simplexIntegrand Y05.Δ :=
  ⟨sI_meas.aestronglyMeasurable,
    (hasFiniteIntegral_iff_ofReal sI_ae_nn).mpr (lt_top_iff_ne_top.mpr lint_ne_top)⟩

lemma Q0_eq : Z04.Q 0 = -∫ t in Y05.Δ, S02.simplexIntegrand t := by
  have k := congrArg ENNReal.toReal key
  rw [ENNReal.toReal_ofReal (setIntegral_nonneg measurableSet_Ioi
      (fun y _ => pow_nonneg (Z04.h_nn 0 y) 3)), ENNReal.toReal_mul,
    ENNReal.toReal_ofReal (Real.sqrt_nonneg _),
    ← integral_eq_lintegral_of_nonneg_ae sI_ae_nn sI_meas.aestronglyMeasurable] at k
  unfold Z04.Q
  rw [k, show (1 : ℝ) + 2 * 0 = 1 by norm_num, Real.Gamma_one, mul_one]
  have : 0 < Real.sqrt π := Real.sqrt_pos.mpr Real.pi_pos
  field_simp

theorem main : S06.stmt_1 :=
  ⟨fun x₁ x₂ x₃ h1 h2 h3 => S02.hlawka3 x₁ x₂ x₃ h1 h2 h3, sI_int, Q0_eq⟩

end S06.PA

/-! ### Merged from work/S06/PartB.lean -/
/-!
# S06, part B: `stmt_1 → stmt_2` (Lemma `g3:lem:step1`)

Route (no regularisation needed). On the open simplex, with `σ = (t₁, t₂, 1-t₁-t₂)`,
`A_{ij} = √(σᵢ²+σⱼ²)`, `ρ = √(σ₁²+σ₂²+σ₃²)`, one has the exact pointwise splitting
`D(σ)/(σ₁σ₂σ₃) = T₁ + T₂ + T₃`, `T_k = (2/σ_k)(1/(σᵢ+σⱼ+A_{ij}) - 1/(1+ρ))` (`{i,j,k} = {1,2,3}`),
with every `T_k ≥ 0` (hence integrable, being dominated by the integrable sum).
The coordinates `σ_k = r/(1+r)`, `{σᵢ, σⱼ} = {s/(1+r), (1-s)/(1+r)}` (Jacobian `1/(1+r)³`) turn
each `∫_Δ T_k` into `-2 ∫₀¹ ∫₀^∞ Jintegrand (q(s)) r dr ds`.
-/

namespace S06.PB

open Real MeasureTheory Set Filter Topology

/-! ## 2×2 matrices as continuous linear maps -/

noncomputable def M (a b c d : ℝ) : ℝ × ℝ →L[ℝ] ℝ × ℝ :=
  (Matrix.toLin (Module.Basis.finTwoProd ℝ) (Module.Basis.finTwoProd ℝ) !![a, b; c, d]).toContinuousLinearMap

lemma M_det (a b c d : ℝ) : (M a b c d).det = a * d - b * c := by
  unfold M
  simp only [ContinuousLinearMap.det, LinearMap.coe_toContinuousLinearMap, LinearMap.det_toLin,
    Matrix.det_fin_two_of]

lemma M_apply (a b c d : ℝ) (v : ℝ × ℝ) : M a b c d v = (a * v.1 + b * v.2, c * v.1 + d * v.2) := by
  unfold M
  simp

/-! ## the parametrisations -/

/-- the open domain of `(s, r)` -/
def U : Set (ℝ × ℝ) := Ioo (0 : ℝ) 1 ×ˢ Ioi (0 : ℝ)

lemma U_meas : MeasurableSet U := measurableSet_Ioo.prod measurableSet_Ioi

noncomputable def ψ (p : ℝ × ℝ) : ℝ := (1 + p.2)⁻¹

lemma hasFDerivAt_ψ {p : ℝ × ℝ} (hp : 0 < p.2) :
    HasFDerivAt ψ ((-((1 + p.2) ^ 2)⁻¹) • ContinuousLinearMap.snd ℝ ℝ ℝ) p := by
  have h1 : HasFDerivAt (fun q : ℝ × ℝ => 1 + q.2) (ContinuousLinearMap.snd ℝ ℝ ℝ) p :=
    (hasFDerivAt_snd).const_add 1
  have := (hasDerivAt_inv (show (1 + p.2) ≠ 0 by linarith)).comp_hasFDerivAt p h1
  exact this

/-- `σ₃ = r/(1+r)`: `t = (s/(1+r), (1-s)/(1+r))` -/
noncomputable def Φ3 (p : ℝ × ℝ) : ℝ × ℝ := (p.1 * ψ p, (1 - p.1) * ψ p)
/-- `σ₁ = r/(1+r)`: `t = (r/(1+r), s/(1+r))` -/
noncomputable def Φ1 (p : ℝ × ℝ) : ℝ × ℝ := (p.2 * ψ p, p.1 * ψ p)
/-- `σ₂ = r/(1+r)`: `t = (s/(1+r), r/(1+r))` -/
noncomputable def Φ2 (p : ℝ × ℝ) : ℝ × ℝ := (p.1 * ψ p, p.2 * ψ p)

noncomputable def M3 (p : ℝ × ℝ) : ℝ × ℝ →L[ℝ] ℝ × ℝ :=
  M (1 + p.2)⁻¹ (-(p.1 / (1 + p.2) ^ 2)) (-(1 + p.2)⁻¹) (-((1 - p.1) / (1 + p.2) ^ 2))
noncomputable def M1 (p : ℝ × ℝ) : ℝ × ℝ →L[ℝ] ℝ × ℝ :=
  M 0 (1 / (1 + p.2) ^ 2) (1 + p.2)⁻¹ (-(p.1 / (1 + p.2) ^ 2))
noncomputable def M2 (p : ℝ × ℝ) : ℝ × ℝ →L[ℝ] ℝ × ℝ :=
  M (1 + p.2)⁻¹ (-(p.1 / (1 + p.2) ^ 2)) 0 (1 / (1 + p.2) ^ 2)

lemma hasFDerivAt_Φ3 {p : ℝ × ℝ} (hp : 0 < p.2) : HasFDerivAt Φ3 (M3 p) p := by
  have h1 := (hasFDerivAt_fst : HasFDerivAt (Prod.fst : ℝ × ℝ → ℝ) (ContinuousLinearMap.fst ℝ ℝ ℝ) p).mul (hasFDerivAt_ψ hp)
  have h2 := ((hasFDerivAt_const (1 : ℝ) p).sub (hasFDerivAt_fst : HasFDerivAt (Prod.fst : ℝ × ℝ → ℝ) (ContinuousLinearMap.fst ℝ ℝ ℝ) p)).mul
    (hasFDerivAt_ψ hp)
  have h := h1.prodMk h2
  convert h using 1
  · funext x; simp [Φ3, Pi.mul_apply, Pi.sub_apply]
  have hr : (1 + p.2) ≠ 0 := by linarith
  apply ContinuousLinearMap.ext; intro v
  simp only [M3, M_apply, ContinuousLinearMap.prod_apply, ContinuousLinearMap.add_apply,
    ContinuousLinearMap.smul_apply, ContinuousLinearMap.coe_fst', ContinuousLinearMap.coe_snd',
    ContinuousLinearMap.sub_apply, ContinuousLinearMap.zero_apply, smul_eq_mul, ψ, Pi.sub_apply]
  refine Prod.ext ?_ ?_ <;> simp only <;> field_simp <;> ring

lemma hasFDerivAt_Φ1 {p : ℝ × ℝ} (hp : 0 < p.2) : HasFDerivAt Φ1 (M1 p) p := by
  have h1 := (hasFDerivAt_snd : HasFDerivAt (Prod.snd : ℝ × ℝ → ℝ) (ContinuousLinearMap.snd ℝ ℝ ℝ) p).mul (hasFDerivAt_ψ hp)
  have h2 := (hasFDerivAt_fst : HasFDerivAt (Prod.fst : ℝ × ℝ → ℝ) (ContinuousLinearMap.fst ℝ ℝ ℝ) p).mul (hasFDerivAt_ψ hp)
  have h := h1.prodMk h2
  convert h using 1
  · funext x; simp [Φ1, Pi.mul_apply]
  have hr : (1 + p.2) ≠ 0 := by linarith
  apply ContinuousLinearMap.ext; intro v
  simp only [M1, M_apply, ContinuousLinearMap.prod_apply, ContinuousLinearMap.add_apply,
    ContinuousLinearMap.smul_apply, ContinuousLinearMap.coe_fst', ContinuousLinearMap.coe_snd',
    smul_eq_mul, ψ]
  refine Prod.ext ?_ ?_ <;> simp only <;> field_simp <;> ring

lemma hasFDerivAt_Φ2 {p : ℝ × ℝ} (hp : 0 < p.2) : HasFDerivAt Φ2 (M2 p) p := by
  have h1 := (hasFDerivAt_fst : HasFDerivAt (Prod.fst : ℝ × ℝ → ℝ) (ContinuousLinearMap.fst ℝ ℝ ℝ) p).mul (hasFDerivAt_ψ hp)
  have h2 := (hasFDerivAt_snd : HasFDerivAt (Prod.snd : ℝ × ℝ → ℝ) (ContinuousLinearMap.snd ℝ ℝ ℝ) p).mul (hasFDerivAt_ψ hp)
  have h := h1.prodMk h2
  convert h using 1
  · funext x; simp [Φ2, Pi.mul_apply]
  have hr : (1 + p.2) ≠ 0 := by linarith
  apply ContinuousLinearMap.ext; intro v
  simp only [M2, M_apply, ContinuousLinearMap.prod_apply, ContinuousLinearMap.add_apply,
    ContinuousLinearMap.smul_apply, ContinuousLinearMap.coe_fst', ContinuousLinearMap.coe_snd',
    smul_eq_mul, ψ]
  refine Prod.ext ?_ ?_ <;> simp only <;> field_simp <;> ring

lemma absdet3 {p : ℝ × ℝ} (hp : 0 < p.2) : |(M3 p).det| = ((1 + p.2) ^ 3)⁻¹ := by
  have hr : (1 + p.2) ≠ 0 := by linarith
  rw [M3, M_det, show (1 + p.2)⁻¹ * -((1 - p.1) / (1 + p.2) ^ 2) - -(p.1 / (1 + p.2) ^ 2) * -(1 + p.2)⁻¹
    = -((1 + p.2) ^ 3)⁻¹ by field_simp; ring, abs_neg, abs_of_pos (by positivity)]

lemma absdet1 {p : ℝ × ℝ} (hp : 0 < p.2) : |(M1 p).det| = ((1 + p.2) ^ 3)⁻¹ := by
  have hr : (1 + p.2) ≠ 0 := by linarith
  rw [M1, M_det, show 0 * -(p.1 / (1 + p.2) ^ 2) - 1 / (1 + p.2) ^ 2 * (1 + p.2)⁻¹
    = -((1 + p.2) ^ 3)⁻¹ by field_simp; ring, abs_neg, abs_of_pos (by positivity)]

lemma absdet2 {p : ℝ × ℝ} (hp : 0 < p.2) : |(M2 p).det| = ((1 + p.2) ^ 3)⁻¹ := by
  have hr : (1 + p.2) ≠ 0 := by linarith
  rw [M2, M_det, show (1 + p.2)⁻¹ * (1 / (1 + p.2) ^ 2) - -(p.1 / (1 + p.2) ^ 2) * 0
    = ((1 + p.2) ^ 3)⁻¹ by rw [mul_zero, sub_zero]; field_simp, abs_of_pos (by positivity)]

/-! ## images and injectivity -/

lemma image_of_inv {Φ G : ℝ × ℝ → ℝ × ℝ} (h1 : ∀ p ∈ U, Φ p ∈ Y05.Δo)
    (h2 : ∀ t ∈ Y05.Δo, G t ∈ U ∧ Φ (G t) = t) : Φ '' U = Y05.Δo := by
  ext t
  constructor
  · rintro ⟨p, hp, rfl⟩; exact h1 p hp
  · intro ht; exact ⟨G t, (h2 t ht).1, (h2 t ht).2⟩

lemma injOn_of_inv {Φ G : ℝ × ℝ → ℝ × ℝ} (h : ∀ p ∈ U, G (Φ p) = p) : InjOn Φ U := by
  intro p hp p' hp' e
  rw [← h p hp, ← h p' hp', e]

lemma mem_U {p : ℝ × ℝ} (hp : p ∈ U) : 0 < p.1 ∧ p.1 < 1 ∧ 0 < p.2 := ⟨hp.1.1, hp.1.2, hp.2⟩

noncomputable def G3 (t : ℝ × ℝ) : ℝ × ℝ := (t.1 / (t.1 + t.2), (1 - t.1 - t.2) / (t.1 + t.2))
noncomputable def G1 (t : ℝ × ℝ) : ℝ × ℝ := (t.2 / (1 - t.1), t.1 / (1 - t.1))
noncomputable def G2 (t : ℝ × ℝ) : ℝ × ℝ := (t.1 / (1 - t.2), t.2 / (1 - t.2))

lemma image3 : Φ3 '' U = Y05.Δo := by
  refine image_of_inv (G := G3) (fun p hp => ?_) (fun t ht => ?_)
  · obtain ⟨a, b, c⟩ := mem_U hp
    have hr : (0:ℝ) < 1 + p.2 := by linarith
    refine ⟨by simp only [Φ3, ψ]; positivity, by simp only [Φ3, ψ]; apply mul_pos (by linarith) (by positivity), ?_⟩
    simp only [Φ3, ψ]
    rw [show 1 - p.1 * (1 + p.2)⁻¹ - (1 - p.1) * (1 + p.2)⁻¹ = p.2 / (1 + p.2) by field_simp; ring]
    positivity
  · obtain ⟨a, b, c⟩ := ht
    have hu : 0 < t.1 + t.2 := by linarith
    refine ⟨⟨⟨by simp only [G3]; positivity, ?_⟩, ?_⟩, ?_⟩
    · simp only [G3]; rw [div_lt_one hu]; linarith
    · simp only [G3, mem_Ioi]; apply div_pos (by linarith) hu
    · simp only [Φ3, G3, ψ]
      have : 1 + (1 - t.1 - t.2) / (t.1 + t.2) = 1 / (t.1 + t.2) := by field_simp; ring
      rw [this]
      refine Prod.ext ?_ ?_ <;> simp only <;> field_simp <;> ring

lemma inj3 : InjOn Φ3 U := by
  refine injOn_of_inv (G := G3) (fun p hp => ?_)
  obtain ⟨a, b, c⟩ := mem_U hp
  have hr : (0:ℝ) < 1 + p.2 := by linarith
  simp only [Φ3, G3, ψ]
  have e : p.1 * (1 + p.2)⁻¹ + (1 - p.1) * (1 + p.2)⁻¹ = (1 + p.2)⁻¹ := by ring
  rw [e]
  refine Prod.ext ?_ ?_ <;> simp only <;> field_simp <;> ring

lemma image1 : Φ1 '' U = Y05.Δo := by
  refine image_of_inv (G := G1) (fun p hp => ?_) (fun t ht => ?_)
  · obtain ⟨a, b, c⟩ := mem_U hp
    have hr : (0:ℝ) < 1 + p.2 := by linarith
    refine ⟨by simp only [Φ1, ψ]; positivity, by simp only [Φ1, ψ]; positivity, ?_⟩
    simp only [Φ1, ψ]
    rw [show 1 - p.2 * (1 + p.2)⁻¹ - p.1 * (1 + p.2)⁻¹ = (1 - p.1) / (1 + p.2) by field_simp; ring]
    apply div_pos (by linarith) hr
  · obtain ⟨a, b, c⟩ := ht
    have hw : 0 < 1 - t.1 := by linarith
    refine ⟨⟨⟨by simp only [G1]; positivity, ?_⟩, ?_⟩, ?_⟩
    · simp only [G1]; rw [div_lt_one hw]; linarith
    · simp only [G1, mem_Ioi]; positivity
    · simp only [Φ1, G1, ψ]
      have : 1 + t.1 / (1 - t.1) = 1 / (1 - t.1) := by field_simp; ring
      rw [this]
      refine Prod.ext ?_ ?_ <;> simp only <;> field_simp

lemma inj1 : InjOn Φ1 U := by
  refine injOn_of_inv (G := G1) (fun p hp => ?_)
  obtain ⟨a, b, c⟩ := mem_U hp
  have hr : (0:ℝ) < 1 + p.2 := by linarith
  simp only [Φ1, G1, ψ]
  have e : 1 - p.2 * (1 + p.2)⁻¹ = (1 + p.2)⁻¹ := by field_simp; ring
  rw [e]
  refine Prod.ext ?_ ?_ <;> simp only <;> field_simp

lemma image2 : Φ2 '' U = Y05.Δo := by
  refine image_of_inv (G := G2) (fun p hp => ?_) (fun t ht => ?_)
  · obtain ⟨a, b, c⟩ := mem_U hp
    have hr : (0:ℝ) < 1 + p.2 := by linarith
    refine ⟨by simp only [Φ2, ψ]; positivity, by simp only [Φ2, ψ]; positivity, ?_⟩
    simp only [Φ2, ψ]
    rw [show 1 - p.1 * (1 + p.2)⁻¹ - p.2 * (1 + p.2)⁻¹ = (1 - p.1) / (1 + p.2) by field_simp; ring]
    apply div_pos (by linarith) hr
  · obtain ⟨a, b, c⟩ := ht
    have hw : 0 < 1 - t.2 := by linarith
    refine ⟨⟨⟨by simp only [G2]; positivity, ?_⟩, ?_⟩, ?_⟩
    · simp only [G2]; rw [div_lt_one hw]; linarith
    · simp only [G2, mem_Ioi]; positivity
    · simp only [Φ2, G2, ψ]
      have : 1 + t.2 / (1 - t.2) = 1 / (1 - t.2) := by field_simp; ring
      rw [this]
      refine Prod.ext ?_ ?_ <;> simp only <;> field_simp

lemma inj2 : InjOn Φ2 U := by
  refine injOn_of_inv (G := G2) (fun p hp => ?_)
  obtain ⟨a, b, c⟩ := mem_U hp
  have hr : (0:ℝ) < 1 + p.2 := by linarith
  simp only [Φ2, G2, ψ]
  have e : 1 - p.2 * (1 + p.2)⁻¹ = (1 + p.2)⁻¹ := by field_simp; ring
  rw [e]
  refine Prod.ext ?_ ?_ <;> simp only <;> field_simp

/-! ## the splitting `D/Π = T₁ + T₂ + T₃` -/

/-- one piece: `(2/c)(1/(a+b+√(a²+b²)) - 1/(1+√P))` -/
noncomputable def piece (a b c P : ℝ) : ℝ :=
  2 / c * (1 / (a + b + Real.sqrt (a ^ 2 + b ^ 2)) - 1 / (1 + Real.sqrt P))

noncomputable def P3 (t : ℝ × ℝ) : ℝ := t.1 ^ 2 + t.2 ^ 2 + (1 - t.1 - t.2) ^ 2

noncomputable def T3 (t : ℝ × ℝ) : ℝ := piece t.1 t.2 (1 - t.1 - t.2) (P3 t)
noncomputable def T1 (t : ℝ × ℝ) : ℝ := piece t.2 (1 - t.1 - t.2) t.1 (P3 t)
noncomputable def T2 (t : ℝ × ℝ) : ℝ := piece t.1 (1 - t.1 - t.2) t.2 (P3 t)

lemma inv_pair {x y A : ℝ} (hx : 0 < x) (hy : 0 < y) (hA0 : 0 ≤ A) (hA : A ^ 2 = x ^ 2 + y ^ 2) :
    1 / (x + y + A) = (x + y - A) / (2 * x * y) := by
  rw [div_eq_div_iff (by positivity) (by positivity)]
  linear_combination hA

lemma split {x y : ℝ} (hx : 0 < x) (hy : 0 < y) (hz : 0 < 1 - x - y) :
    S02.Dfun x y (1 - x - y) / (x * y * (1 - x - y)) =
      piece y (1 - x - y) x (x ^ 2 + y ^ 2 + (1 - x - y) ^ 2)
      + piece x (1 - x - y) y (x ^ 2 + y ^ 2 + (1 - x - y) ^ 2)
      + piece x y (1 - x - y) (x ^ 2 + y ^ 2 + (1 - x - y) ^ 2) := by
  set z := 1 - x - y with hzdef
  set A := Real.sqrt (x ^ 2 + y ^ 2) with hAd
  set B := Real.sqrt (y ^ 2 + z ^ 2) with hBd
  set C := Real.sqrt (x ^ 2 + z ^ 2) with hCd
  set R := Real.sqrt (x ^ 2 + y ^ 2 + z ^ 2) with hRd
  have hA : A ^ 2 = x ^ 2 + y ^ 2 := Real.sq_sqrt (by positivity)
  have hB : B ^ 2 = y ^ 2 + z ^ 2 := Real.sq_sqrt (by positivity)
  have hC : C ^ 2 = x ^ 2 + z ^ 2 := Real.sq_sqrt (by positivity)
  have hR : R ^ 2 = x ^ 2 + y ^ 2 + z ^ 2 := Real.sq_sqrt (by positivity)
  have hR0 : 0 ≤ R := Real.sqrt_nonneg _
  have k4 : 1 / (1 + R) = (1 - R) / (2 * (x * y + y * z + x * z)) := by
    rw [div_eq_div_iff (by positivity) (by positivity)]
    rw [hzdef] at hR ⊢
    linear_combination hR
  unfold piece S02.Dfun
  rw [← hAd, ← hBd, ← hCd, ← hRd, k4, inv_pair hx hy (Real.sqrt_nonneg _) hA,
    inv_pair hy hz (Real.sqrt_nonneg _) hB, inv_pair hx hz (Real.sqrt_nonneg _) hC]
  have : y ^ 2 + x ^ 2 = x ^ 2 + y ^ 2 := by ring
  field_simp
  rw [hzdef]
  ring

lemma piece_nonneg {a b c P : ℝ} (ha : 0 < a) (hb : 0 < b) (hc : 0 < c) (hab : a + b ≤ 1)
    (hP : a ^ 2 + b ^ 2 ≤ P) : 0 ≤ piece a b c P := by
  unfold piece
  apply mul_nonneg (by positivity)
  rw [sub_nonneg]
  apply one_div_le_one_div_of_le (by positivity)
  have := Real.sqrt_le_sqrt hP
  linarith

lemma simplex_split {t : ℝ × ℝ} (ht : t ∈ Y05.Δo) :
    S02.simplexIntegrand t = T1 t + T2 t + T3 t := by
  obtain ⟨a, b, c⟩ := ht
  exact split a b c

lemma T_nonneg {t : ℝ × ℝ} (ht : t ∈ Y05.Δo) : 0 ≤ T1 t ∧ 0 ≤ T2 t ∧ 0 ≤ T3 t := by
  obtain ⟨a, b, c⟩ := ht
  refine ⟨piece_nonneg b c a (by linarith) ?_, piece_nonneg a c b (by linarith) ?_,
    piece_nonneg a b c (by linarith) ?_⟩ <;> unfold P3 <;> nlinarith

lemma meas_piece : ∀ (f g h k : ℝ × ℝ → ℝ), Measurable f → Measurable g → Measurable h →
    Measurable k → Measurable (fun t => piece (f t) (g t) (h t) (k t)) := by
  intro f g h k hf hg hh hk
  unfold piece
  fun_prop

lemma meas_T3 : Measurable T3 := by
  unfold T3 P3; apply meas_piece <;> fun_prop
lemma meas_T1 : Measurable T1 := by
  unfold T1 P3; apply meas_piece <;> fun_prop
lemma meas_T2 : Measurable T2 := by
  unfold T2 P3; apply meas_piece <;> fun_prop

/-! ## the pull-back of each piece -/

lemma core {r q R ψ' x y z P : ℝ} (hr : 0 < r) (hq : 0 < q) (hR : 0 < R) (hψ : ψ' = (1 + r)⁻¹)
    (hz : z = r * ψ') (hxy : x + y = ψ') (hA : Real.sqrt (x ^ 2 + y ^ 2) = ψ' * q)
    (hρ : Real.sqrt P = ψ' * R) :
    ((1 + r) ^ 3)⁻¹ * piece x y z P = -2 * ((1 / (1 + r + R) - 1 / (1 + q)) / (r * (1 + r))) := by
  unfold piece
  rw [hA, hρ, hxy, hz, hψ]
  have h1 : (1 + r) ≠ 0 := by linarith
  field_simp
  ring

lemma sqrt_scale {X Y c : ℝ} (hc : 0 ≤ c) (hY : 0 ≤ Y) (h : X = c ^ 2 * Y) :
    Real.sqrt X = c * Real.sqrt Y := by
  rw [h, Real.sqrt_mul (sq_nonneg c), Real.sqrt_sq hc]

lemma qs_sq (s : ℝ) : S02.qs s ^ 2 = s ^ 2 + (1 - s) ^ 2 := Real.sq_sqrt (by positivity)

/-- the common value of the three pull-backs -/
noncomputable def Jp (p : ℝ × ℝ) : ℝ := -2 * S02.Jintegrand (S02.qs p.1) p.2

lemma pull3 {p : ℝ × ℝ} (hp : p ∈ U) : |(M3 p).det| • T3 (Φ3 p) = Jp p := by
  obtain ⟨a, b, c⟩ := mem_U hp
  have hr : (0:ℝ) < 1 + p.2 := by linarith
  rw [absdet3 c, smul_eq_mul, Jp, S02.Jintegrand]
  unfold T3 P3 Φ3 ψ
  simp only
  refine core c (S02.qs_pos p.1) (Real.sqrt_pos.mpr (by positivity)) rfl ?_ (by ring) ?_ ?_
  · field_simp; ring
  · rw [S02.qs]
    exact sqrt_scale (by positivity) (by positivity) (by ring)
  · exact sqrt_scale (by positivity) (by positivity) (by rw [qs_sq]; field_simp; ring)

lemma pull1 {p : ℝ × ℝ} (hp : p ∈ U) : |(M1 p).det| • T1 (Φ1 p) = Jp p := by
  obtain ⟨a, b, c⟩ := mem_U hp
  have hr : (0:ℝ) < 1 + p.2 := by linarith
  rw [absdet1 c, smul_eq_mul, Jp, S02.Jintegrand]
  unfold T1 P3 Φ1 ψ
  simp only
  refine core c (S02.qs_pos p.1) (Real.sqrt_pos.mpr (by positivity)) rfl rfl (by field_simp; ring) ?_ ?_
  · rw [S02.qs]
    exact sqrt_scale (by positivity) (by positivity) (by field_simp; ring)
  · exact sqrt_scale (by positivity) (by positivity) (by rw [qs_sq]; field_simp; ring)

lemma pull2 {p : ℝ × ℝ} (hp : p ∈ U) : |(M2 p).det| • T2 (Φ2 p) = Jp p := by
  obtain ⟨a, b, c⟩ := mem_U hp
  have hr : (0:ℝ) < 1 + p.2 := by linarith
  rw [absdet2 c, smul_eq_mul, Jp, S02.Jintegrand]
  unfold T2 P3 Φ2 ψ
  simp only
  refine core c (S02.qs_pos p.1) (Real.sqrt_pos.mpr (by positivity)) rfl rfl (by field_simp; ring) ?_ ?_
  · rw [S02.qs]
    exact sqrt_scale (by positivity) (by positivity) (by field_simp; ring)
  · exact sqrt_scale (by positivity) (by positivity) (by rw [qs_sq]; field_simp; ring)

/-! ## integrals -/

lemma int3 : ∫ t in Y05.Δo, T3 t = ∫ p in U, Jp p := by
  have hV := integral_image_eq_integral_abs_det_fderiv_smul (μ := volume) U_meas
    (fun p hp => (hasFDerivAt_Φ3 hp.2).hasFDerivWithinAt) inj3 T3
  rw [image3] at hV
  rw [hV]
  exact setIntegral_congr_fun U_meas (fun p hp => pull3 hp)

lemma int1 : ∫ t in Y05.Δo, T1 t = ∫ p in U, Jp p := by
  have hV := integral_image_eq_integral_abs_det_fderiv_smul (μ := volume) U_meas
    (fun p hp => (hasFDerivAt_Φ1 hp.2).hasFDerivWithinAt) inj1 T1
  rw [image1] at hV
  rw [hV]
  exact setIntegral_congr_fun U_meas (fun p hp => pull1 hp)

lemma int2 : ∫ t in Y05.Δo, T2 t = ∫ p in U, Jp p := by
  have hV := integral_image_eq_integral_abs_det_fderiv_smul (μ := volume) U_meas
    (fun p hp => (hasFDerivAt_Φ2 hp.2).hasFDerivWithinAt) inj2 T2
  rw [image2] at hV
  rw [hV]
  exact setIntegral_congr_fun U_meas (fun p hp => pull2 hp)

lemma Jp_integrable (h3 : IntegrableOn T3 Y05.Δo) : IntegrableOn Jp U := by
  have hI := integrableOn_image_iff_integrableOn_abs_det_fderiv_smul (μ := volume) U_meas
    (fun p hp => (hasFDerivAt_Φ3 hp.2).hasFDerivWithinAt) inj3 T3
  rw [image3] at hI
  exact (hI.mp h3).congr_fun (fun p hp => pull3 hp) U_meas

lemma T_integrable (hint : IntegrableOn S02.simplexIntegrand Y05.Δ) :
    IntegrableOn T1 Y05.Δo ∧ IntegrableOn T2 Y05.Δo ∧ IntegrableOn T3 Y05.Δo := by
  have h : IntegrableOn S02.simplexIntegrand Y05.Δo :=
    hint.mono_set (fun t ht => ⟨ht.1.le, ht.2.1.le, by linarith [ht.2.2]⟩)
  have hb : ∀ (T : ℝ × ℝ → ℝ), Measurable T → (∀ t ∈ Y05.Δo, 0 ≤ T t ∧ T t ≤ S02.simplexIntegrand t) →
      IntegrableOn T Y05.Δo := by
    intro T hT hle
    refine Integrable.mono' h hT.aestronglyMeasurable ?_
    refine (ae_restrict_iff' Y05.Δo_meas).2 (Eventually.of_forall (fun t ht => ?_))
    rw [Real.norm_eq_abs, abs_of_nonneg (hle t ht).1]
    exact (hle t ht).2
  refine ⟨hb T1 meas_T1 (fun t ht => ?_), hb T2 meas_T2 (fun t ht => ?_),
    hb T3 meas_T3 (fun t ht => ?_)⟩ <;>
  · obtain ⟨n1, n2, n3⟩ := T_nonneg ht
    rw [simplex_split ht]
    refine ⟨by assumption, by linarith⟩

lemma int_simplex (hint : IntegrableOn S02.simplexIntegrand Y05.Δ) :
    ∫ t in Y05.Δ, S02.simplexIntegrand t = 3 * ∫ p in U, Jp p := by
  obtain ⟨i1, i2, i3⟩ := T_integrable hint
  rw [← setIntegral_congr_set Y05.Δo_ae,
    setIntegral_congr_fun Y05.Δo_meas (fun t ht => simplex_split ht),
    integral_add (f := fun x => T1 x + T2 x) (g := T3) (i1.add i2) i3, integral_add i1 i2, int1, int2, int3]
  ring

lemma int_U (h : IntegrableOn Jp U) :
    ∫ p in U, Jp p = -2 * ∫ s in (0 : ℝ)..1, S02.J (S02.qs s) := by
  have e : (fun p : ℝ × ℝ => S02.Jintegrand (S02.qs p.1) p.2) = fun p => (-1 / 2) * Jp p := by
    funext p; unfold Jp; ring
  have hf : IntegrableOn (fun p : ℝ × ℝ => S02.Jintegrand (S02.qs p.1) p.2)
      (Ioo (0 : ℝ) 1 ×ˢ Ioi (0 : ℝ)) (volume.prod volume) := by
    rw [e, ← Measure.volume_eq_prod]; exact h.const_mul _
  have hF := setIntegral_prod (fun p : ℝ × ℝ => S02.Jintegrand (S02.qs p.1) p.2) hf
  rw [← Measure.volume_eq_prod] at hF
  have hJp : (fun p : ℝ × ℝ => Jp p) = fun p => -2 * S02.Jintegrand (S02.qs p.1) p.2 := rfl
  rw [hJp, integral_const_mul, U, hF, intervalIntegral.integral_of_le zero_le_one,
    integral_Ioc_eq_integral_Ioo]
  rfl

theorem main : S06.stmt_1 → S06.stmt_2 := by
  intro h1
  obtain ⟨-, hint, hQ⟩ := h1
  refine ⟨fun _ hq _ => S02.integrableOn_Jintegrand (le_trans (by positivity) hq), S02.ii_J, ?_⟩
  have hJ := Jp_integrable (T_integrable hint).2.2
  rw [hQ, int_simplex hint, int_U hJ]
  ring

end S06.PB

/-! ### Merged from work/S06/PartC.lean -/
/-!
# S06 Part C: the closed form of `P₁` (`S06.stmt_6`)

Route: substitute `s ↦ v = 2s - 1 + √2 q(s)` (so `q = (v²+1)/(2√2 v)`, `v ∈ [b, a]`), write down an
explicit real primitive `G` on `[b, a]` (real dilogarithms and logarithms, plus
`T(v) = Re Li₂((a-v)/(a+i))` from the complex power series), and evaluate `G(a) - G(b)` with the
dilogarithm identities of D00/D02 and one complex reflection at `(1-ib)/√2`.
-/

namespace S06.PC

open Real MeasureTheory Set Filter Topology

/-! ## The complex dilogarithm series on the unit disc -/

/-- `Li₂(z) = Σ z^{n+1}/(n+1)²`. -/
noncomputable def Lc (z : ℂ) : ℂ := ∑' n : ℕ, z ^ (n + 1) / ((n : ℂ) + 1) ^ 2

/-- `Σ zⁿ/(n+1)`, the derivative of `Lc`. -/
noncomputable def Sc (z : ℂ) : ℂ := ∑' n : ℕ, z ^ n / ((n : ℂ) + 1)

lemma norm_natp1 (n : ℕ) : ‖((n : ℂ) + 1)‖ = (n : ℝ) + 1 := by
  rw [show ((n : ℂ) + 1) = (((n : ℝ) + 1 : ℝ) : ℂ) by push_cast; ring, Complex.norm_real,
    Real.norm_of_nonneg (by positivity)]

lemma summable_Sc {z : ℂ} (hz : ‖z‖ < 1) : Summable (fun n : ℕ => z ^ n / ((n : ℂ) + 1)) := by
  refine Summable.of_norm_bounded (summable_geometric_of_lt_one (norm_nonneg z) hz) (fun n => ?_)
  rw [norm_div, norm_pow, norm_natp1, div_le_iff₀ (by positivity)]
  have : (1 : ℝ) ≤ (n : ℝ) + 1 := by linarith [(Nat.cast_nonneg n : (0:ℝ) ≤ n)]
  nlinarith [pow_nonneg (norm_nonneg z) n]

lemma summable_Lc {z : ℂ} (hz : ‖z‖ ≤ 1) :
    Summable (fun n : ℕ => z ^ (n + 1) / ((n : ℂ) + 1) ^ 2) := by
  refine Summable.of_norm_bounded D00.summable_inv_sq (fun n => ?_)
  rw [norm_div, norm_pow, norm_pow, norm_natp1]
  apply div_le_div_of_nonneg_right _ (by positivity)
  exact pow_le_one₀ (norm_nonneg z) hz

lemma mul_Sc {z : ℂ} (hz : ‖z‖ < 1) : z * Sc z = -Complex.log (1 - z) := by
  have h1 := (summable_Sc hz).hasSum.mul_left z
  have h2 := Complex.hasSum_taylorSeries_neg_log' hz
  refine h1.unique (h2.congr_fun fun n => ?_)
  rw [pow_succ]; ring

lemma hasDerivAt_term (n : ℕ) (y : ℂ) :
    HasDerivAt (fun z : ℂ => z ^ (n + 1) / ((n : ℂ) + 1) ^ 2) (y ^ n / ((n : ℂ) + 1)) y := by
  have h := (hasDerivAt_pow (n + 1) y).div_const (((n : ℂ) + 1) ^ 2)
  convert h using 1
  have : ((n : ℂ) + 1) ≠ 0 := by exact_mod_cast Nat.succ_ne_zero n
  rw [Nat.add_sub_cancel, Nat.cast_add, Nat.cast_one, pow_two, mul_div_mul_comm, div_self this,
    one_mul]

lemma hasDerivAt_Lc {z : ℂ} (hz : ‖z‖ < 1) : HasDerivAt Lc (Sc z) z := by
  set r := (‖z‖ + 1) / 2 with hr
  have hr1 : r < 1 := by rw [hr]; linarith
  have hr0 : 0 ≤ r := by rw [hr]; positivity
  have hzr : z ∈ Metric.ball (0 : ℂ) r := by
    rw [Metric.mem_ball, dist_zero_right, hr]; linarith
  have h0r : (0 : ℂ) ∈ Metric.ball (0 : ℂ) r := Metric.mem_ball_self (by rw [hr]; positivity)
  have hb : ∀ (n : ℕ) (y : ℂ), y ∈ Metric.ball (0 : ℂ) r → ‖y ^ n / ((n : ℂ) + 1)‖ ≤ r ^ n := by
    intro n y hy
    rw [Metric.mem_ball, dist_zero_right] at hy
    rw [norm_div, norm_pow, norm_natp1, div_le_iff₀ (by positivity)]
    have : (1 : ℝ) ≤ (n : ℝ) + 1 := by linarith [(Nat.cast_nonneg n : (0:ℝ) ≤ n)]
    have := pow_le_pow_left₀ (norm_nonneg y) hy.le n
    nlinarith [pow_nonneg hr0 n]
  have hs0 : Summable (fun n : ℕ => (0 : ℂ) ^ (n + 1) / ((n : ℂ) + 1) ^ 2) := by
    simp only [zero_pow (Nat.succ_ne_zero _), zero_div]; exact summable_zero
  exact hasDerivAt_tsum_of_isPreconnected (g := fun (n : ℕ) (z : ℂ) => z ^ (n + 1) / ((n : ℂ) + 1) ^ 2)
    (g' := fun (n : ℕ) (z : ℂ) => z ^ n / ((n : ℂ) + 1)) (u := fun n : ℕ => r ^ n)
    (summable_geometric_of_lt_one hr0 hr1) Metric.isOpen_ball (convex_ball 0 r).isPreconnected
    (fun n y _ => hasDerivAt_term n y) hb h0r hs0 hzr

lemma Lc_ofReal {x : ℝ} (hx : |x| ≤ 1) : Lc (x : ℂ) = ((D00.Li2 x : ℝ) : ℂ) := by
  have h1 := Complex.ofRealCLM.hasSum (D00.proof_6 x hx)
  have hz : ‖(x : ℂ)‖ ≤ 1 := by rw [Complex.norm_real, Real.norm_eq_abs]; exact hx
  refine (summable_Lc hz).hasSum.unique (h1.congr_fun fun n => ?_)
  simp only [Complex.ofRealCLM_apply]; push_cast; ring

lemma slit_of_re {z : ℂ} (h : 0 < z.re) : z ∈ Complex.slitPlane :=
  Complex.mem_slitPlane_iff.2 (Or.inl h)

/-- The reflection function `Li₂(z) + Li₂(1-z) + log z log(1-z)`. -/
noncomputable def hRefl (z : ℂ) : ℂ := Lc z + Lc (1 - z) + Complex.log z * Complex.log (1 - z)

lemma ball_facts {z : ℂ} (hz : z ∈ Metric.ball ((1 : ℂ) / 2) (2 / 5)) :
    ‖z‖ < 1 ∧ ‖1 - z‖ < 1 ∧ 0 < z.re ∧ 0 < (1 - z).re := by
  rw [Metric.mem_ball, dist_eq_norm] at hz
  have hre : |(z - 1 / 2).re| ≤ ‖z - 1 / 2‖ := Complex.abs_re_le_norm _
  have e1 : ‖z‖ ≤ ‖z - 1 / 2‖ + 1 / 2 := by
    have := norm_add_le (z - 1 / 2) (1 / 2 : ℂ)
    rw [sub_add_cancel] at this
    have h2 : ‖(1 / 2 : ℂ)‖ = 1 / 2 := by norm_num
    linarith
  have e2 : ‖1 - z‖ ≤ ‖z - 1 / 2‖ + 1 / 2 := by
    have := norm_add_le (-(z - 1 / 2)) (1 / 2 : ℂ)
    rw [norm_neg] at this
    have h2 : ‖(1 / 2 : ℂ)‖ = 1 / 2 := by norm_num
    have h3 : -(z - 1 / 2) + 1 / 2 = 1 - z := by ring
    rw [h3] at this
    linarith
  have hr : (z - 1 / 2).re = z.re - 1 / 2 := by simp
  rw [hr, abs_le] at hre
  refine ⟨by linarith, by linarith, by linarith [hre.1], ?_⟩
  simp only [Complex.sub_re, Complex.one_re]; linarith [hre.2]

lemma hasDerivAt_hRefl {z : ℂ} (hz : z ∈ Metric.ball ((1 : ℂ) / 2) (2 / 5)) :
    HasDerivAt hRefl 0 z := by
  obtain ⟨n1, n2, r1, r2⟩ := ball_facts hz
  have d1 := hasDerivAt_Lc n1
  have d2 := (hasDerivAt_Lc n2).comp z ((hasDerivAt_id z).const_sub 1)
  have d3 := (hasDerivAt_id z).clog (slit_of_re r1)
  have d4 := ((hasDerivAt_id z).const_sub 1).clog (slit_of_re r2)
  have h := (d1.add d2).add (d3.mul d4)
  convert h using 1
  · funext w; simp only [hRefl, Pi.add_apply, Pi.mul_apply, Function.comp, id]
  have z0 : z ≠ 0 := by intro h0; rw [h0] at r1; simp at r1
  have z1 : (1 : ℂ) - z ≠ 0 := by intro h0; rw [h0] at r2; simp at r2
  have m1 := mul_Sc n1
  have m2 := mul_Sc n2
  rw [sub_sub_cancel] at m2
  have s1 : Sc z = -Complex.log (1 - z) / z := by field_simp; linear_combination m1
  have s2 : Sc (1 - z) = -Complex.log z / (1 - z) := by field_simp; linear_combination m2
  simp only [id]
  rw [s1, s2]
  field_simp
  ring

lemma hRefl_const {z : ℂ} (hz : z ∈ Metric.ball ((1 : ℂ) / 2) (2 / 5)) :
    hRefl z = hRefl (1 / 2) := by
  have hc : (1 / 2 : ℂ) ∈ Metric.ball ((1 : ℂ) / 2) (2 / 5) := Metric.mem_ball_self (by norm_num)
  refine Metric.isOpen_ball.is_const_of_deriv_eq_zero (convex_ball _ _).isPreconnected
    (fun w hw => (hasDerivAt_hRefl hw).differentiableAt.differentiableWithinAt)
    (fun w hw => (hasDerivAt_hRefl hw).deriv) hz hc

lemma hRefl_half : hRefl (1 / 2) = ((π ^ 2 / 6 : ℝ) : ℂ) := by
  unfold hRefl
  have e : (1 : ℂ) - 1 / 2 = ((1 / 2 : ℝ) : ℂ) := by push_cast; ring
  have e2 : (1 / 2 : ℂ) = ((1 / 2 : ℝ) : ℂ) := by push_cast; ring
  rw [e, e2, Lc_ofReal (by rw [abs_of_pos (by norm_num)]; norm_num), ← Complex.ofReal_log (by norm_num),
    D00.proof_14, one_div, Real.log_inv]
  push_cast
  ring

/-- `Re Li₂(ρ e^{iθ}) = Σ ρ^{n+1} cos((n+1)θ)/(n+1)²` for `0 ≤ ρ ≤ 1`. -/
lemma re_Lc_polar {ρ θ : ℝ} (h0 : 0 ≤ ρ) (h1 : ρ ≤ 1) :
    HasSum (fun n : ℕ => ρ ^ (n + 1) * Real.cos ((n + 1 : ℝ) * θ) / ((n : ℝ) + 1) ^ 2)
      (Lc ((ρ : ℂ) * (Complex.cos θ + Complex.sin θ * Complex.I))).re := by
  set z : ℂ := (ρ : ℂ) * (Complex.cos θ + Complex.sin θ * Complex.I) with hz
  have hzn : ‖z‖ ≤ 1 := by
    rw [hz, norm_mul, Complex.norm_cos_add_sin_mul_I, mul_one, Complex.norm_real,
      Real.norm_of_nonneg h0]
    exact h1
  have hS := Complex.reCLM.hasSum (summable_Lc hzn).hasSum
  have hterm : ∀ n : ℕ, z ^ (n + 1) / ((n : ℂ) + 1) ^ 2 =
      ((ρ ^ (n + 1) * Real.cos ((n + 1 : ℝ) * θ) / ((n : ℝ) + 1) ^ 2 : ℝ) : ℂ)
        + ((ρ ^ (n + 1) * Real.sin ((n + 1 : ℝ) * θ) / ((n : ℝ) + 1) ^ 2 : ℝ) : ℂ) * Complex.I := by
    intro n
    rw [hz, mul_pow, Complex.cos_add_sin_mul_I_pow]
    push_cast
    have : ((n : ℂ) + 1) ≠ 0 := by exact_mod_cast (Nat.succ_ne_zero n)
    field_simp
  unfold Lc
  refine hS.congr_fun (fun n => ?_)
  simp only [Complex.reCLM_apply]
  rw [hterm n]
  simp only [Complex.add_re, Complex.ofReal_re, Complex.mul_re, Complex.ofReal_im, Complex.I_re,
    Complex.I_im]
  ring

/-! ## Constants `cos(π/8)`, `sin(π/8)` -/

noncomputable def c8 : ℝ := Real.cos (π / 8)
noncomputable def s8 : ℝ := Real.sin (π / 8)

lemma sqrt2_sq : Real.sqrt 2 ^ 2 = 2 := Real.sq_sqrt (by norm_num)

lemma s8_pos : 0 < s8 := by
  unfold s8; apply Real.sin_pos_of_pos_of_lt_pi <;> linarith [Real.pi_pos]

lemma c8_pos : 0 < c8 := by
  unfold c8; apply Real.cos_pos_of_mem_Ioo; constructor <;> linarith [Real.pi_pos]

lemma s8_sq : s8 ^ 2 = (2 - Real.sqrt 2) / 4 := by
  unfold s8; rw [Real.sin_pi_div_eight, div_pow, Real.sq_sqrt]; · norm_num
  nlinarith [Real.sqrt_nonneg 2, sqrt2_sq, show Real.sqrt 2 < 2 by
    rw [Real.sqrt_lt' (by norm_num)]; norm_num]

lemma c8_sq : c8 ^ 2 = (2 + Real.sqrt 2) / 4 := by
  unfold c8; rw [Real.cos_pi_div_eight, div_pow, Real.sq_sqrt (by positivity)]; norm_num

lemma s8c8 : s8 * c8 = Real.sqrt 2 / 4 := by
  have h1 : (s8 * c8) ^ 2 = (Real.sqrt 2 / 4) ^ 2 := by
    rw [mul_pow, s8_sq, c8_sq, div_pow, sqrt2_sq]; ring_nf; rw [sqrt2_sq]; ring
  have := mul_pos s8_pos c8_pos
  nlinarith [Real.sqrt_nonneg 2, sq_nonneg (s8 * c8 - Real.sqrt 2 / 4),
    sq_nonneg (s8 * c8 + Real.sqrt 2 / 4)]

/-! ## `T(v) = Re Li₂((a-v)/(a+i))` -/

/-- `(a-v)/(a+i) = (a-v) sin(π/8) e^{-iπ/8}`. -/
noncomputable def zT (v : ℝ) : ℂ := (((D02.a - v) * s8 : ℝ) : ℂ) * ((c8 : ℂ) - (s8 : ℂ) * Complex.I)

noncomputable def T (v : ℝ) : ℝ := (Lc (zT v)).re

lemma norm_omega : ‖((c8 : ℂ) - (s8 : ℂ) * Complex.I)‖ = 1 := by
  have h : ((c8 : ℂ) - (s8 : ℂ) * Complex.I) = Complex.cos ((-(π / 8) : ℝ) : ℂ)
      + Complex.sin ((-(π / 8) : ℝ) : ℂ) * Complex.I := by
    rw [← Complex.ofReal_cos, ← Complex.ofReal_sin, Real.cos_neg, Real.sin_neg]
    unfold c8 s8; push_cast; ring
  rw [h, Complex.norm_cos_add_sin_mul_I]

lemma norm_zT (v : ℝ) : ‖zT v‖ = |D02.a - v| * s8 := by
  unfold zT
  rw [norm_mul, norm_omega, mul_one, Complex.norm_real, Real.norm_eq_abs, abs_mul,
    abs_of_pos s8_pos]

lemma two_s8_lt : 2 * s8 < 1 := by
  have h := s8_sq
  have hs := D02.hs1
  by_contra hc
  push Not at hc
  have : 1 ≤ (2 * s8) ^ 2 := by nlinarith
  rw [mul_pow, h] at this
  linarith

lemma norm_zT_lt {v : ℝ} (hv : |D02.a - v| ≤ 2) : ‖zT v‖ < 1 := by
  rw [norm_zT]
  nlinarith [two_s8_lt, s8_pos, abs_nonneg (D02.a - v)]

lemma c8s8_one : c8 ^ 2 + s8 ^ 2 = 1 := by unfold c8 s8; rw [add_comm]; exact Real.sin_sq_add_cos_sq _

lemma normSq_one_sub_zT (v : ℝ) : ‖1 - zT v‖ ^ 2 = s8 ^ 2 * (1 + v ^ 2) := by
  have hre : (1 - zT v).re = 1 - (D02.a - v) * s8 * c8 := by unfold zT; simp
  have him : (1 - zT v).im = (D02.a - v) * s8 * s8 := by unfold zT; simp
  rw [Complex.sq_norm, Complex.normSq_apply, hre, him]
  have ha : D02.a = 1 + Real.sqrt 2 := rfl
  calc (1 - (D02.a - v) * s8 * c8) * (1 - (D02.a - v) * s8 * c8)
        + (D02.a - v) * s8 * s8 * ((D02.a - v) * s8 * s8)
      = 1 - 2 * (D02.a - v) * (s8 * c8) + (D02.a - v) ^ 2 * s8 ^ 2 * (c8 ^ 2 + s8 ^ 2) := by ring
    _ = 1 - 2 * (1 + √2 - v) * (√2 / 4) + (1 + √2 - v) ^ 2 * ((2 - √2) / 4) * 1 := by
        rw [c8s8_one, s8c8, s8_sq, ha]
    _ = s8 ^ 2 * (1 + v ^ 2) := by
        rw [s8_sq]; linear_combination (-√2 / 4 + v / 2 - 1 / 2) * sqrt2_sq

/-- `ω = e^{-iπ/8}`. -/
noncomputable def om : ℂ := (c8 : ℂ) - (s8 : ℂ) * Complex.I

lemma zT_eq (v : ℝ) : zT v = (((D02.a - v : ℝ) : ℂ) * (s8 : ℂ)) * om := by
  unfold zT om; push_cast; ring

lemma log_norm_one_sub_zT (v : ℝ) :
    Real.log ‖1 - zT v‖ = Real.log (s8 ^ 2 * (1 + v ^ 2)) / 2 := by
  rw [← normSq_one_sub_zT, Real.log_pow]; push_cast; ring

lemma hasDerivAt_T {v : ℝ} (hv : |D02.a - v| ≤ 2) (hva : v ≠ D02.a) :
    HasDerivAt T (Real.log (s8 ^ 2 * (1 + v ^ 2)) / (2 * (D02.a - v))) v := by
  have hz := norm_zT_lt hv
  have hin : HasDerivAt (fun w : ℂ => ((D02.a : ℂ) - w) * (s8 : ℂ) * om) (-1 * (s8 : ℂ) * om) (v : ℂ) :=
    (((hasDerivAt_id (v : ℂ)).const_sub (D02.a : ℂ)).mul_const (s8 : ℂ)).mul_const om
  have e1 : ((D02.a : ℂ) - (v : ℂ)) * (s8 : ℂ) * om = zT v := by rw [zT_eq]; push_cast; ring
  have hL := hasDerivAt_Lc hz
  rw [← e1] at hL
  have hc := (hL.comp (v : ℂ) hin).real_of_complex
  convert hc using 1
  · funext x; simp only [T, Function.comp]; rw [zT_eq]; push_cast; ring_nf
  · rw [e1]
    have hav : ((D02.a - v : ℝ) : ℂ) ≠ 0 := by
      exact_mod_cast sub_ne_zero.mpr (Ne.symm hva)
    have hm := mul_Sc hz
    have key : Sc (zT v) * (-1 * (s8 : ℂ) * om) = Complex.log (1 - zT v) / ((D02.a - v : ℝ) : ℂ) := by
      rw [eq_div_iff hav, ← neg_neg (Complex.log (1 - zT v)), ← hm, zT_eq]; ring
    rw [key, Complex.div_ofReal_re, Complex.log_re, log_norm_one_sub_zT]
    field_simp

lemma continuousAt_T {v : ℝ} (hv : |D02.a - v| ≤ 2) : ContinuousAt T v := by
  have hz := norm_zT_lt hv
  have hzc : Continuous zT := by unfold zT; fun_prop
  have h1 : ContinuousAt (fun x => Lc (zT x)) v :=
    ContinuousAt.comp (f := zT) (g := Lc) (hasDerivAt_Lc hz).continuousAt hzc.continuousAt
  exact Complex.continuous_re.continuousAt.comp h1

lemma T_a : T D02.a = 0 := by
  unfold T; rw [zT_eq, sub_self]; simp [Lc]

lemma b_def : D02.b = Real.sqrt 2 - 1 := rfl
lemma a_def : D02.a = 1 + Real.sqrt 2 := rfl

lemma zTb_re : (zT D02.b).re = Real.sqrt 2 / 2 := by
  unfold zT; simp only [Complex.mul_re, Complex.ofReal_re, Complex.ofReal_im, Complex.sub_re,
    Complex.sub_im, Complex.mul_im, Complex.I_re, Complex.I_im]
  rw [a_def, b_def]; linear_combination 2 * s8c8

lemma zTb_im : (zT D02.b).im = -(2 - Real.sqrt 2) / 2 := by
  unfold zT; simp only [Complex.mul_re, Complex.ofReal_re, Complex.ofReal_im, Complex.sub_re,
    Complex.sub_im, Complex.mul_im, Complex.I_re, Complex.I_im]
  rw [a_def, b_def]; linear_combination (-2) * s8_sq

lemma zTb_ball : zT D02.b ∈ Metric.ball ((1 : ℂ) / 2) (2 / 5) := by
  rw [Metric.mem_ball, dist_eq_norm]
  have h14 : (1.4 : ℝ) < Real.sqrt 2 := by rw [Real.lt_sqrt (by norm_num)]; norm_num
  have hs2 : Real.sqrt 2 < 1.5 := by rw [Real.sqrt_lt' (by norm_num)]; norm_num
  have hsq : ‖zT D02.b - 1 / 2‖ ^ 2 < (2 / 5) ^ 2 := by
    rw [Complex.sq_norm, Complex.normSq_apply]
    have r1 : (zT D02.b - 1 / 2).re = Real.sqrt 2 / 2 - 1 / 2 := by
      rw [Complex.sub_re, zTb_re]; norm_num
    have r2 : (zT D02.b - 1 / 2).im = -(2 - Real.sqrt 2) / 2 := by
      rw [Complex.sub_im, zTb_im]; norm_num
    rw [r1, r2]
    nlinarith [sqrt2_sq]
  exact abs_lt_of_sq_lt_sq' hsq (by norm_num) |>.2

lemma one_sub_zTb : 1 - zT D02.b = (D02.b : ℂ) * (Complex.cos ((π / 4 : ℝ) : ℂ)
    + Complex.sin ((π / 4 : ℝ) : ℂ) * Complex.I) := by
  rw [← Complex.ofReal_cos, ← Complex.ofReal_sin, Real.cos_pi_div_four, Real.sin_pi_div_four]
  apply Complex.ext
  · rw [Complex.sub_re, zTb_re]; simp only [Complex.one_re, Complex.mul_re, Complex.ofReal_re,
      Complex.ofReal_im, Complex.add_re, Complex.add_im, Complex.I_re, Complex.I_im]
    rw [b_def]; linear_combination (-1 / 2) * sqrt2_sq
  · rw [Complex.sub_im, zTb_im]; simp only [Complex.one_im, Complex.mul_im, Complex.ofReal_re,
      Complex.ofReal_im, Complex.add_re, Complex.add_im, Complex.I_re, Complex.I_im]
    rw [b_def]; linear_combination (-1 / 2) * sqrt2_sq

lemma re_Lc_one_sub_zTb : (Lc (1 - zT D02.b)).re = (1 / 2) * D00.Li2 D02.b + π ^ 2 / 192 := by
  rw [one_sub_zTb]
  exact (re_Lc_polar D02.b_pos.le D02.b_lt.le).unique S02.lemmaA_pi4

lemma zTb_polar : zT D02.b = ((2 * s8 : ℝ) : ℂ) * (Complex.cos ((-(π / 8) : ℝ) : ℂ)
    + Complex.sin ((-(π / 8) : ℝ) : ℂ) * Complex.I) := by
  rw [zT_eq, ← Complex.ofReal_cos, ← Complex.ofReal_sin, Real.cos_neg, Real.sin_neg]
  unfold om c8 s8
  rw [a_def, b_def]; push_cast; ring

lemma log_zTb : Complex.log (zT D02.b) = ((Real.log (2 * s8) : ℝ) : ℂ) + ((-(π / 8) : ℝ) : ℂ) * Complex.I := by
  have hr : 0 < 2 * s8 := by linarith [s8_pos]
  have hθ : -(π / 8) ∈ Set.Ioc (-π) π := by constructor <;> linarith [Real.pi_pos]
  apply Complex.ext
  · rw [Complex.log_re, zTb_polar, norm_mul, Complex.norm_cos_add_sin_mul_I, mul_one, Complex.norm_real,
      Real.norm_of_nonneg hr.le]; simp
  · rw [Complex.log_im, zTb_polar, Complex.arg_mul_cos_add_sin_mul_I hr hθ]; simp

lemma log_one_sub_zTb : Complex.log (1 - zT D02.b) = ((Real.log D02.b : ℝ) : ℂ) + ((π / 4 : ℝ) : ℂ) * Complex.I := by
  have hr := D02.b_pos
  have hθ : π / 4 ∈ Set.Ioc (-π) π := by constructor <;> linarith [Real.pi_pos]
  apply Complex.ext
  · rw [Complex.log_re, one_sub_zTb, norm_mul, Complex.norm_cos_add_sin_mul_I, mul_one, Complex.norm_real,
      Real.norm_of_nonneg hr.le]; simp
  · rw [Complex.log_im, one_sub_zTb, Complex.arg_mul_cos_add_sin_mul_I hr hθ]; simp

lemma log_two_s8 : Real.log (2 * s8) = D02.lam / 4 - D02.L / 2 := by
  have h1 : (2 * s8) ^ 2 = Real.sqrt 2 * D02.b := by
    rw [mul_pow, s8_sq, b_def]; linear_combination (-1) * sqrt2_sq
  have h2 : Real.log ((2 * s8) ^ 2) = D02.lam / 2 - D02.L := by
    rw [h1, Real.log_mul (by positivity) D02.b_pos.ne', D02.log_s, D02.log_b]; ring
  rw [Real.log_pow] at h2; push_cast at h2; linarith

lemma T_b : T D02.b = 25 * π ^ 2 / 192 - D00.Li2 D02.b / 2 + D02.lam * D02.L / 4 - D02.L ^ 2 / 2 := by
  have h := congrArg Complex.re ((hRefl_const zTb_ball).trans hRefl_half)
  unfold hRefl at h
  rw [Complex.add_re, Complex.add_re, re_Lc_one_sub_zTb, log_zTb, log_one_sub_zTb] at h
  simp only [Complex.mul_re, Complex.add_re, Complex.add_im, Complex.ofReal_re, Complex.ofReal_im,
    Complex.mul_im, Complex.I_re, Complex.I_im] at h
  rw [log_two_s8, D02.log_b] at h
  unfold T
  linarith

/-! ## Real primitives -/

lemma dL1 {c p v : ℝ} (hcp : 0 < c + p) (hvp : 0 < v + p) (hvc : v ≠ c) :
    HasDerivAt (fun v => D00.Li2 ((c - v) / (c + p)))
      ((Real.log (v + p) - Real.log (c + p)) / (c - v)) v := by
  have hf : HasDerivAt (fun v => (c - v) / (c + p)) (-1 / (c + p)) v := by
    simpa using ((hasDerivAt_id v).const_sub c).div_const (c + p)
  have h1 : (c - v) / (c + p) ≠ 1 := by
    rw [Ne, div_eq_one_iff_eq hcp.ne']; intro h; linarith
  have h0 : (c - v) / (c + p) ≠ 0 := div_ne_zero (sub_ne_zero.mpr (Ne.symm hvc)) hcp.ne'
  convert D00.Li2_comp hf h1 h0 using 1
  have e : 1 - (c - v) / (c + p) = (v + p) / (c + p) := by field_simp; ring
  rw [e, Real.log_div hvp.ne' hcp.ne']
  have : c - v ≠ 0 := sub_ne_zero.mpr (Ne.symm hvc)
  field_simp

lemma dL3 {w e v : ℝ} (he : 0 < e) (hev : e < v + w) :
    HasDerivAt (fun v => D00.Li2 (e / (v + w)))
      ((Real.log (v + w - e) - Real.log (v + w)) / (v + w)) v := by
  have hvw : 0 < v + w := by linarith
  have hf : HasDerivAt (fun v => e / (v + w)) (-e / (v + w) ^ 2) v := by
    have := ((hasDerivAt_id v).add_const w).inv hvw.ne'
    convert this.const_mul e using 1
    · funext x; simp [div_eq_mul_inv]
    · simp only [id]; field_simp
  have h1 : e / (v + w) ≠ 1 := by
    rw [Ne, div_eq_one_iff_eq hvw.ne']; intro h; linarith
  have h0 : e / (v + w) ≠ 0 := div_ne_zero he.ne' hvw.ne'
  convert D00.Li2_comp hf h1 h0 using 1
  have e1 : 1 - e / (v + w) = (v + w - e) / (v + w) := by field_simp
  rw [e1, Real.log_div (by linarith) hvw.ne']
  field_simp

lemma dL4 {v : ℝ} (hv : v ≠ 0) :
    HasDerivAt (fun v => D00.Li2 (-v ^ 2)) (-2 * Real.log (v ^ 2 + 1) / v) v := by
  have hf : HasDerivAt (fun v : ℝ => -v ^ 2) (-(2 * v)) v := by
    have := (hasDerivAt_pow 2 v).neg
    convert this using 1; simp
  have hv2 : 0 < v ^ 2 := by positivity
  have h1 : -v ^ 2 ≠ 1 := by intro h; linarith
  have h0 : -v ^ 2 ≠ 0 := by intro h; linarith
  convert D00.Li2_comp hf h1 h0 using 1
  rw [sub_neg_eq_add, add_comm]
  field_simp

set_option maxRecDepth 100000 in
set_option maxHeartbeats 4000000 in
lemma key_identity (t v A B C D l L : ℝ) (ht : t ^ 2 = 2) (hv : 0 < v) (h1 : 1 + t - v ≠ 0)
    (h2 : t - 1 - v ≠ 0) (h3 : (1 + t) * v - 1 ≠ 0) (h4 : v + (1 + t) ≠ 0) (h5 : v + (t - 1) ≠ 0)
    (h6 : -(t - 1) - v ≠ 0) (h7 : v - (1 + t) ≠ 0) (h8 : v - (t - 1) ≠ 0) (h9 : t ≠ 0) :
    (A - (l + L)) / (1 + t - v) + (B - 3 * l / 2) / (1 + t - v)
      - 2 * ((D - (3 * l / 2 + L)) / (2 * (1 + t - v)))
      - (A - 3 * l / 2) / (t - 1 - v) - (B - (l - L)) / (t - 1 - v)
      + 2 * (-(D - 2 * C - (3 * l / 2 + L)) / (2 * v * ((1 + t) * v - 1)))
      + 2 * ((C + L) / (t - 1 - v)) + (1 / 2 - t / 4) * (-2 * D / v) + (1 - t / 4) * (2 * C / v)
      + (3 * l / 2 + L - t * l / 4) / v + A / (v + (1 + t)) + (B - A) / (v + (1 + t))
      - (C - A) / (v + (1 + t)) - 3 / 2 * l / (v + (1 + t)) + (A - l) / (-(t - 1) - v)
      + (C - B) / (v + (t - 1)) + l / 2 / (v + (t - 1))
    = -2 * (v ^ 2 + 1) / ((v - (1 + t)) * (v - (t - 1))) *
        (2 * (A + B - C - 3 * l / 2) / ((v + (1 + t)) * (v + (t - 1))) - (D - C - l / 2) / (2 * t * v)) := by
  have hv0 : v ≠ 0 := hv.ne'
  have t3 : t ^ 3 = 2 * t := by rw [pow_succ, ht]
  have t4 : t ^ 4 = 4 := by rw [show t ^ 4 = (t ^ 2) ^ 2 by ring, ht]; norm_num
  have t5 : t ^ 5 = 4 * t := by rw [pow_succ, t4]
  have t6 : t ^ 6 = 8 := by rw [show t ^ 6 = (t ^ 2) ^ 3 by ring, ht]; norm_num
  have t7 : t ^ 7 = 8 * t := by rw [pow_succ, t6]
  have t8 : t ^ 8 = 16 := by rw [show t ^ 8 = (t ^ 2) ^ 4 by ring, ht]; norm_num
  have t9 : t ^ 9 = 16 * t := by rw [pow_succ, t8]
  have t10 : t ^ 10 = 32 := by rw [show t ^ 10 = (t ^ 2) ^ 5 by ring, ht]; norm_num
  field_simp
  ring_nf
  simp only [ht, t3, t4, t5, t6, t7, t8, t9, t10]
  ring_nf

/-! ## The primitive `G` on `[b, a]` -/

/-- `K = log(a² + 1) = 3λ/2 + L`. -/
noncomputable def K : ℝ := 3 * D02.lam / 2 + D02.L

lemma sqrt2_pos : 0 < Real.sqrt 2 := by positivity

lemma log_a_a : Real.log (D02.a + D02.a) = D02.lam + D02.L := by
  rw [show D02.a + D02.a = 2 * D02.a by ring, Real.log_mul (by norm_num) D02.a_pos.ne']; rfl

lemma log_a_b : Real.log (D02.a + D02.b) = 3 * D02.lam / 2 := by
  rw [show D02.a + D02.b = 2 * Real.sqrt 2 by rw [a_def, b_def]; ring,
    Real.log_mul (by norm_num) sqrt2_pos.ne', D02.log_s]; unfold D02.lam; ring

lemma log_b_a : Real.log (D02.b + D02.a) = 3 * D02.lam / 2 := by rw [add_comm]; exact log_a_b

lemma log_b_b : Real.log (D02.b + D02.b) = D02.lam - D02.L := by
  rw [show D02.b + D02.b = 2 * D02.b by ring, Real.log_mul (by norm_num) D02.b_pos.ne', D02.log_b]
  unfold D02.lam; ring

lemma log_b_0 : Real.log (D02.b + 0) = -D02.L := by rw [add_zero, D02.log_b]

lemma log_mb_a : Real.log (-D02.b + D02.a) = D02.lam := by
  rw [show -D02.b + D02.a = 2 by rw [a_def, b_def]; ring]; rfl

lemma s8_sq_eq : s8 ^ 2 * (2 * Real.sqrt 2 * D02.a) = 1 := by
  rw [s8_sq, a_def]; linear_combination (1 / 2 - Real.sqrt 2 / 2) * sqrt2_sq

lemma log_s8_sq : Real.log (s8 ^ 2) = -K := by
  have h : s8 ^ 2 = (2 * Real.sqrt 2 * D02.a)⁻¹ := by
    exact eq_inv_of_mul_eq_one_left s8_sq_eq
  rw [h, Real.log_inv, Real.log_mul (by positivity) D02.a_pos.ne',
    Real.log_mul (by norm_num) sqrt2_pos.ne', D02.log_s]
  unfold K D02.lam D02.L; ring

lemma log_T_arg (v : ℝ) : Real.log (s8 ^ 2 * (1 + v ^ 2)) = Real.log (v ^ 2 + 1) - K := by
  rw [Real.log_mul (by positivity [s8_pos]) (by positivity), log_s8_sq,
    show (1 : ℝ) + v ^ 2 = v ^ 2 + 1 by ring]; ring

lemma Tinv_val {v : ℝ} (hv : 0 < v) (hav : D02.a * v - 1 ≠ 0) :
    Real.log (s8 ^ 2 * (1 + (v⁻¹) ^ 2)) / (2 * (D02.a - v⁻¹)) * -(v ^ 2)⁻¹
      = -(Real.log (v ^ 2 + 1) - 2 * Real.log v - K) / (2 * v * (D02.a * v - 1)) := by
  have e : (v⁻¹) ^ 2 + 1 = (v ^ 2 + 1) / v ^ 2 := by field_simp; ring
  rw [log_T_arg, e, Real.log_div (by positivity) (by positivity), Real.log_pow]
  have : D02.a - v⁻¹ = (D02.a * v - 1) / v := by field_simp
  rw [this]
  field_simp
  push_cast
  ring

lemma dlog {w v : ℝ} (h : 0 < v + w) : HasDerivAt (fun v => Real.log (v + w)) (1 / (v + w)) v :=
  ((hasDerivAt_id v).add_const w).log h.ne'

lemma dlog2 {w v : ℝ} (h : 0 < v + w) :
    HasDerivAt (fun v => Real.log (v + w) ^ 2) (2 * Real.log (v + w) / (v + w)) v := by
  have h2 := (dlog h).mul (dlog h)
  have e : (fun v => Real.log (v + w)) * (fun v => Real.log (v + w)) = fun v => Real.log (v + w) ^ 2 := by
    funext x; simp [pow_two]
  rw [e] at h2; convert h2 using 1; ring

lemma dlogv {v : ℝ} (h : 0 < v) : HasDerivAt (fun v => Real.log v) (1 / v) v := by
  rw [one_div]; exact Real.hasDerivAt_log h.ne'

lemma dlogv2 {v : ℝ} (h : 0 < v) : HasDerivAt (fun v => Real.log v ^ 2) (2 * Real.log v / v) v := by
  have h2 := (dlogv h).mul (dlogv h)
  have e : (fun v => Real.log v) * (fun v => Real.log v) = fun v => Real.log v ^ 2 := by
    funext x; simp [pow_two]
  rw [e] at h2; convert h2 using 1; ring

noncomputable def G (v : ℝ) : ℝ :=
  D00.Li2 ((D02.a - v) / (D02.a + D02.a)) + D00.Li2 ((D02.a - v) / (D02.a + D02.b)) - 2 * T v
    - D00.Li2 ((D02.b - v) / (D02.b + D02.a)) - D00.Li2 ((D02.b - v) / (D02.b + D02.b))
    + 2 * T v⁻¹ + 2 * D00.Li2 ((D02.b - v) / (D02.b + 0))
    + (1 / 2 - Real.sqrt 2 / 4) * D00.Li2 (-v ^ 2) + (1 - Real.sqrt 2 / 4) * Real.log v ^ 2
    + (3 * D02.lam / 2 + D02.L - Real.sqrt 2 * D02.lam / 4) * Real.log v
    + Real.log (v + D02.a) ^ 2 / 2 + D00.Li2 (2 / (v + D02.a)) - D00.Li2 (D02.a / (v + D02.a))
    - 3 / 2 * D02.lam * Real.log (v + D02.a) + D00.Li2 ((-D02.b - v) / (-D02.b + D02.a))
    + D00.Li2 (D02.b / (v + D02.b)) + D02.lam / 2 * Real.log (v + D02.b)

/-- The integrand in the variable `v`, in factored form. -/
noncomputable def gform (v : ℝ) : ℝ :=
  -2 * (v ^ 2 + 1) / ((v - D02.a) * (v - D02.b)) *
    (2 * (Real.log (v + D02.a) + Real.log (v + D02.b) - Real.log v - 3 * D02.lam / 2)
        / ((v + D02.a) * (v + D02.b))
      - (Real.log (v ^ 2 + 1) - Real.log v - D02.lam / 2) / (2 * Real.sqrt 2 * v))

lemma hasDerivAt_G {v : ℝ} (hv : v ∈ Ioo D02.b D02.a) : HasDerivAt G (gform v) v := by
  obtain ⟨hb, ha⟩ := hv
  have b0 := D02.b_pos
  have a0 := D02.a_pos
  have hab : D02.a = D02.b + 2 := by rw [a_def, b_def]; ring
  have v0 : 0 < v := by linarith
  have ab1 := D02.ab
  have h1 := dL1 (c := D02.a) (p := D02.a) (v := v) (by linarith) (by linarith) ha.ne
  have h2 := dL1 (c := D02.a) (p := D02.b) (v := v) (by linarith) (by linarith) ha.ne
  have h3 := hasDerivAt_T (v := v) (by rw [abs_le]; constructor <;> linarith) ha.ne
  have h4 := dL1 (c := D02.b) (p := D02.a) (v := v) (by linarith) (by linarith) hb.ne'
  have h5 := dL1 (c := D02.b) (p := D02.b) (v := v) (by linarith) (by linarith) hb.ne'
  have vi1 : D02.b < v⁻¹ := by
    rw [lt_inv_comm₀ b0 v0, show D02.b⁻¹ = D02.a by exact inv_eq_of_mul_eq_one_right (by rw [mul_comm]; exact ab1)]
    exact ha
  have vi2 : v⁻¹ < D02.a := by
    rw [inv_lt_comm₀ v0 a0, show D02.a⁻¹ = D02.b by exact inv_eq_of_mul_eq_one_right ab1]
    exact hb
  have h6 := (hasDerivAt_T (v := v⁻¹) (by rw [abs_le]; constructor <;> linarith) vi2.ne).comp v
    (hasDerivAt_inv v0.ne')
  have h7 := dL1 (c := D02.b) (p := 0) (v := v) (by linarith) (by linarith) hb.ne'
  have h8 := dL4 v0.ne'
  have h9 := dlogv2 v0
  have h10 := dlogv v0
  have h11 := dlog2 (w := D02.a) (v := v) (by linarith)
  have h12 := dL3 (w := D02.a) (e := 2) (v := v) (by norm_num) (by linarith)
  have h13 := dL3 (w := D02.a) (e := D02.a) (v := v) a0 (by linarith)
  have h14 := dlog (w := D02.a) (v := v) (by linarith)
  have h15 := dL1 (c := -D02.b) (p := D02.a) (v := v) (by linarith) (by linarith) (by linarith)
  have h16 := dL3 (w := D02.b) (e := D02.b) (v := v) b0 (by linarith)
  have h17 := dlog (w := D02.b) (v := v) (by linarith)
  have g1 := (h1.add h2).sub (h3.const_mul 2)
  have g2 := ((g1.sub h4).sub h5).add (h6.const_mul 2)
  have g3 := (g2.add (h7.const_mul 2)).add (h8.const_mul (1 / 2 - Real.sqrt 2 / 4))
  have g4 := (g3.add (h9.const_mul (1 - Real.sqrt 2 / 4))).add
    (h10.const_mul (3 * D02.lam / 2 + D02.L - Real.sqrt 2 * D02.lam / 4))
  have g5 := ((g4.add (h11.div_const 2)).add h12).sub h13
  have g6 := ((g5.sub (h14.const_mul (3 / 2 * D02.lam))).add h15).add h16
  have h := g6.add (h17.const_mul (D02.lam / 2))
  convert h using 1
  · funext x; rfl
  have hav : D02.a * v - 1 ≠ 0 := by
    have : D02.a * D02.b < D02.a * v := mul_lt_mul_of_pos_left hb a0
    exact (sub_pos.mpr (by linarith)).ne'
  rw [Tinv_val v0 hav, log_T_arg, log_a_a, log_a_b, log_b_a, log_b_b, log_b_0, log_mb_a,
    show v + D02.a - 2 = v + D02.b by rw [hab]; ring, show v + D02.a - D02.a = v by ring,
    show v + D02.b - D02.b = v by ring, add_zero]
  have hb' : Real.sqrt 2 - 1 < v := hb
  have ha' : v < 1 + Real.sqrt 2 := ha
  have hav' : (1 + Real.sqrt 2) * v - 1 ≠ 0 := hav
  have hs1 := D02.hs1
  have key := key_identity (Real.sqrt 2) v (Real.log (v + (1 + Real.sqrt 2)))
    (Real.log (v + (Real.sqrt 2 - 1))) (Real.log v) (Real.log (v ^ 2 + 1)) D02.lam D02.L sqrt2_sq v0
    (by linarith) (by linarith) hav' (by linarith) (by linarith) (by linarith) (by linarith)
    (by linarith) sqrt2_pos.ne'
  unfold gform K
  rw [a_def, b_def]
  linear_combination -key

lemma continuousAt_G {v : ℝ} (hv : v ∈ Icc D02.b D02.a) : ContinuousAt G v := by
  obtain ⟨hb, ha⟩ := hv
  have b0 := D02.b_pos
  have a0 := D02.a_pos
  have hab : D02.a = D02.b + 2 := by rw [a_def, b_def]; ring
  have v0 : 0 < v := by linarith
  have ab1 := D02.ab
  have hT1 : ContinuousAt T v := continuousAt_T (by rw [abs_le]; constructor <;> linarith)
  have vi1 : D02.b ≤ v⁻¹ := by
    rw [le_inv_comm₀ b0 v0, show D02.b⁻¹ = D02.a from
      inv_eq_of_mul_eq_one_right (by rw [mul_comm]; exact ab1)]
    exact ha
  have vi2 : v⁻¹ ≤ D02.a := by
    rw [inv_le_comm₀ v0 a0, show D02.a⁻¹ = D02.b from inv_eq_of_mul_eq_one_right ab1]
    exact hb
  have hT2 : ContinuousAt (fun x => T x⁻¹) v :=
    ContinuousAt.comp (g := T) (f := fun x => x⁻¹)
      (continuousAt_T (by rw [abs_le]; constructor <;> linarith)) (continuousAt_inv₀ v0.ne')
  have hL := D00.continuous_Li2
  have p1 : v + D02.a ≠ 0 := by linarith
  have p2 : v + D02.b ≠ 0 := by linarith
  have p3 : v ≠ 0 := v0.ne'
  have cL : ∀ f : ℝ → ℝ, ContinuousAt f v → ContinuousAt (fun x => D00.Li2 (f x)) v :=
    fun f hf => hL.continuousAt.comp (f := f) hf
  have lA : ContinuousAt (fun x => Real.log (x + D02.a)) v :=
    (Real.continuousAt_log p1).comp (f := fun x => x + D02.a) (continuousAt_id.add continuousAt_const)
  have lB : ContinuousAt (fun x => Real.log (x + D02.b)) v :=
    (Real.continuousAt_log p2).comp (f := fun x => x + D02.b) (continuousAt_id.add continuousAt_const)
  have lC : ContinuousAt (fun x => Real.log x) v := Real.continuousAt_log p3
  have t1 := cL (fun x => (D02.a - x) / (D02.a + D02.a)) (by fun_prop)
  have t2 := cL (fun x => (D02.a - x) / (D02.a + D02.b)) (by fun_prop)
  have t4 := cL (fun x => (D02.b - x) / (D02.b + D02.a)) (by fun_prop)
  have t5 := cL (fun x => (D02.b - x) / (D02.b + D02.b)) (by fun_prop)
  have t7 := cL (fun x => (D02.b - x) / (D02.b + 0)) (by fun_prop)
  have t8 := cL (fun x => -x ^ 2) (by fun_prop)
  have t12 := cL (fun x => 2 / (x + D02.a))
    (continuousAt_const.div (continuousAt_id.add continuousAt_const) p1)
  have t13 := cL (fun x => D02.a / (x + D02.a))
    (continuousAt_const.div (continuousAt_id.add continuousAt_const) p1)
  have t15 := cL (fun x => (-D02.b - x) / (-D02.b + D02.a)) (by fun_prop)
  have t16 := cL (fun x => D02.b / (x + D02.b))
    (continuousAt_const.div (continuousAt_id.add continuousAt_const) p2)
  have cm : ∀ (c : ℝ) (f : ℝ → ℝ), ContinuousAt f v → ContinuousAt (fun x => c * f x) v :=
    fun c f hf => continuousAt_const.mul hf
  have k := (((((((((((((((t1.add t2).sub (cm 2 _ hT1)).sub t4).sub t5).add
    (cm 2 _ hT2)).add (cm 2 _ t7)).add
    (cm (1 / 2 - Real.sqrt 2 / 4) _ t8)).add (cm (1 - Real.sqrt 2 / 4) _ (lC.pow 2))).add
    (cm (3 * D02.lam / 2 + D02.L - Real.sqrt 2 * D02.lam / 4) _ lC)).add ((lA.pow 2).div_const 2)).add
    t12).sub t13).sub (cm (3 / 2 * D02.lam) _ lA)).add t15).add t16).add (cm (D02.lam / 2) _ lB)
  exact k

lemma continuousOn_G : ContinuousOn G (Icc D02.b D02.a) :=
  fun _ hv => (continuousAt_G hv).continuousWithinAt

/-! ## The substitution `v = 2s - 1 + √2 q(s)` -/

noncomputable def Vs (s : ℝ) : ℝ := 2 * s - 1 + Real.sqrt 2 * S02.qs s

lemma qs_sq (s : ℝ) : S02.qs s ^ 2 = s ^ 2 + (1 - s) ^ 2 := Real.sq_sqrt (by positivity)

lemma hasDerivAt_qs (s : ℝ) : HasDerivAt S02.qs ((2 * s - 1) / S02.qs s) s := by
  have hf : HasDerivAt (fun s : ℝ => s ^ 2 + (1 - s) ^ 2) (2 * s - 2 * (1 - s)) s := by
    have h := ((hasDerivAt_id s).mul (hasDerivAt_id s)).add
      (((hasDerivAt_id s).const_sub 1).mul ((hasDerivAt_id s).const_sub 1))
    convert h using 1
    · funext x; simp only [Pi.add_apply, Pi.mul_apply, id]; ring
    · simp only [id]; ring
  have hq := S02.qs_pos s
  have h := hf.sqrt (by have := qs_sq s; nlinarith)
  unfold S02.qs at hq ⊢
  convert h using 1
  field_simp
  ring

lemma W_mul (s : ℝ) : Vs s * (Real.sqrt 2 * S02.qs s - (2 * s - 1)) = 1 := by
  unfold Vs
  linear_combination (S02.qs s ^ 2) * sqrt2_sq + 2 * qs_sq s

lemma Vs_quad (s : ℝ) : Vs s ^ 2 - (4 * s - 2) * Vs s - 1 = 0 := by
  have h := W_mul s
  have e : Real.sqrt 2 * S02.qs s - (2 * s - 1) = Vs s - (4 * s - 2) := by unfold Vs; ring
  rw [e] at h
  linear_combination h

lemma Vs_pos (s : ℝ) : 0 < Vs s := by
  have h1 : (Real.sqrt 2 * S02.qs s) ^ 2 = (2 * s - 1) ^ 2 + 1 := by
    rw [mul_pow, sqrt2_sq, qs_sq]; ring
  have h2 : 0 ≤ Real.sqrt 2 * S02.qs s := by have := S02.qs_pos s; positivity
  unfold Vs
  nlinarith [sq_nonneg (2 * s - 1 + Real.sqrt 2 * S02.qs s)]

lemma Vs_gt {s : ℝ} (hs : 0 < s) : D02.b < Vs s := by
  have hq := Vs_quad s
  have hp := Vs_pos s
  have hs1 := D02.hs1
  rw [b_def]
  by_contra hc
  push Not at hc
  have : (Vs s + 1) ^ 2 ≤ 2 := by
    have := sqrt2_sq
    nlinarith
  nlinarith

lemma Vs_lt {s : ℝ} (hs : s < 1) : Vs s < D02.a := by
  have hq := Vs_quad s
  have hp := Vs_pos s
  have hs1 := D02.hs1
  rw [a_def]
  by_contra hc
  push Not at hc
  have : 2 ≤ (Vs s - 1) ^ 2 := by
    have := sqrt2_sq
    nlinarith
  nlinarith

lemma Vs_ge {s : ℝ} (hs : 0 ≤ s) : D02.b ≤ Vs s := by
  have hq := Vs_quad s
  have hp := Vs_pos s
  have hs1 := D02.hs1
  rw [b_def]
  by_contra hc
  push Not at hc
  have : (Vs s + 1) ^ 2 < 2 := by
    have := sqrt2_sq
    nlinarith
  nlinarith

lemma Vs_le {s : ℝ} (hs : s ≤ 1) : Vs s ≤ D02.a := by
  have hq := Vs_quad s
  have hp := Vs_pos s
  have hs1 := D02.hs1
  rw [a_def]
  by_contra hc
  push Not at hc
  have : 2 < (Vs s - 1) ^ 2 := by
    have := sqrt2_sq
    nlinarith
  nlinarith

lemma Vs_zero : Vs 0 = D02.b := by
  unfold Vs S02.qs; rw [b_def]; norm_num; ring

lemma Vs_one : Vs 1 = D02.a := by
  unfold Vs S02.qs; rw [a_def]; norm_num

lemma continuous_Vs : Continuous Vs := by
  unfold Vs; exact (by fun_prop : Continuous fun s : ℝ => 2 * s - 1).add
    (continuous_const.mul S02.continuous_qs)

lemma hasDerivAt_Vs (s : ℝ) :
    HasDerivAt Vs (2 + Real.sqrt 2 * ((2 * s - 1) / S02.qs s)) s := by
  have h := ((((hasDerivAt_id s).const_mul 2).sub_const 1).add ((hasDerivAt_qs s).const_mul (Real.sqrt 2)))
  convert h using 1
  · funext x; rfl
  · ring

lemma gform_Vs {s : ℝ} (hs : s ∈ Ioo (0 : ℝ) 1) :
    gform (Vs s) * (2 + Real.sqrt 2 * ((2 * s - 1) / S02.qs s)) = S02.f1 (S02.qs s) := by
  obtain ⟨h0, h1⟩ := hs
  have hq0 := S02.qs_pos s
  have hq1 := S02.qs_lt_one h0 h1
  have hV := Vs_pos s
  have hW := W_mul s
  have t2 := sqrt2_sq
  have st := sqrt2_pos
  set q := S02.qs s with hqdef
  set V := Vs s with hVdef
  have hVd : V = 2 * s - 1 + Real.sqrt 2 * q := rfl
  have e2 : V ^ 2 + 1 = 2 * Real.sqrt 2 * q * V := by linear_combination V * hVd - hW
  have e1 : (V + D02.a) * (V + D02.b) = 2 * Real.sqrt 2 * V * (1 + q) := by
    rw [a_def, b_def]; linear_combination e2 + t2
  have e3 : (V - D02.a) * (V - D02.b) = 2 * Real.sqrt 2 * V * (q - 1) := by
    rw [a_def, b_def]; linear_combination e2 + t2
  have hl22 : Real.log (2 * Real.sqrt 2) = 3 * D02.lam / 2 := by
    rw [Real.log_mul (by norm_num) st.ne', D02.log_s]; unfold D02.lam; ring
  have hX : Real.log (V + D02.a) + Real.log (V + D02.b) - Real.log V - 3 * D02.lam / 2
      = Real.log (1 + q) := by
    have := congrArg Real.log e1
    rw [Real.log_mul (by linarith [D02.a_pos]) (by linarith [D02.b_pos]),
      Real.log_mul (by positivity) (by linarith), Real.log_mul (by positivity) hV.ne', hl22] at this
    linarith
  have hY : Real.log (V ^ 2 + 1) - Real.log V - D02.lam / 2 = Real.log (2 * q) := by
    rw [e2, show 2 * Real.sqrt 2 * q * V = Real.sqrt 2 * ((2 * q) * V) by ring,
      Real.log_mul st.ne' (by positivity), Real.log_mul (by positivity) hV.ne', D02.log_s]
    ring
  have hVq : 2 + Real.sqrt 2 * ((2 * s - 1) / q) = Real.sqrt 2 * V / q := by
    rw [hVd]; field_simp; linear_combination (-q) * t2
  unfold gform S02.f1
  rw [hX, hY, e2, e1, e3, hVq]
  have q1 : q - 1 ≠ 0 := by linarith
  have q2 : 1 + q ≠ 0 := by linarith
  have q3 : 1 - q ^ 2 ≠ 0 := by nlinarith
  field_simp
  ring

lemma P1_eq_G : S02.P1 = G D02.a - G D02.b := by
  have hcont : ContinuousOn (fun s => G (Vs s)) (Icc 0 1) :=
    continuousOn_G.comp continuous_Vs.continuousOn (fun s hs => ⟨Vs_ge hs.1, Vs_le hs.2⟩)
  have hderiv : ∀ x ∈ Ioo (0 : ℝ) 1, HasDerivAt (fun s => G (Vs s)) (S02.f1 (S02.qs x)) x := by
    intro x hx
    have h := (hasDerivAt_G ⟨Vs_gt hx.1, Vs_lt hx.2⟩).comp x (hasDerivAt_Vs x)
    rw [gform_Vs hx] at h
    exact h
  unfold S02.P1
  rw [intervalIntegral.integral_eq_sub_of_hasDerivAt_of_le zero_le_one hcont hderiv S02.ii_f1,
    Vs_one, Vs_zero]

/-! ## Dilogarithm values -/

lemma s_ne : Real.sqrt 2 ≠ 0 := sqrt2_pos.ne'

lemma log_a_eq : Real.log D02.a = D02.L := rfl

lemma V1 : D00.Li2 (-(1 / Real.sqrt 2)) = -D00.Li2 D02.b - (1 / 2) * (D02.L - D02.lam / 2) ^ 2 := by
  have hs1 := D02.hs1
  have h := D00.landen' (x := -(1 / Real.sqrt 2)) (by
    have : 0 < 1 / Real.sqrt 2 := by positivity
    linarith)
  have e1 : -(1 / Real.sqrt 2) / (-(1 / Real.sqrt 2) - 1) = D02.b := by
    have hne : -(1 / Real.sqrt 2) - 1 ≠ 0 := by
      have : 0 < 1 / Real.sqrt 2 := by positivity
      linarith
    rw [div_eq_iff hne, b_def]; field_simp; linear_combination (1 : ℝ) * sqrt2_sq
  have e2 : Real.log (1 - -(1 / Real.sqrt 2)) = D02.L - D02.lam / 2 := by
    rw [show 1 - -(1 / Real.sqrt 2) = D02.a / Real.sqrt 2 by rw [a_def]; field_simp; ring,
      Real.log_div D02.a_pos.ne' s_ne, D02.log_s, log_a_eq]
  rw [e1, e2] at h
  linarith

lemma V3 : D00.Li2 (1 / Real.sqrt 2) = (1 / 2) * D00.Li2 (1 / 2) - D00.Li2 (-(1 / Real.sqrt 2)) := by
  have hs1 := D02.hs1
  have h := D00.duplication' (x := 1 / Real.sqrt 2) (by
    rw [abs_of_pos (by positivity), div_le_one sqrt2_pos]; linarith)
  rw [div_pow, one_pow, sqrt2_sq] at h
  linarith

lemma V2 : D00.Li2 (-D02.a) = -(1 / 2) * (D02.lam / 2 + D02.L) ^ 2 - D00.Li2 (1 / Real.sqrt 2) := by
  have hs1 := D02.hs1
  have h := D00.landen' (x := -D02.a) (by linarith [D02.a_pos])
  have e1 : -D02.a / (-D02.a - 1) = 1 / Real.sqrt 2 := by
    rw [a_def]
    have : -(1 + Real.sqrt 2) - 1 ≠ 0 := by linarith
    field_simp; linear_combination (-1 : ℝ) * sqrt2_sq
  have e2 : Real.log (1 - -D02.a) = D02.lam / 2 + D02.L := by
    rw [show 1 - -D02.a = Real.sqrt 2 * D02.a by rw [a_def]; linear_combination -sqrt2_sq,
      Real.log_mul s_ne D02.a_pos.ne', D02.log_s, log_a_eq]
  rw [e1, e2] at h
  linarith

lemma b_sq : D02.b ^ 2 = 3 - 2 * Real.sqrt 2 := by rw [b_def]; linear_combination sqrt2_sq

lemma V4a : D00.Li2 (2 * D02.b) = π ^ 2 / 6 + 2 * D02.L * (D02.lam - D02.L) - D00.Li2 (D02.b ^ 2) := by
  have hs1 := D02.hs1
  have hs2 := D02.hs2
  have hb2 := b_sq
  have h := D00.reflection' (x := D02.b ^ 2) (by positivity [D02.b_pos]) (by rw [hb2]; linarith)
  have e1 : 1 - D02.b ^ 2 = 2 * D02.b := by rw [hb2, b_def]; ring
  have e2 : Real.log (D02.b ^ 2) = -2 * D02.L := by rw [Real.log_pow, D02.log_b]; push_cast; ring
  have e3 : Real.log (2 * D02.b) = D02.lam - D02.L := by
    rw [Real.log_mul (by norm_num) D02.b_pos.ne', D02.log_b]; rfl
  rw [e1, e2, e3] at h
  linarith

lemma V4 : D00.Li2 (-(2 * D02.a)) = -2 * D02.L ^ 2 - D00.Li2 (2 * D02.b) := by
  have hs1 := D02.hs1
  have hs2 := D02.hs2
  have h := D00.landen' (x := 2 * D02.b) (by rw [b_def]; linarith)
  have e1 : 2 * D02.b / (2 * D02.b - 1) = -(2 * D02.a) := by
    rw [a_def, b_def]
    have : 2 * (Real.sqrt 2 - 1) - 1 ≠ 0 := by linarith
    field_simp; linear_combination (2 : ℝ) * sqrt2_sq
  have e2 : Real.log (1 - 2 * D02.b) = -2 * D02.L := by
    rw [show 1 - 2 * D02.b = D02.b ^ 2 by rw [b_sq, b_def]; ring, Real.log_pow, D02.log_b]
    push_cast; ring
  rw [e1, e2] at h
  linarith

lemma V5 : D00.Li2 (-D02.a ^ 2) = -(π ^ 2 / 6) - 2 * D02.L ^ 2 - D00.Li2 (-D02.b ^ 2) := by
  have h := D00.inversion_neg' (x := -D02.a ^ 2) (by have := D02.a_pos; nlinarith)
  have e1 : 1 / -D02.a ^ 2 = -D02.b ^ 2 := by
    have := D02.ab
    have ha := D02.a_pos.ne'
    field_simp
    linear_combination (D02.a * D02.b + 1) * this
  have e2 : Real.log (- -D02.a ^ 2) = 2 * D02.L := by
    rw [neg_neg, Real.log_pow, log_a_eq]; push_cast; ring
  rw [e1, e2] at h
  linarith

lemma V8 : D00.Li2 (D02.b / (2 * Real.sqrt 2)) = π ^ 2 / 6
    - (-D02.L - 3 * D02.lam / 2) * (D02.L - 3 * D02.lam / 2) - D00.Li2 ((2 + Real.sqrt 2) / 4) := by
  have hs1 := D02.hs1
  have hs2 := D02.hs2
  have h := D00.reflection' (x := D02.b / (2 * Real.sqrt 2)) (by positivity [D02.b_pos])
    (by rw [div_lt_one (by positivity), b_def]; linarith)
  have e1 : 1 - D02.b / (2 * Real.sqrt 2) = (2 + Real.sqrt 2) / 4 := by
    rw [b_def]; field_simp; linear_combination (-2 : ℝ) * sqrt2_sq
  have hl22 : Real.log (2 * Real.sqrt 2) = 3 * D02.lam / 2 := by
    rw [Real.log_mul (by norm_num) s_ne, D02.log_s]; unfold D02.lam; ring
  have e2 : Real.log (D02.b / (2 * Real.sqrt 2)) = -D02.L - 3 * D02.lam / 2 := by
    rw [Real.log_div D02.b_pos.ne' (by positivity), D02.log_b, hl22]
  have e3 : Real.log ((2 + Real.sqrt 2) / 4) = D02.L - 3 * D02.lam / 2 := by
    rw [show (2 + Real.sqrt 2) / 4 = Real.sqrt 2 * D02.a / 2 ^ 2 by
        rw [a_def]; linear_combination (-1 / 4 : ℝ) * sqrt2_sq,
      Real.log_div (by positivity [D02.a_pos]) (by norm_num), Real.log_mul s_ne D02.a_pos.ne',
      D02.log_s, log_a_eq, Real.log_pow]
    unfold D02.lam; push_cast; ring
  rw [e1, e2, e3] at h
  linarith

/-! ## `G(a)` and `G(b)` -/

lemma a_inv : D02.a⁻¹ = D02.b := inv_eq_of_mul_eq_one_right D02.ab
lemma b_inv : D02.b⁻¹ = D02.a := inv_eq_of_mul_eq_one_right (by rw [mul_comm]; exact D02.ab)

lemma arg1 : (D02.b - D02.a) / (D02.b + D02.a) = -(1 / Real.sqrt 2) := by
  have := s_ne
  rw [show D02.b - D02.a = -2 by rw [a_def, b_def]; ring,
    show D02.b + D02.a = 2 * Real.sqrt 2 by rw [a_def, b_def]; ring]
  field_simp
lemma arg2 : (D02.b - D02.a) / (D02.b + D02.b) = -D02.a := by
  have hb := D02.b_pos.ne'
  have := D02.ab
  rw [show D02.b - D02.a = -2 by rw [a_def, b_def]; ring]
  field_simp
  linear_combination (2 : ℝ) * this
lemma arg3 : (D02.b - D02.a) / (D02.b + 0) = -(2 * D02.a) := by
  have hb := D02.b_pos.ne'
  have := D02.ab
  rw [show D02.b - D02.a = -2 by rw [a_def, b_def]; ring, add_zero]
  field_simp
  linear_combination (1 : ℝ) * this
lemma arg4 : 2 / (D02.a + D02.a) = D02.b := by
  rw [show D02.a + D02.a = 2 * D02.a by ring, ← a_inv]; field_simp
lemma arg5 : D02.a / (D02.a + D02.a) = 1 / 2 := by
  have := D02.a_pos.ne'; field_simp; ring
lemma arg6 : (-D02.b - D02.a) / (-D02.b + D02.a) = -Real.sqrt 2 := by
  rw [a_def, b_def]; ring
lemma arg7 : D02.b / (D02.a + D02.b) = D02.b / (2 * Real.sqrt 2) := by
  rw [a_def, b_def]; ring_nf
lemma arg8 : (D02.a - D02.b) / (D02.a + D02.a) = D02.b := by
  rw [show D02.a - D02.b = 2 by rw [a_def, b_def]; ring]; exact arg4
lemma arg9 : (D02.a - D02.b) / (D02.a + D02.b) = 1 / Real.sqrt 2 := by
  have := s_ne
  rw [show D02.a - D02.b = 2 by rw [a_def, b_def]; ring,
    show D02.a + D02.b = 2 * Real.sqrt 2 by rw [a_def, b_def]; ring]
  field_simp
lemma arg10 : 2 / (D02.b + D02.a) = 1 / Real.sqrt 2 := by
  have := s_ne
  rw [show D02.b + D02.a = 2 * Real.sqrt 2 by rw [a_def, b_def]; ring]
  field_simp
lemma arg11 : D02.a / (D02.b + D02.a) = (2 + Real.sqrt 2) / 4 := by
  have := s_ne
  rw [show D02.b + D02.a = 2 * Real.sqrt 2 by rw [a_def, b_def]; ring, a_def]
  field_simp
  linear_combination (-2 : ℝ) * sqrt2_sq
lemma arg12 : (-D02.b - D02.b) / (-D02.b + D02.a) = -D02.b := by
  rw [a_def, b_def]; ring
lemma arg13 : D02.b / (D02.b + D02.b) = 1 / 2 := by
  have := D02.b_pos.ne'; field_simp; ring

lemma G_a : G D02.a = -D00.Li2 (-(1 / Real.sqrt 2)) - D00.Li2 (-D02.a) + 2 * T D02.b
    + 2 * D00.Li2 (-(2 * D02.a)) + (1 / 2 - Real.sqrt 2 / 4) * D00.Li2 (-D02.a ^ 2)
    + (1 - Real.sqrt 2 / 4) * D02.L ^ 2 + (3 * D02.lam / 2 + D02.L - Real.sqrt 2 * D02.lam / 4) * D02.L
    + (D02.lam + D02.L) ^ 2 / 2 + D00.Li2 D02.b - D00.Li2 (1 / 2) - 3 / 2 * D02.lam * (D02.lam + D02.L)
    + D00.Li2 (-Real.sqrt 2) + D00.Li2 (D02.b / (2 * Real.sqrt 2)) + D02.lam / 2 * (3 * D02.lam / 2) := by
  unfold G
  rw [T_a, arg1, arg2, a_inv, arg3, log_a_eq, log_a_a, arg4, arg5, arg6, arg7,
    add_comm D02.a D02.b, log_b_a]
  simp only [sub_self, zero_div, D00.Li2_zero']
  ring

lemma G_b : G D02.b = D00.Li2 D02.b + D00.Li2 (1 / Real.sqrt 2) - 2 * T D02.b
    + (1 / 2 - Real.sqrt 2 / 4) * D00.Li2 (-D02.b ^ 2)
    + (1 - Real.sqrt 2 / 4) * D02.L ^ 2 + (3 * D02.lam / 2 + D02.L - Real.sqrt 2 * D02.lam / 4) * (-D02.L)
    + (3 * D02.lam / 2) ^ 2 / 2 + D00.Li2 (1 / Real.sqrt 2) - D00.Li2 ((2 + Real.sqrt 2) / 4)
    - 3 / 2 * D02.lam * (3 * D02.lam / 2)
    + D00.Li2 (-D02.b) + D00.Li2 (1 / 2) + D02.lam / 2 * (D02.lam - D02.L) := by
  unfold G
  rw [b_inv, T_a, arg8, arg9, D02.log_b, log_b_a, arg10, arg11, arg12, arg13, log_b_b]
  simp only [sub_self, zero_div, D00.Li2_zero']
  ring

/-! ## The closed form of `P₁` -/

lemma P1_closed : S02.P1 = 3 * Real.sqrt 2 * D00.Li2 D02.b - D02.lam ^ 2 / 2
      - Real.sqrt 2 / 2 * D02.lam * D02.L + (1 + 3 * Real.sqrt 2 / 2) * D02.L ^ 2
      - (1 / 48 + 3 * Real.sqrt 2 / 16) * π ^ 2 := by
  rw [P1_eq_G, G_a, G_b, T_b]
  have v1 := V1
  have v2 := V2
  have v3 := V3
  have v4 := V4
  have v4a := V4a
  have v5 := V5
  have v8 := V8
  have d1 := D02.proof_1
  have d2 := D02.proof_2
  have d4 := D02.proof_4
  have d5 := D02.proof_5
  have d6 := D02.proof_6
  have h12 : D00.Li2 (1 / 2) = π ^ 2 / 12 - (1 / 2) * D02.lam ^ 2 := D00.proof_14
  unfold D02.stmt_1 at d1
  unfold D02.stmt_2 at d2
  unfold D02.stmt_4 at d4
  unfold D02.stmt_5 at d5
  unfold D02.stmt_6 at d6
  linear_combination (1 / 2 - Real.sqrt 2 / 4) * v5 + v8 + 0 * d6 + 2 * v4 - 2 * v4a + 2 * d2
    - 2 * (1 / 2 - Real.sqrt 2 / 4) * d5 - v2 - v3 + 0 * v1 - 5 / 2 * h12 + d4 + d1

/-- **Part C.** `S06.stmt_6`: the integrand of `P₁` is integrable and `P₁` has the closed form
of Lemma `g3:lem:P`. -/
theorem main : S06.stmt_6 := ⟨S02.ii_f1, P1_closed⟩

end S06.PC

/-! ### Merged from work/S06/PartD.lean -/
/-!
# S06 Part D: the closed form of `P_{2a}` (Lemma `g3:lem:P`, second line)

Route (all real). Substitute `s = S(z) = (1 + T(z))/2`, `T(z) = -(z²-4z+1)/(2√C(z))`,
`C(z) = (z-2)(z²-1)`, `z ∈ [a, 5+3√2]` (`S(a) = 1`, `S(5+3√2) = 0`). Then
`q(S(z)) = (z²-2z+3)/(2√2√C)`, `√(1+q²) = (z²+2z-5)/(2√2√C)`, `arsinh q = ½ log((z²-1)/(2(z-2)))`,
and the pulled-back integrand is
`R(z) = -(√2/8)(log(z-1)+log(z+1)-log2-log(z-2))(1/(z+1)-1/(z-1)+1/(z-2))`,
with an explicit primitive `Φ` (logs, `Li₂(-(z-1)/2)`, `Li₂((z-2)/(z+1))`). The end values give
the paper's bracket `Li₂(-1/√2)-Li₂(1/√2)+Li₂(b²/√2)-Li₂(-a²/√2)`, reduced by Landen,
duplication, the five-term relation, inversion, reflection and D02.
-/

namespace S06.PD

open Real MeasureTheory Set

noncomputable def Cf (z : ℝ) : ℝ := (z - 2) * (z ^ 2 - 1)

noncomputable def Tf (z : ℝ) : ℝ := -(z ^ 2 - 4 * z + 1) / (2 * Real.sqrt (Cf z))

noncomputable def Sf (z : ℝ) : ℝ := (1 + Tf z) / 2

noncomputable def Sd (z : ℝ) : ℝ :=
  -((z ^ 2 - 2 * z + 3) * (z ^ 2 + 2 * z - 5)) / (8 * Real.sqrt (Cf z) * Cf z)

noncomputable def gf (s : ℝ) : ℝ := Real.arsinh (S02.qs s) / Real.sqrt (1 + S02.qs s ^ 2)

noncomputable def ell (z : ℝ) : ℝ :=
  Real.log (z - 1) + Real.log (z + 1) - Real.log 2 - Real.log (z - 2)

noncomputable def Rf (z : ℝ) : ℝ :=
  -(Real.sqrt 2 / 8) * ell z * (1 / (z + 1) - 1 / (z - 1) + 1 / (z - 2))

lemma Cf_pos {z : ℝ} (hz : 2 < z) : 0 < Cf z := by
  unfold Cf
  have : 0 < z ^ 2 - 1 := by nlinarith
  exact mul_pos (by linarith) this

lemma hasDerivAt_Cf (z : ℝ) : HasDerivAt Cf (3 * z ^ 2 - 4 * z - 1) z := by
  have h := ((hasDerivAt_id z).sub_const 2).mul ((hasDerivAt_pow 2 z).sub_const 1)
  exact h.congr_deriv (by simp; ring)

lemma hasDerivAt_Sf {z : ℝ} (hz : 2 < z) : HasDerivAt Sf (Sd z) z := by
  have hC := Cf_pos hz
  have hr := (hasDerivAt_Cf z).sqrt hC.ne'
  have hr0 : 0 < Real.sqrt (Cf z) := Real.sqrt_pos.mpr hC
  have hN : HasDerivAt (fun z : ℝ => -(z ^ 2 - 4 * z + 1)) (-(2 * z - 4)) z := by
    have := (((hasDerivAt_pow 2 z).sub ((hasDerivAt_id z).const_mul 4)).add_const 1).neg
    exact this.congr_deriv (by simp)
  have hT := hN.div (hr.const_mul 2) (by positivity)
  have hS := (hT.const_add 1).div_const 2
  refine hS.congr_deriv ?_
  have hr2 : Real.sqrt (Cf z) ^ 2 = Cf z := Real.sq_sqrt hC.le
  set r := Real.sqrt (Cf z) with hrdef
  have e1 : ((-(2 * z - 4) * (2 * r) - -(z ^ 2 - 4 * z + 1) * (2 * ((3 * z ^ 2 - 4 * z - 1) / (2 * r))))
      / (2 * r) ^ 2) / 2
      = (2 * (-(2 * z - 4)) * r ^ 2 - (-(z ^ 2 - 4 * z + 1)) * (3 * z ^ 2 - 4 * z - 1)) / (8 * r ^ 3) := by
    field_simp; ring
  rw [e1, show r ^ 3 = r ^ 2 * r by ring, hr2]
  unfold Sd
  rw [← hrdef]
  unfold Cf at hC ⊢
  have : z - 2 ≠ 0 := by linarith
  have : z ^ 2 - 1 ≠ 0 := by nlinarith
  field_simp
  ring

lemma s2_facts : Real.sqrt 2 ^ 2 = 2 ∧ 0 < Real.sqrt 2 :=
  ⟨Real.sq_sqrt (by norm_num), Real.sqrt_pos.mpr (by norm_num)⟩

lemma qs_Sf {z : ℝ} (hz : 2 < z) :
    S02.qs (Sf z) = (z ^ 2 - 2 * z + 3) / (2 * Real.sqrt 2 * Real.sqrt (Cf z)) := by
  have hC := Cf_pos hz
  have hr2 : Real.sqrt (Cf z) ^ 2 = Cf z := Real.sq_sqrt hC.le
  have hr0 : 0 < Real.sqrt (Cf z) := Real.sqrt_pos.mpr hC
  obtain ⟨hs2, hs0⟩ := s2_facts
  have hP : 0 < z ^ 2 - 2 * z + 3 := by nlinarith [sq_nonneg (z - 1)]
  unfold S02.qs
  rw [show Sf z ^ 2 + (1 - Sf z) ^ 2 = ((z ^ 2 - 2 * z + 3) / (2 * Real.sqrt 2 * Real.sqrt (Cf z))) ^ 2 from ?_]
  · exact Real.sqrt_sq (by positivity)
  have e1 : Sf z ^ 2 + (1 - Sf z) ^ 2 = (1 + Tf z ^ 2) / 2 := by unfold Sf; ring
  have e2 : Tf z ^ 2 = (z ^ 2 - 4 * z + 1) ^ 2 / (4 * Cf z) := by
    unfold Tf; rw [div_pow, mul_pow, hr2, neg_sq]; norm_num
  rw [e1, e2, div_pow, mul_pow, mul_pow, hr2, hs2]
  unfold Cf at hC ⊢
  have : z - 2 ≠ 0 := by linarith
  have : z ^ 2 - 1 ≠ 0 := by nlinarith
  field_simp
  ring

lemma sq1_Sf {z : ℝ} (hz : 2 < z) :
    Real.sqrt (1 + S02.qs (Sf z) ^ 2) = (z ^ 2 + 2 * z - 5) / (2 * Real.sqrt 2 * Real.sqrt (Cf z)) := by
  have hC := Cf_pos hz
  have hr2 : Real.sqrt (Cf z) ^ 2 = Cf z := Real.sq_sqrt hC.le
  have hr0 : 0 < Real.sqrt (Cf z) := Real.sqrt_pos.mpr hC
  obtain ⟨hs2, hs0⟩ := s2_facts
  have hB : 0 < z ^ 2 + 2 * z - 5 := by nlinarith
  rw [qs_Sf hz]
  rw [show 1 + ((z ^ 2 - 2 * z + 3) / (2 * Real.sqrt 2 * Real.sqrt (Cf z))) ^ 2
      = ((z ^ 2 + 2 * z - 5) / (2 * Real.sqrt 2 * Real.sqrt (Cf z))) ^ 2 from ?_]
  · exact Real.sqrt_sq (by positivity)
  rw [div_pow, div_pow, mul_pow, mul_pow, hr2, hs2]
  unfold Cf at hC ⊢
  have : z - 2 ≠ 0 := by linarith
  have : z ^ 2 - 1 ≠ 0 := by nlinarith
  field_simp
  ring

lemma arsinh_Sf {z : ℝ} (hz : 2 < z) : Real.arsinh (S02.qs (Sf z)) = ell z / 2 := by
  have hC := Cf_pos hz
  have hr2 : Real.sqrt (Cf z) ^ 2 = Cf z := Real.sq_sqrt hC.le
  have hr0 : 0 < Real.sqrt (Cf z) := Real.sqrt_pos.mpr hC
  obtain ⟨hs2, hs0⟩ := s2_facts
  unfold Real.arsinh
  rw [sq1_Sf hz, qs_Sf hz]
  set y := (z ^ 2 - 2 * z + 3) / (2 * Real.sqrt 2 * Real.sqrt (Cf z))
    + (z ^ 2 + 2 * z - 5) / (2 * Real.sqrt 2 * Real.sqrt (Cf z)) with hy
  have hy0 : 0 < y := by
    rw [hy, ← add_div]; apply div_pos _ (by positivity); nlinarith
  have hy2 : y ^ 2 = (z - 1) * (z + 1) / (2 * (z - 2)) := by
    rw [hy, ← add_div, div_pow, mul_pow, mul_pow, hr2, hs2]
    unfold Cf
    have : z - 2 ≠ 0 := by linarith
    have : z ^ 2 - 1 ≠ 0 := by nlinarith
    field_simp
    ring
  have hl : Real.log (y ^ 2) = 2 * Real.log y := by
    rw [Real.log_pow]; norm_num
  have h1 : (0 : ℝ) < z - 1 := by linarith
  have h2 : (0 : ℝ) < z + 1 := by linarith
  have h3 : (0 : ℝ) < z - 2 := by linarith
  rw [hy2, Real.log_div (by positivity) (by positivity), Real.log_mul h1.ne' h2.ne',
    Real.log_mul (by norm_num) h3.ne'] at hl
  unfold ell
  linarith

lemma gf_Sf_mul {z : ℝ} (hz : 2 < z) : gf (Sf z) * Sd z = Rf z := by
  have hC := Cf_pos hz
  have hr0 : 0 < Real.sqrt (Cf z) := Real.sqrt_pos.mpr hC
  obtain ⟨hs2, hs0⟩ := s2_facts
  have hB : 0 < z ^ 2 + 2 * z - 5 := by nlinarith
  unfold gf
  rw [arsinh_Sf hz, sq1_Sf hz]
  have cancel : ∀ X B P k m : ℝ, B ≠ 0 → X / (B / k) * (-(P * B) / m) = X * k * (-P) / m := by
    intro X B P k m hB0
    rw [div_div_eq_mul_div, mul_div_assoc]
    rw [show X * (k / B) * (-(P * B) / m) = X * k * -P / m * (B / B) by ring, div_self hB0, mul_one]
  unfold Sd
  rw [cancel _ _ _ _ _ hB.ne']
  unfold Rf
  unfold Cf at hC ⊢
  have : z - 2 ≠ 0 := by linarith
  have : z - 1 ≠ 0 := by linarith
  have : z + 1 ≠ 0 := by linarith
  have : z ^ 2 - 1 ≠ 0 := by nlinarith
  have : z ^ 2 + 2 * z - 5 ≠ 0 := hB.ne'
  field_simp
  ring

/-- The primitive of `Rf`. -/
noncomputable def Phi (z : ℝ) : ℝ :=
  Real.sqrt 2 * (Real.log (z - 1) ^ 2 / 16 - Real.log (z - 1) * Real.log (z + 1) / 8
    - Real.log (z - 1) * Real.log (z - 2) / 8 - 3 * Real.log (z + 1) ^ 2 / 16
    + Real.log (z + 1) * Real.log (z - 2) / 8 + Real.log (z - 2) ^ 2 / 16
    + Real.log 2 * Real.log (z - 1) / 8 + Real.log 2 * Real.log (z + 1) / 8
    + Real.log 2 * Real.log (z - 2) / 8 + Real.log 3 * Real.log (z + 1) / 4
    - Real.log 3 * Real.log (z - 2) / 4 - D00.Li2 (-(z - 1) / 2) / 4
    - D00.Li2 ((z - 2) / (z + 1)) / 4)

lemma hasDerivAt_Phi {z : ℝ} (hz : 2 < z) : HasDerivAt Phi (Rf z) z := by
  have p1 : (0 : ℝ) < z - 1 := by linarith
  have p2 : (0 : ℝ) < z + 1 := by linarith
  have p3 : (0 : ℝ) < z - 2 := by linarith
  have h1 : HasDerivAt (fun y : ℝ => Real.log (y - 1)) (1 / (z - 1)) z :=
    (((hasDerivAt_id z).sub_const 1).log (by simp; linarith)).congr_deriv (by simp)
  have hp : HasDerivAt (fun y : ℝ => Real.log (y + 1)) (1 / (z + 1)) z :=
    (((hasDerivAt_id z).add_const 1).log (by simp; linarith)).congr_deriv (by simp)
  have h2 : HasDerivAt (fun y : ℝ => Real.log (y - 2)) (1 / (z - 2)) z :=
    (((hasDerivAt_id z).sub_const 2).log (by simp; linarith)).congr_deriv (by simp)
  have hfA : HasDerivAt (fun y : ℝ => -(y - 1) / 2) (-1 / 2) z :=
    ((((hasDerivAt_id z).sub_const 1).neg).div_const 2).congr_deriv (by simp)
  have hA : HasDerivAt (fun y : ℝ => D00.Li2 (-(y - 1) / 2))
      (-Real.log (1 - -(z - 1) / 2) / (-(z - 1) / 2) * (-1 / 2)) z :=
    D00.Li2_comp hfA (by intro h; linarith) (by intro h; linarith)
  have hfB : HasDerivAt (fun y : ℝ => (y - 2) / (y + 1)) (3 / (z + 1) ^ 2) z :=
    (((hasDerivAt_id z).sub_const 2).div ((hasDerivAt_id z).add_const 1) (by simp; linarith)).congr_deriv
      (by simp; field_simp; ring)
  have hB1 : (z - 2) / (z + 1) ≠ 1 := by
    rw [Ne, div_eq_one_iff_eq p2.ne']; intro h; linarith
  have hB : HasDerivAt (fun y : ℝ => D00.Li2 ((y - 2) / (y + 1)))
      (-Real.log (1 - (z - 2) / (z + 1)) / ((z - 2) / (z + 1)) * (3 / (z + 1) ^ 2)) z :=
    D00.Li2_comp hfB hB1 (div_pos p3 p2).ne'
  have H := ((((((((((((((h1.pow 2).div_const 16).sub ((h1.mul hp).div_const 8)).sub
    ((h1.mul h2).div_const 8)).sub (((hp.pow 2).const_mul 3).div_const 16)).add
    ((hp.mul h2).div_const 8)).add ((h2.pow 2).div_const 16)).add
    ((h1.const_mul (Real.log 2)).div_const 8)).add ((hp.const_mul (Real.log 2)).div_const 8)).add
    ((h2.const_mul (Real.log 2)).div_const 8)).add ((hp.const_mul (Real.log 3)).div_const 4)).sub
    ((h2.const_mul (Real.log 3)).div_const 4)).sub (hA.div_const 4)).sub (hB.div_const 4)).const_mul
    (Real.sqrt 2)
  refine H.congr_deriv ?_
  have e1 : Real.log (1 - -(z - 1) / 2) = Real.log (z + 1) - Real.log 2 := by
    rw [show 1 - -(z - 1) / 2 = (z + 1) / 2 by ring, Real.log_div p2.ne' (by norm_num)]
  have e2 : Real.log (1 - (z - 2) / (z + 1)) = Real.log 3 - Real.log (z + 1) := by
    rw [show 1 - (z - 2) / (z + 1) = 3 / (z + 1) by field_simp; ring, Real.log_div (by norm_num) p2.ne']
  rw [e1, e2]
  unfold Rf ell
  field_simp
  ring

/-- The upper end point `z₁ = 5 + 3√2` (where `S = 0`). -/
noncomputable def z1 : ℝ := 5 + 3 * Real.sqrt 2

lemma a_gt_two : 2 < D02.a := by unfold D02.a; linarith [D02.hs1]

lemma a_le_z1 : D02.a ≤ z1 := by unfold D02.a z1; linarith [D02.hs1]

lemma Sf_a : Sf D02.a = 1 := by
  obtain ⟨hs2, hs0⟩ := s2_facts
  have hC : Cf D02.a = Real.sqrt 2 ^ 2 := by unfold Cf D02.a; ring_nf; rw [hs2]; nlinarith [hs2]
  unfold Sf Tf
  rw [hC, Real.sqrt_sq hs0.le]
  unfold D02.a
  field_simp
  ring_nf
  rw [hs2]
  ring

lemma Sf_z1 : Sf z1 = 0 := by
  obtain ⟨hs2, hs0⟩ := s2_facts
  have h3 : Real.sqrt 2 ^ 3 = 2 * Real.sqrt 2 := by rw [pow_succ, hs2]
  have h4 : Real.sqrt 2 ^ 4 = 4 := by rw [show (4 : ℕ) = 2 * 2 by rfl, pow_mul, hs2]; norm_num
  have hC : Cf z1 = (12 + 9 * Real.sqrt 2) ^ 2 := by
    unfold Cf z1; ring_nf; rw [hs2, h3]; ring
  unfold Sf Tf
  rw [hC, Real.sqrt_sq (by positivity)]
  unfold z1
  field_simp
  ring_nf
  rw [hs2]; ring

lemma continuous_gf : Continuous gf := by
  unfold gf
  refine (Real.continuous_arsinh.comp S02.continuous_qs).div (by have := S02.continuous_qs; fun_prop)
    (fun s => ?_)
  exact (Real.sqrt_pos.mpr (by positivity)).ne'

lemma mem_I {x : ℝ} (hx : x ∈ uIcc z1 D02.a) : 2 < x := by
  rw [uIcc_of_ge a_le_z1] at hx
  linarith [hx.1, a_gt_two]

lemma continuousOn_Sd : ContinuousOn Sd (uIcc z1 D02.a) := by
  intro x hx
  have h2 := mem_I hx
  have hC := Cf_pos h2
  apply ContinuousAt.continuousWithinAt
  show ContinuousAt (fun x => -((x ^ 2 - 2 * x + 3) * (x ^ 2 + 2 * x - 5))
    / (8 * Real.sqrt (Cf x) * Cf x)) x
  have hd : ContinuousAt (fun x => 8 * Real.sqrt (Cf x) * Cf x) x := by unfold Cf; fun_prop
  have hn : Continuous (fun x : ℝ => -((x ^ 2 - 2 * x + 3) * (x ^ 2 + 2 * x - 5))) := by fun_prop
  have := Real.sqrt_pos.mpr hC
  exact hn.continuousAt.div hd (mul_pos (mul_pos (by norm_num) this) hC).ne'

lemma P2a_eq_Phi : S02.P2a = Phi D02.a - Phi z1 := by
  have hderS : ∀ x ∈ uIcc z1 D02.a, HasDerivAt Sf (Sd x) x := fun x hx => hasDerivAt_Sf (mem_I hx)
  have hcv := intervalIntegral.integral_comp_mul_deriv hderS continuousOn_Sd continuous_gf
  rw [Sf_z1, Sf_a] at hcv
  have hderP : ∀ x ∈ uIcc z1 D02.a, HasDerivAt Phi ((gf ∘ Sf) x * Sd x) x := fun x hx =>
    (hasDerivAt_Phi (mem_I hx)).congr_deriv (gf_Sf_mul (mem_I hx)).symm
  have hSc : ContinuousOn Sf (uIcc z1 D02.a) := fun x hx =>
    (hderS x hx).continuousAt.continuousWithinAt
  have hint : IntervalIntegrable (fun x => (gf ∘ Sf) x * Sd x) volume z1 D02.a :=
    ((continuous_gf.comp_continuousOn hSc).mul continuousOn_Sd).intervalIntegrable
  have hFTC := intervalIntegral.integral_eq_sub_of_hasDerivAt hderP hint
  rw [hcv] at hFTC
  exact hFTC

/-! ## End values -/

lemma s2_pows : Real.sqrt 2 ^ 2 = 2 ∧ Real.sqrt 2 ^ 3 = 2 * Real.sqrt 2 ∧ Real.sqrt 2 ^ 4 = 4 := by
  have hs2 : Real.sqrt 2 ^ 2 = 2 := Real.sq_sqrt (by norm_num)
  refine ⟨hs2, by rw [pow_succ, hs2], by rw [show (4 : ℕ) = 2 * 2 by rfl, pow_mul, hs2]; norm_num⟩

lemma e_am1 : D02.a - 1 = Real.sqrt 2 := by unfold D02.a; ring
lemma e_ap1 : D02.a + 1 = Real.sqrt 2 * D02.a := by
  obtain ⟨h2, -, -⟩ := s2_pows; unfold D02.a; ring_nf; rw [h2]; ring
lemma e_am2 : D02.a - 2 = D02.b := by unfold D02.a D02.b; ring
lemma e_zm1 : z1 - 1 = Real.sqrt 2 * D02.a ^ 2 := by
  obtain ⟨h2, h3, -⟩ := s2_pows; unfold D02.a z1; ring_nf; rw [h2, h3]; ring
lemma e_zp1 : z1 + 1 = 3 * (Real.sqrt 2 * D02.a) := by
  obtain ⟨h2, -, -⟩ := s2_pows; unfold D02.a z1; ring_nf; rw [h2]; ring
lemma e_zm2 : z1 - 2 = 3 * D02.a := by unfold D02.a z1; ring

lemma log_am1 : Real.log (D02.a - 1) = D02.lam / 2 := by rw [e_am1, D02.log_s]
lemma log_ap1 : Real.log (D02.a + 1) = D02.lam / 2 + D02.L := by
  rw [e_ap1, Real.log_mul (by positivity) D02.a_pos.ne', D02.log_s]; rfl
lemma log_am2 : Real.log (D02.a - 2) = -D02.L := by rw [e_am2, D02.log_b]
lemma log_zm1 : Real.log (z1 - 1) = D02.lam / 2 + 2 * D02.L := by
  rw [e_zm1, Real.log_mul (by positivity) (by have := D02.a_pos; positivity), D02.log_s,
    Real.log_pow]; unfold D02.L; push_cast; ring
lemma log_zp1 : Real.log (z1 + 1) = Real.log 3 + (D02.lam / 2 + D02.L) := by
  rw [e_zp1, Real.log_mul (by norm_num) (by have := D02.a_pos; positivity),
    Real.log_mul (by positivity) D02.a_pos.ne', D02.log_s]; rfl
lemma log_zm2 : Real.log (z1 - 2) = Real.log 3 + D02.L := by
  rw [e_zm2, Real.log_mul (by norm_num) D02.a_pos.ne']; rfl

lemma ab' : D02.a * D02.b = 1 := D02.ab

lemma b_sq : D02.b ^ 2 = 1 - 2 * D02.b := by
  obtain ⟨h2, -, -⟩ := s2_pows; unfold D02.b; ring_nf; rw [h2]; ring

lemma a_sq : 1 + D02.a ^ 2 = 2 * Real.sqrt 2 * D02.a := by
  obtain ⟨h2, -, -⟩ := s2_pows; unfold D02.a; ring_nf; rw [h2]; ring

lemma arg1 : (-(D02.a - 1) / 2) / (-(D02.a - 1) / 2 - 1) = D02.b := by
  obtain ⟨h2, h3, -⟩ := s2_pows
  rw [e_am1, div_eq_iff (by nlinarith [D02.hs1])]
  unfold D02.b; ring_nf; rw [h2]; ring

lemma log1_arg1 : Real.log (1 - -(D02.a - 1) / 2) = D02.L - D02.lam / 2 := by
  obtain ⟨h2, -, -⟩ := s2_pows
  have : 1 - -(D02.a - 1) / 2 = D02.a / Real.sqrt 2 := by
    rw [e_am1, eq_div_iff (by positivity)]; unfold D02.a; ring_nf; rw [h2]; ring
  rw [this, Real.log_div D02.a_pos.ne' (by positivity), D02.log_s]; rfl

lemma yZ_eq : (z1 - 2) / (z1 + 1) = Real.sqrt 2 / 2 := by
  obtain ⟨h2, -, -⟩ := s2_pows
  rw [e_zm2, e_zp1, div_eq_iff (by have := D02.a_pos; positivity)]
  ring_nf; rw [h2]; ring

lemma arg_five1 : -D02.a ^ 2 / (1 - 2 * D02.b) = -D02.a ^ 4 := by
  rw [← b_sq, div_eq_iff (by have := D02.b_pos; positivity)]
  have := ab'
  linear_combination (D02.a ^ 2 * (D02.a * D02.b + 1)) * this

lemma arg_five2 : 2 * D02.b / (1 - -D02.a ^ 2) = (D02.a - 2) / (D02.a + 1) := by
  obtain ⟨h2, h3, -⟩ := s2_pows
  rw [e_am2, e_ap1, sub_neg_eq_add, a_sq]
  have := D02.a_pos
  field_simp

lemma arg_five3 : -D02.a ^ 2 * (2 * D02.b) / ((1 - -D02.a ^ 2) * (1 - 2 * D02.b)) = -(z1 - 1) / 2 := by
  obtain ⟨h2, h3, -⟩ := s2_pows
  rw [e_zm1, sub_neg_eq_add, a_sq, ← b_sq]
  have := D02.a_pos
  have hb := D02.b_pos
  have hab := ab'
  rw [div_eq_div_iff (by positivity) (by norm_num)]
  linear_combination (4 * D02.a ^ 2 * D02.b) * hab + (2 * D02.a ^ 3 * D02.b ^ 2) * h2

lemma log_five1 : Real.log (1 - -D02.a ^ 2) = 3 * D02.lam / 2 + D02.L := by
  rw [sub_neg_eq_add, a_sq, Real.log_mul (by positivity) D02.a_pos.ne',
    Real.log_mul (by norm_num) (by positivity), D02.log_s]
  unfold D02.lam D02.L; ring

lemma log_five2 : Real.log (1 - 2 * D02.b) = -2 * D02.L := by
  rw [← b_sq, Real.log_pow, D02.log_b]; push_cast; ring

lemma inv_a4 : 1 / -D02.a ^ 4 = -D02.b ^ 4 := by
  have h := ab'
  have : D02.a ^ 4 * D02.b ^ 4 = 1 := by rw [← mul_pow, h, one_pow]
  have ha := D02.a_pos
  field_simp
  linear_combination this

lemma inv_a2 : 1 / -D02.a ^ 2 = -D02.b ^ 2 := by
  have h := ab'
  have : D02.a ^ 2 * D02.b ^ 2 = 1 := by rw [← mul_pow, h, one_pow]
  have ha := D02.a_pos
  field_simp
  linear_combination this

lemma bracket :
    D00.Li2 (-(D02.a - 1) / 2) - D00.Li2 ((z1 - 2) / (z1 + 1)) + D00.Li2 ((D02.a - 2) / (D02.a + 1))
      - D00.Li2 (-(z1 - 1) / 2)
      = -12 * D00.Li2 D02.b + D00.Li2 (-D02.b ^ 4) - 2 * D02.L ^ 2 + 5 * π ^ 2 / 6 := by
  have hs1 := D02.hs1
  have hs2' := D02.hs2
  have ha := D02.a_pos
  have hb := D02.b_pos
  have hb1 := D02.b_lt
  -- (2) Landen at x = -1/√2
  have I1 := D00.landen (x := -(D02.a - 1) / 2) (by rw [e_am1]; linarith)
  rw [arg1, log1_arg1] at I1
  -- (3) duplication at x = 1/√2
  have I2 := D00.duplication (x := (z1 - 2) / (z1 + 1))
    (by rw [yZ_eq, abs_le]; constructor <;> linarith)
  rw [show -((z1 - 2) / (z1 + 1)) = -(D02.a - 1) / 2 by rw [yZ_eq, e_am1]; ring,
    show ((z1 - 2) / (z1 + 1)) ^ 2 = 1 / 2 by rw [yZ_eq, div_pow, s2_pows.1]; norm_num] at I2
  -- (15) five-term at x = -a², y = 2b
  have hbh : D02.b < 1 / 2 := by unfold D02.b; linarith
  have I3 := D00.five_term (x := -D02.a ^ 2) (y := 2 * D02.b) (by nlinarith) (by linarith)
    (by nlinarith)
  rw [arg_five1, arg_five2, arg_five3, log_five1, log_five2] at I3
  -- (1) inversion at -a⁴
  have I4 := D00.inversion_neg (x := -D02.a ^ 4) (by have := pow_pos ha 4; linarith)
  rw [inv_a4, neg_neg, Real.log_pow] at I4
  -- inversion at -a²
  have I5 := D00.inversion_neg (x := -D02.a ^ 2) (by have := pow_pos ha 2; linarith)
  rw [inv_a2, neg_neg, Real.log_pow] at I5
  have I6 := D02.proof_5
  -- reflection at b²
  have I7 := D00.reflection (x := D02.b ^ 2) (by positivity) (by nlinarith)
  rw [show 1 - D02.b ^ 2 = 2 * D02.b by rw [b_sq]; ring, Real.log_pow, D02.log_b,
    Real.log_mul (by norm_num) hb.ne', D02.log_b] at I7
  have I8 := D02.proof_2
  have I9 := D00.Li2_half
  unfold D02.stmt_5 at I6
  unfold D02.stmt_2 at I8
  unfold D02.L D02.lam at *
  push_cast at I4 I5 I7
  linear_combination 2 * I1 - I2 + I3 - I4 + I5 - I6 + I7 - I8 - (1 / 2) * I9

lemma Phi_diff : Phi D02.a - Phi z1 = Real.sqrt 2 * (3 * D00.Li2 D02.b - (1 / 4) * D00.Li2 (-D02.b ^ 4)
      - D02.lam * D02.L / 2 + D02.L * Real.log 3 + D02.L ^ 2 / 2 - 5 * π ^ 2 / 24) := by
  have B := bracket
  unfold Phi
  rw [log_am1, log_ap1, log_am2, log_zm1, log_zp1, log_zm2]
  rw [show Real.log 2 = D02.lam from rfl]
  linear_combination (-(Real.sqrt 2 / 4)) * B

/-- Part D: Lemma `g3:lem:P`, second line (`S06.stmt_7`). -/
theorem main : S06.stmt_7 := ⟨S02.ii_2a, by rw [P2a_eq_Phi, Phi_diff]⟩

end S06.PD

/-! ### Merged from work/S06/PartE.lean -/
/-!
# S06, part E: the closed form of `P_{2c}` (`g3:lem:P`, third line)

Route (real-variable reduction plus a complex Landen identity proved by a derivative argument).
* `χ(z) = ½ Σ (z^m - (-z)^m)/m²` (complex power series, closed unit disc); `χ'(z) = ℓ(z)/(2z)`,
  `ℓ(z) = log(1+z) - log(1-z)`.
* `x(s) = (2s - 1 + 2√(s²-s+1))/√3` maps `[0,½]` onto `[1/√3, 1]`, and
  `d/ds Re χ(x(s) ω) = φ(s)/√2` (`ω = e^{iα_c}`), so `P_{2c} = 2√2 (Re χ(ω) - Re χ(ω/√3))`.
* `Re χ(ω) = π²/8 - πα_c/4` (Fourier series of `B₂`).
* Landen for `χ` on `{Re z > 0, |z| < 1}` (derivative zero, value at `z = b`), at `t = ω/√3`,
  with `(1-t)/(1+t) = b e^{-iπ/4}`, and `Re χ(b e^{-iπ/4}) = ½Li₂(b) + π²/192 - Li₂(-b⁴)/16`
  (Lemma A plus the even part of the series).
-/

namespace S06.PE

open Real MeasureTheory Set Filter Topology

/-! ## 1. The complex `χ` series -/

noncomputable def chiT (z : ℂ) (n : ℕ) : ℂ :=
  (z ^ (n + 1) - (-z) ^ (n + 1)) / (2 * ((n : ℂ) + 1) ^ 2)

noncomputable def chi (z : ℂ) : ℂ := ∑' n, chiT z n

noncomputable def chiD (z : ℂ) (n : ℕ) : ℂ := (z ^ n + (-z) ^ n) / (2 * ((n : ℂ) + 1))

noncomputable def ell (z : ℂ) : ℂ := Complex.log (1 + z) - Complex.log (1 - z)

lemma norm_chiT_le {z : ℂ} (hz : ‖z‖ ≤ 1) (n : ℕ) : ‖chiT z n‖ ≤ 1 / ((n : ℝ) + 1) ^ 2 := by
  unfold chiT
  have h1 : ‖z ^ (n + 1) - (-z) ^ (n + 1)‖ ≤ 2 := by
    calc ‖z ^ (n + 1) - (-z) ^ (n + 1)‖ ≤ ‖z ^ (n + 1)‖ + ‖(-z) ^ (n + 1)‖ := norm_sub_le _ _
      _ = ‖z‖ ^ (n + 1) + ‖z‖ ^ (n + 1) := by rw [norm_pow, norm_pow, norm_neg]
      _ ≤ 1 + 1 := by gcongr <;> exact pow_le_one₀ (norm_nonneg _) hz
      _ = 2 := by norm_num
  have h2 : ‖(2 * ((n : ℂ) + 1) ^ 2)‖ = 2 * ((n : ℝ) + 1) ^ 2 := by
    rw [norm_mul, norm_pow]
    have : ‖((n : ℂ) + 1)‖ = (n : ℝ) + 1 := by
      rw [show ((n : ℂ) + 1) = (((n : ℝ) + 1 : ℝ) : ℂ) by push_cast; ring, Complex.norm_real,
        Real.norm_of_nonneg (by positivity)]
    rw [this]; simp
  rw [norm_div, h2, div_le_iff₀ (by positivity)]
  have : 1 / ((n : ℝ) + 1) ^ 2 * (2 * ((n : ℝ) + 1) ^ 2) = 2 := by field_simp
  rw [this]; exact h1

lemma summable_chiT {z : ℂ} (hz : ‖z‖ ≤ 1) : Summable (chiT z) :=
  Summable.of_norm_bounded D00.summable_inv_sq (norm_chiT_le hz)

lemma hasDerivAt_chiT (n : ℕ) (z : ℂ) : HasDerivAt (fun w => chiT w n) (chiD z n) z := by
  have h1 := hasDerivAt_pow (n + 1) z
  have h2 := (hasDerivAt_neg z).pow (n + 1)
  have h := (h1.sub h2).div_const (2 * ((n : ℂ) + 1) ^ 2)
  convert h using 1
  · funext w; simp [chiT]
  unfold chiD
  have : ((n : ℂ) + 1) ≠ 0 := by exact_mod_cast Nat.succ_ne_zero n
  push_cast
  field_simp
  ring

lemma hasSum_chiD {z : ℂ} (hz : ‖z‖ < 1) (h0 : z ≠ 0) : HasSum (chiD z) (ell z / (2 * z)) := by
  have S1 := Complex.hasSum_taylorSeries_neg_log' hz
  have S2 := Complex.hasSum_taylorSeries_neg_log' (z := -z) (by rwa [norm_neg])
  have S := (S1.sub S2).div_const (2 * z)
  rw [sub_neg_eq_add] at S
  convert S using 1
  · funext n
    unfold chiD
    have : ((n : ℂ) + 1) ≠ 0 := by exact_mod_cast Nat.succ_ne_zero n
    rw [pow_succ, pow_succ, ← sub_div,
      show z ^ n * z - (-z) ^ n * -z = z * (z ^ n + (-z) ^ n) by ring]
    field_simp
  · unfold ell; rw [sub_neg_eq_add]; ring

lemma norm_chiD_le {y : ℂ} {R : ℝ} (hy : ‖y‖ ≤ R) (n : ℕ) : ‖chiD y n‖ ≤ R ^ n := by
  unfold chiD
  have hR : 0 ≤ R := le_trans (norm_nonneg _) hy
  have h1 : ‖y ^ n + (-y) ^ n‖ ≤ 2 * R ^ n := by
    calc ‖y ^ n + (-y) ^ n‖ ≤ ‖y ^ n‖ + ‖(-y) ^ n‖ := norm_add_le _ _
      _ = ‖y‖ ^ n + ‖y‖ ^ n := by rw [norm_pow, norm_pow, norm_neg]
      _ ≤ R ^ n + R ^ n := by gcongr
      _ = 2 * R ^ n := by ring
  have h2 : ‖(2 * ((n : ℂ) + 1))‖ = 2 * ((n : ℝ) + 1) := by
    rw [norm_mul]
    have : ‖((n : ℂ) + 1)‖ = (n : ℝ) + 1 := by
      rw [show ((n : ℂ) + 1) = (((n : ℝ) + 1 : ℝ) : ℂ) by push_cast; ring, Complex.norm_real,
        Real.norm_of_nonneg (by positivity)]
    rw [this]; simp
  rw [norm_div, h2, div_le_iff₀ (by positivity)]
  have : 0 ≤ R ^ n := pow_nonneg hR n
  have hn : (0 : ℝ) ≤ n := Nat.cast_nonneg n
  nlinarith

lemma hasDerivAt_chi {z : ℂ} (hz : ‖z‖ < 1) : HasDerivAt chi (∑' n, chiD z n) z := by
  set R := (1 + ‖z‖) / 2 with hRdef
  have hR1 : R < 1 := by rw [hRdef]; linarith
  have hzR : ‖z‖ < R := by rw [hRdef]; linarith
  have hR0 : 0 ≤ R := le_trans (norm_nonneg _) hzR.le
  have hu : Summable (fun n : ℕ => R ^ n) := summable_geometric_of_lt_one hR0 hR1
  have hzt : z ∈ Metric.ball (0 : ℂ) R := by simpa using hzR
  have h0t : (0 : ℂ) ∈ Metric.ball (0 : ℂ) R := by simpa using lt_of_le_of_lt (norm_nonneg z) hzR
  exact hasDerivAt_tsum_of_isPreconnected (g := fun n w => chiT w n) (g' := fun n w => chiD w n)
    hu Metric.isOpen_ball (convex_ball _ _).isPreconnected (fun n y _ => hasDerivAt_chiT n y)
    (fun n y hy => norm_chiD_le (by simpa using (le_of_lt (Metric.mem_ball.mp hy))) n) h0t
    (summable_chiT (by simp)) hzt

lemma hasDerivAt_chi' {z : ℂ} (hz : ‖z‖ < 1) (h0 : z ≠ 0) :
    HasDerivAt chi (ell z / (2 * z)) z := by
  have := hasDerivAt_chi hz
  rwa [(hasSum_chiD hz h0).tsum_eq] at this

lemma continuousOn_chi : ContinuousOn chi (Metric.closedBall 0 1) := by
  refine continuousOn_tsum (u := fun n : ℕ => 1 / ((n : ℝ) + 1) ^ 2) (fun n => ?_)
    D00.summable_inv_sq (fun n x hx => norm_chiT_le (by simpa using hx) n)
  exact (continuous_iff_continuousAt.2 (fun w => (hasDerivAt_chiT n w).continuousAt)).continuousOn

/-- `χ` at a real point is the real `χ₂ = ½(Li₂ x - Li₂(-x))`. -/
lemma chi_ofReal {x : ℝ} (hx : |x| ≤ 1) :
    chi (x : ℂ) = (((D00.Li2 x - D00.Li2 (-x)) / 2 : ℝ) : ℂ) := by
  have hx' : |-x| ≤ 1 := by rwa [abs_neg]
  have hR : HasSum (fun n : ℕ => (x ^ (n + 1) - (-x) ^ (n + 1)) / (2 * ((n : ℝ) + 1) ^ 2))
      ((D00.Li2 x - D00.Li2 (-x)) / 2) := by
    have := ((D00.hasSum_Li2 hx).sub (D00.hasSum_Li2 hx')).div_const 2
    convert this using 1
    funext n
    have : ((n : ℝ) + 1) ≠ 0 := by positivity
    field_simp
  have hC := Complex.hasSum_ofReal.mpr hR
  unfold chi
  rw [← hC.tsum_eq]
  congr 1
  funext n
  simp only [chiT]
  push_cast
  ring

/-! ## 2. Landen's identity for `χ` on `U = {Re z > 0, |z| < 1}` -/

noncomputable def U : Set ℂ := {z | 0 < z.re} ∩ Metric.ball 0 1

lemma mem_U {z : ℂ} : z ∈ U ↔ 0 < z.re ∧ ‖z‖ < 1 := by
  simp [U]

lemma norm_lt_one_iff' (z : ℂ) : ‖z‖ < 1 ↔ z.re ^ 2 + z.im ^ 2 < 1 := by
  rw [← sq_lt_one_iff₀ (norm_nonneg z), Complex.sq_norm, Complex.normSq_apply]
  constructor <;> intro h <;> nlinarith

lemma isOpen_U : IsOpen U :=
  (isOpen_lt continuous_const Complex.continuous_re).inter Metric.isOpen_ball

lemma convex_U : Convex ℝ U := (convex_halfSpace_re_gt 0).inter (convex_ball 0 1)

noncomputable def wz (z : ℂ) : ℂ := (1 - z) / (1 + z)

lemma U_facts {z : ℂ} (hz : z ∈ U) : ‖z‖ < 1 ∧ z ≠ 0 ∧ 1 + z ≠ 0 ∧ 1 - z ≠ 0 ∧
    z ∈ Complex.slitPlane ∧ 1 + z ∈ Complex.slitPlane ∧ 1 - z ∈ Complex.slitPlane := by
  obtain ⟨h1, hn⟩ := mem_U.mp hz
  have h2 := (norm_lt_one_iff' z).mp hn
  have hre : z.re < 1 := by nlinarith [sq_nonneg z.im]
  refine ⟨hn, ?_, ?_, ?_, ?_, ?_, ?_⟩
  · intro h; rw [h] at h1; simp at h1
  · intro h; have := congrArg Complex.re h; simp at this; linarith
  · intro h; have := congrArg Complex.re h; simp at this; linarith
  · exact Complex.mem_slitPlane_iff.mpr (Or.inl h1)
  · exact Complex.mem_slitPlane_iff.mpr (Or.inl (by simp; linarith))
  · exact Complex.mem_slitPlane_iff.mpr (Or.inl (by simp; linarith))

lemma wz_mem {z : ℂ} (hz : z ∈ U) : wz z ∈ U := by
  obtain ⟨h1, hn⟩ := mem_U.mp hz
  obtain ⟨-, -, hp, -, -, -, -⟩ := U_facts hz
  have h2 := (norm_lt_one_iff' z).mp hn
  have hN : 0 < Complex.normSq (1 + z) := Complex.normSq_pos.mpr hp
  rw [mem_U]
  constructor
  · unfold wz
    rw [Complex.div_re, ← add_div]
    apply div_pos _ hN
    simp
    nlinarith
  · unfold wz
    rw [norm_div, div_lt_one (norm_pos_iff.mpr hp)]
    have e1 : ‖1 - z‖ ^ 2 < ‖1 + z‖ ^ 2 := by
      rw [Complex.sq_norm, Complex.sq_norm, Complex.normSq_apply, Complex.normSq_apply]
      simp; nlinarith
    by_contra h
    push Not at h
    nlinarith [norm_nonneg (1 - z), norm_nonneg (1 + z)]

lemma wz_wz {z : ℂ} (hz : 1 + z ≠ 0) : wz (wz z) = z := by
  unfold wz
  have h2 : (1 + (1 - z) / (1 + z)) = 2 / (1 + z) := by field_simp; ring
  have h3 : (1 - (1 - z) / (1 + z)) = 2 * z / (1 + z) := by field_simp; ring
  rw [h2, h3]
  field_simp

lemma hasDerivAt_wz {z : ℂ} (hz : 1 + z ≠ 0) : HasDerivAt wz (-2 / (1 + z) ^ 2) z := by
  have := ((hasDerivAt_id z).const_sub 1).div ((hasDerivAt_id z).const_add 1) hz
  convert this using 1
  · funext x; simp [wz]
  · simp only [id]
    field_simp
    ring

lemma hasDerivAt_ell {u : ℂ} (h1 : 1 + u ∈ Complex.slitPlane) (h2 : 1 - u ∈ Complex.slitPlane) :
    HasDerivAt ell (1 / (1 + u) + 1 / (1 - u)) u := by
  have e1 := ((hasDerivAt_id u).const_add 1).clog h1
  have e2 := ((hasDerivAt_id u).const_sub 1).clog h2
  have := e1.sub e2
  convert this using 1
  · funext x; simp [ell]
  · simp only [id]
    ring

lemma ell_ofReal {x : ℝ} (h : |x| < 1) :
    ell (x : ℂ) = ((Real.log (1 + x) - Real.log (1 - x) : ℝ) : ℂ) := by
  rw [abs_lt] at h
  unfold ell
  push_cast
  rw [Complex.ofReal_log (by linarith), Complex.ofReal_log (by linarith)]
  push_cast; ring

lemma b_mem_U : ((D02.b : ℝ) : ℂ) ∈ U := by
  rw [mem_U]
  simp only [Complex.ofReal_re, Complex.norm_real, Real.norm_of_nonneg D02.b_pos.le]
  exact ⟨D02.b_pos, D02.b_lt⟩

lemma b_fix : (1 - D02.b) / (1 + D02.b) = D02.b := by
  have hb := D02.b_pos
  rw [div_eq_iff (by linarith)]
  unfold D02.b
  have := D02.hs
  nlinarith

lemma wz_b : wz (D02.b : ℂ) = (D02.b : ℂ) := by
  unfold wz
  calc (1 - (D02.b : ℂ)) / (1 + (D02.b : ℂ)) = (((1 - D02.b) / (1 + D02.b) : ℝ) : ℂ) := by
        push_cast; ring
    _ = (D02.b : ℂ) := by rw [b_fix]

lemma ell_b : ell (D02.b : ℂ) = (D02.L : ℂ) := by
  rw [ell_ofReal (by rw [abs_lt]; constructor <;> linarith [D02.b_pos, D02.b_lt]),
    D02.log_1pb, D02.log_1mb]
  congr 1; ring

lemma log_b_C : Complex.log (D02.b : ℂ) = ((-D02.L : ℝ) : ℂ) := by
  rw [← Complex.ofReal_log D02.b_pos.le, D02.log_b]

noncomputable def Gf (z : ℂ) : ℂ := ell (wz z) + Complex.log z

lemma hasDerivAt_Gf {z : ℂ} (hz : z ∈ U) : HasDerivAt Gf 0 z := by
  obtain ⟨-, h0, hp, hm, hs, -, -⟩ := U_facts hz
  obtain ⟨-, -, -, -, -, hs1, hs2⟩ := U_facts (wz_mem hz)
  have e := ((hasDerivAt_ell hs1 hs2).comp z (hasDerivAt_wz hp)).add (Complex.hasDerivAt_log hs)
  convert e using 1
  · funext x; simp [Gf]
  have e1 : 1 + wz z = 2 / (1 + z) := by unfold wz; field_simp; ring
  have e2 : 1 - wz z = 2 * z / (1 + z) := by unfold wz; field_simp; ring
  rw [e1, e2]
  field_simp
  ring

lemma Gf_eq_zero {z : ℂ} (hz : z ∈ U) : Gf z = 0 := by
  have hd : DifferentiableOn ℂ Gf U := fun x hx => (hasDerivAt_Gf hx).differentiableAt.differentiableWithinAt
  have key := isOpen_U.is_const_of_deriv_eq_zero convex_U.isPreconnected hd
    (fun x hx => (hasDerivAt_Gf hx).deriv) hz b_mem_U
  rw [key]
  unfold Gf
  rw [wz_b, ell_b, log_b_C]
  push_cast; ring

lemma ell_wz {z : ℂ} (hz : z ∈ U) : ell (wz z) = -Complex.log z := by
  have := Gf_eq_zero hz
  unfold Gf at this
  linear_combination this

end S06.PE

namespace S06.PE

open Real MeasureTheory Set Filter Topology

noncomputable def Phi (z : ℂ) : ℂ := chi z + chi (wz z) - (1 / 2) * Complex.log z * ell z

lemma hasDerivAt_Phi {z : ℂ} (hz : z ∈ U) : HasDerivAt Phi 0 z := by
  obtain ⟨hn, h0, hp, hm, hs, hs1, hs2⟩ := U_facts hz
  have hw := wz_mem hz
  obtain ⟨hwn, hw0, -, -, -, -, -⟩ := U_facts hw
  have d1 := hasDerivAt_chi' hn h0
  have d2 := (hasDerivAt_chi' hwn hw0).comp z (hasDerivAt_wz hp)
  have d3 := ((Complex.hasDerivAt_log hs).const_mul (1 / 2 : ℂ)).mul (hasDerivAt_ell hs1 hs2)
  have e := (d1.add d2).sub d3
  convert e using 1
  · funext x; simp only [Phi, Function.comp_apply, Pi.add_apply, Pi.sub_apply, Pi.mul_apply]
  rw [ell_wz hz]
  have e2 : wz z = (1 - z) / (1 + z) := rfl
  rw [e2]
  field_simp
  ring

lemma Phi_b : Phi (D02.b : ℂ) = ((π ^ 2 / 8 : ℝ) : ℂ) := by
  unfold Phi
  rw [wz_b, ell_b, log_b_C, chi_ofReal (by rw [abs_le]; constructor <;> linarith [D02.b_pos, D02.b_lt]),
    D02.proof_1]
  push_cast; ring

lemma landen {z : ℂ} (hz : z ∈ U) : Phi z = ((π ^ 2 / 8 : ℝ) : ℂ) := by
  have hd : DifferentiableOn ℂ Phi U := fun x hx =>
    (hasDerivAt_Phi hx).differentiableAt.differentiableWithinAt
  rw [isOpen_U.is_const_of_deriv_eq_zero convex_U.isPreconnected hd
    (fun x hx => (hasDerivAt_Phi hx).deriv) hz b_mem_U, Phi_b]

end S06.PE

/-! ## 3. Real parts of `χ` in polar form; the Fourier value `Re χ(e^{iα})` -/

namespace S06.PE

open Real MeasureTheory Set Filter Topology

noncomputable def pol (r θ : ℝ) : ℂ :=
  (r : ℂ) * (Complex.cos (θ : ℂ) + Complex.sin (θ : ℂ) * Complex.I)

lemma norm_pol (r θ : ℝ) : ‖pol r θ‖ = |r| := by
  unfold pol
  rw [norm_mul, Complex.norm_cos_add_sin_mul_I, mul_one, Complex.norm_real, Real.norm_eq_abs]

lemma pow_pol (r θ : ℝ) (m : ℕ) : pol r θ ^ m = pol (r ^ m) (m * θ) := by
  unfold pol
  rw [mul_pow, Complex.cos_add_sin_mul_I_pow]
  push_cast; ring_nf

lemma neg_pol (r θ : ℝ) : -pol r θ = pol (-r) θ := by
  unfold pol; push_cast; ring

lemma re_pol (r θ : ℝ) : (pol r θ).re = r * Real.cos θ := by
  unfold pol
  rw [Complex.re_ofReal_mul]
  simp [← Complex.ofReal_cos, ← Complex.ofReal_sin]

lemma re_chiT_pol (r θ : ℝ) (n : ℕ) :
    (chiT (pol r θ) n).re
      = (r ^ (n + 1) - (-r) ^ (n + 1)) * Real.cos ((n + 1 : ℝ) * θ) / (2 * ((n : ℝ) + 1) ^ 2) := by
  unfold chiT
  rw [neg_pol, pow_pol, pow_pol]
  have : (2 * ((n : ℂ) + 1) ^ 2) = ((2 * ((n : ℝ) + 1) ^ 2 : ℝ) : ℂ) := by push_cast; ring
  rw [this, Complex.div_ofReal_re, Complex.sub_re, re_pol, re_pol]
  push_cast
  ring

lemma hasSum_re_chi_pol {r : ℝ} (θ : ℝ) (hr : |r| ≤ 1) :
    HasSum (fun n : ℕ => (r ^ (n + 1) - (-r) ^ (n + 1)) * Real.cos ((n + 1 : ℝ) * θ)
      / (2 * ((n : ℝ) + 1) ^ 2)) (chi (pol r θ)).re := by
  have h := Complex.hasSum_re (summable_chiT (z := pol r θ) (by rw [norm_pol]; exact hr)).hasSum
  simp only [re_chiT_pol] at h
  exact h

/-- `Σ_{m≥1} cos(mθ)/m² = θ²/4 - πθ/2 + π²/6` on `[0, 2π]`. -/
lemma hasSum_cos_sq {θ : ℝ} (h0 : 0 ≤ θ) (h1 : θ ≤ 2 * π) :
    HasSum (fun n : ℕ => Real.cos ((n + 1 : ℝ) * θ) / ((n : ℝ) + 1) ^ 2)
      (θ ^ 2 / 4 - π * θ / 2 + π ^ 2 / 6) := by
  have hpi := Real.pi_pos
  have hx : θ / (2 * π) ∈ Icc (0 : ℝ) 1 := by
    constructor
    · positivity
    · rw [div_le_one (by positivity)]; exact h1
  have h := hasSum_one_div_nat_pow_mul_cos (k := 1) one_ne_zero hx
  rw [← hasSum_nat_add_iff' 1] at h
  simp only [Finset.range_one, Finset.sum_singleton, Nat.cast_zero] at h
  change HasSum _ (_ * bernoulliFun 2 _ - _) at h
  rw [bernoulliFun_two] at h
  convert h using 1
  · funext n
    push_cast
    have : 2 * π * ((n : ℝ) + 1) * (θ / (2 * π)) = ((n : ℝ) + 1) * θ := by field_simp
    rw [this]; ring
  · simp
    field_simp
    ring

end S06.PE

/-! ## 4. The point values -/

namespace S06.PE

open Real MeasureTheory Set Filter Topology

lemma alpha_pos : 0 < S02.alphac := by
  unfold S02.alphac; exact Real.arctan_pos.mpr (by positivity)

lemma alpha_lt : S02.alphac < π / 2 := by
  unfold S02.alphac; exact Real.arctan_lt_pi_div_two _

lemma re_chi_omega : (chi (pol 1 S02.alphac)).re = π ^ 2 / 8 - π * S02.alphac / 4 := by
  have hpi := Real.pi_pos
  have ha0 := alpha_pos
  have ha1 := alpha_lt
  have h := hasSum_re_chi_pol S02.alphac (r := 1) (by norm_num)
  have F1 := hasSum_cos_sq (θ := S02.alphac) ha0.le (by linarith)
  have F2 := hasSum_cos_sq (θ := S02.alphac + π) (by linarith) (by linarith)
  have F := (F1.sub F2).div_const 2
  have F' : HasSum (fun n : ℕ => (1 ^ (n + 1) - (-1) ^ (n + 1)) * Real.cos (((n : ℝ) + 1) * S02.alphac)
      / (2 * ((n : ℝ) + 1) ^ 2)) (π ^ 2 / 8 - π * S02.alphac / 4) := by
    convert F using 1
    · funext n
      have hc : Real.cos (((n : ℝ) + 1) * (S02.alphac + π))
          = (-1) ^ (n + 1) * Real.cos (((n : ℝ) + 1) * S02.alphac) := by
        rw [show ((n : ℝ) + 1) * (S02.alphac + π) = ((n : ℝ) + 1) * S02.alphac + ((n + 1 : ℕ) : ℝ) * π by
          push_cast; ring, Real.cos_add_nat_mul_pi]
      rw [hc]
      have : ((n : ℝ) + 1) ≠ 0 := by positivity
      field_simp
      ring
    · ring
  exact h.unique F'


/-- The even part: `Σ_m ½(x^m + (-x)^m) cos(mπ/4)/m² = Li₂(-x⁴)/16`. -/
noncomputable def evT (x : ℝ) (n : ℕ) : ℝ :=
  (x ^ (n + 1) + (-x) ^ (n + 1)) * Real.cos (((n : ℝ) + 1) * (π / 4)) / (2 * ((n : ℝ) + 1) ^ 2)

lemma evT_zero (x : ℝ) (n : ℕ) (hn : n ∉ Set.range (fun j : ℕ => 4 * j + 3)) : evT x n = 0 := by
  unfold evT
  have hdm := Nat.div_add_mod n 4
  have hlt := Nat.mod_lt n (show 4 > 0 by norm_num)
  set k := n / 4
  interval_cases h : n % 4
  · have ho : Odd (n + 1) := ⟨2 * k, by omega⟩
    rw [ho.neg_pow]; simp
  · have hc : Real.cos (((n : ℝ) + 1) * (π / 4)) = 0 := by
      have : ((n : ℝ) + 1) * (π / 4) = π / 2 + (k : ℕ) * π := by
        have : (n : ℝ) = 4 * k + 1 := by exact_mod_cast (by omega : n = 4 * k + 1)
        rw [this]; ring
      rw [this, Real.cos_add_nat_mul_pi, Real.cos_pi_div_two, mul_zero]
    rw [hc]; simp
  · have ho : Odd (n + 1) := ⟨2 * k + 1, by omega⟩
    rw [ho.neg_pow]; simp
  · exact absurd ⟨k, by omega⟩ hn

lemma evT_g (x : ℝ) (j : ℕ) :
    evT x (4 * j + 3) = (-x ^ 4) ^ (j + 1) / ((j : ℝ) + 1) ^ 2 / 16 := by
  unfold evT
  have he : Even (4 * j + 3 + 1) := ⟨2 * j + 2, by omega⟩
  rw [he.neg_pow]
  have hc : Real.cos ((((4 * j + 3 : ℕ) : ℝ) + 1) * (π / 4)) = (-1) ^ (j + 1) := by
    rw [show (((4 * j + 3 : ℕ) : ℝ) + 1) * (π / 4) = ((j + 1 : ℕ) : ℝ) * π by push_cast; ring,
      Real.cos_nat_mul_pi]
  rw [hc]
  have : ((j : ℝ) + 1) ≠ 0 := by positivity
  rw [show 4 * j + 3 + 1 = 4 * (j + 1) by ring, pow_mul, neg_pow (x ^ 4)]
  push_cast
  field_simp
  ring

lemma hasSum_evT {x : ℝ} (hx : |x| ≤ 1) : HasSum (evT x) (D00.Li2 (-x ^ 4) / 16) := by
  have hx4 : |-x ^ 4| ≤ 1 := by
    rw [abs_neg, abs_pow]; exact pow_le_one₀ (abs_nonneg _) hx
  have h := (D00.hasSum_Li2 hx4).div_const 16
  have hinj : Function.Injective (fun j : ℕ => 4 * j + 3) := fun a b h => by simp at h; omega
  rw [← hinj.hasSum_iff (evT_zero x)]
  convert h using 1
  funext j
  simp only [Function.comp_apply]
  rw [evT_g]

lemma re_chi_y : (chi (pol D02.b (-(π / 4)))).re
    = (1 / 2) * D00.Li2 D02.b + π ^ 2 / 192 - D00.Li2 (-D02.b ^ 4) / 16 := by
  have hb : |D02.b| ≤ 1 := by rw [abs_le]; constructor <;> linarith [D02.b_pos, D02.b_lt]
  have h := hasSum_re_chi_pol (-(π / 4)) hb
  have F := S02.lemmaA_pi4.sub (hasSum_evT hb)
  refine h.unique ?_
  convert F using 1
  funext n
  unfold evT
  rw [mul_neg, Real.cos_neg]
  have : ((n : ℝ) + 1) ≠ 0 := by positivity
  field_simp
  ring


lemma im_pol (r θ : ℝ) : (pol r θ).im = r * Real.sin θ := by
  unfold pol
  rw [Complex.im_ofReal_mul]
  simp [← Complex.ofReal_cos, ← Complex.ofReal_sin]

lemma sqrt32 : Real.sqrt (1 + (1 / Real.sqrt 2) ^ 2) = Real.sqrt 3 / Real.sqrt 2 := by
  rw [div_pow, Real.sq_sqrt (by norm_num), ← Real.sqrt_div (by norm_num)]
  norm_num

lemma cos_alpha : Real.cos S02.alphac = Real.sqrt 2 / Real.sqrt 3 := by
  unfold S02.alphac
  rw [Real.cos_arctan, sqrt32]
  have : Real.sqrt 3 ≠ 0 := by positivity
  have : Real.sqrt 2 ≠ 0 := by positivity
  field_simp

lemma sin_alpha : Real.sin S02.alphac = 1 / Real.sqrt 3 := by
  unfold S02.alphac
  rw [Real.sin_arctan, sqrt32]
  have : Real.sqrt 3 ≠ 0 := by positivity
  have : Real.sqrt 2 ≠ 0 := by positivity
  field_simp

lemma s3sq : Real.sqrt 3 * Real.sqrt 3 = 3 := Real.mul_self_sqrt (by norm_num)

/-- `t = ω/√3`. -/
noncomputable def tP : ℂ := pol (1 / Real.sqrt 3) S02.alphac

/-- `y = b e^{-iπ/4}`. -/
noncomputable def yP : ℂ := pol D02.b (-(π / 4))

lemma tP_re : tP.re = Real.sqrt 2 / 3 := by
  unfold tP; rw [re_pol, cos_alpha]
  have : Real.sqrt 3 ≠ 0 := by positivity
  field_simp; rw [Real.sq_sqrt (by norm_num)]

lemma tP_im : tP.im = 1 / 3 := by
  unfold tP; rw [im_pol, sin_alpha]
  have : Real.sqrt 3 ≠ 0 := by positivity
  field_simp; rw [Real.sq_sqrt (by norm_num)]

lemma yP_re : yP.re = D02.b * (Real.sqrt 2 / 2) := by
  unfold yP; rw [re_pol, Real.cos_neg, Real.cos_pi_div_four]

lemma yP_im : yP.im = -(D02.b * (Real.sqrt 2 / 2)) := by
  unfold yP; rw [im_pol, Real.sin_neg, Real.sin_pi_div_four]; ring

lemma tP_mem : tP ∈ U := by
  rw [mem_U]
  constructor
  · rw [tP_re]; positivity
  · rw [norm_lt_one_iff', tP_re, tP_im, div_pow, Real.sq_sqrt (by norm_num)]; norm_num

lemma yP_mem : yP ∈ U := by
  rw [mem_U]
  constructor
  · rw [yP_re]; have := D02.b_pos; positivity
  · unfold yP; rw [norm_pol, abs_of_pos D02.b_pos]; exact D02.b_lt

lemma wz_tP : wz tP = yP := by
  obtain ⟨-, -, hp, -, -, -, -⟩ := U_facts tP_mem
  unfold wz
  rw [div_eq_iff hp]
  have h2 := D02.hs
  apply Complex.ext
  · simp only [Complex.sub_re, Complex.one_re, Complex.mul_re, Complex.add_re, Complex.add_im,
      Complex.one_im, tP_re, tP_im, yP_re, yP_im]
    unfold D02.b
    linear_combination (-(1:ℝ) / 2 - Real.sqrt 2 / 6) * h2
  · simp only [Complex.sub_im, Complex.one_im, Complex.mul_im, Complex.add_re, Complex.add_im,
      Complex.one_re, tP_re, tP_im, yP_re, yP_im]
    unfold D02.b
    linear_combination ((1:ℝ) / 6 + Real.sqrt 2 / 6) * h2

lemma wz_yP : wz yP = tP := by
  obtain ⟨-, -, hp, -, -, -, -⟩ := U_facts tP_mem
  rw [← wz_tP, wz_wz hp]

lemma log_pol {r θ : ℝ} (hr : 0 < r) (hθ : θ ∈ Ioc (-π) π) :
    Complex.log (pol r θ) = ((Real.log r : ℝ) : ℂ) + (θ : ℂ) * Complex.I := by
  apply Complex.ext
  · rw [Complex.log_re, norm_pol, abs_of_pos hr]; simp
  · rw [Complex.log_im]
    unfold pol
    rw [Complex.arg_real_mul _ hr, Complex.arg_cos_add_sin_mul_I hθ]
    simp

lemma log_tP : Complex.log tP = ((-(Real.log 3 / 2) : ℝ) : ℂ) + (S02.alphac : ℂ) * Complex.I := by
  have hpi := Real.pi_pos
  unfold tP
  rw [log_pol (by positivity) ⟨by linarith [alpha_pos], by linarith [alpha_lt]⟩]
  congr 2
  rw [one_div, Real.log_inv, Real.log_sqrt (by norm_num)]

lemma log_yP : Complex.log yP = ((-D02.L : ℝ) : ℂ) + ((-(π / 4) : ℝ) : ℂ) * Complex.I := by
  have hpi := Real.pi_pos
  unfold yP
  rw [log_pol D02.b_pos ⟨by linarith, by linarith⟩, D02.log_b]

lemma ell_tP : ell tP = ((D02.L : ℝ) : ℂ) + ((π / 4 : ℝ) : ℂ) * Complex.I := by
  have := ell_wz yP_mem
  rw [wz_yP, log_yP] at this
  rw [this]; push_cast; ring

/-- Landen at `t`, real part. -/
lemma re_chi_tP : (chi tP).re = π ^ 2 / 8 - (chi yP).re - D02.L * Real.log 3 / 4
    - π * S02.alphac / 8 := by
  have h := congrArg Complex.re (landen tP_mem)
  unfold Phi at h
  rw [wz_tP, log_tP, ell_tP] at h
  simp only [Complex.sub_re, Complex.add_re, Complex.mul_re, Complex.ofReal_re, Complex.ofReal_im,
    Complex.I_re, Complex.I_im, Complex.add_im, Complex.mul_im] at h
  norm_num at h
  linarith


end S06.PE

/-! ## 5. The real function `T(x) = Re χ(x e^{iα})` and the substitution `x = x(s)` -/

namespace S06.PE

open Real MeasureTheory Set Filter Topology

noncomputable def T (x : ℝ) : ℝ := (chi (pol x S02.alphac)).re

lemma pol_eq_mul (x : ℝ) : pol x S02.alphac = (x : ℂ) * pol 1 S02.alphac := by
  unfold pol; push_cast; ring

lemma norm_one_pol : ‖pol 1 S02.alphac‖ = 1 := by rw [norm_pol]; simp

lemma hasDerivAt_T {x : ℝ} (h0 : 0 < x) (h1 : x < 1) :
    HasDerivAt T ((Real.log ‖1 + pol x S02.alphac‖ - Real.log ‖1 - pol x S02.alphac‖) / (2 * x)) x := by
  set ω := pol 1 S02.alphac with hω
  have hω0 : ω ≠ 0 := by
    intro h; have := norm_one_pol; rw [← hω, h, norm_zero] at this; norm_num at this
  have hz : ‖(x : ℂ) * ω‖ < 1 := by
    rw [norm_mul, norm_one_pol, mul_one, Complex.norm_real, Real.norm_of_nonneg h0.le]; exact h1
  have hne : (x : ℂ) * ω ≠ 0 := mul_ne_zero (by exact_mod_cast h0.ne') hω0
  have d := ((hasDerivAt_chi' hz hne).comp (x : ℂ) (hasDerivAt_mul_const ω)).real_of_complex
  convert d using 1
  · funext y; simp only [T, Function.comp_apply]; rw [pol_eq_mul]
  · have e : ell ((x : ℂ) * ω) / (2 * ((x : ℂ) * ω)) * ω = ell ((x : ℂ) * ω) / ((2 * x : ℝ) : ℂ) := by
      have : (x : ℂ) ≠ 0 := by exact_mod_cast h0.ne'
      push_cast
      field_simp
    rw [e, Complex.div_ofReal_re, ← pol_eq_mul]
    unfold ell
    rw [Complex.sub_re, Complex.log_re, Complex.log_re]

/-- `S(s) = √(s² - s + 1)`. -/
noncomputable def Sq (s : ℝ) : ℝ := Real.sqrt (s ^ 2 - s + 1)

/-- `x(s) = (2s - 1 + 2S(s))/√3`. -/
noncomputable def xs (s : ℝ) : ℝ := (2 * s - 1 + 2 * Sq s) / Real.sqrt 3

lemma Sq_arg_pos (s : ℝ) : 0 < s ^ 2 - s + 1 := by nlinarith [sq_nonneg (s - 1 / 2)]

lemma Sq_pos (s : ℝ) : 0 < Sq s := Real.sqrt_pos.mpr (Sq_arg_pos s)

lemma Sq_sq (s : ℝ) : Sq s ^ 2 = s ^ 2 - s + 1 := Real.sq_sqrt (Sq_arg_pos s).le

lemma qs_sq (s : ℝ) : S02.qs s ^ 2 = 2 * Sq s ^ 2 - 1 := by
  unfold S02.qs; rw [Real.sq_sqrt (by positivity), Sq_sq]; ring

lemma beta_eq (s : ℝ) : Real.sqrt (1 + S02.qs s ^ 2) = Real.sqrt 2 * Sq s := by
  rw [qs_sq, show 1 + (2 * Sq s ^ 2 - 1) = 2 * Sq s ^ 2 by ring, Real.sqrt_mul (by norm_num),
    Real.sqrt_sq (Sq_pos s).le]

lemma xs_key (s : ℝ) : xs s ^ 2 + 1 = 4 * Sq s / Real.sqrt 3 * xs s := by
  unfold xs
  have h3 : Real.sqrt 3 ^ 2 = 3 := Real.sq_sqrt (by norm_num)
  have h3' : Real.sqrt 3 ≠ 0 := by positivity
  have hS := Sq_sq s
  field_simp
  rw [h3]
  linear_combination (-4) * hS

lemma xs_pos {s : ℝ} : 0 < xs s := by
  unfold xs
  apply div_pos _ (by positivity)
  have hS := Sq_sq s
  have h0 := Sq_pos s
  nlinarith [sq_nonneg (2 * s - 1)]

lemma xs_lt_one {s : ℝ} (h : s < 1 / 2) : xs s < 1 := by
  unfold xs
  rw [div_lt_one (by positivity)]
  have hS := Sq_sq s
  have h0 := Sq_pos s
  have h3 : Real.sqrt 3 ^ 2 = 3 := Real.sq_sqrt (by norm_num)
  have h3p : 0 < Real.sqrt 3 := by positivity
  -- 2S < √3 + 1 - 2s
  have hpos : 0 < Real.sqrt 3 + 1 - 2 * s := by linarith
  have : (2 * Sq s) ^ 2 < (Real.sqrt 3 + 1 - 2 * s) ^ 2 := by nlinarith
  have h2S : 2 * Sq s < Real.sqrt 3 + 1 - 2 * s := by
    by_contra hc; push Not at hc; nlinarith
  linarith

lemma xs_half : xs (1 / 2) = 1 := by
  unfold xs Sq
  rw [show (1 / 2 : ℝ) ^ 2 - 1 / 2 + 1 = 3 / 4 by norm_num, Real.sqrt_div' _ (by norm_num : (0:ℝ) ≤ 4),
    show (4 : ℝ) = 2 ^ 2 by norm_num, Real.sqrt_sq (by norm_num)]
  have h3' : Real.sqrt 3 ≠ 0 := by positivity
  field_simp
  ring

lemma xs_zero : xs 0 = 1 / Real.sqrt 3 := by
  unfold xs Sq; norm_num

lemma hasDerivAt_xs (s : ℝ) :
    HasDerivAt xs ((2 + (2 * s - 1) / Sq s) / Real.sqrt 3) s := by
  have h1 : HasDerivAt (fun s : ℝ => s ^ 2 - s + 1) (2 * s - 1) s := by
    have := ((hasDerivAt_pow 2 s).sub (hasDerivAt_id s)).add_const 1
    convert this using 1 <;> norm_num
  have h2 := h1.sqrt (Sq_arg_pos s).ne'
  have h3 := ((((hasDerivAt_id s).const_mul 2).sub_const 1).add (h2.const_mul 2)).div_const (Real.sqrt 3)
  convert h3 using 1
  · funext y; simp [xs, Sq]
  · unfold Sq
    have := (Sq_pos s).ne'
    unfold Sq at this
    field_simp

end S06.PE

/-! ## 6. The derivative of `s ↦ √2 T(x(s))` is the integrand of `P_{2c}` -/

namespace S06.PE

open Real MeasureTheory Set Filter Topology

lemma nsq_add (x : ℝ) : ‖1 + pol x S02.alphac‖ ^ 2 = 1 + 2 * x * Real.cos S02.alphac + x ^ 2 := by
  rw [Complex.sq_norm, Complex.normSq_apply]
  simp only [Complex.add_re, Complex.add_im, Complex.one_re, Complex.one_im, re_pol, im_pol]
  have := Real.sin_sq_add_cos_sq S02.alphac
  linear_combination x ^ 2 * this

lemma nsq_sub (x : ℝ) : ‖1 - pol x S02.alphac‖ ^ 2 = 1 - 2 * x * Real.cos S02.alphac + x ^ 2 := by
  rw [Complex.sq_norm, Complex.normSq_apply]
  simp only [Complex.sub_re, Complex.sub_im, Complex.one_re, Complex.one_im, re_pol, im_pol]
  have := Real.sin_sq_add_cos_sq S02.alphac
  linear_combination x ^ 2 * this

lemma alg (x q β c : ℝ) (hβ : β ^ 2 = 1 + q ^ 2)
    (key : x ^ 2 + 1 = 2 * c * β * x) :
    (1 + 2 * x * c + x ^ 2) * q ^ 2 = (1 + β) ^ 2 * (1 - 2 * x * c + x ^ 2) := by
  linear_combination (q ^ 2 - (1 + β) ^ 2) * key + (-2 * c * x * (1 + β)) * hβ

lemma cos_alpha_lt : Real.cos S02.alphac < 1 := by
  rw [cos_alpha, div_lt_one (by positivity)]
  exact Real.sqrt_lt_sqrt (by norm_num) (by norm_num)

lemma log_diff (s : ℝ) :
    Real.log ‖1 + pol (xs s) S02.alphac‖ - Real.log ‖1 - pol (xs s) S02.alphac‖
      = Real.arsinh (1 / S02.qs s) := by
  set x := xs s with hx
  set c := Real.cos S02.alphac with hc
  set q := S02.qs s with hq
  set β := Real.sqrt (1 + q ^ 2) with hβdef
  have q0 : 0 < q := S02.qs_pos s
  have hβ : β ^ 2 = 1 + q ^ 2 := Real.sq_sqrt (by positivity)
  have β0 : 0 ≤ β := Real.sqrt_nonneg _
  have hc1 := cos_alpha_lt
  have hc0 : 0 < c := by rw [hc, cos_alpha]; positivity
  have key : x ^ 2 + 1 = 2 * c * β * x := by
    rw [hx, xs_key, hβdef, hq, beta_eq, hc, cos_alpha]
    have h2 : Real.sqrt 2 * Real.sqrt 2 = 2 := D02.hs
    have h3' : Real.sqrt 3 ≠ 0 := by positivity
    field_simp
    linear_combination (-(2:ℝ)) * Sq s * xs s * h2
  have hP := nsq_add x
  have hM := nsq_sub x
  rw [← hc] at hP hM
  have hA := alg x q β c hβ key
  have Np : 0 < ‖1 + pol x S02.alphac‖ := by
    have : 0 < ‖1 + pol x S02.alphac‖ ^ 2 := by rw [hP]; nlinarith [sq_nonneg (x + c)]
    exact lt_of_le_of_ne (norm_nonneg _) (by intro h; rw [← h] at this; norm_num at this)
  have Nm : 0 < ‖1 - pol x S02.alphac‖ := by
    have : 0 < ‖1 - pol x S02.alphac‖ ^ 2 := by rw [hM]; nlinarith [sq_nonneg (x - c)]
    exact lt_of_le_of_ne (norm_nonneg _) (by intro h; rw [← h] at this; norm_num at this)
  have heq : ‖1 + pol x S02.alphac‖ * q = (1 + β) * ‖1 - pol x S02.alphac‖ := by
    have hsq : (‖1 + pol x S02.alphac‖ * q) ^ 2 = ((1 + β) * ‖1 - pol x S02.alphac‖) ^ 2 := by
      rw [mul_pow, mul_pow, hP, hM]; linarith [hA]
    exact (sq_eq_sq₀ (by positivity) (by positivity)).mp hsq
  have hs1 : Real.sqrt (1 + (1 / q) ^ 2) = β / q := by
    have e : (β / q) ^ 2 = 1 + (1 / q) ^ 2 := by rw [div_pow, hβ]; field_simp; ring
    rw [← e, Real.sqrt_sq (by positivity)]
  unfold Real.arsinh
  rw [hs1, ← Real.log_div Np.ne' Nm.ne']
  congr 1
  rw [div_eq_iff Nm.ne']
  field_simp
  linarith [heq]

lemma deriv_val (s D : ℝ) :
    Real.sqrt 2 * (D / (2 * xs s) * ((2 + (2 * s - 1) / Sq s) / Real.sqrt 3))
      = D / (Real.sqrt 2 * Sq s) := by
  have h2 : Real.sqrt 2 * Real.sqrt 2 = 2 := D02.hs
  have hS := (Sq_pos s).ne'
  have hx := (xs_pos (s := s)).ne'
  have h3' : Real.sqrt 3 ≠ 0 := by positivity
  have h2' : Real.sqrt 2 ≠ 0 := by positivity
  have hxs : xs s * Real.sqrt 3 = 2 * s - 1 + 2 * Sq s := by unfold xs; field_simp
  field_simp
  linear_combination (-2 * D) * hxs + D * (2 * Sq s + (2 * s - 1)) * h2

end S06.PE

/-! ## 7. The fundamental theorem of calculus and the closed form -/

namespace S06.PE

open Real MeasureTheory Set Filter Topology

noncomputable def phi (s : ℝ) : ℝ := Real.arsinh (1 / S02.qs s) / Real.sqrt (1 + S02.qs s ^ 2)

noncomputable def H (s : ℝ) : ℝ := Real.sqrt 2 * T (xs s)

lemma hasDerivAt_H {s : ℝ} (h1 : s < 1 / 2) : HasDerivAt H (phi s) s := by
  have hT := hasDerivAt_T (xs_pos (s := s)) (xs_lt_one h1)
  have d := (hT.comp s (hasDerivAt_xs s)).const_mul (Real.sqrt 2)
  convert d using 1
  · funext y; simp [H]
  rw [log_diff, deriv_val]
  unfold phi; rw [beta_eq]

lemma continuous_pol : Continuous (fun x : ℝ => pol x S02.alphac) := by
  unfold pol; fun_prop

lemma continuousOn_T : ContinuousOn T (Icc 0 1) := by
  have hm : MapsTo (fun x : ℝ => pol x S02.alphac) (Icc 0 1) (Metric.closedBall 0 1) := by
    intro x hx
    rw [Metric.mem_closedBall, dist_zero_right, norm_pol, abs_le]
    constructor <;> linarith [hx.1, hx.2]
  exact Complex.continuous_re.comp_continuousOn (continuousOn_chi.comp continuous_pol.continuousOn hm)

lemma continuous_xs : Continuous xs := by unfold xs Sq; fun_prop

lemma continuousOn_H : ContinuousOn H (Icc 0 (1 / 2)) := by
  have hm : MapsTo xs (Icc 0 (1 / 2)) (Icc 0 1) := by
    intro s hs
    refine ⟨xs_pos.le, ?_⟩
    rcases eq_or_lt_of_le hs.2 with h | h
    · rw [h, xs_half]
    · exact (xs_lt_one h).le
  exact continuousOn_const.mul (continuousOn_T.comp continuous_xs.continuousOn hm)

lemma ii_phi_half : IntervalIntegrable phi volume 0 (1 / 2) := by
  refine S02.ii_2c.mono_set ?_
  rw [uIcc_of_le (by norm_num), uIcc_of_le (by norm_num)]
  exact Icc_subset_Icc le_rfl (by norm_num)

lemma ii_phi_half' : IntervalIntegrable phi volume (1 / 2) 1 := by
  refine S02.ii_2c.mono_set ?_
  rw [uIcc_of_le (by norm_num), uIcc_of_le (by norm_num)]
  exact Icc_subset_Icc (by norm_num) le_rfl

lemma int_half : ∫ s in (0 : ℝ)..1 / 2, phi s = H (1 / 2) - H 0 :=
  intervalIntegral.integral_eq_sub_of_hasDerivAt_of_le (by norm_num) continuousOn_H
    (fun s hs => hasDerivAt_H hs.2) ii_phi_half

lemma qs_symm (s : ℝ) : S02.qs (1 - s) = S02.qs s := by
  unfold S02.qs; congr 1; ring

lemma phi_symm (s : ℝ) : phi (1 - s) = phi s := by
  unfold phi; rw [qs_symm]

lemma P2c_eq : S02.P2c = 2 * ∫ s in (0 : ℝ)..1 / 2, phi s := by
  have e : S02.P2c = ∫ s in (0 : ℝ)..1, phi s := rfl
  rw [e, ← intervalIntegral.integral_add_adjacent_intervals ii_phi_half ii_phi_half']
  have h2 : ∫ s in (1 / 2 : ℝ)..1, phi s = ∫ s in (0 : ℝ)..1 / 2, phi s := by
    have := intervalIntegral.integral_comp_sub_left phi (1 : ℝ) (a := 1 / 2) (b := 1)
    simp only [phi_symm] at this
    rw [this]; norm_num
  rw [h2]; ring

lemma P2c_closed : S02.P2c = Real.sqrt 2 * (D00.Li2 D02.b - (1 / 8) * D00.Li2 (-D02.b ^ 4)
      + D02.L * Real.log 3 / 2 + π ^ 2 / 96 - π * S02.alphac / 4) := by
  rw [P2c_eq, int_half]
  unfold H
  rw [xs_half, xs_zero]
  have h1 : T 1 = π ^ 2 / 8 - π * S02.alphac / 4 := re_chi_omega
  have h2 : T (1 / Real.sqrt 3) = (chi tP).re := rfl
  have h3 : (chi yP).re = (1 / 2) * D00.Li2 D02.b + π ^ 2 / 192 - D00.Li2 (-D02.b ^ 4) / 16 :=
    re_chi_y
  rw [h1, h2, re_chi_tP, h3]
  ring

theorem main : S06.stmt_8 := ⟨S02.ii_2c, P2c_closed⟩

end S06.PE

theorem S06.proof_1 : S06.stmt_1 := S06.PA.main

theorem S06.proof_2 : S06.stmt_2 := S06.PB.main S06.PA.main

theorem S06.proof_4 : S06.stmt_4 :=
  ⟨S02.ii_f1, S02.ii_2a, S02.ii_2c, S02.Q0_P_of_J S06.proof_2.2.2⟩

theorem S06.proof_6 : S06.stmt_6 := S06.PC.main

theorem S06.proof_7 : S06.stmt_7 := S06.PD.main

theorem S06.proof_8 : S06.stmt_8 := S06.PE.main

theorem S06.proof_9 : S06.stmt_9 :=
  S02.Q0_closed_of S06.proof_4.2.2.2 S06.proof_6.2 S06.proof_7.2 S06.proof_8.2

theorem S06.proof_10 : S06.stmt_10 := S02.s2_of_Q0 S06.proof_9

theorem S06.proof_11 : S06.stmt_11 := by
  rw [S06.stmt_11, S06.proof_10]; exact S02.s2closed_bounds
