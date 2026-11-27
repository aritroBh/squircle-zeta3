import Mathlib
import Zeta3.Wave10.S05_Hexagon
/-!
# S07 -- continuation of S05: the formula of `g4:lem:precise` and the five sentences of `thm:hexagon`

Target id: S07, the continuation of `Zeta3.Wave10.S05_Hexagon` (planner split, 2026-11-26).
Paper label: `thm:hexagon` (main.tex 1709-1725, Theorem "P'''(1) of every hexagon"); second sentence in
the precise form of Lemma `g4:lem:precise` (sec_g4_hexagon.tex 39-52).
`stmt_10` and `stmt_14` keep their S05 numbers; the notation (`S05.Hex`, `S05.Ee`, `S05.closedForm`, …)
and `S05.stmt_9, 11, 12, 13` are those of S05.
-/
-- BEGIN STATEMENTS
namespace S07

/-- Lemma g4:lem:precise (second sentence of thm:hexagon): there are integer polynomials
`c_4, c_5, c_6, c_7, c_10, c_11` in `B, r, q_0, V_{02}, L², l_0, l_1, l_B` such that for every
hexagon, all weights and every edge (so every `r > 0`), `L⁵E_e` is the eleven-term combination
with the five printed blocks `c_1, c_2, c_3, c_8, c_9`. -/
def stmt_10 : Prop :=
  ∃ c4 c5 c6 c7 c10 c11 : MvPolynomial (Fin 8) ℤ,
    ∀ H : S05.Hex, S05.IsHex H → ∀ k : Fin 3,
      S05.Len H k ^ 5 * S05.Ee H k = S05.closedForm H k c4 c5 c6 c7 c10 c11

/-- Theorem thm:hexagon, all five sentences as formalised: sentence 1 (stmt_9), sentence 2 in
the precise form of g4:lem:precise (stmt_10), sentence 3 (stmt_11), sentence 4 (stmt_12),
sentence 5 at `r = 1`, `w = 1` (stmt_13). -/
def stmt_14 : Prop := S05.stmt_9 ∧ stmt_10 ∧ S05.stmt_11 ∧ S05.stmt_12 ∧ S05.stmt_13

end S07
-- END STATEMENTS

/-!
# S07 Part A: `S07.stmt_10` (the closed form of `E_e` for every `r > 0`)

Plan (all inside `S07.PA`):
1. generating functions on `|x| < 1` and the Beta values `𝓘_N(a,b)` for `a + b ≤ 2`;
2. the 20 non-local `J(n;a,b,c)(x)`, `n ≤ 1`, on `0 < |x| < 1` (series), in the basis
   `Li₃(x), Li₃(x/(x-1)), Li₂(x), 𝓛₁ = -log(1-x), π², ζ(3)`;
3. the same on `x < 1/2` by `t ↦ 1-t` (Landen), hence on every `x < 1`, `x ≠ 0`;
4. the general-`r` third `p`-derivative of the length element, the exact `L`-part, the `h²`-part
   as 40 monomials, and the assembly against the table.
-/

set_option linter.unusedSimpArgs false
set_option linter.unusedTactic false
set_option linter.unreachableTactic false
set_option linter.style.multiGoal false
set_option linter.unusedVariables false

namespace S07.PA

open MeasureTheory Set Filter Real Interval Topology Finset

/-! ## 1. generating functions -/

/-- `Li₃(x/(x-1))`. -/
noncomputable def L3p (x : ℝ) : ℝ := Z03.Li3 (x / (x - 1))

/-- the Nielsen--Kölbig value of `A(x)`. -/
noncomputable def Avl (x : ℝ) : ℝ :=
  2 * Z03.Li3 x + L3p x + D00.Li2 x * S05.cL1 x + S05.cL1 x ^ 3 / 6

lemma hs_A {x : ℝ} (hx : |x| < 1) : HasSum (Z03.Aterm x) (Avl x) := by
  have h := Z03.proof_11 x (by linarith [neg_abs_le x]) (by linarith [le_abs_self x])
  convert h using 1
  unfold Avl L3p S05.cL1; ring

lemma A_eq {x : ℝ} (hx : |x| < 1) : Z03.A x = Avl x := (hs_A hx).tsum_eq

noncomputable def g1 (q : ℕ) : ℝ := 1 / ((q : ℝ) + 1)
noncomputable def g2 (q : ℕ) : ℝ := 1 / ((q : ℝ) + 1) ^ 2
noncomputable def g3 (q : ℕ) : ℝ := 1 / ((q : ℝ) + 1) ^ 3
noncomputable def gH1 (q : ℕ) : ℝ := (harmonic (q + 1) : ℝ) / ((q : ℝ) + 1)
noncomputable def gH2 (q : ℕ) : ℝ := (harmonic (q + 1) : ℝ) / ((q : ℝ) + 1) ^ 2
noncomputable def gHH (q : ℕ) : ℝ := (harmonic (q + 1) : ℝ) ^ 2 / ((q : ℝ) + 1)
noncomputable def gHt (q : ℕ) : ℝ := S05.H2 (q + 1) / ((q : ℝ) + 1)

lemma f_g1 {x : ℝ} (hx : |x| < 1) : HasSum (fun q => g1 q * x ^ (q + 1)) (S05.cL1 x) := by
  convert S05.hs_L1 hx using 1; funext q; unfold g1; ring

lemma f_g2 {x : ℝ} (hx : |x| < 1) : HasSum (fun q => g2 q * x ^ (q + 1)) (D00.Li2 x) := by
  convert D00.hasSum_Li2 hx.le using 1; funext q; unfold g2; ring

lemma f_g3 {x : ℝ} (hx : |x| < 1) : HasSum (fun q => g3 q * x ^ (q + 1)) (Z03.Li3 x) := by
  convert Z03.hasSum_Li3 hx.le using 1; funext q; unfold g3 Z03.Li3term; ring

lemma f_gH1 {x : ℝ} (hx : |x| < 1) :
    HasSum (fun q => gH1 q * x ^ (q + 1)) (D00.Li2 x + 1 / 2 * S05.cL1 x ^ 2) := by
  convert S05.euler4 hx using 1
  funext q; rfl

lemma f_gH2 {x : ℝ} (hx : |x| < 1) : HasSum (fun q => gH2 q * x ^ (q + 1)) (Avl x) := by
  convert hs_A hx using 1; funext q; unfold gH2 Z03.Aterm; ring

lemma f_gHH {x : ℝ} (hx : |x| < 1) :
    HasSum (fun q => gHH q * x ^ (q + 1))
      (Z03.Li3 x + S05.cL1 x * D00.Li2 x + 1 / 3 * S05.cL1 x ^ 3) := by
  convert S05.euler5 hx using 1
  funext q; rfl

lemma f_gHt {x : ℝ} (hx : |x| < 1) :
    HasSum (fun q => gHt q * x ^ (q + 1))
      (3 * Z03.Li3 x - 2 * Avl x + S05.cL1 x * D00.Li2 x) := by
  convert S05.euler6 hx using 1
  · funext q; rfl
  · rw [A_eq hx]

lemma sh1 (g : ℕ → ℝ) {x F : ℝ} (h0 : x ≠ 0) (h : HasSum (fun q => g q * x ^ (q + 1)) F) :
    HasSum (fun q => g (q + 1) * x ^ (q + 1)) ((F - g 0 * x) / x) := S05.shift_div g h0 h

lemma sh2 (g : ℕ → ℝ) {x F : ℝ} (h0 : x ≠ 0) (h : HasSum (fun q => g q * x ^ (q + 1)) F) :
    HasSum (fun q => g (q + 1 + 1) * x ^ (q + 1)) (((F - g 0 * x) / x - g 1 * x) / x) :=
  S05.shift_div (fun q => g (q + 1)) h0 (sh1 g h0 h)

/-! ## harmonic normal forms -/

lemma hA1 (q : ℕ) : (harmonic (q + 1) : ℝ) = harmonic q + 1 / ((q : ℝ) + 1) := S05.hR_succ q
lemma hA2 (q : ℕ) : (harmonic (q + 1 + 1) : ℝ) = harmonic q + 1 / ((q : ℝ) + 1) + 1 / ((q : ℝ) + 2) := by
  rw [S05.hR_succ (q + 1), S05.hR_succ q]; push_cast; ring
lemma hA3 (q : ℕ) : (harmonic (q + 1 + 1 + 1) : ℝ) =
    harmonic q + 1 / ((q : ℝ) + 1) + 1 / ((q : ℝ) + 2) + 1 / ((q : ℝ) + 3) := by
  rw [S05.hR_succ (q + 1 + 1), hA2]; push_cast; ring
lemma tA1 (q : ℕ) : S05.H2 (q + 1) = S05.H2 q + 1 / ((q : ℝ) + 1) ^ 2 := S05.H2_succ q
lemma tA2 (q : ℕ) : S05.H2 (q + 1 + 1) = S05.H2 q + 1 / ((q : ℝ) + 1) ^ 2 + 1 / ((q : ℝ) + 2) ^ 2 := by
  rw [S05.H2_succ (q + 1), S05.H2_succ q]; push_cast; ring
lemma tA3 (q : ℕ) : S05.H2 (q + 1 + 1 + 1) =
    S05.H2 q + 1 / ((q : ℝ) + 1) ^ 2 + 1 / ((q : ℝ) + 2) ^ 2 + 1 / ((q : ℝ) + 3) ^ 2 := by
  rw [S05.H2_succ (q + 1 + 1), tA2]; push_cast; ring

/-! ## Beta values `𝓘_N(0,1)` and `𝓘_N(1,1)` -/

lemma Iint01 (N : ℕ) : S05.Iint N 0 1 = -(harmonic (N + 1) : ℝ) / ((N : ℝ) + 1) := by
  have hsub : S05.Iint N 0 1 = ∫ u in (0:ℝ)..1, (1 - u) ^ N * Real.log u := by
    unfold S05.Iint
    have := intervalIntegral.integral_comp_sub_left
      (fun u : ℝ => (1 - u) ^ N * Real.log u) (a := 0) (b := 1) 1
    simp only [sub_sub_cancel, sub_zero, sub_self] at this
    rw [← this]
    congr 1; funext s; ring
  rw [hsub, S05.intP]

lemma H2_range (M : ℕ) : S05.H2 M = ∑ i ∈ range M, 1 / ((i : ℝ) + 1) ^ 2 := by
  induction M with
  | zero => simp [S05.H2_zero]
  | succ M ih => rw [S05.H2_succ, ih, sum_range_succ]

lemma harm_range (M : ℕ) : (harmonic M : ℝ) = ∑ i ∈ range M, 1 / ((i : ℝ) + 1) := S05.hR_sum M

lemma tele_gen (f : ℕ → ℝ) (hf : Tendsto f atTop (𝓝 0)) (hmono : ∀ q, f (q + 1) ≤ f q) :
    HasSum (fun q => f q - f (q + 1)) (f 0) := by
  have hnn : ∀ q, 0 ≤ f q - f (q + 1) := fun q => sub_nonneg.2 (hmono q)
  rw [hasSum_iff_tendsto_nat_of_nonneg hnn]
  have e : ∀ n, ∑ i ∈ range n, (f i - f (i + 1)) = f 0 - f n := fun n => by
    rw [Finset.sum_range_sub']
  simp_rw [e]
  have := (tendsto_const_nhds : Tendsto (fun _ : ℕ => f 0) atTop (𝓝 (f 0))).sub hf
  simpa using this

lemma inv_shift_tendsto (M : ℕ) : Tendsto (fun q : ℕ => 1 / ((q : ℝ) + 1 + M)) atTop (𝓝 0) := by
  have h := (tendsto_one_div_add_atTop_nhds_zero_nat (𝕜 := ℝ)).comp (tendsto_add_atTop_nat M)
  refine h.congr fun q => ?_
  simp only [Function.comp]; push_cast; ring_nf

lemma tele_M (M : ℕ) :
    HasSum (fun q : ℕ => 1 / ((q : ℝ) + 1) - 1 / ((q : ℝ) + 1 + M)) (harmonic M : ℝ) := by
  induction M with
  | zero => simp
  | succ M ih =>
    have h2 := tele_gen (fun q : ℕ => 1 / ((q : ℝ) + 1 + M)) (inv_shift_tendsto M) (fun q => by
      push_cast
      gcongr; linarith)
    have := ih.add h2
    convert this using 1
    · funext q; push_cast; ring
    · rw [S05.hR_succ]; push_cast; ring

lemma zeta2_shift (M : ℕ) :
    HasSum (fun q : ℕ => 1 / ((q : ℝ) + 1 + M) ^ 2) (π ^ 2 / 6 - S05.H2 M) := by
  have h0 : HasSum (fun n : ℕ => 1 / ((n : ℝ) + 1) ^ 2) (π ^ 2 / 6) := by
    have := D00.hasSum_Li2 (x := 1) (by norm_num)
    rw [D00.proof_7] at this
    simpa using this
  have h := (hasSum_nat_add_iff' M).2 h0
  rw [← H2_range] at h
  convert h using 1
  funext q; push_cast; ring

lemma Iint11 (N : ℕ) : S05.Iint N 1 1 =
    (harmonic (N + 1) : ℝ) / ((N : ℝ) + 1) ^ 2 + S05.H2 (N + 1) / ((N : ℝ) + 1)
      - (π ^ 2 / 6) / ((N : ℝ) + 1) := by
  refine (S05.Iint_a1 N 1).unique ?_
  have hM : (0:ℝ) < (N : ℝ) + 1 := by positivity
  have h := ((tele_M (N + 1)).mul_left (1 / ((N : ℝ) + 1) ^ 2)).sub
    ((zeta2_shift (N + 1)).mul_left (1 / ((N : ℝ) + 1)))
  convert h using 1
  · funext q
    rw [S05.Iint10]
    push_cast
    field_simp
    ring
  · push_cast
    field_simp
    ring

lemma hv0 : (harmonic 0 : ℝ) = 0 := by simp
lemma hv1 : (harmonic 1 : ℝ) = 1 := by simp
lemma hv2 : (harmonic (1 + 1) : ℝ) = 3 / 2 := by rw [S05.hR_succ, hv1]; norm_num
lemma tv1 : S05.H2 1 = 1 := by simp [S05.H2]
lemma tv2 : S05.H2 (1 + 1) = 5 / 4 := by rw [S05.H2_succ, tv1]; norm_num

/-! ## 2. the twenty non-local `J` on `0 < |x| < 1` -/

noncomputable def CF0001 (x : ℝ) : ℝ :=
  (S05.cL1 x + ((-1) * x) + ((-1) * S05.cL1 x * x)) / (x)

noncomputable def CF0002 (x : ℝ) : ℝ :=
  (((-1) * (S05.cL1 x) ^ 2) + ((-2) * S05.cL1 x) + (2 * x) + (x * (S05.cL1 x) ^ 2) + (2 * S05.cL1 x * x)) / (x)

noncomputable def CF0003 (x : ℝ) : ℝ :=
  ((S05.cL1 x) ^ 3 + ((-6) * x) + (3 * (S05.cL1 x) ^ 2) + (6 * S05.cL1 x) + ((-1) * x * (S05.cL1 x) ^ 3) + ((-6) * S05.cL1 x * x) + ((-3) * x * (S05.cL1 x) ^ 2)) / (x)

noncomputable def CF0011 (x : ℝ) : ℝ :=
  (((-1) * (S05.cL1 x) ^ 2) + ((-2) * S05.cL1 x) + ((-2) * D00.Li2 x) + (4 * x) + (x * (S05.cL1 x) ^ 2) + (2 * S05.cL1 x * x) + (2 * D00.Li2 x * x)) / ((2 * x))

noncomputable def CF0012 (x : ℝ) : ℝ :=
  ((S05.cL1 x) ^ 3 + ((-18) * x) + ((-6) * Z03.Li3 x) + ((-6) * L3p x) + (6 * D00.Li2 x) + (6 * (S05.cL1 x) ^ 2) + (12 * S05.cL1 x) + ((-1) * x * (S05.cL1 x) ^ 3) + ((-12) * S05.cL1 x * x) + ((-6) * D00.Li2 x * x) + ((-6) * x * (S05.cL1 x) ^ 2) + (6 * Z03.Li3 x * x) + (6 * L3p x * x)) / ((3 * x))

noncomputable def CF0021 (x : ℝ) : ℝ :=
  ((S05.cL1 x) ^ 2 + ((-6) * x) + ((-2) * L3p x) + (2 * S05.cL1 x) + (2 * D00.Li2 x) + ((-1) * x * (S05.cL1 x) ^ 2) + ((-2) * S05.cL1 x * x) + ((-2) * D00.Li2 x * x) + (2 * L3p x * x)) / (x)

noncomputable def CF0101 (x : ℝ) : ℝ :=
  (((-1) * S05.cL1 x) + ((-1) * D00.Li2 x) + (2 * x) + (S05.cL1 x * x)) / (x)

noncomputable def CF0102 (x : ℝ) : ℝ :=
  ((S05.cL1 x) ^ 3 + ((-18) * x) + (3 * (S05.cL1 x) ^ 2) + (6 * D00.Li2 x) + (6 * Z03.Li3 x) + (6 * L3p x) + (12 * S05.cL1 x) + ((-12) * S05.cL1 x * x) + ((-3) * x * (S05.cL1 x) ^ 2) + (6 * S05.cL1 x * D00.Li2 x)) / ((3 * x))

noncomputable def CF0111 (x : ℝ) : ℝ :=
  (((-1) * (S05.cL1 x) ^ 3) + ((-36) * x) + ((-6) * L3p x) + (3 * (S05.cL1 x) ^ 2) + (6 * Z03.Li3 x) + (12 * S05.cL1 x) + (12 * D00.Li2 x) + ((-12) * S05.cL1 x * x) + ((-6) * S05.cL1 x * (π ^ 2 / 6)) + ((-6) * D00.Li2 x * x) + ((-3) * x * (S05.cL1 x) ^ 2) + (2 * x * (S05.cL1 x) ^ 3) + (6 * Z03.Li3 x * x) + (6 * x * (π ^ 2 / 6)) + (12 * L3p x * x) + (6 * S05.cL1 x * D00.Li2 x * x) + (6 * S05.cL1 x * x * (π ^ 2 / 6))) / ((6 * x))

noncomputable def CF0201 (x : ℝ) : ℝ :=
  (((-6) * x) + (2 * S05.cL1 x) + (2 * D00.Li2 x) + (2 * Z03.Li3 x) + ((-2) * S05.cL1 x * x)) / (x)

noncomputable def CF1001 (x : ℝ) : ℝ :=
  (((-1) * (x) ^ 2) + ((-2) * x) + (2 * S05.cL1 x) + ((-2) * S05.cL1 x * (x) ^ 2)) / ((4 * (x) ^ 2))

noncomputable def CF1002 (x : ℝ) : ℝ :=
  ((x) ^ 2 + ((-6) * S05.cL1 x) + ((-2) * (S05.cL1 x) ^ 2) + (6 * x) + (2 * S05.cL1 x * (x) ^ 2) + (2 * (S05.cL1 x) ^ 2 * (x) ^ 2) + (4 * S05.cL1 x * x)) / ((4 * (x) ^ 2))

noncomputable def CF1003 (x : ℝ) : ℝ :=
  (((-42) * x) + ((-3) * (x) ^ 2) + (4 * (S05.cL1 x) ^ 3) + (18 * (S05.cL1 x) ^ 2) + (42 * S05.cL1 x) + ((-36) * S05.cL1 x * x) + ((-12) * x * (S05.cL1 x) ^ 2) + ((-6) * S05.cL1 x * (x) ^ 2) + ((-6) * (S05.cL1 x) ^ 2 * (x) ^ 2) + ((-4) * (S05.cL1 x) ^ 3 * (x) ^ 2)) / ((8 * (x) ^ 2))

noncomputable def CF1011 (x : ℝ) : ℝ :=
  (((-1) * S05.cL1 x) + ((-1) * (S05.cL1 x) ^ 2) + ((-2) * D00.Li2 x) + (3 * x) + (4 * (x) ^ 2) + ((S05.cL1 x) ^ 2 * (x) ^ 2) + ((-2) * S05.cL1 x * x) + (2 * D00.Li2 x * (x) ^ 2) + (3 * S05.cL1 x * (x) ^ 2)) / ((4 * (x) ^ 2))

noncomputable def CF1012 (x : ℝ) : ℝ :=
  (((-84) * x) + ((-51) * (x) ^ 2) + ((-24) * Z03.Li3 x) + ((-24) * L3p x) + (4 * (S05.cL1 x) ^ 3) + (24 * (S05.cL1 x) ^ 2) + (36 * D00.Li2 x) + (48 * S05.cL1 x) + ((-48) * S05.cL1 x * (x) ^ 2) + ((-24) * D00.Li2 x * x) + ((-24) * (S05.cL1 x) ^ 2 * (x) ^ 2) + ((-12) * D00.Li2 x * (x) ^ 2) + ((-4) * (S05.cL1 x) ^ 3 * (x) ^ 2) + (24 * Z03.Li3 x * (x) ^ 2) + (24 * L3p x * (x) ^ 2)) / ((24 * (x) ^ 2))

noncomputable def CF1021 (x : ℝ) : ℝ :=
  (((-31) * (x) ^ 2) + ((-14) * x) + ((-8) * L3p x) + (2 * S05.cL1 x) + (2 * (S05.cL1 x) ^ 2) + (4 * D00.Li2 x) + ((-14) * S05.cL1 x * (x) ^ 2) + ((-12) * D00.Li2 x * (x) ^ 2) + ((-6) * (S05.cL1 x) ^ 2 * (x) ^ 2) + (4 * x * (S05.cL1 x) ^ 2) + (8 * D00.Li2 x * x) + (8 * L3p x * (x) ^ 2) + (12 * S05.cL1 x * x)) / ((8 * (x) ^ 2))

noncomputable def CF1101 (x : ℝ) : ℝ :=
  ((x) ^ 2 + ((-1) * S05.cL1 x) + ((-2) * D00.Li2 x) + (3 * x) + (S05.cL1 x * (x) ^ 2)) / ((4 * (x) ^ 2))

noncomputable def CF1102 (x : ℝ) : ℝ :=
  (((-84) * x) + ((-9) * (x) ^ 2) + (4 * (S05.cL1 x) ^ 3) + (6 * (S05.cL1 x) ^ 2) + (24 * Z03.Li3 x) + (24 * L3p x) + (36 * D00.Li2 x) + (48 * S05.cL1 x) + ((-36) * S05.cL1 x * x) + ((-12) * S05.cL1 x * (x) ^ 2) + ((-6) * (S05.cL1 x) ^ 2 * (x) ^ 2) + (24 * S05.cL1 x * D00.Li2 x)) / ((24 * (x) ^ 2))

noncomputable def CF1111 (x : ℝ) : ℝ :=
  (((-51) * (x) ^ 2) + ((-42) * x) + ((-12) * L3p x) + ((-2) * (S05.cL1 x) ^ 3) + (3 * (S05.cL1 x) ^ 2) + (6 * S05.cL1 x) + (12 * D00.Li2 x) + (12 * Z03.Li3 x) + ((-24) * S05.cL1 x * (x) ^ 2) + ((-12) * S05.cL1 x * (π ^ 2 / 6)) + ((-6) * D00.Li2 x * (x) ^ 2) + ((-3) * (S05.cL1 x) ^ 2 * (x) ^ 2) + (4 * (S05.cL1 x) ^ 3 * (x) ^ 2) + (6 * (π ^ 2 / 6) * (x) ^ 2) + (12 * D00.Li2 x * x) + (12 * Z03.Li3 x * (x) ^ 2) + (12 * x * (π ^ 2 / 6)) + (18 * S05.cL1 x * x) + (24 * L3p x * (x) ^ 2) + (12 * S05.cL1 x * D00.Li2 x * (x) ^ 2) + (12 * S05.cL1 x * (π ^ 2 / 6) * (x) ^ 2)) / ((24 * (x) ^ 2))

noncomputable def CF1201 (x : ℝ) : ℝ :=
  (((-14) * x) + ((-3) * (x) ^ 2) + (2 * S05.cL1 x) + (4 * D00.Li2 x) + (8 * Z03.Li3 x) + ((-2) * S05.cL1 x * (x) ^ 2)) / ((8 * (x) ^ 2))


lemma Ja0001 {x : ℝ} (hx : |x| < 1) (h0 : x ≠ 0) : S05.Jint 0 0 0 1 x = CF0001 x := by
  have h := (S05.proof_3 x hx 0 0 0 1 (by norm_num) (by norm_num)).2
  have hc := (((f_g1 hx).mul_left (-1 : ℝ)).add ((sh1 g1 h0 (f_g1 hx)).mul_left (1 : ℝ)))
  refine h.unique ?_
  convert hc using 1
  · funext q
    rw [show 0 + q + 1 = q + 1 by ring]
    simp only [S05.lam, S05.Iint00, g1, g2, g3, gH1, gH2, gHH, gHt, Nat.add_sub_cancel, hA1 q, hA2 q,
      hA3 q, tA1 q, tA2 q, tA3 q]
    push_cast
    field_simp
    ring
  · simp only [g1, g2, g3, gH1, gH2, gHH, gHt, CF0001, Avl, hv0, hv1, hv2, tv1, tv2, Nat.cast_zero, Nat.cast_one, zero_add]
    field_simp
    ring

lemma Ja0002 {x : ℝ} (hx : |x| < 1) (h0 : x ≠ 0) : S05.Jint 0 0 0 2 x = CF0002 x := by
  have h := (S05.proof_3 x hx 0 0 0 2 (by norm_num) (by norm_num)).2
  have hc := (((((((f_g1 hx).mul_left (2 : ℝ)).add ((sh1 g1 h0 (f_g1 hx)).mul_left (-2 : ℝ))).add ((f_g2 hx).mul_left (-2 : ℝ))).add ((sh1 g2 h0 (f_g2 hx)).mul_left (2 : ℝ))).add ((f_gH1 hx).mul_left (2 : ℝ))).add ((sh1 gH1 h0 (f_gH1 hx)).mul_left (-2 : ℝ)))
  refine h.unique ?_
  convert hc using 1
  · funext q
    rw [show 0 + q + 1 = q + 1 by ring]
    simp only [S05.lam, S05.Iint00, g1, g2, g3, gH1, gH2, gHH, gHt, Nat.add_sub_cancel, hA1 q, hA2 q,
      hA3 q, tA1 q, tA2 q, tA3 q]
    push_cast
    field_simp
    ring
  · simp only [g1, g2, g3, gH1, gH2, gHH, gHt, CF0002, Avl, hv0, hv1, hv2, tv1, tv2, Nat.cast_zero, Nat.cast_one, zero_add]
    field_simp
    ring

lemma Ja0003 {x : ℝ} (hx : |x| < 1) (h0 : x ≠ 0) : S05.Jint 0 0 0 3 x = CF0003 x := by
  have h := (S05.proof_3 x hx 0 0 0 3 (by norm_num) (by norm_num)).2
  have hc := (((((((((((((((f_g1 hx).mul_left (-6 : ℝ)).add ((sh1 g1 h0 (f_g1 hx)).mul_left (6 : ℝ))).add ((f_g2 hx).mul_left (6 : ℝ))).add ((sh1 g2 h0 (f_g2 hx)).mul_left (-6 : ℝ))).add ((f_g3 hx).mul_left (-6 : ℝ))).add ((sh1 g3 h0 (f_g3 hx)).mul_left (6 : ℝ))).add ((f_gH1 hx).mul_left (-6 : ℝ))).add ((sh1 gH1 h0 (f_gH1 hx)).mul_left (6 : ℝ))).add ((f_gH2 hx).mul_left (6 : ℝ))).add ((sh1 gH2 h0 (f_gH2 hx)).mul_left (-6 : ℝ))).add ((f_gHH hx).mul_left (-3 : ℝ))).add ((sh1 gHH h0 (f_gHH hx)).mul_left (3 : ℝ))).add ((f_gHt hx).mul_left (3 : ℝ))).add ((sh1 gHt h0 (f_gHt hx)).mul_left (-3 : ℝ)))
  refine h.unique ?_
  convert hc using 1
  · funext q
    rw [show 0 + q + 1 = q + 1 by ring]
    simp only [S05.lam, S05.Iint00, g1, g2, g3, gH1, gH2, gHH, gHt, Nat.add_sub_cancel, hA1 q, hA2 q,
      hA3 q, tA1 q, tA2 q, tA3 q]
    push_cast
    field_simp
    ring
  · simp only [g1, g2, g3, gH1, gH2, gHH, gHt, CF0003, Avl, hv0, hv1, hv2, tv1, tv2, Nat.cast_zero, Nat.cast_one, zero_add]
    field_simp
    ring

lemma Ja0011 {x : ℝ} (hx : |x| < 1) (h0 : x ≠ 0) : S05.Jint 0 0 1 1 x = CF0011 x := by
  have h := (S05.proof_3 x hx 0 0 1 1 (by norm_num) (by norm_num)).2
  have hc := (((((f_g1 hx).mul_left (1 : ℝ)).add ((sh1 g1 h0 (f_g1 hx)).mul_left (-1 : ℝ))).add ((f_gH1 hx).mul_left (1 : ℝ))).add ((sh1 gH1 h0 (f_gH1 hx)).mul_left (-1 : ℝ)))
  refine h.unique ?_
  convert hc using 1
  · funext q
    rw [show 0 + q + 1 = q + 1 by ring]
    simp only [S05.lam, Iint01, g1, g2, g3, gH1, gH2, gHH, gHt, Nat.add_sub_cancel, hA1 q, hA2 q,
      hA3 q, tA1 q, tA2 q, tA3 q]
    push_cast
    field_simp
    ring
  · simp only [g1, g2, g3, gH1, gH2, gHH, gHt, CF0011, Avl, hv0, hv1, hv2, tv1, tv2, Nat.cast_zero, Nat.cast_one, zero_add]
    field_simp
    ring

lemma Ja0012 {x : ℝ} (hx : |x| < 1) (h0 : x ≠ 0) : S05.Jint 0 0 1 2 x = CF0012 x := by
  have h := (S05.proof_3 x hx 0 0 1 2 (by norm_num) (by norm_num)).2
  have hc := (((((((((((f_g1 hx).mul_left (-4 : ℝ)).add ((sh1 g1 h0 (f_g1 hx)).mul_left (4 : ℝ))).add ((f_g2 hx).mul_left (2 : ℝ))).add ((sh1 g2 h0 (f_g2 hx)).mul_left (-2 : ℝ))).add ((f_gH1 hx).mul_left (-4 : ℝ))).add ((sh1 gH1 h0 (f_gH1 hx)).mul_left (4 : ℝ))).add ((f_gH2 hx).mul_left (2 : ℝ))).add ((sh1 gH2 h0 (f_gH2 hx)).mul_left (-2 : ℝ))).add ((f_gHH hx).mul_left (-2 : ℝ))).add ((sh1 gHH h0 (f_gHH hx)).mul_left (2 : ℝ)))
  refine h.unique ?_
  convert hc using 1
  · funext q
    rw [show 0 + q + 1 = q + 1 by ring]
    simp only [S05.lam, Iint01, g1, g2, g3, gH1, gH2, gHH, gHt, Nat.add_sub_cancel, hA1 q, hA2 q,
      hA3 q, tA1 q, tA2 q, tA3 q]
    push_cast
    field_simp
    ring
  · simp only [g1, g2, g3, gH1, gH2, gHH, gHt, CF0012, Avl, hv0, hv1, hv2, tv1, tv2, Nat.cast_zero, Nat.cast_one, zero_add]
    field_simp
    ring

lemma Ja0021 {x : ℝ} (hx : |x| < 1) (h0 : x ≠ 0) : S05.Jint 0 0 2 1 x = CF0021 x := by
  have h := (S05.proof_3 x hx 0 0 2 1 (by norm_num) (by norm_num)).2
  have hc := (((((((((f_g1 hx).mul_left (-2 : ℝ)).add ((sh1 g1 h0 (f_g1 hx)).mul_left (2 : ℝ))).add ((f_gH1 hx).mul_left (-2 : ℝ))).add ((sh1 gH1 h0 (f_gH1 hx)).mul_left (2 : ℝ))).add ((f_gHH hx).mul_left (-1 : ℝ))).add ((sh1 gHH h0 (f_gHH hx)).mul_left (1 : ℝ))).add ((f_gHt hx).mul_left (-1 : ℝ))).add ((sh1 gHt h0 (f_gHt hx)).mul_left (1 : ℝ)))
  refine h.unique ?_
  convert hc using 1
  · funext q
    rw [show 0 + q + 1 = q + 1 by ring]
    simp only [S05.lam, S05.Iint02, g1, g2, g3, gH1, gH2, gHH, gHt, Nat.add_sub_cancel, hA1 q, hA2 q,
      hA3 q, tA1 q, tA2 q, tA3 q]
    push_cast
    field_simp
    ring
  · simp only [g1, g2, g3, gH1, gH2, gHH, gHt, CF0021, Avl, hv0, hv1, hv2, tv1, tv2, Nat.cast_zero, Nat.cast_one, zero_add]
    field_simp
    ring

lemma Ja0101 {x : ℝ} (hx : |x| < 1) (h0 : x ≠ 0) : S05.Jint 0 1 0 1 x = CF0101 x := by
  have h := (S05.proof_3 x hx 0 1 0 1 (by norm_num) (by norm_num)).2
  have hc := ((((f_g1 hx).mul_left (1 : ℝ)).add ((sh1 g1 h0 (f_g1 hx)).mul_left (-1 : ℝ))).add ((sh1 g2 h0 (f_g2 hx)).mul_left (-1 : ℝ)))
  refine h.unique ?_
  convert hc using 1
  · funext q
    rw [show 0 + q + 1 = q + 1 by ring]
    simp only [S05.lam, S05.Iint10, g1, g2, g3, gH1, gH2, gHH, gHt, Nat.add_sub_cancel, hA1 q, hA2 q,
      hA3 q, tA1 q, tA2 q, tA3 q]
    push_cast
    field_simp
    ring
  · simp only [g1, g2, g3, gH1, gH2, gHH, gHt, CF0101, Avl, hv0, hv1, hv2, tv1, tv2, Nat.cast_zero, Nat.cast_one, zero_add]
    field_simp
    ring

lemma Ja0102 {x : ℝ} (hx : |x| < 1) (h0 : x ≠ 0) : S05.Jint 0 1 0 2 x = CF0102 x := by
  have h := (S05.proof_3 x hx 0 1 0 2 (by norm_num) (by norm_num)).2
  have hc := ((((((((f_g1 hx).mul_left (-4 : ℝ)).add ((sh1 g1 h0 (f_g1 hx)).mul_left (4 : ℝ))).add ((f_g2 hx).mul_left (2 : ℝ))).add ((sh1 g3 h0 (f_g3 hx)).mul_left (-2 : ℝ))).add ((f_gH1 hx).mul_left (-2 : ℝ))).add ((sh1 gH1 h0 (f_gH1 hx)).mul_left (2 : ℝ))).add ((sh1 gH2 h0 (f_gH2 hx)).mul_left (2 : ℝ)))
  refine h.unique ?_
  convert hc using 1
  · funext q
    rw [show 0 + q + 1 = q + 1 by ring]
    simp only [S05.lam, S05.Iint10, g1, g2, g3, gH1, gH2, gHH, gHt, Nat.add_sub_cancel, hA1 q, hA2 q,
      hA3 q, tA1 q, tA2 q, tA3 q]
    push_cast
    field_simp
    ring
  · simp only [g1, g2, g3, gH1, gH2, gHH, gHt, CF0102, Avl, hv0, hv1, hv2, tv1, tv2, Nat.cast_zero, Nat.cast_one, zero_add]
    field_simp
    ring

lemma Ja0111 {x : ℝ} (hx : |x| < 1) (h0 : x ≠ 0) : S05.Jint 0 1 1 1 x = CF0111 x := by
  have h := (S05.proof_3 x hx 0 1 1 1 (by norm_num) (by norm_num)).2
  have hc := (((((((((((f_g1 hx).mul_left (-2 : ℝ)).add ((sh1 g1 h0 (f_g1 hx)).mul_left (2 : ℝ))).add ((sh1 g2 h0 (f_g2 hx)).mul_left (1 : ℝ))).add ((f_gH1 hx).mul_left (-1 : ℝ))).add ((sh1 gH1 h0 (f_gH1 hx)).mul_left (1 : ℝ))).add ((sh1 gH2 h0 (f_gH2 hx)).mul_left (1 : ℝ))).add ((f_gHt hx).mul_left (-1 : ℝ))).add ((sh1 gHt h0 (f_gHt hx)).mul_left (1 : ℝ))).add ((f_g1 hx).mul_left ((1 : ℝ) * (π ^ 2 / 6)))).add ((sh1 g1 h0 (f_g1 hx)).mul_left ((-1 : ℝ) * (π ^ 2 / 6))))
  refine h.unique ?_
  convert hc using 1
  · funext q
    rw [show 0 + q + 1 = q + 1 by ring]
    simp only [S05.lam, Iint11, g1, g2, g3, gH1, gH2, gHH, gHt, Nat.add_sub_cancel, hA1 q, hA2 q,
      hA3 q, tA1 q, tA2 q, tA3 q]
    push_cast
    field_simp
    ring
  · simp only [g1, g2, g3, gH1, gH2, gHH, gHt, CF0111, Avl, hv0, hv1, hv2, tv1, tv2, Nat.cast_zero, Nat.cast_one, zero_add]
    field_simp
    ring

lemma Ja0201 {x : ℝ} (hx : |x| < 1) (h0 : x ≠ 0) : S05.Jint 0 2 0 1 x = CF0201 x := by
  have h := (S05.proof_3 x hx 0 2 0 1 (by norm_num) (by norm_num)).2
  have hc := (((((f_g1 hx).mul_left (-2 : ℝ)).add ((sh1 g1 h0 (f_g1 hx)).mul_left (2 : ℝ))).add ((sh1 g2 h0 (f_g2 hx)).mul_left (2 : ℝ))).add ((sh1 g3 h0 (f_g3 hx)).mul_left (2 : ℝ)))
  refine h.unique ?_
  convert hc using 1
  · funext q
    rw [show 0 + q + 1 = q + 1 by ring]
    simp only [S05.lam, S05.Iint20, g1, g2, g3, gH1, gH2, gHH, gHt, Nat.add_sub_cancel, hA1 q, hA2 q,
      hA3 q, tA1 q, tA2 q, tA3 q]
    push_cast
    field_simp
    ring
  · simp only [g1, g2, g3, gH1, gH2, gHH, gHt, CF0201, Avl, hv0, hv1, hv2, tv1, tv2, Nat.cast_zero, Nat.cast_one, zero_add]
    field_simp
    ring

lemma Ja1001 {x : ℝ} (hx : |x| < 1) (h0 : x ≠ 0) : S05.Jint 1 0 0 1 x = CF1001 x := by
  have h := (S05.proof_3 x hx 1 0 0 1 (by norm_num) (by norm_num)).2
  have hc := (((f_g1 hx).mul_left (-1 / 2 : ℝ)).add ((sh2 g1 h0 (f_g1 hx)).mul_left (1 / 2 : ℝ)))
  refine h.unique ?_
  convert hc using 1
  · funext q
    rw [show 1 + q + 1 = q + 1 + 1 by ring]
    simp only [S05.lam, S05.Iint00, g1, g2, g3, gH1, gH2, gHH, gHt, Nat.add_sub_cancel, hA1 q, hA2 q,
      hA3 q, tA1 q, tA2 q, tA3 q]
    push_cast
    field_simp
    ring
  · simp only [g1, g2, g3, gH1, gH2, gHH, gHt, CF1001, Avl, hv0, hv1, hv2, tv1, tv2, Nat.cast_zero, Nat.cast_one, zero_add]
    field_simp
    ring

lemma Ja1002 {x : ℝ} (hx : |x| < 1) (h0 : x ≠ 0) : S05.Jint 1 0 0 2 x = CF1002 x := by
  have h := (S05.proof_3 x hx 1 0 0 2 (by norm_num) (by norm_num)).2
  have hc := ((((((((f_g1 hx).mul_left (1 / 2 : ℝ)).add ((sh1 g1 h0 (f_g1 hx)).mul_left (1 : ℝ))).add ((sh2 g1 h0 (f_g1 hx)).mul_left (-3 / 2 : ℝ))).add ((f_g2 hx).mul_left (-1 : ℝ))).add ((sh2 g2 h0 (f_g2 hx)).mul_left (1 : ℝ))).add ((f_gH1 hx).mul_left (1 : ℝ))).add ((sh2 gH1 h0 (f_gH1 hx)).mul_left (-1 : ℝ)))
  refine h.unique ?_
  convert hc using 1
  · funext q
    rw [show 1 + q + 1 = q + 1 + 1 by ring]
    simp only [S05.lam, S05.Iint00, g1, g2, g3, gH1, gH2, gHH, gHt, Nat.add_sub_cancel, hA1 q, hA2 q,
      hA3 q, tA1 q, tA2 q, tA3 q]
    push_cast
    field_simp
    ring
  · simp only [g1, g2, g3, gH1, gH2, gHH, gHt, CF1002, Avl, hv0, hv1, hv2, tv1, tv2, Nat.cast_zero, Nat.cast_one, zero_add]
    field_simp
    ring

lemma Ja1003 {x : ℝ} (hx : |x| < 1) (h0 : x ≠ 0) : S05.Jint 1 0 0 3 x = CF1003 x := by
  have h := (S05.proof_3 x hx 1 0 0 3 (by norm_num) (by norm_num)).2
  have hc := ((((((((((((((((((f_g1 hx).mul_left (-3 / 4 : ℝ)).add ((sh1 g1 h0 (f_g1 hx)).mul_left (-9 / 2 : ℝ))).add ((sh2 g1 h0 (f_g1 hx)).mul_left (21 / 4 : ℝ))).add ((f_g2 hx).mul_left (3 / 2 : ℝ))).add ((sh1 g2 h0 (f_g2 hx)).mul_left (3 : ℝ))).add ((sh2 g2 h0 (f_g2 hx)).mul_left (-9 / 2 : ℝ))).add ((f_g3 hx).mul_left (-3 : ℝ))).add ((sh2 g3 h0 (f_g3 hx)).mul_left (3 : ℝ))).add ((f_gH1 hx).mul_left (-3 / 2 : ℝ))).add ((sh1 gH1 h0 (f_gH1 hx)).mul_left (-3 : ℝ))).add ((sh2 gH1 h0 (f_gH1 hx)).mul_left (9 / 2 : ℝ))).add ((f_gH2 hx).mul_left (3 : ℝ))).add ((sh2 gH2 h0 (f_gH2 hx)).mul_left (-3 : ℝ))).add ((f_gHH hx).mul_left (-3 / 2 : ℝ))).add ((sh2 gHH h0 (f_gHH hx)).mul_left (3 / 2 : ℝ))).add ((f_gHt hx).mul_left (3 / 2 : ℝ))).add ((sh2 gHt h0 (f_gHt hx)).mul_left (-3 / 2 : ℝ)))
  refine h.unique ?_
  convert hc using 1
  · funext q
    rw [show 1 + q + 1 = q + 1 + 1 by ring]
    simp only [S05.lam, S05.Iint00, g1, g2, g3, gH1, gH2, gHH, gHt, Nat.add_sub_cancel, hA1 q, hA2 q,
      hA3 q, tA1 q, tA2 q, tA3 q]
    push_cast
    field_simp
    ring
  · simp only [g1, g2, g3, gH1, gH2, gHH, gHt, CF1003, Avl, hv0, hv1, hv2, tv1, tv2, Nat.cast_zero, Nat.cast_one, zero_add]
    field_simp
    ring

lemma Ja1011 {x : ℝ} (hx : |x| < 1) (h0 : x ≠ 0) : S05.Jint 1 0 1 1 x = CF1011 x := by
  have h := (S05.proof_3 x hx 1 0 1 1 (by norm_num) (by norm_num)).2
  have hc := ((((((f_g1 hx).mul_left (3 / 4 : ℝ)).add ((sh1 g1 h0 (f_g1 hx)).mul_left (-1 / 2 : ℝ))).add ((sh2 g1 h0 (f_g1 hx)).mul_left (-1 / 4 : ℝ))).add ((f_gH1 hx).mul_left (1 / 2 : ℝ))).add ((sh2 gH1 h0 (f_gH1 hx)).mul_left (-1 / 2 : ℝ)))
  refine h.unique ?_
  convert hc using 1
  · funext q
    rw [show 1 + q + 1 = q + 1 + 1 by ring]
    simp only [S05.lam, Iint01, g1, g2, g3, gH1, gH2, gHH, gHt, Nat.add_sub_cancel, hA1 q, hA2 q,
      hA3 q, tA1 q, tA2 q, tA3 q]
    push_cast
    field_simp
    ring
  · simp only [g1, g2, g3, gH1, gH2, gHH, gHt, CF1011, Avl, hv0, hv1, hv2, tv1, tv2, Nat.cast_zero, Nat.cast_one, zero_add]
    field_simp
    ring

lemma Ja1012 {x : ℝ} (hx : |x| < 1) (h0 : x ≠ 0) : S05.Jint 1 0 1 2 x = CF1012 x := by
  have h := (S05.proof_3 x hx 1 0 1 2 (by norm_num) (by norm_num)).2
  have hc := ((((((((((((f_g1 hx).mul_left (-2 : ℝ)).add ((sh2 g1 h0 (f_g1 hx)).mul_left (2 : ℝ))).add ((f_g2 hx).mul_left (3 / 2 : ℝ))).add ((sh1 g2 h0 (f_g2 hx)).mul_left (-1 : ℝ))).add ((sh2 g2 h0 (f_g2 hx)).mul_left (-1 / 2 : ℝ))).add ((f_gH1 hx).mul_left (-2 : ℝ))).add ((sh2 gH1 h0 (f_gH1 hx)).mul_left (2 : ℝ))).add ((f_gH2 hx).mul_left (1 : ℝ))).add ((sh2 gH2 h0 (f_gH2 hx)).mul_left (-1 : ℝ))).add ((f_gHH hx).mul_left (-1 : ℝ))).add ((sh2 gHH h0 (f_gHH hx)).mul_left (1 : ℝ)))
  refine h.unique ?_
  convert hc using 1
  · funext q
    rw [show 1 + q + 1 = q + 1 + 1 by ring]
    simp only [S05.lam, Iint01, g1, g2, g3, gH1, gH2, gHH, gHt, Nat.add_sub_cancel, hA1 q, hA2 q,
      hA3 q, tA1 q, tA2 q, tA3 q]
    push_cast
    field_simp
    ring
  · simp only [g1, g2, g3, gH1, gH2, gHH, gHt, CF1012, Avl, hv0, hv1, hv2, tv1, tv2, Nat.cast_zero, Nat.cast_one, zero_add]
    field_simp
    ring

lemma Ja1021 {x : ℝ} (hx : |x| < 1) (h0 : x ≠ 0) : S05.Jint 1 0 2 1 x = CF1021 x := by
  have h := (S05.proof_3 x hx 1 0 2 1 (by norm_num) (by norm_num)).2
  have hc := (((((((((((f_g1 hx).mul_left (-7 / 4 : ℝ)).add ((sh1 g1 h0 (f_g1 hx)).mul_left (3 / 2 : ℝ))).add ((sh2 g1 h0 (f_g1 hx)).mul_left (1 / 4 : ℝ))).add ((f_gH1 hx).mul_left (-3 / 2 : ℝ))).add ((sh1 gH1 h0 (f_gH1 hx)).mul_left (1 : ℝ))).add ((sh2 gH1 h0 (f_gH1 hx)).mul_left (1 / 2 : ℝ))).add ((f_gHH hx).mul_left (-1 / 2 : ℝ))).add ((sh2 gHH h0 (f_gHH hx)).mul_left (1 / 2 : ℝ))).add ((f_gHt hx).mul_left (-1 / 2 : ℝ))).add ((sh2 gHt h0 (f_gHt hx)).mul_left (1 / 2 : ℝ)))
  refine h.unique ?_
  convert hc using 1
  · funext q
    rw [show 1 + q + 1 = q + 1 + 1 by ring]
    simp only [S05.lam, S05.Iint02, g1, g2, g3, gH1, gH2, gHH, gHt, Nat.add_sub_cancel, hA1 q, hA2 q,
      hA3 q, tA1 q, tA2 q, tA3 q]
    push_cast
    field_simp
    ring
  · simp only [g1, g2, g3, gH1, gH2, gHH, gHt, CF1021, Avl, hv0, hv1, hv2, tv1, tv2, Nat.cast_zero, Nat.cast_one, zero_add]
    field_simp
    ring

lemma Ja1101 {x : ℝ} (hx : |x| < 1) (h0 : x ≠ 0) : S05.Jint 1 1 0 1 x = CF1101 x := by
  have h := (S05.proof_3 x hx 1 1 0 1 (by norm_num) (by norm_num)).2
  have hc := ((((f_g1 hx).mul_left (1 / 4 : ℝ)).add ((sh2 g1 h0 (f_g1 hx)).mul_left (-1 / 4 : ℝ))).add ((sh2 g2 h0 (f_g2 hx)).mul_left (-1 / 2 : ℝ)))
  refine h.unique ?_
  convert hc using 1
  · funext q
    rw [show 1 + q + 1 = q + 1 + 1 by ring]
    simp only [S05.lam, S05.Iint10, g1, g2, g3, gH1, gH2, gHH, gHt, Nat.add_sub_cancel, hA1 q, hA2 q,
      hA3 q, tA1 q, tA2 q, tA3 q]
    push_cast
    field_simp
    ring
  · simp only [g1, g2, g3, gH1, gH2, gHH, gHt, CF1101, Avl, hv0, hv1, hv2, tv1, tv2, Nat.cast_zero, Nat.cast_one, zero_add]
    field_simp
    ring

lemma Ja1102 {x : ℝ} (hx : |x| < 1) (h0 : x ≠ 0) : S05.Jint 1 1 0 2 x = CF1102 x := by
  have h := (S05.proof_3 x hx 1 1 0 2 (by norm_num) (by norm_num)).2
  have hc := ((((((((((f_g1 hx).mul_left (-1 / 2 : ℝ)).add ((sh1 g1 h0 (f_g1 hx)).mul_left (-3 / 2 : ℝ))).add ((sh2 g1 h0 (f_g1 hx)).mul_left (2 : ℝ))).add ((f_g2 hx).mul_left (1 / 2 : ℝ))).add ((sh2 g2 h0 (f_g2 hx)).mul_left (1 : ℝ))).add ((sh2 g3 h0 (f_g3 hx)).mul_left (-1 : ℝ))).add ((f_gH1 hx).mul_left (-1 / 2 : ℝ))).add ((sh2 gH1 h0 (f_gH1 hx)).mul_left (1 / 2 : ℝ))).add ((sh2 gH2 h0 (f_gH2 hx)).mul_left (1 : ℝ)))
  refine h.unique ?_
  convert hc using 1
  · funext q
    rw [show 1 + q + 1 = q + 1 + 1 by ring]
    simp only [S05.lam, S05.Iint10, g1, g2, g3, gH1, gH2, gHH, gHt, Nat.add_sub_cancel, hA1 q, hA2 q,
      hA3 q, tA1 q, tA2 q, tA3 q]
    push_cast
    field_simp
    ring
  · simp only [g1, g2, g3, gH1, gH2, gHH, gHt, CF1102, Avl, hv0, hv1, hv2, tv1, tv2, Nat.cast_zero, Nat.cast_one, zero_add]
    field_simp
    ring

lemma Ja1111 {x : ℝ} (hx : |x| < 1) (h0 : x ≠ 0) : S05.Jint 1 1 1 1 x = CF1111 x := by
  have h := (S05.proof_3 x hx 1 1 1 1 (by norm_num) (by norm_num)).2
  have hc := (((((((((((((f_g1 hx).mul_left (-1 : ℝ)).add ((sh1 g1 h0 (f_g1 hx)).mul_left (3 / 4 : ℝ))).add ((sh2 g1 h0 (f_g1 hx)).mul_left (1 / 4 : ℝ))).add ((sh1 g2 h0 (f_g2 hx)).mul_left (1 / 2 : ℝ))).add ((sh2 g2 h0 (f_g2 hx)).mul_left (1 / 4 : ℝ))).add ((f_gH1 hx).mul_left (-1 / 4 : ℝ))).add ((sh2 gH1 h0 (f_gH1 hx)).mul_left (1 / 4 : ℝ))).add ((sh2 gH2 h0 (f_gH2 hx)).mul_left (1 / 2 : ℝ))).add ((f_gHt hx).mul_left (-1 / 2 : ℝ))).add ((sh2 gHt h0 (f_gHt hx)).mul_left (1 / 2 : ℝ))).add ((f_g1 hx).mul_left ((1 / 2 : ℝ) * (π ^ 2 / 6)))).add ((sh2 g1 h0 (f_g1 hx)).mul_left ((-1 / 2 : ℝ) * (π ^ 2 / 6))))
  refine h.unique ?_
  convert hc using 1
  · funext q
    rw [show 1 + q + 1 = q + 1 + 1 by ring]
    simp only [S05.lam, Iint11, g1, g2, g3, gH1, gH2, gHH, gHt, Nat.add_sub_cancel, hA1 q, hA2 q,
      hA3 q, tA1 q, tA2 q, tA3 q]
    push_cast
    field_simp
    ring
  · simp only [g1, g2, g3, gH1, gH2, gHH, gHt, CF1111, Avl, hv0, hv1, hv2, tv1, tv2, Nat.cast_zero, Nat.cast_one, zero_add]
    field_simp
    ring

lemma Ja1201 {x : ℝ} (hx : |x| < 1) (h0 : x ≠ 0) : S05.Jint 1 2 0 1 x = CF1201 x := by
  have h := (S05.proof_3 x hx 1 2 0 1 (by norm_num) (by norm_num)).2
  have hc := (((((f_g1 hx).mul_left (-1 / 4 : ℝ)).add ((sh2 g1 h0 (f_g1 hx)).mul_left (1 / 4 : ℝ))).add ((sh2 g2 h0 (f_g2 hx)).mul_left (1 / 2 : ℝ))).add ((sh2 g3 h0 (f_g3 hx)).mul_left (1 : ℝ)))
  refine h.unique ?_
  convert hc using 1
  · funext q
    rw [show 1 + q + 1 = q + 1 + 1 by ring]
    simp only [S05.lam, S05.Iint20, g1, g2, g3, gH1, gH2, gHH, gHt, Nat.add_sub_cancel, hA1 q, hA2 q,
      hA3 q, tA1 q, tA2 q, tA3 q]
    push_cast
    field_simp
    ring
  · simp only [g1, g2, g3, gH1, gH2, gHH, gHt, CF1201, Avl, hv0, hv1, hv2, tv1, tv2, Nat.cast_zero, Nat.cast_one, zero_add]
    field_simp
    ring


/-! ## 3. integrability for every `y < 1`, linear combinations, `t ↦ 1 - t` -/

/-- the monomial `t^n log^a t log^b(1-t) log^c(1-yt)`. -/
noncomputable def jmon (n a b c : ℕ) (y t : ℝ) : ℝ :=
  t ^ n * Real.log t ^ a * Real.log (1 - t) ^ b * Real.log (1 - y * t) ^ c

lemma one_sub_mul_pos {y s : ℝ} (hy : y < 1) (hs0 : 0 ≤ s) (hs1 : s ≤ 1) : 0 < 1 - y * s := by
  rcases le_or_gt y 0 with h | h
  · nlinarith
  · nlinarith

lemma ii_jmon {y : ℝ} (hy : y < 1) (n a b c : ℕ) : IntervalIntegrable (jmon n a b c y) volume 0 1 := by
  have hc : ContinuousOn (fun t : ℝ => Real.log (1 - y * t) ^ c) (uIcc 0 1) := by
    rw [Set.uIcc_of_le zero_le_one]
    refine ContinuousOn.pow (ContinuousOn.log (by fun_prop) fun t ht => ?_) c
    exact (one_sub_mul_pos hy ht.1 ht.2).ne'
  exact (S05.ii_g n a b).mul_continuousOn hc

/-- a finite linear combination of the `jmon`. -/
noncomputable def jcomb {m : ℕ} (kv : Fin m → ℝ) (nv av bv cv : Fin m → ℕ) (y t : ℝ) : ℝ :=
  ∑ i, kv i * jmon (nv i) (av i) (bv i) (cv i) y t

lemma int_jcomb {y : ℝ} (hy : y < 1) {m : ℕ} (kv : Fin m → ℝ) (nv av bv cv : Fin m → ℕ) :
    (∫ t in (0:ℝ)..1, jcomb kv nv av bv cv y t)
      = ∑ i, kv i * S05.Jint (nv i) (av i) (bv i) (cv i) y := by
  unfold jcomb
  rw [intervalIntegral.integral_finsetSum fun i _ => (ii_jmon hy _ _ _ _).const_mul _]
  refine Finset.sum_congr rfl fun i _ => ?_
  rw [intervalIntegral.integral_const_mul]; rfl

lemma Jint_zero (n a b : ℕ) (y : ℝ) : S05.Jint n a b 0 y = S05.Iint n a b := by
  simp [S05.Jint, S05.Iint]

lemma Jint_at_zero (n a b c : ℕ) (hc : c ≠ 0) : S05.Jint n a b c 0 = 0 := by
  simp [S05.Jint, zero_pow hc]

/-- the substitution `t ↦ 1 - t` with `1 - x(1-s) = (1-x)(1 - ys)`, `y = x/(x-1)`. -/
lemma landen_sub (n a b c : ℕ) {x : ℝ} (hx : x < 1) :
    S05.Jint n a b c x = ∫ s in (0:ℝ)..1, (1 - s) ^ n * Real.log (1 - s) ^ a * Real.log s ^ b *
      (Real.log (1 - x) + Real.log (1 - x / (x - 1) * s)) ^ c := by
  unfold S05.Jint
  have h := intervalIntegral.integral_comp_sub_left
    (fun t : ℝ => t ^ n * Real.log t ^ a * Real.log (1 - t) ^ b * Real.log (1 - x * t) ^ c)
    (a := 0) (b := 1) 1
  simp only [sub_zero, sub_self] at h
  rw [← h]
  refine intervalIntegral.integral_congr fun s hs => ?_
  rw [Set.uIcc_of_le zero_le_one] at hs
  have hy : x / (x - 1) < 1 := by
    rw [div_lt_one_of_neg (by linarith)]; linarith
  have h1 : 0 < 1 - x := by linarith
  have h2 := one_sub_mul_pos hy hs.1 hs.2
  have e : 1 - x * (1 - s) = (1 - x) * (1 - x / (x - 1) * s) := by
    have : x - 1 ≠ 0 := by linarith
    field_simp; ring
  simp only [sub_sub_cancel]
  rw [e, Real.log_mul h1.ne' h2.ne']

/-! ## the basis at `y = x/(x-1)` -/

lemma y_lt {x : ℝ} (hx : x < 1) : x / (x - 1) < 1 := by
  rw [div_lt_one_of_neg (by linarith)]; linarith

lemma y_abs {x : ℝ} (hx : x < 1 / 2) : |x / (x - 1)| < 1 := by
  rw [abs_lt]; constructor
  · rw [lt_div_iff_of_neg (by linarith)]; linarith
  · exact y_lt (by linarith)

lemma y_ne {x : ℝ} (hx : x < 1) (h0 : x ≠ 0) : x / (x - 1) ≠ 0 :=
  div_ne_zero h0 (by linarith)

lemma cL1_y {x : ℝ} (hx : x < 1) : S05.cL1 (x / (x - 1)) = -S05.cL1 x := by
  unfold S05.cL1
  have : 1 - x / (x - 1) = (1 - x)⁻¹ := by
    have : x - 1 ≠ 0 := by linarith
    have : 1 - x ≠ 0 := by linarith
    field_simp; ring
  rw [this, Real.log_inv]

lemma Li2_y {x : ℝ} (hx : x < 1) :
    D00.Li2 (x / (x - 1)) = -D00.Li2 x - (1 / 2) * S05.cL1 x ^ 2 := by
  have := D00.landen hx
  unfold S05.cL1
  linear_combination this

lemma Li3_y (x : ℝ) : Z03.Li3 (x / (x - 1)) = L3p x := rfl

lemma L3p_y {x : ℝ} (hx : x < 1) : L3p (x / (x - 1)) = Z03.Li3 x := by
  unfold L3p
  congr 1
  have : x - 1 ≠ 0 := by linarith
  have : x / (x - 1) - 1 ≠ 0 := by
    rw [div_sub_one this]; exact div_ne_zero (by norm_num) this
  field_simp
  ring

lemma log_one_sub {x : ℝ} : Real.log (1 - x) = -S05.cL1 x := by unfold S05.cL1; ring


/-! ## the local values -/

lemma Iv000 : S05.Iint 0 0 0 = 1 := by
  rw [S05.Iint00]; norm_num [hv1, hv2, tv1, tv2]

lemma Iv001 : S05.Iint 0 0 1 = (-1) := by
  rw [Iint01]; norm_num [hv1, hv2, tv1, tv2]

lemma Iv002 : S05.Iint 0 0 2 = 2 := by
  rw [S05.Iint02]; norm_num [hv1, hv2, tv1, tv2]

lemma Iv003 : S05.Iint 0 0 3 = (-6) := by
  rw [S05.I003] <;> ring

lemma Iv010 : S05.Iint 0 1 0 = (-1) := by
  rw [S05.Iint10]; norm_num [hv1, hv2, tv1, tv2]

lemma Iv011 : S05.Iint 0 1 1 = (2 + ((-1) * (π ^ 2 / 6))) := by
  rw [Iint11]; norm_num [hv1, hv2, tv1, tv2] <;> ring

lemma Iv012 : S05.Iint 0 1 2 = ((-6) + (2 * (π ^ 2 / 6)) + (2 * Z03.z3)) := by
  rw [S05.I012] <;> ring

lemma Iv020 : S05.Iint 0 2 0 = 2 := by
  rw [S05.Iint20]; norm_num [hv1, hv2, tv1, tv2]

lemma Iv021 : S05.Iint 0 2 1 = ((-6) + (2 * (π ^ 2 / 6)) + (2 * Z03.z3)) := by
  rw [S05.I021] <;> ring

lemma Iv030 : S05.Iint 0 3 0 = (-6) := by
  rw [S05.I030] <;> ring

lemma Iv100 : S05.Iint 1 0 0 = (1 / 2) := by
  rw [S05.Iint00]; norm_num [hv1, hv2, tv1, tv2]

lemma Iv101 : S05.Iint 1 0 1 = ((-3) / 4) := by
  rw [Iint01]; norm_num [hv1, hv2, tv1, tv2]

lemma Iv102 : S05.Iint 1 0 2 = (7 / 4) := by
  rw [S05.Iint02]; norm_num [hv1, hv2, tv1, tv2]

lemma Iv103 : S05.Iint 1 0 3 = ((-45) / 8) := by
  rw [S05.I103] <;> ring

lemma Iv110 : S05.Iint 1 1 0 = ((-1) / 4) := by
  rw [S05.Iint10]; norm_num [hv1, hv2, tv1, tv2]

lemma Iv111 : S05.Iint 1 1 1 = (1 + (((-1) * (π ^ 2 / 6)) / (2))) := by
  rw [Iint11]; norm_num [hv1, hv2, tv1, tv2] <;> ring

lemma Iv112 : S05.Iint 1 1 2 = (((-31) / 8) + Z03.z3 + ((3 * (π ^ 2 / 6)) / (2))) := by
  rw [S05.I112] <;> ring

lemma Iv120 : S05.Iint 1 2 0 = (1 / 4) := by
  rw [S05.Iint20]; norm_num [hv1, hv2, tv1, tv2]

lemma Iv121 : S05.Iint 1 2 1 = (((-17) / 8) + Z03.z3 + ((1 * (π ^ 2 / 6)) / (2))) := by
  rw [S05.I121] <;> ring

lemma Iv130 : S05.Iint 1 3 0 = ((-3) / 8) := by
  rw [S05.I130] <;> ring

/-! ## the twenty `J` on `x < 1/2` -/

lemma Jb0001 {x : ℝ} (hx : x < 1 / 2) (h0 : x ≠ 0) : S05.Jint 0 0 0 1 x = CF0001 x := by
  have hx1 : x < 1 := by linarith
  have hy := y_abs hx
  have hy0 := y_ne hx1 h0
  rw [landen_sub 0 0 0 1 hx1]
  have e : (fun s : ℝ => (1 - s) ^ 0 * Real.log (1 - s) ^ 0 * Real.log s ^ 0 *
      (Real.log (1 - x) + Real.log (1 - x / (x - 1) * s)) ^ 1) =
      jcomb ![1, Real.log (1 - x)] ![0, 0] ![0, 0] ![0, 0] ![1, 0] (x / (x - 1)) := by
    funext s
    simp [jcomb, jmon, Fin.sum_univ_succ]
    ring
  rw [e, int_jcomb (y_lt hx1)]
  simp only [Fin.sum_univ_succ, Fin.sum_univ_zero, Matrix.cons_val_zero, Matrix.cons_val_succ,
    Fin.succ_zero_eq_one, Fin.succ_one_eq_two, Matrix.cons_val_one, Matrix.head_cons, add_zero, Jint_zero]
  simp only [Ja0001 hy hy0, Iv000]
  simp only [CF0001, CF0001, cL1_y hx1, Li2_y hx1, Li3_y, L3p_y hx1, log_one_sub]
  have hx2 : x - 1 ≠ 0 := by linarith
  field_simp
  ring

lemma Jb0002 {x : ℝ} (hx : x < 1 / 2) (h0 : x ≠ 0) : S05.Jint 0 0 0 2 x = CF0002 x := by
  have hx1 : x < 1 := by linarith
  have hy := y_abs hx
  have hy0 := y_ne hx1 h0
  rw [landen_sub 0 0 0 2 hx1]
  have e : (fun s : ℝ => (1 - s) ^ 0 * Real.log (1 - s) ^ 0 * Real.log s ^ 0 *
      (Real.log (1 - x) + Real.log (1 - x / (x - 1) * s)) ^ 2) =
      jcomb ![1, (2 * Real.log (1 - x)), (Real.log (1 - x)) ^ 2] ![0, 0, 0] ![0, 0, 0] ![0, 0, 0] ![2, 1, 0] (x / (x - 1)) := by
    funext s
    simp [jcomb, jmon, Fin.sum_univ_succ]
    ring
  rw [e, int_jcomb (y_lt hx1)]
  simp only [Fin.sum_univ_succ, Fin.sum_univ_zero, Matrix.cons_val_zero, Matrix.cons_val_succ,
    Fin.succ_zero_eq_one, Fin.succ_one_eq_two, Matrix.cons_val_one, Matrix.head_cons, add_zero, Jint_zero]
  simp only [Ja0001 hy hy0, Ja0002 hy hy0, Iv000]
  simp only [CF0001, CF0002, CF0002, cL1_y hx1, Li2_y hx1, Li3_y, L3p_y hx1, log_one_sub]
  have hx2 : x - 1 ≠ 0 := by linarith
  field_simp
  ring

lemma Jb0003 {x : ℝ} (hx : x < 1 / 2) (h0 : x ≠ 0) : S05.Jint 0 0 0 3 x = CF0003 x := by
  have hx1 : x < 1 := by linarith
  have hy := y_abs hx
  have hy0 := y_ne hx1 h0
  rw [landen_sub 0 0 0 3 hx1]
  have e : (fun s : ℝ => (1 - s) ^ 0 * Real.log (1 - s) ^ 0 * Real.log s ^ 0 *
      (Real.log (1 - x) + Real.log (1 - x / (x - 1) * s)) ^ 3) =
      jcomb ![1, (3 * Real.log (1 - x)), (3 * (Real.log (1 - x)) ^ 2), (Real.log (1 - x)) ^ 3] ![0, 0, 0, 0] ![0, 0, 0, 0] ![0, 0, 0, 0] ![3, 2, 1, 0] (x / (x - 1)) := by
    funext s
    simp [jcomb, jmon, Fin.sum_univ_succ]
    ring
  rw [e, int_jcomb (y_lt hx1)]
  simp only [Fin.sum_univ_succ, Fin.sum_univ_zero, Matrix.cons_val_zero, Matrix.cons_val_succ,
    Fin.succ_zero_eq_one, Fin.succ_one_eq_two, Matrix.cons_val_one, Matrix.head_cons, add_zero, Jint_zero]
  simp only [Ja0001 hy hy0, Ja0002 hy hy0, Ja0003 hy hy0, Iv000]
  simp only [CF0001, CF0002, CF0003, CF0003, cL1_y hx1, Li2_y hx1, Li3_y, L3p_y hx1, log_one_sub]
  have hx2 : x - 1 ≠ 0 := by linarith
  field_simp
  ring

lemma Jb0011 {x : ℝ} (hx : x < 1 / 2) (h0 : x ≠ 0) : S05.Jint 0 0 1 1 x = CF0011 x := by
  have hx1 : x < 1 := by linarith
  have hy := y_abs hx
  have hy0 := y_ne hx1 h0
  rw [landen_sub 0 0 1 1 hx1]
  have e : (fun s : ℝ => (1 - s) ^ 0 * Real.log (1 - s) ^ 0 * Real.log s ^ 1 *
      (Real.log (1 - x) + Real.log (1 - x / (x - 1) * s)) ^ 1) =
      jcomb ![1, Real.log (1 - x)] ![0, 0] ![1, 1] ![0, 0] ![1, 0] (x / (x - 1)) := by
    funext s
    simp [jcomb, jmon, Fin.sum_univ_succ]
    ring
  rw [e, int_jcomb (y_lt hx1)]
  simp only [Fin.sum_univ_succ, Fin.sum_univ_zero, Matrix.cons_val_zero, Matrix.cons_val_succ,
    Fin.succ_zero_eq_one, Fin.succ_one_eq_two, Matrix.cons_val_one, Matrix.head_cons, add_zero, Jint_zero]
  simp only [Ja0101 hy hy0, Iv010]
  simp only [CF0101, CF0011, cL1_y hx1, Li2_y hx1, Li3_y, L3p_y hx1, log_one_sub]
  have hx2 : x - 1 ≠ 0 := by linarith
  field_simp
  ring

lemma Jb0012 {x : ℝ} (hx : x < 1 / 2) (h0 : x ≠ 0) : S05.Jint 0 0 1 2 x = CF0012 x := by
  have hx1 : x < 1 := by linarith
  have hy := y_abs hx
  have hy0 := y_ne hx1 h0
  rw [landen_sub 0 0 1 2 hx1]
  have e : (fun s : ℝ => (1 - s) ^ 0 * Real.log (1 - s) ^ 0 * Real.log s ^ 1 *
      (Real.log (1 - x) + Real.log (1 - x / (x - 1) * s)) ^ 2) =
      jcomb ![1, (2 * Real.log (1 - x)), (Real.log (1 - x)) ^ 2] ![0, 0, 0] ![1, 1, 1] ![0, 0, 0] ![2, 1, 0] (x / (x - 1)) := by
    funext s
    simp [jcomb, jmon, Fin.sum_univ_succ]
    ring
  rw [e, int_jcomb (y_lt hx1)]
  simp only [Fin.sum_univ_succ, Fin.sum_univ_zero, Matrix.cons_val_zero, Matrix.cons_val_succ,
    Fin.succ_zero_eq_one, Fin.succ_one_eq_two, Matrix.cons_val_one, Matrix.head_cons, add_zero, Jint_zero]
  simp only [Ja0101 hy hy0, Ja0102 hy hy0, Iv010]
  simp only [CF0101, CF0102, CF0012, cL1_y hx1, Li2_y hx1, Li3_y, L3p_y hx1, log_one_sub]
  have hx2 : x - 1 ≠ 0 := by linarith
  field_simp
  ring

lemma Jb0021 {x : ℝ} (hx : x < 1 / 2) (h0 : x ≠ 0) : S05.Jint 0 0 2 1 x = CF0021 x := by
  have hx1 : x < 1 := by linarith
  have hy := y_abs hx
  have hy0 := y_ne hx1 h0
  rw [landen_sub 0 0 2 1 hx1]
  have e : (fun s : ℝ => (1 - s) ^ 0 * Real.log (1 - s) ^ 0 * Real.log s ^ 2 *
      (Real.log (1 - x) + Real.log (1 - x / (x - 1) * s)) ^ 1) =
      jcomb ![1, Real.log (1 - x)] ![0, 0] ![2, 2] ![0, 0] ![1, 0] (x / (x - 1)) := by
    funext s
    simp [jcomb, jmon, Fin.sum_univ_succ]
    ring
  rw [e, int_jcomb (y_lt hx1)]
  simp only [Fin.sum_univ_succ, Fin.sum_univ_zero, Matrix.cons_val_zero, Matrix.cons_val_succ,
    Fin.succ_zero_eq_one, Fin.succ_one_eq_two, Matrix.cons_val_one, Matrix.head_cons, add_zero, Jint_zero]
  simp only [Ja0201 hy hy0, Iv020]
  simp only [CF0201, CF0021, cL1_y hx1, Li2_y hx1, Li3_y, L3p_y hx1, log_one_sub]
  have hx2 : x - 1 ≠ 0 := by linarith
  field_simp
  ring

lemma Jb0101 {x : ℝ} (hx : x < 1 / 2) (h0 : x ≠ 0) : S05.Jint 0 1 0 1 x = CF0101 x := by
  have hx1 : x < 1 := by linarith
  have hy := y_abs hx
  have hy0 := y_ne hx1 h0
  rw [landen_sub 0 1 0 1 hx1]
  have e : (fun s : ℝ => (1 - s) ^ 0 * Real.log (1 - s) ^ 1 * Real.log s ^ 0 *
      (Real.log (1 - x) + Real.log (1 - x / (x - 1) * s)) ^ 1) =
      jcomb ![1, Real.log (1 - x)] ![0, 0] ![0, 0] ![1, 1] ![1, 0] (x / (x - 1)) := by
    funext s
    simp [jcomb, jmon, Fin.sum_univ_succ]
    ring
  rw [e, int_jcomb (y_lt hx1)]
  simp only [Fin.sum_univ_succ, Fin.sum_univ_zero, Matrix.cons_val_zero, Matrix.cons_val_succ,
    Fin.succ_zero_eq_one, Fin.succ_one_eq_two, Matrix.cons_val_one, Matrix.head_cons, add_zero, Jint_zero]
  simp only [Ja0011 hy hy0, Iv001]
  simp only [CF0011, CF0101, cL1_y hx1, Li2_y hx1, Li3_y, L3p_y hx1, log_one_sub]
  have hx2 : x - 1 ≠ 0 := by linarith
  field_simp
  ring

lemma Jb0102 {x : ℝ} (hx : x < 1 / 2) (h0 : x ≠ 0) : S05.Jint 0 1 0 2 x = CF0102 x := by
  have hx1 : x < 1 := by linarith
  have hy := y_abs hx
  have hy0 := y_ne hx1 h0
  rw [landen_sub 0 1 0 2 hx1]
  have e : (fun s : ℝ => (1 - s) ^ 0 * Real.log (1 - s) ^ 1 * Real.log s ^ 0 *
      (Real.log (1 - x) + Real.log (1 - x / (x - 1) * s)) ^ 2) =
      jcomb ![1, (2 * Real.log (1 - x)), (Real.log (1 - x)) ^ 2] ![0, 0, 0] ![0, 0, 0] ![1, 1, 1] ![2, 1, 0] (x / (x - 1)) := by
    funext s
    simp [jcomb, jmon, Fin.sum_univ_succ]
    ring
  rw [e, int_jcomb (y_lt hx1)]
  simp only [Fin.sum_univ_succ, Fin.sum_univ_zero, Matrix.cons_val_zero, Matrix.cons_val_succ,
    Fin.succ_zero_eq_one, Fin.succ_one_eq_two, Matrix.cons_val_one, Matrix.head_cons, add_zero, Jint_zero]
  simp only [Ja0011 hy hy0, Ja0012 hy hy0, Iv001]
  simp only [CF0011, CF0012, CF0102, cL1_y hx1, Li2_y hx1, Li3_y, L3p_y hx1, log_one_sub]
  have hx2 : x - 1 ≠ 0 := by linarith
  field_simp
  ring

lemma Jb0111 {x : ℝ} (hx : x < 1 / 2) (h0 : x ≠ 0) : S05.Jint 0 1 1 1 x = CF0111 x := by
  have hx1 : x < 1 := by linarith
  have hy := y_abs hx
  have hy0 := y_ne hx1 h0
  rw [landen_sub 0 1 1 1 hx1]
  have e : (fun s : ℝ => (1 - s) ^ 0 * Real.log (1 - s) ^ 1 * Real.log s ^ 1 *
      (Real.log (1 - x) + Real.log (1 - x / (x - 1) * s)) ^ 1) =
      jcomb ![1, Real.log (1 - x)] ![0, 0] ![1, 1] ![1, 1] ![1, 0] (x / (x - 1)) := by
    funext s
    simp [jcomb, jmon, Fin.sum_univ_succ]
    ring
  rw [e, int_jcomb (y_lt hx1)]
  simp only [Fin.sum_univ_succ, Fin.sum_univ_zero, Matrix.cons_val_zero, Matrix.cons_val_succ,
    Fin.succ_zero_eq_one, Fin.succ_one_eq_two, Matrix.cons_val_one, Matrix.head_cons, add_zero, Jint_zero]
  simp only [Ja0111 hy hy0, Iv011]
  simp only [CF0111, CF0111, cL1_y hx1, Li2_y hx1, Li3_y, L3p_y hx1, log_one_sub]
  have hx2 : x - 1 ≠ 0 := by linarith
  field_simp
  ring

lemma Jb0201 {x : ℝ} (hx : x < 1 / 2) (h0 : x ≠ 0) : S05.Jint 0 2 0 1 x = CF0201 x := by
  have hx1 : x < 1 := by linarith
  have hy := y_abs hx
  have hy0 := y_ne hx1 h0
  rw [landen_sub 0 2 0 1 hx1]
  have e : (fun s : ℝ => (1 - s) ^ 0 * Real.log (1 - s) ^ 2 * Real.log s ^ 0 *
      (Real.log (1 - x) + Real.log (1 - x / (x - 1) * s)) ^ 1) =
      jcomb ![1, Real.log (1 - x)] ![0, 0] ![0, 0] ![2, 2] ![1, 0] (x / (x - 1)) := by
    funext s
    simp [jcomb, jmon, Fin.sum_univ_succ]
    ring
  rw [e, int_jcomb (y_lt hx1)]
  simp only [Fin.sum_univ_succ, Fin.sum_univ_zero, Matrix.cons_val_zero, Matrix.cons_val_succ,
    Fin.succ_zero_eq_one, Fin.succ_one_eq_two, Matrix.cons_val_one, Matrix.head_cons, add_zero, Jint_zero]
  simp only [Ja0021 hy hy0, Iv002]
  simp only [CF0021, CF0201, cL1_y hx1, Li2_y hx1, Li3_y, L3p_y hx1, log_one_sub]
  have hx2 : x - 1 ≠ 0 := by linarith
  field_simp
  ring

lemma Jb1001 {x : ℝ} (hx : x < 1 / 2) (h0 : x ≠ 0) : S05.Jint 1 0 0 1 x = CF1001 x := by
  have hx1 : x < 1 := by linarith
  have hy := y_abs hx
  have hy0 := y_ne hx1 h0
  rw [landen_sub 1 0 0 1 hx1]
  have e : (fun s : ℝ => (1 - s) ^ 1 * Real.log (1 - s) ^ 0 * Real.log s ^ 0 *
      (Real.log (1 - x) + Real.log (1 - x / (x - 1) * s)) ^ 1) =
      jcomb ![(-1), ((-1) * Real.log (1 - x)), 1, Real.log (1 - x)] ![1, 1, 0, 0] ![0, 0, 0, 0] ![0, 0, 0, 0] ![1, 0, 1, 0] (x / (x - 1)) := by
    funext s
    simp [jcomb, jmon, Fin.sum_univ_succ]
    ring
  rw [e, int_jcomb (y_lt hx1)]
  simp only [Fin.sum_univ_succ, Fin.sum_univ_zero, Matrix.cons_val_zero, Matrix.cons_val_succ,
    Fin.succ_zero_eq_one, Fin.succ_one_eq_two, Matrix.cons_val_one, Matrix.head_cons, add_zero, Jint_zero]
  simp only [Ja0001 hy hy0, Ja1001 hy hy0, Iv000, Iv100]
  simp only [CF0001, CF1001, CF1001, cL1_y hx1, Li2_y hx1, Li3_y, L3p_y hx1, log_one_sub]
  have hx2 : x - 1 ≠ 0 := by linarith
  field_simp
  ring

lemma Jb1002 {x : ℝ} (hx : x < 1 / 2) (h0 : x ≠ 0) : S05.Jint 1 0 0 2 x = CF1002 x := by
  have hx1 : x < 1 := by linarith
  have hy := y_abs hx
  have hy0 := y_ne hx1 h0
  rw [landen_sub 1 0 0 2 hx1]
  have e : (fun s : ℝ => (1 - s) ^ 1 * Real.log (1 - s) ^ 0 * Real.log s ^ 0 *
      (Real.log (1 - x) + Real.log (1 - x / (x - 1) * s)) ^ 2) =
      jcomb ![(-1), ((-2) * Real.log (1 - x)), ((-1) * (Real.log (1 - x)) ^ 2), 1, (2 * Real.log (1 - x)), (Real.log (1 - x)) ^ 2] ![1, 1, 1, 0, 0, 0] ![0, 0, 0, 0, 0, 0] ![0, 0, 0, 0, 0, 0] ![2, 1, 0, 2, 1, 0] (x / (x - 1)) := by
    funext s
    simp [jcomb, jmon, Fin.sum_univ_succ]
    ring
  rw [e, int_jcomb (y_lt hx1)]
  simp only [Fin.sum_univ_succ, Fin.sum_univ_zero, Matrix.cons_val_zero, Matrix.cons_val_succ,
    Fin.succ_zero_eq_one, Fin.succ_one_eq_two, Matrix.cons_val_one, Matrix.head_cons, add_zero, Jint_zero]
  simp only [Ja0001 hy hy0, Ja0002 hy hy0, Ja1001 hy hy0, Ja1002 hy hy0, Iv000, Iv100]
  simp only [CF0001, CF0002, CF1001, CF1002, CF1002, cL1_y hx1, Li2_y hx1, Li3_y, L3p_y hx1, log_one_sub]
  have hx2 : x - 1 ≠ 0 := by linarith
  field_simp
  ring

lemma Jb1003 {x : ℝ} (hx : x < 1 / 2) (h0 : x ≠ 0) : S05.Jint 1 0 0 3 x = CF1003 x := by
  have hx1 : x < 1 := by linarith
  have hy := y_abs hx
  have hy0 := y_ne hx1 h0
  rw [landen_sub 1 0 0 3 hx1]
  have e : (fun s : ℝ => (1 - s) ^ 1 * Real.log (1 - s) ^ 0 * Real.log s ^ 0 *
      (Real.log (1 - x) + Real.log (1 - x / (x - 1) * s)) ^ 3) =
      jcomb ![(-1), ((-3) * Real.log (1 - x)), ((-3) * (Real.log (1 - x)) ^ 2), ((-1) * (Real.log (1 - x)) ^ 3), 1, (3 * Real.log (1 - x)), (3 * (Real.log (1 - x)) ^ 2), (Real.log (1 - x)) ^ 3] ![1, 1, 1, 1, 0, 0, 0, 0] ![0, 0, 0, 0, 0, 0, 0, 0] ![0, 0, 0, 0, 0, 0, 0, 0] ![3, 2, 1, 0, 3, 2, 1, 0] (x / (x - 1)) := by
    funext s
    simp [jcomb, jmon, Fin.sum_univ_succ]
    ring
  rw [e, int_jcomb (y_lt hx1)]
  simp only [Fin.sum_univ_succ, Fin.sum_univ_zero, Matrix.cons_val_zero, Matrix.cons_val_succ,
    Fin.succ_zero_eq_one, Fin.succ_one_eq_two, Matrix.cons_val_one, Matrix.head_cons, add_zero, Jint_zero]
  simp only [Ja0001 hy hy0, Ja0002 hy hy0, Ja0003 hy hy0, Ja1001 hy hy0, Ja1002 hy hy0, Ja1003 hy hy0, Iv000, Iv100]
  simp only [CF0001, CF0002, CF0003, CF1001, CF1002, CF1003, CF1003, cL1_y hx1, Li2_y hx1, Li3_y, L3p_y hx1, log_one_sub]
  have hx2 : x - 1 ≠ 0 := by linarith
  field_simp
  ring

lemma Jb1011 {x : ℝ} (hx : x < 1 / 2) (h0 : x ≠ 0) : S05.Jint 1 0 1 1 x = CF1011 x := by
  have hx1 : x < 1 := by linarith
  have hy := y_abs hx
  have hy0 := y_ne hx1 h0
  rw [landen_sub 1 0 1 1 hx1]
  have e : (fun s : ℝ => (1 - s) ^ 1 * Real.log (1 - s) ^ 0 * Real.log s ^ 1 *
      (Real.log (1 - x) + Real.log (1 - x / (x - 1) * s)) ^ 1) =
      jcomb ![(-1), ((-1) * Real.log (1 - x)), 1, Real.log (1 - x)] ![1, 1, 0, 0] ![1, 1, 1, 1] ![0, 0, 0, 0] ![1, 0, 1, 0] (x / (x - 1)) := by
    funext s
    simp [jcomb, jmon, Fin.sum_univ_succ]
    ring
  rw [e, int_jcomb (y_lt hx1)]
  simp only [Fin.sum_univ_succ, Fin.sum_univ_zero, Matrix.cons_val_zero, Matrix.cons_val_succ,
    Fin.succ_zero_eq_one, Fin.succ_one_eq_two, Matrix.cons_val_one, Matrix.head_cons, add_zero, Jint_zero]
  simp only [Ja0101 hy hy0, Ja1101 hy hy0, Iv010, Iv110]
  simp only [CF0101, CF1101, CF1011, cL1_y hx1, Li2_y hx1, Li3_y, L3p_y hx1, log_one_sub]
  have hx2 : x - 1 ≠ 0 := by linarith
  field_simp
  ring

lemma Jb1012 {x : ℝ} (hx : x < 1 / 2) (h0 : x ≠ 0) : S05.Jint 1 0 1 2 x = CF1012 x := by
  have hx1 : x < 1 := by linarith
  have hy := y_abs hx
  have hy0 := y_ne hx1 h0
  rw [landen_sub 1 0 1 2 hx1]
  have e : (fun s : ℝ => (1 - s) ^ 1 * Real.log (1 - s) ^ 0 * Real.log s ^ 1 *
      (Real.log (1 - x) + Real.log (1 - x / (x - 1) * s)) ^ 2) =
      jcomb ![(-1), ((-2) * Real.log (1 - x)), ((-1) * (Real.log (1 - x)) ^ 2), 1, (2 * Real.log (1 - x)), (Real.log (1 - x)) ^ 2] ![1, 1, 1, 0, 0, 0] ![1, 1, 1, 1, 1, 1] ![0, 0, 0, 0, 0, 0] ![2, 1, 0, 2, 1, 0] (x / (x - 1)) := by
    funext s
    simp [jcomb, jmon, Fin.sum_univ_succ]
    ring
  rw [e, int_jcomb (y_lt hx1)]
  simp only [Fin.sum_univ_succ, Fin.sum_univ_zero, Matrix.cons_val_zero, Matrix.cons_val_succ,
    Fin.succ_zero_eq_one, Fin.succ_one_eq_two, Matrix.cons_val_one, Matrix.head_cons, add_zero, Jint_zero]
  simp only [Ja0101 hy hy0, Ja0102 hy hy0, Ja1101 hy hy0, Ja1102 hy hy0, Iv010, Iv110]
  simp only [CF0101, CF0102, CF1101, CF1102, CF1012, cL1_y hx1, Li2_y hx1, Li3_y, L3p_y hx1, log_one_sub]
  have hx2 : x - 1 ≠ 0 := by linarith
  field_simp
  ring

lemma Jb1021 {x : ℝ} (hx : x < 1 / 2) (h0 : x ≠ 0) : S05.Jint 1 0 2 1 x = CF1021 x := by
  have hx1 : x < 1 := by linarith
  have hy := y_abs hx
  have hy0 := y_ne hx1 h0
  rw [landen_sub 1 0 2 1 hx1]
  have e : (fun s : ℝ => (1 - s) ^ 1 * Real.log (1 - s) ^ 0 * Real.log s ^ 2 *
      (Real.log (1 - x) + Real.log (1 - x / (x - 1) * s)) ^ 1) =
      jcomb ![(-1), ((-1) * Real.log (1 - x)), 1, Real.log (1 - x)] ![1, 1, 0, 0] ![2, 2, 2, 2] ![0, 0, 0, 0] ![1, 0, 1, 0] (x / (x - 1)) := by
    funext s
    simp [jcomb, jmon, Fin.sum_univ_succ]
    ring
  rw [e, int_jcomb (y_lt hx1)]
  simp only [Fin.sum_univ_succ, Fin.sum_univ_zero, Matrix.cons_val_zero, Matrix.cons_val_succ,
    Fin.succ_zero_eq_one, Fin.succ_one_eq_two, Matrix.cons_val_one, Matrix.head_cons, add_zero, Jint_zero]
  simp only [Ja0201 hy hy0, Ja1201 hy hy0, Iv020, Iv120]
  simp only [CF0201, CF1201, CF1021, cL1_y hx1, Li2_y hx1, Li3_y, L3p_y hx1, log_one_sub]
  have hx2 : x - 1 ≠ 0 := by linarith
  field_simp
  ring

lemma Jb1101 {x : ℝ} (hx : x < 1 / 2) (h0 : x ≠ 0) : S05.Jint 1 1 0 1 x = CF1101 x := by
  have hx1 : x < 1 := by linarith
  have hy := y_abs hx
  have hy0 := y_ne hx1 h0
  rw [landen_sub 1 1 0 1 hx1]
  have e : (fun s : ℝ => (1 - s) ^ 1 * Real.log (1 - s) ^ 1 * Real.log s ^ 0 *
      (Real.log (1 - x) + Real.log (1 - x / (x - 1) * s)) ^ 1) =
      jcomb ![(-1), ((-1) * Real.log (1 - x)), 1, Real.log (1 - x)] ![1, 1, 0, 0] ![0, 0, 0, 0] ![1, 1, 1, 1] ![1, 0, 1, 0] (x / (x - 1)) := by
    funext s
    simp [jcomb, jmon, Fin.sum_univ_succ]
    ring
  rw [e, int_jcomb (y_lt hx1)]
  simp only [Fin.sum_univ_succ, Fin.sum_univ_zero, Matrix.cons_val_zero, Matrix.cons_val_succ,
    Fin.succ_zero_eq_one, Fin.succ_one_eq_two, Matrix.cons_val_one, Matrix.head_cons, add_zero, Jint_zero]
  simp only [Ja0011 hy hy0, Ja1011 hy hy0, Iv001, Iv101]
  simp only [CF0011, CF1011, CF1101, cL1_y hx1, Li2_y hx1, Li3_y, L3p_y hx1, log_one_sub]
  have hx2 : x - 1 ≠ 0 := by linarith
  field_simp
  ring

lemma Jb1102 {x : ℝ} (hx : x < 1 / 2) (h0 : x ≠ 0) : S05.Jint 1 1 0 2 x = CF1102 x := by
  have hx1 : x < 1 := by linarith
  have hy := y_abs hx
  have hy0 := y_ne hx1 h0
  rw [landen_sub 1 1 0 2 hx1]
  have e : (fun s : ℝ => (1 - s) ^ 1 * Real.log (1 - s) ^ 1 * Real.log s ^ 0 *
      (Real.log (1 - x) + Real.log (1 - x / (x - 1) * s)) ^ 2) =
      jcomb ![(-1), ((-2) * Real.log (1 - x)), ((-1) * (Real.log (1 - x)) ^ 2), 1, (2 * Real.log (1 - x)), (Real.log (1 - x)) ^ 2] ![1, 1, 1, 0, 0, 0] ![0, 0, 0, 0, 0, 0] ![1, 1, 1, 1, 1, 1] ![2, 1, 0, 2, 1, 0] (x / (x - 1)) := by
    funext s
    simp [jcomb, jmon, Fin.sum_univ_succ]
    ring
  rw [e, int_jcomb (y_lt hx1)]
  simp only [Fin.sum_univ_succ, Fin.sum_univ_zero, Matrix.cons_val_zero, Matrix.cons_val_succ,
    Fin.succ_zero_eq_one, Fin.succ_one_eq_two, Matrix.cons_val_one, Matrix.head_cons, add_zero, Jint_zero]
  simp only [Ja0011 hy hy0, Ja0012 hy hy0, Ja1011 hy hy0, Ja1012 hy hy0, Iv001, Iv101]
  simp only [CF0011, CF0012, CF1011, CF1012, CF1102, cL1_y hx1, Li2_y hx1, Li3_y, L3p_y hx1, log_one_sub]
  have hx2 : x - 1 ≠ 0 := by linarith
  field_simp
  ring

lemma Jb1111 {x : ℝ} (hx : x < 1 / 2) (h0 : x ≠ 0) : S05.Jint 1 1 1 1 x = CF1111 x := by
  have hx1 : x < 1 := by linarith
  have hy := y_abs hx
  have hy0 := y_ne hx1 h0
  rw [landen_sub 1 1 1 1 hx1]
  have e : (fun s : ℝ => (1 - s) ^ 1 * Real.log (1 - s) ^ 1 * Real.log s ^ 1 *
      (Real.log (1 - x) + Real.log (1 - x / (x - 1) * s)) ^ 1) =
      jcomb ![(-1), ((-1) * Real.log (1 - x)), 1, Real.log (1 - x)] ![1, 1, 0, 0] ![1, 1, 1, 1] ![1, 1, 1, 1] ![1, 0, 1, 0] (x / (x - 1)) := by
    funext s
    simp [jcomb, jmon, Fin.sum_univ_succ]
    ring
  rw [e, int_jcomb (y_lt hx1)]
  simp only [Fin.sum_univ_succ, Fin.sum_univ_zero, Matrix.cons_val_zero, Matrix.cons_val_succ,
    Fin.succ_zero_eq_one, Fin.succ_one_eq_two, Matrix.cons_val_one, Matrix.head_cons, add_zero, Jint_zero]
  simp only [Ja0111 hy hy0, Ja1111 hy hy0, Iv011, Iv111]
  simp only [CF0111, CF1111, CF1111, cL1_y hx1, Li2_y hx1, Li3_y, L3p_y hx1, log_one_sub]
  have hx2 : x - 1 ≠ 0 := by linarith
  field_simp
  ring

lemma Jb1201 {x : ℝ} (hx : x < 1 / 2) (h0 : x ≠ 0) : S05.Jint 1 2 0 1 x = CF1201 x := by
  have hx1 : x < 1 := by linarith
  have hy := y_abs hx
  have hy0 := y_ne hx1 h0
  rw [landen_sub 1 2 0 1 hx1]
  have e : (fun s : ℝ => (1 - s) ^ 1 * Real.log (1 - s) ^ 2 * Real.log s ^ 0 *
      (Real.log (1 - x) + Real.log (1 - x / (x - 1) * s)) ^ 1) =
      jcomb ![(-1), ((-1) * Real.log (1 - x)), 1, Real.log (1 - x)] ![1, 1, 0, 0] ![0, 0, 0, 0] ![2, 2, 2, 2] ![1, 0, 1, 0] (x / (x - 1)) := by
    funext s
    simp [jcomb, jmon, Fin.sum_univ_succ]
    ring
  rw [e, int_jcomb (y_lt hx1)]
  simp only [Fin.sum_univ_succ, Fin.sum_univ_zero, Matrix.cons_val_zero, Matrix.cons_val_succ,
    Fin.succ_zero_eq_one, Fin.succ_one_eq_two, Matrix.cons_val_one, Matrix.head_cons, add_zero, Jint_zero]
  simp only [Ja0021 hy hy0, Ja1021 hy hy0, Iv002, Iv102]
  simp only [CF0021, CF1021, CF1201, cL1_y hx1, Li2_y hx1, Li3_y, L3p_y hx1, log_one_sub]
  have hx2 : x - 1 ≠ 0 := by linarith
  field_simp
  ring


/-! ## the twenty `J` on every `x < 1`, `x ≠ 0` -/

lemma Jx0001 {x : ℝ} (hx : x < 1) (h0 : x ≠ 0) : S05.Jint 0 0 0 1 x = CF0001 x := by
  rcases lt_or_ge x (1 / 2) with h | h
  · exact Jb0001 h h0
  · exact Ja0001 (abs_lt.2 ⟨by linarith, hx⟩) h0

lemma Jx0002 {x : ℝ} (hx : x < 1) (h0 : x ≠ 0) : S05.Jint 0 0 0 2 x = CF0002 x := by
  rcases lt_or_ge x (1 / 2) with h | h
  · exact Jb0002 h h0
  · exact Ja0002 (abs_lt.2 ⟨by linarith, hx⟩) h0

lemma Jx0003 {x : ℝ} (hx : x < 1) (h0 : x ≠ 0) : S05.Jint 0 0 0 3 x = CF0003 x := by
  rcases lt_or_ge x (1 / 2) with h | h
  · exact Jb0003 h h0
  · exact Ja0003 (abs_lt.2 ⟨by linarith, hx⟩) h0

lemma Jx0011 {x : ℝ} (hx : x < 1) (h0 : x ≠ 0) : S05.Jint 0 0 1 1 x = CF0011 x := by
  rcases lt_or_ge x (1 / 2) with h | h
  · exact Jb0011 h h0
  · exact Ja0011 (abs_lt.2 ⟨by linarith, hx⟩) h0

lemma Jx0012 {x : ℝ} (hx : x < 1) (h0 : x ≠ 0) : S05.Jint 0 0 1 2 x = CF0012 x := by
  rcases lt_or_ge x (1 / 2) with h | h
  · exact Jb0012 h h0
  · exact Ja0012 (abs_lt.2 ⟨by linarith, hx⟩) h0

lemma Jx0021 {x : ℝ} (hx : x < 1) (h0 : x ≠ 0) : S05.Jint 0 0 2 1 x = CF0021 x := by
  rcases lt_or_ge x (1 / 2) with h | h
  · exact Jb0021 h h0
  · exact Ja0021 (abs_lt.2 ⟨by linarith, hx⟩) h0

lemma Jx0101 {x : ℝ} (hx : x < 1) (h0 : x ≠ 0) : S05.Jint 0 1 0 1 x = CF0101 x := by
  rcases lt_or_ge x (1 / 2) with h | h
  · exact Jb0101 h h0
  · exact Ja0101 (abs_lt.2 ⟨by linarith, hx⟩) h0

lemma Jx0102 {x : ℝ} (hx : x < 1) (h0 : x ≠ 0) : S05.Jint 0 1 0 2 x = CF0102 x := by
  rcases lt_or_ge x (1 / 2) with h | h
  · exact Jb0102 h h0
  · exact Ja0102 (abs_lt.2 ⟨by linarith, hx⟩) h0

lemma Jx0111 {x : ℝ} (hx : x < 1) (h0 : x ≠ 0) : S05.Jint 0 1 1 1 x = CF0111 x := by
  rcases lt_or_ge x (1 / 2) with h | h
  · exact Jb0111 h h0
  · exact Ja0111 (abs_lt.2 ⟨by linarith, hx⟩) h0

lemma Jx0201 {x : ℝ} (hx : x < 1) (h0 : x ≠ 0) : S05.Jint 0 2 0 1 x = CF0201 x := by
  rcases lt_or_ge x (1 / 2) with h | h
  · exact Jb0201 h h0
  · exact Ja0201 (abs_lt.2 ⟨by linarith, hx⟩) h0

lemma Jx1001 {x : ℝ} (hx : x < 1) (h0 : x ≠ 0) : S05.Jint 1 0 0 1 x = CF1001 x := by
  rcases lt_or_ge x (1 / 2) with h | h
  · exact Jb1001 h h0
  · exact Ja1001 (abs_lt.2 ⟨by linarith, hx⟩) h0

lemma Jx1002 {x : ℝ} (hx : x < 1) (h0 : x ≠ 0) : S05.Jint 1 0 0 2 x = CF1002 x := by
  rcases lt_or_ge x (1 / 2) with h | h
  · exact Jb1002 h h0
  · exact Ja1002 (abs_lt.2 ⟨by linarith, hx⟩) h0

lemma Jx1003 {x : ℝ} (hx : x < 1) (h0 : x ≠ 0) : S05.Jint 1 0 0 3 x = CF1003 x := by
  rcases lt_or_ge x (1 / 2) with h | h
  · exact Jb1003 h h0
  · exact Ja1003 (abs_lt.2 ⟨by linarith, hx⟩) h0

lemma Jx1011 {x : ℝ} (hx : x < 1) (h0 : x ≠ 0) : S05.Jint 1 0 1 1 x = CF1011 x := by
  rcases lt_or_ge x (1 / 2) with h | h
  · exact Jb1011 h h0
  · exact Ja1011 (abs_lt.2 ⟨by linarith, hx⟩) h0

lemma Jx1012 {x : ℝ} (hx : x < 1) (h0 : x ≠ 0) : S05.Jint 1 0 1 2 x = CF1012 x := by
  rcases lt_or_ge x (1 / 2) with h | h
  · exact Jb1012 h h0
  · exact Ja1012 (abs_lt.2 ⟨by linarith, hx⟩) h0

lemma Jx1021 {x : ℝ} (hx : x < 1) (h0 : x ≠ 0) : S05.Jint 1 0 2 1 x = CF1021 x := by
  rcases lt_or_ge x (1 / 2) with h | h
  · exact Jb1021 h h0
  · exact Ja1021 (abs_lt.2 ⟨by linarith, hx⟩) h0

lemma Jx1101 {x : ℝ} (hx : x < 1) (h0 : x ≠ 0) : S05.Jint 1 1 0 1 x = CF1101 x := by
  rcases lt_or_ge x (1 / 2) with h | h
  · exact Jb1101 h h0
  · exact Ja1101 (abs_lt.2 ⟨by linarith, hx⟩) h0

lemma Jx1102 {x : ℝ} (hx : x < 1) (h0 : x ≠ 0) : S05.Jint 1 1 0 2 x = CF1102 x := by
  rcases lt_or_ge x (1 / 2) with h | h
  · exact Jb1102 h h0
  · exact Ja1102 (abs_lt.2 ⟨by linarith, hx⟩) h0

lemma Jx1111 {x : ℝ} (hx : x < 1) (h0 : x ≠ 0) : S05.Jint 1 1 1 1 x = CF1111 x := by
  rcases lt_or_ge x (1 / 2) with h | h
  · exact Jb1111 h h0
  · exact Ja1111 (abs_lt.2 ⟨by linarith, hx⟩) h0

lemma Jx1201 {x : ℝ} (hx : x < 1) (h0 : x ≠ 0) : S05.Jint 1 2 0 1 x = CF1201 x := by
  rcases lt_or_ge x (1 / 2) with h | h
  · exact Jb1201 h h0
  · exact Ja1201 (abs_lt.2 ⟨by linarith, hx⟩) h0

/-! ## 4. the third `p`-derivative of the length element for every `r` -/

/-- `F_m(t) = Σ_j w_jα_j log^m α_j` at `p = 1` on an edge with data `B, r, l_0, l_1, l_B`
(`A_0 = 1 - Br`, `A_1 = 1 - B`, `1 + κt = 1 - (1-r)t`). -/
noncomputable def gF (B r l0 l1 lB : ℝ) (m : ℕ) (t : ℝ) : ℝ :=
  (1 - B * r) * t * (l0 + Real.log t) ^ m + (1 - B) * (1 - t) * (l1 + Real.log (1 - t)) ^ m
    + B * (1 - (1 - r) * t) * (lB + Real.log (1 - (1 - r) * t)) ^ m

/-- `S_m(t) = Σ_j w_jβ_j log^m α_j` at `p = 1`. -/
noncomputable def gS (B r l0 l1 lB : ℝ) (m : ℕ) (t : ℝ) : ℝ :=
  (1 - B * r) * (l0 + Real.log t) ^ m - (1 - B) * (l1 + Real.log (1 - t)) ^ m
    + B * (r - 1) * (lB + Real.log (1 - (1 - r) * t)) ^ m

/-- the edge parameters in the order `B, r, l_0, l_1, l_B`. -/
noncomputable def eF (H : S05.Hex) (k : Fin 3) (m : ℕ) (t : ℝ) : ℝ :=
  gF (S05.Bv H k) (S05.rr H k) (S05.l0 H k) (S05.l1 H k) (S05.lB H k) m t
noncomputable def eS (H : S05.Hex) (k : Fin 3) (m : ℕ) (t : ℝ) : ℝ :=
  gS (S05.Bv H k) (S05.rr H k) (S05.l0 H k) (S05.l1 H k) (S05.lB H k) m t

lemma edge_facts (H : S05.Hex) (hH : S05.IsHex H) (k : Fin 3) :
    S05.A0 H k = 1 - S05.Bv H k * S05.rr H k ∧ S05.A1 H k = 1 - S05.Bv H k ∧
    0 < S05.Bv H k ∧ S05.Bv H k < 1 ∧ 0 < S05.A0 H k ∧ 0 < S05.A1 H k ∧ 0 < S05.rr H k := by
  obtain ⟨-, h1, h2, h3, h4, h5, h6, h7⟩ := S05.proof_8 H hH k
  exact ⟨by linarith, by linarith, h3, h4, h5, h6, h7⟩

lemma aff_pos {r t : ℝ} (hr : 0 < r) (ht0 : 0 ≤ t) (ht1 : t ≤ 1) : 0 < 1 - (1 - r) * t := by
  rcases eq_or_lt_of_le ht0 with h | h
  · subst h; norm_num
  · nlinarith [mul_pos hr h]

lemma eα_val (H : S05.Hex) (k : Fin 3) (t : ℝ) :
    S05.eα H k t 0 = S05.A0 H k * t / H.w k ∧ S05.eα H k t 1 = S05.A1 H k * (1 - t) / H.w (k + 1) ∧
    S05.eα H k t 2 = S05.Bv H k * (1 + (S05.rr H k - 1) * t) / H.w (k + 2) := by
  simp [S05.eα]

lemma eβ_val (H : S05.Hex) (k : Fin 3) :
    S05.eβ H k 0 = S05.A0 H k / H.w k ∧ S05.eβ H k 1 = -(S05.A1 H k / H.w (k + 1)) ∧
    S05.eβ H k 2 = S05.Bv H k * (S05.rr H k - 1) / H.w (k + 2) := by
  simp [S05.eβ]

lemma ew_val (H : S05.Hex) (k : Fin 3) :
    S05.ew H k 0 = H.w k ∧ S05.ew H k 1 = H.w (k + 1) ∧ S05.ew H k 2 = H.w (k + 2) := by
  simp [S05.ew]

lemma log_α (H : S05.Hex) (hH : S05.IsHex H) (k : Fin 3) {t : ℝ} (ht : t ∈ Ioo (0:ℝ) 1) :
    Real.log (S05.eα H k t 0) = S05.l0 H k + Real.log t ∧
    Real.log (S05.eα H k t 1) = S05.l1 H k + Real.log (1 - t) ∧
    Real.log (S05.eα H k t 2) = S05.lB H k + Real.log (1 - (1 - S05.rr H k) * t) := by
  obtain ⟨-, -, hB, -, hA0, hA1, hr⟩ := edge_facts H hH k
  obtain ⟨a0, a1, a2⟩ := eα_val H k t
  have hw := hH.1
  have ht0 := ht.1
  have ht1 : 0 < 1 - t := by linarith [ht.2]
  have h2 := aff_pos hr ht.1.le ht.2.le
  refine ⟨?_, ?_, ?_⟩
  · rw [a0, show S05.A0 H k * t / H.w k = (S05.A0 H k / H.w k) * t by ring,
      Real.log_mul (div_pos hA0 (hw k)).ne' ht0.ne']; rfl
  · rw [a1, show S05.A1 H k * (1 - t) / H.w (k + 1) = (S05.A1 H k / H.w (k + 1)) * (1 - t) by ring,
      Real.log_mul (div_pos hA1 (hw _)).ne' ht1.ne']; rfl
  · rw [a2, show S05.Bv H k * (1 + (S05.rr H k - 1) * t) / H.w (k + 2) =
      (S05.Bv H k / H.w (k + 2)) * (1 - (1 - S05.rr H k) * t) by ring,
      Real.log_mul (div_pos hB (hw _)).ne' h2.ne']; rfl

lemma Ed_F (H : S05.Hex) (hH : S05.IsHex H) (k : Fin 3) {t : ℝ} (ht : t ∈ Ioo (0:ℝ) 1) (m : ℕ) :
    S05.Ed (S05.ew H k) (S05.eα H k t) m 1 = eF H k m t := by
  obtain ⟨eA0, eA1, -, -, -, -, -⟩ := edge_facts H hH k
  have hα := S05.eα_pos H hH k ht
  obtain ⟨l0', l1', l2'⟩ := log_α H hH k ht
  obtain ⟨a0, a1, a2⟩ := eα_val H k t
  obtain ⟨w0, w1, w2⟩ := ew_val H k
  have hw := hH.1
  simp only [S05.Ed, Fin.sum_univ_three, mul_one]
  rw [Real.exp_log (hα 0), Real.exp_log (hα 1), Real.exp_log (hα 2), l0', l1', l2', w0, w1, w2,
    a0, a1, a2, eA0, eA1]
  unfold eF gF
  have := (hw k).ne'; have := (hw (k + 1)).ne'; have := (hw (k + 2)).ne'
  field_simp
  ring

lemma Ed_S (H : S05.Hex) (hH : S05.IsHex H) (k : Fin 3) {t : ℝ} (ht : t ∈ Ioo (0:ℝ) 1) (m : ℕ) :
    S05.Ed (fun j => S05.ew H k j * S05.eβ H k j / S05.eα H k t j) (S05.eα H k t) m 1 = eS H k m t := by
  obtain ⟨eA0, eA1, -, -, -, -, -⟩ := edge_facts H hH k
  have hα := S05.eα_pos H hH k ht
  obtain ⟨l0', l1', l2'⟩ := log_α H hH k ht
  obtain ⟨b0, b1, b2⟩ := eβ_val H k
  obtain ⟨w0, w1, w2⟩ := ew_val H k
  have hw := hH.1
  simp only [S05.Ed, Fin.sum_univ_three, mul_one]
  rw [Real.exp_log (hα 0), Real.exp_log (hα 1), Real.exp_log (hα 2), l0', l1', l2', w0, w1, w2,
    b0, b1, b2]
  have := (hα 0).ne'; have := (hα 1).ne'; have := (hα 2).ne'
  obtain ⟨-, -, hB, -, hA0, hA1, hr⟩ := edge_facts H hH k
  have := (hw k).ne'; have := (hw (k + 1)).ne'; have := (hw (k + 2)).ne'
  have ht0 : t ≠ 0 := ht.1.ne'
  have ht1 : 1 - t ≠ 0 := by linarith [ht.2]
  have h2 : 1 + (S05.rr H k - 1) * t ≠ 0 := by
    have := aff_pos hr ht.1.le ht.2.le; intro h; apply this.ne'; linarith
  have := hA0.ne'; have := hA1.ne'; have := hB.ne'
  unfold eS gS
  field_simp
  simp only [eA0, eA1]
  ring

lemma Gw3_gen (H : S05.Hex) (hH : S05.IsHex H) (k : Fin 3) {t : ℝ} (ht : t ∈ Ioo (0:ℝ) 1) :
    S05.Gw H k 3 1 t = S05.Len H k * (S05.gg3 (eF H k 1 t) (eF H k 2 t) (eF H k 3 t) -
        (S05.q0 H k / S05.Lsq H k + t) *
          S05.WW (eF H k 1 t) (eF H k 2 t) (eS H k 1 t) (eS H k 2 t) (eS H k 3 t)) +
      (S05.V02 H k - S05.q0 H k ^ 2 / S05.Lsq H k) / S05.Len H k *
        S05.HH (S05.q0 H k / S05.Lsq H k + t) (eF H k 1 t) (eS H k 1 t) (eS H k 2 t) := by
  obtain ⟨eA0, eA1, -, -, -, -, -⟩ := edge_facts H hH k
  have hα := S05.eα_pos H hH k ht
  obtain ⟨a0, a1, a2⟩ := eα_val H k t
  obtain ⟨b0, b1, b2⟩ := eβ_val H k
  obtain ⟨w0, w1, w2⟩ := ew_val H k
  have hw := hH.1
  have hL := S05.Lsq_pos H hH k
  have hLn : S05.Len H k ^ 2 = S05.Lsq H k := Real.sq_sqrt hL.le
  have hLp : 0 < S05.Len H k := Real.sqrt_pos.2 hL
  have hQ := fun q => S05.cQ_pos (α := S05.eα H k t) (β := S05.eβ H k) (Y0 := S05.ept H k t 0)
    (Y1 := S05.ept H k t 1) (E0 := S05.eE H k 0) (E1 := S05.eE H k 1) (S05.ew_pos H hH k)
    (S05.edge_det_ne H hH k t) q
  have hFp := fun q => S05.cF_pos (α := S05.eα H k t) (S05.ew_pos H hH k) q
  have hE : S05.eE H k 0 ^ 2 + S05.eE H k 1 ^ 2 = S05.Lsq H k := by
    simp [S05.Lsq, S05.eE, X04.dot]; ring
  have hF0 : S05.cF (S05.ew H k) (S05.eα H k t) 1 = 1 := by
    simp only [S05.cF, Fin.sum_univ_three, mul_one]
    rw [Real.exp_log (hα 0), Real.exp_log (hα 1), Real.exp_log (hα 2), w0, w1, w2, a0, a1, a2,
      eA0, eA1]
    have := (hw k).ne'; have := (hw (k + 1)).ne'; have := (hw (k + 2)).ne'
    field_simp
    ring
  have hS0 : S05.cS (S05.ew H k) (S05.eα H k t) (S05.eβ H k) 1 = 0 := by
    simp only [S05.cS, Fin.sum_univ_three, mul_one]
    rw [Real.exp_log (hα 0), Real.exp_log (hα 1), Real.exp_log (hα 2),
      div_mul_cancel₀ _ (hα 0).ne', div_mul_cancel₀ _ (hα 1).ne', div_mul_cancel₀ _ (hα 2).ne',
      w0, w1, w2, b0, b1, b2, eA0, eA1]
    have := (hw k).ne'; have := (hw (k + 1)).ne'; have := (hw (k + 2)).ne'
    field_simp
    ring
  have key := S05.G_iter (S05.cF_J3 _ _) (S05.cS_J3 _ _ _) hF0 hS0 hFp hQ (by rw [hE]; exact hL)
  unfold S05.Gw
  rw [iteratedDerivWithin_of_isOpen isOpen_Ioi (by norm_num : (1:ℝ) ∈ Set.Ioi 0)]
  have hG : (fun q => S05.Gf H k q t) = S05.cG (S05.ew H k) (S05.eα H k t) (S05.eβ H k)
      (S05.eE H k 0) (S05.eE H k 1) (S05.ept H k t 0) (S05.ept H k t 1) := by
    funext q; rfl
  rw [hG, key]
  rw [Ed_F H hH k ht 1, Ed_F H hH k ht 2, Ed_F H hH k ht 3, Ed_S H hH k ht 1, Ed_S H hH k ht 2,
    Ed_S H hH k ht 3, hE]
  have hq : S05.eE H k 0 * S05.ept H k t 0 + S05.eE H k 1 * S05.ept H k t 1 =
      S05.q0 H k + t * S05.Lsq H k := by
    simp [S05.eE, S05.ept, S05.q0, S05.Lsq, X04.dot]; ring
  have hY : S05.ept H k t 0 ^ 2 + S05.ept H k t 1 ^ 2 =
      S05.V02 H k + 2 * t * S05.q0 H k + t ^ 2 * S05.Lsq H k := by
    simp [S05.ept, S05.q0, S05.Lsq, S05.V02, X04.dot]; ring
  rw [hq, hY]
  have hsq : Real.sqrt (S05.Lsq H k) = S05.Len H k := rfl
  rw [hsq]
  have hL0 : S05.Lsq H k ≠ 0 := hL.ne'
  have hLn0 : S05.Len H k ≠ 0 := hLp.ne'
  rw [← hLn]
  field_simp
  ring

/-! ## 5. the `L`-part is an exact derivative -/

lemma hasDerivAt_gF (B r l0 l1 lB : ℝ) (hr : 0 < r) (m : ℕ) {t : ℝ} (ht : t ∈ Ioo (0:ℝ) 1) :
    HasDerivAt (gF B r l0 l1 lB m) (gS B r l0 l1 lB m t + m * gS B r l0 l1 lB (m - 1) t) t := by
  have ht0 : t ≠ 0 := ht.1.ne'
  have ht1 : 1 - t ≠ 0 := by linarith [ht.2]
  have h2 : 1 - (1 - r) * t ≠ 0 := (aff_pos hr ht.1.le ht.2.le).ne'
  have dP : HasDerivAt (fun s => l0 + Real.log s) (1 / t) t := by
    have := (Real.hasDerivAt_log ht0).const_add l0
    refine this.congr_deriv ?_; field_simp
  have dM : HasDerivAt (fun s => l1 + Real.log (1 - s)) (-(1 / (1 - t))) t := by
    have := (((hasDerivAt_id' t).const_sub 1).log ht1).const_add l1
    refine this.congr_deriv ?_
    field_simp
  have dW : HasDerivAt (fun s => lB + Real.log (1 - (1 - r) * s)) (-(1 - r) / (1 - (1 - r) * t)) t := by
    have := ((((hasDerivAt_id' t).const_mul (1 - r)).const_sub 1).log h2).const_add lB
    refine this.congr_deriv ?_
    field_simp
  have := (((((hasDerivAt_id' t).const_mul (1 - B * r)).mul (dP.pow m))).add
    ((((hasDerivAt_id' t).const_sub 1).const_mul (1 - B)).mul (dM.pow m))).add
    (((((hasDerivAt_id' t).const_mul (1 - r)).const_sub 1).const_mul B).mul (dW.pow m))
  unfold gF
  refine this.congr_deriv ?_
  unfold gS
  simp only [Pi.pow_apply, id]
  rcases m with _ | m
  · norm_num
  · simp only [Nat.add_sub_cancel, Nat.cast_add, Nat.cast_one, pow_succ]
    generalize 1 - (1 - r) * t = u at *
    generalize 1 - t = v at *
    field_simp
    ring

lemma gS_zero (B r l0 l1 lB t : ℝ) : gS B r l0 l1 lB 0 t = 0 := by
  unfold gS; ring

lemma cont_tlog (c : ℝ) (m : ℕ) : ContinuousOn (fun t : ℝ => t * (c + Real.log t) ^ m) (Icc 0 1) := by
  have e : (fun t : ℝ => t * (c + Real.log t) ^ m) = fun t =>
      ∑ j ∈ range (m + 1), (c ^ j * (m.choose j : ℝ)) * (t * Real.log t ^ (m - j)) := by
    funext t
    rw [add_pow, Finset.mul_sum]
    refine Finset.sum_congr rfl fun j _ => ?_
    ring
  rw [e]
  refine continuousOn_finsetSum _ fun j _ => ?_
  have hj : ContinuousOn (fun t : ℝ => t * Real.log t ^ (m - j)) (Icc 0 1) := by
    rcases Nat.eq_zero_or_pos (m - j) with h | h
    · simp only [h, pow_zero, mul_one]; exact continuousOn_id
    · have := S05.G_cont 0 (m - j) h.ne'
      simp only [zero_add, pow_one] at this
      exact this
  exact ContinuousOn.mul (f := fun _ : ℝ => (c ^ j * (m.choose j : ℝ))) continuousOn_const hj

lemma cont_tlog1 (c : ℝ) (m : ℕ) :
    ContinuousOn (fun t : ℝ => (1 - t) * (c + Real.log (1 - t)) ^ m) (Icc 0 1) := by
  have hmaps : Set.MapsTo (fun t : ℝ => 1 - t) (Icc 0 1) (Icc 0 1) :=
    fun t ht => ⟨by linarith [ht.2], by linarith [ht.1]⟩
  exact (cont_tlog c m).comp (continuousOn_const.sub continuousOn_id) hmaps

lemma cont_gF (B r l0 l1 lB : ℝ) (hr : 0 < r) (m : ℕ) : ContinuousOn (gF B r l0 l1 lB m) (Icc 0 1) := by
  have hu : ContinuousOn (fun t : ℝ => 1 - (1 - r) * t) (Icc 0 1) :=
    continuousOn_const.sub (continuousOn_const.mul continuousOn_id)
  have hl : ContinuousOn (fun t => Real.log (1 - (1 - r) * t)) (Icc 0 1) :=
    hu.log (fun t ht => (aff_pos hr ht.1 ht.2).ne')
  have h3 : ContinuousOn (fun t : ℝ => B * (1 - (1 - r) * t) * (lB + Real.log (1 - (1 - r) * t)) ^ m)
      (Icc 0 1) := (continuousOn_const.mul hu).mul ((continuousOn_const.add hl).pow m)
  have h1 := (cont_tlog l0 m).const_smul (1 - B * r)
  have h2 := (cont_tlog1 l1 m).const_smul (1 - B)
  have := (h1.add h2).add h3
  refine this.congr fun t _ => ?_
  simp only [gF, Pi.add_apply, Pi.smul_apply, smul_eq_mul]
  ring

lemma gF_one (B r l0 l1 lB : ℝ) (m : ℕ) :
    gF B r l0 l1 lB m 1 = (1 - B * r) * l0 ^ m + B * r * (lB + Real.log r) ^ m := by
  unfold gF; simp

lemma gF_zero (B r l0 l1 lB : ℝ) (m : ℕ) :
    gF B r l0 l1 lB m 0 = (1 - B) * l1 ^ m + B * lB ^ m := by
  unfold gF; simp

/-- the `L`-part: `∫_0^1 (g''' - aW) dt = [(a_0+t) g'''(F(t))]_0^1`. -/
lemma Lpart_gen (B r l0 l1 lB a0 : ℝ) (hr : 0 < r)
    (hi : IntervalIntegrable (fun t => S05.gg3 (gF B r l0 l1 lB 1 t) (gF B r l0 l1 lB 2 t)
      (gF B r l0 l1 lB 3 t) - (a0 + t) * S05.WW (gF B r l0 l1 lB 1 t) (gF B r l0 l1 lB 2 t)
        (gS B r l0 l1 lB 1 t) (gS B r l0 l1 lB 2 t) (gS B r l0 l1 lB 3 t)) volume 0 1) :
    (∫ t in (0:ℝ)..1, (S05.gg3 (gF B r l0 l1 lB 1 t) (gF B r l0 l1 lB 2 t)
      (gF B r l0 l1 lB 3 t) - (a0 + t) * S05.WW (gF B r l0 l1 lB 1 t) (gF B r l0 l1 lB 2 t)
        (gS B r l0 l1 lB 1 t) (gS B r l0 l1 lB 2 t) (gS B r l0 l1 lB 3 t))) =
      (a0 + 1) * S05.gg3 (gF B r l0 l1 lB 1 1) (gF B r l0 l1 lB 2 1) (gF B r l0 l1 lB 3 1)
        - a0 * S05.gg3 (gF B r l0 l1 lB 1 0) (gF B r l0 l1 lB 2 0) (gF B r l0 l1 lB 3 0) := by
  set F := gF B r l0 l1 lB
  set τ : ℝ → ℝ := fun t => (a0 + t) * (-F 3 t + 6 * (F 1 t * F 2 t) - 6 * F 1 t ^ 3 +
    3 * F 2 t - 9 * F 1 t ^ 2 - 6 * F 1 t) with hτ
  have hc : ContinuousOn τ (Icc 0 1) := by
    have c1 := cont_gF B r l0 l1 lB hr 1
    have c2 := cont_gF B r l0 l1 lB hr 2
    have c3 := cont_gF B r l0 l1 lB hr 3
    have e6 : ContinuousOn (fun t : ℝ => a0 + t) (Icc 0 1) :=
      (continuous_const.add continuous_id).continuousOn
    exact e6.mul ((((((c3.neg.add ((c1.mul c2).const_smul (6:ℝ))).sub
      ((c1.pow 3).const_smul (6:ℝ))).add (c2.const_smul (3:ℝ))).sub
      ((c1.pow 2).const_smul (9:ℝ))).sub (c1.const_smul (6:ℝ))))
  have hd : ∀ t ∈ Ioo (0:ℝ) 1, HasDerivAt τ (S05.gg3 (F 1 t) (F 2 t) (F 3 t) -
      (a0 + t) * S05.WW (F 1 t) (F 2 t) (gS B r l0 l1 lB 1 t) (gS B r l0 l1 lB 2 t)
        (gS B r l0 l1 lB 3 t)) t := by
    intro t ht
    have d1 := hasDerivAt_gF B r l0 l1 lB hr 1 ht
    have d2 := hasDerivAt_gF B r l0 l1 lB hr 2 ht
    have d3 := hasDerivAt_gF B r l0 l1 lB hr 3 ht
    have := ((hasDerivAt_id' t).const_add a0).mul
      ((((((d3.neg.add ((d1.mul d2).const_mul 6)).sub ((d1.pow 3).const_mul 6)).add
        (d2.const_mul 3)).sub ((d1.pow 2).const_mul 9)).sub (d1.const_mul 6)))
    refine this.congr_deriv ?_
    rw [S05.gg3_eq]
    have e0 : gS B r l0 l1 lB 0 t = 0 := gS_zero _ _ _ _ _ _
    simp only [S05.WW, S05.gg2, S05.gg1, S05.kp1, S05.kp2, S05.uu2, S05.uu3, Nat.cast_ofNat,
      Nat.cast_one, e0]
    norm_num
    ring
  rw [intervalIntegral.integral_eq_sub_of_hasDerivAt_of_le zero_le_one hc hd hi]
  simp only [hτ, S05.gg3_eq]
  ring

/-! ## 6. the `h²`-part as forty monomials -/

/- the coefficients of the forty monomials `t^n log^a t log^b(1-t) log^c(1-(1-r)t)` of
`HH(a_0+t, F_1, S_1, S_2)` (computed by expansion; checked by `ring` below). -/
set_option maxHeartbeats 2000000 in
noncomputable def kc0 (B r l0 l1 lB a0 : ℝ) : ℝ :=
  (((-27) * (B) ^ 3 * (l1) ^ 2) + ((-27) * (B) ^ 3 * (lB) ^ 2) + ((-24) * B * (l1) ^ 2) + ((-12) * B * (l0) ^ 2) + (9 * (B) ^ 2 * (lB) ^ 2) + (51 * (B) ^ 2 * (l1) ^ 2) + ((-60) * l1 * lB * (B) ^ 2) + ((-36) * l0 * l1 * (B) ^ 2) + ((-33) * r * (B) ^ 2 * (l1) ^ 2) + ((-27) * (B) ^ 3 * (lB) ^ 2 * (r) ^ 2) + ((-18) * r * (B) ^ 2 * (lB) ^ 2) + ((-9) * B * a0 * (l0) ^ 2) + ((-9) * B * a0 * (l1) ^ 2) + ((-9) * a0 * (B) ^ 3 * (l1) ^ 2) + ((-9) * a0 * (B) ^ 3 * (lB) ^ 2) + ((-9) * (B) ^ 3 * (l0) ^ 2 * (r) ^ 2) + ((-6) * B * l0 * lB) + ((-3) * (B) ^ 2 * (l0) ^ 2 * (r) ^ 2) + (3 * B * r * (l0) ^ 2) + (6 * B * l1 * lB) + (9 * (B) ^ 2 * (lB) ^ 2 * (r) ^ 2) + (15 * B * r * (l1) ^ 2) + (18 * a0 * (B) ^ 2 * (l1) ^ 2) + (18 * r * (B) ^ 3 * (l1) ^ 2) + (21 * r * (B) ^ 2 * (l0) ^ 2) + (36 * B * l0 * l1) + (36 * l0 * lB * (B) ^ 2) + (54 * l1 * lB * (B) ^ 3) + (54 * r * (B) ^ 3 * (lB) ^ 2) + ((-72) * l1 * lB * r * (B) ^ 3) + ((-36) * l0 * lB * r * (B) ^ 3) + ((-30) * l0 * lB * r * (B) ^ 2) + ((-27) * a0 * (B) ^ 3 * (lB) ^ 2 * (r) ^ 2) + ((-18) * B * l0 * l1 * r) + ((-18) * a0 * l0 * l1 * (B) ^ 2) + ((-18) * a0 * l1 * lB * (B) ^ 2) + ((-18) * a0 * r * (B) ^ 2 * (l1) ^ 2) + ((-18) * a0 * (B) ^ 2 * (l0) ^ 2 * (r) ^ 2) + ((-18) * l0 * l1 * r * (B) ^ 2) + ((-18) * l0 * l1 * (B) ^ 3 * (r) ^ 2) + ((-18) * l1 * lB * (B) ^ 2 * (r) ^ 2) + ((-9) * a0 * (B) ^ 3 * (l0) ^ 2 * (r) ^ 2) + ((-6) * B * l1 * lB * r) + ((-6) * l0 * lB * (B) ^ 2 * (r) ^ 2) + (6 * B * l0 * lB * r) + (9 * B * a0 * r * (l0) ^ 2) + (9 * B * a0 * r * (l1) ^ 2) + (9 * a0 * r * (B) ^ 3 * (l1) ^ 2) + (9 * a0 * (B) ^ 3 * (l0) ^ 2 * (r) ^ 3) + (9 * a0 * (B) ^ 3 * (lB) ^ 2 * (r) ^ 3) + (18 * B * a0 * l0 * l1) + (18 * a0 * l0 * lB * (B) ^ 2) + (18 * a0 * l1 * lB * (B) ^ 3) + (18 * a0 * r * (B) ^ 2 * (l0) ^ 2) + (18 * l0 * l1 * (B) ^ 2 * (r) ^ 2) + (18 * l1 * lB * (B) ^ 3 * (r) ^ 2) + (27 * a0 * r * (B) ^ 3 * (lB) ^ 2) + (36 * l0 * l1 * r * (B) ^ 3) + (36 * l0 * lB * (B) ^ 3 * (r) ^ 2) + (78 * l1 * lB * r * (B) ^ 2) + ((-36) * a0 * l0 * lB * r * (B) ^ 2) + ((-36) * a0 * l1 * lB * r * (B) ^ 3) + ((-18) * B * a0 * l0 * l1 * r) + ((-18) * a0 * l0 * l1 * (B) ^ 3 * (r) ^ 2) + ((-18) * a0 * l0 * lB * r * (B) ^ 3) + ((-18) * a0 * l0 * lB * (B) ^ 3 * (r) ^ 3) + ((-18) * a0 * l1 * lB * (B) ^ 2 * (r) ^ 2) + (18 * a0 * l0 * l1 * r * (B) ^ 3) + (18 * a0 * l0 * l1 * (B) ^ 2 * (r) ^ 2) + (18 * a0 * l0 * lB * (B) ^ 2 * (r) ^ 2) + (18 * a0 * l1 * lB * (B) ^ 3 * (r) ^ 2) + (36 * a0 * l0 * lB * (B) ^ 3 * (r) ^ 2) + (36 * a0 * l1 * lB * r * (B) ^ 2))

set_option maxHeartbeats 2000000 in
noncomputable def kc1 (B r l0 l1 lB a0 : ℝ) : ℝ :=
  (((-30) * l1 * (B) ^ 2) + ((-27) * lB * (B) ^ 3) + ((-3) * B * l0) + (3 * B * l1) + (9 * lB * (B) ^ 2) + (18 * l0 * (B) ^ 2) + (27 * l1 * (B) ^ 3) + ((-36) * l1 * r * (B) ^ 3) + ((-27) * lB * (B) ^ 3 * (r) ^ 2) + ((-18) * l0 * r * (B) ^ 3) + ((-18) * lB * r * (B) ^ 2) + ((-15) * l0 * r * (B) ^ 2) + ((-9) * a0 * l1 * (B) ^ 2) + ((-9) * a0 * lB * (B) ^ 3) + ((-9) * l1 * (B) ^ 2 * (r) ^ 2) + ((-3) * B * l1 * r) + ((-3) * l0 * (B) ^ 2 * (r) ^ 2) + (3 * B * l0 * r) + (9 * a0 * l0 * (B) ^ 2) + (9 * a0 * l1 * (B) ^ 3) + (9 * l1 * (B) ^ 3 * (r) ^ 2) + (9 * lB * (B) ^ 2 * (r) ^ 2) + (18 * l0 * (B) ^ 3 * (r) ^ 2) + (39 * l1 * r * (B) ^ 2) + (54 * lB * r * (B) ^ 3) + ((-27) * a0 * lB * (B) ^ 3 * (r) ^ 2) + ((-18) * a0 * l0 * r * (B) ^ 2) + ((-18) * a0 * l1 * r * (B) ^ 3) + ((-9) * a0 * l0 * r * (B) ^ 3) + ((-9) * a0 * l0 * (B) ^ 3 * (r) ^ 3) + ((-9) * a0 * l1 * (B) ^ 2 * (r) ^ 2) + (9 * a0 * l0 * (B) ^ 2 * (r) ^ 2) + (9 * a0 * l1 * (B) ^ 3 * (r) ^ 2) + (9 * a0 * lB * (B) ^ 3 * (r) ^ 3) + (18 * a0 * l0 * (B) ^ 3 * (r) ^ 2) + (18 * a0 * l1 * r * (B) ^ 2) + (27 * a0 * lB * r * (B) ^ 3))

set_option maxHeartbeats 2000000 in
noncomputable def kc2 (B r l0 l1 lB a0 : ℝ) : ℝ :=
  (((-9) * (B) ^ 3) + (3 * (B) ^ 2) + ((-9) * (B) ^ 3 * (r) ^ 2) + ((-6) * r * (B) ^ 2) + ((-3) * a0 * (B) ^ 3) + (3 * (B) ^ 2 * (r) ^ 2) + (18 * r * (B) ^ 3) + ((-9) * a0 * (B) ^ 3 * (r) ^ 2) + (3 * a0 * (B) ^ 3 * (r) ^ 3) + (9 * a0 * r * (B) ^ 3))

set_option maxHeartbeats 2000000 in
noncomputable def kc3 (B r l0 l1 lB a0 : ℝ) : ℝ :=
  (((-60) * lB * (B) ^ 2) + ((-54) * l1 * (B) ^ 3) + ((-48) * B * l1) + ((-36) * l0 * (B) ^ 2) + (6 * B * lB) + (36 * B * l0) + (54 * lB * (B) ^ 3) + (102 * l1 * (B) ^ 2) + ((-72) * lB * r * (B) ^ 3) + ((-66) * l1 * r * (B) ^ 2) + ((-18) * B * a0 * l1) + ((-18) * B * l0 * r) + ((-18) * a0 * l0 * (B) ^ 2) + ((-18) * a0 * l1 * (B) ^ 3) + ((-18) * a0 * lB * (B) ^ 2) + ((-18) * l0 * r * (B) ^ 2) + ((-18) * l0 * (B) ^ 3 * (r) ^ 2) + ((-18) * lB * (B) ^ 2 * (r) ^ 2) + ((-6) * B * lB * r) + (18 * B * a0 * l0) + (18 * a0 * lB * (B) ^ 3) + (18 * l0 * (B) ^ 2 * (r) ^ 2) + (18 * lB * (B) ^ 3 * (r) ^ 2) + (30 * B * l1 * r) + (36 * a0 * l1 * (B) ^ 2) + (36 * l0 * r * (B) ^ 3) + (36 * l1 * r * (B) ^ 3) + (78 * lB * r * (B) ^ 2) + ((-36) * a0 * l1 * r * (B) ^ 2) + ((-36) * a0 * lB * r * (B) ^ 3) + ((-18) * B * a0 * l0 * r) + ((-18) * a0 * l0 * (B) ^ 3 * (r) ^ 2) + ((-18) * a0 * lB * (B) ^ 2 * (r) ^ 2) + (18 * B * a0 * l1 * r) + (18 * a0 * l0 * r * (B) ^ 3) + (18 * a0 * l0 * (B) ^ 2 * (r) ^ 2) + (18 * a0 * l1 * r * (B) ^ 3) + (18 * a0 * lB * (B) ^ 3 * (r) ^ 2) + (36 * a0 * lB * r * (B) ^ 2))

set_option maxHeartbeats 2000000 in
noncomputable def kc4 (B r l0 l1 lB a0 : ℝ) : ℝ :=
  (((-30) * (B) ^ 2) + (3 * B) + (27 * (B) ^ 3) + ((-36) * r * (B) ^ 3) + ((-9) * a0 * (B) ^ 2) + ((-9) * (B) ^ 2 * (r) ^ 2) + ((-3) * B * r) + (9 * a0 * (B) ^ 3) + (9 * (B) ^ 3 * (r) ^ 2) + (39 * r * (B) ^ 2) + ((-18) * a0 * r * (B) ^ 3) + ((-9) * a0 * (B) ^ 2 * (r) ^ 2) + (9 * a0 * (B) ^ 3 * (r) ^ 2) + (18 * a0 * r * (B) ^ 2))

set_option maxHeartbeats 2000000 in
noncomputable def kc5 (B r l0 l1 lB a0 : ℝ) : ℝ :=
  (((-27) * (B) ^ 3) + ((-24) * B) + (51 * (B) ^ 2) + ((-33) * r * (B) ^ 2) + ((-9) * B * a0) + ((-9) * a0 * (B) ^ 3) + (15 * B * r) + (18 * a0 * (B) ^ 2) + (18 * r * (B) ^ 3) + ((-18) * a0 * r * (B) ^ 2) + (9 * B * a0 * r) + (9 * a0 * r * (B) ^ 3))

set_option maxHeartbeats 2000000 in
noncomputable def kc6 (B r l0 l1 lB a0 : ℝ) : ℝ :=
  (((-36) * l1 * (B) ^ 2) + ((-24) * B * l0) + ((-6) * B * lB) + (36 * B * l1) + (36 * lB * (B) ^ 2) + ((-36) * lB * r * (B) ^ 3) + ((-30) * lB * r * (B) ^ 2) + ((-18) * B * a0 * l0) + ((-18) * B * l1 * r) + ((-18) * a0 * l1 * (B) ^ 2) + ((-18) * l0 * (B) ^ 3 * (r) ^ 2) + ((-18) * l1 * r * (B) ^ 2) + ((-18) * l1 * (B) ^ 3 * (r) ^ 2) + ((-6) * l0 * (B) ^ 2 * (r) ^ 2) + ((-6) * lB * (B) ^ 2 * (r) ^ 2) + (6 * B * l0 * r) + (6 * B * lB * r) + (18 * B * a0 * l1) + (18 * a0 * lB * (B) ^ 2) + (18 * l1 * (B) ^ 2 * (r) ^ 2) + (36 * l1 * r * (B) ^ 3) + (36 * lB * (B) ^ 3 * (r) ^ 2) + (42 * l0 * r * (B) ^ 2) + ((-36) * a0 * l0 * (B) ^ 2 * (r) ^ 2) + ((-36) * a0 * lB * r * (B) ^ 2) + ((-18) * B * a0 * l1 * r) + ((-18) * a0 * l0 * (B) ^ 3 * (r) ^ 2) + ((-18) * a0 * l1 * (B) ^ 3 * (r) ^ 2) + ((-18) * a0 * lB * r * (B) ^ 3) + ((-18) * a0 * lB * (B) ^ 3 * (r) ^ 3) + (18 * B * a0 * l0 * r) + (18 * a0 * l0 * (B) ^ 3 * (r) ^ 3) + (18 * a0 * l1 * r * (B) ^ 3) + (18 * a0 * l1 * (B) ^ 2 * (r) ^ 2) + (18 * a0 * lB * (B) ^ 2 * (r) ^ 2) + (36 * a0 * l0 * r * (B) ^ 2) + (36 * a0 * lB * (B) ^ 3 * (r) ^ 2))

set_option maxHeartbeats 2000000 in
noncomputable def kc7 (B r l0 l1 lB a0 : ℝ) : ℝ :=
  (((-3) * B) + (18 * (B) ^ 2) + ((-18) * r * (B) ^ 3) + ((-15) * r * (B) ^ 2) + ((-3) * (B) ^ 2 * (r) ^ 2) + (3 * B * r) + (9 * a0 * (B) ^ 2) + (18 * (B) ^ 3 * (r) ^ 2) + ((-18) * a0 * r * (B) ^ 2) + ((-9) * a0 * r * (B) ^ 3) + ((-9) * a0 * (B) ^ 3 * (r) ^ 3) + (9 * a0 * (B) ^ 2 * (r) ^ 2) + (18 * a0 * (B) ^ 3 * (r) ^ 2))

set_option maxHeartbeats 2000000 in
noncomputable def kc8 (B r l0 l1 lB a0 : ℝ) : ℝ :=
  (((-36) * (B) ^ 2) + (36 * B) + ((-18) * B * r) + ((-18) * a0 * (B) ^ 2) + ((-18) * r * (B) ^ 2) + ((-18) * (B) ^ 3 * (r) ^ 2) + (18 * B * a0) + (18 * (B) ^ 2 * (r) ^ 2) + (36 * r * (B) ^ 3) + ((-18) * B * a0 * r) + ((-18) * a0 * (B) ^ 3 * (r) ^ 2) + (18 * a0 * r * (B) ^ 3) + (18 * a0 * (B) ^ 2 * (r) ^ 2))

set_option maxHeartbeats 2000000 in
noncomputable def kc9 (B r l0 l1 lB a0 : ℝ) : ℝ :=
  (((-12) * B) + ((-9) * B * a0) + ((-9) * (B) ^ 3 * (r) ^ 2) + ((-3) * (B) ^ 2 * (r) ^ 2) + (3 * B * r) + (21 * r * (B) ^ 2) + ((-18) * a0 * (B) ^ 2 * (r) ^ 2) + ((-9) * a0 * (B) ^ 3 * (r) ^ 2) + (9 * B * a0 * r) + (9 * a0 * (B) ^ 3 * (r) ^ 3) + (18 * a0 * r * (B) ^ 2))

set_option maxHeartbeats 2000000 in
noncomputable def kc10 (B r l0 l1 lB a0 : ℝ) : ℝ :=
  (((-36) * (B) ^ 2 * (l1) ^ 2) + (18 * B * (l0) ^ 2) + (18 * B * (l1) ^ 2) + (18 * (B) ^ 3 * (l1) ^ 2) + (18 * (B) ^ 3 * (lB) ^ 2) + ((-54) * r * (B) ^ 3 * (lB) ^ 2) + ((-36) * B * l0 * l1) + ((-36) * l0 * lB * (B) ^ 2) + ((-36) * l1 * lB * (B) ^ 3) + ((-36) * r * (B) ^ 2 * (l0) ^ 2) + ((-18) * B * r * (l0) ^ 2) + ((-18) * B * r * (l1) ^ 2) + ((-18) * r * (B) ^ 3 * (l1) ^ 2) + ((-18) * (B) ^ 3 * (l0) ^ 2 * (r) ^ 3) + ((-18) * (B) ^ 3 * (lB) ^ 2 * (r) ^ 3) + (18 * (B) ^ 3 * (l0) ^ 2 * (r) ^ 2) + (36 * l0 * l1 * (B) ^ 2) + (36 * l1 * lB * (B) ^ 2) + (36 * r * (B) ^ 2 * (l1) ^ 2) + (36 * (B) ^ 2 * (l0) ^ 2 * (r) ^ 2) + (54 * (B) ^ 3 * (lB) ^ 2 * (r) ^ 2) + ((-72) * l0 * lB * (B) ^ 3 * (r) ^ 2) + ((-72) * l1 * lB * r * (B) ^ 2) + ((-36) * l0 * l1 * r * (B) ^ 3) + ((-36) * l0 * l1 * (B) ^ 2 * (r) ^ 2) + ((-36) * l0 * lB * (B) ^ 2 * (r) ^ 2) + ((-36) * l1 * lB * (B) ^ 3 * (r) ^ 2) + (36 * B * l0 * l1 * r) + (36 * l0 * l1 * (B) ^ 3 * (r) ^ 2) + (36 * l0 * lB * r * (B) ^ 3) + (36 * l0 * lB * (B) ^ 3 * (r) ^ 3) + (36 * l1 * lB * (B) ^ 2 * (r) ^ 2) + (72 * l0 * lB * r * (B) ^ 2) + (72 * l1 * lB * r * (B) ^ 3))

set_option maxHeartbeats 2000000 in
noncomputable def kc11 (B r l0 l1 lB a0 : ℝ) : ℝ :=
  (((-18) * l0 * (B) ^ 2) + ((-18) * l1 * (B) ^ 3) + (18 * l1 * (B) ^ 2) + (18 * lB * (B) ^ 3) + ((-54) * lB * r * (B) ^ 3) + ((-36) * l0 * (B) ^ 3 * (r) ^ 2) + ((-36) * l1 * r * (B) ^ 2) + ((-18) * l0 * (B) ^ 2 * (r) ^ 2) + ((-18) * l1 * (B) ^ 3 * (r) ^ 2) + ((-18) * lB * (B) ^ 3 * (r) ^ 3) + (18 * l0 * r * (B) ^ 3) + (18 * l0 * (B) ^ 3 * (r) ^ 3) + (18 * l1 * (B) ^ 2 * (r) ^ 2) + (36 * l0 * r * (B) ^ 2) + (36 * l1 * r * (B) ^ 3) + (54 * lB * (B) ^ 3 * (r) ^ 2))

set_option maxHeartbeats 2000000 in
noncomputable def kc12 (B r l0 l1 lB a0 : ℝ) : ℝ :=
  ((6 * (B) ^ 3) + ((-18) * r * (B) ^ 3) + ((-6) * (B) ^ 3 * (r) ^ 3) + (18 * (B) ^ 3 * (r) ^ 2))

set_option maxHeartbeats 2000000 in
noncomputable def kc13 (B r l0 l1 lB a0 : ℝ) : ℝ :=
  (((-72) * l1 * (B) ^ 2) + ((-36) * B * l0) + ((-36) * lB * (B) ^ 3) + (36 * B * l1) + (36 * l0 * (B) ^ 2) + (36 * l1 * (B) ^ 3) + (36 * lB * (B) ^ 2) + ((-72) * lB * r * (B) ^ 2) + ((-36) * B * l1 * r) + ((-36) * l0 * r * (B) ^ 3) + ((-36) * l0 * (B) ^ 2 * (r) ^ 2) + ((-36) * l1 * r * (B) ^ 3) + ((-36) * lB * (B) ^ 3 * (r) ^ 2) + (36 * B * l0 * r) + (36 * l0 * (B) ^ 3 * (r) ^ 2) + (36 * lB * (B) ^ 2 * (r) ^ 2) + (72 * l1 * r * (B) ^ 2) + (72 * lB * r * (B) ^ 3))

set_option maxHeartbeats 2000000 in
noncomputable def kc14 (B r l0 l1 lB a0 : ℝ) : ℝ :=
  (((-18) * (B) ^ 3) + (18 * (B) ^ 2) + ((-36) * r * (B) ^ 2) + ((-18) * (B) ^ 3 * (r) ^ 2) + (18 * (B) ^ 2 * (r) ^ 2) + (36 * r * (B) ^ 3))

set_option maxHeartbeats 2000000 in
noncomputable def kc15 (B r l0 l1 lB a0 : ℝ) : ℝ :=
  (((-36) * (B) ^ 2) + (18 * B) + (18 * (B) ^ 3) + ((-18) * B * r) + ((-18) * r * (B) ^ 3) + (36 * r * (B) ^ 2))

set_option maxHeartbeats 2000000 in
noncomputable def kc16 (B r l0 l1 lB a0 : ℝ) : ℝ :=
  (((-36) * B * l1) + ((-36) * lB * (B) ^ 2) + (36 * B * l0) + (36 * l1 * (B) ^ 2) + ((-72) * l0 * r * (B) ^ 2) + ((-72) * lB * (B) ^ 3 * (r) ^ 2) + ((-36) * B * l0 * r) + ((-36) * l0 * (B) ^ 3 * (r) ^ 3) + ((-36) * l1 * r * (B) ^ 3) + ((-36) * l1 * (B) ^ 2 * (r) ^ 2) + ((-36) * lB * (B) ^ 2 * (r) ^ 2) + (36 * B * l1 * r) + (36 * l0 * (B) ^ 3 * (r) ^ 2) + (36 * l1 * (B) ^ 3 * (r) ^ 2) + (36 * lB * r * (B) ^ 3) + (36 * lB * (B) ^ 3 * (r) ^ 3) + (72 * l0 * (B) ^ 2 * (r) ^ 2) + (72 * lB * r * (B) ^ 2))

set_option maxHeartbeats 2000000 in
noncomputable def kc17 (B r l0 l1 lB a0 : ℝ) : ℝ :=
  (((-18) * (B) ^ 2) + ((-36) * (B) ^ 3 * (r) ^ 2) + ((-18) * (B) ^ 2 * (r) ^ 2) + (18 * r * (B) ^ 3) + (18 * (B) ^ 3 * (r) ^ 3) + (36 * r * (B) ^ 2))

set_option maxHeartbeats 2000000 in
noncomputable def kc18 (B r l0 l1 lB a0 : ℝ) : ℝ :=
  (((-36) * B) + (36 * (B) ^ 2) + ((-36) * r * (B) ^ 3) + ((-36) * (B) ^ 2 * (r) ^ 2) + (36 * B * r) + (36 * (B) ^ 3 * (r) ^ 2))

set_option maxHeartbeats 2000000 in
noncomputable def kc19 (B r l0 l1 lB a0 : ℝ) : ℝ :=
  ((18 * B) + ((-36) * r * (B) ^ 2) + ((-18) * B * r) + ((-18) * (B) ^ 3 * (r) ^ 3) + (18 * (B) ^ 3 * (r) ^ 2) + (36 * (B) ^ 2 * (r) ^ 2))

set_option maxHeartbeats 2000000 in
noncomputable def kc20 (B r l0 l1 lB a0 : ℝ) : ℝ :=
  (((-6) * (l1) ^ 3) + (3 * (l0) ^ 3) + ((-24) * (B) ^ 2 * (l1) ^ 3) + ((-12) * l1 * (l0) ^ 2) + ((-9) * (B) ^ 3 * (lB) ^ 3) + ((-3) * a0 * (l1) ^ 3) + (3 * a0 * (l0) ^ 3) + (3 * (B) ^ 2 * (lB) ^ 3) + (9 * (B) ^ 3 * (l1) ^ 3) + (15 * l0 * (l1) ^ 2) + (21 * B * (l1) ^ 3) + ((-33) * B * l0 * (l1) ^ 2) + ((-30) * l1 * (B) ^ 2 * (lB) ^ 2) + ((-27) * lB * (B) ^ 3 * (l1) ^ 2) + ((-24) * B * lB * (l1) ^ 2) + ((-12) * B * lB * (l0) ^ 2) + ((-9) * a0 * l1 * (l0) ^ 2) + ((-9) * a0 * (B) ^ 2 * (l1) ^ 3) + ((-9) * (B) ^ 3 * (lB) ^ 3 * (r) ^ 2) + ((-6) * B * r * (l0) ^ 3) + ((-6) * r * (B) ^ 2 * (lB) ^ 3) + ((-3) * B * l0 * (lB) ^ 2) + ((-3) * a0 * (B) ^ 3 * (lB) ^ 3) + (3 * B * l1 * (lB) ^ 2) + (3 * a0 * (B) ^ 3 * (l1) ^ 3) + (3 * (B) ^ 2 * (l0) ^ 3 * (r) ^ 2) + (3 * (B) ^ 2 * (lB) ^ 3 * (r) ^ 2) + (9 * B * a0 * (l1) ^ 3) + (9 * a0 * l0 * (l1) ^ 2) + (12 * B * l1 * (l0) ^ 2) + (18 * l0 * (B) ^ 2 * (l1) ^ 2) + (18 * l0 * (B) ^ 2 * (lB) ^ 2) + (18 * r * (B) ^ 3 * (lB) ^ 3) + (27 * l1 * (B) ^ 3 * (lB) ^ 2) + (51 * lB * (B) ^ 2 * (l1) ^ 2) + ((-36) * l0 * l1 * lB * (B) ^ 2) + ((-36) * l1 * r * (B) ^ 3 * (lB) ^ 2) + ((-33) * lB * r * (B) ^ 2 * (l1) ^ 2) + ((-21) * l1 * r * (B) ^ 2 * (l0) ^ 2) + ((-18) * B * a0 * l0 * (l1) ^ 2) + ((-18) * l0 * r * (B) ^ 3 * (l1) ^ 2) + ((-18) * l0 * r * (B) ^ 3 * (lB) ^ 2) + ((-15) * B * l0 * r * (l1) ^ 2) + ((-15) * l0 * r * (B) ^ 2 * (lB) ^ 2) + ((-9) * B * a0 * lB * (l0) ^ 2) + ((-9) * B * a0 * lB * (l1) ^ 2) + ((-9) * B * a0 * r * (l0) ^ 3) + ((-9) * a0 * l1 * (B) ^ 2 * (lB) ^ 2) + ((-9) * a0 * lB * (B) ^ 3 * (l1) ^ 2) + ((-9) * a0 * (B) ^ 3 * (lB) ^ 3 * (r) ^ 2) + ((-9) * l1 * (B) ^ 2 * (l0) ^ 2 * (r) ^ 2) + ((-9) * l1 * (B) ^ 2 * (lB) ^ 2 * (r) ^ 2) + ((-9) * lB * (B) ^ 3 * (l0) ^ 2 * (r) ^ 2) + ((-3) * B * l1 * r * (lB) ^ 2) + ((-3) * a0 * (B) ^ 3 * (l0) ^ 3 * (r) ^ 3) + ((-3) * l0 * (B) ^ 2 * (lB) ^ 2 * (r) ^ 2) + ((-3) * lB * (B) ^ 2 * (l0) ^ 2 * (r) ^ 2) + (3 * B * l0 * r * (lB) ^ 2) + (3 * B * lB * r * (l0) ^ 2) + (3 * a0 * (B) ^ 3 * (lB) ^ 3 * (r) ^ 3) + (9 * B * a0 * l1 * (l0) ^ 2) + (9 * a0 * l0 * (B) ^ 2 * (l1) ^ 2) + (9 * a0 * l0 * (B) ^ 2 * (lB) ^ 2) + (9 * a0 * l1 * (B) ^ 3 * (lB) ^ 2) + (9 * a0 * r * (B) ^ 3 * (lB) ^ 3) + (9 * a0 * (B) ^ 2 * (l0) ^ 3 * (r) ^ 2) + (9 * l1 * (B) ^ 3 * (l0) ^ 2 * (r) ^ 2) + (9 * l1 * (B) ^ 3 * (lB) ^ 2 * (r) ^ 2) + (15 * B * lB * r * (l1) ^ 2) + (18 * a0 * lB * (B) ^ 2 * (l1) ^ 2) + (18 * l0 * (B) ^ 3 * (lB) ^ 2 * (r) ^ 2) + (18 * lB * r * (B) ^ 3 * (l1) ^ 2) + (21 * B * l1 * r * (l0) ^ 2) + (21 * lB * r * (B) ^ 2 * (l0) ^ 2) + (33 * l0 * r * (B) ^ 2 * (l1) ^ 2) + (36 * B * l0 * l1 * lB) + (39 * l1 * r * (B) ^ 2 * (lB) ^ 2) + ((-18) * B * l0 * l1 * lB * r) + ((-18) * a0 * l0 * l1 * lB * (B) ^ 2) + ((-18) * a0 * l0 * r * (B) ^ 2 * (lB) ^ 2) + ((-18) * a0 * l1 * r * (B) ^ 2 * (l0) ^ 2) + ((-18) * a0 * l1 * r * (B) ^ 3 * (lB) ^ 2) + ((-18) * a0 * lB * r * (B) ^ 2 * (l1) ^ 2) + ((-18) * a0 * lB * (B) ^ 2 * (l0) ^ 2 * (r) ^ 2) + ((-18) * l0 * l1 * lB * r * (B) ^ 2) + ((-18) * l0 * l1 * lB * (B) ^ 3 * (r) ^ 2) + ((-9) * B * a0 * l0 * r * (l1) ^ 2) + ((-9) * a0 * l0 * r * (B) ^ 3 * (l1) ^ 2) + ((-9) * a0 * l0 * r * (B) ^ 3 * (lB) ^ 2) + ((-9) * a0 * l0 * (B) ^ 3 * (lB) ^ 2 * (r) ^ 3) + ((-9) * a0 * l1 * (B) ^ 2 * (l0) ^ 2 * (r) ^ 2) + ((-9) * a0 * l1 * (B) ^ 2 * (lB) ^ 2 * (r) ^ 2) + ((-9) * a0 * lB * (B) ^ 3 * (l0) ^ 2 * (r) ^ 2) + (9 * B * a0 * lB * r * (l0) ^ 2) + (9 * B * a0 * lB * r * (l1) ^ 2) + (9 * a0 * l0 * (B) ^ 2 * (lB) ^ 2 * (r) ^ 2) + (9 * a0 * l1 * (B) ^ 3 * (l0) ^ 2 * (r) ^ 2) + (9 * a0 * l1 * (B) ^ 3 * (lB) ^ 2 * (r) ^ 2) + (9 * a0 * lB * r * (B) ^ 3 * (l1) ^ 2) + (9 * a0 * lB * (B) ^ 3 * (l0) ^ 2 * (r) ^ 3) + (18 * B * a0 * l0 * l1 * lB) + (18 * B * a0 * l1 * r * (l0) ^ 2) + (18 * a0 * l0 * r * (B) ^ 2 * (l1) ^ 2) + (18 * a0 * l0 * (B) ^ 3 * (lB) ^ 2 * (r) ^ 2) + (18 * a0 * l1 * r * (B) ^ 2 * (lB) ^ 2) + (18 * a0 * lB * r * (B) ^ 2 * (l0) ^ 2) + (18 * l0 * l1 * lB * (B) ^ 2 * (r) ^ 2) + (36 * l0 * l1 * lB * r * (B) ^ 3) + ((-18) * B * a0 * l0 * l1 * lB * r) + ((-18) * a0 * l0 * l1 * lB * (B) ^ 3 * (r) ^ 2) + (18 * a0 * l0 * l1 * lB * r * (B) ^ 3) + (18 * a0 * l0 * l1 * lB * (B) ^ 2 * (r) ^ 2))

set_option maxHeartbeats 2000000 in
noncomputable def kc21 (B r l0 l1 lB a0 : ℝ) : ℝ :=
  (((-18) * (l1) ^ 2) + ((-12) * (l0) ^ 2) + ((-72) * (B) ^ 2 * (l1) ^ 2) + ((-30) * (B) ^ 2 * (lB) ^ 2) + ((-9) * a0 * (l0) ^ 2) + ((-9) * a0 * (l1) ^ 2) + (3 * B * (lB) ^ 2) + (12 * B * (l0) ^ 2) + (27 * (B) ^ 3 * (l1) ^ 2) + (27 * (B) ^ 3 * (lB) ^ 2) + (30 * l0 * l1) + (63 * B * (l1) ^ 2) + ((-66) * B * l0 * l1) + ((-54) * l1 * lB * (B) ^ 3) + ((-48) * B * l1 * lB) + ((-36) * l0 * lB * (B) ^ 2) + ((-36) * r * (B) ^ 3 * (lB) ^ 2) + ((-27) * a0 * (B) ^ 2 * (l1) ^ 2) + ((-21) * r * (B) ^ 2 * (l0) ^ 2) + ((-9) * a0 * (B) ^ 2 * (lB) ^ 2) + ((-9) * (B) ^ 2 * (l0) ^ 2 * (r) ^ 2) + ((-9) * (B) ^ 2 * (lB) ^ 2 * (r) ^ 2) + ((-3) * B * r * (lB) ^ 2) + (9 * B * a0 * (l0) ^ 2) + (9 * a0 * (B) ^ 3 * (l1) ^ 2) + (9 * a0 * (B) ^ 3 * (lB) ^ 2) + (9 * (B) ^ 3 * (l0) ^ 2 * (r) ^ 2) + (9 * (B) ^ 3 * (lB) ^ 2 * (r) ^ 2) + (18 * a0 * l0 * l1) + (21 * B * r * (l0) ^ 2) + (27 * B * a0 * (l1) ^ 2) + (36 * B * l0 * lB) + (36 * l0 * l1 * (B) ^ 2) + (39 * r * (B) ^ 2 * (lB) ^ 2) + (102 * l1 * lB * (B) ^ 2) + ((-66) * l1 * lB * r * (B) ^ 2) + ((-36) * B * a0 * l0 * l1) + ((-36) * l0 * l1 * r * (B) ^ 3) + ((-30) * B * l0 * l1 * r) + ((-18) * B * a0 * l1 * lB) + ((-18) * B * l0 * lB * r) + ((-18) * a0 * l0 * lB * (B) ^ 2) + ((-18) * a0 * l1 * lB * (B) ^ 3) + ((-18) * a0 * r * (B) ^ 2 * (l0) ^ 2) + ((-18) * a0 * r * (B) ^ 3 * (lB) ^ 2) + ((-18) * l0 * lB * r * (B) ^ 2) + ((-18) * l0 * lB * (B) ^ 3 * (r) ^ 2) + ((-9) * a0 * (B) ^ 2 * (l0) ^ 2 * (r) ^ 2) + ((-9) * a0 * (B) ^ 2 * (lB) ^ 2 * (r) ^ 2) + (9 * a0 * (B) ^ 3 * (l0) ^ 2 * (r) ^ 2) + (9 * a0 * (B) ^ 3 * (lB) ^ 2 * (r) ^ 2) + (18 * B * a0 * l0 * lB) + (18 * B * a0 * r * (l0) ^ 2) + (18 * a0 * l0 * l1 * (B) ^ 2) + (18 * a0 * r * (B) ^ 2 * (lB) ^ 2) + (18 * l0 * lB * (B) ^ 2 * (r) ^ 2) + (30 * B * l1 * lB * r) + (36 * a0 * l1 * lB * (B) ^ 2) + (36 * l0 * lB * r * (B) ^ 3) + (36 * l1 * lB * r * (B) ^ 3) + (66 * l0 * l1 * r * (B) ^ 2) + ((-36) * a0 * l1 * lB * r * (B) ^ 2) + ((-18) * B * a0 * l0 * l1 * r) + ((-18) * B * a0 * l0 * lB * r) + ((-18) * a0 * l0 * l1 * r * (B) ^ 3) + ((-18) * a0 * l0 * lB * (B) ^ 3 * (r) ^ 2) + (18 * B * a0 * l1 * lB * r) + (18 * a0 * l0 * lB * r * (B) ^ 3) + (18 * a0 * l0 * lB * (B) ^ 2 * (r) ^ 2) + (18 * a0 * l1 * lB * r * (B) ^ 3) + (36 * a0 * l0 * l1 * r * (B) ^ 2))

set_option maxHeartbeats 2000000 in
noncomputable def kc22 (B r l0 l1 lB a0 : ℝ) : ℝ :=
  (((-18) * l1) + (15 * l0) + ((-72) * l1 * (B) ^ 2) + ((-33) * B * l0) + ((-27) * lB * (B) ^ 3) + ((-24) * B * lB) + ((-9) * a0 * l1) + (9 * a0 * l0) + (18 * l0 * (B) ^ 2) + (27 * l1 * (B) ^ 3) + (51 * lB * (B) ^ 2) + (63 * B * l1) + ((-33) * lB * r * (B) ^ 2) + ((-27) * a0 * l1 * (B) ^ 2) + ((-18) * B * a0 * l0) + ((-18) * l0 * r * (B) ^ 3) + ((-15) * B * l0 * r) + ((-9) * B * a0 * lB) + ((-9) * a0 * lB * (B) ^ 3) + (9 * a0 * l0 * (B) ^ 2) + (9 * a0 * l1 * (B) ^ 3) + (15 * B * lB * r) + (18 * a0 * lB * (B) ^ 2) + (18 * lB * r * (B) ^ 3) + (27 * B * a0 * l1) + (33 * l0 * r * (B) ^ 2) + ((-18) * a0 * lB * r * (B) ^ 2) + ((-9) * B * a0 * l0 * r) + ((-9) * a0 * l0 * r * (B) ^ 3) + (9 * B * a0 * lB * r) + (9 * a0 * lB * r * (B) ^ 3) + (18 * a0 * l0 * r * (B) ^ 2))

set_option maxHeartbeats 2000000 in
noncomputable def kc23 (B r l0 l1 lB a0 : ℝ) : ℝ :=
  ((-6) + ((-24) * (B) ^ 2) + ((-3) * a0) + (9 * (B) ^ 3) + (21 * B) + ((-9) * a0 * (B) ^ 2) + (3 * a0 * (B) ^ 3) + (9 * B * a0))

set_option maxHeartbeats 2000000 in
noncomputable def kc24 (B r l0 l1 lB a0 : ℝ) : ℝ :=
  ((9 * (l0) ^ 2) + (15 * (l1) ^ 2) + ((-33) * B * (l1) ^ 2) + ((-24) * l0 * l1) + ((-3) * B * (lB) ^ 2) + (9 * a0 * (l0) ^ 2) + (9 * a0 * (l1) ^ 2) + (18 * (B) ^ 2 * (l1) ^ 2) + (18 * (B) ^ 2 * (lB) ^ 2) + ((-36) * l1 * lB * (B) ^ 2) + ((-24) * B * l0 * lB) + ((-18) * B * a0 * (l1) ^ 2) + ((-18) * B * r * (l0) ^ 2) + ((-18) * a0 * l0 * l1) + ((-18) * r * (B) ^ 3 * (l1) ^ 2) + ((-18) * r * (B) ^ 3 * (lB) ^ 2) + ((-15) * B * r * (l1) ^ 2) + ((-15) * r * (B) ^ 2 * (lB) ^ 2) + ((-3) * (B) ^ 2 * (lB) ^ 2 * (r) ^ 2) + (3 * B * r * (lB) ^ 2) + (9 * a0 * (B) ^ 2 * (l1) ^ 2) + (9 * a0 * (B) ^ 2 * (lB) ^ 2) + (9 * (B) ^ 2 * (l0) ^ 2 * (r) ^ 2) + (18 * (B) ^ 3 * (lB) ^ 2 * (r) ^ 2) + (24 * B * l0 * l1) + (33 * r * (B) ^ 2 * (l1) ^ 2) + (36 * B * l1 * lB) + ((-42) * l0 * l1 * r * (B) ^ 2) + ((-27) * B * a0 * r * (l0) ^ 2) + ((-18) * B * a0 * l0 * lB) + ((-18) * B * l1 * lB * r) + ((-18) * a0 * l1 * lB * (B) ^ 2) + ((-18) * a0 * r * (B) ^ 2 * (lB) ^ 2) + ((-18) * l0 * l1 * (B) ^ 2 * (r) ^ 2) + ((-18) * l0 * lB * (B) ^ 3 * (r) ^ 2) + ((-18) * l1 * lB * r * (B) ^ 2) + ((-18) * l1 * lB * (B) ^ 3 * (r) ^ 2) + ((-9) * B * a0 * r * (l1) ^ 2) + ((-9) * a0 * r * (B) ^ 3 * (l1) ^ 2) + ((-9) * a0 * r * (B) ^ 3 * (lB) ^ 2) + ((-9) * a0 * (B) ^ 3 * (l0) ^ 2 * (r) ^ 3) + ((-9) * a0 * (B) ^ 3 * (lB) ^ 2 * (r) ^ 3) + ((-6) * l0 * lB * (B) ^ 2 * (r) ^ 2) + (6 * B * l0 * lB * r) + (9 * a0 * (B) ^ 2 * (lB) ^ 2 * (r) ^ 2) + (18 * B * a0 * l0 * l1) + (18 * B * a0 * l1 * lB) + (18 * a0 * r * (B) ^ 2 * (l1) ^ 2) + (18 * a0 * (B) ^ 3 * (lB) ^ 2 * (r) ^ 2) + (18 * l0 * l1 * (B) ^ 3 * (r) ^ 2) + (18 * l1 * lB * (B) ^ 2 * (r) ^ 2) + (27 * a0 * (B) ^ 2 * (l0) ^ 2 * (r) ^ 2) + (36 * l1 * lB * r * (B) ^ 3) + (42 * B * l0 * l1 * r) + (42 * l0 * lB * r * (B) ^ 2) + ((-36) * a0 * l0 * l1 * r * (B) ^ 2) + ((-36) * a0 * l0 * lB * (B) ^ 2 * (r) ^ 2) + ((-18) * B * a0 * l1 * lB * r) + ((-18) * a0 * l0 * l1 * (B) ^ 2 * (r) ^ 2) + ((-18) * a0 * l0 * lB * (B) ^ 3 * (r) ^ 2) + ((-18) * a0 * l1 * lB * (B) ^ 3 * (r) ^ 2) + (18 * B * a0 * l0 * lB * r) + (18 * a0 * l0 * l1 * (B) ^ 3 * (r) ^ 2) + (18 * a0 * l0 * lB * (B) ^ 3 * (r) ^ 3) + (18 * a0 * l1 * lB * r * (B) ^ 3) + (18 * a0 * l1 * lB * (B) ^ 2 * (r) ^ 2) + (36 * B * a0 * l0 * l1 * r) + (36 * a0 * l0 * lB * r * (B) ^ 2))

set_option maxHeartbeats 2000000 in
noncomputable def kc25 (B r l0 l1 lB a0 : ℝ) : ℝ :=
  (((-24) * l0) + (30 * l1) + ((-66) * B * l1) + ((-36) * lB * (B) ^ 2) + ((-18) * a0 * l0) + (18 * a0 * l1) + (24 * B * l0) + (36 * B * lB) + (36 * l1 * (B) ^ 2) + ((-42) * l0 * r * (B) ^ 2) + ((-36) * B * a0 * l1) + ((-36) * l1 * r * (B) ^ 3) + ((-30) * B * l1 * r) + ((-18) * B * lB * r) + ((-18) * a0 * lB * (B) ^ 2) + ((-18) * l0 * (B) ^ 2 * (r) ^ 2) + ((-18) * lB * r * (B) ^ 2) + ((-18) * lB * (B) ^ 3 * (r) ^ 2) + (18 * B * a0 * l0) + (18 * B * a0 * lB) + (18 * a0 * l1 * (B) ^ 2) + (18 * l0 * (B) ^ 3 * (r) ^ 2) + (18 * lB * (B) ^ 2 * (r) ^ 2) + (36 * lB * r * (B) ^ 3) + (42 * B * l0 * r) + (66 * l1 * r * (B) ^ 2) + ((-36) * a0 * l0 * r * (B) ^ 2) + ((-18) * B * a0 * l1 * r) + ((-18) * B * a0 * lB * r) + ((-18) * a0 * l0 * (B) ^ 2 * (r) ^ 2) + ((-18) * a0 * l1 * r * (B) ^ 3) + ((-18) * a0 * lB * (B) ^ 3 * (r) ^ 2) + (18 * a0 * l0 * (B) ^ 3 * (r) ^ 2) + (18 * a0 * lB * r * (B) ^ 3) + (18 * a0 * lB * (B) ^ 2 * (r) ^ 2) + (36 * B * a0 * l0 * r) + (36 * a0 * l1 * r * (B) ^ 2))

set_option maxHeartbeats 2000000 in
noncomputable def kc26 (B r l0 l1 lB a0 : ℝ) : ℝ :=
  (15 + ((-33) * B) + (9 * a0) + (18 * (B) ^ 2) + ((-18) * B * a0) + ((-18) * r * (B) ^ 3) + ((-15) * B * r) + (9 * a0 * (B) ^ 2) + (33 * r * (B) ^ 2) + ((-9) * B * a0 * r) + ((-9) * a0 * r * (B) ^ 3) + (18 * a0 * r * (B) ^ 2))

set_option maxHeartbeats 2000000 in
noncomputable def kc27 (B r l0 l1 lB a0 : ℝ) : ℝ :=
  (((-12) * l1) + (9 * l0) + ((-12) * B * lB) + ((-9) * a0 * l1) + (9 * a0 * l0) + (12 * B * l1) + ((-21) * l1 * r * (B) ^ 2) + ((-18) * B * l0 * r) + ((-9) * B * a0 * lB) + ((-9) * l1 * (B) ^ 2 * (r) ^ 2) + ((-9) * lB * (B) ^ 3 * (r) ^ 2) + ((-3) * lB * (B) ^ 2 * (r) ^ 2) + (3 * B * lB * r) + (9 * B * a0 * l1) + (9 * l0 * (B) ^ 2 * (r) ^ 2) + (9 * l1 * (B) ^ 3 * (r) ^ 2) + (21 * B * l1 * r) + (21 * lB * r * (B) ^ 2) + ((-27) * B * a0 * l0 * r) + ((-18) * a0 * l1 * r * (B) ^ 2) + ((-18) * a0 * lB * (B) ^ 2 * (r) ^ 2) + ((-9) * a0 * l0 * (B) ^ 3 * (r) ^ 3) + ((-9) * a0 * l1 * (B) ^ 2 * (r) ^ 2) + ((-9) * a0 * lB * (B) ^ 3 * (r) ^ 2) + (9 * B * a0 * lB * r) + (9 * a0 * l1 * (B) ^ 3 * (r) ^ 2) + (9 * a0 * lB * (B) ^ 3 * (r) ^ 3) + (18 * B * a0 * l1 * r) + (18 * a0 * lB * r * (B) ^ 2) + (27 * a0 * l0 * (B) ^ 2 * (r) ^ 2))

set_option maxHeartbeats 2000000 in
noncomputable def kc28 (B r l0 l1 lB a0 : ℝ) : ℝ :=
  ((-12) + ((-9) * a0) + (12 * B) + ((-21) * r * (B) ^ 2) + ((-9) * (B) ^ 2 * (r) ^ 2) + (9 * B * a0) + (9 * (B) ^ 3 * (r) ^ 2) + (21 * B * r) + ((-18) * a0 * r * (B) ^ 2) + ((-9) * a0 * (B) ^ 2 * (r) ^ 2) + (9 * a0 * (B) ^ 3 * (r) ^ 2) + (18 * B * a0 * r))

set_option maxHeartbeats 2000000 in
noncomputable def kc29 (B r l0 l1 lB a0 : ℝ) : ℝ :=
  (3 + (3 * a0) + ((-6) * B * r) + (3 * (B) ^ 2 * (r) ^ 2) + ((-9) * B * a0 * r) + ((-3) * a0 * (B) ^ 3 * (r) ^ 3) + (9 * a0 * (B) ^ 2 * (r) ^ 2))

set_option maxHeartbeats 2000000 in
noncomputable def kc30 (B r l0 l1 lB a0 : ℝ) : ℝ :=
  (((-6) * (l0) ^ 3) + (6 * (l1) ^ 3) + ((-18) * B * (l1) ^ 3) + ((-18) * l0 * (l1) ^ 2) + ((-6) * (B) ^ 3 * (l1) ^ 3) + (6 * (B) ^ 3 * (lB) ^ 3) + (18 * l1 * (l0) ^ 2) + (18 * (B) ^ 2 * (l1) ^ 3) + ((-36) * lB * (B) ^ 2 * (l1) ^ 2) + ((-18) * B * l1 * (l0) ^ 2) + ((-18) * l0 * (B) ^ 2 * (l1) ^ 2) + ((-18) * l0 * (B) ^ 2 * (lB) ^ 2) + ((-18) * l1 * (B) ^ 3 * (lB) ^ 2) + ((-18) * r * (B) ^ 3 * (lB) ^ 3) + ((-18) * (B) ^ 2 * (l0) ^ 3 * (r) ^ 2) + ((-6) * (B) ^ 3 * (lB) ^ 3 * (r) ^ 3) + (6 * (B) ^ 3 * (l0) ^ 3 * (r) ^ 3) + (18 * B * lB * (l0) ^ 2) + (18 * B * lB * (l1) ^ 2) + (18 * B * r * (l0) ^ 3) + (18 * l1 * (B) ^ 2 * (lB) ^ 2) + (18 * lB * (B) ^ 3 * (l1) ^ 2) + (18 * (B) ^ 3 * (lB) ^ 3 * (r) ^ 2) + (36 * B * l0 * (l1) ^ 2) + ((-36) * B * l0 * l1 * lB) + ((-36) * B * l1 * r * (l0) ^ 2) + ((-36) * l0 * r * (B) ^ 2 * (l1) ^ 2) + ((-36) * l0 * (B) ^ 3 * (lB) ^ 2 * (r) ^ 2) + ((-36) * l1 * r * (B) ^ 2 * (lB) ^ 2) + ((-36) * lB * r * (B) ^ 2 * (l0) ^ 2) + ((-18) * B * lB * r * (l0) ^ 2) + ((-18) * B * lB * r * (l1) ^ 2) + ((-18) * l0 * (B) ^ 2 * (lB) ^ 2 * (r) ^ 2) + ((-18) * l1 * (B) ^ 3 * (l0) ^ 2 * (r) ^ 2) + ((-18) * l1 * (B) ^ 3 * (lB) ^ 2 * (r) ^ 2) + ((-18) * lB * r * (B) ^ 3 * (l1) ^ 2) + ((-18) * lB * (B) ^ 3 * (l0) ^ 2 * (r) ^ 3) + (18 * B * l0 * r * (l1) ^ 2) + (18 * l0 * r * (B) ^ 3 * (l1) ^ 2) + (18 * l0 * r * (B) ^ 3 * (lB) ^ 2) + (18 * l0 * (B) ^ 3 * (lB) ^ 2 * (r) ^ 3) + (18 * l1 * (B) ^ 2 * (l0) ^ 2 * (r) ^ 2) + (18 * l1 * (B) ^ 2 * (lB) ^ 2 * (r) ^ 2) + (18 * lB * (B) ^ 3 * (l0) ^ 2 * (r) ^ 2) + (36 * l0 * l1 * lB * (B) ^ 2) + (36 * l0 * r * (B) ^ 2 * (lB) ^ 2) + (36 * l1 * r * (B) ^ 2 * (l0) ^ 2) + (36 * l1 * r * (B) ^ 3 * (lB) ^ 2) + (36 * lB * r * (B) ^ 2 * (l1) ^ 2) + (36 * lB * (B) ^ 2 * (l0) ^ 2 * (r) ^ 2) + ((-36) * l0 * l1 * lB * r * (B) ^ 3) + ((-36) * l0 * l1 * lB * (B) ^ 2 * (r) ^ 2) + (36 * B * l0 * l1 * lB * r) + (36 * l0 * l1 * lB * (B) ^ 3 * (r) ^ 2))

set_option maxHeartbeats 2000000 in
noncomputable def kc31 (B r l0 l1 lB a0 : ℝ) : ℝ :=
  ((18 * (l0) ^ 2) + (18 * (l1) ^ 2) + ((-54) * B * (l1) ^ 2) + ((-36) * l0 * l1) + ((-18) * B * (l0) ^ 2) + ((-18) * (B) ^ 3 * (l1) ^ 2) + ((-18) * (B) ^ 3 * (lB) ^ 2) + (18 * (B) ^ 2 * (lB) ^ 2) + (54 * (B) ^ 2 * (l1) ^ 2) + ((-72) * l1 * lB * (B) ^ 2) + ((-36) * B * l0 * lB) + ((-36) * B * r * (l0) ^ 2) + ((-36) * l0 * l1 * (B) ^ 2) + ((-36) * r * (B) ^ 2 * (lB) ^ 2) + ((-18) * (B) ^ 3 * (l0) ^ 2 * (r) ^ 2) + ((-18) * (B) ^ 3 * (lB) ^ 2 * (r) ^ 2) + (18 * (B) ^ 2 * (l0) ^ 2 * (r) ^ 2) + (18 * (B) ^ 2 * (lB) ^ 2 * (r) ^ 2) + (36 * B * l1 * lB) + (36 * l0 * lB * (B) ^ 2) + (36 * l1 * lB * (B) ^ 3) + (36 * r * (B) ^ 2 * (l0) ^ 2) + (36 * r * (B) ^ 3 * (lB) ^ 2) + (72 * B * l0 * l1) + ((-72) * l0 * l1 * r * (B) ^ 2) + ((-36) * B * l1 * lB * r) + ((-36) * l0 * lB * r * (B) ^ 3) + ((-36) * l0 * lB * (B) ^ 2 * (r) ^ 2) + ((-36) * l1 * lB * r * (B) ^ 3) + (36 * B * l0 * l1 * r) + (36 * B * l0 * lB * r) + (36 * l0 * l1 * r * (B) ^ 3) + (36 * l0 * lB * (B) ^ 3 * (r) ^ 2) + (72 * l1 * lB * r * (B) ^ 2))

set_option maxHeartbeats 2000000 in
noncomputable def kc32 (B r l0 l1 lB a0 : ℝ) : ℝ :=
  (((-18) * l0) + (18 * l1) + ((-54) * B * l1) + ((-36) * lB * (B) ^ 2) + ((-18) * l0 * (B) ^ 2) + ((-18) * l1 * (B) ^ 3) + (18 * B * lB) + (18 * lB * (B) ^ 3) + (36 * B * l0) + (54 * l1 * (B) ^ 2) + ((-36) * l0 * r * (B) ^ 2) + ((-18) * B * lB * r) + ((-18) * lB * r * (B) ^ 3) + (18 * B * l0 * r) + (18 * l0 * r * (B) ^ 3) + (36 * lB * r * (B) ^ 2))

set_option maxHeartbeats 2000000 in
noncomputable def kc33 (B r l0 l1 lB a0 : ℝ) : ℝ :=
  (6 + ((-18) * B) + ((-6) * (B) ^ 3) + (18 * (B) ^ 2))

set_option maxHeartbeats 2000000 in
noncomputable def kc34 (B r l0 l1 lB a0 : ℝ) : ℝ :=
  (((-18) * (l0) ^ 2) + ((-18) * (l1) ^ 2) + ((-18) * (B) ^ 2 * (l1) ^ 2) + ((-18) * (B) ^ 2 * (lB) ^ 2) + (36 * B * (l1) ^ 2) + (36 * l0 * l1) + ((-54) * (B) ^ 2 * (l0) ^ 2 * (r) ^ 2) + ((-36) * B * l0 * l1) + ((-36) * B * l1 * lB) + ((-36) * r * (B) ^ 2 * (l1) ^ 2) + ((-36) * (B) ^ 3 * (lB) ^ 2 * (r) ^ 2) + ((-18) * (B) ^ 2 * (lB) ^ 2 * (r) ^ 2) + (18 * B * r * (l1) ^ 2) + (18 * r * (B) ^ 3 * (l1) ^ 2) + (18 * r * (B) ^ 3 * (lB) ^ 2) + (18 * (B) ^ 3 * (l0) ^ 2 * (r) ^ 3) + (18 * (B) ^ 3 * (lB) ^ 2 * (r) ^ 3) + (36 * B * l0 * lB) + (36 * l1 * lB * (B) ^ 2) + (36 * r * (B) ^ 2 * (lB) ^ 2) + (54 * B * r * (l0) ^ 2) + ((-72) * B * l0 * l1 * r) + ((-72) * l0 * lB * r * (B) ^ 2) + ((-36) * B * l0 * lB * r) + ((-36) * l0 * l1 * (B) ^ 3 * (r) ^ 2) + ((-36) * l0 * lB * (B) ^ 3 * (r) ^ 3) + ((-36) * l1 * lB * r * (B) ^ 3) + ((-36) * l1 * lB * (B) ^ 2 * (r) ^ 2) + (36 * B * l1 * lB * r) + (36 * l0 * l1 * (B) ^ 2 * (r) ^ 2) + (36 * l0 * lB * (B) ^ 3 * (r) ^ 2) + (36 * l1 * lB * (B) ^ 3 * (r) ^ 2) + (72 * l0 * l1 * r * (B) ^ 2) + (72 * l0 * lB * (B) ^ 2 * (r) ^ 2))

set_option maxHeartbeats 2000000 in
noncomputable def kc35 (B r l0 l1 lB a0 : ℝ) : ℝ :=
  (((-36) * l1) + (36 * l0) + ((-36) * B * l0) + ((-36) * B * lB) + ((-36) * l1 * (B) ^ 2) + (36 * lB * (B) ^ 2) + (72 * B * l1) + ((-72) * B * l0 * r) + ((-72) * l1 * r * (B) ^ 2) + ((-36) * l0 * (B) ^ 3 * (r) ^ 2) + ((-36) * lB * r * (B) ^ 3) + ((-36) * lB * (B) ^ 2 * (r) ^ 2) + (36 * B * l1 * r) + (36 * B * lB * r) + (36 * l0 * (B) ^ 2 * (r) ^ 2) + (36 * l1 * r * (B) ^ 3) + (36 * lB * (B) ^ 3 * (r) ^ 2) + (72 * l0 * r * (B) ^ 2))

set_option maxHeartbeats 2000000 in
noncomputable def kc36 (B r l0 l1 lB a0 : ℝ) : ℝ :=
  ((-18) + ((-18) * (B) ^ 2) + (36 * B) + ((-36) * r * (B) ^ 2) + (18 * B * r) + (18 * r * (B) ^ 3))

set_option maxHeartbeats 2000000 in
noncomputable def kc37 (B r l0 l1 lB a0 : ℝ) : ℝ :=
  (((-18) * l0) + (18 * l1) + ((-18) * B * l1) + (18 * B * lB) + ((-54) * l0 * (B) ^ 2 * (r) ^ 2) + ((-36) * B * l1 * r) + ((-36) * lB * r * (B) ^ 2) + ((-18) * B * lB * r) + ((-18) * l1 * (B) ^ 3 * (r) ^ 2) + ((-18) * lB * (B) ^ 3 * (r) ^ 3) + (18 * l0 * (B) ^ 3 * (r) ^ 3) + (18 * l1 * (B) ^ 2 * (r) ^ 2) + (18 * lB * (B) ^ 3 * (r) ^ 2) + (36 * l1 * r * (B) ^ 2) + (36 * lB * (B) ^ 2 * (r) ^ 2) + (54 * B * l0 * r))

set_option maxHeartbeats 2000000 in
noncomputable def kc38 (B r l0 l1 lB a0 : ℝ) : ℝ :=
  (18 + ((-18) * B) + ((-36) * B * r) + ((-18) * (B) ^ 3 * (r) ^ 2) + (18 * (B) ^ 2 * (r) ^ 2) + (36 * r * (B) ^ 2))

set_option maxHeartbeats 2000000 in
noncomputable def kc39 (B r l0 l1 lB a0 : ℝ) : ℝ :=
  ((-6) + ((-18) * (B) ^ 2 * (r) ^ 2) + (6 * (B) ^ 3 * (r) ^ 3) + (18 * B * r))

noncomputable def kvK (B r l0 l1 lB a0 : ℝ) : Fin 40 → ℝ :=
  ![kc0 B r l0 l1 lB a0, kc1 B r l0 l1 lB a0, kc2 B r l0 l1 lB a0, kc3 B r l0 l1 lB a0, kc4 B r l0 l1 lB a0, kc5 B r l0 l1 lB a0, kc6 B r l0 l1 lB a0, kc7 B r l0 l1 lB a0, kc8 B r l0 l1 lB a0, kc9 B r l0 l1 lB a0, kc10 B r l0 l1 lB a0, kc11 B r l0 l1 lB a0, kc12 B r l0 l1 lB a0, kc13 B r l0 l1 lB a0, kc14 B r l0 l1 lB a0, kc15 B r l0 l1 lB a0, kc16 B r l0 l1 lB a0, kc17 B r l0 l1 lB a0, kc18 B r l0 l1 lB a0, kc19 B r l0 l1 lB a0, kc20 B r l0 l1 lB a0, kc21 B r l0 l1 lB a0, kc22 B r l0 l1 lB a0, kc23 B r l0 l1 lB a0, kc24 B r l0 l1 lB a0, kc25 B r l0 l1 lB a0, kc26 B r l0 l1 lB a0, kc27 B r l0 l1 lB a0, kc28 B r l0 l1 lB a0, kc29 B r l0 l1 lB a0, kc30 B r l0 l1 lB a0, kc31 B r l0 l1 lB a0, kc32 B r l0 l1 lB a0, kc33 B r l0 l1 lB a0, kc34 B r l0 l1 lB a0, kc35 B r l0 l1 lB a0, kc36 B r l0 l1 lB a0, kc37 B r l0 l1 lB a0, kc38 B r l0 l1 lB a0, kc39 B r l0 l1 lB a0]
def nvK : Fin 40 → ℕ := ![0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
def avK : Fin 40 → ℕ := ![0, 0, 0, 0, 0, 0, 1, 1, 1, 2, 0, 0, 0, 0, 0, 0, 1, 1, 1, 2, 0, 0, 0, 0, 1, 1, 1, 2, 2, 3, 0, 0, 0, 0, 1, 1, 1, 2, 2, 3]
def bvK : Fin 40 → ℕ := ![0, 0, 0, 1, 1, 2, 0, 0, 1, 0, 0, 0, 0, 1, 1, 2, 0, 0, 1, 0, 0, 1, 2, 3, 0, 1, 2, 0, 1, 0, 0, 1, 2, 3, 0, 1, 2, 0, 1, 0]
def cvK : Fin 40 → ℕ := ![1, 2, 3, 1, 2, 1, 1, 2, 1, 1, 1, 2, 3, 1, 2, 1, 1, 2, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]

set_option maxRecDepth 20000 in
set_option maxHeartbeats 2000000 in
lemma HH_expand (B r l0 l1 lB a0 t : ℝ) :
    S05.HH (a0 + t) (gF B r l0 l1 lB 1 t) (gS B r l0 l1 lB 1 t) (gS B r l0 l1 lB 2 t) =
      jcomb (kvK B r l0 l1 lB a0) nvK avK bvK cvK (1 - r) t := by
  simp only [S05.HH, S05.gg1, S05.kp1, S05.uu2, gF, gS, jcomb, jmon, Fin.sum_univ_succ,
    Fin.sum_univ_zero, kvK, kc0, kc1, kc2, kc3, kc4, kc5, kc6, kc7, kc8, kc9, kc10, kc11, kc12, kc13, kc14, kc15, kc16, kc17, kc18, kc19, kc20, kc21, kc22, kc23, kc24, kc25, kc26, kc27, kc28, kc29, kc30, kc31, kc32, kc33, kc34, kc35, kc36, kc37, kc38, kc39, nvK, avK, bvK, cvK, Matrix.cons_val_zero, Matrix.cons_val_succ, add_zero]
  ring

/-! ## 7. `E_e` as the `L`-part plus forty `J`'s -/

lemma ii_jcomb {y : ℝ} (hy : y < 1) {m : ℕ} (kv : Fin m → ℝ) (nv av bv cv : Fin m → ℕ) :
    IntervalIntegrable (jcomb kv nv av bv cv y) volume 0 1 := by
  unfold jcomb
  have : IntervalIntegrable (∑ i, fun t : ℝ => kv i * jmon (nv i) (av i) (bv i) (cv i) y t)
      volume 0 1 :=
    IntervalIntegrable.sum _ fun i _ => (ii_jmon hy _ _ _ _).const_mul _
  convert this using 1
  funext t; simp [Finset.sum_apply]

lemma Ee_gen (H : S05.Hex) (hH : S05.IsHex H) (k : Fin 3) :
    S05.Ee H k = S05.Len H k * ((S05.q0 H k / S05.Lsq H k + 1) *
        S05.gg3 (eF H k 1 1) (eF H k 2 1) (eF H k 3 1) -
        S05.q0 H k / S05.Lsq H k * S05.gg3 (eF H k 1 0) (eF H k 2 0) (eF H k 3 0)) +
      (S05.V02 H k - S05.q0 H k ^ 2 / S05.Lsq H k) / S05.Len H k *
        ∑ i, kvK (S05.Bv H k) (S05.rr H k) (S05.l0 H k) (S05.l1 H k) (S05.lB H k)
          (S05.q0 H k / S05.Lsq H k) i * S05.Jint (nvK i) (avK i) (bvK i) (cvK i) (1 - S05.rr H k) := by
  obtain ⟨-, -, -, -, -, -, hr⟩ := edge_facts H hH k
  set a0 := S05.q0 H k / S05.Lsq H k
  set K := (S05.V02 H k - S05.q0 H k ^ 2 / S05.Lsq H k) / S05.Len H k
  set EL := fun t => S05.gg3 (eF H k 1 t) (eF H k 2 t) (eF H k 3 t) -
    (a0 + t) * S05.WW (eF H k 1 t) (eF H k 2 t) (eS H k 1 t) (eS H k 2 t) (eS H k 3 t) with hEL
  set HL := fun t => S05.HH (a0 + t) (eF H k 1 t) (eS H k 1 t) (eS H k 2 t) with hHL
  have hHLj : HL = jcomb (kvK (S05.Bv H k) (S05.rr H k) (S05.l0 H k) (S05.l1 H k) (S05.lB H k) a0)
      nvK avK bvK cvK (1 - S05.rr H k) := by
    funext t; exact HH_expand _ _ _ _ _ _ _
  have hy : 1 - S05.rr H k < 1 := by linarith
  have hLp : 0 < S05.Len H k := Real.sqrt_pos.2 (S05.Lsq_pos H hH k)
  have hG := S05.Gw_int H hH k (j := 3) le_rfl (p := 1) ⟨by norm_num, by norm_num⟩
  have hae : ∀ᵐ t ∂volume, t ∈ Ι (0:ℝ) 1 → S05.Gw H k 3 1 t = S05.Len H k * EL t + K * HL t :=
    S05.ae_io fun t ht => Gw3_gen H hH k ht
  have hsum : IntervalIntegrable (fun t => S05.Len H k * EL t + K * HL t) volume 0 1 :=
    hG.congr_ae ((ae_restrict_iff' measurableSet_uIoc).2 hae)
  have hHi0 : IntervalIntegrable HL volume 0 1 := by rw [hHLj]; exact ii_jcomb hy _ _ _ _ _
  have hHi : IntervalIntegrable (fun t => K * HL t) volume 0 1 := hHi0.const_mul K
  have hLi : IntervalIntegrable (fun t => S05.Len H k * EL t) volume 0 1 := by
    have := hsum.sub hHi
    refine this.congr_ae (Filter.Eventually.of_forall fun t => ?_)
    simp
  have hEi : IntervalIntegrable EL volume 0 1 := by
    have := hLi.const_mul (S05.Len H k)⁻¹
    refine this.congr_ae (Filter.Eventually.of_forall fun t => ?_)
    simp only
    field_simp
  have hL := Lpart_gen (S05.Bv H k) (S05.rr H k) (S05.l0 H k) (S05.l1 H k) (S05.lB H k) a0 hr hEi
  have hK : (∫ t in (0:ℝ)..1, HL t) = ∑ i, kvK (S05.Bv H k) (S05.rr H k) (S05.l0 H k) (S05.l1 H k)
      (S05.lB H k) a0 i * S05.Jint (nvK i) (avK i) (bvK i) (cvK i) (1 - S05.rr H k) := by
    rw [hHLj, int_jcomb hy]
  rw [S05.Ee_eq H hH, S05.Phi, intervalIntegral.integral_congr_ae hae,
    intervalIntegral.integral_add hLi hHi, intervalIntegral.integral_const_mul,
    intervalIntegral.integral_const_mul, hK]
  have hL' : (∫ t in (0:ℝ)..1, EL t) = (a0 + 1) * S05.gg3 (eF H k 1 1) (eF H k 2 1) (eF H k 3 1) -
      a0 * S05.gg3 (eF H k 1 0) (eF H k 2 0) (eF H k 3 0) := hL
  rw [hL']

/-! ## 8. the six unprinted blocks (from `scripts/g4_formula_table.txt`) -/

noncomputable def pc4_0 : MvPolynomial (Fin 8) ℤ :=
  (MvPolynomial.C (3 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (3 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (-3 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (-3 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (3 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (6 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (2 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (3 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (3 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (2 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (6 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (3 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (-9 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (3 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (3 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (9 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (3 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (-6 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (2 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (-3 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (3 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (2 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (-6 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (3 : ℕ) : MvPolynomial (Fin 8) ℤ)

noncomputable def rc4_0 (B r q V L2 l0 l1 lB : ℝ) : ℝ :=
  ((3 : ℝ) * B ^ (3 : ℕ) * r ^ (3 : ℕ) * q ^ (3 : ℕ) : ℝ) +
    ((-3 : ℝ) * B ^ (3 : ℕ) * r ^ (3 : ℕ) * q ^ (2 : ℕ) * L2 ^ (1 : ℕ) : ℝ) +
    ((-3 : ℝ) * B ^ (3 : ℕ) * r ^ (3 : ℕ) * q ^ (1 : ℕ) * V ^ (1 : ℕ) * L2 ^ (1 : ℕ) : ℝ) +
    ((6 : ℝ) * B ^ (3 : ℕ) * r ^ (3 : ℕ) * q ^ (1 : ℕ) * L2 ^ (2 : ℕ) : ℝ) +
    ((3 : ℝ) * B ^ (3 : ℕ) * r ^ (3 : ℕ) * V ^ (1 : ℕ) * L2 ^ (2 : ℕ) : ℝ) +
    ((6 : ℝ) * B ^ (3 : ℕ) * r ^ (3 : ℕ) * L2 ^ (3 : ℕ) : ℝ) +
    ((-9 : ℝ) * B ^ (2 : ℕ) * r ^ (2 : ℕ) * q ^ (3 : ℕ) : ℝ) +
    ((3 : ℝ) * B ^ (2 : ℕ) * r ^ (2 : ℕ) * q ^ (2 : ℕ) * L2 ^ (1 : ℕ) : ℝ) +
    ((9 : ℝ) * B ^ (2 : ℕ) * r ^ (2 : ℕ) * q ^ (1 : ℕ) * V ^ (1 : ℕ) * L2 ^ (1 : ℕ) : ℝ) +
    ((-6 : ℝ) * B ^ (2 : ℕ) * r ^ (2 : ℕ) * q ^ (1 : ℕ) * L2 ^ (2 : ℕ) : ℝ) +
    ((-3 : ℝ) * B ^ (2 : ℕ) * r ^ (2 : ℕ) * V ^ (1 : ℕ) * L2 ^ (2 : ℕ) : ℝ) +
    ((-6 : ℝ) * B ^ (2 : ℕ) * r ^ (2 : ℕ) * L2 ^ (3 : ℕ) : ℝ)

noncomputable def pc4_1 : MvPolynomial (Fin 8) ℤ :=
  (MvPolynomial.C (3 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (3 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (3 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (-3 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (3 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (-3 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (3 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (2 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (6 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (3 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (-1 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (-6 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (3 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (1 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (2 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (1 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (3 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (2 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (1 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (3 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (-3 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (3 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (-2 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ)

noncomputable def rc4_1 (B r q V L2 l0 l1 lB : ℝ) : ℝ :=
  ((3 : ℝ) * B ^ (2 : ℕ) * r ^ (1 : ℕ) * q ^ (3 : ℕ) : ℝ) +
    ((3 : ℝ) * B ^ (2 : ℕ) * r ^ (1 : ℕ) * q ^ (2 : ℕ) * L2 ^ (1 : ℕ) : ℝ) +
    ((-3 : ℝ) * B ^ (2 : ℕ) * r ^ (1 : ℕ) * q ^ (1 : ℕ) * V ^ (1 : ℕ) * L2 ^ (1 : ℕ) : ℝ) +
    ((-3 : ℝ) * B ^ (2 : ℕ) * r ^ (1 : ℕ) * V ^ (1 : ℕ) * L2 ^ (2 : ℕ) : ℝ) +
    ((6 : ℝ) * B ^ (1 : ℕ) * r ^ (1 : ℕ) * q ^ (3 : ℕ) : ℝ) +
    ((-1 : ℝ) * B ^ (1 : ℕ) * r ^ (1 : ℕ) * q ^ (2 : ℕ) * L2 ^ (1 : ℕ) : ℝ) +
    ((-6 : ℝ) * B ^ (1 : ℕ) * r ^ (1 : ℕ) * q ^ (1 : ℕ) * V ^ (1 : ℕ) * L2 ^ (1 : ℕ) : ℝ) +
    ((1 : ℝ) * B ^ (1 : ℕ) * r ^ (1 : ℕ) * q ^ (1 : ℕ) * L2 ^ (2 : ℕ) : ℝ) +
    ((1 : ℝ) * B ^ (1 : ℕ) * r ^ (1 : ℕ) * V ^ (1 : ℕ) * L2 ^ (2 : ℕ) : ℝ) +
    ((1 : ℝ) * B ^ (1 : ℕ) * r ^ (1 : ℕ) * L2 ^ (3 : ℕ) : ℝ) +
    ((-3 : ℝ) * B ^ (1 : ℕ) * q ^ (3 : ℕ) : ℝ) +
    ((-2 : ℝ) * B ^ (1 : ℕ) * q ^ (2 : ℕ) * L2 ^ (1 : ℕ) : ℝ)

noncomputable def pc4_2 : MvPolynomial (Fin 8) ℤ :=
  (MvPolynomial.C (3 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (3 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (2 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (3 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (2 : ℕ) : MvPolynomial (Fin 8) ℤ)

noncomputable def rc4_2 (B r q V L2 l0 l1 lB : ℝ) : ℝ :=
  ((3 : ℝ) * B ^ (1 : ℕ) * q ^ (1 : ℕ) * V ^ (1 : ℕ) * L2 ^ (1 : ℕ) : ℝ) +
    ((2 : ℝ) * B ^ (1 : ℕ) * V ^ (1 : ℕ) * L2 ^ (2 : ℕ) : ℝ)

/-- block `c4` (coefficient of `ell^3`). -/
noncomputable def pc4 : MvPolynomial (Fin 8) ℤ :=
  pc4_0 + pc4_1 + pc4_2

noncomputable def rc4 (B r q V L2 l0 l1 lB : ℝ) : ℝ :=
  rc4_0 B r q V L2 l0 l1 lB + rc4_1 B r q V L2 l0 l1 lB + rc4_2 B r q V L2 l0 l1 lB

noncomputable def pc5_0 : MvPolynomial (Fin 8) ℤ :=
  (MvPolynomial.C (-18 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (5 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (18 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (18 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (3 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (5 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (-18 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (3 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (18 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (5 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (18 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (6 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (-36 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (6 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (5 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (12 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (6 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (-18 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (-18 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (3 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (5 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (-18 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (3 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (6 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ)

noncomputable def rc5_0 (B r q V L2 l0 l1 lB : ℝ) : ℝ :=
  ((-18 : ℝ) * B ^ (2 : ℕ) * r ^ (2 : ℕ) * q ^ (3 : ℕ) * l0 ^ (1 : ℕ) : ℝ) +
    ((18 : ℝ) * B ^ (2 : ℕ) * r ^ (2 : ℕ) * q ^ (3 : ℕ) * lB ^ (1 : ℕ) : ℝ) +
    ((18 : ℝ) * B ^ (2 : ℕ) * r ^ (2 : ℕ) * q ^ (1 : ℕ) * V ^ (1 : ℕ) * L2 ^ (1 : ℕ) * l0 ^ (1 : ℕ) : ℝ) +
    ((-18 : ℝ) * B ^ (2 : ℕ) * r ^ (2 : ℕ) * q ^ (1 : ℕ) * V ^ (1 : ℕ) * L2 ^ (1 : ℕ) * lB ^ (1 : ℕ) : ℝ) +
    ((18 : ℝ) * B ^ (2 : ℕ) * r ^ (1 : ℕ) * q ^ (3 : ℕ) * l0 ^ (1 : ℕ) : ℝ) +
    ((18 : ℝ) * B ^ (2 : ℕ) * r ^ (1 : ℕ) * q ^ (3 : ℕ) * l1 ^ (1 : ℕ) : ℝ) +
    ((-36 : ℝ) * B ^ (2 : ℕ) * r ^ (1 : ℕ) * q ^ (3 : ℕ) * lB ^ (1 : ℕ) : ℝ) +
    ((6 : ℝ) * B ^ (2 : ℕ) * r ^ (1 : ℕ) * q ^ (2 : ℕ) * L2 ^ (1 : ℕ) * l0 ^ (1 : ℕ) : ℝ) +
    ((12 : ℝ) * B ^ (2 : ℕ) * r ^ (1 : ℕ) * q ^ (2 : ℕ) * L2 ^ (1 : ℕ) * l1 ^ (1 : ℕ) : ℝ) +
    ((-18 : ℝ) * B ^ (2 : ℕ) * r ^ (1 : ℕ) * q ^ (2 : ℕ) * L2 ^ (1 : ℕ) * lB ^ (1 : ℕ) : ℝ) +
    ((-18 : ℝ) * B ^ (2 : ℕ) * r ^ (1 : ℕ) * q ^ (1 : ℕ) * V ^ (1 : ℕ) * L2 ^ (1 : ℕ) * l0 ^ (1 : ℕ) : ℝ) +
    ((-18 : ℝ) * B ^ (2 : ℕ) * r ^ (1 : ℕ) * q ^ (1 : ℕ) * V ^ (1 : ℕ) * L2 ^ (1 : ℕ) * l1 ^ (1 : ℕ) : ℝ)

noncomputable def pc5_1 : MvPolynomial (Fin 8) ℤ :=
  (MvPolynomial.C (36 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (3 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (-6 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (3 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (5 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (-12 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (3 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (6 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (18 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (3 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (-18 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (6 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (18 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (-18 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (6 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (18 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (18 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (3 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (6 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (-18 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (3 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (18 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (3 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (6 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (-18 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (3 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ)

noncomputable def rc5_1 (B r q V L2 l0 l1 lB : ℝ) : ℝ :=
  ((36 : ℝ) * B ^ (2 : ℕ) * r ^ (1 : ℕ) * q ^ (1 : ℕ) * V ^ (1 : ℕ) * L2 ^ (1 : ℕ) * lB ^ (1 : ℕ) : ℝ) +
    ((-6 : ℝ) * B ^ (2 : ℕ) * r ^ (1 : ℕ) * V ^ (1 : ℕ) * L2 ^ (2 : ℕ) * l0 ^ (1 : ℕ) : ℝ) +
    ((-12 : ℝ) * B ^ (2 : ℕ) * r ^ (1 : ℕ) * V ^ (1 : ℕ) * L2 ^ (2 : ℕ) * l1 ^ (1 : ℕ) : ℝ) +
    ((18 : ℝ) * B ^ (2 : ℕ) * r ^ (1 : ℕ) * V ^ (1 : ℕ) * L2 ^ (2 : ℕ) * lB ^ (1 : ℕ) : ℝ) +
    ((-18 : ℝ) * B ^ (2 : ℕ) * q ^ (3 : ℕ) * l1 ^ (1 : ℕ) : ℝ) +
    ((18 : ℝ) * B ^ (2 : ℕ) * q ^ (3 : ℕ) * lB ^ (1 : ℕ) : ℝ) +
    ((-18 : ℝ) * B ^ (2 : ℕ) * q ^ (2 : ℕ) * L2 ^ (1 : ℕ) * l1 ^ (1 : ℕ) : ℝ) +
    ((18 : ℝ) * B ^ (2 : ℕ) * q ^ (2 : ℕ) * L2 ^ (1 : ℕ) * lB ^ (1 : ℕ) : ℝ) +
    ((18 : ℝ) * B ^ (2 : ℕ) * q ^ (1 : ℕ) * V ^ (1 : ℕ) * L2 ^ (1 : ℕ) * l1 ^ (1 : ℕ) : ℝ) +
    ((-18 : ℝ) * B ^ (2 : ℕ) * q ^ (1 : ℕ) * V ^ (1 : ℕ) * L2 ^ (1 : ℕ) * lB ^ (1 : ℕ) : ℝ) +
    ((18 : ℝ) * B ^ (2 : ℕ) * V ^ (1 : ℕ) * L2 ^ (2 : ℕ) * l1 ^ (1 : ℕ) : ℝ) +
    ((-18 : ℝ) * B ^ (2 : ℕ) * V ^ (1 : ℕ) * L2 ^ (2 : ℕ) * lB ^ (1 : ℕ) : ℝ)

noncomputable def pc5_2 : MvPolynomial (Fin 8) ℤ :=
  (MvPolynomial.C (18 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (5 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (-18 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (6 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (-12 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (6 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (6 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (-12 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (-18 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (3 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (5 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (18 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (3 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (6 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (12 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (3 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (6 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (-6 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (3 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (12 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (3 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (2 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (-18 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (5 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (18 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (6 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ)

noncomputable def rc5_2 (B r q V L2 l0 l1 lB : ℝ) : ℝ :=
  ((18 : ℝ) * B ^ (1 : ℕ) * r ^ (1 : ℕ) * q ^ (3 : ℕ) * l0 ^ (1 : ℕ) : ℝ) +
    ((-18 : ℝ) * B ^ (1 : ℕ) * r ^ (1 : ℕ) * q ^ (3 : ℕ) * l1 ^ (1 : ℕ) : ℝ) +
    ((-12 : ℝ) * B ^ (1 : ℕ) * r ^ (1 : ℕ) * q ^ (2 : ℕ) * L2 ^ (1 : ℕ) * l1 ^ (1 : ℕ) : ℝ) +
    ((6 : ℝ) * B ^ (1 : ℕ) * r ^ (1 : ℕ) * q ^ (2 : ℕ) * L2 ^ (1 : ℕ) * lB ^ (1 : ℕ) : ℝ) +
    ((-12 : ℝ) * B ^ (1 : ℕ) * r ^ (1 : ℕ) * q ^ (2 : ℕ) * L2 ^ (1 : ℕ) : ℝ) +
    ((-18 : ℝ) * B ^ (1 : ℕ) * r ^ (1 : ℕ) * q ^ (1 : ℕ) * V ^ (1 : ℕ) * L2 ^ (1 : ℕ) * l0 ^ (1 : ℕ) : ℝ) +
    ((18 : ℝ) * B ^ (1 : ℕ) * r ^ (1 : ℕ) * q ^ (1 : ℕ) * V ^ (1 : ℕ) * L2 ^ (1 : ℕ) * l1 ^ (1 : ℕ) : ℝ) +
    ((12 : ℝ) * B ^ (1 : ℕ) * r ^ (1 : ℕ) * V ^ (1 : ℕ) * L2 ^ (2 : ℕ) * l1 ^ (1 : ℕ) : ℝ) +
    ((-6 : ℝ) * B ^ (1 : ℕ) * r ^ (1 : ℕ) * V ^ (1 : ℕ) * L2 ^ (2 : ℕ) * lB ^ (1 : ℕ) : ℝ) +
    ((12 : ℝ) * B ^ (1 : ℕ) * r ^ (1 : ℕ) * V ^ (1 : ℕ) * L2 ^ (2 : ℕ) : ℝ) +
    ((-18 : ℝ) * B ^ (1 : ℕ) * q ^ (3 : ℕ) * l0 ^ (1 : ℕ) : ℝ) +
    ((18 : ℝ) * B ^ (1 : ℕ) * q ^ (3 : ℕ) * l1 ^ (1 : ℕ) : ℝ)

noncomputable def pc5_3 : MvPolynomial (Fin 8) ℤ :=
  (MvPolynomial.C (-6 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (5 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (18 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (6 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (-6 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (12 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (18 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (3 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (5 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (-18 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (3 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (6 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (6 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (3 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (5 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (-18 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (3 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (6 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (6 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (3 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (-12 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (3 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (2 : ℕ) : MvPolynomial (Fin 8) ℤ)

noncomputable def rc5_3 (B r q V L2 l0 l1 lB : ℝ) : ℝ :=
  ((-6 : ℝ) * B ^ (1 : ℕ) * q ^ (2 : ℕ) * L2 ^ (1 : ℕ) * l0 ^ (1 : ℕ) : ℝ) +
    ((18 : ℝ) * B ^ (1 : ℕ) * q ^ (2 : ℕ) * L2 ^ (1 : ℕ) * l1 ^ (1 : ℕ) : ℝ) +
    ((-6 : ℝ) * B ^ (1 : ℕ) * q ^ (2 : ℕ) * L2 ^ (1 : ℕ) * lB ^ (1 : ℕ) : ℝ) +
    ((12 : ℝ) * B ^ (1 : ℕ) * q ^ (2 : ℕ) * L2 ^ (1 : ℕ) : ℝ) +
    ((18 : ℝ) * B ^ (1 : ℕ) * q ^ (1 : ℕ) * V ^ (1 : ℕ) * L2 ^ (1 : ℕ) * l0 ^ (1 : ℕ) : ℝ) +
    ((-18 : ℝ) * B ^ (1 : ℕ) * q ^ (1 : ℕ) * V ^ (1 : ℕ) * L2 ^ (1 : ℕ) * l1 ^ (1 : ℕ) : ℝ) +
    ((6 : ℝ) * B ^ (1 : ℕ) * V ^ (1 : ℕ) * L2 ^ (2 : ℕ) * l0 ^ (1 : ℕ) : ℝ) +
    ((-18 : ℝ) * B ^ (1 : ℕ) * V ^ (1 : ℕ) * L2 ^ (2 : ℕ) * l1 ^ (1 : ℕ) : ℝ) +
    ((6 : ℝ) * B ^ (1 : ℕ) * V ^ (1 : ℕ) * L2 ^ (2 : ℕ) * lB ^ (1 : ℕ) : ℝ) +
    ((-12 : ℝ) * B ^ (1 : ℕ) * V ^ (1 : ℕ) * L2 ^ (2 : ℕ) : ℝ)

/-- block `c5` (coefficient of `Li2(x)`). -/
noncomputable def pc5 : MvPolynomial (Fin 8) ℤ :=
  pc5_0 + pc5_1 + pc5_2 + pc5_3

noncomputable def rc5 (B r q V L2 l0 l1 lB : ℝ) : ℝ :=
  rc5_0 B r q V L2 l0 l1 lB + rc5_1 B r q V L2 l0 l1 lB + rc5_2 B r q V L2 l0 l1 lB + rc5_3 B r q V L2 l0 l1 lB

noncomputable def pc6_0 : MvPolynomial (Fin 8) ℤ :=
  (MvPolynomial.C (9 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (5 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (-9 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (-9 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (5 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (9 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (-9 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (3 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (5 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (9 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (3 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (18 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (5 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (-18 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (9 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (3 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (5 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (-9 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (3 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (18 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (5 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (-18 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ)

noncomputable def rc6_0 (B r q V L2 l0 l1 lB : ℝ) : ℝ :=
  ((9 : ℝ) * B ^ (3 : ℕ) * r ^ (3 : ℕ) * q ^ (3 : ℕ) * l0 ^ (1 : ℕ) : ℝ) +
    ((-9 : ℝ) * B ^ (3 : ℕ) * r ^ (3 : ℕ) * q ^ (3 : ℕ) * lB ^ (1 : ℕ) : ℝ) +
    ((-9 : ℝ) * B ^ (3 : ℕ) * r ^ (3 : ℕ) * q ^ (2 : ℕ) * L2 ^ (1 : ℕ) * l0 ^ (1 : ℕ) : ℝ) +
    ((9 : ℝ) * B ^ (3 : ℕ) * r ^ (3 : ℕ) * q ^ (2 : ℕ) * L2 ^ (1 : ℕ) * lB ^ (1 : ℕ) : ℝ) +
    ((-9 : ℝ) * B ^ (3 : ℕ) * r ^ (3 : ℕ) * q ^ (1 : ℕ) * V ^ (1 : ℕ) * L2 ^ (1 : ℕ) * l0 ^ (1 : ℕ) : ℝ) +
    ((9 : ℝ) * B ^ (3 : ℕ) * r ^ (3 : ℕ) * q ^ (1 : ℕ) * V ^ (1 : ℕ) * L2 ^ (1 : ℕ) * lB ^ (1 : ℕ) : ℝ) +
    ((18 : ℝ) * B ^ (3 : ℕ) * r ^ (3 : ℕ) * q ^ (1 : ℕ) * L2 ^ (2 : ℕ) * l0 ^ (1 : ℕ) : ℝ) +
    ((-18 : ℝ) * B ^ (3 : ℕ) * r ^ (3 : ℕ) * q ^ (1 : ℕ) * L2 ^ (2 : ℕ) * lB ^ (1 : ℕ) : ℝ) +
    ((9 : ℝ) * B ^ (3 : ℕ) * r ^ (3 : ℕ) * V ^ (1 : ℕ) * L2 ^ (2 : ℕ) * l0 ^ (1 : ℕ) : ℝ) +
    ((-9 : ℝ) * B ^ (3 : ℕ) * r ^ (3 : ℕ) * V ^ (1 : ℕ) * L2 ^ (2 : ℕ) * lB ^ (1 : ℕ) : ℝ) +
    ((18 : ℝ) * B ^ (3 : ℕ) * r ^ (3 : ℕ) * L2 ^ (3 : ℕ) * l0 ^ (1 : ℕ) : ℝ) +
    ((-18 : ℝ) * B ^ (3 : ℕ) * r ^ (3 : ℕ) * L2 ^ (3 : ℕ) * lB ^ (1 : ℕ) : ℝ)

noncomputable def pc6_1 : MvPolynomial (Fin 8) ℤ :=
  (MvPolynomial.C (-9 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (6 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (9 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (9 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (3 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (6 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (-9 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (3 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (-18 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (5 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (9 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (6 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (9 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (12 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (5 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (-9 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (6 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (18 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (3 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (5 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (-9 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (3 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (6 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ)

noncomputable def rc6_1 (B r q V L2 l0 l1 lB : ℝ) : ℝ :=
  ((-9 : ℝ) * B ^ (3 : ℕ) * r ^ (2 : ℕ) * q ^ (3 : ℕ) * l1 ^ (1 : ℕ) : ℝ) +
    ((9 : ℝ) * B ^ (3 : ℕ) * r ^ (2 : ℕ) * q ^ (3 : ℕ) * lB ^ (1 : ℕ) : ℝ) +
    ((9 : ℝ) * B ^ (3 : ℕ) * r ^ (2 : ℕ) * q ^ (1 : ℕ) * V ^ (1 : ℕ) * L2 ^ (1 : ℕ) * l1 ^ (1 : ℕ) : ℝ) +
    ((-9 : ℝ) * B ^ (3 : ℕ) * r ^ (2 : ℕ) * q ^ (1 : ℕ) * V ^ (1 : ℕ) * L2 ^ (1 : ℕ) * lB ^ (1 : ℕ) : ℝ) +
    ((-18 : ℝ) * B ^ (2 : ℕ) * r ^ (2 : ℕ) * q ^ (3 : ℕ) * l0 ^ (1 : ℕ) : ℝ) +
    ((9 : ℝ) * B ^ (2 : ℕ) * r ^ (2 : ℕ) * q ^ (3 : ℕ) * l1 ^ (1 : ℕ) : ℝ) +
    ((9 : ℝ) * B ^ (2 : ℕ) * r ^ (2 : ℕ) * q ^ (3 : ℕ) * lB ^ (1 : ℕ) : ℝ) +
    ((12 : ℝ) * B ^ (2 : ℕ) * r ^ (2 : ℕ) * q ^ (2 : ℕ) * L2 ^ (1 : ℕ) * l0 ^ (1 : ℕ) : ℝ) +
    ((-9 : ℝ) * B ^ (2 : ℕ) * r ^ (2 : ℕ) * q ^ (2 : ℕ) * L2 ^ (1 : ℕ) * lB ^ (1 : ℕ) : ℝ) +
    ((6 : ℝ) * B ^ (2 : ℕ) * r ^ (2 : ℕ) * q ^ (2 : ℕ) * L2 ^ (1 : ℕ) : ℝ) +
    ((18 : ℝ) * B ^ (2 : ℕ) * r ^ (2 : ℕ) * q ^ (1 : ℕ) * V ^ (1 : ℕ) * L2 ^ (1 : ℕ) * l0 ^ (1 : ℕ) : ℝ) +
    ((-9 : ℝ) * B ^ (2 : ℕ) * r ^ (2 : ℕ) * q ^ (1 : ℕ) * V ^ (1 : ℕ) * L2 ^ (1 : ℕ) * l1 ^ (1 : ℕ) : ℝ)

noncomputable def pc6_2 : MvPolynomial (Fin 8) ℤ :=
  (MvPolynomial.C (-9 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (3 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (-24 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (5 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (18 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (-9 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (2 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (-12 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (3 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (5 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (9 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (3 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (-6 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (3 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (2 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (-24 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (5 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (18 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (-9 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (3 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (9 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (6 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (-9 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ)

noncomputable def rc6_2 (B r q V L2 l0 l1 lB : ℝ) : ℝ :=
  ((-9 : ℝ) * B ^ (2 : ℕ) * r ^ (2 : ℕ) * q ^ (1 : ℕ) * V ^ (1 : ℕ) * L2 ^ (1 : ℕ) * lB ^ (1 : ℕ) : ℝ) +
    ((-24 : ℝ) * B ^ (2 : ℕ) * r ^ (2 : ℕ) * q ^ (1 : ℕ) * L2 ^ (2 : ℕ) * l0 ^ (1 : ℕ) : ℝ) +
    ((18 : ℝ) * B ^ (2 : ℕ) * r ^ (2 : ℕ) * q ^ (1 : ℕ) * L2 ^ (2 : ℕ) * lB ^ (1 : ℕ) : ℝ) +
    ((-9 : ℝ) * B ^ (2 : ℕ) * r ^ (2 : ℕ) * q ^ (1 : ℕ) * L2 ^ (2 : ℕ) : ℝ) +
    ((-12 : ℝ) * B ^ (2 : ℕ) * r ^ (2 : ℕ) * V ^ (1 : ℕ) * L2 ^ (2 : ℕ) * l0 ^ (1 : ℕ) : ℝ) +
    ((9 : ℝ) * B ^ (2 : ℕ) * r ^ (2 : ℕ) * V ^ (1 : ℕ) * L2 ^ (2 : ℕ) * lB ^ (1 : ℕ) : ℝ) +
    ((-6 : ℝ) * B ^ (2 : ℕ) * r ^ (2 : ℕ) * V ^ (1 : ℕ) * L2 ^ (2 : ℕ) : ℝ) +
    ((-24 : ℝ) * B ^ (2 : ℕ) * r ^ (2 : ℕ) * L2 ^ (3 : ℕ) * l0 ^ (1 : ℕ) : ℝ) +
    ((18 : ℝ) * B ^ (2 : ℕ) * r ^ (2 : ℕ) * L2 ^ (3 : ℕ) * lB ^ (1 : ℕ) : ℝ) +
    ((-9 : ℝ) * B ^ (2 : ℕ) * r ^ (2 : ℕ) * L2 ^ (3 : ℕ) : ℝ) +
    ((9 : ℝ) * B ^ (2 : ℕ) * r ^ (1 : ℕ) * q ^ (3 : ℕ) * l1 ^ (1 : ℕ) : ℝ) +
    ((-9 : ℝ) * B ^ (2 : ℕ) * r ^ (1 : ℕ) * q ^ (3 : ℕ) * lB ^ (1 : ℕ) : ℝ)

noncomputable def pc6_3 : MvPolynomial (Fin 8) ℤ :=
  (MvPolynomial.C (3 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (6 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (-3 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (-9 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (3 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (6 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (9 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (3 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (-3 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (3 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (6 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (3 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (3 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (9 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (5 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (-9 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (6 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (-3 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (5 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (-3 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (6 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (3 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (-6 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ)

noncomputable def rc6_3 (B r q V L2 l0 l1 lB : ℝ) : ℝ :=
  ((3 : ℝ) * B ^ (2 : ℕ) * r ^ (1 : ℕ) * q ^ (2 : ℕ) * L2 ^ (1 : ℕ) * l1 ^ (1 : ℕ) : ℝ) +
    ((-3 : ℝ) * B ^ (2 : ℕ) * r ^ (1 : ℕ) * q ^ (2 : ℕ) * L2 ^ (1 : ℕ) * lB ^ (1 : ℕ) : ℝ) +
    ((-9 : ℝ) * B ^ (2 : ℕ) * r ^ (1 : ℕ) * q ^ (1 : ℕ) * V ^ (1 : ℕ) * L2 ^ (1 : ℕ) * l1 ^ (1 : ℕ) : ℝ) +
    ((9 : ℝ) * B ^ (2 : ℕ) * r ^ (1 : ℕ) * q ^ (1 : ℕ) * V ^ (1 : ℕ) * L2 ^ (1 : ℕ) * lB ^ (1 : ℕ) : ℝ) +
    ((-3 : ℝ) * B ^ (2 : ℕ) * r ^ (1 : ℕ) * V ^ (1 : ℕ) * L2 ^ (2 : ℕ) * l1 ^ (1 : ℕ) : ℝ) +
    ((3 : ℝ) * B ^ (2 : ℕ) * r ^ (1 : ℕ) * V ^ (1 : ℕ) * L2 ^ (2 : ℕ) * lB ^ (1 : ℕ) : ℝ) +
    ((9 : ℝ) * B ^ (1 : ℕ) * r ^ (1 : ℕ) * q ^ (3 : ℕ) * l0 ^ (1 : ℕ) : ℝ) +
    ((-9 : ℝ) * B ^ (1 : ℕ) * r ^ (1 : ℕ) * q ^ (3 : ℕ) * l1 ^ (1 : ℕ) : ℝ) +
    ((-3 : ℝ) * B ^ (1 : ℕ) * r ^ (1 : ℕ) * q ^ (2 : ℕ) * L2 ^ (1 : ℕ) * l0 ^ (1 : ℕ) : ℝ) +
    ((-3 : ℝ) * B ^ (1 : ℕ) * r ^ (1 : ℕ) * q ^ (2 : ℕ) * L2 ^ (1 : ℕ) * l1 ^ (1 : ℕ) : ℝ) +
    ((3 : ℝ) * B ^ (1 : ℕ) * r ^ (1 : ℕ) * q ^ (2 : ℕ) * L2 ^ (1 : ℕ) * lB ^ (1 : ℕ) : ℝ) +
    ((-6 : ℝ) * B ^ (1 : ℕ) * r ^ (1 : ℕ) * q ^ (2 : ℕ) * L2 ^ (1 : ℕ) : ℝ)

noncomputable def pc6_4 : MvPolynomial (Fin 8) ℤ :=
  (MvPolynomial.C (-9 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (3 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (5 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (9 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (3 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (6 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (6 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (5 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (-3 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (3 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (2 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (3 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (3 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (5 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (3 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (3 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (6 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (-3 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (3 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (6 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (3 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (2 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (6 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (5 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (-3 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (3 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (3 : ℕ) : MvPolynomial (Fin 8) ℤ)

noncomputable def rc6_4 (B r q V L2 l0 l1 lB : ℝ) : ℝ :=
  ((-9 : ℝ) * B ^ (1 : ℕ) * r ^ (1 : ℕ) * q ^ (1 : ℕ) * V ^ (1 : ℕ) * L2 ^ (1 : ℕ) * l0 ^ (1 : ℕ) : ℝ) +
    ((9 : ℝ) * B ^ (1 : ℕ) * r ^ (1 : ℕ) * q ^ (1 : ℕ) * V ^ (1 : ℕ) * L2 ^ (1 : ℕ) * l1 ^ (1 : ℕ) : ℝ) +
    ((6 : ℝ) * B ^ (1 : ℕ) * r ^ (1 : ℕ) * q ^ (1 : ℕ) * L2 ^ (2 : ℕ) * l0 ^ (1 : ℕ) : ℝ) +
    ((-3 : ℝ) * B ^ (1 : ℕ) * r ^ (1 : ℕ) * q ^ (1 : ℕ) * L2 ^ (2 : ℕ) * lB ^ (1 : ℕ) : ℝ) +
    ((3 : ℝ) * B ^ (1 : ℕ) * r ^ (1 : ℕ) * q ^ (1 : ℕ) * L2 ^ (2 : ℕ) : ℝ) +
    ((3 : ℝ) * B ^ (1 : ℕ) * r ^ (1 : ℕ) * V ^ (1 : ℕ) * L2 ^ (2 : ℕ) * l0 ^ (1 : ℕ) : ℝ) +
    ((3 : ℝ) * B ^ (1 : ℕ) * r ^ (1 : ℕ) * V ^ (1 : ℕ) * L2 ^ (2 : ℕ) * l1 ^ (1 : ℕ) : ℝ) +
    ((-3 : ℝ) * B ^ (1 : ℕ) * r ^ (1 : ℕ) * V ^ (1 : ℕ) * L2 ^ (2 : ℕ) * lB ^ (1 : ℕ) : ℝ) +
    ((6 : ℝ) * B ^ (1 : ℕ) * r ^ (1 : ℕ) * V ^ (1 : ℕ) * L2 ^ (2 : ℕ) : ℝ) +
    ((6 : ℝ) * B ^ (1 : ℕ) * r ^ (1 : ℕ) * L2 ^ (3 : ℕ) * l0 ^ (1 : ℕ) : ℝ) +
    ((-3 : ℝ) * B ^ (1 : ℕ) * r ^ (1 : ℕ) * L2 ^ (3 : ℕ) * lB ^ (1 : ℕ) : ℝ) +
    ((3 : ℝ) * B ^ (1 : ℕ) * r ^ (1 : ℕ) * L2 ^ (3 : ℕ) : ℝ)

/-- block `c6` (coefficient of `ell^2`). -/
noncomputable def pc6 : MvPolynomial (Fin 8) ℤ :=
  pc6_0 + pc6_1 + pc6_2 + pc6_3 + pc6_4

noncomputable def rc6 (B r q V L2 l0 l1 lB : ℝ) : ℝ :=
  rc6_0 B r q V L2 l0 l1 lB + rc6_1 B r q V L2 l0 l1 lB + rc6_2 B r q V L2 l0 l1 lB + rc6_3 B r q V L2 l0 l1 lB + rc6_4 B r q V L2 l0 l1 lB

noncomputable def pc7_0 : MvPolynomial (Fin 8) ℤ :=
  (MvPolynomial.C (9 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (5 : Fin 8) ^ (2 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (-18 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (5 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (9 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (2 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (-9 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (5 : Fin 8) ^ (2 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (18 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (5 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (-9 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (2 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (-9 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (3 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (5 : Fin 8) ^ (2 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (18 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (3 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (5 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (-9 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (3 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (2 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (18 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (5 : Fin 8) ^ (2 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (-36 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (5 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (18 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (2 : ℕ) : MvPolynomial (Fin 8) ℤ)

noncomputable def rc7_0 (B r q V L2 l0 l1 lB : ℝ) : ℝ :=
  ((9 : ℝ) * B ^ (3 : ℕ) * r ^ (3 : ℕ) * q ^ (3 : ℕ) * l0 ^ (2 : ℕ) : ℝ) +
    ((-18 : ℝ) * B ^ (3 : ℕ) * r ^ (3 : ℕ) * q ^ (3 : ℕ) * l0 ^ (1 : ℕ) * lB ^ (1 : ℕ) : ℝ) +
    ((9 : ℝ) * B ^ (3 : ℕ) * r ^ (3 : ℕ) * q ^ (3 : ℕ) * lB ^ (2 : ℕ) : ℝ) +
    ((-9 : ℝ) * B ^ (3 : ℕ) * r ^ (3 : ℕ) * q ^ (2 : ℕ) * L2 ^ (1 : ℕ) * l0 ^ (2 : ℕ) : ℝ) +
    ((18 : ℝ) * B ^ (3 : ℕ) * r ^ (3 : ℕ) * q ^ (2 : ℕ) * L2 ^ (1 : ℕ) * l0 ^ (1 : ℕ) * lB ^ (1 : ℕ) : ℝ) +
    ((-9 : ℝ) * B ^ (3 : ℕ) * r ^ (3 : ℕ) * q ^ (2 : ℕ) * L2 ^ (1 : ℕ) * lB ^ (2 : ℕ) : ℝ) +
    ((-9 : ℝ) * B ^ (3 : ℕ) * r ^ (3 : ℕ) * q ^ (1 : ℕ) * V ^ (1 : ℕ) * L2 ^ (1 : ℕ) * l0 ^ (2 : ℕ) : ℝ) +
    ((18 : ℝ) * B ^ (3 : ℕ) * r ^ (3 : ℕ) * q ^ (1 : ℕ) * V ^ (1 : ℕ) * L2 ^ (1 : ℕ) * l0 ^ (1 : ℕ) * lB ^ (1 : ℕ) : ℝ) +
    ((-9 : ℝ) * B ^ (3 : ℕ) * r ^ (3 : ℕ) * q ^ (1 : ℕ) * V ^ (1 : ℕ) * L2 ^ (1 : ℕ) * lB ^ (2 : ℕ) : ℝ) +
    ((18 : ℝ) * B ^ (3 : ℕ) * r ^ (3 : ℕ) * q ^ (1 : ℕ) * L2 ^ (2 : ℕ) * l0 ^ (2 : ℕ) : ℝ) +
    ((-36 : ℝ) * B ^ (3 : ℕ) * r ^ (3 : ℕ) * q ^ (1 : ℕ) * L2 ^ (2 : ℕ) * l0 ^ (1 : ℕ) * lB ^ (1 : ℕ) : ℝ) +
    ((18 : ℝ) * B ^ (3 : ℕ) * r ^ (3 : ℕ) * q ^ (1 : ℕ) * L2 ^ (2 : ℕ) * lB ^ (2 : ℕ) : ℝ)

noncomputable def pc7_1 : MvPolynomial (Fin 8) ℤ :=
  (MvPolynomial.C (9 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (3 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (5 : Fin 8) ^ (2 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (-18 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (3 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (5 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (9 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (3 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (2 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (18 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (5 : Fin 8) ^ (2 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (-36 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (5 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (18 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (2 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (-18 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (5 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (6 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (18 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (5 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (18 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (6 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (-18 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (2 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (18 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (3 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (5 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (6 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (-18 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (3 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (5 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ)

noncomputable def rc7_1 (B r q V L2 l0 l1 lB : ℝ) : ℝ :=
  ((9 : ℝ) * B ^ (3 : ℕ) * r ^ (3 : ℕ) * V ^ (1 : ℕ) * L2 ^ (2 : ℕ) * l0 ^ (2 : ℕ) : ℝ) +
    ((-18 : ℝ) * B ^ (3 : ℕ) * r ^ (3 : ℕ) * V ^ (1 : ℕ) * L2 ^ (2 : ℕ) * l0 ^ (1 : ℕ) * lB ^ (1 : ℕ) : ℝ) +
    ((9 : ℝ) * B ^ (3 : ℕ) * r ^ (3 : ℕ) * V ^ (1 : ℕ) * L2 ^ (2 : ℕ) * lB ^ (2 : ℕ) : ℝ) +
    ((18 : ℝ) * B ^ (3 : ℕ) * r ^ (3 : ℕ) * L2 ^ (3 : ℕ) * l0 ^ (2 : ℕ) : ℝ) +
    ((-36 : ℝ) * B ^ (3 : ℕ) * r ^ (3 : ℕ) * L2 ^ (3 : ℕ) * l0 ^ (1 : ℕ) * lB ^ (1 : ℕ) : ℝ) +
    ((18 : ℝ) * B ^ (3 : ℕ) * r ^ (3 : ℕ) * L2 ^ (3 : ℕ) * lB ^ (2 : ℕ) : ℝ) +
    ((-18 : ℝ) * B ^ (3 : ℕ) * r ^ (2 : ℕ) * q ^ (3 : ℕ) * l0 ^ (1 : ℕ) * l1 ^ (1 : ℕ) : ℝ) +
    ((18 : ℝ) * B ^ (3 : ℕ) * r ^ (2 : ℕ) * q ^ (3 : ℕ) * l0 ^ (1 : ℕ) * lB ^ (1 : ℕ) : ℝ) +
    ((18 : ℝ) * B ^ (3 : ℕ) * r ^ (2 : ℕ) * q ^ (3 : ℕ) * l1 ^ (1 : ℕ) * lB ^ (1 : ℕ) : ℝ) +
    ((-18 : ℝ) * B ^ (3 : ℕ) * r ^ (2 : ℕ) * q ^ (3 : ℕ) * lB ^ (2 : ℕ) : ℝ) +
    ((18 : ℝ) * B ^ (3 : ℕ) * r ^ (2 : ℕ) * q ^ (1 : ℕ) * V ^ (1 : ℕ) * L2 ^ (1 : ℕ) * l0 ^ (1 : ℕ) * l1 ^ (1 : ℕ) : ℝ) +
    ((-18 : ℝ) * B ^ (3 : ℕ) * r ^ (2 : ℕ) * q ^ (1 : ℕ) * V ^ (1 : ℕ) * L2 ^ (1 : ℕ) * l0 ^ (1 : ℕ) * lB ^ (1 : ℕ) : ℝ)

noncomputable def pc7_2 : MvPolynomial (Fin 8) ℤ :=
  (MvPolynomial.C (-18 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (3 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (6 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (18 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (3 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (2 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (9 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (6 : Fin 8) ^ (2 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (-18 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (6 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (9 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (2 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (9 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (6 : Fin 8) ^ (2 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (-18 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (6 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (9 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (2 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (-9 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (3 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (6 : Fin 8) ^ (2 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (18 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (3 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (6 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (-9 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (3 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (2 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (-9 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (3 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (6 : Fin 8) ^ (2 : ℕ) : MvPolynomial (Fin 8) ℤ)

noncomputable def rc7_2 (B r q V L2 l0 l1 lB : ℝ) : ℝ :=
  ((-18 : ℝ) * B ^ (3 : ℕ) * r ^ (2 : ℕ) * q ^ (1 : ℕ) * V ^ (1 : ℕ) * L2 ^ (1 : ℕ) * l1 ^ (1 : ℕ) * lB ^ (1 : ℕ) : ℝ) +
    ((18 : ℝ) * B ^ (3 : ℕ) * r ^ (2 : ℕ) * q ^ (1 : ℕ) * V ^ (1 : ℕ) * L2 ^ (1 : ℕ) * lB ^ (2 : ℕ) : ℝ) +
    ((9 : ℝ) * B ^ (3 : ℕ) * r ^ (1 : ℕ) * q ^ (3 : ℕ) * l1 ^ (2 : ℕ) : ℝ) +
    ((-18 : ℝ) * B ^ (3 : ℕ) * r ^ (1 : ℕ) * q ^ (3 : ℕ) * l1 ^ (1 : ℕ) * lB ^ (1 : ℕ) : ℝ) +
    ((9 : ℝ) * B ^ (3 : ℕ) * r ^ (1 : ℕ) * q ^ (3 : ℕ) * lB ^ (2 : ℕ) : ℝ) +
    ((9 : ℝ) * B ^ (3 : ℕ) * r ^ (1 : ℕ) * q ^ (2 : ℕ) * L2 ^ (1 : ℕ) * l1 ^ (2 : ℕ) : ℝ) +
    ((-18 : ℝ) * B ^ (3 : ℕ) * r ^ (1 : ℕ) * q ^ (2 : ℕ) * L2 ^ (1 : ℕ) * l1 ^ (1 : ℕ) * lB ^ (1 : ℕ) : ℝ) +
    ((9 : ℝ) * B ^ (3 : ℕ) * r ^ (1 : ℕ) * q ^ (2 : ℕ) * L2 ^ (1 : ℕ) * lB ^ (2 : ℕ) : ℝ) +
    ((-9 : ℝ) * B ^ (3 : ℕ) * r ^ (1 : ℕ) * q ^ (1 : ℕ) * V ^ (1 : ℕ) * L2 ^ (1 : ℕ) * l1 ^ (2 : ℕ) : ℝ) +
    ((18 : ℝ) * B ^ (3 : ℕ) * r ^ (1 : ℕ) * q ^ (1 : ℕ) * V ^ (1 : ℕ) * L2 ^ (1 : ℕ) * l1 ^ (1 : ℕ) * lB ^ (1 : ℕ) : ℝ) +
    ((-9 : ℝ) * B ^ (3 : ℕ) * r ^ (1 : ℕ) * q ^ (1 : ℕ) * V ^ (1 : ℕ) * L2 ^ (1 : ℕ) * lB ^ (2 : ℕ) : ℝ) +
    ((-9 : ℝ) * B ^ (3 : ℕ) * r ^ (1 : ℕ) * V ^ (1 : ℕ) * L2 ^ (2 : ℕ) * l1 ^ (2 : ℕ) : ℝ)

noncomputable def pc7_3 : MvPolynomial (Fin 8) ℤ :=
  (MvPolynomial.C (18 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (3 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (6 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (-9 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (3 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (2 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (-18 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (5 : Fin 8) ^ (2 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (18 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (5 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (6 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (18 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (5 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (-18 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (6 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (15 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (5 : Fin 8) ^ (2 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (-24 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (5 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (12 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (5 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (9 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (2 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (-12 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (18 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (3 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (5 : Fin 8) ^ (2 : ℕ) : MvPolynomial (Fin 8) ℤ)

noncomputable def rc7_3 (B r q V L2 l0 l1 lB : ℝ) : ℝ :=
  ((18 : ℝ) * B ^ (3 : ℕ) * r ^ (1 : ℕ) * V ^ (1 : ℕ) * L2 ^ (2 : ℕ) * l1 ^ (1 : ℕ) * lB ^ (1 : ℕ) : ℝ) +
    ((-9 : ℝ) * B ^ (3 : ℕ) * r ^ (1 : ℕ) * V ^ (1 : ℕ) * L2 ^ (2 : ℕ) * lB ^ (2 : ℕ) : ℝ) +
    ((-18 : ℝ) * B ^ (2 : ℕ) * r ^ (2 : ℕ) * q ^ (3 : ℕ) * l0 ^ (2 : ℕ) : ℝ) +
    ((18 : ℝ) * B ^ (2 : ℕ) * r ^ (2 : ℕ) * q ^ (3 : ℕ) * l0 ^ (1 : ℕ) * l1 ^ (1 : ℕ) : ℝ) +
    ((18 : ℝ) * B ^ (2 : ℕ) * r ^ (2 : ℕ) * q ^ (3 : ℕ) * l0 ^ (1 : ℕ) * lB ^ (1 : ℕ) : ℝ) +
    ((-18 : ℝ) * B ^ (2 : ℕ) * r ^ (2 : ℕ) * q ^ (3 : ℕ) * l1 ^ (1 : ℕ) * lB ^ (1 : ℕ) : ℝ) +
    ((15 : ℝ) * B ^ (2 : ℕ) * r ^ (2 : ℕ) * q ^ (2 : ℕ) * L2 ^ (1 : ℕ) * l0 ^ (2 : ℕ) : ℝ) +
    ((-24 : ℝ) * B ^ (2 : ℕ) * r ^ (2 : ℕ) * q ^ (2 : ℕ) * L2 ^ (1 : ℕ) * l0 ^ (1 : ℕ) * lB ^ (1 : ℕ) : ℝ) +
    ((12 : ℝ) * B ^ (2 : ℕ) * r ^ (2 : ℕ) * q ^ (2 : ℕ) * L2 ^ (1 : ℕ) * l0 ^ (1 : ℕ) : ℝ) +
    ((9 : ℝ) * B ^ (2 : ℕ) * r ^ (2 : ℕ) * q ^ (2 : ℕ) * L2 ^ (1 : ℕ) * lB ^ (2 : ℕ) : ℝ) +
    ((-12 : ℝ) * B ^ (2 : ℕ) * r ^ (2 : ℕ) * q ^ (2 : ℕ) * L2 ^ (1 : ℕ) * lB ^ (1 : ℕ) : ℝ) +
    ((18 : ℝ) * B ^ (2 : ℕ) * r ^ (2 : ℕ) * q ^ (1 : ℕ) * V ^ (1 : ℕ) * L2 ^ (1 : ℕ) * l0 ^ (2 : ℕ) : ℝ)

noncomputable def pc7_4 : MvPolynomial (Fin 8) ℤ :=
  (MvPolynomial.C (-18 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (3 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (5 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (6 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (-18 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (3 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (5 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (18 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (3 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (6 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (-30 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (5 : Fin 8) ^ (2 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (48 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (5 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (-18 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (5 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (-18 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (2 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (18 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (-15 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (3 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (5 : Fin 8) ^ (2 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (24 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (3 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (5 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (-12 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (3 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (5 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (-9 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (3 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (2 : ℕ) : MvPolynomial (Fin 8) ℤ)

noncomputable def rc7_4 (B r q V L2 l0 l1 lB : ℝ) : ℝ :=
  ((-18 : ℝ) * B ^ (2 : ℕ) * r ^ (2 : ℕ) * q ^ (1 : ℕ) * V ^ (1 : ℕ) * L2 ^ (1 : ℕ) * l0 ^ (1 : ℕ) * l1 ^ (1 : ℕ) : ℝ) +
    ((-18 : ℝ) * B ^ (2 : ℕ) * r ^ (2 : ℕ) * q ^ (1 : ℕ) * V ^ (1 : ℕ) * L2 ^ (1 : ℕ) * l0 ^ (1 : ℕ) * lB ^ (1 : ℕ) : ℝ) +
    ((18 : ℝ) * B ^ (2 : ℕ) * r ^ (2 : ℕ) * q ^ (1 : ℕ) * V ^ (1 : ℕ) * L2 ^ (1 : ℕ) * l1 ^ (1 : ℕ) * lB ^ (1 : ℕ) : ℝ) +
    ((-30 : ℝ) * B ^ (2 : ℕ) * r ^ (2 : ℕ) * q ^ (1 : ℕ) * L2 ^ (2 : ℕ) * l0 ^ (2 : ℕ) : ℝ) +
    ((48 : ℝ) * B ^ (2 : ℕ) * r ^ (2 : ℕ) * q ^ (1 : ℕ) * L2 ^ (2 : ℕ) * l0 ^ (1 : ℕ) * lB ^ (1 : ℕ) : ℝ) +
    ((-18 : ℝ) * B ^ (2 : ℕ) * r ^ (2 : ℕ) * q ^ (1 : ℕ) * L2 ^ (2 : ℕ) * l0 ^ (1 : ℕ) : ℝ) +
    ((-18 : ℝ) * B ^ (2 : ℕ) * r ^ (2 : ℕ) * q ^ (1 : ℕ) * L2 ^ (2 : ℕ) * lB ^ (2 : ℕ) : ℝ) +
    ((18 : ℝ) * B ^ (2 : ℕ) * r ^ (2 : ℕ) * q ^ (1 : ℕ) * L2 ^ (2 : ℕ) * lB ^ (1 : ℕ) : ℝ) +
    ((-15 : ℝ) * B ^ (2 : ℕ) * r ^ (2 : ℕ) * V ^ (1 : ℕ) * L2 ^ (2 : ℕ) * l0 ^ (2 : ℕ) : ℝ) +
    ((24 : ℝ) * B ^ (2 : ℕ) * r ^ (2 : ℕ) * V ^ (1 : ℕ) * L2 ^ (2 : ℕ) * l0 ^ (1 : ℕ) * lB ^ (1 : ℕ) : ℝ) +
    ((-12 : ℝ) * B ^ (2 : ℕ) * r ^ (2 : ℕ) * V ^ (1 : ℕ) * L2 ^ (2 : ℕ) * l0 ^ (1 : ℕ) : ℝ) +
    ((-9 : ℝ) * B ^ (2 : ℕ) * r ^ (2 : ℕ) * V ^ (1 : ℕ) * L2 ^ (2 : ℕ) * lB ^ (2 : ℕ) : ℝ)

noncomputable def pc7_5 : MvPolynomial (Fin 8) ℤ :=
  (MvPolynomial.C (12 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (3 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (-30 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (5 : Fin 8) ^ (2 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (48 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (5 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (-18 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (5 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (-18 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (2 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (18 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (18 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (5 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (6 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (-18 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (5 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (-18 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (6 : Fin 8) ^ (2 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (18 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (6 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (-15 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (6 : Fin 8) ^ (2 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (24 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (6 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ)

noncomputable def rc7_5 (B r q V L2 l0 l1 lB : ℝ) : ℝ :=
  ((12 : ℝ) * B ^ (2 : ℕ) * r ^ (2 : ℕ) * V ^ (1 : ℕ) * L2 ^ (2 : ℕ) * lB ^ (1 : ℕ) : ℝ) +
    ((-30 : ℝ) * B ^ (2 : ℕ) * r ^ (2 : ℕ) * L2 ^ (3 : ℕ) * l0 ^ (2 : ℕ) : ℝ) +
    ((48 : ℝ) * B ^ (2 : ℕ) * r ^ (2 : ℕ) * L2 ^ (3 : ℕ) * l0 ^ (1 : ℕ) * lB ^ (1 : ℕ) : ℝ) +
    ((-18 : ℝ) * B ^ (2 : ℕ) * r ^ (2 : ℕ) * L2 ^ (3 : ℕ) * l0 ^ (1 : ℕ) : ℝ) +
    ((-18 : ℝ) * B ^ (2 : ℕ) * r ^ (2 : ℕ) * L2 ^ (3 : ℕ) * lB ^ (2 : ℕ) : ℝ) +
    ((18 : ℝ) * B ^ (2 : ℕ) * r ^ (2 : ℕ) * L2 ^ (3 : ℕ) * lB ^ (1 : ℕ) : ℝ) +
    ((18 : ℝ) * B ^ (2 : ℕ) * r ^ (1 : ℕ) * q ^ (3 : ℕ) * l0 ^ (1 : ℕ) * l1 ^ (1 : ℕ) : ℝ) +
    ((-18 : ℝ) * B ^ (2 : ℕ) * r ^ (1 : ℕ) * q ^ (3 : ℕ) * l0 ^ (1 : ℕ) * lB ^ (1 : ℕ) : ℝ) +
    ((-18 : ℝ) * B ^ (2 : ℕ) * r ^ (1 : ℕ) * q ^ (3 : ℕ) * l1 ^ (2 : ℕ) : ℝ) +
    ((18 : ℝ) * B ^ (2 : ℕ) * r ^ (1 : ℕ) * q ^ (3 : ℕ) * l1 ^ (1 : ℕ) * lB ^ (1 : ℕ) : ℝ) +
    ((-15 : ℝ) * B ^ (2 : ℕ) * r ^ (1 : ℕ) * q ^ (2 : ℕ) * L2 ^ (1 : ℕ) * l1 ^ (2 : ℕ) : ℝ) +
    ((24 : ℝ) * B ^ (2 : ℕ) * r ^ (1 : ℕ) * q ^ (2 : ℕ) * L2 ^ (1 : ℕ) * l1 ^ (1 : ℕ) * lB ^ (1 : ℕ) : ℝ)

noncomputable def pc7_6 : MvPolynomial (Fin 8) ℤ :=
  (MvPolynomial.C (-12 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (6 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (-9 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (2 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (12 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (-18 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (3 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (5 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (6 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (18 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (3 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (5 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (18 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (3 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (6 : Fin 8) ^ (2 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (-18 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (3 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (6 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (15 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (3 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (6 : Fin 8) ^ (2 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (-24 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (3 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (6 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (12 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (3 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (6 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (9 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (3 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (2 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (-12 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (3 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ)

noncomputable def rc7_6 (B r q V L2 l0 l1 lB : ℝ) : ℝ :=
  ((-12 : ℝ) * B ^ (2 : ℕ) * r ^ (1 : ℕ) * q ^ (2 : ℕ) * L2 ^ (1 : ℕ) * l1 ^ (1 : ℕ) : ℝ) +
    ((-9 : ℝ) * B ^ (2 : ℕ) * r ^ (1 : ℕ) * q ^ (2 : ℕ) * L2 ^ (1 : ℕ) * lB ^ (2 : ℕ) : ℝ) +
    ((12 : ℝ) * B ^ (2 : ℕ) * r ^ (1 : ℕ) * q ^ (2 : ℕ) * L2 ^ (1 : ℕ) * lB ^ (1 : ℕ) : ℝ) +
    ((-18 : ℝ) * B ^ (2 : ℕ) * r ^ (1 : ℕ) * q ^ (1 : ℕ) * V ^ (1 : ℕ) * L2 ^ (1 : ℕ) * l0 ^ (1 : ℕ) * l1 ^ (1 : ℕ) : ℝ) +
    ((18 : ℝ) * B ^ (2 : ℕ) * r ^ (1 : ℕ) * q ^ (1 : ℕ) * V ^ (1 : ℕ) * L2 ^ (1 : ℕ) * l0 ^ (1 : ℕ) * lB ^ (1 : ℕ) : ℝ) +
    ((18 : ℝ) * B ^ (2 : ℕ) * r ^ (1 : ℕ) * q ^ (1 : ℕ) * V ^ (1 : ℕ) * L2 ^ (1 : ℕ) * l1 ^ (2 : ℕ) : ℝ) +
    ((-18 : ℝ) * B ^ (2 : ℕ) * r ^ (1 : ℕ) * q ^ (1 : ℕ) * V ^ (1 : ℕ) * L2 ^ (1 : ℕ) * l1 ^ (1 : ℕ) * lB ^ (1 : ℕ) : ℝ) +
    ((15 : ℝ) * B ^ (2 : ℕ) * r ^ (1 : ℕ) * V ^ (1 : ℕ) * L2 ^ (2 : ℕ) * l1 ^ (2 : ℕ) : ℝ) +
    ((-24 : ℝ) * B ^ (2 : ℕ) * r ^ (1 : ℕ) * V ^ (1 : ℕ) * L2 ^ (2 : ℕ) * l1 ^ (1 : ℕ) * lB ^ (1 : ℕ) : ℝ) +
    ((12 : ℝ) * B ^ (2 : ℕ) * r ^ (1 : ℕ) * V ^ (1 : ℕ) * L2 ^ (2 : ℕ) * l1 ^ (1 : ℕ) : ℝ) +
    ((9 : ℝ) * B ^ (2 : ℕ) * r ^ (1 : ℕ) * V ^ (1 : ℕ) * L2 ^ (2 : ℕ) * lB ^ (2 : ℕ) : ℝ) +
    ((-12 : ℝ) * B ^ (2 : ℕ) * r ^ (1 : ℕ) * V ^ (1 : ℕ) * L2 ^ (2 : ℕ) * lB ^ (1 : ℕ) : ℝ)

noncomputable def pc7_7 : MvPolynomial (Fin 8) ℤ :=
  (MvPolynomial.C (9 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (5 : Fin 8) ^ (2 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (-18 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (5 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (6 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (9 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (6 : Fin 8) ^ (2 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (-6 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (5 : Fin 8) ^ (2 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (6 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (5 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (-12 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (5 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (6 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (6 : Fin 8) ^ (2 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (-6 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (6 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (12 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (6 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (-9 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (3 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (5 : Fin 8) ^ (2 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (18 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (3 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (5 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (6 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (-9 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (3 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (6 : Fin 8) ^ (2 : ℕ) : MvPolynomial (Fin 8) ℤ)

noncomputable def rc7_7 (B r q V L2 l0 l1 lB : ℝ) : ℝ :=
  ((9 : ℝ) * B ^ (1 : ℕ) * r ^ (1 : ℕ) * q ^ (3 : ℕ) * l0 ^ (2 : ℕ) : ℝ) +
    ((-18 : ℝ) * B ^ (1 : ℕ) * r ^ (1 : ℕ) * q ^ (3 : ℕ) * l0 ^ (1 : ℕ) * l1 ^ (1 : ℕ) : ℝ) +
    ((9 : ℝ) * B ^ (1 : ℕ) * r ^ (1 : ℕ) * q ^ (3 : ℕ) * l1 ^ (2 : ℕ) : ℝ) +
    ((-6 : ℝ) * B ^ (1 : ℕ) * r ^ (1 : ℕ) * q ^ (2 : ℕ) * L2 ^ (1 : ℕ) * l0 ^ (2 : ℕ) : ℝ) +
    ((6 : ℝ) * B ^ (1 : ℕ) * r ^ (1 : ℕ) * q ^ (2 : ℕ) * L2 ^ (1 : ℕ) * l0 ^ (1 : ℕ) * lB ^ (1 : ℕ) : ℝ) +
    ((-12 : ℝ) * B ^ (1 : ℕ) * r ^ (1 : ℕ) * q ^ (2 : ℕ) * L2 ^ (1 : ℕ) * l0 ^ (1 : ℕ) : ℝ) +
    ((6 : ℝ) * B ^ (1 : ℕ) * r ^ (1 : ℕ) * q ^ (2 : ℕ) * L2 ^ (1 : ℕ) * l1 ^ (2 : ℕ) : ℝ) +
    ((-6 : ℝ) * B ^ (1 : ℕ) * r ^ (1 : ℕ) * q ^ (2 : ℕ) * L2 ^ (1 : ℕ) * l1 ^ (1 : ℕ) * lB ^ (1 : ℕ) : ℝ) +
    ((12 : ℝ) * B ^ (1 : ℕ) * r ^ (1 : ℕ) * q ^ (2 : ℕ) * L2 ^ (1 : ℕ) * l1 ^ (1 : ℕ) : ℝ) +
    ((-9 : ℝ) * B ^ (1 : ℕ) * r ^ (1 : ℕ) * q ^ (1 : ℕ) * V ^ (1 : ℕ) * L2 ^ (1 : ℕ) * l0 ^ (2 : ℕ) : ℝ) +
    ((18 : ℝ) * B ^ (1 : ℕ) * r ^ (1 : ℕ) * q ^ (1 : ℕ) * V ^ (1 : ℕ) * L2 ^ (1 : ℕ) * l0 ^ (1 : ℕ) * l1 ^ (1 : ℕ) : ℝ) +
    ((-9 : ℝ) * B ^ (1 : ℕ) * r ^ (1 : ℕ) * q ^ (1 : ℕ) * V ^ (1 : ℕ) * L2 ^ (1 : ℕ) * l1 ^ (2 : ℕ) : ℝ)

noncomputable def pc7_8 : MvPolynomial (Fin 8) ℤ :=
  (MvPolynomial.C (12 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (5 : Fin 8) ^ (2 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (-12 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (5 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (18 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (5 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (3 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (2 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (-6 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (6 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (2 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (6 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (3 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (5 : Fin 8) ^ (2 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (-6 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (3 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (5 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (12 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (3 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (5 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (-6 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (3 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (6 : Fin 8) ^ (2 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (6 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (3 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (6 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (-12 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (3 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (6 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ)

noncomputable def rc7_8 (B r q V L2 l0 l1 lB : ℝ) : ℝ :=
  ((12 : ℝ) * B ^ (1 : ℕ) * r ^ (1 : ℕ) * q ^ (1 : ℕ) * L2 ^ (2 : ℕ) * l0 ^ (2 : ℕ) : ℝ) +
    ((-12 : ℝ) * B ^ (1 : ℕ) * r ^ (1 : ℕ) * q ^ (1 : ℕ) * L2 ^ (2 : ℕ) * l0 ^ (1 : ℕ) * lB ^ (1 : ℕ) : ℝ) +
    ((18 : ℝ) * B ^ (1 : ℕ) * r ^ (1 : ℕ) * q ^ (1 : ℕ) * L2 ^ (2 : ℕ) * l0 ^ (1 : ℕ) : ℝ) +
    ((3 : ℝ) * B ^ (1 : ℕ) * r ^ (1 : ℕ) * q ^ (1 : ℕ) * L2 ^ (2 : ℕ) * lB ^ (2 : ℕ) : ℝ) +
    ((-6 : ℝ) * B ^ (1 : ℕ) * r ^ (1 : ℕ) * q ^ (1 : ℕ) * L2 ^ (2 : ℕ) * lB ^ (1 : ℕ) : ℝ) +
    ((6 : ℝ) * B ^ (1 : ℕ) * r ^ (1 : ℕ) * q ^ (1 : ℕ) * L2 ^ (2 : ℕ) : ℝ) +
    ((6 : ℝ) * B ^ (1 : ℕ) * r ^ (1 : ℕ) * V ^ (1 : ℕ) * L2 ^ (2 : ℕ) * l0 ^ (2 : ℕ) : ℝ) +
    ((-6 : ℝ) * B ^ (1 : ℕ) * r ^ (1 : ℕ) * V ^ (1 : ℕ) * L2 ^ (2 : ℕ) * l0 ^ (1 : ℕ) * lB ^ (1 : ℕ) : ℝ) +
    ((12 : ℝ) * B ^ (1 : ℕ) * r ^ (1 : ℕ) * V ^ (1 : ℕ) * L2 ^ (2 : ℕ) * l0 ^ (1 : ℕ) : ℝ) +
    ((-6 : ℝ) * B ^ (1 : ℕ) * r ^ (1 : ℕ) * V ^ (1 : ℕ) * L2 ^ (2 : ℕ) * l1 ^ (2 : ℕ) : ℝ) +
    ((6 : ℝ) * B ^ (1 : ℕ) * r ^ (1 : ℕ) * V ^ (1 : ℕ) * L2 ^ (2 : ℕ) * l1 ^ (1 : ℕ) * lB ^ (1 : ℕ) : ℝ) +
    ((-12 : ℝ) * B ^ (1 : ℕ) * r ^ (1 : ℕ) * V ^ (1 : ℕ) * L2 ^ (2 : ℕ) * l1 ^ (1 : ℕ) : ℝ)

noncomputable def pc7_9 : MvPolynomial (Fin 8) ℤ :=
  (MvPolynomial.C (12 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (5 : Fin 8) ^ (2 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (-12 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (5 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (18 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (5 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (3 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (2 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (-6 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (6 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (3 : ℕ) : MvPolynomial (Fin 8) ℤ)

noncomputable def rc7_9 (B r q V L2 l0 l1 lB : ℝ) : ℝ :=
  ((12 : ℝ) * B ^ (1 : ℕ) * r ^ (1 : ℕ) * L2 ^ (3 : ℕ) * l0 ^ (2 : ℕ) : ℝ) +
    ((-12 : ℝ) * B ^ (1 : ℕ) * r ^ (1 : ℕ) * L2 ^ (3 : ℕ) * l0 ^ (1 : ℕ) * lB ^ (1 : ℕ) : ℝ) +
    ((18 : ℝ) * B ^ (1 : ℕ) * r ^ (1 : ℕ) * L2 ^ (3 : ℕ) * l0 ^ (1 : ℕ) : ℝ) +
    ((3 : ℝ) * B ^ (1 : ℕ) * r ^ (1 : ℕ) * L2 ^ (3 : ℕ) * lB ^ (2 : ℕ) : ℝ) +
    ((-6 : ℝ) * B ^ (1 : ℕ) * r ^ (1 : ℕ) * L2 ^ (3 : ℕ) * lB ^ (1 : ℕ) : ℝ) +
    ((6 : ℝ) * B ^ (1 : ℕ) * r ^ (1 : ℕ) * L2 ^ (3 : ℕ) : ℝ)

/-- block `c7` (coefficient of `ell`). -/
noncomputable def pc7 : MvPolynomial (Fin 8) ℤ :=
  pc7_0 + pc7_1 + pc7_2 + pc7_3 + pc7_4 + pc7_5 + pc7_6 + pc7_7 + pc7_8 + pc7_9

noncomputable def rc7 (B r q V L2 l0 l1 lB : ℝ) : ℝ :=
  rc7_0 B r q V L2 l0 l1 lB + rc7_1 B r q V L2 l0 l1 lB + rc7_2 B r q V L2 l0 l1 lB + rc7_3 B r q V L2 l0 l1 lB + rc7_4 B r q V L2 l0 l1 lB + rc7_5 B r q V L2 l0 l1 lB + rc7_6 B r q V L2 l0 l1 lB + rc7_7 B r q V L2 l0 l1 lB + rc7_8 B r q V L2 l0 l1 lB + rc7_9 B r q V L2 l0 l1 lB

noncomputable def pc10_0 : MvPolynomial (Fin 8) ℤ :=
  (MvPolynomial.C (18 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (5 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (-18 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (-18 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (3 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (5 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (18 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (3 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (-18 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (6 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (18 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (-18 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (6 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (18 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (18 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (3 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (6 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (-18 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (3 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (18 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (3 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (6 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (-18 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (3 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ)

noncomputable def rc10_0 (B r q V L2 l0 l1 lB : ℝ) : ℝ :=
  ((18 : ℝ) * B ^ (3 : ℕ) * r ^ (2 : ℕ) * q ^ (3 : ℕ) * l0 ^ (1 : ℕ) : ℝ) +
    ((-18 : ℝ) * B ^ (3 : ℕ) * r ^ (2 : ℕ) * q ^ (3 : ℕ) * lB ^ (1 : ℕ) : ℝ) +
    ((-18 : ℝ) * B ^ (3 : ℕ) * r ^ (2 : ℕ) * q ^ (1 : ℕ) * V ^ (1 : ℕ) * L2 ^ (1 : ℕ) * l0 ^ (1 : ℕ) : ℝ) +
    ((18 : ℝ) * B ^ (3 : ℕ) * r ^ (2 : ℕ) * q ^ (1 : ℕ) * V ^ (1 : ℕ) * L2 ^ (1 : ℕ) * lB ^ (1 : ℕ) : ℝ) +
    ((-18 : ℝ) * B ^ (3 : ℕ) * r ^ (1 : ℕ) * q ^ (3 : ℕ) * l1 ^ (1 : ℕ) : ℝ) +
    ((18 : ℝ) * B ^ (3 : ℕ) * r ^ (1 : ℕ) * q ^ (3 : ℕ) * lB ^ (1 : ℕ) : ℝ) +
    ((-18 : ℝ) * B ^ (3 : ℕ) * r ^ (1 : ℕ) * q ^ (2 : ℕ) * L2 ^ (1 : ℕ) * l1 ^ (1 : ℕ) : ℝ) +
    ((18 : ℝ) * B ^ (3 : ℕ) * r ^ (1 : ℕ) * q ^ (2 : ℕ) * L2 ^ (1 : ℕ) * lB ^ (1 : ℕ) : ℝ) +
    ((18 : ℝ) * B ^ (3 : ℕ) * r ^ (1 : ℕ) * q ^ (1 : ℕ) * V ^ (1 : ℕ) * L2 ^ (1 : ℕ) * l1 ^ (1 : ℕ) : ℝ) +
    ((-18 : ℝ) * B ^ (3 : ℕ) * r ^ (1 : ℕ) * q ^ (1 : ℕ) * V ^ (1 : ℕ) * L2 ^ (1 : ℕ) * lB ^ (1 : ℕ) : ℝ) +
    ((18 : ℝ) * B ^ (3 : ℕ) * r ^ (1 : ℕ) * V ^ (1 : ℕ) * L2 ^ (2 : ℕ) * l1 ^ (1 : ℕ) : ℝ) +
    ((-18 : ℝ) * B ^ (3 : ℕ) * r ^ (1 : ℕ) * V ^ (1 : ℕ) * L2 ^ (2 : ℕ) * lB ^ (1 : ℕ) : ℝ)

noncomputable def pc10_1 : MvPolynomial (Fin 8) ℤ :=
  (MvPolynomial.C (-18 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (5 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (18 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (18 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (3 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (5 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (-18 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (3 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (-36 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (5 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (36 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (6 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (-6 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (5 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (30 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (6 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (-18 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (12 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (36 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (3 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (5 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (-36 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (3 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (6 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ)

noncomputable def rc10_1 (B r q V L2 l0 l1 lB : ℝ) : ℝ :=
  ((-18 : ℝ) * B ^ (2 : ℕ) * r ^ (2 : ℕ) * q ^ (3 : ℕ) * l0 ^ (1 : ℕ) : ℝ) +
    ((18 : ℝ) * B ^ (2 : ℕ) * r ^ (2 : ℕ) * q ^ (3 : ℕ) * lB ^ (1 : ℕ) : ℝ) +
    ((18 : ℝ) * B ^ (2 : ℕ) * r ^ (2 : ℕ) * q ^ (1 : ℕ) * V ^ (1 : ℕ) * L2 ^ (1 : ℕ) * l0 ^ (1 : ℕ) : ℝ) +
    ((-18 : ℝ) * B ^ (2 : ℕ) * r ^ (2 : ℕ) * q ^ (1 : ℕ) * V ^ (1 : ℕ) * L2 ^ (1 : ℕ) * lB ^ (1 : ℕ) : ℝ) +
    ((-36 : ℝ) * B ^ (2 : ℕ) * r ^ (1 : ℕ) * q ^ (3 : ℕ) * l0 ^ (1 : ℕ) : ℝ) +
    ((36 : ℝ) * B ^ (2 : ℕ) * r ^ (1 : ℕ) * q ^ (3 : ℕ) * l1 ^ (1 : ℕ) : ℝ) +
    ((-6 : ℝ) * B ^ (2 : ℕ) * r ^ (1 : ℕ) * q ^ (2 : ℕ) * L2 ^ (1 : ℕ) * l0 ^ (1 : ℕ) : ℝ) +
    ((30 : ℝ) * B ^ (2 : ℕ) * r ^ (1 : ℕ) * q ^ (2 : ℕ) * L2 ^ (1 : ℕ) * l1 ^ (1 : ℕ) : ℝ) +
    ((-18 : ℝ) * B ^ (2 : ℕ) * r ^ (1 : ℕ) * q ^ (2 : ℕ) * L2 ^ (1 : ℕ) * lB ^ (1 : ℕ) : ℝ) +
    ((12 : ℝ) * B ^ (2 : ℕ) * r ^ (1 : ℕ) * q ^ (2 : ℕ) * L2 ^ (1 : ℕ) : ℝ) +
    ((36 : ℝ) * B ^ (2 : ℕ) * r ^ (1 : ℕ) * q ^ (1 : ℕ) * V ^ (1 : ℕ) * L2 ^ (1 : ℕ) * l0 ^ (1 : ℕ) : ℝ) +
    ((-36 : ℝ) * B ^ (2 : ℕ) * r ^ (1 : ℕ) * q ^ (1 : ℕ) * V ^ (1 : ℕ) * L2 ^ (1 : ℕ) * l1 ^ (1 : ℕ) : ℝ)

noncomputable def pc10_2 : MvPolynomial (Fin 8) ℤ :=
  (MvPolynomial.C (6 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (3 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (5 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (-30 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (3 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (6 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (18 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (3 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (-12 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (3 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (2 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (18 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (6 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (-18 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (18 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (6 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (-18 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (-18 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (3 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (6 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (18 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (3 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (-18 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (3 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (6 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (18 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (3 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ)

noncomputable def rc10_2 (B r q V L2 l0 l1 lB : ℝ) : ℝ :=
  ((6 : ℝ) * B ^ (2 : ℕ) * r ^ (1 : ℕ) * V ^ (1 : ℕ) * L2 ^ (2 : ℕ) * l0 ^ (1 : ℕ) : ℝ) +
    ((-30 : ℝ) * B ^ (2 : ℕ) * r ^ (1 : ℕ) * V ^ (1 : ℕ) * L2 ^ (2 : ℕ) * l1 ^ (1 : ℕ) : ℝ) +
    ((18 : ℝ) * B ^ (2 : ℕ) * r ^ (1 : ℕ) * V ^ (1 : ℕ) * L2 ^ (2 : ℕ) * lB ^ (1 : ℕ) : ℝ) +
    ((-12 : ℝ) * B ^ (2 : ℕ) * r ^ (1 : ℕ) * V ^ (1 : ℕ) * L2 ^ (2 : ℕ) : ℝ) +
    ((18 : ℝ) * B ^ (2 : ℕ) * q ^ (3 : ℕ) * l1 ^ (1 : ℕ) : ℝ) +
    ((-18 : ℝ) * B ^ (2 : ℕ) * q ^ (3 : ℕ) * lB ^ (1 : ℕ) : ℝ) +
    ((18 : ℝ) * B ^ (2 : ℕ) * q ^ (2 : ℕ) * L2 ^ (1 : ℕ) * l1 ^ (1 : ℕ) : ℝ) +
    ((-18 : ℝ) * B ^ (2 : ℕ) * q ^ (2 : ℕ) * L2 ^ (1 : ℕ) * lB ^ (1 : ℕ) : ℝ) +
    ((-18 : ℝ) * B ^ (2 : ℕ) * q ^ (1 : ℕ) * V ^ (1 : ℕ) * L2 ^ (1 : ℕ) * l1 ^ (1 : ℕ) : ℝ) +
    ((18 : ℝ) * B ^ (2 : ℕ) * q ^ (1 : ℕ) * V ^ (1 : ℕ) * L2 ^ (1 : ℕ) * lB ^ (1 : ℕ) : ℝ) +
    ((-18 : ℝ) * B ^ (2 : ℕ) * V ^ (1 : ℕ) * L2 ^ (2 : ℕ) * l1 ^ (1 : ℕ) : ℝ) +
    ((18 : ℝ) * B ^ (2 : ℕ) * V ^ (1 : ℕ) * L2 ^ (2 : ℕ) * lB ^ (1 : ℕ) : ℝ)

noncomputable def pc10_3 : MvPolynomial (Fin 8) ℤ :=
  (MvPolynomial.C (36 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (5 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (-18 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (6 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (-18 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (6 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (5 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (-12 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (6 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (-12 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (-36 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (3 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (5 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (18 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (3 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (6 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (18 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (3 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (-6 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (3 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (5 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (12 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (3 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (6 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (12 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (3 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (2 : ℕ) : MvPolynomial (Fin 8) ℤ)

noncomputable def rc10_3 (B r q V L2 l0 l1 lB : ℝ) : ℝ :=
  ((36 : ℝ) * B ^ (1 : ℕ) * r ^ (1 : ℕ) * q ^ (3 : ℕ) * l0 ^ (1 : ℕ) : ℝ) +
    ((-18 : ℝ) * B ^ (1 : ℕ) * r ^ (1 : ℕ) * q ^ (3 : ℕ) * l1 ^ (1 : ℕ) : ℝ) +
    ((-18 : ℝ) * B ^ (1 : ℕ) * r ^ (1 : ℕ) * q ^ (3 : ℕ) * lB ^ (1 : ℕ) : ℝ) +
    ((6 : ℝ) * B ^ (1 : ℕ) * r ^ (1 : ℕ) * q ^ (2 : ℕ) * L2 ^ (1 : ℕ) * l0 ^ (1 : ℕ) : ℝ) +
    ((-12 : ℝ) * B ^ (1 : ℕ) * r ^ (1 : ℕ) * q ^ (2 : ℕ) * L2 ^ (1 : ℕ) * l1 ^ (1 : ℕ) : ℝ) +
    ((-12 : ℝ) * B ^ (1 : ℕ) * r ^ (1 : ℕ) * q ^ (2 : ℕ) * L2 ^ (1 : ℕ) : ℝ) +
    ((-36 : ℝ) * B ^ (1 : ℕ) * r ^ (1 : ℕ) * q ^ (1 : ℕ) * V ^ (1 : ℕ) * L2 ^ (1 : ℕ) * l0 ^ (1 : ℕ) : ℝ) +
    ((18 : ℝ) * B ^ (1 : ℕ) * r ^ (1 : ℕ) * q ^ (1 : ℕ) * V ^ (1 : ℕ) * L2 ^ (1 : ℕ) * l1 ^ (1 : ℕ) : ℝ) +
    ((18 : ℝ) * B ^ (1 : ℕ) * r ^ (1 : ℕ) * q ^ (1 : ℕ) * V ^ (1 : ℕ) * L2 ^ (1 : ℕ) * lB ^ (1 : ℕ) : ℝ) +
    ((-6 : ℝ) * B ^ (1 : ℕ) * r ^ (1 : ℕ) * V ^ (1 : ℕ) * L2 ^ (2 : ℕ) * l0 ^ (1 : ℕ) : ℝ) +
    ((12 : ℝ) * B ^ (1 : ℕ) * r ^ (1 : ℕ) * V ^ (1 : ℕ) * L2 ^ (2 : ℕ) * l1 ^ (1 : ℕ) : ℝ) +
    ((12 : ℝ) * B ^ (1 : ℕ) * r ^ (1 : ℕ) * V ^ (1 : ℕ) * L2 ^ (2 : ℕ) : ℝ)

noncomputable def pc10_4 : MvPolynomial (Fin 8) ℤ :=
  (MvPolynomial.C (18 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (5 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (-36 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (6 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (18 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (6 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (5 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (-30 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (6 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (18 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (-12 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (-18 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (3 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (5 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (36 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (3 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (6 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (-18 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (3 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (-6 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (3 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (5 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (30 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (3 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (6 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ)

noncomputable def rc10_4 (B r q V L2 l0 l1 lB : ℝ) : ℝ :=
  ((18 : ℝ) * B ^ (1 : ℕ) * q ^ (3 : ℕ) * l0 ^ (1 : ℕ) : ℝ) +
    ((-36 : ℝ) * B ^ (1 : ℕ) * q ^ (3 : ℕ) * l1 ^ (1 : ℕ) : ℝ) +
    ((18 : ℝ) * B ^ (1 : ℕ) * q ^ (3 : ℕ) * lB ^ (1 : ℕ) : ℝ) +
    ((6 : ℝ) * B ^ (1 : ℕ) * q ^ (2 : ℕ) * L2 ^ (1 : ℕ) * l0 ^ (1 : ℕ) : ℝ) +
    ((-30 : ℝ) * B ^ (1 : ℕ) * q ^ (2 : ℕ) * L2 ^ (1 : ℕ) * l1 ^ (1 : ℕ) : ℝ) +
    ((18 : ℝ) * B ^ (1 : ℕ) * q ^ (2 : ℕ) * L2 ^ (1 : ℕ) * lB ^ (1 : ℕ) : ℝ) +
    ((-12 : ℝ) * B ^ (1 : ℕ) * q ^ (2 : ℕ) * L2 ^ (1 : ℕ) : ℝ) +
    ((-18 : ℝ) * B ^ (1 : ℕ) * q ^ (1 : ℕ) * V ^ (1 : ℕ) * L2 ^ (1 : ℕ) * l0 ^ (1 : ℕ) : ℝ) +
    ((36 : ℝ) * B ^ (1 : ℕ) * q ^ (1 : ℕ) * V ^ (1 : ℕ) * L2 ^ (1 : ℕ) * l1 ^ (1 : ℕ) : ℝ) +
    ((-18 : ℝ) * B ^ (1 : ℕ) * q ^ (1 : ℕ) * V ^ (1 : ℕ) * L2 ^ (1 : ℕ) * lB ^ (1 : ℕ) : ℝ) +
    ((-6 : ℝ) * B ^ (1 : ℕ) * V ^ (1 : ℕ) * L2 ^ (2 : ℕ) * l0 ^ (1 : ℕ) : ℝ) +
    ((30 : ℝ) * B ^ (1 : ℕ) * V ^ (1 : ℕ) * L2 ^ (2 : ℕ) * l1 ^ (1 : ℕ) : ℝ)

noncomputable def pc10_5 : MvPolynomial (Fin 8) ℤ :=
  (MvPolynomial.C (-18 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (3 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (12 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (3 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (2 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (-18 : ℤ) * MvPolynomial.X (2 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (5 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (18 : ℤ) * MvPolynomial.X (2 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (6 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (-6 : ℤ) * MvPolynomial.X (2 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (5 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (12 : ℤ) * MvPolynomial.X (2 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (6 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (12 : ℤ) * MvPolynomial.X (2 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (18 : ℤ) * MvPolynomial.X (2 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (3 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (5 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (-18 : ℤ) * MvPolynomial.X (2 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (3 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (6 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (6 : ℤ) * MvPolynomial.X (3 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (5 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (-12 : ℤ) * MvPolynomial.X (3 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (6 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (-12 : ℤ) * MvPolynomial.X (3 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (2 : ℕ) : MvPolynomial (Fin 8) ℤ)

noncomputable def rc10_5 (B r q V L2 l0 l1 lB : ℝ) : ℝ :=
  ((-18 : ℝ) * B ^ (1 : ℕ) * V ^ (1 : ℕ) * L2 ^ (2 : ℕ) * lB ^ (1 : ℕ) : ℝ) +
    ((12 : ℝ) * B ^ (1 : ℕ) * V ^ (1 : ℕ) * L2 ^ (2 : ℕ) : ℝ) +
    ((-18 : ℝ) * q ^ (3 : ℕ) * l0 ^ (1 : ℕ) : ℝ) +
    ((18 : ℝ) * q ^ (3 : ℕ) * l1 ^ (1 : ℕ) : ℝ) +
    ((-6 : ℝ) * q ^ (2 : ℕ) * L2 ^ (1 : ℕ) * l0 ^ (1 : ℕ) : ℝ) +
    ((12 : ℝ) * q ^ (2 : ℕ) * L2 ^ (1 : ℕ) * l1 ^ (1 : ℕ) : ℝ) +
    ((12 : ℝ) * q ^ (2 : ℕ) * L2 ^ (1 : ℕ) : ℝ) +
    ((18 : ℝ) * q ^ (1 : ℕ) * V ^ (1 : ℕ) * L2 ^ (1 : ℕ) * l0 ^ (1 : ℕ) : ℝ) +
    ((-18 : ℝ) * q ^ (1 : ℕ) * V ^ (1 : ℕ) * L2 ^ (1 : ℕ) * l1 ^ (1 : ℕ) : ℝ) +
    ((6 : ℝ) * V ^ (1 : ℕ) * L2 ^ (2 : ℕ) * l0 ^ (1 : ℕ) : ℝ) +
    ((-12 : ℝ) * V ^ (1 : ℕ) * L2 ^ (2 : ℕ) * l1 ^ (1 : ℕ) : ℝ) +
    ((-12 : ℝ) * V ^ (1 : ℕ) * L2 ^ (2 : ℕ) : ℝ)

/-- block `c10` (coefficient of `zeta2`). -/
noncomputable def pc10 : MvPolynomial (Fin 8) ℤ :=
  pc10_0 + pc10_1 + pc10_2 + pc10_3 + pc10_4 + pc10_5

noncomputable def rc10 (B r q V L2 l0 l1 lB : ℝ) : ℝ :=
  rc10_0 B r q V L2 l0 l1 lB + rc10_1 B r q V L2 l0 l1 lB + rc10_2 B r q V L2 l0 l1 lB + rc10_3 B r q V L2 l0 l1 lB + rc10_4 B r q V L2 l0 l1 lB + rc10_5 B r q V L2 l0 l1 lB

noncomputable def pc11_0 : MvPolynomial (Fin 8) ℤ :=
  (MvPolynomial.C (3 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (5 : Fin 8) ^ (3 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (-9 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (5 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (9 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (5 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (2 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (-3 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (3 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (-3 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (5 : Fin 8) ^ (3 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (9 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (5 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (-9 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (5 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (2 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (3 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (3 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (-3 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (3 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (5 : Fin 8) ^ (3 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (9 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (3 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (5 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (-9 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (3 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (5 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (2 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (3 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (3 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (3 : ℕ) : MvPolynomial (Fin 8) ℤ)

noncomputable def rc11_0 (B r q V L2 l0 l1 lB : ℝ) : ℝ :=
  ((3 : ℝ) * B ^ (3 : ℕ) * r ^ (3 : ℕ) * q ^ (3 : ℕ) * l0 ^ (3 : ℕ) : ℝ) +
    ((-9 : ℝ) * B ^ (3 : ℕ) * r ^ (3 : ℕ) * q ^ (3 : ℕ) * l0 ^ (2 : ℕ) * lB ^ (1 : ℕ) : ℝ) +
    ((9 : ℝ) * B ^ (3 : ℕ) * r ^ (3 : ℕ) * q ^ (3 : ℕ) * l0 ^ (1 : ℕ) * lB ^ (2 : ℕ) : ℝ) +
    ((-3 : ℝ) * B ^ (3 : ℕ) * r ^ (3 : ℕ) * q ^ (3 : ℕ) * lB ^ (3 : ℕ) : ℝ) +
    ((-3 : ℝ) * B ^ (3 : ℕ) * r ^ (3 : ℕ) * q ^ (2 : ℕ) * L2 ^ (1 : ℕ) * l0 ^ (3 : ℕ) : ℝ) +
    ((9 : ℝ) * B ^ (3 : ℕ) * r ^ (3 : ℕ) * q ^ (2 : ℕ) * L2 ^ (1 : ℕ) * l0 ^ (2 : ℕ) * lB ^ (1 : ℕ) : ℝ) +
    ((-9 : ℝ) * B ^ (3 : ℕ) * r ^ (3 : ℕ) * q ^ (2 : ℕ) * L2 ^ (1 : ℕ) * l0 ^ (1 : ℕ) * lB ^ (2 : ℕ) : ℝ) +
    ((3 : ℝ) * B ^ (3 : ℕ) * r ^ (3 : ℕ) * q ^ (2 : ℕ) * L2 ^ (1 : ℕ) * lB ^ (3 : ℕ) : ℝ) +
    ((-3 : ℝ) * B ^ (3 : ℕ) * r ^ (3 : ℕ) * q ^ (1 : ℕ) * V ^ (1 : ℕ) * L2 ^ (1 : ℕ) * l0 ^ (3 : ℕ) : ℝ) +
    ((9 : ℝ) * B ^ (3 : ℕ) * r ^ (3 : ℕ) * q ^ (1 : ℕ) * V ^ (1 : ℕ) * L2 ^ (1 : ℕ) * l0 ^ (2 : ℕ) * lB ^ (1 : ℕ) : ℝ) +
    ((-9 : ℝ) * B ^ (3 : ℕ) * r ^ (3 : ℕ) * q ^ (1 : ℕ) * V ^ (1 : ℕ) * L2 ^ (1 : ℕ) * l0 ^ (1 : ℕ) * lB ^ (2 : ℕ) : ℝ) +
    ((3 : ℝ) * B ^ (3 : ℕ) * r ^ (3 : ℕ) * q ^ (1 : ℕ) * V ^ (1 : ℕ) * L2 ^ (1 : ℕ) * lB ^ (3 : ℕ) : ℝ)

noncomputable def pc11_1 : MvPolynomial (Fin 8) ℤ :=
  (MvPolynomial.C (6 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (5 : Fin 8) ^ (3 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (-18 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (5 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (18 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (5 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (2 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (-6 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (3 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (3 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (3 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (5 : Fin 8) ^ (3 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (-9 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (3 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (5 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (9 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (3 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (5 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (2 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (-3 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (3 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (3 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (6 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (5 : Fin 8) ^ (3 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (-18 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (5 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (18 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (5 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (2 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (-6 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (3 : ℕ) : MvPolynomial (Fin 8) ℤ)

noncomputable def rc11_1 (B r q V L2 l0 l1 lB : ℝ) : ℝ :=
  ((6 : ℝ) * B ^ (3 : ℕ) * r ^ (3 : ℕ) * q ^ (1 : ℕ) * L2 ^ (2 : ℕ) * l0 ^ (3 : ℕ) : ℝ) +
    ((-18 : ℝ) * B ^ (3 : ℕ) * r ^ (3 : ℕ) * q ^ (1 : ℕ) * L2 ^ (2 : ℕ) * l0 ^ (2 : ℕ) * lB ^ (1 : ℕ) : ℝ) +
    ((18 : ℝ) * B ^ (3 : ℕ) * r ^ (3 : ℕ) * q ^ (1 : ℕ) * L2 ^ (2 : ℕ) * l0 ^ (1 : ℕ) * lB ^ (2 : ℕ) : ℝ) +
    ((-6 : ℝ) * B ^ (3 : ℕ) * r ^ (3 : ℕ) * q ^ (1 : ℕ) * L2 ^ (2 : ℕ) * lB ^ (3 : ℕ) : ℝ) +
    ((3 : ℝ) * B ^ (3 : ℕ) * r ^ (3 : ℕ) * V ^ (1 : ℕ) * L2 ^ (2 : ℕ) * l0 ^ (3 : ℕ) : ℝ) +
    ((-9 : ℝ) * B ^ (3 : ℕ) * r ^ (3 : ℕ) * V ^ (1 : ℕ) * L2 ^ (2 : ℕ) * l0 ^ (2 : ℕ) * lB ^ (1 : ℕ) : ℝ) +
    ((9 : ℝ) * B ^ (3 : ℕ) * r ^ (3 : ℕ) * V ^ (1 : ℕ) * L2 ^ (2 : ℕ) * l0 ^ (1 : ℕ) * lB ^ (2 : ℕ) : ℝ) +
    ((-3 : ℝ) * B ^ (3 : ℕ) * r ^ (3 : ℕ) * V ^ (1 : ℕ) * L2 ^ (2 : ℕ) * lB ^ (3 : ℕ) : ℝ) +
    ((6 : ℝ) * B ^ (3 : ℕ) * r ^ (3 : ℕ) * L2 ^ (3 : ℕ) * l0 ^ (3 : ℕ) : ℝ) +
    ((-18 : ℝ) * B ^ (3 : ℕ) * r ^ (3 : ℕ) * L2 ^ (3 : ℕ) * l0 ^ (2 : ℕ) * lB ^ (1 : ℕ) : ℝ) +
    ((18 : ℝ) * B ^ (3 : ℕ) * r ^ (3 : ℕ) * L2 ^ (3 : ℕ) * l0 ^ (1 : ℕ) * lB ^ (2 : ℕ) : ℝ) +
    ((-6 : ℝ) * B ^ (3 : ℕ) * r ^ (3 : ℕ) * L2 ^ (3 : ℕ) * lB ^ (3 : ℕ) : ℝ)

noncomputable def pc11_2 : MvPolynomial (Fin 8) ℤ :=
  (MvPolynomial.C (-9 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (5 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (6 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (9 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (5 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (18 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (5 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (6 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (-18 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (5 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (2 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (-9 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (6 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (2 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (9 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (3 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (9 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (3 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (5 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (6 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (-9 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (3 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (5 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (-18 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (3 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (5 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (6 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (18 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (3 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (5 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (2 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (9 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (3 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (6 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (2 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (-9 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (3 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (3 : ℕ) : MvPolynomial (Fin 8) ℤ)

noncomputable def rc11_2 (B r q V L2 l0 l1 lB : ℝ) : ℝ :=
  ((-9 : ℝ) * B ^ (3 : ℕ) * r ^ (2 : ℕ) * q ^ (3 : ℕ) * l0 ^ (2 : ℕ) * l1 ^ (1 : ℕ) : ℝ) +
    ((9 : ℝ) * B ^ (3 : ℕ) * r ^ (2 : ℕ) * q ^ (3 : ℕ) * l0 ^ (2 : ℕ) * lB ^ (1 : ℕ) : ℝ) +
    ((18 : ℝ) * B ^ (3 : ℕ) * r ^ (2 : ℕ) * q ^ (3 : ℕ) * l0 ^ (1 : ℕ) * l1 ^ (1 : ℕ) * lB ^ (1 : ℕ) : ℝ) +
    ((-18 : ℝ) * B ^ (3 : ℕ) * r ^ (2 : ℕ) * q ^ (3 : ℕ) * l0 ^ (1 : ℕ) * lB ^ (2 : ℕ) : ℝ) +
    ((-9 : ℝ) * B ^ (3 : ℕ) * r ^ (2 : ℕ) * q ^ (3 : ℕ) * l1 ^ (1 : ℕ) * lB ^ (2 : ℕ) : ℝ) +
    ((9 : ℝ) * B ^ (3 : ℕ) * r ^ (2 : ℕ) * q ^ (3 : ℕ) * lB ^ (3 : ℕ) : ℝ) +
    ((9 : ℝ) * B ^ (3 : ℕ) * r ^ (2 : ℕ) * q ^ (1 : ℕ) * V ^ (1 : ℕ) * L2 ^ (1 : ℕ) * l0 ^ (2 : ℕ) * l1 ^ (1 : ℕ) : ℝ) +
    ((-9 : ℝ) * B ^ (3 : ℕ) * r ^ (2 : ℕ) * q ^ (1 : ℕ) * V ^ (1 : ℕ) * L2 ^ (1 : ℕ) * l0 ^ (2 : ℕ) * lB ^ (1 : ℕ) : ℝ) +
    ((-18 : ℝ) * B ^ (3 : ℕ) * r ^ (2 : ℕ) * q ^ (1 : ℕ) * V ^ (1 : ℕ) * L2 ^ (1 : ℕ) * l0 ^ (1 : ℕ) * l1 ^ (1 : ℕ) * lB ^ (1 : ℕ) : ℝ) +
    ((18 : ℝ) * B ^ (3 : ℕ) * r ^ (2 : ℕ) * q ^ (1 : ℕ) * V ^ (1 : ℕ) * L2 ^ (1 : ℕ) * l0 ^ (1 : ℕ) * lB ^ (2 : ℕ) : ℝ) +
    ((9 : ℝ) * B ^ (3 : ℕ) * r ^ (2 : ℕ) * q ^ (1 : ℕ) * V ^ (1 : ℕ) * L2 ^ (1 : ℕ) * l1 ^ (1 : ℕ) * lB ^ (2 : ℕ) : ℝ) +
    ((-9 : ℝ) * B ^ (3 : ℕ) * r ^ (2 : ℕ) * q ^ (1 : ℕ) * V ^ (1 : ℕ) * L2 ^ (1 : ℕ) * lB ^ (3 : ℕ) : ℝ)

noncomputable def pc11_3 : MvPolynomial (Fin 8) ℤ :=
  (MvPolynomial.C (9 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (5 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (6 : Fin 8) ^ (2 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (-18 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (5 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (6 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (9 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (5 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (2 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (-9 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (6 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (18 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (6 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (2 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (-9 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (3 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (9 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (5 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (6 : Fin 8) ^ (2 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (-18 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (5 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (6 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (9 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (5 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (2 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (-9 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (6 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (18 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (6 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (2 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (-9 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (3 : ℕ) : MvPolynomial (Fin 8) ℤ)

noncomputable def rc11_3 (B r q V L2 l0 l1 lB : ℝ) : ℝ :=
  ((9 : ℝ) * B ^ (3 : ℕ) * r ^ (1 : ℕ) * q ^ (3 : ℕ) * l0 ^ (1 : ℕ) * l1 ^ (2 : ℕ) : ℝ) +
    ((-18 : ℝ) * B ^ (3 : ℕ) * r ^ (1 : ℕ) * q ^ (3 : ℕ) * l0 ^ (1 : ℕ) * l1 ^ (1 : ℕ) * lB ^ (1 : ℕ) : ℝ) +
    ((9 : ℝ) * B ^ (3 : ℕ) * r ^ (1 : ℕ) * q ^ (3 : ℕ) * l0 ^ (1 : ℕ) * lB ^ (2 : ℕ) : ℝ) +
    ((-9 : ℝ) * B ^ (3 : ℕ) * r ^ (1 : ℕ) * q ^ (3 : ℕ) * l1 ^ (2 : ℕ) * lB ^ (1 : ℕ) : ℝ) +
    ((18 : ℝ) * B ^ (3 : ℕ) * r ^ (1 : ℕ) * q ^ (3 : ℕ) * l1 ^ (1 : ℕ) * lB ^ (2 : ℕ) : ℝ) +
    ((-9 : ℝ) * B ^ (3 : ℕ) * r ^ (1 : ℕ) * q ^ (3 : ℕ) * lB ^ (3 : ℕ) : ℝ) +
    ((9 : ℝ) * B ^ (3 : ℕ) * r ^ (1 : ℕ) * q ^ (2 : ℕ) * L2 ^ (1 : ℕ) * l0 ^ (1 : ℕ) * l1 ^ (2 : ℕ) : ℝ) +
    ((-18 : ℝ) * B ^ (3 : ℕ) * r ^ (1 : ℕ) * q ^ (2 : ℕ) * L2 ^ (1 : ℕ) * l0 ^ (1 : ℕ) * l1 ^ (1 : ℕ) * lB ^ (1 : ℕ) : ℝ) +
    ((9 : ℝ) * B ^ (3 : ℕ) * r ^ (1 : ℕ) * q ^ (2 : ℕ) * L2 ^ (1 : ℕ) * l0 ^ (1 : ℕ) * lB ^ (2 : ℕ) : ℝ) +
    ((-9 : ℝ) * B ^ (3 : ℕ) * r ^ (1 : ℕ) * q ^ (2 : ℕ) * L2 ^ (1 : ℕ) * l1 ^ (2 : ℕ) * lB ^ (1 : ℕ) : ℝ) +
    ((18 : ℝ) * B ^ (3 : ℕ) * r ^ (1 : ℕ) * q ^ (2 : ℕ) * L2 ^ (1 : ℕ) * l1 ^ (1 : ℕ) * lB ^ (2 : ℕ) : ℝ) +
    ((-9 : ℝ) * B ^ (3 : ℕ) * r ^ (1 : ℕ) * q ^ (2 : ℕ) * L2 ^ (1 : ℕ) * lB ^ (3 : ℕ) : ℝ)

noncomputable def pc11_4 : MvPolynomial (Fin 8) ℤ :=
  (MvPolynomial.C (-9 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (3 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (5 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (6 : Fin 8) ^ (2 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (18 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (3 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (5 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (6 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (-9 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (3 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (5 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (2 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (9 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (3 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (6 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (-18 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (3 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (6 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (2 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (9 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (3 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (3 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (-9 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (3 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (5 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (6 : Fin 8) ^ (2 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (18 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (3 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (5 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (6 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (-9 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (3 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (5 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (2 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (9 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (3 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (6 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (-18 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (3 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (6 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (2 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (9 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (3 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (3 : ℕ) : MvPolynomial (Fin 8) ℤ)

noncomputable def rc11_4 (B r q V L2 l0 l1 lB : ℝ) : ℝ :=
  ((-9 : ℝ) * B ^ (3 : ℕ) * r ^ (1 : ℕ) * q ^ (1 : ℕ) * V ^ (1 : ℕ) * L2 ^ (1 : ℕ) * l0 ^ (1 : ℕ) * l1 ^ (2 : ℕ) : ℝ) +
    ((18 : ℝ) * B ^ (3 : ℕ) * r ^ (1 : ℕ) * q ^ (1 : ℕ) * V ^ (1 : ℕ) * L2 ^ (1 : ℕ) * l0 ^ (1 : ℕ) * l1 ^ (1 : ℕ) * lB ^ (1 : ℕ) : ℝ) +
    ((-9 : ℝ) * B ^ (3 : ℕ) * r ^ (1 : ℕ) * q ^ (1 : ℕ) * V ^ (1 : ℕ) * L2 ^ (1 : ℕ) * l0 ^ (1 : ℕ) * lB ^ (2 : ℕ) : ℝ) +
    ((9 : ℝ) * B ^ (3 : ℕ) * r ^ (1 : ℕ) * q ^ (1 : ℕ) * V ^ (1 : ℕ) * L2 ^ (1 : ℕ) * l1 ^ (2 : ℕ) * lB ^ (1 : ℕ) : ℝ) +
    ((-18 : ℝ) * B ^ (3 : ℕ) * r ^ (1 : ℕ) * q ^ (1 : ℕ) * V ^ (1 : ℕ) * L2 ^ (1 : ℕ) * l1 ^ (1 : ℕ) * lB ^ (2 : ℕ) : ℝ) +
    ((9 : ℝ) * B ^ (3 : ℕ) * r ^ (1 : ℕ) * q ^ (1 : ℕ) * V ^ (1 : ℕ) * L2 ^ (1 : ℕ) * lB ^ (3 : ℕ) : ℝ) +
    ((-9 : ℝ) * B ^ (3 : ℕ) * r ^ (1 : ℕ) * V ^ (1 : ℕ) * L2 ^ (2 : ℕ) * l0 ^ (1 : ℕ) * l1 ^ (2 : ℕ) : ℝ) +
    ((18 : ℝ) * B ^ (3 : ℕ) * r ^ (1 : ℕ) * V ^ (1 : ℕ) * L2 ^ (2 : ℕ) * l0 ^ (1 : ℕ) * l1 ^ (1 : ℕ) * lB ^ (1 : ℕ) : ℝ) +
    ((-9 : ℝ) * B ^ (3 : ℕ) * r ^ (1 : ℕ) * V ^ (1 : ℕ) * L2 ^ (2 : ℕ) * l0 ^ (1 : ℕ) * lB ^ (2 : ℕ) : ℝ) +
    ((9 : ℝ) * B ^ (3 : ℕ) * r ^ (1 : ℕ) * V ^ (1 : ℕ) * L2 ^ (2 : ℕ) * l1 ^ (2 : ℕ) * lB ^ (1 : ℕ) : ℝ) +
    ((-18 : ℝ) * B ^ (3 : ℕ) * r ^ (1 : ℕ) * V ^ (1 : ℕ) * L2 ^ (2 : ℕ) * l1 ^ (1 : ℕ) * lB ^ (2 : ℕ) : ℝ) +
    ((9 : ℝ) * B ^ (3 : ℕ) * r ^ (1 : ℕ) * V ^ (1 : ℕ) * L2 ^ (2 : ℕ) * lB ^ (3 : ℕ) : ℝ)

noncomputable def pc11_5 : MvPolynomial (Fin 8) ℤ :=
  (MvPolynomial.C (-3 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (6 : Fin 8) ^ (3 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (9 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (6 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (-9 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (6 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (2 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (3 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (3 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (-6 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (6 : Fin 8) ^ (3 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (18 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (6 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (-18 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (6 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (2 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (6 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (3 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (3 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (3 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (6 : Fin 8) ^ (3 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (-9 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (3 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (6 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (9 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (3 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (6 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (2 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (-3 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (3 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (3 : ℕ) : MvPolynomial (Fin 8) ℤ)

noncomputable def rc11_5 (B r q V L2 l0 l1 lB : ℝ) : ℝ :=
  ((-3 : ℝ) * B ^ (3 : ℕ) * q ^ (3 : ℕ) * l1 ^ (3 : ℕ) : ℝ) +
    ((9 : ℝ) * B ^ (3 : ℕ) * q ^ (3 : ℕ) * l1 ^ (2 : ℕ) * lB ^ (1 : ℕ) : ℝ) +
    ((-9 : ℝ) * B ^ (3 : ℕ) * q ^ (3 : ℕ) * l1 ^ (1 : ℕ) * lB ^ (2 : ℕ) : ℝ) +
    ((3 : ℝ) * B ^ (3 : ℕ) * q ^ (3 : ℕ) * lB ^ (3 : ℕ) : ℝ) +
    ((-6 : ℝ) * B ^ (3 : ℕ) * q ^ (2 : ℕ) * L2 ^ (1 : ℕ) * l1 ^ (3 : ℕ) : ℝ) +
    ((18 : ℝ) * B ^ (3 : ℕ) * q ^ (2 : ℕ) * L2 ^ (1 : ℕ) * l1 ^ (2 : ℕ) * lB ^ (1 : ℕ) : ℝ) +
    ((-18 : ℝ) * B ^ (3 : ℕ) * q ^ (2 : ℕ) * L2 ^ (1 : ℕ) * l1 ^ (1 : ℕ) * lB ^ (2 : ℕ) : ℝ) +
    ((6 : ℝ) * B ^ (3 : ℕ) * q ^ (2 : ℕ) * L2 ^ (1 : ℕ) * lB ^ (3 : ℕ) : ℝ) +
    ((3 : ℝ) * B ^ (3 : ℕ) * q ^ (1 : ℕ) * V ^ (1 : ℕ) * L2 ^ (1 : ℕ) * l1 ^ (3 : ℕ) : ℝ) +
    ((-9 : ℝ) * B ^ (3 : ℕ) * q ^ (1 : ℕ) * V ^ (1 : ℕ) * L2 ^ (1 : ℕ) * l1 ^ (2 : ℕ) * lB ^ (1 : ℕ) : ℝ) +
    ((9 : ℝ) * B ^ (3 : ℕ) * q ^ (1 : ℕ) * V ^ (1 : ℕ) * L2 ^ (1 : ℕ) * l1 ^ (1 : ℕ) * lB ^ (2 : ℕ) : ℝ) +
    ((-3 : ℝ) * B ^ (3 : ℕ) * q ^ (1 : ℕ) * V ^ (1 : ℕ) * L2 ^ (1 : ℕ) * lB ^ (3 : ℕ) : ℝ)

noncomputable def pc11_6 : MvPolynomial (Fin 8) ℤ :=
  (MvPolynomial.C (-6 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (6 : Fin 8) ^ (3 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (18 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (6 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (-18 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (6 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (2 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (6 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (3 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (6 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (3 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (6 : Fin 8) ^ (3 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (-18 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (3 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (6 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (18 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (3 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (6 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (2 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (-6 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (3 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (3 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (-9 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (5 : Fin 8) ^ (3 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (9 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (5 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (6 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (18 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (5 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (-18 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (5 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (6 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ)

noncomputable def rc11_6 (B r q V L2 l0 l1 lB : ℝ) : ℝ :=
  ((-6 : ℝ) * B ^ (3 : ℕ) * q ^ (1 : ℕ) * L2 ^ (2 : ℕ) * l1 ^ (3 : ℕ) : ℝ) +
    ((18 : ℝ) * B ^ (3 : ℕ) * q ^ (1 : ℕ) * L2 ^ (2 : ℕ) * l1 ^ (2 : ℕ) * lB ^ (1 : ℕ) : ℝ) +
    ((-18 : ℝ) * B ^ (3 : ℕ) * q ^ (1 : ℕ) * L2 ^ (2 : ℕ) * l1 ^ (1 : ℕ) * lB ^ (2 : ℕ) : ℝ) +
    ((6 : ℝ) * B ^ (3 : ℕ) * q ^ (1 : ℕ) * L2 ^ (2 : ℕ) * lB ^ (3 : ℕ) : ℝ) +
    ((6 : ℝ) * B ^ (3 : ℕ) * V ^ (1 : ℕ) * L2 ^ (2 : ℕ) * l1 ^ (3 : ℕ) : ℝ) +
    ((-18 : ℝ) * B ^ (3 : ℕ) * V ^ (1 : ℕ) * L2 ^ (2 : ℕ) * l1 ^ (2 : ℕ) * lB ^ (1 : ℕ) : ℝ) +
    ((18 : ℝ) * B ^ (3 : ℕ) * V ^ (1 : ℕ) * L2 ^ (2 : ℕ) * l1 ^ (1 : ℕ) * lB ^ (2 : ℕ) : ℝ) +
    ((-6 : ℝ) * B ^ (3 : ℕ) * V ^ (1 : ℕ) * L2 ^ (2 : ℕ) * lB ^ (3 : ℕ) : ℝ) +
    ((-9 : ℝ) * B ^ (2 : ℕ) * r ^ (2 : ℕ) * q ^ (3 : ℕ) * l0 ^ (3 : ℕ) : ℝ) +
    ((9 : ℝ) * B ^ (2 : ℕ) * r ^ (2 : ℕ) * q ^ (3 : ℕ) * l0 ^ (2 : ℕ) * l1 ^ (1 : ℕ) : ℝ) +
    ((18 : ℝ) * B ^ (2 : ℕ) * r ^ (2 : ℕ) * q ^ (3 : ℕ) * l0 ^ (2 : ℕ) * lB ^ (1 : ℕ) : ℝ) +
    ((-18 : ℝ) * B ^ (2 : ℕ) * r ^ (2 : ℕ) * q ^ (3 : ℕ) * l0 ^ (1 : ℕ) * l1 ^ (1 : ℕ) * lB ^ (1 : ℕ) : ℝ)

noncomputable def pc11_7 : MvPolynomial (Fin 8) ℤ :=
  (MvPolynomial.C (-9 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (5 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (2 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (9 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (6 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (2 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (6 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (5 : Fin 8) ^ (3 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (-15 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (5 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (6 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (5 : Fin 8) ^ (2 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (12 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (5 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (2 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (-12 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (5 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (-3 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (3 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (6 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (2 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (9 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (3 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (5 : Fin 8) ^ (3 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (-9 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (3 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (5 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (6 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (-18 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (3 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (5 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ)

noncomputable def rc11_7 (B r q V L2 l0 l1 lB : ℝ) : ℝ :=
  ((-9 : ℝ) * B ^ (2 : ℕ) * r ^ (2 : ℕ) * q ^ (3 : ℕ) * l0 ^ (1 : ℕ) * lB ^ (2 : ℕ) : ℝ) +
    ((9 : ℝ) * B ^ (2 : ℕ) * r ^ (2 : ℕ) * q ^ (3 : ℕ) * l1 ^ (1 : ℕ) * lB ^ (2 : ℕ) : ℝ) +
    ((6 : ℝ) * B ^ (2 : ℕ) * r ^ (2 : ℕ) * q ^ (2 : ℕ) * L2 ^ (1 : ℕ) * l0 ^ (3 : ℕ) : ℝ) +
    ((-15 : ℝ) * B ^ (2 : ℕ) * r ^ (2 : ℕ) * q ^ (2 : ℕ) * L2 ^ (1 : ℕ) * l0 ^ (2 : ℕ) * lB ^ (1 : ℕ) : ℝ) +
    ((6 : ℝ) * B ^ (2 : ℕ) * r ^ (2 : ℕ) * q ^ (2 : ℕ) * L2 ^ (1 : ℕ) * l0 ^ (2 : ℕ) : ℝ) +
    ((12 : ℝ) * B ^ (2 : ℕ) * r ^ (2 : ℕ) * q ^ (2 : ℕ) * L2 ^ (1 : ℕ) * l0 ^ (1 : ℕ) * lB ^ (2 : ℕ) : ℝ) +
    ((-12 : ℝ) * B ^ (2 : ℕ) * r ^ (2 : ℕ) * q ^ (2 : ℕ) * L2 ^ (1 : ℕ) * l0 ^ (1 : ℕ) * lB ^ (1 : ℕ) : ℝ) +
    ((-3 : ℝ) * B ^ (2 : ℕ) * r ^ (2 : ℕ) * q ^ (2 : ℕ) * L2 ^ (1 : ℕ) * lB ^ (3 : ℕ) : ℝ) +
    ((6 : ℝ) * B ^ (2 : ℕ) * r ^ (2 : ℕ) * q ^ (2 : ℕ) * L2 ^ (1 : ℕ) * lB ^ (2 : ℕ) : ℝ) +
    ((9 : ℝ) * B ^ (2 : ℕ) * r ^ (2 : ℕ) * q ^ (1 : ℕ) * V ^ (1 : ℕ) * L2 ^ (1 : ℕ) * l0 ^ (3 : ℕ) : ℝ) +
    ((-9 : ℝ) * B ^ (2 : ℕ) * r ^ (2 : ℕ) * q ^ (1 : ℕ) * V ^ (1 : ℕ) * L2 ^ (1 : ℕ) * l0 ^ (2 : ℕ) * l1 ^ (1 : ℕ) : ℝ) +
    ((-18 : ℝ) * B ^ (2 : ℕ) * r ^ (2 : ℕ) * q ^ (1 : ℕ) * V ^ (1 : ℕ) * L2 ^ (1 : ℕ) * l0 ^ (2 : ℕ) * lB ^ (1 : ℕ) : ℝ)

noncomputable def pc11_8 : MvPolynomial (Fin 8) ℤ :=
  (MvPolynomial.C (18 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (3 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (5 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (6 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (9 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (3 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (5 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (2 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (-9 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (3 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (6 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (2 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (-12 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (5 : Fin 8) ^ (3 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (30 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (5 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (-9 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (5 : Fin 8) ^ (2 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (-24 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (5 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (2 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (18 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (5 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (6 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (3 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (-9 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (2 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (-6 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (3 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (5 : Fin 8) ^ (3 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (15 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (3 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (5 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ)

noncomputable def rc11_8 (B r q V L2 l0 l1 lB : ℝ) : ℝ :=
  ((18 : ℝ) * B ^ (2 : ℕ) * r ^ (2 : ℕ) * q ^ (1 : ℕ) * V ^ (1 : ℕ) * L2 ^ (1 : ℕ) * l0 ^ (1 : ℕ) * l1 ^ (1 : ℕ) * lB ^ (1 : ℕ) : ℝ) +
    ((9 : ℝ) * B ^ (2 : ℕ) * r ^ (2 : ℕ) * q ^ (1 : ℕ) * V ^ (1 : ℕ) * L2 ^ (1 : ℕ) * l0 ^ (1 : ℕ) * lB ^ (2 : ℕ) : ℝ) +
    ((-9 : ℝ) * B ^ (2 : ℕ) * r ^ (2 : ℕ) * q ^ (1 : ℕ) * V ^ (1 : ℕ) * L2 ^ (1 : ℕ) * l1 ^ (1 : ℕ) * lB ^ (2 : ℕ) : ℝ) +
    ((-12 : ℝ) * B ^ (2 : ℕ) * r ^ (2 : ℕ) * q ^ (1 : ℕ) * L2 ^ (2 : ℕ) * l0 ^ (3 : ℕ) : ℝ) +
    ((30 : ℝ) * B ^ (2 : ℕ) * r ^ (2 : ℕ) * q ^ (1 : ℕ) * L2 ^ (2 : ℕ) * l0 ^ (2 : ℕ) * lB ^ (1 : ℕ) : ℝ) +
    ((-9 : ℝ) * B ^ (2 : ℕ) * r ^ (2 : ℕ) * q ^ (1 : ℕ) * L2 ^ (2 : ℕ) * l0 ^ (2 : ℕ) : ℝ) +
    ((-24 : ℝ) * B ^ (2 : ℕ) * r ^ (2 : ℕ) * q ^ (1 : ℕ) * L2 ^ (2 : ℕ) * l0 ^ (1 : ℕ) * lB ^ (2 : ℕ) : ℝ) +
    ((18 : ℝ) * B ^ (2 : ℕ) * r ^ (2 : ℕ) * q ^ (1 : ℕ) * L2 ^ (2 : ℕ) * l0 ^ (1 : ℕ) * lB ^ (1 : ℕ) : ℝ) +
    ((6 : ℝ) * B ^ (2 : ℕ) * r ^ (2 : ℕ) * q ^ (1 : ℕ) * L2 ^ (2 : ℕ) * lB ^ (3 : ℕ) : ℝ) +
    ((-9 : ℝ) * B ^ (2 : ℕ) * r ^ (2 : ℕ) * q ^ (1 : ℕ) * L2 ^ (2 : ℕ) * lB ^ (2 : ℕ) : ℝ) +
    ((-6 : ℝ) * B ^ (2 : ℕ) * r ^ (2 : ℕ) * V ^ (1 : ℕ) * L2 ^ (2 : ℕ) * l0 ^ (3 : ℕ) : ℝ) +
    ((15 : ℝ) * B ^ (2 : ℕ) * r ^ (2 : ℕ) * V ^ (1 : ℕ) * L2 ^ (2 : ℕ) * l0 ^ (2 : ℕ) * lB ^ (1 : ℕ) : ℝ)

noncomputable def pc11_9 : MvPolynomial (Fin 8) ℤ :=
  (MvPolynomial.C (-6 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (3 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (5 : Fin 8) ^ (2 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (-12 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (3 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (5 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (2 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (12 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (3 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (5 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (3 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (3 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (3 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (-6 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (3 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (2 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (-12 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (5 : Fin 8) ^ (3 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (30 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (5 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (-9 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (5 : Fin 8) ^ (2 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (-24 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (5 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (2 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (18 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (5 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (6 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (3 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (-9 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (2 : ℕ) : MvPolynomial (Fin 8) ℤ)

noncomputable def rc11_9 (B r q V L2 l0 l1 lB : ℝ) : ℝ :=
  ((-6 : ℝ) * B ^ (2 : ℕ) * r ^ (2 : ℕ) * V ^ (1 : ℕ) * L2 ^ (2 : ℕ) * l0 ^ (2 : ℕ) : ℝ) +
    ((-12 : ℝ) * B ^ (2 : ℕ) * r ^ (2 : ℕ) * V ^ (1 : ℕ) * L2 ^ (2 : ℕ) * l0 ^ (1 : ℕ) * lB ^ (2 : ℕ) : ℝ) +
    ((12 : ℝ) * B ^ (2 : ℕ) * r ^ (2 : ℕ) * V ^ (1 : ℕ) * L2 ^ (2 : ℕ) * l0 ^ (1 : ℕ) * lB ^ (1 : ℕ) : ℝ) +
    ((3 : ℝ) * B ^ (2 : ℕ) * r ^ (2 : ℕ) * V ^ (1 : ℕ) * L2 ^ (2 : ℕ) * lB ^ (3 : ℕ) : ℝ) +
    ((-6 : ℝ) * B ^ (2 : ℕ) * r ^ (2 : ℕ) * V ^ (1 : ℕ) * L2 ^ (2 : ℕ) * lB ^ (2 : ℕ) : ℝ) +
    ((-12 : ℝ) * B ^ (2 : ℕ) * r ^ (2 : ℕ) * L2 ^ (3 : ℕ) * l0 ^ (3 : ℕ) : ℝ) +
    ((30 : ℝ) * B ^ (2 : ℕ) * r ^ (2 : ℕ) * L2 ^ (3 : ℕ) * l0 ^ (2 : ℕ) * lB ^ (1 : ℕ) : ℝ) +
    ((-9 : ℝ) * B ^ (2 : ℕ) * r ^ (2 : ℕ) * L2 ^ (3 : ℕ) * l0 ^ (2 : ℕ) : ℝ) +
    ((-24 : ℝ) * B ^ (2 : ℕ) * r ^ (2 : ℕ) * L2 ^ (3 : ℕ) * l0 ^ (1 : ℕ) * lB ^ (2 : ℕ) : ℝ) +
    ((18 : ℝ) * B ^ (2 : ℕ) * r ^ (2 : ℕ) * L2 ^ (3 : ℕ) * l0 ^ (1 : ℕ) * lB ^ (1 : ℕ) : ℝ) +
    ((6 : ℝ) * B ^ (2 : ℕ) * r ^ (2 : ℕ) * L2 ^ (3 : ℕ) * lB ^ (3 : ℕ) : ℝ) +
    ((-9 : ℝ) * B ^ (2 : ℕ) * r ^ (2 : ℕ) * L2 ^ (3 : ℕ) * lB ^ (2 : ℕ) : ℝ)

noncomputable def pc11_10 : MvPolynomial (Fin 8) ℤ :=
  (MvPolynomial.C (18 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (5 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (6 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (-18 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (5 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (-18 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (5 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (6 : Fin 8) ^ (2 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (18 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (5 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (2 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (18 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (6 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (-18 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (6 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (2 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (3 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (5 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (6 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (-3 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (5 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (-15 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (5 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (6 : Fin 8) ^ (2 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (18 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (5 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (6 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (-12 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (5 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (6 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (-3 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (5 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (2 : ℕ) : MvPolynomial (Fin 8) ℤ)

noncomputable def rc11_10 (B r q V L2 l0 l1 lB : ℝ) : ℝ :=
  ((18 : ℝ) * B ^ (2 : ℕ) * r ^ (1 : ℕ) * q ^ (3 : ℕ) * l0 ^ (2 : ℕ) * l1 ^ (1 : ℕ) : ℝ) +
    ((-18 : ℝ) * B ^ (2 : ℕ) * r ^ (1 : ℕ) * q ^ (3 : ℕ) * l0 ^ (2 : ℕ) * lB ^ (1 : ℕ) : ℝ) +
    ((-18 : ℝ) * B ^ (2 : ℕ) * r ^ (1 : ℕ) * q ^ (3 : ℕ) * l0 ^ (1 : ℕ) * l1 ^ (2 : ℕ) : ℝ) +
    ((18 : ℝ) * B ^ (2 : ℕ) * r ^ (1 : ℕ) * q ^ (3 : ℕ) * l0 ^ (1 : ℕ) * lB ^ (2 : ℕ) : ℝ) +
    ((18 : ℝ) * B ^ (2 : ℕ) * r ^ (1 : ℕ) * q ^ (3 : ℕ) * l1 ^ (2 : ℕ) * lB ^ (1 : ℕ) : ℝ) +
    ((-18 : ℝ) * B ^ (2 : ℕ) * r ^ (1 : ℕ) * q ^ (3 : ℕ) * l1 ^ (1 : ℕ) * lB ^ (2 : ℕ) : ℝ) +
    ((3 : ℝ) * B ^ (2 : ℕ) * r ^ (1 : ℕ) * q ^ (2 : ℕ) * L2 ^ (1 : ℕ) * l0 ^ (2 : ℕ) * l1 ^ (1 : ℕ) : ℝ) +
    ((-3 : ℝ) * B ^ (2 : ℕ) * r ^ (1 : ℕ) * q ^ (2 : ℕ) * L2 ^ (1 : ℕ) * l0 ^ (2 : ℕ) * lB ^ (1 : ℕ) : ℝ) +
    ((-15 : ℝ) * B ^ (2 : ℕ) * r ^ (1 : ℕ) * q ^ (2 : ℕ) * L2 ^ (1 : ℕ) * l0 ^ (1 : ℕ) * l1 ^ (2 : ℕ) : ℝ) +
    ((18 : ℝ) * B ^ (2 : ℕ) * r ^ (1 : ℕ) * q ^ (2 : ℕ) * L2 ^ (1 : ℕ) * l0 ^ (1 : ℕ) * l1 ^ (1 : ℕ) * lB ^ (1 : ℕ) : ℝ) +
    ((-12 : ℝ) * B ^ (2 : ℕ) * r ^ (1 : ℕ) * q ^ (2 : ℕ) * L2 ^ (1 : ℕ) * l0 ^ (1 : ℕ) * l1 ^ (1 : ℕ) : ℝ) +
    ((-3 : ℝ) * B ^ (2 : ℕ) * r ^ (1 : ℕ) * q ^ (2 : ℕ) * L2 ^ (1 : ℕ) * l0 ^ (1 : ℕ) * lB ^ (2 : ℕ) : ℝ)

noncomputable def pc11_11 : MvPolynomial (Fin 8) ℤ :=
  (MvPolynomial.C (12 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (5 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (15 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (6 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (-21 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (6 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (2 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (12 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (6 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (6 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (3 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (-12 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (2 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (-18 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (3 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (5 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (6 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (18 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (3 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (5 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (18 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (3 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (5 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (6 : Fin 8) ^ (2 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (-18 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (3 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (5 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (2 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (-18 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (3 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (6 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (18 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (3 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (6 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (2 : ℕ) : MvPolynomial (Fin 8) ℤ)

noncomputable def rc11_11 (B r q V L2 l0 l1 lB : ℝ) : ℝ :=
  ((12 : ℝ) * B ^ (2 : ℕ) * r ^ (1 : ℕ) * q ^ (2 : ℕ) * L2 ^ (1 : ℕ) * l0 ^ (1 : ℕ) * lB ^ (1 : ℕ) : ℝ) +
    ((15 : ℝ) * B ^ (2 : ℕ) * r ^ (1 : ℕ) * q ^ (2 : ℕ) * L2 ^ (1 : ℕ) * l1 ^ (2 : ℕ) * lB ^ (1 : ℕ) : ℝ) +
    ((-21 : ℝ) * B ^ (2 : ℕ) * r ^ (1 : ℕ) * q ^ (2 : ℕ) * L2 ^ (1 : ℕ) * l1 ^ (1 : ℕ) * lB ^ (2 : ℕ) : ℝ) +
    ((12 : ℝ) * B ^ (2 : ℕ) * r ^ (1 : ℕ) * q ^ (2 : ℕ) * L2 ^ (1 : ℕ) * l1 ^ (1 : ℕ) * lB ^ (1 : ℕ) : ℝ) +
    ((6 : ℝ) * B ^ (2 : ℕ) * r ^ (1 : ℕ) * q ^ (2 : ℕ) * L2 ^ (1 : ℕ) * lB ^ (3 : ℕ) : ℝ) +
    ((-12 : ℝ) * B ^ (2 : ℕ) * r ^ (1 : ℕ) * q ^ (2 : ℕ) * L2 ^ (1 : ℕ) * lB ^ (2 : ℕ) : ℝ) +
    ((-18 : ℝ) * B ^ (2 : ℕ) * r ^ (1 : ℕ) * q ^ (1 : ℕ) * V ^ (1 : ℕ) * L2 ^ (1 : ℕ) * l0 ^ (2 : ℕ) * l1 ^ (1 : ℕ) : ℝ) +
    ((18 : ℝ) * B ^ (2 : ℕ) * r ^ (1 : ℕ) * q ^ (1 : ℕ) * V ^ (1 : ℕ) * L2 ^ (1 : ℕ) * l0 ^ (2 : ℕ) * lB ^ (1 : ℕ) : ℝ) +
    ((18 : ℝ) * B ^ (2 : ℕ) * r ^ (1 : ℕ) * q ^ (1 : ℕ) * V ^ (1 : ℕ) * L2 ^ (1 : ℕ) * l0 ^ (1 : ℕ) * l1 ^ (2 : ℕ) : ℝ) +
    ((-18 : ℝ) * B ^ (2 : ℕ) * r ^ (1 : ℕ) * q ^ (1 : ℕ) * V ^ (1 : ℕ) * L2 ^ (1 : ℕ) * l0 ^ (1 : ℕ) * lB ^ (2 : ℕ) : ℝ) +
    ((-18 : ℝ) * B ^ (2 : ℕ) * r ^ (1 : ℕ) * q ^ (1 : ℕ) * V ^ (1 : ℕ) * L2 ^ (1 : ℕ) * l1 ^ (2 : ℕ) * lB ^ (1 : ℕ) : ℝ) +
    ((18 : ℝ) * B ^ (2 : ℕ) * r ^ (1 : ℕ) * q ^ (1 : ℕ) * V ^ (1 : ℕ) * L2 ^ (1 : ℕ) * l1 ^ (1 : ℕ) * lB ^ (2 : ℕ) : ℝ)

noncomputable def pc11_12 : MvPolynomial (Fin 8) ℤ :=
  (MvPolynomial.C (-3 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (3 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (5 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (6 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (3 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (3 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (5 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (15 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (3 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (5 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (6 : Fin 8) ^ (2 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (-18 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (3 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (5 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (6 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (12 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (3 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (5 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (6 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (3 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (3 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (5 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (2 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (-12 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (3 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (5 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (-15 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (3 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (6 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (21 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (3 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (6 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (2 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (-12 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (3 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (6 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (-6 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (3 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (3 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (12 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (3 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (2 : ℕ) : MvPolynomial (Fin 8) ℤ)

noncomputable def rc11_12 (B r q V L2 l0 l1 lB : ℝ) : ℝ :=
  ((-3 : ℝ) * B ^ (2 : ℕ) * r ^ (1 : ℕ) * V ^ (1 : ℕ) * L2 ^ (2 : ℕ) * l0 ^ (2 : ℕ) * l1 ^ (1 : ℕ) : ℝ) +
    ((3 : ℝ) * B ^ (2 : ℕ) * r ^ (1 : ℕ) * V ^ (1 : ℕ) * L2 ^ (2 : ℕ) * l0 ^ (2 : ℕ) * lB ^ (1 : ℕ) : ℝ) +
    ((15 : ℝ) * B ^ (2 : ℕ) * r ^ (1 : ℕ) * V ^ (1 : ℕ) * L2 ^ (2 : ℕ) * l0 ^ (1 : ℕ) * l1 ^ (2 : ℕ) : ℝ) +
    ((-18 : ℝ) * B ^ (2 : ℕ) * r ^ (1 : ℕ) * V ^ (1 : ℕ) * L2 ^ (2 : ℕ) * l0 ^ (1 : ℕ) * l1 ^ (1 : ℕ) * lB ^ (1 : ℕ) : ℝ) +
    ((12 : ℝ) * B ^ (2 : ℕ) * r ^ (1 : ℕ) * V ^ (1 : ℕ) * L2 ^ (2 : ℕ) * l0 ^ (1 : ℕ) * l1 ^ (1 : ℕ) : ℝ) +
    ((3 : ℝ) * B ^ (2 : ℕ) * r ^ (1 : ℕ) * V ^ (1 : ℕ) * L2 ^ (2 : ℕ) * l0 ^ (1 : ℕ) * lB ^ (2 : ℕ) : ℝ) +
    ((-12 : ℝ) * B ^ (2 : ℕ) * r ^ (1 : ℕ) * V ^ (1 : ℕ) * L2 ^ (2 : ℕ) * l0 ^ (1 : ℕ) * lB ^ (1 : ℕ) : ℝ) +
    ((-15 : ℝ) * B ^ (2 : ℕ) * r ^ (1 : ℕ) * V ^ (1 : ℕ) * L2 ^ (2 : ℕ) * l1 ^ (2 : ℕ) * lB ^ (1 : ℕ) : ℝ) +
    ((21 : ℝ) * B ^ (2 : ℕ) * r ^ (1 : ℕ) * V ^ (1 : ℕ) * L2 ^ (2 : ℕ) * l1 ^ (1 : ℕ) * lB ^ (2 : ℕ) : ℝ) +
    ((-12 : ℝ) * B ^ (2 : ℕ) * r ^ (1 : ℕ) * V ^ (1 : ℕ) * L2 ^ (2 : ℕ) * l1 ^ (1 : ℕ) * lB ^ (1 : ℕ) : ℝ) +
    ((-6 : ℝ) * B ^ (2 : ℕ) * r ^ (1 : ℕ) * V ^ (1 : ℕ) * L2 ^ (2 : ℕ) * lB ^ (3 : ℕ) : ℝ) +
    ((12 : ℝ) * B ^ (2 : ℕ) * r ^ (1 : ℕ) * V ^ (1 : ℕ) * L2 ^ (2 : ℕ) * lB ^ (2 : ℕ) : ℝ)

noncomputable def pc11_13 : MvPolynomial (Fin 8) ℤ :=
  (MvPolynomial.C (-9 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (5 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (6 : Fin 8) ^ (2 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (18 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (5 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (6 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (-9 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (5 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (2 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (9 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (6 : Fin 8) ^ (3 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (-18 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (6 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (9 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (6 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (2 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (-9 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (5 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (6 : Fin 8) ^ (2 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (18 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (5 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (6 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (-9 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (5 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (2 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (15 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (6 : Fin 8) ^ (3 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (-33 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (6 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (6 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (6 : Fin 8) ^ (2 : ℕ) : MvPolynomial (Fin 8) ℤ)

noncomputable def rc11_13 (B r q V L2 l0 l1 lB : ℝ) : ℝ :=
  ((-9 : ℝ) * B ^ (2 : ℕ) * q ^ (3 : ℕ) * l0 ^ (1 : ℕ) * l1 ^ (2 : ℕ) : ℝ) +
    ((18 : ℝ) * B ^ (2 : ℕ) * q ^ (3 : ℕ) * l0 ^ (1 : ℕ) * l1 ^ (1 : ℕ) * lB ^ (1 : ℕ) : ℝ) +
    ((-9 : ℝ) * B ^ (2 : ℕ) * q ^ (3 : ℕ) * l0 ^ (1 : ℕ) * lB ^ (2 : ℕ) : ℝ) +
    ((9 : ℝ) * B ^ (2 : ℕ) * q ^ (3 : ℕ) * l1 ^ (3 : ℕ) : ℝ) +
    ((-18 : ℝ) * B ^ (2 : ℕ) * q ^ (3 : ℕ) * l1 ^ (2 : ℕ) * lB ^ (1 : ℕ) : ℝ) +
    ((9 : ℝ) * B ^ (2 : ℕ) * q ^ (3 : ℕ) * l1 ^ (1 : ℕ) * lB ^ (2 : ℕ) : ℝ) +
    ((-9 : ℝ) * B ^ (2 : ℕ) * q ^ (2 : ℕ) * L2 ^ (1 : ℕ) * l0 ^ (1 : ℕ) * l1 ^ (2 : ℕ) : ℝ) +
    ((18 : ℝ) * B ^ (2 : ℕ) * q ^ (2 : ℕ) * L2 ^ (1 : ℕ) * l0 ^ (1 : ℕ) * l1 ^ (1 : ℕ) * lB ^ (1 : ℕ) : ℝ) +
    ((-9 : ℝ) * B ^ (2 : ℕ) * q ^ (2 : ℕ) * L2 ^ (1 : ℕ) * l0 ^ (1 : ℕ) * lB ^ (2 : ℕ) : ℝ) +
    ((15 : ℝ) * B ^ (2 : ℕ) * q ^ (2 : ℕ) * L2 ^ (1 : ℕ) * l1 ^ (3 : ℕ) : ℝ) +
    ((-33 : ℝ) * B ^ (2 : ℕ) * q ^ (2 : ℕ) * L2 ^ (1 : ℕ) * l1 ^ (2 : ℕ) * lB ^ (1 : ℕ) : ℝ) +
    ((6 : ℝ) * B ^ (2 : ℕ) * q ^ (2 : ℕ) * L2 ^ (1 : ℕ) * l1 ^ (2 : ℕ) : ℝ)

noncomputable def pc11_14 : MvPolynomial (Fin 8) ℤ :=
  (MvPolynomial.C (21 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (6 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (2 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (-12 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (6 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (-3 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (3 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (6 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (2 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (9 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (3 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (5 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (6 : Fin 8) ^ (2 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (-18 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (3 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (5 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (6 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (9 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (3 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (5 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (2 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (-9 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (3 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (6 : Fin 8) ^ (3 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (18 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (3 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (6 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (-9 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (3 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (6 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (2 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (12 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (6 : Fin 8) ^ (3 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (-30 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (6 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ)

noncomputable def rc11_14 (B r q V L2 l0 l1 lB : ℝ) : ℝ :=
  ((21 : ℝ) * B ^ (2 : ℕ) * q ^ (2 : ℕ) * L2 ^ (1 : ℕ) * l1 ^ (1 : ℕ) * lB ^ (2 : ℕ) : ℝ) +
    ((-12 : ℝ) * B ^ (2 : ℕ) * q ^ (2 : ℕ) * L2 ^ (1 : ℕ) * l1 ^ (1 : ℕ) * lB ^ (1 : ℕ) : ℝ) +
    ((-3 : ℝ) * B ^ (2 : ℕ) * q ^ (2 : ℕ) * L2 ^ (1 : ℕ) * lB ^ (3 : ℕ) : ℝ) +
    ((6 : ℝ) * B ^ (2 : ℕ) * q ^ (2 : ℕ) * L2 ^ (1 : ℕ) * lB ^ (2 : ℕ) : ℝ) +
    ((9 : ℝ) * B ^ (2 : ℕ) * q ^ (1 : ℕ) * V ^ (1 : ℕ) * L2 ^ (1 : ℕ) * l0 ^ (1 : ℕ) * l1 ^ (2 : ℕ) : ℝ) +
    ((-18 : ℝ) * B ^ (2 : ℕ) * q ^ (1 : ℕ) * V ^ (1 : ℕ) * L2 ^ (1 : ℕ) * l0 ^ (1 : ℕ) * l1 ^ (1 : ℕ) * lB ^ (1 : ℕ) : ℝ) +
    ((9 : ℝ) * B ^ (2 : ℕ) * q ^ (1 : ℕ) * V ^ (1 : ℕ) * L2 ^ (1 : ℕ) * l0 ^ (1 : ℕ) * lB ^ (2 : ℕ) : ℝ) +
    ((-9 : ℝ) * B ^ (2 : ℕ) * q ^ (1 : ℕ) * V ^ (1 : ℕ) * L2 ^ (1 : ℕ) * l1 ^ (3 : ℕ) : ℝ) +
    ((18 : ℝ) * B ^ (2 : ℕ) * q ^ (1 : ℕ) * V ^ (1 : ℕ) * L2 ^ (1 : ℕ) * l1 ^ (2 : ℕ) * lB ^ (1 : ℕ) : ℝ) +
    ((-9 : ℝ) * B ^ (2 : ℕ) * q ^ (1 : ℕ) * V ^ (1 : ℕ) * L2 ^ (1 : ℕ) * l1 ^ (1 : ℕ) * lB ^ (2 : ℕ) : ℝ) +
    ((12 : ℝ) * B ^ (2 : ℕ) * q ^ (1 : ℕ) * L2 ^ (2 : ℕ) * l1 ^ (3 : ℕ) : ℝ) +
    ((-30 : ℝ) * B ^ (2 : ℕ) * q ^ (1 : ℕ) * L2 ^ (2 : ℕ) * l1 ^ (2 : ℕ) * lB ^ (1 : ℕ) : ℝ)

noncomputable def pc11_15 : MvPolynomial (Fin 8) ℤ :=
  (MvPolynomial.C (9 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (6 : Fin 8) ^ (2 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (24 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (6 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (2 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (-18 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (6 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (-6 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (3 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (9 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (2 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (9 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (3 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (5 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (6 : Fin 8) ^ (2 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (-18 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (3 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (5 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (6 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (9 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (3 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (5 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (2 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (-15 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (3 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (6 : Fin 8) ^ (3 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (33 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (3 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (6 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (-6 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (3 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (6 : Fin 8) ^ (2 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (-21 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (3 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (6 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (2 : ℕ) : MvPolynomial (Fin 8) ℤ)

noncomputable def rc11_15 (B r q V L2 l0 l1 lB : ℝ) : ℝ :=
  ((9 : ℝ) * B ^ (2 : ℕ) * q ^ (1 : ℕ) * L2 ^ (2 : ℕ) * l1 ^ (2 : ℕ) : ℝ) +
    ((24 : ℝ) * B ^ (2 : ℕ) * q ^ (1 : ℕ) * L2 ^ (2 : ℕ) * l1 ^ (1 : ℕ) * lB ^ (2 : ℕ) : ℝ) +
    ((-18 : ℝ) * B ^ (2 : ℕ) * q ^ (1 : ℕ) * L2 ^ (2 : ℕ) * l1 ^ (1 : ℕ) * lB ^ (1 : ℕ) : ℝ) +
    ((-6 : ℝ) * B ^ (2 : ℕ) * q ^ (1 : ℕ) * L2 ^ (2 : ℕ) * lB ^ (3 : ℕ) : ℝ) +
    ((9 : ℝ) * B ^ (2 : ℕ) * q ^ (1 : ℕ) * L2 ^ (2 : ℕ) * lB ^ (2 : ℕ) : ℝ) +
    ((9 : ℝ) * B ^ (2 : ℕ) * V ^ (1 : ℕ) * L2 ^ (2 : ℕ) * l0 ^ (1 : ℕ) * l1 ^ (2 : ℕ) : ℝ) +
    ((-18 : ℝ) * B ^ (2 : ℕ) * V ^ (1 : ℕ) * L2 ^ (2 : ℕ) * l0 ^ (1 : ℕ) * l1 ^ (1 : ℕ) * lB ^ (1 : ℕ) : ℝ) +
    ((9 : ℝ) * B ^ (2 : ℕ) * V ^ (1 : ℕ) * L2 ^ (2 : ℕ) * l0 ^ (1 : ℕ) * lB ^ (2 : ℕ) : ℝ) +
    ((-15 : ℝ) * B ^ (2 : ℕ) * V ^ (1 : ℕ) * L2 ^ (2 : ℕ) * l1 ^ (3 : ℕ) : ℝ) +
    ((33 : ℝ) * B ^ (2 : ℕ) * V ^ (1 : ℕ) * L2 ^ (2 : ℕ) * l1 ^ (2 : ℕ) * lB ^ (1 : ℕ) : ℝ) +
    ((-6 : ℝ) * B ^ (2 : ℕ) * V ^ (1 : ℕ) * L2 ^ (2 : ℕ) * l1 ^ (2 : ℕ) : ℝ) +
    ((-21 : ℝ) * B ^ (2 : ℕ) * V ^ (1 : ℕ) * L2 ^ (2 : ℕ) * l1 ^ (1 : ℕ) * lB ^ (2 : ℕ) : ℝ)

noncomputable def pc11_16 : MvPolynomial (Fin 8) ℤ :=
  (MvPolynomial.C (12 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (3 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (6 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (3 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (3 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (3 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (-6 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (3 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (2 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (9 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (5 : Fin 8) ^ (3 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (-18 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (5 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (6 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (-9 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (5 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (9 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (5 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (6 : Fin 8) ^ (2 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (18 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (5 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (6 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (-9 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (6 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (-3 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (5 : Fin 8) ^ (3 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (-3 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (5 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (6 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (6 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (5 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ)

noncomputable def rc11_16 (B r q V L2 l0 l1 lB : ℝ) : ℝ :=
  ((12 : ℝ) * B ^ (2 : ℕ) * V ^ (1 : ℕ) * L2 ^ (2 : ℕ) * l1 ^ (1 : ℕ) * lB ^ (1 : ℕ) : ℝ) +
    ((3 : ℝ) * B ^ (2 : ℕ) * V ^ (1 : ℕ) * L2 ^ (2 : ℕ) * lB ^ (3 : ℕ) : ℝ) +
    ((-6 : ℝ) * B ^ (2 : ℕ) * V ^ (1 : ℕ) * L2 ^ (2 : ℕ) * lB ^ (2 : ℕ) : ℝ) +
    ((9 : ℝ) * B ^ (1 : ℕ) * r ^ (1 : ℕ) * q ^ (3 : ℕ) * l0 ^ (3 : ℕ) : ℝ) +
    ((-18 : ℝ) * B ^ (1 : ℕ) * r ^ (1 : ℕ) * q ^ (3 : ℕ) * l0 ^ (2 : ℕ) * l1 ^ (1 : ℕ) : ℝ) +
    ((-9 : ℝ) * B ^ (1 : ℕ) * r ^ (1 : ℕ) * q ^ (3 : ℕ) * l0 ^ (2 : ℕ) * lB ^ (1 : ℕ) : ℝ) +
    ((9 : ℝ) * B ^ (1 : ℕ) * r ^ (1 : ℕ) * q ^ (3 : ℕ) * l0 ^ (1 : ℕ) * l1 ^ (2 : ℕ) : ℝ) +
    ((18 : ℝ) * B ^ (1 : ℕ) * r ^ (1 : ℕ) * q ^ (3 : ℕ) * l0 ^ (1 : ℕ) * l1 ^ (1 : ℕ) * lB ^ (1 : ℕ) : ℝ) +
    ((-9 : ℝ) * B ^ (1 : ℕ) * r ^ (1 : ℕ) * q ^ (3 : ℕ) * l1 ^ (2 : ℕ) * lB ^ (1 : ℕ) : ℝ) +
    ((-3 : ℝ) * B ^ (1 : ℕ) * r ^ (1 : ℕ) * q ^ (2 : ℕ) * L2 ^ (1 : ℕ) * l0 ^ (3 : ℕ) : ℝ) +
    ((-3 : ℝ) * B ^ (1 : ℕ) * r ^ (1 : ℕ) * q ^ (2 : ℕ) * L2 ^ (1 : ℕ) * l0 ^ (2 : ℕ) * l1 ^ (1 : ℕ) : ℝ) +
    ((6 : ℝ) * B ^ (1 : ℕ) * r ^ (1 : ℕ) * q ^ (2 : ℕ) * L2 ^ (1 : ℕ) * l0 ^ (2 : ℕ) * lB ^ (1 : ℕ) : ℝ)

noncomputable def pc11_17 : MvPolynomial (Fin 8) ℤ :=
  (MvPolynomial.C (-12 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (5 : Fin 8) ^ (2 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (6 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (5 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (6 : Fin 8) ^ (2 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (12 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (5 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (6 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (-3 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (5 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (2 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (12 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (5 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (-6 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (6 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (3 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (6 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (2 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (-12 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (6 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (-9 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (3 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (5 : Fin 8) ^ (3 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (18 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (3 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (5 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (6 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (9 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (3 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (5 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (-9 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (3 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (5 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (6 : Fin 8) ^ (2 : ℕ) : MvPolynomial (Fin 8) ℤ)

noncomputable def rc11_17 (B r q V L2 l0 l1 lB : ℝ) : ℝ :=
  ((-12 : ℝ) * B ^ (1 : ℕ) * r ^ (1 : ℕ) * q ^ (2 : ℕ) * L2 ^ (1 : ℕ) * l0 ^ (2 : ℕ) : ℝ) +
    ((6 : ℝ) * B ^ (1 : ℕ) * r ^ (1 : ℕ) * q ^ (2 : ℕ) * L2 ^ (1 : ℕ) * l0 ^ (1 : ℕ) * l1 ^ (2 : ℕ) : ℝ) +
    ((12 : ℝ) * B ^ (1 : ℕ) * r ^ (1 : ℕ) * q ^ (2 : ℕ) * L2 ^ (1 : ℕ) * l0 ^ (1 : ℕ) * l1 ^ (1 : ℕ) : ℝ) +
    ((-3 : ℝ) * B ^ (1 : ℕ) * r ^ (1 : ℕ) * q ^ (2 : ℕ) * L2 ^ (1 : ℕ) * l0 ^ (1 : ℕ) * lB ^ (2 : ℕ) : ℝ) +
    ((12 : ℝ) * B ^ (1 : ℕ) * r ^ (1 : ℕ) * q ^ (2 : ℕ) * L2 ^ (1 : ℕ) * l0 ^ (1 : ℕ) * lB ^ (1 : ℕ) : ℝ) +
    ((-6 : ℝ) * B ^ (1 : ℕ) * r ^ (1 : ℕ) * q ^ (2 : ℕ) * L2 ^ (1 : ℕ) * l1 ^ (2 : ℕ) * lB ^ (1 : ℕ) : ℝ) +
    ((3 : ℝ) * B ^ (1 : ℕ) * r ^ (1 : ℕ) * q ^ (2 : ℕ) * L2 ^ (1 : ℕ) * l1 ^ (1 : ℕ) * lB ^ (2 : ℕ) : ℝ) +
    ((-12 : ℝ) * B ^ (1 : ℕ) * r ^ (1 : ℕ) * q ^ (2 : ℕ) * L2 ^ (1 : ℕ) * l1 ^ (1 : ℕ) * lB ^ (1 : ℕ) : ℝ) +
    ((-9 : ℝ) * B ^ (1 : ℕ) * r ^ (1 : ℕ) * q ^ (1 : ℕ) * V ^ (1 : ℕ) * L2 ^ (1 : ℕ) * l0 ^ (3 : ℕ) : ℝ) +
    ((18 : ℝ) * B ^ (1 : ℕ) * r ^ (1 : ℕ) * q ^ (1 : ℕ) * V ^ (1 : ℕ) * L2 ^ (1 : ℕ) * l0 ^ (2 : ℕ) * l1 ^ (1 : ℕ) : ℝ) +
    ((9 : ℝ) * B ^ (1 : ℕ) * r ^ (1 : ℕ) * q ^ (1 : ℕ) * V ^ (1 : ℕ) * L2 ^ (1 : ℕ) * l0 ^ (2 : ℕ) * lB ^ (1 : ℕ) : ℝ) +
    ((-9 : ℝ) * B ^ (1 : ℕ) * r ^ (1 : ℕ) * q ^ (1 : ℕ) * V ^ (1 : ℕ) * L2 ^ (1 : ℕ) * l0 ^ (1 : ℕ) * l1 ^ (2 : ℕ) : ℝ)

noncomputable def pc11_18 : MvPolynomial (Fin 8) ℤ :=
  (MvPolynomial.C (-18 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (3 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (5 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (6 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (9 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (3 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (6 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (7 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (5 : Fin 8) ^ (3 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (-12 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (5 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (15 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (5 : Fin 8) ^ (2 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (6 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (5 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (2 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (-18 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (5 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (6 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (5 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (-1 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (3 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (3 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (2 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (-6 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (3 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (3 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (5 : Fin 8) ^ (3 : ℕ) : MvPolynomial (Fin 8) ℤ)

noncomputable def rc11_18 (B r q V L2 l0 l1 lB : ℝ) : ℝ :=
  ((-18 : ℝ) * B ^ (1 : ℕ) * r ^ (1 : ℕ) * q ^ (1 : ℕ) * V ^ (1 : ℕ) * L2 ^ (1 : ℕ) * l0 ^ (1 : ℕ) * l1 ^ (1 : ℕ) * lB ^ (1 : ℕ) : ℝ) +
    ((9 : ℝ) * B ^ (1 : ℕ) * r ^ (1 : ℕ) * q ^ (1 : ℕ) * V ^ (1 : ℕ) * L2 ^ (1 : ℕ) * l1 ^ (2 : ℕ) * lB ^ (1 : ℕ) : ℝ) +
    ((7 : ℝ) * B ^ (1 : ℕ) * r ^ (1 : ℕ) * q ^ (1 : ℕ) * L2 ^ (2 : ℕ) * l0 ^ (3 : ℕ) : ℝ) +
    ((-12 : ℝ) * B ^ (1 : ℕ) * r ^ (1 : ℕ) * q ^ (1 : ℕ) * L2 ^ (2 : ℕ) * l0 ^ (2 : ℕ) * lB ^ (1 : ℕ) : ℝ) +
    ((15 : ℝ) * B ^ (1 : ℕ) * r ^ (1 : ℕ) * q ^ (1 : ℕ) * L2 ^ (2 : ℕ) * l0 ^ (2 : ℕ) : ℝ) +
    ((6 : ℝ) * B ^ (1 : ℕ) * r ^ (1 : ℕ) * q ^ (1 : ℕ) * L2 ^ (2 : ℕ) * l0 ^ (1 : ℕ) * lB ^ (2 : ℕ) : ℝ) +
    ((-18 : ℝ) * B ^ (1 : ℕ) * r ^ (1 : ℕ) * q ^ (1 : ℕ) * L2 ^ (2 : ℕ) * l0 ^ (1 : ℕ) * lB ^ (1 : ℕ) : ℝ) +
    ((6 : ℝ) * B ^ (1 : ℕ) * r ^ (1 : ℕ) * q ^ (1 : ℕ) * L2 ^ (2 : ℕ) * l0 ^ (1 : ℕ) : ℝ) +
    ((-1 : ℝ) * B ^ (1 : ℕ) * r ^ (1 : ℕ) * q ^ (1 : ℕ) * L2 ^ (2 : ℕ) * lB ^ (3 : ℕ) : ℝ) +
    ((3 : ℝ) * B ^ (1 : ℕ) * r ^ (1 : ℕ) * q ^ (1 : ℕ) * L2 ^ (2 : ℕ) * lB ^ (2 : ℕ) : ℝ) +
    ((-6 : ℝ) * B ^ (1 : ℕ) * r ^ (1 : ℕ) * q ^ (1 : ℕ) * L2 ^ (2 : ℕ) * lB ^ (1 : ℕ) : ℝ) +
    ((3 : ℝ) * B ^ (1 : ℕ) * r ^ (1 : ℕ) * V ^ (1 : ℕ) * L2 ^ (2 : ℕ) * l0 ^ (3 : ℕ) : ℝ)

noncomputable def pc11_19 : MvPolynomial (Fin 8) ℤ :=
  (MvPolynomial.C (3 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (3 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (5 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (6 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (-6 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (3 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (5 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (12 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (3 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (5 : Fin 8) ^ (2 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (-6 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (3 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (5 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (6 : Fin 8) ^ (2 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (-12 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (3 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (5 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (6 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (3 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (3 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (5 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (2 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (-12 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (3 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (5 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (6 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (3 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (6 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (-3 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (3 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (6 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (2 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (12 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (3 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (6 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (7 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (5 : Fin 8) ^ (3 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (-12 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (5 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ)

noncomputable def rc11_19 (B r q V L2 l0 l1 lB : ℝ) : ℝ :=
  ((3 : ℝ) * B ^ (1 : ℕ) * r ^ (1 : ℕ) * V ^ (1 : ℕ) * L2 ^ (2 : ℕ) * l0 ^ (2 : ℕ) * l1 ^ (1 : ℕ) : ℝ) +
    ((-6 : ℝ) * B ^ (1 : ℕ) * r ^ (1 : ℕ) * V ^ (1 : ℕ) * L2 ^ (2 : ℕ) * l0 ^ (2 : ℕ) * lB ^ (1 : ℕ) : ℝ) +
    ((12 : ℝ) * B ^ (1 : ℕ) * r ^ (1 : ℕ) * V ^ (1 : ℕ) * L2 ^ (2 : ℕ) * l0 ^ (2 : ℕ) : ℝ) +
    ((-6 : ℝ) * B ^ (1 : ℕ) * r ^ (1 : ℕ) * V ^ (1 : ℕ) * L2 ^ (2 : ℕ) * l0 ^ (1 : ℕ) * l1 ^ (2 : ℕ) : ℝ) +
    ((-12 : ℝ) * B ^ (1 : ℕ) * r ^ (1 : ℕ) * V ^ (1 : ℕ) * L2 ^ (2 : ℕ) * l0 ^ (1 : ℕ) * l1 ^ (1 : ℕ) : ℝ) +
    ((3 : ℝ) * B ^ (1 : ℕ) * r ^ (1 : ℕ) * V ^ (1 : ℕ) * L2 ^ (2 : ℕ) * l0 ^ (1 : ℕ) * lB ^ (2 : ℕ) : ℝ) +
    ((-12 : ℝ) * B ^ (1 : ℕ) * r ^ (1 : ℕ) * V ^ (1 : ℕ) * L2 ^ (2 : ℕ) * l0 ^ (1 : ℕ) * lB ^ (1 : ℕ) : ℝ) +
    ((6 : ℝ) * B ^ (1 : ℕ) * r ^ (1 : ℕ) * V ^ (1 : ℕ) * L2 ^ (2 : ℕ) * l1 ^ (2 : ℕ) * lB ^ (1 : ℕ) : ℝ) +
    ((-3 : ℝ) * B ^ (1 : ℕ) * r ^ (1 : ℕ) * V ^ (1 : ℕ) * L2 ^ (2 : ℕ) * l1 ^ (1 : ℕ) * lB ^ (2 : ℕ) : ℝ) +
    ((12 : ℝ) * B ^ (1 : ℕ) * r ^ (1 : ℕ) * V ^ (1 : ℕ) * L2 ^ (2 : ℕ) * l1 ^ (1 : ℕ) * lB ^ (1 : ℕ) : ℝ) +
    ((7 : ℝ) * B ^ (1 : ℕ) * r ^ (1 : ℕ) * L2 ^ (3 : ℕ) * l0 ^ (3 : ℕ) : ℝ) +
    ((-12 : ℝ) * B ^ (1 : ℕ) * r ^ (1 : ℕ) * L2 ^ (3 : ℕ) * l0 ^ (2 : ℕ) * lB ^ (1 : ℕ) : ℝ)

noncomputable def pc11_20 : MvPolynomial (Fin 8) ℤ :=
  (MvPolynomial.C (15 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (5 : Fin 8) ^ (2 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (6 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (5 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (2 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (-18 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (5 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (6 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (5 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (-1 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (3 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (3 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (2 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (-6 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (1 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (-9 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (5 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (6 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (9 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (5 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (18 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (5 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (6 : Fin 8) ^ (2 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (-18 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (5 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (6 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (-9 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (6 : Fin 8) ^ (3 : ℕ) : MvPolynomial (Fin 8) ℤ)

noncomputable def rc11_20 (B r q V L2 l0 l1 lB : ℝ) : ℝ :=
  ((15 : ℝ) * B ^ (1 : ℕ) * r ^ (1 : ℕ) * L2 ^ (3 : ℕ) * l0 ^ (2 : ℕ) : ℝ) +
    ((6 : ℝ) * B ^ (1 : ℕ) * r ^ (1 : ℕ) * L2 ^ (3 : ℕ) * l0 ^ (1 : ℕ) * lB ^ (2 : ℕ) : ℝ) +
    ((-18 : ℝ) * B ^ (1 : ℕ) * r ^ (1 : ℕ) * L2 ^ (3 : ℕ) * l0 ^ (1 : ℕ) * lB ^ (1 : ℕ) : ℝ) +
    ((6 : ℝ) * B ^ (1 : ℕ) * r ^ (1 : ℕ) * L2 ^ (3 : ℕ) * l0 ^ (1 : ℕ) : ℝ) +
    ((-1 : ℝ) * B ^ (1 : ℕ) * r ^ (1 : ℕ) * L2 ^ (3 : ℕ) * lB ^ (3 : ℕ) : ℝ) +
    ((3 : ℝ) * B ^ (1 : ℕ) * r ^ (1 : ℕ) * L2 ^ (3 : ℕ) * lB ^ (2 : ℕ) : ℝ) +
    ((-6 : ℝ) * B ^ (1 : ℕ) * r ^ (1 : ℕ) * L2 ^ (3 : ℕ) * lB ^ (1 : ℕ) : ℝ) +
    ((-9 : ℝ) * B ^ (1 : ℕ) * q ^ (3 : ℕ) * l0 ^ (2 : ℕ) * l1 ^ (1 : ℕ) : ℝ) +
    ((9 : ℝ) * B ^ (1 : ℕ) * q ^ (3 : ℕ) * l0 ^ (2 : ℕ) * lB ^ (1 : ℕ) : ℝ) +
    ((18 : ℝ) * B ^ (1 : ℕ) * q ^ (3 : ℕ) * l0 ^ (1 : ℕ) * l1 ^ (2 : ℕ) : ℝ) +
    ((-18 : ℝ) * B ^ (1 : ℕ) * q ^ (3 : ℕ) * l0 ^ (1 : ℕ) * l1 ^ (1 : ℕ) * lB ^ (1 : ℕ) : ℝ) +
    ((-9 : ℝ) * B ^ (1 : ℕ) * q ^ (3 : ℕ) * l1 ^ (3 : ℕ) : ℝ)

noncomputable def pc11_21 : MvPolynomial (Fin 8) ℤ :=
  (MvPolynomial.C (9 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (6 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (-3 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (5 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (6 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (3 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (5 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (15 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (5 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (6 : Fin 8) ^ (2 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (-18 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (5 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (6 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (12 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (5 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (6 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (3 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (5 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (2 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (-12 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (5 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (-12 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (6 : Fin 8) ^ (3 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (15 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (6 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (-12 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (6 : Fin 8) ^ (2 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (-3 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (6 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (2 : ℕ) : MvPolynomial (Fin 8) ℤ)

noncomputable def rc11_21 (B r q V L2 l0 l1 lB : ℝ) : ℝ :=
  ((9 : ℝ) * B ^ (1 : ℕ) * q ^ (3 : ℕ) * l1 ^ (2 : ℕ) * lB ^ (1 : ℕ) : ℝ) +
    ((-3 : ℝ) * B ^ (1 : ℕ) * q ^ (2 : ℕ) * L2 ^ (1 : ℕ) * l0 ^ (2 : ℕ) * l1 ^ (1 : ℕ) : ℝ) +
    ((3 : ℝ) * B ^ (1 : ℕ) * q ^ (2 : ℕ) * L2 ^ (1 : ℕ) * l0 ^ (2 : ℕ) * lB ^ (1 : ℕ) : ℝ) +
    ((15 : ℝ) * B ^ (1 : ℕ) * q ^ (2 : ℕ) * L2 ^ (1 : ℕ) * l0 ^ (1 : ℕ) * l1 ^ (2 : ℕ) : ℝ) +
    ((-18 : ℝ) * B ^ (1 : ℕ) * q ^ (2 : ℕ) * L2 ^ (1 : ℕ) * l0 ^ (1 : ℕ) * l1 ^ (1 : ℕ) * lB ^ (1 : ℕ) : ℝ) +
    ((12 : ℝ) * B ^ (1 : ℕ) * q ^ (2 : ℕ) * L2 ^ (1 : ℕ) * l0 ^ (1 : ℕ) * l1 ^ (1 : ℕ) : ℝ) +
    ((3 : ℝ) * B ^ (1 : ℕ) * q ^ (2 : ℕ) * L2 ^ (1 : ℕ) * l0 ^ (1 : ℕ) * lB ^ (2 : ℕ) : ℝ) +
    ((-12 : ℝ) * B ^ (1 : ℕ) * q ^ (2 : ℕ) * L2 ^ (1 : ℕ) * l0 ^ (1 : ℕ) * lB ^ (1 : ℕ) : ℝ) +
    ((-12 : ℝ) * B ^ (1 : ℕ) * q ^ (2 : ℕ) * L2 ^ (1 : ℕ) * l1 ^ (3 : ℕ) : ℝ) +
    ((15 : ℝ) * B ^ (1 : ℕ) * q ^ (2 : ℕ) * L2 ^ (1 : ℕ) * l1 ^ (2 : ℕ) * lB ^ (1 : ℕ) : ℝ) +
    ((-12 : ℝ) * B ^ (1 : ℕ) * q ^ (2 : ℕ) * L2 ^ (1 : ℕ) * l1 ^ (2 : ℕ) : ℝ) +
    ((-3 : ℝ) * B ^ (1 : ℕ) * q ^ (2 : ℕ) * L2 ^ (1 : ℕ) * l1 ^ (1 : ℕ) * lB ^ (2 : ℕ) : ℝ)

noncomputable def pc11_22 : MvPolynomial (Fin 8) ℤ :=
  (MvPolynomial.C (12 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (6 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (9 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (3 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (5 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (6 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (-9 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (3 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (5 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (-18 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (3 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (5 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (6 : Fin 8) ^ (2 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (18 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (3 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (5 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (6 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (9 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (3 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (6 : Fin 8) ^ (3 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (-9 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (3 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (6 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (-7 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (6 : Fin 8) ^ (3 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (12 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (6 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (-15 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (6 : Fin 8) ^ (2 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (-6 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (6 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (2 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (18 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (6 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ)

noncomputable def rc11_22 (B r q V L2 l0 l1 lB : ℝ) : ℝ :=
  ((12 : ℝ) * B ^ (1 : ℕ) * q ^ (2 : ℕ) * L2 ^ (1 : ℕ) * l1 ^ (1 : ℕ) * lB ^ (1 : ℕ) : ℝ) +
    ((9 : ℝ) * B ^ (1 : ℕ) * q ^ (1 : ℕ) * V ^ (1 : ℕ) * L2 ^ (1 : ℕ) * l0 ^ (2 : ℕ) * l1 ^ (1 : ℕ) : ℝ) +
    ((-9 : ℝ) * B ^ (1 : ℕ) * q ^ (1 : ℕ) * V ^ (1 : ℕ) * L2 ^ (1 : ℕ) * l0 ^ (2 : ℕ) * lB ^ (1 : ℕ) : ℝ) +
    ((-18 : ℝ) * B ^ (1 : ℕ) * q ^ (1 : ℕ) * V ^ (1 : ℕ) * L2 ^ (1 : ℕ) * l0 ^ (1 : ℕ) * l1 ^ (2 : ℕ) : ℝ) +
    ((18 : ℝ) * B ^ (1 : ℕ) * q ^ (1 : ℕ) * V ^ (1 : ℕ) * L2 ^ (1 : ℕ) * l0 ^ (1 : ℕ) * l1 ^ (1 : ℕ) * lB ^ (1 : ℕ) : ℝ) +
    ((9 : ℝ) * B ^ (1 : ℕ) * q ^ (1 : ℕ) * V ^ (1 : ℕ) * L2 ^ (1 : ℕ) * l1 ^ (3 : ℕ) : ℝ) +
    ((-9 : ℝ) * B ^ (1 : ℕ) * q ^ (1 : ℕ) * V ^ (1 : ℕ) * L2 ^ (1 : ℕ) * l1 ^ (2 : ℕ) * lB ^ (1 : ℕ) : ℝ) +
    ((-7 : ℝ) * B ^ (1 : ℕ) * q ^ (1 : ℕ) * L2 ^ (2 : ℕ) * l1 ^ (3 : ℕ) : ℝ) +
    ((12 : ℝ) * B ^ (1 : ℕ) * q ^ (1 : ℕ) * L2 ^ (2 : ℕ) * l1 ^ (2 : ℕ) * lB ^ (1 : ℕ) : ℝ) +
    ((-15 : ℝ) * B ^ (1 : ℕ) * q ^ (1 : ℕ) * L2 ^ (2 : ℕ) * l1 ^ (2 : ℕ) : ℝ) +
    ((-6 : ℝ) * B ^ (1 : ℕ) * q ^ (1 : ℕ) * L2 ^ (2 : ℕ) * l1 ^ (1 : ℕ) * lB ^ (2 : ℕ) : ℝ) +
    ((18 : ℝ) * B ^ (1 : ℕ) * q ^ (1 : ℕ) * L2 ^ (2 : ℕ) * l1 ^ (1 : ℕ) * lB ^ (1 : ℕ) : ℝ)

noncomputable def pc11_23 : MvPolynomial (Fin 8) ℤ :=
  (MvPolynomial.C (-6 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (6 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (1 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (3 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (-3 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (2 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (6 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (2 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (3 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (3 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (5 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (6 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (-3 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (3 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (5 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (-15 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (3 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (5 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (6 : Fin 8) ^ (2 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (18 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (3 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (5 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (6 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (-12 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (3 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (5 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (6 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (-3 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (3 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (5 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (2 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (12 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (3 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (5 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (12 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (3 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (6 : Fin 8) ^ (3 : ℕ) : MvPolynomial (Fin 8) ℤ)

noncomputable def rc11_23 (B r q V L2 l0 l1 lB : ℝ) : ℝ :=
  ((-6 : ℝ) * B ^ (1 : ℕ) * q ^ (1 : ℕ) * L2 ^ (2 : ℕ) * l1 ^ (1 : ℕ) : ℝ) +
    ((1 : ℝ) * B ^ (1 : ℕ) * q ^ (1 : ℕ) * L2 ^ (2 : ℕ) * lB ^ (3 : ℕ) : ℝ) +
    ((-3 : ℝ) * B ^ (1 : ℕ) * q ^ (1 : ℕ) * L2 ^ (2 : ℕ) * lB ^ (2 : ℕ) : ℝ) +
    ((6 : ℝ) * B ^ (1 : ℕ) * q ^ (1 : ℕ) * L2 ^ (2 : ℕ) * lB ^ (1 : ℕ) : ℝ) +
    ((3 : ℝ) * B ^ (1 : ℕ) * V ^ (1 : ℕ) * L2 ^ (2 : ℕ) * l0 ^ (2 : ℕ) * l1 ^ (1 : ℕ) : ℝ) +
    ((-3 : ℝ) * B ^ (1 : ℕ) * V ^ (1 : ℕ) * L2 ^ (2 : ℕ) * l0 ^ (2 : ℕ) * lB ^ (1 : ℕ) : ℝ) +
    ((-15 : ℝ) * B ^ (1 : ℕ) * V ^ (1 : ℕ) * L2 ^ (2 : ℕ) * l0 ^ (1 : ℕ) * l1 ^ (2 : ℕ) : ℝ) +
    ((18 : ℝ) * B ^ (1 : ℕ) * V ^ (1 : ℕ) * L2 ^ (2 : ℕ) * l0 ^ (1 : ℕ) * l1 ^ (1 : ℕ) * lB ^ (1 : ℕ) : ℝ) +
    ((-12 : ℝ) * B ^ (1 : ℕ) * V ^ (1 : ℕ) * L2 ^ (2 : ℕ) * l0 ^ (1 : ℕ) * l1 ^ (1 : ℕ) : ℝ) +
    ((-3 : ℝ) * B ^ (1 : ℕ) * V ^ (1 : ℕ) * L2 ^ (2 : ℕ) * l0 ^ (1 : ℕ) * lB ^ (2 : ℕ) : ℝ) +
    ((12 : ℝ) * B ^ (1 : ℕ) * V ^ (1 : ℕ) * L2 ^ (2 : ℕ) * l0 ^ (1 : ℕ) * lB ^ (1 : ℕ) : ℝ) +
    ((12 : ℝ) * B ^ (1 : ℕ) * V ^ (1 : ℕ) * L2 ^ (2 : ℕ) * l1 ^ (3 : ℕ) : ℝ)

noncomputable def pc11_24 : MvPolynomial (Fin 8) ℤ :=
  (MvPolynomial.C (-15 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (3 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (6 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (12 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (3 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (6 : Fin 8) ^ (2 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (3 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (3 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (6 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (2 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (-12 : ℤ) * MvPolynomial.X (0 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (3 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (6 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (7 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (-3 : ℤ) * MvPolynomial.X (2 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (5 : Fin 8) ^ (3 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (9 : ℤ) * MvPolynomial.X (2 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (5 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (6 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (-9 : ℤ) * MvPolynomial.X (2 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (5 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (6 : Fin 8) ^ (2 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (3 : ℤ) * MvPolynomial.X (2 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (6 : Fin 8) ^ (3 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (3 : ℤ) * MvPolynomial.X (2 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (5 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (6 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (6 : ℤ) * MvPolynomial.X (2 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (5 : Fin 8) ^ (2 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (-6 : ℤ) * MvPolynomial.X (2 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (5 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (6 : Fin 8) ^ (2 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (-12 : ℤ) * MvPolynomial.X (2 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (5 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (6 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ)

noncomputable def rc11_24 (B r q V L2 l0 l1 lB : ℝ) : ℝ :=
  ((-15 : ℝ) * B ^ (1 : ℕ) * V ^ (1 : ℕ) * L2 ^ (2 : ℕ) * l1 ^ (2 : ℕ) * lB ^ (1 : ℕ) : ℝ) +
    ((12 : ℝ) * B ^ (1 : ℕ) * V ^ (1 : ℕ) * L2 ^ (2 : ℕ) * l1 ^ (2 : ℕ) : ℝ) +
    ((3 : ℝ) * B ^ (1 : ℕ) * V ^ (1 : ℕ) * L2 ^ (2 : ℕ) * l1 ^ (1 : ℕ) * lB ^ (2 : ℕ) : ℝ) +
    ((-12 : ℝ) * B ^ (1 : ℕ) * V ^ (1 : ℕ) * L2 ^ (2 : ℕ) * l1 ^ (1 : ℕ) * lB ^ (1 : ℕ) : ℝ) +
    ((-3 : ℝ) * q ^ (3 : ℕ) * l0 ^ (3 : ℕ) : ℝ) +
    ((9 : ℝ) * q ^ (3 : ℕ) * l0 ^ (2 : ℕ) * l1 ^ (1 : ℕ) : ℝ) +
    ((-9 : ℝ) * q ^ (3 : ℕ) * l0 ^ (1 : ℕ) * l1 ^ (2 : ℕ) : ℝ) +
    ((3 : ℝ) * q ^ (3 : ℕ) * l1 ^ (3 : ℕ) : ℝ) +
    ((3 : ℝ) * q ^ (2 : ℕ) * L2 ^ (1 : ℕ) * l0 ^ (2 : ℕ) * l1 ^ (1 : ℕ) : ℝ) +
    ((6 : ℝ) * q ^ (2 : ℕ) * L2 ^ (1 : ℕ) * l0 ^ (2 : ℕ) : ℝ) +
    ((-6 : ℝ) * q ^ (2 : ℕ) * L2 ^ (1 : ℕ) * l0 ^ (1 : ℕ) * l1 ^ (2 : ℕ) : ℝ) +
    ((-12 : ℝ) * q ^ (2 : ℕ) * L2 ^ (1 : ℕ) * l0 ^ (1 : ℕ) * l1 ^ (1 : ℕ) : ℝ)

noncomputable def pc11_25 : MvPolynomial (Fin 8) ℤ :=
  (MvPolynomial.C (3 : ℤ) * MvPolynomial.X (2 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (6 : Fin 8) ^ (3 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (6 : ℤ) * MvPolynomial.X (2 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (6 : Fin 8) ^ (2 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (3 : ℤ) * MvPolynomial.X (2 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (3 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (5 : Fin 8) ^ (3 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (-9 : ℤ) * MvPolynomial.X (2 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (3 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (5 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (6 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (9 : ℤ) * MvPolynomial.X (2 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (3 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (5 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (6 : Fin 8) ^ (2 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (-3 : ℤ) * MvPolynomial.X (2 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (3 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (6 : Fin 8) ^ (3 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (-1 : ℤ) * MvPolynomial.X (2 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (5 : Fin 8) ^ (3 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (-6 : ℤ) * MvPolynomial.X (2 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (5 : Fin 8) ^ (2 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (-6 : ℤ) * MvPolynomial.X (2 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (5 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (1 : ℤ) * MvPolynomial.X (2 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (6 : Fin 8) ^ (3 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (6 : ℤ) * MvPolynomial.X (2 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (6 : Fin 8) ^ (2 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (6 : ℤ) * MvPolynomial.X (2 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (6 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ)

noncomputable def rc11_25 (B r q V L2 l0 l1 lB : ℝ) : ℝ :=
  ((3 : ℝ) * q ^ (2 : ℕ) * L2 ^ (1 : ℕ) * l1 ^ (3 : ℕ) : ℝ) +
    ((6 : ℝ) * q ^ (2 : ℕ) * L2 ^ (1 : ℕ) * l1 ^ (2 : ℕ) : ℝ) +
    ((3 : ℝ) * q ^ (1 : ℕ) * V ^ (1 : ℕ) * L2 ^ (1 : ℕ) * l0 ^ (3 : ℕ) : ℝ) +
    ((-9 : ℝ) * q ^ (1 : ℕ) * V ^ (1 : ℕ) * L2 ^ (1 : ℕ) * l0 ^ (2 : ℕ) * l1 ^ (1 : ℕ) : ℝ) +
    ((9 : ℝ) * q ^ (1 : ℕ) * V ^ (1 : ℕ) * L2 ^ (1 : ℕ) * l0 ^ (1 : ℕ) * l1 ^ (2 : ℕ) : ℝ) +
    ((-3 : ℝ) * q ^ (1 : ℕ) * V ^ (1 : ℕ) * L2 ^ (1 : ℕ) * l1 ^ (3 : ℕ) : ℝ) +
    ((-1 : ℝ) * q ^ (1 : ℕ) * L2 ^ (2 : ℕ) * l0 ^ (3 : ℕ) : ℝ) +
    ((-6 : ℝ) * q ^ (1 : ℕ) * L2 ^ (2 : ℕ) * l0 ^ (2 : ℕ) : ℝ) +
    ((-6 : ℝ) * q ^ (1 : ℕ) * L2 ^ (2 : ℕ) * l0 ^ (1 : ℕ) : ℝ) +
    ((1 : ℝ) * q ^ (1 : ℕ) * L2 ^ (2 : ℕ) * l1 ^ (3 : ℕ) : ℝ) +
    ((6 : ℝ) * q ^ (1 : ℕ) * L2 ^ (2 : ℕ) * l1 ^ (2 : ℕ) : ℝ) +
    ((6 : ℝ) * q ^ (1 : ℕ) * L2 ^ (2 : ℕ) * l1 ^ (1 : ℕ) : ℝ)

noncomputable def pc11_26 : MvPolynomial (Fin 8) ℤ :=
  (MvPolynomial.C (-3 : ℤ) * MvPolynomial.X (3 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (5 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (6 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (-6 : ℤ) * MvPolynomial.X (3 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (5 : Fin 8) ^ (2 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (6 : ℤ) * MvPolynomial.X (3 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (5 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (6 : Fin 8) ^ (2 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (12 : ℤ) * MvPolynomial.X (3 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (5 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (6 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (-3 : ℤ) * MvPolynomial.X (3 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (6 : Fin 8) ^ (3 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (-6 : ℤ) * MvPolynomial.X (3 : Fin 8) ^ (1 : ℕ) * MvPolynomial.X (4 : Fin 8) ^ (2 : ℕ) * MvPolynomial.X (6 : Fin 8) ^ (2 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (-1 : ℤ) * MvPolynomial.X (4 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (5 : Fin 8) ^ (3 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (-6 : ℤ) * MvPolynomial.X (4 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (5 : Fin 8) ^ (2 : ℕ) : MvPolynomial (Fin 8) ℤ) +
    (MvPolynomial.C (-6 : ℤ) * MvPolynomial.X (4 : Fin 8) ^ (3 : ℕ) * MvPolynomial.X (5 : Fin 8) ^ (1 : ℕ) : MvPolynomial (Fin 8) ℤ)

noncomputable def rc11_26 (B r q V L2 l0 l1 lB : ℝ) : ℝ :=
  ((-3 : ℝ) * V ^ (1 : ℕ) * L2 ^ (2 : ℕ) * l0 ^ (2 : ℕ) * l1 ^ (1 : ℕ) : ℝ) +
    ((-6 : ℝ) * V ^ (1 : ℕ) * L2 ^ (2 : ℕ) * l0 ^ (2 : ℕ) : ℝ) +
    ((6 : ℝ) * V ^ (1 : ℕ) * L2 ^ (2 : ℕ) * l0 ^ (1 : ℕ) * l1 ^ (2 : ℕ) : ℝ) +
    ((12 : ℝ) * V ^ (1 : ℕ) * L2 ^ (2 : ℕ) * l0 ^ (1 : ℕ) * l1 ^ (1 : ℕ) : ℝ) +
    ((-3 : ℝ) * V ^ (1 : ℕ) * L2 ^ (2 : ℕ) * l1 ^ (3 : ℕ) : ℝ) +
    ((-6 : ℝ) * V ^ (1 : ℕ) * L2 ^ (2 : ℕ) * l1 ^ (2 : ℕ) : ℝ) +
    ((-1 : ℝ) * L2 ^ (3 : ℕ) * l0 ^ (3 : ℕ) : ℝ) +
    ((-6 : ℝ) * L2 ^ (3 : ℕ) * l0 ^ (2 : ℕ) : ℝ) +
    ((-6 : ℝ) * L2 ^ (3 : ℕ) * l0 ^ (1 : ℕ) : ℝ)

/-- block `c11` (coefficient of `1`). -/
noncomputable def pc11 : MvPolynomial (Fin 8) ℤ :=
  pc11_0 + pc11_1 + pc11_2 + pc11_3 + pc11_4 + pc11_5 + pc11_6 + pc11_7 + pc11_8 + pc11_9 + pc11_10 + pc11_11 + pc11_12 + pc11_13 + pc11_14 + pc11_15 + pc11_16 + pc11_17 + pc11_18 + pc11_19 + pc11_20 + pc11_21 + pc11_22 + pc11_23 + pc11_24 + pc11_25 + pc11_26

noncomputable def rc11 (B r q V L2 l0 l1 lB : ℝ) : ℝ :=
  rc11_0 B r q V L2 l0 l1 lB + rc11_1 B r q V L2 l0 l1 lB + rc11_2 B r q V L2 l0 l1 lB + rc11_3 B r q V L2 l0 l1 lB + rc11_4 B r q V L2 l0 l1 lB + rc11_5 B r q V L2 l0 l1 lB + rc11_6 B r q V L2 l0 l1 lB + rc11_7 B r q V L2 l0 l1 lB + rc11_8 B r q V L2 l0 l1 lB + rc11_9 B r q V L2 l0 l1 lB + rc11_10 B r q V L2 l0 l1 lB + rc11_11 B r q V L2 l0 l1 lB + rc11_12 B r q V L2 l0 l1 lB + rc11_13 B r q V L2 l0 l1 lB + rc11_14 B r q V L2 l0 l1 lB + rc11_15 B r q V L2 l0 l1 lB + rc11_16 B r q V L2 l0 l1 lB + rc11_17 B r q V L2 l0 l1 lB + rc11_18 B r q V L2 l0 l1 lB + rc11_19 B r q V L2 l0 l1 lB + rc11_20 B r q V L2 l0 l1 lB + rc11_21 B r q V L2 l0 l1 lB + rc11_22 B r q V L2 l0 l1 lB + rc11_23 B r q V L2 l0 l1 lB + rc11_24 B r q V L2 l0 l1 lB + rc11_25 B r q V L2 l0 l1 lB + rc11_26 B r q V L2 l0 l1 lB

lemma evr0 (H : S05.Hex) (k : Fin 3) : S05.evars H k 0 = S05.Bv H k := rfl
lemma evr1 (H : S05.Hex) (k : Fin 3) : S05.evars H k 1 = S05.rr H k := rfl
lemma evr2 (H : S05.Hex) (k : Fin 3) : S05.evars H k 2 = S05.q0 H k := rfl
lemma evr3 (H : S05.Hex) (k : Fin 3) : S05.evars H k 3 = S05.V02 H k := rfl
lemma evr4 (H : S05.Hex) (k : Fin 3) : S05.evars H k 4 = S05.Lsq H k := rfl
lemma evr5 (H : S05.Hex) (k : Fin 3) : S05.evars H k 5 = S05.l0 H k := rfl
lemma evr6 (H : S05.Hex) (k : Fin 3) : S05.evars H k 6 = S05.l1 H k := rfl
lemma evr7 (H : S05.Hex) (k : Fin 3) : S05.evars H k 7 = S05.lB H k := rfl

lemma ev_add (H : S05.Hex) (k : Fin 3) (p q : MvPolynomial (Fin 8) ℤ) :
    S05.ev H k (p + q) = S05.ev H k p + S05.ev H k q := MvPolynomial.eval₂_add _ _

lemma ev_c4_0 (H : S05.Hex) (k : Fin 3) : S05.ev H k pc4_0 = rc4_0 (S05.Bv H k) (S05.rr H k) (S05.q0 H k)
    (S05.V02 H k) (S05.Lsq H k) (S05.l0 H k) (S05.l1 H k) (S05.lB H k) := by
  rw [show S05.ev H k pc4_0 = MvPolynomial.eval₂Hom (Int.castRingHom ℝ) (S05.evars H k) pc4_0 from rfl]
  simp only [pc4_0, rc4_0, map_add, map_mul, map_pow, MvPolynomial.eval₂Hom_X', MvPolynomial.eval₂Hom_C,
    map_intCast, map_neg, eq_intCast, evr0, evr1, evr2, evr3, evr4, evr5, evr6, evr7]
  push_cast
  ring

lemma ev_c4_1 (H : S05.Hex) (k : Fin 3) : S05.ev H k pc4_1 = rc4_1 (S05.Bv H k) (S05.rr H k) (S05.q0 H k)
    (S05.V02 H k) (S05.Lsq H k) (S05.l0 H k) (S05.l1 H k) (S05.lB H k) := by
  rw [show S05.ev H k pc4_1 = MvPolynomial.eval₂Hom (Int.castRingHom ℝ) (S05.evars H k) pc4_1 from rfl]
  simp only [pc4_1, rc4_1, map_add, map_mul, map_pow, MvPolynomial.eval₂Hom_X', MvPolynomial.eval₂Hom_C,
    map_intCast, map_neg, eq_intCast, evr0, evr1, evr2, evr3, evr4, evr5, evr6, evr7]
  push_cast
  ring

lemma ev_c4_2 (H : S05.Hex) (k : Fin 3) : S05.ev H k pc4_2 = rc4_2 (S05.Bv H k) (S05.rr H k) (S05.q0 H k)
    (S05.V02 H k) (S05.Lsq H k) (S05.l0 H k) (S05.l1 H k) (S05.lB H k) := by
  rw [show S05.ev H k pc4_2 = MvPolynomial.eval₂Hom (Int.castRingHom ℝ) (S05.evars H k) pc4_2 from rfl]
  simp only [pc4_2, rc4_2, map_add, map_mul, map_pow, MvPolynomial.eval₂Hom_X', MvPolynomial.eval₂Hom_C,
    map_intCast, map_neg, eq_intCast, evr0, evr1, evr2, evr3, evr4, evr5, evr6, evr7]
  push_cast
  ring

lemma ev_c4 (H : S05.Hex) (k : Fin 3) : S05.ev H k pc4 = rc4 (S05.Bv H k) (S05.rr H k) (S05.q0 H k)
    (S05.V02 H k) (S05.Lsq H k) (S05.l0 H k) (S05.l1 H k) (S05.lB H k) := by
  simp only [pc4, rc4, ev_add, ev_c4_0, ev_c4_1, ev_c4_2]

lemma ev_c5_0 (H : S05.Hex) (k : Fin 3) : S05.ev H k pc5_0 = rc5_0 (S05.Bv H k) (S05.rr H k) (S05.q0 H k)
    (S05.V02 H k) (S05.Lsq H k) (S05.l0 H k) (S05.l1 H k) (S05.lB H k) := by
  rw [show S05.ev H k pc5_0 = MvPolynomial.eval₂Hom (Int.castRingHom ℝ) (S05.evars H k) pc5_0 from rfl]
  simp only [pc5_0, rc5_0, map_add, map_mul, map_pow, MvPolynomial.eval₂Hom_X', MvPolynomial.eval₂Hom_C,
    map_intCast, map_neg, eq_intCast, evr0, evr1, evr2, evr3, evr4, evr5, evr6, evr7]
  push_cast
  ring

lemma ev_c5_1 (H : S05.Hex) (k : Fin 3) : S05.ev H k pc5_1 = rc5_1 (S05.Bv H k) (S05.rr H k) (S05.q0 H k)
    (S05.V02 H k) (S05.Lsq H k) (S05.l0 H k) (S05.l1 H k) (S05.lB H k) := by
  rw [show S05.ev H k pc5_1 = MvPolynomial.eval₂Hom (Int.castRingHom ℝ) (S05.evars H k) pc5_1 from rfl]
  simp only [pc5_1, rc5_1, map_add, map_mul, map_pow, MvPolynomial.eval₂Hom_X', MvPolynomial.eval₂Hom_C,
    map_intCast, map_neg, eq_intCast, evr0, evr1, evr2, evr3, evr4, evr5, evr6, evr7]
  push_cast
  ring

lemma ev_c5_2 (H : S05.Hex) (k : Fin 3) : S05.ev H k pc5_2 = rc5_2 (S05.Bv H k) (S05.rr H k) (S05.q0 H k)
    (S05.V02 H k) (S05.Lsq H k) (S05.l0 H k) (S05.l1 H k) (S05.lB H k) := by
  rw [show S05.ev H k pc5_2 = MvPolynomial.eval₂Hom (Int.castRingHom ℝ) (S05.evars H k) pc5_2 from rfl]
  simp only [pc5_2, rc5_2, map_add, map_mul, map_pow, MvPolynomial.eval₂Hom_X', MvPolynomial.eval₂Hom_C,
    map_intCast, map_neg, eq_intCast, evr0, evr1, evr2, evr3, evr4, evr5, evr6, evr7]
  push_cast
  ring

lemma ev_c5_3 (H : S05.Hex) (k : Fin 3) : S05.ev H k pc5_3 = rc5_3 (S05.Bv H k) (S05.rr H k) (S05.q0 H k)
    (S05.V02 H k) (S05.Lsq H k) (S05.l0 H k) (S05.l1 H k) (S05.lB H k) := by
  rw [show S05.ev H k pc5_3 = MvPolynomial.eval₂Hom (Int.castRingHom ℝ) (S05.evars H k) pc5_3 from rfl]
  simp only [pc5_3, rc5_3, map_add, map_mul, map_pow, MvPolynomial.eval₂Hom_X', MvPolynomial.eval₂Hom_C,
    map_intCast, map_neg, eq_intCast, evr0, evr1, evr2, evr3, evr4, evr5, evr6, evr7]
  push_cast
  ring

lemma ev_c5 (H : S05.Hex) (k : Fin 3) : S05.ev H k pc5 = rc5 (S05.Bv H k) (S05.rr H k) (S05.q0 H k)
    (S05.V02 H k) (S05.Lsq H k) (S05.l0 H k) (S05.l1 H k) (S05.lB H k) := by
  simp only [pc5, rc5, ev_add, ev_c5_0, ev_c5_1, ev_c5_2, ev_c5_3]

lemma ev_c6_0 (H : S05.Hex) (k : Fin 3) : S05.ev H k pc6_0 = rc6_0 (S05.Bv H k) (S05.rr H k) (S05.q0 H k)
    (S05.V02 H k) (S05.Lsq H k) (S05.l0 H k) (S05.l1 H k) (S05.lB H k) := by
  rw [show S05.ev H k pc6_0 = MvPolynomial.eval₂Hom (Int.castRingHom ℝ) (S05.evars H k) pc6_0 from rfl]
  simp only [pc6_0, rc6_0, map_add, map_mul, map_pow, MvPolynomial.eval₂Hom_X', MvPolynomial.eval₂Hom_C,
    map_intCast, map_neg, eq_intCast, evr0, evr1, evr2, evr3, evr4, evr5, evr6, evr7]
  push_cast
  ring

lemma ev_c6_1 (H : S05.Hex) (k : Fin 3) : S05.ev H k pc6_1 = rc6_1 (S05.Bv H k) (S05.rr H k) (S05.q0 H k)
    (S05.V02 H k) (S05.Lsq H k) (S05.l0 H k) (S05.l1 H k) (S05.lB H k) := by
  rw [show S05.ev H k pc6_1 = MvPolynomial.eval₂Hom (Int.castRingHom ℝ) (S05.evars H k) pc6_1 from rfl]
  simp only [pc6_1, rc6_1, map_add, map_mul, map_pow, MvPolynomial.eval₂Hom_X', MvPolynomial.eval₂Hom_C,
    map_intCast, map_neg, eq_intCast, evr0, evr1, evr2, evr3, evr4, evr5, evr6, evr7]
  push_cast
  ring

lemma ev_c6_2 (H : S05.Hex) (k : Fin 3) : S05.ev H k pc6_2 = rc6_2 (S05.Bv H k) (S05.rr H k) (S05.q0 H k)
    (S05.V02 H k) (S05.Lsq H k) (S05.l0 H k) (S05.l1 H k) (S05.lB H k) := by
  rw [show S05.ev H k pc6_2 = MvPolynomial.eval₂Hom (Int.castRingHom ℝ) (S05.evars H k) pc6_2 from rfl]
  simp only [pc6_2, rc6_2, map_add, map_mul, map_pow, MvPolynomial.eval₂Hom_X', MvPolynomial.eval₂Hom_C,
    map_intCast, map_neg, eq_intCast, evr0, evr1, evr2, evr3, evr4, evr5, evr6, evr7]
  push_cast
  ring

lemma ev_c6_3 (H : S05.Hex) (k : Fin 3) : S05.ev H k pc6_3 = rc6_3 (S05.Bv H k) (S05.rr H k) (S05.q0 H k)
    (S05.V02 H k) (S05.Lsq H k) (S05.l0 H k) (S05.l1 H k) (S05.lB H k) := by
  rw [show S05.ev H k pc6_3 = MvPolynomial.eval₂Hom (Int.castRingHom ℝ) (S05.evars H k) pc6_3 from rfl]
  simp only [pc6_3, rc6_3, map_add, map_mul, map_pow, MvPolynomial.eval₂Hom_X', MvPolynomial.eval₂Hom_C,
    map_intCast, map_neg, eq_intCast, evr0, evr1, evr2, evr3, evr4, evr5, evr6, evr7]
  push_cast
  ring

lemma ev_c6_4 (H : S05.Hex) (k : Fin 3) : S05.ev H k pc6_4 = rc6_4 (S05.Bv H k) (S05.rr H k) (S05.q0 H k)
    (S05.V02 H k) (S05.Lsq H k) (S05.l0 H k) (S05.l1 H k) (S05.lB H k) := by
  rw [show S05.ev H k pc6_4 = MvPolynomial.eval₂Hom (Int.castRingHom ℝ) (S05.evars H k) pc6_4 from rfl]
  simp only [pc6_4, rc6_4, map_add, map_mul, map_pow, MvPolynomial.eval₂Hom_X', MvPolynomial.eval₂Hom_C,
    map_intCast, map_neg, eq_intCast, evr0, evr1, evr2, evr3, evr4, evr5, evr6, evr7]
  push_cast
  ring

lemma ev_c6 (H : S05.Hex) (k : Fin 3) : S05.ev H k pc6 = rc6 (S05.Bv H k) (S05.rr H k) (S05.q0 H k)
    (S05.V02 H k) (S05.Lsq H k) (S05.l0 H k) (S05.l1 H k) (S05.lB H k) := by
  simp only [pc6, rc6, ev_add, ev_c6_0, ev_c6_1, ev_c6_2, ev_c6_3, ev_c6_4]

lemma ev_c7_0 (H : S05.Hex) (k : Fin 3) : S05.ev H k pc7_0 = rc7_0 (S05.Bv H k) (S05.rr H k) (S05.q0 H k)
    (S05.V02 H k) (S05.Lsq H k) (S05.l0 H k) (S05.l1 H k) (S05.lB H k) := by
  rw [show S05.ev H k pc7_0 = MvPolynomial.eval₂Hom (Int.castRingHom ℝ) (S05.evars H k) pc7_0 from rfl]
  simp only [pc7_0, rc7_0, map_add, map_mul, map_pow, MvPolynomial.eval₂Hom_X', MvPolynomial.eval₂Hom_C,
    map_intCast, map_neg, eq_intCast, evr0, evr1, evr2, evr3, evr4, evr5, evr6, evr7]
  push_cast
  ring

lemma ev_c7_1 (H : S05.Hex) (k : Fin 3) : S05.ev H k pc7_1 = rc7_1 (S05.Bv H k) (S05.rr H k) (S05.q0 H k)
    (S05.V02 H k) (S05.Lsq H k) (S05.l0 H k) (S05.l1 H k) (S05.lB H k) := by
  rw [show S05.ev H k pc7_1 = MvPolynomial.eval₂Hom (Int.castRingHom ℝ) (S05.evars H k) pc7_1 from rfl]
  simp only [pc7_1, rc7_1, map_add, map_mul, map_pow, MvPolynomial.eval₂Hom_X', MvPolynomial.eval₂Hom_C,
    map_intCast, map_neg, eq_intCast, evr0, evr1, evr2, evr3, evr4, evr5, evr6, evr7]
  push_cast
  ring

lemma ev_c7_2 (H : S05.Hex) (k : Fin 3) : S05.ev H k pc7_2 = rc7_2 (S05.Bv H k) (S05.rr H k) (S05.q0 H k)
    (S05.V02 H k) (S05.Lsq H k) (S05.l0 H k) (S05.l1 H k) (S05.lB H k) := by
  rw [show S05.ev H k pc7_2 = MvPolynomial.eval₂Hom (Int.castRingHom ℝ) (S05.evars H k) pc7_2 from rfl]
  simp only [pc7_2, rc7_2, map_add, map_mul, map_pow, MvPolynomial.eval₂Hom_X', MvPolynomial.eval₂Hom_C,
    map_intCast, map_neg, eq_intCast, evr0, evr1, evr2, evr3, evr4, evr5, evr6, evr7]
  push_cast
  ring

lemma ev_c7_3 (H : S05.Hex) (k : Fin 3) : S05.ev H k pc7_3 = rc7_3 (S05.Bv H k) (S05.rr H k) (S05.q0 H k)
    (S05.V02 H k) (S05.Lsq H k) (S05.l0 H k) (S05.l1 H k) (S05.lB H k) := by
  rw [show S05.ev H k pc7_3 = MvPolynomial.eval₂Hom (Int.castRingHom ℝ) (S05.evars H k) pc7_3 from rfl]
  simp only [pc7_3, rc7_3, map_add, map_mul, map_pow, MvPolynomial.eval₂Hom_X', MvPolynomial.eval₂Hom_C,
    map_intCast, map_neg, eq_intCast, evr0, evr1, evr2, evr3, evr4, evr5, evr6, evr7]
  push_cast
  ring

lemma ev_c7_4 (H : S05.Hex) (k : Fin 3) : S05.ev H k pc7_4 = rc7_4 (S05.Bv H k) (S05.rr H k) (S05.q0 H k)
    (S05.V02 H k) (S05.Lsq H k) (S05.l0 H k) (S05.l1 H k) (S05.lB H k) := by
  rw [show S05.ev H k pc7_4 = MvPolynomial.eval₂Hom (Int.castRingHom ℝ) (S05.evars H k) pc7_4 from rfl]
  simp only [pc7_4, rc7_4, map_add, map_mul, map_pow, MvPolynomial.eval₂Hom_X', MvPolynomial.eval₂Hom_C,
    map_intCast, map_neg, eq_intCast, evr0, evr1, evr2, evr3, evr4, evr5, evr6, evr7]
  push_cast
  ring

lemma ev_c7_5 (H : S05.Hex) (k : Fin 3) : S05.ev H k pc7_5 = rc7_5 (S05.Bv H k) (S05.rr H k) (S05.q0 H k)
    (S05.V02 H k) (S05.Lsq H k) (S05.l0 H k) (S05.l1 H k) (S05.lB H k) := by
  rw [show S05.ev H k pc7_5 = MvPolynomial.eval₂Hom (Int.castRingHom ℝ) (S05.evars H k) pc7_5 from rfl]
  simp only [pc7_5, rc7_5, map_add, map_mul, map_pow, MvPolynomial.eval₂Hom_X', MvPolynomial.eval₂Hom_C,
    map_intCast, map_neg, eq_intCast, evr0, evr1, evr2, evr3, evr4, evr5, evr6, evr7]
  push_cast
  ring

lemma ev_c7_6 (H : S05.Hex) (k : Fin 3) : S05.ev H k pc7_6 = rc7_6 (S05.Bv H k) (S05.rr H k) (S05.q0 H k)
    (S05.V02 H k) (S05.Lsq H k) (S05.l0 H k) (S05.l1 H k) (S05.lB H k) := by
  rw [show S05.ev H k pc7_6 = MvPolynomial.eval₂Hom (Int.castRingHom ℝ) (S05.evars H k) pc7_6 from rfl]
  simp only [pc7_6, rc7_6, map_add, map_mul, map_pow, MvPolynomial.eval₂Hom_X', MvPolynomial.eval₂Hom_C,
    map_intCast, map_neg, eq_intCast, evr0, evr1, evr2, evr3, evr4, evr5, evr6, evr7]
  push_cast
  ring

lemma ev_c7_7 (H : S05.Hex) (k : Fin 3) : S05.ev H k pc7_7 = rc7_7 (S05.Bv H k) (S05.rr H k) (S05.q0 H k)
    (S05.V02 H k) (S05.Lsq H k) (S05.l0 H k) (S05.l1 H k) (S05.lB H k) := by
  rw [show S05.ev H k pc7_7 = MvPolynomial.eval₂Hom (Int.castRingHom ℝ) (S05.evars H k) pc7_7 from rfl]
  simp only [pc7_7, rc7_7, map_add, map_mul, map_pow, MvPolynomial.eval₂Hom_X', MvPolynomial.eval₂Hom_C,
    map_intCast, map_neg, eq_intCast, evr0, evr1, evr2, evr3, evr4, evr5, evr6, evr7]
  push_cast
  ring

lemma ev_c7_8 (H : S05.Hex) (k : Fin 3) : S05.ev H k pc7_8 = rc7_8 (S05.Bv H k) (S05.rr H k) (S05.q0 H k)
    (S05.V02 H k) (S05.Lsq H k) (S05.l0 H k) (S05.l1 H k) (S05.lB H k) := by
  rw [show S05.ev H k pc7_8 = MvPolynomial.eval₂Hom (Int.castRingHom ℝ) (S05.evars H k) pc7_8 from rfl]
  simp only [pc7_8, rc7_8, map_add, map_mul, map_pow, MvPolynomial.eval₂Hom_X', MvPolynomial.eval₂Hom_C,
    map_intCast, map_neg, eq_intCast, evr0, evr1, evr2, evr3, evr4, evr5, evr6, evr7]
  push_cast
  ring

lemma ev_c7_9 (H : S05.Hex) (k : Fin 3) : S05.ev H k pc7_9 = rc7_9 (S05.Bv H k) (S05.rr H k) (S05.q0 H k)
    (S05.V02 H k) (S05.Lsq H k) (S05.l0 H k) (S05.l1 H k) (S05.lB H k) := by
  rw [show S05.ev H k pc7_9 = MvPolynomial.eval₂Hom (Int.castRingHom ℝ) (S05.evars H k) pc7_9 from rfl]
  simp only [pc7_9, rc7_9, map_add, map_mul, map_pow, MvPolynomial.eval₂Hom_X', MvPolynomial.eval₂Hom_C,
    map_intCast, map_neg, eq_intCast, evr0, evr1, evr2, evr3, evr4, evr5, evr6, evr7]
  push_cast
  ring

lemma ev_c7 (H : S05.Hex) (k : Fin 3) : S05.ev H k pc7 = rc7 (S05.Bv H k) (S05.rr H k) (S05.q0 H k)
    (S05.V02 H k) (S05.Lsq H k) (S05.l0 H k) (S05.l1 H k) (S05.lB H k) := by
  simp only [pc7, rc7, ev_add, ev_c7_0, ev_c7_1, ev_c7_2, ev_c7_3, ev_c7_4, ev_c7_5, ev_c7_6, ev_c7_7, ev_c7_8, ev_c7_9]

lemma ev_c10_0 (H : S05.Hex) (k : Fin 3) : S05.ev H k pc10_0 = rc10_0 (S05.Bv H k) (S05.rr H k) (S05.q0 H k)
    (S05.V02 H k) (S05.Lsq H k) (S05.l0 H k) (S05.l1 H k) (S05.lB H k) := by
  rw [show S05.ev H k pc10_0 = MvPolynomial.eval₂Hom (Int.castRingHom ℝ) (S05.evars H k) pc10_0 from rfl]
  simp only [pc10_0, rc10_0, map_add, map_mul, map_pow, MvPolynomial.eval₂Hom_X', MvPolynomial.eval₂Hom_C,
    map_intCast, map_neg, eq_intCast, evr0, evr1, evr2, evr3, evr4, evr5, evr6, evr7]
  push_cast
  ring

lemma ev_c10_1 (H : S05.Hex) (k : Fin 3) : S05.ev H k pc10_1 = rc10_1 (S05.Bv H k) (S05.rr H k) (S05.q0 H k)
    (S05.V02 H k) (S05.Lsq H k) (S05.l0 H k) (S05.l1 H k) (S05.lB H k) := by
  rw [show S05.ev H k pc10_1 = MvPolynomial.eval₂Hom (Int.castRingHom ℝ) (S05.evars H k) pc10_1 from rfl]
  simp only [pc10_1, rc10_1, map_add, map_mul, map_pow, MvPolynomial.eval₂Hom_X', MvPolynomial.eval₂Hom_C,
    map_intCast, map_neg, eq_intCast, evr0, evr1, evr2, evr3, evr4, evr5, evr6, evr7]
  push_cast
  ring

lemma ev_c10_2 (H : S05.Hex) (k : Fin 3) : S05.ev H k pc10_2 = rc10_2 (S05.Bv H k) (S05.rr H k) (S05.q0 H k)
    (S05.V02 H k) (S05.Lsq H k) (S05.l0 H k) (S05.l1 H k) (S05.lB H k) := by
  rw [show S05.ev H k pc10_2 = MvPolynomial.eval₂Hom (Int.castRingHom ℝ) (S05.evars H k) pc10_2 from rfl]
  simp only [pc10_2, rc10_2, map_add, map_mul, map_pow, MvPolynomial.eval₂Hom_X', MvPolynomial.eval₂Hom_C,
    map_intCast, map_neg, eq_intCast, evr0, evr1, evr2, evr3, evr4, evr5, evr6, evr7]
  push_cast
  ring

lemma ev_c10_3 (H : S05.Hex) (k : Fin 3) : S05.ev H k pc10_3 = rc10_3 (S05.Bv H k) (S05.rr H k) (S05.q0 H k)
    (S05.V02 H k) (S05.Lsq H k) (S05.l0 H k) (S05.l1 H k) (S05.lB H k) := by
  rw [show S05.ev H k pc10_3 = MvPolynomial.eval₂Hom (Int.castRingHom ℝ) (S05.evars H k) pc10_3 from rfl]
  simp only [pc10_3, rc10_3, map_add, map_mul, map_pow, MvPolynomial.eval₂Hom_X', MvPolynomial.eval₂Hom_C,
    map_intCast, map_neg, eq_intCast, evr0, evr1, evr2, evr3, evr4, evr5, evr6, evr7]
  push_cast
  ring

lemma ev_c10_4 (H : S05.Hex) (k : Fin 3) : S05.ev H k pc10_4 = rc10_4 (S05.Bv H k) (S05.rr H k) (S05.q0 H k)
    (S05.V02 H k) (S05.Lsq H k) (S05.l0 H k) (S05.l1 H k) (S05.lB H k) := by
  rw [show S05.ev H k pc10_4 = MvPolynomial.eval₂Hom (Int.castRingHom ℝ) (S05.evars H k) pc10_4 from rfl]
  simp only [pc10_4, rc10_4, map_add, map_mul, map_pow, MvPolynomial.eval₂Hom_X', MvPolynomial.eval₂Hom_C,
    map_intCast, map_neg, eq_intCast, evr0, evr1, evr2, evr3, evr4, evr5, evr6, evr7]
  push_cast
  ring

lemma ev_c10_5 (H : S05.Hex) (k : Fin 3) : S05.ev H k pc10_5 = rc10_5 (S05.Bv H k) (S05.rr H k) (S05.q0 H k)
    (S05.V02 H k) (S05.Lsq H k) (S05.l0 H k) (S05.l1 H k) (S05.lB H k) := by
  rw [show S05.ev H k pc10_5 = MvPolynomial.eval₂Hom (Int.castRingHom ℝ) (S05.evars H k) pc10_5 from rfl]
  simp only [pc10_5, rc10_5, map_add, map_mul, map_pow, MvPolynomial.eval₂Hom_X', MvPolynomial.eval₂Hom_C,
    map_intCast, map_neg, eq_intCast, evr0, evr1, evr2, evr3, evr4, evr5, evr6, evr7]
  push_cast
  ring

lemma ev_c10 (H : S05.Hex) (k : Fin 3) : S05.ev H k pc10 = rc10 (S05.Bv H k) (S05.rr H k) (S05.q0 H k)
    (S05.V02 H k) (S05.Lsq H k) (S05.l0 H k) (S05.l1 H k) (S05.lB H k) := by
  simp only [pc10, rc10, ev_add, ev_c10_0, ev_c10_1, ev_c10_2, ev_c10_3, ev_c10_4, ev_c10_5]

lemma ev_c11_0 (H : S05.Hex) (k : Fin 3) : S05.ev H k pc11_0 = rc11_0 (S05.Bv H k) (S05.rr H k) (S05.q0 H k)
    (S05.V02 H k) (S05.Lsq H k) (S05.l0 H k) (S05.l1 H k) (S05.lB H k) := by
  rw [show S05.ev H k pc11_0 = MvPolynomial.eval₂Hom (Int.castRingHom ℝ) (S05.evars H k) pc11_0 from rfl]
  simp only [pc11_0, rc11_0, map_add, map_mul, map_pow, MvPolynomial.eval₂Hom_X', MvPolynomial.eval₂Hom_C,
    map_intCast, map_neg, eq_intCast, evr0, evr1, evr2, evr3, evr4, evr5, evr6, evr7]
  push_cast
  ring

lemma ev_c11_1 (H : S05.Hex) (k : Fin 3) : S05.ev H k pc11_1 = rc11_1 (S05.Bv H k) (S05.rr H k) (S05.q0 H k)
    (S05.V02 H k) (S05.Lsq H k) (S05.l0 H k) (S05.l1 H k) (S05.lB H k) := by
  rw [show S05.ev H k pc11_1 = MvPolynomial.eval₂Hom (Int.castRingHom ℝ) (S05.evars H k) pc11_1 from rfl]
  simp only [pc11_1, rc11_1, map_add, map_mul, map_pow, MvPolynomial.eval₂Hom_X', MvPolynomial.eval₂Hom_C,
    map_intCast, map_neg, eq_intCast, evr0, evr1, evr2, evr3, evr4, evr5, evr6, evr7]
  push_cast
  ring

lemma ev_c11_2 (H : S05.Hex) (k : Fin 3) : S05.ev H k pc11_2 = rc11_2 (S05.Bv H k) (S05.rr H k) (S05.q0 H k)
    (S05.V02 H k) (S05.Lsq H k) (S05.l0 H k) (S05.l1 H k) (S05.lB H k) := by
  rw [show S05.ev H k pc11_2 = MvPolynomial.eval₂Hom (Int.castRingHom ℝ) (S05.evars H k) pc11_2 from rfl]
  simp only [pc11_2, rc11_2, map_add, map_mul, map_pow, MvPolynomial.eval₂Hom_X', MvPolynomial.eval₂Hom_C,
    map_intCast, map_neg, eq_intCast, evr0, evr1, evr2, evr3, evr4, evr5, evr6, evr7]
  push_cast
  ring

lemma ev_c11_3 (H : S05.Hex) (k : Fin 3) : S05.ev H k pc11_3 = rc11_3 (S05.Bv H k) (S05.rr H k) (S05.q0 H k)
    (S05.V02 H k) (S05.Lsq H k) (S05.l0 H k) (S05.l1 H k) (S05.lB H k) := by
  rw [show S05.ev H k pc11_3 = MvPolynomial.eval₂Hom (Int.castRingHom ℝ) (S05.evars H k) pc11_3 from rfl]
  simp only [pc11_3, rc11_3, map_add, map_mul, map_pow, MvPolynomial.eval₂Hom_X', MvPolynomial.eval₂Hom_C,
    map_intCast, map_neg, eq_intCast, evr0, evr1, evr2, evr3, evr4, evr5, evr6, evr7]
  push_cast
  ring

lemma ev_c11_4 (H : S05.Hex) (k : Fin 3) : S05.ev H k pc11_4 = rc11_4 (S05.Bv H k) (S05.rr H k) (S05.q0 H k)
    (S05.V02 H k) (S05.Lsq H k) (S05.l0 H k) (S05.l1 H k) (S05.lB H k) := by
  rw [show S05.ev H k pc11_4 = MvPolynomial.eval₂Hom (Int.castRingHom ℝ) (S05.evars H k) pc11_4 from rfl]
  simp only [pc11_4, rc11_4, map_add, map_mul, map_pow, MvPolynomial.eval₂Hom_X', MvPolynomial.eval₂Hom_C,
    map_intCast, map_neg, eq_intCast, evr0, evr1, evr2, evr3, evr4, evr5, evr6, evr7]
  push_cast
  ring

lemma ev_c11_5 (H : S05.Hex) (k : Fin 3) : S05.ev H k pc11_5 = rc11_5 (S05.Bv H k) (S05.rr H k) (S05.q0 H k)
    (S05.V02 H k) (S05.Lsq H k) (S05.l0 H k) (S05.l1 H k) (S05.lB H k) := by
  rw [show S05.ev H k pc11_5 = MvPolynomial.eval₂Hom (Int.castRingHom ℝ) (S05.evars H k) pc11_5 from rfl]
  simp only [pc11_5, rc11_5, map_add, map_mul, map_pow, MvPolynomial.eval₂Hom_X', MvPolynomial.eval₂Hom_C,
    map_intCast, map_neg, eq_intCast, evr0, evr1, evr2, evr3, evr4, evr5, evr6, evr7]
  push_cast
  ring

lemma ev_c11_6 (H : S05.Hex) (k : Fin 3) : S05.ev H k pc11_6 = rc11_6 (S05.Bv H k) (S05.rr H k) (S05.q0 H k)
    (S05.V02 H k) (S05.Lsq H k) (S05.l0 H k) (S05.l1 H k) (S05.lB H k) := by
  rw [show S05.ev H k pc11_6 = MvPolynomial.eval₂Hom (Int.castRingHom ℝ) (S05.evars H k) pc11_6 from rfl]
  simp only [pc11_6, rc11_6, map_add, map_mul, map_pow, MvPolynomial.eval₂Hom_X', MvPolynomial.eval₂Hom_C,
    map_intCast, map_neg, eq_intCast, evr0, evr1, evr2, evr3, evr4, evr5, evr6, evr7]
  push_cast
  ring

lemma ev_c11_7 (H : S05.Hex) (k : Fin 3) : S05.ev H k pc11_7 = rc11_7 (S05.Bv H k) (S05.rr H k) (S05.q0 H k)
    (S05.V02 H k) (S05.Lsq H k) (S05.l0 H k) (S05.l1 H k) (S05.lB H k) := by
  rw [show S05.ev H k pc11_7 = MvPolynomial.eval₂Hom (Int.castRingHom ℝ) (S05.evars H k) pc11_7 from rfl]
  simp only [pc11_7, rc11_7, map_add, map_mul, map_pow, MvPolynomial.eval₂Hom_X', MvPolynomial.eval₂Hom_C,
    map_intCast, map_neg, eq_intCast, evr0, evr1, evr2, evr3, evr4, evr5, evr6, evr7]
  push_cast
  ring

lemma ev_c11_8 (H : S05.Hex) (k : Fin 3) : S05.ev H k pc11_8 = rc11_8 (S05.Bv H k) (S05.rr H k) (S05.q0 H k)
    (S05.V02 H k) (S05.Lsq H k) (S05.l0 H k) (S05.l1 H k) (S05.lB H k) := by
  rw [show S05.ev H k pc11_8 = MvPolynomial.eval₂Hom (Int.castRingHom ℝ) (S05.evars H k) pc11_8 from rfl]
  simp only [pc11_8, rc11_8, map_add, map_mul, map_pow, MvPolynomial.eval₂Hom_X', MvPolynomial.eval₂Hom_C,
    map_intCast, map_neg, eq_intCast, evr0, evr1, evr2, evr3, evr4, evr5, evr6, evr7]
  push_cast
  ring

lemma ev_c11_9 (H : S05.Hex) (k : Fin 3) : S05.ev H k pc11_9 = rc11_9 (S05.Bv H k) (S05.rr H k) (S05.q0 H k)
    (S05.V02 H k) (S05.Lsq H k) (S05.l0 H k) (S05.l1 H k) (S05.lB H k) := by
  rw [show S05.ev H k pc11_9 = MvPolynomial.eval₂Hom (Int.castRingHom ℝ) (S05.evars H k) pc11_9 from rfl]
  simp only [pc11_9, rc11_9, map_add, map_mul, map_pow, MvPolynomial.eval₂Hom_X', MvPolynomial.eval₂Hom_C,
    map_intCast, map_neg, eq_intCast, evr0, evr1, evr2, evr3, evr4, evr5, evr6, evr7]
  push_cast
  ring

lemma ev_c11_10 (H : S05.Hex) (k : Fin 3) : S05.ev H k pc11_10 = rc11_10 (S05.Bv H k) (S05.rr H k) (S05.q0 H k)
    (S05.V02 H k) (S05.Lsq H k) (S05.l0 H k) (S05.l1 H k) (S05.lB H k) := by
  rw [show S05.ev H k pc11_10 = MvPolynomial.eval₂Hom (Int.castRingHom ℝ) (S05.evars H k) pc11_10 from rfl]
  simp only [pc11_10, rc11_10, map_add, map_mul, map_pow, MvPolynomial.eval₂Hom_X', MvPolynomial.eval₂Hom_C,
    map_intCast, map_neg, eq_intCast, evr0, evr1, evr2, evr3, evr4, evr5, evr6, evr7]
  push_cast
  ring

lemma ev_c11_11 (H : S05.Hex) (k : Fin 3) : S05.ev H k pc11_11 = rc11_11 (S05.Bv H k) (S05.rr H k) (S05.q0 H k)
    (S05.V02 H k) (S05.Lsq H k) (S05.l0 H k) (S05.l1 H k) (S05.lB H k) := by
  rw [show S05.ev H k pc11_11 = MvPolynomial.eval₂Hom (Int.castRingHom ℝ) (S05.evars H k) pc11_11 from rfl]
  simp only [pc11_11, rc11_11, map_add, map_mul, map_pow, MvPolynomial.eval₂Hom_X', MvPolynomial.eval₂Hom_C,
    map_intCast, map_neg, eq_intCast, evr0, evr1, evr2, evr3, evr4, evr5, evr6, evr7]
  push_cast
  ring

lemma ev_c11_12 (H : S05.Hex) (k : Fin 3) : S05.ev H k pc11_12 = rc11_12 (S05.Bv H k) (S05.rr H k) (S05.q0 H k)
    (S05.V02 H k) (S05.Lsq H k) (S05.l0 H k) (S05.l1 H k) (S05.lB H k) := by
  rw [show S05.ev H k pc11_12 = MvPolynomial.eval₂Hom (Int.castRingHom ℝ) (S05.evars H k) pc11_12 from rfl]
  simp only [pc11_12, rc11_12, map_add, map_mul, map_pow, MvPolynomial.eval₂Hom_X', MvPolynomial.eval₂Hom_C,
    map_intCast, map_neg, eq_intCast, evr0, evr1, evr2, evr3, evr4, evr5, evr6, evr7]
  push_cast
  ring

lemma ev_c11_13 (H : S05.Hex) (k : Fin 3) : S05.ev H k pc11_13 = rc11_13 (S05.Bv H k) (S05.rr H k) (S05.q0 H k)
    (S05.V02 H k) (S05.Lsq H k) (S05.l0 H k) (S05.l1 H k) (S05.lB H k) := by
  rw [show S05.ev H k pc11_13 = MvPolynomial.eval₂Hom (Int.castRingHom ℝ) (S05.evars H k) pc11_13 from rfl]
  simp only [pc11_13, rc11_13, map_add, map_mul, map_pow, MvPolynomial.eval₂Hom_X', MvPolynomial.eval₂Hom_C,
    map_intCast, map_neg, eq_intCast, evr0, evr1, evr2, evr3, evr4, evr5, evr6, evr7]
  push_cast
  ring

lemma ev_c11_14 (H : S05.Hex) (k : Fin 3) : S05.ev H k pc11_14 = rc11_14 (S05.Bv H k) (S05.rr H k) (S05.q0 H k)
    (S05.V02 H k) (S05.Lsq H k) (S05.l0 H k) (S05.l1 H k) (S05.lB H k) := by
  rw [show S05.ev H k pc11_14 = MvPolynomial.eval₂Hom (Int.castRingHom ℝ) (S05.evars H k) pc11_14 from rfl]
  simp only [pc11_14, rc11_14, map_add, map_mul, map_pow, MvPolynomial.eval₂Hom_X', MvPolynomial.eval₂Hom_C,
    map_intCast, map_neg, eq_intCast, evr0, evr1, evr2, evr3, evr4, evr5, evr6, evr7]
  push_cast
  ring

lemma ev_c11_15 (H : S05.Hex) (k : Fin 3) : S05.ev H k pc11_15 = rc11_15 (S05.Bv H k) (S05.rr H k) (S05.q0 H k)
    (S05.V02 H k) (S05.Lsq H k) (S05.l0 H k) (S05.l1 H k) (S05.lB H k) := by
  rw [show S05.ev H k pc11_15 = MvPolynomial.eval₂Hom (Int.castRingHom ℝ) (S05.evars H k) pc11_15 from rfl]
  simp only [pc11_15, rc11_15, map_add, map_mul, map_pow, MvPolynomial.eval₂Hom_X', MvPolynomial.eval₂Hom_C,
    map_intCast, map_neg, eq_intCast, evr0, evr1, evr2, evr3, evr4, evr5, evr6, evr7]
  push_cast
  ring

lemma ev_c11_16 (H : S05.Hex) (k : Fin 3) : S05.ev H k pc11_16 = rc11_16 (S05.Bv H k) (S05.rr H k) (S05.q0 H k)
    (S05.V02 H k) (S05.Lsq H k) (S05.l0 H k) (S05.l1 H k) (S05.lB H k) := by
  rw [show S05.ev H k pc11_16 = MvPolynomial.eval₂Hom (Int.castRingHom ℝ) (S05.evars H k) pc11_16 from rfl]
  simp only [pc11_16, rc11_16, map_add, map_mul, map_pow, MvPolynomial.eval₂Hom_X', MvPolynomial.eval₂Hom_C,
    map_intCast, map_neg, eq_intCast, evr0, evr1, evr2, evr3, evr4, evr5, evr6, evr7]
  push_cast
  ring

lemma ev_c11_17 (H : S05.Hex) (k : Fin 3) : S05.ev H k pc11_17 = rc11_17 (S05.Bv H k) (S05.rr H k) (S05.q0 H k)
    (S05.V02 H k) (S05.Lsq H k) (S05.l0 H k) (S05.l1 H k) (S05.lB H k) := by
  rw [show S05.ev H k pc11_17 = MvPolynomial.eval₂Hom (Int.castRingHom ℝ) (S05.evars H k) pc11_17 from rfl]
  simp only [pc11_17, rc11_17, map_add, map_mul, map_pow, MvPolynomial.eval₂Hom_X', MvPolynomial.eval₂Hom_C,
    map_intCast, map_neg, eq_intCast, evr0, evr1, evr2, evr3, evr4, evr5, evr6, evr7]
  push_cast
  ring

lemma ev_c11_18 (H : S05.Hex) (k : Fin 3) : S05.ev H k pc11_18 = rc11_18 (S05.Bv H k) (S05.rr H k) (S05.q0 H k)
    (S05.V02 H k) (S05.Lsq H k) (S05.l0 H k) (S05.l1 H k) (S05.lB H k) := by
  rw [show S05.ev H k pc11_18 = MvPolynomial.eval₂Hom (Int.castRingHom ℝ) (S05.evars H k) pc11_18 from rfl]
  simp only [pc11_18, rc11_18, map_add, map_mul, map_pow, MvPolynomial.eval₂Hom_X', MvPolynomial.eval₂Hom_C,
    map_intCast, map_neg, eq_intCast, evr0, evr1, evr2, evr3, evr4, evr5, evr6, evr7]
  push_cast
  ring

lemma ev_c11_19 (H : S05.Hex) (k : Fin 3) : S05.ev H k pc11_19 = rc11_19 (S05.Bv H k) (S05.rr H k) (S05.q0 H k)
    (S05.V02 H k) (S05.Lsq H k) (S05.l0 H k) (S05.l1 H k) (S05.lB H k) := by
  rw [show S05.ev H k pc11_19 = MvPolynomial.eval₂Hom (Int.castRingHom ℝ) (S05.evars H k) pc11_19 from rfl]
  simp only [pc11_19, rc11_19, map_add, map_mul, map_pow, MvPolynomial.eval₂Hom_X', MvPolynomial.eval₂Hom_C,
    map_intCast, map_neg, eq_intCast, evr0, evr1, evr2, evr3, evr4, evr5, evr6, evr7]
  push_cast
  ring

lemma ev_c11_20 (H : S05.Hex) (k : Fin 3) : S05.ev H k pc11_20 = rc11_20 (S05.Bv H k) (S05.rr H k) (S05.q0 H k)
    (S05.V02 H k) (S05.Lsq H k) (S05.l0 H k) (S05.l1 H k) (S05.lB H k) := by
  rw [show S05.ev H k pc11_20 = MvPolynomial.eval₂Hom (Int.castRingHom ℝ) (S05.evars H k) pc11_20 from rfl]
  simp only [pc11_20, rc11_20, map_add, map_mul, map_pow, MvPolynomial.eval₂Hom_X', MvPolynomial.eval₂Hom_C,
    map_intCast, map_neg, eq_intCast, evr0, evr1, evr2, evr3, evr4, evr5, evr6, evr7]
  push_cast
  ring

lemma ev_c11_21 (H : S05.Hex) (k : Fin 3) : S05.ev H k pc11_21 = rc11_21 (S05.Bv H k) (S05.rr H k) (S05.q0 H k)
    (S05.V02 H k) (S05.Lsq H k) (S05.l0 H k) (S05.l1 H k) (S05.lB H k) := by
  rw [show S05.ev H k pc11_21 = MvPolynomial.eval₂Hom (Int.castRingHom ℝ) (S05.evars H k) pc11_21 from rfl]
  simp only [pc11_21, rc11_21, map_add, map_mul, map_pow, MvPolynomial.eval₂Hom_X', MvPolynomial.eval₂Hom_C,
    map_intCast, map_neg, eq_intCast, evr0, evr1, evr2, evr3, evr4, evr5, evr6, evr7]
  push_cast
  ring

lemma ev_c11_22 (H : S05.Hex) (k : Fin 3) : S05.ev H k pc11_22 = rc11_22 (S05.Bv H k) (S05.rr H k) (S05.q0 H k)
    (S05.V02 H k) (S05.Lsq H k) (S05.l0 H k) (S05.l1 H k) (S05.lB H k) := by
  rw [show S05.ev H k pc11_22 = MvPolynomial.eval₂Hom (Int.castRingHom ℝ) (S05.evars H k) pc11_22 from rfl]
  simp only [pc11_22, rc11_22, map_add, map_mul, map_pow, MvPolynomial.eval₂Hom_X', MvPolynomial.eval₂Hom_C,
    map_intCast, map_neg, eq_intCast, evr0, evr1, evr2, evr3, evr4, evr5, evr6, evr7]
  push_cast
  ring

lemma ev_c11_23 (H : S05.Hex) (k : Fin 3) : S05.ev H k pc11_23 = rc11_23 (S05.Bv H k) (S05.rr H k) (S05.q0 H k)
    (S05.V02 H k) (S05.Lsq H k) (S05.l0 H k) (S05.l1 H k) (S05.lB H k) := by
  rw [show S05.ev H k pc11_23 = MvPolynomial.eval₂Hom (Int.castRingHom ℝ) (S05.evars H k) pc11_23 from rfl]
  simp only [pc11_23, rc11_23, map_add, map_mul, map_pow, MvPolynomial.eval₂Hom_X', MvPolynomial.eval₂Hom_C,
    map_intCast, map_neg, eq_intCast, evr0, evr1, evr2, evr3, evr4, evr5, evr6, evr7]
  push_cast
  ring

lemma ev_c11_24 (H : S05.Hex) (k : Fin 3) : S05.ev H k pc11_24 = rc11_24 (S05.Bv H k) (S05.rr H k) (S05.q0 H k)
    (S05.V02 H k) (S05.Lsq H k) (S05.l0 H k) (S05.l1 H k) (S05.lB H k) := by
  rw [show S05.ev H k pc11_24 = MvPolynomial.eval₂Hom (Int.castRingHom ℝ) (S05.evars H k) pc11_24 from rfl]
  simp only [pc11_24, rc11_24, map_add, map_mul, map_pow, MvPolynomial.eval₂Hom_X', MvPolynomial.eval₂Hom_C,
    map_intCast, map_neg, eq_intCast, evr0, evr1, evr2, evr3, evr4, evr5, evr6, evr7]
  push_cast
  ring

lemma ev_c11_25 (H : S05.Hex) (k : Fin 3) : S05.ev H k pc11_25 = rc11_25 (S05.Bv H k) (S05.rr H k) (S05.q0 H k)
    (S05.V02 H k) (S05.Lsq H k) (S05.l0 H k) (S05.l1 H k) (S05.lB H k) := by
  rw [show S05.ev H k pc11_25 = MvPolynomial.eval₂Hom (Int.castRingHom ℝ) (S05.evars H k) pc11_25 from rfl]
  simp only [pc11_25, rc11_25, map_add, map_mul, map_pow, MvPolynomial.eval₂Hom_X', MvPolynomial.eval₂Hom_C,
    map_intCast, map_neg, eq_intCast, evr0, evr1, evr2, evr3, evr4, evr5, evr6, evr7]
  push_cast
  ring

lemma ev_c11_26 (H : S05.Hex) (k : Fin 3) : S05.ev H k pc11_26 = rc11_26 (S05.Bv H k) (S05.rr H k) (S05.q0 H k)
    (S05.V02 H k) (S05.Lsq H k) (S05.l0 H k) (S05.l1 H k) (S05.lB H k) := by
  rw [show S05.ev H k pc11_26 = MvPolynomial.eval₂Hom (Int.castRingHom ℝ) (S05.evars H k) pc11_26 from rfl]
  simp only [pc11_26, rc11_26, map_add, map_mul, map_pow, MvPolynomial.eval₂Hom_X', MvPolynomial.eval₂Hom_C,
    map_intCast, map_neg, eq_intCast, evr0, evr1, evr2, evr3, evr4, evr5, evr6, evr7]
  push_cast
  ring

lemma ev_c11 (H : S05.Hex) (k : Fin 3) : S05.ev H k pc11 = rc11 (S05.Bv H k) (S05.rr H k) (S05.q0 H k)
    (S05.V02 H k) (S05.Lsq H k) (S05.l0 H k) (S05.l1 H k) (S05.lB H k) := by
  simp only [pc11, rc11, ev_add, ev_c11_0, ev_c11_1, ev_c11_2, ev_c11_3, ev_c11_4, ev_c11_5, ev_c11_6, ev_c11_7, ev_c11_8, ev_c11_9, ev_c11_10, ev_c11_11, ev_c11_12, ev_c11_13, ev_c11_14, ev_c11_15, ev_c11_16, ev_c11_17, ev_c11_18, ev_c11_19, ev_c11_20, ev_c11_21, ev_c11_22, ev_c11_23, ev_c11_24, ev_c11_25, ev_c11_26]

/-! ## 9. assembly -/

noncomputable def CFa0001 (x T1 T2 T3 T3p P : ℝ) : ℝ :=
  (T1 + ((-1) * x) + ((-1) * T1 * x)) / (x)

lemma CF0001_eq (x : ℝ) : CF0001 x = CFa0001 x (S05.cL1 x) (D00.Li2 x) (Z03.Li3 x) (L3p x) (π ^ 2 / 6) := rfl

noncomputable def CFa0002 (x T1 T2 T3 T3p P : ℝ) : ℝ :=
  (((-1) * (T1) ^ 2) + ((-2) * T1) + (2 * x) + (x * (T1) ^ 2) + (2 * T1 * x)) / (x)

lemma CF0002_eq (x : ℝ) : CF0002 x = CFa0002 x (S05.cL1 x) (D00.Li2 x) (Z03.Li3 x) (L3p x) (π ^ 2 / 6) := rfl

noncomputable def CFa0003 (x T1 T2 T3 T3p P : ℝ) : ℝ :=
  ((T1) ^ 3 + ((-6) * x) + (3 * (T1) ^ 2) + (6 * T1) + ((-1) * x * (T1) ^ 3) + ((-6) * T1 * x) + ((-3) * x * (T1) ^ 2)) / (x)

lemma CF0003_eq (x : ℝ) : CF0003 x = CFa0003 x (S05.cL1 x) (D00.Li2 x) (Z03.Li3 x) (L3p x) (π ^ 2 / 6) := rfl

noncomputable def CFa0011 (x T1 T2 T3 T3p P : ℝ) : ℝ :=
  (((-1) * (T1) ^ 2) + ((-2) * T1) + ((-2) * T2) + (4 * x) + (x * (T1) ^ 2) + (2 * T1 * x) + (2 * T2 * x)) / ((2 * x))

lemma CF0011_eq (x : ℝ) : CF0011 x = CFa0011 x (S05.cL1 x) (D00.Li2 x) (Z03.Li3 x) (L3p x) (π ^ 2 / 6) := rfl

noncomputable def CFa0012 (x T1 T2 T3 T3p P : ℝ) : ℝ :=
  ((T1) ^ 3 + ((-18) * x) + ((-6) * T3) + ((-6) * T3p) + (6 * T2) + (6 * (T1) ^ 2) + (12 * T1) + ((-1) * x * (T1) ^ 3) + ((-12) * T1 * x) + ((-6) * T2 * x) + ((-6) * x * (T1) ^ 2) + (6 * T3 * x) + (6 * T3p * x)) / ((3 * x))

lemma CF0012_eq (x : ℝ) : CF0012 x = CFa0012 x (S05.cL1 x) (D00.Li2 x) (Z03.Li3 x) (L3p x) (π ^ 2 / 6) := rfl

noncomputable def CFa0021 (x T1 T2 T3 T3p P : ℝ) : ℝ :=
  ((T1) ^ 2 + ((-6) * x) + ((-2) * T3p) + (2 * T1) + (2 * T2) + ((-1) * x * (T1) ^ 2) + ((-2) * T1 * x) + ((-2) * T2 * x) + (2 * T3p * x)) / (x)

lemma CF0021_eq (x : ℝ) : CF0021 x = CFa0021 x (S05.cL1 x) (D00.Li2 x) (Z03.Li3 x) (L3p x) (π ^ 2 / 6) := rfl

noncomputable def CFa0101 (x T1 T2 T3 T3p P : ℝ) : ℝ :=
  (((-1) * T1) + ((-1) * T2) + (2 * x) + (T1 * x)) / (x)

lemma CF0101_eq (x : ℝ) : CF0101 x = CFa0101 x (S05.cL1 x) (D00.Li2 x) (Z03.Li3 x) (L3p x) (π ^ 2 / 6) := rfl

noncomputable def CFa0102 (x T1 T2 T3 T3p P : ℝ) : ℝ :=
  ((T1) ^ 3 + ((-18) * x) + (3 * (T1) ^ 2) + (6 * T2) + (6 * T3) + (6 * T3p) + (12 * T1) + ((-12) * T1 * x) + ((-3) * x * (T1) ^ 2) + (6 * T1 * T2)) / ((3 * x))

lemma CF0102_eq (x : ℝ) : CF0102 x = CFa0102 x (S05.cL1 x) (D00.Li2 x) (Z03.Li3 x) (L3p x) (π ^ 2 / 6) := rfl

noncomputable def CFa0111 (x T1 T2 T3 T3p P : ℝ) : ℝ :=
  (((-1) * (T1) ^ 3) + ((-36) * x) + ((-6) * T3p) + (3 * (T1) ^ 2) + (6 * T3) + (12 * T1) + (12 * T2) + ((-12) * T1 * x) + ((-6) * T1 * P) + ((-6) * T2 * x) + ((-3) * x * (T1) ^ 2) + (2 * x * (T1) ^ 3) + (6 * T3 * x) + (6 * x * P) + (12 * T3p * x) + (6 * T1 * T2 * x) + (6 * T1 * x * P)) / ((6 * x))

lemma CF0111_eq (x : ℝ) : CF0111 x = CFa0111 x (S05.cL1 x) (D00.Li2 x) (Z03.Li3 x) (L3p x) (π ^ 2 / 6) := rfl

noncomputable def CFa0201 (x T1 T2 T3 T3p P : ℝ) : ℝ :=
  (((-6) * x) + (2 * T1) + (2 * T2) + (2 * T3) + ((-2) * T1 * x)) / (x)

lemma CF0201_eq (x : ℝ) : CF0201 x = CFa0201 x (S05.cL1 x) (D00.Li2 x) (Z03.Li3 x) (L3p x) (π ^ 2 / 6) := rfl

noncomputable def CFa1001 (x T1 T2 T3 T3p P : ℝ) : ℝ :=
  (((-1) * (x) ^ 2) + ((-2) * x) + (2 * T1) + ((-2) * T1 * (x) ^ 2)) / ((4 * (x) ^ 2))

lemma CF1001_eq (x : ℝ) : CF1001 x = CFa1001 x (S05.cL1 x) (D00.Li2 x) (Z03.Li3 x) (L3p x) (π ^ 2 / 6) := rfl

noncomputable def CFa1002 (x T1 T2 T3 T3p P : ℝ) : ℝ :=
  ((x) ^ 2 + ((-6) * T1) + ((-2) * (T1) ^ 2) + (6 * x) + (2 * T1 * (x) ^ 2) + (2 * (T1) ^ 2 * (x) ^ 2) + (4 * T1 * x)) / ((4 * (x) ^ 2))

lemma CF1002_eq (x : ℝ) : CF1002 x = CFa1002 x (S05.cL1 x) (D00.Li2 x) (Z03.Li3 x) (L3p x) (π ^ 2 / 6) := rfl

noncomputable def CFa1003 (x T1 T2 T3 T3p P : ℝ) : ℝ :=
  (((-42) * x) + ((-3) * (x) ^ 2) + (4 * (T1) ^ 3) + (18 * (T1) ^ 2) + (42 * T1) + ((-36) * T1 * x) + ((-12) * x * (T1) ^ 2) + ((-6) * T1 * (x) ^ 2) + ((-6) * (T1) ^ 2 * (x) ^ 2) + ((-4) * (T1) ^ 3 * (x) ^ 2)) / ((8 * (x) ^ 2))

lemma CF1003_eq (x : ℝ) : CF1003 x = CFa1003 x (S05.cL1 x) (D00.Li2 x) (Z03.Li3 x) (L3p x) (π ^ 2 / 6) := rfl

noncomputable def CFa1011 (x T1 T2 T3 T3p P : ℝ) : ℝ :=
  (((-1) * T1) + ((-1) * (T1) ^ 2) + ((-2) * T2) + (3 * x) + (4 * (x) ^ 2) + ((T1) ^ 2 * (x) ^ 2) + ((-2) * T1 * x) + (2 * T2 * (x) ^ 2) + (3 * T1 * (x) ^ 2)) / ((4 * (x) ^ 2))

lemma CF1011_eq (x : ℝ) : CF1011 x = CFa1011 x (S05.cL1 x) (D00.Li2 x) (Z03.Li3 x) (L3p x) (π ^ 2 / 6) := rfl

noncomputable def CFa1012 (x T1 T2 T3 T3p P : ℝ) : ℝ :=
  (((-84) * x) + ((-51) * (x) ^ 2) + ((-24) * T3) + ((-24) * T3p) + (4 * (T1) ^ 3) + (24 * (T1) ^ 2) + (36 * T2) + (48 * T1) + ((-48) * T1 * (x) ^ 2) + ((-24) * T2 * x) + ((-24) * (T1) ^ 2 * (x) ^ 2) + ((-12) * T2 * (x) ^ 2) + ((-4) * (T1) ^ 3 * (x) ^ 2) + (24 * T3 * (x) ^ 2) + (24 * T3p * (x) ^ 2)) / ((24 * (x) ^ 2))

lemma CF1012_eq (x : ℝ) : CF1012 x = CFa1012 x (S05.cL1 x) (D00.Li2 x) (Z03.Li3 x) (L3p x) (π ^ 2 / 6) := rfl

noncomputable def CFa1021 (x T1 T2 T3 T3p P : ℝ) : ℝ :=
  (((-31) * (x) ^ 2) + ((-14) * x) + ((-8) * T3p) + (2 * T1) + (2 * (T1) ^ 2) + (4 * T2) + ((-14) * T1 * (x) ^ 2) + ((-12) * T2 * (x) ^ 2) + ((-6) * (T1) ^ 2 * (x) ^ 2) + (4 * x * (T1) ^ 2) + (8 * T2 * x) + (8 * T3p * (x) ^ 2) + (12 * T1 * x)) / ((8 * (x) ^ 2))

lemma CF1021_eq (x : ℝ) : CF1021 x = CFa1021 x (S05.cL1 x) (D00.Li2 x) (Z03.Li3 x) (L3p x) (π ^ 2 / 6) := rfl

noncomputable def CFa1101 (x T1 T2 T3 T3p P : ℝ) : ℝ :=
  ((x) ^ 2 + ((-1) * T1) + ((-2) * T2) + (3 * x) + (T1 * (x) ^ 2)) / ((4 * (x) ^ 2))

lemma CF1101_eq (x : ℝ) : CF1101 x = CFa1101 x (S05.cL1 x) (D00.Li2 x) (Z03.Li3 x) (L3p x) (π ^ 2 / 6) := rfl

noncomputable def CFa1102 (x T1 T2 T3 T3p P : ℝ) : ℝ :=
  (((-84) * x) + ((-9) * (x) ^ 2) + (4 * (T1) ^ 3) + (6 * (T1) ^ 2) + (24 * T3) + (24 * T3p) + (36 * T2) + (48 * T1) + ((-36) * T1 * x) + ((-12) * T1 * (x) ^ 2) + ((-6) * (T1) ^ 2 * (x) ^ 2) + (24 * T1 * T2)) / ((24 * (x) ^ 2))

lemma CF1102_eq (x : ℝ) : CF1102 x = CFa1102 x (S05.cL1 x) (D00.Li2 x) (Z03.Li3 x) (L3p x) (π ^ 2 / 6) := rfl

noncomputable def CFa1111 (x T1 T2 T3 T3p P : ℝ) : ℝ :=
  (((-51) * (x) ^ 2) + ((-42) * x) + ((-12) * T3p) + ((-2) * (T1) ^ 3) + (3 * (T1) ^ 2) + (6 * T1) + (12 * T2) + (12 * T3) + ((-24) * T1 * (x) ^ 2) + ((-12) * T1 * P) + ((-6) * T2 * (x) ^ 2) + ((-3) * (T1) ^ 2 * (x) ^ 2) + (4 * (T1) ^ 3 * (x) ^ 2) + (6 * P * (x) ^ 2) + (12 * T2 * x) + (12 * T3 * (x) ^ 2) + (12 * x * P) + (18 * T1 * x) + (24 * T3p * (x) ^ 2) + (12 * T1 * T2 * (x) ^ 2) + (12 * T1 * P * (x) ^ 2)) / ((24 * (x) ^ 2))

lemma CF1111_eq (x : ℝ) : CF1111 x = CFa1111 x (S05.cL1 x) (D00.Li2 x) (Z03.Li3 x) (L3p x) (π ^ 2 / 6) := rfl

noncomputable def CFa1201 (x T1 T2 T3 T3p P : ℝ) : ℝ :=
  (((-14) * x) + ((-3) * (x) ^ 2) + (2 * T1) + (4 * T2) + (8 * T3) + ((-2) * T1 * (x) ^ 2)) / ((8 * (x) ^ 2))

lemma CF1201_eq (x : ℝ) : CF1201 x = CFa1201 x (S05.cL1 x) (D00.Li2 x) (Z03.Li3 x) (L3p x) (π ^ 2 / 6) := rfl

lemma cL1_one_sub (r : ℝ) : S05.cL1 (1 - r) = -Real.log r := by
  unfold S05.cL1; ring_nf

lemma len_alg (L L2 V q X Y : ℝ) (hL : 0 < L) (hLn : L ^ 2 = L2) :
    L ^ 5 * (L * X + (V - q ^ 2 / L2) / L * Y) = L2 ^ 3 * X + L2 ^ 2 * (V - q ^ 2 / L2) * Y := by
  rw [← hLn]; field_simp

set_option maxRecDepth 40000 in
set_option maxHeartbeats 100000000 in
lemma edge_ne (H : S05.Hex) (hH : S05.IsHex H) (k : Fin 3) (hr1 : S05.rr H k ≠ 1) :
    S05.Len H k ^ 5 * S05.Ee H k = S05.closedForm H k pc4 pc5 pc6 pc7 pc10 pc11 := by
  obtain ⟨-, -, -, -, -, -, hr⟩ := edge_facts H hH k
  have hx : 1 - S05.rr H k ≠ 0 := fun h => hr1 (by linarith)
  have hx1 : 1 - S05.rr H k < 1 := by linarith
  have hL := S05.Lsq_pos H hH k
  have hLn : S05.Len H k ^ 2 = S05.Lsq H k := Real.sq_sqrt hL.le
  have hLp : 0 < S05.Len H k := Real.sqrt_pos.2 hL
  rw [Ee_gen H hH k, len_alg _ _ _ _ _ _ hLp hLn]
  simp only [Fin.sum_univ_succ, Fin.sum_univ_zero, kvK, nvK, avK, bvK, cvK, Matrix.cons_val_zero,
    Matrix.cons_val_succ, add_zero]
  simp only [Jx0001 hx1 hx, Jx0002 hx1 hx, Jx0003 hx1 hx, Jx0011 hx1 hx, Jx0012 hx1 hx, Jx0021 hx1 hx, Jx0101 hx1 hx, Jx0102 hx1 hx, Jx0111 hx1 hx, Jx0201 hx1 hx, Jx1001 hx1 hx, Jx1002 hx1 hx, Jx1003 hx1 hx, Jx1011 hx1 hx, Jx1012 hx1 hx, Jx1021 hx1 hx, Jx1101 hx1 hx, Jx1102 hx1 hx, Jx1111 hx1 hx, Jx1201 hx1 hx, Jint_zero, Iv000, Iv001, Iv002, Iv003, Iv010, Iv011, Iv012, Iv020, Iv021, Iv030, Iv100, Iv101, Iv102, Iv103, Iv110, Iv111, Iv112, Iv120, Iv121, Iv130]
  simp only [CF0001_eq, CF0002_eq, CF0003_eq, CF0011_eq, CF0012_eq, CF0021_eq, CF0101_eq, CF0102_eq, CF0111_eq, CF0201_eq, CF1001_eq, CF1002_eq, CF1003_eq, CF1011_eq, CF1012_eq, CF1021_eq, CF1101_eq, CF1102_eq, CF1111_eq, CF1201_eq, cL1_one_sub, eF, gF_one, gF_zero]
  simp only [S05.closedForm, ev_c4, ev_c5, ev_c6, ev_c7, ev_c10, ev_c11, Li3_y]
  simp only [CFa0001, CFa0002, CFa0003, CFa0011, CFa0012, CFa0021, CFa0101, CFa0102, CFa0111, CFa0201, CFa1001, CFa1002, CFa1003, CFa1011, CFa1012, CFa1021, CFa1101, CFa1102, CFa1111, CFa1201, kc0, kc1, kc2, kc3, kc4, kc5, kc6, kc7, kc8, kc9, kc10, kc11, kc12, kc13, kc14, kc15, kc16, kc17, kc18, kc19, kc20, kc21, kc22, kc23, kc24, kc25, kc26, kc27, kc28, kc29, kc30, kc31, kc32, kc33, kc34, kc35, kc36, kc37, kc38, kc39, rc4, rc5, rc6, rc7, rc10, rc11, rc4_0, rc4_1, rc4_2, rc5_0, rc5_1, rc5_2, rc5_3, rc6_0, rc6_1, rc6_2, rc6_3, rc6_4, rc7_0, rc7_1, rc7_2, rc7_3, rc7_4, rc7_5, rc7_6, rc7_7, rc7_8, rc7_9, rc10_0, rc10_1, rc10_2, rc10_3, rc10_4, rc10_5, rc11_0, rc11_1, rc11_2, rc11_3, rc11_4, rc11_5, rc11_6, rc11_7, rc11_8, rc11_9, rc11_10, rc11_11, rc11_12, rc11_13, rc11_14, rc11_15, rc11_16, rc11_17, rc11_18, rc11_19, rc11_20, rc11_21, rc11_22, rc11_23, rc11_24, rc11_25, rc11_26, S05.c1, S05.c2, S05.c3, S05.c8, S05.c9, S05.Dl, S05.gg3_eq]
  have hL0 : S05.Lsq H k ≠ 0 := hL.ne'
  field_simp
  ring

lemma Li3_zero : Z03.Li3 0 = 0 := by simp [Z03.Li3]

lemma Li2_zero : D00.Li2 0 = 0 := D00.proof_2

set_option maxRecDepth 40000 in
set_option maxHeartbeats 100000000 in
lemma edge_one (H : S05.Hex) (hH : S05.IsHex H) (k : Fin 3) (hr1 : S05.rr H k = 1) :
    S05.Len H k ^ 5 * S05.Ee H k = S05.closedForm H k pc4 pc5 pc6 pc7 pc10 pc11 := by
  have hL := S05.Lsq_pos H hH k
  have hLn : S05.Len H k ^ 2 = S05.Lsq H k := Real.sq_sqrt hL.le
  have hLp : 0 < S05.Len H k := Real.sqrt_pos.2 hL
  rw [Ee_gen H hH k, len_alg _ _ _ _ _ _ hLp hLn]
  simp only [Fin.sum_univ_succ, Fin.sum_univ_zero, kvK, nvK, avK, bvK, cvK, Matrix.cons_val_zero,
    Matrix.cons_val_succ, add_zero]
  simp only [S05.closedForm, ev_c4, ev_c5, ev_c6, ev_c7, ev_c10, ev_c11, eF, gF_one, gF_zero]
  simp only [hr1, sub_self, zero_div, Real.log_one, neg_zero, Li2_zero, Li3_zero, Jint_zero]
  simp only [Jint_at_zero 0 0 0 1 (by norm_num), Jint_at_zero 0 0 0 2 (by norm_num), Jint_at_zero 0 0 0 3 (by norm_num), Jint_at_zero 0 0 1 1 (by norm_num), Jint_at_zero 0 0 1 2 (by norm_num), Jint_at_zero 0 0 2 1 (by norm_num), Jint_at_zero 0 1 0 1 (by norm_num), Jint_at_zero 0 1 0 2 (by norm_num), Jint_at_zero 0 1 1 1 (by norm_num), Jint_at_zero 0 2 0 1 (by norm_num), Jint_at_zero 1 0 0 1 (by norm_num), Jint_at_zero 1 0 0 2 (by norm_num), Jint_at_zero 1 0 0 3 (by norm_num), Jint_at_zero 1 0 1 1 (by norm_num), Jint_at_zero 1 0 1 2 (by norm_num), Jint_at_zero 1 0 2 1 (by norm_num), Jint_at_zero 1 1 0 1 (by norm_num), Jint_at_zero 1 1 0 2 (by norm_num), Jint_at_zero 1 1 1 1 (by norm_num), Jint_at_zero 1 2 0 1 (by norm_num), Iv000, Iv001, Iv002, Iv003, Iv010, Iv011, Iv012, Iv020, Iv021, Iv030, Iv100, Iv101, Iv102, Iv103, Iv110, Iv111, Iv112, Iv120, Iv121, Iv130]
  simp only [kc0, kc1, kc2, kc3, kc4, kc5, kc6, kc7, kc8, kc9, kc10, kc11, kc12, kc13, kc14, kc15, kc16, kc17, kc18, kc19, kc20, kc21, kc22, kc23, kc24, kc25, kc26, kc27, kc28, kc29, kc30, kc31, kc32, kc33, kc34, kc35, kc36, kc37, kc38, kc39, rc4, rc5, rc6, rc7, rc10, rc11, rc4_0, rc4_1, rc4_2, rc5_0, rc5_1, rc5_2, rc5_3, rc6_0, rc6_1, rc6_2, rc6_3, rc6_4, rc7_0, rc7_1, rc7_2, rc7_3, rc7_4, rc7_5, rc7_6, rc7_7, rc7_8, rc7_9, rc10_0, rc10_1, rc10_2, rc10_3, rc10_4, rc10_5, rc11_0, rc11_1, rc11_2, rc11_3, rc11_4, rc11_5, rc11_6, rc11_7, rc11_8, rc11_9, rc11_10, rc11_11, rc11_12, rc11_13, rc11_14, rc11_15, rc11_16, rc11_17, rc11_18, rc11_19, rc11_20, rc11_21, rc11_22, rc11_23, rc11_24, rc11_25, rc11_26, S05.c1, S05.c2, S05.c3, S05.c8, S05.c9, S05.Dl, S05.gg3_eq]
  have hL0 : S05.Lsq H k ≠ 0 := hL.ne'
  field_simp
  ring

/-- `S07.stmt_10`, with the six blocks `c_4, c_5, c_6, c_7, c_10, c_11` of the table. -/
theorem main : S07.stmt_10 := by
  refine ⟨pc4, pc5, pc6, pc7, pc10, pc11, fun H hH k => ?_⟩
  by_cases h : S05.rr H k = 1
  · exact edge_one H hH k h
  · exact edge_ne H hH k h

end S07.PA

theorem S07.proof_10 : S07.stmt_10 := S07.PA.main

theorem S07.proof_14 : S07.stmt_14 :=
  ⟨S05.proof_9, S07.proof_10, S05.proof_11, S05.proof_12, S05.proof_13⟩
