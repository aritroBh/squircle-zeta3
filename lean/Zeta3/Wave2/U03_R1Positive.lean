import Mathlib
/-!
# U03 — `R_1 > 0` for `m ≥ 3`, so `P'(1) = 0` only for `m = 2` (paper label `thm:poly`)

Target id: U03.
Source: `paper/ieee/main.tex`, Subsection "Polygons at the Diamond End", lines 1645–1653
(detailed proof: `paper/ieee/proofs/sec_g4_poly.tex`, label `g4:thm:poly`, line 156).

Setup quoted verbatim (main.tex line 1645):
"With $n_i=(\cos\frac{i\pi}m,\sin\frac{i\pi}m)$, let $K_p=\{\sum_{i=1}^m|\langle n_i,x\rangle|^p\le1\}$.
At $p=1$ this is a regular $2m$-gon of perimeter $P(1)=4m(\sec\frac{\pi}{2m}-\cos\frac{\pi}{2m})$,
and $m=2$ is the squircle."

Theorem quoted verbatim (main.tex lines 1647–1653):
"\begin{theorem}[First derivative at the regular polygon$^\dagger$]\label{thm:poly}
$P'(1)=P(1)\,R_1$, where $V=\cot\frac{\pi}{2m}$ and
\[
R_1=\frac1V\sum_{i=1}^{m-1}\sin\frac{i\pi}{m}\,\log\frac{V}{\sin(i\pi/m)} .
\]
In particular $P'(1)=0$ only for $m=2$, and for the hexagon ($m=3$) $P'(1)=2\sqrt3\log2$.
\end{theorem}"

Proof step quoted verbatim (sec_g4_poly.tex line 156):
"For $m=2$: $V=1$, $\sin\frac\pi2=1$, so $R_1=0$. For $m\ge3$: $V=\cot\frac\pi{2m}\ge\cot\frac\pi6=\sqrt3>1\ge\sin\frac{i\pi}m$,
so every term of $R_1$ is $\ge0$ and the term $i=1$ is $>0$, hence $R_1>0$ and $P'(1)>0$."

What is formalized: `V`, `R1`, `P1` are defined exactly by the paper's formulas (same
definitions as Wave1 `T06`; `\sec` written as `1 / cos`; the sum `∑_{i=1}^{m-1}` is
`Finset.Ico 1 m`). The range is `m ≥ 2` (a `2m`-gon needs `m ≥ 2`; this is the paper's
setting, and for `m = 1` the formula degenerates: `V = cot(π/2) = 0` and the sum is empty).
For `m ≥ 2` all divisions are by strictly positive reals and every `log` argument is
positive, so no Lean convention (`x/0 = 0`, `log` of a nonpositive) is involved.
"P'(1)=0 only for m=2" is formalized through the identity `P'(1) = P(1) R_1` as the
statement `P(1) R_1 = 0 ↔ m = 2` for `m ≥ 2`, together with its ingredients
`R_1 > 0` (`m ≥ 3`), `R_1 = 0` (`m = 2`) and `P(1) > 0`.

NOT COVERED:
* The identity `P'(1) = P(1) R_1` itself and the fact that `4m(sec − cos)` is the
  perimeter of `K_1` (Mathlib has no convenient perimeter of a planar convex body as a
  function of `p`). Here `P'(1)` appears only as the product `P(1) R_1`.
* The hexagon value `P'(1) = 2√3 log 2` (already in Wave1 `T06`).
-/
-- BEGIN STATEMENTS
namespace U03

open Real

/-- `V = cot(π/(2m))`. -/
noncomputable def V (m : ℕ) : ℝ := Real.cot (π / (2 * (m : ℝ)))

/-- `R_1 = (1/V) ∑_{i=1}^{m-1} sin(iπ/m) log(V / sin(iπ/m))`. -/
noncomputable def R1 (m : ℕ) : ℝ :=
  (1 / V m) * ∑ i ∈ Finset.Ico 1 m,
    Real.sin ((i : ℝ) * π / (m : ℝ)) * Real.log (V m / Real.sin ((i : ℝ) * π / (m : ℝ)))

/-- The paper's perimeter formula `P(1) = 4m(sec(π/2m) − cos(π/2m))` of the regular `2m`-gon. -/
noncomputable def P1 (m : ℕ) : ℝ :=
  4 * (m : ℝ) * (1 / Real.cos (π / (2 * (m : ℝ))) - Real.cos (π / (2 * (m : ℝ))))

/-- Proof of `thm:poly` (sec_g4_poly.tex l.156): "For $m\ge3$ ... hence $R_1>0$". -/
def stmt_1 : Prop := ∀ m : ℕ, 3 ≤ m → 0 < R1 m

/-- Proof of `thm:poly` (sec_g4_poly.tex l.156): "For $m=2$: ... so $R_1=0$". -/
def stmt_2 : Prop := R1 2 = 0

/-- Combined: for every `m ≥ 2`, `R_1 = 0` exactly when `m = 2`. -/
def stmt_3 : Prop := ∀ m : ℕ, 2 ≤ m → (R1 m = 0 ↔ m = 2)

/-- The perimeter `P(1) = 4m(sec(π/2m) − cos(π/2m))` is positive for every `m ≥ 2`. -/
def stmt_4 : Prop := ∀ m : ℕ, 2 ≤ m → 0 < P1 m

/-- `thm:poly`: "In particular $P'(1)=0$ only for $m=2$", with `P'(1) = P(1) R_1`:
for `m ≥ 2`, `P(1) R_1 = 0 ↔ m = 2`. -/
def stmt_5 : Prop := ∀ m : ℕ, 2 ≤ m → (P1 m * R1 m = 0 ↔ m = 2)

/-- Proof of `thm:poly` (sec_g4_poly.tex l.156): "hence $R_1>0$ and $P'(1)>0$" for `m ≥ 3`,
with `P'(1) = P(1) R_1`. -/
def stmt_6 : Prop := ∀ m : ℕ, 3 ≤ m → 0 < P1 m * R1 m

end U03
-- END STATEMENTS

lemma U03.V_gt_one (m : ℕ) (hm : 3 ≤ m) : 1 < U03.V m := by
  have hm' : (3:ℝ) ≤ m := by exact_mod_cast hm
  have hpos : (0:ℝ) < m := by linarith
  set x := Real.pi / (2 * (m:ℝ)) with hx
  have hx0 : 0 < x := by positivity
  have hx1 : x < Real.pi / 4 := by
    rw [hx, div_lt_div_iff₀ (by positivity) (by norm_num)]
    nlinarith [Real.pi_pos]
  have hsin : 0 < Real.sin x := Real.sin_pos_of_pos_of_lt_pi hx0 (by linarith [Real.pi_pos])
  have h1 : Real.sin x < Real.sin (Real.pi / 4) :=
    Real.strictMonoOn_sin ⟨by linarith [Real.pi_pos], by linarith [Real.pi_pos]⟩
      ⟨by linarith [Real.pi_pos], by linarith [Real.pi_pos]⟩ hx1
  have h2 : Real.cos (Real.pi / 4) < Real.cos x :=
    Real.cos_lt_cos_of_nonneg_of_le_pi_div_two hx0.le (by linarith [Real.pi_pos]) hx1
  rw [Real.cos_pi_div_four] at h2
  rw [Real.sin_pi_div_four] at h1
  unfold U03.V
  rw [Real.cot_eq_cos_div_sin, one_lt_div hsin]
  linarith

lemma U03.R1_pos (m : ℕ) (hm : 3 ≤ m) : 0 < U03.R1 m := by
  have hV := U03.V_gt_one m hm
  have hm' : (3:ℝ) ≤ m := by exact_mod_cast hm
  have hs : ∀ i ∈ Finset.Ico 1 m, 0 < Real.sin ((i:ℝ) * Real.pi / (m:ℝ)) ∧
      Real.sin ((i:ℝ) * Real.pi / (m:ℝ)) ≤ 1 ∧ 
      (i = 1 → Real.sin ((i:ℝ) * Real.pi / (m:ℝ)) < U03.V m) := by
    intro i hi
    rw [Finset.mem_Ico] at hi
    have hi1 : (1:ℝ) ≤ i := by exact_mod_cast hi.1
    have hi2 : (i:ℝ) + 1 ≤ m := by exact_mod_cast hi.2
    refine ⟨?_, Real.sin_le_one _, fun _ => lt_of_le_of_lt (Real.sin_le_one _) hV⟩
    apply Real.sin_pos_of_pos_of_lt_pi
    · positivity
    · rw [div_lt_iff₀ (by linarith)]
      nlinarith [Real.pi_pos]
  unfold U03.R1
  apply mul_pos (by positivity)
  apply Finset.sum_pos'
  · intro i hi
    obtain ⟨h0, h1, _⟩ := hs i hi
    apply mul_nonneg h0.le
    apply Real.log_nonneg
    rw [le_div_iff₀ h0]
    linarith
  · refine ⟨1, Finset.mem_Ico.2 ⟨le_refl _, by omega⟩, ?_⟩
    obtain ⟨h0, h1, h2⟩ := hs 1 (Finset.mem_Ico.2 ⟨le_refl _, by omega⟩)
    apply mul_pos h0
    apply Real.log_pos
    rw [lt_div_iff₀ h0]
    linarith

lemma U03.R1_two : U03.R1 2 = 0 := by
  have : U03.V 2 = 1 := by
    unfold U03.V
    have : Real.pi / (2 * ((2:ℕ):ℝ)) = Real.pi / 4 := by norm_num
    rw [this, Real.cot_eq_cos_div_sin, Real.cos_pi_div_four, Real.sin_pi_div_four]
    have : Real.sqrt 2 / 2 ≠ 0 := by positivity
    exact div_self this
  unfold U03.R1
  rw [this]
  simp [Finset.sum_singleton]
  
lemma U03.P1_pos (m : ℕ) (hm : 2 ≤ m) : 0 < U03.P1 m := by
  have hm' : (2:ℝ) ≤ m := by exact_mod_cast hm
  set x := Real.pi / (2 * (m:ℝ)) with hx
  have hx0 : 0 < x := by positivity
  have hx1 : x < Real.pi / 2 := by
    rw [hx, div_lt_div_iff₀ (by positivity) (by norm_num)]
    nlinarith [Real.pi_pos]
  have hc : 0 < Real.cos x := Real.cos_pos_of_mem_Ioo ⟨by linarith, hx1⟩
  have hs : 0 < Real.sin x := Real.sin_pos_of_pos_of_lt_pi hx0 (by linarith [Real.pi_pos])
  have hc1 : Real.cos x < 1 := by
    have := Real.sin_sq_add_cos_sq x
    nlinarith
  unfold U03.P1
  rw [← hx]
  apply mul_pos (by positivity)
  rw [sub_pos, lt_div_iff₀ hc]
  nlinarith

theorem U03.proof_1 : U03.stmt_1 := fun m hm => U03.R1_pos m hm

theorem U03.proof_2 : U03.stmt_2 := U03.R1_two

theorem U03.proof_3 : U03.stmt_3 := by
  intro m hm
  constructor
  · intro h
    by_contra hne
    have := U03.R1_pos m (by omega)
    linarith
  · rintro rfl; exact U03.R1_two

theorem U03.proof_4 : U03.stmt_4 := fun m hm => U03.P1_pos m hm

theorem U03.proof_5 : U03.stmt_5 := by
  intro m hm
  have hp := U03.P1_pos m hm
  rw [mul_eq_zero]
  constructor
  · rintro (h | h)
    · linarith
    · by_contra hne
      have := U03.R1_pos m (by omega)
      linarith
  · rintro rfl; exact Or.inr U03.R1_two

theorem U03.proof_6 : U03.stmt_6 := fun m hm =>
  mul_pos (U03.P1_pos m (by omega)) (U03.R1_pos m hm)
