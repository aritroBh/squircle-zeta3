import Mathlib
/-!
# T06 — the ratio `R_1` at the regular polygon (paper label `thm:poly`)

Source: `paper/ieee/main.tex`, Subsection "Polygons at the Diamond End", lines 1645–1653
(detailed proof: `paper/ieee/proofs/sec_g4_poly.tex`, label `g4:thm:poly`).

Setup quoted verbatim (line 1645):
"With $n_i=(\cos\frac{i\pi}m,\sin\frac{i\pi}m)$, let $K_p=\{\sum_{i=1}^m|\langle n_i,x\rangle|^p\le1\}$.
At $p=1$ this is a regular $2m$-gon of perimeter $P(1)=4m(\sec\frac{\pi}{2m}-\cos\frac{\pi}{2m})$,
and $m=2$ is the squircle."

Theorem quoted verbatim (lines 1647–1653):
"\begin{theorem}[First derivative at the regular polygon$^\dagger$]\label{thm:poly}
$P'(1)=P(1)\,R_1$, where $V=\cot\frac{\pi}{2m}$ and
\[
R_1=\frac1V\sum_{i=1}^{m-1}\sin\frac{i\pi}{m}\,\log\frac{V}{\sin(i\pi/m)} .
\]
In particular $P'(1)=0$ only for $m=2$, and for the hexagon ($m=3$) $P'(1)=2\sqrt3\log2$.
\end{theorem}"

Also used (line 1686): "The regular hexagon has $P(1)=2\sqrt3$ ...".

What is formalized: `V`, `R1` and `P1` are defined exactly by the paper's formulas
(`\sec` written as `1 / cos`; all divisions are by strictly positive reals for the `m` used,
so no `x/0` or `log` of a nonpositive number is involved). The statements are the
closed-form evaluations at `m = 2` and `m = 3`, and the hexagon value of `P(1) R_1`, i.e. the
value `2√3 log 2` that the theorem asserts for `P'(1)` once `P'(1) = P(1) R_1` is granted.

NOT COVERED:
* The main identity `P'(1) = P(1) R_1` itself: it needs the perimeter of the convex body
  `K_p` as a function of `p` and its derivative at `p = 1`; Mathlib has no convenient
  perimeter (boundary length) of a planar convex body, so this is not stated.
* That the formula `4m(sec(π/2m) − cos(π/2m))` is the perimeter of `K_1` (geometric fact).
  Here `P1` is the paper's stated formula only.
* "P'(1)=0 only for m=2" (excluded by instruction).
-/
-- BEGIN STATEMENTS
namespace T06

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

/-- Theorem `thm:poly`, squircle case: `R_1 = 0` for `m = 2` (so `P'(1) = P(1) R_1 = 0`). -/
def stmt_1 : Prop := R1 2 = 0

/-- Theorem `thm:poly`, hexagon case: `R_1 = log 2` for `m = 3`. -/
def stmt_2 : Prop := R1 3 = Real.log 2

/-- Line 1686: the regular hexagon (`m = 3`) has `P(1) = 2√3`, from the formula of line 1645. -/
def stmt_3 : Prop := P1 3 = 2 * Real.sqrt 3

/-- Theorem `thm:poly`: "for the hexagon ($m=3$) $P'(1)=2\sqrt3\log2$", i.e. `P(1) R_1 = 2√3 log 2`. -/
def stmt_4 : Prop := P1 3 * R1 3 = 2 * Real.sqrt 3 * Real.log 2

end T06
-- END STATEMENTS
open Real in
lemma T06.V2 : T06.V 2 = 1 := by
  unfold T06.V
  have : π / (2 * ((2:ℕ):ℝ)) = π / 4 := by push_cast; ring
  rw [this, Real.cot_eq_cos_div_sin, Real.cos_pi_div_four, Real.sin_pi_div_four]
  have : (0:ℝ) < √2 / 2 := by positivity
  field_simp

open Real in
lemma T06.V3 : T06.V 3 = √3 := by
  unfold T06.V
  have : π / (2 * ((3:ℕ):ℝ)) = π / 6 := by push_cast; ring
  rw [this, Real.cot_eq_cos_div_sin, Real.cos_pi_div_six, Real.sin_pi_div_six]
  have : (0:ℝ) < √3 := by positivity
  field_simp

open Real in
lemma T06.sin3 : Real.sin (π / 3) = √3 / 2 := Real.sin_pi_div_three

open Real in
lemma T06.sin3' : Real.sin (2 * π / 3) = √3 / 2 := by
  have : 2 * π / 3 = π - π / 3 := by ring
  rw [this, Real.sin_pi_sub, Real.sin_pi_div_three]

open Real in
theorem T06.proof_1 : T06.stmt_1 := by
  unfold T06.stmt_1 T06.R1
  rw [T06.V2]
  have : Finset.Ico 1 2 = {1} := by decide
  rw [this]
  have h : (1:ℝ) * π / ((2:ℕ):ℝ) = π / 2 := by push_cast; ring
  simp

open Real in
theorem T06.proof_2 : T06.stmt_2 := by
  unfold T06.stmt_2 T06.R1
  rw [T06.V3]
  have : Finset.Ico 1 3 = {1, 2} := by decide
  rw [this]
  have h1 : ((1:ℕ):ℝ) * π / ((3:ℕ):ℝ) = π / 3 := by push_cast; ring
  have h2 : ((2:ℕ):ℝ) * π / ((3:ℕ):ℝ) = 2 * π / 3 := by push_cast; ring
  rw [Finset.sum_pair (by decide), h1, h2, T06.sin3, T06.sin3']
  have h3 : (0:ℝ) < √3 := by positivity
  have h4 : √3 / (√3 / 2) = 2 := by field_simp
  rw [h4]
  field_simp
  ring

open Real in
theorem T06.proof_3 : T06.stmt_3 := by
  unfold T06.stmt_3 T06.P1
  have : π / (2 * ((3:ℕ):ℝ)) = π / 6 := by push_cast; ring
  rw [this, Real.cos_pi_div_six]
  have h3 : (0:ℝ) < √3 := by positivity
  have h4 : √3 * √3 = 3 := Real.mul_self_sqrt (by norm_num)
  push_cast
  field_simp
  nlinarith [h4]

theorem T06.proof_4 : T06.stmt_4 := by
  unfold T06.stmt_4
  have h3 := T06.proof_3
  have h2 := T06.proof_2
  unfold T06.stmt_3 at h3
  unfold T06.stmt_2 at h2
  rw [h3, h2]
