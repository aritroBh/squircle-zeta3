import Mathlib
/-!
# W04 -- the perimeter integral `eq:Pint` and its three checks

Target id: W04.
Paper labels: `eq:Pint` and the paragraph around it, "An Integral for the Perimeter"
(paper/ieee/main.tex, lines 823-829), with the squircle defined at main.tex lines 785-788 and the
landmarks listed at main.tex lines 802-805.

Quoted verbatim (main.tex 785-788):

  Let's deform a shape. The shape is the squircle
  \[
  |x|^p+|y|^p=1,\qquad p\ge1 .
  \]

Quoted verbatim (main.tex 802-805):

  The perimeter $P(p)$ has no such closed form. We only know the landmarks
  \[
  P(1)=4\sqrt2,\qquad P(2)=2\pi,\qquad \lim_{p\to\infty}P(p)=8 ,
  \]

Quoted verbatim (main.tex 824-829):

  Following Ahmed~\cite{ahmed}, parametrize the quarter arc in the first quadrant by
  $x=t^{1/p}$, $y=(1-t)^{1/p}$, $0\le t\le1$. Then $x^p+y^p=1$, and the arc-length formula gives
  \begin{equation}\label{eq:Pint}
  P(p)=\frac4p\int_0^1\sqrt{t^{2/p-2}+(1-t)^{2/p-2}}\;dt .
  \end{equation}
  Let's check the shapes we know. \textbf{Case 1:} $p=1$. The integrand is $\sqrt2$, so
  $P(1)=4\sqrt2$: four sides of length $\sqrt2$. \textbf{Case 2:} $p=2$. Here
  $t^{-1}+(1-t)^{-1}=1/(t(1-t))$, so $P(2)=2\int_0^1dt/\sqrt{t(1-t)}=2B(\frac12,\frac12)=2\pi$.
  \textbf{Case 3:} $p\to\infty$ gives $P\to8$; exact hypergeometric series are in~\cite{MRL,RM}.

## Formalization

Points of the plane are `EuclideanSpace ℝ (Fin 2)` (`Plane`), so `dist`/`edist`/`‖·‖` are the
Euclidean ones (NOT the sup metric of `ℝ × ℝ`); `v 0`, `v 1` are the coordinates. Lengths are
total variations `eVariationOn` (supremum over partitions of sums of Euclidean distances), so
"length" is the true arc length, with no smoothness assumed. Every statement is for real `p ≥ 1`,
the paper's range; `p = 1` is included throughout (there the quarter arc is the segment
`(t, 1-t)` and the integrand is the constant `√2`, not singular).

* `pIntegrand`, `Pint` are literally `W02.pIntegrand`, `W02.Pint` (Zeta3/Wave3/W02_Holo.lean):
  `Pint p = 4 / p * ∫ t in 0..1, √(t^(2/p-2) + (1-t)^(2/p-2))`, i.e. the right side of eq:Pint.
  For `p > 1` the integrand blows up like `t^(1/p-1)` at the endpoints; its integrability is a
  separate statement (stmt_3) so that no claim rests on Lean's "non-integrable ↦ 0" convention.
  `Real.rpow` gives `0^x = 0` for `x ≠ 0` and `0^0 = 1`; this only affects the null set `{0,1}`.
* `quarterArc p t = (t^(1/p), (1-t)^(1/p))`, the paper's parametrization (used on `[0,1]`, where
  the bases are nonnegative).
* `squircle p = {v | |v 0|^p + |v 1|^p = 1}`, `quarterSet p` = its points with `v 0, v 1 ≥ 0`.
* `loop p` on `[0,4]`: the quarter arc on `[0,1]` followed by its reflections in the other three
  quadrants (same layout as `T14.loop`, but valued in the Euclidean plane).

Statements:
* stmt_1: `quarterArc p` is a bijection of `[0,1]` onto the closed first-quadrant part of the
  curve ("Then $x^p+y^p=1$").
* stmt_2: on `(0,1)` the map is differentiable with speed `‖γ'(t)‖ = (1/p)·√(t^(2/p-2)+(1-t)^(2/p-2))`
  (the integrand of the arc-length formula).
* stmt_3: the integrand of eq:Pint is interval integrable on `[0,1]`.
* stmt_4: the length of the quarter arc is `(1/p) ∫_0^1 √(t^(2/p-2)+(1-t)^(2/p-2)) dt`
  (in particular finite).
* stmt_5: `loop p` is a simple closed parametrization of the whole curve `|x|^p+|y|^p=1`
  (continuous on `[0,4]`, closed, injective on `[0,4)`, image = the curve); so the four
  reflected quarter arcs cover the curve and meet only at their shared endpoints, the four axis
  points `(±1,0), (0,±1)`.
* stmt_6: symmetry: the perimeter (length of `loop p` on `[0,4]`) is four times the length of
  the quarter arc.
* stmt_7: eq:Pint: the perimeter equals `Pint p`.
* stmt_8: Case 1: `Pint 1 = 4√2`.
* stmt_9: Case 2 intermediate: `∫_0^1 dt/√(t(1-t)) = π` (= `B(1/2,1/2)`), with integrability.
* stmt_10: Case 2: `Pint 2 = 2π`.
* stmt_11: Case 3: `Pint p → 8` as `p → ∞`.

## NOT COVERED

* The citation remark "exact hypergeometric series are in~\cite{MRL,RM}" (a reference, not a
  claim of this paper).
* The identification `B(1/2,1/2)` as Euler's Beta function is not stated through Mathlib's
  `Complex.betaIntegral`; stmt_9 states its value `π` directly as the real integral.
* "Perimeter" is formalized as the total variation of the explicit loop `loop p`; no
  parametrization-free notion (e.g. the 1-dimensional Hausdorff measure of the curve) is used.
-/
-- BEGIN STATEMENTS
namespace W04

open MeasureTheory Set Filter

/-- The Euclidean plane. -/
abbrev Plane := EuclideanSpace ℝ (Fin 2)

/-- The integrand of eq:Pint: `√(t^(2/p-2) + (1-t)^(2/p-2))` (literally `W02.pIntegrand`). -/
noncomputable def pIntegrand (p t : ℝ) : ℝ :=
  Real.sqrt (t ^ (2 / p - 2) + (1 - t) ^ (2 / p - 2))

/-- eq:Pint: `P(p) = (4/p) ∫_0^1 √(t^(2/p-2) + (1-t)^(2/p-2)) dt` (literally `W02.Pint`). -/
noncomputable def Pint (p : ℝ) : ℝ := 4 / p * ∫ t in (0 : ℝ)..1, pIntegrand p t

/-- The squircle `|x|^p + |y|^p = 1`. -/
def squircle (p : ℝ) : Set Plane := {v | |v 0| ^ p + |v 1| ^ p = 1}

/-- The part of the squircle in the closed first quadrant. -/
def quarterSet (p : ℝ) : Set Plane := {v | v ∈ squircle p ∧ 0 ≤ v 0 ∧ 0 ≤ v 1}

/-- The paper's quarter arc `x = t^{1/p}`, `y = (1-t)^{1/p}`, `0 ≤ t ≤ 1`. -/
noncomputable def quarterArc (p t : ℝ) : Plane := !₂[t ^ (1 / p), (1 - t) ^ (1 / p)]

/-- The whole curve traversed once on `[0,4]`: the quarter arc on `[0,1]` (from `(0,1)` to
`(1,0)`), then its reflections into the fourth, third and second quadrants. -/
noncomputable def loop (p t : ℝ) : Plane :=
  if t ≤ 1 then !₂[t ^ (1 / p), (1 - t) ^ (1 / p)]
  else if t ≤ 2 then !₂[(2 - t) ^ (1 / p), -((t - 1) ^ (1 / p))]
  else if t ≤ 3 then !₂[-((t - 2) ^ (1 / p)), -((3 - t) ^ (1 / p))]
  else !₂[-((4 - t) ^ (1 / p)), (t - 3) ^ (1 / p)]

/-- The perimeter of the squircle: the Euclidean length of `loop p` over `[0,4]`. -/
noncomputable def perimeter (p : ℝ) : ENNReal := eVariationOn (loop p) (Icc (0 : ℝ) 4)

/-- "parametrize the quarter arc in the first quadrant by $x=t^{1/p}$, $y=(1-t)^{1/p}$,
$0\le t\le1$. Then $x^p+y^p=1$": a bijection of `[0,1]` onto the first-quadrant arc. -/
def stmt_1 : Prop :=
  ∀ p : ℝ, 1 ≤ p → BijOn (quarterArc p) (Icc 0 1) (quarterSet p)

/-- "the arc-length formula": on `(0,1)` the speed of the quarter arc is
`(1/p) √(t^(2/p-2)+(1-t)^(2/p-2))`. -/
def stmt_2 : Prop :=
  ∀ p : ℝ, 1 ≤ p → ∀ t ∈ Ioo (0 : ℝ) 1,
    ∃ v : Plane, HasDerivAt (quarterArc p) v t ∧ ‖v‖ = 1 / p * pIntegrand p t

/-- The integrand of eq:Pint (singular at the endpoints for `p > 1`) is integrable on `[0,1]`. -/
def stmt_3 : Prop :=
  ∀ p : ℝ, 1 ≤ p → IntervalIntegrable (pIntegrand p) volume 0 1

/-- The arc-length formula: the (finite) length of the quarter arc is
`(1/p) ∫_0^1 √(t^(2/p-2)+(1-t)^(2/p-2)) dt`. -/
def stmt_4 : Prop :=
  ∀ p : ℝ, 1 ≤ p →
    eVariationOn (quarterArc p) (Icc 0 1) =
      ENNReal.ofReal (1 / p * ∫ t in (0 : ℝ)..1, pIntegrand p t)

/-- The whole curve is the union of the quarter arc and its three reflections, traversed once:
`loop p` is a simple closed parametrization of `|x|^p+|y|^p=1`. -/
def stmt_5 : Prop :=
  ∀ p : ℝ, 1 ≤ p →
    ContinuousOn (loop p) (Icc 0 4) ∧ loop p 4 = loop p 0 ∧ InjOn (loop p) (Ico 0 4) ∧
    loop p '' Icc 0 4 = squircle p

/-- Symmetry: the perimeter is four times the length of the quarter arc. -/
def stmt_6 : Prop :=
  ∀ p : ℝ, 1 ≤ p → perimeter p = 4 * eVariationOn (quarterArc p) (Icc 0 1)

/-- eq:Pint: `P(p)=\frac4p\int_0^1\sqrt{t^{2/p-2}+(1-t)^{2/p-2}}\;dt` for the true perimeter. -/
def stmt_7 : Prop :=
  ∀ p : ℝ, 1 ≤ p → perimeter p = ENNReal.ofReal (Pint p)

/-- Case 1: "The integrand is $\sqrt2$, so $P(1)=4\sqrt2$". -/
def stmt_8 : Prop :=
  Pint 1 = 4 * Real.sqrt 2

/-- Case 2 (intermediate): `\int_0^1dt/\sqrt{t(1-t)}=B(\frac12,\frac12)=\pi`. -/
def stmt_9 : Prop :=
  IntervalIntegrable (fun t : ℝ => 1 / Real.sqrt (t * (1 - t))) volume 0 1 ∧
    ∫ t in (0 : ℝ)..1, 1 / Real.sqrt (t * (1 - t)) = Real.pi

/-- Case 2: "$P(2)=2\int_0^1dt/\sqrt{t(1-t)}=2B(\frac12,\frac12)=2\pi$". -/
def stmt_10 : Prop :=
  Pint 2 = 2 * Real.pi

/-- Case 3: "$p\to\infty$ gives $P\to8$". -/
def stmt_11 : Prop :=
  Tendsto Pint atTop (nhds 8)

end W04
-- END STATEMENTS
open MeasureTheory Set Filter
open scoped Topology

namespace W04


theorem pw1 {p : ℝ} (hp : 1 ≤ p) {x : ℝ} (hx : 0 ≤ x) : (x ^ (1 / p)) ^ p = x := by
  rw [← Real.rpow_mul hx, one_div_mul_cancel (by linarith), Real.rpow_one]

theorem pw2 {p : ℝ} (hp : 1 ≤ p) {x : ℝ} (hx : 0 ≤ x) : (x ^ p) ^ (1 / p) = x := by
  rw [← Real.rpow_mul hx, mul_one_div_cancel (by linarith), Real.rpow_one]

@[simp] theorem vec0 (a b : ℝ) : (!₂[a, b] : Plane) 0 = a := rfl
@[simp] theorem vec1 (a b : ℝ) : (!₂[a, b] : Plane) 1 = b := rfl

theorem plane_ext {v w : Plane} (h0 : v 0 = w 0) (h1 : v 1 = w 1) : v = w := by
  ext i; fin_cases i <;> simp [h0, h1]

theorem plane_eta (v : Plane) : (!₂[v 0, v 1] : Plane) = v := plane_ext rfl rfl

theorem norm_vec (a b : ℝ) : ‖(!₂[a, b] : Plane)‖ = Real.sqrt (a ^ 2 + b ^ 2) := by
  rw [EuclideanSpace.norm_eq, Fin.sum_univ_two]; simp

theorem dist_vec (a b c d : ℝ) :
    dist (!₂[a, b] : Plane) !₂[c, d] = Real.sqrt ((a - c) ^ 2 + (b - d) ^ 2) := by
  rw [EuclideanSpace.dist_eq, Fin.sum_univ_two]; simp [Real.dist_eq, sq_abs]

theorem dist_plane (v w : Plane) : dist v w = Real.sqrt ((v 0 - w 0) ^ 2 + (v 1 - w 1) ^ 2) := by
  rw [← dist_vec, plane_eta, plane_eta]

/-- Embedding of `ℝ × ℝ` into the plane. -/
noncomputable def emb (z : ℝ × ℝ) : Plane := !₂[z.1, z.2]

theorem emb_inj : Function.Injective emb := by
  intro z w h
  have h0 := congrArg (fun v : Plane => v 0) h
  have h1 := congrArg (fun v : Plane => v 1) h
  simp [emb] at h0 h1
  exact Prod.ext h0 h1

theorem emb_cont : Continuous emb := by
  unfold emb
  exact (PiLp.continuousLinearEquiv 2 ℝ (fun _ : Fin 2 => ℝ)).symm.continuous.comp
    (continuous_pi fun i => by fin_cases i <;> simp <;> fun_prop)

noncomputable def loopP (p : ℝ) (t : ℝ) : ℝ × ℝ :=
  if t ≤ 1 then (t ^ (1 / p), (1 - t) ^ (1 / p))
  else if t ≤ 2 then ((2 - t) ^ (1 / p), -((t - 1) ^ (1 / p)))
  else if t ≤ 3 then (-((t - 2) ^ (1 / p)), -((3 - t) ^ (1 / p)))
  else (-((4 - t) ^ (1 / p)), (t - 3) ^ (1 / p))

theorem loop_eq (p t : ℝ) : loop p t = emb (loopP p t) := by
  unfold loop loopP emb; split_ifs <;> rfl

noncomputable def arc (p s : ℝ) : ℝ × ℝ := (s ^ (1 / p), (1 - s) ^ (1 / p))

theorem qa_eq (p s : ℝ) : quarterArc p s = emb (arc p s) := rfl

theorem loop_Q1 {p : ℝ} (t : ℝ) (ht : t ∈ Icc (0:ℝ) 1) :
    loopP p t = arc p (t - 0) := by
  simp [loopP, arc, ht.2]

theorem loop_Q2 {p : ℝ} (hp : 1 ≤ p) (t : ℝ) (ht : t ∈ Icc (1:ℝ) 2) :
    loopP p t = ((arc p (t - 1)).2, -(arc p (t - 1)).1) := by
  have hr : (1 / p) ≠ 0 := by positivity
  have h0 : (0:ℝ) ^ p⁻¹ = 0 := by simpa using Real.zero_rpow hr
  rcases ht.1.eq_or_lt with rfl | h1
  · simp [loopP, arc, h0]
  · have h : ¬ t ≤ 1 := not_le.2 h1
    simp only [loopP, arc, h, ite_false, ht.2, ite_true]
    rw [show 1 - (t - 1) = 2 - t by ring]

theorem loop_Q3 {p : ℝ} (hp : 1 ≤ p) (t : ℝ) (ht : t ∈ Icc (2:ℝ) 3) :
    loopP p t = (-(arc p (t - 2)).1, -(arc p (t - 2)).2) := by
  have hr : (1 / p) ≠ 0 := by positivity
  have h0 : (0:ℝ) ^ p⁻¹ = 0 := by simpa using Real.zero_rpow hr
  have h : ¬ t ≤ 1 := by have := ht.1; intro; linarith
  rcases ht.1.eq_or_lt with rfl | h1
  · norm_num [loopP, arc, h0]
  · have h' : ¬ t ≤ 2 := not_le.2 h1
    simp only [loopP, arc, h, h', ite_false, ht.2, ite_true]
    rw [show 1 - (t - 2) = 3 - t by ring]

theorem loop_Q4 {p : ℝ} (hp : 1 ≤ p) (t : ℝ) (ht : t ∈ Icc (3:ℝ) 4) :
    loopP p t = (-(arc p (t - 3)).2, (arc p (t - 3)).1) := by
  have hr : (1 / p) ≠ 0 := by positivity
  have h0 : (0:ℝ) ^ p⁻¹ = 0 := by simpa using Real.zero_rpow hr
  have h : ¬ t ≤ 1 := by have := ht.1; intro; linarith
  have h' : ¬ t ≤ 2 := by have := ht.1; intro; linarith
  rcases ht.1.eq_or_lt with rfl | h1
  · norm_num [loopP, arc, h0]
  · have h'' : ¬ t ≤ 3 := not_le.2 h1
    simp only [loopP, arc, h, h', h'', ite_false]
    rw [show 1 - (t - 3) = 4 - t by ring]

theorem arc_cont {p : ℝ} (hp : 1 ≤ p) : Continuous (arc p) := by
  have hr : 0 ≤ 1 / p := by positivity
  exact ((Real.continuous_rpow_const hr).comp continuous_id).prodMk
    ((Real.continuous_rpow_const hr).comp (continuous_const.sub continuous_id))

theorem loopP_cont {p : ℝ} (hp : 1 ≤ p) : ContinuousOn (loopP p) (Icc 0 4) := by
  have hc := arc_cont hp
  have c1 : ContinuousOn (loopP p) (Icc 0 1) :=
    (hc.comp (continuous_id.sub continuous_const) |>.continuousOn).congr
      (fun t ht => loop_Q1 t ht)
  have c2 : ContinuousOn (loopP p) (Icc 1 2) :=
    ((hc.comp (continuous_id.sub continuous_const)).snd.prodMk
      (hc.comp (continuous_id.sub continuous_const)).fst.neg |>.continuousOn).congr
      (fun t ht => loop_Q2 hp t ht)
  have c3 : ContinuousOn (loopP p) (Icc 2 3) :=
    ((hc.comp (continuous_id.sub continuous_const)).fst.neg.prodMk
      (hc.comp (continuous_id.sub continuous_const)).snd.neg |>.continuousOn).congr
      (fun t ht => loop_Q3 hp t ht)
  have c4 : ContinuousOn (loopP p) (Icc 3 4) :=
    ((hc.comp (continuous_id.sub continuous_const)).snd.neg.prodMk
      (hc.comp (continuous_id.sub continuous_const)).fst |>.continuousOn).congr
      (fun t ht => loop_Q4 hp t ht)
  rw [← Icc_union_Icc_eq_Icc (a := 0) (b := 1) (c := 4) (by norm_num) (by norm_num),
    ← Icc_union_Icc_eq_Icc (a := 1) (b := 2) (c := 4) (by norm_num) (by norm_num),
    ← Icc_union_Icc_eq_Icc (a := 2) (b := 3) (c := 4) (by norm_num) (by norm_num)]
  exact c1.union_of_isClosed (c2.union_of_isClosed (c3.union_of_isClosed c4 isClosed_Icc isClosed_Icc)
    isClosed_Icc (isClosed_Icc.union isClosed_Icc)) isClosed_Icc
    (isClosed_Icc.union (isClosed_Icc.union isClosed_Icc))

noncomputable def ginv (p : ℝ) (z : ℝ × ℝ) : ℝ :=
  if 0 ≤ z.1 ∧ 0 < z.2 then |z.1| ^ p
  else if 0 < z.1 ∧ z.2 ≤ 0 then 1 + |z.2| ^ p
  else if z.1 ≤ 0 ∧ z.2 < 0 then 2 + |z.1| ^ p
  else 3 + |z.2| ^ p

theorem ginv_loop {p : ℝ} (hp : 1 ≤ p) (t : ℝ) (ht : t ∈ Ico (0:ℝ) 4) :
    ginv p (loopP p t) = t := by
  have hr : 0 < 1 / p := by positivity
  have nn : ∀ x : ℝ, 0 ≤ x → 0 ≤ x ^ (1 / p) := fun x hx => Real.rpow_nonneg hx _
  have ps : ∀ x : ℝ, 0 < x → 0 < x ^ (1 / p) := fun x hx => Real.rpow_pos_of_pos hx _
  have ht2 := ht.2
  by_cases h1 : t < 1
  · have hz := loop_Q1 (p := p) t ⟨ht.1, h1.le⟩
    rw [hz]
    simp only [arc, sub_zero]
    have a1 := nn t ht.1
    have a2 := ps (1 - t) (by linarith)
    rw [ginv, if_pos ⟨a1, a2⟩, abs_of_nonneg a1, pw1 hp ht.1]
  by_cases h2 : t < 2
  · have hz := loop_Q2 hp t ⟨not_lt.1 h1, h2.le⟩
    rw [hz]
    simp only [arc]
    have a1 := nn (t - 1) (by linarith)
    have a2 := ps (1 - (t - 1)) (by linarith)
    rw [ginv, if_neg (by intro h; linarith [h.2]), if_pos ⟨a2, by linarith⟩, abs_neg,
      abs_of_nonneg a1, pw1 hp (by linarith)]
    ring
  by_cases h3 : t < 3
  · have hz := loop_Q3 hp t ⟨not_lt.1 h2, h3.le⟩
    rw [hz]
    simp only [arc]
    have a1 := nn (t - 2) (by linarith)
    have a2 := ps (1 - (t - 2)) (by linarith)
    rw [ginv, if_neg (by intro h; linarith [h.2]), if_neg (by intro h; linarith [h.1]),
      if_pos ⟨by linarith, by linarith⟩, abs_neg, abs_of_nonneg a1, pw1 hp (by linarith)]
    ring
  · have hz := loop_Q4 hp t ⟨not_lt.1 h3, ht.2.le⟩
    rw [hz]
    simp only [arc]
    have a1 := nn (t - 3) (by linarith)
    have a2 := ps (1 - (t - 3)) (by linarith)
    rw [ginv, if_neg (by intro h; linarith [h.1]), if_neg (by intro h; linarith [h.1]),
      if_neg (by intro h; linarith [h.2]), abs_of_nonneg a1, pw1 hp (by linarith)]
    ring

theorem norm_loop {p : ℝ} (hp : 1 ≤ p) (t : ℝ) (ht : t ∈ Icc (0:ℝ) 4) :
    |(loopP p t).1| ^ p + |(loopP p t).2| ^ p = 1 := by
  have nn : ∀ x : ℝ, 0 ≤ x → 0 ≤ x ^ (1 / p) := fun x hx => Real.rpow_nonneg hx _
  have key : ∀ s : ℝ, 0 ≤ s → s ≤ 1 →
      |(arc p s).1| ^ p + |(arc p s).2| ^ p = 1 := by
    intro s h0 h1
    simp only [arc]
    rw [abs_of_nonneg (nn _ h0), abs_of_nonneg (nn _ (by linarith)), pw1 hp h0,
      pw1 hp (by linarith)]
    ring
  by_cases h1 : t ≤ 1
  · rw [loop_Q1 t ⟨ht.1, h1⟩]; exact key _ (by linarith [ht.1]) (by linarith)
  by_cases h2 : t ≤ 2
  · rw [loop_Q2 hp t ⟨(not_le.1 h1).le, h2⟩]
    have := key (t - 1) (by linarith [not_le.1 h1]) (by linarith)
    simpa [abs_neg, add_comm] using this
  by_cases h3 : t ≤ 3
  · rw [loop_Q3 hp t ⟨(not_le.1 h2).le, h3⟩]
    have := key (t - 2) (by linarith [not_le.1 h2]) (by linarith)
    simpa [abs_neg] using this
  · rw [loop_Q4 hp t ⟨(not_le.1 h3).le, ht.2⟩]
    have := key (t - 3) (by linarith [not_le.1 h3]) (by linarith [ht.2])
    simpa [abs_neg, add_comm] using this

theorem exists_t {p : ℝ} (hp : 1 ≤ p) (z : ℝ × ℝ) (hz : |z.1| ^ p + |z.2| ^ p = 1) :
    ∃ t ∈ Icc (0:ℝ) 4, loopP p t = z := by
  have hp0 : 0 < p := by linarith
  set u := |z.1| ^ p with hu
  set v := |z.2| ^ p with hv
  have u0 : 0 ≤ u := Real.rpow_nonneg (abs_nonneg _) _
  have v0 : 0 ≤ v := Real.rpow_nonneg (abs_nonneg _) _
  have ur : u ^ (1 / p) = |z.1| := pw2 hp (abs_nonneg _)
  have vr : v ^ (1 / p) = |z.2| := pw2 hp (abs_nonneg _)
  have u1 : (1 - u) = v := by linarith
  have v1 : (1 - v) = u := by linarith
  by_cases c1 : 0 ≤ z.1 ∧ 0 < z.2
  · refine ⟨u, ⟨u0, by linarith⟩, ?_⟩
    rw [loop_Q1 u ⟨u0, by linarith⟩]
    simp only [arc, sub_zero, u1, ur, vr]
    ext <;> simp [abs_of_nonneg c1.1, abs_of_pos c1.2]
  by_cases c2 : 0 < z.1 ∧ z.2 ≤ 0
  · refine ⟨1 + v, ⟨by linarith, by linarith⟩, ?_⟩
    rw [loop_Q2 hp (1 + v) ⟨by linarith, by linarith⟩]
    simp only [arc, add_sub_cancel_left, v1, ur, vr]
    ext <;> simp [abs_of_pos c2.1, abs_of_nonpos c2.2]
  by_cases c3 : z.1 ≤ 0 ∧ z.2 < 0
  · refine ⟨2 + u, ⟨by linarith, by linarith⟩, ?_⟩
    rw [loop_Q3 hp (2 + u) ⟨by linarith, by linarith⟩]
    simp only [arc, add_sub_cancel_left, u1, ur, vr]
    ext <;> simp [abs_of_nonpos c3.1, abs_of_neg c3.2]
  · have h1 : z.1 ≤ 0 := by
      by_contra h
      have hx := not_le.1 h
      have : 0 < z.2 := by
        by_contra h'
        exact c2 ⟨hx, not_lt.1 h'⟩
      exact c1 ⟨hx.le, this⟩
    have h2 : 0 ≤ z.2 := by
      by_contra h
      exact c3 ⟨h1, not_le.1 h⟩
    refine ⟨3 + v, ⟨by linarith, by linarith⟩, ?_⟩
    rw [loop_Q4 hp (3 + v) ⟨by linarith, by linarith⟩]
    simp only [arc, add_sub_cancel_left, v1, ur, vr]
    ext <;> simp [abs_of_nonpos h1, abs_of_nonneg h2]

/-! ### variation helpers -/

theorem var_congr {E F : Type*} [PseudoEMetricSpace E] [PseudoEMetricSpace F]
    (f : ℝ → E) (h : ℝ → F) (s : Set ℝ)
    (hd : ∀ a ∈ s, ∀ b ∈ s, edist (f a) (f b) = edist (h a) (h b)) :
    eVariationOn f s = eVariationOn h s := by
  unfold eVariationOn
  refine iSup_congr fun p => Finset.sum_congr rfl fun i _ => ?_
  exact hd _ (p.2.2.2 _) _ (p.2.2.2 _)

theorem var_shift {E F : Type*} [PseudoEMetricSpace E] [PseudoEMetricSpace F]
    (f : ℝ → E) (g : ℝ → F) (k : ℝ)
    (hd : ∀ a ∈ Icc k (k+1), ∀ b ∈ Icc k (k+1), edist (f a) (f b) = edist (g (a-k)) (g (b-k))) :
    eVariationOn f (Icc k (k+1)) = eVariationOn g (Icc 0 1) := by
  rw [var_congr f (g ∘ fun t => t - k) _ hd,
    eVariationOn.comp_eq_of_monotoneOn g (fun t => t - k) (fun a _ b _ h => by simpa using h)]
  congr 1
  ext x
  simp only [mem_image, mem_Icc]
  constructor
  · rintro ⟨y, ⟨h1, h2⟩, rfl⟩; constructor <;> linarith
  · rintro ⟨h1, h2⟩; exact ⟨x + k, ⟨by linarith, by linarith⟩, by ring⟩

theorem var_split {E : Type*} [PseudoEMetricSpace E] (f : ℝ → E) {a b c : ℝ}
    (hab : a ≤ b) (hbc : b ≤ c) :
    eVariationOn f (Icc a b) + eVariationOn f (Icc b c) = eVariationOn f (Icc a c) := by
  have := eVariationOn.Icc_add_Icc f (s := Icc a c) hab hbc ⟨hab, hbc⟩
  rwa [Set.inter_eq_right.2 (Icc_subset_Icc le_rfl hbc),
    Set.inter_eq_right.2 (Icc_subset_Icc hab le_rfl), Set.inter_self] at this

theorem loop_split {E : Type*} [PseudoEMetricSpace E] (f : ℝ → E) :
    eVariationOn f (Icc 0 4) = eVariationOn f (Icc 0 1) + eVariationOn f (Icc 1 2) +
      eVariationOn f (Icc 2 3) + eVariationOn f (Icc 3 4) := by
  rw [← var_split f (b := 1) (by norm_num) (by norm_num),
    ← var_split f (a := 1) (b := 2) (by norm_num) (by norm_num),
    ← var_split f (a := 2) (b := 3) (by norm_num) (by norm_num)]
  ring

theorem piece {p : ℝ} (k m : ℝ) (hm : m = k + 1) (A : ℝ × ℝ → ℝ × ℝ)
    (hA : ∀ z w, edist (emb (A z)) (emb (A w)) = edist (emb z) (emb w))
    (hf : ∀ t ∈ Icc k m, loopP p t = A (arc p (t - k))) :
    eVariationOn (loop p) (Icc k m) = eVariationOn (quarterArc p) (Icc 0 1) := by
  subst hm
  refine var_shift _ _ k fun a ha b hb => ?_
  rw [loop_eq, loop_eq, hf a ha, hf b hb, hA]
  rfl

theorem edist_emb (z w : ℝ × ℝ) :
    edist (emb z) (emb w) = ENNReal.ofReal (Real.sqrt ((z.1 - w.1) ^ 2 + (z.2 - w.2) ^ 2)) := by
  rw [edist_dist, emb, emb, dist_vec]



theorem rpow_sq' {t p : ℝ} (ht : 0 ≤ t) : t ^ (2 / p - 2) = (t ^ (1 / p - 1)) ^ 2 := by
  rw [← Real.rpow_mul_natCast ht]; congr 1; push_cast; ring

theorem sqrt_bounds {A B : ℝ} (hA : 0 ≤ A) (hB : 0 ≤ B) :
    max A B ≤ Real.sqrt (A ^ 2 + B ^ 2) ∧ Real.sqrt (A ^ 2 + B ^ 2) ≤ A + B := by
  constructor
  · apply Real.le_sqrt_of_sq_le
    rcases max_choice A B with h | h <;> rw [h] <;> nlinarith
  · rw [Real.sqrt_le_left (by positivity)]; nlinarith

theorem pInt_eq {p t : ℝ} (ht0 : 0 ≤ t) (ht1 : t ≤ 1) :
    pIntegrand p t = Real.sqrt ((t ^ (1 / p - 1)) ^ 2 + ((1 - t) ^ (1 / p - 1)) ^ 2) := by
  rw [pIntegrand, rpow_sq' ht0, rpow_sq' (by linarith)]

theorem int_rpow {s : ℝ} (hs : -1 < s) :
    IntervalIntegrable (fun x : ℝ => x ^ s) volume 0 1 :=
  intervalIntegral.intervalIntegrable_rpow' hs

theorem int_rpow' {s : ℝ} (hs : -1 < s) :
    IntervalIntegrable (fun x : ℝ => (1 - x) ^ s) volume 0 1 := by
  have := (int_rpow hs).comp_sub_left 1
  simpa using this.symm

theorem int_rpow_val {p : ℝ} (hp : 1 ≤ p) :
    ∫ t in (0:ℝ)..1, t ^ (1 / p - 1) = p := by
  have hr : 0 < 1 / p := by positivity
  rw [integral_rpow (Or.inl (by linarith))]
  rw [show 1 / p - 1 + 1 = 1 / p by ring, Real.zero_rpow hr.ne', Real.one_rpow, sub_zero]
  field_simp

theorem int_rpow_val' {p : ℝ} (hp : 1 ≤ p) :
    ∫ t in (0:ℝ)..1, (1 - t) ^ (1 / p - 1) = p := by
  rw [intervalIntegral.integral_comp_sub_left (fun x : ℝ => x ^ (1 / p - 1)) 1]
  simpa using int_rpow_val hp

theorem int_rpow_half {p : ℝ} (hp : 1 ≤ p) :
    ∫ t in (0:ℝ)..(1 / 2), t ^ (1 / p - 1) = p * (2 : ℝ) ^ (-1 / p) := by
  have hr : 0 < 1 / p := by positivity
  have hc : (1 / 2 : ℝ) ^ (1 / p) = (2 : ℝ) ^ (-1 / p) := by
    rw [one_div, Real.inv_rpow (by norm_num), ← Real.rpow_neg (by norm_num)]
    congr 1; ring
  rw [integral_rpow (Or.inl (by linarith))]
  rw [show 1 / p - 1 + 1 = 1 / p by ring, Real.zero_rpow hr.ne', hc, sub_zero]
  field_simp

theorem pInt_meas (p : ℝ) : Measurable (pIntegrand p) := by
  unfold pIntegrand; fun_prop

/-- Upper/lower bounds of the integrand. -/
theorem pInt_bounds {p t : ℝ} (ht0 : 0 ≤ t) (ht1 : t ≤ 1) :
    max (t ^ (1 / p - 1)) ((1 - t) ^ (1 / p - 1)) ≤ pIntegrand p t ∧
      pIntegrand p t ≤ t ^ (1 / p - 1) + (1 - t) ^ (1 / p - 1) := by
  rw [pInt_eq ht0 ht1]
  exact sqrt_bounds (Real.rpow_nonneg ht0 _) (Real.rpow_nonneg (by linarith) _)

theorem maxInt {p : ℝ} (hp : 1 ≤ p) :
    IntervalIntegrable (fun t => max (t ^ (1 / p - 1)) ((1 - t) ^ (1 / p - 1))) volume 0 1 ∧
    ∫ t in (0 : ℝ)..1, max (t ^ (1 / p - 1)) ((1 - t) ^ (1 / p - 1)) = 2 * p * (2 : ℝ) ^ (-1 / p) := by
  have hr : 0 < 1 / p := by positivity
  set s : ℝ := 1 / p - 1 with hsdef
  have hs : -1 < s := by rw [hsdef]; linarith
  have hr1 : 1 / p ≤ 1 := by rw [div_le_one (by linarith)]; exact hp
  have hs0 : s ≤ 0 := by rw [hsdef]; linarith
  set F : ℝ → ℝ := fun t => max (t ^ s) ((1 - t) ^ s) with hF
  have hFi : IntervalIntegrable F volume 0 1 := by
    refine ((int_rpow hs).add (int_rpow' hs)).mono_fun' ?_ ?_
    · exact (by fun_prop : Measurable F).aestronglyMeasurable
    · refine ae_restrict_of_forall_mem measurableSet_uIoc fun x hx => ?_
      rw [Set.uIoc_of_le zero_le_one] at hx
      have a1 : 0 ≤ x ^ s := Real.rpow_nonneg hx.1.le _
      have a2 : 0 ≤ (1 - x) ^ s := Real.rpow_nonneg (by linarith [hx.2]) _
      simp only [hF, Real.norm_eq_abs]
      rw [abs_of_nonneg (le_max_of_le_left a1)]
      rcases max_choice (x ^ s) ((1 - x) ^ s) with h | h <;> rw [h] <;> linarith
  refine ⟨hFi, ?_⟩
  have sub : ∀ a b : ℝ, 0 ≤ a → a ≤ b → b ≤ 1 → IntervalIntegrable F volume a b := by
    intro a b ha hab hb
    refine hFi.mono_set ?_
    rw [Set.uIcc_of_le zero_le_one, Set.uIcc_of_le hab]
    exact Icc_subset_Icc ha hb
  show ∫ t in (0:ℝ)..1, F t = _
  rw [← intervalIntegral.integral_add_adjacent_intervals (b := 1 / 2)
    (sub 0 (1 / 2) le_rfl (by norm_num) (by norm_num)) (sub (1 / 2) 1 (by norm_num) (by norm_num) le_rfl)]
  have L : ∫ t in (0:ℝ)..(1 / 2), F t = p * (2 : ℝ) ^ (-1 / p) := by
    rw [← int_rpow_half hp]
    refine intervalIntegral.integral_congr_Ioo_of_le (by norm_num) fun t ht => ?_
    simp only [hF]
    apply max_eq_left
    exact Real.rpow_le_rpow_of_nonpos ht.1 (by linarith [ht.2]) hs0
  have R : ∫ t in (1 / 2 : ℝ)..1, F t = p * (2 : ℝ) ^ (-1 / p) := by
    rw [← int_rpow_half hp]
    have : ∫ t in (1 / 2 : ℝ)..1, F t = ∫ t in (1 / 2 : ℝ)..1, (1 - t) ^ s := by
      refine intervalIntegral.integral_congr_Ioo_of_le (by norm_num) fun t ht => ?_
      simp only [hF]
      apply max_eq_right
      exact Real.rpow_le_rpow_of_nonpos (by linarith [ht.2]) (by linarith [ht.1]) hs0
    rw [this, intervalIntegral.integral_comp_sub_left (fun x : ℝ => x ^ s) 1,
      show (1:ℝ) - 1 = 0 by norm_num, show (1:ℝ) - 1 / 2 = 1 / 2 by norm_num]
  rw [L, R]; ring

theorem pInt_iInt {p : ℝ} (hp : 1 ≤ p) : IntervalIntegrable (pIntegrand p) volume 0 1 := by
  have hs : -1 < 1 / p - 1 := by have : 0 < 1 / p := by positivity
                                 linarith
  refine ((int_rpow hs).add (int_rpow' hs)).mono_fun' (pInt_meas p).aestronglyMeasurable ?_
  refine ae_restrict_of_forall_mem measurableSet_uIoc fun x hx => ?_
  rw [Set.uIoc_of_le zero_le_one] at hx
  show ‖pIntegrand p x‖ ≤ _
  rw [Real.norm_eq_abs, abs_of_nonneg (by unfold pIntegrand; exact Real.sqrt_nonneg _)]
  exact (pInt_bounds hx.1.le hx.2).2

theorem Pint_bounds {p : ℝ} (hp : 1 ≤ p) :
    8 * (2 : ℝ) ^ (-1 / p) ≤ Pint p ∧ Pint p ≤ 8 := by
  have hp0 : 0 < p := by linarith
  have hs : -1 < 1 / p - 1 := by have : 0 < 1 / p := by positivity
                                 linarith
  obtain ⟨hmi, hmv⟩ := maxInt hp
  have hsum := (int_rpow hs).add (int_rpow' hs)
  have h1 : ∫ t in (0:ℝ)..1, max (t ^ (1 / p - 1)) ((1 - t) ^ (1 / p - 1)) ≤
      ∫ t in (0:ℝ)..1, pIntegrand p t :=
    intervalIntegral.integral_mono_on zero_le_one hmi (pInt_iInt hp)
      (fun t ht => (pInt_bounds ht.1 ht.2).1)
  have h2 : ∫ t in (0:ℝ)..1, pIntegrand p t ≤
      ∫ t in (0:ℝ)..1, (t ^ (1 / p - 1) + (1 - t) ^ (1 / p - 1)) :=
    intervalIntegral.integral_mono_on zero_le_one (pInt_iInt hp) hsum
      (fun t ht => (pInt_bounds ht.1 ht.2).2)
  rw [intervalIntegral.integral_add (int_rpow hs) (int_rpow' hs), int_rpow_val hp,
    int_rpow_val' hp] at h2
  rw [hmv] at h1
  unfold Pint
  have hq : 0 < 4 / p := by positivity
  constructor
  · calc 8 * (2 : ℝ) ^ (-1 / p) = 4 / p * (2 * p * (2 : ℝ) ^ (-1 / p)) := by field_simp; ring
      _ ≤ _ := mul_le_mul_of_nonneg_left h1 hq.le
  · calc 4 / p * ∫ t in (0:ℝ)..1, pIntegrand p t ≤ 4 / p * (p + p) :=
          mul_le_mul_of_nonneg_left h2 hq.le
      _ = 8 := by field_simp; ring


theorem qa_deriv {p : ℝ} (hp : 1 ≤ p) {t : ℝ} (ht : t ∈ Ioo (0 : ℝ) 1) :
    ∃ v : Plane, HasDerivAt (quarterArc p) v t ∧ ‖v‖ = 1 / p * pIntegrand p t := by
  have hp0 : 0 < p := by linarith
  have ht1 : 0 < 1 - t := by linarith [ht.2]
  set a := 1 / p * t ^ (1 / p - 1)
  set b := -(1 / p * (1 - t) ^ (1 / p - 1))
  refine ⟨!₂[a, b], ?_, ?_⟩
  · have hA : HasDerivAt (fun s : ℝ => s ^ (1 / p)) a t :=
      Real.hasDerivAt_rpow_const (Or.inl ht.1.ne')
    have hB : HasDerivAt (fun s : ℝ => (1 - s) ^ (1 / p)) b t := by
      have := ((hasDerivAt_id' t).const_sub 1).rpow_const (p := 1 / p) (Or.inl ht1.ne')
      convert this using 1
      show -(1 / p * (1 - t) ^ (1 / p - 1)) = _
      ring
    have hpi : HasDerivAt (fun s : ℝ => (![s ^ (1 / p), (1 - s) ^ (1 / p)] : Fin 2 → ℝ)) ![a, b] t := by
      rw [hasDerivAt_pi]
      intro i; fin_cases i
      · simpa using hA
      · simpa using hB
    exact (PiLp.continuousLinearEquiv 2 ℝ (fun _ : Fin 2 => ℝ)).symm.hasFDerivAt.comp_hasDerivAt t hpi
  · rw [norm_vec, pIntegrand]
    have e1 : a ^ 2 = (1 / p) ^ 2 * t ^ (2 / p - 2) := by
      simp only [a]; rw [mul_pow, ← Real.rpow_mul_natCast ht.1.le]; congr 2; push_cast; ring
    have e2 : b ^ 2 = (1 / p) ^ 2 * (1 - t) ^ (2 / p - 2) := by
      simp only [b]; rw [neg_sq, mul_pow, ← Real.rpow_mul_natCast ht1.le]; congr 2; push_cast; ring
    rw [e1, e2, ← mul_add, Real.sqrt_mul (by positivity), Real.sqrt_sq (by positivity)]


theorem qa_cont {p : ℝ} (hp : 1 ≤ p) : Continuous (quarterArc p) := by
  have : quarterArc p = emb ∘ arc p := funext fun s => qa_eq p s
  rw [this]; exact emb_cont.comp (arc_cont hp)

/-- The speed of the quarter arc. -/
noncomputable def spd (p t : ℝ) : ℝ := 1 / p * pIntegrand p t

/-- Arc length from `0` to `x`, as an integral of the speed. -/
noncomputable def Gp (p x : ℝ) : ℝ := ∫ t in (0:ℝ)..x, spd p t

theorem spd_iInt {p : ℝ} (hp : 1 ≤ p) : IntervalIntegrable (spd p) volume 0 1 :=
  (pInt_iInt hp).const_mul _

theorem spd_nonneg {p : ℝ} (_hp : 1 ≤ p) (t : ℝ) : 0 ≤ spd p t := by
  unfold spd pIntegrand
  have : 0 < 1 / p := by positivity
  positivity

theorem spd_sub {p : ℝ} (hp : 1 ≤ p) {a b : ℝ} (ha : 0 ≤ a) (hab : a ≤ b) (hb : b ≤ 1) :
    IntervalIntegrable (spd p) volume a b := by
  refine (spd_iInt hp).mono_set ?_
  rw [uIcc_of_le zero_le_one, uIcc_of_le hab]
  exact Icc_subset_Icc ha hb

theorem Gp_cont {p : ℝ} (hp : 1 ≤ p) : ContinuousOn (Gp p) (Icc 0 1) := by
  have := intervalIntegral.continuousOn_primitive_interval' (spd_iInt hp) left_mem_uIcc
  rw [uIcc_of_le zero_le_one] at this
  exact this

theorem spd_contAt {p : ℝ} (hp : 1 ≤ p) {t : ℝ} (ht : t ∈ Ioo (0:ℝ) 1) :
    ContinuousAt (spd p) t := by
  have h1 : 0 < 1 - t := by linarith [ht.2]
  unfold spd pIntegrand
  have c1 : ContinuousAt (fun s : ℝ => s ^ (2 / p - 2)) t :=
    Real.continuousAt_rpow_const _ _ (Or.inl ht.1.ne')
  have c2 : ContinuousAt (fun s : ℝ => (1 - s) ^ (2 / p - 2)) t :=
    (Real.continuousAt_rpow_const _ _ (Or.inl h1.ne')).comp (continuousAt_const.sub continuousAt_id)
  exact continuousAt_const.mul (c1.add c2).sqrt

theorem Gp_deriv {p : ℝ} (hp : 1 ≤ p) {t : ℝ} (ht : t ∈ Ioo (0:ℝ) 1) :
    HasDerivAt (Gp p) (spd p t) t := by
  have hm : Measurable (spd p) := measurable_const.mul (pInt_meas p)
  exact intervalIntegral.integral_hasDerivAt_right (spd_sub hp le_rfl ht.1.le ht.2.le)
    hm.stronglyMeasurable.stronglyMeasurableAtFilter (spd_contAt hp ht)

theorem Gp_sub {p : ℝ} (hp : 1 ≤ p) {a b : ℝ} (ha : 0 ≤ a) (hab : a ≤ b) (hb : b ≤ 1) :
    Gp p b - Gp p a = ∫ t in a..b, spd p t := by
  unfold Gp
  rw [intervalIntegral.integral_interval_sub_left (spd_sub hp le_rfl (ha.trans hab) hb)
    (spd_sub hp le_rfl ha (hab.trans hb))]

theorem Gp_mono {p : ℝ} (hp : 1 ≤ p) {a b : ℝ} (ha : 0 ≤ a) (hab : a ≤ b) (hb : b ≤ 1) :
    Gp p a ≤ Gp p b := by
  have := Gp_sub hp ha hab hb
  have h2 : 0 ≤ ∫ t in a..b, spd p t :=
    intervalIntegral.integral_nonneg hab fun u _ => spd_nonneg hp u
  linarith

theorem qa_hasDeriv {p : ℝ} (hp : 1 ≤ p) {t : ℝ} (ht : t ∈ Ioo (0:ℝ) 1) :
    HasDerivAt (quarterArc p) (deriv (quarterArc p) t) t ∧ ‖deriv (quarterArc p) t‖ = spd p t := by
  obtain ⟨v, hv, hn⟩ := qa_deriv hp ht
  rw [hv.deriv]; exact ⟨hv, hn⟩

theorem up_open {p : ℝ} (hp : 1 ≤ p) {x y : ℝ} (hx : 0 < x) (hxy : x ≤ y) (hy : y < 1) :
    ‖quarterArc p y - quarterArc p x‖ ≤ Gp p y - Gp p x := by
  have hI : ∀ z ∈ Ico x y, z ∈ Ioo (0:ℝ) 1 := fun z hz => ⟨by linarith [hz.1], by linarith [hz.2]⟩
  refine image_norm_le_of_norm_deriv_right_le_deriv_boundary' (a := x) (b := y)
    (f := fun z => quarterArc p z - quarterArc p x) (f' := deriv (quarterArc p))
    ((qa_cont hp).continuousOn.sub continuousOn_const)
    (fun z hz => ((qa_hasDeriv hp (hI z hz)).1.sub_const _).hasDerivWithinAt)
    (B := fun z => Gp p z - Gp p x) (B' := spd p) (by simp)
    (((Gp_cont hp).mono (Icc_subset_Icc hx.le hy.le)).sub continuousOn_const)
    (fun z hz => ((Gp_deriv hp (hI z hz)).sub_const _).hasDerivWithinAt)
    (fun z hz => (qa_hasDeriv hp (hI z hz)).2.le) ⟨hxy, le_rfl⟩

theorem eps_tend {x : ℝ} (hx : x ∈ Icc (0:ℝ) 1) :
    Tendsto (fun e : ℝ => e + (1 - 2 * e) * x) (𝓝[>] 0) (𝓝[Icc 0 1] x) := by
  apply tendsto_nhdsWithin_iff.2 ⟨?_, ?_⟩
  · have : Tendsto (fun e : ℝ => e + (1 - 2 * e) * x) (𝓝 0) (𝓝 (0 + (1 - 2 * 0) * x)) :=
      ((continuous_id.add ((continuous_const.sub (continuous_const.mul continuous_id)).mul
        continuous_const)).tendsto 0)
    simpa using this.mono_left nhdsWithin_le_nhds
  · filter_upwards [Ioo_mem_nhdsGT (show (0:ℝ) < 1/2 by norm_num)] with e he
    constructor <;> nlinarith [he.1, he.2, hx.1, hx.2]

theorem up_closed {p : ℝ} (hp : 1 ≤ p) {x y : ℝ} (hx : 0 ≤ x) (hxy : x ≤ y) (hy : y ≤ 1) :
    ‖quarterArc p y - quarterArc p x‖ ≤ Gp p y - Gp p x := by
  have tx := eps_tend ⟨hx, hxy.trans hy⟩
  have ty := eps_tend ⟨hx.trans hxy, hy⟩
  have L : Tendsto (fun e : ℝ => ‖quarterArc p (e + (1 - 2 * e) * y) - quarterArc p (e + (1 - 2 * e) * x)‖)
      (𝓝[>] 0) (𝓝 ‖quarterArc p y - quarterArc p x‖) :=
    (((qa_cont hp).tendsto y).comp (ty.mono_right nhdsWithin_le_nhds)).sub
      (((qa_cont hp).tendsto x).comp (tx.mono_right nhdsWithin_le_nhds)) |>.norm
  have R : Tendsto (fun e : ℝ => Gp p (e + (1 - 2 * e) * y) - Gp p (e + (1 - 2 * e) * x))
      (𝓝[>] 0) (𝓝 (Gp p y - Gp p x)) :=
    ((Gp_cont hp y ⟨hx.trans hxy, hy⟩).tendsto.comp ty).sub
      ((Gp_cont hp x ⟨hx, hxy.trans hy⟩).tendsto.comp tx)
  refine le_of_tendsto_of_tendsto L R ?_
  filter_upwards [Ioo_mem_nhdsGT (show (0:ℝ) < 1/2 by norm_num)] with e he
  exact up_open hp (by nlinarith [he.1, he.2]) (by nlinarith [he.1, he.2]) (by nlinarith [he.1, he.2])

theorem var_le {E F : Type*} [PseudoEMetricSpace E] [PseudoEMetricSpace F]
    (f : ℝ → E) (h : ℝ → F) (s : Set ℝ)
    (hd : ∀ a ∈ s, ∀ b ∈ s, edist (f a) (f b) ≤ edist (h a) (h b)) :
    eVariationOn f s ≤ eVariationOn h s := by
  unfold eVariationOn
  exact iSup_mono fun p => Finset.sum_le_sum fun i _ => hd _ (p.2.2.2 _) _ (p.2.2.2 _)

theorem Gp_zero (p : ℝ) : Gp p 0 = 0 := intervalIntegral.integral_same

theorem var_up {p : ℝ} (hp : 1 ≤ p) :
    eVariationOn (quarterArc p) (Icc 0 1) ≤ ENNReal.ofReal (Gp p 1) := by
  have key : ∀ a ∈ Icc (0:ℝ) 1, ∀ b ∈ Icc (0:ℝ) 1, a ≤ b →
      edist (quarterArc p a) (quarterArc p b) ≤ edist (Gp p a) (Gp p b) := by
    intro a ha b hb hab
    rw [edist_dist, edist_dist, dist_eq_norm, norm_sub_rev, Real.dist_eq, abs_sub_comm,
      abs_of_nonneg (by linarith [Gp_mono hp ha.1 hab hb.2])]
    exact ENNReal.ofReal_le_ofReal (up_closed hp ha.1 hab hb.2)
  calc eVariationOn (quarterArc p) (Icc 0 1) ≤ eVariationOn (Gp p) (Icc 0 1) := by
        refine var_le _ _ _ fun a ha b hb => ?_
        rcases le_total a b with h | h
        · exact key a ha b hb h
        · rw [edist_comm, edist_comm (Gp p a)]; exact key b hb a ha h
    _ = ENNReal.ofReal (Gp p 1) := by
        have hm : MonotoneOn (Gp p) (Icc 0 1) := fun a ha b hb hab => Gp_mono hp ha.1 hab hb.2
        have := hm.eVariationOn_eq (a := 0) (b := 1) ⟨le_rfl, zero_le_one⟩ ⟨zero_le_one, le_rfl⟩
        simpa [Gp_zero] using this

theorem V_inc {p : ℝ} (hbv : BoundedVariationOn (quarterArc p) (Icc 0 1)) {x z : ℝ}
    (hx : 0 ≤ x) (hxz : x ≤ z) (hz : z ≤ 1) :
    ‖quarterArc p z - quarterArc p x‖ ≤
      variationOnFromTo (quarterArc p) (Icc 0 1) 0 z - variationOnFromTo (quarterArc p) (Icc 0 1) 0 x := by
  rw [variationOnFromTo.sub_right hbv.locallyBoundedVariationOn ⟨le_rfl, zero_le_one⟩
    ⟨hx.trans hxz, hz⟩ ⟨hx, hxz.trans hz⟩, variationOnFromTo.eq_of_le _ _ hxz, ← dist_eq_norm,
    dist_comm, dist_edist]
  refine ENNReal.toReal_mono (hbv.mono inter_subset_left) ?_
  exact eVariationOn.edist_le _ ⟨⟨hx, hxz.trans hz⟩, le_rfl, hxz⟩ ⟨⟨hx.trans hxz, hz⟩, hxz, le_rfl⟩

theorem V_cont {p : ℝ} (hp : 1 ≤ p) (hbv : BoundedVariationOn (quarterArc p) (Icc 0 1)) :
    ContinuousOn (variationOnFromTo (quarterArc p) (Icc 0 1) 0) (Icc 0 1) := fun x hx =>
  (hbv.continuousWithinAt_variationOnFromTo_iff ⟨le_rfl, zero_le_one⟩ hx).2
    (qa_cont hp).continuousWithinAt

theorem low_open {p : ℝ} (hp : 1 ≤ p) (hbv : BoundedVariationOn (quarterArc p) (Icc 0 1))
    {a b : ℝ} (ha : 0 < a) (hab : a ≤ b) (hb : b < 1) :
    Gp p b - Gp p a ≤
      variationOnFromTo (quarterArc p) (Icc 0 1) 0 b - variationOnFromTo (quarterArc p) (Icc 0 1) 0 a := by
  set V := variationOnFromTo (quarterArc p) (Icc 0 1) 0 with hV
  set f : ℝ → ℝ := fun z => (Gp p z - Gp p a) - (V z - V a) with hf
  have hsub : Icc a b ⊆ Icc 0 1 := Icc_subset_Icc ha.le hb.le
  have hfc : ContinuousOn f (Icc a b) :=
    (((Gp_cont hp).mono hsub).sub continuousOn_const).sub (((V_cont hp hbv).mono hsub).sub continuousOn_const)
  have key := image_le_of_liminf_slope_right_le_deriv_boundary (a := a) (b := b) hfc
    (B := fun _ => 0) (B' := fun _ => 0) (by simp [f]) continuousOn_const
    (fun _ _ => hasDerivWithinAt_const _ _ _) ?_ ⟨hab, le_rfl⟩
  · simp only [f] at key; linarith
  intro x hx r hr
  have hx01 : x ∈ Ioo (0:ℝ) 1 := ⟨by linarith [hx.1], by linarith [hx.2]⟩
  obtain ⟨hd, hn⟩ := qa_hasDeriv hp hx01
  have t1 : Tendsto (slope (quarterArc p) x) (𝓝[>] x) (𝓝 (deriv (quarterArc p) x)) :=
    (hasDerivAt_iff_tendsto_slope.1 hd).mono_left (nhdsWithin_mono _ fun z hz => ne_of_gt hz)
  have t2 : Tendsto (slope (Gp p) x) (𝓝[>] x) (𝓝 (spd p x)) :=
    (hasDerivAt_iff_tendsto_slope.1 (Gp_deriv hp hx01)).mono_left
      (nhdsWithin_mono _ fun z hz => ne_of_gt hz)
  have t3 : Tendsto (fun z => slope (Gp p) x z - ‖slope (quarterArc p) x z‖) (𝓝[>] x)
      (𝓝 (spd p x - ‖deriv (quarterArc p) x‖)) := t2.sub t1.norm
  rw [hn, sub_self] at t3
  refine Filter.Eventually.frequently ?_
  filter_upwards [t3 (Iio_mem_nhds hr), self_mem_nhdsWithin,
    nhdsWithin_le_nhds (Iio_mem_nhds hx01.2)] with z hz hzx hz1
  simp only [mem_preimage, mem_Iio] at hz
  have hzx' : x < z := hzx
  have hpos : 0 < z - x := by linarith
  have hVi := V_inc hbv hx01.1.le hzx'.le (le_of_lt hz1)
  refine lt_of_le_of_lt ?_ hz
  rw [slope_def_field, slope_def_field, slope_def_module, norm_smul, Real.norm_eq_abs,
    abs_of_pos (inv_pos.2 hpos)]
  simp only [f]
  rw [div_le_iff₀ hpos]
  have : (z - x)⁻¹ * ‖quarterArc p z - quarterArc p x‖ * (z - x) = ‖quarterArc p z - quarterArc p x‖ := by
    field_simp
  rw [sub_mul, div_mul_cancel₀ _ hpos.ne', this]
  linarith

theorem var_eq {p : ℝ} (hp : 1 ≤ p) :
    eVariationOn (quarterArc p) (Icc 0 1) = ENNReal.ofReal (Gp p 1) := by
  have hup := var_up hp
  have hne : eVariationOn (quarterArc p) (Icc 0 1) ≠ ⊤ := ne_top_of_le_ne_top ENNReal.ofReal_ne_top hup
  have hbv : BoundedVariationOn (quarterArc p) (Icc 0 1) := hne
  set W := (eVariationOn (quarterArc p) (Icc 0 1)).toReal with hW
  have hab : ∀ a b : ℝ, 0 < a → a ≤ b → b < 1 → Gp p b - Gp p a ≤ W := by
    intro a b ha hab hb
    have h1 := low_open hp hbv ha hab hb
    rw [variationOnFromTo.sub_right hbv.locallyBoundedVariationOn ⟨le_rfl, zero_le_one⟩
      ⟨by linarith, hb.le⟩ ⟨ha.le, by linarith⟩, variationOnFromTo.eq_of_le _ _ hab] at h1
    exact h1.trans (ENNReal.toReal_mono hne (eVariationOn.mono _ inter_subset_left))
  have t0 := eps_tend (x := 0) ⟨le_rfl, zero_le_one⟩
  have t1 := eps_tend (x := 1) ⟨zero_le_one, le_rfl⟩
  have R : Tendsto (fun e : ℝ => Gp p (e + (1 - 2 * e) * 1) - Gp p (e + (1 - 2 * e) * 0))
      (𝓝[>] 0) (𝓝 (Gp p 1 - Gp p 0)) :=
    ((Gp_cont hp 1 ⟨zero_le_one, le_rfl⟩).tendsto.comp t1).sub
      ((Gp_cont hp 0 ⟨le_rfl, zero_le_one⟩).tendsto.comp t0)
  have hG : Gp p 1 ≤ W := by
    have := le_of_tendsto (b := W) R ?_
    · simpa [Gp_zero] using this
    filter_upwards [Ioo_mem_nhdsGT (show (0:ℝ) < 1/2 by norm_num)] with e he
    exact hab _ _ (by nlinarith [he.1, he.2]) (by nlinarith [he.1, he.2]) (by nlinarith [he.1, he.2])
  refine le_antisymm hup ?_
  rw [← ENNReal.ofReal_toReal hne]
  exact ENNReal.ofReal_le_ofReal hG


end W04

theorem W04.proof_1 : W04.stmt_1 := by
  intro p hp
  have hr : 0 < 1 / p := by positivity
  refine ⟨?_, ?_, ?_⟩
  · intro t ht
    have a1 : 0 ≤ t ^ (1 / p) := Real.rpow_nonneg ht.1 _
    have a2 : 0 ≤ (1 - t) ^ (1 / p) := Real.rpow_nonneg (by linarith [ht.2]) _
    refine ⟨?_, a1, a2⟩
    show |t ^ (1 / p)| ^ p + |(1 - t) ^ (1 / p)| ^ p = 1
    rw [abs_of_nonneg a1, abs_of_nonneg a2, W04.pw1 hp ht.1, W04.pw1 hp (by linarith [ht.2])]
    ring
  · intro s hs t ht hst
    have := congrArg (fun v : W04.Plane => (v 0) ^ p) hst
    simp only [W04.quarterArc, W04.vec0] at this
    rwa [W04.pw1 hp hs.1, W04.pw1 hp ht.1] at this
  · intro v hv
    obtain ⟨hsq, h0, h1⟩ := hv
    have hsq' : |v 0| ^ p + |v 1| ^ p = 1 := hsq
    rw [abs_of_nonneg h0, abs_of_nonneg h1] at hsq'
    have u0 : 0 ≤ v 0 ^ p := Real.rpow_nonneg h0 _
    have u1 : 0 ≤ v 1 ^ p := Real.rpow_nonneg h1 _
    refine ⟨v 0 ^ p, ⟨u0, by linarith⟩, ?_⟩
    apply W04.plane_ext
    · simp only [W04.quarterArc, W04.vec0]; exact W04.pw2 hp h0
    · simp only [W04.quarterArc, W04.vec1]
      rw [show 1 - v 0 ^ p = v 1 ^ p by linarith, W04.pw2 hp h1]

theorem W04.proof_2 : W04.stmt_2 := fun _ hp _ ht => W04.qa_deriv hp ht

theorem W04.proof_3 : W04.stmt_3 := fun _ hp => W04.pInt_iInt hp

theorem W04.proof_4 : W04.stmt_4 := by
  intro p hp
  rw [W04.var_eq hp, W04.Gp]
  unfold W04.spd
  rw [intervalIntegral.integral_const_mul]

theorem W04.proof_5 : W04.stmt_5 := by
  intro p hp
  have hr : (1 / p) ≠ 0 := by positivity
  have h0 : (0:ℝ) ^ p⁻¹ = 0 := by simpa using Real.zero_rpow hr
  have hl : W04.loop p = W04.emb ∘ W04.loopP p := funext fun t => W04.loop_eq p t
  refine ⟨?_, ?_, ?_, ?_⟩
  · rw [hl]; exact W04.emb_cont.comp_continuousOn (W04.loopP_cont hp)
  · rw [W04.loop_eq, W04.loop_eq]; norm_num [W04.loopP, h0]
  · rw [hl]
    exact W04.emb_inj.comp_injOn
      (Set.LeftInvOn.injOn (f₁' := W04.ginv p) (fun t ht => W04.ginv_loop hp t ht))
  · ext v
    constructor
    · rintro ⟨t, ht, rfl⟩
      have := W04.norm_loop hp t ht
      rw [W04.loop_eq]
      exact this
    · intro hv
      have hv' : |v 0| ^ p + |v 1| ^ p = 1 := hv
      obtain ⟨t, ht, he⟩ := W04.exists_t hp (v 0, v 1) hv'
      refine ⟨t, ht, ?_⟩
      rw [W04.loop_eq, he]
      exact W04.plane_eta v

theorem W04.proof_6 : W04.stmt_6 := by
  intro p hp
  unfold W04.perimeter
  rw [W04.loop_split]
  rw [W04.piece 0 1 (by norm_num) (fun z => z) (fun z w => rfl) (fun t ht => W04.loop_Q1 t ht),
    W04.piece 1 2 (by norm_num) (fun z => (z.2, -z.1))
      (fun z w => by
        rw [W04.edist_emb, W04.edist_emb]
        congr 2; ring)
      (fun t ht => W04.loop_Q2 hp t ht),
    W04.piece 2 3 (by norm_num) (fun z => (-z.1, -z.2))
      (fun z w => by
        rw [W04.edist_emb, W04.edist_emb]
        congr 2; ring)
      (fun t ht => W04.loop_Q3 hp t ht),
    W04.piece 3 4 (by norm_num) (fun z => (-z.2, z.1))
      (fun z w => by
        rw [W04.edist_emb, W04.edist_emb]
        congr 2; ring)
      (fun t ht => W04.loop_Q4 hp t ht)]
  ring

theorem W04.proof_7 : W04.stmt_7 := by
  intro p hp
  rw [W04.proof_6 p hp, W04.proof_4 p hp, W04.Pint]
  rw [show (4:ENNReal) = ENNReal.ofReal 4 by norm_num, ← ENNReal.ofReal_mul (by norm_num)]
  congr 1; ring

theorem W04.proof_8 : W04.stmt_8 := by
  unfold W04.stmt_8 W04.Pint W04.pIntegrand
  norm_num

theorem W04.proof_9 : W04.stmt_9 := by
  have hd : ∀ x ∈ Ioo (0:ℝ) 1,
      HasDerivAt (fun t => Real.arcsin (2 * t - 1)) (1 / Real.sqrt (x * (1 - x))) x := by
    intro x hx
    have a1 : 2 * x - 1 ≠ -1 := by intro h; linarith [hx.1]
    have a2 : 2 * x - 1 ≠ 1 := by intro h; linarith [hx.2]
    have h := (Real.hasDerivAt_arcsin a1 a2).comp x
      (((hasDerivAt_id' x).const_mul 2).sub_const 1)
    convert h using 1
    · rfl
    have e : 1 - (2 * x - 1) ^ 2 = 2 ^ 2 * (x * (1 - x)) := by ring
    rw [e, Real.sqrt_mul' _ (mul_pos hx.1 (by linarith [hx.2])).le, Real.sqrt_sq (by norm_num)]
    have : 0 < Real.sqrt (x * (1 - x)) := Real.sqrt_pos.2 (mul_pos hx.1 (by linarith [hx.2]))
    field_simp
  have hc : ContinuousOn (fun t : ℝ => Real.arcsin (2 * t - 1)) (Icc 0 1) := by fun_prop
  have hi : IntervalIntegrable (fun t : ℝ => 1 / Real.sqrt (t * (1 - t))) volume 0 1 := by
    refine intervalIntegral.intervalIntegrable_deriv_of_nonneg (g := fun t => Real.arcsin (2 * t - 1))
      (by rw [uIcc_of_le zero_le_one]; exact hc) ?_ ?_
    · simpa using hd
    · intro x _; positivity
  refine ⟨hi, ?_⟩
  rw [intervalIntegral.integral_eq_sub_of_hasDerivAt_of_le zero_le_one hc hd hi]
  norm_num

theorem W04.proof_10 : W04.stmt_10 := by
  unfold W04.stmt_10 W04.Pint
  rw [show (4:ℝ) / 2 = 2 by norm_num]
  rw [← W04.proof_9.2]
  congr 1
  refine intervalIntegral.integral_congr_Ioo_of_le zero_le_one fun t ht => ?_
  have h1 : 0 < 1 - t := by linarith [ht.2]
  have h0 : t ≠ 0 := ht.1.ne'
  have h1' : 1 - t ≠ 0 := h1.ne'
  unfold W04.pIntegrand
  rw [show (2:ℝ) / 2 - 2 = -1 by norm_num, Real.rpow_neg_one, Real.rpow_neg_one,
    one_div, ← Real.sqrt_inv]
  congr 1
  field_simp
  ring

theorem W04.proof_11 : W04.stmt_11 := by
  have h2 : Tendsto (fun p : ℝ => 8 * (2 : ℝ) ^ (-1 / p)) atTop (nhds 8) := by
    have h0 : Tendsto (fun q : ℝ => -1 / q) atTop (nhds 0) := by
      have e : (fun q : ℝ => -1 / q) = fun q => -q⁻¹ := by funext q; ring
      rw [e]; simpa using (tendsto_inv_atTop_zero (𝕜 := ℝ)).neg
    have := (tendsto_const_nhds (x := (2 : ℝ))).rpow h0 (Or.inl two_ne_zero)
    simpa using this.const_mul 8
  refine tendsto_of_tendsto_of_tendsto_of_le_of_le' h2 tendsto_const_nhds ?_ ?_
  · filter_upwards [eventually_ge_atTop (1 : ℝ)] with p hp
    exact (W04.Pint_bounds hp).1
  · filter_upwards [eventually_ge_atTop (1 : ℝ)] with p hp
    exact (W04.Pint_bounds hp).2

