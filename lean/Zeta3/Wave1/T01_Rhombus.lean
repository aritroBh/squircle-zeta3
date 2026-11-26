import Mathlib
/-!
# T01 -- the rhombus perimeter `P(1) = 4r`

Target id: T01.
Paper label: `thm:rhombus` (paper/ieee/main.tex, lines 1136-1145), only the claim `P(1)=4r`.
Supporting setup: main.tex lines 824-829 (the perimeter integral `eq:Pint` and "Case 1").

Quoted verbatim (main.tex 1137, 1139-1145):

  The same four steps work for $|x/A|^p+|y/B|^p=1$ with semi-axes $A,B>0$.
  Put $r=\sqrt{A^2+B^2}$ ...

  \begin{theorem}[Derivatives at the rhombus]\label{thm:rhombus}
  For $|x/A|^p+|y/B|^p=1$ with $r=\sqrt{A^2+B^2}$,
  \begin{gather*}
  P(1)=4r,\qquad P'(1)=0,\qquad P''(1)=\frac{4\pi^2A^2B^2}{3r^3},\\
  P'''(1)=\frac{24A^2B^2}{r^3}\Bigl(\zeta(3)-\frac{\pi^2}{3}\Bigr).
  \end{gather*}
  \end{theorem}

Quoted verbatim (main.tex 825-829, the case A = B = 1):

  parametrize the quarter arc in the first quadrant by $x=t^{1/p}$, $y=(1-t)^{1/p}$,
  $0\le t\le1$. Then $x^p+y^p=1$, and the arc-length formula gives
  P(p)=\frac4p\int_0^1\sqrt{t^{2/p-2}+(1-t)^{2/p-2}}\;dt .
  ... \textbf{Case 1:} $p=1$. The integrand is $\sqrt2$, so $P(1)=4\sqrt2$:
  four sides of length $\sqrt2$.

Formalization (points of the plane are `EuclideanSpace ℝ (Fin 2)`, so `dist` is the Euclidean
distance; `x 0`, `x 1` are the coordinates). Throughout `A B > 0`, as in the paper.

* stmt_1: the curve `|x/A| + |y/B| = 1` is exactly the union of the four closed segments
  joining consecutive vertices `(A,0), (0,B), (-A,0), (0,-B)`.
* stmt_2: each of the four sides has Euclidean length `r = sqrt(A^2+B^2)`.
* stmt_3: the sum of the four side lengths is `4r`.
* stmt_4: the paper's quarter-arc parametrization at `p = 1`, `t ↦ (A t, B (1-t))`, `t ∈ [0,1]`,
  lies on the curve, and its arc length (total variation, `eVariationOn`) is `r`; so the
  perimeter, four congruent quarter arcs, is `4r`.
* stmt_5: the perimeter integral. The paper writes `eq:Pint` only for `A = B = 1`; for general
  semi-axes the same parametrization `x = A t^{1/p}`, `y = B (1-t)^{1/p}` gives
  `P(p) = (4/p) ∫_0^1 sqrt(A^2 t^{2/p-2} + B^2 (1-t)^{2/p-2}) dt` (our extension, stated by the
  paper only as "the same four steps work"). We state that this formula gives `4r` at `p = 1`.
  `Real.rpow` has `0^0 = 1`, which only affects the measure-zero endpoints.

NOT COVERED: `P'(1) = 0`, `P''(1)`, `P'''(1)` of `thm:rhombus` (excluded from this target by
instruction). No general Mathlib notion of "perimeter of a convex set" is used; the perimeter is
expressed as the sum of the side lengths (stmt_3) and as four times the arc length of the
paper's quarter-arc parametrization (stmt_4, stmt_5).
-/
-- BEGIN STATEMENTS
namespace T01

open scoped ENNReal

/-- The Euclidean plane. -/
abbrev Plane := EuclideanSpace ℝ (Fin 2)

/-- `r = sqrt(A^2 + B^2)`. -/
noncomputable def r (A B : ℝ) : ℝ := Real.sqrt (A ^ 2 + B ^ 2)

/-- The curve `|x/A| + |y/B| = 1` (the case `p = 1` of `|x/A|^p+|y/B|^p=1`). -/
def curve (A B : ℝ) : Set Plane := {v | |v 0 / A| + |v 1 / B| = 1}

/-- The four vertices `(A,0), (0,B), (-A,0), (0,-B)`. -/
noncomputable def V1 (A : ℝ) : Plane := !₂[A, 0]
noncomputable def V2 (B : ℝ) : Plane := !₂[0, B]
noncomputable def V3 (A : ℝ) : Plane := !₂[-A, 0]
noncomputable def V4 (B : ℝ) : Plane := !₂[0, -B]

/-- The quarter arc at `p = 1`: `x = A t^{1/p}`, `y = B (1-t)^{1/p}` with `p = 1`. -/
noncomputable def quarterArc (A B : ℝ) (t : ℝ) : Plane := !₂[A * t, B * (1 - t)]

/-- The perimeter integral `P(p) = (4/p) ∫_0^1 sqrt(A^2 t^{2/p-2} + B^2 (1-t)^{2/p-2}) dt`
(paper's `eq:Pint`, with semi-axes `A, B`). -/
noncomputable def P (A B p : ℝ) : ℝ :=
  4 / p * ∫ t in (0 : ℝ)..1,
    Real.sqrt (A ^ 2 * t ^ (2 / p - 2) + B ^ 2 * (1 - t) ^ (2 / p - 2))

/-- "the curve |x/A|+|y/B|=1 is the rhombus with vertices (±A,0),(0,±B)": it is the union of
the four sides. -/
def stmt_1 : Prop :=
  ∀ A B : ℝ, 0 < A → 0 < B →
    curve A B = segment ℝ (V1 A) (V2 B) ∪ segment ℝ (V2 B) (V3 A) ∪
      segment ℝ (V3 A) (V4 B) ∪ segment ℝ (V4 B) (V1 A)

/-- Each side of the rhombus has length `r` ("four sides of length ..."). -/
def stmt_2 : Prop :=
  ∀ A B : ℝ, 0 < A → 0 < B →
    dist (V1 A) (V2 B) = r A B ∧ dist (V2 B) (V3 A) = r A B ∧
    dist (V3 A) (V4 B) = r A B ∧ dist (V4 B) (V1 A) = r A B

/-- `P(1) = 4r`, as the sum of the four side lengths. -/
def stmt_3 : Prop :=
  ∀ A B : ℝ, 0 < A → 0 < B →
    dist (V1 A) (V2 B) + dist (V2 B) (V3 A) + dist (V3 A) (V4 B) + dist (V4 B) (V1 A)
      = 4 * r A B

/-- The quarter-arc parametrization at `p = 1` lies on the curve and has arc length `r`,
so `P(1) = 4 r`. -/
def stmt_4 : Prop :=
  ∀ A B : ℝ, 0 < A → 0 < B →
    (∀ t ∈ Set.Icc (0 : ℝ) 1, quarterArc A B t ∈ curve A B) ∧
    eVariationOn (quarterArc A B) (Set.Icc 0 1) = ENNReal.ofReal (r A B)

/-- `P(1) = 4r` from the perimeter integral. -/
def stmt_5 : Prop :=
  ∀ A B : ℝ, 0 < A → 0 < B → P A B 1 = 4 * r A B

end T01
-- END STATEMENTS
namespace T01

lemma dist2 (a b c d : ℝ) : dist (!₂[a, b] : Plane) !₂[c, d] = Real.sqrt ((a - c) ^ 2 + (b - d) ^ 2) := by
  rw [EuclideanSpace.dist_eq]; simp [Fin.sum_univ_two, Real.dist_eq, sq_abs]

lemma pl_ext (v w : Plane) (h0 : v 0 = w 0) (h1 : v 1 = w 1) : v = w := by
  ext i; fin_cases i <;> simpa

end T01

theorem T01.proof_1 : T01.stmt_1 := by
  intro A B hA hB
  ext v
  simp only [T01.curve, Set.mem_setOf_eq, Set.mem_union, segment, T01.V1, T01.V2, T01.V3, T01.V4]
  constructor
  · intro h
    rcases le_or_gt 0 (v 0) with h0 | h0 <;> rcases le_or_gt 0 (v 1) with h1 | h1
    · rw [abs_of_nonneg (div_nonneg h0 hA.le), abs_of_nonneg (div_nonneg h1 hB.le)] at h
      refine Or.inl (Or.inl (Or.inl ⟨v 0 / A, v 1 / B, by positivity, by positivity, h, ?_⟩))
      apply T01.pl_ext <;> simp <;> field_simp
    · rw [abs_of_nonneg (div_nonneg h0 hA.le), abs_of_neg (div_neg_of_neg_of_pos h1 hB)] at h
      refine Or.inr ⟨-(v 1) / B, v 0 / A, by apply div_nonneg <;> linarith, by positivity, by rw [neg_div] at *; linarith, ?_⟩
      apply T01.pl_ext <;> simp <;> field_simp
    · rw [abs_of_neg (div_neg_of_neg_of_pos h0 hA), abs_of_nonneg (div_nonneg h1 hB.le)] at h
      refine Or.inl (Or.inl (Or.inr ⟨v 1 / B, -(v 0) / A, by positivity, by apply div_nonneg <;> linarith, by rw [neg_div] at *; linarith, ?_⟩))
      apply T01.pl_ext <;> simp <;> field_simp
    · rw [abs_of_neg (div_neg_of_neg_of_pos h0 hA), abs_of_neg (div_neg_of_neg_of_pos h1 hB)] at h
      refine Or.inl (Or.inr ⟨-(v 0) / A, -(v 1) / B, by apply div_nonneg <;> linarith, by apply div_nonneg <;> linarith, by rw [neg_div, neg_div] at *; linarith, ?_⟩)
      apply T01.pl_ext <;> simp <;> field_simp
  · rintro (((⟨a, b, ha, hb, hab, rfl⟩ | ⟨a, b, ha, hb, hab, rfl⟩) | ⟨a, b, ha, hb, hab, rfl⟩) | ⟨a, b, ha, hb, hab, rfl⟩) <;>
    · simp [abs_div, abs_mul, abs_of_pos hA, abs_of_pos hB, abs_of_nonneg ha, abs_of_nonneg hb, hA.ne', hB.ne']
      linarith
theorem T01.proof_2 : T01.stmt_2 := by
  intro A B hA hB
  simp only [T01.V1, T01.V2, T01.V3, T01.V4, T01.dist2, T01.r]
  refine ⟨?_, ?_, ?_, ?_⟩ <;> congr 1 <;> ring
theorem T01.proof_3 : T01.stmt_3 := by
  intro A B hA hB
  obtain ⟨h1, h2, h3, h4⟩ := T01.proof_2 A B hA hB
  rw [h1, h2, h3, h4]; ring
theorem T01.proof_4 : T01.stmt_4 := by
  intro A B hA hB
  refine ⟨?_, ?_⟩
  · intro t ht
    simp only [T01.curve, T01.quarterArc, Set.mem_setOf_eq]
    simp
    rw [mul_div_cancel_left₀ _ hA.ne', mul_div_cancel_left₀ _ hB.ne', abs_of_nonneg ht.1,
      abs_of_nonneg (by linarith [ht.2])]
    ring
  · apply le_antisymm
    · have hL : LipschitzOnWith (Real.toNNReal (T01.r A B)) (T01.quarterArc A B) (Set.Icc 0 1) := by
        apply LipschitzOnWith.of_dist_le_mul
        intro x _ y _
        simp only [T01.quarterArc, T01.dist2, Real.coe_toNNReal', T01.r]
        rw [max_eq_left (Real.sqrt_nonneg _), Real.dist_eq, ← Real.sqrt_sq_eq_abs,
          ← Real.sqrt_mul (by positivity)]
        apply le_of_eq; congr 1; ring
      have h1 := hL.comp_eVariationOn_le (g := id) (s := Set.Icc (0:ℝ) 1) (Set.mapsTo_id _)
      have h2 : eVariationOn (id : ℝ → ℝ) (Set.Icc 0 1) ≤ ENNReal.ofReal 1 := by
        have := (monotone_id.monotoneOn (Set.Icc (0:ℝ) 1)).eVariationOn_le (a := 0) (b := 1)
          (by simp) (by simp)
        simpa using this
      calc eVariationOn (T01.quarterArc A B) (Set.Icc 0 1)
          ≤ (Real.toNNReal (T01.r A B) : ENNReal) * eVariationOn (id : ℝ → ℝ) (Set.Icc 0 1) := h1
        _ ≤ (Real.toNNReal (T01.r A B) : ENNReal) * ENNReal.ofReal 1 := by gcongr
        _ = ENNReal.ofReal (T01.r A B) := by simp [ENNReal.ofReal]
    · have := eVariationOn.edist_le (T01.quarterArc A B) (s := Set.Icc (0:ℝ) 1) (x := 1) (y := 0)
        (by simp) (by simp)
      refine le_trans (le_of_eq ?_) this
      rw [edist_dist, T01.quarterArc, T01.quarterArc, T01.dist2]
      congr 1
      simp [T01.r]
theorem T01.proof_5 : T01.stmt_5 := by
  intro A B hA hB
  simp [T01.P, T01.r]
