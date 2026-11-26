import Mathlib
/-!
# T03 — the coefficients `ω_n` and `V_n` of Proposition `prop:an`

Target id: T03.
Paper label: `prop:an` (only the `ω_n`, `V_n` part), paper/ieee/main.tex lines 1163–1178.

Claim quoted verbatim (main.tex lines 1169–1173):

  For odd $n$ its coefficient in $a_n$ is
  \[
  \omega_n=\frac{4(-1)^n}{n}\bigl(1-2^{1-n}V_n\bigr),\qquad
  V_n=\frac{d^n}{dt^n}\sqrt{\tfrac{1+e^{4t}}{2}}\,\Big|_{t=0},
  \]
  so $\omega_3=1$, $\omega_5=-\frac34$, and $\omega_n\ne0$ for every odd $n\ge3$.

and from the following paragraph (line 1178): "$V_3,V_5,V_7=7,1,967$ give
$\omega_3,\omega_5,\omega_7=1,-\frac34,\frac{129}{16}$".

Formalization choices:
* `V n` is `iteratedDeriv n f 0` with `f t = Real.sqrt ((1 + Real.exp (4 t)) / 2)`.
  The radicand is `> 0` for all real `t`, so `Real.sqrt`'s convention at negatives
  plays no role and `f` is smooth; `iteratedDeriv` is the honest `n`-th derivative.
* `ω n` uses the integer power `(2:ℝ) ^ ((1:ℤ) - n)` so that `2^{1-n}` is the real
  number `2^{1-n}` (no truncated ℕ subtraction), and division by `(n:ℝ)`, which is
  nonzero for the cases `n = 3, 5` stated here.

NOT COVERED:
* `ω_n ≠ 0` for every odd `n ≥ 3` (excluded by the task brief), and the parity
  argument that `V_n` is odd for every `n`.
* The identification of `ω_n` as the coefficient of `ζ(n)` in the Taylor coefficient
  `a_n = P^{(n)}(1)/(√2 n!)` of the diamond-end perimeter (requires the full
  perimeter expansion; not part of this target).
* `V_7 = 967`, `ω_7 = 129/16` (stated in the follow-up paragraph, not in the
  proposition itself).
-/
-- BEGIN STATEMENTS
namespace T03

/-- The function `t ↦ √((1 + e^{4t})/2)` whose derivatives at `0` define `V_n`. -/
noncomputable def f (t : ℝ) : ℝ := Real.sqrt ((1 + Real.exp (4 * t)) / 2)

/-- `V_n = (d^n/dt^n) √((1+e^{4t})/2) |_{t=0}`. -/
noncomputable def V (n : ℕ) : ℝ := iteratedDeriv n f 0

/-- `ω_n = (4(-1)^n / n) (1 - 2^{1-n} V_n)`. -/
noncomputable def ω (n : ℕ) : ℝ :=
  4 * (-1 : ℝ) ^ n / (n : ℝ) * (1 - (2 : ℝ) ^ ((1 : ℤ) - (n : ℤ)) * V n)

/-- Line 1178: "$V_3,V_5,\dots=7,1,\dots$" — the third derivative at `0` is `7`. -/
def stmt_1 : Prop := V 3 = 7

/-- Line 1178: the fifth derivative at `0` is `1`. -/
def stmt_2 : Prop := V 5 = 1

/-- Line 1173: "so $\omega_3=1$". -/
def stmt_3 : Prop := ω 3 = 1

/-- Line 1173: "$\omega_5=-\frac34$". -/
def stmt_4 : Prop := ω 5 = -3 / 4

end T03
-- END STATEMENTS
namespace T03

lemma g_pos (t : ℝ) : 0 < (1 + Real.exp (4 * t)) / 2 := by positivity

lemma f_contDiff : ContDiff ℝ (⊤ : ℕ∞) f := by
  unfold f
  apply ContDiff.sqrt
  · fun_prop
  · intro t; exact (g_pos t).ne'

lemma f_sq (t : ℝ) : f t ^ 2 = (1 + Real.exp (4 * t)) / 2 := by
  unfold f; rw [Real.sq_sqrt (g_pos t).le]

lemma f_pos (t : ℝ) : 0 < f t := Real.sqrt_pos.mpr (g_pos t)

lemma f_deriv (t : ℝ) : f t * deriv f t = Real.exp (4 * t) := by
  have hd : HasDerivAt f (deriv f t) t :=
    ((f_contDiff.differentiable (by simp)) t).hasDerivAt
  have h2 : HasDerivAt (fun s => f s ^ 2) (2 * f t * deriv f t) t := by
    have := hd.mul hd
    have e : (fun s => f s ^ 2) = fun s => f s * f s := funext fun s => sq (f s)
    rw [e]; convert this using 1; ring
  have h3 : HasDerivAt (fun s => (1 + Real.exp (4 * s)) / 2) (Real.exp (4 * t) * 4 / 2) t := by
    have := ((hasDerivAt_id t).const_mul 4).exp
    simpa using (this.const_add 1).div_const 2
  have h4 : (fun s => f s ^ 2) = fun s => (1 + Real.exp (4 * s)) / 2 := funext f_sq
  rw [h4] at h2
  have := h2.unique h3
  linarith

lemma V_rec (n : ℕ) :
    ∑ i ∈ Finset.range (n + 1), (n.choose i : ℝ) * V i * V (n - i + 1) = 4 ^ n := by
  have hf : ContDiff ℝ n f := f_contDiff.of_le (by exact_mod_cast le_top)
  have hdf : ContDiff ℝ n (deriv f) := by
    exact (contDiff_infty_iff_deriv.mp (by simpa using f_contDiff)).2.of_le
      (by exact_mod_cast le_top)
  have h : f * deriv f = fun t => Real.exp (4 * t) := funext f_deriv
  have := iteratedDeriv_mul (n := n) (x := 0) hf.contDiffAt hdf.contDiffAt
  rw [h, iteratedDeriv_exp_const_mul] at this
  simp at this
  rw [this]
  apply Finset.sum_congr rfl
  intro i hi
  have hi' : i ≤ n := Nat.lt_succ_iff.mp (Finset.mem_range.mp hi)
  simp only [V, iteratedDeriv_succ']

lemma V0 : V 0 = 1 := by
  simp [V, f]

lemma V1 : V 1 = 1 := by
  have h := V_rec 0
  simpa [V0] using h

lemma V2 : V 2 = 3 := by
  have h := V_rec 1
  simp [Finset.sum_range_succ, V0, V1] at h
  linarith

lemma V3 : V 3 = 7 := by
  have h := V_rec 2
  simp [Finset.sum_range_succ, Nat.choose, V0, V1, V2] at h
  linarith

lemma V4 : V 4 = 9 := by
  have h := V_rec 3
  simp [Finset.sum_range_succ, Nat.choose, V0, V1, V2, V3] at h
  linarith

lemma V5 : V 5 = 1 := by
  have h := V_rec 4
  simp [Finset.sum_range_succ, Nat.choose, V0, V1, V2, V3, V4] at h
  linarith

end T03

theorem T03.proof_1 : T03.stmt_1 := T03.V3
theorem T03.proof_2 : T03.stmt_2 := T03.V5
theorem T03.proof_3 : T03.stmt_3 := by
  unfold T03.stmt_3 T03.ω
  rw [T03.V3]; norm_num
theorem T03.proof_4 : T03.stmt_4 := by
  unfold T03.stmt_4 T03.ω
  rw [T03.V5]; norm_num
