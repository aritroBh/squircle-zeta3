import Mathlib
/-!
# T11 — printed decimals of the constants of `thm:xnyn`

Source: `paper/ieee/main.tex`, Theorem `thm:xnyn` (lines 1880-1891) and the
"Numerically" display right after it (lines 1893-1898).

**Definitions, `thm:xnyn`, lines 1881 and 1887-1888:**

> Let $\ell=\log n$ and $\lambda=\log2$.

> \beta_1=\frac{\lambda}{2}-1,\qquad \beta_0=-(1+\sqrt2)+\frac{\lambda}{2}-\frac{\pi^2}{12}-\frac{\lambda^2}{8},\\
> b_3=\frac{1-\lambda}{2},\qquad c_3=1-\frac{\lambda}{2}+\frac{\lambda^2}{4}+\frac{\pi^2}{6},

**Numerical claim, lines 1893-1897:**

> Numerically
> \beta_1=-0.6534264097\ldots,\qquad\beta_0=-2.9501636322\ldots,\\
> d_3=2.4722752293\ldots

`β₁`, `β₀` below are the paper's closed forms with `λ = Real.log 2`, `√2 = Real.sqrt 2`,
`π = Real.pi` (all arguments positive, so no Lean convention is involved).
A printed value `-0.6534264097…` (truncated decimal expansion of a negative number) is read as
`-0.6534264098 < β₁ < -0.6534264097` (strict on both sides; the true value
`-0.65342640972…` is not a 10-digit decimal). Likewise for `β₀`.

* `stmt_1`: confirms all 10 printed decimals of `β₁` (margin about 2·10⁻¹¹ on the upper side,
  so a proof needs `log 2` to within about 4·10⁻¹¹; Mathlib's `Real.log_two_gt_d9`/`lt_d9`
  alone give only 9 decimals).
* `stmt_2`: confirms all 10 printed decimals of `β₀` (margins about 4·10⁻¹¹ and 6·10⁻¹¹).

NOT COVERED:
* `b_3` and `c_3`: the paper prints no decimals for them (only for `β₁`, `β₀`, `d_3`), so
  there are no printed digits to confirm. (Their values are `b_3 = 0.1534264097…`,
  `c_3 = 2.4184737300…`, not claimed in the paper.)
* `d_3 = 2.4722752293…`: needs `ζ(3)` to about 10⁻¹⁰; no honest bound on `riemannZeta 3`
  that sharp is available in Mathlib, so it is left out.
* The sentence "at $n=10^5$ the theorem gives $R(n)+R_3(n)/n=-76.7411$" (needs `d_3`) and the
  asymptotic expansion of `A(n)` itself (a separate target).
-/
-- BEGIN STATEMENTS
namespace T11

/-- `λ = log 2` (`thm:xnyn`, line 1881). -/
noncomputable def lam : ℝ := Real.log 2

/-- `β₁ = λ/2 - 1` (`thm:xnyn`, line 1887). -/
noncomputable def beta1 : ℝ := lam / 2 - 1

/-- `β₀ = -(1+√2) + λ/2 - π²/12 - λ²/8` (`thm:xnyn`, line 1887). -/
noncomputable def beta0 : ℝ :=
  -(1 + Real.sqrt 2) + lam / 2 - Real.pi ^ 2 / 12 - lam ^ 2 / 8

/-- "Numerically $\beta_1=-0.6534264097\ldots$" (line 1895): all 10 printed decimals. -/
def stmt_1 : Prop := (-0.6534264098 : ℝ) < beta1 ∧ beta1 < (-0.6534264097 : ℝ)

/-- "Numerically ... $\beta_0=-2.9501636322\ldots$" (line 1895): all 10 printed decimals. -/
def stmt_2 : Prop := (-2.9501636323 : ℝ) < beta0 ∧ beta0 < (-2.9501636322 : ℝ)

end T11
-- END STATEMENTS
theorem T11.lam_bounds : (0.69314718055 : ℝ) < T11.lam ∧ T11.lam < (0.69314718057 : ℝ) := by
  have h := Real.abs_log_sub_add_sum_range_le (x := (1/2 : ℝ)) (by norm_num [abs_of_pos]) 40
  have hl : Real.log (1 - 1/2 : ℝ) = - Real.log 2 := by
    rw [show (1 - 1/2 : ℝ) = (2:ℝ)⁻¹ by norm_num, Real.log_inv]
  rw [hl] at h
  simp only [Finset.sum_range_succ, Finset.sum_range_zero] at h
  norm_num [abs_le] at h
  unfold T11.lam
  constructor <;> linarith [h.1, h.2]

theorem T11.sqrt2_bounds : (1.41421356237 : ℝ) < Real.sqrt 2 ∧ Real.sqrt 2 < (1.41421356238 : ℝ) := by
  constructor
  · apply Real.lt_sqrt_of_sq_lt; norm_num
  · rw [Real.sqrt_lt' (by norm_num)]; norm_num

theorem T11.proof_1 : T11.stmt_1 := by
  obtain ⟨h1, h2⟩ := T11.lam_bounds
  unfold T11.stmt_1 T11.beta1
  constructor <;> linarith

theorem T11.proof_2 : T11.stmt_2 := by
  obtain ⟨h1, h2⟩ := T11.lam_bounds
  obtain ⟨s1, s2⟩ := T11.sqrt2_bounds
  have p1 := Real.pi_gt_d20
  have p2 := Real.pi_lt_d20
  unfold T11.stmt_2 T11.beta0
  have q1 : (3.14159265358979323846 : ℝ)^2 < Real.pi^2 := by gcongr
  have q2 : Real.pi^2 < (3.14159265358979323847 : ℝ)^2 := by gcongr
  have l1 : (0.69314718055 : ℝ)^2 < T11.lam^2 := by gcongr
  have l2 : T11.lam^2 < (0.69314718057 : ℝ)^2 := by gcongr
  constructor <;> nlinarith
