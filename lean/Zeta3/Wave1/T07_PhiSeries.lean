import Mathlib
/-!
# T07 — closed form of `Φ` (paper label `thm:flat`)

Source: `paper/ieee/main.tex`, Theorem `thm:flat` ("Expanding about the flat angle"),
lines 1501–1508. The claim formalized here is the closed form of `Φ` in line 1507:

> So, with $\Phi(c)=\sum_{n\ge1}\frac{c^{2n}}{4n^2-1}=\frac12-\frac{1-c^2}{2c}\operatorname{artanh}c$,
> the functions $M_2-\zeta(2)\Phi$ and $M_3+\zeta(3)\Phi-\zeta(2)\Psi$ are power series in $c^2$
> with rational coefficients, for some $\Psi\in\mathbb Q[[c^2]]$.

In the paper `c = \cos\frac\alpha2` with `0<\alpha<\pi` (line 1502), so `0 < c < 1`.

* `artanh c` is defined here as `(1/2) * log((1+c)/(1-c))`, the standard real inverse
  hyperbolic tangent on `(-1,1)` (Mathlib has no `Real.artanh`).
* The sum is stated with `HasSum` (summability + value), indexed by `n : ℕ` with the term
  `c^(2(n+1)) / (4(n+1)^2 - 1)` so that it runs over `n ≥ 1`; the denominator is computed
  in `ℝ`, so no natural-number truncation occurs.
* `stmt_1` states the identity for `0 < |c| < 1`. This is slightly more general than the
  paper's `0 < c < 1`; both sides are even in `c`, so it is the same claim. No hypothesis is
  strengthened.
* `stmt_2` is the paper's own parametrization: `0 < α < π`, `c = cos(α/2)`.

NOT COVERED: the rest of `thm:flat` (existence of the expansion `M_k = Σ μ_{k,n} c^{2n}`,
the weight/rationality structure of `μ_{k,n}`, the forms of `μ_{2,n}`, `μ_{3,n}`, and the
rationality of `M_2 - ζ(2)Φ`, `M_3 + ζ(3)Φ - ζ(2)Ψ`). These need the corner integrals `M_k(α)`,
which are defined via the paper's multi-section setup and are outside this target.
-/
-- BEGIN STATEMENTS
namespace T07

/-- The real inverse hyperbolic tangent on `(-1,1)`: `artanh c = (1/2) log((1+c)/(1-c))`. -/
noncomputable def artanh (c : ℝ) : ℝ := (1 / 2) * Real.log ((1 + c) / (1 - c))

/-- The `n`-th term (shifted so `n = 0` is the paper's `n = 1`) of
`Φ(c) = Σ_{n≥1} c^{2n} / (4n^2 - 1)`. -/
noncomputable def phiTerm (c : ℝ) (n : ℕ) : ℝ :=
  c ^ (2 * (n + 1)) / (4 * ((n : ℝ) + 1) ^ 2 - 1)

/-- thm:flat, line 1507: for `0 < |c| < 1`,
`Σ_{n≥1} c^{2n}/(4n^2-1) = 1/2 - ((1-c^2)/(2c)) artanh c` (the series converges). -/
def stmt_1 : Prop :=
  ∀ c : ℝ, 0 < |c| → |c| < 1 →
    HasSum (phiTerm c) (1 / 2 - (1 - c ^ 2) / (2 * c) * artanh c)

/-- thm:flat, lines 1502 and 1507, in the paper's variables: for `0 < α < π` and
`c = cos(α/2)`, `Φ(c) = 1/2 - ((1-c^2)/(2c)) artanh c`. -/
def stmt_2 : Prop :=
  ∀ α : ℝ, 0 < α → α < Real.pi →
    HasSum (phiTerm (Real.cos (α / 2)))
      (1 / 2 - (1 - Real.cos (α / 2) ^ 2) / (2 * Real.cos (α / 2)) * artanh (Real.cos (α / 2)))

end T07
-- END STATEMENTS
theorem T07.proof_1 : T07.stmt_1 := by
  intro c hc0 hc1
  have hc : c ≠ 0 := abs_pos.mp hc0
  have h := Real.hasSum_log_sub_log_of_abs_lt_one hc1
  have hA : HasSum (fun k : ℕ => c ^ (2 * k + 1) / (2 * (k : ℝ) + 1)) (T07.artanh c) := by
    have := h.mul_left (1 / 2)
    convert this using 1
    · ext k; field_simp
    · unfold T07.artanh
      rw [← Real.log_div (by linarith [abs_lt.mp hc1]) (by linarith [abs_lt.mp hc1])]
  have hB := (hasSum_nat_add_iff' 1).mpr hA
  simp only [Finset.sum_range_one, pow_one, Nat.cast_zero, mul_zero, zero_add, div_one] at hB
  have hC := (hA.mul_left c).sub (hB.mul_left c⁻¹)
  have hD := hC.mul_left (1 / 2)
  convert hD using 1
  · ext k
    unfold T07.phiTerm
    push_cast
    have h1 : (2 * (k : ℝ) + 1) ≠ 0 := by positivity
    have h2 : (2 * ((k : ℝ) + 1) + 1) ≠ 0 := by positivity
    have h3 : (4 * ((k : ℝ) + 1) ^ 2 - 1) ≠ 0 := by nlinarith [(Nat.cast_nonneg k : (0:ℝ) ≤ k)]
    have e : c ^ (2 * (k + 1)) = c * c ^ (2 * k + 1) := by ring
    have e2 : c ^ (2 * (k + 1) + 1) = c ^ 2 * c ^ (2 * k + 1) := by ring
    have e3 : (4 * ((k : ℝ) + 1) ^ 2 - 1) = (2 * k + 1) * (2 * ((k : ℝ) + 1) + 1) := by ring
    rw [e, e2, e3]
    field_simp
    ring
  · field_simp
    ring
theorem T07.proof_2 : T07.stmt_2 := by
  intro α h0 h1
  have hp : 0 < Real.cos (α / 2) :=
    Real.cos_pos_of_mem_Ioo ⟨by linarith, by linarith [Real.pi_pos]⟩
  have hl : Real.cos (α / 2) < 1 := by
    have := Real.cos_lt_cos_of_nonneg_of_le_pi (le_refl 0) (by linarith : α / 2 ≤ Real.pi)
      (by linarith : (0:ℝ) < α / 2)
    simpa using this
  exact T07.proof_1 _ (by rw [abs_of_pos hp]; exact hp) (by rw [abs_of_pos hp]; exact hl)
