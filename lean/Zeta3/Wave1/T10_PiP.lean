import Mathlib
/-!
# T10 — the term `4·2^{-1/p}` of `π_p` (paper label `thm:pip`)

Source: `paper/ieee/main.tex`, Theorem `thm:pip` ("Expansion of $\pi_p$"), lines 1336–1350,
with the splitting of `paper/ieee/proofs/sec_pip.tex`, Lemma `lem:pip:split` (lines 52–66)
and line 132. The theorem reads (lines 1337–1344):

> $\pi_p=4\cdot2^{-1/p}+E(p)$ with an explicit integral $E(p)$ (given in the proof), and
> \[ \pi_p=4\Bigl(1+\frac{N_0}{p}+\frac{N_1}{p^2}+\frac{N_2}{p^3}+\frac{N_3}{p^4}+O(p^{-5})\Bigr), \]
> where $N_0=-\log2$, $N_1=\frac12\log^22$, $N_2=\frac{\pi^2}{24}-\frac{\log^32}{6}$ and
> \[ N_3=\frac{\log^42}{24}+\frac{\pi^2}{24}(1-\log2)+\frac{\zeta(3)}{4} . \]

and the proof (sec_pip.tex line 65 and line 132):

> The first part gives $\frac4p\int_0^{1/2}t^{-1+h}dt=\frac4p\,p\,2^{-h}=4\cdot2^{-1/p}$.
> ... Finally $2^{-h}=\sum_k\frac{(-\log2)^{k}}{k!}h^k$ and $\pi_p/4=2^{-h}+\hat E(h)$.

What is covered here is exactly the contribution of the term `4·2^{-1/p}`, with `h = ε = 1/p`:

* `stmt_1`: the main part of the splitting, `(4/p)∫_0^{1/2} t^{-1+1/p} dt = 4·2^{-1/p}` for `p > 1`
  (sec_pip.tex line 65).
* `stmt_2`: the Taylor coefficients of `ε ↦ 2^{-ε}` at `0` are `(-log 2)^k` (so the
  coefficient of `ε^k` is `(-log 2)^k/k!`, sec_pip.tex line 132).
* `stmt_3`: `2^{-ε} = Σ_{k≤4} (-log 2)^k ε^k/k! + O(ε^5)` on `0 ≤ ε ≤ 1`, with the explicit
  constant `(log 2)^5/120` (the Lagrange bound). The paper prints no constant for this
  `O(ε^5)`; the constant is our choice (it is the sharp Lagrange constant, valid for all ε ≥ 0).
* `stmt_4`: the printed values match the Taylor coefficients: `N_0 = -log 2 = (-log 2)^1/1!`,
  `N_1 = ½ log²2 = (-log 2)^2/2!`, and the pure-`log 2` parts of `N_2`, `N_3`,
  `-log³2/6 = (-log 2)^3/3!`, `log⁴2/24 = (-log 2)^4/4!`.
* `stmt_5`: in the variable `p` (`p ≥ 1`), with the paper's printed `N_0`, `N_1` and the
  pure-log parts of `N_2`, `N_3`:
  `|4·2^{-1/p} - 4(1 + N_0/p + N_1/p² + N_2^{log}/p³ + N_3^{log}/p⁴)| ≤ 4 (log 2)^5/(120 p^5)`.

`2^{-ε}` is `Real.rpow 2 (-ε)`; all divisions are by nonzero reals under the hypotheses.

NOT COVERED:
* The integral `E(p)` and its expansion `E(p) = (π²/6)p^{-3} + (π²/6 + ζ(3) - π² log2/6)p^{-4}
  + O(p^{-5})`, i.e. the `π²/24` part of `N_2` and the `π²/24 (1 - log 2) + ζ(3)/4` part of
  `N_3`; hence also the full expansion of `π_p` itself (the definition of `π_p` as an integral
  and the identity `π_p = 4·2^{-1/p} + E(p)` are not stated here; only its main part is).
* `N_4`, `N_5`, the expansion to every order, its divergence, and the growth law of `N_k`.
-/
-- BEGIN STATEMENTS
namespace T10

/-- The function `ε ↦ 2^{-ε}` (real power). -/
noncomputable def twoPow (ε : ℝ) : ℝ := (2 : ℝ) ^ (-ε)

/-- The degree-4 Taylor polynomial `Σ_{k≤4} (-log 2)^k ε^k / k!`. -/
noncomputable def taylor4 (ε : ℝ) : ℝ :=
  ∑ k ∈ Finset.range 5, (-Real.log 2) ^ k * ε ^ k / (Nat.factorial k : ℝ)

/-- Printed `N_0 = -log 2` (main.tex line 1341). -/
noncomputable def N0 : ℝ := -Real.log 2
/-- Printed `N_1 = ½ log² 2` (main.tex line 1341). -/
noncomputable def N1 : ℝ := (1 / 2) * Real.log 2 ^ 2
/-- Pure-`log 2` part of the printed `N_2 = π²/24 - log³2/6` (main.tex line 1341). -/
noncomputable def N2log : ℝ := -(Real.log 2 ^ 3 / 6)
/-- Pure-`log 2` part of the printed `N_3 = log⁴2/24 + π²/24(1-log 2) + ζ(3)/4`
(main.tex line 1343). -/
noncomputable def N3log : ℝ := Real.log 2 ^ 4 / 24

/-- sec_pip.tex line 65: for `p > 1`, `(4/p) ∫_0^{1/2} t^{-1+1/p} dt = 4·2^{-1/p}`. -/
def stmt_1 : Prop :=
  ∀ p : ℝ, 1 < p →
    (4 / p) * ∫ t in (0 : ℝ)..(1 / 2), t ^ (-1 + 1 / p) = 4 * (2 : ℝ) ^ (-(1 / p))

/-- sec_pip.tex line 132: the `k`-th derivative of `ε ↦ 2^{-ε}` at `0` is `(-log 2)^k`,
so `2^{-ε} = Σ_k (-log 2)^k ε^k / k!` term by term. -/
def stmt_2 : Prop :=
  ∀ k : ℕ, iteratedDeriv k twoPow 0 = (-Real.log 2) ^ k

/-- `2^{-ε} = Σ_{k≤4} (-log 2)^k ε^k/k! + O(ε^5)` on `0 ≤ ε ≤ 1`, with the explicit constant
`(log 2)^5/120` (constant not printed in the paper). -/
def stmt_3 : Prop :=
  ∀ ε : ℝ, 0 ≤ ε → ε ≤ 1 → |twoPow ε - taylor4 ε| ≤ Real.log 2 ^ 5 / 120 * ε ^ 5

/-- main.tex line 1341–1343: `N_0`, `N_1` and the pure-log parts of `N_2`, `N_3` are the
Taylor coefficients `(-log 2)^{k+1}/(k+1)!` of `2^{-h}`. -/
def stmt_4 : Prop :=
  N0 = (-Real.log 2) ^ 1 / (Nat.factorial 1 : ℝ) ∧
  N1 = (-Real.log 2) ^ 2 / (Nat.factorial 2 : ℝ) ∧
  N2log = (-Real.log 2) ^ 3 / (Nat.factorial 3 : ℝ) ∧
  N3log = (-Real.log 2) ^ 4 / (Nat.factorial 4 : ℝ)

/-- The `4·2^{-1/p}` part of thm:pip: for `p ≥ 1`,
`4·2^{-1/p} = 4(1 + N_0/p + N_1/p² + N_2^{log}/p³ + N_3^{log}/p⁴) + O(p^{-5})`,
with the explicit bound `4 (log 2)^5 / (120 p^5)`. -/
def stmt_5 : Prop :=
  ∀ p : ℝ, 1 ≤ p →
    |4 * (2 : ℝ) ^ (-(1 / p)) - 4 * (1 + N0 / p + N1 / p ^ 2 + N2log / p ^ 3 + N3log / p ^ 4)|
      ≤ 4 * Real.log 2 ^ 5 / (120 * p ^ 5)

end T10
-- END STATEMENTS
namespace T10

lemma twoPow_eq : twoPow = fun x => Real.exp ((-Real.log 2) * x) := by
  funext x
  unfold twoPow
  rw [Real.rpow_def_of_pos (by norm_num : (0:ℝ) < 2)]
  congr 1; ring

lemma iter_twoPow (k : ℕ) : iteratedDeriv k twoPow = fun x => (-Real.log 2) ^ k * Real.exp ((-Real.log 2) * x) := by
  rw [twoPow_eq]; exact iteratedDeriv_exp_const_mul k _

lemma twoPow_contDiff (n : WithTop ℕ∞) : ContDiff ℝ n twoPow := by
  rw [twoPow_eq]; fun_prop

lemma taylor4_eq (ε : ℝ) (hε : 0 < ε) :
    taylorWithinEval twoPow 4 (Set.Icc 0 ε) 0 ε = taylor4 ε := by
  rw [taylor_within_apply]
  unfold taylor4
  apply Finset.sum_congr rfl
  intro k _
  rw [iteratedDerivWithin_eq_iteratedDeriv (uniqueDiffOn_Icc hε) (twoPow_contDiff _).contDiffAt
    (Set.left_mem_Icc.2 hε.le), iter_twoPow]
  simp only [smul_eq_mul, mul_zero, Real.exp_zero, mul_one, sub_zero]
  field_simp

end T10

theorem T10.proof_1 : T10.stmt_1 := by
  intro p hp
  have hp0 : p ≠ 0 := by positivity
  have hr : (-1 + 1 / p) + 1 = 1 / p := by ring
  have hpos : 0 < 1 / p := by positivity
  rw [integral_rpow (Or.inl (by linarith))]
  rw [hr, Real.zero_rpow (by positivity)]
  have : ((1:ℝ) / 2) ^ (1 / p) = (2:ℝ) ^ (-(1 / p)) := by
    rw [Real.rpow_neg (by norm_num), one_div, Real.inv_rpow (by norm_num)]
  rw [this]
  field_simp
  simp
theorem T10.proof_2 : T10.stmt_2 := by
  intro k
  rw [T10.iter_twoPow]; simp
theorem T10.proof_3 : T10.stmt_3 := by
  intro ε h0 h1
  rcases h0.eq_or_lt with h | hε
  · subst h
    simp [T10.taylor4, T10.twoPow, Finset.sum_range_succ]
  obtain ⟨x', hx', h⟩ := taylor_mean_remainder_lagrange_iteratedDeriv (f := T10.twoPow) (x := ε) (x₀ := 0) (n := 4)
    hε.ne ((T10.twoPow_contDiff _).contDiffOn)
  rw [Set.uIcc_of_le hε.le, T10.taylor4_eq ε hε] at h
  rw [T10.iter_twoPow] at h
  have hx : x' ∈ Set.Ioo 0 ε := by
    simpa [Set.uIoo, Set.uIcc_of_le hε.le, hε.le] using hx'
  have hL : 0 < Real.log 2 := Real.log_pos (by norm_num)
  have he : Real.exp ((-Real.log 2) * x') ≤ 1 := by
    rw [Real.exp_le_one_iff]; nlinarith [hx.1]
  have he0 := Real.exp_pos ((-Real.log 2) * x')
  rw [h]
  simp only [sub_zero]
  have : ((4:ℕ) + 1).factorial = 120 := by decide
  rw [abs_div, abs_mul, abs_mul, abs_pow, abs_pow, abs_neg, abs_of_pos hL, abs_of_pos hε, abs_of_pos he0]
  norm_num [Nat.factorial]
  rw [div_le_iff₀ (by norm_num)]
  have h5 : 0 ≤ Real.log 2 ^ 5 * ε ^ 5 := by positivity
  simp only [neg_mul] at he he0
  nlinarith [mul_le_mul_of_nonneg_left he h5]
theorem T10.proof_4 : T10.stmt_4 := by
  refine ⟨?_, ?_, ?_, ?_⟩ <;> simp [T10.N0, T10.N1, T10.N2log, T10.N3log, Nat.factorial] <;> ring
theorem T10.proof_5 : T10.stmt_5 := by
  intro p hp
  have hp0 : 0 < p := by linarith
  have h := T10.proof_3 (1 / p) (by positivity) (by rw [div_le_one hp0]; exact hp)
  have e : 4 * (2 : ℝ) ^ (-(1 / p)) - 4 * (1 + T10.N0 / p + T10.N1 / p ^ 2 + T10.N2log / p ^ 3 + T10.N3log / p ^ 4)
      = 4 * (T10.twoPow (1 / p) - T10.taylor4 (1 / p)) := by
    unfold T10.twoPow T10.taylor4 T10.N0 T10.N1 T10.N2log T10.N3log
    simp only [Finset.sum_range_succ, Finset.sum_range_zero, Nat.factorial]
    push_cast
    field_simp
    ring
  rw [e, abs_mul, abs_of_pos (by norm_num : (0:ℝ) < 4)]
  have : 4 * Real.log 2 ^ 5 / (120 * p ^ 5) = 4 * (Real.log 2 ^ 5 / 120 * (1 / p) ^ 5) := by
    field_simp
  rw [this]
  exact mul_le_mul_of_nonneg_left h (by norm_num)
