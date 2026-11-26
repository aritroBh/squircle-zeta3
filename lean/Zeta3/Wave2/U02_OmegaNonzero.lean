import Mathlib
/-!
# U02 — `ω_n ≠ 0` for every odd `n ≥ 3` (last sentence of Proposition `prop:an`)

Target id: U02.
Paper label: `prop:an`, paper/ieee/main.tex lines 1163–1174 (statement) and the
argument in the paragraph at line 1178 (the proofs supplement has no separate proof).

Claim quoted verbatim (main.tex lines 1170–1173):

  \[
  \omega_n=\frac{4(-1)^n}{n}\bigl(1-2^{1-n}V_n\bigr),\qquad
  V_n=\frac{d^n}{dt^n}\sqrt{\tfrac{1+e^{4t}}{2}}\,\Big|_{t=0},
  \]
  so $\omega_3=1$, $\omega_5=-\frac34$, and $\omega_n\ne0$ for every odd $n\ge3$.

Argument quoted verbatim (main.tex line 1178):

  Is $\omega_n$ ever $0$? No. Write $\sqrt{(1+e^{4t})/2}=e^t\sqrt{\cosh2t}$ and
  $y=(\cosh2t-1)/4$. Then $y(0)=0$ and the derivatives of $y$ at $0$ are integers
  ($2^{m-2}$ for even $m\ge2$, and $0$ otherwise). Also
  $\sqrt{\cosh2t}=\sqrt{1+4y}=1+\sum_{j\ge1}(-1)^{j-1}\frac2j\binom{2j-2}{j-1}y^j$ has even
  integer coefficients. Sums and products of functions with integer derivatives at $0$
  have the same property (Leibniz rule), and $y^j$ vanishes to order $2j$, so each
  derivative at $0$ involves only finitely many $j$. Hence every derivative $\nu_k$ of
  $\sqrt{\cosh2t}$ at $0$ is an integer, with $\nu_0=1$ and $\nu_k$ even for $k\ge1$.
  Then $V_n=\sum_k\binom nk\nu_k$ is odd, so $V_n\ne2^{n-1}$ and $\omega_n\ne0$ for every
  odd $n\ge3$. (The formula gives $\omega_1=0$, but $\zeta(1)$ is not a value.)

Formalization choices:
* `f`, `V`, `ω` are defined exactly as in Wave1/T03_Omega.lean (re-declared here, not
  imported): `V n = iteratedDeriv n f 0` with `f t = √((1 + e^{4t})/2)`; the radicand is
  `> 0`, so `Real.sqrt`'s convention plays no role and `f` is smooth.
  `ω n` uses the integer power `(2:ℝ) ^ ((1:ℤ) - n)` (no truncated ℕ subtraction);
  division by `(n:ℝ)` is only used for `n ≥ 1`.
* `g t = √(cosh 2t)` (radicand `≥ 1`), `ν k = iteratedDeriv k g 0`.
* The binomial series is stated with the reindexing `j ↦ j+1`, so the coefficient
  `(-1)^{j-1}(2/j)C(2j-2, j-1)` for `j ≥ 1` becomes `(-1)^j (2/(j+1)) C(2j, j)` for
  `j ≥ 0` (avoids ℕ subtraction). The paper gives no domain for the expansion; it is
  stated as a `HasSum` on `|y| < 1/4`, its disc of convergence (an added hypothesis,
  needed for the identity to be meaningful).
* "Integer" is `∃ m : ℤ, x = m`; "odd" is `∃ m : ℤ, x = 2m+1`; "even" is `∃ m : ℤ, x = 2m`.

NOT COVERED:
* "Sums and products of functions with integer derivatives at 0 have the same property"
  and "y^j vanishes to order 2j" as general lemmas (proof devices, not claims about the
  paper's objects; their consequence, `ν_k` integer and even, is stmt_5).
* The identification of `ω_n` as the coefficient of `ζ(n)` in `a_n` (needs the full
  perimeter expansion; outside this target).
* `ω_3 = 1`, `ω_5 = -3/4` (covered by Wave1 T03).
-/
-- BEGIN STATEMENTS
namespace U02

/-- The function `t ↦ √((1 + e^{4t})/2)` whose derivatives at `0` define `V_n`. -/
noncomputable def f (t : ℝ) : ℝ := Real.sqrt ((1 + Real.exp (4 * t)) / 2)

/-- `V_n = (d^n/dt^n) √((1+e^{4t})/2) |_{t=0}`. -/
noncomputable def V (n : ℕ) : ℝ := iteratedDeriv n f 0

/-- `ω_n = (4(-1)^n / n) (1 - 2^{1-n} V_n)`. -/
noncomputable def ω (n : ℕ) : ℝ :=
  4 * (-1 : ℝ) ^ n / (n : ℝ) * (1 - (2 : ℝ) ^ ((1 : ℤ) - (n : ℤ)) * V n)

/-- `y = (cosh 2t - 1)/4`. -/
noncomputable def y (t : ℝ) : ℝ := (Real.cosh (2 * t) - 1) / 4

/-- `√(cosh 2t)`. -/
noncomputable def g (t : ℝ) : ℝ := Real.sqrt (Real.cosh (2 * t))

/-- `ν_k`, the `k`-th derivative of `√(cosh 2t)` at `0`. -/
noncomputable def ν (k : ℕ) : ℝ := iteratedDeriv k g 0

/-- The coefficient of `y^{j+1}` in `√(1+4y)`, i.e. the paper's
`(-1)^{j-1} (2/j) C(2j-2, j-1)` at index `j+1`. -/
noncomputable def c (j : ℕ) : ℝ := (-1 : ℝ) ^ j * (2 / ((j : ℝ) + 1)) * (Nat.choose (2 * j) j : ℝ)

/-- Line 1178: "Write $\sqrt{(1+e^{4t})/2}=e^t\sqrt{\cosh2t}$". -/
def stmt_1 : Prop := ∀ t : ℝ, f t = Real.exp t * g t

/-- Line 1178: "$y(0)=0$ and the derivatives of $y$ at $0$ are integers
($2^{m-2}$ for even $m\ge2$, and $0$ otherwise)". -/
def stmt_2 : Prop :=
  y 0 = 0 ∧ ∀ m : ℕ, iteratedDeriv m y 0 = if Even m ∧ 2 ≤ m then (2 : ℝ) ^ (m - 2) else 0

/-- Line 1178: "$\sqrt{1+4y}=1+\sum_{j\ge1}(-1)^{j-1}\frac2j\binom{2j-2}{j-1}y^j$"
(on the disc of convergence `|y| < 1/4`, reindexed `j ↦ j+1`). -/
def stmt_3 : Prop :=
  ∀ x : ℝ, |x| < 1 / 4 →
    HasSum (fun j : ℕ => c j * x ^ (j + 1)) (Real.sqrt (1 + 4 * x) - 1)

/-- Line 1178: the coefficients of that series "has even integer coefficients". -/
def stmt_4 : Prop := ∀ j : ℕ, ∃ m : ℤ, c j = 2 * m

/-- Line 1178: "every derivative $\nu_k$ of $\sqrt{\cosh2t}$ at $0$ is an integer, with
$\nu_0=1$ and $\nu_k$ even for $k\ge1$". -/
def stmt_5 : Prop :=
  ν 0 = 1 ∧ ∀ k : ℕ, 1 ≤ k → ∃ m : ℤ, ν k = 2 * m

/-- Line 1178: "$V_n=\sum_k\binom nk\nu_k$". -/
def stmt_6 : Prop :=
  ∀ n : ℕ, V n = ∑ k ∈ Finset.range (n + 1), (Nat.choose n k : ℝ) * ν k

/-- Line 1178 / line 2150 ("$V_n$ odd for every $n$"): `V_n` is an odd integer. -/
def stmt_7 : Prop := ∀ n : ℕ, ∃ m : ℤ, V n = 2 * m + 1

/-- Line 1178: "so $V_n\ne2^{n-1}$ ... for every odd $n\ge3$". -/
def stmt_8 : Prop := ∀ n : ℕ, Odd n → 3 ≤ n → V n ≠ (2 : ℝ) ^ (n - 1)

/-- Line 1173 (and 1178): "$\omega_n\ne0$ for every odd $n\ge3$". -/
def stmt_9 : Prop := ∀ n : ℕ, Odd n → 3 ≤ n → ω n ≠ 0

/-- Line 1178: "(The formula gives $\omega_1=0$ ...)". -/
def stmt_10 : Prop := ω 1 = 0

end U02
-- END STATEMENTS

namespace U02

lemma cosh2_deriv (m : ℕ) :
    iteratedDeriv m (fun t : ℝ => Real.cosh (2 * t)) 0 = if Even m then (2 : ℝ) ^ m else 0 := by
  rw [iteratedDeriv_comp_const_mul Real.contDiff_cosh]
  simp only [mul_zero]
  rcases Nat.even_or_odd' m with ⟨n, rfl | rfl⟩
  · simp
  · simp

lemma g_sq (t : ℝ) : g t ^ 2 = Real.cosh (2 * t) := by
  unfold g
  rw [Real.sq_sqrt (by have := Real.cosh_pos (2*t); linarith)]

lemma g_contDiff : ContDiff ℝ ⊤ g := by
  unfold g
  apply ContDiff.sqrt
  · exact Real.contDiff_cosh.comp (contDiff_const.mul contDiff_id)
  · intro t; exact (Real.cosh_pos _).ne'

end U02

theorem U02.proof_1 : U02.stmt_1 := by
  intro t
  unfold U02.f U02.g
  have h1 : (1 + Real.exp (4 * t)) / 2 = Real.exp t ^ 2 * Real.cosh (2 * t) := by
    rw [Real.cosh_eq, sq, ← Real.exp_add]
    have e4 : Real.exp (4 * t) = Real.exp (2*t) * Real.exp (2*t) := by rw [← Real.exp_add]; ring_nf
    have e2 : Real.exp (t + t) = Real.exp (2*t) := by ring_nf
    rw [e2, e4, Real.exp_neg]
    have := Real.exp_pos (2*t)
    field_simp
    ring
  rw [h1, Real.sqrt_mul (by positivity), Real.sqrt_sq (Real.exp_pos t).le]

theorem U02.proof_2 : U02.stmt_2 := by
  refine ⟨by simp [U02.y], fun m => ?_⟩
  have hy : U02.y = fun t => (Real.cosh (2 * t) - 1) / 4 := rfl
  rw [hy, iteratedDeriv_div_const]
  have hc : ContDiff ℝ ⊤ (fun t : ℝ => Real.cosh (2 * t)) :=
    Real.contDiff_cosh.comp (contDiff_const.mul contDiff_id)
  rcases Nat.eq_zero_or_pos m with rfl | hm
  · simp
  · rw [iteratedDeriv_fun_sub (hc.contDiffAt.of_le le_top) contDiffAt_const, cosh2_deriv,
      iteratedDeriv_const]
    have hm0 : m ≠ 0 := hm.ne'
    by_cases he : Even m
    · have h2 : 2 ≤ m := by
        obtain ⟨k, hk⟩ := he; omega
      simp only [he, h2, true_and, ite_true, hm0, ite_false, sub_zero]
      obtain ⟨k, hk⟩ : ∃ k, m = k + 2 := ⟨m - 2, by omega⟩
      subst hk
      simp [pow_add]
      ring
    · simp [he, hm0]

namespace U02

lemma nu_zero : ν 0 = 1 := by
  simp [ν, g]

lemma nu_conv (n : ℕ) :
    ∑ i ∈ Finset.range (n + 1), (n.choose i : ℝ) * ν i * ν (n - i)
      = if Even n then (2 : ℝ) ^ n else 0 := by
  have hg : ContDiffAt ℝ n g 0 := (g_contDiff.of_le le_top).contDiffAt
  have := iteratedDeriv_mul (n := n) (x := (0:ℝ)) hg hg
  have e : g * g = fun t : ℝ => Real.cosh (2 * t) := by
    funext t; simp only [Pi.mul_apply]; rw [← sq, g_sq]
  rw [e, cosh2_deriv] at this
  exact this.symm

lemma sum_four {ι : Type*} (s : Finset ι) (f : ι → ℝ) (h : ∀ i ∈ s, ∃ N : ℤ, f i = 4 * N) :
    ∃ N : ℤ, ∑ i ∈ s, f i = 4 * N := by
  classical
  induction s using Finset.induction_on with
  | empty => exact ⟨0, by simp⟩
  | insert a s ha ih =>
    obtain ⟨N, hN⟩ := h a (Finset.mem_insert_self _ _)
    obtain ⟨M, hM⟩ := ih (fun i hi => h i (Finset.mem_insert_of_mem hi))
    exact ⟨N + M, by rw [Finset.sum_insert ha, hN, hM]; push_cast; ring⟩

lemma rhs4 (n : ℕ) (hn : 1 ≤ n) : ∃ R : ℤ, (if Even n then (2 : ℝ) ^ n else 0) = 4 * R := by
  by_cases he : Even n
  · obtain ⟨k, hk⟩ : ∃ k, n = k + 2 := by
      obtain ⟨j, hj⟩ := he; exact ⟨n - 2, by omega⟩
    refine ⟨2 ^ k, ?_⟩
    subst hk
    simp [he, pow_add]; ring
  · exact ⟨0, by simp [he]⟩

lemma nu_even : ∀ k : ℕ, 1 ≤ k → ∃ m : ℤ, ν k = 2 * m := by
  intro k
  induction k using Nat.strong_induction_on with
  | _ k ih =>
    intro hk
    obtain ⟨k, rfl⟩ : ∃ j, k = j + 1 := ⟨k - 1, by omega⟩
    have h := nu_conv (k + 1)
    rw [Finset.sum_range_succ, Finset.sum_range_succ'] at h
    simp only [Nat.choose_zero_right, Nat.choose_self, Nat.sub_zero, Nat.sub_self, nu_zero,
      Nat.cast_one] at h
    obtain ⟨N, hN⟩ := sum_four (Finset.range k)
      (fun i => ((k + 1).choose (i + 1) : ℝ) * ν (i + 1) * ν (k + 1 - (i + 1))) (by
        intro i hi
        rw [Finset.mem_range] at hi
        obtain ⟨a, ha⟩ := ih (i + 1) (by omega) (by omega)
        obtain ⟨b, hb⟩ := ih (k + 1 - (i + 1)) (by omega) (by omega)
        refine ⟨(k + 1).choose (i + 1) * a * b, ?_⟩
        rw [ha, hb]; push_cast; ring)
    obtain ⟨R, hR⟩ := rhs4 (k + 1) hk
    refine ⟨R - N, ?_⟩
    push_cast
    rw [hN, hR] at h
    linarith

lemma nu_one : ν 1 = 0 := by
  have h := nu_conv 1
  simp [Finset.sum_range_succ, nu_zero] at h
  linarith

end U02

theorem U02.proof_5 : U02.stmt_5 := ⟨U02.nu_zero, U02.nu_even⟩

theorem U02.proof_4 : U02.stmt_4 := by
  intro j
  obtain ⟨q, hq⟩ := Nat.succ_dvd_centralBinom j
  refine ⟨(-1) ^ j * q, ?_⟩
  unfold U02.c
  have h1 : (Nat.choose (2 * j) j : ℝ) = ((j : ℝ) + 1) * q := by
    have : Nat.choose (2 * j) j = (j + 1) * q := by
      rw [← Nat.centralBinom_eq_two_mul_choose]; exact hq
    rw [this]; push_cast; ring
  rw [h1]
  have : ((j : ℝ) + 1) ≠ 0 := by positivity
  push_cast
  field_simp

namespace U02

lemma V_eq : V = fun n => iteratedDeriv n (g * Real.exp) 0 := by
  funext n
  unfold V
  congr 1
  funext t
  simp only [Pi.mul_apply]
  rw [mul_comm]
  exact proof_1 t

lemma V_sum (n : ℕ) : V n = ∑ k ∈ Finset.range (n + 1), (Nat.choose n k : ℝ) * ν k := by
  rw [V_eq]
  simp only
  rw [iteratedDeriv_mul (g_contDiff.of_le le_top).contDiffAt Real.contDiff_exp.contDiffAt]
  refine Finset.sum_congr rfl fun i _ => ?_
  have : iteratedDeriv (n - i) Real.exp 0 = 1 := by
    rw [iteratedDeriv_eq_iterate, Real.iter_deriv_exp]; simp
  rw [this]; unfold ν; ring

lemma sum_two {ι : Type*} (s : Finset ι) (f : ι → ℝ) (h : ∀ i ∈ s, ∃ N : ℤ, f i = 2 * N) :
    ∃ N : ℤ, ∑ i ∈ s, f i = 2 * N := by
  classical
  induction s using Finset.induction_on with
  | empty => exact ⟨0, by simp⟩
  | insert a s ha ih =>
    obtain ⟨N, hN⟩ := h a (Finset.mem_insert_self _ _)
    obtain ⟨M, hM⟩ := ih (fun i hi => h i (Finset.mem_insert_of_mem hi))
    exact ⟨N + M, by rw [Finset.sum_insert ha, hN, hM]; push_cast; ring⟩

lemma V_odd (n : ℕ) : ∃ m : ℤ, V n = 2 * m + 1 := by
  rw [V_sum, Finset.sum_range_succ']
  obtain ⟨N, hN⟩ := sum_two (Finset.range n) (fun i => ((n.choose (i + 1) : ℕ) : ℝ) * ν (i + 1)) (by
    intro i _
    obtain ⟨a, ha⟩ := nu_even (i + 1) (by omega)
    exact ⟨n.choose (i + 1) * a, by rw [ha]; push_cast; ring⟩)
  exact ⟨N, by rw [hN]; simp [nu_zero]⟩

end U02

theorem U02.proof_6 : U02.stmt_6 := U02.V_sum
theorem U02.proof_7 : U02.stmt_7 := U02.V_odd

theorem U02.proof_8 : U02.stmt_8 := by
  intro n hn h3 h
  obtain ⟨m, hm⟩ := U02.V_odd n
  have h1 : ((2:ℤ) ^ (n - 1) : ℤ) = 2 * m + 1 := by
    have : ((2:ℝ) ^ (n - 1)) = 2 * m + 1 := by rw [← h, hm]
    exact_mod_cast this
  have h2 : (2:ℤ) ∣ 2 ^ (n - 1) := dvd_pow_self 2 (by omega)
  omega

theorem U02.proof_9 : U02.stmt_9 := by
  intro n hn h3
  have h8 := U02.proof_8 n hn h3
  have hn0 : (n : ℝ) ≠ 0 := by positivity
  unfold U02.ω
  refine mul_ne_zero (div_ne_zero (mul_ne_zero (by norm_num) (pow_ne_zero _ (by norm_num))) hn0) ?_
  intro h0
  apply h8
  have e : (2 : ℝ) ^ ((1 : ℤ) - (n : ℤ)) = ((2 : ℝ) ^ (n - 1))⁻¹ := by
    have : (1 : ℤ) - (n : ℤ) = -((n - 1 : ℕ) : ℤ) := by omega
    rw [this, zpow_neg, zpow_natCast]
  rw [e] at h0
  have hp : (2 : ℝ) ^ (n - 1) ≠ 0 := by positivity
  field_simp at h0
  linarith

theorem U02.proof_10 : U02.stmt_10 := by
  have h : U02.V 1 = 1 := by
    rw [U02.V_sum]; simp [Finset.sum_range_succ, U02.nu_zero, U02.nu_one]
  show U02.ω 1 = 0
  unfold U02.ω
  rw [h]; norm_num

namespace U02

lemma choose_ratio (r : ℝ) (k : ℕ) :
    ((k : ℝ) + 1) * Ring.choose r (k + 1) = Ring.choose r k * (r - k) := by
  have h := Ring.choose_smul_choose r (n := k + 1) (k := k) (by omega)
  rw [Nat.choose_succ_self_right, Nat.add_sub_cancel_left] at h
  simp only [Ring.choose_one_right, nsmul_eq_mul] at h
  rw [← h]; push_cast; ring

lemma c_rec (j : ℕ) : ((j : ℝ) + 2) * c (j + 1) = -(4 * j + 2) * c j := by
  have h := Nat.succ_mul_centralBinom_succ j
  have h' : ((j : ℝ) + 1) * (Nat.choose (2 * (j + 1)) (j + 1) : ℝ)
      = 2 * (2 * j + 1) * (Nat.choose (2 * j) j : ℝ) := by
    have := congrArg (fun x : ℕ => (x : ℝ)) h
    simpa [Nat.centralBinom_eq_two_mul_choose] using this
  unfold c
  have h1 : ((j : ℝ) + 1) ≠ 0 := by positivity
  have h2 : ((j : ℝ) + 1 + 1) ≠ 0 := by positivity
  push_cast
  field_simp
  rw [pow_succ]
  linear_combination (-((j : ℝ) + 2) * (-1 : ℝ) ^ j) * h'

lemma choose_half (j : ℕ) : Ring.choose (1 / 2 : ℝ) (j + 1) * 4 ^ (j + 1) = c j := by
  induction j with
  | zero => simp [c]; norm_num
  | succ j ih =>
    have hr := choose_ratio (1 / 2 : ℝ) (j + 1)
    have hc := c_rec j
    have h2 : ((j : ℝ) + 2) ≠ 0 := by positivity
    apply mul_left_cancel₀ h2
    rw [hc, ← ih]
    push_cast at hr
    have : ((j : ℝ) + 2) * Ring.choose (1 / 2 : ℝ) (j + 1 + 1) = Ring.choose (1 / 2 : ℝ) (j + 1) * (1 / 2 - (j + 1)) := by
      rw [← hr]; ring
    calc ((j : ℝ) + 2) * (Ring.choose (1 / 2 : ℝ) (j + 1 + 1) * 4 ^ (j + 1 + 1))
        = (((j : ℝ) + 2) * Ring.choose (1 / 2 : ℝ) (j + 1 + 1)) * 4 ^ (j + 1 + 1) := by ring
      _ = _ := by rw [this]; ring

end U02

theorem U02.proof_3 : U02.stmt_3 := by
  intro x hx
  have hx' : |4 * x| < 1 := by
    rw [abs_lt] at hx ⊢; constructor <;> linarith [hx.1, hx.2]
  have hb : (4 * x) ∈ Metric.eball (0 : ℝ) 1 := by
    rw [Metric.mem_eball, edist_dist, Real.dist_eq, sub_zero]
    exact ENNReal.ofReal_lt_one.mpr hx'
  have h := (Real.one_add_rpow_hasFPowerSeriesOnBall_zero (a := (1 / 2 : ℝ))).hasSum hb
  simp only [binomialSeries_apply, List.ofFn_const, List.prod_replicate, zero_add,
    smul_eq_mul] at h
  have h1 := (hasSum_nat_add_iff' 1).mpr h
  simp only [Finset.sum_range_one, Ring.choose_zero_right, pow_zero, mul_one] at h1
  have e : (1 + 4 * x) ^ (1 / 2 : ℝ) = Real.sqrt (1 + 4 * x) := (Real.sqrt_eq_rpow _).symm
  rw [e] at h1
  convert h1 using 1
  funext j
  rw [mul_pow, ← U02.choose_half]
  ring
