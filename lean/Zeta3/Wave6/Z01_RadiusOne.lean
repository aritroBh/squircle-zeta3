import Mathlib
import Zeta3.Wave3.W02_Holo
import Zeta3.Wave4.X01_RhombusDerivs
import Zeta3.Wave5.Y06_RadiusZero
/-!
# Z01 -- radius exactly one at the square end (`thm:one`)

Target id: Z01.
Paper labels: `thm:one` with its proof (paper/ieee/main.tex, lines 1255-1263), which uses
`lem:holo` (main.tex lines 845-855) and `thm:zero` (main.tex lines 1237-1240). `thm:one` carries
no dagger; the proofs supplement has no separate proof of it, so the intermediate statements
follow the five-line proof printed in main.tex.

Quoted verbatim (main.tex, Theorem thm:one):

  \begin{theorem}[Radius one; uses Theorem~\ref{thm:zero}]\label{thm:one}
  $\limsup_{k\to\infty}|M_k|^{1/k}=1$, and $\varepsilon=1$ (the diamond) is the only singularity of
  $M(\varepsilon)$ on $|\varepsilon|=1$.
  \end{theorem}

Quoted verbatim (main.tex, proof of thm:one):

  Lemma~\ref{lem:holo} gives radius at least $1$, and $\varepsilon=1$ is the only point of
  $|\varepsilon|=1$ outside $\{\operatorname{Re}\varepsilon<1\}$. If the radius were larger than
  $1$, the sum of the series would continue $M$ holomorphically across $\varepsilon=1$, and
  $P(p)=8(1+M(1/p)/p)$ would agree for $p\ge1$ with a function analytic at $p=1$. The Taylor
  series of $P$ at $p=1$ would then converge, which contradicts Theorem~\ref{thm:zero}.

Quoted verbatim (main.tex, Lemma lem:holo):

  $F$ extends to a holomorphic function on $\{\operatorname{Re}\varepsilon<1\}$. Hence
  $P(p)/8=1+\varepsilon M(\varepsilon)$ with $M=-F/2$ holomorphic in $|\varepsilon|<1$, the expansion
  P(p)=8\Bigl(1+\frac{M_0}{p}+\frac{M_1}{p^2}+\frac{M_2}{p^3}+\frac{M_3}{p^4}+\cdots\Bigr)
  converges for every $p>1$, and $M_k=-F^{(k)}(0)/(2\,k!)$.

Quoted verbatim (main.tex, Theorem thm:zero, last sentence):

  In particular the Taylor series of $P$ at $p=1$, in $\eta$ or in $d$, has radius of
  convergence $0$: $P$ is infinitely differentiable but not analytic at $p=1$.

## Conventions and imported objects (nothing is redefined)

* `F` is `W02.FC` (the extension `F(ε) = 2∫_0^{1/2} gC ε t dt` of the proof of `lem:holo`;
  `W02.stmt_5` is its holomorphy on `{Re ε < 1}`), and `M_k` is `W02.Mk k = -F^{(k)}(0)/(2 k!)`
  (a complex number). The only new definition is the function `M = -F/2` of `lem:holo`,
  `Z01.M ε = -W02.FC ε / 2`.
* The square perimeter appears in two imported forms: `W02.Pint p` (eq:Pint verbatim) and
  `X01.P 1 1 p` (the rhombus perimeter with `A = B = 1`, the object of `thm:zero` / `Y06`).
  `stmt_1` proves they agree for every `p`.
* "the series" is `∑_k M_k ε^k`, i.e. `FormalMultilinearSeries.ofScalars ℂ W02.Mk`; its radius is
  stated with Mathlib's `FormalMultilinearSeries.radius`.
* `|M_k|^{1/k}` is the real sequence `‖W02.Mk k‖ ^ (1 / (k : ℝ))` (real `rpow` of a nonnegative
  base; at `k = 0` it is `‖M_0‖^0 = 1`, which does not affect a limsup). As in `Y06.stmt_9`, the
  limsup claim comes with eventual boundedness above, so that `Filter.limsup` is the genuine
  limsup and not a junk value.
* "the only singularity of `M` on `|ε| = 1`" is stated as: (a) at every `e₀` with `‖e₀‖ = 1`,
  `e₀ ≠ 1`, there are `r > 0` and a function `g` holomorphic on the ball `B(e₀, r)` with
  `∑ M_k ε^k = g ε` (as a `HasSum`) for every `ε` in `B(e₀, r) ∩ B(0,1)`; and (b) no such pair
  `(r, g)` exists at `e₀ = 1`.
* `p ≥ 1` link: the paper's "$P(p)=8(1+M(1/p)/p)$ ... for $p\ge1$" is stated over `ℂ` for the
  real perimeter `W02.Pint p` cast to `ℂ` (at `p = 1` it reads `4√2 = 8(1 - (2-√2)/2)`, using the
  integral value `W02.FC 1 = 2 - √2`; `ε = 1` is outside `{Re ε < 1}` but `FC 1` is still the
  honest integral there).

## NOT COVERED

Nothing of `thm:one` is left out. `thm:zero` itself is not restated; its consequence for the square
used by the proof is `stmt_6` (no real-analytic function agrees with `P` on `p ≥ 1` near `p = 1`),
which is to be derived from the imported `Y06.stmt_11` (`A = B = 1`).
-/
-- BEGIN STATEMENTS
namespace Z01

open Filter Metric

/-- The function `M = -F/2` of `lem:holo` (`F` = the extension `W02.FC`). -/
noncomputable def M (ε : ℂ) : ℂ := -W02.FC ε / 2

/-- Link of the imported perimeters: eq:Pint (`W02.Pint`) is the rhombus perimeter `X01.P` with
`A = B = 1`, for every `p`. -/
def stmt_1 : Prop :=
  ∀ p : ℝ, W02.Pint p = X01.P 1 1 p

/-- Proof of `thm:one`: "$P(p)=8(1+M(1/p)/p)$ ... for $p\ge1$", with `M = -F/2`. -/
def stmt_2 : Prop :=
  ∀ p : ℝ, 1 ≤ p → ((W02.Pint p : ℝ) : ℂ) = 8 * (1 + M (1 / (p : ℂ)) / (p : ℂ))

/-- `lem:holo`: "$M=-F/2$ holomorphic in $|\varepsilon|<1$", with Taylor coefficients `M_k`:
`M` is holomorphic on the open unit disk and `∑ M_k ε^k = M(ε)` there. -/
def stmt_3 : Prop :=
  DifferentiableOn ℂ M (ball 0 1) ∧
    ∀ ε : ℂ, ‖ε‖ < 1 → HasSum (fun k : ℕ => W02.Mk k * ε ^ k) (M ε)

/-- Proof of `thm:one`: "Lemma~\ref{lem:holo} gives radius at least $1$", as the radius of
`∑ M_k ε^k` and as `limsup |M_k|^{1/k} ≤ 1` (with eventual boundedness). -/
def stmt_4 : Prop :=
  1 ≤ (FormalMultilinearSeries.ofScalars ℂ W02.Mk).radius ∧
    IsBoundedUnder (· ≤ ·) atTop (fun k : ℕ => ‖W02.Mk k‖ ^ (1 / (k : ℝ))) ∧
    limsup (fun k : ℕ => ‖W02.Mk k‖ ^ (1 / (k : ℝ))) atTop ≤ 1

/-- Proof of `thm:one`: "$\varepsilon=1$ is the only point of $|\varepsilon|=1$ outside
$\{\operatorname{Re}\varepsilon<1\}$", so the sum of the series extends holomorphically to a
neighbourhood of every other point of the unit circle. -/
def stmt_5 : Prop :=
  ∀ e₀ : ℂ, ‖e₀‖ = 1 → e₀ ≠ 1 →
    ∃ r : ℝ, 0 < r ∧ ∃ g : ℂ → ℂ, DifferentiableOn ℂ g (ball e₀ r) ∧
      ∀ ε ∈ ball e₀ r ∩ ball (0 : ℂ) 1, HasSum (fun k : ℕ => W02.Mk k * ε ^ k) (g ε)

/-- Proof of `thm:one`, the step from `thm:zero`: no function analytic at `p = 1` agrees with the
square perimeter `P` for `p ≥ 1` near `1` (otherwise "the Taylor series of $P$ at $p=1$ would then
converge"). -/
def stmt_6 : Prop :=
  ¬ ∃ G : ℝ → ℝ, AnalyticAt ℝ G 1 ∧ ∀ᶠ p in nhdsWithin (1 : ℝ) (Set.Ici 1), G p = X01.P 1 1 p

/-- Proof of `thm:one`: the sum of the series cannot be continued holomorphically across
`ε = 1` ("$\varepsilon=1$ ... is a singularity of $M$"). -/
def stmt_7 : Prop :=
  ¬ ∃ r : ℝ, 0 < r ∧ ∃ g : ℂ → ℂ, DifferentiableOn ℂ g (ball 1 r) ∧
      ∀ ε ∈ ball (1 : ℂ) r ∩ ball (0 : ℂ) 1, HasSum (fun k : ℕ => W02.Mk k * ε ^ k) (g ε)

/-- `thm:one` in full: "$\limsup_{k\to\infty}|M_k|^{1/k}=1$" (with eventual boundedness), i.e.
the series `∑ M_k ε^k` has radius of convergence exactly `1`, "and $\varepsilon=1$ (the diamond)
is the only singularity of $M(\varepsilon)$ on $|\varepsilon|=1$": the sum extends holomorphically
across every other point of the unit circle, and not across `1`. -/
def stmt_8 : Prop :=
  IsBoundedUnder (· ≤ ·) atTop (fun k : ℕ => ‖W02.Mk k‖ ^ (1 / (k : ℝ))) ∧
    limsup (fun k : ℕ => ‖W02.Mk k‖ ^ (1 / (k : ℝ))) atTop = 1 ∧
    (FormalMultilinearSeries.ofScalars ℂ W02.Mk).radius = 1 ∧
    (∀ e₀ : ℂ, ‖e₀‖ = 1 → e₀ ≠ 1 →
      ∃ r : ℝ, 0 < r ∧ ∃ g : ℂ → ℂ, DifferentiableOn ℂ g (ball e₀ r) ∧
        ∀ ε ∈ ball e₀ r ∩ ball (0 : ℂ) 1, HasSum (fun k : ℕ => W02.Mk k * ε ^ k) (g ε)) ∧
    (¬ ∃ r : ℝ, 0 < r ∧ ∃ g : ℂ → ℂ, DifferentiableOn ℂ g (ball 1 r) ∧
      ∀ ε ∈ ball (1 : ℂ) r ∩ ball (0 : ℂ) 1, HasSum (fun k : ℕ => W02.Mk k * ε ^ k) (g ε))

end Z01
-- END STATEMENTS

open Filter Metric Topology Set

theorem Z01.proof_1 : Z01.stmt_1 := by
  intro p
  simp [W02.Pint, X01.P, W02.pIntegrand, X01.integrand]

lemma Z01.FC_one : W02.FC 1 = 2 - ((Real.sqrt 2 : ℝ) : ℂ) := by
  have hg : ∀ t : ℝ, W02.gC 1 t = 2 - ((Real.sqrt 2 : ℝ) : ℂ) := by
    intro t
    have h2 : (2 : ℂ) ^ (1 / 2 : ℂ) = ((Real.sqrt 2 : ℝ) : ℂ) := by
      rw [Real.sqrt_eq_rpow, Complex.ofReal_cpow (by norm_num)]
      push_cast; rfl
    simp only [W02.gC, W02.aC, W02.bC, W02.zC, sub_self, Complex.cpow_zero]
    norm_num
    rw [h2]
  unfold W02.FC
  simp only [hg, intervalIntegral.integral_const, Complex.real_smul]
  push_cast
  ring

theorem Z01.proof_2 : Z01.stmt_2 := by
  intro p hp
  rcases hp.lt_or_eq with hp | rfl
  · obtain ⟨-, -, h⟩ := W02.proof_2 p hp
    have hp0 : (0 : ℝ) < p := by linarith
    obtain ⟨-, h4⟩ := W02.proof_4 (1 / p) (by rw [div_lt_one hp0]; exact hp)
    have hP : W02.Pint p = 8 - 4 * (1 / p) * W02.FR (1 / p) := by linarith
    rw [hP, Z01.M]
    have : (1 / (p : ℂ)) = ((1 / p : ℝ) : ℂ) := by push_cast; ring
    rw [this, h4]
    have hpc : (p : ℂ) ≠ 0 := by exact_mod_cast hp0.ne'
    push_cast
    field_simp
    ring
  · have hP : W02.Pint 1 = 4 * Real.sqrt 2 := by
      simp [W02.Pint, W02.pIntegrand]
      norm_num
    rw [hP, Z01.M]
    simp only [Complex.ofReal_one, div_one, Z01.FC_one]
    push_cast
    ring

lemma Z01.re_sub {e : ℂ} (h : ‖e‖ < 1) : e.re < 1 :=
  lt_of_le_of_lt (Complex.re_le_norm e) h

lemma Z01.M_diff : DifferentiableOn ℂ Z01.M {ε : ℂ | ε.re < 1} := by
  have := W02.proof_5.2
  exact this.neg.div_const 2

theorem Z01.proof_3 : Z01.stmt_3 := by
  refine ⟨Z01.M_diff.mono (fun z hz => Z01.re_sub (by simpa using hz)), fun ε hε => ?_⟩
  have hF : DifferentiableOn ℂ W02.FC (ball 0 1) :=
    W02.proof_5.2.mono (fun z hz => Z01.re_sub (by simpa using hz))
  have h := (Complex.hasSum_taylorSeries_on_ball hF (by simpa using hε)).mul_left (-1 / 2)
  have hM : Z01.M ε = -1 / 2 * W02.FC ε := by rw [Z01.M]; ring
  rw [hM]
  convert h using 1
  funext k
  have hk : (k.factorial : ℂ) ≠ 0 := by exact_mod_cast (Nat.factorial_pos k).ne'
  simp only [W02.Mk, sub_zero, smul_eq_mul]
  field_simp

lemma Z01.radius_ge : 1 ≤ (FormalMultilinearSeries.ofScalars ℂ W02.Mk).radius := by
  refine ENNReal.le_of_forall_nnreal_lt (fun r hr => ?_)
  have hr1 : (r : ℝ) < 1 := by exact_mod_cast hr
  have hs := (Z01.proof_3.2 ((r : ℝ) : ℂ) (by simpa [abs_of_nonneg r.2] using hr1)).summable
  obtain ⟨C, hC⟩ := isBounded_iff_forall_norm_le.mp
    (Metric.isBounded_range_of_tendsto _ hs.tendsto_atTop_zero)
  refine FormalMultilinearSeries.le_radius_of_bound _ C (fun n => ?_)
  rw [FormalMultilinearSeries.ofScalars_norm]
  have := hC _ ⟨n, rfl⟩
  simpa [norm_mul, norm_pow, abs_of_nonneg r.2] using this

lemma Z01.UB_Mk {R : ℝ} (hR : (R.toNNReal : ENNReal) ≤ (FormalMultilinearSeries.ofScalars ℂ W02.Mk).radius) :
    Y06.UB (fun k => ‖W02.Mk k‖) R := by
  intro σ hσ hσR
  have hlt : ((σ.toNNReal : ENNReal)) < (FormalMultilinearSeries.ofScalars ℂ W02.Mk).radius := by
    refine lt_of_lt_of_le ?_ hR
    exact_mod_cast (Real.toNNReal_lt_toNNReal_iff (by linarith)).mpr hσR
  obtain ⟨C, -, hC⟩ := (FormalMultilinearSeries.ofScalars ℂ W02.Mk).norm_mul_pow_le_of_lt_radius hlt
  refine ⟨C, fun m => ?_⟩
  have h := hC m
  rw [FormalMultilinearSeries.ofScalars_norm, Real.coe_toNNReal σ hσ.le] at h
  have hp : 0 < σ ^ m := by positivity
  rw [abs_norm, inv_pow, ← div_eq_mul_inv, le_div_iff₀ hp]
  exact h

lemma Z01.limsup_le :
    IsBoundedUnder (· ≤ ·) atTop (fun k : ℕ => ‖W02.Mk k‖ ^ (1 / (k : ℝ))) ∧
    limsup (fun k : ℕ => ‖W02.Mk k‖ ^ (1 / (k : ℝ))) atTop ≤ 1 := by
  have hUB : Y06.UB (fun k => ‖W02.Mk k‖) 1 := Z01.UB_Mk (by simpa using Z01.radius_ge)
  have hfun : (fun k : ℕ => ‖W02.Mk k‖ ^ (1 / (k : ℝ))) =
      fun k : ℕ => |(fun k => ‖W02.Mk k‖) k| ^ (1 / (k : ℝ)) := by
    funext k; simp
  rw [hfun]
  have hcob : IsCoboundedUnder (· ≤ ·) atTop
      (fun m : ℕ => |(fun k => ‖W02.Mk k‖) m| ^ (1 / (m : ℝ))) :=
    isCoboundedUnder_le_of_le atTop (Y06.root_nonneg _)
  obtain ⟨C, hC⟩ := hUB (1 / 2) (by norm_num) (by norm_num)
  refine ⟨isBoundedUnder_of_eventually_le (Y06.ev_root_le (by norm_num) hC
    (lt_add_one _)), ?_⟩
  apply le_of_forall_gt_imp_ge_of_dense
  intro y hy
  have hy0 : 0 < y := lt_trans one_pos hy
  obtain ⟨σ, hσ1, hσ2⟩ : ∃ σ, y⁻¹ < σ ∧ σ < 1 := exists_between (inv_lt_one_of_one_lt₀ hy)
  have hσ0 : 0 < σ := lt_trans (inv_pos.mpr hy0) hσ1
  obtain ⟨C, hC⟩ := hUB σ hσ0 hσ2
  exact limsup_le_of_le hcob (Y06.ev_root_le hσ0 hC ((inv_lt_comm₀ hy0 hσ0).mp hσ1))

theorem Z01.proof_4 : Z01.stmt_4 :=
  ⟨Z01.radius_ge, Z01.limsup_le.1, Z01.limsup_le.2⟩

lemma Z01.re_lt_one {e₀ : ℂ} (h1 : ‖e₀‖ = 1) (hne : e₀ ≠ 1) : e₀.re < 1 := by
  have hle : e₀.re ≤ 1 := h1 ▸ Complex.re_le_norm e₀
  rcases hle.lt_or_eq with h | h
  · exact h
  · exfalso; apply hne
    have hsq : e₀.re ^ 2 + e₀.im ^ 2 = 1 := by
      have := Complex.sq_norm e₀
      rw [h1, Complex.normSq_apply] at this
      nlinarith [this]
    have him : e₀.im = 0 := by rw [h] at hsq; nlinarith [sq_nonneg e₀.im]
    exact Complex.ext h him

theorem Z01.proof_5 : Z01.stmt_5 := by
  intro e₀ h1 hne
  have hre := Z01.re_lt_one h1 hne
  refine ⟨1 - e₀.re, by linarith, Z01.M, Z01.M_diff.mono (fun z hz => ?_), fun ε hε => ?_⟩
  · have hz' : ‖z - e₀‖ < 1 - e₀.re := by simpa [dist_eq_norm] using hz
    have := Complex.re_le_norm (z - e₀)
    simp only [Complex.sub_re] at this
    show z.re < 1
    linarith
  · exact Z01.proof_3.2 ε (by simpa using hε.2)

theorem Z01.proof_6 : Z01.stmt_6 := by
  rintro ⟨G, hG, hEq⟩
  have hPcd := Y06.P_cd (A := 1) (B := 1) one_pos one_pos
  have hiter : ∀ n : ℕ, iteratedDeriv n G 1 = iteratedDeriv n (X01.P 1 1) 1 := by
    intro n
    have hG1 : G 1 = X01.P 1 1 1 := hEq.self_of_nhdsWithin self_mem_Ici
    rw [← iteratedDerivWithin_eq_iteratedDeriv (uniqueDiffOn_Ici 1) hG.contDiffAt self_mem_Ici,
      ← iteratedDerivWithin_eq_iteratedDeriv (uniqueDiffOn_Ici 1)
        (hPcd.of_le (mod_cast le_top)) self_mem_Ici]
    exact Filter.EventuallyEq.iteratedDerivWithin_eq (f := G) (g := X01.P 1 1) hEq hG1
  obtain ⟨C, s, hs, hC⟩ := Y06.analytic_coeff hG
  have hd : (1 / (2 * s) : ℝ) ≠ 0 := by positivity
  apply Y06.P_not_summable (A := 1) (B := 1) one_pos one_pos hd
  refine Summable.of_norm_bounded ((summable_geometric_two).mul_left C) (fun n => ?_)
  rw [Real.norm_eq_abs, abs_mul, ← hiter n, abs_pow, abs_of_pos (by positivity : (0:ℝ) < 1 / (2 * s))]
  calc |iteratedDeriv n G 1 / (n.factorial : ℝ)| * (1 / (2 * s)) ^ n
      ≤ C * s ^ n * (1 / (2 * s)) ^ n := by gcongr; exact hC n
    _ = C * (1 / 2) ^ n := by
      rw [mul_assoc, ← mul_pow]; congr 2; field_simp

theorem Z01.proof_7 : Z01.stmt_7 := by
  rintro ⟨r, hr, g, hg, hsum⟩
  apply Z01.proof_6
  have hgA : AnalyticAt ℂ g 1 := hg.analyticAt (isOpen_ball.mem_nhds (mem_ball_self hr))
  set f : ℂ → ℂ := fun w => 8 * (1 + g (1 / w) / w) with hf
  have hinv : AnalyticAt ℂ (fun w : ℂ => 1 / w) 1 := by
    simpa [one_div] using analyticAt_inv (𝕜 := ℂ) (z := (1 : ℂ)) one_ne_zero
  have hfA : AnalyticAt ℂ f 1 := by
    have h1 : AnalyticAt ℂ g ((fun w : ℂ => 1 / w) 1) := by simpa using hgA
    have h2 : AnalyticAt ℂ (fun w : ℂ => g (1 / w)) 1 := h1.comp hinv
    have h3 : AnalyticAt ℂ (fun w : ℂ => g (1 / w) / w) 1 :=
      h2.div analyticAt_id one_ne_zero
    exact analyticAt_const.mul (analyticAt_const.add h3)
  have hGA : AnalyticAt ℝ (fun p : ℝ => (f p).re) 1 := by
    simpa using hfA.re_ofReal (x := 1)
  have hright : ∀ p : ℝ, 1 < p → p < 1 + r → (f p).re = X01.P 1 1 p := by
    intro p hp1 hp2
    have hp0 : 0 < p := by linarith
    have hn1 : ‖(1 / (p : ℂ))‖ < 1 := by
      rw [norm_div, norm_one, Complex.norm_real, Real.norm_eq_abs, abs_of_pos hp0]
      exact (div_lt_one hp0).mpr hp1
    have hb : 1 / (p : ℂ) ∈ ball (1 : ℂ) r := by
      rw [mem_ball, dist_eq_norm]
      have : 1 / (p : ℂ) - 1 = ((1 / p - 1 : ℝ) : ℂ) := by push_cast; ring
      rw [this, Complex.norm_real, Real.norm_eq_abs]
      have h1 : 1 / p - 1 < 0 := by rw [sub_neg, div_lt_one hp0]; exact hp1
      rw [abs_of_neg h1]
      have h1p : 1 / p ≤ 1 := by rw [div_le_one hp0]; exact hp1.le
      have : (p - 1) / p ≤ p - 1 := div_le_self (by linarith) hp1.le
      have : 1 - 1 / p = (p - 1) / p := by field_simp
      linarith
    have hg_eq : g (1 / (p : ℂ)) = Z01.M (1 / (p : ℂ)) :=
      (hsum _ ⟨hb, by simpa using hn1⟩).unique (Z01.proof_3.2 _ hn1)
    have h2 := Z01.proof_2 p hp1.le
    simp only [hf, hg_eq, ← h2, Complex.ofReal_re]
    exact Z01.proof_1 p
  have hPc : ContinuousAt (X01.P 1 1) 1 :=
    (Y06.P_cd (A := 1) (B := 1) one_pos one_pos).continuousAt
  have hev : (fun p : ℝ => (f p).re) =ᶠ[𝓝[>] (1 : ℝ)] X01.P 1 1 := by
    filter_upwards [Ioo_mem_nhdsGT (by linarith : (1 : ℝ) < 1 + r)] with p hp
    exact hright p hp.1 hp.2
  have h1 : (f 1).re = X01.P 1 1 1 :=
    tendsto_nhds_unique (hGA.continuousAt.tendsto.mono_left nhdsWithin_le_nhds)
      ((hPc.tendsto.mono_left nhdsWithin_le_nhds).congr' hev.symm)
  refine ⟨fun p => (f p).re, hGA, ?_⟩
  filter_upwards [Ico_mem_nhdsGE (by linarith : (1 : ℝ) < 1 + r)] with p hp
  rcases hp.1.lt_or_eq with h | h
  · exact hright p h hp.2
  · subst h; simpa using h1

lemma Z01.radius_le : (FormalMultilinearSeries.ofScalars ℂ W02.Mk).radius ≤ 1 := by
  by_contra hlt
  push Not at hlt
  set ps := FormalMultilinearSeries.ofScalars ℂ W02.Mk with hps
  obtain ⟨R, hR1, hR2⟩ := ENNReal.lt_iff_exists_nnreal_btwn.mp hlt
  have hR1' : (1 : ℝ) < R := by exact_mod_cast hR1
  have hpos : 0 < ps.radius := lt_trans (by norm_num) hlt
  have hP := ps.hasFPowerSeriesOnBall hpos
  apply Z01.proof_7
  have hsub : ∀ ε ∈ ball (1 : ℂ) (R - 1), ε ∈ Metric.eball (0 : ℂ) ps.radius := by
    intro ε hε
    have h1 : ‖ε‖ < R := by
      have := norm_sub_norm_le ε 1
      rw [mem_ball, dist_eq_norm] at hε
      simp at this
      linarith
    have : ε ∈ Metric.eball (0 : ℂ) (R : ENNReal) := by
      rw [Metric.eball_coe, mem_ball, dist_zero_right]; exact h1
    exact Metric.eball_subset_eball hR2.le this
  refine ⟨R - 1, by linarith, ps.sum, hP.differentiableOn.mono hsub, fun ε hε => ?_⟩
  have h := hP.hasSum (hsub ε hε.1)
  rw [zero_add] at h
  convert h using 1
  funext k
  rw [hps, FormalMultilinearSeries.ofScalars_apply_eq, smul_eq_mul]

lemma Z01.noBetter : Y06.NoBetter (fun k => ‖W02.Mk k‖) 1 := by
  rintro ⟨K, τ, hτ0, hτ1, hK⟩
  have hb : ((τ⁻¹.toNNReal : NNReal) : ENNReal) ≤ (FormalMultilinearSeries.ofScalars ℂ W02.Mk).radius := by
    refine FormalMultilinearSeries.le_radius_of_bound _ K (fun n => ?_)
    rw [FormalMultilinearSeries.ofScalars_norm, Real.coe_toNNReal _ (inv_pos.mpr hτ0).le]
    have h := hK n
    rw [abs_norm] at h
    calc ‖W02.Mk n‖ * τ⁻¹ ^ n ≤ K * τ ^ n * τ⁻¹ ^ n := by gcongr
      _ = K := by rw [mul_assoc, ← mul_pow, mul_inv_cancel₀ hτ0.ne', one_pow, mul_one]
  have h1 : (1 : ℝ) < τ⁻¹ := by rw [one_div] at hτ1; simpa using one_lt_inv₀ hτ0 |>.mpr (by simpa using hτ1)
  have : (1 : ENNReal) < ((τ⁻¹.toNNReal : NNReal) : ENNReal) := by
    have : (1 : NNReal) < τ⁻¹.toNNReal := by
      rw [← Real.toNNReal_one]; exact (Real.toNNReal_lt_toNNReal_iff (by linarith)).mpr h1
    exact_mod_cast this
  exact absurd (lt_of_lt_of_le this hb) (not_lt.mpr Z01.radius_le)

theorem Z01.proof_8 : Z01.stmt_8 := by
  have hrad : (FormalMultilinearSeries.ofScalars ℂ W02.Mk).radius = 1 :=
    le_antisymm Z01.radius_le Z01.radius_ge
  have hUB : Y06.UB (fun k => ‖W02.Mk k‖) 1 := Z01.UB_Mk (by simp [hrad])
  have hL := Y06.limsup_of one_pos hUB Z01.noBetter
  have hfun : (fun k : ℕ => ‖W02.Mk k‖ ^ (1 / (k : ℝ))) =
      fun k : ℕ => |(fun k => ‖W02.Mk k‖) k| ^ (1 / (k : ℝ)) := by
    funext k; simp
  unfold Z01.stmt_8
  rw [hfun]
  exact ⟨hL.1, by rw [hL.2, div_one], hrad, Z01.proof_5, Z01.proof_7⟩
