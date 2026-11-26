import Mathlib
/-!
# W01 -- the length of the graph of `f^n` (Ahmed's Open Problem 3)

Target id: W01.
Paper label: `thm:fn` (paper/ieee/main.tex, lines 1811-1816).
Detailed proof: paper/ieee/proofs/sec_g3_fn.tex (setting line 5-6, theorem lines 9-14, bound
`g3:eq:phibound` lines 18-23, Lemma `g3:lem:layer` lines 27-33, Lemma `g3:lem:moments`
lines 35-48, localisation lines 50-53, Lemma `g3:lem:shape` lines 55-68, Lemma `g3:lem:window`
lines 70-81, proof of the theorem lines 83-97).

Quoted verbatim (main.tex 1811-1816):

  \begin{theorem}[Length of the graph of $f^n$$^\dagger$]\label{thm:fn}
  Let $f\in C^1[0,1]$ be nondecreasing, $f(0)=0$, $f(1)=1$, $f'(1)>0$, and $f\in C^2$ near $x=1$. Then
  \[
  2-L_n=\frac{\log n+\log f'(1)+1-\log2}{f'(1)\,(n-1)}+O\Bigl(\frac{\log^2n}{n^2}\Bigr).
  \]
  \end{theorem}

Quoted verbatim (sec_g3_fn.tex 6, the definition of `L_n`):

  For $f$ as below and $n>1$, $L_n$ is the length of the graph of $y_n=f^n$ on $[0,1]$,
  $L_n=\int_0^1\sqrt{1+y_n'(x)^2}\dd x$. The $O(\cdot)$ refers to $n\to\infty$, with a constant
  that depends on $f$.

Quoted verbatim (sec_g3_fn.tex 18-23, the setup paragraph and `g3:eq:phibound`):

  $a=f'(1)>0$ and $\mu=\frac1{n-1}$. Let $\varphi(s)=1-\frac s{\sqrt{1+s^2}}$ ... for $s>0$ ...
  \[0<\varphi(s)\le\min\Bigl(1,\frac1{2s^2}\Bigr),\]
  ... So $\varphi$ is a probability density on $(0,\infty)$. ...
  For $y\in C^1[0,1]$ nondecreasing, $m_y(s)=|\{x\in[0,1]:y'(x)>s\}|$.

Quoted verbatim (sec_g3_fn.tex 27-29, Lemma `g3:lem:layer`):

  Let $y\in C^1[0,1]$ be nondecreasing with $y(1)-y(0)=1$. Then
  $2-\int_0^1\sqrt{1+y'^2}\dd x=\int_0^\infty\varphi(s)\,m_y(s)\dd s$.

Quoted verbatim (sec_g3_fn.tex 35-37, Lemma `g3:lem:moments`):

  $\E[\log S]=\log2-1$, and $\E[\log^2S]<\infty$.

Quoted verbatim (sec_g3_fn.tex 55-66, Lemma `g3:lem:shape`):

  There is $n_0$ such that for $n\ge n_0$, with $\sigma_n=nF_1\theta_0^{\,n-1}$:
  (a) $y_n'\le\sigma_n<na$ on $[0,1-\eta]$, and $y_n'$ is continuous and strictly increasing on
  $[1-\eta,1]$ with $y_n'(1)=na$;
  (b) $m_{y_n}(s)=0$ for $s\ge na$; $0\le m_{y_n}(s)\le1$ for all $s$; and for $\sigma_n\le s<na$,
  $m_{y_n}(s)=\tau_n(s)$, the unique $\tau\in(0,\eta]$ with $y_n'(1-\tau)=s$, equivalently with
  $A(\tau)=\mu\bigl(\lambda+B(\tau)\bigr)$, $\lambda=\log\frac ns$.
  (Here $\eta\in(0,1)$ is chosen with $f\in C^2[1-\eta,1]$, $f'\ge\frac a2$ on $[1-\eta,1]$,
  $f(1-\eta)\ge\frac12$; $\theta_0=f(1-\eta)<1$ and $F_1=\max_{[0,1]}f'$.)

Quoted verbatim (sec_g3_fn.tex 70-75, Lemma `g3:lem:window`):

  There are constants $C_1,n_0$ (depending on $f$) such that for $n\ge n_0$ and
  $\sigma_n\le s<na$, with $\lambda=\log\frac ns$,
  \[\Bigl|\tau_n(s)-\frac{\mu\,(\lambda+\log a)}a\Bigr|\le C_1\mu^2(1+|\lambda|)^2 .\]

## Formalization choices

* `n` is a natural number. The paper writes `L_n` as a sequence and lets `n → ∞` ("the graph of
  $x^n$", "$n=10^3,\dots,10^6$"); the supplement only says `n > 1`. The proof works verbatim for
  real `n > 1`, but the sequence reading is the one the symbol `L_n` and Ahmed's problem use, so
  `f x ^ n` is the natural-number power (`Monoid.npow`), no `rpow` is involved. Statements that
  are "for `n ≥ n_0`" are `∃ N : ℕ, ∀ n ≥ N`.
* `f : ℝ → ℝ` is arbitrary outside `[0,1]` and nothing depends on those values: every derivative is
  `derivWithin _ (Icc 0 1)`, which on `[0,1]` (a set of unique differentiability) is the two-sided
  derivative at interior points and the one-sided derivative at `0` and `1`.
* "`f ∈ C^1[0,1]`" is `ContDiffOn ℝ 1 f (Icc 0 1)` (derivative within `[0,1]` exists everywhere on
  `[0,1]` and is continuous on `[0,1]`, one-sided at the endpoints).
* "`f'(1)`" is `derivWithin f (Icc 0 1) 1`, the left derivative at `1`; "`f'(1) > 0`" is
  `0 < derivWithin f (Icc 0 1) 1`.
* "`f ∈ C^2` near `x = 1`" is `∃ η > 0, ContDiffOn ℝ 2 f (Icc (1 - η) 1)`, exactly the supplement's
  use of it ("$f\in C^2[1-\eta,1]$"); since `f` lives on `[0,1]`, "near 1" is a one-sided
  neighbourhood.
* "nondecreasing" is `MonotoneOn f (Icc 0 1)`.
* Length. `L f n = ∫ x in 0..1, √(1 + (y_n'(x))²)` with `y_n' = derivWithin (fun x ↦ f x ^ n)
  (Icc 0 1)`, as the supplement defines it. This is the arc length of the graph: for a `C^1`
  function `y` on `[0,1]` the curve `γ(x) = (x, y(x))` is `C^1`, and the length of a `C^1` curve
  (the supremum over partitions of the polygonal lengths) equals `∫ ‖γ'‖ = ∫ √(1 + y'²)`. Here
  `y_n = f^n` is `C^1` on `[0,1]` (product of `C^1` functions), so the integrand is continuous on
  `[0,1]`, hence integrable: the value is never the "non-integrable ⇒ 0" convention.
* `m_y(s)` is `(volume {x ∈ [0,1] | s < y'(x)}).toReal`; the set has finite measure (`≤ 1`), so
  `toReal` loses nothing.
* All integrals over `(0,∞)` are stated together with the integrability of their integrand, so no
  statement is true because of the "non-integrable ⇒ 0" convention. `Real.log` is only applied to
  positive arguments in the statements (`s > 0`, `n ≥ 2`, `a > 0`).
* The `O(log² n / n²)` is stated explicitly: `∃ C N, ∀ n ≥ N, |error| ≤ C * (log n)^2 / n^2`.

## Statements

* stmt_1: the bound `g3:eq:phibound`, `0 < φ(s) ≤ min(1, 1/(2s²))` for `s > 0`.
* stmt_2: Lemma `g3:lem:moments` (with "`φ` is a probability density"): `φ`, `φ·log`, `φ·log²`
  are integrable on `(0,∞)`, `∫ φ = 1` and `∫ φ(s) log s ds = log 2 - 1`.
* stmt_3: Lemma `g3:lem:layer` (layer cake): for `y ∈ C^1[0,1]` nondecreasing with
  `y(1) - y(0) = 1`, `s ↦ φ(s) m_y(s)` is integrable on `(0,∞)` and
  `2 - ∫_0^1 √(1+y'²) = ∫_0^∞ φ(s) m_y(s) ds`.
* stmt_4: Lemma `g3:lem:shape` for `y_n = f^n`, with the constants `η, θ_0, F_1` of the
  supplement existentially quantified (with the properties the supplement gives them).
* stmt_5: Lemma `g3:lem:window` (combined with the localisation of Lemma `g3:lem:shape`, which
  fixes `σ_n`): the window `m_{y_n}(s)` is `μ (λ + log a)/a` up to `C_1 μ² (1+|λ|)²` for
  `σ_n ≤ s < na`, with `σ_n = n F_1 θ_0^{n-1}`, `0 < θ_0 < 1`.
* stmt_6: Theorem `thm:fn` in full.

## NOT COVERED

* The higher-order expansion `2-L_n=\sum_k\varepsilon^k\mathcal J_k(\Xi)+O(\cdots)`, the cumulants
  `k_2,k_3,k_4`, the `2/e` slope law beyond `E log S = log 2 - 1`, and the higher-contact regime
  stated after the theorem in main.tex (not part of `thm:fn`).
* The identification "integral formula = arc length (supremum of polygonal lengths)" is explained
  above, not stated as a separate Lean statement.
-/
-- BEGIN STATEMENTS
namespace W01

open MeasureTheory Set Filter

/-- The density `φ(s) = 1 - s/√(1+s²)` of the random slope `S` (sec_g3_fn.tex, lines 18-19). -/
noncomputable def φ (s : ℝ) : ℝ := 1 - s / Real.sqrt (1 + s ^ 2)

/-- `m_y(s) = |{x ∈ [0,1] : y'(x) > s}|`, for a given derivative function `y'`. -/
noncomputable def levelMeasure (y' : ℝ → ℝ) (s : ℝ) : ℝ :=
  (volume {x : ℝ | x ∈ Icc (0 : ℝ) 1 ∧ s < y' x}).toReal

/-- The derivative on `[0,1]` (one-sided at the endpoints) of a function given on `[0,1]`. -/
noncomputable def D (y : ℝ → ℝ) (x : ℝ) : ℝ := derivWithin y (Icc (0 : ℝ) 1) x

/-- `y_n' = (f^n)'`, the derivative on `[0,1]` of `y_n = f^n` (`n` a natural number). -/
noncomputable def yd (f : ℝ → ℝ) (n : ℕ) (x : ℝ) : ℝ := D (fun t => f t ^ n) x

/-- The length `L_n = ∫_0^1 √(1 + y_n'(x)²) dx` of the graph of `y_n = f^n` on `[0,1]`. -/
noncomputable def L (f : ℝ → ℝ) (n : ℕ) : ℝ :=
  ∫ x in (0 : ℝ)..1, Real.sqrt (1 + (yd f n x) ^ 2)

/-- The hypotheses of `thm:fn`: `f ∈ C^1[0,1]` nondecreasing, `f(0)=0`, `f(1)=1`, `f'(1)>0`
(left derivative), and `f ∈ C^2[1-η,1]` for some `η > 0`. -/
def Hyp (f : ℝ → ℝ) : Prop :=
  ContDiffOn ℝ 1 f (Icc 0 1) ∧ MonotoneOn f (Icc 0 1) ∧ f 0 = 0 ∧ f 1 = 1 ∧
    0 < D f 1 ∧ ∃ η > 0, ContDiffOn ℝ 2 f (Icc (1 - η) 1)

/-- `g3:eq:phibound`: `0 < φ(s) ≤ min(1, 1/(2s²))` for `s > 0`. -/
def stmt_1 : Prop :=
  ∀ s : ℝ, 0 < s → 0 < φ s ∧ φ s ≤ min 1 (1 / (2 * s ^ 2))

/-- Lemma `g3:lem:moments` and "`φ` is a probability density on `(0,∞)`":
`∫ φ = 1`, `E[log S] = log 2 - 1`, `E[log² S] < ∞`. -/
def stmt_2 : Prop :=
  IntegrableOn φ (Ioi 0) ∧ ∫ s in Ioi (0 : ℝ), φ s = 1 ∧
    IntegrableOn (fun s => φ s * Real.log s) (Ioi 0) ∧
    IntegrableOn (fun s => φ s * Real.log s ^ 2) (Ioi 0) ∧
    ∫ s in Ioi (0 : ℝ), φ s * Real.log s = Real.log 2 - 1

/-- Lemma `g3:lem:layer` (layer cake): for `y ∈ C^1[0,1]` nondecreasing with `y(1)-y(0)=1`,
`2 - ∫_0^1 √(1+y'²) dx = ∫_0^∞ φ(s) m_y(s) ds` (the right side is an honest integral). -/
def stmt_3 : Prop :=
  ∀ y : ℝ → ℝ, ContDiffOn ℝ 1 y (Icc 0 1) → MonotoneOn y (Icc 0 1) → y 1 - y 0 = 1 →
    IntegrableOn (fun s => φ s * levelMeasure (D y) s) (Ioi 0) ∧
    2 - ∫ x in (0 : ℝ)..1, Real.sqrt (1 + (D y x) ^ 2) =
      ∫ s in Ioi (0 : ℝ), φ s * levelMeasure (D y) s

/-- Lemma `g3:lem:shape`, for `y_n = f^n` and `a = f'(1)`: there are `η ∈ (0,1)`,
`θ_0 = f(1-η) ∈ [1/2,1)`, `F_1 ≥ a` and `n_0` such that for `n ≥ n_0`, with
`σ_n = n F_1 θ_0^{n-1}`: (a) `y_n' ≤ σ_n < na` on `[0,1-η]`, `y_n'` is continuous and strictly
increasing on `[1-η,1]`, `y_n'(1) = na`; (b) `m_{y_n}(s) = 0` for `s ≥ na`, `0 ≤ m_{y_n} ≤ 1`,
and for `σ_n ≤ s < na`, `m_{y_n}(s) = τ` for the unique `τ ∈ (0,η]` with `y_n'(1-τ) = s`. -/
def stmt_4 : Prop :=
  ∀ f : ℝ → ℝ, Hyp f →
    ∃ η θ₀ F₁ : ℝ, ∃ n₀ : ℕ, 0 < η ∧ η < 1 ∧ θ₀ = f (1 - η) ∧ 1 / 2 ≤ θ₀ ∧ θ₀ < 1 ∧
      D f 1 ≤ F₁ ∧
      ∀ n : ℕ, n₀ ≤ n →
        let σ : ℝ := n * F₁ * θ₀ ^ (n - 1)
        let a : ℝ := D f 1
        (∀ x ∈ Icc (0 : ℝ) (1 - η), yd f n x ≤ σ) ∧ σ < n * a ∧
        ContinuousOn (yd f n) (Icc (1 - η) 1) ∧ StrictMonoOn (yd f n) (Icc (1 - η) 1) ∧
        yd f n 1 = n * a ∧
        (∀ s : ℝ, n * a ≤ s → levelMeasure (yd f n) s = 0) ∧
        (∀ s : ℝ, 0 ≤ levelMeasure (yd f n) s ∧ levelMeasure (yd f n) s ≤ 1) ∧
        (∀ s : ℝ, σ ≤ s → s < n * a →
          ∃! τ : ℝ, τ ∈ Ioc 0 η ∧ yd f n (1 - τ) = s ∧ levelMeasure (yd f n) s = τ)

/-- Lemma `g3:lem:window` (with the `σ_n` of Lemma `g3:lem:shape`): with `a = f'(1)`,
`μ = 1/(n-1)`, `λ = log(n/s)`, for `n ≥ n_0` and `σ_n ≤ s < na`,
`|m_{y_n}(s) - μ(λ + log a)/a| ≤ C_1 μ² (1+|λ|)²`. Here `σ_n = n F_1 θ_0^{n-1}` with
`0 < θ_0 < 1` and `F_1 > 0`, so `σ_n > 0` is exponentially small. -/
def stmt_5 : Prop :=
  ∀ f : ℝ → ℝ, Hyp f →
    ∃ θ₀ F₁ C₁ : ℝ, ∃ n₀ : ℕ, 0 < θ₀ ∧ θ₀ < 1 ∧ 0 < F₁ ∧
      ∀ n : ℕ, n₀ ≤ n → ∀ s : ℝ,
        let σ : ℝ := n * F₁ * θ₀ ^ (n - 1)
        let a : ℝ := D f 1
        let μ : ℝ := 1 / ((n : ℝ) - 1)
        let lam : ℝ := Real.log (n / s)
        σ ≤ s → s < n * a →
          |levelMeasure (yd f n) s - μ * (lam + Real.log a) / a| ≤ C₁ * μ ^ 2 * (1 + |lam|) ^ 2

/-- Theorem `thm:fn`: `2 - L_n = (log n + log f'(1) + 1 - log 2)/(f'(1)(n-1)) + O(log² n / n²)`,
with the `O` made explicit. -/
def stmt_6 : Prop :=
  ∀ f : ℝ → ℝ, Hyp f →
    ∃ C : ℝ, ∃ N : ℕ, ∀ n : ℕ, N ≤ n →
      |(2 - L f n) - (Real.log n + Real.log (D f 1) + 1 - Real.log 2) / (D f 1 * ((n : ℝ) - 1))|
        ≤ C * Real.log n ^ 2 / (n : ℝ) ^ 2

end W01
-- END STATEMENTS

namespace W01
open MeasureTheory Set Filter

lemma phi_eq (s : ℝ) (hs : 0 ≤ s) :
    φ s = 1 / (Real.sqrt (1 + s ^ 2) * (Real.sqrt (1 + s ^ 2) + s)) := by
  have hr : 0 < Real.sqrt (1 + s ^ 2) := Real.sqrt_pos.2 (by positivity)
  have hr2 : Real.sqrt (1 + s ^ 2) ^ 2 = 1 + s ^ 2 := Real.sq_sqrt (by positivity)
  unfold φ
  field_simp
  nlinarith

lemma phi_pos (s : ℝ) (hs : 0 ≤ s) : 0 < φ s := by
  rw [phi_eq s hs]; have : 0 < Real.sqrt (1 + s ^ 2) := Real.sqrt_pos.2 (by positivity)
  positivity

lemma phi_le_one (s : ℝ) (hs : 0 ≤ s) : φ s ≤ 1 := by
  unfold φ; have : 0 ≤ s / Real.sqrt (1 + s ^ 2) := by positivity
  linarith

lemma phi_le_inv (s : ℝ) (hs : 0 < s) : φ s ≤ 1 / (2 * s ^ 2) := by
  rw [phi_eq s hs.le]
  have hr : s ≤ Real.sqrt (1 + s ^ 2) := Real.le_sqrt_of_sq_le (by linarith)
  apply one_div_le_one_div_of_le (by positivity)
  nlinarith

end W01

namespace W01
open MeasureTheory Set Filter Topology

lemma hasDerivAt_sqrt1 (s : ℝ) :
    HasDerivAt (fun t => Real.sqrt (1 + t ^ 2)) (s / Real.sqrt (1 + s ^ 2)) s := by
  have h1 : HasDerivAt (fun t : ℝ => 1 + t ^ 2) (2 * s) s := by
    simpa using (hasDerivAt_pow 2 s).const_add 1
  have := h1.sqrt (by positivity)
  convert this using 1
  field_simp

lemma hasDerivAt_G0 (s : ℝ) :
    HasDerivAt (fun t => t - Real.sqrt (1 + t ^ 2)) (φ s) s := by
  have := (hasDerivAt_id' s).sub (hasDerivAt_sqrt1 s)
  convert this using 1
  rfl

lemma tendsto_G0 : Tendsto (fun t : ℝ => t - Real.sqrt (1 + t ^ 2)) atTop (𝓝 0) := by
  have h : ∀ t : ℝ, 0 ≤ t → t - Real.sqrt (1 + t ^ 2) = -(t + Real.sqrt (1 + t ^ 2))⁻¹ := by
    intro t ht
    have hr : 0 < Real.sqrt (1 + t ^ 2) := Real.sqrt_pos.2 (by positivity)
    have hr2 : Real.sqrt (1 + t ^ 2) ^ 2 = 1 + t ^ 2 := Real.sq_sqrt (by positivity)
    field_simp
    nlinarith
  have h2 : Tendsto (fun t : ℝ => t + Real.sqrt (1 + t ^ 2)) atTop atTop :=
    tendsto_atTop_mono (fun t => le_add_of_nonneg_right (Real.sqrt_nonneg (1 + t ^ 2))) tendsto_id
  have := (h2.inv_tendsto_atTop).neg
  simp only [neg_zero] at this
  refine this.congr' ?_
  filter_upwards [eventually_ge_atTop 0] with t ht
  exact (h t ht).symm

lemma phi_integrable : IntegrableOn φ (Ioi 0) :=
  integrableOn_Ioi_deriv_of_nonneg' (fun x _ => hasDerivAt_G0 x)
    (fun x hx => (phi_pos x (le_of_lt (mem_Ioi.1 hx))).le) tendsto_G0

lemma phi_integral : ∫ s in Ioi (0 : ℝ), φ s = 1 := by
  rw [integral_Ioi_of_hasDerivAt_of_nonneg' (fun x _ => hasDerivAt_G0 x)
    (fun x hx => (phi_pos x (le_of_lt (mem_Ioi.1 hx))).le) tendsto_G0]
  simp

end W01

namespace W01
open MeasureTheory Set Filter Topology

lemma rpow_quarter_facts (s : ℝ) (hs : 0 < s) :
    let u := s ^ (1 / 4 : ℝ)
    0 < u ∧ s = u ^ 4 ∧ s ^ (-1 / 2 : ℝ) = (u ^ 2)⁻¹ ∧ s ^ (-3 / 2 : ℝ) = (u ^ 6)⁻¹ := by
  intro u
  have hu : 0 < u := Real.rpow_pos_of_pos hs _
  have key : ∀ k : ℕ, u ^ k = s ^ ((k : ℝ) / 4) := by
    intro k
    rw [← Real.rpow_natCast, ← Real.rpow_mul hs.le]; ring_nf
  refine ⟨hu, ?_, ?_, ?_⟩
  · rw [key 4]; norm_num
  · rw [key 2, ← Real.rpow_neg hs.le]; norm_num
  · rw [key 6, ← Real.rpow_neg hs.le]; norm_num

lemma log_sq_le_small (s : ℝ) (hs : 0 < s) (hs1 : s ≤ 1) :
    Real.log s ^ 2 ≤ 16 * s ^ (-1 / 2 : ℝ) := by
  obtain ⟨hu, hsu, h2, -⟩ := rpow_quarter_facts s hs
  have h := Real.abs_log_mul_self_rpow_lt s (1 / 4) hs hs1 (by norm_num)
  rw [h2]
  set u := s ^ (1 / 4 : ℝ)
  have h' : (Real.log s * u) ^ 2 ≤ 16 := by
    have := sq_le_sq' (abs_lt.1 h).1.le (abs_lt.1 h).2.le
    norm_num at this ⊢; linarith
  rw [le_mul_inv_iff₀ (by positivity)]
  nlinarith

lemma log_sq_le_large (s : ℝ) (hs : 1 ≤ s) :
    Real.log s ^ 2 ≤ 16 * (s ^ (1 / 4 : ℝ)) ^ 2 := by
  have h := Real.log_le_rpow_div (x := s) (by linarith) (ε := 1 / 4) (by norm_num)
  have h0 := Real.log_nonneg hs
  have : Real.log s ≤ 4 * s ^ (1 / 4 : ℝ) := by
    calc Real.log s ≤ s ^ (1 / 4 : ℝ) / (1 / 4) := h
      _ = 4 * s ^ (1 / 4 : ℝ) := by ring
  nlinarith

lemma phi_continuous : Continuous φ := by
  show Continuous (fun s : ℝ => 1 - s / Real.sqrt (1 + s ^ 2))
  apply continuous_const.sub
  apply continuous_id.div (by fun_prop)
  intro x; positivity

lemma phi_log_cont : ContinuousOn (fun s => φ s * Real.log s ^ 2) (Ioi 0) := by
  apply ContinuousOn.mul
  · exact phi_continuous.continuousOn
  · exact (Real.continuousOn_log.mono (fun x hx => ne_of_gt hx)).pow 2

lemma phi_log2_integrable : IntegrableOn (fun s => φ s * Real.log s ^ 2) (Ioi 0) := by
  rw [← Ioc_union_Ioi_eq_Ioi (zero_le_one' ℝ)]
  apply IntegrableOn.union
  · have hint : IntegrableOn (fun s : ℝ => 16 * s ^ (-1 / 2 : ℝ)) (Ioc 0 1) := by
      apply Integrable.const_mul
      have := intervalIntegral.intervalIntegrable_rpow' (a := 0) (b := 1) (r := -1 / 2) (by norm_num)
      rwa [intervalIntegrable_iff_integrableOn_Ioc_of_le zero_le_one] at this
    refine hint.mono' ((phi_log_cont.mono Ioc_subset_Ioi_self).aestronglyMeasurable
      measurableSet_Ioc) ?_
    refine (ae_restrict_iff' measurableSet_Ioc).2 (Eventually.of_forall fun s hs => ?_)
    have h1 := phi_pos s hs.1.le
    have h2 := phi_le_one s hs.1.le
    rw [Real.norm_eq_abs, abs_of_nonneg (by positivity)]
    have := log_sq_le_small s hs.1 hs.2
    nlinarith [sq_nonneg (Real.log s)]
  · have hint : IntegrableOn (fun s : ℝ => 8 * s ^ (-3 / 2 : ℝ)) (Ioi 1) :=
      (integrableOn_Ioi_rpow_of_lt (by norm_num) zero_lt_one).const_mul 8
    refine hint.mono' ((phi_log_cont.mono (Ioi_subset_Ioi zero_le_one)).aestronglyMeasurable
      measurableSet_Ioi) ?_
    refine (ae_restrict_iff' measurableSet_Ioi).2 (Eventually.of_forall fun s hs => ?_)
    have hs' : (0 : ℝ) < s := lt_trans zero_lt_one hs
    have h1 := phi_pos s hs'.le
    have h2 := phi_le_inv s hs'
    rw [Real.norm_eq_abs, abs_of_nonneg (by positivity)]
    have := log_sq_le_large s hs.le
    obtain ⟨hu, hsu, -, h3⟩ := rpow_quarter_facts s hs'
    rw [h3]
    set u := s ^ (1 / 4 : ℝ)
    calc φ s * Real.log s ^ 2 ≤ 1 / (2 * s ^ 2) * (16 * u ^ 2) :=
          mul_le_mul h2 this (sq_nonneg _) (by positivity)
      _ = 8 * (u ^ 6)⁻¹ := by rw [hsu]; field_simp; ring

lemma phi_log_integrable : IntegrableOn (fun s => φ s * Real.log s) (Ioi 0) := by
  refine (phi_integrable.add phi_log2_integrable).mono' ?_ ?_
  · apply ContinuousOn.aestronglyMeasurable _ measurableSet_Ioi
    apply ContinuousOn.mul
    · exact phi_continuous.continuousOn
    · exact Real.continuousOn_log.mono (fun x hx => ne_of_gt hx)
  · refine (ae_restrict_iff' measurableSet_Ioi).2 (Eventually.of_forall fun s hs => ?_)
    have h1 := phi_pos s (le_of_lt hs)
    rw [Real.norm_eq_abs, abs_mul, abs_of_pos h1]
    have : |Real.log s| ≤ 1 + Real.log s ^ 2 := by
      rcases abs_cases (Real.log s) with ⟨h, -⟩ | ⟨h, -⟩ <;> rw [h] <;>
        nlinarith [sq_nonneg (Real.log s - 1), sq_nonneg (Real.log s + 1)]
    show _ ≤ φ s + φ s * Real.log s ^ 2
    nlinarith

end W01

namespace W01
open MeasureTheory Set Filter Topology

/-- antiderivative of `φ log`. -/
noncomputable def Glog (s : ℝ) : ℝ :=
  (s - Real.sqrt (1 + s ^ 2) + 1) * Real.log s - s + Real.sqrt (1 + s ^ 2)
    - Real.log (1 + Real.sqrt (1 + s ^ 2))

lemma hasDerivAt_Glog (s : ℝ) (hs : 0 < s) : HasDerivAt Glog (φ s * Real.log s) s := by
  set r := Real.sqrt (1 + s ^ 2) with hr_def
  have hr : 1 ≤ r := Real.le_sqrt_of_sq_le (by nlinarith)
  have hr2 : r ^ 2 = 1 + s ^ 2 := Real.sq_sqrt (by positivity)
  have hR := hasDerivAt_sqrt1 s
  have h1 := (((hasDerivAt_id' s).sub hR).add_const 1).mul (Real.hasDerivAt_log hs.ne')
  have h2 := ((h1.sub (hasDerivAt_id' s)).add hR).sub ((hR.const_add 1).log (by positivity))
  convert h2 using 1
  case e'_8 => funext x; simp only [Glog, Pi.sub_apply, Pi.add_apply, Pi.mul_apply]
  show φ s * Real.log s = _
  simp only [Pi.sub_apply]
  unfold φ
  rw [← hr_def]
  field_simp
  linear_combination r * hr2

lemma Glog_zero : Glog 0 = 1 - Real.log 2 := by
  unfold Glog; norm_num

end W01

namespace W01
open MeasureTheory Set Filter Topology

lemma Glog_cont0 : ContinuousWithinAt Glog (Ici 0) 0 := by
  have hA : ContinuousWithinAt
      (fun s : ℝ => (s - Real.sqrt (1 + s ^ 2) + 1) * Real.log s) (Ici 0) 0 := by
    show Tendsto _ _ _
    simp only [Real.log_zero, mul_zero]
    have hB : Tendsto (fun s : ℝ => s * Real.log s) (𝓝[Ici 0] 0) (𝓝 0) := by
      have := (Real.continuous_mul_log.tendsto 0).mono_left
        (nhdsWithin_le_nhds (s := Ici (0 : ℝ)))
      simpa using this
    refine squeeze_zero_norm' ?_ hB.norm
    filter_upwards [self_mem_nhdsWithin] with s hs
    have hs0 : (0 : ℝ) ≤ s := hs
    have hr : 1 ≤ Real.sqrt (1 + s ^ 2) := Real.le_sqrt_of_sq_le (by nlinarith)
    have hr' : Real.sqrt (1 + s ^ 2) ≤ s + 1 := Real.sqrt_le_iff.2 ⟨by linarith, by nlinarith⟩
    simp only [norm_mul, Real.norm_eq_abs]
    rw [abs_of_nonneg (by linarith : (0:ℝ) ≤ s - Real.sqrt (1 + s ^ 2) + 1)]
    exact mul_le_mul_of_nonneg_right (le_trans (by linarith) (Real.le_norm_self s))
      (abs_nonneg _)
  have hC : Continuous (fun s : ℝ => - s + Real.sqrt (1 + s ^ 2)
      - Real.log (1 + Real.sqrt (1 + s ^ 2))) := by
    exact (continuous_neg.add (by fun_prop : Continuous fun s : ℝ => Real.sqrt (1 + s ^ 2))).sub
      ((by fun_prop : Continuous fun s : ℝ => 1 + Real.sqrt (1 + s ^ 2)).log
        (fun x => by positivity))
  have := hA.add hC.continuousWithinAt
  convert this using 1
  funext s; simp only [Glog, Pi.add_apply]; ring

lemma Glog_tendsto : Tendsto Glog atTop (𝓝 0) := by
  -- three pieces
  have h1 : Tendsto (fun s : ℝ => (s - Real.sqrt (1 + s ^ 2)) * Real.log s) atTop (𝓝 0) := by
    have hL : Tendsto (fun s : ℝ => Real.log s ^ 1 / (1 * s + 0)) atTop (𝓝 0) :=
      Real.tendsto_pow_log_div_mul_add_atTop 1 0 1 one_ne_zero
    refine squeeze_zero_norm' ?_ hL
    filter_upwards [eventually_ge_atTop 1] with s hs
    have hs0 : 0 < s := by linarith
    have hr : s ≤ Real.sqrt (1 + s ^ 2) := Real.le_sqrt_of_sq_le (by nlinarith)
    have hr2 : Real.sqrt (1 + s ^ 2) ^ 2 = 1 + s ^ 2 := Real.sq_sqrt (by positivity)
    have hl := Real.log_nonneg hs
    have hd : Real.sqrt (1 + s ^ 2) - s ≤ 1 / s := by
      rw [le_div_iff₀ hs0]; nlinarith
    rw [norm_mul, Real.norm_eq_abs, Real.norm_eq_abs, abs_of_nonpos (by linarith),
      abs_of_nonneg hl]
    simp only [pow_one, one_mul, add_zero]
    calc -(s - Real.sqrt (1 + s ^ 2)) * Real.log s ≤ 1 / s * Real.log s :=
          mul_le_mul_of_nonneg_right (by linarith) hl
      _ = Real.log s / s := by ring
  have h2 : Tendsto (fun s : ℝ => Real.log s - Real.log (1 + Real.sqrt (1 + s ^ 2))) atTop
      (𝓝 0) := by
    have hq : Tendsto (fun s : ℝ => s / (1 + Real.sqrt (1 + s ^ 2))) atTop (𝓝 1) := by
      have hlow : Tendsto (fun s : ℝ => 1 - 2 / s) atTop (𝓝 1) := by
        simpa using (tendsto_const_nhds (x := (1 : ℝ))).sub
          ((tendsto_const_nhds (x := (2 : ℝ))).div_atTop tendsto_id)
      refine tendsto_of_tendsto_of_tendsto_of_le_of_le' hlow tendsto_const_nhds ?_ ?_
      · filter_upwards [eventually_ge_atTop 1] with s hs
        have hs0 : 0 < s := by linarith
        have hr : s ≤ Real.sqrt (1 + s ^ 2) := Real.le_sqrt_of_sq_le (by nlinarith)
        have hr' : Real.sqrt (1 + s ^ 2) ≤ s + 1 :=
          Real.sqrt_le_iff.2 ⟨by linarith, by nlinarith⟩
        rw [le_div_iff₀ (by linarith)]
        have : 2 / s * (1 + Real.sqrt (1 + s ^ 2)) ≥ 2 := by
          rw [ge_iff_le, div_mul_eq_mul_div, le_div_iff₀ hs0]; nlinarith
        nlinarith
      · filter_upwards [eventually_ge_atTop 1] with s hs
        have hr : s ≤ Real.sqrt (1 + s ^ 2) := Real.le_sqrt_of_sq_le (by nlinarith)
        rw [div_le_one (by linarith [Real.sqrt_nonneg (1 + s ^ 2)])]; linarith
    have := (Real.continuousAt_log one_ne_zero).tendsto.comp hq
    rw [Real.log_one] at this
    refine this.congr' ?_
    filter_upwards [eventually_gt_atTop 0] with s hs
    simp only [Function.comp]
    rw [Real.log_div hs.ne' (by positivity)]
  have h3 := tendsto_G0.neg
  rw [neg_zero] at h3
  have := (h1.add h2).add h3
  rw [add_zero, add_zero] at this
  refine this.congr' (Eventually.of_forall fun s => ?_)
  simp only [Glog]; ring

lemma phi_log_integral : ∫ s in Ioi (0 : ℝ), φ s * Real.log s = Real.log 2 - 1 := by
  rw [integral_Ioi_of_hasDerivAt_of_tendsto Glog_cont0 (fun x hx => hasDerivAt_Glog x hx)
    phi_log_integrable Glog_tendsto, Glog_zero]
  ring

end W01


namespace W01
open MeasureTheory Set Filter Topology

lemma levelMeasure_le_one (g : ℝ → ℝ) (s : ℝ) : 0 ≤ levelMeasure g s ∧ levelMeasure g s ≤ 1 := by
  refine ⟨ENNReal.toReal_nonneg, ?_⟩
  unfold levelMeasure
  have h : volume {x : ℝ | x ∈ Icc (0 : ℝ) 1 ∧ s < g x} ≤ volume (Icc (0 : ℝ) 1) :=
    measure_mono (fun x hx => hx.1)
  rw [Real.volume_Icc] at h
  have := ENNReal.toReal_mono ENNReal.ofReal_ne_top h
  simpa using this

lemma levelMeasure_antitone (g : ℝ → ℝ) : Antitone (levelMeasure g) := by
  intro s t hst
  unfold levelMeasure
  apply ENNReal.toReal_mono
  · exact ne_top_of_le_ne_top (b := volume (Icc (0 : ℝ) 1)) (by simp)
      (measure_mono (fun x hx => hx.1))
  · exact measure_mono (fun x hx => ⟨hx.1, lt_of_le_of_lt hst hx.2⟩)

lemma phi_level_integrable (g : ℝ → ℝ) :
    IntegrableOn (fun s => φ s * levelMeasure g s) (Ioi 0) := by
  refine phi_integrable.mono' ?_ ?_
  · exact (phi_continuous.measurable.mul
      (levelMeasure_antitone g).measurable).aestronglyMeasurable
  · refine (ae_restrict_iff' measurableSet_Ioi).2 (Eventually.of_forall fun s hs => ?_)
    have h1 := phi_pos s (le_of_lt hs)
    obtain ⟨h2, h3⟩ := levelMeasure_le_one g s
    rw [Real.norm_eq_abs, abs_of_nonneg (by positivity)]
    nlinarith

lemma integral_phi_eq (t : ℝ) : ∫ s in (0 : ℝ)..t, φ s = 1 + t - Real.sqrt (1 + t ^ 2) := by
  rw [intervalIntegral.integral_eq_sub_of_hasDerivAt (fun x _ => hasDerivAt_G0 x)
    (phi_continuous.intervalIntegrable _ _)]
  norm_num; ring

lemma layer (y : ℝ → ℝ) (hy : ContDiffOn ℝ 1 y (Icc 0 1)) (hm : MonotoneOn y (Icc 0 1))
    (h1 : y 1 - y 0 = 1) :
    2 - ∫ x in (0 : ℝ)..1, Real.sqrt (1 + (D y x) ^ 2) =
      ∫ s in Ioi (0 : ℝ), φ s * levelMeasure (D y) s := by
  set g := D y with hg_def
  have hgc : ContinuousOn g (Icc 0 1) :=
    hy.continuousOn_derivWithin (uniqueDiffOn_Icc zero_lt_one) le_rfl
  have hg0 : ∀ x, 0 ≤ g x := fun x => hm.derivWithin_nonneg
  have hint1 : ∫ x in (0 : ℝ)..1, g x = 1 := by
    have := intervalIntegral.integral_derivWithin_Icc_of_contDiffOn_Icc (a := (0 : ℝ)) (b := 1)
      hy zero_le_one
    rw [h1] at this
    exact this
  have hsc : ContinuousOn (fun x => Real.sqrt (1 + g x ^ 2)) (Icc 0 1) :=
    (continuousOn_const.add (hgc.pow 2)).sqrt
  -- rewrite LHS as ∫ F
  set F : ℝ → ℝ := fun x => 1 + g x - Real.sqrt (1 + g x ^ 2) with hF
  have hF0 : ∀ x, 0 ≤ F x := by
    intro x
    have : Real.sqrt (1 + g x ^ 2) ≤ 1 + g x :=
      Real.sqrt_le_iff.2 ⟨by linarith [hg0 x], by nlinarith [hg0 x]⟩
    simp only [hF]; linarith
  have hFc : ContinuousOn F (Icc 0 1) := (continuousOn_const.add hgc).sub hsc
  have hL : 2 - ∫ x in (0 : ℝ)..1, Real.sqrt (1 + g x ^ 2) = ∫ x in Icc (0 : ℝ) 1, F x := by
    have hgi : IntervalIntegrable g volume 0 1 :=
      (hgc.mono (by rw [uIcc_of_le zero_le_one])).intervalIntegrable
    have hsi : IntervalIntegrable (fun x => Real.sqrt (1 + g x ^ 2)) volume 0 1 :=
      (hsc.mono (by rw [uIcc_of_le zero_le_one])).intervalIntegrable
    have : ∫ x in (0 : ℝ)..1, F x = 2 - ∫ x in (0 : ℝ)..1, Real.sqrt (1 + g x ^ 2) := by
      simp only [hF]
      rw [intervalIntegral.integral_sub ((intervalIntegrable_const).add hgi) hsi,
        intervalIntegral.integral_add intervalIntegrable_const hgi, hint1]
      norm_num
    rw [← this, intervalIntegral.integral_of_le zero_le_one, integral_Icc_eq_integral_Ioc]
  rw [hL]
  -- layer cake
  set μ : Measure ℝ := volume.restrict (Icc (0 : ℝ) 1)
  have hcake := lintegral_comp_eq_lintegral_meas_lt_mul μ (f := g) (g := φ)
    (Eventually.of_forall hg0) (hgc.aemeasurable measurableSet_Icc)
    (fun t _ => phi_continuous.intervalIntegrable _ _)
    ((ae_restrict_iff' measurableSet_Ioi).2
      (Eventually.of_forall fun t ht => (phi_pos t (le_of_lt ht)).le))
  have hFeq : ∀ x, F x = ∫ s in (0 : ℝ)..g x, φ s := by
    intro x; rw [integral_phi_eq]
  rw [integral_eq_lintegral_of_nonneg_ae (Eventually.of_forall hF0)
    (hFc.aestronglyMeasurable measurableSet_Icc)]
  rw [integral_eq_lintegral_of_nonneg_ae
    ((ae_restrict_iff' measurableSet_Ioi).2 (Eventually.of_forall fun t ht =>
      mul_nonneg (phi_pos t (le_of_lt (mem_Ioi.1 ht))).le (levelMeasure_le_one g t).1))
    (phi_level_integrable g).aestronglyMeasurable]
  congr 1
  simp_rw [hFeq]
  rw [hcake]
  apply setLIntegral_congr_fun measurableSet_Ioi
  intro t ht
  have hphi := (phi_pos t (le_of_lt (mem_Ioi.1 ht))).le
  show μ {a | t < g a} * ENNReal.ofReal (φ t) = ENNReal.ofReal (φ t * levelMeasure g t)
  rw [ENNReal.ofReal_mul hphi, mul_comm]
  congr 1
  unfold levelMeasure
  rw [ENNReal.ofReal_toReal (ne_top_of_le_ne_top (b := volume (Icc (0 : ℝ) 1)) (by simp)
      (measure_mono (s := {x : ℝ | x ∈ Icc (0 : ℝ) 1 ∧ t < g x}) (fun x hx => hx.1))),
    Measure.restrict_apply' measurableSet_Icc]
  congr 1
  ext x; simp [and_comm]

end W01

namespace W01
open MeasureTheory Set Filter Topology

/-- The localisation data of the supplement. -/
def Setup (f : ℝ → ℝ) (η θ₀ F₁ β : ℝ) (g2 : ℝ → ℝ) : Prop :=
      0 < D f 1 ∧ f 1 = 1 ∧ MonotoneOn f (Icc 0 1) ∧
      0 < η ∧ η < 1 ∧ θ₀ = f (1 - η) ∧ 1 / 2 ≤ θ₀ ∧ θ₀ < 1 ∧ D f 1 ≤ F₁ ∧ 0 ≤ β ∧
      (∀ x ∈ Icc (0:ℝ) 1, 0 ≤ f x ∧ f x ≤ 1) ∧
      (∀ x ∈ Icc (0:ℝ) 1, 0 ≤ D f x ∧ D f x ≤ F₁) ∧
      (∀ x ∈ Icc (0:ℝ) 1, HasDerivWithinAt f (D f x) (Icc 0 1) x) ∧
      ContinuousOn (D f) (Icc 0 1) ∧
      (∀ x ∈ Icc (1 - η) 1, D f 1 / 2 ≤ D f x) ∧
      (∀ x ∈ Icc (1 - η) 1, HasDerivWithinAt (D f) (g2 x) (Icc (1 - η) 1) x) ∧
      (∀ x ∈ Icc (1 - η) 1, |g2 x| ≤ β)

lemma setup (f : ℝ → ℝ) (hf : Hyp f) :
    ∃ η θ₀ F₁ β : ℝ, ∃ g2 : ℝ → ℝ,
      0 < D f 1 ∧ f 1 = 1 ∧ MonotoneOn f (Icc 0 1) ∧
      0 < η ∧ η < 1 ∧ θ₀ = f (1 - η) ∧ 1 / 2 ≤ θ₀ ∧ θ₀ < 1 ∧ D f 1 ≤ F₁ ∧ 0 ≤ β ∧
      (∀ x ∈ Icc (0:ℝ) 1, 0 ≤ f x ∧ f x ≤ 1) ∧
      (∀ x ∈ Icc (0:ℝ) 1, 0 ≤ D f x ∧ D f x ≤ F₁) ∧
      (∀ x ∈ Icc (0:ℝ) 1, HasDerivWithinAt f (D f x) (Icc 0 1) x) ∧
      ContinuousOn (D f) (Icc 0 1) ∧
      (∀ x ∈ Icc (1 - η) 1, D f 1 / 2 ≤ D f x) ∧
      (∀ x ∈ Icc (1 - η) 1, HasDerivWithinAt (D f) (g2 x) (Icc (1 - η) 1) x) ∧
      (∀ x ∈ Icc (1 - η) 1, |g2 x| ≤ β) := by
  obtain ⟨hC1, hmono, h0, h1, ha, η₀, hη₀, hC2⟩ := hf
  set a := D f 1 with ha_def
  have hU : UniqueDiffOn ℝ (Icc (0:ℝ) 1) := uniqueDiffOn_Icc zero_lt_one
  have hDc : ContinuousOn (D f) (Icc 0 1) := hC1.continuousOn_derivWithin hU le_rfl
  have hfc : ContinuousOn f (Icc 0 1) := hC1.continuousOn
  have hderiv : ∀ x ∈ Icc (0:ℝ) 1, HasDerivWithinAt f (D f x) (Icc 0 1) x :=
    fun x hx => ((hC1.differentiableOn one_ne_zero) x hx).hasDerivWithinAt
  have hf01 : ∀ x ∈ Icc (0:ℝ) 1, 0 ≤ f x ∧ f x ≤ 1 := fun x hx =>
    ⟨h0 ▸ hmono ⟨le_refl 0, zero_le_one⟩ hx hx.1, h1 ▸ hmono hx ⟨zero_le_one, le_refl 1⟩ hx.2⟩
  have h1mem : (1:ℝ) ∈ Icc (0:ℝ) 1 := ⟨zero_le_one, le_rfl⟩
  obtain ⟨δ1, hδ1, hδ1'⟩ := Metric.continuousWithinAt_iff.1 (hDc 1 h1mem) (a / 2) (by positivity)
  obtain ⟨δ2, hδ2, hδ2'⟩ := Metric.continuousWithinAt_iff.1 (hfc 1 h1mem) (1 / 2) (by norm_num)
  set η := min (η₀ / 2) (min (δ1 / 2) (min (δ2 / 2) (1 / 2))) with hη
  have hηpos : 0 < η := lt_min (by positivity) (lt_min (by positivity) (lt_min (by positivity)
    (by norm_num)))
  have hη1 : η ≤ η₀ / 2 := min_le_left _ _
  have hη2 : η ≤ δ1 / 2 := le_trans (min_le_right _ _) (min_le_left _ _)
  have hη3 : η ≤ δ2 / 2 :=
    le_trans (min_le_right _ _) (le_trans (min_le_right _ _) (min_le_left _ _))
  have hη4 : η ≤ 1 / 2 :=
    le_trans (min_le_right _ _) (le_trans (min_le_right _ _) (min_le_right _ _))
  have hsub : Icc (1 - η) 1 ⊆ Icc (0:ℝ) 1 := Icc_subset_Icc (by linarith) le_rfl
  have hlow : ∀ x ∈ Icc (1 - η) 1, a / 2 ≤ D f x := by
    intro x hx
    have hd : dist x 1 < δ1 := by
      rw [Real.dist_eq, abs_of_nonpos (by linarith [hx.2])]; linarith [hx.1]
    have := hδ1' (hsub hx) hd
    rw [Real.dist_eq] at this
    linarith [(abs_lt.1 this).1]
  have hθhalf : 1 / 2 ≤ f (1 - η) := by
    have hd : dist (1 - η) 1 < δ2 := by
      rw [Real.dist_eq, abs_of_nonpos (by linarith)]; linarith
    have := hδ2' (hsub ⟨le_rfl, by linarith⟩) hd
    rw [Real.dist_eq, h1] at this
    linarith [(abs_lt.1 this).1]
  have hθlt : f (1 - η) < 1 := by
    have hsm : StrictMonoOn f (Icc (1 - η) 1) := by
      apply strictMonoOn_of_hasDerivWithinAt_pos (convex_Icc _ _) (hfc.mono hsub)
        (f' := D f)
      · intro x hx
        rw [interior_Icc] at hx ⊢
        exact (hderiv x (hsub (Ioo_subset_Icc_self hx))).mono
          (Ioo_subset_Icc_self.trans hsub)
      · intro x hx
        rw [interior_Icc] at hx
        linarith [hlow x (Ioo_subset_Icc_self hx)]
    have := hsm ⟨le_rfl, by linarith⟩ ⟨by linarith, le_rfl⟩ (by linarith : 1 - η < 1)
    rwa [h1] at this
  obtain ⟨C, hC⟩ := isCompact_Icc.exists_bound_of_continuousOn hDc
  -- second derivative
  set I' := Icc (1 - η₀) (1:ℝ)
  have hU' : UniqueDiffOn ℝ I' := uniqueDiffOn_Icc (by linarith)
  set g1 := derivWithin f I'
  have hg1 : ContDiffOn ℝ 1 g1 I' := hC2.derivWithin hU' (by norm_num)
  set g2 := derivWithin g1 I'
  have hg2c : ContinuousOn g2 I' := hg1.continuousOn_derivWithin hU' le_rfl
  obtain ⟨C2, hC2b⟩ := isCompact_Icc.exists_bound_of_continuousOn hg2c
  have hsub' : Icc (1 - η) 1 ⊆ I' := Icc_subset_Icc (by linarith) le_rfl
  have hDg1 : ∀ x ∈ Icc (1 - η) 1, D f x = g1 x := by
    intro x hx
    have h2 : HasDerivWithinAt f (g1 x) I' x :=
      ((hC2.differentiableOn (by norm_num)) x (hsub' hx)).hasDerivWithinAt
    have h3 : HasDerivWithinAt f (g1 x) (Icc 0 1) x := by
      apply h2.mono_of_mem_nhdsWithin
      refine mem_nhdsWithin_iff_exists_mem_nhds_inter.2
        ⟨Ioi (1 - η₀), Ioi_mem_nhds (by linarith [hx.1]), ?_⟩
      rintro t ⟨ht1, ht2⟩
      exact ⟨le_of_lt ht1, ht2.2⟩
    exact h3.derivWithin (hU x (hsub hx))
  refine ⟨η, f (1 - η), max C a, max C2 0, g2, ha, h1, hmono, hηpos, by linarith, rfl, hθhalf, hθlt,
    le_max_right _ _, le_max_right _ _, hf01, ?_, hderiv, hDc, hlow, ?_, ?_⟩
  · intro x hx
    refine ⟨hmono.derivWithin_nonneg, ?_⟩
    have := hC x hx
    rw [Real.norm_eq_abs] at this
    exact le_trans (le_abs_self _) (le_trans this (le_max_left _ _))
  · intro x hx
    have h4 : HasDerivWithinAt g1 (g2 x) I' x :=
      ((hg1.differentiableOn one_ne_zero) x (hsub' hx)).hasDerivWithinAt
    exact (h4.mono hsub').congr (fun t ht => hDg1 t ht) (hDg1 x hx)
  · intro x hx
    have := hC2b x (hsub' hx)
    rw [Real.norm_eq_abs] at this
    exact le_trans this (le_max_left _ _)

lemma setup' (f : ℝ → ℝ) (hf : Hyp f) : ∃ η θ₀ F₁ β : ℝ, ∃ g2 : ℝ → ℝ, Setup f η θ₀ F₁ β g2 :=
  setup f hf

lemma yd_eq (f : ℝ → ℝ) (hderiv : ∀ x ∈ Icc (0:ℝ) 1, HasDerivWithinAt f (D f x) (Icc 0 1) x)
    (n : ℕ) (x : ℝ) (hx : x ∈ Icc (0:ℝ) 1) : yd f n x = n * f x ^ (n - 1) * D f x :=
  ((hderiv x hx).fun_pow n).derivWithin (uniqueDiffOn_Icc zero_lt_one x hx)

/-- The body of Lemma `g3:lem:shape`. -/
def ShapeBody (f : ℝ → ℝ) (η θ₀ F₁ : ℝ) (n : ℕ) : Prop :=
  let σ : ℝ := n * F₁ * θ₀ ^ (n - 1)
  let a : ℝ := D f 1
  (∀ x ∈ Icc (0 : ℝ) (1 - η), yd f n x ≤ σ) ∧ σ < n * a ∧
  ContinuousOn (yd f n) (Icc (1 - η) 1) ∧ StrictMonoOn (yd f n) (Icc (1 - η) 1) ∧
  yd f n 1 = n * a ∧
  (∀ s : ℝ, n * a ≤ s → levelMeasure (yd f n) s = 0) ∧
  (∀ s : ℝ, 0 ≤ levelMeasure (yd f n) s ∧ levelMeasure (yd f n) s ≤ 1) ∧
  (∀ s : ℝ, σ ≤ s → s < n * a →
    ∃! τ : ℝ, τ ∈ Ioc 0 η ∧ yd f n (1 - τ) = s ∧ levelMeasure (yd f n) s = τ)

lemma shape (f : ℝ → ℝ) (η θ₀ F₁ β : ℝ) (g2 : ℝ → ℝ) (hS : Setup f η θ₀ F₁ β g2) :
    ∃ n₀ : ℕ, 2 ≤ n₀ ∧ ∀ n : ℕ, n₀ ≤ n → ShapeBody f η θ₀ F₁ n := by
  obtain ⟨ha, h1f, hmono, hη, hη1, hθ, hθh, hθ1, haF, hβ, hf01, hDf, hderiv, hDc, hlow, hg2, hg2b⟩ := hS
  set a := D f 1 with ha_def
  have h1mem : (1:ℝ) ∈ Icc (0:ℝ) 1 := ⟨zero_le_one, le_rfl⟩
  have hF : 0 < F₁ := lt_of_lt_of_le ha haF
  have hsub : Icc (1 - η) 1 ⊆ Icc (0:ℝ) 1 := Icc_subset_Icc (by linarith) le_rfl
  obtain ⟨N1, hN1⟩ := exists_pow_lt_of_lt_one (div_pos ha hF) hθ1
  set N2 := ⌈4 * β / a ^ 2⌉₊
  refine ⟨N1 + N2 + 2, by omega, fun n hn => ?_⟩
  obtain ⟨m, rfl⟩ : ∃ m, n = m + 2 := ⟨n - 2, by omega⟩
  have hm1 : N1 ≤ m + 1 := by omega
  have hm2 : (N2 : ℝ) ≤ m := by exact_mod_cast (by omega : N2 ≤ m)
  have hβm : β < (m + 1) * (a ^ 2 / 4) := by
    have := Nat.le_ceil (4 * β / a ^ 2)
    have h4 : 4 * β / a ^ 2 < m + 1 := by linarith
    rw [div_lt_iff₀ (by positivity)] at h4
    linarith
  have hn2 : m + 2 - 1 = m + 1 := by omega
  have hsub : Icc (1 - η) 1 ⊆ Icc (0:ℝ) 1 := Icc_subset_Icc (by linarith) le_rfl
  have hfc : ContinuousOn f (Icc 0 1) := fun x hx => (hderiv x hx).continuousWithinAt
  set h : ℝ → ℝ := fun x => ((m + 2 : ℕ) : ℝ) * f x ^ (m + 1) * D f x with hh
  have hyd : ∀ x ∈ Icc (0:ℝ) 1, yd f (m + 2) x = h x := by
    intro x hx; rw [yd_eq f hderiv _ x hx, hn2]
  have hθ0 : 0 ≤ θ₀ := by linarith
  have hfθ : ∀ x ∈ Icc (1 - η) 1, θ₀ ≤ f x := fun x hx =>
    hθ ▸ hmono (hsub ⟨le_rfl, by linarith⟩) (hsub hx) hx.1
  have hnpos : (0:ℝ) < ((m + 2 : ℕ) : ℝ) := by positivity
  -- strict monotonicity of h
  have hsm : StrictMonoOn h (Icc (1 - η) 1) := by
    have hd : ∀ x ∈ Icc (1 - η) 1, HasDerivWithinAt h
        (((m + 2 : ℕ) : ℝ) * (((m + 1 : ℕ) : ℝ) * f x ^ m * D f x * D f x + f x ^ (m + 1) * g2 x))
        (Icc (1 - η) 1) x := by
      intro x hx
      have h1 := ((hderiv x (hsub hx)).mono hsub).fun_pow (m + 1)
      have h2 := (h1.mul (hg2 x hx)).const_mul (((m + 2 : ℕ) : ℝ))
      convert h2 using 1
      · funext t; simp only [hh, Pi.mul_apply]; ring
      · simp
    apply strictMonoOn_of_hasDerivWithinAt_pos (convex_Icc _ _)
      (fun x hx => (hd x hx).continuousWithinAt)
    · intro x hx
      have hx' := interior_subset hx
      exact (hd x hx').mono interior_subset
    · intro x hx
      have hx' : x ∈ Icc (1 - η) 1 := interior_subset hx
      have hfx := hfθ x hx'
      have hf1 := (hf01 x (hsub hx')).2
      have hD := hlow x hx'
      have hg := hg2b x hx'
      have hfm : 0 < f x ^ m := by apply pow_pos; linarith
      have hfg : -β ≤ f x * g2 x := by
        have : |f x * g2 x| ≤ β := by
          rw [abs_mul, abs_of_nonneg (by linarith)]
          nlinarith [abs_nonneg (g2 x)]
        linarith [neg_abs_le (f x * g2 x)]
      have hD2 : a ^ 2 / 4 ≤ D f x * D f x := by nlinarith
      have key : 0 < ((m + 1 : ℕ) : ℝ) * (D f x * D f x) + f x * g2 x := by
        push_cast
        nlinarith
      have : ((m + 1 : ℕ) : ℝ) * f x ^ m * D f x * D f x + f x ^ (m + 1) * g2 x =
          f x ^ m * (((m + 1 : ℕ) : ℝ) * (D f x * D f x) + f x * g2 x) := by ring
      rw [this]
      positivity
  have hsmy : StrictMonoOn (yd f (m + 2)) (Icc (1 - η) 1) := by
    intro x hx y hy hxy
    rw [hyd x (hsub hx), hyd y (hsub hy)]; exact hsm hx hy hxy
  have hyd1 : yd f (m + 2) 1 = ((m + 2 : ℕ) : ℝ) * a := by
    rw [hyd 1 h1mem]; simp only [hh, h1f, one_pow, mul_one, ha_def]
  have hσle : ∀ x ∈ Icc (0 : ℝ) (1 - η),
      yd f (m + 2) x ≤ ((m + 2 : ℕ) : ℝ) * F₁ * θ₀ ^ (m + 2 - 1) := by
    intro x hx
    have hx01 : x ∈ Icc (0:ℝ) 1 := ⟨hx.1, by linarith [hx.2]⟩
    rw [hyd x hx01, hn2]
    have hfx : f x ≤ θ₀ := hθ ▸ hmono hx01 (hsub ⟨le_rfl, by linarith⟩) hx.2
    have h0 := (hf01 x hx01).1
    have hDx := hDf x hx01
    have : f x ^ (m + 1) ≤ θ₀ ^ (m + 1) := pow_le_pow_left₀ h0 hfx _
    simp only [hh]
    have := mul_le_mul this hDx.2 hDx.1 (by positivity)
    nlinarith
  have hσlt : ((m + 2 : ℕ) : ℝ) * F₁ * θ₀ ^ (m + 2 - 1) < ((m + 2 : ℕ) : ℝ) * a := by
    rw [hn2, mul_assoc]
    apply mul_lt_mul_of_pos_left _ hnpos
    have : θ₀ ^ (m + 1) ≤ θ₀ ^ N1 := pow_le_pow_of_le_one hθ0 hθ1.le hm1
    rw [lt_div_iff₀ hF] at hN1
    nlinarith
  have hle_na : ∀ x ∈ Icc (0:ℝ) 1, yd f (m + 2) x ≤ ((m + 2 : ℕ) : ℝ) * a := by
    intro x hx
    by_cases hx' : x ≤ 1 - η
    · exact le_trans (hσle x ⟨hx.1, hx'⟩) hσlt.le
    · push Not at hx'
      rw [← hyd1]
      exact hsmy.monotoneOn ⟨hx'.le, hx.2⟩ ⟨by linarith, le_rfl⟩ hx.2
  refine ⟨hσle, hσlt, ?_, hsmy, hyd1, ?_, fun s => levelMeasure_le_one _ s, ?_⟩
  · exact (ContinuousOn.mul (ContinuousOn.mul continuousOn_const ((hfc.mono hsub).pow _))
      (hDc.mono hsub)).congr (fun x hx => hyd x (hsub hx))
  · intro s hs
    unfold levelMeasure
    have : {x : ℝ | x ∈ Icc (0:ℝ) 1 ∧ s < yd f (m + 2) x} = ∅ := by
      ext x; simp only [mem_ofPred_eq, mem_empty_iff_false, iff_false, not_and, not_lt]
      intro hx; exact le_trans (hle_na x hx) hs
    rw [this, measure_empty, ENNReal.toReal_zero]
  · intro s hs1 hs2
    have hcont : ContinuousOn (yd f (m + 2)) (Icc (1 - η) 1) :=
      (ContinuousOn.mul (ContinuousOn.mul continuousOn_const ((hfc.mono hsub).pow _))
        (hDc.mono hsub)).congr (fun x hx => hyd x (hsub hx))
    have hlo : yd f (m + 2) (1 - η) ≤ s :=
      le_trans (hσle (1 - η) ⟨by linarith, le_rfl⟩) hs1
    obtain ⟨x, hx, hxs⟩ := intermediate_value_Icc (by linarith : 1 - η ≤ 1) hcont
      ⟨hlo, by rw [hyd1]; exact hs2.le⟩
    have hx1 : x < 1 := by
      rcases eq_or_lt_of_le hx.2 with h | h
      · rw [h, hyd1] at hxs; linarith
      · exact h
    refine ⟨1 - x, ⟨⟨by linarith, by linarith [hx.1]⟩, by rw [sub_sub_cancel]; exact hxs, ?_⟩, ?_⟩
    · unfold levelMeasure
      have : {t : ℝ | t ∈ Icc (0:ℝ) 1 ∧ s < yd f (m + 2) t} = Ioc x 1 := by
        ext t
        simp only [mem_ofPred_eq, mem_Ioc, mem_Icc]
        constructor
        · rintro ⟨⟨ht0, ht1⟩, hts⟩
          by_cases ht : t ≤ 1 - η
          · have := hσle t ⟨ht0, ht⟩; linarith
          · push Not at ht
            refine ⟨?_, ht1⟩
            rw [← hxs] at hts
            exact (hsmy.lt_iff_lt hx ⟨ht.le, ht1⟩).1 hts
        · rintro ⟨htx, ht1⟩
          have ht : t ∈ Icc (1 - η) 1 := ⟨by linarith [hx.1], ht1⟩
          refine ⟨⟨by linarith [hx.1], ht1⟩, ?_⟩
          rw [← hxs]; exact hsmy hx ht htx
      rw [this, Real.volume_Ioc, ENNReal.toReal_ofReal (by linarith)]
    · rintro τ ⟨hτ, hτs, -⟩
      have h1τ : 1 - τ ∈ Icc (1 - η) 1 := ⟨by linarith [hτ.2], by linarith [hτ.1]⟩
      have := hsmy.injOn h1τ hx (hτs.trans hxs.symm)
      linarith

end W01

namespace W01
open MeasureTheory Set Filter Topology

lemma neglog_bounds (w : ℝ) (hw : 1 / 2 ≤ w) :
    1 - w ≤ -Real.log w ∧ -Real.log w ≤ (1 - w) + 2 * (1 - w) ^ 2 := by
  have hw0 : 0 < w := by linarith
  constructor
  · have := Real.log_le_sub_one_of_pos hw0; linarith
  · have h := Real.log_le_sub_one_of_pos (inv_pos.2 hw0)
    rw [Real.log_inv] at h
    have : w⁻¹ - 1 ≤ (1 - w) + 2 * (1 - w) ^ 2 := by
      rw [inv_eq_one_div, div_sub_one hw0.ne', div_le_iff₀ hw0]
      nlinarith [sq_nonneg (1 - w)]
    linarith

lemma log_sub_log_le (d a : ℝ) (ha : 0 < a) (hd : a / 2 ≤ d) :
    |Real.log d - Real.log a| ≤ 2 / a * |d - a| := by
  have hd0 : 0 < d := by linarith
  rw [abs_le]
  constructor
  · -- log a - log d ≤ (a - d)/d ≤ 2|d-a|/a
    have h := Real.log_le_sub_one_of_pos (div_pos ha hd0)
    rw [Real.log_div ha.ne' hd0.ne'] at h
    have h2 : a / d - 1 ≤ 2 / a * |d - a| := by
      rw [div_sub_one hd0.ne', div_le_iff₀ hd0]
      have : a - d ≤ |d - a| := by rw [abs_sub_comm]; exact le_abs_self _
      have h3 : 2 / a * d ≥ 1 := by rw [ge_iff_le, div_mul_eq_mul_div, le_div_iff₀ ha]; linarith
      nlinarith [abs_nonneg (d - a)]
    linarith
  · have h := Real.log_le_sub_one_of_pos (div_pos hd0 ha)
    rw [Real.log_div hd0.ne' ha.ne'] at h
    have h2 : d / a - 1 ≤ 2 / a * |d - a| := by
      rw [div_sub_one ha.ne', div_le_iff₀ ha]
      have : d - a ≤ |d - a| := le_abs_self _
      have h3 : 2 / a * a = 2 := by field_simp
      nlinarith [abs_nonneg (d - a)]
    linarith

lemma window_alg (a β F₁ μ τ w d lam : ℝ) (ha : 0 < a) (hβ : 0 ≤ β) (hμ : 0 < μ)
    (hτ0 : 0 < τ) (hτ1 : τ ≤ 1) (hw : 1 / 2 ≤ w) (hw1 : w ≤ 1) (hd : a / 2 ≤ d) (hdF : d ≤ F₁)
    (E1 : |d - a| ≤ β * τ) (E2 : |w - (1 - a * τ)| ≤ β * τ ^ 2) (E3 : w ≤ 1 - a / 2 * τ)
    (heq : -Real.log w = μ * (lam + Real.log d)) :
    let K := 2 / a * (1 + |Real.log F₁|)
    |τ - μ * (lam + Real.log a) / a| ≤
      1 / a * ((β + 2 * (a + β) ^ 2) * K ^ 2 + 2 * β / a * K) * μ ^ 2 * (1 + |lam|) ^ 2 := by
  intro K
  obtain ⟨hA1, hA2⟩ := neglog_bounds w hw
  have hd0 : 0 < d := by linarith
  have hlogd : Real.log d ≤ |Real.log F₁| :=
    le_trans (Real.log_le_log hd0 hdF) (le_abs_self _)
  -- a priori bound on τ
  have hτK : τ ≤ K * μ * (1 + |lam|) := by
    have h1 : a / 2 * τ ≤ μ * (lam + Real.log d) := by linarith
    have h2 : μ * (lam + Real.log d) ≤ μ * (|lam| + |Real.log F₁|) := by
      apply mul_le_mul_of_nonneg_left _ hμ.le; linarith [le_abs_self lam]
    have h3 : |lam| + |Real.log F₁| ≤ (1 + |Real.log F₁|) * (1 + |lam|) := by
      nlinarith [abs_nonneg lam, abs_nonneg (Real.log F₁)]
    have : τ ≤ 2 / a * (μ * (|lam| + |Real.log F₁|)) := by
      rw [div_mul_eq_mul_div, le_div_iff₀ ha]; linarith
    calc τ ≤ 2 / a * (μ * (|lam| + |Real.log F₁|)) := this
      _ ≤ 2 / a * (μ * ((1 + |Real.log F₁|) * (1 + |lam|))) := by gcongr
      _ = K * μ * (1 + |lam|) := by ring
  have hK0 : 0 ≤ K := by positivity
  set X := μ * (1 + |lam|) with hX
  have hX0 : 0 ≤ X := by positivity
  have hX1 : μ ≤ X := by
    rw [hX]; nlinarith [abs_nonneg lam]
  -- 1 - w ≤ (a + β) τ
  have h1w : 1 - w ≤ (a + β) * τ := by
    have := (abs_le.1 E2).1; nlinarith [mul_nonneg hβ (mul_nonneg hτ0.le (sub_nonneg.2 hτ1))]
  have h1w0 : 0 ≤ 1 - w := by linarith
  -- main identity
  have hmain : a * τ - μ * (lam + Real.log a) =
      (a * τ - (1 - w)) + ((1 - w) - (-Real.log w)) + μ * (Real.log d - Real.log a) := by
    rw [heq]; ring
  have hB : |Real.log d - Real.log a| ≤ 2 / a * (β * τ) :=
    le_trans (log_sub_log_le d a ha hd) (by gcongr)
  have hT1 : |a * τ - (1 - w)| ≤ β * τ ^ 2 := by
    rw [show a * τ - (1 - w) = w - (1 - a * τ) by ring]; exact E2
  have hT2 : |(1 - w) - (-Real.log w)| ≤ 2 * ((a + β) * τ) ^ 2 := by
    rw [abs_le]; constructor <;> nlinarith
  have hT3 : |μ * (Real.log d - Real.log a)| ≤ μ * (2 / a * (β * τ)) := by
    rw [abs_mul, abs_of_pos hμ]; exact mul_le_mul_of_nonneg_left hB hμ.le
  have htot : |a * τ - μ * (lam + Real.log a)| ≤
      (β + 2 * (a + β) ^ 2) * τ ^ 2 + 2 * β / a * (μ * τ) := by
    rw [hmain]
    calc _ ≤ |a * τ - (1 - w)| + |(1 - w) - (-Real.log w)| + |μ * (Real.log d - Real.log a)| :=
          abs_add_three _ _ _
      _ ≤ β * τ ^ 2 + 2 * ((a + β) * τ) ^ 2 + μ * (2 / a * (β * τ)) := by gcongr
      _ = _ := by ring
  have hτ2 : τ ^ 2 ≤ K ^ 2 * X ^ 2 := by
    have : τ ≤ K * X := by rw [hX]; linarith
    calc τ ^ 2 ≤ (K * X) ^ 2 := by gcongr
      _ = _ := by ring
  have hμτ : μ * τ ≤ K * X ^ 2 := by
    have : τ ≤ K * X := by rw [hX]; linarith
    calc μ * τ ≤ X * (K * X) := by gcongr
      _ = _ := by ring
  have hfin : |a * τ - μ * (lam + Real.log a)| ≤
      ((β + 2 * (a + β) ^ 2) * K ^ 2 + 2 * β / a * K) * X ^ 2 := by
    calc _ ≤ _ := htot
      _ ≤ (β + 2 * (a + β) ^ 2) * (K ^ 2 * X ^ 2) + 2 * β / a * (K * X ^ 2) := by
          gcongr
      _ = _ := by ring
  have : τ - μ * (lam + Real.log a) / a = (a * τ - μ * (lam + Real.log a)) / a := by
    field_simp
  rw [this, abs_div, abs_of_pos ha, div_le_iff₀ ha]
  calc _ ≤ _ := hfin
    _ = _ := by rw [hX]; field_simp

lemma local_est (f : ℝ → ℝ) (η θ₀ F₁ β : ℝ) (g2 : ℝ → ℝ) (hS : Setup f η θ₀ F₁ β g2)
    (u : ℝ) (hu : u ∈ Icc (1 - η) 1) :
    |D f u - D f 1| ≤ β * (1 - u) ∧ |f u - (1 - D f 1 * (1 - u))| ≤ β * (1 - u) ^ 2 ∧
      f u ≤ 1 - D f 1 / 2 * (1 - u) := by
  obtain ⟨ha, h1f, hmono, hη, hη1, hθ, hθh, hθ1, haF, hβ, hf01, hDf, hderiv, hDc, hlow, hg2,
    hg2b⟩ := hS
  set a := D f 1
  have hsub : Icc (1 - η) 1 ⊆ Icc (0:ℝ) 1 := Icc_subset_Icc (by linarith) le_rfl
  have hE1 : ∀ t ∈ Icc (1 - η) 1, |D f t - a| ≤ β * (1 - t) := by
    intro t ht
    have := Convex.norm_image_sub_le_of_norm_hasDerivWithin_le hg2
      (fun x hx => by rw [Real.norm_eq_abs]; exact hg2b x hx) (convex_Icc _ _) ht
      (⟨by linarith, le_rfl⟩ : (1:ℝ) ∈ Icc (1 - η) 1)
    rw [Real.norm_eq_abs, Real.norm_eq_abs, abs_of_nonneg (by linarith [ht.2] : (0:ℝ) ≤ 1 - t),
      abs_sub_comm] at this
    exact this
  have huI : Icc u 1 ⊆ Icc (1 - η) 1 := Icc_subset_Icc hu.1 le_rfl
  have hu1 : u ∈ Icc u 1 := ⟨le_rfl, hu.2⟩
  have h11 : (1:ℝ) ∈ Icc u 1 := ⟨hu.2, le_rfl⟩
  refine ⟨hE1 u hu, ?_, ?_⟩
  · have hd : ∀ t ∈ Icc u 1, HasDerivWithinAt (fun t => f t - 1 + a * (1 - t)) (D f t - a)
        (Icc u 1) t := by
      intro t ht
      have h1 := (hderiv t (hsub (huI ht))).mono (huI.trans hsub)
      have h2 := ((hasDerivWithinAt_id t (Icc u 1)).const_sub 1).const_mul a
      convert (h1.sub_const 1).add h2 using 1
      · funext x; simp
      · simp; ring
    have := Convex.norm_image_sub_le_of_norm_hasDerivWithin_le (C := β * (1 - u)) hd
      (fun x hx => by
        rw [Real.norm_eq_abs]
        exact le_trans (hE1 x (huI hx))
          (mul_le_mul_of_nonneg_left (by linarith [hx.1] : 1 - x ≤ 1 - u) hβ))
      (convex_Icc _ _) hu1 h11
    simp only [h1f, sub_self, mul_zero, add_zero, zero_sub, norm_neg, Real.norm_eq_abs] at this
    rw [abs_of_nonneg (by linarith [hu.2] : (0:ℝ) ≤ 1 - u)] at this
    rw [show f u - (1 - a * (1 - u)) = f u - 1 + a * (1 - u) by ring]
    nlinarith
  · have hmon : MonotoneOn (fun t => f t + a / 2 * (1 - t)) (Icc u 1) := by
      have hd : ∀ t ∈ Icc u 1, HasDerivWithinAt (fun t => f t + a / 2 * (1 - t))
          (D f t - a / 2) (Icc u 1) t := by
        intro t ht
        have h1 := (hderiv t (hsub (huI ht))).mono (huI.trans hsub)
        have h2 := ((hasDerivWithinAt_id t (Icc u 1)).const_sub 1).const_mul (a / 2)
        convert h1.add h2 using 1
        · funext x; simp
        · simp; ring
      apply monotoneOn_of_hasDerivWithinAt_nonneg (convex_Icc _ _)
        (fun x hx => (hd x hx).continuousWithinAt)
      · intro x hx; exact (hd x (interior_subset hx)).mono interior_subset
      · intro x hx; linarith [hlow x (huI (interior_subset hx))]
    have := hmon hu1 h11 hu.2
    simp only [h1f, sub_self, mul_zero, add_zero] at this
    linarith

end W01

namespace W01
open MeasureTheory Set Filter Topology

/-- The window estimate, in terms of the localisation data. -/
lemma window (f : ℝ → ℝ) (η θ₀ F₁ β : ℝ) (g2 : ℝ → ℝ) (hS : Setup f η θ₀ F₁ β g2)
    (n₀ : ℕ) (hn₀ : 2 ≤ n₀) (hshape : ∀ n : ℕ, n₀ ≤ n → ShapeBody f η θ₀ F₁ n)
    (n : ℕ) (hn : n₀ ≤ n) (s : ℝ)
    (hs1 : (n : ℝ) * F₁ * θ₀ ^ (n - 1) ≤ s) (hs2 : s < n * D f 1) :
    |levelMeasure (yd f n) s - 1 / ((n : ℝ) - 1) * (Real.log (n / s) + Real.log (D f 1)) / D f 1|
      ≤ 1 / D f 1 * ((β + 2 * (D f 1 + β) ^ 2) * (2 / D f 1 * (1 + |Real.log F₁|)) ^ 2 +
          2 * β / D f 1 * (2 / D f 1 * (1 + |Real.log F₁|))) *
        (1 / ((n : ℝ) - 1)) ^ 2 * (1 + |Real.log (n / s)|) ^ 2 := by
  have hS' := hS
  obtain ⟨ha, h1f, hmono, hη, hη1, hθ, hθh, hθ1, haF, hβ, hf01, hDf, hderiv, hDc, hlow, hg2,
    hg2b⟩ := hS'
  set a := D f 1 with ha_def
  obtain ⟨-, -, -, -, -, -, -, hex⟩ := hshape n hn
  obtain ⟨τ, ⟨hτ, hyds, hm⟩, -⟩ := hex s hs1 hs2
  rw [hm]
  have hsub : Icc (1 - η) 1 ⊆ Icc (0:ℝ) 1 := Icc_subset_Icc (by linarith) le_rfl
  have hu : 1 - τ ∈ Icc (1 - η) 1 := ⟨by linarith [hτ.2], by linarith [hτ.1]⟩
  obtain ⟨E1, E2, E3⟩ := local_est f η θ₀ F₁ β g2 hS (1 - τ) hu
  rw [sub_sub_cancel] at E1 E2 E3
  set w := f (1 - τ) with hw_def
  set d := D f (1 - τ) with hd_def
  have hw : 1 / 2 ≤ w := le_trans hθh (hθ ▸ hmono (hsub ⟨le_rfl, by linarith⟩) (hsub hu) hu.1)
  have hw1 : w ≤ 1 := (hf01 _ (hsub hu)).2
  have hd : a / 2 ≤ d := hlow _ hu
  have hdF : d ≤ F₁ := (hDf _ (hsub hu)).2
  have hn2 : 2 ≤ n := le_trans hn₀ hn
  have hnR : (2:ℝ) ≤ n := by exact_mod_cast hn2
  have hμ : 0 < 1 / ((n : ℝ) - 1) := by apply div_pos one_pos; linarith
  have hw0 : 0 < w := by linarith
  have hd0 : 0 < d := by linarith
  have hyd := yd_eq f hderiv n (1 - τ) (hsub hu)
  rw [hyds] at hyd
  have hcast : ((n - 1 : ℕ) : ℝ) = (n : ℝ) - 1 := by
    rw [Nat.cast_sub (by omega)]; simp
  have heq : -Real.log w = 1 / ((n : ℝ) - 1) * (Real.log (n / s) + Real.log d) := by
    have hs0 : 0 < s := by rw [hyd]; positivity
    have hl : Real.log s = Real.log n + ((n : ℝ) - 1) * Real.log w + Real.log d := by
      rw [hyd, Real.log_mul (by positivity) hd0.ne', Real.log_mul (by positivity)
        (by positivity), Real.log_pow, hcast]
    rw [Real.log_div (by positivity) hs0.ne', hl]
    have hn1 : (n:ℝ) - 1 ≠ 0 := by linarith
    field_simp
    ring
  exact window_alg a β F₁ (1 / ((n : ℝ) - 1)) τ w d (Real.log (n / s)) ha hβ hμ hτ.1
    (by linarith [hτ.2]) hw hw1 hd hdF E1 E2 E3 heq

end W01

namespace W01
open MeasureTheory Set Filter Topology

lemma tail_bound (A : ℝ) (hA : 0 < A) :
    IntegrableOn (fun s => φ s * (Real.log s - Real.log A)) (Ioi A) ∧
      ∫ s in Ioi A, φ s * (Real.log s - Real.log A) ≤ 2 / A := by
  have hsub : Ioi A ⊆ Ioi (0:ℝ) := Ioi_subset_Ioi hA.le
  have i1 : IntegrableOn (fun s => φ s * Real.log s) (Ioi A) := phi_log_integrable.mono_set hsub
  have i2 : IntegrableOn (fun s => φ s * Real.log A) (Ioi A) :=
    (phi_integrable.mono_set hsub).mul_const _
  have hint : IntegrableOn (fun s => φ s * (Real.log s - Real.log A)) (Ioi A) := by
    have := i1.sub i2
    refine this.congr_fun (fun s _ => ?_) measurableSet_Ioi
    simp only [Pi.sub_apply]; ring
  refine ⟨hint, ?_⟩
  have e1 : ∫ s in Ioi A, φ s * Real.log s = 0 - Glog A :=
    integral_Ioi_of_hasDerivAt_of_tendsto' (fun x hx => hasDerivAt_Glog x (lt_of_lt_of_le hA hx))
      i1 Glog_tendsto
  have e2 : ∫ s in Ioi A, φ s = 0 - (A - Real.sqrt (1 + A ^ 2)) :=
    integral_Ioi_of_hasDerivAt_of_tendsto' (f := fun t => t - Real.sqrt (1 + t ^ 2))
      (fun x _ => hasDerivAt_G0 x) (phi_integrable.mono_set hsub) tendsto_G0
  have e3 : ∫ s in Ioi A, φ s * (Real.log s - Real.log A) =
      (∫ s in Ioi A, φ s * Real.log s) - (∫ s in Ioi A, φ s) * Real.log A := by
    rw [← integral_mul_const, ← integral_sub i1 i2]
    congr 1; funext s; ring
  rw [e3, e1, e2]
  set r := Real.sqrt (1 + A ^ 2) with hr
  have hr1 : A ≤ r := Real.le_sqrt_of_sq_le (by nlinarith)
  have hr2 : r ≤ A + 1 := Real.sqrt_le_iff.2 ⟨by linarith, by nlinarith⟩
  have hr0 : 0 < 1 + r := by linarith
  have hlog : Real.log (1 + r) - Real.log A ≤ (1 + r) / A - 1 := by
    rw [← Real.log_div hr0.ne' hA.ne']
    exact Real.log_le_sub_one_of_pos (div_pos hr0 hA)
  have : (1 + r) / A - 1 ≤ 2 / A := by
    rw [div_sub_one hA.ne']; exact div_le_div_of_nonneg_right (by linarith) hA.le
  unfold Glog
  rw [← hr]
  nlinarith

lemma small_bound (s σ : ℝ) (hs : 0 < s) (hsσ : s ≤ σ) (hσ1 : σ ≤ 1) :
    1 + |Real.log s| ≤ 5 * σ ^ (1 / 4 : ℝ) * s ^ (-1 / 2 : ℝ) := by
  obtain ⟨hu, hsu, h2, -⟩ := rpow_quarter_facts s hs
  have h := Real.abs_log_mul_self_rpow_lt s (1 / 4) hs (le_trans hsσ hσ1) (by norm_num)
  rw [h2]
  set u := s ^ (1 / 4 : ℝ)
  have huσ : u ≤ σ ^ (1 / 4 : ℝ) := Real.rpow_le_rpow hs.le hsσ (by norm_num)
  have hu1 : u ≤ 1 := Real.rpow_le_one hs.le (le_trans hsσ hσ1) (by norm_num)
  have hl : |Real.log s| * u ≤ 4 := by
    rw [abs_mul, abs_of_pos hu] at h; norm_num at h; linarith
  rw [← div_eq_mul_inv, le_div_iff₀ (by positivity)]
  have : (1 + |Real.log s|) * u ≤ 5 := by nlinarith
  nlinarith [abs_nonneg (Real.log s)]

end W01

namespace W01
open MeasureTheory Set Filter Topology

/-- The `T` identity: `(log n + log a + 1 - log 2)/(a(n-1)) = ∫ φ(s) μ(log(n/s) + log a)/a`. -/
lemma T_eq (n a : ℝ) (hn : 1 < n) (ha : 0 < a) :
    IntegrableOn (fun s => φ s * (1 / (n - 1) * (Real.log (n / s) + Real.log a) / a)) (Ioi 0) ∧
    (Real.log n + Real.log a + 1 - Real.log 2) / (a * (n - 1)) =
      ∫ s in Ioi (0:ℝ), φ s * (1 / (n - 1) * (Real.log (n / s) + Real.log a) / a) := by
  set c1 := 1 / (n - 1) * (Real.log n + Real.log a) / a
  set c2 := 1 / (n - 1) / a
  have hcongr : EqOn (fun s => φ s * (1 / (n - 1) * (Real.log (n / s) + Real.log a) / a))
      (fun s => c1 * φ s - c2 * (φ s * Real.log s)) (Ioi 0) := by
    intro s hs
    simp only
    rw [Real.log_div (by linarith) (ne_of_gt hs)]
    simp only [c1, c2]; ring
  have hi : IntegrableOn (fun s => c1 * φ s - c2 * (φ s * Real.log s)) (Ioi 0) :=
    (phi_integrable.const_mul c1).sub (phi_log_integrable.const_mul c2)
  refine ⟨hi.congr_fun hcongr.symm measurableSet_Ioi, ?_⟩
  rw [setIntegral_congr_fun measurableSet_Ioi hcongr, integral_sub
    (phi_integrable.const_mul c1) (phi_log_integrable.const_mul c2), integral_const_mul,
    integral_const_mul, phi_integral, phi_log_integral]
  have : n - 1 ≠ 0 := by linarith
  simp only [c1, c2]
  field_simp
  ring

end W01

namespace W01
open MeasureTheory Set Filter Topology

/-- The three dominating functions. -/
noncomputable def B1 (n C1 μ : ℝ) (s : ℝ) : ℝ :=
  3 * C1 * μ ^ 2 * (φ s * (1 + Real.log n ^ 2) + φ s * Real.log s ^ 2)
noncomputable def B2 (n K1 σ : ℝ) (s : ℝ) : ℝ :=
  (Ioc (0:ℝ) 1).indicator
    (fun s => 5 * K1 * (1 + Real.log n) * σ ^ (1 / 4 : ℝ) * s ^ (-1 / 2 : ℝ)) s
noncomputable def B3 (n a μ : ℝ) (s : ℝ) : ℝ :=
  (Ioi (n * a)).indicator (fun s => μ / a * (φ s * (Real.log s - Real.log (n * a)))) s

lemma pt_bound (n a μ C1 σ mv s : ℝ) (hs : 0 < s) (ha : 0 < a) (hμ0 : 0 < μ) (hμ1 : μ ≤ 1)
    (hn1 : 1 ≤ n) (hC1 : 0 ≤ C1) (hσ0 : 0 < σ) (hσ1 : σ ≤ 1) (hm0 : 0 ≤ mv) (hm1 : mv ≤ 1)
    (hmA : n * a ≤ s → mv = 0)
    (hwin : σ ≤ s → s < n * a →
      |mv - μ * (Real.log (n / s) + Real.log a) / a| ≤ C1 * μ ^ 2 * (1 + |Real.log (n / s)|) ^ 2) :
    ‖φ s * (mv - μ * (Real.log (n / s) + Real.log a) / a)‖ ≤
      B1 n C1 μ s + B2 n (1 + (1 + |Real.log a|) / a) σ s + B3 n a μ s := by
  have hφ0 := phi_pos s hs.le
  have hφ1 := phi_le_one s hs.le
  have hln : 0 ≤ Real.log n := Real.log_nonneg hn1
  have hlam : Real.log (n / s) = Real.log n - Real.log s :=
    Real.log_div (by linarith) hs.ne'
  have hB1 : 0 ≤ B1 n C1 μ s := by unfold B1; positivity
  have hB2 : 0 ≤ B2 n (1 + (1 + |Real.log a|) / a) σ s := by
    unfold B2; apply indicator_nonneg; intro x hx
    have : 0 < x := hx.1
    positivity
  have hB3 : 0 ≤ B3 n a μ s := by
    unfold B3; apply indicator_nonneg; intro x hx
    have hx' : n * a < x := hx
    have : Real.log (n * a) ≤ Real.log x := Real.log_le_log (by positivity) hx'.le
    have : 0 < φ x := phi_pos x (by nlinarith)
    have : 0 ≤ Real.log x - Real.log (n * a) := by linarith
    positivity
  rw [Real.norm_eq_abs, abs_mul, abs_of_pos hφ0]
  rcases lt_or_ge s σ with h1 | h1
  · -- small s
    have hsI : s ∈ Ioc (0:ℝ) 1 := ⟨hs, by linarith⟩
    have hB2v : B2 n (1 + (1 + |Real.log a|) / a) σ s = 5 * (1 + (1 + |Real.log a|) / a) *
        (1 + Real.log n) * σ ^ (1 / 4 : ℝ) * s ^ (-1 / 2 : ℝ) := by
      unfold B2; rw [indicator_of_mem hsI]
    have hsm := small_bound s σ hs h1.le hσ1
    set K1 := 1 + (1 + |Real.log a|) / a
    set ls := |Real.log s|
    have hh : |μ * (Real.log (n / s) + Real.log a) / a| ≤ (Real.log n + ls + |Real.log a|) / a := by
      rw [abs_div, abs_of_pos ha, abs_mul, abs_of_pos hμ0]
      apply div_le_div_of_nonneg_right _ ha.le
      have : |Real.log (n / s) + Real.log a| ≤ Real.log n + ls + |Real.log a| := by
        rw [hlam]
        calc _ ≤ |Real.log n - Real.log s| + |Real.log a| := abs_add_le _ _
          _ ≤ |Real.log n| + |Real.log s| + |Real.log a| := by gcongr; exact abs_sub _ _
          _ = _ := by rw [abs_of_nonneg hln]
      calc μ * |Real.log (n / s) + Real.log a| ≤ 1 * (Real.log n + ls + |Real.log a|) :=
            mul_le_mul hμ1 this (abs_nonneg _) zero_le_one
        _ = _ := one_mul _
    have hdiff : |mv - μ * (Real.log (n / s) + Real.log a) / a| ≤
        1 + (Real.log n + ls + |Real.log a|) / a := by
      calc _ ≤ |mv| + |μ * (Real.log (n / s) + Real.log a) / a| := abs_sub _ _
        _ ≤ _ := by rw [abs_of_nonneg hm0]; linarith
    have hK : 1 + (Real.log n + ls + |Real.log a|) / a ≤ K1 * (1 + Real.log n) * (1 + ls) := by
      have hls : 0 ≤ ls := abs_nonneg _
      have hla : 0 ≤ |Real.log a| := abs_nonneg _
      have e : K1 * (1 + Real.log n) * (1 + ls) = (1 + Real.log n) * (1 + ls) +
          ((1 + |Real.log a|) * (1 + Real.log n) * (1 + ls)) / a := by
        simp only [K1]; field_simp
      rw [e]
      have h2 : Real.log n + ls + |Real.log a| ≤ (1 + |Real.log a|) * (1 + Real.log n) * (1 + ls) := by
        nlinarith [mul_nonneg hln hls, mul_nonneg hla hln, mul_nonneg hla hls,
          mul_nonneg (mul_nonneg hla hln) hls]
      have h3 : (Real.log n + ls + |Real.log a|) / a ≤
          ((1 + |Real.log a|) * (1 + Real.log n) * (1 + ls)) / a :=
        div_le_div_of_nonneg_right h2 ha.le
      nlinarith [mul_nonneg hln hls]
    have hK1 : 0 ≤ K1 := by positivity
    calc φ s * |mv - μ * (Real.log (n / s) + Real.log a) / a| ≤
          1 * (K1 * (1 + Real.log n) * (1 + ls)) :=
          mul_le_mul hφ1 (hdiff.trans hK) (abs_nonneg _) zero_le_one
      _ ≤ K1 * (1 + Real.log n) * (5 * σ ^ (1 / 4 : ℝ) * s ^ (-1 / 2 : ℝ)) := by
          rw [one_mul]; gcongr
      _ = B2 n K1 σ s := by rw [hB2v]; ring
      _ ≤ _ := by linarith
  rcases lt_or_ge s (n * a) with h2 | h2
  · -- window
    have hw := hwin h1 h2
    have hl : |Real.log (n / s)| ≤ Real.log n + |Real.log s| := by
      rw [hlam]
      calc _ ≤ |Real.log n| + |Real.log s| := abs_sub _ _
        _ = _ := by rw [abs_of_nonneg hln]
    have hsq : (1 + |Real.log (n / s)|) ^ 2 ≤ 3 * (1 + Real.log n ^ 2 + Real.log s ^ 2) := by
      have h0 : 0 ≤ 1 + |Real.log (n / s)| := by positivity
      have : (1 + |Real.log (n / s)|) ^ 2 ≤ (1 + Real.log n + |Real.log s|) ^ 2 :=
        pow_le_pow_left₀ h0 (by linarith) 2
      have h3 : (1 + Real.log n + |Real.log s|) ^ 2 ≤ 3 * (1 + Real.log n ^ 2 + |Real.log s| ^ 2) := by
        nlinarith [sq_nonneg (1 - Real.log n), sq_nonneg (1 - |Real.log s|),
          sq_nonneg (Real.log n - |Real.log s|)]
      rw [sq_abs] at h3
      linarith
    have : φ s * |mv - μ * (Real.log (n / s) + Real.log a) / a| ≤ B1 n C1 μ s := by
      calc _ ≤ φ s * (C1 * μ ^ 2 * (1 + |Real.log (n / s)|) ^ 2) :=
            mul_le_mul_of_nonneg_left hw hφ0.le
        _ ≤ φ s * (C1 * μ ^ 2 * (3 * (1 + Real.log n ^ 2 + Real.log s ^ 2))) := by gcongr
        _ = B1 n C1 μ s := by unfold B1; ring
    linarith
  · -- large s
    have hm := hmA h2
    rw [hm, zero_sub, abs_neg]
    have hlog : Real.log (n / s) + Real.log a = -(Real.log s - Real.log (n * a)) := by
      rw [hlam, Real.log_mul (by linarith) ha.ne']; ring
    rw [hlog]
    rcases eq_or_lt_of_le h2 with h3 | h3
    · have h0 : Real.log s - Real.log (n * a) = 0 := by rw [h3, sub_self]
      rw [h0, neg_zero, mul_zero, zero_div, abs_zero, mul_zero]
      linarith
    · have hB3v : B3 n a μ s = μ / a * (φ s * (Real.log s - Real.log (n * a))) := by
        unfold B3; rw [indicator_of_mem (show s ∈ Ioi (n * a) from h3)]
      have hpos : 0 ≤ Real.log s - Real.log (n * a) := by
        have := Real.log_le_log (by positivity) h3.le; linarith
      rw [abs_div, abs_of_pos ha, abs_mul, abs_of_pos hμ0, abs_neg, abs_of_nonneg hpos]
      have : φ s * (μ * (Real.log s - Real.log (n * a)) / a) = B3 n a μ s := by
        rw [hB3v]; field_simp
      linarith

end W01

namespace W01
open MeasureTheory Set Filter Topology

lemma rpow_half_integrable : IntegrableOn (fun s : ℝ => s ^ (-1 / 2 : ℝ)) (Ioc 0 1) := by
  have := intervalIntegral.intervalIntegrable_rpow' (a := 0) (b := 1) (r := -1 / 2) (by norm_num)
  rwa [intervalIntegrable_iff_integrableOn_Ioc_of_le zero_le_one] at this

lemma B_integrals (n C1 μ K1 σ a : ℝ) (hn : 1 ≤ n) (ha : 0 < a) (hμ : 0 ≤ μ) :
    IntegrableOn (fun s => B1 n C1 μ s + B2 n K1 σ s + B3 n a μ s) (Ioi 0) ∧
    ∫ s in Ioi (0:ℝ), (B1 n C1 μ s + B2 n K1 σ s + B3 n a μ s) ≤
      3 * C1 * μ ^ 2 * (1 + Real.log n ^ 2 + ∫ s in Ioi (0:ℝ), φ s * Real.log s ^ 2) +
      5 * K1 * (1 + Real.log n) * σ ^ (1 / 4 : ℝ) * (∫ s in Ioc (0:ℝ) 1, s ^ (-1 / 2 : ℝ)) +
      μ / a * (2 / (n * a)) := by
  have hA : 0 < n * a := by nlinarith
  obtain ⟨hti, htb⟩ := tail_bound (n * a) hA
  have iB1 : IntegrableOn (B1 n C1 μ) (Ioi 0) := by
    have := ((phi_integrable.mul_const (1 + Real.log n ^ 2)).add phi_log2_integrable).const_mul
      (3 * C1 * μ ^ 2)
    exact IntegrableOn.congr_fun this (fun s _ => by simp [B1]) measurableSet_Ioi
  have ig2 : IntegrableOn (fun s : ℝ => 5 * K1 * (1 + Real.log n) * σ ^ (1 / 4 : ℝ) *
      s ^ (-1 / 2 : ℝ)) (Ioc 0 1) := rpow_half_integrable.const_mul _
  have iB2 : IntegrableOn (B2 n K1 σ) (Ioi 0) :=
    (ig2.integrable_indicator measurableSet_Ioc).integrableOn
  have ig3 : IntegrableOn (fun s => μ / a * (φ s * (Real.log s - Real.log (n * a)))) (Ioi (n * a)) :=
    hti.const_mul _
  have iB3 : IntegrableOn (B3 n a μ) (Ioi 0) :=
    (ig3.integrable_indicator measurableSet_Ioi).integrableOn
  refine ⟨(iB1.add iB2).add iB3, ?_⟩
  rw [integral_add (f := fun s => B1 n C1 μ s + B2 n K1 σ s) (g := B3 n a μ) (iB1.add iB2) iB3,
    integral_add (f := B1 n C1 μ) (g := B2 n K1 σ) iB1 iB2]
  have e1 : ∫ s in Ioi (0:ℝ), B1 n C1 μ s =
      3 * C1 * μ ^ 2 * (1 + Real.log n ^ 2 + ∫ s in Ioi (0:ℝ), φ s * Real.log s ^ 2) := by
    unfold B1
    rw [integral_const_mul, integral_add (phi_integrable.mul_const _) phi_log2_integrable,
      integral_mul_const, phi_integral]
    ring
  have e2 : ∫ s in Ioi (0:ℝ), B2 n K1 σ s =
      5 * K1 * (1 + Real.log n) * σ ^ (1 / 4 : ℝ) * (∫ s in Ioc (0:ℝ) 1, s ^ (-1 / 2 : ℝ)) := by
    unfold B2
    rw [setIntegral_indicator measurableSet_Ioc,
      inter_eq_right.2 Ioc_subset_Ioi_self, integral_const_mul]
  have e3 : ∫ s in Ioi (0:ℝ), B3 n a μ s ≤ μ / a * (2 / (n * a)) := by
    unfold B3
    rw [setIntegral_indicator measurableSet_Ioi, inter_eq_right.2 (Ioi_subset_Ioi hA.le),
      integral_const_mul]
    exact mul_le_mul_of_nonneg_left htb (by positivity)
  rw [e1, e2]
  linarith

end W01

namespace W01
open MeasureTheory Set Filter Topology

lemma final_arith (n ln μ C1 M2 K1 q I0 a : ℝ) (hn : 3 ≤ n) (hln : 1 ≤ ln) (hμ0 : 0 ≤ μ)
    (hμ : μ ≤ 2 / n) (hC1 : 0 ≤ C1) (hM2 : 0 ≤ M2) (hK1 : 0 ≤ K1) (hq0 : 0 ≤ q)
    (hq : q ≤ 1 / n ^ 2) (hI0 : 0 ≤ I0) (ha : 0 < a) :
    3 * C1 * μ ^ 2 * (1 + ln ^ 2 + M2) + 5 * K1 * (1 + ln) * q * I0 + μ / a * (2 / (n * a)) ≤
      (12 * C1 * (2 + M2) + 10 * K1 * I0 + 4 / a ^ 2) * ln ^ 2 / n ^ 2 := by
  have hn0 : 0 < n := by linarith
  set t := ln ^ 2 with ht
  have ht1 : 1 ≤ t := by nlinarith
  have hlt : ln ≤ t := by nlinarith
  have hμ2 : μ ^ 2 ≤ 4 / n ^ 2 := by
    calc μ ^ 2 ≤ (2 / n) ^ 2 := pow_le_pow_left₀ hμ0 hμ 2
      _ = 4 / n ^ 2 := by ring
  have T1 : 3 * C1 * μ ^ 2 * (1 + t + M2) ≤ 12 * C1 * (2 + M2) * t / n ^ 2 := by
    have h1 : 1 + t + M2 ≤ (2 + M2) * t := by nlinarith
    calc 3 * C1 * μ ^ 2 * (1 + t + M2) ≤ 3 * C1 * (4 / n ^ 2) * ((2 + M2) * t) := by
          gcongr
      _ = 12 * C1 * (2 + M2) * t / n ^ 2 := by ring
  have T2 : 5 * K1 * (1 + ln) * q * I0 ≤ 10 * K1 * I0 * t / n ^ 2 := by
    calc 5 * K1 * (1 + ln) * q * I0 ≤ 5 * K1 * (2 * t) * (1 / n ^ 2) * I0 := by
          gcongr; linarith
      _ = 10 * K1 * I0 * t / n ^ 2 := by ring
  have T3 : μ / a * (2 / (n * a)) ≤ 4 / a ^ 2 * t / n ^ 2 := by
    calc μ / a * (2 / (n * a)) ≤ (2 / n) / a * (2 / (n * a)) := by gcongr
      _ = 4 / a ^ 2 * 1 / n ^ 2 := by field_simp; norm_num
      _ ≤ 4 / a ^ 2 * t / n ^ 2 := by gcongr
  calc _ ≤ 12 * C1 * (2 + M2) * t / n ^ 2 + 10 * K1 * I0 * t / n ^ 2 + 4 / a ^ 2 * t / n ^ 2 := by
        linarith
    _ = _ := by ring

lemma main_thm (f : ℝ → ℝ) (hf : Hyp f) :
    ∃ C : ℝ, ∃ N : ℕ, ∀ n : ℕ, N ≤ n →
      |(2 - L f n) - (Real.log n + Real.log (D f 1) + 1 - Real.log 2) / (D f 1 * ((n : ℝ) - 1))|
        ≤ C * Real.log n ^ 2 / (n : ℝ) ^ 2 := by
  have hf' := hf
  obtain ⟨hC1, hmono0, h0, h1, -, -⟩ := hf'
  obtain ⟨η, θ₀, F₁, β, g2, hS⟩ := setup' f hf
  obtain ⟨n₀, hn₀, hshape⟩ := shape f η θ₀ F₁ β g2 hS
  have hS' := hS
  obtain ⟨ha, -, -, -, -, -, hθh, hθ1, haF, hβ, hf01, -⟩ := hS'
  set a := D f 1 with ha_def
  have hF : 0 < F₁ := lt_of_lt_of_le ha haF
  have hθ0 : 0 < θ₀ := by linarith
  set K := 2 / a * (1 + |Real.log F₁|) with hK
  set C1 := 1 / a * ((β + 2 * (a + β) ^ 2) * K ^ 2 + 2 * β / a * K) with hC1def
  have hK0 : 0 ≤ K := by positivity
  have hC1nn : 0 ≤ C1 := by
    apply mul_nonneg (by positivity)
    apply add_nonneg
    · exact mul_nonneg (add_nonneg hβ (by positivity)) (by positivity)
    · exact mul_nonneg (div_nonneg (by linarith) ha.le) hK0
  set K1 := 1 + (1 + |Real.log a|) / a with hK1
  have hK1nn : 0 ≤ K1 := by positivity
  set M2 := ∫ s in Ioi (0:ℝ), φ s * Real.log s ^ 2 with hM2def
  have hM2 : 0 ≤ M2 := setIntegral_nonneg measurableSet_Ioi
    (fun s hs => mul_nonneg (phi_pos s (le_of_lt hs)).le (sq_nonneg _))
  set I0 := ∫ s in Ioc (0:ℝ) 1, s ^ (-1 / 2 : ℝ) with hI0def
  have hI0 : 0 ≤ I0 := setIntegral_nonneg measurableSet_Ioc
    (fun s hs => Real.rpow_nonneg hs.1.le _)
  have hθabs : |θ₀| < 1 := by rw [abs_of_pos hθ0]; exact hθ1
  obtain ⟨N3, hN3⟩ := eventually_atTop.1
    ((tendsto_pow_const_mul_const_pow_of_abs_lt_one 9 hθabs).eventually
      (gt_mem_nhds (div_pos hθ0 hF)))
  refine ⟨12 * C1 * (2 + M2) + 10 * K1 * I0 + 4 / a ^ 2, max n₀ (max N3 3), fun n hn => ?_⟩
  have hnn₀ : n₀ ≤ n := le_trans (le_max_left _ _) hn
  have hnN3 : N3 ≤ n := le_trans (le_trans (le_max_left _ _) (le_max_right _ _)) hn
  have hn3 : 3 ≤ n := le_trans (le_trans (le_max_right _ _) (le_max_right _ _)) hn
  have hnR : (3:ℝ) ≤ n := by exact_mod_cast hn3
  have hn0 : (0:ℝ) < n := by linarith
  set μ := 1 / ((n : ℝ) - 1) with hμdef
  have hμ0 : 0 < μ := by apply div_pos one_pos; linarith
  have hμ1 : μ ≤ 1 := by rw [div_le_one (by linarith)]; linarith
  have hμ2 : μ ≤ 2 / n := by
    rw [hμdef, div_le_div_iff₀ (by linarith) hn0]; linarith
  have hln : 1 ≤ Real.log n := by
    rw [← Real.log_exp 1]
    apply Real.log_le_log (Real.exp_pos 1)
    have := Real.exp_one_lt_d9; linarith
  -- σ
  set σ := (n : ℝ) * F₁ * θ₀ ^ (n - 1) with hσdef
  have hσ0 : 0 < σ := by positivity
  have hσn : σ * (n : ℝ) ^ 8 < 1 := by
    have h9 := hN3 n hnN3
    have hp : θ₀ ^ n = θ₀ ^ (n - 1) * θ₀ := by
      rw [← pow_succ]; congr 1; omega
    rw [hp, lt_div_iff₀ hF] at h9
    have : (n : ℝ) ^ 9 * θ₀ ^ (n - 1) * F₁ < 1 := by
      have h' : (n : ℝ) ^ 9 * θ₀ ^ (n - 1) * F₁ * θ₀ < 1 * θ₀ := by nlinarith
      exact lt_of_mul_lt_mul_right h' hθ0.le
    calc σ * (n : ℝ) ^ 8 = (n : ℝ) ^ 9 * θ₀ ^ (n - 1) * F₁ := by rw [hσdef]; ring
      _ < 1 := this
  have hn8 : (1:ℝ) ≤ (n : ℝ) ^ 8 := one_le_pow₀ (by linarith)
  have hσ1 : σ ≤ 1 := by nlinarith
  have hσq : σ ^ (1 / 4 : ℝ) ≤ 1 / (n : ℝ) ^ 2 := by
    have hle : σ ≤ (1 / (n : ℝ) ^ 2) ^ 4 := by
      rw [div_pow, one_pow, ← pow_mul, le_div_iff₀ (by positivity)]; norm_num; linarith
    calc σ ^ (1 / 4 : ℝ) ≤ ((1 / (n : ℝ) ^ 2) ^ 4) ^ (1 / 4 : ℝ) :=
          Real.rpow_le_rpow hσ0.le hle (by norm_num)
      _ = 1 / (n : ℝ) ^ 2 := by
          have := Real.pow_rpow_inv_natCast (x := 1 / (n : ℝ) ^ 2) (n := 4) (by positivity)
            (by norm_num)
          have e : (1 / 4 : ℝ) = ((4 : ℕ) : ℝ)⁻¹ := by norm_num
          rw [e, this]
  have hσq0 : 0 ≤ σ ^ (1 / 4 : ℝ) := Real.rpow_nonneg hσ0.le _
  -- layer cake
  have hlayer : 2 - L f n = ∫ s in Ioi (0:ℝ), φ s * levelMeasure (yd f n) s := by
    have hy : ContDiffOn ℝ 1 (fun t => f t ^ n) (Icc 0 1) := hC1.pow n
    have hm : MonotoneOn (fun t => f t ^ n) (Icc 0 1) := by
      intro x hx y hy hxy
      exact pow_le_pow_left₀ (hf01 x hx).1 (hmono0 hx hy hxy) n
    have h10 : (fun t => f t ^ n) 1 - (fun t => f t ^ n) 0 = 1 := by
      simp only [h0, h1, one_pow]; rw [zero_pow (by omega)]; ring
    exact layer (fun t => f t ^ n) hy hm h10
  obtain ⟨hTi, hT⟩ := T_eq (n : ℝ) a (by linarith) ha
  rw [hlayer, hT, ← integral_sub (phi_level_integrable _) hTi]
  obtain ⟨hBi, hBint⟩ := B_integrals (n : ℝ) C1 μ K1 σ a (by linarith) ha hμ0.le
  obtain ⟨-, -, -, -, -, hmA, -, -⟩ := hshape n hnn₀
  have hbound : ‖∫ s in Ioi (0:ℝ), (φ s * levelMeasure (yd f n) s -
      φ s * (1 / ((n : ℝ) - 1) * (Real.log (n / s) + Real.log a) / a))‖ ≤
      ∫ s in Ioi (0:ℝ), (B1 n C1 μ s + B2 n K1 σ s + B3 n a μ s) := by
    apply norm_integral_le_of_norm_le hBi
    refine (ae_restrict_iff' measurableSet_Ioi).2 (Eventually.of_forall fun s hs => ?_)
    rw [← mul_sub]
    exact pt_bound (n : ℝ) a μ C1 σ (levelMeasure (yd f n) s) s hs ha hμ0 hμ1 (by linarith)
      hC1nn hσ0 hσ1 (levelMeasure_le_one _ s).1 (levelMeasure_le_one _ s).2 (hmA s)
      (fun h1 h2 => window f η θ₀ F₁ β g2 hS n₀ hn₀ hshape n hnn₀ s h1 h2)
  rw [Real.norm_eq_abs] at hbound
  refine le_trans hbound (le_trans hBint ?_)
  exact final_arith (n : ℝ) (Real.log n) μ C1 M2 K1 (σ ^ (1 / 4 : ℝ)) I0 a hnR hln hμ0.le hμ2
    hC1nn hM2 hK1nn hσq0 hσq hI0 ha

end W01

theorem W01.proof_1 : W01.stmt_1 := fun s hs =>
  ⟨W01.phi_pos s hs.le, le_min (W01.phi_le_one s hs.le) (W01.phi_le_inv s hs)⟩
theorem W01.proof_2 : W01.stmt_2 :=
  ⟨W01.phi_integrable, W01.phi_integral, W01.phi_log_integrable, W01.phi_log2_integrable,
    W01.phi_log_integral⟩
theorem W01.proof_3 : W01.stmt_3 := fun y hy hm h1 =>
  ⟨W01.phi_level_integrable _, W01.layer y hy hm h1⟩
theorem W01.proof_4 : W01.stmt_4 := by
  intro f hf
  obtain ⟨η, θ₀, F₁, β, g2, hS⟩ := W01.setup' f hf
  obtain ⟨n₀, -, hn⟩ := W01.shape f η θ₀ F₁ β g2 hS
  obtain ⟨-, -, -, hη, hη1, hθ, hθh, hθ1, haF, -⟩ := hS
  exact ⟨η, θ₀, F₁, n₀, hη, hη1, hθ, hθh, hθ1, haF, fun n h => hn n h⟩
theorem W01.proof_5 : W01.stmt_5 := by
  intro f hf
  obtain ⟨η, θ₀, F₁, β, g2, hS⟩ := W01.setup' f hf
  obtain ⟨n₀, hn₀, hshape⟩ := W01.shape f η θ₀ F₁ β g2 hS
  have hS' := hS
  obtain ⟨ha, -, -, -, -, -, hθh, hθ1, haF, -⟩ := hS'
  refine ⟨θ₀, F₁, 1 / D f 1 * ((β + 2 * (D f 1 + β) ^ 2) * (2 / D f 1 * (1 + |Real.log F₁|)) ^ 2 +
    2 * β / D f 1 * (2 / D f 1 * (1 + |Real.log F₁|))), n₀, by linarith, hθ1, by linarith,
    fun n hn s => ?_⟩
  intro σ a μ lam hs1 hs2
  exact W01.window f η θ₀ F₁ β g2 hS n₀ hn₀ hshape n hn s hs1 hs2
theorem W01.proof_6 : W01.stmt_6 := W01.main_thm
