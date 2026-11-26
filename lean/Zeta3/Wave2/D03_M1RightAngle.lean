import Mathlib
import Zeta3.Wave2.D00_Dilog
/-!
# D03 -- the four-dilogarithm formula at the right angle reproduces `M_1` of the squircle

Target id: D03.
Paper labels: `prop:M1short` (main.tex 1623-1629), `thm:M1` (main.tex 889-895), the definition
of `M_0(α)` and of the corner family `M(u;α)` (main.tex 1466-1474), and the integral
`eq:M1t` for the squircle's `M_1` (main.tex 858-861).

Quoted verbatim (main.tex 1466-1474):

  The coefficients $M_k$ of the squircle come from an integral over one right-angled corner.
  ... its length is $2+2u\,M(u;\alpha)$, where
  \begin{multline*}
  M(u;\alpha)=\int_0^1\Bigl[(1+g)^{-1-u}\sqrt{1+g^{2-2u}-2\cos\alpha\,g^{1-u}}\\
  -1\Bigr]g^{u-1}\dd g=\sum_{k\ge0}M_k(\alpha)u^k,
  \end{multline*}
  which gives back the squircle's $M_k$ at $\alpha=\pi/2$. The first coefficient is elementary:
  with $c=\cos\frac\alpha2$,
  \[
  M_0(\alpha)=-(1+c)\log(1+c)-(1-c)\log(1-c).
  \]

Quoted verbatim (main.tex 1623-1629, `prop:M1short`):

  For $0<\alpha<\pi$, with $\lambda=\log2$ and $M_0(\alpha)$ as above,
  \begin{multline*}
  M_1(\alpha)=\zeta(2)-\lambda^2-\lambda M_0(\alpha)-(1-c)\bigl[\Li_2(c)+\Li_2(\tfrac{1-c}2)\bigr]\\
  -(1+c)\bigl[\Li_2(-c)+\Li_2(\tfrac{1+c}2)\bigr].
  \end{multline*}

Quoted verbatim (main.tex 879-895, `M_0` and `thm:M1`; `a=1+\sqrt2`, `L=\log a`, `λ=\log 2`):

  M_0=\lambda-\sqrt2\,L=-0.55330329972051571737\ldots,
  \begin{theorem}[$M_1$; answers Problem~1 of Ahmed~\cite{ahmed}]\label{thm:M1}
  \begin{multline*}
  M_1=\Bigl(\frac{19\sqrt2}{48}-\frac1{24}\Bigr)\pi^2+\frac{\lambda^2}{2}+2\sqrt2\,\lambda L\\
  -\Bigl(1+\frac{5\sqrt2}{2}\Bigr)L^2-5\sqrt2\,\Li_2(\sqrt2-1).
  \end{multline*}
  Numerically $\kappa=8M_1=2.02095276965635835695\ldots$
  \end{theorem}

Quoted verbatim (main.tex 858-861, eq. `eq:M1t`):

  M_1=\int_0^1\frac{\log t}{t}\Bigl(\frac{1-t}{q}-1\Bigr)dt,\qquad q=\sqrt{t^2+(1-t)^2},

## Normalisation

The paper says (line 1471) that the corner family "gives back the squircle's `M_k` at
`α = π/2`", with no factor or shift: `M_k(π/2) = M_k`. (Checked: at `α = π/2` the corner
integrand is the squircle integrand, and both `M_0` and `M_1` agree numerically to 40 digits,
see below.) So the right-angle claim is `M_1(π/2) = M_1`, i.e. the four-dilogarithm expression
at `c = cos(π/4) = √2/2` equals the closed form of `thm:M1`, and likewise
`M_0(π/2) = λ - √2 L`.

## Conventions

* `Li_2` is `D00.Li2`. All four arguments `c, (1-c)/2, -c, (1+c)/2` lie in `(-1,1)` at
  `c = √2/2`, as does `√2-1`, so only the genuine dilogarithm is used.
* `ζ(2)` is written `π²/6` (Mathlib: `riemannZeta_two`).
* All logarithms are of positive numbers (`1 ± c > 0`, `1+√2 > 0`, `2 > 0`).
* `c` is defined as the paper does, `c = cos(α/2)`, and the formula is evaluated at `α = π/2`.
* stmt_4 uses `eq:M1t` as the definition of the squircle's `M_1` (an integral defined via a
  power-series coefficient is not used). The integrand is integrable on `(0,1)` (it is
  `O(t log t)` at `0` and bounded near `1`); stmt_4 asserts that integrability, so the value is
  not Lean's junk `0`. At `t = 0` Lean's `log 0 / 0 = 0` is a single point, irrelevant.

## Statements

* stmt_1: `M_0(π/2) = λ - √2 L` (the corner `M_0` at the right angle is the squircle `M_0`).
* stmt_2: the four-dilogarithm right side of `prop:M1short` at `α = π/2` equals the closed form
  of `thm:M1`.
* stmt_3: `8 ×` that value is `2.02095276965635835695…` to within `10⁻²⁰` (the paper's
  numerical value of `κ`).
* stmt_4: the integral `eq:M1t` is integrable on `(0,1)` and equals the four-dilogarithm value
  at `α = π/2`.

Numerical check (mpmath, 40 digits): four-dilog value = 0.25261909620704479461836352096565486053,
`thm:M1` closed form agrees to 4·10⁻⁴⁰, `eq:M1t` integral agrees to 30 digits,
`8M_1 = 2.0209527696563583569469…`; `M_0(π/2) = λ-√2L = -0.55330329972051571737…`.

NOT COVERED:
* `prop:M1short` for general `0 < α < π`, and the corner family `M(u;α)` itself (defining
  `M_1(α)` as a Taylor coefficient of a parameter integral); only its value at `α = π/2` is
  formalized, through the four-dilogarithm right side, the closed form, and `eq:M1t`.
* The statement "`M(u;π/2)` gives back the squircle's `M_k`" for `k ≥ 2`.
-/
-- BEGIN STATEMENTS
namespace D03

open Real

/-- `λ = log 2`. -/
noncomputable def lam : ℝ := Real.log 2

/-- `L = log(1+√2)`. -/
noncomputable def L : ℝ := Real.log (1 + Real.sqrt 2)

/-- `c = cos(α/2)` (main.tex 1472). -/
noncomputable def c (α : ℝ) : ℝ := Real.cos (α / 2)

/-- `M_0(α) = -(1+c)log(1+c) - (1-c)log(1-c)` (main.tex 1473). -/
noncomputable def M0 (α : ℝ) : ℝ :=
  -(1 + c α) * Real.log (1 + c α) - (1 - c α) * Real.log (1 - c α)

/-- Right side of `prop:M1short`:
`ζ(2) - λ² - λM_0(α) - (1-c)[Li₂(c)+Li₂((1-c)/2)] - (1+c)[Li₂(-c)+Li₂((1+c)/2)]`. -/
noncomputable def M1four (α : ℝ) : ℝ :=
  π ^ 2 / 6 - lam ^ 2 - lam * M0 α
    - (1 - c α) * (D00.Li2 (c α) + D00.Li2 ((1 - c α) / 2))
    - (1 + c α) * (D00.Li2 (-(c α)) + D00.Li2 ((1 + c α) / 2))

/-- Closed form of `thm:M1`. -/
noncomputable def M1closed : ℝ :=
  (19 * Real.sqrt 2 / 48 - 1 / 24) * π ^ 2 + lam ^ 2 / 2 + 2 * Real.sqrt 2 * lam * L
    - (1 + 5 * Real.sqrt 2 / 2) * L ^ 2 - 5 * Real.sqrt 2 * D00.Li2 (Real.sqrt 2 - 1)

/-- The integrand of `eq:M1t`: `(log t / t)((1-t)/q - 1)`, `q = √(t²+(1-t)²)`. -/
noncomputable def M1tIntegrand (t : ℝ) : ℝ :=
  Real.log t / t * ((1 - t) / Real.sqrt (t ^ 2 + (1 - t) ^ 2) - 1)

/-- main.tex 1471-1473 with 882: at the right angle the corner `M_0(π/2)` is the squircle
`M_0 = λ - √2 L`. -/
def stmt_1 : Prop := M0 (π / 2) = lam - Real.sqrt 2 * L

/-- `prop:M1short` at `α = π/2` reproduces `thm:M1`: `M_1(π/2) = M_1`. -/
def stmt_2 : Prop := M1four (π / 2) = M1closed

/-- `thm:M1`: "Numerically `κ = 8M_1 = 2.02095276965635835695…`", for the value at `α = π/2`. -/
def stmt_3 : Prop := |8 * M1four (π / 2) - 2.02095276965635835695| < (1 : ℝ) / 10 ^ 20

/-- `eq:M1t` (the squircle's `M_1` as an integral) equals the four-dilogarithm value at the
right angle; the integrand is integrable on `(0,1)`. -/
def stmt_4 : Prop :=
  MeasureTheory.IntegrableOn M1tIntegrand (Set.Ioo 0 1) ∧
    ∫ t in (0 : ℝ)..1, M1tIntegrand t = M1four (π / 2)

end D03
-- END STATEMENTS

namespace D03
open Real


lemma s_sq : Real.sqrt 2 ^ 2 = 2 := Real.sq_sqrt (by norm_num)
lemma s_gt : (1.41 : ℝ) < Real.sqrt 2 := by
  rw [Real.lt_sqrt (by norm_num)]; norm_num
lemma s_lt : Real.sqrt 2 < 1.42 := by
  rw [Real.sqrt_lt' (by norm_num)]; norm_num
lemma s_pos : 0 < Real.sqrt 2 := by linarith [s_gt]

lemma log_s : Real.log (Real.sqrt 2) = lam / 2 := by
  rw [Real.log_sqrt (by norm_num)]; rfl

lemma ab : (Real.sqrt 2 - 1) * (1 + Real.sqrt 2) = 1 := by nlinarith [s_sq]

lemma log_b : Real.log (Real.sqrt 2 - 1) = -L := by
  have h := Real.log_mul (x := Real.sqrt 2 - 1) (y := 1 + Real.sqrt 2) (by linarith [s_gt]) (by linarith [s_gt])
  rw [ab, Real.log_one] at h
  unfold L; linarith

lemma L_pos : 0 < L := by
  unfold L; apply Real.log_pos; linarith [s_gt]

/-- K1 -/
lemma K1 : D00.Li2 (Real.sqrt 2 - 1) - D00.Li2 (-(Real.sqrt 2 - 1)) = π ^ 2 / 8 - L ^ 2 / 2 := by
  have hb0 : 0 < Real.sqrt 2 - 1 := by linarith [s_gt]
  have hb1 : Real.sqrt 2 - 1 < 1 := by linarith [s_lt]
  have h := D00.chi2_landen hb0 hb1
  have e1 : (1 - (Real.sqrt 2 - 1)) / (1 + (Real.sqrt 2 - 1)) = Real.sqrt 2 - 1 := by
    field_simp; nlinarith [s_sq]
  have e2 : (1 + (Real.sqrt 2 - 1)) / (1 - (Real.sqrt 2 - 1)) = 1 + Real.sqrt 2 := by
    have : (2 - Real.sqrt 2) ≠ 0 := by linarith [s_lt]
    have h2 : 1 - (Real.sqrt 2 - 1) = 2 - Real.sqrt 2 := by ring
    rw [h2, show 1 + (Real.sqrt 2 - 1) = Real.sqrt 2 by ring, div_eq_iff this]
    nlinarith [s_sq]
  rw [e1, e2, log_b] at h
  have hL : Real.log (1 + Real.sqrt 2) = L := rfl
  rw [hL] at h
  unfold D00.chi2 at h
  linarith


lemma log_a : Real.log (1 + Real.sqrt 2) = L := rfl

lemma K2 : D00.Li2 ((Real.sqrt 2 - 1) ^ 2) = 4 * D00.Li2 (Real.sqrt 2 - 1) - π ^ 2 / 4 + L ^ 2 := by
  have h := D00.duplication (x := Real.sqrt 2 - 1) (by
    rw [abs_le]; constructor <;> linarith [s_gt, s_lt])
  have := K1
  linarith

lemma K3 : D00.Li2 (2 - Real.sqrt 2) = -D00.Li2 (Real.sqrt 2 - 1) + π ^ 2 / 6 + (1/2) * lam * L - L ^ 2 := by
  have hb0 : 0 < Real.sqrt 2 - 1 := by linarith [s_gt]
  have hb1 : Real.sqrt 2 - 1 < 1 := by linarith [s_lt]
  have h := D00.reflection hb0 hb1
  have e : 1 - (Real.sqrt 2 - 1) = 2 - Real.sqrt 2 := by ring
  have e2 : 2 - Real.sqrt 2 = (Real.sqrt 2 - 1) * Real.sqrt 2 := by nlinarith [s_sq]
  have hl : Real.log (2 - Real.sqrt 2) = -L + lam / 2 := by
    rw [e2, Real.log_mul (by linarith [s_gt]) (by linarith [s_gt]), log_b, log_s]
  rw [e, hl, log_b] at h
  linarith

lemma K4 : D00.Li2 (-Real.sqrt 2) = D00.Li2 (Real.sqrt 2 - 1) - π ^ 2 / 6 - (1/2) * lam * L + (1/2) * L ^ 2 := by
  have h := D00.landen (x := -Real.sqrt 2) (by linarith [s_pos])
  have e : -Real.sqrt 2 / (-Real.sqrt 2 - 1) = 2 - Real.sqrt 2 := by
    have : (-Real.sqrt 2 - 1) ≠ 0 := by linarith [s_pos]
    rw [div_eq_iff this]; nlinarith [s_sq]
  rw [e, K3, show 1 - -Real.sqrt 2 = 1 + Real.sqrt 2 by ring, log_a] at h
  linarith

lemma K5 : D00.Li2 (-(Real.sqrt 2 - 1) ^ 2) = 6 * D00.Li2 (Real.sqrt 2 - 1) - 11 * π ^ 2 / 24 + 2 * L ^ 2 := by
  have h := D00.five_term (x := -Real.sqrt 2) (y := -(Real.sqrt 2 - 1)) (by linarith [s_pos])
    (by linarith [s_gt]) (by linarith [s_gt])
  have e1 : -Real.sqrt 2 / (1 - -(Real.sqrt 2 - 1)) = -1 := by
    have : (1 - -(Real.sqrt 2 - 1)) = Real.sqrt 2 := by ring
    rw [this]; field_simp
  have e2 : -(Real.sqrt 2 - 1) / (1 - -Real.sqrt 2) = -(Real.sqrt 2 - 1) ^ 2 := by
    have : (1 - -Real.sqrt 2) ≠ 0 := by linarith [s_pos]
    rw [div_eq_iff this]; nlinarith [s_sq]
  have e3 : -Real.sqrt 2 * -(Real.sqrt 2 - 1) / ((1 - -Real.sqrt 2) * (1 - -(Real.sqrt 2 - 1))) = (Real.sqrt 2 - 1) ^ 2 := by
    have h1 : (1 - -Real.sqrt 2) * (1 - -(Real.sqrt 2 - 1)) ≠ 0 := by
      have : 0 < (1 - -Real.sqrt 2) * (1 - -(Real.sqrt 2 - 1)) := by
        apply mul_pos <;> linarith [s_gt]
      exact this.ne'
    rw [div_eq_iff h1]; nlinarith [s_sq, s_gt]
  rw [e1, e2, e3, show 1 - -Real.sqrt 2 = 1 + Real.sqrt 2 by ring,
    show 1 - -(Real.sqrt 2 - 1) = Real.sqrt 2 by ring, log_a, log_s] at h
  have := K1
  have := K2
  have := K4
  rw [D00.Li2_neg_one] at h
  nlinarith


lemma log_two_s : Real.log (2 - Real.sqrt 2) = -L + lam / 2 := by
  have e2 : 2 - Real.sqrt 2 = (Real.sqrt 2 - 1) * Real.sqrt 2 := by nlinarith [s_sq]
  rw [e2, Real.log_mul (by linarith [s_gt]) (by linarith [s_gt]), log_b, log_s]

lemma log_four : Real.log 4 = 2 * lam := by
  rw [show (4 : ℝ) = 2 ^ 2 by norm_num, Real.log_pow]; unfold lam; push_cast; ring

lemma K6 : D00.Li2 ((2 + Real.sqrt 2) / 4) = 6 * D00.Li2 (Real.sqrt 2 - 1) - 7 * π ^ 2 / 24
    - (9/8) * lam ^ 2 - (3/2) * lam * L + (7/2) * L ^ 2 := by
  have h := D00.landen (x := (2 + Real.sqrt 2) / 4) (by linarith [s_lt])
  have e : (2 + Real.sqrt 2) / 4 / ((2 + Real.sqrt 2) / 4 - 1) = -(1 + Real.sqrt 2) ^ 2 := by
    have : ((2 + Real.sqrt 2) / 4 - 1) ≠ 0 := by linarith [s_lt]
    rw [div_eq_iff this]; nlinarith [s_sq]
  have hl : Real.log (1 - (2 + Real.sqrt 2) / 4) = -L - 3 * lam / 2 := by
    have : 1 - (2 + Real.sqrt 2) / 4 = (2 - Real.sqrt 2) / 4 := by ring
    rw [this, Real.log_div (by linarith [s_lt]) (by norm_num), log_two_s, log_four]; ring
  have hi := D00.inversion_neg (x := -(1 + Real.sqrt 2) ^ 2) (by nlinarith [s_pos])
  have e3 : 1 / -(1 + Real.sqrt 2) ^ 2 = -(Real.sqrt 2 - 1) ^ 2 := by
    have : (1 + Real.sqrt 2) ^ 2 ≠ 0 := by positivity
    field_simp; nlinarith [s_sq]
  have hl2 : Real.log (-(-(1 + Real.sqrt 2) ^ 2)) = 2 * L := by
    rw [neg_neg, Real.log_pow]; unfold L; push_cast; ring
  rw [e, hl] at h
  rw [e3, hl2, K5] at hi
  nlinarith

lemma K7a : D00.Li2 (-(Real.sqrt 2 / 2)) = -(1/2) * (L - lam / 2) ^ 2 - D00.Li2 (Real.sqrt 2 - 1) := by
  have h := D00.landen (x := -(Real.sqrt 2 / 2)) (by linarith [s_pos])
  have e : -(Real.sqrt 2 / 2) / (-(Real.sqrt 2 / 2) - 1) = Real.sqrt 2 - 1 := by
    have : (-(Real.sqrt 2 / 2) - 1) ≠ 0 := by linarith [s_pos]
    rw [div_eq_iff this]; nlinarith [s_sq]
  have hl : Real.log (1 - -(Real.sqrt 2 / 2)) = L - lam / 2 := by
    have : 1 - -(Real.sqrt 2 / 2) = (1 + Real.sqrt 2) * Real.sqrt 2 / 2 := by nlinarith [s_sq]
    rw [this, Real.log_div (by positivity) (by norm_num), Real.log_mul (by linarith [s_gt]) (by linarith [s_gt]), log_s, log_a]
    unfold lam; ring
  rw [e, hl] at h
  linarith

lemma K7b : D00.Li2 (Real.sqrt 2 / 2) = -(1/2) * (L + lam / 2) ^ 2 + π ^ 2 / 24 + L ^ 2 + D00.Li2 (Real.sqrt 2 - 1) := by
  have h := D00.landen (x := Real.sqrt 2 / 2) (by linarith [s_lt])
  have e : Real.sqrt 2 / 2 / (Real.sqrt 2 / 2 - 1) = -(1 + Real.sqrt 2) := by
    have : (Real.sqrt 2 / 2 - 1) ≠ 0 := by linarith [s_lt]
    rw [div_eq_iff this]; nlinarith [s_sq]
  have hl : Real.log (1 - Real.sqrt 2 / 2) = -L - lam / 2 := by
    have : 1 - Real.sqrt 2 / 2 = (Real.sqrt 2 - 1) * Real.sqrt 2 / 2 := by nlinarith [s_sq]
    rw [this, Real.log_div (mul_pos (by linarith [s_gt]) s_pos).ne' (by norm_num), Real.log_mul (by linarith [s_gt]) (by linarith [s_gt]), log_s, log_b]
    unfold lam; ring
  have hi := D00.inversion_neg (x := -(1 + Real.sqrt 2)) (by linarith [s_pos])
  have e3 : 1 / -(1 + Real.sqrt 2) = -(Real.sqrt 2 - 1) := by
    have : (1 + Real.sqrt 2) ≠ 0 := by positivity
    field_simp; nlinarith [s_sq]
  rw [neg_neg, log_a, e3] at hi
  rw [e, hl] at h
  have := K1
  nlinarith

lemma K8 : D00.Li2 ((2 - Real.sqrt 2) / 4) = π ^ 2 / 6 - (L - 3 * lam / 2) * (-L - 3 * lam / 2)
    - D00.Li2 ((2 + Real.sqrt 2) / 4) := by
  have h := D00.reflection (x := (2 + Real.sqrt 2) / 4) (by positivity) (by linarith [s_lt])
  have e : 1 - (2 + Real.sqrt 2) / 4 = (2 - Real.sqrt 2) / 4 := by ring
  have hl1 : Real.log (1 - (2 + Real.sqrt 2) / 4) = -L - 3 * lam / 2 := by
    rw [e, Real.log_div (by linarith [s_lt]) (by norm_num), log_two_s, log_four]; ring
  have hl2 : Real.log ((2 + Real.sqrt 2) / 4) = L - 3 * lam / 2 := by
    have : 2 + Real.sqrt 2 = (1 + Real.sqrt 2) * Real.sqrt 2 := by nlinarith [s_sq]
    rw [this, Real.log_div (by positivity) (by norm_num), Real.log_mul (by linarith [s_gt]) (by linarith [s_gt]), log_s, log_a, log_four]; ring
  rw [hl1, hl2, e] at h
  linarith



lemma hasSum_bounds {f : ℕ → ℝ} {S : ℝ} (hS : HasSum f S) (hf0 : ∀ n, 0 ≤ f n) (N : ℕ)
    (K r : ℝ) (hr0 : 0 ≤ r) (hr1 : r < 1) (hbd : ∀ n, f (n + N) ≤ K * r ^ n) :
    ∑ n ∈ Finset.range N, f n ≤ S ∧ S ≤ ∑ n ∈ Finset.range N, f n + K / (1 - r) := by
  have h1 := (hasSum_nat_add_iff' N).mpr hS
  have h2 : HasSum (fun n : ℕ => K * r ^ n) (K * (1 - r)⁻¹) :=
    (hasSum_geometric_of_lt_one hr0 hr1).mul_left K
  constructor
  · have := h1.nonneg (fun n => hf0 (n + N))
    linarith
  · have := hasSum_le hbd h1 h2
    rw [div_eq_mul_inv]
    linarith

lemma li2_bounds {x lo hi : ℝ} (h0 : 0 ≤ lo) (hlo : lo ≤ x) (hhi : x ≤ hi) (hhi1 : hi < 1) (N : ℕ) :
    ∑ n ∈ Finset.range N, lo ^ (n + 1) / ((n : ℝ) + 1) ^ 2 ≤ D00.Li2 x ∧
    D00.Li2 x ≤ ∑ n ∈ Finset.range N, hi ^ (n + 1) / ((n : ℝ) + 1) ^ 2
      + hi ^ (N + 1) / ((N : ℝ) + 1) ^ 2 / (1 - hi) := by
  have hx0 : 0 ≤ x := h0.trans hlo
  have hS := D00.hasSum_Li2 (x := x) (by rw [abs_le]; constructor <;> linarith)
  have hf0 : ∀ n : ℕ, 0 ≤ x ^ (n + 1) / ((n : ℝ) + 1) ^ 2 := fun n => by positivity
  have hb := hasSum_bounds hS hf0 N (hi ^ (N + 1) / ((N : ℝ) + 1) ^ 2) hi (h0.trans (hlo.trans hhi))
    hhi1 (fun n => by
      have e1 : x ^ (n + N + 1) ≤ hi ^ (n + N + 1) := pow_le_pow_left₀ hx0 hhi _
      have e2 : ((N : ℝ) + 1) ^ 2 ≤ (((n + N : ℕ) : ℝ) + 1) ^ 2 := by
        push_cast; apply pow_le_pow_left₀ (by positivity); linarith [(Nat.cast_nonneg n : (0:ℝ) ≤ n)]
      calc x ^ (n + N + 1) / (((n + N : ℕ) : ℝ) + 1) ^ 2
          ≤ hi ^ (n + N + 1) / ((N : ℝ) + 1) ^ 2 := by
            apply div_le_div₀ (pow_nonneg (h0.trans (hlo.trans hhi)) _) e1 (by positivity) e2
        _ = hi ^ (N + 1) / ((N : ℝ) + 1) ^ 2 * hi ^ n := by ring)
  refine ⟨?_, ?_⟩
  · refine le_trans ?_ hb.1
    apply Finset.sum_le_sum
    intro n _
    apply div_le_div_of_nonneg_right (pow_le_pow_left₀ h0 hlo _) (by positivity)
  · refine le_trans hb.2 ?_
    gcongr with n hn

lemma art_bounds {x lo hi : ℝ} (h0 : 0 ≤ lo) (hlo : lo ≤ x) (hhi : x ≤ hi) (hhi1 : hi < 1) (N : ℕ) :
    ∑ k ∈ Finset.range N, 2 * (1 / (2 * (k : ℝ) + 1)) * lo ^ (2 * k + 1) ≤ Real.log (1 + x) - Real.log (1 - x) ∧
    Real.log (1 + x) - Real.log (1 - x) ≤ ∑ k ∈ Finset.range N, 2 * (1 / (2 * (k : ℝ) + 1)) * hi ^ (2 * k + 1)
      + 2 * hi ^ (2 * N + 1) / (2 * (N : ℝ) + 1) / (1 - hi ^ 2) := by
  have hx0 : 0 ≤ x := h0.trans hlo
  have hS := Real.hasSum_log_sub_log_of_abs_lt_one (x := x) (by rw [abs_lt]; constructor <;> linarith)
  have hf0 : ∀ k : ℕ, 0 ≤ 2 * (1 / (2 * (k : ℝ) + 1)) * x ^ (2 * k + 1) := fun k => by positivity
  have hb := hasSum_bounds hS hf0 N (2 * hi ^ (2 * N + 1) / (2 * (N : ℝ) + 1)) (hi ^ 2) (by positivity)
    (by nlinarith [h0.trans (hlo.trans hhi)]) (fun n => by
      have e1 : x ^ (2 * (n + N) + 1) ≤ hi ^ (2 * (n + N) + 1) := pow_le_pow_left₀ hx0 hhi _
      have e2 : 1 / (2 * (((n + N : ℕ)) : ℝ) + 1) ≤ 1 / (2 * (N : ℝ) + 1) := by
        apply one_div_le_one_div_of_le (by positivity)
        push_cast; linarith [(Nat.cast_nonneg n : (0:ℝ) ≤ n)]
      have e3 : 2 * (1 / (2 * (((n + N : ℕ)) : ℝ) + 1)) * x ^ (2 * (n + N) + 1)
          ≤ 2 * (1 / (2 * (N : ℝ) + 1)) * hi ^ (2 * (n + N) + 1) := by
        apply mul_le_mul (by linarith) e1 (by positivity) (by positivity)
      calc _ ≤ _ := e3
        _ = 2 * hi ^ (2 * N + 1) / (2 * (N : ℝ) + 1) * (hi ^ 2) ^ n := by ring)
  refine ⟨?_, ?_⟩
  · refine le_trans ?_ hb.1
    apply Finset.sum_le_sum
    intro n _
    apply mul_le_mul_of_nonneg_left (pow_le_pow_left₀ h0 hlo _) (by positivity)
  · refine le_trans hb.2 ?_
    have : ∑ k ∈ Finset.range N, 2 * (1 / (2 * (k : ℝ) + 1)) * x ^ (2 * k + 1) ≤
        ∑ k ∈ Finset.range N, 2 * (1 / (2 * (k : ℝ) + 1)) * hi ^ (2 * k + 1) :=
      Finset.sum_le_sum (fun n _ => mul_le_mul_of_nonneg_left (pow_le_pow_left₀ hx0 hhi _) (by positivity))
    linarith

lemma log_1p : Real.log (1 + Real.sqrt 2 / 2) = L - lam / 2 := by
  have : 1 + Real.sqrt 2 / 2 = (1 + Real.sqrt 2) * Real.sqrt 2 / 2 := by nlinarith [s_sq]
  rw [this, Real.log_div (by positivity) (by norm_num), Real.log_mul (by linarith [s_gt]) (by linarith [s_gt]), log_s, log_a]
  unfold lam; ring

lemma log_1m : Real.log (1 - Real.sqrt 2 / 2) = -L - lam / 2 := by
  have : 1 - Real.sqrt 2 / 2 = (Real.sqrt 2 - 1) * Real.sqrt 2 / 2 := by nlinarith [s_sq]
  rw [this, Real.log_div (mul_pos (by linarith [s_gt]) s_pos).ne' (by norm_num), Real.log_mul (by linarith [s_gt]) (by linarith [s_gt]), log_s, log_b]
  unfold lam; ring

lemma c_half_pi : c (π / 2) = Real.sqrt 2 / 2 := by
  unfold c
  rw [show π / 2 / 2 = π / 4 by ring, Real.cos_pi_div_four]

lemma proof1 : stmt_1 := by
  unfold stmt_1 M0
  rw [c_half_pi, log_1p, log_1m]
  ring

lemma proof2 : stmt_2 := by
  unfold stmt_2 M1four M1closed M0
  rw [c_half_pi, log_1p, log_1m, K7a, K7b,
    show (1 - Real.sqrt 2 / 2) / 2 = (2 - Real.sqrt 2) / 4 by ring,
    show (1 + Real.sqrt 2 / 2) / 2 = (2 + Real.sqrt 2) / 4 by ring, K8, K6]
  ring

lemma s_bounds : (1.41421356237309504880168872 : ℝ) ≤ Real.sqrt 2 ∧ Real.sqrt 2 ≤ (1.41421356237309504880168873 : ℝ) := by
  constructor
  · apply le_of_lt; rw [Real.lt_sqrt (by norm_num)]; norm_num
  · apply le_of_lt; rw [Real.sqrt_lt' (by norm_num)]; norm_num

lemma u_bounds : (0.46753399702300461762696479 : ℝ) ≤ D00.Li2 (Real.sqrt 2 - 1) ∧ D00.Li2 (Real.sqrt 2 - 1) ≤ (0.46753399702300461762696481 : ℝ) := by
  have hs := s_bounds
  have h := li2_bounds (x := Real.sqrt 2 - 1) (lo := (0.41421356237309504880168872 : ℝ)) (hi := (0.41421356237309504880168873 : ℝ)) (by norm_num)
    (by linarith [hs.1]) (by linarith [hs.2]) (by norm_num) 60
  constructor
  · refine le_trans ?_ h.1
    norm_num [Finset.sum_range_succ]
  · refine le_trans h.2 ?_
    norm_num [Finset.sum_range_succ]

lemma w_bounds : (0.58224052646501250590265631 : ℝ) ≤ D00.Li2 (1/2) ∧ D00.Li2 (1/2) ≤ (0.58224052646501250590265633 : ℝ) := by
  have h := li2_bounds (x := 1/2) (lo := 1/2) (hi := 1/2) (by norm_num) le_rfl le_rfl (by norm_num) 75
  constructor
  · refine le_trans ?_ h.1
    norm_num [Finset.sum_range_succ]
  · refine le_trans h.2 ?_
    norm_num [Finset.sum_range_succ]

lemma lam_bounds : (0.69314718055994530941723212 : ℝ) ≤ lam ∧ lam ≤ (0.69314718055994530941723213 : ℝ) := by
  have h := art_bounds (x := 1/3) (lo := 1/3) (hi := 1/3) (by norm_num) le_rfl le_rfl (by norm_num) 28
  have e : Real.log (1 + 1/3) - Real.log (1 - 1/3) = lam := by
    rw [← Real.log_div (by norm_num) (by norm_num)]; unfold lam; norm_num
  rw [e] at h
  constructor
  · refine le_trans ?_ h.1
    norm_num [Finset.sum_range_succ]
  · refine le_trans h.2 ?_
    norm_num [Finset.sum_range_succ]

lemma L_bounds : (0.88137358701954302523260931 : ℝ) ≤ L ∧ L ≤ (0.88137358701954302523260934 : ℝ) := by
  have hs := s_bounds
  have h := art_bounds (x := Real.sqrt 2 - 1) (lo := (0.41421356237309504880168872 : ℝ)) (hi := (0.41421356237309504880168873 : ℝ)) (by norm_num)
    (by linarith [hs.1]) (by linarith [hs.2]) (by norm_num) 36
  have e : Real.log (1 + (Real.sqrt 2 - 1)) - Real.log (1 - (Real.sqrt 2 - 1)) = L := by
    rw [show 1 + (Real.sqrt 2 - 1) = Real.sqrt 2 by ring, show 1 - (Real.sqrt 2 - 1) = 2 - Real.sqrt 2 by ring,
      log_s, log_two_s]; ring
  rw [e] at h
  constructor
  · refine le_trans ?_ h.1
    norm_num [Finset.sum_range_succ]
  · refine le_trans h.2 ?_
    norm_num [Finset.sum_range_succ]

lemma M_abs (s u w l Lg : ℝ) (hs : (1.41421356237309504880168872 : ℝ) ≤ s ∧ s ≤ (1.41421356237309504880168873 : ℝ))
    (hu : (0.46753399702300461762696479 : ℝ) ≤ u ∧ u ≤ (0.46753399702300461762696481 : ℝ)) (hw : (0.58224052646501250590265631 : ℝ) ≤ w ∧ w ≤ (0.58224052646501250590265633 : ℝ))
    (hl : (0.69314718055994530941723212 : ℝ) ≤ l ∧ l ≤ (0.69314718055994530941723213 : ℝ)) (hL : (0.88137358701954302523260931 : ℝ) ≤ Lg ∧ Lg ≤ (0.88137358701954302523260934 : ℝ)) :
    |8 * ((19/4) * s * w - w / 2 + (19/8 * s + 1/4) * l ^ 2 + 2 * s * l * Lg - Lg ^ 2
      - (5/2) * s * Lg ^ 2 - 5 * s * u) - 2.02095276965635835695| < (1 : ℝ) / 10 ^ 20 := by
  obtain ⟨hs1, hs2⟩ := hs
  obtain ⟨hu1, hu2⟩ := hu
  obtain ⟨hw1, hw2⟩ := hw
  obtain ⟨hl1, hl2⟩ := hl
  obtain ⟨hL1, hL2⟩ := hL
  have hs0 : 0 ≤ s := by linarith
  have hu0 : 0 ≤ u := by linarith
  have hw0 : 0 ≤ w := by linarith
  have hl0 : 0 ≤ l := by linarith
  have hL0 : 0 ≤ Lg := by linarith
  have p1 : s * w ≤ (1.41421356237309504880168873 : ℝ) * (0.58224052646501250590265633 : ℝ) := mul_le_mul hs2 hw2 hw0 (by norm_num)
  have p1' : (1.41421356237309504880168872 : ℝ) * (0.58224052646501250590265631 : ℝ) ≤ s * w := mul_le_mul hs1 hw1 (by norm_num) hs0
  have q : l ^ 2 ≤ (0.69314718055994530941723213 : ℝ) ^ 2 := pow_le_pow_left₀ hl0 hl2 2
  have q' : (0.69314718055994530941723212 : ℝ) ^ 2 ≤ l ^ 2 := pow_le_pow_left₀ (by norm_num) hl1 2
  have p2 : s * l ^ 2 ≤ (1.41421356237309504880168873 : ℝ) * (0.69314718055994530941723213 : ℝ) ^ 2 := mul_le_mul hs2 q (by positivity) (by norm_num)
  have p2' : (1.41421356237309504880168872 : ℝ) * (0.69314718055994530941723212 : ℝ) ^ 2 ≤ s * l ^ 2 := mul_le_mul hs1 q' (by norm_num) hs0
  have r1 : s * l ≤ (1.41421356237309504880168873 : ℝ) * (0.69314718055994530941723213 : ℝ) := mul_le_mul hs2 hl2 hl0 (by norm_num)
  have r1' : (1.41421356237309504880168872 : ℝ) * (0.69314718055994530941723212 : ℝ) ≤ s * l := mul_le_mul hs1 hl1 (by norm_num) hs0
  have p3 : s * l * Lg ≤ (1.41421356237309504880168873 : ℝ) * (0.69314718055994530941723213 : ℝ) * (0.88137358701954302523260934 : ℝ) := mul_le_mul r1 hL2 hL0 (by norm_num)
  have p3' : (1.41421356237309504880168872 : ℝ) * (0.69314718055994530941723212 : ℝ) * (0.88137358701954302523260931 : ℝ) ≤ s * l * Lg := mul_le_mul r1' hL1 (by norm_num) (by positivity)
  have t : Lg ^ 2 ≤ (0.88137358701954302523260934 : ℝ) ^ 2 := pow_le_pow_left₀ hL0 hL2 2
  have t' : (0.88137358701954302523260931 : ℝ) ^ 2 ≤ Lg ^ 2 := pow_le_pow_left₀ (by norm_num) hL1 2
  have p4 : s * Lg ^ 2 ≤ (1.41421356237309504880168873 : ℝ) * (0.88137358701954302523260934 : ℝ) ^ 2 := mul_le_mul hs2 t (by positivity) (by norm_num)
  have p4' : (1.41421356237309504880168872 : ℝ) * (0.88137358701954302523260931 : ℝ) ^ 2 ≤ s * Lg ^ 2 := mul_le_mul hs1 t' (by norm_num) hs0
  have p5 : s * u ≤ (1.41421356237309504880168873 : ℝ) * (0.46753399702300461762696481 : ℝ) := mul_le_mul hs2 hu2 hu0 (by norm_num)
  have p5' : (1.41421356237309504880168872 : ℝ) * (0.46753399702300461762696479 : ℝ) ≤ s * u := mul_le_mul hs1 hu1 (by norm_num) hs0
  rw [abs_lt]
  constructor
  · nlinarith
  · nlinarith

lemma M1closed_expand : M1closed = (19/4) * Real.sqrt 2 * D00.Li2 (1/2) - D00.Li2 (1/2) / 2
      + (19/8 * Real.sqrt 2 + 1/4) * lam ^ 2 + 2 * Real.sqrt 2 * lam * L - L ^ 2
      - (5/2) * Real.sqrt 2 * L ^ 2 - 5 * Real.sqrt 2 * D00.Li2 (Real.sqrt 2 - 1) := by
  have hpi : π ^ 2 = 12 * D00.Li2 (1/2) + 6 * lam ^ 2 := by
    have := D00.Li2_half
    unfold lam; linarith
  unfold M1closed; rw [hpi]; ring

lemma proof3 : stmt_3 := by
  unfold stmt_3
  rw [proof2, M1closed_expand]
  exact M_abs _ _ _ _ _ s_bounds u_bounds w_bounds lam_bounds L_bounds
/-- the Euler substitution -/
noncomputable def tt (v : ℝ) : ℝ := (v ^ 2 + 2 * v - 1) / (4 * v)

lemma tt_eq (v : ℝ) (hv : v ≠ 0) : tt v = (v - (Real.sqrt 2 - 1)) * (v + (1 + Real.sqrt 2)) / (4 * v) := by
  unfold tt; congr 1; nlinarith [s_sq]

lemma sqrt_tt {v : ℝ} (hv : 0 < v) :
    Real.sqrt (tt v ^ 2 + (1 - tt v) ^ 2) = Real.sqrt 2 * (v ^ 2 + 1) / (4 * v) := by
  have h : tt v ^ 2 + (1 - tt v) ^ 2 = (Real.sqrt 2 * (v ^ 2 + 1) / (4 * v)) ^ 2 := by
    unfold tt
    field_simp
    nlinarith [s_sq]
  rw [h, Real.sqrt_sq (by have := s_pos; positivity)]

lemma key {v : ℝ} (hv : 0 < v) (hb : Real.sqrt 2 - 1 < v) :
    1 / tt v * ((1 - tt v) / Real.sqrt (tt v ^ 2 + (1 - tt v) ^ 2) - 1) * ((v ^ 2 + 1) / (4 * v ^ 2))
      = (1 - Real.sqrt 2 / 2) / v - 2 / (v + (1 + Real.sqrt 2)) := by
  rw [sqrt_tt hv, tt_eq v hv.ne']
  have h1 : v - (Real.sqrt 2 - 1) ≠ 0 := by linarith
  have h2 : v + (1 + Real.sqrt 2) ≠ 0 := by have := s_pos; linarith
  have h3 : Real.sqrt 2 ≠ 0 := s_pos.ne'
  have h4 : v ^ 2 + 1 ≠ 0 := by positivity
  field_simp
  linear_combination (v ^ 2 - 2 * v - Real.sqrt 2 ^ 2 + 2 * Real.sqrt 2 + 1) * s_sq

lemma h_bound {t : ℝ} (h0 : 0 ≤ t) (h1 : t ≤ 1) :
    |(1 - t) / Real.sqrt (t ^ 2 + (1 - t) ^ 2) - 1| ≤ 2 * t ^ 2 := by
  set q := Real.sqrt (t ^ 2 + (1 - t) ^ 2) with hq
  have hq2 : q ^ 2 = t ^ 2 + (1 - t) ^ 2 := Real.sq_sqrt (by positivity)
  have hqge : 1 / 2 ≤ q ^ 2 := by nlinarith [sq_nonneg (t - 1/2)]
  have hq0 : 0 < q := by
    rw [hq]; apply Real.sqrt_pos.mpr; nlinarith [sq_nonneg (t - 1/2)]
  have hq1 : 1 - t ≤ q := by
    rw [hq]; apply Real.le_sqrt_of_sq_le; nlinarith
  have hd : (q - (1 - t)) * (q + (1 - t)) = t ^ 2 := by nlinarith
  have e : (1 - t) / q - 1 = -((q - (1 - t)) / q) := by field_simp; ring
  rw [e, abs_neg, abs_of_nonneg (div_nonneg (by linarith) hq0.le), div_le_iff₀ hq0]
  have hd0 : 0 ≤ q - (1 - t) := by linarith
  have : (q - (1 - t)) * q ≤ t ^ 2 := by nlinarith
  nlinarith [mul_nonneg hd0 hq0.le, sq_nonneg t]

lemma g_bound {t : ℝ} (h0 : 0 < t) (h1 : t ≤ 1) : |M1tIntegrand t| ≤ 2 * |t * Real.log t| := by
  unfold M1tIntegrand
  rw [abs_mul, abs_div, abs_of_pos h0, abs_mul, abs_of_pos h0]
  have := h_bound h0.le h1
  calc |Real.log t| / t * |(1 - t) / Real.sqrt (t ^ 2 + (1 - t) ^ 2) - 1|
      ≤ |Real.log t| / t * (2 * t ^ 2) := by gcongr
    _ = 2 * (t * |Real.log t|) := by field_simp

lemma g_zero : M1tIntegrand 0 = 0 := by simp [M1tIntegrand]

lemma g_cont : ContinuousOn M1tIntegrand (Set.Icc 0 1) := by
  intro t ht
  rcases eq_or_lt_of_le ht.1 with h | h
  · subst h
    rw [ContinuousWithinAt, g_zero]
    have hc : Filter.Tendsto (fun t : ℝ => 2 * |t * Real.log t|) (nhdsWithin 0 (Set.Icc 0 1)) (nhds 0) := by
      have := ((Real.continuous_mul_log.abs).const_mul 2).tendsto 0
      simpa using this.mono_left nhdsWithin_le_nhds
    refine squeeze_zero_norm' ?_ hc
    filter_upwards [self_mem_nhdsWithin] with t ht'
    rcases eq_or_lt_of_le ht'.1 with h' | h'
    · subst h'; simp [g_zero]
    · rw [Real.norm_eq_abs]; exact g_bound h' ht'.2
  · have hq : (t ^ 2 + (1 - t) ^ 2) ≠ 0 :=
      ne_of_gt (by nlinarith [pow_pos h 2, sq_nonneg (1 - t)])
    have hsq : Real.sqrt (t ^ 2 + (1 - t) ^ 2) ≠ 0 := by
      rw [Ne, Real.sqrt_eq_zero (by positivity)]; exact hq
    have hsc : ContinuousAt (fun t : ℝ => Real.sqrt (t ^ 2 + (1 - t) ^ 2)) t :=
      Real.continuous_sqrt.continuousAt.comp (by fun_prop)
    have h1 : ContinuousAt (fun t : ℝ => (1 - t) / Real.sqrt (t ^ 2 + (1 - t) ^ 2) - 1) t :=
      (((continuousAt_const.sub continuousAt_id).div hsc hsq)).sub continuousAt_const
    have h2 : ContinuousAt (fun t : ℝ => Real.log t / t) t :=
      (Real.continuousAt_log h.ne').div continuousAt_id h.ne'
    exact (h2.mul h1).continuousWithinAt


noncomputable def PsiG (x y v : ℝ) : ℝ :=
  Real.log (x - y) * Real.log ((v - y) / (x - y)) - D00.Li2 ((v - y) / (x - y))

lemma hasDerivAt_PsiG {x y v : ℝ} (hxy : x ≠ y) (hv : v ≠ y) (hvx : v ≠ x) :
    HasDerivAt (PsiG x y) (Real.log (v - x) / (v - y)) v := by
  have hxy' : x - y ≠ 0 := sub_ne_zero.mpr hxy
  have hv' : v - y ≠ 0 := sub_ne_zero.mpr hv
  have hw : HasDerivAt (fun v => (v - y) / (x - y)) (1 / (x - y)) v := by
    have := ((hasDerivAt_id v).sub_const y).div_const (x - y)
    simpa using this
  have hw0 : (v - y) / (x - y) ≠ 0 := div_ne_zero hv' hxy'
  have hw1 : (v - y) / (x - y) ≠ 1 := by
    intro h
    apply hvx
    rw [div_eq_one_iff_eq hxy'] at h
    linarith
  have hlog := hw.log hw0
  have hli := D00.Li2_comp hw hw1 hw0
  have h := (hlog.const_mul (Real.log (x - y))).sub hli
  have hx : Real.log (v - x) = Real.log (x - y) + Real.log (1 - (v - y) / (x - y)) := by
    have : v - x = -((x - y) * (1 - (v - y) / (x - y))) := by field_simp; ring
    rw [this, Real.log_neg_eq_log, Real.log_mul hxy' ?_]
    intro h0
    apply hw1
    linarith
  unfold PsiG
  convert h using 1
  rw [hx]
  field_simp
  ring

lemma hasDerivAt_logsq {x v : ℝ} (hv : v ≠ x) :
    HasDerivAt (fun v => Real.log (v - x) ^ 2 / 2) (Real.log (v - x) / (v - x)) v := by
  have hv' : v - x ≠ 0 := sub_ne_zero.mpr hv
  have h1 : HasDerivAt (fun v => v - x) 1 v := (hasDerivAt_id v).sub_const x
  have h2 := h1.log hv'
  have h3 := (h2.pow 2).div_const 2
  convert h3 using 1
  field_simp
  simp
  ring


noncomputable def Phi (v : ℝ) : ℝ :=
  (1 - Real.sqrt 2 / 2) * (PsiG (Real.sqrt 2 - 1) 0 v + PsiG (-(1 + Real.sqrt 2)) 0 v
      - Real.log (v - 0) ^ 2 / 2)
    - 2 * (PsiG (Real.sqrt 2 - 1) (-(1 + Real.sqrt 2)) v + Real.log (v - -(1 + Real.sqrt 2)) ^ 2 / 2
      - PsiG 0 (-(1 + Real.sqrt 2)) v)
    - 2 * lam * ((1 - Real.sqrt 2 / 2) * Real.log v - 2 * Real.log (v + (1 + Real.sqrt 2)))

lemma log_tt {v : ℝ} (hv : 0 < v) (hb : Real.sqrt 2 - 1 < v) :
    Real.log (tt v) = Real.log (v - (Real.sqrt 2 - 1)) + Real.log (v + (1 + Real.sqrt 2))
      - Real.log v - 2 * lam := by
  have h1 : v - (Real.sqrt 2 - 1) ≠ 0 := by linarith
  have h2 : v + (1 + Real.sqrt 2) ≠ 0 := by have := s_pos; linarith
  rw [tt_eq v hv.ne', Real.log_div (mul_ne_zero h1 h2) (by positivity), Real.log_mul h1 h2,
    Real.log_mul (by norm_num) hv.ne', show (4 : ℝ) = 4 from rfl, log_four]
  ring

lemma hasDerivAt_Phi {v : ℝ} (hb : Real.sqrt 2 - 1 < v) :
    HasDerivAt Phi (M1tIntegrand (tt v) * ((v ^ 2 + 1) / (4 * v ^ 2))) v := by
  have hs := s_pos
  have hs1 := s_gt
  have hv : 0 < v := by linarith
  have hv2 : v ≠ Real.sqrt 2 - 1 := hb.ne'
  have hva : v ≠ -(1 + Real.sqrt 2) := by intro h; linarith
  have d1 := hasDerivAt_PsiG (x := Real.sqrt 2 - 1) (y := 0) (v := v) (by linarith) hv.ne' hv2
  have d2 := hasDerivAt_PsiG (x := -(1 + Real.sqrt 2)) (y := 0) (v := v) (by linarith) hv.ne' hva
  have d3 := hasDerivAt_logsq (x := 0) (v := v) hv.ne'
  have d4 := hasDerivAt_PsiG (x := Real.sqrt 2 - 1) (y := -(1 + Real.sqrt 2)) (v := v) (by linarith) hva hv2
  have d5 := hasDerivAt_logsq (x := -(1 + Real.sqrt 2)) (v := v) hva
  have d6 := hasDerivAt_PsiG (x := 0) (y := -(1 + Real.sqrt 2)) (v := v) (by linarith) hva hv.ne'
  have d7 : HasDerivAt (fun v => Real.log v) (1 / v) v := by
    simpa using (hasDerivAt_id v).log hv.ne'
  have d8 : HasDerivAt (fun v => Real.log (v + (1 + Real.sqrt 2))) (1 / (v + (1 + Real.sqrt 2))) v := by
    have := ((hasDerivAt_id v).add_const (1 + Real.sqrt 2)).log (by simp only [id]; exact (by linarith : 0 < v + (1 + Real.sqrt 2)).ne')
    simpa using this
  have hsum := (((d1.add d2).sub d3).const_mul (1 - Real.sqrt 2 / 2)).sub
    ((((d4.add d5).sub d6)).const_mul 2) |>.sub
    (((d7.const_mul (1 - Real.sqrt 2 / 2)).sub (d8.const_mul 2)).const_mul (2 * lam))
  have hfun : Phi = _ := rfl
  unfold Phi
  convert hsum using 1
  have hk := key hv hb
  have hlt := log_tt hv hb
  unfold M1tIntegrand
  have e : Real.log (tt v) / tt v * ((1 - tt v) / Real.sqrt (tt v ^ 2 + (1 - tt v) ^ 2) - 1) *
      ((v ^ 2 + 1) / (4 * v ^ 2)) = Real.log (tt v) *
      (1 / tt v * ((1 - tt v) / Real.sqrt (tt v ^ 2 + (1 - tt v) ^ 2) - 1) * ((v ^ 2 + 1) / (4 * v ^ 2))) := by
    ring
  rw [e, hk, hlt]
  have h2 : v - -(1 + Real.sqrt 2) = v + (1 + Real.sqrt 2) := by ring
  simp only [h2, sub_zero]
  field_simp

lemma log_2s : Real.log (2 * Real.sqrt 2) = 3 * lam / 2 := by
  rw [Real.log_mul (by norm_num) s_pos.ne', log_s]; unfold lam; ring

lemma log_aa_s : Real.log ((1 + Real.sqrt 2) / Real.sqrt 2) = L - lam / 2 := by
  rw [Real.log_div (by have := s_pos; positivity) s_pos.ne', log_s, log_a]

lemma log_2aa : Real.log (2 * (1 + Real.sqrt 2)) = lam + L := by
  rw [Real.log_mul (by norm_num) (by have := s_pos; positivity), log_a]; unfold lam; ring

lemma log_aa_sq : Real.log ((1 + Real.sqrt 2) ^ 2) = 2 * L := by
  rw [Real.log_pow]; push_cast; rw [log_a]

lemma log_bb_sq : Real.log ((Real.sqrt 2 - 1) ^ 2) = -(2 * L) := by
  rw [Real.log_pow, log_b]; push_cast; ring

lemma W1 : PsiG (Real.sqrt 2 - 1) 0 (1 + Real.sqrt 2) - PsiG (Real.sqrt 2 - 1) 0 (Real.sqrt 2 - 1)
    = D00.Li2 ((Real.sqrt 2 - 1) ^ 2) - π ^ 2 / 6 := by
  have hb : 0 < Real.sqrt 2 - 1 := by linarith [s_gt]
  have e1 : ((1 + Real.sqrt 2) - 0) / ((Real.sqrt 2 - 1) - 0) = (1 + Real.sqrt 2) ^ 2 := by
    rw [sub_zero, sub_zero, div_eq_iff hb.ne']; nlinarith [s_sq]
  have e2 : ((Real.sqrt 2 - 1) - 0) / ((Real.sqrt 2 - 1) - 0) = 1 := by
    rw [sub_zero]; exact div_self hb.ne'
  have hinv := D00.inversion_gt_one (x := (1 + Real.sqrt 2) ^ 2) (by nlinarith [s_gt])
  have e3 : 1 / ((1 + Real.sqrt 2) ^ 2) = (Real.sqrt 2 - 1) ^ 2 := by
    have : (1 + Real.sqrt 2) ^ 2 ≠ 0 := by have := s_pos; positivity
    field_simp; nlinarith [s_sq]
  rw [e3, log_aa_sq] at hinv
  unfold PsiG
  rw [e1, e2, sub_zero, log_b, log_aa_sq, Real.log_one, D00.Li2_one]
  nlinarith [hinv]

lemma W2 : PsiG (-(1 + Real.sqrt 2)) 0 (1 + Real.sqrt 2) - PsiG (-(1 + Real.sqrt 2)) 0 (Real.sqrt 2 - 1)
    = π ^ 2 / 12 + 2 * L ^ 2 + D00.Li2 (-(Real.sqrt 2 - 1) ^ 2) := by
  have ha : 0 < 1 + Real.sqrt 2 := by linarith [s_pos]
  have e1 : ((1 + Real.sqrt 2) - 0) / (-(1 + Real.sqrt 2) - 0) = -1 := by
    rw [sub_zero, sub_zero]; field_simp
  have e2 : ((Real.sqrt 2 - 1) - 0) / (-(1 + Real.sqrt 2) - 0) = -(Real.sqrt 2 - 1) ^ 2 := by
    rw [sub_zero, sub_zero, div_eq_iff (by linarith [s_pos])]; nlinarith [s_sq]
  have hn : Real.log (-(Real.sqrt 2 - 1) ^ 2) = -(2 * L) := by
    rw [Real.log_neg_eq_log, log_bb_sq]
  unfold PsiG
  rw [e1, e2, sub_zero, Real.log_neg_eq_log, log_a, hn, Real.log_neg_eq_log, Real.log_one,
    D00.Li2_neg_one]
  ring

lemma W4 : PsiG (Real.sqrt 2 - 1) (-(1 + Real.sqrt 2)) (1 + Real.sqrt 2)
      - PsiG (Real.sqrt 2 - 1) (-(1 + Real.sqrt 2)) (Real.sqrt 2 - 1)
    = D00.Li2 (2 - Real.sqrt 2) - π ^ 2 / 6 + lam * L + L ^ 2 / 2 - 5 * lam ^ 2 / 8 := by
  have hs := s_pos
  have e0 : (Real.sqrt 2 - 1) - -(1 + Real.sqrt 2) = 2 * Real.sqrt 2 := by ring
  have e1 : (1 + Real.sqrt 2) - -(1 + Real.sqrt 2) = 2 * (1 + Real.sqrt 2) := by ring
  have e2 : (Real.sqrt 2 - 1) - -(1 + Real.sqrt 2) = 2 * Real.sqrt 2 := e0
  have r1 : 2 * (1 + Real.sqrt 2) / (2 * Real.sqrt 2) = (2 + Real.sqrt 2) / 2 := by
    rw [div_eq_div_iff (by positivity) (by norm_num)]; nlinarith [s_sq]
  have r2 : 2 * Real.sqrt 2 / (2 * Real.sqrt 2) = 1 := div_self (by positivity)
  have hr1 : 2 * (1 + Real.sqrt 2) / (2 * Real.sqrt 2) = (1 + Real.sqrt 2) / Real.sqrt 2 := by
    field_simp
  have hinv := D00.inversion_gt_one (x := (2 + Real.sqrt 2) / 2) (by linarith [s_gt])
  have e3 : 1 / ((2 + Real.sqrt 2) / 2) = 2 - Real.sqrt 2 := by
    have : (2 + Real.sqrt 2) ≠ 0 := by positivity
    field_simp; nlinarith [s_sq]
  have hl : Real.log ((2 + Real.sqrt 2) / 2) = L - lam / 2 := by
    rw [← r1, hr1, log_aa_s]
  rw [e3, hl] at hinv
  unfold PsiG
  rw [e1, e2, r1, r2, hl, Real.log_one, log_2s, D00.Li2_one]
  nlinarith [hinv]

lemma W5 : Real.log (((1 + Real.sqrt 2)) - -(1 + Real.sqrt 2)) ^ 2 / 2
      - Real.log ((Real.sqrt 2 - 1) - -(1 + Real.sqrt 2)) ^ 2 / 2
    = (lam + L) ^ 2 / 2 - (3 * lam / 2) ^ 2 / 2 := by
  have e1 : (1 + Real.sqrt 2) - -(1 + Real.sqrt 2) = 2 * (1 + Real.sqrt 2) := by ring
  have e2 : (Real.sqrt 2 - 1) - -(1 + Real.sqrt 2) = 2 * Real.sqrt 2 := by ring
  rw [e1, e2, log_2aa, log_2s]

lemma W6 : PsiG 0 (-(1 + Real.sqrt 2)) (1 + Real.sqrt 2) - PsiG 0 (-(1 + Real.sqrt 2)) (Real.sqrt 2 - 1)
    = π ^ 2 / 12 + lam * L + L ^ 2 / 2 - 9 * lam ^ 2 / 8 - D00.Li2 ((2 + Real.sqrt 2) / 4) := by
  have hs := s_pos
  have ha : 0 < 1 + Real.sqrt 2 := by linarith
  have e0 : (0 : ℝ) - -(1 + Real.sqrt 2) = 1 + Real.sqrt 2 := by ring
  have e1 : (1 + Real.sqrt 2) - -(1 + Real.sqrt 2) = 2 * (1 + Real.sqrt 2) := by ring
  have e2 : (Real.sqrt 2 - 1) - -(1 + Real.sqrt 2) = 2 * Real.sqrt 2 := by ring
  have r1 : 2 * (1 + Real.sqrt 2) / (1 + Real.sqrt 2) = 2 := by field_simp
  have r2 : 2 * Real.sqrt 2 / (1 + Real.sqrt 2) = 4 - 2 * Real.sqrt 2 := by
    rw [div_eq_iff ha.ne']; nlinarith [s_sq]
  have hl2 : Real.log (4 - 2 * Real.sqrt 2) = 3 * lam / 2 - L := by
    have : 4 - 2 * Real.sqrt 2 = 2 * (2 - Real.sqrt 2) := by ring
    rw [this, Real.log_mul (by norm_num) (by linarith [s_lt]), log_two_s]
    unfold lam; ring
  have hi2 := D00.inversion_gt_one (x := 2) (by norm_num)
  have hi4 := D00.inversion_gt_one (x := 4 - 2 * Real.sqrt 2) (by linarith [s_lt])
  have e3 : 1 / (4 - 2 * Real.sqrt 2) = (2 + Real.sqrt 2) / 4 := by
    have : (4 - 2 * Real.sqrt 2) ≠ 0 := by linarith [s_lt]
    field_simp; nlinarith [s_sq]
  rw [e3, hl2] at hi4
  have hl' : Real.log 2 = lam := rfl
  rw [D00.Li2_half, hl'] at hi2
  unfold PsiG
  rw [e0, e1, e2, r1, r2, log_a, hl2, hl']
  nlinarith [hi2, hi4]

lemma Phi_diff : Phi (1 + Real.sqrt 2) - Phi (Real.sqrt 2 - 1) = M1closed := by
  have w1 := W1
  have w2 := W2
  have w4 := W4
  have w5 := W5
  have w6 := W6
  have k2 := K2
  have k3 := K3
  have k5 := K5
  have k6 := K6
  have e1 : 1 + Real.sqrt 2 + (1 + Real.sqrt 2) = 2 * (1 + Real.sqrt 2) := by ring
  have e2 : Real.sqrt 2 - 1 + (1 + Real.sqrt 2) = 2 * Real.sqrt 2 := by ring
  unfold Phi M1closed
  rw [e1, e2, log_2aa, log_2s, log_a, log_b]
  simp only [sub_zero] at w5 ⊢
  rw [log_a, log_b] at *
  linear_combination (1 - Real.sqrt 2 / 2) * w1 + (1 - Real.sqrt 2 / 2) * w2 - 2 * w4 - 2 * w5
    + 2 * w6 + (1 - Real.sqrt 2 / 2) * k2 + (1 - Real.sqrt 2 / 2) * k5 - 2 * k3 - 2 * k6

lemma continuousAt_PsiG {x y v : ℝ} (hxy : x ≠ y) (hv : v ≠ y) : ContinuousAt (PsiG x y) v := by
  have hxy' : x - y ≠ 0 := sub_ne_zero.mpr hxy
  have hw : ContinuousAt (fun v : ℝ => (v - y) / (x - y)) v :=
    ((continuousAt_id.sub continuousAt_const).div_const _)
  have hw0 : (v - y) / (x - y) ≠ 0 := div_ne_zero (sub_ne_zero.mpr hv) hxy'
  unfold PsiG
  exact (continuousAt_const.mul (hw.log hw0)).sub (D00.continuous_Li2.continuousAt.comp (f := fun v : ℝ => (v - y) / (x - y)) hw)

lemma continuousAt_logsq {x v : ℝ} (hv : v ≠ x) :
    ContinuousAt (fun v : ℝ => Real.log (v - x) ^ 2 / 2) v := by
  have h : ContinuousAt (fun v : ℝ => Real.log (v - x)) v :=
    (continuousAt_id.sub continuousAt_const).log (sub_ne_zero.mpr hv)
  exact (h.pow 2).div_const 2

lemma continuousOn_Phi : ContinuousOn Phi (Set.Icc (Real.sqrt 2 - 1) (1 + Real.sqrt 2)) := by
  intro v hv
  have hs := s_gt
  have hv0 : 0 < v := by linarith [hv.1]
  have hva : v ≠ -(1 + Real.sqrt 2) := by intro h; linarith
  have c1 := continuousAt_PsiG (x := Real.sqrt 2 - 1) (y := 0) (v := v) (by linarith) hv0.ne'
  have c2 := continuousAt_PsiG (x := -(1 + Real.sqrt 2)) (y := 0) (v := v) (by linarith) hv0.ne'
  have c3 := continuousAt_logsq (x := 0) (v := v) hv0.ne'
  have c4 := continuousAt_PsiG (x := Real.sqrt 2 - 1) (y := -(1 + Real.sqrt 2)) (v := v) (by linarith) hva
  have c5 := continuousAt_logsq (x := -(1 + Real.sqrt 2)) (v := v) hva
  have c6 := continuousAt_PsiG (x := 0) (y := -(1 + Real.sqrt 2)) (v := v) (by linarith) hva
  have c7 : ContinuousAt (fun v : ℝ => Real.log v) v := Real.continuousAt_log hv0.ne'
  have c8 : ContinuousAt (fun v : ℝ => Real.log (v + (1 + Real.sqrt 2))) v :=
    (continuousAt_id.add continuousAt_const).log
      (by show v + (1 + Real.sqrt 2) ≠ 0; exact (by linarith : 0 < v + (1 + Real.sqrt 2)).ne')
  have := (((c1.add c2).sub c3).const_mul (1 - Real.sqrt 2 / 2)).sub
    ((((c4.add c5).sub c6)).const_mul 2) |>.sub
    (((c7.const_mul (1 - Real.sqrt 2 / 2)).sub (c8.const_mul 2)).const_mul (2 * lam))
  unfold Phi
  exact this.continuousWithinAt

lemma tt_b : tt (Real.sqrt 2 - 1) = 0 := by
  unfold tt
  have : (Real.sqrt 2 - 1) ^ 2 + 2 * (Real.sqrt 2 - 1) - 1 = 0 := by nlinarith [s_sq]
  rw [this]; simp

lemma tt_a : tt (1 + Real.sqrt 2) = 1 := by
  unfold tt
  have hs := s_pos
  rw [div_eq_one_iff_eq (by positivity)]
  nlinarith [s_sq]

lemma tt_mem {v : ℝ} (h1 : Real.sqrt 2 - 1 ≤ v) (h2 : v ≤ 1 + Real.sqrt 2) : tt v ∈ Set.Icc (0 : ℝ) 1 := by
  have hs := s_gt
  have hv : 0 < v := by linarith
  unfold tt
  constructor
  · apply div_nonneg _ (by linarith)
    nlinarith [s_sq]
  · rw [div_le_one (by linarith)]
    nlinarith [s_sq]

lemma tt_hasDerivAt {v : ℝ} (hv : v ≠ 0) : HasDerivAt tt ((v ^ 2 + 1) / (4 * v ^ 2)) v := by
  have h1 : HasDerivAt (fun v : ℝ => v ^ 2 + 2 * v - 1) (2 * v + 2) v := by
    have := (((hasDerivAt_id v).pow 2).add ((hasDerivAt_id v).const_mul 2)).sub_const 1
    simpa using this
  have h2 : HasDerivAt (fun v : ℝ => 4 * v) 4 v := by
    simpa using (hasDerivAt_id v).const_mul 4
  have := h1.div h2 (by positivity)
  unfold tt
  convert this using 1
  field_simp
  ring

lemma integral_eq : ∫ t in (0 : ℝ)..1, M1tIntegrand t = M1closed := by
  have hs := s_gt
  have hle : Real.sqrt 2 - 1 ≤ 1 + Real.sqrt 2 := by linarith
  have h : ∀ x ∈ Set.uIcc (Real.sqrt 2 - 1) (1 + Real.sqrt 2),
      HasDerivAt tt ((x ^ 2 + 1) / (4 * x ^ 2)) x := by
    intro x hx
    rw [Set.uIcc_of_le hle] at hx
    exact tt_hasDerivAt (by linarith [hx.1])
  have h' : ContinuousOn (fun x : ℝ => (x ^ 2 + 1) / (4 * x ^ 2))
      (Set.uIcc (Real.sqrt 2 - 1) (1 + Real.sqrt 2)) := by
    intro x hx
    rw [Set.uIcc_of_le hle] at hx
    have hx0' : 0 < x := by linarith [hx.1]
    have hx0 : x ≠ 0 := hx0'.ne'
    have : ContinuousAt (fun x : ℝ => (x ^ 2 + 1) / (4 * x ^ 2)) x :=
      (by fun_prop : Continuous fun x : ℝ => x ^ 2 + 1).continuousAt.div
        (by fun_prop : Continuous fun x : ℝ => 4 * x ^ 2).continuousAt (by positivity)
    exact this.continuousWithinAt
  have himg : tt '' Set.uIcc (Real.sqrt 2 - 1) (1 + Real.sqrt 2) ⊆ Set.Icc 0 1 := by
    rintro _ ⟨v, hv, rfl⟩
    rw [Set.uIcc_of_le hle] at hv
    exact tt_mem hv.1 hv.2
  have hg : ContinuousOn M1tIntegrand (tt '' Set.uIcc (Real.sqrt 2 - 1) (1 + Real.sqrt 2)) :=
    g_cont.mono himg
  have hsub := intervalIntegral.integral_comp_mul_deriv' h h' hg
  rw [tt_b, tt_a] at hsub
  rw [← hsub, ← Phi_diff]
  apply intervalIntegral.integral_eq_sub_of_hasDerivAt_of_le hle continuousOn_Phi
  · intro x hx
    exact hasDerivAt_Phi hx.1
  · apply ContinuousOn.intervalIntegrable
    rw [Set.uIcc_of_le hle]
    have hc1 : ContinuousOn tt (Set.Icc (Real.sqrt 2 - 1) (1 + Real.sqrt 2)) := by
      intro x hx
      have hx0' : 0 < x := by linarith [hx.1]
      exact (tt_hasDerivAt hx0'.ne').continuousAt.continuousWithinAt
    have hc2 : ContinuousOn (fun x => M1tIntegrand (tt x)) (Set.Icc (Real.sqrt 2 - 1) (1 + Real.sqrt 2)) :=
      g_cont.comp hc1 (fun x hx => tt_mem hx.1 hx.2)
    have h' := h'
    rw [Set.uIcc_of_le hle] at h'
    exact hc2.mul h'

lemma proof4 : stmt_4 := by
  refine ⟨?_, ?_⟩
  · exact (g_cont.integrableOn_Icc).mono_set Set.Ioo_subset_Icc_self
  · rw [integral_eq, ← proof2]

end D03

theorem D03.proof_1 : D03.stmt_1 := D03.proof1

theorem D03.proof_2 : D03.stmt_2 := D03.proof2

theorem D03.proof_3 : D03.stmt_3 := D03.proof3

theorem D03.proof_4 : D03.stmt_4 := D03.proof4
