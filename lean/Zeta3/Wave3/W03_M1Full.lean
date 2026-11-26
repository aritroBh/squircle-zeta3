import Mathlib
import Zeta3.Wave3.W02_Holo
import Zeta3.Wave2.D03_M1RightAngle
/-!
# W03 -- capstone: `M_0` and `M_1` of the squircle expansion in closed form (Theorem `thm:M1`)

Target id: W03.
Paper labels: `thm:M1` (main.tex 889-895), `lem:holo` with `eq:expansion` (main.tex 846-852),
`eq:M1t` (main.tex 858-861), and the unlabelled display `M_0 = λ - √2 L` (main.tex 880-883).
Constants (main.tex 870-875): `a = 1+√2`, `L = log a`, `λ = log 2`.

Quoted verbatim (main.tex, Lemma lem:holo with eq:expansion):

  $F$ extends to a holomorphic function on $\{\operatorname{Re}\varepsilon<1\}$. Hence
  $P(p)/8=1+\varepsilon M(\varepsilon)$ with $M=-F/2$ holomorphic in $|\varepsilon|<1$, the expansion
  P(p)=8\Bigl(1+\frac{M_0}{p}+\frac{M_1}{p^2}+\frac{M_2}{p^3}+\frac{M_3}{p^4}+\cdots\Bigr)
  converges for every $p>1$, and $M_k=-F^{(k)}(0)/(2\,k!)$.

Quoted verbatim (main.tex, eq:M1t):

  M_1=\int_0^1\frac{\log t}{t}\Bigl(\frac{1-t}{q}-1\Bigr)dt,\qquad q=\sqrt{t^2+(1-t)^2},

Quoted verbatim (main.tex 880-883, subsection "M_0: Logarithms Only"):

  Ahmed~\cite{ahmed} proved $P(p)=8-4K/p+O(p^{-2})$ with $K=2\sqrt2\log(1+\sqrt2)-2\log2$.
  In our notation
  \[
  M_0=\lambda-\sqrt2\,L=-0.55330329972051571737\ldots,
  \]

Quoted verbatim (main.tex 889-895, Theorem thm:M1):

  \begin{theorem}[$M_1$; answers Problem~1 of Ahmed~\cite{ahmed}]\label{thm:M1}
  \begin{multline*}
  M_1=\Bigl(\frac{19\sqrt2}{48}-\frac1{24}\Bigr)\pi^2+\frac{\lambda^2}{2}+2\sqrt2\,\lambda L\\
  -\Bigl(1+\frac{5\sqrt2}{2}\Bigr)L^2-5\sqrt2\,\Li_2(\sqrt2-1).
  \end{multline*}
  Numerically $\kappa=8M_1=2.02095276965635835695\ldots$
  \end{theorem}

## Objects used (imported, not redefined)

* `W02.Pint p = (4/p) ∫_0^1 √(t^(2/p-2) + (1-t)^(2/p-2)) dt` (eq:Pint), `W02.gR`, `W02.FR`
  (the real `F` of eq:F), `W02.FC` (its holomorphic extension of lem:holo) and
  `W02.Mk k = -F^{(k)}(0)/(2 k!)` (complex derivatives of `W02.FC`). W02 proves
  (`W02.proof_6`) that `Σ_{k≥0} M_k / p^(k+1) = P(p)/8 - 1` for every real `p > 1` (`HasSum`),
  and (`W02.proof_7`) that `M_1` equals the integral eq:M1t.
* `D03.lam = log 2`, `D03.L = log(1+√2)`, `D03.M1closed` = the right side of thm:M1 with the
  dilogarithm `D00.Li2` (`Li_2(x) = -∫_0^x log(1-t)/t dt`; its argument `√2-1` is in `(0,1)`,
  where this is the series `Σ x^m/m²`). D03 proves that the integral eq:M1t equals `M1closed`
  (`D03.proof_4`, `D03.proof_2`) and that `|8·M1closed - 2.02095276965635835695| < 10⁻²⁰`
  (via `D03.proof_3`).

## How the statements fit together

`M_0 = -F(0)/2`, and `F(0) = ∫_0^1 (1/t + 1/(1-t) - √(1/t² + 1/(1-t)²)) dt` (eq:F at `ε = 0`;
in Lean `W02.FR 0`, with `t^(0-1) = t⁻¹` for `t > 0`). The paper states `M_0 = λ - √2 L`
without proof (it credits Ahmed and [MRL]); the evaluation `F(0) = 2√2 L - 2λ` is stmt_1, and
`W02.proof_4` (`FC ε = FR ε` for real `ε < 1`) turns it into stmt_2.
`M_1` goes through `W02.proof_7` (eq:M1t) and `D03` (eq:M1t = closed form), giving stmt_4.

Numerical check of stmt_1 (40 digits, from D03's notes and supplement sec_g5_decay):
`λ - √2 L = -0.55330329972051571737…`, so `F(0) = 1.10660659944103143474…`.

## Statements

* stmt_1: `F(0)` in closed form: the integrand of `F` at `ε = 0` is integrable on `[0,1]` and
  `F(0) = 2√2 L - 2λ`.
* stmt_2: the display `M_0 = λ - √2 L`, for `M_0 = W02.Mk 0` (as a complex number).
* stmt_3: the printed decimals: `‖M_0 - (-0.55330329972051571737)‖ < 10⁻²⁰` (the printed
  digits are all correct; the ellipsis means truncation, and the true value lies within
  `10⁻²⁰` of the printed decimal either way).
* stmt_4: Theorem thm:M1: `M_1 = W02.Mk 1` equals the closed form `D03.M1closed`.
* stmt_5: thm:M1, "Numerically `κ = 8M_1 = 2.02095276965635835695…`":
  `‖8 M_1 - 2.02095276965635835695‖ < 10⁻²⁰`, for `M_1 = W02.Mk 1`.
* stmt_6 (final, the combined statement): for every real `p > 1`,
  `P(p)/8 = 1 + M_0/p + M_1/p² + Σ_{k≥2} M_k/p^(k+1)` with `M_0 = λ - √2 L` and `M_1` the closed
  form of thm:M1, in `HasSum` form (the tail series converges, to the stated value).

## NOT COVERED

* That `W02.Pint p` (the integral eq:Pint) is the perimeter of the squircle `|x|^p + |y|^p = 1`
  (the arc-length computation "the arc-length formula gives") is target W04, not this file.
  Here `P(p)` always means the integral eq:Pint.
* The coefficients `M_k` for `k ≥ 2` are not evaluated (they appear only through `W02.Mk` in the
  tail of stmt_6); `M_2`, `M_3` are other targets.
* Ahmed's form `P(p) = 8 - 4K/p + O(p⁻²)` with `K = 2√2 log(1+√2) - 2 log 2` (a citation) is not
  stated separately; it is the case `M_0 = -K/2` of stmt_6 plus a tail bound.
* The identification `κ = 8 M_1` with Ahmed's numerically found constant and "It matches Ahmed's
  ten digits" (citations) are not formalized; stmt_5 states the decimal for `8 M_1` itself.
* The `O`-rate claims of Fig. `fig:error` (each coefficient buys one power of `p`) are not stated.
-/
-- BEGIN STATEMENTS
namespace W03

open Real

/-- `M_0 = λ - √2 L` (main.tex 882), the claimed closed form of `M_0`. -/
noncomputable def M0closed : ℝ := D03.lam - Real.sqrt 2 * D03.L

/-- `M_0 = -F(0)/2` with `F(0) = ∫_0^1 (1/t + 1/(1-t) - √(1/t² + 1/(1-t)²)) dt` (eq:F at
`ε = 0`): the integrand is integrable on `[0,1]` and `F(0) = 2√2 L - 2λ`. This evaluation is
used, not printed, by the display `M_0 = λ - √2 L` (main.tex 882). -/
def stmt_1 : Prop :=
  IntervalIntegrable (W02.gR 0) MeasureTheory.volume 0 1 ∧
    W02.FR 0 = 2 * Real.sqrt 2 * D03.L - 2 * D03.lam

/-- main.tex 882: "$M_0=\lambda-\sqrt2\,L$", for `M_0 = -F(0)/2 = W02.Mk 0`. -/
def stmt_2 : Prop := W02.Mk 0 = ((M0closed : ℝ) : ℂ)

/-- main.tex 882: "$M_0=\dots=-0.55330329972051571737\ldots$", to all printed digits. -/
def stmt_3 : Prop :=
  ‖W02.Mk 0 - ((-0.55330329972051571737 : ℝ) : ℂ)‖ < (1 : ℝ) / 10 ^ 20

/-- Theorem thm:M1: the coefficient `M_1` of `1/p²` in eq:expansion,
`M_1 = -F'(0)/2 = W02.Mk 1`, equals
`(19√2/48 - 1/24)π² + λ²/2 + 2√2 λL - (1 + 5√2/2)L² - 5√2 Li₂(√2-1)`. -/
def stmt_4 : Prop := W02.Mk 1 = ((D03.M1closed : ℝ) : ℂ)

/-- Theorem thm:M1: "Numerically $\kappa=8M_1=2.02095276965635835695\ldots$". -/
def stmt_5 : Prop :=
  ‖8 * W02.Mk 1 - ((2.02095276965635835695 : ℝ) : ℂ)‖ < (1 : ℝ) / 10 ^ 20

/-- eq:expansion with the two closed forms: for every real `p > 1`,
`P(p)/8 = 1 + M_0/p + M_1/p² + Σ_{k≥2} M_k/p^(k+1)`, `M_0 = λ - √2 L`, `M_1` as in thm:M1;
the tail `Σ_{j≥0} M_{j+2}/p^(j+3)` converges to `P(p)/8 - 1 - M_0/p - M_1/p²`. -/
def stmt_6 : Prop :=
  ∀ p : ℝ, 1 < p →
    HasSum (fun j : ℕ => W02.Mk (j + 2) / (p : ℂ) ^ (j + 3))
      (((W02.Pint p / 8 : ℝ) : ℂ) - 1 - ((M0closed : ℝ) : ℂ) / (p : ℂ)
        - ((D03.M1closed : ℝ) : ℂ) / (p : ℂ) ^ 2)

end W03
-- END STATEMENTS

namespace W03
open Real

lemma s2_sq : Real.sqrt 2 ^ 2 = 2 := Real.sq_sqrt (by norm_num)
lemma s2_pos : 0 < Real.sqrt 2 := by positivity

lemma A_pos (t : ℝ) : 0 < Real.sqrt 2 * W02.q t + 2 * t - 1 := by
  have hq := W02.q_pos t
  have hq2 := W02.q_sq t
  have hx : 0 < Real.sqrt 2 * W02.q t := mul_pos s2_pos hq
  have hx2 : (Real.sqrt 2 * W02.q t) ^ 2 = (2 * t - 1) ^ 2 + 1 := by
    rw [mul_pow, s2_sq, hq2]; ring
  by_contra h0
  have h := not_lt.mp h0
  nlinarith [mul_nonneg (by linarith : (0:ℝ) ≤ 1 - 2 * t - Real.sqrt 2 * W02.q t)
    (by linarith : (0:ℝ) ≤ 1 - 2 * t + Real.sqrt 2 * W02.q t)]

lemma B_pos (t : ℝ) : 0 < W02.q t - t + 1 := by
  have hq := W02.q_pos t
  have hq2 := W02.q_sq t
  by_contra h0
  have h := not_lt.mp h0
  nlinarith [mul_nonneg (by linarith : (0:ℝ) ≤ t - 1 - W02.q t) (by linarith : (0:ℝ) ≤ t - 1 + W02.q t)]

noncomputable def G (t : ℝ) : ℝ :=
  -W02.q t + (1 / Real.sqrt 2) * Real.log (2 * Real.sqrt 2 * W02.q t + 4 * t - 2)
    + Real.log (2 * W02.q t - 2 * t + 2)

lemma q_hasDerivAt (t : ℝ) : HasDerivAt W02.q ((2 * t - 1) / W02.q t) t := by
  have h1 : HasDerivAt (fun y : ℝ => y ^ 2 + (1 - y) ^ 2) (4 * t - 2) t := by
    have := ((hasDerivAt_id t).pow 2).add (((hasDerivAt_id t).const_sub 1).pow 2)
    refine this.congr_deriv ?_
    simp; ring
  have h2 := h1.sqrt (by nlinarith [sq_nonneg (t - 1/2)])
  refine h2.congr_deriv ?_
  have e : √(t^2+(1-t)^2) = W02.q t := rfl
  rw [e]; field_simp; ring

lemma G_hasDerivAt (t : ℝ) :
    HasDerivAt G ((W02.q t + 1 - 2 * t) / (W02.q t - t + 1)) t := by
  have hq := W02.q_pos t
  have hq2 := W02.q_sq t
  have hA := A_pos t
  have hB := B_pos t
  have hd := q_hasDerivAt t
  have hA' : HasDerivAt (fun y => 2 * Real.sqrt 2 * W02.q y + 4 * y - 2)
      (2 * Real.sqrt 2 * ((2 * t - 1) / W02.q t) + 4) t := by
    have := ((hd.const_mul (2 * Real.sqrt 2)).add ((hasDerivAt_id t).const_mul 4)).sub_const 2
    exact this.congr_deriv (by ring)
  have hB' : HasDerivAt (fun y => 2 * W02.q y - 2 * y + 2)
      (2 * ((2 * t - 1) / W02.q t) - 2) t := by
    have := ((hd.const_mul 2).sub ((hasDerivAt_id t).const_mul 2)).add_const 2
    exact this.congr_deriv (by ring)
  have hAne : 2 * Real.sqrt 2 * W02.q t + 4 * t - 2 ≠ 0 := by nlinarith
  have hBne : 2 * W02.q t - 2 * t + 2 ≠ 0 := by nlinarith
  have := (hd.neg.add ((hA'.log hAne).const_mul (1 / Real.sqrt 2))).add (hB'.log hBne)
  have hs := s2_pos
  have hs2 := s2_sq
  have e1 : 1 / √2 * ((2 * √2 * ((2 * t - 1) / W02.q t) + 4) / (2 * √2 * W02.q t + 4 * t - 2))
      = 1 / W02.q t := by
    field_simp
    rw [div_eq_iff (by nlinarith)]
    linear_combination (-2 * W02.q t) * hs2
  have e2 : (2 * ((2 * t - 1) / W02.q t) - 2) / (2 * W02.q t - 2 * t + 2)
      = (2 * t - 1 - W02.q t) / (W02.q t * (W02.q t - t + 1)) := by
    field_simp
  refine this.congr_deriv ?_
  rw [e1, e2]
  field_simp
  linear_combination (-1) * hq2


lemma G_deriv_eq {t : ℝ} (ht : t ≠ 0) :
    (W02.q t + 1 - 2 * t) / (W02.q t - t + 1) = (1 - W02.q t) / t := by
  have hB := B_pos t
  have hq2 := W02.q_sq t
  rw [div_eq_div_iff hB.ne' ht]
  linear_combination hq2

noncomputable def H (t : ℝ) : ℝ := G t - G (1 - t)

lemma H_hasDerivAt (t : ℝ) :
    HasDerivAt H ((W02.q t + 1 - 2 * t) / (W02.q t - t + 1)
      + (W02.q (1 - t) + 1 - 2 * (1 - t)) / (W02.q (1 - t) - (1 - t) + 1)) t := by
  have h1 := G_hasDerivAt t
  have h2 := (G_hasDerivAt (1 - t)).comp t ((hasDerivAt_id t).const_sub 1)
  exact (h1.sub h2).congr_deriv (by simp)

lemma gR0_eq {t : ℝ} (ht0 : 0 < t) (ht1 : t < 1) :
    W02.gR 0 t = (1 - W02.q t) / t + (1 - W02.q t) / (1 - t) := by
  have h1 : 0 < 1 - t := by linarith
  have hq := W02.q_pos t
  have hq2 := W02.q_sq t
  have e : (t⁻¹) ^ 2 + ((1 - t)⁻¹) ^ 2 = (W02.q t / (t * (1 - t))) ^ 2 := by
    rw [div_pow, hq2]; field_simp; ring
  simp only [W02.gR, W02.aR, W02.bR, zero_sub, Real.rpow_neg_one]
  rw [e, Real.sqrt_sq (by positivity)]
  field_simp
  ring

lemma H_one_sub_zero : H 1 - H 0 = 2 * Real.sqrt 2 * D03.L - 2 * D03.lam := by
  have q1 : W02.q 1 = 1 := by simp [W02.q]
  have q0 : W02.q 0 = 1 := by simp [W02.q]
  have hs := s2_pos
  have hs2 := s2_sq
  have e1 : Real.log (2 * Real.sqrt 2 * 1 + 4 * 1 - 2) = D03.lam + D03.L := by
    rw [show 2 * Real.sqrt 2 * 1 + 4 * 1 - 2 = 2 * (1 + Real.sqrt 2) by ring,
      Real.log_mul (by norm_num) (by positivity)]
    rfl
  have e0 : Real.log (2 * Real.sqrt 2 * 1 + 4 * 0 - 2) = D03.lam - D03.L := by
    rw [show 2 * Real.sqrt 2 * 1 + 4 * 0 - 2 = 2 * (Real.sqrt 2 - 1) by ring,
      Real.log_mul (by norm_num) (by nlinarith), D03.log_b]
    rfl
  have f1 : Real.log (0 + 2) = D03.lam := by norm_num [D03.lam]
  have f0 : Real.log (2 * 1 - 2 * 0 + 2) = 2 * D03.lam := by
    rw [show (2 : ℝ) * 1 - 2 * 0 + 2 = 4 by norm_num]; exact D03.log_four
  simp only [H, G, sub_zero, sub_self, q1, q0, e1, e0, f1, f0]
  field_simp
  linear_combination (-2 * D03.L) * hs2

lemma FR0 : W02.FR 0 = 2 * Real.sqrt 2 * D03.L - 2 * D03.lam := by
  have hint := (W02.proof_4 0 (by norm_num)).1
  rw [W02.FR, intervalIntegral.integral_eq_sub_of_hasDerivAt_of_le zero_le_one
    (fun x _ => (H_hasDerivAt x).continuousAt.continuousWithinAt) _ hint, H_one_sub_zero]
  intro x hx
  refine (H_hasDerivAt x).congr_deriv ?_
  rw [gR0_eq hx.1 hx.2, W02.q_symm, G_deriv_eq hx.1.ne', ← W02.q_symm x,
    G_deriv_eq (by linarith [hx.2] : (1 : ℝ) - x ≠ 0), W02.q_symm]


lemma proof1 : stmt_1 := ⟨(W02.proof_4 0 (by norm_num)).1, FR0⟩

lemma proof2 : stmt_2 := by
  have h := (W02.proof_4 0 (by norm_num)).2
  simp only [Complex.ofReal_zero] at h
  unfold stmt_2 W02.Mk M0closed
  rw [iteratedDeriv_zero, h, FR0]
  push_cast
  ring

lemma proof3 : stmt_3 := by
  unfold stmt_3
  rw [proof2, ← Complex.ofReal_sub, Complex.norm_real, Real.norm_eq_abs, M0closed]
  obtain ⟨hs1, hs2⟩ := D03.s_bounds
  obtain ⟨hl1, hl2⟩ := D03.lam_bounds
  obtain ⟨hL1, hL2⟩ := D03.L_bounds
  have p : Real.sqrt 2 * D03.L ≤ (1.41421356237309504880168873 : ℝ) * (0.88137358701954302523260934 : ℝ) :=
    mul_le_mul hs2 hL2 (by linarith) (by norm_num)
  have p' : (1.41421356237309504880168872 : ℝ) * (0.88137358701954302523260931 : ℝ) ≤ Real.sqrt 2 * D03.L :=
    mul_le_mul hs1 hL1 (by norm_num) (by linarith)
  rw [abs_lt]
  constructor <;> norm_num at p p' ⊢ <;> linarith

lemma proof4 : stmt_4 := by
  unfold stmt_4
  rw [W02.proof_7.2]
  congr 1
  exact D03.integral_eq

lemma proof5 : stmt_5 := by
  unfold stmt_5
  rw [proof4, show (8 : ℂ) * ((D03.M1closed : ℝ) : ℂ) - ((2.02095276965635835695 : ℝ) : ℂ)
      = ((8 * D03.M1closed - 2.02095276965635835695 : ℝ) : ℂ) by push_cast; ring,
    Complex.norm_real, Real.norm_eq_abs, ← D03.proof_2]
  exact D03.proof_3

lemma proof6 : stmt_6 := by
  intro p hp
  have h := (hasSum_nat_add_iff' 2).mpr (W02.proof_6 p hp)
  convert h using 1
  rw [Finset.sum_range_succ, Finset.sum_range_succ, Finset.sum_range_zero, proof2, proof4]
  push_cast
  ring

end W03

theorem W03.proof_1 : W03.stmt_1 := W03.proof1

theorem W03.proof_2 : W03.stmt_2 := W03.proof2

theorem W03.proof_3 : W03.stmt_3 := W03.proof3

theorem W03.proof_4 : W03.stmt_4 := W03.proof4

theorem W03.proof_5 : W03.stmt_5 := W03.proof5

theorem W03.proof_6 : W03.stmt_6 := W03.proof6
