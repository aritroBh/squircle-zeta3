import Mathlib
import Zeta3.Wave2.D00_Dilog
/-!
# D02 -- the dilogarithm values at `√2`-arguments (Lemma `lem:dilog`)

Target id: D02.
Paper label: `lem:dilog` (paper/ieee/main.tex, lines 2332-2343), Appendix "Proof of Theorem
thm:M1". Setup: main.tex line 875 (`L = log a`, `λ = log 2`, with `a = 1+√2`, see the caption at
main.tex 918: "$a=1+\sqrt2$, $b=\sqrt2-1$, $L=\log a$, $\lambda=\log2$"), and main.tex 2325
("We use the dilogarithm on $(-\infty,1]$ ...").

Quoted verbatim (main.tex 2332-2343):

  \begin{lemma}\label{lem:dilog}
  With $b=\sqrt2-1$,
  \begin{align*}
  \Li_2(b)-\Li_2(-b)&=\tfrac{\pi^2}{8}-\tfrac12L^2,\\
  \Li_2(b^2)&=4\Li_2(b)-\tfrac{\pi^2}{4}+L^2,\\
  \Li_2(2-\sqrt2)&=-\Li_2(b)+\tfrac{\pi^2}{6}+\tfrac12\lambda L-L^2,\\
  \Li_2(-\sqrt2)&=\Li_2(b)-\tfrac{\pi^2}{6}-\tfrac12\lambda L+\tfrac12L^2,\\
  \Li_2(-b^2)&=6\Li_2(b)-\tfrac{11\pi^2}{24}+2L^2,\\
  \Li_2\bigl(\tfrac{2+\sqrt2}{4}\bigr)&=6\Li_2(b)-\tfrac{7\pi^2}{24}-\tfrac98\lambda^2\\
  &\quad-\tfrac32\lambda L+\tfrac72L^2.
  \end{align*}
  \end{lemma}

Quoted verbatim (main.tex 875): L=\log a,\qquad \lambda=\log2 .

The align environment has seven printed lines but six identities: the last identity spans the
last two printed lines (the `\quad` continuation). Each identity is its own statement
(stmt_1 .. stmt_6). `Li_2` is `D00.Li2` throughout.

Argument check: the arguments are `b ≈ 0.414`, `-b`, `b² ≈ 0.172`, `2-√2 ≈ 0.586`, `-√2`,
`-b² `, `(2+√2)/4 ≈ 0.854`; all lie in `(-∞,1]` (indeed in `(-∞,1)`), where `D00.Li2` is the
genuine real dilogarithm. stmt_7 records this as a formal claim.

Numerical sanity check (Simpson quadrature of the defining integral, double precision): all six
differences are below `4e-14`.

NOT COVERED: nothing of the lemma's statement is omitted. (The proof's intermediate identities,
e.g. the five-term instance, are not stated separately; they belong to the proof.)
-/
-- BEGIN STATEMENTS
namespace D02

open Real

/-- `b = √2 - 1` (main.tex 2333). -/
noncomputable def b : ℝ := Real.sqrt 2 - 1

/-- `a = 1 + √2` (main.tex 918). -/
noncomputable def a : ℝ := 1 + Real.sqrt 2

/-- `L = log a` (main.tex 875). -/
noncomputable def L : ℝ := Real.log a

/-- `λ = log 2` (main.tex 875). -/
noncomputable def lam : ℝ := Real.log 2

/-- Line 1: `Li₂(b) - Li₂(-b) = π²/8 - ½L²`. -/
def stmt_1 : Prop :=
  D00.Li2 b - D00.Li2 (-b) = π ^ 2 / 8 - (1 / 2) * L ^ 2

/-- Line 2: `Li₂(b²) = 4Li₂(b) - π²/4 + L²`. -/
def stmt_2 : Prop :=
  D00.Li2 (b ^ 2) = 4 * D00.Li2 b - π ^ 2 / 4 + L ^ 2

/-- Line 3: `Li₂(2-√2) = -Li₂(b) + π²/6 + ½λL - L²`. -/
def stmt_3 : Prop :=
  D00.Li2 (2 - Real.sqrt 2) = -D00.Li2 b + π ^ 2 / 6 + (1 / 2) * lam * L - L ^ 2

/-- Line 4: `Li₂(-√2) = Li₂(b) - π²/6 - ½λL + ½L²`. -/
def stmt_4 : Prop :=
  D00.Li2 (-Real.sqrt 2) = D00.Li2 b - π ^ 2 / 6 - (1 / 2) * lam * L + (1 / 2) * L ^ 2

/-- Line 5: `Li₂(-b²) = 6Li₂(b) - 11π²/24 + 2L²`. -/
def stmt_5 : Prop :=
  D00.Li2 (-b ^ 2) = 6 * D00.Li2 b - 11 * π ^ 2 / 24 + 2 * L ^ 2

/-- Lines 6-7 (one identity): `Li₂((2+√2)/4) = 6Li₂(b) - 7π²/24 - 9/8 λ² - 3/2 λL + 7/2 L²`. -/
def stmt_6 : Prop :=
  D00.Li2 ((2 + Real.sqrt 2) / 4)
    = 6 * D00.Li2 b - 7 * π ^ 2 / 24 - (9 / 8) * lam ^ 2 - (3 / 2) * lam * L
      + (7 / 2) * L ^ 2

/-- Every argument of `Li₂` in the lemma lies in `(-∞,1]`, the paper's domain of the real
dilogarithm (main.tex 2325). -/
def stmt_7 : Prop :=
  b ≤ 1 ∧ -b ≤ 1 ∧ b ^ 2 ≤ 1 ∧ 2 - Real.sqrt 2 ≤ 1 ∧ -Real.sqrt 2 ≤ 1 ∧ -b ^ 2 ≤ 1 ∧
    (2 + Real.sqrt 2) / 4 ≤ 1

end D02
-- END STATEMENTS

namespace D02

open Real

lemma hs : Real.sqrt 2 * Real.sqrt 2 = 2 := Real.mul_self_sqrt (by norm_num)

lemma hs1 : 1 < Real.sqrt 2 := by
  nlinarith [hs, Real.sqrt_nonneg 2]

lemma hs2 : Real.sqrt 2 < 3 / 2 := by
  nlinarith [hs, Real.sqrt_nonneg 2]

lemma b_pos : 0 < b := by unfold b; linarith [hs1]

lemma b_lt : b < 1 := by unfold b; linarith [hs2]

lemma a_pos : 0 < a := by unfold a; linarith [hs1]

lemma ab : a * b = 1 := by unfold a b; nlinarith [hs]

lemma log_b : Real.log b = -L := by
  have : b = a⁻¹ := by
    exact eq_inv_of_mul_eq_one_right ab
  rw [this, Real.log_inv]; rfl

lemma log_s : Real.log (Real.sqrt 2) = lam / 2 := by
  rw [Real.log_sqrt (by norm_num)]; rfl

lemma log_1mb : Real.log (1 - b) = lam / 2 - L := by
  have h : 1 - b = Real.sqrt 2 * b := by unfold b; nlinarith [hs]
  rw [h, Real.log_mul (by linarith [hs1]) b_pos.ne', log_s, log_b]; ring

lemma log_1pb : Real.log (1 + b) = lam / 2 := by
  have h : 1 + b = Real.sqrt 2 := by unfold b; ring
  rw [h, log_s]

lemma p1 : stmt_1 := by
  have h := D00.chi2_landen b_pos b_lt
  have e1 : (1 - b) / (1 + b) = b := by
    have : (1 + b) ≠ 0 := by linarith [b_pos]
    rw [div_eq_iff this]; unfold b; nlinarith [hs]
  have e2 : (1 + b) / (1 - b) = a := by
    have : (1 - b) ≠ 0 := by linarith [b_lt]
    rw [div_eq_iff this]; unfold b a; nlinarith [hs]
  rw [e1, e2, log_b] at h
  unfold D00.chi2 at h
  unfold stmt_1
  have hL : Real.log a = L := rfl
  rw [hL] at h
  nlinarith [h]

lemma p2 : stmt_2 := by
  have h1 := p1
  unfold stmt_1 at h1
  have hd := D00.duplication (x := b) (by rw [abs_of_pos b_pos]; exact b_lt.le)
  unfold stmt_2
  linarith

lemma p3 : stmt_3 := by
  have h := D00.reflection b_pos b_lt
  rw [log_b, log_1mb] at h
  have e : (2 - Real.sqrt 2) = 1 - b := by unfold b; ring
  unfold stmt_3
  rw [e]
  linear_combination h

lemma p4 : stmt_4 := by
  have h := D00.landen (x := -Real.sqrt 2) (by linarith [hs1])
  have e : -Real.sqrt 2 / (-Real.sqrt 2 - 1) = 2 - Real.sqrt 2 := by
    have : (-Real.sqrt 2 - 1) ≠ 0 := by linarith [hs1]
    rw [div_eq_iff this]; nlinarith [hs]
  have e2 : (1 - -Real.sqrt 2) = a := by unfold a; ring
  rw [e, e2] at h
  have h3 := p3
  unfold stmt_3 at h3
  unfold stmt_4
  have hL : Real.log a = L := rfl
  rw [hL] at h
  nlinarith [h, h3]

lemma p5 : stmt_5 := by
  have h := D00.five_term (x := -Real.sqrt 2) (y := 1 - Real.sqrt 2) (by linarith [hs1])
    (by linarith [hs1]) (by linarith [hs1])
  have hbb : b = Real.sqrt 2 - 1 := rfl
  have e1 : -Real.sqrt 2 / (1 - (1 - Real.sqrt 2)) = -1 := by
    have : Real.sqrt 2 ≠ 0 := by linarith [hs1]
    field_simp
    ring
  have e2 : (1 - Real.sqrt 2) / (1 - -Real.sqrt 2) = -b ^ 2 := by
    have : (1 - -Real.sqrt 2) ≠ 0 := by linarith [hs1]
    rw [div_eq_iff this]; unfold b; nlinarith [hs]
  have e3 : -Real.sqrt 2 * (1 - Real.sqrt 2) / ((1 - -Real.sqrt 2) * (1 - (1 - Real.sqrt 2)))
      = b ^ 2 := by
    have : ((1 - -Real.sqrt 2) * (1 - (1 - Real.sqrt 2))) ≠ 0 := by
      have : 0 < (1 - -Real.sqrt 2) * (1 - (1 - Real.sqrt 2)) := by
        apply mul_pos <;> linarith [hs1]
      exact this.ne'
    rw [div_eq_iff this]; unfold b; nlinarith [hs]
  have e4 : 1 - Real.sqrt 2 = -b := by unfold b; ring
  have e5 : (1 - -Real.sqrt 2) = a := by unfold a; ring
  have e6 : (1 - (1 - Real.sqrt 2)) = Real.sqrt 2 := by ring
  rw [e1, e2, e3, e5, e6] at h
  rw [e4] at h
  rw [log_s] at h
  have hL : Real.log a = L := rfl
  rw [hL] at h
  have h2 := p2
  have h1 := p1
  have h4 := p4
  unfold stmt_2 at h2
  unfold stmt_1 at h1
  unfold stmt_4 at h4
  unfold stmt_5
  rw [D00.Li2_neg_one] at h
  nlinarith [h, h1, h2, h4]

lemma p6 : stmt_6 := by
  have hc : (2 + Real.sqrt 2) / 4 < 1 := by linarith [hs2]
  have h := D00.landen hc
  have e : (2 + Real.sqrt 2) / 4 / ((2 + Real.sqrt 2) / 4 - 1) = -a ^ 2 := by
    have : ((2 + Real.sqrt 2) / 4 - 1) ≠ 0 := by linarith [hs2]
    rw [div_eq_iff this]; unfold a; nlinarith [hs]
  have e2 : (1 - (2 + Real.sqrt 2) / 4) = Real.sqrt 2 * b / 4 := by unfold b; nlinarith [hs]
  have hl4 : Real.log (1 - (2 + Real.sqrt 2) / 4) = -(3 / 2) * lam - L := by
    rw [e2, Real.log_div (mul_ne_zero (by linarith [hs1]) b_pos.ne') (by norm_num),
      Real.log_mul (by linarith [hs1]) b_pos.ne', log_s, log_b,
      show (4:ℝ) = 2 ^ 2 by norm_num, Real.log_pow]
    unfold lam; push_cast; ring
  have hi := D00.inversion_neg (x := -a ^ 2) (by nlinarith [a_pos])
  have e3 : 1 / (-a ^ 2) = -b ^ 2 := by
    have : b = a⁻¹ := eq_inv_of_mul_eq_one_right ab
    rw [this]; field_simp
  have hl5 : Real.log (- -a ^ 2) = 2 * L := by
    rw [neg_neg, Real.log_pow]; push_cast; rfl
  rw [e3, hl5] at hi
  rw [e, hl4] at h
  have h5 := p5
  unfold stmt_5 at h5
  unfold stmt_6
  nlinarith [h, hi, h5]

end D02

theorem D02.proof_1 : D02.stmt_1 := D02.p1
theorem D02.proof_2 : D02.stmt_2 := D02.p2
theorem D02.proof_3 : D02.stmt_3 := D02.p3
theorem D02.proof_4 : D02.stmt_4 := D02.p4
theorem D02.proof_5 : D02.stmt_5 := D02.p5
theorem D02.proof_6 : D02.stmt_6 := D02.p6
theorem D02.proof_7 : D02.stmt_7 := by
  have := D02.hs1
  have := D02.hs2
  have := D02.b_pos
  have := D02.b_lt
  unfold D02.stmt_7
  refine ⟨D02.b_lt.le, by linarith, ?_, by linarith, by linarith, ?_, by linarith⟩
  · nlinarith
  · nlinarith [sq_nonneg D02.b]
