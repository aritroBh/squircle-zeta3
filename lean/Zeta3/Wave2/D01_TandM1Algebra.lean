import Mathlib
import Zeta3.Wave2.D00_Dilog
/-!
# D01 -- the weight-two term `T` of `thm:M2` and the closing algebra of the `M_1` proof

Target id: D01.
Paper labels: `thm:M2` equation `(T)` (label `eq:T`, main.tex lines 901-907); the last paragraph
of the proof of `thm:M1` in Appendix `app:M1` (main.tex lines 2349-2368), together with the
statement of `thm:M1` (main.tex lines 889-896) and the section constants (main.tex lines 873-876).

Quoted verbatim (main.tex 873-876):

  \begin{gather*}
  a=1+\sqrt2,\qquad b=\sqrt2-1=\tfrac1a=\tan\tfrac\pi8,\\
  L=\log a,\qquad \lambda=\log2 .
  \end{gather*}

Quoted verbatim (main.tex 903-906, inside thm:M2):

  \begin{equation}\label{eq:T}
  T=\frac{\pi^2+6\lambda^2}{24}=\frac{\zeta(2)-\Li_2(\frac12)}{2}
  \end{equation}

Quoted verbatim (main.tex 889-894, thm:M1):

  M_1=\Bigl(\frac{19\sqrt2}{48}-\frac1{24}\Bigr)\pi^2+\frac{\lambda^2}{2}+2\sqrt2\,\lambda L\\
  -\Bigl(1+\frac{5\sqrt2}{2}\Bigr)L^2-5\sqrt2\,\Li_2(\sqrt2-1).

Quoted verbatim (main.tex 2350-2368, end of the proof of thm:M1):

  with $c_*=1-1/\sqrt2$ ...
  I(b,0)&=\Li_2(b^2)-\tfrac{\pi^2}{6},\\
  I(-a,0)&=\Li_2(-b^2)+\tfrac{\pi^2}{12}+2L^2,\\
  I(0,0)&=0,\qquad I(-a,-a)=-\tfrac58\lambda^2+\lambda L+\tfrac12L^2,\\
  I(b,-a)&=\Li_2(2-\sqrt2)-\tfrac{\pi^2}{6}-\tfrac58\lambda^2\\
  &\quad+\lambda L+\tfrac12L^2,\\
  I(0,-a)&=-\Li_2\bigl(\tfrac{2+\sqrt2}{4}\bigr)+\tfrac{\pi^2}{12}-\tfrac98\lambda^2\\
  &\quad+\lambda L+\tfrac12L^2.
  ...
  M_1=c_*\bigl[I(b,0)+I(-a,0)-I(0,0)\bigr]\\
  -2\bigl[I(b,-a)+I(-a,-a)-I(0,-a)\bigr]-2\lambda^2+2\sqrt2\,\lambda L .
  By Lemma~\ref{lem:dilog} the bracket multiplying $c_*$ is $10\Li_2(b)-\frac{19}{24}\pi^2+5L^2$,
  and the rest is $\frac34\pi^2+\frac12\lambda^2+2\sqrt2\lambda L-6L^2-10\Li_2(b)$. Adding gives
  Theorem~\ref{thm:M1}.

## Conventions

* `Li_2` is the library dilogarithm `D00.Li2` throughout. All its arguments here
  (`1/2`, `b`, `b^2`, `-b^2`, `2-√2`, `(2+√2)/4`) lie in `(-∞,1)`, where `D00.Li2` is the
  genuine dilogarithm, so no branch convention is involved.
* `ζ(2)` is Mathlib's `riemannZeta 2 : ℂ`; statement 1 is therefore an equation in `ℂ`, with the
  real quantities coerced.
* The integrals `I(x,y)` are NOT formalized here: statements 2, 3 and 5 take the paper's six
  displayed evaluations of `I(x,y)` as the inputs (written out as the paper's right-hand sides)
  and check the claimed combination. Their content is the `lem:dilog` evaluations (lines 2 and 5
  for the `c_*` bracket, lines 3 and 6 for the rest). Statement 4 is the final "Adding" step,
  pure algebra with `Li_2(b)` as an arbitrary real number.

NOT COVERED: the integral definition of `T` (main.tex 2372) and its evaluation; the integral
representation `M_1 = ∫_b^a log t (c_*/v - 2/(v+a)) dv` and the six evaluations of `I(x,y)`
as integrals; `∫_b^a dv/v = 2L`, `∫_b^a dv/(v+a) = L - λ/2`; the rest of `thm:M2` (the
weight-three part `(T1)`) and the numerical values.
-/
-- BEGIN STATEMENTS
namespace D01

open Real

/-- `a = 1 + √2` (main.tex 874). -/
noncomputable def a : ℝ := 1 + Real.sqrt 2
/-- `b = √2 - 1` (main.tex 874). -/
noncomputable def b : ℝ := Real.sqrt 2 - 1
/-- `L = log a` (main.tex 875). -/
noncomputable def L : ℝ := Real.log a
/-- `λ = log 2` (main.tex 875). -/
noncomputable def lam : ℝ := Real.log 2
/-- `c_* = 1 - 1/√2` (main.tex 2350). -/
noncomputable def cstar : ℝ := 1 - 1 / Real.sqrt 2

/-- The paper's displayed values of `I(b,0)`, `I(-a,0)`, `I(0,0)` (main.tex 2353-2355). -/
noncomputable def I_b_0 : ℝ := D00.Li2 (b ^ 2) - π ^ 2 / 6
noncomputable def I_ma_0 : ℝ := D00.Li2 (-(b ^ 2)) + π ^ 2 / 12 + 2 * L ^ 2
noncomputable def I_0_0 : ℝ := 0
/-- The paper's displayed values of `I(-a,-a)`, `I(b,-a)`, `I(0,-a)` (main.tex 2355-2359). -/
noncomputable def I_ma_ma : ℝ := -(5 / 8) * lam ^ 2 + lam * L + (1 / 2) * L ^ 2
noncomputable def I_b_ma : ℝ :=
  D00.Li2 (2 - Real.sqrt 2) - π ^ 2 / 6 - (5 / 8) * lam ^ 2 + lam * L + (1 / 2) * L ^ 2
noncomputable def I_0_ma : ℝ :=
  -D00.Li2 ((2 + Real.sqrt 2) / 4) + π ^ 2 / 12 - (9 / 8) * lam ^ 2 + lam * L + (1 / 2) * L ^ 2

/-- The closed form of Theorem `thm:M1` (main.tex 891-893). -/
noncomputable def M1closed : ℝ :=
  (19 * Real.sqrt 2 / 48 - 1 / 24) * π ^ 2 + lam ^ 2 / 2 + 2 * Real.sqrt 2 * lam * L
    - (1 + 5 * Real.sqrt 2 / 2) * L ^ 2 - 5 * Real.sqrt 2 * D00.Li2 (Real.sqrt 2 - 1)

/-- eq. `(T)` of `thm:M2`: `(π² + 6λ²)/24 = (ζ(2) - Li₂(1/2))/2`, with `ζ(2)` Mathlib's
`riemannZeta 2` (an equation in `ℂ`). -/
def stmt_1 : Prop :=
  ((((π ^ 2 + 6 * lam ^ 2) / 24 : ℝ)) : ℂ) = (riemannZeta 2 - ((D00.Li2 (1 / 2) : ℝ) : ℂ)) / 2

/-- "By Lemma lem:dilog the bracket multiplying `c_*` is `10 Li₂(b) - 19/24 π² + 5L²`":
`I(b,0) + I(-a,0) - I(0,0) = 10 Li₂(b) - 19π²/24 + 5L²` with the paper's displayed `I` values. -/
def stmt_2 : Prop :=
  I_b_0 + I_ma_0 - I_0_0 = 10 * D00.Li2 b - 19 / 24 * π ^ 2 + 5 * L ^ 2

/-- "... and the rest is `3/4 π² + 1/2 λ² + 2√2 λL - 6L² - 10 Li₂(b)`":
`-2[I(b,-a) + I(-a,-a) - I(0,-a)] - 2λ² + 2√2 λL` equals that, with the paper's `I` values. -/
def stmt_3 : Prop :=
  -2 * (I_b_ma + I_ma_ma - I_0_ma) - 2 * lam ^ 2 + 2 * Real.sqrt 2 * lam * L
    = 3 / 4 * π ^ 2 + 1 / 2 * lam ^ 2 + 2 * Real.sqrt 2 * lam * L - 6 * L ^ 2 - 10 * D00.Li2 b

/-- "Adding gives Theorem thm:M1": `c_*` times the bracket plus the rest is the closed form of
`thm:M1`. -/
def stmt_4 : Prop :=
  cstar * (10 * D00.Li2 b - 19 / 24 * π ^ 2 + 5 * L ^ 2)
    + (3 / 4 * π ^ 2 + 1 / 2 * lam ^ 2 + 2 * Real.sqrt 2 * lam * L - 6 * L ^ 2 - 10 * D00.Li2 b)
    = M1closed

/-- The whole closing computation: the paper's formula
`M_1 = c_*[I(b,0)+I(-a,0)-I(0,0)] - 2[I(b,-a)+I(-a,-a)-I(0,-a)] - 2λ² + 2√2 λL`, evaluated with
the paper's displayed `I` values, equals the closed form of `thm:M1`. -/
def stmt_5 : Prop :=
  cstar * (I_b_0 + I_ma_0 - I_0_0) - 2 * (I_b_ma + I_ma_ma - I_0_ma) - 2 * lam ^ 2
    + 2 * Real.sqrt 2 * lam * L = M1closed

end D01
-- END STATEMENTS

namespace D01

open Real

lemma sq2 : Real.sqrt 2 ^ 2 = 2 := Real.sq_sqrt (by norm_num)

lemma s_gt : 1 < Real.sqrt 2 := by
  have := sq2
  nlinarith [Real.sqrt_nonneg 2]

lemma s_lt : Real.sqrt 2 < 2 := by
  have := sq2
  nlinarith [Real.sqrt_nonneg 2]

lemma ab : a * b = 1 := by
  unfold a b; nlinarith [sq2]

lemma a_pos : 0 < a := by unfold a; positivity

lemma b_pos : 0 < b := by unfold b; linarith [s_gt]

lemma b_lt : b < 1 := by unfold b; linarith [s_lt]

lemma logb : Real.log b = -L := by
  have : b = a⁻¹ := (inv_eq_of_mul_eq_one_right ab).symm
  rw [this, Real.log_inv]; rfl

lemma logs : Real.log (Real.sqrt 2) = lam / 2 := by
  rw [Real.log_sqrt (by norm_num)]; rfl

lemma log_2ms : Real.log (2 - Real.sqrt 2) = lam / 2 - L := by
  have h : 2 - Real.sqrt 2 = Real.sqrt 2 * b := by unfold b; nlinarith [sq2]
  rw [h, Real.log_mul (by linarith [s_gt]) b_pos.ne', logs, logb]; ring

lemma log_a2 : Real.log (a ^ 2) = 2 * L := by
  rw [Real.log_pow]; simp [L]

lemma log_4 : Real.log 4 = 2 * lam := by
  have : (4:ℝ) = 2 ^ 2 := by norm_num
  rw [this, Real.log_pow]; simp [lam]

lemma fa : a = 1 + Real.sqrt 2 := rfl
lemma fb : b = Real.sqrt 2 - 1 := rfl

/-- (A) chi Landen at the fixed point b. -/
lemma relA : D00.Li2 b - D00.Li2 (-b) = π ^ 2 / 8 - L ^ 2 / 2 := by
  have h := D00.chi2_landen b_pos b_lt
  have s0 := s_gt
  have e1 : (1 - b) / (1 + b) = b := by
    rw [div_eq_iff (by linarith [b_pos])]; unfold b; nlinarith [sq2]
  have e2 : (1 + b) / (1 - b) = a := by
    rw [div_eq_iff (by linarith [b_lt])]; unfold a b; nlinarith [sq2]
  rw [e1, e2] at h
  unfold D00.chi2 at h
  have : Real.log a = L := rfl
  rw [this, logb] at h
  linarith

lemma relB : D00.Li2 b + D00.Li2 (-b) = (1 / 2) * D00.Li2 (b ^ 2) :=
  D00.duplication (by rw [abs_of_pos b_pos]; exact b_lt.le)

lemma relC : D00.Li2 (2 - Real.sqrt 2) + D00.Li2 b
    = π ^ 2 / 6 + L * (lam / 2 - L) := by
  have h0 : 0 < 2 - Real.sqrt 2 := by linarith [s_lt]
  have h1 : 2 - Real.sqrt 2 < 1 := by linarith [s_gt]
  have h := D00.reflection h0 h1
  have e : 1 - (2 - Real.sqrt 2) = b := by unfold b; ring
  rw [e, log_2ms, logb] at h
  linarith

lemma relD : D00.Li2 (-Real.sqrt 2) + D00.Li2 (2 - Real.sqrt 2) = -(1 / 2) * L ^ 2 := by
  have h := D00.landen (x := -Real.sqrt 2) (by linarith [s_gt])
  have e : -Real.sqrt 2 / (-Real.sqrt 2 - 1) = 2 - Real.sqrt 2 := by
    rw [div_eq_iff (by linarith [s_gt])]; nlinarith [sq2]
  have e2 : 1 - -Real.sqrt 2 = a := by unfold a; ring
  rw [e, e2] at h
  exact h

lemma relE : -(π ^ 2 / 12) + D00.Li2 (-(b ^ 2)) - D00.Li2 (b ^ 2)
    = D00.Li2 (-Real.sqrt 2) + D00.Li2 (-b) + L * (lam / 2) := by
  have s0 := s_gt
  have hx : -Real.sqrt 2 < 1 := by linarith
  have hy : 1 - Real.sqrt 2 < 1 := by linarith
  have hxy : -Real.sqrt 2 + (1 - Real.sqrt 2) < 1 := by linarith
  have h := D00.five_term hx hy hxy
  have e1 : -Real.sqrt 2 / (1 - (1 - Real.sqrt 2)) = -1 := by
    rw [div_eq_iff (by simp)]; ring
  have e2 : (1 - Real.sqrt 2) / (1 - -Real.sqrt 2) = -(b ^ 2) := by
    rw [div_eq_iff (by linarith)]; unfold b; nlinarith [sq2]
  have e3 : -Real.sqrt 2 * (1 - Real.sqrt 2) / ((1 - -Real.sqrt 2) * (1 - (1 - Real.sqrt 2)))
      = b ^ 2 := by
    rw [div_eq_iff (by nlinarith [sq2])]; unfold b; nlinarith [sq2]
  have e4 : 1 - Real.sqrt 2 = -b := by unfold b; ring
  have e5 : 1 - -Real.sqrt 2 = a := by unfold a; ring
  have e6 : 1 - (1 - Real.sqrt 2) = Real.sqrt 2 := by ring
  rw [e1, e2, e3, e5, e6, logs] at h
  rw [e4] at h
  have : Real.log a = L := rfl
  rw [this] at h
  rw [D00.Li2_neg_one] at h
  linarith

lemma relF : D00.Li2 (-(a ^ 2)) + D00.Li2 (-(b ^ 2)) = -(π ^ 2 / 6) - (1 / 2) * (2 * L) ^ 2 := by
  have h := D00.inversion_neg (x := -(a ^ 2)) (by nlinarith [a_pos])
  have e : 1 / -(a ^ 2) = -(b ^ 2) := by
    have := ab
    have ha := a_pos.ne'
    field_simp
    nlinarith
  rw [e, neg_neg, log_a2] at h
  exact h

lemma relG : D00.Li2 ((2 + Real.sqrt 2) / 4) + D00.Li2 (-(a ^ 2))
    = -(1 / 2) * (lam / 2 - L - 2 * lam) ^ 2 := by
  have s0 := s_gt
  have s1 := s_lt
  have h := D00.landen (x := (2 + Real.sqrt 2) / 4) (by linarith)
  have e : (2 + Real.sqrt 2) / 4 / ((2 + Real.sqrt 2) / 4 - 1) = -(a ^ 2) := by
    rw [div_eq_iff (by linarith)]; unfold a; nlinarith [sq2]
  have e2 : Real.log (1 - (2 + Real.sqrt 2) / 4) = lam / 2 - L - 2 * lam := by
    have : 1 - (2 + Real.sqrt 2) / 4 = (2 - Real.sqrt 2) / 4 := by ring
    rw [this, Real.log_div (by linarith) (by norm_num), log_2ms, log_4]
  rw [e, e2] at h
  exact h

end D01

theorem D01.proof_1 : D01.stmt_1 := by
  unfold D01.stmt_1
  rw [D00.Li2_half, riemannZeta_two]
  unfold D01.lam
  push_cast
  ring

theorem D01.proof_2 : D01.stmt_2 := by
  unfold D01.stmt_2 D01.I_b_0 D01.I_ma_0 D01.I_0_0
  have A := D01.relA
  have B := D01.relB
  have C := D01.relC
  have D := D01.relD
  have E := D01.relE
  linarith

theorem D01.proof_3 : D01.stmt_3 := by
  unfold D01.stmt_3 D01.I_b_ma D01.I_ma_ma D01.I_0_ma
  have A := D01.relA
  have B := D01.relB
  have C := D01.relC
  have D := D01.relD
  have E := D01.relE
  have F := D01.relF
  have G := D01.relG
  nlinarith

theorem D01.proof_4 : D01.stmt_4 := by
  unfold D01.stmt_4 D01.cstar D01.M1closed
  have h : 1 / Real.sqrt 2 = Real.sqrt 2 / 2 := by
    have := D01.s_gt
    field_simp
    nlinarith [D01.sq2]
  rw [h, ← D01.fb]
  ring

theorem D01.proof_5 : D01.stmt_5 := by
  have h2 := D01.proof_2
  have h3 := D01.proof_3
  have h4 := D01.proof_4
  unfold D01.stmt_2 at h2
  unfold D01.stmt_3 at h3
  unfold D01.stmt_4 at h4
  unfold D01.stmt_5
  linear_combination D01.cstar * h2 + h3 + h4
