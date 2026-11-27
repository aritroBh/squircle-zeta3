import Mathlib
import Zeta3.Wave9.R04_SuperellipsoidAllK
import Zeta3.Wave2.D00_Dilog
import Zeta3.Wave2.D02_DilogSqrt2
import Zeta3.Wave3.W03_M1Full
/-!
# S02 — the vertex constant `Q₀` of the superellipsoid in closed form

Target id: S02.
Paper label: `thm:S` (paper/ieee/main.tex lines 1921–1943), the closed form of `Q₀`; the sentence
"Then $s_2$ is \eqref{eq:s2} in Table~\ref{tab:closed}" (main.tex 1944) with `eq:s2`
(main.tex 949–950) and the decimals of `s₂` (main.tex 1944–1947).
Detailed proof: paper/ieee/proofs/sec_g3_S.tex, subsection "The vertex constant $Q_0$"
(lines 123–485): Lemmas `g3:lem:Q0simplex`, `g3:lem:step1`, `g3:lem:J`, `g3:lem:A`, `g3:lem:P`
and "End of the proof of Theorem `g3:thm:S`".

Earlier Lean used here (imported, not restated):
* Z04 (`Z04.Q`, `Z04.h`): `Q(ε) = -(1/(√π Γ(1+2ε))) ∫₀^∞ h_ε(y)³ dy`, so `Q₀ = Z04.Q 0`.
* R04 (`R04.s`): the coefficients `s_k` of `S(p) = Σ s_k p^{-k} + O(p^{-K-1})`; R04 proves
  `s₂ = 48M₁ - 4π² - 8Q(0)`.
* D00 (`D00.Li2`): the real dilogarithm `Li₂(x) = -∫₀ˣ log(1-t)/t dt` (genuine on `(-∞,1]`).
  The `Li₂` arguments below are `b = √2-1 ≈ 0.414`, `-b⁴ ≈ -0.029` and, in `stmt_5`,
  `2cos θ - 1 ∈ [0,1]`; all lie in `(-∞,1]`.
* D02 (`D02.a`, `D02.b`, `D02.L`, `D02.lam`): `a = 1+√2`, `b = √2-1`, `L = log a`, `λ = log 2`
  (main.tex 875, 918, Section `sec:square`).
* W03 (imported for the proof of the corollary only): `M₁ = W02.Mk 1 = D03.M1closed`.

Quoted verbatim (main.tex 1936–1941):

  where $Q_0=Q(0)$ and $Q_1=Q'(0)$ ($Q$ is smooth on $[0,\frac14]$, one-sided at $0$) and
  \begin{multline*}
  Q_0=-\frac{\pi^2}{8}+\frac{\sqrt2}{16}\pi^2+\frac{3\sqrt2}{2}\pi\alpha_c-3\lambda^2-9\sqrt2\,L\log3\\
  +(6+6\sqrt2)L^2-6\sqrt2\,\Li_2(b)+\frac{9\sqrt2}{4}\Li_2(-b^4),
  \end{multline*}
  with $\alpha_c=\arctan\frac1{\sqrt2}$ and $a,b,L,\lambda$ as in Section~\ref{sec:square}.

Quoted verbatim (main.tex 949–950, Table `tab:closed`):

  s_2&=-5\pi^2+\frac{37\sqrt2}{2}\pi^2-12\sqrt2\,\pi\alpha_c+48\lambda^2-96L^2-168\sqrt2\,L^2+96\sqrt2\,\lambda L\notag\\
  &\quad+72\sqrt2\,L\log3-192\sqrt2\,\Li_2(b)-18\sqrt2\,\Li_2(-b^4),\qquad \alpha_c=\arctan\tfrac{1}{\sqrt2}\tag{T6}\label{eq:s2}

Quoted verbatim (main.tex 1944–1946): "Numerically
  \begin{gather*} s_1=-26.5585583865\ldots,\qquad s_2=-4.6892155184\ldots, ..."

Quoted verbatim (sec_g3_S.tex 127–129, Lemma `g3:lem:Q0simplex`):

  Let $D(x)=x_1+x_2+x_3-\sum_{i<j}\sqrt{x_i^2+x_j^2}+\sqrt{x_1^2+x_2^2+x_3^2}$. Then $D\ge0$ and
  $Q_0=-\int_\Delta\frac{D(\sigma)}{\sigma_1\sigma_2\sigma_3}\dd\sigma$, that is, with
  $\rho=\sqrt{\sigma_1^2+\sigma_2^2+\sigma_3^2}$,
  $Q_0=\int_\Delta\frac{\sum_{i<j}\sqrt{\sigma_i^2+\sigma_j^2}-1-\rho}{\sigma_1\sigma_2\sigma_3}\dd\sigma$.

Quoted verbatim (sec_g3_S.tex 137–145, Lemma `g3:lem:step1`):

  Let $q(s)=\sqrt{s^2+(1-s)^2}$, $s\in[0,1]$ ...
  $Q_0=6\int_0^1J\bigl(q(s)\bigr)\dd s$, where for $q\in[\frac1{\sqrt2},1]$
  \[ J(q)=\int_0^\infty\Bigl[\frac1{1+r+\sqrt{q^2+r^2}}-\frac1{1+q}\Bigr]\frac{\dd r}{r(1+r)} . \]

Quoted verbatim (sec_g3_S.tex 168): Put
  $P_1=\int_0^1\frac{2\log(1+q)-(1+q)\log(2q)}{1-q^2}\dd s$,
  $P_{2a}=\int_0^1\frac{\arsinh q}{\sqrt{1+q^2}}\dd s$,
  $P_{2c}=\int_0^1\frac{\arsinh(1/q)}{\sqrt{1+q^2}}\dd s$ ($q=q(s)$). The inner integral $J$ is
  elementary (Lemma~\ref{g3:lem:J}), and so $Q_0=6(P_1-P_{2a}-P_{2c})$.

Quoted verbatim (sec_g3_S.tex 170–174, Lemma `g3:lem:J`):

  For $0<q<1$, with $\beta=\sqrt{1+q^2}$,
  \[ J(q)=\frac{2\log(1+q)-(1+q)\log(2q)}{1-q^2}-\frac{\operatorname{arsinh}q+\operatorname{arsinh}(1/q)}{\beta}. \]

Quoted verbatim (sec_g3_S.tex 211–213, Lemma `g3:lem:A`):

  For $0\le\theta\le\pi/3$ and $\zeta_\theta=(2\cos\theta-1)e^{i\theta}$:
  $\operatorname{Re}\Li_2(\zeta_\theta)=\frac12\Li_2(2\cos\theta-1)+\frac34\theta^2-\frac\pi2\theta+\frac{\pi^2}{12}$.
  In particular $\operatorname{Re}\Li_2\bigl(b\,e^{\pm i\pi/4}\bigr)=\frac12\Li_2(b)+\frac{\pi^2}{192}$
  (the case $\theta=\pi/4$, $2\cos\frac\pi4-1=b$).

Quoted verbatim (sec_g3_S.tex 219–226, Lemma `g3:lem:P`):

  With $\lambda=\log2$, $L=\log(1+\sqrt2)$, $b=\sqrt2-1$, $\alpha_c=\arctan\frac1{\sqrt2}$:
  \begin{align*}
  P_1&=3\sqrt2\Li_2(b)-\tfrac{\lambda^2}2-\tfrac{\sqrt2}2\lambda L+\bigl(1+\tfrac{3\sqrt2}2\bigr)L^2-\bigl(\tfrac1{48}+\tfrac{3\sqrt2}{16}\bigr)\pi^2,\\
  P_{2a}&=\sqrt2\Bigl[3\Li_2(b)-\tfrac14\Li_2(-b^4)-\tfrac{\lambda L}2+L\log3+\tfrac{L^2}2-\tfrac{5\pi^2}{24}\Bigr],\\
  P_{2c}&=\sqrt2\Bigl[\Li_2(b)-\tfrac18\Li_2(-b^4)+\tfrac{L\log3}2+\tfrac{\pi^2}{96}-\tfrac{\pi\alpha_c}4\Bigr].
  \end{align*}

Conventions.
* `Δ` is `Y05.Δ = {t ∈ ℝ² | t₁ ≥ 0, t₂ ≥ 0, t₁+t₂ ≤ 1}` with `σ = (t₁, t₂, 1-t₁-t₂)` and
  `dσ = dt₁ dt₂` (Lebesgue measure on `ℝ²`); this is the normalisation of the supplement's
  Dirichlet integrals (`∫_Δ Π^{δ-1} dσ = Γ(δ)³/Γ(3δ)`). On the boundary of `Δ` the integrand
  has the junk value `x/0 = 0`; the boundary is a null set.
* The integrand of `P₁` has a removable singularity at `q = 1` (`s ∈ {0,1}`); Lean's value
  there is `x/0 = 0`, on a null set. Interval integrals are stated with `IntervalIntegrable`.
* `Re Li₂(z)` for `|z| ≤ 1` is written by the power series: `Re Σ zⁿ/n² = Σ rⁿ cos(nθ)/n²`
  for `z = r e^{iθ}`, `0 ≤ r ≤ 1` (the principal dilogarithm equals its power series on the
  closed unit disc). In Lemma A, `r = 2cos θ - 1 ∈ [0,1]` for `θ ∈ [0, π/3]`. Since cosine is
  even, the case `b e^{-iπ/4}` has the same series as `b e^{iπ/4}`.
* `D ≥ 0` is stated for `x ∈ [0,∞)³` (the paper prints `D ≥ 0` with no domain; for
  `x = (-1,0,0)`, `D = -2`. The paper applies it only on `Δ`).
* The integrability of `h_0³` (so that `Z04.Q 0` is a genuine integral) is R04 `stmt_6`; not
  restated.

Numerical check (mpmath, 34–40 digits, before freezing): `Q₀` from the definition
`-(1/√π)∫h₀³` is `-2.832935683494847760462057251515769`, equal to the closed form to 33 digits;
the simplex integral of stmt_1 (in the coordinates `(s,r)` of the proof) agrees to 25 digits;
`J` against its closed form at `q = 0.72, 0.85, 0.97`: `≤ 3·10⁻³⁹`; `P₁, P_{2a}, P_{2c}` against
their closed forms: `≤ 2·10⁻⁴⁰`; `6(P₁-P_{2a}-P_{2c})` and `6∫J(q(s))ds` against the closed form:
`5·10⁻⁴¹`; Lemma A at `θ = 0.01, 0.3, π/4, 0.9`: `≤ 3·10⁻⁴¹`; `48M₁-4π²-8Q₀` (with the closed
form of `thm:M1`) against (T6): `7·10⁻³⁹`, value `-4.689215518460502249960...`.

NOT COVERED:
* The second form of the simplex integrand in `g3:lem:Q0simplex` (it is the first one with
  `σ₁+σ₂+σ₃ = 1`, i.e. the same function on `Δ`).
* Inner steps of the proofs: the Gaussian representation of `D`, the regularised integral
  `I(δ)` and its Dirichlet evaluation, the substitution (`g3:eq:sq`), the formula
  `P_{2a} = √2(L log 3/2 - W)` and the integral `W`, the complex `χ₂` identities (a)–(e) for
  `P_{2c}`, the values `Re Li₂(Z₁), Re Li₂(Z₂), Re Li₂(Z₃)`, the eighteen rows of Table
  `g3:tab:P1` and the fifteen dilogarithm instances of Table `g3:tab:dilog` (several of them are
  already D02 statements).
* The claim "fifteen instances of the classical dilogarithm identities" as a count, and the
  computer-algebra remarks.
* `Q₁` and its decimals; `s₁`, `s₃` decimals; the `n`-dimensional `S_n(p)` formula.

Statements 1, 2, 4, 6, 7, 8, 9, 10, 11 moved to `Zeta3.Wave10.S06_VertexQ0Rest` (continuation of S02); this file keeps 3 and 5.
-/
-- BEGIN STATEMENTS
namespace S02

open Real MeasureTheory Set

/-- The cube angle `α_c = arctan(1/√2)` (main.tex 1941). -/
noncomputable def alphac : ℝ := Real.arctan (1 / Real.sqrt 2)

/-- `q(s) = √(s² + (1-s)²)` (sec_g3_S.tex 137). -/
noncomputable def qs (s : ℝ) : ℝ := Real.sqrt (s ^ 2 + (1 - s) ^ 2)

/-- `D(x) = x₁+x₂+x₃ - Σ_{i<j} √(x_i²+x_j²) + √(x₁²+x₂²+x₃²)` (`g3:lem:Q0simplex`). -/
noncomputable def Dfun (x₁ x₂ x₃ : ℝ) : ℝ :=
  x₁ + x₂ + x₃ - Real.sqrt (x₁ ^ 2 + x₂ ^ 2) - Real.sqrt (x₂ ^ 2 + x₃ ^ 2)
    - Real.sqrt (x₁ ^ 2 + x₃ ^ 2) + Real.sqrt (x₁ ^ 2 + x₂ ^ 2 + x₃ ^ 2)

/-- `D(σ)/(σ₁σ₂σ₃)` on the simplex, with `σ = (t₁, t₂, 1-t₁-t₂)`. -/
noncomputable def simplexIntegrand (t : ℝ × ℝ) : ℝ :=
  Dfun t.1 t.2 (1 - t.1 - t.2) / (t.1 * t.2 * (1 - t.1 - t.2))

/-- The integrand of `J(q)`: `[1/(1+r+√(q²+r²)) - 1/(1+q)] / (r(1+r))`. -/
noncomputable def Jintegrand (q r : ℝ) : ℝ :=
  (1 / (1 + r + Real.sqrt (q ^ 2 + r ^ 2)) - 1 / (1 + q)) / (r * (1 + r))

/-- `J(q) = ∫₀^∞ [1/(1+r+√(q²+r²)) - 1/(1+q)] dr/(r(1+r))` (`g3:lem:step1`). -/
noncomputable def J (q : ℝ) : ℝ := ∫ r in Ioi (0 : ℝ), Jintegrand q r

/-- `(2log(1+q) - (1+q)log(2q))/(1-q²)`, the integrand of `P₁` and the first part of `J`. -/
noncomputable def f1 (q : ℝ) : ℝ :=
  (2 * Real.log (1 + q) - (1 + q) * Real.log (2 * q)) / (1 - q ^ 2)

/-- `P₁ = ∫₀¹ (2log(1+q)-(1+q)log(2q))/(1-q²) ds`, `q = q(s)`. -/
noncomputable def P1 : ℝ := ∫ s in (0 : ℝ)..1, f1 (qs s)

/-- `P_{2a} = ∫₀¹ arsinh q / √(1+q²) ds`, `q = q(s)`. -/
noncomputable def P2a : ℝ :=
  ∫ s in (0 : ℝ)..1, Real.arsinh (qs s) / Real.sqrt (1 + qs s ^ 2)

/-- `P_{2c} = ∫₀¹ arsinh(1/q) / √(1+q²) ds`, `q = q(s)`. -/
noncomputable def P2c : ℝ :=
  ∫ s in (0 : ℝ)..1, Real.arsinh (1 / qs s) / Real.sqrt (1 + qs s ^ 2)

/-- The closed form of `Q₀` printed in `thm:S`. -/
noncomputable def Q0closed : ℝ :=
  -(π ^ 2 / 8) + Real.sqrt 2 / 16 * π ^ 2 + 3 * Real.sqrt 2 / 2 * π * alphac
    - 3 * D02.lam ^ 2 - 9 * Real.sqrt 2 * D02.L * Real.log 3
    + (6 + 6 * Real.sqrt 2) * D02.L ^ 2 - 6 * Real.sqrt 2 * D00.Li2 D02.b
    + 9 * Real.sqrt 2 / 4 * D00.Li2 (-D02.b ^ 4)

/-- The closed form (T6) of `s₂` printed in Table `tab:closed` (`eq:s2`). -/
noncomputable def s2closed : ℝ :=
  -5 * π ^ 2 + 37 * Real.sqrt 2 / 2 * π ^ 2 - 12 * Real.sqrt 2 * π * alphac
    + 48 * D02.lam ^ 2 - 96 * D02.L ^ 2 - 168 * Real.sqrt 2 * D02.L ^ 2
    + 96 * Real.sqrt 2 * D02.lam * D02.L + 72 * Real.sqrt 2 * D02.L * Real.log 3
    - 192 * Real.sqrt 2 * D00.Li2 D02.b - 18 * Real.sqrt 2 * D00.Li2 (-D02.b ^ 4)

/-- Lemma `g3:lem:J`: for `0 < q < 1`, with `β = √(1+q²)`,
`J(q) = (2log(1+q)-(1+q)log(2q))/(1-q²) - (arsinh q + arsinh(1/q))/β`. -/
def stmt_3 : Prop :=
  ∀ q : ℝ, 0 < q → q < 1 →
    J q = f1 q - (Real.arsinh q + Real.arsinh (1 / q)) / Real.sqrt (1 + q ^ 2)

/-- Lemma `g3:lem:A` (Lemma A): for `0 ≤ θ ≤ π/3`, with `r = 2cos θ - 1`,
`Re Li₂(r e^{iθ}) = Σ_{n≥1} rⁿ cos(nθ)/n² = ½Li₂(r) + ¾θ² - (π/2)θ + π²/12`; and the case
`θ = π/4` (`r = b`): `Re Li₂(b e^{±iπ/4}) = ½Li₂(b) + π²/192`. -/
def stmt_5 : Prop :=
  (∀ θ : ℝ, 0 ≤ θ → θ ≤ π / 3 →
    HasSum (fun n : ℕ => (2 * Real.cos θ - 1) ^ (n + 1) * Real.cos ((n + 1 : ℝ) * θ)
        / ((n : ℝ) + 1) ^ 2)
      ((1 / 2) * D00.Li2 (2 * Real.cos θ - 1) + (3 / 4) * θ ^ 2 - π / 2 * θ + π ^ 2 / 12)) ∧
  HasSum (fun n : ℕ => D02.b ^ (n + 1) * Real.cos ((n + 1 : ℝ) * (π / 4)) / ((n : ℝ) + 1) ^ 2)
    ((1 / 2) * D00.Li2 D02.b + π ^ 2 / 192)

end S02
-- END STATEMENTS

namespace S02

open Real MeasureTheory Set Filter Topology


lemma hlawka3 (x₁ x₂ x₃ : ℝ) (h1 : 0 ≤ x₁) (h2 : 0 ≤ x₂) (h3 : 0 ≤ x₃) :
    0 ≤ x₁ + x₂ + x₃ - Real.sqrt (x₁ ^ 2 + x₂ ^ 2) - Real.sqrt (x₂ ^ 2 + x₃ ^ 2)
    - Real.sqrt (x₁ ^ 2 + x₃ ^ 2) + Real.sqrt (x₁ ^ 2 + x₂ ^ 2 + x₃ ^ 2) := by
  set a := Real.sqrt (x₁ ^ 2 + x₂ ^ 2) with ha
  set b := Real.sqrt (x₂ ^ 2 + x₃ ^ 2) with hb
  set c := Real.sqrt (x₁ ^ 2 + x₃ ^ 2) with hc
  set r := Real.sqrt (x₁ ^ 2 + x₂ ^ 2 + x₃ ^ 2) with hr
  have a2 : a ^ 2 = x₁ ^ 2 + x₂ ^ 2 := Real.sq_sqrt (by positivity)
  have b2 : b ^ 2 = x₂ ^ 2 + x₃ ^ 2 := Real.sq_sqrt (by positivity)
  have c2 : c ^ 2 = x₁ ^ 2 + x₃ ^ 2 := Real.sq_sqrt (by positivity)
  have r2 : r ^ 2 = x₁ ^ 2 + x₂ ^ 2 + x₃ ^ 2 := Real.sq_sqrt (by positivity)
  have a0 : 0 ≤ a := Real.sqrt_nonneg _
  have b0 : 0 ≤ b := Real.sqrt_nonneg _
  have c0 : 0 ≤ c := Real.sqrt_nonneg _
  have r0 : 0 ≤ r := Real.sqrt_nonneg _
  have ta : a ≤ x₁ + x₂ := by
    rw [ha]; apply Real.sqrt_le_iff.mpr; constructor
    · positivity
    · nlinarith [mul_nonneg h1 h2]
  have tb : b ≤ x₂ + x₃ := by
    rw [hb]; apply Real.sqrt_le_iff.mpr; constructor
    · positivity
    · nlinarith [mul_nonneg h2 h3]
  have tc : c ≤ x₁ + x₃ := by
    rw [hc]; apply Real.sqrt_le_iff.mpr; constructor
    · positivity
    · nlinarith [mul_nonneg h1 h3]
  have ra : a ≤ r := by rw [ha, hr]; apply Real.sqrt_le_sqrt; nlinarith
  have rb : b ≤ r := by rw [hb, hr]; apply Real.sqrt_le_sqrt; nlinarith
  have rc : c ≤ r := by rw [hc, hr]; apply Real.sqrt_le_sqrt; nlinarith
  have key : (x₁ + x₂ + x₃ - a - b - c + r) * (x₁ + x₂ + x₃ + r)
      = (x₁ + x₂ - a) * (x₃ - a + r) + (x₂ + x₃ - b) * (x₁ - b + r)
        + (x₁ + x₃ - c) * (x₂ - c + r) := by
    linear_combination r2 - a2 - b2 - c2
  have hpos : 0 ≤ (x₁ + x₂ + x₃ - a - b - c + r) * (x₁ + x₂ + x₃ + r) := by
    rw [key]
    have := mul_nonneg (sub_nonneg.mpr ta) (show 0 ≤ x₃ - a + r by linarith)
    have := mul_nonneg (sub_nonneg.mpr tb) (show 0 ≤ x₁ - b + r by linarith)
    have := mul_nonneg (sub_nonneg.mpr tc) (show 0 ≤ x₂ - c + r by linarith)
    linarith
  rcases eq_or_lt_of_le (show 0 ≤ x₁ + x₂ + x₃ + r by linarith) with h | h
  · have e1 : x₁ = 0 := by linarith
    have e2 : x₂ = 0 := by linarith
    have e3 : x₃ = 0 := by linarith
    have : r = 0 := by linarith
    have : a = 0 := by linarith
    have : b = 0 := by linarith
    have : c = 0 := by linarith
    linarith
  · exact nonneg_of_mul_nonneg_left hpos h


lemma tail_abs {f : ℕ → ℝ} {S : ℝ} (hS : HasSum f S) (N : ℕ) (K r : ℝ) (hr0 : 0 ≤ r)
    (hr1 : r < 1) (hbd : ∀ n, |f (n + N)| ≤ K * r ^ n) :
    |S - ∑ n ∈ Finset.range N, f n| ≤ K / (1 - r) := by
  have h1 := (hasSum_nat_add_iff' N).mpr hS
  have h2 : HasSum (fun n : ℕ => K * r ^ n) (K * (1 - r)⁻¹) :=
    (hasSum_geometric_of_lt_one hr0 hr1).mul_left K
  have hs : Summable (fun n => |f (n + N)|) :=
    Summable.of_nonneg_of_le (fun n => abs_nonneg _) hbd h2.summable
  have e1 : |∑' n, f (n + N)| ≤ ∑' n, |f (n + N)| := by
    simpa [Real.norm_eq_abs] using norm_tsum_le_tsum_norm (f := fun n => f (n + N)) (by simpa [Real.norm_eq_abs] using hs)
  have e2 : ∑' n, |f (n + N)| ≤ K * (1 - r)⁻¹ := by
    rw [← h2.tsum_eq]; exact hs.tsum_le_tsum hbd h2.summable
  rw [← h1.tsum_eq, div_eq_mul_inv]
  linarith

lemma alpha_bounds : (0.615479708670387341 : ℝ) ≤ Real.arctan (1 / Real.sqrt 2) ∧
    Real.arctan (1 / Real.sqrt 2) ≤ (0.615479708670387342 : ℝ) := by
  have hs := D03.s_bounds
  set x := 1 / Real.sqrt 2 with hx
  have hx2 : x ^ 2 = 1 / 2 := by
    rw [hx, div_pow, Real.sq_sqrt (by norm_num)]; norm_num
  have hx0 : 0 < x := by positivity
  have hx1 : x < 1 := by
    rw [hx, div_lt_one (by positivity)]; linarith
  have hS := Real.hasSum_arctan (x := x) (by rw [Real.norm_eq_abs, abs_lt]; constructor <;> linarith)
  have hpow : ∀ n : ℕ, x ^ (2 * n + 1) = x * (1 / 2) ^ n := by
    intro n; rw [pow_succ, pow_mul, hx2]; ring
  have ht := tail_abs hS 70 (x * (1 / 2) ^ 70) (1 / 2) (by norm_num) (by norm_num) (fun n => by
    rw [abs_div, abs_mul, abs_pow, abs_neg, abs_one, one_pow, one_mul, abs_of_pos (by positivity : (0:ℝ) < x ^ (2 * (n + 70) + 1)),
      Nat.abs_cast, hpow]
    rw [div_le_iff₀ (by positivity)]
    have : (1 : ℝ) ≤ ((2 * (n + 70) + 1 : ℕ) : ℝ) := by norm_cast; omega
    calc x * (1 / 2) ^ (n + 70) = x * (1 / 2) ^ 70 * (1 / 2) ^ n * 1 := by ring
      _ ≤ x * (1 / 2) ^ 70 * (1 / 2) ^ n * ((2 * (n + 70) + 1 : ℕ) : ℝ) := by gcongr)
  have hsum : ∑ n ∈ Finset.range 70, (-1) ^ n * x ^ (2 * n + 1) / ((2 * n + 1 : ℕ) : ℝ)
      = x * ∑ n ∈ Finset.range 70, (-1) ^ n * (1 / 2 : ℝ) ^ n / ((2 * n + 1 : ℕ) : ℝ) := by
    rw [Finset.mul_sum]; refine Finset.sum_congr rfl (fun n _ => ?_); rw [hpow]; ring
  rw [hsum] at ht
  have hxl : (0.7071067811865475244008443 : ℝ) ≤ x := by
    rw [hx, le_div_iff₀ (by positivity)]; nlinarith
  have hxh : x ≤ (0.7071067811865475244008444 : ℝ) := by
    rw [hx, div_le_iff₀ (by positivity)]; nlinarith
  set T := ∑ n ∈ Finset.range 70, (-1) ^ n * (1 / 2 : ℝ) ^ n / ((2 * n + 1 : ℕ) : ℝ) with hTdef
  have hTv : (0.870419751367103197473551 : ℝ) ≤ T ∧ T ≤ (0.870419751367103197473552 : ℝ) := by
    rw [hTdef]; constructor <;> norm_num [Finset.sum_range_succ]
  rw [abs_le] at ht
  have p1 : x * T ≤ (0.7071067811865475244008444 : ℝ) * (0.870419751367103197473552 : ℝ) :=
    mul_le_mul hxh hTv.2 (le_trans (by norm_num) hTv.1) (by norm_num)
  have p2 : (0.7071067811865475244008443 : ℝ) * (0.870419751367103197473551 : ℝ) ≤ x * T :=
    mul_le_mul hxl hTv.1 (by norm_num) (by linarith)
  have p3 : x * (1 / 2) ^ 70 / (1 - 1 / 2) ≤ 1 / 10 ^ 20 := by
    rw [div_le_iff₀ (by norm_num)]; nlinarith
  constructor <;> norm_num at p1 p2 ⊢ <;> linarith [ht.1, ht.2]


lemma l3_bounds : (1.09861228866810969139 : ℝ) ≤ Real.log 3 ∧ Real.log 3 ≤ (1.0986122886681096914 : ℝ) := by
  have h := D03.art_bounds (x := 1/2) (lo := 1/2) (hi := 1/2) (by norm_num) le_rfl le_rfl (by norm_num) 36
  have e : Real.log (1 + 1/2) - Real.log (1 - 1/2) = Real.log 3 := by
    rw [← Real.log_div (by norm_num) (by norm_num)]; norm_num
  rw [e] at h
  constructor
  · refine le_trans ?_ h.1
    norm_num [Finset.sum_range_succ]
  · refine le_trans h.2 ?_
    norm_num [Finset.sum_range_succ]

lemma b4_eq : D02.b ^ 4 = 17 - 12 * Real.sqrt 2 := by
  have h := Real.sq_sqrt (show (0:ℝ) ≤ 2 by norm_num)
  unfold D02.b
  have : (Real.sqrt 2 - 1) ^ 4 = 17 - 12 * Real.sqrt 2 + 0 * (Real.sqrt 2 ^ 2 - 2) := by
    have e : Real.sqrt 2 ^ 4 = (Real.sqrt 2 ^ 2) ^ 2 := by ring
    have e3 : Real.sqrt 2 ^ 3 = Real.sqrt 2 ^ 2 * Real.sqrt 2 := by ring
    ring_nf; rw [e3, e, h]; ring
  rw [this, h]; ring

lemma m_bounds : (-0.0292234018355678025 : ℝ) ≤ D00.Li2 (-D02.b ^ 4) ∧
    D00.Li2 (-D02.b ^ 4) ≤ (-0.0292234018355678024 : ℝ) := by
  have hs := D03.s_bounds
  have hb4 := b4_eq
  have hb8 : (D02.b ^ 4) ^ 2 = 577 - 408 * Real.sqrt 2 := by
    rw [hb4]; ring_nf; rw [Real.sq_sqrt (by norm_num)]; ring
  have hd := D00.duplication' (x := D02.b ^ 4) (by rw [hb4, abs_le]; constructor <;> linarith)
  have h1 := D03.li2_bounds (x := D02.b ^ 4) (lo := (0.02943725152285941437973 : ℝ))
    (hi := (0.02943725152285941437974 : ℝ)) (by norm_num) (by rw [hb4]; linarith) (by rw [hb4]; linarith)
    (by norm_num) 14
  have h2 := D03.li2_bounds (x := (D02.b ^ 4) ^ 2) (lo := (0.000866551777220088910998 : ℝ))
    (hi := (0.000866551777220088911003 : ℝ)) (by norm_num) (by rw [hb8]; linarith) (by rw [hb8]; linarith)
    (by norm_num) 8
  have a1 : D00.Li2 (D02.b ^ 4) ≤ (0.0296567716243435315273 : ℝ) := by
    refine le_trans h1.2 ?_; norm_num [Finset.sum_range_succ]
  have a2 : (0.0296567716243435315271 : ℝ) ≤ D00.Li2 (D02.b ^ 4) := by
    refine le_trans ?_ h1.1; norm_num [Finset.sum_range_succ]
  have a3 : D00.Li2 ((D02.b ^ 4) ^ 2) ≤ (0.000866739577551458099 : ℝ) := by
    refine le_trans h2.2 ?_; norm_num [Finset.sum_range_succ]
  have a4 : (0.000866739577551458097 : ℝ) ≤ D00.Li2 ((D02.b ^ 4) ^ 2) := by
    refine le_trans ?_ h2.1; norm_num [Finset.sum_range_succ]
  constructor <;> linarith


lemma mb {x y a b c d : ℝ} (hx : a ≤ x ∧ x ≤ b) (hy : c ≤ y ∧ y ≤ d) (ha : 0 ≤ a) (hc : 0 ≤ c) :
    a * c ≤ x * y ∧ x * y ≤ b * d :=
  ⟨mul_le_mul hx.1 hy.1 hc (ha.trans hx.1), mul_le_mul hx.2 hy.2 (hc.trans hy.1) (ha.trans (hx.1.trans hx.2))⟩

lemma s2_num (s p al l L l3 u m : ℝ)
    (hs : (1.41421356237309504880168872 : ℝ) ≤ s ∧ s ≤ (1.41421356237309504880168873 : ℝ))
    (hp : (3.14159265358979323846 : ℝ) ≤ p ∧ p ≤ (3.14159265358979323847 : ℝ))
    (ha : (0.615479708670387341 : ℝ) ≤ al ∧ al ≤ (0.615479708670387342 : ℝ))
    (hl : (0.69314718055994530941723212 : ℝ) ≤ l ∧ l ≤ (0.69314718055994530941723213 : ℝ))
    (hL : (0.88137358701954302523260931 : ℝ) ≤ L ∧ L ≤ (0.88137358701954302523260934 : ℝ))
    (h3 : (1.09861228866810969139 : ℝ) ≤ l3 ∧ l3 ≤ (1.0986122886681096914 : ℝ))
    (hu : (0.46753399702300461762696479 : ℝ) ≤ u ∧ u ≤ (0.46753399702300461762696481 : ℝ))
    (hm : (-0.0292234018355678025 : ℝ) ≤ m ∧ m ≤ (-0.0292234018355678024 : ℝ)) (v : ℝ)
    (hv : v = -5 * p ^ 2 + 37 * s / 2 * p ^ 2 - 12 * s * p * al
      + 48 * l ^ 2 - 96 * L ^ 2 - 168 * s * L ^ 2
      + 96 * s * l * L + 72 * s * L * l3 - 192 * s * u - 18 * s * m) :
    -4.6892155185 < v ∧ v ≤ -4.6892155184 := by
  have hm' : (0.0292234018355678024 : ℝ) ≤ -m ∧ -m ≤ (0.0292234018355678025 : ℝ) := ⟨by linarith [hm.2], by linarith [hm.1]⟩
  have pp := mb hp hp (by norm_num) (by norm_num)
  have spp := mb hs pp (by norm_num) (by norm_num)
  have sp := mb hs hp (by norm_num) (by norm_num)
  have spa := mb sp ha (by norm_num) (by norm_num)
  have ll := mb hl hl (by norm_num) (by norm_num)
  have LL := mb hL hL (by norm_num) (by norm_num)
  have sLL := mb hs LL (by norm_num) (by norm_num)
  have sl := mb hs hl (by norm_num) (by norm_num)
  have slL := mb sl hL (by norm_num) (by norm_num)
  have sL := mb hs hL (by norm_num) (by norm_num)
  have sL3 := mb sL h3 (by norm_num) (by norm_num)
  have su := mb hs hu (by norm_num) (by norm_num)
  have sm := mb hs hm' (by norm_num) (by norm_num)
  have ev : v = -5 * (p * p) + 37 / 2 * (s * (p * p)) - 12 * (s * p * al) + 48 * (l * l)
      - 96 * (L * L) - 168 * (s * (L * L)) + 96 * (s * l * L) + 72 * (s * L * l3)
      - 192 * (s * u) + 18 * (s * -m) := by rw [hv]; ring
  rw [ev]
  generalize p * p = P at pp spp
  generalize s * (P) = SP at spp
  generalize s * p * al = SPA at spa
  generalize l * l = LL2 at ll
  generalize L * L = L2 at LL sLL
  generalize s * L2 = SL2 at sLL
  generalize s * l * L = SLL at slL
  generalize s * L * l3 = SL3 at sL3
  generalize s * u = SU at su
  generalize s * -m = SM at sm
  norm_num at pp spp spa ll LL sLL slL sL3 su sm ⊢
  constructor <;> linarith [pp.1, pp.2, spp.1, spp.2, spa.1, spa.2, ll.1, ll.2, LL.1, LL.2, sLL.1, sLL.2,
    slL.1, slL.2, sL3.1, sL3.2, su.1, su.2, sm.1, sm.2]


/-- `u(r) = r + √(q²+r²)`. -/
noncomputable def uf (q r : ℝ) : ℝ := r + Real.sqrt (q ^ 2 + r ^ 2)

lemma uf_ge {q r : ℝ} (_hq : 0 ≤ q) (hr : 0 ≤ r) : q ≤ uf q r := by
  unfold uf
  have : q ≤ Real.sqrt (q ^ 2 + r ^ 2) := Real.le_sqrt_of_sq_le (by nlinarith)
  linarith

lemma uf_ge_r {q r : ℝ} : r ≤ uf q r := by
  unfold uf; linarith [Real.sqrt_nonneg (q ^ 2 + r ^ 2)]

lemma hasDerivAt_uf {q : ℝ} (hq : 0 < q) (r : ℝ) :
    HasDerivAt (uf q) (1 + r / Real.sqrt (q ^ 2 + r ^ 2)) r := by
  have hpos : 0 < q ^ 2 + r ^ 2 := by positivity
  have h1 : HasDerivAt (fun r => q ^ 2 + r ^ 2) (2 * r) r := by
    simpa using ((hasDerivAt_id r).pow 2).const_add (q ^ 2)
  have h2 := h1.sqrt hpos.ne'
  have h3 := (hasDerivAt_id' r).add h2
  show HasDerivAt (fun r => r + Real.sqrt (q ^ 2 + r ^ 2)) _ r
  convert h3 using 1
  rw [mul_div_mul_left _ _ (by norm_num : (2:ℝ) ≠ 0)]

/-- The four residues. -/
noncomputable def cA (q : ℝ) : ℝ := 2 / (q ^ 2 - 1)
noncomputable def cB (q : ℝ) : ℝ := -(2 * q) / (q ^ 2 - 1)
noncomputable def cC (q : ℝ) : ℝ := 1 / (1 + q) - 1 / Real.sqrt (1 + q ^ 2)
noncomputable def cD (q : ℝ) : ℝ := 1 / Real.sqrt (1 + q ^ 2) + 1 / (1 + q)

/-- The primitive `F(u) = Σ c_k log(u - u_k)`. -/
noncomputable def Fp (q u : ℝ) : ℝ :=
  cA q * Real.log (u + 1) + cB q * Real.log (u + q)
    + cC q * Real.log (u + 1 - Real.sqrt (1 + q ^ 2))
    + cD q * Real.log (u + 1 + Real.sqrt (1 + q ^ 2))

lemma tendsto_log_one_add (d : ℝ) :
    Tendsto (fun u : ℝ => Real.log (1 + d / u)) atTop (𝓝 0) := by
  have h1 : Tendsto (fun u : ℝ => 1 + d / u) atTop (𝓝 1) := by
    simpa using (tendsto_const_nhds (x := (1:ℝ))).add
      (Filter.Tendsto.div_atTop (tendsto_const_nhds (x := d)) tendsto_id)
  rw [← Real.log_one]; exact ((Real.continuousAt_log one_ne_zero).tendsto).comp h1

section
variable {q : ℝ} (hq : 0 < q) (hq1 : q < 1)
include hq hq1

lemma beta_facts : 1 < Real.sqrt (1 + q ^ 2) ∧ Real.sqrt (1 + q ^ 2) < 1 + q ∧
    Real.sqrt (1 + q ^ 2) ^ 2 = 1 + q ^ 2 := by
  have h2 : Real.sqrt (1 + q ^ 2) ^ 2 = 1 + q ^ 2 := Real.sq_sqrt (by positivity)
  refine ⟨?_, ?_, h2⟩
  · rw [Real.lt_sqrt (by norm_num)]; nlinarith
  · rw [Real.sqrt_lt' (by linarith)]; nlinarith

lemma args_pos {u : ℝ} (hu : q ≤ u) :
    0 < u + 1 ∧ 0 < u + q ∧ 0 < u + 1 - Real.sqrt (1 + q ^ 2) ∧
      0 < u + 1 + Real.sqrt (1 + q ^ 2) := by
  obtain ⟨b1, b2, -⟩ := beta_facts hq hq1
  refine ⟨by linarith, by linarith, by linarith, by linarith⟩

lemma hasDerivAt_Fu {r : ℝ} (hr : 0 ≤ r) :
    HasDerivAt (fun r => Fp q (uf q r))
      ((cA q / (uf q r + 1) + cB q / (uf q r + q) + cC q / (uf q r + 1 - Real.sqrt (1 + q ^ 2))
        + cD q / (uf q r + 1 + Real.sqrt (1 + q ^ 2))) * (1 + r / Real.sqrt (q ^ 2 + r ^ 2))) r := by
  have h := hasDerivAt_uf hq r
  obtain ⟨p1, p2, p3, p4⟩ := args_pos hq hq1 (uf_ge hq.le hr)
  have e1 := ((h.add_const 1).log p1.ne').const_mul (cA q)
  have e2 := ((h.add_const q).log p2.ne').const_mul (cB q)
  have e3 := (((h.add_const 1).sub_const (Real.sqrt (1 + q ^ 2))).log p3.ne').const_mul (cC q)
  have e4 := (((h.add_const 1).add_const (Real.sqrt (1 + q ^ 2))).log p4.ne').const_mul (cD q)
  convert ((e1.add e2).add e3).add e4 using 1
  · funext y; simp only [Fp, Pi.add_apply]
  · ring

lemma Dv_eq {r : ℝ} (hr : 0 < r) :
    (cA q / (uf q r + 1) + cB q / (uf q r + q) + cC q / (uf q r + 1 - Real.sqrt (1 + q ^ 2))
        + cD q / (uf q r + 1 + Real.sqrt (1 + q ^ 2))) * (1 + r / Real.sqrt (q ^ 2 + r ^ 2))
      = Jintegrand q r := by
  obtain ⟨b1, b2, hβ⟩ := beta_facts hq hq1
  have hu := uf_ge hq.le hr.le
  obtain ⟨p1, p2, p3, p4⟩ := args_pos hq hq1 hu
  set β := Real.sqrt (1 + q ^ 2) with hβdef
  set R := Real.sqrt (q ^ 2 + r ^ 2) with hRdef
  have hR : R ^ 2 = q ^ 2 + r ^ 2 := Real.sq_sqrt (by positivity)
  have hR0 : 0 < R := Real.sqrt_pos.mpr (by positivity)
  set u := uf q r with hudef
  have hur : u = r + R := rfl
  have hu0 : 0 < u := by linarith
  have hq1' : q ^ 2 - 1 ≠ 0 := by nlinarith
  have hb0 : β ≠ 0 := by linarith
  -- (a)
  have ha : (u + 1 - β) * (u + 1 + β) = 2 * u * (1 + r) := by
    rw [hur]; linear_combination hR - hβ
  -- (b)
  have hb : cC q * (u + 1 + β) + cD q * (u + 1 - β) = 2 * (u + 1) / (1 + q) - 2 := by
    unfold cC cD; rw [← hβdef]; field_simp; ring
  -- (c)
  have hc : cC q / (u + 1 - β) + cD q / (u + 1 + β) = (u - q) / ((1 + q) * u * (1 + r)) := by
    rw [div_add_div _ _ p3.ne' p4.ne', show cC q * (u + 1 + β) + (u + 1 - β) * cD q
      = cC q * (u + 1 + β) + cD q * (u + 1 - β) by ring, hb, ha]
    field_simp; ring
  -- (d)
  have hd : cA q / (u + 1) + cB q / (u + q) = -(2 * u) / ((1 + q) * (u + 1) * (u + q)) := by
    unfold cA cB
    have : (1 + q) ≠ 0 := by linarith
    field_simp; ring
  have hrel : (u - q) * (u + q) = 2 * r * u := by rw [hur]; linear_combination hR
  have huq : u - q = 2 * r * u / (u + q) := by rw [eq_div_iff p2.ne']; exact hrel
  have huR : 1 + r / R = u / R := by rw [hur]; field_simp; ring
  rw [show cA q / (u + 1) + cB q / (u + q) + cC q / (u + 1 - β) + cD q / (u + 1 + β)
      = (cA q / (u + 1) + cB q / (u + q)) + (cC q / (u + 1 - β) + cD q / (u + 1 + β)) by ring,
    hc, hd, huR, huq]
  unfold Jintegrand
  rw [← hRdef, show 1 + r + R = u + 1 by rw [hur]; ring]
  have hqu : 1 / (u + 1) - 1 / (1 + q) = -(u - q) / ((u + 1) * (1 + q)) := by
    field_simp; ring
  rw [hqu, huq]
  have : R = u - r := by rw [hur]; ring
  rw [this]
  have hur0 : u - r ≠ 0 := by rw [← this]; exact hR0.ne'
  field_simp
  ring

lemma csum : cA q + cB q + cC q + cD q = 0 := by
  obtain ⟨b1, -, -⟩ := beta_facts hq hq1
  have hq1' : q ^ 2 - 1 ≠ 0 := by nlinarith
  have : (1 + q) ≠ 0 := by linarith
  have : Real.sqrt (1 + q ^ 2) ≠ 0 := by linarith
  unfold cA cB cC cD
  field_simp; ring

lemma tendsto_Fu : Tendsto (fun r => Fp q (uf q r)) atTop (𝓝 0) := by
  set β := Real.sqrt (1 + q ^ 2) with hβdef
  have hH : Tendsto (fun u : ℝ => cA q * Real.log (1 + 1 / u) + cB q * Real.log (1 + q / u)
      + cC q * Real.log (1 + (1 - β) / u) + cD q * Real.log (1 + (1 + β) / u)) atTop (𝓝 0) := by
    have := ((((tendsto_log_one_add 1).const_mul (cA q)).add
      ((tendsto_log_one_add q).const_mul (cB q))).add
      ((tendsto_log_one_add (1 - β)).const_mul (cC q))).add
      ((tendsto_log_one_add (1 + β)).const_mul (cD q))
    simpa using this
  have hu : Tendsto (uf q) atTop atTop := tendsto_atTop_mono (fun r => uf_ge_r) tendsto_id
  refine (hH.comp hu).congr' ?_
  filter_upwards [eventually_ge_atTop (0 : ℝ)] with r hr
  have hge := uf_ge hq.le hr
  obtain ⟨p1, p2, p3, p4⟩ := args_pos hq hq1 hge
  have hu0 : 0 < uf q r := by linarith
  have key : ∀ d : ℝ, 0 < uf q r + d →
      Real.log (1 + d / uf q r) = Real.log (uf q r + d) - Real.log (uf q r) := by
    intro d hd
    rw [← Real.log_div hd.ne' hu0.ne']
    congr 1; field_simp
  simp only [Function.comp]
  rw [key 1 p1, key q p2, key (1 - β) (by linarith), key (1 + β) (by linarith)]
  unfold Fp
  rw [← hβdef]
  have hc := csum hq hq1
  have e3 : uf q r + (1 - β) = uf q r + 1 - β := by ring
  have e4 : uf q r + (1 + β) = uf q r + 1 + β := by ring
  rw [e3, e4]
  linear_combination (-Real.log (uf q r)) * hc

omit hq1 in
lemma uf_zero : uf q 0 = q := by
  unfold uf; simp [Real.sqrt_sq hq.le]

lemma J_eq_neg_F : J q = -Fp q q := by
  unfold J
  have hcont : ContinuousWithinAt (fun r => Fp q (uf q r)) (Ici 0) 0 :=
    (hasDerivAt_Fu hq hq1 le_rfl).continuousAt.continuousWithinAt
  have hderiv : ∀ r ∈ Ioi (0 : ℝ), HasDerivAt (fun r => Fp q (uf q r)) (Jintegrand q r) r := by
    intro r hr
    exact (hasDerivAt_Fu hq hq1 (le_of_lt hr)).congr_deriv (Dv_eq hq hq1 hr)
  have hneg : ∀ r ∈ Ioi (0 : ℝ), Jintegrand q r ≤ 0 := by
    intro r hr
    have hr0 : (0:ℝ) < r := hr
    unfold Jintegrand
    have hge := uf_ge hq.le hr0.le
    unfold uf at hge
    apply div_nonpos_of_nonpos_of_nonneg _ (by positivity)
    rw [sub_nonpos]
    apply one_div_le_one_div_of_le (by linarith)
    linarith
  rw [integral_Ioi_of_hasDerivAt_of_nonpos hcont hderiv hneg (tendsto_Fu hq hq1), uf_zero hq]
  ring

lemma arsinh_sum : Real.arsinh q + Real.arsinh (1 / q) =
    Real.log (q + 1 + Real.sqrt (1 + q ^ 2)) - Real.log (q + 1 - Real.sqrt (1 + q ^ 2)) := by
  obtain ⟨b1, b2, hβ⟩ := beta_facts hq hq1
  set β := Real.sqrt (1 + q ^ 2) with hβdef
  have h1 : Real.sqrt (1 + (1 / q) ^ 2) = β / q := by
    rw [show 1 + (1 / q) ^ 2 = (1 + q ^ 2) / q ^ 2 by field_simp; ring, Real.sqrt_div' _ (by positivity),
      Real.sqrt_sq hq.le]
  unfold Real.arsinh
  rw [h1, ← hβdef, ← Real.log_mul (by positivity) (by positivity),
    ← Real.log_div (by linarith) (by linarith)]
  congr 1
  rw [eq_div_iff (by linarith)]
  field_simp
  linear_combination (-β) * hβ

lemma prod_two_q : Real.log (q + 1 + Real.sqrt (1 + q ^ 2)) + Real.log (q + 1 - Real.sqrt (1 + q ^ 2))
    = Real.log (2 * q) := by
  obtain ⟨b1, b2, hβ⟩ := beta_facts hq hq1
  rw [← Real.log_mul (by linarith) (by linarith)]
  congr 1
  linear_combination -hβ

lemma J_closed : J q = f1 q - (Real.arsinh q + Real.arsinh (1 / q)) / Real.sqrt (1 + q ^ 2) := by
  rw [J_eq_neg_F hq hq1, arsinh_sum hq hq1]
  have hp := prod_two_q hq hq1
  obtain ⟨b1, b2, hβ⟩ := beta_facts hq hq1
  set β := Real.sqrt (1 + q ^ 2) with hβdef
  set X := Real.log (q + 1 + β)
  set Y := Real.log (q + 1 - β)
  have hY : Y = Real.log (2 * q) - X := by linarith
  unfold Fp f1 cA cB cC cD
  rw [← hβdef, show q + 1 - β = q + 1 - β from rfl, show Real.log (q + q) = Real.log (2 * q) by ring_nf,
    show Real.log (q + 1) = Real.log (1 + q) by ring_nf]
  change _ = _ - (X - Y) / β
  rw [show Real.log (q + 1 - β) = Y from rfl, show Real.log (q + 1 + β) = X from rfl, hY]
  have hq1' : q ^ 2 - 1 ≠ 0 := by nlinarith
  have hq1'' : 1 - q ^ 2 ≠ 0 := by nlinarith
  have : (1 + q) ≠ 0 := by linarith
  have : β ≠ 0 := by linarith
  field_simp
  ring

end


/-- The two real log series at `z = ρ e^{iθ}`, `0 ≤ ρ < 1`. -/
lemma log_series_parts {ρ θ : ℝ} (h0 : 0 ≤ ρ) (h1 : ρ < 1) :
    HasSum (fun n : ℕ => ρ ^ (n + 1) * Real.cos ((n + 1 : ℝ) * θ) / ((n : ℝ) + 1))
        (-(Complex.log (1 - (ρ : ℂ) * (Complex.cos θ + Complex.sin θ * Complex.I)))).re ∧
      HasSum (fun n : ℕ => ρ ^ (n + 1) * Real.sin ((n + 1 : ℝ) * θ) / ((n : ℝ) + 1))
        (-(Complex.log (1 - (ρ : ℂ) * (Complex.cos θ + Complex.sin θ * Complex.I)))).im := by
  set z : ℂ := (ρ : ℂ) * (Complex.cos θ + Complex.sin θ * Complex.I) with hz
  have hzn : ‖z‖ < 1 := by
    rw [hz, norm_mul, Complex.norm_cos_add_sin_mul_I, mul_one, Complex.norm_real,
      Real.norm_of_nonneg h0]
    exact h1
  have hS := Complex.hasSum_taylorSeries_neg_log' hzn
  have hterm : ∀ n : ℕ, z ^ (n + 1) / ((n : ℂ) + 1) =
      ((ρ ^ (n + 1) * Real.cos ((n + 1 : ℝ) * θ) / ((n : ℝ) + 1) : ℝ) : ℂ)
        + ((ρ ^ (n + 1) * Real.sin ((n + 1 : ℝ) * θ) / ((n : ℝ) + 1) : ℝ) : ℂ) * Complex.I := by
    intro n
    rw [hz, mul_pow, Complex.cos_add_sin_mul_I_pow]
    push_cast
    have : ((n : ℂ) + 1) ≠ 0 := by exact_mod_cast (Nat.succ_ne_zero n)
    field_simp
  simp_rw [hterm] at hS
  constructor
  · convert Complex.reCLM.hasSum hS using 1
    swap; · rfl
    funext n
    simp only [Complex.reCLM_apply, Complex.add_re, Complex.ofReal_re, Complex.mul_re,
      Complex.ofReal_im, Complex.I_re, Complex.I_im]
    ring
  · convert Complex.imCLM.hasSum hS using 1
    swap; · rfl
    funext n
    simp only [Complex.imCLM_apply, Complex.add_im, Complex.ofReal_re, Complex.mul_im,
      Complex.ofReal_im, Complex.I_re, Complex.I_im]
    ring

lemma trig_aux (x : ℝ) :
    1 - (2 * Real.cos (2 * x) - 1) * Real.cos (2 * x) = 2 * Real.sin x * Real.sin (3 * x) ∧
    (2 * Real.cos (2 * x) - 1) * Real.sin (2 * x) = 2 * Real.sin x * Real.cos (3 * x) := by
  have h := Real.sin_sq_add_cos_sq x
  rw [Real.cos_two_mul, Real.sin_two_mul, Real.sin_three_mul, Real.cos_three_mul]
  constructor
  · linear_combination (8 * Real.sin x ^ 2 + 2 - 8 * Real.cos x ^ 2) * h
  · ring

/-- `1 - ζ_θ = 2 sin(θ/2) e^{i(3θ/2 - π/2)}`. -/
lemma one_sub_zeta (θ : ℝ) :
    1 - ((2 * Real.cos θ - 1 : ℝ) : ℂ) * (Complex.cos θ + Complex.sin θ * Complex.I)
      = ((2 * Real.sin (θ / 2) : ℝ) : ℂ) * (Complex.cos ((3 * θ / 2 - π / 2 : ℝ) : ℂ)
        + Complex.sin ((3 * θ / 2 - π / 2 : ℝ) : ℂ) * Complex.I) := by
  obtain ⟨h1, h2⟩ := trig_aux (θ / 2)
  rw [show 2 * (θ / 2) = θ by ring, show 3 * (θ / 2) = 3 * θ / 2 by ring] at h1 h2
  have c1 : Real.cos (3 * θ / 2 - π / 2) = Real.sin (3 * θ / 2) := by
    rw [Real.cos_sub_pi_div_two]
  have s1 : Real.sin (3 * θ / 2 - π / 2) = -Real.cos (3 * θ / 2) := by
    rw [Real.sin_sub_pi_div_two]
  apply Complex.ext
  · simp only [Complex.sub_re, Complex.one_re, Complex.mul_re, Complex.ofReal_re, Complex.ofReal_im,
      Complex.add_re, Complex.add_im, Complex.I_re, Complex.I_im, ← Complex.ofReal_cos,
      ← Complex.ofReal_sin]
    rw [c1, s1]; linear_combination h1
  · simp only [Complex.sub_im, Complex.one_im, Complex.mul_im, Complex.ofReal_re, Complex.ofReal_im,
      Complex.add_re, Complex.add_im, Complex.I_re, Complex.I_im, ← Complex.ofReal_cos,
      ← Complex.ofReal_sin]
    rw [c1, s1]; linear_combination -h2

lemma log_zeta {θ : ℝ} (h0 : 0 < θ) (h1 : θ < π / 3) :
    (-(Complex.log (1 - ((2 * Real.cos θ - 1 : ℝ) : ℂ) * (Complex.cos θ + Complex.sin θ * Complex.I)))).re
        = -(1 / 2) * Real.log (1 - (2 * Real.cos θ - 1)) ∧
    (-(Complex.log (1 - ((2 * Real.cos θ - 1 : ℝ) : ℂ) * (Complex.cos θ + Complex.sin θ * Complex.I)))).im
        = π / 2 - 3 * θ / 2 := by
  have hpi := Real.pi_pos
  have hs : 0 < Real.sin (θ / 2) := Real.sin_pos_of_pos_of_lt_pi (by linarith) (by linarith)
  have hr : (0:ℝ) < 2 * Real.sin (θ / 2) := by linarith
  rw [one_sub_zeta]
  constructor
  · rw [Complex.neg_re, Complex.log_re, norm_mul, Complex.norm_real, Real.norm_of_nonneg hr.le,
      Complex.norm_cos_add_sin_mul_I, mul_one]
    have : 1 - (2 * Real.cos θ - 1) = (2 * Real.sin (θ / 2)) ^ 2 := by
      have h3 := Real.cos_two_mul (θ / 2)
      rw [show 2 * (θ / 2) = θ by ring] at h3
      have h2 := Real.sin_sq_add_cos_sq (θ / 2)
      linear_combination (-2) * h3 - 4 * h2
    rw [this, Real.log_pow]; push_cast; ring
  · rw [Complex.neg_im, Complex.log_im, Complex.arg_mul_cos_add_sin_mul_I hr]
    · ring
    · constructor <;> linarith

/-- `r(θ) = 2cos θ - 1`. -/
noncomputable def rr (θ : ℝ) : ℝ := 2 * Real.cos θ - 1

/-- The terms of `Re Li₂(ζ_θ)`. -/
noncomputable def aT (n : ℕ) (θ : ℝ) : ℝ :=
  rr θ ^ (n + 1) * Real.cos ((n + 1 : ℝ) * θ) / ((n : ℝ) + 1) ^ 2

/-- Their `θ`-derivatives. -/
noncomputable def aD (n : ℕ) (θ : ℝ) : ℝ :=
  (rr θ ^ n * (-2 * Real.sin θ) * Real.cos ((n + 1 : ℝ) * θ)
    - rr θ ^ (n + 1) * Real.sin ((n + 1 : ℝ) * θ)) / ((n : ℝ) + 1)

lemma hasDerivAt_rr (θ : ℝ) : HasDerivAt rr (-2 * Real.sin θ) θ := by
  have := ((Real.hasDerivAt_cos θ).const_mul 2).sub_const 1
  show HasDerivAt (fun x => 2 * Real.cos x - 1) _ θ
  convert this using 1; ring

lemma hasDerivAt_aT (n : ℕ) (θ : ℝ) : HasDerivAt (aT n) (aD n θ) θ := by
  have hp := (hasDerivAt_rr θ).pow (n + 1)
  have hc := ((hasDerivAt_id θ).const_mul ((n : ℝ) + 1)).cos
  have h := (hp.mul hc).div_const (((n : ℝ) + 1) ^ 2)
  convert h using 1
  · funext x; simp only [aT, Pi.pow_apply, Pi.mul_apply, id]
  unfold aD
  simp only [id, Nat.add_sub_cancel, mul_one, Pi.pow_apply]
  have : ((n : ℝ) + 1) ≠ 0 := by positivity
  push_cast
  field_simp
  ring

lemma abs_rr_le {θ : ℝ} (h0 : 0 ≤ θ) (h1 : θ ≤ π / 3) : |rr θ| ≤ 1 := by
  have hpi := Real.pi_pos
  have c1 : Real.cos θ ≤ 1 := Real.cos_le_one θ
  have c2 : Real.cos (π / 3) ≤ Real.cos θ :=
    Real.cos_le_cos_of_nonneg_of_le_pi h0 (by linarith) h1
  rw [Real.cos_pi_div_three] at c2
  unfold rr; rw [abs_le]; constructor <;> linarith

lemma abs_aT_le {θ : ℝ} (h0 : 0 ≤ θ) (h1 : θ ≤ π / 3) (n : ℕ) :
    ‖aT n θ‖ ≤ 1 / ((n : ℝ) + 1) ^ 2 := by
  unfold aT
  rw [Real.norm_eq_abs, abs_div, abs_mul, abs_pow, abs_of_pos (by positivity : (0:ℝ) < ((n:ℝ)+1)^2)]
  apply div_le_div_of_nonneg_right _ (by positivity)
  have := abs_rr_le h0 h1
  calc |rr θ| ^ (n + 1) * |Real.cos ((n + 1 : ℝ) * θ)| ≤ 1 * 1 := by
        apply mul_le_mul (pow_le_one₀ (abs_nonneg _) this) (Real.abs_cos_le_one _) (abs_nonneg _)
          (by norm_num)
    _ = 1 := by ring

lemma rr_pos_lt {θ : ℝ} (h0 : 0 < θ) (h1 : θ < π / 3) : 0 < rr θ ∧ rr θ < 1 := by
  have hpi := Real.pi_pos
  have c1 : Real.cos θ < Real.cos 0 :=
    Real.cos_lt_cos_of_nonneg_of_le_pi le_rfl (by linarith) h0
  have c2 : Real.cos (π / 3) < Real.cos θ :=
    Real.cos_lt_cos_of_nonneg_of_le_pi h0.le (by linarith) h1
  rw [Real.cos_pi_div_three] at c2
  rw [Real.cos_zero] at c1
  unfold rr; constructor <;> linarith

lemma hasSum_aD {θ : ℝ} (h0 : 0 < θ) (h1 : θ < π / 3) :
    HasSum (fun n => aD n θ)
      ((-2 * Real.sin θ) / rr θ * (-(1 / 2) * Real.log (1 - rr θ)) - (π / 2 - 3 * θ / 2)) := by
  obtain ⟨r0, r1⟩ := rr_pos_lt h0 h1
  obtain ⟨hc, hs⟩ := log_series_parts (θ := θ) r0.le r1
  obtain ⟨e1, e2⟩ := log_zeta h0 h1
  rw [show (2 * Real.cos θ - 1) = rr θ from rfl] at e1 e2
  rw [e1] at hc; rw [e2] at hs
  have := (hc.mul_left ((-2 * Real.sin θ) / rr θ)).sub hs
  convert this using 1
  funext n
  unfold aD
  have : ((n : ℝ) + 1) ≠ 0 := by positivity
  field_simp
  ring

lemma hasDerivAt_tsum_aT {θ : ℝ} (h0 : 0 < θ) (h1 : θ < π / 3) :
    HasDerivAt (fun z => ∑' n, aT n z) (∑' n, aD n θ) θ := by
  have hpi := Real.pi_pos
  set t := Ioo (θ / 2) ((θ + π / 2) / 2) with ht
  set ρ0 := max (2 * Real.cos (θ / 2) - 1) (1 - 2 * Real.cos ((θ + π / 2) / 2)) with hρ0
  have hc1 : Real.cos (θ / 2) < 1 := by
    have := Real.cos_lt_cos_of_nonneg_of_le_pi le_rfl (by linarith) (by linarith : 0 < θ / 2)
    rwa [Real.cos_zero] at this
  have hc2 : 0 < Real.cos ((θ + π / 2) / 2) :=
    Real.cos_pos_of_mem_Ioo ⟨by linarith, by linarith⟩
  have hρ1 : ρ0 < 1 := max_lt (by linarith) (by linarith)
  have hbd : ∀ y ∈ t, |rr y| ≤ ρ0 := by
    intro y hy
    obtain ⟨y1, y2⟩ := hy
    have a1 : Real.cos y < Real.cos (θ / 2) :=
      Real.cos_lt_cos_of_nonneg_of_le_pi (by linarith) (by linarith) y1
    have a2 : Real.cos ((θ + π / 2) / 2) < Real.cos y :=
      Real.cos_lt_cos_of_nonneg_of_le_pi (by linarith) (by linarith) y2
    unfold rr; rw [abs_le]
    constructor
    · have := le_max_right (2 * Real.cos (θ / 2) - 1) (1 - 2 * Real.cos ((θ + π / 2) / 2))
      linarith
    · have := le_max_left (2 * Real.cos (θ / 2) - 1) (1 - 2 * Real.cos ((θ + π / 2) / 2))
      linarith
  have hρ0 : 0 ≤ ρ0 := le_trans (abs_nonneg _) (hbd θ ⟨by linarith, by linarith⟩)
  have hu : Summable (fun n : ℕ => 3 * ρ0 ^ n) :=
    (summable_geometric_of_lt_one hρ0 hρ1).mul_left 3
  have hθt : θ ∈ t := ⟨by linarith, by linarith⟩
  refine hasDerivAt_tsum_of_isPreconnected (t := t) hu isOpen_Ioo isPreconnected_Ioo
    (fun n y _ => hasDerivAt_aT n y) ?_ hθt ?_ hθt
  · intro n y hy
    have hr := hbd y hy
    have hn : (1 : ℝ) ≤ (n : ℝ) + 1 := by linarith [(Nat.cast_nonneg n : (0:ℝ) ≤ n)]
    unfold aD
    rw [Real.norm_eq_abs, abs_div, abs_of_pos (by positivity : (0:ℝ) < (n:ℝ) + 1),
      div_le_iff₀ (by positivity)]
    have p1 : |rr y ^ n * (-2 * Real.sin y) * Real.cos ((n + 1 : ℝ) * y)| ≤ ρ0 ^ n * 2 := by
      rw [abs_mul, abs_mul, abs_pow, abs_mul, abs_neg, abs_two]
      have := pow_le_pow_left₀ (abs_nonneg _) hr n
      have s1 := Real.abs_sin_le_one y
      have s2 := Real.abs_cos_le_one ((n + 1 : ℝ) * y)
      calc |rr y| ^ n * (2 * |Real.sin y|) * |Real.cos ((n + 1 : ℝ) * y)|
          ≤ ρ0 ^ n * (2 * 1) * 1 := by gcongr
        _ = ρ0 ^ n * 2 := by ring
    have p2 : |rr y ^ (n + 1) * Real.sin ((n + 1 : ℝ) * y)| ≤ ρ0 ^ n := by
      rw [abs_mul, abs_pow]
      have := pow_le_pow_left₀ (abs_nonneg _) hr (n + 1)
      have s2 := Real.abs_sin_le_one ((n + 1 : ℝ) * y)
      have : ρ0 ^ (n + 1) ≤ ρ0 ^ n := pow_le_pow_of_le_one hρ0 hρ1.le (by omega)
      calc |rr y| ^ (n + 1) * |Real.sin ((n + 1 : ℝ) * y)| ≤ ρ0 ^ (n + 1) * 1 := by gcongr
        _ ≤ ρ0 ^ n := by linarith
    have := abs_sub _ _ |>.trans (add_le_add p1 p2)
    have h3 : 0 ≤ ρ0 ^ n := pow_nonneg hρ0 n
    nlinarith
  · exact Summable.of_norm_bounded D00.summable_inv_sq (abs_aT_le h0.le h1.le)

/-- The right side of Lemma A. -/
noncomputable def RA (θ : ℝ) : ℝ :=
  (1 / 2) * D00.Li2 (rr θ) + (3 / 4) * θ ^ 2 - π / 2 * θ + π ^ 2 / 12

lemma tsum_aT_eq {θ : ℝ} (h0 : 0 ≤ θ) (h1 : θ ≤ π / 3) : ∑' n, aT n θ = RA θ := by
  have hpi := Real.pi_pos
  have hrr : Continuous rr := (continuous_const.mul Real.continuous_cos).sub continuous_const
  have hcont : ContinuousOn (fun z => ∑' n, aT n z - RA z) (Icc θ (π / 3)) := by
    have c1 : ContinuousOn (fun z => ∑' n, aT n z) (Icc θ (π / 3)) := by
      refine continuousOn_tsum (u := fun n : ℕ => 1 / ((n : ℝ) + 1) ^ 2) (fun n => ?_)
        D00.summable_inv_sq (fun n z hz => abs_aT_le (by linarith [hz.1]) hz.2 n)
      exact (continuous_iff_continuousAt.2 (fun z => (hasDerivAt_aT n z).continuousAt)).continuousOn
    have c2 : Continuous RA := by
      unfold RA
      exact (((continuous_const.mul (D00.continuous_Li2.comp hrr)).add
        (continuous_const.mul (continuous_pow 2))).sub (continuous_const.mul continuous_id)).add
        continuous_const
    exact c1.sub c2.continuousOn
  have hd : ∀ x ∈ Ioo θ (π / 3), HasDerivAt (fun z => ∑' n, aT n z - RA z) 0 x := by
    intro x hx
    have hx0 : 0 < x := by linarith [hx.1]
    have hx1 : x < π / 3 := hx.2
    obtain ⟨r0, r1⟩ := rr_pos_lt hx0 hx1
    have e1 := hasDerivAt_tsum_aT hx0 hx1
    rw [(hasSum_aD hx0 hx1).tsum_eq] at e1
    have e2 : HasDerivAt RA ((1 / 2) * (-Real.log (1 - rr x) / rr x * (-2 * Real.sin x))
        + (3 / 4) * (2 * x) - π / 2) x := by
      have hl := (D00.hasDerivAt_Li2 r1 r0.ne').comp x (hasDerivAt_rr x)
      have hq : HasDerivAt (fun z : ℝ => z ^ 2) (2 * x) x := by
        simpa using hasDerivAt_pow 2 x
      have := (((hl.const_mul (1 / 2)).add (hq.const_mul (3 / 4))).sub
        ((hasDerivAt_id x).const_mul (π / 2))).add_const (π ^ 2 / 12)
      convert this using 1
      · funext z; simp only [RA, Function.comp, id, Pi.add_apply, Pi.sub_apply]
      · ring
    have := e1.sub e2
    convert this using 1
    field_simp
    ring
  have key := D00.eq_of_deriv_zero h1 hcont hd
  have hend : ∑' n, aT n (π / 3) - RA (π / 3) = 0 := by
    have hz : rr (π / 3) = 0 := by unfold rr; rw [Real.cos_pi_div_three]; norm_num
    simp only [RA, aT, hz, D00.Li2_zero]
    simp
    ring
  linarith [key, hend]

/-- Lemma A, the general case. -/
lemma lemmaA {θ : ℝ} (h0 : 0 ≤ θ) (h1 : θ ≤ π / 3) :
    HasSum (fun n : ℕ => (2 * Real.cos θ - 1) ^ (n + 1) * Real.cos ((n + 1 : ℝ) * θ)
        / ((n : ℝ) + 1) ^ 2)
      ((1 / 2) * D00.Li2 (2 * Real.cos θ - 1) + (3 / 4) * θ ^ 2 - π / 2 * θ + π ^ 2 / 12) := by
  have hs : Summable (fun n => aT n θ) :=
    Summable.of_norm_bounded D00.summable_inv_sq (abs_aT_le h0 h1)
  have := hs.hasSum
  rw [tsum_aT_eq h0 h1] at this
  exact this

lemma lemmaA_pi4 :
    HasSum (fun n : ℕ => D02.b ^ (n + 1) * Real.cos ((n + 1 : ℝ) * (π / 4)) / ((n : ℝ) + 1) ^ 2)
      ((1 / 2) * D00.Li2 D02.b + π ^ 2 / 192) := by
  have hpi := Real.pi_pos
  have h := lemmaA (θ := π / 4) (by positivity) (by linarith)
  have hb : 2 * Real.cos (π / 4) - 1 = D02.b := by
    rw [Real.cos_pi_div_four]; unfold D02.b; ring
  rw [hb] at h
  convert h using 1
  ring



lemma continuous_qs : Continuous qs := by unfold qs; fun_prop

lemma qs_pos (s : ℝ) : 0 < qs s := by
  unfold qs; apply Real.sqrt_pos.mpr
  rcases le_or_gt s (1/2) with h | h <;> nlinarith [sq_nonneg s, sq_nonneg (1 - s)]

lemma qs_bounds {s : ℝ} (h0 : 0 ≤ s) (h1 : s ≤ 1) : 1 / 2 ≤ qs s ∧ qs s ≤ 1 := by
  unfold qs
  constructor
  · apply Real.le_sqrt_of_sq_le; nlinarith [sq_nonneg (s - 1/2)]
  · rw [Real.sqrt_le_one]; nlinarith

lemma f1_abs_le {q : ℝ} (h0 : 1 / 2 ≤ q) (h1 : q ≤ 1) : |f1 q| ≤ 7 := by
  rcases eq_or_lt_of_le h1 with rfl | h1
  · simp [f1]
  have hq0 : 0 < q := by linarith
  set A := Real.log ((1 + q) / 2) with hA
  set B := Real.log q with hB
  have a1 : A ≤ (1 + q) / 2 - 1 := Real.log_le_sub_one_of_pos (by linarith)
  have a2 : 1 - ((1 + q) / 2)⁻¹ ≤ A := Real.one_sub_inv_le_log_of_pos (by linarith)
  have b1 : B ≤ q - 1 := Real.log_le_sub_one_of_pos hq0
  have b2 : 1 - q⁻¹ ≤ B := Real.one_sub_inv_le_log_of_pos hq0
  have a2' : q - 1 ≤ A := by
    have : 1 - ((1 + q) / 2)⁻¹ = (q - 1) / (1 + q) := by field_simp; ring
    rw [this] at a2
    have : q - 1 ≤ (q - 1) / (1 + q) := by
      rw [le_div_iff₀ (by linarith)]; nlinarith
    linarith
  have b2' : 2 * (q - 1) ≤ B := by
    have : 1 - q⁻¹ = (q - 1) / q := by field_simp
    rw [this] at b2
    have : 2 * (q - 1) ≤ (q - 1) / q := by
      rw [le_div_iff₀ hq0]; nlinarith
    linarith
  have hl2 : 0 < Real.log 2 ∧ Real.log 2 < 1 := ⟨Real.log_pos (by norm_num),
    by have := Real.log_two_lt_d9; linarith⟩
  have hN : 2 * Real.log (1 + q) - (1 + q) * Real.log (2 * q)
      = 2 * A + (1 - q) * Real.log 2 - (1 + q) * B := by
    rw [hA, hB, Real.log_div (by linarith) (by norm_num), Real.log_mul (by norm_num) hq0.ne']
    ring
  unfold f1
  rw [hN, abs_div, abs_of_pos (by nlinarith : (0:ℝ) < 1 - q ^ 2), div_le_iff₀ (by nlinarith)]
  rw [abs_le]
  constructor <;> nlinarith [mul_le_mul_of_nonneg_left b1 (by linarith : (0:ℝ) ≤ 1 + q),
    mul_le_mul_of_nonneg_left b2' (by linarith : (0:ℝ) ≤ 1 + q),
    mul_le_mul_of_nonneg_left hl2.1.le (by linarith : (0:ℝ) ≤ 1 - q),
    mul_le_mul_of_nonneg_left hl2.2.le (by linarith : (0:ℝ) ≤ 1 - q)]

lemma ii_f1 : IntervalIntegrable (fun s => f1 (qs s)) volume 0 1 := by
  rw [intervalIntegrable_iff_integrableOn_Ioc_of_le zero_le_one]
  refine Measure.integrableOn_of_bounded (M := 7) (by simp) ?_ ?_
  · have : Measurable (fun s => f1 (qs s)) := by unfold f1 qs; fun_prop
    exact this.aestronglyMeasurable
  · refine (ae_restrict_iff' measurableSet_Ioc).2 (Eventually.of_forall (fun s hs => ?_))
    obtain ⟨h0, h1⟩ := qs_bounds hs.1.le hs.2
    rw [Real.norm_eq_abs]; exact f1_abs_le h0 h1

lemma ii_2a : IntervalIntegrable (fun s => Real.arsinh (qs s) / Real.sqrt (1 + qs s ^ 2)) volume 0 1 := by
  apply Continuous.intervalIntegrable
  have := continuous_qs
  refine (Real.continuous_arsinh.comp this).div (by fun_prop) (fun s => ?_)
  exact (Real.sqrt_pos.mpr (by positivity)).ne'

lemma ii_2c : IntervalIntegrable (fun s => Real.arsinh (1 / qs s) / Real.sqrt (1 + qs s ^ 2)) volume 0 1 := by
  apply Continuous.intervalIntegrable
  have := continuous_qs
  refine (Real.continuous_arsinh.comp (continuous_const.div this (fun s => (qs_pos s).ne'))).div
    (by fun_prop) (fun s => ?_)
  exact (Real.sqrt_pos.mpr (by positivity)).ne'


lemma integrableOn_inv_sq : IntegrableOn (fun r : ℝ => 2 / (1 + r) ^ 2) (Ioi 0) := by
  refine integrableOn_Ioi_deriv_of_nonneg' (g := fun r => -2 / (1 + r)) (l := 0) ?_ ?_ ?_
  · intro r hr
    have hr0 : (0:ℝ) ≤ r := hr
    have h1 : (1 + r) ≠ 0 := by linarith
    have := ((hasDerivAt_id r).const_add 1).inv h1
    convert this.const_mul (-2) using 1
    · funext x; simp [div_eq_mul_inv]
    · simp only [id]; field_simp
  · intro r hr; positivity
  · have : Tendsto (fun r : ℝ => 1 + r) atTop atTop := tendsto_atTop_add_const_left _ 1 tendsto_id
    simpa using (tendsto_const_nhds (x := (-2 : ℝ))).div_atTop this

lemma abs_Jintegrand_le {q r : ℝ} (hq : 0 ≤ q) (hr : 0 < r) :
    |Jintegrand q r| ≤ 2 / (1 + r) ^ 2 := by
  unfold Jintegrand
  set R := Real.sqrt (q ^ 2 + r ^ 2) with hR
  have R1 : q ≤ R := Real.le_sqrt_of_sq_le (by nlinarith)
  have R2 : R ≤ q + r := by
    rw [hR, Real.sqrt_le_left (by linarith)]; nlinarith
  have n1 : (1 + r + R) ≠ 0 := by linarith
  have n2 : (1 + q) ≠ 0 := by linarith
  have hb : |1 / (1 + r + R) - 1 / (1 + q)| ≤ 2 * r / (1 + r) := by
    have e : 1 / (1 + r + R) - 1 / (1 + q) = -((r + R - q) / ((1 + r + R) * (1 + q))) := by
      field_simp; ring
    rw [e, abs_neg, abs_of_nonneg (by apply div_nonneg <;> nlinarith)]
    rw [div_le_div_iff₀ (by nlinarith) (by linarith)]
    nlinarith [mul_nonneg hq hr.le, mul_nonneg hr.le hr.le, mul_nonneg hq hq]
  rw [abs_div, abs_of_pos (by positivity : (0:ℝ) < r * (1 + r)), div_le_iff₀ (by positivity)]
  calc |1 / (1 + r + R) - 1 / (1 + q)| ≤ 2 * r / (1 + r) := hb
    _ = 2 / (1 + r) ^ 2 * (r * (1 + r)) := by field_simp

lemma integrableOn_Jintegrand {q : ℝ} (hq : 0 ≤ q) : IntegrableOn (Jintegrand q) (Ioi 0) := by
  refine integrableOn_inv_sq.mono' ?_ ?_
  · have : Measurable (Jintegrand q) := by unfold Jintegrand; fun_prop
    exact this.aestronglyMeasurable
  · refine (ae_restrict_iff' measurableSet_Ioi).2 (Eventually.of_forall (fun r hr => ?_))
    rw [Real.norm_eq_abs]; exact abs_Jintegrand_le hq hr

lemma qs_lt_one {s : ℝ} (h0 : 0 < s) (h1 : s < 1) : qs s < 1 := by
  unfold qs; rw [Real.sqrt_lt' one_pos]; nlinarith [mul_pos h0 (sub_pos.mpr h1)]

lemma J_qs_eq {s : ℝ} (hs : s ∈ Ioo (0 : ℝ) 1) :
    J (qs s) = f1 (qs s) - Real.arsinh (qs s) / Real.sqrt (1 + qs s ^ 2)
      - Real.arsinh (1 / qs s) / Real.sqrt (1 + qs s ^ 2) := by
  rw [J_closed (qs_pos s) (qs_lt_one hs.1 hs.2), add_div]; ring

lemma ii_J : IntervalIntegrable (fun s => J (qs s)) volume 0 1 := by
  have h := (ii_f1.sub ii_2a).sub ii_2c
  rw [intervalIntegrable_iff_integrableOn_Ioo_of_le zero_le_one] at h ⊢
  exact h.congr_fun (fun s hs => (J_qs_eq hs).symm) measurableSet_Ioo

lemma int_J : ∫ s in (0 : ℝ)..1, J (qs s) = P1 - P2a - P2c := by
  have h1 : ∀ᵐ x ∂(volume : Measure ℝ), x ≠ 1 := by
    simp [ae_iff, measure_singleton]
  rw [intervalIntegral.integral_congr_ae (g := fun s => f1 (qs s)
      - Real.arsinh (qs s) / Real.sqrt (1 + qs s ^ 2)
      - Real.arsinh (1 / qs s) / Real.sqrt (1 + qs s ^ 2)) ?_]
  · rw [intervalIntegral.integral_sub (ii_f1.sub ii_2a) ii_2c,
      intervalIntegral.integral_sub ii_f1 ii_2a]
    rfl
  · filter_upwards [h1] with x hx hxI
    rw [Set.uIoc_of_le zero_le_one] at hxI
    exact J_qs_eq ⟨hxI.1, lt_of_le_of_ne hxI.2 hx⟩

/-- `Q₀ = 6∫J(q(s))ds` gives `Q₀ = 6(P₁ - P_{2a} - P_{2c})`. -/
lemma Q0_P_of_J (h : Z04.Q 0 = 6 * ∫ s in (0 : ℝ)..1, J (qs s)) :
    Z04.Q 0 = 6 * (P1 - P2a - P2c) := by
  rw [h, int_J]

/-- End of the proof of `thm:S`: adding the three closed forms of Lemma `g3:lem:P`. -/
lemma Q0_closed_of (h4 : Z04.Q 0 = 6 * (P1 - P2a - P2c))
    (h6 : P1 = 3 * Real.sqrt 2 * D00.Li2 D02.b - D02.lam ^ 2 / 2
      - Real.sqrt 2 / 2 * D02.lam * D02.L + (1 + 3 * Real.sqrt 2 / 2) * D02.L ^ 2
      - (1 / 48 + 3 * Real.sqrt 2 / 16) * π ^ 2)
    (h7 : P2a = Real.sqrt 2 * (3 * D00.Li2 D02.b - (1 / 4) * D00.Li2 (-D02.b ^ 4)
      - D02.lam * D02.L / 2 + D02.L * Real.log 3 + D02.L ^ 2 / 2 - 5 * π ^ 2 / 24))
    (h8 : P2c = Real.sqrt 2 * (D00.Li2 D02.b - (1 / 8) * D00.Li2 (-D02.b ^ 4)
      + D02.L * Real.log 3 / 2 + π ^ 2 / 96 - π * alphac / 4)) :
    Z04.Q 0 = Q0closed := by
  rw [h4, h6, h7, h8]; unfold Q0closed; ring

/-- `s₂ = 48M₁ - 4π² - 8Q₀` (R04) with `M₁ = M1closed` (W03): the closed form of `Q₀` gives (T6). -/
lemma s2_of_Q0 (h : Z04.Q 0 = Q0closed) : R04.s 2 = s2closed := by
  have h9 := R04.proof_9.2.2.1
  have hM : R04.M 1 = D03.M1closed := by
    unfold R04.M; rw [W03.proof_4, Complex.ofReal_re]
  rw [h9, hM, h]
  unfold D03.M1closed Q0closed s2closed
  simp only [D03.lam, D02.lam, D03.L, D02.L, D02.a, D02.b]
  ring

lemma s2closed_bounds : -4.6892155185 < s2closed ∧ s2closed ≤ -4.6892155184 := by
  refine s2_num (Real.sqrt 2) π alphac D02.lam D02.L (Real.log 3) (D00.Li2 D02.b)
    (D00.Li2 (-D02.b ^ 4)) D03.s_bounds ⟨Real.pi_gt_d20.le, Real.pi_lt_d20.le⟩ alpha_bounds
    D03.lam_bounds D03.L_bounds l3_bounds D03.u_bounds m_bounds s2closed ?_
  unfold s2closed; ring

end S02

theorem S02.proof_3 : S02.stmt_3 := fun _ hq hq1 => S02.J_closed hq hq1

theorem S02.proof_5 : S02.stmt_5 := ⟨fun _ h0 h1 => S02.lemmaA h0 h1, S02.lemmaA_pi4⟩
