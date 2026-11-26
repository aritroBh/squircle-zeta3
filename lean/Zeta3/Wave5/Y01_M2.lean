import Mathlib
import Zeta3.Wave3.W02_Holo
import Zeta3.Wave2.D00_Dilog
/-!
# Y01 -- the weight-three coefficient `M_2` of the squircle expansion (Theorem `thm:M2`)

Target id: Y01.
Paper labels: `thm:M2` (main.tex 901-907), `eq:T` (main.tex 903-905), `eq:M2` = (T1) of Table
`tab:closed` (main.tex 917-931); constants of Section `sec:square` (main.tex 870-877).
Detailed proof: paper/ieee/proofs/sec_g1_M2.tex (Lemmas `g1:lem:Mk`, `g1:lem:euler`,
`g1:lem:k2`, `g1:lem:red2`, `g1:lem:T`, the table of functional equations in
`g1:sub:M2:feq`, Lemma `g1:lem:K`, and the proof in `g1:sub:M2:proof`).

Quoted verbatim (main.tex 870-877):

  \begin{gather*}
  a=1+\sqrt2,\qquad b=\sqrt2-1=\tfrac1a=\tan\tfrac\pi8,\\
  L=\log a,\qquad \lambda=\log2 .
  \end{gather*}
  ... We use $\Li_n(z)=\sum_{k\ge1}z^k/k^n$, so that $\zeta(n)=\Li_n(1)$ for every $n\ge2$,
  even or odd. The series converges for $|z|\le1$ only. Outside the disk we use the analytic
  continuation, which is real for real $z<1$.

Quoted verbatim (main.tex 901-907, Theorem thm:M2 with eq:T):

  \begin{theorem}[$M_2$$^\dagger$]\label{thm:M2}
  $M_2$ is given by \eqref{eq:M2} in Table~\ref{tab:closed}. It is the sum of the weight-two term
  \begin{equation}\label{eq:T}
  T=\frac{\pi^2+6\lambda^2}{24}=\frac{\zeta(2)-\Li_2(\frac12)}{2}
  \end{equation}
  and a combination of weight three of $\zeta(3)$, $\Li_3(\pm b)$, $\Li_3(1/\sqrt2)$ and products
  of lower-weight constants. Numerically $M_2=-0.02216191108726939942\ldots$
  \end{theorem}

Quoted verbatim (main.tex 924-930, eq:M2 = (T1) of Table tab:closed):

  M_2&=\frac{\pi^2+6\lambda^2}{24}
  +\Bigl(\frac{189\sqrt2}{64}-\frac{19}{16}\Bigr)\zeta(3)
  +(11+15\sqrt2)\Li_3(b)+(11-15\sqrt2)\Li_3(-b)-27\sqrt2\,\Li_3\bigl(\tfrac{1}{\sqrt2}\bigr)
  +\Bigl(12L-\frac{7\sqrt2}{2}\lambda\Bigr)\Li_2(b)\notag\\
  &\quad+\Bigl(\frac{15\sqrt2}{8}-\frac7{12}\Bigr)\pi^2L-\Bigl(\frac1{24}+\frac{65\sqrt2}{48}\Bigr)\pi^2\lambda
  +\Bigl(\frac16+\frac{9\sqrt2}{16}\Bigr)\lambda^3+\frac{11\sqrt2}{8}\lambda^2L-\Bigl(1+\frac{7\sqrt2}{4}\Bigr)\lambda L^2
  +\Bigl(\frac{25}{6}-\frac{7\sqrt2}{3}\Bigr)L^3\tag{T1}\label{eq:M2}

Quoted verbatim (supplement, Lemma g1:lem:Mk):

  M_k=-\frac12\int_0^1\varphi_k(t)\,dt,\qquad \varphi_k(t)=[\varepsilon^k]\Phi_\varepsilon(t),

and from the proof of Lemma g1:lem:k2:

  So $\varphi_2=\tfrac12Sc(1-2c+c^3)\ell_1^2+\tfrac12Sc'(1-2c'+c'^3)\ell_2^2+Sc^2c'^2\ell_1\ell_2$,
  and $Sc=1/t$, $Sc'=1/(1-t)$, $Sc^2c'^2=t(1-t)/q^3$.

(with $c=(1-t)/q$, $c'=t/q$, $q=\sqrt{t^2+(1-t)^2}$, $\ell_1=\log t$, $\ell_2=\log(1-t)$).

Quoted verbatim (supplement, Lemma g1:lem:k2):

  M_2=\int_b^a\bigl[\ell_1^2R_{20}+\ell_1\ell_2R_{11}\bigr]dv,
  \qquad
  R_{20}=g_{20}'-\frac1{v+a}+\frac{4-3\sqrt2}{8v},\quad R_{11}=-g_{20}'+\frac{\sqrt2}{8v},
  where $g_{20}'(v)=\sqrt2\,v/(v^2+1)^2$.

(with $t(v)=\frac{(v-b)(v+a)}{4v}=\frac{v^2+2v-1}{4v}$, Lemma g1:lem:euler).

Quoted verbatim (supplement, Lemma g1:lem:red2):

  $M_2=W_3+T$, where
  W_3=\int_b^a\Bigl[\ell_1^2\Bigl(\frac{4-3\sqrt2}{8v}-\frac1{v+a}\Bigr)+\frac{\sqrt2}8\frac{\ell_1\ell_2}v\Bigr]dv,\qquad
  T=\int_b^a\ell_1\,\Rwt\,dv,
  \Rwt=\frac1{2v}+\frac1{2(v-a)}-\frac1{v+a}-\frac1{2(v+b)} .
  All integrals converge.

Quoted verbatim (supplement, Lemma g1:lem:T):

  $T=\dfrac{\pi^2+6\lambda^2}{24}=\dfrac{\zeta(2)-\Li_2(\frac12)}2$.

Quoted verbatim (supplement, table of g1:sub:M2:feq, the trilogarithm rows):

  (I3) & $\Li_3(z)-\Li_3(1/z)=-\frac{\pi^2}6\log(-z)-\frac16\log^3(-z)$ & $z<0$\\
  (Tr) & $\Li_3(z)+\Li_3(1-z)+\Li_3\bigl(1-\frac1z\bigr)=\zeta(3)+\frac16\log^3z+\frac{\pi^2}6\log z-\frac12\log^2z\log(1-z)$ & $0<z<1$\\
  (D3) & $\Li_3(w^2)=4\Li_3(w)+4\Li_3(-w)$ & $|w|\le1$\\

and the Spence--Kummer equation (g1:K), proved on $0<x<y<1$ (Lemma g1:lem:K):

  \Li_3\Bigl(\tfrac{x(1-y)^2}{y(1-x)^2}\Bigr)+\Li_3(xy)+\Li_3\bigl(\tfrac xy\bigr)
  -2\Li_3\Bigl(\tfrac{x(1-y)}{y(1-x)}\Bigr)-2\Li_3\Bigl(\tfrac{x(1-y)}{x-1}\Bigr)-2\Li_3\Bigl(\tfrac{1-y}{1-x}\Bigr)\\
  -2\Li_3\Bigl(\tfrac{1-y}{y(x-1)}\Bigr)-2\Li_3(x)-2\Li_3(y)+2\zeta(3)
  =\log^2y\,\log\tfrac{1-y}{1-x}-\tfrac{\pi^2}{3}\log y-\tfrac13\log^3y .

## Conventions

* `M_2` is `W02.Mk 2 = -F''(0)/(2·2!)` (complex, imported from Zeta3.Wave3.W02_Holo, not
  redefined). The real numbers below are cast to `ℂ`; `Mk 2 = ↑x` with `x : ℝ` also says that
  `M_2` is real.
* `Li₂ = D00.Li2` (imported: `Li₂(x) = -∫_0^x log(1-t)/t dt`, the genuine dilogarithm on
  `(-∞,1]`).
* The trilogarithm is defined as `Li3 x = ∫_0^x Li₂(t)/t dt`. This is the analytic continuation
  of the paper's series on `(-∞,1]` (all arguments used here are real and `≤ 1`); stmt_1 to
  stmt_4 record the properties relied on: integrability of the integrand on `[0,x]` for
  `x ≤ 1` (so the definition never uses the junk value of a non-integrable integral), the
  paper's series `Σ_{k≥1} x^k/k³` for `|x| ≤ 1`, the derivative `Li₂(x)/x`, and `Li3 1 = ζ(3)`.
* `ζ(3)` is `z3 := (riemannZeta 3).re`, the real part of Mathlib's complex `riemannZeta 3`.
  stmt_4 also states that `riemannZeta 3 = ↑z3` (so it is real) and that
  `HasSum (fun n => 1/(n+1)^3) z3` (the relation to Mathlib is via `zeta_nat_eq_tsum_of_gt_one`
  / `hasSum_zeta_nat`-style lemmas in the proof). `ζ(2)` is written `π²/6`.
* `a = 1+√2`, `b = √2-1`, `L = log a`, `λ = log 2` are defined here (`a`, `b`, `L`, `lam`).
* `t(v) = (v²+2v-1)/(4v)` is `tv`; on `(b,a)` it lies in `(0,1)`, so `Real.log (tv v)` and
  `Real.log (1 - tv v)` are genuine logarithms there; the endpoints are a null set.
  `1/(v-a)` at `v = a` is Lean's `0`, again on a null set only.
* Every integral whose value is asserted comes with an explicit `IntervalIntegrable` claim.

## Statements

* stmt_1: `t ↦ Li₂(t)/t` is integrable on `[0,x]` for every `x ≤ 1`.
* stmt_2: `Li3 x = Σ_{n≥0} x^(n+1)/(n+1)³` for `|x| ≤ 1` (the paper's definition of `Li_3`).
* stmt_3: `Li3' (x) = Li₂(x)/x` for `x < 1`, `x ≠ 0`.
* stmt_4: `riemannZeta 3 = ↑z3`, `Σ 1/(n+1)³ = z3`, and `Li3 1 = z3` ("ζ(n) = Li_n(1)").
* stmt_5: Lemma g1:lem:Mk at `k = 2` with the explicit `φ₂` (the analogue of eq:M1t):
  `M_2 = -½ ∫_0^1 φ₂(t) dt`.
* stmt_6: Lemma g1:lem:k2: `M_2 = ∫_b^a [ℓ₁² R₂₀ + ℓ₁ℓ₂ R₁₁] dv`.
* stmt_7: Lemma g1:lem:red2: `M_2 = W_3 + T` with the two integrals `W_3`, `T` over `[b,a]`.
* stmt_8: Lemma g1:lem:T and eq:T: `T = (π²+6λ²)/24 = (ζ(2) - Li₂(½))/2`.
* stmt_9: the trilogarithm functional equations (I3), (Tr), (D3) in the supplement's ranges.
* stmt_10: the Spence--Kummer equation (g1:K) on `0 < x < y < 1`.
* stmt_11: the reduction of the weight-three part: `W_3 = (T1) - (π²+6λ²)/24`.
* stmt_12: Theorem thm:M2, eq:M2: `M_2 = (T1)`.
* stmt_13: Theorem thm:M2: `|M_2 - (-0.02216191108726939942)| < 10⁻²⁰` (the printed digits).

Numerical check before freezing (mpmath, 40 digits): (T1) = -0.022161911087269399423754937739...;
the integrals of stmt_5, stmt_6 and `W_3 + T` of stmt_7 agree with (T1) to better than `10⁻³⁷`;
`T` agrees with `(π²+6λ²)/24` to `10⁻³⁹`; (I3), (Tr), (D3), (K) at test points and at the three
Kummer points of the supplement vanish to `10⁻⁴⁰`; the integral definition of `Li3` agrees with
mpmath's `polylog(3, ·)` at `-2.5`. The printed decimal differs from (T1) by `3.75·10⁻²¹`.

## NOT COVERED

* "a combination of weight three of ζ(3), Li₃(±b), Li₃(1/√2) and products of lower-weight
  constants": the notion of weight is a property of the displayed expression (T1), not a
  statement about the number `M_2`; it is not formalized (it is visible in the definition
  `M2closed`).
* Lemma g1:lem:euler (the substitution identities) and the folding step are used inside the
  proofs of stmt_6/stmt_7, not stated separately.
* Lemma g1:lem:prim, g1:lem:elem and g1:lem:raw (the 96 elementary integrals and the raw form
  of `W_3` with 15 trilogarithms and 17 dilogarithms, Table g1:tab:rawW3), and the two
  reduction tables of Lemmas g1:lem:dilog, g1:lem:trilog: their content is a computer-generated
  table not reproduced here; stmt_11 states their combined conclusion.
* The paper's remark that `Li₃(½)` can replace `ζ(3)`, the 290-digit search, eq:chi3 and the
  duplication remark after the theorem are not part of thm:M2 and are not formalized
  (the duplication relation itself is (D3) in stmt_9).
* The dilogarithm functional equations of the table (R), (I), (La), (D), (Ln), (A) are already
  statements of Zeta3.Wave2.D00_Dilog and are not restated.
-/
-- BEGIN STATEMENTS
namespace Y01

open MeasureTheory Real

/-- `a = 1 + √2`. -/
noncomputable def a : ℝ := 1 + Real.sqrt 2

/-- `b = √2 - 1`. -/
noncomputable def b : ℝ := Real.sqrt 2 - 1

/-- `L = log a`. -/
noncomputable def L : ℝ := Real.log a

/-- `λ = log 2`. -/
noncomputable def lam : ℝ := Real.log 2

/-- `ζ(3)` as a real number: the real part of Mathlib's `riemannZeta 3`. -/
noncomputable def z3 : ℝ := (riemannZeta 3).re

/-- The real trilogarithm `Li₃(x) = ∫_0^x Li₂(t)/t dt` (genuine on `(-∞,1]`). -/
noncomputable def Li3 (x : ℝ) : ℝ := ∫ t in (0 : ℝ)..x, D00.Li2 t / t

/-- The right-hand side of eq:M2 = (T1) of Table tab:closed. -/
noncomputable def M2closed : ℝ :=
  (π ^ 2 + 6 * lam ^ 2) / 24
    + (189 * Real.sqrt 2 / 64 - 19 / 16) * z3
    + (11 + 15 * Real.sqrt 2) * Li3 b + (11 - 15 * Real.sqrt 2) * Li3 (-b)
    - 27 * Real.sqrt 2 * Li3 (1 / Real.sqrt 2)
    + (12 * L - 7 * Real.sqrt 2 / 2 * lam) * D00.Li2 b
    + (15 * Real.sqrt 2 / 8 - 7 / 12) * π ^ 2 * L
    - (1 / 24 + 65 * Real.sqrt 2 / 48) * π ^ 2 * lam
    + (1 / 6 + 9 * Real.sqrt 2 / 16) * lam ^ 3
    + 11 * Real.sqrt 2 / 8 * lam ^ 2 * L
    - (1 + 7 * Real.sqrt 2 / 4) * lam * L ^ 2
    + (25 / 6 - 7 * Real.sqrt 2 / 3) * L ^ 3

/-- `c = (1-t)/q` (supplement, Lemma g1:lem:k2), `q = W02.q`. -/
noncomputable def c (t : ℝ) : ℝ := (1 - t) / W02.q t

/-- `c' = t/q`. -/
noncomputable def c' (t : ℝ) : ℝ := t / W02.q t

/-- `φ₂ = [ε²]Φ_ε = (1/(2t))(1-2c+c³)ℓ₁² + (1/(2(1-t)))(1-2c'+c'³)ℓ₂² + t(1-t)/q³ ℓ₁ℓ₂`. -/
noncomputable def phi2 (t : ℝ) : ℝ :=
  1 / (2 * t) * (1 - 2 * c t + c t ^ 3) * Real.log t ^ 2
    + 1 / (2 * (1 - t)) * (1 - 2 * c' t + c' t ^ 3) * Real.log (1 - t) ^ 2
    + t * (1 - t) / W02.q t ^ 3 * Real.log t * Real.log (1 - t)

/-- The Euler substitution `t(v) = (v² + 2v - 1)/(4v)` (Lemma g1:lem:euler). -/
noncomputable def tv (v : ℝ) : ℝ := (v ^ 2 + 2 * v - 1) / (4 * v)

/-- `ℓ₁ = log t(v)`. -/
noncomputable def l1 (v : ℝ) : ℝ := Real.log (tv v)

/-- `ℓ₂ = log (1 - t(v))`. -/
noncomputable def l2 (v : ℝ) : ℝ := Real.log (1 - tv v)

/-- `g₂₀'(v) = √2 v/(v²+1)²`. -/
noncomputable def g20' (v : ℝ) : ℝ := Real.sqrt 2 * v / (v ^ 2 + 1) ^ 2

/-- `R₂₀ = g₂₀' - 1/(v+a) + (4-3√2)/(8v)`. -/
noncomputable def R20 (v : ℝ) : ℝ := g20' v - 1 / (v + a) + (4 - 3 * Real.sqrt 2) / (8 * v)

/-- `R₁₁ = -g₂₀' + √2/(8v)`. -/
noncomputable def R11 (v : ℝ) : ℝ := -g20' v + Real.sqrt 2 / (8 * v)

/-- The integrand of `W_3` (Lemma g1:lem:red2). -/
noncomputable def W3integrand (v : ℝ) : ℝ :=
  l1 v ^ 2 * ((4 - 3 * Real.sqrt 2) / (8 * v) - 1 / (v + a))
    + Real.sqrt 2 / 8 * (l1 v * l2 v / v)

/-- `R_wt = 1/(2v) + 1/(2(v-a)) - 1/(v+a) - 1/(2(v+b))`. -/
noncomputable def Rwt (v : ℝ) : ℝ :=
  1 / (2 * v) + 1 / (2 * (v - a)) - 1 / (v + a) - 1 / (2 * (v + b))

/-- `W_3 = ∫_b^a [ℓ₁²((4-3√2)/(8v) - 1/(v+a)) + (√2/8) ℓ₁ℓ₂/v] dv`. -/
noncomputable def W3 : ℝ := ∫ v in b..a, W3integrand v

/-- `T = ∫_b^a ℓ₁ R_wt dv`. -/
noncomputable def Tint : ℝ := ∫ v in b..a, l1 v * Rwt v

/-- Properties relied on (1): the integrand `Li₂(t)/t` of `Li3` is integrable on `[0,x]`
for every `x ≤ 1`. -/
def stmt_1 : Prop :=
  ∀ x : ℝ, x ≤ 1 → IntervalIntegrable (fun t : ℝ => D00.Li2 t / t) volume 0 x

/-- Properties relied on (2), main.tex 877: "$\Li_n(z)=\sum_{k\ge1}z^k/k^n$" at `n = 3`, for
`|x| ≤ 1`. -/
def stmt_2 : Prop :=
  ∀ x : ℝ, |x| ≤ 1 →
    HasSum (fun n : ℕ => x ^ (n + 1) / ((n : ℝ) + 1) ^ 3) (Li3 x)

/-- Properties relied on (3), supplement g1:lem:prim: "$\frac d{ds}\Li_n(s)=\Li_{n-1}(s)/s$" at
`n = 3`, for `x < 1`, `x ≠ 0`. -/
def stmt_3 : Prop :=
  ∀ x : ℝ, x < 1 → x ≠ 0 → HasDerivAt Li3 (D00.Li2 x / x) x

/-- Properties relied on (4), main.tex 877: "$\zeta(n)=\Li_n(1)$" at `n = 3`: `riemannZeta 3`
is the real number `z3 = Σ_{n≥1} 1/n³`, and `Li3 1 = z3`. -/
def stmt_4 : Prop :=
  riemannZeta 3 = ((z3 : ℝ) : ℂ) ∧
    HasSum (fun n : ℕ => 1 / ((n : ℝ) + 1) ^ 3) z3 ∧ Li3 1 = z3

/-- Supplement, Lemma g1:lem:Mk at `k = 2` with `φ₂` from the proof of g1:lem:k2 (the analogue
of eq:M1t): `φ₂` is integrable on `[0,1]` and `M_2 = -½ ∫_0^1 φ₂(t) dt`. -/
def stmt_5 : Prop :=
  IntervalIntegrable phi2 volume 0 1 ∧
    W02.Mk 2 = (((-(1 / 2) * ∫ t in (0 : ℝ)..1, phi2 t) : ℝ) : ℂ)

/-- Supplement, Lemma g1:lem:k2: "$M_2=\int_b^a[\ell_1^2R_{20}+\ell_1\ell_2R_{11}]dv$", with the
integrability. -/
def stmt_6 : Prop :=
  IntervalIntegrable (fun v => l1 v ^ 2 * R20 v + l1 v * l2 v * R11 v) volume b a ∧
    W02.Mk 2 = (((∫ v in b..a, (l1 v ^ 2 * R20 v + l1 v * l2 v * R11 v)) : ℝ) : ℂ)

/-- Supplement, Lemma g1:lem:red2: "$M_2=W_3+T$ ... All integrals converge." -/
def stmt_7 : Prop :=
  IntervalIntegrable W3integrand volume b a ∧
    IntervalIntegrable (fun v => l1 v * Rwt v) volume b a ∧
    W02.Mk 2 = (((W3 + Tint) : ℝ) : ℂ)

/-- Supplement, Lemma g1:lem:T, and eq:T:
"$T=\dfrac{\pi^2+6\lambda^2}{24}=\dfrac{\zeta(2)-\Li_2(\frac12)}2$". -/
def stmt_8 : Prop :=
  Tint = (π ^ 2 + 6 * lam ^ 2) / 24 ∧
    (π ^ 2 + 6 * lam ^ 2) / 24 = (π ^ 2 / 6 - D00.Li2 (1 / 2)) / 2

/-- Supplement, table of g1:sub:M2:feq: (I3) for `z < 0`, (Tr) for `0 < z < 1`, (D3) for
`|w| ≤ 1`. -/
def stmt_9 : Prop :=
  (∀ z : ℝ, z < 0 →
      Li3 z - Li3 (1 / z) = -(π ^ 2 / 6) * Real.log (-z) - (1 / 6) * Real.log (-z) ^ 3) ∧
    (∀ z : ℝ, 0 < z → z < 1 →
      Li3 z + Li3 (1 - z) + Li3 (1 - 1 / z)
        = z3 + (1 / 6) * Real.log z ^ 3 + π ^ 2 / 6 * Real.log z
          - (1 / 2) * Real.log z ^ 2 * Real.log (1 - z)) ∧
    (∀ w : ℝ, |w| ≤ 1 → Li3 (w ^ 2) = 4 * Li3 w + 4 * Li3 (-w))

/-- Supplement, Lemma g1:lem:K: the Spence--Kummer equation (g1:K) for `0 < x < y < 1`. -/
def stmt_10 : Prop :=
  ∀ x y : ℝ, 0 < x → x < y → y < 1 →
    Li3 (x * (1 - y) ^ 2 / (y * (1 - x) ^ 2)) + Li3 (x * y) + Li3 (x / y)
      - 2 * Li3 (x * (1 - y) / (y * (1 - x))) - 2 * Li3 (x * (1 - y) / (x - 1))
      - 2 * Li3 ((1 - y) / (1 - x))
      - 2 * Li3 ((1 - y) / (y * (x - 1))) - 2 * Li3 x - 2 * Li3 y + 2 * z3
      = Real.log y ^ 2 * Real.log ((1 - y) / (1 - x)) - π ^ 2 / 3 * Real.log y
        - (1 / 3) * Real.log y ^ 3

/-- Supplement, proof of Theorem g1:thm:M2 (Lemmas g1:lem:raw, g1:lem:dilog, g1:lem:trilog
combined): the weight-three integral `W_3` equals (T1) minus `(π²+6λ²)/24`. -/
def stmt_11 : Prop := W3 = M2closed - (π ^ 2 + 6 * lam ^ 2) / 24

/-- Theorem thm:M2: "$M_2$ is given by \eqref{eq:M2} in Table~\ref{tab:closed}", for
`M_2 = W02.Mk 2`, the coefficient of `1/p³` in eq:expansion. -/
def stmt_12 : Prop := W02.Mk 2 = ((M2closed : ℝ) : ℂ)

/-- Theorem thm:M2: "Numerically $M_2=-0.02216191108726939942\ldots$", to all printed digits. -/
def stmt_13 : Prop :=
  ‖W02.Mk 2 - ((-0.02216191108726939942 : ℝ) : ℂ)‖ < (1 : ℝ) / 10 ^ 20

end Y01
-- END STATEMENTS

namespace Y01

open MeasureTheory Real Filter Topology Set

/-- The continuous extension of `Li₂(t)/t` (value `1` at `0`). -/
noncomputable def H (t : ℝ) : ℝ := if t = 0 then 1 else D00.Li2 t / t

lemma H_of_ne {t : ℝ} (h : t ≠ 0) : H t = D00.Li2 t / t := by simp [H, h]

lemma continuousAt_H_zero : ContinuousAt H 0 := by
  rw [← continuousWithinAt_compl_self, ContinuousWithinAt]
  have hd := D00.hasDerivAt_Li2_zero
  rw [hasDerivAt_iff_tendsto_slope] at hd
  have h0 : H 0 = 1 := by simp [H]
  rw [h0]
  apply hd.congr'
  filter_upwards [self_mem_nhdsWithin] with t ht
  rw [H_of_ne ht, slope_def_field]
  simp

lemma continuous_H : Continuous H := by
  rw [continuous_iff_continuousAt]
  intro t
  rcases eq_or_ne t 0 with rfl | h0
  · exact continuousAt_H_zero
  · have hc : ContinuousAt (fun s => D00.Li2 s / s) t :=
      (D00.continuous_Li2.continuousAt).div continuousAt_id h0
    apply hc.congr
    filter_upwards [isOpen_ne.mem_nhds h0] with s hs
    exact (H_of_ne hs).symm

lemma Li3_eq_H (x : ℝ) : Li3 x = ∫ t in (0:ℝ)..x, H t := by
  unfold Li3
  apply intervalIntegral.integral_congr_ae
  filter_upwards [(Set.countable_singleton (0:ℝ)).ae_notMem volume] with t ht _
  exact (H_of_ne ht).symm

lemma hasDerivAt_Li3_H (x : ℝ) : HasDerivAt Li3 (H x) x := by
  have e : Li3 = fun x => ∫ t in (0:ℝ)..x, H t := funext Li3_eq_H
  rw [e]
  exact intervalIntegral.integral_hasDerivAt_right (continuous_H.intervalIntegrable _ _)
    (continuous_H.stronglyMeasurableAtFilter _ _) continuous_H.continuousAt

lemma continuous_Li3 : Continuous Li3 :=
  continuous_iff_continuousAt.2 fun x => (hasDerivAt_Li3_H x).continuousAt

lemma hasDerivAt_Li3 {x : ℝ} (h0 : x ≠ 0) : HasDerivAt Li3 (D00.Li2 x / x) x := by
  have := hasDerivAt_Li3_H x
  rwa [H_of_ne h0] at this

lemma Li3_zero : Li3 0 = 0 := by simp [Li3]

lemma Li3_comp {f : ℝ → ℝ} {f' x : ℝ} (hf : HasDerivAt f f' x) (h0 : f x ≠ 0) :
    HasDerivAt (fun y => Li3 (f y)) (D00.Li2 (f x) / f x * f') x :=
  (hasDerivAt_Li3 h0).comp x hf

/-- The power series `Σ x^(n+1)/(n+1)³`. -/
noncomputable def S3 (x : ℝ) : ℝ := ∑' n : ℕ, x ^ (n + 1) / ((n : ℝ) + 1) ^ 3

lemma summable_inv_cube : Summable (fun n : ℕ => 1 / ((n : ℝ) + 1) ^ 3) := by
  have := (summable_nat_add_iff 1).mpr (Real.summable_one_div_nat_pow.mpr (by norm_num : 1 < 3))
  simpa using this

lemma term_bound3 {x : ℝ} (hx : |x| ≤ 1) (n : ℕ) :
    ‖x ^ (n + 1) / ((n : ℝ) + 1) ^ 3‖ ≤ 1 / ((n : ℝ) + 1) ^ 3 := by
  rw [Real.norm_eq_abs, abs_div, abs_pow, abs_of_pos (by positivity : (0:ℝ) < ((n:ℝ) + 1) ^ 3)]
  apply div_le_div_of_nonneg_right _ (by positivity)
  exact pow_le_one₀ (abs_nonneg x) hx

lemma summable_S3 {x : ℝ} (hx : |x| ≤ 1) :
    Summable (fun n : ℕ => x ^ (n + 1) / ((n : ℝ) + 1) ^ 3) :=
  Summable.of_norm_bounded summable_inv_cube (term_bound3 hx)

lemma continuousOn_S3 : ContinuousOn S3 (Icc (-1) 1) := by
  apply continuousOn_tsum (fun n => (continuous_pow (n + 1)).continuousOn.div_const _)
    summable_inv_cube
  intro n x hx
  exact term_bound3 (abs_le.mpr hx) n

lemma S3_zero : S3 0 = 0 := by simp [S3]

lemma hasDerivAt_S3 {y : ℝ} (hy : |y| < 1) (h0 : y ≠ 0) : HasDerivAt S3 (D00.Li2 y / y) y := by
  set r := (|y| + 1) / 2 with hr
  have hr1 : r < 1 := by rw [hr]; linarith
  have hyr : |y| < r := by rw [hr]; linarith
  have hr0 : 0 ≤ r := by rw [hr]; positivity
  have ht : y ∈ Ioo (-r) r := ⟨by linarith [neg_abs_le y], by linarith [le_abs_self y]⟩
  have hd := hasDerivAt_tsum_of_isPreconnected (t := Ioo (-r) r) (y₀ := 0)
    (g := fun (n : ℕ) (z : ℝ) => z ^ (n + 1) / ((n : ℝ) + 1) ^ 3)
    (g' := fun (n : ℕ) (z : ℝ) => z ^ n / ((n : ℝ) + 1) ^ 2) (u := fun n : ℕ => r ^ n)
    (summable_geometric_of_lt_one hr0 hr1) isOpen_Ioo isPreconnected_Ioo
    (by
      intro n z _
      have := (hasDerivAt_pow (n + 1) z).div_const (((n : ℝ) + 1) ^ 3)
      refine this.congr_deriv ?_
      push_cast
      field_simp)
    (by
      intro n z hz
      rw [Real.norm_eq_abs, abs_div, abs_pow,
        abs_of_pos (by positivity : (0:ℝ) < ((n:ℝ) + 1) ^ 2)]
      have hz' : |z| ≤ r := abs_le.mpr ⟨hz.1.le, hz.2.le⟩
      calc |z| ^ n / ((n:ℝ) + 1) ^ 2 ≤ |z| ^ n := by
            apply div_le_self (by positivity)
            have : (0:ℝ) ≤ n := Nat.cast_nonneg n
            nlinarith
        _ ≤ r ^ n := pow_le_pow_left₀ (abs_nonneg z) hz' n)
    ⟨by linarith [abs_nonneg y, hr], by linarith [abs_nonneg y, hr]⟩
    (by simp)
    ht
  have hsum : HasSum (fun n : ℕ => y ^ n / ((n : ℝ) + 1) ^ 2) (D00.Li2 y / y) := by
    have := (D00.hasSum_Li2 hy.le).div_const y
    convert this using 1
    funext n
    field_simp
    ring
  rw [hsum.tsum_eq] at hd
  exact hd

lemma Li3_eq_S3 {x : ℝ} (hx : |x| ≤ 1) : Li3 x = S3 x := by
  have key : (fun z => Li3 z - S3 z) 0 = (fun z => Li3 z - S3 z) x := by
    apply D00.eq_of_deriv_zero' (f := fun z => Li3 z - S3 z)
    · apply continuous_Li3.continuousOn.sub
      apply continuousOn_S3.mono
      intro z hz
      rw [mem_uIcc] at hz
      have hx' := abs_le.mp hx
      constructor <;> rcases hz with h | h <;> linarith [h.1, h.2]
    · intro z hz
      obtain ⟨hz0, hzx⟩ := D00.mem_uIoo_zero hz
      have := (hasDerivAt_Li3 hz0).sub (hasDerivAt_S3 (lt_of_lt_of_le hzx hx) hz0)
      rw [sub_self] at this
      exact this
  simp [S3_zero, Li3_zero] at key
  linarith

lemma Li3_hasSum {x : ℝ} (hx : |x| ≤ 1) :
    HasSum (fun n : ℕ => x ^ (n + 1) / ((n : ℝ) + 1) ^ 3) (Li3 x) := by
  rw [Li3_eq_S3 hx]
  exact (summable_S3 hx).hasSum

lemma zeta3_eq : riemannZeta 3 = ((∑' n : ℕ, 1 / (n : ℝ) ^ 3 : ℝ) : ℂ) := by
  rw [show (3:ℂ) = ((3:ℕ):ℂ) by norm_num, zeta_nat_eq_tsum_of_gt_one (by norm_num),
    Complex.ofReal_tsum]
  push_cast
  rfl

lemma z3_eq : z3 = ∑' n : ℕ, 1 / (n : ℝ) ^ 3 := by
  rw [z3, zeta3_eq, Complex.ofReal_re]

lemma z3_hasSum : HasSum (fun n : ℕ => 1 / ((n : ℝ) + 1) ^ 3) z3 := by
  have h := (Real.summable_one_div_nat_pow.mpr (by norm_num : 1 < 3)).hasSum
  rw [← z3_eq] at h
  have := (hasSum_nat_add_iff' 1).mpr h
  simpa using this

lemma Li3_one : Li3 1 = z3 := by
  have h := Li3_hasSum (x := 1) (by norm_num)
  simp only [one_pow] at h
  exact h.unique z3_hasSum

end Y01


namespace Y01

open MeasureTheory Real Filter Topology Set

lemma I3 {z : ℝ} (hz : z < 0) :
    Li3 z - Li3 (1 / z) = -(π ^ 2 / 6) * Real.log (-z) - (1 / 6) * Real.log (-z) ^ 3 := by
  set F : ℝ → ℝ := fun u => Li3 u - Li3 u⁻¹ + (π ^ 2 / 6) * Real.log (-u)
    + (1 / 6) * Real.log (-u) ^ 3 with hF
  have hd : ∀ u : ℝ, u < 0 → HasDerivAt F 0 u := by
    intro u hu
    have hu0 : u ≠ 0 := hu.ne
    have hnu : -u ≠ 0 := by linarith
    have d1 := hasDerivAt_Li3 hu0
    have d2 := Li3_comp (hasDerivAt_inv hu0) (inv_ne_zero hu0)
    have d3 := (hasDerivAt_neg u).log hnu
    have d4 := d3.pow 3
    have := ((d1.sub d2).add (d3.const_mul (π ^ 2 / 6))).add (d4.const_mul (1 / 6))
    refine this.congr_deriv ?_
    have inv := D00.inversion_neg hu
    rw [one_div] at inv
    have e : D00.Li2 u⁻¹ = -(π ^ 2 / 6) - (1 / 2) * Real.log (-u) ^ 2 - D00.Li2 u := by linarith
    rw [e]
    field_simp
    ring
  have key := D00.eq_of_deriv_zero' (f := F) (a := -1) (b := z)
    (fun u hu => by
      have : u < 0 := by
        rw [mem_uIcc] at hu; rcases hu with h | h <;> linarith [h.1, h.2]
      exact (hd u this).continuousAt.continuousWithinAt)
    (fun u hu => by
      have : u < 0 := by
        simp only [uIoo, mem_Ioo] at hu
        rcases le_total (-1) z with h | h
        · rw [max_eq_right h] at hu; linarith [hu.2]
        · rw [max_eq_left h] at hu; linarith [hu.2]
      exact hd u this)
  simp only [hF] at key
  norm_num at key
  rw [one_div, Real.log_neg_eq_log]
  linarith

lemma D3 {w : ℝ} (hw : |w| ≤ 1) : Li3 (w ^ 2) = 4 * Li3 w + 4 * Li3 (-w) := by
  have key := D00.eq_of_deriv_zero' (f := fun u => Li3 (u ^ 2) - 4 * Li3 u - 4 * Li3 (-u))
    (a := 0) (b := w)
    ((continuous_Li3.comp (continuous_pow 2)).sub (continuous_const.mul continuous_Li3)
      |>.sub (continuous_const.mul (continuous_Li3.comp continuous_neg))).continuousOn
    (by
      intro u hu
      obtain ⟨hu0, huw⟩ := D00.mem_uIoo_zero hu
      have hu1 : |u| ≤ 1 := (lt_of_lt_of_le huw hw).le
      have d1 := Li3_comp (hasDerivAt_pow 2 u) (pow_ne_zero 2 hu0)
      have d2 := hasDerivAt_Li3 hu0
      have d3 := Li3_comp (hasDerivAt_neg u) (neg_ne_zero.mpr hu0)
      have := (d1.sub (d2.const_mul 4)).sub (d3.const_mul 4)
      refine this.congr_deriv ?_
      have dup := D00.duplication hu1
      have e : D00.Li2 (u ^ 2) = 2 * D00.Li2 u + 2 * D00.Li2 (-u) := by linarith
      simp only [e]
      field_simp
      ring)
  simp [Li3_zero] at key
  linarith

lemma g_tendsto (z : ℝ) :
    Tendsto (fun u => Real.log u ^ 2 * Real.log (1 - u)) (𝓝[Icc z 1] 1) (𝓝 0) := by
  have hc : Continuous (fun u : ℝ => 4 * (1 - u) * |(1 - u) * Real.log (1 - u)|) :=
    (continuous_const.mul (continuous_const.sub continuous_id)).mul
      ((Real.continuous_mul_log.comp (continuous_const.sub continuous_id)).abs)
  have ht : Tendsto (fun u : ℝ => 4 * (1 - u) * |(1 - u) * Real.log (1 - u)|) (𝓝[Icc z 1] 1)
      (𝓝 0) := by
    have := (hc.tendsto 1).mono_left (nhdsWithin_le_nhds (s := Icc z 1))
    simpa using this
  apply squeeze_zero_norm' _ ht
  have h12 : ∀ᶠ u in 𝓝[Icc z 1] (1:ℝ), (1:ℝ) / 2 < u :=
    nhdsWithin_le_nhds (lt_mem_nhds (by norm_num))
  filter_upwards [self_mem_nhdsWithin, h12] with u hu h2
  have hu0 : 0 < u := by linarith
  have hu1 : u ≤ 1 := hu.2
  have hl1 : Real.log u ≤ 0 := Real.log_nonpos hu0.le hu1
  have hl2 : 1 - u⁻¹ ≤ Real.log u := Real.one_sub_inv_le_log_of_pos hu0
  have hinv : u⁻¹ ≤ 2 := by
    rw [inv_le_comm₀ hu0 (by norm_num)]; linarith
  have hl3 : -(2 * (1 - u)) ≤ Real.log u := by
    have : 1 - u⁻¹ = -((1 - u) * u⁻¹) := by
      rw [mul_comm, mul_sub, inv_mul_cancel₀ hu0.ne']; ring
    nlinarith
  have hsq : Real.log u ^ 2 ≤ 4 * (1 - u) ^ 2 := by nlinarith
  rw [Real.norm_eq_abs, abs_mul, abs_of_nonneg (sq_nonneg _), abs_mul,
    abs_of_nonneg (by linarith : (0:ℝ) ≤ 1 - u)]
  have := abs_nonneg (Real.log (1 - u))
  nlinarith

lemma Tr {z : ℝ} (hz0 : 0 < z) (hz1 : z < 1) :
    Li3 z + Li3 (1 - z) + Li3 (1 - 1 / z)
      = z3 + (1 / 6) * Real.log z ^ 3 + π ^ 2 / 6 * Real.log z
        - (1 / 2) * Real.log z ^ 2 * Real.log (1 - z) := by
  set F : ℝ → ℝ := fun u => Li3 u + Li3 (1 - u) + Li3 (1 - u⁻¹)
    - ((1 / 6) * Real.log u ^ 3 + π ^ 2 / 6 * Real.log u
        - (1 / 2) * (Real.log u ^ 2 * Real.log (1 - u))) with hF
  have hd : ∀ u : ℝ, 0 < u → u < 1 → HasDerivAt F 0 u := by
    intro u hu0 hu1
    have hu0' : u ≠ 0 := hu0.ne'
    have h1u : 1 - u ≠ 0 := by linarith
    have hA : 1 - u⁻¹ ≠ 0 := by
      intro h
      have : u⁻¹ = 1 := by linarith
      rw [inv_eq_one] at this; linarith
    have d1 := hasDerivAt_Li3 hu0'
    have d2 := Li3_comp ((hasDerivAt_id u).const_sub 1) h1u
    have d3 := Li3_comp ((hasDerivAt_inv hu0').const_sub 1) hA
    have dl := Real.hasDerivAt_log hu0'
    have dl1 := ((hasDerivAt_id u).const_sub 1).log h1u
    have := ((d1.add d2).add d3).sub
      (((dl.pow 3).const_mul (1 / 6) |>.add (dl.const_mul (π ^ 2 / 6))).sub
        (((dl.pow 2).mul dl1).const_mul (1 / 2)))
    refine this.congr_deriv ?_
    have hrefl := D00.reflection hu0 hu1
    have land := D00.landen (x := 1 - u) (by linarith)
    have e1 : (1 - u) / (1 - u - 1) = 1 - u⁻¹ := by
      rw [show (1:ℝ) - u - 1 = -u by ring, div_neg, sub_div, div_self hu0']; ring
    rw [e1, show 1 - (1 - u) = u by ring] at land
    have eA : D00.Li2 (1 - u⁻¹) = -(1 / 2) * Real.log u ^ 2 - D00.Li2 (1 - u) := by linarith
    have eB : D00.Li2 (1 - u) = π ^ 2 / 6 - Real.log u * Real.log (1 - u) - D00.Li2 u := by
      linarith
    simp only [id, Pi.pow_apply]
    rw [eA, eB]
    have hu1' : u - 1 ≠ 0 := by linarith
    field_simp
    ring
  have hc1 : ContinuousWithinAt F (Icc z 1) 1 := by
    have hA : ContinuousAt (fun u : ℝ => Li3 u + Li3 (1 - u) + Li3 (1 - u⁻¹)
        - ((1 / 6) * Real.log u ^ 3 + π ^ 2 / 6 * Real.log u)) 1 := by
      have hl : ContinuousAt Real.log 1 := Real.continuousAt_log one_ne_zero
      have hi : ContinuousAt (fun u : ℝ => u⁻¹) 1 := continuousAt_inv₀ one_ne_zero
      exact ((continuous_Li3.continuousAt.add
        (continuous_Li3.continuousAt.comp (continuousAt_const.sub continuousAt_id))).add
        (continuous_Li3.continuousAt.comp (continuousAt_const.sub hi))).sub
        (((hl.pow 3).const_mul _).add (hl.const_mul _))
    have hB : ContinuousWithinAt (fun u : ℝ => (1 / 2) * (Real.log u ^ 2 * Real.log (1 - u)))
        (Icc z 1) 1 := by
      rw [ContinuousWithinAt]
      have := (g_tendsto z).const_mul (1 / 2)
      simpa using this
    have := hA.continuousWithinAt.add hB
    convert this using 1
    funext u
    simp only [hF, Pi.add_apply]
    ring
  have key := D00.eq_of_deriv_zero (f := F) (a := z) (b := 1) hz1.le
    (fun u hu => by
      rcases eq_or_lt_of_le hu.2 with h | h
      · rw [h]; exact hc1
      · exact (hd u (by linarith [hu.1]) h).continuousAt.continuousWithinAt)
    (fun u hu => hd u (by linarith [hu.1]) hu.2)
  simp only [hF] at key
  norm_num [Li3_zero, Li3_one] at key
  rw [one_div]
  linarith

end Y01


namespace Y01

open MeasureTheory Real Filter Topology Set

lemma Li3_comp' {f : ℝ → ℝ} {f' x r : ℝ} (hf : HasDerivAt f f' x) (h0 : f x ≠ 0)
    (hr : f' = f x * r) : HasDerivAt (fun y => Li3 (f y)) (D00.Li2 (f x) * r) x := by
  refine (Li3_comp hf h0).congr_deriv ?_
  rw [hr]; field_simp

lemma Li2_comp' {f : ℝ → ℝ} {f' x r : ℝ} (hf : HasDerivAt f f' x) (h0 : f x ≠ 0)
    (h1 : f x ≠ 1) (hr : f' = f x * r) :
    HasDerivAt (fun y => D00.Li2 (f y)) (-Real.log (1 - f x) * r) x := by
  refine (D00.Li2_comp hf h1 h0).congr_deriv ?_
  rw [hr]; field_simp

/-! The nine arguments of the Spence--Kummer equation (the last two are `x`, `y`). -/
noncomputable def KA1 (x y : ℝ) : ℝ := x * (1 - y) ^ 2 / (y * (1 - x) ^ 2)
noncomputable def KA2 (x y : ℝ) : ℝ := x * y
noncomputable def KA3 (x y : ℝ) : ℝ := x / y
noncomputable def KA4 (x y : ℝ) : ℝ := x * (1 - y) / (y * (1 - x))
noncomputable def KA5 (x y : ℝ) : ℝ := x * (1 - y) / (x - 1)
noncomputable def KA6 (x y : ℝ) : ℝ := (1 - y) / (1 - x)
noncomputable def KA7 (x y : ℝ) : ℝ := (1 - y) / (y * (x - 1))

section KummerDerivs

variable {x y u v : ℝ}

/-! Derivatives in the first variable at `u` (logarithmic form). -/
lemma dx1 (hy : y ≠ 0) (hu0 : u ≠ 0) (hu1 : 1 - u ≠ 0) :
    HasDerivAt (fun s => KA1 s y) (KA1 u y * (1 / u + 2 / (1 - u))) u := by
  have h := ((hasDerivAt_id u).mul_const ((1 - y) ^ 2)).div
    ((((hasDerivAt_id u).const_sub 1).pow 2).const_mul y) (by simp [hy, hu1])
  refine h.congr_deriv ?_
  unfold KA1; try simp only [id, Pi.pow_apply]
  field_simp
  try ring

lemma dx2 (hu0 : u ≠ 0) : HasDerivAt (fun s => KA2 s y) (KA2 u y * (1 / u)) u := by
  refine ((hasDerivAt_id u).mul_const y).congr_deriv ?_
  unfold KA2; field_simp

lemma dx3 (hy : y ≠ 0) (hu0 : u ≠ 0) : HasDerivAt (fun s => KA3 s y) (KA3 u y * (1 / u)) u := by
  refine ((hasDerivAt_id u).div_const y).congr_deriv ?_
  unfold KA3; field_simp

lemma dx4 (hy : y ≠ 0) (hu0 : u ≠ 0) (hu1 : 1 - u ≠ 0) :
    HasDerivAt (fun s => KA4 s y) (KA4 u y * (1 / u + 1 / (1 - u))) u := by
  have h := ((hasDerivAt_id u).mul_const (1 - y)).div
    (((hasDerivAt_id u).const_sub 1).const_mul y) (by simp [hy, hu1])
  refine h.congr_deriv ?_
  unfold KA4; try simp only [id, Pi.pow_apply]
  field_simp
  try ring

lemma dx5 (hu0 : u ≠ 0) (hu1 : 1 - u ≠ 0) :
    HasDerivAt (fun s => KA5 s y) (KA5 u y * (1 / u + 1 / (1 - u))) u := by
  have hu1' : u - 1 ≠ 0 := by intro h; apply hu1; linarith
  have h := ((hasDerivAt_id u).mul_const (1 - y)).div
    ((hasDerivAt_id u).sub_const 1) (by simpa using hu1')
  refine h.congr_deriv ?_
  unfold KA5; try simp only [id, Pi.pow_apply]
  field_simp
  try ring

lemma dx6 (hu1 : 1 - u ≠ 0) :
    HasDerivAt (fun s => KA6 s y) (KA6 u y * (1 / (1 - u))) u := by
  have h := (hasDerivAt_const u (1 - y)).div
    ((hasDerivAt_id u).const_sub 1) (by simpa using hu1)
  refine h.congr_deriv ?_
  unfold KA6; try simp only [id, Pi.pow_apply]
  field_simp
  try ring

lemma dx7 (hy : y ≠ 0) (hu1 : 1 - u ≠ 0) :
    HasDerivAt (fun s => KA7 s y) (KA7 u y * (1 / (1 - u))) u := by
  have hu1' : u - 1 ≠ 0 := by intro h; apply hu1; linarith
  have h := (hasDerivAt_const u (1 - y)).div
    (((hasDerivAt_id u).sub_const 1).const_mul y) (by simp [hy, hu1'])
  refine h.congr_deriv ?_
  unfold KA7; try simp only [id, Pi.pow_apply]
  field_simp
  try ring

/-! Derivatives in the second variable at `v`. -/
lemma dy1 (hx : 1 - x ≠ 0) (hv0 : v ≠ 0) (hv1 : 1 - v ≠ 0) :
    HasDerivAt (fun s => KA1 x s) (KA1 x v * (-2 / (1 - v) - 1 / v)) v := by
  have h := ((((hasDerivAt_id v).const_sub 1).pow 2).const_mul x).div
    ((hasDerivAt_id v).mul_const ((1 - x) ^ 2)) (by simp [hv0, hx])
  refine h.congr_deriv ?_
  unfold KA1; try simp only [id, Pi.pow_apply]
  field_simp
  try ring

lemma dy2 (hv0 : v ≠ 0) : HasDerivAt (fun s => KA2 x s) (KA2 x v * (1 / v)) v := by
  refine ((hasDerivAt_id v).const_mul x).congr_deriv ?_
  unfold KA2; try simp only [id, Pi.pow_apply]
  field_simp

lemma dy3 (hv0 : v ≠ 0) : HasDerivAt (fun s => KA3 x s) (KA3 x v * (-1 / v)) v := by
  refine ((hasDerivAt_inv hv0).const_mul x).congr_deriv ?_
  unfold KA3; field_simp

lemma dy4 (hx : 1 - x ≠ 0) (hv0 : v ≠ 0) (hv1 : 1 - v ≠ 0) :
    HasDerivAt (fun s => KA4 x s) (KA4 x v * (-1 / (1 - v) - 1 / v)) v := by
  have h := (((hasDerivAt_id v).const_sub 1).const_mul x).div
    ((hasDerivAt_id v).mul_const (1 - x)) (by simp [hv0, hx])
  refine h.congr_deriv ?_
  unfold KA4; try simp only [id, Pi.pow_apply]
  field_simp
  try ring

lemma dy5 (hx : 1 - x ≠ 0) (hv1 : 1 - v ≠ 0) :
    HasDerivAt (fun s => KA5 x s) (KA5 x v * (-1 / (1 - v))) v := by
  have hx' : x - 1 ≠ 0 := by intro h; apply hx; linarith
  have h := (((hasDerivAt_id v).const_sub 1).const_mul x).div_const (x - 1)
  refine h.congr_deriv ?_
  unfold KA5; try simp only [id, Pi.pow_apply]
  field_simp

lemma dy6 (hx : 1 - x ≠ 0) (hv1 : 1 - v ≠ 0) :
    HasDerivAt (fun s => KA6 x s) (KA6 x v * (-1 / (1 - v))) v := by
  have h := ((hasDerivAt_id v).const_sub 1).div_const (1 - x)
  refine h.congr_deriv ?_
  unfold KA6; try simp only [id, Pi.pow_apply]
  field_simp

lemma dy7 (hx : 1 - x ≠ 0) (hv0 : v ≠ 0) (hv1 : 1 - v ≠ 0) :
    HasDerivAt (fun s => KA7 x s) (KA7 x v * (-1 / (1 - v) - 1 / v)) v := by
  have hx' : x - 1 ≠ 0 := by intro h; apply hx; linarith
  have h := ((hasDerivAt_id v).const_sub 1).div
    ((hasDerivAt_id v).mul_const (x - 1)) (by simp [hv0, hx'])
  refine h.congr_deriv ?_
  unfold KA7; try simp only [id, Pi.pow_apply]
  field_simp
  try ring

end KummerDerivs

end Y01

namespace Y01

open MeasureTheory Real Filter Topology Set

section KummerMain

variable {x y : ℝ}

lemma oneSub (hy : y ≠ 0) (hx : 1 - x ≠ 0) :
    1 - KA1 x y = (y - x) * (1 - x * y) / (y * (1 - x) ^ 2) ∧
    1 - KA2 x y = 1 - x * y ∧
    1 - KA3 x y = (y - x) / y ∧
    1 - KA4 x y = (y - x) / (y * (1 - x)) ∧
    1 - KA5 x y = (1 - x * y) / (1 - x) ∧
    1 - KA6 x y = (y - x) / (1 - x) ∧
    1 - KA7 x y = (1 - x * y) / (y * (1 - x)) := by
  have hx' : x - 1 ≠ 0 := by intro h; apply hx; linarith
  unfold KA1 KA2 KA3 KA4 KA5 KA6 KA7
  refine ⟨?_, rfl, ?_, ?_, ?_, ?_, ?_⟩ <;> field_simp <;> ring

lemma logs (hx0 : 0 < x) (hxy : x < y) (hy1 : y < 1) :
    Real.log (1 - KA1 x y) = Real.log (y - x) + Real.log (1 - x * y) - Real.log y
        - 2 * Real.log (1 - x) ∧
    Real.log (1 - KA2 x y) = Real.log (1 - x * y) ∧
    Real.log (1 - KA3 x y) = Real.log (y - x) - Real.log y ∧
    Real.log (1 - KA4 x y) = Real.log (y - x) - Real.log y - Real.log (1 - x) ∧
    Real.log (1 - KA5 x y) = Real.log (1 - x * y) - Real.log (1 - x) ∧
    Real.log (1 - KA6 x y) = Real.log (y - x) - Real.log (1 - x) ∧
    Real.log (1 - KA7 x y) = Real.log (1 - x * y) - Real.log y - Real.log (1 - x) := by
  have pa : 0 < y - x := by linarith
  have pb : 0 < 1 - x * y := by nlinarith
  have py : 0 < y := by linarith
  have px : 0 < 1 - x := by linarith
  obtain ⟨e1, e2, e3, e4, e5, e6, e7⟩ := oneSub py.ne' px.ne'
  rw [e1, e2, e3, e4, e5, e6, e7]
  refine ⟨?_, rfl, ?_, ?_, ?_, ?_, ?_⟩
  · rw [Real.log_div (mul_pos pa pb).ne' (mul_pos py (pow_pos px 2)).ne',
      Real.log_mul pa.ne' pb.ne', Real.log_mul py.ne' (pow_pos px 2).ne', Real.log_pow]
    push_cast; ring
  · rw [Real.log_div pa.ne' py.ne']
  · rw [Real.log_div pa.ne' (mul_pos py px).ne', Real.log_mul py.ne' px.ne']; ring
  · rw [Real.log_div pb.ne' px.ne']
  · rw [Real.log_div pa.ne' px.ne']
  · rw [Real.log_div pb.ne' (mul_pos py px).ne', Real.log_mul py.ne' px.ne']; ring

/-- The arguments `A_1,...,A_7` are nonzero and `≠ 1` on `0 < x < y < 1`. -/
lemma args_ne (hx0 : 0 < x) (hxy : x < y) (hy1 : y < 1) :
    (KA1 x y ≠ 0 ∧ KA2 x y ≠ 0 ∧ KA3 x y ≠ 0 ∧ KA4 x y ≠ 0 ∧ KA5 x y ≠ 0 ∧ KA6 x y ≠ 0 ∧
      KA7 x y ≠ 0) ∧
    (KA1 x y ≠ 1 ∧ KA2 x y ≠ 1 ∧ KA3 x y ≠ 1 ∧ KA4 x y ≠ 1 ∧ KA5 x y ≠ 1 ∧ KA6 x y ≠ 1 ∧
      KA7 x y ≠ 1) := by
  have pa : 0 < y - x := by linarith
  have pb : 0 < 1 - x * y := by nlinarith
  have py : 0 < y := by linarith
  have px : 0 < 1 - x := by linarith
  have p1y : 0 < 1 - y := by linarith
  have hx' : x - 1 ≠ 0 := by linarith
  obtain ⟨e1, e2, e3, e4, e5, e6, e7⟩ := oneSub py.ne' px.ne'
  refine ⟨⟨?_, ?_, ?_, ?_, ?_, ?_, ?_⟩, ?_, ?_, ?_, ?_, ?_, ?_, ?_⟩
  · unfold KA1; positivity
  · unfold KA2; positivity
  · unfold KA3; positivity
  · unfold KA4; positivity
  · unfold KA5; exact div_ne_zero (by positivity) hx'
  · unfold KA6; positivity
  · unfold KA7; exact div_ne_zero p1y.ne' (mul_ne_zero py.ne' hx')
  all_goals intro h
  · have : (0:ℝ) < 1 - KA1 x y := by rw [e1]; positivity
    linarith
  · nlinarith [h, show KA2 x y = x * y from rfl]
  · have : (0:ℝ) < 1 - KA3 x y := by rw [e3]; positivity
    linarith
  · have : (0:ℝ) < 1 - KA4 x y := by rw [e4]; positivity
    linarith
  · have : (0:ℝ) < 1 - KA5 x y := by rw [e5]; positivity
    linarith
  · have : (0:ℝ) < 1 - KA6 x y := by rw [e6]; positivity
    linarith
  · have : (0:ℝ) < 1 - KA7 x y := by rw [e7]; positivity
    linarith

/-- `B₁` of the proof of Lemma g1:lem:K. -/
noncomputable def KB1 (x y : ℝ) : ℝ :=
  D00.Li2 (KA1 x y) + D00.Li2 (KA2 x y) + D00.Li2 (KA3 x y) - 2 * D00.Li2 (KA4 x y)
    - 2 * D00.Li2 (KA5 x y) - 2 * D00.Li2 x

/-- `B₂` of the proof of Lemma g1:lem:K. -/
noncomputable def KB2 (x y : ℝ) : ℝ :=
  2 * (D00.Li2 (KA1 x y) - D00.Li2 (KA4 x y) - D00.Li2 (KA5 x y) - D00.Li2 (KA6 x y)
    - D00.Li2 (KA7 x y))

lemma diag (hx0 : 0 < x) (hx1 : x < 1) :
    KA1 x x = 1 ∧ KA2 x x = x ^ 2 ∧ KA3 x x = 1 ∧ KA4 x x = 1 ∧ KA5 x x = -x ∧ KA6 x x = 1 ∧
      KA7 x x = 1 / (-x) := by
  have h1 : 1 - x ≠ 0 := by linarith
  have h1' : x - 1 ≠ 0 := by linarith
  have h0 : x ≠ 0 := hx0.ne'
  unfold KA1 KA2 KA3 KA4 KA5 KA6 KA7
  refine ⟨?_, by ring, ?_, ?_, ?_, ?_, ?_⟩ <;> field_simp <;> ring

lemma KB_cont (hx1 : x < 1) {v : ℝ} (hv : 0 < v) :
    ContinuousAt (fun s => KB1 x s) v ∧ ContinuousAt (fun s => KB2 x s - Real.log s ^ 2) v := by
  have hL := D00.continuous_Li2
  have h1 : 1 - x ≠ 0 := by linarith
  have h1' : x - 1 ≠ 0 := by linarith
  have hv0 : v ≠ 0 := hv.ne'
  have a1 : v * (1 - x) ^ 2 ≠ 0 := by positivity
  have a2 : v * (1 - x) ≠ 0 := mul_ne_zero hv0 h1
  have a3 : v * (x - 1) ≠ 0 := mul_ne_zero hv0 h1'
  unfold KB1 KB2 KA1 KA2 KA3 KA4 KA5 KA6 KA7
  constructor <;> fun_prop (disch := assumption)

lemma KB1_zero (hx0 : 0 < x) (hxy : x < y) (hy1 : y < 1) : KB1 x y = 0 := by
  have hx1 : x < 1 := by linarith
  have key := D00.eq_of_deriv_zero (f := fun s => KB1 x s) (a := x) (b := y) hxy.le
    (fun v hv => ((KB_cont hx1 (by linarith [hv.1] : 0 < v)).1).continuousWithinAt)
    (by
      intro v hv
      have hv0 : v ≠ 0 := by linarith [hv.1]
      have hv1 : 1 - v ≠ 0 := by linarith [hv.2]
      have hx : 1 - x ≠ 0 := by linarith
      obtain ⟨⟨n1, n2, n3, n4, n5, -, -⟩, ⟨o1, o2, o3, o4, o5, -, -⟩⟩ :=
        args_ne hx0 hv.1 (by linarith [hv.2])
      obtain ⟨l1, l2, l3, l4, l5, -, -⟩ := logs hx0 hv.1 (by linarith [hv.2])
      have d1 := Li2_comp' (dy1 hx hv0 hv1) n1 o1 rfl
      have d2 := Li2_comp' (dy2 (x := x) hv0) n2 o2 rfl
      have d3 := Li2_comp' (dy3 (x := x) hv0) n3 o3 rfl
      have d4 := Li2_comp' (dy4 hx hv0 hv1) n4 o4 rfl
      have d5 := Li2_comp' (dy5 hx hv1) n5 o5 rfl
      have := ((((d1.add d2).add d3).sub (d4.const_mul 2)).sub (d5.const_mul 2)).sub
        (hasDerivAt_const v (2 * D00.Li2 x))
      refine this.congr_deriv ?_
      rw [l1, l2, l3, l4, l5]
      ring)
  have hd := diag hx0 hx1
  obtain ⟨e1, e2, e3, e4, e5, -, -⟩ := hd
  simp only [KB1, e1, e2, e3, e4, e5] at key
  have dup := D00.duplication (x := x) (by rw [abs_le]; constructor <;> linarith)
  unfold KB1
  linarith

lemma KB2_eq (hx0 : 0 < x) (hxy : x < y) (hy1 : y < 1) : KB2 x y = Real.log y ^ 2 := by
  have hx1 : x < 1 := by linarith
  have key := D00.eq_of_deriv_zero (f := fun s => KB2 x s - Real.log s ^ 2) (a := x) (b := y)
    hxy.le
    (fun v hv => ((KB_cont hx1 (by linarith [hv.1] : 0 < v)).2).continuousWithinAt)
    (by
      intro v hv
      have hv0 : v ≠ 0 := by linarith [hv.1]
      have hv1 : 1 - v ≠ 0 := by linarith [hv.2]
      have hx : 1 - x ≠ 0 := by linarith
      obtain ⟨⟨n1, -, -, n4, n5, n6, n7⟩, ⟨o1, -, -, o4, o5, o6, o7⟩⟩ :=
        args_ne hx0 hv.1 (by linarith [hv.2])
      obtain ⟨l1, -, -, l4, l5, l6, l7⟩ := logs hx0 hv.1 (by linarith [hv.2])
      have d1 := Li2_comp' (dy1 hx hv0 hv1) n1 o1 rfl
      have d4 := Li2_comp' (dy4 hx hv0 hv1) n4 o4 rfl
      have d5 := Li2_comp' (dy5 hx hv1) n5 o5 rfl
      have d6 := Li2_comp' (dy6 hx hv1) n6 o6 rfl
      have d7 := Li2_comp' (dy7 hx hv0 hv1) n7 o7 rfl
      have dl := (Real.hasDerivAt_log hv0).pow 2
      have := ((((((d1.sub d4).sub d5).sub d6).sub d7).const_mul 2).sub dl)
      refine this.congr_deriv ?_
      rw [l1, l4, l5, l6, l7]
      field_simp
      ring)
  have hd := diag hx0 hx1
  obtain ⟨e1, -, -, e4, e5, e6, e7⟩ := hd
  simp only [KB2, e1, e4, e5, e6, e7] at key
  have inv := D00.inversion_neg (x := -x) (by linarith)
  rw [neg_neg] at inv
  unfold KB2
  simp only [D00.Li2_one] at key
  linarith

/-- `Δ(x,y)` of Lemma g1:lem:K, as a function of `x` for fixed `y`. -/
noncomputable def KD (y s : ℝ) : ℝ :=
  Li3 (KA1 s y) + Li3 (KA2 s y) + Li3 (KA3 s y) - 2 * Li3 (KA4 s y) - 2 * Li3 (KA5 s y)
    - 2 * Li3 (KA6 s y) - 2 * Li3 (KA7 s y) - 2 * Li3 s - 2 * Li3 y + 2 * z3
    - (Real.log y ^ 2 * Real.log (KA6 s y) - π ^ 2 / 3 * Real.log y - (1 / 3) * Real.log y ^ 3)

lemma KD_cont (hy0 : 0 < y) (hy1 : y < 1) {u : ℝ} (hu : u < 1) : ContinuousAt (KD y) u := by
  have hL := continuous_Li3
  have h1 : 1 - u ≠ 0 := by linarith
  have h1' : u - 1 ≠ 0 := by linarith
  have hy : y ≠ 0 := hy0.ne'
  have a1 : y * (1 - u) ^ 2 ≠ 0 := by positivity
  have a2 : y * (1 - u) ≠ 0 := mul_ne_zero hy h1
  have a3 : y * (u - 1) ≠ 0 := mul_ne_zero hy h1'
  have a4 : (1 - y) / (1 - u) ≠ 0 := div_ne_zero (by linarith) h1
  unfold KD KA1 KA2 KA3 KA4 KA5 KA6 KA7
  fun_prop (disch := assumption)

lemma KD_zero_at (hy0 : 0 < y) (hy1 : y < 1) : KD y 0 = 0 := by
  have e7 : KA7 0 y = 1 - 1 / y := by unfold KA7; field_simp; ring
  have e6 : KA6 0 y = 1 - y := by unfold KA6; simp
  have z1 : KA1 0 y = 0 := by simp [KA1]
  have z2 : KA2 0 y = 0 := by simp [KA2]
  have z3' : KA3 0 y = 0 := by simp [KA3]
  have z4 : KA4 0 y = 0 := by simp [KA4]
  have z5 : KA5 0 y = 0 := by simp [KA5]
  have tr := Tr hy0 hy1
  simp only [KD, e6, e7, z1, z2, z3', z4, z5, Li3_zero]
  linarith

lemma KD_deriv (hy1 : y < 1) {u : ℝ} (hu0 : 0 < u) (huy : u < y) : HasDerivAt (KD y) 0 u := by
  have hy0 : y ≠ 0 := by linarith
  have hu0' : u ≠ 0 := hu0.ne'
  have hu1 : 1 - u ≠ 0 := by linarith
  obtain ⟨⟨n1, n2, n3, n4, n5, n6, n7⟩, -⟩ := args_ne hu0 huy hy1
  have d1 := Li3_comp' (dx1 hy0 hu0' hu1) n1 rfl
  have d2 := Li3_comp' (dx2 (y := y) hu0') n2 rfl
  have d3 := Li3_comp' (dx3 hy0 hu0') n3 rfl
  have d4 := Li3_comp' (dx4 hy0 hu0' hu1) n4 rfl
  have d5 := Li3_comp' (dx5 (y := y) hu0' hu1) n5 rfl
  have d6 := Li3_comp' (dx6 (y := y) hu1) n6 rfl
  have d7 := Li3_comp' (dx7 hy0 hu1) n7 rfl
  have d8 := hasDerivAt_Li3 hu0'
  have dl := (dx6 (y := y) hu1).log n6
  have := (((((((((d1.add d2).add d3).sub (d4.const_mul 2)).sub (d5.const_mul 2)).sub
    (d6.const_mul 2)).sub (d7.const_mul 2)).sub (d8.const_mul 2)).sub
    (hasDerivAt_const u (2 * Li3 y))).add (hasDerivAt_const u (2 * z3))).sub
    (((dl.const_mul (Real.log y ^ 2)).sub (hasDerivAt_const u (π ^ 2 / 3 * Real.log y))).sub
      (hasDerivAt_const u ((1 / 3) * Real.log y ^ 3)))
  refine this.congr_deriv ?_
  have hB1 := KB1_zero hu0 huy hy1
  have hB2 := KB2_eq hu0 huy hy1
  unfold KB1 at hB1
  unfold KB2 at hB2
  field_simp
  linear_combination (1 - u) * hB1 + u * hB2

lemma kummer (hx0 : 0 < x) (hxy : x < y) (hy1 : y < 1) : KD y x = 0 := by
  have hy0 : 0 < y := by linarith
  have key := D00.eq_of_deriv_zero (f := KD y) (a := 0) (b := x) hx0.le
    (fun u hu => (KD_cont hy0 hy1 (by linarith [hu.2])).continuousWithinAt)
    (fun u hu => KD_deriv hy1 hu.1 (by linarith [hu.2]))
  rw [← key]
  exact KD_zero_at hy0 hy1

end KummerMain

end Y01


/-!
# Part A of Y01: `stmt_5` (Lemma g1:lem:Mk at `k = 2`)

`M_2 = -F''(0)/4`, `F = 2 ∫_0^{1/2} gC`. We differentiate `F' = 2 ∫_0^{1/2} gCd` once more under
the integral sign. The domination comes from the Cauchy estimate applied to `gCd` (which is itself
bounded near `ε = 0` by `W02.gCd_bound`). At `ε = 0` the second `ε`-derivative of the integrand
is `2 φ₂(t)` on `(0, 1/2]`, and `φ₂` is symmetric under `t ↦ 1 - t`.
-/

namespace Y01.PA

open MeasureTheory Real Set Filter Topology

/-- explicit `ε`-derivative of `W02.gCd`. -/
noncomputable def gCdd (e : ℂ) (t : ℝ) : ℂ :=
  W02.aC e t * Complex.log (t : ℂ) * Complex.log (t : ℂ)
    + W02.bC e t * Complex.log ((1 - t : ℝ) : ℂ) * Complex.log ((1 - t : ℝ) : ℂ) -
   (W02.aC e t * Complex.log (t : ℂ) * Complex.log (t : ℂ) * (1 + W02.zC e t) ^ (1 / 2 : ℂ)
    + W02.aC e t * Complex.log (t : ℂ) * ((1 / 2 : ℂ) * (1 + W02.zC e t) ^ ((1 / 2 : ℂ) - 1) *
        (W02.zC e t * Complex.log ((t / (1 - t) : ℝ) : ℂ) * (-2)))
    + W02.aC e t * Complex.log (t : ℂ) * ((1 / 2 : ℂ) * (1 + W02.zC e t) ^ ((1 / 2 : ℂ) - 1) *
        (W02.zC e t * Complex.log ((t / (1 - t) : ℝ) : ℂ) * (-2)))
    + W02.aC e t * ((1 / 2 : ℂ) *
        (((1 / 2 : ℂ) - 1) * (1 + W02.zC e t) ^ ((1 / 2 : ℂ) - 1 - 1) *
            (W02.zC e t * Complex.log ((t / (1 - t) : ℝ) : ℂ) * (-2)) *
            (W02.zC e t * Complex.log ((t / (1 - t) : ℝ) : ℂ) * (-2))
          + (1 + W02.zC e t) ^ ((1 / 2 : ℂ) - 1) *
            (W02.zC e t * Complex.log ((t / (1 - t) : ℝ) : ℂ) * (-2) *
              Complex.log ((t / (1 - t) : ℝ) : ℂ) * (-2)))))

lemma gCdd_meas (ε : ℂ) : Measurable (gCdd ε) := by
  unfold gCdd W02.aC W02.bC W02.zC
  have := Complex.measurable_log
  fun_prop

lemma gCd_hasDerivAt {e : ℂ} (he : e.re < 1) {t : ℝ} (ht0 : 0 < t) (ht : t ≤ 1 / 2) :
    HasDerivAt (fun x => W02.gCd x t) (gCdd e t) e := by
  have h1 : 0 < 1 - t := by linarith
  have ha : HasDerivAt (fun x => W02.aC x t) (W02.aC e t * Complex.log (t : ℂ) * 1) e :=
    ((hasDerivAt_id e).sub_const 1).const_cpow (Or.inl (by exact_mod_cast ht0.ne'))
  have hb : HasDerivAt (fun x => W02.bC x t)
      (W02.bC e t * Complex.log ((1 - t : ℝ) : ℂ) * 1) e :=
    ((hasDerivAt_id e).sub_const 1).const_cpow (Or.inl (by exact_mod_cast h1.ne'))
  have hz : HasDerivAt (fun x => W02.zC x t)
      (W02.zC e t * Complex.log ((t / (1 - t) : ℝ) : ℂ) * (0 - 2 * 1)) e :=
    ((hasDerivAt_const e (2 : ℂ)).sub ((hasDerivAt_id e).const_mul 2)).const_cpow
      (Or.inl (by exact_mod_cast (div_pos ht0 h1).ne'))
  have hS := ((hasDerivAt_const e (1 : ℂ)).add hz).cpow_const (c := (1 / 2 : ℂ))
    (W02.zC_slit he ht0 ht)
  have hT := ((hasDerivAt_const e (1 : ℂ)).add hz).cpow_const (c := (1 / 2 : ℂ) - 1)
    (W02.zC_slit he ht0 ht)
  have h := ((ha.mul_const (Complex.log (t : ℂ))).add
      (hb.mul_const (Complex.log ((1 - t : ℝ) : ℂ)))).sub
    (((ha.mul_const (Complex.log (t : ℂ))).mul hS).add
      (ha.mul ((hT.const_mul (1 / 2 : ℂ)).mul
        ((hz.mul_const (Complex.log ((t / (1 - t) : ℝ) : ℂ))).mul_const (-2)))))
  unfold W02.gCd
  convert h using 1
  unfold gCdd
  simp only [Pi.add_apply, Pi.mul_apply]
  ring

/-- the uniform bound of `W02.gCd` on the ball of radius `1/3` around `0`. -/
noncomputable def Bd : ℝ :=
  (2 : ℝ) ^ (2 - (0 : ℂ).re + 2 * ((1 - (0 : ℂ).re) / 3)) / ((1 - (0 : ℂ).re) / 3)

lemma re_lt_of_ball {y : ℂ} (hy : y ∈ Metric.ball (0 : ℂ) (1 / 3)) : y.re < 1 := by
  rw [mem_ball_zero_iff] at hy
  have := Complex.re_le_norm y
  linarith

lemma gCd_bound0 {y : ℂ} (hy : y ∈ Metric.ball (0 : ℂ) (1 / 3)) {t : ℝ} (ht0 : 0 < t)
    (ht : t ≤ 1 / 2) : ‖W02.gCd y t‖ ≤ Bd := by
  have hy3 : y ∈ Metric.ball (0 : ℂ) ((1 - (0 : ℂ).re) / 3) := by
    simpa using hy
  exact W02.gCd_bound (x₀ := 0) (by simp) hy3 ht0 ht

lemma gCdd_bound {x : ℂ} (hx : x ∈ Metric.ball (0 : ℂ) (1 / 6)) {t : ℝ} (ht0 : 0 < t)
    (ht : t ≤ 1 / 2) : ‖gCdd x t‖ ≤ Bd / (1 / 6) := by
  have hsub : Metric.closedBall x (1 / 6) ⊆ Metric.ball (0 : ℂ) (1 / 3) := by
    intro y hy
    rw [Metric.mem_closedBall, dist_eq_norm] at hy
    rw [mem_ball_zero_iff] at hx ⊢
    calc ‖y‖ = ‖(y - x) + x‖ := by ring_nf
      _ ≤ ‖y - x‖ + ‖x‖ := norm_add_le _ _
      _ < 1 / 3 := by linarith
  have hdiff : DifferentiableOn ℂ (fun e => W02.gCd e t) (Metric.ball (0 : ℂ) (1 / 3)) :=
    fun e he => (gCd_hasDerivAt (re_lt_of_ball he) ht0 ht).differentiableAt.differentiableWithinAt
  have hx3 : x ∈ Metric.ball (0 : ℂ) (1 / 3) := hsub (Metric.mem_closedBall_self (by norm_num))
  rw [← (gCd_hasDerivAt (re_lt_of_ball hx3) ht0 ht).deriv]
  apply Complex.norm_deriv_le_of_forall_mem_sphere_norm_le (by norm_num)
    (hdiff.diffContOnCl_ball hsub)
  intro y hy
  exact gCd_bound0 (hsub (Metric.sphere_subset_closedBall hy)) ht0 ht

lemma gCd_iInt0 : IntervalIntegrable (W02.gCd 0) volume 0 (1 / 2) := by
  refine IntervalIntegrable.mono_fun' (g := fun _ => Bd)
    intervalIntegrable_const (W02.gCd_meas 0).aestronglyMeasurable ?_
  rw [uIoc_of_le (by norm_num)]
  exact ae_restrict_of_forall_mem measurableSet_Ioc fun t ht =>
    gCd_bound0 (by simp) ht.1 ht.2

lemma key2 : IntervalIntegrable (gCdd 0) volume 0 (1 / 2) ∧
    HasDerivAt (fun x => ∫ t in (0 : ℝ)..(1 / 2), W02.gCd x t)
      (∫ t in (0 : ℝ)..(1 / 2), gCdd 0 t) 0 := by
  have h6 : Metric.ball (0 : ℂ) (1 / 6) ⊆ Metric.ball (0 : ℂ) (1 / 3) :=
    Metric.ball_subset_ball (by norm_num)
  exact intervalIntegral.hasDerivAt_integral_of_dominated_loc_of_deriv_le
    (μ := volume) (a := 0) (b := 1 / 2) (F := fun e t => W02.gCd e t) (F' := gCdd)
    (x₀ := 0) (s := Metric.ball (0 : ℂ) (1 / 6)) (bound := fun _ => Bd / (1 / 6))
    (Metric.ball_mem_nhds 0 (by norm_num))
    (Filter.Eventually.of_forall fun x => (W02.gCd_meas x).aestronglyMeasurable)
    gCd_iInt0 ((gCdd_meas 0).aestronglyMeasurable)
    (Filter.Eventually.of_forall fun t ht x hx => by
      rw [uIoc_of_le (by norm_num)] at ht
      exact gCdd_bound hx ht.1 ht.2)
    intervalIntegrable_const
    (Filter.Eventually.of_forall fun t ht x hx => by
      rw [uIoc_of_le (by norm_num)] at ht
      exact gCd_hasDerivAt (re_lt_of_ball (h6 hx)) ht.1 ht.2)

lemma iter2 : iteratedDeriv 2 W02.FC 0 = 2 * ∫ t in (0 : ℝ)..(1 / 2), gCdd 0 t := by
  have hev : deriv W02.FC =ᶠ[𝓝 0] fun x => 2 * ∫ t in (0 : ℝ)..(1 / 2), W02.gCd x t := by
    have ho : IsOpen {e : ℂ | e.re < 1} := isOpen_lt Complex.continuous_re continuous_const
    filter_upwards [ho.mem_nhds (by simp : (0 : ℂ) ∈ {e : ℂ | e.re < 1})] with x hx
    exact (W02.FC_hasDerivAt hx).deriv
  rw [iteratedDeriv_succ, iteratedDeriv_one, hev.deriv_eq]
  exact (key2.2.const_mul (2 : ℂ)).deriv

lemma alg (t s q x y : ℝ) (ht : t ≠ 0) (hs : s ≠ 0) (hq : q ≠ 0)
    (hq2 : q ^ 2 = t ^ 2 + s ^ 2) :
    t⁻¹ * x * x + s⁻¹ * y * y -
      (t⁻¹ * x * x * (q / s)
        + t⁻¹ * x * ((1 / 2) * (q / s)⁻¹ * ((t / s) ^ 2 * (x - y) * (-2)))
        + t⁻¹ * x * ((1 / 2) * (q / s)⁻¹ * ((t / s) ^ 2 * (x - y) * (-2)))
        + t⁻¹ * ((1 / 2) *
          (((1 / 2) - 1) * ((q / s) ^ 3)⁻¹ * ((t / s) ^ 2 * (x - y) * (-2)) *
              ((t / s) ^ 2 * (x - y) * (-2))
            + (q / s)⁻¹ * ((t / s) ^ 2 * (x - y) * (-2) * (x - y) * (-2))))) =
    2 * (1 / (2 * t) * (1 - 2 * (s / q) + (s / q) ^ 3) * x ^ 2
      + 1 / (2 * s) * (1 - 2 * (t / q) + (t / q) ^ 3) * y ^ 2
      + t * s / q ^ 3 * x * y) := by
  have key : t⁻¹ * x * x + s⁻¹ * y * y -
      (t⁻¹ * x * x * (q / s)
        + t⁻¹ * x * ((1 / 2) * (q / s)⁻¹ * ((t / s) ^ 2 * (x - y) * (-2)))
        + t⁻¹ * x * ((1 / 2) * (q / s)⁻¹ * ((t / s) ^ 2 * (x - y) * (-2)))
        + t⁻¹ * ((1 / 2) *
          (((1 / 2) - 1) * ((q / s) ^ 3)⁻¹ * ((t / s) ^ 2 * (x - y) * (-2)) *
              ((t / s) ^ 2 * (x - y) * (-2))
            + (q / s)⁻¹ * ((t / s) ^ 2 * (x - y) * (-2) * (x - y) * (-2))))) -
    2 * (1 / (2 * t) * (1 - 2 * (s / q) + (s / q) ^ 3) * x ^ 2
      + 1 / (2 * s) * (1 - 2 * (t / q) + (t / q) ^ 3) * y ^ 2
      + t * s / q ^ 3 * x * y) =
      ((-(x ^ 2 * q ^ 2) / (t * s) + s * x ^ 2 / t + 2 * t * x * y / s - t * x ^ 2 / s) / q ^ 3) *
        (q ^ 2 - (t ^ 2 + s ^ 2)) := by
    field_simp
    ring
  linear_combination key +
    ((-(x ^ 2 * q ^ 2) / (t * s) + s * x ^ 2 / t + 2 * t * x * y / s - t * x ^ 2 / s) / q ^ 3) * hq2

lemma cpow_m32 (z : ℝ) (hz : 0 ≤ z) :
    (1 + (z : ℂ)) ^ ((1 / 2 : ℂ) - 1 - 1) = ((((Real.sqrt (1 + z)) ^ 3)⁻¹ : ℝ) : ℂ) := by
  have h1 : 0 ≤ 1 + z := by linarith
  rw [show (1 / 2 : ℂ) - 1 - 1 = ((-(3 / 2)) : ℝ) by push_cast; ring,
    ← Complex.ofReal_one, ← Complex.ofReal_add, ← Complex.ofReal_cpow h1]
  congr 1
  rw [Real.rpow_neg h1, Real.sqrt_eq_rpow, ← Real.rpow_natCast, ← Real.rpow_mul h1]
  norm_num

lemma rc {t : ℝ} (ht0 : 0 < t) (ht1 : t < 1) :
    W02.aR 0 t * Real.log t * Real.log t + W02.bR 0 t * Real.log (1 - t) * Real.log (1 - t) -
      (W02.aR 0 t * Real.log t * Real.log t *
          Real.sqrt (1 + (t / (1 - t)) ^ (2 - 2 * (0:ℝ)))
        + W02.aR 0 t * Real.log t * ((1 / 2) *
            (Real.sqrt (1 + (t / (1 - t)) ^ (2 - 2 * (0:ℝ))))⁻¹ *
            ((t / (1 - t)) ^ (2 - 2 * (0:ℝ)) * Real.log (t / (1 - t)) * (-2)))
        + W02.aR 0 t * Real.log t * ((1 / 2) *
            (Real.sqrt (1 + (t / (1 - t)) ^ (2 - 2 * (0:ℝ))))⁻¹ *
            ((t / (1 - t)) ^ (2 - 2 * (0:ℝ)) * Real.log (t / (1 - t)) * (-2)))
        + W02.aR 0 t * ((1 / 2) *
          (((1 / 2) - 1) * ((Real.sqrt (1 + (t / (1 - t)) ^ (2 - 2 * (0:ℝ)))) ^ 3)⁻¹ *
              ((t / (1 - t)) ^ (2 - 2 * (0:ℝ)) * Real.log (t / (1 - t)) * (-2)) *
              ((t / (1 - t)) ^ (2 - 2 * (0:ℝ)) * Real.log (t / (1 - t)) * (-2))
            + (Real.sqrt (1 + (t / (1 - t)) ^ (2 - 2 * (0:ℝ))))⁻¹ *
              ((t / (1 - t)) ^ (2 - 2 * (0:ℝ)) * Real.log (t / (1 - t)) * (-2) *
                Real.log (t / (1 - t)) * (-2))))) =
      2 * phi2 t := by
  have h1 : 0 < 1 - t := by linarith
  have hq := W02.q_pos t
  have hq2 := W02.q_sq t
  have ha : W02.aR 0 t = t⁻¹ := by
    unfold W02.aR; rw [show (0:ℝ) - 1 = -1 by ring, Real.rpow_neg_one]
  have hb : W02.bR 0 t = (1 - t)⁻¹ := by
    unfold W02.bR; rw [show (0:ℝ) - 1 = -1 by ring, Real.rpow_neg_one]
  have hz : (t / (1 - t)) ^ (2 - 2 * (0:ℝ)) = (t / (1 - t)) ^ 2 := by
    rw [show (2:ℝ) - 2 * 0 = ((2:ℕ):ℝ) by norm_num, Real.rpow_natCast]
  have hs : Real.sqrt (1 + (t / (1 - t)) ^ 2) = W02.q t / (1 - t) := by
    rw [show 1 + (t / (1 - t)) ^ 2 = (W02.q t / (1 - t)) ^ 2 by
      rw [div_pow, div_pow, hq2]; field_simp; ring]
    exact Real.sqrt_sq (by positivity)
  have hl : Real.log (t / (1 - t)) = Real.log t - Real.log (1 - t) :=
    Real.log_div ht0.ne' h1.ne'
  have hA := alg t (1 - t) (W02.q t) (Real.log t) (Real.log (1 - t)) ht0.ne' h1.ne' hq.ne' hq2
  rw [ha, hb, hz, hs, hl]
  unfold phi2 c c'
  linear_combination hA

lemma gCdd_zero {t : ℝ} (ht0 : 0 < t) (ht : t ≤ 1 / 2) :
    gCdd 0 t = (((2 * phi2 t : ℝ)) : ℂ) := by
  have h1 : t < 1 := by linarith
  have h1' : 0 < 1 - t := by linarith
  have hz0 : 0 ≤ (t / (1 - t)) ^ (2 - 2 * (0:ℝ)) :=
    Real.rpow_nonneg (div_nonneg ht0.le h1'.le) _
  have hc := rc ht0 h1
  rw [← Complex.ofReal_zero]
  unfold gCdd
  rw [W02.aC_ofReal 0 ht0.le, W02.bC_ofReal 0 h1.le, W02.zC_ofReal 0 ht0.le h1,
    W02.sqrt_cpow _ hz0, W02.sqrt_inv_cpow _ hz0, cpow_m32 _ hz0,
    ← Complex.ofReal_log ht0.le, ← Complex.ofReal_log h1'.le,
    ← Complex.ofReal_log (div_nonneg ht0.le h1'.le)]
  have hc' := congrArg Complex.ofReal hc
  push_cast at hc' ⊢
  linear_combination hc'

lemma phi2_symm (t : ℝ) : phi2 (1 - t) = phi2 t := by
  unfold phi2 c c'
  rw [W02.q_symm, sub_sub_cancel]
  ring

lemma phi2_meas : Measurable phi2 := by
  unfold phi2 c c' W02.q
  fun_prop

lemma phi2_bound_half {t : ℝ} (ht0 : 0 < t) (ht : t ≤ 1 / 2) :
    ‖phi2 t‖ ≤ Bd / (1 / 6) / 2 := by
  have h := gCdd_bound (x := 0) (by simp) ht0 ht
  rw [gCdd_zero ht0 ht, Complex.norm_real, Real.norm_eq_abs, abs_mul] at h
  rw [Real.norm_eq_abs]
  norm_num at h ⊢
  linarith

lemma phi2_bound {t : ℝ} (ht0 : 0 < t) (ht1 : t < 1) : ‖phi2 t‖ ≤ Bd / (1 / 6) / 2 := by
  rcases le_total t (1 / 2) with h | h
  · exact phi2_bound_half ht0 h
  · rw [← phi2_symm t]
    exact phi2_bound_half (by linarith) (by linarith)

lemma phi2_iInt : IntervalIntegrable phi2 volume 0 1 := by
  refine IntervalIntegrable.mono_fun' (g := fun _ => Bd / (1 / 6) / 2)
    intervalIntegrable_const phi2_meas.aestronglyMeasurable ?_
  rw [uIoc_of_le zero_le_one, ← restrict_Ioo_eq_restrict_Ioc]
  exact ae_restrict_of_forall_mem measurableSet_Ioo fun t ht => phi2_bound ht.1 ht.2

lemma phi2_half : ∫ t in (0:ℝ)..1, phi2 t = 2 * ∫ t in (0:ℝ)..(1/2), phi2 t := by
  have hi := phi2_iInt
  have h1 : IntervalIntegrable phi2 volume 0 (1/2) :=
    hi.mono_set (by rw [uIcc_of_le (by norm_num), uIcc_of_le (by norm_num)]
                    exact Icc_subset_Icc le_rfl (by norm_num))
  have h2 : IntervalIntegrable phi2 volume (1/2) 1 :=
    hi.mono_set (by rw [uIcc_of_le (by norm_num), uIcc_of_le (by norm_num)]
                    exact Icc_subset_Icc (by norm_num) le_rfl)
  rw [← intervalIntegral.integral_add_adjacent_intervals h1 h2]
  have : ∫ t in (1/2:ℝ)..1, phi2 t = ∫ t in (0:ℝ)..(1/2), phi2 t := by
    have := intervalIntegral.integral_comp_sub_left phi2 (a := 0) (b := 1/2) 1
    simp only [phi2_symm] at this
    rw [this]; norm_num
  rw [this]; ring

theorem main : Y01.stmt_5 := by
  refine ⟨phi2_iInt, ?_⟩
  unfold W02.Mk
  rw [iter2]
  have : ∫ t in (0:ℝ)..(1 / 2), gCdd 0 t = ∫ t in (0:ℝ)..(1 / 2), (((2 * phi2 t : ℝ)) : ℂ) := by
    refine intervalIntegral.integral_congr_ae (Filter.Eventually.of_forall fun t ht => ?_)
    rw [uIoc_of_le (by norm_num)] at ht
    exact gCdd_zero ht.1 ht.2
  rw [this, intervalIntegral.integral_ofReal, intervalIntegral.integral_const_mul, phi2_half]
  simp only [Nat.factorial_two, Nat.cast_ofNat]
  push_cast
  ring

end Y01.PA

namespace Y01.PB

open MeasureTheory Real Set

lemma s_sq : Real.sqrt 2 ^ 2 = 2 := Real.sq_sqrt (by norm_num)
lemma s_gt : (1.41 : ℝ) < Real.sqrt 2 := by
  rw [Real.lt_sqrt (by norm_num)]; norm_num
lemma s_pos : 0 < Real.sqrt 2 := by linarith [s_gt]

/-- the `ℓ₁²` part of `φ₂`. -/
noncomputable def f (t : ℝ) : ℝ := 1 / (2 * t) * (1 - 2 * c t + c t ^ 3) * Real.log t ^ 2

/-- the `ℓ₁ℓ₂` part of `φ₂`. -/
noncomputable def g (t : ℝ) : ℝ := t * (1 - t) / W02.q t ^ 3 * Real.log t * Real.log (1 - t)

lemma c_one_sub (t : ℝ) : c (1 - t) = c' t := by
  unfold c c'; rw [W02.q_symm]; ring

lemma phi2_eq (t : ℝ) : phi2 t = f t + f (1 - t) + g t := by
  unfold phi2 f g
  rw [c_one_sub]

lemma q_cont : Continuous W02.q := by
  unfold W02.q; fun_prop

lemma g_eq (t : ℝ) : g t = (t * Real.log t) * ((1 - t) * Real.log (1 - t)) / W02.q t ^ 3 := by
  unfold g; ring

lemma g_cont : Continuous g := by
  have h : g = fun t => (t * Real.log t) * ((1 - t) * Real.log (1 - t)) / W02.q t ^ 3 :=
    funext g_eq
  rw [h]
  refine Continuous.div ?_ (q_cont.pow 3) (fun t => (pow_pos (W02.q_pos t) 3).ne')
  exact Real.continuous_mul_log.mul (Real.continuous_mul_log.comp (continuous_const.sub continuous_id))

lemma den_pos (t : ℝ) : 0 < 1 - t + W02.q t := by
  have hq := W02.q_pos t
  have hq2 := W02.q_sq t
  by_contra h
  push Not at h
  nlinarith [sq_nonneg t]

lemma c_sub_one (t : ℝ) : c t - 1 = -t ^ 2 / (W02.q t * (1 - t + W02.q t)) := by
  have hq := W02.q_pos t
  have hq2 := W02.q_sq t
  have hd := den_pos t
  unfold c
  rw [div_sub_one hq.ne', div_eq_div_iff hq.ne' (by positivity)]
  linear_combination (-W02.q t) * hq2

lemma f_eq {t : ℝ} (ht : 0 ≤ t) :
    f t = -(4 * (Real.sqrt t * Real.log (Real.sqrt t)) ^ 2) * (c t ^ 2 + c t - 1)
      / (2 * W02.q t * (1 - t + W02.q t)) := by
  rcases eq_or_lt_of_le ht with h | h
  · subst h; simp [f]
  have hq := W02.q_pos t
  have hd := den_pos t
  have h1 : 1 - 2 * c t + c t ^ 3 = (c t - 1) * (c t ^ 2 + c t - 1) := by ring
  unfold f
  rw [h1, c_sub_one, Real.log_sqrt ht, mul_pow, Real.sq_sqrt ht]
  field_simp
  ring

lemma c_cont : Continuous c := by
  unfold c
  exact Continuous.div (continuous_const.sub continuous_id) q_cont (fun t => (W02.q_pos t).ne')

lemma f_cont : ContinuousOn f (Icc 0 1) := by
  have hc : Continuous fun t => -(4 * (Real.sqrt t * Real.log (Real.sqrt t)) ^ 2) *
      (c t ^ 2 + c t - 1) / (2 * W02.q t * (1 - t + W02.q t)) := by
    refine Continuous.div ?_ ?_ ?_
    · exact ((((Real.continuous_mul_log.comp Real.continuous_sqrt).pow 2).const_mul 4).neg).mul
        ((c_cont.pow 2).add c_cont |>.sub continuous_const)
    · exact (continuous_const.mul q_cont).mul ((continuous_const.sub continuous_id).add q_cont)
    · intro t
      exact (mul_pos (mul_pos two_pos (W02.q_pos t)) (den_pos t)).ne'
  exact hc.continuousOn.congr (fun t ht => f_eq ht.1)

lemma tpoly_pos {v r : ℝ} (hr : r ^ 2 = 2) (hr1 : 1 < r) (hb : r - 1 < v) :
    0 < v ^ 2 + 2 * v - 1 := by
  have : v ^ 2 + 2 * v - 1 = (v - (r - 1)) * (v + (1 + r)) := by linear_combination hr
  rw [this]; apply mul_pos <;> linarith

lemma idA {v r : ℝ} (hr : r ^ 2 = 2) (hr1 : 1 < r) (hv : 0 < v) (hb : r - 1 < v) :
    r * v / (v ^ 2 + 1) ^ 2 - 1 / (v + (1 + r)) + (4 - 3 * r) / (8 * v)
      = -(1 / (2 * ((v ^ 2 + 2 * v - 1) / (4 * v))) *
          (1 - 2 * (r * (1 + 2 * v - v ^ 2) / (2 * (v ^ 2 + 1)))
            + (r * (1 + 2 * v - v ^ 2) / (2 * (v ^ 2 + 1))) ^ 3)) * ((v ^ 2 + 1) / (4 * v ^ 2)) := by
  have h1 := tpoly_pos hr hr1 hb
  have h2 : 0 < v + (1 + r) := by linarith
  have h3 : 0 < v ^ 2 + 1 := by positivity
  have h1' : v * (v + 2) - 1 ≠ 0 := by nlinarith
  field_simp
  linear_combination ((8:ℝ)*r^1 + (8:ℝ)*r^2 + (-128:ℝ)*v^1 + (56:ℝ)*v^1*r^1 + (48:ℝ)*v^1*r^2 + (120:ℝ)*v^2*r^1 + (72:ℝ)*v^2*r^2 + (-256:ℝ)*v^3 + (40:ℝ)*v^3*r^1 + (-32:ℝ)*v^3*r^2 + (-104:ℝ)*v^4*r^1 + (-72:ℝ)*v^4*r^2 + (-128:ℝ)*v^5 + (-24:ℝ)*v^5*r^1 + (48:ℝ)*v^5*r^2 + (40:ℝ)*v^6*r^1 + (-8:ℝ)*v^6*r^2 + (-8:ℝ)*v^7*r^1) * hr

lemma idB {v r : ℝ} (hr : r ^ 2 = 2) (hr1 : 1 < r) (hv : 0 < v) :
    -(r * v / (v ^ 2 + 1) ^ 2) + r / (8 * v)
      = -(1 / 2) * (((v ^ 2 + 2 * v - 1) / (4 * v)) * (1 - (v ^ 2 + 2 * v - 1) / (4 * v))
          / (r * (v ^ 2 + 1) / (4 * v)) ^ 3) * ((v ^ 2 + 1) / (4 * v ^ 2)) := by
  have h3 : 0 < v ^ 2 + 1 := by positivity
  have h4 : 0 < r := by linarith
  field_simp
  linear_combination ((4:ℝ) + (2:ℝ)*r^2 + (-24:ℝ)*v^2 + (-12:ℝ)*v^2*r^2 + (4:ℝ)*v^4 + (2:ℝ)*v^4*r^2) * hr

lemma q_tv {v : ℝ} (hv : 0 < v) : W02.q (tv v) = Real.sqrt 2 * (v ^ 2 + 1) / (4 * v) := by
  have h : tv v ^ 2 + (1 - tv v) ^ 2 = (Real.sqrt 2 * (v ^ 2 + 1) / (4 * v)) ^ 2 := by
    unfold tv
    field_simp
    nlinarith [s_sq]
  unfold W02.q
  rw [h, Real.sqrt_sq (by have := s_pos; positivity)]

lemma c_tv {v : ℝ} (hv : 0 < v) :
    c (tv v) = Real.sqrt 2 * (1 + 2 * v - v ^ 2) / (2 * (v ^ 2 + 1)) := by
  have hs := s_pos
  unfold c
  rw [q_tv hv]
  unfold tv
  have h3 : 0 < v ^ 2 + 1 := by positivity
  field_simp
  rw [s_sq]; ring

lemma b_eq : b = Real.sqrt 2 - 1 := rfl
lemma a_eq : a = 1 + Real.sqrt 2 := rfl

lemma hr1 : (1 : ℝ) < Real.sqrt 2 := by linarith [s_gt]

lemma R20_eq {v : ℝ} (hv : 0 < v) (hb : b < v) :
    R20 v = -(1 / (2 * tv v) * (1 - 2 * c (tv v) + c (tv v) ^ 3)) * ((v ^ 2 + 1) / (4 * v ^ 2)) := by
  rw [c_tv hv]
  unfold R20 g20' a tv
  exact idA s_sq hr1 hv hb

lemma R11_eq {v : ℝ} (hv : 0 < v) :
    R11 v = -(1 / 2) * (tv v * (1 - tv v) / W02.q (tv v) ^ 3) * ((v ^ 2 + 1) / (4 * v ^ 2)) := by
  rw [q_tv hv]
  unfold R11 g20' tv
  exact idB s_sq hr1 hv

/-- the folded integrand in `t`. -/
noncomputable def G (t : ℝ) : ℝ := 2 * f t + g t

/-- the derivative of `tv`. -/
noncomputable def tv' (v : ℝ) : ℝ := (v ^ 2 + 1) / (4 * v ^ 2)

lemma F_eq_inner {v : ℝ} (hv : 0 < v) (hb : b < v) :
    l1 v ^ 2 * R20 v + l1 v * l2 v * R11 v = -(1 / 2) * (G (tv v) * tv' v) := by
  rw [R20_eq hv hb, R11_eq hv]
  unfold l1 l2 G f g tv'
  ring

lemma tv_b : tv b = 0 := by
  unfold tv b
  have : (Real.sqrt 2 - 1) ^ 2 + 2 * (Real.sqrt 2 - 1) - 1 = 0 := by nlinarith [s_sq]
  rw [this]; simp

lemma tv_a : tv a = 1 := by
  unfold tv a
  have hs := s_pos
  rw [div_eq_one_iff_eq (by positivity)]
  nlinarith [s_sq]

lemma b_pos : 0 < b := by unfold b; linarith [s_gt]

lemma b_le_a : b ≤ a := by unfold a b; linarith

lemma F_eq {v : ℝ} (hv : v ∈ Icc b a) :
    l1 v ^ 2 * R20 v + l1 v * l2 v * R11 v = -(1 / 2) * (G (tv v) * tv' v) := by
  rcases eq_or_lt_of_le hv.1 with h | h
  · subst h
    simp [l1, l2, tv_b, G, f, g]
  rcases eq_or_lt_of_le hv.2 with h' | h'
  · subst h'
    simp [l1, l2, tv_a, G, f, g]
  exact F_eq_inner (b_pos.trans h) h

lemma tv_mem {v : ℝ} (hv : v ∈ Icc b a) : tv v ∈ Icc (0 : ℝ) 1 := by
  have hs := s_gt
  have h1 := hv.1
  have h2 := hv.2
  unfold b at h1
  unfold a at h2
  have hv0 : 0 < v := by linarith
  unfold tv
  constructor
  · apply div_nonneg _ (by linarith)
    nlinarith [s_sq]
  · rw [div_le_one (by linarith)]
    nlinarith [s_sq]

lemma tv_hasDerivAt {v : ℝ} (hv : v ≠ 0) : HasDerivAt tv (tv' v) v := by
  have h1 : HasDerivAt (fun v : ℝ => v ^ 2 + 2 * v - 1) (2 * v + 2) v := by
    have := (((hasDerivAt_id v).pow 2).add ((hasDerivAt_id v).const_mul 2)).sub_const 1
    simpa using this
  have h2 : HasDerivAt (fun v : ℝ => 4 * v) 4 v := by
    simpa using (hasDerivAt_id v).const_mul 4
  have := h1.div h2 (by positivity)
  unfold tv tv'
  convert this using 1
  field_simp
  ring

lemma G_cont : ContinuousOn G (Icc 0 1) := by
  unfold G
  exact (continuousOn_const.mul f_cont).add g_cont.continuousOn

lemma tv'_cont : ContinuousOn tv' (Icc b a) := by
  intro x hx
  have hx0 : 0 < x := b_pos.trans_le hx.1
  have : ContinuousAt tv' x := by
    unfold tv'
    exact (by fun_prop : Continuous fun x : ℝ => x ^ 2 + 1).continuousAt.div
      (by fun_prop : Continuous fun x : ℝ => 4 * x ^ 2).continuousAt (by positivity)
  exact this.continuousWithinAt

lemma tv_cont : ContinuousOn tv (Icc b a) := fun _ hx =>
  (tv_hasDerivAt (b_pos.trans_le hx.1).ne').continuousAt.continuousWithinAt

lemma inner_cont : ContinuousOn (fun v => -(1 / 2) * (G (tv v) * tv' v)) (Icc b a) :=
  continuousOn_const.mul ((G_cont.comp tv_cont (fun _ hx => tv_mem hx)).mul tv'_cont)

lemma subst : ∫ v in b..a, G (tv v) * tv' v = ∫ t in (0 : ℝ)..1, G t := by
  have h : ∀ x ∈ uIcc b a, HasDerivAt tv (tv' x) x := by
    intro x hx
    rw [uIcc_of_le b_le_a] at hx
    exact tv_hasDerivAt (b_pos.trans_le hx.1).ne'
  have h' : ContinuousOn tv' (uIcc b a) := by rw [uIcc_of_le b_le_a]; exact tv'_cont
  have hg : ContinuousOn G (tv '' uIcc b a) := by
    apply G_cont.mono
    rintro _ ⟨v, hv, rfl⟩
    rw [uIcc_of_le b_le_a] at hv
    exact tv_mem hv
  have := intervalIntegral.integral_comp_mul_deriv' h h' hg
  rw [tv_b, tv_a] at this
  exact this

lemma f_int : IntervalIntegrable f volume 0 1 :=
  (f_cont.mono (by rw [uIcc_of_le zero_le_one])).intervalIntegrable

lemma f1_int : IntervalIntegrable (fun t => f (1 - t)) volume 0 1 := by
  apply ContinuousOn.intervalIntegrable
  rw [uIcc_of_le zero_le_one]
  refine f_cont.comp (continuousOn_const.sub continuousOn_id) ?_
  intro x hx
  exact ⟨by linarith [hx.2], by linarith [hx.1]⟩

lemma g_int : IntervalIntegrable g volume 0 1 := g_cont.intervalIntegrable _ _

lemma phi2_int_eq : ∫ t in (0 : ℝ)..1, phi2 t = ∫ t in (0 : ℝ)..1, G t := by
  have e1 : (fun t => phi2 t) = fun t => f t + f (1 - t) + g t := funext phi2_eq
  have e2 : ∫ t in (0 : ℝ)..1, f (1 - t) = ∫ t in (0 : ℝ)..1, f t := by
    rw [intervalIntegral.integral_comp_sub_left f 1]; norm_num
  rw [e1, intervalIntegral.integral_add (f_int.add f1_int) g_int,
    intervalIntegral.integral_add f_int f1_int, e2]
  unfold G
  rw [intervalIntegral.integral_add (f_int.const_mul 2) g_int, intervalIntegral.integral_const_mul]
  ring

lemma F_int_eq : ∫ v in b..a, (l1 v ^ 2 * R20 v + l1 v * l2 v * R11 v)
    = -(1 / 2) * ∫ t in (0 : ℝ)..1, phi2 t := by
  rw [phi2_int_eq, ← subst, ← intervalIntegral.integral_const_mul]
  apply intervalIntegral.integral_congr
  intro v hv
  rw [uIcc_of_le b_le_a] at hv
  exact F_eq hv

lemma F_intble :
    IntervalIntegrable (fun v => l1 v ^ 2 * R20 v + l1 v * l2 v * R11 v) volume b a := by
  apply ContinuousOn.intervalIntegrable
  rw [uIcc_of_le b_le_a]
  exact inner_cont.congr (fun v hv => F_eq hv)

theorem main : Y01.stmt_5 → Y01.stmt_6 := by
  rintro ⟨_, hM⟩
  refine ⟨F_intble, ?_⟩
  rw [F_int_eq]
  exact hM

end Y01.PB

namespace Y01.PC

open MeasureTheory Real Filter Topology Set

/-! ### Constants -/

lemma s2sq : Real.sqrt 2 ^ 2 = 2 := Real.sq_sqrt (by norm_num)
lemma s2pos : 0 < Real.sqrt 2 := by positivity
lemma s2gt : 1 < Real.sqrt 2 := by
  rw [show (1:ℝ) = Real.sqrt 1 by simp]
  exact Real.sqrt_lt_sqrt (by norm_num) (by norm_num)
lemma s2lt : Real.sqrt 2 < 3 / 2 := by
  rw [Real.sqrt_lt' (by norm_num)]; norm_num

lemma hab : a * b = 1 := by
  unfold a b; nlinarith [s2sq]
lemma hamb : a - b = 2 := by unfold a b; ring
lemma hapb : a + b = 2 * Real.sqrt 2 := by unfold a b; ring
lemma b_pos : 0 < b := by unfold b; linarith [s2gt]
lemma b_lt_one : b < 1 := by unfold b; linarith [s2lt]
lemma one_lt_a : 1 < a := by unfold a; linarith [s2pos]
lemma a_pos : 0 < a := by linarith [one_lt_a]
lemma b_lt_a : b < a := by linarith [b_lt_one, one_lt_a]

/-! ### The substitution `t(v)` -/

lemma tv_eq {v : ℝ} (_hv : v ≠ 0) : tv v = (v - b) * (v + a) / (4 * v) := by
  unfold tv
  congr 1
  linear_combination (-v) * hamb + hab

lemma one_sub_tv {v : ℝ} (hv : v ≠ 0) : 1 - tv v = (a - v) * (v + b) / (4 * v) := by
  unfold tv
  field_simp
  linear_combination (-v) * hamb - hab

lemma tv_pos {v : ℝ} (hv : v ∈ Ioo b a) : 0 < tv v := by
  have hv0 : 0 < v := b_pos.trans hv.1
  rw [tv_eq hv0.ne']
  have := hv.1; have := a_pos
  apply div_pos (mul_pos (by linarith) (by linarith)) (by linarith)

lemma one_sub_tv_pos {v : ℝ} (hv : v ∈ Ioo b a) : 0 < 1 - tv v := by
  have hv0 : 0 < v := b_pos.trans hv.1
  rw [one_sub_tv hv0.ne']
  have := hv.2; have := b_pos
  apply div_pos (mul_pos (by linarith) (by linarith)) (by linarith)

lemma tv_lt_one {v : ℝ} (hv : v ∈ Ioo b a) : tv v < 1 := by
  linarith [one_sub_tv_pos hv]

lemma tv_b : tv b = 0 := by
  rw [tv_eq b_pos.ne']; simp

lemma tv_a : tv a = 1 := by
  have := one_sub_tv a_pos.ne'
  simp at this; linarith

lemma tv_inv {v : ℝ} (hv : v ≠ 0) : tv v⁻¹ = 1 - tv v := by
  unfold tv
  field_simp
  ring

lemma l2_inv {v : ℝ} (hv : v ≠ 0) : l2 v⁻¹ = l1 v := by
  unfold l1 l2; rw [tv_inv hv]; ring_nf

/-! ### Derivatives -/

noncomputable def D1 (v : ℝ) : ℝ := 1 / (v - b) + 1 / (v + a) - 1 / v
noncomputable def D2 (v : ℝ) : ℝ := 1 / (v - a) + 1 / (v + b) - 1 / v
noncomputable def g20 (v : ℝ) : ℝ := a / 4 - Real.sqrt 2 / (2 * (v ^ 2 + 1))

lemma l1_eq {v : ℝ} (hv : v ∈ Ioo b a) :
    l1 v = Real.log (v - b) + Real.log (v + a) - Real.log (4 * v) := by
  have hv0 : 0 < v := b_pos.trans hv.1
  have h1 : 0 < v - b := by linarith [hv.1]
  have h2 : 0 < v + a := by linarith [a_pos]
  unfold l1; rw [tv_eq hv0.ne', Real.log_div (by positivity) (by positivity),
    Real.log_mul h1.ne' h2.ne']

lemma l2_eq {v : ℝ} (hv : v ∈ Ioo b a) :
    l2 v = Real.log (a - v) + Real.log (v + b) - Real.log (4 * v) := by
  have hv0 : 0 < v := b_pos.trans hv.1
  have h1 : 0 < a - v := by linarith [hv.2]
  have h2 : 0 < v + b := by linarith [b_pos]
  unfold l2; rw [one_sub_tv hv0.ne', Real.log_div (by positivity) (by positivity),
    Real.log_mul h1.ne' h2.ne']

lemma hasDerivAt_l1 {v : ℝ} (hv : v ∈ Ioo b a) : HasDerivAt l1 (D1 v) v := by
  have hv0 : 0 < v := b_pos.trans hv.1
  have h1 : 0 < v - b := by linarith [hv.1]
  have h2 : 0 < v + a := by linarith [a_pos]
  have hev : (fun w => Real.log (w - b) + Real.log (w + a) - Real.log (4 * w)) =ᶠ[𝓝 v] l1 :=
    Filter.eventually_of_mem (isOpen_Ioo.mem_nhds hv) (fun w hw => (l1_eq hw).symm)
  have hd : HasDerivAt (fun w => Real.log (w - b) + Real.log (w + a) - Real.log (4 * w))
      (1 / (v - b) + 1 / (v + a) - 4 / (4 * v)) v := by
    have e1 := ((hasDerivAt_id' v).sub_const b).log h1.ne'
    have e2 := ((hasDerivAt_id' v).add_const a).log h2.ne'
    have e3 := ((hasDerivAt_id' v).const_mul 4).log (by positivity)
    convert (e1.add e2).sub e3 using 1; simp
  refine (hd.congr_of_eventuallyEq hev.symm).congr_deriv ?_
  unfold D1; field_simp

lemma hasDerivAt_l2 {v : ℝ} (hv : v ∈ Ioo b a) : HasDerivAt l2 (D2 v) v := by
  have hv0 : 0 < v := b_pos.trans hv.1
  have h1 : 0 < a - v := by linarith [hv.2]
  have h2 : 0 < v + b := by linarith [b_pos]
  have hev : (fun w => Real.log (a - w) + Real.log (w + b) - Real.log (4 * w)) =ᶠ[𝓝 v] l2 :=
    Filter.eventually_of_mem (isOpen_Ioo.mem_nhds hv) (fun w hw => (l2_eq hw).symm)
  have hd : HasDerivAt (fun w => Real.log (a - w) + Real.log (w + b) - Real.log (4 * w))
      (-1 / (a - v) + 1 / (v + b) - 4 / (4 * v)) v := by
    have e1 := ((hasDerivAt_id' v).const_sub a).log h1.ne'
    have e2 := ((hasDerivAt_id' v).add_const b).log h2.ne'
    have e3 := ((hasDerivAt_id' v).const_mul 4).log (by positivity)
    convert (e1.add e2).sub e3 using 1; simp
  refine (hd.congr_of_eventuallyEq hev.symm).congr_deriv ?_
  have h3 : v - a ≠ 0 := by linarith [hv.2]
  unfold D2; field_simp; ring

lemma hasDerivAt_g20 (v : ℝ) : HasDerivAt g20 (g20' v) v := by
  have hp : 0 < 2 * (v ^ 2 + 1) := by positivity
  have e := ((hasDerivAt_pow 2 v).add_const 1).const_mul 2
  have e2 := (e.inv hp.ne').const_mul (Real.sqrt 2)
  have e3 := e2.const_sub (a / 4)
  refine (e3.congr_deriv ?_).congr_of_eventuallyEq ?_
  · unfold g20'; field_simp; ring
  · exact Filter.Eventually.of_forall (fun w => by unfold g20; simp [div_eq_mul_inv])

/-! ### Rational identities -/

noncomputable def r1 (v : ℝ) : ℝ := a * (v + b) / (4 * v * (v + a))
noncomputable def r2 (v : ℝ) : ℝ := a * (v - b) / (4 * v * (v - a))
noncomputable def kc (v : ℝ) : ℝ := b * (v + a) / (4 * v * (v + b))

lemma g20_eq (v : ℝ) : g20 v = a * (v - b) * (v + b) / (4 * (v ^ 2 + 1)) := by
  have hp : 0 < v ^ 2 + 1 := by positivity
  unfold g20
  field_simp
  linear_combination 2 * b * hab + 2 * hapb

lemma pf1 {v : ℝ} (hv : 0 < v) : a * (v + b) / (4 * v * (v + a)) = b / (4 * v) + 1 / (2 * (v + a)) := by
  have h2 : 0 < v + a := by linarith [a_pos]
  field_simp
  linear_combination 2 * v * hamb

lemma pf2 {v : ℝ} (hv : 0 < v) (hva : v ≠ a) :
    a * (v - b) / (4 * v * (v - a)) = b / (4 * v) + 1 / (2 * (v - a)) := by
  have h2 : v - a ≠ 0 := sub_ne_zero.mpr hva
  field_simp
  linear_combination 2 * v * hamb

lemma pf3 {v : ℝ} (hv : 0 < v) :
    b * (v + a) / (4 * v * (v + b)) = a / (4 * v) - 1 / (2 * (v + b)) := by
  have h2 : 0 < v + b := by linarith [b_pos]
  field_simp
  linear_combination (-2 * v) * hamb

lemma c1 {v : ℝ} (hv : v ∈ Ioo b a) : g20 v * D1 v = r1 v := by
  have hv0 : 0 < v := b_pos.trans hv.1
  have h1 : v - b ≠ 0 := by linarith [hv.1]
  have h2 : 0 < v + a := by linarith [a_pos]
  have hp : 0 < v ^ 2 + 1 := by positivity
  rw [g20_eq]; unfold D1 r1
  field_simp
  linear_combination (a * v + a * b) * hab

lemma c2 {v : ℝ} (hv : v ∈ Ioo b a) : g20 v * D2 v = r2 v := by
  have hv0 : 0 < v := b_pos.trans hv.1
  have h1 : v - a ≠ 0 := by linarith [hv.2]
  have h2 : 0 < v + b := by linarith [b_pos]
  have hp : 0 < v ^ 2 + 1 := by positivity
  rw [g20_eq]; unfold D2 r2
  field_simp
  linear_combination (a * v - a * b) * hab

lemma c3 {v : ℝ} (hv : v ∈ Ioo b a) : Rwt v = -2 * r1 v + r2 v + kc v := by
  have hv0 : 0 < v := b_pos.trans hv.1
  have h1 : v - a ≠ 0 := by linarith [hv.2]
  have h2 : 0 < v + b := by linarith [b_pos]
  have h3 : 0 < v + a := by linarith [a_pos]
  unfold r1 r2 kc
  rw [pf1 hv0, pf2 hv0 (by linarith [hv.2]), pf3 hv0]
  unfold Rwt
  field_simp
  linear_combination (-2 * v ^ 3 - 2 * b * v ^ 2 + 2 * a ^ 2 * v + 2 * a ^ 2 * b) * hamb

/-! ### The functions of the integration by parts -/

noncomputable def hfun (v : ℝ) : ℝ := g20' v * (l1 v * (l1 v - l2 v))
noncomputable def k (v : ℝ) : ℝ := r1 v * l2 v
noncomputable def kf (v : ℝ) : ℝ := kc v * l1 v
noncomputable def G (v : ℝ) : ℝ := g20 v * (l1 v * (l1 v - l2 v))
noncomputable def Gd (v : ℝ) : ℝ := hfun v - l1 v * Rwt v - k v + kf v

lemma hasDerivAt_G {v : ℝ} (hv : v ∈ Ioo b a) : HasDerivAt G (Gd v) v := by
  have e := (hasDerivAt_g20 v).mul ((hasDerivAt_l1 hv).mul ((hasDerivAt_l1 hv).sub (hasDerivAt_l2 hv)))
  have hG : G = g20 * (l1 * (l1 - l2)) := rfl
  rw [hG]
  refine e.congr_deriv ?_
  simp only [Pi.mul_apply, Pi.sub_apply]
  unfold Gd hfun k kf
  linear_combination (2 * l1 v - l2 v) * c1 hv - l1 v * c2 hv + l1 v * c3 hv

/-! ### Integrability -/

lemma ii_of_bound {x y : ℝ} (hxy : x ≤ y) {f g : ℝ → ℝ} (hg : IntervalIntegrable g volume x y)
    (hf : AEStronglyMeasurable f volume) (hb : ∀ v ∈ Ioo x y, ‖f v‖ ≤ g v) :
    IntervalIntegrable f volume x y := by
  rw [intervalIntegrable_iff_integrableOn_Ioo_of_le hxy] at hg ⊢
  exact hg.mono' hf.restrict ((ae_restrict_iff' measurableSet_Ioo).mpr (Eventually.of_forall hb))

lemma logsq_le {x y : ℝ} (hy : 0 < y) (hyx : y ≤ x) (hx1 : x ≤ 1) :
    Real.log x ^ 2 ≤ 16 * y ^ (-(1 / 2 : ℝ)) := by
  have hx0 : 0 < x := hy.trans_le hyx
  have h := neg_rpow_div_le_log hx0.le (ε := 1 / 4) (by norm_num)
  have hl : Real.log x ≤ 0 := Real.log_nonpos hx0.le hx1
  have hp : 0 ≤ x ^ (-(1 / 4 : ℝ)) := Real.rpow_nonneg hx0.le _
  have hsq : (x ^ (-(1 / 4 : ℝ))) ^ 2 = x ^ (-(1 / 2 : ℝ)) := by
    rw [← Real.rpow_natCast, ← Real.rpow_mul hx0.le]; norm_num
  have hmono : x ^ (-(1 / 2 : ℝ)) ≤ y ^ (-(1 / 2 : ℝ)) :=
    Real.rpow_le_rpow_of_nonpos hy hyx (by norm_num)
  have h4 : -Real.log x ≤ 4 * x ^ (-(1 / 4 : ℝ)) := by
    have : -x ^ (-(1 / 4 : ℝ)) / (1 / 4) = -(4 * x ^ (-(1 / 4 : ℝ))) := by ring
    linarith
  have : Real.log x ^ 2 ≤ (4 * x ^ (-(1 / 4 : ℝ))) ^ 2 := by nlinarith
  nlinarith

lemma tv_ge {v : ℝ} (hv : v ∈ Ioo b a) : (v - b) / 4 ≤ tv v := by
  have hv0 : 0 < v := b_pos.trans hv.1
  rw [tv_eq hv0.ne', div_le_div_iff₀ (by norm_num) (by positivity)]
  have := hv.1; have := a_pos
  nlinarith

lemma one_sub_tv_ge {v : ℝ} (hv : v ∈ Ioo b a) : (a - v) / 4 ≤ 1 - tv v := by
  have hv0 : 0 < v := b_pos.trans hv.1
  rw [one_sub_tv hv0.ne', div_le_div_iff₀ (by norm_num) (by positivity)]
  have := hv.2; have := b_pos
  nlinarith

noncomputable def K : ℝ := 16 / (4 : ℝ) ^ (-(1 / 2 : ℝ))

lemma K_pos : 0 < K := by unfold K; positivity

lemma l1sq_le {v : ℝ} (hv : v ∈ Ioo b a) : l1 v ^ 2 ≤ K * (v - b) ^ (-(1 / 2 : ℝ)) := by
  have h1 : 0 < v - b := by linarith [hv.1]
  have := logsq_le (by positivity) (tv_ge hv) (tv_lt_one hv).le
  rw [Real.div_rpow h1.le (by norm_num)] at this
  unfold l1 K; calc _ ≤ _ := this
    _ = _ := by ring

lemma l2sq_le {v : ℝ} (hv : v ∈ Ioo b a) : l2 v ^ 2 ≤ K * (a - v) ^ (-(1 / 2 : ℝ)) := by
  have h1 : 0 < a - v := by linarith [hv.2]
  have := logsq_le (by positivity) (one_sub_tv_ge hv) (by linarith [tv_pos hv])
  rw [Real.div_rpow h1.le (by norm_num)] at this
  unfold l2 K; calc _ ≤ _ := this
    _ = _ := by ring

noncomputable def Phi (v : ℝ) : ℝ :=
  K * (v - b) ^ (-(1 / 2 : ℝ)) + K * (a - v) ^ (-(1 / 2 : ℝ)) + 1

lemma ii_Phi : IntervalIntegrable Phi volume b a := by
  have h1 := (intervalIntegral.intervalIntegrable_rpow' (r := -(1 / 2 : ℝ)) (by norm_num)
    (a := 0) (b := a - b)).comp_sub_right b
  have h2 := (intervalIntegral.intervalIntegrable_rpow' (r := -(1 / 2 : ℝ)) (by norm_num)
    (a := a - b) (b := 0)).comp_sub_left a
  simp only [zero_add, sub_add_cancel, sub_sub_cancel, sub_zero] at h1 h2
  exact ((h1.const_mul K).add (h2.const_mul K)).add intervalIntegrable_const

lemma Phi_ge {v : ℝ} (hv : v ∈ Ioo b a) :
    l1 v ^ 2 + 1 ≤ Phi v ∧ l2 v ^ 2 + 1 ≤ Phi v ∧ l1 v ^ 2 + l2 v ^ 2 + 1 ≤ Phi v := by
  have h1 : 0 < v - b := by linarith [hv.1]
  have h2 : 0 < a - v := by linarith [hv.2]
  have p1 : 0 ≤ K * (v - b) ^ (-(1 / 2 : ℝ)) := by have := K_pos; positivity
  have p2 : 0 ≤ K * (a - v) ^ (-(1 / 2 : ℝ)) := by have := K_pos; positivity
  have := l1sq_le hv; have := l2sq_le hv
  unfold Phi
  refine ⟨by linarith, by linarith, by linarith⟩

lemma meas_l1 : Measurable l1 := by unfold l1 tv; fun_prop
lemma meas_l2 : Measurable l2 := by unfold l2 tv; fun_prop

lemma ii_l1sq : IntervalIntegrable (fun v => l1 v ^ 2) volume b a :=
  ii_of_bound b_lt_a.le ii_Phi (meas_l1.pow_const 2).aestronglyMeasurable (fun v hv => by
    have := Phi_ge hv; rw [Real.norm_eq_abs, abs_of_nonneg (sq_nonneg _)]; linarith [this.1])

lemma ii_l1l2 : IntervalIntegrable (fun v => l1 v * l2 v) volume b a :=
  ii_of_bound b_lt_a.le ii_Phi (meas_l1.mul meas_l2).aestronglyMeasurable (fun v hv => by
    have := Phi_ge hv; rw [Real.norm_eq_abs]
    apply abs_le.mpr
    constructor <;> nlinarith [sq_nonneg (l1 v + l2 v), sq_nonneg (l1 v - l2 v)])

lemma ii_l1 : IntervalIntegrable l1 volume b a :=
  ii_of_bound b_lt_a.le ii_Phi meas_l1.aestronglyMeasurable (fun v hv => by
    have := Phi_ge hv; rw [Real.norm_eq_abs]
    apply abs_le.mpr
    constructor <;> nlinarith [sq_nonneg (l1 v + 1), sq_nonneg (l1 v - 1)])

lemma ii_l2 : IntervalIntegrable l2 volume b a :=
  ii_of_bound b_lt_a.le ii_Phi meas_l2.aestronglyMeasurable (fun v hv => by
    have := Phi_ge hv; rw [Real.norm_eq_abs]
    apply abs_le.mpr
    constructor <;> nlinarith [sq_nonneg (l2 v + 1), sq_nonneg (l2 v - 1)])

lemma cont_g20' : Continuous g20' := by
  rw [show g20' = fun v => Real.sqrt 2 * v / (v ^ 2 + 1) ^ 2 from rfl]
  exact Continuous.div (f := fun v => Real.sqrt 2 * v) (g := fun v => (v ^ 2 + 1) ^ 2)
    (by fun_prop) (by fun_prop) (fun x => by positivity)

lemma ii_hfun : IntervalIntegrable hfun volume b a := by
  have e : hfun = fun v => g20' v * (l1 v ^ 2 - l1 v * l2 v) := by
    funext v; unfold hfun; ring
  rw [e]
  exact (ii_l1sq.sub ii_l1l2).continuousOn_mul cont_g20'.continuousOn

lemma pos_of_mem {v : ℝ} (hv : v ∈ uIcc b a) : 0 < v := by
  rw [uIcc_of_le b_lt_a.le] at hv; exact b_pos.trans_le hv.1

lemma ii_k : IntervalIntegrable k volume b a := by
  unfold k
  refine ii_l2.continuousOn_mul ?_
  intro x hx
  have := pos_of_mem hx; have := a_pos; have := b_pos
  have : x + a ≠ 0 := by positivity
  apply ContinuousAt.continuousWithinAt
  unfold r1
  fun_prop (disch := positivity)

lemma ii_kf : IntervalIntegrable kf volume b a := by
  unfold kf
  refine ii_l1.continuousOn_mul ?_
  intro x hx
  have := pos_of_mem hx; have := a_pos; have := b_pos
  have : x + b ≠ 0 := by positivity
  apply ContinuousAt.continuousWithinAt
  unfold kc
  fun_prop (disch := positivity)

noncomputable def Rr (v : ℝ) : ℝ := 1 / (2 * v) - 1 / (v + a) - 1 / (2 * (v + b))

lemma Rwt_split (v : ℝ) : Rwt v = Rr v + 1 / (2 * (v - a)) := by
  unfold Rwt Rr; ring

lemma tv_ge_half {v : ℝ} (hv : 1 ≤ v) : 1 / 2 ≤ tv v := by
  unfold tv
  rw [div_le_div_iff₀ (by norm_num) (by linarith)]
  nlinarith

lemma l1_abs_le {v : ℝ} (hv : v ∈ Ioo 1 a) : |l1 v| ≤ 2 * (1 - tv v) := by
  have hv' : v ∈ Ioo b a := ⟨b_lt_one.trans hv.1, hv.2⟩
  have ht0 := tv_pos hv'
  have ht1 := tv_lt_one hv'
  have hth := tv_ge_half hv.1.le
  have hl : l1 v ≤ 0 := Real.log_nonpos ht0.le ht1.le
  have hlog : -l1 v ≤ 1 / tv v - 1 := by
    have := Real.log_le_sub_one_of_pos (one_div_pos.mpr ht0)
    rw [one_div, Real.log_inv] at this
    unfold l1; rw [one_div]; linarith
  have h2 : 1 / tv v - 1 ≤ 2 * (1 - tv v) := by
    rw [div_sub_one ht0.ne', div_le_iff₀ ht0]
    nlinarith
  rw [abs_of_nonpos hl]; linarith

lemma l1_pole_bound {v : ℝ} (hv : v ∈ Ioo 1 a) : ‖l1 v * (1 / (2 * (v - a)))‖ ≤ 1 := by
  have hv' : v ∈ Ioo b a := ⟨b_lt_one.trans hv.1, hv.2⟩
  have hv0 : 0 < v := by linarith [hv.1]
  have ht0 := tv_pos hv'
  have ht1 := tv_lt_one hv'
  have hth := tv_ge_half hv.1.le
  have hl : l1 v ≤ 0 := Real.log_nonpos ht0.le ht1.le
  have hlog : -l1 v ≤ 1 / tv v - 1 := by
    have := Real.log_le_sub_one_of_pos (one_div_pos.mpr ht0)
    rw [one_div, Real.log_inv] at this
    unfold l1; rw [one_div]; linarith
  have h2 : 1 / tv v - 1 ≤ 2 * (1 - tv v) := by
    rw [div_sub_one ht0.ne', div_le_iff₀ ht0]
    nlinarith
  have h3 : 2 * (1 - tv v) ≤ 2 * (a - v) := by
    rw [one_sub_tv hv0.ne']
    have : (a - v) * (v + b) / (4 * v) ≤ a - v := by
      rw [div_le_iff₀ (by positivity)]
      nlinarith [mul_le_mul_of_nonneg_left (show v + b ≤ 4 * v by linarith [b_lt_one, hv.1])
        (show 0 ≤ a - v by linarith [hv.2])]
    linarith
  have hneg : 2 * (v - a) < 0 := by linarith [hv.2]
  rw [Real.norm_eq_abs, abs_mul, abs_of_nonpos hl, abs_div, abs_one, abs_of_neg hneg,
    mul_one_div, div_le_one (by linarith)]
  linarith

lemma ii_l1Rwt : IntervalIntegrable (fun v => l1 v * Rwt v) volume b a := by
  have hsub1 : uIcc b 1 ⊆ uIcc b a := by
    rw [uIcc_of_le b_lt_one.le, uIcc_of_le b_lt_a.le]
    exact Icc_subset_Icc le_rfl one_lt_a.le
  have hsub2 : uIcc 1 a ⊆ uIcc b a := by
    rw [uIcc_of_le one_lt_a.le, uIcc_of_le b_lt_a.le]
    exact Icc_subset_Icc b_lt_one.le le_rfl
  have i1 : IntervalIntegrable (fun v => l1 v * Rwt v) volume b 1 := by
    refine (ii_l1.mono_set hsub1).mul_continuousOn ?_
    rw [uIcc_of_le b_lt_one.le]
    intro x hx
    have hx0 : 0 < x := b_pos.trans_le hx.1
    have hxa : x - a ≠ 0 := by linarith [hx.2, one_lt_a]
    have := a_pos; have := b_pos
    unfold Rwt
    apply ContinuousAt.continuousWithinAt
    have : x + a ≠ 0 := by positivity
    have : x + b ≠ 0 := by positivity
    fun_prop (disch := positivity)
  have i2 : IntervalIntegrable (fun v => l1 v * Rr v) volume 1 a := by
    refine (ii_l1.mono_set hsub2).mul_continuousOn ?_
    rw [uIcc_of_le one_lt_a.le]
    intro x hx
    have hx0 : 0 < x := by linarith [hx.1]
    have := a_pos; have := b_pos
    unfold Rr
    apply ContinuousAt.continuousWithinAt
    have : x + a ≠ 0 := by positivity
    have : x + b ≠ 0 := by positivity
    fun_prop (disch := positivity)
  have i3 : IntervalIntegrable (fun v => l1 v * (1 / (2 * (v - a)))) volume 1 a :=
    ii_of_bound one_lt_a.le intervalIntegrable_const
      (meas_l1.mul (by fun_prop)).aestronglyMeasurable (fun v hv => l1_pole_bound hv)
  have i4 : IntervalIntegrable (fun v => l1 v * Rwt v) volume 1 a := by
    have e : (fun v => l1 v * Rwt v) = fun v => l1 v * Rr v + l1 v * (1 / (2 * (v - a))) := by
      funext v; rw [Rwt_split]; ring
    rw [e]; exact i2.add i3
  exact i1.trans i4

lemma ii_Gd : IntervalIntegrable Gd volume b a := by
  have e : Gd = fun v => hfun v - l1 v * Rwt v - k v + kf v := rfl
  rw [e]
  exact ((ii_hfun.sub ii_l1Rwt).sub ii_k).add ii_kf

/-! ### Continuity of `G` on `[b,a]` and the fundamental theorem -/

lemma g20_b : g20 b = 0 := by rw [g20_eq]; simp

lemma l1_a : l1 a = 0 := by unfold l1; rw [tv_a]; simp

lemma G_b : G b = 0 := by unfold G; rw [g20_b]; ring

lemma G_a : G a = 0 := by unfold G; rw [l1_a]; ring

lemma contAt_tv {v : ℝ} (hv : v ≠ 0) : ContinuousAt tv v := by
  unfold tv
  have : 4 * v ≠ 0 := by positivity
  fun_prop (disch := assumption)

noncomputable def p (v : ℝ) : ℝ := a * (v + b) * v / ((v ^ 2 + 1) * (v + a))

lemma g20_eq_p {v : ℝ} (hv : 0 < v) : g20 v = p v * tv v := by
  have h2 : 0 < v + a := by linarith [a_pos]
  have hp : 0 < v ^ 2 + 1 := by positivity
  rw [g20_eq, tv_eq hv.ne']; unfold p
  field_simp

lemma contAt_l2_b : ContinuousAt l2 b := by
  unfold l2
  apply ContinuousAt.log ((continuousAt_const).sub (contAt_tv b_pos.ne'))
  show (1:ℝ) - tv b ≠ 0
  rw [tv_b]; norm_num

lemma l2_b : l2 b = 0 := by unfold l2; rw [tv_b]; simp

lemma tendsto_G_b : Tendsto G (𝓝[>] b) (𝓝 (G b)) := by
  have hT : Tendsto tv (𝓝[>] b) (𝓝[>] 0) := by
    rw [tendsto_nhdsWithin_iff]
    refine ⟨?_, ?_⟩
    · have := (contAt_tv b_pos.ne').tendsto
      rw [tv_b] at this
      exact this.mono_left nhdsWithin_le_nhds
    · filter_upwards [Ioo_mem_nhdsGT b_lt_a] with v hv using tv_pos hv
  have hTc : Tendsto tv (𝓝[>] b) (𝓝 0) := tendsto_nhdsWithin_iff.1 hT |>.1
  have h1 := (tendsto_log_mul_rpow_nhdsGT_zero (r := 1 / 2) (by norm_num)).comp hT
  have h2 : Tendsto (fun v => tv v * Real.log (tv v)) (𝓝[>] b) (𝓝 0) := by
    have := (Real.continuous_mul_log.tendsto 0).comp hTc
    simp only [Function.comp_def, zero_mul] at this
    exact this
  have hpc : ContinuousAt p b := by
    unfold p
    have := a_pos; have := b_pos
    have : (b ^ 2 + 1) * (b + a) ≠ 0 := by positivity
    fun_prop (disch := assumption)
  have hp := hpc.tendsto.mono_left (nhdsWithin_le_nhds (s := Ioi b))
  have hl2 := contAt_l2_b.tendsto.mono_left (nhdsWithin_le_nhds (s := Ioi b))
  rw [l2_b] at hl2
  have hlim := (hp.mul (h1.pow 2)).sub ((hp.mul h2).mul hl2)
  simp only [Function.comp_def] at hlim
  rw [G_b]
  have e0 : p b * (0:ℝ) ^ 2 - p b * 0 * 0 = 0 := by ring
  rw [e0] at hlim
  refine hlim.congr' ?_
  filter_upwards [Ioo_mem_nhdsGT b_lt_a] with v hv
  have hv0 : 0 < v := b_pos.trans hv.1
  have ht := tv_pos hv
  have hsq : (tv v ^ (1 / 2 : ℝ)) ^ 2 = tv v := by
    rw [← Real.rpow_natCast, ← Real.rpow_mul ht.le]; norm_num
  unfold G; rw [g20_eq_p hv0]; unfold l1
  linear_combination (p v * Real.log (tv v) ^ 2) * hsq

lemma contAt_l1_a : ContinuousAt l1 a := by
  unfold l1
  apply ContinuousAt.log (contAt_tv a_pos.ne')
  rw [tv_a]; norm_num

lemma tendsto_G_a : Tendsto G (𝓝[<] a) (𝓝 (G a)) := by
  have hg := (hasDerivAt_g20 a).continuousAt.tendsto.mono_left (nhdsWithin_le_nhds (s := Iio a))
  have hl1 := contAt_l1_a.tendsto.mono_left (nhdsWithin_le_nhds (s := Iio a))
  rw [l1_a] at hl1
  have hS : Tendsto (fun v => (1 - tv v) * Real.log (1 - tv v)) (𝓝[<] a) (𝓝 0) := by
    have hc : Tendsto (fun v => 1 - tv v) (𝓝 a) (𝓝 0) := by
      have h : Tendsto (fun v => 1 - tv v) (𝓝 a) (𝓝 (1 - tv a)) :=
        tendsto_const_nhds.sub (contAt_tv a_pos.ne').tendsto
      rwa [tv_a, sub_self] at h
    have := (Real.continuous_mul_log.tendsto 0).comp hc
    simp only [Function.comp_def, zero_mul] at this
    exact this.mono_left nhdsWithin_le_nhds
  have hprod : Tendsto (fun v => l1 v * l2 v) (𝓝[<] a) (𝓝 0) := by
    refine squeeze_zero_norm' ?_ (by simpa using (hS.norm.const_mul 2))
    filter_upwards [Ioo_mem_nhdsLT one_lt_a] with v hv
    have hv' : v ∈ Ioo b a := ⟨b_lt_one.trans hv.1, hv.2⟩
    have h1 := l1_abs_le hv
    have h0 : 0 ≤ 1 - tv v := (one_sub_tv_pos hv').le
    simp only [Real.norm_eq_abs, abs_mul]
    rw [abs_of_nonneg h0]
    unfold l2
    nlinarith [abs_nonneg (Real.log (1 - tv v))]
  have hlim := (hg.mul (hl1.pow 2)).sub (hg.mul hprod)
  rw [G_a]
  have e0 : g20 a * (0:ℝ) ^ 2 - g20 a * 0 = 0 := by ring
  rw [e0] at hlim
  refine hlim.congr' (Eventually.of_forall fun v => ?_)
  unfold G; ring

lemma contOn_G : ContinuousOn G (Icc b a) := by
  intro x hx
  rcases eq_or_lt_of_le hx.1 with h | h
  · subst h
    exact (continuousWithinAt_Icc_iff_Ici b_lt_a).2
      (continuousWithinAt_Ioi_iff_Ici.1 tendsto_G_b)
  rcases eq_or_lt_of_le hx.2 with h' | h'
  · subst h'
    exact (continuousWithinAt_Icc_iff_Iic b_lt_a).2
      (continuousWithinAt_Iio_iff_Iic.1 tendsto_G_a)
  exact (hasDerivAt_G ⟨h, h'⟩).continuousAt.continuousWithinAt

lemma int_Gd : ∫ v in b..a, Gd v = 0 := by
  rw [intervalIntegral.integral_eq_sub_of_hasDerivAt_of_le b_lt_a.le contOn_G
    (fun x hx => hasDerivAt_G hx) ii_Gd, G_a, G_b, sub_zero]

/-! ### Folding `v ↦ 1/v` -/

lemma b_eq_inv : b = a⁻¹ := eq_inv_of_mul_eq_one_right hab
lemma a_eq_inv : a = b⁻¹ := eq_inv_of_mul_eq_one_left hab

lemma image_inv : (fun u : ℝ => u⁻¹) '' Ioo b a = Ioo b a := by
  ext x
  constructor
  · rintro ⟨u, hu, rfl⟩
    have hu0 : 0 < u := b_pos.trans hu.1
    refine ⟨?_, ?_⟩
    · rw [b_eq_inv]; exact (inv_lt_inv₀ a_pos hu0).2 hu.2
    · conv_rhs => rw [a_eq_inv]
      exact (inv_lt_inv₀ hu0 b_pos).2 hu.1
  · intro hx
    have hx0 : 0 < x := b_pos.trans hx.1
    refine ⟨x⁻¹, ⟨?_, ?_⟩, inv_inv x⟩
    · rw [b_eq_inv]; exact (inv_lt_inv₀ a_pos hx0).2 hx.2
    · conv_rhs => rw [a_eq_inv]
      exact (inv_lt_inv₀ hx0 b_pos).2 hx.1

lemma fold_pt {u : ℝ} (hu : u ∈ Ioo b a) : |-(u ^ 2)⁻¹| • k u⁻¹ = kf u := by
  have hu0 : 0 < u := b_pos.trans hu.1
  have h2 : 0 < u + b := by linarith [b_pos]
  have h3 : 0 < u + a := by linarith [a_pos]
  have := a_pos; have := b_pos
  unfold k kf
  rw [l2_inv hu0.ne', abs_neg, abs_of_pos (by positivity), smul_eq_mul]
  unfold r1 kc
  have h4 : 0 < u⁻¹ + a := by positivity
  field_simp
  linear_combination (u * l1 u * (b - a)) * hab

lemma fold_k : ∫ v in b..a, k v = ∫ v in b..a, kf v := by
  have himg := integral_image_eq_integral_abs_deriv_smul (s := Ioo b a)
    (f := fun u : ℝ => u⁻¹) (f' := fun u => -(u ^ 2)⁻¹) measurableSet_Ioo
    (fun x hx => (hasDerivAt_inv (b_pos.trans hx.1).ne').hasDerivWithinAt)
    inv_injective.injOn k
  rw [image_inv] at himg
  rw [intervalIntegral.integral_of_le b_lt_a.le, intervalIntegral.integral_of_le b_lt_a.le,
    integral_Ioc_eq_integral_Ioo, integral_Ioc_eq_integral_Ioo, himg]
  exact setIntegral_congr_fun measurableSet_Ioo (fun u hu => fold_pt hu)

/-! ### Assembly -/

lemma int_hfun : ∫ v in b..a, hfun v = Tint := by
  have h := int_Gd
  have e : ∫ v in b..a, Gd v = (∫ v in b..a, hfun v) - (∫ v in b..a, l1 v * Rwt v)
      - (∫ v in b..a, k v) + ∫ v in b..a, kf v := by
    show ∫ v in b..a, (hfun v - l1 v * Rwt v - k v + kf v) = _
    rw [intervalIntegral.integral_add ((ii_hfun.sub ii_l1Rwt).sub ii_k) ii_kf,
      intervalIntegral.integral_sub (ii_hfun.sub ii_l1Rwt) ii_k,
      intervalIntegral.integral_sub ii_hfun ii_l1Rwt]
  have hT : Tint = ∫ v in b..a, l1 v * Rwt v := rfl
  linarith [fold_k]

/-- Lemma g1:lem:red2 from Lemma g1:lem:k2. -/
theorem main : stmt_6 → stmt_7 := by
  rintro ⟨hS, hM⟩
  have hpt : ∀ v, l1 v ^ 2 * R20 v + l1 v * l2 v * R11 v = W3integrand v + hfun v := by
    intro v; unfold R20 R11 W3integrand hfun; ring
  have hW : IntervalIntegrable W3integrand volume b a := by
    have e : W3integrand = fun v => (l1 v ^ 2 * R20 v + l1 v * l2 v * R11 v) - hfun v := by
      funext v; rw [hpt]; ring
    rw [e]; exact hS.sub ii_hfun
  refine ⟨hW, ii_l1Rwt, ?_⟩
  rw [hM]
  congr 1
  simp_rw [hpt]
  rw [intervalIntegral.integral_add hW ii_hfun, int_hfun]
  rfl

end Y01.PC

namespace Y01.PD

open MeasureTheory Real Filter Topology Set

lemma s_sq : Real.sqrt 2 ^ 2 = 2 := Real.sq_sqrt (by norm_num)
lemma s_pos : 0 < Real.sqrt 2 := by positivity
lemma s_gt_one : 1 < Real.sqrt 2 := by
  rw [show (1:ℝ) = Real.sqrt 1 by simp]
  exact Real.sqrt_lt_sqrt (by norm_num) (by norm_num)
lemma s_lt : Real.sqrt 2 < 3 / 2 := by nlinarith [s_sq, s_pos]

lemma b_pos : 0 < b := by unfold b; linarith [s_gt_one]
lemma a_pos : 0 < a := by unfold a; linarith [s_gt_one]
lemma b_lt_a : b < a := by unfold a b; linarith
lemma ab : a * b = 1 := by unfold a b; nlinarith [s_sq]

noncomputable def rho (v : ℝ) : ℝ := (a - v) * (v - b) / (4 * Real.sqrt 2 * v)
noncomputable def w (v : ℝ) : ℝ := (a - v) / (a + v)
noncomputable def mm (v : ℝ) : ℝ := Real.log (1 - rho v)

lemma tv_eq {v : ℝ} (hv : v ≠ 0) : tv v = (v - b) * (v + a) / (4 * v) := by
  unfold tv a b; field_simp; ring_nf; rw [s_sq]; ring

lemma one_sub_tv {v : ℝ} (hv : v ≠ 0) : 1 - tv v = (a - v) * (v + b) / (4 * v) := by
  unfold tv a b; field_simp; ring_nf; rw [s_sq]; ring

lemma rho_eq {v : ℝ} (hv : v ≠ 0) (hav : a + v ≠ 0) :
    rho v = tv v * (w v / Real.sqrt 2) := by
  rw [tv_eq hv]; unfold rho w
  have := s_pos.ne'
  field_simp; ring

lemma s3 : Real.sqrt 2 ^ 3 = 2 * Real.sqrt 2 := by rw [pow_succ, s_sq]
lemma s4 : Real.sqrt 2 ^ 4 = 4 := by rw [show (4:ℕ) = 2 * 2 from rfl, pow_mul, s_sq]; norm_num

lemma mem_facts {v : ℝ} (hv : v ∈ Ioo b a) :
    0 < v ∧ 0 < a - v ∧ 0 < a + v ∧ 0 < v - b ∧ 0 < v + b := by
  obtain ⟨h1, h2⟩ := hv
  have := b_pos; have := a_pos
  refine ⟨by linarith, by linarith, by linarith, by linarith, by linarith⟩

lemma one_sub_rho {v : ℝ} (hv : v ∈ Ioo b a) :
    1 - rho v = (1 - w v / Real.sqrt 2) / ((1 - w v) * (1 + w v)) := by
  obtain ⟨h0, h1, h2, h3, h4⟩ := mem_facts hv
  have hs := s_pos.ne'
  have e1 : 1 - w v = 2 * v / (a + v) := by unfold w; field_simp; ring
  have e2 : 1 + w v = 2 * a / (a + v) := by unfold w; field_simp; ring
  rw [e1, e2]
  have ha := a_pos
  unfold rho w
  field_simp
  unfold a b
  ring_nf

lemma w_pos {v : ℝ} (hv : v ∈ Ioo b a) : 0 < w v ∧ w v < 1 := by
  obtain ⟨h0, h1, h2, h3, h4⟩ := mem_facts hv
  unfold w
  exact ⟨div_pos h1 h2, (div_lt_one h2).2 (by linarith)⟩

lemma tv_mem {v : ℝ} (hv : v ∈ Ioo b a) : 0 < tv v ∧ tv v < 1 := by
  obtain ⟨h0, h1, h2, h3, h4⟩ := mem_facts hv
  have e1 := tv_eq h0.ne'
  have e2 := one_sub_tv h0.ne'
  constructor
  · rw [e1]; exact div_pos (mul_pos h3 (by linarith)) (by linarith)
  · have : 0 < 1 - tv v := by rw [e2]; exact div_pos (mul_pos h1 h4) (by linarith)
    linarith

lemma mm_eq {v : ℝ} (hv : v ∈ Ioo b a) :
    mm v = Real.log (1 - w v / Real.sqrt 2) - Real.log (1 - w v) - Real.log (1 + w v) := by
  obtain ⟨h0, h1, h2, h3, h4⟩ := mem_facts hv
  obtain ⟨w0, w1⟩ := w_pos hv
  have hs := s_gt_one
  have p1 : 0 < 1 - w v / Real.sqrt 2 := by
    rw [sub_pos, div_lt_one (by linarith)]; linarith
  unfold mm
  rw [one_sub_rho hv, Real.log_div p1.ne'
    (by nlinarith), Real.log_mul (by linarith) (by linarith)]
  ring

lemma d_tv {v : ℝ} (hv : v ≠ 0) : HasDerivAt tv ((v ^ 2 + 1) / (4 * v ^ 2)) v := by
  have h := ((((hasDerivAt_id v).pow 2).add ((hasDerivAt_id v).const_mul 2)).sub_const 1).div
    ((hasDerivAt_id v).const_mul 4) (by simpa using hv)
  convert h using 1
  · funext x; simp [tv]
  · simp; field_simp; ring

lemma d_rho {v : ℝ} (hv : v ≠ 0) :
    HasDerivAt rho ((1 - v ^ 2) / (4 * Real.sqrt 2 * v ^ 2)) v := by
  have hs := s_pos.ne'
  have h := (((hasDerivAt_id v).const_sub a).mul ((hasDerivAt_id v).sub_const b)).div
    ((hasDerivAt_id v).const_mul (4 * Real.sqrt 2)) (mul_ne_zero (mul_ne_zero (by norm_num) hs) hv)
  convert h using 1
  · funext x; simp [rho]
  · simp; field_simp; unfold a b; ring_nf; simp only [s_sq]; ring

lemma d_w {v : ℝ} (hav : a + v ≠ 0) : HasDerivAt w (-2 * a / (a + v) ^ 2) v := by
  have h := ((hasDerivAt_id v).const_sub a).div ((hasDerivAt_id v).const_add a) (by simpa using hav)
  convert h using 1
  · funext x; simp [w]
  · simp; field_simp; ring

lemma one_sub_rho' {v : ℝ} (hv : v ≠ 0) : 1 - rho v = (v + a) * (v + b) / (4 * Real.sqrt 2 * v) := by
  have hs := s_pos.ne'
  unfold rho; field_simp; unfold a b; ring_nf

lemma hX {v : ℝ} (hv : v ∈ Ioo b a) :
    -(1 / 2) * (((v ^ 2 + 1) / (4 * v ^ 2)) / (1 - tv v))
      + ((1 - v ^ 2) / (4 * Real.sqrt 2 * v ^ 2)) / (1 - rho v) = Rwt v := by
  obtain ⟨h0, h1, h2, h3, h4⟩ := mem_facts hv
  have hs := s_pos.ne'
  rw [one_sub_tv h0.ne', one_sub_rho' h0.ne']
  have h5 : v - a ≠ 0 := by linarith
  have h6 : v + a ≠ 0 := by linarith
  unfold Rwt
  field_simp
  unfold a b; ring_nf; simp only [s_sq, s3, s4]; ring_nf

lemma hY {v : ℝ} (hv : v ∈ Ioo b a) :
    -(((v ^ 2 + 1) / (4 * v ^ 2)) / tv v) + ((1 - v ^ 2) / (4 * Real.sqrt 2 * v ^ 2)) / rho v
      - (-2 * a / (a + v) ^ 2) / w v = 0 := by
  obtain ⟨h0, h1, h2, h3, h4⟩ := mem_facts hv
  have hs := s_pos.ne'
  rw [tv_eq h0.ne']
  have h5 : a - v ≠ 0 := by linarith
  have h6 : v + a ≠ 0 := by linarith
  have h7 : v - b ≠ 0 := by linarith
  unfold rho w
  field_simp
  unfold a b; ring_nf; simp only [s_sq, s3]; ring_nf

/-- The primitive of `l1 * Rwt` on `(b,a)`. -/
noncomputable def F (v : ℝ) : ℝ :=
  -(1 / 2) * D00.Li2 (1 - tv v) - l1 v * mm v - D00.Li2 (rho v)
    - (-D00.Li2 (w v / Real.sqrt 2) + D00.Li2 (w v) + D00.Li2 (-w v))

lemma rho_mem {v : ℝ} (hv : v ∈ Ioo b a) : 0 < rho v ∧ rho v < 1 := by
  obtain ⟨h0, h1, h2, h3, h4⟩ := mem_facts hv
  obtain ⟨t0, t1⟩ := tv_mem hv
  obtain ⟨w0, w1⟩ := w_pos hv
  have hs := s_gt_one
  have q0 : 0 < w v / Real.sqrt 2 := div_pos w0 (by linarith)
  have q1 : w v / Real.sqrt 2 < 1 := by rw [div_lt_one (by linarith)]; linarith
  rw [rho_eq h0.ne' h2.ne']
  exact ⟨mul_pos t0 q0, by nlinarith⟩

lemma hasDerivAt_F {v : ℝ} (hv : v ∈ Ioo b a) : HasDerivAt F (l1 v * Rwt v) v := by
  obtain ⟨h0, h1, h2, h3, h4⟩ := mem_facts hv
  obtain ⟨t0, t1⟩ := tv_mem hv
  obtain ⟨w0, w1⟩ := w_pos hv
  obtain ⟨r0, r1⟩ := rho_mem hv
  have hs := s_gt_one
  have hs0 : Real.sqrt 2 ≠ 0 := by linarith
  have q0 : 0 < w v / Real.sqrt 2 := div_pos w0 (by linarith)
  have q1 : w v / Real.sqrt 2 < 1 := by rw [div_lt_one (by linarith)]; linarith
  have ht := d_tv h0.ne'
  have hr := d_rho h0.ne'
  have hw := d_w h2.ne'
  have u1 : 1 - tv v < 1 := by linarith
  have u0 : 1 - tv v ≠ 0 := by linarith
  have P1 := HasDerivAt.comp (h := fun x => 1 - tv x) v (D00.hasDerivAt_Li2 u1 u0) (ht.const_sub 1)
  have P2 := ht.log t0.ne'
  have r2 : 1 - rho v ≠ 0 := by linarith
  have P2m := (hr.const_sub 1).log r2
  have P3 := (D00.hasDerivAt_Li2 r1 r0.ne').comp v hr
  have P4 := HasDerivAt.comp (h := fun x => w x / Real.sqrt 2) v (D00.hasDerivAt_Li2 q1 q0.ne')
    (hw.div_const (Real.sqrt 2))
  have P5 := (D00.hasDerivAt_Li2 w1 w0.ne').comp v hw
  have n1 : -w v < 1 := by linarith
  have n0 : -w v ≠ 0 := by linarith
  have P6 := HasDerivAt.comp (h := fun x => -w x) v (D00.hasDerivAt_Li2 n1 n0) hw.neg
  have H := (((P1.const_mul (-(1 / 2))).sub (P2.mul P2m)).sub P3).sub ((P4.neg.add P5).add P6)
  convert H using 1
  · funext x; simp only [F, l1, mm, Function.comp_def, Pi.sub_apply, Pi.add_apply, Pi.neg_apply,
      Pi.mul_apply]
  · have key := mm_eq hv
    unfold mm at key
    rw [sub_sub_cancel, sub_neg_eq_add]
    have e4 : -Real.log (1 - w v / Real.sqrt 2) / (w v / Real.sqrt 2)
        * (-2 * a / (a + v) ^ 2 / Real.sqrt 2)
        = -Real.log (1 - w v / Real.sqrt 2) * ((-2 * a / (a + v) ^ 2) / w v) := by
      field_simp
    rw [e4]
    unfold l1
    linear_combination (-Real.log (tv v)) * hX hv
      - (-(((v ^ 2 + 1) / (4 * v ^ 2)) / tv v) + ((1 - v ^ 2) / (4 * Real.sqrt 2 * v ^ 2)) / rho v)
        * key
      - (Real.log (1 - w v / Real.sqrt 2) - Real.log (1 - w v) - Real.log (1 + w v)) * hY hv

lemma Rwt_eq {v : ℝ} (hv : v ∈ Ioo b a) :
    Rwt v = (v - b) * (a ^ 2 + 2 * Real.sqrt 2 * a * v - v ^ 2)
      / (2 * v * (v + b) * (v + a) * (v - a)) := by
  obtain ⟨h0, h1, h2, h3, h4⟩ := mem_facts hv
  have h5 : v - a ≠ 0 := by linarith
  have h6 : v + a ≠ 0 := by linarith
  unfold Rwt
  field_simp
  unfold a b; ring_nf; simp only [s_sq, s3]; ring_nf

lemma Rwt_nonpos {v : ℝ} (hv : v ∈ Ioo b a) : Rwt v ≤ 0 := by
  obtain ⟨h0, h1, h2, h3, h4⟩ := mem_facts hv
  rw [Rwt_eq hv]
  have hs := s_pos
  have ha := a_pos
  apply div_nonpos_of_nonneg_of_nonpos
  · apply mul_nonneg h3.le
    nlinarith [mul_pos (mul_pos hs ha) h0]
  · have : 0 < 2 * v * (v + b) * (v + a) := by positivity
    nlinarith

lemma l1_nonpos {v : ℝ} (hv : v ∈ Ioo b a) : l1 v ≤ 0 := by
  obtain ⟨t0, t1⟩ := tv_mem hv
  exact Real.log_nonpos t0.le t1.le

lemma deriv_nonneg {v : ℝ} (hv : v ∈ Ioo b a) : 0 ≤ l1 v * Rwt v :=
  mul_nonneg_of_nonpos_of_nonpos (l1_nonpos hv) (Rwt_nonpos hv)

lemma memI {v : ℝ} (hv : v ∈ Icc b a) :
    0 < v ∧ 0 ≤ a - v ∧ 0 < a + v ∧ 0 ≤ v - b ∧ 0 < v + b := by
  obtain ⟨h1, h2⟩ := hv
  have := b_pos; have := a_pos
  refine ⟨by linarith, by linarith, by linarith, by linarith, by linarith⟩

lemma tvI {v : ℝ} (hv : v ∈ Icc b a) : 0 ≤ tv v ∧ tv v ≤ 1 := by
  obtain ⟨h0, h1, h2, h3, h4⟩ := memI hv
  have e1 := tv_eq h0.ne'
  have e2 := one_sub_tv h0.ne'
  constructor
  · rw [e1]; exact div_nonneg (mul_nonneg h3 (by linarith)) (by linarith)
  · have : 0 ≤ 1 - tv v := by rw [e2]; exact div_nonneg (mul_nonneg h1 h4.le) (by linarith)
    linarith

lemma wI {v : ℝ} (hv : v ∈ Icc b a) : 0 ≤ w v / Real.sqrt 2 ∧ w v / Real.sqrt 2 ≤ 1 / 2 := by
  obtain ⟨h0, h1, h2, h3, h4⟩ := memI hv
  have hs := s_pos
  have hs1 := s_gt_one
  unfold w
  constructor
  · positivity
  · rw [div_div, div_le_div_iff₀ (by positivity) (by norm_num)]
    have hb : b ≤ v := hv.1
    unfold a b at *
    nlinarith [s_sq]

lemma rhoI {v : ℝ} (hv : v ∈ Icc b a) : 0 ≤ rho v ∧ 2 * rho v ≤ tv v := by
  obtain ⟨h0, h1, h2, h3, h4⟩ := memI hv
  obtain ⟨t0, t1⟩ := tvI hv
  obtain ⟨q0, q1⟩ := wI hv
  rw [rho_eq h0.ne' h2.ne']
  exact ⟨mul_nonneg t0 q0, by nlinarith⟩

lemma mm_bound {v : ℝ} (hv : v ∈ Icc b a) : |mm v| ≤ tv v := by
  obtain ⟨r0, r1⟩ := rhoI hv
  obtain ⟨t0, t1⟩ := tvI hv
  have hp : 0 < 1 - rho v := by linarith
  have hle : mm v ≤ 0 := Real.log_nonpos hp.le (by linarith)
  have hge := Real.one_sub_inv_le_log_of_pos hp
  have : -(2 * rho v) ≤ 1 - (1 - rho v)⁻¹ := by
    rw [show 1 - (1 - rho v)⁻¹ = -(rho v / (1 - rho v)) by field_simp; ring]
    rw [neg_le_neg_iff, div_le_iff₀ hp]
    nlinarith
  unfold mm at hle ⊢
  rw [abs_of_nonpos hle]
  linarith

lemma tv_b : tv b = 0 := by
  unfold tv b; ring_nf; simp only [s_sq]; ring

lemma tv_a : tv a = 1 := by
  have := a_pos
  unfold tv a; field_simp; ring_nf; simp only [s_sq]; ring

lemma rho_b : rho b = 0 := by unfold rho; simp
lemma rho_a : rho a = 0 := by unfold rho; simp
lemma w_a : w a = 0 := by unfold w; simp

lemma w_b : w b = 1 / Real.sqrt 2 := by
  have hs := s_pos.ne'
  unfold w a b; field_simp; ring_nf; exact mul_inv_cancel₀ hs

lemma tvC : ContinuousOn tv (Icc b a) := by
  have h : ContinuousOn (fun v : ℝ => (v ^ 2 + 2 * v - 1) / (4 * v)) (Icc b a) :=
    ContinuousOn.div (f := fun v : ℝ => v ^ 2 + 2 * v - 1) (g := fun v : ℝ => 4 * v) (by fun_prop) (by fun_prop) (fun v hv => by
      have := (memI hv).1; positivity)
  exact h

lemma rhoC : ContinuousOn rho (Icc b a) := by
  have := s_pos
  have h : ContinuousOn (fun v : ℝ => (a - v) * (v - b) / (4 * Real.sqrt 2 * v)) (Icc b a) :=
    ContinuousOn.div (f := fun v : ℝ => (a - v) * (v - b)) (g := fun v : ℝ => 4 * Real.sqrt 2 * v)
      (by fun_prop) (by fun_prop) (fun v hv => by
      have := (memI hv).1; positivity)
  exact h

lemma wC : ContinuousOn w (Icc b a) := by
  have h : ContinuousOn (fun v : ℝ => (a - v) / (a + v)) (Icc b a) :=
    ContinuousOn.div (f := fun v : ℝ => a - v) (g := fun v : ℝ => a + v) (by fun_prop) (by fun_prop) (fun v hv => by
      have := (memI hv).2.2.1; linarith)
  exact h

lemma mmC : ContinuousOn mm (Icc b a) := by
  unfold mm
  refine ContinuousOn.log (continuousOn_const.sub rhoC) (fun v hv => ?_)
  have := rhoI hv; have := tvI hv
  linarith

lemma tv_pos_Ioc {v : ℝ} (hv : v ∈ Icc b a) (hvb : v ≠ b) : 0 < tv v := by
  obtain ⟨h0, h1, h2, h3, h4⟩ := memI hv
  have h3' : 0 < v - b := lt_of_le_of_ne h3 (fun h => hvb (by linarith))
  rw [tv_eq h0.ne']
  exact div_pos (mul_pos h3' (by linarith)) (by linarith)

lemma prod_at_b : ContinuousWithinAt (fun v => l1 v * mm v) (Icc b a) b := by
  rw [ContinuousWithinAt, show l1 b * mm b = 0 by simp [l1, tv_b]]
  have hb : b ∈ Icc b a := ⟨le_refl _, b_lt_a.le⟩
  have hc : ContinuousWithinAt ((fun x => |x * Real.log x|) ∘ tv) (Icc b a) b :=
    ((continuous_abs.comp Real.continuous_mul_log).continuousAt).comp_continuousWithinAt
      (tvC b hb)
  have ht := hc.tendsto
  simp only [Function.comp_def, tv_b, zero_mul, abs_zero] at ht
  apply squeeze_zero_norm' _ ht
  filter_upwards [self_mem_nhdsWithin] with v hv
  rw [Real.norm_eq_abs, abs_mul, abs_mul, abs_of_nonneg (tvI hv).1, mul_comm (tv v)]
  exact mul_le_mul_of_nonneg_left (mm_bound hv) (abs_nonneg _)

lemma prodC : ContinuousOn (fun v => l1 v * mm v) (Icc b a) := by
  intro v hv
  by_cases hvb : v = b
  · subst hvb; exact prod_at_b
  · apply ContinuousWithinAt.mul _ (mmC v hv)
    have h0 := (memI hv).1
    have htc : ContinuousAt tv v := by
      have h : ContinuousAt (fun v : ℝ => (v ^ 2 + 2 * v - 1) / (4 * v)) v :=
        ContinuousAt.div (f := fun v : ℝ => v ^ 2 + 2 * v - 1) (g := fun v : ℝ => 4 * v)
          (by fun_prop) (by fun_prop) (by positivity)
      exact h
    exact ((Real.continuousAt_log (tv_pos_Ioc hv hvb).ne').comp htc).continuousWithinAt

lemma FC : ContinuousOn F (Icc b a) := by
  have L := D00.continuous_Li2
  have c1 : ContinuousOn (fun v => D00.Li2 (1 - tv v)) (Icc b a) :=
    L.comp_continuousOn ((continuousOn_const (c := (1:ℝ))).sub tvC)
  have c3 : ContinuousOn (fun v => D00.Li2 (rho v)) (Icc b a) := L.comp_continuousOn rhoC
  have c4 : ContinuousOn (fun v => D00.Li2 (w v / Real.sqrt 2)) (Icc b a) :=
    L.comp_continuousOn (wC.div_const (Real.sqrt 2))
  have c5 : ContinuousOn (fun v => D00.Li2 (w v)) (Icc b a) := L.comp_continuousOn wC
  have c6 : ContinuousOn (fun v => D00.Li2 (-w v)) (Icc b a) := L.comp_continuousOn wC.neg
  have h := ((((continuousOn_const (c := -(1 / 2 : ℝ))).mul c1).sub prodC).sub c3).sub
    ((c4.neg.add c5).add c6)
  exact h

lemma F_a : F a = 0 := by
  unfold F l1
  rw [tv_a, rho_a, w_a]
  simp

lemma F_b : F b = -(π ^ 2 / 12) + (1 / 2) * D00.Li2 (1 / 2) := by
  have hs := s_pos
  have h2 : (1 / Real.sqrt 2) / Real.sqrt 2 = 1 / 2 := by
    rw [div_div, ← sq, s_sq]
  have hd := D00.duplication (x := 1 / Real.sqrt 2) (by
    rw [abs_of_pos (by positivity), div_le_one hs]; exact s_gt_one.le)
  have h3 : (1 / Real.sqrt 2) ^ 2 = 1 / 2 := by rw [div_pow, s_sq]; norm_num
  rw [h3] at hd
  unfold F l1
  rw [tv_b, rho_b, w_b, h2]
  simp only [sub_zero, D00.Li2_one, Real.log_zero, zero_mul, D00.Li2_zero]
  linear_combination (-1 : ℝ) * hd

lemma Tint_eq : Tint = π ^ 2 / 12 - (1 / 2) * D00.Li2 (1 / 2) := by
  have hab := b_lt_a.le
  have hint : IntervalIntegrable (fun v => l1 v * Rwt v) volume b a :=
    (intervalIntegrable_iff_integrableOn_Ioc_of_le hab).2
      (intervalIntegral.integrableOn_deriv_of_nonneg FC (fun x hx => hasDerivAt_F hx) (fun x hx => deriv_nonneg hx))
  unfold Tint
  rw [intervalIntegral.integral_eq_sub_of_hasDerivAt_of_le hab FC (fun x hx => hasDerivAt_F hx)
    hint, F_a, F_b]
  ring

theorem main : Y01.stmt_8 := by
  have h := D00.Li2_half
  refine ⟨?_, ?_⟩
  · rw [Tint_eq, h]; unfold lam; ring
  · rw [h]; unfold lam; ring

end Y01.PD

namespace Y01

noncomputable def rawW3 : ℝ :=
  -- u = -a²
  (0 : ℝ) * Li3 (-(a ^ 2))
    + (1 / 4 * Real.sqrt 2 * lam + 1 / 8 * Real.sqrt 2 * L) * D00.Li2 (-(a ^ 2))
  -- u = -2a
  + (-1 + 5 / 8 * Real.sqrt 2) * Li3 (-(2 * a))
    + ((-1 + 5 / 8 * Real.sqrt 2) * lam + 0 * L) * D00.Li2 (-(2 * a))
  -- u = -a
  + (1 / 2 - 1 / 4 * Real.sqrt 2) * Li3 (-a)
    + (0 * lam + (1 / 2 - 1 / 4 * Real.sqrt 2) * L) * D00.Li2 (-a)
  -- u = -1
  + (-2) * Li3 (-1)
    + ((-1 / 2 - 3 / 2 * Real.sqrt 2) * lam + (2 - Real.sqrt 2) * L) * D00.Li2 (-1)
  -- u = -1/√2
  + (2 + 1 / 4 * Real.sqrt 2) * Li3 (-(1 / Real.sqrt 2))
    + ((2 + 1 / 8 * Real.sqrt 2) * lam + 0 * L) * D00.Li2 (-(1 / Real.sqrt 2))
  -- u = -b
  + (1 / 2 - 1 / 2 * Real.sqrt 2) * Li3 (-b)
    + (-3 / 2 * lam + (1 / 2 + 1 / 2 * Real.sqrt 2) * L) * D00.Li2 (-b)
  -- u = -b/2
  + (-1 / 8 * Real.sqrt 2) * Li3 (-(b / 2))
    + (-1 / 8 * Real.sqrt 2 * lam - 1 / 8 * Real.sqrt 2 * L) * D00.Li2 (-(b / 2))
  -- u = -b²
  + 2 * Li3 (-(b ^ 2))
    + ((2 + 5 / 4 * Real.sqrt 2) * lam + (3 - 5 / 8 * Real.sqrt 2) * L) * D00.Li2 (-(b ^ 2))
  -- u = (2-√2)/4
  + (1 + 1 / 2 * Real.sqrt 2) * Li3 ((2 - Real.sqrt 2) / 4)
    + ((3 / 2 + 3 / 4 * Real.sqrt 2) * lam + (1 + 1 / 2 * Real.sqrt 2) * L)
      * D00.Li2 ((2 - Real.sqrt 2) / 4)
  -- u = b²
  + (0 : ℝ) * Li3 (b ^ 2)
    + (-1 / 4 * Real.sqrt 2 * lam + 1 / 8 * Real.sqrt 2 * L) * D00.Li2 (b ^ 2)
  -- u = 1 - 1/√2
  + (1 / 4 * Real.sqrt 2) * Li3 (1 - 1 / Real.sqrt 2)
    + ((3 / 2 + 1 / 8 * Real.sqrt 2) * lam + (-1 + 1 / 4 * Real.sqrt 2) * L)
      * D00.Li2 (1 - 1 / Real.sqrt 2)
  -- u = b
  + (1 + 3 / 4 * Real.sqrt 2) * Li3 b
    + (0 * lam + (1 + 1 / 2 * Real.sqrt 2) * L) * D00.Li2 b
  -- u = 1/2
  + (-2) * Li3 (1 / 2)
    + (-1 / 2 * lam - L) * D00.Li2 (1 / 2)
  -- u = 2 - √2
  + (-1 / 4 * Real.sqrt 2) * Li3 (2 - Real.sqrt 2)
    + (1 / 8 * Real.sqrt 2 * lam - 1 / 4 * Real.sqrt 2 * L) * D00.Li2 (2 - Real.sqrt 2)
  -- u = 1/√2
  + (2 - 1 / 4 * Real.sqrt 2) * Li3 (1 / Real.sqrt 2)
    + ((1 - 1 / 8 * Real.sqrt 2) * lam + 1 / 4 * Real.sqrt 2 * L) * D00.Li2 (1 / Real.sqrt 2)
  -- u = 2b
  + (-1 - 1 / 2 * Real.sqrt 2) * Li3 (2 * b)
    + ((1 + 1 / 2 * Real.sqrt 2) * lam + (-1 - 1 / 2 * Real.sqrt 2) * L) * D00.Li2 (2 * b)
  -- u = (2+√2)/4
  + (1 - 1 / 2 * Real.sqrt 2) * Li3 ((2 + Real.sqrt 2) / 4)
    + ((3 / 2 - 3 / 4 * Real.sqrt 2) * lam + (-1 + 1 / 2 * Real.sqrt 2) * L)
      * D00.Li2 ((2 + Real.sqrt 2) / 4)
  -- terms without Li_n
  + (1 / 3 - 1 / 48 * Real.sqrt 2) * Real.pi ^ 2 * L
  + (-1 / 2 + 1 / 24 * Real.sqrt 2) * Real.pi ^ 2 * lam
  + (4 / 3 - 3 / 4 * Real.sqrt 2) * L ^ 3
  + (-2 + 31 / 8 * Real.sqrt 2) * lam * L ^ 2
  + (-41 / 8 - 21 / 32 * Real.sqrt 2) * lam ^ 2 * L
  + (143 / 48 - 1 / 32 * Real.sqrt 2) * lam ^ 3

end Y01

/-!
# Part E1 of Y01: `W_3 = rawW3` (supplement Lemmas g1:lem:prim, g1:lem:elem, g1:lem:raw)

`W_3 = ∫_b^a W3integrand` is evaluated exactly as in `scripts/g1_M2_raw.py`:
`ℓ₁ = ℓ_b + ℓ_{-a} - ℓ_0 - 2λ`, `ℓ₂ = ℓ_a + ℓ_{-b} - ℓ_0 - 2λ` (`ℓ_x = log|v - x|`), every product
`ℓ_x ℓ_y` is polarised with the endpoint letter in the numerator of `w = (v-x)/(v-y)`, and `[b,a]`
is split at `1` where the two pieces use different orientations.  This gives 43 elementary
integrals of the forms `ℓ_x²/(v-z)`, `log²|w|/(v-z)`, `ℓ_x/(v-z)`, `1/(v-z)` (`z ∈ {0,-a}`).

1. `P1n, P1p, P2a, P2g` are the primitives `P₁, P₂` of g1:lem:prim; `Phi`, `Psi` combine them
   with the constant `κ`; `FA, FC, FG, FIV, FIII, FL, FH, FK` are the primitives of
   g1:lem:elem in the variable `v`, with `HasDerivAt` lemmas `dFA ... dFK`.
2. `ftc_pos` (nonnegative integrand, continuous primitive on the closed piece) and the lemmas
   `intA ... intK` turn a primitive into the value of the integral (one-sided limits are
   handled by continuity at `s = 0`: `contAt_P1n`, `contAt_P2a`).
3. `blkN` are the 43 instances, `valNq`, `valNp` their endpoint values in the atoms of
   Table g1:tab:rawW3.
4. `pw1`, `pw2` decompose the integrand on `(b,1)` and `(1,a)`; `main` sums everything.
-/

set_option linter.unusedSimpArgs false
set_option linter.unusedTactic false
set_option linter.unusedVariables false

namespace Y01
namespace PE1

open Real MeasureTheory Set Filter Topology

/-! ## The primitives of Lemma g1:lem:prim -/

noncomputable def P1n (s : ℝ) : ℝ := Real.log s * Real.log (1 - s) + D00.Li2 s
noncomputable def P1p (s : ℝ) : ℝ := -D00.Li2 (1 - s)
noncomputable def P2a (s : ℝ) : ℝ :=
  Real.log s ^ 2 * Real.log (1 - s) + 2 * Real.log s * D00.Li2 s - 2 * Li3 s
noncomputable def E3 (r : ℝ) : ℝ :=
  -(1 / 3) * Real.log r ^ 3 + Real.log r ^ 2 * Real.log (1 - r)
    + 2 * Real.log r * D00.Li2 r - 2 * Li3 r
noncomputable def P2g (s : ℝ) : ℝ := E3 s⁻¹
noncomputable def Phi (P1 P2 : ℝ → ℝ) (k s : ℝ) : ℝ :=
  k ^ 2 * Real.log (1 - s) + 2 * k * P1 s + P2 s
noncomputable def Psi (P1 : ℝ → ℝ) (k s : ℝ) : ℝ := k * Real.log (1 - s) + P1 s

lemma hasDerivAt_log1 {s : ℝ} (h1 : s ≠ 1) :
    HasDerivAt (fun s => Real.log (1 - s)) (1 / (s - 1)) s := by
  have h : (1 : ℝ) - s ≠ 0 := sub_ne_zero.mpr (Ne.symm h1)
  have := ((hasDerivAt_id s).const_sub 1).log h
  refine this.congr_deriv ?_
  have h' : s - 1 ≠ 0 := sub_ne_zero.mpr h1
  simp only [id]
  field_simp
  ring

lemma hasDerivAt_P1n {s : ℝ} (h0 : s ≠ 0) (h1 : s ≠ 1) :
    HasDerivAt P1n (Real.log s / (s - 1)) s := by
  have hl := Real.hasDerivAt_log h0
  have hl1 := hasDerivAt_log1 h1
  have hL := D00.hasDerivAt_Li2_ne h1 h0
  have := (hl.mul hl1).add hL
  refine this.congr_deriv ?_
  have h' : s - 1 ≠ 0 := sub_ne_zero.mpr h1
  have h'' : 1 - s ≠ 0 := sub_ne_zero.mpr (Ne.symm h1)
  field_simp
  ring

lemma hasDerivAt_P1p {s : ℝ} (h0 : s ≠ 0) (h1 : s ≠ 1) :
    HasDerivAt P1p (Real.log s / (s - 1)) s := by
  have ha : (1 : ℝ) - s ≠ 1 := by intro h; apply h0; linarith
  have hb : (1 : ℝ) - s ≠ 0 := sub_ne_zero.mpr (Ne.symm h1)
  have := (D00.Li2_comp ((hasDerivAt_id s).const_sub 1) ha hb).neg
  unfold P1p
  refine this.congr_deriv ?_
  simp only [id, sub_sub_cancel]
  have h' : s - 1 ≠ 0 := sub_ne_zero.mpr h1
  field_simp
  ring

lemma hasDerivAt_P2a {s : ℝ} (h0 : s ≠ 0) (h1 : s ≠ 1) :
    HasDerivAt P2a (Real.log s ^ 2 / (s - 1)) s := by
  have hl := Real.hasDerivAt_log h0
  have hl1 := hasDerivAt_log1 h1
  have hL := D00.hasDerivAt_Li2_ne h1 h0
  have h3 := hasDerivAt_Li3 h0
  have := (((hl.pow 2).mul hl1).add ((hl.const_mul 2).mul hL)).sub (h3.const_mul 2)
  unfold P2a
  refine this.congr_deriv ?_
  have h' : s - 1 ≠ 0 := sub_ne_zero.mpr h1
  have h'' : 1 - s ≠ 0 := sub_ne_zero.mpr (Ne.symm h1)
  simp only [Pi.pow_apply]
  field_simp
  ring

lemma hasDerivAt_E3 {r : ℝ} (h0 : r ≠ 0) (h1 : r ≠ 1) :
    HasDerivAt E3 (-(Real.log r ^ 2) / (r * (1 - r))) r := by
  have hl := Real.hasDerivAt_log h0
  have hl1 := hasDerivAt_log1 h1
  have hL := D00.hasDerivAt_Li2_ne h1 h0
  have h3 := hasDerivAt_Li3 h0
  have := ((((hl.pow 3).const_mul (-(1 / 3))).add ((hl.pow 2).mul hl1)).add
    ((hl.const_mul 2).mul hL)).sub (h3.const_mul 2)
  unfold E3
  refine this.congr_deriv ?_
  have h' : r - 1 ≠ 0 := sub_ne_zero.mpr h1
  have h'' : 1 - r ≠ 0 := sub_ne_zero.mpr (Ne.symm h1)
  simp only [Pi.pow_apply]
  field_simp
  ring

lemma hasDerivAt_P2g {s : ℝ} (h0 : s ≠ 0) (h1 : s ≠ 1) :
    HasDerivAt P2g (Real.log s ^ 2 / (s - 1)) s := by
  have hi0 : s⁻¹ ≠ 0 := inv_ne_zero h0
  have hi1 : s⁻¹ ≠ 1 := by
    intro h; apply h1; rw [inv_eq_one] at h; exact h
  have := (hasDerivAt_E3 hi0 hi1).comp s (hasDerivAt_inv h0)
  unfold P2g
  refine this.congr_deriv ?_
  rw [Real.log_inv]
  have h' : s - 1 ≠ 0 := sub_ne_zero.mpr h1
  field_simp

lemma hasDerivAt_Phi {P1 P2 : ℝ → ℝ} {k s : ℝ} (h1 : s ≠ 1)
    (hP1 : HasDerivAt P1 (Real.log s / (s - 1)) s)
    (hP2 : HasDerivAt P2 (Real.log s ^ 2 / (s - 1)) s) :
    HasDerivAt (Phi P1 P2 k) ((k + Real.log s) ^ 2 / (s - 1)) s := by
  have := (((hasDerivAt_log1 h1).const_mul (k ^ 2)).add (hP1.const_mul (2 * k))).add hP2
  unfold Phi
  refine this.congr_deriv ?_
  have h' : s - 1 ≠ 0 := sub_ne_zero.mpr h1
  field_simp
  ring

lemma hasDerivAt_Psi {P1 : ℝ → ℝ} {k s : ℝ} (h1 : s ≠ 1)
    (hP1 : HasDerivAt P1 (Real.log s / (s - 1)) s) :
    HasDerivAt (Psi P1 k) ((k + Real.log s) / (s - 1)) s := by
  have := ((hasDerivAt_log1 h1).const_mul k).add hP1
  unfold Psi
  refine this.congr_deriv ?_
  have h' : s - 1 ≠ 0 := sub_ne_zero.mpr h1
  field_simp


/-! ## Continuity of the primitives -/

lemma xlog2_eq (s : ℝ) :
    |s * Real.log s ^ 2| = (2 * (Real.sqrt |s| * Real.log (Real.sqrt |s|))) ^ 2 := by
  rw [Real.log_sqrt (abs_nonneg s), abs_mul, abs_pow, sq_abs, ← Real.log_abs s]
  have := Real.sq_sqrt (abs_nonneg s)
  rw [mul_pow, mul_pow, this]
  ring

lemma contAt_xlog2 : ContinuousAt (fun s : ℝ => s * Real.log s ^ 2) 0 := by
  have hg : Continuous (fun s : ℝ => (2 * (Real.sqrt |s| * Real.log (Real.sqrt |s|))) ^ 2) :=
    (continuous_const.mul (Real.continuous_mul_log.comp
      (Real.continuous_sqrt.comp continuous_abs))).pow 2
  have h0 : Tendsto (fun s : ℝ => (2 * (Real.sqrt |s| * Real.log (Real.sqrt |s|))) ^ 2)
      (𝓝 0) (𝓝 0) := by
    have := hg.tendsto 0
    simpa using this
  rw [ContinuousAt]
  simp only [zero_mul]
  apply squeeze_zero_norm _ h0
  intro s
  rw [Real.norm_eq_abs, xlog2_eq]

lemma contAt_xlog : ContinuousAt (fun s : ℝ => s * Real.log s) 0 :=
  Real.continuous_mul_log.continuousAt

lemma log1_eq_G (s : ℝ) : Real.log (1 - s) = s * D00.G s := by
  rcases eq_or_ne s 0 with rfl | h
  · simp
  · rw [D00.G_of_ne h]; field_simp

lemma Li2_eq_H (s : ℝ) : D00.Li2 s = s * H s := by
  rcases eq_or_ne s 0 with rfl | h
  · simp
  · rw [H_of_ne h]; field_simp

lemma P1n_eq : P1n = fun s => (s * Real.log s) * D00.G s + D00.Li2 s := by
  funext s; rw [P1n, log1_eq_G]; ring

lemma P2a_eq : P2a = fun s => (s * Real.log s ^ 2) * D00.G s
    + 2 * ((s * Real.log s) * H s) - 2 * Li3 s := by
  funext s; rw [P2a, log1_eq_G, Li2_eq_H]; ring

lemma contAt_log1 {s : ℝ} (h1 : s ≠ 1) : ContinuousAt (fun s => Real.log (1 - s)) s :=
  (hasDerivAt_log1 h1).continuousAt

lemma contAt_P1n {s : ℝ} (h1 : s ≠ 1) : ContinuousAt P1n s := by
  rcases eq_or_ne s 0 with rfl | h0
  · rw [P1n_eq]
    exact (contAt_xlog.mul D00.continuousAt_G_zero).add D00.continuous_Li2.continuousAt
  · exact (hasDerivAt_P1n h0 h1).continuousAt

lemma contAt_P1p (s : ℝ) : ContinuousAt P1p s := by
  unfold P1p
  exact (D00.continuous_Li2.comp (continuous_const.sub continuous_id)).neg.continuousAt

lemma contAt_P2a {s : ℝ} (h1 : s ≠ 1) : ContinuousAt P2a s := by
  rcases eq_or_ne s 0 with rfl | h0
  · rw [P2a_eq]
    exact ((contAt_xlog2.mul D00.continuousAt_G_zero).add
      (continuousAt_const.mul (contAt_xlog.mul continuous_H.continuousAt))).sub
      (continuousAt_const.mul continuous_Li3.continuousAt)
  · exact (hasDerivAt_P2a h0 h1).continuousAt

lemma contAt_P2g {s : ℝ} (h0 : s ≠ 0) (h1 : s ≠ 1) : ContinuousAt P2g s :=
  (hasDerivAt_P2g h0 h1).continuousAt

lemma contAt_Phi {P1 P2 : ℝ → ℝ} {k s : ℝ} (h1 : s ≠ 1) (c1 : ContinuousAt P1 s)
    (c2 : ContinuousAt P2 s) : ContinuousAt (Phi P1 P2 k) s := by
  unfold Phi
  exact ((continuousAt_const.mul (contAt_log1 h1)).add (continuousAt_const.mul c1)).add c2

lemma contAt_Psi {P1 : ℝ → ℝ} {k s : ℝ} (h1 : s ≠ 1) (c1 : ContinuousAt P1 s) :
    ContinuousAt (Psi P1 k) s := by
  unfold Psi
  exact (continuousAt_const.mul (contAt_log1 h1)).add c1


/-! ## Elementary integrands (Lemma g1:lem:elem) and their primitives -/

noncomputable def fsq (x z v : ℝ) : ℝ := Real.log (v - x) ^ 2 / (v - z)
noncomputable def fpr (x y z v : ℝ) : ℝ := (Real.log (v - x) - Real.log (v - y)) ^ 2 / (v - z)
noncomputable def flin (x z v : ℝ) : ℝ := Real.log (v - x) / (v - z)
noncomputable def fcst (z v : ℝ) : ℝ := 1 / (v - z)

def IsP1 (P : ℝ → ℝ) : Prop := ∀ s, s ≠ 0 → s ≠ 1 → HasDerivAt P (Real.log s / (s - 1)) s
def IsP2 (P : ℝ → ℝ) : Prop := ∀ s, s ≠ 0 → s ≠ 1 → HasDerivAt P (Real.log s ^ 2 / (s - 1)) s

lemma isP1n : IsP1 P1n := fun _ h0 h1 => hasDerivAt_P1n h0 h1
lemma isP1p : IsP1 P1p := fun _ h0 h1 => hasDerivAt_P1p h0 h1
lemma isP2a : IsP2 P2a := fun _ h0 h1 => hasDerivAt_P2a h0 h1
lemma isP2g : IsP2 P2g := fun _ h0 h1 => hasDerivAt_P2g h0 h1

noncomputable def sA (x z v : ℝ) : ℝ := (v - x) / (z - x)
noncomputable def FA (P1 P2 : ℝ → ℝ) (x z v : ℝ) : ℝ := Phi P1 P2 (Real.log (z - x)) (sA x z v)
noncomputable def FC (x v : ℝ) : ℝ := Real.log (v - x) ^ 3 / 3
noncomputable def s1 (x y z v : ℝ) : ℝ := (v - x) * (z - y) / ((v - y) * (z - x))
noncomputable def mb (x y v : ℝ) : ℝ := (v - x) / (v - y)
noncomputable def FG (P1 P2 Q : ℝ → ℝ) (x y z v : ℝ) : ℝ :=
  Phi P1 P2 (Real.log ((z - x) / (z - y))) (s1 x y z v) - Q (mb x y v)
noncomputable def FIV (Q : ℝ → ℝ) (x y v : ℝ) : ℝ := -Q (mb x y v)
noncomputable def FIII (Q : ℝ → ℝ) (x y v : ℝ) : ℝ :=
  (Real.log (v - x) - Real.log (v - y)) ^ 3 / 3 - Q (mb x y v)
noncomputable def FL (P1 : ℝ → ℝ) (x z v : ℝ) : ℝ := Psi P1 (Real.log (z - x)) (sA x z v)
noncomputable def FH (x v : ℝ) : ℝ := Real.log (v - x) ^ 2 / 2
noncomputable def FK (z v : ℝ) : ℝ := Real.log (v - z)

lemma sA_ne_zero {x z v : ℝ} (hxz : z - x ≠ 0) (hvx : v - x ≠ 0) : sA x z v ≠ 0 :=
  div_ne_zero hvx hxz

lemma sA_ne_one {x z v : ℝ} (hxz : z - x ≠ 0) (hvz : v - z ≠ 0) : sA x z v ≠ 1 := by
  intro h
  rw [sA, div_eq_one_iff_eq hxz] at h
  apply hvz; linarith

lemma s1_ne_zero {x y z v : ℝ} (hzx : z - x ≠ 0) (hzy : z - y ≠ 0) (hvx : v - x ≠ 0)
    (hvy : v - y ≠ 0) : s1 x y z v ≠ 0 :=
  div_ne_zero (mul_ne_zero hvx hzy) (mul_ne_zero hvy hzx)

lemma s1_ne_one {x y z v : ℝ} (hxy : x - y ≠ 0) (hzx : z - x ≠ 0) (hvy : v - y ≠ 0)
    (hvz : v - z ≠ 0) : s1 x y z v ≠ 1 := by
  intro h
  rw [s1, div_eq_one_iff_eq (mul_ne_zero hvy hzx)] at h
  have : (v - z) * (x - y) = 0 := by linear_combination h
  rcases mul_eq_zero.mp this with h' | h'
  · exact hvz h'
  · exact hxy h'

lemma mb_ne_zero {x y v : ℝ} (hvx : v - x ≠ 0) (hvy : v - y ≠ 0) : mb x y v ≠ 0 :=
  div_ne_zero hvx hvy

lemma mb_ne_one {x y v : ℝ} (hxy : x - y ≠ 0) (hvy : v - y ≠ 0) : mb x y v ≠ 1 := by
  intro h
  rw [mb, div_eq_one_iff_eq hvy] at h
  apply hxy; linarith

lemma hasDerivAt_sA (x z v : ℝ) : HasDerivAt (sA x z) (1 / (z - x)) v := by
  have := ((hasDerivAt_id v).sub_const x).div_const (z - x)
  unfold sA
  exact this

lemma hasDerivAt_mb {x y v : ℝ} (hvy : v - y ≠ 0) :
    HasDerivAt (mb x y) ((x - y) / (v - y) ^ 2) v := by
  have := ((hasDerivAt_id v).sub_const x).div ((hasDerivAt_id v).sub_const y) hvy
  refine this.congr_deriv ?_
  simp only [id]
  ring

lemma hasDerivAt_s1 {x y z v : ℝ} (hzx : z - x ≠ 0) (hvy : v - y ≠ 0) :
    HasDerivAt (s1 x y z) ((z - y) * (x - y) / ((v - y) ^ 2 * (z - x))) v := by
  have := (((hasDerivAt_id v).sub_const x).mul_const (z - y)).div
    (((hasDerivAt_id v).sub_const y).mul_const (z - x)) (mul_ne_zero hvy hzx)
  refine this.congr_deriv ?_
  simp only [id]
  field_simp
  ring

lemma log_mb {x y v : ℝ} (hvx : v - x ≠ 0) (hvy : v - y ≠ 0) :
    Real.log (mb x y v) = Real.log (v - x) - Real.log (v - y) := Real.log_div hvx hvy

lemma log_s1 {x y z v : ℝ} (hzx : z - x ≠ 0) (hzy : z - y ≠ 0) (hvx : v - x ≠ 0)
    (hvy : v - y ≠ 0) : Real.log (s1 x y z v)
      = Real.log (v - x) + Real.log (z - y) - (Real.log (v - y) + Real.log (z - x)) := by
  rw [s1, Real.log_div (mul_ne_zero hvx hzy) (mul_ne_zero hvy hzx), Real.log_mul hvx hzy,
    Real.log_mul hvy hzx]

lemma dFA {P1 P2 : ℝ → ℝ} (h1 : IsP1 P1) (h2 : IsP2 P2) {x z v : ℝ} (hxz : z - x ≠ 0)
    (hvx : v - x ≠ 0) (hvz : v - z ≠ 0) : HasDerivAt (FA P1 P2 x z) (fsq x z v) v := by
  have n0 := sA_ne_zero hxz hvx
  have n1 := sA_ne_one hxz hvz
  have := (hasDerivAt_Phi (k := Real.log (z - x)) n1 (h1 _ n0 n1) (h2 _ n0 n1)).comp v
    (hasDerivAt_sA x z v)
  refine this.congr_deriv ?_
  rw [sA, Real.log_div hvx hxz]
  unfold fsq
  have e : (v - x) / (z - x) - 1 = (v - z) / (z - x) := by field_simp; ring
  rw [e]
  field_simp
  ring

lemma dFC {x v : ℝ} (hvx : v - x ≠ 0) : HasDerivAt (FC x) (fsq x x v) v := by
  have := ((((hasDerivAt_id v).sub_const x).log hvx).pow 3).div_const 3
  unfold FC fsq
  refine this.congr_deriv ?_
  simp only [id]
  field_simp
  ring

lemma dFG {P1 P2 Q : ℝ → ℝ} (h1 : IsP1 P1) (h2 : IsP2 P2) (hQ : IsP2 Q) {x y z v : ℝ}
    (hxy : x - y ≠ 0) (hzx : z - x ≠ 0) (hzy : z - y ≠ 0) (hvx : v - x ≠ 0) (hvy : v - y ≠ 0)
    (hvz : v - z ≠ 0) : HasDerivAt (FG P1 P2 Q x y z) (fpr x y z v) v := by
  have n0 := s1_ne_zero hzx hzy hvx hvy
  have n1 := s1_ne_one hxy hzx hvy hvz
  have m0 := mb_ne_zero hvx hvy
  have m1 := mb_ne_one hxy hvy
  have dA := (hasDerivAt_Phi (k := Real.log ((z - x) / (z - y))) n1 (h1 _ n0 n1)
    (h2 _ n0 n1)).comp v (hasDerivAt_s1 hzx hvy)
  have dB := (hQ _ m0 m1).comp v (hasDerivAt_mb (x := x) hvy)
  have := dA.sub dB
  refine this.congr_deriv ?_
  rw [log_s1 hzx hzy hvx hvy, log_mb hvx hvy, Real.log_div hzx hzy]
  unfold fpr
  have e1 : s1 x y z v - 1 = (v - z) * (x - y) / ((v - y) * (z - x)) := by
    rw [s1]; field_simp; ring
  have e2 : mb x y v - 1 = (y - x) / (v - y) := by
    rw [mb]; field_simp; ring
  rw [e1, e2]
  have hyx : y - x ≠ 0 := by intro h; apply hxy; linarith
  field_simp
  ring

lemma dFIV {Q : ℝ → ℝ} (hQ : IsP2 Q) {x y v : ℝ} (hxy : x - y ≠ 0) (hvx : v - x ≠ 0)
    (hvy : v - y ≠ 0) : HasDerivAt (FIV Q x y) (fpr x y y v) v := by
  have m0 := mb_ne_zero hvx hvy
  have m1 := mb_ne_one hxy hvy
  have := ((hQ _ m0 m1).comp v (hasDerivAt_mb (x := x) hvy)).neg
  refine this.congr_deriv ?_
  rw [log_mb hvx hvy]
  unfold fpr
  have e2 : mb x y v - 1 = (y - x) / (v - y) := by
    rw [mb]; field_simp; ring
  rw [e2]
  have hyx : y - x ≠ 0 := by intro h; apply hxy; linarith
  field_simp
  ring

lemma dFIII {Q : ℝ → ℝ} (hQ : IsP2 Q) {x y v : ℝ} (hxy : x - y ≠ 0) (hvx : v - x ≠ 0)
    (hvy : v - y ≠ 0) : HasDerivAt (FIII Q x y) (fpr x y x v) v := by
  have m0 := mb_ne_zero hvx hvy
  have m1 := mb_ne_one hxy hvy
  have dl := ((((hasDerivAt_id v).sub_const x).log hvx).sub
    (((hasDerivAt_id v).sub_const y).log hvy))
  have := ((dl.pow 3).div_const 3).sub ((hQ _ m0 m1).comp v (hasDerivAt_mb (x := x) hvy))
  refine this.congr_deriv ?_
  rw [log_mb hvx hvy]
  unfold fpr
  have e2 : mb x y v - 1 = (y - x) / (v - y) := by
    rw [mb]; field_simp; ring
  rw [e2]
  have hyx : y - x ≠ 0 := by intro h; apply hxy; linarith
  simp only [id, Pi.sub_apply]
  field_simp
  ring

lemma dFL {P1 : ℝ → ℝ} (h1 : IsP1 P1) {x z v : ℝ} (hxz : z - x ≠ 0)
    (hvx : v - x ≠ 0) (hvz : v - z ≠ 0) : HasDerivAt (FL P1 x z) (flin x z v) v := by
  have n0 := sA_ne_zero hxz hvx
  have n1 := sA_ne_one hxz hvz
  have := (hasDerivAt_Psi (k := Real.log (z - x)) n1 (h1 _ n0 n1)).comp v (hasDerivAt_sA x z v)
  refine this.congr_deriv ?_
  rw [sA, Real.log_div hvx hxz]
  unfold flin
  have e : (v - x) / (z - x) - 1 = (v - z) / (z - x) := by field_simp; ring
  rw [e]
  field_simp
  ring

lemma dFH {x v : ℝ} (hvx : v - x ≠ 0) : HasDerivAt (FH x) (flin x x v) v := by
  have := ((((hasDerivAt_id v).sub_const x).log hvx).pow 2).div_const 2
  unfold FH flin
  refine this.congr_deriv ?_
  simp only [id]
  norm_num
  field_simp

lemma dFK {z v : ℝ} (hvz : v - z ≠ 0) : HasDerivAt (FK z) (fcst z v) v := by
  have := ((hasDerivAt_id v).sub_const z).log hvz
  unfold FK fcst
  refine this.congr_deriv ?_
  simp


/-! ## From a primitive to the value of the integral -/

lemma ftc_pos {F f : ℝ → ℝ} {p q : ℝ} (hpq : p ≤ q) (hc : ContinuousOn F (Icc p q))
    (hd : ∀ v ∈ Ioo p q, HasDerivAt F (f v) v) (hf : ∀ v ∈ Ioo p q, 0 ≤ f v) :
    IntervalIntegrable f volume p q ∧ ∫ v in p..q, f v = F q - F p := by
  have hi : IntervalIntegrable f volume p q := by
    apply intervalIntegral.intervalIntegrable_deriv_of_nonneg (g := F)
    · rwa [uIcc_of_le hpq]
    · rw [min_eq_left hpq, max_eq_right hpq]; exact hd
    · rw [min_eq_left hpq, max_eq_right hpq]; exact hf
  exact ⟨hi, intervalIntegral.integral_eq_sub_of_hasDerivAt_of_le hpq hc hd hi⟩

lemma contOn_comp {g σ : ℝ → ℝ} {p q : ℝ} (hσ : ContinuousOn σ (Icc p q))
    (hg : ∀ v ∈ Icc p q, ContinuousAt g (σ v)) : ContinuousOn (fun v => g (σ v)) (Icc p q) :=
  fun v hv => (hg v hv).comp_continuousWithinAt (hσ v hv)

lemma continuous_sA (x z : ℝ) : Continuous (sA x z) :=
  (continuous_id.sub continuous_const).div_const _

lemma contOn_mb {x y p q : ℝ} (hy : ∀ v ∈ Icc p q, v - y ≠ 0) :
    ContinuousOn (mb x y) (Icc p q) := fun v hv =>
  ((continuousAt_id.sub continuousAt_const).div (continuousAt_id.sub continuousAt_const)
    (hy v hv)).continuousWithinAt

lemma contOn_s1 {x y z p q : ℝ} (hzx : z - x ≠ 0) (hy : ∀ v ∈ Icc p q, v - y ≠ 0) :
    ContinuousOn (s1 x y z) (Icc p q) := fun v hv =>
  (((continuousAt_id.sub continuousAt_const).mul continuousAt_const).div
    ((continuousAt_id.sub continuousAt_const).mul continuousAt_const)
    (mul_ne_zero (hy v hv) hzx)).continuousWithinAt

lemma contOn_logsub {x p q : ℝ} (hx : ∀ v ∈ Icc p q, v - x ≠ 0) :
    ContinuousOn (fun v => Real.log (v - x)) (Icc p q) := fun v hv =>
  ((continuousAt_id.sub continuousAt_const).log (hx v hv)).continuousWithinAt

lemma cPhiN {k s : ℝ} (h1 : s ≠ 1) : ContinuousAt (Phi P1n P2a k) s :=
  contAt_Phi h1 (contAt_P1n h1) (contAt_P2a h1)
lemma cPhiM {k s : ℝ} (h1 : s ≠ 1) : ContinuousAt (Phi P1p P2a k) s :=
  contAt_Phi h1 (contAt_P1p s) (contAt_P2a h1)
lemma cPhiG {k s : ℝ} (h0 : s ≠ 0) (h1 : s ≠ 1) : ContinuousAt (Phi P1p P2g k) s :=
  contAt_Phi h1 (contAt_P1p s) (contAt_P2g h0 h1)
lemma cPsiN {k s : ℝ} (h1 : s ≠ 1) : ContinuousAt (Psi P1n k) s :=
  contAt_Psi h1 (contAt_P1n h1)
lemma cPsiM {k s : ℝ} (h1 : s ≠ 1) : ContinuousAt (Psi P1p k) s :=
  contAt_Psi h1 (contAt_P1p s)

lemma intA {P1 P2 : ℝ → ℝ} (h1 : IsP1 P1) (h2 : IsP2 P2) {x z p q : ℝ} (hpq : p ≤ q)
    (hxz : z - x ≠ 0) (hz : ∀ v ∈ Icc p q, 0 < v - z) (hx : ∀ v ∈ Ioo p q, v - x ≠ 0)
    (hc : ∀ v ∈ Icc p q, ContinuousAt (Phi P1 P2 (Real.log (z - x))) (sA x z v)) :
    IntervalIntegrable (fsq x z) volume p q ∧
      ∫ v in p..q, fsq x z v = FA P1 P2 x z q - FA P1 P2 x z p :=
  ftc_pos hpq (contOn_comp (σ := sA x z) (continuous_sA x z).continuousOn hc)
    (fun v hv => dFA h1 h2 hxz (hx v hv) (hz v (Ioo_subset_Icc_self hv)).ne')
    (fun v hv => div_nonneg (sq_nonneg _) (hz v (Ioo_subset_Icc_self hv)).le)

lemma intC {x p q : ℝ} (hpq : p ≤ q) (hx : ∀ v ∈ Icc p q, 0 < v - x) :
    IntervalIntegrable (fsq x x) volume p q ∧
      ∫ v in p..q, fsq x x v = FC x q - FC x p := by
  have hc : ContinuousOn (FC x) (Icc p q) := by
    have := (contOn_logsub (fun v hv => (hx v hv).ne')).pow 3
    exact this.div_const 3
  exact ftc_pos hpq hc (fun v hv => dFC (hx v (Ioo_subset_Icc_self hv)).ne')
    (fun v hv => div_nonneg (sq_nonneg _) (hx v (Ioo_subset_Icc_self hv)).le)

lemma intG {P1 P2 Q : ℝ → ℝ} (h1 : IsP1 P1) (h2 : IsP2 P2) (hQ : IsP2 Q) {x y z p q : ℝ}
    (hpq : p ≤ q) (hxy : x - y ≠ 0) (hzx : z - x ≠ 0) (hzy : z - y ≠ 0)
    (hz : ∀ v ∈ Icc p q, 0 < v - z) (hy : ∀ v ∈ Icc p q, v - y ≠ 0)
    (hx : ∀ v ∈ Ioo p q, v - x ≠ 0)
    (hc1 : ∀ v ∈ Icc p q, ContinuousAt (Phi P1 P2 (Real.log ((z - x) / (z - y)))) (s1 x y z v))
    (hc2 : ∀ v ∈ Icc p q, ContinuousAt Q (mb x y v)) :
    IntervalIntegrable (fpr x y z) volume p q ∧
      ∫ v in p..q, fpr x y z v = FG P1 P2 Q x y z q - FG P1 P2 Q x y z p := by
  have hc : ContinuousOn (FG P1 P2 Q x y z) (Icc p q) :=
    (contOn_comp (contOn_s1 hzx hy) hc1).sub (contOn_comp (contOn_mb hy) hc2)
  exact ftc_pos hpq hc
    (fun v hv => dFG h1 h2 hQ hxy hzx hzy (hx v hv) (hy v (Ioo_subset_Icc_self hv))
      (hz v (Ioo_subset_Icc_self hv)).ne')
    (fun v hv => div_nonneg (sq_nonneg _) (hz v (Ioo_subset_Icc_self hv)).le)

lemma intIV {Q : ℝ → ℝ} (hQ : IsP2 Q) {x y p q : ℝ} (hpq : p ≤ q) (hxy : x - y ≠ 0)
    (hy : ∀ v ∈ Icc p q, 0 < v - y) (hx : ∀ v ∈ Ioo p q, v - x ≠ 0)
    (hc2 : ∀ v ∈ Icc p q, ContinuousAt Q (mb x y v)) :
    IntervalIntegrable (fpr x y y) volume p q ∧
      ∫ v in p..q, fpr x y y v = FIV Q x y q - FIV Q x y p := by
  have hc : ContinuousOn (FIV Q x y) (Icc p q) :=
    (contOn_comp (contOn_mb (fun v hv => (hy v hv).ne')) hc2).neg
  exact ftc_pos hpq hc
    (fun v hv => dFIV hQ hxy (hx v hv) (hy v (Ioo_subset_Icc_self hv)).ne')
    (fun v hv => div_nonneg (sq_nonneg _) (hy v (Ioo_subset_Icc_self hv)).le)

lemma intIII {Q : ℝ → ℝ} (hQ : IsP2 Q) {x y p q : ℝ} (hpq : p ≤ q) (hxy : x - y ≠ 0)
    (hx : ∀ v ∈ Icc p q, 0 < v - x) (hy : ∀ v ∈ Icc p q, v - y ≠ 0)
    (hc2 : ∀ v ∈ Icc p q, ContinuousAt Q (mb x y v)) :
    IntervalIntegrable (fpr x y x) volume p q ∧
      ∫ v in p..q, fpr x y x v = FIII Q x y q - FIII Q x y p := by
  have hc : ContinuousOn (FIII Q x y) (Icc p q) :=
    ((((contOn_logsub (fun v hv => (hx v hv).ne')).sub (contOn_logsub hy)).pow 3).div_const 3).sub
      (contOn_comp (contOn_mb hy) hc2)
  exact ftc_pos hpq hc
    (fun v hv => dFIII hQ hxy (hx v (Ioo_subset_Icc_self hv)).ne' (hy v (Ioo_subset_Icc_self hv)))
    (fun v hv => div_nonneg (sq_nonneg _) (hx v (Ioo_subset_Icc_self hv)).le)

lemma intble_flin {x z p q : ℝ} (hpq : p ≤ q) (hz : ∀ v ∈ Icc p q, 0 < v - z) :
    IntervalIntegrable (flin x z) volume p q := by
  have h1 := (intervalIntegral.intervalIntegrable_log' (a := p - x) (b := q - x)).comp_sub_right x
  simp only [sub_add_cancel] at h1
  have h2 := h1.mul_continuousOn (g := fun v => (v - z)⁻¹) (by
    rw [uIcc_of_le hpq]
    intro v hv
    exact ((continuousAt_id.sub continuousAt_const).inv₀ (hz v hv).ne').continuousWithinAt)
  have e : flin x z = fun v => Real.log (v - x) * (v - z)⁻¹ := by
    funext v; rw [flin, div_eq_mul_inv]
  rw [e]
  exact h2

lemma intL {P1 : ℝ → ℝ} (h1 : IsP1 P1) {x z p q : ℝ} (hpq : p ≤ q)
    (hxz : z - x ≠ 0) (hz : ∀ v ∈ Icc p q, 0 < v - z) (hx : ∀ v ∈ Ioo p q, v - x ≠ 0)
    (hc : ∀ v ∈ Icc p q, ContinuousAt (Psi P1 (Real.log (z - x))) (sA x z v)) :
    IntervalIntegrable (flin x z) volume p q ∧
      ∫ v in p..q, flin x z v = FL P1 x z q - FL P1 x z p := by
  have hi := intble_flin (x := x) hpq hz
  exact ⟨hi, intervalIntegral.integral_eq_sub_of_hasDerivAt_of_le hpq
    (contOn_comp (σ := sA x z) (continuous_sA x z).continuousOn hc)
    (fun v hv => dFL h1 hxz (hx v hv) (hz v (Ioo_subset_Icc_self hv)).ne') hi⟩

lemma intH {x p q : ℝ} (hpq : p ≤ q) (hx : ∀ v ∈ Icc p q, 0 < v - x) :
    IntervalIntegrable (flin x x) volume p q ∧
      ∫ v in p..q, flin x x v = FH x q - FH x p := by
  have hi := intble_flin (x := x) hpq hx
  have hc : ContinuousOn (FH x) (Icc p q) :=
    ((contOn_logsub (fun v hv => (hx v hv).ne')).pow 2).div_const 2
  exact ⟨hi, intervalIntegral.integral_eq_sub_of_hasDerivAt_of_le hpq hc
    (fun v hv => dFH (hx v (Ioo_subset_Icc_self hv)).ne') hi⟩

lemma intK {z p q : ℝ} (hpq : p ≤ q) (hz : ∀ v ∈ Icc p q, 0 < v - z) :
    IntervalIntegrable (fcst z) volume p q ∧
      ∫ v in p..q, fcst z v = FK z q - FK z p := by
  have hc : ContinuousOn (fcst z) (Icc p q) := fun v hv =>
    (continuousAt_const.div (continuousAt_id.sub continuousAt_const)
      (hz v hv).ne').continuousWithinAt
  have hi : IntervalIntegrable (fcst z) volume p q := by
    apply ContinuousOn.intervalIntegrable; rwa [uIcc_of_le hpq]
  have hK : ContinuousOn (FK z) (Icc p q) := contOn_logsub (fun v hv => (hz v hv).ne')
  exact ⟨hi, intervalIntegral.integral_eq_sub_of_hasDerivAt_of_le hpq hK
    (fun v hv => dFK (hz v (Ioo_subset_Icc_self hv)).ne') hi⟩


/-! ## Numerical facts about `√2`, `a`, `b` -/

lemma s2sq : Real.sqrt 2 ^ 2 = 2 := Real.sq_sqrt (by norm_num)
lemma s2pos : 0 < Real.sqrt 2 := Real.sqrt_pos.mpr (by norm_num)
lemma s2lo : (1.414 : ℝ) < Real.sqrt 2 := by
  rw [Real.lt_sqrt (by norm_num)]; norm_num
lemma s2hi : Real.sqrt 2 < 1.415 := by
  rw [Real.sqrt_lt' (by norm_num)]; norm_num
lemma a_pos : 0 < a := by simp only [Y01.a]; linarith [s2lo]
lemma hba : b ≤ a := by simp only [Y01.a, Y01.b]; linarith
lemma hb1 : b ≤ 1 := by simp only [Y01.b]; linarith [s2hi]
lemma h1a : 1 ≤ a := by simp only [Y01.a]; linarith [s2lo]

lemma s2_2 : Real.sqrt 2 ^ 2 = 2 := by
  rw [show Real.sqrt 2 ^ 2 = (Real.sqrt 2 ^ 2) ^ 1 by ring, s2sq]; norm_num
lemma s2_3 : Real.sqrt 2 ^ 3 = 2 * Real.sqrt 2 := by
  rw [show Real.sqrt 2 ^ 3 = (Real.sqrt 2 ^ 2) ^ 1 * Real.sqrt 2 by ring, s2sq]; norm_num
lemma s2_4 : Real.sqrt 2 ^ 4 = 4 := by
  rw [show Real.sqrt 2 ^ 4 = (Real.sqrt 2 ^ 2) ^ 2 by ring, s2sq]; norm_num
lemma s2_5 : Real.sqrt 2 ^ 5 = 4 * Real.sqrt 2 := by
  rw [show Real.sqrt 2 ^ 5 = (Real.sqrt 2 ^ 2) ^ 2 * Real.sqrt 2 by ring, s2sq]; norm_num
lemma s2_6 : Real.sqrt 2 ^ 6 = 8 := by
  rw [show Real.sqrt 2 ^ 6 = (Real.sqrt 2 ^ 2) ^ 3 by ring, s2sq]; norm_num
lemma s2_7 : Real.sqrt 2 ^ 7 = 8 * Real.sqrt 2 := by
  rw [show Real.sqrt 2 ^ 7 = (Real.sqrt 2 ^ 2) ^ 3 * Real.sqrt 2 by ring, s2sq]; norm_num
lemma s2_8 : Real.sqrt 2 ^ 8 = 16 := by
  rw [show Real.sqrt 2 ^ 8 = (Real.sqrt 2 ^ 2) ^ 4 by ring, s2sq]; norm_num
lemma s2_9 : Real.sqrt 2 ^ 9 = 16 * Real.sqrt 2 := by
  rw [show Real.sqrt 2 ^ 9 = (Real.sqrt 2 ^ 2) ^ 4 * Real.sqrt 2 by ring, s2sq]; norm_num
lemma s2_10 : Real.sqrt 2 ^ 10 = 32 := by
  rw [show Real.sqrt 2 ^ 10 = (Real.sqrt 2 ^ 2) ^ 5 by ring, s2sq]; norm_num
lemma s2_11 : Real.sqrt 2 ^ 11 = 32 * Real.sqrt 2 := by
  rw [show Real.sqrt 2 ^ 11 = (Real.sqrt 2 ^ 2) ^ 5 * Real.sqrt 2 by ring, s2sq]; norm_num
lemma s2_12 : Real.sqrt 2 ^ 12 = 64 := by
  rw [show Real.sqrt 2 ^ 12 = (Real.sqrt 2 ^ 2) ^ 6 by ring, s2sq]; norm_num
lemma s2_13 : Real.sqrt 2 ^ 13 = 64 * Real.sqrt 2 := by
  rw [show Real.sqrt 2 ^ 13 = (Real.sqrt 2 ^ 2) ^ 6 * Real.sqrt 2 by ring, s2sq]; norm_num
lemma s2_14 : Real.sqrt 2 ^ 14 = 128 := by
  rw [show Real.sqrt 2 ^ 14 = (Real.sqrt 2 ^ 2) ^ 7 by ring, s2sq]; norm_num
lemma s2_15 : Real.sqrt 2 ^ 15 = 128 * Real.sqrt 2 := by
  rw [show Real.sqrt 2 ^ 15 = (Real.sqrt 2 ^ 2) ^ 7 * Real.sqrt 2 by ring, s2sq]; norm_num
lemma s2_16 : Real.sqrt 2 ^ 16 = 256 := by
  rw [show Real.sqrt 2 ^ 16 = (Real.sqrt 2 ^ 2) ^ 8 by ring, s2sq]; norm_num
lemma s2_17 : Real.sqrt 2 ^ 17 = 256 * Real.sqrt 2 := by
  rw [show Real.sqrt 2 ^ 17 = (Real.sqrt 2 ^ 2) ^ 8 * Real.sqrt 2 by ring, s2sq]; norm_num
lemma s2_18 : Real.sqrt 2 ^ 18 = 512 := by
  rw [show Real.sqrt 2 ^ 18 = (Real.sqrt 2 ^ 2) ^ 9 by ring, s2sq]; norm_num
lemma s2_19 : Real.sqrt 2 ^ 19 = 512 * Real.sqrt 2 := by
  rw [show Real.sqrt 2 ^ 19 = (Real.sqrt 2 ^ 2) ^ 9 * Real.sqrt 2 by ring, s2sq]; norm_num
lemma s2_20 : Real.sqrt 2 ^ 20 = 1024 := by
  rw [show Real.sqrt 2 ^ 20 = (Real.sqrt 2 ^ 2) ^ 10 by ring, s2sq]; norm_num
lemma s2_21 : Real.sqrt 2 ^ 21 = 1024 * Real.sqrt 2 := by
  rw [show Real.sqrt 2 ^ 21 = (Real.sqrt 2 ^ 2) ^ 10 * Real.sqrt 2 by ring, s2sq]; norm_num
lemma s2_22 : Real.sqrt 2 ^ 22 = 2048 := by
  rw [show Real.sqrt 2 ^ 22 = (Real.sqrt 2 ^ 2) ^ 11 by ring, s2sq]; norm_num
lemma s2_23 : Real.sqrt 2 ^ 23 = 2048 * Real.sqrt 2 := by
  rw [show Real.sqrt 2 ^ 23 = (Real.sqrt 2 ^ 2) ^ 11 * Real.sqrt 2 by ring, s2sq]; norm_num
lemma s2_24 : Real.sqrt 2 ^ 24 = 4096 := by
  rw [show Real.sqrt 2 ^ 24 = (Real.sqrt 2 ^ 2) ^ 12 by ring, s2sq]; norm_num
lemma nza : a ≠ 0 := a_pos.ne'
lemma nzb : b ≠ 0 := ne_of_gt (by simp only [Y01.b]; linarith [s2lo])
lemma nzs : Real.sqrt 2 ≠ 0 := s2pos.ne'
lemma nzd0 : (b : ℝ) - (-a) ≠ 0 := ne_of_gt (by (try simp only [Y01.a, Y01.b]); linarith [s2lo, s2hi])
lemma nzd1 : (b : ℝ) - 0 ≠ 0 := ne_of_gt (by (try simp only [Y01.a, Y01.b]); linarith [s2lo, s2hi])
lemma nzd2 : (b : ℝ) - a ≠ 0 := ne_of_lt (by (try simp only [Y01.a, Y01.b]); linarith [s2lo, s2hi])
lemma nzd3 : (b : ℝ) - (-b) ≠ 0 := ne_of_gt (by (try simp only [Y01.a, Y01.b]); linarith [s2lo, s2hi])
lemma nzd4 : (b : ℝ) - 1 ≠ 0 := ne_of_lt (by (try simp only [Y01.a, Y01.b]); linarith [s2lo, s2hi])
lemma nzd5 : ((-a) : ℝ) - b ≠ 0 := ne_of_lt (by (try simp only [Y01.a, Y01.b]); linarith [s2lo, s2hi])
lemma nzd6 : ((-a) : ℝ) - 0 ≠ 0 := ne_of_lt (by (try simp only [Y01.a, Y01.b]); linarith [s2lo, s2hi])
lemma nzd7 : ((-a) : ℝ) - a ≠ 0 := ne_of_lt (by (try simp only [Y01.a, Y01.b]); linarith [s2lo, s2hi])
lemma nzd8 : ((-a) : ℝ) - (-b) ≠ 0 := ne_of_lt (by (try simp only [Y01.a, Y01.b]); linarith [s2lo, s2hi])
lemma nzd9 : ((-a) : ℝ) - 1 ≠ 0 := ne_of_lt (by (try simp only [Y01.a, Y01.b]); linarith [s2lo, s2hi])
lemma nzd10 : (0 : ℝ) - b ≠ 0 := ne_of_lt (by (try simp only [Y01.a, Y01.b]); linarith [s2lo, s2hi])
lemma nzd11 : (0 : ℝ) - (-a) ≠ 0 := ne_of_gt (by (try simp only [Y01.a, Y01.b]); linarith [s2lo, s2hi])
lemma nzd12 : (0 : ℝ) - a ≠ 0 := ne_of_lt (by (try simp only [Y01.a, Y01.b]); linarith [s2lo, s2hi])
lemma nzd13 : (0 : ℝ) - (-b) ≠ 0 := ne_of_gt (by (try simp only [Y01.a, Y01.b]); linarith [s2lo, s2hi])
lemma nzd14 : (0 : ℝ) - 1 ≠ 0 := ne_of_lt (by (try simp only [Y01.a, Y01.b]); linarith [s2lo, s2hi])
lemma nzd15 : (a : ℝ) - b ≠ 0 := ne_of_gt (by (try simp only [Y01.a, Y01.b]); linarith [s2lo, s2hi])
lemma nzd16 : (a : ℝ) - (-a) ≠ 0 := ne_of_gt (by (try simp only [Y01.a, Y01.b]); linarith [s2lo, s2hi])
lemma nzd17 : (a : ℝ) - 0 ≠ 0 := ne_of_gt (by (try simp only [Y01.a, Y01.b]); linarith [s2lo, s2hi])
lemma nzd18 : (a : ℝ) - (-b) ≠ 0 := ne_of_gt (by (try simp only [Y01.a, Y01.b]); linarith [s2lo, s2hi])
lemma nzd19 : (a : ℝ) - 1 ≠ 0 := ne_of_gt (by (try simp only [Y01.a, Y01.b]); linarith [s2lo, s2hi])
lemma nzd20 : ((-b) : ℝ) - b ≠ 0 := ne_of_lt (by (try simp only [Y01.a, Y01.b]); linarith [s2lo, s2hi])
lemma nzd21 : ((-b) : ℝ) - (-a) ≠ 0 := ne_of_gt (by (try simp only [Y01.a, Y01.b]); linarith [s2lo, s2hi])
lemma nzd22 : ((-b) : ℝ) - 0 ≠ 0 := ne_of_lt (by (try simp only [Y01.a, Y01.b]); linarith [s2lo, s2hi])
lemma nzd23 : ((-b) : ℝ) - a ≠ 0 := ne_of_lt (by (try simp only [Y01.a, Y01.b]); linarith [s2lo, s2hi])
lemma nzd24 : ((-b) : ℝ) - 1 ≠ 0 := ne_of_lt (by (try simp only [Y01.a, Y01.b]); linarith [s2lo, s2hi])
lemma nzd25 : (1 : ℝ) - b ≠ 0 := ne_of_gt (by (try simp only [Y01.a, Y01.b]); linarith [s2lo, s2hi])
lemma nzd26 : (1 : ℝ) - (-a) ≠ 0 := ne_of_gt (by (try simp only [Y01.a, Y01.b]); linarith [s2lo, s2hi])
lemma nzd27 : (1 : ℝ) - 0 ≠ 0 := ne_of_gt (by (try simp only [Y01.a, Y01.b]); linarith [s2lo, s2hi])
lemma nzd28 : (1 : ℝ) - a ≠ 0 := ne_of_lt (by (try simp only [Y01.a, Y01.b]); linarith [s2lo, s2hi])
lemma nzd29 : (1 : ℝ) - (-b) ≠ 0 := ne_of_gt (by (try simp only [Y01.a, Y01.b]); linarith [s2lo, s2hi])

lemma logP {u : ℝ} (j1 m1 j2 m2 : ℕ) (c1 c2 : ℝ) (hc1 : c1 = ((j1 : ℝ) - j2) / 2)
    (hc2 : c2 = (m1 : ℝ) - m2) (h : u * (Real.sqrt 2 ^ j2 * a ^ m2) = Real.sqrt 2 ^ j1 * a ^ m1) :
    Real.log u = c1 * lam + c2 * L := by
  have hD : 0 < Real.sqrt 2 ^ j2 * a ^ m2 := mul_pos (pow_pos s2pos _) (pow_pos a_pos _)
  have hN : 0 < Real.sqrt 2 ^ j1 * a ^ m1 := mul_pos (pow_pos s2pos _) (pow_pos a_pos _)
  have hu : u = (Real.sqrt 2 ^ j1 * a ^ m1) / (Real.sqrt 2 ^ j2 * a ^ m2) := by
    rw [eq_div_iff hD.ne']; exact h
  rw [hu, Real.log_div hN.ne' hD.ne', Real.log_mul (pow_pos s2pos _).ne' (pow_pos a_pos _).ne',
    Real.log_mul (pow_pos s2pos _).ne' (pow_pos a_pos _).ne', Real.log_pow, Real.log_pow,
    Real.log_pow, Real.log_pow, Real.log_sqrt (by norm_num), hc1, hc2]
  unfold lam L
  ring

lemma logN {u : ℝ} (j1 m1 j2 m2 : ℕ) (c1 c2 : ℝ) (hc1 : c1 = ((j1 : ℝ) - j2) / 2)
    (hc2 : c2 = (m1 : ℝ) - m2) (h : -u * (Real.sqrt 2 ^ j2 * a ^ m2) = Real.sqrt 2 ^ j1 * a ^ m1) :
    Real.log u = c1 * lam + c2 * L := by
  rw [← Real.log_neg_eq_log]; exact logP j1 m1 j2 m2 c1 c2 hc1 hc2 h

lemma blk0 : IntervalIntegrable (fcst 0) volume b a ∧
    ∫ v in b..a, fcst 0 v = FK 0 a - FK 0 b := by
  have hz : ∀ v ∈ Icc b a, 0 < v - 0 := fun v hv => by
    have h1 := hv.1; have h2 := hv.2; simp only [Y01.a, Y01.b] at *; linarith [s2lo, s2hi]
  exact intK hba hz

lemma val0q : FK 0 a = (1 : ℝ) * L := by
  have l0 : Real.log (a - 0) = (0 : ℝ) * lam + (1 : ℝ) * L :=
    logP 0 1 0 0 _ _ (by norm_num) (by norm_num)
      (by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf))
  simp only [FK, l0]
  try simp only [Real.log_zero, sub_zero, Real.log_one, D00.Li2_zero', Li3_zero, D00.Li2_one']
  ring

lemma val0p : FK 0 b = (-1 : ℝ) * L := by
  have l0 : Real.log (b - 0) = (0 : ℝ) * lam + (-1 : ℝ) * L :=
    logP 0 0 0 1 _ _ (by norm_num) (by norm_num)
      (by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf))
  simp only [FK, l0]
  try simp only [Real.log_zero, sub_zero, Real.log_one, D00.Li2_zero', Li3_zero, D00.Li2_one']
  ring

lemma blk1 : IntervalIntegrable (flin b 0) volume b a ∧
    ∫ v in b..a, flin b 0 v = FL P1n b 0 a - FL P1n b 0 b := by
  have hxz : (0 : ℝ) - b ≠ 0 := ne_of_lt (by (try simp only [Y01.a, Y01.b]); linarith [s2lo, s2hi])
  have hz : ∀ v ∈ Icc b a, 0 < v - 0 := fun v hv => by
    have h1 := hv.1; have h2 := hv.2; simp only [Y01.a, Y01.b] at *; linarith [s2lo, s2hi]
  have hx : ∀ v ∈ Ioo b a, v - b ≠ 0 := fun v hv => ne_of_gt (by
    have h1 := hv.1; have h2 := hv.2; simp only [Y01.a, Y01.b] at *; linarith [s2lo, s2hi])
  exact intL isP1n hba hxz hz hx (fun v hv => cPsiN (sA_ne_one hxz (hz v hv).ne'))

lemma val1q : FL P1n b 0 a = (2 : ℝ) * L * lam + (1 : ℝ) * D00.Li2 (-(2 * a)) := by
  have e0 : sA b 0 a = -(2 * a) := by unfold sA; have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf)
  have l0 : Real.log (0 - b) = (0 : ℝ) * lam + (-1 : ℝ) * L :=
    logN 0 0 0 1 _ _ (by norm_num) (by norm_num)
      (by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf))
  have l1 : Real.log (-(2 * a)) = (1 : ℝ) * lam + (1 : ℝ) * L :=
    logN 2 1 0 0 _ _ (by norm_num) (by norm_num)
      (by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf))
  have l2 : Real.log (1 - (-(2 * a))) = (0 : ℝ) * lam + (2 : ℝ) * L :=
    logP 0 2 0 0 _ _ (by norm_num) (by norm_num)
      (by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf))
  simp only [FL, e0, l0, Psi, P1n, l1, l2]
  try simp only [Real.log_zero, sub_zero, Real.log_one, D00.Li2_zero', Li3_zero, D00.Li2_one']
  ring

lemma val1p : FL P1n b 0 b = 0 := by
  have e0 : sA b 0 b = 0 := by unfold sA; have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf)
  have l0 : Real.log (0 - b) = (0 : ℝ) * lam + (-1 : ℝ) * L :=
    logN 0 0 0 1 _ _ (by norm_num) (by norm_num)
      (by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf))
  simp only [FL, e0, l0, Psi, P1n]
  try simp only [Real.log_zero, sub_zero, Real.log_one, D00.Li2_zero', Li3_zero, D00.Li2_one']
  ring

lemma blk2 : IntervalIntegrable (flin (-a) 0) volume b a ∧
    ∫ v in b..a, flin (-a) 0 v = FL P1p (-a) 0 a - FL P1p (-a) 0 b := by
  have hxz : (0 : ℝ) - (-a) ≠ 0 := ne_of_gt (by (try simp only [Y01.a, Y01.b]); linarith [s2lo, s2hi])
  have hz : ∀ v ∈ Icc b a, 0 < v - 0 := fun v hv => by
    have h1 := hv.1; have h2 := hv.2; simp only [Y01.a, Y01.b] at *; linarith [s2lo, s2hi]
  have hx : ∀ v ∈ Ioo b a, v - (-a) ≠ 0 := fun v hv => ne_of_gt (by
    have h1 := hv.1; have h2 := hv.2; simp only [Y01.a, Y01.b] at *; linarith [s2lo, s2hi])
  exact intL isP1p hba hxz hz hx (fun v hv => cPsiM (sA_ne_one hxz (hz v hv).ne'))

lemma val2q : FL P1p (-a) 0 a = (-1 : ℝ) * D00.Li2 (-1) := by
  have e0 : sA (-a) 0 a = (Real.sqrt 2 ^ 2 * a ^ 0 / (Real.sqrt 2 ^ 0 * a ^ 0)) := by unfold sA; have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf)
  have l0 : Real.log (0 - (-a)) = (0 : ℝ) * lam + (1 : ℝ) * L :=
    logP 0 1 0 0 _ _ (by norm_num) (by norm_num)
      (by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf))
  have e1 : (1:ℝ) - ((Real.sqrt 2 ^ 2 * a ^ 0 / (Real.sqrt 2 ^ 0 * a ^ 0))) = -1 := by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf)
  have l1 : Real.log (-1) = (0 : ℝ) * lam + (0 : ℝ) * L :=
    logN 0 0 0 0 _ _ (by norm_num) (by norm_num)
      (by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf))
  simp only [FL, e0, l0, Psi, P1p, e1, l1]
  try simp only [Real.log_zero, sub_zero, Real.log_one, D00.Li2_zero', Li3_zero, D00.Li2_one']
  ring

lemma val2p : FL P1p (-a) 0 b = (-2 : ℝ) * L ^ 2 + (-1 : ℝ) * D00.Li2 (-(b ^ 2)) := by
  have e0 : sA (-a) 0 b = (Real.sqrt 2 ^ 3 * a ^ 0 / (Real.sqrt 2 ^ 0 * a ^ 1)) := by unfold sA; have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf)
  have l0 : Real.log (0 - (-a)) = (0 : ℝ) * lam + (1 : ℝ) * L :=
    logP 0 1 0 0 _ _ (by norm_num) (by norm_num)
      (by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf))
  have e1 : (1:ℝ) - ((Real.sqrt 2 ^ 3 * a ^ 0 / (Real.sqrt 2 ^ 0 * a ^ 1))) = -(b ^ 2) := by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf)
  have l1 : Real.log (-(b ^ 2)) = (0 : ℝ) * lam + (-2 : ℝ) * L :=
    logN 0 0 0 2 _ _ (by norm_num) (by norm_num)
      (by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf))
  simp only [FL, e0, l0, Psi, P1p, e1, l1]
  try simp only [Real.log_zero, sub_zero, Real.log_one, D00.Li2_zero', Li3_zero, D00.Li2_one']
  ring

lemma blk3 : IntervalIntegrable (flin 0 0) volume b a ∧
    ∫ v in b..a, flin 0 0 v = FH 0 a - FH 0 b := by
  have hx : ∀ v ∈ Icc b a, 0 < v - 0 := fun v hv => by
    have h1 := hv.1; have h2 := hv.2; simp only [Y01.a, Y01.b] at *; linarith [s2lo, s2hi]
  exact intH hba hx

lemma val3q : FH 0 a = (1 / 2 : ℝ) * L ^ 2 := by
  have l0 : Real.log (a - 0) = (0 : ℝ) * lam + (1 : ℝ) * L :=
    logP 0 1 0 0 _ _ (by norm_num) (by norm_num)
      (by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf))
  simp only [FH, l0]
  try simp only [Real.log_zero, sub_zero, Real.log_one, D00.Li2_zero', Li3_zero, D00.Li2_one']
  ring

lemma val3p : FH 0 b = (1 / 2 : ℝ) * L ^ 2 := by
  have l0 : Real.log (b - 0) = (0 : ℝ) * lam + (-1 : ℝ) * L :=
    logP 0 0 0 1 _ _ (by norm_num) (by norm_num)
      (by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf))
  simp only [FH, l0]
  try simp only [Real.log_zero, sub_zero, Real.log_one, D00.Li2_zero', Li3_zero, D00.Li2_one']
  ring

lemma blk4 : IntervalIntegrable (fsq b 0) volume b a ∧
    ∫ v in b..a, fsq b 0 v = FA P1n P2a b 0 a - FA P1n P2a b 0 b := by
  have hxz : (0 : ℝ) - b ≠ 0 := ne_of_lt (by (try simp only [Y01.a, Y01.b]); linarith [s2lo, s2hi])
  have hz : ∀ v ∈ Icc b a, 0 < v - 0 := fun v hv => by
    have h1 := hv.1; have h2 := hv.2; simp only [Y01.a, Y01.b] at *; linarith [s2lo, s2hi]
  have hx : ∀ v ∈ Ioo b a, v - b ≠ 0 := fun v hv => ne_of_gt (by
    have h1 := hv.1; have h2 := hv.2; simp only [Y01.a, Y01.b] at *; linarith [s2lo, s2hi])
  exact intA isP1n isP2a hba hxz hz hx (fun v hv => cPhiN (sA_ne_one hxz (hz v hv).ne'))

lemma val4q : FA P1n P2a b 0 a = (2 : ℝ) * L * lam ^ 2 + (2 : ℝ) * D00.Li2 (-(2 * a)) * lam + (-2 : ℝ) * Li3 (-(2 * a)) := by
  have e0 : sA b 0 a = -(2 * a) := by unfold sA; have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf)
  have l0 : Real.log (0 - b) = (0 : ℝ) * lam + (-1 : ℝ) * L :=
    logN 0 0 0 1 _ _ (by norm_num) (by norm_num)
      (by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf))
  have l1 : Real.log (-(2 * a)) = (1 : ℝ) * lam + (1 : ℝ) * L :=
    logN 2 1 0 0 _ _ (by norm_num) (by norm_num)
      (by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf))
  have l2 : Real.log (1 - (-(2 * a))) = (0 : ℝ) * lam + (2 : ℝ) * L :=
    logP 0 2 0 0 _ _ (by norm_num) (by norm_num)
      (by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf))
  simp only [FA, e0, l0, Phi, P1n, P2a, l1, l2]
  try simp only [Real.log_zero, sub_zero, Real.log_one, D00.Li2_zero', Li3_zero, D00.Li2_one']
  ring

lemma val4p : FA P1n P2a b 0 b = 0 := by
  have e0 : sA b 0 b = 0 := by unfold sA; have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf)
  have l0 : Real.log (0 - b) = (0 : ℝ) * lam + (-1 : ℝ) * L :=
    logN 0 0 0 1 _ _ (by norm_num) (by norm_num)
      (by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf))
  simp only [FA, e0, l0, Phi, P1n, P2a]
  try simp only [Real.log_zero, sub_zero, Real.log_one, D00.Li2_zero', Li3_zero, D00.Li2_one']
  ring

lemma blk5 : IntervalIntegrable (fsq (-a) 0) volume b a ∧
    ∫ v in b..a, fsq (-a) 0 v = FA P1p P2g (-a) 0 a - FA P1p P2g (-a) 0 b := by
  have hxz : (0 : ℝ) - (-a) ≠ 0 := ne_of_gt (by (try simp only [Y01.a, Y01.b]); linarith [s2lo, s2hi])
  have hz : ∀ v ∈ Icc b a, 0 < v - 0 := fun v hv => by
    have h1 := hv.1; have h2 := hv.2; simp only [Y01.a, Y01.b] at *; linarith [s2lo, s2hi]
  have hx : ∀ v ∈ Ioo b a, v - (-a) ≠ 0 := fun v hv => ne_of_gt (by
    have h1 := hv.1; have h2 := hv.2; simp only [Y01.a, Y01.b] at *; linarith [s2lo, s2hi])
  have hx2 : ∀ v ∈ Icc b a, v - (-a) ≠ 0 := fun v hv => ne_of_gt (by
    have h1 := hv.1; have h2 := hv.2; simp only [Y01.a, Y01.b] at *; linarith [s2lo, s2hi])
  exact intA isP1p isP2g hba hxz hz hx (fun v hv => cPhiG (sA_ne_zero hxz (hx2 v hv)) (sA_ne_one hxz (hz v hv).ne'))

lemma val5q : FA P1p P2g (-a) 0 a = (-2 : ℝ) * L * D00.Li2 (-1) + (-2 : ℝ) * D00.Li2 (1 / 2) * lam + (-2 : ℝ) * Li3 (1 / 2) + (-2 / 3 : ℝ) * lam ^ 3 := by
  have e0 : sA (-a) 0 a = (Real.sqrt 2 ^ 2 * a ^ 0 / (Real.sqrt 2 ^ 0 * a ^ 0)) := by unfold sA; have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf)
  have l0 : Real.log (0 - (-a)) = (0 : ℝ) * lam + (1 : ℝ) * L :=
    logP 0 1 0 0 _ _ (by norm_num) (by norm_num)
      (by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf))
  have e1 : (1:ℝ) - ((Real.sqrt 2 ^ 2 * a ^ 0 / (Real.sqrt 2 ^ 0 * a ^ 0))) = -1 := by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf)
  have e2 : ((Real.sqrt 2 ^ 2 * a ^ 0 / (Real.sqrt 2 ^ 0 * a ^ 0)))⁻¹ = 1 / 2 := by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf)
  have l1 : Real.log (-1) = (0 : ℝ) * lam + (0 : ℝ) * L :=
    logN 0 0 0 0 _ _ (by norm_num) (by norm_num)
      (by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf))
  have l2 : Real.log (1 / 2) = (-1 : ℝ) * lam + (0 : ℝ) * L :=
    logP 0 0 2 0 _ _ (by norm_num) (by norm_num)
      (by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf))
  have l3 : Real.log (1 - (1 / 2)) = (-1 : ℝ) * lam + (0 : ℝ) * L :=
    logP 0 0 2 0 _ _ (by norm_num) (by norm_num)
      (by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf))
  simp only [FA, e0, l0, Phi, P1p, P2g, e1, e2, E3, l1, l2, l3]
  try simp only [Real.log_zero, sub_zero, Real.log_one, D00.Li2_zero', Li3_zero, D00.Li2_one']
  ring

lemma val5p : FA P1p P2g (-a) 0 b = (-10 / 3 : ℝ) * L ^ 3 + (3 : ℝ) * L ^ 2 * lam + (-2 : ℝ) * L * D00.Li2 (-(b ^ 2)) + (2 : ℝ) * L * D00.Li2 ((2 + Real.sqrt 2) / 4) + (-3 : ℝ) * D00.Li2 ((2 + Real.sqrt 2) / 4) * lam + (-2 : ℝ) * Li3 ((2 + Real.sqrt 2) / 4) + (-9 / 4 : ℝ) * lam ^ 3 := by
  have e0 : sA (-a) 0 b = (Real.sqrt 2 ^ 3 * a ^ 0 / (Real.sqrt 2 ^ 0 * a ^ 1)) := by unfold sA; have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf)
  have l0 : Real.log (0 - (-a)) = (0 : ℝ) * lam + (1 : ℝ) * L :=
    logP 0 1 0 0 _ _ (by norm_num) (by norm_num)
      (by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf))
  have e1 : (1:ℝ) - ((Real.sqrt 2 ^ 3 * a ^ 0 / (Real.sqrt 2 ^ 0 * a ^ 1))) = -(b ^ 2) := by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf)
  have e2 : ((Real.sqrt 2 ^ 3 * a ^ 0 / (Real.sqrt 2 ^ 0 * a ^ 1)))⁻¹ = (2 + Real.sqrt 2) / 4 := by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf)
  have l1 : Real.log (-(b ^ 2)) = (0 : ℝ) * lam + (-2 : ℝ) * L :=
    logN 0 0 0 2 _ _ (by norm_num) (by norm_num)
      (by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf))
  have l2 : Real.log ((2 + Real.sqrt 2) / 4) = (-3 / 2 : ℝ) * lam + (1 : ℝ) * L :=
    logP 0 1 3 0 _ _ (by norm_num) (by norm_num)
      (by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf))
  have l3 : Real.log (1 - ((2 + Real.sqrt 2) / 4)) = (-3 / 2 : ℝ) * lam + (-1 : ℝ) * L :=
    logP 0 0 3 1 _ _ (by norm_num) (by norm_num)
      (by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf))
  simp only [FA, e0, l0, Phi, P1p, P2g, e1, e2, E3, l1, l2, l3]
  try simp only [Real.log_zero, sub_zero, Real.log_one, D00.Li2_zero', Li3_zero, D00.Li2_one']
  ring

lemma blk6 : IntervalIntegrable (fpr b (-a) 0) volume b 1 ∧
    ∫ v in b..1, fpr b (-a) 0 v = FG P1n P2a P2a b (-a) 0 1 - FG P1n P2a P2a b (-a) 0 b := by
  have hxy : (b : ℝ) - (-a) ≠ 0 := ne_of_gt (by (try simp only [Y01.a, Y01.b]); linarith [s2lo, s2hi])
  have hzx : (0 : ℝ) - b ≠ 0 := ne_of_lt (by (try simp only [Y01.a, Y01.b]); linarith [s2lo, s2hi])
  have hzy : (0 : ℝ) - (-a) ≠ 0 := ne_of_gt (by (try simp only [Y01.a, Y01.b]); linarith [s2lo, s2hi])
  have hz : ∀ v ∈ Icc b 1, 0 < v - 0 := fun v hv => by
    have h1 := hv.1; have h2 := hv.2; simp only [Y01.a, Y01.b] at *; linarith [s2lo, s2hi]
  have hy : ∀ v ∈ Icc b 1, v - (-a) ≠ 0 := fun v hv => ne_of_gt (by
    have h1 := hv.1; have h2 := hv.2; simp only [Y01.a, Y01.b] at *; linarith [s2lo, s2hi])
  have hx : ∀ v ∈ Ioo b 1, v - b ≠ 0 := fun v hv => ne_of_gt (by
    have h1 := hv.1; have h2 := hv.2; simp only [Y01.a, Y01.b] at *; linarith [s2lo, s2hi])
  exact intG isP1n isP2a isP2a hb1 hxy hzx hzy hz hy hx
    (fun v hv => cPhiN (s1_ne_one hxy hzx (hy v hv) (hz v hv).ne'))
    (fun v hv => contAt_P2a (mb_ne_one hxy (hy v hv)))

lemma val6q : FG P1n P2a P2a b (-a) 0 1 = (4 : ℝ) * L ^ 3 + (-4 : ℝ) * L * D00.Li2 (-1) + (4 : ℝ) * L * D00.Li2 (b ^ 2) + (-2 : ℝ) * Li3 (-1) + (2 : ℝ) * Li3 (b ^ 2) := by
  have e0 : mb b (-a) 1 = b ^ 2 := by unfold mb; have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf)
  have e1 : s1 b (-a) 0 1 = -1 := by unfold s1; have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf)
  have l0 : Real.log ((0 - b) / (0 - (-a))) = (0 : ℝ) * lam + (-2 : ℝ) * L :=
    logN 0 0 0 2 _ _ (by norm_num) (by norm_num)
      (by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf))
  have l1 : Real.log (-1) = (0 : ℝ) * lam + (0 : ℝ) * L :=
    logN 0 0 0 0 _ _ (by norm_num) (by norm_num)
      (by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf))
  have l2 : Real.log (1 - (-1)) = (1 : ℝ) * lam + (0 : ℝ) * L :=
    logP 2 0 0 0 _ _ (by norm_num) (by norm_num)
      (by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf))
  have l3 : Real.log (b ^ 2) = (0 : ℝ) * lam + (-2 : ℝ) * L :=
    logP 0 0 0 2 _ _ (by norm_num) (by norm_num)
      (by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf))
  have l4 : Real.log (1 - (b ^ 2)) = (1 : ℝ) * lam + (-1 : ℝ) * L :=
    logP 2 0 0 1 _ _ (by norm_num) (by norm_num)
      (by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf))
  simp only [FG, e0, e1, l0, Phi, P1n, P2a, l1, l2, l3, l4]
  try simp only [Real.log_zero, sub_zero, Real.log_one, D00.Li2_zero', Li3_zero, D00.Li2_one']
  ring

lemma val6p : FG P1n P2a P2a b (-a) 0 b = 0 := by
  have e0 : mb b (-a) b = 0 := by unfold mb; have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf)
  have e1 : s1 b (-a) 0 b = 0 := by unfold s1; have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf)
  have l0 : Real.log ((0 - b) / (0 - (-a))) = (0 : ℝ) * lam + (-2 : ℝ) * L :=
    logN 0 0 0 2 _ _ (by norm_num) (by norm_num)
      (by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf))
  simp only [FG, e0, e1, l0, Phi, P1n, P2a]
  try simp only [Real.log_zero, sub_zero, Real.log_one, D00.Li2_zero', Li3_zero, D00.Li2_one']
  ring

lemma blk7 : IntervalIntegrable (fsq 0 0) volume b a ∧
    ∫ v in b..a, fsq 0 0 v = FC 0 a - FC 0 b := by
  have hx : ∀ v ∈ Icc b a, 0 < v - 0 := fun v hv => by
    have h1 := hv.1; have h2 := hv.2; simp only [Y01.a, Y01.b] at *; linarith [s2lo, s2hi]
  exact intC hba hx

lemma val7q : FC 0 a = (1 / 3 : ℝ) * L ^ 3 := by
  have l0 : Real.log (a - 0) = (0 : ℝ) * lam + (1 : ℝ) * L :=
    logP 0 1 0 0 _ _ (by norm_num) (by norm_num)
      (by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf))
  simp only [FC, l0]
  try simp only [Real.log_zero, sub_zero, Real.log_one, D00.Li2_zero', Li3_zero, D00.Li2_one']
  ring

lemma val7p : FC 0 b = (-1 / 3 : ℝ) * L ^ 3 := by
  have l0 : Real.log (b - 0) = (0 : ℝ) * lam + (-1 : ℝ) * L :=
    logP 0 0 0 1 _ _ (by norm_num) (by norm_num)
      (by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf))
  simp only [FC, l0]
  try simp only [Real.log_zero, sub_zero, Real.log_one, D00.Li2_zero', Li3_zero, D00.Li2_one']
  ring

lemma blk8 : IntervalIntegrable (fpr b 0 0) volume b 1 ∧
    ∫ v in b..1, fpr b 0 0 v = FIV P2a b 0 1 - FIV P2a b 0 b := by
  have hxy : (b : ℝ) - 0 ≠ 0 := ne_of_gt (by (try simp only [Y01.a, Y01.b]); linarith [s2lo, s2hi])
  have hy : ∀ v ∈ Icc b 1, 0 < v - 0 := fun v hv => by
    have h1 := hv.1; have h2 := hv.2; simp only [Y01.a, Y01.b] at *; linarith [s2lo, s2hi]
  have hx : ∀ v ∈ Ioo b 1, v - b ≠ 0 := fun v hv => ne_of_gt (by
    have h1 := hv.1; have h2 := hv.2; simp only [Y01.a, Y01.b] at *; linarith [s2lo, s2hi])
  exact intIV isP2a hb1 hxy hy hx (fun v hv => contAt_P2a (mb_ne_one hxy (hy v hv).ne'))

lemma val8q : FIV P2a b 0 1 = (1 : ℝ) * L ^ 3 + (-1 : ℝ) * L ^ 2 * lam + (2 : ℝ) * L * D00.Li2 (2 - Real.sqrt 2) + (1 / 4 : ℝ) * L * lam ^ 2 + (-1 : ℝ) * D00.Li2 (2 - Real.sqrt 2) * lam + (2 : ℝ) * Li3 (2 - Real.sqrt 2) := by
  have e0 : mb b 0 1 = 2 - Real.sqrt 2 := by unfold mb; have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf)
  have l0 : Real.log (2 - Real.sqrt 2) = (1 / 2 : ℝ) * lam + (-1 : ℝ) * L :=
    logP 1 0 0 1 _ _ (by norm_num) (by norm_num)
      (by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf))
  have l1 : Real.log (1 - (2 - Real.sqrt 2)) = (0 : ℝ) * lam + (-1 : ℝ) * L :=
    logP 0 0 0 1 _ _ (by norm_num) (by norm_num)
      (by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf))
  simp only [FIV, e0, P2a, l0, l1]
  try simp only [Real.log_zero, sub_zero, Real.log_one, D00.Li2_zero', Li3_zero, D00.Li2_one']
  ring

lemma val8p : FIV P2a b 0 b = 0 := by
  have e0 : mb b 0 b = 0 := by unfold mb; have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf)
  simp only [FIV, e0, P2a]
  try simp only [Real.log_zero, sub_zero, Real.log_one, D00.Li2_zero', Li3_zero, D00.Li2_one']
  ring

lemma blk9 : IntervalIntegrable (fpr (-a) 0 0) volume b a ∧
    ∫ v in b..a, fpr (-a) 0 0 v = FIV P2g (-a) 0 a - FIV P2g (-a) 0 b := by
  have hxy : ((-a) : ℝ) - 0 ≠ 0 := ne_of_lt (by (try simp only [Y01.a, Y01.b]); linarith [s2lo, s2hi])
  have hy : ∀ v ∈ Icc b a, 0 < v - 0 := fun v hv => by
    have h1 := hv.1; have h2 := hv.2; simp only [Y01.a, Y01.b] at *; linarith [s2lo, s2hi]
  have hx : ∀ v ∈ Ioo b a, v - (-a) ≠ 0 := fun v hv => ne_of_gt (by
    have h1 := hv.1; have h2 := hv.2; simp only [Y01.a, Y01.b] at *; linarith [s2lo, s2hi])
  have hx2 : ∀ v ∈ Icc b a, v - (-a) ≠ 0 := fun v hv => ne_of_gt (by
    have h1 := hv.1; have h2 := hv.2; simp only [Y01.a, Y01.b] at *; linarith [s2lo, s2hi])
  exact intIV isP2g hba hxy hy hx (fun v hv => contAt_P2g (mb_ne_zero (hx2 v hv) (hy v hv).ne') (mb_ne_one hxy (hy v hv).ne'))

lemma val9q : FIV P2g (-a) 0 a = (2 : ℝ) * D00.Li2 (1 / 2) * lam + (2 : ℝ) * Li3 (1 / 2) + (2 / 3 : ℝ) * lam ^ 3 := by
  have e0 : mb (-a) 0 a = (Real.sqrt 2 ^ 2 * a ^ 0 / (Real.sqrt 2 ^ 0 * a ^ 0)) := by unfold mb; have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf)
  have e1 : ((Real.sqrt 2 ^ 2 * a ^ 0 / (Real.sqrt 2 ^ 0 * a ^ 0)))⁻¹ = 1 / 2 := by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf)
  have l0 : Real.log (1 / 2) = (-1 : ℝ) * lam + (0 : ℝ) * L :=
    logP 0 0 2 0 _ _ (by norm_num) (by norm_num)
      (by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf))
  have l1 : Real.log (1 - (1 / 2)) = (-1 : ℝ) * lam + (0 : ℝ) * L :=
    logP 0 0 2 0 _ _ (by norm_num) (by norm_num)
      (by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf))
  simp only [FIV, e0, P2g, E3, e1, l0, l1]
  try simp only [Real.log_zero, sub_zero, Real.log_one, D00.Li2_zero', Li3_zero, D00.Li2_one']
  ring

lemma val9p : FIV P2g (-a) 0 b = (-4 / 3 : ℝ) * L ^ 3 + (-3 : ℝ) * L ^ 2 * lam + (2 : ℝ) * L * D00.Li2 ((2 - Real.sqrt 2) / 4) + (3 : ℝ) * D00.Li2 ((2 - Real.sqrt 2) / 4) * lam + (2 : ℝ) * Li3 ((2 - Real.sqrt 2) / 4) + (9 / 4 : ℝ) * lam ^ 3 := by
  have e0 : mb (-a) 0 b = (Real.sqrt 2 ^ 3 * a ^ 1 / (Real.sqrt 2 ^ 0 * a ^ 0)) := by unfold mb; have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf)
  have e1 : ((Real.sqrt 2 ^ 3 * a ^ 1 / (Real.sqrt 2 ^ 0 * a ^ 0)))⁻¹ = (2 - Real.sqrt 2) / 4 := by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf)
  have l0 : Real.log ((2 - Real.sqrt 2) / 4) = (-3 / 2 : ℝ) * lam + (-1 : ℝ) * L :=
    logP 0 0 3 1 _ _ (by norm_num) (by norm_num)
      (by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf))
  have l1 : Real.log (1 - ((2 - Real.sqrt 2) / 4)) = (-3 / 2 : ℝ) * lam + (1 : ℝ) * L :=
    logP 0 1 3 0 _ _ (by norm_num) (by norm_num)
      (by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf))
  simp only [FIV, e0, P2g, E3, e1, l0, l1]
  try simp only [Real.log_zero, sub_zero, Real.log_one, D00.Li2_zero', Li3_zero, D00.Li2_one']
  ring

lemma blk10 : IntervalIntegrable (fpr 0 (-a) 0) volume b a ∧
    ∫ v in b..a, fpr 0 (-a) 0 v = FIII P2a 0 (-a) a - FIII P2a 0 (-a) b := by
  have hxy : (0 : ℝ) - (-a) ≠ 0 := ne_of_gt (by (try simp only [Y01.a, Y01.b]); linarith [s2lo, s2hi])
  have hx : ∀ v ∈ Icc b a, 0 < v - 0 := fun v hv => by
    have h1 := hv.1; have h2 := hv.2; simp only [Y01.a, Y01.b] at *; linarith [s2lo, s2hi]
  have hy : ∀ v ∈ Icc b a, v - (-a) ≠ 0 := fun v hv => ne_of_gt (by
    have h1 := hv.1; have h2 := hv.2; simp only [Y01.a, Y01.b] at *; linarith [s2lo, s2hi])
  exact intIII isP2a hba hxy hx hy (fun v hv => contAt_P2a (mb_ne_one hxy (hy v hv)))

lemma val10q : FIII P2a 0 (-a) a = (2 : ℝ) * D00.Li2 (1 / 2) * lam + (2 : ℝ) * Li3 (1 / 2) + (2 / 3 : ℝ) * lam ^ 3 := by
  have e0 : mb 0 (-a) a = 1 / 2 := by unfold mb; have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf)
  have l0 : Real.log (a - 0) = (0 : ℝ) * lam + (1 : ℝ) * L :=
    logP 0 1 0 0 _ _ (by norm_num) (by norm_num)
      (by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf))
  have l1 : Real.log (a - (-a)) = (1 : ℝ) * lam + (1 : ℝ) * L :=
    logP 2 1 0 0 _ _ (by norm_num) (by norm_num)
      (by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf))
  have l2 : Real.log (1 / 2) = (-1 : ℝ) * lam + (0 : ℝ) * L :=
    logP 0 0 2 0 _ _ (by norm_num) (by norm_num)
      (by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf))
  have l3 : Real.log (1 - (1 / 2)) = (-1 : ℝ) * lam + (0 : ℝ) * L :=
    logP 0 0 2 0 _ _ (by norm_num) (by norm_num)
      (by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf))
  simp only [FIII, e0, l0, l1, P2a, l2, l3]
  try simp only [Real.log_zero, sub_zero, Real.log_one, D00.Li2_zero', Li3_zero, D00.Li2_one']
  ring

lemma val10p : FIII P2a 0 (-a) b = (-4 / 3 : ℝ) * L ^ 3 + (-3 : ℝ) * L ^ 2 * lam + (2 : ℝ) * L * D00.Li2 ((2 - Real.sqrt 2) / 4) + (3 : ℝ) * D00.Li2 ((2 - Real.sqrt 2) / 4) * lam + (2 : ℝ) * Li3 ((2 - Real.sqrt 2) / 4) + (9 / 4 : ℝ) * lam ^ 3 := by
  have e0 : mb 0 (-a) b = (2 - Real.sqrt 2) / 4 := by unfold mb; have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf)
  have l0 : Real.log (b - 0) = (0 : ℝ) * lam + (-1 : ℝ) * L :=
    logP 0 0 0 1 _ _ (by norm_num) (by norm_num)
      (by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf))
  have l1 : Real.log (b - (-a)) = (3 / 2 : ℝ) * lam + (0 : ℝ) * L :=
    logP 3 0 0 0 _ _ (by norm_num) (by norm_num)
      (by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf))
  have l2 : Real.log ((2 - Real.sqrt 2) / 4) = (-3 / 2 : ℝ) * lam + (-1 : ℝ) * L :=
    logP 0 0 3 1 _ _ (by norm_num) (by norm_num)
      (by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf))
  have l3 : Real.log (1 - ((2 - Real.sqrt 2) / 4)) = (-3 / 2 : ℝ) * lam + (1 : ℝ) * L :=
    logP 0 1 3 0 _ _ (by norm_num) (by norm_num)
      (by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf))
  simp only [FIII, e0, l0, l1, P2a, l2, l3]
  try simp only [Real.log_zero, sub_zero, Real.log_one, D00.Li2_zero', Li3_zero, D00.Li2_one']
  ring

lemma blk11 : IntervalIntegrable (fpr b (-a) 0) volume 1 a ∧
    ∫ v in 1..a, fpr b (-a) 0 v = FG P1n P2a P2a b (-a) 0 a - FG P1n P2a P2a b (-a) 0 1 := by
  have hxy : (b : ℝ) - (-a) ≠ 0 := ne_of_gt (by (try simp only [Y01.a, Y01.b]); linarith [s2lo, s2hi])
  have hzx : (0 : ℝ) - b ≠ 0 := ne_of_lt (by (try simp only [Y01.a, Y01.b]); linarith [s2lo, s2hi])
  have hzy : (0 : ℝ) - (-a) ≠ 0 := ne_of_gt (by (try simp only [Y01.a, Y01.b]); linarith [s2lo, s2hi])
  have hz : ∀ v ∈ Icc 1 a, 0 < v - 0 := fun v hv => by
    have h1 := hv.1; have h2 := hv.2; simp only [Y01.a, Y01.b] at *; linarith [s2lo, s2hi]
  have hy : ∀ v ∈ Icc 1 a, v - (-a) ≠ 0 := fun v hv => ne_of_gt (by
    have h1 := hv.1; have h2 := hv.2; simp only [Y01.a, Y01.b] at *; linarith [s2lo, s2hi])
  have hx : ∀ v ∈ Ioo 1 a, v - b ≠ 0 := fun v hv => ne_of_gt (by
    have h1 := hv.1; have h2 := hv.2; simp only [Y01.a, Y01.b] at *; linarith [s2lo, s2hi])
  exact intG isP1n isP2a isP2a h1a hxy hzx hzy hz hy hx
    (fun v hv => cPhiN (s1_ne_one hxy hzx (hy v hv) (hz v hv).ne'))
    (fun v hv => contAt_P2a (mb_ne_one hxy (hy v hv)))

lemma val11q : FG P1n P2a P2a b (-a) 0 a = (2 : ℝ) * L ^ 3 + (2 : ℝ) * L * D00.Li2 (b) + (-2 : ℝ) * L * D00.Li2 (-a) + (2 : ℝ) * Li3 (b) + (-2 : ℝ) * Li3 (-a) := by
  have e0 : mb b (-a) a = b := by unfold mb; have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf)
  have e1 : s1 b (-a) 0 a = -a := by unfold s1; have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf)
  have l0 : Real.log ((0 - b) / (0 - (-a))) = (0 : ℝ) * lam + (-2 : ℝ) * L :=
    logN 0 0 0 2 _ _ (by norm_num) (by norm_num)
      (by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf))
  have l1 : Real.log (-a) = (0 : ℝ) * lam + (1 : ℝ) * L :=
    logN 0 1 0 0 _ _ (by norm_num) (by norm_num)
      (by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf))
  have l2 : Real.log (1 - (-a)) = (1 / 2 : ℝ) * lam + (1 : ℝ) * L :=
    logP 1 1 0 0 _ _ (by norm_num) (by norm_num)
      (by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf))
  have l3 : Real.log (b) = (0 : ℝ) * lam + (-1 : ℝ) * L :=
    logP 0 0 0 1 _ _ (by norm_num) (by norm_num)
      (by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf))
  have l4 : Real.log (1 - (b)) = (1 / 2 : ℝ) * lam + (-1 : ℝ) * L :=
    logP 1 0 0 1 _ _ (by norm_num) (by norm_num)
      (by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf))
  simp only [FG, e0, e1, l0, Phi, P1n, P2a, l1, l2, l3, l4]
  try simp only [Real.log_zero, sub_zero, Real.log_one, D00.Li2_zero', Li3_zero, D00.Li2_one']
  ring

lemma val11p : FG P1n P2a P2a b (-a) 0 1 = (4 : ℝ) * L ^ 3 + (-4 : ℝ) * L * D00.Li2 (-1) + (4 : ℝ) * L * D00.Li2 (b ^ 2) + (-2 : ℝ) * Li3 (-1) + (2 : ℝ) * Li3 (b ^ 2) := by
  have e0 : mb b (-a) 1 = b ^ 2 := by unfold mb; have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf)
  have e1 : s1 b (-a) 0 1 = -1 := by unfold s1; have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf)
  have l0 : Real.log ((0 - b) / (0 - (-a))) = (0 : ℝ) * lam + (-2 : ℝ) * L :=
    logN 0 0 0 2 _ _ (by norm_num) (by norm_num)
      (by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf))
  have l1 : Real.log (-1) = (0 : ℝ) * lam + (0 : ℝ) * L :=
    logN 0 0 0 0 _ _ (by norm_num) (by norm_num)
      (by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf))
  have l2 : Real.log (1 - (-1)) = (1 : ℝ) * lam + (0 : ℝ) * L :=
    logP 2 0 0 0 _ _ (by norm_num) (by norm_num)
      (by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf))
  have l3 : Real.log (b ^ 2) = (0 : ℝ) * lam + (-2 : ℝ) * L :=
    logP 0 0 0 2 _ _ (by norm_num) (by norm_num)
      (by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf))
  have l4 : Real.log (1 - (b ^ 2)) = (1 : ℝ) * lam + (-1 : ℝ) * L :=
    logP 2 0 0 1 _ _ (by norm_num) (by norm_num)
      (by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf))
  simp only [FG, e0, e1, l0, Phi, P1n, P2a, l1, l2, l3, l4]
  try simp only [Real.log_zero, sub_zero, Real.log_one, D00.Li2_zero', Li3_zero, D00.Li2_one']
  ring

lemma blk12 : IntervalIntegrable (fpr b 0 0) volume 1 a ∧
    ∫ v in 1..a, fpr b 0 0 v = FIV P2a b 0 a - FIV P2a b 0 1 := by
  have hxy : (b : ℝ) - 0 ≠ 0 := ne_of_gt (by (try simp only [Y01.a, Y01.b]); linarith [s2lo, s2hi])
  have hy : ∀ v ∈ Icc 1 a, 0 < v - 0 := fun v hv => by
    have h1 := hv.1; have h2 := hv.2; simp only [Y01.a, Y01.b] at *; linarith [s2lo, s2hi]
  have hx : ∀ v ∈ Ioo 1 a, v - b ≠ 0 := fun v hv => ne_of_gt (by
    have h1 := hv.1; have h2 := hv.2; simp only [Y01.a, Y01.b] at *; linarith [s2lo, s2hi])
  exact intIV isP2a h1a hxy hy hx (fun v hv => contAt_P2a (mb_ne_one hxy (hy v hv).ne'))

lemma val12q : FIV P2a b 0 a = (2 : ℝ) * L ^ 3 + (-4 : ℝ) * L ^ 2 * lam + (2 : ℝ) * L * D00.Li2 (2 * b) + (2 : ℝ) * L * lam ^ 2 + (-2 : ℝ) * D00.Li2 (2 * b) * lam + (2 : ℝ) * Li3 (2 * b) := by
  have e0 : mb b 0 a = 2 * b := by unfold mb; have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf)
  have l0 : Real.log (2 * b) = (1 : ℝ) * lam + (-1 : ℝ) * L :=
    logP 2 0 0 1 _ _ (by norm_num) (by norm_num)
      (by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf))
  have l1 : Real.log (1 - (2 * b)) = (0 : ℝ) * lam + (-2 : ℝ) * L :=
    logP 0 0 0 2 _ _ (by norm_num) (by norm_num)
      (by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf))
  simp only [FIV, e0, P2a, l0, l1]
  try simp only [Real.log_zero, sub_zero, Real.log_one, D00.Li2_zero', Li3_zero, D00.Li2_one']
  ring

lemma val12p : FIV P2a b 0 1 = (1 : ℝ) * L ^ 3 + (-1 : ℝ) * L ^ 2 * lam + (2 : ℝ) * L * D00.Li2 (2 - Real.sqrt 2) + (1 / 4 : ℝ) * L * lam ^ 2 + (-1 : ℝ) * D00.Li2 (2 - Real.sqrt 2) * lam + (2 : ℝ) * Li3 (2 - Real.sqrt 2) := by
  have e0 : mb b 0 1 = 2 - Real.sqrt 2 := by unfold mb; have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf)
  have l0 : Real.log (2 - Real.sqrt 2) = (1 / 2 : ℝ) * lam + (-1 : ℝ) * L :=
    logP 1 0 0 1 _ _ (by norm_num) (by norm_num)
      (by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf))
  have l1 : Real.log (1 - (2 - Real.sqrt 2)) = (0 : ℝ) * lam + (-1 : ℝ) * L :=
    logP 0 0 0 1 _ _ (by norm_num) (by norm_num)
      (by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf))
  simp only [FIV, e0, P2a, l0, l1]
  try simp only [Real.log_zero, sub_zero, Real.log_one, D00.Li2_zero', Li3_zero, D00.Li2_one']
  ring

lemma blk13 : IntervalIntegrable (fpr (-a) b 0) volume 1 a ∧
    ∫ v in 1..a, fpr (-a) b 0 v = FG P1n P2a P2g (-a) b 0 a - FG P1n P2a P2g (-a) b 0 1 := by
  have hxy : ((-a) : ℝ) - b ≠ 0 := ne_of_lt (by (try simp only [Y01.a, Y01.b]); linarith [s2lo, s2hi])
  have hzx : (0 : ℝ) - (-a) ≠ 0 := ne_of_gt (by (try simp only [Y01.a, Y01.b]); linarith [s2lo, s2hi])
  have hzy : (0 : ℝ) - b ≠ 0 := ne_of_lt (by (try simp only [Y01.a, Y01.b]); linarith [s2lo, s2hi])
  have hz : ∀ v ∈ Icc 1 a, 0 < v - 0 := fun v hv => by
    have h1 := hv.1; have h2 := hv.2; simp only [Y01.a, Y01.b] at *; linarith [s2lo, s2hi]
  have hy : ∀ v ∈ Icc 1 a, v - b ≠ 0 := fun v hv => ne_of_gt (by
    have h1 := hv.1; have h2 := hv.2; simp only [Y01.a, Y01.b] at *; linarith [s2lo, s2hi])
  have hx : ∀ v ∈ Ioo 1 a, v - (-a) ≠ 0 := fun v hv => ne_of_gt (by
    have h1 := hv.1; have h2 := hv.2; simp only [Y01.a, Y01.b] at *; linarith [s2lo, s2hi])
  have hx2 : ∀ v ∈ Icc 1 a, v - (-a) ≠ 0 := fun v hv => ne_of_gt (by
    have h1 := hv.1; have h2 := hv.2; simp only [Y01.a, Y01.b] at *; linarith [s2lo, s2hi])
  exact intG isP1n isP2a isP2g h1a hxy hzx hzy hz hy hx
    (fun v hv => cPhiN (s1_ne_one hxy hzx (hy v hv) (hz v hv).ne'))
    (fun v hv => contAt_P2g (mb_ne_zero (hx2 v hv) (hy v hv)) (mb_ne_one hxy (hy v hv)))

lemma val13q : FG P1n P2a P2g (-a) b 0 a = (2 / 3 : ℝ) * L ^ 3 + (2 : ℝ) * L * D00.Li2 (b) + (2 : ℝ) * L * D00.Li2 (-b) + (2 : ℝ) * Li3 (b) + (-2 : ℝ) * Li3 (-b) := by
  have e0 : mb (-a) b a = (Real.sqrt 2 ^ 0 * a ^ 1 / (Real.sqrt 2 ^ 0 * a ^ 0)) := by unfold mb; have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf)
  have e1 : s1 (-a) b 0 a = -b := by unfold s1; have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf)
  have l0 : Real.log ((0 - (-a)) / (0 - b)) = (0 : ℝ) * lam + (2 : ℝ) * L :=
    logN 0 2 0 0 _ _ (by norm_num) (by norm_num)
      (by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf))
  have l1 : Real.log (-b) = (0 : ℝ) * lam + (-1 : ℝ) * L :=
    logN 0 0 0 1 _ _ (by norm_num) (by norm_num)
      (by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf))
  have l2 : Real.log (1 - (-b)) = (1 / 2 : ℝ) * lam + (0 : ℝ) * L :=
    logP 1 0 0 0 _ _ (by norm_num) (by norm_num)
      (by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf))
  have e2 : ((Real.sqrt 2 ^ 0 * a ^ 1 / (Real.sqrt 2 ^ 0 * a ^ 0)))⁻¹ = b := by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf)
  have l3 : Real.log (b) = (0 : ℝ) * lam + (-1 : ℝ) * L :=
    logP 0 0 0 1 _ _ (by norm_num) (by norm_num)
      (by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf))
  have l4 : Real.log (1 - (b)) = (1 / 2 : ℝ) * lam + (-1 : ℝ) * L :=
    logP 1 0 0 1 _ _ (by norm_num) (by norm_num)
      (by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf))
  simp only [FG, e0, e1, l0, Phi, P1n, P2a, l1, l2, P2g, E3, e2, l3, l4]
  try simp only [Real.log_zero, sub_zero, Real.log_one, D00.Li2_zero', Li3_zero, D00.Li2_one']
  ring

lemma val13p : FG P1n P2a P2g (-a) b 0 1 = (4 / 3 : ℝ) * L ^ 3 + (4 : ℝ) * L * D00.Li2 (-1) + (4 : ℝ) * L * D00.Li2 (b ^ 2) + (-2 : ℝ) * Li3 (-1) + (2 : ℝ) * Li3 (b ^ 2) := by
  have e0 : mb (-a) b 1 = (Real.sqrt 2 ^ 0 * a ^ 2 / (Real.sqrt 2 ^ 0 * a ^ 0)) := by unfold mb; have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf)
  have e1 : s1 (-a) b 0 1 = -1 := by unfold s1; have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf)
  have l0 : Real.log ((0 - (-a)) / (0 - b)) = (0 : ℝ) * lam + (2 : ℝ) * L :=
    logN 0 2 0 0 _ _ (by norm_num) (by norm_num)
      (by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf))
  have l1 : Real.log (-1) = (0 : ℝ) * lam + (0 : ℝ) * L :=
    logN 0 0 0 0 _ _ (by norm_num) (by norm_num)
      (by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf))
  have l2 : Real.log (1 - (-1)) = (1 : ℝ) * lam + (0 : ℝ) * L :=
    logP 2 0 0 0 _ _ (by norm_num) (by norm_num)
      (by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf))
  have e2 : ((Real.sqrt 2 ^ 0 * a ^ 2 / (Real.sqrt 2 ^ 0 * a ^ 0)))⁻¹ = b ^ 2 := by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf)
  have l3 : Real.log (b ^ 2) = (0 : ℝ) * lam + (-2 : ℝ) * L :=
    logP 0 0 0 2 _ _ (by norm_num) (by norm_num)
      (by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf))
  have l4 : Real.log (1 - (b ^ 2)) = (1 : ℝ) * lam + (-1 : ℝ) * L :=
    logP 2 0 0 1 _ _ (by norm_num) (by norm_num)
      (by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf))
  simp only [FG, e0, e1, l0, Phi, P1n, P2a, l1, l2, P2g, E3, e2, l3, l4]
  try simp only [Real.log_zero, sub_zero, Real.log_one, D00.Li2_zero', Li3_zero, D00.Li2_one']
  ring

lemma blk14 : IntervalIntegrable (fpr 0 b 0) volume 1 a ∧
    ∫ v in 1..a, fpr 0 b 0 v = FIII P2g 0 b a - FIII P2g 0 b 1 := by
  have hxy : (0 : ℝ) - b ≠ 0 := ne_of_lt (by (try simp only [Y01.a, Y01.b]); linarith [s2lo, s2hi])
  have hx : ∀ v ∈ Icc 1 a, 0 < v - 0 := fun v hv => by
    have h1 := hv.1; have h2 := hv.2; simp only [Y01.a, Y01.b] at *; linarith [s2lo, s2hi]
  have hy : ∀ v ∈ Icc 1 a, v - b ≠ 0 := fun v hv => ne_of_gt (by
    have h1 := hv.1; have h2 := hv.2; simp only [Y01.a, Y01.b] at *; linarith [s2lo, s2hi])
  exact intIII isP2g h1a hxy hx hy (fun v hv => contAt_P2g (mb_ne_zero (hx v hv).ne' (hy v hv)) (mb_ne_one hxy (hy v hv)))

lemma val14q : FIII P2g 0 b a = (2 : ℝ) * L ^ 3 + (-4 : ℝ) * L ^ 2 * lam + (2 : ℝ) * L * D00.Li2 (2 * b) + (2 : ℝ) * L * lam ^ 2 + (-2 : ℝ) * D00.Li2 (2 * b) * lam + (2 : ℝ) * Li3 (2 * b) := by
  have e0 : mb 0 b a = (Real.sqrt 2 ^ 0 * a ^ 1 / (Real.sqrt 2 ^ 2 * a ^ 0)) := by unfold mb; have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf)
  have l0 : Real.log (a - 0) = (0 : ℝ) * lam + (1 : ℝ) * L :=
    logP 0 1 0 0 _ _ (by norm_num) (by norm_num)
      (by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf))
  have l1 : Real.log (a - b) = (1 : ℝ) * lam + (0 : ℝ) * L :=
    logP 2 0 0 0 _ _ (by norm_num) (by norm_num)
      (by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf))
  have e1 : ((Real.sqrt 2 ^ 0 * a ^ 1 / (Real.sqrt 2 ^ 2 * a ^ 0)))⁻¹ = 2 * b := by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf)
  have l2 : Real.log (2 * b) = (1 : ℝ) * lam + (-1 : ℝ) * L :=
    logP 2 0 0 1 _ _ (by norm_num) (by norm_num)
      (by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf))
  have l3 : Real.log (1 - (2 * b)) = (0 : ℝ) * lam + (-2 : ℝ) * L :=
    logP 0 0 0 2 _ _ (by norm_num) (by norm_num)
      (by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf))
  simp only [FIII, e0, l0, l1, P2g, E3, e1, l2, l3]
  try simp only [Real.log_zero, sub_zero, Real.log_one, D00.Li2_zero', Li3_zero, D00.Li2_one']
  ring

lemma val14p : FIII P2g 0 b 1 = (1 : ℝ) * L ^ 3 + (-1 : ℝ) * L ^ 2 * lam + (2 : ℝ) * L * D00.Li2 (2 - Real.sqrt 2) + (1 / 4 : ℝ) * L * lam ^ 2 + (-1 : ℝ) * D00.Li2 (2 - Real.sqrt 2) * lam + (2 : ℝ) * Li3 (2 - Real.sqrt 2) := by
  have e0 : mb 0 b 1 = (Real.sqrt 2 ^ 0 * a ^ 1 / (Real.sqrt 2 ^ 1 * a ^ 0)) := by unfold mb; have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf)
  have l0 : Real.log (1 - 0) = (0 : ℝ) * lam + (0 : ℝ) * L :=
    logP 0 0 0 0 _ _ (by norm_num) (by norm_num)
      (by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf))
  have l1 : Real.log (1 - b) = (1 / 2 : ℝ) * lam + (-1 : ℝ) * L :=
    logP 1 0 0 1 _ _ (by norm_num) (by norm_num)
      (by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf))
  have e1 : ((Real.sqrt 2 ^ 0 * a ^ 1 / (Real.sqrt 2 ^ 1 * a ^ 0)))⁻¹ = 2 - Real.sqrt 2 := by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf)
  have l2 : Real.log (2 - Real.sqrt 2) = (1 / 2 : ℝ) * lam + (-1 : ℝ) * L :=
    logP 1 0 0 1 _ _ (by norm_num) (by norm_num)
      (by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf))
  have l3 : Real.log (1 - (2 - Real.sqrt 2)) = (0 : ℝ) * lam + (-1 : ℝ) * L :=
    logP 0 0 0 1 _ _ (by norm_num) (by norm_num)
      (by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf))
  simp only [FIII, e0, l0, l1, P2g, E3, e1, l2, l3]
  try simp only [Real.log_zero, sub_zero, Real.log_one, D00.Li2_zero', Li3_zero, D00.Li2_one']
  ring

lemma blk15 : IntervalIntegrable (fcst (-a)) volume b a ∧
    ∫ v in b..a, fcst (-a) v = FK (-a) a - FK (-a) b := by
  have hz : ∀ v ∈ Icc b a, 0 < v - (-a) := fun v hv => by
    have h1 := hv.1; have h2 := hv.2; simp only [Y01.a, Y01.b] at *; linarith [s2lo, s2hi]
  exact intK hba hz

lemma val15q : FK (-a) a = (1 : ℝ) * L + (1 : ℝ) * lam := by
  have l0 : Real.log (a - (-a)) = (1 : ℝ) * lam + (1 : ℝ) * L :=
    logP 2 1 0 0 _ _ (by norm_num) (by norm_num)
      (by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf))
  simp only [FK, l0]
  try simp only [Real.log_zero, sub_zero, Real.log_one, D00.Li2_zero', Li3_zero, D00.Li2_one']
  ring

lemma val15p : FK (-a) b = (3 / 2 : ℝ) * lam := by
  have l0 : Real.log (b - (-a)) = (3 / 2 : ℝ) * lam + (0 : ℝ) * L :=
    logP 3 0 0 0 _ _ (by norm_num) (by norm_num)
      (by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf))
  simp only [FK, l0]
  try simp only [Real.log_zero, sub_zero, Real.log_one, D00.Li2_zero', Li3_zero, D00.Li2_one']
  ring

lemma blk16 : IntervalIntegrable (flin b (-a)) volume b a ∧
    ∫ v in b..a, flin b (-a) v = FL P1n b (-a) a - FL P1n b (-a) b := by
  have hxz : ((-a) : ℝ) - b ≠ 0 := ne_of_lt (by (try simp only [Y01.a, Y01.b]); linarith [s2lo, s2hi])
  have hz : ∀ v ∈ Icc b a, 0 < v - (-a) := fun v hv => by
    have h1 := hv.1; have h2 := hv.2; simp only [Y01.a, Y01.b] at *; linarith [s2lo, s2hi]
  have hx : ∀ v ∈ Ioo b a, v - b ≠ 0 := fun v hv => ne_of_gt (by
    have h1 := hv.1; have h2 := hv.2; simp only [Y01.a, Y01.b] at *; linarith [s2lo, s2hi])
  exact intL isP1n hba hxz hz hx (fun v hv => cPsiN (sA_ne_one hxz (hz v hv).ne'))

lemma val16q : FL P1n b (-a) a = (1 : ℝ) * L * lam + (1 : ℝ) * D00.Li2 (-(1 / Real.sqrt 2)) + (-1 / 2 : ℝ) * lam ^ 2 := by
  have e0 : sA b (-a) a = -(1 / Real.sqrt 2) := by unfold sA; have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf)
  have l0 : Real.log ((-a) - b) = (3 / 2 : ℝ) * lam + (0 : ℝ) * L :=
    logN 3 0 0 0 _ _ (by norm_num) (by norm_num)
      (by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf))
  have l1 : Real.log (-(1 / Real.sqrt 2)) = (-1 / 2 : ℝ) * lam + (0 : ℝ) * L :=
    logN 0 0 1 0 _ _ (by norm_num) (by norm_num)
      (by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf))
  have l2 : Real.log (1 - (-(1 / Real.sqrt 2))) = (-1 / 2 : ℝ) * lam + (1 : ℝ) * L :=
    logP 0 1 1 0 _ _ (by norm_num) (by norm_num)
      (by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf))
  simp only [FL, e0, l0, Psi, P1n, l1, l2]
  try simp only [Real.log_zero, sub_zero, Real.log_one, D00.Li2_zero', Li3_zero, D00.Li2_one']
  ring

lemma val16p : FL P1n b (-a) b = 0 := by
  have e0 : sA b (-a) b = 0 := by unfold sA; have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf)
  have l0 : Real.log ((-a) - b) = (3 / 2 : ℝ) * lam + (0 : ℝ) * L :=
    logN 3 0 0 0 _ _ (by norm_num) (by norm_num)
      (by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf))
  simp only [FL, e0, l0, Psi, P1n]
  try simp only [Real.log_zero, sub_zero, Real.log_one, D00.Li2_zero', Li3_zero, D00.Li2_one']
  ring

lemma blk17 : IntervalIntegrable (flin (-a) (-a)) volume b a ∧
    ∫ v in b..a, flin (-a) (-a) v = FH (-a) a - FH (-a) b := by
  have hx : ∀ v ∈ Icc b a, 0 < v - (-a) := fun v hv => by
    have h1 := hv.1; have h2 := hv.2; simp only [Y01.a, Y01.b] at *; linarith [s2lo, s2hi]
  exact intH hba hx

lemma val17q : FH (-a) a = (1 / 2 : ℝ) * L ^ 2 + (1 : ℝ) * L * lam + (1 / 2 : ℝ) * lam ^ 2 := by
  have l0 : Real.log (a - (-a)) = (1 : ℝ) * lam + (1 : ℝ) * L :=
    logP 2 1 0 0 _ _ (by norm_num) (by norm_num)
      (by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf))
  simp only [FH, l0]
  try simp only [Real.log_zero, sub_zero, Real.log_one, D00.Li2_zero', Li3_zero, D00.Li2_one']
  ring

lemma val17p : FH (-a) b = (9 / 8 : ℝ) * lam ^ 2 := by
  have l0 : Real.log (b - (-a)) = (3 / 2 : ℝ) * lam + (0 : ℝ) * L :=
    logP 3 0 0 0 _ _ (by norm_num) (by norm_num)
      (by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf))
  simp only [FH, l0]
  try simp only [Real.log_zero, sub_zero, Real.log_one, D00.Li2_zero', Li3_zero, D00.Li2_one']
  ring

lemma blk18 : IntervalIntegrable (flin 0 (-a)) volume b a ∧
    ∫ v in b..a, flin 0 (-a) v = FL P1n 0 (-a) a - FL P1n 0 (-a) b := by
  have hxz : ((-a) : ℝ) - 0 ≠ 0 := ne_of_lt (by (try simp only [Y01.a, Y01.b]); linarith [s2lo, s2hi])
  have hz : ∀ v ∈ Icc b a, 0 < v - (-a) := fun v hv => by
    have h1 := hv.1; have h2 := hv.2; simp only [Y01.a, Y01.b] at *; linarith [s2lo, s2hi]
  have hx : ∀ v ∈ Ioo b a, v - 0 ≠ 0 := fun v hv => ne_of_gt (by
    have h1 := hv.1; have h2 := hv.2; simp only [Y01.a, Y01.b] at *; linarith [s2lo, s2hi])
  exact intL isP1n hba hxz hz hx (fun v hv => cPsiN (sA_ne_one hxz (hz v hv).ne'))

lemma val18q : FL P1n 0 (-a) a = (1 : ℝ) * L * lam + (1 : ℝ) * D00.Li2 (-1) := by
  have e0 : sA 0 (-a) a = -1 := by unfold sA; have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf)
  have l0 : Real.log ((-a) - 0) = (0 : ℝ) * lam + (1 : ℝ) * L :=
    logN 0 1 0 0 _ _ (by norm_num) (by norm_num)
      (by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf))
  have l1 : Real.log (-1) = (0 : ℝ) * lam + (0 : ℝ) * L :=
    logN 0 0 0 0 _ _ (by norm_num) (by norm_num)
      (by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf))
  have l2 : Real.log (1 - (-1)) = (1 : ℝ) * lam + (0 : ℝ) * L :=
    logP 2 0 0 0 _ _ (by norm_num) (by norm_num)
      (by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf))
  simp only [FL, e0, l0, Psi, P1n, l1, l2]
  try simp only [Real.log_zero, sub_zero, Real.log_one, D00.Li2_zero', Li3_zero, D00.Li2_one']
  ring

lemma val18p : FL P1n 0 (-a) b = (1 : ℝ) * L ^ 2 + (-3 / 2 : ℝ) * L * lam + (1 : ℝ) * D00.Li2 (-(b ^ 2)) := by
  have e0 : sA 0 (-a) b = -(b ^ 2) := by unfold sA; have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf)
  have l0 : Real.log ((-a) - 0) = (0 : ℝ) * lam + (1 : ℝ) * L :=
    logN 0 1 0 0 _ _ (by norm_num) (by norm_num)
      (by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf))
  have l1 : Real.log (-(b ^ 2)) = (0 : ℝ) * lam + (-2 : ℝ) * L :=
    logN 0 0 0 2 _ _ (by norm_num) (by norm_num)
      (by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf))
  have l2 : Real.log (1 - (-(b ^ 2))) = (3 / 2 : ℝ) * lam + (-1 : ℝ) * L :=
    logP 3 0 0 1 _ _ (by norm_num) (by norm_num)
      (by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf))
  simp only [FL, e0, l0, Psi, P1n, l1, l2]
  try simp only [Real.log_zero, sub_zero, Real.log_one, D00.Li2_zero', Li3_zero, D00.Li2_one']
  ring

lemma blk19 : IntervalIntegrable (fsq b (-a)) volume b a ∧
    ∫ v in b..a, fsq b (-a) v = FA P1n P2a b (-a) a - FA P1n P2a b (-a) b := by
  have hxz : ((-a) : ℝ) - b ≠ 0 := ne_of_lt (by (try simp only [Y01.a, Y01.b]); linarith [s2lo, s2hi])
  have hz : ∀ v ∈ Icc b a, 0 < v - (-a) := fun v hv => by
    have h1 := hv.1; have h2 := hv.2; simp only [Y01.a, Y01.b] at *; linarith [s2lo, s2hi]
  have hx : ∀ v ∈ Ioo b a, v - b ≠ 0 := fun v hv => ne_of_gt (by
    have h1 := hv.1; have h2 := hv.2; simp only [Y01.a, Y01.b] at *; linarith [s2lo, s2hi])
  exact intA isP1n isP2a hba hxz hz hx (fun v hv => cPhiN (sA_ne_one hxz (hz v hv).ne'))

lemma val19q : FA P1n P2a b (-a) a = (1 : ℝ) * L * lam ^ 2 + (2 : ℝ) * D00.Li2 (-(1 / Real.sqrt 2)) * lam + (-2 : ℝ) * Li3 (-(1 / Real.sqrt 2)) + (-1 / 2 : ℝ) * lam ^ 3 := by
  have e0 : sA b (-a) a = -(1 / Real.sqrt 2) := by unfold sA; have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf)
  have l0 : Real.log ((-a) - b) = (3 / 2 : ℝ) * lam + (0 : ℝ) * L :=
    logN 3 0 0 0 _ _ (by norm_num) (by norm_num)
      (by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf))
  have l1 : Real.log (-(1 / Real.sqrt 2)) = (-1 / 2 : ℝ) * lam + (0 : ℝ) * L :=
    logN 0 0 1 0 _ _ (by norm_num) (by norm_num)
      (by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf))
  have l2 : Real.log (1 - (-(1 / Real.sqrt 2))) = (-1 / 2 : ℝ) * lam + (1 : ℝ) * L :=
    logP 0 1 1 0 _ _ (by norm_num) (by norm_num)
      (by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf))
  simp only [FA, e0, l0, Phi, P1n, P2a, l1, l2]
  try simp only [Real.log_zero, sub_zero, Real.log_one, D00.Li2_zero', Li3_zero, D00.Li2_one']
  ring

lemma val19p : FA P1n P2a b (-a) b = 0 := by
  have e0 : sA b (-a) b = 0 := by unfold sA; have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf)
  have l0 : Real.log ((-a) - b) = (3 / 2 : ℝ) * lam + (0 : ℝ) * L :=
    logN 3 0 0 0 _ _ (by norm_num) (by norm_num)
      (by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf))
  simp only [FA, e0, l0, Phi, P1n, P2a]
  try simp only [Real.log_zero, sub_zero, Real.log_one, D00.Li2_zero', Li3_zero, D00.Li2_one']
  ring

lemma blk20 : IntervalIntegrable (fsq (-a) (-a)) volume b a ∧
    ∫ v in b..a, fsq (-a) (-a) v = FC (-a) a - FC (-a) b := by
  have hx : ∀ v ∈ Icc b a, 0 < v - (-a) := fun v hv => by
    have h1 := hv.1; have h2 := hv.2; simp only [Y01.a, Y01.b] at *; linarith [s2lo, s2hi]
  exact intC hba hx

lemma val20q : FC (-a) a = (1 / 3 : ℝ) * L ^ 3 + (1 : ℝ) * L ^ 2 * lam + (1 : ℝ) * L * lam ^ 2 + (1 / 3 : ℝ) * lam ^ 3 := by
  have l0 : Real.log (a - (-a)) = (1 : ℝ) * lam + (1 : ℝ) * L :=
    logP 2 1 0 0 _ _ (by norm_num) (by norm_num)
      (by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf))
  simp only [FC, l0]
  try simp only [Real.log_zero, sub_zero, Real.log_one, D00.Li2_zero', Li3_zero, D00.Li2_one']
  ring

lemma val20p : FC (-a) b = (9 / 8 : ℝ) * lam ^ 3 := by
  have l0 : Real.log (b - (-a)) = (3 / 2 : ℝ) * lam + (0 : ℝ) * L :=
    logP 3 0 0 0 _ _ (by norm_num) (by norm_num)
      (by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf))
  simp only [FC, l0]
  try simp only [Real.log_zero, sub_zero, Real.log_one, D00.Li2_zero', Li3_zero, D00.Li2_one']
  ring

lemma blk21 : IntervalIntegrable (fpr b (-a) (-a)) volume b 1 ∧
    ∫ v in b..1, fpr b (-a) (-a) v = FIV P2a b (-a) 1 - FIV P2a b (-a) b := by
  have hxy : (b : ℝ) - (-a) ≠ 0 := ne_of_gt (by (try simp only [Y01.a, Y01.b]); linarith [s2lo, s2hi])
  have hy : ∀ v ∈ Icc b 1, 0 < v - (-a) := fun v hv => by
    have h1 := hv.1; have h2 := hv.2; simp only [Y01.a, Y01.b] at *; linarith [s2lo, s2hi]
  have hx : ∀ v ∈ Ioo b 1, v - b ≠ 0 := fun v hv => ne_of_gt (by
    have h1 := hv.1; have h2 := hv.2; simp only [Y01.a, Y01.b] at *; linarith [s2lo, s2hi])
  exact intIV isP2a hb1 hxy hy hx (fun v hv => contAt_P2a (mb_ne_one hxy (hy v hv).ne'))

lemma val21q : FIV P2a b (-a) 1 = (4 : ℝ) * L ^ 3 + (-4 : ℝ) * L ^ 2 * lam + (4 : ℝ) * L * D00.Li2 (b ^ 2) + (2 : ℝ) * Li3 (b ^ 2) := by
  have e0 : mb b (-a) 1 = b ^ 2 := by unfold mb; have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf)
  have l0 : Real.log (b ^ 2) = (0 : ℝ) * lam + (-2 : ℝ) * L :=
    logP 0 0 0 2 _ _ (by norm_num) (by norm_num)
      (by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf))
  have l1 : Real.log (1 - (b ^ 2)) = (1 : ℝ) * lam + (-1 : ℝ) * L :=
    logP 2 0 0 1 _ _ (by norm_num) (by norm_num)
      (by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf))
  simp only [FIV, e0, P2a, l0, l1]
  try simp only [Real.log_zero, sub_zero, Real.log_one, D00.Li2_zero', Li3_zero, D00.Li2_one']
  ring

lemma val21p : FIV P2a b (-a) b = 0 := by
  have e0 : mb b (-a) b = 0 := by unfold mb; have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf)
  simp only [FIV, e0, P2a]
  try simp only [Real.log_zero, sub_zero, Real.log_one, D00.Li2_zero', Li3_zero, D00.Li2_one']
  ring

lemma blk22 : IntervalIntegrable (fsq 0 (-a)) volume b a ∧
    ∫ v in b..a, fsq 0 (-a) v = FA P1n P2a 0 (-a) a - FA P1n P2a 0 (-a) b := by
  have hxz : ((-a) : ℝ) - 0 ≠ 0 := ne_of_lt (by (try simp only [Y01.a, Y01.b]); linarith [s2lo, s2hi])
  have hz : ∀ v ∈ Icc b a, 0 < v - (-a) := fun v hv => by
    have h1 := hv.1; have h2 := hv.2; simp only [Y01.a, Y01.b] at *; linarith [s2lo, s2hi]
  have hx : ∀ v ∈ Ioo b a, v - 0 ≠ 0 := fun v hv => ne_of_gt (by
    have h1 := hv.1; have h2 := hv.2; simp only [Y01.a, Y01.b] at *; linarith [s2lo, s2hi])
  exact intA isP1n isP2a hba hxz hz hx (fun v hv => cPhiN (sA_ne_one hxz (hz v hv).ne'))

lemma val22q : FA P1n P2a 0 (-a) a = (1 : ℝ) * L ^ 2 * lam + (2 : ℝ) * L * D00.Li2 (-1) + (-2 : ℝ) * Li3 (-1) := by
  have e0 : sA 0 (-a) a = -1 := by unfold sA; have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf)
  have l0 : Real.log ((-a) - 0) = (0 : ℝ) * lam + (1 : ℝ) * L :=
    logN 0 1 0 0 _ _ (by norm_num) (by norm_num)
      (by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf))
  have l1 : Real.log (-1) = (0 : ℝ) * lam + (0 : ℝ) * L :=
    logN 0 0 0 0 _ _ (by norm_num) (by norm_num)
      (by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf))
  have l2 : Real.log (1 - (-1)) = (1 : ℝ) * lam + (0 : ℝ) * L :=
    logP 2 0 0 0 _ _ (by norm_num) (by norm_num)
      (by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf))
  simp only [FA, e0, l0, Phi, P1n, P2a, l1, l2]
  try simp only [Real.log_zero, sub_zero, Real.log_one, D00.Li2_zero', Li3_zero, D00.Li2_one']
  ring

lemma val22p : FA P1n P2a 0 (-a) b = (-1 : ℝ) * L ^ 3 + (3 / 2 : ℝ) * L ^ 2 * lam + (-2 : ℝ) * L * D00.Li2 (-(b ^ 2)) + (-2 : ℝ) * Li3 (-(b ^ 2)) := by
  have e0 : sA 0 (-a) b = -(b ^ 2) := by unfold sA; have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf)
  have l0 : Real.log ((-a) - 0) = (0 : ℝ) * lam + (1 : ℝ) * L :=
    logN 0 1 0 0 _ _ (by norm_num) (by norm_num)
      (by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf))
  have l1 : Real.log (-(b ^ 2)) = (0 : ℝ) * lam + (-2 : ℝ) * L :=
    logN 0 0 0 2 _ _ (by norm_num) (by norm_num)
      (by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf))
  have l2 : Real.log (1 - (-(b ^ 2))) = (3 / 2 : ℝ) * lam + (-1 : ℝ) * L :=
    logP 3 0 0 1 _ _ (by norm_num) (by norm_num)
      (by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf))
  simp only [FA, e0, l0, Phi, P1n, P2a, l1, l2]
  try simp only [Real.log_zero, sub_zero, Real.log_one, D00.Li2_zero', Li3_zero, D00.Li2_one']
  ring

lemma blk23 : IntervalIntegrable (fpr b 0 (-a)) volume b 1 ∧
    ∫ v in b..1, fpr b 0 (-a) v = FG P1p P2a P2a b 0 (-a) 1 - FG P1p P2a P2a b 0 (-a) b := by
  have hxy : (b : ℝ) - 0 ≠ 0 := ne_of_gt (by (try simp only [Y01.a, Y01.b]); linarith [s2lo, s2hi])
  have hzx : ((-a) : ℝ) - b ≠ 0 := ne_of_lt (by (try simp only [Y01.a, Y01.b]); linarith [s2lo, s2hi])
  have hzy : ((-a) : ℝ) - 0 ≠ 0 := ne_of_lt (by (try simp only [Y01.a, Y01.b]); linarith [s2lo, s2hi])
  have hz : ∀ v ∈ Icc b 1, 0 < v - (-a) := fun v hv => by
    have h1 := hv.1; have h2 := hv.2; simp only [Y01.a, Y01.b] at *; linarith [s2lo, s2hi]
  have hy : ∀ v ∈ Icc b 1, v - 0 ≠ 0 := fun v hv => ne_of_gt (by
    have h1 := hv.1; have h2 := hv.2; simp only [Y01.a, Y01.b] at *; linarith [s2lo, s2hi])
  have hx : ∀ v ∈ Ioo b 1, v - b ≠ 0 := fun v hv => ne_of_gt (by
    have h1 := hv.1; have h2 := hv.2; simp only [Y01.a, Y01.b] at *; linarith [s2lo, s2hi])
  exact intG isP1p isP2a isP2a hb1 hxy hzx hzy hz hy hx
    (fun v hv => cPhiM (s1_ne_one hxy hzx (hy v hv) (hz v hv).ne'))
    (fun v hv => contAt_P2a (mb_ne_one hxy (hy v hv)))

lemma val23q : FG P1p P2a P2a b 0 (-a) 1 = (1 : ℝ) * L ^ 3 + (-2 : ℝ) * L ^ 2 * lam + (2 : ℝ) * L * D00.Li2 (1 / 2) + (2 : ℝ) * L * D00.Li2 (2 - Real.sqrt 2) + (13 / 4 : ℝ) * L * lam ^ 2 + (-5 : ℝ) * D00.Li2 (1 / 2) * lam + (-1 : ℝ) * D00.Li2 (2 - Real.sqrt 2) * lam + (-2 : ℝ) * Li3 (1 / 2) + (2 : ℝ) * Li3 (2 - Real.sqrt 2) + (-13 / 4 : ℝ) * lam ^ 3 := by
  have e0 : mb b 0 1 = 2 - Real.sqrt 2 := by unfold mb; have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf)
  have e1 : s1 b 0 (-a) 1 = 1 / 2 := by unfold s1; have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf)
  have l0 : Real.log (((-a) - b) / ((-a) - 0)) = (3 / 2 : ℝ) * lam + (-1 : ℝ) * L :=
    logP 3 0 0 1 _ _ (by norm_num) (by norm_num)
      (by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf))
  have e2 : (1:ℝ) - (1 / 2) = 1 / 2 := by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf)
  have l1 : Real.log (1 / 2) = (-1 : ℝ) * lam + (0 : ℝ) * L :=
    logP 0 0 2 0 _ _ (by norm_num) (by norm_num)
      (by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf))
  have l2 : Real.log (2 - Real.sqrt 2) = (1 / 2 : ℝ) * lam + (-1 : ℝ) * L :=
    logP 1 0 0 1 _ _ (by norm_num) (by norm_num)
      (by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf))
  have l3 : Real.log (1 - (2 - Real.sqrt 2)) = (0 : ℝ) * lam + (-1 : ℝ) * L :=
    logP 0 0 0 1 _ _ (by norm_num) (by norm_num)
      (by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf))
  simp only [FG, e0, e1, l0, Phi, P1p, P2a, e2, l1, l2, l3]
  try simp only [Real.log_zero, sub_zero, Real.log_one, D00.Li2_zero', Li3_zero, D00.Li2_one']
  ring

lemma val23p : FG P1p P2a P2a b 0 (-a) b = (1 / 3 : ℝ) * L * (Real.pi ^ 2) + (-1 / 2 : ℝ) * lam * (Real.pi ^ 2) := by
  have e0 : mb b 0 b = 0 := by unfold mb; have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf)
  have e1 : s1 b 0 (-a) b = 0 := by unfold s1; have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf)
  have l0 : Real.log (((-a) - b) / ((-a) - 0)) = (3 / 2 : ℝ) * lam + (-1 : ℝ) * L :=
    logP 3 0 0 1 _ _ (by norm_num) (by norm_num)
      (by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf))
  simp only [FG, e0, e1, l0, Phi, P1p, P2a]
  try simp only [Real.log_zero, sub_zero, Real.log_one, D00.Li2_zero', Li3_zero, D00.Li2_one']
  ring

lemma blk24 : IntervalIntegrable (fpr (-a) 0 (-a)) volume b a ∧
    ∫ v in b..a, fpr (-a) 0 (-a) v = FIII P2g (-a) 0 a - FIII P2g (-a) 0 b := by
  have hxy : ((-a) : ℝ) - 0 ≠ 0 := ne_of_lt (by (try simp only [Y01.a, Y01.b]); linarith [s2lo, s2hi])
  have hx : ∀ v ∈ Icc b a, 0 < v - (-a) := fun v hv => by
    have h1 := hv.1; have h2 := hv.2; simp only [Y01.a, Y01.b] at *; linarith [s2lo, s2hi]
  have hy : ∀ v ∈ Icc b a, v - 0 ≠ 0 := fun v hv => ne_of_gt (by
    have h1 := hv.1; have h2 := hv.2; simp only [Y01.a, Y01.b] at *; linarith [s2lo, s2hi])
  exact intIII isP2g hba hxy hx hy (fun v hv => contAt_P2g (mb_ne_zero (hx v hv).ne' (hy v hv)) (mb_ne_one hxy (hy v hv)))

lemma val24q : FIII P2g (-a) 0 a = (2 : ℝ) * D00.Li2 (1 / 2) * lam + (2 : ℝ) * Li3 (1 / 2) + (1 : ℝ) * lam ^ 3 := by
  have e0 : mb (-a) 0 a = (Real.sqrt 2 ^ 2 * a ^ 0 / (Real.sqrt 2 ^ 0 * a ^ 0)) := by unfold mb; have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf)
  have l0 : Real.log (a - (-a)) = (1 : ℝ) * lam + (1 : ℝ) * L :=
    logP 2 1 0 0 _ _ (by norm_num) (by norm_num)
      (by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf))
  have l1 : Real.log (a - 0) = (0 : ℝ) * lam + (1 : ℝ) * L :=
    logP 0 1 0 0 _ _ (by norm_num) (by norm_num)
      (by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf))
  have e1 : ((Real.sqrt 2 ^ 2 * a ^ 0 / (Real.sqrt 2 ^ 0 * a ^ 0)))⁻¹ = 1 / 2 := by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf)
  have l2 : Real.log (1 / 2) = (-1 : ℝ) * lam + (0 : ℝ) * L :=
    logP 0 0 2 0 _ _ (by norm_num) (by norm_num)
      (by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf))
  have l3 : Real.log (1 - (1 / 2)) = (-1 : ℝ) * lam + (0 : ℝ) * L :=
    logP 0 0 2 0 _ _ (by norm_num) (by norm_num)
      (by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf))
  simp only [FIII, e0, l0, l1, P2g, E3, e1, l2, l3]
  try simp only [Real.log_zero, sub_zero, Real.log_one, D00.Li2_zero', Li3_zero, D00.Li2_one']
  ring

lemma val24p : FIII P2g (-a) 0 b = (-1 : ℝ) * L ^ 3 + (-3 / 2 : ℝ) * L ^ 2 * lam + (2 : ℝ) * L * D00.Li2 ((2 - Real.sqrt 2) / 4) + (9 / 4 : ℝ) * L * lam ^ 2 + (3 : ℝ) * D00.Li2 ((2 - Real.sqrt 2) / 4) * lam + (2 : ℝ) * Li3 ((2 - Real.sqrt 2) / 4) + (27 / 8 : ℝ) * lam ^ 3 := by
  have e0 : mb (-a) 0 b = (Real.sqrt 2 ^ 3 * a ^ 1 / (Real.sqrt 2 ^ 0 * a ^ 0)) := by unfold mb; have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf)
  have l0 : Real.log (b - (-a)) = (3 / 2 : ℝ) * lam + (0 : ℝ) * L :=
    logP 3 0 0 0 _ _ (by norm_num) (by norm_num)
      (by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf))
  have l1 : Real.log (b - 0) = (0 : ℝ) * lam + (-1 : ℝ) * L :=
    logP 0 0 0 1 _ _ (by norm_num) (by norm_num)
      (by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf))
  have e1 : ((Real.sqrt 2 ^ 3 * a ^ 1 / (Real.sqrt 2 ^ 0 * a ^ 0)))⁻¹ = (2 - Real.sqrt 2) / 4 := by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf)
  have l2 : Real.log ((2 - Real.sqrt 2) / 4) = (-3 / 2 : ℝ) * lam + (-1 : ℝ) * L :=
    logP 0 0 3 1 _ _ (by norm_num) (by norm_num)
      (by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf))
  have l3 : Real.log (1 - ((2 - Real.sqrt 2) / 4)) = (-3 / 2 : ℝ) * lam + (1 : ℝ) * L :=
    logP 0 1 3 0 _ _ (by norm_num) (by norm_num)
      (by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf))
  simp only [FIII, e0, l0, l1, P2g, E3, e1, l2, l3]
  try simp only [Real.log_zero, sub_zero, Real.log_one, D00.Li2_zero', Li3_zero, D00.Li2_one']
  ring

lemma blk25 : IntervalIntegrable (fpr 0 (-a) (-a)) volume b a ∧
    ∫ v in b..a, fpr 0 (-a) (-a) v = FIV P2a 0 (-a) a - FIV P2a 0 (-a) b := by
  have hxy : (0 : ℝ) - (-a) ≠ 0 := ne_of_gt (by (try simp only [Y01.a, Y01.b]); linarith [s2lo, s2hi])
  have hy : ∀ v ∈ Icc b a, 0 < v - (-a) := fun v hv => by
    have h1 := hv.1; have h2 := hv.2; simp only [Y01.a, Y01.b] at *; linarith [s2lo, s2hi]
  have hx : ∀ v ∈ Ioo b a, v - 0 ≠ 0 := fun v hv => ne_of_gt (by
    have h1 := hv.1; have h2 := hv.2; simp only [Y01.a, Y01.b] at *; linarith [s2lo, s2hi])
  exact intIV isP2a hba hxy hy hx (fun v hv => contAt_P2a (mb_ne_one hxy (hy v hv).ne'))

lemma val25q : FIV P2a 0 (-a) a = (2 : ℝ) * D00.Li2 (1 / 2) * lam + (2 : ℝ) * Li3 (1 / 2) + (1 : ℝ) * lam ^ 3 := by
  have e0 : mb 0 (-a) a = 1 / 2 := by unfold mb; have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf)
  have l0 : Real.log (1 / 2) = (-1 : ℝ) * lam + (0 : ℝ) * L :=
    logP 0 0 2 0 _ _ (by norm_num) (by norm_num)
      (by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf))
  have l1 : Real.log (1 - (1 / 2)) = (-1 : ℝ) * lam + (0 : ℝ) * L :=
    logP 0 0 2 0 _ _ (by norm_num) (by norm_num)
      (by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf))
  simp only [FIV, e0, P2a, l0, l1]
  try simp only [Real.log_zero, sub_zero, Real.log_one, D00.Li2_zero', Li3_zero, D00.Li2_one']
  ring

lemma val25p : FIV P2a 0 (-a) b = (-1 : ℝ) * L ^ 3 + (-3 / 2 : ℝ) * L ^ 2 * lam + (2 : ℝ) * L * D00.Li2 ((2 - Real.sqrt 2) / 4) + (9 / 4 : ℝ) * L * lam ^ 2 + (3 : ℝ) * D00.Li2 ((2 - Real.sqrt 2) / 4) * lam + (2 : ℝ) * Li3 ((2 - Real.sqrt 2) / 4) + (27 / 8 : ℝ) * lam ^ 3 := by
  have e0 : mb 0 (-a) b = (2 - Real.sqrt 2) / 4 := by unfold mb; have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf)
  have l0 : Real.log ((2 - Real.sqrt 2) / 4) = (-3 / 2 : ℝ) * lam + (-1 : ℝ) * L :=
    logP 0 0 3 1 _ _ (by norm_num) (by norm_num)
      (by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf))
  have l1 : Real.log (1 - ((2 - Real.sqrt 2) / 4)) = (-3 / 2 : ℝ) * lam + (1 : ℝ) * L :=
    logP 0 1 3 0 _ _ (by norm_num) (by norm_num)
      (by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf))
  simp only [FIV, e0, P2a, l0, l1]
  try simp only [Real.log_zero, sub_zero, Real.log_one, D00.Li2_zero', Li3_zero, D00.Li2_one']
  ring

lemma blk26 : IntervalIntegrable (fpr b (-a) (-a)) volume 1 a ∧
    ∫ v in 1..a, fpr b (-a) (-a) v = FIV P2a b (-a) a - FIV P2a b (-a) 1 := by
  have hxy : (b : ℝ) - (-a) ≠ 0 := ne_of_gt (by (try simp only [Y01.a, Y01.b]); linarith [s2lo, s2hi])
  have hy : ∀ v ∈ Icc 1 a, 0 < v - (-a) := fun v hv => by
    have h1 := hv.1; have h2 := hv.2; simp only [Y01.a, Y01.b] at *; linarith [s2lo, s2hi]
  have hx : ∀ v ∈ Ioo 1 a, v - b ≠ 0 := fun v hv => ne_of_gt (by
    have h1 := hv.1; have h2 := hv.2; simp only [Y01.a, Y01.b] at *; linarith [s2lo, s2hi])
  exact intIV isP2a h1a hxy hy hx (fun v hv => contAt_P2a (mb_ne_one hxy (hy v hv).ne'))

lemma val26q : FIV P2a b (-a) a = (1 : ℝ) * L ^ 3 + (-1 / 2 : ℝ) * L ^ 2 * lam + (2 : ℝ) * L * D00.Li2 (b) + (2 : ℝ) * Li3 (b) := by
  have e0 : mb b (-a) a = b := by unfold mb; have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf)
  have l0 : Real.log (b) = (0 : ℝ) * lam + (-1 : ℝ) * L :=
    logP 0 0 0 1 _ _ (by norm_num) (by norm_num)
      (by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf))
  have l1 : Real.log (1 - (b)) = (1 / 2 : ℝ) * lam + (-1 : ℝ) * L :=
    logP 1 0 0 1 _ _ (by norm_num) (by norm_num)
      (by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf))
  simp only [FIV, e0, P2a, l0, l1]
  try simp only [Real.log_zero, sub_zero, Real.log_one, D00.Li2_zero', Li3_zero, D00.Li2_one']
  ring

lemma val26p : FIV P2a b (-a) 1 = (4 : ℝ) * L ^ 3 + (-4 : ℝ) * L ^ 2 * lam + (4 : ℝ) * L * D00.Li2 (b ^ 2) + (2 : ℝ) * Li3 (b ^ 2) := by
  have e0 : mb b (-a) 1 = b ^ 2 := by unfold mb; have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf)
  have l0 : Real.log (b ^ 2) = (0 : ℝ) * lam + (-2 : ℝ) * L :=
    logP 0 0 0 2 _ _ (by norm_num) (by norm_num)
      (by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf))
  have l1 : Real.log (1 - (b ^ 2)) = (1 : ℝ) * lam + (-1 : ℝ) * L :=
    logP 2 0 0 1 _ _ (by norm_num) (by norm_num)
      (by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf))
  simp only [FIV, e0, P2a, l0, l1]
  try simp only [Real.log_zero, sub_zero, Real.log_one, D00.Li2_zero', Li3_zero, D00.Li2_one']
  ring

lemma blk27 : IntervalIntegrable (fpr b 0 (-a)) volume 1 a ∧
    ∫ v in 1..a, fpr b 0 (-a) v = FG P1p P2a P2a b 0 (-a) a - FG P1p P2a P2a b 0 (-a) 1 := by
  have hxy : (b : ℝ) - 0 ≠ 0 := ne_of_gt (by (try simp only [Y01.a, Y01.b]); linarith [s2lo, s2hi])
  have hzx : ((-a) : ℝ) - b ≠ 0 := ne_of_lt (by (try simp only [Y01.a, Y01.b]); linarith [s2lo, s2hi])
  have hzy : ((-a) : ℝ) - 0 ≠ 0 := ne_of_lt (by (try simp only [Y01.a, Y01.b]); linarith [s2lo, s2hi])
  have hz : ∀ v ∈ Icc 1 a, 0 < v - (-a) := fun v hv => by
    have h1 := hv.1; have h2 := hv.2; simp only [Y01.a, Y01.b] at *; linarith [s2lo, s2hi]
  have hy : ∀ v ∈ Icc 1 a, v - 0 ≠ 0 := fun v hv => ne_of_gt (by
    have h1 := hv.1; have h2 := hv.2; simp only [Y01.a, Y01.b] at *; linarith [s2lo, s2hi])
  have hx : ∀ v ∈ Ioo 1 a, v - b ≠ 0 := fun v hv => ne_of_gt (by
    have h1 := hv.1; have h2 := hv.2; simp only [Y01.a, Y01.b] at *; linarith [s2lo, s2hi])
  exact intG isP1p isP2a isP2a h1a hxy hzx hzy hz hy hx
    (fun v hv => cPhiM (s1_ne_one hxy hzx (hy v hv) (hz v hv).ne'))
    (fun v hv => contAt_P2a (mb_ne_one hxy (hy v hv)))

lemma val27q : FG P1p P2a P2a b 0 (-a) a = (1 : ℝ) * L ^ 3 + (-3 / 2 : ℝ) * L ^ 2 * lam + (2 : ℝ) * L * D00.Li2 (2 * b) + (2 : ℝ) * L * D00.Li2 (1 - 1 / Real.sqrt 2) + (1 : ℝ) * L * lam ^ 2 + (-2 : ℝ) * D00.Li2 (2 * b) * lam + (-3 : ℝ) * D00.Li2 (1 - 1 / Real.sqrt 2) * lam + (-1 : ℝ) * D00.Li2 (1 / Real.sqrt 2) * lam + (2 : ℝ) * Li3 (2 * b) + (-2 : ℝ) * Li3 (1 / Real.sqrt 2) + (-5 / 4 : ℝ) * lam ^ 3 := by
  have e0 : mb b 0 a = 2 * b := by unfold mb; have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf)
  have e1 : s1 b 0 (-a) a = 1 / Real.sqrt 2 := by unfold s1; have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf)
  have l0 : Real.log (((-a) - b) / ((-a) - 0)) = (3 / 2 : ℝ) * lam + (-1 : ℝ) * L :=
    logP 3 0 0 1 _ _ (by norm_num) (by norm_num)
      (by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf))
  have e2 : (1:ℝ) - (1 / Real.sqrt 2) = 1 - 1 / Real.sqrt 2 := by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf)
  have l1 : Real.log (1 / Real.sqrt 2) = (-1 / 2 : ℝ) * lam + (0 : ℝ) * L :=
    logP 0 0 1 0 _ _ (by norm_num) (by norm_num)
      (by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf))
  have l2 : Real.log (1 - 1 / Real.sqrt 2) = (-1 / 2 : ℝ) * lam + (-1 : ℝ) * L :=
    logP 0 0 1 1 _ _ (by norm_num) (by norm_num)
      (by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf))
  have l3 : Real.log (2 * b) = (1 : ℝ) * lam + (-1 : ℝ) * L :=
    logP 2 0 0 1 _ _ (by norm_num) (by norm_num)
      (by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf))
  have l4 : Real.log (1 - (2 * b)) = (0 : ℝ) * lam + (-2 : ℝ) * L :=
    logP 0 0 0 2 _ _ (by norm_num) (by norm_num)
      (by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf))
  simp only [FG, e0, e1, l0, Phi, P1p, P2a, e2, l1, l2, l3, l4]
  try simp only [Real.log_zero, sub_zero, Real.log_one, D00.Li2_zero', Li3_zero, D00.Li2_one']
  ring

lemma val27p : FG P1p P2a P2a b 0 (-a) 1 = (1 : ℝ) * L ^ 3 + (-2 : ℝ) * L ^ 2 * lam + (2 : ℝ) * L * D00.Li2 (1 / 2) + (2 : ℝ) * L * D00.Li2 (2 - Real.sqrt 2) + (13 / 4 : ℝ) * L * lam ^ 2 + (-5 : ℝ) * D00.Li2 (1 / 2) * lam + (-1 : ℝ) * D00.Li2 (2 - Real.sqrt 2) * lam + (-2 : ℝ) * Li3 (1 / 2) + (2 : ℝ) * Li3 (2 - Real.sqrt 2) + (-13 / 4 : ℝ) * lam ^ 3 := by
  have e0 : mb b 0 1 = 2 - Real.sqrt 2 := by unfold mb; have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf)
  have e1 : s1 b 0 (-a) 1 = 1 / 2 := by unfold s1; have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf)
  have l0 : Real.log (((-a) - b) / ((-a) - 0)) = (3 / 2 : ℝ) * lam + (-1 : ℝ) * L :=
    logP 3 0 0 1 _ _ (by norm_num) (by norm_num)
      (by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf))
  have e2 : (1:ℝ) - (1 / 2) = 1 / 2 := by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf)
  have l1 : Real.log (1 / 2) = (-1 : ℝ) * lam + (0 : ℝ) * L :=
    logP 0 0 2 0 _ _ (by norm_num) (by norm_num)
      (by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf))
  have l2 : Real.log (2 - Real.sqrt 2) = (1 / 2 : ℝ) * lam + (-1 : ℝ) * L :=
    logP 1 0 0 1 _ _ (by norm_num) (by norm_num)
      (by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf))
  have l3 : Real.log (1 - (2 - Real.sqrt 2)) = (0 : ℝ) * lam + (-1 : ℝ) * L :=
    logP 0 0 0 1 _ _ (by norm_num) (by norm_num)
      (by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf))
  simp only [FG, e0, e1, l0, Phi, P1p, P2a, e2, l1, l2, l3]
  try simp only [Real.log_zero, sub_zero, Real.log_one, D00.Li2_zero', Li3_zero, D00.Li2_one']
  ring

lemma blk28 : IntervalIntegrable (fpr (-a) b (-a)) volume 1 a ∧
    ∫ v in 1..a, fpr (-a) b (-a) v = FIII P2g (-a) b a - FIII P2g (-a) b 1 := by
  have hxy : ((-a) : ℝ) - b ≠ 0 := ne_of_lt (by (try simp only [Y01.a, Y01.b]); linarith [s2lo, s2hi])
  have hx : ∀ v ∈ Icc 1 a, 0 < v - (-a) := fun v hv => by
    have h1 := hv.1; have h2 := hv.2; simp only [Y01.a, Y01.b] at *; linarith [s2lo, s2hi]
  have hy : ∀ v ∈ Icc 1 a, v - b ≠ 0 := fun v hv => ne_of_gt (by
    have h1 := hv.1; have h2 := hv.2; simp only [Y01.a, Y01.b] at *; linarith [s2lo, s2hi])
  exact intIII isP2g h1a hxy hx hy (fun v hv => contAt_P2g (mb_ne_zero (hx v hv).ne' (hy v hv)) (mb_ne_one hxy (hy v hv)))

lemma val28q : FIII P2g (-a) b a = (1 : ℝ) * L ^ 3 + (-1 / 2 : ℝ) * L ^ 2 * lam + (2 : ℝ) * L * D00.Li2 (b) + (2 : ℝ) * Li3 (b) := by
  have e0 : mb (-a) b a = (Real.sqrt 2 ^ 0 * a ^ 1 / (Real.sqrt 2 ^ 0 * a ^ 0)) := by unfold mb; have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf)
  have l0 : Real.log (a - (-a)) = (1 : ℝ) * lam + (1 : ℝ) * L :=
    logP 2 1 0 0 _ _ (by norm_num) (by norm_num)
      (by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf))
  have l1 : Real.log (a - b) = (1 : ℝ) * lam + (0 : ℝ) * L :=
    logP 2 0 0 0 _ _ (by norm_num) (by norm_num)
      (by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf))
  have e1 : ((Real.sqrt 2 ^ 0 * a ^ 1 / (Real.sqrt 2 ^ 0 * a ^ 0)))⁻¹ = b := by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf)
  have l2 : Real.log (b) = (0 : ℝ) * lam + (-1 : ℝ) * L :=
    logP 0 0 0 1 _ _ (by norm_num) (by norm_num)
      (by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf))
  have l3 : Real.log (1 - (b)) = (1 / 2 : ℝ) * lam + (-1 : ℝ) * L :=
    logP 1 0 0 1 _ _ (by norm_num) (by norm_num)
      (by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf))
  simp only [FIII, e0, l0, l1, P2g, E3, e1, l2, l3]
  try simp only [Real.log_zero, sub_zero, Real.log_one, D00.Li2_zero', Li3_zero, D00.Li2_one']
  ring

lemma val28p : FIII P2g (-a) b 1 = (4 : ℝ) * L ^ 3 + (-4 : ℝ) * L ^ 2 * lam + (4 : ℝ) * L * D00.Li2 (b ^ 2) + (2 : ℝ) * Li3 (b ^ 2) := by
  have e0 : mb (-a) b 1 = (Real.sqrt 2 ^ 0 * a ^ 2 / (Real.sqrt 2 ^ 0 * a ^ 0)) := by unfold mb; have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf)
  have l0 : Real.log (1 - (-a)) = (1 / 2 : ℝ) * lam + (1 : ℝ) * L :=
    logP 1 1 0 0 _ _ (by norm_num) (by norm_num)
      (by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf))
  have l1 : Real.log (1 - b) = (1 / 2 : ℝ) * lam + (-1 : ℝ) * L :=
    logP 1 0 0 1 _ _ (by norm_num) (by norm_num)
      (by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf))
  have e1 : ((Real.sqrt 2 ^ 0 * a ^ 2 / (Real.sqrt 2 ^ 0 * a ^ 0)))⁻¹ = b ^ 2 := by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf)
  have l2 : Real.log (b ^ 2) = (0 : ℝ) * lam + (-2 : ℝ) * L :=
    logP 0 0 0 2 _ _ (by norm_num) (by norm_num)
      (by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf))
  have l3 : Real.log (1 - (b ^ 2)) = (1 : ℝ) * lam + (-1 : ℝ) * L :=
    logP 2 0 0 1 _ _ (by norm_num) (by norm_num)
      (by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf))
  simp only [FIII, e0, l0, l1, P2g, E3, e1, l2, l3]
  try simp only [Real.log_zero, sub_zero, Real.log_one, D00.Li2_zero', Li3_zero, D00.Li2_one']
  ring

lemma blk29 : IntervalIntegrable (fpr 0 b (-a)) volume 1 a ∧
    ∫ v in 1..a, fpr 0 b (-a) v = FG P1p P2g P2g 0 b (-a) a - FG P1p P2g P2g 0 b (-a) 1 := by
  have hxy : (0 : ℝ) - b ≠ 0 := ne_of_lt (by (try simp only [Y01.a, Y01.b]); linarith [s2lo, s2hi])
  have hzx : ((-a) : ℝ) - 0 ≠ 0 := ne_of_lt (by (try simp only [Y01.a, Y01.b]); linarith [s2lo, s2hi])
  have hzy : ((-a) : ℝ) - b ≠ 0 := ne_of_lt (by (try simp only [Y01.a, Y01.b]); linarith [s2lo, s2hi])
  have hz : ∀ v ∈ Icc 1 a, 0 < v - (-a) := fun v hv => by
    have h1 := hv.1; have h2 := hv.2; simp only [Y01.a, Y01.b] at *; linarith [s2lo, s2hi]
  have hy : ∀ v ∈ Icc 1 a, v - b ≠ 0 := fun v hv => ne_of_gt (by
    have h1 := hv.1; have h2 := hv.2; simp only [Y01.a, Y01.b] at *; linarith [s2lo, s2hi])
  have hx : ∀ v ∈ Ioo 1 a, v - 0 ≠ 0 := fun v hv => ne_of_gt (by
    have h1 := hv.1; have h2 := hv.2; simp only [Y01.a, Y01.b] at *; linarith [s2lo, s2hi])
  have hx2 : ∀ v ∈ Icc 1 a, v - 0 ≠ 0 := fun v hv => ne_of_gt (by
    have h1 := hv.1; have h2 := hv.2; simp only [Y01.a, Y01.b] at *; linarith [s2lo, s2hi])
  exact intG isP1p isP2g isP2g h1a hxy hzx hzy hz hy hx
    (fun v hv => cPhiG (s1_ne_zero hzx hzy (hx2 v hv) (hy v hv)) (s1_ne_one hxy hzx (hy v hv) (hz v hv).ne'))
    (fun v hv => contAt_P2g (mb_ne_zero (hx2 v hv) (hy v hv)) (mb_ne_one hxy (hy v hv)))

lemma val29q : FG P1p P2g P2g 0 b (-a) a = (2 / 3 : ℝ) * L ^ 3 + (2 : ℝ) * L * D00.Li2 (2 * b) + (-2 : ℝ) * L * D00.Li2 (-b) + (-3 / 2 : ℝ) * L * lam ^ 2 + (-2 : ℝ) * D00.Li2 (2 * b) * lam + (3 : ℝ) * D00.Li2 (-b) * lam + (-1 : ℝ) * D00.Li2 (1 / Real.sqrt 2) * lam + (2 : ℝ) * Li3 (2 * b) + (-2 : ℝ) * Li3 (1 / Real.sqrt 2) + (1 / 4 : ℝ) * lam ^ 3 := by
  have e0 : mb 0 b a = (Real.sqrt 2 ^ 0 * a ^ 1 / (Real.sqrt 2 ^ 2 * a ^ 0)) := by unfold mb; have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf)
  have e1 : s1 0 b (-a) a = (Real.sqrt 2 ^ 1 * a ^ 0 / (Real.sqrt 2 ^ 0 * a ^ 0)) := by unfold s1; have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf)
  have l0 : Real.log (((-a) - 0) / ((-a) - b)) = (-3 / 2 : ℝ) * lam + (1 : ℝ) * L :=
    logP 0 1 3 0 _ _ (by norm_num) (by norm_num)
      (by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf))
  have e2 : (1:ℝ) - ((Real.sqrt 2 ^ 1 * a ^ 0 / (Real.sqrt 2 ^ 0 * a ^ 0))) = -b := by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf)
  have e3 : ((Real.sqrt 2 ^ 1 * a ^ 0 / (Real.sqrt 2 ^ 0 * a ^ 0)))⁻¹ = 1 / Real.sqrt 2 := by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf)
  have l1 : Real.log (-b) = (0 : ℝ) * lam + (-1 : ℝ) * L :=
    logN 0 0 0 1 _ _ (by norm_num) (by norm_num)
      (by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf))
  have l2 : Real.log (1 / Real.sqrt 2) = (-1 / 2 : ℝ) * lam + (0 : ℝ) * L :=
    logP 0 0 1 0 _ _ (by norm_num) (by norm_num)
      (by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf))
  have l3 : Real.log (1 - (1 / Real.sqrt 2)) = (-1 / 2 : ℝ) * lam + (-1 : ℝ) * L :=
    logP 0 0 1 1 _ _ (by norm_num) (by norm_num)
      (by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf))
  have e4 : ((Real.sqrt 2 ^ 0 * a ^ 1 / (Real.sqrt 2 ^ 2 * a ^ 0)))⁻¹ = 2 * b := by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf)
  have l4 : Real.log (2 * b) = (1 : ℝ) * lam + (-1 : ℝ) * L :=
    logP 2 0 0 1 _ _ (by norm_num) (by norm_num)
      (by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf))
  have l5 : Real.log (1 - (2 * b)) = (0 : ℝ) * lam + (-2 : ℝ) * L :=
    logP 0 0 0 2 _ _ (by norm_num) (by norm_num)
      (by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf))
  simp only [FG, e0, e1, l0, Phi, P1p, P2g, e2, e3, E3, l1, l2, l3, e4, l4, l5]
  try simp only [Real.log_zero, sub_zero, Real.log_one, D00.Li2_zero', Li3_zero, D00.Li2_one']
  ring

lemma val29p : FG P1p P2g P2g 0 b (-a) 1 = (2 / 3 : ℝ) * L ^ 3 + (-1 / 2 : ℝ) * L ^ 2 * lam + (-2 : ℝ) * L * D00.Li2 (-1) + (2 : ℝ) * L * D00.Li2 (2 - Real.sqrt 2) + (3 : ℝ) * D00.Li2 (-1) * lam + (-2 : ℝ) * D00.Li2 (1 / 2) * lam + (-1 : ℝ) * D00.Li2 (2 - Real.sqrt 2) * lam + (-2 : ℝ) * Li3 (1 / 2) + (2 : ℝ) * Li3 (2 - Real.sqrt 2) + (-5 / 8 : ℝ) * lam ^ 3 := by
  have e0 : mb 0 b 1 = (Real.sqrt 2 ^ 0 * a ^ 1 / (Real.sqrt 2 ^ 1 * a ^ 0)) := by unfold mb; have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf)
  have e1 : s1 0 b (-a) 1 = (Real.sqrt 2 ^ 2 * a ^ 0 / (Real.sqrt 2 ^ 0 * a ^ 0)) := by unfold s1; have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf)
  have l0 : Real.log (((-a) - 0) / ((-a) - b)) = (-3 / 2 : ℝ) * lam + (1 : ℝ) * L :=
    logP 0 1 3 0 _ _ (by norm_num) (by norm_num)
      (by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf))
  have e2 : (1:ℝ) - ((Real.sqrt 2 ^ 2 * a ^ 0 / (Real.sqrt 2 ^ 0 * a ^ 0))) = -1 := by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf)
  have e3 : ((Real.sqrt 2 ^ 2 * a ^ 0 / (Real.sqrt 2 ^ 0 * a ^ 0)))⁻¹ = 1 / 2 := by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf)
  have l1 : Real.log (-1) = (0 : ℝ) * lam + (0 : ℝ) * L :=
    logN 0 0 0 0 _ _ (by norm_num) (by norm_num)
      (by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf))
  have l2 : Real.log (1 / 2) = (-1 : ℝ) * lam + (0 : ℝ) * L :=
    logP 0 0 2 0 _ _ (by norm_num) (by norm_num)
      (by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf))
  have l3 : Real.log (1 - (1 / 2)) = (-1 : ℝ) * lam + (0 : ℝ) * L :=
    logP 0 0 2 0 _ _ (by norm_num) (by norm_num)
      (by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf))
  have e4 : ((Real.sqrt 2 ^ 0 * a ^ 1 / (Real.sqrt 2 ^ 1 * a ^ 0)))⁻¹ = 2 - Real.sqrt 2 := by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf)
  have l4 : Real.log (2 - Real.sqrt 2) = (1 / 2 : ℝ) * lam + (-1 : ℝ) * L :=
    logP 1 0 0 1 _ _ (by norm_num) (by norm_num)
      (by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf))
  have l5 : Real.log (1 - (2 - Real.sqrt 2)) = (0 : ℝ) * lam + (-1 : ℝ) * L :=
    logP 0 0 0 1 _ _ (by norm_num) (by norm_num)
      (by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf))
  simp only [FG, e0, e1, l0, Phi, P1p, P2g, e2, e3, E3, l1, l2, l3, e4, l4, l5]
  try simp only [Real.log_zero, sub_zero, Real.log_one, D00.Li2_zero', Li3_zero, D00.Li2_one']
  ring

lemma blk30 : IntervalIntegrable (flin a 0) volume b a ∧
    ∫ v in b..a, flin a 0 v = FL P1p a 0 a - FL P1p a 0 b := by
  have hxz : (0 : ℝ) - a ≠ 0 := ne_of_lt (by (try simp only [Y01.a, Y01.b]); linarith [s2lo, s2hi])
  have hz : ∀ v ∈ Icc b a, 0 < v - 0 := fun v hv => by
    have h1 := hv.1; have h2 := hv.2; simp only [Y01.a, Y01.b] at *; linarith [s2lo, s2hi]
  have hx : ∀ v ∈ Ioo b a, v - a ≠ 0 := fun v hv => ne_of_lt (by
    have h1 := hv.1; have h2 := hv.2; simp only [Y01.a, Y01.b] at *; linarith [s2lo, s2hi])
  exact intL isP1p hba hxz hz hx (fun v hv => cPsiM (sA_ne_one hxz (hz v hv).ne'))

lemma val30q : FL P1p a 0 a = (-1 / 6 : ℝ) * (Real.pi ^ 2) := by
  have e0 : sA a 0 a = 0 := by unfold sA; have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf)
  have l0 : Real.log (0 - a) = (0 : ℝ) * lam + (1 : ℝ) * L :=
    logN 0 1 0 0 _ _ (by norm_num) (by norm_num)
      (by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf))
  simp only [FL, e0, l0, Psi, P1p]
  try simp only [Real.log_zero, sub_zero, Real.log_one, D00.Li2_zero', Li3_zero, D00.Li2_one']
  ring

lemma val30p : FL P1p a 0 b = (-2 : ℝ) * L ^ 2 + (-1 : ℝ) * D00.Li2 (b ^ 2) := by
  have e0 : sA a 0 b = 2 * b := by unfold sA; have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf)
  have l0 : Real.log (0 - a) = (0 : ℝ) * lam + (1 : ℝ) * L :=
    logN 0 1 0 0 _ _ (by norm_num) (by norm_num)
      (by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf))
  have e1 : (1:ℝ) - (2 * b) = b ^ 2 := by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf)
  have l1 : Real.log (b ^ 2) = (0 : ℝ) * lam + (-2 : ℝ) * L :=
    logP 0 0 0 2 _ _ (by norm_num) (by norm_num)
      (by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf))
  simp only [FL, e0, l0, Psi, P1p, e1, l1]
  try simp only [Real.log_zero, sub_zero, Real.log_one, D00.Li2_zero', Li3_zero, D00.Li2_one']
  ring

lemma blk31 : IntervalIntegrable (flin (-b) 0) volume b a ∧
    ∫ v in b..a, flin (-b) 0 v = FL P1p (-b) 0 a - FL P1p (-b) 0 b := by
  have hxz : (0 : ℝ) - (-b) ≠ 0 := ne_of_gt (by (try simp only [Y01.a, Y01.b]); linarith [s2lo, s2hi])
  have hz : ∀ v ∈ Icc b a, 0 < v - 0 := fun v hv => by
    have h1 := hv.1; have h2 := hv.2; simp only [Y01.a, Y01.b] at *; linarith [s2lo, s2hi]
  have hx : ∀ v ∈ Ioo b a, v - (-b) ≠ 0 := fun v hv => ne_of_gt (by
    have h1 := hv.1; have h2 := hv.2; simp only [Y01.a, Y01.b] at *; linarith [s2lo, s2hi])
  exact intL isP1p hba hxz hz hx (fun v hv => cPsiM (sA_ne_one hxz (hz v hv).ne'))

lemma val31q : FL P1p (-b) 0 a = (-2 : ℝ) * L ^ 2 + (-1 : ℝ) * D00.Li2 (-(a ^ 2)) := by
  have e0 : sA (-b) 0 a = (Real.sqrt 2 ^ 3 * a ^ 1 / (Real.sqrt 2 ^ 0 * a ^ 0)) := by unfold sA; have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf)
  have l0 : Real.log (0 - (-b)) = (0 : ℝ) * lam + (-1 : ℝ) * L :=
    logP 0 0 0 1 _ _ (by norm_num) (by norm_num)
      (by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf))
  have e1 : (1:ℝ) - ((Real.sqrt 2 ^ 3 * a ^ 1 / (Real.sqrt 2 ^ 0 * a ^ 0))) = -(a ^ 2) := by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf)
  have l1 : Real.log (-(a ^ 2)) = (0 : ℝ) * lam + (2 : ℝ) * L :=
    logN 0 2 0 0 _ _ (by norm_num) (by norm_num)
      (by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf))
  simp only [FL, e0, l0, Psi, P1p, e1, l1]
  try simp only [Real.log_zero, sub_zero, Real.log_one, D00.Li2_zero', Li3_zero, D00.Li2_one']
  ring

lemma val31p : FL P1p (-b) 0 b = (-1 : ℝ) * D00.Li2 (-1) := by
  have e0 : sA (-b) 0 b = (Real.sqrt 2 ^ 2 * a ^ 0 / (Real.sqrt 2 ^ 0 * a ^ 0)) := by unfold sA; have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf)
  have l0 : Real.log (0 - (-b)) = (0 : ℝ) * lam + (-1 : ℝ) * L :=
    logP 0 0 0 1 _ _ (by norm_num) (by norm_num)
      (by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf))
  have e1 : (1:ℝ) - ((Real.sqrt 2 ^ 2 * a ^ 0 / (Real.sqrt 2 ^ 0 * a ^ 0))) = -1 := by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf)
  have l1 : Real.log (-1) = (0 : ℝ) * lam + (0 : ℝ) * L :=
    logN 0 0 0 0 _ _ (by norm_num) (by norm_num)
      (by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf))
  simp only [FL, e0, l0, Psi, P1p, e1, l1]
  try simp only [Real.log_zero, sub_zero, Real.log_one, D00.Li2_zero', Li3_zero, D00.Li2_one']
  ring

lemma blk32 : IntervalIntegrable (fsq a 0) volume b a ∧
    ∫ v in b..a, fsq a 0 v = FA P1p P2a a 0 a - FA P1p P2a a 0 b := by
  have hxz : (0 : ℝ) - a ≠ 0 := ne_of_lt (by (try simp only [Y01.a, Y01.b]); linarith [s2lo, s2hi])
  have hz : ∀ v ∈ Icc b a, 0 < v - 0 := fun v hv => by
    have h1 := hv.1; have h2 := hv.2; simp only [Y01.a, Y01.b] at *; linarith [s2lo, s2hi]
  have hx : ∀ v ∈ Ioo b a, v - a ≠ 0 := fun v hv => ne_of_lt (by
    have h1 := hv.1; have h2 := hv.2; simp only [Y01.a, Y01.b] at *; linarith [s2lo, s2hi])
  exact intA isP1p isP2a hba hxz hz hx (fun v hv => cPhiM (sA_ne_one hxz (hz v hv).ne'))

lemma val32q : FA P1p P2a a 0 a = (-1 / 3 : ℝ) * L * (Real.pi ^ 2) := by
  have e0 : sA a 0 a = 0 := by unfold sA; have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf)
  have l0 : Real.log (0 - a) = (0 : ℝ) * lam + (1 : ℝ) * L :=
    logN 0 1 0 0 _ _ (by norm_num) (by norm_num)
      (by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf))
  simp only [FA, e0, l0, Phi, P1p, P2a]
  try simp only [Real.log_zero, sub_zero, Real.log_one, D00.Li2_zero', Li3_zero, D00.Li2_one']
  ring

lemma val32p : FA P1p P2a a 0 b = (-4 : ℝ) * L ^ 3 + (4 : ℝ) * L ^ 2 * lam + (-2 : ℝ) * L * D00.Li2 (2 * b) + (-2 : ℝ) * L * D00.Li2 (b ^ 2) + (-2 : ℝ) * L * lam ^ 2 + (2 : ℝ) * D00.Li2 (2 * b) * lam + (-2 : ℝ) * Li3 (2 * b) := by
  have e0 : sA a 0 b = 2 * b := by unfold sA; have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf)
  have l0 : Real.log (0 - a) = (0 : ℝ) * lam + (1 : ℝ) * L :=
    logN 0 1 0 0 _ _ (by norm_num) (by norm_num)
      (by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf))
  have e1 : (1:ℝ) - (2 * b) = b ^ 2 := by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf)
  have l1 : Real.log (2 * b) = (1 : ℝ) * lam + (-1 : ℝ) * L :=
    logP 2 0 0 1 _ _ (by norm_num) (by norm_num)
      (by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf))
  have l2 : Real.log (b ^ 2) = (0 : ℝ) * lam + (-2 : ℝ) * L :=
    logP 0 0 0 2 _ _ (by norm_num) (by norm_num)
      (by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf))
  simp only [FA, e0, l0, Phi, P1p, P2a, e1, l1, l2]
  try simp only [Real.log_zero, sub_zero, Real.log_one, D00.Li2_zero', Li3_zero, D00.Li2_one']
  ring

lemma blk33 : IntervalIntegrable (fpr b a 0) volume b 1 ∧
    ∫ v in b..1, fpr b a 0 v = FG P1n P2a P2a b a 0 1 - FG P1n P2a P2a b a 0 b := by
  have hxy : (b : ℝ) - a ≠ 0 := ne_of_lt (by (try simp only [Y01.a, Y01.b]); linarith [s2lo, s2hi])
  have hzx : (0 : ℝ) - b ≠ 0 := ne_of_lt (by (try simp only [Y01.a, Y01.b]); linarith [s2lo, s2hi])
  have hzy : (0 : ℝ) - a ≠ 0 := ne_of_lt (by (try simp only [Y01.a, Y01.b]); linarith [s2lo, s2hi])
  have hz : ∀ v ∈ Icc b 1, 0 < v - 0 := fun v hv => by
    have h1 := hv.1; have h2 := hv.2; simp only [Y01.a, Y01.b] at *; linarith [s2lo, s2hi]
  have hy : ∀ v ∈ Icc b 1, v - a ≠ 0 := fun v hv => ne_of_lt (by
    have h1 := hv.1; have h2 := hv.2; simp only [Y01.a, Y01.b] at *; linarith [s2lo, s2hi])
  have hx : ∀ v ∈ Ioo b 1, v - b ≠ 0 := fun v hv => ne_of_gt (by
    have h1 := hv.1; have h2 := hv.2; simp only [Y01.a, Y01.b] at *; linarith [s2lo, s2hi])
  exact intG isP1n isP2a isP2a hb1 hxy hzx hzy hz hy hx
    (fun v hv => cPhiN (s1_ne_one hxy hzx (hy v hv) (hz v hv).ne'))
    (fun v hv => contAt_P2a (mb_ne_one hxy (hy v hv)))

lemma val33q : FG P1n P2a P2a b a 0 1 = (1 : ℝ) * L ^ 3 + (-2 : ℝ) * L * D00.Li2 (-a) + (2 : ℝ) * L * D00.Li2 (-b) + (-2 : ℝ) * Li3 (-a) + (2 : ℝ) * Li3 (-b) := by
  have e0 : mb b a 1 = -b := by unfold mb; have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf)
  have e1 : s1 b a 0 1 = -a := by unfold s1; have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf)
  have l0 : Real.log ((0 - b) / (0 - a)) = (0 : ℝ) * lam + (-2 : ℝ) * L :=
    logP 0 0 0 2 _ _ (by norm_num) (by norm_num)
      (by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf))
  have l1 : Real.log (-a) = (0 : ℝ) * lam + (1 : ℝ) * L :=
    logN 0 1 0 0 _ _ (by norm_num) (by norm_num)
      (by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf))
  have l2 : Real.log (1 - (-a)) = (1 / 2 : ℝ) * lam + (1 : ℝ) * L :=
    logP 1 1 0 0 _ _ (by norm_num) (by norm_num)
      (by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf))
  have l3 : Real.log (-b) = (0 : ℝ) * lam + (-1 : ℝ) * L :=
    logN 0 0 0 1 _ _ (by norm_num) (by norm_num)
      (by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf))
  have l4 : Real.log (1 - (-b)) = (1 / 2 : ℝ) * lam + (0 : ℝ) * L :=
    logP 1 0 0 0 _ _ (by norm_num) (by norm_num)
      (by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf))
  simp only [FG, e0, e1, l0, Phi, P1n, P2a, l1, l2, l3, l4]
  try simp only [Real.log_zero, sub_zero, Real.log_one, D00.Li2_zero', Li3_zero, D00.Li2_one']
  ring

lemma val33p : FG P1n P2a P2a b a 0 b = 0 := by
  have e0 : mb b a b = 0 := by unfold mb; have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf)
  have e1 : s1 b a 0 b = 0 := by unfold s1; have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf)
  have l0 : Real.log ((0 - b) / (0 - a)) = (0 : ℝ) * lam + (-2 : ℝ) * L :=
    logP 0 0 0 2 _ _ (by norm_num) (by norm_num)
      (by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf))
  simp only [FG, e0, e1, l0, Phi, P1n, P2a]
  try simp only [Real.log_zero, sub_zero, Real.log_one, D00.Li2_zero', Li3_zero, D00.Li2_one']
  ring

lemma blk34 : IntervalIntegrable (fsq (-b) 0) volume b a ∧
    ∫ v in b..a, fsq (-b) 0 v = FA P1p P2g (-b) 0 a - FA P1p P2g (-b) 0 b := by
  have hxz : (0 : ℝ) - (-b) ≠ 0 := ne_of_gt (by (try simp only [Y01.a, Y01.b]); linarith [s2lo, s2hi])
  have hz : ∀ v ∈ Icc b a, 0 < v - 0 := fun v hv => by
    have h1 := hv.1; have h2 := hv.2; simp only [Y01.a, Y01.b] at *; linarith [s2lo, s2hi]
  have hx : ∀ v ∈ Ioo b a, v - (-b) ≠ 0 := fun v hv => ne_of_gt (by
    have h1 := hv.1; have h2 := hv.2; simp only [Y01.a, Y01.b] at *; linarith [s2lo, s2hi])
  have hx2 : ∀ v ∈ Icc b a, v - (-b) ≠ 0 := fun v hv => ne_of_gt (by
    have h1 := hv.1; have h2 := hv.2; simp only [Y01.a, Y01.b] at *; linarith [s2lo, s2hi])
  exact intA isP1p isP2g hba hxz hz hx (fun v hv => cPhiG (sA_ne_zero hxz (hx2 v hv)) (sA_ne_one hxz (hz v hv).ne'))

lemma val34q : FA P1p P2g (-b) 0 a = (10 / 3 : ℝ) * L ^ 3 + (3 : ℝ) * L ^ 2 * lam + (2 : ℝ) * L * D00.Li2 (-(a ^ 2)) + (-2 : ℝ) * L * D00.Li2 ((2 - Real.sqrt 2) / 4) + (-3 : ℝ) * D00.Li2 ((2 - Real.sqrt 2) / 4) * lam + (-2 : ℝ) * Li3 ((2 - Real.sqrt 2) / 4) + (-9 / 4 : ℝ) * lam ^ 3 := by
  have e0 : sA (-b) 0 a = (Real.sqrt 2 ^ 3 * a ^ 1 / (Real.sqrt 2 ^ 0 * a ^ 0)) := by unfold sA; have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf)
  have l0 : Real.log (0 - (-b)) = (0 : ℝ) * lam + (-1 : ℝ) * L :=
    logP 0 0 0 1 _ _ (by norm_num) (by norm_num)
      (by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf))
  have e1 : (1:ℝ) - ((Real.sqrt 2 ^ 3 * a ^ 1 / (Real.sqrt 2 ^ 0 * a ^ 0))) = -(a ^ 2) := by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf)
  have e2 : ((Real.sqrt 2 ^ 3 * a ^ 1 / (Real.sqrt 2 ^ 0 * a ^ 0)))⁻¹ = (2 - Real.sqrt 2) / 4 := by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf)
  have l1 : Real.log (-(a ^ 2)) = (0 : ℝ) * lam + (2 : ℝ) * L :=
    logN 0 2 0 0 _ _ (by norm_num) (by norm_num)
      (by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf))
  have l2 : Real.log ((2 - Real.sqrt 2) / 4) = (-3 / 2 : ℝ) * lam + (-1 : ℝ) * L :=
    logP 0 0 3 1 _ _ (by norm_num) (by norm_num)
      (by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf))
  have l3 : Real.log (1 - ((2 - Real.sqrt 2) / 4)) = (-3 / 2 : ℝ) * lam + (1 : ℝ) * L :=
    logP 0 1 3 0 _ _ (by norm_num) (by norm_num)
      (by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf))
  simp only [FA, e0, l0, Phi, P1p, P2g, e1, e2, E3, l1, l2, l3]
  try simp only [Real.log_zero, sub_zero, Real.log_one, D00.Li2_zero', Li3_zero, D00.Li2_one']
  ring

lemma val34p : FA P1p P2g (-b) 0 b = (2 : ℝ) * L * D00.Li2 (-1) + (-2 : ℝ) * D00.Li2 (1 / 2) * lam + (-2 : ℝ) * Li3 (1 / 2) + (-2 / 3 : ℝ) * lam ^ 3 := by
  have e0 : sA (-b) 0 b = (Real.sqrt 2 ^ 2 * a ^ 0 / (Real.sqrt 2 ^ 0 * a ^ 0)) := by unfold sA; have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf)
  have l0 : Real.log (0 - (-b)) = (0 : ℝ) * lam + (-1 : ℝ) * L :=
    logP 0 0 0 1 _ _ (by norm_num) (by norm_num)
      (by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf))
  have e1 : (1:ℝ) - ((Real.sqrt 2 ^ 2 * a ^ 0 / (Real.sqrt 2 ^ 0 * a ^ 0))) = -1 := by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf)
  have e2 : ((Real.sqrt 2 ^ 2 * a ^ 0 / (Real.sqrt 2 ^ 0 * a ^ 0)))⁻¹ = 1 / 2 := by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf)
  have l1 : Real.log (-1) = (0 : ℝ) * lam + (0 : ℝ) * L :=
    logN 0 0 0 0 _ _ (by norm_num) (by norm_num)
      (by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf))
  have l2 : Real.log (1 / 2) = (-1 : ℝ) * lam + (0 : ℝ) * L :=
    logP 0 0 2 0 _ _ (by norm_num) (by norm_num)
      (by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf))
  have l3 : Real.log (1 - (1 / 2)) = (-1 : ℝ) * lam + (0 : ℝ) * L :=
    logP 0 0 2 0 _ _ (by norm_num) (by norm_num)
      (by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf))
  simp only [FA, e0, l0, Phi, P1p, P2g, e1, e2, E3, l1, l2, l3]
  try simp only [Real.log_zero, sub_zero, Real.log_one, D00.Li2_zero', Li3_zero, D00.Li2_one']
  ring

lemma blk35 : IntervalIntegrable (fpr b (-b) 0) volume b a ∧
    ∫ v in b..a, fpr b (-b) 0 v = FG P1n P2a P2a b (-b) 0 a - FG P1n P2a P2a b (-b) 0 b := by
  have hxy : (b : ℝ) - (-b) ≠ 0 := ne_of_gt (by (try simp only [Y01.a, Y01.b]); linarith [s2lo, s2hi])
  have hzx : (0 : ℝ) - b ≠ 0 := ne_of_lt (by (try simp only [Y01.a, Y01.b]); linarith [s2lo, s2hi])
  have hzy : (0 : ℝ) - (-b) ≠ 0 := ne_of_gt (by (try simp only [Y01.a, Y01.b]); linarith [s2lo, s2hi])
  have hz : ∀ v ∈ Icc b a, 0 < v - 0 := fun v hv => by
    have h1 := hv.1; have h2 := hv.2; simp only [Y01.a, Y01.b] at *; linarith [s2lo, s2hi]
  have hy : ∀ v ∈ Icc b a, v - (-b) ≠ 0 := fun v hv => ne_of_gt (by
    have h1 := hv.1; have h2 := hv.2; simp only [Y01.a, Y01.b] at *; linarith [s2lo, s2hi])
  have hx : ∀ v ∈ Ioo b a, v - b ≠ 0 := fun v hv => ne_of_gt (by
    have h1 := hv.1; have h2 := hv.2; simp only [Y01.a, Y01.b] at *; linarith [s2lo, s2hi])
  exact intG isP1n isP2a isP2a hba hxy hzx hzy hz hy hx
    (fun v hv => cPhiN (s1_ne_one hxy hzx (hy v hv) (hz v hv).ne'))
    (fun v hv => contAt_P2a (mb_ne_one hxy (hy v hv)))

lemma val35q : FG P1n P2a P2a b (-b) 0 a = (1 / 2 : ℝ) * L * lam ^ 2 + (-1 : ℝ) * D00.Li2 (-(1 / Real.sqrt 2)) * lam + (1 : ℝ) * D00.Li2 (1 / Real.sqrt 2) * lam + (-2 : ℝ) * Li3 (-(1 / Real.sqrt 2)) + (2 : ℝ) * Li3 (1 / Real.sqrt 2) := by
  have e0 : mb b (-b) a = 1 / Real.sqrt 2 := by unfold mb; have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf)
  have e1 : s1 b (-b) 0 a = -(1 / Real.sqrt 2) := by unfold s1; have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf)
  have l0 : Real.log ((0 - b) / (0 - (-b))) = (0 : ℝ) * lam + (0 : ℝ) * L :=
    logN 0 0 0 0 _ _ (by norm_num) (by norm_num)
      (by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf))
  have l1 : Real.log (-(1 / Real.sqrt 2)) = (-1 / 2 : ℝ) * lam + (0 : ℝ) * L :=
    logN 0 0 1 0 _ _ (by norm_num) (by norm_num)
      (by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf))
  have l2 : Real.log (1 - (-(1 / Real.sqrt 2))) = (-1 / 2 : ℝ) * lam + (1 : ℝ) * L :=
    logP 0 1 1 0 _ _ (by norm_num) (by norm_num)
      (by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf))
  have l3 : Real.log (1 / Real.sqrt 2) = (-1 / 2 : ℝ) * lam + (0 : ℝ) * L :=
    logP 0 0 1 0 _ _ (by norm_num) (by norm_num)
      (by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf))
  have l4 : Real.log (1 - (1 / Real.sqrt 2)) = (-1 / 2 : ℝ) * lam + (-1 : ℝ) * L :=
    logP 0 0 1 1 _ _ (by norm_num) (by norm_num)
      (by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf))
  simp only [FG, e0, e1, l0, Phi, P1n, P2a, l1, l2, l3, l4]
  try simp only [Real.log_zero, sub_zero, Real.log_one, D00.Li2_zero', Li3_zero, D00.Li2_one']
  ring

lemma val35p : FG P1n P2a P2a b (-b) 0 b = 0 := by
  have e0 : mb b (-b) b = 0 := by unfold mb; have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf)
  have e1 : s1 b (-b) 0 b = 0 := by unfold s1; have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf)
  have l0 : Real.log ((0 - b) / (0 - (-b))) = (0 : ℝ) * lam + (0 : ℝ) * L :=
    logN 0 0 0 0 _ _ (by norm_num) (by norm_num)
      (by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf))
  simp only [FG, e0, e1, l0, Phi, P1n, P2a]
  try simp only [Real.log_zero, sub_zero, Real.log_one, D00.Li2_zero', Li3_zero, D00.Li2_one']
  ring

lemma blk36 : IntervalIntegrable (fpr (-a) a 0) volume b 1 ∧
    ∫ v in b..1, fpr (-a) a 0 v = FG P1p P2g P2a (-a) a 0 1 - FG P1p P2g P2a (-a) a 0 b := by
  have hxy : ((-a) : ℝ) - a ≠ 0 := ne_of_lt (by (try simp only [Y01.a, Y01.b]); linarith [s2lo, s2hi])
  have hzx : (0 : ℝ) - (-a) ≠ 0 := ne_of_gt (by (try simp only [Y01.a, Y01.b]); linarith [s2lo, s2hi])
  have hzy : (0 : ℝ) - a ≠ 0 := ne_of_lt (by (try simp only [Y01.a, Y01.b]); linarith [s2lo, s2hi])
  have hz : ∀ v ∈ Icc b 1, 0 < v - 0 := fun v hv => by
    have h1 := hv.1; have h2 := hv.2; simp only [Y01.a, Y01.b] at *; linarith [s2lo, s2hi]
  have hy : ∀ v ∈ Icc b 1, v - a ≠ 0 := fun v hv => ne_of_lt (by
    have h1 := hv.1; have h2 := hv.2; simp only [Y01.a, Y01.b] at *; linarith [s2lo, s2hi])
  have hx : ∀ v ∈ Ioo b 1, v - (-a) ≠ 0 := fun v hv => ne_of_gt (by
    have h1 := hv.1; have h2 := hv.2; simp only [Y01.a, Y01.b] at *; linarith [s2lo, s2hi])
  have hx2 : ∀ v ∈ Icc b 1, v - (-a) ≠ 0 := fun v hv => ne_of_gt (by
    have h1 := hv.1; have h2 := hv.2; simp only [Y01.a, Y01.b] at *; linarith [s2lo, s2hi])
  exact intG isP1p isP2g isP2a hb1 hxy hzx hzy hz hy hx
    (fun v hv => cPhiG (s1_ne_zero hzx hzy (hx2 v hv) (hy v hv)) (s1_ne_one hxy hzx (hy v hv) (hz v hv).ne'))
    (fun v hv => contAt_P2a (mb_ne_one hxy (hy v hv)))

lemma val36q : FG P1p P2g P2a (-a) a 0 1 = (-5 / 3 : ℝ) * L ^ 3 + (-2 : ℝ) * L * D00.Li2 (b) + (-2 : ℝ) * L * D00.Li2 (-a) + (-2 : ℝ) * Li3 (b) + (2 : ℝ) * Li3 (-a) := by
  have e0 : mb (-a) a 1 = -a := by unfold mb; have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf)
  have e1 : s1 (-a) a 0 1 = (Real.sqrt 2 ^ 0 * a ^ 1 / (Real.sqrt 2 ^ 0 * a ^ 0)) := by unfold s1; have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf)
  have l0 : Real.log ((0 - (-a)) / (0 - a)) = (0 : ℝ) * lam + (0 : ℝ) * L :=
    logN 0 0 0 0 _ _ (by norm_num) (by norm_num)
      (by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf))
  have e2 : (1:ℝ) - ((Real.sqrt 2 ^ 0 * a ^ 1 / (Real.sqrt 2 ^ 0 * a ^ 0))) = -Real.sqrt 2 := by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf)
  have e3 : ((Real.sqrt 2 ^ 0 * a ^ 1 / (Real.sqrt 2 ^ 0 * a ^ 0)))⁻¹ = b := by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf)
  have l1 : Real.log (-Real.sqrt 2) = (1 / 2 : ℝ) * lam + (0 : ℝ) * L :=
    logN 1 0 0 0 _ _ (by norm_num) (by norm_num)
      (by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf))
  have l2 : Real.log (b) = (0 : ℝ) * lam + (-1 : ℝ) * L :=
    logP 0 0 0 1 _ _ (by norm_num) (by norm_num)
      (by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf))
  have l3 : Real.log (1 - (b)) = (1 / 2 : ℝ) * lam + (-1 : ℝ) * L :=
    logP 1 0 0 1 _ _ (by norm_num) (by norm_num)
      (by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf))
  have l4 : Real.log (-a) = (0 : ℝ) * lam + (1 : ℝ) * L :=
    logN 0 1 0 0 _ _ (by norm_num) (by norm_num)
      (by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf))
  have l5 : Real.log (1 - (-a)) = (1 / 2 : ℝ) * lam + (1 : ℝ) * L :=
    logP 1 1 0 0 _ _ (by norm_num) (by norm_num)
      (by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf))
  simp only [FG, e0, e1, l0, Phi, P1p, P2g, e2, e3, E3, l1, l2, l3, P2a, l4, l5]
  try simp only [Real.log_zero, sub_zero, Real.log_one, D00.Li2_zero', Li3_zero, D00.Li2_one']
  ring

lemma val36p : FG P1p P2g P2a (-a) a 0 b = (-1 / 2 : ℝ) * L * lam ^ 2 + (-1 : ℝ) * D00.Li2 (-Real.sqrt 2) * lam + (-1 : ℝ) * D00.Li2 (1 / Real.sqrt 2) * lam + (2 : ℝ) * Li3 (-Real.sqrt 2) + (-2 : ℝ) * Li3 (1 / Real.sqrt 2) + (-1 / 12 : ℝ) * lam ^ 3 := by
  have e0 : mb (-a) a b = -Real.sqrt 2 := by unfold mb; have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf)
  have e1 : s1 (-a) a 0 b = (Real.sqrt 2 ^ 1 * a ^ 0 / (Real.sqrt 2 ^ 0 * a ^ 0)) := by unfold s1; have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf)
  have l0 : Real.log ((0 - (-a)) / (0 - a)) = (0 : ℝ) * lam + (0 : ℝ) * L :=
    logN 0 0 0 0 _ _ (by norm_num) (by norm_num)
      (by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf))
  have e2 : (1:ℝ) - ((Real.sqrt 2 ^ 1 * a ^ 0 / (Real.sqrt 2 ^ 0 * a ^ 0))) = -b := by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf)
  have e3 : ((Real.sqrt 2 ^ 1 * a ^ 0 / (Real.sqrt 2 ^ 0 * a ^ 0)))⁻¹ = 1 / Real.sqrt 2 := by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf)
  have l1 : Real.log (-b) = (0 : ℝ) * lam + (-1 : ℝ) * L :=
    logN 0 0 0 1 _ _ (by norm_num) (by norm_num)
      (by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf))
  have l2 : Real.log (1 / Real.sqrt 2) = (-1 / 2 : ℝ) * lam + (0 : ℝ) * L :=
    logP 0 0 1 0 _ _ (by norm_num) (by norm_num)
      (by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf))
  have l3 : Real.log (1 - (1 / Real.sqrt 2)) = (-1 / 2 : ℝ) * lam + (-1 : ℝ) * L :=
    logP 0 0 1 1 _ _ (by norm_num) (by norm_num)
      (by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf))
  have l4 : Real.log (-Real.sqrt 2) = (1 / 2 : ℝ) * lam + (0 : ℝ) * L :=
    logN 1 0 0 0 _ _ (by norm_num) (by norm_num)
      (by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf))
  have l5 : Real.log (1 - (-Real.sqrt 2)) = (0 : ℝ) * lam + (1 : ℝ) * L :=
    logP 0 1 0 0 _ _ (by norm_num) (by norm_num)
      (by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf))
  simp only [FG, e0, e1, l0, Phi, P1p, P2g, e2, e3, E3, l1, l2, l3, P2a, l4, l5]
  try simp only [Real.log_zero, sub_zero, Real.log_one, D00.Li2_zero', Li3_zero, D00.Li2_one']
  ring

lemma blk37 : IntervalIntegrable (fpr (-a) (-b) 0) volume b a ∧
    ∫ v in b..a, fpr (-a) (-b) 0 v = FG P1p P2a P2g (-a) (-b) 0 a - FG P1p P2a P2g (-a) (-b) 0 b := by
  have hxy : ((-a) : ℝ) - (-b) ≠ 0 := ne_of_lt (by (try simp only [Y01.a, Y01.b]); linarith [s2lo, s2hi])
  have hzx : (0 : ℝ) - (-a) ≠ 0 := ne_of_gt (by (try simp only [Y01.a, Y01.b]); linarith [s2lo, s2hi])
  have hzy : (0 : ℝ) - (-b) ≠ 0 := ne_of_gt (by (try simp only [Y01.a, Y01.b]); linarith [s2lo, s2hi])
  have hz : ∀ v ∈ Icc b a, 0 < v - 0 := fun v hv => by
    have h1 := hv.1; have h2 := hv.2; simp only [Y01.a, Y01.b] at *; linarith [s2lo, s2hi]
  have hy : ∀ v ∈ Icc b a, v - (-b) ≠ 0 := fun v hv => ne_of_gt (by
    have h1 := hv.1; have h2 := hv.2; simp only [Y01.a, Y01.b] at *; linarith [s2lo, s2hi])
  have hx : ∀ v ∈ Ioo b a, v - (-a) ≠ 0 := fun v hv => ne_of_gt (by
    have h1 := hv.1; have h2 := hv.2; simp only [Y01.a, Y01.b] at *; linarith [s2lo, s2hi])
  have hx2 : ∀ v ∈ Icc b a, v - (-a) ≠ 0 := fun v hv => ne_of_gt (by
    have h1 := hv.1; have h2 := hv.2; simp only [Y01.a, Y01.b] at *; linarith [s2lo, s2hi])
  exact intG isP1p isP2a isP2g hba hxy hzx hzy hz hy hx
    (fun v hv => cPhiM (s1_ne_one hxy hzx (hy v hv) (hz v hv).ne'))
    (fun v hv => contAt_P2g (mb_ne_zero (hx2 v hv) (hy v hv)) (mb_ne_one hxy (hy v hv)))

lemma val37q : FG P1p P2a P2g (-a) (-b) 0 a = (2 / 3 : ℝ) * L ^ 3 + (-3 : ℝ) * L ^ 2 * lam + (-2 : ℝ) * L * D00.Li2 (1 - 1 / Real.sqrt 2) + (2 : ℝ) * L * D00.Li2 (2 - Real.sqrt 2) + (-4 : ℝ) * L * D00.Li2 (1 / Real.sqrt 2) + (-1 / 2 : ℝ) * L * lam ^ 2 + (-1 : ℝ) * D00.Li2 (1 - 1 / Real.sqrt 2) * lam + (-1 : ℝ) * D00.Li2 (2 - Real.sqrt 2) * lam + (-2 : ℝ) * Li3 (1 - 1 / Real.sqrt 2) + (2 : ℝ) * Li3 (2 - Real.sqrt 2) + (-1 / 12 : ℝ) * lam ^ 3 := by
  have e0 : mb (-a) (-b) a = (Real.sqrt 2 ^ 0 * a ^ 1 / (Real.sqrt 2 ^ 1 * a ^ 0)) := by unfold mb; have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf)
  have e1 : s1 (-a) (-b) 0 a = 1 - 1 / Real.sqrt 2 := by unfold s1; have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf)
  have l0 : Real.log ((0 - (-a)) / (0 - (-b))) = (0 : ℝ) * lam + (2 : ℝ) * L :=
    logP 0 2 0 0 _ _ (by norm_num) (by norm_num)
      (by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf))
  have e2 : (1:ℝ) - (1 - 1 / Real.sqrt 2) = 1 / Real.sqrt 2 := by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf)
  have l1 : Real.log (1 - 1 / Real.sqrt 2) = (-1 / 2 : ℝ) * lam + (-1 : ℝ) * L :=
    logP 0 0 1 1 _ _ (by norm_num) (by norm_num)
      (by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf))
  have l2 : Real.log (1 / Real.sqrt 2) = (-1 / 2 : ℝ) * lam + (0 : ℝ) * L :=
    logP 0 0 1 0 _ _ (by norm_num) (by norm_num)
      (by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf))
  have e3 : ((Real.sqrt 2 ^ 0 * a ^ 1 / (Real.sqrt 2 ^ 1 * a ^ 0)))⁻¹ = 2 - Real.sqrt 2 := by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf)
  have l3 : Real.log (2 - Real.sqrt 2) = (1 / 2 : ℝ) * lam + (-1 : ℝ) * L :=
    logP 1 0 0 1 _ _ (by norm_num) (by norm_num)
      (by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf))
  have l4 : Real.log (1 - (2 - Real.sqrt 2)) = (0 : ℝ) * lam + (-1 : ℝ) * L :=
    logP 0 0 0 1 _ _ (by norm_num) (by norm_num)
      (by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf))
  simp only [FG, e0, e1, l0, Phi, P1p, P2a, e2, l1, l2, P2g, E3, e3, l3, l4]
  try simp only [Real.log_zero, sub_zero, Real.log_one, D00.Li2_zero', Li3_zero, D00.Li2_one']
  ring

lemma val37p : FG P1p P2a P2g (-a) (-b) 0 b = (-16 / 3 : ℝ) * L ^ 3 + (1 : ℝ) * L ^ 2 * lam + (-4 : ℝ) * L * D00.Li2 (b) + (2 : ℝ) * L * D00.Li2 (1 - 1 / Real.sqrt 2) + (-2 : ℝ) * L * D00.Li2 (2 - Real.sqrt 2) + (1 : ℝ) * D00.Li2 (1 - 1 / Real.sqrt 2) * lam + (1 : ℝ) * D00.Li2 (2 - Real.sqrt 2) * lam + (2 : ℝ) * Li3 (1 - 1 / Real.sqrt 2) + (-2 : ℝ) * Li3 (2 - Real.sqrt 2) + (1 / 12 : ℝ) * lam ^ 3 := by
  have e0 : mb (-a) (-b) b = (Real.sqrt 2 ^ 1 * a ^ 1 / (Real.sqrt 2 ^ 0 * a ^ 0)) := by unfold mb; have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf)
  have e1 : s1 (-a) (-b) 0 b = 2 - Real.sqrt 2 := by unfold s1; have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf)
  have l0 : Real.log ((0 - (-a)) / (0 - (-b))) = (0 : ℝ) * lam + (2 : ℝ) * L :=
    logP 0 2 0 0 _ _ (by norm_num) (by norm_num)
      (by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf))
  have e2 : (1:ℝ) - (2 - Real.sqrt 2) = b := by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf)
  have l1 : Real.log (2 - Real.sqrt 2) = (1 / 2 : ℝ) * lam + (-1 : ℝ) * L :=
    logP 1 0 0 1 _ _ (by norm_num) (by norm_num)
      (by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf))
  have l2 : Real.log (b) = (0 : ℝ) * lam + (-1 : ℝ) * L :=
    logP 0 0 0 1 _ _ (by norm_num) (by norm_num)
      (by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf))
  have e3 : ((Real.sqrt 2 ^ 1 * a ^ 1 / (Real.sqrt 2 ^ 0 * a ^ 0)))⁻¹ = 1 - 1 / Real.sqrt 2 := by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf)
  have l3 : Real.log (1 - 1 / Real.sqrt 2) = (-1 / 2 : ℝ) * lam + (-1 : ℝ) * L :=
    logP 0 0 1 1 _ _ (by norm_num) (by norm_num)
      (by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf))
  have l4 : Real.log (1 - (1 - 1 / Real.sqrt 2)) = (-1 / 2 : ℝ) * lam + (0 : ℝ) * L :=
    logP 0 0 1 0 _ _ (by norm_num) (by norm_num)
      (by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf))
  simp only [FG, e0, e1, l0, Phi, P1p, P2a, e2, l1, l2, P2g, E3, e3, l3, l4]
  try simp only [Real.log_zero, sub_zero, Real.log_one, D00.Li2_zero', Li3_zero, D00.Li2_one']
  ring

lemma blk38 : IntervalIntegrable (fpr 0 a 0) volume b 1 ∧
    ∫ v in b..1, fpr 0 a 0 v = FIII P2a 0 a 1 - FIII P2a 0 a b := by
  have hxy : (0 : ℝ) - a ≠ 0 := ne_of_lt (by (try simp only [Y01.a, Y01.b]); linarith [s2lo, s2hi])
  have hx : ∀ v ∈ Icc b 1, 0 < v - 0 := fun v hv => by
    have h1 := hv.1; have h2 := hv.2; simp only [Y01.a, Y01.b] at *; linarith [s2lo, s2hi]
  have hy : ∀ v ∈ Icc b 1, v - a ≠ 0 := fun v hv => ne_of_lt (by
    have h1 := hv.1; have h2 := hv.2; simp only [Y01.a, Y01.b] at *; linarith [s2lo, s2hi])
  exact intIII isP2a hb1 hxy hx hy (fun v hv => contAt_P2a (mb_ne_one hxy (hy v hv)))

lemma val38q : FIII P2a 0 a 1 = (-1 / 4 : ℝ) * L * lam ^ 2 + (1 : ℝ) * D00.Li2 (-(1 / Real.sqrt 2)) * lam + (2 : ℝ) * Li3 (-(1 / Real.sqrt 2)) + (1 / 12 : ℝ) * lam ^ 3 := by
  have e0 : mb 0 a 1 = -(1 / Real.sqrt 2) := by unfold mb; have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf)
  have l0 : Real.log (1 - 0) = (0 : ℝ) * lam + (0 : ℝ) * L :=
    logP 0 0 0 0 _ _ (by norm_num) (by norm_num)
      (by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf))
  have l1 : Real.log (1 - a) = (1 / 2 : ℝ) * lam + (0 : ℝ) * L :=
    logN 1 0 0 0 _ _ (by norm_num) (by norm_num)
      (by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf))
  have l2 : Real.log (-(1 / Real.sqrt 2)) = (-1 / 2 : ℝ) * lam + (0 : ℝ) * L :=
    logN 0 0 1 0 _ _ (by norm_num) (by norm_num)
      (by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf))
  have l3 : Real.log (1 - (-(1 / Real.sqrt 2))) = (-1 / 2 : ℝ) * lam + (1 : ℝ) * L :=
    logP 0 1 1 0 _ _ (by norm_num) (by norm_num)
      (by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf))
  simp only [FIII, e0, l0, l1, P2a, l2, l3]
  try simp only [Real.log_zero, sub_zero, Real.log_one, D00.Li2_zero', Li3_zero, D00.Li2_one']
  ring

lemma val38p : FIII P2a 0 a b = (-4 / 3 : ℝ) * L ^ 3 + (-2 : ℝ) * L ^ 2 * lam + (2 : ℝ) * L * D00.Li2 (-(b / 2)) + (2 : ℝ) * D00.Li2 (-(b / 2)) * lam + (2 : ℝ) * Li3 (-(b / 2)) + (2 / 3 : ℝ) * lam ^ 3 := by
  have e0 : mb 0 a b = -(b / 2) := by unfold mb; have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf)
  have l0 : Real.log (b - 0) = (0 : ℝ) * lam + (-1 : ℝ) * L :=
    logP 0 0 0 1 _ _ (by norm_num) (by norm_num)
      (by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf))
  have l1 : Real.log (b - a) = (1 : ℝ) * lam + (0 : ℝ) * L :=
    logN 2 0 0 0 _ _ (by norm_num) (by norm_num)
      (by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf))
  have l2 : Real.log (-(b / 2)) = (-1 : ℝ) * lam + (-1 : ℝ) * L :=
    logN 0 0 2 1 _ _ (by norm_num) (by norm_num)
      (by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf))
  have l3 : Real.log (1 - (-(b / 2))) = (-1 : ℝ) * lam + (1 : ℝ) * L :=
    logP 0 1 2 0 _ _ (by norm_num) (by norm_num)
      (by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf))
  simp only [FIII, e0, l0, l1, P2a, l2, l3]
  try simp only [Real.log_zero, sub_zero, Real.log_one, D00.Li2_zero', Li3_zero, D00.Li2_one']
  ring

lemma blk39 : IntervalIntegrable (fpr 0 (-b) 0) volume b a ∧
    ∫ v in b..a, fpr 0 (-b) 0 v = FIII P2a 0 (-b) a - FIII P2a 0 (-b) b := by
  have hxy : (0 : ℝ) - (-b) ≠ 0 := ne_of_gt (by (try simp only [Y01.a, Y01.b]); linarith [s2lo, s2hi])
  have hx : ∀ v ∈ Icc b a, 0 < v - 0 := fun v hv => by
    have h1 := hv.1; have h2 := hv.2; simp only [Y01.a, Y01.b] at *; linarith [s2lo, s2hi]
  have hy : ∀ v ∈ Icc b a, v - (-b) ≠ 0 := fun v hv => ne_of_gt (by
    have h1 := hv.1; have h2 := hv.2; simp only [Y01.a, Y01.b] at *; linarith [s2lo, s2hi])
  exact intIII isP2a hba hxy hx hy (fun v hv => contAt_P2a (mb_ne_one hxy (hy v hv)))

lemma val39q : FIII P2a 0 (-b) a = (4 / 3 : ℝ) * L ^ 3 + (-3 : ℝ) * L ^ 2 * lam + (-2 : ℝ) * L * D00.Li2 ((2 + Real.sqrt 2) / 4) + (3 : ℝ) * D00.Li2 ((2 + Real.sqrt 2) / 4) * lam + (2 : ℝ) * Li3 ((2 + Real.sqrt 2) / 4) + (9 / 4 : ℝ) * lam ^ 3 := by
  have e0 : mb 0 (-b) a = (2 + Real.sqrt 2) / 4 := by unfold mb; have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf)
  have l0 : Real.log (a - 0) = (0 : ℝ) * lam + (1 : ℝ) * L :=
    logP 0 1 0 0 _ _ (by norm_num) (by norm_num)
      (by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf))
  have l1 : Real.log (a - (-b)) = (3 / 2 : ℝ) * lam + (0 : ℝ) * L :=
    logP 3 0 0 0 _ _ (by norm_num) (by norm_num)
      (by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf))
  have l2 : Real.log ((2 + Real.sqrt 2) / 4) = (-3 / 2 : ℝ) * lam + (1 : ℝ) * L :=
    logP 0 1 3 0 _ _ (by norm_num) (by norm_num)
      (by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf))
  have l3 : Real.log (1 - ((2 + Real.sqrt 2) / 4)) = (-3 / 2 : ℝ) * lam + (-1 : ℝ) * L :=
    logP 0 0 3 1 _ _ (by norm_num) (by norm_num)
      (by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf))
  simp only [FIII, e0, l0, l1, P2a, l2, l3]
  try simp only [Real.log_zero, sub_zero, Real.log_one, D00.Li2_zero', Li3_zero, D00.Li2_one']
  ring

lemma val39p : FIII P2a 0 (-b) b = (2 : ℝ) * D00.Li2 (1 / 2) * lam + (2 : ℝ) * Li3 (1 / 2) + (2 / 3 : ℝ) * lam ^ 3 := by
  have e0 : mb 0 (-b) b = 1 / 2 := by unfold mb; have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf)
  have l0 : Real.log (b - 0) = (0 : ℝ) * lam + (-1 : ℝ) * L :=
    logP 0 0 0 1 _ _ (by norm_num) (by norm_num)
      (by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf))
  have l1 : Real.log (b - (-b)) = (1 : ℝ) * lam + (-1 : ℝ) * L :=
    logP 2 0 0 1 _ _ (by norm_num) (by norm_num)
      (by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf))
  have l2 : Real.log (1 / 2) = (-1 : ℝ) * lam + (0 : ℝ) * L :=
    logP 0 0 2 0 _ _ (by norm_num) (by norm_num)
      (by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf))
  have l3 : Real.log (1 - (1 / 2)) = (-1 : ℝ) * lam + (0 : ℝ) * L :=
    logP 0 0 2 0 _ _ (by norm_num) (by norm_num)
      (by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf))
  simp only [FIII, e0, l0, l1, P2a, l2, l3]
  try simp only [Real.log_zero, sub_zero, Real.log_one, D00.Li2_zero', Li3_zero, D00.Li2_one']
  ring

lemma blk40 : IntervalIntegrable (fpr a b 0) volume 1 a ∧
    ∫ v in 1..a, fpr a b 0 v = FG P1n P2a P2a a b 0 a - FG P1n P2a P2a a b 0 1 := by
  have hxy : (a : ℝ) - b ≠ 0 := ne_of_gt (by (try simp only [Y01.a, Y01.b]); linarith [s2lo, s2hi])
  have hzx : (0 : ℝ) - a ≠ 0 := ne_of_lt (by (try simp only [Y01.a, Y01.b]); linarith [s2lo, s2hi])
  have hzy : (0 : ℝ) - b ≠ 0 := ne_of_lt (by (try simp only [Y01.a, Y01.b]); linarith [s2lo, s2hi])
  have hz : ∀ v ∈ Icc 1 a, 0 < v - 0 := fun v hv => by
    have h1 := hv.1; have h2 := hv.2; simp only [Y01.a, Y01.b] at *; linarith [s2lo, s2hi]
  have hy : ∀ v ∈ Icc 1 a, v - b ≠ 0 := fun v hv => ne_of_gt (by
    have h1 := hv.1; have h2 := hv.2; simp only [Y01.a, Y01.b] at *; linarith [s2lo, s2hi])
  have hx : ∀ v ∈ Ioo 1 a, v - a ≠ 0 := fun v hv => ne_of_lt (by
    have h1 := hv.1; have h2 := hv.2; simp only [Y01.a, Y01.b] at *; linarith [s2lo, s2hi])
  exact intG isP1n isP2a isP2a h1a hxy hzx hzy hz hy hx
    (fun v hv => cPhiN (s1_ne_one hxy hzx (hy v hv) (hz v hv).ne'))
    (fun v hv => contAt_P2a (mb_ne_one hxy (hy v hv)))

lemma val40q : FG P1n P2a P2a a b 0 a = 0 := by
  have e0 : mb a b a = 0 := by unfold mb; have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf)
  have e1 : s1 a b 0 a = 0 := by unfold s1; have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf)
  have l0 : Real.log ((0 - a) / (0 - b)) = (0 : ℝ) * lam + (2 : ℝ) * L :=
    logP 0 2 0 0 _ _ (by norm_num) (by norm_num)
      (by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf))
  simp only [FG, e0, e1, l0, Phi, P1n, P2a]
  try simp only [Real.log_zero, sub_zero, Real.log_one, D00.Li2_zero', Li3_zero, D00.Li2_one']
  ring

lemma val40p : FG P1n P2a P2a a b 0 1 = (-1 : ℝ) * L ^ 3 + (-2 : ℝ) * L * D00.Li2 (-a) + (2 : ℝ) * L * D00.Li2 (-b) + (2 : ℝ) * Li3 (-a) + (-2 : ℝ) * Li3 (-b) := by
  have e0 : mb a b 1 = -a := by unfold mb; have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf)
  have e1 : s1 a b 0 1 = -b := by unfold s1; have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf)
  have l0 : Real.log ((0 - a) / (0 - b)) = (0 : ℝ) * lam + (2 : ℝ) * L :=
    logP 0 2 0 0 _ _ (by norm_num) (by norm_num)
      (by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf))
  have l1 : Real.log (-b) = (0 : ℝ) * lam + (-1 : ℝ) * L :=
    logN 0 0 0 1 _ _ (by norm_num) (by norm_num)
      (by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf))
  have l2 : Real.log (1 - (-b)) = (1 / 2 : ℝ) * lam + (0 : ℝ) * L :=
    logP 1 0 0 0 _ _ (by norm_num) (by norm_num)
      (by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf))
  have l3 : Real.log (-a) = (0 : ℝ) * lam + (1 : ℝ) * L :=
    logN 0 1 0 0 _ _ (by norm_num) (by norm_num)
      (by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf))
  have l4 : Real.log (1 - (-a)) = (1 / 2 : ℝ) * lam + (1 : ℝ) * L :=
    logP 1 1 0 0 _ _ (by norm_num) (by norm_num)
      (by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf))
  simp only [FG, e0, e1, l0, Phi, P1n, P2a, l1, l2, l3, l4]
  try simp only [Real.log_zero, sub_zero, Real.log_one, D00.Li2_zero', Li3_zero, D00.Li2_one']
  ring

lemma blk41 : IntervalIntegrable (fpr a (-a) 0) volume 1 a ∧
    ∫ v in 1..a, fpr a (-a) 0 v = FG P1p P2a P2a a (-a) 0 a - FG P1p P2a P2a a (-a) 0 1 := by
  have hxy : (a : ℝ) - (-a) ≠ 0 := ne_of_gt (by (try simp only [Y01.a, Y01.b]); linarith [s2lo, s2hi])
  have hzx : (0 : ℝ) - a ≠ 0 := ne_of_lt (by (try simp only [Y01.a, Y01.b]); linarith [s2lo, s2hi])
  have hzy : (0 : ℝ) - (-a) ≠ 0 := ne_of_gt (by (try simp only [Y01.a, Y01.b]); linarith [s2lo, s2hi])
  have hz : ∀ v ∈ Icc 1 a, 0 < v - 0 := fun v hv => by
    have h1 := hv.1; have h2 := hv.2; simp only [Y01.a, Y01.b] at *; linarith [s2lo, s2hi]
  have hy : ∀ v ∈ Icc 1 a, v - (-a) ≠ 0 := fun v hv => ne_of_gt (by
    have h1 := hv.1; have h2 := hv.2; simp only [Y01.a, Y01.b] at *; linarith [s2lo, s2hi])
  have hx : ∀ v ∈ Ioo 1 a, v - a ≠ 0 := fun v hv => ne_of_lt (by
    have h1 := hv.1; have h2 := hv.2; simp only [Y01.a, Y01.b] at *; linarith [s2lo, s2hi])
  exact intG isP1p isP2a isP2a h1a hxy hzx hzy hz hy hx
    (fun v hv => cPhiM (s1_ne_one hxy hzx (hy v hv) (hz v hv).ne'))
    (fun v hv => contAt_P2a (mb_ne_one hxy (hy v hv)))

lemma val41q : FG P1p P2a P2a a (-a) 0 a = 0 := by
  have e0 : mb a (-a) a = 0 := by unfold mb; have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf)
  have e1 : s1 a (-a) 0 a = 0 := by unfold s1; have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf)
  have l0 : Real.log ((0 - a) / (0 - (-a))) = (0 : ℝ) * lam + (0 : ℝ) * L :=
    logN 0 0 0 0 _ _ (by norm_num) (by norm_num)
      (by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf))
  simp only [FG, e0, e1, l0, Phi, P1p, P2a]
  try simp only [Real.log_zero, sub_zero, Real.log_one, D00.Li2_zero', Li3_zero, D00.Li2_one']
  ring

lemma val41p : FG P1p P2a P2a a (-a) 0 1 = (-1 : ℝ) * L ^ 3 + (-2 : ℝ) * L * D00.Li2 (b) + (2 : ℝ) * L * D00.Li2 (-b) + (-2 : ℝ) * Li3 (b) + (2 : ℝ) * Li3 (-b) := by
  have e0 : mb a (-a) 1 = -b := by unfold mb; have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf)
  have e1 : s1 a (-a) 0 1 = b := by unfold s1; have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf)
  have l0 : Real.log ((0 - a) / (0 - (-a))) = (0 : ℝ) * lam + (0 : ℝ) * L :=
    logN 0 0 0 0 _ _ (by norm_num) (by norm_num)
      (by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf))
  have e2 : (1:ℝ) - (b) = 2 - Real.sqrt 2 := by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf)
  have l1 : Real.log (b) = (0 : ℝ) * lam + (-1 : ℝ) * L :=
    logP 0 0 0 1 _ _ (by norm_num) (by norm_num)
      (by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf))
  have l2 : Real.log (2 - Real.sqrt 2) = (1 / 2 : ℝ) * lam + (-1 : ℝ) * L :=
    logP 1 0 0 1 _ _ (by norm_num) (by norm_num)
      (by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf))
  have l3 : Real.log (-b) = (0 : ℝ) * lam + (-1 : ℝ) * L :=
    logN 0 0 0 1 _ _ (by norm_num) (by norm_num)
      (by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf))
  have l4 : Real.log (1 - (-b)) = (1 / 2 : ℝ) * lam + (0 : ℝ) * L :=
    logP 1 0 0 0 _ _ (by norm_num) (by norm_num)
      (by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf))
  simp only [FG, e0, e1, l0, Phi, P1p, P2a, e2, l1, l2, l3, l4]
  try simp only [Real.log_zero, sub_zero, Real.log_one, D00.Li2_zero', Li3_zero, D00.Li2_one']
  ring

lemma blk42 : IntervalIntegrable (fpr a 0 0) volume 1 a ∧
    ∫ v in 1..a, fpr a 0 0 v = FIV P2a a 0 a - FIV P2a a 0 1 := by
  have hxy : (a : ℝ) - 0 ≠ 0 := ne_of_gt (by (try simp only [Y01.a, Y01.b]); linarith [s2lo, s2hi])
  have hy : ∀ v ∈ Icc 1 a, 0 < v - 0 := fun v hv => by
    have h1 := hv.1; have h2 := hv.2; simp only [Y01.a, Y01.b] at *; linarith [s2lo, s2hi]
  have hx : ∀ v ∈ Ioo 1 a, v - a ≠ 0 := fun v hv => ne_of_lt (by
    have h1 := hv.1; have h2 := hv.2; simp only [Y01.a, Y01.b] at *; linarith [s2lo, s2hi])
  exact intIV isP2a h1a hxy hy hx (fun v hv => contAt_P2a (mb_ne_one hxy (hy v hv).ne'))

lemma val42q : FIV P2a a 0 a = 0 := by
  have e0 : mb a 0 a = 0 := by unfold mb; have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf)
  simp only [FIV, e0, P2a]
  try simp only [Real.log_zero, sub_zero, Real.log_one, D00.Li2_zero', Li3_zero, D00.Li2_one']
  ring

lemma val42p : FIV P2a a 0 1 = (-1 / 4 : ℝ) * L * lam ^ 2 + (-1 : ℝ) * D00.Li2 (-Real.sqrt 2) * lam + (2 : ℝ) * Li3 (-Real.sqrt 2) := by
  have e0 : mb a 0 1 = -Real.sqrt 2 := by unfold mb; have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf)
  have l0 : Real.log (-Real.sqrt 2) = (1 / 2 : ℝ) * lam + (0 : ℝ) * L :=
    logN 1 0 0 0 _ _ (by norm_num) (by norm_num)
      (by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf))
  have l1 : Real.log (1 - (-Real.sqrt 2)) = (0 : ℝ) * lam + (1 : ℝ) * L :=
    logP 0 1 0 0 _ _ (by norm_num) (by norm_num)
      (by have := nza; have := nzb; have := nzs; have := nzd0; have := nzd1; have := nzd2; have := nzd3; have := nzd4; have := nzd5; have := nzd6; have := nzd7; have := nzd8; have := nzd9; have := nzd10; have := nzd11; have := nzd12; have := nzd13; have := nzd14; have := nzd15; have := nzd16; have := nzd17; have := nzd18; have := nzd19; have := nzd20; have := nzd21; have := nzd22; have := nzd23; have := nzd24; have := nzd25; have := nzd26; have := nzd27; have := nzd28; have := nzd29; (try field_simp); (try simp only [Y01.a, Y01.b]); (try ring_nf); (try simp only [s2_2, s2_3, s2_4, s2_5, s2_6, s2_7, s2_8, s2_9, s2_10, s2_11, s2_12, s2_13, s2_14, s2_15, s2_16, s2_17, s2_18, s2_19, s2_20, s2_21, s2_22, s2_23, s2_24]); (try ring_nf))
  simp only [FIV, e0, P2a, l0, l1]
  try simp only [Real.log_zero, sub_zero, Real.log_one, D00.Li2_zero', Li3_zero, D00.Li2_one']
  ring


/-! ## Assembly: linearity, the pointwise decomposition, and the sum -/

lemma lin_one {f : ℝ → ℝ} {p q A : ℝ} (c : ℝ)
    (h : IntervalIntegrable f volume p q ∧ ∫ v in p..q, f v = A) :
    IntervalIntegrable (fun v => c * f v) volume p q ∧ ∫ v in p..q, c * f v = c * A :=
  ⟨h.1.const_mul c, by rw [intervalIntegral.integral_const_mul, h.2]⟩

lemma lin_cons {f g : ℝ → ℝ} {p q A B : ℝ} (c : ℝ)
    (hf : IntervalIntegrable f volume p q ∧ ∫ v in p..q, f v = A)
    (hg : IntervalIntegrable g volume p q ∧ ∫ v in p..q, g v = B) :
    IntervalIntegrable (fun v => c * f v + g v) volume p q ∧
      ∫ v in p..q, (c * f v + g v) = c * A + B :=
  ⟨(hf.1.const_mul c).add hg.1, by
    rw [intervalIntegral.integral_add (hf.1.const_mul c) hg.1, intervalIntegral.integral_const_mul,
      hf.2, hg.2]⟩

lemma congr_Ioo {f g : ℝ → ℝ} {p q A : ℝ} (hpq : p ≤ q) (h : ∀ v ∈ Ioo p q, f v = g v)
    (hg : IntervalIntegrable g volume p q ∧ ∫ v in p..q, g v = A) :
    IntervalIntegrable f volume p q ∧ ∫ v in p..q, f v = A := by
  have hae : ∀ᵐ v ∂volume, v ∈ Set.uIoc p q → f v = g v := by
    filter_upwards [(Set.countable_singleton q).ae_notMem volume] with v hv hvI
    rw [uIoc_of_le hpq] at hvI
    exact h v ⟨hvI.1, lt_of_le_of_ne hvI.2 hv⟩
  have hae' : f =ᵐ[volume.restrict (Set.uIoc p q)] g := (ae_restrict_iff' measurableSet_uIoc).mpr hae
  exact ⟨hg.1.congr_ae hae'.symm, by rw [intervalIntegral.integral_congr_ae hae, hg.2]⟩

lemma tv_eq {v : ℝ} (hv0 : v ≠ 0) : tv v = (v - b) * (v - -a) / (4 * (v - 0)) := by
  unfold tv
  rw [sub_zero]
  field_simp
  simp only [Y01.a, Y01.b]
  ring_nf
  (try simp only [s2_2])
  (try ring_nf)

lemma one_sub_tv_eq {v : ℝ} (hv0 : v ≠ 0) :
    1 - tv v = -((v - a) * (v - -b) / (4 * (v - 0))) := by
  unfold tv
  rw [sub_zero]
  field_simp
  simp only [Y01.a, Y01.b]
  ring_nf
  (try simp only [s2_2])
  (try ring_nf)

lemma log_four : Real.log 4 = 2 * lam := by
  unfold lam
  rw [show (4:ℝ) = 2 ^ 2 by norm_num, Real.log_pow]
  norm_num

lemma log_tv {v : ℝ} (hv : v ∈ Ioo b a) :
    Real.log (tv v) = Real.log (v - b) + Real.log (v - -a) - Real.log (v - 0) - 2 * lam := by
  have hb0 : 0 < b := by simp only [Y01.b]; linarith [s2lo]
  have h1 : 0 < v - b := by linarith [hv.1]
  have h2 : 0 < v - -a := by linarith [hv.1, a_pos]
  have h3 : 0 < v - 0 := by linarith [hv.1]
  rw [tv_eq (ne_of_gt (by linarith [hv.1] : (0:ℝ) < v)), Real.log_div (mul_pos h1 h2).ne'
    (mul_pos (by norm_num) h3).ne',
    Real.log_mul h1.ne' h2.ne', Real.log_mul (by norm_num) h3.ne', log_four]
  ring

lemma log_one_sub_tv {v : ℝ} (hv : v ∈ Ioo b a) :
    Real.log (1 - tv v) = Real.log (v - a) + Real.log (v - -b) - Real.log (v - 0) - 2 * lam := by
  have hb0 : 0 < b := by simp only [Y01.b]; linarith [s2lo]
  have h1 : v - a < 0 := by linarith [hv.2]
  have h2 : 0 < v - -b := by linarith [hv.1]
  have h3 : 0 < v - 0 := by linarith [hv.1]
  rw [one_sub_tv_eq (ne_of_gt (by linarith [hv.1] : (0:ℝ) < v)), Real.log_neg_eq_log,
    Real.log_div (mul_ne_zero h1.ne h2.ne') (mul_pos (by norm_num) h3).ne',
    Real.log_mul h1.ne h2.ne', Real.log_mul (by norm_num) h3.ne', log_four]
  ring

noncomputable def S0 (v : ℝ) : ℝ :=
    (((2 : ℝ) + (-1 : ℝ) * Real.sqrt 2) * lam ^ 2) * fcst 0 v +
    ((((-2 : ℝ) + (5 / 4 : ℝ) * Real.sqrt 2) * lam) * flin b 0 v +
    ((((-2 : ℝ) + (5 / 4 : ℝ) * Real.sqrt 2) * lam) * flin (-a) 0 v +
    ((((2 : ℝ) + (-1 : ℝ) * Real.sqrt 2) * lam) * flin 0 0 v +
    ((((1 / 2 : ℝ) + (-5 / 16 : ℝ) * Real.sqrt 2)) * fsq b 0 v +
    ((((1 / 2 : ℝ) + (-5 / 16 : ℝ) * Real.sqrt 2)) * fsq (-a) 0 v +
    ((((-1 / 2 : ℝ) + (1 / 4 : ℝ) * Real.sqrt 2)) * fsq 0 0 v +
    ((((1 / 4 : ℝ) + (-1 / 8 : ℝ) * Real.sqrt 2)) * fpr (-a) 0 0 v +
    ((((1 / 4 : ℝ) + (-3 / 16 : ℝ) * Real.sqrt 2)) * fpr 0 (-a) 0 v +
    ((((-4 : ℝ) + (0 : ℝ) * Real.sqrt 2) * lam ^ 2) * fcst (-a) v +
    ((((4 : ℝ) + (0 : ℝ) * Real.sqrt 2) * lam) * flin b (-a) v +
    ((((4 : ℝ) + (0 : ℝ) * Real.sqrt 2) * lam) * flin (-a) (-a) v +
    ((((-4 : ℝ) + (0 : ℝ) * Real.sqrt 2) * lam) * flin 0 (-a) v +
    ((((-1 : ℝ) + (0 : ℝ) * Real.sqrt 2)) * fsq b (-a) v +
    ((((-1 : ℝ) + (0 : ℝ) * Real.sqrt 2)) * fsq (-a) (-a) v +
    ((((1 : ℝ) + (0 : ℝ) * Real.sqrt 2)) * fsq 0 (-a) v +
    ((((-1 / 2 : ℝ) + (0 : ℝ) * Real.sqrt 2)) * fpr (-a) 0 (-a) v +
    ((((-1 / 2 : ℝ) + (0 : ℝ) * Real.sqrt 2)) * fpr 0 (-a) (-a) v +
    ((((0 : ℝ) + (-1 / 4 : ℝ) * Real.sqrt 2) * lam) * flin a 0 v +
    ((((0 : ℝ) + (-1 / 4 : ℝ) * Real.sqrt 2) * lam) * flin (-b) 0 v +
    ((((0 : ℝ) + (1 / 16 : ℝ) * Real.sqrt 2)) * fsq a 0 v +
    ((((0 : ℝ) + (1 / 16 : ℝ) * Real.sqrt 2)) * fsq (-b) 0 v +
    ((((0 : ℝ) + (-1 / 16 : ℝ) * Real.sqrt 2)) * fpr b (-b) 0 v +
    ((((0 : ℝ) + (-1 / 16 : ℝ) * Real.sqrt 2)) * fpr (-a) (-b) 0 v +
    ((((0 : ℝ) + (1 / 16 : ℝ) * Real.sqrt 2)) * fpr 0 (-b) 0 v))))))))))))))))))))))))

noncomputable def V0 : ℝ :=
    (((2 : ℝ) + (-1 : ℝ) * Real.sqrt 2) * lam ^ 2) * (FK 0 a - FK 0 b) +
    ((((-2 : ℝ) + (5 / 4 : ℝ) * Real.sqrt 2) * lam) * (FL P1n b 0 a - FL P1n b 0 b) +
    ((((-2 : ℝ) + (5 / 4 : ℝ) * Real.sqrt 2) * lam) * (FL P1p (-a) 0 a - FL P1p (-a) 0 b) +
    ((((2 : ℝ) + (-1 : ℝ) * Real.sqrt 2) * lam) * (FH 0 a - FH 0 b) +
    ((((1 / 2 : ℝ) + (-5 / 16 : ℝ) * Real.sqrt 2)) * (FA P1n P2a b 0 a - FA P1n P2a b 0 b) +
    ((((1 / 2 : ℝ) + (-5 / 16 : ℝ) * Real.sqrt 2)) * (FA P1p P2g (-a) 0 a - FA P1p P2g (-a) 0 b) +
    ((((-1 / 2 : ℝ) + (1 / 4 : ℝ) * Real.sqrt 2)) * (FC 0 a - FC 0 b) +
    ((((1 / 4 : ℝ) + (-1 / 8 : ℝ) * Real.sqrt 2)) * (FIV P2g (-a) 0 a - FIV P2g (-a) 0 b) +
    ((((1 / 4 : ℝ) + (-3 / 16 : ℝ) * Real.sqrt 2)) * (FIII P2a 0 (-a) a - FIII P2a 0 (-a) b) +
    ((((-4 : ℝ) + (0 : ℝ) * Real.sqrt 2) * lam ^ 2) * (FK (-a) a - FK (-a) b) +
    ((((4 : ℝ) + (0 : ℝ) * Real.sqrt 2) * lam) * (FL P1n b (-a) a - FL P1n b (-a) b) +
    ((((4 : ℝ) + (0 : ℝ) * Real.sqrt 2) * lam) * (FH (-a) a - FH (-a) b) +
    ((((-4 : ℝ) + (0 : ℝ) * Real.sqrt 2) * lam) * (FL P1n 0 (-a) a - FL P1n 0 (-a) b) +
    ((((-1 : ℝ) + (0 : ℝ) * Real.sqrt 2)) * (FA P1n P2a b (-a) a - FA P1n P2a b (-a) b) +
    ((((-1 : ℝ) + (0 : ℝ) * Real.sqrt 2)) * (FC (-a) a - FC (-a) b) +
    ((((1 : ℝ) + (0 : ℝ) * Real.sqrt 2)) * (FA P1n P2a 0 (-a) a - FA P1n P2a 0 (-a) b) +
    ((((-1 / 2 : ℝ) + (0 : ℝ) * Real.sqrt 2)) * (FIII P2g (-a) 0 a - FIII P2g (-a) 0 b) +
    ((((-1 / 2 : ℝ) + (0 : ℝ) * Real.sqrt 2)) * (FIV P2a 0 (-a) a - FIV P2a 0 (-a) b) +
    ((((0 : ℝ) + (-1 / 4 : ℝ) * Real.sqrt 2) * lam) * (FL P1p a 0 a - FL P1p a 0 b) +
    ((((0 : ℝ) + (-1 / 4 : ℝ) * Real.sqrt 2) * lam) * (FL P1p (-b) 0 a - FL P1p (-b) 0 b) +
    ((((0 : ℝ) + (1 / 16 : ℝ) * Real.sqrt 2)) * (FA P1p P2a a 0 a - FA P1p P2a a 0 b) +
    ((((0 : ℝ) + (1 / 16 : ℝ) * Real.sqrt 2)) * (FA P1p P2g (-b) 0 a - FA P1p P2g (-b) 0 b) +
    ((((0 : ℝ) + (-1 / 16 : ℝ) * Real.sqrt 2)) * (FG P1n P2a P2a b (-b) 0 a - FG P1n P2a P2a b (-b) 0 b) +
    ((((0 : ℝ) + (-1 / 16 : ℝ) * Real.sqrt 2)) * (FG P1p P2a P2g (-a) (-b) 0 a - FG P1p P2a P2g (-a) (-b) 0 b) +
    ((((0 : ℝ) + (1 / 16 : ℝ) * Real.sqrt 2)) * (FIII P2a 0 (-b) a - FIII P2a 0 (-b) b)))))))))))))))))))))))))

lemma S0_int : IntervalIntegrable S0 volume b a ∧ ∫ v in b..a, S0 v = V0 := by
  have h := lin_cons (((2 : ℝ) + (-1 : ℝ) * Real.sqrt 2) * lam ^ 2) blk0
    (lin_cons (((-2 : ℝ) + (5 / 4 : ℝ) * Real.sqrt 2) * lam) blk1
    (lin_cons (((-2 : ℝ) + (5 / 4 : ℝ) * Real.sqrt 2) * lam) blk2
    (lin_cons (((2 : ℝ) + (-1 : ℝ) * Real.sqrt 2) * lam) blk3
    (lin_cons (((1 / 2 : ℝ) + (-5 / 16 : ℝ) * Real.sqrt 2)) blk4
    (lin_cons (((1 / 2 : ℝ) + (-5 / 16 : ℝ) * Real.sqrt 2)) blk5
    (lin_cons (((-1 / 2 : ℝ) + (1 / 4 : ℝ) * Real.sqrt 2)) blk7
    (lin_cons (((1 / 4 : ℝ) + (-1 / 8 : ℝ) * Real.sqrt 2)) blk9
    (lin_cons (((1 / 4 : ℝ) + (-3 / 16 : ℝ) * Real.sqrt 2)) blk10
    (lin_cons (((-4 : ℝ) + (0 : ℝ) * Real.sqrt 2) * lam ^ 2) blk15
    (lin_cons (((4 : ℝ) + (0 : ℝ) * Real.sqrt 2) * lam) blk16
    (lin_cons (((4 : ℝ) + (0 : ℝ) * Real.sqrt 2) * lam) blk17
    (lin_cons (((-4 : ℝ) + (0 : ℝ) * Real.sqrt 2) * lam) blk18
    (lin_cons (((-1 : ℝ) + (0 : ℝ) * Real.sqrt 2)) blk19
    (lin_cons (((-1 : ℝ) + (0 : ℝ) * Real.sqrt 2)) blk20
    (lin_cons (((1 : ℝ) + (0 : ℝ) * Real.sqrt 2)) blk22
    (lin_cons (((-1 / 2 : ℝ) + (0 : ℝ) * Real.sqrt 2)) blk24
    (lin_cons (((-1 / 2 : ℝ) + (0 : ℝ) * Real.sqrt 2)) blk25
    (lin_cons (((0 : ℝ) + (-1 / 4 : ℝ) * Real.sqrt 2) * lam) blk30
    (lin_cons (((0 : ℝ) + (-1 / 4 : ℝ) * Real.sqrt 2) * lam) blk31
    (lin_cons (((0 : ℝ) + (1 / 16 : ℝ) * Real.sqrt 2)) blk32
    (lin_cons (((0 : ℝ) + (1 / 16 : ℝ) * Real.sqrt 2)) blk34
    (lin_cons (((0 : ℝ) + (-1 / 16 : ℝ) * Real.sqrt 2)) blk35
    (lin_cons (((0 : ℝ) + (-1 / 16 : ℝ) * Real.sqrt 2)) blk37
    (lin_one (((0 : ℝ) + (1 / 16 : ℝ) * Real.sqrt 2)) blk39))))))))))))))))))))))))
  exact h

noncomputable def S1 (v : ℝ) : ℝ :=
    (((-1 / 2 : ℝ) + (3 / 8 : ℝ) * Real.sqrt 2)) * fpr b (-a) 0 v +
    ((((1 / 2 : ℝ) + (-5 / 16 : ℝ) * Real.sqrt 2)) * fpr b 0 0 v +
    ((((1 : ℝ) + (0 : ℝ) * Real.sqrt 2)) * fpr b (-a) (-a) v +
    ((((-1 : ℝ) + (0 : ℝ) * Real.sqrt 2)) * fpr b 0 (-a) v +
    ((((0 : ℝ) + (-1 / 16 : ℝ) * Real.sqrt 2)) * fpr b a 0 v +
    ((((0 : ℝ) + (-1 / 16 : ℝ) * Real.sqrt 2)) * fpr (-a) a 0 v +
    ((((0 : ℝ) + (1 / 16 : ℝ) * Real.sqrt 2)) * fpr 0 a 0 v))))))

noncomputable def V1 : ℝ :=
    (((-1 / 2 : ℝ) + (3 / 8 : ℝ) * Real.sqrt 2)) * (FG P1n P2a P2a b (-a) 0 1 - FG P1n P2a P2a b (-a) 0 b) +
    ((((1 / 2 : ℝ) + (-5 / 16 : ℝ) * Real.sqrt 2)) * (FIV P2a b 0 1 - FIV P2a b 0 b) +
    ((((1 : ℝ) + (0 : ℝ) * Real.sqrt 2)) * (FIV P2a b (-a) 1 - FIV P2a b (-a) b) +
    ((((-1 : ℝ) + (0 : ℝ) * Real.sqrt 2)) * (FG P1p P2a P2a b 0 (-a) 1 - FG P1p P2a P2a b 0 (-a) b) +
    ((((0 : ℝ) + (-1 / 16 : ℝ) * Real.sqrt 2)) * (FG P1n P2a P2a b a 0 1 - FG P1n P2a P2a b a 0 b) +
    ((((0 : ℝ) + (-1 / 16 : ℝ) * Real.sqrt 2)) * (FG P1p P2g P2a (-a) a 0 1 - FG P1p P2g P2a (-a) a 0 b) +
    ((((0 : ℝ) + (1 / 16 : ℝ) * Real.sqrt 2)) * (FIII P2a 0 a 1 - FIII P2a 0 a b)))))))

lemma S1_int : IntervalIntegrable S1 volume b 1 ∧ ∫ v in b..1, S1 v = V1 := by
  have h := lin_cons (((-1 / 2 : ℝ) + (3 / 8 : ℝ) * Real.sqrt 2)) blk6
    (lin_cons (((1 / 2 : ℝ) + (-5 / 16 : ℝ) * Real.sqrt 2)) blk8
    (lin_cons (((1 : ℝ) + (0 : ℝ) * Real.sqrt 2)) blk21
    (lin_cons (((-1 : ℝ) + (0 : ℝ) * Real.sqrt 2)) blk23
    (lin_cons (((0 : ℝ) + (-1 / 16 : ℝ) * Real.sqrt 2)) blk33
    (lin_cons (((0 : ℝ) + (-1 / 16 : ℝ) * Real.sqrt 2)) blk36
    (lin_one (((0 : ℝ) + (1 / 16 : ℝ) * Real.sqrt 2)) blk38))))))
  exact h

noncomputable def S2 (v : ℝ) : ℝ :=
    (((-1 / 4 : ℝ) + (3 / 16 : ℝ) * Real.sqrt 2)) * fpr b (-a) 0 v +
    ((((1 / 4 : ℝ) + (-1 / 8 : ℝ) * Real.sqrt 2)) * fpr b 0 0 v +
    ((((-1 / 4 : ℝ) + (3 / 16 : ℝ) * Real.sqrt 2)) * fpr (-a) b 0 v +
    ((((1 / 4 : ℝ) + (-3 / 16 : ℝ) * Real.sqrt 2)) * fpr 0 b 0 v +
    ((((1 / 2 : ℝ) + (0 : ℝ) * Real.sqrt 2)) * fpr b (-a) (-a) v +
    ((((-1 / 2 : ℝ) + (0 : ℝ) * Real.sqrt 2)) * fpr b 0 (-a) v +
    ((((1 / 2 : ℝ) + (0 : ℝ) * Real.sqrt 2)) * fpr (-a) b (-a) v +
    ((((-1 / 2 : ℝ) + (0 : ℝ) * Real.sqrt 2)) * fpr 0 b (-a) v +
    ((((0 : ℝ) + (-1 / 16 : ℝ) * Real.sqrt 2)) * fpr a b 0 v +
    ((((0 : ℝ) + (-1 / 16 : ℝ) * Real.sqrt 2)) * fpr a (-a) 0 v +
    ((((0 : ℝ) + (1 / 16 : ℝ) * Real.sqrt 2)) * fpr a 0 0 v))))))))))

noncomputable def V2 : ℝ :=
    (((-1 / 4 : ℝ) + (3 / 16 : ℝ) * Real.sqrt 2)) * (FG P1n P2a P2a b (-a) 0 a - FG P1n P2a P2a b (-a) 0 1) +
    ((((1 / 4 : ℝ) + (-1 / 8 : ℝ) * Real.sqrt 2)) * (FIV P2a b 0 a - FIV P2a b 0 1) +
    ((((-1 / 4 : ℝ) + (3 / 16 : ℝ) * Real.sqrt 2)) * (FG P1n P2a P2g (-a) b 0 a - FG P1n P2a P2g (-a) b 0 1) +
    ((((1 / 4 : ℝ) + (-3 / 16 : ℝ) * Real.sqrt 2)) * (FIII P2g 0 b a - FIII P2g 0 b 1) +
    ((((1 / 2 : ℝ) + (0 : ℝ) * Real.sqrt 2)) * (FIV P2a b (-a) a - FIV P2a b (-a) 1) +
    ((((-1 / 2 : ℝ) + (0 : ℝ) * Real.sqrt 2)) * (FG P1p P2a P2a b 0 (-a) a - FG P1p P2a P2a b 0 (-a) 1) +
    ((((1 / 2 : ℝ) + (0 : ℝ) * Real.sqrt 2)) * (FIII P2g (-a) b a - FIII P2g (-a) b 1) +
    ((((-1 / 2 : ℝ) + (0 : ℝ) * Real.sqrt 2)) * (FG P1p P2g P2g 0 b (-a) a - FG P1p P2g P2g 0 b (-a) 1) +
    ((((0 : ℝ) + (-1 / 16 : ℝ) * Real.sqrt 2)) * (FG P1n P2a P2a a b 0 a - FG P1n P2a P2a a b 0 1) +
    ((((0 : ℝ) + (-1 / 16 : ℝ) * Real.sqrt 2)) * (FG P1p P2a P2a a (-a) 0 a - FG P1p P2a P2a a (-a) 0 1) +
    ((((0 : ℝ) + (1 / 16 : ℝ) * Real.sqrt 2)) * (FIV P2a a 0 a - FIV P2a a 0 1)))))))))))

lemma S2_int : IntervalIntegrable S2 volume 1 a ∧ ∫ v in 1..a, S2 v = V2 := by
  have h := lin_cons (((-1 / 4 : ℝ) + (3 / 16 : ℝ) * Real.sqrt 2)) blk11
    (lin_cons (((1 / 4 : ℝ) + (-1 / 8 : ℝ) * Real.sqrt 2)) blk12
    (lin_cons (((-1 / 4 : ℝ) + (3 / 16 : ℝ) * Real.sqrt 2)) blk13
    (lin_cons (((1 / 4 : ℝ) + (-3 / 16 : ℝ) * Real.sqrt 2)) blk14
    (lin_cons (((1 / 2 : ℝ) + (0 : ℝ) * Real.sqrt 2)) blk26
    (lin_cons (((-1 / 2 : ℝ) + (0 : ℝ) * Real.sqrt 2)) blk27
    (lin_cons (((1 / 2 : ℝ) + (0 : ℝ) * Real.sqrt 2)) blk28
    (lin_cons (((-1 / 2 : ℝ) + (0 : ℝ) * Real.sqrt 2)) blk29
    (lin_cons (((0 : ℝ) + (-1 / 16 : ℝ) * Real.sqrt 2)) blk40
    (lin_cons (((0 : ℝ) + (-1 / 16 : ℝ) * Real.sqrt 2)) blk41
    (lin_one (((0 : ℝ) + (1 / 16 : ℝ) * Real.sqrt 2)) blk42))))))))))
  exact h


lemma pw1 : ∀ v ∈ Ioo b 1, W3integrand v = S0 v + S1 v := by
  intro v hv
  have hv' : v ∈ Ioo b a := ⟨hv.1, lt_of_lt_of_le hv.2 h1a⟩
  unfold W3integrand l1 l2
  rw [log_tv hv', log_one_sub_tv hv']
  unfold S0 S1 fsq fpr flin fcst
  simp only [sub_zero, sub_neg_eq_add]
  ring


lemma pw2 : ∀ v ∈ Ioo 1 a, W3integrand v = S0 v + S2 v := by
  intro v hv
  have hv' : v ∈ Ioo b a := ⟨lt_of_le_of_lt hb1 hv.1, hv.2⟩
  unfold W3integrand l1 l2
  rw [log_tv hv', log_one_sub_tv hv']
  unfold S0 S2 fsq fpr flin fcst
  simp only [sub_zero, sub_neg_eq_add]
  ring

set_option maxHeartbeats 4000000 in
lemma final_eq : V0 + V1 + V2 = rawW3 := by
  unfold V0 V1 V2
  simp only [val0q, val0p, val1q, val1p, val2q, val2p, val3q, val3p, val4q, val4p, val5q, val5p, val6q, val6p, val7q, val7p, val8q, val8p, val9q, val9p, val10q, val10p, val11q, val11p, val12q, val12p, val13q, val13p, val14q, val14p, val15q, val15p, val16q, val16p, val17q, val17p, val18q, val18p, val19q, val19p, val20q, val20p, val21q, val21p, val22q, val22p, val23q, val23p, val24q, val24p, val25q, val25p, val26q, val26p, val27q, val27p, val28q, val28p, val29q, val29p, val30q, val30p, val31q, val31p, val32q, val32p, val33q, val33p, val34q, val34p, val35q, val35p, val36q, val36p, val37q, val37p, val38q, val38p, val39q, val39p, val40q, val40p, val41q, val41p, val42q, val42p]
  unfold rawW3
  ring

theorem main : W3 = rawW3 := by
  obtain ⟨i0, e0⟩ := S0_int
  obtain ⟨i1, e1⟩ := S1_int
  obtain ⟨i2, e2⟩ := S2_int
  have i0a : IntervalIntegrable S0 volume b 1 := i0.mono_set (by
    rw [uIcc_of_le hba, uIcc_of_le hb1]; exact Icc_subset_Icc le_rfl h1a)
  have i0b : IntervalIntegrable S0 volume 1 a := i0.mono_set (by
    rw [uIcc_of_le hba, uIcc_of_le h1a]; exact Icc_subset_Icc hb1 le_rfl)
  have k1 := congr_Ioo (A := (∫ v in b..1, S0 v) + V1) hb1 pw1
    ⟨i0a.add i1, by rw [intervalIntegral.integral_add i0a i1, e1]⟩
  have k2 := congr_Ioo (A := (∫ v in (1:ℝ)..a, S0 v) + V2) h1a pw2
    ⟨i0b.add i2, by rw [intervalIntegral.integral_add i0b i2, e2]⟩
  have hs := intervalIntegral.integral_add_adjacent_intervals i0a i0b
  have hw := intervalIntegral.integral_add_adjacent_intervals k1.1 k2.1
  unfold W3
  rw [← hw, k1.2, k2.2]
  have h3 : (∫ v in b..1, S0 v) + V1 + ((∫ v in (1:ℝ)..a, S0 v) + V2) = V0 + V1 + V2 := by
    rw [← e0, ← hs]; ring
  rw [h3]
  exact final_eq

end PE1
end Y01

/-! # Part E2: `rawW3 = M2closed - (π² + 6λ²)/24` (Lemmas g1:lem:dilog, g1:lem:trilog) -/


namespace Y01.PE2

open Real

/-! ## Basic facts about `√2`, `a`, `b` -/

lemma hs : Real.sqrt 2 ^ 2 = 2 := Real.sq_sqrt (by norm_num)
lemma hs3 : Real.sqrt 2 ^ 3 = 2 * Real.sqrt 2 := by rw [pow_succ, hs]
lemma hs4 : Real.sqrt 2 ^ 4 = 4 := by
  rw [show (4:ℕ) = 2 * 2 by norm_num, pow_mul, hs]; norm_num
lemma hs5 : Real.sqrt 2 ^ 5 = 4 * Real.sqrt 2 := by rw [pow_succ, hs4]
lemma hs6 : Real.sqrt 2 ^ 6 = 8 := by
  rw [show Real.sqrt 2 ^ 6 = (Real.sqrt 2 ^ 2) ^ 3 by ring, hs]; norm_num
lemma hs0 : 0 < Real.sqrt 2 := by positivity
lemma hslo : 1.41 < Real.sqrt 2 := by
  rw [Real.lt_sqrt (by norm_num)]; norm_num
lemma hshi : Real.sqrt 2 < 1.42 := by
  rw [Real.sqrt_lt' (by norm_num)]; norm_num
lemma hsne : Real.sqrt 2 ≠ 0 := hs0.ne'

lemma ha : a = 1 + Real.sqrt 2 := rfl
lemma hb : b = Real.sqrt 2 - 1 := rfl
lemma a_pos : 0 < a := by rw [ha]; linarith [hs0]
lemma b_pos : 0 < b := by rw [hb]; linarith [hslo]
lemma b_lt : b < 1 := by rw [hb]; linarith [hshi]
lemma ab : a * b = 1 := by rw [ha, hb]; linear_combination hs
lemma b_eq_inv : b = a⁻¹ := eq_inv_of_mul_eq_one_right ab

/-! ## Logarithms -/

lemma log_a : Real.log a = L := rfl
lemma log_two : Real.log 2 = lam := rfl
lemma log_b : Real.log b = -L := by rw [b_eq_inv, Real.log_inv]; rfl
lemma log_s : Real.log (Real.sqrt 2) = lam / 2 := by
  rw [Real.log_sqrt (by norm_num)]; rfl
lemma log_is : Real.log (1 / Real.sqrt 2) = -(lam / 2) := by
  rw [one_div, Real.log_inv, log_s]
lemma log_half : Real.log (1 / 2) = -lam := by
  rw [one_div, Real.log_inv]; rfl
lemma log_mul' {x y : ℝ} (hx : 0 < x) (hy : 0 < y) : Real.log (x * y) = Real.log x + Real.log y :=
  Real.log_mul hx.ne' hy.ne'

lemma log_2ms : Real.log (2 - Real.sqrt 2) = lam / 2 - L := by
  rw [show 2 - Real.sqrt 2 = Real.sqrt 2 * b by rw [hb]; linear_combination -hs,
    log_mul' hs0 b_pos, log_s, log_b]; ring
lemma log_1mis : Real.log (1 - 1 / Real.sqrt 2) = -L - lam / 2 := by
  rw [show 1 - 1 / Real.sqrt 2 = b * (1 / Real.sqrt 2) by
      rw [hb]; field_simp,
    log_mul' b_pos (by positivity), log_b, log_is]; ring
lemma log_2b : Real.log (2 * b) = lam - L := by
  rw [log_mul' (by norm_num) b_pos, log_b]; rfl
lemma log_b2 : Real.log (b ^ 2) = -2 * L := by
  rw [Real.log_pow, log_b]; push_cast; ring
lemma log_a2 : Real.log (a ^ 2) = 2 * L := by
  rw [Real.log_pow]; push_cast; rfl
lemma log_2a : Real.log (2 * a) = lam + L := by
  rw [log_mul' (by norm_num) a_pos]; rfl
lemma log_p4 : Real.log ((2 + Real.sqrt 2) / 4) = L - 3 * lam / 2 := by
  rw [show (2 + Real.sqrt 2) / 4 = a * (1 / Real.sqrt 2) * (1 / 2) by
      rw [ha]; field_simp; linear_combination 2 * hs,
    log_mul' (mul_pos a_pos (by positivity)) (by norm_num), log_mul' a_pos (by positivity),
    log_is, log_half, log_a]; ring
lemma log_m4 : Real.log ((2 - Real.sqrt 2) / 4) = -L - 3 * lam / 2 := by
  rw [show (2 - Real.sqrt 2) / 4 = b * (1 / Real.sqrt 2) * (1 / 2) by
      rw [hb]; field_simp; linear_combination -2 * hs,
    log_mul' (mul_pos b_pos (by positivity)) (by norm_num), log_mul' b_pos (by positivity),
    log_is, log_half, log_b]; ring

end Y01.PE2

namespace Y01.PE2

open Real

/-! ## Wrappers around the dilogarithm relations of D00 -/

lemma is_eq : 1 / Real.sqrt 2 = Real.sqrt 2 / 2 := by
  field_simp; rw [hs]

lemma div_eq_of {p q r : ℝ} (hq : q ≠ 0) (h : p = r * q) : p / q = r := by
  rw [div_eq_iff hq]; exact h

lemma W_refl {x y c d : ℝ} (h0 : 0 < x) (h1 : x < 1) (hy : 1 - x = y) (hc : Real.log x = c)
    (hd : Real.log y = d) : D00.Li2 x + D00.Li2 y = π ^ 2 / 6 - c * d := by
  have h := D00.reflection h0 h1
  rw [hy, hc, hd] at h; exact h

lemma W_landen {x y c : ℝ} (hx : x < 1) (hy : x / (x - 1) = y) (hc : Real.log (1 - x) = c) :
    D00.Li2 x + D00.Li2 y = -(1 / 2) * c ^ 2 := by
  have h := D00.landen hx
  rw [hy, hc] at h; exact h

lemma W_inv {x y c : ℝ} (hx : x < 0) (hy : 1 / x = y) (hc : Real.log (-x) = c) :
    D00.Li2 x + D00.Li2 y = -(π ^ 2 / 6) - (1 / 2) * c ^ 2 := by
  have h := D00.inversion_neg hx
  rw [hy, hc] at h; exact h

lemma W_five {x y u v w c d : ℝ} (hx : x < 1) (hy : y < 1) (hxy : x + y < 1)
    (hu : x / (1 - y) = u) (hv : y / (1 - x) = v) (hw : x * y / ((1 - x) * (1 - y)) = w)
    (hc : Real.log (1 - x) = c) (hd : Real.log (1 - y) = d) :
    D00.Li2 u + D00.Li2 v - D00.Li2 w = D00.Li2 x + D00.Li2 y + c * d := by
  have h := D00.five_term hx hy hxy
  rw [hu, hv, hw, hc, hd] at h; exact h

/-! ## Table g1:tab:dilog, one lemma per line -/

lemma d_mb : D00.Li2 (-b) = D00.Li2 b - π ^ 2 / 8 + L ^ 2 / 2 := by
  have h := D00.chi2_landen b_pos b_lt
  rw [div_eq_of (p := 1 - b) (q := 1 + b) (r := b) (by linarith [b_pos]) (by rw [hb]; linear_combination (-1:ℝ) * hs),
    div_eq_of (p := 1 + b) (q := 1 - b) (r := a) (by linarith [b_lt]) (by rw [ha, hb]; linear_combination hs),
    log_b, log_a] at h
  unfold D00.chi2 at h
  linarith

lemma d_b2 : D00.Li2 (b ^ 2) = 4 * D00.Li2 b - π ^ 2 / 4 + L ^ 2 := by
  have h := D00.duplication (x := b) (by rw [abs_of_pos b_pos]; exact b_lt.le)
  linarith [d_mb]

lemma d_2ms : D00.Li2 (2 - Real.sqrt 2) = -D00.Li2 b + π ^ 2 / 6 - L ^ 2 + lam * L / 2 := by
  have h := W_refl (x := 2 - Real.sqrt 2) (y := b) (by linarith [hshi]) (by linarith [hslo])
    (by rw [hb]; ring) log_2ms log_b
  linarith

lemma d_ms : D00.Li2 (-Real.sqrt 2)
    = D00.Li2 b - π ^ 2 / 6 + L ^ 2 / 2 - lam * L / 2 := by
  have h := W_landen (x := -Real.sqrt 2) (y := 2 - Real.sqrt 2) (c := L) (by linarith [hs0])
    (div_eq_of (by linarith [hs0]) (by linear_combination (-1:ℝ) * hs))
    (by rw [show 1 - -Real.sqrt 2 = a by rw [ha]; ring, log_a])
  linarith [d_2ms]

lemma d_mb2 : D00.Li2 (-(b ^ 2)) = 6 * D00.Li2 b - 11 / 24 * π ^ 2 + 2 * L ^ 2 := by
  have h := W_five (x := -Real.sqrt 2) (y := -b) (u := -1) (v := -(b ^ 2)) (w := b ^ 2)
    (c := L) (d := lam / 2) (by linarith [hs0]) (by linarith [b_pos])
    (by linarith [hs0, b_pos])
    (div_eq_of (by rw [hb]; linarith [hs0]) (by rw [hb]; linear_combination (0:ℝ) * hs))
    (div_eq_of (by linarith [hs0]) (by rw [hb]; linear_combination (Real.sqrt 2 - 1) * hs))
    (div_eq_of (by rw [hb]; exact (mul_pos (by linarith [hs0]) (by linarith [hs0])).ne')
      (by rw [hb]; linear_combination (Real.sqrt 2 - Real.sqrt 2 ^ 2) * hs))
    (by rw [show 1 - -Real.sqrt 2 = a by rw [ha]; ring, log_a])
    (by rw [show 1 - -b = Real.sqrt 2 by rw [hb]; ring]; exact log_s)
  rw [D00.Li2_neg_one] at h
  linarith [d_ms, d_mb, d_b2]

lemma d_ma : D00.Li2 (-a) = -D00.Li2 b - π ^ 2 / 24 - L ^ 2 := by
  have h := W_inv (x := -a) (y := -b) (c := L) (by linarith [a_pos])
    (by rw [one_div, inv_neg, ← b_eq_inv]) (by rw [neg_neg]; rfl)
  linarith [d_mb]

lemma d_ma2 : D00.Li2 (-(a ^ 2)) = -6 * D00.Li2 b + 7 / 24 * π ^ 2 - 4 * L ^ 2 := by
  have h := W_inv (x := -(a ^ 2)) (y := -(b ^ 2)) (c := 2 * L) (by nlinarith [a_pos])
    (by rw [one_div, inv_neg, ← inv_pow, ← b_eq_inv]) (by rw [neg_neg, log_a2])
  linarith [d_mb2]

lemma d_is : D00.Li2 (1 / Real.sqrt 2)
    = D00.Li2 b + π ^ 2 / 24 + L ^ 2 / 2 - lam * L / 2 - lam ^ 2 / 8 := by
  have h := W_landen (x := -a) (y := 1 / Real.sqrt 2) (c := lam / 2 + L) (by linarith [a_pos])
    (div_eq_of (by linarith [a_pos]) (by rw [is_eq, ha]; linear_combination (1 / 2 : ℝ) * hs))
    (by rw [show 1 - -a = Real.sqrt 2 * a by rw [ha]; linear_combination (-1:ℝ) * hs,
      log_mul' hs0 a_pos, log_s]; rfl)
  linarith [d_ma]

lemma d_mis : D00.Li2 (-(1 / Real.sqrt 2))
    = -D00.Li2 b - L ^ 2 / 2 + lam * L / 2 - lam ^ 2 / 8 := by
  have h := W_landen (x := -(1 / Real.sqrt 2)) (y := b) (c := L - lam / 2)
    (by have : 0 < 1 / Real.sqrt 2 := by positivity
        linarith)
    (div_eq_of (by have : 0 < 1 / Real.sqrt 2 := by positivity
                   linarith)
      (by rw [is_eq, hb]; linear_combination (1 / 2 : ℝ) * hs))
    (by rw [show 1 - -(1 / Real.sqrt 2) = a * (1 / Real.sqrt 2) by
          rw [is_eq, ha]; linear_combination (-1 / 2 : ℝ) * hs,
        log_mul' a_pos (by positivity), log_is]; rw [log_a]; ring)
  linarith

lemma d_1mis : D00.Li2 (1 - 1 / Real.sqrt 2)
    = -D00.Li2 b + π ^ 2 / 8 - L ^ 2 / 2 - lam ^ 2 / 8 := by
  have h := W_refl (x := 1 / Real.sqrt 2) (y := 1 - 1 / Real.sqrt 2) (by positivity)
    (by rw [div_lt_one hs0]; linarith [hslo]) rfl log_is log_1mis
  linarith [d_is]

lemma d_2b : D00.Li2 (2 * b) = -4 * D00.Li2 b + 5 / 12 * π ^ 2 - 3 * L ^ 2 + 2 * lam * L := by
  have h := W_refl (x := 2 * b) (y := b ^ 2) (by linarith [b_pos]) (by rw [hb]; linarith [hshi])
    (by rw [hb]; linear_combination (-1:ℝ) * hs) log_2b log_b2
  linarith [d_b2]

lemma d_mbh : D00.Li2 (-(b / 2))
    = -4 * D00.Li2 b + π ^ 2 / 4 - 3 / 2 * L ^ 2 + lam * L - lam ^ 2 / 2 := by
  have h := W_landen (x := -(b / 2)) (y := b ^ 2) (c := L - lam) (by linarith [b_pos])
    (div_eq_of (by linarith [b_pos]) (by rw [hb]; linear_combination (Real.sqrt 2 - 1) / 2 * hs))
    (by rw [show 1 - -(b / 2) = a * (1 / 2) by rw [ha, hb]; ring,
        log_mul' a_pos (by norm_num), log_half]; rw [log_a]; ring)
  linarith [d_b2]

lemma d_m2a : D00.Li2 (-(2 * a)) = 4 * D00.Li2 b - 5 / 12 * π ^ 2 + L ^ 2 - 2 * lam * L := by
  have h := W_inv (x := -(2 * a)) (y := -(b / 2)) (c := lam + L) (by linarith [a_pos])
    (by rw [b_eq_inv]; field_simp) (by rw [neg_neg, log_2a])
  linarith [d_mbh]

lemma d_p4 : D00.Li2 ((2 + Real.sqrt 2) / 4)
    = 6 * D00.Li2 b - 7 / 24 * π ^ 2 + 7 / 2 * L ^ 2 - 3 / 2 * lam * L - 9 / 8 * lam ^ 2 := by
  have h := W_landen (x := (2 + Real.sqrt 2) / 4) (y := -(a ^ 2)) (c := -L - 3 * lam / 2)
    (by linarith [hshi])
    (div_eq_of (by linarith [hshi]) (by rw [ha]; linear_combination Real.sqrt 2 / 4 * hs))
    (by rw [show 1 - (2 + Real.sqrt 2) / 4 = (2 - Real.sqrt 2) / 4 by ring, log_m4])
  linarith [d_ma2]

lemma d_m4 : D00.Li2 ((2 - Real.sqrt 2) / 4)
    = -6 * D00.Li2 b + 11 / 24 * π ^ 2 - 5 / 2 * L ^ 2 + 3 / 2 * lam * L - 9 / 8 * lam ^ 2 := by
  have h := W_refl (x := (2 + Real.sqrt 2) / 4) (y := (2 - Real.sqrt 2) / 4)
    (by positivity) (by linarith [hshi]) (by ring) log_p4 log_m4
  linarith [d_p4]

lemma d_m1 : D00.Li2 (-1) = -(π ^ 2 / 12) := D00.Li2_neg_one

lemma d_h : D00.Li2 (1 / 2) = π ^ 2 / 12 - lam ^ 2 / 2 := by
  rw [D00.Li2_half, log_two]; ring

end Y01.PE2

namespace Y01.PE2

open Real

/-! ## Wrappers around the trilogarithm relations (I3), (Tr), (D3), (K) -/

lemma one_sub_inv_eq {z w : ℝ} (hz : z ≠ 0) (h : z - 1 = w * z) : 1 - 1 / z = w := by
  rw [show 1 - 1 / z = (z - 1) / z by field_simp]
  exact div_eq_of hz h

lemma W_I3 {z w c : ℝ} (hz : z < 0) (hw : 1 / z = w) (hc : Real.log (-z) = c) :
    Li3 z - Li3 w = -(π ^ 2 / 6) * c - (1 / 6) * c ^ 3 := by
  have h := I3 hz
  rw [hw, hc] at h; exact h

lemma W_Tr {z u w c d : ℝ} (h0 : 0 < z) (h1 : z < 1) (hu : 1 - z = u) (hw : 1 - 1 / z = w)
    (hc : Real.log z = c) (hd : Real.log (1 - z) = d) :
    Li3 z + Li3 u + Li3 w = z3 + (1 / 6) * c ^ 3 + π ^ 2 / 6 * c - (1 / 2) * c ^ 2 * d := by
  have h := Tr h0 h1
  rw [hd, hu, hw, hc] at h; exact h

lemma W_D3 {w v : ℝ} (hw : |w| ≤ 1) (hv : w ^ 2 = v) : Li3 v = 4 * Li3 w + 4 * Li3 (-w) := by
  have h := D3 hw
  rw [hv] at h; exact h

lemma W_K {x y k1 k2 k3 k4 k5 k6 k7 ly lr : ℝ} (hx0 : 0 < x) (hxy : x < y) (hy1 : y < 1)
    (h1 : x * (1 - y) ^ 2 / (y * (1 - x) ^ 2) = k1) (h2 : x * y = k2) (h3 : x / y = k3)
    (h4 : x * (1 - y) / (y * (1 - x)) = k4) (h5 : x * (1 - y) / (x - 1) = k5)
    (h6 : (1 - y) / (1 - x) = k6) (h7 : (1 - y) / (y * (x - 1)) = k7)
    (hly : Real.log y = ly) (hlr : Real.log k6 = lr) :
    Li3 k1 + Li3 k2 + Li3 k3 - 2 * Li3 k4 - 2 * Li3 k5 - 2 * Li3 k6 - 2 * Li3 k7
      - 2 * Li3 x - 2 * Li3 y + 2 * z3
      = ly ^ 2 * lr - π ^ 2 / 3 * ly - (1 / 3) * ly ^ 3 := by
  have h := Y01.kummer hx0 hxy hy1
  simp only [Y01.KD, Y01.KA1, Y01.KA2, Y01.KA3, Y01.KA4, Y01.KA5, Y01.KA6, Y01.KA7] at h
  replace h : Li3 (x * (1 - y) ^ 2 / (y * (1 - x) ^ 2)) + Li3 (x * y) + Li3 (x / y)
      - 2 * Li3 (x * (1 - y) / (y * (1 - x))) - 2 * Li3 (x * (1 - y) / (x - 1))
      - 2 * Li3 ((1 - y) / (1 - x))
      - 2 * Li3 ((1 - y) / (y * (x - 1))) - 2 * Li3 x - 2 * Li3 y + 2 * z3
      = Real.log y ^ 2 * Real.log ((1 - y) / (1 - x)) - π ^ 2 / 3 * Real.log y
        - (1 / 3) * Real.log y ^ 3 := by linarith
  rw [h1, h2, h3, h4, h5, h6, h7, hly, hlr] at h; exact h

lemma is_pos : 0 < 1 / Real.sqrt 2 := by positivity
lemma is_lt : 1 / Real.sqrt 2 < 1 := by rw [div_lt_one hs0]; linarith [hslo]
lemma is_lo : 0.7 < 1 / Real.sqrt 2 := by rw [is_eq]; linarith [hslo]
lemma is_hi : 1 / Real.sqrt 2 < 0.71 := by rw [is_eq]; linarith [hshi]

/-! ## Table g1:tab:trilog, one lemma per line (plus `Li₃(-1)`, `Li₃(½)`, `Li₃(b²)`) -/

lemma t_m1 : Li3 (-1) = -(3 / 4) * z3 := by
  have h := W_D3 (w := 1) (v := 1) (by norm_num) (by norm_num)
  rw [Li3_one] at h; linarith

lemma t_h : Li3 (1 / 2) = 7 / 8 * z3 - π ^ 2 * lam / 12 + lam ^ 3 / 6 := by
  have h := W_Tr (z := 1 / 2) (u := 1 / 2) (w := -1) (c := -lam) (d := -lam)
    (by norm_num) (by norm_num) (by norm_num) (by norm_num) log_half
    (by rw [show (1:ℝ) - 1 / 2 = 1 / 2 by norm_num, log_half])
  linarith [t_m1]

lemma t_b2 : Li3 (b ^ 2) = 4 * Li3 b + 4 * Li3 (-b) :=
  W_D3 (by rw [abs_of_pos b_pos]; exact b_lt.le) rfl

lemma t_ma : Li3 (-a) = Li3 (-b) - π ^ 2 * L / 6 - L ^ 3 / 6 := by
  have h := W_I3 (z := -a) (w := -b) (c := L) (by linarith [a_pos])
    (by rw [one_div, inv_neg, ← b_eq_inv]) (by rw [neg_neg, log_a])
  linarith

lemma t_mis : Li3 (-(1 / Real.sqrt 2))
    = -Li3 (1 / Real.sqrt 2) - π ^ 2 * lam / 48 + 7 / 32 * z3 + lam ^ 3 / 24 := by
  have h := W_D3 (w := 1 / Real.sqrt 2) (v := 1 / 2)
    (by rw [abs_of_pos is_pos]; exact is_lt.le)
    (by rw [is_eq]; linear_combination (1 / 4 : ℝ) * hs)
  linarith [t_h]

lemma t_1mis : Li3 (1 - 1 / Real.sqrt 2)
    = -Li3 (-b) - Li3 (1 / Real.sqrt 2) - π ^ 2 * lam / 12 + z3 + lam ^ 2 * L / 8
      + lam ^ 3 / 24 := by
  have h := W_Tr (z := 1 / Real.sqrt 2) (u := 1 - 1 / Real.sqrt 2) (w := -b)
    (c := -(lam / 2)) (d := -L - lam / 2) is_pos is_lt rfl
    (by rw [one_div_one_div, hb]; ring) log_is log_1mis
  linarith

lemma t_2ms : Li3 (2 - Real.sqrt 2)
    = -Li3 b + Li3 (1 / Real.sqrt 2) - π ^ 2 * L / 6 + 5 / 48 * π ^ 2 * lam + 25 / 32 * z3
      + L ^ 3 / 3 - lam * L ^ 2 / 4 - lam ^ 3 / 48 := by
  have h := W_Tr (z := 2 - Real.sqrt 2) (u := b) (w := -(1 / Real.sqrt 2))
    (c := lam / 2 - L) (d := -L) (by linarith [hshi]) (by linarith [hslo])
    (by rw [hb]; ring)
    (one_sub_inv_eq (by linarith [hshi]) (by rw [is_eq]; linear_combination (-1 / 2 : ℝ) * hs))
    log_2ms (by rw [show 1 - (2 - Real.sqrt 2) = b by rw [hb]; ring, log_b])
  linarith [t_mis]

lemma t_2b : Li3 (2 * b)
    = -6 * Li3 b + 2 * Li3 (-b) + 8 * Li3 (1 / Real.sqrt 2) - 2 / 3 * π ^ 2 * L
      + 7 / 12 * π ^ 2 * lam - 3 / 8 * z3 + 4 / 3 * L ^ 3 - lam * L ^ 2 - lam ^ 3 / 6 := by
  have h := W_K (x := 2 - Real.sqrt 2) (y := 1 / Real.sqrt 2) (k1 := b) (k2 := b) (k3 := 2 * b)
    (k4 := 2 - Real.sqrt 2) (k5 := -b) (k6 := 1 / Real.sqrt 2) (k7 := -1)
    (ly := -(lam / 2)) (lr := -(lam / 2))
    (by linarith [hshi]) (by linarith [hslo, is_lo]) is_lt
    (div_eq_of (by rw [is_eq]; intro h; nlinarith [hs, hslo, hshi])
      (by rw [is_eq, hb]; linear_combination (-Real.sqrt 2 ^ 2 / 2 + 5 * Real.sqrt 2 / 4 - 1) * hs))
    (by rw [is_eq, hb]; linear_combination (-1 / 2 : ℝ) * hs)
    (div_eq_of is_pos.ne' (by rw [is_eq, hb]; linear_combination (-1 : ℝ) * hs))
    (div_eq_of (by rw [is_eq]; intro h; nlinarith [hs, hslo, hshi])
      (by rw [is_eq]; linear_combination (Real.sqrt 2 / 2 - 1) * hs))
    (div_eq_of (by intro h; nlinarith [hs, hslo, hshi])
      (by rw [is_eq, hb]; linear_combination (-1 / 2 : ℝ) * hs))
    (div_eq_of (by intro h; nlinarith [hs, hslo, hshi])
      (by rw [is_eq]; linear_combination (-1 / 2 : ℝ) * hs))
    (div_eq_of (by rw [is_eq]; intro h; nlinarith [hs, hslo, hshi])
      (by rw [is_eq]; linear_combination (-1 / 2 : ℝ) * hs))
    log_is log_is
  linarith [t_2ms, t_m1]

lemma t_m2a : Li3 (-(2 * a))
    = 2 * Li3 b - 6 * Li3 (-b) - 8 * Li3 (1 / Real.sqrt 2) + π ^ 2 * L / 3
      - 7 / 12 * π ^ 2 * lam + 11 / 8 * z3 - 2 / 3 * L ^ 3 - lam * L ^ 2 + lam ^ 3 / 6 := by
  have h := W_Tr (z := b ^ 2) (u := 2 * b) (w := -(2 * a)) (c := -2 * L) (d := lam - L)
    (pow_pos b_pos 2) (by nlinarith [b_pos, b_lt])
    (by rw [hb]; linear_combination (-1 : ℝ) * hs)
    (by rw [one_div, b_eq_inv, inv_pow, inv_inv, ha]; linear_combination (-1 : ℝ) * hs)
    log_b2 (by rw [show 1 - b ^ 2 = 2 * b by rw [hb]; linear_combination (-1 : ℝ) * hs, log_2b])
  linarith [t_b2, t_2b]

lemma t_mbh : Li3 (-(b / 2))
    = 2 * Li3 b - 6 * Li3 (-b) - 8 * Li3 (1 / Real.sqrt 2) + π ^ 2 * L / 2
      - 5 / 12 * π ^ 2 * lam + 11 / 8 * z3 - L ^ 3 / 2 - lam * L ^ 2 / 2 + lam ^ 2 * L / 2
      + lam ^ 3 / 3 := by
  have h := W_I3 (z := -(2 * a)) (w := -(b / 2)) (c := lam + L) (by linarith [a_pos])
    (by rw [b_eq_inv]; field_simp) (by rw [neg_neg, log_2a])
  linarith [t_m2a]

lemma t_m4 : Li3 ((2 - Real.sqrt 2) / 4)
    = 7 * Li3 b - 13 * Li3 (-b) - 18 * Li3 (1 / Real.sqrt 2) + 7 / 6 * π ^ 2 * L
      - 17 / 16 * π ^ 2 * lam + 101 / 32 * z3 - 4 / 3 * L ^ 3 - 3 / 4 * lam * L ^ 2
      + 9 / 8 * lam ^ 2 * L + 15 / 16 * lam ^ 3 := by
  have h := W_K (x := b) (y := 1 / Real.sqrt 2) (k1 := (2 - Real.sqrt 2) / 4)
    (k2 := 1 - 1 / Real.sqrt 2) (k3 := 2 - Real.sqrt 2) (k4 := 1 - 1 / Real.sqrt 2)
    (k5 := -(b / 2)) (k6 := 1 / 2) (k7 := -(1 / Real.sqrt 2))
    (ly := -(lam / 2)) (lr := -lam)
    b_pos (by rw [hb]; linarith [hshi, is_lo]) is_lt
    (div_eq_of (by rw [is_eq, hb]; intro h; nlinarith [hs, hslo, hshi])
      (by rw [is_eq, hb]; linear_combination (Real.sqrt 2 ^ 2 / 8 - Real.sqrt 2 / 2 + 1 / 2) * hs))
    (by rw [is_eq, hb]; linear_combination (1 / 2 : ℝ) * hs)
    (div_eq_of is_pos.ne' (by rw [is_eq, hb]; linear_combination (1 / 2 : ℝ) * hs))
    (div_eq_of (by rw [is_eq, hb]; intro h; nlinarith [hs, hslo, hshi])
      (by rw [is_eq, hb]; linear_combination (1 / 2 - Real.sqrt 2 / 4) * hs))
    (div_eq_of (by rw [hb]; intro h; nlinarith [hs, hslo, hshi])
      (by rw [is_eq, hb]; ring))
    (div_eq_of (by rw [hb]; intro h; nlinarith [hs, hslo, hshi])
      (by rw [is_eq, hb]; ring))
    (div_eq_of (by rw [is_eq, hb]; intro h; nlinarith [hs, hslo, hshi])
      (by rw [is_eq, hb]; linear_combination (Real.sqrt 2 / 4 - 1 / 2) * hs))
    log_is log_half
  linarith [t_1mis, t_2ms, t_mbh, t_h, t_mis]

lemma t_p4 : Li3 ((2 + Real.sqrt 2) / 4)
    = -13 * Li3 b + 7 * Li3 (-b) + 18 * Li3 (1 / Real.sqrt 2) - 7 / 6 * π ^ 2 * L
      + 13 / 16 * π ^ 2 * lam - 25 / 32 * z3 + 7 / 3 * L ^ 3 - 3 / 4 * lam * L ^ 2
      - 9 / 8 * lam ^ 2 * L + 3 / 16 * lam ^ 3 := by
  have h := W_K (x := 1 / Real.sqrt 2) (y := 2 * b) (k1 := 1 - 1 / Real.sqrt 2)
    (k2 := 2 - Real.sqrt 2) (k3 := (2 + Real.sqrt 2) / 4) (k4 := 1 / 2)
    (k5 := -b) (k6 := 2 - Real.sqrt 2) (k7 := -(1 / Real.sqrt 2))
    (ly := lam - L) (lr := lam / 2 - L)
    is_pos (by rw [hb]; linarith [hslo, is_hi]) (by rw [hb]; linarith [hshi])
    (div_eq_of (by rw [is_eq, hb]; intro h; nlinarith [hs, hslo, hshi])
      (by rw [is_eq, hb]; linear_combination (Real.sqrt 2 ^ 2 / 4 + Real.sqrt 2 / 4 - 1) * hs))
    (by rw [is_eq, hb]; linear_combination (1 : ℝ) * hs)
    (div_eq_of (by rw [hb]; intro h; nlinarith [hs, hslo, hshi])
      (by rw [is_eq, hb]; linear_combination (-1 / 2 : ℝ) * hs))
    (div_eq_of (by rw [is_eq, hb]; intro h; nlinarith [hs, hslo, hshi])
      (by rw [is_eq, hb]; linear_combination (-1 / 2 : ℝ) * hs))
    (div_eq_of (by rw [is_eq]; intro h; nlinarith [hs, hslo, hshi])
      (by rw [is_eq, hb]; linear_combination (-1 / 2 : ℝ) * hs))
    (div_eq_of (by rw [is_eq]; intro h; nlinarith [hs, hslo, hshi])
      (by rw [is_eq, hb]; linear_combination (-1 / 2 : ℝ) * hs))
    (div_eq_of (by rw [is_eq, hb]; intro h; nlinarith [hs, hslo, hshi])
      (by rw [is_eq, hb]; linear_combination (Real.sqrt 2 / 2 - 3 / 2) * hs))
    log_2b log_2ms
  linarith [t_1mis, t_2ms, t_h, t_2b, t_mis]

lemma t_mb2 : Li3 (-(b ^ 2))
    = 6 * Li3 b + 6 * Li3 (-b) + π ^ 2 * L / 6 - 11 / 8 * z3 - L ^ 3 / 3 := by
  have h := W_Tr (z := (2 + Real.sqrt 2) / 4) (u := (2 - Real.sqrt 2) / 4) (w := -(b ^ 2))
    (c := L - 3 * lam / 2) (d := -L - 3 * lam / 2) (by positivity) (by linarith [hshi])
    (by ring)
    (one_sub_inv_eq (by positivity) (by rw [hb]; linear_combination Real.sqrt 2 / 4 * hs))
    log_p4 (by rw [show 1 - (2 + Real.sqrt 2) / 4 = (2 - Real.sqrt 2) / 4 by ring, log_m4])
  linarith [t_p4, t_m4]

end Y01.PE2

namespace Y01.PE2

open Real

/-! ## The `ℚ(√2)` bookkeeping: substitute every line and compare with (T1) -/

theorem main : Y01.rawW3 = Y01.M2closed - (Real.pi ^ 2 + 6 * Y01.lam ^ 2) / 24 := by
  unfold rawW3 M2closed
  rw [t_m2a, t_ma, t_m1, t_mis, t_mbh, t_mb2, t_m4, t_1mis, t_h, t_2ms, t_2b, t_p4,
    d_ma2, d_m2a, d_ma, d_m1, d_mis, d_mb, d_mbh, d_mb2, d_m4, d_b2, d_1mis, d_h, d_2ms, d_is,
    d_2b, d_p4]
  ring

end Y01.PE2

/-!
# Y01 Part F: numerical value of `M_2` (stmt_12 → stmt_13)

Enclosures: `√2`, `λ = log 2`, `L = log(1+√2)`, `w = Li₂(1/2)`, `u = Li₂(b)` (copied from D03),
`T = Li₃(1/2)`, `P = Li₃(b)`, `Q = Li₃(b²)`, `O = Σ_k 2^{-k}/(2k+1)³` (geometric tails).
Reductions: `π² = 12w + 6λ²`, `ζ(3) = (8/7)T + (8/7)wλ + (8/21)λ³` (Tr at 1/2 and D3 at 1),
`Li₃(-b) = Q/4 - P` (D3 at b), `Li₃(1/√2) = T/8 + √2·O/2` (D3 at 1/√2 and the odd part).
-/

namespace Y01.PF

open Real

lemma s_sq : Real.sqrt 2 ^ 2 = 2 := Real.sq_sqrt (by norm_num)
lemma s_gt : (1.41 : ℝ) < Real.sqrt 2 := by
  rw [Real.lt_sqrt (by norm_num)]; norm_num
lemma s_lt : Real.sqrt 2 < 1.42 := by
  rw [Real.sqrt_lt' (by norm_num)]; norm_num
lemma s_pos : 0 < Real.sqrt 2 := by linarith [s_gt]

lemma log_s : Real.log (Real.sqrt 2) = Y01.lam / 2 := by
  rw [Real.log_sqrt (by norm_num)]; rfl

lemma ab : (Real.sqrt 2 - 1) * (1 + Real.sqrt 2) = 1 := by nlinarith [s_sq]

lemma log_b : Real.log (Real.sqrt 2 - 1) = -Y01.L := by
  have h := Real.log_mul (x := Real.sqrt 2 - 1) (y := 1 + Real.sqrt 2)
    (by linarith [s_gt]) (by linarith [s_gt])
  rw [ab, Real.log_one] at h
  unfold Y01.L Y01.a; linarith

lemma log_two_s : Real.log (2 - Real.sqrt 2) = -Y01.L + Y01.lam / 2 := by
  have e2 : 2 - Real.sqrt 2 = (Real.sqrt 2 - 1) * Real.sqrt 2 := by nlinarith [s_sq]
  rw [e2, Real.log_mul (by linarith [s_gt]) (by linarith [s_gt]), log_b, log_s]

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

lemma li2_bounds {x lo hi : ℝ} (h0 : 0 ≤ lo) (hlo : lo ≤ x) (hhi : x ≤ hi) (hhi1 : hi < 1)
    (N : ℕ) :
    ∑ n ∈ Finset.range N, lo ^ (n + 1) / ((n : ℝ) + 1) ^ 2 ≤ D00.Li2 x ∧
    D00.Li2 x ≤ ∑ n ∈ Finset.range N, hi ^ (n + 1) / ((n : ℝ) + 1) ^ 2
      + hi ^ (N + 1) / ((N : ℝ) + 1) ^ 2 / (1 - hi) := by
  have hx0 : 0 ≤ x := h0.trans hlo
  have hS := D00.hasSum_Li2 (x := x) (by rw [abs_le]; constructor <;> linarith)
  have hf0 : ∀ n : ℕ, 0 ≤ x ^ (n + 1) / ((n : ℝ) + 1) ^ 2 := fun n => by positivity
  have hb := hasSum_bounds hS hf0 N (hi ^ (N + 1) / ((N : ℝ) + 1) ^ 2) hi
    (h0.trans (hlo.trans hhi)) hhi1 (fun n => by
      have e1 : x ^ (n + N + 1) ≤ hi ^ (n + N + 1) := pow_le_pow_left₀ hx0 hhi _
      have e2 : ((N : ℝ) + 1) ^ 2 ≤ (((n + N : ℕ) : ℝ) + 1) ^ 2 := by
        push_cast; apply pow_le_pow_left₀ (by positivity)
        linarith [(Nat.cast_nonneg n : (0:ℝ) ≤ n)]
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

lemma li3_bounds {x lo hi : ℝ} (h0 : 0 ≤ lo) (hlo : lo ≤ x) (hhi : x ≤ hi) (hhi1 : hi < 1)
    (N : ℕ) :
    ∑ n ∈ Finset.range N, lo ^ (n + 1) / ((n : ℝ) + 1) ^ 3 ≤ Y01.Li3 x ∧
    Y01.Li3 x ≤ ∑ n ∈ Finset.range N, hi ^ (n + 1) / ((n : ℝ) + 1) ^ 3
      + hi ^ (N + 1) / ((N : ℝ) + 1) ^ 3 / (1 - hi) := by
  have hx0 : 0 ≤ x := h0.trans hlo
  have hS := Y01.Li3_hasSum (x := x) (by rw [abs_le]; constructor <;> linarith)
  have hf0 : ∀ n : ℕ, 0 ≤ x ^ (n + 1) / ((n : ℝ) + 1) ^ 3 := fun n => by positivity
  have hb := hasSum_bounds hS hf0 N (hi ^ (N + 1) / ((N : ℝ) + 1) ^ 3) hi
    (h0.trans (hlo.trans hhi)) hhi1 (fun n => by
      have e1 : x ^ (n + N + 1) ≤ hi ^ (n + N + 1) := pow_le_pow_left₀ hx0 hhi _
      have e2 : ((N : ℝ) + 1) ^ 3 ≤ (((n + N : ℕ) : ℝ) + 1) ^ 3 := by
        push_cast; apply pow_le_pow_left₀ (by positivity)
        linarith [(Nat.cast_nonneg n : (0:ℝ) ≤ n)]
      calc x ^ (n + N + 1) / (((n + N : ℕ) : ℝ) + 1) ^ 3
          ≤ hi ^ (n + N + 1) / ((N : ℝ) + 1) ^ 3 := by
            apply div_le_div₀ (pow_nonneg (h0.trans (hlo.trans hhi)) _) e1 (by positivity) e2
        _ = hi ^ (N + 1) / ((N : ℝ) + 1) ^ 3 * hi ^ n := by ring)
  refine ⟨?_, ?_⟩
  · refine le_trans ?_ hb.1
    apply Finset.sum_le_sum
    intro n _
    apply div_le_div_of_nonneg_right (pow_le_pow_left₀ h0 hlo _) (by positivity)
  · refine le_trans hb.2 ?_
    gcongr with n hn

lemma art_bounds {x lo hi : ℝ} (h0 : 0 ≤ lo) (hlo : lo ≤ x) (hhi : x ≤ hi) (hhi1 : hi < 1)
    (N : ℕ) :
    ∑ k ∈ Finset.range N, 2 * (1 / (2 * (k : ℝ) + 1)) * lo ^ (2 * k + 1)
      ≤ Real.log (1 + x) - Real.log (1 - x) ∧
    Real.log (1 + x) - Real.log (1 - x)
      ≤ ∑ k ∈ Finset.range N, 2 * (1 / (2 * (k : ℝ) + 1)) * hi ^ (2 * k + 1)
        + 2 * hi ^ (2 * N + 1) / (2 * (N : ℝ) + 1) / (1 - hi ^ 2) := by
  have hx0 : 0 ≤ x := h0.trans hlo
  have hS := Real.hasSum_log_sub_log_of_abs_lt_one (x := x)
    (by rw [abs_lt]; constructor <;> linarith)
  have hf0 : ∀ k : ℕ, 0 ≤ 2 * (1 / (2 * (k : ℝ) + 1)) * x ^ (2 * k + 1) := fun k => by
    positivity
  have hb := hasSum_bounds hS hf0 N (2 * hi ^ (2 * N + 1) / (2 * (N : ℝ) + 1)) (hi ^ 2)
    (by positivity) (by nlinarith [h0.trans (hlo.trans hhi)]) (fun n => by
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
      Finset.sum_le_sum (fun n _ =>
        mul_le_mul_of_nonneg_left (pow_le_pow_left₀ hx0 hhi _) (by positivity))
    linarith

lemma s_bounds : (1.41421356237309504880168872 : ℝ) ≤ Real.sqrt 2 ∧
    Real.sqrt 2 ≤ (1.41421356237309504880168873 : ℝ) := by
  constructor
  · apply le_of_lt; rw [Real.lt_sqrt (by norm_num)]; norm_num
  · apply le_of_lt; rw [Real.sqrt_lt' (by norm_num)]; norm_num

lemma u_bounds : (0.46753399702300461762696479 : ℝ) ≤ D00.Li2 (Real.sqrt 2 - 1) ∧
    D00.Li2 (Real.sqrt 2 - 1) ≤ (0.46753399702300461762696481 : ℝ) := by
  have hs := s_bounds
  have h := li2_bounds (x := Real.sqrt 2 - 1) (lo := (0.41421356237309504880168872 : ℝ))
    (hi := (0.41421356237309504880168873 : ℝ)) (by norm_num)
    (by linarith [hs.1]) (by linarith [hs.2]) (by norm_num) 60
  constructor
  · refine le_trans ?_ h.1
    norm_num [Finset.sum_range_succ]
  · refine le_trans h.2 ?_
    norm_num [Finset.sum_range_succ]

lemma w_bounds : (0.58224052646501250590265631 : ℝ) ≤ D00.Li2 (1/2) ∧
    D00.Li2 (1/2) ≤ (0.58224052646501250590265633 : ℝ) := by
  have h := li2_bounds (x := 1/2) (lo := 1/2) (hi := 1/2) (by norm_num) le_rfl le_rfl
    (by norm_num) 75
  constructor
  · refine le_trans ?_ h.1
    norm_num [Finset.sum_range_succ]
  · refine le_trans h.2 ?_
    norm_num [Finset.sum_range_succ]

lemma lam_bounds : (0.69314718055994530941723212 : ℝ) ≤ Y01.lam ∧
    Y01.lam ≤ (0.69314718055994530941723213 : ℝ) := by
  have h := art_bounds (x := 1/3) (lo := 1/3) (hi := 1/3) (by norm_num) le_rfl le_rfl
    (by norm_num) 28
  have e : Real.log (1 + 1/3) - Real.log (1 - 1/3) = Y01.lam := by
    rw [← Real.log_div (by norm_num) (by norm_num)]; unfold Y01.lam; norm_num
  rw [e] at h
  constructor
  · refine le_trans ?_ h.1
    norm_num [Finset.sum_range_succ]
  · refine le_trans h.2 ?_
    norm_num [Finset.sum_range_succ]

lemma L_bounds : (0.88137358701954302523260931 : ℝ) ≤ Y01.L ∧
    Y01.L ≤ (0.88137358701954302523260934 : ℝ) := by
  have hs := s_bounds
  have h := art_bounds (x := Real.sqrt 2 - 1) (lo := (0.41421356237309504880168872 : ℝ))
    (hi := (0.41421356237309504880168873 : ℝ)) (by norm_num)
    (by linarith [hs.1]) (by linarith [hs.2]) (by norm_num) 36
  have e : Real.log (1 + (Real.sqrt 2 - 1)) - Real.log (1 - (Real.sqrt 2 - 1)) = Y01.L := by
    rw [show 1 + (Real.sqrt 2 - 1) = Real.sqrt 2 by ring,
      show 1 - (Real.sqrt 2 - 1) = 2 - Real.sqrt 2 by ring, log_s, log_two_s]; ring
  rw [e] at h
  constructor
  · refine le_trans ?_ h.1
    norm_num [Finset.sum_range_succ]
  · refine le_trans h.2 ?_
    norm_num [Finset.sum_range_succ]

lemma T_bounds : (0.537213193608040200940623182388 : ℝ) ≤ Y01.Li3 (1/2) ∧
    Y01.Li3 (1/2) ≤ (0.537213193608040200940623227450 : ℝ) := by
  have h := li3_bounds (x := 1/2) (lo := 1/2) (hi := 1/2) (by norm_num) le_rfl le_rfl
    (by norm_num) 66
  constructor
  · refine le_trans ?_ h.1
    norm_num [Finset.sum_range_succ]
  · refine le_trans h.2 ?_
    norm_num [Finset.sum_range_succ]

lemma P_bounds : (0.438881700078672395657262704418 : ℝ) ≤ Y01.Li3 (Real.sqrt 2 - 1) ∧
    Y01.Li3 (Real.sqrt 2 - 1) ≤ (0.438881700078672395657263102955 : ℝ) := by
  have hs := s_bounds
  have h := li3_bounds (x := Real.sqrt 2 - 1) (lo := (0.41421356237309504880168872 : ℝ))
    (hi := (0.41421356237309504880168873 : ℝ)) (by norm_num)
    (by linarith [hs.1]) (by linarith [hs.2]) (by norm_num) 50
  constructor
  · refine le_trans ?_ h.1
    norm_num [Finset.sum_range_succ]
  · refine le_trans h.2 ?_
    norm_num [Finset.sum_range_succ]

lemma b_sq : (Real.sqrt 2 - 1) ^ 2 = 3 - 2 * Real.sqrt 2 := by nlinarith [s_sq]

lemma Q_bounds : (0.175454454019511266772526112476 : ℝ) ≤ Y01.Li3 ((Real.sqrt 2 - 1) ^ 2) ∧
    Y01.Li3 ((Real.sqrt 2 - 1) ^ 2) ≤ (0.175454454019511266772526153583 : ℝ) := by
  have hs := s_bounds
  have h := li3_bounds (x := (Real.sqrt 2 - 1) ^ 2)
    (lo := (0.17157287525380990239662254 : ℝ))
    (hi := (0.17157287525380990239662256 : ℝ)) (by norm_num)
    (by rw [b_sq]; linarith [hs.2]) (by rw [b_sq]; linarith [hs.1]) (by norm_num) 27
  constructor
  · refine le_trans ?_ h.1
    norm_num [Finset.sum_range_succ]
  · refine le_trans h.2 ?_
    norm_num [Finset.sum_range_succ]

end Y01.PF

namespace Y01.PF

open Real

/-! ## The odd series `O = Σ_k 2^{-k}/(2k+1)³` -/

noncomputable def Og (k : ℕ) : ℝ := (1 / 2 : ℝ) ^ k / (2 * (k : ℝ) + 1) ^ 3

lemma Og_nonneg (k : ℕ) : 0 ≤ Og k := by unfold Og; positivity

lemma Og_le (k : ℕ) : Og k ≤ (1 / 2 : ℝ) ^ k := by
  unfold Og
  apply div_le_self (by positivity)
  have : (1:ℝ) ≤ 2 * (k : ℝ) + 1 := by linarith [(Nat.cast_nonneg k : (0:ℝ) ≤ k)]
  exact one_le_pow₀ this

lemma Og_summable : Summable Og :=
  Summable.of_nonneg_of_le Og_nonneg Og_le
    (summable_geometric_of_lt_one (by norm_num) (by norm_num))

noncomputable def O : ℝ := ∑' k, Og k

lemma O_hasSum : HasSum Og O := Og_summable.hasSum

lemma Og_tail (N n : ℕ) :
    Og (n + N) ≤ (1 / 2 : ℝ) ^ N / (2 * (N : ℝ) + 1) ^ 3 * (1 / 2 : ℝ) ^ n := by
  unfold Og
  rw [pow_add, mul_comm ((1 / 2 : ℝ) ^ n), div_mul_eq_mul_div]
  apply div_le_div_of_nonneg_left (by positivity) (by positivity)
  apply pow_le_pow_left₀ (by positivity)
  push_cast; linarith [(Nat.cast_nonneg n : (0:ℝ) ≤ n)]

lemma O_bounds : (1.021002840767543829163857617042 : ℝ) ≤ O ∧
    O ≤ (1.021002840767543829163857839087 : ℝ) := by
  have h := hasSum_bounds O_hasSum Og_nonneg 62 _ (1 / 2) (by norm_num) (by norm_num)
    (Og_tail 62)
  constructor
  · refine le_trans ?_ h.1
    norm_num [Finset.sum_range_succ, Og]
  · refine le_trans h.2 ?_
    norm_num [Finset.sum_range_succ, Og]

/-! ## Reductions -/

lemma pi2_eq : π ^ 2 = 12 * D00.Li2 (1 / 2) + 6 * Y01.lam ^ 2 := by
  have h := D00.Li2_half
  unfold Y01.lam
  linear_combination (-12) * h

lemma z3_id : Y01.z3 = 8 / 7 * Y01.Li3 (1 / 2) + 8 / 7 * D00.Li2 (1 / 2) * Y01.lam
    + 8 / 21 * Y01.lam ^ 3 := by
  have hTr := Y01.Tr (z := 1 / 2) (by norm_num) (by norm_num)
  have hl : Real.log (1 / 2) = -Y01.lam := by rw [one_div, Real.log_inv]; rfl
  rw [show (1:ℝ) - 1 / 2 = 1 / 2 by norm_num, show (1:ℝ) - 1 / (1 / 2) = -1 by norm_num,
    hl] at hTr
  have hD : Y01.Li3 ((1:ℝ) ^ 2) = 4 * Y01.Li3 1 + 4 * Y01.Li3 (-1) := Y01.D3 (by norm_num)
  rw [one_pow, Y01.Li3_one] at hD
  have hpi := pi2_eq
  linear_combination (-4 / 7) * hTr + (-1 / 7) * hD + (2 / 21 * Y01.lam) * hpi

lemma Li3_negb : Y01.Li3 (-(Real.sqrt 2 - 1)) =
    Y01.Li3 ((Real.sqrt 2 - 1) ^ 2) / 4 - Y01.Li3 (Real.sqrt 2 - 1) := by
  have h := Y01.D3 (w := Real.sqrt 2 - 1) (by
    rw [abs_le]; constructor <;> linarith [s_gt, s_lt])
  linarith

lemma inv_s : 1 / Real.sqrt 2 = Real.sqrt 2 / 2 := by
  have h := s_pos.ne'
  field_simp
  rw [s_sq]

lemma Li3_half_sqrt : Y01.Li3 (Real.sqrt 2 / 2) =
    Y01.Li3 (1 / 2) / 8 + Real.sqrt 2 * O / 2 := by
  set x : ℝ := Real.sqrt 2 / 2 with hxdef
  have hx2 : x ^ 2 = 1 / 2 := by rw [hxdef, div_pow, s_sq]; norm_num
  have hx0 : 0 ≤ x := by rw [hxdef]; positivity
  have hx1 : x ≤ 1 := by rw [hxdef]; linarith [s_lt]
  have habs : |x| ≤ 1 := by rw [abs_le]; constructor <;> linarith
  have habs' : |-x| ≤ 1 := by rw [abs_neg]; exact habs
  have hD := Y01.D3 habs
  rw [hx2] at hD
  have hf := (Y01.Li3_hasSum habs).sub (Y01.Li3_hasSum habs')
  have hinj : Function.Injective (fun k : ℕ => 2 * k) := fun a b h => by
    simp only at h; omega
  have hzero : ∀ n : ℕ, n ∉ Set.range (fun k : ℕ => 2 * k) →
      x ^ (n + 1) / ((n : ℝ) + 1) ^ 3 - (-x) ^ (n + 1) / ((n : ℝ) + 1) ^ 3 = 0 := by
    intro n hn
    rcases Nat.even_or_odd' n with ⟨k, rfl | rfl⟩
    · exact absurd ⟨k, rfl⟩ hn
    · rw [show 2 * k + 1 + 1 = 2 * (k + 1) by ring, pow_mul, pow_mul, neg_sq, sub_self]
  have hc := (hinj.hasSum_iff hzero).mpr hf
  have heq : ((fun n : ℕ => x ^ (n + 1) / ((n : ℝ) + 1) ^ 3
      - (-x) ^ (n + 1) / ((n : ℝ) + 1) ^ 3) ∘ (fun k : ℕ => 2 * k))
      = fun k => Real.sqrt 2 * Og k := by
    funext k
    simp only [Function.comp, Og]
    have e1 : x ^ (2 * k + 1) = x * (1 / 2) ^ k := by
      rw [pow_succ, pow_mul, hx2, mul_comm]
    have e2 : (-x) ^ (2 * k + 1) = -(x * (1 / 2) ^ k) := by
      rw [Odd.neg_pow ⟨k, rfl⟩, e1]
    rw [e1, e2]
    push_cast
    rw [hxdef]
    ring
  rw [heq] at hc
  have hu := hc.unique (O_hasSum.mul_left (Real.sqrt 2))
  linarith

end Y01.PF

namespace Y01.PF

open Real

lemma mb {x y xl xh yl yh : ℝ} (hx : xl ≤ x ∧ x ≤ xh) (hy : yl ≤ y ∧ y ≤ yh)
    (h1 : 0 ≤ xl) (h2 : 0 ≤ yl) : xl * yl ≤ x * y ∧ x * y ≤ xh * yh :=
  ⟨mul_le_mul hx.1 hy.1 h2 (h1.trans hx.1),
    mul_le_mul hx.2 hy.2 (h2.trans hy.1) (h1.trans (hx.1.trans hx.2))⟩

set_option maxHeartbeats 4000000 in
lemma final_abs (s l L w u T P Q O : ℝ) (hs2 : s ^ 2 = 2)
    (hs : (1.41421356237309504880168872 : ℝ) ≤ s ∧ s ≤ (1.41421356237309504880168873 : ℝ))
    (hu : (0.46753399702300461762696479 : ℝ) ≤ u ∧ u ≤ (0.46753399702300461762696481 : ℝ))
    (hw : (0.58224052646501250590265631 : ℝ) ≤ w ∧ w ≤ (0.58224052646501250590265633 : ℝ))
    (hl : (0.69314718055994530941723212 : ℝ) ≤ l ∧ l ≤ (0.69314718055994530941723213 : ℝ))
    (hL : (0.88137358701954302523260931 : ℝ) ≤ L ∧ L ≤ (0.88137358701954302523260934 : ℝ))
    (hT : (0.537213193608040200940623182388 : ℝ) ≤ T ∧ T ≤ (0.537213193608040200940623227450 : ℝ))
    (hP : (0.438881700078672395657262704418 : ℝ) ≤ P ∧ P ≤ (0.438881700078672395657263102955 : ℝ))
    (hQ : (0.175454454019511266772526112476 : ℝ) ≤ Q ∧ Q ≤ (0.175454454019511266772526153583 : ℝ))
    (hO : (1.021002840767543829163857617042 : ℝ) ≤ O ∧ O ≤ (1.021002840767543829163857839087 : ℝ)) :
    |(((12 * w + 6 * l ^ 2) + 6 * l ^ 2) / 24
      + (189 * s / 64 - 19 / 16) * (8 / 7 * T + 8 / 7 * w * l + 8 / 21 * l ^ 3)
      + (11 + 15 * s) * P + (11 - 15 * s) * (Q / 4 - P)
      - 27 * s * (T / 8 + s * O / 2)
      + (12 * L - 7 * s / 2 * l) * u
      + (15 * s / 8 - 7 / 12) * (12 * w + 6 * l ^ 2) * L
      - (1 / 24 + 65 * s / 48) * (12 * w + 6 * l ^ 2) * l
      + (1 / 6 + 9 * s / 16) * l ^ 3
      + 11 * s / 8 * l ^ 2 * L
      - (1 + 7 * s / 4) * l * L ^ 2
      + (25 / 6 - 7 * s / 3) * L ^ 3)
      - (-0.02216191108726939942 : ℝ)| < (1 : ℝ) / 10 ^ 20 := by
  have hsO : s ^ 2 * O = 2 * O := by rw [hs2]
  have m1 := (mb hL hu (by norm_num) (by norm_num))
  have m2 := (mb hL hw (by norm_num) (by norm_num))
  have m3 := (mb (mb hL hL (by norm_num) (by norm_num)) hL (by norm_num) (by norm_num))
  have m4 := (mb hl hw (by norm_num) (by norm_num))
  have m5 := (mb (mb hl hL (by norm_num) (by norm_num)) hL (by norm_num) (by norm_num))
  have m6 := (mb hl hl (by norm_num) (by norm_num))
  have m7 := (mb (mb hl hl (by norm_num) (by norm_num)) hL (by norm_num) (by norm_num))
  have m8 := (mb (mb hl hl (by norm_num) (by norm_num)) hl (by norm_num) (by norm_num))
  have m9 := (mb hs hQ (by norm_num) (by norm_num))
  have m10 := (mb hs hP (by norm_num) (by norm_num))
  have m11 := (mb (mb hs hL (by norm_num) (by norm_num)) hw (by norm_num) (by norm_num))
  have m12 := (mb (mb (mb hs hL (by norm_num) (by norm_num)) hL (by norm_num) (by norm_num)) hL (by norm_num) (by norm_num))
  have m13 := (mb (mb hs hl (by norm_num) (by norm_num)) hu (by norm_num) (by norm_num))
  have m14 := (mb (mb hs hl (by norm_num) (by norm_num)) hw (by norm_num) (by norm_num))
  have m15 := (mb (mb (mb hs hl (by norm_num) (by norm_num)) hL (by norm_num) (by norm_num)) hL (by norm_num) (by norm_num))
  have m16 := (mb (mb (mb hs hl (by norm_num) (by norm_num)) hl (by norm_num) (by norm_num)) hL (by norm_num) (by norm_num))
  have m17 := (mb (mb (mb hs hl (by norm_num) (by norm_num)) hl (by norm_num) (by norm_num)) hl (by norm_num) (by norm_num))

  rw [abs_lt]
  constructor
  · linarith
  · linarith

lemma M2closed_eq : Y01.M2closed =
    ((12 * D00.Li2 (1 / 2) + 6 * Y01.lam ^ 2) + 6 * Y01.lam ^ 2) / 24
      + (189 * Real.sqrt 2 / 64 - 19 / 16) * (8 / 7 * Y01.Li3 (1 / 2) + 8 / 7 * D00.Li2 (1 / 2) * Y01.lam + 8 / 21 * Y01.lam ^ 3)
      + (11 + 15 * Real.sqrt 2) * Y01.Li3 (Real.sqrt 2 - 1) + (11 - 15 * Real.sqrt 2) * (Y01.Li3 ((Real.sqrt 2 - 1) ^ 2) / 4 - Y01.Li3 (Real.sqrt 2 - 1))
      - 27 * Real.sqrt 2 * (Y01.Li3 (1 / 2) / 8 + Real.sqrt 2 * O / 2)
      + (12 * Y01.L - 7 * Real.sqrt 2 / 2 * Y01.lam) * D00.Li2 (Real.sqrt 2 - 1)
      + (15 * Real.sqrt 2 / 8 - 7 / 12) * (12 * D00.Li2 (1 / 2) + 6 * Y01.lam ^ 2) * Y01.L
      - (1 / 24 + 65 * Real.sqrt 2 / 48) * (12 * D00.Li2 (1 / 2) + 6 * Y01.lam ^ 2) * Y01.lam
      + (1 / 6 + 9 * Real.sqrt 2 / 16) * Y01.lam ^ 3
      + 11 * Real.sqrt 2 / 8 * Y01.lam ^ 2 * Y01.L
      - (1 + 7 * Real.sqrt 2 / 4) * Y01.lam * Y01.L ^ 2
      + (25 / 6 - 7 * Real.sqrt 2 / 3) * Y01.L ^ 3 := by
  unfold Y01.M2closed Y01.b
  rw [Li3_negb, inv_s, Li3_half_sqrt, z3_id, pi2_eq]

theorem main : Y01.stmt_12 → Y01.stmt_13 := by
  intro h
  unfold Y01.stmt_13
  rw [h, ← Complex.ofReal_sub, Complex.norm_real, Real.norm_eq_abs, M2closed_eq]
  exact final_abs _ _ _ _ _ _ _ _ _ s_sq s_bounds u_bounds w_bounds lam_bounds L_bounds
    T_bounds P_bounds Q_bounds O_bounds

end Y01.PF


theorem Y01.proof_1 : Y01.stmt_1 := fun x _ => by
  apply (Y01.continuous_H.intervalIntegrable 0 x).congr_ae
  apply MeasureTheory.ae_restrict_of_ae
  filter_upwards [(Set.countable_singleton (0:ℝ)).ae_notMem MeasureTheory.volume] with t ht
  exact Y01.H_of_ne ht
theorem Y01.proof_2 : Y01.stmt_2 := fun _ hx => Y01.Li3_hasSum hx
theorem Y01.proof_3 : Y01.stmt_3 := fun _ _ h0 => Y01.hasDerivAt_Li3 h0
theorem Y01.proof_4 : Y01.stmt_4 := by
  refine ⟨?_, Y01.z3_hasSum, Y01.Li3_one⟩
  rw [Y01.z3_eq]; exact Y01.zeta3_eq
theorem Y01.proof_5 : Y01.stmt_5 := Y01.PA.main
theorem Y01.proof_6 : Y01.stmt_6 := Y01.PB.main Y01.proof_5
theorem Y01.proof_7 : Y01.stmt_7 := Y01.PC.main Y01.proof_6
theorem Y01.proof_8 : Y01.stmt_8 := Y01.PD.main
theorem Y01.proof_9 : Y01.stmt_9 :=
  ⟨fun _ hz => Y01.I3 hz, fun _ h0 h1 => Y01.Tr h0 h1, fun _ hw => Y01.D3 hw⟩
theorem Y01.proof_10 : Y01.stmt_10 := by
  intro x y hx0 hxy hy1
  have h := Y01.kummer hx0 hxy hy1
  simp only [Y01.KD, Y01.KA1, Y01.KA2, Y01.KA3, Y01.KA4, Y01.KA5, Y01.KA6, Y01.KA7] at h
  linarith
theorem Y01.proof_11 : Y01.stmt_11 := Y01.PE1.main.trans Y01.PE2.main
theorem Y01.proof_12 : Y01.stmt_12 := by
  unfold Y01.stmt_12
  rw [Y01.proof_7.2.2, Y01.proof_11, Y01.proof_8.1, sub_add_cancel]
theorem Y01.proof_13 : Y01.stmt_13 := Y01.PF.main Y01.proof_12
