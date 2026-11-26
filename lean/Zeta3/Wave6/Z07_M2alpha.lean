import Mathlib
import Zeta3.Wave2.D00_Dilog
import Zeta3.Wave2.D04_MalphaForms
/-!
# Z07 -- `M_2` at every angle (paper labels `thm:M2alpha`, `thm:flat`)

Source: `paper/ieee/main.tex`, Subsection "Smoothed Corners of Any Angle" (line 1465);
detailed proofs in `paper/ieee/proofs/sec_g8_M2alpha.tex`, `sec_g8_W3int.tex`, `sec_g6_flat.tex`.

**Setup, lines 1466-1471** (the corner integral; `D04.Mint` is this integral):

> its length is $2+2u\,M(u;\alpha)$, where
> \begin{multline*}
> M(u;\alpha)=\int_0^1\Bigl[(1+g)^{-1-u}\sqrt{1+g^{2-2u}-2\cos\alpha\,g^{1-u}}\\
> -1\Bigr]g^{u-1}\dd g=\sum_{k\ge0}M_k(\alpha)u^k,
> \end{multline*}

**`thm:flat`, lines 1501-1508** (the parts used here):

> For every $k\ge0$ and $0<\alpha<\pi$, $M_k(\alpha)=\sum_{n\ge1}\mu_{k,n}c^{2n}$, where each
> $\mu_{k,n}$ is a polynomial with rational coefficients in $\zeta(2),\dots,\zeta(k)$ of weight at
> most $k$. In particular, with $r_n,r_n',q_n'\in\mathbb Q$,
> \[ \mu_{2,n}=r_n+\frac{\zeta(2)}{4n^2-1},\qquad ... \]
> So, with $\Phi(c)=\sum_{n\ge1}\frac{c^{2n}}{4n^2-1}=\frac12-\frac{1-c^2}{2c}\operatorname{artanh}c$, ...

**Harmonic polylogarithms, lines 1516-1517**:

> for words $w$ in the letters $0,\pm1$ not ending in $0$, let $H_1(x)=-\log(1-x)$,
> $H_{-1}(x)=\log(1+x)$ and $H_{a,w}(x)=\int_0^xf_a(t)H_w(t)\dd t$, with $f_0=\frac1t$,
> $f_1=\frac1{1-t}$, $f_{-1}=\frac1{1+t}$; for example $H_{0,0,1}=\Li_3$.

**`thm:M2alpha`, lines 1518-1530** (verbatim):

> For $0<\alpha<\pi$, $M_2(\alpha)=W_2+W_3$, where $W_2=\frac12\zeta(2)-\frac12\Li_2(c^2)$ and
> $W_3$ is a combination, with coefficients rational in $c$, of iterated integrals of weight three
> whose letters are all real. Explicitly,
> \[ M_2(\alpha)=\zeta(2)\,\Phi(c)-\tfrac12\Li_2(c^2)+Y(c)+Y(-c), \]
> where
> \begin{multline*}
> Y(x)=(1-x)\bigl[H_{0,0,1}+H_{0,1,-1}-H_{1,-1,1}\bigr](x)\\
> +\frac{1-x^2}{2x}\,H_{1,0,1}(x)-\frac{(1-x)^2}{2x}\,H_{1,0,-1}(x).
> \end{multline*}

**After the theorem, lines 1536-1537**: "at the square corner $W_2$ is the term $T$ of
\eqref{eq:T}", where (line 903) $T=\frac{\pi^2+6\lambda^2}{24}=\frac{\zeta(2)-\Li_2(\frac12)}{2}$.

**Supplement, `sec_g8_W3int.tex`, Lemma `g8c:lem:M2int` (lines 99-112)**: with
$R(h)=\sqrt{1+h^2-2kh}$, $k=\cos\alpha$, $L=\log g$, $m=\log(1+g)$,
> f_A=\frac1g\Bigl[\frac{R-gR'+g^2R''}{2(1+g)}-\frac12\Bigr],\qquad
> f_B=\frac{gR'-R}{g(1+g)},\qquad f_C=\frac{R}{2g(1+g)},
> ... Then $f_AL^2+f_BLm+f_Cm^2$ is integrable on $(0,1]$ and
> M_2(\alpha)=\int_0^1\bigl[f_A(g)L^2+f_B(g)Lm+f_C(g)m^2\bigr]\dd g .

**Supplement, `sec_g8_M2alpha.tex`, eq. `g8:eq:M2rat` (line 100) and the output of M1
(lines 116-119)**:
> \sum_{n\ge1}r_nc^{2n}=-\tfrac12\Li_2(c^2)+Y(c)+Y(-c)\qquad(0<c<1),
> r_n=-\frac{(n-1)^2}{n^3(2n-1)}-\frac{n-1}{n^2(2n-1)}Z_{-1}(2n)-\frac{2}{(2n-1)(2n+1)}Z_{2}(2n)
> +\frac{Z_{-2}(2n)}{n(2n-1)(2n+1)}-\frac{Z_{-1,-1}(2n)}{n(2n-1)}
with the nested sums $Z_{a,v}(M)=\sum_{i=1}^M\operatorname{sgn}(a)^i\,i^{-|a|}Z_v(i-1)$,
$Z_\emptyset=1$ (proof of Lemma `g8:lem:shift`, line 59). For example $r_1=-\frac7{12}$.

## Formalization

* `c = cos(α/2)` is `D04.c`; `Li_2` is `D00.Li2`; the corner integral is `D04.Mint`.
* `ζ(2)` is written `π^2/6` (Mathlib `riemannZeta_two : riemannZeta 2 = π^2/6`), as in D04.
  No `ζ(3)` occurs in the statements (the theorem is `ζ(3)`-free, as the paper remarks).
* **`M_2(α)`.** The paper defines `M_k(α)` as the coefficient of `u^k` in the expansion of
  `M(u;α)`. We do not assume the power series; we use the Peano form of "coefficient of `u^2`":
  `IsCoeff2 α m` says `M(u;α) = a + b u + m u^2 + o(u^2)` as `u → 0` for some `a, b`. Such an `m`
  is unique (`stmt_6`), and when `M(·;α)` is analytic at `0` it is the Taylor coefficient. This
  matches D04, where `M_1(α)` is `∂_u M(0;α)` (the coefficient `b`), and the supplement, where
  `M_2 = ½ M''(0;α)`. `D04.Mint u α` is a genuine (integrable) integral for `|u| < 1` (D04
  docstring), so no junk value enters near `u = 0`.
* **Harmonic polylogarithms.** Letters are the inductive type `Letter` (`zero`, `one`, `negOne`
  for `0, 1, -1`); `H [] = 1` and `H (a :: w) x = ∫_0^x f_a(t) H_w(t) dt`. With this recursion
  `H [one] = ∫_0^x dt/(1-t)` and `H [negOne] = ∫_0^x dt/(1+t)`; `stmt_1` says these are the
  paper's `-log(1-x)` and `log(1+x)` on `(-1,1)`. `stmt_2` says every integrand that occurs in
  the eleven words used is interval integrable on `[0,x]` for `-1 < x < 1`, so no
  `intervalIntegral` junk value is used. (For `a = 0` the integrand `H_w(t)/t` takes the value
  `0` at the single point `t = 0` by Lean's `1/0 = 0`; this does not affect the integral.)
* `Y(x)` is evaluated at `x = ±c` with `0 < c < 1`, so `x ≠ 0` and no `x/0` convention is used.
* `Φ` is written in its closed form `½ - (1-c²)/(2c) artanh c` (`Real.artanh`, the genuine
  inverse of `tanh` on `(-1,1)`); `stmt_4` is the paper's identity with the series. In the final
  statement `Φ(c)` is this closed form.
* The integrand of the supplement's integral formula uses `R = √(1+g²-2kg) > 0` for `0 < α < π`,
  `R' = (g-k)/R`, `R'' = (1-k²)/R³`; `stmt_7` states its integrability on `(0,1]`.
* `n ≥ 1` is written as `n + 1` with `n : ℕ`.

Numerical check (mpmath, 40 digits; definition = `½ ∂_u² ` of the integrand under the integral,
HPLs by their Taylor series to order 900): the right side of `stmt_10` equals `M_2(α)` from the
definition, and from the integral of `stmt_7`, to `< 2·10^{-40}` at `α = π/3`
(`-0.04287175497685288253`), `π/2` (`-0.02216191108726939942`, the squircle's `M_2` of
`thm:M2`), `2π/3` (`-0.00956789565295638477`) and `1.8`; `stmt_9` (with the `r_n` above,
summed to `n = 400`) and `stmt_4` agree to `< 2·10^{-40}`; `r_1 = -7/12`; `stmt_5` holds exactly.

NOT COVERED:
* The first sentence of `thm:M2alpha`: "$W_3$ is a combination, with coefficients rational in
  $c$, of iterated integrals of weight three whose letters are all real". With `W_3 := M_2 - W_2`
  the decomposition `M_2 = W_2 + W_3` is a tautology; the content is the explicit formula
  (`stmt_10`) and the iterated-integral representation of `sec_g8_W3int.tex` (105 regularised
  iterated integrals), which is not formalized.
* The full expansion `M(u;α) = Σ_k M_k(α) u^k` (analyticity in `u`), and `thm:flat` for `k ≠ 2`
  (the statements for `M_3`, `Ψ`, and the weight claims for general `k`); only `k = 2` is used.
* The rationality `r_n ∈ ℚ` (evident from the explicit formula for `r_n`).
* The classical form with eight trilogarithms and the coefficient `7/4 + 2c` of `ζ(3)`
  (Corollary `g8:thm:classical` of the supplement; mentioned after the theorem in the paper).
* The decimal value of `M_2(π/2)` (checked numerically only, see above).
-/
-- BEGIN STATEMENTS
namespace Z07

open Real Filter Topology MeasureTheory

/-- The letters `0, 1, -1` of the harmonic polylogarithms. -/
inductive Letter
  | zero
  | one
  | negOne

/-- The kernels `f_0 = 1/t`, `f_1 = 1/(1-t)`, `f_{-1} = 1/(1+t)` (line 1517). -/
noncomputable def f : Letter → ℝ → ℝ
  | .zero, t => 1 / t
  | .one, t => 1 / (1 - t)
  | .negOne, t => 1 / (1 + t)

/-- Harmonic polylogarithm `H_w(x)`: `H_∅ = 1`, `H_{a,w}(x) = ∫_0^x f_a(t) H_w(t) dt` (line 1517;
the base cases `H_{±1}` agree with the paper's logarithms by `stmt_1`). -/
noncomputable def H : List Letter → ℝ → ℝ
  | [], _ => 1
  | a :: w, x => ∫ t in (0 : ℝ)..x, f a t * H w t

/-- The eleven words occurring in `Y` (weight three) and their tails. -/
def usedWords : List (List Letter) :=
  [ [.one], [.negOne],
    [.zero, .one], [.one, .negOne], [.negOne, .one], [.zero, .negOne],
    [.zero, .zero, .one], [.zero, .one, .negOne], [.one, .negOne, .one],
    [.one, .zero, .one], [.one, .zero, .negOne] ]

/-- `Y(x)` of `thm:M2alpha` (lines 1525-1528). -/
noncomputable def Y (x : ℝ) : ℝ :=
  (1 - x) * (H [.zero, .zero, .one] x + H [.zero, .one, .negOne] x - H [.one, .negOne, .one] x)
    + (1 - x ^ 2) / (2 * x) * H [.one, .zero, .one] x
    - (1 - x) ^ 2 / (2 * x) * H [.one, .zero, .negOne] x

/-- The closed form `Φ(c) = ½ - (1-c²)/(2c) artanh c` (line 1506). -/
noncomputable def PhiClosed (c : ℝ) : ℝ := 1 / 2 - (1 - c ^ 2) / (2 * c) * Real.artanh c

/-- "`m` is the coefficient of `u^2` of `M(u;α)`": `M(u;α) = a + b u + m u² + o(u²)` as `u → 0`. -/
def IsCoeff2 (α m : ℝ) : Prop :=
  ∃ a b : ℝ, (fun u : ℝ => D04.Mint u α - (a + b * u + m * u ^ 2)) =o[𝓝 0] (fun u => u ^ 2)

/-- The integrand `f_A L² + f_B L m + f_C m²` of Lemma `g8c:lem:M2int` (supplement). -/
noncomputable def M2integrand (α g : ℝ) : ℝ :=
  let k := Real.cos α
  let R := Real.sqrt (1 + g ^ 2 - 2 * k * g)
  let R1 := (g - k) / R
  let R2 := (1 - k ^ 2) / R ^ 3
  let L := Real.log g
  let m := Real.log (1 + g)
  let fA := ((R - g * R1 + g ^ 2 * R2) / (2 * (1 + g)) - 1 / 2) / g
  let fB := (g * R1 - R) / (g * (1 + g))
  let fC := R / (2 * g * (1 + g))
  fA * L ^ 2 + fB * L * m + fC * m ^ 2

/-- `Z_{-1}(M) = Σ_{i=1}^M (-1)^i / i`. -/
noncomputable def Zm1 (M : ℕ) : ℝ := ∑ i ∈ Finset.Icc 1 M, (-1 : ℝ) ^ i / (i : ℝ)

/-- `Z_{2}(M) = Σ_{i=1}^M 1 / i²`. -/
noncomputable def Z2 (M : ℕ) : ℝ := ∑ i ∈ Finset.Icc 1 M, 1 / (i : ℝ) ^ 2

/-- `Z_{-2}(M) = Σ_{i=1}^M (-1)^i / i²`. -/
noncomputable def Zm2 (M : ℕ) : ℝ := ∑ i ∈ Finset.Icc 1 M, (-1 : ℝ) ^ i / (i : ℝ) ^ 2

/-- `Z_{-1,-1}(M) = Σ_{i=1}^M (-1)^i / i · Z_{-1}(i-1)`. -/
noncomputable def Zm1m1 (M : ℕ) : ℝ :=
  ∑ i ∈ Finset.Icc 1 M, (-1 : ℝ) ^ i / (i : ℝ) * Zm1 (i - 1)

/-- The rational part `r_n` of `μ_{2,n}` (supplement, output of M1, lines 116-119); used for
`n ≥ 1`. -/
noncomputable def r (n : ℕ) : ℝ :=
  let N : ℝ := n;
  -(N - 1) ^ 2 / (N ^ 3 * (2 * N - 1)) - (N - 1) / (N ^ 2 * (2 * N - 1)) * Zm1 (2 * n)
    - 2 / ((2 * N - 1) * (2 * N + 1)) * Z2 (2 * n)
    + Zm2 (2 * n) / (N * (2 * N - 1) * (2 * N + 1)) - Zm1m1 (2 * n) / (N * (2 * N - 1))

/-- Line 1516: the base cases are the paper's `H_1(x) = -log(1-x)`, `H_{-1}(x) = log(1+x)`. -/
def stmt_1 : Prop :=
  ∀ x : ℝ, -1 < x → x < 1 →
    H [.one] x = -Real.log (1 - x) ∧ H [.negOne] x = Real.log (1 + x)

/-- Line 1517: the iterated integrals defining the harmonic polylogarithms used are genuine:
each integrand `f_a · H_w` is interval integrable on `[0, x]` for `-1 < x < 1`. -/
def stmt_2 : Prop :=
  ∀ (a : Letter) (w : List Letter), (a :: w) ∈ usedWords → ∀ x : ℝ, -1 < x → x < 1 →
    IntervalIntegrable (fun t => f a t * H w t) volume 0 x

/-- Line 1517: "for example `H_{0,0,1} = Li_3`", as the series `Σ_{n≥1} x^n/n³`, `|x| < 1`. -/
def stmt_3 : Prop :=
  ∀ x : ℝ, |x| < 1 →
    HasSum (fun n : ℕ => x ^ (n + 1) / ((n : ℝ) + 1) ^ 3) (H [.zero, .zero, .one] x)

/-- `thm:flat`, line 1506: `Φ(c) = Σ_{n≥1} c^{2n}/(4n²-1) = ½ - (1-c²)/(2c) artanh c`. -/
def stmt_4 : Prop :=
  ∀ c : ℝ, |c| < 1 → c ≠ 0 →
    HasSum (fun n : ℕ => c ^ (2 * (n + 1)) / (4 * ((n : ℝ) + 1) ^ 2 - 1)) (PhiClosed c)

/-- Lines 1536-1537 with eq. `eq:T`: at the square corner (`α = π/2`),
`W_2 = ½ζ(2) - ½Li_2(c²) = T = (π² + 6 log² 2)/24`. -/
def stmt_5 : Prop :=
  (π ^ 2 / 6) / 2 - D00.Li2 (D04.c (π / 2) ^ 2) / 2 = (π ^ 2 + 6 * Real.log 2 ^ 2) / 24

/-- The coefficient of `u²` is well defined: `IsCoeff2 α` has at most one solution. -/
def stmt_6 : Prop :=
  ∀ α m m' : ℝ, IsCoeff2 α m → IsCoeff2 α m' → m = m'

/-- Supplement, Lemma `g8c:lem:M2int`: the integrand is integrable on `(0,1]` and `M_2(α)`
(the coefficient of `u²` of `M(u;α)`) is `∫_0^1 [f_A L² + f_B L m + f_C m²] dg`. -/
def stmt_7 : Prop :=
  ∀ α : ℝ, 0 < α → α < π →
    IntegrableOn (M2integrand α) (Set.Ioc 0 1) ∧
      IsCoeff2 α (∫ g in (0 : ℝ)..1, M2integrand α g)

/-- `thm:flat` for `k = 2` (with `r_n` explicit, supplement line 117):
`M_2(α) = Σ_{n≥1} (r_n + ζ(2)/(4n²-1)) c^{2n}`, the series converging. -/
def stmt_8 : Prop :=
  ∀ α : ℝ, 0 < α → α < π →
    ∃ S : ℝ,
      HasSum (fun n : ℕ => (r (n + 1) + (π ^ 2 / 6) / (4 * ((n : ℝ) + 1) ^ 2 - 1))
        * D04.c α ^ (2 * (n + 1))) S ∧ IsCoeff2 α S

/-- Supplement, eq. `g8:eq:M2rat`: `Σ_{n≥1} r_n c^{2n} = -½ Li_2(c²) + Y(c) + Y(-c)` for
`0 < c < 1`. -/
def stmt_9 : Prop :=
  ∀ c : ℝ, 0 < c → c < 1 →
    HasSum (fun n : ℕ => r (n + 1) * c ^ (2 * (n + 1)))
      (-(1 / 2) * D00.Li2 (c ^ 2) + Y c + Y (-c))

/-- `thm:M2alpha`: for `0 < α < π`, `c = cos(α/2)`,
`M_2(α) = ζ(2) Φ(c) - ½ Li_2(c²) + Y(c) + Y(-c)`, with `M_2(α)` the coefficient of `u²` of
`M(u;α)`. -/
def stmt_10 : Prop :=
  ∀ α : ℝ, 0 < α → α < π →
    IsCoeff2 α ((π ^ 2 / 6) * PhiClosed (D04.c α) - (1 / 2) * D00.Li2 (D04.c α ^ 2)
      + Y (D04.c α) + Y (-(D04.c α)))

end Z07
-- END STATEMENTS

namespace Z07
open Real Filter Topology MeasureTheory Set

lemma H_one_eq {x : ℝ} (h0 : -1 < x) (h1 : x < 1) : H [.one] x = -Real.log (1 - x) := by
  simp only [H, f, mul_one]
  have : ∀ t ∈ uIcc (0:ℝ) x, HasDerivAt (fun t => -Real.log (1 - t)) (1 / (1 - t)) t := by
    intro t ht
    have ht1 : t < 1 := by
      rcases le_total 0 x with h | h
      · rw [uIcc_of_le h] at ht; linarith [ht.2]
      · rw [uIcc_of_ge h] at ht; linarith [ht.2]
    have h2 : HasDerivAt (fun y : ℝ => 1 - y) (-1) t := by
      simpa using (hasDerivAt_id t).const_sub 1
    have := (h2.log (by linarith)).neg
    convert this using 1; field_simp
  rw [intervalIntegral.integral_eq_sub_of_hasDerivAt this]
  · simp
  · apply ContinuousOn.intervalIntegrable
    apply ContinuousOn.div continuousOn_const (by fun_prop : ContinuousOn (fun t : ℝ => 1 - t) (uIcc 0 x))
    intro t ht
    rcases le_total 0 x with h | h
    · rw [uIcc_of_le h] at ht; have := ht.2; intro h'; linarith
    · rw [uIcc_of_ge h] at ht; have := ht.2; intro h'; linarith

lemma H_negOne_eq {x : ℝ} (h0 : -1 < x) (h1 : x < 1) : H [.negOne] x = Real.log (1 + x) := by
  simp only [H, f, mul_one]
  have : ∀ t ∈ uIcc (0:ℝ) x, HasDerivAt (fun t => Real.log (1 + t)) (1 / (1 + t)) t := by
    intro t ht
    have ht1 : -1 < t := by
      rcases le_total 0 x with h | h
      · rw [uIcc_of_le h] at ht; linarith [ht.1]
      · rw [uIcc_of_ge h] at ht; linarith [ht.1]
    have h2 : HasDerivAt (fun y : ℝ => 1 + y) 1 t := by
      simpa using (hasDerivAt_id t).const_add 1
    exact h2.log (by linarith)
  rw [intervalIntegral.integral_eq_sub_of_hasDerivAt this]
  · simp
  · apply ContinuousOn.intervalIntegrable
    apply ContinuousOn.div continuousOn_const (by fun_prop : ContinuousOn (fun t : ℝ => 1 + t) (uIcc 0 x))
    intro t ht
    rcases le_total 0 x with h | h
    · rw [uIcc_of_le h] at ht; have := ht.1; intro h'; linarith
    · rw [uIcc_of_ge h] at ht; have := ht.1; intro h'; linarith

/-- `H w` is continuous on `[-r, r]` and `O(t)` there, for every `r < 1`. -/
def Good (w : List Letter) : Prop :=
  ∀ r : ℝ, 0 ≤ r → r < 1 →
    ContinuousOn (H w) (Icc (-r) r) ∧ ∃ C : ℝ, ∀ t ∈ Icc (-r) r, |H w t| ≤ C * |t|

lemma measurable_f (a : Letter) : Measurable (f a) := by
  cases a
  · exact measurable_const.div measurable_id
  · exact measurable_const.div (measurable_const.sub measurable_id)
  · exact measurable_const.div (measurable_const.add measurable_id)

lemma abs_f_le {a : Letter} (ha : a ≠ .zero) {r t : ℝ} (hr : r < 1) (ht : t ∈ Icc (-r) r) :
    |f a t| ≤ 1 / (1 - r) := by
  have h1 : 0 < 1 - r := by linarith
  cases a
  · exact absurd rfl ha
  · simp only [f]
    have : 1 - r ≤ 1 - t := by linarith [ht.2]
    rw [abs_of_pos (by apply one_div_pos.mpr; linarith)]
    exact one_div_le_one_div_of_le h1 this
  · simp only [f]
    have : 1 - r ≤ 1 + t := by linarith [ht.1]
    rw [abs_of_pos (by apply one_div_pos.mpr; linarith)]
    exact one_div_le_one_div_of_le h1 this

/-- The step: from a bounded measurable integrand on `[-r, r]` to `Good (a :: w)`. -/
lemma good_of_bdd {a : Letter} {w : List Letter}
    (hb : ∀ r : ℝ, 0 ≤ r → r < 1 → AEStronglyMeasurable (fun t => f a t * H w t)
        (volume.restrict (Icc (-r) r)) ∧ ∃ D : ℝ, ∀ t ∈ Icc (-r) r, |f a t * H w t| ≤ D) :
    (∀ r : ℝ, 0 ≤ r → r < 1 → IntegrableOn (fun t => f a t * H w t) (Icc (-r) r)) ∧
      Good (a :: w) := by
  have hint : ∀ r : ℝ, 0 ≤ r → r < 1 → IntegrableOn (fun t => f a t * H w t) (Icc (-r) r) := by
    intro r hr0 hr1
    obtain ⟨hm, D, hD⟩ := hb r hr0 hr1
    have : IntegrableOn (fun t => f a t * H w t) univ (volume.restrict (Icc (-r) r)) := by
      apply Measure.integrableOn_of_bounded (M := D)
      · simp
      · exact hm
      · rw [Measure.restrict_univ]
        exact (ae_restrict_iff' measurableSet_Icc).mpr (Filter.Eventually.of_forall hD)
    simpa [IntegrableOn] using this
  refine ⟨hint, fun r hr0 hr1 => ⟨?_, ?_⟩⟩
  · have hii : IntervalIntegrable (fun t => f a t * H w t) volume (-r) r := by
      apply IntegrableOn.intervalIntegrable
      rw [uIcc_of_le (by linarith)]; exact hint r hr0 hr1
    have := intervalIntegral.continuousOn_primitive_interval' hii
      (a := 0) (by rw [uIcc_of_le (by linarith)]; constructor <;> linarith)
    rw [uIcc_of_le (by linarith)] at this
    exact this
  · obtain ⟨_, D, hD⟩ := hb r hr0 hr1
    refine ⟨D, fun t ht => ?_⟩
    simp only [H]
    have := intervalIntegral.norm_integral_le_of_norm_le_const (a := 0) (b := t)
      (f := fun t => f a t * H w t) (C := D) (fun s hs => by
        apply hD
        have hs' := uIoc_subset_uIcc hs
        rcases le_total 0 t with h | h
        · rw [uIcc_of_le h] at hs'; constructor <;> linarith [hs'.1, hs'.2, ht.1, ht.2]
        · rw [uIcc_of_ge h] at hs'; constructor <;> linarith [hs'.1, hs'.2, ht.1, ht.2])
    simpa using this

lemma bdd_nil {a : Letter} (ha : a ≠ .zero) (r : ℝ) (hr0 : 0 ≤ r) (hr1 : r < 1) :
    AEStronglyMeasurable (fun t => f a t * H [] t) (volume.restrict (Icc (-r) r)) ∧
      ∃ D : ℝ, ∀ t ∈ Icc (-r) r, |f a t * H [] t| ≤ D := by
  refine ⟨?_, 1 / (1 - r), fun t ht => ?_⟩
  · simp only [H, mul_one]; exact (measurable_f a).aestronglyMeasurable
  · simp only [H, mul_one]; exact abs_f_le ha hr1 ht

lemma bdd_cons (a : Letter) {w : List Letter} (hw : Good w) (r : ℝ) (hr0 : 0 ≤ r) (hr1 : r < 1) :
    AEStronglyMeasurable (fun t => f a t * H w t) (volume.restrict (Icc (-r) r)) ∧
      ∃ D : ℝ, ∀ t ∈ Icc (-r) r, |f a t * H w t| ≤ D := by
  obtain ⟨hc, C, hC⟩ := hw r hr0 hr1
  refine ⟨(measurable_f a).aestronglyMeasurable.mul (hc.aestronglyMeasurable measurableSet_Icc),
    ?_⟩
  have hC0 : ∀ t ∈ Icc (-r) r, |t| ≤ r := fun t ht => abs_le.mpr ⟨ht.1, ht.2⟩
  by_cases ha : a = .zero
  · subst ha
    refine ⟨|C|, fun t ht => ?_⟩
    simp only [f]
    rcases eq_or_ne t 0 with h | h
    · simp [h]
    · rw [abs_mul, abs_div, abs_one, div_mul_eq_mul_div, one_mul, div_le_iff₀ (abs_pos.mpr h)]
      exact (hC t ht).trans (mul_le_mul_of_nonneg_right (le_abs_self C) (abs_nonneg t))
  · refine ⟨1 / (1 - r) * (|C| * r), fun t ht => ?_⟩
    rw [abs_mul]
    apply mul_le_mul (abs_f_le ha hr1 ht) _ (abs_nonneg _) (le_of_lt (one_div_pos.mpr (by linarith)))
    exact (hC t ht).trans (mul_le_mul (le_abs_self C) (hC0 t ht) (abs_nonneg _) (abs_nonneg _))

lemma good_one : Good [.one] := (good_of_bdd (bdd_nil (by simp))).2
lemma good_negOne : Good [.negOne] := (good_of_bdd (bdd_nil (by simp))).2
lemma good_cons (a : Letter) {w : List Letter} (hw : Good w) : Good (a :: w) :=
  (good_of_bdd (bdd_cons a hw)).2

lemma ii_of_good (a : Letter) {w : List Letter} (hw : Good w ∨ (w = [] ∧ a ≠ .zero)) {x : ℝ}
    (h0 : -1 < x) (h1 : x < 1) : IntervalIntegrable (fun t => f a t * H w t) volume 0 x := by
  have hint : ∀ r : ℝ, 0 ≤ r → r < 1 → IntegrableOn (fun t => f a t * H w t) (Icc (-r) r) := by
    rcases hw with hw | ⟨rfl, ha⟩
    · exact (good_of_bdd (bdd_cons a hw)).1
    · exact (good_of_bdd (bdd_nil ha)).1
  apply IntegrableOn.intervalIntegrable
  apply (hint |x| (abs_nonneg x) (abs_lt.mpr ⟨h0, h1⟩)).mono_set
  intro t ht
  rcases le_total 0 x with h | h
  · rw [uIcc_of_le h] at ht; rw [abs_of_nonneg h]; constructor <;> linarith [ht.1, ht.2]
  · rw [uIcc_of_ge h] at ht; rw [abs_of_nonpos h]; constructor <;> linarith [ht.1, ht.2]

end Z07

theorem Z07.proof_1 : Z07.stmt_1 := fun _ h0 h1 => ⟨Z07.H_one_eq h0 h1, Z07.H_negOne_eq h0 h1⟩

theorem Z07.proof_2 : Z07.stmt_2 := by
  intro a w hmem x h0 h1
  open Z07 in
  have g1 := good_one; have gm := good_negOne
  simp only [usedWords, List.mem_cons, List.cons.injEq, List.not_mem_nil, or_false] at hmem
  rcases hmem with ⟨rfl, rfl⟩ | ⟨rfl, rfl⟩ | ⟨rfl, rfl⟩ | ⟨rfl, rfl⟩ | ⟨rfl, rfl⟩ | ⟨rfl, rfl⟩ |
      ⟨rfl, rfl⟩ | ⟨rfl, rfl⟩ | ⟨rfl, rfl⟩ | ⟨rfl, rfl⟩ | ⟨rfl, rfl⟩
  all_goals first
    | exact Z07.ii_of_good _ (Or.inr ⟨rfl, by simp⟩) h0 h1
    | exact Z07.ii_of_good _ (Or.inl (by assumption)) h0 h1
    | exact Z07.ii_of_good _ (Or.inl (Z07.good_cons _ (by assumption))) h0 h1

namespace Z07
open Real Filter Topology MeasureTheory Set

lemma ae_ne_zero : ∀ᵐ t : ℝ ∂volume, t ≠ 0 := by
  simpa [ae_iff] using measure_singleton (0 : ℝ)

lemma abs_le_of_mem_uIoc {x t : ℝ} (ht : t ∈ uIoc 0 x) : |t| ≤ |x| := by
  rcases le_total 0 x with h | h
  · rw [uIoc_of_le h] at ht
    rw [abs_of_nonneg h, abs_of_pos ht.1]; exact ht.2
  · rw [uIoc_of_ge h] at ht
    rw [abs_of_nonpos h, abs_of_nonpos ht.2]; linarith [ht.1]

/-- Termwise integration of `∫_0^x (1/t) G(t) dt` for a power series `G`. -/
lemma hasSum_H_zero {k : ℕ} {w : List Letter}
    (hG : ∀ t : ℝ, |t| < 1 → HasSum (fun n : ℕ => t ^ (n + 1) / ((n : ℝ) + 1) ^ k) (H w t))
    {x : ℝ} (hx : |x| < 1) :
    HasSum (fun n : ℕ => x ^ (n + 1) / ((n : ℝ) + 1) ^ (k + 1)) (H (.zero :: w) x) := by
  have key := intervalIntegral.hasSum_integral_of_dominated_convergence (μ := volume)
    (a := 0) (b := x) (F := fun (n : ℕ) (t : ℝ) => t ^ n / ((n : ℝ) + 1) ^ k)
    (f := fun t => f .zero t * H w t) (fun n _ => |x| ^ n)
    (fun n => by
      apply Continuous.aestronglyMeasurable
      fun_prop)
    (fun n => Eventually.of_forall fun t ht => by
      have h1 : (1 : ℝ) ≤ ((n : ℝ) + 1) ^ k := one_le_pow₀ (by linarith [n.cast_nonneg (α := ℝ)])
      rw [norm_div, norm_pow, Real.norm_eq_abs, Real.norm_eq_abs,
        abs_of_pos (by positivity : (0:ℝ) < ((n : ℝ) + 1) ^ k)]
      calc |t| ^ n / ((n : ℝ) + 1) ^ k ≤ |t| ^ n := div_le_self (by positivity) h1
        _ ≤ |x| ^ n := pow_le_pow_left₀ (abs_nonneg _) (abs_le_of_mem_uIoc ht) n)
    (Eventually.of_forall fun t _ => summable_geometric_of_lt_one (abs_nonneg _) hx)
    intervalIntegrable_const
    (by
      filter_upwards [ae_ne_zero] with t h0 ht
      have ht1 : |t| < 1 := (abs_le_of_mem_uIoc ht).trans_lt hx
      have := (hG t ht1).mul_left (1 / t)
      convert this using 1
      · funext n
        field_simp
        ring
      · rfl)
  have e : ∀ n : ℕ, (∫ t in (0:ℝ)..x, t ^ n / ((n : ℝ) + 1) ^ k)
      = x ^ (n + 1) / ((n : ℝ) + 1) ^ (k + 1) := by
    intro n
    rw [intervalIntegral.integral_div, integral_pow]
    have : ((n : ℝ) + 1) ≠ 0 := by positivity
    field_simp
    ring
  simp only [e] at key
  exact key

lemma hasSum_Li3 {x : ℝ} (hx : |x| < 1) :
    HasSum (fun n : ℕ => x ^ (n + 1) / ((n : ℝ) + 1) ^ 3) (H [.zero, .zero, .one] x) := by
  have h1 : ∀ t : ℝ, |t| < 1 →
      HasSum (fun n : ℕ => t ^ (n + 1) / ((n : ℝ) + 1) ^ 1) (H [.one] t) := by
    intro t ht
    rw [H_one_eq (abs_lt.mp ht).1 (abs_lt.mp ht).2]
    simpa using hasSum_pow_div_log_of_abs_lt_one ht
  exact hasSum_H_zero (fun t ht => hasSum_H_zero h1 ht) hx

lemma hasSum_artanh {c : ℝ} (hc : |c| < 1) :
    HasSum (fun k : ℕ => c ^ (2 * k + 1) / (2 * (k : ℝ) + 1)) (Real.artanh c) := by
  have h := (hasSum_log_sub_log_of_abs_lt_one hc).mul_left (1 / 2)
  have hc' := abs_lt.mp hc
  have e : Real.artanh c = 1 / 2 * (Real.log (1 + c) - Real.log (1 - c)) := by
    rw [Real.artanh_eq_half_log ⟨by linarith, by linarith⟩,
      Real.log_div (by linarith) (by linarith)]
  rw [e]
  convert h using 1
  funext k
  ring

lemma hasSum_Phi {c : ℝ} (hc : |c| < 1) (h0 : c ≠ 0) :
    HasSum (fun n : ℕ => c ^ (2 * (n + 1)) / (4 * ((n : ℝ) + 1) ^ 2 - 1)) (PhiClosed c) := by
  have hA := hasSum_artanh hc
  have hA1 := (hasSum_nat_add_iff' 1).mpr hA
  simp only [Finset.sum_range_one, Nat.cast_zero, mul_zero, zero_add, pow_one] at hA1
  have := ((hA.mul_left c).sub (hA1.mul_left (1 / c))).mul_left (1 / 2)
  convert this using 1
  · funext n
    have e1 : (4 * ((n : ℝ) + 1) ^ 2 - 1) = (2 * n + 1) * (2 * n + 3) := by ring
    have p1 : (2 * (n : ℝ) + 1) ≠ 0 := by positivity
    have p2 : (2 * (n : ℝ) + 3) ≠ 0 := by positivity
    rw [e1]
    push_cast
    have q1 : c ^ (2 * (n + 1)) = c ^ (2 * n + 1) * c := by ring
    have q2 : c ^ (2 * (n + 1) + 1) = c ^ (2 * n + 1) * c ^ 2 := by ring
    rw [q1, q2]
    field_simp
    ring
  · simp only [PhiClosed]
    field_simp
    ring

lemma stmt5_aux : (π ^ 2 / 6) / 2 - D00.Li2 (D04.c (π / 2) ^ 2) / 2
    = (π ^ 2 + 6 * Real.log 2 ^ 2) / 24 := by
  have hc : D04.c (π / 2) ^ 2 = 1 / 2 := by
    simp only [D04.c]
    rw [show π / 2 / 2 = π / 4 by ring, Real.cos_pi_div_four, div_pow, Real.sq_sqrt (by norm_num)]
    norm_num
  rw [hc, D00.Li2_half]
  ring

lemma coeff2_unique {α m m' : ℝ} (h : IsCoeff2 α m) (h' : IsCoeff2 α m') : m = m' := by
  obtain ⟨a, b, hab⟩ := h
  obtain ⟨a', b', hab'⟩ := h'
  have hd := hab'.sub hab
  by_contra hne
  have hC : 0 < |m - m'| := abs_pos.mpr (sub_ne_zero.mpr hne)
  obtain ⟨δ, hδ, hball⟩ := Metric.eventually_nhds_iff.mp (hd.def (half_pos hC))
  have key : ∀ u : ℝ, |u| < δ →
      |(a - a') + (b - b') * u + (m - m') * u ^ 2| ≤ |m - m'| / 2 * u ^ 2 := by
    intro u hu
    have := hball (by simpa [Real.dist_eq] using hu)
    simp only [Real.norm_eq_abs, abs_pow, sq_abs] at this
    convert this using 2
    ring
  have h0 := key 0 (by simpa using hδ)
  simp only [mul_zero, add_zero, ne_eq, OfNat.ofNat_ne_zero, not_false_eq_true, zero_pow,
    abs_nonpos_iff, sub_eq_zero] at h0
  have hu : |δ / 2| < δ := by rw [abs_of_pos (half_pos hδ)]; linarith
  have hu' : |-(δ / 2)| < δ := by rwa [abs_neg]
  have k1 := key _ hu
  have k2 := key _ hu'
  rw [h0, sub_self, zero_add] at k1 k2
  have : |2 * ((m - m') * (δ / 2) ^ 2)| ≤ 2 * (|m - m'| / 2 * (δ / 2) ^ 2) := by
    calc |2 * ((m - m') * (δ / 2) ^ 2)|
        = |((b - b') * (δ / 2) + (m - m') * (δ / 2) ^ 2)
            + ((b - b') * -(δ / 2) + (m - m') * (-(δ / 2)) ^ 2)| := by ring_nf
      _ ≤ _ := (abs_add_le _ _).trans (by linarith)
  rw [abs_mul, abs_mul, abs_two, abs_of_pos (by positivity : (0:ℝ) < (δ / 2) ^ 2)] at this
  have hδ2 : 0 < (δ / 2) ^ 2 := by positivity
  nlinarith

end Z07

theorem Z07.proof_3 : Z07.stmt_3 := fun _ hx => Z07.hasSum_Li3 hx

theorem Z07.proof_4 : Z07.stmt_4 := fun _ hc h0 => Z07.hasSum_Phi hc h0

theorem Z07.proof_5 : Z07.stmt_5 := Z07.stmt5_aux

theorem Z07.proof_6 : Z07.stmt_6 := fun _ _ _ h h' => Z07.coeff2_unique h h'

namespace Z07
open Real

/-- `thm:M2alpha` from `thm:flat` (`stmt_8`), eq. `g8:eq:M2rat` (`stmt_9`) and the closed form
of `Φ` (`stmt_4`). -/
lemma stmt10_of (h8 : stmt_8) (h9 : stmt_9) : stmt_10 := by
  intro α h0 h1
  obtain ⟨hc0, hc1, -⟩ := D04.cs_facts h0 h1
  obtain ⟨S, hS, hcoef⟩ := h8 α h0 h1
  have hc : |D04.c α| < 1 := abs_lt.mpr ⟨by linarith, hc1⟩
  have hsum := (h9 _ hc0 hc1).add ((hasSum_Phi hc hc0.ne').mul_left (π ^ 2 / 6))
  have : S = (π ^ 2 / 6) * PhiClosed (D04.c α) - (1 / 2) * D00.Li2 (D04.c α ^ 2)
      + Y (D04.c α) + Y (-(D04.c α)) := by
    refine hS.unique ?_
    convert hsum using 1
    · funext n; ring
    · ring
  rw [← this]; exact hcoef

end Z07

namespace Z07
open Real Filter Topology Asymptotics

noncomputable def Fdd (α u g : ℝ) : ℝ :=
  Real.log (1 + g) ^ 2 * (1 + g) ^ (-1 - u)
      * Real.sqrt (1 + (g ^ (1 - u)) ^ 2 - 2 * Real.cos α * g ^ (1 - u)) / g ^ (1 - u)
    + 2 * Real.log (1 + g) * Real.log g * (1 + g) ^ (-1 - u) * (g ^ (1 - u) - Real.cos α)
      / Real.sqrt (1 + (g ^ (1 - u)) ^ 2 - 2 * Real.cos α * g ^ (1 - u))
    - 2 * Real.log (1 + g) * Real.log g * (1 + g) ^ (-1 - u)
      * Real.sqrt (1 + (g ^ (1 - u)) ^ 2 - 2 * Real.cos α * g ^ (1 - u)) / g ^ (1 - u)
    + Real.log g ^ 2 * ((1 + g) ^ (-1 - u) * g ^ (1 - u)
        / Real.sqrt (1 + (g ^ (1 - u)) ^ 2 - 2 * Real.cos α * g ^ (1 - u))
      - (1 + g) ^ (-1 - u) * g ^ (1 - u) * (g ^ (1 - u) - Real.cos α) ^ 2
        / Real.sqrt (1 + (g ^ (1 - u)) ^ 2 - 2 * Real.cos α * g ^ (1 - u)) ^ 3
      - (1 + g) ^ (-1 - u) * (g ^ (1 - u) - Real.cos α)
        / Real.sqrt (1 + (g ^ (1 - u)) ^ 2 - 2 * Real.cos α * g ^ (1 - u))
      + ((1 + g) ^ (-1 - u)
        * Real.sqrt (1 + (g ^ (1 - u)) ^ 2 - 2 * Real.cos α * g ^ (1 - u)) - 1) / g ^ (1 - u))

lemma hasDerivAt_Fd {α u g : ℝ} (hα0 : 0 < α) (hα1 : α < π) (hg : 0 < g) :
    HasDerivAt (fun u => D04.Fd α u g) (Fdd α u g) u := by
  have hh : HasDerivAt (fun u : ℝ => g ^ (1 - u)) (Real.log g * (-1) * g ^ (1 - u)) u :=
    ((hasDerivAt_id u).const_sub 1).const_rpow hg |>.congr_deriv (by simp)
  have hp : HasDerivAt (fun u : ℝ => (1 + g) ^ (-1 - u))
      (Real.log (1 + g) * (-1) * (1 + g) ^ (-1 - u)) u :=
    ((hasDerivAt_id u).const_sub (-1)).const_rpow (by linarith) |>.congr_deriv (by simp)
  have hr := ((hh.pow 2).const_add 1).sub (hh.const_mul (2 * Real.cos α))
  have hpos := D04.rad_pos (x := g ^ (1 - u)) hα0 hα1
  have hsq := hr.sqrt (by exact hpos.ne')
  have hgu : 0 < g ^ (1 - u) := Real.rpow_pos_of_pos hg _
  have hsr : 0 < Real.sqrt (1 + (g ^ (1 - u)) ^ 2 - 2 * Real.cos α * g ^ (1 - u)) :=
    Real.sqrt_pos.mpr hpos
  have t1 := ((hp.mul hsq).div hh hgu.ne').const_mul (-Real.log (1 + g))
  have t2 := ((hp.mul (hh.sub_const (Real.cos α))).div hsq hsr.ne').const_mul (Real.log g)
  have t3 := (((hp.mul hsq).sub_const 1).div hh hgu.ne').const_mul (Real.log g)
  have := (t1.sub t2).add t3
  have hfun : (fun u => D04.Fd α u g) = fun u =>
      -Real.log (1 + g) * (((1 + g) ^ (-1 - u)
        * Real.sqrt (1 + (g ^ (1 - u)) ^ 2 - 2 * Real.cos α * g ^ (1 - u))) / g ^ (1 - u))
      - Real.log g * (((1 + g) ^ (-1 - u) * (g ^ (1 - u) - Real.cos α))
        / Real.sqrt (1 + (g ^ (1 - u)) ^ 2 - 2 * Real.cos α * g ^ (1 - u)))
      + Real.log g * (((1 + g) ^ (-1 - u)
        * Real.sqrt (1 + (g ^ (1 - u)) ^ 2 - 2 * Real.cos α * g ^ (1 - u)) - 1) / g ^ (1 - u)) := by
    funext u; unfold D04.Fd; ring
  rw [hfun]
  refine this.congr_deriv ?_
  unfold Fdd
  simp only [Pi.sub_apply, Pi.mul_apply, Pi.pow_apply, Nat.cast_ofNat]
  set S := Real.sqrt (1 + (g ^ (1 - u)) ^ 2 - 2 * Real.cos α * g ^ (1 - u)) with hS
  set h := g ^ (1 - u) with hhd
  set p := (1 + g) ^ (-1 - u) with hpd
  have hS0 : S ≠ 0 := hsr.ne'
  have hh0 : h ≠ 0 := hgu.ne'
  field_simp
  ring


lemma fdd_alg1 {m L p h S w κ g : ℝ} (hh0 : 0 < h) (hp0 : 0 < p) (hp1 : p ≤ 1)
    (hg0 : 0 < g) (hg1 : g ≤ 1) (hgw : g ≤ h * w) (hw1 : 1 ≤ w) (hlog : |L| ≤ 4 * w)
    (hL0 : 0 ≤ m) (hL1 : m ≤ g) (hS0 : 0 < S) (hS1 : S ≤ 1 + h) (hh1 : h ≤ 1)
    (hS3 : |h - κ| ≤ S) :
    |m ^ 2 * p * S / h + 2 * m * L * p * (h - κ) / S - 2 * m * L * p * S / h|
      ≤ 26 * (w * w) := by
  have hw0 : 0 ≤ w := by linarith
  have hww1 : w ≤ w * w := le_mul_of_one_le_left hw0 hw1
  have U1 : |m ^ 2 * p * S / h| ≤ 2 * w := by
    rw [abs_div, abs_of_pos hh0, div_le_iff₀ hh0,
      abs_of_nonneg (by positivity : 0 ≤ m ^ 2 * p * S)]
    have a1 : m ^ 2 ≤ g := by
      rw [sq]; calc m * m ≤ g * 1 := mul_le_mul hL1 (hL1.trans hg1) hL0 hg0.le
        _ = g := mul_one g
    have a2 : m ^ 2 * p * S ≤ g * 1 * 2 :=
      mul_le_mul (mul_le_mul a1 hp1 hp0.le hg0.le) (by linarith) hS0.le (by positivity)
    have : 2 * g ≤ 2 * w * h := by linarith
    linarith
  have U2 : |2 * m * L * p * (h - κ) / S| ≤ 8 * w := by
    rw [abs_div, abs_of_pos hS0, div_le_iff₀ hS0]
    have : |2 * m * L * p * (h - κ)| = 2 * m * |L| * p * |h - κ| := by
      rw [abs_mul, abs_mul, abs_mul, abs_mul, abs_of_pos hp0, abs_of_nonneg hL0]; norm_num
    rw [this]
    have b1 : 2 * m * |L| * p ≤ 2 * 1 * (4 * w) * 1 := by
      apply mul_le_mul _ hp1 hp0.le (by positivity)
      apply mul_le_mul _ hlog (abs_nonneg _) (by positivity)
      linarith
    have b2 := mul_le_mul b1 hS3 (abs_nonneg _) (by positivity)
    linarith
  have U3 : |2 * m * L * p * S / h| ≤ 16 * (w * w) := by
    rw [abs_div, abs_of_pos hh0, div_le_iff₀ hh0]
    have : |2 * m * L * p * S| = 2 * m * |L| * p * S := by
      rw [abs_mul, abs_mul, abs_mul, abs_mul, abs_of_pos hp0, abs_of_nonneg hL0,
        abs_of_pos hS0]; norm_num
    rw [this]
    have b1 : 2 * m * |L| * p * S ≤ 2 * (h * w) * (4 * w) * 1 * 2 := by
      apply mul_le_mul _ (by linarith) hS0.le (by positivity)
      apply mul_le_mul _ hp1 hp0.le (by positivity)
      apply mul_le_mul _ hlog (abs_nonneg _) (by positivity)
      linarith
    have : 2 * (h * w) * (4 * w) * 1 * 2 = 16 * (w * w) * h := by ring
    linarith
  have := abs_sub (m ^ 2 * p * S / h + 2 * m * L * p * (h - κ) / S) (2 * m * L * p * S / h)
  have := abs_add_le (m ^ 2 * p * S / h) (2 * m * L * p * (h - κ) / S)
  linarith

lemma fdd_alg2 {L p h S w κ g σ : ℝ} (hh0 : 0 < h) (hp0 : 0 < p) (hp1 : p ≤ 1)
    (hp2 : 1 - 2 * g ≤ p) (hg0 : 0 < g) (hgw : g ≤ h * w) (hw1 : 1 ≤ w)
    (hlog : |L| ≤ 4 * w) (hS0 : 0 < S) (hS1 : S ≤ 1 + h) (hS2 : 1 - h ≤ S) (hh1 : h ≤ 1)
    (hS3 : |h - κ| ≤ S) (hσ0 : 0 < σ) (hSσ : σ ≤ S) :
    |L ^ 2 * (p * h / S - p * h * (h - κ) ^ 2 / S ^ 3 - p * (h - κ) / S
      + (p * S - 1) / h)| ≤ 32 / σ * (w * w) + 64 * (w * w * w) := by
  have hw0 : 0 ≤ w := by linarith
  have hlog2 : L ^ 2 ≤ 16 * (w * w) := by
    have := mul_le_mul hlog hlog (abs_nonneg _) (by positivity)
    rw [← abs_mul, abs_of_nonneg (mul_self_nonneg _)] at this
    rw [sq]; linarith
  have hph : p * h ≤ 1 := by
    calc p * h ≤ 1 * 1 := mul_le_mul hp1 hh1 hh0.le zero_le_one
      _ = 1 := one_mul 1
  have V1 : |p * h / S| ≤ 1 / σ := by
    rw [abs_div, abs_of_pos hS0, abs_of_pos (by positivity), div_le_div_iff₀ hS0 hσ0]
    have := mul_le_mul hph hSσ hσ0.le zero_le_one
    linarith
  have V2 : |p * h * (h - κ) ^ 2 / S ^ 3| ≤ 1 / σ := by
    rw [abs_div, abs_of_nonneg (by positivity), abs_of_pos (by positivity),
      div_le_div_iff₀ (by positivity) hσ0]
    have c1 : (h - κ) ^ 2 ≤ S ^ 2 := by
      have := sq_le_sq' (abs_le.mp hS3).1 (abs_le.mp hS3).2; linarith
    have c3 : p * h * (h - κ) ^ 2 ≤ 1 * S ^ 2 :=
      mul_le_mul hph c1 (sq_nonneg _) zero_le_one
    have c4 : p * h * (h - κ) ^ 2 * σ ≤ (1 * S ^ 2) * S :=
      mul_le_mul c3 hSσ hσ0.le (by positivity)
    have : (1 * S ^ 2) * S = 1 * S ^ 3 := by ring
    linarith
  have V3 : |p * (h - κ) / S| ≤ 1 := by
    rw [abs_div, abs_of_pos hS0, div_le_iff₀ hS0, abs_mul, abs_of_pos hp0]
    have := mul_le_mul hp1 hS3 (abs_nonneg _) zero_le_one
    linarith
  have hpS : |p * S - 1| ≤ h + 2 * g := by
    rw [abs_le]; constructor
    · have b1 : p * (1 - h) ≤ p * S := mul_le_mul_of_nonneg_left hS2 hp0.le
      have b2 : (1 - 2 * g) * (1 - h) ≤ p * (1 - h) :=
        mul_le_mul_of_nonneg_right hp2 (by linarith)
      have b3 : (1 - 2 * g) * (1 - h) = 1 - h - 2 * g + 2 * (g * h) := by ring
      have b4 : 0 ≤ g * h := mul_nonneg hg0.le hh0.le
      linarith
    · have : p * S ≤ 1 * (1 + h) := mul_le_mul hp1 hS1 hS0.le (by norm_num)
      linarith
  have V4 : |(p * S - 1) / h| ≤ 1 + 2 * w := by
    rw [abs_div, abs_of_pos hh0, div_le_iff₀ hh0]
    have : (1 + 2 * w) * h = h + 2 * (h * w) := by ring
    linarith
  have hinner : |p * h / S - p * h * (h - κ) ^ 2 / S ^ 3 - p * (h - κ) / S + (p * S - 1) / h|
      ≤ 2 / σ + 4 * w := by
    have := abs_add_le (p * h / S - p * h * (h - κ) ^ 2 / S ^ 3 - p * (h - κ) / S)
      ((p * S - 1) / h)
    have := abs_sub (p * h / S - p * h * (h - κ) ^ 2 / S ^ 3) (p * (h - κ) / S)
    have := abs_sub (p * h / S) (p * h * (h - κ) ^ 2 / S ^ 3)
    have : 2 / σ = 1 / σ + 1 / σ := by ring
    linarith
  rw [abs_mul, abs_of_nonneg (sq_nonneg _)]
  have := mul_le_mul hlog2 hinner (abs_nonneg _) (by positivity)
  have e : 16 * (w * w) * (2 / σ + 4 * w) = 32 / σ * (w * w) + 64 * (w * w * w) := by ring
  linarith

lemma Fdd_bound {α u g : ℝ} (hα0 : 0 < α) (hα1 : α < π) (hg0 : 0 < g) (hg1 : g ≤ 1)
    (hu0 : -(1 / 4) < u) (hu1 : u < 1 / 4) :
    ‖Fdd α u g‖ ≤ (90 + 32 / Real.sin α) * g ^ (-(3 / 4) : ℝ) := by
  set κ := Real.cos α with hκ
  set σ := Real.sin α with hσ
  have hσ0 : 0 < σ := Real.sin_pos_of_pos_of_lt_pi hα0 hα1
  have hκσ : κ ^ 2 + σ ^ 2 = 1 := by rw [hκ, hσ, add_comm]; exact Real.sin_sq_add_cos_sq α
  set h := g ^ (1 - u) with hh
  set p := (1 + g) ^ (-1 - u) with hp
  set w := g ^ (-(1 / 4) : ℝ) with hw
  have hh0 : 0 < h := Real.rpow_pos_of_pos hg0 _
  have hh1 : h ≤ 1 := Real.rpow_le_one hg0.le hg1 (by linarith)
  have hw1 : 1 ≤ w := Real.one_le_rpow_of_pos_of_le_one_of_nonpos hg0 hg1 (by norm_num)
  have hgw : g ≤ h * w := by
    have e : g = h * g ^ u := by
      rw [hh, ← Real.rpow_add hg0]; norm_num
    have : g ^ u ≤ w := Real.rpow_le_rpow_of_exponent_ge hg0 hg1 hu0.le
    rw [e]; exact mul_le_mul_of_nonneg_left this hh0.le
  have hwww : w * w * w = g ^ (-(3 / 4) : ℝ) := by
    rw [hw, ← Real.rpow_add hg0, ← Real.rpow_add hg0]; norm_num
  have hlog : |Real.log g| ≤ 4 * w := by
    have := Real.abs_log_mul_self_rpow_lt g (1 / 4) hg0 hg1 (by norm_num)
    have e : g ^ (1 / 4 : ℝ) * w = 1 := by
      rw [hw, ← Real.rpow_add hg0]; norm_num
    have hq : 0 < g ^ (1 / 4 : ℝ) := Real.rpow_pos_of_pos hg0 _
    rw [abs_mul, abs_of_pos hq] at this
    nlinarith [abs_nonneg (Real.log g)]
  have hp0 : 0 < p := Real.rpow_pos_of_pos (by linarith) _
  have hp1 : p ≤ 1 := Real.rpow_le_one_of_one_le_of_nonpos (by linarith) (by linarith)
  have hp2 : 1 - 2 * g ≤ p := by
    have : (1 + g) ^ (-2 : ℝ) ≤ p := Real.rpow_le_rpow_of_exponent_le (by linarith) (by linarith)
    have e : (1 + g) ^ (-2 : ℝ) = 1 / (1 + g) ^ 2 := by
      rw [Real.rpow_neg (by linarith), one_div]; norm_num
    rw [e] at this
    have : 1 - 2 * g ≤ 1 / (1 + g) ^ 2 := by
      rw [le_div_iff₀ (by positivity)]; nlinarith [sq_nonneg g]
    linarith
  set r := 1 + h ^ 2 - 2 * κ * h with hr
  have hr0 : 0 < r := D04.rad_pos hα0 hα1
  set S := Real.sqrt r with hS
  have hS0 : 0 < S := Real.sqrt_pos.mpr hr0
  have hrσ : r = (h - κ) ^ 2 + σ ^ 2 := by rw [hr]; linear_combination (-1 : ℝ) * hκσ
  have hκ1 : -1 ≤ κ := Real.neg_one_le_cos α
  have hκ2 : κ ≤ 1 := Real.cos_le_one α
  have hS1 : S ≤ 1 + h := by
    rw [hS, show 1 + h = Real.sqrt ((1 + h) ^ 2) by rw [Real.sqrt_sq (by linarith)]]
    exact Real.sqrt_le_sqrt (by rw [hr]; nlinarith)
  have hS2 : 1 - h ≤ S := by
    rw [hS, show 1 - h = Real.sqrt ((1 - h) ^ 2) by rw [Real.sqrt_sq (by linarith)]]
    exact Real.sqrt_le_sqrt (by rw [hr]; nlinarith)
  have hS3 : |h - κ| ≤ S := Real.abs_le_sqrt (by rw [hrσ]; linarith [sq_nonneg σ])
  have hSσ : σ ≤ S := Real.le_sqrt_of_sq_le (by rw [hrσ]; linarith [sq_nonneg (h - κ)])
  have hL0 : 0 ≤ Real.log (1 + g) := Real.log_nonneg (by linarith)
  have hL1 : Real.log (1 + g) ≤ g := by
    have := Real.log_le_sub_one_of_pos (show 0 < 1 + g by linarith); linarith
  have A := fdd_alg1 (κ := κ) hh0 hp0 hp1 hg0 hg1 hgw hw1 hlog hL0 hL1 hS0 hS1 hh1 hS3
  have B := fdd_alg2 (κ := κ) hh0 hp0 hp1 hp2 hg0 hgw hw1 hlog hS0 hS1 hS2 hh1 hS3 hσ0 hSσ
  unfold Fdd
  rw [← hκ, ← hh, ← hp, ← hr, ← hS, Real.norm_eq_abs, ← hwww]
  have hw0 : 0 ≤ w := by linarith
  have hww1 : w ≤ w * w := le_mul_of_one_le_left hw0 hw1
  have hww2 : w * w ≤ w * w * w := le_mul_of_one_le_right (mul_nonneg hw0 hw0) hw1
  have k1 : 32 / σ * (w * w) ≤ 32 / σ * (w * w * w) :=
    mul_le_mul_of_nonneg_left hww2 (by positivity)
  have e : (90 + 32 / σ) * (w * w * w) = 90 * (w * w * w) + 32 / σ * (w * w * w) := by ring
  refine (abs_add_le _ _).trans ?_
  linarith


lemma Fdd_zero {α g : ℝ} (hα0 : 0 < α) (hα1 : α < π) (hg : 0 < g) :
    Fdd α 0 g = 2 * M2integrand α g := by
  unfold Fdd M2integrand
  simp only [sub_zero, Real.rpow_one, Real.rpow_neg_one]
  have hr0 := D04.rad_pos (x := g) hα0 hα1
  set κ := Real.cos α
  set S := Real.sqrt (1 + g ^ 2 - 2 * κ * g) with hS
  have hS0 : S ≠ 0 := (Real.sqrt_pos.mpr hr0).ne'
  have hSS : S ^ 2 = 1 + g ^ 2 - 2 * κ * g := Real.sq_sqrt hr0.le
  have e : (g - κ) ^ 2 = S ^ 2 - (1 - κ ^ 2) := by rw [hSS]; ring
  rw [e]
  have hg' : (1 + g) ≠ 0 := by linarith
  have hg0 : g ≠ 0 := hg.ne'
  field_simp
  ring

lemma taylor2 {f D : ℝ → ℝ} {d δ : ℝ} (hδ : 0 < δ)
    (hf : ∀ u, |u| < δ → HasDerivAt f (D u) u) (hD : HasDerivAt D d 0) :
    (fun u => f u - (f 0 + D 0 * u + d / 2 * u ^ 2)) =o[𝓝 0] (fun u => u ^ 2) := by
  have hD' := hD.isLittleO
  simp only [sub_zero, smul_eq_mul] at hD'
  rw [Asymptotics.isLittleO_iff] at hD' ⊢
  intro ε hε
  obtain ⟨η, hη, hball⟩ := Metric.eventually_nhds_iff.mp (hD' hε)
  rw [Metric.eventually_nhds_iff]
  refine ⟨min η δ, lt_min hη hδ, fun u hu => ?_⟩
  rw [Real.dist_eq, sub_zero] at hu
  set φ := fun v => f v - (f 0 + D 0 * v + d / 2 * v ^ 2) with hφ
  have hder : ∀ v ∈ Metric.closedBall (0:ℝ) |u|,
      HasDerivWithinAt φ (D v - D 0 - v * d) (Metric.closedBall (0:ℝ) |u|) v := by
    intro v hv
    rw [Metric.mem_closedBall, Real.dist_eq, sub_zero] at hv
    have hv' : |v| < δ := lt_of_le_of_lt hv (lt_of_lt_of_le hu (min_le_right _ _))
    have := ((hf v hv').sub (((hasDerivAt_id v).const_mul (D 0)).const_add (f 0) |>.add
      (((hasDerivAt_id v).pow 2).const_mul (d / 2))))
    refine (this.congr_deriv ?_).hasDerivWithinAt
    simp; ring
  have hbd : ∀ v ∈ Metric.closedBall (0:ℝ) |u|, ‖D v - D 0 - v * d‖ ≤ ε * |u| := by
    intro v hv
    rw [Metric.mem_closedBall, Real.dist_eq, sub_zero] at hv
    have := hball (show dist v 0 < η by
      rw [Real.dist_eq, sub_zero]; exact lt_of_le_of_lt hv (lt_of_lt_of_le hu (min_le_left _ _)))
    rw [Real.norm_eq_abs] at this ⊢
    refine this.trans ?_
    exact mul_le_mul_of_nonneg_left hv hε.le
  have := (convex_closedBall (0:ℝ) |u|).norm_image_sub_le_of_norm_hasDerivWithin_le hder hbd
    (Metric.mem_closedBall_self (abs_nonneg u)) (show u ∈ Metric.closedBall (0:ℝ) |u| by simp)
  have h0 : φ 0 = 0 := by simp [hφ]
  rw [h0, sub_zero, sub_zero, Real.norm_eq_abs] at this
  rw [Real.norm_eq_abs, Real.norm_eq_abs, abs_pow, sq_abs]
  calc |φ u| ≤ ε * |u| * |u| := this
    _ = ε * u ^ 2 := by rw [mul_assoc, ← sq, sq_abs]

lemma Fi_bound {α u g : ℝ} (hα0 : 0 < α) (hα1 : α < π) (hg0 : 0 < g) (hg1 : g ≤ 1)
    (hu : |u| < 1 / 4) : ‖D04.Fi α u g‖ ≤ 3 + 18 * g ^ (-(1 / 2) : ℝ) := by
  have hmv := (convex_ball (0:ℝ) (1 / 4)).norm_image_sub_le_of_norm_hasDerivWithin_le
    (f := fun u => D04.Fi α u g) (f' := fun u => D04.Fd α u g)
    (fun x _ => (D04.hasDerivAt_Fi hα0 hα1 hg0).hasDerivWithinAt)
    (fun x hx => by
      rw [Metric.mem_ball, Real.dist_eq, sub_zero, abs_lt] at hx
      exact D04.Fd_bound hα0 hα1 hg0 hg1 hx.1 hx.2)
    (Metric.mem_ball_self (by norm_num)) (show u ∈ Metric.ball (0:ℝ) (1 / 4) by
      rwa [Metric.mem_ball, Real.dist_eq, sub_zero])
  have h0 := D04.Fi0_bound hα0 hα1 hg0 hg1
  have hq : 0 ≤ 18 * g ^ (-(1 / 2) : ℝ) := by positivity
  simp only [sub_zero, Real.norm_eq_abs] at hmv h0 ⊢
  have : 18 * g ^ (-(1 / 2) : ℝ) * |u| ≤ 18 * g ^ (-(1 / 2) : ℝ) * 1 :=
    mul_le_mul_of_nonneg_left (by linarith) hq
  have := abs_sub_abs_le_abs_sub (D04.Fi α u g) (D04.Fi α 0 g)
  linarith

lemma meas_Fi (α u : ℝ) : Measurable (D04.Fi α u) := by unfold D04.Fi; fun_prop
lemma meas_Fd (α u : ℝ) : Measurable (D04.Fd α u) := by unfold D04.Fd; fun_prop
lemma meas_Fdd (α u : ℝ) : Measurable (Fdd α u) := by unfold Fdd; fun_prop

lemma ii_Fi {α u : ℝ} (hα0 : 0 < α) (hα1 : α < π) (hu : |u| < 1 / 4) :
    IntervalIntegrable (D04.Fi α u) MeasureTheory.volume 0 1 := by
  refine ((intervalIntegrable_const (c := (3 : ℝ))).add
    ((intervalIntegral.intervalIntegrable_rpow' (r := -(1 / 2)) (by norm_num)).const_mul 18)
    |>.mono_fun' (meas_Fi α u).aestronglyMeasurable ?_)
  refine (MeasureTheory.ae_restrict_iff' measurableSet_uIoc).mpr
    (Filter.Eventually.of_forall fun g hg => ?_)
  rw [Set.uIoc_of_le zero_le_one] at hg
  exact Fi_bound hα0 hα1 hg.1 hg.2 hu

lemma hasDerivAt_Mint_u {α u : ℝ} (hα0 : 0 < α) (hα1 : α < π) (hu : |u| < 1 / 8) :
    HasDerivAt (fun u => D04.Mint u α) (∫ g in (0 : ℝ)..1, D04.Fd α u g) u := by
  have hsub : Metric.ball u (1 / 8) ⊆ Metric.ball (0 : ℝ) (1 / 4) := by
    intro x hx
    rw [Metric.mem_ball, Real.dist_eq] at hx ⊢
    rw [sub_zero]
    calc |x| = |(x - u) + u| := by ring_nf
      _ ≤ |x - u| + |u| := abs_add_le _ _
      _ < 1 / 4 := by linarith
  have key := intervalIntegral.hasDerivAt_integral_of_dominated_loc_of_deriv_le
    (μ := MeasureTheory.volume) (a := 0) (b := 1)
    (F := fun u g => D04.Fi α u g) (F' := fun u g => D04.Fd α u g) (x₀ := u)
    (bound := fun g => 18 * g ^ (-(1 / 2) : ℝ)) (s := Metric.ball u (1 / 8))
    (Metric.ball_mem_nhds _ (by norm_num))
    (Filter.Eventually.of_forall (fun u => (meas_Fi α u).aestronglyMeasurable))
    (ii_Fi hα0 hα1 (by linarith)) (meas_Fd α u).aestronglyMeasurable
    (Filter.Eventually.of_forall fun g hg x hx => by
      rw [Set.uIoc_of_le zero_le_one] at hg
      have := hsub hx
      rw [Metric.mem_ball, Real.dist_eq, sub_zero, abs_lt] at this
      exact D04.Fd_bound hα0 hα1 hg.1 hg.2 this.1 this.2)
    ((intervalIntegral.intervalIntegrable_rpow' (by norm_num)).const_mul 18)
    (Filter.Eventually.of_forall fun g hg x _ => by
      rw [Set.uIoc_of_le zero_le_one] at hg
      exact D04.hasDerivAt_Fi hα0 hα1 hg.1)
  exact key.2

lemma hasDerivAt_D {α : ℝ} (hα0 : 0 < α) (hα1 : α < π) :
    IntervalIntegrable (Fdd α 0) MeasureTheory.volume 0 1 ∧
    HasDerivAt (fun u => ∫ g in (0 : ℝ)..1, D04.Fd α u g) (∫ g in (0 : ℝ)..1, Fdd α 0 g) 0 := by
  have key := intervalIntegral.hasDerivAt_integral_of_dominated_loc_of_deriv_le
    (μ := MeasureTheory.volume) (a := 0) (b := 1)
    (F := fun u g => D04.Fd α u g) (F' := fun u g => Fdd α u g) (x₀ := 0)
    (bound := fun g => (90 + 32 / Real.sin α) * g ^ (-(3 / 4) : ℝ))
    (s := Metric.ball (0 : ℝ) (1 / 4))
    (Metric.ball_mem_nhds _ (by norm_num))
    (Filter.Eventually.of_forall (fun u => (meas_Fd α u).aestronglyMeasurable))
    (D04.hasDerivAt_Mint hα0 hα1).1 (meas_Fdd α 0).aestronglyMeasurable
    (Filter.Eventually.of_forall fun g hg x hx => by
      rw [Set.uIoc_of_le zero_le_one] at hg
      rw [Metric.mem_ball, Real.dist_eq, sub_zero, abs_lt] at hx
      exact Fdd_bound hα0 hα1 hg.1 hg.2 hx.1 hx.2)
    ((intervalIntegral.intervalIntegrable_rpow' (by norm_num)).const_mul _)
    (Filter.Eventually.of_forall fun g hg x _ => by
      rw [Set.uIoc_of_le zero_le_one] at hg
      exact hasDerivAt_Fd hα0 hα1 hg.1)
  exact key

lemma stmt7_aux {α : ℝ} (hα0 : 0 < α) (hα1 : α < π) :
    MeasureTheory.IntegrableOn (M2integrand α) (Set.Ioc 0 1) ∧
      IsCoeff2 α (∫ g in (0 : ℝ)..1, M2integrand α g) := by
  obtain ⟨hii, hD⟩ := hasDerivAt_D hα0 hα1
  have hT := taylor2 (f := fun u => D04.Mint u α) (D := fun u => ∫ g in (0 : ℝ)..1, D04.Fd α u g)
    (δ := 1 / 8) (by norm_num) (fun u hu => hasDerivAt_Mint_u hα0 hα1 hu) hD
  have heq : Set.EqOn (fun g => Fdd α 0 g / 2) (M2integrand α) (Set.Ioc 0 1) := by
    intro g hg
    simp only
    rw [Fdd_zero hα0 hα1 hg.1]; ring
  have hint : ∫ g in (0 : ℝ)..1, M2integrand α g = (∫ g in (0 : ℝ)..1, Fdd α 0 g) / 2 := by
    rw [← intervalIntegral.integral_div]
    refine intervalIntegral.integral_congr_ae (Filter.Eventually.of_forall fun g hg => ?_)
    rw [Set.uIoc_of_le zero_le_one] at hg
    exact (heq hg).symm
  refine ⟨?_, ?_⟩
  · have := ((intervalIntegrable_iff_integrableOn_Ioc_of_le zero_le_one).mp hii).div_const 2
    exact MeasureTheory.IntegrableOn.congr_fun this heq measurableSet_Ioc
  · rw [hint]
    exact ⟨_, _, hT⟩
end Z07

theorem Z07.proof_7 : Z07.stmt_7 := fun _ h0 h1 => Z07.stmt7_aux h0 h1

namespace Z07
open Real Filter Topology MeasureTheory Set

/-- `H w t = Σ a m t^m` on `(-1,1)`, with `a 0 = 0` and `|a m| ≤ B`. -/
def Ser (w : List Letter) (a : ℕ → ℝ) (B : ℝ) : Prop :=
  (∀ t : ℝ, |t| < 1 → HasSum (fun m => a m * t ^ m) (H w t)) ∧ a 0 = 0 ∧ ∀ m, |a m| ≤ B

lemma summable_lin {c : ℕ → ℝ} {B : ℝ} (hc : ∀ m, |c m| ≤ B * ((m : ℝ) + 1)) {t : ℝ}
    (ht : |t| < 1) : Summable (fun m => c m * t ^ m) := by
  have h1 : Summable (fun m : ℕ => (m : ℝ) ^ 1 * |t| ^ m) :=
    summable_pow_mul_geometric_of_norm_lt_one 1 (by simpa using ht)
  have h0 : Summable (fun m : ℕ => |t| ^ m) := summable_geometric_of_lt_one (abs_nonneg _) ht
  refine Summable.of_norm_bounded ((h1.add h0).mul_left B) (fun m => ?_)
  rw [Real.norm_eq_abs, abs_mul, abs_pow]
  have := hc m
  have hp : 0 ≤ |t| ^ m := by positivity
  calc |c m| * |t| ^ m ≤ B * ((m : ℝ) + 1) * |t| ^ m := mul_le_mul_of_nonneg_right this hp
    _ = _ := by ring

/-- Termwise integration of `∫_0^x f_a(t) H_w(t) dt` from a power series of the integrand. -/
lemma hasSum_H_cons {a : Letter} {w : List Letter} {c : ℕ → ℝ} {B : ℝ}
    (hc : ∀ m, |c m| ≤ B * ((m : ℝ) + 1))
    (hG : ∀ t : ℝ, |t| < 1 → t ≠ 0 → HasSum (fun m => c m * t ^ m) (f a t * H w t))
    {x : ℝ} (hx : |x| < 1) :
    HasSum (fun m : ℕ => c m / ((m : ℝ) + 1) * x ^ (m + 1)) (H (a :: w) x) := by
  have hsx : Summable (fun m : ℕ => B * ((m : ℝ) + 1) * |x| ^ m) := by
    have h1 : Summable (fun m : ℕ => (m : ℝ) ^ 1 * |x| ^ m) :=
      summable_pow_mul_geometric_of_norm_lt_one 1 (by simpa using hx)
    have h0 : Summable (fun m : ℕ => |x| ^ m) := summable_geometric_of_lt_one (abs_nonneg _) hx
    convert (h1.add h0).mul_left B using 1
    funext m; ring
  have key := intervalIntegral.hasSum_integral_of_dominated_convergence (μ := volume)
    (a := 0) (b := x) (F := fun (n : ℕ) (t : ℝ) => c n * t ^ n)
    (f := fun t => f a t * H w t) (fun n _ => B * ((n : ℝ) + 1) * |x| ^ n)
    (fun n => (continuous_const.mul (continuous_pow n)).aestronglyMeasurable)
    (fun n => Eventually.of_forall fun t ht => by
      rw [Real.norm_eq_abs, abs_mul, abs_pow]
      exact mul_le_mul (hc n) (pow_le_pow_left₀ (abs_nonneg _) (abs_le_of_mem_uIoc ht) n)
        (by positivity) (le_trans (abs_nonneg _) (hc n)))
    (Eventually.of_forall fun t _ => hsx)
    intervalIntegrable_const
    (by
      filter_upwards [ae_ne_zero] with t h0 ht
      exact hG t ((abs_le_of_mem_uIoc ht).trans_lt hx) h0)
  have e : ∀ n : ℕ, (∫ t in (0:ℝ)..x, c n * t ^ n) = c n / ((n : ℝ) + 1) * x ^ (n + 1) := by
    intro n
    rw [intervalIntegral.integral_const_mul, integral_pow]
    have : ((n : ℝ) + 1) ≠ 0 := by positivity
    field_simp
    ring
  simp only [e] at key
  have hH : H (a :: w) x = ∫ t in (0:ℝ)..x, f a t * H w t := by rw [H]
  rw [hH]
  exact key

lemma hasSum_shift {a : ℕ → ℝ} {x S : ℝ} (h : HasSum (fun m => a (m + 1) * x ^ (m + 1)) S)
    (h0 : a 0 = 0) : HasSum (fun m => a m * x ^ m) S := by
  have := (hasSum_nat_add_iff' 1 (f := fun m => a m * x ^ m) (g := S)).mp
  simp only [Finset.sum_range_one, h0, zero_mul, sub_zero] at this
  exact this h

lemma hasSum_unshift {a : ℕ → ℝ} {x S : ℝ} (h : HasSum (fun m => a m * x ^ m) S)
    (h0 : a 0 = 0) : HasSum (fun m => a (m + 1) * x ^ (m + 1)) S := by
  have := (hasSum_nat_add_iff' 1 (f := fun m => a m * x ^ m) (g := S)).mpr h
  simpa [h0] using this

lemma ser_zero {w : List Letter} {p q : ℕ → ℝ} {B : ℝ} (hw : Ser w p B) (hq0 : q 0 = 0)
    (hq : ∀ m : ℕ, q (m + 1) = p (m + 1) / ((m : ℝ) + 1)) : Ser (.zero :: w) q B := by
  obtain ⟨hs, hp0, hb⟩ := hw
  have hB : 0 ≤ B := le_trans (abs_nonneg _) (hb 0)
  refine ⟨fun x hx => ?_, hq0, fun m => ?_⟩
  · apply hasSum_shift _ hq0
    have := hasSum_H_cons (a := .zero) (w := w) (c := fun m => p (m + 1)) (B := B)
      (fun m => le_trans (hb _) (le_mul_of_one_le_right hB (by linarith [m.cast_nonneg (α := ℝ)])))
      (fun t ht h0 => by
        have h1 := (hasSum_unshift (hs t ht) hp0).mul_left (1 / t)
        convert h1 using 1
        · funext m; rw [pow_succ]; field_simp
        · simp [f]) hx
    convert this using 1
    funext m; rw [hq]
  · rcases m with _ | m
    · rw [hq0, abs_zero]; exact hB
    · rw [hq, abs_div, abs_of_pos (by positivity : (0:ℝ) < (m : ℝ) + 1)]
      exact (div_le_self (abs_nonneg _) (by linarith [m.cast_nonneg (α := ℝ)])).trans (hb _)

lemma ser_pm {a : Letter} {w : List Letter} {p q : ℕ → ℝ} {B s : ℝ} (hw : Ser w p B)
    (hs1 : |s| ≤ 1) (hf : ∀ t, f a t = 1 / (1 - s * t)) (hq0 : q 0 = 0)
    (hq : ∀ m : ℕ, ((m : ℝ) + 1) * q (m + 1) = p m + s * m * q m) : Ser (a :: w) q B := by
  obtain ⟨hs, hp0, hb⟩ := hw
  have hB : 0 ≤ B := le_trans (abs_nonneg _) (hb 0)
  have hqb : ∀ m, |q m| ≤ B := by
    intro m
    induction m with
    | zero => rw [hq0, abs_zero]; exact hB
    | succ m ih =>
      have hm : (0:ℝ) < (m : ℝ) + 1 := by positivity
      have h1 : |((m : ℝ) + 1) * q (m + 1)| ≤ B * ((m : ℝ) + 1) := by
        rw [hq]
        calc |p m + s * m * q m| ≤ |p m| + |s| * m * |q m| := by
              refine (abs_add_le _ _).trans (le_of_eq ?_)
              rw [abs_mul, abs_mul, Nat.abs_cast]
          _ ≤ B + 1 * m * B := by
              gcongr
              · exact hb m
          _ = B * ((m : ℝ) + 1) := by ring
      rw [abs_mul, abs_of_pos hm] at h1
      nlinarith
  set c : ℕ → ℝ := fun m => ((m : ℝ) + 1) * q (m + 1) with hcdef
  have hcb : ∀ m, |c m| ≤ B * ((m : ℝ) + 1) := by
    intro m
    simp only [hcdef, abs_mul, abs_of_pos (by positivity : (0:ℝ) < (m : ℝ) + 1)]
    rw [mul_comm]; exact mul_le_mul_of_nonneg_right (hqb _) (by positivity)
  refine ⟨fun x hx => ?_, hq0, hqb⟩
  apply hasSum_shift _ hq0
  have := hasSum_H_cons (a := a) (w := w) (c := c) (B := B) hcb
    (fun t ht h0 => by
      have hst : |s * t| < 1 := by
        rw [abs_mul]; exact lt_of_le_of_lt (mul_le_of_le_one_left (abs_nonneg _) hs1) ht
      have hne : 1 - s * t ≠ 0 := by
        intro h; have := (abs_lt.mp hst).2; linarith
      obtain ⟨S, hS⟩ := summable_lin hcb ht
      have hS1 := (hasSum_nat_add_iff' 1 (f := fun m => c m * t ^ m) (g := S)).mpr hS
      have hP1 := hasSum_unshift (hs t ht) hp0
      have hsum := hP1.add (hS.mul_left (s * t))
      have e : (fun m => c (m + 1) * t ^ (m + 1))
          = (fun m => p (m + 1) * t ^ (m + 1) + s * t * (c m * t ^ m)) := by
        funext m
        simp only [hcdef]
        have := hq (m + 1)
        push_cast at this ⊢
        rw [pow_succ, this]; ring
      simp only [Finset.sum_range_one] at hS1
      rw [e] at hS1
      have hu := hS1.unique hsum
      have hc0 : c 0 = p 0 := by simp only [hcdef]; have := hq 0; simpa using this
      rw [hc0, hp0] at hu
      have : S = f a t * H w t := by
        rw [hf]; field_simp; linarith
      rwa [← this]) hx
  convert this using 1
  funext m
  simp only [hcdef]
  have : ((m : ℝ) + 1) ≠ 0 := by positivity
  push_cast
  field_simp

end Z07

namespace Z07
open Real Filter Topology MeasureTheory Set

lemma Zm1_succ (k : ℕ) : Zm1 (k + 1) = Zm1 k + (-1 : ℝ) ^ (k + 1) / ((k : ℝ) + 1) := by
  simp only [Zm1]; rw [Finset.sum_Icc_succ_top (by omega)]; push_cast; ring

lemma Z2_succ (k : ℕ) : Z2 (k + 1) = Z2 k + 1 / ((k : ℝ) + 1) ^ 2 := by
  simp only [Z2]; rw [Finset.sum_Icc_succ_top (by omega)]; push_cast; ring

lemma Zm2_succ (k : ℕ) : Zm2 (k + 1) = Zm2 k + (-1 : ℝ) ^ (k + 1) / ((k : ℝ) + 1) ^ 2 := by
  simp only [Zm2]; rw [Finset.sum_Icc_succ_top (by omega)]; push_cast; ring

lemma Zm1m1_succ (k : ℕ) :
    Zm1m1 (k + 1) = Zm1m1 k + (-1 : ℝ) ^ (k + 1) / ((k : ℝ) + 1) * Zm1 k := by
  simp only [Zm1m1]; rw [Finset.sum_Icc_succ_top (by omega)]; push_cast; simp

@[simp] lemma Zm1_zero : Zm1 0 = 0 := by simp [Zm1]
@[simp] lemma Z2_zero : Z2 0 = 0 := by simp [Z2]
@[simp] lemma Zm2_zero : Zm2 0 = 0 := by simp [Zm2]
@[simp] lemma Zm1m1_zero : Zm1m1 0 = 0 := by simp [Zm1m1]

lemma ser_one : Ser [.one] (fun m => 1 / (m : ℝ)) 1 := by
  refine ⟨fun t ht => ?_, by simp, fun m => ?_⟩
  · apply hasSum_shift _ (by simp)
    rw [H_one_eq (abs_lt.mp ht).1 (abs_lt.mp ht).2]
    convert hasSum_pow_div_log_of_abs_lt_one ht using 1
    funext m; push_cast; ring
  · rcases m with _ | m
    · simp
    · push_cast
      rw [abs_of_pos (by positivity)]
      exact div_le_one_of_le₀ (by linarith [m.cast_nonneg (α := ℝ)]) (by positivity)

lemma ser_negOne : Ser [.negOne] (fun m => (-1 : ℝ) ^ (m + 1) / (m : ℝ)) 1 := by
  refine ⟨fun t ht => ?_, by simp, fun m => ?_⟩
  · apply hasSum_shift _ (by simp)
    rw [H_negOne_eq (abs_lt.mp ht).1 (abs_lt.mp ht).2]
    have ht' : |-t| < 1 := by rwa [abs_neg]
    have := (hasSum_pow_div_log_of_abs_lt_one ht').neg
    rw [sub_neg_eq_add, neg_neg] at this
    convert this using 1
    funext m; push_cast; rw [neg_pow t]; ring
  · rw [abs_div, abs_pow, abs_neg, abs_one, one_pow, Nat.abs_cast]
    rcases m with _ | m
    · simp
    · push_cast
      exact div_le_one_of_le₀ (by linarith [m.cast_nonneg (α := ℝ)]) (by positivity)

lemma f_one (t : ℝ) : f .one t = 1 / (1 - 1 * t) := by simp [f]
lemma f_negOne (t : ℝ) : f .negOne t = 1 / (1 - (-1) * t) := by simp [f]

lemma ser_01 : Ser [.zero, .one] (fun m => 1 / (m : ℝ) ^ 2) 1 :=
  ser_zero ser_one (by simp) (fun m => by push_cast; field_simp)

lemma ser_001 : Ser [.zero, .zero, .one] (fun m => 1 / (m : ℝ) ^ 3) 1 :=
  ser_zero ser_01 (by simp) (fun m => by push_cast; field_simp)

lemma ser_0m1 : Ser [.zero, .negOne] (fun m => (-1 : ℝ) ^ (m + 1) / (m : ℝ) ^ 2) 1 :=
  ser_zero ser_negOne (by simp) (fun m => by push_cast; field_simp)

lemma ser_1m1 : Ser [.one, .negOne] (fun m => -Zm1 (m - 1) / (m : ℝ)) 1 :=
  ser_pm ser_negOne (by norm_num) f_one (by simp) (fun m => by
    rcases m with _ | m
    · simp
    · simp only [Nat.add_sub_cancel]
      rw [Zm1_succ]
      have : ((m : ℝ) + 1) ≠ 0 := by positivity
      push_cast
      field_simp
      ring)

lemma ser_01m1 : Ser [.zero, .one, .negOne] (fun m => -Zm1 (m - 1) / (m : ℝ) ^ 2) 1 :=
  ser_zero ser_1m1 (by simp) (fun m => by simp only [Nat.add_sub_cancel]; push_cast; field_simp)

lemma ser_m11 : Ser [.negOne, .one] (fun m => (-1 : ℝ) ^ (m + 1) * Zm1 (m - 1) / (m : ℝ)) 1 :=
  ser_pm ser_one (by norm_num) f_negOne (by simp) (fun m => by
    rcases m with _ | m
    · simp
    · simp only [Nat.add_sub_cancel]
      rw [Zm1_succ]
      have : ((m : ℝ) + 1) ≠ 0 := by positivity
      push_cast
      rcases neg_one_pow_eq_or ℝ m with h | h <;> simp only [pow_succ, h] <;> field_simp <;> ring)

lemma ser_1m11 : Ser [.one, .negOne, .one] (fun m => -Zm1m1 (m - 1) / (m : ℝ)) 1 :=
  ser_pm ser_m11 (by norm_num) f_one (by simp) (fun m => by
    rcases m with _ | m
    · simp
    · simp only [Nat.add_sub_cancel]
      rw [Zm1m1_succ]
      have : ((m : ℝ) + 1) ≠ 0 := by positivity
      push_cast
      field_simp
      ring)

lemma ser_101 : Ser [.one, .zero, .one] (fun m => Z2 (m - 1) / (m : ℝ)) 1 :=
  ser_pm ser_01 (by norm_num) f_one (by simp) (fun m => by
    rcases m with _ | m
    · simp
    · simp only [Nat.add_sub_cancel]
      rw [Z2_succ]
      have : ((m : ℝ) + 1) ≠ 0 := by positivity
      push_cast
      field_simp
      ring)

lemma ser_10m1 : Ser [.one, .zero, .negOne] (fun m => -Zm2 (m - 1) / (m : ℝ)) 1 :=
  ser_pm ser_0m1 (by norm_num) f_one (by simp) (fun m => by
    rcases m with _ | m
    · simp
    · simp only [Nat.add_sub_cancel]
      rw [Zm2_succ]
      have : ((m : ℝ) + 1) ≠ 0 := by positivity
      push_cast
      field_simp
      ring)

end Z07

namespace Z07
open Real Filter Topology MeasureTheory Set

/-- Coefficients of `H_{0,0,1} + H_{0,1,-1} - H_{1,-1,1}`, of `H_{1,0,1}` and of `H_{1,0,-1}`. -/
noncomputable def cA (m : ℕ) : ℝ := 1 / (m : ℝ) ^ 3 + -Zm1 (m - 1) / (m : ℝ) ^ 2 - -Zm1m1 (m - 1) / (m : ℝ)
noncomputable def cQ (m : ℕ) : ℝ := Z2 (m - 1) / (m : ℝ)
noncomputable def cR (m : ℕ) : ℝ := -Zm2 (m - 1) / (m : ℝ)

/-- Coefficients of `2x·Y(x)`. -/
noncomputable def cz (m : ℕ) : ℝ :=
  2 * (cA (m - 1) - cA (m - 1 - 1)) + (cQ m - cQ (m - 1 - 1))
    - (cR m - 2 * cR (m - 1) + cR (m - 1 - 1))

lemma hasSum_mulX {a : ℕ → ℝ} {x S : ℝ} (h : HasSum (fun m => a m * x ^ m) S) (h0 : a 0 = 0) :
    HasSum (fun m => a (m - 1) * x ^ m) (x * S) := by
  apply hasSum_shift _ (by simpa using h0)
  convert h.mul_left x using 1
  funext m; simp only [Nat.add_sub_cancel]; ring

lemma hasSum_Y {x : ℝ} (hx : |x| < 1) (hx0 : x ≠ 0) :
    HasSum (fun m => cz (m + 1) / 2 * x ^ m) (Y x) := by
  have hP : HasSum (fun m => cA m * x ^ m)
      (H [.zero, .zero, .one] x + H [.zero, .one, .negOne] x - H [.one, .negOne, .one] x) := by
    have := ((ser_001.1 x hx).add (ser_01m1.1 x hx)).sub (ser_1m11.1 x hx)
    convert this using 1
    funext m; simp only [cA]; ring
  have hQ : HasSum (fun m => cQ m * x ^ m) (H [.one, .zero, .one] x) := ser_101.1 x hx
  have hR : HasSum (fun m => cR m * x ^ m) (H [.one, .zero, .negOne] x) := ser_10m1.1 x hx
  have hA0 : cA 0 = 0 := by simp [cA]
  have hQ0 : cQ 0 = 0 := by simp [cQ]
  have hR0 : cR 0 = 0 := by simp [cR]
  have T1 := hasSum_mulX hP hA0
  have T2 := hasSum_mulX T1 (by simpa using hA0)
  have Q2 := hasSum_mulX (hasSum_mulX hQ hQ0) (by simpa using hQ0)
  have R1 := hasSum_mulX hR hR0
  have R2 := hasSum_mulX R1 (by simpa using hR0)
  have tot := (((T1.sub T2).mul_left 2).add (hQ.sub Q2)).sub ((hR.sub (R1.mul_left 2)).add R2)
  have h2 : HasSum (fun m => cz m * x ^ m) (2 * x * Y x) := by
    convert tot using 1
    · funext m; simp only [cz]; ring
    · simp only [Y]; field_simp; ring
  have hz0 : cz 0 = 0 := by simp [cz, hA0, hQ0, hR0]
  have h3 := (hasSum_unshift h2 hz0).mul_left (1 / (2 * x))
  convert h3 using 1
  · funext m; rw [pow_succ]; field_simp
  · field_simp

lemma cz_one : cz 1 = 0 := by simp [cz, cA, cQ, cR]

lemma coef_id (n : ℕ) : r (n + 1) = -1 / (2 * ((n : ℝ) + 1) ^ 2) + cz (2 * n + 3) := by
  have e1 : cz (2 * n + 3) = 2 * (cA (2 * n + 1 + 1) - cA (2 * n + 1)) + (cQ (2 * n + 2 + 1)
      - cQ (2 * n + 1)) - (cR (2 * n + 2 + 1) - 2 * cR (2 * n + 1 + 1) + cR (2 * n + 1)) := by
    simp only [cz]
    congr 3
  rw [e1]
  simp only [cA, cQ, cR, r, Nat.add_sub_cancel, show 2 * (n + 1) = 2 * n + 1 + 1 by ring]
  rw [Zm1m1_succ (2 * n + 1), Zm1m1_succ (2 * n), Zm1_succ (2 * n + 1), Zm1_succ (2 * n),
    Z2_succ (2 * n + 1), Z2_succ (2 * n), Zm2_succ (2 * n + 1), Zm2_succ (2 * n)]
  have o1 : (-1 : ℝ) ^ (2 * n + 1) = -1 := by rw [pow_succ, pow_mul]; simp
  have o2 : (-1 : ℝ) ^ (2 * n + 1 + 1) = 1 := by rw [pow_succ, o1]; simp
  rw [o1, o2]
  push_cast
  have p1 : (2 * (n : ℝ) + 1) ≠ 0 := by positivity
  have p2 : (2 * (n : ℝ) + 2) ≠ 0 := by positivity
  have p3 : (2 * (n : ℝ) + 3) ≠ 0 := by positivity
  have p4 : ((n : ℝ) + 1) ≠ 0 := by positivity
  rw [show (2 * ((n : ℝ) + 1) - 1) = 2 * n + 1 by ring, show (2 * ((n : ℝ) + 1) + 1) = 2 * n + 3 by ring,
    show (2 * (n : ℝ) + 1 + 1) = 2 * n + 2 by ring, show (2 * (n : ℝ) + 2 + 1) = 2 * n + 3 by ring,
    show ((n : ℝ) + 1 - 1) = n by ring]
  field_simp
  ring

lemma hasSum_r {c : ℝ} (h0 : 0 < c) (h1 : c < 1) :
    HasSum (fun n : ℕ => r (n + 1) * c ^ (2 * (n + 1)))
      (-(1 / 2) * D00.Li2 (c ^ 2) + Y c + Y (-c)) := by
  have hc : |c| < 1 := abs_lt.mpr ⟨by linarith, h1⟩
  have hc' : |-c| < 1 := by rwa [abs_neg]
  have hYY := (hasSum_Y hc h0.ne').add (hasSum_Y hc' (neg_ne_zero.mpr h0.ne'))
  have hinj : Function.Injective (fun n : ℕ => 2 * n + 2) := fun a b h => by simp at h; omega
  have hE := (hinj.hasSum_iff (f := fun m => cz (m + 1) / 2 * c ^ m + cz (m + 1) / 2 * (-c) ^ m)
    (a := Y c + Y (-c)) (fun m hm => by
      rcases Nat.even_or_odd m with ⟨k, hk⟩ | ho
      · rcases k with _ | k
        · simp only [hk]; simp [cz_one]
        · exact absurd ⟨k, by show 2 * k + 2 = m; omega⟩ hm
      · simp only [ho.neg_pow]; ring)).mpr hYY
  have hL := (D00.hasSum_Li2 (x := c ^ 2) (by rw [abs_of_pos (by positivity)]; nlinarith)).mul_left
    (-(1 / 2))
  convert hL.add hE using 1
  · funext n
    simp only [Function.comp]
    have hev : (-c) ^ (2 * n + 2) = c ^ (2 * n + 2) := by
      rw [show 2 * n + 2 = 2 * (n + 1) by ring, pow_mul, pow_mul, neg_sq]
    rw [show 2 * n + 2 + 1 = 2 * n + 3 by ring, hev, coef_id n, ← pow_mul]
    have p4 : ((n : ℝ) + 1) ≠ 0 := by positivity
    rw [show 2 * (n + 1) = 2 * n + 2 by ring]
    field_simp
    ring
  · ring

end Z07

namespace Z07
open Real Filter Topology MeasureTheory Set Interval

/-! ## Toward `stmt_8`: the coefficients `μ_{2,n}` (supplement `sec_g6_flat.tex`) -/

/-! ### Beta-log integrals with integer parameters -/

lemma bl_ae_Ioo01 : ∀ᵐ t ∂(volume.restrict (Ι (0:ℝ) 1)), t ∈ Ioo (0:ℝ) 1 := by
  rw [uIoc_of_le zero_le_one, ← Measure.restrict_congr_set Ioo_ae_eq_Ioc]
  exact ae_restrict_mem measurableSet_Ioo

lemma bl_ae_Ioo01' : ∀ᵐ t ∂(volume : Measure ℝ), t ∈ Ι (0:ℝ) 1 → t ∈ Ioo (0:ℝ) 1 :=
  (ae_restrict_iff' measurableSet_uIoc).mp bl_ae_Ioo01

lemma bl_ii_tt (c : ℝ) (hc : -1 < c) :
    IntervalIntegrable (fun t : ℝ => (t * (1 - t)) ^ c) volume 0 1 := by
  have h1 : IntervalIntegrable (fun t : ℝ => (t * (1 - t)) ^ c) volume 0 (1/2) := by
    have : IntervalIntegrable (fun t : ℝ => t ^ c * (1 - t) ^ c) volume 0 (1/2) := by
      apply IntervalIntegrable.mul_continuousOn (intervalIntegral.intervalIntegrable_rpow' hc)
      apply ContinuousOn.rpow_const (by fun_prop)
      intro x hx; left; rw [uIcc_of_le (by norm_num)] at hx; linarith [hx.2]
    refine this.congr ?_
    intro t ht
    rw [uIoc_of_le (by norm_num)] at ht
    exact (mul_rpow ht.1.le (by linarith [ht.2])).symm
  have h2 := (h1.comp_sub_left 1).symm
  refine h1.trans ?_
  convert h2 using 1
  · funext t; ring_nf
  · norm_num
  · norm_num

lemma bl_abs_log_le {x θ : ℝ} (h0 : 0 < x) (h1 : x ≤ 1) (hθ : 0 < θ) :
    |log x| ≤ x ^ (-θ) / θ := by
  have := abs_log_mul_self_rpow_lt x θ h0 h1 hθ
  have hp := rpow_pos_of_pos h0 θ
  rw [abs_mul, abs_of_pos hp] at this
  rw [rpow_neg h0.le, show (x ^ θ)⁻¹ / θ = (1 / θ) / x ^ θ by field_simp, le_div_iff₀ hp]
  linarith

lemma bl_logsum_le {t θ : ℝ} (ht : t ∈ Ioo (0:ℝ) 1) (hθ : 0 < θ) :
    |log t| + |log (1 - t)| ≤ 2 / θ * (t * (1 - t)) ^ (-θ) := by
  obtain ⟨h0, h1⟩ := ht
  have hu : 0 < 1 - t := by linarith
  have hp : 0 < t * (1 - t) := mul_pos h0 hu
  have e1 := bl_abs_log_le h0 h1.le hθ
  have e2 := bl_abs_log_le hu (by linarith) hθ
  have f1 : t ^ (-θ) ≤ (t * (1 - t)) ^ (-θ) :=
    rpow_le_rpow_of_nonpos hp (by nlinarith) (by linarith)
  have f2 : (1 - t) ^ (-θ) ≤ (t * (1 - t)) ^ (-θ) :=
    rpow_le_rpow_of_nonpos hp (by nlinarith) (by linarith)
  have : t ^ (-θ) / θ + (1 - t) ^ (-θ) / θ ≤ 2 / θ * (t * (1 - t)) ^ (-θ) := by
    rw [← add_div, div_le_iff₀ hθ]
    have : 2 / θ * (t * (1 - t)) ^ (-θ) * θ = 2 * (t * (1 - t)) ^ (-θ) := by field_simp
    rw [this]; linarith
  linarith

/-- the log-sum `|log t| + |log(1-t)|` -/
noncomputable def bls (t : ℝ) : ℝ := |log t| + |log (1 - t)|

lemma bls_nonneg (t : ℝ) : 0 ≤ bls t := by unfold bls; positivity

lemma bls_pow_le {t : ℝ} (ht : t ∈ Ioo (0:ℝ) 1) (m : ℕ) :
    bls t ^ m ≤ (4 * (m + 1)) ^ m * (t * (1 - t)) ^ (-(1/2:ℝ)) := by
  have hp : 0 < t * (1 - t) := mul_pos ht.1 (by linarith [ht.2])
  have hθ : (0:ℝ) < 1 / (2 * (m + 1)) := by positivity
  have h := bl_logsum_le ht hθ
  have h' : bls t ^ m ≤ (2 / (1 / (2 * (m + 1))) * (t * (1 - t)) ^ (-(1 / (2 * ((m:ℝ) + 1))))) ^ m :=
    pow_le_pow_left₀ (bls_nonneg t) h m
  rw [mul_pow, ← rpow_mul_natCast hp.le] at h'
  have hle1 : t * (1 - t) ≤ 1 := by nlinarith [ht.1, ht.2]
  have : (t * (1 - t)) ^ (-(1 / (2 * ((m:ℝ) + 1))) * (m:ℝ)) ≤ (t * (1 - t)) ^ (-(1/2:ℝ)) := by
    apply rpow_le_rpow_of_exponent_ge hp hle1
    rw [neg_mul, neg_le_neg_iff, div_mul_eq_mul_div, one_mul, div_le_iff₀ (by positivity)]
    linarith
  have e : (2 / (1 / (2 * ((m:ℝ) + 1)))) = 4 * (m + 1) := by field_simp; ring
  rw [e] at h'
  calc bls t ^ m ≤ _ := h'
    _ ≤ _ := by gcongr

lemma bl_ii_of_ls {f : ℝ → ℝ} (hf : AEStronglyMeasurable f (volume.restrict (Ι (0:ℝ) 1)))
    (K : ℝ) (m : ℕ) (h : ∀ t ∈ Ioo (0:ℝ) 1, |f t| ≤ K * bls t ^ m) :
    IntervalIntegrable f volume 0 1 := by
  refine IntervalIntegrable.mono_fun'
    (((bl_ii_tt (-(1/2)) (by norm_num)).const_mul (|K| * (4 * (m + 1)) ^ m))) hf ?_
  filter_upwards [bl_ae_Ioo01] with t ht
  have h1 := h t ht
  have h2 := bls_pow_le ht m
  have hK : K ≤ |K| := le_abs_self K
  have hl := pow_nonneg (bls_nonneg t) m
  simp only [Real.norm_eq_abs]
  calc |f t| ≤ K * bls t ^ m := h1
    _ ≤ |K| * bls t ^ m := by gcongr
    _ ≤ |K| * ((4 * (m + 1)) ^ m * (t * (1 - t)) ^ (-(1/2:ℝ))) := by gcongr
    _ = _ := by ring

/-- `Kb i m P Q = ∫_0^1 y^P (1-y)^Q log^i y log^m (1-y) dy`. -/
noncomputable def Kb (i m P Q : ℕ) : ℝ :=
  ∫ y in (0:ℝ)..1, y ^ P * (1 - y) ^ Q * (log y ^ i * log (1 - y) ^ m)

lemma ii_Kb (i m P Q : ℕ) : IntervalIntegrable
    (fun y : ℝ => y ^ P * (1 - y) ^ Q * (log y ^ i * log (1 - y) ^ m)) volume 0 1 := by
  refine bl_ii_of_ls (by fun_prop : Measurable (fun y : ℝ => y ^ P * (1 - y) ^ Q *
    (log y ^ i * log (1 - y) ^ m))).aestronglyMeasurable 1 (i + m) fun t ht => ?_
  have h0 : 0 < t := ht.1
  have h1 : 0 < 1 - t := by linarith [ht.2]
  rw [abs_mul, abs_mul, abs_mul, abs_pow, abs_pow, abs_pow, abs_pow, one_mul, pow_add,
    abs_of_pos h0, abs_of_pos h1]
  have a1 : |log t| ≤ bls t := by unfold bls; linarith [abs_nonneg (log (1 - t))]
  have a2 : |log (1 - t)| ≤ bls t := by unfold bls; linarith [abs_nonneg (log t)]
  have hb := bls_nonneg t
  have p1 : t ^ P ≤ 1 := pow_le_one₀ h0.le (by linarith [ht.2])
  have p2 : (1 - t) ^ Q ≤ 1 := pow_le_one₀ h1.le (by linarith)
  have q1 : |log t| ^ i * |log (1 - t)| ^ m ≤ bls t ^ i * bls t ^ m := by gcongr
  have q0 : 0 ≤ |log t| ^ i * |log (1 - t)| ^ m := by positivity
  have : t ^ P * (1 - t) ^ Q ≤ 1 := by
    calc t ^ P * (1 - t) ^ Q ≤ 1 * 1 := by gcongr
      _ = 1 := by ring
  calc t ^ P * (1 - t) ^ Q * (|log t| ^ i * |log (1 - t)| ^ m)
      ≤ 1 * (|log t| ^ i * |log (1 - t)| ^ m) := by gcongr
    _ ≤ _ := by linarith

lemma Kb_succ (i m P Q : ℕ) : Kb i m P (Q + 1) = Kb i m P Q - Kb i m (P + 1) Q := by
  unfold Kb
  rw [← intervalIntegral.integral_sub (ii_Kb i m P Q) (ii_Kb i m (P + 1) Q)]
  congr 1; funext y; ring

lemma Kb_symm (i m P Q : ℕ) : Kb i m P Q = Kb m i Q P := by
  unfold Kb
  have e := intervalIntegral.integral_comp_sub_left
    (fun y : ℝ => y ^ Q * (1 - y) ^ P * (log y ^ m * log (1 - y) ^ i)) (a := 0) (b := 1) 1
  norm_num at e
  rw [← e]
  congr 1; funext y; ring

/-- the Beta values `Bt P Q = P! Q! / (P+Q+1)!`. -/
lemma Kb00 (P Q : ℕ) : Kb 0 0 P Q = (P.factorial : ℝ) * Q.factorial / (P + Q + 1).factorial := by
  induction Q generalizing P with
  | zero =>
    simp only [Kb, pow_zero, mul_one, integral_pow]
    simp [Nat.factorial_succ]
    field_simp
  | succ Q ih =>
    rw [Kb_succ, ih, ih]
    rw [show P + 1 + Q + 1 = (P + Q + 1) + 1 by ring, show P + (Q + 1) + 1 = (P + Q + 1) + 1 by ring]
    simp only [Nat.factorial_succ]
    push_cast
    have : ((P + Q + 1).factorial : ℝ) ≠ 0 := by positivity
    field_simp
    ring

lemma Kb00_pos (P Q : ℕ) : 0 < Kb 0 0 P Q := by rw [Kb00]; positivity

lemma Kb00_succR (P Q : ℕ) : Kb 0 0 P (Q + 1) = Kb 0 0 P Q * ((Q : ℝ) + 1) / ((P : ℝ) + Q + 2) := by
  rw [Kb00, Kb00, show P + (Q + 1) + 1 = (P + Q + 1) + 1 by ring]
  simp only [Nat.factorial_succ]
  push_cast
  have : ((P + Q + 1).factorial : ℝ) ≠ 0 := by positivity
  have : (P : ℝ) + Q + 2 ≠ 0 := by positivity
  field_simp
  ring

lemma Kb00_succL (P Q : ℕ) : Kb 0 0 (P + 1) Q = Kb 0 0 P Q * ((P : ℝ) + 1) / ((P : ℝ) + Q + 2) := by
  rw [Kb00, Kb00, show P + 1 + Q + 1 = (P + Q + 1) + 1 by ring]
  simp only [Nat.factorial_succ]
  push_cast
  have : ((P + Q + 1).factorial : ℝ) ≠ 0 := by positivity
  have : (P : ℝ) + Q + 2 ≠ 0 := by positivity
  field_simp
  ring

/-- harmonic numbers -/
noncomputable def hh (N : ℕ) : ℝ := ∑ i ∈ Finset.range N, 1 / ((i : ℝ) + 1)
noncomputable def hh2 (N : ℕ) : ℝ := ∑ i ∈ Finset.range N, 1 / ((i : ℝ) + 1) ^ 2

lemma hh_succ (N : ℕ) : hh (N + 1) = hh N + 1 / ((N : ℝ) + 1) := by
  simp [hh, Finset.sum_range_succ]
lemma hh2_succ (N : ℕ) : hh2 (N + 1) = hh2 N + 1 / ((N : ℝ) + 1) ^ 2 := by
  simp [hh2, Finset.sum_range_succ]

lemma bl_uIoo01 : uIoo (0:ℝ) 1 = Ioo 0 1 := by simp [uIoo]

lemma bl_t_log_pow_tendsto (k : ℕ) (hk : k ≠ 0) :
    Tendsto (fun t : ℝ => t * log t ^ k) (𝓝[>] 0) (𝓝 0) := by
  have hr : (0:ℝ) < (k:ℝ)⁻¹ := by positivity
  have h := (tendsto_log_mul_rpow_nhdsGT_zero hr).pow k
  rw [zero_pow hk] at h
  refine h.congr' ?_
  filter_upwards [self_mem_nhdsWithin] with t ht
  rw [mul_pow, rpow_inv_natCast_pow (le_of_lt ht) hk, mul_comm]

lemma bl_G_cont (n k : ℕ) (hk : k ≠ 0) :
    ContinuousOn (fun t : ℝ => t ^ (n+1) * log t ^ k) (Icc 0 1) := by
  intro x hx
  rcases eq_or_lt_of_le hx.1 with h | h
  · subst h
    apply ContinuousWithinAt.mono _ Icc_subset_Ici_self
    rw [← continuousWithinAt_Ioi_iff_Ici, ContinuousWithinAt]
    have := (((continuous_pow n).tendsto 0).mono_left nhdsWithin_le_nhds).mul
      (bl_t_log_pow_tendsto k hk)
    rw [mul_zero] at this
    simp only [ne_eq, Nat.add_eq_zero_iff, one_ne_zero, and_false, not_false_eq_true, zero_pow,
      zero_mul]
    refine this.congr fun t => ?_
    ring
  · exact ((continuous_pow _).continuousAt.mul ((continuousAt_log h.ne').pow k)).continuousWithinAt

/-- `bJ n k = ∫_0^1 t^n log^k t` -/
noncomputable def bJ (n k : ℕ) : ℝ := ∫ t in (0:ℝ)..1, t ^ n * log t ^ k

lemma ii_bJ (n k : ℕ) : IntervalIntegrable (fun t : ℝ => t ^ n * log t ^ k) volume 0 1 := by
  have := ii_Kb k 0 n 0
  simpa using this

lemma bJ_rec (n k : ℕ) : (n + 1 : ℝ) * bJ n (k+1) + (k + 1) * bJ n k = 0 := by
  have hderiv : ∀ x ∈ Ioo (0:ℝ) 1, HasDerivAt (fun t : ℝ => t ^ (n+1) * log t ^ (k+1))
      ((n + 1 : ℝ) * (x ^ n * log x ^ (k+1)) + (k + 1) * (x ^ n * log x ^ k)) x := by
    intro x hx
    have hx0 : x ≠ 0 := hx.1.ne'
    have h2 : HasDerivAt (fun t => log t ^ (k+1)) (((k+1:ℕ):ℝ) * log x ^ k * x⁻¹) x := by
      have := (hasDerivAt_pow (k+1) (log x)).comp x (hasDerivAt_log hx0)
      simpa [Function.comp_def] using this
    have h1 := (hasDerivAt_pow (n+1) x).mul h2
    convert h1 using 1
    simp only [Nat.add_sub_cancel]
    field_simp
    push_cast; ring
  have hint : IntervalIntegrable (fun x : ℝ => (n + 1 : ℝ) * (x ^ n * log x ^ (k+1)) +
      (k + 1) * (x ^ n * log x ^ k)) volume 0 1 :=
    ((ii_bJ n (k+1)).const_mul _).add ((ii_bJ n k).const_mul _)
  have := intervalIntegral.integral_eq_sub_of_hasDerivAt_of_le zero_le_one
    (bl_G_cont n (k+1) (by omega)) hderiv hint
  rw [intervalIntegral.integral_add ((ii_bJ n (k+1)).const_mul _) ((ii_bJ n k).const_mul _),
    intervalIntegral.integral_const_mul, intervalIntegral.integral_const_mul] at this
  simpa [bJ] using this

lemma bJ0 (n : ℕ) : bJ n 0 = 1 / (n + 1) := by
  simp [bJ, integral_pow]

lemma bJ1 (n : ℕ) : bJ n 1 = -1 / (n + 1) ^ 2 := by
  have := bJ_rec n 0; rw [bJ0] at this; push_cast at this
  have hn : (n:ℝ) + 1 ≠ 0 := by positivity
  field_simp; field_simp at this; linarith

lemma bJ2 (n : ℕ) : bJ n 2 = 2 / (n + 1) ^ 3 := by
  have := bJ_rec n 1; rw [bJ1] at this; push_cast at this
  have hn : (n:ℝ) + 1 ≠ 0 := by positivity
  field_simp; field_simp at this; linarith

/-- `∫ t^P log t log(1-t)` as a series. -/
lemma hasSum_logmulP (P : ℕ) :
    HasSum (fun n : ℕ => -(bJ (P + (n+1)) 1) / (n+1)) (Kb 1 1 P 0) := by
  have key : ∀ t ∈ Ioo (0:ℝ) 1, HasSum (fun n : ℕ => -(t ^ P * log t * (t ^ (n+1) / (n+1))))
      (t ^ P * (1 - t) ^ 0 * (log t ^ 1 * log (1 - t) ^ 1)) := by
    intro t ht
    have h := (hasSum_pow_div_log_of_abs_lt_one
      (show |t| < 1 by rw [abs_of_pos ht.1]; exact ht.2)).mul_left (t ^ P * log t)
    convert h.neg using 1
    ring
  have hint : ∀ n : ℕ, ∫ t in (0:ℝ)..1, -(t ^ P * log t * (t ^ (n+1) / (n+1)))
      = -(bJ (P + (n+1)) 1) / (n+1) := by
    intro n
    rw [intervalIntegral.integral_neg]
    have : (fun t : ℝ => t ^ P * log t * (t ^ (n+1) / (n+1))) =
        fun t => (1 / (n+1 : ℝ)) * (t ^ (P + (n+1)) * log t ^ 1) := by funext t; rw [pow_add]; ring
    rw [this, intervalIntegral.integral_const_mul, bJ]; ring
  simp_rw [← hint]
  have hg : IntervalIntegrable (fun t : ℝ => |log t| * (-log (1 - t))) volume 0 1 := by
    refine bl_ii_of_ls (by fun_prop : Measurable (fun t : ℝ => |log t| * (-log (1 - t)))).aestronglyMeasurable
      1 2 fun t ht => ?_
    rw [abs_mul, abs_abs, abs_neg, one_mul, sq]
    have a1 : |log t| ≤ bls t := by unfold bls; linarith [abs_nonneg (log (1 - t))]
    have a2 : |log (1 - t)| ≤ bls t := by unfold bls; linarith [abs_nonneg (log t)]
    have := bls_nonneg t
    gcongr
  refine intervalIntegral.hasSum_integral_of_dominated_convergence
    (fun n t => |log t| * (t ^ (n+1) / (n+1))) (fun n => ?_) (fun n => ?_) ?_ ?_ ?_
  · exact (by fun_prop : Measurable (fun t : ℝ => -(t ^ P * log t * (t ^ (n+1) / (n+1))))).aestronglyMeasurable
  · filter_upwards [bl_ae_Ioo01'] with t ht' ht0
    have ht := ht' ht0
    have hp : 0 < t ^ (n+1) / ((n:ℝ)+1) := by have := ht.1; positivity
    have hP : 0 < t ^ P := pow_pos ht.1 P
    have hP1 : t ^ P ≤ 1 := pow_le_one₀ ht.1.le ht.2.le
    rw [Real.norm_eq_abs, abs_neg, abs_mul, abs_mul, abs_of_pos hp, abs_of_pos hP]
    have := abs_nonneg (log t)
    have : t ^ P * |log t| ≤ |log t| := by nlinarith
    nlinarith
  · filter_upwards [bl_ae_Ioo01'] with t ht' ht0
    have ht := ht' ht0
    exact ((hasSum_pow_div_log_of_abs_lt_one
      (show |t| < 1 by rw [abs_of_pos ht.1]; exact ht.2)).mul_left _).summable
  · refine hg.congr_uIoo ?_
    rw [bl_uIoo01]
    intro t ht
    dsimp only
    rw [tsum_mul_left, (hasSum_pow_div_log_of_abs_lt_one
      (show |t| < 1 by rw [abs_of_pos ht.1]; exact ht.2)).tsum_eq]
  · filter_upwards [bl_ae_Ioo01'] with t ht' ht0
    have ht := ht' ht0
    exact key t ht

/-- telescoping -/
lemma bl_tele {a : ℝ} (ha : 0 < a) :
    HasSum (fun n : ℕ => 1 / ((n:ℝ) + a) - 1 / ((n:ℝ) + a + 1)) (1 / a) := by
  rw [hasSum_iff_tendsto_nat_of_nonneg]
  · have e : ∀ N : ℕ, ∑ i ∈ Finset.range N, (1 / ((i:ℝ) + a) - 1 / ((i:ℝ) + a + 1)) =
        1 / a - 1 / ((N:ℝ) + a) := by
      intro N
      have := Finset.sum_range_sub' (fun i : ℕ => 1 / ((i:ℝ) + a)) N
      simp only [Nat.cast_add, Nat.cast_one, Nat.cast_zero, zero_add] at this
      rw [← this]
      refine Finset.sum_congr rfl fun i _ => ?_
      ring_nf
    simp_rw [e]
    have : Tendsto (fun N : ℕ => 1 / ((N:ℝ) + a)) atTop (𝓝 0) := by
      have := tendsto_one_div_add_atTop_nhds_zero_nat (𝕜 := ℝ)
      refine tendsto_of_tendsto_of_tendsto_of_le_of_le' tendsto_const_nhds
        (tendsto_one_div_atTop_nhds_zero_nat (𝕜 := ℝ)) ?_ ?_
      · filter_upwards with N; positivity
      · filter_upwards [eventually_gt_atTop 0] with N hN
        apply one_div_le_one_div_of_le (by positivity); linarith
    simpa using this.const_sub (1 / a)
  · intro n
    rw [sub_nonneg]
    apply one_div_le_one_div_of_le (by positivity); linarith

lemma bl_harm (p : ℕ) :
    HasSum (fun n : ℕ => 1 / ((n:ℝ) + 1) - 1 / ((n:ℝ) + 1 + p)) (hh p) := by
  induction p with
  | zero => simp [hh]
  | succ p ih =>
    have h2 := bl_tele (a := (p:ℝ) + 1) (by positivity)
    rw [hh_succ]
    convert ih.add h2 using 1
    · funext n; push_cast; ring_nf

lemma bl_sq_tail (p : ℕ) :
    HasSum (fun n : ℕ => 1 / ((n:ℝ) + 1 + p) ^ 2) (π ^ 2 / 6 - hh2 p) := by
  have := (hasSum_nat_add_iff' (p + 1)).mpr hasSum_zeta_two
  rw [Finset.sum_range_succ'] at this
  simp only [Nat.cast_zero, ne_eq, OfNat.ofNat_ne_zero, not_false_eq_true, zero_pow, div_zero,
    add_zero] at this
  convert this using 1
  · funext n; push_cast; ring_nf
  · simp [hh2]

/-- `∫ y^P log y log(1-y) = (hh p/p + hh2 p - ζ(2))/p`, `p = P+1`. -/
lemma Kb11_base (P : ℕ) : Kb 1 1 P 0 =
    (hh (P + 1) / ((P:ℝ) + 1) + hh2 (P + 1) - π ^ 2 / 6) / ((P:ℝ) + 1) := by
  have h := hasSum_logmulP P
  refine h.unique ?_
  have hp : (0:ℝ) < (P:ℝ) + 1 := by positivity
  have := ((bl_harm (P + 1)).mul_left (1 / ((P:ℝ) + 1) ^ 2)).sub
    ((bl_sq_tail (P + 1)).mul_left (1 / ((P:ℝ) + 1)))
  convert this using 1
  · funext n; rw [bJ1]; push_cast
    have : (n:ℝ) + 1 ≠ 0 := by positivity
    have : (P:ℝ) + 1 ≠ 0 := by positivity
    have : ((P + (n + 1) : ℕ) : ℝ) + 1 ≠ 0 := by positivity
    push_cast at this
    field_simp
    ring
  · field_simp
    ring

noncomputable def F20 (P Q : ℕ) : ℝ := (hh P - hh (P + Q + 1)) ^ 2 + hh2 (P + Q + 1) - hh2 P
noncomputable def F11 (P Q : ℕ) : ℝ :=
  (hh P - hh (P + Q + 1)) * (hh Q - hh (P + Q + 1)) + hh2 (P + Q + 1) - π ^ 2 / 6

lemma Kb00_zero (P : ℕ) : Kb 0 0 P 0 = 1 / ((P : ℝ) + 1) := by
  rw [Kb00]; simp [Nat.factorial_succ]; field_simp

lemma Kb20_eq (P Q : ℕ) : Kb 2 0 P Q = Kb 0 0 P Q * F20 P Q := by
  induction Q generalizing P with
  | zero =>
    have e : Kb 2 0 P 0 = bJ P 2 := by simp [Kb, bJ]
    rw [e, bJ2, Kb00_zero, F20, show P + 0 + 1 = P + 1 by ring, hh_succ, hh2_succ]
    have : (P : ℝ) + 1 ≠ 0 := by positivity
    field_simp; ring
  | succ Q ih =>
    rw [Kb_succ, ih, ih, Kb00_succR, Kb00_succL, F20, F20, F20,
      show P + (Q + 1) + 1 = (P + Q + 1) + 1 by ring, show P + 1 + Q + 1 = (P + Q + 1) + 1 by ring,
      hh_succ (P + Q + 1), hh2_succ (P + Q + 1), hh_succ P, hh2_succ P]
    push_cast
    have : (P : ℝ) + 1 ≠ 0 := by positivity
    have : (P : ℝ) + Q + 2 ≠ 0 := by positivity
    have : (P : ℝ) + Q + 1 + 1 ≠ 0 := by positivity
    field_simp
    ring

lemma Kb11_eq (P Q : ℕ) : Kb 1 1 P Q = Kb 0 0 P Q * F11 P Q := by
  induction Q generalizing P with
  | zero =>
    rw [Kb11_base, Kb00_zero, F11, show P + 0 + 1 = P + 1 by ring, hh_succ P]
    simp only [hh, Finset.range_zero, Finset.sum_empty]
    have : (P : ℝ) + 1 ≠ 0 := by positivity
    field_simp; ring
  | succ Q ih =>
    rw [Kb_succ, ih, ih, Kb00_succR, Kb00_succL, F11, F11, F11,
      show P + (Q + 1) + 1 = (P + Q + 1) + 1 by ring, show P + 1 + Q + 1 = (P + Q + 1) + 1 by ring,
      hh_succ (P + Q + 1), hh2_succ (P + Q + 1), hh_succ P, hh_succ Q]
    push_cast
    have : (P : ℝ) + 1 ≠ 0 := by positivity
    have : (Q : ℝ) + 1 ≠ 0 := by positivity
    have : (P : ℝ) + Q + 2 ≠ 0 := by positivity
    have : (P : ℝ) + Q + 1 + 1 ≠ 0 := by positivity
    field_simp
    ring

lemma Kb02_eq (P Q : ℕ) : Kb 0 2 P Q = Kb 0 0 Q P * F20 Q P := by
  rw [Kb_symm, Kb20_eq]


lemma Zm1_even (n : ℕ) : Zm1 (2 * n) = hh n - hh (2 * n) := by
  induction n with
  | zero => simp [hh]
  | succ n ih =>
    rw [show 2 * (n + 1) = 2 * n + 1 + 1 by ring, Zm1_succ, Zm1_succ, ih, hh_succ, hh_succ, hh_succ]
    have o1 : (-1 : ℝ) ^ (2 * n + 1) = -1 := by rw [pow_succ, pow_mul]; simp
    have o2 : (-1 : ℝ) ^ (2 * n + 1 + 1) = 1 := by rw [pow_succ, o1]; simp
    rw [o1, o2]; push_cast
    have : (n : ℝ) + 1 ≠ 0 := by positivity
    have : 2 * (n : ℝ) + 1 ≠ 0 := by positivity
    have : 2 * (n : ℝ) + 1 + 1 ≠ 0 := by positivity
    field_simp; ring

lemma Zm1_even' (n : ℕ) : Zm1 (2 * n + 1 + 1) = hh (n + 1) - hh (2 * n + 1 + 1) := by
  have := Zm1_even (n + 1); rwa [show 2 * (n + 1) = 2 * n + 1 + 1 by ring] at this

lemma Z2_eq (M : ℕ) : Z2 M = hh2 M := by
  induction M with
  | zero => simp [hh2]
  | succ M ih => rw [Z2_succ, hh2_succ, ih]

lemma Zm2_even (n : ℕ) : Zm2 (2 * n) = hh2 n / 2 - hh2 (2 * n) := by
  induction n with
  | zero => simp [hh2]
  | succ n ih =>
    rw [show 2 * (n + 1) = 2 * n + 1 + 1 by ring, Zm2_succ, Zm2_succ, ih, hh2_succ, hh2_succ,
      hh2_succ]
    have o1 : (-1 : ℝ) ^ (2 * n + 1) = -1 := by rw [pow_succ, pow_mul]; simp
    have o2 : (-1 : ℝ) ^ (2 * n + 1 + 1) = 1 := by rw [pow_succ, o1]; simp
    rw [o1, o2]; push_cast
    have : (n : ℝ) + 1 ≠ 0 := by positivity
    have : 2 * (n : ℝ) + 1 ≠ 0 := by positivity
    have : 2 * (n : ℝ) + 1 + 1 ≠ 0 := by positivity
    field_simp; ring

lemma Zm2_even' (n : ℕ) : Zm2 (2 * n + 1 + 1) = hh2 (n + 1) / 2 - hh2 (2 * n + 1 + 1) := by
  have := Zm2_even (n + 1); rwa [show 2 * (n + 1) = 2 * n + 1 + 1 by ring] at this

lemma Zm1m1_eq (M : ℕ) : Zm1m1 M = (Zm1 M ^ 2 - Z2 M) / 2 := by
  induction M with
  | zero => simp
  | succ M ih =>
    rw [Zm1m1_succ, Zm1_succ, Z2_succ, ih]
    have hs : ((-1 : ℝ) ^ (M + 1)) ^ 2 = 1 := by rw [← pow_mul, pow_mul']; simp
    have : (M : ℝ) + 1 ≠ 0 := by positivity
    rw [add_sq, div_pow, hs]
    ring

lemma coef_alg (N : ℕ) :
    ((N:ℝ) ^ 2 / 2) * Kb 2 0 N N + (-(2 * (N:ℝ) + 1) ^ 2 / 2) * Kb 2 0 (N + 1) N
      + ((2 * (N:ℝ) + 1) * ((N:ℝ) + 1)) * Kb 2 0 (N + 2) N
      + (-(N:ℝ) * ((N:ℝ) + 1)) * Kb 1 1 N N + (2 * (2 * (N:ℝ) + 1) * ((N:ℝ) + 1)) * Kb 1 1 (N + 1) N
      + (-2 * (2 * (N:ℝ) + 1) * ((N:ℝ) + 1)) * Kb 1 1 (N + 2) N
      + (((N:ℝ) + 1) ^ 2 / 2) * Kb 0 2 N N + (-(2 * (N:ℝ) + 1) * (2 * (N:ℝ) + 3) / 2) * Kb 0 2 (N + 1) N
      + ((2 * (N:ℝ) + 1) * ((N:ℝ) + 1)) * Kb 0 2 (N + 2) N
    = -(((N:ℝ) + 1) * (2 * (N:ℝ) + 1)) * (r (N + 1) + (π ^ 2 / 6) / (4 * ((N:ℝ) + 1) ^ 2 - 1))
      * Kb 0 0 N N := by
  rw [Kb20_eq, Kb20_eq, Kb20_eq, Kb11_eq, Kb11_eq, Kb11_eq, Kb02_eq, Kb02_eq, Kb02_eq]
  rw [show N + 2 = N + 1 + 1 by ring, Kb00_succL (N + 1) N, Kb00_succL N N,
    Kb00_succR N (N + 1), Kb00_succR N N]
  simp only [F20, F11]
  rw [show N + N + 1 = 2 * N + 1 by ring, show N + 1 + N + 1 = 2 * N + 1 + 1 by ring,
    show N + 1 + 1 + N + 1 = 2 * N + 1 + 1 + 1 by ring, show N + (N + 1) + 1 = 2 * N + 1 + 1 by ring,
    show N + (N + 1 + 1) + 1 = 2 * N + 1 + 1 + 1 by ring]
  simp only [r]
  rw [show 2 * (N + 1) = 2 * N + 1 + 1 by ring, Zm1m1_eq, Zm1_even', Z2_eq, Zm2_even']
  simp only [hh_succ, hh2_succ]
  push_cast
  set A := hh N
  set C := hh (2 * N)
  set A2 := hh2 N
  set C2 := hh2 (2 * N)
  set B := Kb 0 0 N N
  have h1 : (N:ℝ) + 1 ≠ 0 := by positivity
  have h2 : (N:ℝ) + 1 + 1 ≠ 0 := by positivity
  have h3 : 2 * (N:ℝ) + 1 ≠ 0 := by positivity
  have h4 : 2 * (N:ℝ) + 1 + 1 ≠ 0 := by positivity
  have h5 : 2 * (N:ℝ) + 1 + 1 + 1 ≠ 0 := by positivity
  have h6 : (N:ℝ) + N + 2 ≠ 0 := by positivity
  have h7 : (N:ℝ) + 1 + N + 2 ≠ 0 := by positivity
  have h8 : (N:ℝ) + (N + 1) + 2 ≠ 0 := by positivity
  have h0 : (0:ℝ) ≤ N := N.cast_nonneg
  have h9 : 4 * ((N:ℝ) + 1) ^ 2 - 1 ≠ 0 := by nlinarith
  have h10 : 2 * ((N:ℝ) + 1) - 1 ≠ 0 := by nlinarith
  have h11 : 2 * ((N:ℝ) + 1) + 1 ≠ 0 := by nlinarith
  have h12 : 2 * (N:ℝ) + 2 + 1 ≠ 0 := by positivity
  have h13 : 2 * (N:ℝ) + 3 ≠ 0 := by positivity
  rw [show 4 * ((N:ℝ) + 1) ^ 2 - 1 = (2 * N + 1) * (2 * N + 3) by ring]
  field_simp
  ring

/-- `E_{n,2}(y)` of Lemma `g6:lem:En` for `n = N+1`: `½X² + γℓ²`, `X = α₁ log y + β₁ log(1-y)`. -/
noncomputable def Ey (N : ℕ) (y : ℝ) : ℝ :=
  ((-(N:ℝ) + (2 * N + 1) * y) * log y + (((N:ℝ) + 1) - (2 * N + 1) * y) * log (1 - y)) ^ 2 / 2
    - (2 * (N:ℝ) + 1) * y * (1 - y) / 2 * (log (1 - y) - log y) ^ 2

/-- the integrand of `Kb`. -/
noncomputable def kf (i m P Q : ℕ) (y : ℝ) : ℝ := y ^ P * (1 - y) ^ Q * (log y ^ i * log (1 - y) ^ m)

lemma ii_kf (i m P Q : ℕ) (c : ℝ) : IntervalIntegrable (fun y => c * kf i m P Q y) volume 0 1 :=
  (ii_Kb i m P Q).const_mul c

lemma psi_decomp (N : ℕ) : (fun y => (y * (1 - y)) ^ N * Ey N y) = fun y =>
    ((N:ℝ) ^ 2 / 2) * kf 2 0 N N y + (-(2 * (N:ℝ) + 1) ^ 2 / 2) * kf 2 0 (N + 1) N y
    + ((2 * (N:ℝ) + 1) * ((N:ℝ) + 1)) * kf 2 0 (N + 2) N y
    + (-(N:ℝ) * ((N:ℝ) + 1)) * kf 1 1 N N y + (2 * (2 * (N:ℝ) + 1) * ((N:ℝ) + 1)) * kf 1 1 (N + 1) N y
    + (-2 * (2 * (N:ℝ) + 1) * ((N:ℝ) + 1)) * kf 1 1 (N + 2) N y
    + (((N:ℝ) + 1) ^ 2 / 2) * kf 0 2 N N y + (-(2 * (N:ℝ) + 1) * (2 * (N:ℝ) + 3) / 2) * kf 0 2 (N + 1) N y
    + ((2 * (N:ℝ) + 1) * ((N:ℝ) + 1)) * kf 0 2 (N + 2) N y := by
  funext y; simp only [Ey, kf]; rw [mul_pow]; ring

lemma ii_psi (N : ℕ) : IntervalIntegrable (fun y => (y * (1 - y)) ^ N * Ey N y) volume 0 1 := by
  rw [psi_decomp N]
  repeat apply IntervalIntegrable.add
  all_goals exact ii_kf _ _ _ _ _

/-- The coefficient identity: `½∫_0^1 [y(1-y)]^{n-1} E_{n,2} dy`, normalised by `B(n,n)`, is
`-n(2n-1)(r_n + ζ(2)/(4n²-1))` (`n = N+1`). -/
lemma coef_main (N : ℕ) : ∫ y in (0:ℝ)..1, (y * (1 - y)) ^ N * Ey N y
    = -(((N:ℝ) + 1) * (2 * (N:ℝ) + 1)) * (r (N + 1) + (π ^ 2 / 6) / (4 * ((N:ℝ) + 1) ^ 2 - 1))
      * Kb 0 0 N N := by
  rw [← coef_alg N]
  have e := psi_decomp N
  rw [e]
  rw [intervalIntegral.integral_add, intervalIntegral.integral_add, intervalIntegral.integral_add,
    intervalIntegral.integral_add, intervalIntegral.integral_add, intervalIntegral.integral_add,
    intervalIntegral.integral_add, intervalIntegral.integral_add]
  · simp only [intervalIntegral.integral_const_mul, Kb, kf]
  all_goals first
    | exact ii_kf _ _ _ _ _
    | (repeat apply IntervalIntegrable.add) <;> exact ii_kf _ _ _ _ _

/-- `ρ_n(g) = g^n (1+g)^{1-2n}`, `n = N+1`, and its first two derivatives. -/
noncomputable def rho (N : ℕ) (g : ℝ) : ℝ := g ^ (N + 1) / (1 + g) ^ (2 * N + 1)
noncomputable def rho1 (N : ℕ) (g : ℝ) : ℝ := rho N g * (((N:ℝ) + 1) - N * g) / (g * (1 + g))
noncomputable def rho2 (N : ℕ) (g : ℝ) : ℝ :=
  rho N g * ((((N:ℝ) + 1) / g - (2 * N + 1) / (1 + g)) ^ 2 - ((N:ℝ) + 1) / g ^ 2
    + (2 * (N:ℝ) + 1) / (1 + g) ^ 2)

/-- the linear operator of `M2integrand` applied to `ρ_n` (without the constant `-1/(2g)`). -/
noncomputable def Lam (N : ℕ) (g : ℝ) : ℝ :=
  (rho N g - g * rho1 N g + g ^ 2 * rho2 N g) / (2 * g * (1 + g)) * log g ^ 2
    + (g * rho1 N g - rho N g) / (g * (1 + g)) * log g * log (1 + g)
    + rho N g / (2 * g * (1 + g)) * log (1 + g) ^ 2

lemma psi_eq (N : ℕ) {y : ℝ} (h0 : 0 < y) (h1 : y < 1) :
    Lam N (y / (1 - y)) * (1 / (1 - y) ^ 2) = (y * (1 - y)) ^ N * Ey N y := by
  have hy : 0 < 1 - y := by linarith
  have e1 : 1 + y / (1 - y) = 1 / (1 - y) := by field_simp; ring
  have eL : log (y / (1 - y)) = log y - log (1 - y) := log_div h0.ne' hy.ne'
  have em : log (1 + y / (1 - y)) = -log (1 - y) := by rw [e1, one_div, log_inv]
  have er : rho N (y / (1 - y)) = y ^ (N + 1) * (1 - y) ^ N := by
    rw [rho, e1, div_pow, one_div_pow, div_div_eq_mul_div, div_one]
    rw [show 2 * N + 1 = (N + 1) + N by ring, pow_add]
    have : (1 - y) ^ (N + 1) ≠ 0 := by positivity
    field_simp
    rw [pow_add]
  have em' : log (1 / (1 - y)) = -log (1 - y) := by rw [one_div, log_inv]
  simp only [Lam, rho1, rho2, er, eL, e1, em', Ey]
  rw [mul_pow]
  field_simp
  ring

lemma Ey_symm (N : ℕ) (y : ℝ) : Ey N (1 - y) = Ey N y := by
  simp only [Ey, sub_sub_cancel]; ring


lemma int_Lam (N : ℕ) :
    ∫ g in (0:ℝ)..1, Lam N g = (1 / 2) * ∫ y in (0:ℝ)..1, (y * (1 - y)) ^ N * Ey N y := by
  set ψ := fun y : ℝ => (y * (1 - y)) ^ N * Ey N y with hψ
  have hc : ContinuousOn (fun y : ℝ => y / (1 - y)) (uIcc 0 (1 / 2)) := by
    have h : ∀ x ∈ uIcc (0:ℝ) (1 / 2), (1 - x) ≠ 0 := by
      intro x hx; rw [uIcc_of_le (by norm_num)] at hx; have := hx.2; intro h; linarith
    exact ContinuousOn.div continuousOn_id (continuousOn_const.sub continuousOn_id) h
  have hd : ∀ x ∈ Ioo (min (0:ℝ) (1 / 2)) (max 0 (1 / 2)),
      HasDerivAt (fun y : ℝ => y / (1 - y)) (1 / (1 - x) ^ 2) x := by
    intro x hx
    rw [min_eq_left (by norm_num), max_eq_right (by norm_num)] at hx
    have h1 : (1 - x) ≠ 0 := by linarith [hx.2]
    have := (hasDerivAt_id' x).div ((hasDerivAt_id' x).const_sub 1) h1
    convert this using 1
    field_simp; ring
  have hp : ∀ x ∈ Ioo (min (0:ℝ) (1 / 2)) (max 0 (1 / 2)), 0 ≤ 1 / (1 - x) ^ 2 :=
    fun x _ => by positivity
  have hcv := intervalIntegral.integral_comp_mul_deriv_of_deriv_nonneg (g := Lam N) hc hd hp
  norm_num at hcv
  have h1 : ∫ y in (0:ℝ)..1 / 2, ψ y = ∫ g in (0:ℝ)..1, Lam N g := by
    rw [← hcv]
    refine intervalIntegral.integral_congr fun y hy => ?_
    rw [uIcc_of_le (by norm_num)] at hy
    rcases eq_or_lt_of_le hy.1 with h | h
    · subst h; simp [ψ, Ey, Lam]
    · simp only [ψ]
      rw [← psi_eq N h (by linarith [hy.2])]
      simp
  have h2 : ∫ y in (1/2:ℝ)..1, ψ y = ∫ y in (0:ℝ)..1 / 2, ψ y := by
    have := intervalIntegral.integral_comp_sub_left ψ (a := 0) (b := 1 / 2) 1
    rw [show (1:ℝ) - 1 / 2 = 1 / 2 by norm_num, sub_zero] at this
    rw [← this]
    congr 1; funext y; simp only [ψ, Ey_symm]; ring
  have hi := ii_psi N
  have hs := intervalIntegral.integral_add_adjacent_intervals (a := 0) (b := 1 / 2) (c := 1)
    (hi.mono_set (by
      rw [uIcc_of_le (by norm_num), uIcc_of_le (by norm_num)]
      exact Icc_subset_Icc le_rfl (by norm_num)))
    (hi.mono_set (by
      rw [uIcc_of_le (by norm_num), uIcc_of_le (by norm_num)]
      exact Icc_subset_Icc (by norm_num) le_rfl))
  rw [← hs, h2, h1]; ring

lemma hasSum_sqrt_one_sub {x : ℝ} (h0 : 0 ≤ x) (h1 : x < 1) :
    HasSum (fun n : ℕ => Ring.choose (1 / 2 : ℝ) n * (-x) ^ n) (Real.sqrt (1 - x)) := by
  have hb := (Complex.one_add_cpow_hasFPowerSeriesOnBall_zero (a := (1 / 2 : ℂ))).hasSum
    (y := ((-x : ℝ) : ℂ)) (by
      rw [Metric.mem_eball, edist_zero_right]
      simp only [Complex.ofReal_neg, eNorm_neg]
      rw [← ofReal_norm, ENNReal.ofReal_lt_one, Complex.norm_real, Real.norm_eq_abs,
        abs_of_nonneg h0]
      exact h1)
  simp only [binomialSeries_apply, List.ofFn_const, List.prod_replicate, zero_add] at hb
  rw [← Complex.hasSum_ofReal]
  convert hb using 1
  · funext n
    rw [smul_eq_mul]
    push_cast [Complex.ofReal_choose]
    rfl
  · rw [Real.sqrt_eq_rpow, Complex.ofReal_cpow (by linarith)]
    push_cast
    ring_nf

lemma choose_half_succ (n : ℕ) :
    Ring.choose (1 / 2 : ℝ) (n + 1) = Ring.choose (1 / 2 : ℝ) n * ((1 / 2 - n) / (n + 1)) := by
  rw [Ring.choose_eq_smul, Ring.choose_eq_smul, descPochhammer_succ_right]
  rw [Polynomial.smeval_mul]
  simp [Nat.factorial_succ, Polynomial.smeval_natCast]
  have : (n : ℝ) + 1 ≠ 0 := by positivity
  field_simp

/-- `κ_n = a_n 4^{-n}` (`n = N+1`, `a_n = binom(2n,n)/(2n-1)`), written through the Beta value
`B(n,n)`: `a_n B(n,n)/2 = 1/(n(2n-1))`. -/
noncomputable def kap (N : ℕ) : ℝ :=
  2 / (4 ^ (N + 1) * (((N:ℝ) + 1) * (2 * (N:ℝ) + 1)) * Kb 0 0 N N)

lemma choose_half_eq (N : ℕ) :
    Ring.choose (1 / 2 : ℝ) (N + 1) * (-1) ^ (N + 1) = -kap N := by
  induction N with
  | zero =>
    rw [Ring.choose_one_right]; simp [kap, Kb00]; norm_num
  | succ N ih =>
    rw [choose_half_succ, pow_succ, show
      Ring.choose (1 / 2 : ℝ) (N + 1) * ((1 / 2 - ((N + 1 : ℕ) : ℝ)) / ((N + 1 : ℕ) + 1)) *
        ((-1) ^ (N + 1) * -1) = (Ring.choose (1 / 2 : ℝ) (N + 1) * (-1) ^ (N + 1)) *
        (-(1 / 2 - ((N + 1 : ℕ) : ℝ)) / ((N + 1 : ℕ) + 1)) by ring, ih]
    simp only [kap]
    rw [Kb00_succR (N + 1) N, Kb00_succL N N]
    push_cast
    have hB := Kb00_pos N N
    have : (N : ℝ) + 1 ≠ 0 := by positivity
    have : (N : ℝ) + 2 ≠ 0 := by positivity
    have : 2 * (N : ℝ) + 3 ≠ 0 := by positivity
    have : (N : ℝ) + N + 2 ≠ 0 := by positivity
    have : (N : ℝ) + 1 + N + 2 ≠ 0 := by positivity
    have : (2 : ℝ) * (N + 1) + 1 ≠ 0 := by positivity
    rw [pow_succ]
    field_simp
    ring

/-- the binomial series `1 - √(1-x) = Σ_{n≥1} κ_n x^n`, `0 ≤ x < 1`. -/
lemma hasSum_kap {x : ℝ} (h0 : 0 ≤ x) (h1 : x < 1) :
    HasSum (fun N : ℕ => kap N * x ^ (N + 1)) (1 - Real.sqrt (1 - x)) := by
  have h := (hasSum_nat_add_iff' 1).mpr (hasSum_sqrt_one_sub h0 h1)
  simp only [Finset.range_one, Finset.sum_singleton, Ring.choose_zero_right, pow_zero,
    mul_one] at h
  convert h.neg using 1
  · funext N
    rw [neg_pow, ← mul_assoc, choose_half_eq]; ring
  · ring

lemma hasDerivAt_rho (N : ℕ) {g : ℝ} (hg : 0 < g) : HasDerivAt (rho N) (rho1 N g) g := by
  have h1 : (1 + g) ≠ 0 := by positivity
  have := (hasDerivAt_pow (N + 1) g).div ((hasDerivAt_id' g).const_add 1 |>.pow (2 * N + 1))
    (pow_ne_zero _ h1)
  have e : rho N = fun z => z ^ (N + 1) / (1 + z) ^ (2 * N + 1) := rfl
  rw [e]
  convert this using 1
  simp only [rho1, rho, Nat.add_sub_cancel]
  have hg' : g ≠ 0 := hg.ne'
  simp only [Pi.pow_apply]
  rw [pow_succ (1 + g) (2 * N)]
  push_cast
  field_simp
  ring

lemma rho1_eq (N : ℕ) {g : ℝ} (hg : 0 < g) :
    rho1 N g = rho N g * (((N:ℝ) + 1) / g - (2 * N + 1) / (1 + g)) := by
  simp only [rho1]
  have : (1 + g) ≠ 0 := by positivity
  field_simp; ring

lemma hasDerivAt_rho1 (N : ℕ) {g : ℝ} (hg : 0 < g) : HasDerivAt (rho1 N) (rho2 N g) g := by
  have h1 : (1 + g) ≠ 0 := by positivity
  have hev : rho1 N =ᶠ[𝓝 g] fun z => rho N z * (((N:ℝ) + 1) / z - (2 * N + 1) / (1 + z)) := by
    filter_upwards [Ioi_mem_nhds hg] with z hz using rho1_eq N hz
  refine HasDerivAt.congr_of_eventuallyEq ?_ hev
  have ht : HasDerivAt (fun z : ℝ => ((N:ℝ) + 1) / z - (2 * N + 1) / (1 + z))
      (-((N:ℝ) + 1) / g ^ 2 + (2 * N + 1) / (1 + g) ^ 2) g := by
    have a := (hasDerivAt_inv hg.ne').const_mul ((N:ℝ) + 1)
    have b := ((hasDerivAt_id' g).const_add 1).inv h1 |>.const_mul (2 * (N:ℝ) + 1)
    convert a.sub b using 1
    · funext z; simp [div_eq_mul_inv]
    · field_simp; ring
  convert (hasDerivAt_rho N hg).mul ht using 1
  rw [rho2, rho1_eq N hg]; ring

lemma q_le {g : ℝ} (hg : 0 < g) : g / (1 + g) ^ 2 ≤ 1 / 4 := by
  rw [div_le_div_iff₀ (by positivity) (by norm_num)]; nlinarith [sq_nonneg (g - 1)]

lemma q_nonneg {g : ℝ} (hg : 0 < g) : 0 ≤ g / (1 + g) ^ 2 := by positivity

lemma rho_eq_q (N : ℕ) {g : ℝ} (hg : 0 < g) : rho N g = (g / (1 + g) ^ 2) ^ (N + 1) * (1 + g) := by
  have : (1 + g) ≠ 0 := by positivity
  rw [rho, div_pow, ← pow_mul, show 2 * (N + 1) = (2 * N + 1) + 1 by ring, pow_succ]
  field_simp
  ring

lemma rho1_cf (N : ℕ) {g : ℝ} (hg : 0 < g) :
    rho1 N g = (g / (1 + g) ^ 2) ^ N * (((N:ℝ) + 1) - N * g) / (1 + g) ^ 2 := by
  have : (1 + g) ≠ 0 := by positivity
  rw [rho1, rho_eq_q N hg, pow_succ]
  field_simp

lemma rho2_cf (N : ℕ) {g : ℝ} (hg : 0 < g) :
    rho2 N g * (g * (1 + g) ^ 3) = (g / (1 + g) ^ 2) ^ N *
      ((N:ℝ) * (N + 1) * (1 + g ^ 2) - 2 * ((N:ℝ) + 1) ^ 2 * g) := by
  have : (1 + g) ≠ 0 := by positivity
  rw [rho2, rho_eq_q N hg, pow_succ]
  field_simp
  ring

lemma pow_q_le (N : ℕ) {g : ℝ} (hg : 0 < g) : (g / (1 + g) ^ 2) ^ N ≤ 1 / 4 ^ N := by
  rw [one_div, ← inv_pow]
  exact pow_le_pow_left₀ (q_nonneg hg) (by have := q_le hg; linarith) N

lemma abs_rho1_le (N : ℕ) {g : ℝ} (hg : 0 < g) : |rho1 N g| ≤ ((N:ℝ) + 1) / 4 ^ N := by
  rw [rho1_cf N hg, abs_div, abs_mul, abs_of_nonneg (pow_nonneg (q_nonneg hg) N),
    abs_of_pos (by positivity : (0:ℝ) < (1 + g) ^ 2)]
  have h1 : |((N:ℝ) + 1) - N * g| ≤ ((N:ℝ) + 1) * (1 + g) := by
    rw [abs_le]; constructor <;> nlinarith [N.cast_nonneg (α := ℝ)]
  have h2 := pow_q_le N hg
  have h3 : (1:ℝ) + g ≤ (1 + g) ^ 2 := by nlinarith
  rw [div_le_iff₀ (by positivity)]
  calc (g / (1 + g) ^ 2) ^ N * |((N:ℝ) + 1) - N * g|
      ≤ (1 / 4 ^ N) * (((N:ℝ) + 1) * (1 + g)) := by
        gcongr
    _ ≤ (1 / 4 ^ N) * (((N:ℝ) + 1) * (1 + g) ^ 2) := by gcongr
    _ = ((N:ℝ) + 1) / 4 ^ N * (1 + g) ^ 2 := by ring

lemma abs_rho2_le (N : ℕ) {g : ℝ} (hg : 0 < g) : |rho2 N g| ≤ 8 * ((N:ℝ) + 1) ^ 2 / 4 ^ N := by
  have hc := rho2_cf N hg
  have hp : 0 < g * (1 + g) ^ 3 := by positivity
  have e : rho2 N g = (g / (1 + g) ^ 2) ^ N *
      ((N:ℝ) * (N + 1) * (1 + g ^ 2) - 2 * ((N:ℝ) + 1) ^ 2 * g) / (g * (1 + g) ^ 3) := by
    rw [← hc]; field_simp
  rcases N with _ | M
  · rw [e, abs_div, abs_of_pos hp, div_le_iff₀ hp]
    have h3 : |(g / (1 + g) ^ 2) ^ 0 * (((0:ℕ):ℝ) * (((0:ℕ):ℝ) + 1) * (1 + g ^ 2)
        - 2 * (((0:ℕ):ℝ) + 1) ^ 2 * g)| = 2 * g := by
      rw [show (g / (1 + g) ^ 2) ^ 0 * (((0:ℕ):ℝ) * (((0:ℕ):ℝ) + 1) * (1 + g ^ 2)
        - 2 * (((0:ℕ):ℝ) + 1) ^ 2 * g) = -(2 * g) by push_cast; ring, abs_neg,
        abs_of_pos (by linarith)]
    rw [h3]
    have : (1:ℝ) ≤ (1 + g) ^ 3 := one_le_pow₀ (by linarith)
    norm_num
    nlinarith
  · have e2 : rho2 (M + 1) g = (g / (1 + g) ^ 2) ^ M *
        (((M + 1 : ℕ) : ℝ) * ((M + 1 : ℕ) + 1) * (1 + g ^ 2) - 2 * (((M + 1 : ℕ) : ℝ) + 1) ^ 2 * g)
          / (1 + g) ^ 5 := by
      rw [e, pow_succ]; field_simp
    rw [e2, abs_div, abs_mul, abs_of_nonneg (pow_nonneg (q_nonneg hg) M),
      abs_of_pos (by positivity : (0:ℝ) < (1 + g) ^ 5)]
    push_cast
    set x : ℝ := (M : ℝ) + 1 with hx
    have hx0 : 0 < x := by positivity
    have hnum : |x * (x + 1) * (1 + g ^ 2) - 2 * (x + 1) ^ 2 * g| ≤ 2 * (x + 1) ^ 2 * (1 + g) ^ 2 := by
      rw [abs_le]; constructor <;> nlinarith [sq_nonneg g, mul_pos hx0 hg, sq_nonneg (x + 1)]
    have h2 := pow_q_le M hg
    have h5 : (1 + g) ^ 2 ≤ (1 + g) ^ 5 := pow_le_pow_right₀ (by linarith) (by norm_num)
    rw [div_le_iff₀ (by positivity)]
    have : (4:ℝ) ^ (M + 1) = 4 ^ M * 4 := pow_succ 4 M
    rw [this]
    calc (g / (1 + g) ^ 2) ^ M * |x * (x + 1) * (1 + g ^ 2) - 2 * (x + 1) ^ 2 * g|
        ≤ (1 / 4 ^ M) * (2 * (x + 1) ^ 2 * (1 + g) ^ 2) := by gcongr
      _ ≤ (1 / 4 ^ M) * (2 * (x + 1) ^ 2 * (1 + g) ^ 5) := by gcongr
      _ = 8 * (x + 1) ^ 2 / (4 ^ M * 4) * (1 + g) ^ 5 := by field_simp; ring


lemma kap_succ (N : ℕ) : kap (N + 1) = kap N * ((2 * (N:ℝ) + 1) / (2 * ((N:ℝ) + 2))) := by
  simp only [kap]
  rw [Kb00_succR (N + 1) N, Kb00_succL N N]
  push_cast
  have hB := Kb00_pos N N
  have : (N : ℝ) + 1 ≠ 0 := by positivity
  have : (N : ℝ) + 2 ≠ 0 := by positivity
  have : 2 * (N : ℝ) + 3 ≠ 0 := by positivity
  have : (N : ℝ) + N + 2 ≠ 0 := by positivity
  have : (N : ℝ) + 1 + N + 2 ≠ 0 := by positivity
  rw [pow_succ]
  field_simp
  ring

lemma kap_pos (N : ℕ) : 0 < kap N := by
  have := Kb00_pos N N; unfold kap; positivity

lemma kap_le (N : ℕ) : kap N ≤ 1 / 2 := by
  induction N with
  | zero => simp [kap, Kb00]; norm_num
  | succ N ih =>
    rw [kap_succ]
    have h := kap_pos N
    have : (2 * (N:ℝ) + 1) / (2 * ((N:ℝ) + 2)) ≤ 1 := by
      rw [div_le_one (by positivity)]; linarith
    calc kap N * ((2 * (N:ℝ) + 1) / (2 * ((N:ℝ) + 2))) ≤ kap N * 1 := by gcongr
      _ ≤ 1 / 2 := by linarith

/-- `R(g) = √(1 + g² - 2kg)` with `k = 2c² - 1`. -/
noncomputable def Rc (c g : ℝ) : ℝ := Real.sqrt (1 + g ^ 2 - 2 * (2 * c ^ 2 - 1) * g)

/-- the weights `a_n c^{2n}` of the expansion. -/
noncomputable def wt (c : ℝ) (N : ℕ) : ℝ := 4 ^ (N + 1) * kap N * c ^ (2 * N + 2)

lemma Rc_arg_pos {c : ℝ} (hc0 : 0 < c) (hc1 : c < 1) (g : ℝ) :
    0 < 1 + g ^ 2 - 2 * (2 * c ^ 2 - 1) * g := by
  have h1 : 0 < 1 - (2 * c ^ 2 - 1) ^ 2 := by
    have e : 1 - (2 * c ^ 2 - 1) ^ 2 = 4 * c ^ 2 * (1 - c ^ 2) := by ring
    have : c ^ 2 < 1 := by nlinarith
    rw [e]; exact mul_pos (by positivity) (by linarith)
  nlinarith [sq_nonneg (g - (2 * c ^ 2 - 1))]

lemma Rc_eq {c g : ℝ} (hg : 0 < g) :
    Rc c g = (1 + g) * Real.sqrt (1 - 4 * c ^ 2 * (g / (1 + g) ^ 2)) := by
  have h1 : (1 + g) ≠ 0 := by positivity
  have e : 1 + g ^ 2 - 2 * (2 * c ^ 2 - 1) * g = (1 + g) ^ 2 * (1 - 4 * c ^ 2 * (g / (1 + g) ^ 2)) := by
    field_simp; ring
  rw [Rc, e, Real.sqrt_mul (by positivity), Real.sqrt_sq (by linarith)]

lemma hasSum_G {c g : ℝ} (hc0 : 0 < c) (hc1 : c < 1) (hg : 0 < g) :
    HasSum (fun N => wt c N * rho N g) ((1 + g) - Rc c g) := by
  have hq := q_le hg
  have hq0 := q_nonneg hg
  have hx0 : 0 ≤ 4 * c ^ 2 * (g / (1 + g) ^ 2) := by positivity
  have hx1 : 4 * c ^ 2 * (g / (1 + g) ^ 2) < 1 := by
    have : 4 * c ^ 2 * (g / (1 + g) ^ 2) ≤ c ^ 2 := by nlinarith [sq_nonneg c]
    nlinarith
  have h := (hasSum_kap hx0 hx1).mul_left (1 + g)
  convert h using 1
  · funext N
    rw [wt, rho_eq_q N hg, mul_pow, mul_pow, ← pow_mul, show 2 * (N + 1) = 2 * N + 2 by ring]
    ring
  · rw [Rc_eq hg]; ring

lemma summable_bd {c : ℝ} (hc0 : 0 < c) (hc1 : c < 1) :
    Summable (fun N : ℕ => 16 * ((N:ℝ) + 1) ^ 2 * c ^ (2 * N + 2)) := by
  have hr : ‖c ^ 2‖ < 1 := by
    rw [Real.norm_eq_abs, abs_of_pos (by positivity)]; nlinarith
  have h0 : Summable (fun n : ℕ => (n:ℝ) ^ 2 * (c ^ 2) ^ n) :=
    summable_pow_mul_geometric_of_norm_lt_one 2 hr
  have h1 : Summable (fun n : ℕ => ((n + 1 : ℕ) : ℝ) ^ 2 * (c ^ 2) ^ (n + 1)) :=
    (summable_nat_add_iff (f := fun n : ℕ => (n:ℝ) ^ 2 * (c ^ 2) ^ n) 1).mpr h0
  refine (h1.mul_left 16).congr fun N => ?_
  rw [← pow_mul, show 2 * (N + 1) = 2 * N + 2 by ring]; push_cast; ring

lemma wt_nonneg (c : ℝ) (N : ℕ) : 0 ≤ wt c N := by
  have := kap_pos N; unfold wt
  have : 0 ≤ c ^ (2 * N + 2) := by rw [show 2 * N + 2 = 2 * (N + 1) by ring, pow_mul]; positivity
  positivity

lemma wt_le (c : ℝ) (N : ℕ) : wt c N ≤ 4 ^ (N + 1) / 2 * c ^ (2 * N + 2) := by
  have := kap_le N
  have : 0 ≤ c ^ (2 * N + 2) := by rw [show 2 * N + 2 = 2 * (N + 1) by ring, pow_mul]; positivity
  unfold wt
  calc 4 ^ (N + 1) * kap N * c ^ (2 * N + 2) ≤ 4 ^ (N + 1) * (1 / 2) * c ^ (2 * N + 2) := by gcongr
    _ = _ := by ring

lemma abs_wt_rho1_le (c : ℝ) (N : ℕ) {g : ℝ} (hg : 0 < g) :
    ‖wt c N * rho1 N g‖ ≤ 16 * ((N:ℝ) + 1) ^ 2 * c ^ (2 * N + 2) := by
  rw [Real.norm_eq_abs, abs_mul, abs_of_nonneg (wt_nonneg c N)]
  have h1 := abs_rho1_le N hg
  have h2 := wt_le c N
  have hc : 0 ≤ c ^ (2 * N + 2) := by rw [show 2 * N + 2 = 2 * (N + 1) by ring, pow_mul]; positivity
  calc wt c N * |rho1 N g| ≤ (4 ^ (N + 1) / 2 * c ^ (2 * N + 2)) * (((N:ℝ) + 1) / 4 ^ N) := by
        gcongr
    _ = 2 * ((N:ℝ) + 1) * c ^ (2 * N + 2) := by rw [pow_succ]; field_simp; ring
    _ ≤ _ := by
        have : (0:ℝ) ≤ N := N.cast_nonneg
        have h3 : (N:ℝ) + 1 ≤ ((N:ℝ) + 1) ^ 2 := by nlinarith
        nlinarith [mul_le_mul_of_nonneg_right h3 hc]

lemma abs_wt_rho2_le (c : ℝ) (N : ℕ) {g : ℝ} (hg : 0 < g) :
    ‖wt c N * rho2 N g‖ ≤ 16 * ((N:ℝ) + 1) ^ 2 * c ^ (2 * N + 2) := by
  rw [Real.norm_eq_abs, abs_mul, abs_of_nonneg (wt_nonneg c N)]
  have h1 := abs_rho2_le N hg
  have h2 := wt_le c N
  have hc : 0 ≤ c ^ (2 * N + 2) := by rw [show 2 * N + 2 = 2 * (N + 1) by ring, pow_mul]; positivity
  calc wt c N * |rho2 N g| ≤ (4 ^ (N + 1) / 2 * c ^ (2 * N + 2)) * (8 * ((N:ℝ) + 1) ^ 2 / 4 ^ N) := by
        gcongr
    _ = _ := by rw [pow_succ]; field_simp; ring

lemma Rc_pos {c : ℝ} (hc0 : 0 < c) (hc1 : c < 1) (g : ℝ) : 0 < Rc c g :=
  Real.sqrt_pos.mpr (Rc_arg_pos hc0 hc1 g)

lemma Rc_sq {c : ℝ} (hc0 : 0 < c) (hc1 : c < 1) (g : ℝ) :
    Rc c g ^ 2 = 1 + g ^ 2 - 2 * (2 * c ^ 2 - 1) * g :=
  Real.sq_sqrt (Rc_arg_pos hc0 hc1 g).le

lemma hasDerivAt_Rc {c : ℝ} (hc0 : 0 < c) (hc1 : c < 1) (g : ℝ) :
    HasDerivAt (Rc c) ((g - (2 * c ^ 2 - 1)) / Rc c g) g := by
  have h1 : HasDerivAt (fun z : ℝ => 1 + z ^ 2 - 2 * (2 * c ^ 2 - 1) * z)
      (2 * g - 2 * (2 * c ^ 2 - 1)) g := by
    have := ((hasDerivAt_pow 2 g).const_add 1).sub ((hasDerivAt_id' g).const_mul (2 * (2 * c ^ 2 - 1)))
    convert this using 1; ring
  have h2 := h1.sqrt (Rc_arg_pos hc0 hc1 g).ne'
  have e : Rc c = fun z => Real.sqrt (1 + z ^ 2 - 2 * (2 * c ^ 2 - 1) * z) := rfl
  rw [e]
  convert h2 using 1
  have := Rc_pos hc0 hc1 g
  unfold Rc at this
  field_simp

lemma hasDerivAt_R1 {c : ℝ} (hc0 : 0 < c) (hc1 : c < 1) (g : ℝ) :
    HasDerivAt (fun z => (z - (2 * c ^ 2 - 1)) / Rc c z)
      ((1 - (2 * c ^ 2 - 1) ^ 2) / Rc c g ^ 3) g := by
  have hR := Rc_pos hc0 hc1 g
  have := ((hasDerivAt_id' g).sub_const (2 * c ^ 2 - 1)).div (hasDerivAt_Rc hc0 hc1 g) hR.ne'
  convert this using 1
  have hs := Rc_sq hc0 hc1 g
  field_simp
  rw [hs]; ring

lemma hasSum_G1 {c g : ℝ} (hc0 : 0 < c) (hc1 : c < 1) (hg : 0 < g) :
    HasSum (fun N => wt c N * rho1 N g) (1 - (g - (2 * c ^ 2 - 1)) / Rc c g) := by
  have hsum := summable_bd hc0 hc1
  have hd := hasDerivAt_tsum_of_isPreconnected (t := Ioi (0:ℝ)) (y₀ := g) (y := g) hsum isOpen_Ioi
    isPreconnected_Ioi (g := fun N z => wt c N * rho N z) (g' := fun N z => wt c N * rho1 N z)
    (fun N z hz => (hasDerivAt_rho N hz).const_mul _) (fun N z hz => abs_wt_rho1_le c N hz)
    hg (hasSum_G hc0 hc1 hg).summable hg
  have hev : (fun z => ∑' N, wt c N * rho N z) =ᶠ[𝓝 g] fun z => (1 + z) - Rc c z := by
    filter_upwards [Ioi_mem_nhds hg] with z hz using (hasSum_G hc0 hc1 hz).tsum_eq
  have h2 : HasDerivAt (fun z => (1 + z) - Rc c z) (1 - (g - (2 * c ^ 2 - 1)) / Rc c g) g := by
    have := ((hasDerivAt_id' g).const_add 1).sub (hasDerivAt_Rc hc0 hc1 g)
    convert this using 1
  have hu := (hd.congr_of_eventuallyEq hev.symm).unique h2
  have hs : Summable (fun N => wt c N * rho1 N g) :=
    Summable.of_norm_bounded hsum (fun N => abs_wt_rho1_le c N hg)
  rw [← hu]; exact hs.hasSum

lemma hasSum_G2 {c g : ℝ} (hc0 : 0 < c) (hc1 : c < 1) (hg : 0 < g) :
    HasSum (fun N => wt c N * rho2 N g) (-((1 - (2 * c ^ 2 - 1) ^ 2) / Rc c g ^ 3)) := by
  have hsum := summable_bd hc0 hc1
  have hd := hasDerivAt_tsum_of_isPreconnected (t := Ioi (0:ℝ)) (y₀ := g) (y := g) hsum isOpen_Ioi
    isPreconnected_Ioi (g := fun N z => wt c N * rho1 N z) (g' := fun N z => wt c N * rho2 N z)
    (fun N z hz => (hasDerivAt_rho1 N hz).const_mul _) (fun N z hz => abs_wt_rho2_le c N hz)
    hg (hasSum_G1 hc0 hc1 hg).summable hg
  have hev : (fun z => ∑' N, wt c N * rho1 N z) =ᶠ[𝓝 g]
      fun z => 1 - (z - (2 * c ^ 2 - 1)) / Rc c z := by
    filter_upwards [Ioi_mem_nhds hg] with z hz using (hasSum_G1 hc0 hc1 hz).tsum_eq
  have h2 : HasDerivAt (fun z => 1 - (z - (2 * c ^ 2 - 1)) / Rc c z)
      (-((1 - (2 * c ^ 2 - 1) ^ 2) / Rc c g ^ 3)) g := by
    have := (hasDerivAt_R1 hc0 hc1 g).const_sub 1
    convert this using 1
  have hu := (hd.congr_of_eventuallyEq hev.symm).unique h2
  have hs : Summable (fun N => wt c N * rho2 N g) :=
    Summable.of_norm_bounded hsum (fun N => abs_wt_rho2_le c N hg)
  rw [← hu]; exact hs.hasSum

/-- the linear form of `M2integrand` in `(R, R', R'')`. -/
noncomputable def LL (ρ ρ1 ρ2 g : ℝ) : ℝ :=
  (ρ - g * ρ1 + g ^ 2 * ρ2) / (2 * g * (1 + g)) * log g ^ 2
    + (g * ρ1 - ρ) / (g * (1 + g)) * log g * log (1 + g)
    + ρ / (2 * g * (1 + g)) * log (1 + g) ^ 2

lemma LL_lin (ρ ρ1 ρ2 g : ℝ) : LL ρ ρ1 ρ2 g =
    (log g ^ 2 / (2 * g * (1 + g)) - log g * log (1 + g) / (g * (1 + g))
      + log (1 + g) ^ 2 / (2 * g * (1 + g))) * ρ
    + (-(g * log g ^ 2) / (2 * g * (1 + g)) + g * log g * log (1 + g) / (g * (1 + g))) * ρ1
    + (g ^ 2 * log g ^ 2 / (2 * g * (1 + g))) * ρ2 := by
  unfold LL; ring

lemma Lam_eq (N : ℕ) (g : ℝ) : Lam N g = LL (rho N g) (rho1 N g) (rho2 N g) g := rfl

/-- the `c`-independent part `Z₀ = LL(1+g, 1, 0) - L²/(2g)`. -/
noncomputable def Z0 (g : ℝ) : ℝ := LL (1 + g) 1 0 g - log g ^ 2 / (2 * g)

lemma M2_eq {α c g : ℝ} (hc0 : 0 < c) (hc1 : c < 1) (hg : 0 < g) (hk : Real.cos α = 2 * c ^ 2 - 1) :
    M2integrand α g = LL (Rc c g) ((g - (2 * c ^ 2 - 1)) / Rc c g)
      ((1 - (2 * c ^ 2 - 1) ^ 2) / Rc c g ^ 3) g - log g ^ 2 / (2 * g) := by
  simp only [M2integrand, hk, LL]
  rw [← Rc]
  have : (1 + g) ≠ 0 := by positivity
  have : g ≠ 0 := hg.ne'
  have hR := (Rc_pos hc0 hc1 g).ne'
  generalize Rc c g = R at hR ⊢
  field_simp
  ring

lemma hasSum_Lam {α c g : ℝ} (hc0 : 0 < c) (hc1 : c < 1) (hg : 0 < g)
    (hk : Real.cos α = 2 * c ^ 2 - 1) :
    HasSum (fun N => wt c N * Lam N g) (Z0 g - M2integrand α g) := by
  have h := (((hasSum_G hc0 hc1 hg).mul_left (log g ^ 2 / (2 * g * (1 + g))
      - log g * log (1 + g) / (g * (1 + g)) + log (1 + g) ^ 2 / (2 * g * (1 + g)))).add
    ((hasSum_G1 hc0 hc1 hg).mul_left (-(g * log g ^ 2) / (2 * g * (1 + g))
      + g * log g * log (1 + g) / (g * (1 + g))))).add
    ((hasSum_G2 hc0 hc1 hg).mul_left (g ^ 2 * log g ^ 2 / (2 * g * (1 + g))))
  convert h using 1
  · funext N; rw [Lam_eq, LL_lin]; ring
  · rw [Z0, M2_eq hc0 hc1 hg hk, LL_lin, LL_lin]; ring

lemma abs_rho_le (N : ℕ) {g : ℝ} (hg : 0 < g) : |rho N g| ≤ g / 4 ^ N := by
  have e : rho N g = (g / (1 + g) ^ 2) ^ N * (g / (1 + g)) := by
    rw [rho_eq_q N hg, pow_succ]; field_simp
  rw [e, abs_mul, abs_of_nonneg (pow_nonneg (q_nonneg hg) N), abs_of_pos (by positivity)]
  have h1 := pow_q_le N hg
  have h2 : g / (1 + g) ≤ g := div_le_self hg.le (by linarith)
  calc (g / (1 + g) ^ 2) ^ N * (g / (1 + g)) ≤ (1 / 4 ^ N) * g := by
        gcongr
    _ = g / 4 ^ N := by ring

lemma bd_A {L m g ρ P : ℝ} (hg : 0 < g) (hm0 : 0 ≤ m) (hm1 : m ≤ 1) (hP : 0 ≤ P)
    (hρ : |ρ| ≤ g * P) : |(L - m) ^ 2 / (2 * g * (1 + g)) * ρ| ≤ (2 * L ^ 2 + 2) * P := by
  rw [abs_mul, abs_div, abs_of_nonneg (sq_nonneg _), abs_of_pos (by positivity)]
  have h1 : (L - m) ^ 2 / (2 * g * (1 + g)) * |ρ| ≤ (L - m) ^ 2 / (2 * g * (1 + g)) * (g * P) := by
    gcongr
  have h2 : (L - m) ^ 2 / (2 * g * (1 + g)) * (g * P) = (L - m) ^ 2 / (2 * (1 + g)) * P := by
    field_simp
  have h3 : (L - m) ^ 2 / (2 * (1 + g)) ≤ (L - m) ^ 2 := div_le_self (sq_nonneg _) (by linarith)
  have h4 : (L - m) ^ 2 ≤ 2 * L ^ 2 + 2 := by nlinarith [sq_nonneg (L + m)]
  have h5 : (L - m) ^ 2 / (2 * (1 + g)) * P ≤ (2 * L ^ 2 + 2) * P :=
    mul_le_mul_of_nonneg_right (h3.trans h4) hP
  linarith

lemma bd_B0 {L m : ℝ} (hm0 : 0 ≤ m) (hm1 : m ≤ 1) : |2 * L * m - L ^ 2| ≤ 2 * L ^ 2 + 1 := by
  have h1 : |L| * m ≤ |L| := by nlinarith [abs_nonneg L]
  have h2 : |L| ≤ (L ^ 2 + 1) / 2 := by nlinarith [sq_abs L, sq_nonneg (|L| - 1)]
  have h3 : |2 * L * m - L ^ 2| ≤ |2 * L * m| + |L ^ 2| := abs_sub _ _
  have h4 : |2 * L * m| = 2 * (|L| * m) := by
    rw [abs_mul, abs_mul, abs_of_nonneg hm0]; norm_num; ring
  rw [h4, abs_of_nonneg (sq_nonneg L)] at h3
  linarith

lemma bd_B {L m g ρ Q : ℝ} (hg : 0 < g) (hm0 : 0 ≤ m) (hm1 : m ≤ 1)
    (hρ : |ρ| ≤ Q) : |(2 * L * m - L ^ 2) / (2 * (1 + g)) * ρ| ≤ (2 * L ^ 2 + 1) * Q := by
  rw [abs_mul, abs_div, abs_of_pos (by positivity : (0:ℝ) < 2 * (1 + g))]
  have hb := bd_B0 (L := L) hm0 hm1
  have h1 : |2 * L * m - L ^ 2| / (2 * (1 + g)) ≤ |2 * L * m - L ^ 2| :=
    div_le_self (abs_nonneg _) (by linarith)
  have h2 : |2 * L * m - L ^ 2| / (2 * (1 + g)) * |ρ| ≤ |2 * L * m - L ^ 2| * |ρ| :=
    mul_le_mul_of_nonneg_right h1 (abs_nonneg _)
  have h3 : |2 * L * m - L ^ 2| * |ρ| ≤ (2 * L ^ 2 + 1) * Q :=
    mul_le_mul hb hρ (abs_nonneg _) (by positivity)
  linarith

lemma bd_C {L g ρ Q : ℝ} (hg : 0 < g) (hρ : |ρ| ≤ Q) :
    |g * L ^ 2 / (2 * (1 + g)) * ρ| ≤ L ^ 2 * Q := by
  rw [abs_mul, abs_div, abs_mul, abs_of_pos hg, abs_of_nonneg (sq_nonneg L),
    abs_of_pos (by positivity : (0:ℝ) < 2 * (1 + g))]
  have h1 : g * L ^ 2 / (2 * (1 + g)) ≤ L ^ 2 := by
    rw [div_le_iff₀ (by positivity)]; nlinarith [sq_nonneg L]
  have h2 := mul_le_mul h1 hρ (abs_nonneg _) (sq_nonneg L)
  linarith

lemma bd_sum {L x P : ℝ} (hx : 1 ≤ x) (hP : 0 ≤ P) :
    (2 * L ^ 2 + 2) * P + (2 * L ^ 2 + 1) * (x * P) + L ^ 2 * (8 * x ^ 2 * P)
      ≤ 12 * x ^ 2 * P * (L ^ 2 + 1) := by
  have hL := sq_nonneg L
  have a : 0 ≤ 2 * L ^ 2 * (2 * x + 1) * (x - 1) + (12 * x ^ 2 - x - 2) := by
    have := mul_nonneg (mul_nonneg (mul_nonneg (by norm_num : (0:ℝ) ≤ 2) hL)
      (by linarith : 0 ≤ 2 * x + 1)) (by linarith : 0 ≤ x - 1)
    nlinarith
  nlinarith [mul_nonneg hP a]

lemma abs_Lam_le (N : ℕ) {g : ℝ} (hg : 0 < g) (hg1 : g ≤ 1) :
    |Lam N g| ≤ 12 * ((N:ℝ) + 1) ^ 2 / 4 ^ N * (log g ^ 2 + 1) := by
  have h1p : (0:ℝ) < 1 + g := by linarith
  have hm0 : 0 ≤ log (1 + g) := log_nonneg (by linarith)
  have hm1 : log (1 + g) ≤ 1 := by have := log_le_sub_one_of_pos h1p; linarith
  have eA : log g ^ 2 / (2 * g * (1 + g)) - log g * log (1 + g) / (g * (1 + g))
      + log (1 + g) ^ 2 / (2 * g * (1 + g)) = (log g - log (1 + g)) ^ 2 / (2 * g * (1 + g)) := by
    field_simp; ring
  have eB : -(g * log g ^ 2) / (2 * g * (1 + g)) + g * log g * log (1 + g) / (g * (1 + g))
      = (2 * log g * log (1 + g) - log g ^ 2) / (2 * (1 + g)) := by field_simp; ring
  have eC : g ^ 2 * log g ^ 2 / (2 * g * (1 + g)) = g * log g ^ 2 / (2 * (1 + g)) := by
    field_simp
  rw [Lam_eq, LL_lin, eA, eB, eC]
  have hP0 : (0:ℝ) ≤ 1 / 4 ^ N := by positivity
  have r0 : |rho N g| ≤ g * (1 / 4 ^ N) := by have := abs_rho_le N hg; rwa [← mul_one_div] at this
  have r1 : |rho1 N g| ≤ ((N:ℝ) + 1) * (1 / 4 ^ N) := by
    have := abs_rho1_le N hg; rwa [← mul_one_div] at this
  have r2 : |rho2 N g| ≤ 8 * ((N:ℝ) + 1) ^ 2 * (1 / 4 ^ N) := by
    have := abs_rho2_le N hg; rwa [← mul_one_div] at this
  have tA := bd_A (L := log g) hg hm0 hm1 hP0 r0
  have tB := bd_B (L := log g) hg hm0 hm1 r1
  have tC := bd_C (L := log g) hg r2
  have hs := bd_sum (L := log g) (x := (N:ℝ) + 1) (by linarith [N.cast_nonneg (α := ℝ)]) hP0
  have h3 := abs_add_three ((log g - log (1 + g)) ^ 2 / (2 * g * (1 + g)) * rho N g)
    ((2 * log g * log (1 + g) - log g ^ 2) / (2 * (1 + g)) * rho1 N g)
    (g * log g ^ 2 / (2 * (1 + g)) * rho2 N g)
  have e : 12 * ((N:ℝ) + 1) ^ 2 / 4 ^ N * (log g ^ 2 + 1)
      = 12 * ((N:ℝ) + 1) ^ 2 * (1 / 4 ^ N) * (log g ^ 2 + 1) := by ring
  rw [e]
  linarith

lemma abs_wt_Lam_le {c : ℝ} (N : ℕ) {g : ℝ} (hg : 0 < g) (hg1 : g ≤ 1) :
    ‖wt c N * Lam N g‖ ≤ 16 * ((N:ℝ) + 1) ^ 2 * c ^ (2 * N + 2) * (3 / 2 * (log g ^ 2 + 1)) := by
  rw [Real.norm_eq_abs, abs_mul, abs_of_nonneg (wt_nonneg c N)]
  have h1 := abs_Lam_le N hg hg1
  have h2 := wt_le c N
  have hc : 0 ≤ c ^ (2 * N + 2) := by rw [show 2 * N + 2 = 2 * (N + 1) by ring, pow_mul]; positivity
  have hK : 0 ≤ log g ^ 2 + 1 := by positivity
  calc wt c N * |Lam N g| ≤ (4 ^ (N + 1) / 2 * c ^ (2 * N + 2)) *
        (12 * ((N:ℝ) + 1) ^ 2 / 4 ^ N * (log g ^ 2 + 1)) := by
        gcongr
    _ = _ := by rw [pow_succ]; field_simp; ring

lemma ii_logsq : IntervalIntegrable (fun t : ℝ => log t ^ 2) volume 0 1 := by
  have := ii_Kb 2 0 0 0; simpa using this

lemma meas_Lam (N : ℕ) : Measurable (Lam N) := by
  unfold Lam rho1 rho2 rho; fun_prop

lemma hasSum_intLam {α c : ℝ} (hc0 : 0 < c) (hc1 : c < 1) (hk : Real.cos α = 2 * c ^ 2 - 1) :
    HasSum (fun N => wt c N * ∫ g in (0:ℝ)..1, Lam N g)
      (∫ g in (0:ℝ)..1, (Z0 g - M2integrand α g)) := by
  have hs := summable_bd hc0 hc1
  have e : (fun N => wt c N * ∫ g in (0:ℝ)..1, Lam N g)
      = fun N => ∫ g in (0:ℝ)..1, wt c N * Lam N g := by
    funext N; rw [intervalIntegral.integral_const_mul]
  rw [e]
  refine intervalIntegral.hasSum_integral_of_dominated_convergence
    (fun N g => 16 * ((N:ℝ) + 1) ^ 2 * c ^ (2 * N + 2) * (3 / 2 * (log g ^ 2 + 1)))
    (fun N => ((meas_Lam N).const_mul _).aestronglyMeasurable) (fun N => ?_) ?_ ?_ ?_
  · refine Filter.Eventually.of_forall fun g hg => ?_
    rw [uIoc_of_le zero_le_one] at hg
    exact abs_wt_Lam_le N hg.1 hg.2
  · exact Filter.Eventually.of_forall fun g _ => hs.mul_right _
  · have e2 : (fun g : ℝ => ∑' N : ℕ, 16 * ((N:ℝ) + 1) ^ 2 * c ^ (2 * N + 2) * (3 / 2 * (log g ^ 2 + 1)))
        = fun g => (∑' N : ℕ, 16 * ((N:ℝ) + 1) ^ 2 * c ^ (2 * N + 2)) * (3 / 2 * (log g ^ 2 + 1)) := by
      funext g; rw [tsum_mul_right]
    rw [e2]
    exact ((ii_logsq.add intervalIntegrable_const).const_mul (3 / 2)).const_mul _
  · refine Filter.Eventually.of_forall fun g hg => ?_
    rw [uIoc_of_le zero_le_one] at hg
    exact hasSum_Lam hc0 hc1 hg.1 hk

lemma Z0_eq {g : ℝ} (hg : 0 < g) : Z0 g = log (1 + g) ^ 2 / (2 * g)
    - log g * log (1 + g) / (g * (1 + g)) - log g ^ 2 / (2 * (1 + g)) := by
  have : (1 + g) ≠ 0 := by positivity
  have : g ≠ 0 := hg.ne'
  unfold Z0 LL; field_simp; ring

/-- primitive of `Z₀`. -/
noncomputable def P0 (g : ℝ) : ℝ := (log g * log (1 + g) ^ 2 - log g ^ 2 * log (1 + g)) / 2

lemma hasDerivAt_P0 {g : ℝ} (hg : 0 < g) : HasDerivAt P0 (Z0 g) g := by
  have h1p : (0:ℝ) < 1 + g := by linarith
  have hL := hasDerivAt_log hg.ne'
  have hm : HasDerivAt (fun z => log (1 + z)) (1 / (1 + g)) g := by
    have := ((hasDerivAt_id' g).const_add 1).log h1p.ne'
    simpa using this
  have := ((hL.mul (hm.pow 2)).sub ((hL.pow 2).mul hm)).div_const 2
  have e : P0 = fun z => (log z * log (1 + z) ^ 2 - log z ^ 2 * log (1 + z)) / 2 := rfl
  rw [e]
  convert this using 1
  rw [Z0_eq hg]
  simp only [Pi.pow_apply]
  field_simp
  ring

lemma P0_cont : ContinuousOn P0 (Icc 0 1) := by
  intro x hx
  rcases eq_or_lt_of_le hx.1 with h | h
  · subst h
    apply ContinuousWithinAt.mono _ Icc_subset_Ici_self
    rw [← continuousWithinAt_Ioi_iff_Ici, ContinuousWithinAt]
    have hP0 : P0 0 = 0 := by simp [P0]
    rw [hP0]
    have t1 := bl_t_log_pow_tendsto 1 one_ne_zero
    have t2 := bl_t_log_pow_tendsto 2 two_ne_zero
    have tb : Tendsto (fun t : ℝ => (t * |t * log t ^ 1| + |t * log t ^ 2|) / 2) (𝓝[>] 0) (𝓝 0) := by
      have := ((((tendsto_nhdsWithin_of_tendsto_nhds (continuous_id.tendsto 0)).mul
        t1.abs).add t2.abs).div_const 2)
      simpa using this
    refine squeeze_zero_norm' ?_ tb
    filter_upwards [Ioo_mem_nhdsGT (show (0:ℝ) < 1 by norm_num)] with t ht
    have ht0 : 0 < t := ht.1
    have hm0 : 0 ≤ log (1 + t) := log_nonneg (by linarith)
    have hm1 : log (1 + t) ≤ t := by have := log_le_sub_one_of_pos (by linarith : (0:ℝ) < 1 + t); linarith
    rw [Real.norm_eq_abs, P0, abs_div, abs_of_pos (by norm_num : (0:ℝ) < 2)]
    gcongr
    calc |log t * log (1 + t) ^ 2 - log t ^ 2 * log (1 + t)|
        ≤ |log t * log (1 + t) ^ 2| + |log t ^ 2 * log (1 + t)| := abs_sub _ _
      _ ≤ t * |t * log t ^ 1| + |t * log t ^ 2| := by
          simp only [abs_mul, abs_pow, abs_of_nonneg hm0, abs_of_pos ht0, pow_one, sq_abs]
          have a1 : log (1 + t) ^ 2 ≤ t ^ 2 := pow_le_pow_left₀ hm0 hm1 2
          have a2 : |log t| * log (1 + t) ^ 2 ≤ |log t| * t ^ 2 := by gcongr
          have a3 : log t ^ 2 * log (1 + t) ≤ log t ^ 2 * t := by gcongr
          nlinarith [abs_nonneg (log t)]
  · have hx0 : x ≠ 0 := h.ne'
    exact (hasDerivAt_P0 h).continuousAt.continuousWithinAt

lemma meas_Z0 : Measurable Z0 := by unfold Z0 LL; fun_prop

lemma ii_Z0 : IntervalIntegrable Z0 volume 0 1 := by
  have hb : IntervalIntegrable (fun t : ℝ => log t ^ 2 + 1) volume 0 1 :=
    ii_logsq.add intervalIntegrable_const
  refine IntervalIntegrable.mono_fun' hb meas_Z0.aestronglyMeasurable ?_
  have hae : ∀ᵐ t ∂(volume.restrict (Ι (0:ℝ) 1)), t ∈ Ioc (0:ℝ) 1 := by
    rw [uIoc_of_le zero_le_one]; exact ae_restrict_mem measurableSet_Ioc
  filter_upwards [hae] with t ht
  have ht0 : 0 < t := ht.1
  have h1p : (0:ℝ) < 1 + t := by linarith
  have hm0 : 0 ≤ log (1 + t) := log_nonneg (by linarith)
  have hm1 : log (1 + t) ≤ t := by have := log_le_sub_one_of_pos h1p; linarith
  rw [Real.norm_eq_abs, Z0_eq ht0]
  have a1 : log (1 + t) ^ 2 / (2 * t) ≤ 1 / 2 := by
    rw [div_le_iff₀ (by positivity)]
    have := pow_le_pow_left₀ hm0 hm1 2
    have : t ^ 2 ≤ t := by nlinarith [ht.2]
    nlinarith
  have a1' : 0 ≤ log (1 + t) ^ 2 / (2 * t) := by positivity
  have a2 : |log t * log (1 + t) / (t * (1 + t))| ≤ |log t| := by
    rw [abs_div, abs_mul, abs_of_nonneg hm0, abs_of_pos (by positivity : (0:ℝ) < t * (1 + t))]
    rw [div_le_iff₀ (by positivity)]
    have : log (1 + t) ≤ t * (1 + t) := by nlinarith
    calc |log t| * log (1 + t) ≤ |log t| * (t * (1 + t)) := by gcongr
      _ = _ := by ring
  have a3 : 0 ≤ log t ^ 2 / (2 * (1 + t)) := by positivity
  have a4 : log t ^ 2 / (2 * (1 + t)) ≤ log t ^ 2 / 2 := by
    gcongr; linarith
  have a5 : |log t| ≤ (log t ^ 2 + 1) / 2 := by nlinarith [sq_abs (log t), sq_nonneg (|log t| - 1)]
  have := abs_sub (log (1 + t) ^ 2 / (2 * t) - log t * log (1 + t) / (t * (1 + t)))
    (log t ^ 2 / (2 * (1 + t)))
  have := abs_sub (log (1 + t) ^ 2 / (2 * t)) (log t * log (1 + t) / (t * (1 + t)))
  rw [abs_of_nonneg a1', abs_of_nonneg a3] at *
  linarith

lemma int_Z0 : ∫ g in (0:ℝ)..1, Z0 g = 0 := by
  rw [intervalIntegral.integral_eq_sub_of_hasDerivAt_of_le zero_le_one P0_cont
    (fun x hx => hasDerivAt_P0 hx.1) ii_Z0]
  simp [P0]

/-- `thm:flat` at `k = 2`, with the coefficient `M_2 = ∫_0^1 M2integrand` of `stmt_7`:
`∫_0^1 M2integrand α = Σ_{n≥1} (r_n + ζ(2)/(4n²-1)) c^{2n}`. -/
lemma hasSum_M2 {α c : ℝ} (hc0 : 0 < c) (hc1 : c < 1) (hk : Real.cos α = 2 * c ^ 2 - 1)
    (hM : IntervalIntegrable (M2integrand α) volume 0 1) :
    HasSum (fun N : ℕ => (r (N + 1) + (π ^ 2 / 6) / (4 * ((N:ℝ) + 1) ^ 2 - 1)) * c ^ (2 * (N + 1)))
      (∫ g in (0:ℝ)..1, M2integrand α g) := by
  have h := (hasSum_intLam hc0 hc1 hk).neg
  rw [intervalIntegral.integral_sub ii_Z0 hM, int_Z0, zero_sub, neg_neg] at h
  convert h using 1
  funext N
  rw [int_Lam, coef_main, wt, kap, show 2 * (N + 1) = 2 * N + 2 by ring]
  have hB := (Kb00_pos N N).ne'
  have : (N:ℝ) + 1 ≠ 0 := by positivity
  have : 2 * (N:ℝ) + 1 ≠ 0 := by positivity
  field_simp

end Z07

theorem Z07.proof_8 : Z07.stmt_8 := by
  intro α h0 h1
  have hs : 0 < Real.sin (α / 2) := Real.sin_pos_of_pos_of_lt_pi (by linarith) (by linarith)
  have hc0 : 0 < D04.c α := by
    unfold D04.c
    exact Real.cos_pos_of_mem_Ioo ⟨by linarith [Real.pi_pos], by linarith⟩
  have hc1 : D04.c α < 1 := by
    have h2 := Real.sin_sq_add_cos_sq (α / 2)
    unfold D04.c
    nlinarith [Real.cos_le_one (α / 2)]
  have hk : Real.cos α = 2 * D04.c α ^ 2 - 1 := by
    unfold D04.c
    rw [← Real.cos_two_mul]; ring_nf
  obtain ⟨hI, hC⟩ := Z07.proof_7 α h0 h1
  exact ⟨_, Z07.hasSum_M2 hc0 hc1 hk
    ((intervalIntegrable_iff_integrableOn_Ioc_of_le zero_le_one).mpr hI), hC⟩

theorem Z07.proof_9 : Z07.stmt_9 := fun _ h0 h1 => Z07.hasSum_r h0 h1

theorem Z07.proof_10 : Z07.stmt_10 := Z07.stmt10_of Z07.proof_8 Z07.proof_9
