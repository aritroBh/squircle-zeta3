import Mathlib
import Zeta3.Wave8.Q01_PureWeight
import Zeta3.Wave5.Y05_Octahedron
import Zeta3.Wave1.T04_SmallFacts
/-!
# R03 — pure weight in every dimension `N` (`thm:pw`, last sentence) and the weight claims
of `thm:oct`

Target id: R03.
Paper labels: `thm:pw`, last sentence (paper/ieee/main.tex lines 1183–1191, `τ = 1 - 1/p` from
line 1181); `thm:oct`, the two weight claims (main.tex lines 1973–1990, the sentence on `N`
dimensions line 1992). Detailed proofs: paper/ieee/proofs/sec_pureweight.tex, subsection
"Dimension $n$" (lines 157–244: `eq:pw:Sn`, `lem:pw:dirichlet`, the centring and
`eq:pw:factor`, the flat-face factor, `lem:pw:X`, `lem:pw:moments`, `prop:pw:n`; the setting
bullet lines 25–28) and paper/ieee/proofs/sec_g7_octa.tex (`g7:lem:both` lines 43–62,
`g7:lem:X2` 64–75, `g7:lem:struct` 95–104, `g7:lem:w23` 112–132, `g7:lem:zeta3n` 134–140,
`g7:lem:dcoef` 142–152).

Earlier files used (imported, not restated): `Q01.IsZetaPoly` (the meaning of weight,
sec_pureweight.tex line 39), `Q01` (the two-dimensional theorem), `Y05.S` (the face integral
`eq:pw:Sn` for `n = 3`, with its integrability, the Dirichlet–Gamma factorisation for `n = 3`
and the explicit series through `τ⁶`), `T04.Y` (`Y_N`), `T04.coeffFormula`
(`-(N²-1)((N-1)³-2)/(3N²)`), `X01.P` (the perimeter integral of the rhombus, `A = B = 1` is
the squircle), `X05.z n = (riemannZeta n).re` (the real number `ζ(n)`; all zeta values in this
file are `X05.z`, i.e. real parts of Mathlib's `riemannZeta` at natural arguments).

Quoted verbatim.

main.tex 1190 (`thm:pw`, last sentence):

  The same holds for every rhombus (with polynomials in $q^2$ as coefficients) and, after
  division by the rational factor
  $Y_N(\tau)=(1-\tau)^{N-1}\big/\prod_{i=1}^{N-1}\bigl(1-\frac{(N-1)\tau}{i}\bigr)$, so
  $Y_2=1$ and $Y_3=1/(2-p)$, for the surface area of $|x_1|^p+\dots+|x_N|^p=1$ in every
  dimension $N$.

main.tex 1984 and 1989 (`thm:oct`, the two weight claims):

  and for every $k$, $[d^k]S$ is $4\sqrt3$ times a polynomial with rational coefficients in
  $\zeta(2),\dots,\zeta(k)$ of weight at most $k$.

  where the coefficient of $\tau^k$ is a polynomial in $\zeta(2),\dots,\zeta(k)$, homogeneous
  of weight $k$.

(`thm:oct` concerns `\log\frac{(2-p)\,S(p)}{4\sqrt3}` and `S(1+d)`, `d = p - 1`, "as $d\to0$
from either side".)

sec_pureweight.tex 26 (setting):

  \emph{Dimension $n\ge2$.} $S_n(p)$ is the $(n-1)$-dimensional area of
  $\{x\in\RR^n:\sum|x_i|^p=1\}$, and $S_n(1)=2^n\sqrt n/(n-1)!$ is the area of the
  cross-polytope. ... We prove $\log\frac{S_n(p)}{S_n(1)Y_n(\tau)}=\sum_{k\ge2}w_k(n)\tau^k$
  with $w_k(n)\in\QQ[\zeta(2),\dots,\zeta(k)]$ of weight $k$. For $n=2$, $S_2=P$,
  $S_2(1)=4\sqrt2$ and $Y_2=\frac{1-\tau}{1-\tau}=1$.

sec_pureweight.tex 163–166 (`eq:pw:Sn`):

  S_n(p)=2^n\varepsilon^{n-1}\int_\Delta\Bigl(\prod_{i=1}^nt_i\Bigr)^{-\tau}
  \Bigl(\sum_{i=1}^nt_i^{2\tau}\Bigr)^{1/2}dt_1\cdots dt_{n-1},
  $\Delta=\{t_i\ge0,\ \sum_{i<n}t_i\le1\}$. The integral converges since $\tau<1$.

sec_pureweight.tex 168–172 (`lem:pw:dirichlet`):

  Let $a>-1$ and let $F\ge0$ be measurable on $(0,\infty)^n$ and homogeneous of degree $0$.
  Then \int_{(0,\infty)^n}\prod g_i^{\,a}\,\eu^{-\sum g_i}F(g)\,dg=\Gamma(n(1+a))
  \int_\Delta\Bigl(\prod t_i\Bigr)^aF(t)\,dt_1\cdots dt_{n-1}.

sec_pureweight.tex 182–183 (centring):

  S_n(p)=2^n\sqrt n\,\varepsilon^{n-1}\frac{\Gamma(1+a)^n}{\Gamma(n-(n-1)\tau)}\,B_n(\tau),
  \qquad B_n(\tau)=\E\bigl[f(\tau,m(G))\bigr],

sec_pureweight.tex 197 (flat-face factor):

  So $\Gamma(1+a)^n/\Gamma(1-(n-1)\tau)=\exp\bigl(\sum_{k\ge2}\lambda_k\zeta(k)\tau^k\bigr)$
  with $\lambda_k=-\frac{(n-1)^k(1-n^{1-k})}k\in\QQ$.

sec_pureweight.tex 202–203 (`lem:pw:X`):

  For each $J\ge0$ there are $C_J,\tau_J>0$ with $\E|X|^{J+1}\le C_J\tau^{J+1}$ and
  $\bigl|B_n-\sum_{j\le J}\binom{1/2}j\E[X^j]\bigr|\le C_J\,\tau^{J+1}$ for
  $0<\tau\le\tau_J$.

sec_pureweight.tex 212–213 (`lem:pw:moments`):

  For each $j\ge0$ there is a pure series $\mathbf E_j\in\mathcal R$ with $\mathbf E_0=1$ whose
  evaluation at $Z_k=\zeta(k)$ is the Taylor series of $\E[X^j]$ at $\tau=0$. Moreover,
  $\E[X^j]$ is real analytic near $\tau=0$, and $\mathbf E_j\in\tau^j\mathcal R$.

sec_pureweight.tex 234–238 (`prop:pw:n`):

  Let $n\ge2$. There is a pure series $\sum_{k\ge2}w_k(n)\tau^k\in\mathcal R$ (so
  $w_k(n)\in\QQ[Z_2,\dots,Z_k]$ is homogeneous of weight $k$) such that, for every $K$, as
  $\tau\to0^+$, \log\frac{S_n(p)}{S_n(1)Y_n(\tau)}=\sum_{k=2}^Kw_k(n)(\zeta)\,\tau^k
  +O(\tau^{K+1}).

sec_g7_octa.tex 44–48 (`g7:lem:both`):

  For every $K$ there are $C,\tau_0>0$ such that for $0<|\tau|<\tau_0$ (both signs,
  $p=\frac1{1-\tau}$) \Bigl|\log\frac{S_n(p)}{S_n(1)Y_n(\tau)}-\sum_{k=2}^Kw_k(n)(\zeta)\,
  \tau^k\Bigr|\le C|\tau|^{K+1}.

sec_g7_octa.tex 113–116 (`g7:lem:w23`):

  For every $n\ge2$, w_2(n)=-\frac{(n-1)(n^2-2n-1)}{2n}\,Z_2,\qquad
  w_3(n)=-\frac{(n^2-1)\bigl((n-1)^3-2\bigr)}{3n^2}\,Z_3 .

sec_g7_octa.tex 96–100 (`g7:lem:struct`):

  Then, as $\tau\to0$ through both signs, for every $K$, \log\frac{(2-p)S(p)}{4\sqrt3}=
  \sum_{k=2}^Kw_k(\zeta)\tau^k+O(|\tau|^{K+1}), and for every $k$ the number
  $[d^k]\frac S{4\sqrt3}$ is the value at $Z_j=\zeta(j)$ of a polynomial
  $P_k\in\QQ[Z_2,\dots,Z_k]$ all of whose monomials have weight at most $k$.

## Conventions and reading

* Dimension `N ≥ 2` (the supplement's range "$n\ge2$"). The simplex `Δ` lives in
  `Fin (N-1) → ℝ` with Lebesgue measure, and `t_N = 1 - Σ_{i<N} t_i` is appended by `ext`
  (`Fin.snoc`). `S N p` is the face integral `eq:pw:Sn` with `ε = 1/p`, `τ = 1 - 1/p`, exactly
  as Y05 does for `n = 3`; that it is the Hausdorff area is not claimed (as in Y05).
* The Gamma variables: `aPar N τ = -(N-1)τ/N`, the density `∏ g_i^a e^{-Σ g_i}` on the open
  orthant `(0,∞)^N`, and `gammaE N τ F = Γ(1+a)^{-N} ∫ dens · F` (the expectation of the
  supplement). `m_i = log g_i - (1/N) Σ_j log g_j`, `Y = (1/N) Σ_i e^{2τ m_i}`, `X = Y - 1`,
  `B_N(τ) = E[√Y]`. Bochner integrals of non-integrable functions are `0`; every statement that
  uses `gammaE` therefore also asserts the integrability it needs.
* `lem:pw:dirichlet` is stated with lower Lebesgue integrals (`F ≥ 0`, both sides may be `∞`),
  which is the paper's statement for nonnegative measurable `F`.
* Weight is at the level of real numbers, as in Q01: `Q01.IsZetaPoly k x` says `x` is the value
  of a rational polynomial in `ζ(2), ζ(3), …` homogeneous of weight `k`. "Weight at most `k`"
  (`IsZetaPolyLE k`) is a finite sum of `IsZetaPoly j` numbers, `j ≤ k`. The formal statements
  (`𝐄_j ∈ τ^j𝓡`) are recorded through their values (Taylor coefficients of order `< j` vanish).
* `binom(1/2, j)` is written out as `∏_{i<j}(1/2 - i)/j!` (`binomHalf`).
* `lem:pw:X` is stated for both signs of `τ` (`|τ| ≤ τ_J`), as `g7:lem:both` proves it.
* The expansions of `thm:oct` are stated for `Y05.S`; stmt_3 shows `S 3 = Y05.S`.

## Statements

* stmt_1: `eq:pw:Sn` converges: the face integrand is integrable on `Δ` for every `N ≥ 2`,
  `p > 0`.
* stmt_2: `S_N(1) = 2^N √N/(N-1)!` (the cross-polytope).
* stmt_3: `S_2 = X01.P 1 1` (the squircle perimeter) for `p > 0`, and `S_3 = Y05.S`.
* stmt_4: `lem:pw:dirichlet` for every `n ≥ 2`.
* stmt_5: the centring: `S_N(p) = 2^N √N ε^{N-1} Γ(1+a)^N/Γ(N-(N-1)τ) · B_N(τ)` for `p > 0`.
* stmt_6: the flat-face factor `log(Γ(1+a)^N/Γ(1-(N-1)τ)) = Σ_{k≥2} λ_k ζ(k) τ^k`,
  `|τ| < 1/(N-1)`.
* stmt_7: `lem:pw:X` (both signs).
* stmt_8: `lem:pw:moments`.
* stmt_9: `prop:pw:n` = the last sentence of `thm:pw` (`τ → 0⁺`, every `N ≥ 2`, every order).
* stmt_10: `g7:lem:both` with `g7:lem:w23`: both signs, and `w_2(N)`, `w_3(N)` explicit; the
  `ζ(3)` coefficient is `T04.coeffFormula N`.
* stmt_11: `thm:oct`, the two weight claims for every `k` (`= g7:lem:struct`), the final claim.

## Numerical check (before freezing)

Direct tanh–sinh quadrature of `eq:pw:Sn` (Duffy map of the simplex, double precision) at
16 values `τ ∈ [-0.04, 0.04]`, least-squares fit of `c_2 τ² + … + c_8 τ⁸` to
`log(S_N/(S_N(1)Y_N))`: `N=2`: `c_2 = 0.41123351627`, `c_3 = 0.300514026` (formula `ζ(2)/4 =
0.41123351671`, `ζ(3)/4 = 0.300514226`); `N=3`: `-1.09662271123`, `-2.13698991` (formula
`-2ζ(2)/3 = -1.09662271123`, `-16ζ(3)/9 = -2.13699005`); `N=4`: `-4.31795192`, `-9.3910715`
(formula `-21ζ(2)/8 = -4.31795193`, `-125ζ(3)/16 = -9.3910696`). Also `(N-1)!·∫_Δ 1·√N = √N`
at `τ = 0` for `N = 2,3,4`.

## NOT COVERED

* That the face integral `eq:pw:Sn` is the `(N-1)`-dimensional Hausdorff area of
  `|x_1|^p+…+|x_N|^p = 1` (symmetry, graph parametrisation); as in Y05 the area is defined by
  `eq:pw:Sn`.
* The formal-series level of `lem:pw:pure`/`lem:pw:moments` (`𝐄_j` as elements of
  `ℚ[Z][[τ]]`, uniqueness of the formal polynomial, `ζ(2)² = 5ζ(4)/2` non-uniqueness): weight is
  "has a representation of that weight" for the real number, as in Q01.
* `g7:lem:X2` and `g7:lem:Eformal` (the second-order size of `X`; Y05 has the case `n=3`,
  `J=2`), `g7:lem:zeta3n` as a statement on formal coefficients of `[d^3]` (its numerical
  content `w_3(N) = coeffFormula N · ζ(3)` is in stmt_10), the explicit coefficients through
  `τ⁶`/`d⁴` (Y05), the divergence of the series.
* The rhombus part of `thm:pw` (Q01).
-/
-- BEGIN STATEMENTS
namespace R03

open Real

/-- The simplex `Δ = {t_i ≥ 0, Σ_{i<n} t_i ≤ 1}` in `ℝ^m`, `m = n - 1` (sec_pureweight.tex 166). -/
def simplex (m : ℕ) : Set (Fin m → ℝ) := {t | (∀ i, 0 ≤ t i) ∧ ∑ i, t i ≤ 1}

/-- The point `(t_1, …, t_m, t_{m+1})` with `t_{m+1} = 1 - Σ_{i≤m} t_i`. -/
noncomputable def ext {m : ℕ} (t : Fin m → ℝ) : Fin (m + 1) → ℝ :=
  Fin.snoc (α := fun _ => ℝ) t (1 - ∑ i, t i)

/-- The integrand of `eq:pw:Sn`: `(∏_{i=1}^n t_i)^{-τ} (Σ_{i=1}^n t_i^{2τ})^{1/2}`, `n = m + 1`. -/
noncomputable def faceIntegrand (m : ℕ) (τ : ℝ) (t : Fin m → ℝ) : ℝ :=
  (∏ i, ext t i) ^ (-τ) * Real.sqrt (∑ i, ext t i ^ (2 * τ))

/-- `S_N(p)`, the surface area of `|x_1|^p+…+|x_N|^p = 1`, DEFINED by the paper's face integral
`eq:pw:Sn`: `S_N(p) = 2^N ε^{N-1} ∫_Δ (∏ t_i)^{-τ} (Σ t_i^{2τ})^{1/2}`, `ε = 1/p`,
`τ = 1 - 1/p`. -/
noncomputable def S (N : ℕ) (p : ℝ) : ℝ :=
  2 ^ N * (1 / p) ^ (N - 1) * ∫ t in simplex (N - 1), faceIntegrand (N - 1) (1 - 1 / p) t

/-- The open orthant `(0,∞)^n`. -/
def orth (n : ℕ) : Set (Fin n → ℝ) := Set.pi Set.univ (fun _ => Set.Ioi (0 : ℝ))

/-- `a = -(N-1)τ/N` (sec_pureweight.tex 159). -/
noncomputable def aPar (N : ℕ) (τ : ℝ) : ℝ := -(((N : ℝ) - 1) * τ) / N

/-- The unnormalised Gamma density `∏ g_i^a e^{-Σ g_i}`, `a = aPar N τ`. -/
noncomputable def dens (N : ℕ) (τ : ℝ) (g : Fin N → ℝ) : ℝ :=
  (∏ i, g i ^ aPar N τ) * Real.exp (-∑ i, g i)

/-- `E[F(G)]` for independent `G_1,…,G_N` with density `g^a e^{-g}/Γ(1+a)` on `(0,∞)`. -/
noncomputable def gammaE (N : ℕ) (τ : ℝ) (F : (Fin N → ℝ) → ℝ) : ℝ :=
  (Real.Gamma (1 + aPar N τ) ^ N)⁻¹ * ∫ g in orth N, dens N τ g * F g

/-- `Y = (1/N) Σ_i e^{2τ m_i}` with `m_i = log g_i - (1/N) Σ_j log g_j`. -/
noncomputable def Ymean (N : ℕ) (τ : ℝ) (g : Fin N → ℝ) : ℝ :=
  (∑ i, Real.exp (2 * τ * (Real.log (g i) - (∑ j, Real.log (g j)) / N))) / N

/-- The bending factor `B_N(τ) = E[√Y]`. -/
noncomputable def B (N : ℕ) (τ : ℝ) : ℝ := gammaE N τ (fun g => Real.sqrt (Ymean N τ g))

/-- The moment `E[X^j]`, `X = Y - 1`. -/
noncomputable def momX (N j : ℕ) (τ : ℝ) : ℝ := gammaE N τ (fun g => (Ymean N τ g - 1) ^ j)

/-- `binom(1/2, j) = ∏_{i<j} (1/2 - i) / j!`. -/
noncomputable def binomHalf (j : ℕ) : ℝ :=
  (∏ i ∈ Finset.range j, ((1 : ℝ) / 2 - i)) / (j.factorial : ℝ)

/-- "a polynomial with rational coefficients in `ζ(2),…,ζ(k)` of weight at most `k`"
(sec_g7_octa.tex line 24): a finite sum of numbers of pure weight `j ≤ k`. -/
def IsZetaPolyLE (k : ℕ) (x : ℝ) : Prop :=
  ∃ f : ℕ → ℝ, (∀ j ≤ k, Q01.IsZetaPoly j (f j)) ∧ x = ∑ j ∈ Finset.range (k + 1), f j

/-- sec_pureweight.tex 166 ("The integral converges since $\tau<1$"; both signs in
`g7:lem:both`): for every `N ≥ 2` and `p > 0` the integrand of `eq:pw:Sn` is integrable on `Δ`. -/
def stmt_1 : Prop :=
  ∀ N : ℕ, 2 ≤ N → ∀ p : ℝ, 0 < p →
    MeasureTheory.IntegrableOn (faceIntegrand (N - 1) (1 - 1 / p)) (simplex (N - 1))

/-- sec_pureweight.tex 26 and 183: `S_N(1) = 2^N √N/(N-1)!`, the area of the cross-polytope. -/
def stmt_2 : Prop :=
  ∀ N : ℕ, 2 ≤ N → S N 1 = 2 ^ N * Real.sqrt N / ((N - 1).factorial : ℝ)

/-- sec_pureweight.tex 26 ("For $n=2$, $S_2=P$") and sec_g7_octa.tex 6 ("$S=S_3$"): the face
integral for `N = 2` is the perimeter integral of the squircle `X01.P 1 1` (`p > 0`), and for
`N = 3` it is Y05's `S`. -/
def stmt_3 : Prop :=
  (∀ p : ℝ, 0 < p → S 2 p = X01.P 1 1 p) ∧ (∀ p : ℝ, S 3 p = Y05.S p)

/-- `lem:pw:dirichlet` (sec_pureweight.tex 168–172), `n = m + 1 ≥ 2`: for `a > -1` and `F ≥ 0`
measurable, homogeneous of degree `0` on `(0,∞)^n`,
`∫_{(0,∞)^n} ∏ g_i^a e^{-Σ g_i} F(g) dg = Γ(n(1+a)) ∫_Δ (∏ t_i)^a F(t) dt_1⋯dt_{n-1}`
(lower Lebesgue integrals; `t_n = 1 - Σ_{i<n} t_i`). -/
def stmt_4 : Prop :=
  ∀ m : ℕ, 1 ≤ m → ∀ a : ℝ, -1 < a → ∀ F : (Fin (m + 1) → ℝ) → ℝ, Measurable F →
    (∀ g ∈ orth (m + 1), 0 ≤ F g) →
    (∀ c : ℝ, 0 < c → ∀ g ∈ orth (m + 1), F (c • g) = F g) →
    ∫⁻ g in orth (m + 1), ENNReal.ofReal ((∏ i, g i ^ a) * Real.exp (-∑ i, g i) * F g) =
      ENNReal.ofReal (Real.Gamma (((m : ℝ) + 1) * (1 + a))) *
        ∫⁻ t in simplex m, ENNReal.ofReal ((∏ i, ext t i ^ a) * F (ext t))

/-- The centring (sec_pureweight.tex 175–183), every `N ≥ 2`, every `p > 0` (`τ = 1 - 1/p`,
`a = -(N-1)τ/N`): `√Y` is integrable against the Gamma density and
`S_N(p) = 2^N √N (1/p)^{N-1} Γ(1+a)^N / Γ(N-(N-1)τ) · B_N(τ)`. -/
def stmt_5 : Prop :=
  ∀ N : ℕ, 2 ≤ N → ∀ p : ℝ, 0 < p →
    MeasureTheory.IntegrableOn
        (fun g => dens N (1 - 1 / p) g * Real.sqrt (Ymean N (1 - 1 / p) g)) (orth N) ∧
      S N p = 2 ^ N * Real.sqrt N * (1 / p) ^ (N - 1) *
        Real.Gamma (1 + aPar N (1 - 1 / p)) ^ N /
          Real.Gamma (N - ((N : ℝ) - 1) * (1 - 1 / p)) * B N (1 - 1 / p)

/-- The flat-face factor (sec_pureweight.tex 190–197), every `N ≥ 2`, `|τ| < 1/(N-1)`:
`log(Γ(1+a)^N/Γ(1-(N-1)τ)) = Σ_{k≥2} λ_k ζ(k) τ^k`, `λ_k = -(N-1)^k(1-N^{1-k})/k`
(Euler's `γ` cancels). -/
def stmt_6 : Prop :=
  ∀ N : ℕ, 2 ≤ N → ∀ τ : ℝ, |τ| < 1 / ((N : ℝ) - 1) →
    HasSum (fun k : ℕ => -(((N : ℝ) - 1) ^ (k + 2) * (1 - 1 / (N : ℝ) ^ (k + 1)))
          / ((k + 2 : ℕ) : ℝ) * X05.z (k + 2) * τ ^ (k + 2))
      (Real.log (Real.Gamma (1 + aPar N τ) ^ N / Real.Gamma (1 - ((N : ℝ) - 1) * τ)))

/-- `lem:pw:X` (sec_pureweight.tex 201–208), for both signs of `τ` (`g7:lem:both`): for every
`N ≥ 2`, `J ≥ 0` there are `C`, `τ_J > 0` such that for `|τ| ≤ τ_J` the moments `X^j`,
`j ≤ J+1`, are integrable against the Gamma density, `E|X|^{J+1} ≤ C|τ|^{J+1}` and
`|B_N(τ) - Σ_{j≤J} binom(1/2,j) E[X^j]| ≤ C|τ|^{J+1}`. -/
def stmt_7 : Prop :=
  ∀ N : ℕ, 2 ≤ N → ∀ J : ℕ, ∃ C τJ : ℝ, 0 < τJ ∧ ∀ τ : ℝ, |τ| ≤ τJ →
    (∀ j ≤ J + 1, MeasureTheory.IntegrableOn
        (fun g => dens N τ g * (Ymean N τ g - 1) ^ j) (orth N)) ∧
    gammaE N τ (fun g => |Ymean N τ g - 1| ^ (J + 1)) ≤ C * |τ| ^ (J + 1) ∧
    |B N τ - ∑ j ∈ Finset.range (J + 1), binomHalf j * momX N j τ| ≤ C * |τ| ^ (J + 1)

/-- `lem:pw:moments` (sec_pureweight.tex 211–231), every `N ≥ 2`, every `j`: near `τ = 0`
the moment `X^j` is integrable and `E[X^j]` is real analytic at `0`; its Taylor coefficient of
`τ^k` is a rational zeta polynomial homogeneous of weight `k`, and it is `0` for `k < j`
(`𝐄_j ∈ τ^j 𝓡`, at the level of values). -/
def stmt_8 : Prop :=
  ∀ N : ℕ, 2 ≤ N → ∀ j : ℕ,
    (∃ r : ℝ, 0 < r ∧ ∀ τ : ℝ, |τ| < r →
      MeasureTheory.IntegrableOn (fun g => dens N τ g * (Ymean N τ g - 1) ^ j) (orth N)) ∧
    AnalyticAt ℝ (momX N j) 0 ∧
    (∀ k : ℕ, Q01.IsZetaPoly k (iteratedDeriv k (momX N j) 0 / (k.factorial : ℝ))) ∧
    (∀ k < j, iteratedDeriv k (momX N j) 0 = 0)

/-- `prop:pw:n` (sec_pureweight.tex 233–238) = `thm:pw`, last sentence (main.tex 1190): for
every dimension `N ≥ 2` there is one sequence `w_k(N)`, each a rational polynomial in
`ζ(2),…,ζ(k)` homogeneous of weight `k`, such that for every `K` there are `C, τ₀ > 0` with
`|log(S_N(p)/(S_N(1) Y_N(τ))) - Σ_{k=2}^K w_k(N) τ^k| ≤ C τ^{K+1}` for `p > 1`,
`τ = 1 - 1/p ≤ τ₀`. -/
def stmt_9 : Prop :=
  ∀ N : ℕ, 2 ≤ N → ∃ w : ℕ → ℝ, (∀ k, Q01.IsZetaPoly k (w k)) ∧
    ∀ K : ℕ, ∃ C τ₀ : ℝ, 0 < τ₀ ∧ ∀ p : ℝ, 1 < p → 1 - 1 / p ≤ τ₀ →
      |Real.log (S N p / (S N 1 * T04.Y N (1 - 1 / p))) -
          ∑ k ∈ Finset.Icc 2 K, w k * (1 - 1 / p) ^ k|
        ≤ C * (1 - 1 / p) ^ (K + 1)

/-- `g7:lem:both` (sec_g7_octa.tex 43–48) with `g7:lem:w23` (112–116): for every `N ≥ 2` the
expansion of stmt_9 holds for both signs of `τ` (`p > 0`, `|τ| ≤ τ₀`), with
`w_2(N) = -(N-1)(N²-2N-1)/(2N) ζ(2)` and `w_3(N) = -(N²-1)((N-1)³-2)/(3N²) ζ(3)`
(`= T04.coeffFormula N · ζ(3)`, the `ζ(3)` coefficient of main.tex 1992). -/
def stmt_10 : Prop :=
  ∀ N : ℕ, 2 ≤ N → ∃ w : ℕ → ℝ, (∀ k, Q01.IsZetaPoly k (w k)) ∧
    w 2 = -(((N : ℝ) - 1) * ((N : ℝ) ^ 2 - 2 * N - 1)) / (2 * N) * X05.z 2 ∧
    w 3 = (T04.coeffFormula N : ℝ) * X05.z 3 ∧
    ∀ K : ℕ, ∃ C τ₀ : ℝ, 0 < τ₀ ∧ ∀ p : ℝ, 0 < p → |1 - 1 / p| ≤ τ₀ →
      |Real.log (S N p / (S N 1 * T04.Y N (1 - 1 / p))) -
          ∑ k ∈ Finset.Icc 2 K, w k * (1 - 1 / p) ^ k|
        ≤ C * |1 - 1 / p| ^ (K + 1)

/-- `thm:oct`, the two weight claims (main.tex 1984, 1989; `g7:lem:struct`), for `N = 3` and
every `k`, both signs: (i) there is one sequence `w_k`, each a rational zeta polynomial
homogeneous of weight `k`, with `|log((2-p)S(p)/(4√3)) - Σ_{k=2}^K w_k τ^k| ≤ C|τ|^{K+1}` for
`p > 0`, `|τ| ≤ δ`, every `K`; (ii) there is one sequence `c_k` (`c_k = [d^k]S/(4√3)`), each a
rational polynomial in `ζ(2),…,ζ(k)` of weight at most `k`, with
`|S(1+d) - 4√3 Σ_{k≤K} c_k d^k| ≤ C|d|^{K+1}` for `|d| ≤ δ`, every `K`. -/
def stmt_11 : Prop :=
  (∃ w : ℕ → ℝ, (∀ k, Q01.IsZetaPoly k (w k)) ∧
    ∀ K : ℕ, ∃ C δ : ℝ, 0 < δ ∧ ∀ p : ℝ, 0 < p → |1 - 1 / p| ≤ δ →
      |Real.log ((2 - p) * Y05.S p / (4 * Real.sqrt 3)) -
          ∑ k ∈ Finset.Icc 2 K, w k * (1 - 1 / p) ^ k|
        ≤ C * |1 - 1 / p| ^ (K + 1)) ∧
  (∃ c : ℕ → ℝ, (∀ k, IsZetaPolyLE k (c k)) ∧
    ∀ K : ℕ, ∃ C δ : ℝ, 0 < δ ∧ ∀ d : ℝ, |d| ≤ δ →
      |Y05.S (1 + d) - 4 * Real.sqrt 3 * ∑ k ∈ Finset.range (K + 1), c k * d ^ k|
        ≤ C * |d| ^ (K + 1))

end R03
-- END STATEMENTS

open MeasureTheory Set

namespace R03

/-- the telescoped log-Gamma series: `log(Γ(1+u)^M/Γ(1+Mu))`. -/
lemma logGammaM {u r : ℝ} (hr : r < 1) (hu : |u| ≤ r) :
    ∀ M : ℕ, 1 ≤ M → |(M : ℝ) * u| ≤ r →
    HasSum (fun k : ℕ => (-1 : ℝ) ^ (k + 2) * X05.z (k + 2) / ((k + 2 : ℕ) : ℝ) *
        ((M : ℝ) * u ^ (k + 2) - ((M : ℝ) * u) ^ (k + 2)))
      (Real.log (Real.Gamma (1 + u) ^ M / Real.Gamma (1 + (M : ℝ) * u))) := by
  intro M hM
  induction M, hM using Nat.le_induction with
  | base => intro _; simp
  | succ M hM ih =>
    intro hM1
    have hMu : |(M : ℝ) * u| ≤ r := by
      have h1 : |(M : ℝ) * u| = M * |u| := by rw [abs_mul, abs_of_nonneg (by positivity)]
      have h2 : |((M + 1 : ℕ) : ℝ) * u| = (M + 1) * |u| := by
        rw [abs_mul, abs_of_nonneg (by positivity)]; push_cast; ring
      rw [h1]; rw [h2] at hM1; nlinarith [abs_nonneg u]
    have A := ih hMu
    have hsum : (M : ℝ) * u + u = ((M + 1 : ℕ) : ℝ) * u := by push_cast; ring
    have B := Y05.logGamma2 hr hMu hu (by rw [hsum]; exact hM1)
    rw [hsum] at B
    have hG : ∀ x : ℝ, |x| ≤ r → 0 < Real.Gamma (1 + x) := fun x hx =>
      Real.Gamma_pos_of_pos (by cases abs_lt.mp (lt_of_le_of_lt hx hr); linarith)
    have g0 := hG u hu; have g1 := hG _ hMu; have g2 := hG _ hM1
    have hlog : Real.log (Real.Gamma (1 + u) ^ (M + 1) / Real.Gamma (1 + ((M + 1 : ℕ) : ℝ) * u)) =
        Real.log (Real.Gamma (1 + u) ^ M / Real.Gamma (1 + (M : ℝ) * u)) +
        Real.log (Real.Gamma (1 + (M : ℝ) * u) * Real.Gamma (1 + u) /
          Real.Gamma (1 + ((M + 1 : ℕ) : ℝ) * u)) := by
      rw [← Real.log_mul (by positivity) (by positivity)]
      congr 1
      generalize Real.Gamma (1 + u) = a at *
      generalize Real.Gamma (1 + (M : ℝ) * u) = b at *
      generalize Real.Gamma (1 + ((M + 1 : ℕ) : ℝ) * u) = c at *
      field_simp
      ring
    rw [hlog]
    convert A.add B using 1
    ext k
    push_cast
    ring

end R03

theorem R03.proof_6 : R03.stmt_6 := by
  intro N hN τ hτ
  have hN1 : (1 : ℝ) ≤ (N : ℝ) - 1 := by
    have : (2 : ℝ) ≤ N := by exact_mod_cast hN
    linarith
  have hNpos : (0 : ℝ) < N := by linarith
  set u := R03.aPar N τ with hu
  set r := ((N : ℝ) - 1) * |τ| with hrdef
  have hr : r < 1 := by
    rw [hrdef]; rw [lt_div_iff₀ (by linarith)] at hτ; linarith
  have hNu : (N : ℝ) * u = -(((N : ℝ) - 1) * τ) := by
    rw [hu, R03.aPar]; field_simp
  have hNur : |(N : ℝ) * u| ≤ r := by
    rw [hNu, abs_neg, abs_mul, abs_of_nonneg (by linarith : (0 : ℝ) ≤ (N : ℝ) - 1)]
  have hur : |u| ≤ r := by
    have : |(N : ℝ) * u| = N * |u| := by rw [abs_mul, abs_of_pos hNpos]
    rw [this] at hNur
    have h2 : (1 : ℝ) ≤ N := by linarith
    nlinarith [abs_nonneg u]
  have H := R03.logGammaM hr hur N (by omega) hNur
  rw [hNu, ← sub_eq_add_neg] at H
  convert H using 1
  ext k
  set n := k + 2
  have hn : ((n : ℕ) : ℝ) ≠ 0 := by positivity
  have hN0 : (N : ℝ) ≠ 0 := hNpos.ne'
  have e0 : u = -((((N : ℝ) - 1) * τ) / N) := by rw [hu, R03.aPar]; ring
  have key : ∀ x : ℝ, (-1 : ℝ) ^ n * (-x) ^ n = x ^ n := fun x => by
    rw [← mul_pow, neg_one_mul, neg_neg]
  set w := ((N : ℝ) - 1) * τ with hw
  have hNn : (N : ℝ) ^ n = N * N ^ (k + 1) := by ring
  have hNk : (N : ℝ) ^ (k + 1) ≠ 0 := by positivity
  calc -(((N : ℝ) - 1) ^ n * (1 - 1 / (N : ℝ) ^ (k + 1))) / ((n : ℕ) : ℝ) * X05.z n * τ ^ n
      = X05.z n / ((n : ℕ) : ℝ) * (N * (w / N) ^ n - w ^ n) := by
        rw [hw, div_pow, mul_pow, hNn]
        field_simp
        ring
    _ = X05.z n / ((n : ℕ) : ℝ) * (N * ((-1 : ℝ) ^ n * u ^ n) - (-1 : ℝ) ^ n * (-w) ^ n) := by
        rw [e0, key, key]
    _ = _ := by ring

namespace R03

lemma face3 (τ : ℝ) (x : ℝ × ℝ) :
    faceIntegrand 2 τ ((MeasurableEquiv.finTwoArrow (α := ℝ)).symm x) = Y05.faceIntegrand τ x := by
  simp [faceIntegrand, Y05.faceIntegrand, ext, Fin.prod_univ_three, Fin.sum_univ_three,
    MeasurableEquiv.finTwoArrow, Fin.snoc, Fin.sum_univ_two, sub_add_eq_sub_sub]

lemma pre3 : (MeasurableEquiv.finTwoArrow (α := ℝ)).symm ⁻¹' simplex 2 = Y05.Δ := by
  ext x
  simp [simplex, Y05.Δ, MeasurableEquiv.finTwoArrow, Fin.forall_fin_two, Fin.sum_univ_two, and_assoc]

lemma face2 (τ : ℝ) (x : ℝ) :
    faceIntegrand 1 τ ((MeasurableEquiv.funUnique (Fin 1) ℝ).symm x) =
      (x * (1 - x)) ^ (-τ) * Real.sqrt (x ^ (2 * τ) + (1 - x) ^ (2 * τ)) := by
  simp [faceIntegrand, ext, Fin.prod_univ_two, Fin.sum_univ_two, Fin.snoc]

lemma pre2 : (MeasurableEquiv.funUnique (Fin 1) ℝ).symm ⁻¹' simplex 1 = Icc 0 1 := by
  ext x
  simp [simplex]

end R03

namespace R03

lemma int3 (τ : ℝ) : ∫ t in simplex 2, faceIntegrand 2 τ t = ∫ t in Y05.Δ, Y05.faceIntegrand τ t := by
  have h := ((volume_preserving_finTwoArrow ℝ).symm).setIntegral_preimage_emb
    (MeasurableEquiv.finTwoArrow (α := ℝ)).symm.measurableEmbedding (faceIntegrand 2 τ) (simplex 2)
  rw [← h, pre3]
  exact setIntegral_congr_fun Y05.Δ_meas (fun x _ => face3 τ x)

lemma int2 (τ : ℝ) : ∫ t in simplex 1, faceIntegrand 1 τ t =
    ∫ x in Icc (0 : ℝ) 1, (x * (1 - x)) ^ (-τ) * Real.sqrt (x ^ (2 * τ) + (1 - x) ^ (2 * τ)) := by
  have h := ((volume_preserving_funUnique (Fin 1) ℝ).symm).setIntegral_preimage_emb
    (MeasurableEquiv.funUnique (Fin 1) ℝ).symm.measurableEmbedding (faceIntegrand 1 τ) (simplex 1)
  refine Eq.trans h.symm ?_
  rw [pre2]
  exact setIntegral_congr_fun measurableSet_Icc (fun x _ => face2 τ x)

lemma pt2 {p x : ℝ} (hx0 : 0 < x) (hx1 : x < 1) :
    (x * (1 - x)) ^ (-(1 - 1 / p)) * Real.sqrt (x ^ (2 * (1 - 1 / p)) + (1 - x) ^ (2 * (1 - 1 / p)))
      = X01.integrand 1 1 p x := by
  set τ := 1 - 1 / p
  have hy : 0 < 1 - x := by linarith
  have e : 2 / p - 2 = -τ + -τ := by simp only [τ]; ring
  have ha : x ^ (-τ) * x ^ (-τ) * x ^ (2 * τ) = 1 := by
    rw [← Real.rpow_add hx0, ← Real.rpow_add hx0]; ring_nf; simp
  have hb : (1 - x) ^ (-τ) * (1 - x) ^ (-τ) * (1 - x) ^ (2 * τ) = 1 := by
    rw [← Real.rpow_add hy, ← Real.rpow_add hy]; ring_nf; simp
  have hA : 0 < (x * (1 - x)) ^ (-τ) := by positivity
  rw [X01.integrand, e, Real.rpow_add hx0, Real.rpow_add hy]
  rw [← Real.sqrt_sq hA.le, ← Real.sqrt_mul (sq_nonneg _)]
  congr 1
  rw [Real.mul_rpow hx0.le hy.le]
  linear_combination ((1 - x) ^ (-τ) * (1 - x) ^ (-τ)) * ha + (x ^ (-τ) * x ^ (-τ)) * hb

end R03

theorem R03.proof_3 : R03.stmt_3 := by
  refine ⟨fun p hp => ?_, fun p => ?_⟩
  · rw [R03.S, X01.P, intervalIntegral.integral_of_le zero_le_one, show 2 - 1 = 1 from rfl,
      R03.int2, integral_Icc_eq_integral_Ioo, integral_Ioc_eq_integral_Ioo]
    rw [setIntegral_congr_fun measurableSet_Ioo (g := X01.integrand 1 1 p)]
    · ring
    · exact fun x hx => R03.pt2 hx.1 hx.2
  · rw [R03.S, Y05.S, R03.int3]; norm_num

namespace R03

/-- the scaled simplex `{t ≥ 0, Σ t ≤ s}`. -/
def T (m : ℕ) (s : ℝ) : Set (Fin m → ℝ) := {t | (∀ i, 0 ≤ t i) ∧ ∑ i, t i ≤ s}

lemma T_closed (m : ℕ) (s : ℝ) : IsClosed (T m s) := by
  have : T m s = (⋂ i, {t : Fin m → ℝ | 0 ≤ t i}) ∩ {t | ∑ i, t i ≤ s} := by
    ext t; simp [T]
  rw [this]
  exact (isClosed_iInter fun i => isClosed_le continuous_const (continuous_apply i)).inter
    (isClosed_le (by fun_prop) continuous_const)

lemma T_empty {m : ℕ} {s : ℝ} (hs : s < 0) : T m s = ∅ := by
  ext t
  simp only [T, Set.mem_ofPred_eq, Set.mem_empty_iff_false, iff_false, not_and, not_le]
  intro h
  have := Finset.sum_nonneg (fun i (_ : i ∈ Finset.univ) => h i)
  linarith

lemma slice (m : ℕ) (s x : ℝ) :
    Prod.mk x ⁻¹' ((MeasurableEquiv.piFinSuccAbove (fun _ : Fin (m + 1) => ℝ) 0).symm ⁻¹'
      T (m + 1) s) = if 0 ≤ x then T m (s - x) else ∅ := by
  ext u
  split_ifs with hx
  · simp [T, MeasurableEquiv.piFinSuccAbove, Fin.forall_fin_succ, Fin.sum_univ_succ, hx,
      le_sub_iff_add_le']
  · simp [T, MeasurableEquiv.piFinSuccAbove, Fin.forall_fin_succ, hx]

lemma volT : ∀ m : ℕ, ∀ s : ℝ, 0 ≤ s →
    volume (T m s) = ENNReal.ofReal (s ^ m / (m.factorial : ℝ)) := by
  intro m
  induction m with
  | zero =>
    intro s hs
    have : T 0 s = Set.univ := by ext t; simp [T, hs]
    rw [this, volume_pi, Measure.pi_univ]; simp
  | succ m ih =>
    intro s hs
    have hP := (volume_preserving_piFinSuccAbove (fun _ : Fin (m + 1) => ℝ) 0).symm
    have hmeas : MeasurableSet ((MeasurableEquiv.piFinSuccAbove (fun _ : Fin (m + 1) => ℝ) 0).symm
        ⁻¹' T (m + 1) s) :=
      (MeasurableEquiv.measurable _) (T_closed _ _).measurableSet
    rw [← hP.measure_preimage (T_closed _ _).measurableSet.nullMeasurableSet, Measure.volume_eq_prod,
      Measure.prod_apply hmeas]
    simp_rw [slice]
    have hfun : (fun x : ℝ => volume (if 0 ≤ x then T m (s - x) else ∅)) =
        (Icc 0 s).indicator (fun x => ENNReal.ofReal ((s - x) ^ m / (m.factorial : ℝ))) := by
      ext x
      by_cases h0 : 0 ≤ x
      · by_cases h1 : x ≤ s
        · simp only [h0, ↓reduceIte]
          rw [Set.indicator_of_mem (show x ∈ Icc 0 s from ⟨h0, h1⟩), ih _ (by linarith)]
        · simp only [h0, ↓reduceIte]
          rw [Set.indicator_of_notMem (show x ∉ Icc 0 s from fun h => h1 h.2),
            T_empty (by linarith)]
          simp
      · simp only [h0, ↓reduceIte]
        rw [Set.indicator_of_notMem (show x ∉ Icc 0 s from fun h => h0 h.1)]; simp
    rw [hfun, lintegral_indicator measurableSet_Icc]
    rw [← ofReal_integral_eq_lintegral_ofReal]
    · congr 1
      rw [integral_Icc_eq_integral_Ioc, ← intervalIntegral.integral_of_le hs,
        intervalIntegral.integral_div, intervalIntegral.integral_comp_sub_left
          (fun x => x ^ m), integral_pow]
      push_cast [Nat.factorial_succ]
      field_simp
      ring
    · exact (Continuous.integrableOn_Icc (by fun_prop))
    · refine ae_restrict_of_forall_mem measurableSet_Icc (fun x hx => ?_)
      have : 0 ≤ s - x := by linarith [hx.2]
      simp only [Pi.zero_apply]
      positivity

end R03

theorem R03.proof_2 : R03.stmt_2 := by
  intro N hN
  have hS : R03.simplex (N - 1) = R03.T (N - 1) 1 := rfl
  have hc : ∀ t, R03.faceIntegrand (N - 1) (1 - 1 / 1) t = Real.sqrt N := by
    intro t
    simp only [R03.faceIntegrand, div_one, sub_self, neg_zero, Real.rpow_zero, mul_zero,
      Finset.sum_const, Finset.card_univ, Fintype.card_fin, one_mul, nsmul_eq_mul, mul_one]
    congr 1
    have : N - 1 + 1 = N := by omega
    rw [this]
  rw [R03.S, hS, setIntegral_congr_fun (R03.T_closed _ _).measurableSet (fun t _ => hc t),
    setIntegral_const, Measure.real, R03.volT _ _ zero_le_one,
    ENNReal.toReal_ofReal (by positivity)]
  simp
  field_simp


namespace R03

/-- `x^{-τ}` on `[0,1]`, `0` elsewhere. -/
noncomputable def f1 (τ : ℝ) (x : ℝ) : ℝ := (Icc (0 : ℝ) 1).indicator (fun x => x ^ (-τ)) x

lemma f1_meas (τ : ℝ) : Measurable (f1 τ) :=
  (measurable_id.pow_const _).indicator measurableSet_Icc

lemma f1_int {τ : ℝ} (hτ : τ < 1) : Integrable (f1 τ) :=
  (integrable_indicator_iff measurableSet_Icc).mpr (Y05.pow_int hτ)

lemma f1_nn (τ x : ℝ) : 0 ≤ f1 τ x := by
  unfold f1
  by_cases h : x ∈ Icc (0 : ℝ) 1
  · rw [indicator_of_mem h]; exact Real.rpow_nonneg h.1 _
  · rw [indicator_of_notMem h]

lemma f1_eq {τ x : ℝ} (h0 : 0 ≤ x) (h1 : x ≤ 1) : f1 τ x = x ^ (-τ) := by
  unfold f1; rw [indicator_of_mem (show x ∈ Icc (0 : ℝ) 1 from ⟨h0, h1⟩)]

/-- the dominating function in which coordinate `j` is replaced by `1 - Σ t`. -/
noncomputable def phi (τ : ℝ) (n : ℕ) (j : Fin (n + 1)) (t : Fin (n + 1) → ℝ) : ℝ :=
  (∏ i : Fin n, f1 τ (t (j.succAbove i))) * f1 τ (1 - ∑ i, t i)

lemma phi_int {τ : ℝ} (hτ : τ < 1) (n : ℕ) (j : Fin (n + 1)) : Integrable (phi τ n j) := by
  have hP := (volume_preserving_piFinSuccAbove (fun _ : Fin (n + 1) => ℝ) j).symm
  rw [← hP.integrable_comp_emb (MeasurableEquiv.measurableEmbedding _)]
  have hc : (phi τ n j ∘ (MeasurableEquiv.piFinSuccAbove (fun _ : Fin (n + 1) => ℝ) j).symm) =
      fun q : ℝ × (Fin n → ℝ) => (∏ i, f1 τ (q.2 i)) * f1 τ (1 - (q.1 + ∑ i, q.2 i)) := by
    ext q
    simp [phi, MeasurableEquiv.piFinSuccAbove_symm_apply, Fin.insertNthEquiv,
      Fin.sum_univ_succAbove _ j]
  rw [hc]
  have hm : Measurable fun q : ℝ × (Fin n → ℝ) => (∏ i, f1 τ (q.2 i)) * f1 τ (1 - (q.1 + ∑ i, q.2 i)) := by
    have := f1_meas τ
    fun_prop
  rw [Measure.volume_eq_prod, integrable_prod_iff' hm.aestronglyMeasurable]
  constructor
  · refine Filter.Eventually.of_forall (fun u => ?_)
    exact ((f1_int hτ).comp_sub_left (1 - ∑ i, u i)).const_mul (∏ i, f1 τ (u i)) |>.congr
      (Filter.Eventually.of_forall fun x => by
        simp only [show 1 - (x + ∑ i, u i) = 1 - ∑ i, u i - x by ring])
  · have hF : Integrable (fun u : Fin n → ℝ => ∏ i, f1 τ (u i)) :=
      Integrable.fintype_prod (fun _ => f1_int hτ)
    refine (hF.norm.mul_const (∫ x, ‖f1 τ x‖)).congr (Filter.Eventually.of_forall fun u => ?_)
    simp only [norm_mul]
    rw [integral_const_mul]
    congr 1
    have := integral_sub_left_eq_self (fun x => ‖f1 τ x‖) (volume : Measure ℝ) (1 - ∑ i, u i)
    rw [← this]
    congr 1; ext x; ring_nf

end R03

namespace R03

lemma ext_cs {m : ℕ} (t : Fin m → ℝ) (i : Fin m) : ext t (Fin.castSucc i) = t i := by
  simp [ext]

lemma ext_last {m : ℕ} (t : Fin m → ℝ) : ext t (Fin.last m) = 1 - ∑ i, t i := by
  simp [ext]

lemma ext_bd {m : ℕ} {t : Fin m → ℝ} (ht : t ∈ simplex m) (i : Fin (m + 1)) :
    0 ≤ ext t i ∧ ext t i ≤ 1 := by
  have hs : 0 ≤ ∑ i, t i := Finset.sum_nonneg (fun i _ => ht.1 i)
  induction i using Fin.lastCases with
  | last => rw [ext_last]; constructor <;> linarith [ht.2]
  | cast j =>
    rw [ext_cs]
    refine ⟨ht.1 j, ?_⟩
    have := Finset.single_le_sum (fun i (_ : i ∈ Finset.univ) => ht.1 i) (Finset.mem_univ j)
    linarith [ht.2]

lemma ext_sum {m : ℕ} (t : Fin m → ℝ) : ∑ i, ext t i = 1 := by
  rw [Fin.sum_univ_castSucc]; simp only [ext_cs, ext_last]; ring

lemma ext_cont {m : ℕ} (i : Fin (m + 1)) : Continuous fun t : Fin m → ℝ => ext t i := by
  induction i using Fin.lastCases with
  | last => simp only [ext_last]; fun_prop
  | cast j => simp only [ext_cs]; fun_prop

lemma face_meas (m : ℕ) (τ : ℝ) : Measurable (faceIntegrand m τ) := by
  change Measurable fun t => (∏ i, ext t i) ^ (-τ) * Real.sqrt (∑ i, ext t i ^ (2 * τ))
  have h : ∀ i, Measurable fun t : Fin m → ℝ => ext t i := fun i => (ext_cont i).measurable
  have h1 : Measurable fun t : Fin m → ℝ => (∏ i, ext t i) ^ (-τ) :=
    (Finset.measurable_prod (f := fun i t => ext t i) Finset.univ (fun i _ => h i)).pow_const (-τ)
  have h2 : Measurable fun t : Fin m → ℝ => Real.sqrt (∑ i, ext t i ^ (2 * τ)) :=
    (Finset.measurable_sum (f := fun i t => ext t i ^ (2 * τ)) Finset.univ
      (fun i _ => (h i).pow_const _)).sqrt
  exact h1.mul h2

lemma simplex_meas (m : ℕ) : MeasurableSet (simplex m) := (T_closed m 1).measurableSet

lemma sqrt_bd {m : ℕ} {τ : ℝ} (hτ : 0 ≤ τ) {t : Fin m → ℝ} (ht : t ∈ simplex m) :
    Real.sqrt (∑ i, ext t i ^ (2 * τ)) ≤ Real.sqrt ((m : ℝ) + 1) := by
  apply Real.sqrt_le_sqrt
  calc ∑ i, ext t i ^ (2 * τ) ≤ ∑ _i : Fin (m + 1), (1 : ℝ) :=
        Finset.sum_le_sum (fun i _ => Real.rpow_le_one (ext_bd ht i).1 (ext_bd ht i).2 (by linarith))
    _ = (m : ℝ) + 1 := by simp

/-- the bound for `τ ≤ 0`. -/
lemma face_bd_neg {m : ℕ} {τ : ℝ} (hτ : τ ≤ 0) {t : Fin m → ℝ} (ht : t ∈ simplex m) :
    ‖faceIntegrand m τ t‖ ≤ Real.sqrt ((m : ℝ) + 1) := by
  have hb := ext_bd ht
  have hP0 : 0 ≤ ∏ i, ext t i := Finset.prod_nonneg (fun i _ => (hb i).1)
  have hA : 0 ≤ (∏ i, ext t i) ^ (-τ) := Real.rpow_nonneg hP0 _
  have hS : 0 ≤ ∑ i, ext t i ^ (2 * τ) := Finset.sum_nonneg (fun i _ => Real.rpow_nonneg (hb i).1 _)
  unfold faceIntegrand
  rw [Real.norm_eq_abs, abs_of_nonneg (mul_nonneg hA (Real.sqrt_nonneg _))]
  rw [← Real.sqrt_sq hA, ← Real.sqrt_mul (sq_nonneg _)]
  apply Real.sqrt_le_sqrt
  rw [Finset.mul_sum]
  calc ∑ i, ((∏ i, ext t i) ^ (-τ)) ^ 2 * ext t i ^ (2 * τ) ≤ ∑ _i : Fin (m + 1), (1 : ℝ) := by
        apply Finset.sum_le_sum
        intro i _
        rcases (hb i).1.lt_or_eq with hpos | hzero
        · have hle : ∏ j, ext t j ≤ ext t i := by
            rw [← Finset.mul_prod_erase _ _ (Finset.mem_univ i)]
            exact mul_le_of_le_one_right (hb i).1
              (Finset.prod_le_one₀ (fun j _ => (hb j).1) (fun j _ => (hb j).2))
          have h1 : (∏ j, ext t j) ^ (-τ) ≤ ext t i ^ (-τ) :=
            Real.rpow_le_rpow hP0 hle (by linarith)
          have h2 : (ext t i ^ (-τ)) ^ 2 * ext t i ^ (2 * τ) = 1 := by
            rw [← Real.rpow_natCast, ← Real.rpow_mul hpos.le, ← Real.rpow_add hpos]
            ring_nf; simp
          calc ((∏ i, ext t i) ^ (-τ)) ^ 2 * ext t i ^ (2 * τ)
              ≤ (ext t i ^ (-τ)) ^ 2 * ext t i ^ (2 * τ) := by
                gcongr
            _ = 1 := h2
        · rcases hτ.lt_or_eq with hneg | h0
          · have hP : ∏ j, ext t j = 0 := Finset.prod_eq_zero (Finset.mem_univ i) hzero.symm
            rw [hP, Real.zero_rpow (by linarith)]; norm_num
          · rw [h0]; simp
    _ = (m : ℝ) + 1 := by simp

end R03

namespace R03

lemma f1_le {τ x M : ℝ} (hτ : 0 ≤ τ) (hM : 0 < M) (hx : 1 / M ≤ x) (hx1 : x ≤ 1) :
    f1 τ x ≤ M ^ τ := by
  have hx0 : 0 < x := lt_of_lt_of_le (by positivity) hx
  rw [f1_eq hx0.le hx1]
  calc x ^ (-τ) ≤ (1 / M) ^ (-τ) := Real.rpow_le_rpow_of_nonpos (by positivity) hx (by linarith)
    _ = M ^ τ := by
      rw [Real.rpow_neg (by positivity), Real.div_rpow zero_le_one hM.le, Real.one_rpow,
        inv_div, div_one]

lemma face_bd_pos {n : ℕ} {τ : ℝ} (hτ0 : 0 < τ) {t : Fin (n + 1) → ℝ}
    (ht : t ∈ simplex (n + 1)) :
    ‖faceIntegrand (n + 1) τ t‖ ≤ Real.sqrt (((n + 1 : ℕ) : ℝ) + 1) *
      ((((n + 1 : ℕ) : ℝ) + 1) ^ τ * ((∏ i, f1 τ (t i)) + ∑ j, phi τ n j t)) := by
  set M : ℝ := ((n + 1 : ℕ) : ℝ) + 1 with hMdef
  have hM : 0 < M := by positivity
  have hb := ext_bd ht
  have hP0 : 0 ≤ ∏ i, ext t i := Finset.prod_nonneg (fun i _ => (hb i).1)
  have hA : 0 ≤ (∏ i, ext t i) ^ (-τ) := Real.rpow_nonneg hP0 _
  have hs : 0 ≤ ∑ i, t i := Finset.sum_nonneg (fun i _ => ht.1 i)
  have hQ : (∏ i, ext t i) ^ (-τ) = (∏ i, f1 τ (t i)) * f1 τ (1 - ∑ i, t i) := by
    rw [← Real.finsetProd_rpow _ _ (fun i _ => (hb i).1), Fin.prod_univ_castSucc]
    simp only [ext_cs, ext_last]
    congr 1
    · exact Finset.prod_congr rfl (fun i _ => (f1_eq (ht.1 i) (by
        have := Finset.single_le_sum (fun i (_ : i ∈ Finset.univ) => ht.1 i) (Finset.mem_univ i)
        linarith [ht.2])).symm)
    · exact (f1_eq (by linarith [ht.2]) (by linarith)).symm
  have hF0 : 0 ≤ ∏ i, f1 τ (t i) := Finset.prod_nonneg (fun i _ => f1_nn _ _)
  have hphi0 : ∀ j, 0 ≤ phi τ n j t := fun j =>
    mul_nonneg (Finset.prod_nonneg (fun i _ => f1_nn _ _)) (f1_nn _ _)
  have hsum0 : 0 ≤ ∑ j, phi τ n j t := Finset.sum_nonneg (fun j _ => hphi0 j)
  have hMt : 0 ≤ M ^ τ := by positivity
  -- some coordinate is at least `1/M`
  obtain ⟨k, -, hk⟩ := Finset.exists_le_of_sum_le (s := Finset.univ) Finset.univ_nonempty
    (f := fun _ => 1 / M) (g := fun i => ext t i) (by
      rw [ext_sum]; simp [hMdef]; field_simp; exact le_rfl)
  have hmain : (∏ i, ext t i) ^ (-τ) ≤ M ^ τ * ((∏ i, f1 τ (t i)) + ∑ j, phi τ n j t) := by
    rw [hQ]
    induction k using Fin.lastCases with
    | last =>
      rw [ext_last] at hk
      have := f1_le hτ0.le hM hk (by linarith)
      calc (∏ i, f1 τ (t i)) * f1 τ (1 - ∑ i, t i) ≤ (∏ i, f1 τ (t i)) * M ^ τ := by gcongr
        _ ≤ M ^ τ * ((∏ i, f1 τ (t i)) + ∑ j, phi τ n j t) := by nlinarith
    | cast j =>
      rw [ext_cs] at hk
      have := f1_le hτ0.le hM hk (by have := (hb (Fin.castSucc j)).2; rwa [ext_cs] at this)
      have hsplit : (∏ i, f1 τ (t i)) * f1 τ (1 - ∑ i, t i) = f1 τ (t j) * phi τ n j t := by
        rw [Fin.prod_univ_succAbove _ j, phi]; ring
      rw [hsplit]
      have h1 : phi τ n j t ≤ ∑ j, phi τ n j t :=
        Finset.single_le_sum (fun j _ => hphi0 j) (Finset.mem_univ j)
      calc f1 τ (t j) * phi τ n j t ≤ M ^ τ * phi τ n j t := by gcongr; exact hphi0 j
        _ ≤ M ^ τ * ((∏ i, f1 τ (t i)) + ∑ j, phi τ n j t) := by gcongr; linarith
  have hsq := sqrt_bd hτ0.le ht
  unfold faceIntegrand
  rw [Real.norm_eq_abs, abs_of_nonneg (mul_nonneg hA (Real.sqrt_nonneg _))]
  calc (∏ i, ext t i) ^ (-τ) * Real.sqrt (∑ i, ext t i ^ (2 * τ))
      ≤ (M ^ τ * ((∏ i, f1 τ (t i)) + ∑ j, phi τ n j t)) * Real.sqrt M := by
        rw [hQ] at hA
        exact mul_le_mul hmain (by simpa [hMdef] using hsq) (Real.sqrt_nonneg _) (by positivity)
    _ = _ := by ring

end R03

theorem R03.proof_1 : R03.stmt_1 := by
  intro N hN p hp
  obtain ⟨n, hn⟩ : ∃ n, N - 1 = n + 1 := ⟨N - 2, by omega⟩
  rw [hn]
  set τ := 1 - 1 / p with hτdef
  have hτ1 : τ < 1 := by have : 0 < 1 / p := by positivity
                         linarith
  have hfin : volume (R03.simplex (n + 1)) ≠ ⊤ := by
    show volume (R03.T (n + 1) 1) ≠ ⊤
    rw [R03.volT _ _ zero_le_one]; exact ENNReal.ofReal_ne_top
  have hmeas := (R03.face_meas (n + 1) τ).aestronglyMeasurable
    (μ := volume.restrict (R03.simplex (n + 1)))
  by_cases hτ : τ ≤ 0
  · exact Integrable.mono' (integrableOn_const hfin) hmeas
      (ae_restrict_of_forall_mem (R03.simplex_meas _) (fun t ht => R03.face_bd_neg hτ ht))
  · replace hτ : 0 < τ := lt_of_not_ge hτ
    have hg : Integrable (fun t : Fin (n + 1) → ℝ => Real.sqrt (((n + 1 : ℕ) : ℝ) + 1) *
        ((((n + 1 : ℕ) : ℝ) + 1) ^ τ * ((∏ i, R03.f1 τ (t i)) + ∑ j, R03.phi τ n j t))) :=
      (((Integrable.fintype_prod (fun _ => R03.f1_int hτ1)).add
        (integrable_finsetSum _ (fun j _ => R03.phi_int hτ1 n j))).const_mul _).const_mul _
    exact Integrable.mono' hg.integrableOn hmeas
      (ae_restrict_of_forall_mem (R03.simplex_meas _) (fun t ht => R03.face_bd_pos hτ ht))


namespace R03

open scoped ENNReal

/-- graphs over the other coordinates are Lebesgue-null. -/
lemma graph_null {k : ℕ} (j : Fin (k + 1)) (h : (Fin k → ℝ) → ℝ) (hh : Measurable h) :
    volume {t : Fin (k + 1) → ℝ | t j = h (fun i => t (j.succAbove i))} = 0 := by
  have hP := volume_preserving_piFinSuccAbove (fun _ : Fin (k + 1) => ℝ) j
  have hA : MeasurableSet {q : ℝ × (Fin k → ℝ) | q.1 = h q.2} :=
    measurableSet_eq_fun measurable_fst (hh.comp measurable_snd)
  have : {t : Fin (k + 1) → ℝ | t j = h (fun i => t (j.succAbove i))} =
      MeasurableEquiv.piFinSuccAbove (fun _ => ℝ) j ⁻¹' {q | q.1 = h q.2} := by
    ext t; simp [MeasurableEquiv.piFinSuccAbove_apply, Fin.insertNthEquiv]; rfl
  rw [this, hP.measure_preimage hA.nullMeasurableSet, Measure.volume_eq_prod,
    Measure.prod_apply_symm hA]
  simp

/-- the open simplex. -/
def opS (m : ℕ) : Set (Fin m → ℝ) := {t | ∀ i, 0 < ext t i}

lemma opS_meas (m : ℕ) : MeasurableSet (opS m) := by
  have : opS m = ⋂ i, {t | 0 < ext t i} := by ext; simp [opS]
  rw [this]
  exact MeasurableSet.iInter fun i => measurableSet_lt measurable_const (ext_cont i).measurable

lemma opS_sub (m : ℕ) : opS m ⊆ simplex m := by
  intro t ht
  refine ⟨fun i => ?_, ?_⟩
  · have := ht (Fin.castSucc i); rw [ext_cs] at this; exact this.le
  · have := ht (Fin.last m); rw [ext_last] at this; linarith

lemma simplex_ae (k : ℕ) : opS (k + 1) =ᵐ[volume] simplex (k + 1) := by
  refine (ae_eq_set).mpr ⟨by rw [sdiff_eq_empty.mpr (opS_sub _)]; simp, ?_⟩
  have h1 : ∀ i : Fin (k + 1), volume {t : Fin (k + 1) → ℝ | t i = 0} = 0 := fun i => by
    simpa using graph_null i (fun _ => (0 : ℝ)) measurable_const
  have h2 : volume {t : Fin (k + 1) → ℝ | t 0 = 1 - ∑ i : Fin k, t i.succ} = 0 := by
    simpa using graph_null 0 (fun u : Fin k → ℝ => 1 - ∑ i, u i) (by fun_prop)
  apply measure_mono_null (t := (⋃ i, {t : Fin (k + 1) → ℝ | t i = 0}) ∪
      {t | t 0 = 1 - ∑ i : Fin k, t i.succ})
  · rintro t ⟨ht, hno⟩
    simp only [opS, mem_ofPred_eq, not_forall, not_lt] at hno
    obtain ⟨i, hi⟩ := hno
    induction i using Fin.lastCases with
    | last =>
      right
      rw [ext_last] at hi
      simp only [mem_ofPred_eq]
      have := ht.2
      rw [Fin.sum_univ_succ] at this hi
      linarith
    | cast j =>
      left
      rw [ext_cs] at hi
      simp only [mem_iUnion, mem_ofPred_eq]
      exact ⟨j, le_antisymm hi (ht.1 j)⟩
  · exact measure_union_null (measure_iUnion_null h1) h2

/-- the integrand of `lem:pw:dirichlet`, as an indicator on the whole space. -/
noncomputable def dirPhi (m : ℕ) (a : ℝ) (F : (Fin (m + 1) → ℝ) → ℝ) (g : Fin (m + 1) → ℝ) :
    ℝ≥0∞ :=
  (orth (m + 1)).indicator
    (fun g => ENNReal.ofReal ((∏ i, g i ^ a) * Real.exp (-∑ i, g i) * F g)) g

lemma orth_meas (n : ℕ) : MeasurableSet (orth n) :=
  MeasurableSet.univ_pi (fun _ => measurableSet_Ioi)

lemma dirPhi_meas (m : ℕ) (a : ℝ) {F : (Fin (m + 1) → ℝ) → ℝ} (hF : Measurable F) :
    Measurable (dirPhi m a F) := by
  unfold dirPhi
  refine Measurable.indicator ?_ (orth_meas _)
  have h1 : Measurable fun g : Fin (m + 1) → ℝ => ∏ i, g i ^ a :=
    Finset.measurable_prod (f := fun i (g : Fin (m + 1) → ℝ) => g i ^ a) _
      (fun i _ => (measurable_pi_apply i).pow_const a)
  have h2 : Measurable fun g : Fin (m + 1) → ℝ => Real.exp (-∑ i, g i) :=
    Real.measurable_exp.comp (Finset.measurable_sum (f := fun i (g : Fin (m + 1) → ℝ) => g i) _
      (fun i _ => measurable_pi_apply i)).neg
  exact ENNReal.measurable_ofReal.comp ((h1.mul h2).mul hF)

lemma snoc_sum {m : ℕ} (u : Fin m → ℝ) (y : ℝ) :
    ∑ i, (Fin.snoc (α := fun _ => ℝ) u y : Fin (m + 1) → ℝ) i = ∑ i, u i + y := by
  rw [Fin.sum_univ_castSucc]; simp

lemma snoc_smul {m : ℕ} {s : ℝ} (t : Fin m → ℝ) :
    (Fin.snoc (α := fun _ => ℝ) (s • t) (s - ∑ i, (s • t) i) : Fin (m + 1) → ℝ) = s • ext t := by
  funext i
  induction i using Fin.lastCases with
  | last => simp [ext]; rw [← Finset.mul_sum]; ring
  | cast j => simp [ext]

lemma dirPhi_neg {m : ℕ} {a : ℝ} {F : (Fin (m + 1) → ℝ) → ℝ} {x : ℝ} (hx : x ≤ 0)
    (u : Fin m → ℝ) :
    dirPhi m a F (Fin.snoc (α := fun _ => ℝ) u (x - ∑ i, u i)) = 0 := by
  unfold dirPhi
  rw [indicator_of_notMem]
  intro h
  have hpos : 0 < ∑ i, (Fin.snoc (α := fun _ => ℝ) u (x - ∑ i, u i) : Fin (m + 1) → ℝ) i :=
    Finset.sum_pos (fun i _ => h i (mem_univ i)) Finset.univ_nonempty
  rw [snoc_sum] at hpos
  linarith

lemma dirPhi_pos {m : ℕ} {a : ℝ} {F : (Fin (m + 1) → ℝ) → ℝ}
    (hF : ∀ c : ℝ, 0 < c → ∀ g ∈ orth (m + 1), F (c • g) = F g) {s : ℝ} (hs : 0 < s)
    (t : Fin m → ℝ) :
    dirPhi m a F (s • ext t) = (opS m).indicator (fun t =>
      ENNReal.ofReal ((s ^ a) ^ (m + 1) * Real.exp (-s)) *
        ENNReal.ofReal ((∏ i, ext t i ^ a) * F (ext t))) t := by
  unfold dirPhi
  by_cases ht : t ∈ opS m
  · have hmem : ext t ∈ orth (m + 1) := fun i _ => ht i
    have hmem' : s • ext t ∈ orth (m + 1) := fun i _ => by
      simp only [Pi.smul_apply, smul_eq_mul, mem_Ioi]; exact mul_pos hs (ht i)
    rw [indicator_of_mem hmem', indicator_of_mem ht, ← ENNReal.ofReal_mul (by positivity),
      hF s hs _ hmem]
    congr 1
    simp only [Pi.smul_apply, smul_eq_mul]
    rw [← Finset.mul_sum, ext_sum, mul_one,
      Finset.prod_congr rfl (fun i _ => Real.mul_rpow hs.le (ht i).le), Finset.prod_mul_distrib,
      Finset.prod_const, Finset.card_univ, Fintype.card_fin]
    ring
  · rw [indicator_of_notMem ht, indicator_of_notMem]
    intro h
    apply ht
    intro i
    have := h i (mem_univ i)
    simp only [Pi.smul_apply, smul_eq_mul, mem_Ioi] at this
    exact pos_of_mul_pos_right this hs.le

lemma lint_smul {m : ℕ} {s : ℝ} (hs : 0 < s) {f : (Fin m → ℝ) → ℝ≥0∞} (hf : Measurable f) :
    ∫⁻ u, f u = ENNReal.ofReal (s ^ m) * ∫⁻ t, f (s • t) := by
  have h1 : ∫⁻ t, f (s • t) = ∫⁻ y, f y ∂(Measure.map (s • ·) (volume : Measure (Fin m → ℝ))) :=
    (lintegral_map hf (measurable_const_smul s)).symm
  rw [h1, Measure.map_addHaar_smul _ hs.ne', lintegral_smul_measure, Module.finrank_fin_fun,
    smul_eq_mul, ← mul_assoc, abs_of_pos (by positivity), ← ENNReal.ofReal_mul (by positivity),
    mul_inv_cancel₀ (by positivity), ENNReal.ofReal_one, one_mul]

lemma dirichlet (n : ℕ) (hae : opS n =ᵐ[volume] simplex n) (a : ℝ) (ha : -1 < a)
    (F : (Fin (n + 1) → ℝ) → ℝ) (hFm : Measurable F)
    (hF : ∀ c : ℝ, 0 < c → ∀ g ∈ orth (n + 1), F (c • g) = F g) :
    ∫⁻ g in orth (n + 1), ENNReal.ofReal ((∏ i, g i ^ a) * Real.exp (-∑ i, g i) * F g) =
      ENNReal.ofReal (Real.Gamma (((n : ℝ) + 1) * (1 + a))) *
        ∫⁻ t in simplex n, ENNReal.ofReal ((∏ i, ext t i ^ a) * F (ext t)) := by
  set I := ∫⁻ t in simplex n, ENNReal.ofReal ((∏ i, ext t i ^ a) * F (ext t)) with hI
  set e := MeasurableEquiv.piFinSuccAbove (fun _ : Fin (n + 1) => ℝ) (Fin.last n) with he
  have hP : MeasurePreserving e.symm (volume.prod volume) volume :=
    (volume_preserving_piFinSuccAbove (fun _ : Fin (n + 1) => ℝ) (Fin.last n)).symm
  have hes : ∀ x u, e.symm (x, u) = Fin.snoc (α := fun _ => ℝ) u x := by
    intro x u
    simp [he, MeasurableEquiv.piFinSuccAbove_symm_apply, Fin.insertNthEquiv]
  have hΦ := dirPhi_meas n a hFm
  have step1 : ∫⁻ g in orth (n + 1), ENNReal.ofReal ((∏ i, g i ^ a) * Real.exp (-∑ i, g i) * F g)
      = ∫⁻ u, ∫⁻ x, dirPhi n a F (Fin.snoc (α := fun _ => ℝ) u (x - ∑ i, u i)) := by
    rw [← lintegral_indicator (orth_meas _)]
    change ∫⁻ g, dirPhi n a F g = _
    rw [← hP.lintegral_comp_emb e.symm.measurableEmbedding,
      lintegral_prod_symm (fun z => dirPhi n a F (e.symm z)) (hΦ.comp e.symm.measurable).aemeasurable]
    refine lintegral_congr (fun u => ?_)
    simp only [hes]
    exact (lintegral_sub_right_eq_self
      (fun x => dirPhi n a F (Fin.snoc (α := fun _ => ℝ) u x)) (∑ i, u i)).symm
  have hmeas2 : Measurable (Function.uncurry fun (u : Fin n → ℝ) (x : ℝ) =>
      dirPhi n a F (Fin.snoc (α := fun _ => ℝ) u (x - ∑ i, u i))) := by
    have : (Function.uncurry fun (u : Fin n → ℝ) (x : ℝ) =>
        dirPhi n a F (Fin.snoc (α := fun _ => ℝ) u (x - ∑ i, u i))) =
        fun q : (Fin n → ℝ) × ℝ => dirPhi n a F (e.symm (q.2 - ∑ i, q.1 i, q.1)) := by
      funext ⟨u, x⟩; simp [hes]
    rw [this]
    exact hΦ.comp (e.symm.measurable.comp (by fun_prop))
  rw [step1, lintegral_lintegral_swap hmeas2.aemeasurable]
  have step2 : ∀ x : ℝ, ∫⁻ u, dirPhi n a F (Fin.snoc (α := fun _ => ℝ) u (x - ∑ i, u i)) =
      (Ioi (0 : ℝ)).indicator (fun s => ENNReal.ofReal (s ^ n * ((s ^ a) ^ (n + 1) *
        Real.exp (-s)))) x * I := by
    intro x
    by_cases hx : 0 < x
    · rw [indicator_of_mem (show x ∈ Ioi (0 : ℝ) from hx)]
      have hm : Measurable fun u : Fin n → ℝ =>
          dirPhi n a F (Fin.snoc (α := fun _ => ℝ) u (x - ∑ i, u i)) := by
        have : (fun u : Fin n → ℝ => dirPhi n a F (Fin.snoc (α := fun _ => ℝ) u (x - ∑ i, u i)))
            = fun u => dirPhi n a F (e.symm (x - ∑ i, u i, u)) := by funext u; simp [hes]
        rw [this]; exact hΦ.comp (e.symm.measurable.comp (by fun_prop))
      rw [lint_smul hx hm]
      simp_rw [snoc_smul, dirPhi_pos hF hx]
      rw [lintegral_indicator (opS_meas _), lintegral_const_mul' _ _ ENNReal.ofReal_ne_top,
        setLIntegral_congr hae, ← mul_assoc, ← ENNReal.ofReal_mul (by positivity)]
    · rw [indicator_of_notMem (show x ∉ Ioi (0 : ℝ) from hx), zero_mul]
      simp only [dirPhi_neg (not_lt.mp hx), lintegral_zero]
  simp_rw [step2]
  rw [lintegral_mul_const _ ((Measurable.ennreal_ofReal (by fun_prop)).indicator
      measurableSet_Ioi), lintegral_indicator measurableSet_Ioi]
  congr 1
  have hs : 0 < (((n : ℕ) : ℝ) + 1) * (1 + a) := by
    have : (0 : ℝ) < (n : ℝ) + 1 := by positivity
    have : 0 < 1 + a := by linarith
    positivity
  rw [Real.Gamma_eq_integral hs, ofReal_integral_eq_lintegral_ofReal
    (Real.GammaIntegral_convergent hs)
    (ae_restrict_of_forall_mem measurableSet_Ioi (fun x hx => by
      have : (0 : ℝ) < x := hx
      simp only [Pi.zero_apply]; positivity))]
  refine setLIntegral_congr_fun measurableSet_Ioi (fun x hx => ?_)
  have hx0 : (0 : ℝ) < x := hx
  congr 1
  rw [← Real.rpow_mul_natCast hx0.le, ← Real.rpow_natCast, ← mul_assoc, ← Real.rpow_add hx0,
    mul_comm]
  congr 2
  push_cast
  ring

end R03

theorem R03.proof_4 : R03.stmt_4 := by
  intro m hm a ha F hFm hF0 hF
  obtain ⟨k, rfl⟩ : ∃ k, m = k + 1 := ⟨m - 1, by omega⟩
  exact R03.dirichlet (k + 1) (R03.simplex_ae k) a ha F hFm hF

namespace R03

lemma Ymean_hom {N : ℕ} (hN : 0 < N) (τ : ℝ) {c : ℝ} (hc : 0 < c) {g : Fin N → ℝ}
    (hg : g ∈ orth N) : Ymean N τ (c • g) = Ymean N τ g := by
  have hN' : (N : ℝ) ≠ 0 := by positivity
  have hgi : ∀ i, 0 < g i := fun i => hg i (mem_univ i)
  unfold Ymean
  congr 1
  refine Finset.sum_congr rfl (fun i _ => ?_)
  congr 2
  simp only [Pi.smul_apply, smul_eq_mul]
  rw [Finset.sum_congr rfl (fun j _ => Real.log_mul hc.ne' (hgi j).ne'), Finset.sum_add_distrib,
    Finset.sum_const, Finset.card_univ, Fintype.card_fin, nsmul_eq_mul, Real.log_mul hc.ne'
    (hgi i).ne']
  field_simp
  ring

lemma Ymean_meas (N : ℕ) (τ : ℝ) : Measurable (Ymean N τ) := by
  unfold Ymean
  refine Measurable.div_const ?_ _
  refine Finset.measurable_sum _ (fun i _ => ?_)
  refine Real.measurable_exp.comp (Measurable.const_mul ?_ _)
  refine (Real.measurable_log.comp (measurable_pi_apply i)).sub (Measurable.div_const ?_ _)
  exact Finset.measurable_sum _ (fun j _ => Real.measurable_log.comp (measurable_pi_apply j))

/-- the centring identity, pointwise. -/
lemma face_key {N : ℕ} (hN : 0 < N) (τ : ℝ) {x : Fin N → ℝ} (hx : ∀ i, 0 < x i) :
    (∏ i, x i) ^ (-τ) * Real.sqrt (∑ i, x i ^ (2 * τ)) =
      Real.sqrt N * ((∏ i, x i ^ aPar N τ) * Real.sqrt (Ymean N τ x)) := by
  have hN' : (N : ℝ) ≠ 0 := by positivity
  set L := ∑ i, Real.log (x i) with hL
  have hP : ∏ i, x i = Real.exp L := by
    rw [hL, Real.exp_sum]; exact Finset.prod_congr rfl (fun i _ => (Real.exp_log (hx i)).symm)
  have hPa : ∏ i, x i ^ aPar N τ = Real.exp (L * aPar N τ) := by
    rw [Real.finsetProd_rpow _ _ (fun i _ => (hx i).le), hP, ← Real.exp_mul]
  have hY : Ymean N τ x = Real.exp (-(τ * L / N)) ^ 2 * (∑ i, x i ^ (2 * τ)) / N := by
    unfold Ymean
    rw [← hL, Finset.mul_sum]
    congr 1
    refine Finset.sum_congr rfl (fun i _ => ?_)
    rw [Real.rpow_def_of_pos (hx i), ← Real.exp_nat_mul, ← Real.exp_add]
    congr 1
    push_cast
    ring
  have hS : 0 ≤ ∑ i, x i ^ (2 * τ) := Finset.sum_nonneg (fun i _ => (Real.rpow_pos_of_pos (hx i) _).le)
  have hsN : 0 < Real.sqrt N := Real.sqrt_pos.mpr (by positivity)
  rw [hPa, hY, Real.sqrt_div' _ (Nat.cast_nonneg N), Real.sqrt_mul (sq_nonneg _),
    Real.sqrt_sq (Real.exp_pos _).le, hP, ← Real.exp_mul]
  have hE : Real.exp (L * aPar N τ) * Real.exp (-(τ * L / N)) = Real.exp (L * -τ) := by
    rw [← Real.exp_add]; congr 1; unfold aPar; field_simp; ring
  rw [← hE]
  field_simp

end R03

theorem R03.proof_5 : R03.stmt_5 := by
  intro N hN p hp
  obtain ⟨k, rfl⟩ : ∃ k, N = k + 1 + 1 := ⟨N - 2, by omega⟩
  set n := k + 1 with hn
  set τ := 1 - 1 / p with hτ
  have hτ1 : τ < 1 := by have : 0 < 1 / p := by positivity
                         linarith
  set a := R03.aPar (n + 1) τ with ha
  have ha1 : -1 < a := by
    rw [ha, R03.aPar]
    have h1 : (0 : ℝ) < ((n + 1 : ℕ) : ℝ) := by positivity
    rw [lt_div_iff₀ h1]
    push_cast
    nlinarith
  have hNpos : 0 < n + 1 := by omega
  set F : (Fin (n + 1) → ℝ) → ℝ := fun g => Real.sqrt (R03.Ymean (n + 1) τ g) with hFdef
  have hFm : Measurable F := (R03.Ymean_meas _ _).sqrt
  have hFh : ∀ c : ℝ, 0 < c → ∀ g ∈ R03.orth (n + 1), F (c • g) = F g := by
    intro c hc g hg; simp only [hFdef, R03.Ymean_hom hNpos τ hc hg]
  have hD := R03.dirichlet n (R03.simplex_ae k) a ha1 F hFm hFh
  -- the simplex side
  set G : (Fin n → ℝ) → ℝ := fun t => (∏ i, R03.ext t i ^ a) * F (R03.ext t) with hG
  have hface : ∀ t ∈ R03.opS n, R03.faceIntegrand n τ t = Real.sqrt ((n + 1 : ℕ) : ℝ) * G t := by
    intro t ht
    exact R03.face_key hNpos τ ht
  have hint1 : IntegrableOn (R03.faceIntegrand n τ) (R03.simplex n) :=
    R03.proof_1 (n + 1) (by omega) p hp
  have hsN : 0 < Real.sqrt ((n + 1 : ℕ) : ℝ) := Real.sqrt_pos.mpr (by positivity)
  have hGint : IntegrableOn G (R03.simplex n) := by
    have : IntegrableOn (fun t => (Real.sqrt ((n + 1 : ℕ) : ℝ))⁻¹ * R03.faceIntegrand n τ t)
        (R03.opS n) := (hint1.mono_set (R03.opS_sub n)).const_mul _
    refine (this.congr_fun (fun t ht => ?_) (R03.opS_meas n)).congr_set_ae (R03.simplex_ae k).symm
    rw [hface t ht]; field_simp
  have hG0 : ∀ t ∈ R03.simplex n, 0 ≤ G t := by
    intro t ht
    have hb := R03.ext_bd ht
    exact mul_nonneg (Finset.prod_nonneg (fun i _ => Real.rpow_nonneg (hb i).1 _))
      (Real.sqrt_nonneg _)
  have hGlint : ∫⁻ t in R03.simplex n, ENNReal.ofReal (G t) =
      ENNReal.ofReal (∫ t in R03.simplex n, G t) :=
    (ofReal_integral_eq_lintegral_ofReal hGint
      (ae_restrict_of_forall_mem (R03.simplex_meas n) hG0)).symm
  have hfaceint : ∫ t in R03.simplex n, R03.faceIntegrand n τ t =
      Real.sqrt ((n + 1 : ℕ) : ℝ) * ∫ t in R03.simplex n, G t := by
    rw [← setIntegral_congr_set (R03.simplex_ae k), ← setIntegral_congr_set (R03.simplex_ae k),
      ← integral_const_mul]
    exact setIntegral_congr_fun (R03.opS_meas n) hface
  -- the orthant side
  have hdF0 : ∀ g ∈ R03.orth (n + 1), 0 ≤ R03.dens (n + 1) τ g * F g := by
    intro g hg
    have hgi : ∀ i, 0 < g i := fun i => hg i (mem_univ i)
    exact mul_nonneg (mul_nonneg (Finset.prod_nonneg (fun i _ => (Real.rpow_pos_of_pos (hgi i) _).le))
      (Real.exp_pos _).le) (Real.sqrt_nonneg _)
  have hdFm : Measurable fun g => R03.dens (n + 1) τ g * F g := by
    have h1 : Measurable fun g : Fin (n + 1) → ℝ => ∏ i, g i ^ a :=
      Finset.measurable_prod (f := fun i (g : Fin (n + 1) → ℝ) => g i ^ a) _
        (fun i _ => (measurable_pi_apply i).pow_const a)
    have h2 : Measurable fun g : Fin (n + 1) → ℝ => Real.exp (-∑ i, g i) :=
      Real.measurable_exp.comp (Finset.measurable_sum (f := fun i (g : Fin (n + 1) → ℝ) => g i) _
        (fun i _ => measurable_pi_apply i)).neg
    exact (h1.mul h2).mul hFm
  have hlint : ∫⁻ g in R03.orth (n + 1), ENNReal.ofReal (R03.dens (n + 1) τ g * F g) =
      ENNReal.ofReal (Real.Gamma (((n : ℝ) + 1) * (1 + a)) * ∫ t in R03.simplex n, G t) := by
    rw [ENNReal.ofReal_mul (Real.Gamma_nonneg_of_nonneg (by
      have : (0 : ℝ) < (n : ℝ) + 1 := by positivity
      have : 0 < 1 + a := by linarith
      positivity)), ← hGlint]
    exact hD
  have hint : IntegrableOn (fun g => R03.dens (n + 1) τ g * F g) (R03.orth (n + 1)) := by
    refine ⟨hdFm.aestronglyMeasurable, ?_⟩
    rw [hasFiniteIntegral_iff_ofReal (ae_restrict_of_forall_mem (R03.orth_meas _) hdF0), hlint]
    exact ENNReal.ofReal_lt_top
  refine ⟨hint, ?_⟩
  have hI : ∫ g in R03.orth (n + 1), R03.dens (n + 1) τ g * F g =
      Real.Gamma (((n : ℝ) + 1) * (1 + a)) * ∫ t in R03.simplex n, G t := by
    rw [integral_eq_lintegral_of_nonneg_ae (ae_restrict_of_forall_mem (R03.orth_meas _) hdF0)
      hdFm.aestronglyMeasurable, hlint, ENNReal.toReal_ofReal]
    have : (0 : ℝ) < (n : ℝ) + 1 := by positivity
    have : 0 < 1 + a := by linarith
    exact mul_nonneg (Real.Gamma_nonneg_of_nonneg (by positivity))
      (setIntegral_nonneg (R03.simplex_meas n) hG0)
  have hΓa : 0 < Real.Gamma (1 + a) := Real.Gamma_pos_of_pos (by linarith)
  have hNa : ((n : ℝ) + 1) * (1 + a) = ((n + 1 : ℕ) : ℝ) - (((n + 1 : ℕ) : ℝ) - 1) * τ := by
    rw [ha, R03.aPar]; push_cast; field_simp; ring
  have hΓN : 0 < Real.Gamma (((n : ℝ) + 1) * (1 + a)) := Real.Gamma_pos_of_pos (by
    have : (0 : ℝ) < (n : ℝ) + 1 := by positivity
    have : 0 < 1 + a := by linarith
    positivity)
  unfold R03.B R03.gammaE
  rw [← ha, ← hNa]
  change R03.S (n + 1) p = _ * (_ * ∫ g in R03.orth (n + 1), R03.dens (n + 1) τ g * F g)
  rw [hI, R03.S]
  simp only [Nat.add_sub_cancel]
  rw [← hτ]
  change 2 ^ (n + 1) * (1 / p) ^ n * ∫ t in R03.simplex n, R03.faceIntegrand n τ t = _
  rw [hfaceint]
  have hΓN' : Real.Gamma ((1 + a) * ((n : ℝ) + 1)) ≠ 0 := by rw [mul_comm]; exact hΓN.ne'
  field_simp


namespace R03

lemma choose_half (k : ℕ) : Ring.choose (1 / 2 : ℝ) k = binomHalf k := by
  have h := Ring.descPochhammer_eq_factorial_smul_choose (1 / 2 : ℝ) k
  rw [← Polynomial.aeval_eq_smeval, Polynomial.aeval_def, Polynomial.eval₂_eq_eval_map,
    descPochhammer_map, descPochhammer_eval_eq_prod_range, nsmul_eq_mul] at h
  rw [binomHalf, h]; field_simp

lemma pow_le_of_ge {δ x : ℝ} (hδ0 : 0 ≤ δ) (hδ1 : δ ≤ 1) (hx : δ ≤ x) {j J : ℕ} (hj : j ≤ J) :
    x ^ j * δ ^ J ≤ x ^ J := by
  have hx0 : 0 ≤ x := hδ0.trans hx
  calc x ^ j * δ ^ J ≤ x ^ j * δ ^ (J - j) :=
        mul_le_mul_of_nonneg_left (pow_le_pow_of_le_one hδ0 hδ1 (Nat.sub_le J j)) (pow_nonneg hx0 _)
    _ ≤ x ^ j * x ^ (J - j) := by gcongr
    _ = x ^ J := by rw [← pow_add, Nat.add_sub_cancel' hj]

/-- the binomial remainder for `√(1+x)`, `x ≥ -1` (proof of `lem:pw:X`). -/
lemma sqrt_rem (J : ℕ) : ∃ C : ℝ, 0 ≤ C ∧ ∀ x : ℝ, -1 ≤ x →
    |Real.sqrt (1 + x) - ∑ j ∈ Finset.range (J + 1), binomHalf j * x ^ j| ≤ C * |x| ^ (J + 1) := by
  have hB := (Real.one_add_rpow_hasFPowerSeriesOnBall_zero (a := 1 / 2)).hasFPowerSeriesAt.isBigO_sub_partialSum_pow (J + 1)
  simp only [zero_add, FormalMultilinearSeries.partialSum, binomialSeries_apply] at hB
  obtain ⟨C0, hC0⟩ := hB.bound
  rw [Metric.eventually_nhds_iff] at hC0
  obtain ⟨ε, hε, hεb⟩ := hC0
  set δ := min ε 1 with hδ
  have hδ0 : 0 < δ := lt_min hε one_pos
  have hδ1 : δ ≤ 1 := min_le_right _ _
  set C1 := (2 + ∑ j ∈ Finset.range (J + 1), |binomHalf j|) / δ ^ (J + 1) with hC1
  have hC1n : 0 ≤ C1 := by
    have : 0 ≤ ∑ j ∈ Finset.range (J + 1), |binomHalf j| := Finset.sum_nonneg (fun j _ => abs_nonneg _)
    positivity
  refine ⟨max C0 C1, le_max_of_le_right hC1n, fun x hx => ?_⟩
  by_cases hsmall : |x| < ε
  · have := hεb (y := x) (by simpa [Real.dist_eq] using hsmall)
    simp only [List.prod_ofFn, Fin.prod_const, smul_eq_mul, choose_half, Real.norm_eq_abs,
      abs_pow, abs_abs] at this
    rw [Real.sqrt_eq_rpow]
    calc _ ≤ C0 * |x| ^ (J + 1) := this
      _ ≤ max C0 C1 * |x| ^ (J + 1) := by gcongr; exact le_max_left _ _
  · have hxδ : δ ≤ |x| := le_trans (min_le_left _ _) (not_lt.mp hsmall)
    have hsq : Real.sqrt (1 + x) ≤ 1 + |x| := by
      rw [Real.sqrt_le_left (by positivity)]
      nlinarith [le_abs_self x, abs_nonneg x]
    have hP : |∑ j ∈ Finset.range (J + 1), binomHalf j * x ^ j| ≤
        ∑ j ∈ Finset.range (J + 1), |binomHalf j| * |x| ^ j := by
      refine (Finset.abs_sum_le_sum_abs _ _).trans (le_of_eq ?_)
      exact Finset.sum_congr rfl (fun j _ => by rw [abs_mul, abs_pow])
    have hpow : 0 < δ ^ (J + 1) := by positivity
    have h0 : 1 * δ ^ (J + 1) ≤ |x| ^ (J + 1) := by
      simpa using pow_le_of_ge hδ0.le hδ1 hxδ (j := 0) (J := J + 1) (by omega)
    have h1 : |x| * δ ^ (J + 1) ≤ |x| ^ (J + 1) := by
      simpa using pow_le_of_ge hδ0.le hδ1 hxδ (j := 1) (J := J + 1) (by omega)
    have h2 : (∑ j ∈ Finset.range (J + 1), |binomHalf j| * |x| ^ j) * δ ^ (J + 1) ≤
        (∑ j ∈ Finset.range (J + 1), |binomHalf j|) * |x| ^ (J + 1) := by
      rw [Finset.sum_mul, Finset.sum_mul]
      refine Finset.sum_le_sum (fun j hj => ?_)
      rw [mul_assoc]
      gcongr
      exact pow_le_of_ge hδ0.le hδ1 hxδ (by simp at hj; omega)
    have hmain : |Real.sqrt (1 + x) - ∑ j ∈ Finset.range (J + 1), binomHalf j * x ^ j| ≤
        C1 * |x| ^ (J + 1) := by
      rw [hC1, div_mul_eq_mul_div, le_div_iff₀ hpow]
      calc |Real.sqrt (1 + x) - ∑ j ∈ Finset.range (J + 1), binomHalf j * x ^ j| * δ ^ (J + 1)
          ≤ ((1 + |x|) + ∑ j ∈ Finset.range (J + 1), |binomHalf j| * |x| ^ j) * δ ^ (J + 1) := by
            gcongr
            refine (abs_sub _ _).trans ?_
            rw [abs_of_nonneg (Real.sqrt_nonneg _)]
            linarith
        _ ≤ _ := by nlinarith
    calc _ ≤ C1 * |x| ^ (J + 1) := hmain
      _ ≤ max C0 C1 * |x| ^ (J + 1) := by gcongr; exact le_max_right _ _

lemma exp_sub_one_abs (u : ℝ) : |Real.exp u - 1| ≤ |u| * Real.exp |u| := by
  rcases le_or_gt 0 u with hu | hu
  · rw [abs_of_nonneg hu, abs_of_nonneg (by linarith [Real.add_one_le_exp u])]
    have h1 := Real.add_one_le_exp (-u)
    have h2 : Real.exp u * Real.exp (-u) = 1 := by rw [← Real.exp_add]; simp
    nlinarith [Real.exp_pos u]
  · rw [abs_of_neg hu, abs_of_nonpos (by linarith [(Real.exp_le_exp.mpr hu.le).trans_eq Real.exp_zero])]
    have h1 := Real.add_one_le_exp u
    have h2 : 1 ≤ Real.exp (-u) := Real.one_le_exp (by linarith)
    nlinarith

/-- `Ξ = Σ_j |log g_j|`. -/
noncomputable def Xi {N : ℕ} (g : Fin N → ℝ) : ℝ := ∑ j, |Real.log (g j)|

lemma Xi_nn {N : ℕ} (g : Fin N → ℝ) : 0 ≤ Xi g := Finset.sum_nonneg (fun _ _ => abs_nonneg _)

lemma X_bd {N : ℕ} (hN : 1 ≤ N) (τ : ℝ) (g : Fin N → ℝ) :
    |Ymean N τ g - 1| ≤ 4 * |τ| * Xi g * Real.exp (4 * |τ| * Xi g) := by
  have hN' : (0 : ℝ) < N := by exact_mod_cast hN
  have hN1 : (1 : ℝ) ≤ N := by exact_mod_cast hN
  set L := ∑ j, Real.log (g j)
  have hL : |L| ≤ Xi g := Finset.abs_sum_le_sum_abs _ _
  have hu : ∀ i, |2 * τ * (Real.log (g i) - L / N)| ≤ 4 * |τ| * Xi g := by
    intro i
    have hi : |Real.log (g i)| ≤ Xi g :=
      Finset.single_le_sum (f := fun j => |Real.log (g j)|) (fun j _ => abs_nonneg _)
        (Finset.mem_univ i)
    have hLN : |L / N| ≤ Xi g := by
      rw [abs_div, abs_of_pos hN', div_le_iff₀ hN']; nlinarith [abs_nonneg L, Xi_nn g]
    rw [abs_mul, abs_mul, abs_two]
    have := abs_sub (Real.log (g i)) (L / N)
    have : 0 ≤ |τ| := abs_nonneg τ
    nlinarith
  have hmono : ∀ i, |2 * τ * (Real.log (g i) - L / N)| *
      Real.exp |2 * τ * (Real.log (g i) - L / N)| ≤ 4 * |τ| * Xi g * Real.exp (4 * |τ| * Xi g) :=
    fun i => mul_le_mul (hu i) (Real.exp_le_exp.mpr (hu i)) (Real.exp_pos _).le
      (by have := Xi_nn g; positivity)
  have hY : Ymean N τ g - 1 = (∑ i, (Real.exp (2 * τ * (Real.log (g i) - L / N)) - 1)) / N := by
    unfold Ymean
    rw [Finset.sum_sub_distrib]
    simp only [Finset.sum_const, Finset.card_univ, Fintype.card_fin, nsmul_eq_mul, mul_one]
    rw [sub_div, div_self hN'.ne']
  rw [hY, abs_div, abs_of_pos hN', div_le_iff₀ hN']
  calc |∑ i, (Real.exp (2 * τ * (Real.log (g i) - L / N)) - 1)|
      ≤ ∑ i, |Real.exp (2 * τ * (Real.log (g i) - L / N)) - 1| := Finset.abs_sum_le_sum_abs _ _
    _ ≤ ∑ _i : Fin N, 4 * |τ| * Xi g * Real.exp (4 * |τ| * Xi g) :=
        Finset.sum_le_sum (fun i _ => (exp_sub_one_abs _).trans (hmono i))
    _ = _ := by simp; ring

lemma Xpow_bd {N : ℕ} (hN : 1 ≤ N) (J : ℕ) {τ : ℝ} (hτ : |τ| ≤ 1 / (16 * ((J : ℝ) + 1)))
    (g : Fin N → ℝ) :
    |Ymean N τ g - 1| ^ (J + 1) ≤
      (4 * |τ|) ^ (J + 1) * (4 ^ (J + 1) * ((J + 1).factorial : ℝ)) * Real.exp (Xi g / 2) := by
  have hX := Xi_nn g
  have hJ : (0 : ℝ) < (J : ℝ) + 1 := by positivity
  have h1 := pow_le_pow_left₀ (abs_nonneg _) (X_bd hN τ g) (J + 1)
  have h2 : (Xi g / 4) ^ (J + 1) / ((J + 1).factorial : ℝ) ≤ Real.exp (Xi g / 4) :=
    Real.pow_div_factorial_le_exp (Xi g / 4) (by positivity) (J + 1)
  have h2' : Xi g ^ (J + 1) ≤ 4 ^ (J + 1) * ((J + 1).factorial : ℝ) * Real.exp (Xi g / 4) := by
    rw [div_pow, div_div, div_le_iff₀ (by positivity)] at h2
    nlinarith
  have h3 : Real.exp (4 * |τ| * Xi g) ^ (J + 1) ≤ Real.exp (Xi g / 4) := by
    rw [← Real.exp_nat_mul, Real.exp_le_exp]
    have : ((J + 1 : ℕ) : ℝ) * (4 * |τ|) ≤ 1 / 4 := by
      push_cast
      rw [le_div_iff₀ (by positivity)] at hτ
      nlinarith [abs_nonneg τ]
    nlinarith
  have he : Real.exp (Xi g / 2) = Real.exp (Xi g / 4) * Real.exp (Xi g / 4) := by
    rw [← Real.exp_add]; ring_nf
  calc |Ymean N τ g - 1| ^ (J + 1) ≤ (4 * |τ| * Xi g * Real.exp (4 * |τ| * Xi g)) ^ (J + 1) := h1
    _ = (4 * |τ|) ^ (J + 1) * Xi g ^ (J + 1) * Real.exp (4 * |τ| * Xi g) ^ (J + 1) := by
        rw [mul_pow, mul_pow]
    _ ≤ (4 * |τ|) ^ (J + 1) * (4 ^ (J + 1) * ((J + 1).factorial : ℝ) * Real.exp (Xi g / 4)) *
          Real.exp (Xi g / 4) := by gcongr
    _ = _ := by rw [he]; ring

lemma exp_half_log {x : ℝ} (hx : 0 < x) :
    Real.exp (|Real.log x| / 2) ≤ x ^ (1 / 2 : ℝ) + x ^ (-(1 / 2) : ℝ) := by
  rw [Real.rpow_def_of_pos hx, Real.rpow_def_of_pos hx]
  rcases le_or_gt 0 (Real.log x) with h | h
  · rw [abs_of_nonneg h]
    have : Real.exp (|Real.log x| / 2) = Real.exp (Real.log x * (1 / 2)) := by
      rw [abs_of_nonneg h]; ring_nf
    rw [abs_of_nonneg h] at this
    rw [this]; linarith [Real.exp_pos (Real.log x * -(1 / 2))]
  · rw [abs_of_neg h]
    have : Real.exp (-Real.log x / 2) = Real.exp (Real.log x * -(1 / 2)) := by ring_nf
    rw [this]; linarith [Real.exp_pos (Real.log x * (1 / 2))]

end R03

namespace R03

lemma orth_ind {N : ℕ} (φ : ℝ → ℝ) :
    (orth N).indicator (fun g : Fin N → ℝ => ∏ j, φ (g j)) =
      fun g => ∏ j, (Ioi (0 : ℝ)).indicator φ (g j) := by
  funext g
  by_cases hg : g ∈ orth N
  · rw [indicator_of_mem hg]
    exact Finset.prod_congr rfl (fun j _ => (indicator_of_mem (hg j (mem_univ j)) _).symm)
  · rw [indicator_of_notMem hg]
    have : ∃ j, g j ∉ Ioi (0 : ℝ) := by
      by_contra h; push Not at h; exact hg (fun j _ => h j)
    obtain ⟨j, hj⟩ := this
    exact (Finset.prod_eq_zero (Finset.mem_univ j) (indicator_of_notMem hj _)).symm

lemma prod_intOn {N : ℕ} {φ : ℝ → ℝ} (h : IntegrableOn φ (Ioi 0)) :
    IntegrableOn (fun g : Fin N → ℝ => ∏ j, φ (g j)) (orth N) := by
  rw [← integrable_indicator_iff (orth_meas N), orth_ind]
  exact Integrable.fintype_prod (f := fun _ => (Ioi (0 : ℝ)).indicator φ)
    (fun _ => h.integrable_indicator measurableSet_Ioi)

lemma prod_integral {N : ℕ} (φ : ℝ → ℝ) :
    ∫ g in orth N, ∏ j, φ (g j) = (∫ x in Ioi 0, φ x) ^ N := by
  rw [← integral_indicator (orth_meas N), orth_ind, volume_pi,
    integral_fintype_prod_eq_prod (fun _ => (Ioi (0 : ℝ)).indicator φ),
    integral_indicator measurableSet_Ioi, Finset.prod_const, Finset.card_univ, Fintype.card_fin]

lemma dens_eq (N : ℕ) (τ : ℝ) (g : Fin N → ℝ) :
    dens N τ g = ∏ j, (g j ^ aPar N τ * Real.exp (-g j)) := by
  rw [dens, Finset.prod_mul_distrib, ← Real.exp_sum, Finset.sum_neg_distrib]

lemma dens_meas (N : ℕ) (τ : ℝ) : Measurable (dens N τ) := by
  have h1 : Measurable fun g : Fin N → ℝ => ∏ i, g i ^ aPar N τ :=
    Finset.measurable_prod (f := fun i (g : Fin N → ℝ) => g i ^ aPar N τ) _
      (fun i _ => (measurable_pi_apply i).pow_const _)
  have h2 : Measurable fun g : Fin N → ℝ => Real.exp (-∑ i, g i) :=
    Real.measurable_exp.comp (Finset.measurable_sum (f := fun i (g : Fin N → ℝ) => g i) _
      (fun i _ => measurable_pi_apply i)).neg
  exact h1.mul h2

lemma dens_nn {N : ℕ} (τ : ℝ) {g : Fin N → ℝ} (hg : g ∈ orth N) : 0 ≤ dens N τ g :=
  mul_nonneg (Finset.prod_nonneg (fun i _ => (Real.rpow_pos_of_pos (hg i (mem_univ i)) _).le))
    (Real.exp_pos _).le

lemma gam_int {a : ℝ} (ha : -1 < a) :
    IntegrableOn (fun x : ℝ => x ^ a * Real.exp (-x)) (Ioi 0) ∧
      ∫ x in Ioi 0, x ^ a * Real.exp (-x) = Real.Gamma (1 + a) := by
  have hs : 0 < 1 + a := by linarith
  have he : EqOn (fun x : ℝ => Real.exp (-x) * x ^ (1 + a - 1)) (fun x => x ^ a * Real.exp (-x))
      (Ioi 0) := fun x _ => by simp only; rw [show 1 + a - 1 = a by ring, mul_comm]
  refine ⟨(Real.GammaIntegral_convergent hs).congr_fun he measurableSet_Ioi, ?_⟩
  rw [Real.Gamma_eq_integral hs]
  exact (setIntegral_congr_fun measurableSet_Ioi he).symm

/-- `x^a e^{-x}(x^{1/2}+x^{-1/2})` integrates to `Γ(a+3/2)+Γ(a+1/2)`. -/
noncomputable def phiH (a x : ℝ) : ℝ := x ^ a * Real.exp (-x) * (x ^ (1 / 2 : ℝ) + x ^ (-(1 / 2) : ℝ))

lemma phiH_int {a : ℝ} (ha : -(1 / 2) < a) :
    IntegrableOn (phiH a) (Ioi 0) ∧
      ∫ x in Ioi 0, phiH a x = Real.Gamma (a + 3 / 2) + Real.Gamma (a + 1 / 2) := by
  have h1 := gam_int (a := a + 1 / 2) (by linarith)
  have h2 := gam_int (a := a - 1 / 2) (by linarith)
  have he : EqOn (fun x => x ^ (a + 1 / 2) * Real.exp (-x) + x ^ (a - 1 / 2) * Real.exp (-x))
      (phiH a) (Ioi 0) := by
    intro x hx
    have hx0 : (0 : ℝ) < x := hx
    simp only [phiH]
    rw [Real.rpow_add hx0, show a - 1 / 2 = a + -(1 / 2) by ring, Real.rpow_add hx0]
    ring
  refine ⟨(h1.1.add h2.1).congr_fun he measurableSet_Ioi, ?_⟩
  rw [← setIntegral_congr_fun measurableSet_Ioi he, integral_add h1.1 h2.1, h1.2, h2.2]
  ring_nf

lemma gamma_bds : ∃ m M : ℝ, 0 < m ∧ ∀ x ∈ Icc (1 / 4 : ℝ) 2, m ≤ Real.Gamma x ∧ Real.Gamma x ≤ M := by
  have hc : ContinuousOn Real.Gamma (Icc (1 / 4 : ℝ) 2) := by
    intro x hx
    have hx0 : 0 < x := lt_of_lt_of_le (by norm_num) hx.1
    refine (Real.differentiableAt_Gamma (fun m => ?_)).continuousAt.continuousWithinAt
    have : (0 : ℝ) ≤ m := Nat.cast_nonneg m
    linarith
  have hne : (Icc (1 / 4 : ℝ) 2).Nonempty := ⟨1, by norm_num, by norm_num⟩
  obtain ⟨x0, hx0, hmin⟩ := isCompact_Icc.exists_isMinOn hne hc
  obtain ⟨x1, -, hmax⟩ := isCompact_Icc.exists_isMaxOn hne hc
  refine ⟨Real.Gamma x0, Real.Gamma x1, Real.Gamma_pos_of_pos (lt_of_lt_of_le (by norm_num) hx0.1),
    fun x hx => ⟨hmin hx, hmax hx⟩⟩

lemma aPar_abs {N : ℕ} (hN : 1 ≤ N) (τ : ℝ) : |aPar N τ| ≤ |τ| := by
  have hN' : (0 : ℝ) < N := by exact_mod_cast hN
  have hN1 : (1 : ℝ) ≤ N := by exact_mod_cast hN
  rw [aPar, abs_div, abs_neg, abs_mul, abs_of_pos hN', abs_of_nonneg (by linarith),
    div_le_iff₀ hN']
  nlinarith [abs_nonneg τ]

end R03

theorem R03.proof_7 : R03.stmt_7 := by
  intro N hN J
  have hN1 : 1 ≤ N := by omega
  obtain ⟨CJ, hCJ0, hCJ⟩ := R03.sqrt_rem J
  obtain ⟨m, M, hm, hmM⟩ := R03.gamma_bds
  have hMpos : 0 < M := lt_of_lt_of_le hm ((hmM 1 ⟨by norm_num, by norm_num⟩).1.trans
    (hmM 1 ⟨by norm_num, by norm_num⟩).2)
  set K : ℝ := 4 ^ (J + 1) * ((J + 1).factorial : ℝ) with hK
  have hK0 : 0 < K := by positivity
  set C0 : ℝ := 4 ^ (J + 1) * K * (2 * M) ^ N / m ^ N with hC0
  have hC00 : 0 ≤ C0 := by positivity
  refine ⟨(1 + CJ) * C0, 1 / (16 * ((J : ℝ) + 1)), by positivity, fun τ hτ => ?_⟩
  have hτ16 : |τ| ≤ 1 / 16 := by
    refine hτ.trans ?_
    rw [div_le_div_iff₀ (by positivity) (by norm_num)]
    have : (0 : ℝ) ≤ J := Nat.cast_nonneg J
    linarith
  set a := R03.aPar N τ with ha
  have hab := R03.aPar_abs hN1 τ
  rw [← ha] at hab
  have hab' := abs_le.mp (hab.trans hτ16)
  -- the Gamma factors
  have hΓa := hmM (1 + a) ⟨by linarith, by linarith⟩
  have hΓ3 := hmM (a + 3 / 2) ⟨by linarith, by linarith⟩
  have hΓ1 := hmM (a + 1 / 2) ⟨by linarith, by linarith⟩
  have hΓpos : 0 < Real.Gamma (1 + a) := lt_of_lt_of_le hm hΓa.1
  -- the integrable majorants
  obtain ⟨hψi, -⟩ := R03.gam_int (a := a) (by linarith)
  obtain ⟨hφi, hφv⟩ := R03.phiH_int (a := a) (by linarith)
  have hdi : IntegrableOn (R03.dens N τ) (R03.orth N) := by
    have := R03.prod_intOn (N := N) hψi
    refine this.congr_fun (fun g _ => ?_) (R03.orth_meas N)
    rw [R03.dens_eq]
  set D : (Fin N → ℝ) → ℝ := fun g => (4 * |τ|) ^ (J + 1) * K * ∏ j, R03.phiH a (g j) with hD
  have hDi : IntegrableOn D (R03.orth N) := (R03.prod_intOn hφi).const_mul _
  have hDv : ∫ g in R03.orth N, D g = (4 * |τ|) ^ (J + 1) * K *
      (Real.Gamma (a + 3 / 2) + Real.Gamma (a + 1 / 2)) ^ N := by
    rw [hD, integral_const_mul, R03.prod_integral, hφv]
  have hpt : ∀ g ∈ R03.orth N, R03.dens N τ g * |R03.Ymean N τ g - 1| ^ (J + 1) ≤ D g := by
    intro g hg
    have hgi : ∀ i, 0 < g i := fun i => hg i (mem_univ i)
    have h1 := R03.Xpow_bd hN1 J hτ g
    have h2 : Real.exp (R03.Xi g / 2) ≤ ∏ j, (g j ^ (1 / 2 : ℝ) + g j ^ (-(1 / 2) : ℝ)) := by
      rw [R03.Xi, Finset.sum_div, Real.exp_sum]
      exact Finset.prod_le_prod₀ (fun j _ => (Real.exp_pos _).le)
        (fun j _ => R03.exp_half_log (hgi j))
    have h3 : R03.dens N τ g * ∏ j, (g j ^ (1 / 2 : ℝ) + g j ^ (-(1 / 2) : ℝ)) =
        ∏ j, R03.phiH a (g j) := by
      rw [R03.dens_eq, ← Finset.prod_mul_distrib]; rfl
    have hd0 := R03.dens_nn τ hg
    calc R03.dens N τ g * |R03.Ymean N τ g - 1| ^ (J + 1)
        ≤ R03.dens N τ g * ((4 * |τ|) ^ (J + 1) * K * Real.exp (R03.Xi g / 2)) := by gcongr
      _ ≤ R03.dens N τ g * ((4 * |τ|) ^ (J + 1) * K *
            ∏ j, (g j ^ (1 / 2 : ℝ) + g j ^ (-(1 / 2) : ℝ))) := by gcongr
      _ = D g := by simp only [hD]; rw [← h3]; ring
  have hD0 : ∀ g ∈ R03.orth N, 0 ≤ D g := fun g hg =>
    (mul_nonneg (R03.dens_nn τ hg) (by positivity)).trans (hpt g hg)
  have hYm := R03.Ymean_meas N τ
  have hdm := R03.dens_meas N τ
  have hpowi : ∀ j ≤ J + 1, IntegrableOn
      (fun g => R03.dens N τ g * (R03.Ymean N τ g - 1) ^ j) (R03.orth N) := by
    intro j hj
    refine Integrable.mono' (hdi.add hDi) ((hdm.mul ((hYm.sub_const 1).pow_const j)).aestronglyMeasurable)
      (ae_restrict_of_forall_mem (R03.orth_meas N) (fun g hg => ?_))
    have hd0 := R03.dens_nn τ hg
    rw [Real.norm_eq_abs, abs_mul, abs_of_nonneg hd0, abs_pow]
    have : |R03.Ymean N τ g - 1| ^ j ≤ 1 + |R03.Ymean N τ g - 1| ^ (J + 1) := by
      rcases le_or_gt |R03.Ymean N τ g - 1| 1 with h | h
      · have := pow_le_one₀ (abs_nonneg _) h (n := j)
        have := pow_nonneg (abs_nonneg (R03.Ymean N τ g - 1)) (J + 1)
        linarith
      · have := pow_le_pow_right₀ h.le hj
        linarith
    calc R03.dens N τ g * |R03.Ymean N τ g - 1| ^ j
        ≤ R03.dens N τ g * (1 + |R03.Ymean N τ g - 1| ^ (J + 1)) := by gcongr
      _ ≤ R03.dens N τ g + D g := by have := hpt g hg; nlinarith
  have habsi : IntegrableOn (fun g => R03.dens N τ g * |R03.Ymean N τ g - 1| ^ (J + 1))
      (R03.orth N) := by
    refine Integrable.mono' hDi ((hdm.mul ((hYm.sub_const 1).abs.pow_const _)).aestronglyMeasurable)
      (ae_restrict_of_forall_mem (R03.orth_meas N) (fun g hg => ?_))
    rw [Real.norm_eq_abs, abs_of_nonneg (mul_nonneg (R03.dens_nn τ hg) (by positivity))]
    exact hpt g hg
  -- the moment bound
  have hE : R03.gammaE N τ (fun g => |R03.Ymean N τ g - 1| ^ (J + 1)) ≤ C0 * |τ| ^ (J + 1) := by
    unfold R03.gammaE
    rw [← ha]
    have hI : ∫ g in R03.orth N, R03.dens N τ g * |R03.Ymean N τ g - 1| ^ (J + 1) ≤
        (4 * |τ|) ^ (J + 1) * K * (2 * M) ^ N := by
      refine (setIntegral_mono_on habsi hDi (R03.orth_meas N) hpt).trans ?_
      rw [hDv]
      gcongr
      · linarith [hΓ3.1, hm]
      · linarith [hΓ3.2, hΓ1.2]
    have hI0 : 0 ≤ ∫ g in R03.orth N, R03.dens N τ g * |R03.Ymean N τ g - 1| ^ (J + 1) :=
      setIntegral_nonneg (R03.orth_meas N) (fun g hg => mul_nonneg (R03.dens_nn τ hg) (by positivity))
    have hinv : (Real.Gamma (1 + a) ^ N)⁻¹ ≤ (m ^ N)⁻¹ :=
      inv_anti₀ (by positivity) (pow_le_pow_left₀ hm.le hΓa.1 N)
    calc (Real.Gamma (1 + a) ^ N)⁻¹ *
          ∫ g in R03.orth N, R03.dens N τ g * |R03.Ymean N τ g - 1| ^ (J + 1)
        ≤ (m ^ N)⁻¹ * ((4 * |τ|) ^ (J + 1) * K * (2 * M) ^ N) :=
          mul_le_mul hinv hI hI0 (by positivity)
      _ = C0 * |τ| ^ (J + 1) := by rw [hC0, mul_pow]; field_simp
  refine ⟨hpowi, ?_, ?_⟩
  · calc _ ≤ C0 * |τ| ^ (J + 1) := hE
      _ ≤ (1 + CJ) * C0 * |τ| ^ (J + 1) := by
        have : 0 ≤ C0 * |τ| ^ (J + 1) := by positivity
        nlinarith
  · -- the binomial remainder
    have hsqi : IntegrableOn (fun g => R03.dens N τ g * Real.sqrt (R03.Ymean N τ g)) (R03.orth N) := by
      refine Integrable.mono' ((hdi.const_mul 2).add hDi) ((hdm.mul hYm.sqrt).aestronglyMeasurable)
        (ae_restrict_of_forall_mem (R03.orth_meas N) (fun g hg => ?_))
      have hd0 := R03.dens_nn τ hg
      rw [Real.norm_eq_abs, abs_of_nonneg (mul_nonneg hd0 (Real.sqrt_nonneg _))]
      set X := R03.Ymean N τ g - 1
      have hsq : Real.sqrt (R03.Ymean N τ g) ≤ 1 + |X| := by
        rw [Real.sqrt_le_left (by positivity)]
        nlinarith [le_abs_self X, abs_nonneg X]
      have hX1 : |X| ≤ 1 + |X| ^ (J + 1) := by
        rcases le_or_gt |X| 1 with h | h
        · have := pow_nonneg (abs_nonneg X) (J + 1); linarith
        · have := pow_le_pow_right₀ h.le (show 1 ≤ J + 1 by omega)
          rw [pow_one] at this; linarith
      calc R03.dens N τ g * Real.sqrt (R03.Ymean N τ g) ≤ R03.dens N τ g * (2 + |X| ^ (J + 1)) := by
            gcongr; linarith
        _ ≤ 2 * R03.dens N τ g + D g := by have := hpt g hg; nlinarith
    have hlin : R03.B N τ - ∑ j ∈ Finset.range (J + 1), R03.binomHalf j * R03.momX N j τ =
        (Real.Gamma (1 + a) ^ N)⁻¹ * ∫ g in R03.orth N, (R03.dens N τ g * Real.sqrt (R03.Ymean N τ g) -
          ∑ j ∈ Finset.range (J + 1), R03.binomHalf j * (R03.dens N τ g * (R03.Ymean N τ g - 1) ^ j)) := by
      rw [integral_sub hsqi (integrable_finsetSum _ (fun j hj => (hpowi j (by
          simp at hj; omega)).const_mul _)),
        integral_finsetSum _ (fun j hj => (hpowi j (by simp at hj; omega)).const_mul _)]
      unfold R03.B R03.momX R03.gammaE
      rw [← ha, mul_sub, Finset.mul_sum]
      congr 1
      refine Finset.sum_congr rfl (fun j _ => ?_)
      rw [integral_const_mul]; ring
    rw [hlin, abs_mul, abs_of_pos (by positivity)]
    have hbd : ‖∫ g in R03.orth N, (R03.dens N τ g * Real.sqrt (R03.Ymean N τ g) -
          ∑ j ∈ Finset.range (J + 1), R03.binomHalf j * (R03.dens N τ g * (R03.Ymean N τ g - 1) ^ j))‖
        ≤ ∫ g in R03.orth N, CJ * (R03.dens N τ g * |R03.Ymean N τ g - 1| ^ (J + 1)) := by
      refine norm_integral_le_of_norm_le (habsi.const_mul CJ)
        (ae_restrict_of_forall_mem (R03.orth_meas N) (fun g hg => ?_))
      have hd0 := R03.dens_nn τ hg
      have hY0 : 0 ≤ R03.Ymean N τ g := by
        unfold R03.Ymean; exact div_nonneg (Finset.sum_nonneg (fun _ _ => (Real.exp_pos _).le))
          (Nat.cast_nonneg N)
      have := hCJ (R03.Ymean N τ g - 1) (by linarith)
      rw [show 1 + (R03.Ymean N τ g - 1) = R03.Ymean N τ g by ring] at this
      rw [Real.norm_eq_abs]
      calc |R03.dens N τ g * Real.sqrt (R03.Ymean N τ g) -
            ∑ j ∈ Finset.range (J + 1), R03.binomHalf j * (R03.dens N τ g * (R03.Ymean N τ g - 1) ^ j)|
          = R03.dens N τ g * |Real.sqrt (R03.Ymean N τ g) -
            ∑ j ∈ Finset.range (J + 1), R03.binomHalf j * (R03.Ymean N τ g - 1) ^ j| := by
            rw [← abs_of_nonneg hd0, ← abs_mul, abs_of_nonneg hd0, mul_sub, Finset.mul_sum]
            congr 2
            refine Finset.sum_congr rfl (fun j _ => ?_)
            ring
        _ ≤ R03.dens N τ g * (CJ * |R03.Ymean N τ g - 1| ^ (J + 1)) := by gcongr
        _ = _ := by ring
    rw [Real.norm_eq_abs, integral_const_mul] at hbd
    have hE' : (Real.Gamma (1 + a) ^ N)⁻¹ *
        ∫ g in R03.orth N, R03.dens N τ g * |R03.Ymean N τ g - 1| ^ (J + 1) ≤ C0 * |τ| ^ (J + 1) := by
      have := hE; unfold R03.gammaE at this; rwa [← ha] at this
    calc (Real.Gamma (1 + a) ^ N)⁻¹ * |∫ g in R03.orth N, (R03.dens N τ g * Real.sqrt (R03.Ymean N τ g) -
          ∑ j ∈ Finset.range (J + 1), R03.binomHalf j * (R03.dens N τ g * (R03.Ymean N τ g - 1) ^ j))|
        ≤ (Real.Gamma (1 + a) ^ N)⁻¹ *
          (CJ * ∫ g in R03.orth N, R03.dens N τ g * |R03.Ymean N τ g - 1| ^ (J + 1)) := by
          gcongr
      _ = CJ * ((Real.Gamma (1 + a) ^ N)⁻¹ *
          ∫ g in R03.orth N, R03.dens N τ g * |R03.Ymean N τ g - 1| ^ (J + 1)) := by ring
      _ ≤ CJ * (C0 * |τ| ^ (J + 1)) := by gcongr
      _ ≤ (1 + CJ) * C0 * |τ| ^ (J + 1) := by
        have : 0 ≤ C0 * |τ| ^ (J + 1) := by positivity
        nlinarith


namespace R03

open Filter Topology

lemma small_le {c K ε : ℝ} (hε : 0 < ε) (h : ∀ z : ℝ, 0 < z → z < ε → c ≤ K * z) : c ≤ 0 := by
  by_contra hc
  push Not at hc
  rcases le_or_gt K 0 with hK | hK
  · have := h (ε / 2) (by positivity) (by linarith)
    nlinarith
  · set z := min (ε / 2) (c / (2 * K))
    have hz0 : 0 < z := lt_min (by positivity) (by positivity)
    have hz1 : z < ε := lt_of_le_of_lt (min_le_left _ _) (by linarith)
    have := h z hz0 hz1
    have hz2 : K * z ≤ c / 2 := by
      have := min_le_right (ε / 2) (c / (2 * K))
      calc K * z ≤ K * (c / (2 * K)) := by gcongr
        _ = c / 2 := by field_simp
    linarith

/-- an analytic germ with `|f x| ≤ C|x|^j` has vanishing Taylor coefficients below `j`. -/
lemma coeff_vanish {f : ℝ → ℝ} (hf : AnalyticAt ℝ f 0) {j : ℕ} {C r : ℝ} (hr : 0 < r)
    (hb : ∀ x, |x| < r → |f x| ≤ C * |x| ^ j) : ∀ k < j, iteratedDeriv k f 0 = 0 := by
  intro k
  induction k using Nat.strong_induction_on with
  | _ k ih =>
  intro hk
  obtain ⟨F, hFa, hF⟩ := hf.exists_eventuallyEq_sum_add_pow_mul (k + 1)
  have hF' : ∀ᶠ z in 𝓝 (0 : ℝ), f z = z ^ k / (k.factorial : ℝ) * iteratedDeriv k f 0 +
      z ^ (k + 1) * F z := by
    filter_upwards [hF] with z hz
    rw [hz, Finset.sum_range_succ, Finset.sum_eq_zero (fun i hi => by
      rw [ih i (Finset.mem_range.mp hi) (by have := Finset.mem_range.mp hi; omega), smul_zero])]
    simp [smul_eq_mul]
  have hFb : ∀ᶠ z in 𝓝 (0 : ℝ), |F z| < |F 0| + 1 :=
    (hFa.continuousAt.abs).eventually_lt continuousAt_const (by linarith)
  have hball : ∀ᶠ z in 𝓝 (0 : ℝ), |z| < min r 1 := by
    have := Metric.ball_mem_nhds (0 : ℝ) (lt_min hr one_pos)
    filter_upwards [this] with z hz
    simpa [Real.dist_eq] using hz
  obtain ⟨ε, hε, hall⟩ := Metric.eventually_nhds_iff.mp (hF'.and (hFb.and hball))
  set d := iteratedDeriv k f 0
  have hk! : (0 : ℝ) < k.factorial := by positivity
  have key : |d| / k.factorial ≤ 0 := by
    refine small_le (K := |C| + |F 0| + 1) hε (fun z hz0 hzε => ?_)
    obtain ⟨h1, h2, h3⟩ := hall (y := z) (by simpa [Real.dist_eq, abs_of_pos hz0] using hzε)
    have hz1 : z < 1 := by have := lt_of_lt_of_le h3 (min_le_right _ _); rwa [abs_of_pos hz0] at this
    have hzr : |z| < r := lt_of_lt_of_le h3 (min_le_left _ _)
    have hfb := hb z hzr
    rw [abs_of_pos hz0] at hfb
    have hzj : z ^ j ≤ z ^ (k + 1) := pow_le_pow_of_le_one hz0.le hz1.le hk
    have hzk : 0 < z ^ k := by positivity
    have e1 : z ^ k / (k.factorial : ℝ) * d = f z - z ^ (k + 1) * F z := by rw [h1]; ring
    have e2 : |z ^ k / (k.factorial : ℝ) * d| = z ^ k * (|d| / k.factorial) := by
      rw [abs_mul, abs_div, abs_of_pos hzk, abs_of_pos hk!]; ring
    have e3 : |f z - z ^ (k + 1) * F z| ≤ |C| * z ^ (k + 1) + z ^ (k + 1) * (|F 0| + 1) := by
      refine (abs_sub _ _).trans ?_
      rw [abs_mul, abs_of_pos (by positivity : (0 : ℝ) < z ^ (k + 1))]
      have : |f z| ≤ |C| * z ^ (k + 1) := hfb.trans (by
        calc C * z ^ j ≤ |C| * z ^ j := by gcongr; exact le_abs_self C
          _ ≤ |C| * z ^ (k + 1) := by gcongr)
      nlinarith [abs_nonneg (F z), (by positivity : (0 : ℝ) < z ^ (k + 1))]
    rw [← e1, e2] at e3
    have : z ^ k * (|d| / k.factorial) ≤ z ^ k * ((|C| + |F 0| + 1) * z) := by
      calc z ^ k * (|d| / k.factorial) ≤ |C| * z ^ (k + 1) + z ^ (k + 1) * (|F 0| + 1) := e3
        _ = z ^ k * ((|C| + |F 0| + 1) * z) := by ring
    exact le_of_mul_le_mul_left this hzk
  have : |d| ≤ 0 := by
    have := (div_le_iff₀ hk!).mp key; simpa using this
  exact abs_nonpos_iff.mp this

/-- `PF` is stable under `x ↦ q x`, `q ∈ ℚ`. -/
lemma pf_scale {f : ℝ → ℝ} (hf : Q01.PF 0 f) (q : ℚ) : Q01.PF 0 (fun x => f ((q : ℝ) * x)) := by
  set u : ℝ →L[ℝ] ℝ := (q : ℝ) • ContinuousLinearMap.id ℝ ℝ with hu
  have hg : (fun x => f ((q : ℝ) * x)) = f ∘ u := by funext x; simp [hu]
  have hA : AnalyticAt ℝ (fun x => f ((q : ℝ) * x)) 0 := by
    rw [hg]; exact (by simpa [hu] using hf.1 : AnalyticAt ℝ f (u 0)).comp (u.analyticAt 0)
  refine ⟨hA, fun a => ?_⟩
  have h1 : HasFPowerSeriesAt (f ∘ u) ((FormalMultilinearSeries.ofScalars ℝ
      (fun n => iteratedDeriv n f 0 / (n.factorial : ℝ))).compContinuousLinearMap u) 0 :=
    HasFPowerSeriesAt.compContinuousLinearMap (by simpa [hu] using hf.1.hasFPowerSeriesAt)
  rw [← hg] at h1
  have h2 := hA.hasFPowerSeriesAt.eq_formalMultilinearSeries h1
  have h3 := congrArg (fun p : FormalMultilinearSeries ℝ ℝ ℝ => p a (fun _ => 1)) h2
  simp only [FormalMultilinearSeries.ofScalars_apply_eq, one_pow, smul_eq_mul, mul_one,
    FormalMultilinearSeries.compContinuousLinearMap_apply] at h3
  rw [h3]
  have h4 : (fun _ : Fin a => u 1) = fun _ => (q : ℝ) := by funext; simp [hu]
  rw [Function.comp_def, h4, FormalMultilinearSeries.ofScalars_apply_eq, smul_eq_mul]
  have := Q01.zp_smul (q ^ a) (hf.2 a)
  rw [add_zero] at this ⊢
  convert this using 1
  push_cast; ring

lemma Gf_Hf (x : ℝ) : Q01.Gf x * Q01.Hf x = 1 := by
  unfold Q01.Gf Q01.Hf; rw [← Real.exp_add]; simp

lemma pf_prod {ι : Type*} (s : Finset ι) (f : ι → ℝ → ℝ) (hf : ∀ i ∈ s, Q01.PF 0 (f i)) :
    Q01.PF 0 (fun x => ∏ i ∈ s, f i x) := by
  classical
  induction s using Finset.induction_on with
  | empty =>
    simp only [Finset.prod_empty]
    exact Q01.pf_cast (by simp) (Q01.pf_congr (Eventually.of_forall fun x => Gf_Hf x)
      (Q01.pf_mul Q01.pf_G Q01.pf_H))
  | insert i s hi ih =>
    simp only [Finset.prod_insert hi]
    exact Q01.pf_cast (by simp) (Q01.pf_mul (hf i (Finset.mem_insert_self _ _))
      (ih fun j hj => hf j (Finset.mem_insert_of_mem hj)))

/-- `log(Γ(1+u)Γ(1+v)/Γ(1+u+v)) = ℓ(u)+ℓ(v)-ℓ(u+v)`. -/
lemma gam2 {u v : ℝ} (hu : |u| ≤ 1 / 2) (hv : |v| ≤ 1 / 2) (huv : |u + v| ≤ 1 / 2) :
    Real.Gamma (1 + u) * Real.Gamma (1 + v) =
      Q01.Gf u * Q01.Gf v * Q01.Hf (u + v) * Real.Gamma (1 + (u + v)) := by
  have hA := Y05.logGamma2 (by norm_num : (1 / 2 : ℝ) < 1) hu hv huv
  have hs : ∀ x : ℝ, |x| ≤ 1 / 2 → HasSum (fun k : ℕ => (-1 : ℝ) ^ (k + 2) * X05.z (k + 2) /
      ((k + 2 : ℕ) : ℝ) * x ^ (k + 2)) (Q01.ell x) := by
    intro x hx
    have h0 := Q01.ell_hasSum (x := x) (by linarith)
    have h := (hasSum_nat_add_iff' 2).mpr h0
    simp [Finset.sum_range_succ, Q01.cz] at h
    convert h using 1
    funext k; push_cast; ring
  have hB := ((hs u hu).add (hs v hv)).sub (hs (u + v) huv)
  have heq := hA.unique (hB.congr_fun (fun k => by ring))
  have hG : ∀ x : ℝ, |x| ≤ 1 / 2 → 0 < Real.Gamma (1 + x) := fun x hx =>
    Real.Gamma_pos_of_pos (by cases abs_le.mp hx; linarith)
  have g1 := hG u hu; have g2 := hG v hv; have g3 := hG _ huv
  have hE := congrArg Real.exp heq
  rw [Real.exp_log (by positivity), Real.exp_sub, Real.exp_add] at hE
  unfold Q01.Gf Q01.Hf
  rw [Real.exp_neg]
  field_simp at hE ⊢
  linarith

/-- `∏ Γ(1+x_v) = ∏ G(x_v) · H(Σx) · Γ(1+Σx)` for small `x`. -/
lemma gam_prod : ∀ (n : ℕ) (x : Fin n → ℝ) (ρ : ℝ), (∀ v, |x v| ≤ ρ) → (n : ℝ) * ρ ≤ 1 / 2 →
    ∏ v, Real.Gamma (1 + x v) =
      (∏ v, Q01.Gf (x v)) * Q01.Hf (∑ v, x v) * Real.Gamma (1 + ∑ v, x v) := by
  intro n
  induction n with
  | zero =>
    intro x ρ _ _
    simp only [Finset.univ_eq_empty, Finset.prod_empty, Finset.sum_empty, one_mul, add_zero]
    rw [Real.Gamma_one, mul_one]
    unfold Q01.Hf
    simp [Q01.ell_zero]
  | succ n ih =>
    intro x ρ hx hρ
    have hρ0 : 0 ≤ ρ := (abs_nonneg _).trans (hx 0)
    have hS : ∀ m : ℕ, m ≤ n + 1 → ∀ y : Fin m → ℝ, (∀ v, |y v| ≤ ρ) → |∑ v, y v| ≤ m * ρ := by
      intro m _ y hy
      calc |∑ v, y v| ≤ ∑ v, |y v| := Finset.abs_sum_le_sum_abs _ _
        _ ≤ ∑ _v : Fin m, ρ := Finset.sum_le_sum (fun v _ => hy v)
        _ = m * ρ := by simp
    have hn : (n : ℝ) * ρ ≤ 1 / 2 := by
      have : (n : ℝ) * ρ ≤ ((n + 1 : ℕ) : ℝ) * ρ := by push_cast; nlinarith
      linarith
    rw [Fin.prod_univ_castSucc, Fin.prod_univ_castSucc, Fin.sum_univ_castSucc,
      ih (fun v => x (Fin.castSucc v)) ρ (fun v => hx _) hn]
    set S := ∑ v : Fin n, x (Fin.castSucc v)
    have h1 : |S| ≤ 1 / 2 := (hS n (by omega) _ (fun v => hx _)).trans hn
    have h2 : |x (Fin.last n)| ≤ 1 / 2 := (hx _).trans (by
      have : (1 : ℝ) * ρ ≤ ((n + 1 : ℕ) : ℝ) * ρ := by
        have : (1 : ℝ) ≤ ((n + 1 : ℕ) : ℝ) := by push_cast; linarith [(Nat.cast_nonneg n : (0:ℝ) ≤ n)]
        nlinarith
      linarith)
    have h3 : |S + x (Fin.last n)| ≤ 1 / 2 := by
      have := hS (n + 1) le_rfl x hx
      rw [Fin.sum_univ_castSucc] at this
      exact this.trans hρ
    rw [mul_assoc (_ * _), mul_assoc, gam2 h1 h2 h3]
    have := Gf_Hf S
    calc (∏ v : Fin n, Q01.Gf (x (Fin.castSucc v))) * (Q01.Hf S * (Q01.Gf S * Q01.Gf (x (Fin.last n)) *
          Q01.Hf (S + x (Fin.last n)) * Real.Gamma (1 + (S + x (Fin.last n)))))
        = (∏ v : Fin n, Q01.Gf (x (Fin.castSucc v))) * Q01.Gf (x (Fin.last n)) *
          Q01.Hf (S + x (Fin.last n)) * Real.Gamma (1 + (S + x (Fin.last n))) * (Q01.Gf S * Q01.Hf S) := by
          ring
      _ = _ := by rw [this, mul_one]

end R03

namespace R03

open Filter Topology

lemma orth_ind' {N : ℕ} (φ : Fin N → ℝ → ℝ) :
    (orth N).indicator (fun g : Fin N → ℝ => ∏ j, φ j (g j)) =
      fun g => ∏ j, (Ioi (0 : ℝ)).indicator (φ j) (g j) := by
  funext g
  by_cases hg : g ∈ orth N
  · rw [indicator_of_mem hg]
    exact Finset.prod_congr rfl (fun j _ => (indicator_of_mem (hg j (mem_univ j)) _).symm)
  · rw [indicator_of_notMem hg]
    have : ∃ j, g j ∉ Ioi (0 : ℝ) := by
      by_contra h; push Not at h; exact hg (fun j _ => h j)
    obtain ⟨j, hj⟩ := this
    exact (Finset.prod_eq_zero (Finset.mem_univ j) (indicator_of_notMem hj _)).symm

lemma prod_intOn' {N : ℕ} {φ : Fin N → ℝ → ℝ} (h : ∀ j, IntegrableOn (φ j) (Ioi 0)) :
    IntegrableOn (fun g : Fin N → ℝ => ∏ j, φ j (g j)) (orth N) := by
  rw [← integrable_indicator_iff (orth_meas N), orth_ind']
  exact Integrable.fintype_prod (f := fun j => (Ioi (0 : ℝ)).indicator (φ j))
    (fun j => (h j).integrable_indicator measurableSet_Ioi)

lemma prod_integral' {N : ℕ} (φ : Fin N → ℝ → ℝ) :
    ∫ g in orth N, ∏ j, φ j (g j) = ∏ j, ∫ x in Ioi 0, φ j x := by
  rw [← integral_indicator (orth_meas N), orth_ind', volume_pi,
    integral_fintype_prod_eq_prod (fun j => (Ioi (0 : ℝ)).indicator (φ j))]
  exact Finset.prod_congr rfl (fun j _ => integral_indicator measurableSet_Ioi)

/-- `ν_v = #{l : ι l = v}`. -/
def cnt {k N : ℕ} (ι : Fin k → Fin N) (v : Fin N) : ℕ := (Finset.univ.filter (fun l => ι l = v)).card

lemma cnt_sum {k N : ℕ} (ι : Fin k → Fin N) : ∑ v, (cnt ι v : ℝ) = k := by
  have := Finset.card_eq_sum_card_fiberwise (f := ι) (s := Finset.univ) (t := Finset.univ)
    (fun x _ => Finset.mem_univ _)
  simp only [Finset.card_univ, Fintype.card_fin] at this
  unfold cnt; exact_mod_cast this.symm

lemma cnt_le {k N : ℕ} (ι : Fin k → Fin N) (v : Fin N) : cnt ι v ≤ k := by
  unfold cnt
  calc _ ≤ (Finset.univ : Finset (Fin k)).card := Finset.card_filter_le _ _
    _ = k := by simp

/-- the rational multiplier `q_v = -(N-1)/N + 2(ν_v - k/N)`. -/
noncomputable def qc (N k : ℕ) (ι : Fin k → Fin N) (v : Fin N) : ℚ :=
  -((N : ℚ) - 1) / N + 2 * ((cnt ι v : ℚ) - k / N)

/-- `β = -(N-1)/N`, so `a = βτ`. -/
noncomputable def bq (N : ℕ) : ℚ := -((N : ℚ) - 1) / N

lemma aPar_eq (N : ℕ) (τ : ℝ) : aPar N τ = (bq N : ℝ) * τ := by
  unfold aPar bq; push_cast; ring

lemma qc_eq (N k : ℕ) (ι : Fin k → Fin N) (v : Fin N) (τ : ℝ) :
    (qc N k ι v : ℝ) * τ = aPar N τ + 2 * τ * ((cnt ι v : ℝ) - k / N) := by
  unfold qc aPar; push_cast; ring

lemma qc_sum {N k : ℕ} (hN : 1 ≤ N) (ι : Fin k → Fin N) (τ : ℝ) :
    ∑ v, (qc N k ι v : ℝ) * τ = N * aPar N τ := by
  have hN' : (N : ℝ) ≠ 0 := by have : (0:ℝ) < N := by exact_mod_cast hN
                               exact this.ne'
  simp_rw [qc_eq]
  rw [Finset.sum_add_distrib, ← Finset.mul_sum, Finset.sum_sub_distrib, cnt_sum]
  simp only [Finset.sum_const, Finset.card_univ, Fintype.card_fin, nsmul_eq_mul]
  field_simp
  ring

lemma qc_abs {N k : ℕ} (hN : 1 ≤ N) (ι : Fin k → Fin N) (v : Fin N) :
    |(qc N k ι v : ℝ)| ≤ 1 + 2 * k := by
  have hN' : (0 : ℝ) < N := by exact_mod_cast hN
  have hN1 : (1 : ℝ) ≤ N := by exact_mod_cast hN
  have hc := cnt_le ι v
  have hc' : (cnt ι v : ℝ) ≤ k := by exact_mod_cast hc
  have hc0 : (0 : ℝ) ≤ cnt ι v := Nat.cast_nonneg _
  have hk0 : (0 : ℝ) ≤ k := Nat.cast_nonneg _
  have h1 : 0 ≤ (k : ℝ) / N := by positivity
  have h2 : (k : ℝ) / N ≤ k := div_le_self hk0 hN1
  have h3 : 0 ≤ ((N : ℝ) - 1) / N := div_nonneg (by linarith) hN'.le
  have h4 : ((N : ℝ) - 1) / N ≤ 1 := by rw [div_le_one hN']; linarith
  unfold qc; push_cast
  rw [neg_div, abs_le]; constructor <;> nlinarith

/-- pointwise expansion of the product of exponentials. -/
lemma prod_e {N k : ℕ} (_hN : 1 ≤ N) (τ : ℝ) (ι : Fin k → Fin N) {g : Fin N → ℝ}
    (hg : g ∈ orth N) :
    ∏ l, Real.exp (2 * τ * (Real.log (g (ι l)) - (∑ j, Real.log (g j)) / N)) =
      ∏ v, g v ^ (2 * τ * ((cnt ι v : ℝ) - k / N)) := by
  have hgi : ∀ i, 0 < g i := fun i => hg i (mem_univ i)
  rw [← Real.exp_sum]
  rw [Finset.prod_congr rfl (fun v _ => Real.rpow_def_of_pos (hgi v) _), ← Real.exp_sum]
  congr 1
  have hfib : ∑ l, Real.log (g (ι l)) = ∑ v, (cnt ι v : ℝ) * Real.log (g v) := by
    rw [← Finset.sum_fiberwise Finset.univ ι (fun l => Real.log (g (ι l)))]
    refine Finset.sum_congr rfl (fun v _ => ?_)
    rw [Finset.sum_congr rfl (fun l hl => by rw [(Finset.mem_filter.mp hl).2]), Finset.sum_const,
      nsmul_eq_mul]
    rfl
  rw [← Finset.mul_sum, Finset.sum_sub_distrib, hfib]
  simp only [Finset.sum_const, Finset.card_univ, Fintype.card_fin, nsmul_eq_mul]
  have : ∑ v, Real.log (g v) * (2 * τ * ((cnt ι v : ℝ) - k / N)) =
      2 * τ * (∑ v, (cnt ι v : ℝ) * Real.log (g v)) - 2 * τ * (k / N) * ∑ j, Real.log (g j) := by
    rw [Finset.mul_sum, Finset.mul_sum, ← Finset.sum_sub_distrib]
    exact Finset.sum_congr rfl (fun v _ => by ring)
  rw [this]
  ring

/-- the Gamma integrand of coordinate `v`. -/
noncomputable def phiQ (b : ℝ) (x : ℝ) : ℝ := x ^ b * Real.exp (-x)

/-- coefficient `binom(j,k)(-1)^{j-k}N^{-k}`. -/
noncomputable def coefQ (N j k : ℕ) : ℚ := (j.choose k : ℚ) * (-1) ^ (j - k) * (1 / (N : ℚ)) ^ k

lemma dens_X_pow {N : ℕ} (hN : 1 ≤ N) (τ : ℝ) (j : ℕ) {g : Fin N → ℝ} (hg : g ∈ orth N) :
    dens N τ g * (Ymean N τ g - 1) ^ j =
      ∑ k ∈ Finset.range (j + 1), (coefQ N j k : ℝ) *
        ∑ ι : Fin k → Fin N, ∏ v, phiQ ((qc N k ι v : ℝ) * τ) (g v) := by
  have hgi : ∀ i, 0 < g i := fun i => hg i (mem_univ i)
  have hY : ∀ k : ℕ, Ymean N τ g ^ k = (1 / (N : ℝ)) ^ k *
      ∑ ι : Fin k → Fin N, ∏ v, g v ^ (2 * τ * ((cnt ι v : ℝ) - k / N)) := by
    intro k
    unfold Ymean
    rw [div_pow, Fintype.sum_pow, one_div_pow, div_eq_inv_mul, one_div]
    congr 1
    exact Finset.sum_congr rfl (fun ι _ => prod_e hN τ ι hg)
  have hp : ∀ k (ι : Fin k → Fin N), dens N τ g * ∏ v, g v ^ (2 * τ * ((cnt ι v : ℝ) - k / N)) =
      ∏ v, phiQ ((qc N k ι v : ℝ) * τ) (g v) := by
    intro k ι
    rw [dens_eq, ← Finset.prod_mul_distrib]
    refine Finset.prod_congr rfl (fun v _ => ?_)
    rw [phiQ, qc_eq, Real.rpow_add (hgi v)]
    ring
  rw [sub_eq_add_neg, add_pow, Finset.mul_sum]
  refine Finset.sum_congr rfl (fun k _ => ?_)
  have hs : ∑ ι : Fin k → Fin N, ∏ v, phiQ ((qc N k ι v : ℝ) * τ) (g v) =
      dens N τ g * ∑ ι : Fin k → Fin N, ∏ v, g v ^ (2 * τ * ((cnt ι v : ℝ) - k / N)) := by
    rw [Finset.mul_sum]; exact Finset.sum_congr rfl (fun ι _ => (hp k ι).symm)
  rw [hs, hY]
  unfold coefQ
  push_cast
  ring

lemma Ymean_cont_meas (N : ℕ) (τ : ℝ) (j : ℕ) :
    Measurable fun g => dens N τ g * (Ymean N τ g - 1) ^ j :=
  (dens_meas N τ).mul ((Ymean_meas N τ).sub_const 1 |>.pow_const j)

/-- the closed form of `E[X^j]` for small `τ`. -/
noncomputable def Fm (N j : ℕ) (τ : ℝ) : ℝ :=
  ∑ k ∈ Finset.range (j + 1), (coefQ N j k : ℝ) *
    ∑ ι : Fin k → Fin N, ∏ v, (Q01.Gf ((qc N k ι v : ℝ) * τ) * Q01.Hf ((bq N : ℝ) * τ))

lemma pf_Fm (N j : ℕ) : Q01.PF 0 (Fm N j) := by
  unfold Fm
  refine Q01.pf_sum _ _ (fun k _ => ?_)
  refine Q01.pf_smul (coefQ N j k) ?_
  refine Q01.pf_sum _ _ (fun ι _ => ?_)
  refine pf_prod _ _ (fun v _ => ?_)
  exact Q01.pf_cast (by simp) (Q01.pf_mul (pf_scale Q01.pf_G _) (pf_scale Q01.pf_H _))

lemma momX_eq {N : ℕ} (hN : 1 ≤ N) (j : ℕ) {τ : ℝ}
    (hτ : |τ| ≤ 1 / (2 * N * (1 + 2 * j))) : momX N j τ = Fm N j τ := by
  have hN' : (0 : ℝ) < N := by exact_mod_cast hN
  have hj0 : (0 : ℝ) ≤ j := Nat.cast_nonneg j
  set ρ := (1 + 2 * (j : ℝ)) * |τ| with hρ
  have hNρ : (N : ℝ) * ρ ≤ 1 / 2 := by
    rw [le_div_iff₀ (by positivity)] at hτ
    rw [hρ]; nlinarith
  have hρ0 : 0 ≤ ρ := by positivity
  have hρ2 : ρ ≤ 1 / 2 := by
    have : (1 : ℝ) ≤ N := by exact_mod_cast hN
    have := mul_le_mul_of_nonneg_right this hρ0
    linarith
  have hxv : ∀ k ∈ Finset.range (j + 1), ∀ (ι : Fin k → Fin N) (v : Fin N),
      |(qc N k ι v : ℝ) * τ| ≤ ρ := by
    intro k hk ι v
    have hk' : (k : ℝ) ≤ j := by exact_mod_cast Nat.lt_succ_iff.mp (Finset.mem_range.mp hk)
    rw [abs_mul, hρ]
    gcongr
    exact (qc_abs hN ι v).trans (by linarith)
  have ha : |aPar N τ| ≤ ρ := (aPar_abs hN τ).trans (by rw [hρ]; nlinarith [abs_nonneg τ])
  -- integrability of each coordinate factor
  have hφ : ∀ k ∈ Finset.range (j + 1), ∀ (ι : Fin k → Fin N) (v : Fin N),
      IntegrableOn (phiQ ((qc N k ι v : ℝ) * τ)) (Ioi 0) ∧
        ∫ x in Ioi 0, phiQ ((qc N k ι v : ℝ) * τ) x = Real.Gamma (1 + (qc N k ι v : ℝ) * τ) := by
    intro k hk ι v
    have := hxv k hk ι v
    exact gam_int (by cases abs_le.mp this; linarith)
  unfold momX gammaE
  rw [setIntegral_congr_fun (orth_meas N) (fun g hg => dens_X_pow hN τ j hg)]
  rw [integral_finsetSum _ (fun k hk => (integrable_finsetSum _ (fun ι _ =>
      prod_intOn' (fun v => (hφ k hk ι v).1))).const_mul _)]
  rw [Finset.mul_sum]
  unfold Fm
  refine Finset.sum_congr rfl (fun k hk => ?_)
  rw [integral_const_mul, integral_finsetSum _ (fun ι _ => prod_intOn' (fun v => (hφ k hk ι v).1)),
    Finset.mul_sum, Finset.mul_sum, Finset.mul_sum]
  refine Finset.sum_congr rfl (fun ι _ => ?_)
  rw [prod_integral', Finset.prod_congr rfl (fun v _ => (hφ k hk ι v).2)]
  -- the Gamma ratio
  have h1 := gam_prod N (fun v => (qc N k ι v : ℝ) * τ) ρ (hxv k hk ι) hNρ
  have h2 := gam_prod N (fun _ => aPar N τ) ρ (fun _ => ha) hNρ
  rw [qc_sum hN] at h1
  simp only [Finset.sum_const, Finset.card_univ, Fintype.card_fin, nsmul_eq_mul,
    Finset.prod_const] at h2
  rw [h1, h2, Finset.prod_mul_distrib, Finset.prod_const, Finset.card_univ, Fintype.card_fin,
    ← aPar_eq]
  have hNa : |(N : ℝ) * aPar N τ| ≤ 1 / 2 := by
    rw [abs_mul, abs_of_pos hN']
    calc (N : ℝ) * |aPar N τ| ≤ N * ρ := by gcongr
      _ ≤ 1 / 2 := hNρ
  have hG : 0 < Real.Gamma (1 + (N : ℝ) * aPar N τ) :=
    Real.Gamma_pos_of_pos (by linarith [(abs_le.mp hNa).1])
  have hH : 0 < Q01.Hf ((N : ℝ) * aPar N τ) := Real.exp_pos _
  have hGH : Q01.Gf (aPar N τ) ^ N * Q01.Hf (aPar N τ) ^ N = 1 := by
    rw [← mul_pow, Gf_Hf, one_pow]
  have hGa : 0 < Q01.Gf (aPar N τ) := Real.exp_pos _
  rw [eq_inv_of_mul_eq_one_right hGH]
  field_simp
  simp only [mul_comm τ]

end R03

theorem R03.proof_8 : R03.stmt_8 := by
  intro N hN j
  have hN1 : 1 ≤ N := by omega
  have hN' : (0 : ℝ) < N := by exact_mod_cast hN1
  set r : ℝ := 1 / (2 * N * (1 + 2 * j)) with hr
  have hr0 : 0 < r := by positivity
  have hev : R03.Fm N j =ᶠ[nhds 0] R03.momX N j := by
    have := Metric.ball_mem_nhds (0 : ℝ) hr0
    filter_upwards [this] with τ hτ
    rw [Metric.mem_ball, Real.dist_eq, sub_zero] at hτ
    exact (R03.momX_eq hN1 j hτ.le).symm
  have hPF := Q01.pf_congr hev (R03.pf_Fm N j)
  obtain ⟨C, τJ, hτJ, h7⟩ := R03.proof_7 N hN j
  refine ⟨⟨τJ, hτJ, fun τ hτ => (h7 τ hτ.le).1 j (by omega)⟩, hPF.1,
    fun k => by simpa using hPF.2 k, ?_⟩
  -- vanishing below order `j`
  rcases Nat.eq_zero_or_pos j with hj | hj
  · intro k hk; omega
  obtain ⟨C', τ', hτ', h7'⟩ := R03.proof_7 N hN (j - 1)
  rw [show j - 1 + 1 = j by omega] at h7'
  refine R03.coeff_vanish hPF.1 (C := C') (r := min τ' (1 / 2)) (lt_min hτ' (by norm_num))
    (fun τ hτ => ?_)
  have hτa : |τ| ≤ τ' := (lt_of_lt_of_le hτ (min_le_left _ _)).le
  have hτb : |τ| < 1 / 2 := lt_of_lt_of_le hτ (min_le_right _ _)
  obtain ⟨hint, hE, -⟩ := h7' τ hτa
  unfold R03.momX R03.gammaE
  unfold R03.gammaE at hE
  set a := R03.aPar N τ
  have hΓ : 0 < Real.Gamma (1 + a) := Real.Gamma_pos_of_pos (by
    have := abs_le.mp ((R03.aPar_abs hN1 τ).trans hτb.le)
    linarith [this.1])
  rw [abs_mul, abs_of_pos (inv_pos.mpr (pow_pos hΓ N))]
  refine le_trans ?_ hE
  gcongr
  refine (abs_integral_le_integral_abs).trans (le_of_eq ?_)
  refine setIntegral_congr_fun (R03.orth_meas N) (fun g hg => ?_)
  simp only
  rw [abs_mul, abs_of_nonneg (R03.dens_nn τ hg), abs_pow]


namespace R03

open Filter Topology

lemma fact_prod (c : ℝ) : ∀ m : ℕ,
    (m.factorial : ℝ) * ∏ i ∈ Finset.Icc 1 m, (1 - c / (i : ℝ)) = ∏ i ∈ Finset.Icc 1 m, ((i : ℝ) - c)
  | 0 => by simp
  | m + 1 => by
    rw [Finset.prod_Icc_succ_top (by omega), Finset.prod_Icc_succ_top (by omega), ← fact_prod c m,
      Nat.factorial_succ]
    push_cast
    have : ((m : ℝ) + 1) ≠ 0 := by positivity
    field_simp

lemma gam_shift {x : ℝ} (hx : 0 < x) : ∀ m : ℕ,
    Real.Gamma (x + m) = Real.Gamma x * ∏ i ∈ Finset.Icc 1 m, (x + i - 1)
  | 0 => by simp
  | m + 1 => by
    rw [Finset.prod_Icc_succ_top (by omega), ← mul_assoc, ← gam_shift hx m]
    have : x + ((m + 1 : ℕ) : ℝ) = (x + m) + 1 := by push_cast; ring
    rw [this, Real.Gamma_add_one (by positivity)]
    ring

/-- the flat-face factor `Φ = Γ(1+a)^N/Γ(1-(N-1)τ)` as `∏ G(βτ) · H(Nβτ)`. -/
noncomputable def Phif (N : ℕ) (τ : ℝ) : ℝ :=
  (∏ _v : Fin N, Q01.Gf ((bq N : ℝ) * τ)) * Q01.Hf (((N * bq N : ℚ) : ℝ) * τ)

lemma pf_Phif (N : ℕ) : Q01.PF 0 (Phif N) := by
  unfold Phif
  exact Q01.pf_cast (by simp) (Q01.pf_mul (pf_prod _ _ (fun v _ => pf_scale Q01.pf_G _))
    (pf_scale Q01.pf_H _))

lemma Phif_eq {N : ℕ} (hN : 1 ≤ N) {τ : ℝ} (hτ : (N : ℝ) * |τ| ≤ 1 / 2) :
    Phif N τ = Real.Gamma (1 + aPar N τ) ^ N / Real.Gamma (1 - ((N : ℝ) - 1) * τ) := by
  have hN' : (0 : ℝ) < N := by exact_mod_cast hN
  have ha := aPar_abs hN τ
  have h2 := gam_prod N (fun _ => aPar N τ) |τ| (fun _ => ha) hτ
  simp only [Finset.sum_const, Finset.card_univ, Fintype.card_fin, nsmul_eq_mul,
    Finset.prod_const] at h2
  have hNa : (N : ℝ) * aPar N τ = -(((N : ℝ) - 1) * τ) := by unfold aPar; field_simp
  rw [hNa] at h2
  have hG : 0 < Real.Gamma (1 + -(((N : ℝ) - 1) * τ)) := Real.Gamma_pos_of_pos (by
    have : ((N : ℝ) - 1) * τ ≤ 1 / 2 := by
      have h1 : ((N : ℝ) - 1) * τ ≤ ((N : ℝ) - 1) * |τ| :=
        mul_le_mul_of_nonneg_left (le_abs_self τ) (by
          have : (1 : ℝ) ≤ N := by exact_mod_cast hN
          linarith)
      nlinarith [abs_nonneg τ]
    linarith)
  rw [h2, show 1 - ((N : ℝ) - 1) * τ = 1 + -(((N : ℝ) - 1) * τ) by ring]
  unfold Phif
  rw [Finset.prod_const, Finset.card_univ, Fintype.card_fin, ← aPar_eq]
  have : (((N * bq N : ℚ) : ℝ) * τ) = -(((N : ℝ) - 1) * τ) := by
    push_cast; rw [mul_assoc, ← aPar_eq, hNa]
  rw [this]
  field_simp

/-- `binom(1/2,j)` is rational. -/
noncomputable def bhq (j : ℕ) : ℚ := (∏ i ∈ Finset.range j, ((1 : ℚ) / 2 - i)) / (j.factorial : ℚ)

lemma bhq_eq (j : ℕ) : (bhq j : ℝ) = binomHalf j := by
  unfold bhq binomHalf; push_cast; rfl

lemma pf_momX {N : ℕ} (hN : 2 ≤ N) (j : ℕ) : Q01.PF 0 (momX N j) :=
  ⟨(proof_8 N hN j).2.1, fun a => by simpa using (proof_8 N hN j).2.2.1 a⟩

/-- `B_K = Σ_{j≤K} binom(1/2,j) E[X^j]`. -/
noncomputable def Bpoly (N K : ℕ) (τ : ℝ) : ℝ :=
  ∑ j ∈ Finset.range (K + 1), binomHalf j * momX N j τ

/-- the analytic model `G_K = Φ · B_K`. -/
noncomputable def Gm (N K : ℕ) (τ : ℝ) : ℝ := Phif N τ * Bpoly N K τ

lemma pf_Gm {N : ℕ} (hN : 2 ≤ N) (K : ℕ) : Q01.PF 0 (Gm N K) := by
  unfold Gm Bpoly
  refine Q01.pf_cast (zero_add 0) (Q01.pf_mul (pf_Phif N) ?_)
  refine Q01.pf_sum _ _ (fun j _ => ?_)
  simp_rw [← bhq_eq]
  exact Q01.pf_smul _ (pf_momX hN j)

lemma momX_bd {N : ℕ} (hN : 2 ≤ N) {j : ℕ} (hj : 1 ≤ j) :
    ∃ C r : ℝ, 0 < r ∧ ∀ τ : ℝ, |τ| < r → |momX N j τ| ≤ C * |τ| ^ j := by
  have hN1 : 1 ≤ N := by omega
  obtain ⟨C', τ', hτ', h7'⟩ := proof_7 N hN (j - 1)
  rw [show j - 1 + 1 = j by omega] at h7'
  refine ⟨C', min τ' (1 / 2), lt_min hτ' (by norm_num), fun τ hτ => ?_⟩
  have hτa : |τ| ≤ τ' := (lt_of_lt_of_le hτ (min_le_left _ _)).le
  have hτb : |τ| < 1 / 2 := lt_of_lt_of_le hτ (min_le_right _ _)
  obtain ⟨-, hE, -⟩ := h7' τ hτa
  unfold momX gammaE
  unfold gammaE at hE
  set a := aPar N τ
  have hΓ : 0 < Real.Gamma (1 + a) := Real.Gamma_pos_of_pos (by
    have := abs_le.mp ((aPar_abs hN1 τ).trans hτb.le)
    linarith [this.1])
  rw [abs_mul, abs_of_pos (inv_pos.mpr (pow_pos hΓ N))]
  refine le_trans ?_ hE
  gcongr
  refine (abs_integral_le_integral_abs).trans (le_of_eq ?_)
  refine setIntegral_congr_fun (orth_meas N) (fun g hg => ?_)
  simp only
  rw [abs_mul, abs_of_nonneg (dens_nn τ hg), abs_pow]

/-- a sequence of analytic germs whose consecutive differences are `O(|x|^{K+1})` has stable
Taylor coefficients. -/
lemma stable_coeff {f : ℕ → ℝ → ℝ} (hf : ∀ K, AnalyticAt ℝ (f K) 0)
    (hd : ∀ K, ∃ C r : ℝ, 0 < r ∧ ∀ x : ℝ, |x| < r → |f (K + 1) x - f K x| ≤ C * |x| ^ (K + 1)) :
    ∀ K k, k ≤ K → iteratedDeriv k (f K) 0 = iteratedDeriv k (f k) 0 := by
  intro K
  induction K with
  | zero => intro k hk; rw [Nat.le_zero.mp hk]
  | succ K ih =>
    intro k hk
    rcases Nat.lt_or_ge k (K + 1) with hlt | hge
    · obtain ⟨C, r, hr, hb⟩ := hd K
      have hv := coeff_vanish ((hf (K + 1)).sub (hf K)) hr (fun x hx => hb x hx) k hlt
      rw [iteratedDeriv_sub (hf (K + 1)).contDiffAt (hf K).contDiffAt, sub_eq_zero] at hv
      rw [hv, ih k (by omega)]
    · rw [show k = K + 1 by omega]

lemma Gm_step {N : ℕ} (hN : 2 ≤ N) (K : ℕ) :
    ∃ C r : ℝ, 0 < r ∧ ∀ x : ℝ, |x| < r → |Gm N (K + 1) x - Gm N K x| ≤ C * |x| ^ (K + 1) := by
  obtain ⟨C, r, hr, hb⟩ := momX_bd hN (j := K + 1) (by omega)
  have hc : ContinuousAt (Phif N) 0 := (pf_Phif N).1.continuousAt
  obtain ⟨δ, hδ, hcont⟩ := Metric.continuousAt_iff.mp hc 1 one_pos
  refine ⟨(|Phif N 0| + 1) * |binomHalf (K + 1)| * |C|, min r δ, lt_min hr hδ, fun x hx => ?_⟩
  have hxr : |x| < r := lt_of_lt_of_le hx (min_le_left _ _)
  have hxδ : |x| < δ := lt_of_lt_of_le hx (min_le_right _ _)
  have hP : |Phif N x| ≤ |Phif N 0| + 1 := by
    have := hcont (x := x) (by rwa [Real.dist_eq, sub_zero])
    rw [Real.dist_eq] at this
    have := abs_sub_abs_le_abs_sub (Phif N x) (Phif N 0)
    linarith
  have hdiff : Gm N (K + 1) x - Gm N K x = Phif N x * (binomHalf (K + 1) * momX N (K + 1) x) := by
    unfold Gm Bpoly; rw [Finset.sum_range_succ _ (K + 1)]; ring
  rw [hdiff, abs_mul, abs_mul]
  have hm := hb x hxr
  have hm' : |momX N (K + 1) x| ≤ |C| * |x| ^ (K + 1) := hm.trans (by gcongr; exact le_abs_self C)
  calc |Phif N x| * (|binomHalf (K + 1)| * |momX N (K + 1) x|)
      ≤ (|Phif N 0| + 1) * (|binomHalf (K + 1)| * (|C| * |x| ^ (K + 1))) := by gcongr
    _ = _ := by ring

lemma Gm_zero {N : ℕ} (hN : 2 ≤ N) (K : ℕ) : Gm N K 0 = 1 := by
  have hN1 : 1 ≤ N := by omega
  have hP : Phif N 0 = 1 := by
    unfold Phif; simp [Q01.Gf, Q01.Hf, Q01.ell_zero]
  have hm : ∀ j ∈ Finset.range (K + 1), binomHalf j * momX N j 0 = if j = 0 then 1 else 0 := by
    intro j _
    rcases Nat.eq_zero_or_pos j with rfl | hj
    · simp only [↓reduceIte]
      rw [binomHalf]
      simp only [Finset.range_zero, Finset.prod_empty, Nat.factorial_zero, Nat.cast_one, div_one,
        one_mul]
      rw [momX_eq hN1 0 (by simp)]
      unfold Fm
      simp [coefQ, Q01.Gf, Q01.Hf, Q01.ell_zero]
    · have := (proof_8 N hN j).2.2.2 0 hj
      simp only [iteratedDeriv_zero] at this
      simp [this, show j ≠ 0 by omega]
  unfold Gm Bpoly
  rw [hP, one_mul, Finset.sum_congr rfl hm, Finset.sum_ite_eq' (Finset.range (K + 1)) 0]
  simp

end R03

namespace R03

open Filter Topology

/-- `eq:pw:factor`: `S_N(p)/(S_N(1) Y_N(τ)) = Φ(τ) B_N(τ)` for small `|τ|`. -/
lemma ratio_eq {N : ℕ} (hN : 2 ≤ N) {p : ℝ} (hp : 0 < p) (hτ : (N : ℝ) * |1 - 1 / p| ≤ 1 / 2) :
    S N p / (S N 1 * T04.Y N (1 - 1 / p)) = Phif N (1 - 1 / p) * B N (1 - 1 / p) := by
  have hN1 : 1 ≤ N := by omega
  set τ := 1 - 1 / p with hτdef
  have hN' : (1 : ℝ) ≤ N := by exact_mod_cast hN1
  set c := ((N : ℝ) - 1) * τ with hc
  have hc2 : c ≤ 1 / 2 := by
    have h1 : c ≤ ((N : ℝ) - 1) * |τ| := mul_le_mul_of_nonneg_left (le_abs_self τ) (by linarith)
    nlinarith [abs_nonneg τ]
  rw [proof_5 N hN p hp |>.2, proof_2 N hN, Phif_eq hN1 (by simpa [hτdef] using hτ)]
  rw [← hτdef]
  unfold T04.Y
  rw [show 1 - τ = 1 / p by rw [hτdef]; ring]
  have hcast : ((N - 1 : ℕ) : ℝ) = (N : ℝ) - 1 := by rw [Nat.cast_sub hN1]; simp
  have hx : 0 < 1 - c := by linarith
  have hshift := gam_shift hx (N - 1)
  rw [hcast, show 1 - c + ((N : ℝ) - 1) = (N : ℝ) - ((N : ℝ) - 1) * τ by rw [hc]; ring] at hshift
  have hfp := fact_prod c (N - 1)
  have hP : ∏ i ∈ Finset.Icc 1 (N - 1), (1 - c + (i : ℝ) - 1) = ∏ i ∈ Finset.Icc 1 (N - 1), ((i : ℝ) - c) :=
    Finset.prod_congr rfl (fun i _ => by ring)
  rw [hP] at hshift
  have hPpos : 0 < ∏ i ∈ Finset.Icc 1 (N - 1), ((i : ℝ) - c) := by
    refine Finset.prod_pos (fun i hi => ?_)
    have : (1 : ℝ) ≤ i := by exact_mod_cast (Finset.mem_Icc.mp hi).1
    linarith
  have hQ : ∏ i ∈ Finset.Icc 1 (N - 1), (1 - ((N : ℝ) - 1) * τ / (i : ℝ)) =
      ∏ i ∈ Finset.Icc 1 (N - 1), (1 - c / (i : ℝ)) := rfl
  rw [hQ, hshift]
  have hQne : ∏ i ∈ Finset.Icc 1 (N - 1), (1 - c / (i : ℝ)) ≠ 0 := by
    intro h; rw [h, mul_zero] at hfp; linarith
  have hG1 : 0 < Real.Gamma (1 - c) := Real.Gamma_pos_of_pos hx
  have hF : (0 : ℝ) < ((N - 1).factorial : ℝ) := by positivity
  rw [← hfp, ← hc]
  have hp' : (1 / p) ^ (N - 1) ≠ 0 := by positivity
  have hsN : 0 < Real.sqrt N := Real.sqrt_pos.mpr (by positivity)
  field_simp

lemma PB_bd {N : ℕ} (hN : 2 ≤ N) (K : ℕ) :
    ∃ C r : ℝ, 0 < r ∧ ∀ τ : ℝ, |τ| < r →
      |Phif N τ * B N τ - Gm N K τ| ≤ C * |τ| ^ (K + 1) := by
  obtain ⟨C7, τ7, hτ7, h7⟩ := proof_7 N hN K
  have hc : ContinuousAt (Phif N) 0 := (pf_Phif N).1.continuousAt
  obtain ⟨δ, hδ, hcont⟩ := Metric.continuousAt_iff.mp hc 1 one_pos
  refine ⟨(|Phif N 0| + 1) * |C7|, min τ7 δ, lt_min hτ7 hδ, fun τ hτ => ?_⟩
  have hτa : |τ| ≤ τ7 := (lt_of_lt_of_le hτ (min_le_left _ _)).le
  have hτb : |τ| < δ := lt_of_lt_of_le hτ (min_le_right _ _)
  have hP : |Phif N τ| ≤ |Phif N 0| + 1 := by
    have := hcont (x := τ) (by rwa [Real.dist_eq, sub_zero])
    rw [Real.dist_eq] at this
    have := abs_sub_abs_le_abs_sub (Phif N τ) (Phif N 0)
    linarith
  have hB := (h7 τ hτa).2.2
  have hB' : |B N τ - Bpoly N K τ| ≤ |C7| * |τ| ^ (K + 1) :=
    hB.trans (by gcongr; exact le_abs_self C7)
  have : Phif N τ * B N τ - Gm N K τ = Phif N τ * (B N τ - Bpoly N K τ) := by unfold Gm; ring
  rw [this, abs_mul]
  calc |Phif N τ| * |B N τ - Bpoly N K τ| ≤ (|Phif N 0| + 1) * (|C7| * |τ| ^ (K + 1)) := by gcongr
    _ = _ := by ring

/-- `log` of a pure germ with value `1` at `0` is pure (`lem:pw:pure`). -/
lemma pf_log {G : ℝ → ℝ} (hG : Q01.PF 0 G) (h0 : G 0 = 1) :
    Q01.PF 0 (fun x => Real.log (G x)) := by
  have hLan : AnalyticAt ℝ (fun x => Real.log (G x)) 0 := hG.1.log (by rw [h0]; norm_num)
  refine ⟨hLan, fun a => ?_⟩
  have hpos : ∀ᶠ x in 𝓝 (0 : ℝ), 0 < G x :=
    continuousAt_const.eventually_lt hG.1.continuousAt (by rw [h0]; norm_num)
  have hexp : (fun x => Real.exp (Real.log (G x))) =ᶠ[𝓝 0] G := by
    filter_upwards [hpos] with x hx
    exact Real.exp_log hx
  set c : ℕ → ℝ := fun n => iteratedDeriv n (fun x => Real.exp (Real.log (G x))) 0 /
    (n.factorial : ℝ)
  set w : ℕ → ℝ := fun n => iteratedDeriv n (fun x => Real.log (G x)) 0 / (n.factorial : ℝ)
  have hc0 : c 0 = 1 := by simp [c, h0]
  have hw0 : w 0 = 0 := by simp [w, h0]
  have hrec : ∀ n : ℕ, 1 ≤ n →
      (n : ℝ) * c n = ∑ k ∈ Finset.Icc 1 n, (k : ℝ) * w k * c (n - k) :=
    fun n hn => Q01.exp_rec hLan n hn
  have hcp : ∀ n, Q01.IsZetaPoly n (c n) := by
    intro n
    have := hG.2 n
    simp only [add_zero] at this
    simp only [c]
    rw [hexp.iteratedDeriv_eq]
    exact this
  have := ((Q01.proof_1).2.2.2 c w hc0 hw0 hrec).1 hcp a
  simpa using this

lemma log_step {N : ℕ} (hN : 2 ≤ N) (K : ℕ) :
    ∃ C r : ℝ, 0 < r ∧ ∀ x : ℝ, |x| < r →
      |Real.log (Gm N (K + 1) x) - Real.log (Gm N K x)| ≤ C * |x| ^ (K + 1) := by
  obtain ⟨C, r, hr, hb⟩ := Gm_step hN K
  have hc1 : ContinuousAt (Gm N (K + 1)) 0 := (pf_Gm hN (K + 1)).1.continuousAt
  have hc0 : ContinuousAt (Gm N K) 0 := (pf_Gm hN K).1.continuousAt
  obtain ⟨δ1, hδ1, h1⟩ := Metric.continuousAt_iff.mp hc1 (1 / 2) (by norm_num)
  obtain ⟨δ0, hδ0, h0⟩ := Metric.continuousAt_iff.mp hc0 (1 / 2) (by norm_num)
  refine ⟨2 * C, min r (min δ1 δ0), lt_min hr (lt_min hδ1 hδ0), fun x hx => ?_⟩
  have hxr : |x| < r := lt_of_lt_of_le hx (min_le_left _ _)
  have hx1 : |x| < δ1 := lt_of_lt_of_le hx ((min_le_right _ _).trans (min_le_left _ _))
  have hx0 : |x| < δ0 := lt_of_lt_of_le hx ((min_le_right _ _).trans (min_le_right _ _))
  have e1 := h1 (x := x) (by rwa [Real.dist_eq, sub_zero])
  have e0 := h0 (x := x) (by rwa [Real.dist_eq, sub_zero])
  rw [Real.dist_eq, Gm_zero hN] at e1 e0
  have l1 : 1 / 2 ≤ Gm N (K + 1) x := by rw [abs_lt] at e1; linarith
  have l0 : 1 / 2 ≤ Gm N K x := by rw [abs_lt] at e0; linarith
  calc _ ≤ 2 * |Gm N (K + 1) x - Gm N K x| := Q01.log_lip l1 l0
    _ ≤ 2 * (C * |x| ^ (K + 1)) := by gcongr; exact hb x hxr
    _ = _ := by ring

/-- the logarithmic coefficients `w_k(N)`: the `τ^k` Taylor coefficient of `log G_k`. -/
noncomputable def wN (N k : ℕ) : ℝ :=
  iteratedDeriv k (fun x => Real.log (Gm N k x)) 0 / (k.factorial : ℝ)

/-- `prop:pw:n` and `g7:lem:both` without the explicit `w₂, w₃`: one pure sequence for both signs. -/
lemma main_exp {N : ℕ} (hN : 2 ≤ N) :
    ∃ w : ℕ → ℝ, w = wN N ∧ (∀ k, Q01.IsZetaPoly k (w k)) ∧ ∀ K : ℕ, ∃ C τ₀ : ℝ, 0 < τ₀ ∧
      ∀ τ : ℝ, |τ| ≤ τ₀ →
        |Real.log (Phif N τ * B N τ) - ∑ k ∈ Finset.Icc 2 K, w k * τ ^ k| ≤ C * |τ| ^ (K + 1) := by
  set Lf : ℕ → ℝ → ℝ := fun K x => Real.log (Gm N K x) with hLf
  have hpf : ∀ K, Q01.PF 0 (Lf K) := fun K => pf_log (pf_Gm hN K) (Gm_zero hN K)
  have hst := stable_coeff (f := Lf) (fun K => (hpf K).1) (fun K => log_step hN K)
  set w : ℕ → ℝ := fun k => iteratedDeriv k (Lf k) 0 / (k.factorial : ℝ) with hw
  have hwp : ∀ k, Q01.IsZetaPoly k (w k) := fun k => by simpa using (hpf k).2 k
  have hw0 : w 0 = 0 := by simp [hw, hLf, Gm_zero hN]
  have hw1 : w 1 = 0 := Q01.zp_one_eq_zero (hwp 1)
  refine ⟨w, rfl, hwp, fun K => ?_⟩
  obtain ⟨C1, r1, hr1, hb1⟩ := PB_bd hN K
  obtain ⟨δ, hδ, D, hD⟩ := Q01.taylor_bd (hpf K).1 (K + 1)
  have hcG : ContinuousAt (Gm N K) 0 := (pf_Gm hN K).1.continuousAt
  obtain ⟨δ2, hδ2, hG⟩ := Metric.continuousAt_iff.mp hcG (1 / 4) (by norm_num)
  set τ₀ := min (min (r1 / 2) (δ / 2)) (min (δ2 / 2) (min 1 (1 / (4 * (|C1| + 1))))) with hτ₀
  refine ⟨2 * |C1| + |D|, τ₀, by positivity, fun τ hτ => ?_⟩
  have hτr : |τ| < r1 := by
    have := hτ.trans ((min_le_left _ _).trans (min_le_left _ _)); linarith
  have hτδ : |τ| < δ := by
    have := hτ.trans ((min_le_left _ _).trans (min_le_right _ _)); linarith
  have hτδ2 : |τ| < δ2 := by
    have := hτ.trans ((min_le_right _ _).trans (min_le_left _ _)); linarith
  have hτ1 : |τ| ≤ 1 :=
    hτ.trans ((min_le_right _ _).trans ((min_le_right _ _).trans (min_le_left _ _)))
  have hτc : |τ| ≤ 1 / (4 * (|C1| + 1)) :=
    hτ.trans ((min_le_right _ _).trans ((min_le_right _ _).trans (min_le_right _ _)))
  have hpow : |τ| ^ (K + 1) ≤ |τ| := by
    calc |τ| ^ (K + 1) ≤ |τ| ^ 1 := pow_le_pow_of_le_one (abs_nonneg τ) hτ1 (by omega)
      _ = |τ| := pow_one _
  have hb := hb1 τ hτr
  have hb' : |Phif N τ * B N τ - Gm N K τ| ≤ |C1| * |τ| ^ (K + 1) :=
    hb.trans (by gcongr; exact le_abs_self C1)
  have hsmall : |C1| * |τ| ^ (K + 1) ≤ 1 / 4 := by
    calc |C1| * |τ| ^ (K + 1) ≤ |C1| * |τ| := mul_le_mul_of_nonneg_left hpow (abs_nonneg _)
      _ ≤ |C1| * (1 / (4 * (|C1| + 1))) := mul_le_mul_of_nonneg_left hτc (abs_nonneg _)
      _ ≤ 1 / 4 := by
          rw [mul_one_div, div_le_div_iff₀ (by positivity) (by norm_num)]
          nlinarith [abs_nonneg C1]
  have hGτ := hG (x := τ) (by rwa [Real.dist_eq, sub_zero])
  rw [Real.dist_eq, Gm_zero hN] at hGτ
  have hGlo' : 3 / 4 < Gm N K τ := by rw [abs_lt] at hGτ; linarith
  have hGlo : 1 / 2 ≤ Gm N K τ := by linarith
  have hPlo : 1 / 2 ≤ Phif N τ * B N τ := by rw [abs_le] at hb'; linarith
  have hlip := Q01.log_lip hPlo hGlo
  have hT := hD τ hτδ
  have hsw : ∑ k ∈ Finset.Icc 2 K, w k * τ ^ k =
      ∑ k ∈ Finset.range (K + 1), iteratedDeriv k (Lf K) 0 / (k.factorial : ℝ) * τ ^ k := by
    have h1 : ∑ k ∈ Finset.range (K + 1), iteratedDeriv k (Lf K) 0 / (k.factorial : ℝ) * τ ^ k =
        ∑ k ∈ Finset.range (K + 1), w k * τ ^ k := by
      refine Finset.sum_congr rfl (fun k hk => ?_)
      rw [hw, hst K k (Nat.lt_succ_iff.mp (Finset.mem_range.mp hk))]
    rw [h1]
    refine Finset.sum_subset (fun k hk => by
      have := (Finset.mem_Icc.mp hk).2; simp; omega) (fun k hk hk' => ?_)
    have hk2 : k < 2 := by
      by_contra h; push Not at h
      exact hk' (Finset.mem_Icc.mpr ⟨h, Nat.lt_succ_iff.mp (Finset.mem_range.mp hk)⟩)
    interval_cases k
    · rw [hw0, zero_mul]
    · rw [hw1, zero_mul]
  rw [hsw]
  have hT' : D * |τ| ^ (K + 1) ≤ |D| * |τ| ^ (K + 1) :=
    mul_le_mul_of_nonneg_right (le_abs_self _) (by positivity)
  calc |Real.log (Phif N τ * B N τ) -
        ∑ k ∈ Finset.range (K + 1), iteratedDeriv k (Lf K) 0 / (k.factorial : ℝ) * τ ^ k|
      ≤ |Real.log (Phif N τ * B N τ) - Real.log (Gm N K τ)| +
        |Lf K τ - ∑ k ∈ Finset.range (K + 1), iteratedDeriv k (Lf K) 0 / (k.factorial : ℝ) * τ ^ k| := by
        have := abs_add_le (Real.log (Phif N τ * B N τ) - Real.log (Gm N K τ))
          (Lf K τ - ∑ k ∈ Finset.range (K + 1), iteratedDeriv k (Lf K) 0 / (k.factorial : ℝ) * τ ^ k)
        simpa [hLf] using this
    _ ≤ 2 * (|C1| * |τ| ^ (K + 1)) + |D| * |τ| ^ (K + 1) := by
        have := mul_le_mul_of_nonneg_left hb' (by norm_num : (0:ℝ) ≤ 2)
        linarith
    _ = _ := by ring

/-- the expansion in terms of `p`. -/
lemma main_p {N : ℕ} (hN : 2 ≤ N) :
    ∃ w : ℕ → ℝ, w = wN N ∧ (∀ k, Q01.IsZetaPoly k (w k)) ∧ ∀ K : ℕ, ∃ C τ₀ : ℝ, 0 < τ₀ ∧
      ∀ p : ℝ, 0 < p → |1 - 1 / p| ≤ τ₀ →
        |Real.log (S N p / (S N 1 * T04.Y N (1 - 1 / p))) -
          ∑ k ∈ Finset.Icc 2 K, w k * (1 - 1 / p) ^ k| ≤ C * |1 - 1 / p| ^ (K + 1) := by
  obtain ⟨w, hw0, hw, hK⟩ := main_exp hN
  refine ⟨w, hw0, hw, fun K => ?_⟩
  obtain ⟨C, τ₀, hτ₀, hb⟩ := hK K
  have hN' : (0 : ℝ) < N := by have : (2 : ℝ) ≤ N := by exact_mod_cast hN
                               linarith
  refine ⟨C, min τ₀ (1 / (2 * N)), lt_min hτ₀ (by positivity), fun p hp hτ => ?_⟩
  have h1 : (N : ℝ) * |1 - 1 / p| ≤ 1 / 2 := by
    have := hτ.trans (min_le_right _ _)
    rw [le_div_iff₀ (by positivity)] at this
    linarith
  rw [ratio_eq hN hp h1]
  exact hb _ (hτ.trans (min_le_left _ _))

end R03

theorem R03.proof_9 : R03.stmt_9 := by
  intro N hN
  obtain ⟨w, -, hw, hK⟩ := R03.main_p hN
  refine ⟨w, hw, fun K => ?_⟩
  obtain ⟨C, τ₀, hτ₀, hb⟩ := hK K
  refine ⟨C, τ₀, hτ₀, fun p hp hτ => ?_⟩
  have hτ0 : 0 ≤ 1 - 1 / p := by
    have : 1 / p < 1 := by rw [div_lt_one (by linarith)]; exact hp
    linarith
  have := hb p (by linarith) (by rwa [abs_of_nonneg hτ0])
  rwa [abs_of_nonneg hτ0] at this


namespace R03

open Filter Topology

lemma exp_sub_one_sub_abs (u : ℝ) : |Real.exp u - 1 - u| ≤ u ^ 2 * Real.exp |u| := by
  have h := Complex.norm_exp_sub_sum_le_norm_mul_exp (u : ℂ) 2
  simp only [Finset.sum_range_succ, Finset.sum_range_zero, pow_zero, Nat.factorial_zero,
    Nat.cast_one, div_one, zero_add, pow_one, Nat.factorial_one] at h
  rw [show Complex.exp (u : ℂ) - (1 + (u : ℂ)) = ((Real.exp u - 1 - u : ℝ) : ℂ) by
    push_cast; ring, Complex.norm_real, Complex.norm_real, Real.norm_eq_abs, Real.norm_eq_abs,
    sq_abs] at h
  exact h

/-- `X = O(τ²)`: `|X| ≤ 16τ²Ξ² e^{4|τ|Ξ}` (the sum of the `m_i` vanishes). -/
lemma X_bd2 {N : ℕ} (hN : 1 ≤ N) (τ : ℝ) (g : Fin N → ℝ) :
    |Ymean N τ g - 1| ≤ 16 * τ ^ 2 * Xi g ^ 2 * Real.exp (4 * |τ| * Xi g) := by
  have hN' : (0 : ℝ) < N := by exact_mod_cast hN
  have hN1 : (1 : ℝ) ≤ N := by exact_mod_cast hN
  set L := ∑ j, Real.log (g j)
  set u : Fin N → ℝ := fun i => 2 * τ * (Real.log (g i) - L / N) with hu
  have hL : |L| ≤ Xi g := Finset.abs_sum_le_sum_abs _ _
  have hub : ∀ i, |u i| ≤ 4 * |τ| * Xi g := by
    intro i
    have hi : |Real.log (g i)| ≤ Xi g :=
      Finset.single_le_sum (f := fun j => |Real.log (g j)|) (fun j _ => abs_nonneg _)
        (Finset.mem_univ i)
    have hLN : |L / N| ≤ Xi g := by
      rw [abs_div, abs_of_pos hN', div_le_iff₀ hN']; nlinarith [abs_nonneg L, Xi_nn g]
    simp only [hu]
    rw [abs_mul, abs_mul, abs_two]
    have := abs_sub (Real.log (g i)) (L / N)
    have : 0 ≤ |τ| := abs_nonneg τ
    nlinarith
  have hsum0 : ∑ i, u i = 0 := by
    simp only [hu, ← Finset.mul_sum, Finset.sum_sub_distrib, Finset.sum_const, Finset.card_univ,
      Fintype.card_fin, nsmul_eq_mul]
    field_simp; ring
  have hY : Ymean N τ g - 1 = (∑ i, (Real.exp (u i) - 1 - u i)) / N := by
    rw [Finset.sum_sub_distrib, hsum0, sub_zero, Finset.sum_sub_distrib]
    unfold Ymean
    simp only [Finset.sum_const, Finset.card_univ, Fintype.card_fin, nsmul_eq_mul, mul_one]
    rw [sub_div, div_self hN'.ne']
  have hmono : ∀ i, u i ^ 2 * Real.exp |u i| ≤ 16 * τ ^ 2 * Xi g ^ 2 * Real.exp (4 * |τ| * Xi g) := by
    intro i
    have h1 : u i ^ 2 ≤ (4 * |τ| * Xi g) ^ 2 := by
      rw [← sq_abs]; exact pow_le_pow_left₀ (abs_nonneg _) (hub i) 2
    have h2 : Real.exp |u i| ≤ Real.exp (4 * |τ| * Xi g) := Real.exp_le_exp.mpr (hub i)
    calc u i ^ 2 * Real.exp |u i| ≤ (4 * |τ| * Xi g) ^ 2 * Real.exp (4 * |τ| * Xi g) :=
          mul_le_mul h1 h2 (Real.exp_pos _).le (sq_nonneg _)
      _ = _ := by rw [mul_pow, mul_pow, sq_abs]; ring
  rw [hY, abs_div, abs_of_pos hN', div_le_iff₀ hN']
  calc |∑ i, (Real.exp (u i) - 1 - u i)| ≤ ∑ i, |Real.exp (u i) - 1 - u i| :=
        Finset.abs_sum_le_sum_abs _ _
    _ ≤ ∑ _i : Fin N, 16 * τ ^ 2 * Xi g ^ 2 * Real.exp (4 * |τ| * Xi g) :=
        Finset.sum_le_sum (fun i _ => (exp_sub_one_sub_abs _).trans (hmono i))
    _ = _ := by simp; ring

/-- `E[X^j] = O(τ^{2j})`. -/
lemma momX_bd2 {N : ℕ} (hN : 2 ≤ N) (j : ℕ) :
    ∃ C r : ℝ, 0 < r ∧ ∀ τ : ℝ, |τ| < r → |momX N j τ| ≤ C * (τ ^ 2) ^ j := by
  have hN1 : 1 ≤ N := by omega
  obtain ⟨m, M, hm, hmM⟩ := gamma_bds
  set K : ℝ := 4 ^ (2 * j) * ((2 * j).factorial : ℝ) with hK
  refine ⟨(m ^ N)⁻¹ * (16 ^ j * K * (2 * M) ^ N), 1 / (32 * ((j : ℝ) + 1)), by positivity,
    fun τ hτ => ?_⟩
  have hj0 : (0 : ℝ) ≤ j := Nat.cast_nonneg j
  have hτ16 : |τ| ≤ 1 / 16 := by
    refine hτ.le.trans ?_
    rw [div_le_div_iff₀ (by positivity) (by norm_num)]; linarith
  have hτj : 4 * (j : ℝ) * |τ| ≤ 1 / 4 := by
    have := hτ.le
    rw [le_div_iff₀ (by positivity)] at this
    nlinarith [abs_nonneg τ]
  set a := aPar N τ with ha
  have hab' := abs_le.mp ((aPar_abs hN1 τ).trans hτ16)
  have hΓa := hmM (1 + a) ⟨by linarith, by linarith⟩
  have hΓ3 := hmM (a + 3 / 2) ⟨by linarith, by linarith⟩
  have hΓ1 := hmM (a + 1 / 2) ⟨by linarith, by linarith⟩
  have hΓpos : 0 < Real.Gamma (1 + a) := lt_of_lt_of_le hm hΓa.1
  have hMpos : 0 < M := lt_of_lt_of_le hm (hΓa.1.trans hΓa.2)
  obtain ⟨hφi, hφv⟩ := phiH_int (a := a) (by linarith)
  set D : (Fin N → ℝ) → ℝ := fun g => (16 * τ ^ 2) ^ j * K * ∏ v, phiH a (g v) with hD
  have hDi : IntegrableOn D (orth N) := (prod_intOn hφi).const_mul _
  have hDv : ∫ g in orth N, D g = (16 * τ ^ 2) ^ j * K *
      (Real.Gamma (a + 3 / 2) + Real.Gamma (a + 1 / 2)) ^ N := by
    rw [hD, integral_const_mul, prod_integral, hφv]
  have hpt : ∀ g ∈ orth N, ‖dens N τ g * (Ymean N τ g - 1) ^ j‖ ≤ D g := by
    intro g hg
    have hgi : ∀ i, 0 < g i := fun i => hg i (mem_univ i)
    have hd0 := dens_nn τ hg
    have hX := Xi_nn g
    rw [Real.norm_eq_abs, abs_mul, abs_of_nonneg hd0, abs_pow]
    have h1 := pow_le_pow_left₀ (abs_nonneg _) (X_bd2 hN1 τ g) j
    have h2 : (Xi g / 4) ^ (2 * j) / ((2 * j).factorial : ℝ) ≤ Real.exp (Xi g / 4) :=
      Real.pow_div_factorial_le_exp (Xi g / 4) (by positivity) (2 * j)
    have h2' : (Xi g ^ 2) ^ j ≤ K * Real.exp (Xi g / 4) := by
      rw [div_pow, div_div, div_le_iff₀ (by positivity)] at h2
      rw [← pow_mul, hK]; nlinarith
    have h3 : Real.exp (4 * |τ| * Xi g) ^ j ≤ Real.exp (Xi g / 4) := by
      rw [← Real.exp_nat_mul, Real.exp_le_exp]; nlinarith
    have h4 : Real.exp (Xi g / 2) ≤ ∏ v, (g v ^ (1 / 2 : ℝ) + g v ^ (-(1 / 2) : ℝ)) := by
      rw [Xi, Finset.sum_div, Real.exp_sum]
      exact Finset.prod_le_prod₀ (fun j _ => (Real.exp_pos _).le)
        (fun j _ => exp_half_log (hgi j))
    have h5 : dens N τ g * ∏ v, (g v ^ (1 / 2 : ℝ) + g v ^ (-(1 / 2) : ℝ)) =
        ∏ v, phiH a (g v) := by
      rw [dens_eq, ← Finset.prod_mul_distrib]; rfl
    have he : Real.exp (Xi g / 2) = Real.exp (Xi g / 4) * Real.exp (Xi g / 4) := by
      rw [← Real.exp_add]; ring_nf
    calc dens N τ g * |Ymean N τ g - 1| ^ j
        ≤ dens N τ g * ((16 * τ ^ 2) ^ j * (Xi g ^ 2) ^ j * Real.exp (4 * |τ| * Xi g) ^ j) := by
          gcongr; rw [← mul_pow, ← mul_pow]; exact h1
      _ ≤ dens N τ g * ((16 * τ ^ 2) ^ j * (K * Real.exp (Xi g / 4)) * Real.exp (Xi g / 4)) := by
          gcongr
      _ = dens N τ g * ((16 * τ ^ 2) ^ j * K * Real.exp (Xi g / 2)) := by rw [he]; ring
      _ ≤ dens N τ g * ((16 * τ ^ 2) ^ j * K * ∏ v, (g v ^ (1 / 2 : ℝ) + g v ^ (-(1 / 2) : ℝ))) := by
          gcongr
      _ = D g := by simp only [hD]; rw [← h5]; ring
  unfold momX gammaE
  rw [← ha, abs_mul, abs_of_pos (inv_pos.mpr (pow_pos hΓpos N))]
  have hI := norm_integral_le_of_norm_le hDi (ae_restrict_of_forall_mem (orth_meas N) hpt)
  rw [Real.norm_eq_abs, hDv] at hI
  have hinv : (Real.Gamma (1 + a) ^ N)⁻¹ ≤ (m ^ N)⁻¹ :=
    inv_anti₀ (by positivity) (pow_le_pow_left₀ hm.le hΓa.1 N)
  have hI2 : (16 * τ ^ 2) ^ j * K * (Real.Gamma (a + 3 / 2) + Real.Gamma (a + 1 / 2)) ^ N ≤
      (16 * τ ^ 2) ^ j * K * (2 * M) ^ N := by
    gcongr
    · linarith [hΓ3.1, hΓ1.1]
    · linarith [hΓ3.2, hΓ1.2]
  calc (Real.Gamma (1 + a) ^ N)⁻¹ * |∫ g in orth N, dens N τ g * (Ymean N τ g - 1) ^ j|
      ≤ (m ^ N)⁻¹ * ((16 * τ ^ 2) ^ j * K * (2 * M) ^ N) :=
        mul_le_mul hinv (hI.trans hI2) (abs_nonneg _) (by positivity)
    _ = _ := by rw [mul_pow]; ring

end R03

namespace R03

open Filter Topology

/-- the `k`-th Taylor coefficient at `0`. -/
noncomputable def cf (f : ℝ → ℝ) (k : ℕ) : ℝ := iteratedDeriv k f 0 / (k.factorial : ℝ)

lemma cf_add {f g : ℝ → ℝ} (hf : AnalyticAt ℝ f 0) (hg : AnalyticAt ℝ g 0) (k : ℕ) :
    cf (fun x => f x + g x) k = cf f k + cf g k := by
  unfold cf
  rw [show (fun x => f x + g x) = f + g from rfl,
    iteratedDeriv_add hf.contDiffAt hg.contDiffAt, add_div]

lemma cf_sub {f g : ℝ → ℝ} (hf : AnalyticAt ℝ f 0) (hg : AnalyticAt ℝ g 0) (k : ℕ) :
    cf (fun x => f x - g x) k = cf f k - cf g k := by
  unfold cf
  rw [show (fun x => f x - g x) = f - g from rfl,
    iteratedDeriv_sub hf.contDiffAt hg.contDiffAt, sub_div]

lemma cf_smul (c : ℝ) (f : ℝ → ℝ) (k : ℕ) : cf (fun x => c * f x) k = c * cf f k := by
  unfold cf; rw [iteratedDeriv_const_mul_field, mul_div_assoc]

lemma cf_scale {f : ℝ → ℝ} (hf : AnalyticAt ℝ f 0) (q : ℝ) (a : ℕ) :
    cf (fun x => f (q * x)) a = q ^ a * cf f a := by
  unfold cf
  set u : ℝ →L[ℝ] ℝ := q • ContinuousLinearMap.id ℝ ℝ with hu
  have hg : (fun x => f (q * x)) = f ∘ u := by funext x; simp [hu]
  have hA : AnalyticAt ℝ (fun x => f (q * x)) 0 := by
    rw [hg]; exact (by simpa [hu] using hf : AnalyticAt ℝ f (u 0)).comp (u.analyticAt 0)
  have h1 : HasFPowerSeriesAt (f ∘ u) ((FormalMultilinearSeries.ofScalars ℝ
      (fun n => iteratedDeriv n f 0 / (n.factorial : ℝ))).compContinuousLinearMap u) 0 :=
    HasFPowerSeriesAt.compContinuousLinearMap (by simpa [hu] using hf.hasFPowerSeriesAt)
  rw [← hg] at h1
  have h2 := hA.hasFPowerSeriesAt.eq_formalMultilinearSeries h1
  have h3 := congrArg (fun p : FormalMultilinearSeries ℝ ℝ ℝ => p a (fun _ => 1)) h2
  simp only [FormalMultilinearSeries.ofScalars_apply_eq, one_pow, smul_eq_mul, mul_one,
    FormalMultilinearSeries.compContinuousLinearMap_apply] at h3
  rw [h3]
  have h4 : (fun _ : Fin a => u 1) = fun _ => q := by funext; simp [hu]
  rw [Function.comp_def, h4, FormalMultilinearSeries.ofScalars_apply_eq, smul_eq_mul]
  ring

lemma ellq_an (q : ℝ) : AnalyticAt ℝ (fun x => Q01.ell (q * x)) 0 :=
  (by simpa using Q01.ell_analytic (x := q * 0) (by simp) : AnalyticAt ℝ Q01.ell (q * 0)).comp
    (analyticAt_const.mul analyticAt_id)

lemma ellq_cf (q : ℝ) (k : ℕ) : cf (fun x => Q01.ell (q * x)) k = q ^ k * Q01.cz k := by
  rw [cf_scale (Q01.ell_analytic (by simp)) q k]
  unfold cf; rw [Q01.ell_coef]

/-- `log Φ = N ℓ(βτ) - ℓ(Nβτ)`. -/
noncomputable def Afun (N : ℕ) (x : ℝ) : ℝ :=
  (N : ℝ) * Q01.ell ((bq N : ℝ) * x) - Q01.ell (((N : ℝ) * bq N) * x)

/-- `log E[e^{2τ m_1}] = ℓ(q₊τ) + (N-1)ℓ(q₋τ) - Nℓ(βτ)`. -/
noncomputable def Efun (N : ℕ) (x : ℝ) : ℝ :=
  Q01.ell ((((N : ℝ) - 1) / N) * x) + (((N : ℝ) - 1) * Q01.ell ((-((N : ℝ) + 1) / N) * x) -
    (N : ℝ) * Q01.ell ((bq N : ℝ) * x))

lemma Afun_an (N : ℕ) : AnalyticAt ℝ (Afun N) 0 :=
  (analyticAt_const.mul (ellq_an _)).sub (ellq_an _)

lemma Efun_an (N : ℕ) : AnalyticAt ℝ (Efun N) 0 :=
  (ellq_an _).add ((analyticAt_const.mul (ellq_an _)).sub (analyticAt_const.mul (ellq_an _)))

lemma Afun_cf (N k : ℕ) : cf (Afun N) k =
    ((N : ℝ) * (bq N : ℝ) ^ k - ((N : ℝ) * bq N) ^ k) * Q01.cz k := by
  unfold Afun
  rw [cf_sub (f := fun x => (N : ℝ) * Q01.ell ((bq N : ℝ) * x))
      (g := fun x => Q01.ell (((N : ℝ) * bq N) * x)) (analyticAt_const.mul (ellq_an _)) (ellq_an _),
    cf_smul (N : ℝ) (fun x => Q01.ell ((bq N : ℝ) * x)), ellq_cf, ellq_cf]; ring

lemma Efun_cf (N k : ℕ) : cf (Efun N) k =
    ((((N : ℝ) - 1) / N) ^ k + ((N : ℝ) - 1) * (-((N : ℝ) + 1) / N) ^ k -
      (N : ℝ) * (bq N : ℝ) ^ k) * Q01.cz k := by
  unfold Efun
  rw [cf_add (f := fun x => Q01.ell ((((N : ℝ) - 1) / N) * x))
      (g := fun x => ((N : ℝ) - 1) * Q01.ell ((-((N : ℝ) + 1) / N) * x) -
        (N : ℝ) * Q01.ell ((bq N : ℝ) * x)) (ellq_an _)
      ((analyticAt_const.mul (ellq_an _)).sub (analyticAt_const.mul (ellq_an _))),
    cf_sub (f := fun x => ((N : ℝ) - 1) * Q01.ell ((-((N : ℝ) + 1) / N) * x))
      (g := fun x => (N : ℝ) * Q01.ell ((bq N : ℝ) * x))
      (analyticAt_const.mul (ellq_an _)) (analyticAt_const.mul (ellq_an _)),
    cf_smul ((N : ℝ) - 1) (fun x => Q01.ell ((-((N : ℝ) + 1) / N) * x)),
    cf_smul (N : ℝ) (fun x => Q01.ell ((bq N : ℝ) * x)),
    ellq_cf, ellq_cf, ellq_cf]; ring

lemma Phif_exp (N : ℕ) (x : ℝ) : Phif N x = Real.exp (Afun N x) := by
  unfold Phif Afun Q01.Gf Q01.Hf
  rw [Finset.prod_const, Finset.card_univ, Fintype.card_fin, ← Real.exp_nat_mul, ← Real.exp_add]
  push_cast; ring_nf

lemma Fm0 (N : ℕ) (x : ℝ) : Fm N 0 x = 1 := by
  unfold Fm
  simp only [zero_add, Finset.sum_range_one, coefQ, Nat.choose_self, Nat.cast_one, Nat.sub_self,
    pow_zero, mul_one, Rat.cast_one, one_mul]
  have : ∀ (ι : Fin 0 → Fin N) (v : Fin N), qc N 0 ι v = bq N := by
    intro ι v; unfold qc bq cnt; simp
  simp only [this, Gf_Hf, Finset.prod_const_one, Finset.sum_const, Finset.card_univ]
  simp

lemma cnt1 {N : ℕ} (ι : Fin 1 → Fin N) (v : Fin N) : cnt ι v = if v = ι 0 then 1 else 0 := by
  unfold cnt
  rw [Finset.univ_unique, Finset.filter_singleton]
  by_cases h : v = ι 0
  · simp [h, Fin.default_eq_zero]
  · simp [h]; intro h'; exact h h'.symm

lemma Fm1 {N : ℕ} (hN : 1 ≤ N) (x : ℝ) : Fm N 1 x = Real.exp (Efun N x) - 1 := by
  have hN' : (0 : ℝ) < N := by exact_mod_cast hN
  have hprod : ∀ ι : Fin 1 → Fin N,
      ∏ v, (Q01.Gf ((qc N 1 ι v : ℝ) * x) * Q01.Hf ((bq N : ℝ) * x)) = Real.exp (Efun N x) := by
    intro ι
    unfold Q01.Gf Q01.Hf
    simp_rw [← Real.exp_add]
    rw [← Real.exp_sum]
    congr 1
    rw [Finset.sum_add_distrib, Finset.sum_const, Finset.card_univ, Fintype.card_fin,
      ← Finset.add_sum_erase _ _ (Finset.mem_univ (ι 0))]
    have h1 : (qc N 1 ι (ι 0) : ℝ) = ((N : ℝ) - 1) / N := by
      unfold qc; rw [cnt1]; simp only [↓reduceIte]; push_cast; field_simp; ring
    have h2 : ∀ v ∈ Finset.univ.erase (ι 0), (qc N 1 ι v : ℝ) = -((N : ℝ) + 1) / N := by
      intro v hv
      have : v ≠ ι 0 := Finset.ne_of_mem_erase hv
      unfold qc; rw [cnt1]; simp only [this, ↓reduceIte]; push_cast; field_simp; ring
    rw [h1, Finset.sum_congr rfl (fun v hv => by rw [h2 v hv]), Finset.sum_const,
      Finset.card_erase_of_mem (Finset.mem_univ _), Finset.card_univ, Fintype.card_fin,
      nsmul_eq_mul, Nat.cast_sub hN]
    unfold Efun
    push_cast
    ring
  unfold Fm
  rw [Finset.sum_range_succ, Finset.sum_range_one]
  simp only [hprod, Finset.sum_const, Finset.card_univ]
  have h0 : ∑ ι : Fin 0 → Fin N, ∏ v, (Q01.Gf ((qc N 0 ι v : ℝ) * x) * Q01.Hf ((bq N : ℝ) * x)) = 1 := by
    have : ∀ (ι : Fin 0 → Fin N) (v : Fin N), qc N 0 ι v = bq N := by
      intro ι v; unfold qc bq cnt; simp
    simp [this, Gf_Hf]
  rw [h0]
  simp [coefQ]
  field_simp
  ring

lemma logcosh (y : ℝ) : |Real.log ((1 + Real.exp y) / 2) - y / 2| ≤ y ^ 2 / 8 := by
  have hc : (1 + Real.exp y) / 2 = Real.exp (y / 2) * Real.cosh (y / 2) := by
    rw [Real.cosh_eq]
    have : Real.exp (y / 2) * Real.exp (y / 2) = Real.exp y := by rw [← Real.exp_add]; ring_nf
    have h2 : Real.exp (y / 2) * Real.exp (-(y / 2)) = 1 := by rw [← Real.exp_add]; simp
    linear_combination (-1 / 2 : ℝ) * this - (1 / 2 : ℝ) * h2
  have hcp : 0 < Real.cosh (y / 2) := Real.cosh_pos _
  rw [hc, Real.log_mul (Real.exp_pos _).ne' hcp.ne', Real.log_exp, add_sub_cancel_left]
  have h1 : 0 ≤ Real.log (Real.cosh (y / 2)) := Real.log_nonneg (Real.one_le_cosh _)
  have h2 : Real.log (Real.cosh (y / 2)) ≤ y ^ 2 / 8 := by
    have := Real.cosh_le_exp_half_sq (y / 2)
    have := Real.log_le_log hcp this
    rw [Real.log_exp] at this
    linarith [show (y / 2) ^ 2 / 2 = y ^ 2 / 8 by ring]
  rw [abs_of_nonneg h1]; exact h2

end R03

namespace R03

open Filter Topology

lemma Efun_zero (N : ℕ) : Efun N 0 = 0 := by simp [Efun, Q01.ell_zero]

noncomputable def Lg (N : ℕ) (x : ℝ) : ℝ := Real.log ((1 + Real.exp (Efun N x)) / 2)

lemma Lg_an (N : ℕ) : AnalyticAt ℝ (Lg N) 0 := by
  have h : AnalyticAt ℝ (fun x => (1 + Real.exp (Efun N x)) / 2) 0 :=
    (analyticAt_const.add (Efun_an N).rexp').div_const
  exact h.log (by simp [Efun_zero])

lemma Efun_sq (N : ℕ) : ∃ D δ : ℝ, 0 < δ ∧ ∀ x : ℝ, |x| < δ → |Efun N x| ≤ D * x ^ 2 := by
  obtain ⟨δ, hδ, D, hD⟩ := Q01.taylor_bd (Efun_an N) 2
  refine ⟨D, δ, hδ, fun x hx => ?_⟩
  have h := hD x hx
  have h0 : iteratedDeriv 0 (Efun N) 0 / (Nat.factorial 0 : ℝ) = 0 := by simp [Efun_zero]
  have h1 : iteratedDeriv 1 (Efun N) 0 / (Nat.factorial 1 : ℝ) = 0 := by
    have := Efun_cf N 1
    unfold cf at this
    rw [this]; simp [Q01.cz]
  simp only [Finset.sum_range_succ, Finset.sum_range_zero, h0, h1, zero_mul, zero_add,
    sub_zero] at h
  rwa [sq_abs] at h

lemma Lg_cf {N : ℕ} (k : ℕ) (hk : k < 4) : cf (Lg N) k = cf (Efun N) k / 2 := by
  obtain ⟨D, δ, hδ, hE⟩ := Efun_sq N
  have han : AnalyticAt ℝ (fun x => Lg N x - (1 / 2 : ℝ) * Efun N x) 0 :=
    (Lg_an N).sub (analyticAt_const.mul (Efun_an N))
  have hv := coeff_vanish han (C := D ^ 2 / 8) hδ (fun x hx => by
    have h1 := logcosh (Efun N x)
    have h2 := hE x hx
    have h3 : Efun N x ^ 2 ≤ (D * x ^ 2) ^ 2 := by
      rw [← sq_abs]; exact pow_le_pow_left₀ (abs_nonneg _) h2 2
    unfold Lg
    calc |Real.log ((1 + Real.exp (Efun N x)) / 2) - 1 / 2 * Efun N x|
        = |Real.log ((1 + Real.exp (Efun N x)) / 2) - Efun N x / 2| := by ring_nf
      _ ≤ Efun N x ^ 2 / 8 := h1
      _ ≤ (D * x ^ 2) ^ 2 / 8 := by gcongr
      _ = D ^ 2 / 8 * |x| ^ 4 := by rw [← sq_abs x]; ring) k hk
  have : cf (fun x => Lg N x - (1 / 2 : ℝ) * Efun N x) k = 0 := by unfold cf; rw [hv]; simp
  rw [cf_sub (f := Lg N) (g := fun x => (1 / 2 : ℝ) * Efun N x) (Lg_an N)
    (analyticAt_const.mul (Efun_an N)), cf_smul (1 / 2 : ℝ) (Efun N)] at this
  linarith

lemma binomHalf01 : binomHalf 0 = 1 ∧ binomHalf 1 = 1 / 2 := by
  constructor <;> simp [binomHalf]

/-- `log G_3 = log Φ + log((1+e^E)/2) + O(τ⁴)`. -/
lemma logG3_bd {N : ℕ} (hN : 2 ≤ N) :
    ∃ C r : ℝ, 0 < r ∧ ∀ x : ℝ, |x| < r →
      |Real.log (Gm N 3 x) - (Afun N x + Lg N x)| ≤ C * |x| ^ 4 := by
  have hN1 : 1 ≤ N := by omega
  have hN' : (0 : ℝ) < N := by exact_mod_cast hN1
  obtain ⟨C2, r2, hr2, h2⟩ := momX_bd2 hN 2
  obtain ⟨C3, r3, hr3, h3⟩ := momX_bd2 hN 3
  obtain ⟨D, δ, hδ, hE⟩ := Efun_sq N
  set Cr := |binomHalf 2| * |C2| + |binomHalf 3| * |C3| with hCr
  set r0 : ℝ := 1 / (2 * N * (1 + 2 * 1)) with hr0
  set ρ := min (min (min r0 r2) (min r3 δ)) (min 1 (min (1 / (2 * (|D| + 1))) (1 / (4 * (Cr + 1)))))
    with hρ
  have hρ0 : 0 < ρ := by
    have : 0 ≤ Cr := by positivity
    positivity
  refine ⟨2 * Cr, ρ, hρ0, fun x hx => ?_⟩
  have hx0 : |x| < r0 := lt_of_lt_of_le hx ((min_le_left _ _).trans ((min_le_left _ _).trans
    (min_le_left _ _)))
  have hx2 : |x| < r2 := lt_of_lt_of_le hx ((min_le_left _ _).trans ((min_le_left _ _).trans
    (min_le_right _ _)))
  have hx3 : |x| < r3 := lt_of_lt_of_le hx ((min_le_left _ _).trans ((min_le_right _ _).trans
    (min_le_left _ _)))
  have hxδ : |x| < δ := lt_of_lt_of_le hx ((min_le_left _ _).trans ((min_le_right _ _).trans
    (min_le_right _ _)))
  have hx1 : |x| ≤ 1 := hx.le.trans ((min_le_right _ _).trans (min_le_left _ _))
  have hxD : |x| ≤ 1 / (2 * (|D| + 1)) :=
    hx.le.trans ((min_le_right _ _).trans ((min_le_right _ _).trans (min_le_left _ _)))
  have hxC : |x| ≤ 1 / (4 * (Cr + 1)) :=
    hx.le.trans ((min_le_right _ _).trans ((min_le_right _ _).trans (min_le_right _ _)))
  -- the moments
  have hm0 : momX N 0 x = 1 := by
    rw [momX_eq hN1 0 (by
      refine hx0.le.trans ?_
      rw [hr0]; gcongr; norm_num), Fm0]
  have hm1 : momX N 1 x = Real.exp (Efun N x) - 1 := by
    rw [momX_eq hN1 1 (by rw [hr0] at hx0; push_cast; linarith), Fm1 hN1]
  have hx4 : |x| ^ 4 ≤ |x| := by
    calc |x| ^ 4 ≤ |x| ^ 1 := pow_le_pow_of_le_one (abs_nonneg x) hx1 (by norm_num)
      _ = |x| := pow_one _
  have hx2' : x ^ 2 ≤ |x| := by
    rw [← sq_abs]; nlinarith [abs_nonneg x]
  set R := binomHalf 2 * momX N 2 x + binomHalf 3 * momX N 3 x with hR
  have hRb : |R| ≤ Cr * |x| ^ 4 := by
    have e2 := h2 x hx2
    have e3 := h3 x hx3
    have e2' : |momX N 2 x| ≤ |C2| * |x| ^ 4 := by
      refine e2.trans ?_
      rw [show (x ^ 2) ^ 2 = |x| ^ 4 by rw [← sq_abs x]; ring]
      gcongr; exact le_abs_self _
    have e3' : |momX N 3 x| ≤ |C3| * |x| ^ 4 := by
      have h6 : (x ^ 2) ^ 3 = |x| ^ 6 := by rw [← sq_abs x]; ring
      have h64 : |x| ^ 6 ≤ |x| ^ 4 := pow_le_pow_of_le_one (abs_nonneg x) hx1 (by norm_num)
      rw [h6] at e3
      calc |momX N 3 x| ≤ C3 * |x| ^ 6 := e3
        _ ≤ |C3| * |x| ^ 6 := mul_le_mul_of_nonneg_right (le_abs_self _) (by positivity)
        _ ≤ |C3| * |x| ^ 4 := mul_le_mul_of_nonneg_left h64 (abs_nonneg _)
    calc |R| ≤ |binomHalf 2 * momX N 2 x| + |binomHalf 3 * momX N 3 x| := abs_add_le _ _
      _ = |binomHalf 2| * |momX N 2 x| + |binomHalf 3| * |momX N 3 x| := by
          rw [abs_mul, abs_mul]
      _ ≤ |binomHalf 2| * (|C2| * |x| ^ 4) + |binomHalf 3| * (|C3| * |x| ^ 4) := by gcongr
      _ = Cr * |x| ^ 4 := by rw [hCr]; ring
  have hB : Bpoly N 3 x = (1 + Real.exp (Efun N x)) / 2 + R := by
    unfold Bpoly
    simp only [Finset.sum_range_succ, Finset.sum_range_zero, zero_add]
    rw [hm0, hm1, binomHalf01.1, binomHalf01.2, hR]
    ring
  -- lower bounds
  have hEx := hE x hxδ
  have hEs : |Efun N x| ≤ 1 / 2 := by
    refine hEx.trans ?_
    have : D * x ^ 2 ≤ |D| * |x| := by
      calc D * x ^ 2 ≤ |D| * x ^ 2 := by gcongr; exact le_abs_self _
        _ ≤ |D| * |x| := by gcongr
    refine this.trans ?_
    have : |D| * |x| ≤ |D| * (1 / (2 * (|D| + 1))) := by gcongr
    refine this.trans ?_
    rw [mul_one_div, div_le_div_iff₀ (by positivity) (by norm_num)]
    nlinarith [abs_nonneg D]
  have hexp : 1 / 2 ≤ Real.exp (Efun N x) := by
    have := Real.add_one_le_exp (Efun N x)
    have := (abs_le.mp hEs).1
    linarith
  have hRs : |R| ≤ 1 / 4 := by
    refine hRb.trans ?_
    have : Cr * |x| ^ 4 ≤ Cr * |x| := by gcongr
    refine this.trans ?_
    have hCr0 : 0 ≤ Cr := by positivity
    have : Cr * |x| ≤ Cr * (1 / (4 * (Cr + 1))) := by gcongr
    refine this.trans ?_
    rw [mul_one_div, div_le_div_iff₀ (by positivity) (by norm_num)]
    nlinarith
  have hL1 : 1 / 2 ≤ (1 + Real.exp (Efun N x)) / 2 := by linarith [Real.exp_pos (Efun N x)]
  have hL2 : 1 / 2 ≤ Bpoly N 3 x := by rw [hB]; linarith [(abs_le.mp hRs).1]
  have hG : Real.log (Gm N 3 x) = Afun N x + Real.log (Bpoly N 3 x) := by
    unfold Gm
    rw [Phif_exp, Real.log_mul (Real.exp_pos _).ne' (by linarith), Real.log_exp]
  rw [hG, show Afun N x + Real.log (Bpoly N 3 x) - (Afun N x + Lg N x) =
    Real.log (Bpoly N 3 x) - Lg N x by ring]
  unfold Lg
  calc |Real.log (Bpoly N 3 x) - Real.log ((1 + Real.exp (Efun N x)) / 2)|
      ≤ 2 * |Bpoly N 3 x - (1 + Real.exp (Efun N x)) / 2| := Q01.log_lip hL2 hL1
    _ = 2 * |R| := by rw [hB]; ring_nf
    _ ≤ 2 * (Cr * |x| ^ 4) := by gcongr
    _ = _ := by ring

lemma wN_cf {N : ℕ} (hN : 2 ≤ N) {k : ℕ} (hk : k ≤ 3) :
    wN N k = cf (Afun N) k + cf (Efun N) k / 2 := by
  have hst := stable_coeff (f := fun K x => Real.log (Gm N K x))
    (fun K => (pf_log (pf_Gm hN K) (Gm_zero hN K)).1) (fun K => log_step hN K) 3 k hk
  have hG3 : AnalyticAt ℝ (fun x => Real.log (Gm N 3 x)) 0 :=
    (pf_log (pf_Gm hN 3) (Gm_zero hN 3)).1
  have hM : AnalyticAt ℝ (fun x => Afun N x + Lg N x) 0 := (Afun_an N).add (Lg_an N)
  obtain ⟨C, r, hr, hb⟩ := logG3_bd hN
  have hv := coeff_vanish (hG3.sub hM) hr (fun x hx => hb x hx) k (by omega)
  rw [iteratedDeriv_sub hG3.contDiffAt hM.contDiffAt, sub_eq_zero] at hv
  have : wN N k = cf (fun x => Afun N x + Lg N x) k := by
    unfold wN cf; rw [← hst, hv]
  rw [this, cf_add (Afun_an N) (Lg_an N), Lg_cf k (by omega)]

end R03

theorem R03.proof_10 : R03.stmt_10 := by
  intro N hN
  obtain ⟨w, hw0, hw, hK⟩ := R03.main_p hN
  have hN' : (N : ℝ) ≠ 0 := by
    have : (2 : ℝ) ≤ N := by exact_mod_cast hN
    linarith
  have hz2 : Q01.cz 2 = X05.z 2 / 2 := by simp [Q01.cz]
  have hz3 : Q01.cz 3 = -X05.z 3 / 3 := by simp [Q01.cz]; ring
  refine ⟨w, hw, ?_, ?_, hK⟩
  · rw [hw0, R03.wN_cf hN (by norm_num), R03.Afun_cf, R03.Efun_cf, hz2]
    simp only [R03.bq]
    push_cast
    field_simp
    ring
  · rw [hw0, R03.wN_cf hN (by norm_num), R03.Afun_cf, R03.Efun_cf, hz3]
    simp only [R03.bq, T04.coeffFormula]
    push_cast
    field_simp
    ring


namespace R03

open Filter Topology

/-- analytic germs with rational Taylor coefficients. -/
def RQ (f : ℝ → ℝ) : Prop :=
  AnalyticAt ℝ f 0 ∧ ∀ a : ℕ, ∃ q : ℚ, iteratedDeriv a f 0 / (a.factorial : ℝ) = q

lemma choose_cast (r : ℝ) (k : ℕ) :
    Ring.choose r k = (∏ i ∈ Finset.range k, (r - i)) / (k.factorial : ℝ) := by
  have h := Ring.descPochhammer_eq_factorial_smul_choose r k
  rw [← Polynomial.aeval_eq_smeval, Polynomial.aeval_def, Polynomial.eval₂_eq_eval_map,
    descPochhammer_map, descPochhammer_eval_eq_prod_range, nsmul_eq_mul] at h
  rw [h]; field_simp

lemma rq_binom (r : ℚ) : RQ (fun x => (1 + x) ^ (r : ℝ)) := by
  have hB := Real.one_add_rpow_hasFPowerSeriesOnBall_zero (a := (r : ℝ))
  have hA : AnalyticAt ℝ (fun x : ℝ => (1 + x) ^ (r : ℝ)) 0 := hB.analyticAt
  refine ⟨hA, fun a => ?_⟩
  have h2 := hA.hasFPowerSeriesAt.eq_formalMultilinearSeries hB.hasFPowerSeriesAt
  have h3 := congrArg (fun p : FormalMultilinearSeries ℝ ℝ ℝ => p a (fun _ => 1)) h2
  simp only [FormalMultilinearSeries.ofScalars_apply_eq, one_pow, smul_eq_mul, mul_one,
    binomialSeries_apply, List.ofFn_const, List.prod_replicate] at h3
  refine ⟨(∏ i ∈ Finset.range a, (r - i)) / (a.factorial : ℚ), ?_⟩
  rw [h3, choose_cast]
  push_cast; ring

lemma rq_scale {f : ℝ → ℝ} (hf : RQ f) (q : ℚ) : RQ (fun x => f ((q : ℝ) * x)) := by
  set u : ℝ →L[ℝ] ℝ := (q : ℝ) • ContinuousLinearMap.id ℝ ℝ with hu
  have hg : (fun x => f ((q : ℝ) * x)) = f ∘ u := by funext x; simp [hu]
  have hA : AnalyticAt ℝ (fun x => f ((q : ℝ) * x)) 0 := by
    rw [hg]; exact (by simpa [hu] using hf.1 : AnalyticAt ℝ f (u 0)).comp (u.analyticAt 0)
  refine ⟨hA, fun a => ?_⟩
  have h1 : HasFPowerSeriesAt (f ∘ u) ((FormalMultilinearSeries.ofScalars ℝ
      (fun n => iteratedDeriv n f 0 / (n.factorial : ℝ))).compContinuousLinearMap u) 0 :=
    HasFPowerSeriesAt.compContinuousLinearMap (by simpa [hu] using hf.1.hasFPowerSeriesAt)
  rw [← hg] at h1
  have h2 := hA.hasFPowerSeriesAt.eq_formalMultilinearSeries h1
  have h3 := congrArg (fun p : FormalMultilinearSeries ℝ ℝ ℝ => p a (fun _ => 1)) h2
  simp only [FormalMultilinearSeries.ofScalars_apply_eq, one_pow, smul_eq_mul, mul_one,
    FormalMultilinearSeries.compContinuousLinearMap_apply] at h3
  rw [h3]
  have h4 : (fun _ : Fin a => u 1) = fun _ => (q : ℝ) := by funext; simp [hu]
  rw [Function.comp_def, h4, FormalMultilinearSeries.ofScalars_apply_eq, smul_eq_mul]
  obtain ⟨c, hc⟩ := hf.2 a
  refine ⟨c * q ^ a, ?_⟩
  rw [hc]; push_cast; ring

lemma rq_add {f g : ℝ → ℝ} (hf : RQ f) (hg : RQ g) : RQ (fun x => f x + g x) := by
  refine ⟨hf.1.add hg.1, fun a => ?_⟩
  have := iteratedDeriv_add (n := a) hf.1.contDiffAt hg.1.contDiffAt
  rw [show (fun x => f x + g x) = f + g from rfl, this, add_div]
  obtain ⟨c, hc⟩ := hf.2 a; obtain ⟨d, hd⟩ := hg.2 a
  exact ⟨c + d, by rw [hc, hd]; push_cast; ring⟩

lemma rq_smul {f : ℝ → ℝ} (c : ℚ) (hf : RQ f) : RQ (fun x => (c : ℝ) * f x) := by
  refine ⟨analyticAt_const.mul hf.1, fun a => ?_⟩
  rw [iteratedDeriv_const_mul_field, mul_div_assoc]
  obtain ⟨d, hd⟩ := hf.2 a
  exact ⟨c * d, by rw [hd]; push_cast; ring⟩

lemma rq_mul {f g : ℝ → ℝ} (hf : RQ f) (hg : RQ g) : RQ (fun x => f x * g x) := by
  refine ⟨hf.1.mul hg.1, fun a => ?_⟩
  have := iteratedDeriv_mul (n := a) hf.1.contDiffAt hg.1.contDiffAt
  rw [show (fun x => f x * g x) = f * g from rfl, this, Finset.sum_div]
  choose qf hqf using hf.2
  choose qg hqg using hg.2
  refine ⟨∑ i ∈ Finset.range (a + 1), qf i * qg (a - i), ?_⟩
  push_cast
  refine Finset.sum_congr rfl (fun i hi => ?_)
  have hi' : i ≤ a := Nat.lt_succ_iff.1 (Finset.mem_range.1 hi)
  rw [← hqf i, ← hqg (a - i), Nat.cast_choose ℝ hi']
  have h1 : (i.factorial : ℝ) ≠ 0 := by positivity
  have h2 : ((a - i).factorial : ℝ) ≠ 0 := by positivity
  have h3 : (a.factorial : ℝ) ≠ 0 := by positivity
  field_simp

lemma rq_pow {f : ℝ → ℝ} (hf : RQ f) : ∀ n : ℕ, RQ (fun x => f x ^ n)
  | 0 => by
    have := rq_binom 0
    simpa using this
  | n + 1 => by
    simp only [pow_succ]
    exact rq_mul (rq_pow hf n) hf

/-- `Ψ_n(d) = (1-d)^{-1} (1 - (1+d)^{-1})^n`, i.e. `(d/(1+d))^n/(1-d)`. -/
noncomputable def Psi (n : ℕ) (d : ℝ) : ℝ :=
  (1 + ((-1 : ℚ) : ℝ) * d) ^ ((-1 : ℚ) : ℝ) *
    ((1 + d) ^ ((0 : ℚ) : ℝ) + ((-1 : ℚ) : ℝ) * (1 + d) ^ ((-1 : ℚ) : ℝ)) ^ n

lemma rq_Psi (n : ℕ) : RQ (Psi n) := by
  unfold Psi
  refine rq_mul (rq_scale (rq_binom (-1)) (-1)) (rq_pow ?_ n)
  exact rq_add (rq_binom 0) (rq_smul (-1) (rq_binom (-1)))

lemma Psi_eq (n : ℕ) {d : ℝ} (hd : |d| < 1) :
    Psi n d = (1 / (1 - d)) * (d / (1 + d)) ^ n := by
  have h1 : 0 < 1 + d := by cases abs_lt.mp hd; linarith
  have h2 : 0 < 1 - d := by cases abs_lt.mp hd; linarith
  unfold Psi
  push_cast
  rw [Real.rpow_zero, Real.rpow_neg_one, Real.rpow_neg_one]
  congr 1
  · rw [show 1 + -1 * d = 1 - d by ring, one_div]
  · congr 1; field_simp; ring

/-- the `d`-model `H_K(d) = Σ_{n≤K} r_n Ψ_n(d)`. -/
noncomputable def Hd (r : ℕ → ℝ) (K : ℕ) (d : ℝ) : ℝ :=
  ∑ n ∈ Finset.range (K + 1), r n * Psi n d

lemma Hd_an (r : ℕ → ℝ) (K : ℕ) : AnalyticAt ℝ (Hd r K) 0 := by
  unfold Hd
  have : (fun d => ∑ n ∈ Finset.range (K + 1), r n * Psi n d) =
      ∑ n ∈ Finset.range (K + 1), fun d => r n * Psi n d := by
    funext d; simp
  rw [this]
  exact Finset.analyticAt_sum _ (fun n _ => analyticAt_const.mul (rq_Psi n).1)

lemma Hd_coeff (r : ℕ → ℝ) (hr : ∀ n, Q01.IsZetaPoly n (r n)) (m : ℕ) :
    IsZetaPolyLE m (iteratedDeriv m (Hd r m) 0 / (m.factorial : ℝ)) := by
  choose q hq using fun n => (rq_Psi n).2 m
  refine ⟨fun n => if n ≤ m then (q n : ℝ) * r n else 0, fun j hj => ?_, ?_⟩
  · simp only [hj, ↓reduceIte]; exact Q01.zp_smul _ (hr j)
  · have hsum : iteratedDeriv m (Hd r m) 0 =
        ∑ n ∈ Finset.range (m + 1), r n * iteratedDeriv m (Psi n) 0 := by
      unfold Hd
      rw [iteratedDeriv_fun_sum (f := fun n d => r n * Psi n d)
        (fun n _ => contDiffAt_const.mul (rq_Psi n).1.contDiffAt)]
      refine Finset.sum_congr rfl (fun n _ => ?_)
      exact iteratedDeriv_const_mul_field _ _
    rw [hsum, Finset.sum_div]
    refine Finset.sum_congr rfl (fun n hn => ?_)
    have : n ≤ m := Nat.lt_succ_iff.mp (Finset.mem_range.mp hn)
    simp only [this, ↓reduceIte]
    rw [mul_div_assoc, hq n]; ring

lemma Hd_step (r : ℕ → ℝ) (K : ℕ) :
    ∃ C ρ : ℝ, 0 < ρ ∧ ∀ x : ℝ, |x| < ρ → |Hd r (K + 1) x - Hd r K x| ≤ C * |x| ^ (K + 1) := by
  refine ⟨2 * 2 ^ (K + 1) * |r (K + 1)|, 1 / 2, by norm_num, fun x hx => ?_⟩
  have hx1 : |x| < 1 := by linarith
  have h1 : 0 < 1 + x := by cases abs_lt.mp hx1; linarith
  have h2 : 1 / 2 < 1 - x := by cases abs_lt.mp hx; linarith
  have h3 : 1 / 2 < 1 + x := by cases abs_lt.mp hx; linarith
  have hdiff : Hd r (K + 1) x - Hd r K x = r (K + 1) * Psi (K + 1) x := by
    unfold Hd; rw [Finset.sum_range_succ _ (K + 1)]; ring
  rw [hdiff, Psi_eq _ hx1, abs_mul, abs_mul, abs_pow, abs_div, abs_div, abs_of_pos h1,
    abs_of_pos (by linarith : (0 : ℝ) < 1 - x), abs_one]
  have e1 : 1 / (1 - x) ≤ 2 := by rw [div_le_iff₀ (by linarith)]; linarith
  have e2 : |x| / (1 + x) ≤ 2 * |x| := by rw [div_le_iff₀ h1]; nlinarith [abs_nonneg x]
  calc |r (K + 1)| * (1 / (1 - x) * (|x| / (1 + x)) ^ (K + 1))
      ≤ |r (K + 1)| * (2 * (2 * |x|) ^ (K + 1)) := by
        gcongr
    _ = _ := by rw [mul_pow]; ring

end R03

theorem R03.proof_11 : R03.stmt_11 := by
  have h3 : (2 : ℕ) ≤ 3 := by norm_num
  have hS1 : R03.S 3 1 = 4 * Real.sqrt 3 := by
    rw [R03.proof_2 3 h3]; norm_num; ring
  have hY : ∀ p : ℝ, 0 < p → p ≠ 2 → T04.Y 3 (1 - 1 / p) = 1 / (2 - p) :=
    fun p hp hp2 => T04.proof_2 p hp.ne' hp2
  constructor
  · -- (i)
    obtain ⟨w, -, hw, hK⟩ := R03.main_p h3
    refine ⟨w, hw, fun K => ?_⟩
    obtain ⟨C, τ₀, hτ₀, hb⟩ := hK K
    refine ⟨C, min τ₀ (1 / 4), lt_min hτ₀ (by norm_num), fun p hp hτ => ?_⟩
    have hp2 : p ≠ 2 := by
      intro h; rw [h] at hτ
      have := hτ.trans (min_le_right _ _); norm_num [abs_of_pos] at this
    have := hb p hp (hτ.trans (min_le_left _ _))
    rw [hY p hp hp2, hS1, (R03.proof_3).2 p] at this
    convert this using 3
    have h2p : (2 : ℝ) - p ≠ 0 := sub_ne_zero.mpr (Ne.symm hp2)
    field_simp
  · -- (ii)
    set G : ℕ → ℝ → ℝ := fun K => R03.Gm 3 K with hGdef
    have hstG := R03.stable_coeff (f := G) (fun K => (R03.pf_Gm h3 K).1) (fun K => R03.Gm_step h3 K)
    set r : ℕ → ℝ := fun n => iteratedDeriv n (G n) 0 / (n.factorial : ℝ) with hrdef
    have hr : ∀ n, Q01.IsZetaPoly n (r n) := fun n => by simpa using (R03.pf_Gm h3 n).2 n
    have hstH := R03.stable_coeff (f := R03.Hd r) (fun K => R03.Hd_an r K) (fun K => R03.Hd_step r K)
    refine ⟨fun m => iteratedDeriv m (R03.Hd r m) 0 / (m.factorial : ℝ),
      fun m => R03.Hd_coeff r hr m, fun K => ?_⟩
    obtain ⟨C1, r1, hr1, hb1⟩ := R03.PB_bd h3 K
    obtain ⟨δ, hδ, D, hD⟩ := Q01.taylor_bd (R03.pf_Gm h3 K).1 (K + 1)
    obtain ⟨δ', hδ', D', hD'⟩ := Q01.taylor_bd (R03.Hd_an r K) (K + 1)
    set δ0 := min (min (r1 / 4) (δ / 4)) (min (δ' / 2) (1 / 16)) with hδ0
    refine ⟨4 * Real.sqrt 3 * (2 * 2 ^ (K + 1) * (|C1| + |D|) + |D'|), δ0, by positivity,
      fun d hd => ?_⟩
    have hd1 : |d| ≤ 1 / 16 := hd.trans ((min_le_right _ _).trans (min_le_right _ _))
    have hdr : |d| ≤ r1 / 4 := hd.trans ((min_le_left _ _).trans (min_le_left _ _))
    have hdδ : |d| ≤ δ / 4 := hd.trans ((min_le_left _ _).trans (min_le_right _ _))
    have hdδ' : |d| < δ' := by
      have := hd.trans ((min_le_right _ _).trans (min_le_left _ _)); linarith
    have hdab := abs_le.mp hd1
    have hp : 0 < 1 + d := by linarith
    set τ := 1 - 1 / (1 + d) with hτdef
    have hτd : τ = d / (1 + d) := by rw [hτdef]; field_simp; ring
    have hτabs : |τ| ≤ 2 * |d| := by
      rw [hτd, abs_div, abs_of_pos hp, div_le_iff₀ hp]; nlinarith [abs_nonneg d]
    have hτr : |τ| < r1 := by
      have : 0 < r1 := hr1; linarith
    have hτδ : |τ| < δ := by linarith
    have h3τ : ((3 : ℕ) : ℝ) * |τ| ≤ 1 / 2 := by push_cast; linarith
    have hp2 : 1 + d ≠ 2 := by intro h; linarith
    -- the factorisation
    have hSd : Y05.S (1 + d) = 4 * Real.sqrt 3 * (1 / (1 - d)) * (R03.Phif 3 τ * R03.B 3 τ) := by
      have hr := R03.ratio_eq h3 hp h3τ
      rw [← hτdef, hY _ hp hp2, hS1, (R03.proof_3).2] at hr
      have h2p : (2 : ℝ) - (1 + d) ≠ 0 := by intro h; linarith
      have h1d : (1 : ℝ) - d ≠ 0 := by intro h; linarith
      have hs3 : 0 < Real.sqrt 3 := by positivity
      rw [← hr, show (2 : ℝ) - (1 + d) = 1 - d by ring]
      field_simp
    -- `G_K(τ)` and its Taylor polynomial
    have hTay := hD τ hτδ
    have hsum : ∑ k ∈ Finset.range (K + 1), iteratedDeriv k (R03.Gm 3 K) 0 / (k.factorial : ℝ) * τ ^ k
        = ∑ k ∈ Finset.range (K + 1), r k * τ ^ k := by
      refine Finset.sum_congr rfl (fun k hk => ?_)
      have := hstG K k (Nat.lt_succ_iff.mp (Finset.mem_range.mp hk))
      simp only [hGdef] at this
      rw [hrdef]; simp only [hGdef]; rw [this]
    rw [hsum] at hTay
    have hHd : R03.Hd r K d = (1 / (1 - d)) * ∑ k ∈ Finset.range (K + 1), r k * τ ^ k := by
      unfold R03.Hd
      rw [Finset.mul_sum]
      refine Finset.sum_congr rfl (fun k _ => ?_)
      rw [R03.Psi_eq k (by linarith [abs_nonneg d] : |d| < 1), hτd]; ring
    have hHsum : ∑ k ∈ Finset.range (K + 1), iteratedDeriv k (R03.Hd r K) 0 / (k.factorial : ℝ) * d ^ k
        = ∑ k ∈ Finset.range (K + 1), iteratedDeriv k (R03.Hd r k) 0 / (k.factorial : ℝ) * d ^ k := by
      refine Finset.sum_congr rfl (fun k hk => ?_)
      rw [hstH K k (Nat.lt_succ_iff.mp (Finset.mem_range.mp hk))]
    have hT2 := hD' d hdδ'
    rw [hHsum] at hT2
    set P := ∑ k ∈ Finset.range (K + 1), iteratedDeriv k (R03.Hd r k) 0 / (k.factorial : ℝ) * d ^ k
    have hpre : |R03.Phif 3 τ * R03.B 3 τ - ∑ k ∈ Finset.range (K + 1), r k * τ ^ k| ≤
        (|C1| + |D|) * |τ| ^ (K + 1) := by
      have e1 := hb1 τ hτr
      calc _ ≤ |R03.Phif 3 τ * R03.B 3 τ - R03.Gm 3 K τ| +
            |R03.Gm 3 K τ - ∑ k ∈ Finset.range (K + 1), r k * τ ^ k| := abs_sub_le _ _ _
        _ ≤ |C1| * |τ| ^ (K + 1) + |D| * |τ| ^ (K + 1) := by
            gcongr
            · exact e1.trans (by gcongr; exact le_abs_self _)
            · exact hTay.trans (by gcongr; exact le_abs_self _)
        _ = _ := by ring
    have hτpow : |τ| ^ (K + 1) ≤ 2 ^ (K + 1) * |d| ^ (K + 1) := by
      rw [← mul_pow]; exact pow_le_pow_left₀ (abs_nonneg _) hτabs _
    have h1d : 0 < 1 - d := by linarith
    have hinv : |1 / (1 - d)| ≤ 2 := by
      rw [abs_of_pos (by positivity), div_le_iff₀ h1d]; linarith
    have hmain : |Y05.S (1 + d) / (4 * Real.sqrt 3) - R03.Hd r K d| ≤
        2 * 2 ^ (K + 1) * (|C1| + |D|) * |d| ^ (K + 1) := by
      have hs3 : 0 < Real.sqrt 3 := by positivity
      rw [hSd, hHd, show 4 * Real.sqrt 3 * (1 / (1 - d)) * (R03.Phif 3 τ * R03.B 3 τ) /
          (4 * Real.sqrt 3) = 1 / (1 - d) * (R03.Phif 3 τ * R03.B 3 τ) by field_simp,
        ← mul_sub, abs_mul]
      calc |1 / (1 - d)| * |R03.Phif 3 τ * R03.B 3 τ - ∑ k ∈ Finset.range (K + 1), r k * τ ^ k|
          ≤ 2 * ((|C1| + |D|) * (2 ^ (K + 1) * |d| ^ (K + 1))) := by
            gcongr
            exact hpre.trans (by gcongr)
        _ = _ := by ring
    have hs3 : 0 < Real.sqrt 3 := by positivity
    have hfin : |Y05.S (1 + d) / (4 * Real.sqrt 3) - P| ≤
        (2 * 2 ^ (K + 1) * (|C1| + |D|) + |D'|) * |d| ^ (K + 1) := by
      calc _ ≤ |Y05.S (1 + d) / (4 * Real.sqrt 3) - R03.Hd r K d| + |R03.Hd r K d - P| :=
            abs_sub_le _ _ _
        _ ≤ 2 * 2 ^ (K + 1) * (|C1| + |D|) * |d| ^ (K + 1) + |D'| * |d| ^ (K + 1) := by
            gcongr
            exact hT2.trans (by gcongr; exact le_abs_self _)
        _ = _ := by ring
    have : Y05.S (1 + d) - 4 * Real.sqrt 3 * P =
        4 * Real.sqrt 3 * (Y05.S (1 + d) / (4 * Real.sqrt 3) - P) := by field_simp
    rw [this, abs_mul, abs_of_pos (by positivity : (0 : ℝ) < 4 * Real.sqrt 3)]
    calc 4 * Real.sqrt 3 * |Y05.S (1 + d) / (4 * Real.sqrt 3) - P|
        ≤ 4 * Real.sqrt 3 * ((2 * 2 ^ (K + 1) * (|C1| + |D|) + |D'|) * |d| ^ (K + 1)) :=
          mul_le_mul_of_nonneg_left hfin (by positivity)
      _ = _ := by ring
