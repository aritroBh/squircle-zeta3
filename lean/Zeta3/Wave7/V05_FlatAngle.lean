import Mathlib
import Zeta3.Wave2.D04_MalphaForms
/-!
# V05 -- expanding the corner coefficients `M_k(α)` about the flat angle (first orders)

Target id: V05.
Paper label: `thm:flat` (paper/ieee/main.tex, lines 1501-1511; proof idea 1513-1515).
Supplement: paper/ieee/proofs/sec_g6_flat.tex (`g6:thm:flat`, lines 13-20; Lemmas
`g6:lem:series` 36-48, `g6:lem:En` 82-88, `g6:lem:gamma` 110-112, Definition `g6:def:formal`
122-127, Lemma `g6:prop:formula` 129-131, Lemma `g6:lem:first` 139-146; proof of the theorem
180-186). The corner integral is defined in main.tex lines 1466-1471 (reused from
`Zeta3.Wave2.D04_MalphaForms`, `D04.Mint`, and `D04.c α = cos(α/2)`).

Quoted verbatim (main.tex 1466-1471):

  its length is $2+2u\,M(u;\alpha)$, where
  \begin{multline*}
  M(u;\alpha)=\int_0^1\Bigl[(1+g)^{-1-u}\sqrt{1+g^{2-2u}-2\cos\alpha\,g^{1-u}}\\
  -1\Bigr]g^{u-1}\dd g=\sum_{k\ge0}M_k(\alpha)u^k,
  \end{multline*}

Quoted verbatim (main.tex 1501-1511, `thm:flat`):

  \begin{theorem}[Expanding about the flat angle$^\dagger$]\label{thm:flat}
  For every $k\ge0$ and $0<\alpha<\pi$, $M_k(\alpha)=\sum_{n\ge1}\mu_{k,n}c^{2n}$, where each
  $\mu_{k,n}$ is a polynomial with rational coefficients in $\zeta(2),\dots,\zeta(k)$ of weight at
  most $k$. In particular, with $r_n,r_n',q_n'\in\mathbb Q$,
  \[
  \mu_{2,n}=r_n+\frac{\zeta(2)}{4n^2-1},\qquad
  \mu_{3,n}=r_n'+q_n'\,\zeta(2)-\frac{\zeta(3)}{4n^2-1}.
  \]
  So, with $\Phi(c)=\sum_{n\ge1}\frac{c^{2n}}{4n^2-1}=\frac12-\frac{1-c^2}{2c}\operatorname{artanh}c$,
  the functions $M_2-\zeta(2)\Phi$ and $M_3+\zeta(3)\Phi-\zeta(2)\Psi$ are power series in $c^2$
  with rational coefficients, for some $\Psi\in\mathbb Q[[c^2]]$.
  \end{theorem}

Quoted verbatim (supplement, line 10): "The numbers $M_k(\alpha)$ are its Taylor coefficients,
$M_k(\alpha)=M^{(k)}(0;\alpha)/k!$. (We never use that the series converges.)"

Quoted verbatim (supplement, Lemma `g6:lem:series`, lines 30-47):

  Put $a_n=\binom{2n}{n}/(2n-1)$ for $n\ge1$. For $0<y<1$ write $a=\log y$, $b=\log(1-y)$
  (both negative), $\ell=b-a$, $A=-(n-1)a+nb$, and for real $u$
  \[ f(u,y)=1-y+y\,\eu^{u\ell},\qquad E_n(u,y)=\eu^{uA}\,f(u,y)^{1-2n}. \]
  \item[(a)] For $|u|<1$, $M(u;\alpha)=-\sum_{n\ge1}a_nc^{2n}J_n(u)$ with
  $J_n(u)=\int_0^1(1+g)^{-1-u}\,h^{n-1}(1+h)^{1-2n}\dd g$, $h=g^{1-u}$, and all terms are positive.
  \item[(b)] $J_n(u)=\frac12\int_0^1[y(1-y)]^{n-1}E_n(u,y)\dd y$.
  \item[(c)] There is $\delta\in(0,\frac1{12}]$, depending only on $c$, such that
  $M(\cdot\,;\alpha)$ is $C^\infty$ on $[-\delta,\delta]$ and for every $k\ge0$
  \[ M_k(\alpha)=-\sum_{n\ge1}a_nc^{2n}\,\frac{J_n^{(k)}(0)}{k!},\qquad
  \frac{J_n^{(k)}(0)}{k!}=\frac12\int_0^1[y(1-y)]^{n-1}\,\frac{\partial_u^kE_n(0,y)}{k!}\dd y , \]
  where the series converges absolutely.

Quoted verbatim (supplement, Lemma `g6:lem:En`, lines 83-87):

  For $n\ge1$ and $k\ge0$ there is a polynomial $E_{n,k}(y,a,b)=\sum e^{(k)}_{j,i,m}(n)\,y^ja^ib^m$,
  with coefficients $e^{(k)}_{j,i,m}\in\QQ[n]$, such that $\partial_u^kE_n(0,y)/k!=
  E_{n,k}(y,\log y,\log(1-y))$. Here $j\le k$ and $i+m=k$ in every monomial. Moreover, with
  $X=\alpha_1a+\beta_1b$, $\alpha_1=-(n-1)+(2n-1)y$, $\beta_1=n-(2n-1)y$, and
  $\gamma=-\frac12(2n-1)y(1-y)$,
  E_{n,0}=1,\quad E_{n,1}=X,\quad E_{n,2}=\tfrac12X^2+\gamma\,\ell^2,\quad
  E_{n,3}=\tfrac16X^3+\gamma\,X\ell^2+\gamma_3\,\ell^3,
  where $\ell=b-a$ and $\gamma_3=-\frac16(2n-1)y(1-y)(1-2y)$.

Quoted verbatim (supplement, Lemma `g6:lem:gamma`, line 111):

  For integers $p,q\ge1$ and $i,m\ge0$, $D_{i,m}(p,q):=\partial_p^i\partial_q^m\gBeta(p,q)/\gBeta(p,q)$
  is the value of $D^\sharp_{i,m}(p,q)$ at $\gvs w=\zeta(w)$.

Quoted verbatim (supplement, Definition `g6:def:formal` and Lemma `g6:prop:formula`, 123-130):

  For $n\ge1$, $k\ge0$ put $\rho_j(n)=(n)_j/(2n)_j$, $(t)_j=t(t+1)\cdots(t+j-1)$, and
  P_{k,n}=-\frac1{n(2n-1)}\sum_{j=0}^k\ \sum_{i+m=k}e^{(k)}_{j,i,m}(n)\,\rho_j(n)\,
  D^\sharp_{i,m}(n+j,n)\ \in\QQ[\gvs2,\dots,\gvs k].
  For $0<\alpha<\pi$ and $k\ge0$, $M_k(\alpha)=\sum_{n\ge1}\mu_{k,n}c^{2n}$ with
  $\mu_{k,n}=P_{k,n}(\zeta(2),\dots,\zeta(k))$, absolutely convergent.

Quoted verbatim (supplement, Lemma `g6:lem:first`, lines 140-145):

  For $n\ge1$ put $D_n=H_{2n-1}-H_{n-1}$ (harmonic numbers $H_N=H^{(1)}_N$). Then
  \item[(a)] $P_{0,n}=-\dfrac1{n(2n-1)}$ and $P_{1,n}=\dfrac{D_n}{n(2n-1)}-\dfrac1{2n^2}$.
  \item[(b)] $P_{2,n}=r_n+\dfrac{\gvs2}{4n^2-1}$ with $r_n\in\QQ$.
  \item[(c)] $P_{3,n}=r'_n+q'_n\gvs2-\dfrac{\gvs3}{4n^2-1}$ with $r'_n,q'_n\in\QQ$.

## Formalization choices

* `ζ(2)` is `z2 := π^2/6` (Mathlib `riemannZeta_two`); `ζ(3)` is `z3 := (riemannZeta 3).re`,
  the real part of Mathlib's complex `riemannZeta 3` (it is real, and
  `zeta_nat_eq_tsum_of_gt_one` identifies `riemannZeta 3` with `∑ 1/n^3`).
* `M(u;α)` is `D04.Mint u α` (imported, not redefined), `c = D04.c α = cos(α/2)`.
  `M_k(α)` is DEFINED as the paper's Taylor coefficient `M^{(k)}(0;α)/k!`
  (`Mk k α := iteratedDeriv k (M(·;α)) 0 / k!`), never as the series. Smoothness of `M(·;α)`
  near `0` is asserted (stmt_4, stmt_12), so `iteratedDeriv` is the true derivative.
* The formal coefficient `P_{k,n}` is a polynomial in formal zeta values; we represent it as a
  real function `P k n x2 x3` of two real variables `x2, x3` standing for `𝔷2, 𝔷3` (enough for
  `k ≤ 3`). "Polynomial with rational coefficients" for `k = 2, 3` is stated as: there are
  rationals `r_n` (resp. `r'_n, q'_n`) with `P 2 n x2 x3 = r_n + x2/(4n²-1)` for ALL real
  `x2, x3` (resp. `P 3 n x2 x3 = r'_n + q'_n x2 - x3/(4n²-1)`). This is the formal (identity in
  the variables) statement of the supplement; nothing is claimed about independence of
  `1, π², ζ(3)`.
* `D^♯_{i,m}(p,q)` (`Dsharp`) is written out explicitly for `i + m ≤ 3`, read off from
  `exp(Δ_⋆)` up to degree 3 (with `A = h₁-h₃`, `B = h₂-h₃`, `P = h₁⁽²⁾-h₃⁽²⁾`,
  `Q = h₂⁽²⁾-h₃⁽²⁾`, `R = h₃⁽²⁾-𝔷2`, `S = h₁⁽³⁾-h₃⁽³⁾`, `T = h₂⁽³⁾-h₃⁽³⁾`, `W = 𝔷3-h₃⁽³⁾`,
  `h₁ = H_{p-1}`, `h₂ = H_{q-1}`, `h₃ = H_{p+q-1}`):
  `D10 = A, D01 = B, D20 = A²-P, D02 = B²-Q, D11 = AB+R, D30 = 2S-3AP+A³, D03 = 2T-3BQ+B³,
  D21 = 2W+2AR-BP+A²B, D12 = 2W+2BR-AQ+AB²`. It is `0` (unused) for `i + m > 3`.
* The coefficients `e^{(k)}_{j,i,m}(n)` are the `y^j`-coefficients of the polynomials
  `ecoef k n i m ∈ ℝ[y]`, which are the coefficients of `a^i b^m` in the displayed `E_{n,k}`
  (expanded from `X`, `γ`, `γ₃`, `ℓ` exactly as displayed; for `k ≤ 3`).
* `B(p,q)` for integers `p, q ≥ 1` is `(p-1)!(q-1)!/(p+q-1)!`; the Beta-log integrals of
  `g6:lem:gamma` are stated with integrability (so no junk value of the interval integral).
* Sums over `n ≥ 1` are written over `n : ℕ` with index `n + 1`. Natural subtraction `n - 1`
  only occurs with `n ≥ 1`. All `rpow` bases are positive on `(0,1]` (resp. `(0,1)`) except the
  endpoint `g = 0`, a null set.
* `C^∞ on [-δ,δ]` (`g6:lem:series`(c)) is stated as `C^k` on the open interval `(-δ,δ)` for
  every `k` (enough for the Taylor coefficients at `0`).

## Statements

* stmt_1: `Φ(c) = Σ_{n≥1} c^{2n}/(4n²-1) = 1/2 - (1-c²)/(2c) artanh c` for `0 < |c| < 1`.
* stmt_2: `g6:lem:series`(a): for `0<α<π`, `|u|<1`, the corner integrand is integrable, every
  `J_n(u) > 0`, and `M(u;α) = -Σ_{n≥1} a_n c^{2n} J_n(u)`.
* stmt_3: `g6:lem:series`(b): `J_n(u) = ½∫_0^1 [y(1-y)]^{n-1} E_n(u,y) dy`, `|u| < 1`.
* stmt_4: `g6:lem:series`(c), first half: `M(·;α)` is `C^∞` on `(-δ,δ)` for some
  `0 < δ ≤ 1/12`, and `M_k(α) = -Σ a_n c^{2n} J_n^{(k)}(0)/k!` for every `k`.
* stmt_5: `g6:lem:series`(c), second half: `J_n` is `C^k` at `0` and
  `J_n^{(k)}(0)/k! = ½∫ [y(1-y)]^{n-1} ∂_u^k E_n(0,y)/k! dy`, every `k`.
* stmt_6: `g6:lem:En` for `k ≤ 3`: `∂_u^k E_n(0,y)/k! = Σ_{i+m=k} e^{(k)}_{i,m}(y) a^i b^m`.
* stmt_7: `g6:lem:gamma` for `i + m ≤ 3`: `∫_0^1 y^{p-1}(1-y)^{q-1} log^i y log^m(1-y) dy
  = B(p,q) D^♯_{i,m}(p,q)` at `𝔷2 = ζ(2)`, `𝔷3 = ζ(3)` (integers `p, q ≥ 1`).
* stmt_8: `g6:prop:formula` for `k ≤ 3`: `M_k(α) = Σ_{n≥1} P_{k,n}(ζ(2),ζ(3)) c^{2n}`.
* stmt_9: `g6:lem:first`(a): `P_{0,n} = -1/(n(2n-1))`, `P_{1,n} = D_n/(n(2n-1)) - 1/(2n²)`.
* stmt_10: `g6:lem:first`(b): `P_{2,n} = r_n + 𝔷2/(4n²-1)`, `r_n ∈ ℚ`.
* stmt_11: `g6:lem:first`(c): `P_{3,n} = r'_n + q'_n 𝔷2 - 𝔷3/(4n²-1)`, `r'_n, q'_n ∈ ℚ`.
* stmt_12: `thm:flat` for `k = 0,1,2,3`: rationals `r_n, r'_n, q'_n` (independent of `α`, given
  by the formal coefficients) such that for every `0<α<π`, `M(·;α)` is smooth at `0`,
  `M_0 = Σ -c^{2n}/(n(2n-1))`, `M_1 = Σ (D_n/(n(2n-1)) - 1/(2n²)) c^{2n}`,
  `M_2 = Σ (r_n + ζ(2)/(4n²-1)) c^{2n}`, `M_3 = Σ (r'_n + q'_n ζ(2) - ζ(3)/(4n²-1)) c^{2n}`,
  `M_2 - ζ(2)Φ(c) = Σ r_n c^{2n}` and, with `Ψ = Σ q'_n c^{2n}`,
  `M_3 + ζ(3)Φ(c) - ζ(2)Ψ = Σ r'_n c^{2n}`.

Numerical check before freezing (mpmath, 30 digits, script mirroring these definitions with
exact rationals): `P_{0,n}`, `P_{1,n}` closed forms and the `𝔷2`/`𝔷3` coefficients
`±1/(4n²-1)` hold exactly for `n ≤ 7`; `r_1 = -7/12`, `r_2 = -103/864`, `r_3 = -3139/54000`;
`P_{3,1} = 35/72 - 𝔷2/36 - 𝔷3/3` (as in the supplement). `-a_n J_n^{(k)}(0)/k!` computed from
the `g`-integral by a Cauchy integral in `u` agrees with `P_{k,n}(ζ(2),ζ(3))` for `n = 1,2,3`,
`k ≤ 3` to ~1e-21; `M_k(α)` from the defining integral (Cauchy integral, radius 0.2) agrees with
`Σ_{n<90} μ_{k,n} c^{2n}` to 15 digits at `α = π/2` (`M_2 = -0.0221619110872694`) and
`α = 2π/3` (`M_2 = -0.00956789565295638`), `k ≤ 3`.

## NOT COVERED

* `k ≥ 4`, and the general claim "for every `k`, `μ_{k,n}` is a polynomial with rational
  coefficients in `ζ(2),…,ζ(k)` of weight at most `k`" (the general `D^♯_{i,m}` via
  `exp(Δ_⋆)`, the weight grading `g6:eq:iso`, and `E_{n,k}` for general `k`).
* `Ψ ∈ ℚ[[c²]]` is recorded only through its rational coefficients `q'_n` and the convergence
  of `Σ q'_n c^{2n}`; the uniqueness remark (independence of `1, π², ζ(3)`) is not stated.
* "`δ` depends only on `c`" and the absolute convergence in `g6:lem:series`(c) (only `HasSum`).
* The closed form of `M_0(α)` (main.tex line 1474) and the bound
  `|coeff P_{k,n}| ≤ C_k n^k (1+log(2n+k))^k`.
-/
-- BEGIN STATEMENTS
namespace V05

open Real MeasureTheory Set Filter

/-- `ζ(2) = π²/6`. -/
noncomputable def z2 : ℝ := π ^ 2 / 6

/-- `ζ(3)`, the real part of Mathlib's `riemannZeta 3`. -/
noncomputable def z3 : ℝ := (riemannZeta 3).re

/-- `M_k(α) = M^{(k)}(0;α)/k!`, the Taylor coefficient at `u = 0` of the corner integral
`M(u;α)` (`D04.Mint`). -/
noncomputable def Mk (k : ℕ) (α : ℝ) : ℝ :=
  iteratedDeriv k (fun u => D04.Mint u α) 0 / (k.factorial : ℝ)

/-- `Φ(c) = 1/2 - (1-c²)/(2c) artanh c`. -/
noncomputable def Phi (c : ℝ) : ℝ := 1 / 2 - (1 - c ^ 2) / (2 * c) * Real.artanh c

/-- `a_n = binom(2n,n)/(2n-1)`. -/
noncomputable def an (n : ℕ) : ℝ := (Nat.choose (2 * n) n : ℝ) / (2 * (n : ℝ) - 1)

/-- `J_n(u) = ∫_0^1 (1+g)^{-1-u} h^{n-1} (1+h)^{1-2n} dg`, `h = g^{1-u}`. -/
noncomputable def J (n : ℕ) (u : ℝ) : ℝ :=
  ∫ g in (0 : ℝ)..1,
    (1 + g) ^ (-1 - u) * (g ^ (1 - u)) ^ (n - 1) * (1 + g ^ (1 - u)) ^ (1 - 2 * (n : ℝ))

/-- `f(u,y) = 1 - y + y e^{uℓ}`, `ℓ = log(1-y) - log y`. -/
noncomputable def f (u y : ℝ) : ℝ := 1 - y + y * Real.exp (u * (Real.log (1 - y) - Real.log y))

/-- `E_n(u,y) = e^{uA} f(u,y)^{1-2n}`, `A = -(n-1) log y + n log(1-y)`. -/
noncomputable def E (n : ℕ) (u y : ℝ) : ℝ :=
  Real.exp (u * (-((n : ℝ) - 1) * Real.log y + (n : ℝ) * Real.log (1 - y)))
    * f u y ^ (1 - 2 * (n : ℝ))

/-- Generalized harmonic number `H^{(r)}_N = Σ_{l=1}^N l^{-r}`. -/
noncomputable def Hr (r N : ℕ) : ℝ := ∑ l ∈ Finset.range N, 1 / ((l : ℝ) + 1) ^ r

/-- `D^♯_{i,m}(p,q)` for `i + m ≤ 3` (degree `≤ 3` part of `exp(Δ_⋆)`), with the formal zeta
values `𝔷2 = x2`, `𝔷3 = x3`; `0` (unused) for `i + m > 3`. -/
noncomputable def Dsharp (i m p q : ℕ) (x2 x3 : ℝ) : ℝ :=
  let A := Hr 1 (p - 1) - Hr 1 (p + q - 1)
  let B := Hr 1 (q - 1) - Hr 1 (p + q - 1)
  let P := Hr 2 (p - 1) - Hr 2 (p + q - 1)
  let Q := Hr 2 (q - 1) - Hr 2 (p + q - 1)
  let R := Hr 2 (p + q - 1) - x2
  let S := Hr 3 (p - 1) - Hr 3 (p + q - 1)
  let T := Hr 3 (q - 1) - Hr 3 (p + q - 1)
  let W := x3 - Hr 3 (p + q - 1)
  match i, m with
  | 0, 0 => 1
  | 1, 0 => A
  | 0, 1 => B
  | 2, 0 => A ^ 2 - P
  | 0, 2 => B ^ 2 - Q
  | 1, 1 => A * B + R
  | 3, 0 => 2 * S - 3 * A * P + A ^ 3
  | 0, 3 => 2 * T - 3 * B * Q + B ^ 3
  | 2, 1 => 2 * W + 2 * A * R - B * P + A ^ 2 * B
  | 1, 2 => 2 * W + 2 * B * R - A * Q + A * B ^ 2
  | _, _ => 0

/-- `α₁ = -(n-1) + (2n-1) y`. -/
noncomputable def al1 (n : ℕ) : Polynomial ℝ :=
  Polynomial.C (-((n : ℝ) - 1)) + Polynomial.C (2 * (n : ℝ) - 1) * Polynomial.X

/-- `β₁ = n - (2n-1) y`. -/
noncomputable def be1 (n : ℕ) : Polynomial ℝ :=
  Polynomial.C (n : ℝ) - Polynomial.C (2 * (n : ℝ) - 1) * Polynomial.X

/-- `γ = -½(2n-1) y(1-y)`. -/
noncomputable def ga (n : ℕ) : Polynomial ℝ :=
  Polynomial.C (-(2 * (n : ℝ) - 1) / 2) * Polynomial.X * (1 - Polynomial.X)

/-- `γ₃ = -⅙(2n-1) y(1-y)(1-2y)`. -/
noncomputable def ga3 (n : ℕ) : Polynomial ℝ :=
  Polynomial.C (-(2 * (n : ℝ) - 1) / 6) * Polynomial.X * (1 - Polynomial.X)
    * (1 - 2 * Polynomial.X)

/-- The coefficient of `a^i b^m` in `E_{n,k}` (`k ≤ 3`), a polynomial in `y` (its `y^j`
coefficient is `e^{(k)}_{j,i,m}(n)`), from `E_{n,0} = 1`, `E_{n,1} = X`,
`E_{n,2} = ½X² + γℓ²`, `E_{n,3} = ⅙X³ + γXℓ² + γ₃ℓ³`, `X = α₁a + β₁b`, `ℓ = b - a`. -/
noncomputable def ecoef (k n i m : ℕ) : Polynomial ℝ :=
  match k, i, m with
  | 0, 0, 0 => 1
  | 1, 1, 0 => al1 n
  | 1, 0, 1 => be1 n
  | 2, 2, 0 => Polynomial.C (1 / 2) * al1 n ^ 2 + ga n
  | 2, 0, 2 => Polynomial.C (1 / 2) * be1 n ^ 2 + ga n
  | 2, 1, 1 => al1 n * be1 n - 2 * ga n
  | 3, 3, 0 => Polynomial.C (1 / 6) * al1 n ^ 3 + ga n * al1 n - ga3 n
  | 3, 2, 1 => Polynomial.C (1 / 2) * al1 n ^ 2 * be1 n + ga n * (be1 n - 2 * al1 n) + 3 * ga3 n
  | 3, 1, 2 => Polynomial.C (1 / 2) * al1 n * be1 n ^ 2 + ga n * (al1 n - 2 * be1 n) - 3 * ga3 n
  | 3, 0, 3 => Polynomial.C (1 / 6) * be1 n ^ 3 + ga n * be1 n + ga3 n
  | _, _, _ => 0

/-- `ρ_j(n) = (n)_j/(2n)_j`. -/
noncomputable def rho (j n : ℕ) : ℝ :=
  ∏ t ∈ Finset.range j, ((n : ℝ) + t) / (2 * (n : ℝ) + t)

/-- The formal coefficient `P_{k,n}` (Definition `g6:def:formal`), `k ≤ 3`, as a function of the
formal zeta values `x2 = 𝔷2`, `x3 = 𝔷3`. -/
noncomputable def P (k n : ℕ) (x2 x3 : ℝ) : ℝ :=
  -(1 / ((n : ℝ) * (2 * (n : ℝ) - 1))) *
    ∑ j ∈ Finset.range (k + 1), ∑ i ∈ Finset.range (k + 1),
      (ecoef k n i (k - i)).coeff j * rho j n * Dsharp i (k - i) (n + j) n x2 x3

/-- `B(p,q) = (p-1)!(q-1)!/(p+q-1)!` for integers `p, q ≥ 1`. -/
noncomputable def Bet (p q : ℕ) : ℝ :=
  ((p - 1).factorial * (q - 1).factorial : ℝ) / ((p + q - 1).factorial : ℝ)

/-- The paper's `Φ(c) = Σ_{n≥1} c^{2n}/(4n²-1) = 1/2 - (1-c²)/(2c) artanh c`, `0 < |c| < 1`. -/
def stmt_1 : Prop :=
  ∀ c : ℝ, 0 < |c| → |c| < 1 →
    HasSum (fun n : ℕ => c ^ (2 * (n + 1)) / (4 * ((n : ℝ) + 1) ^ 2 - 1)) (Phi c)

/-- `g6:lem:series`(a): for `|u| < 1` the integrand of `M(u;α)` is integrable, all `J_n(u) > 0`,
and `M(u;α) = -Σ_{n≥1} a_n c^{2n} J_n(u)`. -/
def stmt_2 : Prop :=
  ∀ α : ℝ, 0 < α → α < π → ∀ u : ℝ, |u| < 1 →
    IntervalIntegrable (fun g : ℝ =>
      ((1 + g) ^ (-1 - u) * Real.sqrt (1 + g ^ (2 - 2 * u) - 2 * Real.cos α * g ^ (1 - u)) - 1)
        * g ^ (u - 1)) volume 0 1 ∧
    (∀ n : ℕ, 1 ≤ n → 0 < J n u) ∧
    HasSum (fun n : ℕ => an (n + 1) * D04.c α ^ (2 * (n + 1)) * J (n + 1) u) (-(D04.Mint u α))

/-- `g6:lem:series`(b): `J_n(u) = ½ ∫_0^1 [y(1-y)]^{n-1} E_n(u,y) dy` for `|u| < 1`. -/
def stmt_3 : Prop :=
  ∀ n : ℕ, 1 ≤ n → ∀ u : ℝ, |u| < 1 →
    J n u = 1 / 2 * ∫ y in (0 : ℝ)..1, (y * (1 - y)) ^ (n - 1) * E n u y

/-- `g6:lem:series`(c), first part: `M(·;α)` is `C^∞` on `(-δ,δ)` for some `0 < δ ≤ 1/12`, and
`M_k(α) = -Σ_{n≥1} a_n c^{2n} J_n^{(k)}(0)/k!` for every `k`. -/
def stmt_4 : Prop :=
  ∀ α : ℝ, 0 < α → α < π → ∃ δ : ℝ, 0 < δ ∧ δ ≤ 1 / 12 ∧
    (∀ k : ℕ, ContDiffOn ℝ k (fun u => D04.Mint u α) (Ioo (-δ) δ)) ∧
    ∀ k : ℕ, HasSum (fun n : ℕ => an (n + 1) * D04.c α ^ (2 * (n + 1)) *
      (iteratedDeriv k (J (n + 1)) 0 / (k.factorial : ℝ))) (-(Mk k α))

/-- `g6:lem:series`(c), second part: `J_n` is smooth at `0` and
`J_n^{(k)}(0)/k! = ½ ∫_0^1 [y(1-y)]^{n-1} ∂_u^k E_n(0,y)/k! dy`. -/
def stmt_5 : Prop :=
  ∀ n : ℕ, 1 ≤ n → ∀ k : ℕ, ContDiffAt ℝ k (J n) 0 ∧
    iteratedDeriv k (J n) 0 / (k.factorial : ℝ) =
      1 / 2 * ∫ y in (0 : ℝ)..1, (y * (1 - y)) ^ (n - 1) *
        (iteratedDeriv k (fun u => E n u y) 0 / (k.factorial : ℝ))

/-- `g6:lem:En` for `k ≤ 3`: `∂_u^k E_n(0,y)/k! = E_{n,k}(y, log y, log(1-y))`, written as
`Σ_{i+m=k} e^{(k)}_{i,m}(y) a^i b^m`. -/
def stmt_6 : Prop :=
  ∀ n : ℕ, 1 ≤ n → ∀ y : ℝ, 0 < y → y < 1 → ∀ k : ℕ, k ≤ 3 →
    iteratedDeriv k (fun u => E n u y) 0 / (k.factorial : ℝ) =
      ∑ i ∈ Finset.range (k + 1),
        (ecoef k n i (k - i)).eval y * Real.log y ^ i * Real.log (1 - y) ^ (k - i)

/-- `g6:lem:gamma` for `i + m ≤ 3`: for integers `p, q ≥ 1`,
`∫_0^1 y^{p-1}(1-y)^{q-1} log^i y log^m(1-y) dy = B(p,q) D^♯_{i,m}(p,q)` at `𝔷w = ζ(w)`. -/
def stmt_7 : Prop :=
  ∀ p q : ℕ, 1 ≤ p → 1 ≤ q → ∀ i m : ℕ, i + m ≤ 3 →
    IntervalIntegrable (fun y : ℝ =>
      y ^ (p - 1) * (1 - y) ^ (q - 1) * Real.log y ^ i * Real.log (1 - y) ^ m) volume 0 1 ∧
    (∫ y in (0 : ℝ)..1, y ^ (p - 1) * (1 - y) ^ (q - 1) * Real.log y ^ i * Real.log (1 - y) ^ m)
      = Bet p q * Dsharp i m p q z2 z3

/-- `g6:prop:formula` for `k ≤ 3`: `M_k(α) = Σ_{n≥1} P_{k,n}(ζ(2),ζ(3)) c^{2n}`. -/
def stmt_8 : Prop :=
  ∀ α : ℝ, 0 < α → α < π → ∀ k : ℕ, k ≤ 3 →
    HasSum (fun n : ℕ => P k (n + 1) z2 z3 * D04.c α ^ (2 * (n + 1))) (Mk k α)

/-- `g6:lem:first`(a): `P_{0,n} = -1/(n(2n-1))` and `P_{1,n} = D_n/(n(2n-1)) - 1/(2n²)`,
`D_n = H_{2n-1} - H_{n-1}`. -/
def stmt_9 : Prop :=
  ∀ n : ℕ, 1 ≤ n → ∀ x2 x3 : ℝ,
    P 0 n x2 x3 = -(1 / ((n : ℝ) * (2 * (n : ℝ) - 1))) ∧
    P 1 n x2 x3 = (Hr 1 (2 * n - 1) - Hr 1 (n - 1)) / ((n : ℝ) * (2 * (n : ℝ) - 1))
      - 1 / (2 * (n : ℝ) ^ 2)

/-- `g6:lem:first`(b): `P_{2,n} = r_n + 𝔷2/(4n²-1)` with `r_n ∈ ℚ` (identity in `𝔷2, 𝔷3`). -/
def stmt_10 : Prop :=
  ∃ r : ℕ → ℚ, ∀ n : ℕ, 1 ≤ n → ∀ x2 x3 : ℝ,
    P 2 n x2 x3 = (r n : ℝ) + x2 / (4 * (n : ℝ) ^ 2 - 1)

/-- `g6:lem:first`(c): `P_{3,n} = r'_n + q'_n 𝔷2 - 𝔷3/(4n²-1)` with `r'_n, q'_n ∈ ℚ`. -/
def stmt_11 : Prop :=
  ∃ r' q' : ℕ → ℚ, ∀ n : ℕ, 1 ≤ n → ∀ x2 x3 : ℝ,
    P 3 n x2 x3 = (r' n : ℝ) + (q' n : ℝ) * x2 - x3 / (4 * (n : ℝ) ^ 2 - 1)

/-- `thm:flat` for `k = 0, 1, 2, 3`: with rationals `r_n, r'_n, q'_n` independent of the angle
(the evaluated formal coefficients), for every `0 < α < π`: `M(·;α)` is smooth at `0`,
`M_k(α) = Σ_{n≥1} μ_{k,n} c^{2n}` with `μ_{0,n} = -1/(n(2n-1))`,
`μ_{1,n} = D_n/(n(2n-1)) - 1/(2n²)`, `μ_{2,n} = r_n + ζ(2)/(4n²-1)`,
`μ_{3,n} = r'_n + q'_n ζ(2) - ζ(3)/(4n²-1)`; `M_2 - ζ(2)Φ = Σ r_n c^{2n}`; and with
`Ψ = Σ q'_n c^{2n}`, `M_3 + ζ(3)Φ - ζ(2)Ψ = Σ r'_n c^{2n}`. -/
def stmt_12 : Prop :=
  ∃ r r' q' : ℕ → ℚ,
    (∀ n : ℕ, 1 ≤ n → (r n : ℝ) = P 2 n 0 0 ∧ (r' n : ℝ) = P 3 n 0 0 ∧
      (q' n : ℝ) = P 3 n 1 0 - P 3 n 0 0) ∧
    ∀ α : ℝ, 0 < α → α < π →
      (∀ k : ℕ, ContDiffAt ℝ k (fun u => D04.Mint u α) 0) ∧
      HasSum (fun n : ℕ => -(1 / (((n : ℝ) + 1) * (2 * ((n : ℝ) + 1) - 1)))
        * D04.c α ^ (2 * (n + 1))) (Mk 0 α) ∧
      HasSum (fun n : ℕ => ((Hr 1 (2 * (n + 1) - 1) - Hr 1 n) / (((n : ℝ) + 1) * (2 * ((n : ℝ) + 1) - 1))
        - 1 / (2 * ((n : ℝ) + 1) ^ 2)) * D04.c α ^ (2 * (n + 1))) (Mk 1 α) ∧
      HasSum (fun n : ℕ => ((r (n + 1) : ℝ) + z2 / (4 * ((n : ℝ) + 1) ^ 2 - 1))
        * D04.c α ^ (2 * (n + 1))) (Mk 2 α) ∧
      HasSum (fun n : ℕ => ((r' (n + 1) : ℝ) + (q' (n + 1) : ℝ) * z2
        - z3 / (4 * ((n : ℝ) + 1) ^ 2 - 1)) * D04.c α ^ (2 * (n + 1))) (Mk 3 α) ∧
      HasSum (fun n : ℕ => (r (n + 1) : ℝ) * D04.c α ^ (2 * (n + 1)))
        (Mk 2 α - z2 * Phi (D04.c α)) ∧
      ∃ Ψ : ℝ, HasSum (fun n : ℕ => (q' (n + 1) : ℝ) * D04.c α ^ (2 * (n + 1))) Ψ ∧
        HasSum (fun n : ℕ => (r' (n + 1) : ℝ) * D04.c α ^ (2 * (n + 1)))
          (Mk 3 α + z3 * Phi (D04.c α) - z2 * Ψ)

end V05
-- END STATEMENTS

namespace V05

open Real MeasureTheory Set Filter

noncomputable def HrQ (r N : ℕ) : ℚ := ∑ l ∈ Finset.range N, 1 / ((l : ℚ) + 1) ^ r

lemma HrQ_cast (r N : ℕ) : ((HrQ r N : ℚ) : ℝ) = Hr r N := by
  simp [HrQ, Hr]

noncomputable def DsharpQ (i m p q : ℕ) (x2 x3 : ℚ) : ℚ :=
  let A := HrQ 1 (p - 1) - HrQ 1 (p + q - 1)
  let B := HrQ 1 (q - 1) - HrQ 1 (p + q - 1)
  let P := HrQ 2 (p - 1) - HrQ 2 (p + q - 1)
  let Q := HrQ 2 (q - 1) - HrQ 2 (p + q - 1)
  let R := HrQ 2 (p + q - 1) - x2
  let S := HrQ 3 (p - 1) - HrQ 3 (p + q - 1)
  let T := HrQ 3 (q - 1) - HrQ 3 (p + q - 1)
  let W := x3 - HrQ 3 (p + q - 1)
  match i, m with
  | 0, 0 => 1
  | 1, 0 => A
  | 0, 1 => B
  | 2, 0 => A ^ 2 - P
  | 0, 2 => B ^ 2 - Q
  | 1, 1 => A * B + R
  | 3, 0 => 2 * S - 3 * A * P + A ^ 3
  | 0, 3 => 2 * T - 3 * B * Q + B ^ 3
  | 2, 1 => 2 * W + 2 * A * R - B * P + A ^ 2 * B
  | 1, 2 => 2 * W + 2 * B * R - A * Q + A * B ^ 2
  | _, _ => 0

lemma DsharpQ_cast (i m p q : ℕ) (x2 x3 : ℚ) :
    ((DsharpQ i m p q x2 x3 : ℚ) : ℝ) = Dsharp i m p q x2 x3 := by
  rcases i with _ | _ | _ | _ | i <;> rcases m with _ | _ | _ | _ | m <;>
    simp [DsharpQ, Dsharp, HrQ_cast]

noncomputable def al1Q (n : ℕ) : Polynomial ℚ :=
  Polynomial.C (-((n : ℚ) - 1)) + Polynomial.C (2 * (n : ℚ) - 1) * Polynomial.X
noncomputable def be1Q (n : ℕ) : Polynomial ℚ :=
  Polynomial.C (n : ℚ) - Polynomial.C (2 * (n : ℚ) - 1) * Polynomial.X
noncomputable def gaQ (n : ℕ) : Polynomial ℚ :=
  Polynomial.C (-(2 * (n : ℚ) - 1) / 2) * Polynomial.X * (1 - Polynomial.X)
noncomputable def ga3Q (n : ℕ) : Polynomial ℚ :=
  Polynomial.C (-(2 * (n : ℚ) - 1) / 6) * Polynomial.X * (1 - Polynomial.X)
    * (1 - 2 * Polynomial.X)

noncomputable def ecoefQ (k n i m : ℕ) : Polynomial ℚ :=
  match k, i, m with
  | 0, 0, 0 => 1
  | 1, 1, 0 => al1Q n
  | 1, 0, 1 => be1Q n
  | 2, 2, 0 => Polynomial.C (1 / 2) * al1Q n ^ 2 + gaQ n
  | 2, 0, 2 => Polynomial.C (1 / 2) * be1Q n ^ 2 + gaQ n
  | 2, 1, 1 => al1Q n * be1Q n - 2 * gaQ n
  | 3, 3, 0 => Polynomial.C (1 / 6) * al1Q n ^ 3 + gaQ n * al1Q n - ga3Q n
  | 3, 2, 1 => Polynomial.C (1 / 2) * al1Q n ^ 2 * be1Q n + gaQ n * (be1Q n - 2 * al1Q n) + 3 * ga3Q n
  | 3, 1, 2 => Polynomial.C (1 / 2) * al1Q n * be1Q n ^ 2 + gaQ n * (al1Q n - 2 * be1Q n) - 3 * ga3Q n
  | 3, 0, 3 => Polynomial.C (1 / 6) * be1Q n ^ 3 + gaQ n * be1Q n + ga3Q n
  | _, _, _ => 0

lemma ecoef_map (k n i m : ℕ) :
    ecoef k n i m = (ecoefQ k n i m).map (algebraMap ℚ ℝ) := by
  rcases k with _ | _ | _ | _ | k <;> rcases i with _ | _ | _ | _ | i <;>
    rcases m with _ | _ | _ | _ | m <;>
    simp [ecoef, ecoefQ, al1, al1Q, be1, be1Q, ga, gaQ, ga3, ga3Q, Polynomial.map_sub]

noncomputable def rhoQ (j n : ℕ) : ℚ :=
  ∏ t ∈ Finset.range j, ((n : ℚ) + t) / (2 * (n : ℚ) + t)

noncomputable def PQ (k n : ℕ) (x2 x3 : ℚ) : ℚ :=
  -(1 / ((n : ℚ) * (2 * (n : ℚ) - 1))) *
    ∑ j ∈ Finset.range (k + 1), ∑ i ∈ Finset.range (k + 1),
      (ecoefQ k n i (k - i)).coeff j * rhoQ j n * DsharpQ i (k - i) (n + j) n x2 x3

lemma PQ_cast (k n : ℕ) (x2 x3 : ℚ) : ((PQ k n x2 x3 : ℚ) : ℝ) = P k n x2 x3 := by
  simp [PQ, P, rhoQ, rho, ecoef_map, Polynomial.coeff_map, DsharpQ_cast]

/-- the `𝔷2`-part of `D^♯`. -/
noncomputable def D2 (i m p q : ℕ) : ℝ :=
  match i, m with
  | 1, 1 => -1
  | 2, 1 => -2 * (Hr 1 (p - 1) - Hr 1 (p + q - 1))
  | 1, 2 => -2 * (Hr 1 (q - 1) - Hr 1 (p + q - 1))
  | _, _ => 0

/-- the `𝔷3`-part of `D^♯`. -/
noncomputable def D3 (i m : ℕ) : ℝ :=
  match i, m with
  | 2, 1 => 2
  | 1, 2 => 2
  | _, _ => 0

lemma Dsharp_lin (i m p q : ℕ) (x2 x3 : ℝ) :
    Dsharp i m p q x2 x3 = Dsharp i m p q 0 0 + x2 * D2 i m p q + x3 * D3 i m := by
  rcases i with _ | _ | _ | _ | i <;> rcases m with _ | _ | _ | _ | m <;>
    simp [Dsharp, D2, D3] <;> ring

noncomputable def L2 (k n : ℕ) : ℝ :=
  -(1 / ((n : ℝ) * (2 * (n : ℝ) - 1))) *
    ∑ j ∈ Finset.range (k + 1), ∑ i ∈ Finset.range (k + 1),
      (ecoef k n i (k - i)).coeff j * rho j n * D2 i (k - i) (n + j) n

noncomputable def L3 (k n : ℕ) : ℝ :=
  -(1 / ((n : ℝ) * (2 * (n : ℝ) - 1))) *
    ∑ j ∈ Finset.range (k + 1), ∑ i ∈ Finset.range (k + 1),
      (ecoef k n i (k - i)).coeff j * rho j n * D3 i (k - i)

lemma P_lin (k n : ℕ) (x2 x3 : ℝ) :
    P k n x2 x3 = P k n 0 0 + x2 * L2 k n + x3 * L3 k n := by
  simp only [P, L2, L3, Finset.mul_sum, ← Finset.sum_add_distrib]
  refine Finset.sum_congr rfl fun j _ => Finset.sum_congr rfl fun i _ => ?_
  rw [Dsharp_lin _ _ _ _ x2 x3]; ring

lemma e211_eq (n : ℕ) : ecoef 2 n 1 1 = Polynomial.C (-((n : ℝ) - 1) * n)
    + Polynomial.C (2 * (n : ℝ) * (2 * n - 1)) * Polynomial.X
    + Polynomial.C (-2 * (n : ℝ) * (2 * n - 1)) * Polynomial.X ^ 2 := by
  apply Polynomial.funext; intro y
  simp [ecoef, al1, be1, ga]; ring

lemma e3sum_eq (n : ℕ) : ecoef 3 n 2 1 + ecoef 3 n 1 2 = Polynomial.C (1 / 2) * ecoef 2 n 1 1 := by
  apply Polynomial.funext; intro y
  simp [ecoef, al1, be1, ga, ga3]; ring

lemma coeff3 (a0 a1 a2 a3 : ℝ) (j : ℕ) :
    (Polynomial.C a0 + Polynomial.C a1 * Polynomial.X + Polynomial.C a2 * Polynomial.X ^ 2
      + Polynomial.C a3 * Polynomial.X ^ 3).coeff j =
      if j = 0 then a0 else if j = 1 then a1 else if j = 2 then a2 else if j = 3 then a3 else 0 := by
  simp only [Polynomial.coeff_add, Polynomial.coeff_C_mul, Polynomial.coeff_X_pow,
    Polynomial.coeff_X, Polynomial.coeff_C]
  rcases j with _ | _ | _ | _ | j <;> simp

lemma e211_eq3 (n : ℕ) : ecoef 2 n 1 1 = Polynomial.C (-((n : ℝ) - 1) * n)
    + Polynomial.C (2 * (n : ℝ) * (2 * n - 1)) * Polynomial.X
    + Polynomial.C (-2 * (n : ℝ) * (2 * n - 1)) * Polynomial.X ^ 2
    + Polynomial.C 0 * Polynomial.X ^ 3 := by
  rw [e211_eq]; simp

lemma L2_two {n : ℕ} (hn : 1 ≤ n) : L2 2 n = 1 / (4 * (n : ℝ) ^ 2 - 1) := by
  have hn' : (1 : ℝ) ≤ n := by exact_mod_cast hn
  simp only [L2, Finset.sum_range_succ, Finset.sum_range_zero, D2, rho, Finset.prod_range_succ,
    Finset.prod_range_zero]
  norm_num
  simp only [e211_eq3, coeff3]
  norm_num
  have h1 : (n : ℝ) ≠ 0 := by positivity
  have h2 : (2 * (n : ℝ) - 1) ≠ 0 := by linarith
  have h3 : (2 * (n : ℝ) + 1) ≠ 0 := by linarith
  have h4 : (4 * (n : ℝ) ^ 2 - 1) ≠ 0 := by nlinarith
  have h5 : (-1 + (n : ℝ) ^ 2 * 4) ≠ 0 := by nlinarith
  field_simp
  ring_nf
  field_simp
  ring

lemma L3_two (n : ℕ) : L3 2 n = 0 := by
  simp [L3, Finset.sum_range_succ, D3]

lemma L3_three {n : ℕ} (hn : 1 ≤ n) : L3 3 n = -(1 / (4 * (n : ℝ) ^ 2 - 1)) := by
  have hn' : (1 : ℝ) ≤ n := by exact_mod_cast hn
  have hc : ∀ j, (ecoef 3 n 2 1).coeff j = 1 / 2 * (ecoef 2 n 1 1).coeff j - (ecoef 3 n 1 2).coeff j := by
    intro j
    have := congrArg (fun p : Polynomial ℝ => p.coeff j) (e3sum_eq n)
    simp only [Polynomial.coeff_add, Polynomial.coeff_C_mul] at this
    linarith
  simp only [L3, Finset.sum_range_succ, Finset.sum_range_zero, D3, rho, Finset.prod_range_succ,
    Finset.prod_range_zero]
  norm_num
  simp only [hc, e211_eq3, coeff3]
  norm_num
  have h1 : (n : ℝ) ≠ 0 := by positivity
  have h2 : (2 * (n : ℝ) - 1) ≠ 0 := by linarith
  have h3 : (2 * (n : ℝ) + 1) ≠ 0 := by linarith
  have h4 : (4 * (n : ℝ) ^ 2 - 1) ≠ 0 := by nlinarith
  have h5 : (-1 + (n : ℝ) ^ 2 * 4) ≠ 0 := by nlinarith
  field_simp
  ring_nf
  field_simp

lemma Hr_succ (r N : ℕ) : Hr r (N + 1) = Hr r N + 1 / ((N : ℝ) + 1) ^ r := by
  simp [Hr, Finset.sum_range_succ]

lemma P0_eq (n : ℕ) (x2 x3 : ℝ) : P 0 n x2 x3 = -(1 / ((n : ℝ) * (2 * (n : ℝ) - 1))) := by
  simp [P, ecoef, rho, Dsharp]

lemma P1_eq {n : ℕ} (hn : 1 ≤ n) (x2 x3 : ℝ) : P 1 n x2 x3 =
    (Hr 1 (2 * n - 1) - Hr 1 (n - 1)) / ((n : ℝ) * (2 * (n : ℝ) - 1)) - 1 / (2 * (n : ℝ) ^ 2) := by
  obtain ⟨m, rfl⟩ : ∃ m, n = m + 1 := ⟨n - 1, by omega⟩
  simp only [P, Finset.sum_range_succ, Finset.sum_range_zero, rho, Finset.prod_range_succ,
    Finset.prod_range_zero]
  norm_num
  simp only [ecoef, al1, be1, Polynomial.coeff_add, Polynomial.coeff_sub, Polynomial.coeff_C_mul,
    Polynomial.coeff_C, Polynomial.coeff_X, Dsharp]
  norm_num
  rw [show m + 1 + m = 2 * m + 1 by ring, show m + 1 + 1 + m = 2 * m + 1 + 1 by ring,
    show 2 * (m + 1) - 1 = 2 * m + 1 by omega, Hr_succ 1 (2 * m + 1), Hr_succ 1 m]
  have h1 : ((m : ℝ) + 1) ≠ 0 := by positivity
  have h2 : (2 * ((m : ℝ) + 1) - 1) ≠ 0 := by have : (0:ℝ) ≤ m := by positivity
                                              linarith
  push_cast
  field_simp
  ring

lemma PQ_cast00 (k n : ℕ) : ((PQ k n 0 0 : ℚ) : ℝ) = P k n 0 0 := by
  have := PQ_cast k n 0 0; push_cast at this; exact this

lemma PQ_cast10 (k n : ℕ) : ((PQ k n 1 0 : ℚ) : ℝ) = P k n 1 0 := by
  have := PQ_cast k n 1 0; push_cast at this; exact this

lemma P3_x2 (n : ℕ) : P 3 n 1 0 - P 3 n 0 0 = L2 3 n := by
  rw [P_lin 3 n 1 0]; ring

end V05

theorem V05.proof_9 : V05.stmt_9 := fun _ hn x2 x3 => ⟨V05.P0_eq _ x2 x3, V05.P1_eq hn x2 x3⟩

theorem V05.proof_10 : V05.stmt_10 := by
  refine ⟨fun n => V05.PQ 2 n 0 0, fun n hn x2 x3 => ?_⟩
  rw [V05.P_lin, V05.L2_two hn, V05.L3_two, V05.PQ_cast00]; ring

theorem V05.proof_11 : V05.stmt_11 := by
  refine ⟨fun n => V05.PQ 3 n 0 0, fun n => V05.PQ 3 n 1 0 - V05.PQ 3 n 0 0, fun n hn x2 x3 => ?_⟩
  push_cast
  rw [V05.P_lin, V05.L3_three hn, V05.PQ_cast00, V05.PQ_cast10, V05.P3_x2]; ring

namespace V05
open Real

lemma artanh_hasSum {c : ℝ} (hc : |c| < 1) :
    HasSum (fun k : ℕ => c ^ (2 * k + 1) / (2 * (k : ℝ) + 1)) (Real.artanh c) := by
  have h := (Real.hasSum_log_sub_log_of_abs_lt_one hc).mul_left (1 / 2)
  have h1 : 0 < 1 + c := by linarith [neg_abs_le c]
  have h2 : 0 < 1 - c := by linarith [le_abs_self c]
  rw [Real.artanh_eq_half_log ⟨by linarith, by linarith⟩, Real.log_div h1.ne' h2.ne']
  convert h using 1
  funext k; field_simp

end V05

theorem V05.proof_1 : V05.stmt_1 := by
  intro c hc0 hc1
  have hc : c ≠ 0 := abs_pos.mp hc0
  have h := V05.artanh_hasSum hc1
  have hA := h.mul_left c
  have hB := ((hasSum_nat_add_iff' 1).mpr h).mul_left (1 / c)
  have hC := (hA.sub hB).mul_left (1 / 2)
  simp only [Finset.range_one, Finset.sum_singleton] at hC
  unfold V05.Phi
  convert hC using 1
  · funext n
    have h4 : (4 * ((n : ℝ) + 1) ^ 2 - 1) ≠ 0 := by
      have : (0 : ℝ) ≤ n := by positivity
      nlinarith
    have h5 : (2 * (n : ℝ) + 1) ≠ 0 := by positivity
    have h6 : (2 * ((n + 1 : ℕ) : ℝ) + 1) ≠ 0 := by positivity
    push_cast at h6 ⊢
    field_simp
    ring
  · norm_num
    field_simp
    ring

namespace V05
open Real

noncomputable def Ff (y l u : ℝ) : ℝ := 1 - y + y * Real.exp (u * l)

lemma Ff_pos {y : ℝ} (hy0 : 0 < y) (hy1 : y < 1) (l u : ℝ) : 0 < Ff y l u := by
  unfold Ff; have := Real.exp_pos (u * l); nlinarith

lemma hasDerivAt_Ff (y l u : ℝ) :
    HasDerivAt (Ff y l) (y * l * Real.exp (u * l)) u := by
  have h := (((hasDerivAt_id' u).mul_const l).exp).const_mul y |>.const_add (1 - y)
  unfold Ff
  convert h using 1; simp; ring

noncomputable def G1 (y l A s u : ℝ) : ℝ := A + s * y * l * Real.exp (u * l) / Ff y l u
noncomputable def G2 (y l s u : ℝ) : ℝ :=
  s * y * (1 - y) * l ^ 2 * Real.exp (u * l) / Ff y l u ^ 2
noncomputable def G3 (y l s u : ℝ) : ℝ :=
  s * y * (1 - y) * l ^ 3 * Real.exp (u * l) * (Ff y l u - 2 * y * Real.exp (u * l)) / Ff y l u ^ 3

noncomputable def EE (y l A s u : ℝ) : ℝ := Real.exp (u * A) * Ff y l u ^ s

section
variable {y : ℝ} (hy0 : 0 < y) (hy1 : y < 1) (l A s : ℝ)
include hy0 hy1

lemma hasDerivAt_EE (u : ℝ) :
    HasDerivAt (EE y l A s) (G1 y l A s u * EE y l A s u) u := by
  have hF := Ff_pos hy0 hy1 l u
  have h := (((hasDerivAt_id' u).mul_const A).exp).mul ((hasDerivAt_Ff y l u).rpow_const
    (p := s) (Or.inl hF.ne'))
  refine h.congr_deriv ?_
  unfold G1 EE
  rw [Real.rpow_sub_one hF.ne']
  field_simp

lemma hasDerivAt_G1 (u : ℝ) : HasDerivAt (G1 y l A s) (G2 y l s u) u := by
  have hF := Ff_pos hy0 hy1 l u
  have h := ((((hasDerivAt_id' u).mul_const l).exp).const_mul (s * y * l)).div
    (hasDerivAt_Ff y l u) hF.ne' |>.const_add A
  refine h.congr_deriv ?_
  unfold G2 Ff
  field_simp
  ring

lemma hasDerivAt_G2 (u : ℝ) : HasDerivAt (G2 y l s) (G3 y l s u) u := by
  have hF := Ff_pos hy0 hy1 l u
  have h := ((((hasDerivAt_id' u).mul_const l).exp).const_mul (s * y * (1 - y) * l ^ 2)).div
    ((hasDerivAt_Ff y l u).pow 2) (pow_ne_zero 2 hF.ne')
  refine h.congr_deriv ?_
  unfold G3
  simp only [Pi.pow_apply]
  field_simp
  ring

lemma deriv1_EE : deriv (EE y l A s) = fun u => G1 y l A s u * EE y l A s u := by
  funext u; exact (hasDerivAt_EE hy0 hy1 l A s u).deriv

lemma deriv2_EE : deriv (fun u => G1 y l A s u * EE y l A s u) =
    fun u => (G2 y l s u + G1 y l A s u ^ 2) * EE y l A s u := by
  funext u
  refine ((hasDerivAt_G1 hy0 hy1 l A s u).mul (hasDerivAt_EE hy0 hy1 l A s u)).deriv.trans ?_
  ring

lemma deriv3_EE : deriv (fun u => (G2 y l s u + G1 y l A s u ^ 2) * EE y l A s u) =
    fun u => (G3 y l s u + 3 * G1 y l A s u * G2 y l s u + G1 y l A s u ^ 3) * EE y l A s u := by
  funext u
  refine (((hasDerivAt_G2 hy0 hy1 l s u).add ((hasDerivAt_G1 hy0 hy1 l A s u).pow 2)).mul
    (hasDerivAt_EE hy0 hy1 l A s u)).deriv.trans ?_
  simp; ring

end

lemma EE_zero (y l A s : ℝ) : EE y l A s 0 = 1 := by simp [EE, Ff]

end V05

theorem V05.proof_6 : V05.stmt_6 := by
  intro n hn y hy0 hy1 k hk
  have e : (fun u => V05.E n u y) = V05.EE y (Real.log (1 - y) - Real.log y)
      (-((n : ℝ) - 1) * Real.log y + (n : ℝ) * Real.log (1 - y)) (1 - 2 * (n : ℝ)) := by
    funext u; rfl
  rw [e]
  have hF0 : V05.Ff y (Real.log (1 - y) - Real.log y) 0 = 1 := by simp [V05.Ff]
  interval_cases k
  · simp [V05.EE_zero, V05.ecoef]
  · simp only [iteratedDeriv_one, V05.deriv1_EE hy0 hy1, V05.EE_zero]
    simp [V05.G1, hF0, V05.ecoef, V05.al1, V05.be1, Finset.sum_range_succ]
    ring
  · simp only [iteratedDeriv_succ, iteratedDeriv_zero, V05.deriv1_EE hy0 hy1,
      V05.deriv2_EE hy0 hy1, V05.EE_zero]
    simp [V05.G1, V05.G2, hF0, V05.ecoef, V05.al1, V05.be1, V05.ga, Finset.sum_range_succ]
    ring
  · simp only [iteratedDeriv_succ, iteratedDeriv_zero, V05.deriv1_EE hy0 hy1,
      V05.deriv2_EE hy0 hy1, V05.deriv3_EE hy0 hy1, V05.EE_zero]
    simp [V05.G1, V05.G2, V05.G3, hF0, V05.ecoef, V05.al1, V05.be1, V05.ga, V05.ga3,
      Finset.sum_range_succ, Nat.factorial]
    ring

namespace V05
open Real MeasureTheory Set

noncomputable def phi (n : ℕ) (u g : ℝ) : ℝ :=
  (1 + g) ^ (-1 - u) * (g ^ (1 - u)) ^ (n - 1) * (1 + g ^ (1 - u)) ^ (1 - 2 * (n : ℝ))

noncomputable def psi (n : ℕ) (u y : ℝ) : ℝ := (y * (1 - y)) ^ (n - 1) * E n u y

lemma exp4 {a b c d e f : ℝ} (K : ℝ) (h : a + (b + (c + d)) = e + f) :
    Real.exp a * (Real.exp b * (Real.exp c * (Real.exp d * K))) = Real.exp e * Real.exp f * K := by
  rw [← mul_assoc, ← mul_assoc, ← mul_assoc, ← Real.exp_add, ← Real.exp_add, ← Real.exp_add,
    ← Real.exp_add, ← h]; ring_nf

lemma pt1 {n : ℕ} (hn : 1 ≤ n) (u : ℝ) {g : ℝ} (hg : 0 < g) :
    |1 / (1 + g) ^ 2| * psi n u (g / (1 + g)) = phi n u g := by
  have h1 : 0 < 1 + g := by linarith
  set X := Real.log g with hX
  set Y := Real.log (1 + g) with hY
  have eg : g = Real.exp X := (Real.exp_log hg).symm
  have e1g : 1 + g = Real.exp Y := (Real.exp_log h1).symm
  have hly : Real.log (g / (1 + g)) = X - Y := Real.log_div hg.ne' h1.ne'
  have h1y : 1 - g / (1 + g) = Real.exp (-Y) := by
    rw [Real.exp_neg, ← e1g]; field_simp; ring
  have hl1y : Real.log (1 - g / (1 + g)) = -Y := by rw [h1y, Real.log_exp]
  have hyy : g / (1 + g) = Real.exp (X - Y) := by rw [Real.exp_sub, ← eg, ← e1g]
  have hh : g ^ (1 - u) = Real.exp (X * (1 - u)) := Real.rpow_def_of_pos hg _
  have hF : f u (g / (1 + g)) = Real.exp (-Y) * (1 + Real.exp (X * (1 - u))) := by
    unfold f
    rw [hly, hl1y, h1y, hyy, ← Real.exp_add]
    rw [mul_add, ← Real.exp_add]; ring_nf
  unfold psi E phi
  rw [hly, hl1y, hF, hh, Real.mul_rpow (Real.exp_pos _).le (by positivity), ← Real.exp_mul,
    Real.rpow_def_of_pos h1, ← hY, h1y, hyy, ← Real.exp_add]
  rw [e1g, abs_of_pos (by positivity), ← Real.exp_nat_mul, ← Real.exp_nat_mul, one_div,
    ← Real.exp_nat_mul, ← Real.exp_neg]
  have hc : ((n - 1 : ℕ) : ℝ) = (n : ℝ) - 1 := by rw [Nat.cast_sub hn]; simp
  rw [hc]
  apply exp4
  push_cast
  ring

lemma pt2 {n : ℕ} (hn : 1 ≤ n) (u : ℝ) {g : ℝ} (hg : 0 < g) :
    |-(1 / (1 + g) ^ 2)| * psi n u (1 / (1 + g)) = phi n u g := by
  have h1 : 0 < 1 + g := by linarith
  set X := Real.log g with hX
  set Y := Real.log (1 + g) with hY
  have eg : g = Real.exp X := (Real.exp_log hg).symm
  have e1g : 1 + g = Real.exp Y := (Real.exp_log h1).symm
  have hly : Real.log (1 / (1 + g)) = -Y := by rw [one_div, Real.log_inv]
  have h1y : 1 - 1 / (1 + g) = Real.exp (X - Y) := by
    rw [Real.exp_sub, ← eg, ← e1g]; field_simp; ring
  have hl1y : Real.log (1 - 1 / (1 + g)) = X - Y := by rw [h1y, Real.log_exp]
  have hyy : 1 / (1 + g) = Real.exp (-Y) := by rw [Real.exp_neg, ← e1g, one_div]
  have hh : g ^ (1 - u) = Real.exp (X * (1 - u)) := Real.rpow_def_of_pos hg _
  have hF : f u (1 / (1 + g)) = Real.exp (-Y + u * X) * (1 + Real.exp (X * (1 - u))) := by
    unfold f
    rw [hly, hl1y, h1y, hyy, ← Real.exp_add, mul_add, ← Real.exp_add]; ring_nf
  unfold psi E phi
  rw [hly, hl1y, hF, hh, Real.mul_rpow (Real.exp_pos _).le (by positivity), ← Real.exp_mul,
    Real.rpow_def_of_pos h1, ← hY, h1y, hyy, ← Real.exp_add]
  rw [e1g, abs_neg, abs_of_pos (by positivity), ← Real.exp_nat_mul, ← Real.exp_nat_mul, one_div,
    ← Real.exp_nat_mul, ← Real.exp_neg]
  have hc : ((n - 1 : ℕ) : ℝ) = (n : ℝ) - 1 := by rw [Nat.cast_sub hn]; simp
  rw [hc]
  apply exp4
  push_cast
  ring

lemma phi_bound {n : ℕ} (hn : 1 ≤ n) {u : ℝ} (hu : |u| < 1) {g : ℝ} (hg0 : 0 < g) (hg1 : g < 1) :
    ‖phi n u g‖ ≤ 1 := by
  have hu' := abs_lt.mp hu
  have hh0 : 0 < g ^ (1 - u) := Real.rpow_pos_of_pos hg0 _
  have hh1 : g ^ (1 - u) ≤ 1 := Real.rpow_le_one hg0.le hg1.le (by linarith)
  have a1 : (1 + g) ^ (-1 - u) ≤ 1 := Real.rpow_le_one_of_one_le_of_nonpos (by linarith) (by linarith)
  have a2 : (g ^ (1 - u)) ^ (n - 1) ≤ 1 := pow_le_one₀ hh0.le hh1
  have hn' : (1 : ℝ) ≤ n := by exact_mod_cast hn
  have a3 : (1 + g ^ (1 - u)) ^ (1 - 2 * (n : ℝ)) ≤ 1 :=
    Real.rpow_le_one_of_one_le_of_nonpos (by linarith) (by linarith)
  have b1 : 0 ≤ (1 + g) ^ (-1 - u) := Real.rpow_nonneg (by linarith) _
  have b2 : 0 ≤ (g ^ (1 - u)) ^ (n - 1) := pow_nonneg hh0.le _
  have b3 : 0 ≤ (1 + g ^ (1 - u)) ^ (1 - 2 * (n : ℝ)) := Real.rpow_nonneg (by linarith) _
  unfold phi
  rw [Real.norm_eq_abs, abs_of_nonneg (by positivity)]
  simpa using mul_le_mul (mul_le_mul a1 a2 b2 zero_le_one) a3 b3 (by norm_num)

lemma phi_int {n : ℕ} (hn : 1 ≤ n) {u : ℝ} (hu : |u| < 1) : IntegrableOn (phi n u) (Ioo 0 1) := by
  refine IntegrableOn.of_bound (by simp) ?_ 1 ?_
  · exact (by unfold phi; fun_prop : Measurable (phi n u)).aestronglyMeasurable
  · exact (ae_restrict_iff' measurableSet_Ioo).mpr
      (Filter.Eventually.of_forall fun g hg => phi_bound hn hu hg.1 hg.2)

lemma img1 : (fun g : ℝ => g / (1 + g)) '' Ioo 0 1 = Ioo 0 (1 / 2) := by
  ext y; constructor
  · rintro ⟨g, ⟨h0, h1⟩, rfl⟩
    refine ⟨by positivity, ?_⟩
    rw [div_lt_iff₀ (by linarith)]; linarith
  · rintro ⟨h0, h1⟩
    refine ⟨y / (1 - y), ⟨by apply div_pos h0; linarith, ?_⟩, ?_⟩
    · rw [div_lt_one (by linarith)]; linarith
    · have : 1 - y ≠ 0 := by linarith
      field_simp; ring

lemma img2 : (fun g : ℝ => 1 / (1 + g)) '' Ioo 0 1 = Ioo (1 / 2) 1 := by
  ext y; constructor
  · rintro ⟨g, ⟨h0, h1⟩, rfl⟩
    refine ⟨?_, ?_⟩
    · rw [lt_div_iff₀ (by linarith)]; linarith
    · rw [div_lt_one (by linarith)]; linarith
  · rintro ⟨h0, h1⟩
    refine ⟨(1 - y) / y, ⟨by apply div_pos <;> linarith, ?_⟩, ?_⟩
    · rw [div_lt_one (by linarith)]; linarith
    · have : y ≠ 0 := by linarith
      field_simp; ring

lemma inj1 : InjOn (fun g : ℝ => g / (1 + g)) (Ioo 0 1) := by
  intro a ha b hb h
  have h1 : (1 + a) ≠ 0 := by linarith [ha.1]
  have h2 : (1 + b) ≠ 0 := by linarith [hb.1]
  simp only at h
  field_simp at h
  linarith

lemma inj2 : InjOn (fun g : ℝ => 1 / (1 + g)) (Ioo 0 1) := by
  intro a ha b hb h
  have h1 : (1 + a) ≠ 0 := by linarith [ha.1]
  have h2 : (1 + b) ≠ 0 := by linarith [hb.1]
  simp only at h
  field_simp at h
  linarith

lemma der1 (g : ℝ) (hg : 0 < g) :
    HasDerivAt (fun g : ℝ => g / (1 + g)) (1 / (1 + g) ^ 2) g := by
  have h1 : (1 + g) ≠ 0 := by linarith
  have := (hasDerivAt_id' g).div ((hasDerivAt_id' g).const_add 1) h1
  refine this.congr_deriv ?_
  field_simp; ring

lemma der2 (g : ℝ) (hg : 0 < g) :
    HasDerivAt (fun g : ℝ => 1 / (1 + g)) (-(1 / (1 + g) ^ 2)) g := by
  have h1 : (1 + g) ≠ 0 := by linarith
  have := (hasDerivAt_const g (1 : ℝ)).div ((hasDerivAt_id' g).const_add 1) h1
  refine this.congr_deriv ?_
  field_simp; ring

lemma half1 {n : ℕ} (hn : 1 ≤ n) {u : ℝ} (hu : |u| < 1) :
    IntegrableOn (psi n u) (Ioo 0 (1 / 2)) ∧
      ∫ y in Ioo 0 (1 / 2), psi n u y = ∫ g in Ioo 0 1, phi n u g := by
  have hd : ∀ x ∈ Ioo (0 : ℝ) 1, HasDerivWithinAt (fun g : ℝ => g / (1 + g)) (1 / (1 + x) ^ 2)
      (Ioo 0 1) x := fun x hx => (der1 x hx.1).hasDerivWithinAt
  have hc : EqOn (fun x => |1 / (1 + x) ^ 2| • psi n u (x / (1 + x))) (phi n u) (Ioo 0 1) :=
    fun x hx => pt1 hn u hx.1
  constructor
  · rw [← img1, integrableOn_image_iff_integrableOn_abs_deriv_smul measurableSet_Ioo hd inj1]
    exact (phi_int hn hu).congr_fun hc.symm measurableSet_Ioo
  · rw [← img1, integral_image_eq_integral_abs_deriv_smul measurableSet_Ioo hd inj1]
    exact setIntegral_congr_fun measurableSet_Ioo hc

lemma half2 {n : ℕ} (hn : 1 ≤ n) {u : ℝ} (hu : |u| < 1) :
    IntegrableOn (psi n u) (Ioo (1 / 2) 1) ∧
      ∫ y in Ioo (1 / 2) 1, psi n u y = ∫ g in Ioo 0 1, phi n u g := by
  have hd : ∀ x ∈ Ioo (0 : ℝ) 1, HasDerivWithinAt (fun g : ℝ => 1 / (1 + g)) (-(1 / (1 + x) ^ 2))
      (Ioo 0 1) x := fun x hx => (der2 x hx.1).hasDerivWithinAt
  have hc : EqOn (fun x => |-(1 / (1 + x) ^ 2)| • psi n u (1 / (1 + x))) (phi n u) (Ioo 0 1) :=
    fun x hx => pt2 hn u hx.1
  constructor
  · rw [← img2, integrableOn_image_iff_integrableOn_abs_deriv_smul measurableSet_Ioo hd inj2]
    exact (phi_int hn hu).congr_fun hc.symm measurableSet_Ioo
  · rw [← img2, integral_image_eq_integral_abs_deriv_smul measurableSet_Ioo hd inj2]
    exact setIntegral_congr_fun measurableSet_Ioo hc

lemma psi_int {n : ℕ} (hn : 1 ≤ n) {u : ℝ} (hu : |u| < 1) :
    IntegrableOn (psi n u) (Ioo 0 1) ∧ ∫ y in Ioo 0 1, psi n u y = 2 * J n u := by
  obtain ⟨i1, e1⟩ := half1 hn hu
  obtain ⟨i2, e2⟩ := half2 hn hu
  have i2' : IntegrableOn (psi n u) (Ico (1 / 2) 1) := (integrableOn_Ico_iff_integrableOn_Ioo (f := psi n u) (a := 1 / 2) (b := 1)).mpr i2
  have hU : Ioo (0 : ℝ) 1 = Ioo 0 (1 / 2) ∪ Ico (1 / 2) 1 := by
    rw [Ioo_union_Ico_eq_Ioo] <;> norm_num
  have hD : Disjoint (Ioo (0 : ℝ) (1 / 2)) (Ico (1 / 2) 1) :=
    Set.disjoint_left.mpr fun x h1 h2 => by linarith [h1.2, h2.1]
  refine ⟨hU ▸ i1.union i2', ?_⟩
  rw [hU, setIntegral_union hD measurableSet_Ico i1 i2', integral_Ico_eq_integral_Ioo, e1, e2]
  unfold J
  rw [intervalIntegral.integral_of_le zero_le_one, integral_Ioc_eq_integral_Ioo, two_mul]
  rfl

end V05

theorem V05.proof_3 : V05.stmt_3 := by
  intro n hn u hu
  have h := (V05.psi_int hn hu).2
  rw [intervalIntegral.integral_of_le zero_le_one, MeasureTheory.integral_Ioc_eq_integral_Ioo]
  change _ = 1 / 2 * ∫ y in Set.Ioo 0 1, V05.psi n u y
  rw [h]; ring

namespace V05
open Real MeasureTheory Set Filter Topology

lemma choose_ratio (r : ℝ) (k : ℕ) :
    ((k : ℝ) + 1) * Ring.choose r (k + 1) = Ring.choose r k * (r - k) := by
  have h := Ring.choose_smul_choose r (n := k + 1) (k := k) (by omega)
  rw [Nat.choose_succ_self_right, Nat.add_sub_cancel_left] at h
  simp only [Ring.choose_one_right, nsmul_eq_mul] at h
  rw [← h]; push_cast; ring

lemma an_rec (j : ℕ) : ((j : ℝ) + 2) * an (j + 1 + 1) = 2 * (2 * j + 1) * an (j + 1) := by
  have h := Nat.succ_mul_centralBinom_succ (j + 1)
  have h' : ((j : ℝ) + 2) * (Nat.choose (2 * (j + 2)) (j + 2) : ℝ)
      = 2 * (2 * ((j : ℝ) + 1) + 1) * (Nat.choose (2 * (j + 1)) (j + 1) : ℝ) := by
    have := congrArg (fun x : ℕ => (x : ℝ)) h
    simp only [Nat.centralBinom_eq_two_mul_choose] at this
    push_cast at this
    linarith
  unfold an
  have h1 : (2 * ((j + 2 : ℕ) : ℝ) - 1) ≠ 0 := by push_cast; linarith [(Nat.cast_nonneg j : (0:ℝ) ≤ j)]
  have h2 : (2 * ((j + 1 : ℕ) : ℝ) - 1) ≠ 0 := by push_cast; linarith [(Nat.cast_nonneg j : (0:ℝ) ≤ j)]
  rw [div_eq_mul_inv, div_eq_mul_inv]
  field_simp
  push_cast at h' ⊢
  linear_combination (2 * (j:ℝ) + 1) * h'

lemma choose_half (j : ℕ) : Ring.choose (1 / 2 : ℝ) (j + 1) * (-4) ^ (j + 1) = -an (j + 1) := by
  induction j with
  | zero => simp [an]; norm_num
  | succ j ih =>
    have hr := choose_ratio (1 / 2 : ℝ) (j + 1)
    have hc := an_rec j
    have h2 : ((j : ℝ) + 2) ≠ 0 := by positivity
    apply mul_left_cancel₀ h2
    push_cast at hr
    have : ((j : ℝ) + 2) * Ring.choose (1 / 2 : ℝ) (j + 1 + 1)
        = Ring.choose (1 / 2 : ℝ) (j + 1) * (1 / 2 - (j + 1)) := by
      rw [← hr]; ring
    have key : ((j : ℝ) + 2) * (Ring.choose (1 / 2 : ℝ) (j + 1 + 1) * (-4) ^ (j + 1 + 1))
        = 2 * (2 * j + 1) * (Ring.choose (1 / 2 : ℝ) (j + 1) * (-4) ^ (j + 1)) := by
      calc _ = (((j : ℝ) + 2) * Ring.choose (1 / 2 : ℝ) (j + 1 + 1)) * (-4) ^ (j + 1 + 1) := by ring
        _ = _ := by rw [this]; ring
    rw [key, ih]; linarith

lemma sqrt_series {x : ℝ} (hx0 : 0 ≤ x) (hx1 : x < 1) :
    HasSum (fun n : ℕ => an (n + 1) / 4 ^ (n + 1) * x ^ (n + 1)) (1 - Real.sqrt (1 - x)) := by
  have hb : (-x) ∈ Metric.eball (0 : ℝ) 1 := by
    rw [Metric.mem_eball, edist_dist, Real.dist_eq, sub_zero, abs_neg, abs_of_nonneg hx0]
    exact ENNReal.ofReal_lt_one.mpr hx1
  have h := (Real.one_add_rpow_hasFPowerSeriesOnBall_zero (a := (1 / 2 : ℝ))).hasSum hb
  simp only [binomialSeries_apply, List.ofFn_const, List.prod_replicate, zero_add,
    smul_eq_mul] at h
  have h1 := (hasSum_nat_add_iff' 1).mpr h
  simp only [Finset.sum_range_one, Ring.choose_zero_right, pow_zero, mul_one] at h1
  have e : (1 + -x) ^ (1 / 2 : ℝ) = Real.sqrt (1 - x) := by
    rw [Real.sqrt_eq_rpow]; ring_nf
  rw [e] at h1
  have h2 := h1.neg
  convert h2 using 1
  · funext j
    have hc := choose_half j
    have : (-x) ^ (j + 1) = (-4) ^ (j + 1) * x ^ (j + 1) / 4 ^ (j + 1) := by
      rw [show (-x) = (-4) * x / 4 by ring, div_pow, mul_pow]
    rw [this, show an (j + 1) = -(Ring.choose (1 / 2 : ℝ) (j + 1) * (-4) ^ (j + 1)) by linarith]
    field_simp
  · ring


noncomputable def I0 (u g : ℝ) : ℝ := ((1 + g) ^ (-1 - u) * (1 + g ^ (1 - u)) - 1) * g ^ (u - 1)

lemma pw {α u g : ℝ} (hα0 : 0 < α) (hα1 : α < π) (hu : |u| < 1) (hg0 : 0 < g) (hg1 : g ≤ 1) :
    HasSum (fun n : ℕ => an (n + 1) * D04.c α ^ (2 * (n + 1)) * phi (n + 1) u g)
      (I0 u g - D04.Fi α u g) := by
  obtain ⟨hc0, hc1, -, -, -⟩ := D04.cs_facts hα0 hα1
  have hu' := abs_lt.mp hu
  set c := D04.c α with hc
  set h := g ^ (1 - u) with hh
  set P := (1 + g) ^ (-1 - u) with hP
  have hh0 : 0 < h := Real.rpow_pos_of_pos hg0 _
  have hh1 : h ≤ 1 := Real.rpow_le_one hg0.le hg1 (by linarith)
  have hw : g ^ (u - 1) = h⁻¹ := by rw [hh, ← Real.rpow_neg hg0.le]; congr 1; ring
  have h2 : g ^ (2 - 2 * u) = h ^ 2 := by
    rw [hh, ← Real.rpow_natCast, ← Real.rpow_mul hg0.le]; congr 1; push_cast; ring
  have hcos : Real.cos α = 2 * c ^ 2 - 1 := by rw [D04.cos_eq_kap]; rfl
  set x := 4 * c ^ 2 * h / (1 + h) ^ 2 with hx
  have h1h : 0 < 1 + h := by linarith
  have hx0 : 0 ≤ x := by positivity
  have hx1 : x < 1 := by
    rw [hx, div_lt_one (by positivity)]
    nlinarith [sq_nonneg (1 - h), mul_pos hh0 (by nlinarith : (0:ℝ) < 1 - c ^ 2)]
  have hrad : 1 + g ^ (2 - 2 * u) - 2 * Real.cos α * g ^ (1 - u) = (1 + h) ^ 2 * (1 - x) := by
    rw [h2, hcos, ← hh, hx]; field_simp; ring
  have hsq : Real.sqrt (1 + g ^ (2 - 2 * u) - 2 * Real.cos α * g ^ (1 - u))
      = (1 + h) * Real.sqrt (1 - x) := by
    rw [hrad, Real.sqrt_mul (by positivity), Real.sqrt_sq h1h.le]
  have hdiff : I0 u g - D04.Fi α u g = P * (1 + h) * h⁻¹ * (1 - Real.sqrt (1 - x)) := by
    unfold I0 D04.Fi
    rw [hsq, hw, ← hh, ← hP]; ring
  rw [hdiff]
  convert (sqrt_series hx0 hx1).mul_left (P * (1 + h) * h⁻¹) using 1
  funext n
  unfold phi
  rw [← hh, ← hP]
  have he : (1 - 2 * ((n + 1 : ℕ) : ℝ)) = -((2 * n + 1 : ℕ) : ℝ) := by push_cast; ring
  rw [he, Real.rpow_neg h1h.le, Real.rpow_natCast, hx]
  simp only [Nat.add_sub_cancel]
  have hne : (1 + h) ≠ 0 := h1h.ne'
  have hne' : h ≠ 0 := hh0.ne'
  rw [div_pow, mul_pow, mul_pow, ← pow_mul, ← pow_mul, pow_succ h]
  field_simp
  ring


lemma I0_eq {u g : ℝ} (hg : 0 < g) :
    I0 u g = ((1 + g) ^ (-1 - u) - 1) * g ^ (u - 1) + (1 + g) ^ (-1 - u) := by
  have : g ^ (1 - u) * g ^ (u - 1) = 1 := by
    rw [← Real.rpow_add hg]; simp
  unfold I0; linear_combination (1 + g) ^ (-1 - u) * this

lemma one_sub_rpow_le {r g : ℝ} (hr : 0 ≤ r) (hg : 0 ≤ g) :
    0 ≤ 1 - (1 + g) ^ (-r) ∧ 1 - (1 + g) ^ (-r) ≤ r * g := by
  have h1 : 0 < 1 + g := by linarith
  have hL0 : 0 ≤ Real.log (1 + g) := Real.log_nonneg (by linarith)
  have hL1 : Real.log (1 + g) ≤ g := by
    have := Real.log_le_sub_one_of_pos h1; linarith
  rw [Real.rpow_def_of_pos h1]
  constructor
  · have : Real.exp (Real.log (1 + g) * -r) ≤ 1 := Real.exp_le_one_iff.mpr (by nlinarith)
    linarith
  · have := Real.add_one_le_exp (Real.log (1 + g) * -r)
    nlinarith

lemma I0_bound {u g : ℝ} (hu : |u| < 1) (hg0 : 0 < g) (_hg1 : g ≤ 1) :
    ‖I0 u g‖ ≤ 2 * g ^ u + 1 := by
  have hu' := abs_lt.mp hu
  rw [I0_eq hg0]
  obtain ⟨a1, a2⟩ := one_sub_rpow_le (r := 1 + u) (by linarith) hg0.le
  have e : -(1 + u) = -1 - u := by ring
  rw [e] at a1 a2
  have hP0 : 0 ≤ (1 + g) ^ (-1 - u) := Real.rpow_nonneg (by linarith) _
  have hP1 : (1 + g) ^ (-1 - u) ≤ 1 := by linarith
  have hw : 0 ≤ g ^ (u - 1) := Real.rpow_nonneg hg0.le _
  have hgw : g * g ^ (u - 1) = g ^ u := by
    rw [mul_comm, Real.rpow_sub_one hg0.ne']; field_simp
  rw [Real.norm_eq_abs]
  have : |((1 + g) ^ (-1 - u) - 1) * g ^ (u - 1)| ≤ 2 * g ^ u := by
    rw [abs_mul, abs_of_nonneg hw, abs_sub_comm, abs_of_nonneg a1, ← hgw]
    have : (1 + u) * g ≤ 2 * g := by nlinarith
    nlinarith
  calc _ ≤ |((1 + g) ^ (-1 - u) - 1) * g ^ (u - 1)| + |(1 + g) ^ (-1 - u)| := abs_add_le _ _
    _ ≤ _ := by rw [abs_of_nonneg hP0]; linarith

lemma I0_meas (u : ℝ) : Measurable (I0 u) := by unfold I0; fun_prop

lemma I0_int {u : ℝ} (hu : |u| < 1) : IntervalIntegrable (I0 u) volume 0 1 := by
  have hu' := abs_lt.mp hu
  have hb : IntervalIntegrable (fun g : ℝ => 2 * g ^ u + 1) volume 0 1 :=
    ((intervalIntegral.intervalIntegrable_rpow' (by linarith)).const_mul 2).add
      intervalIntegrable_const
  refine hb.mono_fun' (I0_meas u).aestronglyMeasurable ?_
  refine (ae_restrict_iff' measurableSet_uIoc).mpr (Eventually.of_forall fun g hg => ?_)
  rw [uIoc_of_le zero_le_one] at hg
  exact I0_bound hu hg.1 hg.2

noncomputable def Qp (u g : ℝ) : ℝ := u⁻¹ * (g ^ u * (1 + g) ^ (-u) - g ^ u - (1 + g) ^ (-u))

lemma Qp_deriv {u g : ℝ} (hu0 : u ≠ 0) (hg : 0 < g) : HasDerivAt (Qp u) (I0 u g) g := by
  have h1 : 0 < 1 + g := by linarith
  have dA : HasDerivAt (fun x : ℝ => x ^ u) (u * g ^ (u - 1)) g :=
    Real.hasDerivAt_rpow_const (Or.inl hg.ne')
  have dB : HasDerivAt (fun x : ℝ => (1 + x) ^ (-u)) (1 * (-u) * (1 + g) ^ (-u - 1)) g :=
    ((hasDerivAt_id' g).const_add 1).rpow_const (Or.inl h1.ne')
  have h := ((dA.mul dB).sub dA).sub dB |>.const_mul u⁻¹
  unfold Qp
  refine h.congr_deriv ?_
  unfold I0
  rw [Real.rpow_sub_one hg.ne', Real.rpow_sub_one h1.ne', show (-1 - u) = -u - 1 by ring,
    Real.rpow_sub_one h1.ne', show (1 - u) = -u + 1 by ring, Real.rpow_add_one hg.ne',
    Real.rpow_neg hg.le]
  have hA : 0 < g ^ u := Real.rpow_pos_of_pos hg _
  field_simp
  ring

lemma K_tendsto {u : ℝ} (hu : |u| < 1) :
    Tendsto (fun g : ℝ => g ^ u * ((1 + g) ^ (-u) - 1)) (𝓝[>] 0) (𝓝 0) := by
  have hu' := abs_lt.mp hu
  have hc : Tendsto (fun g : ℝ => 2 * g ^ (u + 1)) (𝓝[>] 0) (𝓝 0) := by
    have : ContinuousAt (fun g : ℝ => 2 * g ^ (u + 1)) 0 :=
      (Real.continuousAt_rpow_const _ _ (Or.inr (by linarith))).const_mul 2
    have h0 := this.tendsto
    simp only [Real.zero_rpow (show u + 1 ≠ 0 by linarith), mul_zero] at h0
    exact h0.mono_left nhdsWithin_le_nhds
  refine squeeze_zero_norm' ?_ hc
  have hmem : Ioo (0 : ℝ) 1 ∈ 𝓝[>] 0 := Ioo_mem_nhdsGT (by norm_num)
  filter_upwards [hmem] with g hg
  have h1 : 0 < 1 + g := by linarith [hg.1]
  have hL0 : 0 ≤ Real.log (1 + g) := Real.log_nonneg (by linarith [hg.1])
  have hL1 : Real.log (1 + g) ≤ g := by
    have := Real.log_le_sub_one_of_pos h1; linarith
  have ht : |Real.log (1 + g) * -u| ≤ g := by
    rw [abs_mul, abs_neg, abs_of_nonneg hL0]
    calc _ ≤ g * 1 := mul_le_mul hL1 hu.le (abs_nonneg _) hg.1.le
      _ = g := mul_one g
  have he := Real.abs_exp_sub_one_le (ht.trans hg.2.le)
  rw [← Real.rpow_def_of_pos h1] at he
  have hA : 0 < g ^ u := Real.rpow_pos_of_pos hg.1 _
  rw [Real.norm_eq_abs, abs_mul, abs_of_pos hA, Real.rpow_add hg.1, Real.rpow_one]
  calc g ^ u * |(1 + g) ^ (-u) - 1| ≤ g ^ u * (2 * g) := by
        apply mul_le_mul_of_nonneg_left _ hA.le; linarith
    _ = _ := by ring

lemma I0_integral {u : ℝ} (hu : |u| < 1) : ∫ g in (0 : ℝ)..1, I0 u g = 0 := by
  by_cases hu0 : u = 0
  · subst hu0
    have : ∀ g ∈ uIcc (0 : ℝ) 1, I0 0 g = 0 := by
      intro g hg
      rw [uIcc_of_le zero_le_one] at hg
      unfold I0
      have h1 : (1 + g) ≠ 0 := by linarith [hg.1]
      simp [Real.rpow_neg_one, h1]
    rw [intervalIntegral.integral_congr this]; simp
  · have hlim0 : Tendsto (Qp u) (𝓝[>] 0) (𝓝 (u⁻¹ * (0 - 1))) := by
      have hB : Tendsto (fun g : ℝ => (1 + g) ^ (-u)) (𝓝[>] 0) (𝓝 1) := by
        have : ContinuousAt (fun g : ℝ => (1 + g) ^ (-u)) 0 :=
          (continuousAt_const.add continuousAt_id).rpow_const (Or.inl (by norm_num : (1:ℝ) + 0 ≠ 0))
        have h0 := this.tendsto
        simp only [add_zero, Real.one_rpow] at h0
        exact h0.mono_left nhdsWithin_le_nhds
      have := ((K_tendsto hu).sub hB).const_mul u⁻¹
      refine this.congr' (Eventually.of_forall fun g => ?_)
      unfold Qp; ring
    have hlim1 : Tendsto (Qp u) (𝓝[<] 1) (𝓝 (Qp u 1)) := by
      have : ContinuousAt (Qp u) 1 := by
        unfold Qp
        have c1 : ContinuousAt (fun g : ℝ => g ^ u) 1 :=
          Real.continuousAt_rpow_const _ _ (Or.inl (by norm_num))
        have c2 : ContinuousAt (fun g : ℝ => (1 + g) ^ (-u)) 1 :=
          (continuousAt_const.add continuousAt_id).rpow_const (Or.inl (by norm_num : (1:ℝ) + 1 ≠ 0))
        exact ((c1.mul c2).sub c1 |>.sub c2).const_mul _
      exact this.tendsto.mono_left nhdsWithin_le_nhds
    rw [intervalIntegral.integral_eq_sub_of_hasDerivAt_of_tendsto zero_lt_one
      (fun g hg => Qp_deriv hu0 hg.1) (I0_int hu) hlim0 hlim1]
    unfold Qp; simp

lemma diff_bound {α u g : ℝ} (hα0 : 0 < α) (hα1 : α < π) (hu : |u| < 1) (hg0 : 0 < g)
    (_hg1 : g ≤ 1) : I0 u g - D04.Fi α u g ≤ 4 := by
  obtain ⟨hc0, hc1, -, -, -⟩ := D04.cs_facts hα0 hα1
  have hu' := abs_lt.mp hu
  set c := D04.c α with hc
  set h := g ^ (1 - u) with hh
  set P := (1 + g) ^ (-1 - u) with hP
  have hh0 : 0 < h := Real.rpow_pos_of_pos hg0 _
  have hw : g ^ (u - 1) = h⁻¹ := by rw [hh, ← Real.rpow_neg hg0.le]; congr 1; ring
  have h2 : g ^ (2 - 2 * u) = h ^ 2 := by
    rw [hh, ← Real.rpow_natCast, ← Real.rpow_mul hg0.le]; congr 1; push_cast; ring
  have hcos : Real.cos α = 2 * c ^ 2 - 1 := by rw [D04.cos_eq_kap]; rfl
  set x := 4 * c ^ 2 * h / (1 + h) ^ 2 with hx
  have h1h : 0 < 1 + h := by linarith
  have hx0 : 0 ≤ x := by positivity
  have hx1 : x ≤ 1 := by
    rw [hx, div_le_one (by positivity)]
    nlinarith [sq_nonneg (1 - h), mul_pos hh0 (by nlinarith : (0:ℝ) < 1 - c ^ 2)]
  have hrad : 1 + g ^ (2 - 2 * u) - 2 * Real.cos α * g ^ (1 - u) = (1 + h) ^ 2 * (1 - x) := by
    rw [h2, hcos, ← hh, hx]; field_simp; ring
  have hsq : Real.sqrt (1 + g ^ (2 - 2 * u) - 2 * Real.cos α * g ^ (1 - u))
      = (1 + h) * Real.sqrt (1 - x) := by
    rw [hrad, Real.sqrt_mul (by positivity), Real.sqrt_sq h1h.le]
  have hdiff : I0 u g - D04.Fi α u g = P * (1 + h) * h⁻¹ * (1 - Real.sqrt (1 - x)) := by
    unfold I0 D04.Fi
    rw [hsq, hw, ← hh, ← hP]; ring
  have hs : 1 - x ≤ Real.sqrt (1 - x) :=
    Real.le_sqrt_of_sq_le (by nlinarith)
  have hP0 : 0 ≤ P := Real.rpow_nonneg (by linarith) _
  have hP1 : P ≤ 1 := Real.rpow_le_one_of_one_le_of_nonpos (by linarith) (by linarith)
  have hk : 0 ≤ P * (1 + h) * h⁻¹ := by positivity
  rw [hdiff]
  calc P * (1 + h) * h⁻¹ * (1 - Real.sqrt (1 - x)) ≤ P * (1 + h) * h⁻¹ * x :=
        mul_le_mul_of_nonneg_left (by linarith) hk
    _ = 4 * c ^ 2 * P / (1 + h) := by rw [hx]; field_simp
    _ ≤ 4 := by
        rw [div_le_iff₀ h1h]
        have : c ^ 2 ≤ 1 := by nlinarith
        nlinarith [mul_le_mul this hP1 hP0 zero_le_one]

lemma an_nonneg (n : ℕ) : 0 ≤ an (n + 1) := by
  unfold an; apply div_nonneg (by positivity); push_cast; linarith [(Nat.cast_nonneg n : (0:ℝ) ≤ n)]

lemma phi_nonneg (n : ℕ) (u : ℝ) {g : ℝ} (hg : 0 < g) : 0 < phi n u g := by
  unfold phi
  have : 0 < g ^ (1 - u) := Real.rpow_pos_of_pos hg _
  have h1 : 0 < (1 + g) ^ (-1 - u) := Real.rpow_pos_of_pos (by linarith) _
  have h2 : 0 < (1 + g ^ (1 - u)) ^ (1 - 2 * (n : ℝ)) := Real.rpow_pos_of_pos (by linarith) _
  positivity

lemma term_nonneg {α : ℝ} (n : ℕ) (u : ℝ) {g : ℝ} (hg : 0 < g) :
    0 ≤ an (n + 1) * D04.c α ^ (2 * (n + 1)) * phi (n + 1) u g := by
  have := an_nonneg n
  have := phi_nonneg (n + 1) u hg
  have : 0 ≤ D04.c α ^ (2 * (n + 1)) := by rw [pow_mul]; positivity
  positivity

lemma phi_intIoc {n : ℕ} (hn : 1 ≤ n) {u : ℝ} (hu : |u| < 1) :
    IntegrableOn (phi n u) (Ioc 0 1) :=
  (integrableOn_Ioc_iff_integrableOn_Ioo).mpr (phi_int hn hu)

lemma J_eq (n : ℕ) (u : ℝ) : J n u = ∫ g in Ioc (0 : ℝ) 1, phi n u g := by
  unfold J; rw [intervalIntegral.integral_of_le zero_le_one]; rfl

lemma stmt2_all {α : ℝ} (hα0 : 0 < α) (hα1 : α < π) {u : ℝ} (hu : |u| < 1) :
    IntervalIntegrable (D04.Fi α u) volume 0 1 ∧
    (∀ n : ℕ, 1 ≤ n → 0 < J n u) ∧
    HasSum (fun n : ℕ => an (n + 1) * D04.c α ^ (2 * (n + 1)) * J (n + 1) u) (-(D04.Mint u α)) := by
  have hFm : Measurable (D04.Fi α u) := by unfold D04.Fi; fun_prop
  have hD : IntervalIntegrable (fun g => I0 u g - D04.Fi α u g) volume 0 1 := by
    refine (intervalIntegrable_const (c := (4 : ℝ))).mono_fun'
      ((I0_meas u).sub hFm).aestronglyMeasurable ?_
    refine (ae_restrict_iff' measurableSet_uIoc).mpr (Eventually.of_forall fun g hg => ?_)
    rw [uIoc_of_le zero_le_one] at hg
    have hs := pw hα0 hα1 hu hg.1 hg.2
    have h0 : 0 ≤ I0 u g - D04.Fi α u g := hs.nonneg fun n => term_nonneg n u hg.1
    show ‖I0 u g - D04.Fi α u g‖ ≤ 4
    rw [Real.norm_eq_abs, abs_of_nonneg h0]
    exact diff_bound hα0 hα1 hu hg.1 hg.2
  have hFi : IntervalIntegrable (D04.Fi α u) volume 0 1 := by
    have e : D04.Fi α u = fun g => I0 u g - (I0 u g - D04.Fi α u g) := by funext g; ring
    rw [e]; exact (I0_int hu).sub hD
  refine ⟨hFi, ?_, ?_⟩
  · intro n hn
    unfold J
    exact intervalIntegral.intervalIntegral_pos_of_pos_on
      ((intervalIntegrable_iff_integrableOn_Ioo_of_le zero_le_one).mpr (phi_int hn hu))
      (fun g hg => phi_nonneg n u hg.1) zero_lt_one
  · have hmct := integral_tendsto_of_tendsto_of_monotone (μ := volume.restrict (Ioc (0 : ℝ) 1))
      (f := fun N g => ∑ n ∈ Finset.range N, an (n + 1) * D04.c α ^ (2 * (n + 1)) * phi (n + 1) u g)
      (F := fun g => I0 u g - D04.Fi α u g)
      (fun N => integrable_finsetSum _ fun n _ => (phi_intIoc (by omega) hu).const_mul _)
      hD.1
      ((ae_restrict_iff' measurableSet_Ioc).mpr (Eventually.of_forall fun g hg =>
        monotone_nat_of_le_succ fun N => by
          simp only [Finset.sum_range_succ]
          linarith [term_nonneg (α := α) N u hg.1]))
      ((ae_restrict_iff' measurableSet_Ioc).mpr (Eventually.of_forall fun g hg =>
        (pw hα0 hα1 hu hg.1 hg.2).tendsto_sum_nat))
    have e1 : ∀ N, (∫ g in Ioc (0 : ℝ) 1, ∑ n ∈ Finset.range N,
        an (n + 1) * D04.c α ^ (2 * (n + 1)) * phi (n + 1) u g)
        = ∑ n ∈ Finset.range N, an (n + 1) * D04.c α ^ (2 * (n + 1)) * J (n + 1) u := by
      intro N
      rw [integral_finsetSum _ fun n _ => (phi_intIoc (by omega) hu).const_mul _]
      refine Finset.sum_congr rfl fun n _ => ?_
      rw [integral_const_mul, J_eq]
    have e2 : (∫ g in Ioc (0 : ℝ) 1, I0 u g - D04.Fi α u g) = -(D04.Mint u α) := by
      rw [integral_sub (I0_int hu).1 hFi.1]
      have h1 := I0_integral hu
      rw [intervalIntegral.integral_of_le zero_le_one] at h1
      rw [h1]
      unfold D04.Mint
      rw [intervalIntegral.integral_of_le zero_le_one, zero_sub]
      rfl
    simp only [e1, e2] at hmct
    refine (hasSum_iff_tendsto_nat_of_nonneg (fun n => ?_) _).mpr hmct
    have := an_nonneg n
    have : 0 < J (n + 1) u := by
      unfold J
      exact intervalIntegral.intervalIntegral_pos_of_pos_on
        ((intervalIntegrable_iff_integrableOn_Ioo_of_le zero_le_one).mpr (phi_int (by omega) hu))
        (fun g hg => phi_nonneg (n + 1) u hg.1) zero_lt_one
    have : 0 ≤ D04.c α ^ (2 * (n + 1)) := by rw [pow_mul]; positivity
    positivity

end V05

theorem V05.proof_2 : V05.stmt_2 := by
  intro α hα0 hα1 u hu
  exact V05.stmt2_all hα0 hα1 hu

namespace V05
open Real MeasureTheory Set Filter Topology

noncomputable def ell (y : ℝ) : ℝ := Real.log (1 - y) - Real.log y
noncomputable def AA (n : ℕ) (y : ℝ) : ℝ := -((n : ℝ) - 1) * Real.log y + (n : ℝ) * Real.log (1 - y)
noncomputable def pif (u y : ℝ) : ℝ := y * Real.exp (u * ell y) / f u y
noncomputable def fam (n : ℕ) (z : ℕ × ℕ × ℕ) (u y : ℝ) : ℝ :=
  AA n y ^ z.1 * ell y ^ z.2.1 * pif u y ^ z.2.2 * E n u y

noncomputable def step (n : ℕ) (z : ℕ × ℕ × ℕ) : List (ℝ × (ℕ × ℕ × ℕ)) :=
  [(1, (z.1 + 1, z.2.1, z.2.2)), ((1 - 2 * (n : ℝ)) - z.2.2, (z.1, z.2.1 + 1, z.2.2 + 1)),
    ((z.2.2 : ℝ), (z.1, z.2.1 + 1, z.2.2))]

noncomputable def Dl (n : ℕ) (L : List (ℝ × (ℕ × ℕ × ℕ))) : List (ℝ × (ℕ × ℕ × ℕ)) :=
  L.flatMap (fun p => (step n p.2).map (fun q => (p.1 * q.1, q.2)))

noncomputable def comb (G : ℕ × ℕ × ℕ → ℝ → ℝ) (L : List (ℝ × (ℕ × ℕ × ℕ))) (u : ℝ) : ℝ :=
  (L.map (fun p => p.1 * G p.2 u)).sum

lemma f_eq (u y : ℝ) : f u y = 1 - y + y * Real.exp (u * ell y) := rfl

lemma E_eq (n : ℕ) (u y : ℝ) : E n u y = Real.exp (u * AA n y) * f u y ^ (1 - 2 * (n : ℝ)) := rfl

lemma f_bounds {u y : ℝ} (hy0 : 0 < y) (hy1 : y < 1) :
    Real.exp (-|u * ell y|) ≤ f u y ∧ f u y ≤ Real.exp |u * ell y| := by
  rw [f_eq]
  set v := u * ell y
  have a1 : Real.exp (-|v|) ≤ 1 := Real.exp_le_one_iff.mpr (by linarith [abs_nonneg v])
  have a2 : Real.exp (-|v|) ≤ Real.exp v := Real.exp_le_exp.mpr (by linarith [neg_abs_le v])
  have b1 : 1 ≤ Real.exp |v| := Real.one_le_exp (abs_nonneg v)
  have b2 : Real.exp v ≤ Real.exp |v| := Real.exp_le_exp.mpr (le_abs_self v)
  constructor <;> nlinarith

lemma f_pos {u y : ℝ} (hy0 : 0 < y) (hy1 : y < 1) : 0 < f u y :=
  lt_of_lt_of_le (Real.exp_pos _) (f_bounds hy0 hy1).1

lemma pif_pos {u y : ℝ} (hy0 : 0 < y) (hy1 : y < 1) : 0 < pif u y := by
  unfold pif; have := f_pos (u := u) hy0 hy1; positivity

lemma pif_lt_one {u y : ℝ} (hy0 : 0 < y) (hy1 : y < 1) : pif u y < 1 := by
  unfold pif; have := f_pos (u := u) hy0 hy1
  rw [div_lt_one this, f_eq]; linarith

lemma fam_hasDerivAt (n : ℕ) (z : ℕ × ℕ × ℕ) {y : ℝ} (hy0 : 0 < y) (hy1 : y < 1) (u : ℝ) :
    HasDerivAt (fun u => fam n z u y) (comb (fun z u => fam n z u y) (step n z) u) u := by
  obtain ⟨a, b, j⟩ := z
  have hF := f_pos (u := u) hy0 hy1
  have hX : HasDerivAt (fun u => Real.exp (u * ell y)) (Real.exp (u * ell y) * ell y) u := by
    simpa using ((hasDerivAt_id' u).mul_const (ell y)).exp
  have hf : HasDerivAt (fun u => f u y) (y * (Real.exp (u * ell y) * ell y)) u := by
    have := (hX.const_mul y).const_add (1 - y)
    exact this.congr_of_eventuallyEq (Eventually.of_forall fun v => rfl)
  have hP' : pif u y = y * Real.exp (u * ell y) / f u y := rfl
  have hF' : f u y ≠ 0 := hF.ne'
  have hS : HasDerivAt (fun u => f u y ^ (1 - 2 * (n : ℝ)))
      ((1 - 2 * (n : ℝ)) * ell y * pif u y * f u y ^ (1 - 2 * (n : ℝ))) u := by
    refine (hf.rpow_const (p := 1 - 2 * (n : ℝ)) (Or.inl hF')).congr_deriv ?_
    rw [Real.rpow_sub_one hF', hP']; field_simp
  have hEx : HasDerivAt (fun u => Real.exp (u * AA n y)) (Real.exp (u * AA n y) * AA n y) u := by
    simpa using ((hasDerivAt_id' u).mul_const (AA n y)).exp
  have hP : HasDerivAt (fun u => pif u y) (ell y * pif u y * (1 - pif u y)) u := by
    have := (hX.const_mul y).div hf hF'
    refine (this.congr_of_eventuallyEq (Eventually.of_forall fun v => rfl)).congr_deriv ?_
    rw [hP']; field_simp
  have hPj : HasDerivAt (fun u => pif u y ^ j)
      ((j : ℝ) * pif u y ^ (j - 1) * (ell y * pif u y * (1 - pif u y))) u :=
    hP.fun_pow j
  have h := ((hPj.mul (hEx.mul hS)).const_mul (AA n y ^ a * ell y ^ b))
  have heq : (fun u => fam n (a, b, j) u y) =
      fun u => AA n y ^ a * ell y ^ b * (pif u y ^ j *
        (Real.exp (u * AA n y) * f u y ^ (1 - 2 * (n : ℝ)))) := by
    funext v; simp only [fam, E_eq]; ring
  rw [heq]
  refine h.congr_deriv ?_
  simp only [Pi.mul_apply]
  simp only [comb, step, fam, E_eq, List.map_cons, List.map_nil, List.sum_cons, List.sum_nil]
  rcases j with _ | j
  · simp only [pow_zero, Nat.cast_zero, zero_mul, zero_add, add_zero, sub_zero, one_mul, mul_one]
    ring
  · simp only [Nat.add_sub_cancel, Nat.cast_add, Nat.cast_one]
    ring

lemma comb_nil (G : ℕ × ℕ × ℕ → ℝ → ℝ) : comb G [] = fun _ => 0 := by
  funext u; simp [comb]

lemma comb_cons (G : ℕ × ℕ × ℕ → ℝ → ℝ) (p : ℝ × (ℕ × ℕ × ℕ)) (L : List (ℝ × (ℕ × ℕ × ℕ))) :
    comb G (p :: L) = fun u => p.1 * G p.2 u + comb G L u := by
  funext u; simp [comb]

lemma comb_Dl_cons (n : ℕ) (G : ℕ × ℕ × ℕ → ℝ → ℝ) (p : ℝ × (ℕ × ℕ × ℕ))
    (L : List (ℝ × (ℕ × ℕ × ℕ))) (u : ℝ) :
    comb G (Dl n (p :: L)) u = p.1 * comb G (step n p.2) u + comb G (Dl n L) u := by
  simp only [Dl, List.flatMap_cons, comb, List.map_append, List.sum_append, List.map_map]
  congr 1
  rw [← List.sum_map_mul_left]
  congr 1
  apply List.map_congr_left
  intro q _
  simp only [Function.comp]; ring

lemma comb_hasDerivAt (n : ℕ) {G : ℕ × ℕ × ℕ → ℝ → ℝ} {u : ℝ}
    (hG : ∀ z, HasDerivAt (G z) (comb G (step n z) u) u) :
    ∀ L, HasDerivAt (comb G L) (comb G (Dl n L) u) u := by
  intro L
  induction L with
  | nil => rw [comb_nil]; simpa [Dl, comb] using hasDerivAt_const u (0 : ℝ)
  | cons p L ih =>
    rw [comb_cons, comb_Dl_cons]
    exact ((hG p.2).const_mul p.1).add ih

lemma comb_iter_global (n : ℕ) {G : ℕ × ℕ × ℕ → ℝ → ℝ}
    (hG : ∀ z u, HasDerivAt (G z) (comb G (step n z) u) u) :
    ∀ k L, iteratedDeriv k (comb G L) = comb G ((Dl n)^[k] L) := by
  intro k
  induction k with
  | zero => intro L; simp
  | succ k ih =>
    intro L
    rw [iteratedDeriv_succ']
    have : deriv (comb G L) = comb G (Dl n L) := by
      funext u; exact (comb_hasDerivAt n (fun z => hG z u) L).deriv
    rw [this, ih, Function.iterate_succ_apply]

lemma comb_iter_local (n : ℕ) {G : ℕ × ℕ × ℕ → ℝ → ℝ} {S : Set ℝ} (hS : IsOpen S)
    (hG : ∀ z, ∀ u ∈ S, HasDerivAt (G z) (comb G (step n z) u) u) (L : List (ℝ × (ℕ × ℕ × ℕ))) :
    ∀ k, ∀ u ∈ S, iteratedDeriv k (comb G L) u = comb G ((Dl n)^[k] L) u := by
  intro k
  induction k with
  | zero => intro u _; simp
  | succ k ih =>
    intro u hu
    rw [iteratedDeriv_succ]
    have hev : iteratedDeriv k (comb G L) =ᶠ[𝓝 u] comb G ((Dl n)^[k] L) :=
      Filter.eventually_of_mem (hS.mem_nhds hu) fun v hv => ih v hv
    rw [hev.deriv_eq, (comb_hasDerivAt n (fun z => hG z u hu) _).deriv,
      Function.iterate_succ_apply']

lemma comb_contDiffOn (n : ℕ) {G : ℕ × ℕ × ℕ → ℝ → ℝ} {S : Set ℝ} (hS : IsOpen S)
    (hG : ∀ z, ∀ u ∈ S, HasDerivAt (G z) (comb G (step n z) u) u) :
    ∀ m : ℕ, ∀ L, ContDiffOn ℝ m (comb G L) S := by
  intro m
  induction m with
  | zero =>
    intro L
    rw [Nat.cast_zero, contDiffOn_zero]
    intro u hu
    exact (comb_hasDerivAt n (fun z => hG z u hu) L).continuousAt.continuousWithinAt
  | succ m ih =>
    intro L
    rw [show ((m + 1 : ℕ) : WithTop ℕ∞) = (m : WithTop ℕ∞) + 1 by push_cast; rfl,
      contDiffOn_succ_iff_deriv_of_isOpen hS]
    refine ⟨fun u hu => (comb_hasDerivAt n (fun z => hG z u hu) L).differentiableAt.differentiableWithinAt,
      fun h => absurd h (by simp), ?_⟩
    refine (ih (Dl n L)).congr fun u hu => ?_
    exact (comb_hasDerivAt n (fun z => hG z u hu) L).deriv

noncomputable def dl (_n : ℕ) : ℝ := 1 / 12
noncomputable def ww (n : ℕ) (y : ℝ) : ℝ := (y * (1 - y)) ^ (n - 1)
noncomputable def TT (y : ℝ) : ℝ := |Real.log y| + |Real.log (1 - y)|
noncomputable def qq (t : ℝ) : ℝ := (t * (1 - t)) ^ (-(5 / 8 : ℝ))

lemma dl_pos {n : ℕ} (_hn : 1 ≤ n) : 0 < dl n := by unfold dl; norm_num

lemma TT_eq {y : ℝ} (hy0 : 0 < y) (hy1 : y < 1) : TT y = -Real.log (y * (1 - y)) := by
  have h1 : 0 < 1 - y := by linarith
  unfold TT
  rw [abs_of_nonpos (Real.log_nonpos hy0.le hy1.le),
    abs_of_nonpos (Real.log_nonpos h1.le (by linarith)), Real.log_mul hy0.ne' h1.ne']
  ring

lemma E_bdG {n : ℕ} (hn : 1 ≤ n) {u y d : ℝ} (hy0 : 0 < y) (hy1 : y < 1) (_hd : 0 ≤ d)
    (hu : |u| ≤ d) : 0 < E n u y ∧ E n u y ≤ (y * (1 - y)) ^ (-((3 * (n : ℝ) - 1) * d)) := by
  have h1 : 0 < 1 - y := by linarith
  have hy : 0 < y * (1 - y) := mul_pos hy0 h1
  have hF := f_pos (u := u) hy0 hy1
  obtain ⟨fb1, fb2⟩ := f_bounds (u := u) hy0 hy1
  have hn' : (1 : ℝ) ≤ n := by exact_mod_cast hn
  have hlf : |Real.log (f u y)| ≤ |u| * TT y := by
    have a1 : -|u * ell y| ≤ Real.log (f u y) := (Real.le_log_iff_exp_le hF).mpr fb1
    have a2 : Real.log (f u y) ≤ |u * ell y| := (Real.log_le_iff_le_exp hF).mpr fb2
    have hl : |ell y| ≤ TT y := by
      unfold ell TT; exact (abs_sub _ _).trans (le_of_eq (add_comm _ _))
    have : |u * ell y| ≤ |u| * TT y := by
      rw [abs_mul]; exact mul_le_mul_of_nonneg_left hl (abs_nonneg u)
    rw [abs_le]; constructor <;> linarith
  have hA : |AA n y| ≤ n * TT y := by
    unfold AA TT
    calc |-((n : ℝ) - 1) * Real.log y + n * Real.log (1 - y)|
        ≤ |-((n : ℝ) - 1) * Real.log y| + |(n : ℝ) * Real.log (1 - y)| := abs_add_le _ _
      _ = ((n : ℝ) - 1) * |Real.log y| + n * |Real.log (1 - y)| := by
        rw [abs_mul, abs_mul, abs_neg, abs_of_nonneg (show (0:ℝ) ≤ (n : ℝ) - 1 by linarith),
          Nat.abs_cast]
      _ ≤ (n : ℝ) * (|Real.log y| + |Real.log (1 - y)|) := by
        nlinarith [abs_nonneg (Real.log y), abs_nonneg (Real.log (1 - y))]
  have hT0 : 0 ≤ TT y := by unfold TT; positivity
  have he : E n u y = Real.exp (u * AA n y + Real.log (f u y) * (1 - 2 * (n : ℝ))) := by
    rw [E_eq, Real.rpow_def_of_pos hF, ← Real.exp_add]
  refine ⟨by rw [he]; exact Real.exp_pos _, ?_⟩
  rw [he, Real.rpow_def_of_pos hy]
  apply Real.exp_le_exp.mpr
  have k1 : u * AA n y ≤ |u| * (n * TT y) := by
    calc u * AA n y ≤ |u * AA n y| := le_abs_self _
      _ = |u| * |AA n y| := abs_mul _ _
      _ ≤ _ := mul_le_mul_of_nonneg_left hA (abs_nonneg u)
  have k2 : Real.log (f u y) * (1 - 2 * (n : ℝ)) ≤ |u| * TT y * (2 * n - 1) := by
    have : Real.log (f u y) * (1 - 2 * (n : ℝ)) = -(Real.log (f u y)) * (2 * n - 1) := by ring
    rw [this]
    exact mul_le_mul_of_nonneg_right ((neg_le_abs _).trans hlf) (by linarith)
  have k3 : |u| * (3 * n - 1) ≤ d * (3 * n - 1) :=
    mul_le_mul_of_nonneg_right hu (show (0:ℝ) ≤ 3 * n - 1 by linarith)
  rw [TT_eq hy0 hy1] at k1 k2 hT0
  nlinarith [abs_nonneg u]

lemma wwE_bd {n : ℕ} (hn : 1 ≤ n) {u y d : ℝ} (hy0 : 0 < y) (hy1 : y < 1) (hd : 0 ≤ d)
    (hu : |u| ≤ d) : 0 ≤ ww n y * E n u y ∧
      ww n y * E n u y ≤ (y * (1 - y)) ^ ((n : ℝ) - 1 - (3 * (n : ℝ) - 1) * d) := by
  have h1 : 0 < 1 - y := by linarith
  have hy : 0 < y * (1 - y) := mul_pos hy0 h1
  obtain ⟨hE0, hE1⟩ := E_bdG hn hy0 hy1 hd hu
  have hw : ww n y = (y * (1 - y)) ^ ((n : ℝ) - 1) := by
    unfold ww
    rw [← Real.rpow_natCast, Nat.cast_sub hn, Nat.cast_one]
  have hw0 : 0 ≤ ww n y := by unfold ww; positivity
  refine ⟨mul_nonneg hw0 hE0.le, ?_⟩
  calc ww n y * E n u y ≤ ww n y * (y * (1 - y)) ^ (-((3 * (n : ℝ) - 1) * d)) :=
        mul_le_mul_of_nonneg_left hE1 hw0
    _ = _ := by rw [hw, ← Real.rpow_add hy]; ring_nf

lemma log_bd_eps {x ε : ℝ} (hx0 : 0 < x) (hx1 : x ≤ 1) (hε : 0 < ε) :
    |Real.log x| ≤ x ^ (-ε) / ε := by
  have h := Real.neg_rpow_div_le_log hx0.le hε
  have hl : Real.log x ≤ 0 := Real.log_nonpos hx0.le hx1
  rw [abs_of_nonpos hl]
  have : -x ^ (-ε) / ε = -(x ^ (-ε) / ε) := by ring
  linarith

lemma TT_pow_bd (N : ℕ) : ∃ K : ℝ, 0 ≤ K ∧ ∀ t : ℝ, 0 < t → t < 1 →
    TT t ^ N ≤ K * (t * (1 - t)) ^ (-(1 / 8 : ℝ)) := by
  set ε : ℝ := 1 / (8 * ((N : ℝ) + 1)) with hε
  have hε0 : 0 < ε := by positivity
  refine ⟨(2 / ε) ^ N, by positivity, fun t ht0 ht1 => ?_⟩
  have h1 : 0 < 1 - t := by linarith
  have hy : 0 < t * (1 - t) := mul_pos ht0 h1
  have hy1 : t * (1 - t) ≤ 1 := by nlinarith
  have a := log_bd_eps ht0 ht1.le hε0
  have b := log_bd_eps h1 (by linarith) hε0
  have a' : t ^ (-ε) ≤ (t * (1 - t)) ^ (-ε) :=
    Real.rpow_le_rpow_of_nonpos hy (by nlinarith) (by linarith)
  have b' : (1 - t) ^ (-ε) ≤ (t * (1 - t)) ^ (-ε) :=
    Real.rpow_le_rpow_of_nonpos hy (by nlinarith) (by linarith)
  have hsum : TT t ≤ 2 / ε * (t * (1 - t)) ^ (-ε) := by
    have ha'' : |Real.log t| ≤ (t * (1 - t)) ^ (-ε) / ε :=
      a.trans (div_le_div_of_nonneg_right a' hε0.le)
    have hb'' : |Real.log (1 - t)| ≤ (t * (1 - t)) ^ (-ε) / ε :=
      b.trans (div_le_div_of_nonneg_right b' hε0.le)
    have : 2 / ε * (t * (1 - t)) ^ (-ε) = (t * (1 - t)) ^ (-ε) / ε + (t * (1 - t)) ^ (-ε) / ε := by
      ring
    unfold TT; linarith
  have hpos : 0 ≤ TT t := by unfold TT; positivity
  calc TT t ^ N ≤ (2 / ε * (t * (1 - t)) ^ (-ε)) ^ N := pow_le_pow_left₀ hpos hsum N
    _ = (2 / ε) ^ N * (t * (1 - t)) ^ (-ε * N) := by
        rw [mul_pow, ← Real.rpow_mul_natCast hy.le]
    _ ≤ (2 / ε) ^ N * (t * (1 - t)) ^ (-(1 / 8 : ℝ)) := by
        refine mul_le_mul_of_nonneg_left (Real.rpow_le_rpow_of_exponent_ge hy hy1 ?_) (by positivity)
        rw [hε]
        have hN : (0 : ℝ) ≤ N := N.cast_nonneg
        rw [neg_mul, neg_le_neg_iff, div_mul_eq_mul_div, one_mul, div_le_div_iff₀ (by positivity)
          (by norm_num)]
        nlinarith

lemma qq_int : IntervalIntegrable qq volume 0 1 := by
  have ha : IntervalIntegrable (fun t : ℝ => t ^ (-(5 / 8 : ℝ))) volume 0 1 :=
    intervalIntegral.intervalIntegrable_rpow' (by norm_num)
  have hb : IntervalIntegrable (fun t : ℝ => (1 - t) ^ (-(5 / 8 : ℝ))) volume 0 1 := by
    have := ha.comp_sub_left 1
    simpa using this.symm
  have hab := (ha.add hb).const_mul 2
  refine IntervalIntegrable.mono_fun' hab
    (by unfold qq; exact (by fun_prop : Measurable fun t : ℝ => (t * (1 - t)) ^ (-(5 / 8 : ℝ))).aestronglyMeasurable) ?_
  rw [uIoc_of_le zero_le_one, ← restrict_Ioo_eq_restrict_Ioc]
  refine ae_restrict_of_forall_mem measurableSet_Ioo fun t ht => ?_
  have ht0 := ht.1; have ht1 := ht.2
  have h1 : 0 < 1 - t := by linarith
  have hq : 0 ≤ qq t := by unfold qq; positivity
  dsimp only
  rw [Real.norm_eq_abs, abs_of_nonneg hq, qq, Real.mul_rpow ht0.le h1.le]
  have hA : 0 < t ^ (-(5 / 8 : ℝ)) := Real.rpow_pos_of_pos ht0 _
  have hB : 0 < (1 - t) ^ (-(5 / 8 : ℝ)) := Real.rpow_pos_of_pos h1 _
  have h2 : (1 / 2 : ℝ) ^ (-(5 / 8 : ℝ)) ≤ 2 := by
    have := Real.rpow_le_rpow_of_exponent_ge (show (0:ℝ) < 1 / 2 by norm_num)
      (show (1 / 2 : ℝ) ≤ 1 by norm_num) (show (-1 : ℝ) ≤ -(5 / 8) by norm_num)
    rw [Real.rpow_neg_one] at this; norm_num at this ⊢; linarith
  rcases le_total t (1 / 2) with h | h
  · have : (1 - t) ^ (-(5 / 8 : ℝ)) ≤ 2 :=
      (Real.rpow_le_rpow_of_nonpos (by norm_num) (by linarith) (by norm_num)).trans h2
    nlinarith
  · have : t ^ (-(5 / 8 : ℝ)) ≤ 2 :=
      (Real.rpow_le_rpow_of_nonpos (by norm_num) h (by norm_num)).trans h2
    nlinarith

lemma fam_bdG (N : ℕ) : ∃ K : ℝ, 0 ≤ K ∧ ∀ n : ℕ, 1 ≤ n → ∀ z : ℕ × ℕ × ℕ, z.1 + z.2.1 = N →
    ∀ d u y : ℝ, 0 ≤ d → d ≤ 1 / 12 → |u| ≤ d → 0 < y → y < 1 →
      |ww n y * fam n z u y| ≤
        (n : ℝ) ^ N * K * qq y * (1 / 4 : ℝ) ^ ((n : ℝ) - 1 / 2 - (3 * (n : ℝ) - 1) * d) := by
  obtain ⟨K, hK0, hK⟩ := TT_pow_bd N
  refine ⟨K, hK0, fun n hn z hz d u y hd0 hd1 hu hy0 hy1 => ?_⟩
  obtain ⟨a, b, j⟩ := z
  simp only at hz
  have h1 : 0 < 1 - y := by linarith
  have hy : 0 < y * (1 - y) := mul_pos hy0 h1
  have hy4 : y * (1 - y) ≤ 1 / 4 := by nlinarith [sq_nonneg (y - 1 / 2)]
  obtain ⟨hwE0, hwE1⟩ := wwE_bd hn hy0 hy1 hd0 hu
  have hP0 := pif_pos (u := u) hy0 hy1
  have hP1 := pif_lt_one (u := u) hy0 hy1
  have hn' : (1 : ℝ) ≤ n := by exact_mod_cast hn
  have hT0 : 0 ≤ TT y := by unfold TT; positivity
  have hl : |ell y| ≤ TT y := by
    unfold ell TT; exact (abs_sub _ _).trans (le_of_eq (add_comm _ _))
  have hA : |AA n y| ≤ n * TT y := by
    unfold AA TT
    calc |-((n : ℝ) - 1) * Real.log y + n * Real.log (1 - y)|
        ≤ |-((n : ℝ) - 1) * Real.log y| + |(n : ℝ) * Real.log (1 - y)| := abs_add_le _ _
      _ = ((n : ℝ) - 1) * |Real.log y| + n * |Real.log (1 - y)| := by
        rw [abs_mul, abs_mul, abs_neg, abs_of_nonneg (show (0:ℝ) ≤ (n : ℝ) - 1 by linarith),
          Nat.abs_cast]
      _ ≤ (n : ℝ) * (|Real.log y| + |Real.log (1 - y)|) := by
        nlinarith [abs_nonneg (Real.log y), abs_nonneg (Real.log (1 - y))]
  have hpj : |pif u y ^ j| ≤ 1 := by
    rw [abs_pow, abs_of_pos hP0]; exact pow_le_one₀ hP0.le hP1.le
  have hpoly : |AA n y ^ a * ell y ^ b| ≤ (n : ℝ) ^ N * TT y ^ N := by
    rw [abs_mul, abs_pow, abs_pow, ← hz, pow_add, pow_add]
    have e1 : |AA n y| ^ a ≤ (n * TT y) ^ a := pow_le_pow_left₀ (abs_nonneg _) hA a
    have e2 : |ell y| ^ b ≤ TT y ^ b := pow_le_pow_left₀ (abs_nonneg _) hl b
    have e3 : (1 : ℝ) ≤ (n : ℝ) ^ b := one_le_pow₀ hn'
    rw [mul_pow] at e1
    calc |AA n y| ^ a * |ell y| ^ b ≤ ((n : ℝ) ^ a * TT y ^ a) * TT y ^ b :=
          mul_le_mul e1 e2 (by positivity) (by positivity)
      _ ≤ ((n : ℝ) ^ a * TT y ^ a) * TT y ^ b * (n : ℝ) ^ b := by
          have : 0 ≤ ((n : ℝ) ^ a * TT y ^ a) * TT y ^ b := by positivity
          nlinarith
      _ = _ := by ring
  set e : ℝ := (n : ℝ) - 1 / 2 - (3 * (n : ℝ) - 1) * d with he
  have he0 : 0 ≤ e := by rw [he]; nlinarith
  have hX : (y * (1 - y)) ^ (-(1 / 8 : ℝ)) * (y * (1 - y)) ^ ((n : ℝ) - 1 - (3 * (n : ℝ) - 1) * d)
      = qq y * (y * (1 - y)) ^ e := by
    rw [qq, ← Real.rpow_add hy, ← Real.rpow_add hy, he]; ring_nf
  have hXe : (y * (1 - y)) ^ e ≤ (1 / 4 : ℝ) ^ e := Real.rpow_le_rpow hy.le hy4 he0
  have hq0 : 0 ≤ qq y := by unfold qq; positivity
  have hfam : ww n y * fam n (a, b, j) u y = (AA n y ^ a * ell y ^ b) * pif u y ^ j * (ww n y * E n u y) := by
    unfold fam; ring
  rw [hfam, abs_mul, abs_mul, abs_of_nonneg hwE0]
  calc |AA n y ^ a * ell y ^ b| * |pif u y ^ j| * (ww n y * E n u y)
      ≤ ((n : ℝ) ^ N * TT y ^ N) * 1 * (y * (1 - y)) ^ ((n : ℝ) - 1 - (3 * (n : ℝ) - 1) * d) :=
        mul_le_mul (mul_le_mul hpoly hpj (abs_nonneg _) (by positivity)) hwE1 hwE0 (by positivity)
    _ ≤ (n : ℝ) ^ N * (K * (y * (1 - y)) ^ (-(1 / 8 : ℝ)))
          * (y * (1 - y)) ^ ((n : ℝ) - 1 - (3 * (n : ℝ) - 1) * d) := by
        rw [mul_one]
        gcongr
        exact hK y hy0 hy1
    _ = (n : ℝ) ^ N * K * (qq y * (y * (1 - y)) ^ e) := by rw [← hX]; ring
    _ ≤ (n : ℝ) ^ N * K * (qq y * (1 / 4 : ℝ) ^ e) := by gcongr
    _ = _ := by ring

lemma fam_bd {n : ℕ} (hn : 1 ≤ n) (z : ℕ × ℕ × ℕ) : ∃ B : ℝ, 0 ≤ B ∧ ∀ u y : ℝ, |u| ≤ dl n →
    0 < y → y < 1 → |ww n y * fam n z u y| ≤ B * qq y := by
  obtain ⟨K, hK0, hK⟩ := fam_bdG (z.1 + z.2.1)
  refine ⟨(n : ℝ) ^ (z.1 + z.2.1) * K, by positivity, fun u y hu hy0 hy1 => ?_⟩
  have hn' : (1 : ℝ) ≤ n := by exact_mod_cast hn
  have h := hK n hn z rfl (1 / 12) u y (by norm_num) le_rfl hu hy0 hy1
  have he0 : 0 ≤ (n : ℝ) - 1 / 2 - (3 * (n : ℝ) - 1) * (1 / 12) := by nlinarith
  have h4 : (1 / 4 : ℝ) ^ ((n : ℝ) - 1 / 2 - (3 * (n : ℝ) - 1) * (1 / 12)) ≤ 1 :=
    Real.rpow_le_one (by norm_num) (by norm_num) he0
  have hq0 : 0 ≤ qq y := by unfold qq; positivity
  calc _ ≤ _ := h
    _ ≤ (n : ℝ) ^ (z.1 + z.2.1) * K * qq y * 1 :=
        mul_le_mul_of_nonneg_left h4 (by positivity)
    _ = _ := by ring

lemma f_meas (u : ℝ) : Measurable (fun y => f u y) := by unfold f; fun_prop

lemma fam_meas (n : ℕ) (z : ℕ × ℕ × ℕ) (u : ℝ) : Measurable (fun y => ww n y * fam n z u y) := by
  have h1 : Measurable (ww n) := by unfold ww; fun_prop
  have h2 : Measurable (AA n) := by unfold AA; fun_prop
  have h3 : Measurable ell := by unfold ell; fun_prop
  have h4 : Measurable (fun y => pif u y) := by
    unfold pif; exact (measurable_id.mul (h3.const_mul u).exp).div (f_meas u)
  have h5 : Measurable (fun y => E n u y) := by
    simp only [E_eq]; exact (h2.const_mul u).exp.mul ((f_meas u).pow_const _)
  unfold fam
  exact h1.mul ((((h2.pow_const _).mul (h3.pow_const _)).mul (h4.pow_const _)).mul h5)

noncomputable def PhiF (n : ℕ) (z : ℕ × ℕ × ℕ) (u : ℝ) : ℝ :=
  ∫ y in (0 : ℝ)..1, ww n y * fam n z u y

lemma ae_ne_one : ∀ᵐ t ∂(volume : Measure ℝ), t ≠ 1 := by
  rw [ae_iff]; simp

lemma ww_comb_cons (n : ℕ) (p : ℝ × (ℕ × ℕ × ℕ)) (L : List (ℝ × (ℕ × ℕ × ℕ))) (u : ℝ) :
    (fun y => ww n y * comb (fun z u => fam n z u y) (p :: L) u)
      = fun y => p.1 * (ww n y * fam n p.2 u y) + ww n y * comb (fun z u => fam n z u y) L u := by
  funext y; rw [comb_cons]; ring

lemma comb_meas (n : ℕ) (L : List (ℝ × (ℕ × ℕ × ℕ))) (u : ℝ) :
    Measurable (fun y => ww n y * comb (fun z u => fam n z u y) L u) := by
  induction L with
  | nil => simp [comb_nil]
  | cons p L ih => rw [ww_comb_cons]; exact ((fam_meas n p.2 u).const_mul _).add ih

lemma comb_bd {n : ℕ} (hn : 1 ≤ n) (L : List (ℝ × (ℕ × ℕ × ℕ))) : ∃ B : ℝ, 0 ≤ B ∧
    ∀ u y : ℝ, |u| ≤ dl n → 0 < y → y < 1 →
      |ww n y * comb (fun z u => fam n z u y) L u| ≤ B * qq y := by
  induction L with
  | nil => exact ⟨0, le_rfl, fun u y _ _ _ => by simp [comb_nil]⟩
  | cons p L ih =>
    obtain ⟨B1, hB1, h1⟩ := fam_bd hn p.2
    obtain ⟨B2, hB2, h2⟩ := ih
    refine ⟨|p.1| * B1 + B2, by positivity, fun u y hu hy0 hy1 => ?_⟩
    have e := congrFun (ww_comb_cons n p L u) y
    rw [e]
    have hq : 0 ≤ qq y := by unfold qq; positivity
    calc _ ≤ |p.1 * (ww n y * fam n p.2 u y)| + |ww n y * comb (fun z u => fam n z u y) L u| :=
          abs_add_le _ _
      _ ≤ |p.1| * (B1 * qq y) + B2 * qq y := by
          rw [abs_mul]
          exact add_le_add (mul_le_mul_of_nonneg_left (h1 u y hu hy0 hy1) (abs_nonneg _))
            (h2 u y hu hy0 hy1)
      _ = _ := by ring

lemma comb_intg {n : ℕ} (hn : 1 ≤ n) {u : ℝ} (hu : |u| ≤ dl n) (L : List (ℝ × (ℕ × ℕ × ℕ))) :
    IntervalIntegrable (fun y => ww n y * comb (fun z u => fam n z u y) L u) volume 0 1 := by
  obtain ⟨B, hB0, hB⟩ := comb_bd hn L
  refine IntervalIntegrable.mono_fun' (qq_int.const_mul B) (comb_meas n L u).aestronglyMeasurable ?_
  rw [uIoc_of_le zero_le_one, ← restrict_Ioo_eq_restrict_Ioc]
  exact ae_restrict_of_forall_mem measurableSet_Ioo fun t ht => by
    dsimp only; rw [Real.norm_eq_abs]; exact hB u t hu ht.1 ht.2

lemma famw_int {n : ℕ} (hn : 1 ≤ n) {u : ℝ} (hu : |u| ≤ dl n) (z : ℕ × ℕ × ℕ) :
    IntervalIntegrable (fun y => ww n y * fam n z u y) volume 0 1 := by
  have := comb_intg hn hu [(1, z)]
  simpa [comb] using this

lemma comb_integral {n : ℕ} (hn : 1 ≤ n) {u : ℝ} (hu : |u| ≤ dl n) (L : List (ℝ × (ℕ × ℕ × ℕ))) :
    ∫ y in (0 : ℝ)..1, ww n y * comb (fun z u => fam n z u y) L u = comb (PhiF n) L u := by
  induction L with
  | nil => simp [comb_nil]
  | cons p L ih =>
    rw [ww_comb_cons, intervalIntegral.integral_add ((famw_int hn hu p.2).const_mul _)
      (comb_intg hn hu L), intervalIntegral.integral_const_mul, ih, comb_cons]
    rfl

lemma PhiF_hasDerivAt {n : ℕ} (hn : 1 ≤ n) (z : ℕ × ℕ × ℕ) {u₀ : ℝ}
    (h : u₀ ∈ Ioo (-dl n) (dl n)) :
    HasDerivAt (PhiF n z) (comb (PhiF n) (step n z) u₀) u₀ := by
  obtain ⟨B, hB0, hB⟩ := comb_bd hn (step n z)
  have hS : Ioo (-dl n) (dl n) ∈ 𝓝 u₀ := Ioo_mem_nhds h.1 h.2
  have hu₀ : |u₀| ≤ dl n := abs_le.mpr ⟨h.1.le, h.2.le⟩
  have key := intervalIntegral.hasDerivAt_integral_of_dominated_loc_of_deriv_le
    (μ := volume) (a := 0) (b := 1)
    (F := fun u y => ww n y * fam n z u y)
    (F' := fun u y => ww n y * comb (fun z u => fam n z u y) (step n z) u)
    (bound := fun y => B * qq y) hS
    (Eventually.of_forall fun u => (fam_meas n z u).aestronglyMeasurable)
    (famw_int hn hu₀ z) (comb_meas n _ u₀).aestronglyMeasurable ?_
    (qq_int.const_mul B) ?_
  · rw [comb_integral hn hu₀] at key
    exact key.2
  · filter_upwards [ae_ne_one] with t ht1 htm x hx
    rw [uIoc_of_le zero_le_one] at htm
    have hx' : |x| ≤ dl n := abs_le.mpr ⟨hx.1.le, hx.2.le⟩
    rw [Real.norm_eq_abs]
    exact hB x t hx' htm.1 (lt_of_le_of_ne htm.2 ht1)
  · filter_upwards [ae_ne_one] with t ht1 htm x _
    rw [uIoc_of_le zero_le_one] at htm
    exact (fam_hasDerivAt n z htm.1 (lt_of_le_of_ne htm.2 ht1) x).const_mul (ww n t)

lemma comb_L0 (n : ℕ) (u : ℝ) :
    comb (PhiF n) [((1 : ℝ) / 2, (0, 0, 0))] u = 1 / 2 * ∫ y in (0 : ℝ)..1, ww n y * E n u y := by
  simp [comb, PhiF, fam]

lemma stmt5_core {n : ℕ} (hn : 1 ≤ n) (k : ℕ) :
    ContDiffAt ℝ k (comb (PhiF n) [((1 : ℝ) / 2, (0, 0, 0))]) 0 ∧
    iteratedDeriv k (comb (PhiF n) [((1 : ℝ) / 2, (0, 0, 0))]) 0 =
      1 / 2 * ∫ y in (0 : ℝ)..1, ww n y * iteratedDeriv k (fun u => E n u y) 0 := by
  have hd := dl_pos hn
  have hSo : IsOpen (Ioo (-dl n) (dl n)) := isOpen_Ioo
  have h0 : (0 : ℝ) ∈ Ioo (-dl n) (dl n) := ⟨by linarith, hd⟩
  have hG : ∀ z, ∀ u ∈ Ioo (-dl n) (dl n), HasDerivAt (PhiF n z) (comb (PhiF n) (step n z) u) u :=
    fun z u hu => PhiF_hasDerivAt hn z hu
  refine ⟨(comb_contDiffOn n hSo hG k _).contDiffAt (hSo.mem_nhds h0), ?_⟩
  rw [comb_iter_local n hSo hG _ k 0 h0,
    ← comb_integral hn (by rw [abs_zero]; exact hd.le), ← intervalIntegral.integral_const_mul]
  refine intervalIntegral.integral_congr_uIoo fun y hy => ?_
  rw [uIoo_of_le zero_le_one] at hy
  have hG' : ∀ z u, HasDerivAt ((fun z u => fam n z u y) z)
      (comb (fun z u => fam n z u y) (step n z) u) u :=
    fun z u => fam_hasDerivAt n z hy.1 hy.2 u
  have e1 := comb_iter_global n hG' k [((1 : ℝ) / 2, (0, 0, 0))]
  have e2 : comb (fun z u => fam n z u y) [((1 : ℝ) / 2, (0, 0, 0))] = fun u => 1 / 2 * E n u y := by
    funext u; simp [comb, fam]
  rw [← e1, e2, iteratedDeriv_const_mul_field]
  ring

noncomputable def LL (n k : ℕ) : List (ℝ × (ℕ × ℕ × ℕ)) := (Dl n)^[k] [((1 : ℝ) / 2, (0, 0, 0))]

noncomputable def wsum (L : List (ℝ × (ℕ × ℕ × ℕ))) : ℝ := (L.map (fun p => |p.1|)).sum

lemma wsum_nonneg (L : List (ℝ × (ℕ × ℕ × ℕ))) : 0 ≤ wsum L := by
  unfold wsum; exact List.sum_nonneg (fun x hx => by
    obtain ⟨p, _, rfl⟩ := List.mem_map.mp hx; exact abs_nonneg _)

lemma wsum_Dl {n : ℕ} (hn : 1 ≤ n) (K : ℕ) (L : List (ℝ × (ℕ × ℕ × ℕ)))
    (hL : ∀ p ∈ L, p.2.2.2 ≤ K) : wsum (Dl n L) ≤ (2 * n + 2 * K) * wsum L := by
  induction L with
  | nil => simp [wsum, Dl]
  | cons p L ih =>
    have hp := hL p (by simp)
    have ih' := ih (fun q hq => hL q (by simp [hq]))
    have hn' : (1 : ℝ) ≤ n := by exact_mod_cast hn
    have hj : ((p.2.2.2 : ℕ) : ℝ) ≤ K := by exact_mod_cast hp
    have h0 : (0 : ℝ) ≤ p.2.2.2 := Nat.cast_nonneg _
    have hD : wsum (Dl n (p :: L)) = |p.1| * (1 + |(1 - 2 * (n : ℝ)) - p.2.2.2| + |((p.2.2.2 : ℕ) : ℝ)|)
        + wsum (Dl n L) := by
      simp only [wsum, Dl, List.flatMap_cons, List.map_append, List.sum_append, step,
        List.map_cons, List.map_nil, List.sum_cons, List.sum_nil, abs_mul, abs_one]
      ring
    have hW : wsum (p :: L) = |p.1| + wsum L := by simp [wsum]
    rw [hD, hW, abs_of_nonpos (show (1 - 2 * (n : ℝ)) - p.2.2.2 ≤ 0 by linarith),
      abs_of_nonneg h0]
    have : |p.1| * (1 + -(1 - 2 * (n : ℝ) - p.2.2.2) + p.2.2.2) ≤ |p.1| * (2 * n + 2 * K) :=
      mul_le_mul_of_nonneg_left (by linarith) (abs_nonneg _)
    have hK0 : (0 : ℝ) ≤ 2 * n + 2 * K := by positivity
    nlinarith [abs_nonneg p.1]

lemma LL_succ (n k : ℕ) : LL n (k + 1) = Dl n (LL n k) := by
  unfold LL; rw [Function.iterate_succ_apply']

lemma LL_struct {n : ℕ} (hn : 1 ≤ n) (k : ℕ) :
    (∀ p ∈ LL n k, p.2.1 + p.2.2.1 = k ∧ p.2.2.2 ≤ k) ∧ wsum (LL n k) ≤ (2 * n + 2 * k) ^ k := by
  induction k with
  | zero => refine ⟨fun p hp => ?_, ?_⟩
            · simp [LL] at hp; subst hp; simp
            · simp [LL, wsum]; norm_num
  | succ k ih =>
    obtain ⟨h1, h2⟩ := ih
    refine ⟨fun p hp => ?_, ?_⟩
    · rw [LL_succ] at hp
      simp only [Dl, List.mem_flatMap, List.mem_map] at hp
      obtain ⟨q, hq, r, hr, rfl⟩ := hp
      obtain ⟨e1, e2⟩ := h1 q hq
      simp only [step, List.mem_cons, List.not_mem_nil, or_false] at hr
      rcases hr with rfl | rfl | rfl <;> simp <;> omega
    · rw [LL_succ]
      have hw := wsum_Dl hn k (LL n k) (fun p hp => (h1 p hp).2)
      have hK0 : (0 : ℝ) ≤ 2 * n + 2 * k := by positivity
      have hmono : (2 * (n : ℝ) + 2 * k) ^ k ≤ (2 * n + 2 * ((k + 1 : ℕ) : ℝ)) ^ k :=
        pow_le_pow_left₀ hK0 (by push_cast; linarith) k
      calc wsum (Dl n (LL n k)) ≤ (2 * n + 2 * k) * wsum (LL n k) := hw
        _ ≤ (2 * n + 2 * ((k + 1 : ℕ) : ℝ)) * (2 * n + 2 * ((k + 1 : ℕ) : ℝ)) ^ k := by
          apply mul_le_mul (by push_cast; linarith) (h2.trans hmono) (wsum_nonneg _)
          push_cast; positivity
        _ = _ := by rw [pow_succ]; ring

lemma comb_pt_bd (N : ℕ) : ∃ K : ℝ, 0 ≤ K ∧ ∀ n : ℕ, 1 ≤ n → ∀ L : List (ℝ × (ℕ × ℕ × ℕ)),
    (∀ p ∈ L, p.2.1 + p.2.2.1 = N) → ∀ d u y : ℝ, 0 ≤ d → d ≤ 1 / 12 → |u| ≤ d → 0 < y → y < 1 →
      |ww n y * comb (fun z u => fam n z u y) L u| ≤
        wsum L * ((n : ℝ) ^ N * K * qq y * (1 / 4 : ℝ) ^ ((n : ℝ) - 1 / 2 - (3 * (n : ℝ) - 1) * d)) := by
  obtain ⟨K, hK0, hK⟩ := fam_bdG N
  refine ⟨K, hK0, fun n hn L hL d u y hd0 hd1 hu hy0 hy1 => ?_⟩
  induction L with
  | nil => simp [comb_nil, wsum]
  | cons p L ih =>
    have ih' := ih (fun q hq => hL q (by simp [hq]))
    have hp := hK n hn p.2 (hL p (by simp)) d u y hd0 hd1 hu hy0 hy1
    have e := congrFun (ww_comb_cons n p L u) y
    rw [e]
    have hW : wsum (p :: L) = |p.1| + wsum L := by simp [wsum]
    rw [hW]
    calc _ ≤ |p.1 * (ww n y * fam n p.2 u y)| + |ww n y * comb (fun z u => fam n z u y) L u| :=
          abs_add_le _ _
      _ ≤ _ := by
          rw [abs_mul, add_mul]
          exact add_le_add (mul_le_mul_of_nonneg_left hp (abs_nonneg _)) ih'

noncomputable def Qint : ℝ := ∫ y in (0 : ℝ)..1, qq y

lemma comb_decay (k : ℕ) : ∃ K : ℝ, 0 ≤ K ∧ ∀ n : ℕ, 1 ≤ n → ∀ d u : ℝ, 0 ≤ d → d ≤ 1 / 12 →
    |u| ≤ d → |comb (PhiF n) (LL n k) u| ≤
      (2 * (n : ℝ) + 2 * k) ^ k * ((n : ℝ) ^ k * K * Qint *
        (1 / 4 : ℝ) ^ ((n : ℝ) - 1 / 2 - (3 * (n : ℝ) - 1) * d)) := by
  obtain ⟨K, hK0, hK⟩ := comb_pt_bd k
  refine ⟨K, hK0, fun n hn d u hd0 hd1 hu => ?_⟩
  obtain ⟨hs1, hs2⟩ := LL_struct hn k
  have hu' : |u| ≤ dl n := by unfold dl; linarith
  set e : ℝ := (n : ℝ) - 1 / 2 - (3 * (n : ℝ) - 1) * d
  rw [← comb_integral hn hu']
  have hQ : 0 ≤ Qint := by
    unfold Qint
    exact intervalIntegral.integral_nonneg zero_le_one fun y hy => by
      unfold qq; exact Real.rpow_nonneg (mul_nonneg hy.1 (by linarith [hy.2])) _
  have hb := intervalIntegral.norm_integral_le_of_norm_le (μ := volume)
    (f := fun y => ww n y * comb (fun z u => fam n z u y) (LL n k) u)
    (g := fun y => wsum (LL n k) * ((n : ℝ) ^ k * K * (1 / 4 : ℝ) ^ e) * qq y) zero_le_one ?_
    (qq_int.const_mul _)
  · rw [intervalIntegral.integral_const_mul, Real.norm_eq_abs] at hb
    refine hb.trans ?_
    have h4 : 0 ≤ (1 / 4 : ℝ) ^ e := by positivity
    calc wsum (LL n k) * ((n : ℝ) ^ k * K * (1 / 4 : ℝ) ^ e) * (∫ y in (0 : ℝ)..1, qq y)
        = wsum (LL n k) * ((n : ℝ) ^ k * K * Qint * (1 / 4 : ℝ) ^ e) := by rw [Qint]; ring
      _ ≤ _ := mul_le_mul_of_nonneg_right hs2
          (mul_nonneg (mul_nonneg (mul_nonneg (by positivity) hK0) hQ) h4)
  · filter_upwards [ae_ne_one] with t ht1 htm
    have := hK n hn (LL n k) (fun p hp => (hs1 p hp).1) d u t hd0 hd1 hu htm.1
      (lt_of_le_of_ne htm.2 ht1)
    rw [Real.norm_eq_abs]
    calc _ ≤ _ := this
      _ = _ := by ring

noncomputable def rho4 (δ : ℝ) : ℝ := 4 * (1 / 4 : ℝ) ^ (1 - 3 * δ)

lemma an_le (m : ℕ) (hm : 1 ≤ m) : 0 ≤ an m ∧ an m ≤ 4 ^ m := by
  have hm' : (1 : ℝ) ≤ m := by exact_mod_cast hm
  have h1 : (1 : ℝ) ≤ 2 * m - 1 := by linarith
  have hc : ((Nat.choose (2 * m) m : ℕ) : ℝ) ≤ 4 ^ m := by
    have := Nat.choose_le_two_pow (2 * m) m
    have h : ((Nat.choose (2 * m) m : ℕ) : ℝ) ≤ ((2 ^ (2 * m) : ℕ) : ℝ) := by exact_mod_cast this
    rw [Nat.cast_pow, pow_mul] at h; norm_num at h; exact h
  unfold an
  refine ⟨div_nonneg (Nat.cast_nonneg _) (by linarith), ?_⟩
  rw [div_le_iff₀ (by linarith)]
  nlinarith [Nat.cast_nonneg (α := ℝ) (Nat.choose (2 * m) m), pow_pos (by norm_num : (0:ℝ) < 4) m]

noncomputable def trm (c : ℝ) (k n : ℕ) (u : ℝ) : ℝ :=
  an (n + 1) * c ^ (2 * (n + 1)) * comb (PhiF (n + 1)) (LL (n + 1) k) u

lemma trm_bd (k : ℕ) : ∃ C : ℝ, 0 ≤ C ∧ ∀ c δ : ℝ, 0 < c → 0 ≤ δ → δ ≤ 1 / 12 → ∀ (n : ℕ) (u : ℝ),
    |u| ≤ δ → |trm c k n u| ≤ C * ((n : ℝ) + 1) ^ (2 * k) * (c ^ 2 * rho4 δ) ^ (n + 1) := by
  obtain ⟨K, hK0, hK⟩ := comb_decay k
  have hQ : 0 ≤ Qint := by
    unfold Qint
    exact intervalIntegral.integral_nonneg zero_le_one fun y hy => by
      unfold qq; exact Real.rpow_nonneg (mul_nonneg hy.1 (by linarith [hy.2])) _
  refine ⟨K * Qint * 4 * (2 + 2 * k) ^ k, by positivity, fun c δ hc hδ0 hδ1 n u hu => ?_⟩
  set m : ℕ := n + 1 with hm
  have hm1 : 1 ≤ m := by omega
  have hm' : (1 : ℝ) ≤ m := by exact_mod_cast hm1
  have hmc : (m : ℝ) = (n : ℝ) + 1 := by rw [hm]; push_cast; ring
  have hb := hK m hm1 δ u hδ0 hδ1 hu
  obtain ⟨ha0, ha1⟩ := an_le m hm1
  have hpoly : (2 * (m : ℝ) + 2 * k) ^ k * (m : ℝ) ^ k ≤ (2 + 2 * k) ^ k * (m : ℝ) ^ (2 * k) := by
    have : (2 * (m : ℝ) + 2 * k) ≤ (2 + 2 * k) * m := by nlinarith
    calc (2 * (m : ℝ) + 2 * k) ^ k * (m : ℝ) ^ k ≤ ((2 + 2 * k) * m) ^ k * (m : ℝ) ^ k :=
          mul_le_mul_of_nonneg_right (pow_le_pow_left₀ (by positivity) this k) (by positivity)
      _ = _ := by rw [mul_pow, mul_assoc, ← pow_add, two_mul k]
  -- the exponential factor
  have hexp : (1 / 4 : ℝ) ^ ((m : ℝ) - 1 / 2 - (3 * (m : ℝ) - 1) * δ)
      ≤ 4 * ((1 / 4 : ℝ) ^ (1 - 3 * δ)) ^ m := by
    have e1 : (m : ℝ) - 1 / 2 - (3 * (m : ℝ) - 1) * δ = (1 - 3 * δ) * m + (δ - 1 / 2) := by ring
    rw [e1, Real.rpow_add (by norm_num), Real.rpow_mul (by norm_num), Real.rpow_natCast]
    have : (1 / 4 : ℝ) ^ (δ - 1 / 2) ≤ 4 := by
      have h := Real.rpow_le_rpow_of_exponent_ge (show (0:ℝ) < 1 / 4 by norm_num)
        (show (1 / 4 : ℝ) ≤ 1 by norm_num) (show (-1 : ℝ) ≤ δ - 1 / 2 by linarith)
      rw [Real.rpow_neg_one] at h; norm_num at h ⊢; linarith
    have h0 : 0 ≤ ((1 / 4 : ℝ) ^ (1 - 3 * δ)) ^ m := by positivity
    nlinarith
  have hcomb0 : 0 ≤ |comb (PhiF m) (LL m k) u| := abs_nonneg _
  unfold trm
  rw [← hm, abs_mul, abs_mul, abs_of_nonneg ha0, abs_of_nonneg (by positivity : (0:ℝ) ≤ c ^ (2 * m))]
  have hX0 : 0 ≤ (1 / 4 : ℝ) ^ ((m : ℝ) - 1 / 2 - (3 * (m : ℝ) - 1) * δ) := by positivity
  calc an m * c ^ (2 * m) * |comb (PhiF m) (LL m k) u|
      ≤ 4 ^ m * c ^ (2 * m) * ((2 * (m : ℝ) + 2 * k) ^ k * ((m : ℝ) ^ k * K * Qint *
          (1 / 4 : ℝ) ^ ((m : ℝ) - 1 / 2 - (3 * (m : ℝ) - 1) * δ))) :=
        mul_le_mul (mul_le_mul_of_nonneg_right ha1 (by positivity)) hb hcomb0 (by positivity)
    _ = ((2 * (m : ℝ) + 2 * k) ^ k * (m : ℝ) ^ k) * (K * Qint) * (4 ^ m * c ^ (2 * m)) *
          (1 / 4 : ℝ) ^ ((m : ℝ) - 1 / 2 - (3 * (m : ℝ) - 1) * δ) := by ring
    _ ≤ ((2 + 2 * k) ^ k * (m : ℝ) ^ (2 * k)) * (K * Qint) * (4 ^ m * c ^ (2 * m)) *
          (4 * ((1 / 4 : ℝ) ^ (1 - 3 * δ)) ^ m) := by
        gcongr
    _ = K * Qint * 4 * (2 + 2 * k) ^ k * ((n : ℝ) + 1) ^ (2 * k) * (c ^ 2 * rho4 δ) ^ m := by
        rw [rho4, ← hmc, mul_pow, mul_pow, ← pow_mul]; ring

lemma rho4_eq (δ : ℝ) : rho4 δ = Real.exp (3 * δ * Real.log 4) := by
  unfold rho4
  rw [Real.rpow_def_of_pos (by norm_num), one_div, Real.log_inv,
    show (4 : ℝ) = Real.exp (Real.log 4) from (Real.exp_log (by norm_num)).symm, Real.log_exp,
    ← Real.exp_add]
  congr 1; ring

lemma exists_delta {c : ℝ} (hc0 : 0 < c) (hc1 : c < 1) :
    ∃ δ : ℝ, 0 < δ ∧ δ ≤ 1 / 12 ∧ c ^ 2 * rho4 δ < 1 := by
  have hl4 : 0 < Real.log 4 := Real.log_pos (by norm_num)
  have hlc : 0 < Real.log (1 / c) := Real.log_pos (by rw [lt_div_iff₀ hc0]; linarith)
  refine ⟨min (1 / 12) (Real.log (1 / c) / (3 * Real.log 4)), lt_min (by norm_num) (by positivity),
    min_le_left _ _, ?_⟩
  rw [rho4_eq]
  have h1 : 3 * min (1 / 12) (Real.log (1 / c) / (3 * Real.log 4)) * Real.log 4 ≤ Real.log (1 / c) := by
    have := min_le_right (1 / 12 : ℝ) (Real.log (1 / c) / (3 * Real.log 4))
    have e : 3 * (Real.log (1 / c) / (3 * Real.log 4)) * Real.log 4 = Real.log (1 / c) := by
      field_simp
    nlinarith
  have h2 : Real.exp (3 * min (1 / 12) (Real.log (1 / c) / (3 * Real.log 4)) * Real.log 4) ≤ 1 / c := by
    calc _ ≤ Real.exp (Real.log (1 / c)) := Real.exp_le_exp.mpr h1
      _ = 1 / c := Real.exp_log (by positivity)
  have h3 : c ^ 2 * (1 / c) = c := by field_simp
  have h4 : 0 ≤ c ^ 2 := by positivity
  nlinarith [mul_le_mul_of_nonneg_left h2 h4]

lemma maj_summable (k : ℕ) {r : ℝ} (hr0 : 0 ≤ r) (hr : r < 1) (C : ℝ) :
    Summable (fun n : ℕ => C * ((n : ℝ) + 1) ^ (2 * k) * r ^ (n + 1)) := by
  have h := summable_pow_mul_geometric_of_norm_lt_one (R := ℝ) (2 * k)
    (by rw [Real.norm_eq_abs, abs_of_nonneg hr0]; exact hr)
  have h1 := (summable_nat_add_iff 1).mpr h
  have h2 := h1.mul_left C
  refine h2.congr fun n => ?_
  push_cast; ring

noncomputable def Hs (c : ℝ) (k : ℕ) (u : ℝ) : ℝ := ∑' n, trm c k n u

lemma trm_hasDerivAt (c : ℝ) (k n : ℕ) {u : ℝ} (hu : u ∈ Ioo (-(1 / 12 : ℝ)) (1 / 12)) :
    HasDerivAt (trm c k n) (trm c (k + 1) n u) u := by
  unfold trm
  rw [LL_succ]
  have hu' : u ∈ Ioo (-dl (n + 1)) (dl (n + 1)) := by unfold dl; exact hu
  exact (comb_hasDerivAt (n + 1) (fun z => PhiF_hasDerivAt (by omega) z hu') _).const_mul _

lemma Hs_hasDerivAt {c δ : ℝ} (hc0 : 0 < c) (hδ0 : 0 < δ) (hδ1 : δ ≤ 1 / 12)
    (hr : c ^ 2 * rho4 δ < 1) (k : ℕ) {u : ℝ} (hu : u ∈ Ioo (-δ) δ) :
    HasDerivAt (Hs c k) (Hs c (k + 1) u) u := by
  obtain ⟨C1, hC1, hb1⟩ := trm_bd (k + 1)
  obtain ⟨C0, hC0, hb0⟩ := trm_bd k
  have hr0 : 0 ≤ c ^ 2 * rho4 δ := by unfold rho4; positivity
  have hsub : ∀ y ∈ Ioo (-δ) δ, y ∈ Ioo (-(1 / 12 : ℝ)) (1 / 12) := fun y hy =>
    ⟨by linarith [hy.1], by linarith [hy.2]⟩
  have habs : ∀ y ∈ Ioo (-δ) δ, |y| ≤ δ := fun y hy => abs_le.mpr ⟨hy.1.le, hy.2.le⟩
  unfold Hs
  exact hasDerivAt_tsum_of_isPreconnected (maj_summable (k + 1) hr0 hr C1) isOpen_Ioo
    isPreconnected_Ioo (fun n y hy => trm_hasDerivAt c k n (hsub y hy))
    (fun n y hy => by rw [Real.norm_eq_abs]; exact hb1 c δ hc0 hδ0.le hδ1 n y (habs y hy))
    (show (0 : ℝ) ∈ Ioo (-δ) δ from ⟨by linarith, hδ0⟩)
    ((maj_summable k hr0 hr C0).of_norm_bounded fun n => by
      rw [Real.norm_eq_abs]; exact hb0 c δ hc0 hδ0.le hδ1 n 0 (by rw [abs_zero]; exact hδ0.le))
    hu

lemma Hs_summable {c δ : ℝ} (hc0 : 0 < c) (hδ0 : 0 < δ) (hδ1 : δ ≤ 1 / 12)
    (hr : c ^ 2 * rho4 δ < 1) (k : ℕ) {u : ℝ} (hu : |u| ≤ δ) : Summable (fun n => trm c k n u) := by
  obtain ⟨C0, hC0, hb0⟩ := trm_bd k
  have hr0 : 0 ≤ c ^ 2 * rho4 δ := by unfold rho4; positivity
  exact (maj_summable k hr0 hr C0).of_norm_bounded fun n => by
    rw [Real.norm_eq_abs]; exact hb0 c δ hc0 hδ0.le hδ1 n u hu

lemma Hs_iter {c δ : ℝ} (hc0 : 0 < c) (hδ0 : 0 < δ) (hδ1 : δ ≤ 1 / 12)
    (hr : c ^ 2 * rho4 δ < 1) :
    ∀ k, ∀ u ∈ Ioo (-δ) δ, iteratedDeriv k (Hs c 0) u = Hs c k u := by
  intro k
  induction k with
  | zero => intro u _; simp
  | succ k ih =>
    intro u hu
    rw [iteratedDeriv_succ]
    have hev : iteratedDeriv k (Hs c 0) =ᶠ[𝓝 u] Hs c k :=
      Filter.eventually_of_mem (isOpen_Ioo.mem_nhds hu) fun v hv => ih v hv
    rw [hev.deriv_eq, (Hs_hasDerivAt hc0 hδ0 hδ1 hr k hu).deriv]

lemma Hs_contDiffOn {c δ : ℝ} (hc0 : 0 < c) (hδ0 : 0 < δ) (hδ1 : δ ≤ 1 / 12)
    (hr : c ^ 2 * rho4 δ < 1) : ∀ m : ℕ, ∀ k, ContDiffOn ℝ m (Hs c k) (Ioo (-δ) δ) := by
  intro m
  induction m with
  | zero =>
    intro k
    rw [Nat.cast_zero, contDiffOn_zero]
    intro u hu
    exact (Hs_hasDerivAt hc0 hδ0 hδ1 hr k hu).continuousAt.continuousWithinAt
  | succ m ih =>
    intro k
    rw [show ((m + 1 : ℕ) : WithTop ℕ∞) = (m : WithTop ℕ∞) + 1 by push_cast; rfl,
      contDiffOn_succ_iff_deriv_of_isOpen isOpen_Ioo]
    refine ⟨fun u hu => (Hs_hasDerivAt hc0 hδ0 hδ1 hr k hu).differentiableAt.differentiableWithinAt,
      fun h => absurd h (by simp), ?_⟩
    refine (ih (k + 1)).congr fun u hu => ?_
    exact (Hs_hasDerivAt hc0 hδ0 hδ1 hr k hu).deriv

end V05

theorem V05.proof_5 : V05.stmt_5 := by
  intro n hn k
  obtain ⟨h1, h2⟩ := V05.stmt5_core hn k
  have hev : V05.J n =ᶠ[nhds 0] V05.comb (V05.PhiF n) [((1 : ℝ) / 2, (0, 0, 0))] := by
    filter_upwards [Ioo_mem_nhds (show (-1 : ℝ) < 0 by norm_num) (show (0 : ℝ) < 1 by norm_num)]
      with u hu
    rw [V05.comb_L0, V05.proof_3 n hn u (abs_lt.mpr hu)]
    rfl
  refine ⟨h1.congr_of_eventuallyEq hev, ?_⟩
  rw [hev.iteratedDeriv_eq, h2, mul_div_assoc, ← intervalIntegral.integral_div]
  congr 1
  refine intervalIntegral.integral_congr fun y _ => ?_
  simp only [V05.ww]
  ring

theorem V05.proof_4 : V05.stmt_4 := by
  intro α hα0 hα1
  obtain ⟨hc0, hc1, -, -, -⟩ := D04.cs_facts hα0 hα1
  obtain ⟨δ, hδ0, hδ1, hr⟩ := V05.exists_delta hc0 hc1
  have hJ : ∀ (n : ℕ) (u : ℝ), |u| < 1 →
      V05.J (n + 1) u = V05.comb (V05.PhiF (n + 1)) (V05.LL (n + 1) 0) u := by
    intro n u hu
    rw [show V05.LL (n + 1) 0 = [((1 : ℝ) / 2, (0, 0, 0))] from rfl, V05.comb_L0,
      V05.proof_3 (n + 1) (by omega) u hu]
    rfl
  have hM : ∀ u ∈ Set.Ioo (-δ) δ, D04.Mint u α = -(V05.Hs (D04.c α) 0 u) := by
    intro u hu
    have hu1 : |u| < 1 := by rw [abs_lt]; constructor <;> linarith [hu.1, hu.2]
    have h2 := (V05.proof_2 α hα0 hα1 u hu1).2.2
    have e : (fun n : ℕ => V05.an (n + 1) * D04.c α ^ (2 * (n + 1)) * V05.J (n + 1) u)
        = fun n => V05.trm (D04.c α) 0 n u := by
      funext n; unfold V05.trm; rw [hJ n u hu1]
    rw [e] at h2
    unfold V05.Hs; rw [h2.tsum_eq]; ring
  have h0 : (0 : ℝ) ∈ Set.Ioo (-δ) δ := ⟨by linarith, hδ0⟩
  refine ⟨δ, hδ0, hδ1, fun k => ?_, fun k => ?_⟩
  · exact ((V05.Hs_contDiffOn hc0 hδ0 hδ1 hr k 0).neg).congr (fun u hu => hM u hu)
  · have hev : (fun u => D04.Mint u α) =ᶠ[nhds 0] fun u => (-1) * V05.Hs (D04.c α) 0 u :=
      Filter.eventually_of_mem (isOpen_Ioo.mem_nhds h0) fun u hu => by
        show D04.Mint u α = _; rw [hM u hu]; ring
    have hMk : V05.Mk k α = -(V05.Hs (D04.c α) k 0) / (k.factorial : ℝ) := by
      unfold V05.Mk
      rw [hev.iteratedDeriv_eq, iteratedDeriv_const_mul_field,
        V05.Hs_iter hc0 hδ0 hδ1 hr k 0 h0]
      ring
    have hJk : ∀ n : ℕ, iteratedDeriv k (V05.J (n + 1)) 0
        = V05.comb (V05.PhiF (n + 1)) (V05.LL (n + 1) k) 0 := by
      intro n
      have hevJ : V05.J (n + 1) =ᶠ[nhds 0] V05.comb (V05.PhiF (n + 1)) (V05.LL (n + 1) 0) :=
        Filter.eventually_of_mem (Ioo_mem_nhds (show (-1 : ℝ) < 0 by norm_num)
          (show (0 : ℝ) < 1 by norm_num)) fun u hu => hJ n u (abs_lt.mpr hu)
      rw [hevJ.iteratedDeriv_eq]
      have h0' : (0 : ℝ) ∈ Set.Ioo (-V05.dl (n + 1)) (V05.dl (n + 1)) := by
        unfold V05.dl; constructor <;> norm_num
      exact V05.comb_iter_local (n + 1) isOpen_Ioo
        (fun z u hu => V05.PhiF_hasDerivAt (by omega) z hu) _ k 0 h0'
    have hs := (V05.Hs_summable hc0 hδ0 hδ1 hr k (u := 0)
      (by rw [abs_zero]; exact hδ0.le)).hasSum
    rw [hMk, neg_div, neg_neg]
    unfold V05.Hs
    convert hs.div_const (k.factorial : ℝ) using 1
    funext n
    rw [hJk n]
    unfold V05.trm
    ring

namespace V05
open MeasureTheory Set Filter Real Interval Topology

lemma ae_Ioo01 : ∀ᵐ t ∂(volume.restrict (Ι (0:ℝ) 1)), t ∈ Ioo (0:ℝ) 1 := by
  rw [uIoc_of_le zero_le_one, ← Measure.restrict_congr_set Ioo_ae_eq_Ioc]
  exact ae_restrict_mem measurableSet_Ioo

lemma ii_tt (c : ℝ) (hc : -1 < c) :
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

lemma abs_log_le {x θ : ℝ} (h0 : 0 < x) (h1 : x ≤ 1) (hθ : 0 < θ) :
    |log x| ≤ x ^ (-θ) / θ := by
  have := abs_log_mul_self_rpow_lt x θ h0 h1 hθ
  have hp := rpow_pos_of_pos h0 θ
  rw [abs_mul, abs_of_pos hp] at this
  rw [rpow_neg h0.le, show (x ^ θ)⁻¹ / θ = (1 / θ) / x ^ θ by field_simp, le_div_iff₀ hp]
  linarith

lemma logsum_le {t θ : ℝ} (ht : t ∈ Ioo (0:ℝ) 1) (hθ : 0 < θ) :
    |log t| + |log (1 - t)| ≤ 2 / θ * (t * (1 - t)) ^ (-θ) := by
  obtain ⟨h0, h1⟩ := ht
  have hu : 0 < 1 - t := by linarith
  have hp : 0 < t * (1 - t) := mul_pos h0 hu
  have e1 := abs_log_le h0 h1.le hθ
  have e2 := abs_log_le hu (by linarith) hθ
  have f1 : t ^ (-θ) ≤ (t * (1 - t)) ^ (-θ) :=
    rpow_le_rpow_of_nonpos hp (by nlinarith) (by linarith)
  have f2 : (1 - t) ^ (-θ) ≤ (t * (1 - t)) ^ (-θ) :=
    rpow_le_rpow_of_nonpos hp (by nlinarith) (by linarith)
  have : t ^ (-θ) / θ + (1 - t) ^ (-θ) / θ ≤ 2 / θ * (t * (1 - t)) ^ (-θ) := by
    rw [← add_div, div_le_iff₀ hθ]
    have : 2 / θ * (t * (1 - t)) ^ (-θ) * θ = 2 * (t * (1 - t)) ^ (-θ) := by field_simp
    rw [this]; linarith
  linarith

/-- the log-sum `s(t) = |log t| + |log(1-t)|` -/
noncomputable def ls (t : ℝ) : ℝ := |log t| + |log (1 - t)|

lemma ls_nonneg (t : ℝ) : 0 ≤ ls t := by unfold ls; positivity

lemma ls_pow_le {t : ℝ} (ht : t ∈ Ioo (0:ℝ) 1) (m : ℕ) :
    ls t ^ m ≤ (4 * (m + 1)) ^ m * (t * (1 - t)) ^ (-(1/2:ℝ)) := by
  have hp : 0 < t * (1 - t) := mul_pos ht.1 (by linarith [ht.2])
  have hθ : (0:ℝ) < 1 / (2 * (m + 1)) := by positivity
  have h := logsum_le ht hθ
  have h' : ls t ^ m ≤ (2 / (1 / (2 * (m + 1))) * (t * (1 - t)) ^ (-(1 / (2 * ((m:ℝ) + 1))))) ^ m :=
    pow_le_pow_left₀ (ls_nonneg t) h m
  rw [mul_pow, ← rpow_mul_natCast hp.le] at h'
  have hle1 : t * (1 - t) ≤ 1 := by nlinarith [ht.1, ht.2]
  have : (t * (1 - t)) ^ (-(1 / (2 * ((m:ℝ) + 1))) * (m:ℝ)) ≤ (t * (1 - t)) ^ (-(1/2:ℝ)) := by
    apply rpow_le_rpow_of_exponent_ge hp hle1
    rw [neg_mul, neg_le_neg_iff, div_mul_eq_mul_div, one_mul, div_le_iff₀ (by positivity)]
    linarith
  have e : (2 / (1 / (2 * ((m:ℝ) + 1)))) = 4 * (m + 1) := by field_simp; ring
  rw [e] at h'
  calc ls t ^ m ≤ _ := h'
    _ ≤ _ := by gcongr

/-- Integrability of anything bounded by `K * s(t)^m` -/
lemma ii_of_ls {f : ℝ → ℝ} (hf : AEStronglyMeasurable f (volume.restrict (Ι (0:ℝ) 1)))
    (K : ℝ) (m : ℕ) (h : ∀ t ∈ Ioo (0:ℝ) 1, |f t| ≤ K * ls t ^ m) :
    IntervalIntegrable f volume 0 1 := by
  refine IntervalIntegrable.mono_fun'
    (((ii_tt (-(1/2)) (by norm_num)).const_mul (|K| * (4 * (m + 1)) ^ m))) hf ?_
  filter_upwards [ae_Ioo01] with t ht
  have h1 := h t ht
  have h2 := ls_pow_le ht m
  have hK : K ≤ |K| := le_abs_self K
  have hl := pow_nonneg (ls_nonneg t) m
  simp only [Real.norm_eq_abs]
  calc |f t| ≤ K * ls t ^ m := h1
    _ ≤ |K| * ls t ^ m := by gcongr
    _ ≤ |K| * ((4 * (m + 1)) ^ m * (t * (1 - t)) ^ (-(1/2:ℝ))) := by gcongr
    _ = _ := by ring

lemma meas_mono (i j : ℕ) :
    AEStronglyMeasurable (fun t : ℝ => log t ^ i * log (1 - t) ^ j) (volume.restrict (Ι (0:ℝ) 1)) :=
  (by fun_prop : Measurable (fun t : ℝ => log t ^ i * log (1 - t) ^ j)).aestronglyMeasurable

lemma ii_mono (i j : ℕ) : IntervalIntegrable (fun t : ℝ => log t ^ i * log (1 - t) ^ j) volume 0 1 := by
  refine ii_of_ls (meas_mono i j) 1 (i + j) fun t _ => ?_
  rw [abs_mul, abs_pow, abs_pow, one_mul, pow_add]
  have a1 : |log t| ≤ ls t := by unfold ls; linarith [abs_nonneg (log (1 - t))]
  have a2 : |log (1 - t)| ≤ ls t := by unfold ls; linarith [abs_nonneg (log t)]
  have := ls_nonneg t
  gcongr


lemma ae_Ioo01' : ∀ᵐ t ∂(volume : Measure ℝ), t ∈ Ι (0:ℝ) 1 → t ∈ Ioo (0:ℝ) 1 :=
  (ae_restrict_iff' measurableSet_uIoc).mp ae_Ioo01

lemma uIoo01 : uIoo (0:ℝ) 1 = Ioo 0 1 := by simp [uIoo]

lemma t_log_pow_tendsto (k : ℕ) (hk : k ≠ 0) :
    Tendsto (fun t : ℝ => t * log t ^ k) (𝓝[>] 0) (𝓝 0) := by
  have hr : (0:ℝ) < (k:ℝ)⁻¹ := by positivity
  have h := (tendsto_log_mul_rpow_nhdsGT_zero hr).pow k
  rw [zero_pow hk] at h
  refine h.congr' ?_
  filter_upwards [self_mem_nhdsWithin] with t ht
  rw [mul_pow, rpow_inv_natCast_pow (le_of_lt ht) hk, mul_comm]

lemma G_cont (n k : ℕ) (hk : k ≠ 0) :
    ContinuousOn (fun t : ℝ => t ^ (n+1) * log t ^ k) (Icc 0 1) := by
  intro x hx
  rcases eq_or_lt_of_le hx.1 with h | h
  · subst h
    apply ContinuousWithinAt.mono _ Icc_subset_Ici_self
    rw [← continuousWithinAt_Ioi_iff_Ici, ContinuousWithinAt]
    have := (((continuous_pow n).tendsto 0).mono_left nhdsWithin_le_nhds).mul
      (t_log_pow_tendsto k hk)
    rw [mul_zero] at this
    simp only [ne_eq, Nat.add_eq_zero_iff, one_ne_zero, and_false, not_false_eq_true, zero_pow,
      zero_mul]
    refine this.congr fun t => ?_
    ring
  · exact ((continuous_pow _).continuousAt.mul ((continuousAt_log h.ne').pow k)).continuousWithinAt

/-- `Jl n k = ∫_0^1 t^n log^k t` -/
noncomputable def Jl (n k : ℕ) : ℝ := ∫ t in (0:ℝ)..1, t ^ n * log t ^ k

lemma ii_J (n k : ℕ) : IntervalIntegrable (fun t : ℝ => t ^ n * log t ^ k) volume 0 1 := by
  refine ii_of_ls (by fun_prop : Measurable (fun t : ℝ => t ^ n * log t ^ k)).aestronglyMeasurable
    1 k fun t ht => ?_
  rw [abs_mul, abs_pow, abs_pow, one_mul]
  have h1 : |t| ^ n ≤ 1 := pow_le_one₀ (abs_nonneg t) (by rw [abs_of_pos ht.1]; exact ht.2.le)
  have a1 : |log t| ≤ ls t := by unfold ls; linarith [abs_nonneg (log (1 - t))]
  have := pow_le_pow_left₀ (abs_nonneg _) a1 k
  have := pow_nonneg (abs_nonneg (log t)) k
  nlinarith

lemma J_rec (n k : ℕ) : (n + 1 : ℝ) * Jl n (k+1) + (k + 1) * Jl n k = 0 := by
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
    ((ii_J n (k+1)).const_mul _).add ((ii_J n k).const_mul _)
  have := intervalIntegral.integral_eq_sub_of_hasDerivAt_of_le zero_le_one
    (G_cont n (k+1) (by omega)) hderiv hint
  rw [intervalIntegral.integral_add ((ii_J n (k+1)).const_mul _) ((ii_J n k).const_mul _),
    intervalIntegral.integral_const_mul, intervalIntegral.integral_const_mul] at this
  simpa [Jl] using this

lemma J0 (n : ℕ) : Jl n 0 = 1 / (n + 1) := by
  simp [Jl, integral_pow]

lemma J1 (n : ℕ) : Jl n 1 = -1 / (n + 1) ^ 2 := by
  have := J_rec n 0; rw [J0] at this; push_cast at this
  have hn : (n:ℝ) + 1 ≠ 0 := by positivity
  field_simp; field_simp at this; linarith

lemma J2 (n : ℕ) : Jl n 2 = 2 / (n + 1) ^ 3 := by
  have := J_rec n 1; rw [J1] at this; push_cast at this
  have hn : (n:ℝ) + 1 ≠ 0 := by positivity
  field_simp; field_simp at this; linarith

lemma J3 (n : ℕ) : Jl n 3 = -6 / (n + 1) ^ 4 := by
  have := J_rec n 2; rw [J2] at this; push_cast at this
  have hn : (n:ℝ) + 1 ≠ 0 := by positivity
  field_simp; field_simp at this; linarith

lemma hasSum_logmul (k : ℕ) :
    HasSum (fun n : ℕ => -(Jl (n+1) k) / (n+1)) (∫ t in (0:ℝ)..1, log t ^ k * log (1 - t)) := by
  have key : ∀ t ∈ Ioo (0:ℝ) 1, HasSum (fun n : ℕ => -(log t ^ k * (t ^ (n+1) / (n+1))))
      (log t ^ k * log (1 - t)) := by
    intro t ht
    have h := (hasSum_pow_div_log_of_abs_lt_one
      (show |t| < 1 by rw [abs_of_pos ht.1]; exact ht.2)).mul_left (log t ^ k)
    convert h.neg using 1; ring
  have hint : ∀ n : ℕ, ∫ t in (0:ℝ)..1, -(log t ^ k * (t ^ (n+1) / (n+1))) = -(Jl (n+1) k) / (n+1) := by
    intro n
    rw [intervalIntegral.integral_neg]
    have : (fun t : ℝ => log t ^ k * (t ^ (n+1) / (n+1))) =
        fun t => (1 / (n+1 : ℝ)) * (t ^ (n+1) * log t ^ k) := by funext t; ring
    rw [this, intervalIntegral.integral_const_mul, Jl]; ring
  simp_rw [← hint]
  have hg : IntervalIntegrable (fun t : ℝ => |log t| ^ k * (-log (1 - t))) volume 0 1 := by
    refine ii_of_ls (by fun_prop : Measurable (fun t : ℝ => |log t| ^ k * (-log (1 - t)))).aestronglyMeasurable
      1 (k+1) fun t ht => ?_
    rw [abs_mul, abs_pow, abs_abs, abs_neg, one_mul, pow_succ]
    have a1 : |log t| ≤ ls t := by unfold ls; linarith [abs_nonneg (log (1 - t))]
    have a2 : |log (1 - t)| ≤ ls t := by unfold ls; linarith [abs_nonneg (log t)]
    have := ls_nonneg t
    gcongr
  refine intervalIntegral.hasSum_integral_of_dominated_convergence
    (fun n t => |log t| ^ k * (t ^ (n+1) / (n+1))) (fun n => ?_) (fun n => ?_) ?_ ?_ ?_
  · exact (by fun_prop : Measurable (fun t : ℝ => -(log t ^ k * (t ^ (n+1) / (n+1))))).aestronglyMeasurable
  · filter_upwards [ae_Ioo01'] with t ht' ht0
    have ht := ht' ht0
    have hp : 0 < t ^ (n+1) / ((n:ℝ)+1) := by have := ht.1; positivity
    rw [Real.norm_eq_abs, abs_neg, abs_mul, abs_pow, abs_of_pos hp]
  · filter_upwards [ae_Ioo01'] with t ht' ht0
    have ht := ht' ht0
    exact ((hasSum_pow_div_log_of_abs_lt_one
      (show |t| < 1 by rw [abs_of_pos ht.1]; exact ht.2)).mul_left _).summable
  · refine hg.congr_uIoo ?_
    rw [uIoo01]
    intro t ht
    dsimp only
    rw [tsum_mul_left, (hasSum_pow_div_log_of_abs_lt_one
      (show |t| < 1 by rw [abs_of_pos ht.1]; exact ht.2)).tsum_eq]
  · filter_upwards [ae_Ioo01'] with t ht' ht0
    have ht := ht' ht0
    exact key t ht


lemma hs_tele : HasSum (fun n : ℕ => 1 / ((n:ℝ) + 1) - 1 / ((n:ℝ) + 2)) 1 := by
  rw [hasSum_iff_tendsto_nat_of_nonneg]
  · have e : ∀ N : ℕ, ∑ i ∈ Finset.range N, (1 / ((i:ℝ) + 1) - 1 / ((i:ℝ) + 2)) =
        1 - 1 / ((N:ℝ) + 1) := by
      intro N
      have := Finset.sum_range_sub' (fun i : ℕ => 1 / ((i:ℝ) + 1)) N
      simp only [Nat.cast_add, Nat.cast_one, Nat.cast_zero, zero_add, div_one] at this
      rw [← this]
      refine Finset.sum_congr rfl fun i _ => ?_
      ring
    simp_rw [e]
    simpa using (tendsto_one_div_add_atTop_nhds_zero_nat).const_sub (1:ℝ)
  · intro n
    rw [sub_nonneg]
    apply one_div_le_one_div_of_le (by positivity); linarith

lemma hs_sq : HasSum (fun n : ℕ => 1 / ((n:ℝ) + 2) ^ 2) (π ^ 2 / 6 - 1) := by
  have := (hasSum_nat_add_iff' 2).mpr hasSum_zeta_two
  simp [Finset.sum_range_succ] at this
  convert this using 2 with n
  push_cast; ring_nf

lemma hz3 : HasSum (fun n : ℕ => 1 / (n:ℝ) ^ 3) z3 := by
  have hs : Summable (fun n : ℕ => 1 / (n:ℝ) ^ 3) := Real.summable_one_div_nat_pow.mpr (by norm_num)
  have h := zeta_nat_eq_tsum_of_gt_one (k := 3) (by norm_num)
  have h2 : (∑' n : ℕ, 1 / (n:ℂ) ^ 3) = ((∑' n : ℕ, 1 / (n:ℝ) ^ 3 : ℝ) : ℂ) := by
    rw [Complex.ofReal_tsum]; push_cast; rfl
  have : z3 = ∑' n : ℕ, 1 / (n:ℝ) ^ 3 := by
    unfold z3
    rw [show (3:ℂ) = ((3:ℕ):ℂ) by norm_num, h, h2, Complex.ofReal_re]
  rw [this]; exact hs.hasSum

lemma hs_cube : HasSum (fun n : ℕ => 1 / ((n:ℝ) + 2) ^ 3) (z3 - 1) := by
  have := (hasSum_nat_add_iff' 2).mpr hz3
  simp [Finset.sum_range_succ] at this
  convert this using 2 with n
  push_cast; ring_nf

/-- `∫ log t log(1-t) = 2 - π²/6` -/
lemma int_ab : ∫ t in (0:ℝ)..1, log t * log (1 - t) = 2 - π ^ 2 / 6 := by
  have h := hasSum_logmul 1
  simp only [pow_one] at h
  refine h.unique ?_
  convert hs_tele.sub hs_sq using 1
  · funext n; rw [J1]; push_cast
    field_simp; ring
  · ring

/-- `∫ log² t log(1-t) = 2ζ(3) - 6 + π²/3` -/
lemma int_aab : ∫ t in (0:ℝ)..1, log t ^ 2 * log (1 - t) = 2 * z3 - 6 + π ^ 2 / 3 := by
  have h := hasSum_logmul 2
  refine h.unique ?_
  convert ((hs_tele.sub hs_sq).sub hs_cube).mul_left (-2) using 1
  · funext n; rw [J2]; push_cast
    field_simp; ring
  · ring

lemma int_a (k : ℕ) : ∫ t in (0:ℝ)..1, log t ^ k = Jl 0 k := by simp [Jl]

lemma int_symm (f : ℝ → ℝ) : ∫ t in (0:ℝ)..1, f (1 - t) = ∫ t in (0:ℝ)..1, f t := by
  rw [intervalIntegral.integral_comp_sub_left]; norm_num

lemma int_b (k : ℕ) : ∫ t in (0:ℝ)..1, log (1 - t) ^ k = Jl 0 k := by
  rw [int_symm (fun t => log t ^ k), int_a]

lemma int_abb : ∫ t in (0:ℝ)..1, log t * log (1 - t) ^ 2 = 2 * z3 - 6 + π ^ 2 / 3 := by
  rw [← int_aab, ← int_symm (fun t => log t ^ 2 * log (1 - t))]
  congr 1; funext t; simp only [sub_sub_cancel]; ring


noncomputable def ibf (p q i m : ℕ) (y : ℝ) : ℝ :=
  y ^ (p - 1) * (1 - y) ^ (q - 1) * Real.log y ^ i * Real.log (1 - y) ^ m

noncomputable def Ib (p q i m : ℕ) : ℝ := ∫ y in (0 : ℝ)..1, ibf p q i m y

lemma ibf_int (p q i m : ℕ) : IntervalIntegrable (ibf p q i m) volume 0 1 := by
  refine ii_of_ls (by unfold ibf; fun_prop : Measurable (ibf p q i m)).aestronglyMeasurable
    1 (i + m) fun t ht => ?_
  have h0 := ht.1; have h1 := ht.2
  unfold ibf
  rw [abs_mul, abs_mul, abs_mul, abs_pow, abs_pow, abs_pow, abs_pow, one_mul, pow_add]
  have a0 : |t| ^ (p - 1) ≤ 1 := pow_le_one₀ (abs_nonneg t) (by rw [abs_of_pos h0]; exact h1.le)
  have b0 : |1 - t| ^ (q - 1) ≤ 1 :=
    pow_le_one₀ (abs_nonneg _) (by rw [abs_of_pos (by linarith)]; linarith)
  have a1 : |Real.log t| ≤ ls t := by unfold ls; linarith [abs_nonneg (Real.log (1 - t))]
  have a2 : |Real.log (1 - t)| ≤ ls t := by unfold ls; linarith [abs_nonneg (Real.log t)]
  have c1 := pow_le_pow_left₀ (abs_nonneg _) a1 i
  have c2 := pow_le_pow_left₀ (abs_nonneg _) a2 m
  have := ls_nonneg t
  calc |t| ^ (p - 1) * |1 - t| ^ (q - 1) * |Real.log t| ^ i * |Real.log (1 - t)| ^ m
      ≤ 1 * 1 * ls t ^ i * ls t ^ m := by gcongr
    _ = _ := by ring

lemma Ib_symm (p q i m : ℕ) : Ib p q i m = Ib q p m i := by
  unfold Ib
  rw [← int_symm (ibf q p m i)]
  congr 1; funext y; unfold ibf; simp only [sub_sub_cancel]; ring

lemma Ib_pascal {p q : ℕ} (hp : 1 ≤ p) (hq : 1 ≤ q) (i m : ℕ) :
    Ib p q i m = Ib (p + 1) q i m + Ib p (q + 1) i m := by
  unfold Ib
  rw [← intervalIntegral.integral_add (ibf_int _ _ _ _) (ibf_int _ _ _ _)]
  congr 1; funext y
  obtain ⟨a, rfl⟩ : ∃ a, p = a + 1 := ⟨p - 1, by omega⟩
  obtain ⟨b, rfl⟩ : ∃ b, q = b + 1 := ⟨q - 1, by omega⟩
  unfold ibf; simp only [Nat.add_sub_cancel]; ring

lemma pow_log_tendsto (p i : ℕ) (hp : 1 ≤ p) :
    Tendsto (fun y : ℝ => y ^ p * Real.log y ^ i) (𝓝[>] 0) (𝓝 0) := by
  obtain ⟨a, rfl⟩ : ∃ a, p = a + 1 := ⟨p - 1, by omega⟩
  rcases Nat.eq_zero_or_pos i with hi | hi
  · subst hi
    have : ContinuousAt (fun y : ℝ => y ^ (a + 1) * Real.log y ^ 0) 0 := by
      simp only [pow_zero, mul_one]; exact (continuous_pow _).continuousAt
    have h := this.tendsto
    simp only [pow_zero, mul_one, ne_eq, Nat.add_eq_zero_iff, one_ne_zero, and_false,
      not_false_eq_true, zero_pow] at h ⊢
    exact h.mono_left nhdsWithin_le_nhds
  · have h := (((continuous_pow a).tendsto 0).mono_left nhdsWithin_le_nhds).mul
      (t_log_pow_tendsto i hi.ne')
    rw [mul_zero] at h
    refine h.congr fun y => ?_
    ring

lemma G_tendsto0 (p q i m : ℕ) (hp : 1 ≤ p) :
    Tendsto (fun y : ℝ => y ^ p * (1 - y) ^ q * Real.log y ^ i * Real.log (1 - y) ^ m)
      (𝓝[>] 0) (𝓝 0) := by
  have h1 := pow_log_tendsto p i hp
  have h2 : ContinuousAt (fun y : ℝ => (1 - y) ^ q * Real.log (1 - y) ^ m) 0 := by
    have : ContinuousAt (fun y : ℝ => Real.log (1 - y)) 0 :=
      (continuousAt_log (by norm_num)).comp (continuousAt_const.sub continuousAt_id)
    exact ((continuousAt_const.sub continuousAt_id).pow q).mul (this.pow m)
  have h := h1.mul (h2.tendsto.mono_left nhdsWithin_le_nhds)
  rw [zero_mul] at h
  refine h.congr fun y => ?_
  ring

lemma G_tendsto1 (p q i m : ℕ) (hq : 1 ≤ q) :
    Tendsto (fun y : ℝ => y ^ p * (1 - y) ^ q * Real.log y ^ i * Real.log (1 - y) ^ m)
      (𝓝[<] 1) (𝓝 0) := by
  have hs : Tendsto (fun y : ℝ => 1 - y) (𝓝[<] 1) (𝓝[>] 0) := by
    refine tendsto_nhdsWithin_of_tendsto_nhds_of_eventually_within _ ?_ ?_
    · have : Tendsto (fun y : ℝ => 1 - y) (𝓝 1) (𝓝 (1 - 1)) :=
        (continuous_const.sub continuous_id).tendsto 1
      rw [sub_self] at this
      exact this.mono_left nhdsWithin_le_nhds
    · filter_upwards [self_mem_nhdsWithin] with y hy
      simp only [mem_Ioi, sub_pos]; exact hy
  have h0 := (G_tendsto0 q p m i hq).comp hs
  refine h0.congr fun y => ?_
  simp only [Function.comp, sub_sub_cancel]; ring

lemma Ib_ibp {p q : ℕ} (hp : 1 ≤ p) (hq : 1 ≤ q) (i m : ℕ) :
    (p : ℝ) * Ib p (q + 1) i m - q * Ib (p + 1) q i m + i * Ib p (q + 1) (i - 1) m
      - m * Ib (p + 1) q i (m - 1) = 0 := by
  have hderiv : ∀ y ∈ Ioo (0 : ℝ) 1, HasDerivAt
      (fun y : ℝ => y ^ p * (1 - y) ^ q * Real.log y ^ i * Real.log (1 - y) ^ m)
      ((p : ℝ) * ibf p (q + 1) i m y - q * ibf (p + 1) q i m y + i * ibf p (q + 1) (i - 1) m y
        - m * ibf (p + 1) q i (m - 1) y) y := by
    intro y hy
    have hy0 : y ≠ 0 := hy.1.ne'
    have hy1 : 1 - y ≠ 0 := by linarith [hy.2]
    have d1 : HasDerivAt (fun y : ℝ => y ^ p) ((p : ℝ) * y ^ (p - 1)) y := hasDerivAt_pow p y
    have d2 : HasDerivAt (fun y : ℝ => (1 - y) ^ q) ((q : ℝ) * (1 - y) ^ (q - 1) * (-1)) y := by
      have := (hasDerivAt_pow q (1 - y)).comp y ((hasDerivAt_id' y).const_sub 1)
      simpa [Function.comp_def] using this
    have d3 : HasDerivAt (fun y : ℝ => Real.log y ^ i) ((i : ℝ) * Real.log y ^ (i - 1) * y⁻¹) y := by
      have := (hasDerivAt_pow i (Real.log y)).comp y (Real.hasDerivAt_log hy0)
      simpa [Function.comp_def] using this
    have d4 : HasDerivAt (fun y : ℝ => Real.log (1 - y) ^ m)
        ((m : ℝ) * Real.log (1 - y) ^ (m - 1) * ((1 - y)⁻¹ * (-1))) y := by
      have h4 := ((Real.hasDerivAt_log hy1).comp y ((hasDerivAt_id' y).const_sub 1))
      have := (hasDerivAt_pow m (Real.log (1 - y))).comp y h4
      simpa [Function.comp_def] using this
    have h := ((d1.mul d2).mul d3).mul d4
    refine h.congr_deriv ?_
    obtain ⟨a, rfl⟩ : ∃ a, p = a + 1 := ⟨p - 1, by omega⟩
    obtain ⟨b, rfl⟩ : ∃ b, q = b + 1 := ⟨q - 1, by omega⟩
    simp only [Pi.mul_apply]
    unfold ibf
    simp only [Nat.add_sub_cancel]
    rw [pow_succ y a, pow_succ (1 - y) b]
    field_simp
    push_cast
    ring
  have hint : IntervalIntegrable (fun y => (p : ℝ) * ibf p (q + 1) i m y - q * ibf (p + 1) q i m y
      + i * ibf p (q + 1) (i - 1) m y - m * ibf (p + 1) q i (m - 1) y) volume 0 1 :=
    ((((ibf_int _ _ _ _).const_mul _).sub ((ibf_int _ _ _ _).const_mul _)).add
      ((ibf_int _ _ _ _).const_mul _)).sub ((ibf_int _ _ _ _).const_mul _)
  have key := intervalIntegral.integral_eq_sub_of_hasDerivAt_of_tendsto zero_lt_one hderiv hint
    (G_tendsto0 p q i m hp) (G_tendsto1 p q i m hq)
  rw [sub_zero, intervalIntegral.integral_sub _ ((ibf_int _ _ _ _).const_mul _),
    intervalIntegral.integral_add _ ((ibf_int _ _ _ _).const_mul _),
    intervalIntegral.integral_sub ((ibf_int _ _ _ _).const_mul _) ((ibf_int _ _ _ _).const_mul _),
    intervalIntegral.integral_const_mul, intervalIntegral.integral_const_mul,
    intervalIntegral.integral_const_mul, intervalIntegral.integral_const_mul] at key
  · unfold Ib; linarith
  · exact ((ibf_int _ _ _ _).const_mul _).sub ((ibf_int _ _ _ _).const_mul _)
  · exact (((ibf_int _ _ _ _).const_mul _).sub ((ibf_int _ _ _ _).const_mul _)).add
      ((ibf_int _ _ _ _).const_mul _)

lemma Ib_rec {p q : ℕ} (hp : 1 ≤ p) (hq : 1 ≤ q) (i m : ℕ) :
    ((p : ℝ) + q) * Ib p (q + 1) i m = q * Ib p q i m - i * Ib p (q + 1) (i - 1) m
      + m * (Ib p q i (m - 1) - Ib p (q + 1) i (m - 1)) := by
  have h1 := Ib_ibp hp hq i m
  have h2 := Ib_pascal hp hq i m
  have h3 := Ib_pascal hp hq i (m - 1)
  linear_combination h1 - (q : ℝ) * h2 - (m : ℝ) * h3

noncomputable def DF (i m : ℕ) (a1 b1 c1 a2 b2 c2 a3 b3 c3 x2 x3 : ℝ) : ℝ :=
  let A := a1 - c1
  let B := b1 - c1
  let P := a2 - c2
  let Q := b2 - c2
  let R := c2 - x2
  let S := a3 - c3
  let T := b3 - c3
  let W := x3 - c3
  match i, m with
  | 0, 0 => 1
  | 1, 0 => A
  | 0, 1 => B
  | 2, 0 => A ^ 2 - P
  | 0, 2 => B ^ 2 - Q
  | 1, 1 => A * B + R
  | 3, 0 => 2 * S - 3 * A * P + A ^ 3
  | 0, 3 => 2 * T - 3 * B * Q + B ^ 3
  | 2, 1 => 2 * W + 2 * A * R - B * P + A ^ 2 * B
  | 1, 2 => 2 * W + 2 * B * R - A * Q + A * B ^ 2
  | _, _ => 0

lemma Dsharp_DF (i m p q : ℕ) (x2 x3 : ℝ) : Dsharp i m p q x2 x3 =
    DF i m (Hr 1 (p - 1)) (Hr 1 (q - 1)) (Hr 1 (p + q - 1)) (Hr 2 (p - 1)) (Hr 2 (q - 1))
      (Hr 2 (p + q - 1)) (Hr 3 (p - 1)) (Hr 3 (q - 1)) (Hr 3 (p + q - 1)) x2 x3 := rfl

noncomputable def Fm (p q i m : ℕ) : ℝ := Bet p q * Dsharp i m p q z2 z3

lemma Bet_pos (p q : ℕ) : 0 < Bet p q := by unfold Bet; positivity

lemma Bet_succ (a b : ℕ) : Bet (a + 1) (b + 1 + 1) = Bet (a + 1) (b + 1) * ((b : ℝ) + 1) / ((a : ℝ) + b + 2) := by
  unfold Bet
  have e1 : a + 1 + (b + 1 + 1) - 1 = (a + b + 1) + 1 := by omega
  have e2 : a + 1 + (b + 1) - 1 = a + b + 1 := by omega
  simp only [e1, e2, Nat.add_sub_cancel, Nat.factorial_succ]
  push_cast
  have : ((a + b + 1).factorial : ℝ) ≠ 0 := by positivity
  field_simp
  ring

lemma Fm_rec {p q : ℕ} (hp : 1 ≤ p) (hq : 1 ≤ q) {i m : ℕ} (h : i + m ≤ 3) :
    ((p : ℝ) + q) * Fm p (q + 1) i m = q * Fm p q i m - i * Fm p (q + 1) (i - 1) m
      + m * (Fm p q i (m - 1) - Fm p (q + 1) i (m - 1)) := by
  obtain ⟨a, rfl⟩ : ∃ a, p = a + 1 := ⟨p - 1, by omega⟩
  obtain ⟨b, rfl⟩ : ∃ b, q = b + 1 := ⟨q - 1, by omega⟩
  unfold Fm
  simp only [Dsharp_DF]
  have e1 : a + 1 + (b + 1 + 1) - 1 = (a + b + 1) + 1 := by omega
  have e2 : a + 1 + (b + 1) - 1 = a + b + 1 := by omega
  have e3 : b + 1 + 1 - 1 = b + 1 := by omega
  simp only [e1, e2, e3, Nat.add_sub_cancel, Bet_succ]
  have h1 : ∀ r, Hr r (a + b + 1 + 1) = Hr r (a + b + 1) + 1 / ((a : ℝ) + b + 2) ^ r := by
    intro r; rw [Hr_succ]; push_cast; ring_nf
  have h2 : ∀ r, Hr r (b + 1) = Hr r b + 1 / ((b : ℝ) + 1) ^ r := fun r => Hr_succ r b
  simp only [h1, h2]
  have hB := Bet_pos (a + 1) (b + 1)
  generalize Bet (a + 1) (b + 1) = B at hB ⊢
  generalize Hr 1 a = a1
  generalize Hr 1 b = b1
  generalize Hr 1 (a + b + 1) = c1
  generalize Hr 2 a = a2
  generalize Hr 2 b = b2
  generalize Hr 2 (a + b + 1) = c2
  generalize Hr 3 a = a3
  generalize Hr 3 b = b3
  generalize Hr 3 (a + b + 1) = c3
  have hb : (b : ℝ) + 1 ≠ 0 := by positivity
  have hab : (a : ℝ) + b + 2 ≠ 0 := by positivity
  push_cast
  have hi : i ≤ 3 := by omega
  have hm : m ≤ 3 := by omega
  interval_cases i <;> interval_cases m <;> first | omega | (simp only [DF]; field_simp; ring)

lemma Fm_symm (p q i m : ℕ) (h : i + m ≤ 3) : Fm p q i m = Fm q p m i := by
  unfold Fm
  simp only [Dsharp_DF]
  have e : q + p - 1 = p + q - 1 := by omega
  have hB : Bet p q = Bet q p := by unfold Bet; rw [e]; ring
  rw [e, hB]
  have hi : i ≤ 3 := by omega
  have hm : m ≤ 3 := by omega
  interval_cases i <;> interval_cases m <;> first | omega | (simp only [DF] <;> ring)

lemma Hr_zero (r : ℕ) : Hr r 0 = 0 := by simp [Hr]
lemma Hr_one (r : ℕ) : Hr r 1 = 1 := by simp [Hr]

lemma Ib_base (i m : ℕ) (h : i + m ≤ 3) : Ib 1 1 i m = Fm 1 1 i m := by
  have hB : Bet 1 1 = 1 := by simp [Bet]
  have hI : Ib 1 1 i m = ∫ y in (0 : ℝ)..1, Real.log y ^ i * Real.log (1 - y) ^ m := by
    unfold Ib ibf; congr 1; funext y; simp
  rw [hI]
  unfold Fm
  rw [hB, one_mul, Dsharp_DF]
  simp only [Nat.sub_self, show 1 + 1 - 1 = 1 by rfl, Hr_zero, Hr_one]
  have i00 : ∫ y in (0 : ℝ)..1, (1 : ℝ) = 1 := by simp
  have i10 : ∫ y in (0 : ℝ)..1, Real.log y = -1 := by
    have := int_a 1; simp only [pow_one] at this; rw [this, J1]; norm_num
  have i20 : ∫ y in (0 : ℝ)..1, Real.log y ^ 2 = 2 := by rw [int_a, J2]; norm_num
  have i30 : ∫ y in (0 : ℝ)..1, Real.log y ^ 3 = -6 := by rw [int_a, J3]; norm_num
  have i01 : ∫ y in (0 : ℝ)..1, Real.log (1 - y) = -1 := by
    have := int_b 1; simp only [pow_one] at this; rw [this, J1]; norm_num
  have i02 : ∫ y in (0 : ℝ)..1, Real.log (1 - y) ^ 2 = 2 := by rw [int_b, J2]; norm_num
  have i03 : ∫ y in (0 : ℝ)..1, Real.log (1 - y) ^ 3 = -6 := by rw [int_b, J3]; norm_num
  have i11 := int_ab
  have i21 := int_aab
  have i12 := int_abb
  have hi : i ≤ 3 := by omega
  have hm : m ≤ 3 := by omega
  interval_cases i <;> interval_cases m <;> first | omega |
    (simp only [DF, pow_zero, mul_one, one_mul, pow_one, i00, i10, i20, i30, i01, i02, i03, i11,
      i21, i12, z2] <;> ring)

lemma mul_eq_of {i : ℕ} {X Y : ℝ} (h : i = 0 ∨ X = Y) : (i : ℝ) * X = i * Y := by
  rcases h with h | h
  · subst h; simp
  · rw [h]

lemma Ib_eq_col {p : ℕ} (hp : 1 ≤ p) (hbase : ∀ i m, i + m ≤ 3 → Ib p 1 i m = Fm p 1 i m) :
    ∀ q, 1 ≤ q → ∀ i m, i + m ≤ 3 → Ib p q i m = Fm p q i m := by
  intro q hq
  induction q, hq using Nat.le_induction with
  | base => exact hbase
  | succ q hq ih =>
    have key : ∀ s, ∀ i m, i + m ≤ s → s ≤ 3 → Ib p (q + 1) i m = Fm p (q + 1) i m := by
      intro s
      induction s with
      | zero =>
        intro i m h _
        have hi : i = 0 := by omega
        have hm : m = 0 := by omega
        subst hi hm
        have r1 := Ib_rec hp hq 0 0
        have r2 := Fm_rec hp hq (i := 0) (m := 0) (by norm_num)
        rw [ih 0 0 (by norm_num)] at r1
        have hpq : (p : ℝ) + q ≠ 0 := by positivity
        apply mul_left_cancel₀ hpq
        rw [r1, r2]; simp
      | succ s ihs =>
        intro i m h hs
        have r1 := Ib_rec hp hq i m
        have r2 := Fm_rec hp hq (i := i) (m := m) (by omega)
        have e1 : (i : ℝ) * Ib p (q + 1) (i - 1) m = i * Fm p (q + 1) (i - 1) m := by
          apply mul_eq_of
          rcases Nat.eq_zero_or_pos i with h0 | h0
          · left; exact h0
          · right; exact ihs (i - 1) m (by omega) (by omega)
        have e2 : (m : ℝ) * (Ib p q i (m - 1) - Ib p (q + 1) i (m - 1))
            = m * (Fm p q i (m - 1) - Fm p (q + 1) i (m - 1)) := by
          apply mul_eq_of
          rcases Nat.eq_zero_or_pos m with h0 | h0
          · left; exact h0
          · right; rw [ih i (m - 1) (by omega), ihs i (m - 1) (by omega) (by omega)]
        rw [ih i m (by omega), e1, e2] at r1
        have hpq : (p : ℝ) + q ≠ 0 := by positivity
        apply mul_left_cancel₀ hpq
        rw [r1, r2]
    intro i m h
    exact key 3 i m h le_rfl

lemma Ib_eq_Fm {p q : ℕ} (hp : 1 ≤ p) (hq : 1 ≤ q) {i m : ℕ} (h : i + m ≤ 3) :
    Ib p q i m = Fm p q i m := by
  have col1 := Ib_eq_col (le_refl 1) Ib_base
  have hbase : ∀ i m, i + m ≤ 3 → Ib p 1 i m = Fm p 1 i m := by
    intro i m h
    rw [Ib_symm, col1 p hp m i (by omega), Fm_symm 1 p m i (by omega)]
  exact Ib_eq_col hp hbase q hq i m h

end V05

theorem V05.proof_7 : V05.stmt_7 := by
  intro p q hp hq i m h
  exact ⟨V05.ibf_int p q i m, V05.Ib_eq_Fm hp hq h⟩

namespace V05
open MeasureTheory Set Filter Real Interval Topology

lemma ecoef_deg {k : ℕ} (hk : k ≤ 3) (n i : ℕ) (hi : i ≤ k) :
    (ecoef k n i (k - i)).natDegree < k + 1 := by
  apply Nat.lt_succ_of_le
  interval_cases k <;> interval_cases i <;> simp only [ecoef, al1, be1, ga, ga3] <;> compute_degree!

lemma Bet_comm (p q : ℕ) : Bet p q = Bet q p := by
  unfold Bet; rw [show q + p - 1 = p + q - 1 by omega]; ring

lemma Bet_succ_left (a b : ℕ) :
    Bet (a + 1 + 1) (b + 1) = Bet (a + 1) (b + 1) * ((a : ℝ) + 1) / ((a : ℝ) + b + 2) := by
  rw [Bet_comm, Bet_succ b a, Bet_comm]; ring

lemma Bet_rho {N : ℕ} (hN : 1 ≤ N) (j : ℕ) : Bet (N + j) N = Bet N N * rho j N := by
  induction j with
  | zero => simp [rho]
  | succ j ih =>
    obtain ⟨M, rfl⟩ : ∃ M, N = M + 1 := ⟨N - 1, by omega⟩
    have e : M + 1 + (j + 1) = (M + j) + 1 + 1 := by omega
    rw [e, Bet_succ_left (M + j) M, show M + j + 1 = M + 1 + j by omega, ih]
    simp only [rho, Finset.prod_range_succ]
    have h1 : (2 * ((M + 1 : ℕ) : ℝ) + j) ≠ 0 := by positivity
    have h2 : ((M + j : ℕ) : ℝ) + M + 2 ≠ 0 := by positivity
    push_cast at h1 h2 ⊢
    field_simp
    ring

lemma an_Bet {N : ℕ} (hN : 1 ≤ N) : an N * Bet N N / 2 = 1 / ((N : ℝ) * (2 * (N : ℝ) - 1)) := by
  obtain ⟨M, rfl⟩ : ∃ M, N = M + 1 := ⟨N - 1, by omega⟩
  unfold an Bet
  have e1 : M + 1 + (M + 1) - 1 = 2 * M + 1 := by omega
  rw [e1, Nat.add_sub_cancel]
  have hc := Nat.choose_mul_factorial_mul_factorial (show M + 1 ≤ 2 * (M + 1) by omega)
  rw [show 2 * (M + 1) - (M + 1) = M + 1 by omega] at hc
  have hc' : ((Nat.choose (2 * (M + 1)) (M + 1) : ℕ) : ℝ) * ((M + 1).factorial : ℝ) * ((M + 1).factorial : ℝ)
      = ((2 * (M + 1)).factorial : ℝ) := by exact_mod_cast hc
  rw [show 2 * (M + 1) = (2 * M + 1) + 1 by omega, Nat.factorial_succ (2 * M + 1)] at hc'
  rw [show 2 * (M + 1) = (2 * M + 1) + 1 by omega]
  simp only [Nat.factorial_succ M] at hc'
  push_cast at hc' ⊢
  have h1 : ((2 * M + 1).factorial : ℝ) ≠ 0 := by positivity
  have h2 : (M.factorial : ℝ) ≠ 0 := by positivity
  have h3 : (2 * ((M : ℝ) + 1) - 1) ≠ 0 := by
    have : (0 : ℝ) ≤ M := Nat.cast_nonneg M
    linarith
  have h4 : ((M : ℝ) + 1) ≠ 0 := by positivity
  field_simp
  nlinarith [hc']

lemma Jk_formula {N : ℕ} (hN : 1 ≤ N) {k : ℕ} (hk : k ≤ 3) :
    iteratedDeriv k (J N) 0 / (k.factorial : ℝ) = 1 / 2 * ∑ j ∈ Finset.range (k + 1),
      ∑ i ∈ Finset.range (k + 1), (ecoef k N i (k - i)).coeff j * Bet (N + j) N
        * Dsharp i (k - i) (N + j) N z2 z3 := by
  rw [(V05.proof_5 N hN k).2]
  congr 1
  have hpt : ∀ y ∈ uIoo (0 : ℝ) 1, (y * (1 - y)) ^ (N - 1) *
      (iteratedDeriv k (fun u => E N u y) 0 / (k.factorial : ℝ)) =
      ∑ j ∈ Finset.range (k + 1), ∑ i ∈ Finset.range (k + 1),
        (ecoef k N i (k - i)).coeff j * ibf (N + j) N i (k - i) y := by
    intro y hy
    rw [uIoo_of_le zero_le_one] at hy
    rw [V05.proof_6 N hN y hy.1 hy.2 k hk, Finset.sum_comm, Finset.mul_sum]
    refine Finset.sum_congr rfl fun i hi => ?_
    have hi' : i ≤ k := by simp at hi; omega
    rw [Polynomial.eval_eq_sum_range' (ecoef_deg hk N i hi') y, Finset.sum_mul, Finset.sum_mul,
      Finset.mul_sum]
    refine Finset.sum_congr rfl fun j _ => ?_
    unfold ibf
    rw [show N + j - 1 = (N - 1) + j by omega, pow_add, mul_pow]
    ring
  rw [intervalIntegral.integral_congr_uIoo hpt]
  have hint : ∀ j i : ℕ, IntervalIntegrable
      (fun y => (ecoef k N i (k - i)).coeff j * ibf (N + j) N i (k - i) y) volume 0 1 :=
    fun j i => (ibf_int _ _ _ _).const_mul _
  rw [intervalIntegral.integral_finsetSum]
  · refine Finset.sum_congr rfl fun j _ => ?_
    rw [intervalIntegral.integral_finsetSum (fun i _ => hint j i)]
    refine Finset.sum_congr rfl fun i hi => ?_
    have hi' : i ≤ k := by simp at hi; omega
    rw [intervalIntegral.integral_const_mul]
    have h7 := (V05.proof_7 (N + j) N (by omega) hN i (k - i) (by omega)).2
    change (ecoef k N i (k - i)).coeff j * Ib (N + j) N i (k - i) = _
    unfold Ib
    rw [show (∫ y in (0 : ℝ)..1, ibf (N + j) N i (k - i) y) = Bet (N + j) N
      * Dsharp i (k - i) (N + j) N z2 z3 from h7]
    ring
  · intro j _
    have := IntervalIntegrable.sum (μ := volume) (a := 0) (b := 1) (Finset.range (k + 1))
      (f := fun i y => (ecoef k N i (k - i)).coeff j * ibf (N + j) N i (k - i) y)
      (fun i _ => hint j i)
    refine this.congr ?_
    intro y _
    simp only [Finset.sum_apply]

end V05

theorem V05.proof_8 : V05.stmt_8 := by
  intro α hα0 hα1 k hk
  obtain ⟨δ, -, -, -, hs⟩ := V05.proof_4 α hα0 hα1
  have h := (hs k).neg
  rw [neg_neg] at h
  convert h using 1
  funext n
  have hN : 1 ≤ n + 1 := by omega
  rw [V05.Jk_formula hN hk]
  unfold V05.P
  simp_rw [V05.Bet_rho hN]
  have e : ∀ j ∈ Finset.range (k + 1), ∑ i ∈ Finset.range (k + 1),
      (V05.ecoef k (n + 1) i (k - i)).coeff j * (V05.Bet (n + 1) (n + 1) * V05.rho j (n + 1))
        * V05.Dsharp i (k - i) (n + 1 + j) (n + 1) V05.z2 V05.z3
      = V05.Bet (n + 1) (n + 1) * ∑ i ∈ Finset.range (k + 1),
        (V05.ecoef k (n + 1) i (k - i)).coeff j * V05.rho j (n + 1)
          * V05.Dsharp i (k - i) (n + 1 + j) (n + 1) V05.z2 V05.z3 := by
    intro j _
    rw [Finset.mul_sum]
    exact Finset.sum_congr rfl fun i _ => by ring
  rw [Finset.sum_congr rfl e, ← Finset.mul_sum, ← V05.an_Bet hN]
  ring

namespace V05
open MeasureTheory Set Filter Real Interval Topology

lemma e321_eq (n : ℕ) : ecoef 3 n 2 1 =
    Polynomial.C ((n : ℝ) / 2 - (n : ℝ) ^ 2 + (n : ℝ) ^ 3 / 2)
    + Polynomial.C (-(n : ℝ) / 2 + 5 * (n : ℝ) ^ 2 / 2 - 3 * (n : ℝ) ^ 3) * Polynomial.X
    + Polynomial.C (-(n : ℝ) - (n : ℝ) ^ 2 + 6 * (n : ℝ) ^ 3) * Polynomial.X ^ 2
    + Polynomial.C ((n : ℝ) - 4 * (n : ℝ) ^ 3) * Polynomial.X ^ 3 := by
  apply Polynomial.funext; intro y
  simp [ecoef, al1, be1, ga, ga3]; ring

lemma e312_eq (n : ℕ) : ecoef 3 n 1 2 =
    Polynomial.C ((n : ℝ) ^ 2 / 2 - (n : ℝ) ^ 3 / 2)
    + Polynomial.C (-(n : ℝ) / 2 - (n : ℝ) ^ 2 / 2 + 3 * (n : ℝ) ^ 3) * Polynomial.X
    + Polynomial.C (2 * (n : ℝ) - (n : ℝ) ^ 2 - 6 * (n : ℝ) ^ 3) * Polynomial.X ^ 2
    + Polynomial.C (-(n : ℝ) + 4 * (n : ℝ) ^ 3) * Polynomial.X ^ 3 := by
  apply Polynomial.funext; intro y
  simp [ecoef, al1, be1, ga, ga3]; ring

lemma coeff_bd3 {a0 a1 a2 a3 B : ℝ} (hB : 0 ≤ B) (h0 : |a0| ≤ B) (h1 : |a1| ≤ B) (h2 : |a2| ≤ B)
    (h3 : |a3| ≤ B) (j : ℕ) :
    |(Polynomial.C a0 + Polynomial.C a1 * Polynomial.X + Polynomial.C a2 * Polynomial.X ^ 2
      + Polynomial.C a3 * Polynomial.X ^ 3).coeff j| ≤ B := by
  rw [coeff3]
  rcases j with _ | _ | _ | _ | j <;> simp [h0, h1, h2, h3, hB]

lemma ecoef3_coeff_bd {n : ℕ} (hn : 1 ≤ n) {i : ℕ} (hi : i = 1 ∨ i = 2) (j : ℕ) :
    |(ecoef 3 n i (3 - i)).coeff j| ≤ 8 * (n : ℝ) ^ 3 := by
  have hn' : (1 : ℝ) ≤ n := by exact_mod_cast hn
  have h2 : (n : ℝ) ≤ (n : ℝ) ^ 2 := by nlinarith
  have h3 : (n : ℝ) ^ 2 ≤ (n : ℝ) ^ 3 := by nlinarith
  have hB : (0 : ℝ) ≤ 8 * (n : ℝ) ^ 3 := by positivity
  rcases hi with rfl | rfl
  · rw [show 3 - 1 = 2 from rfl, e312_eq]
    apply coeff_bd3 hB <;> rw [abs_le] <;> constructor <;> nlinarith
  · rw [show 3 - 2 = 1 from rfl, e321_eq]
    apply coeff_bd3 hB <;> rw [abs_le] <;> constructor <;> nlinarith

lemma rho_bd (j n : ℕ) : 0 ≤ rho j n ∧ rho j n ≤ 1 := by
  unfold rho
  refine ⟨Finset.prod_nonneg fun t _ => by positivity, Finset.prod_le_one₀ (fun t _ => by positivity)
    fun t _ => ?_⟩
  rcases Nat.eq_zero_or_pos n with h | h
  · subst h; simp only [Nat.cast_zero, zero_add, mul_zero]; exact div_self_le_one _
  · have : (0 : ℝ) < 2 * n + t := by positivity
    rw [div_le_one this]
    have : (0 : ℝ) ≤ n := Nat.cast_nonneg n
    linarith

lemma Hr1_bd (N : ℕ) : 0 ≤ Hr 1 N ∧ Hr 1 N ≤ N := by
  unfold Hr
  refine ⟨Finset.sum_nonneg fun l _ => by positivity, ?_⟩
  calc ∑ l ∈ Finset.range N, 1 / ((l : ℝ) + 1) ^ 1 ≤ ∑ l ∈ Finset.range N, (1 : ℝ) := by
        refine Finset.sum_le_sum fun l _ => ?_
        rw [pow_one, div_le_one (by positivity)]
        have : (0 : ℝ) ≤ l := Nat.cast_nonneg l
        linarith
    _ = N := by simp

lemma D2_bd (i m p q : ℕ) : |D2 i m p q| ≤ 2 * ((p : ℝ) + q + 1) := by
  have hp : (0 : ℝ) ≤ p := Nat.cast_nonneg p
  have hq : (0 : ℝ) ≤ q := Nat.cast_nonneg q
  have hdiff : ∀ a b : ℕ, a ≤ p + q → b ≤ p + q → |Hr 1 a - Hr 1 b| ≤ (p : ℝ) + q := by
    intro a b ha hb
    obtain ⟨a0, a1⟩ := Hr1_bd a
    obtain ⟨b0, b1⟩ := Hr1_bd b
    have ha' : (a : ℝ) ≤ p + q := by exact_mod_cast ha
    have hb' : (b : ℝ) ≤ p + q := by exact_mod_cast hb
    rw [abs_le]; constructor <;> linarith
  rcases i with _ | _ | _ | _ | i <;> rcases m with _ | _ | _ | _ | m <;> simp only [D2] <;>
    first
    | (rw [abs_zero]; positivity)
    | (rw [abs_neg, abs_one]; linarith)
    | (rw [abs_mul, abs_neg, abs_two]
       have := hdiff (p - 1) (p + q - 1) (by omega) (by omega)
       linarith)
    | (rw [abs_mul, abs_neg, abs_two]
       have := hdiff (q - 1) (p + q - 1) (by omega) (by omega)
       linarith)

lemma L23_bd {n : ℕ} (hn : 1 ≤ n) : |L2 3 n| ≤ 1536 * (n : ℝ) ^ 4 := by
  have hn' : (1 : ℝ) ≤ n := by exact_mod_cast hn
  have hpre : |(1 / ((n : ℝ) * (2 * (n : ℝ) - 1)))| ≤ 1 := by
    rw [abs_of_pos (by apply div_pos one_pos; nlinarith), div_le_one (by nlinarith)]
    nlinarith
  have hterm : ∀ j ∈ Finset.range 4, ∀ i ∈ Finset.range 4,
      |(ecoef 3 n i (3 - i)).coeff j * rho j n * D2 i (3 - i) (n + j) n| ≤ 96 * (n : ℝ) ^ 4 := by
    intro j hj i hi
    simp only [Finset.mem_range] at hj hi
    obtain ⟨r0, r1⟩ := rho_bd j n
    have hD := D2_bd i (3 - i) (n + j) n
    have hj' : (j : ℝ) ≤ 3 := by exact_mod_cast (by omega : j ≤ 3)
    have hD' : |D2 i (3 - i) (n + j) n| ≤ 12 * n := by
      refine hD.trans ?_; push_cast; nlinarith
    rw [abs_mul, abs_mul, abs_of_nonneg r0]
    have hn4 : (0 : ℝ) ≤ 96 * (n : ℝ) ^ 4 := by positivity
    interval_cases i
    · simp [D2]
    · have hc := ecoef3_coeff_bd hn (i := 1) (Or.inl rfl) j
      calc |(ecoef 3 n 1 (3 - 1)).coeff j| * rho j n * |D2 1 (3 - 1) (n + j) n|
          ≤ (8 * (n : ℝ) ^ 3) * 1 * (12 * n) := by
            apply mul_le_mul (mul_le_mul hc r1 r0 (by positivity)) hD' (abs_nonneg _) (by positivity)
        _ = 96 * (n : ℝ) ^ 4 := by ring
    · have hc := ecoef3_coeff_bd hn (i := 2) (Or.inr rfl) j
      calc |(ecoef 3 n 2 (3 - 2)).coeff j| * rho j n * |D2 2 (3 - 2) (n + j) n|
          ≤ (8 * (n : ℝ) ^ 3) * 1 * (12 * n) := by
            apply mul_le_mul (mul_le_mul hc r1 r0 (by positivity)) hD' (abs_nonneg _) (by positivity)
        _ = 96 * (n : ℝ) ^ 4 := by ring
    · simp [D2]
  unfold L2
  rw [abs_mul, abs_neg]
  have hsum : |∑ j ∈ Finset.range (3 + 1), ∑ i ∈ Finset.range (3 + 1),
      (ecoef 3 n i (3 - i)).coeff j * rho j n * D2 i (3 - i) (n + j) n| ≤ 16 * (96 * (n : ℝ) ^ 4) := by
    refine (Finset.abs_sum_le_sum_abs _ _).trans ?_
    calc ∑ j ∈ Finset.range (3 + 1), |∑ i ∈ Finset.range (3 + 1),
          (ecoef 3 n i (3 - i)).coeff j * rho j n * D2 i (3 - i) (n + j) n|
        ≤ ∑ j ∈ Finset.range 4, ∑ i ∈ Finset.range 4, 96 * (n : ℝ) ^ 4 := by
          refine Finset.sum_le_sum fun j hj => (Finset.abs_sum_le_sum_abs _ _).trans ?_
          exact Finset.sum_le_sum fun i hi => hterm j hj i hi
      _ = 16 * (96 * (n : ℝ) ^ 4) := by simp; ring
  calc _ ≤ 1 * (16 * (96 * (n : ℝ) ^ 4)) :=
        mul_le_mul hpre hsum (abs_nonneg _) zero_le_one
    _ = _ := by ring

end V05

theorem V05.proof_12 : V05.stmt_12 := by
  refine ⟨fun n => V05.PQ 2 n 0 0, fun n => V05.PQ 3 n 0 0,
    fun n => V05.PQ 3 n 1 0 - V05.PQ 3 n 0 0, fun n hn => ?_, fun α hα0 hα1 => ?_⟩
  · refine ⟨V05.PQ_cast00 2 n, V05.PQ_cast00 3 n, ?_⟩
    push_cast; rw [V05.PQ_cast00, V05.PQ_cast10]
  obtain ⟨hc0, hc1, -, -, -⟩ := D04.cs_facts hα0 hα1
  obtain ⟨δ, hδ0, -, hcd, -⟩ := V05.proof_4 α hα0 hα1
  have h0 := V05.proof_8 α hα0 hα1 0 (by norm_num)
  have h1 := V05.proof_8 α hα0 hα1 1 (by norm_num)
  have h2 := V05.proof_8 α hα0 hα1 2 (by norm_num)
  have h3 := V05.proof_8 α hα0 hα1 3 (by norm_num)
  have hPhi := V05.proof_1 (D04.c α) (by rw [abs_of_pos hc0]; exact hc0)
    (by rw [abs_of_pos hc0]; exact hc1)
  have e2 : ∀ n : ℕ, V05.P 2 (n + 1) V05.z2 V05.z3 =
      ((V05.PQ 2 (n + 1) 0 0 : ℚ) : ℝ) + V05.z2 / (4 * ((n : ℝ) + 1) ^ 2 - 1) := by
    intro n
    rw [V05.P_lin, V05.L2_two (by omega), V05.L3_two, V05.PQ_cast00]; push_cast; ring
  have e3 : ∀ n : ℕ, V05.P 3 (n + 1) V05.z2 V05.z3 =
      ((V05.PQ 3 (n + 1) 0 0 : ℚ) : ℝ) + (((V05.PQ 3 (n + 1) 1 0 - V05.PQ 3 (n + 1) 0 0 : ℚ) : ℝ))
        * V05.z2 - V05.z3 / (4 * ((n : ℝ) + 1) ^ 2 - 1) := by
    intro n
    push_cast
    rw [V05.P_lin, V05.L3_three (by omega), V05.PQ_cast00, V05.PQ_cast10, V05.P3_x2]; push_cast; ring
  have hM2 : HasSum (fun n : ℕ => (((V05.PQ 2 (n + 1) 0 0 : ℚ) : ℝ)
      + V05.z2 / (4 * ((n : ℝ) + 1) ^ 2 - 1)) * D04.c α ^ (2 * (n + 1))) (V05.Mk 2 α) := by
    refine h2.congr_fun ?_
    intro n; rw [e2]
  have hM3 : HasSum (fun n : ℕ => (((V05.PQ 3 (n + 1) 0 0 : ℚ) : ℝ)
      + (((V05.PQ 3 (n + 1) 1 0 - V05.PQ 3 (n + 1) 0 0 : ℚ) : ℝ)) * V05.z2
      - V05.z3 / (4 * ((n : ℝ) + 1) ^ 2 - 1)) * D04.c α ^ (2 * (n + 1))) (V05.Mk 3 α) := by
    refine h3.congr_fun ?_
    intro n; rw [e3]
  -- summability of the q'-series
  have hq : Summable (fun n : ℕ => (((V05.PQ 3 (n + 1) 1 0 - V05.PQ 3 (n + 1) 0 0 : ℚ) : ℝ))
      * D04.c α ^ (2 * (n + 1))) := by
    have hr0 : 0 ≤ D04.c α ^ 2 := by positivity
    have hr1 : D04.c α ^ 2 < 1 := by nlinarith
    refine (V05.maj_summable 2 hr0 hr1 1536).of_norm_bounded fun n => ?_
    have hq' : (((V05.PQ 3 (n + 1) 1 0 - V05.PQ 3 (n + 1) 0 0 : ℚ) : ℝ)) = V05.L2 3 (n + 1) := by
      push_cast; rw [V05.PQ_cast00, V05.PQ_cast10, V05.P3_x2]
    rw [hq', Real.norm_eq_abs, abs_mul, abs_of_nonneg (by positivity : (0:ℝ) ≤ D04.c α ^ (2 * (n + 1))),
      pow_mul]
    have hb := V05.L23_bd (n := n + 1) (by omega)
    push_cast at hb
    exact mul_le_mul_of_nonneg_right (by rw [show 2 * 2 = 4 from rfl]; exact hb) (by positivity)
  refine ⟨fun k => (hcd k).contDiffAt (Ioo_mem_nhds (by linarith) hδ0), ?_, ?_, hM2, hM3, ?_, ?_⟩
  · refine h0.congr_fun ?_
    intro n; rw [V05.P0_eq]; push_cast; ring
  · refine h1.congr_fun ?_
    intro n; rw [V05.P1_eq (by omega)]; push_cast; rfl
  · have := hM2.sub (hPhi.mul_left V05.z2)
    refine this.congr_fun ?_
    intro n
    have h4 : (4 * ((n : ℝ) + 1) ^ 2 - 1) ≠ 0 := by
      have : (0 : ℝ) ≤ n := Nat.cast_nonneg n
      nlinarith
    rw [pow_mul] ; field_simp; ring
  · refine ⟨∑' n, (((V05.PQ 3 (n + 1) 1 0 - V05.PQ 3 (n + 1) 0 0 : ℚ) : ℝ))
      * D04.c α ^ (2 * (n + 1)), hq.hasSum, ?_⟩
    have := (hM3.sub (hq.hasSum.mul_left V05.z2)).add (hPhi.mul_left V05.z3)
    rw [show V05.Mk 3 α + V05.z3 * V05.Phi (D04.c α) - V05.z2 * ∑' n, (((V05.PQ 3 (n + 1) 1 0
      - V05.PQ 3 (n + 1) 0 0 : ℚ) : ℝ)) * D04.c α ^ (2 * (n + 1)) = V05.Mk 3 α - V05.z2 * ∑' n,
      (((V05.PQ 3 (n + 1) 1 0 - V05.PQ 3 (n + 1) 0 0 : ℚ) : ℝ)) * D04.c α ^ (2 * (n + 1))
      + V05.z3 * V05.Phi (D04.c α) by ring]
    refine this.congr_fun ?_
    intro n
    have h4 : (4 * ((n : ℝ) + 1) ^ 2 - 1) ≠ 0 := by
      have : (0 : ℝ) ≤ n := Nat.cast_nonneg n
      nlinarith
    dsimp only
    rw [pow_mul]; field_simp; ring


