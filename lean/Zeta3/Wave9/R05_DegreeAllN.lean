import Mathlib
import Zeta3.Wave8.Q03_Theta7
/-!
# R05 — `prop:degree` for every `n ≥ 2`: the degree and leading coefficient of `Θ_n`

Target id: R05.
Paper label: `prop:degree` (main.tex lines 1214–1216).
Detailed proofs: paper/ieee/proofs/sec_g2_degree.tex (Proposition `g2:prop:degree` line 11,
Lemmas `g2:lem:d:coef` line 20, `g2:lem:d:e` line 37, `g2:lem:d:theta` line 51,
Definition `g2:def:amzv` line 66, Lemma `g2:lem:d:atoms` line 74, proof of Proposition 2) and
paper/ieee/proofs/sec_g2_fourier.tex (Lemma `g2:lem:f:red` line 34, Lemma `g2:lem:f:watson`
line 76 for every `N`, Lemma `g2:lem:f:atoms` line 124).
Earlier Lean: `Zeta3.Wave5.Y04_Fourier` (`h_p`, `a_k`, `ℓ`, `θ_x`, `G_k`, `I_{m,j}`, Watson at
`N = 4`), `Zeta3.Wave6.Z02_Theta5` (`N = 5`, `Θ₂…Θ₅`), `Zeta3.Wave7.V02_Theta6` (`N = 6`, `Θ₆`),
`Zeta3.Wave8.Q03_Theta7` (`N = 7`, `Θ₇`). Those objects are imported, not redefined. The
numbers `e_{r,m}`, the polynomials `ρ_m`, `Θ_n`, the Gudermannian `gd` and `H` were never defined
as Lean objects in the earlier files (they only appeared as explicit numbers), so they are
defined here, by the supplement's formulas.

Quoted verbatim.

main.tex 1214–1216 (`prop:degree`):

  For $n\ge2$, the coefficient of $\tau^n$ in $(-1)^k\frac{\pi}{\sqrt2}a_k(h_p)$ is a polynomial
  $\Theta_n(k^2)$ of degree exactly $\lfloor n/2\rfloor-1$, whose coefficients are alternating
  multiple zeta values of weight $n$. For even $n=2m$ its leading coefficient is
  $(-1)^{m-1}(2^{2m-1}-1)\zeta(2m)$.

sec_g2_fourier.tex, Lemma `g2:lem:f:watson`:

  Let $G:[0,\infty)\to\RR$ be bounded and of class $C^\infty$ on $[0,\infty)$, with
  $g_m=G^{(m)}(0)/m!$. For every integer $N\ge2$ there is $C$ (depending on $G$ and $N$) such
  that for $0<\tau<1$
  \Bigl|\,\tau\int_0^\infty\bigl[(1+\eu^{-s})^\tau-1\bigr]G(\tau s)\,ds-\sum_{n=2}^{N}\tau^n
  \sum_{j=1}^{n-1}\frac{g_{n-1-j}}{j!}\,I_{n-1-j,\,j}\Bigr|\le C\tau^{N+1}.

sec_g2_fourier.tex, Lemma `g2:lem:f:red`(c) (the part not already in Y04):

  Moreover $\gd(0)=0$, $\gd'(x)=1/\cosh x$, and $G_k(x)=(-1)^k\cos\bigl(2k\gd(x)\bigr)\,H(x)$.
  (with $\gd(x)=\frac\pi2-2\theta_x$, $H(x)=\eu^{-x}\bigl(1+\eu^{-2x}\bigr)^{-3/2}$)

sec_g2_fourier.tex, Lemma `g2:lem:f:atoms` (general part):

  More generally $I_{m,1}=m!\,\eta(m+2)$ with $\eta(s)=\sum_{l\ge1}(-1)^{l+1}l^{-s}=(1-2^{1-s})\zeta(s)$.

sec_g2_degree.tex, Lemma `g2:lem:d:coef`:

  For $k\ge0$, $N\ge2$, as $\tau\to0^+$,
  a_k(h_p)=A_k+\sum_{n=2}^Nc_n(k)\tau^n+O_{k,N}(\tau^{N+1}),\qquad
  c_n(k)=\frac8\pi\sum_{j=1}^{n-1}\frac{g_{k,n-1-j}}{j!}\,I_{n-1-j,\,j},
  with $g_{k,m}=G_k^{(m)}(0)/m!$ as in Section~\ref{g2:sec:fourier}.

sec_g2_degree.tex, before Lemma `g2:lem:d:e`:

  For integers $r,m\ge0$ put $e_{r,m}=2^{3/2}\,[y^m]\bigl(\gd(y)^{2r}H(y)\bigr)$, and for $m\ge0$
  put the polynomial
  \rho_m(u)=\sum_{r=0}^{\lfloor m/2\rfloor}\frac{(-4u)^r}{(2r)!}\,e_{r,m}.

sec_g2_degree.tex, Lemma `g2:lem:d:e`:

  (a) All $e_{r,m}$ are rational. $e_{r,m}=0$ for $m<2r$, $e_{r,2r}=1$ and $e_{r,2r+1}=\frac12$.
  (b) $g_{k,m}=(-1)^k\frac{\sqrt2}4\,\rho_m(k^2)$ for all $k,m\ge0$.

sec_g2_degree.tex, Lemma `g2:lem:d:theta`:

  For every $n\ge2$ the coefficient of $\tau^n$ in $(-1)^k\frac\pi{\sqrt2}a_k(h_p)$ is
  $\Theta_n(k^2)$, where
  \Theta_n(u)=2\sum_{j=1}^{n-1}\frac{\rho_{n-1-j}(u)}{j!}\,I_{n-1-j,\,j}.

sec_g2_degree.tex, proof of Proposition 2, case $n=2M+1$:

  The coefficient of $u^{M-1}$ is
  \frac{(-4)^{M-1}}{(2M-2)!}\bigl(I_{2M-1,1}+I_{2M-2,2}\bigr).
  Both integrals are positive (positive integrands), so the coefficient is not $0$ and has the
  sign $(-1)^{M-1}$.

sec_g2_degree.tex, Definition `g2:def:amzv` and Lemma `g2:lem:d:atoms`:

  A \emph{word} of length $w$ is a sequence $(a_1,\dots,a_w)$ with $a_i\in\{0,1,-1\}$. It is
  \emph{convergent} if $a_1\ne1$ and $a_w\ne0$. For a convergent word put
  \Gdim(a_1,\dots,a_w)=\int_{1>t_1>\dots>t_w>0}\ \prod_{i=1}^w\frac{dt_i}{t_i-a_i}.
  For $m\ge0$ and $j\ge1$,
  I_{m,j}=m!\,j!\,\Gdim\bigl(\underbrace{0,\dots,0}_{m+1},\underbrace{-1,\dots,-1}_{j}\bigr).

Conventions.
* `ζ(2m)` is `(riemannZeta (2m)).re` (Mathlib's complex zeta at the natural number `2m`, real
  part); in the consistency statement the values `ζ(2), ζ(4), ζ(6)` are written `π²/6, π⁴/90,
  π⁶/945`, and `ζ(3), ζ(5), ζ(7)` are `Y04.z3`, `Z02.z5`, `Q03.z7` (each the real part of
  `riemannZeta`), as in the earlier files.
* `[y^m] f` is written `iteratedDeriv m f 0 / m!`; the functions involved are `C^∞` (this is
  part of `stmt_4`), so this is the Taylor coefficient, not a junk value.
* `2^{3/2}` is written `2 √2`.
* `g_{k,m}` is `iteratedDeriv m (Y04.G k) 0 / m!`, `A_k` is `Y04.a k Y04.hDiamond`.
* Watson's lemma (`stmt_1`) is stated, as in Y04/Z02/V02/Q03, for `G` that is `C^∞` on all of `ℝ`
  and bounded on `[0,∞)`, so that `G^{(m)}(0)` is the ordinary derivative; the paper says
  "`C^∞` on `[0,∞)`" (one-sided derivatives at `0`). Every such function extends to a `C^∞`
  function on `ℝ` (Seeley/Borel), so this is the same class; it is the form the earlier files
  use and the only one needed (for `G = G_k`).
* Integrability of the integrals that occur is stated explicitly.

NOT COVERED.
* "whose coefficients are alternating multiple zeta values of weight `n`" is covered in two
  halves: `stmt_6` says each coefficient of `Θ_n` is a rational combination of the `I_{r,j}`
  with `r + j = n - 1`, and `stmt_10` (Lemma `g2:lem:d:atoms`) says each such `I_{r,j}` is
  `r! j!` times the iterated integral of the convergent word `(0^{r+1}, (-1)^j)`. The general
  notion "an alternating MZV of weight `w`" (any convergent word over `{0,1,-1}`) is not defined
  as a predicate; only the specific words that occur are.
* The remark on the odd case (identification of the top coefficient with
  `(-4)^{M-1}[ζ(2M+1) + 2Σ η(2j)ζ(2M+1-2j)]`, not needed by the proof).
* Uniqueness of asymptotic coefficients (used only to make "the coefficient of `τ^n`"
  meaningful; the statements give the expansion itself).
-/
-- BEGIN STATEMENTS
namespace R05

open Real MeasureTheory Set Polynomial

/-- The Gudermannian function in the supplement's normalisation, `gd(x) = π/2 - 2θ_x`,
`θ_x = arctan e^{-x}`. -/
noncomputable def gd (x : ℝ) : ℝ := π / 2 - 2 * Y04.thetaX x

/-- `H(x) = e^{-x} (1 + e^{-2x})^{-3/2}` (positive base, so `rpow` is the usual power). -/
noncomputable def H (x : ℝ) : ℝ := Real.exp (-x) * (1 + Real.exp (-2 * x)) ^ (-(3 : ℝ) / 2)

/-- `g_{k,m} = G_k^{(m)}(0)/m!`. -/
noncomputable def g (k m : ℕ) : ℝ := iteratedDeriv m (Y04.G k) 0 / (m.factorial : ℝ)

/-- `e_{r,m} = 2^{3/2} [y^m](gd(y)^{2r} H(y))`. -/
noncomputable def e (r m : ℕ) : ℝ :=
  2 * Real.sqrt 2 *
    (iteratedDeriv m (fun y : ℝ => gd y ^ (2 * r) * H y) 0 / (m.factorial : ℝ))

/-- `ρ_m(u) = Σ_{r=0}^{⌊m/2⌋} (-4u)^r/(2r)! e_{r,m}`, a polynomial in `u`. -/
noncomputable def rho (m : ℕ) : ℝ[X] :=
  ∑ r ∈ Finset.range (m / 2 + 1),
    Polynomial.C ((-4 : ℝ) ^ r / ((2 * r).factorial : ℝ) * e r m) * X ^ r

/-- `Θ_n(u) = 2 Σ_{j=1}^{n-1} ρ_{n-1-j}(u)/j! I_{n-1-j,j}` (Lemma `g2:lem:d:theta`). -/
noncomputable def Theta (n : ℕ) : ℝ[X] :=
  2 * ∑ j ∈ Finset.Icc 1 (n - 1),
    Polynomial.C (Y04.I (n - 1 - j) j / (j.factorial : ℝ)) * rho (n - 1 - j)

/-- `c_n(k) = (8/π) Σ_{j=1}^{n-1} g_{k,n-1-j}/j! I_{n-1-j,j}` (Lemma `g2:lem:d:coef`). -/
noncomputable def c (n k : ℕ) : ℝ :=
  8 / π * ∑ j ∈ Finset.Icc 1 (n - 1), g k (n - 1 - j) / (j.factorial : ℝ) * Y04.I (n - 1 - j) j

/-- The convergent word `(0^{m+1}, (-1)^j)` of length `m+1+j`, as a function of the position. -/
noncomputable def word (m j : ℕ) (i : Fin (m + 1 + j)) : ℝ := if i.val < m + 1 then 0 else -1

/-- The simplex `1 > t_1 > ⋯ > t_w > 0` in `ℝ^w` (positions `0, …, w-1`). -/
def simplex (w : ℕ) : Set (Fin w → ℝ) := {t | (∀ i, 0 < t i ∧ t i < 1) ∧ StrictAnti t}

/-- `Γ(a_1,…,a_w) = ∫_{1>t_1>⋯>t_w>0} Π dt_i/(t_i - a_i)` (Definition `g2:def:amzv`). -/
noncomputable def iterInt (w : ℕ) (a : Fin w → ℝ) : ℝ :=
  ∫ t in simplex w, ∏ i, 1 / (t i - a i)

/-- Lemma `g2:lem:f:watson` for every `N ≥ 2`: for `G` smooth and bounded on `[0,∞)`, with
`g_m = G^{(m)}(0)/m!`, there is `C` with, for `0 < τ < 1`, the integrand integrable and
`|τ ∫_0^∞ [(1+e^{-s})^τ - 1] G(τ s) ds - Σ_{n=2}^{N} τ^n Σ_{j=1}^{n-1} g_{n-1-j}/j! I_{n-1-j,j}|
  ≤ C τ^{N+1}`. -/
def stmt_1 : Prop :=
  ∀ Gf : ℝ → ℝ, ContDiff ℝ (⊤ : ℕ∞) Gf → (∃ M : ℝ, ∀ x : ℝ, 0 ≤ x → |Gf x| ≤ M) →
    ∀ N : ℕ, 2 ≤ N → ∃ C : ℝ, ∀ τ : ℝ, 0 < τ → τ < 1 →
      IntegrableOn (fun s : ℝ => ((1 + Real.exp (-s)) ^ τ - 1) * Gf (τ * s)) (Ioi 0) ∧
      |τ * (∫ s in Ioi (0 : ℝ), ((1 + Real.exp (-s)) ^ τ - 1) * Gf (τ * s)) -
          ∑ n ∈ Finset.Icc 2 N, τ ^ n * ∑ j ∈ Finset.Icc 1 (n - 1),
            iteratedDeriv (n - 1 - j) Gf 0 / ((n - 1 - j).factorial : ℝ) /
              (j.factorial : ℝ) * Y04.I (n - 1 - j) j| ≤ C * τ ^ (N + 1)

/-- Lemma `g2:lem:f:red`(c), the Gudermannian part: `gd(0) = 0`, `gd'(x) = 1/cosh x`, and
`G_k(x) = (-1)^k cos(2k gd(x)) H(x)` for all real `x` and all `k`. -/
def stmt_2 : Prop :=
  gd 0 = 0 ∧ (∀ x : ℝ, HasDerivAt gd (1 / Real.cosh x) x) ∧
    ∀ k : ℕ, ∀ x : ℝ, Y04.G k x = (-1 : ℝ) ^ k * Real.cos (2 * k * gd x) * H x

/-- Lemma `g2:lem:d:coef`: for every `k ≥ 0` and `N ≥ 2` there are `C`, `τ₀ > 0` with
`|a_k(h_p) - A_k - Σ_{n=2}^N c_n(k) τ^n| ≤ C τ^{N+1}` for `p > 1`, `τ = 1 - 1/p ≤ τ₀`. -/
def stmt_3 : Prop :=
  ∀ k N : ℕ, 2 ≤ N → ∃ C τ₀ : ℝ, 0 < τ₀ ∧ ∀ p : ℝ, 1 < p → Y04.tau p ≤ τ₀ →
    |Y04.a k (Y04.h p) - Y04.a k Y04.hDiamond - ∑ n ∈ Finset.Icc 2 N, c n k * Y04.tau p ^ n|
      ≤ C * Y04.tau p ^ (N + 1)

/-- Lemma `g2:lem:d:e`(a), with the smoothness that makes `e_{r,m}` a Taylor coefficient: each
`gd^{2r} H` is `C^∞`; every `e_{r,m}` is rational; `e_{r,m} = 0` for `m < 2r`;
`e_{r,2r} = 1`; `e_{r,2r+1} = 1/2`. -/
def stmt_4 : Prop :=
  (∀ r : ℕ, ContDiff ℝ (⊤ : ℕ∞) (fun y : ℝ => gd y ^ (2 * r) * H y)) ∧
    (∀ r m : ℕ, ∃ q : ℚ, e r m = (q : ℝ)) ∧
    (∀ r m : ℕ, m < 2 * r → e r m = 0) ∧
    (∀ r : ℕ, e r (2 * r) = 1) ∧
    (∀ r : ℕ, e r (2 * r + 1) = 1 / 2)

/-- Lemma `g2:lem:d:e`(b): `g_{k,m} = (-1)^k (√2/4) ρ_m(k²)` for all `k, m ≥ 0`. -/
def stmt_5 : Prop :=
  ∀ k m : ℕ, g k m = (-1 : ℝ) ^ k * Real.sqrt 2 / 4 * (rho m).eval ((k : ℝ) ^ 2)

/-- Lemma `g2:lem:d:theta`, with the coefficient description from the proof of Proposition 2:
for `n ≥ 2` and every `k`, `(-1)^k (π/√2) c_n(k) = Θ_n(k²)`; and every coefficient of `Θ_n` is a
rational combination of the integrals `I_{n-1-j,j}`, `1 ≤ j ≤ n-1` (so `r + j = n - 1`). -/
def stmt_6 : Prop :=
  ∀ n : ℕ, 2 ≤ n →
    (∀ k : ℕ, (-1 : ℝ) ^ k * (π / Real.sqrt 2) * c n k = (Theta n).eval ((k : ℝ) ^ 2)) ∧
    ∀ i : ℕ, ∃ q : ℕ → ℚ,
      (Theta n).coeff i = ∑ j ∈ Finset.Icc 1 (n - 1), (q j : ℝ) * Y04.I (n - 1 - j) j

/-- Lemma `g2:lem:f:atoms`, general part: `I_{m,1} = m! η(m+2)`, with
`η(m+2) = Σ_{l≥1} (-1)^{l+1} l^{-(m+2)}` (as a `HasSum`) and `η(m+2) = (1 - 2^{-(m+1)}) ζ(m+2)`,
`ζ(m+2) = (riemannZeta (m+2)).re`. -/
def stmt_7 : Prop :=
  ∀ m : ℕ,
    HasSum (fun l : ℕ => (-1 : ℝ) ^ l / ((l : ℝ) + 1) ^ (m + 2))
      (Y04.I m 1 / (m.factorial : ℝ)) ∧
    Y04.I m 1 = (m.factorial : ℝ) * (1 - (2 : ℝ)⁻¹ ^ (m + 1)) *
      (riemannZeta ((m + 2 : ℕ) : ℂ)).re

/-- Consistency with Y04, Z02, V02, Q03: with `u = k²`, the polynomials `Θ_n` defined by the
supplement's formula are `Θ₂ = ζ(2)`, `Θ₃ = ζ(3)`, `Θ₄ = -(33/16 + 7u) ζ(4)`,
`Θ₅ = -(3/4 + 4u) ζ(5) - (7/4 + 4u) ζ(2)ζ(3)`,
`Θ₆ = (4039/256 + (733/8)u + 31u²) ζ(6) - (7/8 + 2u) ζ(3)²`,
`Θ₇ = (129/16 + 46u + 16u²) ζ(7) + (141/16 + 50u + 16u²) ζ(2)ζ(5)
  + (1071/160 + (181/5)u + (56/5)u²) ζ(2)²ζ(3)`. -/
def stmt_8 : Prop :=
  Theta 2 = Polynomial.C (π ^ 2 / 6) ∧
  Theta 3 = Polynomial.C Y04.z3 ∧
  Theta 4 = Polynomial.C (-(33 / 16) * (π ^ 4 / 90)) + Polynomial.C (-7 * (π ^ 4 / 90)) * X ∧
  Theta 5 = Polynomial.C (-(3 / 4) * Z02.z5 - 7 / 4 * ((π ^ 2 / 6) * Y04.z3))
      + Polynomial.C (-4 * Z02.z5 - 4 * ((π ^ 2 / 6) * Y04.z3)) * X ∧
  Theta 6 = Polynomial.C (4039 / 256 * (π ^ 6 / 945) - 7 / 8 * Y04.z3 ^ 2)
      + Polynomial.C (733 / 8 * (π ^ 6 / 945) - 2 * Y04.z3 ^ 2) * X
      + Polynomial.C (31 * (π ^ 6 / 945)) * X ^ 2 ∧
  Theta 7 = Polynomial.C (129 / 16 * Q03.z7 + 141 / 16 * ((π ^ 2 / 6) * Z02.z5)
        + 1071 / 160 * ((π ^ 2 / 6) ^ 2 * Y04.z3))
      + Polynomial.C (46 * Q03.z7 + 50 * ((π ^ 2 / 6) * Z02.z5)
        + 181 / 5 * ((π ^ 2 / 6) ^ 2 * Y04.z3)) * X
      + Polynomial.C (16 * Q03.z7 + 16 * ((π ^ 2 / 6) * Z02.z5)
        + 56 / 5 * ((π ^ 2 / 6) ^ 2 * Y04.z3)) * X ^ 2

/-- `prop:degree` for every `n ≥ 2`: (i) `Θ_n(k²)` is the coefficient of `τ^n` in
`(-1)^k (π/√2) a_k(h_p)` in the asymptotic sense: for every `k` and `N ≥ 2` there are `C`,
`τ₀ > 0` with `|a_k(h_p) - (-1)^k (√2/π)[-4/(16k²-1) + Σ_{n=2}^N Θ_n(k²) τ^n]| ≤ C τ^{N+1}`;
(ii) `Θ_n ≠ 0` and `natDegree Θ_n = ⌊n/2⌋ - 1` (degree exactly); (iii) for `n = 2m` the
leading coefficient is `(-1)^{m-1}(2^{2m-1}-1) ζ(2m)`; (iv) for `n = 2M+1` the leading
coefficient is `(-4)^{M-1}/(2M-2)! (I_{2M-1,1} + I_{2M-2,2})`, with `I_{2M-1,1}`,
`I_{2M-2,2} > 0` (so it is non-zero, of sign `(-1)^{M-1}`). -/
def stmt_9 : Prop :=
  (∀ k N : ℕ, 2 ≤ N → ∃ C τ₀ : ℝ, 0 < τ₀ ∧ ∀ p : ℝ, 1 < p → Y04.tau p ≤ τ₀ →
    |Y04.a k (Y04.h p) - (-1 : ℝ) ^ k * (Real.sqrt 2 / π) *
        (-(4 / (16 * (k : ℝ) ^ 2 - 1)) +
          ∑ n ∈ Finset.Icc 2 N, (Theta n).eval ((k : ℝ) ^ 2) * Y04.tau p ^ n)|
      ≤ C * Y04.tau p ^ (N + 1)) ∧
  (∀ n : ℕ, 2 ≤ n → Theta n ≠ 0 ∧ (Theta n).natDegree = n / 2 - 1) ∧
  (∀ m : ℕ, 1 ≤ m → (Theta (2 * m)).leadingCoeff =
      (-1 : ℝ) ^ (m - 1) * (2 ^ (2 * m - 1) - 1) * (riemannZeta ((2 * m : ℕ) : ℂ)).re) ∧
  (∀ M : ℕ, 1 ≤ M →
    0 < Y04.I (2 * M - 1) 1 ∧ 0 < Y04.I (2 * M - 2) 2 ∧
    (Theta (2 * M + 1)).leadingCoeff =
      (-4 : ℝ) ^ (M - 1) / ((2 * M - 2).factorial : ℝ) *
        (Y04.I (2 * M - 1) 1 + Y04.I (2 * M - 2) 2))

/-- Lemma `g2:lem:d:atoms`: for `m ≥ 0`, `j ≥ 1`, the integrand of the iterated integral of the
convergent word `(0^{m+1}, (-1)^j)` is integrable on the simplex `1 > t_1 > ⋯ > t_{m+1+j} > 0`,
and `I_{m,j} = m! j! Γ(0^{m+1}, (-1)^j)`; so `I_{m,j}` is `m! j!` times an alternating multiple
zeta value of weight `m + j + 1`. -/
def stmt_10 : Prop :=
  ∀ m j : ℕ, 1 ≤ j →
    IntegrableOn (fun t : Fin (m + 1 + j) → ℝ => ∏ i, 1 / (t i - word m j i))
      (simplex (m + 1 + j)) ∧
    Y04.I m j = (m.factorial : ℝ) * (j.factorial : ℝ) * iterInt (m + 1 + j) (word m j)

end R05
-- END STATEMENTS

open Real MeasureTheory Set

namespace R05

lemma eta_m (m : ℕ) :
    HasSum (fun l : ℕ => (-1 : ℝ) ^ l / ((l : ℝ) + 1) ^ (m + 2))
      (Y04.I m 1 / (m.factorial : ℝ)) ∧
    Y04.I m 1 = (m.factorial : ℝ) * (1 - (2 : ℝ)⁻¹ ^ (m + 1)) *
      (riemannZeta ((m + 2 : ℕ) : ℂ)).re := by
  have hz := Y04.z_hasSum (k := m + 2) (by omega)
  have ha := Y04.hasSum_alt (m + 2) hz
  have hf : (fun l : ℕ => (-1 : ℝ) ^ l / ((l : ℝ) + 1) ^ (m + 2)) =
      fun n : ℕ => (-1) ^ n * (1 / ((n : ℝ) + 1) ^ (m + 2)) := by
    funext n; ring
  have hI := Y04.I_m1 m
  have hfac : (m.factorial : ℝ) ≠ 0 := by positivity
  have ht : ∑' n : ℕ, (-1 : ℝ) ^ n / ((n : ℝ) + 1) ^ (m + 2) =
      Y04.zn (m + 2) - 2 * (Y04.zn (m + 2) / 2 ^ (m + 2)) := by
    rw [hf]; exact ha.tsum_eq
  refine ⟨?_, ?_⟩
  · have e2 : Y04.I m 1 / (m.factorial : ℝ) = Y04.zn (m + 2) - 2 * (Y04.zn (m + 2) / 2 ^ (m + 2)) := by
      rw [hI, ht]; field_simp
    rw [hf, e2]; exact ha
  · rw [hI, ht]
    have : Y04.zn (m + 2) = (riemannZeta ((m + 2 : ℕ) : ℂ)).re := rfl
    rw [this, inv_pow, pow_succ]
    field_simp

lemma gd_zero : gd 0 = 0 := by
  simp [gd, Y04.thetaX_zero]; ring

lemma gd_hasDeriv (x : ℝ) : HasDerivAt gd (1 / Real.cosh x) x := by
  have h := ((Y04.theta_hasDeriv x).const_mul 2).const_sub (π / 2)
  have e : gd = fun x => π / 2 - 2 * Y04.thetaX x := rfl
  rw [e]
  refine h.congr_deriv ?_
  simp only [Y04.thetaX, Real.sin_arctan, Real.cos_arctan, Real.cosh_eq]
  have hs : Real.sqrt (1 + Real.exp (-x) ^ 2) ^ 2 = 1 + Real.exp (-x) ^ 2 :=
    Real.sq_sqrt (by positivity)
  have hp : 0 < Real.sqrt (1 + Real.exp (-x) ^ 2) := Real.sqrt_pos.2 (by positivity)
  have hx : Real.exp x = (Real.exp (-x))⁻¹ := by rw [Real.exp_neg, inv_inv]
  rw [hx]
  have ht : 0 < Real.exp (-x) := Real.exp_pos _
  field_simp
  rw [hs]

lemma cos_sub_pi (k : ℕ) (y : ℝ) : Real.cos (k * π - y) = (-1 : ℝ) ^ k * Real.cos y := by
  rw [Real.cos_sub, Real.cos_nat_mul_pi, Real.sin_nat_mul_pi]; ring

lemma H_eq_trig (x : ℝ) :
    Real.cos (Y04.thetaX x) ^ 2 * Real.sin (Y04.thetaX x) = H x := by
  simp only [Y04.thetaX, Real.sin_arctan, Real.cos_arctan, H]
  set t := Real.exp (-x) with htdef
  have ht : 0 < t := Real.exp_pos _
  have h2 : Real.exp (-2 * x) = t ^ 2 := by
    rw [htdef, ← Real.exp_nat_mul]; congr 1; push_cast; ring
  rw [h2]
  have ha : 0 < 1 + t ^ 2 := by positivity
  have hs : Real.sqrt (1 + t ^ 2) ^ 2 = 1 + t ^ 2 := Real.sq_sqrt ha.le
  have hp : 0 < Real.sqrt (1 + t ^ 2) := Real.sqrt_pos.2 ha
  have hr : (1 + t ^ 2) ^ (-(3 : ℝ) / 2) = 1 / ((1 + t ^ 2) * Real.sqrt (1 + t ^ 2)) := by
    rw [show (-(3 : ℝ) / 2) = -(1 + 1 / 2) by norm_num, Real.rpow_neg ha.le,
      Real.rpow_add ha, Real.rpow_one, ← Real.sqrt_eq_rpow, one_div]
  rw [hr]
  field_simp
  rw [hs]

lemma G_eq_gd (k : ℕ) (x : ℝ) :
    Y04.G k x = (-1 : ℝ) ^ k * Real.cos (2 * k * gd x) * H x := by
  have e : 4 * (k : ℝ) * Y04.thetaX x = k * π - 2 * k * gd x := by
    simp only [gd]; ring
  simp only [Y04.G]
  rw [e, cos_sub_pi, ← H_eq_trig]; ring

end R05

theorem R05.proof_2 : R05.stmt_2 := ⟨R05.gd_zero, R05.gd_hasDeriv, R05.G_eq_gd⟩

theorem R05.proof_7 : R05.stmt_7 := R05.eta_m

open Real MeasureTheory Set

namespace R05

/-- Taylor coefficient at `0`. -/
noncomputable def tc (f : ℝ → ℝ) (i : ℕ) : ℝ := iteratedDeriv i f 0 / (i.factorial : ℝ)

/-- The formal Taylor series at `0`. -/
noncomputable def ts (f : ℝ → ℝ) : PowerSeries ℝ := PowerSeries.mk (tc f)

lemma coeff_ts (f : ℝ → ℝ) (i : ℕ) : PowerSeries.coeff i (ts f) = tc f i := by
  simp [ts]

lemma sm_at {f : ℝ → ℝ} (hf : ContDiff ℝ (⊤ : ℕ∞) f) (m : ℕ) (x : ℝ) : ContDiffAt ℝ m f x :=
  (hf.of_le (by exact_mod_cast le_top)).contDiffAt

lemma ts_mul {f g : ℝ → ℝ} (hf : ContDiff ℝ (⊤ : ℕ∞) f) (hg : ContDiff ℝ (⊤ : ℕ∞) g) :
    ts (fun y => f y * g y) = ts f * ts g := by
  ext m
  rw [PowerSeries.coeff_mul, Finset.Nat.sum_antidiagonal_eq_sum_range_succ_mk, coeff_ts]
  simp only [coeff_ts, tc]
  rw [iteratedDeriv_fun_mul (sm_at hf m 0) (sm_at hg m 0), Finset.sum_div]
  refine Finset.sum_congr rfl (fun i hi => ?_)
  have hi' : i ≤ m := Nat.lt_succ_iff.mp (Finset.mem_range.mp hi)
  rw [Nat.cast_choose ℝ hi']
  have h1 : (i.factorial : ℝ) ≠ 0 := by positivity
  have h2 : ((m - i).factorial : ℝ) ≠ 0 := by positivity
  have h3 : (m.factorial : ℝ) ≠ 0 := by positivity
  field_simp

lemma ts_add {f g : ℝ → ℝ} (hf : ContDiff ℝ (⊤ : ℕ∞) f) (hg : ContDiff ℝ (⊤ : ℕ∞) g) :
    ts (fun y => f y + g y) = ts f + ts g := by
  ext m
  simp only [map_add, coeff_ts, tc]
  rw [iteratedDeriv_fun_add (sm_at hf m 0) (sm_at hg m 0), add_div]

lemma ts_sub {f g : ℝ → ℝ} (hf : ContDiff ℝ (⊤ : ℕ∞) f) (hg : ContDiff ℝ (⊤ : ℕ∞) g) :
    ts (fun y => f y - g y) = ts f - ts g := by
  ext m
  simp only [map_sub, coeff_ts, tc]
  rw [iteratedDeriv_fun_sub (sm_at hf m 0) (sm_at hg m 0), sub_div]

lemma ts_cmul (c : ℝ) (f : ℝ → ℝ) : ts (fun y => c * f y) = PowerSeries.C c * ts f := by
  ext m
  simp only [PowerSeries.coeff_C_mul, coeff_ts, tc]
  rw [iteratedDeriv_const_mul_field]; ring

lemma ts_monomial (i : ℕ) : ts (fun y => y ^ i) = PowerSeries.X ^ i := by
  ext m
  simp only [coeff_ts, tc, PowerSeries.coeff_X_pow]
  rw [iteratedDeriv_fun_pow_zero]
  split_ifs with h
  · subst h; field_simp
  · simp

lemma ts_const (c : ℝ) : ts (fun _ => c) = PowerSeries.C c := by
  have h1 : ts (fun _ : ℝ => (1 : ℝ)) = 1 := by simpa using ts_monomial 0
  have := ts_cmul c (fun _ => (1 : ℝ))
  simp only [mul_one] at this
  rw [this, h1, mul_one]

lemma ts_pow {f : ℝ → ℝ} (hf : ContDiff ℝ (⊤ : ℕ∞) f) (n : ℕ) :
    ts (fun y => f y ^ n) = ts f ^ n := by
  induction n with
  | zero => simp only [pow_zero]; rw [ts_const]; rfl
  | succ n ih =>
    simp only [pow_succ]
    rw [ts_mul (hf.pow n) hf, ih]

lemma ts_sum {ι : Type*} (s : Finset ι) (f : ι → ℝ → ℝ)
    (hf : ∀ i ∈ s, ContDiff ℝ (⊤ : ℕ∞) (f i)) :
    ts (fun y => ∑ i ∈ s, f i y) = ∑ i ∈ s, ts (f i) := by
  classical
  induction s using Finset.induction_on with
  | empty => simp only [Finset.sum_empty]; rw [ts_const]; simp
  | insert a s ha ih =>
    simp only [Finset.sum_insert ha]
    have hs : ContDiff ℝ (⊤ : ℕ∞) (fun y => ∑ i ∈ s, f i y) :=
      ContDiff.sum (fun i hi => hf i (Finset.mem_insert_of_mem hi))
    rw [ts_add (hf a (Finset.mem_insert_self a s)) hs,
      ih (fun i hi => hf i (Finset.mem_insert_of_mem hi))]

/-- A smooth function whose Taylor coefficients below `n` vanish is `O(|y|^n)` on `[-1,1]`. -/
lemma flat_bound (n : ℕ) : ∀ R : ℝ → ℝ, ContDiff ℝ (⊤ : ℕ∞) R →
    (∀ i < n, iteratedDeriv i R 0 = 0) → ∃ C, 0 ≤ C ∧ ∀ y, |y| ≤ 1 → |R y| ≤ C * |y| ^ n := by
  induction n with
  | zero =>
    intro R hR _
    obtain ⟨C, hC⟩ := (isCompact_Icc (a := (-1 : ℝ)) (b := 1)).exists_bound_of_continuousOn
      hR.continuous.continuousOn
    refine ⟨max C 0, le_max_right _ _, fun y hy => ?_⟩
    have := hC y (abs_le.mp hy)
    rw [Real.norm_eq_abs] at this
    simp only [pow_zero, mul_one]
    exact this.trans (le_max_left _ _)
  | succ n ih =>
    intro R hR h0
    have hd := contDiff_infty_iff_deriv.mp hR
    obtain ⟨C, hC0, hC⟩ := ih (deriv R) hd.2 (fun i hi => by
      rw [← iteratedDeriv_succ']; exact h0 (i + 1) (by omega))
    refine ⟨C, hC0, fun y hy => ?_⟩
    have hR0 : R 0 = 0 := by have := h0 0 (by omega); simpa using this
    have hb : ∀ z ∈ Set.uIcc 0 y, ‖deriv R z‖ ≤ C * |y| ^ n := by
      intro z hz
      have hzy : |z| ≤ |y| := by
        rcases Set.mem_uIcc.mp hz with ⟨h1, h2⟩ | ⟨h1, h2⟩
        · rw [abs_of_nonneg h1, abs_of_nonneg (h1.trans h2)]; exact h2
        · rw [abs_of_nonpos h2, abs_of_nonpos (h1.trans h2)]; linarith
      rw [Real.norm_eq_abs]
      calc |deriv R z| ≤ C * |z| ^ n := hC z (hzy.trans hy)
        _ ≤ C * |y| ^ n := by gcongr
    have := (convex_uIcc (0 : ℝ) y).norm_image_sub_le_of_norm_deriv_le
      (fun z _ => hd.1 z) hb Set.left_mem_uIcc Set.right_mem_uIcc
    simp only [hR0, sub_zero, Real.norm_eq_abs] at this
    calc |R y| ≤ C * |y| ^ n * |y| := this
      _ = C * |y| ^ (n + 1) := by ring

/-- Conversely, a smooth function that is `O(|y|^n)` near `0` has vanishing Taylor coefficients
below `n`. -/
lemma flat_unique (n : ℕ) (R : ℝ → ℝ) (hR : ContDiff ℝ (⊤ : ℕ∞) R) {δ C : ℝ} (hδ : 0 < δ)
    (hb : ∀ y, |y| ≤ δ → |R y| ≤ C * |y| ^ n) : ∀ i < n, iteratedDeriv i R 0 = 0 := by
  intro i
  induction i using Nat.strong_induction_on with
  | _ i IH =>
  intro hin
  set a := iteratedDeriv i R 0 with ha
  have hfi : (0 : ℝ) < i.factorial := by positivity
  set Q : ℝ → ℝ := fun y => R y - a / (i.factorial : ℝ) * y ^ i with hQ
  have hp : ContDiff ℝ (⊤ : ℕ∞) (fun y : ℝ => a / (i.factorial : ℝ) * y ^ i) := by fun_prop
  have hQs : ContDiff ℝ (⊤ : ℕ∞) Q := hR.sub hp
  have hQ0 : ∀ j < i + 1, iteratedDeriv j Q 0 = 0 := by
    intro j hj
    rw [hQ, iteratedDeriv_fun_sub (sm_at hR j 0) (sm_at hp j 0)]
    rw [iteratedDeriv_const_mul_field, iteratedDeriv_fun_pow_zero]
    rcases Nat.lt_succ_iff_lt_or_eq.mp hj with h | h
    · rw [IH j h (by omega)]; simp [show j ≠ i by omega]
    · subst h; simp only [ite_true, ← ha]; field_simp; ring
  obtain ⟨C', hC'0, hC'⟩ := flat_bound (i + 1) Q hQs hQ0
  by_contra hne
  set b := |a| / (i.factorial : ℝ) with hbdef
  have hb0 : 0 < b := div_pos (abs_pos.mpr hne) hfi
  set K := |C| + C' + 1 with hK
  have hK0 : 0 < K := by positivity
  set y := min (min δ 1) (b / (2 * K)) with hy
  have hy0 : 0 < y := lt_min (lt_min hδ one_pos) (by positivity)
  have hyδ : y ≤ δ := (min_le_left _ _).trans (min_le_left _ _)
  have hy1 : y ≤ 1 := (min_le_left _ _).trans (min_le_right _ _)
  have hyb : y ≤ b / (2 * K) := min_le_right _ _
  have hay : |y| = y := abs_of_pos hy0
  have e1 := hb y (by rw [hay]; exact hyδ)
  have e2 := hC' y (by rw [hay]; exact hy1)
  rw [hay] at e1 e2
  have hpow : y ^ n ≤ y ^ (i + 1) := pow_le_pow_of_le_one hy0.le hy1 (by omega)
  have hsplit : a / (i.factorial : ℝ) * y ^ i = R y - Q y := by simp only [hQ]; ring
  have habs : b * y ^ i = |R y - Q y| := by
    rw [← hsplit, abs_mul, abs_div, abs_of_pos hfi, abs_of_nonneg (pow_nonneg hy0.le i)]
  have hyi : 0 < y ^ i := pow_pos hy0 i
  have h3 : b * y ^ i ≤ K * y ^ (i + 1) := by
    rw [habs]
    calc |R y - Q y| ≤ |R y| + |Q y| := abs_sub _ _
      _ ≤ |C| * y ^ (i + 1) + C' * y ^ (i + 1) := by
          have : C * y ^ n ≤ |C| * y ^ (i + 1) :=
            (mul_le_mul_of_nonneg_right (le_abs_self C) (pow_nonneg hy0.le n)).trans
              (mul_le_mul_of_nonneg_left hpow (abs_nonneg C))
          linarith
      _ ≤ K * y ^ (i + 1) := by
          have := pow_nonneg hy0.le (i + 1); rw [hK]; nlinarith
  have h4 : b ≤ K * y := by
    rw [pow_succ] at h3
    have : b * y ^ i ≤ (K * y) * y ^ i := by linarith
    exact le_of_mul_le_mul_right this hyi
  have h5 : K * y ≤ b / 2 := by
    calc K * y ≤ K * (b / (2 * K)) := mul_le_mul_of_nonneg_left hyb hK0.le
      _ = b / 2 := by field_simp
  linarith

end R05

namespace R05

/-- Global Taylor bound on `[0,∞)` for a smooth function bounded there. -/
lemma taylor_glob (Gf : ℝ → ℝ) (hG : ContDiff ℝ (⊤ : ℕ∞) Gf) {M : ℝ}
    (hM : ∀ x : ℝ, 0 ≤ x → |Gf x| ≤ M) (n : ℕ) :
    ∃ C, 0 ≤ C ∧ ∀ v : ℝ, 0 ≤ v → |Gf v - ∑ i ∈ Finset.range n, tc Gf i * v ^ i| ≤ C * v ^ n := by
  set Q : ℝ → ℝ := fun y => Gf y - ∑ i ∈ Finset.range n, tc Gf i * y ^ i with hQ
  have hP : ContDiff ℝ (⊤ : ℕ∞) (fun y : ℝ => ∑ i ∈ Finset.range n, tc Gf i * y ^ i) := by
    fun_prop
  have hQs : ContDiff ℝ (⊤ : ℕ∞) Q := hG.sub hP
  have hts : ts Q = ts Gf - ∑ i ∈ Finset.range n, PowerSeries.C (tc Gf i) * PowerSeries.X ^ i := by
    rw [hQ, ts_sub hG hP, ts_sum _ _ (fun i _ => by fun_prop)]
    congr 1
    refine Finset.sum_congr rfl (fun i _ => ?_)
    rw [ts_cmul, ts_monomial]
  have hQ0 : ∀ j < n, iteratedDeriv j Q 0 = 0 := by
    intro j hj
    have hc : PowerSeries.coeff j (ts Q) = 0 := by
      rw [hts, map_sub, map_sum, coeff_ts]
      simp only [PowerSeries.coeff_C_mul, PowerSeries.coeff_X_pow, mul_ite, mul_one, mul_zero]
      rw [Finset.sum_ite_eq]
      simp [Finset.mem_range.mpr hj]
    rw [coeff_ts, tc, div_eq_zero_iff] at hc
    rcases hc with hc | hc
    · exact hc
    · exact absurd hc (by positivity)
  obtain ⟨C1, hC10, hC1⟩ := flat_bound n Q hQs hQ0
  have hM0 : 0 ≤ M := (abs_nonneg _).trans (hM 0 le_rfl)
  refine ⟨C1 + M + ∑ i ∈ Finset.range n, |tc Gf i|, by
    have := Finset.sum_nonneg (fun i (_ : i ∈ Finset.range n) => abs_nonneg (tc Gf i))
    linarith, fun v hv => ?_⟩
  have hS0 : 0 ≤ ∑ i ∈ Finset.range n, |tc Gf i| :=
    Finset.sum_nonneg (fun i _ => abs_nonneg _)
  rcases le_or_gt v 1 with h1 | h1
  · have := hC1 v (by rw [abs_of_nonneg hv]; exact h1)
    rw [abs_of_nonneg hv] at this
    have hvn := pow_nonneg hv n
    calc |Gf v - ∑ i ∈ Finset.range n, tc Gf i * v ^ i| = |Q v| := rfl
      _ ≤ C1 * v ^ n := this
      _ ≤ _ := by nlinarith
  · have hvn : 1 ≤ v ^ n := one_le_pow₀ h1.le
    have hsum : |∑ i ∈ Finset.range n, tc Gf i * v ^ i| ≤ (∑ i ∈ Finset.range n, |tc Gf i|) * v ^ n := by
      rw [Finset.sum_mul]
      refine (Finset.abs_sum_le_sum_abs _ _).trans (Finset.sum_le_sum (fun i hi => ?_))
      rw [abs_mul, abs_of_nonneg (pow_nonneg hv i)]
      exact mul_le_mul_of_nonneg_left
        (pow_le_pow_right₀ h1.le (Finset.mem_range.mp hi).le) (abs_nonneg _)
    calc |Gf v - ∑ i ∈ Finset.range n, tc Gf i * v ^ i|
        ≤ |Gf v| + |∑ i ∈ Finset.range n, tc Gf i * v ^ i| := abs_sub _ _
      _ ≤ M * v ^ n + (∑ i ∈ Finset.range n, |tc Gf i|) * v ^ n := by
          have := hM v hv; nlinarith
      _ ≤ _ := by have := pow_nonneg hv n; nlinarith

lemma reindex (N : ℕ) (F : ℕ → ℕ → ℝ) :
    ∑ n ∈ Finset.Icc 2 N, ∑ j ∈ Finset.Icc 1 (n - 1), F n j =
      ∑ j ∈ Finset.Icc 1 (N - 1), ∑ i ∈ Finset.range (N - j), F (i + j + 1) j := by
  rw [Finset.sum_comm' (t' := Finset.Icc 1 (N - 1)) (s' := fun j => Finset.Icc (j + 1) N)]
  · refine Finset.sum_congr rfl (fun j hj => ?_)
    refine Finset.sum_bij' (fun n _ => n - j - 1) (fun i _ => i + j + 1) ?_ ?_ ?_ ?_ ?_
    · intro n hn; simp only [Finset.mem_Icc, Finset.mem_range] at hn ⊢; omega
    · intro i hi; simp only [Finset.mem_Icc, Finset.mem_range] at hi hj ⊢; omega
    · intro n hn; simp only [Finset.mem_Icc] at hn; omega
    · intro i _; omega
    · intro n hn; simp only [Finset.mem_Icc] at hn; congr 1; omega
  · intro n j; simp only [Finset.mem_Icc]; omega

lemma ell_le1 {s : ℝ} (hs : 0 ≤ s) : Y04.ell s ≤ 1 := Y04.ell_le_one s hs

/-- Lemma `g2:lem:f:watson` for every `N ≥ 2`. -/
lemma watsonN (Gf : ℝ → ℝ) (hG : ContDiff ℝ (⊤ : ℕ∞) Gf) (hb : ∃ M : ℝ, ∀ x : ℝ, 0 ≤ x → |Gf x| ≤ M)
    (N : ℕ) (hN : 2 ≤ N) : ∃ C : ℝ, ∀ τ : ℝ, 0 < τ → τ < 1 →
      IntegrableOn (fun s : ℝ => ((1 + Real.exp (-s)) ^ τ - 1) * Gf (τ * s)) (Ioi 0) ∧
      |τ * (∫ s in Ioi (0 : ℝ), ((1 + Real.exp (-s)) ^ τ - 1) * Gf (τ * s)) -
          ∑ n ∈ Finset.Icc 2 N, τ ^ n * ∑ j ∈ Finset.Icc 1 (n - 1),
            iteratedDeriv (n - 1 - j) Gf 0 / ((n - 1 - j).factorial : ℝ) /
              (j.factorial : ℝ) * Y04.I (n - 1 - j) j| ≤ C * τ ^ (N + 1) := by
  obtain ⟨M, hM⟩ := hb
  have hM0 : 0 ≤ M := (abs_nonneg _).trans (hM 0 le_rfl)
  obtain ⟨_, hW⟩ := Y04.watson Gf hG ⟨M, hM⟩
  have hT := fun m => taylor_glob Gf hG hM m
  choose Cf hCf0 hCf using hT
  set B : ℝ → ℝ := fun s => 2 * M * (s ^ 0 * Y04.ell s ^ N) +
    ∑ j ∈ Finset.Icc 1 (N - 1), Cf (N - j) * (s ^ (N - j) * Y04.ell s ^ j) with hB
  have hBint : IntegrableOn B (Ioi 0) := by
    refine ((Y04.I_int 0 N (by omega)).const_mul _).add ?_
    exact integrable_finsetSum _ (fun j hj => (Y04.I_int _ j (Finset.mem_Icc.mp hj).1).const_mul _)
  refine ⟨∫ s in Ioi (0 : ℝ), B s, fun τ hτ0 hτ1 => ⟨(hW τ hτ0 hτ1).1, ?_⟩⟩
  have hFint := (hW τ hτ0 hτ1).1
  set P : ℝ → ℝ := fun s => ∑ n ∈ Finset.Icc 2 N, τ ^ n * ∑ j ∈ Finset.Icc 1 (n - 1),
      (tc Gf (n - 1 - j) / (j.factorial : ℝ) * (s ^ (n - 1 - j) * Y04.ell s ^ j)) with hP
  have hPint : ∫ s in Ioi (0 : ℝ), P s = ∑ n ∈ Finset.Icc 2 N, τ ^ n * ∑ j ∈ Finset.Icc 1 (n - 1),
      iteratedDeriv (n - 1 - j) Gf 0 / ((n - 1 - j).factorial : ℝ) /
        (j.factorial : ℝ) * Y04.I (n - 1 - j) j := by
    rw [hP, integral_finsetSum]
    · refine Finset.sum_congr rfl (fun n _ => ?_)
      rw [integral_const_mul, integral_finsetSum]
      · congr 1
        refine Finset.sum_congr rfl (fun j _ => ?_)
        rw [integral_const_mul]; rfl
      · intro j hj
        exact (Y04.I_int _ j (Finset.mem_Icc.mp hj).1).const_mul _
    · intro n _
      refine Integrable.const_mul ?_ _
      exact integrable_finsetSum _ (fun j hj => (Y04.I_int _ j (Finset.mem_Icc.mp hj).1).const_mul _)
  have hPi : IntegrableOn P (Ioi 0) := by
    refine integrable_finsetSum _ (fun n _ => Integrable.const_mul ?_ _)
    exact integrable_finsetSum _ (fun j hj => (Y04.I_int _ j (Finset.mem_Icc.mp hj).1).const_mul _)
  rw [← hPint, ← integral_const_mul, ← integral_sub (hFint.const_mul τ) hPi]
  have key : ∀ s : ℝ, 0 < s →
      ‖τ * (((1 + Real.exp (-s)) ^ τ - 1) * Gf (τ * s)) - P s‖ ≤ τ ^ (N + 1) * B s := by
    intro s hs
    set u := τ * Y04.ell s with hu
    set v := τ * s with hv
    have hl0 := Y04.ell_nonneg s
    have hl1 := ell_le1 hs.le
    have hu0 : 0 ≤ u := mul_nonneg hτ0.le hl0
    have hu1 : u ≤ 1 := by rw [hu]; nlinarith
    have hv0 : 0 ≤ v := mul_nonneg hτ0.le hs.le
    -- rewrite P
    have hP' : P s = τ * ∑ j ∈ Finset.Icc 1 (N - 1), u ^ j / (j.factorial : ℝ) *
        ∑ i ∈ Finset.range (N - j), tc Gf i * v ^ i := by
      rw [hP]
      simp only [Finset.mul_sum]
      rw [reindex N (fun n j => τ ^ n *
        (tc Gf (n - 1 - j) / (j.factorial : ℝ) * (s ^ (n - 1 - j) * Y04.ell s ^ j)))]
      refine Finset.sum_congr rfl (fun j _ => Finset.sum_congr rfl (fun i _ => ?_))
      rw [show i + j + 1 - 1 - j = i by omega, hu, hv]
      ring
    have hexp : (1 + Real.exp (-s)) ^ τ = Real.exp u := by rw [Y04.rpow_eq_exp]
    set A := Real.exp u - ∑ m ∈ Finset.range N, u ^ m / (m.factorial : ℝ) with hA
    have hsplit : ∑ m ∈ Finset.range N, u ^ m / (m.factorial : ℝ) =
        1 + ∑ j ∈ Finset.Icc 1 (N - 1), u ^ j / (j.factorial : ℝ) := by
      rw [show Finset.Icc 1 (N - 1) = Finset.Ico 1 N by ext x; simp; omega,
        Finset.sum_Ico_eq_sub _ (by omega : 1 ≤ N)]
      simp
    have hid : τ * (((1 + Real.exp (-s)) ^ τ - 1) * Gf (τ * s)) - P s =
        τ * (A * Gf v + ∑ j ∈ Finset.Icc 1 (N - 1), u ^ j / (j.factorial : ℝ) *
          (Gf v - ∑ i ∈ Finset.range (N - j), tc Gf i * v ^ i)) := by
      rw [hP', hexp, hA, hsplit]
      simp only [mul_sub, Finset.sum_sub_distrib, ← Finset.sum_mul]
      ring
    rw [hid, Real.norm_eq_abs, abs_mul, abs_of_pos hτ0]
    have hA' : |A| ≤ 2 * u ^ N := by
      have := Real.exp_bound (x := u) (by rw [abs_of_nonneg hu0]; exact hu1) (n := N) (by omega)
      rw [abs_of_nonneg hu0] at this
      have hfr : ((N.succ : ℕ) : ℝ) / ((N.factorial : ℝ) * N) ≤ 2 := by
        have h1 : (1 : ℝ) ≤ N.factorial := by exact_mod_cast Nat.one_le_iff_ne_zero.mpr (Nat.factorial_ne_zero N)
        have h2 : (2 : ℝ) ≤ N := by exact_mod_cast hN
        rw [div_le_iff₀ (by positivity)]
        push_cast
        nlinarith
      calc |A| ≤ u ^ N * ((N.succ : ℕ) / ((N.factorial : ℝ) * N)) := this
        _ ≤ u ^ N * 2 := mul_le_mul_of_nonneg_left hfr (pow_nonneg hu0 N)
        _ = 2 * u ^ N := by ring
    have hterm : ∀ j ∈ Finset.Icc 1 (N - 1), |u ^ j / (j.factorial : ℝ) *
        (Gf v - ∑ i ∈ Finset.range (N - j), tc Gf i * v ^ i)| ≤
          τ ^ N * (Cf (N - j) * (s ^ (N - j) * Y04.ell s ^ j)) := by
      intro j hj
      have hjN : j ≤ N := by simp only [Finset.mem_Icc] at hj; omega
      have hf1 : (1 : ℝ) ≤ j.factorial := by exact_mod_cast Nat.one_le_iff_ne_zero.mpr (Nat.factorial_ne_zero j)
      rw [abs_mul, abs_div, abs_of_nonneg (pow_nonneg hu0 j), abs_of_pos (by positivity : (0:ℝ) < j.factorial)]
      have h2 := hCf (N - j) v hv0
      have hcj : u ^ j / (j.factorial : ℝ) ≤ u ^ j := div_le_self (pow_nonneg hu0 j) hf1
      calc u ^ j / (j.factorial : ℝ) * |Gf v - ∑ i ∈ Finset.range (N - j), tc Gf i * v ^ i|
          ≤ u ^ j * (Cf (N - j) * v ^ (N - j)) :=
            mul_le_mul hcj h2 (abs_nonneg _) (pow_nonneg hu0 j)
        _ = τ ^ (j + (N - j)) * (Cf (N - j) * (s ^ (N - j) * Y04.ell s ^ j)) := by
            rw [hu, hv]; ring
        _ = τ ^ N * (Cf (N - j) * (s ^ (N - j) * Y04.ell s ^ j)) := by
            rw [show j + (N - j) = N by omega]
    have hGv := hM v hv0
    calc τ * |A * Gf v + ∑ j ∈ Finset.Icc 1 (N - 1), u ^ j / (j.factorial : ℝ) *
          (Gf v - ∑ i ∈ Finset.range (N - j), tc Gf i * v ^ i)|
        ≤ τ * (|A| * |Gf v| + ∑ j ∈ Finset.Icc 1 (N - 1), τ ^ N * (Cf (N - j) * (s ^ (N - j) * Y04.ell s ^ j))) := by
          refine mul_le_mul_of_nonneg_left ?_ hτ0.le
          refine (abs_add_le _ _).trans ?_
          rw [abs_mul]
          exact add_le_add le_rfl ((Finset.abs_sum_le_sum_abs _ _).trans (Finset.sum_le_sum hterm))
      _ ≤ τ * (2 * u ^ N * M + ∑ j ∈ Finset.Icc 1 (N - 1), τ ^ N * (Cf (N - j) * (s ^ (N - j) * Y04.ell s ^ j))) := by
          have := mul_le_mul hA' hGv (abs_nonneg _) (by positivity : (0:ℝ) ≤ 2 * u ^ N)
          gcongr
      _ = τ ^ (N + 1) * B s := by
          rw [hB, hu, ← Finset.mul_sum]; ring
  have hle := norm_integral_le_of_norm_le (hBint.const_mul (τ ^ (N + 1)))
    ((ae_restrict_iff' measurableSet_Ioi).2 (Filter.Eventually.of_forall key))
  rw [Real.norm_eq_abs, integral_const_mul] at hle
  calc |∫ s in Ioi (0 : ℝ), (τ * (((1 + Real.exp (-s)) ^ τ - 1) * Gf (τ * s)) - P s)|
      ≤ τ ^ (N + 1) * ∫ s in Ioi (0 : ℝ), B s := hle
    _ = (∫ s in Ioi (0 : ℝ), B s) * τ ^ (N + 1) := by ring

end R05

theorem R05.proof_1 : R05.stmt_1 := fun Gf hG hb N hN => R05.watsonN Gf hG hb N hN

namespace R05

lemma coef_expansion (k N : ℕ) (hN : 2 ≤ N) : ∃ C τ₀ : ℝ, 0 < τ₀ ∧ ∀ p : ℝ, 1 < p → Y04.tau p ≤ τ₀ →
    |Y04.a k (Y04.h p) - Y04.a k Y04.hDiamond - ∑ n ∈ Finset.Icc 2 N, c n k * Y04.tau p ^ n|
      ≤ C * Y04.tau p ^ (N + 1) := by
  obtain ⟨CW, hW⟩ := watsonN (Y04.G k) (Y04.G_contDiff k) ⟨1, fun x _ => Y04.G_bdd k x⟩ N hN
  refine ⟨8 / π * CW, 1 / 2, by norm_num, fun p hp _ => ?_⟩
  have hpi := Real.pi_pos
  obtain ⟨_, hb⟩ := hW (Y04.tau p) (Y04.tau_pos hp) (Y04.tau_lt_one hp)
  have h2 := Y04.layer_formula k hp
  set τ := Y04.tau p
  set X := τ * ∫ s in Set.Ioi (0 : ℝ), ((1 + Real.exp (-s)) ^ τ - 1) * Y04.G k (τ * s) with hX
  set S := ∑ n ∈ Finset.Icc 2 N, τ ^ n * ∑ j ∈ Finset.Icc 1 (n - 1),
            iteratedDeriv (n - 1 - j) (Y04.G k) 0 / ((n - 1 - j).factorial : ℝ) /
              (j.factorial : ℝ) * Y04.I (n - 1 - j) j with hS
  have hc : ∑ n ∈ Finset.Icc 2 N, c n k * τ ^ n = 8 / π * S := by
    rw [hS, Finset.mul_sum]
    refine Finset.sum_congr rfl (fun n _ => ?_)
    simp only [c, g, Finset.mul_sum, Finset.sum_mul]
    refine Finset.sum_congr rfl (fun j _ => ?_)
    ring
  have key : Y04.a k (Y04.h p) - Y04.a k Y04.hDiamond - ∑ n ∈ Finset.Icc 2 N, c n k * τ ^ n
      = 8 / π * (X - S) := by
    rw [hc, h2, hX]; ring
  rw [key, abs_mul, abs_of_pos (by positivity)]
  calc 8 / π * |X - S| ≤ 8 / π * (CW * τ ^ (N + 1)) := mul_le_mul_of_nonneg_left hb (by positivity)
    _ = _ := by ring

end R05

theorem R05.proof_3 : R05.stmt_3 := fun k N hN => R05.coef_expansion k N hN

namespace R05

/-- `√2 cos θ_x`. -/
noncomputable def cA (x : ℝ) : ℝ := Real.sqrt 2 * Real.cos (Y04.thetaX x)
/-- `√2 sin θ_x`. -/
noncomputable def sA (x : ℝ) : ℝ := Real.sqrt 2 * Real.sin (Y04.thetaX x)
/-- `W = 2√2 H = cA² sA`. -/
noncomputable def W (x : ℝ) : ℝ := cA x ^ 2 * sA x

lemma sqrt2_sq : Real.sqrt 2 ^ 2 = 2 := Real.sq_sqrt (by norm_num)

lemma thetaX_contDiff : ContDiff ℝ (⊤ : ℕ∞) Y04.thetaX :=
  Real.contDiff_arctan.comp (Real.contDiff_exp.comp contDiff_neg)

lemma gd_contDiff : ContDiff ℝ (⊤ : ℕ∞) gd := by
  have := thetaX_contDiff
  unfold gd; fun_prop

lemma cA_contDiff : ContDiff ℝ (⊤ : ℕ∞) cA := by
  have := thetaX_contDiff
  unfold cA; exact contDiff_const.mul (Real.contDiff_cos.comp this)

lemma sA_contDiff : ContDiff ℝ (⊤ : ℕ∞) sA := by
  have := thetaX_contDiff
  unfold sA; exact contDiff_const.mul (Real.contDiff_sin.comp this)

lemma W_contDiff : ContDiff ℝ (⊤ : ℕ∞) W := by
  unfold W; exact (cA_contDiff.pow 2).mul sA_contDiff

lemma W_eq (x : ℝ) : W x = 2 * Real.sqrt 2 * H x := by
  rw [← H_eq_trig]; simp only [W, cA, sA]
  linear_combination (2 * Real.cos (Y04.thetaX x) ^ 2 * Real.sin (Y04.thetaX x) * Real.sqrt 2) * sqrt2_sq
    - (Real.cos (Y04.thetaX x) ^ 2 * Real.sin (Y04.thetaX x) * Real.sqrt 2) * sqrt2_sq

lemma H_contDiff : ContDiff ℝ (⊤ : ℕ∞) H := by
  have : H = fun x => (2 * Real.sqrt 2)⁻¹ * W x := by
    funext x; rw [W_eq]; field_simp
  rw [this]; exact contDiff_const.mul W_contDiff

lemma cA_hasDeriv (x : ℝ) : HasDerivAt cA (sA x ^ 2 * cA x / 2) x := by
  have h := ((Real.hasDerivAt_cos (Y04.thetaX x)).comp x (Y04.theta_hasDeriv x)).const_mul (Real.sqrt 2)
  refine h.congr_deriv ?_
  simp only [sA, cA]
  linear_combination (-(Real.sqrt 2 * Real.sin (Y04.thetaX x) ^ 2 * Real.cos (Y04.thetaX x)) / 2) * sqrt2_sq

lemma sA_hasDeriv (x : ℝ) : HasDerivAt sA (-(sA x * cA x ^ 2) / 2) x := by
  have h := ((Real.hasDerivAt_sin (Y04.thetaX x)).comp x (Y04.theta_hasDeriv x)).const_mul (Real.sqrt 2)
  refine h.congr_deriv ?_
  simp only [sA, cA]
  linear_combination ((Real.sqrt 2 * Real.sin (Y04.thetaX x) * Real.cos (Y04.thetaX x) ^ 2) / 2) * sqrt2_sq

lemma gd_hasDeriv' (x : ℝ) : HasDerivAt gd (sA x * cA x) x := by
  have h := ((Y04.theta_hasDeriv x).const_mul 2).const_sub (π / 2)
  have e : gd = fun x => π / 2 - 2 * Y04.thetaX x := rfl
  rw [e]
  refine h.congr_deriv ?_
  simp only [sA, cA]
  linear_combination (-(Real.sin (Y04.thetaX x) * Real.cos (Y04.thetaX x))) * sqrt2_sq

lemma deriv_cA : deriv cA = fun x => sA x ^ 2 * cA x / 2 := funext fun x => (cA_hasDeriv x).deriv
lemma deriv_sA : deriv sA = fun x => -(sA x * cA x ^ 2) / 2 := funext fun x => (sA_hasDeriv x).deriv
lemma deriv_gd : deriv gd = fun x => sA x * cA x := funext fun x => (gd_hasDeriv' x).deriv

lemma cA_zero : cA 0 = 1 := by
  simp only [cA, Y04.thetaX_zero, Real.cos_pi_div_four]
  linear_combination sqrt2_sq / 2

lemma sA_zero : sA 0 = 1 := by
  simp only [sA, Y04.thetaX_zero, Real.sin_pi_div_four]
  linear_combination sqrt2_sq / 2

/-- The `ℚ`-algebra of functions generated by `gd`, `√2 cos θ`, `√2 sin θ`. -/
noncomputable def Alg : Subalgebra ℚ (ℝ → ℝ) := Algebra.adjoin ℚ {gd, cA, sA}

lemma gd_mem : gd ∈ Alg := Algebra.subset_adjoin (by simp)
lemma cA_mem : cA ∈ Alg := Algebra.subset_adjoin (by simp)
lemma sA_mem : sA ∈ Alg := Algebra.subset_adjoin (by simp)

lemma half_mem (f : ℝ → ℝ) (hf : f ∈ Alg) : (fun x => f x / 2) ∈ Alg := by
  have h := Alg.mul_mem (Alg.algebraMap_mem (1 / 2 : ℚ)) hf
  convert h using 1
  funext x
  simp only [Pi.mul_apply, Pi.algebraMap_apply]
  rw [show (algebraMap ℚ ℝ) (1 / 2 : ℚ) = (1 / 2 : ℝ) by simp]
  ring

lemma alg_deriv : ∀ f ∈ Alg, Differentiable ℝ f ∧ deriv f ∈ Alg := by
  intro f hf
  induction hf using Algebra.adjoin_induction with
  | mem x hx =>
    simp only [Set.mem_insert_iff, Set.mem_singleton_iff] at hx
    rcases hx with rfl | rfl | rfl
    · refine ⟨fun x => (gd_hasDeriv' x).differentiableAt, ?_⟩
      rw [deriv_gd]; exact Alg.mul_mem sA_mem cA_mem
    · refine ⟨fun x => (cA_hasDeriv x).differentiableAt, ?_⟩
      rw [deriv_cA]
      exact half_mem (fun x => sA x ^ 2 * cA x) (Alg.mul_mem (Alg.pow_mem sA_mem 2) cA_mem)
    · refine ⟨fun x => (sA_hasDeriv x).differentiableAt, ?_⟩
      rw [deriv_sA]
      have := half_mem (fun x => -(sA x * cA x ^ 2)) (Alg.neg_mem (Alg.mul_mem sA_mem (Alg.pow_mem cA_mem 2)))
      exact this
  | algebraMap r =>
    have e : (algebraMap ℚ (ℝ → ℝ) r) = fun _ => (algebraMap ℚ ℝ r) := by
      funext x; simp [Pi.algebraMap_apply]
    rw [e]
    refine ⟨differentiable_const _, ?_⟩
    rw [deriv_const']
    exact Alg.zero_mem
  | add x y hx hy ihx ihy =>
    refine ⟨ihx.1.add ihy.1, ?_⟩
    have : deriv (x + y) = deriv x + deriv y := by
      funext t; exact deriv_add (ihx.1 t) (ihy.1 t)
    rw [this]; exact Alg.add_mem ihx.2 ihy.2
  | mul x y hx hy ihx ihy =>
    refine ⟨ihx.1.mul ihy.1, ?_⟩
    have : deriv (x * y) = deriv x * y + x * deriv y := by
      funext t; exact deriv_mul (ihx.1 t) (ihy.1 t)
    rw [this]; exact Alg.add_mem (Alg.mul_mem ihx.2 hy) (Alg.mul_mem hx ihy.2)

lemma alg_iter (f : ℝ → ℝ) (hf : f ∈ Alg) (m : ℕ) : iteratedDeriv m f ∈ Alg := by
  induction m with
  | zero => simpa using hf
  | succ m ih => rw [iteratedDeriv_succ]; exact (alg_deriv _ ih).2

lemma alg_rat : ∀ f ∈ Alg, ∃ q : ℚ, f 0 = q := by
  intro f hf
  induction hf using Algebra.adjoin_induction with
  | mem x hx =>
    simp only [Set.mem_insert_iff, Set.mem_singleton_iff] at hx
    rcases hx with rfl | rfl | rfl
    · exact ⟨0, by simp [gd_zero]⟩
    · exact ⟨1, by simp [cA_zero]⟩
    · exact ⟨1, by simp [sA_zero]⟩
  | algebraMap r => exact ⟨r, by simp [Pi.algebraMap_apply]⟩
  | add x y hx hy ihx ihy =>
    obtain ⟨a, ha⟩ := ihx; obtain ⟨b, hb⟩ := ihy
    exact ⟨a + b, by simp [ha, hb]⟩
  | mul x y hx hy ihx ihy =>
    obtain ⟨a, ha⟩ := ihx; obtain ⟨b, hb⟩ := ihy
    exact ⟨a * b, by simp [ha, hb]⟩

/-- `e_{r,m}` is the Taylor coefficient of `gd^{2r} W`. -/
lemma e_eq (r m : ℕ) : e r m = tc (fun y => gd y ^ (2 * r) * W y) m := by
  have hf : (fun y => gd y ^ (2 * r) * W y) = fun y => (2 * Real.sqrt 2) * (gd y ^ (2 * r) * H y) := by
    funext y; rw [W_eq]; ring
  rw [e, tc, hf, iteratedDeriv_const_mul_field]; ring

lemma e_rat (r m : ℕ) : ∃ q : ℚ, e r m = (q : ℝ) := by
  have hmem : (fun y => gd y ^ (2 * r) * W y) ∈ Alg := by
    have : (fun y => gd y ^ (2 * r) * W y) = gd ^ (2 * r) * (cA ^ 2 * sA) := by
      funext y; simp [W]
    rw [this]
    exact Alg.mul_mem (Alg.pow_mem gd_mem _) (Alg.mul_mem (Alg.pow_mem cA_mem 2) sA_mem)
  obtain ⟨q, hq⟩ := alg_rat _ (alg_iter _ hmem m)
  refine ⟨q / (m.factorial : ℚ), ?_⟩
  rw [e_eq, tc, hq]; push_cast; ring

end R05

namespace R05

lemma tc_gd0 : tc gd 0 = 0 := by simp [tc, gd_zero]
lemma tc_gd1 : tc gd 1 = 1 := by simp [tc, deriv_gd, sA_zero, cA_zero]
lemma tc_gd2 : tc gd 2 = 0 := by
  have h : HasDerivAt (fun x => sA x * cA x) (-(sA 0 * cA 0 ^ 2) / 2 * cA 0 + sA 0 * (sA 0 ^ 2 * cA 0 / 2)) 0 :=
    (sA_hasDeriv 0).mul (cA_hasDeriv 0)
  have : iteratedDeriv 2 gd 0 = 0 := by
    rw [iteratedDeriv_succ, iteratedDeriv_one, deriv_gd, h.deriv, sA_zero, cA_zero]; ring
  simp [tc, this]
lemma tc_W0 : tc W 0 = 1 := by simp [tc, W, sA_zero, cA_zero]
lemma tc_W1 : tc W 1 = 1 / 2 := by
  have h : HasDerivAt W (2 * cA 0 * (sA 0 ^ 2 * cA 0 / 2) * sA 0 + cA 0 ^ 2 * (-(sA 0 * cA 0 ^ 2) / 2)) 0 := by
    have := ((cA_hasDeriv 0).pow 2).mul (sA_hasDeriv 0)
    convert this using 1
    · funext x; simp [W]
    · simp only [Pi.pow_apply]; norm_num
  simp only [tc, iteratedDeriv_one, h.deriv, sA_zero, cA_zero]; norm_num

/-- `ts gd = X · U`. -/
noncomputable def U : PowerSeries ℝ := PowerSeries.mk (fun i => tc gd (i + 1))

lemma ts_gd : ts gd = PowerSeries.X * U := by
  ext n
  cases n with
  | zero => simp [coeff_ts, tc_gd0]
  | succ n => rw [PowerSeries.coeff_succ_X_mul]; simp [U, coeff_ts]

lemma U0 : PowerSeries.constantCoeff U = 1 := by
  rw [← PowerSeries.coeff_zero_eq_constantCoeff_apply]; simp [U, tc_gd1]

lemma U1 : PowerSeries.coeff 1 U = 0 := by simp [U, tc_gd2]

lemma Upow1 (n : ℕ) : PowerSeries.coeff 1 (U ^ n) = 0 := by
  induction n with
  | zero => simp
  | succ n ih =>
    rw [pow_succ, PowerSeries.coeff_mul, Finset.Nat.antidiagonal_succ]
    simp [ih, U1]

lemma ts_F (r : ℕ) : ts (fun y => gd y ^ (2 * r) * W y) =
    PowerSeries.X ^ (2 * r) * (U ^ (2 * r) * ts W) := by
  rw [ts_mul (gd_contDiff.pow _) W_contDiff, ts_pow gd_contDiff, ts_gd]; ring

lemma e_coeff (r m : ℕ) : e r m = PowerSeries.coeff m
    (PowerSeries.X ^ (2 * r) * (U ^ (2 * r) * ts W)) := by
  rw [e_eq, ← ts_F, coeff_ts]

lemma e_low (r m : ℕ) (h : m < 2 * r) : e r m = 0 := by
  rw [e_coeff, PowerSeries.coeff_X_pow_mul']
  simp only [show ¬ 2 * r ≤ m by omega, ↓reduceIte]

lemma e_diag (r : ℕ) : e r (2 * r) = 1 := by
  rw [e_coeff, PowerSeries.coeff_X_pow_mul']
  simp only [le_refl, ↓reduceIte]
  rw [Nat.sub_self,
    PowerSeries.coeff_zero_eq_constantCoeff_apply, map_mul, map_pow, U0,
    ← PowerSeries.coeff_zero_eq_constantCoeff_apply, coeff_ts, tc_W0]
  simp

lemma e_next (r : ℕ) : e r (2 * r + 1) = 1 / 2 := by
  rw [e_coeff, PowerSeries.coeff_X_pow_mul']
  simp only [show 2 * r ≤ 2 * r + 1 by omega, ↓reduceIte]
  rw [show 2 * r + 1 - 2 * r = 1 by omega, PowerSeries.coeff_mul, Finset.Nat.antidiagonal_succ]
  simp only [Finset.Nat.antidiagonal_zero, Finset.map_singleton, Finset.sum_cons,
    Finset.sum_singleton]
  simp [Upow1, coeff_ts, tc_W1, PowerSeries.coeff_zero_eq_constantCoeff_apply, U0]

lemma e_props : (∀ r : ℕ, ContDiff ℝ (⊤ : ℕ∞) (fun y : ℝ => gd y ^ (2 * r) * H y)) ∧
    (∀ r m : ℕ, ∃ q : ℚ, e r m = (q : ℝ)) ∧
    (∀ r m : ℕ, m < 2 * r → e r m = 0) ∧
    (∀ r : ℕ, e r (2 * r) = 1) ∧
    (∀ r : ℕ, e r (2 * r + 1) = 1 / 2) :=
  ⟨fun _ => (gd_contDiff.pow _).mul H_contDiff, e_rat, e_low, e_diag, e_next⟩

end R05

theorem R05.proof_4 : R05.stmt_4 := R05.e_props

namespace R05

lemma two_pow_le_fac (n R : ℕ) : (2 : ℝ) ^ n ≤ ((2 * (n + R + 1)).factorial : ℝ) := by
  have h1 := @Nat.factorial_mul_pow_le_factorial 1 n
  have h2 : (1 + n).factorial ≤ (2 * (n + R + 1)).factorial := Nat.factorial_le (by omega)
  have h1' : 2 ^ n ≤ (1 + n).factorial := by simpa using h1
  have : 2 ^ n ≤ (2 * (n + R + 1)).factorial := h1'.trans h2
  exact_mod_cast this

lemma cos_tail {z : ℝ} (hz : |z| ≤ 1) (R : ℕ) :
    |Real.cos z - ∑ n ∈ Finset.range (R + 1), (-1 : ℝ) ^ n * z ^ (2 * n) / ((2 * n).factorial : ℝ)|
      ≤ |z| ^ (2 * R + 2) * 2 := by
  have hs := (hasSum_nat_add_iff' (R + 1)).mpr (Real.hasSum_cos z)
  have hg := hasSum_geometric_two.mul_left (|z| ^ (2 * R + 2))
  have := hs.norm_le_of_bounded hg (fun n => ?_)
  · simpa [Real.norm_eq_abs] using this
  · rw [Real.norm_eq_abs, abs_div, abs_mul, abs_pow, abs_neg, abs_one, one_pow, one_mul, abs_pow,
      Nat.abs_cast]
    have hf : (0 : ℝ) < ((2 * (n + (R + 1))).factorial : ℝ) := by positivity
    rw [div_le_iff₀ hf]
    have h2 := two_pow_le_fac n R
    rw [show n + (R + 1) = n + R + 1 by ring] at hf ⊢
    have hz0 := abs_nonneg z
    have hzn : |z| ^ (2 * n) ≤ 1 := pow_le_one₀ hz0 hz
    have e1 : |z| ^ (2 * (n + R + 1)) = |z| ^ (2 * R + 2) * |z| ^ (2 * n) := by ring
    rw [e1]
    have hzR : 0 ≤ |z| ^ (2 * R + 2) := pow_nonneg hz0 _
    have : (1 / 2 : ℝ) ^ n * ((2 * (n + R + 1)).factorial : ℝ) ≥ 1 := by
      rw [one_div, inv_pow, ge_iff_le, le_inv_mul_iff₀ (by positivity)]; simpa using h2
    calc |z| ^ (2 * R + 2) * |z| ^ (2 * n) ≤ |z| ^ (2 * R + 2) * 1 := mul_le_mul_of_nonneg_left hzn hzR
      _ ≤ |z| ^ (2 * R + 2) * ((1 / 2 : ℝ) ^ n * ((2 * (n + R + 1)).factorial : ℝ)) :=
          mul_le_mul_of_nonneg_left this hzR
      _ = _ := by ring

lemma gd_abs_le (y : ℝ) : |gd y| ≤ |y| := by
  have hb : ∀ z ∈ Set.uIcc 0 y, ‖deriv gd z‖ ≤ 1 := by
    intro z _
    rw [(gd_hasDeriv z).deriv, Real.norm_eq_abs, abs_div, abs_one,
      abs_of_pos (Real.cosh_pos z)]
    exact div_le_one_of_le₀ (Real.one_le_cosh z) (Real.cosh_pos z).le
  have := (convex_uIcc (0 : ℝ) y).norm_image_sub_le_of_norm_deriv_le
    (fun z _ => (gd_hasDeriv z).differentiableAt) hb Set.left_mem_uIcc Set.right_mem_uIcc
  simpa [gd_zero, Real.norm_eq_abs] using this

/-- The truncated cosine series in `gd`. -/
noncomputable def Pc (k R : ℕ) (x : ℝ) : ℝ :=
  ∑ n ∈ Finset.range (R + 1), (-1 : ℝ) ^ n * (2 * k * gd x) ^ (2 * n) / ((2 * n).factorial : ℝ)

lemma Pc_contDiff (k R : ℕ) : ContDiff ℝ (⊤ : ℕ∞) (Pc k R) := by
  unfold Pc; have := gd_contDiff; fun_prop

lemma cosgd_contDiff (k : ℕ) : ContDiff ℝ (⊤ : ℕ∞) (fun x => Real.cos (2 * k * gd x)) := by
  have := gd_contDiff; exact Real.contDiff_cos.comp (contDiff_const.mul this)

lemma E_flat (k R : ℕ) : ∀ i < 2 * R + 2,
    iteratedDeriv i (fun x => Real.cos (2 * k * gd x) - Pc k R x) 0 = 0 := by
  refine flat_unique (2 * R + 2) _ ((cosgd_contDiff k).sub (Pc_contDiff k R))
    (δ := 1 / (2 * k + 1)) (C := 2 * (2 * k) ^ (2 * R + 2)) (by positivity) (fun y hy => ?_)
  have hk : (0 : ℝ) ≤ k := Nat.cast_nonneg k
  have hz : |2 * k * gd y| ≤ 2 * k * |y| := by
    rw [abs_mul, abs_of_nonneg (by positivity : (0:ℝ) ≤ 2 * k)]
    exact mul_le_mul_of_nonneg_left (gd_abs_le y) (by positivity)
  have hz1 : |2 * k * gd y| ≤ 1 := by
    refine hz.trans ?_
    have : 2 * (k : ℝ) * |y| ≤ 2 * k * (1 / (2 * k + 1)) := mul_le_mul_of_nonneg_left hy (by positivity)
    refine this.trans ?_
    rw [mul_one_div, div_le_one (by positivity)]; linarith
  have := cos_tail hz1 R
  simp only [Pc]
  calc _ ≤ |2 * k * gd y| ^ (2 * R + 2) * 2 := this
    _ ≤ (2 * k * |y|) ^ (2 * R + 2) * 2 := by gcongr
    _ = 2 * (2 * k) ^ (2 * R + 2) * |y| ^ (2 * R + 2) := by ring

lemma tc_gdH (r m : ℕ) : tc (fun y => gd y ^ (2 * r) * H y) m = e r m / (2 * Real.sqrt 2) := by
  rw [e, tc]; field_simp

lemma g_formula (k m : ℕ) :
    g k m = (-1 : ℝ) ^ k * Real.sqrt 2 / 4 * (rho m).eval ((k : ℝ) ^ 2) := by
  have hG : Y04.G k = fun x => (-1 : ℝ) ^ k * (Real.cos (2 * k * gd x) * H x) := by
    funext x; rw [G_eq_gd]; ring
  set E : ℝ → ℝ := fun x => Real.cos (2 * k * gd x) - Pc k m x with hE
  have hEs : ContDiff ℝ (⊤ : ℕ∞) E := (cosgd_contDiff k).sub (Pc_contDiff k m)
  have hcos : (fun x => Real.cos (2 * k * gd x)) = fun x => Pc k m x + E x := by
    funext x; simp [hE]
  have hts : ts (Y04.G k) = PowerSeries.C ((-1 : ℝ) ^ k) *
      (ts (Pc k m) * ts H + ts E * ts H) := by
    rw [hG, ts_cmul, ts_mul (cosgd_contDiff k) H_contDiff, hcos, ts_add (Pc_contDiff k m) hEs]
    ring
  have hPc : ts (Pc k m) * ts H = ∑ n ∈ Finset.range (m + 1),
      PowerSeries.C ((-1 : ℝ) ^ n * (2 * k) ^ (2 * n) / ((2 * n).factorial : ℝ)) *
        ts (fun y => gd y ^ (2 * n) * H y) := by
    have e1 : Pc k m = fun x => ∑ n ∈ Finset.range (m + 1),
        (fun n x => (-1 : ℝ) ^ n * (2 * k) ^ (2 * n) / ((2 * n).factorial : ℝ) * gd x ^ (2 * n)) n x := by
      funext x; simp only [Pc]
      refine Finset.sum_congr rfl (fun n _ => ?_); ring
    rw [e1, ts_sum _ _ (fun n _ => by have := gd_contDiff; fun_prop), Finset.sum_mul]
    refine Finset.sum_congr rfl (fun n _ => ?_)
    rw [ts_cmul, ts_mul (gd_contDiff.pow _) H_contDiff, ts_pow gd_contDiff]; ring
  have hEz : PowerSeries.coeff m (ts E * ts H) = 0 := by
    rw [PowerSeries.coeff_mul]
    refine Finset.sum_eq_zero (fun p hp => ?_)
    have hp1 : p.1 ≤ m := by have h := hp; simp at h; omega
    rw [coeff_ts, tc, E_flat k m p.1 (by omega)]; simp
  have hg : g k m = PowerSeries.coeff m (ts (Y04.G k)) := by rw [coeff_ts]; rfl
  rw [hg, hts, PowerSeries.coeff_C_mul, map_add, hEz, add_zero, hPc, map_sum]
  simp only [PowerSeries.coeff_C_mul, coeff_ts, tc_gdH]
  -- the polynomial side
  have hrho : (rho m).eval ((k : ℝ) ^ 2) = ∑ r ∈ Finset.range (m + 1),
      (-4 * (k : ℝ) ^ 2) ^ r / ((2 * r).factorial : ℝ) * e r m := by
    calc (rho m).eval ((k : ℝ) ^ 2) = ∑ r ∈ Finset.range (m / 2 + 1),
          (-4 * (k : ℝ) ^ 2) ^ r / ((2 * r).factorial : ℝ) * e r m := by
          rw [rho, Polynomial.eval_finsetSum]
          refine Finset.sum_congr rfl (fun r _ => ?_)
          simp only [Polynomial.eval_mul, Polynomial.eval_C, Polynomial.eval_pow, Polynomial.eval_X]
          ring
      _ = _ := by
          refine Finset.sum_subset (Finset.range_subset_range.mpr (by omega)) (fun r hr hr' => ?_)
          have h1 : ¬ r < m / 2 + 1 := fun h => hr' (Finset.mem_range.mpr h)
          rw [e_low r m (by omega)]; ring
  rw [hrho, Finset.mul_sum, Finset.mul_sum]
  refine Finset.sum_congr rfl (fun r _ => ?_)
  have hs : Real.sqrt 2 ^ 2 = 2 := sqrt2_sq
  have hs0 : Real.sqrt 2 ≠ 0 := by positivity
  rw [show ((-4 : ℝ) * (k : ℝ) ^ 2) ^ r = (-1) ^ r * 4 ^ r * ((k : ℝ) ^ 2) ^ r by
    rw [← mul_pow, ← mul_pow]; norm_num]
  rw [mul_pow, pow_mul, pow_mul]
  field_simp
  rw [hs]
  ring

end R05

theorem R05.proof_5 : R05.stmt_5 := R05.g_formula

namespace R05

lemma ite_f {c : Prop} [Decidable c] {a b : ℝ} (h : ¬c) : (if c then a else b) = b := by simp [h]
lemma ite_t {c : Prop} [Decidable c] {a b : ℝ} (h : c) : (if c then a else b) = a := by simp [h]

lemma neg_one_sq (k : ℕ) : (-1 : ℝ) ^ k * (-1) ^ k = 1 := by
  rw [← mul_pow]; norm_num

lemma theta_eval (n k : ℕ) :
    (-1 : ℝ) ^ k * (π / Real.sqrt 2) * c n k = (Theta n).eval ((k : ℝ) ^ 2) := by
  have hpi := Real.pi_pos
  have hs0 : Real.sqrt 2 ≠ 0 := by positivity
  simp only [Theta, c, g_formula, Polynomial.eval_mul, Polynomial.eval_finsetSum, Polynomial.eval_C,
    Polynomial.eval_ofNat]
  rw [Finset.mul_sum, Finset.mul_sum, Finset.mul_sum]
  refine Finset.sum_congr rfl (fun j _ => ?_)
  have hk := neg_one_sq k
  have hj : (j.factorial : ℝ) ≠ 0 := by positivity
  field_simp
  linear_combination (8 * (Polynomial.eval ((k:ℝ) ^ 2) (rho (n - 1 - j))) * Y04.I (n - 1 - j) j) * hk

lemma rho_coeff (m i : ℕ) : (rho m).coeff i =
    if i < m / 2 + 1 then (-4 : ℝ) ^ i / ((2 * i).factorial : ℝ) * e i m else 0 := by
  rw [rho, Polynomial.finsetSum_coeff]
  simp only [Polynomial.coeff_C_mul, Polynomial.coeff_X_pow, mul_ite, mul_one, mul_zero]
  rw [Finset.sum_ite_eq]
  simp [Finset.mem_range]

lemma theta_coeff (n i : ℕ) : (Theta n).coeff i = ∑ j ∈ Finset.Icc 1 (n - 1),
    2 * (Y04.I (n - 1 - j) j / (j.factorial : ℝ)) * (rho (n - 1 - j)).coeff i := by
  rw [Theta, ← Polynomial.C_ofNat, Polynomial.coeff_C_mul,
    Polynomial.finsetSum_coeff, Finset.mul_sum]
  refine Finset.sum_congr rfl (fun j _ => ?_)
  rw [Polynomial.coeff_C_mul]; ring

/-- A rational witness for `e_{r,m}`. -/
noncomputable def eQ (r m : ℕ) : ℚ := (e_rat r m).choose

lemma eQ_spec (r m : ℕ) : e r m = (eQ r m : ℝ) := (e_rat r m).choose_spec

lemma theta_rat (n i : ℕ) : ∃ q : ℕ → ℚ,
    (Theta n).coeff i = ∑ j ∈ Finset.Icc 1 (n - 1), (q j : ℝ) * Y04.I (n - 1 - j) j := by
  refine ⟨fun j => 2 / (j.factorial : ℚ) *
    (if i < (n - 1 - j) / 2 + 1 then (-4 : ℚ) ^ i / ((2 * i).factorial : ℚ) * eQ i (n - 1 - j) else 0), ?_⟩
  rw [theta_coeff]
  refine Finset.sum_congr rfl (fun j _ => ?_)
  rw [rho_coeff]
  by_cases h : i < (n - 1 - j) / 2 + 1
  · simp only [h, ↓reduceIte]; rw [eQ_spec]; push_cast; ring
  · simp only [h, ↓reduceIte]; push_cast; ring

lemma ell_pos (s : ℝ) : 0 < Y04.ell s := Real.log_pos (by linarith [Real.exp_pos (-s)])

lemma I_pos (m j : ℕ) (hj : 1 ≤ j) : 0 < Y04.I m j := by
  unfold Y04.I
  rw [setIntegral_pos_iff_support_of_nonneg_ae _ (Y04.I_int m j hj)]
  · have hsub : Ioi (0 : ℝ) ⊆ Function.support (fun s : ℝ => s ^ m * Y04.ell s ^ j) ∩ Ioi 0 := by
      intro s hs
      refine ⟨?_, hs⟩
      rw [Function.mem_support]
      exact (mul_pos (pow_pos hs m) (pow_pos (ell_pos s) j)).ne'
    refine lt_of_lt_of_le ?_ (measure_mono hsub)
    rw [Real.volume_Ioi]; exact ENNReal.zero_lt_top
  · exact (ae_restrict_iff' measurableSet_Ioi).2 (Filter.Eventually.of_forall
      (fun s (hs : 0 < s) => (mul_pos (pow_pos hs m) (pow_pos (ell_pos s) j)).le))

lemma theta_coeff_high (n i : ℕ) (hn : 2 ≤ n) (hi : n / 2 - 1 < i) : (Theta n).coeff i = 0 := by
  rw [theta_coeff]
  refine Finset.sum_eq_zero (fun j hj => ?_)
  simp only [Finset.mem_Icc] at hj
  rw [rho_coeff, ite_f (by omega), mul_zero]

lemma theta_even_top (m : ℕ) (hm : 1 ≤ m) : (Theta (2 * m)).coeff (2 * m / 2 - 1) =
    (-1 : ℝ) ^ (m - 1) * (2 ^ (2 * m - 1) - 1) * (riemannZeta ((2 * m : ℕ) : ℂ)).re := by
  rw [theta_coeff, Finset.sum_eq_single 1]
  · rw [rho_coeff, ite_t (by omega)]
    have hd : 2 * m / 2 - 1 = m - 1 := by omega
    have h21 : 2 * m - 1 - 1 = 2 * (m - 1) := by omega
    rw [hd, h21, e_diag]
    have hI := (eta_m (2 * (m - 1))).2
    rw [show 2 * (m - 1) + 2 = 2 * m by omega] at hI
    rw [hI]
    have hf : ((2 * (m - 1)).factorial : ℝ) ≠ 0 := by positivity
    obtain ⟨m', rfl⟩ : ∃ m', m = m' + 1 := ⟨m - 1, by omega⟩
    simp only [Nat.add_sub_cancel, show 2 * (m' + 1) - 1 = 2 * m' + 1 by omega]
    field_simp
    rw [show ((-4 : ℝ)) ^ m' = (-1) ^ m' * 4 ^ m' by rw [← mul_pow]; norm_num,
      show (4 : ℝ) ^ m' = 2 ^ (2 * m') by rw [pow_mul]; norm_num]
    have h2 : (2 : ℝ)⁻¹ ^ (2 * m') * 2 ^ (2 * m') = 1 := by rw [← mul_pow]; norm_num
    linear_combination (-((riemannZeta ((2 * (m' + 1) : ℕ) : ℂ)).re * (-1 : ℝ) ^ m')) * h2
  · intro j hj hj1
    simp only [Finset.mem_Icc] at hj
    rw [rho_coeff, ite_f (by omega), mul_zero]
  · intro h; exact absurd (Finset.mem_Icc.mpr ⟨le_rfl, by omega⟩) h

lemma theta_odd_top (M : ℕ) (hM : 1 ≤ M) : (Theta (2 * M + 1)).coeff ((2 * M + 1) / 2 - 1) =
    (-4 : ℝ) ^ (M - 1) / ((2 * M - 2).factorial : ℝ) *
      (Y04.I (2 * M - 1) 1 + Y04.I (2 * M - 2) 2) := by
  rw [theta_coeff, ← Finset.sum_subset (s₁ := {1, 2})]
  · rw [Finset.sum_pair (by norm_num)]
    have hd : (2 * M + 1) / 2 - 1 = M - 1 := by omega
    rw [hd, rho_coeff, rho_coeff, ite_t (by omega), ite_t (by omega)]
    rw [show 2 * M + 1 - 1 - 1 = 2 * (M - 1) + 1 by omega,
      show 2 * M + 1 - 1 - 2 = 2 * (M - 1) by omega, e_next, e_diag,
      show 2 * (M - 1) + 1 = 2 * M - 1 by omega, show 2 * (M - 1) = 2 * M - 2 by omega]
    simp only [Nat.factorial_one, Nat.factorial_two, Nat.cast_one, Nat.cast_ofNat]
    ring
  · intro j hj
    simp only [Finset.mem_insert, Finset.mem_singleton] at hj
    simp only [Finset.mem_Icc]; omega
  · intro j hj hj'
    simp only [Finset.mem_Icc] at hj
    simp only [Finset.mem_insert, Finset.mem_singleton] at hj'
    rw [rho_coeff, ite_f (by omega), mul_zero]

lemma natDeg_of (p : Polynomial ℝ) (d : ℕ) (hhi : ∀ i, d < i → p.coeff i = 0) (htop : p.coeff d ≠ 0) :
    p ≠ 0 ∧ p.natDegree = d ∧ p.leadingCoeff = p.coeff d := by
  have hne : p ≠ 0 := fun h => htop (by rw [h, Polynomial.coeff_zero])
  have hle : p.natDegree ≤ d := by
    rw [Polynomial.natDegree_le_iff_coeff_eq_zero]
    intro i hi; exact hhi i (by exact_mod_cast hi)
  have hdeg : p.natDegree = d := le_antisymm hle (Polynomial.le_natDegree_of_ne_zero htop)
  exact ⟨hne, hdeg, by rw [Polynomial.leadingCoeff, hdeg]⟩

lemma zeta_pos (n : ℕ) (hn : 1 < n) : 0 < (riemannZeta ((n : ℕ) : ℂ)).re := Q03.zn_pos hn

lemma theta_degree (n : ℕ) (hn : 2 ≤ n) :
    Theta n ≠ 0 ∧ (Theta n).natDegree = n / 2 - 1 ∧
      (Theta n).leadingCoeff = (Theta n).coeff (n / 2 - 1) := by
  refine natDeg_of _ _ (fun i hi => theta_coeff_high n i hn hi) ?_
  obtain ⟨m, rfl | rfl⟩ := Nat.even_or_odd' n
  · rw [theta_even_top m (by omega)]
    have h1 : (0 : ℝ) < 2 ^ (2 * m - 1) - 1 := by
      have : (2 : ℝ) ≤ 2 ^ (2 * m - 1) := by
        calc (2 : ℝ) = 2 ^ 1 := by norm_num
          _ ≤ 2 ^ (2 * m - 1) := pow_le_pow_right₀ (by norm_num) (by omega)
      linarith
    have h2 := zeta_pos (2 * m) (by omega)
    have h3 : (-1 : ℝ) ^ (m - 1) ≠ 0 := by positivity
    exact mul_ne_zero (mul_ne_zero h3 h1.ne') h2.ne'
  · rw [theta_odd_top m (by omega)]
    have h1 := I_pos (2 * m - 1) 1 le_rfl
    have h2 := I_pos (2 * m - 2) 2 (by norm_num)
    have h3 : (-4 : ℝ) ^ (m - 1) ≠ 0 := by positivity
    have h4 : ((2 * m - 2).factorial : ℝ) ≠ 0 := by positivity
    exact mul_ne_zero (div_ne_zero h3 h4) (by linarith)

end R05

theorem R05.proof_6 : R05.stmt_6 := fun n _ => ⟨fun k => R05.theta_eval n k, fun i => R05.theta_rat n i⟩

namespace R05

lemma expansion_theta (k N : ℕ) (hN : 2 ≤ N) : ∃ C τ₀ : ℝ, 0 < τ₀ ∧ ∀ p : ℝ, 1 < p → Y04.tau p ≤ τ₀ →
    |Y04.a k (Y04.h p) - (-1 : ℝ) ^ k * (Real.sqrt 2 / π) *
        (-(4 / (16 * (k : ℝ) ^ 2 - 1)) +
          ∑ n ∈ Finset.Icc 2 N, (Theta n).eval ((k : ℝ) ^ 2) * Y04.tau p ^ n)|
      ≤ C * Y04.tau p ^ (N + 1) := by
  obtain ⟨C, τ₀, hτ₀, h⟩ := coef_expansion k N hN
  refine ⟨C, τ₀, hτ₀, fun p hp hτ => ?_⟩
  have hpi := Real.pi_pos
  have hs0 : Real.sqrt 2 ≠ 0 := by positivity
  have hss : Real.sqrt 2 ^ 2 = 2 := sqrt2_sq
  have a1 := Y04.four_k_ne k; have b1 := Y04.four_k_ne' k
  have c1 : (16 * (k : ℝ) ^ 2 - 1) ≠ 0 := by
    rw [show (16 * (k : ℝ) ^ 2 - 1) = (4 * k - 1) * (4 * k + 1) by ring]; exact mul_ne_zero a1 b1
  have hk := neg_one_sq k
  have e1 : (-1 : ℝ) ^ k * (Real.sqrt 2 / π) *
        (-(4 / (16 * (k : ℝ) ^ 2 - 1)) +
          ∑ n ∈ Finset.Icc 2 N, (Theta n).eval ((k : ℝ) ^ 2) * Y04.tau p ^ n) =
      Y04.a k Y04.hDiamond + ∑ n ∈ Finset.Icc 2 N, c n k * Y04.tau p ^ n := by
    rw [Y04.A_formula, mul_add, Finset.mul_sum]
    congr 1
    · field_simp
    · refine Finset.sum_congr rfl (fun n _ => ?_)
      rw [← theta_eval]
      field_simp
      linear_combination (c n k * Y04.tau p ^ n) * hk
  rw [e1, ← sub_sub]
  exact h p hp hτ

lemma degree_all : (∀ k N : ℕ, 2 ≤ N → ∃ C τ₀ : ℝ, 0 < τ₀ ∧ ∀ p : ℝ, 1 < p → Y04.tau p ≤ τ₀ →
    |Y04.a k (Y04.h p) - (-1 : ℝ) ^ k * (Real.sqrt 2 / π) *
        (-(4 / (16 * (k : ℝ) ^ 2 - 1)) +
          ∑ n ∈ Finset.Icc 2 N, (Theta n).eval ((k : ℝ) ^ 2) * Y04.tau p ^ n)|
      ≤ C * Y04.tau p ^ (N + 1)) ∧
  (∀ n : ℕ, 2 ≤ n → Theta n ≠ 0 ∧ (Theta n).natDegree = n / 2 - 1) ∧
  (∀ m : ℕ, 1 ≤ m → (Theta (2 * m)).leadingCoeff =
      (-1 : ℝ) ^ (m - 1) * (2 ^ (2 * m - 1) - 1) * (riemannZeta ((2 * m : ℕ) : ℂ)).re) ∧
  (∀ M : ℕ, 1 ≤ M →
    0 < Y04.I (2 * M - 1) 1 ∧ 0 < Y04.I (2 * M - 2) 2 ∧
    (Theta (2 * M + 1)).leadingCoeff =
      (-4 : ℝ) ^ (M - 1) / ((2 * M - 2).factorial : ℝ) *
        (Y04.I (2 * M - 1) 1 + Y04.I (2 * M - 2) 2)) := by
  refine ⟨expansion_theta, fun n hn => ⟨(theta_degree n hn).1, (theta_degree n hn).2.1⟩,
    fun m hm => ?_, fun M hM => ⟨I_pos _ 1 le_rfl, I_pos _ 2 (by norm_num), ?_⟩⟩
  · rw [(theta_degree (2 * m) (by omega)).2.2, theta_even_top m hm]
  · rw [(theta_degree (2 * M + 1) (by omega)).2.2, theta_odd_top M hM]

end R05

theorem R05.proof_9 : R05.stmt_9 := R05.degree_all

namespace R05

lemma poly_eq_of_sq (p q : Polynomial ℝ) (h : ∀ k : ℕ, p.eval ((k : ℝ) ^ 2) = q.eval ((k : ℝ) ^ 2)) :
    p = q := by
  apply Polynomial.eq_of_infinite_eval_eq
  have hinj : Function.Injective (fun k : ℕ => (k : ℝ) ^ 2) := by
    intro a b hab
    have : ((a ^ 2 : ℕ) : ℝ) = ((b ^ 2 : ℕ) : ℝ) := by push_cast; exact hab
    exact Nat.pow_left_injective (by norm_num) (Nat.cast_injective this)
  refine (Set.infinite_range_of_injective hinj).mono ?_
  rintro x ⟨k, rfl⟩
  exact h k

lemma common (k : ℕ) : ((-1 : ℝ) ^ k) ^ 2 = 1 ∧ Real.sqrt 2 ^ 2 = 2 :=
  ⟨by rw [← pow_mul, mul_comm, pow_mul]; norm_num, sqrt2_sq⟩

lemma th2 (k : ℕ) : (Theta 2).eval ((k : ℝ) ^ 2) = π ^ 2 / 6 := by
  rw [← theta_eval]
  simp only [c, g, show Finset.Icc 1 (2 - 1) = {1} by decide, Finset.sum_singleton]
  simp only [show 2 - 1 - 1 = 0 by rfl, iteratedDeriv_zero, Y04.G_zero, Y04.I01]
  have hpi := Real.pi_pos
  have hs0 : Real.sqrt 2 ≠ 0 := by positivity
  obtain ⟨h1, h2⟩ := common k
  field_simp
  simp only [Nat.factorial_zero, Nat.factorial_one, Nat.cast_one]
  linear_combination 8 * h1

lemma th3 (k : ℕ) : (Theta 3).eval ((k : ℝ) ^ 2) = Y04.z3 := by
  rw [← theta_eval]
  simp only [c, g, show Finset.Icc 1 (3 - 1) = {1, 2} by decide]
  rw [Finset.sum_pair (by norm_num)]
  simp only [show 3 - 1 - 1 = 1 by rfl, show 3 - 1 - 2 = 0 by rfl, iteratedDeriv_zero,
    iteratedDeriv_one, Y04.G_zero, Y04.G_d1, Y04.I11, Y04.I02]
  have hpi := Real.pi_pos
  have hs0 : Real.sqrt 2 ≠ 0 := by positivity
  obtain ⟨h1, h2⟩ := common k
  simp only [Nat.factorial_zero, Nat.factorial_one, Nat.factorial_two, Nat.cast_one, Nat.cast_ofNat]
  field_simp
  linear_combination (32 * Y04.z3) * h1

lemma th4 (k : ℕ) : (Theta 4).eval ((k : ℝ) ^ 2) =
    -(33 / 16) * (π ^ 4 / 90) + -7 * (π ^ 4 / 90) * (k : ℝ) ^ 2 := by
  rw [← theta_eval]
  simp only [c, g, show Finset.Icc 1 (4 - 1) = {1, 2, 3} by decide]
  rw [Finset.sum_insert (by decide), Finset.sum_pair (by norm_num)]
  simp only [show 4 - 1 - 1 = 2 by rfl, show 4 - 1 - 2 = 1 by rfl, show 4 - 1 - 3 = 0 by rfl,
    iteratedDeriv_zero, iteratedDeriv_one, Y04.G_zero, Y04.G_d1]
  have hd2 := Y04.G_d2 k
  have hI03 : Y04.I 0 3 = (3 / 4 * (π ^ 4 / 90) - 3 * Y04.I 1 2) / 2 := by linarith [Y04.I0312]
  have hd2' : iteratedDeriv 2 (Y04.G k) 0 = 2 * -((-1 : ℝ) ^ k * Real.sqrt 2 * ((k : ℝ) ^ 2 / 2 + 5 / 32)) := by
    linarith
  rw [hd2', hI03, Y04.I21]
  have hpi := Real.pi_pos
  have hs0 : Real.sqrt 2 ≠ 0 := by positivity
  obtain ⟨h1, h2⟩ := common k
  simp only [Nat.factorial, Nat.cast_one, Nat.succ_eq_add_one]
  norm_num
  field_simp
  ring_nf
  linear_combination (-(688128 * π ^ 4 * (k : ℝ) ^ 2 + 202752 * π ^ 4)) * h1

end R05

namespace R05

lemma th5 (k : ℕ) : (Theta 5).eval ((k : ℝ) ^ 2) =
    -(3 / 4 + 4 * (k : ℝ) ^ 2) * Z02.z5 - (7 / 4 + 4 * (k : ℝ) ^ 2) * ((π ^ 2 / 6) * Y04.z3) := by
  rw [← theta_eval]; exact Z02.theta5 k

lemma th6 (k : ℕ) : (Theta 6).eval ((k : ℝ) ^ 2) =
    (4039 / 256 + 733 / 8 * (k : ℝ) ^ 2 + 31 * (k : ℝ) ^ 4) * (π ^ 6 / 945)
      - (7 / 8 + 2 * (k : ℝ) ^ 2) * Y04.z3 ^ 2 := by
  rw [← theta_eval]; exact V02.theta6 k

lemma th7 (k : ℕ) : (Theta 7).eval ((k : ℝ) ^ 2) =
    (129 / 16 + 46 * (k : ℝ) ^ 2 + 16 * (k : ℝ) ^ 4) * Q03.z7
      + (141 / 16 + 50 * (k : ℝ) ^ 2 + 16 * (k : ℝ) ^ 4) * ((π ^ 2 / 6) * Z02.z5)
      + (1071 / 160 + 181 / 5 * (k : ℝ) ^ 2 + 56 / 5 * (k : ℝ) ^ 4) * ((π ^ 2 / 6) ^ 2 * Y04.z3) := by
  rw [← theta_eval]; exact Q03.theta7 k

lemma theta_values : Theta 2 = Polynomial.C (π ^ 2 / 6) ∧
  Theta 3 = Polynomial.C Y04.z3 ∧
  Theta 4 = Polynomial.C (-(33 / 16) * (π ^ 4 / 90)) + Polynomial.C (-7 * (π ^ 4 / 90)) * Polynomial.X ∧
  Theta 5 = Polynomial.C (-(3 / 4) * Z02.z5 - 7 / 4 * ((π ^ 2 / 6) * Y04.z3))
      + Polynomial.C (-4 * Z02.z5 - 4 * ((π ^ 2 / 6) * Y04.z3)) * Polynomial.X ∧
  Theta 6 = Polynomial.C (4039 / 256 * (π ^ 6 / 945) - 7 / 8 * Y04.z3 ^ 2)
      + Polynomial.C (733 / 8 * (π ^ 6 / 945) - 2 * Y04.z3 ^ 2) * Polynomial.X
      + Polynomial.C (31 * (π ^ 6 / 945)) * Polynomial.X ^ 2 ∧
  Theta 7 = Polynomial.C (129 / 16 * Q03.z7 + 141 / 16 * ((π ^ 2 / 6) * Z02.z5)
        + 1071 / 160 * ((π ^ 2 / 6) ^ 2 * Y04.z3))
      + Polynomial.C (46 * Q03.z7 + 50 * ((π ^ 2 / 6) * Z02.z5)
        + 181 / 5 * ((π ^ 2 / 6) ^ 2 * Y04.z3)) * Polynomial.X
      + Polynomial.C (16 * Q03.z7 + 16 * ((π ^ 2 / 6) * Z02.z5)
        + 56 / 5 * ((π ^ 2 / 6) ^ 2 * Y04.z3)) * Polynomial.X ^ 2 := by
  refine ⟨poly_eq_of_sq _ _ (fun k => ?_), poly_eq_of_sq _ _ (fun k => ?_),
    poly_eq_of_sq _ _ (fun k => ?_), poly_eq_of_sq _ _ (fun k => ?_),
    poly_eq_of_sq _ _ (fun k => ?_), poly_eq_of_sq _ _ (fun k => ?_)⟩
  · rw [th2]; simp
  · rw [th3]; simp
  · rw [th4]; simp
  · rw [th5]; simp; ring
  · rw [th6]; simp; ring
  · rw [th7]; simp; ring

end R05

theorem R05.proof_8 : R05.stmt_8 := R05.theta_values

open MeasureTheory Set
open scoped ENNReal

namespace R05

/-- The simplex `x > t_0 > ⋯ > t_{n-1} > 0`. -/
def Sx (x : ℝ) (n : ℕ) : Set (Fin n → ℝ) := {t | (∀ i, 0 < t i ∧ t i < x) ∧ StrictAnti t}

lemma Sx_meas (x : ℝ) (n : ℕ) : MeasurableSet (Sx x n) := by
  have e : Sx x n = (⋂ i, ({t : Fin n → ℝ | 0 < t i} ∩ {t | t i < x})) ∩
      ⋂ i, ⋂ k, ⋂ (_ : i < k), {t : Fin n → ℝ | t k < t i} := by
    ext t; simp [Sx, StrictAnti, forall_and]
  rw [e]
  refine MeasurableSet.inter (MeasurableSet.iInter fun i => (measurableSet_lt measurable_const
    (measurable_pi_apply i)).inter (measurableSet_lt (measurable_pi_apply i) measurable_const))
    (MeasurableSet.iInter fun i => MeasurableSet.iInter fun k => MeasurableSet.iInter fun _ =>
      measurableSet_lt (measurable_pi_apply k) (measurable_pi_apply i))

lemma mem_cons (n : ℕ) (x a : ℝ) (s : Fin n → ℝ) :
    (Fin.cons a s : Fin (n + 1) → ℝ) ∈ Sx x (n + 1) ↔ a ∈ Ioo 0 x ∧ s ∈ Sx a n := by
  constructor
  · rintro ⟨h1, h2⟩
    have ha := h1 0
    simp only [Fin.cons_zero] at ha
    refine ⟨ha, fun i => ?_, fun i k hik => ?_⟩
    · have := h1 i.succ
      have h3 := h2 (Fin.succ_pos i)
      simp only [Fin.cons_succ, Fin.cons_zero] at this h3
      exact ⟨this.1, h3⟩
    · have := h2 (Fin.succ_lt_succ_iff.mpr hik)
      simpa using this
  · rintro ⟨ha, hs1, hs2⟩
    refine ⟨fun i => ?_, fun i k hik => ?_⟩
    · induction i using Fin.cases with
      | zero => simpa using ha
      | succ i => simp only [Fin.cons_succ]; exact ⟨(hs1 i).1, (hs1 i).2.trans ha.2⟩
    · induction k using Fin.cases with
      | zero => exact absurd hik (Fin.not_lt_zero _)
      | succ k =>
        induction i using Fin.cases with
        | zero => simp only [Fin.cons_succ, Fin.cons_zero]; exact (hs1 k).2
        | succ i =>
          simp only [Fin.cons_succ]
          exact hs2 (Fin.succ_lt_succ_iff.mp hik)

/-- Ordered integral with position weights. -/
noncomputable def Om (n : ℕ) (ws : ℕ → ℝ → ℝ≥0∞) (x : ℝ) : ℝ≥0∞ :=
  ∫⁻ t in Sx x n, ∏ i : Fin n, ws i.val (t i)

lemma Om_meas_integrand (n : ℕ) (ws : ℕ → ℝ → ℝ≥0∞) (hws : ∀ p, Measurable (ws p)) :
    Measurable (fun t : Fin n → ℝ => ∏ i : Fin n, ws i.val (t i)) :=
  Finset.measurable_prod _ (fun i _ => (hws i.val).comp (measurable_pi_apply i))

lemma Om_succ (n : ℕ) (ws : ℕ → ℝ → ℝ≥0∞) (hws : ∀ p, Measurable (ws p)) (x : ℝ) :
    Om (n + 1) ws x = ∫⁻ a in Ioo 0 x, ws 0 a * Om n (fun p => ws (p + 1)) a := by
  have hmp := (volume_preserving_piFinSuccAbove (fun _ : Fin (n + 1) => ℝ) 0).symm
  set F : (Fin (n + 1) → ℝ) → ℝ≥0∞ :=
    (Sx x (n + 1)).indicator (fun t => ∏ i : Fin (n + 1), ws i.val (t i)) with hF
  have hFm : Measurable F := (Om_meas_integrand _ ws hws).indicator (Sx_meas _ _)
  have h1 : Om (n + 1) ws x = ∫⁻ t, F t := by
    rw [Om, hF, lintegral_indicator (Sx_meas _ _)]
  rw [h1, ← hmp.lintegral_comp_emb (MeasurableEquiv.measurableEmbedding _)]
  rw [show (volume : Measure (ℝ × (Fin n → ℝ))) = volume.prod volume from rfl,
    lintegral_prod (fun z => F ((MeasurableEquiv.piFinSuccAbove (fun _ : Fin (n + 1) => ℝ) 0).symm z))
      (hFm.comp (MeasurableEquiv.measurable _)).aemeasurable]
  rw [← lintegral_indicator measurableSet_Ioo]
  refine lintegral_congr (fun a => ?_)
  simp only [MeasurableEquiv.piFinSuccAbove_symm_apply, Fin.insertNthEquiv,
    Equiv.coe_fn_mk, Fin.insertNth_zero']
  have hpt : ∀ s : Fin n → ℝ, F (Fin.cons a s) = (Ioo 0 x).indicator (fun a => ws 0 a *
      (Sx a n).indicator (fun s => ∏ i : Fin n, ws (i.val + 1) (s i)) s) a := by
    intro s
    rw [hF]
    by_cases hc : (Fin.cons a s : Fin (n + 1) → ℝ) ∈ Sx x (n + 1)
    · rw [indicator_of_mem hc]
      obtain ⟨ha, hs⟩ := (mem_cons n x a s).mp hc
      rw [indicator_of_mem ha, indicator_of_mem hs, Fin.prod_univ_succ]
      simp
    · rw [indicator_of_notMem hc]
      by_cases ha : a ∈ Ioo 0 x
      · rw [indicator_of_mem ha]
        have hs : s ∉ Sx a n := fun hs => hc ((mem_cons n x a s).mpr ⟨ha, hs⟩)
        rw [indicator_of_notMem hs, mul_zero]
      · rw [indicator_of_notMem ha]
  simp_rw [hpt]
  by_cases ha : a ∈ Ioo 0 x
  · simp only [indicator_of_mem ha]
    rw [lintegral_const_mul _ ((Om_meas_integrand n (fun p => ws (p + 1)) (fun p => hws _)).indicator
      (Sx_meas _ _)), lintegral_indicator (Sx_meas _ _)]
    rfl
  · simp [indicator_of_notMem ha]


lemma Om_zero (ws : ℕ → ℝ → ℝ≥0∞) (x : ℝ) : Om 0 ws x = 1 := by
  have : Sx x 0 = univ := by
    ext t; simp only [Sx, mem_univ, iff_true]
    exact ⟨fun i => i.elim0, fun i => i.elim0⟩
  rw [Om, this]; simp [volume_pi]

lemma ftc_l {g G : ℝ → ℝ} {lo hi : ℝ} (hle : lo ≤ hi)
    (hG : ∀ a ∈ Icc lo hi, HasDerivAt G (g a) a) (hg : ContinuousOn g (Icc lo hi))
    (hnn : ∀ a ∈ Icc lo hi, 0 ≤ g a) :
    ∫⁻ a in Ioo lo hi, ENNReal.ofReal (g a) = ENNReal.ofReal (G hi - G lo) := by
  rw [setLIntegral_congr Ioo_ae_eq_Ioc]
  have hint : IntegrableOn g (Ioc lo hi) := (hg.integrableOn_Icc).mono_set Ioc_subset_Icc_self
  rw [← ofReal_integral_eq_lintegral_ofReal hint
    ((ae_restrict_iff' measurableSet_Ioc).2 (Filter.Eventually.of_forall
      (fun a ha => hnn a (Ioc_subset_Icc_self ha))))]
  congr 1
  rw [← intervalIntegral.integral_of_le hle]
  refine intervalIntegral.integral_eq_sub_of_hasDerivAt (fun a ha => hG a ?_)
    (hg.intervalIntegrable_of_Icc hle)
  rwa [uIcc_of_le hle] at ha

/-- Weights of the word `(0^k, (-1)^∞)`, by position. -/
noncomputable def Wd (k p : ℕ) (t : ℝ) : ℝ≥0∞ :=
  ENNReal.ofReal (1 / (t - if p < k then 0 else -1))

lemma Wd_meas (k p : ℕ) : Measurable (Wd k p) := by
  unfold Wd; fun_prop

lemma Wd_shift (k : ℕ) : (fun p => Wd (k + 1) (p + 1)) = Wd k := by
  funext p t; simp [Wd]

lemma Wd_zero_shift : (fun p => Wd 0 (p + 1)) = Wd 0 := by
  funext p t; simp [Wd]

/-- `L_j(y) = log(1+y)^j / j!`. -/
noncomputable def Lj (j : ℕ) (y : ℝ) : ℝ := Real.log (1 + y) ^ j / (j.factorial : ℝ)

lemma Lj_nonneg (j : ℕ) {y : ℝ} (hy : 0 ≤ y) : 0 ≤ Lj j y := by
  unfold Lj; have := Real.log_nonneg (by linarith : (1 : ℝ) ≤ 1 + y); positivity

lemma inner_val (j : ℕ) : ∀ y : ℝ, 0 < y → Om j (Wd 0) y = ENNReal.ofReal (Lj j y) := by
  induction j with
  | zero => intro y _; rw [Om_zero]; simp [Lj]
  | succ j ih =>
    intro y hy
    rw [Om_succ j (Wd 0) (Wd_meas 0) y, Wd_zero_shift]
    rw [setLIntegral_congr_fun measurableSet_Ioo (g := fun a => ENNReal.ofReal
      (1 / (a + 1) * Lj j a)) (fun a ha => by
        rw [ih a ha.1, Wd]; simp only [Nat.not_lt_zero, ↓reduceIte, sub_neg_eq_add]
        rw [ENNReal.ofReal_mul (by have := ha.1; positivity)])]
    have hG : ∀ a ∈ Icc 0 y, HasDerivAt (Lj (j + 1)) (1 / (a + 1) * Lj j a) a := by
      intro a ha
      have h1 : (1 + a) ≠ 0 := by linarith [ha.1]
      have h := ((((hasDerivAt_id' a).const_add 1).log h1).pow (j + 1)).div_const
        (((j + 1).factorial : ℕ) : ℝ)
      show HasDerivAt (fun y => Real.log (1 + y) ^ (j + 1) / (((j + 1).factorial : ℕ) : ℝ)) _ a
      refine h.congr_deriv ?_
      simp only [Nat.factorial_succ, Nat.cast_mul, Nat.add_sub_cancel, Lj]
      have h2 : ((j.factorial : ℕ) : ℝ) ≠ 0 := by positivity
      have h3 : a + 1 ≠ 0 := by linarith [ha.1]
      field_simp
      ring
    rw [ftc_l hy.le hG ?_ ?_]
    · simp [Lj]
    · intro a ha
      have : ContinuousAt (fun a : ℝ => 1 / (a + 1) * Lj j a) a := by
        unfold Lj
        have h1 : (1 + a) ≠ 0 := by linarith [ha.1]
        have h2 : (a + 1) ≠ 0 := by linarith [ha.1]
        fun_prop (disch := assumption)
      exact this.continuousWithinAt
    · intro a ha
      have := Lj_nonneg j ha.1
      have : 0 < a + 1 := by linarith [ha.1]
      positivity


lemma tri_swap (x : ℝ) (f : ℝ → ℝ → ℝ≥0∞) (hf : Measurable (Function.uncurry f)) :
    ∫⁻ a in Ioo 0 x, ∫⁻ y in Ioo 0 a, f a y = ∫⁻ y in Ioo 0 x, ∫⁻ a in Ioo y x, f a y := by
  set T : Set (ℝ × ℝ) := {p | 0 < p.2 ∧ p.2 < p.1 ∧ p.1 < x} with hT
  have hTm : MeasurableSet T := by
    have e : T = {p : ℝ × ℝ | 0 < p.2} ∩ ({p | p.2 < p.1} ∩ {p | p.1 < x}) := by
      ext p; simp [hT]
    rw [e]
    exact (measurableSet_lt measurable_const measurable_snd).inter
      ((measurableSet_lt measurable_snd measurable_fst).inter (measurableSet_lt measurable_fst measurable_const))
  set G : ℝ → ℝ → ℝ≥0∞ := fun a y => T.indicator (Function.uncurry f) (a, y) with hG
  have hGm : Measurable (Function.uncurry G) := hf.indicator hTm
  have hL : ∫⁻ a in Ioo 0 x, ∫⁻ y in Ioo 0 a, f a y = ∫⁻ a, ∫⁻ y, G a y := by
    rw [← lintegral_indicator measurableSet_Ioo]
    refine lintegral_congr (fun a => ?_)
    by_cases ha : a ∈ Ioo 0 x
    · rw [indicator_of_mem ha, ← lintegral_indicator measurableSet_Ioo]
      refine lintegral_congr (fun y => ?_)
      by_cases hy : y ∈ Ioo 0 a
      · rw [indicator_of_mem hy, hG]; simp only
        rw [indicator_of_mem (show (a, y) ∈ T from ⟨hy.1, hy.2, ha.2⟩)]; rfl
      · rw [indicator_of_notMem hy, hG]; simp only
        rw [indicator_of_notMem (show (a, y) ∉ T from fun h => hy ⟨h.1, h.2.1⟩)]
    · rw [indicator_of_notMem ha]
      symm; refine lintegral_eq_zero_of_ae_eq_zero (Filter.Eventually.of_forall (fun y => ?_))
      rw [hG]; simp only
      rw [indicator_of_notMem (show (a, y) ∉ T from fun h => ha ⟨h.1.trans h.2.1, h.2.2⟩)]; rfl
  have hR : ∫⁻ y in Ioo 0 x, ∫⁻ a in Ioo y x, f a y = ∫⁻ y, ∫⁻ a, G a y := by
    rw [← lintegral_indicator measurableSet_Ioo]
    refine lintegral_congr (fun y => ?_)
    by_cases hy : y ∈ Ioo 0 x
    · rw [indicator_of_mem hy, ← lintegral_indicator measurableSet_Ioo]
      refine lintegral_congr (fun a => ?_)
      by_cases ha : a ∈ Ioo y x
      · rw [indicator_of_mem ha, hG]; simp only
        rw [indicator_of_mem (show (a, y) ∈ T from ⟨hy.1, ha.1, ha.2⟩)]; rfl
      · rw [indicator_of_notMem ha, hG]; simp only
        rw [indicator_of_notMem (show (a, y) ∉ T from fun h => ha ⟨h.2.1, h.2.2⟩)]
    · rw [indicator_of_notMem hy]
      symm; refine lintegral_eq_zero_of_ae_eq_zero (Filter.Eventually.of_forall (fun a => ?_))
      rw [hG]; simp only
      rw [indicator_of_notMem (show (a, y) ∉ T from fun h => hy ⟨h.1, h.2.1.trans h.2.2⟩)]; rfl
  rw [hL, hR, lintegral_lintegral_swap hGm.aemeasurable]

lemma zeros_val (j k : ℕ) : ∀ x : ℝ, 0 < x → Om (k + 1 + j) (Wd (k + 1)) x =
    ∫⁻ y in Ioo 0 x, ENNReal.ofReal (Real.log (x / y) ^ k / (k.factorial : ℝ) * (1 / y) * Lj j y) := by
  induction k with
  | zero =>
    intro x _
    rw [show 0 + 1 + j = j + 1 by omega, Om_succ j _ (Wd_meas 1) x, Wd_shift 0]
    refine setLIntegral_congr_fun measurableSet_Ioo (fun a ha => ?_)
    rw [inner_val j a ha.1, Wd]
    simp only [zero_lt_one, ↓reduceIte, sub_zero, pow_zero, Nat.factorial_zero, Nat.cast_one,
      div_one, one_mul]
    rw [ENNReal.ofReal_mul (by have := ha.1; positivity)]
  | succ k ih =>
    intro x hx
    rw [show k + 1 + 1 + j = (k + 1 + j) + 1 by omega, Om_succ _ _ (Wd_meas _) x, Wd_shift (k + 1)]
    have e1 : ∀ a ∈ Ioo 0 x, Wd (k + 1 + 1) 0 a * Om (k + 1 + j) (Wd (k + 1)) a =
        ∫⁻ y in Ioo 0 a, ENNReal.ofReal (1 / a) *
          ENNReal.ofReal (Real.log (a / y) ^ k / (k.factorial : ℝ) * (1 / y) * Lj j y) := by
      intro a ha
      rw [ih a ha.1, Wd, ← lintegral_const_mul' _ _ ENNReal.ofReal_ne_top]
      simp
    rw [setLIntegral_congr_fun measurableSet_Ioo e1]
    rw [tri_swap x (fun a y => ENNReal.ofReal (1 / a) *
          ENNReal.ofReal (Real.log (a / y) ^ k / (k.factorial : ℝ) * (1 / y) * Lj j y)) (by
      unfold Lj; fun_prop)]
    refine setLIntegral_congr_fun measurableSet_Ioo (fun y hy => ?_)
    have hy0 := hy.1
    have hcy : 0 ≤ 1 / y * Lj j y := by have := Lj_nonneg j hy0.le; positivity
    have e2 : ∀ a ∈ Ioo y x, ENNReal.ofReal (1 / a) *
          ENNReal.ofReal (Real.log (a / y) ^ k / (k.factorial : ℝ) * (1 / y) * Lj j y) =
        ENNReal.ofReal (1 / y * Lj j y) *
          ENNReal.ofReal (Real.log (a / y) ^ k / (k.factorial : ℝ) * (1 / a)) := by
      intro a ha
      have ha0 : 0 < a := hy0.trans ha.1
      rw [← ENNReal.ofReal_mul (by positivity), ← ENNReal.ofReal_mul hcy]
      congr 1; ring
    rw [setLIntegral_congr_fun measurableSet_Ioo e2, lintegral_const_mul' _ _ ENNReal.ofReal_ne_top]
    have hG : ∀ a ∈ Icc y x, HasDerivAt (fun a => Real.log (a / y) ^ (k + 1) / ((k + 1).factorial : ℝ))
        (Real.log (a / y) ^ k / (k.factorial : ℝ) * (1 / a)) a := by
      intro a ha
      have ha0 : 0 < a := hy0.trans_le ha.1
      have h1 : a / y ≠ 0 := by positivity
      have h := ((((hasDerivAt_id' a).div_const y).log h1).pow (k + 1)).div_const
        (((k + 1).factorial : ℕ) : ℝ)
      refine h.congr_deriv ?_
      simp only [Nat.factorial_succ, Nat.cast_mul, Nat.add_sub_cancel]
      have h2 : ((k.factorial : ℕ) : ℝ) ≠ 0 := by positivity
      field_simp
    rw [ftc_l hy.2.le hG ?_ ?_]
    · rw [← ENNReal.ofReal_mul hcy]
      congr 1
      simp only [div_self hy0.ne', Real.log_one]
      ring
    · intro a ha
      have ha0 : 0 < a := hy0.trans_le ha.1
      have : ContinuousAt (fun a : ℝ => Real.log (a / y) ^ k / (k.factorial : ℝ) * (1 / a)) a := by
        have h1 : a / y ≠ 0 := by positivity
        have h2 : a ≠ 0 := ha0.ne'
        fun_prop (disch := assumption)
      exact this.continuousWithinAt
    · intro a ha
      have ha0 : 0 < a := hy0.trans_le ha.1
      have : 0 ≤ Real.log (a / y) := Real.log_nonneg (by rw [le_div_iff₀ hy0]; linarith [ha.1])
      positivity


lemma cov (m j : ℕ) :
    ∫⁻ y in Ioo 0 1, ENNReal.ofReal (Real.log (1 / y) ^ m / (m.factorial : ℝ) * (1 / y) * Lj j y) =
      ∫⁻ s in Ioi 0, ENNReal.ofReal (s ^ m * Y04.ell s ^ j / ((m.factorial : ℝ) * (j.factorial : ℝ))) := by
  have himg : (fun s : ℝ => Real.exp (-s)) '' Ioi 0 = Ioo 0 1 := by
    ext y
    constructor
    · rintro ⟨s, hs, rfl⟩
      exact ⟨Real.exp_pos _, by simp only [mem_Ioi] at hs; show Real.exp (-s) < 1; rw [Real.exp_lt_one_iff]; linarith⟩
    · rintro ⟨h0, h1⟩
      refine ⟨-Real.log y, ?_, ?_⟩
      · simp only [mem_Ioi]; linarith [Real.log_neg h0 h1]
      · simp only [neg_neg, Real.exp_log h0]
  rw [← himg, lintegral_image_eq_lintegral_abs_deriv_mul measurableSet_Ioi
    (f' := fun s => -Real.exp (-s))
    (fun s _ => by
      have := (hasDerivAt_neg s).exp
      simpa using this.hasDerivWithinAt)
    (fun a _ b _ h => by simpa using h)]
  refine setLIntegral_congr_fun measurableSet_Ioi (fun s hs => ?_)
  rw [← ENNReal.ofReal_mul (abs_nonneg _)]
  congr 1
  have he : Real.log (1 / Real.exp (-s)) = s := by
    rw [one_div, Real.log_inv, Real.log_exp, neg_neg]
  rw [he, abs_neg, abs_of_pos (Real.exp_pos _)]
  simp only [Lj, Y04.ell]
  have h1 : Real.exp (-s) ≠ 0 := (Real.exp_pos _).ne'
  field_simp

lemma lint_I (m j : ℕ) (hj : 1 ≤ j) :
    ∫⁻ s in Ioi 0, ENNReal.ofReal (s ^ m * Y04.ell s ^ j / ((m.factorial : ℝ) * (j.factorial : ℝ))) =
      ENNReal.ofReal (Y04.I m j / ((m.factorial : ℝ) * (j.factorial : ℝ))) := by
  rw [← ofReal_integral_eq_lintegral_ofReal ((Y04.I_int m j hj).div_const _)
    ((ae_restrict_iff' measurableSet_Ioi).2 (Filter.Eventually.of_forall (fun s (hs : 0 < s) => by
      have := Y04.ell_nonneg s; have := hs.le; positivity)))]
  rw [integral_div]; rfl

lemma word_pos (m j : ℕ) (t : Fin (m + 1 + j) → ℝ) (ht : t ∈ simplex (m + 1 + j)) (i : Fin (m + 1 + j)) :
    0 < 1 / (t i - word m j i) := by
  have := (ht.1 i).1
  unfold word
  split_ifs <;> (apply one_div_pos.mpr; linarith)

lemma atoms (m j : ℕ) (hj : 1 ≤ j) :
    IntegrableOn (fun t : Fin (m + 1 + j) → ℝ => ∏ i, 1 / (t i - word m j i))
      (simplex (m + 1 + j)) ∧
    Y04.I m j = (m.factorial : ℝ) * (j.factorial : ℝ) * iterInt (m + 1 + j) (word m j) := by
  set n := m + 1 + j
  set f : (Fin (m + 1 + j) → ℝ) → ℝ := fun t => ∏ i, 1 / (t i - word m j i) with hf
  have hfm : Measurable f := by
    rw [hf]
    exact Finset.measurable_prod _ (fun i _ => by fun_prop)
  have hsm : MeasurableSet (simplex (m + 1 + j)) := Sx_meas 1 _
  have hnn : 0 ≤ᵐ[volume.restrict (simplex (m + 1 + j))] f :=
    (ae_restrict_iff' hsm).2 (Filter.Eventually.of_forall (fun t ht =>
      Finset.prod_nonneg (fun i _ => (word_pos m j t ht i).le)))
  have hlin : ∫⁻ t in simplex (m + 1 + j), ENNReal.ofReal (f t) = Om (m + 1 + j) (Wd (m + 1)) 1 := by
    refine setLIntegral_congr_fun hsm (fun t ht => ?_)
    rw [hf, ENNReal.ofReal_prod_of_nonneg (fun i _ => (word_pos m j t ht i).le)]
    rfl
  have hval : Om (m + 1 + j) (Wd (m + 1)) 1 =
      ENNReal.ofReal (Y04.I m j / ((m.factorial : ℝ) * (j.factorial : ℝ))) := by
    rw [zeros_val j m 1 one_pos, cov, lint_I m j hj]
  refine ⟨⟨hfm.aestronglyMeasurable, (hasFiniteIntegral_iff_ofReal hnn).mpr (by
    rw [hlin, hval]; exact ENNReal.ofReal_lt_top)⟩, ?_⟩
  have hI : 0 ≤ Y04.I m j := (I_pos m j hj).le
  have hint : iterInt (m + 1 + j) (word m j) = Y04.I m j / ((m.factorial : ℝ) * (j.factorial : ℝ)) := by
    unfold iterInt
    rw [integral_eq_lintegral_of_nonneg_ae hnn hfm.aestronglyMeasurable, hlin, hval,
      ENNReal.toReal_ofReal (by positivity)]
  rw [hint]
  have h1 : (m.factorial : ℝ) ≠ 0 := by positivity
  have h2 : (j.factorial : ℝ) ≠ 0 := by positivity
  field_simp

end R05

theorem R05.proof_10 : R05.stmt_10 := fun m j hj => R05.atoms m j hj
