import Mathlib
import Zeta3.Wave8.Q01_PureWeight
import Zeta3.Wave9.R02_Theta8
import Zeta3.Wave9.R05_DegreeAllN
/-!
# S01 — `Θ₉` in closed form, and `prop:theta` for all `5 ≤ n ≤ 9`

Target id: S01.
Paper labels: `prop:theta` and `prop:degree` (main.tex lines 1214–1225).
Detailed proofs: paper/ieee/proofs/sec_g2_theta.tex (the display of `Θ₉` at lines 21–23,
Lemma `g2:lem:t:red` at line 33, Table `g2:tab:D` rows `D̃_{9,0}, …, D̃_{9,3}` at lines 59–62, the
printed values of `e_{r,m}` in "Output") and paper/ieee/proofs/sec_g2_degree.tex
(`e_{r,m}`, `ρ_m`, Lemmas `g2:lem:d:e`, `g2:lem:d:theta`).
Earlier Lean (imported, not redefined): `Zeta3.Wave5.Y04_Fourier` (`Y04.tau`, `Y04.h`, `Y04.a`,
`Y04.G`, `Y04.I`, `Y04.z3`), `Zeta3.Wave6.Z02_Theta5` (`Z02.z5`), `Zeta3.Wave8.Q03_Theta7`
(`Q03.z7`), `Zeta3.Wave4.X05_DiamondA4A5` (`X05.z n = (riemannZeta n).re`, used for `ζ(9)`),
`Zeta3.Wave8.Q01_PureWeight` (`Q01.IsZetaPoly`), `Zeta3.Wave9.R05_DegreeAllN` (`R05.e`, `R05.rho`,
`R05.Theta`; R05 proves Watson for every `N`, `g_{k,m} = (-1)^k(√2/4)ρ_m(k²)`, the expansion of
`a_k(h_p)` to every order with coefficients `R05.Theta n`, `natDegree Θ_n = ⌊n/2⌋ - 1` for every
`n ≥ 2` (so "degree exactly 3" for `Θ₉` is `R05.proof_9` and is not restated here), `I_{m,1}`,
and `Θ₂, …, Θ₇` in closed form), `Zeta3.Wave9.R02_Theta8` (`Θ₈` in closed form, rows `D̃_{8,*}`).

Quoted verbatim.

main.tex 1214–1216 (`prop:degree`):

  For $n\ge2$, the coefficient of $\tau^n$ in $(-1)^k\frac{\pi}{\sqrt2}a_k(h_p)$ is a polynomial
  $\Theta_n(k^2)$ of degree exactly $\lfloor n/2\rfloor-1$, whose coefficients are alternating
  multiple zeta values of weight $n$. For even $n=2m$ its leading coefficient is
  $(-1)^{m-1}(2^{2m-1}-1)\zeta(2m)$.

main.tex 1220–1225 (`prop:theta`):

  For $5\le n\le9$, every coefficient of $\Theta_n$ is a rational combination of products of
  single zeta values of total weight $n$. With $u=k^2$,
  \Theta_5=-\bigl(\tfrac34+4u\bigr)\zeta(5)-\bigl(\tfrac74+4u\bigr)\zeta(2)\zeta(3).

Supplement sec_g2_theta.tex lines 21–23:

  \Theta_9&=-\bigl(\tfrac{22805}{192}+\tfrac{3907}4u+\tfrac{1916}3u^2+64u^3\bigr)\zeta(9)-\bigl(\tfrac{8343}{64}+\tfrac{4163}4u+660u^2+64u^3\bigr)\zeta(2)\zeta(7)\\
  &\qquad-\bigl(\tfrac{60093}{640}+\tfrac{29741}{40}u+\tfrac{2334}5u^2+\tfrac{224}5u^3\bigr)\zeta(2)^2\zeta(5)\\
  &\qquad-\bigl(\tfrac{263713}{4480}+\tfrac{130973}{280}u+\tfrac{2062}7u^2+\tfrac{992}{35}u^3\bigr)\zeta(2)^3\zeta(3)+\bigl(\tfrac{169}{96}+9u+\tfrac83u^2\bigr)\zeta(3)^3 .

Supplement, Lemma `g2:lem:t:red`:

  Let $\wD_{n,m}=2\sum_{j=1}^{n-1-2m}\dfrac{e_{m,\,n-1-j}}{j!}\,I_{n-1-j,\,j}$, for $n\ge2$ and
  $0\le m\le\lfloor n/2\rfloor-1$, with $e_{r,m}$ as in Section~\ref{g2:sec:degree}.
  $\Theta_n(u)=\sum_{m=0}^{\lfloor n/2\rfloor-1}\wD_{n,m}\,\dfrac{(-4u)^m}{(2m)!}$.

Supplement, Table `g2:tab:D` (products `ζ(9), ζ(2)ζ(7), ζ(2)²ζ(5), ζ(2)³ζ(3), ζ(3)³`, in this
order):

  $\wD_{9,0}$ & $-\tfrac{22805}{192},\ -\tfrac{8343}{64},\ -\tfrac{60093}{640},\ -\tfrac{263713}{4480},\ \tfrac{169}{96}$ \\
  $\wD_{9,1}$ & $\tfrac{3907}8,\ \tfrac{4163}8,\ \tfrac{29741}{80},\ \tfrac{130973}{560},\ -\tfrac92$ \\
  $\wD_{9,2}$ & $-958,\ -990,\ -\tfrac{3501}5,\ -\tfrac{3093}7,\ 4$ \\
  $\wD_{9,3}$ & $720,\ 720,\ 504,\ \tfrac{2232}7,\ 0$

Supplement, "Output" (values of `e_{m,r}`):

  $m=0$: $1,\frac12,-\frac58,-\frac{17}{48},\frac{121}{384}$; $m=1$ ($r\ge2$):
  $1,\frac12,-\frac{23}{24},-\frac{25}{48},\frac{731}{1152}$; $m=2$ ($r\ge4$):
  $1,\frac12,-\frac{31}{24},-\frac{11}{16},\frac{409}{384}$

The supplement does not print `e_{0,5}, e_{0,6}, e_{0,7}, e_{1,7}`. They were computed here
exactly from the definition `e_{r,m} = 2^{3/2}[y^m](gd(y)^{2r}H(y))` (SymPy series of
`2^{3/2}H(y) = e^{-y}((1+e^{-2y})/2)^{-3/2}` and of `gd(y) = π/2 - 2 arctan e^{-y}`; the printed
values were reproduced): `e_{0,5} = 721/3840`, `e_{0,6} = -1369/9216` (as in Q03, R02),
`e_{0,7} = -58337/645120`, `e_{1,7} = 4163/11520`. `e_{3,6} = 1`, `e_{3,7} = 1/2` are
`e_{r,2r} = 1`, `e_{r,2r+1} = 1/2` (Lemma `g2:lem:d:e`(a), `R05.proof_4`). So
`ρ₇(u) = e_{0,7} - 2u e_{1,7} + (2/3)u² e_{2,7} - (4/45)u³ e_{3,7}
  = -58337/645120 - (4163/5760)u - (11/24)u² - (2/45)u³`.

Conventions.
* `ζ(3) = Y04.z3`, `ζ(5) = Z02.z5`, `ζ(7) = Q03.z7`, `ζ(9) = X05.z 9`: each is the real part of
  Mathlib's `riemannZeta` at that natural number (by `zeta_nat_eq_tsum_of_gt_one` the real series
  `∑ 1/n^s`). `ζ(2) = π²/6`, `ζ(4) = π⁴/90`, `ζ(6) = π⁶/945`, `ζ(8) = π⁸/9450` are written out as
  in Y04, V02, R02 (Mathlib `riemannZeta_two`, …; `R02.proof_1` for `ζ(8)`). So
  `ζ(2)ζ(7) = (π²/6)·Q03.z7`, `ζ(2)²ζ(5) = (π²/6)²·Z02.z5`, `ζ(2)³ζ(3) = (π²/6)³·Y04.z3`,
  `ζ(3)³ = Y04.z3³`.
* "Rational combination of products of single zeta values of total weight `n`" is
  `Q01.IsZetaPoly n` (a finite list of monomials `c·ζ(n_1)⋯ζ(n_r)`, `c ∈ ℚ`, `n_i ≥ 2`,
  `Σ n_i = n`, with `ζ(j) = X05.z j`). Weight is a property of the representation, as in Q01.
* The paper's `Θ_n` is `R05.Theta n`, defined by the supplement's formula
  `Θ_n = 2 Σ_j ρ_{n-1-j} I_{n-1-j,j}/j!`; `R05.proof_9` (i) proves that `Θ_n(k²)` is the `τ^n`
  coefficient of `(-1)^k(π/√2)a_k(h_p)` (expansion with explicit `C τ^{N+1}` error, every `N`).
* `g_{k,7} = G_k^{(7)}(0)/7!` is `iteratedDeriv 7 (Y04.G k) 0 / 5040`; `G_k` is `C^∞` (Y04).
* All `I_{m,j}`, `j ≥ 1`, have integrable integrands (Y04), so no statement holds by the
  "non-integrable ↦ 0" convention.
* The `O(τ^{10})` claim is explicit: `∃ C τ₀, 0 < τ₀ ∧ ∀ p > 1, τ ≤ τ₀ → |error| ≤ C τ^{10}`.
* Numerical check before freezing (mpmath, 40 digits; SymPy for `e_{r,m}` and for Taylor
  coefficients of `G_k`): the four rows `D̃_{9,0}, …, D̃_{9,3}` (`stmt_5`–`stmt_8`, values
  `-910.5617…, 3639.0415…, -6923.6550…, 5035.7251…`) agree with the table's combinations to
  `≤ 1.1·10⁻³²`; the table substituted into Lemma `g2:lem:t:red` gives exactly the displayed `Θ₉`
  (exact rational arithmetic); `g_{k,7}` (`stmt_3`) agrees exactly with the SymPy Taylor coefficient
  of `G_k` for `k = 0, 1, 2, 3`; the `τ⁹` coefficient `Σ_m D̃_{9,m}(-4k²)^m/(2m)!` minus the closed
  form `Θ₉(k²)` is `≤ 1.3·10⁻³⁰` in absolute value for `k = 0, 1, 2, 3` (values `-910.5617…,
  -13252.0349…, -132522.8953…, -766605.6703…`). Direct quadrature of `a_k(h_p)` was not redone.

NOT COVERED.
* The clause "whose coefficients are alternating multiple zeta values of weight `n`" of
  `prop:degree` beyond what `R05.stmt_6`, `R05.stmt_10` state (the general AMZV predicate is not
  defined).
* The certificate route of the supplement (shuffle, stuffle, duality relations and the sixteen
  exact certificates, Lemma `g2:lem:t:euler`): the sixteen identities are stated directly
  (`R02`, earlier files, and `stmt_5`–`stmt_8` here for `n = 9`), whatever the proof route.
* The display of `Θ₅` in `prop:theta` and `Θ₆, Θ₇` are `R05.proof_8`, `Θ₈` is `R02.proof_8`;
  "degree exactly `⌊n/2⌋ - 1`" for `n = 9` is `R05.proof_9`; they are not restated, only used.
* The remark that for larger `n` the question is open; uniqueness of asymptotic coefficients.
-/
-- BEGIN STATEMENTS
namespace S01

open Real MeasureTheory Set Polynomial

/-- `D̃_{n,m} = 2 Σ_{j=1}^{n-1-2m} e_{m,n-1-j}/j! I_{n-1-j,j}` (supplement before Lemma
`g2:lem:t:red`). -/
noncomputable def Dt (n m : ℕ) : ℝ :=
  2 * ∑ j ∈ Finset.Icc 1 (n - 1 - 2 * m),
    R05.e m (n - 1 - j) / (j.factorial : ℝ) * Y04.I (n - 1 - j) j

/-- The values of `e_{r,m}` printed in the supplement ("Output") that `D̃_{9,*}` needs, besides
`e_{r,2r} = 1`, `e_{r,2r+1} = 1/2` (`R05.stmt_4`): `e_{0,2} = -5/8`, `e_{0,3} = -17/48`,
`e_{0,4} = 121/384`, `e_{1,4} = -23/24`, `e_{1,5} = -25/48`, `e_{1,6} = 731/1152`,
`e_{2,6} = -31/24`, `e_{2,7} = -11/16`. -/
def stmt_1 : Prop :=
  R05.e 0 2 = -(5 / 8) ∧ R05.e 0 3 = -(17 / 48) ∧ R05.e 0 4 = 121 / 384 ∧
  R05.e 1 4 = -(23 / 24) ∧ R05.e 1 5 = -(25 / 48) ∧ R05.e 1 6 = 731 / 1152 ∧
  R05.e 2 6 = -(31 / 24) ∧ R05.e 2 7 = -(11 / 16)

/-- The values of `e_{r,m}`, `m ≤ 7`, that the supplement does not print (computed here from the
definition): `e_{0,5} = 721/3840`, `e_{0,6} = -1369/9216`, `e_{0,7} = -58337/645120`,
`e_{1,7} = 4163/11520`. -/
def stmt_2 : Prop :=
  R05.e 0 5 = 721 / 3840 ∧ R05.e 0 6 = -(1369 / 9216) ∧
  R05.e 0 7 = -(58337 / 645120) ∧ R05.e 1 7 = 4163 / 11520

/-- `ρ₇(u) = -58337/645120 - (4163/5760)u - (11/24)u² - (2/45)u³`, and so (Lemma `g2:lem:d:e`(b)
at `m = 7`) the Taylor coefficient `g_{k,7} = G_k^{(7)}(0)/7! = (-1)^k (√2/4) ρ₇(k²)`. -/
def stmt_3 : Prop :=
  R05.rho 7 = Polynomial.C (-(58337 / 645120)) + Polynomial.C (-(4163 / 5760)) * X
      + Polynomial.C (-(11 / 24)) * X ^ 2 + Polynomial.C (-(2 / 45)) * X ^ 3 ∧
  ∀ k : ℕ, iteratedDeriv 7 (Y04.G k) 0 / 5040 =
    (-1 : ℝ) ^ k * Real.sqrt 2 / 4 *
      (-(58337 / 645120) - 4163 / 5760 * (k : ℝ) ^ 2 - 11 / 24 * (k : ℝ) ^ 4
        - 2 / 45 * (k : ℝ) ^ 6)

/-- Lemma `g2:lem:t:red`, for every `n ≥ 2`:
`Θ_n(u) = Σ_{m=0}^{⌊n/2⌋-1} D̃_{n,m} (-4u)^m/(2m)!` as polynomials in `u`. -/
def stmt_4 : Prop :=
  ∀ n : ℕ, 2 ≤ n →
    R05.Theta n = ∑ m ∈ Finset.range (n / 2),
      Polynomial.C (Dt n m * (-4 : ℝ) ^ m / ((2 * m).factorial : ℝ)) * X ^ m

/-- Table `g2:tab:D`, row `D̃_{9,3}` (with `e_{3,7} = 1/2`, `e_{3,6} = 1`):
`D̃_{9,3} = 2((1/2) I_{7,1} + I_{6,2}/2!) = 720 ζ(9) + 720 ζ(2)ζ(7) + 504 ζ(2)²ζ(5)
  + (2232/7) ζ(2)³ζ(3)`. -/
def stmt_5 : Prop :=
  2 * ((1 / 2) * Y04.I 7 1 + Y04.I 6 2 / 2) =
    720 * X05.z 9 + 720 * ((π ^ 2 / 6) * Q03.z7) + 504 * ((π ^ 2 / 6) ^ 2 * Z02.z5)
      + 2232 / 7 * ((π ^ 2 / 6) ^ 3 * Y04.z3)

/-- Table `g2:tab:D`, row `D̃_{9,2}` (with `e_{2,4..7} = 1, 1/2, -31/24, -11/16`):
`D̃_{9,2} = 2(e_{2,7} I_{7,1} + e_{2,6} I_{6,2}/2! + e_{2,5} I_{5,3}/3! + e_{2,4} I_{4,4}/4!)
  = -958 ζ(9) - 990 ζ(2)ζ(7) - (3501/5) ζ(2)²ζ(5) - (3093/7) ζ(2)³ζ(3) + 4 ζ(3)³`. -/
def stmt_6 : Prop :=
  2 * (-(11 / 16) * Y04.I 7 1 + (-(31 / 24)) * Y04.I 6 2 / 2 + (1 / 2) * Y04.I 5 3 / 6
      + Y04.I 4 4 / 24) =
    -958 * X05.z 9 - 990 * ((π ^ 2 / 6) * Q03.z7) - 3501 / 5 * ((π ^ 2 / 6) ^ 2 * Z02.z5)
      - 3093 / 7 * ((π ^ 2 / 6) ^ 3 * Y04.z3) + 4 * Y04.z3 ^ 3

/-- Table `g2:tab:D`, row `D̃_{9,1}` (with `e_{1,2..7} = 1, 1/2, -23/24, -25/48, 731/1152,
4163/11520`): `D̃_{9,1} = 2(e_{1,7} I_{7,1} + e_{1,6} I_{6,2}/2! + e_{1,5} I_{5,3}/3!
  + e_{1,4} I_{4,4}/4! + e_{1,3} I_{3,5}/5! + e_{1,2} I_{2,6}/6!)
  = (3907/8) ζ(9) + (4163/8) ζ(2)ζ(7) + (29741/80) ζ(2)²ζ(5) + (130973/560) ζ(2)³ζ(3)
    - (9/2) ζ(3)³`. -/
def stmt_7 : Prop :=
  2 * (4163 / 11520 * Y04.I 7 1 + 731 / 1152 * Y04.I 6 2 / 2 + (-(25 / 48)) * Y04.I 5 3 / 6
      + (-(23 / 24)) * Y04.I 4 4 / 24 + (1 / 2) * Y04.I 3 5 / 120 + Y04.I 2 6 / 720) =
    3907 / 8 * X05.z 9 + 4163 / 8 * ((π ^ 2 / 6) * Q03.z7)
      + 29741 / 80 * ((π ^ 2 / 6) ^ 2 * Z02.z5) + 130973 / 560 * ((π ^ 2 / 6) ^ 3 * Y04.z3)
      - 9 / 2 * Y04.z3 ^ 3

/-- Table `g2:tab:D`, row `D̃_{9,0}` (with `e_{0,0..7} = 1, 1/2, -5/8, -17/48, 121/384,
721/3840, -1369/9216, -58337/645120`): `D̃_{9,0} = 2 Σ_{j=1}^{8} e_{0,8-j} I_{8-j,j}/j!
  = -(22805/192) ζ(9) - (8343/64) ζ(2)ζ(7) - (60093/640) ζ(2)²ζ(5)
    - (263713/4480) ζ(2)³ζ(3) + (169/96) ζ(3)³`. -/
def stmt_8 : Prop :=
  2 * (-(58337 / 645120) * Y04.I 7 1 + (-(1369 / 9216)) * Y04.I 6 2 / 2
      + 721 / 3840 * Y04.I 5 3 / 6 + 121 / 384 * Y04.I 4 4 / 24
      + (-(17 / 48)) * Y04.I 3 5 / 120 + (-(5 / 8)) * Y04.I 2 6 / 720
      + (1 / 2) * Y04.I 1 7 / 5040 + Y04.I 0 8 / 40320) =
    -(22805 / 192) * X05.z 9 - 8343 / 64 * ((π ^ 2 / 6) * Q03.z7)
      - 60093 / 640 * ((π ^ 2 / 6) ^ 2 * Z02.z5) - 263713 / 4480 * ((π ^ 2 / 6) ^ 3 * Y04.z3)
      + 169 / 96 * Y04.z3 ^ 3

/-- The supplement's display of `Θ₉` (sec_g2_theta.tex lines 21–23), as an identity of
polynomials in `u` for the `Θ₉` of Lemma `g2:lem:d:theta` (`R05.Theta 9`):
`Θ₉ = -(22805/192 + (3907/4)u + (1916/3)u² + 64u³) ζ(9)
  - (8343/64 + (4163/4)u + 660u² + 64u³) ζ(2)ζ(7)
  - (60093/640 + (29741/40)u + (2334/5)u² + (224/5)u³) ζ(2)²ζ(5)
  - (263713/4480 + (130973/280)u + (2062/7)u² + (992/35)u³) ζ(2)³ζ(3)
  + (169/96 + 9u + (8/3)u²) ζ(3)³`. -/
def stmt_9 : Prop :=
  R05.Theta 9 =
    Polynomial.C (-(22805 / 192) * X05.z 9 - 8343 / 64 * ((π ^ 2 / 6) * Q03.z7)
        - 60093 / 640 * ((π ^ 2 / 6) ^ 2 * Z02.z5)
        - 263713 / 4480 * ((π ^ 2 / 6) ^ 3 * Y04.z3) + 169 / 96 * Y04.z3 ^ 3)
    + Polynomial.C (-(3907 / 4) * X05.z 9 - 4163 / 4 * ((π ^ 2 / 6) * Q03.z7)
        - 29741 / 40 * ((π ^ 2 / 6) ^ 2 * Z02.z5)
        - 130973 / 280 * ((π ^ 2 / 6) ^ 3 * Y04.z3) + 9 * Y04.z3 ^ 3) * X
    + Polynomial.C (-(1916 / 3) * X05.z 9 - 660 * ((π ^ 2 / 6) * Q03.z7)
        - 2334 / 5 * ((π ^ 2 / 6) ^ 2 * Z02.z5)
        - 2062 / 7 * ((π ^ 2 / 6) ^ 3 * Y04.z3) + 8 / 3 * Y04.z3 ^ 3) * X ^ 2
    + Polynomial.C (-64 * X05.z 9 - 64 * ((π ^ 2 / 6) * Q03.z7)
        - 224 / 5 * ((π ^ 2 / 6) ^ 2 * Z02.z5)
        - 992 / 35 * ((π ^ 2 / 6) ^ 3 * Y04.z3)) * X ^ 3

/-- `thm:fourier` extended through `τ⁹` with `prop:theta`'s `Θ₅` and the supplement's `Θ₆, …, Θ₉`:
for every `k ≥ 0` there are `C_k` and `τ₀ > 0` such that for every `p > 1` with
`τ = 1 - 1/p ≤ τ₀`, with `u = k²`,
`|a_k(h_p) - (-1)^k (√2/π)[-4/(16k²-1) + ζ(2)τ² + ζ(3)τ³ - (33/16 + 7k²)ζ(4)τ⁴ + Θ₅(u)τ⁵
   + Θ₆(u)τ⁶ + Θ₇(u)τ⁷ + Θ₈(u)τ⁸ + Θ₉(u)τ⁹]| ≤ C_k τ¹⁰`, where `Θ₅, …, Θ₈` are as in
`R02.stmt_10` and `Θ₉` is the closed form of `stmt_9`, all written out. -/
def stmt_10 : Prop :=
  ∀ k : ℕ, ∃ C τ₀ : ℝ, 0 < τ₀ ∧ ∀ p : ℝ, 1 < p → Y04.tau p ≤ τ₀ →
    |Y04.a k (Y04.h p) - (-1 : ℝ) ^ k * (Real.sqrt 2 / π) *
        (-(4 / (16 * (k : ℝ) ^ 2 - 1)) + (π ^ 2 / 6) * Y04.tau p ^ 2 + Y04.z3 * Y04.tau p ^ 3
          - (33 / 16 + 7 * (k : ℝ) ^ 2) * (π ^ 4 / 90) * Y04.tau p ^ 4
          + (-(3 / 4 + 4 * (k : ℝ) ^ 2) * Z02.z5
              - (7 / 4 + 4 * (k : ℝ) ^ 2) * ((π ^ 2 / 6) * Y04.z3)) * Y04.tau p ^ 5
          + ((4039 / 256 + 733 / 8 * (k : ℝ) ^ 2 + 31 * (k : ℝ) ^ 4) * (π ^ 6 / 945)
              - (7 / 8 + 2 * (k : ℝ) ^ 2) * Y04.z3 ^ 2) * Y04.tau p ^ 6
          + ((129 / 16 + 46 * (k : ℝ) ^ 2 + 16 * (k : ℝ) ^ 4) * Q03.z7
              + (141 / 16 + 50 * (k : ℝ) ^ 2 + 16 * (k : ℝ) ^ 4) * ((π ^ 2 / 6) * Z02.z5)
              + (1071 / 160 + 181 / 5 * (k : ℝ) ^ 2 + 56 / 5 * (k : ℝ) ^ 4)
                  * ((π ^ 2 / 6) ^ 2 * Y04.z3)) * Y04.tau p ^ 7
          + (-(2891923 / 12288 + 495965 / 256 * (k : ℝ) ^ 2 + 60889 / 48 * (k : ℝ) ^ 4
                + 127 * (k : ℝ) ^ 6) * (π ^ 8 / 9450)
              + (141 / 16 + 50 * (k : ℝ) ^ 2 + 16 * (k : ℝ) ^ 4) * (Y04.z3 * Z02.z5)
              + (169 / 32 + 27 * (k : ℝ) ^ 2 + 8 * (k : ℝ) ^ 4)
                  * ((π ^ 2 / 6) * Y04.z3 ^ 2)) * Y04.tau p ^ 8
          + (-(22805 / 192 + 3907 / 4 * (k : ℝ) ^ 2 + 1916 / 3 * (k : ℝ) ^ 4
                + 64 * (k : ℝ) ^ 6) * X05.z 9
              - (8343 / 64 + 4163 / 4 * (k : ℝ) ^ 2 + 660 * (k : ℝ) ^ 4
                + 64 * (k : ℝ) ^ 6) * ((π ^ 2 / 6) * Q03.z7)
              - (60093 / 640 + 29741 / 40 * (k : ℝ) ^ 2 + 2334 / 5 * (k : ℝ) ^ 4
                + 224 / 5 * (k : ℝ) ^ 6) * ((π ^ 2 / 6) ^ 2 * Z02.z5)
              - (263713 / 4480 + 130973 / 280 * (k : ℝ) ^ 2 + 2062 / 7 * (k : ℝ) ^ 4
                + 992 / 35 * (k : ℝ) ^ 6) * ((π ^ 2 / 6) ^ 3 * Y04.z3)
              + (169 / 96 + 9 * (k : ℝ) ^ 2 + 8 / 3 * (k : ℝ) ^ 4) * Y04.z3 ^ 3)
              * Y04.tau p ^ 9)|
      ≤ C * Y04.tau p ^ 10

/-- `prop:theta`, first sentence, for all five orders: for `5 ≤ n ≤ 9`, every coefficient of
`Θ_n` (the polynomial of `prop:degree`, `R05.Theta n`) is a rational combination of products of
single zeta values of total weight `n` (`Q01.IsZetaPoly n`). -/
def stmt_11 : Prop :=
  ∀ n : ℕ, 5 ≤ n → n ≤ 9 → ∀ i : ℕ, Q01.IsZetaPoly n ((R05.Theta n).coeff i)

end S01
-- END STATEMENTS

open Real MeasureTheory Set Polynomial

namespace S01

lemma zn9 : Y04.zn 9 = X05.z 9 := rfl

/-! ### The Taylor coefficient `g_{k,7}` and the values `e_{r,m}`, `m ≤ 7` -/

lemma rho_of (m k : ℕ) (c : ℝ)
    (h : iteratedDeriv m (Y04.G k) 0 / (m.factorial : ℝ) = (-1 : ℝ) ^ k * Real.sqrt 2 / 4 * c) :
    (R05.rho m).eval ((k : ℝ) ^ 2) = c := by
  have h2 := R05.proof_5 k m
  unfold R05.g at h2
  rw [h] at h2
  have hne : (-1 : ℝ) ^ k * Real.sqrt 2 / 4 ≠ 0 := by
    have : (-1 : ℝ) ^ k ≠ 0 := pow_ne_zero _ (by norm_num)
    have : Real.sqrt 2 ≠ 0 := by positivity
    positivity
  exact (mul_left_cancel₀ hne h2).symm

lemma rho_eq (m : ℕ) (P : ℝ[X])
    (hP : ∀ k : ℕ, iteratedDeriv m (Y04.G k) 0 / (m.factorial : ℝ) =
      (-1 : ℝ) ^ k * Real.sqrt 2 / 4 * P.eval ((k : ℝ) ^ 2)) : R05.rho m = P :=
  R05.poly_eq_of_sq _ _ (fun k => rho_of m k _ (hP k))

lemma e_of (m i : ℕ) (P : ℝ[X]) (h : R05.rho m = P) (hi : i < m / 2 + 1) :
    (-4 : ℝ) ^ i / ((2 * i).factorial : ℝ) * R05.e i m = P.coeff i := by
  rw [← h, R05.rho_coeff, R05.ite_t hi]

lemma rho2 : R05.rho 2 = C (-(5 / 8)) + C (-2) * X := by
  refine rho_eq 2 _ (fun k => ?_)
  rw [show ((2 : ℕ).factorial : ℝ) = 2 by norm_num [Nat.factorial], Y04.G_d2]
  simp; ring

lemma rho3 : R05.rho 3 = C (-(17 / 48)) + C (-1) * X := by
  refine rho_eq 3 _ (fun k => ?_)
  rw [show ((3 : ℕ).factorial : ℝ) = 6 by norm_num [Nat.factorial], Z02.G_d3]
  simp; ring

lemma rho4 : R05.rho 4 = C (121 / 384) + C (23 / 12) * X + C (2 / 3) * X ^ 2 := by
  refine rho_eq 4 _ (fun k => ?_)
  rw [show ((4 : ℕ).factorial : ℝ) = 24 by norm_num [Nat.factorial], V02.G_d4]
  simp; ring

lemma rho5 : R05.rho 5 = C (721 / 3840) + C (25 / 24) * X + C (1 / 3) * X ^ 2 := by
  refine rho_eq 5 _ (fun k => ?_)
  rw [show ((5 : ℕ).factorial : ℝ) = 120 by norm_num [Nat.factorial], Q03.G_d5]
  simp; ring

lemma rho6 : R05.rho 6 = C (-(1369 / 9216)) + C (-(731 / 576)) * X + C (-(31 / 36)) * X ^ 2
    + C (-(4 / 45)) * X ^ 3 := by
  refine rho_eq 6 _ (fun k => ?_)
  rw [show ((6 : ℕ).factorial : ℝ) = 720 by norm_num [Nat.factorial], R02.G_d6]
  simp; ring

set_option maxRecDepth 100000 in
set_option maxHeartbeats 400000000 in
lemma G_d7 (k : ℕ) : iteratedDeriv 7 (Y04.G k) 0 / 5040 =
    (-1 : ℝ) ^ k * Real.sqrt 2 / 4 *
      (-(58337 / 645120) - 4163 / 5760 * (k : ℝ) ^ 2 - 11 / 24 * (k : ℝ) ^ 4
        - 2 / 45 * (k : ℝ) ^ 6) := by
  rw [V02.iter_G]
  simp (config := {maxSteps := 100000000}) only [Y04.thetaX_zero, Function.iterate_succ_apply', Function.iterate_zero_apply, V02.step, V02.dterm,
    V02.shift, List.flatMap_cons, List.flatMap_nil, List.map_cons, List.map_nil, List.append_nil,
    List.cons_append, List.nil_append, ite_true, ite_false, Bool.false_eq_true]
  simp (config := {maxSteps := 100000000}) only [V02.ev, V02.tm, List.map_cons, List.map_nil, List.sum_cons, List.sum_nil, ite_true, ite_false,
    Bool.false_eq_true, Y04.cos_k, Y04.sin_k, Real.cos_pi_div_four, Real.sin_pi_div_four]
  norm_num (config := {maxSteps := 100000000})
  have : Real.sqrt 2 ^ 2 = 2 := Real.sq_sqrt (by norm_num)
  have h17 : Real.sqrt 2 ^ 17 = 256 * Real.sqrt 2 := by
    rw [show Real.sqrt 2 ^ 17 = (Real.sqrt 2 ^ 2) ^ 8 * Real.sqrt 2 by ring, this]; ring
  ring_nf
  rw [h17]; ring

lemma rho7 : R05.rho 7 = Polynomial.C (-(58337 / 645120)) + Polynomial.C (-(4163 / 5760)) * X
      + Polynomial.C (-(11 / 24)) * X ^ 2 + Polynomial.C (-(2 / 45)) * X ^ 3 := by
  refine rho_eq 7 _ (fun k => ?_)
  rw [show ((7 : ℕ).factorial : ℝ) = 5040 by norm_num [Nat.factorial], G_d7]
  simp; ring

lemma evals1 : R05.e 0 2 = -(5 / 8) ∧ R05.e 0 3 = -(17 / 48) ∧ R05.e 0 4 = 121 / 384 ∧
    R05.e 1 4 = -(23 / 24) ∧ R05.e 1 5 = -(25 / 48) ∧ R05.e 1 6 = 731 / 1152 ∧
    R05.e 2 6 = -(31 / 24) ∧ R05.e 2 7 = -(11 / 16) := by
  have a1 := e_of 2 0 _ rho2 (by norm_num)
  have a2 := e_of 3 0 _ rho3 (by norm_num)
  have a3 := e_of 4 0 _ rho4 (by norm_num)
  have a4 := e_of 4 1 _ rho4 (by norm_num)
  have a5 := e_of 5 1 _ rho5 (by norm_num)
  have a6 := e_of 6 1 _ rho6 (by norm_num)
  have a7 := e_of 6 2 _ rho6 (by norm_num)
  have a8 := e_of 7 2 _ rho7 (by norm_num)
  simp [Polynomial.coeff_X, Polynomial.coeff_X_pow, Polynomial.coeff_C, Nat.factorial] at a1 a2 a3 a4 a5 a6 a7 a8
  refine ⟨by linarith, by linarith, by linarith, by linarith, by linarith, by linarith, by linarith,
    by linarith⟩

lemma evals2 : R05.e 0 5 = 721 / 3840 ∧ R05.e 0 6 = -(1369 / 9216) ∧
    R05.e 0 7 = -(58337 / 645120) ∧ R05.e 1 7 = 4163 / 11520 := by
  have a1 := e_of 5 0 _ rho5 (by norm_num)
  have a2 := e_of 6 0 _ rho6 (by norm_num)
  have a3 := e_of 7 0 _ rho7 (by norm_num)
  have a4 := e_of 7 1 _ rho7 (by norm_num)
  simp [Polynomial.coeff_X, Polynomial.coeff_X_pow, Polynomial.coeff_C, Nat.factorial] at a1 a2 a3 a4
  refine ⟨by linarith, by linarith, by linarith, by linarith⟩

lemma G_d7' (k : ℕ) : iteratedDeriv 7 (Y04.G k) 0 =
    5040 * ((-1 : ℝ) ^ k * Real.sqrt 2 / 4 *
      (-(58337 / 645120) - 4163 / 5760 * (k : ℝ) ^ 2 - 11 / 24 * (k : ℝ) ^ 4
        - 2 / 45 * (k : ℝ) ^ 6)) := by
  have := G_d7 k; linarith


/-! ### The order-8 Beta-identity machinery (one order beyond `R02`) -/

lemma em8_bound {q : ℝ} (hq : 0 ≤ q) :
    |Real.exp (-q) - (1 - q + q ^ 2 / 2 - q ^ 3 / 6 + q ^ 4 / 24 - q ^ 5 / 120 + q ^ 6 / 720 - q ^ 7 / 5040)| ≤ 4 * q ^ 8 := by
  rcases le_or_gt q 1 with h | h
  · have h1 : |-q| ≤ 1 := by rw [abs_neg, abs_of_nonneg hq]; exact h
    have he := Real.exp_bound h1 (n := 8) (by norm_num)
    simp only [Finset.sum_range_succ, Finset.sum_range_zero, Nat.factorial] at he
    norm_num at he
    rw [abs_of_nonneg hq] at he
    have e : Real.exp (-q) - (1 - q + q ^ 2 / 2 - q ^ 3 / 6 + q ^ 4 / 24 - q ^ 5 / 120 + q ^ 6 / 720 - q ^ 7 / 5040)
        = Real.exp (-q) - (1 + -q + q ^ 2 / 2 + (-q) ^ 3 / 6 + q ^ 4 / 24 + (-q) ^ 5 / 120 + q ^ 6 / 720 + (-q) ^ 7 / 5040) := by ring
    rw [e]
    have : 0 ≤ q ^ 8 := by positivity
    nlinarith
  · have he : Real.exp (-q) ≤ 1 := Real.exp_le_one_iff.2 (by linarith)
    have he0 := Real.exp_pos (-q)
    have h2 : q ≤ q ^ 2 := by nlinarith
    have h3 : q ^ 2 ≤ q ^ 3 := by nlinarith
    have h4 : q ^ 3 ≤ q ^ 4 := by nlinarith
    have h5 : q ^ 4 ≤ q ^ 5 := by nlinarith
    have h6 : q ^ 5 ≤ q ^ 6 := by nlinarith
    have h7 : q ^ 6 ≤ q ^ 7 := by nlinarith
    have h8 : q ^ 7 ≤ q ^ 8 := by nlinarith
    rw [abs_le]
    constructor <;> nlinarith

lemma ep9_bound {p : ℝ} (hp0 : 0 ≤ p) (hp1 : p ≤ 1) :
    |Real.exp (-p) - 1 - (-p + p ^ 2 / 2 - p ^ 3 / 6 + p ^ 4 / 24 - p ^ 5 / 120 + p ^ 6 / 720 - p ^ 7 / 5040 + p ^ 8 / 40320)| ≤ p ^ 9 := by
  have h1 : |-p| ≤ 1 := by rw [abs_neg, abs_of_nonneg hp0]; exact hp1
  have he := Real.exp_bound h1 (n := 9) (by norm_num)
  simp only [Finset.sum_range_succ, Finset.sum_range_zero, Nat.factorial] at he
  norm_num at he
  rw [abs_of_nonneg hp0] at he
  have e : Real.exp (-p) - 1 - (-p + p ^ 2 / 2 - p ^ 3 / 6 + p ^ 4 / 24 - p ^ 5 / 120 + p ^ 6 / 720 - p ^ 7 / 5040 + p ^ 8 / 40320)
      = Real.exp (-p) - (1 + -p + p ^ 2 / 2 + (-p) ^ 3 / 6 + p ^ 4 / 24 + (-p) ^ 5 / 120 + p ^ 6 / 720 + (-p) ^ 7 / 5040 + p ^ 8 / 40320) := by ring
  rw [e]
  have : 0 ≤ p ^ 9 := by positivity
  nlinarith

/-- the dominating function of the order-8 remainder -/
noncomputable def B10 (s : ℝ) : ℝ := (1 + s + s ^ 2 + s ^ 3 + s ^ 4 + s ^ 5 + s ^ 6 + s ^ 7 + s ^ 8) * Real.exp (-s)

lemma B10_int : IntegrableOn B10 (Set.Ioi 0) := by
  have h := (((((((((Y04.intOn_pow_exp1 0).add (Y04.intOn_pow_exp1 1)).add (Y04.intOn_pow_exp1 2)).add (Y04.intOn_pow_exp1 3)).add (Y04.intOn_pow_exp1 4)).add (Y04.intOn_pow_exp1 5)).add (Y04.intOn_pow_exp1 6)).add (Y04.intOn_pow_exp1 7)).add (Y04.intOn_pow_exp1 8))
  refine IntegrableOn.congr_fun h (fun s _ => ?_) measurableSet_Ioi
  simp only [B10, Pi.add_apply]; ring

lemma B10_nonneg {s : ℝ} (hs : 0 ≤ s) : 0 ≤ B10 s := by unfold B10; positivity

lemma sj_le_B10 {s : ℝ} (hs : 0 ≤ s) {j : ℕ} (hj : j ≤ 8) : s ^ j * Real.exp (-s) ≤ B10 s := by
  unfold B10
  have h0 : 0 ≤ s ^ 0 * Real.exp (-s) := by positivity
  have h1 : 0 ≤ s ^ 1 * Real.exp (-s) := by positivity
  have h2 : 0 ≤ s ^ 2 * Real.exp (-s) := by positivity
  have h3 : 0 ≤ s ^ 3 * Real.exp (-s) := by positivity
  have h4 : 0 ≤ s ^ 4 * Real.exp (-s) := by positivity
  have h5 : 0 ≤ s ^ 5 * Real.exp (-s) := by positivity
  have h6 : 0 ≤ s ^ 6 * Real.exp (-s) := by positivity
  have h7 : 0 ≤ s ^ 7 * Real.exp (-s) := by positivity
  have h8 : 0 ≤ s ^ 8 * Real.exp (-s) := by positivity
  have w : (1 + s + s ^ 2 + s ^ 3 + s ^ 4 + s ^ 5 + s ^ 6 + s ^ 7 + s ^ 8) * Real.exp (-s)
    = s ^ 0 * Real.exp (-s) + s ^ 1 * Real.exp (-s) + s ^ 2 * Real.exp (-s) + s ^ 3 * Real.exp (-s) + s ^ 4 * Real.exp (-s) + s ^ 5 * Real.exp (-s) + s ^ 6 * Real.exp (-s) + s ^ 7 * Real.exp (-s) + s ^ 8 * Real.exp (-s) := by ring
  rw [w]
  interval_cases j <;> linarith

lemma mono_bd9 {u v t L s : ℝ} (hu0 : 0 ≤ u) (hut : u ≤ t) (hv0 : 0 ≤ v) (hvt : v ≤ t) (ht1 : t ≤ 1)
    (hL0 : 0 ≤ L) (hL1 : L ≤ 1) (hLe : L ≤ Real.exp (-s)) (hs : 0 ≤ s) {i j : ℕ} (hi : 1 ≤ i)
    (hij : 9 ≤ i + j) :
    (u * L) ^ i * (v * s) ^ j ≤ t ^ 9 * (s ^ j * Real.exp (-s)) := by
  have ht0 : 0 ≤ t := hu0.trans hut
  rw [mul_pow, mul_pow]
  have a1 : u ^ i ≤ t ^ i := pow_le_pow_left₀ hu0 hut i
  have a2 : v ^ j ≤ t ^ j := pow_le_pow_left₀ hv0 hvt j
  have a3 : t ^ i * t ^ j ≤ t ^ 9 := by
    rw [← pow_add]; exact pow_le_pow_of_le_one ht0 ht1 hij
  have a4 : L ^ i ≤ Real.exp (-s) := (pow_le_of_le_one hL0 hL1 (by omega)).trans hLe
  have b1 : u ^ i * v ^ j ≤ t ^ 9 :=
    (mul_le_mul a1 a2 (by positivity) (by positivity)).trans a3
  calc u ^ i * L ^ i * (v ^ j * s ^ j) = (u ^ i * v ^ j) * L ^ i * s ^ j := by ring
    _ ≤ t ^ 9 * Real.exp (-s) * s ^ j := by
        apply mul_le_mul_of_nonneg_right _ (by positivity)
        exact mul_le_mul b1 a4 (by positivity) (by positivity)
    _ = _ := by ring

lemma hmono9 {u v t L s p q c : ℝ} (hu0 : 0 ≤ u) (hut : u ≤ t) (hv0 : 0 ≤ v) (hvt : v ≤ t) (ht1 : t ≤ 1)
    (hL0 : 0 ≤ L) (hL1 : L ≤ 1) (hLe : L ≤ Real.exp (-s)) (hs : 0 ≤ s) (hp : p = u * L) (hq : q = v * s)
    (hc : |c| ≤ 1) {i j : ℕ} (hi : 1 ≤ i) (hij : 9 ≤ i + j) (hj : j ≤ 8) :
    |c * (p ^ i * q ^ j)| ≤ t ^ 9 * B10 s := by
  subst hp hq
  have h0 : 0 ≤ (u * L) ^ i * (v * s) ^ j := by positivity
  have h1 := mono_bd9 hu0 hut hv0 hvt ht1 hL0 hL1 hLe hs hi hij
  have h2 := sj_le_B10 hs hj
  have ht0 : 0 ≤ t := hu0.trans hut
  rw [abs_mul, abs_of_nonneg h0]
  calc |c| * ((u * L) ^ i * (v * s) ^ j) ≤ 1 * (t ^ 9 * (s ^ j * Real.exp (-s))) :=
        mul_le_mul hc h1 h0 zero_le_one
    _ ≤ t ^ 9 * B10 s := by rw [one_mul]; exact mul_le_mul_of_nonneg_left h2 (by positivity)

/-- the degree-`≤ 8` part of `(e^{-uℓ} - 1) e^{-vs}` -/
noncomputable def P8 (u v s : ℝ) : ℝ :=
  (-1) * (u ^ 1 * v ^ 0) * (s ^ 0 * Y04.ell s ^ 1)
    + (1) * (u ^ 1 * v ^ 1) * (s ^ 1 * Y04.ell s ^ 1)
    + (-1 / 2) * (u ^ 1 * v ^ 2) * (s ^ 2 * Y04.ell s ^ 1)
    + (1 / 6) * (u ^ 1 * v ^ 3) * (s ^ 3 * Y04.ell s ^ 1)
    + (-1 / 24) * (u ^ 1 * v ^ 4) * (s ^ 4 * Y04.ell s ^ 1)
    + (1 / 120) * (u ^ 1 * v ^ 5) * (s ^ 5 * Y04.ell s ^ 1)
    + (-1 / 720) * (u ^ 1 * v ^ 6) * (s ^ 6 * Y04.ell s ^ 1)
    + (1 / 5040) * (u ^ 1 * v ^ 7) * (s ^ 7 * Y04.ell s ^ 1)
    + (1 / 2) * (u ^ 2 * v ^ 0) * (s ^ 0 * Y04.ell s ^ 2)
    + (-1 / 2) * (u ^ 2 * v ^ 1) * (s ^ 1 * Y04.ell s ^ 2)
    + (1 / 4) * (u ^ 2 * v ^ 2) * (s ^ 2 * Y04.ell s ^ 2)
    + (-1 / 12) * (u ^ 2 * v ^ 3) * (s ^ 3 * Y04.ell s ^ 2)
    + (1 / 48) * (u ^ 2 * v ^ 4) * (s ^ 4 * Y04.ell s ^ 2)
    + (-1 / 240) * (u ^ 2 * v ^ 5) * (s ^ 5 * Y04.ell s ^ 2)
    + (1 / 1440) * (u ^ 2 * v ^ 6) * (s ^ 6 * Y04.ell s ^ 2)
    + (-1 / 6) * (u ^ 3 * v ^ 0) * (s ^ 0 * Y04.ell s ^ 3)
    + (1 / 6) * (u ^ 3 * v ^ 1) * (s ^ 1 * Y04.ell s ^ 3)
    + (-1 / 12) * (u ^ 3 * v ^ 2) * (s ^ 2 * Y04.ell s ^ 3)
    + (1 / 36) * (u ^ 3 * v ^ 3) * (s ^ 3 * Y04.ell s ^ 3)
    + (-1 / 144) * (u ^ 3 * v ^ 4) * (s ^ 4 * Y04.ell s ^ 3)
    + (1 / 720) * (u ^ 3 * v ^ 5) * (s ^ 5 * Y04.ell s ^ 3)
    + (1 / 24) * (u ^ 4 * v ^ 0) * (s ^ 0 * Y04.ell s ^ 4)
    + (-1 / 24) * (u ^ 4 * v ^ 1) * (s ^ 1 * Y04.ell s ^ 4)
    + (1 / 48) * (u ^ 4 * v ^ 2) * (s ^ 2 * Y04.ell s ^ 4)
    + (-1 / 144) * (u ^ 4 * v ^ 3) * (s ^ 3 * Y04.ell s ^ 4)
    + (1 / 576) * (u ^ 4 * v ^ 4) * (s ^ 4 * Y04.ell s ^ 4)
    + (-1 / 120) * (u ^ 5 * v ^ 0) * (s ^ 0 * Y04.ell s ^ 5)
    + (1 / 120) * (u ^ 5 * v ^ 1) * (s ^ 1 * Y04.ell s ^ 5)
    + (-1 / 240) * (u ^ 5 * v ^ 2) * (s ^ 2 * Y04.ell s ^ 5)
    + (1 / 720) * (u ^ 5 * v ^ 3) * (s ^ 3 * Y04.ell s ^ 5)
    + (1 / 720) * (u ^ 6 * v ^ 0) * (s ^ 0 * Y04.ell s ^ 6)
    + (-1 / 720) * (u ^ 6 * v ^ 1) * (s ^ 1 * Y04.ell s ^ 6)
    + (1 / 1440) * (u ^ 6 * v ^ 2) * (s ^ 2 * Y04.ell s ^ 6)
    + (-1 / 5040) * (u ^ 7 * v ^ 0) * (s ^ 0 * Y04.ell s ^ 7)
    + (1 / 5040) * (u ^ 7 * v ^ 1) * (s ^ 1 * Y04.ell s ^ 7)
    + (1 / 40320) * (u ^ 8 * v ^ 0) * (s ^ 0 * Y04.ell s ^ 8)

set_option maxHeartbeats 64000000 in
lemma pt8 {u v t s : ℝ} (hu0 : 0 ≤ u) (hut : u ≤ t) (hv0 : 0 ≤ v) (hvt : v ≤ t) (ht1 : t ≤ 1)
    (hs : 0 ≤ s) :
    |(Real.exp (-u * Y04.ell s) - 1) * Real.exp (-v * s) - P8 u v s| ≤ 37 * t ^ 9 * B10 s := by
  have hl0 := Y04.ell_nonneg s
  have hl1 := Y04.ell_le_one s hs
  have hle := Y04.ell_le s
  have ht0 : 0 ≤ t := hu0.trans hut
  unfold P8
  set L := Y04.ell s
  set p := u * L with hp
  set q := v * s with hq
  have hp0 : 0 ≤ p := by positivity
  have hp1 : p ≤ 1 := by
    have := mul_le_mul (hut.trans ht1) hl1 hl0 zero_le_one; rw [hp]; linarith
  have hq0 : 0 ≤ q := by positivity
  set A := -p + p ^ 2 / 2 - p ^ 3 / 6 + p ^ 4 / 24 - p ^ 5 / 120 + p ^ 6 / 720 - p ^ 7 / 5040 + p ^ 8 / 40320 with hA
  set r1 := Real.exp (-p) - 1 - A with hr1
  set B := 1 - q + q ^ 2 / 2 - q ^ 3 / 6 + q ^ 4 / 24 - q ^ 5 / 120 + q ^ 6 / 720 - q ^ 7 / 5040 with hB
  set r2 := Real.exp (-q) - B with hr2
  have hr1b : |r1| ≤ p ^ 9 := ep9_bound hp0 hp1
  have hr2b : |r2| ≤ 4 * q ^ 8 := em8_bound hq0
  set HT := (-1 / 10080) * (p ^ 2 * q ^ 7) + (-1 / 4320) * (p ^ 3 * q ^ 6) + (1 / 30240) * (p ^ 3 * q ^ 7) + (-1 / 2880) * (p ^ 4 * q ^ 5) + (1 / 17280) * (p ^ 4 * q ^ 6) + (-1 / 120960) * (p ^ 4 * q ^ 7) + (-1 / 2880) * (p ^ 5 * q ^ 4) + (1 / 14400) * (p ^ 5 * q ^ 5) + (-1 / 86400) * (p ^ 5 * q ^ 6) + (1 / 604800) * (p ^ 5 * q ^ 7) + (-1 / 4320) * (p ^ 6 * q ^ 3) + (1 / 17280) * (p ^ 6 * q ^ 4) + (-1 / 86400) * (p ^ 6 * q ^ 5) + (1 / 518400) * (p ^ 6 * q ^ 6) + (-1 / 3628800) * (p ^ 6 * q ^ 7) + (-1 / 10080) * (p ^ 7 * q ^ 2) + (1 / 30240) * (p ^ 7 * q ^ 3) + (-1 / 120960) * (p ^ 7 * q ^ 4) + (1 / 604800) * (p ^ 7 * q ^ 5) + (-1 / 3628800) * (p ^ 7 * q ^ 6) + (1 / 25401600) * (p ^ 7 * q ^ 7) + (-1 / 40320) * (p ^ 8 * q ^ 1) + (1 / 80640) * (p ^ 8 * q ^ 2) + (-1 / 241920) * (p ^ 8 * q ^ 3) + (1 / 967680) * (p ^ 8 * q ^ 4) + (-1 / 4838400) * (p ^ 8 * q ^ 5) + (1 / 29030400) * (p ^ 8 * q ^ 6) + (-1 / 203212800) * (p ^ 8 * q ^ 7) with hHT
  have e1 : -u * L = -p := by rw [hp]; ring
  have e2 : -v * s = -q := by rw [hq]; ring
  rw [e1, e2]
  have key : (Real.exp (-p) - 1) * Real.exp (-q) - ((-1) * (u ^ 1 * v ^ 0) * (s ^ 0 * L ^ 1)
    + (1) * (u ^ 1 * v ^ 1) * (s ^ 1 * L ^ 1)
    + (-1 / 2) * (u ^ 1 * v ^ 2) * (s ^ 2 * L ^ 1)
    + (1 / 6) * (u ^ 1 * v ^ 3) * (s ^ 3 * L ^ 1)
    + (-1 / 24) * (u ^ 1 * v ^ 4) * (s ^ 4 * L ^ 1)
    + (1 / 120) * (u ^ 1 * v ^ 5) * (s ^ 5 * L ^ 1)
    + (-1 / 720) * (u ^ 1 * v ^ 6) * (s ^ 6 * L ^ 1)
    + (1 / 5040) * (u ^ 1 * v ^ 7) * (s ^ 7 * L ^ 1)
    + (1 / 2) * (u ^ 2 * v ^ 0) * (s ^ 0 * L ^ 2)
    + (-1 / 2) * (u ^ 2 * v ^ 1) * (s ^ 1 * L ^ 2)
    + (1 / 4) * (u ^ 2 * v ^ 2) * (s ^ 2 * L ^ 2)
    + (-1 / 12) * (u ^ 2 * v ^ 3) * (s ^ 3 * L ^ 2)
    + (1 / 48) * (u ^ 2 * v ^ 4) * (s ^ 4 * L ^ 2)
    + (-1 / 240) * (u ^ 2 * v ^ 5) * (s ^ 5 * L ^ 2)
    + (1 / 1440) * (u ^ 2 * v ^ 6) * (s ^ 6 * L ^ 2)
    + (-1 / 6) * (u ^ 3 * v ^ 0) * (s ^ 0 * L ^ 3)
    + (1 / 6) * (u ^ 3 * v ^ 1) * (s ^ 1 * L ^ 3)
    + (-1 / 12) * (u ^ 3 * v ^ 2) * (s ^ 2 * L ^ 3)
    + (1 / 36) * (u ^ 3 * v ^ 3) * (s ^ 3 * L ^ 3)
    + (-1 / 144) * (u ^ 3 * v ^ 4) * (s ^ 4 * L ^ 3)
    + (1 / 720) * (u ^ 3 * v ^ 5) * (s ^ 5 * L ^ 3)
    + (1 / 24) * (u ^ 4 * v ^ 0) * (s ^ 0 * L ^ 4)
    + (-1 / 24) * (u ^ 4 * v ^ 1) * (s ^ 1 * L ^ 4)
    + (1 / 48) * (u ^ 4 * v ^ 2) * (s ^ 2 * L ^ 4)
    + (-1 / 144) * (u ^ 4 * v ^ 3) * (s ^ 3 * L ^ 4)
    + (1 / 576) * (u ^ 4 * v ^ 4) * (s ^ 4 * L ^ 4)
    + (-1 / 120) * (u ^ 5 * v ^ 0) * (s ^ 0 * L ^ 5)
    + (1 / 120) * (u ^ 5 * v ^ 1) * (s ^ 1 * L ^ 5)
    + (-1 / 240) * (u ^ 5 * v ^ 2) * (s ^ 2 * L ^ 5)
    + (1 / 720) * (u ^ 5 * v ^ 3) * (s ^ 3 * L ^ 5)
    + (1 / 720) * (u ^ 6 * v ^ 0) * (s ^ 0 * L ^ 6)
    + (-1 / 720) * (u ^ 6 * v ^ 1) * (s ^ 1 * L ^ 6)
    + (1 / 1440) * (u ^ 6 * v ^ 2) * (s ^ 2 * L ^ 6)
    + (-1 / 5040) * (u ^ 7 * v ^ 0) * (s ^ 0 * L ^ 7)
    + (1 / 5040) * (u ^ 7 * v ^ 1) * (s ^ 1 * L ^ 7)
    + (1 / 40320) * (u ^ 8 * v ^ 0) * (s ^ 0 * L ^ 8))
      = HT + A * r2 + r1 * Real.exp (-q) := by
    rw [hHT, hr1, hr2, hA, hB, hp, hq]; ring
  rw [key]
  have hB10 : 0 ≤ B10 s := B10_nonneg hs
  have m0 := hmono9 hu0 hut hv0 hvt ht1 hl0 hl1 hle hs hp hq (c := (-1 / 10080)) (by rw [abs_le]; constructor <;> norm_num) (i := 2) (j := 7) (by norm_num) (by norm_num) (by norm_num)
  have m1 := hmono9 hu0 hut hv0 hvt ht1 hl0 hl1 hle hs hp hq (c := (-1 / 4320)) (by rw [abs_le]; constructor <;> norm_num) (i := 3) (j := 6) (by norm_num) (by norm_num) (by norm_num)
  have m2 := hmono9 hu0 hut hv0 hvt ht1 hl0 hl1 hle hs hp hq (c := (1 / 30240)) (by rw [abs_le]; constructor <;> norm_num) (i := 3) (j := 7) (by norm_num) (by norm_num) (by norm_num)
  have m3 := hmono9 hu0 hut hv0 hvt ht1 hl0 hl1 hle hs hp hq (c := (-1 / 2880)) (by rw [abs_le]; constructor <;> norm_num) (i := 4) (j := 5) (by norm_num) (by norm_num) (by norm_num)
  have m4 := hmono9 hu0 hut hv0 hvt ht1 hl0 hl1 hle hs hp hq (c := (1 / 17280)) (by rw [abs_le]; constructor <;> norm_num) (i := 4) (j := 6) (by norm_num) (by norm_num) (by norm_num)
  have m5 := hmono9 hu0 hut hv0 hvt ht1 hl0 hl1 hle hs hp hq (c := (-1 / 120960)) (by rw [abs_le]; constructor <;> norm_num) (i := 4) (j := 7) (by norm_num) (by norm_num) (by norm_num)
  have m6 := hmono9 hu0 hut hv0 hvt ht1 hl0 hl1 hle hs hp hq (c := (-1 / 2880)) (by rw [abs_le]; constructor <;> norm_num) (i := 5) (j := 4) (by norm_num) (by norm_num) (by norm_num)
  have m7 := hmono9 hu0 hut hv0 hvt ht1 hl0 hl1 hle hs hp hq (c := (1 / 14400)) (by rw [abs_le]; constructor <;> norm_num) (i := 5) (j := 5) (by norm_num) (by norm_num) (by norm_num)
  have m8 := hmono9 hu0 hut hv0 hvt ht1 hl0 hl1 hle hs hp hq (c := (-1 / 86400)) (by rw [abs_le]; constructor <;> norm_num) (i := 5) (j := 6) (by norm_num) (by norm_num) (by norm_num)
  have m9 := hmono9 hu0 hut hv0 hvt ht1 hl0 hl1 hle hs hp hq (c := (1 / 604800)) (by rw [abs_le]; constructor <;> norm_num) (i := 5) (j := 7) (by norm_num) (by norm_num) (by norm_num)
  have m10 := hmono9 hu0 hut hv0 hvt ht1 hl0 hl1 hle hs hp hq (c := (-1 / 4320)) (by rw [abs_le]; constructor <;> norm_num) (i := 6) (j := 3) (by norm_num) (by norm_num) (by norm_num)
  have m11 := hmono9 hu0 hut hv0 hvt ht1 hl0 hl1 hle hs hp hq (c := (1 / 17280)) (by rw [abs_le]; constructor <;> norm_num) (i := 6) (j := 4) (by norm_num) (by norm_num) (by norm_num)
  have m12 := hmono9 hu0 hut hv0 hvt ht1 hl0 hl1 hle hs hp hq (c := (-1 / 86400)) (by rw [abs_le]; constructor <;> norm_num) (i := 6) (j := 5) (by norm_num) (by norm_num) (by norm_num)
  have m13 := hmono9 hu0 hut hv0 hvt ht1 hl0 hl1 hle hs hp hq (c := (1 / 518400)) (by rw [abs_le]; constructor <;> norm_num) (i := 6) (j := 6) (by norm_num) (by norm_num) (by norm_num)
  have m14 := hmono9 hu0 hut hv0 hvt ht1 hl0 hl1 hle hs hp hq (c := (-1 / 3628800)) (by rw [abs_le]; constructor <;> norm_num) (i := 6) (j := 7) (by norm_num) (by norm_num) (by norm_num)
  have m15 := hmono9 hu0 hut hv0 hvt ht1 hl0 hl1 hle hs hp hq (c := (-1 / 10080)) (by rw [abs_le]; constructor <;> norm_num) (i := 7) (j := 2) (by norm_num) (by norm_num) (by norm_num)
  have m16 := hmono9 hu0 hut hv0 hvt ht1 hl0 hl1 hle hs hp hq (c := (1 / 30240)) (by rw [abs_le]; constructor <;> norm_num) (i := 7) (j := 3) (by norm_num) (by norm_num) (by norm_num)
  have m17 := hmono9 hu0 hut hv0 hvt ht1 hl0 hl1 hle hs hp hq (c := (-1 / 120960)) (by rw [abs_le]; constructor <;> norm_num) (i := 7) (j := 4) (by norm_num) (by norm_num) (by norm_num)
  have m18 := hmono9 hu0 hut hv0 hvt ht1 hl0 hl1 hle hs hp hq (c := (1 / 604800)) (by rw [abs_le]; constructor <;> norm_num) (i := 7) (j := 5) (by norm_num) (by norm_num) (by norm_num)
  have m19 := hmono9 hu0 hut hv0 hvt ht1 hl0 hl1 hle hs hp hq (c := (-1 / 3628800)) (by rw [abs_le]; constructor <;> norm_num) (i := 7) (j := 6) (by norm_num) (by norm_num) (by norm_num)
  have m20 := hmono9 hu0 hut hv0 hvt ht1 hl0 hl1 hle hs hp hq (c := (1 / 25401600)) (by rw [abs_le]; constructor <;> norm_num) (i := 7) (j := 7) (by norm_num) (by norm_num) (by norm_num)
  have m21 := hmono9 hu0 hut hv0 hvt ht1 hl0 hl1 hle hs hp hq (c := (-1 / 40320)) (by rw [abs_le]; constructor <;> norm_num) (i := 8) (j := 1) (by norm_num) (by norm_num) (by norm_num)
  have m22 := hmono9 hu0 hut hv0 hvt ht1 hl0 hl1 hle hs hp hq (c := (1 / 80640)) (by rw [abs_le]; constructor <;> norm_num) (i := 8) (j := 2) (by norm_num) (by norm_num) (by norm_num)
  have m23 := hmono9 hu0 hut hv0 hvt ht1 hl0 hl1 hle hs hp hq (c := (-1 / 241920)) (by rw [abs_le]; constructor <;> norm_num) (i := 8) (j := 3) (by norm_num) (by norm_num) (by norm_num)
  have m24 := hmono9 hu0 hut hv0 hvt ht1 hl0 hl1 hle hs hp hq (c := (1 / 967680)) (by rw [abs_le]; constructor <;> norm_num) (i := 8) (j := 4) (by norm_num) (by norm_num) (by norm_num)
  have m25 := hmono9 hu0 hut hv0 hvt ht1 hl0 hl1 hle hs hp hq (c := (-1 / 4838400)) (by rw [abs_le]; constructor <;> norm_num) (i := 8) (j := 5) (by norm_num) (by norm_num) (by norm_num)
  have m26 := hmono9 hu0 hut hv0 hvt ht1 hl0 hl1 hle hs hp hq (c := (1 / 29030400)) (by rw [abs_le]; constructor <;> norm_num) (i := 8) (j := 6) (by norm_num) (by norm_num) (by norm_num)
  have m27 := hmono9 hu0 hut hv0 hvt ht1 hl0 hl1 hle hs hp hq (c := (-1 / 203212800)) (by rw [abs_le]; constructor <;> norm_num) (i := 8) (j := 7) (by norm_num) (by norm_num) (by norm_num)
  have c1 := abs_add_le ((-1 / 10080) * (p ^ 2 * q ^ 7)) ((-1 / 4320) * (p ^ 3 * q ^ 6))
  have c2 := abs_add_le ((-1 / 10080) * (p ^ 2 * q ^ 7) + (-1 / 4320) * (p ^ 3 * q ^ 6)) ((1 / 30240) * (p ^ 3 * q ^ 7))
  have c3 := abs_add_le ((-1 / 10080) * (p ^ 2 * q ^ 7) + (-1 / 4320) * (p ^ 3 * q ^ 6) + (1 / 30240) * (p ^ 3 * q ^ 7)) ((-1 / 2880) * (p ^ 4 * q ^ 5))
  have c4 := abs_add_le ((-1 / 10080) * (p ^ 2 * q ^ 7) + (-1 / 4320) * (p ^ 3 * q ^ 6) + (1 / 30240) * (p ^ 3 * q ^ 7) + (-1 / 2880) * (p ^ 4 * q ^ 5)) ((1 / 17280) * (p ^ 4 * q ^ 6))
  have c5 := abs_add_le ((-1 / 10080) * (p ^ 2 * q ^ 7) + (-1 / 4320) * (p ^ 3 * q ^ 6) + (1 / 30240) * (p ^ 3 * q ^ 7) + (-1 / 2880) * (p ^ 4 * q ^ 5) + (1 / 17280) * (p ^ 4 * q ^ 6)) ((-1 / 120960) * (p ^ 4 * q ^ 7))
  have c6 := abs_add_le ((-1 / 10080) * (p ^ 2 * q ^ 7) + (-1 / 4320) * (p ^ 3 * q ^ 6) + (1 / 30240) * (p ^ 3 * q ^ 7) + (-1 / 2880) * (p ^ 4 * q ^ 5) + (1 / 17280) * (p ^ 4 * q ^ 6) + (-1 / 120960) * (p ^ 4 * q ^ 7)) ((-1 / 2880) * (p ^ 5 * q ^ 4))
  have c7 := abs_add_le ((-1 / 10080) * (p ^ 2 * q ^ 7) + (-1 / 4320) * (p ^ 3 * q ^ 6) + (1 / 30240) * (p ^ 3 * q ^ 7) + (-1 / 2880) * (p ^ 4 * q ^ 5) + (1 / 17280) * (p ^ 4 * q ^ 6) + (-1 / 120960) * (p ^ 4 * q ^ 7) + (-1 / 2880) * (p ^ 5 * q ^ 4)) ((1 / 14400) * (p ^ 5 * q ^ 5))
  have c8 := abs_add_le ((-1 / 10080) * (p ^ 2 * q ^ 7) + (-1 / 4320) * (p ^ 3 * q ^ 6) + (1 / 30240) * (p ^ 3 * q ^ 7) + (-1 / 2880) * (p ^ 4 * q ^ 5) + (1 / 17280) * (p ^ 4 * q ^ 6) + (-1 / 120960) * (p ^ 4 * q ^ 7) + (-1 / 2880) * (p ^ 5 * q ^ 4) + (1 / 14400) * (p ^ 5 * q ^ 5)) ((-1 / 86400) * (p ^ 5 * q ^ 6))
  have c9 := abs_add_le ((-1 / 10080) * (p ^ 2 * q ^ 7) + (-1 / 4320) * (p ^ 3 * q ^ 6) + (1 / 30240) * (p ^ 3 * q ^ 7) + (-1 / 2880) * (p ^ 4 * q ^ 5) + (1 / 17280) * (p ^ 4 * q ^ 6) + (-1 / 120960) * (p ^ 4 * q ^ 7) + (-1 / 2880) * (p ^ 5 * q ^ 4) + (1 / 14400) * (p ^ 5 * q ^ 5) + (-1 / 86400) * (p ^ 5 * q ^ 6)) ((1 / 604800) * (p ^ 5 * q ^ 7))
  have c10 := abs_add_le ((-1 / 10080) * (p ^ 2 * q ^ 7) + (-1 / 4320) * (p ^ 3 * q ^ 6) + (1 / 30240) * (p ^ 3 * q ^ 7) + (-1 / 2880) * (p ^ 4 * q ^ 5) + (1 / 17280) * (p ^ 4 * q ^ 6) + (-1 / 120960) * (p ^ 4 * q ^ 7) + (-1 / 2880) * (p ^ 5 * q ^ 4) + (1 / 14400) * (p ^ 5 * q ^ 5) + (-1 / 86400) * (p ^ 5 * q ^ 6) + (1 / 604800) * (p ^ 5 * q ^ 7)) ((-1 / 4320) * (p ^ 6 * q ^ 3))
  have c11 := abs_add_le ((-1 / 10080) * (p ^ 2 * q ^ 7) + (-1 / 4320) * (p ^ 3 * q ^ 6) + (1 / 30240) * (p ^ 3 * q ^ 7) + (-1 / 2880) * (p ^ 4 * q ^ 5) + (1 / 17280) * (p ^ 4 * q ^ 6) + (-1 / 120960) * (p ^ 4 * q ^ 7) + (-1 / 2880) * (p ^ 5 * q ^ 4) + (1 / 14400) * (p ^ 5 * q ^ 5) + (-1 / 86400) * (p ^ 5 * q ^ 6) + (1 / 604800) * (p ^ 5 * q ^ 7) + (-1 / 4320) * (p ^ 6 * q ^ 3)) ((1 / 17280) * (p ^ 6 * q ^ 4))
  have c12 := abs_add_le ((-1 / 10080) * (p ^ 2 * q ^ 7) + (-1 / 4320) * (p ^ 3 * q ^ 6) + (1 / 30240) * (p ^ 3 * q ^ 7) + (-1 / 2880) * (p ^ 4 * q ^ 5) + (1 / 17280) * (p ^ 4 * q ^ 6) + (-1 / 120960) * (p ^ 4 * q ^ 7) + (-1 / 2880) * (p ^ 5 * q ^ 4) + (1 / 14400) * (p ^ 5 * q ^ 5) + (-1 / 86400) * (p ^ 5 * q ^ 6) + (1 / 604800) * (p ^ 5 * q ^ 7) + (-1 / 4320) * (p ^ 6 * q ^ 3) + (1 / 17280) * (p ^ 6 * q ^ 4)) ((-1 / 86400) * (p ^ 6 * q ^ 5))
  have c13 := abs_add_le ((-1 / 10080) * (p ^ 2 * q ^ 7) + (-1 / 4320) * (p ^ 3 * q ^ 6) + (1 / 30240) * (p ^ 3 * q ^ 7) + (-1 / 2880) * (p ^ 4 * q ^ 5) + (1 / 17280) * (p ^ 4 * q ^ 6) + (-1 / 120960) * (p ^ 4 * q ^ 7) + (-1 / 2880) * (p ^ 5 * q ^ 4) + (1 / 14400) * (p ^ 5 * q ^ 5) + (-1 / 86400) * (p ^ 5 * q ^ 6) + (1 / 604800) * (p ^ 5 * q ^ 7) + (-1 / 4320) * (p ^ 6 * q ^ 3) + (1 / 17280) * (p ^ 6 * q ^ 4) + (-1 / 86400) * (p ^ 6 * q ^ 5)) ((1 / 518400) * (p ^ 6 * q ^ 6))
  have c14 := abs_add_le ((-1 / 10080) * (p ^ 2 * q ^ 7) + (-1 / 4320) * (p ^ 3 * q ^ 6) + (1 / 30240) * (p ^ 3 * q ^ 7) + (-1 / 2880) * (p ^ 4 * q ^ 5) + (1 / 17280) * (p ^ 4 * q ^ 6) + (-1 / 120960) * (p ^ 4 * q ^ 7) + (-1 / 2880) * (p ^ 5 * q ^ 4) + (1 / 14400) * (p ^ 5 * q ^ 5) + (-1 / 86400) * (p ^ 5 * q ^ 6) + (1 / 604800) * (p ^ 5 * q ^ 7) + (-1 / 4320) * (p ^ 6 * q ^ 3) + (1 / 17280) * (p ^ 6 * q ^ 4) + (-1 / 86400) * (p ^ 6 * q ^ 5) + (1 / 518400) * (p ^ 6 * q ^ 6)) ((-1 / 3628800) * (p ^ 6 * q ^ 7))
  have c15 := abs_add_le ((-1 / 10080) * (p ^ 2 * q ^ 7) + (-1 / 4320) * (p ^ 3 * q ^ 6) + (1 / 30240) * (p ^ 3 * q ^ 7) + (-1 / 2880) * (p ^ 4 * q ^ 5) + (1 / 17280) * (p ^ 4 * q ^ 6) + (-1 / 120960) * (p ^ 4 * q ^ 7) + (-1 / 2880) * (p ^ 5 * q ^ 4) + (1 / 14400) * (p ^ 5 * q ^ 5) + (-1 / 86400) * (p ^ 5 * q ^ 6) + (1 / 604800) * (p ^ 5 * q ^ 7) + (-1 / 4320) * (p ^ 6 * q ^ 3) + (1 / 17280) * (p ^ 6 * q ^ 4) + (-1 / 86400) * (p ^ 6 * q ^ 5) + (1 / 518400) * (p ^ 6 * q ^ 6) + (-1 / 3628800) * (p ^ 6 * q ^ 7)) ((-1 / 10080) * (p ^ 7 * q ^ 2))
  have c16 := abs_add_le ((-1 / 10080) * (p ^ 2 * q ^ 7) + (-1 / 4320) * (p ^ 3 * q ^ 6) + (1 / 30240) * (p ^ 3 * q ^ 7) + (-1 / 2880) * (p ^ 4 * q ^ 5) + (1 / 17280) * (p ^ 4 * q ^ 6) + (-1 / 120960) * (p ^ 4 * q ^ 7) + (-1 / 2880) * (p ^ 5 * q ^ 4) + (1 / 14400) * (p ^ 5 * q ^ 5) + (-1 / 86400) * (p ^ 5 * q ^ 6) + (1 / 604800) * (p ^ 5 * q ^ 7) + (-1 / 4320) * (p ^ 6 * q ^ 3) + (1 / 17280) * (p ^ 6 * q ^ 4) + (-1 / 86400) * (p ^ 6 * q ^ 5) + (1 / 518400) * (p ^ 6 * q ^ 6) + (-1 / 3628800) * (p ^ 6 * q ^ 7) + (-1 / 10080) * (p ^ 7 * q ^ 2)) ((1 / 30240) * (p ^ 7 * q ^ 3))
  have c17 := abs_add_le ((-1 / 10080) * (p ^ 2 * q ^ 7) + (-1 / 4320) * (p ^ 3 * q ^ 6) + (1 / 30240) * (p ^ 3 * q ^ 7) + (-1 / 2880) * (p ^ 4 * q ^ 5) + (1 / 17280) * (p ^ 4 * q ^ 6) + (-1 / 120960) * (p ^ 4 * q ^ 7) + (-1 / 2880) * (p ^ 5 * q ^ 4) + (1 / 14400) * (p ^ 5 * q ^ 5) + (-1 / 86400) * (p ^ 5 * q ^ 6) + (1 / 604800) * (p ^ 5 * q ^ 7) + (-1 / 4320) * (p ^ 6 * q ^ 3) + (1 / 17280) * (p ^ 6 * q ^ 4) + (-1 / 86400) * (p ^ 6 * q ^ 5) + (1 / 518400) * (p ^ 6 * q ^ 6) + (-1 / 3628800) * (p ^ 6 * q ^ 7) + (-1 / 10080) * (p ^ 7 * q ^ 2) + (1 / 30240) * (p ^ 7 * q ^ 3)) ((-1 / 120960) * (p ^ 7 * q ^ 4))
  have c18 := abs_add_le ((-1 / 10080) * (p ^ 2 * q ^ 7) + (-1 / 4320) * (p ^ 3 * q ^ 6) + (1 / 30240) * (p ^ 3 * q ^ 7) + (-1 / 2880) * (p ^ 4 * q ^ 5) + (1 / 17280) * (p ^ 4 * q ^ 6) + (-1 / 120960) * (p ^ 4 * q ^ 7) + (-1 / 2880) * (p ^ 5 * q ^ 4) + (1 / 14400) * (p ^ 5 * q ^ 5) + (-1 / 86400) * (p ^ 5 * q ^ 6) + (1 / 604800) * (p ^ 5 * q ^ 7) + (-1 / 4320) * (p ^ 6 * q ^ 3) + (1 / 17280) * (p ^ 6 * q ^ 4) + (-1 / 86400) * (p ^ 6 * q ^ 5) + (1 / 518400) * (p ^ 6 * q ^ 6) + (-1 / 3628800) * (p ^ 6 * q ^ 7) + (-1 / 10080) * (p ^ 7 * q ^ 2) + (1 / 30240) * (p ^ 7 * q ^ 3) + (-1 / 120960) * (p ^ 7 * q ^ 4)) ((1 / 604800) * (p ^ 7 * q ^ 5))
  have c19 := abs_add_le ((-1 / 10080) * (p ^ 2 * q ^ 7) + (-1 / 4320) * (p ^ 3 * q ^ 6) + (1 / 30240) * (p ^ 3 * q ^ 7) + (-1 / 2880) * (p ^ 4 * q ^ 5) + (1 / 17280) * (p ^ 4 * q ^ 6) + (-1 / 120960) * (p ^ 4 * q ^ 7) + (-1 / 2880) * (p ^ 5 * q ^ 4) + (1 / 14400) * (p ^ 5 * q ^ 5) + (-1 / 86400) * (p ^ 5 * q ^ 6) + (1 / 604800) * (p ^ 5 * q ^ 7) + (-1 / 4320) * (p ^ 6 * q ^ 3) + (1 / 17280) * (p ^ 6 * q ^ 4) + (-1 / 86400) * (p ^ 6 * q ^ 5) + (1 / 518400) * (p ^ 6 * q ^ 6) + (-1 / 3628800) * (p ^ 6 * q ^ 7) + (-1 / 10080) * (p ^ 7 * q ^ 2) + (1 / 30240) * (p ^ 7 * q ^ 3) + (-1 / 120960) * (p ^ 7 * q ^ 4) + (1 / 604800) * (p ^ 7 * q ^ 5)) ((-1 / 3628800) * (p ^ 7 * q ^ 6))
  have c20 := abs_add_le ((-1 / 10080) * (p ^ 2 * q ^ 7) + (-1 / 4320) * (p ^ 3 * q ^ 6) + (1 / 30240) * (p ^ 3 * q ^ 7) + (-1 / 2880) * (p ^ 4 * q ^ 5) + (1 / 17280) * (p ^ 4 * q ^ 6) + (-1 / 120960) * (p ^ 4 * q ^ 7) + (-1 / 2880) * (p ^ 5 * q ^ 4) + (1 / 14400) * (p ^ 5 * q ^ 5) + (-1 / 86400) * (p ^ 5 * q ^ 6) + (1 / 604800) * (p ^ 5 * q ^ 7) + (-1 / 4320) * (p ^ 6 * q ^ 3) + (1 / 17280) * (p ^ 6 * q ^ 4) + (-1 / 86400) * (p ^ 6 * q ^ 5) + (1 / 518400) * (p ^ 6 * q ^ 6) + (-1 / 3628800) * (p ^ 6 * q ^ 7) + (-1 / 10080) * (p ^ 7 * q ^ 2) + (1 / 30240) * (p ^ 7 * q ^ 3) + (-1 / 120960) * (p ^ 7 * q ^ 4) + (1 / 604800) * (p ^ 7 * q ^ 5) + (-1 / 3628800) * (p ^ 7 * q ^ 6)) ((1 / 25401600) * (p ^ 7 * q ^ 7))
  have c21 := abs_add_le ((-1 / 10080) * (p ^ 2 * q ^ 7) + (-1 / 4320) * (p ^ 3 * q ^ 6) + (1 / 30240) * (p ^ 3 * q ^ 7) + (-1 / 2880) * (p ^ 4 * q ^ 5) + (1 / 17280) * (p ^ 4 * q ^ 6) + (-1 / 120960) * (p ^ 4 * q ^ 7) + (-1 / 2880) * (p ^ 5 * q ^ 4) + (1 / 14400) * (p ^ 5 * q ^ 5) + (-1 / 86400) * (p ^ 5 * q ^ 6) + (1 / 604800) * (p ^ 5 * q ^ 7) + (-1 / 4320) * (p ^ 6 * q ^ 3) + (1 / 17280) * (p ^ 6 * q ^ 4) + (-1 / 86400) * (p ^ 6 * q ^ 5) + (1 / 518400) * (p ^ 6 * q ^ 6) + (-1 / 3628800) * (p ^ 6 * q ^ 7) + (-1 / 10080) * (p ^ 7 * q ^ 2) + (1 / 30240) * (p ^ 7 * q ^ 3) + (-1 / 120960) * (p ^ 7 * q ^ 4) + (1 / 604800) * (p ^ 7 * q ^ 5) + (-1 / 3628800) * (p ^ 7 * q ^ 6) + (1 / 25401600) * (p ^ 7 * q ^ 7)) ((-1 / 40320) * (p ^ 8 * q ^ 1))
  have c22 := abs_add_le ((-1 / 10080) * (p ^ 2 * q ^ 7) + (-1 / 4320) * (p ^ 3 * q ^ 6) + (1 / 30240) * (p ^ 3 * q ^ 7) + (-1 / 2880) * (p ^ 4 * q ^ 5) + (1 / 17280) * (p ^ 4 * q ^ 6) + (-1 / 120960) * (p ^ 4 * q ^ 7) + (-1 / 2880) * (p ^ 5 * q ^ 4) + (1 / 14400) * (p ^ 5 * q ^ 5) + (-1 / 86400) * (p ^ 5 * q ^ 6) + (1 / 604800) * (p ^ 5 * q ^ 7) + (-1 / 4320) * (p ^ 6 * q ^ 3) + (1 / 17280) * (p ^ 6 * q ^ 4) + (-1 / 86400) * (p ^ 6 * q ^ 5) + (1 / 518400) * (p ^ 6 * q ^ 6) + (-1 / 3628800) * (p ^ 6 * q ^ 7) + (-1 / 10080) * (p ^ 7 * q ^ 2) + (1 / 30240) * (p ^ 7 * q ^ 3) + (-1 / 120960) * (p ^ 7 * q ^ 4) + (1 / 604800) * (p ^ 7 * q ^ 5) + (-1 / 3628800) * (p ^ 7 * q ^ 6) + (1 / 25401600) * (p ^ 7 * q ^ 7) + (-1 / 40320) * (p ^ 8 * q ^ 1)) ((1 / 80640) * (p ^ 8 * q ^ 2))
  have c23 := abs_add_le ((-1 / 10080) * (p ^ 2 * q ^ 7) + (-1 / 4320) * (p ^ 3 * q ^ 6) + (1 / 30240) * (p ^ 3 * q ^ 7) + (-1 / 2880) * (p ^ 4 * q ^ 5) + (1 / 17280) * (p ^ 4 * q ^ 6) + (-1 / 120960) * (p ^ 4 * q ^ 7) + (-1 / 2880) * (p ^ 5 * q ^ 4) + (1 / 14400) * (p ^ 5 * q ^ 5) + (-1 / 86400) * (p ^ 5 * q ^ 6) + (1 / 604800) * (p ^ 5 * q ^ 7) + (-1 / 4320) * (p ^ 6 * q ^ 3) + (1 / 17280) * (p ^ 6 * q ^ 4) + (-1 / 86400) * (p ^ 6 * q ^ 5) + (1 / 518400) * (p ^ 6 * q ^ 6) + (-1 / 3628800) * (p ^ 6 * q ^ 7) + (-1 / 10080) * (p ^ 7 * q ^ 2) + (1 / 30240) * (p ^ 7 * q ^ 3) + (-1 / 120960) * (p ^ 7 * q ^ 4) + (1 / 604800) * (p ^ 7 * q ^ 5) + (-1 / 3628800) * (p ^ 7 * q ^ 6) + (1 / 25401600) * (p ^ 7 * q ^ 7) + (-1 / 40320) * (p ^ 8 * q ^ 1) + (1 / 80640) * (p ^ 8 * q ^ 2)) ((-1 / 241920) * (p ^ 8 * q ^ 3))
  have c24 := abs_add_le ((-1 / 10080) * (p ^ 2 * q ^ 7) + (-1 / 4320) * (p ^ 3 * q ^ 6) + (1 / 30240) * (p ^ 3 * q ^ 7) + (-1 / 2880) * (p ^ 4 * q ^ 5) + (1 / 17280) * (p ^ 4 * q ^ 6) + (-1 / 120960) * (p ^ 4 * q ^ 7) + (-1 / 2880) * (p ^ 5 * q ^ 4) + (1 / 14400) * (p ^ 5 * q ^ 5) + (-1 / 86400) * (p ^ 5 * q ^ 6) + (1 / 604800) * (p ^ 5 * q ^ 7) + (-1 / 4320) * (p ^ 6 * q ^ 3) + (1 / 17280) * (p ^ 6 * q ^ 4) + (-1 / 86400) * (p ^ 6 * q ^ 5) + (1 / 518400) * (p ^ 6 * q ^ 6) + (-1 / 3628800) * (p ^ 6 * q ^ 7) + (-1 / 10080) * (p ^ 7 * q ^ 2) + (1 / 30240) * (p ^ 7 * q ^ 3) + (-1 / 120960) * (p ^ 7 * q ^ 4) + (1 / 604800) * (p ^ 7 * q ^ 5) + (-1 / 3628800) * (p ^ 7 * q ^ 6) + (1 / 25401600) * (p ^ 7 * q ^ 7) + (-1 / 40320) * (p ^ 8 * q ^ 1) + (1 / 80640) * (p ^ 8 * q ^ 2) + (-1 / 241920) * (p ^ 8 * q ^ 3)) ((1 / 967680) * (p ^ 8 * q ^ 4))
  have c25 := abs_add_le ((-1 / 10080) * (p ^ 2 * q ^ 7) + (-1 / 4320) * (p ^ 3 * q ^ 6) + (1 / 30240) * (p ^ 3 * q ^ 7) + (-1 / 2880) * (p ^ 4 * q ^ 5) + (1 / 17280) * (p ^ 4 * q ^ 6) + (-1 / 120960) * (p ^ 4 * q ^ 7) + (-1 / 2880) * (p ^ 5 * q ^ 4) + (1 / 14400) * (p ^ 5 * q ^ 5) + (-1 / 86400) * (p ^ 5 * q ^ 6) + (1 / 604800) * (p ^ 5 * q ^ 7) + (-1 / 4320) * (p ^ 6 * q ^ 3) + (1 / 17280) * (p ^ 6 * q ^ 4) + (-1 / 86400) * (p ^ 6 * q ^ 5) + (1 / 518400) * (p ^ 6 * q ^ 6) + (-1 / 3628800) * (p ^ 6 * q ^ 7) + (-1 / 10080) * (p ^ 7 * q ^ 2) + (1 / 30240) * (p ^ 7 * q ^ 3) + (-1 / 120960) * (p ^ 7 * q ^ 4) + (1 / 604800) * (p ^ 7 * q ^ 5) + (-1 / 3628800) * (p ^ 7 * q ^ 6) + (1 / 25401600) * (p ^ 7 * q ^ 7) + (-1 / 40320) * (p ^ 8 * q ^ 1) + (1 / 80640) * (p ^ 8 * q ^ 2) + (-1 / 241920) * (p ^ 8 * q ^ 3) + (1 / 967680) * (p ^ 8 * q ^ 4)) ((-1 / 4838400) * (p ^ 8 * q ^ 5))
  have c26 := abs_add_le ((-1 / 10080) * (p ^ 2 * q ^ 7) + (-1 / 4320) * (p ^ 3 * q ^ 6) + (1 / 30240) * (p ^ 3 * q ^ 7) + (-1 / 2880) * (p ^ 4 * q ^ 5) + (1 / 17280) * (p ^ 4 * q ^ 6) + (-1 / 120960) * (p ^ 4 * q ^ 7) + (-1 / 2880) * (p ^ 5 * q ^ 4) + (1 / 14400) * (p ^ 5 * q ^ 5) + (-1 / 86400) * (p ^ 5 * q ^ 6) + (1 / 604800) * (p ^ 5 * q ^ 7) + (-1 / 4320) * (p ^ 6 * q ^ 3) + (1 / 17280) * (p ^ 6 * q ^ 4) + (-1 / 86400) * (p ^ 6 * q ^ 5) + (1 / 518400) * (p ^ 6 * q ^ 6) + (-1 / 3628800) * (p ^ 6 * q ^ 7) + (-1 / 10080) * (p ^ 7 * q ^ 2) + (1 / 30240) * (p ^ 7 * q ^ 3) + (-1 / 120960) * (p ^ 7 * q ^ 4) + (1 / 604800) * (p ^ 7 * q ^ 5) + (-1 / 3628800) * (p ^ 7 * q ^ 6) + (1 / 25401600) * (p ^ 7 * q ^ 7) + (-1 / 40320) * (p ^ 8 * q ^ 1) + (1 / 80640) * (p ^ 8 * q ^ 2) + (-1 / 241920) * (p ^ 8 * q ^ 3) + (1 / 967680) * (p ^ 8 * q ^ 4) + (-1 / 4838400) * (p ^ 8 * q ^ 5)) ((1 / 29030400) * (p ^ 8 * q ^ 6))
  have c27 := abs_add_le ((-1 / 10080) * (p ^ 2 * q ^ 7) + (-1 / 4320) * (p ^ 3 * q ^ 6) + (1 / 30240) * (p ^ 3 * q ^ 7) + (-1 / 2880) * (p ^ 4 * q ^ 5) + (1 / 17280) * (p ^ 4 * q ^ 6) + (-1 / 120960) * (p ^ 4 * q ^ 7) + (-1 / 2880) * (p ^ 5 * q ^ 4) + (1 / 14400) * (p ^ 5 * q ^ 5) + (-1 / 86400) * (p ^ 5 * q ^ 6) + (1 / 604800) * (p ^ 5 * q ^ 7) + (-1 / 4320) * (p ^ 6 * q ^ 3) + (1 / 17280) * (p ^ 6 * q ^ 4) + (-1 / 86400) * (p ^ 6 * q ^ 5) + (1 / 518400) * (p ^ 6 * q ^ 6) + (-1 / 3628800) * (p ^ 6 * q ^ 7) + (-1 / 10080) * (p ^ 7 * q ^ 2) + (1 / 30240) * (p ^ 7 * q ^ 3) + (-1 / 120960) * (p ^ 7 * q ^ 4) + (1 / 604800) * (p ^ 7 * q ^ 5) + (-1 / 3628800) * (p ^ 7 * q ^ 6) + (1 / 25401600) * (p ^ 7 * q ^ 7) + (-1 / 40320) * (p ^ 8 * q ^ 1) + (1 / 80640) * (p ^ 8 * q ^ 2) + (-1 / 241920) * (p ^ 8 * q ^ 3) + (1 / 967680) * (p ^ 8 * q ^ 4) + (-1 / 4838400) * (p ^ 8 * q ^ 5) + (1 / 29030400) * (p ^ 8 * q ^ 6)) ((-1 / 203212800) * (p ^ 8 * q ^ 7))
  have bHT : |HT| ≤ 28 * t ^ 9 * B10 s := by rw [hHT]; linarith
  have hAb : |A| ≤ 2 * p := by
    have p2 : p ^ 2 ≤ p := pow_le_of_le_one hp0 hp1 (by norm_num)
    have z2 : 0 ≤ p ^ 2 := by positivity
    have p3 : p ^ 3 ≤ p := pow_le_of_le_one hp0 hp1 (by norm_num)
    have z3 : 0 ≤ p ^ 3 := by positivity
    have p4 : p ^ 4 ≤ p := pow_le_of_le_one hp0 hp1 (by norm_num)
    have z4 : 0 ≤ p ^ 4 := by positivity
    have p5 : p ^ 5 ≤ p := pow_le_of_le_one hp0 hp1 (by norm_num)
    have z5 : 0 ≤ p ^ 5 := by positivity
    have p6 : p ^ 6 ≤ p := pow_le_of_le_one hp0 hp1 (by norm_num)
    have z6 : 0 ≤ p ^ 6 := by positivity
    have p7 : p ^ 7 ≤ p := pow_le_of_le_one hp0 hp1 (by norm_num)
    have z7 : 0 ≤ p ^ 7 := by positivity
    have p8 : p ^ 8 ≤ p := pow_le_of_le_one hp0 hp1 (by norm_num)
    have z8 : 0 ≤ p ^ 8 := by positivity
    rw [abs_le, hA]; constructor <;> linarith
  have bAr : |A * r2| ≤ 8 * t ^ 9 * B10 s := by
    rw [abs_mul]
    have h1 : |A| * |r2| ≤ 2 * p * (4 * q ^ 8) :=
      mul_le_mul hAb hr2b (abs_nonneg _) (by positivity)
    have h2 : |(1 : ℝ) * (p ^ 1 * q ^ 8)| ≤ t ^ 9 * B10 s :=
      hmono9 hu0 hut hv0 hvt ht1 hl0 hl1 hle hs hp hq (by norm_num) (by norm_num) (by norm_num) (by norm_num)
    rw [one_mul, pow_one, abs_of_nonneg (by positivity)] at h2
    calc |A| * |r2| ≤ 2 * p * (4 * q ^ 8) := h1
      _ = 8 * (p * q ^ 8) := by ring
      _ ≤ 8 * (t ^ 9 * B10 s) := by linarith
      _ = _ := by ring
  have br1 : |r1 * Real.exp (-q)| ≤ t ^ 9 * B10 s := by
    rw [abs_mul, abs_of_pos (Real.exp_pos _)]
    have h1 : Real.exp (-q) ≤ 1 := Real.exp_le_one_iff.2 (by linarith)
    have h2 : |(1 : ℝ) * (p ^ 9 * q ^ 0)| ≤ t ^ 9 * B10 s :=
      hmono9 hu0 hut hv0 hvt ht1 hl0 hl1 hle hs hp hq (by norm_num) (by norm_num) (by norm_num) (by norm_num)
    rw [one_mul, pow_zero, mul_one, abs_of_nonneg (by positivity)] at h2
    calc |r1| * Real.exp (-q) ≤ p ^ 9 * 1 :=
          mul_le_mul hr1b h1 (Real.exp_pos _).le (by positivity)
      _ ≤ _ := by linarith
  have t1 := abs_add_le (HT + A * r2) (r1 * Real.exp (-q))
  have t2 := abs_add_le HT (A * r2)
  have : 0 ≤ t ^ 9 * B10 s := by positivity
  linarith

lemma int36 {f1 f2 f3 f4 f5 f6 f7 f8 f9 f10 f11 f12 f13 f14 f15 f16 f17 f18 f19 f20 f21 f22 f23 f24 f25 f26 f27 f28 f29 f30 f31 f32 f33 f34 f35 f36 : ℝ → ℝ}
    (c1 c2 c3 c4 c5 c6 c7 c8 c9 c10 c11 c12 c13 c14 c15 c16 c17 c18 c19 c20 c21 c22 c23 c24 c25 c26 c27 c28 c29 c30 c31 c32 c33 c34 c35 c36 : ℝ)
    (h1 : IntegrableOn f1 (Set.Ioi 0))
    (h2 : IntegrableOn f2 (Set.Ioi 0))
    (h3 : IntegrableOn f3 (Set.Ioi 0))
    (h4 : IntegrableOn f4 (Set.Ioi 0))
    (h5 : IntegrableOn f5 (Set.Ioi 0))
    (h6 : IntegrableOn f6 (Set.Ioi 0))
    (h7 : IntegrableOn f7 (Set.Ioi 0))
    (h8 : IntegrableOn f8 (Set.Ioi 0))
    (h9 : IntegrableOn f9 (Set.Ioi 0))
    (h10 : IntegrableOn f10 (Set.Ioi 0))
    (h11 : IntegrableOn f11 (Set.Ioi 0))
    (h12 : IntegrableOn f12 (Set.Ioi 0))
    (h13 : IntegrableOn f13 (Set.Ioi 0))
    (h14 : IntegrableOn f14 (Set.Ioi 0))
    (h15 : IntegrableOn f15 (Set.Ioi 0))
    (h16 : IntegrableOn f16 (Set.Ioi 0))
    (h17 : IntegrableOn f17 (Set.Ioi 0))
    (h18 : IntegrableOn f18 (Set.Ioi 0))
    (h19 : IntegrableOn f19 (Set.Ioi 0))
    (h20 : IntegrableOn f20 (Set.Ioi 0))
    (h21 : IntegrableOn f21 (Set.Ioi 0))
    (h22 : IntegrableOn f22 (Set.Ioi 0))
    (h23 : IntegrableOn f23 (Set.Ioi 0))
    (h24 : IntegrableOn f24 (Set.Ioi 0))
    (h25 : IntegrableOn f25 (Set.Ioi 0))
    (h26 : IntegrableOn f26 (Set.Ioi 0))
    (h27 : IntegrableOn f27 (Set.Ioi 0))
    (h28 : IntegrableOn f28 (Set.Ioi 0))
    (h29 : IntegrableOn f29 (Set.Ioi 0))
    (h30 : IntegrableOn f30 (Set.Ioi 0))
    (h31 : IntegrableOn f31 (Set.Ioi 0))
    (h32 : IntegrableOn f32 (Set.Ioi 0))
    (h33 : IntegrableOn f33 (Set.Ioi 0))
    (h34 : IntegrableOn f34 (Set.Ioi 0))
    (h35 : IntegrableOn f35 (Set.Ioi 0))
    (h36 : IntegrableOn f36 (Set.Ioi 0))
    :
    IntegrableOn (fun s => c1 * f1 s + c2 * f2 s + c3 * f3 s + c4 * f4 s + c5 * f5 s + c6 * f6 s + c7 * f7 s + c8 * f8 s + c9 * f9 s + c10 * f10 s + c11 * f11 s + c12 * f12 s + c13 * f13 s + c14 * f14 s + c15 * f15 s + c16 * f16 s + c17 * f17 s + c18 * f18 s + c19 * f19 s + c20 * f20 s + c21 * f21 s + c22 * f22 s + c23 * f23 s + c24 * f24 s + c25 * f25 s + c26 * f26 s + c27 * f27 s + c28 * f28 s + c29 * f29 s + c30 * f30 s + c31 * f31 s + c32 * f32 s + c33 * f33 s + c34 * f34 s + c35 * f35 s + c36 * f36 s) (Set.Ioi 0) ∧
    ∫ s in Set.Ioi (0 : ℝ), (c1 * f1 s + c2 * f2 s + c3 * f3 s + c4 * f4 s + c5 * f5 s + c6 * f6 s + c7 * f7 s + c8 * f8 s + c9 * f9 s + c10 * f10 s + c11 * f11 s + c12 * f12 s + c13 * f13 s + c14 * f14 s + c15 * f15 s + c16 * f16 s + c17 * f17 s + c18 * f18 s + c19 * f19 s + c20 * f20 s + c21 * f21 s + c22 * f22 s + c23 * f23 s + c24 * f24 s + c25 * f25 s + c26 * f26 s + c27 * f27 s + c28 * f28 s + c29 * f29 s + c30 * f30 s + c31 * f31 s + c32 * f32 s + c33 * f33 s + c34 * f34 s + c35 * f35 s + c36 * f36 s)
      = c1 * (∫ s in Set.Ioi (0 : ℝ), f1 s) + c2 * (∫ s in Set.Ioi (0 : ℝ), f2 s) + c3 * (∫ s in Set.Ioi (0 : ℝ), f3 s) + c4 * (∫ s in Set.Ioi (0 : ℝ), f4 s) + c5 * (∫ s in Set.Ioi (0 : ℝ), f5 s) + c6 * (∫ s in Set.Ioi (0 : ℝ), f6 s) + c7 * (∫ s in Set.Ioi (0 : ℝ), f7 s) + c8 * (∫ s in Set.Ioi (0 : ℝ), f8 s) + c9 * (∫ s in Set.Ioi (0 : ℝ), f9 s) + c10 * (∫ s in Set.Ioi (0 : ℝ), f10 s) + c11 * (∫ s in Set.Ioi (0 : ℝ), f11 s) + c12 * (∫ s in Set.Ioi (0 : ℝ), f12 s) + c13 * (∫ s in Set.Ioi (0 : ℝ), f13 s) + c14 * (∫ s in Set.Ioi (0 : ℝ), f14 s) + c15 * (∫ s in Set.Ioi (0 : ℝ), f15 s) + c16 * (∫ s in Set.Ioi (0 : ℝ), f16 s) + c17 * (∫ s in Set.Ioi (0 : ℝ), f17 s) + c18 * (∫ s in Set.Ioi (0 : ℝ), f18 s) + c19 * (∫ s in Set.Ioi (0 : ℝ), f19 s) + c20 * (∫ s in Set.Ioi (0 : ℝ), f20 s) + c21 * (∫ s in Set.Ioi (0 : ℝ), f21 s) + c22 * (∫ s in Set.Ioi (0 : ℝ), f22 s) + c23 * (∫ s in Set.Ioi (0 : ℝ), f23 s) + c24 * (∫ s in Set.Ioi (0 : ℝ), f24 s) + c25 * (∫ s in Set.Ioi (0 : ℝ), f25 s) + c26 * (∫ s in Set.Ioi (0 : ℝ), f26 s) + c27 * (∫ s in Set.Ioi (0 : ℝ), f27 s) + c28 * (∫ s in Set.Ioi (0 : ℝ), f28 s) + c29 * (∫ s in Set.Ioi (0 : ℝ), f29 s) + c30 * (∫ s in Set.Ioi (0 : ℝ), f30 s) + c31 * (∫ s in Set.Ioi (0 : ℝ), f31 s) + c32 * (∫ s in Set.Ioi (0 : ℝ), f32 s) + c33 * (∫ s in Set.Ioi (0 : ℝ), f33 s) + c34 * (∫ s in Set.Ioi (0 : ℝ), f34 s) + c35 * (∫ s in Set.Ioi (0 : ℝ), f35 s) + c36 * (∫ s in Set.Ioi (0 : ℝ), f36 s) := by
  obtain ⟨b28, e28⟩ := R02.int28 c1 c2 c3 c4 c5 c6 c7 c8 c9 c10 c11 c12 c13 c14 c15 c16 c17 c18 c19 c20 c21 c22 c23 c24 c25 c26 c27 c28 h1 h2 h3 h4 h5 h6 h7 h8 h9 h10 h11 h12 h13 h14 h15 h16 h17 h18 h19 h20 h21 h22 h23 h24 h25 h26 h27 h28
  have a29 : Integrable (fun s => c29 * f29 s) (volume.restrict (Set.Ioi 0)) := h29.const_mul c29
  have a30 : Integrable (fun s => c30 * f30 s) (volume.restrict (Set.Ioi 0)) := h30.const_mul c30
  have a31 : Integrable (fun s => c31 * f31 s) (volume.restrict (Set.Ioi 0)) := h31.const_mul c31
  have a32 : Integrable (fun s => c32 * f32 s) (volume.restrict (Set.Ioi 0)) := h32.const_mul c32
  have a33 : Integrable (fun s => c33 * f33 s) (volume.restrict (Set.Ioi 0)) := h33.const_mul c33
  have a34 : Integrable (fun s => c34 * f34 s) (volume.restrict (Set.Ioi 0)) := h34.const_mul c34
  have a35 : Integrable (fun s => c35 * f35 s) (volume.restrict (Set.Ioi 0)) := h35.const_mul c35
  have a36 : Integrable (fun s => c36 * f36 s) (volume.restrict (Set.Ioi 0)) := h36.const_mul c36
  have b29 : Integrable (fun s => c1 * f1 s + c2 * f2 s + c3 * f3 s + c4 * f4 s + c5 * f5 s + c6 * f6 s + c7 * f7 s + c8 * f8 s + c9 * f9 s + c10 * f10 s + c11 * f11 s + c12 * f12 s + c13 * f13 s + c14 * f14 s + c15 * f15 s + c16 * f16 s + c17 * f17 s + c18 * f18 s + c19 * f19 s + c20 * f20 s + c21 * f21 s + c22 * f22 s + c23 * f23 s + c24 * f24 s + c25 * f25 s + c26 * f26 s + c27 * f27 s + c28 * f28 s + c29 * f29 s) (volume.restrict (Set.Ioi 0)) := b28.add a29
  have b30 : Integrable (fun s => c1 * f1 s + c2 * f2 s + c3 * f3 s + c4 * f4 s + c5 * f5 s + c6 * f6 s + c7 * f7 s + c8 * f8 s + c9 * f9 s + c10 * f10 s + c11 * f11 s + c12 * f12 s + c13 * f13 s + c14 * f14 s + c15 * f15 s + c16 * f16 s + c17 * f17 s + c18 * f18 s + c19 * f19 s + c20 * f20 s + c21 * f21 s + c22 * f22 s + c23 * f23 s + c24 * f24 s + c25 * f25 s + c26 * f26 s + c27 * f27 s + c28 * f28 s + c29 * f29 s + c30 * f30 s) (volume.restrict (Set.Ioi 0)) := b29.add a30
  have b31 : Integrable (fun s => c1 * f1 s + c2 * f2 s + c3 * f3 s + c4 * f4 s + c5 * f5 s + c6 * f6 s + c7 * f7 s + c8 * f8 s + c9 * f9 s + c10 * f10 s + c11 * f11 s + c12 * f12 s + c13 * f13 s + c14 * f14 s + c15 * f15 s + c16 * f16 s + c17 * f17 s + c18 * f18 s + c19 * f19 s + c20 * f20 s + c21 * f21 s + c22 * f22 s + c23 * f23 s + c24 * f24 s + c25 * f25 s + c26 * f26 s + c27 * f27 s + c28 * f28 s + c29 * f29 s + c30 * f30 s + c31 * f31 s) (volume.restrict (Set.Ioi 0)) := b30.add a31
  have b32 : Integrable (fun s => c1 * f1 s + c2 * f2 s + c3 * f3 s + c4 * f4 s + c5 * f5 s + c6 * f6 s + c7 * f7 s + c8 * f8 s + c9 * f9 s + c10 * f10 s + c11 * f11 s + c12 * f12 s + c13 * f13 s + c14 * f14 s + c15 * f15 s + c16 * f16 s + c17 * f17 s + c18 * f18 s + c19 * f19 s + c20 * f20 s + c21 * f21 s + c22 * f22 s + c23 * f23 s + c24 * f24 s + c25 * f25 s + c26 * f26 s + c27 * f27 s + c28 * f28 s + c29 * f29 s + c30 * f30 s + c31 * f31 s + c32 * f32 s) (volume.restrict (Set.Ioi 0)) := b31.add a32
  have b33 : Integrable (fun s => c1 * f1 s + c2 * f2 s + c3 * f3 s + c4 * f4 s + c5 * f5 s + c6 * f6 s + c7 * f7 s + c8 * f8 s + c9 * f9 s + c10 * f10 s + c11 * f11 s + c12 * f12 s + c13 * f13 s + c14 * f14 s + c15 * f15 s + c16 * f16 s + c17 * f17 s + c18 * f18 s + c19 * f19 s + c20 * f20 s + c21 * f21 s + c22 * f22 s + c23 * f23 s + c24 * f24 s + c25 * f25 s + c26 * f26 s + c27 * f27 s + c28 * f28 s + c29 * f29 s + c30 * f30 s + c31 * f31 s + c32 * f32 s + c33 * f33 s) (volume.restrict (Set.Ioi 0)) := b32.add a33
  have b34 : Integrable (fun s => c1 * f1 s + c2 * f2 s + c3 * f3 s + c4 * f4 s + c5 * f5 s + c6 * f6 s + c7 * f7 s + c8 * f8 s + c9 * f9 s + c10 * f10 s + c11 * f11 s + c12 * f12 s + c13 * f13 s + c14 * f14 s + c15 * f15 s + c16 * f16 s + c17 * f17 s + c18 * f18 s + c19 * f19 s + c20 * f20 s + c21 * f21 s + c22 * f22 s + c23 * f23 s + c24 * f24 s + c25 * f25 s + c26 * f26 s + c27 * f27 s + c28 * f28 s + c29 * f29 s + c30 * f30 s + c31 * f31 s + c32 * f32 s + c33 * f33 s + c34 * f34 s) (volume.restrict (Set.Ioi 0)) := b33.add a34
  have b35 : Integrable (fun s => c1 * f1 s + c2 * f2 s + c3 * f3 s + c4 * f4 s + c5 * f5 s + c6 * f6 s + c7 * f7 s + c8 * f8 s + c9 * f9 s + c10 * f10 s + c11 * f11 s + c12 * f12 s + c13 * f13 s + c14 * f14 s + c15 * f15 s + c16 * f16 s + c17 * f17 s + c18 * f18 s + c19 * f19 s + c20 * f20 s + c21 * f21 s + c22 * f22 s + c23 * f23 s + c24 * f24 s + c25 * f25 s + c26 * f26 s + c27 * f27 s + c28 * f28 s + c29 * f29 s + c30 * f30 s + c31 * f31 s + c32 * f32 s + c33 * f33 s + c34 * f34 s + c35 * f35 s) (volume.restrict (Set.Ioi 0)) := b34.add a35
  have b36 : Integrable (fun s => c1 * f1 s + c2 * f2 s + c3 * f3 s + c4 * f4 s + c5 * f5 s + c6 * f6 s + c7 * f7 s + c8 * f8 s + c9 * f9 s + c10 * f10 s + c11 * f11 s + c12 * f12 s + c13 * f13 s + c14 * f14 s + c15 * f15 s + c16 * f16 s + c17 * f17 s + c18 * f18 s + c19 * f19 s + c20 * f20 s + c21 * f21 s + c22 * f22 s + c23 * f23 s + c24 * f24 s + c25 * f25 s + c26 * f26 s + c27 * f27 s + c28 * f28 s + c29 * f29 s + c30 * f30 s + c31 * f31 s + c32 * f32 s + c33 * f33 s + c34 * f34 s + c35 * f35 s + c36 * f36 s) (volume.restrict (Set.Ioi 0)) := b35.add a36
  refine ⟨b36, ?_⟩
  rw [integral_add b35 a36, integral_add b34 a35, integral_add b33 a34, integral_add b32 a33, integral_add b31 a32, integral_add b30 a31, integral_add b29 a30, integral_add b28 a29, e28, integral_const_mul, integral_const_mul, integral_const_mul, integral_const_mul, integral_const_mul, integral_const_mul, integral_const_mul, integral_const_mul]

/-- the integrated expansion of `(e^{-uℓ} - 1) e^{-vs}` to order `8` -/
lemma int_expand8 {u v t : ℝ} (hu0 : 0 ≤ u) (hut : u ≤ t) (hv0 : 0 ≤ v) (hvt : v ≤ t) (ht1 : t ≤ 1) :
    IntegrableOn (fun s => (Real.exp (-u * Y04.ell s) - 1) * Real.exp (-v * s)) (Set.Ioi 0) ∧
    |(∫ s in Set.Ioi (0 : ℝ), (Real.exp (-u * Y04.ell s) - 1) * Real.exp (-v * s))
      - ((-1) * (u ^ 1 * v ^ 0) * Y04.I 0 1 + (1) * (u ^ 1 * v ^ 1) * Y04.I 1 1 + (-1 / 2) * (u ^ 1 * v ^ 2) * Y04.I 2 1 + (1 / 6) * (u ^ 1 * v ^ 3) * Y04.I 3 1 + (-1 / 24) * (u ^ 1 * v ^ 4) * Y04.I 4 1 + (1 / 120) * (u ^ 1 * v ^ 5) * Y04.I 5 1 + (-1 / 720) * (u ^ 1 * v ^ 6) * Y04.I 6 1 + (1 / 5040) * (u ^ 1 * v ^ 7) * Y04.I 7 1 + (1 / 2) * (u ^ 2 * v ^ 0) * Y04.I 0 2 + (-1 / 2) * (u ^ 2 * v ^ 1) * Y04.I 1 2 + (1 / 4) * (u ^ 2 * v ^ 2) * Y04.I 2 2 + (-1 / 12) * (u ^ 2 * v ^ 3) * Y04.I 3 2 + (1 / 48) * (u ^ 2 * v ^ 4) * Y04.I 4 2 + (-1 / 240) * (u ^ 2 * v ^ 5) * Y04.I 5 2 + (1 / 1440) * (u ^ 2 * v ^ 6) * Y04.I 6 2 + (-1 / 6) * (u ^ 3 * v ^ 0) * Y04.I 0 3 + (1 / 6) * (u ^ 3 * v ^ 1) * Y04.I 1 3 + (-1 / 12) * (u ^ 3 * v ^ 2) * Y04.I 2 3 + (1 / 36) * (u ^ 3 * v ^ 3) * Y04.I 3 3 + (-1 / 144) * (u ^ 3 * v ^ 4) * Y04.I 4 3 + (1 / 720) * (u ^ 3 * v ^ 5) * Y04.I 5 3 + (1 / 24) * (u ^ 4 * v ^ 0) * Y04.I 0 4 + (-1 / 24) * (u ^ 4 * v ^ 1) * Y04.I 1 4 + (1 / 48) * (u ^ 4 * v ^ 2) * Y04.I 2 4 + (-1 / 144) * (u ^ 4 * v ^ 3) * Y04.I 3 4 + (1 / 576) * (u ^ 4 * v ^ 4) * Y04.I 4 4 + (-1 / 120) * (u ^ 5 * v ^ 0) * Y04.I 0 5 + (1 / 120) * (u ^ 5 * v ^ 1) * Y04.I 1 5 + (-1 / 240) * (u ^ 5 * v ^ 2) * Y04.I 2 5 + (1 / 720) * (u ^ 5 * v ^ 3) * Y04.I 3 5 + (1 / 720) * (u ^ 6 * v ^ 0) * Y04.I 0 6 + (-1 / 720) * (u ^ 6 * v ^ 1) * Y04.I 1 6 + (1 / 1440) * (u ^ 6 * v ^ 2) * Y04.I 2 6 + (-1 / 5040) * (u ^ 7 * v ^ 0) * Y04.I 0 7 + (1 / 5040) * (u ^ 7 * v ^ 1) * Y04.I 1 7 + (1 / 40320) * (u ^ 8 * v ^ 0) * Y04.I 0 8)|
      ≤ 37 * t ^ 9 * ∫ s in Set.Ioi (0 : ℝ), B10 s := by
  have hint : IntegrableOn (fun s => (Real.exp (-u * Y04.ell s) - 1) * Real.exp (-v * s)) (Set.Ioi 0) :=
    (Z02.int_expand hu0 hut hv0 hvt ht1).1
  refine ⟨hint, ?_⟩
  obtain ⟨hmi, hmv⟩ := int36 ((-1) * (u ^ 1 * v ^ 0)) ((1) * (u ^ 1 * v ^ 1)) ((-1 / 2) * (u ^ 1 * v ^ 2)) ((1 / 6) * (u ^ 1 * v ^ 3)) ((-1 / 24) * (u ^ 1 * v ^ 4)) ((1 / 120) * (u ^ 1 * v ^ 5)) ((-1 / 720) * (u ^ 1 * v ^ 6)) ((1 / 5040) * (u ^ 1 * v ^ 7)) ((1 / 2) * (u ^ 2 * v ^ 0)) ((-1 / 2) * (u ^ 2 * v ^ 1)) ((1 / 4) * (u ^ 2 * v ^ 2)) ((-1 / 12) * (u ^ 2 * v ^ 3)) ((1 / 48) * (u ^ 2 * v ^ 4)) ((-1 / 240) * (u ^ 2 * v ^ 5)) ((1 / 1440) * (u ^ 2 * v ^ 6)) ((-1 / 6) * (u ^ 3 * v ^ 0)) ((1 / 6) * (u ^ 3 * v ^ 1)) ((-1 / 12) * (u ^ 3 * v ^ 2)) ((1 / 36) * (u ^ 3 * v ^ 3)) ((-1 / 144) * (u ^ 3 * v ^ 4)) ((1 / 720) * (u ^ 3 * v ^ 5)) ((1 / 24) * (u ^ 4 * v ^ 0)) ((-1 / 24) * (u ^ 4 * v ^ 1)) ((1 / 48) * (u ^ 4 * v ^ 2)) ((-1 / 144) * (u ^ 4 * v ^ 3)) ((1 / 576) * (u ^ 4 * v ^ 4)) ((-1 / 120) * (u ^ 5 * v ^ 0)) ((1 / 120) * (u ^ 5 * v ^ 1)) ((-1 / 240) * (u ^ 5 * v ^ 2)) ((1 / 720) * (u ^ 5 * v ^ 3)) ((1 / 720) * (u ^ 6 * v ^ 0)) ((-1 / 720) * (u ^ 6 * v ^ 1)) ((1 / 1440) * (u ^ 6 * v ^ 2)) ((-1 / 5040) * (u ^ 7 * v ^ 0)) ((1 / 5040) * (u ^ 7 * v ^ 1)) ((1 / 40320) * (u ^ 8 * v ^ 0))
    (Y04.I_int 0 1 le_rfl) (Y04.I_int 1 1 le_rfl) (Y04.I_int 2 1 le_rfl) (Y04.I_int 3 1 le_rfl) (Y04.I_int 4 1 le_rfl) (Y04.I_int 5 1 le_rfl) (Y04.I_int 6 1 le_rfl) (Y04.I_int 7 1 le_rfl) (Y04.I_int 0 2 (by norm_num)) (Y04.I_int 1 2 (by norm_num)) (Y04.I_int 2 2 (by norm_num)) (Y04.I_int 3 2 (by norm_num)) (Y04.I_int 4 2 (by norm_num)) (Y04.I_int 5 2 (by norm_num)) (Y04.I_int 6 2 (by norm_num)) (Y04.I_int 0 3 (by norm_num)) (Y04.I_int 1 3 (by norm_num)) (Y04.I_int 2 3 (by norm_num)) (Y04.I_int 3 3 (by norm_num)) (Y04.I_int 4 3 (by norm_num)) (Y04.I_int 5 3 (by norm_num)) (Y04.I_int 0 4 (by norm_num)) (Y04.I_int 1 4 (by norm_num)) (Y04.I_int 2 4 (by norm_num)) (Y04.I_int 3 4 (by norm_num)) (Y04.I_int 4 4 (by norm_num)) (Y04.I_int 0 5 (by norm_num)) (Y04.I_int 1 5 (by norm_num)) (Y04.I_int 2 5 (by norm_num)) (Y04.I_int 3 5 (by norm_num)) (Y04.I_int 0 6 (by norm_num)) (Y04.I_int 1 6 (by norm_num)) (Y04.I_int 2 6 (by norm_num)) (Y04.I_int 0 7 (by norm_num)) (Y04.I_int 1 7 (by norm_num)) (Y04.I_int 0 8 (by norm_num))
  unfold Y04.I
  rw [← hmv, ← integral_sub hint hmi]
  have hB := (B10_int.const_mul 37).const_mul (t ^ 9)
  have := norm_integral_le_of_norm_le hB
    ((ae_restrict_iff' measurableSet_Ioi).2 (Filter.Eventually.of_forall
      (fun s (hs : 0 < s) => by
        rw [Real.norm_eq_abs]
        have := pt8 hu0 hut hv0 hvt ht1 hs.le
        unfold P8 at this
        calc _ = _ := rfl
          _ ≤ _ := this
          _ = t ^ 9 * (37 * B10 s) := by ring)))
  rw [Real.norm_eq_abs, integral_const_mul, integral_const_mul] at this
  calc _ ≤ _ := this
    _ = _ := by ring


/-- eight coefficients vanish if the octic is `O(t⁹)` on `(0, d]` -/
lemma octic_zero {d1 d2 d3 d4 d5 d6 d7 d8 K d : ℝ} (hd : 0 < d) (hd1 : d ≤ 1) (hK : 0 ≤ K)
    (h : ∀ t : ℝ, 0 < t → t ≤ d →
      |d1 * t + d2 * t ^ 2 + d3 * t ^ 3 + d4 * t ^ 4 + d5 * t ^ 5 + d6 * t ^ 6 + d7 * t ^ 7 + d8 * t ^ 8| ≤ K * t ^ 9) :
    d1 = 0 ∧ d2 = 0 ∧ d3 = 0 ∧ d4 = 0 ∧ d5 = 0 ∧ d6 = 0 ∧ d7 = 0 ∧ d8 = 0 := by
  have h1 : d1 = 0 := by
    refine Y04.eq_zero_of_le_lin hd (K := |d2| + |d3| + |d4| + |d5| + |d6| + |d7| + |d8| + K)
      (by have := abs_nonneg d2; have := abs_nonneg d3; have := abs_nonneg d4; have := abs_nonneg d5; have := abs_nonneg d6; have := abs_nonneg d7; have := abs_nonneg d8
          linarith)
      (fun t ht0 htd => ?_)
    have ht1 : t ≤ 1 := htd.trans hd1
    have := h t ht0 htd
    set E := d1 * t + d2 * t ^ 2 + d3 * t ^ 3 + d4 * t ^ 4 + d5 * t ^ 5 + d6 * t ^ 6 + d7 * t ^ 7 + d8 * t ^ 8 with hE
    have e : d1 = E / t - d2 * t - d3 * t ^ 2 - d4 * t ^ 3 - d5 * t ^ 4 - d6 * t ^ 5 - d7 * t ^ 6 - d8 * t ^ 7 := by
      rw [hE]; field_simp; ring
    rw [e]
    have ht2 : t ^ 2 ≤ t := pow_le_of_le_one ht0.le ht1 (by norm_num)
    have ht3 : t ^ 3 ≤ t := pow_le_of_le_one ht0.le ht1 (by norm_num)
    have ht4 : t ^ 4 ≤ t := pow_le_of_le_one ht0.le ht1 (by norm_num)
    have ht5 : t ^ 5 ≤ t := pow_le_of_le_one ht0.le ht1 (by norm_num)
    have ht6 : t ^ 6 ≤ t := pow_le_of_le_one ht0.le ht1 (by norm_num)
    have ht7 : t ^ 7 ≤ t := pow_le_of_le_one ht0.le ht1 (by norm_num)
    have ht8 : t ^ 8 ≤ t := pow_le_of_le_one ht0.le ht1 (by norm_num)
    have b1 : |E / t| ≤ K * t ^ 8 := by
      rw [abs_div, abs_of_pos ht0, div_le_iff₀ ht0]; nlinarith
    have b2 : |d2 * t| ≤ |d2| * t := by rw [abs_mul, abs_of_pos ht0]
    have b3 : |d3 * t ^ 2| ≤ |d3| * t := by
      rw [abs_mul, abs_of_pos (pow_pos ht0 2)]
      exact mul_le_mul_of_nonneg_left ht2 (abs_nonneg _)
    have b4 : |d4 * t ^ 3| ≤ |d4| * t := by
      rw [abs_mul, abs_of_pos (pow_pos ht0 3)]
      exact mul_le_mul_of_nonneg_left ht3 (abs_nonneg _)
    have b5 : |d5 * t ^ 4| ≤ |d5| * t := by
      rw [abs_mul, abs_of_pos (pow_pos ht0 4)]
      exact mul_le_mul_of_nonneg_left ht4 (abs_nonneg _)
    have b6 : |d6 * t ^ 5| ≤ |d6| * t := by
      rw [abs_mul, abs_of_pos (pow_pos ht0 5)]
      exact mul_le_mul_of_nonneg_left ht5 (abs_nonneg _)
    have b7 : |d7 * t ^ 6| ≤ |d7| * t := by
      rw [abs_mul, abs_of_pos (pow_pos ht0 6)]
      exact mul_le_mul_of_nonneg_left ht6 (abs_nonneg _)
    have b8 : |d8 * t ^ 7| ≤ |d8| * t := by
      rw [abs_mul, abs_of_pos (pow_pos ht0 7)]
      exact mul_le_mul_of_nonneg_left ht7 (abs_nonneg _)
    have t8 : K * t ^ 8 ≤ K * t := mul_le_mul_of_nonneg_left ht8 hK
    have := abs_sub (E / t - d2 * t - d3 * t ^ 2 - d4 * t ^ 3 - d5 * t ^ 4 - d6 * t ^ 5 - d7 * t ^ 6) (d8 * t ^ 7)
    have := abs_sub (E / t - d2 * t - d3 * t ^ 2 - d4 * t ^ 3 - d5 * t ^ 4 - d6 * t ^ 5) (d7 * t ^ 6)
    have := abs_sub (E / t - d2 * t - d3 * t ^ 2 - d4 * t ^ 3 - d5 * t ^ 4) (d6 * t ^ 5)
    have := abs_sub (E / t - d2 * t - d3 * t ^ 2 - d4 * t ^ 3) (d5 * t ^ 4)
    have := abs_sub (E / t - d2 * t - d3 * t ^ 2) (d4 * t ^ 3)
    have := abs_sub (E / t - d2 * t) (d3 * t ^ 2)
    have := abs_sub (E / t) (d2 * t)
    have : (|d2| + |d3| + |d4| + |d5| + |d6| + |d7| + |d8| + K) * t
        = |d2| * t + |d3| * t + |d4| * t + |d5| * t + |d6| * t + |d7| * t + |d8| * t + K * t := by ring
    linarith
  subst h1
  obtain ⟨b2, b3, b4, b5, b6, b7, b8⟩ := R02.septic_zero hd hd1 hK (d1 := d2) (d2 := d3) (d3 := d4) (d4 := d5) (d5 := d6) (d6 := d7) (d7 := d8)
    (fun t ht0 htd => by
    have := h t ht0 htd
    have e : 0 * t + d2 * t ^ 2 + d3 * t ^ 3 + d4 * t ^ 4 + d5 * t ^ 5 + d6 * t ^ 6 + d7 * t ^ 7 + d8 * t ^ 8
        = t * (d2 * t + d3 * t ^ 2 + d4 * t ^ 3 + d5 * t ^ 4 + d6 * t ^ 5 + d7 * t ^ 6 + d8 * t ^ 7) := by ring
    rw [e, abs_mul, abs_of_pos ht0, show K * t ^ 9 = t * (K * t ^ 8) by ring] at this
    exact le_of_mul_le_mul_left this ht0)
  exact ⟨rfl, b2, b3, b4, b5, b6, b7, b8⟩

lemma series_trunc9 {α β γ x L : ℝ} (hα : |α| ≤ 3 / 2) (hβ : |β| ≤ 3 / 2) (hγ : |γ| ≤ 3 / 2)
    (hx0 : 0 ≤ x) (hx : x ≤ 1 / 3)
    (hL : HasSum (fun k : ℕ => (-1 : ℝ) ^ (k + 2) * Y04.zn (k + 2) / ((k + 2 : ℕ) : ℝ) *
      ((α * x) ^ (k + 2) + (β * x) ^ (k + 2) - (γ * x) ^ (k + 2))) L) :
    |L - (Y04.cf α β γ 0 * x ^ 2 + Y04.cf α β γ 1 * x ^ 3 + Y04.cf α β γ 2 * x ^ 4 + Y04.cf α β γ 3 * x ^ 5 + Y04.cf α β γ 4 * x ^ 6 + Y04.cf α β γ 5 * x ^ 7 + Y04.cf α β γ 6 * x ^ 8 + Y04.cf α β γ 7 * x ^ 9)| ≤ 346 * x ^ 10 := by
  have hf : (fun k : ℕ => (-1 : ℝ) ^ (k + 2) * Y04.zn (k + 2) / ((k + 2 : ℕ) : ℝ) *
      ((α * x) ^ (k + 2) + (β * x) ^ (k + 2) - (γ * x) ^ (k + 2))) =
      fun k => Y04.cf α β γ k * x ^ (k + 2) := by
    ext k; unfold Y04.cf; rw [mul_pow, mul_pow, mul_pow]; ring
  rw [hf] at hL
  set r := 3 / 2 * x with hr
  have hr0 : 0 ≤ r := by positivity
  have hr1 : r ≤ 1 / 2 := by rw [hr]; linarith
  have hbd : ∀ k, |Y04.cf α β γ k * x ^ (k + 2)| ≤ 3 * r ^ (k + 2) := by
    intro k
    rw [abs_mul, abs_pow, hr, mul_pow, abs_of_nonneg hx0]
    have := Y04.cf_bd hα hβ hγ k
    calc |Y04.cf α β γ k| * x ^ (k + 2) ≤ 3 * (3 / 2) ^ (k + 2) * x ^ (k + 2) :=
          mul_le_mul_of_nonneg_right this (by positivity)
      _ = _ := by ring
  have hg := (hasSum_nat_add_iff' 8).mpr hL
  have hgeo : HasSum (fun k : ℕ => 3 * r ^ 10 * r ^ k) (3 * r ^ 10 * (1 - r)⁻¹) :=
    (hasSum_geometric_of_lt_one hr0 (by linarith)).mul_left _
  have hb1 : ∀ k : ℕ, Y04.cf α β γ (k + 8) * x ^ (k + 8 + 2) ≤ 3 * r ^ 10 * r ^ k := by
    intro k
    have := (abs_le.mp (hbd (k + 8))).2
    calc _ ≤ 3 * r ^ (k + 8 + 2) := this
      _ = _ := by ring
  have hb2 : ∀ k : ℕ, -(Y04.cf α β γ (k + 8) * x ^ (k + 8 + 2)) ≤ 3 * r ^ 10 * r ^ k := by
    intro k
    have := (abs_le.mp (hbd (k + 8))).1
    calc _ ≤ 3 * r ^ (k + 8 + 2) := by linarith
      _ = _ := by ring
  have u1 := hasSum_le hb1 hg hgeo
  have u2 := hasSum_le hb2 hg.neg hgeo
  have hinv : 3 * r ^ 10 * (1 - r)⁻¹ ≤ 6 * r ^ 10 := by
    rw [mul_inv_le_iff₀ (by linarith)]
    nlinarith [pow_nonneg hr0 10]
  have hr9 : 6 * r ^ 10 ≤ 346 * x ^ 10 := by
    rw [hr, mul_pow]
    nlinarith [pow_nonneg hx0 10]
  simp only [Finset.sum_range_succ, Finset.sum_range_zero, zero_add] at u1 u2
  rw [abs_le]
  constructor <;> norm_num at u1 u2 ⊢ <;> linarith

/-- chunk 0 of the order-`≥ 10` part of `P + P²/2 + P³/6 + P⁴/24`, divided by `x^10` -/
noncomputable def Wd0 (a0 a1 a2 a3 a4 a5 a6 a7 x : ℝ) : ℝ :=
  (1 / 6) * (a0 ^ 3 * a2)
    + (1 / 4) * (a0 ^ 2 * a1 ^ 2)
    + (1 / 2) * (a0 ^ 2 * a4)
    + (1) * (a0 * a1 * a3)
    + (1 / 2) * (a0 * a2 ^ 2)
    + (1) * (a0 * a6)
    + (1 / 2) * (a1 ^ 2 * a2)
    + (1) * (a1 * a5)
    + (1) * (a2 * a4)
    + (1 / 2) * (a3 ^ 2)
    + (1 / 6) * (a0 ^ 3 * a3 * x)
    + (1 / 2) * (a0 ^ 2 * a1 * a2 * x)
    + (1 / 2) * (a0 ^ 2 * a5 * x)
    + (1 / 6) * (a0 * a1 ^ 3 * x)
    + (1) * (a0 * a1 * a4 * x)
    + (1) * (a0 * a2 * a3 * x)
    + (1) * (a0 * a7 * x)
    + (1 / 2) * (a1 ^ 2 * a3 * x)
    + (1 / 2) * (a1 * a2 ^ 2 * x)
    + (1) * (a1 * a6 * x)
    + (1) * (a2 * a5 * x)
    + (1) * (a3 * a4 * x)
    + (1 / 6) * (a0 ^ 3 * a4 * x ^ 2)
    + (1 / 2) * (a0 ^ 2 * a1 * a3 * x ^ 2)
    + (1 / 4) * (a0 ^ 2 * a2 ^ 2 * x ^ 2)
    + (1 / 2) * (a0 ^ 2 * a6 * x ^ 2)
    + (1 / 2) * (a0 * a1 ^ 2 * a2 * x ^ 2)
    + (1) * (a0 * a1 * a5 * x ^ 2)
    + (1) * (a0 * a2 * a4 * x ^ 2)
    + (1 / 2) * (a0 * a3 ^ 2 * x ^ 2)
    + (1 / 24) * (a1 ^ 4 * x ^ 2)
    + (1 / 2) * (a1 ^ 2 * a4 * x ^ 2)
    + (1) * (a1 * a2 * a3 * x ^ 2)
    + (1) * (a1 * a7 * x ^ 2)
    + (1 / 6) * (a2 ^ 3 * x ^ 2)
    + (1) * (a2 * a6 * x ^ 2)
    + (1) * (a3 * a5 * x ^ 2)
    + (1 / 2) * (a4 ^ 2 * x ^ 2)
    + (1 / 6) * (a0 ^ 3 * a5 * x ^ 3)
    + (1 / 2) * (a0 ^ 2 * a1 * a4 * x ^ 3)

/-- chunk 1 of the order-`≥ 10` part of `P + P²/2 + P³/6 + P⁴/24`, divided by `x^10` -/
noncomputable def Wd1 (a0 a1 a2 a3 a4 a5 a6 a7 x : ℝ) : ℝ :=
  (1 / 2) * (a0 ^ 2 * a2 * a3 * x ^ 3)
    + (1 / 2) * (a0 ^ 2 * a7 * x ^ 3)
    + (1 / 2) * (a0 * a1 ^ 2 * a3 * x ^ 3)
    + (1 / 2) * (a0 * a1 * a2 ^ 2 * x ^ 3)
    + (1) * (a0 * a1 * a6 * x ^ 3)
    + (1) * (a0 * a2 * a5 * x ^ 3)
    + (1) * (a0 * a3 * a4 * x ^ 3)
    + (1 / 6) * (a1 ^ 3 * a2 * x ^ 3)
    + (1 / 2) * (a1 ^ 2 * a5 * x ^ 3)
    + (1) * (a1 * a2 * a4 * x ^ 3)
    + (1 / 2) * (a1 * a3 ^ 2 * x ^ 3)
    + (1 / 2) * (a2 ^ 2 * a3 * x ^ 3)
    + (1) * (a2 * a7 * x ^ 3)
    + (1) * (a3 * a6 * x ^ 3)
    + (1) * (a4 * a5 * x ^ 3)
    + (1 / 6) * (a0 ^ 3 * a6 * x ^ 4)
    + (1 / 2) * (a0 ^ 2 * a1 * a5 * x ^ 4)
    + (1 / 2) * (a0 ^ 2 * a2 * a4 * x ^ 4)
    + (1 / 4) * (a0 ^ 2 * a3 ^ 2 * x ^ 4)
    + (1 / 2) * (a0 * a1 ^ 2 * a4 * x ^ 4)
    + (1) * (a0 * a1 * a2 * a3 * x ^ 4)
    + (1) * (a0 * a1 * a7 * x ^ 4)
    + (1 / 6) * (a0 * a2 ^ 3 * x ^ 4)
    + (1) * (a0 * a2 * a6 * x ^ 4)
    + (1) * (a0 * a3 * a5 * x ^ 4)
    + (1 / 2) * (a0 * a4 ^ 2 * x ^ 4)
    + (1 / 6) * (a1 ^ 3 * a3 * x ^ 4)
    + (1 / 4) * (a1 ^ 2 * a2 ^ 2 * x ^ 4)
    + (1 / 2) * (a1 ^ 2 * a6 * x ^ 4)
    + (1) * (a1 * a2 * a5 * x ^ 4)
    + (1) * (a1 * a3 * a4 * x ^ 4)
    + (1 / 2) * (a2 ^ 2 * a4 * x ^ 4)
    + (1 / 2) * (a2 * a3 ^ 2 * x ^ 4)
    + (1) * (a3 * a7 * x ^ 4)
    + (1) * (a4 * a6 * x ^ 4)
    + (1 / 2) * (a5 ^ 2 * x ^ 4)
    + (1 / 6) * (a0 ^ 3 * a7 * x ^ 5)
    + (1 / 2) * (a0 ^ 2 * a1 * a6 * x ^ 5)
    + (1 / 2) * (a0 ^ 2 * a2 * a5 * x ^ 5)
    + (1 / 2) * (a0 ^ 2 * a3 * a4 * x ^ 5)

/-- chunk 2 of the order-`≥ 10` part of `P + P²/2 + P³/6 + P⁴/24`, divided by `x^10` -/
noncomputable def Wd2 (a0 a1 a2 a3 a4 a5 a6 a7 x : ℝ) : ℝ :=
  (1 / 2) * (a0 * a1 ^ 2 * a5 * x ^ 5)
    + (1) * (a0 * a1 * a2 * a4 * x ^ 5)
    + (1 / 2) * (a0 * a1 * a3 ^ 2 * x ^ 5)
    + (1 / 2) * (a0 * a2 ^ 2 * a3 * x ^ 5)
    + (1) * (a0 * a2 * a7 * x ^ 5)
    + (1) * (a0 * a3 * a6 * x ^ 5)
    + (1) * (a0 * a4 * a5 * x ^ 5)
    + (1 / 6) * (a1 ^ 3 * a4 * x ^ 5)
    + (1 / 2) * (a1 ^ 2 * a2 * a3 * x ^ 5)
    + (1 / 2) * (a1 ^ 2 * a7 * x ^ 5)
    + (1 / 6) * (a1 * a2 ^ 3 * x ^ 5)
    + (1) * (a1 * a2 * a6 * x ^ 5)
    + (1) * (a1 * a3 * a5 * x ^ 5)
    + (1 / 2) * (a1 * a4 ^ 2 * x ^ 5)
    + (1 / 2) * (a2 ^ 2 * a5 * x ^ 5)
    + (1) * (a2 * a3 * a4 * x ^ 5)
    + (1 / 6) * (a3 ^ 3 * x ^ 5)
    + (1) * (a4 * a7 * x ^ 5)
    + (1) * (a5 * a6 * x ^ 5)
    + (1 / 2) * (a0 ^ 2 * a1 * a7 * x ^ 6)
    + (1 / 2) * (a0 ^ 2 * a2 * a6 * x ^ 6)
    + (1 / 2) * (a0 ^ 2 * a3 * a5 * x ^ 6)
    + (1 / 4) * (a0 ^ 2 * a4 ^ 2 * x ^ 6)
    + (1 / 2) * (a0 * a1 ^ 2 * a6 * x ^ 6)
    + (1) * (a0 * a1 * a2 * a5 * x ^ 6)
    + (1) * (a0 * a1 * a3 * a4 * x ^ 6)
    + (1 / 2) * (a0 * a2 ^ 2 * a4 * x ^ 6)
    + (1 / 2) * (a0 * a2 * a3 ^ 2 * x ^ 6)
    + (1) * (a0 * a3 * a7 * x ^ 6)
    + (1) * (a0 * a4 * a6 * x ^ 6)
    + (1 / 2) * (a0 * a5 ^ 2 * x ^ 6)
    + (1 / 6) * (a1 ^ 3 * a5 * x ^ 6)
    + (1 / 2) * (a1 ^ 2 * a2 * a4 * x ^ 6)
    + (1 / 4) * (a1 ^ 2 * a3 ^ 2 * x ^ 6)
    + (1 / 2) * (a1 * a2 ^ 2 * a3 * x ^ 6)
    + (1) * (a1 * a2 * a7 * x ^ 6)
    + (1) * (a1 * a3 * a6 * x ^ 6)
    + (1) * (a1 * a4 * a5 * x ^ 6)
    + (1 / 24) * (a2 ^ 4 * x ^ 6)
    + (1 / 2) * (a2 ^ 2 * a6 * x ^ 6)

/-- chunk 3 of the order-`≥ 10` part of `P + P²/2 + P³/6 + P⁴/24`, divided by `x^10` -/
noncomputable def Wd3 (a0 a1 a2 a3 a4 a5 a6 a7 x : ℝ) : ℝ :=
  (1) * (a2 * a3 * a5 * x ^ 6)
    + (1 / 2) * (a2 * a4 ^ 2 * x ^ 6)
    + (1 / 2) * (a3 ^ 2 * a4 * x ^ 6)
    + (1) * (a5 * a7 * x ^ 6)
    + (1 / 2) * (a6 ^ 2 * x ^ 6)
    + (1 / 2) * (a0 ^ 2 * a2 * a7 * x ^ 7)
    + (1 / 2) * (a0 ^ 2 * a3 * a6 * x ^ 7)
    + (1 / 2) * (a0 ^ 2 * a4 * a5 * x ^ 7)
    + (1 / 2) * (a0 * a1 ^ 2 * a7 * x ^ 7)
    + (1) * (a0 * a1 * a2 * a6 * x ^ 7)
    + (1) * (a0 * a1 * a3 * a5 * x ^ 7)
    + (1 / 2) * (a0 * a1 * a4 ^ 2 * x ^ 7)
    + (1 / 2) * (a0 * a2 ^ 2 * a5 * x ^ 7)
    + (1) * (a0 * a2 * a3 * a4 * x ^ 7)
    + (1 / 6) * (a0 * a3 ^ 3 * x ^ 7)
    + (1) * (a0 * a4 * a7 * x ^ 7)
    + (1) * (a0 * a5 * a6 * x ^ 7)
    + (1 / 6) * (a1 ^ 3 * a6 * x ^ 7)
    + (1 / 2) * (a1 ^ 2 * a2 * a5 * x ^ 7)
    + (1 / 2) * (a1 ^ 2 * a3 * a4 * x ^ 7)
    + (1 / 2) * (a1 * a2 ^ 2 * a4 * x ^ 7)
    + (1 / 2) * (a1 * a2 * a3 ^ 2 * x ^ 7)
    + (1) * (a1 * a3 * a7 * x ^ 7)
    + (1) * (a1 * a4 * a6 * x ^ 7)
    + (1 / 2) * (a1 * a5 ^ 2 * x ^ 7)
    + (1 / 6) * (a2 ^ 3 * a3 * x ^ 7)
    + (1 / 2) * (a2 ^ 2 * a7 * x ^ 7)
    + (1) * (a2 * a3 * a6 * x ^ 7)
    + (1) * (a2 * a4 * a5 * x ^ 7)
    + (1 / 2) * (a3 ^ 2 * a5 * x ^ 7)
    + (1 / 2) * (a3 * a4 ^ 2 * x ^ 7)
    + (1) * (a6 * a7 * x ^ 7)
    + (1 / 2) * (a0 ^ 2 * a3 * a7 * x ^ 8)
    + (1 / 2) * (a0 ^ 2 * a4 * a6 * x ^ 8)
    + (1 / 4) * (a0 ^ 2 * a5 ^ 2 * x ^ 8)
    + (1) * (a0 * a1 * a2 * a7 * x ^ 8)
    + (1) * (a0 * a1 * a3 * a6 * x ^ 8)
    + (1) * (a0 * a1 * a4 * a5 * x ^ 8)
    + (1 / 2) * (a0 * a2 ^ 2 * a6 * x ^ 8)
    + (1) * (a0 * a2 * a3 * a5 * x ^ 8)

/-- chunk 4 of the order-`≥ 10` part of `P + P²/2 + P³/6 + P⁴/24`, divided by `x^10` -/
noncomputable def Wd4 (a0 a1 a2 a3 a4 a5 a6 a7 x : ℝ) : ℝ :=
  (1 / 2) * (a0 * a2 * a4 ^ 2 * x ^ 8)
    + (1 / 2) * (a0 * a3 ^ 2 * a4 * x ^ 8)
    + (1) * (a0 * a5 * a7 * x ^ 8)
    + (1 / 2) * (a0 * a6 ^ 2 * x ^ 8)
    + (1 / 6) * (a1 ^ 3 * a7 * x ^ 8)
    + (1 / 2) * (a1 ^ 2 * a2 * a6 * x ^ 8)
    + (1 / 2) * (a1 ^ 2 * a3 * a5 * x ^ 8)
    + (1 / 4) * (a1 ^ 2 * a4 ^ 2 * x ^ 8)
    + (1 / 2) * (a1 * a2 ^ 2 * a5 * x ^ 8)
    + (1) * (a1 * a2 * a3 * a4 * x ^ 8)
    + (1 / 6) * (a1 * a3 ^ 3 * x ^ 8)
    + (1) * (a1 * a4 * a7 * x ^ 8)
    + (1) * (a1 * a5 * a6 * x ^ 8)
    + (1 / 6) * (a2 ^ 3 * a4 * x ^ 8)
    + (1 / 4) * (a2 ^ 2 * a3 ^ 2 * x ^ 8)
    + (1) * (a2 * a3 * a7 * x ^ 8)
    + (1) * (a2 * a4 * a6 * x ^ 8)
    + (1 / 2) * (a2 * a5 ^ 2 * x ^ 8)
    + (1 / 2) * (a3 ^ 2 * a6 * x ^ 8)
    + (1) * (a3 * a4 * a5 * x ^ 8)
    + (1 / 6) * (a4 ^ 3 * x ^ 8)
    + (1 / 2) * (a7 ^ 2 * x ^ 8)
    + (1 / 2) * (a0 ^ 2 * a4 * a7 * x ^ 9)
    + (1 / 2) * (a0 ^ 2 * a5 * a6 * x ^ 9)
    + (1) * (a0 * a1 * a3 * a7 * x ^ 9)
    + (1) * (a0 * a1 * a4 * a6 * x ^ 9)
    + (1 / 2) * (a0 * a1 * a5 ^ 2 * x ^ 9)
    + (1 / 2) * (a0 * a2 ^ 2 * a7 * x ^ 9)
    + (1) * (a0 * a2 * a3 * a6 * x ^ 9)
    + (1) * (a0 * a2 * a4 * a5 * x ^ 9)
    + (1 / 2) * (a0 * a3 ^ 2 * a5 * x ^ 9)
    + (1 / 2) * (a0 * a3 * a4 ^ 2 * x ^ 9)
    + (1) * (a0 * a6 * a7 * x ^ 9)
    + (1 / 2) * (a1 ^ 2 * a2 * a7 * x ^ 9)
    + (1 / 2) * (a1 ^ 2 * a3 * a6 * x ^ 9)
    + (1 / 2) * (a1 ^ 2 * a4 * a5 * x ^ 9)
    + (1 / 2) * (a1 * a2 ^ 2 * a6 * x ^ 9)
    + (1) * (a1 * a2 * a3 * a5 * x ^ 9)
    + (1 / 2) * (a1 * a2 * a4 ^ 2 * x ^ 9)
    + (1 / 2) * (a1 * a3 ^ 2 * a4 * x ^ 9)

/-- chunk 5 of the order-`≥ 10` part of `P + P²/2 + P³/6 + P⁴/24`, divided by `x^10` -/
noncomputable def Wd5 (a0 a1 a2 a3 a4 a5 a6 a7 x : ℝ) : ℝ :=
  (1) * (a1 * a5 * a7 * x ^ 9)
    + (1 / 2) * (a1 * a6 ^ 2 * x ^ 9)
    + (1 / 6) * (a2 ^ 3 * a5 * x ^ 9)
    + (1 / 2) * (a2 ^ 2 * a3 * a4 * x ^ 9)
    + (1 / 6) * (a2 * a3 ^ 3 * x ^ 9)
    + (1) * (a2 * a4 * a7 * x ^ 9)
    + (1) * (a2 * a5 * a6 * x ^ 9)
    + (1 / 2) * (a3 ^ 2 * a7 * x ^ 9)
    + (1) * (a3 * a4 * a6 * x ^ 9)
    + (1 / 2) * (a3 * a5 ^ 2 * x ^ 9)
    + (1 / 2) * (a4 ^ 2 * a5 * x ^ 9)
    + (1 / 2) * (a0 ^ 2 * a5 * a7 * x ^ 10)
    + (1 / 4) * (a0 ^ 2 * a6 ^ 2 * x ^ 10)
    + (1) * (a0 * a1 * a4 * a7 * x ^ 10)
    + (1) * (a0 * a1 * a5 * a6 * x ^ 10)
    + (1) * (a0 * a2 * a3 * a7 * x ^ 10)
    + (1) * (a0 * a2 * a4 * a6 * x ^ 10)
    + (1 / 2) * (a0 * a2 * a5 ^ 2 * x ^ 10)
    + (1 / 2) * (a0 * a3 ^ 2 * a6 * x ^ 10)
    + (1) * (a0 * a3 * a4 * a5 * x ^ 10)
    + (1 / 6) * (a0 * a4 ^ 3 * x ^ 10)
    + (1 / 2) * (a0 * a7 ^ 2 * x ^ 10)
    + (1 / 2) * (a1 ^ 2 * a3 * a7 * x ^ 10)
    + (1 / 2) * (a1 ^ 2 * a4 * a6 * x ^ 10)
    + (1 / 4) * (a1 ^ 2 * a5 ^ 2 * x ^ 10)
    + (1 / 2) * (a1 * a2 ^ 2 * a7 * x ^ 10)
    + (1) * (a1 * a2 * a3 * a6 * x ^ 10)
    + (1) * (a1 * a2 * a4 * a5 * x ^ 10)
    + (1 / 2) * (a1 * a3 ^ 2 * a5 * x ^ 10)
    + (1 / 2) * (a1 * a3 * a4 ^ 2 * x ^ 10)
    + (1) * (a1 * a6 * a7 * x ^ 10)
    + (1 / 6) * (a2 ^ 3 * a6 * x ^ 10)
    + (1 / 2) * (a2 ^ 2 * a3 * a5 * x ^ 10)
    + (1 / 4) * (a2 ^ 2 * a4 ^ 2 * x ^ 10)
    + (1 / 2) * (a2 * a3 ^ 2 * a4 * x ^ 10)
    + (1) * (a2 * a5 * a7 * x ^ 10)
    + (1 / 2) * (a2 * a6 ^ 2 * x ^ 10)
    + (1 / 24) * (a3 ^ 4 * x ^ 10)
    + (1) * (a3 * a4 * a7 * x ^ 10)
    + (1) * (a3 * a5 * a6 * x ^ 10)

/-- chunk 6 of the order-`≥ 10` part of `P + P²/2 + P³/6 + P⁴/24`, divided by `x^10` -/
noncomputable def Wd6 (a0 a1 a2 a3 a4 a5 a6 a7 x : ℝ) : ℝ :=
  (1 / 2) * (a4 ^ 2 * a6 * x ^ 10)
    + (1 / 2) * (a4 * a5 ^ 2 * x ^ 10)
    + (1 / 2) * (a0 ^ 2 * a6 * a7 * x ^ 11)
    + (1) * (a0 * a1 * a5 * a7 * x ^ 11)
    + (1 / 2) * (a0 * a1 * a6 ^ 2 * x ^ 11)
    + (1) * (a0 * a2 * a4 * a7 * x ^ 11)
    + (1) * (a0 * a2 * a5 * a6 * x ^ 11)
    + (1 / 2) * (a0 * a3 ^ 2 * a7 * x ^ 11)
    + (1) * (a0 * a3 * a4 * a6 * x ^ 11)
    + (1 / 2) * (a0 * a3 * a5 ^ 2 * x ^ 11)
    + (1 / 2) * (a0 * a4 ^ 2 * a5 * x ^ 11)
    + (1 / 2) * (a1 ^ 2 * a4 * a7 * x ^ 11)
    + (1 / 2) * (a1 ^ 2 * a5 * a6 * x ^ 11)
    + (1) * (a1 * a2 * a3 * a7 * x ^ 11)
    + (1) * (a1 * a2 * a4 * a6 * x ^ 11)
    + (1 / 2) * (a1 * a2 * a5 ^ 2 * x ^ 11)
    + (1 / 2) * (a1 * a3 ^ 2 * a6 * x ^ 11)
    + (1) * (a1 * a3 * a4 * a5 * x ^ 11)
    + (1 / 6) * (a1 * a4 ^ 3 * x ^ 11)
    + (1 / 2) * (a1 * a7 ^ 2 * x ^ 11)
    + (1 / 6) * (a2 ^ 3 * a7 * x ^ 11)
    + (1 / 2) * (a2 ^ 2 * a3 * a6 * x ^ 11)
    + (1 / 2) * (a2 ^ 2 * a4 * a5 * x ^ 11)
    + (1 / 2) * (a2 * a3 ^ 2 * a5 * x ^ 11)
    + (1 / 2) * (a2 * a3 * a4 ^ 2 * x ^ 11)
    + (1) * (a2 * a6 * a7 * x ^ 11)
    + (1 / 6) * (a3 ^ 3 * a4 * x ^ 11)
    + (1) * (a3 * a5 * a7 * x ^ 11)
    + (1 / 2) * (a3 * a6 ^ 2 * x ^ 11)
    + (1 / 2) * (a4 ^ 2 * a7 * x ^ 11)
    + (1) * (a4 * a5 * a6 * x ^ 11)
    + (1 / 6) * (a5 ^ 3 * x ^ 11)
    + (1 / 4) * (a0 ^ 2 * a7 ^ 2 * x ^ 12)
    + (1) * (a0 * a1 * a6 * a7 * x ^ 12)
    + (1) * (a0 * a2 * a5 * a7 * x ^ 12)
    + (1 / 2) * (a0 * a2 * a6 ^ 2 * x ^ 12)
    + (1) * (a0 * a3 * a4 * a7 * x ^ 12)
    + (1) * (a0 * a3 * a5 * a6 * x ^ 12)
    + (1 / 2) * (a0 * a4 ^ 2 * a6 * x ^ 12)
    + (1 / 2) * (a0 * a4 * a5 ^ 2 * x ^ 12)

/-- chunk 7 of the order-`≥ 10` part of `P + P²/2 + P³/6 + P⁴/24`, divided by `x^10` -/
noncomputable def Wd7 (a0 a1 a2 a3 a4 a5 a6 a7 x : ℝ) : ℝ :=
  (1 / 2) * (a1 ^ 2 * a5 * a7 * x ^ 12)
    + (1 / 4) * (a1 ^ 2 * a6 ^ 2 * x ^ 12)
    + (1) * (a1 * a2 * a4 * a7 * x ^ 12)
    + (1) * (a1 * a2 * a5 * a6 * x ^ 12)
    + (1 / 2) * (a1 * a3 ^ 2 * a7 * x ^ 12)
    + (1) * (a1 * a3 * a4 * a6 * x ^ 12)
    + (1 / 2) * (a1 * a3 * a5 ^ 2 * x ^ 12)
    + (1 / 2) * (a1 * a4 ^ 2 * a5 * x ^ 12)
    + (1 / 2) * (a2 ^ 2 * a3 * a7 * x ^ 12)
    + (1 / 2) * (a2 ^ 2 * a4 * a6 * x ^ 12)
    + (1 / 4) * (a2 ^ 2 * a5 ^ 2 * x ^ 12)
    + (1 / 2) * (a2 * a3 ^ 2 * a6 * x ^ 12)
    + (1) * (a2 * a3 * a4 * a5 * x ^ 12)
    + (1 / 6) * (a2 * a4 ^ 3 * x ^ 12)
    + (1 / 2) * (a2 * a7 ^ 2 * x ^ 12)
    + (1 / 6) * (a3 ^ 3 * a5 * x ^ 12)
    + (1 / 4) * (a3 ^ 2 * a4 ^ 2 * x ^ 12)
    + (1) * (a3 * a6 * a7 * x ^ 12)
    + (1) * (a4 * a5 * a7 * x ^ 12)
    + (1 / 2) * (a4 * a6 ^ 2 * x ^ 12)
    + (1 / 2) * (a5 ^ 2 * a6 * x ^ 12)
    + (1 / 2) * (a0 * a1 * a7 ^ 2 * x ^ 13)
    + (1) * (a0 * a2 * a6 * a7 * x ^ 13)
    + (1) * (a0 * a3 * a5 * a7 * x ^ 13)
    + (1 / 2) * (a0 * a3 * a6 ^ 2 * x ^ 13)
    + (1 / 2) * (a0 * a4 ^ 2 * a7 * x ^ 13)
    + (1) * (a0 * a4 * a5 * a6 * x ^ 13)
    + (1 / 6) * (a0 * a5 ^ 3 * x ^ 13)
    + (1 / 2) * (a1 ^ 2 * a6 * a7 * x ^ 13)
    + (1) * (a1 * a2 * a5 * a7 * x ^ 13)
    + (1 / 2) * (a1 * a2 * a6 ^ 2 * x ^ 13)
    + (1) * (a1 * a3 * a4 * a7 * x ^ 13)
    + (1) * (a1 * a3 * a5 * a6 * x ^ 13)
    + (1 / 2) * (a1 * a4 ^ 2 * a6 * x ^ 13)
    + (1 / 2) * (a1 * a4 * a5 ^ 2 * x ^ 13)
    + (1 / 2) * (a2 ^ 2 * a4 * a7 * x ^ 13)
    + (1 / 2) * (a2 ^ 2 * a5 * a6 * x ^ 13)
    + (1 / 2) * (a2 * a3 ^ 2 * a7 * x ^ 13)
    + (1) * (a2 * a3 * a4 * a6 * x ^ 13)
    + (1 / 2) * (a2 * a3 * a5 ^ 2 * x ^ 13)

/-- chunk 8 of the order-`≥ 10` part of `P + P²/2 + P³/6 + P⁴/24`, divided by `x^10` -/
noncomputable def Wd8 (a0 a1 a2 a3 a4 a5 a6 a7 x : ℝ) : ℝ :=
  (1 / 2) * (a2 * a4 ^ 2 * a5 * x ^ 13)
    + (1 / 6) * (a3 ^ 3 * a6 * x ^ 13)
    + (1 / 2) * (a3 ^ 2 * a4 * a5 * x ^ 13)
    + (1 / 6) * (a3 * a4 ^ 3 * x ^ 13)
    + (1 / 2) * (a3 * a7 ^ 2 * x ^ 13)
    + (1) * (a4 * a6 * a7 * x ^ 13)
    + (1 / 2) * (a5 ^ 2 * a7 * x ^ 13)
    + (1 / 2) * (a5 * a6 ^ 2 * x ^ 13)
    + (1 / 2) * (a0 * a2 * a7 ^ 2 * x ^ 14)
    + (1) * (a0 * a3 * a6 * a7 * x ^ 14)
    + (1) * (a0 * a4 * a5 * a7 * x ^ 14)
    + (1 / 2) * (a0 * a4 * a6 ^ 2 * x ^ 14)
    + (1 / 2) * (a0 * a5 ^ 2 * a6 * x ^ 14)
    + (1 / 4) * (a1 ^ 2 * a7 ^ 2 * x ^ 14)
    + (1) * (a1 * a2 * a6 * a7 * x ^ 14)
    + (1) * (a1 * a3 * a5 * a7 * x ^ 14)
    + (1 / 2) * (a1 * a3 * a6 ^ 2 * x ^ 14)
    + (1 / 2) * (a1 * a4 ^ 2 * a7 * x ^ 14)
    + (1) * (a1 * a4 * a5 * a6 * x ^ 14)
    + (1 / 6) * (a1 * a5 ^ 3 * x ^ 14)
    + (1 / 2) * (a2 ^ 2 * a5 * a7 * x ^ 14)
    + (1 / 4) * (a2 ^ 2 * a6 ^ 2 * x ^ 14)
    + (1) * (a2 * a3 * a4 * a7 * x ^ 14)
    + (1) * (a2 * a3 * a5 * a6 * x ^ 14)
    + (1 / 2) * (a2 * a4 ^ 2 * a6 * x ^ 14)
    + (1 / 2) * (a2 * a4 * a5 ^ 2 * x ^ 14)
    + (1 / 6) * (a3 ^ 3 * a7 * x ^ 14)
    + (1 / 2) * (a3 ^ 2 * a4 * a6 * x ^ 14)
    + (1 / 4) * (a3 ^ 2 * a5 ^ 2 * x ^ 14)
    + (1 / 2) * (a3 * a4 ^ 2 * a5 * x ^ 14)
    + (1 / 24) * (a4 ^ 4 * x ^ 14)
    + (1 / 2) * (a4 * a7 ^ 2 * x ^ 14)
    + (1) * (a5 * a6 * a7 * x ^ 14)
    + (1 / 6) * (a6 ^ 3 * x ^ 14)
    + (1 / 2) * (a0 * a3 * a7 ^ 2 * x ^ 15)
    + (1) * (a0 * a4 * a6 * a7 * x ^ 15)
    + (1 / 2) * (a0 * a5 ^ 2 * a7 * x ^ 15)
    + (1 / 2) * (a0 * a5 * a6 ^ 2 * x ^ 15)
    + (1 / 2) * (a1 * a2 * a7 ^ 2 * x ^ 15)
    + (1) * (a1 * a3 * a6 * a7 * x ^ 15)

/-- chunk 9 of the order-`≥ 10` part of `P + P²/2 + P³/6 + P⁴/24`, divided by `x^10` -/
noncomputable def Wd9 (a0 a1 a2 a3 a4 a5 a6 a7 x : ℝ) : ℝ :=
  (1) * (a1 * a4 * a5 * a7 * x ^ 15)
    + (1 / 2) * (a1 * a4 * a6 ^ 2 * x ^ 15)
    + (1 / 2) * (a1 * a5 ^ 2 * a6 * x ^ 15)
    + (1 / 2) * (a2 ^ 2 * a6 * a7 * x ^ 15)
    + (1) * (a2 * a3 * a5 * a7 * x ^ 15)
    + (1 / 2) * (a2 * a3 * a6 ^ 2 * x ^ 15)
    + (1 / 2) * (a2 * a4 ^ 2 * a7 * x ^ 15)
    + (1) * (a2 * a4 * a5 * a6 * x ^ 15)
    + (1 / 6) * (a2 * a5 ^ 3 * x ^ 15)
    + (1 / 2) * (a3 ^ 2 * a4 * a7 * x ^ 15)
    + (1 / 2) * (a3 ^ 2 * a5 * a6 * x ^ 15)
    + (1 / 2) * (a3 * a4 ^ 2 * a6 * x ^ 15)
    + (1 / 2) * (a3 * a4 * a5 ^ 2 * x ^ 15)
    + (1 / 6) * (a4 ^ 3 * a5 * x ^ 15)
    + (1 / 2) * (a5 * a7 ^ 2 * x ^ 15)
    + (1 / 2) * (a6 ^ 2 * a7 * x ^ 15)
    + (1 / 2) * (a0 * a4 * a7 ^ 2 * x ^ 16)
    + (1) * (a0 * a5 * a6 * a7 * x ^ 16)
    + (1 / 6) * (a0 * a6 ^ 3 * x ^ 16)
    + (1 / 2) * (a1 * a3 * a7 ^ 2 * x ^ 16)
    + (1) * (a1 * a4 * a6 * a7 * x ^ 16)
    + (1 / 2) * (a1 * a5 ^ 2 * a7 * x ^ 16)
    + (1 / 2) * (a1 * a5 * a6 ^ 2 * x ^ 16)
    + (1 / 4) * (a2 ^ 2 * a7 ^ 2 * x ^ 16)
    + (1) * (a2 * a3 * a6 * a7 * x ^ 16)
    + (1) * (a2 * a4 * a5 * a7 * x ^ 16)
    + (1 / 2) * (a2 * a4 * a6 ^ 2 * x ^ 16)
    + (1 / 2) * (a2 * a5 ^ 2 * a6 * x ^ 16)
    + (1 / 2) * (a3 ^ 2 * a5 * a7 * x ^ 16)
    + (1 / 4) * (a3 ^ 2 * a6 ^ 2 * x ^ 16)
    + (1 / 2) * (a3 * a4 ^ 2 * a7 * x ^ 16)
    + (1) * (a3 * a4 * a5 * a6 * x ^ 16)
    + (1 / 6) * (a3 * a5 ^ 3 * x ^ 16)
    + (1 / 6) * (a4 ^ 3 * a6 * x ^ 16)
    + (1 / 4) * (a4 ^ 2 * a5 ^ 2 * x ^ 16)
    + (1 / 2) * (a6 * a7 ^ 2 * x ^ 16)
    + (1 / 2) * (a0 * a5 * a7 ^ 2 * x ^ 17)
    + (1 / 2) * (a0 * a6 ^ 2 * a7 * x ^ 17)
    + (1 / 2) * (a1 * a4 * a7 ^ 2 * x ^ 17)
    + (1) * (a1 * a5 * a6 * a7 * x ^ 17)

/-- chunk 10 of the order-`≥ 10` part of `P + P²/2 + P³/6 + P⁴/24`, divided by `x^10` -/
noncomputable def Wd10 (a0 a1 a2 a3 a4 a5 a6 a7 x : ℝ) : ℝ :=
  (1 / 6) * (a1 * a6 ^ 3 * x ^ 17)
    + (1 / 2) * (a2 * a3 * a7 ^ 2 * x ^ 17)
    + (1) * (a2 * a4 * a6 * a7 * x ^ 17)
    + (1 / 2) * (a2 * a5 ^ 2 * a7 * x ^ 17)
    + (1 / 2) * (a2 * a5 * a6 ^ 2 * x ^ 17)
    + (1 / 2) * (a3 ^ 2 * a6 * a7 * x ^ 17)
    + (1) * (a3 * a4 * a5 * a7 * x ^ 17)
    + (1 / 2) * (a3 * a4 * a6 ^ 2 * x ^ 17)
    + (1 / 2) * (a3 * a5 ^ 2 * a6 * x ^ 17)
    + (1 / 6) * (a4 ^ 3 * a7 * x ^ 17)
    + (1 / 2) * (a4 ^ 2 * a5 * a6 * x ^ 17)
    + (1 / 6) * (a4 * a5 ^ 3 * x ^ 17)
    + (1 / 6) * (a7 ^ 3 * x ^ 17)
    + (1 / 2) * (a0 * a6 * a7 ^ 2 * x ^ 18)
    + (1 / 2) * (a1 * a5 * a7 ^ 2 * x ^ 18)
    + (1 / 2) * (a1 * a6 ^ 2 * a7 * x ^ 18)
    + (1 / 2) * (a2 * a4 * a7 ^ 2 * x ^ 18)
    + (1) * (a2 * a5 * a6 * a7 * x ^ 18)
    + (1 / 6) * (a2 * a6 ^ 3 * x ^ 18)
    + (1 / 4) * (a3 ^ 2 * a7 ^ 2 * x ^ 18)
    + (1) * (a3 * a4 * a6 * a7 * x ^ 18)
    + (1 / 2) * (a3 * a5 ^ 2 * a7 * x ^ 18)
    + (1 / 2) * (a3 * a5 * a6 ^ 2 * x ^ 18)
    + (1 / 2) * (a4 ^ 2 * a5 * a7 * x ^ 18)
    + (1 / 4) * (a4 ^ 2 * a6 ^ 2 * x ^ 18)
    + (1 / 2) * (a4 * a5 ^ 2 * a6 * x ^ 18)
    + (1 / 24) * (a5 ^ 4 * x ^ 18)
    + (1 / 6) * (a0 * a7 ^ 3 * x ^ 19)
    + (1 / 2) * (a1 * a6 * a7 ^ 2 * x ^ 19)
    + (1 / 2) * (a2 * a5 * a7 ^ 2 * x ^ 19)
    + (1 / 2) * (a2 * a6 ^ 2 * a7 * x ^ 19)
    + (1 / 2) * (a3 * a4 * a7 ^ 2 * x ^ 19)
    + (1) * (a3 * a5 * a6 * a7 * x ^ 19)
    + (1 / 6) * (a3 * a6 ^ 3 * x ^ 19)
    + (1 / 2) * (a4 ^ 2 * a6 * a7 * x ^ 19)
    + (1 / 2) * (a4 * a5 ^ 2 * a7 * x ^ 19)
    + (1 / 2) * (a4 * a5 * a6 ^ 2 * x ^ 19)
    + (1 / 6) * (a5 ^ 3 * a6 * x ^ 19)
    + (1 / 6) * (a1 * a7 ^ 3 * x ^ 20)
    + (1 / 2) * (a2 * a6 * a7 ^ 2 * x ^ 20)

/-- chunk 11 of the order-`≥ 10` part of `P + P²/2 + P³/6 + P⁴/24`, divided by `x^10` -/
noncomputable def Wd11 (a0 a1 a2 a3 a4 a5 a6 a7 x : ℝ) : ℝ :=
  (1 / 2) * (a3 * a5 * a7 ^ 2 * x ^ 20)
    + (1 / 2) * (a3 * a6 ^ 2 * a7 * x ^ 20)
    + (1 / 4) * (a4 ^ 2 * a7 ^ 2 * x ^ 20)
    + (1) * (a4 * a5 * a6 * a7 * x ^ 20)
    + (1 / 6) * (a4 * a6 ^ 3 * x ^ 20)
    + (1 / 6) * (a5 ^ 3 * a7 * x ^ 20)
    + (1 / 4) * (a5 ^ 2 * a6 ^ 2 * x ^ 20)
    + (1 / 6) * (a2 * a7 ^ 3 * x ^ 21)
    + (1 / 2) * (a3 * a6 * a7 ^ 2 * x ^ 21)
    + (1 / 2) * (a4 * a5 * a7 ^ 2 * x ^ 21)
    + (1 / 2) * (a4 * a6 ^ 2 * a7 * x ^ 21)
    + (1 / 2) * (a5 ^ 2 * a6 * a7 * x ^ 21)
    + (1 / 6) * (a5 * a6 ^ 3 * x ^ 21)
    + (1 / 6) * (a3 * a7 ^ 3 * x ^ 22)
    + (1 / 2) * (a4 * a6 * a7 ^ 2 * x ^ 22)
    + (1 / 4) * (a5 ^ 2 * a7 ^ 2 * x ^ 22)
    + (1 / 2) * (a5 * a6 ^ 2 * a7 * x ^ 22)
    + (1 / 24) * (a6 ^ 4 * x ^ 22)
    + (1 / 6) * (a4 * a7 ^ 3 * x ^ 23)
    + (1 / 2) * (a5 * a6 * a7 ^ 2 * x ^ 23)
    + (1 / 6) * (a6 ^ 3 * a7 * x ^ 23)
    + (1 / 6) * (a5 * a7 ^ 3 * x ^ 24)
    + (1 / 4) * (a6 ^ 2 * a7 ^ 2 * x ^ 24)
    + (1 / 6) * (a6 * a7 ^ 3 * x ^ 25)
    + (1 / 24) * (a7 ^ 4 * x ^ 26)

noncomputable def Wpoly9 (a0 a1 a2 a3 a4 a5 a6 a7 x : ℝ) : ℝ :=
  Wd0 a0 a1 a2 a3 a4 a5 a6 a7 x + Wd1 a0 a1 a2 a3 a4 a5 a6 a7 x + Wd2 a0 a1 a2 a3 a4 a5 a6 a7 x + Wd3 a0 a1 a2 a3 a4 a5 a6 a7 x + Wd4 a0 a1 a2 a3 a4 a5 a6 a7 x + Wd5 a0 a1 a2 a3 a4 a5 a6 a7 x + Wd6 a0 a1 a2 a3 a4 a5 a6 a7 x + Wd7 a0 a1 a2 a3 a4 a5 a6 a7 x + Wd8 a0 a1 a2 a3 a4 a5 a6 a7 x + Wd9 a0 a1 a2 a3 a4 a5 a6 a7 x + Wd10 a0 a1 a2 a3 a4 a5 a6 a7 x + Wd11 a0 a1 a2 a3 a4 a5 a6 a7 x

set_option maxHeartbeats 64000000 in
lemma exp_trunc9 {a0 a1 a2 a3 a4 a5 a6 a7 B K : ℝ} (hB : 0 ≤ B) (_hK : 0 ≤ K)
    (ha0 : |a0| ≤ B) (ha1 : |a1| ≤ B) (ha2 : |a2| ≤ B) (ha3 : |a3| ≤ B) (ha4 : |a4| ≤ B) (ha5 : |a5| ≤ B) (ha6 : |a6| ≤ B) (ha7 : |a7| ≤ B) :
    ∃ KW : ℝ, 0 ≤ KW ∧ ∀ x L : ℝ, 0 ≤ x → x ≤ 1 → 8 * B * x ^ 2 ≤ 1 → K * x ^ 10 ≤ 1 →
      |L - (a0 * x ^ 2 + a1 * x ^ 3 + a2 * x ^ 4 + a3 * x ^ 5 + a4 * x ^ 6 + a5 * x ^ 7 + a6 * x ^ 8 + a7 * x ^ 9)| ≤ K * x ^ 10 →
      |Real.exp L - (1 + a0 * x ^ 2 + a1 * x ^ 3 + (a2 + a0 ^ 2 / 2) * x ^ 4 + (a3 + a0 * a1) * x ^ 5 + (a4 + a0 * a2 + a1 ^ 2 / 2 + a0 ^ 3 / 6) * x ^ 6 + (a5 + a0 * a3 + a1 * a2 + a0 ^ 2 * a1 / 2) * x ^ 7 + (a6 + a0 * a4 + a1 * a3 + a2 ^ 2 / 2 + a0 ^ 2 * a2 / 2 + a0 * a1 ^ 2 / 2 + a0 ^ 4 / 24) * x ^ 8 + (a7 + a0 * a5 + a1 * a4 + a2 * a3 + a0 ^ 2 * a3 / 2 + a0 * a1 * a2 + a1 ^ 3 / 6 + a0 ^ 3 * a1 / 6) * x ^ 9)| ≤ (6 * K + (8 * B) ^ 5 + KW) * x ^ 10 := by
  have hWc : Continuous (fun x : ℝ => Wpoly9 a0 a1 a2 a3 a4 a5 a6 a7 x) := by unfold Wpoly9 Wd0 Wd1 Wd2 Wd3 Wd4 Wd5 Wd6 Wd7 Wd8 Wd9 Wd10 Wd11; fun_prop
  obtain ⟨CW, hCW⟩ := (isCompact_Icc (a := (0 : ℝ)) (b := 1)).exists_bound_of_continuousOn hWc.continuousOn
  refine ⟨max CW 0, le_max_right _ _, fun x L hx0 hx1 hs1 hs2 hL => ?_⟩
  set W := Wpoly9 a0 a1 a2 a3 a4 a5 a6 a7 x with hW
  have hWb : |W| ≤ max CW 0 := by
    have := hCW x ⟨hx0, hx1⟩
    rw [Real.norm_eq_abs] at this
    exact this.trans (le_max_left _ _)
  set P := a0 * x ^ 2 + a1 * x ^ 3 + a2 * x ^ 4 + a3 * x ^ 5 + a4 * x ^ 6 + a5 * x ^ 7 + a6 * x ^ 8 + a7 * x ^ 9 with hPdef
  have hterm : ∀ (a : ℝ) (k : ℕ), |a| ≤ B → |a * x ^ (k + 2)| ≤ B * x ^ 2 := by
    intro a k ha
    have : x ^ (k + 2) ≤ x ^ 2 := pow_le_pow_of_le_one hx0 hx1 (by omega)
    rw [abs_mul, abs_of_nonneg (pow_nonneg hx0 (k + 2))]
    exact mul_le_mul ha this (pow_nonneg hx0 _) hB
  have bP : |P| ≤ 8 * B * x ^ 2 := by
    have e0 : |a0 * x ^ 2| ≤ B * x ^ 2 := hterm a0 0 ha0
    have e1 : |a1 * x ^ 3| ≤ B * x ^ 2 := hterm a1 1 ha1
    have e2 : |a2 * x ^ 4| ≤ B * x ^ 2 := hterm a2 2 ha2
    have e3 : |a3 * x ^ 5| ≤ B * x ^ 2 := hterm a3 3 ha3
    have e4 : |a4 * x ^ 6| ≤ B * x ^ 2 := hterm a4 4 ha4
    have e5 : |a5 * x ^ 7| ≤ B * x ^ 2 := hterm a5 5 ha5
    have e6 : |a6 * x ^ 8| ≤ B * x ^ 2 := hterm a6 6 ha6
    have e7 : |a7 * x ^ 9| ≤ B * x ^ 2 := hterm a7 7 ha7
    have t7 := abs_add_le (a0 * x ^ 2 + a1 * x ^ 3 + a2 * x ^ 4 + a3 * x ^ 5 + a4 * x ^ 6 + a5 * x ^ 7 + a6 * x ^ 8) (a7 * x ^ 9)
    have t6 := abs_add_le (a0 * x ^ 2 + a1 * x ^ 3 + a2 * x ^ 4 + a3 * x ^ 5 + a4 * x ^ 6 + a5 * x ^ 7) (a6 * x ^ 8)
    have t5 := abs_add_le (a0 * x ^ 2 + a1 * x ^ 3 + a2 * x ^ 4 + a3 * x ^ 5 + a4 * x ^ 6) (a5 * x ^ 7)
    have t4 := abs_add_le (a0 * x ^ 2 + a1 * x ^ 3 + a2 * x ^ 4 + a3 * x ^ 5) (a4 * x ^ 6)
    have t3 := abs_add_le (a0 * x ^ 2 + a1 * x ^ 3 + a2 * x ^ 4) (a3 * x ^ 5)
    have t2 := abs_add_le (a0 * x ^ 2 + a1 * x ^ 3) (a2 * x ^ 4)
    have t1 := abs_add_le (a0 * x ^ 2) (a1 * x ^ 3)
    rw [hPdef]; linarith
  have bP1 : |P| ≤ 1 := bP.trans hs1
  set E := L - P with hE
  have bE1 : |E| ≤ 1 := hL.trans hs2
  have hLPE : L = P + E := by rw [hE]; ring
  have hexpP : Real.exp P ≤ 3 := by
    have := Real.exp_le_exp.2 (le_of_abs_le bP1)
    have := Real.exp_one_lt_d9
    linarith
  have hexpE := Real.abs_exp_sub_one_le bE1
  have heP := Real.exp_bound bP1 (n := 5) (by norm_num)
  have hsP : ∑ m ∈ Finset.range 5, P ^ m / (m.factorial : ℝ) = 1 + P + P ^ 2 / 2 + P ^ 3 / 6 + P ^ 4 / 24 := by
    simp [Finset.sum_range_succ, Nat.factorial]
  rw [hsP] at heP
  have hP5 : |P| ^ 5 ≤ (8 * B) ^ 5 * x ^ 10 := by
    calc |P| ^ 5 ≤ (8 * B * x ^ 2) ^ 5 := pow_le_pow_left₀ (abs_nonneg _) bP 5
      _ = (8 * B) ^ 5 * x ^ 10 := by ring
  have key : Real.exp L - (1 + a0 * x ^ 2 + a1 * x ^ 3 + (a2 + a0 ^ 2 / 2) * x ^ 4 + (a3 + a0 * a1) * x ^ 5 + (a4 + a0 * a2 + a1 ^ 2 / 2 + a0 ^ 3 / 6) * x ^ 6 + (a5 + a0 * a3 + a1 * a2 + a0 ^ 2 * a1 / 2) * x ^ 7 + (a6 + a0 * a4 + a1 * a3 + a2 ^ 2 / 2 + a0 ^ 2 * a2 / 2 + a0 * a1 ^ 2 / 2 + a0 ^ 4 / 24) * x ^ 8 + (a7 + a0 * a5 + a1 * a4 + a2 * a3 + a0 ^ 2 * a3 / 2 + a0 * a1 * a2 + a1 ^ 3 / 6 + a0 ^ 3 * a1 / 6) * x ^ 9)
      = Real.exp P * (Real.exp E - 1) + (Real.exp P - (1 + P + P ^ 2 / 2 + P ^ 3 / 6 + P ^ 4 / 24)) + x ^ 10 * W := by
    rw [hLPE, Real.exp_add, hW, hPdef]; unfold Wpoly9 Wd0 Wd1 Wd2 Wd3 Wd4 Wd5 Wd6 Wd7 Wd8 Wd9 Wd10 Wd11; ring
  rw [key]
  have b1 : |Real.exp P * (Real.exp E - 1)| ≤ 3 * (2 * (K * x ^ 10)) := by
    rw [abs_mul, abs_of_pos (Real.exp_pos _)]
    exact mul_le_mul hexpP (hexpE.trans (by linarith)) (abs_nonneg _) (by norm_num)
  have b2 : |Real.exp P - (1 + P + P ^ 2 / 2 + P ^ 3 / 6 + P ^ 4 / 24)| ≤ (8 * B) ^ 5 * x ^ 10 := by
    calc _ ≤ _ := heP
      _ ≤ |P| ^ 5 * 1 := by gcongr; norm_num [Nat.factorial]
      _ ≤ _ := by rw [mul_one]; exact hP5
  have b3 : |x ^ 10 * W| ≤ x ^ 10 * max CW 0 := by
    rw [abs_mul, abs_of_nonneg (by positivity)]
    exact mul_le_mul_of_nonneg_left hWb (by positivity)
  have t1 := abs_add_le (Real.exp P * (Real.exp E - 1) + (Real.exp P - (1 + P + P ^ 2 / 2 + P ^ 3 / 6 + P ^ 4 / 24)))
    (x ^ 10 * W)
  have t2 := abs_add_le (Real.exp P * (Real.exp E - 1)) (Real.exp P - (1 + P + P ^ 2 / 2 + P ^ 3 / 6 + P ^ 4 / 24))
  have : (6 * K + (8 * B) ^ 5 + max CW 0) * x ^ 10 = 3 * (2 * (K * x ^ 10)) + (8 * B) ^ 5 * x ^ 10
    + x ^ 10 * max CW 0 := by ring
  rw [this]
  linarith only [t1, t2, b1, b2, b3]

lemma A_expand9 {α β γ : ℝ} (hα : |α| ≤ 3 / 2) (hβ : |β| ≤ 3 / 2) (hγ : |γ| ≤ 3 / 2) :
    ∃ K9 : ℝ, 0 ≤ K9 ∧ ∀ x L : ℝ, 0 ≤ x → x ≤ 1 / 40 →
    HasSum (fun k : ℕ => (-1 : ℝ) ^ (k + 2) * Y04.zn (k + 2) / ((k + 2 : ℕ) : ℝ) *
      ((α * x) ^ (k + 2) + (β * x) ^ (k + 2) - (γ * x) ^ (k + 2))) L →
    |Real.exp L - (1 + Y04.cf α β γ 0 * x ^ 2 + Y04.cf α β γ 1 * x ^ 3 + (Y04.cf α β γ 2 + Y04.cf α β γ 0 ^ 2 / 2) * x ^ 4 + (Y04.cf α β γ 3 + Y04.cf α β γ 0 * Y04.cf α β γ 1) * x ^ 5 + (Y04.cf α β γ 4 + Y04.cf α β γ 0 * Y04.cf α β γ 2 + Y04.cf α β γ 1 ^ 2 / 2 + Y04.cf α β γ 0 ^ 3 / 6) * x ^ 6 + (Y04.cf α β γ 5 + Y04.cf α β γ 0 * Y04.cf α β γ 3 + Y04.cf α β γ 1 * Y04.cf α β γ 2 + Y04.cf α β γ 0 ^ 2 * Y04.cf α β γ 1 / 2) * x ^ 7 + (Y04.cf α β γ 6 + Y04.cf α β γ 0 * Y04.cf α β γ 4 + Y04.cf α β γ 1 * Y04.cf α β γ 3 + Y04.cf α β γ 2 ^ 2 / 2 + Y04.cf α β γ 0 ^ 2 * Y04.cf α β γ 2 / 2 + Y04.cf α β γ 0 * Y04.cf α β γ 1 ^ 2 / 2 + Y04.cf α β γ 0 ^ 4 / 24) * x ^ 8 + (Y04.cf α β γ 7 + Y04.cf α β γ 0 * Y04.cf α β γ 5 + Y04.cf α β γ 1 * Y04.cf α β γ 4 + Y04.cf α β γ 2 * Y04.cf α β γ 3 + Y04.cf α β γ 0 ^ 2 * Y04.cf α β γ 3 / 2 + Y04.cf α β γ 0 * Y04.cf α β γ 1 * Y04.cf α β γ 2 + Y04.cf α β γ 1 ^ 3 / 6 + Y04.cf α β γ 0 ^ 3 * Y04.cf α β γ 1 / 6) * x ^ 9)| ≤ K9 * x ^ 10 := by
  have hB : ∀ k ≤ 7, |Y04.cf α β γ k| ≤ 116 := by
    intro k hk
    have := Y04.cf_bd hα hβ hγ k
    have : (3 / 2 : ℝ) ^ (k + 2) ≤ (3 / 2) ^ 9 := pow_le_pow_right₀ (by norm_num) (by omega)
    norm_num at this ⊢; linarith
  obtain ⟨KW, hKW, hE⟩ := exp_trunc9 (B := 116) (K := 346) (by norm_num) (by norm_num) (hB 0 (by norm_num)) (hB 1 (by norm_num)) (hB 2 (by norm_num)) (hB 3 (by norm_num)) (hB 4 (by norm_num)) (hB 5 (by norm_num)) (hB 6 (by norm_num)) (hB 7 (by norm_num))
  refine ⟨6 * 346 + (8 * 116) ^ 5 + KW, by positivity, fun x L hx0 hx hL => ?_⟩
  have h := series_trunc9 hα hβ hγ hx0 (by linarith) hL
  have hx2 : x ^ 2 ≤ 1 / 1600 := by
    have := pow_le_pow_left₀ hx0 hx 2
    norm_num at this ⊢; linarith
  have hx9 : x ^ 10 ≤ 1 / 1600 := by
    have : x ^ 10 ≤ x ^ 2 := pow_le_pow_of_le_one hx0 (by linarith) (by norm_num)
    linarith
  exact hE x L hx0 (by linarith) (by nlinarith) (by nlinarith) h


set_option maxHeartbeats 64000000 in
/-- the `x⁸` coefficient of the regularised Beta identity at `a = αx`, `b = βx` -/
theorem coef8 {α β : ℝ} (hα : 0 < α) (hβ : 0 < β) (h3 : α + β ≤ 3) :
    (1 / 5040) * ((α + β) ^ 1 * (α ^ 7 + β ^ 7)) * Y04.I 7 1 + (1 / 1440) * ((α + β) ^ 2 * (α ^ 6 + β ^ 6)) * Y04.I 6 2 + (1 / 720) * ((α + β) ^ 3 * (α ^ 5 + β ^ 5)) * Y04.I 5 3 + (1 / 576) * ((α + β) ^ 4 * (α ^ 4 + β ^ 4)) * Y04.I 4 4 + (1 / 720) * ((α + β) ^ 5 * (α ^ 3 + β ^ 3)) * Y04.I 3 5 + (1 / 1440) * ((α + β) ^ 6 * (α ^ 2 + β ^ 2)) * Y04.I 2 6 + (1 / 5040) * ((α + β) ^ 7 * (α ^ 1 + β ^ 1)) * Y04.I 1 7 + (1 / 40320) * ((α + β) ^ 8 * (α ^ 0 + β ^ 0)) * Y04.I 0 8 = (α + β) / (α * β) * 512 * (Y04.cf (α / 2) (β / 2) ((α + β) / 2) 7 + Y04.cf (α / 2) (β / 2) ((α + β) / 2) 0 * Y04.cf (α / 2) (β / 2) ((α + β) / 2) 5 + Y04.cf (α / 2) (β / 2) ((α + β) / 2) 1 * Y04.cf (α / 2) (β / 2) ((α + β) / 2) 4 + Y04.cf (α / 2) (β / 2) ((α + β) / 2) 2 * Y04.cf (α / 2) (β / 2) ((α + β) / 2) 3 + Y04.cf (α / 2) (β / 2) ((α + β) / 2) 0 ^ 2 * Y04.cf (α / 2) (β / 2) ((α + β) / 2) 3 / 2 + Y04.cf (α / 2) (β / 2) ((α + β) / 2) 0 * Y04.cf (α / 2) (β / 2) ((α + β) / 2) 1 * Y04.cf (α / 2) (β / 2) ((α + β) / 2) 2 + Y04.cf (α / 2) (β / 2) ((α + β) / 2) 1 ^ 3 / 6 + Y04.cf (α / 2) (β / 2) ((α + β) / 2) 0 ^ 3 * Y04.cf (α / 2) (β / 2) ((α + β) / 2) 1 / 6) := by
  set γ := α + β with hγ
  set c0 := Y04.cf (α / 2) (β / 2) (γ / 2) 0 with hc0
  set c1 := Y04.cf (α / 2) (β / 2) (γ / 2) 1 with hc1
  set c2 := Y04.cf (α / 2) (β / 2) (γ / 2) 2 with hc2
  set c3 := Y04.cf (α / 2) (β / 2) (γ / 2) 3 with hc3
  set c4 := Y04.cf (α / 2) (β / 2) (γ / 2) 4 with hc4
  set c5 := Y04.cf (α / 2) (β / 2) (γ / 2) 5 with hc5
  set c6 := Y04.cf (α / 2) (β / 2) (γ / 2) 6 with hc6
  set c7 := Y04.cf (α / 2) (β / 2) (γ / 2) 7 with hc7
  set KB := ∫ s in Set.Ioi (0 : ℝ), B10 s with hKB
  have hKB0 : 0 ≤ KB := setIntegral_nonneg measurableSet_Ioi (fun s hs => B10_nonneg (le_of_lt hs))
  have hα3 : α ≤ 3 := by linarith
  have hβ3 : β ≤ 3 := by linarith
  have hγ0 : 0 < γ := by rw [hγ]; linarith
  obtain ⟨K9, hK9, hA9⟩ := A_expand9 (α := α / 2) (β := β / 2) (γ := γ / 2)
      (by rw [abs_of_pos (by positivity)]; linarith) (by rw [abs_of_pos (by positivity)]; linarith)
      (by rw [abs_of_pos (by positivity)]; linarith)
  set e := γ / (α * β) with he
  have he0 : 0 ≤ e := by positivity
  set d1 := (-1) * (γ ^ 1 * (α ^ 0 + β ^ 0)) * Y04.I 0 1 - e * 4 * (c0) with hd1
  set d2 := (1) * (γ ^ 1 * (α ^ 1 + β ^ 1)) * Y04.I 1 1 + (1 / 2) * (γ ^ 2 * (α ^ 0 + β ^ 0)) * Y04.I 0 2 - e * 8 * (c1) with hd2
  set d3 := (-1 / 2) * (γ ^ 1 * (α ^ 2 + β ^ 2)) * Y04.I 2 1 + (-1 / 2) * (γ ^ 2 * (α ^ 1 + β ^ 1)) * Y04.I 1 2 + (-1 / 6) * (γ ^ 3 * (α ^ 0 + β ^ 0)) * Y04.I 0 3 - e * 16 * ((c2 + c0 ^ 2 / 2)) with hd3
  set d4 := (1 / 6) * (γ ^ 1 * (α ^ 3 + β ^ 3)) * Y04.I 3 1 + (1 / 4) * (γ ^ 2 * (α ^ 2 + β ^ 2)) * Y04.I 2 2 + (1 / 6) * (γ ^ 3 * (α ^ 1 + β ^ 1)) * Y04.I 1 3 + (1 / 24) * (γ ^ 4 * (α ^ 0 + β ^ 0)) * Y04.I 0 4 - e * 32 * ((c3 + c0 * c1)) with hd4
  set d5 := (-1 / 24) * (γ ^ 1 * (α ^ 4 + β ^ 4)) * Y04.I 4 1 + (-1 / 12) * (γ ^ 2 * (α ^ 3 + β ^ 3)) * Y04.I 3 2 + (-1 / 12) * (γ ^ 3 * (α ^ 2 + β ^ 2)) * Y04.I 2 3 + (-1 / 24) * (γ ^ 4 * (α ^ 1 + β ^ 1)) * Y04.I 1 4 + (-1 / 120) * (γ ^ 5 * (α ^ 0 + β ^ 0)) * Y04.I 0 5 - e * 64 * ((c4 + c0 * c2 + c1 ^ 2 / 2 + c0 ^ 3 / 6)) with hd5
  set d6 := (1 / 120) * (γ ^ 1 * (α ^ 5 + β ^ 5)) * Y04.I 5 1 + (1 / 48) * (γ ^ 2 * (α ^ 4 + β ^ 4)) * Y04.I 4 2 + (1 / 36) * (γ ^ 3 * (α ^ 3 + β ^ 3)) * Y04.I 3 3 + (1 / 48) * (γ ^ 4 * (α ^ 2 + β ^ 2)) * Y04.I 2 4 + (1 / 120) * (γ ^ 5 * (α ^ 1 + β ^ 1)) * Y04.I 1 5 + (1 / 720) * (γ ^ 6 * (α ^ 0 + β ^ 0)) * Y04.I 0 6 - e * 128 * ((c5 + c0 * c3 + c1 * c2 + c0 ^ 2 * c1 / 2)) with hd6
  set d7 := (-1 / 720) * (γ ^ 1 * (α ^ 6 + β ^ 6)) * Y04.I 6 1 + (-1 / 240) * (γ ^ 2 * (α ^ 5 + β ^ 5)) * Y04.I 5 2 + (-1 / 144) * (γ ^ 3 * (α ^ 4 + β ^ 4)) * Y04.I 4 3 + (-1 / 144) * (γ ^ 4 * (α ^ 3 + β ^ 3)) * Y04.I 3 4 + (-1 / 240) * (γ ^ 5 * (α ^ 2 + β ^ 2)) * Y04.I 2 5 + (-1 / 720) * (γ ^ 6 * (α ^ 1 + β ^ 1)) * Y04.I 1 6 + (-1 / 5040) * (γ ^ 7 * (α ^ 0 + β ^ 0)) * Y04.I 0 7 - e * 256 * ((c6 + c0 * c4 + c1 * c3 + c2 ^ 2 / 2 + c0 ^ 2 * c2 / 2 + c0 * c1 ^ 2 / 2 + c0 ^ 4 / 24)) with hd7
  set d8 := (1 / 5040) * (γ ^ 1 * (α ^ 7 + β ^ 7)) * Y04.I 7 1 + (1 / 1440) * (γ ^ 2 * (α ^ 6 + β ^ 6)) * Y04.I 6 2 + (1 / 720) * (γ ^ 3 * (α ^ 5 + β ^ 5)) * Y04.I 5 3 + (1 / 576) * (γ ^ 4 * (α ^ 4 + β ^ 4)) * Y04.I 4 4 + (1 / 720) * (γ ^ 5 * (α ^ 3 + β ^ 3)) * Y04.I 3 5 + (1 / 1440) * (γ ^ 6 * (α ^ 2 + β ^ 2)) * Y04.I 2 6 + (1 / 5040) * (γ ^ 7 * (α ^ 1 + β ^ 1)) * Y04.I 1 7 + (1 / 40320) * (γ ^ 8 * (α ^ 0 + β ^ 0)) * Y04.I 0 8 - e * 512 * ((c7 + c0 * c5 + c1 * c4 + c2 * c3 + c0 ^ 2 * c3 / 2 + c0 * c1 * c2 + c1 ^ 3 / 6 + c0 ^ 3 * c1 / 6)) with hd8
  have hmain : ∀ x : ℝ, 0 < x → x ≤ 1 / 80 →
      |d1 * x + d2 * x ^ 2 + d3 * x ^ 3 + d4 * x ^ 4 + d5 * x ^ 5 + d6 * x ^ 6 + d7 * x ^ 7 + d8 * x ^ 8|
        ≤ (74 * 19683 * KB + e * 1024 * K9) * x ^ 9 := by
    intro x hx0 hx1
    have hx3 : 3 * x ≤ 1 := by linarith
    have gx : γ * x ≤ 3 * x := mul_le_mul_of_nonneg_right (by rw [hγ]; exact h3) hx0.le
    have ax : α * x ≤ 3 * x := mul_le_mul_of_nonneg_right hα3 hx0.le
    have bx : β * x ≤ 3 * x := mul_le_mul_of_nonneg_right hβ3 hx0.le
    have dx1 : (α - β) * x ≤ 3 * x := mul_le_mul_of_nonneg_right (by linarith) hx0.le
    have dx2 : -3 * x ≤ (α - β) * x := mul_le_mul_of_nonneg_right (by linarith) hx0.le
    obtain ⟨_, E1⟩ := int_expand8 (u := γ * x) (v := α * x) (t := 3 * x) (by positivity)
      gx (by positivity) ax hx3
    obtain ⟨_, E2⟩ := int_expand8 (u := γ * x) (v := β * x) (t := 3 * x) (by positivity)
      gx (by positivity) bx hx3
    have hC := Z02.lhs_closed (a := α * x) (b := β * x) (by positivity) (by positivity)
    rw [show α * x + β * x = γ * x by rw [hγ]; ring] at hC
    have hη : |γ * x| < 1 / 2 := by
      rw [abs_of_pos (by positivity)]; linarith
    have hs : |(α - β) * x / 2| < 1 / 2 := by
      rw [abs_lt]; constructor <;> linarith
    obtain ⟨hpos, h⟩ := Y04.gamma_log_series (γ * x) ((α - β) * x / 2) hη hs
    rw [show γ * x / 2 + (α - β) * x / 2 = α / 2 * (2 * x) by rw [hγ]; ring,
      show γ * x / 2 - (α - β) * x / 2 = β / 2 * (2 * x) by rw [hγ]; ring] at h hpos
    rw [show γ * x = γ / 2 * (2 * x) by ring] at h hpos
    have hA := hA9 (2 * x) _ (by positivity) (by linarith) h
    rw [Real.exp_log hpos] at hA
    rw [show α / 2 * (2 * x) = α * x by ring, show β / 2 * (2 * x) = β * x by ring,
      show γ / 2 * (2 * x) = γ * x by ring] at hA
    rw [← hc0, ← hc1, ← hc2, ← hc3, ← hc4, ← hc5, ← hc6, ← hc7] at hA
    set Rg := Real.Gamma (1 + α * x) * Real.Gamma (1 + β * x) / Real.Gamma (1 + γ * x) with hRg
    set Iα := ∫ s in Set.Ioi (0 : ℝ), (Real.exp (-(γ * x) * Y04.ell s) - 1) * Real.exp (-(α * x) * s)
    set Iβ := ∫ s in Set.Ioi (0 : ℝ), (Real.exp (-(γ * x) * Y04.ell s) - 1) * Real.exp (-(β * x) * s)
    set u := γ * x with hu
    set Sα := (-1) * (u ^ 1 * (α * x) ^ 0) * Y04.I 0 1 + (1) * (u ^ 1 * (α * x) ^ 1) * Y04.I 1 1 + (-1 / 2) * (u ^ 1 * (α * x) ^ 2) * Y04.I 2 1 + (1 / 6) * (u ^ 1 * (α * x) ^ 3) * Y04.I 3 1 + (-1 / 24) * (u ^ 1 * (α * x) ^ 4) * Y04.I 4 1 + (1 / 120) * (u ^ 1 * (α * x) ^ 5) * Y04.I 5 1 + (-1 / 720) * (u ^ 1 * (α * x) ^ 6) * Y04.I 6 1 + (1 / 5040) * (u ^ 1 * (α * x) ^ 7) * Y04.I 7 1 + (1 / 2) * (u ^ 2 * (α * x) ^ 0) * Y04.I 0 2 + (-1 / 2) * (u ^ 2 * (α * x) ^ 1) * Y04.I 1 2 + (1 / 4) * (u ^ 2 * (α * x) ^ 2) * Y04.I 2 2 + (-1 / 12) * (u ^ 2 * (α * x) ^ 3) * Y04.I 3 2 + (1 / 48) * (u ^ 2 * (α * x) ^ 4) * Y04.I 4 2 + (-1 / 240) * (u ^ 2 * (α * x) ^ 5) * Y04.I 5 2 + (1 / 1440) * (u ^ 2 * (α * x) ^ 6) * Y04.I 6 2 + (-1 / 6) * (u ^ 3 * (α * x) ^ 0) * Y04.I 0 3 + (1 / 6) * (u ^ 3 * (α * x) ^ 1) * Y04.I 1 3 + (-1 / 12) * (u ^ 3 * (α * x) ^ 2) * Y04.I 2 3 + (1 / 36) * (u ^ 3 * (α * x) ^ 3) * Y04.I 3 3 + (-1 / 144) * (u ^ 3 * (α * x) ^ 4) * Y04.I 4 3 + (1 / 720) * (u ^ 3 * (α * x) ^ 5) * Y04.I 5 3 + (1 / 24) * (u ^ 4 * (α * x) ^ 0) * Y04.I 0 4 + (-1 / 24) * (u ^ 4 * (α * x) ^ 1) * Y04.I 1 4 + (1 / 48) * (u ^ 4 * (α * x) ^ 2) * Y04.I 2 4 + (-1 / 144) * (u ^ 4 * (α * x) ^ 3) * Y04.I 3 4 + (1 / 576) * (u ^ 4 * (α * x) ^ 4) * Y04.I 4 4 + (-1 / 120) * (u ^ 5 * (α * x) ^ 0) * Y04.I 0 5 + (1 / 120) * (u ^ 5 * (α * x) ^ 1) * Y04.I 1 5 + (-1 / 240) * (u ^ 5 * (α * x) ^ 2) * Y04.I 2 5 + (1 / 720) * (u ^ 5 * (α * x) ^ 3) * Y04.I 3 5 + (1 / 720) * (u ^ 6 * (α * x) ^ 0) * Y04.I 0 6 + (-1 / 720) * (u ^ 6 * (α * x) ^ 1) * Y04.I 1 6 + (1 / 1440) * (u ^ 6 * (α * x) ^ 2) * Y04.I 2 6 + (-1 / 5040) * (u ^ 7 * (α * x) ^ 0) * Y04.I 0 7 + (1 / 5040) * (u ^ 7 * (α * x) ^ 1) * Y04.I 1 7 + (1 / 40320) * (u ^ 8 * (α * x) ^ 0) * Y04.I 0 8 with hSα
    set Sβ := (-1) * (u ^ 1 * (β * x) ^ 0) * Y04.I 0 1 + (1) * (u ^ 1 * (β * x) ^ 1) * Y04.I 1 1 + (-1 / 2) * (u ^ 1 * (β * x) ^ 2) * Y04.I 2 1 + (1 / 6) * (u ^ 1 * (β * x) ^ 3) * Y04.I 3 1 + (-1 / 24) * (u ^ 1 * (β * x) ^ 4) * Y04.I 4 1 + (1 / 120) * (u ^ 1 * (β * x) ^ 5) * Y04.I 5 1 + (-1 / 720) * (u ^ 1 * (β * x) ^ 6) * Y04.I 6 1 + (1 / 5040) * (u ^ 1 * (β * x) ^ 7) * Y04.I 7 1 + (1 / 2) * (u ^ 2 * (β * x) ^ 0) * Y04.I 0 2 + (-1 / 2) * (u ^ 2 * (β * x) ^ 1) * Y04.I 1 2 + (1 / 4) * (u ^ 2 * (β * x) ^ 2) * Y04.I 2 2 + (-1 / 12) * (u ^ 2 * (β * x) ^ 3) * Y04.I 3 2 + (1 / 48) * (u ^ 2 * (β * x) ^ 4) * Y04.I 4 2 + (-1 / 240) * (u ^ 2 * (β * x) ^ 5) * Y04.I 5 2 + (1 / 1440) * (u ^ 2 * (β * x) ^ 6) * Y04.I 6 2 + (-1 / 6) * (u ^ 3 * (β * x) ^ 0) * Y04.I 0 3 + (1 / 6) * (u ^ 3 * (β * x) ^ 1) * Y04.I 1 3 + (-1 / 12) * (u ^ 3 * (β * x) ^ 2) * Y04.I 2 3 + (1 / 36) * (u ^ 3 * (β * x) ^ 3) * Y04.I 3 3 + (-1 / 144) * (u ^ 3 * (β * x) ^ 4) * Y04.I 4 3 + (1 / 720) * (u ^ 3 * (β * x) ^ 5) * Y04.I 5 3 + (1 / 24) * (u ^ 4 * (β * x) ^ 0) * Y04.I 0 4 + (-1 / 24) * (u ^ 4 * (β * x) ^ 1) * Y04.I 1 4 + (1 / 48) * (u ^ 4 * (β * x) ^ 2) * Y04.I 2 4 + (-1 / 144) * (u ^ 4 * (β * x) ^ 3) * Y04.I 3 4 + (1 / 576) * (u ^ 4 * (β * x) ^ 4) * Y04.I 4 4 + (-1 / 120) * (u ^ 5 * (β * x) ^ 0) * Y04.I 0 5 + (1 / 120) * (u ^ 5 * (β * x) ^ 1) * Y04.I 1 5 + (-1 / 240) * (u ^ 5 * (β * x) ^ 2) * Y04.I 2 5 + (1 / 720) * (u ^ 5 * (β * x) ^ 3) * Y04.I 3 5 + (1 / 720) * (u ^ 6 * (β * x) ^ 0) * Y04.I 0 6 + (-1 / 720) * (u ^ 6 * (β * x) ^ 1) * Y04.I 1 6 + (1 / 1440) * (u ^ 6 * (β * x) ^ 2) * Y04.I 2 6 + (-1 / 5040) * (u ^ 7 * (β * x) ^ 0) * Y04.I 0 7 + (1 / 5040) * (u ^ 7 * (β * x) ^ 1) * Y04.I 1 7 + (1 / 40320) * (u ^ 8 * (β * x) ^ 0) * Y04.I 0 8 with hSβ
    set Q := 1 + c0 * (2 * x) ^ 2 + c1 * (2 * x) ^ 3 + (c2 + c0 ^ 2 / 2) * (2 * x) ^ 4 + (c3 + c0 * c1) * (2 * x) ^ 5 + (c4 + c0 * c2 + c1 ^ 2 / 2 + c0 ^ 3 / 6) * (2 * x) ^ 6 + (c5 + c0 * c3 + c1 * c2 + c0 ^ 2 * c1 / 2) * (2 * x) ^ 7 + (c6 + c0 * c4 + c1 * c3 + c2 ^ 2 / 2 + c0 ^ 2 * c2 / 2 + c0 * c1 ^ 2 / 2 + c0 ^ 4 / 24) * (2 * x) ^ 8 + (c7 + c0 * c5 + c1 * c4 + c2 * c3 + c0 ^ 2 * c3 / 2 + c0 * c1 * c2 + c1 ^ 3 / 6 + c0 ^ 3 * c1 / 6) * (2 * x) ^ 9 with hQ
    have hpoly : d1 * x + d2 * x ^ 2 + d3 * x ^ 3 + d4 * x ^ 4 + d5 * x ^ 5 + d6 * x ^ 6 + d7 * x ^ 7 + d8 * x ^ 8
        = (Sα - Iα) + (Sβ - Iβ) + ((Iα + Iβ) - e * (4 * (c0) * x ^ 1 + 8 * (c1) * x ^ 2 + 16 * ((c2 + c0 ^ 2 / 2)) * x ^ 3 + 32 * ((c3 + c0 * c1)) * x ^ 4 + 64 * ((c4 + c0 * c2 + c1 ^ 2 / 2 + c0 ^ 3 / 6)) * x ^ 5 + 128 * ((c5 + c0 * c3 + c1 * c2 + c0 ^ 2 * c1 / 2)) * x ^ 6 + 256 * ((c6 + c0 * c4 + c1 * c3 + c2 ^ 2 / 2 + c0 ^ 2 * c2 / 2 + c0 * c1 ^ 2 / 2 + c0 ^ 4 / 24)) * x ^ 7 + 512 * ((c7 + c0 * c5 + c1 * c4 + c2 * c3 + c0 ^ 2 * c3 / 2 + c0 * c1 * c2 + c1 ^ 3 / 6 + c0 ^ 3 * c1 / 6)) * x ^ 8)) := by
      rw [hd1, hd2, hd3, hd4, hd5, hd6, hd7, hd8, hSα, hSβ, hu]; ring
    have hlast : (Iα + Iβ) - e * (4 * (c0) * x ^ 1 + 8 * (c1) * x ^ 2 + 16 * ((c2 + c0 ^ 2 / 2)) * x ^ 3 + 32 * ((c3 + c0 * c1)) * x ^ 4 + 64 * ((c4 + c0 * c2 + c1 ^ 2 / 2 + c0 ^ 3 / 6)) * x ^ 5 + 128 * ((c5 + c0 * c3 + c1 * c2 + c0 ^ 2 * c1 / 2)) * x ^ 6 + 256 * ((c6 + c0 * c4 + c1 * c3 + c2 ^ 2 / 2 + c0 ^ 2 * c2 / 2 + c0 * c1 ^ 2 / 2 + c0 ^ 4 / 24)) * x ^ 7 + 512 * ((c7 + c0 * c5 + c1 * c4 + c2 * c3 + c0 ^ 2 * c3 / 2 + c0 * c1 * c2 + c1 ^ 3 / 6 + c0 ^ 3 * c1 / 6)) * x ^ 8) = e / x * (Rg - Q) := by
      rw [hC, he, hQ, hu]
      field_simp
      ring
    rw [hpoly, hlast]
    have b3 : |e / x * (Rg - Q)| ≤ e * 1024 * K9 * x ^ 9 := by
      rw [abs_mul, abs_of_nonneg (by positivity)]
      calc e / x * |Rg - Q| ≤ e / x * (K9 * (2 * x) ^ 10) :=
            mul_le_mul_of_nonneg_left hA (by positivity)
        _ = e * 1024 * K9 * x ^ 9 := by field_simp; ring
    have b1 : |Sα - Iα| ≤ 37 * (3 * x) ^ 9 * KB := by rw [abs_sub_comm]; exact E1
    have b2 : |Sβ - Iβ| ≤ 37 * (3 * x) ^ 9 * KB := by rw [abs_sub_comm]; exact E2
    have t1 := abs_add_le ((Sα - Iα) + (Sβ - Iβ)) (e / x * (Rg - Q))
    have t2 := abs_add_le (Sα - Iα) (Sβ - Iβ)
    have : (74 * 19683 * KB + e * 1024 * K9) * x ^ 9
        = 37 * (3 * x) ^ 9 * KB + 37 * (3 * x) ^ 9 * KB + e * 1024 * K9 * x ^ 9 := by ring
    linarith
  obtain ⟨_, _, _, _, _, _, _, h8⟩ := octic_zero (by norm_num) (by norm_num) (by positivity) hmain
  rw [hd8] at h8
  linarith


/-! ### The four weight-nine rows `D̃_{9,3}, …, D̃_{9,0}` -/

lemma D93 : 2 * ((1 / 2) * Y04.I 7 1 + Y04.I 6 2 / 2) =
    720 * X05.z 9 + 720 * ((π ^ 2 / 6) * Q03.z7) + 504 * ((π ^ 2 / 6) ^ 2 * Z02.z5)
      + 2232 / 7 * ((π ^ 2 / 6) ^ 3 * Y04.z3) := by
  have h1 := coef8 (α := 1) (β := 1) one_pos one_pos (by norm_num)
  have h2 := coef8 (α := 1) (β := 2) one_pos two_pos (by norm_num)
  have h3 := coef8 (α := 1) (β := 1 / 3) one_pos (by norm_num) (by norm_num)
  have h4 := coef8 (α := 2) (β := 1 / 2) two_pos (by norm_num) (by norm_num)
  simp only [Y04.cf, Nat.reduceAdd, Nat.cast_ofNat, Y04.z_two, Y04.zn_three, Y04.z_four, Z02.zn_five,
    V02.z_six, Q03.zn_seven, zn9] at h1 h2 h3 h4
  linear_combination (-18000) * h1 + (12800/7) * h2 + (-13286025/11) * h3 + (1179648/385) * h4

lemma D92 : 2 * (-(11 / 16) * Y04.I 7 1 + (-(31 / 24)) * Y04.I 6 2 / 2 + (1 / 2) * Y04.I 5 3 / 6
      + Y04.I 4 4 / 24) =
    -958 * X05.z 9 - 990 * ((π ^ 2 / 6) * Q03.z7) - 3501 / 5 * ((π ^ 2 / 6) ^ 2 * Z02.z5)
      - 3093 / 7 * ((π ^ 2 / 6) ^ 3 * Y04.z3) + 4 * Y04.z3 ^ 3 := by
  have h1 := coef8 (α := 1) (β := 1) one_pos one_pos (by norm_num)
  have h2 := coef8 (α := 1) (β := 2) one_pos two_pos (by norm_num)
  have h3 := coef8 (α := 1) (β := 1 / 3) one_pos (by norm_num) (by norm_num)
  have h4 := coef8 (α := 2) (β := 1 / 2) two_pos (by norm_num) (by norm_num)
  simp only [Y04.cf, Nat.reduceAdd, Nat.cast_ofNat, Y04.z_two, Y04.zn_three, Y04.z_four, Z02.zn_five,
    V02.z_six, Q03.zn_seven, zn9] at h1 h2 h3 h4
  linear_combination (153649/6) * h1 + (-90992/35) * h2 + (188740287/110) * h3 + (-25120768/5775) * h4

lemma D91 : 2 * (4163 / 11520 * Y04.I 7 1 + 731 / 1152 * Y04.I 6 2 / 2 + (-(25 / 48)) * Y04.I 5 3 / 6
      + (-(23 / 24)) * Y04.I 4 4 / 24 + (1 / 2) * Y04.I 3 5 / 120 + Y04.I 2 6 / 720) =
    3907 / 8 * X05.z 9 + 4163 / 8 * ((π ^ 2 / 6) * Q03.z7)
      + 29741 / 80 * ((π ^ 2 / 6) ^ 2 * Z02.z5) + 130973 / 560 * ((π ^ 2 / 6) ^ 3 * Y04.z3)
      - 9 / 2 * Y04.z3 ^ 3 := by
  have h1 := coef8 (α := 1) (β := 1) one_pos one_pos (by norm_num)
  have h2 := coef8 (α := 1) (β := 2) one_pos two_pos (by norm_num)
  have h3 := coef8 (α := 1) (β := 1 / 3) one_pos (by norm_num) (by norm_num)
  have h4 := coef8 (α := 2) (β := 1 / 2) two_pos (by norm_num) (by norm_num)
  simp only [Y04.cf, Nat.reduceAdd, Nat.cast_ofNat, Y04.z_two, Y04.zn_three, Y04.z_four, Z02.zn_five,
    V02.z_six, Q03.zn_seven, zn9] at h1 h2 h3 h4
  linear_combination (-2005945/144) * h1 + (445181/315) * h2 + (-596467071/640) * h3 + (3718912/1575) * h4

lemma D90 : 2 * (-(58337 / 645120) * Y04.I 7 1 + (-(1369 / 9216)) * Y04.I 6 2 / 2
      + 721 / 3840 * Y04.I 5 3 / 6 + 121 / 384 * Y04.I 4 4 / 24
      + (-(17 / 48)) * Y04.I 3 5 / 120 + (-(5 / 8)) * Y04.I 2 6 / 720
      + (1 / 2) * Y04.I 1 7 / 5040 + Y04.I 0 8 / 40320) =
    -(22805 / 192) * X05.z 9 - 8343 / 64 * ((π ^ 2 / 6) * Q03.z7)
      - 60093 / 640 * ((π ^ 2 / 6) ^ 2 * Z02.z5) - 263713 / 4480 * ((π ^ 2 / 6) ^ 3 * Y04.z3)
      + 169 / 96 * Y04.z3 ^ 3 := by
  have h1 := coef8 (α := 1) (β := 1) one_pos one_pos (by norm_num)
  have h2 := coef8 (α := 1) (β := 2) one_pos two_pos (by norm_num)
  have h3 := coef8 (α := 1) (β := 1 / 3) one_pos (by norm_num) (by norm_num)
  have h4 := coef8 (α := 2) (β := 1 / 2) two_pos (by norm_num) (by norm_num)
  simp only [Y04.cf, Nat.reduceAdd, Nat.cast_ofNat, Y04.z_two, Y04.zn_three, Y04.z_four, Z02.zn_five,
    V02.z_six, Q03.zn_seven, zn9] at h1 h2 h3 h4
  linear_combination (920945/256) * h1 + (-1459/4) * h2 + (1353658959/5632) * h3 + (-33472/55) * h4


/-! ### Lemma `g2:lem:t:red`, `Θ₈`, `Θ₉`, the expansion, the weight claim -/

/-- Lemma `g2:lem:t:red`, coefficientwise. -/
theorem red (n : ℕ) (hn : 2 ≤ n) :
    R05.Theta n = ∑ m ∈ Finset.range (n / 2),
      Polynomial.C (Dt n m * (-4 : ℝ) ^ m / ((2 * m).factorial : ℝ)) * X ^ m := by
  ext i
  rw [R05.theta_coeff, Polynomial.finsetSum_coeff]
  simp only [Polynomial.coeff_C_mul_X_pow]
  rw [Finset.sum_ite_eq]
  by_cases hi : i < n / 2
  · rw [R05.ite_t (Finset.mem_range.2 hi)]
    unfold Dt
    have hs : Finset.Icc 1 (n - 1 - 2 * i) ⊆ Finset.Icc 1 (n - 1) := by
      intro j hj; simp only [Finset.mem_Icc] at hj ⊢; omega
    rw [← Finset.sum_subset hs (f := fun j => 2 * (Y04.I (n - 1 - j) j / (j.factorial : ℝ)) *
      (R05.rho (n - 1 - j)).coeff i)]
    · rw [Finset.mul_sum, Finset.sum_mul, Finset.sum_div]
      refine Finset.sum_congr rfl (fun j hj => ?_)
      simp only [Finset.mem_Icc] at hj
      rw [R05.rho_coeff, R05.ite_t (by omega)]
      ring
    · intro j hj hj'
      simp only [Finset.mem_Icc] at hj hj'
      rw [R05.rho_coeff, R05.ite_f (by omega)]
      ring
  · rw [R05.ite_f (by simpa using hi)]
    refine Finset.sum_eq_zero (fun j hj => ?_)
    simp only [Finset.mem_Icc] at hj
    rw [R05.rho_coeff, R05.ite_f (by omega)]
    ring

lemma theta9 (k : ℕ) :
    (-1 : ℝ) ^ k * (π / Real.sqrt 2) *
        (8 / π * ∑ j ∈ Finset.Icc 1 8,
          iteratedDeriv (8 - j) (Y04.G k) 0 / ((8 - j).factorial : ℝ) / (j.factorial : ℝ) *
            Y04.I (8 - j) j) =
      -(22805 / 192 + 3907 / 4 * (k : ℝ) ^ 2 + 1916 / 3 * (k : ℝ) ^ 4
                + 64 * (k : ℝ) ^ 6) * X05.z 9
              - (8343 / 64 + 4163 / 4 * (k : ℝ) ^ 2 + 660 * (k : ℝ) ^ 4
                + 64 * (k : ℝ) ^ 6) * ((π ^ 2 / 6) * Q03.z7)
              - (60093 / 640 + 29741 / 40 * (k : ℝ) ^ 2 + 2334 / 5 * (k : ℝ) ^ 4
                + 224 / 5 * (k : ℝ) ^ 6) * ((π ^ 2 / 6) ^ 2 * Z02.z5)
              - (263713 / 4480 + 130973 / 280 * (k : ℝ) ^ 2 + 2062 / 7 * (k : ℝ) ^ 4
                + 992 / 35 * (k : ℝ) ^ 6) * ((π ^ 2 / 6) ^ 3 * Y04.z3)
              + (169 / 96 + 9 * (k : ℝ) ^ 2 + 8 / 3 * (k : ℝ) ^ 4) * Y04.z3 ^ 3 := by
  set E := (-1 : ℝ) ^ k with hE
  set Y := (-(58337 / 645120) - 4163 / 5760 * (k : ℝ) ^ 2 - 11 / 24 * (k : ℝ) ^ 4
        - 2 / 45 * (k : ℝ) ^ 6) * Y04.I 7 1
    + (-(1369 / 9216) - 731 / 576 * (k : ℝ) ^ 2 - 31 / 36 * (k : ℝ) ^ 4 - 4 / 45 * (k : ℝ) ^ 6) / 2
      * Y04.I 6 2
    + (721 / 3840 + 25 / 24 * (k : ℝ) ^ 2 + 1 / 3 * (k : ℝ) ^ 4) / 6 * Y04.I 5 3
    + (121 / 384 + 23 / 12 * (k : ℝ) ^ 2 + 2 / 3 * (k : ℝ) ^ 4) / 24 * Y04.I 4 4
    + (-(17 / 48 + (k : ℝ) ^ 2)) / 120 * Y04.I 3 5
    + (-(2 * (k : ℝ) ^ 2 + 5 / 8)) / 720 * Y04.I 2 6
    + (1 / 2) / 5040 * Y04.I 1 7 + 1 / 40320 * Y04.I 0 8 with hY
  have hsum : ∑ j ∈ Finset.Icc 1 8,
          iteratedDeriv (8 - j) (Y04.G k) 0 / ((8 - j).factorial : ℝ) / (j.factorial : ℝ) *
            Y04.I (8 - j) j = E * Real.sqrt 2 / 4 * Y := by
    rw [show Finset.Icc 1 8 = {1, 2, 3, 4, 5, 6, 7, 8} by decide]
    rw [Finset.sum_insert (by decide), Finset.sum_insert (by decide), Finset.sum_insert (by decide),
      Finset.sum_insert (by decide), Finset.sum_insert (by decide), Finset.sum_insert (by decide),
      Finset.sum_insert (by decide), Finset.sum_singleton]
    simp only [show 8 - 1 = 7 from rfl, show 8 - 2 = 6 from rfl, show 8 - 3 = 5 from rfl,
      show 8 - 4 = 4 from rfl, show 8 - 5 = 3 from rfl, show 8 - 6 = 2 from rfl,
      show 8 - 7 = 1 from rfl, show 8 - 8 = 0 from rfl, iteratedDeriv_zero, iteratedDeriv_one]
    rw [G_d7', R02.G_d6', Q03.G_d5', V02.G_d4', Z02.G_d3', Z02.G_d2', Y04.G_d1, Y04.G_zero, hY]
    simp [Nat.factorial]
    ring
  rw [hsum]
  have hE2 : E ^ 2 = 1 := by rw [hE, ← pow_mul, mul_comm, pow_mul, neg_one_sq, one_pow]
  have hs : Real.sqrt 2 ≠ 0 := by positivity
  have hpi : π ≠ 0 := Real.pi_pos.ne'
  have e : E * (π / Real.sqrt 2) * (8 / π * (E * Real.sqrt 2 / 4 * Y)) = 2 * E ^ 2 * Y := by
    field_simp; ring
  rw [e, hE2, hY]
  linear_combination D90 - 2 * (k : ℝ) ^ 2 * D91 + 2 / 3 * (k : ℝ) ^ 4 * D92
    - 4 / 45 * (k : ℝ) ^ 6 * D93

lemma th9 (k : ℕ) : (R05.Theta 9).eval ((k : ℝ) ^ 2) =
      -(22805 / 192 + 3907 / 4 * (k : ℝ) ^ 2 + 1916 / 3 * (k : ℝ) ^ 4
                + 64 * (k : ℝ) ^ 6) * X05.z 9
              - (8343 / 64 + 4163 / 4 * (k : ℝ) ^ 2 + 660 * (k : ℝ) ^ 4
                + 64 * (k : ℝ) ^ 6) * ((π ^ 2 / 6) * Q03.z7)
              - (60093 / 640 + 29741 / 40 * (k : ℝ) ^ 2 + 2334 / 5 * (k : ℝ) ^ 4
                + 224 / 5 * (k : ℝ) ^ 6) * ((π ^ 2 / 6) ^ 2 * Z02.z5)
              - (263713 / 4480 + 130973 / 280 * (k : ℝ) ^ 2 + 2062 / 7 * (k : ℝ) ^ 4
                + 992 / 35 * (k : ℝ) ^ 6) * ((π ^ 2 / 6) ^ 3 * Y04.z3)
              + (169 / 96 + 9 * (k : ℝ) ^ 2 + 8 / 3 * (k : ℝ) ^ 4) * Y04.z3 ^ 3 := by
  rw [← R05.theta_eval]; exact theta9 k

lemma th8 (k : ℕ) : (R05.Theta 8).eval ((k : ℝ) ^ 2) =
      -(2891923 / 12288 + 495965 / 256 * (k : ℝ) ^ 2 + 60889 / 48 * (k : ℝ) ^ 4
          + 127 * (k : ℝ) ^ 6) * (π ^ 8 / 9450)
        + (141 / 16 + 50 * (k : ℝ) ^ 2 + 16 * (k : ℝ) ^ 4) * (Y04.z3 * Z02.z5)
        + (169 / 32 + 27 * (k : ℝ) ^ 2 + 8 * (k : ℝ) ^ 4) * ((π ^ 2 / 6) * Y04.z3 ^ 2) := by
  rw [← R05.theta_eval]; exact R02.theta8 k

lemma theta9_poly : S01.stmt_9 :=
  R05.poly_eq_of_sq _ _ (fun k => by rw [th9 k]; simp; ring)

lemma theta8_poly : R05.Theta 8 =
    Polynomial.C (-(2891923 / 12288) * (π ^ 8 / 9450) + 141 / 16 * (Y04.z3 * Z02.z5)
        + 169 / 32 * ((π ^ 2 / 6) * Y04.z3 ^ 2))
    + Polynomial.C (-(495965 / 256) * (π ^ 8 / 9450) + 50 * (Y04.z3 * Z02.z5)
        + 27 * ((π ^ 2 / 6) * Y04.z3 ^ 2)) * X
    + Polynomial.C (-(60889 / 48) * (π ^ 8 / 9450) + 16 * (Y04.z3 * Z02.z5)
        + 8 * ((π ^ 2 / 6) * Y04.z3 ^ 2)) * X ^ 2
    + Polynomial.C (-127 * (π ^ 8 / 9450)) * X ^ 3 :=
  R05.poly_eq_of_sq _ _ (fun k => by rw [th8 k]; simp; ring)

lemma fourier9 : S01.stmt_10 := by
  intro k
  obtain ⟨C, τ₀, hτ₀, hC⟩ := R05.proof_9.1 k 9 (by norm_num)
  refine ⟨C, τ₀, hτ₀, fun p hp hτ => ?_⟩
  convert hC p hp hτ using 4
  rw [show Finset.Icc 2 9 = {2, 3, 4, 5, 6, 7, 8, 9} by decide]
  rw [Finset.sum_insert (by decide), Finset.sum_insert (by decide), Finset.sum_insert (by decide),
    Finset.sum_insert (by decide), Finset.sum_insert (by decide), Finset.sum_insert (by decide),
    Finset.sum_insert (by decide), Finset.sum_singleton]
  rw [R05.th2, R05.th3, R05.th4, R05.th5, R05.th6, R05.th7, th8, th9]
  ring

lemma zp_zero (n : ℕ) : Q01.IsZetaPoly n 0 := ⟨[], by simp, by simp⟩

lemma cubic_coeff (a0 a1 a2 a3 : ℝ) (i : ℕ) :
    (Polynomial.C a0 + Polynomial.C a1 * X + Polynomial.C a2 * X ^ 2 + Polynomial.C a3 * X ^ 3).coeff i
      = if i = 0 then a0 else if i = 1 then a1 else if i = 2 then a2 else if i = 3 then a3 else 0 := by
  rcases i with _ | _ | _ | _ | i <;> simp [Polynomial.coeff_X, Polynomial.coeff_X_pow, Polynomial.coeff_C]

lemma zp_cubic (n : ℕ) (a0 a1 a2 a3 : ℝ) (h0 : Q01.IsZetaPoly n a0) (h1 : Q01.IsZetaPoly n a1)
    (h2 : Q01.IsZetaPoly n a2) (h3 : Q01.IsZetaPoly n a3) (i : ℕ) :
    Q01.IsZetaPoly n
      ((Polynomial.C a0 + Polynomial.C a1 * X + Polynomial.C a2 * X ^ 2 + Polynomial.C a3 * X ^ 3).coeff i) := by
  rw [cubic_coeff]
  split_ifs
  · exact h0
  · exact h1
  · exact h2
  · exact h3
  · exact zp_zero n

lemma zz2 : X05.z 2 = π ^ 2 / 6 := X05.z_two
lemma zz3 : X05.z 3 = Y04.z3 := Y04.zn_three
lemma zz4 : X05.z 4 = π ^ 4 / 90 := X05.z_four
lemma zz5 : X05.z 5 = Z02.z5 := Z02.zn_five
lemma zz6 : X05.z 6 = π ^ 6 / 945 := V02.z_six
lemma zz7 : X05.z 7 = Q03.z7 := Q03.zn_seven
lemma zz8 : X05.z 8 = π ^ 8 / 9450 := R02.z_eight

lemma zp_of {k : ℕ} {x : ℝ} (L : List (ℚ × List ℕ))
    (hL : ∀ e ∈ L, (∀ n ∈ e.2, 2 ≤ n) ∧ e.2.sum = k)
    (hx : x = (L.map fun e => (e.1 : ℝ) * (e.2.map X05.z).prod).sum) : Q01.IsZetaPoly k x :=
  ⟨L, hL, hx⟩

lemma pure_weight : S01.stmt_11 := by
  intro n h5 h9 i
  obtain ⟨_, _, _, e5, e6, e7⟩ := R05.proof_8
  interval_cases n
  · rw [e5, show Polynomial.C (-(3 / 4) * Z02.z5 - 7 / 4 * ((π ^ 2 / 6) * Y04.z3))
        + Polynomial.C (-4 * Z02.z5 - 4 * ((π ^ 2 / 6) * Y04.z3)) * X
        = Polynomial.C (-(3 / 4) * Z02.z5 - 7 / 4 * ((π ^ 2 / 6) * Y04.z3))
        + Polynomial.C (-4 * Z02.z5 - 4 * ((π ^ 2 / 6) * Y04.z3)) * X
        + Polynomial.C 0 * X ^ 2 + Polynomial.C 0 * X ^ 3 by simp]
    refine zp_cubic _ _ _ _ _ ?_ ?_ (zp_zero _) (zp_zero _) i
    · refine zp_of [(-3/4, [5]), (-7/4, [2, 3])] (by simp) ?_
      simp [zz2, zz3, zz5]; ring
    · refine zp_of [(-4, [5]), (-4, [2, 3])] (by simp) ?_
      simp [zz2, zz3, zz5]; ring
  · rw [e6, show Polynomial.C (4039 / 256 * (π ^ 6 / 945) - 7 / 8 * Y04.z3 ^ 2)
      + Polynomial.C (733 / 8 * (π ^ 6 / 945) - 2 * Y04.z3 ^ 2) * X
      + Polynomial.C (31 * (π ^ 6 / 945)) * X ^ 2
      = Polynomial.C (4039 / 256 * (π ^ 6 / 945) - 7 / 8 * Y04.z3 ^ 2)
      + Polynomial.C (733 / 8 * (π ^ 6 / 945) - 2 * Y04.z3 ^ 2) * X
      + Polynomial.C (31 * (π ^ 6 / 945)) * X ^ 2 + Polynomial.C 0 * X ^ 3 by simp]
    refine zp_cubic _ _ _ _ _ ?_ ?_ ?_ (zp_zero _) i
    · refine zp_of [(4039/256, [6]), (-7/8, [3, 3])] (by simp) ?_
      simp [zz3, zz6]; ring
    · refine zp_of [(733/8, [6]), (-2, [3, 3])] (by simp) ?_
      simp [zz3, zz6]; ring
    · refine zp_of [(31, [6])] (by simp) ?_
      simp [zz6]
  · rw [e7, show Polynomial.C (129 / 16 * Q03.z7 + 141 / 16 * ((π ^ 2 / 6) * Z02.z5)
        + 1071 / 160 * ((π ^ 2 / 6) ^ 2 * Y04.z3))
      + Polynomial.C (46 * Q03.z7 + 50 * ((π ^ 2 / 6) * Z02.z5)
        + 181 / 5 * ((π ^ 2 / 6) ^ 2 * Y04.z3)) * X
      + Polynomial.C (16 * Q03.z7 + 16 * ((π ^ 2 / 6) * Z02.z5)
        + 56 / 5 * ((π ^ 2 / 6) ^ 2 * Y04.z3)) * X ^ 2
      = Polynomial.C (129 / 16 * Q03.z7 + 141 / 16 * ((π ^ 2 / 6) * Z02.z5)
        + 1071 / 160 * ((π ^ 2 / 6) ^ 2 * Y04.z3))
      + Polynomial.C (46 * Q03.z7 + 50 * ((π ^ 2 / 6) * Z02.z5)
        + 181 / 5 * ((π ^ 2 / 6) ^ 2 * Y04.z3)) * X
      + Polynomial.C (16 * Q03.z7 + 16 * ((π ^ 2 / 6) * Z02.z5)
        + 56 / 5 * ((π ^ 2 / 6) ^ 2 * Y04.z3)) * X ^ 2 + Polynomial.C 0 * X ^ 3 by simp]
    refine zp_cubic _ _ _ _ _ ?_ ?_ ?_ (zp_zero _) i
    · refine zp_of [(129/16, [7]), (141/16, [2, 5]), (1071/160, [2, 2, 3])] (by simp) ?_
      simp [zz2, zz3, zz5, zz7]; ring
    · refine zp_of [(46, [7]), (50, [2, 5]), (181/5, [2, 2, 3])] (by simp) ?_
      simp [zz2, zz3, zz5, zz7]; ring
    · refine zp_of [(16, [7]), (16, [2, 5]), (56/5, [2, 2, 3])] (by simp) ?_
      simp [zz2, zz3, zz5, zz7]; ring
  · rw [theta8_poly]
    refine zp_cubic _ _ _ _ _ ?_ ?_ ?_ ?_ i
    · refine zp_of [(-2891923/12288, [8]), (141/16, [3, 5]), (169/32, [2, 3, 3])] (by simp) ?_
      simp [zz2, zz3, zz5, zz8]; ring
    · refine zp_of [(-495965/256, [8]), (50, [3, 5]), (27, [2, 3, 3])] (by simp) ?_
      simp [zz2, zz3, zz5, zz8]; ring
    · refine zp_of [(-60889/48, [8]), (16, [3, 5]), (8, [2, 3, 3])] (by simp) ?_
      simp [zz2, zz3, zz5, zz8]; ring
    · refine zp_of [(-127, [8])] (by simp) ?_
      simp [zz8]
  · rw [theta9_poly]
    refine zp_cubic _ _ _ _ _ ?_ ?_ ?_ ?_ i
    · refine zp_of [(-22805/192, [9]), (-8343/64, [2, 7]), (-60093/640, [2, 2, 5]),
        (-263713/4480, [2, 2, 2, 3]), (169/96, [3, 3, 3])] (by simp) ?_
      simp [zz2, zz3, zz5, zz7]; ring
    · refine zp_of [(-3907/4, [9]), (-4163/4, [2, 7]), (-29741/40, [2, 2, 5]),
        (-130973/280, [2, 2, 2, 3]), (9, [3, 3, 3])] (by simp) ?_
      simp [zz2, zz3, zz5, zz7]; ring
    · refine zp_of [(-1916/3, [9]), (-660, [2, 7]), (-2334/5, [2, 2, 5]),
        (-2062/7, [2, 2, 2, 3]), (8/3, [3, 3, 3])] (by simp) ?_
      simp [zz2, zz3, zz5, zz7]; ring
    · refine zp_of [(-64, [9]), (-64, [2, 7]), (-224/5, [2, 2, 5]),
        (-992/35, [2, 2, 2, 3])] (by simp) ?_
      simp [zz2, zz3, zz5, zz7]; ring


end S01

theorem S01.proof_1 : S01.stmt_1 := S01.evals1

theorem S01.proof_2 : S01.stmt_2 := S01.evals2

theorem S01.proof_3 : S01.stmt_3 := ⟨S01.rho7, S01.G_d7⟩

theorem S01.proof_4 : S01.stmt_4 := S01.red

theorem S01.proof_5 : S01.stmt_5 := S01.D93

theorem S01.proof_6 : S01.stmt_6 := S01.D92

theorem S01.proof_7 : S01.stmt_7 := S01.D91

theorem S01.proof_8 : S01.stmt_8 := S01.D90

theorem S01.proof_9 : S01.stmt_9 := S01.theta9_poly

theorem S01.proof_10 : S01.stmt_10 := S01.fourier9

theorem S01.proof_11 : S01.stmt_11 := S01.pure_weight
