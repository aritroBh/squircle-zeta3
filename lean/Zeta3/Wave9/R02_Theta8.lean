import Mathlib
import Zeta3.Wave5.Y04_Fourier
import Zeta3.Wave6.Z02_Theta5
import Zeta3.Wave7.V02_Theta6
import Zeta3.Wave8.Q03_Theta7
/-!
# R02 — the `τ⁸` coefficient of the Fourier coefficients `a_k(h_p)`: `Θ₈`, and `prop:degree` at `n = 8`

Target id: R02.
Paper labels: `prop:theta` and `prop:degree` (main.tex lines 1214–1225), extending
`thm:fourier` (main.tex 1204–1211, formalised in `Zeta3.Wave5.Y04_Fourier`) and its `τ⁵`, `τ⁶`,
`τ⁷` extensions (`Zeta3.Wave6.Z02_Theta5`, `Zeta3.Wave7.V02_Theta6`, `Zeta3.Wave8.Q03_Theta7`)
by one more order.
Detailed proofs: paper/ieee/proofs/sec_g2_theta.tex (the display of `Θ₈` at lines 19–20,
Lemma `g2:lem:t:red`, Table `g2:tab:D` rows `D̃_{8,0}`, …, `D̃_{8,3}`, and the check of the top
coefficient `-127ζ(8)` after the table) and paper/ieee/proofs/sec_g2_degree.tex (Lemmas
`g2:lem:d:coef`, `g2:lem:d:e`, `g2:lem:d:theta`, proof of Proposition 2, case `n = 2M`); the
Watson expansion is Lemma `g2:lem:f:watson` of sec_g2_fourier.tex.

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

main.tex 1226: "$\Theta_6,\dots,\Theta_9$ are in our notes."

Supplement sec_g2_theta.tex lines 19–20:

  \Theta_8&=-\bigl(\tfrac{2891923}{12288}+\tfrac{495965}{256}u+\tfrac{60889}{48}u^2+127u^3\bigr)\zeta(8)+\bigl(\tfrac{141}{16}+50u+16u^2\bigr)\zeta(3)\zeta(5)\\
  &\qquad+\bigl(\tfrac{169}{32}+27u+8u^2\bigr)\zeta(2)\zeta(3)^2,

Supplement, Lemma `g2:lem:f:watson`:

  For every integer $N\ge2$ there is $C$ (depending on $G$ and $N$) such that for $0<\tau<1$
  \Bigl|\,\tau\int_0^\infty\bigl[(1+\eu^{-s})^\tau-1\bigr]G(\tau s)\,ds-\sum_{n=2}^{N}\tau^n
  \sum_{j=1}^{n-1}\frac{g_{n-1-j}}{j!}\,I_{n-1-j,\,j}\Bigr|\le C\tau^{N+1}.

Supplement, Lemma `g2:lem:d:coef`:

  a_k(h_p)=A_k+\sum_{n=2}^Nc_n(k)\tau^n+O_{k,N}(\tau^{N+1}),\qquad
  c_n(k)=\frac8\pi\sum_{j=1}^{n-1}\frac{g_{k,n-1-j}}{j!}\,I_{n-1-j,\,j},

Supplement, Lemma `g2:lem:d:e`(b) and Lemma `g2:lem:d:theta`:

  $g_{k,m}=(-1)^k\frac{\sqrt2}4\,\rho_m(k^2)$ for all $k,m\ge0$.
  \Theta_n(u)=2\sum_{j=1}^{n-1}\frac{\rho_{n-1-j}(u)}{j!}\,I_{n-1-j,\,j}.

  with $e_{r,m}=2^{3/2}\,[y^m]\bigl(\gd(y)^{2r}H(y)\bigr)$,
  $\rho_m(u)=\sum_{r=0}^{\lfloor m/2\rfloor}\frac{(-4u)^r}{(2r)!}\,e_{r,m}$ and the values
  "$m=0$: $1,\frac12,-\frac58,-\frac{17}{48},\frac{121}{384}$; $m=1$ ($r\ge2$):
  $1,\frac12,-\frac{23}{24},-\frac{25}{48},\frac{731}{1152}$; $m=2$ ($r\ge4$):
  $1,\frac12,-\frac{31}{24},-\frac{11}{16},\frac{409}{384}$", and (Lemma `g2:lem:d:e`(a))
  "$e_{r,2r}=1$". So `e_{1,6} = 731/1152`, `e_{2,6} = -31/24`, `e_{3,6} = 1`.
  The supplement lists `e_{0,m}` only for `m ≤ 4`; `e_{0,5} = 721/3840` (as in Q03) and
  `e_{0,6} = -1369/9216` were computed here from the definition (SymPy series of
  `2^{3/2} H(y) = e^{-y}((1+e^{-2y})/2)^{-3/2}` and of `gd`; the listed values were reproduced).
  So `ρ₆(u) = e_{0,6} - 2u e_{1,6} + (2/3)u² e_{2,6} - (4/45)u³ e_{3,6}
    = -1369/9216 - (731/576)u - (31/36)u² - (4/45)u³`.

Supplement, Lemma `g2:lem:t:red` and Table `g2:tab:D` (products `ζ(8), ζ(3)ζ(5), ζ(2)ζ(3)²`):

  $\wD_{n,m}=2\sum_{j=1}^{n-1-2m}\dfrac{e_{m,\,n-1-j}}{j!}\,I_{n-1-j,\,j}$
  $\Theta_n(u)=\sum_{m=0}^{\lfloor n/2\rfloor-1}\wD_{n,m}\,\dfrac{(-4u)^m}{(2m)!}$.
  $\wD_{8,0}$ & $-\tfrac{2891923}{12288},\ \tfrac{141}{16},\ \tfrac{169}{32}$ \\
  $\wD_{8,1}$ & $\tfrac{495965}{512},\ -25,\ -\tfrac{27}2$ \\
  $\wD_{8,2}$ & $-\tfrac{60889}{32},\ 24,\ 12$ \\
  $\wD_{8,3}$ & $\tfrac{5715}4,\ 0,\ 0$ \\

  "the top coefficients $\frac{(-4)^{M-1}}{(2M-2)!}\wD_{n,M-1}$ at $n=6,8$ are [...]
  $\frac{(-4)^3}{6!}\cdot\frac{5715}4\zeta(8)=-127\zeta(8)$, as in Proposition~2."

Supplement, Lemma `g2:lem:f:atoms`: `I_{m,1} = m! η(m+2)`; at `m = 6`,
`I_{6,1} = 6! (1 - 2⁻⁷) ζ(8) = (5715/8) ζ(8)`.

Reused (not redefined): from `Zeta3.Wave5.Y04_Fourier`: `Y04.tau` (`τ = 1 - 1/p`), `Y04.h`
(the support function `h_p`), `Y04.a` (`a_k(f) = (4/π)∫_0^{π/2} f cos 4kθ`), `Y04.G` (the layer
function `G_k(x) = cos(4kθ_x) cos²θ_x sin θ_x`, `θ_x = arctan e^{-x}`), `Y04.I`
(`I_{m,j} = ∫_0^∞ s^m log^j(1+e^{-s}) ds`), `Y04.z3` (`ζ(3) = (riemannZeta 3).re`); from
`Zeta3.Wave6.Z02_Theta5`: `Z02.z5` (`ζ(5) = (riemannZeta 5).re`); from `Zeta3.Wave8.Q03_Theta7`:
`Q03.z7` (`ζ(7) = (riemannZeta 7).re`). All are statement-block definitions of those files.
Y04, Z02, V02, Q03 prove `g_{k,m}` for `m ≤ 5`, Watson up to `N = 7`, the rows of weight
`≤ 7` and the expansion through `τ⁷` (`Q03.stmt_9`); none is restated here.

Conventions.
* `ζ(3) = Y04.z3`, `ζ(5) = Z02.z5`, `ζ(7) = Q03.z7`, the real parts of Mathlib's complex zeta
  values; by `zeta_nat_eq_tsum_of_gt_one` these are the real numbers `∑ 1/n^3`, `∑ 1/n^5`,
  `∑ 1/n^7`. The even values are written out: `ζ(2) = π²/6`, `ζ(4) = π⁴/90`, `ζ(6) = π⁶/945`
  (as in Y04, V02), and `ζ(8) = π⁸/9450`, which `stmt_1` proves to be Mathlib's
  `riemannZeta 8`. So `ζ(3)ζ(5)` is `Y04.z3 * Z02.z5` and `ζ(2)ζ(3)²` is `(π²/6) * Y04.z3²`.
* `g_{k,m} = G_k^{(m)}(0)/m!` is written `iteratedDeriv m (Y04.G k) 0 / m!`; `G_k` is `C^∞`
  on `ℝ` (proved in Y04), so these are genuine derivatives.
* All half-line integrals `I_{m,j}`, `j ≥ 1`, have integrable integrands (proved in Y04), and
  `stmt_3` asserts the Watson integrand integrable, so no statement holds by the
  "non-integrable ↦ 0" convention.
* As in Y04/Z02/V02/Q03, the Watson lemma (`stmt_3`) asks `G` to be `C^∞` on all of `ℝ` and
  bounded on `[0,∞)` (the supplement: `C^∞` on `[0,∞)`); this covers `G = G_k`.
* "The coefficient of `τ⁸`" in `prop:degree` is the coefficient of an asymptotic expansion
  (the supplement's "Setting"). `stmt_9` states it that way: a polynomial `Θ₈ ∈ ℝ[X]` such
  that `(-1)^k (√2/π)[-4/(16k²-1) + Σ_{n=2}^{7} Θ_n(k²) τ^n + Θ₈(k²) τ⁸]` (with the proved
  `Θ₂, …, Θ₇` of Y04/Z02/V02/Q03 written out) approximates `a_k(h_p)` to `O_k(τ⁹)` for every
  `k`. Such coefficients are unique and a polynomial is determined by its values at `k²`, so
  this pins down the paper's `Θ₈`. "Degree exactly 3" is `Polynomial.degree Θ₈ = 3`, and the
  leading coefficient is the paper's `(-1)^{m-1}(2^{2m-1}-1)ζ(2m)` at `m = 4`, i.e. `-127ζ(8)`.
* The `O(τ⁹)` claims are explicit: `∃ C τ₀, 0 < τ₀ ∧ ∀ p > 1, τ ≤ τ₀ → |error| ≤ C τ⁹`, the
  same form as `Y04.stmt_8`, `Z02.stmt_8`, `V02.stmt_8`, `Q03.stmt_9`.
* Numerical check before freezing (mpmath, 40 digits; SymPy for `e_{r,m}`): `ζ(8) = π⁸/9450`,
  `I_{6,1} = 5715ζ(8)/8`, the four rows `D̃_{8,0}, …, D̃_{8,3}` (`stmt_4`–`stmt_7`) agree to
  `≤ 3·10⁻³⁴`; `g_{k,6}` (`stmt_2`) agrees exactly with the SymPy Taylor coefficient of `G_k` for
  `k = 0, 1, 2, 3`; the `τ⁸` coefficient `2 Σ_j ρ_{7-j}(k²) I_{7-j,j}/j!` minus the closed form
  `Θ₈(k²)` (`stmt_8`) is `≤ 3·10⁻³²` in absolute value for `k = 0, 1, 2, 3` (values
  `-212.768…, -3393.79…, -35404.7…, -209555.7…`). Direct quadrature of `a_k(h_p)` itself was not
  redone at this order.

NOT COVERED: `prop:degree` for `n ≠ 8` beyond what earlier files prove, and its clause
"coefficients are alternating multiple zeta values of weight `n`" (no AMZV/iterated-integral
library is set up; for `n = 8` the explicit closed form in `ζ(8), ζ(3)ζ(5), ζ(2)ζ(3)²` is stated
instead, which is stronger in content but not the same sentence); `Θ₉` and `prop:theta` for
`n = 9`; the iterated-integral form of `I_{m,j}`, the shuffle/stuffle/duality relations and the
exact certificates (the four `n = 8` identities `D̃_{8,0}`, …, `D̃_{8,3}` are stated directly);
the general formula for `e_{r,m}`, `ρ_m` (only `g_{k,6}` is stated); Cauchy's formula
`P = ∫ h_p` and the `k = 0` comparison with `thm:pw`.
-/

-- BEGIN STATEMENTS
namespace R02

open Real MeasureTheory Set Polynomial

/-- The value `ζ(8) = π⁸/9450` (Euler; Mathlib's `riemannZeta` at `8`), used to write `ζ(8)`
as `π⁸/9450` in the statements below. -/
def stmt_1 : Prop := riemannZeta 8 = ((π ^ 8 / 9450 : ℝ) : ℂ)

/-- Lemma `g2:lem:d:e`(b) at `m = 6` (with
`ρ₆(u) = -1369/9216 - (731/576)u - (31/36)u² - (4/45)u³`): the next Taylor coefficient of `G_k`,
`g_{k,6} = G_k^{(6)}(0)/6! = (-1)^k (√2/4)(-1369/9216 - (731/576)k² - (31/36)k⁴ - (4/45)k⁶)`. -/
def stmt_2 : Prop :=
  ∀ k : ℕ, iteratedDeriv 6 (Y04.G k) 0 / 720 =
    (-1 : ℝ) ^ k * Real.sqrt 2 / 4 *
      (-(1369 / 9216) - 731 / 576 * (k : ℝ) ^ 2 - 31 / 36 * (k : ℝ) ^ 4 - 4 / 45 * (k : ℝ) ^ 6)

/-- Lemma `g2:lem:f:watson` with `N = 8`: for `G` smooth and bounded on `[0,∞)`, with
`g_m = G^{(m)}(0)/m!`, there is `C` with, for `0 < τ < 1`, the integrand integrable and
`|τ ∫_0^∞ [(1+e^{-s})^τ - 1] G(τ s) ds - Σ_{n=2}^{8} τ^n Σ_{j=1}^{n-1} g_{n-1-j}/j! I_{n-1-j,j}|
  ≤ C τ⁹`. -/
def stmt_3 : Prop :=
  ∀ Gf : ℝ → ℝ, ContDiff ℝ (⊤ : ℕ∞) Gf → (∃ M : ℝ, ∀ x : ℝ, 0 ≤ x → |Gf x| ≤ M) →
    ∃ C : ℝ, ∀ τ : ℝ, 0 < τ → τ < 1 →
      IntegrableOn (fun s : ℝ => ((1 + Real.exp (-s)) ^ τ - 1) * Gf (τ * s)) (Ioi 0) ∧
      |τ * (∫ s in Ioi (0 : ℝ), ((1 + Real.exp (-s)) ^ τ - 1) * Gf (τ * s)) -
          ∑ n ∈ Finset.Icc 2 8, τ ^ n * ∑ j ∈ Finset.Icc 1 (n - 1),
            iteratedDeriv (n - 1 - j) Gf 0 / ((n - 1 - j).factorial : ℝ) /
              (j.factorial : ℝ) * Y04.I (n - 1 - j) j| ≤ C * τ ^ 9

/-- Table `g2:tab:D`, row `D̃_{8,3}` (with `e_{3,6} = 1`): `D̃_{8,3} = 2 I_{6,1} = 5715 ζ(8)/4`,
i.e. Lemma `g2:lem:f:atoms` at `m = 6`: `I_{6,1} = 6! η(8) = 5715 ζ(8)/8`, `ζ(8) = π⁸/9450`. -/
def stmt_4 : Prop := 2 * Y04.I 6 1 = 5715 / 4 * (π ^ 8 / 9450)

/-- Table `g2:tab:D`, row `D̃_{8,2}` (with `e_{2,4..6} = 1, 1/2, -31/24`):
`D̃_{8,2} = 2(e_{2,6} I_{6,1} + e_{2,5} I_{5,2}/2! + e_{2,4} I_{4,3}/3!)
  = -(60889/32) ζ(8) + 24 ζ(3) ζ(5) + 12 ζ(2) ζ(3)²`. -/
def stmt_5 : Prop :=
  2 * (-(31 / 24) * Y04.I 6 1 + (1 / 2) * Y04.I 5 2 / 2 + Y04.I 4 3 / 6) =
    -(60889 / 32) * (π ^ 8 / 9450) + 24 * (Y04.z3 * Z02.z5)
      + 12 * ((π ^ 2 / 6) * Y04.z3 ^ 2)

/-- Table `g2:tab:D`, row `D̃_{8,1}`: with `e_{1,2..6} = 1, 1/2, -23/24, -25/48, 731/1152`,
`D̃_{8,1} = 2(e_{1,6} I_{6,1} + e_{1,5} I_{5,2}/2! + e_{1,4} I_{4,3}/3! + e_{1,3} I_{3,4}/4!
  + e_{1,2} I_{2,5}/5!) = (495965/512) ζ(8) - 25 ζ(3) ζ(5) - (27/2) ζ(2) ζ(3)²`. -/
def stmt_6 : Prop :=
  2 * (731 / 1152 * Y04.I 6 1 + (-(25 / 48)) * Y04.I 5 2 / 2 + (-(23 / 24)) * Y04.I 4 3 / 6
      + (1 / 2) * Y04.I 3 4 / 24 + Y04.I 2 5 / 120) =
    495965 / 512 * (π ^ 8 / 9450) - 25 * (Y04.z3 * Z02.z5)
      - 27 / 2 * ((π ^ 2 / 6) * Y04.z3 ^ 2)

/-- Table `g2:tab:D`, row `D̃_{8,0}`: with
`e_{0,0..6} = 1, 1/2, -5/8, -17/48, 121/384, 721/3840, -1369/9216`,
`D̃_{8,0} = 2(e_{0,6} I_{6,1} + e_{0,5} I_{5,2}/2! + e_{0,4} I_{4,3}/3! + e_{0,3} I_{3,4}/4!
  + e_{0,2} I_{2,5}/5! + e_{0,1} I_{1,6}/6! + e_{0,0} I_{0,7}/7!)
  = -(2891923/12288) ζ(8) + (141/16) ζ(3) ζ(5) + (169/32) ζ(2) ζ(3)²`. -/
def stmt_7 : Prop :=
  2 * (-(1369 / 9216) * Y04.I 6 1 + 721 / 3840 * Y04.I 5 2 / 2 + 121 / 384 * Y04.I 4 3 / 6
      + (-(17 / 48)) * Y04.I 3 4 / 24 + (-(5 / 8)) * Y04.I 2 5 / 120
      + (1 / 2) * Y04.I 1 6 / 720 + Y04.I 0 7 / 5040) =
    -(2891923 / 12288) * (π ^ 8 / 9450) + 141 / 16 * (Y04.z3 * Z02.z5)
      + 169 / 32 * ((π ^ 2 / 6) * Y04.z3 ^ 2)

/-- Lemmas `g2:lem:d:coef`, `g2:lem:d:theta` and `g2:lem:t:red` at `n = 8`, together with the
closed form of `Θ₈` (sec_g2_theta.tex lines 19–20): the `τ⁸` coefficient
`c_8(k) = (8/π) Σ_{j=1}^{7} g_{k,7-j}/j! I_{7-j,j}` satisfies
`(-1)^k (π/√2) c_8(k) = Θ₈(k²) = -(2891923/12288 + (495965/256)k² + (60889/48)k⁴ + 127k⁶) ζ(8)
  + (141/16 + 50k² + 16k⁴) ζ(3) ζ(5) + (169/32 + 27k² + 8k⁴) ζ(2) ζ(3)²`. -/
def stmt_8 : Prop :=
  ∀ k : ℕ,
    (-1 : ℝ) ^ k * (π / Real.sqrt 2) *
        (8 / π * ∑ j ∈ Finset.Icc 1 7,
          iteratedDeriv (7 - j) (Y04.G k) 0 / ((7 - j).factorial : ℝ) / (j.factorial : ℝ) *
            Y04.I (7 - j) j) =
      -(2891923 / 12288 + 495965 / 256 * (k : ℝ) ^ 2 + 60889 / 48 * (k : ℝ) ^ 4
          + 127 * (k : ℝ) ^ 6) * (π ^ 8 / 9450)
        + (141 / 16 + 50 * (k : ℝ) ^ 2 + 16 * (k : ℝ) ^ 4) * (Y04.z3 * Z02.z5)
        + (169 / 32 + 27 * (k : ℝ) ^ 2 + 8 * (k : ℝ) ^ 4) * ((π ^ 2 / 6) * Y04.z3 ^ 2)

/-- `prop:degree` at `n = 8 = 2m`, `m = 4`: there is a polynomial `Θ₈` of degree exactly
`⌊8/2⌋ - 1 = 3`, with leading coefficient `(-1)^{m-1}(2^{2m-1}-1) ζ(2m) = -127 ζ(8)`, whose value
at `u = k²` is the `τ⁸` coefficient of `(-1)^k (π/√2) a_k(h_p)` in the asymptotic sense
(the lower coefficients `Θ₂, …, Θ₇` written out as proved in Y04, Z02, V02 and Q03): for every
`k` there are `C_k`, `τ₀ > 0` with
`|a_k(h_p) - (-1)^k (√2/π)[-4/(16k²-1) + ζ(2)τ² + ζ(3)τ³ - (33/16 + 7k²)ζ(4)τ⁴ + Θ₅(k²)τ⁵
   + Θ₆(k²)τ⁶ + Θ₇(k²)τ⁷ + Θ₈(k²)τ⁸]| ≤ C_k τ⁹` for `p > 1`, `τ ≤ τ₀`. -/
def stmt_9 : Prop :=
  ∃ Θ₈ : ℝ[X],
    Θ₈.degree = ((8 / 2 - 1 : ℕ) : WithBot ℕ) ∧
    Θ₈.leadingCoeff = (-1 : ℝ) ^ (4 - 1) * (2 ^ (2 * 4 - 1) - 1) * (π ^ 8 / 9450) ∧
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
            + Θ₈.eval ((k : ℝ) ^ 2) * Y04.tau p ^ 8)|
        ≤ C * Y04.tau p ^ 9

/-- `thm:fourier` extended by four orders with `prop:theta`'s `Θ₅` and the supplement's `Θ₆`,
`Θ₇`, `Θ₈`: for every `k ≥ 0` there are `C_k` and `τ₀ > 0` such that for every `p > 1` with
`τ = 1 - 1/p ≤ τ₀`, with `u = k²`,
`|a_k(h_p) - (-1)^k (√2/π)[-4/(16k²-1) + ζ(2)τ² + ζ(3)τ³ - (33/16 + 7k²)ζ(4)τ⁴ + Θ₅(u)τ⁵
   + Θ₆(u)τ⁶ + Θ₇(u)τ⁷ + Θ₈(u)τ⁸]| ≤ C_k τ⁹`, where
`Θ₅(u) = -(3/4 + 4u) ζ(5) - (7/4 + 4u) ζ(2) ζ(3)`,
`Θ₆(u) = (4039/256 + (733/8)u + 31u²) ζ(6) - (7/8 + 2u) ζ(3)²`,
`Θ₇(u) = (129/16 + 46u + 16u²) ζ(7) + (141/16 + 50u + 16u²) ζ(2) ζ(5)
  + (1071/160 + (181/5)u + (56/5)u²) ζ(2)² ζ(3)` and
`Θ₈(u) = -(2891923/12288 + (495965/256)u + (60889/48)u² + 127u³) ζ(8)
  + (141/16 + 50u + 16u²) ζ(3) ζ(5) + (169/32 + 27u + 8u²) ζ(2) ζ(3)²`. -/
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
                  * ((π ^ 2 / 6) * Y04.z3 ^ 2)) * Y04.tau p ^ 8)|
      ≤ C * Y04.tau p ^ 9

end R02
-- END STATEMENTS

open Real MeasureTheory Set

namespace R02


lemma bernoulli'_eight : bernoulli' 8 = -1 / 30 := by
  rw [bernoulli'_def]
  have h7 : bernoulli' 7 = 0 := bernoulli'_eq_zero_of_odd (by decide) (by norm_num)
  have h5 : bernoulli' 5 = 0 := bernoulli'_eq_zero_of_odd (by decide) (by norm_num)
  have h3 : bernoulli' 3 = 0 := bernoulli'_eq_zero_of_odd (by decide) (by norm_num)
  norm_num [Finset.sum_range_succ, Finset.sum_range_zero, h7, h5, h3, V02.bernoulli'_six,
    bernoulli'_four, bernoulli'_two, bernoulli'_one, bernoulli'_zero, Nat.choose]

lemma zeta8 : riemannZeta 8 = ((π ^ 8 / 9450 : ℝ) : ℂ) := by
  have h := riemannZeta_two_mul_nat (k := 4) (by norm_num)
  rw [show (8 : ℂ) = 2 * ((4 : ℕ) : ℂ) by norm_num, h,
    bernoulli_eq_bernoulli'_of_ne_one (by decide), bernoulli'_eight]
  push_cast
  norm_num [Nat.factorial]
  ring

lemma z_eight : Y04.zn 8 = π ^ 8 / 9450 := by
  unfold Y04.zn
  rw [show ((8 : ℕ) : ℂ) = 8 by norm_num, zeta8, Complex.ofReal_re]

lemma I61 : Y04.I 6 1 = 5715 / 8 * (π ^ 8 / 9450) := by
  rw [Y04.I_m1, Y04.eta_val 8 (by norm_num), z_eight]; norm_num [Nat.factorial]; ring

set_option maxRecDepth 40000 in
set_option maxHeartbeats 40000000 in
lemma G_d6 (k : ℕ) : iteratedDeriv 6 (Y04.G k) 0 / 720 =
    (-1 : ℝ) ^ k * Real.sqrt 2 / 4 *
      (-(1369 / 9216) - 731 / 576 * (k : ℝ) ^ 2 - 31 / 36 * (k : ℝ) ^ 4 - 4 / 45 * (k : ℝ) ^ 6) := by
  rw [V02.iter_G]
  simp only [Y04.thetaX_zero, Function.iterate_succ_apply', Function.iterate_zero_apply, V02.step, V02.dterm,
    V02.shift, List.flatMap_cons, List.flatMap_nil, List.map_cons, List.map_nil, List.append_nil,
    List.cons_append, List.nil_append, ite_true, ite_false, Bool.false_eq_true]
  simp only [V02.ev, V02.tm, List.map_cons, List.map_nil, List.sum_cons, List.sum_nil, ite_true, ite_false,
    Bool.false_eq_true, Y04.cos_k, Y04.sin_k, Real.cos_pi_div_four, Real.sin_pi_div_four]
  norm_num
  have : Real.sqrt 2 ^ 2 = 2 := Real.sq_sqrt (by norm_num)
  have h15 : Real.sqrt 2 ^ 15 = 128 * Real.sqrt 2 := by
    rw [show Real.sqrt 2 ^ 15 = (Real.sqrt 2 ^ 2) ^ 7 * Real.sqrt 2 by ring, this]; ring
  ring_nf
  rw [h15]; ring


/-- Taylor's formula at `0` to order 6 with a remainder `O(x^7)` on `[0,∞)`. -/
lemma taylor6 {Gf : ℝ → ℝ} (hG : ContDiff ℝ (⊤ : ℕ∞) Gf) {M : ℝ}
    (hM : ∀ x : ℝ, 0 ≤ x → |Gf x| ≤ M) :
    ∃ C : ℝ, 0 ≤ C ∧ ∀ x : ℝ, 0 ≤ x →
      |Gf x - (iteratedDeriv 0 Gf 0 / (Nat.factorial 0 : ℝ) + iteratedDeriv 1 Gf 0 / (Nat.factorial 1 : ℝ) * x + iteratedDeriv 2 Gf 0 / (Nat.factorial 2 : ℝ) * x ^ 2 + iteratedDeriv 3 Gf 0 / (Nat.factorial 3 : ℝ) * x ^ 3 + iteratedDeriv 4 Gf 0 / (Nat.factorial 4 : ℝ) * x ^ 4 + iteratedDeriv 5 Gf 0 / (Nat.factorial 5 : ℝ) * x ^ 5 + iteratedDeriv 6 Gf 0 / (Nat.factorial 6 : ℝ) * x ^ 6)| ≤ C * x ^ 7 := by
  set g0 := iteratedDeriv 0 Gf 0 / (Nat.factorial 0 : ℝ)
  set g1 := iteratedDeriv 1 Gf 0 / (Nat.factorial 1 : ℝ)
  set g2 := iteratedDeriv 2 Gf 0 / (Nat.factorial 2 : ℝ)
  set g3 := iteratedDeriv 3 Gf 0 / (Nat.factorial 3 : ℝ)
  set g4 := iteratedDeriv 4 Gf 0 / (Nat.factorial 4 : ℝ)
  set g5 := iteratedDeriv 5 Gf 0 / (Nat.factorial 5 : ℝ)
  set g6 := iteratedDeriv 6 Gf 0 / (Nat.factorial 6 : ℝ)
  have hc7 : Continuous (iteratedDeriv 7 Gf) :=
    hG.continuous_iteratedDeriv 7 (WithTop.coe_le_coe.2 le_top)
  obtain ⟨K7, hK7⟩ := (isCompact_Icc (a := (0 : ℝ)) (b := 1)).exists_bound_of_continuousOn
    hc7.continuousOn
  have hM0 : 0 ≤ M := (abs_nonneg _).trans (hM 0 le_rfl)
  refine ⟨max (|K7| / 5040) (M + |g0| + |g1| + |g2| + |g3| + |g4| + |g5| + |g6|), le_max_of_le_left (by positivity), ?_⟩
  intro x hx
  rcases eq_or_lt_of_le hx with h0 | hpos
  · subst h0; simp [g0]
  rcases le_or_gt x 1 with hx1 | hx1
  · have hcd : ContDiffOn ℝ (6 + 1) Gf (uIcc 0 x) := (hG.of_le (WithTop.coe_le_coe.2 le_top)).contDiffOn
    obtain ⟨x', hx', heq⟩ := taylor_mean_remainder_lagrange_iteratedDeriv (n := 6) hpos.ne hcd
    have hu : UniqueDiffOn ℝ (uIcc 0 x) := uniqueDiffOn_uIcc hpos.ne
    have hmem : (0 : ℝ) ∈ uIcc 0 x := left_mem_uIcc
    have hT : taylorWithinEval Gf 6 (uIcc 0 x) 0 x
        = g0 + g1 * x + g2 * x ^ 2 + g3 * x ^ 3 + g4 * x ^ 4 + g5 * x ^ 5 + g6 * x ^ 6 := by
      rw [taylor_within_apply]
      simp only [Finset.sum_range_succ, Finset.sum_range_zero, zero_add, sub_zero, smul_eq_mul]
      rw [iteratedDerivWithin_eq_iteratedDeriv hu (hG.contDiffAt.of_le (WithTop.coe_le_coe.2 le_top)) hmem,
        iteratedDerivWithin_eq_iteratedDeriv hu (hG.contDiffAt.of_le (WithTop.coe_le_coe.2 le_top)) hmem,
        iteratedDerivWithin_eq_iteratedDeriv hu (hG.contDiffAt.of_le (WithTop.coe_le_coe.2 le_top)) hmem,
        iteratedDerivWithin_eq_iteratedDeriv hu (hG.contDiffAt.of_le (WithTop.coe_le_coe.2 le_top)) hmem,
        iteratedDerivWithin_eq_iteratedDeriv hu (hG.contDiffAt.of_le (WithTop.coe_le_coe.2 le_top)) hmem,
        iteratedDerivWithin_eq_iteratedDeriv hu (hG.contDiffAt.of_le (WithTop.coe_le_coe.2 le_top)) hmem,
        iteratedDerivWithin_eq_iteratedDeriv hu (hG.contDiffAt.of_le (WithTop.coe_le_coe.2 le_top)) hmem]
      simp only [g0, g1, g2, g3, g4, g5, g6]; norm_num; ring
    rw [hT] at heq
    rw [heq]
    rw [uIoo_of_le hx] at hx'
    have hb := hK7 x' ⟨hx'.1.le, by linarith [hx'.2]⟩
    rw [Real.norm_eq_abs] at hb
    rw [sub_zero, abs_div, abs_mul, abs_pow, abs_of_pos hpos]
    norm_num [Nat.factorial]
    have : |iteratedDeriv 7 Gf x'| * x ^ 7 / 5040 ≤ |K7| / 5040 * x ^ 7 := by
      have := hb.trans (le_abs_self K7)
      have h3 : 0 ≤ x ^ 7 := by positivity
      nlinarith
    calc |iteratedDeriv 7 Gf x'| * x ^ 7 / 5040 ≤ |K7| / 5040 * x ^ 7 := this
      _ ≤ _ := by gcongr; exact le_max_left _ _
  · have h1 : 1 ≤ x := hx1.le
    have hx2' : 1 ≤ x ^ 2 := one_le_pow₀ h1
    have hx3' : 1 ≤ x ^ 3 := one_le_pow₀ h1
    have hx4' : 1 ≤ x ^ 4 := one_le_pow₀ h1
    have hx5' : 1 ≤ x ^ 5 := one_le_pow₀ h1
    have hx6' : 1 ≤ x ^ 6 := one_le_pow₀ h1
    have hx17 : x ^ 1 ≤ x ^ 7 := pow_le_pow_right₀ h1 (by norm_num)
    have hx27 : x ^ 2 ≤ x ^ 7 := pow_le_pow_right₀ h1 (by norm_num)
    have hx37 : x ^ 3 ≤ x ^ 7 := pow_le_pow_right₀ h1 (by norm_num)
    have hx47 : x ^ 4 ≤ x ^ 7 := pow_le_pow_right₀ h1 (by norm_num)
    have hx57 : x ^ 5 ≤ x ^ 7 := pow_le_pow_right₀ h1 (by norm_num)
    have hx67 : x ^ 6 ≤ x ^ 7 := pow_le_pow_right₀ h1 (by norm_num)
    have h17 : 1 ≤ x ^ 7 := one_le_pow₀ h1
    calc _ ≤ |Gf x| + (|g0| + |g1| * x + |g2| * x ^ 2 + |g3| * x ^ 3 + |g4| * x ^ 4 + |g5| * x ^ 5 + |g6| * x ^ 6) := by
          have := abs_sub (Gf x) (g0 + g1 * x + g2 * x ^ 2 + g3 * x ^ 3 + g4 * x ^ 4 + g5 * x ^ 5 + g6 * x ^ 6)
          have a6 := abs_add_le (g0 + g1 * x + g2 * x ^ 2 + g3 * x ^ 3 + g4 * x ^ 4 + g5 * x ^ 5) (g6 * x ^ 6)
          rw [abs_mul, abs_of_nonneg (by positivity : (0:ℝ) ≤ x ^ 6)] at a6
          have a5 := abs_add_le (g0 + g1 * x + g2 * x ^ 2 + g3 * x ^ 3 + g4 * x ^ 4) (g5 * x ^ 5)
          rw [abs_mul, abs_of_nonneg (by positivity : (0:ℝ) ≤ x ^ 5)] at a5
          have a4 := abs_add_le (g0 + g1 * x + g2 * x ^ 2 + g3 * x ^ 3) (g4 * x ^ 4)
          rw [abs_mul, abs_of_nonneg (by positivity : (0:ℝ) ≤ x ^ 4)] at a4
          have a3 := abs_add_le (g0 + g1 * x + g2 * x ^ 2) (g3 * x ^ 3)
          rw [abs_mul, abs_of_nonneg (by positivity : (0:ℝ) ≤ x ^ 3)] at a3
          have a2 := abs_add_le (g0 + g1 * x) (g2 * x ^ 2)
          rw [abs_mul, abs_of_nonneg (by positivity : (0:ℝ) ≤ x ^ 2)] at a2
          have a1 := abs_add_le (g0) (g1 * x)
          rw [abs_mul, abs_of_pos hpos] at a1
          linarith
      _ ≤ (M + |g0| + |g1| + |g2| + |g3| + |g4| + |g5| + |g6|) * x ^ 7 := by
          have := hM x hx
          have e1 : |g1| * x ≤ |g1| * x ^ 7 := mul_le_mul_of_nonneg_left (by simpa using hx17) (abs_nonneg _)
          have e2 : |g2| * x ^ 2 ≤ |g2| * x ^ 7 := mul_le_mul_of_nonneg_left hx27 (abs_nonneg _)
          have e3 : |g3| * x ^ 3 ≤ |g3| * x ^ 7 := mul_le_mul_of_nonneg_left hx37 (abs_nonneg _)
          have e4 : |g4| * x ^ 4 ≤ |g4| * x ^ 7 := mul_le_mul_of_nonneg_left hx47 (abs_nonneg _)
          have e5 : |g5| * x ^ 5 ≤ |g5| * x ^ 7 := mul_le_mul_of_nonneg_left hx57 (abs_nonneg _)
          have e6 : |g6| * x ^ 6 ≤ |g6| * x ^ 7 := mul_le_mul_of_nonneg_left hx67 (abs_nonneg _)
          have eM : M ≤ M * x ^ 7 := le_mul_of_one_le_right hM0 h17
          have e0 : |g0| ≤ |g0| * x ^ 7 := le_mul_of_one_le_right (abs_nonneg _) h17
          nlinarith
      _ ≤ _ := by gcongr; exact le_max_right _ _

/-- the dominating function of the order-8 Watson remainder -/
noncomputable def B9 (s : ℝ) : ℝ := (1 + s + s ^ 2 + s ^ 3 + s ^ 4 + s ^ 5 + s ^ 6 + s ^ 7) * Real.exp (-s)

lemma B9_int : IntegrableOn B9 (Set.Ioi 0) := by
  have h := (((((((Y04.intOn_pow_exp1 0).add (Y04.intOn_pow_exp1 1)).add (Y04.intOn_pow_exp1 2)).add (Y04.intOn_pow_exp1 3)).add (Y04.intOn_pow_exp1 4)).add (Y04.intOn_pow_exp1 5)).add (Y04.intOn_pow_exp1 6)).add (Y04.intOn_pow_exp1 7)
  refine IntegrableOn.congr_fun h (fun s _ => ?_) measurableSet_Ioi
  simp only [B9, Pi.add_apply]; ring

lemma B9_nonneg {s : ℝ} (hs : 0 ≤ s) : 0 ≤ B9 s := by unfold B9; positivity

lemma sj_le_B9 {s : ℝ} (hs : 0 ≤ s) {j : ℕ} (hj : j ≤ 7) : s ^ j * Real.exp (-s) ≤ B9 s := by
  unfold B9
  have h0 : 0 ≤ s ^ 0 * Real.exp (-s) := by positivity
  have h1 : 0 ≤ s ^ 1 * Real.exp (-s) := by positivity
  have h2 : 0 ≤ s ^ 2 * Real.exp (-s) := by positivity
  have h3 : 0 ≤ s ^ 3 * Real.exp (-s) := by positivity
  have h4 : 0 ≤ s ^ 4 * Real.exp (-s) := by positivity
  have h5 : 0 ≤ s ^ 5 * Real.exp (-s) := by positivity
  have h6 : 0 ≤ s ^ 6 * Real.exp (-s) := by positivity
  have h7 : 0 ≤ s ^ 7 * Real.exp (-s) := by positivity
  have w : (1 + s + s ^ 2 + s ^ 3 + s ^ 4 + s ^ 5 + s ^ 6 + s ^ 7) * Real.exp (-s)
    = s ^ 0 * Real.exp (-s) + s ^ 1 * Real.exp (-s) + s ^ 2 * Real.exp (-s) + s ^ 3 * Real.exp (-s) + s ^ 4 * Real.exp (-s) + s ^ 5 * Real.exp (-s) + s ^ 6 * Real.exp (-s) + s ^ 7 * Real.exp (-s) := by ring
  rw [w]
  interval_cases j <;> linarith


lemma hterm {g s τ L c : ℝ} {m j a : ℕ} (hs : 0 ≤ s) (hτ0 : 0 ≤ τ) (hτ1 : τ ≤ 1) (hL0 : 0 ≤ L)
    (hL1 : L ≤ 1) (hLe : L ≤ Real.exp (-s)) (hj : 1 ≤ j) (hc : 1 ≤ c) :
    |g * s ^ m * (τ ^ a * L ^ j / c)| ≤ |g| * (s ^ m * Real.exp (-s)) := Q03.hterm hs hτ0 hτ1 hL0 hL1 hLe hj hc

set_option maxHeartbeats 16000000 in
/-- the pointwise order-8 Watson remainder -/
lemma watson_pt8 {g0 g1 g2 g3 g4 g5 g6 CR M τ s : ℝ} (hCR : 0 ≤ CR) (hM : 0 ≤ M) (hτ0 : 0 < τ) (hτ1 : τ < 1)
    (hs : 0 ≤ s) {Gv : ℝ} (hGv : |Gv| ≤ M)
    (hR : |Gv - (g0 + g1 * (τ * s) + g2 * (τ * s) ^ 2 + g3 * (τ * s) ^ 3 + g4 * (τ * s) ^ 4 + g5 * (τ * s) ^ 5 + g6 * (τ * s) ^ 6)| ≤ CR * (τ * s) ^ 7) :
    |τ * ((Real.exp (τ * Y04.ell s) - 1) * Gv)
      - (τ ^ 2 * g0 * (s ^ 0 * Y04.ell s ^ 1)
        + τ ^ 3 * g1 * (s ^ 1 * Y04.ell s ^ 1)
        + τ ^ 3 * (g0 / 2) * (s ^ 0 * Y04.ell s ^ 2)
        + τ ^ 4 * g2 * (s ^ 2 * Y04.ell s ^ 1)
        + τ ^ 4 * (g1 / 2) * (s ^ 1 * Y04.ell s ^ 2)
        + τ ^ 4 * (g0 / 6) * (s ^ 0 * Y04.ell s ^ 3)
        + τ ^ 5 * g3 * (s ^ 3 * Y04.ell s ^ 1)
        + τ ^ 5 * (g2 / 2) * (s ^ 2 * Y04.ell s ^ 2)
        + τ ^ 5 * (g1 / 6) * (s ^ 1 * Y04.ell s ^ 3)
        + τ ^ 5 * (g0 / 24) * (s ^ 0 * Y04.ell s ^ 4)
        + τ ^ 6 * g4 * (s ^ 4 * Y04.ell s ^ 1)
        + τ ^ 6 * (g3 / 2) * (s ^ 3 * Y04.ell s ^ 2)
        + τ ^ 6 * (g2 / 6) * (s ^ 2 * Y04.ell s ^ 3)
        + τ ^ 6 * (g1 / 24) * (s ^ 1 * Y04.ell s ^ 4)
        + τ ^ 6 * (g0 / 120) * (s ^ 0 * Y04.ell s ^ 5)
        + τ ^ 7 * g5 * (s ^ 5 * Y04.ell s ^ 1)
        + τ ^ 7 * (g4 / 2) * (s ^ 4 * Y04.ell s ^ 2)
        + τ ^ 7 * (g3 / 6) * (s ^ 3 * Y04.ell s ^ 3)
        + τ ^ 7 * (g2 / 24) * (s ^ 2 * Y04.ell s ^ 4)
        + τ ^ 7 * (g1 / 120) * (s ^ 1 * Y04.ell s ^ 5)
        + τ ^ 7 * (g0 / 720) * (s ^ 0 * Y04.ell s ^ 6)
        + τ ^ 8 * g6 * (s ^ 6 * Y04.ell s ^ 1)
        + τ ^ 8 * (g5 / 2) * (s ^ 5 * Y04.ell s ^ 2)
        + τ ^ 8 * (g4 / 6) * (s ^ 4 * Y04.ell s ^ 3)
        + τ ^ 8 * (g3 / 24) * (s ^ 3 * Y04.ell s ^ 4)
        + τ ^ 8 * (g2 / 120) * (s ^ 2 * Y04.ell s ^ 5)
        + τ ^ 8 * (g1 / 720) * (s ^ 1 * Y04.ell s ^ 6)
        + τ ^ 8 * (g0 / 5040) * (s ^ 0 * Y04.ell s ^ 7))|
      ≤ τ ^ 9 * ((6 * (|g1| + |g2| + |g3| + |g4| + |g5| + |g6|) + 2 * CR + M) * B9 s) := by
  have hl0 := Y04.ell_nonneg s
  have hl1 := Y04.ell_le_one s hs
  have hle := Y04.ell_le s
  set L := Y04.ell s
  have hτL : |τ * L| ≤ 1 := by rw [abs_of_nonneg (by positivity)]; nlinarith
  have he := Real.exp_bound hτL (n := 8) (by norm_num)
  simp only [Finset.sum_range_succ, Finset.sum_range_zero, Nat.factorial] at he
  norm_num at he
  set a1 := τ * L + (τ * L) ^ 2 / 2 + (τ * L) ^ 3 / 6 + (τ * L) ^ 4 / 24 + (τ * L) ^ 5 / 120 + (τ * L) ^ 6 / 720 + (τ * L) ^ 7 / 5040 with ha1
  set r1 := Real.exp (τ * L) - 1 - a1 with hr1
  have hr1b : |r1| ≤ τ ^ 8 * L ^ 8 := by
    have : r1 = Real.exp (τ * L) - (1 + τ * L + (τ * L) ^ 2 / 2 + (τ * L) ^ 3 / 6 + (τ * L) ^ 4 / 24 + (τ * L) ^ 5 / 120 + (τ * L) ^ 6 / 720 + (τ * L) ^ 7 / 5040) := by
      rw [hr1, ha1]; ring
    rw [this]
    have h4 : (|τ| * |L|) ^ 8 = τ ^ 8 * L ^ 8 := by
      rw [abs_of_pos hτ0, abs_of_nonneg hl0]; ring
    rw [h4] at he
    have : 0 ≤ τ ^ 8 * L ^ 8 := by positivity
    linarith
  set q := g0 + g1 * (τ * s) + g2 * (τ * s) ^ 2 + g3 * (τ * s) ^ 3 + g4 * (τ * s) ^ 4 + g5 * (τ * s) ^ 5 + g6 * (τ * s) ^ 6
  set R := Gv - q with hRdef
  set Q := g1 * s ^ 1 * (τ ^ 0 * L ^ 7 / 5040) + g2 * s ^ 2 * (τ ^ 0 * L ^ 6 / 720) + g2 * s ^ 2 * (τ ^ 1 * L ^ 7 / 5040) + g3 * s ^ 3 * (τ ^ 0 * L ^ 5 / 120) + g3 * s ^ 3 * (τ ^ 1 * L ^ 6 / 720) + g3 * s ^ 3 * (τ ^ 2 * L ^ 7 / 5040) + g4 * s ^ 4 * (τ ^ 0 * L ^ 4 / 24) + g4 * s ^ 4 * (τ ^ 1 * L ^ 5 / 120) + g4 * s ^ 4 * (τ ^ 2 * L ^ 6 / 720) + g4 * s ^ 4 * (τ ^ 3 * L ^ 7 / 5040) + g5 * s ^ 5 * (τ ^ 0 * L ^ 3 / 6) + g5 * s ^ 5 * (τ ^ 1 * L ^ 4 / 24) + g5 * s ^ 5 * (τ ^ 2 * L ^ 5 / 120) + g5 * s ^ 5 * (τ ^ 3 * L ^ 6 / 720) + g5 * s ^ 5 * (τ ^ 4 * L ^ 7 / 5040) + g6 * s ^ 6 * (τ ^ 0 * L ^ 2 / 2) + g6 * s ^ 6 * (τ ^ 1 * L ^ 3 / 6) + g6 * s ^ 6 * (τ ^ 2 * L ^ 4 / 24) + g6 * s ^ 6 * (τ ^ 3 * L ^ 5 / 120) + g6 * s ^ 6 * (τ ^ 4 * L ^ 6 / 720) + g6 * s ^ 6 * (τ ^ 5 * L ^ 7 / 5040) with hQ
  have key : τ * ((Real.exp (τ * L) - 1) * Gv)
      - (τ ^ 2 * g0 * (s ^ 0 * L ^ 1)
        + τ ^ 3 * g1 * (s ^ 1 * L ^ 1)
        + τ ^ 3 * (g0 / 2) * (s ^ 0 * L ^ 2)
        + τ ^ 4 * g2 * (s ^ 2 * L ^ 1)
        + τ ^ 4 * (g1 / 2) * (s ^ 1 * L ^ 2)
        + τ ^ 4 * (g0 / 6) * (s ^ 0 * L ^ 3)
        + τ ^ 5 * g3 * (s ^ 3 * L ^ 1)
        + τ ^ 5 * (g2 / 2) * (s ^ 2 * L ^ 2)
        + τ ^ 5 * (g1 / 6) * (s ^ 1 * L ^ 3)
        + τ ^ 5 * (g0 / 24) * (s ^ 0 * L ^ 4)
        + τ ^ 6 * g4 * (s ^ 4 * L ^ 1)
        + τ ^ 6 * (g3 / 2) * (s ^ 3 * L ^ 2)
        + τ ^ 6 * (g2 / 6) * (s ^ 2 * L ^ 3)
        + τ ^ 6 * (g1 / 24) * (s ^ 1 * L ^ 4)
        + τ ^ 6 * (g0 / 120) * (s ^ 0 * L ^ 5)
        + τ ^ 7 * g5 * (s ^ 5 * L ^ 1)
        + τ ^ 7 * (g4 / 2) * (s ^ 4 * L ^ 2)
        + τ ^ 7 * (g3 / 6) * (s ^ 3 * L ^ 3)
        + τ ^ 7 * (g2 / 24) * (s ^ 2 * L ^ 4)
        + τ ^ 7 * (g1 / 120) * (s ^ 1 * L ^ 5)
        + τ ^ 7 * (g0 / 720) * (s ^ 0 * L ^ 6)
        + τ ^ 8 * g6 * (s ^ 6 * L ^ 1)
        + τ ^ 8 * (g5 / 2) * (s ^ 5 * L ^ 2)
        + τ ^ 8 * (g4 / 6) * (s ^ 4 * L ^ 3)
        + τ ^ 8 * (g3 / 24) * (s ^ 3 * L ^ 4)
        + τ ^ 8 * (g2 / 120) * (s ^ 2 * L ^ 5)
        + τ ^ 8 * (g1 / 720) * (s ^ 1 * L ^ 6)
        + τ ^ 8 * (g0 / 5040) * (s ^ 0 * L ^ 7))
      = τ ^ 9 * Q + τ * a1 * R + τ * r1 * Gv := by
    rw [hRdef, hr1, ha1, hQ]; ring
  rw [key]
  have hes : 0 ≤ Real.exp (-s) := (Real.exp_pos _).le
  have hτ' : 0 ≤ τ := hτ0.le
  have q0 : |g1 * s ^ 1 * (τ ^ 0 * L ^ 7 / 5040)| ≤ |g1| * (s ^ 1 * Real.exp (-s)) := hterm hs hτ' hτ1.le hl0 hl1 hle (by norm_num) (by norm_num)
  have q1 : |g2 * s ^ 2 * (τ ^ 0 * L ^ 6 / 720)| ≤ |g2| * (s ^ 2 * Real.exp (-s)) := hterm hs hτ' hτ1.le hl0 hl1 hle (by norm_num) (by norm_num)
  have q2 : |g2 * s ^ 2 * (τ ^ 1 * L ^ 7 / 5040)| ≤ |g2| * (s ^ 2 * Real.exp (-s)) := hterm hs hτ' hτ1.le hl0 hl1 hle (by norm_num) (by norm_num)
  have q3 : |g3 * s ^ 3 * (τ ^ 0 * L ^ 5 / 120)| ≤ |g3| * (s ^ 3 * Real.exp (-s)) := hterm hs hτ' hτ1.le hl0 hl1 hle (by norm_num) (by norm_num)
  have q4 : |g3 * s ^ 3 * (τ ^ 1 * L ^ 6 / 720)| ≤ |g3| * (s ^ 3 * Real.exp (-s)) := hterm hs hτ' hτ1.le hl0 hl1 hle (by norm_num) (by norm_num)
  have q5 : |g3 * s ^ 3 * (τ ^ 2 * L ^ 7 / 5040)| ≤ |g3| * (s ^ 3 * Real.exp (-s)) := hterm hs hτ' hτ1.le hl0 hl1 hle (by norm_num) (by norm_num)
  have q6 : |g4 * s ^ 4 * (τ ^ 0 * L ^ 4 / 24)| ≤ |g4| * (s ^ 4 * Real.exp (-s)) := hterm hs hτ' hτ1.le hl0 hl1 hle (by norm_num) (by norm_num)
  have q7 : |g4 * s ^ 4 * (τ ^ 1 * L ^ 5 / 120)| ≤ |g4| * (s ^ 4 * Real.exp (-s)) := hterm hs hτ' hτ1.le hl0 hl1 hle (by norm_num) (by norm_num)
  have q8 : |g4 * s ^ 4 * (τ ^ 2 * L ^ 6 / 720)| ≤ |g4| * (s ^ 4 * Real.exp (-s)) := hterm hs hτ' hτ1.le hl0 hl1 hle (by norm_num) (by norm_num)
  have q9 : |g4 * s ^ 4 * (τ ^ 3 * L ^ 7 / 5040)| ≤ |g4| * (s ^ 4 * Real.exp (-s)) := hterm hs hτ' hτ1.le hl0 hl1 hle (by norm_num) (by norm_num)
  have q10 : |g5 * s ^ 5 * (τ ^ 0 * L ^ 3 / 6)| ≤ |g5| * (s ^ 5 * Real.exp (-s)) := hterm hs hτ' hτ1.le hl0 hl1 hle (by norm_num) (by norm_num)
  have q11 : |g5 * s ^ 5 * (τ ^ 1 * L ^ 4 / 24)| ≤ |g5| * (s ^ 5 * Real.exp (-s)) := hterm hs hτ' hτ1.le hl0 hl1 hle (by norm_num) (by norm_num)
  have q12 : |g5 * s ^ 5 * (τ ^ 2 * L ^ 5 / 120)| ≤ |g5| * (s ^ 5 * Real.exp (-s)) := hterm hs hτ' hτ1.le hl0 hl1 hle (by norm_num) (by norm_num)
  have q13 : |g5 * s ^ 5 * (τ ^ 3 * L ^ 6 / 720)| ≤ |g5| * (s ^ 5 * Real.exp (-s)) := hterm hs hτ' hτ1.le hl0 hl1 hle (by norm_num) (by norm_num)
  have q14 : |g5 * s ^ 5 * (τ ^ 4 * L ^ 7 / 5040)| ≤ |g5| * (s ^ 5 * Real.exp (-s)) := hterm hs hτ' hτ1.le hl0 hl1 hle (by norm_num) (by norm_num)
  have q15 : |g6 * s ^ 6 * (τ ^ 0 * L ^ 2 / 2)| ≤ |g6| * (s ^ 6 * Real.exp (-s)) := hterm hs hτ' hτ1.le hl0 hl1 hle (by norm_num) (by norm_num)
  have q16 : |g6 * s ^ 6 * (τ ^ 1 * L ^ 3 / 6)| ≤ |g6| * (s ^ 6 * Real.exp (-s)) := hterm hs hτ' hτ1.le hl0 hl1 hle (by norm_num) (by norm_num)
  have q17 : |g6 * s ^ 6 * (τ ^ 2 * L ^ 4 / 24)| ≤ |g6| * (s ^ 6 * Real.exp (-s)) := hterm hs hτ' hτ1.le hl0 hl1 hle (by norm_num) (by norm_num)
  have q18 : |g6 * s ^ 6 * (τ ^ 3 * L ^ 5 / 120)| ≤ |g6| * (s ^ 6 * Real.exp (-s)) := hterm hs hτ' hτ1.le hl0 hl1 hle (by norm_num) (by norm_num)
  have q19 : |g6 * s ^ 6 * (τ ^ 4 * L ^ 6 / 720)| ≤ |g6| * (s ^ 6 * Real.exp (-s)) := hterm hs hτ' hτ1.le hl0 hl1 hle (by norm_num) (by norm_num)
  have q20 : |g6 * s ^ 6 * (τ ^ 5 * L ^ 7 / 5040)| ≤ |g6| * (s ^ 6 * Real.exp (-s)) := hterm hs hτ' hτ1.le hl0 hl1 hle (by norm_num) (by norm_num)
  have c1 := abs_add_le (g1 * s ^ 1 * (τ ^ 0 * L ^ 7 / 5040)) (g2 * s ^ 2 * (τ ^ 0 * L ^ 6 / 720))
  have c2 := abs_add_le (g1 * s ^ 1 * (τ ^ 0 * L ^ 7 / 5040) + g2 * s ^ 2 * (τ ^ 0 * L ^ 6 / 720)) (g2 * s ^ 2 * (τ ^ 1 * L ^ 7 / 5040))
  have c3 := abs_add_le (g1 * s ^ 1 * (τ ^ 0 * L ^ 7 / 5040) + g2 * s ^ 2 * (τ ^ 0 * L ^ 6 / 720) + g2 * s ^ 2 * (τ ^ 1 * L ^ 7 / 5040)) (g3 * s ^ 3 * (τ ^ 0 * L ^ 5 / 120))
  have c4 := abs_add_le (g1 * s ^ 1 * (τ ^ 0 * L ^ 7 / 5040) + g2 * s ^ 2 * (τ ^ 0 * L ^ 6 / 720) + g2 * s ^ 2 * (τ ^ 1 * L ^ 7 / 5040) + g3 * s ^ 3 * (τ ^ 0 * L ^ 5 / 120)) (g3 * s ^ 3 * (τ ^ 1 * L ^ 6 / 720))
  have c5 := abs_add_le (g1 * s ^ 1 * (τ ^ 0 * L ^ 7 / 5040) + g2 * s ^ 2 * (τ ^ 0 * L ^ 6 / 720) + g2 * s ^ 2 * (τ ^ 1 * L ^ 7 / 5040) + g3 * s ^ 3 * (τ ^ 0 * L ^ 5 / 120) + g3 * s ^ 3 * (τ ^ 1 * L ^ 6 / 720)) (g3 * s ^ 3 * (τ ^ 2 * L ^ 7 / 5040))
  have c6 := abs_add_le (g1 * s ^ 1 * (τ ^ 0 * L ^ 7 / 5040) + g2 * s ^ 2 * (τ ^ 0 * L ^ 6 / 720) + g2 * s ^ 2 * (τ ^ 1 * L ^ 7 / 5040) + g3 * s ^ 3 * (τ ^ 0 * L ^ 5 / 120) + g3 * s ^ 3 * (τ ^ 1 * L ^ 6 / 720) + g3 * s ^ 3 * (τ ^ 2 * L ^ 7 / 5040)) (g4 * s ^ 4 * (τ ^ 0 * L ^ 4 / 24))
  have c7 := abs_add_le (g1 * s ^ 1 * (τ ^ 0 * L ^ 7 / 5040) + g2 * s ^ 2 * (τ ^ 0 * L ^ 6 / 720) + g2 * s ^ 2 * (τ ^ 1 * L ^ 7 / 5040) + g3 * s ^ 3 * (τ ^ 0 * L ^ 5 / 120) + g3 * s ^ 3 * (τ ^ 1 * L ^ 6 / 720) + g3 * s ^ 3 * (τ ^ 2 * L ^ 7 / 5040) + g4 * s ^ 4 * (τ ^ 0 * L ^ 4 / 24)) (g4 * s ^ 4 * (τ ^ 1 * L ^ 5 / 120))
  have c8 := abs_add_le (g1 * s ^ 1 * (τ ^ 0 * L ^ 7 / 5040) + g2 * s ^ 2 * (τ ^ 0 * L ^ 6 / 720) + g2 * s ^ 2 * (τ ^ 1 * L ^ 7 / 5040) + g3 * s ^ 3 * (τ ^ 0 * L ^ 5 / 120) + g3 * s ^ 3 * (τ ^ 1 * L ^ 6 / 720) + g3 * s ^ 3 * (τ ^ 2 * L ^ 7 / 5040) + g4 * s ^ 4 * (τ ^ 0 * L ^ 4 / 24) + g4 * s ^ 4 * (τ ^ 1 * L ^ 5 / 120)) (g4 * s ^ 4 * (τ ^ 2 * L ^ 6 / 720))
  have c9 := abs_add_le (g1 * s ^ 1 * (τ ^ 0 * L ^ 7 / 5040) + g2 * s ^ 2 * (τ ^ 0 * L ^ 6 / 720) + g2 * s ^ 2 * (τ ^ 1 * L ^ 7 / 5040) + g3 * s ^ 3 * (τ ^ 0 * L ^ 5 / 120) + g3 * s ^ 3 * (τ ^ 1 * L ^ 6 / 720) + g3 * s ^ 3 * (τ ^ 2 * L ^ 7 / 5040) + g4 * s ^ 4 * (τ ^ 0 * L ^ 4 / 24) + g4 * s ^ 4 * (τ ^ 1 * L ^ 5 / 120) + g4 * s ^ 4 * (τ ^ 2 * L ^ 6 / 720)) (g4 * s ^ 4 * (τ ^ 3 * L ^ 7 / 5040))
  have c10 := abs_add_le (g1 * s ^ 1 * (τ ^ 0 * L ^ 7 / 5040) + g2 * s ^ 2 * (τ ^ 0 * L ^ 6 / 720) + g2 * s ^ 2 * (τ ^ 1 * L ^ 7 / 5040) + g3 * s ^ 3 * (τ ^ 0 * L ^ 5 / 120) + g3 * s ^ 3 * (τ ^ 1 * L ^ 6 / 720) + g3 * s ^ 3 * (τ ^ 2 * L ^ 7 / 5040) + g4 * s ^ 4 * (τ ^ 0 * L ^ 4 / 24) + g4 * s ^ 4 * (τ ^ 1 * L ^ 5 / 120) + g4 * s ^ 4 * (τ ^ 2 * L ^ 6 / 720) + g4 * s ^ 4 * (τ ^ 3 * L ^ 7 / 5040)) (g5 * s ^ 5 * (τ ^ 0 * L ^ 3 / 6))
  have c11 := abs_add_le (g1 * s ^ 1 * (τ ^ 0 * L ^ 7 / 5040) + g2 * s ^ 2 * (τ ^ 0 * L ^ 6 / 720) + g2 * s ^ 2 * (τ ^ 1 * L ^ 7 / 5040) + g3 * s ^ 3 * (τ ^ 0 * L ^ 5 / 120) + g3 * s ^ 3 * (τ ^ 1 * L ^ 6 / 720) + g3 * s ^ 3 * (τ ^ 2 * L ^ 7 / 5040) + g4 * s ^ 4 * (τ ^ 0 * L ^ 4 / 24) + g4 * s ^ 4 * (τ ^ 1 * L ^ 5 / 120) + g4 * s ^ 4 * (τ ^ 2 * L ^ 6 / 720) + g4 * s ^ 4 * (τ ^ 3 * L ^ 7 / 5040) + g5 * s ^ 5 * (τ ^ 0 * L ^ 3 / 6)) (g5 * s ^ 5 * (τ ^ 1 * L ^ 4 / 24))
  have c12 := abs_add_le (g1 * s ^ 1 * (τ ^ 0 * L ^ 7 / 5040) + g2 * s ^ 2 * (τ ^ 0 * L ^ 6 / 720) + g2 * s ^ 2 * (τ ^ 1 * L ^ 7 / 5040) + g3 * s ^ 3 * (τ ^ 0 * L ^ 5 / 120) + g3 * s ^ 3 * (τ ^ 1 * L ^ 6 / 720) + g3 * s ^ 3 * (τ ^ 2 * L ^ 7 / 5040) + g4 * s ^ 4 * (τ ^ 0 * L ^ 4 / 24) + g4 * s ^ 4 * (τ ^ 1 * L ^ 5 / 120) + g4 * s ^ 4 * (τ ^ 2 * L ^ 6 / 720) + g4 * s ^ 4 * (τ ^ 3 * L ^ 7 / 5040) + g5 * s ^ 5 * (τ ^ 0 * L ^ 3 / 6) + g5 * s ^ 5 * (τ ^ 1 * L ^ 4 / 24)) (g5 * s ^ 5 * (τ ^ 2 * L ^ 5 / 120))
  have c13 := abs_add_le (g1 * s ^ 1 * (τ ^ 0 * L ^ 7 / 5040) + g2 * s ^ 2 * (τ ^ 0 * L ^ 6 / 720) + g2 * s ^ 2 * (τ ^ 1 * L ^ 7 / 5040) + g3 * s ^ 3 * (τ ^ 0 * L ^ 5 / 120) + g3 * s ^ 3 * (τ ^ 1 * L ^ 6 / 720) + g3 * s ^ 3 * (τ ^ 2 * L ^ 7 / 5040) + g4 * s ^ 4 * (τ ^ 0 * L ^ 4 / 24) + g4 * s ^ 4 * (τ ^ 1 * L ^ 5 / 120) + g4 * s ^ 4 * (τ ^ 2 * L ^ 6 / 720) + g4 * s ^ 4 * (τ ^ 3 * L ^ 7 / 5040) + g5 * s ^ 5 * (τ ^ 0 * L ^ 3 / 6) + g5 * s ^ 5 * (τ ^ 1 * L ^ 4 / 24) + g5 * s ^ 5 * (τ ^ 2 * L ^ 5 / 120)) (g5 * s ^ 5 * (τ ^ 3 * L ^ 6 / 720))
  have c14 := abs_add_le (g1 * s ^ 1 * (τ ^ 0 * L ^ 7 / 5040) + g2 * s ^ 2 * (τ ^ 0 * L ^ 6 / 720) + g2 * s ^ 2 * (τ ^ 1 * L ^ 7 / 5040) + g3 * s ^ 3 * (τ ^ 0 * L ^ 5 / 120) + g3 * s ^ 3 * (τ ^ 1 * L ^ 6 / 720) + g3 * s ^ 3 * (τ ^ 2 * L ^ 7 / 5040) + g4 * s ^ 4 * (τ ^ 0 * L ^ 4 / 24) + g4 * s ^ 4 * (τ ^ 1 * L ^ 5 / 120) + g4 * s ^ 4 * (τ ^ 2 * L ^ 6 / 720) + g4 * s ^ 4 * (τ ^ 3 * L ^ 7 / 5040) + g5 * s ^ 5 * (τ ^ 0 * L ^ 3 / 6) + g5 * s ^ 5 * (τ ^ 1 * L ^ 4 / 24) + g5 * s ^ 5 * (τ ^ 2 * L ^ 5 / 120) + g5 * s ^ 5 * (τ ^ 3 * L ^ 6 / 720)) (g5 * s ^ 5 * (τ ^ 4 * L ^ 7 / 5040))
  have c15 := abs_add_le (g1 * s ^ 1 * (τ ^ 0 * L ^ 7 / 5040) + g2 * s ^ 2 * (τ ^ 0 * L ^ 6 / 720) + g2 * s ^ 2 * (τ ^ 1 * L ^ 7 / 5040) + g3 * s ^ 3 * (τ ^ 0 * L ^ 5 / 120) + g3 * s ^ 3 * (τ ^ 1 * L ^ 6 / 720) + g3 * s ^ 3 * (τ ^ 2 * L ^ 7 / 5040) + g4 * s ^ 4 * (τ ^ 0 * L ^ 4 / 24) + g4 * s ^ 4 * (τ ^ 1 * L ^ 5 / 120) + g4 * s ^ 4 * (τ ^ 2 * L ^ 6 / 720) + g4 * s ^ 4 * (τ ^ 3 * L ^ 7 / 5040) + g5 * s ^ 5 * (τ ^ 0 * L ^ 3 / 6) + g5 * s ^ 5 * (τ ^ 1 * L ^ 4 / 24) + g5 * s ^ 5 * (τ ^ 2 * L ^ 5 / 120) + g5 * s ^ 5 * (τ ^ 3 * L ^ 6 / 720) + g5 * s ^ 5 * (τ ^ 4 * L ^ 7 / 5040)) (g6 * s ^ 6 * (τ ^ 0 * L ^ 2 / 2))
  have c16 := abs_add_le (g1 * s ^ 1 * (τ ^ 0 * L ^ 7 / 5040) + g2 * s ^ 2 * (τ ^ 0 * L ^ 6 / 720) + g2 * s ^ 2 * (τ ^ 1 * L ^ 7 / 5040) + g3 * s ^ 3 * (τ ^ 0 * L ^ 5 / 120) + g3 * s ^ 3 * (τ ^ 1 * L ^ 6 / 720) + g3 * s ^ 3 * (τ ^ 2 * L ^ 7 / 5040) + g4 * s ^ 4 * (τ ^ 0 * L ^ 4 / 24) + g4 * s ^ 4 * (τ ^ 1 * L ^ 5 / 120) + g4 * s ^ 4 * (τ ^ 2 * L ^ 6 / 720) + g4 * s ^ 4 * (τ ^ 3 * L ^ 7 / 5040) + g5 * s ^ 5 * (τ ^ 0 * L ^ 3 / 6) + g5 * s ^ 5 * (τ ^ 1 * L ^ 4 / 24) + g5 * s ^ 5 * (τ ^ 2 * L ^ 5 / 120) + g5 * s ^ 5 * (τ ^ 3 * L ^ 6 / 720) + g5 * s ^ 5 * (τ ^ 4 * L ^ 7 / 5040) + g6 * s ^ 6 * (τ ^ 0 * L ^ 2 / 2)) (g6 * s ^ 6 * (τ ^ 1 * L ^ 3 / 6))
  have c17 := abs_add_le (g1 * s ^ 1 * (τ ^ 0 * L ^ 7 / 5040) + g2 * s ^ 2 * (τ ^ 0 * L ^ 6 / 720) + g2 * s ^ 2 * (τ ^ 1 * L ^ 7 / 5040) + g3 * s ^ 3 * (τ ^ 0 * L ^ 5 / 120) + g3 * s ^ 3 * (τ ^ 1 * L ^ 6 / 720) + g3 * s ^ 3 * (τ ^ 2 * L ^ 7 / 5040) + g4 * s ^ 4 * (τ ^ 0 * L ^ 4 / 24) + g4 * s ^ 4 * (τ ^ 1 * L ^ 5 / 120) + g4 * s ^ 4 * (τ ^ 2 * L ^ 6 / 720) + g4 * s ^ 4 * (τ ^ 3 * L ^ 7 / 5040) + g5 * s ^ 5 * (τ ^ 0 * L ^ 3 / 6) + g5 * s ^ 5 * (τ ^ 1 * L ^ 4 / 24) + g5 * s ^ 5 * (τ ^ 2 * L ^ 5 / 120) + g5 * s ^ 5 * (τ ^ 3 * L ^ 6 / 720) + g5 * s ^ 5 * (τ ^ 4 * L ^ 7 / 5040) + g6 * s ^ 6 * (τ ^ 0 * L ^ 2 / 2) + g6 * s ^ 6 * (τ ^ 1 * L ^ 3 / 6)) (g6 * s ^ 6 * (τ ^ 2 * L ^ 4 / 24))
  have c18 := abs_add_le (g1 * s ^ 1 * (τ ^ 0 * L ^ 7 / 5040) + g2 * s ^ 2 * (τ ^ 0 * L ^ 6 / 720) + g2 * s ^ 2 * (τ ^ 1 * L ^ 7 / 5040) + g3 * s ^ 3 * (τ ^ 0 * L ^ 5 / 120) + g3 * s ^ 3 * (τ ^ 1 * L ^ 6 / 720) + g3 * s ^ 3 * (τ ^ 2 * L ^ 7 / 5040) + g4 * s ^ 4 * (τ ^ 0 * L ^ 4 / 24) + g4 * s ^ 4 * (τ ^ 1 * L ^ 5 / 120) + g4 * s ^ 4 * (τ ^ 2 * L ^ 6 / 720) + g4 * s ^ 4 * (τ ^ 3 * L ^ 7 / 5040) + g5 * s ^ 5 * (τ ^ 0 * L ^ 3 / 6) + g5 * s ^ 5 * (τ ^ 1 * L ^ 4 / 24) + g5 * s ^ 5 * (τ ^ 2 * L ^ 5 / 120) + g5 * s ^ 5 * (τ ^ 3 * L ^ 6 / 720) + g5 * s ^ 5 * (τ ^ 4 * L ^ 7 / 5040) + g6 * s ^ 6 * (τ ^ 0 * L ^ 2 / 2) + g6 * s ^ 6 * (τ ^ 1 * L ^ 3 / 6) + g6 * s ^ 6 * (τ ^ 2 * L ^ 4 / 24)) (g6 * s ^ 6 * (τ ^ 3 * L ^ 5 / 120))
  have c19 := abs_add_le (g1 * s ^ 1 * (τ ^ 0 * L ^ 7 / 5040) + g2 * s ^ 2 * (τ ^ 0 * L ^ 6 / 720) + g2 * s ^ 2 * (τ ^ 1 * L ^ 7 / 5040) + g3 * s ^ 3 * (τ ^ 0 * L ^ 5 / 120) + g3 * s ^ 3 * (τ ^ 1 * L ^ 6 / 720) + g3 * s ^ 3 * (τ ^ 2 * L ^ 7 / 5040) + g4 * s ^ 4 * (τ ^ 0 * L ^ 4 / 24) + g4 * s ^ 4 * (τ ^ 1 * L ^ 5 / 120) + g4 * s ^ 4 * (τ ^ 2 * L ^ 6 / 720) + g4 * s ^ 4 * (τ ^ 3 * L ^ 7 / 5040) + g5 * s ^ 5 * (τ ^ 0 * L ^ 3 / 6) + g5 * s ^ 5 * (τ ^ 1 * L ^ 4 / 24) + g5 * s ^ 5 * (τ ^ 2 * L ^ 5 / 120) + g5 * s ^ 5 * (τ ^ 3 * L ^ 6 / 720) + g5 * s ^ 5 * (τ ^ 4 * L ^ 7 / 5040) + g6 * s ^ 6 * (τ ^ 0 * L ^ 2 / 2) + g6 * s ^ 6 * (τ ^ 1 * L ^ 3 / 6) + g6 * s ^ 6 * (τ ^ 2 * L ^ 4 / 24) + g6 * s ^ 6 * (τ ^ 3 * L ^ 5 / 120)) (g6 * s ^ 6 * (τ ^ 4 * L ^ 6 / 720))
  have c20 := abs_add_le (g1 * s ^ 1 * (τ ^ 0 * L ^ 7 / 5040) + g2 * s ^ 2 * (τ ^ 0 * L ^ 6 / 720) + g2 * s ^ 2 * (τ ^ 1 * L ^ 7 / 5040) + g3 * s ^ 3 * (τ ^ 0 * L ^ 5 / 120) + g3 * s ^ 3 * (τ ^ 1 * L ^ 6 / 720) + g3 * s ^ 3 * (τ ^ 2 * L ^ 7 / 5040) + g4 * s ^ 4 * (τ ^ 0 * L ^ 4 / 24) + g4 * s ^ 4 * (τ ^ 1 * L ^ 5 / 120) + g4 * s ^ 4 * (τ ^ 2 * L ^ 6 / 720) + g4 * s ^ 4 * (τ ^ 3 * L ^ 7 / 5040) + g5 * s ^ 5 * (τ ^ 0 * L ^ 3 / 6) + g5 * s ^ 5 * (τ ^ 1 * L ^ 4 / 24) + g5 * s ^ 5 * (τ ^ 2 * L ^ 5 / 120) + g5 * s ^ 5 * (τ ^ 3 * L ^ 6 / 720) + g5 * s ^ 5 * (τ ^ 4 * L ^ 7 / 5040) + g6 * s ^ 6 * (τ ^ 0 * L ^ 2 / 2) + g6 * s ^ 6 * (τ ^ 1 * L ^ 3 / 6) + g6 * s ^ 6 * (τ ^ 2 * L ^ 4 / 24) + g6 * s ^ 6 * (τ ^ 3 * L ^ 5 / 120) + g6 * s ^ 6 * (τ ^ 4 * L ^ 6 / 720)) (g6 * s ^ 6 * (τ ^ 5 * L ^ 7 / 5040))
  have hB : ∀ m : ℕ, m ≤ 7 → s ^ m * Real.exp (-s) ≤ B9 s := fun m hm => sj_le_B9 hs hm
  have hB0 : 0 ≤ B9 s := B9_nonneg hs
  have hQb : |Q| ≤ 6 * (|g1| + |g2| + |g3| + |g4| + |g5| + |g6|) * B9 s := by
    have k1 := mul_le_mul_of_nonneg_left (hB 1 (by norm_num)) (abs_nonneg g1)
    have n1 := mul_nonneg (abs_nonneg g1) hB0
    have k2 := mul_le_mul_of_nonneg_left (hB 2 (by norm_num)) (abs_nonneg g2)
    have n2 := mul_nonneg (abs_nonneg g2) hB0
    have k3 := mul_le_mul_of_nonneg_left (hB 3 (by norm_num)) (abs_nonneg g3)
    have n3 := mul_nonneg (abs_nonneg g3) hB0
    have k4 := mul_le_mul_of_nonneg_left (hB 4 (by norm_num)) (abs_nonneg g4)
    have n4 := mul_nonneg (abs_nonneg g4) hB0
    have k5 := mul_le_mul_of_nonneg_left (hB 5 (by norm_num)) (abs_nonneg g5)
    have n5 := mul_nonneg (abs_nonneg g5) hB0
    have k6 := mul_le_mul_of_nonneg_left (hB 6 (by norm_num)) (abs_nonneg g6)
    have n6 := mul_nonneg (abs_nonneg g6) hB0
    rw [hQ]
    linarith
  have b2 : |τ * a1 * R| ≤ τ ^ 9 * (2 * CR * s ^ 7 * Real.exp (-s)) := by
    have ha1b : |a1| ≤ 2 * τ * L := by
      have hu0 : 0 ≤ τ * L := mul_nonneg hτ0.le hl0
      have hu1 : τ * L ≤ 1 := by
        have := mul_le_mul hτ1.le hl1 hl0 zero_le_one; linarith
      have p2 : (τ * L) ^ 2 ≤ τ * L := pow_le_of_le_one hu0 hu1 (by norm_num)
      have z2 : 0 ≤ (τ * L) ^ 2 := pow_nonneg hu0 2
      have p3 : (τ * L) ^ 3 ≤ τ * L := pow_le_of_le_one hu0 hu1 (by norm_num)
      have z3 : 0 ≤ (τ * L) ^ 3 := pow_nonneg hu0 3
      have p4 : (τ * L) ^ 4 ≤ τ * L := pow_le_of_le_one hu0 hu1 (by norm_num)
      have z4 : 0 ≤ (τ * L) ^ 4 := pow_nonneg hu0 4
      have p5 : (τ * L) ^ 5 ≤ τ * L := pow_le_of_le_one hu0 hu1 (by norm_num)
      have z5 : 0 ≤ (τ * L) ^ 5 := pow_nonneg hu0 5
      have p6 : (τ * L) ^ 6 ≤ τ * L := pow_le_of_le_one hu0 hu1 (by norm_num)
      have z6 : 0 ≤ (τ * L) ^ 6 := pow_nonneg hu0 6
      have p7 : (τ * L) ^ 7 ≤ τ * L := pow_le_of_le_one hu0 hu1 (by norm_num)
      have z7 : 0 ≤ (τ * L) ^ 7 := pow_nonneg hu0 7
      have : 0 ≤ a1 := by rw [ha1]; linarith
      rw [abs_of_nonneg this, ha1]
      have : 2 * τ * L = 2 * (τ * L) := by ring
      linarith
    have hRb : |R| ≤ CR * (τ * s) ^ 7 := hR
    rw [abs_mul, abs_mul, abs_of_pos hτ0]
    calc τ * |a1| * |R| ≤ τ * (2 * τ * L) * (CR * (τ * s) ^ 7) :=
          mul_le_mul (mul_le_mul_of_nonneg_left ha1b hτ') hRb (abs_nonneg _) (by positivity)
      _ = τ ^ 9 * (2 * CR * s ^ 7 * L) := by ring
      _ ≤ τ ^ 9 * (2 * CR * s ^ 7 * Real.exp (-s)) :=
          mul_le_mul_of_nonneg_left (mul_le_mul_of_nonneg_left hle (by positivity)) (by positivity)
  have hL8 : L ^ 8 ≤ Real.exp (-s) := (pow_le_of_le_one hl0 hl1 (by norm_num)).trans hle
  have b3 : |τ * r1 * Gv| ≤ τ ^ 9 * (M * Real.exp (-s)) := by
    rw [abs_mul, abs_mul, abs_of_pos hτ0]
    calc τ * |r1| * |Gv| ≤ τ * (τ ^ 8 * L ^ 8) * M :=
          mul_le_mul (mul_le_mul_of_nonneg_left hr1b hτ0.le) hGv (abs_nonneg _) (by positivity)
      _ = τ ^ 9 * (L ^ 8 * M) := by ring
      _ ≤ τ ^ 9 * (Real.exp (-s) * M) :=
          mul_le_mul_of_nonneg_left (mul_le_mul_of_nonneg_right hL8 hM) (by positivity)
      _ = _ := by ring
  have t1 := abs_add_le (τ ^ 9 * Q + τ * a1 * R) (τ * r1 * Gv)
  have t2 := abs_add_le (τ ^ 9 * Q) (τ * a1 * R)
  have b1 : |τ ^ 9 * Q| ≤ τ ^ 9 * (6 * (|g1| + |g2| + |g3| + |g4| + |g5| + |g6|) * B9 s) := by
    rw [abs_mul, abs_of_pos (by positivity)]
    exact mul_le_mul_of_nonneg_left hQb (by positivity)
  have hin : 6 * (|g1| + |g2| + |g3| + |g4| + |g5| + |g6|) * B9 s + 2 * CR * s ^ 7 * Real.exp (-s)
      + M * Real.exp (-s) ≤ (6 * (|g1| + |g2| + |g3| + |g4| + |g5| + |g6|) + 2 * CR + M) * B9 s := by
    have k7 := mul_le_mul_of_nonneg_left (hB 7 (by norm_num)) (by positivity : (0:ℝ) ≤ 2 * CR)
    have k0 := mul_le_mul_of_nonneg_left (hB 0 (by norm_num)) hM
    simp only [pow_zero, one_mul] at k0
    nlinarith
  have hfin := mul_le_mul_of_nonneg_left hin (by positivity : (0:ℝ) ≤ τ ^ 9)
  have e : τ ^ 9 * (6 * (|g1| + |g2| + |g3| + |g4| + |g5| + |g6|) * B9 s + 2 * CR * s ^ 7 * Real.exp (-s)
      + M * Real.exp (-s)) = τ ^ 9 * (6 * (|g1| + |g2| + |g3| + |g4| + |g5| + |g6|) * B9 s)
      + τ ^ 9 * (2 * CR * s ^ 7 * Real.exp (-s)) + τ ^ 9 * (M * Real.exp (-s)) := by ring
  linarith

lemma int28 {f1 f2 f3 f4 f5 f6 f7 f8 f9 f10 f11 f12 f13 f14 f15 f16 f17 f18 f19 f20 f21 f22 f23 f24 f25 f26 f27 f28 : ℝ → ℝ}
    (c1 c2 c3 c4 c5 c6 c7 c8 c9 c10 c11 c12 c13 c14 c15 c16 c17 c18 c19 c20 c21 c22 c23 c24 c25 c26 c27 c28 : ℝ)
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
    :
    IntegrableOn (fun s => c1 * f1 s + c2 * f2 s + c3 * f3 s + c4 * f4 s + c5 * f5 s + c6 * f6 s + c7 * f7 s + c8 * f8 s + c9 * f9 s + c10 * f10 s + c11 * f11 s + c12 * f12 s + c13 * f13 s + c14 * f14 s + c15 * f15 s + c16 * f16 s + c17 * f17 s + c18 * f18 s + c19 * f19 s + c20 * f20 s + c21 * f21 s + c22 * f22 s + c23 * f23 s + c24 * f24 s + c25 * f25 s + c26 * f26 s + c27 * f27 s + c28 * f28 s) (Set.Ioi 0) ∧
    ∫ s in Set.Ioi (0 : ℝ), (c1 * f1 s + c2 * f2 s + c3 * f3 s + c4 * f4 s + c5 * f5 s + c6 * f6 s + c7 * f7 s + c8 * f8 s + c9 * f9 s + c10 * f10 s + c11 * f11 s + c12 * f12 s + c13 * f13 s + c14 * f14 s + c15 * f15 s + c16 * f16 s + c17 * f17 s + c18 * f18 s + c19 * f19 s + c20 * f20 s + c21 * f21 s + c22 * f22 s + c23 * f23 s + c24 * f24 s + c25 * f25 s + c26 * f26 s + c27 * f27 s + c28 * f28 s)
      = c1 * (∫ s in Set.Ioi (0 : ℝ), f1 s) + c2 * (∫ s in Set.Ioi (0 : ℝ), f2 s) + c3 * (∫ s in Set.Ioi (0 : ℝ), f3 s) + c4 * (∫ s in Set.Ioi (0 : ℝ), f4 s) + c5 * (∫ s in Set.Ioi (0 : ℝ), f5 s) + c6 * (∫ s in Set.Ioi (0 : ℝ), f6 s) + c7 * (∫ s in Set.Ioi (0 : ℝ), f7 s) + c8 * (∫ s in Set.Ioi (0 : ℝ), f8 s) + c9 * (∫ s in Set.Ioi (0 : ℝ), f9 s) + c10 * (∫ s in Set.Ioi (0 : ℝ), f10 s) + c11 * (∫ s in Set.Ioi (0 : ℝ), f11 s) + c12 * (∫ s in Set.Ioi (0 : ℝ), f12 s) + c13 * (∫ s in Set.Ioi (0 : ℝ), f13 s) + c14 * (∫ s in Set.Ioi (0 : ℝ), f14 s) + c15 * (∫ s in Set.Ioi (0 : ℝ), f15 s) + c16 * (∫ s in Set.Ioi (0 : ℝ), f16 s) + c17 * (∫ s in Set.Ioi (0 : ℝ), f17 s) + c18 * (∫ s in Set.Ioi (0 : ℝ), f18 s) + c19 * (∫ s in Set.Ioi (0 : ℝ), f19 s) + c20 * (∫ s in Set.Ioi (0 : ℝ), f20 s) + c21 * (∫ s in Set.Ioi (0 : ℝ), f21 s) + c22 * (∫ s in Set.Ioi (0 : ℝ), f22 s) + c23 * (∫ s in Set.Ioi (0 : ℝ), f23 s) + c24 * (∫ s in Set.Ioi (0 : ℝ), f24 s) + c25 * (∫ s in Set.Ioi (0 : ℝ), f25 s) + c26 * (∫ s in Set.Ioi (0 : ℝ), f26 s) + c27 * (∫ s in Set.Ioi (0 : ℝ), f27 s) + c28 * (∫ s in Set.Ioi (0 : ℝ), f28 s) := by
  obtain ⟨b21, e21⟩ := Q03.int21 c1 c2 c3 c4 c5 c6 c7 c8 c9 c10 c11 c12 c13 c14 c15 c16 c17 c18 c19 c20 c21 h1 h2 h3 h4 h5 h6 h7 h8 h9 h10 h11 h12 h13 h14 h15 h16 h17 h18 h19 h20 h21
  have a22 : Integrable (fun s => c22 * f22 s) (volume.restrict (Set.Ioi 0)) := h22.const_mul c22
  have a23 : Integrable (fun s => c23 * f23 s) (volume.restrict (Set.Ioi 0)) := h23.const_mul c23
  have a24 : Integrable (fun s => c24 * f24 s) (volume.restrict (Set.Ioi 0)) := h24.const_mul c24
  have a25 : Integrable (fun s => c25 * f25 s) (volume.restrict (Set.Ioi 0)) := h25.const_mul c25
  have a26 : Integrable (fun s => c26 * f26 s) (volume.restrict (Set.Ioi 0)) := h26.const_mul c26
  have a27 : Integrable (fun s => c27 * f27 s) (volume.restrict (Set.Ioi 0)) := h27.const_mul c27
  have a28 : Integrable (fun s => c28 * f28 s) (volume.restrict (Set.Ioi 0)) := h28.const_mul c28
  have b22 : Integrable (fun s => c1 * f1 s + c2 * f2 s + c3 * f3 s + c4 * f4 s + c5 * f5 s + c6 * f6 s + c7 * f7 s + c8 * f8 s + c9 * f9 s + c10 * f10 s + c11 * f11 s + c12 * f12 s + c13 * f13 s + c14 * f14 s + c15 * f15 s + c16 * f16 s + c17 * f17 s + c18 * f18 s + c19 * f19 s + c20 * f20 s + c21 * f21 s + c22 * f22 s) (volume.restrict (Set.Ioi 0)) := b21.add a22
  have b23 : Integrable (fun s => c1 * f1 s + c2 * f2 s + c3 * f3 s + c4 * f4 s + c5 * f5 s + c6 * f6 s + c7 * f7 s + c8 * f8 s + c9 * f9 s + c10 * f10 s + c11 * f11 s + c12 * f12 s + c13 * f13 s + c14 * f14 s + c15 * f15 s + c16 * f16 s + c17 * f17 s + c18 * f18 s + c19 * f19 s + c20 * f20 s + c21 * f21 s + c22 * f22 s + c23 * f23 s) (volume.restrict (Set.Ioi 0)) := b22.add a23
  have b24 : Integrable (fun s => c1 * f1 s + c2 * f2 s + c3 * f3 s + c4 * f4 s + c5 * f5 s + c6 * f6 s + c7 * f7 s + c8 * f8 s + c9 * f9 s + c10 * f10 s + c11 * f11 s + c12 * f12 s + c13 * f13 s + c14 * f14 s + c15 * f15 s + c16 * f16 s + c17 * f17 s + c18 * f18 s + c19 * f19 s + c20 * f20 s + c21 * f21 s + c22 * f22 s + c23 * f23 s + c24 * f24 s) (volume.restrict (Set.Ioi 0)) := b23.add a24
  have b25 : Integrable (fun s => c1 * f1 s + c2 * f2 s + c3 * f3 s + c4 * f4 s + c5 * f5 s + c6 * f6 s + c7 * f7 s + c8 * f8 s + c9 * f9 s + c10 * f10 s + c11 * f11 s + c12 * f12 s + c13 * f13 s + c14 * f14 s + c15 * f15 s + c16 * f16 s + c17 * f17 s + c18 * f18 s + c19 * f19 s + c20 * f20 s + c21 * f21 s + c22 * f22 s + c23 * f23 s + c24 * f24 s + c25 * f25 s) (volume.restrict (Set.Ioi 0)) := b24.add a25
  have b26 : Integrable (fun s => c1 * f1 s + c2 * f2 s + c3 * f3 s + c4 * f4 s + c5 * f5 s + c6 * f6 s + c7 * f7 s + c8 * f8 s + c9 * f9 s + c10 * f10 s + c11 * f11 s + c12 * f12 s + c13 * f13 s + c14 * f14 s + c15 * f15 s + c16 * f16 s + c17 * f17 s + c18 * f18 s + c19 * f19 s + c20 * f20 s + c21 * f21 s + c22 * f22 s + c23 * f23 s + c24 * f24 s + c25 * f25 s + c26 * f26 s) (volume.restrict (Set.Ioi 0)) := b25.add a26
  have b27 : Integrable (fun s => c1 * f1 s + c2 * f2 s + c3 * f3 s + c4 * f4 s + c5 * f5 s + c6 * f6 s + c7 * f7 s + c8 * f8 s + c9 * f9 s + c10 * f10 s + c11 * f11 s + c12 * f12 s + c13 * f13 s + c14 * f14 s + c15 * f15 s + c16 * f16 s + c17 * f17 s + c18 * f18 s + c19 * f19 s + c20 * f20 s + c21 * f21 s + c22 * f22 s + c23 * f23 s + c24 * f24 s + c25 * f25 s + c26 * f26 s + c27 * f27 s) (volume.restrict (Set.Ioi 0)) := b26.add a27
  have b28 : Integrable (fun s => c1 * f1 s + c2 * f2 s + c3 * f3 s + c4 * f4 s + c5 * f5 s + c6 * f6 s + c7 * f7 s + c8 * f8 s + c9 * f9 s + c10 * f10 s + c11 * f11 s + c12 * f12 s + c13 * f13 s + c14 * f14 s + c15 * f15 s + c16 * f16 s + c17 * f17 s + c18 * f18 s + c19 * f19 s + c20 * f20 s + c21 * f21 s + c22 * f22 s + c23 * f23 s + c24 * f24 s + c25 * f25 s + c26 * f26 s + c27 * f27 s + c28 * f28 s) (volume.restrict (Set.Ioi 0)) := b27.add a28
  refine ⟨b28, ?_⟩
  rw [integral_add b27 a28, integral_add b26 a27, integral_add b25 a26, integral_add b24 a25, integral_add b23 a24, integral_add b22 a23, integral_add b21 a22, e21, integral_const_mul, integral_const_mul, integral_const_mul, integral_const_mul, integral_const_mul, integral_const_mul, integral_const_mul]

lemma watson8 (Gf : ℝ → ℝ) (hG : ContDiff ℝ (⊤ : ℕ∞) Gf) (hb : ∃ M : ℝ, ∀ x : ℝ, 0 ≤ x → |Gf x| ≤ M) :
    ∃ C : ℝ, ∀ τ : ℝ, 0 < τ → τ < 1 →
      IntegrableOn (fun s : ℝ => ((1 + Real.exp (-s)) ^ τ - 1) * Gf (τ * s)) (Ioi 0) ∧
      |τ * (∫ s in Ioi (0 : ℝ), ((1 + Real.exp (-s)) ^ τ - 1) * Gf (τ * s)) -
          ∑ n ∈ Finset.Icc 2 8, τ ^ n * ∑ j ∈ Finset.Icc 1 (n - 1),
            iteratedDeriv (n - 1 - j) Gf 0 / ((n - 1 - j).factorial : ℝ) /
              (j.factorial : ℝ) * Y04.I (n - 1 - j) j| ≤ C * τ ^ 9 := by
  obtain ⟨M, hM⟩ := hb
  have hM0 : 0 ≤ M := (abs_nonneg _).trans (hM 0 le_rfl)
  obtain ⟨CR, hCR, hT⟩ := taylor6 hG hM
  set g0 := iteratedDeriv 0 Gf 0 / (Nat.factorial 0 : ℝ) with hg0
  set g1 := iteratedDeriv 1 Gf 0 / (Nat.factorial 1 : ℝ) with hg1
  set g2 := iteratedDeriv 2 Gf 0 / (Nat.factorial 2 : ℝ) with hg2
  set g3 := iteratedDeriv 3 Gf 0 / (Nat.factorial 3 : ℝ) with hg3
  set g4 := iteratedDeriv 4 Gf 0 / (Nat.factorial 4 : ℝ) with hg4
  set g5 := iteratedDeriv 5 Gf 0 / (Nat.factorial 5 : ℝ) with hg5
  set g6 := iteratedDeriv 6 Gf 0 / (Nat.factorial 6 : ℝ) with hg6
  refine ⟨(6 * (|g1| + |g2| + |g3| + |g4| + |g5| + |g6|) + 2 * CR + M) * ∫ s in Set.Ioi (0 : ℝ), B9 s, fun τ hτ0 hτ1 => ?_⟩
  have e : (fun s : ℝ => ((1 + Real.exp (-s)) ^ τ - 1) * Gf (τ * s))
      = fun s => (Real.exp (τ * Y04.ell s) - 1) * Gf (τ * s) := by
    funext s; rw [Y04.rpow_eq_exp]
  rw [e]
  have hint : IntegrableOn (fun s => (Real.exp (τ * Y04.ell s) - 1) * Gf (τ * s)) (Set.Ioi 0) := by
    have hGc : Continuous Gf := hG.continuous
    refine Y04.intOn_of_bound (by have := Y04.continuous_ell; fun_prop)
      ((Y04.intOn_pow_exp1 0).const_mul (2 * M)) (fun s hs => ?_)
    have hl0 := Y04.ell_nonneg s
    have hl1 := Y04.ell_le_one s hs.le
    have hτL : |τ * Y04.ell s| ≤ 1 := by
      rw [abs_of_nonneg (mul_nonneg hτ0.le hl0)]
      have := mul_le_mul hτ1.le hl1 hl0 zero_le_one; linarith
    have h1 := Real.abs_exp_sub_one_le hτL
    have h2 : |τ * Y04.ell s| ≤ Real.exp (-s) := by
      rw [abs_of_nonneg (mul_nonneg hτ0.le hl0)]
      have := mul_le_mul hτ1.le (Y04.ell_le s) hl0 zero_le_one; linarith
    have h3 := hM (τ * s) (by nlinarith)
    rw [abs_mul]
    simp only [pow_zero, one_mul]
    calc |Real.exp (τ * Y04.ell s) - 1| * |Gf (τ * s)| ≤ (2 * Real.exp (-s)) * M :=
          mul_le_mul (by linarith) h3 (abs_nonneg _) (by positivity)
      _ = 2 * M * Real.exp (-s) := by ring
  refine ⟨hint, ?_⟩
  obtain ⟨hmi, hmv⟩ := int28 (τ ^ 2 * g0) (τ ^ 3 * g1) (τ ^ 3 * (g0 / 2)) (τ ^ 4 * g2) (τ ^ 4 * (g1 / 2)) (τ ^ 4 * (g0 / 6)) (τ ^ 5 * g3) (τ ^ 5 * (g2 / 2)) (τ ^ 5 * (g1 / 6)) (τ ^ 5 * (g0 / 24)) (τ ^ 6 * g4) (τ ^ 6 * (g3 / 2)) (τ ^ 6 * (g2 / 6)) (τ ^ 6 * (g1 / 24)) (τ ^ 6 * (g0 / 120)) (τ ^ 7 * g5) (τ ^ 7 * (g4 / 2)) (τ ^ 7 * (g3 / 6)) (τ ^ 7 * (g2 / 24)) (τ ^ 7 * (g1 / 120)) (τ ^ 7 * (g0 / 720)) (τ ^ 8 * g6) (τ ^ 8 * (g5 / 2)) (τ ^ 8 * (g4 / 6)) (τ ^ 8 * (g3 / 24)) (τ ^ 8 * (g2 / 120)) (τ ^ 8 * (g1 / 720)) (τ ^ 8 * (g0 / 5040))
    (Y04.I_int 0 1 le_rfl) (Y04.I_int 1 1 le_rfl) (Y04.I_int 0 2 (by norm_num)) (Y04.I_int 2 1 le_rfl) (Y04.I_int 1 2 (by norm_num)) (Y04.I_int 0 3 (by norm_num)) (Y04.I_int 3 1 le_rfl) (Y04.I_int 2 2 (by norm_num)) (Y04.I_int 1 3 (by norm_num)) (Y04.I_int 0 4 (by norm_num)) (Y04.I_int 4 1 le_rfl) (Y04.I_int 3 2 (by norm_num)) (Y04.I_int 2 3 (by norm_num)) (Y04.I_int 1 4 (by norm_num)) (Y04.I_int 0 5 (by norm_num)) (Y04.I_int 5 1 le_rfl) (Y04.I_int 4 2 (by norm_num)) (Y04.I_int 3 3 (by norm_num)) (Y04.I_int 2 4 (by norm_num)) (Y04.I_int 1 5 (by norm_num)) (Y04.I_int 0 6 (by norm_num)) (Y04.I_int 6 1 le_rfl) (Y04.I_int 5 2 (by norm_num)) (Y04.I_int 4 3 (by norm_num)) (Y04.I_int 3 4 (by norm_num)) (Y04.I_int 2 5 (by norm_num)) (Y04.I_int 1 6 (by norm_num)) (Y04.I_int 0 7 (by norm_num))
  have hsum : ∑ n ∈ Finset.Icc 2 8, τ ^ n * ∑ j ∈ Finset.Icc 1 (n - 1),
            iteratedDeriv (n - 1 - j) Gf 0 / ((n - 1 - j).factorial : ℝ) /
              (j.factorial : ℝ) * Y04.I (n - 1 - j) j
      = τ ^ 2 * g0 * Y04.I 0 1 + τ ^ 3 * g1 * Y04.I 1 1 + τ ^ 3 * (g0 / 2) * Y04.I 0 2 + τ ^ 4 * g2 * Y04.I 2 1 + τ ^ 4 * (g1 / 2) * Y04.I 1 2 + τ ^ 4 * (g0 / 6) * Y04.I 0 3 + τ ^ 5 * g3 * Y04.I 3 1 + τ ^ 5 * (g2 / 2) * Y04.I 2 2 + τ ^ 5 * (g1 / 6) * Y04.I 1 3 + τ ^ 5 * (g0 / 24) * Y04.I 0 4 + τ ^ 6 * g4 * Y04.I 4 1 + τ ^ 6 * (g3 / 2) * Y04.I 3 2 + τ ^ 6 * (g2 / 6) * Y04.I 2 3 + τ ^ 6 * (g1 / 24) * Y04.I 1 4 + τ ^ 6 * (g0 / 120) * Y04.I 0 5 + τ ^ 7 * g5 * Y04.I 5 1 + τ ^ 7 * (g4 / 2) * Y04.I 4 2 + τ ^ 7 * (g3 / 6) * Y04.I 3 3 + τ ^ 7 * (g2 / 24) * Y04.I 2 4 + τ ^ 7 * (g1 / 120) * Y04.I 1 5 + τ ^ 7 * (g0 / 720) * Y04.I 0 6 + τ ^ 8 * g6 * Y04.I 6 1 + τ ^ 8 * (g5 / 2) * Y04.I 5 2 + τ ^ 8 * (g4 / 6) * Y04.I 4 3 + τ ^ 8 * (g3 / 24) * Y04.I 3 4 + τ ^ 8 * (g2 / 120) * Y04.I 2 5 + τ ^ 8 * (g1 / 720) * Y04.I 1 6 + τ ^ 8 * (g0 / 5040) * Y04.I 0 7 := by
    rw [show Finset.Icc 2 8 = {2, 3, 4, 5, 6, 7, 8} by decide]
    rw [Finset.sum_insert (by decide), Finset.sum_insert (by decide), Finset.sum_insert (by decide), Finset.sum_insert (by decide), Finset.sum_insert (by decide), Finset.sum_insert (by decide), Finset.sum_singleton]
    rw [show Finset.Icc 1 (2 - 1) = {1} by decide, show Finset.Icc 1 (3 - 1) = {1, 2} by decide, show Finset.Icc 1 (4 - 1) = {1, 2, 3} by decide, show Finset.Icc 1 (5 - 1) = {1, 2, 3, 4} by decide, show Finset.Icc 1 (6 - 1) = {1, 2, 3, 4, 5} by decide, show Finset.Icc 1 (7 - 1) = {1, 2, 3, 4, 5, 6} by decide, show Finset.Icc 1 (8 - 1) = {1, 2, 3, 4, 5, 6, 7} by decide]
    simp [Finset.sum_insert, hg0, hg1, hg2, hg3, hg4, hg5, hg6, Nat.factorial]
    ring
  rw [hsum]
  unfold Y04.I
  rw [← hmv, ← integral_const_mul, ← integral_sub (hint.const_mul τ) hmi]
  have hB := (B9_int.const_mul (6 * (|g1| + |g2| + |g3| + |g4| + |g5| + |g6|) + 2 * CR + M)).const_mul (τ ^ 9)
  have := norm_integral_le_of_norm_le hB
    ((ae_restrict_iff' measurableSet_Ioi).2 (Filter.Eventually.of_forall
      (fun s (hs : 0 < s) => by
        rw [Real.norm_eq_abs]
        have hx : 0 ≤ τ * s := by positivity
        have := watson_pt8 (g0 := g0) (g1 := g1) (g2 := g2) (g3 := g3) (g4 := g4) (g5 := g5) (g6 := g6) hCR hM0 hτ0 hτ1 hs.le
          (hM (τ * s) hx)
          (by have := hT (τ * s) hx; simpa [hg0, hg1, hg2, hg3, hg4, hg5, hg6] using this)
        convert this using 2)))
  rw [Real.norm_eq_abs, integral_const_mul, integral_const_mul] at this
  calc _ ≤ _ := this
    _ = _ := by ring


lemma em7_bound {q : ℝ} (hq : 0 ≤ q) :
    |Real.exp (-q) - (1 - q + q ^ 2 / 2 - q ^ 3 / 6 + q ^ 4 / 24 - q ^ 5 / 120 + q ^ 6 / 720)| ≤ 4 * q ^ 7 := by
  rcases le_or_gt q 1 with h | h
  · have h1 : |-q| ≤ 1 := by rw [abs_neg, abs_of_nonneg hq]; exact h
    have he := Real.exp_bound h1 (n := 7) (by norm_num)
    simp only [Finset.sum_range_succ, Finset.sum_range_zero, Nat.factorial] at he
    norm_num at he
    rw [abs_of_nonneg hq] at he
    have e : Real.exp (-q) - (1 - q + q ^ 2 / 2 - q ^ 3 / 6 + q ^ 4 / 24 - q ^ 5 / 120 + q ^ 6 / 720)
        = Real.exp (-q) - (1 + -q + q ^ 2 / 2 + (-q) ^ 3 / 6 + q ^ 4 / 24 + (-q) ^ 5 / 120 + q ^ 6 / 720) := by ring
    rw [e]
    have : 0 ≤ q ^ 7 := by positivity
    nlinarith
  · have he : Real.exp (-q) ≤ 1 := Real.exp_le_one_iff.2 (by linarith)
    have he0 := Real.exp_pos (-q)
    have h2 : q ≤ q ^ 2 := by nlinarith
    have h3 : q ^ 2 ≤ q ^ 3 := by nlinarith
    have h4 : q ^ 3 ≤ q ^ 4 := by nlinarith
    have h5 : q ^ 4 ≤ q ^ 5 := by nlinarith
    have h6 : q ^ 5 ≤ q ^ 6 := by nlinarith
    have h7 : q ^ 6 ≤ q ^ 7 := by nlinarith
    rw [abs_le]
    constructor <;> nlinarith

lemma ep8_bound {p : ℝ} (hp0 : 0 ≤ p) (hp1 : p ≤ 1) :
    |Real.exp (-p) - 1 - (-p + p ^ 2 / 2 - p ^ 3 / 6 + p ^ 4 / 24 - p ^ 5 / 120 + p ^ 6 / 720 - p ^ 7 / 5040)| ≤ p ^ 8 := by
  have h1 : |-p| ≤ 1 := by rw [abs_neg, abs_of_nonneg hp0]; exact hp1
  have he := Real.exp_bound h1 (n := 8) (by norm_num)
  simp only [Finset.sum_range_succ, Finset.sum_range_zero, Nat.factorial] at he
  norm_num at he
  rw [abs_of_nonneg hp0] at he
  have e : Real.exp (-p) - 1 - (-p + p ^ 2 / 2 - p ^ 3 / 6 + p ^ 4 / 24 - p ^ 5 / 120 + p ^ 6 / 720 - p ^ 7 / 5040)
      = Real.exp (-p) - (1 + -p + p ^ 2 / 2 + (-p) ^ 3 / 6 + p ^ 4 / 24 + (-p) ^ 5 / 120
        + p ^ 6 / 720 + (-p) ^ 7 / 5040) := by ring
  rw [e]
  have : 0 ≤ p ^ 8 := by positivity
  nlinarith

lemma mono_bd8 {u v t L s : ℝ} (hu0 : 0 ≤ u) (hut : u ≤ t) (hv0 : 0 ≤ v) (hvt : v ≤ t) (ht1 : t ≤ 1)
    (hL0 : 0 ≤ L) (hL1 : L ≤ 1) (hLe : L ≤ Real.exp (-s)) (hs : 0 ≤ s) {i j : ℕ} (hi : 1 ≤ i)
    (hij : 8 ≤ i + j) :
    (u * L) ^ i * (v * s) ^ j ≤ t ^ 8 * (s ^ j * Real.exp (-s)) := by
  have ht0 : 0 ≤ t := hu0.trans hut
  rw [mul_pow, mul_pow]
  have a1 : u ^ i ≤ t ^ i := pow_le_pow_left₀ hu0 hut i
  have a2 : v ^ j ≤ t ^ j := pow_le_pow_left₀ hv0 hvt j
  have a3 : t ^ i * t ^ j ≤ t ^ 8 := by
    rw [← pow_add]; exact pow_le_pow_of_le_one ht0 ht1 hij
  have a4 : L ^ i ≤ Real.exp (-s) := (pow_le_of_le_one hL0 hL1 (by omega)).trans hLe
  have b1 : u ^ i * v ^ j ≤ t ^ 8 :=
    (mul_le_mul a1 a2 (by positivity) (by positivity)).trans a3
  calc u ^ i * L ^ i * (v ^ j * s ^ j) = (u ^ i * v ^ j) * L ^ i * s ^ j := by ring
    _ ≤ t ^ 8 * Real.exp (-s) * s ^ j := by
        apply mul_le_mul_of_nonneg_right _ (by positivity)
        exact mul_le_mul b1 a4 (by positivity) (by positivity)
    _ = _ := by ring

lemma hmono8 {u v t L s p q c : ℝ} (hu0 : 0 ≤ u) (hut : u ≤ t) (hv0 : 0 ≤ v) (hvt : v ≤ t) (ht1 : t ≤ 1)
    (hL0 : 0 ≤ L) (hL1 : L ≤ 1) (hLe : L ≤ Real.exp (-s)) (hs : 0 ≤ s) (hp : p = u * L) (hq : q = v * s)
    (hc : |c| ≤ 1) {i j : ℕ} (hi : 1 ≤ i) (hij : 8 ≤ i + j) (hj : j ≤ 7) :
    |c * (p ^ i * q ^ j)| ≤ t ^ 8 * B9 s := by
  subst hp hq
  have h0 : 0 ≤ (u * L) ^ i * (v * s) ^ j := by positivity
  have h1 := mono_bd8 hu0 hut hv0 hvt ht1 hL0 hL1 hLe hs hi hij
  have h2 := sj_le_B9 hs hj
  have ht0 : 0 ≤ t := hu0.trans hut
  rw [abs_mul, abs_of_nonneg h0]
  calc |c| * ((u * L) ^ i * (v * s) ^ j) ≤ 1 * (t ^ 8 * (s ^ j * Real.exp (-s))) :=
        mul_le_mul hc h1 h0 zero_le_one
    _ ≤ t ^ 8 * B9 s := by rw [one_mul]; exact mul_le_mul_of_nonneg_left h2 (by positivity)

/-- the degree-`≤ 7` part of `(e^{-uℓ} - 1) e^{-vs}` -/
noncomputable def P7 (u v s : ℝ) : ℝ :=
  (-1) * (u ^ 1 * v ^ 0) * (s ^ 0 * Y04.ell s ^ 1)
    + (1) * (u ^ 1 * v ^ 1) * (s ^ 1 * Y04.ell s ^ 1)
    + (-1 / 2) * (u ^ 1 * v ^ 2) * (s ^ 2 * Y04.ell s ^ 1)
    + (1 / 6) * (u ^ 1 * v ^ 3) * (s ^ 3 * Y04.ell s ^ 1)
    + (-1 / 24) * (u ^ 1 * v ^ 4) * (s ^ 4 * Y04.ell s ^ 1)
    + (1 / 120) * (u ^ 1 * v ^ 5) * (s ^ 5 * Y04.ell s ^ 1)
    + (-1 / 720) * (u ^ 1 * v ^ 6) * (s ^ 6 * Y04.ell s ^ 1)
    + (1 / 2) * (u ^ 2 * v ^ 0) * (s ^ 0 * Y04.ell s ^ 2)
    + (-1 / 2) * (u ^ 2 * v ^ 1) * (s ^ 1 * Y04.ell s ^ 2)
    + (1 / 4) * (u ^ 2 * v ^ 2) * (s ^ 2 * Y04.ell s ^ 2)
    + (-1 / 12) * (u ^ 2 * v ^ 3) * (s ^ 3 * Y04.ell s ^ 2)
    + (1 / 48) * (u ^ 2 * v ^ 4) * (s ^ 4 * Y04.ell s ^ 2)
    + (-1 / 240) * (u ^ 2 * v ^ 5) * (s ^ 5 * Y04.ell s ^ 2)
    + (-1 / 6) * (u ^ 3 * v ^ 0) * (s ^ 0 * Y04.ell s ^ 3)
    + (1 / 6) * (u ^ 3 * v ^ 1) * (s ^ 1 * Y04.ell s ^ 3)
    + (-1 / 12) * (u ^ 3 * v ^ 2) * (s ^ 2 * Y04.ell s ^ 3)
    + (1 / 36) * (u ^ 3 * v ^ 3) * (s ^ 3 * Y04.ell s ^ 3)
    + (-1 / 144) * (u ^ 3 * v ^ 4) * (s ^ 4 * Y04.ell s ^ 3)
    + (1 / 24) * (u ^ 4 * v ^ 0) * (s ^ 0 * Y04.ell s ^ 4)
    + (-1 / 24) * (u ^ 4 * v ^ 1) * (s ^ 1 * Y04.ell s ^ 4)
    + (1 / 48) * (u ^ 4 * v ^ 2) * (s ^ 2 * Y04.ell s ^ 4)
    + (-1 / 144) * (u ^ 4 * v ^ 3) * (s ^ 3 * Y04.ell s ^ 4)
    + (-1 / 120) * (u ^ 5 * v ^ 0) * (s ^ 0 * Y04.ell s ^ 5)
    + (1 / 120) * (u ^ 5 * v ^ 1) * (s ^ 1 * Y04.ell s ^ 5)
    + (-1 / 240) * (u ^ 5 * v ^ 2) * (s ^ 2 * Y04.ell s ^ 5)
    + (1 / 720) * (u ^ 6 * v ^ 0) * (s ^ 0 * Y04.ell s ^ 6)
    + (-1 / 720) * (u ^ 6 * v ^ 1) * (s ^ 1 * Y04.ell s ^ 6)
    + (-1 / 5040) * (u ^ 7 * v ^ 0) * (s ^ 0 * Y04.ell s ^ 7)

set_option maxHeartbeats 16000000 in
lemma pt7 {u v t s : ℝ} (hu0 : 0 ≤ u) (hut : u ≤ t) (hv0 : 0 ≤ v) (hvt : v ≤ t) (ht1 : t ≤ 1)
    (hs : 0 ≤ s) :
    |(Real.exp (-u * Y04.ell s) - 1) * Real.exp (-v * s) - P7 u v s| ≤ 30 * t ^ 8 * B9 s := by
  have hl0 := Y04.ell_nonneg s
  have hl1 := Y04.ell_le_one s hs
  have hle := Y04.ell_le s
  have ht0 : 0 ≤ t := hu0.trans hut
  unfold P7
  set L := Y04.ell s
  set p := u * L with hp
  set q := v * s with hq
  have hp0 : 0 ≤ p := by positivity
  have hp1 : p ≤ 1 := by
    have := mul_le_mul (hut.trans ht1) hl1 hl0 zero_le_one; rw [hp]; linarith
  have hq0 : 0 ≤ q := by positivity
  set A := -p + p ^ 2 / 2 - p ^ 3 / 6 + p ^ 4 / 24 - p ^ 5 / 120 + p ^ 6 / 720 - p ^ 7 / 5040 with hA
  set r1 := Real.exp (-p) - 1 - A with hr1
  set B := 1 - q + q ^ 2 / 2 - q ^ 3 / 6 + q ^ 4 / 24 - q ^ 5 / 120 + q ^ 6 / 720 with hB
  set r2 := Real.exp (-q) - B with hr2
  have hr1b : |r1| ≤ p ^ 8 := ep8_bound hp0 hp1
  have hr2b : |r2| ≤ 4 * q ^ 7 := em7_bound hq0
  set HT := (1 / 1440) * (p ^ 2 * q ^ 6) + (1 / 720) * (p ^ 3 * q ^ 5) + (-1 / 4320) * (p ^ 3 * q ^ 6) + (1 / 576) * (p ^ 4 * q ^ 4) + (-1 / 2880) * (p ^ 4 * q ^ 5) + (1 / 17280) * (p ^ 4 * q ^ 6) + (1 / 720) * (p ^ 5 * q ^ 3) + (-1 / 2880) * (p ^ 5 * q ^ 4) + (1 / 14400) * (p ^ 5 * q ^ 5) + (-1 / 86400) * (p ^ 5 * q ^ 6) + (1 / 1440) * (p ^ 6 * q ^ 2) + (-1 / 4320) * (p ^ 6 * q ^ 3) + (1 / 17280) * (p ^ 6 * q ^ 4) + (-1 / 86400) * (p ^ 6 * q ^ 5) + (1 / 518400) * (p ^ 6 * q ^ 6) + (1 / 5040) * (p ^ 7 * q ^ 1) + (-1 / 10080) * (p ^ 7 * q ^ 2) + (1 / 30240) * (p ^ 7 * q ^ 3) + (-1 / 120960) * (p ^ 7 * q ^ 4) + (1 / 604800) * (p ^ 7 * q ^ 5) + (-1 / 3628800) * (p ^ 7 * q ^ 6) with hHT
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
    + (1 / 2) * (u ^ 2 * v ^ 0) * (s ^ 0 * L ^ 2)
    + (-1 / 2) * (u ^ 2 * v ^ 1) * (s ^ 1 * L ^ 2)
    + (1 / 4) * (u ^ 2 * v ^ 2) * (s ^ 2 * L ^ 2)
    + (-1 / 12) * (u ^ 2 * v ^ 3) * (s ^ 3 * L ^ 2)
    + (1 / 48) * (u ^ 2 * v ^ 4) * (s ^ 4 * L ^ 2)
    + (-1 / 240) * (u ^ 2 * v ^ 5) * (s ^ 5 * L ^ 2)
    + (-1 / 6) * (u ^ 3 * v ^ 0) * (s ^ 0 * L ^ 3)
    + (1 / 6) * (u ^ 3 * v ^ 1) * (s ^ 1 * L ^ 3)
    + (-1 / 12) * (u ^ 3 * v ^ 2) * (s ^ 2 * L ^ 3)
    + (1 / 36) * (u ^ 3 * v ^ 3) * (s ^ 3 * L ^ 3)
    + (-1 / 144) * (u ^ 3 * v ^ 4) * (s ^ 4 * L ^ 3)
    + (1 / 24) * (u ^ 4 * v ^ 0) * (s ^ 0 * L ^ 4)
    + (-1 / 24) * (u ^ 4 * v ^ 1) * (s ^ 1 * L ^ 4)
    + (1 / 48) * (u ^ 4 * v ^ 2) * (s ^ 2 * L ^ 4)
    + (-1 / 144) * (u ^ 4 * v ^ 3) * (s ^ 3 * L ^ 4)
    + (-1 / 120) * (u ^ 5 * v ^ 0) * (s ^ 0 * L ^ 5)
    + (1 / 120) * (u ^ 5 * v ^ 1) * (s ^ 1 * L ^ 5)
    + (-1 / 240) * (u ^ 5 * v ^ 2) * (s ^ 2 * L ^ 5)
    + (1 / 720) * (u ^ 6 * v ^ 0) * (s ^ 0 * L ^ 6)
    + (-1 / 720) * (u ^ 6 * v ^ 1) * (s ^ 1 * L ^ 6)
    + (-1 / 5040) * (u ^ 7 * v ^ 0) * (s ^ 0 * L ^ 7))
      = HT + A * r2 + r1 * Real.exp (-q) := by
    rw [hHT, hr1, hr2, hA, hB, hp, hq]; ring
  rw [key]
  have hB9 : 0 ≤ B9 s := B9_nonneg hs
  have m0 := hmono8 hu0 hut hv0 hvt ht1 hl0 hl1 hle hs hp hq (c := (1 / 1440)) (by rw [abs_le]; constructor <;> norm_num) (i := 2) (j := 6) (by norm_num) (by norm_num) (by norm_num)
  have m1 := hmono8 hu0 hut hv0 hvt ht1 hl0 hl1 hle hs hp hq (c := (1 / 720)) (by rw [abs_le]; constructor <;> norm_num) (i := 3) (j := 5) (by norm_num) (by norm_num) (by norm_num)
  have m2 := hmono8 hu0 hut hv0 hvt ht1 hl0 hl1 hle hs hp hq (c := (-1 / 4320)) (by rw [abs_le]; constructor <;> norm_num) (i := 3) (j := 6) (by norm_num) (by norm_num) (by norm_num)
  have m3 := hmono8 hu0 hut hv0 hvt ht1 hl0 hl1 hle hs hp hq (c := (1 / 576)) (by rw [abs_le]; constructor <;> norm_num) (i := 4) (j := 4) (by norm_num) (by norm_num) (by norm_num)
  have m4 := hmono8 hu0 hut hv0 hvt ht1 hl0 hl1 hle hs hp hq (c := (-1 / 2880)) (by rw [abs_le]; constructor <;> norm_num) (i := 4) (j := 5) (by norm_num) (by norm_num) (by norm_num)
  have m5 := hmono8 hu0 hut hv0 hvt ht1 hl0 hl1 hle hs hp hq (c := (1 / 17280)) (by rw [abs_le]; constructor <;> norm_num) (i := 4) (j := 6) (by norm_num) (by norm_num) (by norm_num)
  have m6 := hmono8 hu0 hut hv0 hvt ht1 hl0 hl1 hle hs hp hq (c := (1 / 720)) (by rw [abs_le]; constructor <;> norm_num) (i := 5) (j := 3) (by norm_num) (by norm_num) (by norm_num)
  have m7 := hmono8 hu0 hut hv0 hvt ht1 hl0 hl1 hle hs hp hq (c := (-1 / 2880)) (by rw [abs_le]; constructor <;> norm_num) (i := 5) (j := 4) (by norm_num) (by norm_num) (by norm_num)
  have m8 := hmono8 hu0 hut hv0 hvt ht1 hl0 hl1 hle hs hp hq (c := (1 / 14400)) (by rw [abs_le]; constructor <;> norm_num) (i := 5) (j := 5) (by norm_num) (by norm_num) (by norm_num)
  have m9 := hmono8 hu0 hut hv0 hvt ht1 hl0 hl1 hle hs hp hq (c := (-1 / 86400)) (by rw [abs_le]; constructor <;> norm_num) (i := 5) (j := 6) (by norm_num) (by norm_num) (by norm_num)
  have m10 := hmono8 hu0 hut hv0 hvt ht1 hl0 hl1 hle hs hp hq (c := (1 / 1440)) (by rw [abs_le]; constructor <;> norm_num) (i := 6) (j := 2) (by norm_num) (by norm_num) (by norm_num)
  have m11 := hmono8 hu0 hut hv0 hvt ht1 hl0 hl1 hle hs hp hq (c := (-1 / 4320)) (by rw [abs_le]; constructor <;> norm_num) (i := 6) (j := 3) (by norm_num) (by norm_num) (by norm_num)
  have m12 := hmono8 hu0 hut hv0 hvt ht1 hl0 hl1 hle hs hp hq (c := (1 / 17280)) (by rw [abs_le]; constructor <;> norm_num) (i := 6) (j := 4) (by norm_num) (by norm_num) (by norm_num)
  have m13 := hmono8 hu0 hut hv0 hvt ht1 hl0 hl1 hle hs hp hq (c := (-1 / 86400)) (by rw [abs_le]; constructor <;> norm_num) (i := 6) (j := 5) (by norm_num) (by norm_num) (by norm_num)
  have m14 := hmono8 hu0 hut hv0 hvt ht1 hl0 hl1 hle hs hp hq (c := (1 / 518400)) (by rw [abs_le]; constructor <;> norm_num) (i := 6) (j := 6) (by norm_num) (by norm_num) (by norm_num)
  have m15 := hmono8 hu0 hut hv0 hvt ht1 hl0 hl1 hle hs hp hq (c := (1 / 5040)) (by rw [abs_le]; constructor <;> norm_num) (i := 7) (j := 1) (by norm_num) (by norm_num) (by norm_num)
  have m16 := hmono8 hu0 hut hv0 hvt ht1 hl0 hl1 hle hs hp hq (c := (-1 / 10080)) (by rw [abs_le]; constructor <;> norm_num) (i := 7) (j := 2) (by norm_num) (by norm_num) (by norm_num)
  have m17 := hmono8 hu0 hut hv0 hvt ht1 hl0 hl1 hle hs hp hq (c := (1 / 30240)) (by rw [abs_le]; constructor <;> norm_num) (i := 7) (j := 3) (by norm_num) (by norm_num) (by norm_num)
  have m18 := hmono8 hu0 hut hv0 hvt ht1 hl0 hl1 hle hs hp hq (c := (-1 / 120960)) (by rw [abs_le]; constructor <;> norm_num) (i := 7) (j := 4) (by norm_num) (by norm_num) (by norm_num)
  have m19 := hmono8 hu0 hut hv0 hvt ht1 hl0 hl1 hle hs hp hq (c := (1 / 604800)) (by rw [abs_le]; constructor <;> norm_num) (i := 7) (j := 5) (by norm_num) (by norm_num) (by norm_num)
  have m20 := hmono8 hu0 hut hv0 hvt ht1 hl0 hl1 hle hs hp hq (c := (-1 / 3628800)) (by rw [abs_le]; constructor <;> norm_num) (i := 7) (j := 6) (by norm_num) (by norm_num) (by norm_num)
  have c1 := abs_add_le ((1 / 1440) * (p ^ 2 * q ^ 6)) ((1 / 720) * (p ^ 3 * q ^ 5))
  have c2 := abs_add_le ((1 / 1440) * (p ^ 2 * q ^ 6) + (1 / 720) * (p ^ 3 * q ^ 5)) ((-1 / 4320) * (p ^ 3 * q ^ 6))
  have c3 := abs_add_le ((1 / 1440) * (p ^ 2 * q ^ 6) + (1 / 720) * (p ^ 3 * q ^ 5) + (-1 / 4320) * (p ^ 3 * q ^ 6)) ((1 / 576) * (p ^ 4 * q ^ 4))
  have c4 := abs_add_le ((1 / 1440) * (p ^ 2 * q ^ 6) + (1 / 720) * (p ^ 3 * q ^ 5) + (-1 / 4320) * (p ^ 3 * q ^ 6) + (1 / 576) * (p ^ 4 * q ^ 4)) ((-1 / 2880) * (p ^ 4 * q ^ 5))
  have c5 := abs_add_le ((1 / 1440) * (p ^ 2 * q ^ 6) + (1 / 720) * (p ^ 3 * q ^ 5) + (-1 / 4320) * (p ^ 3 * q ^ 6) + (1 / 576) * (p ^ 4 * q ^ 4) + (-1 / 2880) * (p ^ 4 * q ^ 5)) ((1 / 17280) * (p ^ 4 * q ^ 6))
  have c6 := abs_add_le ((1 / 1440) * (p ^ 2 * q ^ 6) + (1 / 720) * (p ^ 3 * q ^ 5) + (-1 / 4320) * (p ^ 3 * q ^ 6) + (1 / 576) * (p ^ 4 * q ^ 4) + (-1 / 2880) * (p ^ 4 * q ^ 5) + (1 / 17280) * (p ^ 4 * q ^ 6)) ((1 / 720) * (p ^ 5 * q ^ 3))
  have c7 := abs_add_le ((1 / 1440) * (p ^ 2 * q ^ 6) + (1 / 720) * (p ^ 3 * q ^ 5) + (-1 / 4320) * (p ^ 3 * q ^ 6) + (1 / 576) * (p ^ 4 * q ^ 4) + (-1 / 2880) * (p ^ 4 * q ^ 5) + (1 / 17280) * (p ^ 4 * q ^ 6) + (1 / 720) * (p ^ 5 * q ^ 3)) ((-1 / 2880) * (p ^ 5 * q ^ 4))
  have c8 := abs_add_le ((1 / 1440) * (p ^ 2 * q ^ 6) + (1 / 720) * (p ^ 3 * q ^ 5) + (-1 / 4320) * (p ^ 3 * q ^ 6) + (1 / 576) * (p ^ 4 * q ^ 4) + (-1 / 2880) * (p ^ 4 * q ^ 5) + (1 / 17280) * (p ^ 4 * q ^ 6) + (1 / 720) * (p ^ 5 * q ^ 3) + (-1 / 2880) * (p ^ 5 * q ^ 4)) ((1 / 14400) * (p ^ 5 * q ^ 5))
  have c9 := abs_add_le ((1 / 1440) * (p ^ 2 * q ^ 6) + (1 / 720) * (p ^ 3 * q ^ 5) + (-1 / 4320) * (p ^ 3 * q ^ 6) + (1 / 576) * (p ^ 4 * q ^ 4) + (-1 / 2880) * (p ^ 4 * q ^ 5) + (1 / 17280) * (p ^ 4 * q ^ 6) + (1 / 720) * (p ^ 5 * q ^ 3) + (-1 / 2880) * (p ^ 5 * q ^ 4) + (1 / 14400) * (p ^ 5 * q ^ 5)) ((-1 / 86400) * (p ^ 5 * q ^ 6))
  have c10 := abs_add_le ((1 / 1440) * (p ^ 2 * q ^ 6) + (1 / 720) * (p ^ 3 * q ^ 5) + (-1 / 4320) * (p ^ 3 * q ^ 6) + (1 / 576) * (p ^ 4 * q ^ 4) + (-1 / 2880) * (p ^ 4 * q ^ 5) + (1 / 17280) * (p ^ 4 * q ^ 6) + (1 / 720) * (p ^ 5 * q ^ 3) + (-1 / 2880) * (p ^ 5 * q ^ 4) + (1 / 14400) * (p ^ 5 * q ^ 5) + (-1 / 86400) * (p ^ 5 * q ^ 6)) ((1 / 1440) * (p ^ 6 * q ^ 2))
  have c11 := abs_add_le ((1 / 1440) * (p ^ 2 * q ^ 6) + (1 / 720) * (p ^ 3 * q ^ 5) + (-1 / 4320) * (p ^ 3 * q ^ 6) + (1 / 576) * (p ^ 4 * q ^ 4) + (-1 / 2880) * (p ^ 4 * q ^ 5) + (1 / 17280) * (p ^ 4 * q ^ 6) + (1 / 720) * (p ^ 5 * q ^ 3) + (-1 / 2880) * (p ^ 5 * q ^ 4) + (1 / 14400) * (p ^ 5 * q ^ 5) + (-1 / 86400) * (p ^ 5 * q ^ 6) + (1 / 1440) * (p ^ 6 * q ^ 2)) ((-1 / 4320) * (p ^ 6 * q ^ 3))
  have c12 := abs_add_le ((1 / 1440) * (p ^ 2 * q ^ 6) + (1 / 720) * (p ^ 3 * q ^ 5) + (-1 / 4320) * (p ^ 3 * q ^ 6) + (1 / 576) * (p ^ 4 * q ^ 4) + (-1 / 2880) * (p ^ 4 * q ^ 5) + (1 / 17280) * (p ^ 4 * q ^ 6) + (1 / 720) * (p ^ 5 * q ^ 3) + (-1 / 2880) * (p ^ 5 * q ^ 4) + (1 / 14400) * (p ^ 5 * q ^ 5) + (-1 / 86400) * (p ^ 5 * q ^ 6) + (1 / 1440) * (p ^ 6 * q ^ 2) + (-1 / 4320) * (p ^ 6 * q ^ 3)) ((1 / 17280) * (p ^ 6 * q ^ 4))
  have c13 := abs_add_le ((1 / 1440) * (p ^ 2 * q ^ 6) + (1 / 720) * (p ^ 3 * q ^ 5) + (-1 / 4320) * (p ^ 3 * q ^ 6) + (1 / 576) * (p ^ 4 * q ^ 4) + (-1 / 2880) * (p ^ 4 * q ^ 5) + (1 / 17280) * (p ^ 4 * q ^ 6) + (1 / 720) * (p ^ 5 * q ^ 3) + (-1 / 2880) * (p ^ 5 * q ^ 4) + (1 / 14400) * (p ^ 5 * q ^ 5) + (-1 / 86400) * (p ^ 5 * q ^ 6) + (1 / 1440) * (p ^ 6 * q ^ 2) + (-1 / 4320) * (p ^ 6 * q ^ 3) + (1 / 17280) * (p ^ 6 * q ^ 4)) ((-1 / 86400) * (p ^ 6 * q ^ 5))
  have c14 := abs_add_le ((1 / 1440) * (p ^ 2 * q ^ 6) + (1 / 720) * (p ^ 3 * q ^ 5) + (-1 / 4320) * (p ^ 3 * q ^ 6) + (1 / 576) * (p ^ 4 * q ^ 4) + (-1 / 2880) * (p ^ 4 * q ^ 5) + (1 / 17280) * (p ^ 4 * q ^ 6) + (1 / 720) * (p ^ 5 * q ^ 3) + (-1 / 2880) * (p ^ 5 * q ^ 4) + (1 / 14400) * (p ^ 5 * q ^ 5) + (-1 / 86400) * (p ^ 5 * q ^ 6) + (1 / 1440) * (p ^ 6 * q ^ 2) + (-1 / 4320) * (p ^ 6 * q ^ 3) + (1 / 17280) * (p ^ 6 * q ^ 4) + (-1 / 86400) * (p ^ 6 * q ^ 5)) ((1 / 518400) * (p ^ 6 * q ^ 6))
  have c15 := abs_add_le ((1 / 1440) * (p ^ 2 * q ^ 6) + (1 / 720) * (p ^ 3 * q ^ 5) + (-1 / 4320) * (p ^ 3 * q ^ 6) + (1 / 576) * (p ^ 4 * q ^ 4) + (-1 / 2880) * (p ^ 4 * q ^ 5) + (1 / 17280) * (p ^ 4 * q ^ 6) + (1 / 720) * (p ^ 5 * q ^ 3) + (-1 / 2880) * (p ^ 5 * q ^ 4) + (1 / 14400) * (p ^ 5 * q ^ 5) + (-1 / 86400) * (p ^ 5 * q ^ 6) + (1 / 1440) * (p ^ 6 * q ^ 2) + (-1 / 4320) * (p ^ 6 * q ^ 3) + (1 / 17280) * (p ^ 6 * q ^ 4) + (-1 / 86400) * (p ^ 6 * q ^ 5) + (1 / 518400) * (p ^ 6 * q ^ 6)) ((1 / 5040) * (p ^ 7 * q ^ 1))
  have c16 := abs_add_le ((1 / 1440) * (p ^ 2 * q ^ 6) + (1 / 720) * (p ^ 3 * q ^ 5) + (-1 / 4320) * (p ^ 3 * q ^ 6) + (1 / 576) * (p ^ 4 * q ^ 4) + (-1 / 2880) * (p ^ 4 * q ^ 5) + (1 / 17280) * (p ^ 4 * q ^ 6) + (1 / 720) * (p ^ 5 * q ^ 3) + (-1 / 2880) * (p ^ 5 * q ^ 4) + (1 / 14400) * (p ^ 5 * q ^ 5) + (-1 / 86400) * (p ^ 5 * q ^ 6) + (1 / 1440) * (p ^ 6 * q ^ 2) + (-1 / 4320) * (p ^ 6 * q ^ 3) + (1 / 17280) * (p ^ 6 * q ^ 4) + (-1 / 86400) * (p ^ 6 * q ^ 5) + (1 / 518400) * (p ^ 6 * q ^ 6) + (1 / 5040) * (p ^ 7 * q ^ 1)) ((-1 / 10080) * (p ^ 7 * q ^ 2))
  have c17 := abs_add_le ((1 / 1440) * (p ^ 2 * q ^ 6) + (1 / 720) * (p ^ 3 * q ^ 5) + (-1 / 4320) * (p ^ 3 * q ^ 6) + (1 / 576) * (p ^ 4 * q ^ 4) + (-1 / 2880) * (p ^ 4 * q ^ 5) + (1 / 17280) * (p ^ 4 * q ^ 6) + (1 / 720) * (p ^ 5 * q ^ 3) + (-1 / 2880) * (p ^ 5 * q ^ 4) + (1 / 14400) * (p ^ 5 * q ^ 5) + (-1 / 86400) * (p ^ 5 * q ^ 6) + (1 / 1440) * (p ^ 6 * q ^ 2) + (-1 / 4320) * (p ^ 6 * q ^ 3) + (1 / 17280) * (p ^ 6 * q ^ 4) + (-1 / 86400) * (p ^ 6 * q ^ 5) + (1 / 518400) * (p ^ 6 * q ^ 6) + (1 / 5040) * (p ^ 7 * q ^ 1) + (-1 / 10080) * (p ^ 7 * q ^ 2)) ((1 / 30240) * (p ^ 7 * q ^ 3))
  have c18 := abs_add_le ((1 / 1440) * (p ^ 2 * q ^ 6) + (1 / 720) * (p ^ 3 * q ^ 5) + (-1 / 4320) * (p ^ 3 * q ^ 6) + (1 / 576) * (p ^ 4 * q ^ 4) + (-1 / 2880) * (p ^ 4 * q ^ 5) + (1 / 17280) * (p ^ 4 * q ^ 6) + (1 / 720) * (p ^ 5 * q ^ 3) + (-1 / 2880) * (p ^ 5 * q ^ 4) + (1 / 14400) * (p ^ 5 * q ^ 5) + (-1 / 86400) * (p ^ 5 * q ^ 6) + (1 / 1440) * (p ^ 6 * q ^ 2) + (-1 / 4320) * (p ^ 6 * q ^ 3) + (1 / 17280) * (p ^ 6 * q ^ 4) + (-1 / 86400) * (p ^ 6 * q ^ 5) + (1 / 518400) * (p ^ 6 * q ^ 6) + (1 / 5040) * (p ^ 7 * q ^ 1) + (-1 / 10080) * (p ^ 7 * q ^ 2) + (1 / 30240) * (p ^ 7 * q ^ 3)) ((-1 / 120960) * (p ^ 7 * q ^ 4))
  have c19 := abs_add_le ((1 / 1440) * (p ^ 2 * q ^ 6) + (1 / 720) * (p ^ 3 * q ^ 5) + (-1 / 4320) * (p ^ 3 * q ^ 6) + (1 / 576) * (p ^ 4 * q ^ 4) + (-1 / 2880) * (p ^ 4 * q ^ 5) + (1 / 17280) * (p ^ 4 * q ^ 6) + (1 / 720) * (p ^ 5 * q ^ 3) + (-1 / 2880) * (p ^ 5 * q ^ 4) + (1 / 14400) * (p ^ 5 * q ^ 5) + (-1 / 86400) * (p ^ 5 * q ^ 6) + (1 / 1440) * (p ^ 6 * q ^ 2) + (-1 / 4320) * (p ^ 6 * q ^ 3) + (1 / 17280) * (p ^ 6 * q ^ 4) + (-1 / 86400) * (p ^ 6 * q ^ 5) + (1 / 518400) * (p ^ 6 * q ^ 6) + (1 / 5040) * (p ^ 7 * q ^ 1) + (-1 / 10080) * (p ^ 7 * q ^ 2) + (1 / 30240) * (p ^ 7 * q ^ 3) + (-1 / 120960) * (p ^ 7 * q ^ 4)) ((1 / 604800) * (p ^ 7 * q ^ 5))
  have c20 := abs_add_le ((1 / 1440) * (p ^ 2 * q ^ 6) + (1 / 720) * (p ^ 3 * q ^ 5) + (-1 / 4320) * (p ^ 3 * q ^ 6) + (1 / 576) * (p ^ 4 * q ^ 4) + (-1 / 2880) * (p ^ 4 * q ^ 5) + (1 / 17280) * (p ^ 4 * q ^ 6) + (1 / 720) * (p ^ 5 * q ^ 3) + (-1 / 2880) * (p ^ 5 * q ^ 4) + (1 / 14400) * (p ^ 5 * q ^ 5) + (-1 / 86400) * (p ^ 5 * q ^ 6) + (1 / 1440) * (p ^ 6 * q ^ 2) + (-1 / 4320) * (p ^ 6 * q ^ 3) + (1 / 17280) * (p ^ 6 * q ^ 4) + (-1 / 86400) * (p ^ 6 * q ^ 5) + (1 / 518400) * (p ^ 6 * q ^ 6) + (1 / 5040) * (p ^ 7 * q ^ 1) + (-1 / 10080) * (p ^ 7 * q ^ 2) + (1 / 30240) * (p ^ 7 * q ^ 3) + (-1 / 120960) * (p ^ 7 * q ^ 4) + (1 / 604800) * (p ^ 7 * q ^ 5)) ((-1 / 3628800) * (p ^ 7 * q ^ 6))
  have bHT : |HT| ≤ 21 * t ^ 8 * B9 s := by rw [hHT]; linarith
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
    rw [abs_le, hA]; constructor <;> linarith
  have bAr : |A * r2| ≤ 8 * t ^ 8 * B9 s := by
    rw [abs_mul]
    have h1 : |A| * |r2| ≤ 2 * p * (4 * q ^ 7) :=
      mul_le_mul hAb hr2b (abs_nonneg _) (by positivity)
    have h2 : |(1 : ℝ) * (p ^ 1 * q ^ 7)| ≤ t ^ 8 * B9 s :=
      hmono8 hu0 hut hv0 hvt ht1 hl0 hl1 hle hs hp hq (by norm_num) (by norm_num) (by norm_num) (by norm_num)
    rw [one_mul, pow_one, abs_of_nonneg (by positivity)] at h2
    calc |A| * |r2| ≤ 2 * p * (4 * q ^ 7) := h1
      _ = 8 * (p * q ^ 7) := by ring
      _ ≤ 8 * (t ^ 8 * B9 s) := by linarith
      _ = _ := by ring
  have br1 : |r1 * Real.exp (-q)| ≤ t ^ 8 * B9 s := by
    rw [abs_mul, abs_of_pos (Real.exp_pos _)]
    have h1 : Real.exp (-q) ≤ 1 := Real.exp_le_one_iff.2 (by linarith)
    have h2 : |(1 : ℝ) * (p ^ 8 * q ^ 0)| ≤ t ^ 8 * B9 s :=
      hmono8 hu0 hut hv0 hvt ht1 hl0 hl1 hle hs hp hq (by norm_num) (by norm_num) (by norm_num) (by norm_num)
    rw [one_mul, pow_zero, mul_one, abs_of_nonneg (by positivity)] at h2
    calc |r1| * Real.exp (-q) ≤ p ^ 8 * 1 :=
          mul_le_mul hr1b h1 (Real.exp_pos _).le (by positivity)
      _ ≤ _ := by linarith
  have t1 := abs_add_le (HT + A * r2) (r1 * Real.exp (-q))
  have t2 := abs_add_le HT (A * r2)
  have : 0 ≤ t ^ 8 * B9 s := by positivity
  linarith

/-- the integrated expansion of `(e^{-uℓ} - 1) e^{-vs}` to order `7` -/
lemma int_expand7 {u v t : ℝ} (hu0 : 0 ≤ u) (hut : u ≤ t) (hv0 : 0 ≤ v) (hvt : v ≤ t) (ht1 : t ≤ 1) :
    IntegrableOn (fun s => (Real.exp (-u * Y04.ell s) - 1) * Real.exp (-v * s)) (Set.Ioi 0) ∧
    |(∫ s in Set.Ioi (0 : ℝ), (Real.exp (-u * Y04.ell s) - 1) * Real.exp (-v * s))
      - ((-1) * (u ^ 1 * v ^ 0) * Y04.I 0 1 + (1) * (u ^ 1 * v ^ 1) * Y04.I 1 1 + (-1 / 2) * (u ^ 1 * v ^ 2) * Y04.I 2 1 + (1 / 6) * (u ^ 1 * v ^ 3) * Y04.I 3 1 + (-1 / 24) * (u ^ 1 * v ^ 4) * Y04.I 4 1 + (1 / 120) * (u ^ 1 * v ^ 5) * Y04.I 5 1 + (-1 / 720) * (u ^ 1 * v ^ 6) * Y04.I 6 1 + (1 / 2) * (u ^ 2 * v ^ 0) * Y04.I 0 2 + (-1 / 2) * (u ^ 2 * v ^ 1) * Y04.I 1 2 + (1 / 4) * (u ^ 2 * v ^ 2) * Y04.I 2 2 + (-1 / 12) * (u ^ 2 * v ^ 3) * Y04.I 3 2 + (1 / 48) * (u ^ 2 * v ^ 4) * Y04.I 4 2 + (-1 / 240) * (u ^ 2 * v ^ 5) * Y04.I 5 2 + (-1 / 6) * (u ^ 3 * v ^ 0) * Y04.I 0 3 + (1 / 6) * (u ^ 3 * v ^ 1) * Y04.I 1 3 + (-1 / 12) * (u ^ 3 * v ^ 2) * Y04.I 2 3 + (1 / 36) * (u ^ 3 * v ^ 3) * Y04.I 3 3 + (-1 / 144) * (u ^ 3 * v ^ 4) * Y04.I 4 3 + (1 / 24) * (u ^ 4 * v ^ 0) * Y04.I 0 4 + (-1 / 24) * (u ^ 4 * v ^ 1) * Y04.I 1 4 + (1 / 48) * (u ^ 4 * v ^ 2) * Y04.I 2 4 + (-1 / 144) * (u ^ 4 * v ^ 3) * Y04.I 3 4 + (-1 / 120) * (u ^ 5 * v ^ 0) * Y04.I 0 5 + (1 / 120) * (u ^ 5 * v ^ 1) * Y04.I 1 5 + (-1 / 240) * (u ^ 5 * v ^ 2) * Y04.I 2 5 + (1 / 720) * (u ^ 6 * v ^ 0) * Y04.I 0 6 + (-1 / 720) * (u ^ 6 * v ^ 1) * Y04.I 1 6 + (-1 / 5040) * (u ^ 7 * v ^ 0) * Y04.I 0 7)|
      ≤ 30 * t ^ 8 * ∫ s in Set.Ioi (0 : ℝ), B9 s := by
  have hint : IntegrableOn (fun s => (Real.exp (-u * Y04.ell s) - 1) * Real.exp (-v * s)) (Set.Ioi 0) :=
    (Z02.int_expand hu0 hut hv0 hvt ht1).1
  refine ⟨hint, ?_⟩
  obtain ⟨hmi, hmv⟩ := int28 ((-1) * (u ^ 1 * v ^ 0)) ((1) * (u ^ 1 * v ^ 1)) ((-1 / 2) * (u ^ 1 * v ^ 2)) ((1 / 6) * (u ^ 1 * v ^ 3)) ((-1 / 24) * (u ^ 1 * v ^ 4)) ((1 / 120) * (u ^ 1 * v ^ 5)) ((-1 / 720) * (u ^ 1 * v ^ 6)) ((1 / 2) * (u ^ 2 * v ^ 0)) ((-1 / 2) * (u ^ 2 * v ^ 1)) ((1 / 4) * (u ^ 2 * v ^ 2)) ((-1 / 12) * (u ^ 2 * v ^ 3)) ((1 / 48) * (u ^ 2 * v ^ 4)) ((-1 / 240) * (u ^ 2 * v ^ 5)) ((-1 / 6) * (u ^ 3 * v ^ 0)) ((1 / 6) * (u ^ 3 * v ^ 1)) ((-1 / 12) * (u ^ 3 * v ^ 2)) ((1 / 36) * (u ^ 3 * v ^ 3)) ((-1 / 144) * (u ^ 3 * v ^ 4)) ((1 / 24) * (u ^ 4 * v ^ 0)) ((-1 / 24) * (u ^ 4 * v ^ 1)) ((1 / 48) * (u ^ 4 * v ^ 2)) ((-1 / 144) * (u ^ 4 * v ^ 3)) ((-1 / 120) * (u ^ 5 * v ^ 0)) ((1 / 120) * (u ^ 5 * v ^ 1)) ((-1 / 240) * (u ^ 5 * v ^ 2)) ((1 / 720) * (u ^ 6 * v ^ 0)) ((-1 / 720) * (u ^ 6 * v ^ 1)) ((-1 / 5040) * (u ^ 7 * v ^ 0))
    (Y04.I_int 0 1 le_rfl) (Y04.I_int 1 1 le_rfl) (Y04.I_int 2 1 le_rfl) (Y04.I_int 3 1 le_rfl) (Y04.I_int 4 1 le_rfl) (Y04.I_int 5 1 le_rfl) (Y04.I_int 6 1 le_rfl) (Y04.I_int 0 2 (by norm_num)) (Y04.I_int 1 2 (by norm_num)) (Y04.I_int 2 2 (by norm_num)) (Y04.I_int 3 2 (by norm_num)) (Y04.I_int 4 2 (by norm_num)) (Y04.I_int 5 2 (by norm_num)) (Y04.I_int 0 3 (by norm_num)) (Y04.I_int 1 3 (by norm_num)) (Y04.I_int 2 3 (by norm_num)) (Y04.I_int 3 3 (by norm_num)) (Y04.I_int 4 3 (by norm_num)) (Y04.I_int 0 4 (by norm_num)) (Y04.I_int 1 4 (by norm_num)) (Y04.I_int 2 4 (by norm_num)) (Y04.I_int 3 4 (by norm_num)) (Y04.I_int 0 5 (by norm_num)) (Y04.I_int 1 5 (by norm_num)) (Y04.I_int 2 5 (by norm_num)) (Y04.I_int 0 6 (by norm_num)) (Y04.I_int 1 6 (by norm_num)) (Y04.I_int 0 7 (by norm_num))
  unfold Y04.I
  rw [← hmv, ← integral_sub hint hmi]
  have hB := (B9_int.const_mul 30).const_mul (t ^ 8)
  have := norm_integral_le_of_norm_le hB
    ((ae_restrict_iff' measurableSet_Ioi).2 (Filter.Eventually.of_forall
      (fun s (hs : 0 < s) => by
        rw [Real.norm_eq_abs]
        have := pt7 hu0 hut hv0 hvt ht1 hs.le
        unfold P7 at this
        calc _ = _ := rfl
          _ ≤ _ := this
          _ = t ^ 8 * (30 * B9 s) := by ring)))
  rw [Real.norm_eq_abs, integral_const_mul, integral_const_mul] at this
  calc _ ≤ _ := this
    _ = _ := by ring

/-- seven coefficients vanish if the septic is `O(t⁸)` on `(0, d]` -/
lemma septic_zero {d1 d2 d3 d4 d5 d6 d7 K d : ℝ} (hd : 0 < d) (hd1 : d ≤ 1) (hK : 0 ≤ K)
    (h : ∀ t : ℝ, 0 < t → t ≤ d →
      |d1 * t + d2 * t ^ 2 + d3 * t ^ 3 + d4 * t ^ 4 + d5 * t ^ 5 + d6 * t ^ 6 + d7 * t ^ 7| ≤ K * t ^ 8) :
    d1 = 0 ∧ d2 = 0 ∧ d3 = 0 ∧ d4 = 0 ∧ d5 = 0 ∧ d6 = 0 ∧ d7 = 0 := by
  have h1 : d1 = 0 := by
    refine Y04.eq_zero_of_le_lin hd (K := |d2| + |d3| + |d4| + |d5| + |d6| + |d7| + K)
      (by have := abs_nonneg d2; have := abs_nonneg d3; have := abs_nonneg d4; have := abs_nonneg d5; have := abs_nonneg d6; have := abs_nonneg d7
          linarith)
      (fun t ht0 htd => ?_)
    have ht1 : t ≤ 1 := htd.trans hd1
    have := h t ht0 htd
    set E := d1 * t + d2 * t ^ 2 + d3 * t ^ 3 + d4 * t ^ 4 + d5 * t ^ 5 + d6 * t ^ 6 + d7 * t ^ 7 with hE
    have e : d1 = E / t - d2 * t - d3 * t ^ 2 - d4 * t ^ 3 - d5 * t ^ 4 - d6 * t ^ 5 - d7 * t ^ 6 := by
      rw [hE]; field_simp; ring
    rw [e]
    have ht2 : t ^ 2 ≤ t := pow_le_of_le_one ht0.le ht1 (by norm_num)
    have ht3 : t ^ 3 ≤ t := pow_le_of_le_one ht0.le ht1 (by norm_num)
    have ht4 : t ^ 4 ≤ t := pow_le_of_le_one ht0.le ht1 (by norm_num)
    have ht5 : t ^ 5 ≤ t := pow_le_of_le_one ht0.le ht1 (by norm_num)
    have ht6 : t ^ 6 ≤ t := pow_le_of_le_one ht0.le ht1 (by norm_num)
    have ht7 : t ^ 7 ≤ t := pow_le_of_le_one ht0.le ht1 (by norm_num)
    have b1 : |E / t| ≤ K * t ^ 7 := by
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
    have t7 : K * t ^ 7 ≤ K * t := mul_le_mul_of_nonneg_left ht7 hK
    have := abs_sub (E / t - d2 * t - d3 * t ^ 2 - d4 * t ^ 3 - d5 * t ^ 4 - d6 * t ^ 5) (d7 * t ^ 6)
    have := abs_sub (E / t - d2 * t - d3 * t ^ 2 - d4 * t ^ 3 - d5 * t ^ 4) (d6 * t ^ 5)
    have := abs_sub (E / t - d2 * t - d3 * t ^ 2 - d4 * t ^ 3) (d5 * t ^ 4)
    have := abs_sub (E / t - d2 * t - d3 * t ^ 2) (d4 * t ^ 3)
    have := abs_sub (E / t - d2 * t) (d3 * t ^ 2)
    have := abs_sub (E / t) (d2 * t)
    have : (|d2| + |d3| + |d4| + |d5| + |d6| + |d7| + K) * t
        = |d2| * t + |d3| * t + |d4| * t + |d5| * t + |d6| * t + |d7| * t + K * t := by ring
    linarith
  subst h1
  obtain ⟨a, b, c, e, f, g⟩ := Q03.sextic_zero hd hd1 hK (d1 := d2) (d2 := d3) (d3 := d4) (d4 := d5) (d5 := d6) (d6 := d7)
    (fun t ht0 htd => by
    have := h t ht0 htd
    have e : 0 * t + d2 * t ^ 2 + d3 * t ^ 3 + d4 * t ^ 4 + d5 * t ^ 5 + d6 * t ^ 6 + d7 * t ^ 7
        = t * (d2 * t + d3 * t ^ 2 + d4 * t ^ 3 + d5 * t ^ 4 + d6 * t ^ 5 + d7 * t ^ 6) := by ring
    rw [e, abs_mul, abs_of_pos ht0, show K * t ^ 8 = t * (K * t ^ 7) by ring] at this
    exact le_of_mul_le_mul_left this ht0)
  exact ⟨rfl, a, b, c, e, f, g⟩

lemma series_trunc8 {α β γ x L : ℝ} (hα : |α| ≤ 3 / 2) (hβ : |β| ≤ 3 / 2) (hγ : |γ| ≤ 3 / 2)
    (hx0 : 0 ≤ x) (hx : x ≤ 1 / 3)
    (hL : HasSum (fun k : ℕ => (-1 : ℝ) ^ (k + 2) * Y04.zn (k + 2) / ((k + 2 : ℕ) : ℝ) *
      ((α * x) ^ (k + 2) + (β * x) ^ (k + 2) - (γ * x) ^ (k + 2))) L) :
    |L - (Y04.cf α β γ 0 * x ^ 2 + Y04.cf α β γ 1 * x ^ 3 + Y04.cf α β γ 2 * x ^ 4 + Y04.cf α β γ 3 * x ^ 5 + Y04.cf α β γ 4 * x ^ 6 + Y04.cf α β γ 5 * x ^ 7 + Y04.cf α β γ 6 * x ^ 8)| ≤ 231 * x ^ 9 := by
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
  have hg := (hasSum_nat_add_iff' 7).mpr hL
  have hgeo : HasSum (fun k : ℕ => 3 * r ^ 9 * r ^ k) (3 * r ^ 9 * (1 - r)⁻¹) :=
    (hasSum_geometric_of_lt_one hr0 (by linarith)).mul_left _
  have hb1 : ∀ k : ℕ, Y04.cf α β γ (k + 7) * x ^ (k + 7 + 2) ≤ 3 * r ^ 9 * r ^ k := by
    intro k
    have := (abs_le.mp (hbd (k + 7))).2
    calc _ ≤ 3 * r ^ (k + 7 + 2) := this
      _ = _ := by ring
  have hb2 : ∀ k : ℕ, -(Y04.cf α β γ (k + 7) * x ^ (k + 7 + 2)) ≤ 3 * r ^ 9 * r ^ k := by
    intro k
    have := (abs_le.mp (hbd (k + 7))).1
    calc _ ≤ 3 * r ^ (k + 7 + 2) := by linarith
      _ = _ := by ring
  have u1 := hasSum_le hb1 hg hgeo
  have u2 := hasSum_le hb2 hg.neg hgeo
  have hinv : 3 * r ^ 9 * (1 - r)⁻¹ ≤ 6 * r ^ 9 := by
    rw [mul_inv_le_iff₀ (by linarith)]
    nlinarith [pow_nonneg hr0 9]
  have hr9 : 6 * r ^ 9 ≤ 231 * x ^ 9 := by
    rw [hr, mul_pow]
    nlinarith [pow_nonneg hx0 9]
  simp only [Finset.sum_range_succ, Finset.sum_range_zero, zero_add] at u1 u2
  rw [abs_le]
  constructor <;> norm_num at u1 u2 ⊢ <;> linarith

/-- chunk 0 of the order-`≥ 9` part of `P + P²/2 + P³/6 + P⁴/24`, divided by `x⁹` -/
noncomputable def Wc0 (a0 a1 a2 a3 a4 a5 a6 x : ℝ) : ℝ :=
  (1) * (a2 * a3)
    + (1) * (a1 * a4)
    + (1 / 6) * (a1 ^ 3)
    + (1) * (a0 * a5)
    + (1) * (a0 * a1 * a2)
    + (1 / 2) * (a0 ^ 2 * a3)
    + (1 / 6) * (a0 ^ 3 * a1)
    + (1 / 2) * (a3 ^ 2 * x)
    + (1) * (a2 * a4 * x)
    + (1) * (a1 * a5 * x)
    + (1 / 2) * (a1 ^ 2 * a2 * x)
    + (1) * (a0 * a6 * x)
    + (1 / 2) * (a0 * a2 ^ 2 * x)
    + (1) * (a0 * a1 * a3 * x)
    + (1 / 2) * (a0 ^ 2 * a4 * x)
    + (1 / 4) * (a0 ^ 2 * a1 ^ 2 * x)
    + (1 / 6) * (a0 ^ 3 * a2 * x)
    + (1) * (a3 * a4 * x ^ 2)
    + (1) * (a2 * a5 * x ^ 2)
    + (1) * (a1 * a6 * x ^ 2)
    + (1 / 2) * (a1 * a2 ^ 2 * x ^ 2)
    + (1 / 2) * (a1 ^ 2 * a3 * x ^ 2)
    + (1) * (a0 * a2 * a3 * x ^ 2)
    + (1) * (a0 * a1 * a4 * x ^ 2)
    + (1 / 6) * (a0 * a1 ^ 3 * x ^ 2)
    + (1 / 2) * (a0 ^ 2 * a5 * x ^ 2)
    + (1 / 2) * (a0 ^ 2 * a1 * a2 * x ^ 2)
    + (1 / 6) * (a0 ^ 3 * a3 * x ^ 2)
    + (1 / 2) * (a4 ^ 2 * x ^ 3)
    + (1) * (a3 * a5 * x ^ 3)
    + (1) * (a2 * a6 * x ^ 3)
    + (1 / 6) * (a2 ^ 3 * x ^ 3)
    + (1) * (a1 * a2 * a3 * x ^ 3)
    + (1 / 2) * (a1 ^ 2 * a4 * x ^ 3)
    + (1 / 24) * (a1 ^ 4 * x ^ 3)
    + (1 / 2) * (a0 * a3 ^ 2 * x ^ 3)
    + (1) * (a0 * a2 * a4 * x ^ 3)
    + (1) * (a0 * a1 * a5 * x ^ 3)
    + (1 / 2) * (a0 * a1 ^ 2 * a2 * x ^ 3)
    + (1 / 2) * (a0 ^ 2 * a6 * x ^ 3)
/-- chunk 1 of the order-`≥ 9` part of `P + P²/2 + P³/6 + P⁴/24`, divided by `x⁹` -/
noncomputable def Wc1 (a0 a1 a2 a3 a4 a5 a6 x : ℝ) : ℝ :=
  (1 / 4) * (a0 ^ 2 * a2 ^ 2 * x ^ 3)
    + (1 / 2) * (a0 ^ 2 * a1 * a3 * x ^ 3)
    + (1 / 6) * (a0 ^ 3 * a4 * x ^ 3)
    + (1) * (a4 * a5 * x ^ 4)
    + (1) * (a3 * a6 * x ^ 4)
    + (1 / 2) * (a2 ^ 2 * a3 * x ^ 4)
    + (1 / 2) * (a1 * a3 ^ 2 * x ^ 4)
    + (1) * (a1 * a2 * a4 * x ^ 4)
    + (1 / 2) * (a1 ^ 2 * a5 * x ^ 4)
    + (1 / 6) * (a1 ^ 3 * a2 * x ^ 4)
    + (1) * (a0 * a3 * a4 * x ^ 4)
    + (1) * (a0 * a2 * a5 * x ^ 4)
    + (1) * (a0 * a1 * a6 * x ^ 4)
    + (1 / 2) * (a0 * a1 * a2 ^ 2 * x ^ 4)
    + (1 / 2) * (a0 * a1 ^ 2 * a3 * x ^ 4)
    + (1 / 2) * (a0 ^ 2 * a2 * a3 * x ^ 4)
    + (1 / 2) * (a0 ^ 2 * a1 * a4 * x ^ 4)
    + (1 / 6) * (a0 ^ 3 * a5 * x ^ 4)
    + (1 / 2) * (a5 ^ 2 * x ^ 5)
    + (1) * (a4 * a6 * x ^ 5)
    + (1 / 2) * (a2 * a3 ^ 2 * x ^ 5)
    + (1 / 2) * (a2 ^ 2 * a4 * x ^ 5)
    + (1) * (a1 * a3 * a4 * x ^ 5)
    + (1) * (a1 * a2 * a5 * x ^ 5)
    + (1 / 2) * (a1 ^ 2 * a6 * x ^ 5)
    + (1 / 4) * (a1 ^ 2 * a2 ^ 2 * x ^ 5)
    + (1 / 6) * (a1 ^ 3 * a3 * x ^ 5)
    + (1 / 2) * (a0 * a4 ^ 2 * x ^ 5)
    + (1) * (a0 * a3 * a5 * x ^ 5)
    + (1) * (a0 * a2 * a6 * x ^ 5)
    + (1 / 6) * (a0 * a2 ^ 3 * x ^ 5)
    + (1) * (a0 * a1 * a2 * a3 * x ^ 5)
    + (1 / 2) * (a0 * a1 ^ 2 * a4 * x ^ 5)
    + (1 / 4) * (a0 ^ 2 * a3 ^ 2 * x ^ 5)
    + (1 / 2) * (a0 ^ 2 * a2 * a4 * x ^ 5)
    + (1 / 2) * (a0 ^ 2 * a1 * a5 * x ^ 5)
    + (1 / 6) * (a0 ^ 3 * a6 * x ^ 5)
    + (1) * (a5 * a6 * x ^ 6)
    + (1 / 6) * (a3 ^ 3 * x ^ 6)
    + (1) * (a2 * a3 * a4 * x ^ 6)
/-- chunk 2 of the order-`≥ 9` part of `P + P²/2 + P³/6 + P⁴/24`, divided by `x⁹` -/
noncomputable def Wc2 (a0 a1 a2 a3 a4 a5 a6 x : ℝ) : ℝ :=
  (1 / 2) * (a2 ^ 2 * a5 * x ^ 6)
    + (1 / 2) * (a1 * a4 ^ 2 * x ^ 6)
    + (1) * (a1 * a3 * a5 * x ^ 6)
    + (1) * (a1 * a2 * a6 * x ^ 6)
    + (1 / 6) * (a1 * a2 ^ 3 * x ^ 6)
    + (1 / 2) * (a1 ^ 2 * a2 * a3 * x ^ 6)
    + (1 / 6) * (a1 ^ 3 * a4 * x ^ 6)
    + (1) * (a0 * a4 * a5 * x ^ 6)
    + (1) * (a0 * a3 * a6 * x ^ 6)
    + (1 / 2) * (a0 * a2 ^ 2 * a3 * x ^ 6)
    + (1 / 2) * (a0 * a1 * a3 ^ 2 * x ^ 6)
    + (1) * (a0 * a1 * a2 * a4 * x ^ 6)
    + (1 / 2) * (a0 * a1 ^ 2 * a5 * x ^ 6)
    + (1 / 2) * (a0 ^ 2 * a3 * a4 * x ^ 6)
    + (1 / 2) * (a0 ^ 2 * a2 * a5 * x ^ 6)
    + (1 / 2) * (a0 ^ 2 * a1 * a6 * x ^ 6)
    + (1 / 2) * (a6 ^ 2 * x ^ 7)
    + (1 / 2) * (a3 ^ 2 * a4 * x ^ 7)
    + (1 / 2) * (a2 * a4 ^ 2 * x ^ 7)
    + (1) * (a2 * a3 * a5 * x ^ 7)
    + (1 / 2) * (a2 ^ 2 * a6 * x ^ 7)
    + (1 / 24) * (a2 ^ 4 * x ^ 7)
    + (1) * (a1 * a4 * a5 * x ^ 7)
    + (1) * (a1 * a3 * a6 * x ^ 7)
    + (1 / 2) * (a1 * a2 ^ 2 * a3 * x ^ 7)
    + (1 / 4) * (a1 ^ 2 * a3 ^ 2 * x ^ 7)
    + (1 / 2) * (a1 ^ 2 * a2 * a4 * x ^ 7)
    + (1 / 6) * (a1 ^ 3 * a5 * x ^ 7)
    + (1 / 2) * (a0 * a5 ^ 2 * x ^ 7)
    + (1) * (a0 * a4 * a6 * x ^ 7)
    + (1 / 2) * (a0 * a2 * a3 ^ 2 * x ^ 7)
    + (1 / 2) * (a0 * a2 ^ 2 * a4 * x ^ 7)
    + (1) * (a0 * a1 * a3 * a4 * x ^ 7)
    + (1) * (a0 * a1 * a2 * a5 * x ^ 7)
    + (1 / 2) * (a0 * a1 ^ 2 * a6 * x ^ 7)
    + (1 / 4) * (a0 ^ 2 * a4 ^ 2 * x ^ 7)
    + (1 / 2) * (a0 ^ 2 * a3 * a5 * x ^ 7)
    + (1 / 2) * (a0 ^ 2 * a2 * a6 * x ^ 7)
    + (1 / 2) * (a3 * a4 ^ 2 * x ^ 8)
    + (1 / 2) * (a3 ^ 2 * a5 * x ^ 8)
/-- chunk 3 of the order-`≥ 9` part of `P + P²/2 + P³/6 + P⁴/24`, divided by `x⁹` -/
noncomputable def Wc3 (a0 a1 a2 a3 a4 a5 a6 x : ℝ) : ℝ :=
  (1) * (a2 * a4 * a5 * x ^ 8)
    + (1) * (a2 * a3 * a6 * x ^ 8)
    + (1 / 6) * (a2 ^ 3 * a3 * x ^ 8)
    + (1 / 2) * (a1 * a5 ^ 2 * x ^ 8)
    + (1) * (a1 * a4 * a6 * x ^ 8)
    + (1 / 2) * (a1 * a2 * a3 ^ 2 * x ^ 8)
    + (1 / 2) * (a1 * a2 ^ 2 * a4 * x ^ 8)
    + (1 / 2) * (a1 ^ 2 * a3 * a4 * x ^ 8)
    + (1 / 2) * (a1 ^ 2 * a2 * a5 * x ^ 8)
    + (1 / 6) * (a1 ^ 3 * a6 * x ^ 8)
    + (1) * (a0 * a5 * a6 * x ^ 8)
    + (1 / 6) * (a0 * a3 ^ 3 * x ^ 8)
    + (1) * (a0 * a2 * a3 * a4 * x ^ 8)
    + (1 / 2) * (a0 * a2 ^ 2 * a5 * x ^ 8)
    + (1 / 2) * (a0 * a1 * a4 ^ 2 * x ^ 8)
    + (1) * (a0 * a1 * a3 * a5 * x ^ 8)
    + (1) * (a0 * a1 * a2 * a6 * x ^ 8)
    + (1 / 2) * (a0 ^ 2 * a4 * a5 * x ^ 8)
    + (1 / 2) * (a0 ^ 2 * a3 * a6 * x ^ 8)
    + (1 / 6) * (a4 ^ 3 * x ^ 9)
    + (1) * (a3 * a4 * a5 * x ^ 9)
    + (1 / 2) * (a3 ^ 2 * a6 * x ^ 9)
    + (1 / 2) * (a2 * a5 ^ 2 * x ^ 9)
    + (1) * (a2 * a4 * a6 * x ^ 9)
    + (1 / 4) * (a2 ^ 2 * a3 ^ 2 * x ^ 9)
    + (1 / 6) * (a2 ^ 3 * a4 * x ^ 9)
    + (1) * (a1 * a5 * a6 * x ^ 9)
    + (1 / 6) * (a1 * a3 ^ 3 * x ^ 9)
    + (1) * (a1 * a2 * a3 * a4 * x ^ 9)
    + (1 / 2) * (a1 * a2 ^ 2 * a5 * x ^ 9)
    + (1 / 4) * (a1 ^ 2 * a4 ^ 2 * x ^ 9)
    + (1 / 2) * (a1 ^ 2 * a3 * a5 * x ^ 9)
    + (1 / 2) * (a1 ^ 2 * a2 * a6 * x ^ 9)
    + (1 / 2) * (a0 * a6 ^ 2 * x ^ 9)
    + (1 / 2) * (a0 * a3 ^ 2 * a4 * x ^ 9)
    + (1 / 2) * (a0 * a2 * a4 ^ 2 * x ^ 9)
    + (1) * (a0 * a2 * a3 * a5 * x ^ 9)
    + (1 / 2) * (a0 * a2 ^ 2 * a6 * x ^ 9)
    + (1) * (a0 * a1 * a4 * a5 * x ^ 9)
    + (1) * (a0 * a1 * a3 * a6 * x ^ 9)
/-- chunk 4 of the order-`≥ 9` part of `P + P²/2 + P³/6 + P⁴/24`, divided by `x⁹` -/
noncomputable def Wc4 (a0 a1 a2 a3 a4 a5 a6 x : ℝ) : ℝ :=
  (1 / 4) * (a0 ^ 2 * a5 ^ 2 * x ^ 9)
    + (1 / 2) * (a0 ^ 2 * a4 * a6 * x ^ 9)
    + (1 / 2) * (a4 ^ 2 * a5 * x ^ 10)
    + (1 / 2) * (a3 * a5 ^ 2 * x ^ 10)
    + (1) * (a3 * a4 * a6 * x ^ 10)
    + (1) * (a2 * a5 * a6 * x ^ 10)
    + (1 / 6) * (a2 * a3 ^ 3 * x ^ 10)
    + (1 / 2) * (a2 ^ 2 * a3 * a4 * x ^ 10)
    + (1 / 6) * (a2 ^ 3 * a5 * x ^ 10)
    + (1 / 2) * (a1 * a6 ^ 2 * x ^ 10)
    + (1 / 2) * (a1 * a3 ^ 2 * a4 * x ^ 10)
    + (1 / 2) * (a1 * a2 * a4 ^ 2 * x ^ 10)
    + (1) * (a1 * a2 * a3 * a5 * x ^ 10)
    + (1 / 2) * (a1 * a2 ^ 2 * a6 * x ^ 10)
    + (1 / 2) * (a1 ^ 2 * a4 * a5 * x ^ 10)
    + (1 / 2) * (a1 ^ 2 * a3 * a6 * x ^ 10)
    + (1 / 2) * (a0 * a3 * a4 ^ 2 * x ^ 10)
    + (1 / 2) * (a0 * a3 ^ 2 * a5 * x ^ 10)
    + (1) * (a0 * a2 * a4 * a5 * x ^ 10)
    + (1) * (a0 * a2 * a3 * a6 * x ^ 10)
    + (1 / 2) * (a0 * a1 * a5 ^ 2 * x ^ 10)
    + (1) * (a0 * a1 * a4 * a6 * x ^ 10)
    + (1 / 2) * (a0 ^ 2 * a5 * a6 * x ^ 10)
    + (1 / 2) * (a4 * a5 ^ 2 * x ^ 11)
    + (1 / 2) * (a4 ^ 2 * a6 * x ^ 11)
    + (1) * (a3 * a5 * a6 * x ^ 11)
    + (1 / 24) * (a3 ^ 4 * x ^ 11)
    + (1 / 2) * (a2 * a6 ^ 2 * x ^ 11)
    + (1 / 2) * (a2 * a3 ^ 2 * a4 * x ^ 11)
    + (1 / 4) * (a2 ^ 2 * a4 ^ 2 * x ^ 11)
    + (1 / 2) * (a2 ^ 2 * a3 * a5 * x ^ 11)
    + (1 / 6) * (a2 ^ 3 * a6 * x ^ 11)
    + (1 / 2) * (a1 * a3 * a4 ^ 2 * x ^ 11)
    + (1 / 2) * (a1 * a3 ^ 2 * a5 * x ^ 11)
    + (1) * (a1 * a2 * a4 * a5 * x ^ 11)
    + (1) * (a1 * a2 * a3 * a6 * x ^ 11)
    + (1 / 4) * (a1 ^ 2 * a5 ^ 2 * x ^ 11)
    + (1 / 2) * (a1 ^ 2 * a4 * a6 * x ^ 11)
    + (1 / 6) * (a0 * a4 ^ 3 * x ^ 11)
    + (1) * (a0 * a3 * a4 * a5 * x ^ 11)
/-- chunk 5 of the order-`≥ 9` part of `P + P²/2 + P³/6 + P⁴/24`, divided by `x⁹` -/
noncomputable def Wc5 (a0 a1 a2 a3 a4 a5 a6 x : ℝ) : ℝ :=
  (1 / 2) * (a0 * a3 ^ 2 * a6 * x ^ 11)
    + (1 / 2) * (a0 * a2 * a5 ^ 2 * x ^ 11)
    + (1) * (a0 * a2 * a4 * a6 * x ^ 11)
    + (1) * (a0 * a1 * a5 * a6 * x ^ 11)
    + (1 / 4) * (a0 ^ 2 * a6 ^ 2 * x ^ 11)
    + (1 / 6) * (a5 ^ 3 * x ^ 12)
    + (1) * (a4 * a5 * a6 * x ^ 12)
    + (1 / 2) * (a3 * a6 ^ 2 * x ^ 12)
    + (1 / 6) * (a3 ^ 3 * a4 * x ^ 12)
    + (1 / 2) * (a2 * a3 * a4 ^ 2 * x ^ 12)
    + (1 / 2) * (a2 * a3 ^ 2 * a5 * x ^ 12)
    + (1 / 2) * (a2 ^ 2 * a4 * a5 * x ^ 12)
    + (1 / 2) * (a2 ^ 2 * a3 * a6 * x ^ 12)
    + (1 / 6) * (a1 * a4 ^ 3 * x ^ 12)
    + (1) * (a1 * a3 * a4 * a5 * x ^ 12)
    + (1 / 2) * (a1 * a3 ^ 2 * a6 * x ^ 12)
    + (1 / 2) * (a1 * a2 * a5 ^ 2 * x ^ 12)
    + (1) * (a1 * a2 * a4 * a6 * x ^ 12)
    + (1 / 2) * (a1 ^ 2 * a5 * a6 * x ^ 12)
    + (1 / 2) * (a0 * a4 ^ 2 * a5 * x ^ 12)
    + (1 / 2) * (a0 * a3 * a5 ^ 2 * x ^ 12)
    + (1) * (a0 * a3 * a4 * a6 * x ^ 12)
    + (1) * (a0 * a2 * a5 * a6 * x ^ 12)
    + (1 / 2) * (a0 * a1 * a6 ^ 2 * x ^ 12)
    + (1 / 2) * (a5 ^ 2 * a6 * x ^ 13)
    + (1 / 2) * (a4 * a6 ^ 2 * x ^ 13)
    + (1 / 4) * (a3 ^ 2 * a4 ^ 2 * x ^ 13)
    + (1 / 6) * (a3 ^ 3 * a5 * x ^ 13)
    + (1 / 6) * (a2 * a4 ^ 3 * x ^ 13)
    + (1) * (a2 * a3 * a4 * a5 * x ^ 13)
    + (1 / 2) * (a2 * a3 ^ 2 * a6 * x ^ 13)
    + (1 / 4) * (a2 ^ 2 * a5 ^ 2 * x ^ 13)
    + (1 / 2) * (a2 ^ 2 * a4 * a6 * x ^ 13)
    + (1 / 2) * (a1 * a4 ^ 2 * a5 * x ^ 13)
    + (1 / 2) * (a1 * a3 * a5 ^ 2 * x ^ 13)
    + (1) * (a1 * a3 * a4 * a6 * x ^ 13)
    + (1) * (a1 * a2 * a5 * a6 * x ^ 13)
    + (1 / 4) * (a1 ^ 2 * a6 ^ 2 * x ^ 13)
    + (1 / 2) * (a0 * a4 * a5 ^ 2 * x ^ 13)
    + (1 / 2) * (a0 * a4 ^ 2 * a6 * x ^ 13)
/-- chunk 6 of the order-`≥ 9` part of `P + P²/2 + P³/6 + P⁴/24`, divided by `x⁹` -/
noncomputable def Wc6 (a0 a1 a2 a3 a4 a5 a6 x : ℝ) : ℝ :=
  (1) * (a0 * a3 * a5 * a6 * x ^ 13)
    + (1 / 2) * (a0 * a2 * a6 ^ 2 * x ^ 13)
    + (1 / 2) * (a5 * a6 ^ 2 * x ^ 14)
    + (1 / 6) * (a3 * a4 ^ 3 * x ^ 14)
    + (1 / 2) * (a3 ^ 2 * a4 * a5 * x ^ 14)
    + (1 / 6) * (a3 ^ 3 * a6 * x ^ 14)
    + (1 / 2) * (a2 * a4 ^ 2 * a5 * x ^ 14)
    + (1 / 2) * (a2 * a3 * a5 ^ 2 * x ^ 14)
    + (1) * (a2 * a3 * a4 * a6 * x ^ 14)
    + (1 / 2) * (a2 ^ 2 * a5 * a6 * x ^ 14)
    + (1 / 2) * (a1 * a4 * a5 ^ 2 * x ^ 14)
    + (1 / 2) * (a1 * a4 ^ 2 * a6 * x ^ 14)
    + (1) * (a1 * a3 * a5 * a6 * x ^ 14)
    + (1 / 2) * (a1 * a2 * a6 ^ 2 * x ^ 14)
    + (1 / 6) * (a0 * a5 ^ 3 * x ^ 14)
    + (1) * (a0 * a4 * a5 * a6 * x ^ 14)
    + (1 / 2) * (a0 * a3 * a6 ^ 2 * x ^ 14)
    + (1 / 6) * (a6 ^ 3 * x ^ 15)
    + (1 / 24) * (a4 ^ 4 * x ^ 15)
    + (1 / 2) * (a3 * a4 ^ 2 * a5 * x ^ 15)
    + (1 / 4) * (a3 ^ 2 * a5 ^ 2 * x ^ 15)
    + (1 / 2) * (a3 ^ 2 * a4 * a6 * x ^ 15)
    + (1 / 2) * (a2 * a4 * a5 ^ 2 * x ^ 15)
    + (1 / 2) * (a2 * a4 ^ 2 * a6 * x ^ 15)
    + (1) * (a2 * a3 * a5 * a6 * x ^ 15)
    + (1 / 4) * (a2 ^ 2 * a6 ^ 2 * x ^ 15)
    + (1 / 6) * (a1 * a5 ^ 3 * x ^ 15)
    + (1) * (a1 * a4 * a5 * a6 * x ^ 15)
    + (1 / 2) * (a1 * a3 * a6 ^ 2 * x ^ 15)
    + (1 / 2) * (a0 * a5 ^ 2 * a6 * x ^ 15)
    + (1 / 2) * (a0 * a4 * a6 ^ 2 * x ^ 15)
    + (1 / 6) * (a4 ^ 3 * a5 * x ^ 16)
    + (1 / 2) * (a3 * a4 * a5 ^ 2 * x ^ 16)
    + (1 / 2) * (a3 * a4 ^ 2 * a6 * x ^ 16)
    + (1 / 2) * (a3 ^ 2 * a5 * a6 * x ^ 16)
    + (1 / 6) * (a2 * a5 ^ 3 * x ^ 16)
    + (1) * (a2 * a4 * a5 * a6 * x ^ 16)
    + (1 / 2) * (a2 * a3 * a6 ^ 2 * x ^ 16)
    + (1 / 2) * (a1 * a5 ^ 2 * a6 * x ^ 16)
    + (1 / 2) * (a1 * a4 * a6 ^ 2 * x ^ 16)
/-- chunk 7 of the order-`≥ 9` part of `P + P²/2 + P³/6 + P⁴/24`, divided by `x⁹` -/
noncomputable def Wc7 (a0 a1 a2 a3 a4 a5 a6 x : ℝ) : ℝ :=
  (1 / 2) * (a0 * a5 * a6 ^ 2 * x ^ 16)
    + (1 / 4) * (a4 ^ 2 * a5 ^ 2 * x ^ 17)
    + (1 / 6) * (a4 ^ 3 * a6 * x ^ 17)
    + (1 / 6) * (a3 * a5 ^ 3 * x ^ 17)
    + (1) * (a3 * a4 * a5 * a6 * x ^ 17)
    + (1 / 4) * (a3 ^ 2 * a6 ^ 2 * x ^ 17)
    + (1 / 2) * (a2 * a5 ^ 2 * a6 * x ^ 17)
    + (1 / 2) * (a2 * a4 * a6 ^ 2 * x ^ 17)
    + (1 / 2) * (a1 * a5 * a6 ^ 2 * x ^ 17)
    + (1 / 6) * (a0 * a6 ^ 3 * x ^ 17)
    + (1 / 6) * (a4 * a5 ^ 3 * x ^ 18)
    + (1 / 2) * (a4 ^ 2 * a5 * a6 * x ^ 18)
    + (1 / 2) * (a3 * a5 ^ 2 * a6 * x ^ 18)
    + (1 / 2) * (a3 * a4 * a6 ^ 2 * x ^ 18)
    + (1 / 2) * (a2 * a5 * a6 ^ 2 * x ^ 18)
    + (1 / 6) * (a1 * a6 ^ 3 * x ^ 18)
    + (1 / 24) * (a5 ^ 4 * x ^ 19)
    + (1 / 2) * (a4 * a5 ^ 2 * a6 * x ^ 19)
    + (1 / 4) * (a4 ^ 2 * a6 ^ 2 * x ^ 19)
    + (1 / 2) * (a3 * a5 * a6 ^ 2 * x ^ 19)
    + (1 / 6) * (a2 * a6 ^ 3 * x ^ 19)
    + (1 / 6) * (a5 ^ 3 * a6 * x ^ 20)
    + (1 / 2) * (a4 * a5 * a6 ^ 2 * x ^ 20)
    + (1 / 6) * (a3 * a6 ^ 3 * x ^ 20)
    + (1 / 4) * (a5 ^ 2 * a6 ^ 2 * x ^ 21)
    + (1 / 6) * (a4 * a6 ^ 3 * x ^ 21)
    + (1 / 6) * (a5 * a6 ^ 3 * x ^ 22)
    + (1 / 24) * (a6 ^ 4 * x ^ 23)
/-- the part of `P + P²/2 + P³/6 + P⁴/24` of order `≥ 9`, divided by `x⁹` -/
noncomputable def Wpoly (a0 a1 a2 a3 a4 a5 a6 x : ℝ) : ℝ :=
  Wc0 a0 a1 a2 a3 a4 a5 a6 x + Wc1 a0 a1 a2 a3 a4 a5 a6 x + Wc2 a0 a1 a2 a3 a4 a5 a6 x + Wc3 a0 a1 a2 a3 a4 a5 a6 x + Wc4 a0 a1 a2 a3 a4 a5 a6 x + Wc5 a0 a1 a2 a3 a4 a5 a6 x + Wc6 a0 a1 a2 a3 a4 a5 a6 x + Wc7 a0 a1 a2 a3 a4 a5 a6 x

set_option maxHeartbeats 16000000 in
lemma exp_trunc8 {a0 a1 a2 a3 a4 a5 a6 B K : ℝ} (hB : 0 ≤ B) (_hK : 0 ≤ K)
    (ha0 : |a0| ≤ B) (ha1 : |a1| ≤ B) (ha2 : |a2| ≤ B) (ha3 : |a3| ≤ B) (ha4 : |a4| ≤ B) (ha5 : |a5| ≤ B)
    (ha6 : |a6| ≤ B) :
    ∃ KW : ℝ, 0 ≤ KW ∧ ∀ x L : ℝ, 0 ≤ x → x ≤ 1 → 7 * B * x ^ 2 ≤ 1 → K * x ^ 9 ≤ 1 →
      |L - (a0 * x ^ 2 + a1 * x ^ 3 + a2 * x ^ 4 + a3 * x ^ 5 + a4 * x ^ 6 + a5 * x ^ 7 + a6 * x ^ 8)| ≤ K * x ^ 9 →
      |Real.exp L - (1 + a0 * x ^ 2 + a1 * x ^ 3 + (a2 + a0 ^ 2 / 2) * x ^ 4 + (a3 + a0 * a1) * x ^ 5 + (a4 + a0 * a2 + a1 ^ 2 / 2 + a0 ^ 3 / 6) * x ^ 6 + (a5 + a0 * a3 + a1 * a2 + a0 ^ 2 * a1 / 2) * x ^ 7 + (a6 + a0 * a4 + a1 * a3 + a2 ^ 2 / 2 + a0 ^ 2 * a2 / 2 + a0 * a1 ^ 2 / 2 + a0 ^ 4 / 24) * x ^ 8)| ≤ (6 * K + (7 * B) ^ 5 + KW) * x ^ 9 := by
  have hWc : Continuous (fun x : ℝ => Wpoly a0 a1 a2 a3 a4 a5 a6 x) := by unfold Wpoly Wc0 Wc1 Wc2 Wc3 Wc4 Wc5 Wc6 Wc7; fun_prop
  obtain ⟨CW, hCW⟩ := (isCompact_Icc (a := (0 : ℝ)) (b := 1)).exists_bound_of_continuousOn hWc.continuousOn
  refine ⟨max CW 0, le_max_right _ _, fun x L hx0 hx1 hs1 hs2 hL => ?_⟩
  set W := Wpoly a0 a1 a2 a3 a4 a5 a6 x with hW
  have hWb : |W| ≤ max CW 0 := by
    have := hCW x ⟨hx0, hx1⟩
    rw [Real.norm_eq_abs] at this
    exact this.trans (le_max_left _ _)
  set P := a0 * x ^ 2 + a1 * x ^ 3 + a2 * x ^ 4 + a3 * x ^ 5 + a4 * x ^ 6 + a5 * x ^ 7 + a6 * x ^ 8 with hPdef
  have hterm : ∀ (a : ℝ) (k : ℕ), |a| ≤ B → |a * x ^ (k + 2)| ≤ B * x ^ 2 := by
    intro a k ha
    have : x ^ (k + 2) ≤ x ^ 2 := pow_le_pow_of_le_one hx0 hx1 (by omega)
    rw [abs_mul, abs_of_nonneg (pow_nonneg hx0 (k + 2))]
    exact mul_le_mul ha this (pow_nonneg hx0 _) hB
  have bP : |P| ≤ 7 * B * x ^ 2 := by
    have e0 : |a0 * x ^ 2| ≤ B * x ^ 2 := hterm a0 0 ha0
    have e1 : |a1 * x ^ 3| ≤ B * x ^ 2 := hterm a1 1 ha1
    have e2 : |a2 * x ^ 4| ≤ B * x ^ 2 := hterm a2 2 ha2
    have e3 : |a3 * x ^ 5| ≤ B * x ^ 2 := hterm a3 3 ha3
    have e4 : |a4 * x ^ 6| ≤ B * x ^ 2 := hterm a4 4 ha4
    have e5 : |a5 * x ^ 7| ≤ B * x ^ 2 := hterm a5 5 ha5
    have e6 : |a6 * x ^ 8| ≤ B * x ^ 2 := hterm a6 6 ha6
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
  have hP5 : |P| ^ 5 ≤ (7 * B) ^ 5 * x ^ 9 := by
    calc |P| ^ 5 ≤ (7 * B * x ^ 2) ^ 5 := pow_le_pow_left₀ (abs_nonneg _) bP 5
      _ = (7 * B) ^ 5 * x ^ 10 := by ring
      _ ≤ (7 * B) ^ 5 * x ^ 9 := by
          apply mul_le_mul_of_nonneg_left (pow_le_pow_of_le_one hx0 hx1 (by norm_num)) (by positivity)
  have key : Real.exp L - (1 + a0 * x ^ 2 + a1 * x ^ 3 + (a2 + a0 ^ 2 / 2) * x ^ 4 + (a3 + a0 * a1) * x ^ 5 + (a4 + a0 * a2 + a1 ^ 2 / 2 + a0 ^ 3 / 6) * x ^ 6 + (a5 + a0 * a3 + a1 * a2 + a0 ^ 2 * a1 / 2) * x ^ 7 + (a6 + a0 * a4 + a1 * a3 + a2 ^ 2 / 2 + a0 ^ 2 * a2 / 2 + a0 * a1 ^ 2 / 2 + a0 ^ 4 / 24) * x ^ 8)
      = Real.exp P * (Real.exp E - 1) + (Real.exp P - (1 + P + P ^ 2 / 2 + P ^ 3 / 6 + P ^ 4 / 24)) + x ^ 9 * W := by
    rw [hLPE, Real.exp_add, hW, hPdef]; unfold Wpoly Wc0 Wc1 Wc2 Wc3 Wc4 Wc5 Wc6 Wc7; ring
  rw [key]
  have b1 : |Real.exp P * (Real.exp E - 1)| ≤ 3 * (2 * (K * x ^ 9)) := by
    rw [abs_mul, abs_of_pos (Real.exp_pos _)]
    exact mul_le_mul hexpP (hexpE.trans (by linarith)) (abs_nonneg _) (by norm_num)
  have b2 : |Real.exp P - (1 + P + P ^ 2 / 2 + P ^ 3 / 6 + P ^ 4 / 24)| ≤ (7 * B) ^ 5 * x ^ 9 := by
    calc _ ≤ _ := heP
      _ ≤ |P| ^ 5 * 1 := by gcongr; norm_num [Nat.factorial]
      _ ≤ _ := by rw [mul_one]; exact hP5
  have b3 : |x ^ 9 * W| ≤ x ^ 9 * max CW 0 := by
    rw [abs_mul, abs_of_nonneg (by positivity)]
    exact mul_le_mul_of_nonneg_left hWb (by positivity)
  have t1 := abs_add_le (Real.exp P * (Real.exp E - 1) + (Real.exp P - (1 + P + P ^ 2 / 2 + P ^ 3 / 6 + P ^ 4 / 24)))
    (x ^ 9 * W)
  have t2 := abs_add_le (Real.exp P * (Real.exp E - 1)) (Real.exp P - (1 + P + P ^ 2 / 2 + P ^ 3 / 6 + P ^ 4 / 24))
  have : (6 * K + (7 * B) ^ 5 + max CW 0) * x ^ 9 = 3 * (2 * (K * x ^ 9)) + (7 * B) ^ 5 * x ^ 9
    + x ^ 9 * max CW 0 := by ring
  rw [this]
  linarith only [t1, t2, b1, b2, b3]


lemma A_expand8 {α β γ : ℝ} (hα : |α| ≤ 3 / 2) (hβ : |β| ≤ 3 / 2) (hγ : |γ| ≤ 3 / 2) :
    ∃ K8 : ℝ, 0 ≤ K8 ∧ ∀ x L : ℝ, 0 ≤ x → x ≤ 1 / 40 →
    HasSum (fun k : ℕ => (-1 : ℝ) ^ (k + 2) * Y04.zn (k + 2) / ((k + 2 : ℕ) : ℝ) *
      ((α * x) ^ (k + 2) + (β * x) ^ (k + 2) - (γ * x) ^ (k + 2))) L →
    |Real.exp L - (1 + Y04.cf α β γ 0 * x ^ 2 + Y04.cf α β γ 1 * x ^ 3 + (Y04.cf α β γ 2 + Y04.cf α β γ 0 ^ 2 / 2) * x ^ 4 + (Y04.cf α β γ 3 + Y04.cf α β γ 0 * Y04.cf α β γ 1) * x ^ 5 + (Y04.cf α β γ 4 + Y04.cf α β γ 0 * Y04.cf α β γ 2 + Y04.cf α β γ 1 ^ 2 / 2 + Y04.cf α β γ 0 ^ 3 / 6) * x ^ 6 + (Y04.cf α β γ 5 + Y04.cf α β γ 0 * Y04.cf α β γ 3 + Y04.cf α β γ 1 * Y04.cf α β γ 2 + Y04.cf α β γ 0 ^ 2 * Y04.cf α β γ 1 / 2) * x ^ 7 + (Y04.cf α β γ 6 + Y04.cf α β γ 0 * Y04.cf α β γ 4 + Y04.cf α β γ 1 * Y04.cf α β γ 3 + Y04.cf α β γ 2 ^ 2 / 2 + Y04.cf α β γ 0 ^ 2 * Y04.cf α β γ 2 / 2 + Y04.cf α β γ 0 * Y04.cf α β γ 1 ^ 2 / 2 + Y04.cf α β γ 0 ^ 4 / 24) * x ^ 8)| ≤ K8 * x ^ 9 := by
  have hB : ∀ k ≤ 6, |Y04.cf α β γ k| ≤ 77 := by
    intro k hk
    have := Y04.cf_bd hα hβ hγ k
    have : (3 / 2 : ℝ) ^ (k + 2) ≤ (3 / 2) ^ 8 := pow_le_pow_right₀ (by norm_num) (by omega)
    norm_num at this ⊢; linarith
  obtain ⟨KW, hKW, hE⟩ := exp_trunc8 (B := 77) (K := 231) (by norm_num) (by norm_num) (hB 0 (by norm_num))
    (hB 1 (by norm_num)) (hB 2 (by norm_num)) (hB 3 (by norm_num)) (hB 4 (by norm_num)) (hB 5 (by norm_num))
    (hB 6 (by norm_num))
  refine ⟨6 * 231 + (7 * 77) ^ 5 + KW, by positivity, fun x L hx0 hx hL => ?_⟩
  have h := series_trunc8 hα hβ hγ hx0 (by linarith) hL
  have hx2 : x ^ 2 ≤ 1 / 1600 := by
    have := pow_le_pow_left₀ hx0 hx 2
    norm_num at this ⊢; linarith
  have hx9 : x ^ 9 ≤ 1 / 1600 := by
    have : x ^ 9 ≤ x ^ 2 := pow_le_pow_of_le_one hx0 (by linarith) (by norm_num)
    linarith
  exact hE x L hx0 (by linarith) (by nlinarith) (by nlinarith) h

set_option maxHeartbeats 16000000 in
/-- the `x⁷` coefficient of the regularised Beta identity at `a = αx`, `b = βx` -/
theorem coef7 {α β : ℝ} (hα : 0 < α) (hβ : 0 < β) (h3 : α + β ≤ 3) :
    (-1 / 720) * ((α + β) ^ 1 * (α ^ 6 + β ^ 6)) * Y04.I 6 1 + (-1 / 240) * ((α + β) ^ 2 * (α ^ 5 + β ^ 5)) * Y04.I 5 2 + (-1 / 144) * ((α + β) ^ 3 * (α ^ 4 + β ^ 4)) * Y04.I 4 3 + (-1 / 144) * ((α + β) ^ 4 * (α ^ 3 + β ^ 3)) * Y04.I 3 4 + (-1 / 240) * ((α + β) ^ 5 * (α ^ 2 + β ^ 2)) * Y04.I 2 5 + (-1 / 720) * ((α + β) ^ 6 * (α ^ 1 + β ^ 1)) * Y04.I 1 6 + (-1 / 5040) * ((α + β) ^ 7 * (α ^ 0 + β ^ 0)) * Y04.I 0 7 = (α + β) / (α * β) * 256 * (Y04.cf (α / 2) (β / 2) ((α + β) / 2) 6 + Y04.cf (α / 2) (β / 2) ((α + β) / 2) 0 * Y04.cf (α / 2) (β / 2) ((α + β) / 2) 4 + Y04.cf (α / 2) (β / 2) ((α + β) / 2) 1 * Y04.cf (α / 2) (β / 2) ((α + β) / 2) 3 + Y04.cf (α / 2) (β / 2) ((α + β) / 2) 2 ^ 2 / 2 + Y04.cf (α / 2) (β / 2) ((α + β) / 2) 0 ^ 2 * Y04.cf (α / 2) (β / 2) ((α + β) / 2) 2 / 2 + Y04.cf (α / 2) (β / 2) ((α + β) / 2) 0 * Y04.cf (α / 2) (β / 2) ((α + β) / 2) 1 ^ 2 / 2 + Y04.cf (α / 2) (β / 2) ((α + β) / 2) 0 ^ 4 / 24) := by
  set γ := α + β with hγ
  set c0 := Y04.cf (α / 2) (β / 2) (γ / 2) 0 with hc0
  set c1 := Y04.cf (α / 2) (β / 2) (γ / 2) 1 with hc1
  set c2 := Y04.cf (α / 2) (β / 2) (γ / 2) 2 with hc2
  set c3 := Y04.cf (α / 2) (β / 2) (γ / 2) 3 with hc3
  set c4 := Y04.cf (α / 2) (β / 2) (γ / 2) 4 with hc4
  set c5 := Y04.cf (α / 2) (β / 2) (γ / 2) 5 with hc5
  set c6 := Y04.cf (α / 2) (β / 2) (γ / 2) 6 with hc6
  set KB := ∫ s in Set.Ioi (0 : ℝ), B9 s with hKB
  have hKB0 : 0 ≤ KB := setIntegral_nonneg measurableSet_Ioi (fun s hs => B9_nonneg (le_of_lt hs))
  have hα3 : α ≤ 3 := by linarith
  have hβ3 : β ≤ 3 := by linarith
  have hγ0 : 0 < γ := by rw [hγ]; linarith
  obtain ⟨K8, hK8, hA8⟩ := A_expand8 (α := α / 2) (β := β / 2) (γ := γ / 2)
      (by rw [abs_of_pos (by positivity)]; linarith) (by rw [abs_of_pos (by positivity)]; linarith)
      (by rw [abs_of_pos (by positivity)]; linarith)
  set e := γ / (α * β) with he
  have he0 : 0 ≤ e := by positivity
  set d1 := (-1) * (γ ^ 1 * (α ^ 0 + β ^ 0)) * Y04.I 0 1 - e * 4 * (c0) with hd1
  set d2 := (1) * (γ ^ 1 * (α ^ 1 + β ^ 1)) * Y04.I 1 1 + (1 / 2) * (γ ^ 2 * (α ^ 0 + β ^ 0)) * Y04.I 0 2 - e * 8 * (c1) with hd2
  set d3 := (-1 / 2) * (γ ^ 1 * (α ^ 2 + β ^ 2)) * Y04.I 2 1 + (-1 / 2) * (γ ^ 2 * (α ^ 1 + β ^ 1)) * Y04.I 1 2 + (-1 / 6) * (γ ^ 3 * (α ^ 0 + β ^ 0)) * Y04.I 0 3 - e * 16 * (c2 + c0 ^ 2 / 2) with hd3
  set d4 := (1 / 6) * (γ ^ 1 * (α ^ 3 + β ^ 3)) * Y04.I 3 1 + (1 / 4) * (γ ^ 2 * (α ^ 2 + β ^ 2)) * Y04.I 2 2 + (1 / 6) * (γ ^ 3 * (α ^ 1 + β ^ 1)) * Y04.I 1 3 + (1 / 24) * (γ ^ 4 * (α ^ 0 + β ^ 0)) * Y04.I 0 4 - e * 32 * (c3 + c0 * c1) with hd4
  set d5 := (-1 / 24) * (γ ^ 1 * (α ^ 4 + β ^ 4)) * Y04.I 4 1 + (-1 / 12) * (γ ^ 2 * (α ^ 3 + β ^ 3)) * Y04.I 3 2 + (-1 / 12) * (γ ^ 3 * (α ^ 2 + β ^ 2)) * Y04.I 2 3 + (-1 / 24) * (γ ^ 4 * (α ^ 1 + β ^ 1)) * Y04.I 1 4 + (-1 / 120) * (γ ^ 5 * (α ^ 0 + β ^ 0)) * Y04.I 0 5 - e * 64 * (c4 + c0 * c2 + c1 ^ 2 / 2 + c0 ^ 3 / 6) with hd5
  set d6 := (1 / 120) * (γ ^ 1 * (α ^ 5 + β ^ 5)) * Y04.I 5 1 + (1 / 48) * (γ ^ 2 * (α ^ 4 + β ^ 4)) * Y04.I 4 2 + (1 / 36) * (γ ^ 3 * (α ^ 3 + β ^ 3)) * Y04.I 3 3 + (1 / 48) * (γ ^ 4 * (α ^ 2 + β ^ 2)) * Y04.I 2 4 + (1 / 120) * (γ ^ 5 * (α ^ 1 + β ^ 1)) * Y04.I 1 5 + (1 / 720) * (γ ^ 6 * (α ^ 0 + β ^ 0)) * Y04.I 0 6 - e * 128 * (c5 + c0 * c3 + c1 * c2 + c0 ^ 2 * c1 / 2) with hd6
  set d7 := (-1 / 720) * (γ ^ 1 * (α ^ 6 + β ^ 6)) * Y04.I 6 1 + (-1 / 240) * (γ ^ 2 * (α ^ 5 + β ^ 5)) * Y04.I 5 2 + (-1 / 144) * (γ ^ 3 * (α ^ 4 + β ^ 4)) * Y04.I 4 3 + (-1 / 144) * (γ ^ 4 * (α ^ 3 + β ^ 3)) * Y04.I 3 4 + (-1 / 240) * (γ ^ 5 * (α ^ 2 + β ^ 2)) * Y04.I 2 5 + (-1 / 720) * (γ ^ 6 * (α ^ 1 + β ^ 1)) * Y04.I 1 6 + (-1 / 5040) * (γ ^ 7 * (α ^ 0 + β ^ 0)) * Y04.I 0 7 - e * 256 * (c6 + c0 * c4 + c1 * c3 + c2 ^ 2 / 2 + c0 ^ 2 * c2 / 2 + c0 * c1 ^ 2 / 2 + c0 ^ 4 / 24) with hd7
  have hmain : ∀ x : ℝ, 0 < x → x ≤ 1 / 80 →
      |d1 * x + d2 * x ^ 2 + d3 * x ^ 3 + d4 * x ^ 4 + d5 * x ^ 5 + d6 * x ^ 6 + d7 * x ^ 7|
        ≤ (60 * 6561 * KB + e * 512 * K8) * x ^ 8 := by
    intro x hx0 hx1
    have hx3 : 3 * x ≤ 1 := by linarith
    have gx : γ * x ≤ 3 * x := mul_le_mul_of_nonneg_right (by rw [hγ]; exact h3) hx0.le
    have ax : α * x ≤ 3 * x := mul_le_mul_of_nonneg_right hα3 hx0.le
    have bx : β * x ≤ 3 * x := mul_le_mul_of_nonneg_right hβ3 hx0.le
    have dx1 : (α - β) * x ≤ 3 * x := mul_le_mul_of_nonneg_right (by linarith) hx0.le
    have dx2 : -3 * x ≤ (α - β) * x := mul_le_mul_of_nonneg_right (by linarith) hx0.le
    obtain ⟨_, E1⟩ := int_expand7 (u := γ * x) (v := α * x) (t := 3 * x) (by positivity)
      gx (by positivity) ax hx3
    obtain ⟨_, E2⟩ := int_expand7 (u := γ * x) (v := β * x) (t := 3 * x) (by positivity)
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
    have hA := hA8 (2 * x) _ (by positivity) (by linarith) h
    rw [Real.exp_log hpos] at hA
    rw [show α / 2 * (2 * x) = α * x by ring, show β / 2 * (2 * x) = β * x by ring,
      show γ / 2 * (2 * x) = γ * x by ring] at hA
    rw [← hc0, ← hc1, ← hc2, ← hc3, ← hc4, ← hc5, ← hc6] at hA
    set Rg := Real.Gamma (1 + α * x) * Real.Gamma (1 + β * x) / Real.Gamma (1 + γ * x) with hRg
    set Iα := ∫ s in Set.Ioi (0 : ℝ), (Real.exp (-(γ * x) * Y04.ell s) - 1) * Real.exp (-(α * x) * s)
    set Iβ := ∫ s in Set.Ioi (0 : ℝ), (Real.exp (-(γ * x) * Y04.ell s) - 1) * Real.exp (-(β * x) * s)
    set u := γ * x with hu
    set Sα := (-1) * (u ^ 1 * (α * x) ^ 0) * Y04.I 0 1 + (1) * (u ^ 1 * (α * x) ^ 1) * Y04.I 1 1 + (-1 / 2) * (u ^ 1 * (α * x) ^ 2) * Y04.I 2 1 + (1 / 6) * (u ^ 1 * (α * x) ^ 3) * Y04.I 3 1 + (-1 / 24) * (u ^ 1 * (α * x) ^ 4) * Y04.I 4 1 + (1 / 120) * (u ^ 1 * (α * x) ^ 5) * Y04.I 5 1 + (-1 / 720) * (u ^ 1 * (α * x) ^ 6) * Y04.I 6 1 + (1 / 2) * (u ^ 2 * (α * x) ^ 0) * Y04.I 0 2 + (-1 / 2) * (u ^ 2 * (α * x) ^ 1) * Y04.I 1 2 + (1 / 4) * (u ^ 2 * (α * x) ^ 2) * Y04.I 2 2 + (-1 / 12) * (u ^ 2 * (α * x) ^ 3) * Y04.I 3 2 + (1 / 48) * (u ^ 2 * (α * x) ^ 4) * Y04.I 4 2 + (-1 / 240) * (u ^ 2 * (α * x) ^ 5) * Y04.I 5 2 + (-1 / 6) * (u ^ 3 * (α * x) ^ 0) * Y04.I 0 3 + (1 / 6) * (u ^ 3 * (α * x) ^ 1) * Y04.I 1 3 + (-1 / 12) * (u ^ 3 * (α * x) ^ 2) * Y04.I 2 3 + (1 / 36) * (u ^ 3 * (α * x) ^ 3) * Y04.I 3 3 + (-1 / 144) * (u ^ 3 * (α * x) ^ 4) * Y04.I 4 3 + (1 / 24) * (u ^ 4 * (α * x) ^ 0) * Y04.I 0 4 + (-1 / 24) * (u ^ 4 * (α * x) ^ 1) * Y04.I 1 4 + (1 / 48) * (u ^ 4 * (α * x) ^ 2) * Y04.I 2 4 + (-1 / 144) * (u ^ 4 * (α * x) ^ 3) * Y04.I 3 4 + (-1 / 120) * (u ^ 5 * (α * x) ^ 0) * Y04.I 0 5 + (1 / 120) * (u ^ 5 * (α * x) ^ 1) * Y04.I 1 5 + (-1 / 240) * (u ^ 5 * (α * x) ^ 2) * Y04.I 2 5 + (1 / 720) * (u ^ 6 * (α * x) ^ 0) * Y04.I 0 6 + (-1 / 720) * (u ^ 6 * (α * x) ^ 1) * Y04.I 1 6 + (-1 / 5040) * (u ^ 7 * (α * x) ^ 0) * Y04.I 0 7 with hSα
    set Sβ := (-1) * (u ^ 1 * (β * x) ^ 0) * Y04.I 0 1 + (1) * (u ^ 1 * (β * x) ^ 1) * Y04.I 1 1 + (-1 / 2) * (u ^ 1 * (β * x) ^ 2) * Y04.I 2 1 + (1 / 6) * (u ^ 1 * (β * x) ^ 3) * Y04.I 3 1 + (-1 / 24) * (u ^ 1 * (β * x) ^ 4) * Y04.I 4 1 + (1 / 120) * (u ^ 1 * (β * x) ^ 5) * Y04.I 5 1 + (-1 / 720) * (u ^ 1 * (β * x) ^ 6) * Y04.I 6 1 + (1 / 2) * (u ^ 2 * (β * x) ^ 0) * Y04.I 0 2 + (-1 / 2) * (u ^ 2 * (β * x) ^ 1) * Y04.I 1 2 + (1 / 4) * (u ^ 2 * (β * x) ^ 2) * Y04.I 2 2 + (-1 / 12) * (u ^ 2 * (β * x) ^ 3) * Y04.I 3 2 + (1 / 48) * (u ^ 2 * (β * x) ^ 4) * Y04.I 4 2 + (-1 / 240) * (u ^ 2 * (β * x) ^ 5) * Y04.I 5 2 + (-1 / 6) * (u ^ 3 * (β * x) ^ 0) * Y04.I 0 3 + (1 / 6) * (u ^ 3 * (β * x) ^ 1) * Y04.I 1 3 + (-1 / 12) * (u ^ 3 * (β * x) ^ 2) * Y04.I 2 3 + (1 / 36) * (u ^ 3 * (β * x) ^ 3) * Y04.I 3 3 + (-1 / 144) * (u ^ 3 * (β * x) ^ 4) * Y04.I 4 3 + (1 / 24) * (u ^ 4 * (β * x) ^ 0) * Y04.I 0 4 + (-1 / 24) * (u ^ 4 * (β * x) ^ 1) * Y04.I 1 4 + (1 / 48) * (u ^ 4 * (β * x) ^ 2) * Y04.I 2 4 + (-1 / 144) * (u ^ 4 * (β * x) ^ 3) * Y04.I 3 4 + (-1 / 120) * (u ^ 5 * (β * x) ^ 0) * Y04.I 0 5 + (1 / 120) * (u ^ 5 * (β * x) ^ 1) * Y04.I 1 5 + (-1 / 240) * (u ^ 5 * (β * x) ^ 2) * Y04.I 2 5 + (1 / 720) * (u ^ 6 * (β * x) ^ 0) * Y04.I 0 6 + (-1 / 720) * (u ^ 6 * (β * x) ^ 1) * Y04.I 1 6 + (-1 / 5040) * (u ^ 7 * (β * x) ^ 0) * Y04.I 0 7 with hSβ
    set Q := 1 + c0 * (2 * x) ^ 2 + c1 * (2 * x) ^ 3 + (c2 + c0 ^ 2 / 2) * (2 * x) ^ 4 + (c3 + c0 * c1) * (2 * x) ^ 5 + (c4 + c0 * c2 + c1 ^ 2 / 2 + c0 ^ 3 / 6) * (2 * x) ^ 6 + (c5 + c0 * c3 + c1 * c2 + c0 ^ 2 * c1 / 2) * (2 * x) ^ 7 + (c6 + c0 * c4 + c1 * c3 + c2 ^ 2 / 2 + c0 ^ 2 * c2 / 2 + c0 * c1 ^ 2 / 2 + c0 ^ 4 / 24) * (2 * x) ^ 8 with hQ
    have hpoly : d1 * x + d2 * x ^ 2 + d3 * x ^ 3 + d4 * x ^ 4 + d5 * x ^ 5 + d6 * x ^ 6 + d7 * x ^ 7
        = (Sα - Iα) + (Sβ - Iβ) + ((Iα + Iβ) - e * (4 * (c0) * x ^ 1 + 8 * (c1) * x ^ 2 + 16 * (c2 + c0 ^ 2 / 2) * x ^ 3 + 32 * (c3 + c0 * c1) * x ^ 4 + 64 * (c4 + c0 * c2 + c1 ^ 2 / 2 + c0 ^ 3 / 6) * x ^ 5 + 128 * (c5 + c0 * c3 + c1 * c2 + c0 ^ 2 * c1 / 2) * x ^ 6 + 256 * (c6 + c0 * c4 + c1 * c3 + c2 ^ 2 / 2 + c0 ^ 2 * c2 / 2 + c0 * c1 ^ 2 / 2 + c0 ^ 4 / 24) * x ^ 7)) := by
      rw [hd1, hd2, hd3, hd4, hd5, hd6, hd7, hSα, hSβ, hu]; ring
    have hlast : (Iα + Iβ) - e * (4 * (c0) * x ^ 1 + 8 * (c1) * x ^ 2 + 16 * (c2 + c0 ^ 2 / 2) * x ^ 3 + 32 * (c3 + c0 * c1) * x ^ 4 + 64 * (c4 + c0 * c2 + c1 ^ 2 / 2 + c0 ^ 3 / 6) * x ^ 5 + 128 * (c5 + c0 * c3 + c1 * c2 + c0 ^ 2 * c1 / 2) * x ^ 6 + 256 * (c6 + c0 * c4 + c1 * c3 + c2 ^ 2 / 2 + c0 ^ 2 * c2 / 2 + c0 * c1 ^ 2 / 2 + c0 ^ 4 / 24) * x ^ 7) = e / x * (Rg - Q) := by
      rw [hC, he, hQ, hu]
      field_simp
      ring
    rw [hpoly, hlast]
    have b3 : |e / x * (Rg - Q)| ≤ e * 512 * K8 * x ^ 8 := by
      rw [abs_mul, abs_of_nonneg (by positivity)]
      calc e / x * |Rg - Q| ≤ e / x * (K8 * (2 * x) ^ 9) :=
            mul_le_mul_of_nonneg_left hA (by positivity)
        _ = e * 512 * K8 * x ^ 8 := by field_simp; ring
    have b1 : |Sα - Iα| ≤ 30 * (3 * x) ^ 8 * KB := by rw [abs_sub_comm]; exact E1
    have b2 : |Sβ - Iβ| ≤ 30 * (3 * x) ^ 8 * KB := by rw [abs_sub_comm]; exact E2
    have t1 := abs_add_le ((Sα - Iα) + (Sβ - Iβ)) (e / x * (Rg - Q))
    have t2 := abs_add_le (Sα - Iα) (Sβ - Iβ)
    have : (60 * 6561 * KB + e * 512 * K8) * x ^ 8
        = 30 * (3 * x) ^ 8 * KB + 30 * (3 * x) ^ 8 * KB + e * 512 * K8 * x ^ 8 := by ring
    linarith
  obtain ⟨_, _, _, _, _, _, h7⟩ := septic_zero (by norm_num) (by norm_num) (by positivity) hmain
  rw [hd7] at h7
  linarith


lemma D82 : 2 * (-(31 / 24) * Y04.I 6 1 + (1 / 2) * Y04.I 5 2 / 2 + Y04.I 4 3 / 6) =
    -(60889 / 32) * (π ^ 8 / 9450) + 24 * (Y04.z3 * Z02.z5)
      + 12 * ((π ^ 2 / 6) * Y04.z3 ^ 2) := by
  have h1 := coef7 (α := 1) (β := 1) one_pos one_pos (by norm_num)
  have h2 := coef7 (α := 1) (β := 2) one_pos two_pos (by norm_num)
  have h3 := coef7 (α := 1) (β := 1 / 3) one_pos (by norm_num) (by norm_num)
  have h4 := coef7 (α := 2) (β := 1 / 2) two_pos (by norm_num) (by norm_num)
  simp only [Y04.cf, Nat.reduceAdd, Nat.cast_ofNat, Y04.z_two, Y04.zn_three, Y04.z_four, Z02.zn_five,
    V02.z_six, z_eight] at h1 h2 h3 h4
  linear_combination (-153649 / 3) * h1 + (272976 / 35) * h2 + (-125826858 / 55) * h3
    + (12560384 / 1155) * h4

lemma D81 : 2 * (731 / 1152 * Y04.I 6 1 + (-(25 / 48)) * Y04.I 5 2 / 2 + (-(23 / 24)) * Y04.I 4 3 / 6
      + (1 / 2) * Y04.I 3 4 / 24 + Y04.I 2 5 / 120) =
    495965 / 512 * (π ^ 8 / 9450) - 25 * (Y04.z3 * Z02.z5)
      - 27 / 2 * ((π ^ 2 / 6) * Y04.z3 ^ 2) := by
  have h1 := coef7 (α := 1) (β := 1) one_pos one_pos (by norm_num)
  have h2 := coef7 (α := 1) (β := 2) one_pos two_pos (by norm_num)
  have h3 := coef7 (α := 1) (β := 1 / 3) one_pos (by norm_num) (by norm_num)
  have h4 := coef7 (α := 2) (β := 1 / 2) two_pos (by norm_num) (by norm_num)
  simp only [Y04.cf, Nat.reduceAdd, Nat.cast_ofNat, Y04.z_two, Y04.zn_three, Y04.z_four, Z02.zn_five,
    V02.z_six, z_eight] at h1 h2 h3 h4
  linear_combination (2005945 / 72) * h1 + (-445181 / 105) * h2 + (198822357 / 160) * h3
    + (-1859456 / 315) * h4

lemma D80 : 2 * (-(1369 / 9216) * Y04.I 6 1 + 721 / 3840 * Y04.I 5 2 / 2 + 121 / 384 * Y04.I 4 3 / 6
      + (-(17 / 48)) * Y04.I 3 4 / 24 + (-(5 / 8)) * Y04.I 2 5 / 120
      + (1 / 2) * Y04.I 1 6 / 720 + Y04.I 0 7 / 5040) =
    -(2891923 / 12288) * (π ^ 8 / 9450) + 141 / 16 * (Y04.z3 * Z02.z5)
      + 169 / 32 * ((π ^ 2 / 6) * Y04.z3 ^ 2) := by
  have h1 := coef7 (α := 1) (β := 1) one_pos one_pos (by norm_num)
  have h2 := coef7 (α := 1) (β := 2) one_pos two_pos (by norm_num)
  have h3 := coef7 (α := 1) (β := 1 / 3) one_pos (by norm_num) (by norm_num)
  have h4 := coef7 (α := 2) (β := 1 / 2) two_pos (by norm_num) (by norm_num)
  simp only [Y04.cf, Nat.reduceAdd, Nat.cast_ofNat, Y04.z_two, Y04.zn_three, Y04.z_four, Z02.zn_five,
    V02.z_six, z_eight] at h1 h2 h3 h4
  linear_combination (-920945 / 128) * h1 + (4377 / 4) * h2 + (-451219653 / 1408) * h3
    + (16736 / 11) * h4

lemma D83 : 2 * Y04.I 6 1 = 5715 / 4 * (π ^ 8 / 9450) := by rw [I61]; ring

lemma G_d6' (k : ℕ) : iteratedDeriv 6 (Y04.G k) 0 =
    720 * ((-1 : ℝ) ^ k * Real.sqrt 2 / 4 *
      (-(1369 / 9216) - 731 / 576 * (k : ℝ) ^ 2 - 31 / 36 * (k : ℝ) ^ 4 - 4 / 45 * (k : ℝ) ^ 6)) := by
  have := G_d6 k; linarith

lemma theta8 (k : ℕ) :
    (-1 : ℝ) ^ k * (π / Real.sqrt 2) *
        (8 / π * ∑ j ∈ Finset.Icc 1 7,
          iteratedDeriv (7 - j) (Y04.G k) 0 / ((7 - j).factorial : ℝ) / (j.factorial : ℝ) *
            Y04.I (7 - j) j) =
      -(2891923 / 12288 + 495965 / 256 * (k : ℝ) ^ 2 + 60889 / 48 * (k : ℝ) ^ 4
          + 127 * (k : ℝ) ^ 6) * (π ^ 8 / 9450)
        + (141 / 16 + 50 * (k : ℝ) ^ 2 + 16 * (k : ℝ) ^ 4) * (Y04.z3 * Z02.z5)
        + (169 / 32 + 27 * (k : ℝ) ^ 2 + 8 * (k : ℝ) ^ 4) * ((π ^ 2 / 6) * Y04.z3 ^ 2) := by
  set E := (-1 : ℝ) ^ k with hE
  set Y := (-(1369 / 9216) - 731 / 576 * (k : ℝ) ^ 2 - 31 / 36 * (k : ℝ) ^ 4 - 4 / 45 * (k : ℝ) ^ 6)
      * Y04.I 6 1
    + (721 / 3840 + 25 / 24 * (k : ℝ) ^ 2 + 1 / 3 * (k : ℝ) ^ 4) / 2 * Y04.I 5 2
    + (121 / 384 + 23 / 12 * (k : ℝ) ^ 2 + 2 / 3 * (k : ℝ) ^ 4) / 6 * Y04.I 4 3
    + (-(17 / 48 + (k : ℝ) ^ 2)) / 24 * Y04.I 3 4
    + (-(2 * (k : ℝ) ^ 2 + 5 / 8)) / 120 * Y04.I 2 5
    + (1 / 2) / 720 * Y04.I 1 6 + 1 / 5040 * Y04.I 0 7 with hY
  have hsum : ∑ j ∈ Finset.Icc 1 7,
          iteratedDeriv (7 - j) (Y04.G k) 0 / ((7 - j).factorial : ℝ) / (j.factorial : ℝ) *
            Y04.I (7 - j) j = E * Real.sqrt 2 / 4 * Y := by
    rw [show Finset.Icc 1 7 = {1, 2, 3, 4, 5, 6, 7} by decide]
    rw [Finset.sum_insert (by decide), Finset.sum_insert (by decide), Finset.sum_insert (by decide),
      Finset.sum_insert (by decide), Finset.sum_insert (by decide), Finset.sum_insert (by decide),
      Finset.sum_singleton]
    simp only [show 7 - 1 = 6 from rfl, show 7 - 2 = 5 from rfl, show 7 - 3 = 4 from rfl,
      show 7 - 4 = 3 from rfl, show 7 - 5 = 2 from rfl, show 7 - 6 = 1 from rfl,
      show 7 - 7 = 0 from rfl, iteratedDeriv_zero, iteratedDeriv_one]
    rw [G_d6', Q03.G_d5', V02.G_d4', Z02.G_d3', Z02.G_d2', Y04.G_d1, Y04.G_zero, hY]
    simp [Nat.factorial]
    ring
  rw [hsum]
  have hE2 : E ^ 2 = 1 := by rw [hE, ← pow_mul, mul_comm, pow_mul, neg_one_sq, one_pow]
  have hs : Real.sqrt 2 ≠ 0 := by positivity
  have hpi : π ≠ 0 := Real.pi_pos.ne'
  have e : E * (π / Real.sqrt 2) * (8 / π * (E * Real.sqrt 2 / 4 * Y)) = 2 * E ^ 2 * Y := by
    field_simp; ring
  rw [e, hE2, hY]
  linear_combination D80 - 2 * (k : ℝ) ^ 2 * D81 + 2 / 3 * (k : ℝ) ^ 4 * D82
    - 4 / 45 * (k : ℝ) ^ 6 * D83

lemma sum_eval8 (Gf : ℝ → ℝ) (τ : ℝ) :
    ∑ n ∈ Finset.Icc 2 8, τ ^ n * ∑ j ∈ Finset.Icc 1 (n - 1),
        iteratedDeriv (n - 1 - j) Gf 0 / ((n - 1 - j).factorial : ℝ) / (j.factorial : ℝ) * Y04.I (n - 1 - j) j
      = τ ^ 2 * Gf 0 * Y04.I 0 1 + τ ^ 3 * deriv Gf 0 * Y04.I 1 1 + τ ^ 3 * (Gf 0 / 2) * Y04.I 0 2
        + τ ^ 4 * (iteratedDeriv 2 Gf 0 / 2) * Y04.I 2 1
        + τ ^ 4 * (deriv Gf 0 / 2) * Y04.I 1 2 + τ ^ 4 * (Gf 0 / 6) * Y04.I 0 3
        + τ ^ 5 * (iteratedDeriv 3 Gf 0 / 6) * Y04.I 3 1 + τ ^ 5 * (iteratedDeriv 2 Gf 0 / 4) * Y04.I 2 2
        + τ ^ 5 * (deriv Gf 0 / 6) * Y04.I 1 3 + τ ^ 5 * (Gf 0 / 24) * Y04.I 0 4
        + τ ^ 6 * ∑ j ∈ Finset.Icc 1 5,
          iteratedDeriv (5 - j) Gf 0 / ((5 - j).factorial : ℝ) / (j.factorial : ℝ) * Y04.I (5 - j) j
        + τ ^ 7 * ∑ j ∈ Finset.Icc 1 6,
          iteratedDeriv (6 - j) Gf 0 / ((6 - j).factorial : ℝ) / (j.factorial : ℝ) * Y04.I (6 - j) j
        + τ ^ 8 * ∑ j ∈ Finset.Icc 1 7,
          iteratedDeriv (7 - j) Gf 0 / ((7 - j).factorial : ℝ) / (j.factorial : ℝ) * Y04.I (7 - j) j := by
  rw [Finset.sum_Icc_succ_top (by norm_num : 2 ≤ 7 + 1), Q03.sum_eval7]

set_option maxHeartbeats 8000000 in
lemma fourier8 (k : ℕ) : ∃ C τ₀ : ℝ, 0 < τ₀ ∧ ∀ p : ℝ, 1 < p → Y04.tau p ≤ τ₀ →
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
                  * ((π ^ 2 / 6) * Y04.z3 ^ 2)) * Y04.tau p ^ 8)|
      ≤ C * Y04.tau p ^ 9 := by
  obtain ⟨CW, hW⟩ := watson8 (Y04.G k) (Y04.G_contDiff k) ⟨1, fun x _ => Y04.G_bdd k x⟩
  refine ⟨8 / π * CW, 1 / 2, by norm_num, fun p hp _ => ?_⟩
  have hpi := Real.pi_pos
  have hτ0 := Y04.tau_pos hp
  have hτ1 := Y04.tau_lt_one hp
  obtain ⟨_, hb⟩ := hW (Y04.tau p) hτ0 hτ1
  rw [sum_eval8] at hb
  have h2 := Y04.layer_formula k hp
  have h1 := Y04.A_formula k
  set X := Y04.tau p * ∫ s in Set.Ioi (0 : ℝ), ((1 + Real.exp (-s)) ^ Y04.tau p - 1) * Y04.G k (Y04.tau p * s)
    with hX
  set τ := Y04.tau p with hτ
  set C6 := ∑ j ∈ Finset.Icc 1 5,
    iteratedDeriv (5 - j) (Y04.G k) 0 / ((5 - j).factorial : ℝ) / (j.factorial : ℝ) * Y04.I (5 - j) j
    with hC6
  set C7 := ∑ j ∈ Finset.Icc 1 6,
    iteratedDeriv (6 - j) (Y04.G k) 0 / ((6 - j).factorial : ℝ) / (j.factorial : ℝ) * Y04.I (6 - j) j
    with hC7
  set C8 := ∑ j ∈ Finset.Icc 1 7,
    iteratedDeriv (7 - j) (Y04.G k) 0 / ((7 - j).factorial : ℝ) / (j.factorial : ℝ) * Y04.I (7 - j) j
    with hC8
  set Θ := (4039 / 256 + 733 / 8 * (k : ℝ) ^ 2 + 31 * (k : ℝ) ^ 4) * (π ^ 6 / 945)
        - (7 / 8 + 2 * (k : ℝ) ^ 2) * Y04.z3 ^ 2 with hΘ
  set Θ7 := (129 / 16 + 46 * (k : ℝ) ^ 2 + 16 * (k : ℝ) ^ 4) * Q03.z7
        + (141 / 16 + 50 * (k : ℝ) ^ 2 + 16 * (k : ℝ) ^ 4) * ((π ^ 2 / 6) * Z02.z5)
        + (1071 / 160 + 181 / 5 * (k : ℝ) ^ 2 + 56 / 5 * (k : ℝ) ^ 4)
            * ((π ^ 2 / 6) ^ 2 * Y04.z3) with hΘ7
  set Θ8 := -(2891923 / 12288 + 495965 / 256 * (k : ℝ) ^ 2 + 60889 / 48 * (k : ℝ) ^ 4
          + 127 * (k : ℝ) ^ 6) * (π ^ 8 / 9450)
        + (141 / 16 + 50 * (k : ℝ) ^ 2 + 16 * (k : ℝ) ^ 4) * (Y04.z3 * Z02.z5)
        + (169 / 32 + 27 * (k : ℝ) ^ 2 + 8 * (k : ℝ) ^ 4) * ((π ^ 2 / 6) * Y04.z3 ^ 2) with hΘ8
  have hth := V02.theta6 k
  rw [← hC6, ← hΘ] at hth
  have hth7 := Q03.theta7 k
  rw [← hC7, ← hΘ7] at hth7
  have hth8 := theta8 k
  rw [← hC8, ← hΘ8] at hth8
  have hE2 : ((-1 : ℝ) ^ k) ^ 2 = 1 := by rw [← pow_mul, mul_comm, pow_mul, neg_one_sq, one_pow]
  have hsq : Real.sqrt 2 ^ 2 = 2 := Real.sq_sqrt (by norm_num)
  have hs0 : Real.sqrt 2 ≠ 0 := by positivity
  have hC6v : C6 = (-1 : ℝ) ^ k * Real.sqrt 2 / 8 * Θ := by
    rw [← hth]; field_simp
    linear_combination (-C6) * hE2
  have hC7v : C7 = (-1 : ℝ) ^ k * Real.sqrt 2 / 8 * Θ7 := by
    rw [← hth7]; field_simp
    linear_combination (-C7) * hE2
  have hC8v : C8 = (-1 : ℝ) ^ k * Real.sqrt 2 / 8 * Θ8 := by
    rw [← hth8]; field_simp
    linear_combination (-C8) * hE2
  set S := τ ^ 2 * Y04.G k 0 * Y04.I 0 1 + τ ^ 3 * deriv (Y04.G k) 0 * Y04.I 1 1
        + τ ^ 3 * (Y04.G k 0 / 2) * Y04.I 0 2
        + τ ^ 4 * (iteratedDeriv 2 (Y04.G k) 0 / 2) * Y04.I 2 1
        + τ ^ 4 * (deriv (Y04.G k) 0 / 2) * Y04.I 1 2 + τ ^ 4 * (Y04.G k 0 / 6) * Y04.I 0 3
        + τ ^ 5 * (iteratedDeriv 3 (Y04.G k) 0 / 6) * Y04.I 3 1
        + τ ^ 5 * (iteratedDeriv 2 (Y04.G k) 0 / 4) * Y04.I 2 2
        + τ ^ 5 * (deriv (Y04.G k) 0 / 6) * Y04.I 1 3 + τ ^ 5 * (Y04.G k 0 / 24) * Y04.I 0 4
        + τ ^ 6 * C6 + τ ^ 7 * C7 + τ ^ 8 * C8 with hS
  have ha : Y04.a k (Y04.h p) = Y04.a k Y04.hDiamond + 8 / π * X := by rw [hX, ← mul_assoc, ← h2]; ring
  have a1 := Y04.four_k_ne k; have b1 := Y04.four_k_ne' k
  have c1 : (16 * (k : ℝ) ^ 2 - 1) ≠ 0 := by
    rw [show (16 * (k : ℝ) ^ 2 - 1) = (4 * k - 1) * (4 * k + 1) by ring]; exact mul_ne_zero a1 b1
  have hI03 : Y04.I 0 3 = (3 / 4 * (π ^ 4 / 90) - 3 * Y04.I 1 2) / 2 := by linarith [Y04.I0312]
  have hI22 : Y04.I 2 2 = 2 * Z02.z5 + 2 * (π ^ 2 / 6) * Y04.z3 - 45 / 8 * Z02.z5 := by
    have := Z02.D51; rw [Z02.I31] at this; linarith
  have hI13 : Y04.I 1 3 = 6 * (-(3 / 4) * Z02.z5 - 7 / 4 * ((π ^ 2 / 6) * Y04.z3))
      + 17 / 4 * (45 / 8 * Z02.z5) + 15 / 4 * (2 * Z02.z5 + 2 * (π ^ 2 / 6) * Y04.z3 - 45 / 8 * Z02.z5)
      - Y04.I 0 4 / 2 := by
    have := Z02.D50; rw [Z02.I31, hI22] at this; linarith
  have key : Y04.a k (Y04.h p) - (-1 : ℝ) ^ k * (Real.sqrt 2 / π) *
        (-(4 / (16 * (k : ℝ) ^ 2 - 1)) + (π ^ 2 / 6) * τ ^ 2 + Y04.z3 * τ ^ 3
          - (33 / 16 + 7 * (k : ℝ) ^ 2) * (π ^ 4 / 90) * τ ^ 4
          + (-(3 / 4 + 4 * (k : ℝ) ^ 2) * Z02.z5
              - (7 / 4 + 4 * (k : ℝ) ^ 2) * ((π ^ 2 / 6) * Y04.z3)) * τ ^ 5 + Θ * τ ^ 6 + Θ7 * τ ^ 7
          + Θ8 * τ ^ 8)
        = 8 / π * (X - S) := by
    have hd2 : iteratedDeriv 2 (Y04.G k) 0 = 2 * -((-1 : ℝ) ^ k * Real.sqrt 2 * ((k : ℝ) ^ 2 / 2 + 5 / 32)) :=
      Z02.G_d2' k
    rw [ha, h1, hS, hC6v, hC7v, hC8v, hd2, Z02.G_d3', Y04.G_zero k, Y04.G_d1 k, Y04.I01, Y04.I11, Y04.I02,
      Y04.I21, hI03, Z02.I31, hI22, hI13]
    field_simp
    ring
  rw [key, abs_mul, abs_of_pos (by positivity)]
  calc 8 / π * |X - S| ≤ 8 / π * (CW * τ ^ 9) := mul_le_mul_of_nonneg_left hb (by positivity)
    _ = _ := by ring

end R02

theorem R02.proof_1 : R02.stmt_1 := R02.zeta8

theorem R02.proof_2 : R02.stmt_2 := R02.G_d6

theorem R02.proof_3 : R02.stmt_3 := fun Gf hG hb => R02.watson8 Gf hG hb

theorem R02.proof_4 : R02.stmt_4 := R02.D83

theorem R02.proof_5 : R02.stmt_5 := R02.D82

theorem R02.proof_6 : R02.stmt_6 := R02.D81

theorem R02.proof_7 : R02.stmt_7 := R02.D80

theorem R02.proof_8 : R02.stmt_8 := R02.theta8

theorem R02.proof_9 : R02.stmt_9 := by
  have hpi := Real.pi_pos
  have n8 : (-127 * (π ^ 8 / 9450) : ℝ) ≠ 0 := by
    have : 0 < π ^ 8 / 9450 := by positivity
    linarith
  refine ⟨Polynomial.C (-127 * (π ^ 8 / 9450)) * Polynomial.X ^ 3
      + Polynomial.C (-(60889 / 48) * (π ^ 8 / 9450) + 16 * (Y04.z3 * Z02.z5)
          + 8 * ((π ^ 2 / 6) * Y04.z3 ^ 2)) * Polynomial.X ^ 2
      + Polynomial.C (-(495965 / 256) * (π ^ 8 / 9450) + 50 * (Y04.z3 * Z02.z5)
          + 27 * ((π ^ 2 / 6) * Y04.z3 ^ 2)) * Polynomial.X
      + Polynomial.C (-(2891923 / 12288) * (π ^ 8 / 9450) + 141 / 16 * (Y04.z3 * Z02.z5)
          + 169 / 32 * ((π ^ 2 / 6) * Y04.z3 ^ 2)), ?_, ?_, ?_⟩
  · rw [Polynomial.degree_cubic n8]; rfl
  · rw [Polynomial.leadingCoeff_cubic n8]
    norm_num
  · intro k
    obtain ⟨C, τ₀, hτ₀, hC⟩ := R02.fourier8 k
    refine ⟨C, τ₀, hτ₀, fun p hp hτ => ?_⟩
    have := hC p hp hτ
    simp only [Polynomial.eval_add, Polynomial.eval_mul, Polynomial.eval_C, Polynomial.eval_X,
      Polynomial.eval_pow]
    convert this using 4
    ring

theorem R02.proof_10 : R02.stmt_10 := R02.fourier8
