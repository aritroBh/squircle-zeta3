import Mathlib
import Zeta3.Wave5.Y04_Fourier
import Zeta3.Wave6.Z02_Theta5
import Zeta3.Wave7.V02_Theta6
/-!
# Q03 — the `τ⁷` coefficient of the Fourier coefficients `a_k(h_p)`: `Θ₇`, and `prop:degree` at `n = 7`

Target id: Q03.
Paper labels: `prop:theta` and `prop:degree` (main.tex lines 1214–1225), extending
`thm:fourier` (main.tex 1204–1211, formalised in `Zeta3.Wave5.Y04_Fourier`) and its `τ⁵`, `τ⁶`
extensions (`Zeta3.Wave6.Z02_Theta5`, `Zeta3.Wave7.V02_Theta6`) by one more order.
Detailed proofs: paper/ieee/proofs/sec_g2_theta.tex (the display of `Θ₇` at line 18,
Lemma `g2:lem:t:red`, Table `g2:tab:D` rows `D̃_{7,0}`, `D̃_{7,1}`, `D̃_{7,2}`) and
paper/ieee/proofs/sec_g2_degree.tex (Lemmas `g2:lem:d:coef`, `g2:lem:d:e`, `g2:lem:d:theta`,
proof of Proposition 2, case `n = 2M+1`); the Watson expansion is Lemma `g2:lem:f:watson` of
sec_g2_fourier.tex.

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

Supplement sec_g2_theta.tex line 18:

  \Theta_7&=\bigl(\tfrac{129}{16}+46u+16u^2\bigr)\zeta(7)+\bigl(\tfrac{141}{16}+50u+16u^2\bigr)\zeta(2)\zeta(5)+\bigl(\tfrac{1071}{160}+\tfrac{181}5u+\tfrac{56}5u^2\bigr)\zeta(2)^2\zeta(3),

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
  $1,\frac12,-\frac{31}{24},-\frac{11}{16},\frac{409}{384}$".
  The supplement lists `e_{0,m}` only for `m ≤ 4`; the next value `e_{0,5} = 721/3840` was
  computed here from the definition (SymPy series of `2^{3/2} H`, cross-checked by numerical
  differentiation of `G_k`). So
  `ρ₅(u) = e_{0,5} - 2u e_{1,5} + (2/3)u² e_{2,5} = 721/3840 + (25/24)u + (1/3)u²`.

Supplement, Lemma `g2:lem:t:red` and Table `g2:tab:D` (products `ζ(7), ζ(2)ζ(5), ζ(2)²ζ(3)`):

  $\wD_{n,m}=2\sum_{j=1}^{n-1-2m}\dfrac{e_{m,\,n-1-j}}{j!}\,I_{n-1-j,\,j}$
  $\Theta_n(u)=\sum_{m=0}^{\lfloor n/2\rfloor-1}\wD_{n,m}\,\dfrac{(-4u)^m}{(2m)!}$.
  $\wD_{7,0}$ & $\tfrac{129}{16},\ \tfrac{141}{16},\ \tfrac{1071}{160}$ \\
  $\wD_{7,1}$ & $-23,\ -25,\ -\tfrac{181}{10}$ \\
  $\wD_{7,2}$ & $24,\ 24,\ \tfrac{84}5$ \\

Supplement, proof of Proposition 2, case `n = 2M+1` (here `M = 3`):

  The coefficient of $u^{M-1}$ is
  2\Bigl[\frac{(-4)^{M-1}}{2(2M-2)!}I_{2M-1,1}+\frac{(-4)^{M-1}}{(2M-2)!}\frac{I_{2M-2,2}}{2!}\Bigr]=\frac{(-4)^{M-1}}{(2M-2)!}\bigl(I_{2M-1,1}+I_{2M-2,2}\bigr).
  Both integrals are positive (positive integrands), so the coefficient is not $0$ and has the
  sign $(-1)^{M-1}$. So $\deg\Theta_{2M+1}=M-1$ exactly.

Supplement, Lemma `g2:lem:f:atoms`: `I_{m,1} = m! η(m+2)`; at `m = 5`,
`I_{5,1} = 5! (1 - 2⁻⁶) ζ(7) = (945/8) ζ(7)`.

Reused (not redefined): from `Zeta3.Wave5.Y04_Fourier`: `Y04.tau` (`τ = 1 - 1/p`), `Y04.h`
(the support function `h_p`), `Y04.a` (`a_k(f) = (4/π)∫_0^{π/2} f cos 4kθ`), `Y04.G` (the layer
function `G_k(x) = cos(4kθ_x) cos²θ_x sin θ_x`, `θ_x = arctan e^{-x}`), `Y04.I`
(`I_{m,j} = ∫_0^∞ s^m log^j(1+e^{-s}) ds`), `Y04.z3` (`ζ(3) = (riemannZeta 3).re`); from
`Zeta3.Wave6.Z02_Theta5`: `Z02.z5` (`ζ(5) = (riemannZeta 5).re`). V02 proves `g_{k,4}`, Watson
at `N = 6`, the weight-six rows and the `τ⁶` expansion (`V02.stmt_8`); Z02 the same one order
lower. New here: `Q03.z7`.

Conventions.
* `ζ(3) = Y04.z3`, `ζ(5) = Z02.z5`, `ζ(7) = Q03.z7 = (riemannZeta 7).re`, the real parts of
  Mathlib's complex zeta values; by `zeta_nat_eq_tsum_of_gt_one` these are the real numbers
  `∑ 1/n^3`, `∑ 1/n^5`, `∑ 1/n^7`. The even values are written out: `ζ(2) = π²/6`,
  `ζ(4) = π⁴/90`, `ζ(6) = π⁶/945`; `ζ(2)² ζ(3)` is `(π²/6)² ζ(3)`.
* `g_{k,m} = G_k^{(m)}(0)/m!` is written `iteratedDeriv m (Y04.G k) 0 / m!`; `G_k` is `C^∞`
  on `ℝ` (proved in Y04), so these are genuine derivatives.
* All half-line integrals `I_{m,j}`, `j ≥ 1`, have integrable integrands (proved in Y04), and
  `stmt_2` asserts the Watson integrand integrable, so no statement holds by the
  "non-integrable ↦ 0" convention.
* As in Y04/Z02/V02, the Watson lemma (`stmt_2`) asks `G` to be `C^∞` on all of `ℝ` and bounded
  on `[0,∞)` (the supplement: `C^∞` on `[0,∞)`); this covers `G = G_k`.
* "The coefficient of `τ⁷`" in `prop:degree` is the coefficient of an asymptotic expansion
  (the supplement's "Setting"). `stmt_8` states it that way: a polynomial `Θ₇ ∈ ℝ[X]` such
  that `(-1)^k (√2/π)[-4/(16k²-1) + Σ_{n=2}^{6} Θ_n(k²) τ^n + Θ₇(k²) τ⁷]` (with the proved
  `Θ₂, …, Θ₆` of Y04/Z02/V02 written out) approximates `a_k(h_p)` to `O_k(τ⁸)` for every `k`.
  Such coefficients are unique and a polynomial is determined by its values at `k²`, so this
  pins down the paper's `Θ₇`. "Degree exactly 2" is `Polynomial.degree Θ₇ = 2`. `prop:degree`
  gives no leading coefficient for odd `n`; `stmt_8` uses the supplement's odd-case formula
  `(-4)^{M-1}/(2M-2)! · (I_{2M-1,1} + I_{2M-2,2})` with `M = 3`.
* The `O(τ⁸)` claims are explicit: `∃ C τ₀, 0 < τ₀ ∧ ∀ p > 1, τ ≤ τ₀ → |error| ≤ C τ⁸`, the
  same form as `Y04.stmt_8`, `Z02.stmt_8`, `V02.stmt_8`.
* Numerical check before freezing (mpmath, 40 digits; SymPy for `e_{r,m}`): `I_{5,1} = 945ζ(7)/8`,
  the rows `D̃_{7,0}, D̃_{7,1}, D̃_{7,2}` (`stmt_4`–`stmt_6`), `g_{k,5}` for `k = 0..3` and `stmt_7`
  for `k = 0..3` agree to `≤ 7·10⁻³⁵`; the leading coefficient is `79.8528…`, as in the
  supplement's check. Direct quadrature of `a_k(h_p)` (50 digits) at `τ = 0.02, 0.01, 0.005` for
  `k = 0, 1, 2` gives `(residual after τ⁶)/τ⁷ - Θ₇(k²) ≈ c·τ` with `c` converging
  (≈ -217, -3457, -36015 at `τ = 0.005`), consistent with an `O(τ⁸)` error.

NOT COVERED: `prop:degree` for `n ≥ 8` and in general, and its clause "coefficients are
alternating multiple zeta values of weight `n`" (no AMZV/iterated-integral library is set up;
for `n = 7` the explicit closed form in `ζ(7), ζ(2)ζ(5), ζ(2)²ζ(3)` is stated instead); `Θ₈, Θ₉`
and `prop:theta` for `n = 8, 9`; the iterated-integral form of `I_{m,j}`, the
shuffle/stuffle/duality relations and the exact certificates (the three `n = 7` identities
`D̃_{7,0}`, `D̃_{7,1}`, `D̃_{7,2}` are stated directly); the general formula for `e_{r,m}`,
`ρ_m` (only `g_{k,5}` is stated); the K3 identification of the odd top coefficient; Cauchy's
formula `P = ∫ h_p` and the `k = 0` comparison with `thm:pw`.
-/

-- BEGIN STATEMENTS
namespace Q03

open Real MeasureTheory Set Polynomial

/-- `ζ(7)` as a real number: the real part of Mathlib's `riemannZeta 7`. -/
noncomputable def z7 : ℝ := (riemannZeta 7).re

/-- Lemma `g2:lem:d:e`(b) at `m = 5` (with `ρ₅(u) = 721/3840 + (25/24)u + (1/3)u²`):
the next Taylor coefficient of `G_k`,
`g_{k,5} = G_k^{(5)}(0)/5! = (-1)^k (√2/4)(721/3840 + (25/24)k² + (1/3)k⁴)`. -/
def stmt_1 : Prop :=
  ∀ k : ℕ, iteratedDeriv 5 (Y04.G k) 0 / 120 =
    (-1 : ℝ) ^ k * Real.sqrt 2 / 4 *
      (721 / 3840 + 25 / 24 * (k : ℝ) ^ 2 + 1 / 3 * (k : ℝ) ^ 4)

/-- Lemma `g2:lem:f:watson` with `N = 7`: for `G` smooth and bounded on `[0,∞)`, with
`g_m = G^{(m)}(0)/m!`, there is `C` with, for `0 < τ < 1`, the integrand integrable and
`|τ ∫_0^∞ [(1+e^{-s})^τ - 1] G(τ s) ds - Σ_{n=2}^{7} τ^n Σ_{j=1}^{n-1} g_{n-1-j}/j! I_{n-1-j,j}|
  ≤ C τ⁸`. -/
def stmt_2 : Prop :=
  ∀ Gf : ℝ → ℝ, ContDiff ℝ (⊤ : ℕ∞) Gf → (∃ M : ℝ, ∀ x : ℝ, 0 ≤ x → |Gf x| ≤ M) →
    ∃ C : ℝ, ∀ τ : ℝ, 0 < τ → τ < 1 →
      IntegrableOn (fun s : ℝ => ((1 + Real.exp (-s)) ^ τ - 1) * Gf (τ * s)) (Ioi 0) ∧
      |τ * (∫ s in Ioi (0 : ℝ), ((1 + Real.exp (-s)) ^ τ - 1) * Gf (τ * s)) -
          ∑ n ∈ Finset.Icc 2 7, τ ^ n * ∑ j ∈ Finset.Icc 1 (n - 1),
            iteratedDeriv (n - 1 - j) Gf 0 / ((n - 1 - j).factorial : ℝ) /
              (j.factorial : ℝ) * Y04.I (n - 1 - j) j| ≤ C * τ ^ 8

/-- Lemma `g2:lem:f:atoms` at `m = 5`: `I_{5,1} = 5! η(7) = (945/8) ζ(7)`. -/
def stmt_3 : Prop := Y04.I 5 1 = 945 / 8 * z7

/-- Table `g2:tab:D`, row `D̃_{7,2}` (with `e_{2,4..5} = 1, 1/2`):
`D̃_{7,2} = 2(e_{2,5} I_{5,1} + e_{2,4} I_{4,2}/2!) = 24 ζ(7) + 24 ζ(2) ζ(5) + (84/5) ζ(2)² ζ(3)`. -/
def stmt_4 : Prop :=
  2 * ((1 / 2) * Y04.I 5 1 + Y04.I 4 2 / 2) =
    24 * z7 + 24 * ((π ^ 2 / 6) * Z02.z5) + 84 / 5 * ((π ^ 2 / 6) ^ 2 * Y04.z3)

/-- Table `g2:tab:D`, row `D̃_{7,1}`: with `e_{1,2..5} = 1, 1/2, -23/24, -25/48`,
`D̃_{7,1} = 2(e_{1,5} I_{5,1} + e_{1,4} I_{4,2}/2! + e_{1,3} I_{3,3}/3! + e_{1,2} I_{2,4}/4!)
  = -23 ζ(7) - 25 ζ(2) ζ(5) - (181/10) ζ(2)² ζ(3)`. -/
def stmt_5 : Prop :=
  2 * (-(25 / 48) * Y04.I 5 1 + (-(23 / 24)) * Y04.I 4 2 / 2 + (1 / 2) * Y04.I 3 3 / 6
      + Y04.I 2 4 / 24) =
    -23 * z7 - 25 * ((π ^ 2 / 6) * Z02.z5) - 181 / 10 * ((π ^ 2 / 6) ^ 2 * Y04.z3)

/-- Table `g2:tab:D`, row `D̃_{7,0}`: with `e_{0,0..5} = 1, 1/2, -5/8, -17/48, 121/384, 721/3840`,
`D̃_{7,0} = 2(e_{0,5} I_{5,1} + e_{0,4} I_{4,2}/2! + e_{0,3} I_{3,3}/3! + e_{0,2} I_{2,4}/4!
  + e_{0,1} I_{1,5}/5! + e_{0,0} I_{0,6}/6!)
  = (129/16) ζ(7) + (141/16) ζ(2) ζ(5) + (1071/160) ζ(2)² ζ(3)`. -/
def stmt_6 : Prop :=
  2 * (721 / 3840 * Y04.I 5 1 + 121 / 384 * Y04.I 4 2 / 2 + (-(17 / 48)) * Y04.I 3 3 / 6
      + (-(5 / 8)) * Y04.I 2 4 / 24 + (1 / 2) * Y04.I 1 5 / 120 + Y04.I 0 6 / 720) =
    129 / 16 * z7 + 141 / 16 * ((π ^ 2 / 6) * Z02.z5) + 1071 / 160 * ((π ^ 2 / 6) ^ 2 * Y04.z3)

/-- Lemmas `g2:lem:d:coef`, `g2:lem:d:theta` and `g2:lem:t:red` at `n = 7`, together with the
closed form of `Θ₇` (sec_g2_theta.tex line 18): the `τ⁷` coefficient
`c_7(k) = (8/π) Σ_{j=1}^{6} g_{k,6-j}/j! I_{6-j,j}` satisfies
`(-1)^k (π/√2) c_7(k) = Θ₇(k²) = (129/16 + 46k² + 16k⁴) ζ(7) + (141/16 + 50k² + 16k⁴) ζ(2) ζ(5)
  + (1071/160 + (181/5)k² + (56/5)k⁴) ζ(2)² ζ(3)`. -/
def stmt_7 : Prop :=
  ∀ k : ℕ,
    (-1 : ℝ) ^ k * (π / Real.sqrt 2) *
        (8 / π * ∑ j ∈ Finset.Icc 1 6,
          iteratedDeriv (6 - j) (Y04.G k) 0 / ((6 - j).factorial : ℝ) / (j.factorial : ℝ) *
            Y04.I (6 - j) j) =
      (129 / 16 + 46 * (k : ℝ) ^ 2 + 16 * (k : ℝ) ^ 4) * z7
        + (141 / 16 + 50 * (k : ℝ) ^ 2 + 16 * (k : ℝ) ^ 4) * ((π ^ 2 / 6) * Z02.z5)
        + (1071 / 160 + 181 / 5 * (k : ℝ) ^ 2 + 56 / 5 * (k : ℝ) ^ 4)
            * ((π ^ 2 / 6) ^ 2 * Y04.z3)

/-- `prop:degree` at `n = 7 = 2M+1`, `M = 3`: there is a polynomial `Θ₇` of degree exactly
`⌊7/2⌋ - 1 = 2`, with leading coefficient `(-4)^{M-1}/(2M-2)! · (I_{5,1} + I_{4,2})` (the
supplement's odd case), whose value at `u = k²` is the `τ⁷` coefficient of
`(-1)^k (π/√2) a_k(h_p)` in the asymptotic sense (the lower coefficients `Θ₂, …, Θ₆` written out
as proved in Y04, Z02 and V02): for every `k` there are `C_k`, `τ₀ > 0` with
`|a_k(h_p) - (-1)^k (√2/π)[-4/(16k²-1) + ζ(2)τ² + ζ(3)τ³ - (33/16 + 7k²)ζ(4)τ⁴ + Θ₅(k²)τ⁵
   + Θ₆(k²)τ⁶ + Θ₇(k²)τ⁷]| ≤ C_k τ⁸` for `p > 1`, `τ ≤ τ₀`. -/
def stmt_8 : Prop :=
  ∃ Θ₇ : ℝ[X],
    Θ₇.degree = ((7 / 2 - 1 : ℕ) : WithBot ℕ) ∧
    Θ₇.leadingCoeff =
      (-4 : ℝ) ^ (3 - 1) / ((2 * 3 - 2).factorial : ℝ) * (Y04.I 5 1 + Y04.I 4 2) ∧
    ∀ k : ℕ, ∃ C τ₀ : ℝ, 0 < τ₀ ∧ ∀ p : ℝ, 1 < p → Y04.tau p ≤ τ₀ →
      |Y04.a k (Y04.h p) - (-1 : ℝ) ^ k * (Real.sqrt 2 / π) *
          (-(4 / (16 * (k : ℝ) ^ 2 - 1)) + (π ^ 2 / 6) * Y04.tau p ^ 2 + Y04.z3 * Y04.tau p ^ 3
            - (33 / 16 + 7 * (k : ℝ) ^ 2) * (π ^ 4 / 90) * Y04.tau p ^ 4
            + (-(3 / 4 + 4 * (k : ℝ) ^ 2) * Z02.z5
                - (7 / 4 + 4 * (k : ℝ) ^ 2) * ((π ^ 2 / 6) * Y04.z3)) * Y04.tau p ^ 5
            + ((4039 / 256 + 733 / 8 * (k : ℝ) ^ 2 + 31 * (k : ℝ) ^ 4) * (π ^ 6 / 945)
                - (7 / 8 + 2 * (k : ℝ) ^ 2) * Y04.z3 ^ 2) * Y04.tau p ^ 6
            + Θ₇.eval ((k : ℝ) ^ 2) * Y04.tau p ^ 7)|
        ≤ C * Y04.tau p ^ 8

/-- `thm:fourier` extended by three orders with `prop:theta`'s `Θ₅` and the supplement's `Θ₆`,
`Θ₇`: for every `k ≥ 0` there are `C_k` and `τ₀ > 0` such that for every `p > 1` with
`τ = 1 - 1/p ≤ τ₀`, with `u = k²`,
`|a_k(h_p) - (-1)^k (√2/π)[-4/(16k²-1) + ζ(2)τ² + ζ(3)τ³ - (33/16 + 7k²)ζ(4)τ⁴ + Θ₅(u)τ⁵
   + Θ₆(u)τ⁶ + Θ₇(u)τ⁷]| ≤ C_k τ⁸`, where `Θ₅(u) = -(3/4 + 4u) ζ(5) - (7/4 + 4u) ζ(2) ζ(3)`,
`Θ₆(u) = (4039/256 + (733/8)u + 31u²) ζ(6) - (7/8 + 2u) ζ(3)²` and
`Θ₇(u) = (129/16 + 46u + 16u²) ζ(7) + (141/16 + 50u + 16u²) ζ(2) ζ(5)
  + (1071/160 + (181/5)u + (56/5)u²) ζ(2)² ζ(3)`. -/
def stmt_9 : Prop :=
  ∀ k : ℕ, ∃ C τ₀ : ℝ, 0 < τ₀ ∧ ∀ p : ℝ, 1 < p → Y04.tau p ≤ τ₀ →
    |Y04.a k (Y04.h p) - (-1 : ℝ) ^ k * (Real.sqrt 2 / π) *
        (-(4 / (16 * (k : ℝ) ^ 2 - 1)) + (π ^ 2 / 6) * Y04.tau p ^ 2 + Y04.z3 * Y04.tau p ^ 3
          - (33 / 16 + 7 * (k : ℝ) ^ 2) * (π ^ 4 / 90) * Y04.tau p ^ 4
          + (-(3 / 4 + 4 * (k : ℝ) ^ 2) * Z02.z5
              - (7 / 4 + 4 * (k : ℝ) ^ 2) * ((π ^ 2 / 6) * Y04.z3)) * Y04.tau p ^ 5
          + ((4039 / 256 + 733 / 8 * (k : ℝ) ^ 2 + 31 * (k : ℝ) ^ 4) * (π ^ 6 / 945)
              - (7 / 8 + 2 * (k : ℝ) ^ 2) * Y04.z3 ^ 2) * Y04.tau p ^ 6
          + ((129 / 16 + 46 * (k : ℝ) ^ 2 + 16 * (k : ℝ) ^ 4) * z7
              + (141 / 16 + 50 * (k : ℝ) ^ 2 + 16 * (k : ℝ) ^ 4) * ((π ^ 2 / 6) * Z02.z5)
              + (1071 / 160 + 181 / 5 * (k : ℝ) ^ 2 + 56 / 5 * (k : ℝ) ^ 4)
                  * ((π ^ 2 / 6) ^ 2 * Y04.z3)) * Y04.tau p ^ 7)|
      ≤ C * Y04.tau p ^ 8

end Q03
-- END STATEMENTS

open Real MeasureTheory Set

namespace Q03

set_option maxRecDepth 20000 in
set_option maxHeartbeats 4000000 in
lemma G_d5 (k : ℕ) : iteratedDeriv 5 (Y04.G k) 0 / 120 =
    (-1 : ℝ) ^ k * Real.sqrt 2 / 4 * (721 / 3840 + 25 / 24 * (k : ℝ) ^ 2 + 1 / 3 * (k : ℝ) ^ 4) := by
  rw [V02.iter_G]
  simp only [Y04.thetaX_zero, Function.iterate_succ_apply', Function.iterate_zero_apply, V02.step, V02.dterm,
    V02.shift, List.flatMap_cons, List.flatMap_nil, List.map_cons, List.map_nil, List.append_nil,
    List.cons_append, List.nil_append, ite_true, ite_false, Bool.false_eq_true]
  simp only [V02.ev, V02.tm, List.map_cons, List.map_nil, List.sum_cons, List.sum_nil, ite_true, ite_false,
    Bool.false_eq_true, Y04.cos_k, Y04.sin_k, Real.cos_pi_div_four, Real.sin_pi_div_four]
  norm_num
  have : Real.sqrt 2 ^ 2 = 2 := Real.sq_sqrt (by norm_num)
  have h13 : Real.sqrt 2 ^ 13 = 64 * Real.sqrt 2 := by
    rw [show Real.sqrt 2 ^ 13 = (Real.sqrt 2 ^ 2) ^ 6 * Real.sqrt 2 by ring, this]; ring
  ring_nf
  rw [h13]; ring

lemma zn_seven : Y04.zn 7 = z7 := by simp [Y04.zn, z7]

lemma I51 : Y04.I 5 1 = 945 / 8 * z7 := by
  rw [Y04.I_m1, Y04.eta_val 7 (by norm_num), zn_seven]; norm_num [Nat.factorial]; ring


/-- Taylor's formula at `0` to order 5 with a sextic remainder on `[0,∞)`. -/
lemma taylor5 {Gf : ℝ → ℝ} (hG : ContDiff ℝ (⊤ : ℕ∞) Gf) {M : ℝ}
    (hM : ∀ x : ℝ, 0 ≤ x → |Gf x| ≤ M) :
    ∃ C : ℝ, 0 ≤ C ∧ ∀ x : ℝ, 0 ≤ x →
      |Gf x - (iteratedDeriv 0 Gf 0 / (Nat.factorial 0 : ℝ) + iteratedDeriv 1 Gf 0 / (Nat.factorial 1 : ℝ) * x + iteratedDeriv 2 Gf 0 / (Nat.factorial 2 : ℝ) * x ^ 2 + iteratedDeriv 3 Gf 0 / (Nat.factorial 3 : ℝ) * x ^ 3 + iteratedDeriv 4 Gf 0 / (Nat.factorial 4 : ℝ) * x ^ 4 + iteratedDeriv 5 Gf 0 / (Nat.factorial 5 : ℝ) * x ^ 5)| ≤ C * x ^ 6 := by
  set g0 := iteratedDeriv 0 Gf 0 / (Nat.factorial 0 : ℝ)
  set g1 := iteratedDeriv 1 Gf 0 / (Nat.factorial 1 : ℝ)
  set g2 := iteratedDeriv 2 Gf 0 / (Nat.factorial 2 : ℝ)
  set g3 := iteratedDeriv 3 Gf 0 / (Nat.factorial 3 : ℝ)
  set g4 := iteratedDeriv 4 Gf 0 / (Nat.factorial 4 : ℝ)
  set g5 := iteratedDeriv 5 Gf 0 / (Nat.factorial 5 : ℝ)
  have hc6 : Continuous (iteratedDeriv 6 Gf) :=
    hG.continuous_iteratedDeriv 6 (WithTop.coe_le_coe.2 le_top)
  obtain ⟨K6, hK6⟩ := (isCompact_Icc (a := (0 : ℝ)) (b := 1)).exists_bound_of_continuousOn
    hc6.continuousOn
  have hM0 : 0 ≤ M := (abs_nonneg _).trans (hM 0 le_rfl)
  refine ⟨max (|K6| / 720) (M + |g0| + |g1| + |g2| + |g3| + |g4| + |g5|), le_max_of_le_left (by positivity), ?_⟩
  intro x hx
  rcases eq_or_lt_of_le hx with h0 | hpos
  · subst h0; simp [g0]
  rcases le_or_gt x 1 with hx1 | hx1
  · have hcd : ContDiffOn ℝ (5 + 1) Gf (uIcc 0 x) := (hG.of_le (WithTop.coe_le_coe.2 le_top)).contDiffOn
    obtain ⟨x', hx', heq⟩ := taylor_mean_remainder_lagrange_iteratedDeriv (n := 5) hpos.ne hcd
    have hu : UniqueDiffOn ℝ (uIcc 0 x) := uniqueDiffOn_uIcc hpos.ne
    have hmem : (0 : ℝ) ∈ uIcc 0 x := left_mem_uIcc
    have hT : taylorWithinEval Gf 5 (uIcc 0 x) 0 x
        = g0 + g1 * x + g2 * x ^ 2 + g3 * x ^ 3 + g4 * x ^ 4 + g5 * x ^ 5 := by
      rw [taylor_within_apply]
      simp only [Finset.sum_range_succ, Finset.sum_range_zero, zero_add, sub_zero, smul_eq_mul]
      rw [iteratedDerivWithin_eq_iteratedDeriv hu (hG.contDiffAt.of_le (WithTop.coe_le_coe.2 le_top)) hmem,
        iteratedDerivWithin_eq_iteratedDeriv hu (hG.contDiffAt.of_le (WithTop.coe_le_coe.2 le_top)) hmem,
        iteratedDerivWithin_eq_iteratedDeriv hu (hG.contDiffAt.of_le (WithTop.coe_le_coe.2 le_top)) hmem,
        iteratedDerivWithin_eq_iteratedDeriv hu (hG.contDiffAt.of_le (WithTop.coe_le_coe.2 le_top)) hmem,
        iteratedDerivWithin_eq_iteratedDeriv hu (hG.contDiffAt.of_le (WithTop.coe_le_coe.2 le_top)) hmem,
        iteratedDerivWithin_eq_iteratedDeriv hu (hG.contDiffAt.of_le (WithTop.coe_le_coe.2 le_top)) hmem]
      simp only [g0, g1, g2, g3, g4, g5]; norm_num; ring
    rw [hT] at heq
    rw [heq]
    rw [uIoo_of_le hx] at hx'
    have hb := hK6 x' ⟨hx'.1.le, by linarith [hx'.2]⟩
    rw [Real.norm_eq_abs] at hb
    rw [sub_zero, abs_div, abs_mul, abs_pow, abs_of_pos hpos]
    norm_num [Nat.factorial]
    have : |iteratedDeriv 6 Gf x'| * x ^ 6 / 720 ≤ |K6| / 720 * x ^ 6 := by
      have := hb.trans (le_abs_self K6)
      have h3 : 0 ≤ x ^ 6 := by positivity
      nlinarith
    calc |iteratedDeriv 6 Gf x'| * x ^ 6 / 720 ≤ |K6| / 720 * x ^ 6 := this
      _ ≤ _ := by gcongr; exact le_max_left _ _
  · have h1 : 1 ≤ x := hx1.le
    have hx2' : 1 ≤ x ^ 2 := by nlinarith
    have hx3' : 1 ≤ x ^ 3 := by nlinarith
    have hx4' : 1 ≤ x ^ 4 := by nlinarith
    have hx5' : 1 ≤ x ^ 5 := by nlinarith
    have hx16 : x ≤ x ^ 6 := by
      rw [show x ^ 6 = x * x ^ 5 by ring]; nlinarith
    have hx26 : x ^ 2 ≤ x ^ 6 := by
      rw [show x ^ 6 = x ^ 2 * x ^ 4 by ring]; nlinarith
    have hx36 : x ^ 3 ≤ x ^ 6 := by
      rw [show x ^ 6 = x ^ 3 * x ^ 3 by ring]; nlinarith
    have hx46 : x ^ 4 ≤ x ^ 6 := by
      rw [show x ^ 6 = x ^ 4 * x ^ 2 by ring]; nlinarith
    have hx56 : x ^ 5 ≤ x ^ 6 := by
      rw [show x ^ 6 = x ^ 5 * x by ring]; nlinarith
    have h16 : 1 ≤ x ^ 6 := by linarith
    calc _ ≤ |Gf x| + (|g0| + |g1| * x + |g2| * x ^ 2 + |g3| * x ^ 3 + |g4| * x ^ 4 + |g5| * x ^ 5) := by
          have := abs_sub (Gf x) (g0 + g1 * x + g2 * x ^ 2 + g3 * x ^ 3 + g4 * x ^ 4 + g5 * x ^ 5)
          have a000 := abs_add_le (g0 + g1 * x + g2 * x ^ 2 + g3 * x ^ 3 + g4 * x ^ 4) (g5 * x ^ 5)
          have a00 := abs_add_le (g0 + g1 * x + g2 * x ^ 2 + g3 * x ^ 3) (g4 * x ^ 4)
          have a0 := abs_add_le (g0 + g1 * x + g2 * x ^ 2) (g3 * x ^ 3)
          have a1 := abs_add_le (g0 + g1 * x) (g2 * x ^ 2)
          have a2 := abs_add_le g0 (g1 * x)
          rw [abs_mul, abs_of_pos hpos] at a2
          rw [abs_mul, abs_of_nonneg (by positivity : (0:ℝ) ≤ x ^ 2)] at a1
          rw [abs_mul, abs_of_nonneg (by positivity : (0:ℝ) ≤ x ^ 3)] at a0
          rw [abs_mul, abs_of_nonneg (by positivity : (0:ℝ) ≤ x ^ 4)] at a00
          rw [abs_mul, abs_of_nonneg (by positivity : (0:ℝ) ≤ x ^ 5)] at a000
          linarith
      _ ≤ (M + |g0| + |g1| + |g2| + |g3| + |g4| + |g5|) * x ^ 6 := by
          have := hM x hx
          have e1 : |g1| * x ≤ |g1| * x ^ 6 := mul_le_mul_of_nonneg_left hx16 (abs_nonneg _)
          have e2 : |g2| * x ^ 2 ≤ |g2| * x ^ 6 := mul_le_mul_of_nonneg_left hx26 (abs_nonneg _)
          have e5 : |g3| * x ^ 3 ≤ |g3| * x ^ 6 := mul_le_mul_of_nonneg_left hx36 (abs_nonneg _)
          have e6 : |g4| * x ^ 4 ≤ |g4| * x ^ 6 := mul_le_mul_of_nonneg_left hx46 (abs_nonneg _)
          have e7 : |g5| * x ^ 5 ≤ |g5| * x ^ 6 := mul_le_mul_of_nonneg_left hx56 (abs_nonneg _)
          have e3 : M ≤ M * x ^ 6 := le_mul_of_one_le_right hM0 h16
          have e4 : |g0| ≤ |g0| * x ^ 6 := le_mul_of_one_le_right (abs_nonneg _) h16
          nlinarith
      _ ≤ _ := by gcongr; exact le_max_right _ _

/-- the dominating function of the order-7 Watson remainder -/
noncomputable def B8 (s : ℝ) : ℝ := (1 + s + s ^ 2 + s ^ 3 + s ^ 4 + s ^ 5 + s ^ 6) * Real.exp (-s)

lemma B8_int : IntegrableOn B8 (Set.Ioi 0) := by
  have h := ((((((Y04.intOn_pow_exp1 0).add (Y04.intOn_pow_exp1 1)).add (Y04.intOn_pow_exp1 2)).add
    (Y04.intOn_pow_exp1 3)).add (Y04.intOn_pow_exp1 4)).add (Y04.intOn_pow_exp1 5)).add (Y04.intOn_pow_exp1 6)
  refine IntegrableOn.congr_fun h (fun s _ => ?_) measurableSet_Ioi
  simp only [B8, Pi.add_apply]; ring

lemma B8_nonneg {s : ℝ} (hs : 0 ≤ s) : 0 ≤ B8 s := by unfold B8; positivity

lemma hterm {g s τ L c : ℝ} {m j a : ℕ} (hs : 0 ≤ s) (hτ0 : 0 ≤ τ) (hτ1 : τ ≤ 1) (hL0 : 0 ≤ L)
    (hL1 : L ≤ 1) (hLe : L ≤ Real.exp (-s)) (hj : 1 ≤ j) (hc : 1 ≤ c) :
    |g * s ^ m * (τ ^ a * L ^ j / c)| ≤ |g| * (s ^ m * Real.exp (-s)) := by
  have h1 : τ ^ a ≤ 1 := pow_le_one₀ hτ0 hτ1
  have h2 : L ^ j ≤ Real.exp (-s) := (pow_le_of_le_one hL0 hL1 (by omega)).trans hLe
  have h3 : 0 ≤ τ ^ a * L ^ j / c := div_nonneg (by positivity) (by linarith)
  have h4 : τ ^ a * L ^ j / c ≤ Real.exp (-s) := by
    rw [div_le_iff₀ (by linarith)]
    have : τ ^ a * L ^ j ≤ 1 * Real.exp (-s) :=
      mul_le_mul h1 h2 (by positivity) zero_le_one
    have : Real.exp (-s) ≤ Real.exp (-s) * c := le_mul_of_one_le_right (Real.exp_pos _).le hc
    linarith
  rw [abs_mul, abs_mul, abs_of_nonneg (by positivity : (0:ℝ) ≤ s ^ m), abs_of_nonneg h3, mul_assoc]
  exact mul_le_mul_of_nonneg_left (mul_le_mul_of_nonneg_left h4 (by positivity)) (abs_nonneg _)

set_option maxHeartbeats 8000000 in
/-- the pointwise order-7 Watson remainder -/
lemma watson_pt7 {g0 g1 g2 g3 g4 g5 CR M τ s : ℝ} (hCR : 0 ≤ CR) (hM : 0 ≤ M) (hτ0 : 0 < τ) (hτ1 : τ < 1)
    (hs : 0 ≤ s) {Gv : ℝ} (hGv : |Gv| ≤ M)
    (hR : |Gv - (g0 + g1 * (τ * s) + g2 * (τ * s) ^ 2 + g3 * (τ * s) ^ 3 + g4 * (τ * s) ^ 4 + g5 * (τ * s) ^ 5)| ≤ CR * (τ * s) ^ 6) :
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
        + τ ^ 7 * (g0 / 720) * (s ^ 0 * Y04.ell s ^ 6))|
      ≤ τ ^ 8 * ((5 * (|g1| + |g2| + |g3| + |g4| + |g5|) + 2 * CR + M) * B8 s) := by
  have hl0 := Y04.ell_nonneg s
  have hl1 := Y04.ell_le_one s hs
  have hle := Y04.ell_le s
  set L := Y04.ell s
  have hτL : |τ * L| ≤ 1 := by rw [abs_of_nonneg (by positivity)]; nlinarith
  have he := Real.exp_bound hτL (n := 7) (by norm_num)
  simp only [Finset.sum_range_succ, Finset.sum_range_zero, Nat.factorial] at he
  norm_num at he
  set a1 := τ * L + (τ * L) ^ 2 / 2 + (τ * L) ^ 3 / 6 + (τ * L) ^ 4 / 24 + (τ * L) ^ 5 / 120
    + (τ * L) ^ 6 / 720 with ha1
  set r1 := Real.exp (τ * L) - 1 - a1 with hr1
  have hr1b : |r1| ≤ τ ^ 7 * L ^ 7 := by
    have : r1 = Real.exp (τ * L) - (1 + τ * L + (τ * L) ^ 2 / 2 + (τ * L) ^ 3 / 6 + (τ * L) ^ 4 / 24
        + (τ * L) ^ 5 / 120 + (τ * L) ^ 6 / 720) := by
      rw [hr1, ha1]; ring
    rw [this]
    have h4 : (|τ| * |L|) ^ 7 = τ ^ 7 * L ^ 7 := by
      rw [abs_of_pos hτ0, abs_of_nonneg hl0]; ring
    rw [h4] at he
    have : 0 ≤ τ ^ 7 * L ^ 7 := by positivity
    linarith
  set q := g0 + g1 * (τ * s) + g2 * (τ * s) ^ 2 + g3 * (τ * s) ^ 3 + g4 * (τ * s) ^ 4 + g5 * (τ * s) ^ 5
  set R := Gv - q with hRdef
  set Q := g1 * s ^ 1 * (τ ^ 0 * L ^ 6 / 720) + g2 * s ^ 2 * (τ ^ 0 * L ^ 5 / 120) + g2 * s ^ 2 * (τ ^ 1 * L ^ 6 / 720) + g3 * s ^ 3 * (τ ^ 0 * L ^ 4 / 24) + g3 * s ^ 3 * (τ ^ 1 * L ^ 5 / 120) + g3 * s ^ 3 * (τ ^ 2 * L ^ 6 / 720) + g4 * s ^ 4 * (τ ^ 0 * L ^ 3 / 6) + g4 * s ^ 4 * (τ ^ 1 * L ^ 4 / 24) + g4 * s ^ 4 * (τ ^ 2 * L ^ 5 / 120) + g4 * s ^ 4 * (τ ^ 3 * L ^ 6 / 720) + g5 * s ^ 5 * (τ ^ 0 * L ^ 2 / 2) + g5 * s ^ 5 * (τ ^ 1 * L ^ 3 / 6) + g5 * s ^ 5 * (τ ^ 2 * L ^ 4 / 24) + g5 * s ^ 5 * (τ ^ 3 * L ^ 5 / 120) + g5 * s ^ 5 * (τ ^ 4 * L ^ 6 / 720) with hQ
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
        + τ ^ 7 * (g0 / 720) * (s ^ 0 * L ^ 6))
      = τ ^ 8 * Q + τ * a1 * R + τ * r1 * Gv := by
    rw [hRdef, hr1, ha1, hQ]; ring
  rw [key]
  have hes : 0 ≤ Real.exp (-s) := (Real.exp_pos _).le
  have hτ' : 0 ≤ τ := hτ0.le

  have q0 : |g1 * s ^ 1 * (τ ^ 0 * L ^ 6 / 720)| ≤ |g1| * (s ^ 1 * Real.exp (-s)) := hterm hs hτ' hτ1.le hl0 hl1 hle (by norm_num) (by norm_num)
  have q1 : |g2 * s ^ 2 * (τ ^ 0 * L ^ 5 / 120)| ≤ |g2| * (s ^ 2 * Real.exp (-s)) := hterm hs hτ' hτ1.le hl0 hl1 hle (by norm_num) (by norm_num)
  have q2 : |g2 * s ^ 2 * (τ ^ 1 * L ^ 6 / 720)| ≤ |g2| * (s ^ 2 * Real.exp (-s)) := hterm hs hτ' hτ1.le hl0 hl1 hle (by norm_num) (by norm_num)
  have q3 : |g3 * s ^ 3 * (τ ^ 0 * L ^ 4 / 24)| ≤ |g3| * (s ^ 3 * Real.exp (-s)) := hterm hs hτ' hτ1.le hl0 hl1 hle (by norm_num) (by norm_num)
  have q4 : |g3 * s ^ 3 * (τ ^ 1 * L ^ 5 / 120)| ≤ |g3| * (s ^ 3 * Real.exp (-s)) := hterm hs hτ' hτ1.le hl0 hl1 hle (by norm_num) (by norm_num)
  have q5 : |g3 * s ^ 3 * (τ ^ 2 * L ^ 6 / 720)| ≤ |g3| * (s ^ 3 * Real.exp (-s)) := hterm hs hτ' hτ1.le hl0 hl1 hle (by norm_num) (by norm_num)
  have q6 : |g4 * s ^ 4 * (τ ^ 0 * L ^ 3 / 6)| ≤ |g4| * (s ^ 4 * Real.exp (-s)) := hterm hs hτ' hτ1.le hl0 hl1 hle (by norm_num) (by norm_num)
  have q7 : |g4 * s ^ 4 * (τ ^ 1 * L ^ 4 / 24)| ≤ |g4| * (s ^ 4 * Real.exp (-s)) := hterm hs hτ' hτ1.le hl0 hl1 hle (by norm_num) (by norm_num)
  have q8 : |g4 * s ^ 4 * (τ ^ 2 * L ^ 5 / 120)| ≤ |g4| * (s ^ 4 * Real.exp (-s)) := hterm hs hτ' hτ1.le hl0 hl1 hle (by norm_num) (by norm_num)
  have q9 : |g4 * s ^ 4 * (τ ^ 3 * L ^ 6 / 720)| ≤ |g4| * (s ^ 4 * Real.exp (-s)) := hterm hs hτ' hτ1.le hl0 hl1 hle (by norm_num) (by norm_num)
  have q10 : |g5 * s ^ 5 * (τ ^ 0 * L ^ 2 / 2)| ≤ |g5| * (s ^ 5 * Real.exp (-s)) := hterm hs hτ' hτ1.le hl0 hl1 hle (by norm_num) (by norm_num)
  have q11 : |g5 * s ^ 5 * (τ ^ 1 * L ^ 3 / 6)| ≤ |g5| * (s ^ 5 * Real.exp (-s)) := hterm hs hτ' hτ1.le hl0 hl1 hle (by norm_num) (by norm_num)
  have q12 : |g5 * s ^ 5 * (τ ^ 2 * L ^ 4 / 24)| ≤ |g5| * (s ^ 5 * Real.exp (-s)) := hterm hs hτ' hτ1.le hl0 hl1 hle (by norm_num) (by norm_num)
  have q13 : |g5 * s ^ 5 * (τ ^ 3 * L ^ 5 / 120)| ≤ |g5| * (s ^ 5 * Real.exp (-s)) := hterm hs hτ' hτ1.le hl0 hl1 hle (by norm_num) (by norm_num)
  have q14 : |g5 * s ^ 5 * (τ ^ 4 * L ^ 6 / 720)| ≤ |g5| * (s ^ 5 * Real.exp (-s)) := hterm hs hτ' hτ1.le hl0 hl1 hle (by norm_num) (by norm_num)
  have c1 := abs_add_le (g1 * s ^ 1 * (τ ^ 0 * L ^ 6 / 720)) (g2 * s ^ 2 * (τ ^ 0 * L ^ 5 / 120))
  have c2 := abs_add_le (g1 * s ^ 1 * (τ ^ 0 * L ^ 6 / 720) + g2 * s ^ 2 * (τ ^ 0 * L ^ 5 / 120)) (g2 * s ^ 2 * (τ ^ 1 * L ^ 6 / 720))
  have c3 := abs_add_le (g1 * s ^ 1 * (τ ^ 0 * L ^ 6 / 720) + g2 * s ^ 2 * (τ ^ 0 * L ^ 5 / 120) + g2 * s ^ 2 * (τ ^ 1 * L ^ 6 / 720)) (g3 * s ^ 3 * (τ ^ 0 * L ^ 4 / 24))
  have c4 := abs_add_le (g1 * s ^ 1 * (τ ^ 0 * L ^ 6 / 720) + g2 * s ^ 2 * (τ ^ 0 * L ^ 5 / 120) + g2 * s ^ 2 * (τ ^ 1 * L ^ 6 / 720) + g3 * s ^ 3 * (τ ^ 0 * L ^ 4 / 24)) (g3 * s ^ 3 * (τ ^ 1 * L ^ 5 / 120))
  have c5 := abs_add_le (g1 * s ^ 1 * (τ ^ 0 * L ^ 6 / 720) + g2 * s ^ 2 * (τ ^ 0 * L ^ 5 / 120) + g2 * s ^ 2 * (τ ^ 1 * L ^ 6 / 720) + g3 * s ^ 3 * (τ ^ 0 * L ^ 4 / 24) + g3 * s ^ 3 * (τ ^ 1 * L ^ 5 / 120)) (g3 * s ^ 3 * (τ ^ 2 * L ^ 6 / 720))
  have c6 := abs_add_le (g1 * s ^ 1 * (τ ^ 0 * L ^ 6 / 720) + g2 * s ^ 2 * (τ ^ 0 * L ^ 5 / 120) + g2 * s ^ 2 * (τ ^ 1 * L ^ 6 / 720) + g3 * s ^ 3 * (τ ^ 0 * L ^ 4 / 24) + g3 * s ^ 3 * (τ ^ 1 * L ^ 5 / 120) + g3 * s ^ 3 * (τ ^ 2 * L ^ 6 / 720)) (g4 * s ^ 4 * (τ ^ 0 * L ^ 3 / 6))
  have c7 := abs_add_le (g1 * s ^ 1 * (τ ^ 0 * L ^ 6 / 720) + g2 * s ^ 2 * (τ ^ 0 * L ^ 5 / 120) + g2 * s ^ 2 * (τ ^ 1 * L ^ 6 / 720) + g3 * s ^ 3 * (τ ^ 0 * L ^ 4 / 24) + g3 * s ^ 3 * (τ ^ 1 * L ^ 5 / 120) + g3 * s ^ 3 * (τ ^ 2 * L ^ 6 / 720) + g4 * s ^ 4 * (τ ^ 0 * L ^ 3 / 6)) (g4 * s ^ 4 * (τ ^ 1 * L ^ 4 / 24))
  have c8 := abs_add_le (g1 * s ^ 1 * (τ ^ 0 * L ^ 6 / 720) + g2 * s ^ 2 * (τ ^ 0 * L ^ 5 / 120) + g2 * s ^ 2 * (τ ^ 1 * L ^ 6 / 720) + g3 * s ^ 3 * (τ ^ 0 * L ^ 4 / 24) + g3 * s ^ 3 * (τ ^ 1 * L ^ 5 / 120) + g3 * s ^ 3 * (τ ^ 2 * L ^ 6 / 720) + g4 * s ^ 4 * (τ ^ 0 * L ^ 3 / 6) + g4 * s ^ 4 * (τ ^ 1 * L ^ 4 / 24)) (g4 * s ^ 4 * (τ ^ 2 * L ^ 5 / 120))
  have c9 := abs_add_le (g1 * s ^ 1 * (τ ^ 0 * L ^ 6 / 720) + g2 * s ^ 2 * (τ ^ 0 * L ^ 5 / 120) + g2 * s ^ 2 * (τ ^ 1 * L ^ 6 / 720) + g3 * s ^ 3 * (τ ^ 0 * L ^ 4 / 24) + g3 * s ^ 3 * (τ ^ 1 * L ^ 5 / 120) + g3 * s ^ 3 * (τ ^ 2 * L ^ 6 / 720) + g4 * s ^ 4 * (τ ^ 0 * L ^ 3 / 6) + g4 * s ^ 4 * (τ ^ 1 * L ^ 4 / 24) + g4 * s ^ 4 * (τ ^ 2 * L ^ 5 / 120)) (g4 * s ^ 4 * (τ ^ 3 * L ^ 6 / 720))
  have c10 := abs_add_le (g1 * s ^ 1 * (τ ^ 0 * L ^ 6 / 720) + g2 * s ^ 2 * (τ ^ 0 * L ^ 5 / 120) + g2 * s ^ 2 * (τ ^ 1 * L ^ 6 / 720) + g3 * s ^ 3 * (τ ^ 0 * L ^ 4 / 24) + g3 * s ^ 3 * (τ ^ 1 * L ^ 5 / 120) + g3 * s ^ 3 * (τ ^ 2 * L ^ 6 / 720) + g4 * s ^ 4 * (τ ^ 0 * L ^ 3 / 6) + g4 * s ^ 4 * (τ ^ 1 * L ^ 4 / 24) + g4 * s ^ 4 * (τ ^ 2 * L ^ 5 / 120) + g4 * s ^ 4 * (τ ^ 3 * L ^ 6 / 720)) (g5 * s ^ 5 * (τ ^ 0 * L ^ 2 / 2))
  have c11 := abs_add_le (g1 * s ^ 1 * (τ ^ 0 * L ^ 6 / 720) + g2 * s ^ 2 * (τ ^ 0 * L ^ 5 / 120) + g2 * s ^ 2 * (τ ^ 1 * L ^ 6 / 720) + g3 * s ^ 3 * (τ ^ 0 * L ^ 4 / 24) + g3 * s ^ 3 * (τ ^ 1 * L ^ 5 / 120) + g3 * s ^ 3 * (τ ^ 2 * L ^ 6 / 720) + g4 * s ^ 4 * (τ ^ 0 * L ^ 3 / 6) + g4 * s ^ 4 * (τ ^ 1 * L ^ 4 / 24) + g4 * s ^ 4 * (τ ^ 2 * L ^ 5 / 120) + g4 * s ^ 4 * (τ ^ 3 * L ^ 6 / 720) + g5 * s ^ 5 * (τ ^ 0 * L ^ 2 / 2)) (g5 * s ^ 5 * (τ ^ 1 * L ^ 3 / 6))
  have c12 := abs_add_le (g1 * s ^ 1 * (τ ^ 0 * L ^ 6 / 720) + g2 * s ^ 2 * (τ ^ 0 * L ^ 5 / 120) + g2 * s ^ 2 * (τ ^ 1 * L ^ 6 / 720) + g3 * s ^ 3 * (τ ^ 0 * L ^ 4 / 24) + g3 * s ^ 3 * (τ ^ 1 * L ^ 5 / 120) + g3 * s ^ 3 * (τ ^ 2 * L ^ 6 / 720) + g4 * s ^ 4 * (τ ^ 0 * L ^ 3 / 6) + g4 * s ^ 4 * (τ ^ 1 * L ^ 4 / 24) + g4 * s ^ 4 * (τ ^ 2 * L ^ 5 / 120) + g4 * s ^ 4 * (τ ^ 3 * L ^ 6 / 720) + g5 * s ^ 5 * (τ ^ 0 * L ^ 2 / 2) + g5 * s ^ 5 * (τ ^ 1 * L ^ 3 / 6)) (g5 * s ^ 5 * (τ ^ 2 * L ^ 4 / 24))
  have c13 := abs_add_le (g1 * s ^ 1 * (τ ^ 0 * L ^ 6 / 720) + g2 * s ^ 2 * (τ ^ 0 * L ^ 5 / 120) + g2 * s ^ 2 * (τ ^ 1 * L ^ 6 / 720) + g3 * s ^ 3 * (τ ^ 0 * L ^ 4 / 24) + g3 * s ^ 3 * (τ ^ 1 * L ^ 5 / 120) + g3 * s ^ 3 * (τ ^ 2 * L ^ 6 / 720) + g4 * s ^ 4 * (τ ^ 0 * L ^ 3 / 6) + g4 * s ^ 4 * (τ ^ 1 * L ^ 4 / 24) + g4 * s ^ 4 * (τ ^ 2 * L ^ 5 / 120) + g4 * s ^ 4 * (τ ^ 3 * L ^ 6 / 720) + g5 * s ^ 5 * (τ ^ 0 * L ^ 2 / 2) + g5 * s ^ 5 * (τ ^ 1 * L ^ 3 / 6) + g5 * s ^ 5 * (τ ^ 2 * L ^ 4 / 24)) (g5 * s ^ 5 * (τ ^ 3 * L ^ 5 / 120))
  have c14 := abs_add_le (g1 * s ^ 1 * (τ ^ 0 * L ^ 6 / 720) + g2 * s ^ 2 * (τ ^ 0 * L ^ 5 / 120) + g2 * s ^ 2 * (τ ^ 1 * L ^ 6 / 720) + g3 * s ^ 3 * (τ ^ 0 * L ^ 4 / 24) + g3 * s ^ 3 * (τ ^ 1 * L ^ 5 / 120) + g3 * s ^ 3 * (τ ^ 2 * L ^ 6 / 720) + g4 * s ^ 4 * (τ ^ 0 * L ^ 3 / 6) + g4 * s ^ 4 * (τ ^ 1 * L ^ 4 / 24) + g4 * s ^ 4 * (τ ^ 2 * L ^ 5 / 120) + g4 * s ^ 4 * (τ ^ 3 * L ^ 6 / 720) + g5 * s ^ 5 * (τ ^ 0 * L ^ 2 / 2) + g5 * s ^ 5 * (τ ^ 1 * L ^ 3 / 6) + g5 * s ^ 5 * (τ ^ 2 * L ^ 4 / 24) + g5 * s ^ 5 * (τ ^ 3 * L ^ 5 / 120)) (g5 * s ^ 5 * (τ ^ 4 * L ^ 6 / 720))
  have hB8 : ∀ m : ℕ, m ≤ 6 → s ^ m * Real.exp (-s) ≤ B8 s := by
    intro m hm
    unfold B8
    have h0 : 0 ≤ s ^ 0 * Real.exp (-s) := by positivity
    have h1 : 0 ≤ s ^ 1 * Real.exp (-s) := by positivity
    have h2 : 0 ≤ s ^ 2 * Real.exp (-s) := by positivity
    have h3 : 0 ≤ s ^ 3 * Real.exp (-s) := by positivity
    have h4 : 0 ≤ s ^ 4 * Real.exp (-s) := by positivity
    have h5 : 0 ≤ s ^ 5 * Real.exp (-s) := by positivity
    have h6 : 0 ≤ s ^ 6 * Real.exp (-s) := by positivity
    have w : (1 + s + s ^ 2 + s ^ 3 + s ^ 4 + s ^ 5 + s ^ 6) * Real.exp (-s)
      = s ^ 0 * Real.exp (-s) + s ^ 1 * Real.exp (-s) + s ^ 2 * Real.exp (-s) + s ^ 3 * Real.exp (-s)
        + s ^ 4 * Real.exp (-s) + s ^ 5 * Real.exp (-s) + s ^ 6 * Real.exp (-s) := by ring
    rw [w]
    interval_cases m <;> linarith
  have hB0 : 0 ≤ B8 s := B8_nonneg hs
  have hQb : |Q| ≤ 5 * (|g1| + |g2| + |g3| + |g4| + |g5|) * B8 s := by
    have k1 := mul_le_mul_of_nonneg_left (hB8 1 (by norm_num)) (abs_nonneg g1)
    have k2 := mul_le_mul_of_nonneg_left (hB8 2 (by norm_num)) (abs_nonneg g2)
    have k3 := mul_le_mul_of_nonneg_left (hB8 3 (by norm_num)) (abs_nonneg g3)
    have k4 := mul_le_mul_of_nonneg_left (hB8 4 (by norm_num)) (abs_nonneg g4)
    have k5 := mul_le_mul_of_nonneg_left (hB8 5 (by norm_num)) (abs_nonneg g5)
    have n1 := mul_nonneg (abs_nonneg g1) hB0
    have n2 := mul_nonneg (abs_nonneg g2) hB0
    have n3 := mul_nonneg (abs_nonneg g3) hB0
    have n4 := mul_nonneg (abs_nonneg g4) hB0
    have n5 := mul_nonneg (abs_nonneg g5) hB0
    rw [hQ]
    linarith
  -- term 2
  have b2 : |τ * a1 * R| ≤ τ ^ 8 * (2 * CR * s ^ 6 * Real.exp (-s)) := by
    have ha1b : |a1| ≤ 2 * τ * L := by
      have hu0 : 0 ≤ τ * L := mul_nonneg hτ0.le hl0
      have hu1 : τ * L ≤ 1 := by
        have := mul_le_mul hτ1.le hl1 hl0 zero_le_one; linarith
      have p2 : (τ * L) ^ 2 ≤ τ * L := pow_le_of_le_one hu0 hu1 (by norm_num)
      have p3 : (τ * L) ^ 3 ≤ τ * L := pow_le_of_le_one hu0 hu1 (by norm_num)
      have p4 : (τ * L) ^ 4 ≤ τ * L := pow_le_of_le_one hu0 hu1 (by norm_num)
      have p5 : (τ * L) ^ 5 ≤ τ * L := pow_le_of_le_one hu0 hu1 (by norm_num)
      have p6 : (τ * L) ^ 6 ≤ τ * L := pow_le_of_le_one hu0 hu1 (by norm_num)
      have q2 : 0 ≤ (τ * L) ^ 2 := sq_nonneg _
      have q3 : 0 ≤ (τ * L) ^ 3 := pow_nonneg hu0 3
      have q4 : 0 ≤ (τ * L) ^ 4 := pow_nonneg hu0 4
      have q5 : 0 ≤ (τ * L) ^ 5 := pow_nonneg hu0 5
      have q6 : 0 ≤ (τ * L) ^ 6 := pow_nonneg hu0 6
      have : 0 ≤ a1 := by rw [ha1]; linarith
      rw [abs_of_nonneg this, ha1]
      have : 2 * τ * L = 2 * (τ * L) := by ring
      linarith
    have hRb : |R| ≤ CR * (τ * s) ^ 6 := hR
    rw [abs_mul, abs_mul, abs_of_pos hτ0]
    calc τ * |a1| * |R| ≤ τ * (2 * τ * L) * (CR * (τ * s) ^ 6) :=
          mul_le_mul (mul_le_mul_of_nonneg_left ha1b hτ') hRb (abs_nonneg _) (by positivity)
      _ = τ ^ 8 * (2 * CR * s ^ 6 * L) := by ring
      _ ≤ τ ^ 8 * (2 * CR * s ^ 6 * Real.exp (-s)) :=
          mul_le_mul_of_nonneg_left (mul_le_mul_of_nonneg_left hle (by positivity)) (by positivity)
  -- term 3
  have hL7 : L ^ 7 ≤ Real.exp (-s) := (pow_le_of_le_one hl0 hl1 (by norm_num)).trans hle
  have b3 : |τ * r1 * Gv| ≤ τ ^ 8 * (M * Real.exp (-s)) := by
    rw [abs_mul, abs_mul, abs_of_pos hτ0]
    calc τ * |r1| * |Gv| ≤ τ * (τ ^ 7 * L ^ 7) * M :=
          mul_le_mul (mul_le_mul_of_nonneg_left hr1b hτ0.le) hGv (abs_nonneg _) (by positivity)
      _ = τ ^ 8 * (L ^ 7 * M) := by ring
      _ ≤ τ ^ 8 * (Real.exp (-s) * M) :=
          mul_le_mul_of_nonneg_left (mul_le_mul_of_nonneg_right hL7 hM) (by positivity)
      _ = _ := by ring
  have t1 := abs_add_le (τ ^ 8 * Q + τ * a1 * R) (τ * r1 * Gv)
  have t2 := abs_add_le (τ ^ 8 * Q) (τ * a1 * R)
  have b1 : |τ ^ 8 * Q| ≤ τ ^ 8 * (5 * (|g1| + |g2| + |g3| + |g4| + |g5|) * B8 s) := by
    rw [abs_mul, abs_of_pos (by positivity)]
    exact mul_le_mul_of_nonneg_left hQb (by positivity)
  have hin : 5 * (|g1| + |g2| + |g3| + |g4| + |g5|) * B8 s + 2 * CR * s ^ 6 * Real.exp (-s)
      + M * Real.exp (-s) ≤ (5 * (|g1| + |g2| + |g3| + |g4| + |g5|) + 2 * CR + M) * B8 s := by
    have k6 := mul_le_mul_of_nonneg_left (hB8 6 (by norm_num)) (by positivity : (0:ℝ) ≤ 2 * CR)
    have k0 := mul_le_mul_of_nonneg_left (hB8 0 (by norm_num)) hM
    simp only [pow_zero, one_mul] at k0
    nlinarith
  have hfin := mul_le_mul_of_nonneg_left hin (by positivity : (0:ℝ) ≤ τ ^ 8)
  have e : τ ^ 8 * (5 * (|g1| + |g2| + |g3| + |g4| + |g5|) * B8 s + 2 * CR * s ^ 6 * Real.exp (-s)
      + M * Real.exp (-s)) = τ ^ 8 * (5 * (|g1| + |g2| + |g3| + |g4| + |g5|) * B8 s)
      + τ ^ 8 * (2 * CR * s ^ 6 * Real.exp (-s)) + τ ^ 8 * (M * Real.exp (-s)) := by ring
  linarith


lemma int21 {f1 f2 f3 f4 f5 f6 f7 f8 f9 f10 f11 f12 f13 f14 f15 f16 f17 f18 f19 f20 f21 : ℝ → ℝ}
    (c1 c2 c3 c4 c5 c6 c7 c8 c9 c10 c11 c12 c13 c14 c15 c16 c17 c18 c19 c20 c21 : ℝ)
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
    (h21 : IntegrableOn f21 (Set.Ioi 0)) :
    IntegrableOn (fun s => c1 * f1 s + c2 * f2 s + c3 * f3 s + c4 * f4 s + c5 * f5 s + c6 * f6 s + c7 * f7 s + c8 * f8 s + c9 * f9 s + c10 * f10 s + c11 * f11 s + c12 * f12 s + c13 * f13 s + c14 * f14 s + c15 * f15 s + c16 * f16 s + c17 * f17 s + c18 * f18 s + c19 * f19 s + c20 * f20 s + c21 * f21 s) (Set.Ioi 0) ∧
    ∫ s in Set.Ioi (0 : ℝ), (c1 * f1 s + c2 * f2 s + c3 * f3 s + c4 * f4 s + c5 * f5 s + c6 * f6 s + c7 * f7 s + c8 * f8 s + c9 * f9 s + c10 * f10 s + c11 * f11 s + c12 * f12 s + c13 * f13 s + c14 * f14 s + c15 * f15 s + c16 * f16 s + c17 * f17 s + c18 * f18 s + c19 * f19 s + c20 * f20 s + c21 * f21 s)
      = c1 * (∫ s in Set.Ioi (0 : ℝ), f1 s) + c2 * (∫ s in Set.Ioi (0 : ℝ), f2 s) + c3 * (∫ s in Set.Ioi (0 : ℝ), f3 s) + c4 * (∫ s in Set.Ioi (0 : ℝ), f4 s) + c5 * (∫ s in Set.Ioi (0 : ℝ), f5 s) + c6 * (∫ s in Set.Ioi (0 : ℝ), f6 s) + c7 * (∫ s in Set.Ioi (0 : ℝ), f7 s) + c8 * (∫ s in Set.Ioi (0 : ℝ), f8 s) + c9 * (∫ s in Set.Ioi (0 : ℝ), f9 s) + c10 * (∫ s in Set.Ioi (0 : ℝ), f10 s) + c11 * (∫ s in Set.Ioi (0 : ℝ), f11 s) + c12 * (∫ s in Set.Ioi (0 : ℝ), f12 s) + c13 * (∫ s in Set.Ioi (0 : ℝ), f13 s) + c14 * (∫ s in Set.Ioi (0 : ℝ), f14 s) + c15 * (∫ s in Set.Ioi (0 : ℝ), f15 s) + c16 * (∫ s in Set.Ioi (0 : ℝ), f16 s) + c17 * (∫ s in Set.Ioi (0 : ℝ), f17 s) + c18 * (∫ s in Set.Ioi (0 : ℝ), f18 s) + c19 * (∫ s in Set.Ioi (0 : ℝ), f19 s) + c20 * (∫ s in Set.Ioi (0 : ℝ), f20 s) + c21 * (∫ s in Set.Ioi (0 : ℝ), f21 s) := by
  obtain ⟨b15, e15⟩ := V02.int15 c1 c2 c3 c4 c5 c6 c7 c8 c9 c10 c11 c12 c13 c14 c15 h1 h2 h3 h4 h5 h6 h7 h8 h9 h10 h11 h12 h13 h14 h15

  have a16 : Integrable (fun s => c16 * f16 s) (volume.restrict (Set.Ioi 0)) := h16.const_mul c16
  have a17 : Integrable (fun s => c17 * f17 s) (volume.restrict (Set.Ioi 0)) := h17.const_mul c17
  have a18 : Integrable (fun s => c18 * f18 s) (volume.restrict (Set.Ioi 0)) := h18.const_mul c18
  have a19 : Integrable (fun s => c19 * f19 s) (volume.restrict (Set.Ioi 0)) := h19.const_mul c19
  have a20 : Integrable (fun s => c20 * f20 s) (volume.restrict (Set.Ioi 0)) := h20.const_mul c20
  have a21 : Integrable (fun s => c21 * f21 s) (volume.restrict (Set.Ioi 0)) := h21.const_mul c21
  have b16 : Integrable (fun s => c1 * f1 s + c2 * f2 s + c3 * f3 s + c4 * f4 s + c5 * f5 s + c6 * f6 s + c7 * f7 s + c8 * f8 s + c9 * f9 s + c10 * f10 s + c11 * f11 s + c12 * f12 s + c13 * f13 s + c14 * f14 s + c15 * f15 s + c16 * f16 s) (volume.restrict (Set.Ioi 0)) := b15.add a16
  have b17 : Integrable (fun s => c1 * f1 s + c2 * f2 s + c3 * f3 s + c4 * f4 s + c5 * f5 s + c6 * f6 s + c7 * f7 s + c8 * f8 s + c9 * f9 s + c10 * f10 s + c11 * f11 s + c12 * f12 s + c13 * f13 s + c14 * f14 s + c15 * f15 s + c16 * f16 s + c17 * f17 s) (volume.restrict (Set.Ioi 0)) := b16.add a17
  have b18 : Integrable (fun s => c1 * f1 s + c2 * f2 s + c3 * f3 s + c4 * f4 s + c5 * f5 s + c6 * f6 s + c7 * f7 s + c8 * f8 s + c9 * f9 s + c10 * f10 s + c11 * f11 s + c12 * f12 s + c13 * f13 s + c14 * f14 s + c15 * f15 s + c16 * f16 s + c17 * f17 s + c18 * f18 s) (volume.restrict (Set.Ioi 0)) := b17.add a18
  have b19 : Integrable (fun s => c1 * f1 s + c2 * f2 s + c3 * f3 s + c4 * f4 s + c5 * f5 s + c6 * f6 s + c7 * f7 s + c8 * f8 s + c9 * f9 s + c10 * f10 s + c11 * f11 s + c12 * f12 s + c13 * f13 s + c14 * f14 s + c15 * f15 s + c16 * f16 s + c17 * f17 s + c18 * f18 s + c19 * f19 s) (volume.restrict (Set.Ioi 0)) := b18.add a19
  have b20 : Integrable (fun s => c1 * f1 s + c2 * f2 s + c3 * f3 s + c4 * f4 s + c5 * f5 s + c6 * f6 s + c7 * f7 s + c8 * f8 s + c9 * f9 s + c10 * f10 s + c11 * f11 s + c12 * f12 s + c13 * f13 s + c14 * f14 s + c15 * f15 s + c16 * f16 s + c17 * f17 s + c18 * f18 s + c19 * f19 s + c20 * f20 s) (volume.restrict (Set.Ioi 0)) := b19.add a20
  have b21 : Integrable (fun s => c1 * f1 s + c2 * f2 s + c3 * f3 s + c4 * f4 s + c5 * f5 s + c6 * f6 s + c7 * f7 s + c8 * f8 s + c9 * f9 s + c10 * f10 s + c11 * f11 s + c12 * f12 s + c13 * f13 s + c14 * f14 s + c15 * f15 s + c16 * f16 s + c17 * f17 s + c18 * f18 s + c19 * f19 s + c20 * f20 s + c21 * f21 s) (volume.restrict (Set.Ioi 0)) := b20.add a21
  refine ⟨b21, ?_⟩
  rw [integral_add b20 a21, integral_add b19 a20, integral_add b18 a19, integral_add b17 a18, integral_add b16 a17, integral_add b15 a16, e15, integral_const_mul, integral_const_mul, integral_const_mul, integral_const_mul, integral_const_mul, integral_const_mul]

lemma watson7 (Gf : ℝ → ℝ) (hG : ContDiff ℝ (⊤ : ℕ∞) Gf) (hb : ∃ M : ℝ, ∀ x : ℝ, 0 ≤ x → |Gf x| ≤ M) :
    ∃ C : ℝ, ∀ τ : ℝ, 0 < τ → τ < 1 →
      IntegrableOn (fun s : ℝ => ((1 + Real.exp (-s)) ^ τ - 1) * Gf (τ * s)) (Ioi 0) ∧
      |τ * (∫ s in Ioi (0 : ℝ), ((1 + Real.exp (-s)) ^ τ - 1) * Gf (τ * s)) -
          ∑ n ∈ Finset.Icc 2 7, τ ^ n * ∑ j ∈ Finset.Icc 1 (n - 1),
            iteratedDeriv (n - 1 - j) Gf 0 / ((n - 1 - j).factorial : ℝ) /
              (j.factorial : ℝ) * Y04.I (n - 1 - j) j| ≤ C * τ ^ 8 := by
  obtain ⟨M, hM⟩ := hb
  have hM0 : 0 ≤ M := (abs_nonneg _).trans (hM 0 le_rfl)
  obtain ⟨CR, hCR, hT⟩ := taylor5 hG hM
  set g0 := iteratedDeriv 0 Gf 0 / (Nat.factorial 0 : ℝ) with hg0
  set g1 := iteratedDeriv 1 Gf 0 / (Nat.factorial 1 : ℝ) with hg1
  set g2 := iteratedDeriv 2 Gf 0 / (Nat.factorial 2 : ℝ) with hg2
  set g3 := iteratedDeriv 3 Gf 0 / (Nat.factorial 3 : ℝ) with hg3
  set g4 := iteratedDeriv 4 Gf 0 / (Nat.factorial 4 : ℝ) with hg4
  set g5 := iteratedDeriv 5 Gf 0 / (Nat.factorial 5 : ℝ) with hg5
  refine ⟨(5 * (|g1| + |g2| + |g3| + |g4| + |g5|) + 2 * CR + M) * ∫ s in Set.Ioi (0 : ℝ), B8 s, fun τ hτ0 hτ1 => ?_⟩
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
  obtain ⟨hmi, hmv⟩ := int21 (τ ^ 2 * g0) (τ ^ 3 * g1) (τ ^ 3 * (g0 / 2)) (τ ^ 4 * g2) (τ ^ 4 * (g1 / 2)) (τ ^ 4 * (g0 / 6)) (τ ^ 5 * g3) (τ ^ 5 * (g2 / 2)) (τ ^ 5 * (g1 / 6)) (τ ^ 5 * (g0 / 24)) (τ ^ 6 * g4) (τ ^ 6 * (g3 / 2)) (τ ^ 6 * (g2 / 6)) (τ ^ 6 * (g1 / 24)) (τ ^ 6 * (g0 / 120)) (τ ^ 7 * g5) (τ ^ 7 * (g4 / 2)) (τ ^ 7 * (g3 / 6)) (τ ^ 7 * (g2 / 24)) (τ ^ 7 * (g1 / 120)) (τ ^ 7 * (g0 / 720))
    (Y04.I_int 0 1 le_rfl) (Y04.I_int 1 1 le_rfl) (Y04.I_int 0 2 (by norm_num)) (Y04.I_int 2 1 le_rfl) (Y04.I_int 1 2 (by norm_num)) (Y04.I_int 0 3 (by norm_num)) (Y04.I_int 3 1 le_rfl) (Y04.I_int 2 2 (by norm_num)) (Y04.I_int 1 3 (by norm_num)) (Y04.I_int 0 4 (by norm_num)) (Y04.I_int 4 1 le_rfl) (Y04.I_int 3 2 (by norm_num)) (Y04.I_int 2 3 (by norm_num)) (Y04.I_int 1 4 (by norm_num)) (Y04.I_int 0 5 (by norm_num)) (Y04.I_int 5 1 le_rfl) (Y04.I_int 4 2 (by norm_num)) (Y04.I_int 3 3 (by norm_num)) (Y04.I_int 2 4 (by norm_num)) (Y04.I_int 1 5 (by norm_num)) (Y04.I_int 0 6 (by norm_num))
  have hsum : ∑ n ∈ Finset.Icc 2 7, τ ^ n * ∑ j ∈ Finset.Icc 1 (n - 1),
            iteratedDeriv (n - 1 - j) Gf 0 / ((n - 1 - j).factorial : ℝ) /
              (j.factorial : ℝ) * Y04.I (n - 1 - j) j
      = τ ^ 2 * g0 * Y04.I 0 1 + τ ^ 3 * g1 * Y04.I 1 1 + τ ^ 3 * (g0 / 2) * Y04.I 0 2 + τ ^ 4 * g2 * Y04.I 2 1 + τ ^ 4 * (g1 / 2) * Y04.I 1 2 + τ ^ 4 * (g0 / 6) * Y04.I 0 3 + τ ^ 5 * g3 * Y04.I 3 1 + τ ^ 5 * (g2 / 2) * Y04.I 2 2 + τ ^ 5 * (g1 / 6) * Y04.I 1 3 + τ ^ 5 * (g0 / 24) * Y04.I 0 4 + τ ^ 6 * g4 * Y04.I 4 1 + τ ^ 6 * (g3 / 2) * Y04.I 3 2 + τ ^ 6 * (g2 / 6) * Y04.I 2 3 + τ ^ 6 * (g1 / 24) * Y04.I 1 4 + τ ^ 6 * (g0 / 120) * Y04.I 0 5 + τ ^ 7 * g5 * Y04.I 5 1 + τ ^ 7 * (g4 / 2) * Y04.I 4 2 + τ ^ 7 * (g3 / 6) * Y04.I 3 3 + τ ^ 7 * (g2 / 24) * Y04.I 2 4 + τ ^ 7 * (g1 / 120) * Y04.I 1 5 + τ ^ 7 * (g0 / 720) * Y04.I 0 6 := by
    rw [show Finset.Icc 2 7 = {2, 3, 4, 5, 6, 7} by decide]
    rw [Finset.sum_insert (by decide), Finset.sum_insert (by decide), Finset.sum_insert (by decide),
      Finset.sum_insert (by decide), Finset.sum_insert (by decide), Finset.sum_singleton]
    rw [show Finset.Icc 1 (2 - 1) = {1} by decide, show Finset.Icc 1 (3 - 1) = {1, 2} by decide,
      show Finset.Icc 1 (4 - 1) = {1, 2, 3} by decide, show Finset.Icc 1 (5 - 1) = {1, 2, 3, 4} by decide,
      show Finset.Icc 1 (6 - 1) = {1, 2, 3, 4, 5} by decide,
      show Finset.Icc 1 (7 - 1) = {1, 2, 3, 4, 5, 6} by decide]
    simp [Finset.sum_insert, hg0, hg1, hg2, hg3, hg4, hg5, Nat.factorial]
    ring
  rw [hsum]
  unfold Y04.I
  rw [← hmv, ← integral_const_mul, ← integral_sub (hint.const_mul τ) hmi]
  have hB := (B8_int.const_mul (5 * (|g1| + |g2| + |g3| + |g4| + |g5|) + 2 * CR + M)).const_mul (τ ^ 8)
  have := norm_integral_le_of_norm_le hB
    ((ae_restrict_iff' measurableSet_Ioi).2 (Filter.Eventually.of_forall
      (fun s (hs : 0 < s) => by
        rw [Real.norm_eq_abs]
        have hx : 0 ≤ τ * s := by positivity
        have := watson_pt7 (g0 := g0) (g1 := g1) (g2 := g2) (g3 := g3) (g4 := g4) (g5 := g5) hCR hM0 hτ0 hτ1 hs.le
          (hM (τ * s) hx)
          (by have := hT (τ * s) hx; simpa [hg0, hg1, hg2, hg3, hg4, hg5] using this)
        convert this using 2)))
  rw [Real.norm_eq_abs, integral_const_mul, integral_const_mul] at this
  calc _ ≤ _ := this
    _ = _ := by ring



lemma em6_bound {q : ℝ} (hq : 0 ≤ q) :
    |Real.exp (-q) - (1 - q + q ^ 2 / 2 - q ^ 3 / 6 + q ^ 4 / 24 - q ^ 5 / 120)| ≤ 4 * q ^ 6 := by
  rcases le_or_gt q 1 with h | h
  · have h1 : |-q| ≤ 1 := by rw [abs_neg, abs_of_nonneg hq]; exact h
    have he := Real.exp_bound h1 (n := 6) (by norm_num)
    simp only [Finset.sum_range_succ, Finset.sum_range_zero, Nat.factorial] at he
    norm_num at he
    rw [abs_of_nonneg hq] at he
    have e : Real.exp (-q) - (1 - q + q ^ 2 / 2 - q ^ 3 / 6 + q ^ 4 / 24 - q ^ 5 / 120)
        = Real.exp (-q) - (1 + -q + q ^ 2 / 2 + (-q) ^ 3 / 6 + q ^ 4 / 24 + (-q) ^ 5 / 120) := by ring
    rw [e]
    have : 0 ≤ q ^ 6 := by positivity
    nlinarith
  · have he : Real.exp (-q) ≤ 1 := Real.exp_le_one_iff.2 (by linarith)
    have he0 := Real.exp_pos (-q)
    have h2 : q ≤ q ^ 2 := by nlinarith
    have h3 : q ^ 2 ≤ q ^ 3 := by nlinarith
    have h4 : q ^ 3 ≤ q ^ 4 := by nlinarith
    have h5 : q ^ 4 ≤ q ^ 5 := by nlinarith
    have h6 : q ^ 5 ≤ q ^ 6 := by nlinarith
    rw [abs_le]
    constructor <;> nlinarith

lemma ep7_bound {p : ℝ} (hp0 : 0 ≤ p) (hp1 : p ≤ 1) :
    |Real.exp (-p) - 1 - (-p + p ^ 2 / 2 - p ^ 3 / 6 + p ^ 4 / 24 - p ^ 5 / 120 + p ^ 6 / 720)| ≤ p ^ 7 := by
  have h1 : |-p| ≤ 1 := by rw [abs_neg, abs_of_nonneg hp0]; exact hp1
  have he := Real.exp_bound h1 (n := 7) (by norm_num)
  simp only [Finset.sum_range_succ, Finset.sum_range_zero, Nat.factorial] at he
  norm_num at he
  rw [abs_of_nonneg hp0] at he
  have e : Real.exp (-p) - 1 - (-p + p ^ 2 / 2 - p ^ 3 / 6 + p ^ 4 / 24 - p ^ 5 / 120 + p ^ 6 / 720)
      = Real.exp (-p) - (1 + -p + p ^ 2 / 2 + (-p) ^ 3 / 6 + p ^ 4 / 24 + (-p) ^ 5 / 120
        + p ^ 6 / 720) := by ring
  rw [e]
  have : 0 ≤ p ^ 7 := by positivity
  nlinarith

lemma mono_bd7 {u v t L s : ℝ} (hu0 : 0 ≤ u) (hut : u ≤ t) (hv0 : 0 ≤ v) (hvt : v ≤ t) (ht1 : t ≤ 1)
    (hL0 : 0 ≤ L) (hL1 : L ≤ 1) (hLe : L ≤ Real.exp (-s)) (hs : 0 ≤ s) {i j : ℕ} (hi : 1 ≤ i)
    (hij : 7 ≤ i + j) :
    (u * L) ^ i * (v * s) ^ j ≤ t ^ 7 * (s ^ j * Real.exp (-s)) := by
  have ht0 : 0 ≤ t := hu0.trans hut
  rw [mul_pow, mul_pow]
  have a1 : u ^ i ≤ t ^ i := pow_le_pow_left₀ hu0 hut i
  have a2 : v ^ j ≤ t ^ j := pow_le_pow_left₀ hv0 hvt j
  have a3 : t ^ i * t ^ j ≤ t ^ 7 := by
    rw [← pow_add]; exact pow_le_pow_of_le_one ht0 ht1 hij
  have a4 : L ^ i ≤ Real.exp (-s) := (pow_le_of_le_one hL0 hL1 (by omega)).trans hLe
  have b1 : u ^ i * v ^ j ≤ t ^ 7 :=
    (mul_le_mul a1 a2 (by positivity) (by positivity)).trans a3
  calc u ^ i * L ^ i * (v ^ j * s ^ j) = (u ^ i * v ^ j) * L ^ i * s ^ j := by ring
    _ ≤ t ^ 7 * Real.exp (-s) * s ^ j := by
        apply mul_le_mul_of_nonneg_right _ (by positivity)
        exact mul_le_mul b1 a4 (by positivity) (by positivity)
    _ = _ := by ring

lemma sj_le_B8 {s : ℝ} (hs : 0 ≤ s) {j : ℕ} (hj : j ≤ 6) : s ^ j * Real.exp (-s) ≤ B8 s := by
  unfold B8
  have h0 : 0 ≤ s ^ 0 * Real.exp (-s) := by positivity
  have h1 : 0 ≤ s ^ 1 * Real.exp (-s) := by positivity
  have h2 : 0 ≤ s ^ 2 * Real.exp (-s) := by positivity
  have h3 : 0 ≤ s ^ 3 * Real.exp (-s) := by positivity
  have h4 : 0 ≤ s ^ 4 * Real.exp (-s) := by positivity
  have h5 : 0 ≤ s ^ 5 * Real.exp (-s) := by positivity
  have h6 : 0 ≤ s ^ 6 * Real.exp (-s) := by positivity
  have w : (1 + s + s ^ 2 + s ^ 3 + s ^ 4 + s ^ 5 + s ^ 6) * Real.exp (-s)
    = s ^ 0 * Real.exp (-s) + s ^ 1 * Real.exp (-s) + s ^ 2 * Real.exp (-s) + s ^ 3 * Real.exp (-s)
      + s ^ 4 * Real.exp (-s) + s ^ 5 * Real.exp (-s) + s ^ 6 * Real.exp (-s) := by ring
  rw [w]
  interval_cases j <;> linarith

lemma hmono7 {u v t L s p q c : ℝ} (hu0 : 0 ≤ u) (hut : u ≤ t) (hv0 : 0 ≤ v) (hvt : v ≤ t) (ht1 : t ≤ 1)
    (hL0 : 0 ≤ L) (hL1 : L ≤ 1) (hLe : L ≤ Real.exp (-s)) (hs : 0 ≤ s) (hp : p = u * L) (hq : q = v * s)
    (hc : |c| ≤ 1) {i j : ℕ} (hi : 1 ≤ i) (hij : 7 ≤ i + j) (hj : j ≤ 6) :
    |c * (p ^ i * q ^ j)| ≤ t ^ 7 * B8 s := by
  subst hp hq
  have h0 : 0 ≤ (u * L) ^ i * (v * s) ^ j := by positivity
  have h1 := mono_bd7 hu0 hut hv0 hvt ht1 hL0 hL1 hLe hs hi hij
  have h2 := sj_le_B8 hs hj
  have ht0 : 0 ≤ t := hu0.trans hut
  rw [abs_mul, abs_of_nonneg h0]
  calc |c| * ((u * L) ^ i * (v * s) ^ j) ≤ 1 * (t ^ 7 * (s ^ j * Real.exp (-s))) :=
        mul_le_mul hc h1 h0 zero_le_one
    _ ≤ t ^ 7 * B8 s := by rw [one_mul]; exact mul_le_mul_of_nonneg_left h2 (by positivity)

/-- the degree-`≤ 6` part of `(e^{-uℓ} - 1) e^{-vs}` -/
noncomputable def P6 (u v s : ℝ) : ℝ :=
  (-1) * (u ^ 1 * v ^ 0) * (s ^ 0 * Y04.ell s ^ 1)
    + (1) * (u ^ 1 * v ^ 1) * (s ^ 1 * Y04.ell s ^ 1)
    + (-1 / 2) * (u ^ 1 * v ^ 2) * (s ^ 2 * Y04.ell s ^ 1)
    + (1 / 6) * (u ^ 1 * v ^ 3) * (s ^ 3 * Y04.ell s ^ 1)
    + (-1 / 24) * (u ^ 1 * v ^ 4) * (s ^ 4 * Y04.ell s ^ 1)
    + (1 / 120) * (u ^ 1 * v ^ 5) * (s ^ 5 * Y04.ell s ^ 1)
    + (1 / 2) * (u ^ 2 * v ^ 0) * (s ^ 0 * Y04.ell s ^ 2)
    + (-1 / 2) * (u ^ 2 * v ^ 1) * (s ^ 1 * Y04.ell s ^ 2)
    + (1 / 4) * (u ^ 2 * v ^ 2) * (s ^ 2 * Y04.ell s ^ 2)
    + (-1 / 12) * (u ^ 2 * v ^ 3) * (s ^ 3 * Y04.ell s ^ 2)
    + (1 / 48) * (u ^ 2 * v ^ 4) * (s ^ 4 * Y04.ell s ^ 2)
    + (-1 / 6) * (u ^ 3 * v ^ 0) * (s ^ 0 * Y04.ell s ^ 3)
    + (1 / 6) * (u ^ 3 * v ^ 1) * (s ^ 1 * Y04.ell s ^ 3)
    + (-1 / 12) * (u ^ 3 * v ^ 2) * (s ^ 2 * Y04.ell s ^ 3)
    + (1 / 36) * (u ^ 3 * v ^ 3) * (s ^ 3 * Y04.ell s ^ 3)
    + (1 / 24) * (u ^ 4 * v ^ 0) * (s ^ 0 * Y04.ell s ^ 4)
    + (-1 / 24) * (u ^ 4 * v ^ 1) * (s ^ 1 * Y04.ell s ^ 4)
    + (1 / 48) * (u ^ 4 * v ^ 2) * (s ^ 2 * Y04.ell s ^ 4)
    + (-1 / 120) * (u ^ 5 * v ^ 0) * (s ^ 0 * Y04.ell s ^ 5)
    + (1 / 120) * (u ^ 5 * v ^ 1) * (s ^ 1 * Y04.ell s ^ 5)
    + (1 / 720) * (u ^ 6 * v ^ 0) * (s ^ 0 * Y04.ell s ^ 6)

set_option maxHeartbeats 4000000 in
lemma pt6 {u v t s : ℝ} (hu0 : 0 ≤ u) (hut : u ≤ t) (hv0 : 0 ≤ v) (hvt : v ≤ t) (ht1 : t ≤ 1)
    (hs : 0 ≤ s) :
    |(Real.exp (-u * Y04.ell s) - 1) * Real.exp (-v * s) - P6 u v s| ≤ 24 * t ^ 7 * B8 s := by
  have hl0 := Y04.ell_nonneg s
  have hl1 := Y04.ell_le_one s hs
  have hle := Y04.ell_le s
  have ht0 : 0 ≤ t := hu0.trans hut
  unfold P6
  set L := Y04.ell s
  set p := u * L with hp
  set q := v * s with hq
  have hp0 : 0 ≤ p := by positivity
  have hp1 : p ≤ 1 := by
    have := mul_le_mul (hut.trans ht1) hl1 hl0 zero_le_one; rw [hp]; linarith
  have hq0 : 0 ≤ q := by positivity
  set A := -p + p ^ 2 / 2 - p ^ 3 / 6 + p ^ 4 / 24 - p ^ 5 / 120 + p ^ 6 / 720 with hA
  set r1 := Real.exp (-p) - 1 - A with hr1
  set B := 1 - q + q ^ 2 / 2 - q ^ 3 / 6 + q ^ 4 / 24 - q ^ 5 / 120 with hB
  set r2 := Real.exp (-q) - B with hr2
  have hr1b : |r1| ≤ p ^ 7 := ep7_bound hp0 hp1
  have hr2b : |r2| ≤ 4 * q ^ 6 := em6_bound hq0
  set HT := (-1 / 240) * (p ^ 2 * q ^ 5) + (-1 / 144) * (p ^ 3 * q ^ 4) + (1 / 720) * (p ^ 3 * q ^ 5) + (-1 / 144) * (p ^ 4 * q ^ 3) + (1 / 576) * (p ^ 4 * q ^ 4) + (-1 / 2880) * (p ^ 4 * q ^ 5) + (-1 / 240) * (p ^ 5 * q ^ 2) + (1 / 720) * (p ^ 5 * q ^ 3) + (-1 / 2880) * (p ^ 5 * q ^ 4) + (1 / 14400) * (p ^ 5 * q ^ 5) + (-1 / 720) * (p ^ 6 * q ^ 1) + (1 / 1440) * (p ^ 6 * q ^ 2) + (-1 / 4320) * (p ^ 6 * q ^ 3) + (1 / 17280) * (p ^ 6 * q ^ 4) + (-1 / 86400) * (p ^ 6 * q ^ 5) with hHT
  have e1 : -u * L = -p := by rw [hp]; ring
  have e2 : -v * s = -q := by rw [hq]; ring
  rw [e1, e2]
  have key : (Real.exp (-p) - 1) * Real.exp (-q) - ((-1) * (u ^ 1 * v ^ 0) * (s ^ 0 * L ^ 1)
    + (1) * (u ^ 1 * v ^ 1) * (s ^ 1 * L ^ 1)
    + (-1 / 2) * (u ^ 1 * v ^ 2) * (s ^ 2 * L ^ 1)
    + (1 / 6) * (u ^ 1 * v ^ 3) * (s ^ 3 * L ^ 1)
    + (-1 / 24) * (u ^ 1 * v ^ 4) * (s ^ 4 * L ^ 1)
    + (1 / 120) * (u ^ 1 * v ^ 5) * (s ^ 5 * L ^ 1)
    + (1 / 2) * (u ^ 2 * v ^ 0) * (s ^ 0 * L ^ 2)
    + (-1 / 2) * (u ^ 2 * v ^ 1) * (s ^ 1 * L ^ 2)
    + (1 / 4) * (u ^ 2 * v ^ 2) * (s ^ 2 * L ^ 2)
    + (-1 / 12) * (u ^ 2 * v ^ 3) * (s ^ 3 * L ^ 2)
    + (1 / 48) * (u ^ 2 * v ^ 4) * (s ^ 4 * L ^ 2)
    + (-1 / 6) * (u ^ 3 * v ^ 0) * (s ^ 0 * L ^ 3)
    + (1 / 6) * (u ^ 3 * v ^ 1) * (s ^ 1 * L ^ 3)
    + (-1 / 12) * (u ^ 3 * v ^ 2) * (s ^ 2 * L ^ 3)
    + (1 / 36) * (u ^ 3 * v ^ 3) * (s ^ 3 * L ^ 3)
    + (1 / 24) * (u ^ 4 * v ^ 0) * (s ^ 0 * L ^ 4)
    + (-1 / 24) * (u ^ 4 * v ^ 1) * (s ^ 1 * L ^ 4)
    + (1 / 48) * (u ^ 4 * v ^ 2) * (s ^ 2 * L ^ 4)
    + (-1 / 120) * (u ^ 5 * v ^ 0) * (s ^ 0 * L ^ 5)
    + (1 / 120) * (u ^ 5 * v ^ 1) * (s ^ 1 * L ^ 5)
    + (1 / 720) * (u ^ 6 * v ^ 0) * (s ^ 0 * L ^ 6))
      = HT + A * r2 + r1 * Real.exp (-q) := by
    rw [hHT, hr1, hr2, hA, hB, hp, hq]; ring
  rw [key]
  have hB8 : 0 ≤ B8 s := B8_nonneg hs

  have m0 := hmono7 hu0 hut hv0 hvt ht1 hl0 hl1 hle hs hp hq (c := (-1 / 240)) (by rw [abs_le]; constructor <;> norm_num) (i := 2) (j := 5) (by norm_num) (by norm_num) (by norm_num)
  have m1 := hmono7 hu0 hut hv0 hvt ht1 hl0 hl1 hle hs hp hq (c := (-1 / 144)) (by rw [abs_le]; constructor <;> norm_num) (i := 3) (j := 4) (by norm_num) (by norm_num) (by norm_num)
  have m2 := hmono7 hu0 hut hv0 hvt ht1 hl0 hl1 hle hs hp hq (c := (1 / 720)) (by rw [abs_le]; constructor <;> norm_num) (i := 3) (j := 5) (by norm_num) (by norm_num) (by norm_num)
  have m3 := hmono7 hu0 hut hv0 hvt ht1 hl0 hl1 hle hs hp hq (c := (-1 / 144)) (by rw [abs_le]; constructor <;> norm_num) (i := 4) (j := 3) (by norm_num) (by norm_num) (by norm_num)
  have m4 := hmono7 hu0 hut hv0 hvt ht1 hl0 hl1 hle hs hp hq (c := (1 / 576)) (by rw [abs_le]; constructor <;> norm_num) (i := 4) (j := 4) (by norm_num) (by norm_num) (by norm_num)
  have m5 := hmono7 hu0 hut hv0 hvt ht1 hl0 hl1 hle hs hp hq (c := (-1 / 2880)) (by rw [abs_le]; constructor <;> norm_num) (i := 4) (j := 5) (by norm_num) (by norm_num) (by norm_num)
  have m6 := hmono7 hu0 hut hv0 hvt ht1 hl0 hl1 hle hs hp hq (c := (-1 / 240)) (by rw [abs_le]; constructor <;> norm_num) (i := 5) (j := 2) (by norm_num) (by norm_num) (by norm_num)
  have m7 := hmono7 hu0 hut hv0 hvt ht1 hl0 hl1 hle hs hp hq (c := (1 / 720)) (by rw [abs_le]; constructor <;> norm_num) (i := 5) (j := 3) (by norm_num) (by norm_num) (by norm_num)
  have m8 := hmono7 hu0 hut hv0 hvt ht1 hl0 hl1 hle hs hp hq (c := (-1 / 2880)) (by rw [abs_le]; constructor <;> norm_num) (i := 5) (j := 4) (by norm_num) (by norm_num) (by norm_num)
  have m9 := hmono7 hu0 hut hv0 hvt ht1 hl0 hl1 hle hs hp hq (c := (1 / 14400)) (by rw [abs_le]; constructor <;> norm_num) (i := 5) (j := 5) (by norm_num) (by norm_num) (by norm_num)
  have m10 := hmono7 hu0 hut hv0 hvt ht1 hl0 hl1 hle hs hp hq (c := (-1 / 720)) (by rw [abs_le]; constructor <;> norm_num) (i := 6) (j := 1) (by norm_num) (by norm_num) (by norm_num)
  have m11 := hmono7 hu0 hut hv0 hvt ht1 hl0 hl1 hle hs hp hq (c := (1 / 1440)) (by rw [abs_le]; constructor <;> norm_num) (i := 6) (j := 2) (by norm_num) (by norm_num) (by norm_num)
  have m12 := hmono7 hu0 hut hv0 hvt ht1 hl0 hl1 hle hs hp hq (c := (-1 / 4320)) (by rw [abs_le]; constructor <;> norm_num) (i := 6) (j := 3) (by norm_num) (by norm_num) (by norm_num)
  have m13 := hmono7 hu0 hut hv0 hvt ht1 hl0 hl1 hle hs hp hq (c := (1 / 17280)) (by rw [abs_le]; constructor <;> norm_num) (i := 6) (j := 4) (by norm_num) (by norm_num) (by norm_num)
  have m14 := hmono7 hu0 hut hv0 hvt ht1 hl0 hl1 hle hs hp hq (c := (-1 / 86400)) (by rw [abs_le]; constructor <;> norm_num) (i := 6) (j := 5) (by norm_num) (by norm_num) (by norm_num)
  have c1 := abs_add_le ((-1 / 240) * (p ^ 2 * q ^ 5)) ((-1 / 144) * (p ^ 3 * q ^ 4))
  have c2 := abs_add_le ((-1 / 240) * (p ^ 2 * q ^ 5) + (-1 / 144) * (p ^ 3 * q ^ 4)) ((1 / 720) * (p ^ 3 * q ^ 5))
  have c3 := abs_add_le ((-1 / 240) * (p ^ 2 * q ^ 5) + (-1 / 144) * (p ^ 3 * q ^ 4) + (1 / 720) * (p ^ 3 * q ^ 5)) ((-1 / 144) * (p ^ 4 * q ^ 3))
  have c4 := abs_add_le ((-1 / 240) * (p ^ 2 * q ^ 5) + (-1 / 144) * (p ^ 3 * q ^ 4) + (1 / 720) * (p ^ 3 * q ^ 5) + (-1 / 144) * (p ^ 4 * q ^ 3)) ((1 / 576) * (p ^ 4 * q ^ 4))
  have c5 := abs_add_le ((-1 / 240) * (p ^ 2 * q ^ 5) + (-1 / 144) * (p ^ 3 * q ^ 4) + (1 / 720) * (p ^ 3 * q ^ 5) + (-1 / 144) * (p ^ 4 * q ^ 3) + (1 / 576) * (p ^ 4 * q ^ 4)) ((-1 / 2880) * (p ^ 4 * q ^ 5))
  have c6 := abs_add_le ((-1 / 240) * (p ^ 2 * q ^ 5) + (-1 / 144) * (p ^ 3 * q ^ 4) + (1 / 720) * (p ^ 3 * q ^ 5) + (-1 / 144) * (p ^ 4 * q ^ 3) + (1 / 576) * (p ^ 4 * q ^ 4) + (-1 / 2880) * (p ^ 4 * q ^ 5)) ((-1 / 240) * (p ^ 5 * q ^ 2))
  have c7 := abs_add_le ((-1 / 240) * (p ^ 2 * q ^ 5) + (-1 / 144) * (p ^ 3 * q ^ 4) + (1 / 720) * (p ^ 3 * q ^ 5) + (-1 / 144) * (p ^ 4 * q ^ 3) + (1 / 576) * (p ^ 4 * q ^ 4) + (-1 / 2880) * (p ^ 4 * q ^ 5) + (-1 / 240) * (p ^ 5 * q ^ 2)) ((1 / 720) * (p ^ 5 * q ^ 3))
  have c8 := abs_add_le ((-1 / 240) * (p ^ 2 * q ^ 5) + (-1 / 144) * (p ^ 3 * q ^ 4) + (1 / 720) * (p ^ 3 * q ^ 5) + (-1 / 144) * (p ^ 4 * q ^ 3) + (1 / 576) * (p ^ 4 * q ^ 4) + (-1 / 2880) * (p ^ 4 * q ^ 5) + (-1 / 240) * (p ^ 5 * q ^ 2) + (1 / 720) * (p ^ 5 * q ^ 3)) ((-1 / 2880) * (p ^ 5 * q ^ 4))
  have c9 := abs_add_le ((-1 / 240) * (p ^ 2 * q ^ 5) + (-1 / 144) * (p ^ 3 * q ^ 4) + (1 / 720) * (p ^ 3 * q ^ 5) + (-1 / 144) * (p ^ 4 * q ^ 3) + (1 / 576) * (p ^ 4 * q ^ 4) + (-1 / 2880) * (p ^ 4 * q ^ 5) + (-1 / 240) * (p ^ 5 * q ^ 2) + (1 / 720) * (p ^ 5 * q ^ 3) + (-1 / 2880) * (p ^ 5 * q ^ 4)) ((1 / 14400) * (p ^ 5 * q ^ 5))
  have c10 := abs_add_le ((-1 / 240) * (p ^ 2 * q ^ 5) + (-1 / 144) * (p ^ 3 * q ^ 4) + (1 / 720) * (p ^ 3 * q ^ 5) + (-1 / 144) * (p ^ 4 * q ^ 3) + (1 / 576) * (p ^ 4 * q ^ 4) + (-1 / 2880) * (p ^ 4 * q ^ 5) + (-1 / 240) * (p ^ 5 * q ^ 2) + (1 / 720) * (p ^ 5 * q ^ 3) + (-1 / 2880) * (p ^ 5 * q ^ 4) + (1 / 14400) * (p ^ 5 * q ^ 5)) ((-1 / 720) * (p ^ 6 * q ^ 1))
  have c11 := abs_add_le ((-1 / 240) * (p ^ 2 * q ^ 5) + (-1 / 144) * (p ^ 3 * q ^ 4) + (1 / 720) * (p ^ 3 * q ^ 5) + (-1 / 144) * (p ^ 4 * q ^ 3) + (1 / 576) * (p ^ 4 * q ^ 4) + (-1 / 2880) * (p ^ 4 * q ^ 5) + (-1 / 240) * (p ^ 5 * q ^ 2) + (1 / 720) * (p ^ 5 * q ^ 3) + (-1 / 2880) * (p ^ 5 * q ^ 4) + (1 / 14400) * (p ^ 5 * q ^ 5) + (-1 / 720) * (p ^ 6 * q ^ 1)) ((1 / 1440) * (p ^ 6 * q ^ 2))
  have c12 := abs_add_le ((-1 / 240) * (p ^ 2 * q ^ 5) + (-1 / 144) * (p ^ 3 * q ^ 4) + (1 / 720) * (p ^ 3 * q ^ 5) + (-1 / 144) * (p ^ 4 * q ^ 3) + (1 / 576) * (p ^ 4 * q ^ 4) + (-1 / 2880) * (p ^ 4 * q ^ 5) + (-1 / 240) * (p ^ 5 * q ^ 2) + (1 / 720) * (p ^ 5 * q ^ 3) + (-1 / 2880) * (p ^ 5 * q ^ 4) + (1 / 14400) * (p ^ 5 * q ^ 5) + (-1 / 720) * (p ^ 6 * q ^ 1) + (1 / 1440) * (p ^ 6 * q ^ 2)) ((-1 / 4320) * (p ^ 6 * q ^ 3))
  have c13 := abs_add_le ((-1 / 240) * (p ^ 2 * q ^ 5) + (-1 / 144) * (p ^ 3 * q ^ 4) + (1 / 720) * (p ^ 3 * q ^ 5) + (-1 / 144) * (p ^ 4 * q ^ 3) + (1 / 576) * (p ^ 4 * q ^ 4) + (-1 / 2880) * (p ^ 4 * q ^ 5) + (-1 / 240) * (p ^ 5 * q ^ 2) + (1 / 720) * (p ^ 5 * q ^ 3) + (-1 / 2880) * (p ^ 5 * q ^ 4) + (1 / 14400) * (p ^ 5 * q ^ 5) + (-1 / 720) * (p ^ 6 * q ^ 1) + (1 / 1440) * (p ^ 6 * q ^ 2) + (-1 / 4320) * (p ^ 6 * q ^ 3)) ((1 / 17280) * (p ^ 6 * q ^ 4))
  have c14 := abs_add_le ((-1 / 240) * (p ^ 2 * q ^ 5) + (-1 / 144) * (p ^ 3 * q ^ 4) + (1 / 720) * (p ^ 3 * q ^ 5) + (-1 / 144) * (p ^ 4 * q ^ 3) + (1 / 576) * (p ^ 4 * q ^ 4) + (-1 / 2880) * (p ^ 4 * q ^ 5) + (-1 / 240) * (p ^ 5 * q ^ 2) + (1 / 720) * (p ^ 5 * q ^ 3) + (-1 / 2880) * (p ^ 5 * q ^ 4) + (1 / 14400) * (p ^ 5 * q ^ 5) + (-1 / 720) * (p ^ 6 * q ^ 1) + (1 / 1440) * (p ^ 6 * q ^ 2) + (-1 / 4320) * (p ^ 6 * q ^ 3) + (1 / 17280) * (p ^ 6 * q ^ 4)) ((-1 / 86400) * (p ^ 6 * q ^ 5))
  have bHT : |HT| ≤ 15 * t ^ 7 * B8 s := by rw [hHT]; linarith
  have hAb : |A| ≤ 2 * p := by
    have p2 : p ^ 2 ≤ p := pow_le_of_le_one hp0 hp1 (by norm_num)
    have p3 : p ^ 3 ≤ p := pow_le_of_le_one hp0 hp1 (by norm_num)
    have p4 : p ^ 4 ≤ p := pow_le_of_le_one hp0 hp1 (by norm_num)
    have p5 : p ^ 5 ≤ p := pow_le_of_le_one hp0 hp1 (by norm_num)
    have p6 : p ^ 6 ≤ p := pow_le_of_le_one hp0 hp1 (by norm_num)
    have : 0 ≤ p ^ 2 := by positivity
    have : 0 ≤ p ^ 3 := by positivity
    have : 0 ≤ p ^ 4 := by positivity
    have : 0 ≤ p ^ 5 := by positivity
    have : 0 ≤ p ^ 6 := by positivity
    rw [abs_le, hA]; constructor <;> linarith
  have bAr : |A * r2| ≤ 8 * t ^ 7 * B8 s := by
    rw [abs_mul]
    have h1 : |A| * |r2| ≤ 2 * p * (4 * q ^ 6) :=
      mul_le_mul hAb hr2b (abs_nonneg _) (by positivity)
    have h2 : |(1 : ℝ) * (p ^ 1 * q ^ 6)| ≤ t ^ 7 * B8 s :=
      hmono7 hu0 hut hv0 hvt ht1 hl0 hl1 hle hs hp hq (by norm_num) (by norm_num) (by norm_num) (by norm_num)
    rw [one_mul, pow_one, abs_of_nonneg (by positivity)] at h2
    calc |A| * |r2| ≤ 2 * p * (4 * q ^ 6) := h1
      _ = 8 * (p * q ^ 6) := by ring
      _ ≤ 8 * (t ^ 7 * B8 s) := by linarith
      _ = _ := by ring
  have br1 : |r1 * Real.exp (-q)| ≤ t ^ 7 * B8 s := by
    rw [abs_mul, abs_of_pos (Real.exp_pos _)]
    have h1 : Real.exp (-q) ≤ 1 := Real.exp_le_one_iff.2 (by linarith)
    have h2 : |(1 : ℝ) * (p ^ 7 * q ^ 0)| ≤ t ^ 7 * B8 s :=
      hmono7 hu0 hut hv0 hvt ht1 hl0 hl1 hle hs hp hq (by norm_num) (by norm_num) (by norm_num) (by norm_num)
    rw [one_mul, pow_zero, mul_one, abs_of_nonneg (by positivity)] at h2
    calc |r1| * Real.exp (-q) ≤ p ^ 7 * 1 :=
          mul_le_mul hr1b h1 (Real.exp_pos _).le (by positivity)
      _ ≤ _ := by linarith
  have t1 := abs_add_le (HT + A * r2) (r1 * Real.exp (-q))
  have t2 := abs_add_le HT (A * r2)
  have : 0 ≤ t ^ 7 * B8 s := by positivity
  linarith


/-- the integrated expansion of `(e^{-uℓ} - 1) e^{-vs}` to order `6` -/
lemma int_expand6 {u v t : ℝ} (hu0 : 0 ≤ u) (hut : u ≤ t) (hv0 : 0 ≤ v) (hvt : v ≤ t) (ht1 : t ≤ 1) :
    IntegrableOn (fun s => (Real.exp (-u * Y04.ell s) - 1) * Real.exp (-v * s)) (Set.Ioi 0) ∧
    |(∫ s in Set.Ioi (0 : ℝ), (Real.exp (-u * Y04.ell s) - 1) * Real.exp (-v * s))
      - ((-1) * (u ^ 1 * v ^ 0) * Y04.I 0 1 + (1) * (u ^ 1 * v ^ 1) * Y04.I 1 1 + (-1 / 2) * (u ^ 1 * v ^ 2) * Y04.I 2 1 + (1 / 6) * (u ^ 1 * v ^ 3) * Y04.I 3 1 + (-1 / 24) * (u ^ 1 * v ^ 4) * Y04.I 4 1 + (1 / 120) * (u ^ 1 * v ^ 5) * Y04.I 5 1 + (1 / 2) * (u ^ 2 * v ^ 0) * Y04.I 0 2 + (-1 / 2) * (u ^ 2 * v ^ 1) * Y04.I 1 2 + (1 / 4) * (u ^ 2 * v ^ 2) * Y04.I 2 2 + (-1 / 12) * (u ^ 2 * v ^ 3) * Y04.I 3 2 + (1 / 48) * (u ^ 2 * v ^ 4) * Y04.I 4 2 + (-1 / 6) * (u ^ 3 * v ^ 0) * Y04.I 0 3 + (1 / 6) * (u ^ 3 * v ^ 1) * Y04.I 1 3 + (-1 / 12) * (u ^ 3 * v ^ 2) * Y04.I 2 3 + (1 / 36) * (u ^ 3 * v ^ 3) * Y04.I 3 3 + (1 / 24) * (u ^ 4 * v ^ 0) * Y04.I 0 4 + (-1 / 24) * (u ^ 4 * v ^ 1) * Y04.I 1 4 + (1 / 48) * (u ^ 4 * v ^ 2) * Y04.I 2 4 + (-1 / 120) * (u ^ 5 * v ^ 0) * Y04.I 0 5 + (1 / 120) * (u ^ 5 * v ^ 1) * Y04.I 1 5 + (1 / 720) * (u ^ 6 * v ^ 0) * Y04.I 0 6)|
      ≤ 24 * t ^ 7 * ∫ s in Set.Ioi (0 : ℝ), B8 s := by
  have hint : IntegrableOn (fun s => (Real.exp (-u * Y04.ell s) - 1) * Real.exp (-v * s)) (Set.Ioi 0) :=
    (Z02.int_expand hu0 hut hv0 hvt ht1).1
  refine ⟨hint, ?_⟩
  obtain ⟨hmi, hmv⟩ := int21 ((-1) * (u ^ 1 * v ^ 0)) ((1) * (u ^ 1 * v ^ 1)) ((-1 / 2) * (u ^ 1 * v ^ 2)) ((1 / 6) * (u ^ 1 * v ^ 3)) ((-1 / 24) * (u ^ 1 * v ^ 4)) ((1 / 120) * (u ^ 1 * v ^ 5)) ((1 / 2) * (u ^ 2 * v ^ 0)) ((-1 / 2) * (u ^ 2 * v ^ 1)) ((1 / 4) * (u ^ 2 * v ^ 2)) ((-1 / 12) * (u ^ 2 * v ^ 3)) ((1 / 48) * (u ^ 2 * v ^ 4)) ((-1 / 6) * (u ^ 3 * v ^ 0)) ((1 / 6) * (u ^ 3 * v ^ 1)) ((-1 / 12) * (u ^ 3 * v ^ 2)) ((1 / 36) * (u ^ 3 * v ^ 3)) ((1 / 24) * (u ^ 4 * v ^ 0)) ((-1 / 24) * (u ^ 4 * v ^ 1)) ((1 / 48) * (u ^ 4 * v ^ 2)) ((-1 / 120) * (u ^ 5 * v ^ 0)) ((1 / 120) * (u ^ 5 * v ^ 1)) ((1 / 720) * (u ^ 6 * v ^ 0))
    (Y04.I_int 0 1 le_rfl) (Y04.I_int 1 1 le_rfl) (Y04.I_int 2 1 le_rfl) (Y04.I_int 3 1 le_rfl) (Y04.I_int 4 1 le_rfl) (Y04.I_int 5 1 le_rfl) (Y04.I_int 0 2 (by norm_num)) (Y04.I_int 1 2 (by norm_num)) (Y04.I_int 2 2 (by norm_num)) (Y04.I_int 3 2 (by norm_num)) (Y04.I_int 4 2 (by norm_num)) (Y04.I_int 0 3 (by norm_num)) (Y04.I_int 1 3 (by norm_num)) (Y04.I_int 2 3 (by norm_num)) (Y04.I_int 3 3 (by norm_num)) (Y04.I_int 0 4 (by norm_num)) (Y04.I_int 1 4 (by norm_num)) (Y04.I_int 2 4 (by norm_num)) (Y04.I_int 0 5 (by norm_num)) (Y04.I_int 1 5 (by norm_num)) (Y04.I_int 0 6 (by norm_num))
  unfold Y04.I
  rw [← hmv, ← integral_sub hint hmi]
  have hB := (B8_int.const_mul 24).const_mul (t ^ 7)
  have := norm_integral_le_of_norm_le hB
    ((ae_restrict_iff' measurableSet_Ioi).2 (Filter.Eventually.of_forall
      (fun s (hs : 0 < s) => by
        rw [Real.norm_eq_abs]
        have := pt6 hu0 hut hv0 hvt ht1 hs.le
        unfold P6 at this
        calc _ = _ := rfl
          _ ≤ _ := this
          _ = t ^ 7 * (24 * B8 s) := by ring)))
  rw [Real.norm_eq_abs, integral_const_mul, integral_const_mul] at this
  calc _ ≤ _ := this
    _ = _ := by ring


/-- six coefficients vanish if the sextic is `O(t⁷)` on `(0, d]` -/
lemma sextic_zero {d1 d2 d3 d4 d5 d6 K d : ℝ} (hd : 0 < d) (hd1 : d ≤ 1) (hK : 0 ≤ K)
    (h : ∀ t : ℝ, 0 < t → t ≤ d →
      |d1 * t + d2 * t ^ 2 + d3 * t ^ 3 + d4 * t ^ 4 + d5 * t ^ 5 + d6 * t ^ 6| ≤ K * t ^ 7) :
    d1 = 0 ∧ d2 = 0 ∧ d3 = 0 ∧ d4 = 0 ∧ d5 = 0 ∧ d6 = 0 := by
  have h1 : d1 = 0 := by
    refine Y04.eq_zero_of_le_lin hd (K := |d2| + |d3| + |d4| + |d5| + |d6| + K)
      (by have := abs_nonneg d2; have := abs_nonneg d3; have := abs_nonneg d4; have := abs_nonneg d5
          have := abs_nonneg d6
          linarith)
      (fun t ht0 htd => ?_)
    have ht1 : t ≤ 1 := htd.trans hd1
    have := h t ht0 htd
    set E := d1 * t + d2 * t ^ 2 + d3 * t ^ 3 + d4 * t ^ 4 + d5 * t ^ 5 + d6 * t ^ 6 with hE
    have e : d1 = E / t - d2 * t - d3 * t ^ 2 - d4 * t ^ 3 - d5 * t ^ 4 - d6 * t ^ 5 := by
      rw [hE]; field_simp; ring
    rw [e]
    have ht2 : t ^ 2 ≤ t := by nlinarith
    have ht3 : t ^ 3 ≤ t := by nlinarith [pow_le_one₀ ht0.le ht1 (n := 2)]
    have ht4 : t ^ 4 ≤ t := by nlinarith [pow_le_one₀ ht0.le ht1 (n := 3)]
    have ht5 : t ^ 5 ≤ t := by nlinarith [pow_le_one₀ ht0.le ht1 (n := 4)]
    have ht6 : t ^ 6 ≤ t := by nlinarith [pow_le_one₀ ht0.le ht1 (n := 5)]
    have b1 : |E / t| ≤ K * t ^ 6 := by
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
    have t6 : K * t ^ 6 ≤ K * t := mul_le_mul_of_nonneg_left ht6 hK
    have := abs_sub (E / t - d2 * t - d3 * t ^ 2 - d4 * t ^ 3 - d5 * t ^ 4) (d6 * t ^ 5)
    have := abs_sub (E / t - d2 * t - d3 * t ^ 2 - d4 * t ^ 3) (d5 * t ^ 4)
    have := abs_sub (E / t - d2 * t - d3 * t ^ 2) (d4 * t ^ 3)
    have := abs_sub (E / t - d2 * t) (d3 * t ^ 2)
    have := abs_sub (E / t) (d2 * t)
    have : (|d2| + |d3| + |d4| + |d5| + |d6| + K) * t
        = |d2| * t + |d3| * t + |d4| * t + |d5| * t + |d6| * t + K * t := by ring
    linarith
  subst h1
  obtain ⟨a, b, c, e, f⟩ := V02.quintic_zero hd hd1 hK (d1 := d2) (d2 := d3) (d3 := d4) (d4 := d5) (d5 := d6)
    (fun t ht0 htd => by
    have := h t ht0 htd
    have e : 0 * t + d2 * t ^ 2 + d3 * t ^ 3 + d4 * t ^ 4 + d5 * t ^ 5 + d6 * t ^ 6
        = t * (d2 * t + d3 * t ^ 2 + d4 * t ^ 3 + d5 * t ^ 4 + d6 * t ^ 5) := by ring
    rw [e, abs_mul, abs_of_pos ht0, show K * t ^ 7 = t * (K * t ^ 6) by ring] at this
    exact le_of_mul_le_mul_left this ht0)
  exact ⟨rfl, a, b, c, e, f⟩



lemma series_trunc7 {α β γ x L : ℝ} (hα : |α| ≤ 3 / 2) (hβ : |β| ≤ 3 / 2) (hγ : |γ| ≤ 3 / 2)
    (hx0 : 0 ≤ x) (hx : x ≤ 1 / 3)
    (hL : HasSum (fun k : ℕ => (-1 : ℝ) ^ (k + 2) * Y04.zn (k + 2) / ((k + 2 : ℕ) : ℝ) *
      ((α * x) ^ (k + 2) + (β * x) ^ (k + 2) - (γ * x) ^ (k + 2))) L) :
    |L - (Y04.cf α β γ 0 * x ^ 2 + Y04.cf α β γ 1 * x ^ 3 + Y04.cf α β γ 2 * x ^ 4
      + Y04.cf α β γ 3 * x ^ 5 + Y04.cf α β γ 4 * x ^ 6 + Y04.cf α β γ 5 * x ^ 7)| ≤ 154 * x ^ 8 := by
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
  have hg := (hasSum_nat_add_iff' 6).mpr hL
  have hgeo : HasSum (fun k : ℕ => 3 * r ^ 8 * r ^ k) (3 * r ^ 8 * (1 - r)⁻¹) :=
    (hasSum_geometric_of_lt_one hr0 (by linarith)).mul_left _
  have hb1 : ∀ k : ℕ, Y04.cf α β γ (k + 6) * x ^ (k + 6 + 2) ≤ 3 * r ^ 8 * r ^ k := by
    intro k
    have := (abs_le.mp (hbd (k + 6))).2
    calc _ ≤ 3 * r ^ (k + 6 + 2) := this
      _ = _ := by ring
  have hb2 : ∀ k : ℕ, -(Y04.cf α β γ (k + 6) * x ^ (k + 6 + 2)) ≤ 3 * r ^ 8 * r ^ k := by
    intro k
    have := (abs_le.mp (hbd (k + 6))).1
    calc _ ≤ 3 * r ^ (k + 6 + 2) := by linarith
      _ = _ := by ring
  have u1 := hasSum_le hb1 hg hgeo
  have u2 := hasSum_le hb2 hg.neg hgeo
  have hinv : 3 * r ^ 8 * (1 - r)⁻¹ ≤ 6 * r ^ 8 := by
    rw [mul_inv_le_iff₀ (by linarith)]
    nlinarith [pow_nonneg hr0 8]
  have hr8 : 6 * r ^ 8 ≤ 154 * x ^ 8 := by
    rw [hr, mul_pow]
    nlinarith [pow_nonneg hx0 8]
  simp only [Finset.sum_range_succ, Finset.sum_range_zero, zero_add] at u1 u2
  rw [abs_le]
  constructor <;> norm_num at u1 u2 ⊢ <;> linarith


set_option maxHeartbeats 8000000 in
lemma exp_trunc7 {a0 a1 a2 a3 a4 a5 B K : ℝ} (hB : 0 ≤ B) (_hK : 0 ≤ K)
    (ha0 : |a0| ≤ B) (ha1 : |a1| ≤ B) (ha2 : |a2| ≤ B) (ha3 : |a3| ≤ B) (ha4 : |a4| ≤ B) (ha5 : |a5| ≤ B) :
    ∃ KW : ℝ, 0 ≤ KW ∧ ∀ x L : ℝ, 0 ≤ x → x ≤ 1 → 6 * B * x ^ 2 ≤ 1 → K * x ^ 8 ≤ 1 →
      |L - (a0 * x ^ 2 + a1 * x ^ 3 + a2 * x ^ 4 + a3 * x ^ 5 + a4 * x ^ 6 + a5 * x ^ 7)| ≤ K * x ^ 8 →
      |Real.exp L - (1 + a0 * x ^ 2 + a1 * x ^ 3 + (a2 + a0 ^ 2 / 2) * x ^ 4 + (a3 + a0 * a1) * x ^ 5 + (a4 + a0 * a2 + a1 ^ 2 / 2 + a0 ^ 3 / 6) * x ^ 6 + (a5 + a0 * a3 + a1 * a2 + a0 ^ 2 * a1 / 2) * x ^ 7)| ≤ (6 * K + (6 * B) ^ 4 + KW) * x ^ 8 := by
  have hWc : Continuous (fun x : ℝ => (1 / 2) * (a5 ^ 2 * x ^ 6) + (1 / 6) * (a5 ^ 3 * x ^ 13) + (1) * (a4 * a5 * x ^ 5) + (1 / 2) * (a4 * a5 ^ 2 * x ^ 12) + (1 / 2) * (a4 ^ 2 * x ^ 4) + (1 / 2) * (a4 ^ 2 * a5 * x ^ 11) + (1 / 6) * (a4 ^ 3 * x ^ 10) + (1) * (a3 * a5 * x ^ 4) + (1 / 2) * (a3 * a5 ^ 2 * x ^ 11) + (1) * (a3 * a4 * x ^ 3) + (1) * (a3 * a4 * a5 * x ^ 10) + (1 / 2) * (a3 * a4 ^ 2 * x ^ 9) + (1 / 2) * (a3 ^ 2 * x ^ 2) + (1 / 2) * (a3 ^ 2 * a5 * x ^ 9) + (1 / 2) * (a3 ^ 2 * a4 * x ^ 8) + (1 / 6) * (a3 ^ 3 * x ^ 7) + (1) * (a2 * a5 * x ^ 3) + (1 / 2) * (a2 * a5 ^ 2 * x ^ 10) + (1) * (a2 * a4 * x ^ 2) + (1) * (a2 * a4 * a5 * x ^ 9) + (1 / 2) * (a2 * a4 ^ 2 * x ^ 8) + (1) * (a2 * a3 * x) + (1) * (a2 * a3 * a5 * x ^ 8) + (1) * (a2 * a3 * a4 * x ^ 7) + (1 / 2) * (a2 * a3 ^ 2 * x ^ 6) + (1 / 2) * (a2 ^ 2) + (1 / 2) * (a2 ^ 2 * a5 * x ^ 7) + (1 / 2) * (a2 ^ 2 * a4 * x ^ 6) + (1 / 2) * (a2 ^ 2 * a3 * x ^ 5) + (1 / 6) * (a2 ^ 3 * x ^ 4) + (1) * (a1 * a5 * x ^ 2) + (1 / 2) * (a1 * a5 ^ 2 * x ^ 9) + (1) * (a1 * a4 * x) + (1) * (a1 * a4 * a5 * x ^ 8) + (1 / 2) * (a1 * a4 ^ 2 * x ^ 7) + (1) * (a1 * a3) + (1) * (a1 * a3 * a5 * x ^ 7) + (1) * (a1 * a3 * a4 * x ^ 6) + (1 / 2) * (a1 * a3 ^ 2 * x ^ 5) + (1) * (a1 * a2 * a5 * x ^ 6) + (1) * (a1 * a2 * a4 * x ^ 5) + (1) * (a1 * a2 * a3 * x ^ 4) + (1 / 2) * (a1 * a2 ^ 2 * x ^ 3) + (1 / 2) * (a1 ^ 2 * a5 * x ^ 5) + (1 / 2) * (a1 ^ 2 * a4 * x ^ 4) + (1 / 2) * (a1 ^ 2 * a3 * x ^ 3) + (1 / 2) * (a1 ^ 2 * a2 * x ^ 2) + (1 / 6) * (a1 ^ 3 * x) + (1) * (a0 * a5 * x) + (1 / 2) * (a0 * a5 ^ 2 * x ^ 8) + (1) * (a0 * a4) + (1) * (a0 * a4 * a5 * x ^ 7) + (1 / 2) * (a0 * a4 ^ 2 * x ^ 6) + (1) * (a0 * a3 * a5 * x ^ 6) + (1) * (a0 * a3 * a4 * x ^ 5) + (1 / 2) * (a0 * a3 ^ 2 * x ^ 4) + (1) * (a0 * a2 * a5 * x ^ 5) + (1) * (a0 * a2 * a4 * x ^ 4) + (1) * (a0 * a2 * a3 * x ^ 3) + (1 / 2) * (a0 * a2 ^ 2 * x ^ 2) + (1) * (a0 * a1 * a5 * x ^ 4) + (1) * (a0 * a1 * a4 * x ^ 3) + (1) * (a0 * a1 * a3 * x ^ 2) + (1) * (a0 * a1 * a2 * x) + (1 / 2) * (a0 * a1 ^ 2) + (1 / 2) * (a0 ^ 2 * a5 * x ^ 3) + (1 / 2) * (a0 ^ 2 * a4 * x ^ 2) + (1 / 2) * (a0 ^ 2 * a3 * x) + (1 / 2) * (a0 ^ 2 * a2)) := by fun_prop
  obtain ⟨CW, hCW⟩ := (isCompact_Icc (a := (0 : ℝ)) (b := 1)).exists_bound_of_continuousOn hWc.continuousOn
  refine ⟨max CW 0, le_max_right _ _, fun x L hx0 hx1 hs1 hs2 hL => ?_⟩
  set W := (1 / 2) * (a5 ^ 2 * x ^ 6) + (1 / 6) * (a5 ^ 3 * x ^ 13) + (1) * (a4 * a5 * x ^ 5) + (1 / 2) * (a4 * a5 ^ 2 * x ^ 12) + (1 / 2) * (a4 ^ 2 * x ^ 4) + (1 / 2) * (a4 ^ 2 * a5 * x ^ 11) + (1 / 6) * (a4 ^ 3 * x ^ 10) + (1) * (a3 * a5 * x ^ 4) + (1 / 2) * (a3 * a5 ^ 2 * x ^ 11) + (1) * (a3 * a4 * x ^ 3) + (1) * (a3 * a4 * a5 * x ^ 10) + (1 / 2) * (a3 * a4 ^ 2 * x ^ 9) + (1 / 2) * (a3 ^ 2 * x ^ 2) + (1 / 2) * (a3 ^ 2 * a5 * x ^ 9) + (1 / 2) * (a3 ^ 2 * a4 * x ^ 8) + (1 / 6) * (a3 ^ 3 * x ^ 7) + (1) * (a2 * a5 * x ^ 3) + (1 / 2) * (a2 * a5 ^ 2 * x ^ 10) + (1) * (a2 * a4 * x ^ 2) + (1) * (a2 * a4 * a5 * x ^ 9) + (1 / 2) * (a2 * a4 ^ 2 * x ^ 8) + (1) * (a2 * a3 * x) + (1) * (a2 * a3 * a5 * x ^ 8) + (1) * (a2 * a3 * a4 * x ^ 7) + (1 / 2) * (a2 * a3 ^ 2 * x ^ 6) + (1 / 2) * (a2 ^ 2) + (1 / 2) * (a2 ^ 2 * a5 * x ^ 7) + (1 / 2) * (a2 ^ 2 * a4 * x ^ 6) + (1 / 2) * (a2 ^ 2 * a3 * x ^ 5) + (1 / 6) * (a2 ^ 3 * x ^ 4) + (1) * (a1 * a5 * x ^ 2) + (1 / 2) * (a1 * a5 ^ 2 * x ^ 9) + (1) * (a1 * a4 * x) + (1) * (a1 * a4 * a5 * x ^ 8) + (1 / 2) * (a1 * a4 ^ 2 * x ^ 7) + (1) * (a1 * a3) + (1) * (a1 * a3 * a5 * x ^ 7) + (1) * (a1 * a3 * a4 * x ^ 6) + (1 / 2) * (a1 * a3 ^ 2 * x ^ 5) + (1) * (a1 * a2 * a5 * x ^ 6) + (1) * (a1 * a2 * a4 * x ^ 5) + (1) * (a1 * a2 * a3 * x ^ 4) + (1 / 2) * (a1 * a2 ^ 2 * x ^ 3) + (1 / 2) * (a1 ^ 2 * a5 * x ^ 5) + (1 / 2) * (a1 ^ 2 * a4 * x ^ 4) + (1 / 2) * (a1 ^ 2 * a3 * x ^ 3) + (1 / 2) * (a1 ^ 2 * a2 * x ^ 2) + (1 / 6) * (a1 ^ 3 * x) + (1) * (a0 * a5 * x) + (1 / 2) * (a0 * a5 ^ 2 * x ^ 8) + (1) * (a0 * a4) + (1) * (a0 * a4 * a5 * x ^ 7) + (1 / 2) * (a0 * a4 ^ 2 * x ^ 6) + (1) * (a0 * a3 * a5 * x ^ 6) + (1) * (a0 * a3 * a4 * x ^ 5) + (1 / 2) * (a0 * a3 ^ 2 * x ^ 4) + (1) * (a0 * a2 * a5 * x ^ 5) + (1) * (a0 * a2 * a4 * x ^ 4) + (1) * (a0 * a2 * a3 * x ^ 3) + (1 / 2) * (a0 * a2 ^ 2 * x ^ 2) + (1) * (a0 * a1 * a5 * x ^ 4) + (1) * (a0 * a1 * a4 * x ^ 3) + (1) * (a0 * a1 * a3 * x ^ 2) + (1) * (a0 * a1 * a2 * x) + (1 / 2) * (a0 * a1 ^ 2) + (1 / 2) * (a0 ^ 2 * a5 * x ^ 3) + (1 / 2) * (a0 ^ 2 * a4 * x ^ 2) + (1 / 2) * (a0 ^ 2 * a3 * x) + (1 / 2) * (a0 ^ 2 * a2) with hW
  have hWb : |W| ≤ max CW 0 := by
    have := hCW x ⟨hx0, hx1⟩
    rw [Real.norm_eq_abs] at this
    exact this.trans (le_max_left _ _)
  set P := a0 * x ^ 2 + a1 * x ^ 3 + a2 * x ^ 4 + a3 * x ^ 5 + a4 * x ^ 6 + a5 * x ^ 7 with hPdef
  have hx2 : x ^ 2 ≤ 1 := pow_le_one₀ hx0 hx1
  have hterm : ∀ (a : ℝ) (k : ℕ), |a| ≤ B → |a * x ^ (k + 2)| ≤ B * x ^ 2 := by
    intro a k ha
    have : x ^ (k + 2) ≤ x ^ 2 := pow_le_pow_of_le_one hx0 hx1 (by omega)
    rw [abs_mul, abs_of_nonneg (pow_nonneg hx0 (k + 2))]
    exact mul_le_mul ha this (pow_nonneg hx0 _) hB
  have bP : |P| ≤ 6 * B * x ^ 2 := by
    have e0 : |a0 * x ^ 2| ≤ B * x ^ 2 := hterm a0 0 ha0
    have e1 : |a1 * x ^ 3| ≤ B * x ^ 2 := hterm a1 1 ha1
    have e2 : |a2 * x ^ 4| ≤ B * x ^ 2 := hterm a2 2 ha2
    have e3 : |a3 * x ^ 5| ≤ B * x ^ 2 := hterm a3 3 ha3
    have e4 : |a4 * x ^ 6| ≤ B * x ^ 2 := hterm a4 4 ha4
    have e5 : |a5 * x ^ 7| ≤ B * x ^ 2 := hterm a5 5 ha5
    have t1 := abs_add_le (a0 * x ^ 2 + a1 * x ^ 3 + a2 * x ^ 4 + a3 * x ^ 5 + a4 * x ^ 6) (a5 * x ^ 7)
    have t2 := abs_add_le (a0 * x ^ 2 + a1 * x ^ 3 + a2 * x ^ 4 + a3 * x ^ 5) (a4 * x ^ 6)
    have t3 := abs_add_le (a0 * x ^ 2 + a1 * x ^ 3 + a2 * x ^ 4) (a3 * x ^ 5)
    have t4 := abs_add_le (a0 * x ^ 2 + a1 * x ^ 3) (a2 * x ^ 4)
    have t5 := abs_add_le (a0 * x ^ 2) (a1 * x ^ 3)
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
  have heP := Real.exp_bound bP1 (n := 4) (by norm_num)
  have hsP : ∑ m ∈ Finset.range 4, P ^ m / (m.factorial : ℝ) = 1 + P + P ^ 2 / 2 + P ^ 3 / 6 := by
    simp [Finset.sum_range_succ, Nat.factorial]
  rw [hsP] at heP
  have hP4 : |P| ^ 4 ≤ (6 * B) ^ 4 * x ^ 8 := by
    calc |P| ^ 4 ≤ (6 * B * x ^ 2) ^ 4 := pow_le_pow_left₀ (abs_nonneg _) bP 4
      _ = _ := by ring
  have key : Real.exp L - (1 + a0 * x ^ 2 + a1 * x ^ 3 + (a2 + a0 ^ 2 / 2) * x ^ 4 + (a3 + a0 * a1) * x ^ 5 + (a4 + a0 * a2 + a1 ^ 2 / 2 + a0 ^ 3 / 6) * x ^ 6 + (a5 + a0 * a3 + a1 * a2 + a0 ^ 2 * a1 / 2) * x ^ 7)
      = Real.exp P * (Real.exp E - 1) + (Real.exp P - (1 + P + P ^ 2 / 2 + P ^ 3 / 6)) + x ^ 8 * W := by
    rw [hLPE, Real.exp_add, hW, hPdef]; ring
  rw [key]
  have b1 : |Real.exp P * (Real.exp E - 1)| ≤ 3 * (2 * (K * x ^ 8)) := by
    rw [abs_mul, abs_of_pos (Real.exp_pos _)]
    exact mul_le_mul hexpP (hexpE.trans (by linarith)) (abs_nonneg _) (by norm_num)
  have b2 : |Real.exp P - (1 + P + P ^ 2 / 2 + P ^ 3 / 6)| ≤ (6 * B) ^ 4 * x ^ 8 := by
    calc _ ≤ _ := heP
      _ ≤ |P| ^ 4 * 1 := by gcongr; norm_num [Nat.factorial]
      _ ≤ _ := by rw [mul_one]; exact hP4
  have b3 : |x ^ 8 * W| ≤ x ^ 8 * max CW 0 := by
    rw [abs_mul, abs_of_nonneg (by positivity)]
    exact mul_le_mul_of_nonneg_left hWb (by positivity)
  have t1 := abs_add_le (Real.exp P * (Real.exp E - 1) + (Real.exp P - (1 + P + P ^ 2 / 2 + P ^ 3 / 6)))
    (x ^ 8 * W)
  have t2 := abs_add_le (Real.exp P * (Real.exp E - 1)) (Real.exp P - (1 + P + P ^ 2 / 2 + P ^ 3 / 6))
  have : (6 * K + (6 * B) ^ 4 + max CW 0) * x ^ 8 = 3 * (2 * (K * x ^ 8)) + (6 * B) ^ 4 * x ^ 8
    + x ^ 8 * max CW 0 := by ring
  rw [this]
  linarith only [t1, t2, b1, b2, b3]


lemma A_expand7 {α β γ : ℝ} (hα : |α| ≤ 3 / 2) (hβ : |β| ≤ 3 / 2) (hγ : |γ| ≤ 3 / 2) :
    ∃ K7 : ℝ, 0 ≤ K7 ∧ ∀ x L : ℝ, 0 ≤ x → x ≤ 1 / 20 →
    HasSum (fun k : ℕ => (-1 : ℝ) ^ (k + 2) * Y04.zn (k + 2) / ((k + 2 : ℕ) : ℝ) *
      ((α * x) ^ (k + 2) + (β * x) ^ (k + 2) - (γ * x) ^ (k + 2))) L →
    |Real.exp L - (1 + Y04.cf α β γ 0 * x ^ 2 + Y04.cf α β γ 1 * x ^ 3 + (Y04.cf α β γ 2 + Y04.cf α β γ 0 ^ 2 / 2) * x ^ 4 + (Y04.cf α β γ 3 + Y04.cf α β γ 0 * Y04.cf α β γ 1) * x ^ 5 + (Y04.cf α β γ 4 + Y04.cf α β γ 0 * Y04.cf α β γ 2 + Y04.cf α β γ 1 ^ 2 / 2 + Y04.cf α β γ 0 ^ 3 / 6) * x ^ 6 + (Y04.cf α β γ 5 + Y04.cf α β γ 0 * Y04.cf α β γ 3 + Y04.cf α β γ 1 * Y04.cf α β γ 2 + Y04.cf α β γ 0 ^ 2 * Y04.cf α β γ 1 / 2) * x ^ 7)| ≤ K7 * x ^ 8 := by
  have hB : ∀ k ≤ 5, |Y04.cf α β γ k| ≤ 52 := by
    intro k hk
    have := Y04.cf_bd hα hβ hγ k
    have : (3 / 2 : ℝ) ^ (k + 2) ≤ (3 / 2) ^ 7 := pow_le_pow_right₀ (by norm_num) (by omega)
    norm_num at this ⊢; linarith
  obtain ⟨KW, hKW, hE⟩ := exp_trunc7 (B := 52) (K := 154) (by norm_num) (by norm_num) (hB 0 (by norm_num))
    (hB 1 (by norm_num)) (hB 2 (by norm_num)) (hB 3 (by norm_num)) (hB 4 (by norm_num)) (hB 5 (by norm_num))
  refine ⟨6 * 154 + (6 * 52) ^ 4 + KW, by positivity, fun x L hx0 hx hL => ?_⟩
  have h := series_trunc7 hα hβ hγ hx0 (by linarith) hL
  have hx2 : x ^ 2 ≤ 1 / 400 := by
    have := pow_le_pow_left₀ hx0 hx 2
    norm_num at this ⊢; linarith
  have hx8 : x ^ 8 ≤ 1 / 400 := by
    have : x ^ 8 ≤ x ^ 2 := pow_le_pow_of_le_one hx0 (by linarith) (by norm_num)
    linarith
  exact hE x L hx0 (by linarith) (by nlinarith) (by nlinarith) h


set_option maxHeartbeats 8000000 in
/-- the `x⁶` coefficient of the regularised Beta identity at `a = αx`, `b = βx` -/
theorem coef6 {α β : ℝ} (hα : 0 < α) (hβ : 0 < β) (h3 : α + β ≤ 3) :
    (1 / 120) * ((α + β) ^ 1 * (α ^ 5 + β ^ 5)) * Y04.I 5 1 + (1 / 48) * ((α + β) ^ 2 * (α ^ 4 + β ^ 4)) * Y04.I 4 2 + (1 / 36) * ((α + β) ^ 3 * (α ^ 3 + β ^ 3)) * Y04.I 3 3 + (1 / 48) * ((α + β) ^ 4 * (α ^ 2 + β ^ 2)) * Y04.I 2 4 + (1 / 120) * ((α + β) ^ 5 * (α ^ 1 + β ^ 1)) * Y04.I 1 5 + (1 / 720) * ((α + β) ^ 6 * (α ^ 0 + β ^ 0)) * Y04.I 0 6 = (α + β) / (α * β) * 128 * (Y04.cf (α / 2) (β / 2) ((α + β) / 2) 5 + Y04.cf (α / 2) (β / 2) ((α + β) / 2) 0 * Y04.cf (α / 2) (β / 2) ((α + β) / 2) 3 + Y04.cf (α / 2) (β / 2) ((α + β) / 2) 1 * Y04.cf (α / 2) (β / 2) ((α + β) / 2) 2 + Y04.cf (α / 2) (β / 2) ((α + β) / 2) 0 ^ 2 * Y04.cf (α / 2) (β / 2) ((α + β) / 2) 1 / 2) := by
  set γ := α + β with hγ

  set c0 := Y04.cf (α / 2) (β / 2) (γ / 2) 0 with hc0
  set c1 := Y04.cf (α / 2) (β / 2) (γ / 2) 1 with hc1
  set c2 := Y04.cf (α / 2) (β / 2) (γ / 2) 2 with hc2
  set c3 := Y04.cf (α / 2) (β / 2) (γ / 2) 3 with hc3
  set c4 := Y04.cf (α / 2) (β / 2) (γ / 2) 4 with hc4
  set c5 := Y04.cf (α / 2) (β / 2) (γ / 2) 5 with hc5
  set KB := ∫ s in Set.Ioi (0 : ℝ), B8 s with hKB
  have hKB0 : 0 ≤ KB := setIntegral_nonneg measurableSet_Ioi (fun s hs => B8_nonneg (le_of_lt hs))
  have hα3 : α ≤ 3 := by linarith
  have hβ3 : β ≤ 3 := by linarith
  have hγ0 : 0 < γ := by rw [hγ]; linarith
  obtain ⟨K7, hK7, hA7⟩ := A_expand7 (α := α / 2) (β := β / 2) (γ := γ / 2)
      (by rw [abs_of_pos (by positivity)]; linarith) (by rw [abs_of_pos (by positivity)]; linarith)
      (by rw [abs_of_pos (by positivity)]; linarith)
  set e := γ / (α * β) with he
  have he0 : 0 ≤ e := by positivity
  set d1 := (-1) * (γ ^ 1 * (α ^ 0 + β ^ 0)) * Y04.I 0 1 - e * 4 * c0 with hd1
  set d2 := (1) * (γ ^ 1 * (α ^ 1 + β ^ 1)) * Y04.I 1 1 + (1 / 2) * (γ ^ 2 * (α ^ 0 + β ^ 0)) * Y04.I 0 2 - e * 8 * c1 with hd2
  set d3 := (-1 / 2) * (γ ^ 1 * (α ^ 2 + β ^ 2)) * Y04.I 2 1 + (-1 / 2) * (γ ^ 2 * (α ^ 1 + β ^ 1)) * Y04.I 1 2 + (-1 / 6) * (γ ^ 3 * (α ^ 0 + β ^ 0)) * Y04.I 0 3 - e * 16 * (c2 + c0 ^ 2 / 2) with hd3
  set d4 := (1 / 6) * (γ ^ 1 * (α ^ 3 + β ^ 3)) * Y04.I 3 1 + (1 / 4) * (γ ^ 2 * (α ^ 2 + β ^ 2)) * Y04.I 2 2 + (1 / 6) * (γ ^ 3 * (α ^ 1 + β ^ 1)) * Y04.I 1 3 + (1 / 24) * (γ ^ 4 * (α ^ 0 + β ^ 0)) * Y04.I 0 4 - e * 32 * (c3 + c0 * c1) with hd4
  set d5 := (-1 / 24) * (γ ^ 1 * (α ^ 4 + β ^ 4)) * Y04.I 4 1 + (-1 / 12) * (γ ^ 2 * (α ^ 3 + β ^ 3)) * Y04.I 3 2 + (-1 / 12) * (γ ^ 3 * (α ^ 2 + β ^ 2)) * Y04.I 2 3 + (-1 / 24) * (γ ^ 4 * (α ^ 1 + β ^ 1)) * Y04.I 1 4 + (-1 / 120) * (γ ^ 5 * (α ^ 0 + β ^ 0)) * Y04.I 0 5 - e * 64 * (c4 + c0 * c2 + c1 ^ 2 / 2 + c0 ^ 3 / 6) with hd5
  set d6 := (1 / 120) * (γ ^ 1 * (α ^ 5 + β ^ 5)) * Y04.I 5 1 + (1 / 48) * (γ ^ 2 * (α ^ 4 + β ^ 4)) * Y04.I 4 2 + (1 / 36) * (γ ^ 3 * (α ^ 3 + β ^ 3)) * Y04.I 3 3 + (1 / 48) * (γ ^ 4 * (α ^ 2 + β ^ 2)) * Y04.I 2 4 + (1 / 120) * (γ ^ 5 * (α ^ 1 + β ^ 1)) * Y04.I 1 5 + (1 / 720) * (γ ^ 6 * (α ^ 0 + β ^ 0)) * Y04.I 0 6 - e * 128 * (c5 + c0 * c3 + c1 * c2 + c0 ^ 2 * c1 / 2) with hd6
  have hmain : ∀ x : ℝ, 0 < x → x ≤ 1 / 40 →
      |d1 * x + d2 * x ^ 2 + d3 * x ^ 3 + d4 * x ^ 4 + d5 * x ^ 5 + d6 * x ^ 6|
        ≤ (48 * 2187 * KB + e * 256 * K7) * x ^ 7 := by
    intro x hx0 hx1
    have hx3 : 3 * x ≤ 1 := by linarith
    have gx : γ * x ≤ 3 * x := mul_le_mul_of_nonneg_right (by rw [hγ]; exact h3) hx0.le
    have ax : α * x ≤ 3 * x := mul_le_mul_of_nonneg_right hα3 hx0.le
    have bx : β * x ≤ 3 * x := mul_le_mul_of_nonneg_right hβ3 hx0.le
    have dx1 : (α - β) * x ≤ 3 * x := mul_le_mul_of_nonneg_right (by linarith) hx0.le
    have dx2 : -3 * x ≤ (α - β) * x := mul_le_mul_of_nonneg_right (by linarith) hx0.le
    obtain ⟨_, E1⟩ := int_expand6 (u := γ * x) (v := α * x) (t := 3 * x) (by positivity)
      gx (by positivity) ax hx3
    obtain ⟨_, E2⟩ := int_expand6 (u := γ * x) (v := β * x) (t := 3 * x) (by positivity)
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
    have hA := hA7 (2 * x) _ (by positivity) (by linarith) h
    rw [Real.exp_log hpos] at hA
    rw [show α / 2 * (2 * x) = α * x by ring, show β / 2 * (2 * x) = β * x by ring,
      show γ / 2 * (2 * x) = γ * x by ring] at hA
    rw [← hc0, ← hc1, ← hc2, ← hc3, ← hc4, ← hc5] at hA
    set Rg := Real.Gamma (1 + α * x) * Real.Gamma (1 + β * x) / Real.Gamma (1 + γ * x) with hRg
    set Iα := ∫ s in Set.Ioi (0 : ℝ), (Real.exp (-(γ * x) * Y04.ell s) - 1) * Real.exp (-(α * x) * s)
    set Iβ := ∫ s in Set.Ioi (0 : ℝ), (Real.exp (-(γ * x) * Y04.ell s) - 1) * Real.exp (-(β * x) * s)
    set u := γ * x with hu
    set Sα := (-1) * (u ^ 1 * (α * x) ^ 0) * Y04.I 0 1 + (1) * (u ^ 1 * (α * x) ^ 1) * Y04.I 1 1 + (-1 / 2) * (u ^ 1 * (α * x) ^ 2) * Y04.I 2 1 + (1 / 6) * (u ^ 1 * (α * x) ^ 3) * Y04.I 3 1 + (-1 / 24) * (u ^ 1 * (α * x) ^ 4) * Y04.I 4 1 + (1 / 120) * (u ^ 1 * (α * x) ^ 5) * Y04.I 5 1 + (1 / 2) * (u ^ 2 * (α * x) ^ 0) * Y04.I 0 2 + (-1 / 2) * (u ^ 2 * (α * x) ^ 1) * Y04.I 1 2 + (1 / 4) * (u ^ 2 * (α * x) ^ 2) * Y04.I 2 2 + (-1 / 12) * (u ^ 2 * (α * x) ^ 3) * Y04.I 3 2 + (1 / 48) * (u ^ 2 * (α * x) ^ 4) * Y04.I 4 2 + (-1 / 6) * (u ^ 3 * (α * x) ^ 0) * Y04.I 0 3 + (1 / 6) * (u ^ 3 * (α * x) ^ 1) * Y04.I 1 3 + (-1 / 12) * (u ^ 3 * (α * x) ^ 2) * Y04.I 2 3 + (1 / 36) * (u ^ 3 * (α * x) ^ 3) * Y04.I 3 3 + (1 / 24) * (u ^ 4 * (α * x) ^ 0) * Y04.I 0 4 + (-1 / 24) * (u ^ 4 * (α * x) ^ 1) * Y04.I 1 4 + (1 / 48) * (u ^ 4 * (α * x) ^ 2) * Y04.I 2 4 + (-1 / 120) * (u ^ 5 * (α * x) ^ 0) * Y04.I 0 5 + (1 / 120) * (u ^ 5 * (α * x) ^ 1) * Y04.I 1 5 + (1 / 720) * (u ^ 6 * (α * x) ^ 0) * Y04.I 0 6 with hSα
    set Sβ := (-1) * (u ^ 1 * (β * x) ^ 0) * Y04.I 0 1 + (1) * (u ^ 1 * (β * x) ^ 1) * Y04.I 1 1 + (-1 / 2) * (u ^ 1 * (β * x) ^ 2) * Y04.I 2 1 + (1 / 6) * (u ^ 1 * (β * x) ^ 3) * Y04.I 3 1 + (-1 / 24) * (u ^ 1 * (β * x) ^ 4) * Y04.I 4 1 + (1 / 120) * (u ^ 1 * (β * x) ^ 5) * Y04.I 5 1 + (1 / 2) * (u ^ 2 * (β * x) ^ 0) * Y04.I 0 2 + (-1 / 2) * (u ^ 2 * (β * x) ^ 1) * Y04.I 1 2 + (1 / 4) * (u ^ 2 * (β * x) ^ 2) * Y04.I 2 2 + (-1 / 12) * (u ^ 2 * (β * x) ^ 3) * Y04.I 3 2 + (1 / 48) * (u ^ 2 * (β * x) ^ 4) * Y04.I 4 2 + (-1 / 6) * (u ^ 3 * (β * x) ^ 0) * Y04.I 0 3 + (1 / 6) * (u ^ 3 * (β * x) ^ 1) * Y04.I 1 3 + (-1 / 12) * (u ^ 3 * (β * x) ^ 2) * Y04.I 2 3 + (1 / 36) * (u ^ 3 * (β * x) ^ 3) * Y04.I 3 3 + (1 / 24) * (u ^ 4 * (β * x) ^ 0) * Y04.I 0 4 + (-1 / 24) * (u ^ 4 * (β * x) ^ 1) * Y04.I 1 4 + (1 / 48) * (u ^ 4 * (β * x) ^ 2) * Y04.I 2 4 + (-1 / 120) * (u ^ 5 * (β * x) ^ 0) * Y04.I 0 5 + (1 / 120) * (u ^ 5 * (β * x) ^ 1) * Y04.I 1 5 + (1 / 720) * (u ^ 6 * (β * x) ^ 0) * Y04.I 0 6 with hSβ
    set Q := 1 + c0 * (2 * x) ^ 2 + c1 * (2 * x) ^ 3 + (c2 + c0 ^ 2 / 2) * (2 * x) ^ 4 + (c3 + c0 * c1) * (2 * x) ^ 5 + (c4 + c0 * c2 + c1 ^ 2 / 2 + c0 ^ 3 / 6) * (2 * x) ^ 6 + (c5 + c0 * c3 + c1 * c2 + c0 ^ 2 * c1 / 2) * (2 * x) ^ 7 with hQ
    have hpoly : d1 * x + d2 * x ^ 2 + d3 * x ^ 3 + d4 * x ^ 4 + d5 * x ^ 5 + d6 * x ^ 6
        = (Sα - Iα) + (Sβ - Iβ) + ((Iα + Iβ) - e * (4 * c0 * x ^ 1 + 8 * c1 * x ^ 2 + 16 * (c2 + c0 ^ 2 / 2) * x ^ 3 + 32 * (c3 + c0 * c1) * x ^ 4 + 64 * (c4 + c0 * c2 + c1 ^ 2 / 2 + c0 ^ 3 / 6) * x ^ 5 + 128 * (c5 + c0 * c3 + c1 * c2 + c0 ^ 2 * c1 / 2) * x ^ 6)) := by
      rw [hd1, hd2, hd3, hd4, hd5, hd6, hSα, hSβ, hu]; ring
    have hlast : (Iα + Iβ) - e * (4 * c0 * x ^ 1 + 8 * c1 * x ^ 2 + 16 * (c2 + c0 ^ 2 / 2) * x ^ 3 + 32 * (c3 + c0 * c1) * x ^ 4 + 64 * (c4 + c0 * c2 + c1 ^ 2 / 2 + c0 ^ 3 / 6) * x ^ 5 + 128 * (c5 + c0 * c3 + c1 * c2 + c0 ^ 2 * c1 / 2) * x ^ 6) = e / x * (Rg - Q) := by
      rw [hC, he, hQ, hu]
      field_simp
      ring
    rw [hpoly, hlast]
    have b3 : |e / x * (Rg - Q)| ≤ e * 256 * K7 * x ^ 7 := by
      rw [abs_mul, abs_of_nonneg (by positivity)]
      calc e / x * |Rg - Q| ≤ e / x * (K7 * (2 * x) ^ 8) :=
            mul_le_mul_of_nonneg_left hA (by positivity)
        _ = e * 256 * K7 * x ^ 7 := by field_simp; ring
    have b1 : |Sα - Iα| ≤ 24 * (3 * x) ^ 7 * KB := by rw [abs_sub_comm]; exact E1
    have b2 : |Sβ - Iβ| ≤ 24 * (3 * x) ^ 7 * KB := by rw [abs_sub_comm]; exact E2
    have t1 := abs_add_le ((Sα - Iα) + (Sβ - Iβ)) (e / x * (Rg - Q))
    have t2 := abs_add_le (Sα - Iα) (Sβ - Iβ)
    have : (48 * 2187 * KB + e * 256 * K7) * x ^ 7
        = 24 * (3 * x) ^ 7 * KB + 24 * (3 * x) ^ 7 * KB + e * 256 * K7 * x ^ 7 := by ring
    linarith
  obtain ⟨_, _, _, _, _, h6⟩ := sextic_zero (by norm_num) (by norm_num) (by positivity) hmain
  rw [hd6] at h6
  linarith


lemma D72 : 2 * ((1 / 2) * Y04.I 5 1 + Y04.I 4 2 / 2) =
    24 * z7 + 24 * ((π ^ 2 / 6) * Z02.z5) + 84 / 5 * ((π ^ 2 / 6) ^ 2 * Y04.z3) := by
  have h1 := coef6 (α := 1) (β := 1) one_pos one_pos (by norm_num)
  have h2 := coef6 (α := 1) (β := 2) one_pos two_pos (by norm_num)
  have h3 := coef6 (α := 1) (β := 1 / 3) one_pos (by norm_num) (by norm_num)
  norm_num [Y04.cf, Y04.z_two, Y04.zn_three, Y04.z_four, Z02.zn_five, V02.z_six, zn_seven] at h1 h2 h3
  linear_combination 216 * h1 + (-512 / 15) * h2 + (19683 / 10) * h3

lemma D71 : 2 * (-(25 / 48) * Y04.I 5 1 + (-(23 / 24)) * Y04.I 4 2 / 2 + (1 / 2) * Y04.I 3 3 / 6
      + Y04.I 2 4 / 24) =
    -23 * z7 - 25 * ((π ^ 2 / 6) * Z02.z5) - 181 / 10 * ((π ^ 2 / 6) ^ 2 * Y04.z3) := by
  have h1 := coef6 (α := 1) (β := 1) one_pos one_pos (by norm_num)
  have h2 := coef6 (α := 1) (β := 2) one_pos two_pos (by norm_num)
  have h3 := coef6 (α := 1) (β := 1 / 3) one_pos (by norm_num) (by norm_num)
  norm_num [Y04.cf, Y04.z_two, Y04.zn_three, Y04.z_four, Z02.zn_five, V02.z_six, zn_seven] at h1 h2 h3
  linear_combination (-1885 / 8) * h1 + (1672 / 45) * h2 + (-341901 / 160) * h3

lemma D70 : 2 * (721 / 3840 * Y04.I 5 1 + 121 / 384 * Y04.I 4 2 / 2 + (-(17 / 48)) * Y04.I 3 3 / 6
      + (-(5 / 8)) * Y04.I 2 4 / 24 + (1 / 2) * Y04.I 1 5 / 120 + Y04.I 0 6 / 720) =
    129 / 16 * z7 + 141 / 16 * ((π ^ 2 / 6) * Z02.z5) + 1071 / 160 * ((π ^ 2 / 6) ^ 2 * Y04.z3) := by
  have h1 := coef6 (α := 1) (β := 1) one_pos one_pos (by norm_num)
  have h2 := coef6 (α := 1) (β := 2) one_pos two_pos (by norm_num)
  have h3 := coef6 (α := 1) (β := 1 / 3) one_pos (by norm_num) (by norm_num)
  norm_num [Y04.cf, Y04.z_two, Y04.zn_three, Y04.z_four, Z02.zn_five, V02.z_six, zn_seven] at h1 h2 h3
  linear_combination (5695 / 64) * h1 + (-14) * h2 + (102789 / 128) * h3



lemma G_d5' (k : ℕ) : iteratedDeriv 5 (Y04.G k) 0 =
    120 * ((-1 : ℝ) ^ k * Real.sqrt 2 / 4 * (721 / 3840 + 25 / 24 * (k : ℝ) ^ 2 + 1 / 3 * (k : ℝ) ^ 4)) := by
  have := G_d5 k; linarith

lemma theta7 (k : ℕ) :
    (-1 : ℝ) ^ k * (π / Real.sqrt 2) *
        (8 / π * ∑ j ∈ Finset.Icc 1 6,
          iteratedDeriv (6 - j) (Y04.G k) 0 / ((6 - j).factorial : ℝ) / (j.factorial : ℝ) *
            Y04.I (6 - j) j) =
      (129 / 16 + 46 * (k : ℝ) ^ 2 + 16 * (k : ℝ) ^ 4) * z7
        + (141 / 16 + 50 * (k : ℝ) ^ 2 + 16 * (k : ℝ) ^ 4) * ((π ^ 2 / 6) * Z02.z5)
        + (1071 / 160 + 181 / 5 * (k : ℝ) ^ 2 + 56 / 5 * (k : ℝ) ^ 4)
            * ((π ^ 2 / 6) ^ 2 * Y04.z3) := by
  set E := (-1 : ℝ) ^ k with hE
  set Y := (721 / 3840 + 25 / 24 * (k : ℝ) ^ 2 + 1 / 3 * (k : ℝ) ^ 4) * Y04.I 5 1
    + (121 / 384 + 23 / 12 * (k : ℝ) ^ 2 + 2 / 3 * (k : ℝ) ^ 4) / 2 * Y04.I 4 2
    + (-(17 / 48 + (k : ℝ) ^ 2)) / 6 * Y04.I 3 3
    + (-(2 * (k : ℝ) ^ 2 + 5 / 8)) / 24 * Y04.I 2 4
    + (1 / 2) / 120 * Y04.I 1 5 + 1 / 720 * Y04.I 0 6 with hY
  have hsum : ∑ j ∈ Finset.Icc 1 6,
          iteratedDeriv (6 - j) (Y04.G k) 0 / ((6 - j).factorial : ℝ) / (j.factorial : ℝ) *
            Y04.I (6 - j) j = E * Real.sqrt 2 / 4 * Y := by
    rw [show Finset.Icc 1 6 = {1, 2, 3, 4, 5, 6} by decide]
    rw [Finset.sum_insert (by decide), Finset.sum_insert (by decide), Finset.sum_insert (by decide),
      Finset.sum_insert (by decide), Finset.sum_insert (by decide), Finset.sum_singleton]
    simp only [show 6 - 1 = 5 from rfl, show 6 - 2 = 4 from rfl, show 6 - 3 = 3 from rfl,
      show 6 - 4 = 2 from rfl, show 6 - 5 = 1 from rfl, show 6 - 6 = 0 from rfl,
      iteratedDeriv_zero, iteratedDeriv_one]
    rw [G_d5', V02.G_d4', Z02.G_d3', Z02.G_d2', Y04.G_d1, Y04.G_zero, hY]
    simp [Nat.factorial]
    ring
  rw [hsum]
  have hE2 : E ^ 2 = 1 := by rw [hE, ← pow_mul, mul_comm, pow_mul, neg_one_sq, one_pow]
  have hs : Real.sqrt 2 ≠ 0 := by positivity
  have hpi : π ≠ 0 := Real.pi_pos.ne'
  have e : E * (π / Real.sqrt 2) * (8 / π * (E * Real.sqrt 2 / 4 * Y)) = 2 * E ^ 2 * Y := by
    field_simp; ring
  rw [e, hE2, hY]
  linear_combination D70 - 2 * (k : ℝ) ^ 2 * D71 + 2 / 3 * (k : ℝ) ^ 4 * D72

lemma sum_eval7 (Gf : ℝ → ℝ) (τ : ℝ) :
    ∑ n ∈ Finset.Icc 2 7, τ ^ n * ∑ j ∈ Finset.Icc 1 (n - 1),
        iteratedDeriv (n - 1 - j) Gf 0 / ((n - 1 - j).factorial : ℝ) / (j.factorial : ℝ) * Y04.I (n - 1 - j) j
      = τ ^ 2 * Gf 0 * Y04.I 0 1 + τ ^ 3 * deriv Gf 0 * Y04.I 1 1 + τ ^ 3 * (Gf 0 / 2) * Y04.I 0 2
        + τ ^ 4 * (iteratedDeriv 2 Gf 0 / 2) * Y04.I 2 1
        + τ ^ 4 * (deriv Gf 0 / 2) * Y04.I 1 2 + τ ^ 4 * (Gf 0 / 6) * Y04.I 0 3
        + τ ^ 5 * (iteratedDeriv 3 Gf 0 / 6) * Y04.I 3 1 + τ ^ 5 * (iteratedDeriv 2 Gf 0 / 4) * Y04.I 2 2
        + τ ^ 5 * (deriv Gf 0 / 6) * Y04.I 1 3 + τ ^ 5 * (Gf 0 / 24) * Y04.I 0 4
        + τ ^ 6 * ∑ j ∈ Finset.Icc 1 5,
          iteratedDeriv (5 - j) Gf 0 / ((5 - j).factorial : ℝ) / (j.factorial : ℝ) * Y04.I (5 - j) j
        + τ ^ 7 * ∑ j ∈ Finset.Icc 1 6,
          iteratedDeriv (6 - j) Gf 0 / ((6 - j).factorial : ℝ) / (j.factorial : ℝ) * Y04.I (6 - j) j := by
  rw [Finset.sum_Icc_succ_top (by norm_num : 2 ≤ 6 + 1), V02.sum_eval6]

set_option maxHeartbeats 4000000 in
lemma fourier7 (k : ℕ) : ∃ C τ₀ : ℝ, 0 < τ₀ ∧ ∀ p : ℝ, 1 < p → Y04.tau p ≤ τ₀ →
    |Y04.a k (Y04.h p) - (-1 : ℝ) ^ k * (Real.sqrt 2 / π) *
        (-(4 / (16 * (k : ℝ) ^ 2 - 1)) + (π ^ 2 / 6) * Y04.tau p ^ 2 + Y04.z3 * Y04.tau p ^ 3
          - (33 / 16 + 7 * (k : ℝ) ^ 2) * (π ^ 4 / 90) * Y04.tau p ^ 4
          + (-(3 / 4 + 4 * (k : ℝ) ^ 2) * Z02.z5
              - (7 / 4 + 4 * (k : ℝ) ^ 2) * ((π ^ 2 / 6) * Y04.z3)) * Y04.tau p ^ 5
          + ((4039 / 256 + 733 / 8 * (k : ℝ) ^ 2 + 31 * (k : ℝ) ^ 4) * (π ^ 6 / 945)
              - (7 / 8 + 2 * (k : ℝ) ^ 2) * Y04.z3 ^ 2) * Y04.tau p ^ 6
          + ((129 / 16 + 46 * (k : ℝ) ^ 2 + 16 * (k : ℝ) ^ 4) * z7
        + (141 / 16 + 50 * (k : ℝ) ^ 2 + 16 * (k : ℝ) ^ 4) * ((π ^ 2 / 6) * Z02.z5)
        + (1071 / 160 + 181 / 5 * (k : ℝ) ^ 2 + 56 / 5 * (k : ℝ) ^ 4)
            * ((π ^ 2 / 6) ^ 2 * Y04.z3)) * Y04.tau p ^ 7)|
      ≤ C * Y04.tau p ^ 8 := by
  obtain ⟨CW, hW⟩ := watson7 (Y04.G k) (Y04.G_contDiff k) ⟨1, fun x _ => Y04.G_bdd k x⟩
  refine ⟨8 / π * CW, 1 / 2, by norm_num, fun p hp _ => ?_⟩
  have hpi := Real.pi_pos
  have hτ0 := Y04.tau_pos hp
  have hτ1 := Y04.tau_lt_one hp
  obtain ⟨_, hb⟩ := hW (Y04.tau p) hτ0 hτ1
  rw [sum_eval7] at hb
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
  set Θ := (4039 / 256 + 733 / 8 * (k : ℝ) ^ 2 + 31 * (k : ℝ) ^ 4) * (π ^ 6 / 945)
        - (7 / 8 + 2 * (k : ℝ) ^ 2) * Y04.z3 ^ 2 with hΘ
  set Θ7 := (129 / 16 + 46 * (k : ℝ) ^ 2 + 16 * (k : ℝ) ^ 4) * z7
        + (141 / 16 + 50 * (k : ℝ) ^ 2 + 16 * (k : ℝ) ^ 4) * ((π ^ 2 / 6) * Z02.z5)
        + (1071 / 160 + 181 / 5 * (k : ℝ) ^ 2 + 56 / 5 * (k : ℝ) ^ 4)
            * ((π ^ 2 / 6) ^ 2 * Y04.z3) with hΘ7
  have hth := V02.theta6 k
  rw [← hC6, ← hΘ] at hth
  have hth7 := theta7 k
  rw [← hC7, ← hΘ7] at hth7
  have hE2 : ((-1 : ℝ) ^ k) ^ 2 = 1 := by rw [← pow_mul, mul_comm, pow_mul, neg_one_sq, one_pow]
  have hsq : Real.sqrt 2 ^ 2 = 2 := Real.sq_sqrt (by norm_num)
  have hs0 : Real.sqrt 2 ≠ 0 := by positivity
  have hC6v : C6 = (-1 : ℝ) ^ k * Real.sqrt 2 / 8 * Θ := by
    rw [← hth]; field_simp
    linear_combination (-C6) * hE2
  have hC7v : C7 = (-1 : ℝ) ^ k * Real.sqrt 2 / 8 * Θ7 := by
    rw [← hth7]; field_simp
    linear_combination (-C7) * hE2
  set S := τ ^ 2 * Y04.G k 0 * Y04.I 0 1 + τ ^ 3 * deriv (Y04.G k) 0 * Y04.I 1 1
        + τ ^ 3 * (Y04.G k 0 / 2) * Y04.I 0 2
        + τ ^ 4 * (iteratedDeriv 2 (Y04.G k) 0 / 2) * Y04.I 2 1
        + τ ^ 4 * (deriv (Y04.G k) 0 / 2) * Y04.I 1 2 + τ ^ 4 * (Y04.G k 0 / 6) * Y04.I 0 3
        + τ ^ 5 * (iteratedDeriv 3 (Y04.G k) 0 / 6) * Y04.I 3 1
        + τ ^ 5 * (iteratedDeriv 2 (Y04.G k) 0 / 4) * Y04.I 2 2
        + τ ^ 5 * (deriv (Y04.G k) 0 / 6) * Y04.I 1 3 + τ ^ 5 * (Y04.G k 0 / 24) * Y04.I 0 4
        + τ ^ 6 * C6 + τ ^ 7 * C7 with hS
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
              - (7 / 4 + 4 * (k : ℝ) ^ 2) * ((π ^ 2 / 6) * Y04.z3)) * τ ^ 5 + Θ * τ ^ 6 + Θ7 * τ ^ 7)
        = 8 / π * (X - S) := by
    have hd2 : iteratedDeriv 2 (Y04.G k) 0 = 2 * -((-1 : ℝ) ^ k * Real.sqrt 2 * ((k : ℝ) ^ 2 / 2 + 5 / 32)) :=
      Z02.G_d2' k
    rw [ha, h1, hS, hC6v, hC7v, hd2, Z02.G_d3', Y04.G_zero k, Y04.G_d1 k, Y04.I01, Y04.I11, Y04.I02, Y04.I21,
      hI03, Z02.I31, hI22, hI13]
    field_simp
    ring
  rw [key, abs_mul, abs_of_pos (by positivity)]
  calc 8 / π * |X - S| ≤ 8 / π * (CW * τ ^ 8) := mul_le_mul_of_nonneg_left hb (by positivity)
    _ = _ := by ring

lemma zn_pos {k : ℕ} (hk : 1 < k) : 0 < Y04.zn k :=
  hasSum_lt (f := fun _ => (0 : ℝ)) (i := 0) (fun n => by positivity) (by norm_num) hasSum_zero
    (Y04.z_hasSum hk)


end Q03

theorem Q03.proof_1 : Q03.stmt_1 := Q03.G_d5

theorem Q03.proof_2 : Q03.stmt_2 := fun Gf hG hb => Q03.watson7 Gf hG hb

theorem Q03.proof_3 : Q03.stmt_3 := Q03.I51

theorem Q03.proof_4 : Q03.stmt_4 := Q03.D72

theorem Q03.proof_5 : Q03.stmt_5 := Q03.D71

theorem Q03.proof_6 : Q03.stmt_6 := Q03.D70

theorem Q03.proof_7 : Q03.stmt_7 := Q03.theta7

theorem Q03.proof_8 : Q03.stmt_8 := by
  have hpi := Real.pi_pos
  have h7 : 0 < Q03.z7 := by have := Q03.zn_pos (k := 7) (by norm_num); rwa [Q03.zn_seven] at this
  have h5 : 0 < Z02.z5 := by have := Q03.zn_pos (k := 5) (by norm_num); rwa [Z02.zn_five] at this
  have h3 : 0 < Y04.z3 := by have := Q03.zn_pos (k := 3) (by norm_num); rwa [Y04.zn_three] at this
  have n7 : (16 * Q03.z7 + 16 * ((π ^ 2 / 6) * Z02.z5) + 56 / 5 * ((π ^ 2 / 6) ^ 2 * Y04.z3) : ℝ) ≠ 0 := by
    positivity
  refine ⟨Polynomial.C (16 * Q03.z7 + 16 * ((π ^ 2 / 6) * Z02.z5) + 56 / 5 * ((π ^ 2 / 6) ^ 2 * Y04.z3))
        * Polynomial.X ^ 2
      + Polynomial.C (46 * Q03.z7 + 50 * ((π ^ 2 / 6) * Z02.z5) + 181 / 5 * ((π ^ 2 / 6) ^ 2 * Y04.z3))
        * Polynomial.X
      + Polynomial.C (129 / 16 * Q03.z7 + 141 / 16 * ((π ^ 2 / 6) * Z02.z5)
        + 1071 / 160 * ((π ^ 2 / 6) ^ 2 * Y04.z3)), ?_, ?_, ?_⟩
  · rw [Polynomial.degree_quadratic n7]; rfl
  · rw [Polynomial.leadingCoeff_quadratic n7]
    have := Q03.D72
    norm_num [Nat.factorial]
    linear_combination (-(2 / 3) : ℝ) * this
  · intro k
    obtain ⟨C, τ₀, hτ₀, hC⟩ := Q03.fourier7 k
    refine ⟨C, τ₀, hτ₀, fun p hp hτ => ?_⟩
    have := hC p hp hτ
    simp only [Polynomial.eval_add, Polynomial.eval_mul, Polynomial.eval_C, Polynomial.eval_X,
      Polynomial.eval_pow]
    convert this using 4
    ring

theorem Q03.proof_9 : Q03.stmt_9 := Q03.fourier7
