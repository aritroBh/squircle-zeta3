import Mathlib
import Zeta3.Wave5.Y04_Fourier
import Zeta3.Wave6.Z02_Theta5
/-!
# V02 — the `τ⁶` coefficient of the Fourier coefficients `a_k(h_p)`: `Θ₆`, and `prop:degree` at `n = 6`

Target id: V02.
Paper labels: `prop:theta` and `prop:degree` (main.tex lines 1214–1225), extending
`thm:fourier` (main.tex 1204–1211, formalised in `Zeta3.Wave5.Y04_Fourier`) and its
`τ⁵` extension (`Zeta3.Wave6.Z02_Theta5`) by one more order.
Detailed proofs: paper/ieee/proofs/sec_g2_theta.tex (the display of `Θ₆` at line 17,
Lemma `g2:lem:t:red`, Table `g2:tab:D` rows `D̃_{6,0}`, `D̃_{6,1}`, `D̃_{6,2}`, the remark
"the top coefficients ... at $n=6,8$") and paper/ieee/proofs/sec_g2_degree.tex
(Lemmas `g2:lem:d:coef`, `g2:lem:d:e`, `g2:lem:d:theta`, proof of Proposition 2, case `n = 2M`);
the Watson expansion is Lemma `g2:lem:f:watson` of sec_g2_fourier.tex.

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

Supplement sec_g2_theta.tex line 17:

  \Theta_6&=\bigl(\tfrac{4039}{256}+\tfrac{733}8u+31u^2\bigr)\zeta(6)-\bigl(\tfrac78+2u\bigr)\zeta(3)^2,

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

  with $\rho_m(u)=\sum_{r=0}^{\lfloor m/2\rfloor}\frac{(-4u)^r}{(2r)!}\,e_{r,m}$ and the values
  "$m=0$: $1,\frac12,-\frac58,-\frac{17}{48},\frac{121}{384}$; $m=1$ ($r\ge2$):
  $1,\frac12,-\frac{23}{24},-\frac{25}{48},\frac{731}{1152}$; $m=2$ ($r\ge4$):
  $1,\frac12,-\frac{31}{24},-\frac{11}{16},\frac{409}{384}$"
  (so `ρ₄(u) = e_{0,4} - 2u e_{1,4} + (2/3)u² e_{2,4} = 121/384 + (23/12)u + (2/3)u²`).

Supplement, Lemma `g2:lem:t:red` and Table `g2:tab:D` (products `ζ(6), ζ(3)²`):

  $\wD_{n,m}=2\sum_{j=1}^{n-1-2m}\dfrac{e_{m,\,n-1-j}}{j!}\,I_{n-1-j,\,j}$
  $\Theta_n(u)=\sum_{m=0}^{\lfloor n/2\rfloor-1}\wD_{n,m}\,\dfrac{(-4u)^m}{(2m)!}$.
  $\wD_{6,0}$ & $\tfrac{4039}{256},\ -\tfrac78$ \\
  $\wD_{6,1}$ & $-\tfrac{733}{16},\ 1$ \\
  $\wD_{6,2}$ & $\tfrac{93}2,\ 0$ \\
  "the top coefficients $\frac{(-4)^{M-1}}{(2M-2)!}\wD_{n,M-1}$ at $n=6,8$ are
  $\frac{(-4)^2}{4!}\cdot\frac{93}2\zeta(6)=31\zeta(6)$"

Supplement, proof of Proposition 2, case `n = 2M` (here `M = 3`):

  the coefficient of $u^{M-1}$ is $2\cdot\frac{(-4)^{M-1}}{(2M-2)!}\,I_{2M-2,1}=2(-4)^{M-1}\eta(2M)$,
  ... this is $(-1)^{M-1}(2^{2M-1}-1)\zeta(2M)$, which is not $0$.

Reused (not redefined): from `Zeta3.Wave5.Y04_Fourier`: `Y04.tau` (`τ = 1 - 1/p`), `Y04.h`
(the support function `h_p`), `Y04.a` (`a_k(f) = (4/π)∫_0^{π/2} f cos 4kθ`), `Y04.G` (the layer
function `G_k(x) = cos(4kθ_x) cos²θ_x sin θ_x`, `θ_x = arctan e^{-x}`), `Y04.I`
(`I_{m,j} = ∫_0^∞ s^m log^j(1+e^{-s}) ds`), `Y04.z3` (`ζ(3) = (riemannZeta 3).re`); from
`Zeta3.Wave6.Z02_Theta5`: `Z02.z5` (`ζ(5) = (riemannZeta 5).re`). Z02 proves `g_{k,3}`, Watson
at `N = 5`, the rows `D̃_{5,0}`, `D̃_{5,1}` and the `τ⁵` expansion (`Z02.stmt_8`).

Conventions.
* `ζ(3) = Y04.z3 = (riemannZeta 3).re` and `ζ(5) = Z02.z5 = (riemannZeta 5).re`, the real parts
  of Mathlib's complex zeta values; by `zeta_nat_eq_tsum_of_gt_one` these are the real numbers
  `∑ 1/n^3`, `∑ 1/n^5`. The even values are written out: `ζ(2) = π²/6`, `ζ(4) = π⁴/90`,
  `ζ(6) = π⁶/945`.
* `g_{k,m} = G_k^{(m)}(0)/m!` is written `iteratedDeriv m (Y04.G k) 0 / m!`; `G_k` is `C^∞`
  on `ℝ` (proved in Y04), so these are genuine derivatives.
* All half-line integrals `I_{m,j}`, `j ≥ 1`, have integrable integrands (proved in Y04), and
  `stmt_2` asserts the Watson integrand integrable, so no statement holds by the
  "non-integrable ↦ 0" convention.
* As in Y04/Z02, the Watson lemma (`stmt_2`) asks `G` to be `C^∞` on all of `ℝ` and bounded on
  `[0,∞)` (the supplement: `C^∞` on `[0,∞)`); this covers `G = G_k`.
* "The coefficient of `τ⁶`" in `prop:degree` is the coefficient of an asymptotic expansion
  (the supplement's "Setting"). `stmt_7` states it that way: a polynomial `Θ₆ ∈ ℝ[X]` such
  that `(-1)^k (√2/π)[-4/(16k²-1) + Σ_{n=2}^{5} Θ_n(k²) τ^n + Θ₆(k²) τ⁶]` (with the proved
  `Θ₂, …, Θ₅` of Y04/Z02 written out) approximates `a_k(h_p)` to `O_k(τ⁷)` for every `k`.
  Such coefficients are unique and a polynomial is determined by its values at `k²`, so this
  pins down the paper's `Θ₆`. "Degree exactly 2" is `Polynomial.degree Θ₆ = 2`.
* The `O(τ⁷)` claims are explicit: `∃ C τ₀, 0 < τ₀ ∧ ∀ p > 1, τ ≤ τ₀ → |error| ≤ C τ⁷`, the
  same form as `Y04.stmt_8` and `Z02.stmt_8`.
* Numerical check before freezing (mpmath, 40 digits): `I_{4,1} = 93ζ(6)/4`, the rows
  `D̃_{6,0}, D̃_{6,1}, D̃_{6,2}`, `g_{k,4}` for `k = 0..3` and `stmt_6` for `k = 0..3` agree to
  `≤ 2·10⁻³⁶`; direct quadrature of `a_k(h_p)` at `τ = 0.02, 0.01, 0.005` for `k = 0, 1, 2`
  gives `(residual after τ⁵)/τ⁶ - Θ₆(k²) ≈ c·τ` with `c` bounded (≈ 44, 357, 2140 at
  `τ = 0.005`), consistent with an `O(τ⁷)` error.

NOT COVERED: `prop:degree` for `n ≠ 6` beyond Z02 and in general, and its clause "coefficients
are alternating multiple zeta values of weight `n`" (no AMZV/iterated-integral library is set
up; for `n = 6` the explicit closed form in `ζ(6), ζ(3)²` is stated instead); `Θ₇, Θ₈, Θ₉`
and `prop:theta` for `7 ≤ n ≤ 9`; the iterated-integral form of `I_{m,j}`, the
shuffle/stuffle/duality relations and the exact certificates (the three `n = 6` identities
`D̃_{6,0}`, `D̃_{6,1}`, `D̃_{6,2}` are stated directly); the general formula for `e_{r,m}`,
`ρ_m` (only `g_{k,4}` is stated); Cauchy's formula `P = ∫ h_p` and the `k = 0` comparison
with `thm:pw`.
-/

-- BEGIN STATEMENTS
namespace V02

open Real MeasureTheory Set Polynomial

/-- Lemma `g2:lem:d:e`(b) at `m = 4` (with `ρ₄(u) = 121/384 + (23/12)u + (2/3)u²`):
the next Taylor coefficient of `G_k`,
`g_{k,4} = G_k''''(0)/4! = (-1)^k (√2/4)(121/384 + (23/12)k² + (2/3)k⁴)`. -/
def stmt_1 : Prop :=
  ∀ k : ℕ, iteratedDeriv 4 (Y04.G k) 0 / 24 =
    (-1 : ℝ) ^ k * Real.sqrt 2 / 4 *
      (121 / 384 + 23 / 12 * (k : ℝ) ^ 2 + 2 / 3 * (k : ℝ) ^ 4)

/-- Lemma `g2:lem:f:watson` with `N = 6`: for `G` smooth and bounded on `[0,∞)`, with
`g_m = G^{(m)}(0)/m!`, there is `C` with, for `0 < τ < 1`, the integrand integrable and
`|τ ∫_0^∞ [(1+e^{-s})^τ - 1] G(τ s) ds - Σ_{n=2}^{6} τ^n Σ_{j=1}^{n-1} g_{n-1-j}/j! I_{n-1-j,j}|
  ≤ C τ⁷`. -/
def stmt_2 : Prop :=
  ∀ Gf : ℝ → ℝ, ContDiff ℝ (⊤ : ℕ∞) Gf → (∃ M : ℝ, ∀ x : ℝ, 0 ≤ x → |Gf x| ≤ M) →
    ∃ C : ℝ, ∀ τ : ℝ, 0 < τ → τ < 1 →
      IntegrableOn (fun s : ℝ => ((1 + Real.exp (-s)) ^ τ - 1) * Gf (τ * s)) (Ioi 0) ∧
      |τ * (∫ s in Ioi (0 : ℝ), ((1 + Real.exp (-s)) ^ τ - 1) * Gf (τ * s)) -
          ∑ n ∈ Finset.Icc 2 6, τ ^ n * ∑ j ∈ Finset.Icc 1 (n - 1),
            iteratedDeriv (n - 1 - j) Gf 0 / ((n - 1 - j).factorial : ℝ) /
              (j.factorial : ℝ) * Y04.I (n - 1 - j) j| ≤ C * τ ^ 7

/-- Table `g2:tab:D`, row `D̃_{6,2}` (with `e_{2,4} = 1`): `D̃_{6,2} = 2 I_{4,1} = 93 ζ(6)/2`,
i.e. `I_{4,1} = 4! η(6) = 93 ζ(6)/4`, `ζ(6) = π⁶/945`. -/
def stmt_3 : Prop := 2 * Y04.I 4 1 = 93 / 2 * (π ^ 6 / 945)

/-- Table `g2:tab:D`, row `D̃_{6,1}`: with `e_{1,2..4} = 1, 1/2, -23/24`,
`D̃_{6,1} = 2(e_{1,4} I_{4,1} + e_{1,3} I_{3,2}/2! + e_{1,2} I_{2,3}/3!)
  = -(733/16) ζ(6) + ζ(3)²`. -/
def stmt_4 : Prop :=
  2 * (-(23 / 24) * Y04.I 4 1 + (1 / 2) * Y04.I 3 2 / 2 + Y04.I 2 3 / 6) =
    -(733 / 16) * (π ^ 6 / 945) + Y04.z3 ^ 2

/-- Table `g2:tab:D`, row `D̃_{6,0}`: with `e_{0,0..4} = 1, 1/2, -5/8, -17/48, 121/384`,
`D̃_{6,0} = 2(e_{0,4} I_{4,1} + e_{0,3} I_{3,2}/2! + e_{0,2} I_{2,3}/3! + e_{0,1} I_{1,4}/4!
  + e_{0,0} I_{0,5}/5!) = (4039/256) ζ(6) - (7/8) ζ(3)²`. -/
def stmt_5 : Prop :=
  2 * (121 / 384 * Y04.I 4 1 + (-(17 / 48)) * Y04.I 3 2 / 2 + (-(5 / 8)) * Y04.I 2 3 / 6
      + (1 / 2) * Y04.I 1 4 / 24 + Y04.I 0 5 / 120) =
    4039 / 256 * (π ^ 6 / 945) - 7 / 8 * Y04.z3 ^ 2

/-- Lemmas `g2:lem:d:coef`, `g2:lem:d:theta` and `g2:lem:t:red` at `n = 6`, together with the
closed form of `Θ₆` (sec_g2_theta.tex line 17): the `τ⁶` coefficient
`c_6(k) = (8/π) Σ_{j=1}^{5} g_{k,5-j}/j! I_{5-j,j}` satisfies
`(-1)^k (π/√2) c_6(k) = Θ₆(k²) = (4039/256 + (733/8)k² + 31k⁴) ζ(6) - (7/8 + 2k²) ζ(3)²`. -/
def stmt_6 : Prop :=
  ∀ k : ℕ,
    (-1 : ℝ) ^ k * (π / Real.sqrt 2) *
        (8 / π * ∑ j ∈ Finset.Icc 1 5,
          iteratedDeriv (5 - j) (Y04.G k) 0 / ((5 - j).factorial : ℝ) / (j.factorial : ℝ) *
            Y04.I (5 - j) j) =
      (4039 / 256 + 733 / 8 * (k : ℝ) ^ 2 + 31 * (k : ℝ) ^ 4) * (π ^ 6 / 945)
        - (7 / 8 + 2 * (k : ℝ) ^ 2) * Y04.z3 ^ 2

/-- `prop:degree` at `n = 6 = 2m`, `m = 3`: there is a polynomial `Θ₆` of degree exactly
`⌊6/2⌋ - 1 = 2`, with leading coefficient `(-1)^{m-1}(2^{2m-1}-1) ζ(2m) = 31 ζ(6)`, whose value
at `u = k²` is the `τ⁶` coefficient of `(-1)^k (π/√2) a_k(h_p)` in the asymptotic sense
(the lower coefficients `Θ₂, …, Θ₅` written out as proved in Y04 and Z02): for every `k`
there are `C_k`, `τ₀ > 0` with
`|a_k(h_p) - (-1)^k (√2/π)[-4/(16k²-1) + ζ(2)τ² + ζ(3)τ³ - (33/16 + 7k²)ζ(4)τ⁴ + Θ₅(k²)τ⁵
   + Θ₆(k²)τ⁶]| ≤ C_k τ⁷` for `p > 1`, `τ ≤ τ₀`. -/
def stmt_7 : Prop :=
  ∃ Θ₆ : ℝ[X],
    Θ₆.degree = ((6 / 2 - 1 : ℕ) : WithBot ℕ) ∧
    Θ₆.leadingCoeff = (-1 : ℝ) ^ (3 - 1) * (2 ^ (2 * 3 - 1) - 1) * (π ^ 6 / 945) ∧
    ∀ k : ℕ, ∃ C τ₀ : ℝ, 0 < τ₀ ∧ ∀ p : ℝ, 1 < p → Y04.tau p ≤ τ₀ →
      |Y04.a k (Y04.h p) - (-1 : ℝ) ^ k * (Real.sqrt 2 / π) *
          (-(4 / (16 * (k : ℝ) ^ 2 - 1)) + (π ^ 2 / 6) * Y04.tau p ^ 2 + Y04.z3 * Y04.tau p ^ 3
            - (33 / 16 + 7 * (k : ℝ) ^ 2) * (π ^ 4 / 90) * Y04.tau p ^ 4
            + (-(3 / 4 + 4 * (k : ℝ) ^ 2) * Z02.z5
                - (7 / 4 + 4 * (k : ℝ) ^ 2) * ((π ^ 2 / 6) * Y04.z3)) * Y04.tau p ^ 5
            + Θ₆.eval ((k : ℝ) ^ 2) * Y04.tau p ^ 6)|
        ≤ C * Y04.tau p ^ 7

/-- `thm:fourier` extended by two orders with `prop:theta`'s `Θ₅` and the supplement's `Θ₆`:
for every `k ≥ 0` there are `C_k` and `τ₀ > 0` such that for every `p > 1` with
`τ = 1 - 1/p ≤ τ₀`, with `u = k²`,
`|a_k(h_p) - (-1)^k (√2/π)[-4/(16k²-1) + ζ(2)τ² + ζ(3)τ³ - (33/16 + 7k²)ζ(4)τ⁴ + Θ₅(u)τ⁵
   + Θ₆(u)τ⁶]| ≤ C_k τ⁷`, where `Θ₅(u) = -(3/4 + 4u) ζ(5) - (7/4 + 4u) ζ(2) ζ(3)` and
`Θ₆(u) = (4039/256 + (733/8)u + 31u²) ζ(6) - (7/8 + 2u) ζ(3)²`. -/
def stmt_8 : Prop :=
  ∀ k : ℕ, ∃ C τ₀ : ℝ, 0 < τ₀ ∧ ∀ p : ℝ, 1 < p → Y04.tau p ≤ τ₀ →
    |Y04.a k (Y04.h p) - (-1 : ℝ) ^ k * (Real.sqrt 2 / π) *
        (-(4 / (16 * (k : ℝ) ^ 2 - 1)) + (π ^ 2 / 6) * Y04.tau p ^ 2 + Y04.z3 * Y04.tau p ^ 3
          - (33 / 16 + 7 * (k : ℝ) ^ 2) * (π ^ 4 / 90) * Y04.tau p ^ 4
          + (-(3 / 4 + 4 * (k : ℝ) ^ 2) * Z02.z5
              - (7 / 4 + 4 * (k : ℝ) ^ 2) * ((π ^ 2 / 6) * Y04.z3)) * Y04.tau p ^ 5
          + ((4039 / 256 + 733 / 8 * (k : ℝ) ^ 2 + 31 * (k : ℝ) ^ 4) * (π ^ 6 / 945)
              - (7 / 8 + 2 * (k : ℝ) ^ 2) * Y04.z3 ^ 2) * Y04.tau p ^ 6)|
      ≤ C * Y04.tau p ^ 7

end V02
-- END STATEMENTS

open Real MeasureTheory Set

namespace V02

/-- a trigonometric monomial `r · (cos or sin)(K t) · cos^a t · sin^b t` -/
noncomputable def tm (K t : ℝ) (x : ℝ × Bool × ℕ × ℕ) : ℝ :=
  x.1 * (if x.2.1 then Real.cos (K * t) else Real.sin (K * t)) * Real.cos t ^ x.2.2.1 * Real.sin t ^ x.2.2.2

noncomputable def ev (K : ℝ) (L : List (ℝ × Bool × ℕ × ℕ)) (t : ℝ) : ℝ := (L.map (tm K t)).sum

/-- the derivative of a monomial, as a list of monomials -/
noncomputable def dterm (K : ℝ) (x : ℝ × Bool × ℕ × ℕ) : List (ℝ × Bool × ℕ × ℕ) :=
  if x.2.1 then
    [(-(K * x.1), false, x.2.2.1, x.2.2.2), (-(x.1 * x.2.2.1), true, x.2.2.1 - 1, x.2.2.2 + 1),
      (x.1 * x.2.2.2, true, x.2.2.1 + 1, x.2.2.2 - 1)]
  else
    [(K * x.1, true, x.2.2.1, x.2.2.2), (-(x.1 * x.2.2.1), false, x.2.2.1 - 1, x.2.2.2 + 1),
      (x.1 * x.2.2.2, false, x.2.2.1 + 1, x.2.2.2 - 1)]

/-- multiplication by `θ' = -(sin θ cos θ)` -/
def shift (x : ℝ × Bool × ℕ × ℕ) : ℝ × Bool × ℕ × ℕ := (-x.1, x.2.1, x.2.2.1 + 1, x.2.2.2 + 1)

noncomputable def step (K : ℝ) (L : List (ℝ × Bool × ℕ × ℕ)) : List (ℝ × Bool × ℕ × ℕ) :=
  (L.flatMap (dterm K)).map shift

lemma tm_hasDeriv (K t : ℝ) (x : ℝ × Bool × ℕ × ℕ) :
    HasDerivAt (fun t => tm K t x) (ev K (dterm K x) t) t := by
  obtain ⟨r, e, a, b⟩ := x
  have hS := ((hasDerivAt_id t).const_mul K).sin
  have hC := ((hasDerivAt_id t).const_mul K).cos
  have hc := Real.hasDerivAt_cos t
  have hs := Real.hasDerivAt_sin t
  cases e
  · have h := ((hS.const_mul r).mul (hc.pow a)).mul (hs.pow b)
    convert h using 1
    · funext s; simp [tm, mul_assoc]
    · simp only [ev, dterm, tm, List.map, List.sum_cons, List.sum_nil, id, Bool.false_eq_true,
        ite_false, ite_true, Pi.mul_apply, Pi.pow_apply]
      push_cast; ring
  · have h := ((hC.const_mul r).mul (hc.pow a)).mul (hs.pow b)
    convert h using 1
    · funext s; simp [tm, mul_assoc]
    · simp only [ev, dterm, tm, List.map, List.sum_cons, List.sum_nil, id, ite_true, Pi.mul_apply, Pi.pow_apply]
      push_cast; ring

lemma ev_hasDeriv (K t : ℝ) (L : List (ℝ × Bool × ℕ × ℕ)) :
    HasDerivAt (ev K L) (ev K (L.flatMap (dterm K)) t) t := by
  induction L with
  | nil => simp only [ev, List.map_nil, List.sum_nil, List.flatMap_nil]; exact hasDerivAt_const _ _
  | cons x L ih =>
    have h := (tm_hasDeriv K t x).add ih
    convert h using 1
    · funext s; simp [ev]
    · simp [ev, List.flatMap_cons, List.map_append, List.sum_append]

lemma ev_shift (K t : ℝ) (L : List (ℝ × Bool × ℕ × ℕ)) :
    ev K (L.map shift) t = ev K L t * (-(Real.sin t * Real.cos t)) := by
  induction L with
  | nil => simp [ev]
  | cons x L ih =>
    simp only [ev, List.map_cons, List.sum_cons] at ih ⊢
    rw [ih]
    obtain ⟨r, e, a, b⟩ := x
    cases e <;> simp only [tm, shift, Bool.false_eq_true, ite_false, ite_true] <;> ring

lemma iter_G (k : ℕ) (n : ℕ) :
    iteratedDeriv n (Y04.G k) = fun x => ev (4 * k) ((step (4 * k))^[n] [(1, true, 2, 1)]) (Y04.thetaX x) := by
  induction n with
  | zero =>
    funext x
    simp [ev, tm, Y04.G]
  | succ n ih =>
    rw [iteratedDeriv_succ, ih]
    funext x
    have h := (ev_hasDeriv (4 * k) (Y04.thetaX x) ((step (4 * k))^[n] [(1, true, 2, 1)])).comp x
      (Y04.theta_hasDeriv x)
    rw [show (fun x => ev (4 * (k : ℝ)) ((step (4 * k))^[n] [(1, true, 2, 1)]) (Y04.thetaX x))
      = ev (4 * (k : ℝ)) ((step (4 * k))^[n] [(1, true, 2, 1)]) ∘ Y04.thetaX from rfl, h.deriv,
      Function.iterate_succ_apply', step, ev_shift]


lemma G_d4 (k : ℕ) : iteratedDeriv 4 (Y04.G k) 0 / 24 =
    (-1 : ℝ) ^ k * Real.sqrt 2 / 4 * (121 / 384 + 23 / 12 * (k : ℝ) ^ 2 + 2 / 3 * (k : ℝ) ^ 4) := by
  rw [iter_G]
  simp only [Y04.thetaX_zero, Function.iterate_succ_apply', Function.iterate_zero_apply, step, dterm,
    shift, List.flatMap_cons, List.flatMap_nil, List.map_cons, List.map_nil, List.append_nil,
    List.cons_append, List.nil_append, ite_true, ite_false, Bool.false_eq_true]
  simp only [ev, tm, List.map_cons, List.map_nil, List.sum_cons, List.sum_nil, ite_true, ite_false,
    Bool.false_eq_true, Y04.cos_k, Y04.sin_k, Real.cos_pi_div_four, Real.sin_pi_div_four]
  norm_num
  have : Real.sqrt 2 ^ 2 = 2 := Real.sq_sqrt (by norm_num)
  have h11 : Real.sqrt 2 ^ 11 = 32 * Real.sqrt 2 := by
    rw [show Real.sqrt 2 ^ 11 = (Real.sqrt 2 ^ 2) ^ 5 * Real.sqrt 2 by ring, this]; ring
  ring_nf
  rw [h11]; ring

end V02

namespace V02

lemma bernoulli'_six : bernoulli' 6 = 1 / 42 := by
  have h2 : Nat.choose 6 2 = 15 := by decide
  have h3 : Nat.choose 6 3 = 20 := by decide
  have h4 : Nat.choose 6 4 = 15 := by decide
  have h5 : Nat.choose 5 2 = 10 := by decide
  have h6 : Nat.choose 5 3 = 10 := by decide
  have h7 : Nat.choose 4 2 = 6 := by decide
  have h8 : bernoulli' 5 = 0 := bernoulli'_eq_zero_of_odd (by decide) (by norm_num)
  rw [bernoulli'_def]
  norm_num [Finset.sum_range_succ, Finset.sum_range_zero, h2, h3, h4, h5, h6, h7, h8]

lemma z_six : Y04.zn 6 = π ^ 6 / 945 := by
  have h := riemannZeta_two_mul_nat (k := 3) (by norm_num)
  unfold Y04.zn
  rw [show ((6 : ℕ) : ℂ) = 2 * ((3 : ℕ) : ℂ) by norm_num, h,
    bernoulli_eq_bernoulli'_of_ne_one (by decide), bernoulli'_six]
  have e : (π : ℂ) ^ (2 * 3) = ((π ^ 6 : ℝ) : ℂ) := by push_cast; ring
  rw [e]
  norm_num [Nat.factorial]
  rw [← Complex.ofReal_pow, Complex.ofReal_re]
  ring

lemma I41 : Y04.I 4 1 = 93 / 4 * (π ^ 6 / 945) := by
  rw [Y04.I_m1, Y04.eta_val 6 (by norm_num), z_six]; norm_num [Nat.factorial]; ring

end V02

namespace V02

/-- Taylor's formula at `0` to order 4 with a quintic remainder on `[0,∞)`. -/
lemma taylor4 {Gf : ℝ → ℝ} (hG : ContDiff ℝ (⊤ : ℕ∞) Gf) {M : ℝ}
    (hM : ∀ x : ℝ, 0 ≤ x → |Gf x| ≤ M) :
    ∃ C : ℝ, 0 ≤ C ∧ ∀ x : ℝ, 0 ≤ x →
      |Gf x - (iteratedDeriv 0 Gf 0 / (Nat.factorial 0 : ℝ) + iteratedDeriv 1 Gf 0 / (Nat.factorial 1 : ℝ) * x
        + iteratedDeriv 2 Gf 0 / (Nat.factorial 2 : ℝ) * x ^ 2
        + iteratedDeriv 3 Gf 0 / (Nat.factorial 3 : ℝ) * x ^ 3
        + iteratedDeriv 4 Gf 0 / (Nat.factorial 4 : ℝ) * x ^ 4)| ≤ C * x ^ 5 := by
  set g0 := iteratedDeriv 0 Gf 0 / (Nat.factorial 0 : ℝ)
  set g1 := iteratedDeriv 1 Gf 0 / (Nat.factorial 1 : ℝ)
  set g2 := iteratedDeriv 2 Gf 0 / (Nat.factorial 2 : ℝ)
  set g3 := iteratedDeriv 3 Gf 0 / (Nat.factorial 3 : ℝ)
  set g4 := iteratedDeriv 4 Gf 0 / (Nat.factorial 4 : ℝ)
  have hc5 : Continuous (iteratedDeriv 5 Gf) :=
    hG.continuous_iteratedDeriv 5 (WithTop.coe_le_coe.2 le_top)
  obtain ⟨K5, hK5⟩ := (isCompact_Icc (a := (0 : ℝ)) (b := 1)).exists_bound_of_continuousOn
    hc5.continuousOn
  have hM0 : 0 ≤ M := (abs_nonneg _).trans (hM 0 le_rfl)
  refine ⟨max (|K5| / 120) (M + |g0| + |g1| + |g2| + |g3| + |g4|), le_max_of_le_left (by positivity), ?_⟩
  intro x hx
  rcases eq_or_lt_of_le hx with h0 | hpos
  · subst h0; simp [g0]
  rcases le_or_gt x 1 with hx1 | hx1
  · have hcd : ContDiffOn ℝ (4 + 1) Gf (uIcc 0 x) := (hG.of_le (WithTop.coe_le_coe.2 le_top)).contDiffOn
    obtain ⟨x', hx', heq⟩ := taylor_mean_remainder_lagrange_iteratedDeriv (n := 4) hpos.ne hcd
    have hu : UniqueDiffOn ℝ (uIcc 0 x) := uniqueDiffOn_uIcc hpos.ne
    have hmem : (0 : ℝ) ∈ uIcc 0 x := left_mem_uIcc
    have hT : taylorWithinEval Gf 4 (uIcc 0 x) 0 x
        = g0 + g1 * x + g2 * x ^ 2 + g3 * x ^ 3 + g4 * x ^ 4 := by
      rw [taylor_within_apply]
      simp only [Finset.sum_range_succ, Finset.sum_range_zero, zero_add, sub_zero, smul_eq_mul]
      rw [iteratedDerivWithin_eq_iteratedDeriv hu (hG.contDiffAt.of_le (WithTop.coe_le_coe.2 le_top)) hmem,
        iteratedDerivWithin_eq_iteratedDeriv hu (hG.contDiffAt.of_le (WithTop.coe_le_coe.2 le_top)) hmem,
        iteratedDerivWithin_eq_iteratedDeriv hu (hG.contDiffAt.of_le (WithTop.coe_le_coe.2 le_top)) hmem,
        iteratedDerivWithin_eq_iteratedDeriv hu (hG.contDiffAt.of_le (WithTop.coe_le_coe.2 le_top)) hmem,
        iteratedDerivWithin_eq_iteratedDeriv hu (hG.contDiffAt.of_le (WithTop.coe_le_coe.2 le_top)) hmem]
      simp only [g0, g1, g2, g3, g4]; norm_num; ring
    rw [hT] at heq
    rw [heq]
    rw [uIoo_of_le hx] at hx'
    have hb := hK5 x' ⟨hx'.1.le, by linarith [hx'.2]⟩
    rw [Real.norm_eq_abs] at hb
    rw [sub_zero, abs_div, abs_mul, abs_pow, abs_of_pos hpos]
    norm_num [Nat.factorial]
    have : |iteratedDeriv 5 Gf x'| * x ^ 5 / 120 ≤ |K5| / 120 * x ^ 5 := by
      have := hb.trans (le_abs_self K5)
      have h3 : 0 ≤ x ^ 5 := by positivity
      nlinarith
    calc |iteratedDeriv 5 Gf x'| * x ^ 5 / 120 ≤ |K5| / 120 * x ^ 5 := this
      _ ≤ _ := by gcongr; exact le_max_left _ _
  · have h1 : 1 ≤ x := hx1.le
    have hx2' : 1 ≤ x ^ 2 := by nlinarith
    have hx3' : 1 ≤ x ^ 3 := by nlinarith
    have hx4' : 1 ≤ x ^ 4 := by nlinarith
    have hx15 : x ≤ x ^ 5 := by
      rw [show x ^ 5 = x * x ^ 4 by ring]; nlinarith
    have hx25 : x ^ 2 ≤ x ^ 5 := by
      rw [show x ^ 5 = x ^ 2 * x ^ 3 by ring]; nlinarith
    have hx35 : x ^ 3 ≤ x ^ 5 := by
      rw [show x ^ 5 = x ^ 3 * x ^ 2 by ring]; nlinarith
    have hx45 : x ^ 4 ≤ x ^ 5 := by
      rw [show x ^ 5 = x ^ 4 * x by ring]; nlinarith
    have h15 : 1 ≤ x ^ 5 := by linarith
    calc _ ≤ |Gf x| + (|g0| + |g1| * x + |g2| * x ^ 2 + |g3| * x ^ 3 + |g4| * x ^ 4) := by
          have := abs_sub (Gf x) (g0 + g1 * x + g2 * x ^ 2 + g3 * x ^ 3 + g4 * x ^ 4)
          have a00 := abs_add_le (g0 + g1 * x + g2 * x ^ 2 + g3 * x ^ 3) (g4 * x ^ 4)
          have a0 := abs_add_le (g0 + g1 * x + g2 * x ^ 2) (g3 * x ^ 3)
          have a1 := abs_add_le (g0 + g1 * x) (g2 * x ^ 2)
          have a2 := abs_add_le g0 (g1 * x)
          rw [abs_mul, abs_of_pos hpos] at a2
          rw [abs_mul, abs_of_nonneg (by positivity : (0:ℝ) ≤ x ^ 2)] at a1
          rw [abs_mul, abs_of_nonneg (by positivity : (0:ℝ) ≤ x ^ 3)] at a0
          rw [abs_mul, abs_of_nonneg (by positivity : (0:ℝ) ≤ x ^ 4)] at a00
          linarith
      _ ≤ (M + |g0| + |g1| + |g2| + |g3| + |g4|) * x ^ 5 := by
          have := hM x hx
          have e1 : |g1| * x ≤ |g1| * x ^ 5 := mul_le_mul_of_nonneg_left hx15 (abs_nonneg _)
          have e2 : |g2| * x ^ 2 ≤ |g2| * x ^ 5 := mul_le_mul_of_nonneg_left hx25 (abs_nonneg _)
          have e5 : |g3| * x ^ 3 ≤ |g3| * x ^ 5 := mul_le_mul_of_nonneg_left hx35 (abs_nonneg _)
          have e6 : |g4| * x ^ 4 ≤ |g4| * x ^ 5 := mul_le_mul_of_nonneg_left hx45 (abs_nonneg _)
          have e3 : M ≤ M * x ^ 5 := le_mul_of_one_le_right hM0 h15
          have e4 : |g0| ≤ |g0| * x ^ 5 := le_mul_of_one_le_right (abs_nonneg _) h15
          nlinarith
      _ ≤ _ := by gcongr; exact le_max_right _ _

/-- the dominating function of the order-6 Watson remainder -/
noncomputable def B7 (s : ℝ) : ℝ := (1 + s + s ^ 2 + s ^ 3 + s ^ 4 + s ^ 5) * Real.exp (-s)

lemma B7_int : IntegrableOn B7 (Set.Ioi 0) := by
  have h := (((((Y04.intOn_pow_exp1 0).add (Y04.intOn_pow_exp1 1)).add (Y04.intOn_pow_exp1 2)).add
    (Y04.intOn_pow_exp1 3)).add (Y04.intOn_pow_exp1 4)).add (Y04.intOn_pow_exp1 5)
  refine IntegrableOn.congr_fun h (fun s _ => ?_) measurableSet_Ioi
  simp only [B7, Pi.add_apply]; ring

lemma B7_nonneg {s : ℝ} (hs : 0 ≤ s) : 0 ≤ B7 s := by unfold B7; positivity

/-- the pointwise order-6 Watson remainder -/
lemma watson_pt6 {g0 g1 g2 g3 g4 CR M τ s : ℝ} (hCR : 0 ≤ CR) (hM : 0 ≤ M) (hτ0 : 0 < τ) (hτ1 : τ < 1)
    (hs : 0 ≤ s) {Gv : ℝ} (hGv : |Gv| ≤ M)
    (hR : |Gv - (g0 + g1 * (τ * s) + g2 * (τ * s) ^ 2 + g3 * (τ * s) ^ 3 + g4 * (τ * s) ^ 4)|
      ≤ CR * (τ * s) ^ 5) :
    |τ * ((Real.exp (τ * Y04.ell s) - 1) * Gv)
      - (τ ^ 2 * g0 * (s ^ 0 * Y04.ell s ^ 1) + τ ^ 3 * g1 * (s ^ 1 * Y04.ell s ^ 1)
        + τ ^ 3 * (g0 / 2) * (s ^ 0 * Y04.ell s ^ 2) + τ ^ 4 * g2 * (s ^ 2 * Y04.ell s ^ 1)
        + τ ^ 4 * (g1 / 2) * (s ^ 1 * Y04.ell s ^ 2) + τ ^ 4 * (g0 / 6) * (s ^ 0 * Y04.ell s ^ 3)
        + τ ^ 5 * g3 * (s ^ 3 * Y04.ell s ^ 1) + τ ^ 5 * (g2 / 2) * (s ^ 2 * Y04.ell s ^ 2)
        + τ ^ 5 * (g1 / 6) * (s ^ 1 * Y04.ell s ^ 3) + τ ^ 5 * (g0 / 24) * (s ^ 0 * Y04.ell s ^ 4)
        + τ ^ 6 * g4 * (s ^ 4 * Y04.ell s ^ 1) + τ ^ 6 * (g3 / 2) * (s ^ 3 * Y04.ell s ^ 2)
        + τ ^ 6 * (g2 / 6) * (s ^ 2 * Y04.ell s ^ 3) + τ ^ 6 * (g1 / 24) * (s ^ 1 * Y04.ell s ^ 4)
        + τ ^ 6 * (g0 / 120) * (s ^ 0 * Y04.ell s ^ 5))|
      ≤ τ ^ 7 * ((|g1| + |g2| + |g3| + |g4| + 2 * CR + M) * B7 s) := by
  have hl0 := Y04.ell_nonneg s
  have hl1 := Y04.ell_le_one s hs
  have hle := Y04.ell_le s
  set L := Y04.ell s
  have hτL : |τ * L| ≤ 1 := by rw [abs_of_nonneg (by positivity)]; nlinarith
  have he := Real.exp_bound hτL (n := 6) (by norm_num)
  simp only [Finset.sum_range_succ, Finset.sum_range_zero, Nat.factorial] at he
  norm_num at he
  set a1 := τ * L + (τ * L) ^ 2 / 2 + (τ * L) ^ 3 / 6 + (τ * L) ^ 4 / 24 + (τ * L) ^ 5 / 120 with ha1
  set r1 := Real.exp (τ * L) - 1 - a1 with hr1
  have hr1b : |r1| ≤ τ ^ 6 * L ^ 6 := by
    have : r1 = Real.exp (τ * L) - (1 + τ * L + (τ * L) ^ 2 / 2 + (τ * L) ^ 3 / 6 + (τ * L) ^ 4 / 24
        + (τ * L) ^ 5 / 120) := by
      rw [hr1, ha1]; ring
    rw [this]
    have h4 : (|τ| * |L|) ^ 6 = τ ^ 6 * L ^ 6 := by
      rw [abs_of_pos hτ0, abs_of_nonneg hl0]; ring
    rw [h4] at he
    have : 0 ≤ τ ^ 6 * L ^ 6 := by positivity
    linarith
  set q := g0 + g1 * (τ * s) + g2 * (τ * s) ^ 2 + g3 * (τ * s) ^ 3 + g4 * (τ * s) ^ 4
  set R := Gv - q with hRdef
  set Q := g1 * s * (L ^ 5 / 120) + g2 * s ^ 2 * (L ^ 4 / 24 + τ * L ^ 5 / 120)
    + g3 * s ^ 3 * (L ^ 3 / 6 + τ * L ^ 4 / 24 + τ ^ 2 * L ^ 5 / 120)
    + g4 * s ^ 4 * (L ^ 2 / 2 + τ * L ^ 3 / 6 + τ ^ 2 * L ^ 4 / 24 + τ ^ 3 * L ^ 5 / 120) with hQ
  have key : τ * ((Real.exp (τ * L) - 1) * Gv)
      - (τ ^ 2 * g0 * (s ^ 0 * L ^ 1) + τ ^ 3 * g1 * (s ^ 1 * L ^ 1)
        + τ ^ 3 * (g0 / 2) * (s ^ 0 * L ^ 2) + τ ^ 4 * g2 * (s ^ 2 * L ^ 1)
        + τ ^ 4 * (g1 / 2) * (s ^ 1 * L ^ 2) + τ ^ 4 * (g0 / 6) * (s ^ 0 * L ^ 3)
        + τ ^ 5 * g3 * (s ^ 3 * L ^ 1) + τ ^ 5 * (g2 / 2) * (s ^ 2 * L ^ 2)
        + τ ^ 5 * (g1 / 6) * (s ^ 1 * L ^ 3) + τ ^ 5 * (g0 / 24) * (s ^ 0 * L ^ 4)
        + τ ^ 6 * g4 * (s ^ 4 * L ^ 1) + τ ^ 6 * (g3 / 2) * (s ^ 3 * L ^ 2)
        + τ ^ 6 * (g2 / 6) * (s ^ 2 * L ^ 3) + τ ^ 6 * (g1 / 24) * (s ^ 1 * L ^ 4)
        + τ ^ 6 * (g0 / 120) * (s ^ 0 * L ^ 5))
      = τ ^ 7 * Q + τ * a1 * R + τ * r1 * Gv := by
    rw [hRdef, hr1, ha1, hQ]; ring
  rw [key]
  have hes : 0 ≤ Real.exp (-s) := (Real.exp_pos _).le
  have hL2 : L ^ 2 ≤ Real.exp (-s) := (pow_le_of_le_one hl0 hl1 (by norm_num)).trans hle
  have hL3 : L ^ 3 ≤ Real.exp (-s) := (pow_le_of_le_one hl0 hl1 (by norm_num)).trans hle
  have hL4 : L ^ 4 ≤ Real.exp (-s) := (pow_le_of_le_one hl0 hl1 (by norm_num)).trans hle
  have hL5 : L ^ 5 ≤ Real.exp (-s) := (pow_le_of_le_one hl0 hl1 (by norm_num)).trans hle
  have hL6 : L ^ 6 ≤ Real.exp (-s) := (pow_le_of_le_one hl0 hl1 (by norm_num)).trans hle
  have hτ' : 0 ≤ τ := hτ0.le
  have hτ2 : τ ^ 2 ≤ 1 := pow_le_one₀ hτ' hτ1.le
  have hτ3 : τ ^ 3 ≤ 1 := pow_le_one₀ hτ' hτ1.le
  have hτL3 : τ * L ^ 3 ≤ L ^ 3 := mul_le_of_le_one_left (by positivity) hτ1.le
  have hτL4 : τ * L ^ 4 ≤ L ^ 4 := mul_le_of_le_one_left (by positivity) hτ1.le
  have hτL5 : τ * L ^ 5 ≤ L ^ 5 := mul_le_of_le_one_left (by positivity) hτ1.le
  have hτ2L4 : τ ^ 2 * L ^ 4 ≤ L ^ 4 := mul_le_of_le_one_left (by positivity) hτ2
  have hτ2L5 : τ ^ 2 * L ^ 5 ≤ L ^ 5 := mul_le_of_le_one_left (by positivity) hτ2
  have hτ3L5 : τ ^ 3 * L ^ 5 ≤ L ^ 5 := mul_le_of_le_one_left (by positivity) hτ3
  -- term 1
  have b1 : |τ ^ 7 * Q| ≤ τ ^ 7 * ((|g1| + |g2| + |g3| + |g4|) * (s + s ^ 2 + s ^ 3 + s ^ 4) * Real.exp (-s)) := by
    rw [abs_mul, abs_of_pos (by positivity)]
    apply mul_le_mul_of_nonneg_left _ (by positivity)
    have u1 : 0 ≤ L ^ 5 / 120 := by positivity
    have u1' : L ^ 5 / 120 ≤ Real.exp (-s) := by linarith
    have u2 : 0 ≤ L ^ 4 / 24 + τ * L ^ 5 / 120 := by positivity
    have u2' : L ^ 4 / 24 + τ * L ^ 5 / 120 ≤ Real.exp (-s) := by linarith
    have u3 : 0 ≤ L ^ 3 / 6 + τ * L ^ 4 / 24 + τ ^ 2 * L ^ 5 / 120 := by positivity
    have u3' : L ^ 3 / 6 + τ * L ^ 4 / 24 + τ ^ 2 * L ^ 5 / 120 ≤ Real.exp (-s) := by linarith
    have u4 : 0 ≤ L ^ 2 / 2 + τ * L ^ 3 / 6 + τ ^ 2 * L ^ 4 / 24 + τ ^ 3 * L ^ 5 / 120 := by positivity
    have u4' : L ^ 2 / 2 + τ * L ^ 3 / 6 + τ ^ 2 * L ^ 4 / 24 + τ ^ 3 * L ^ 5 / 120 ≤ Real.exp (-s) := by
      linarith
    have v1 : |g1 * s * (L ^ 5 / 120)| ≤ |g1| * (s * Real.exp (-s)) := by
      rw [abs_mul, abs_mul, abs_of_nonneg hs, abs_of_nonneg u1, mul_assoc]
      exact mul_le_mul_of_nonneg_left (mul_le_mul_of_nonneg_left u1' hs) (abs_nonneg _)
    have v2 : |g2 * s ^ 2 * (L ^ 4 / 24 + τ * L ^ 5 / 120)| ≤ |g2| * (s ^ 2 * Real.exp (-s)) := by
      rw [abs_mul, abs_mul, abs_of_nonneg (by positivity : (0:ℝ) ≤ s ^ 2), abs_of_nonneg u2, mul_assoc]
      exact mul_le_mul_of_nonneg_left (mul_le_mul_of_nonneg_left u2' (by positivity)) (abs_nonneg _)
    have v3 : |g3 * s ^ 3 * (L ^ 3 / 6 + τ * L ^ 4 / 24 + τ ^ 2 * L ^ 5 / 120)|
        ≤ |g3| * (s ^ 3 * Real.exp (-s)) := by
      rw [abs_mul, abs_mul, abs_of_nonneg (by positivity : (0:ℝ) ≤ s ^ 3), abs_of_nonneg u3, mul_assoc]
      exact mul_le_mul_of_nonneg_left (mul_le_mul_of_nonneg_left u3' (by positivity)) (abs_nonneg _)
    have v4 : |g4 * s ^ 4 * (L ^ 2 / 2 + τ * L ^ 3 / 6 + τ ^ 2 * L ^ 4 / 24 + τ ^ 3 * L ^ 5 / 120)|
        ≤ |g4| * (s ^ 4 * Real.exp (-s)) := by
      rw [abs_mul, abs_mul, abs_of_nonneg (by positivity : (0:ℝ) ≤ s ^ 4), abs_of_nonneg u4, mul_assoc]
      exact mul_le_mul_of_nonneg_left (mul_le_mul_of_nonneg_left u4' (by positivity)) (abs_nonneg _)
    have t0 := abs_add_le (g1 * s * (L ^ 5 / 120) + g2 * s ^ 2 * (L ^ 4 / 24 + τ * L ^ 5 / 120)
      + g3 * s ^ 3 * (L ^ 3 / 6 + τ * L ^ 4 / 24 + τ ^ 2 * L ^ 5 / 120))
      (g4 * s ^ 4 * (L ^ 2 / 2 + τ * L ^ 3 / 6 + τ ^ 2 * L ^ 4 / 24 + τ ^ 3 * L ^ 5 / 120))
    have t1 := abs_add_le (g1 * s * (L ^ 5 / 120) + g2 * s ^ 2 * (L ^ 4 / 24 + τ * L ^ 5 / 120))
      (g3 * s ^ 3 * (L ^ 3 / 6 + τ * L ^ 4 / 24 + τ ^ 2 * L ^ 5 / 120))
    have t2 := abs_add_le (g1 * s * (L ^ 5 / 120)) (g2 * s ^ 2 * (L ^ 4 / 24 + τ * L ^ 5 / 120))
    have ag1 := abs_nonneg g1; have ag2 := abs_nonneg g2; have ag3 := abs_nonneg g3
    have ag4 := abs_nonneg g4
    have w : (|g1| + |g2| + |g3| + |g4|) * (s + s ^ 2 + s ^ 3 + s ^ 4) * Real.exp (-s)
        - (|g1| * (s * Real.exp (-s)) + |g2| * (s ^ 2 * Real.exp (-s)) + |g3| * (s ^ 3 * Real.exp (-s))
          + |g4| * (s ^ 4 * Real.exp (-s)))
        = Real.exp (-s) * (|g1| * (s ^ 2 + s ^ 3 + s ^ 4) + |g2| * (s + s ^ 3 + s ^ 4)
          + |g3| * (s + s ^ 2 + s ^ 4) + |g4| * (s + s ^ 2 + s ^ 3)) := by ring
    have : 0 ≤ Real.exp (-s) * (|g1| * (s ^ 2 + s ^ 3 + s ^ 4) + |g2| * (s + s ^ 3 + s ^ 4)
          + |g3| * (s + s ^ 2 + s ^ 4) + |g4| * (s + s ^ 2 + s ^ 3)) := by
      positivity
    linarith
  -- term 2
  have b2 : |τ * a1 * R| ≤ τ ^ 7 * (2 * CR * s ^ 5 * Real.exp (-s)) := by
    have ha1b : |a1| ≤ 2 * τ * L := by
      have hu0 : 0 ≤ τ * L := mul_nonneg hτ0.le hl0
      have hu1 : τ * L ≤ 1 := by
        have := mul_le_mul hτ1.le hl1 hl0 zero_le_one; linarith
      have p2 : (τ * L) ^ 2 ≤ τ * L := pow_le_of_le_one hu0 hu1 (by norm_num)
      have p3 : (τ * L) ^ 3 ≤ τ * L := pow_le_of_le_one hu0 hu1 (by norm_num)
      have p4 : (τ * L) ^ 4 ≤ τ * L := pow_le_of_le_one hu0 hu1 (by norm_num)
      have p5 : (τ * L) ^ 5 ≤ τ * L := pow_le_of_le_one hu0 hu1 (by norm_num)
      have q2 : 0 ≤ (τ * L) ^ 2 := sq_nonneg _
      have q3 : 0 ≤ (τ * L) ^ 3 := pow_nonneg hu0 3
      have q4 : 0 ≤ (τ * L) ^ 4 := pow_nonneg hu0 4
      have q5 : 0 ≤ (τ * L) ^ 5 := pow_nonneg hu0 5
      have : 0 ≤ a1 := by rw [ha1]; linarith
      rw [abs_of_nonneg this, ha1]
      have : 2 * τ * L = 2 * (τ * L) := by ring
      linarith
    have hRb : |R| ≤ CR * (τ * s) ^ 5 := hR
    rw [abs_mul, abs_mul, abs_of_pos hτ0]
    calc τ * |a1| * |R| ≤ τ * (2 * τ * L) * (CR * (τ * s) ^ 5) :=
          mul_le_mul (mul_le_mul_of_nonneg_left ha1b hτ') hRb (abs_nonneg _) (by positivity)
      _ = τ ^ 7 * (2 * CR * s ^ 5 * L) := by ring
      _ ≤ τ ^ 7 * (2 * CR * s ^ 5 * Real.exp (-s)) :=
          mul_le_mul_of_nonneg_left (mul_le_mul_of_nonneg_left hle (by positivity)) (by positivity)
  -- term 3
  have b3 : |τ * r1 * Gv| ≤ τ ^ 7 * (M * Real.exp (-s)) := by
    rw [abs_mul, abs_mul, abs_of_pos hτ0]
    calc τ * |r1| * |Gv| ≤ τ * (τ ^ 6 * L ^ 6) * M :=
          mul_le_mul (mul_le_mul_of_nonneg_left hr1b hτ0.le) hGv (abs_nonneg _) (by positivity)
      _ = τ ^ 7 * (L ^ 6 * M) := by ring
      _ ≤ τ ^ 7 * (Real.exp (-s) * M) :=
          mul_le_mul_of_nonneg_left (mul_le_mul_of_nonneg_right hL6 hM) (by positivity)
      _ = _ := by ring
  have t1 := abs_add_le (τ ^ 7 * Q + τ * a1 * R) (τ * r1 * Gv)
  have t2 := abs_add_le (τ ^ 7 * Q) (τ * a1 * R)
  have hfin : τ ^ 7 * ((|g1| + |g2| + |g3| + |g4|) * (s + s ^ 2 + s ^ 3 + s ^ 4) * Real.exp (-s))
      + τ ^ 7 * (2 * CR * s ^ 5 * Real.exp (-s))
      + τ ^ 7 * (M * Real.exp (-s)) ≤ τ ^ 7 * ((|g1| + |g2| + |g3| + |g4| + 2 * CR + M) * B7 s) := by
    unfold B7
    apply le_of_sub_nonneg
    have : τ ^ 7 * ((|g1| + |g2| + |g3| + |g4| + 2 * CR + M)
          * ((1 + s + s ^ 2 + s ^ 3 + s ^ 4 + s ^ 5) * Real.exp (-s)))
        - (τ ^ 7 * ((|g1| + |g2| + |g3| + |g4|) * (s + s ^ 2 + s ^ 3 + s ^ 4) * Real.exp (-s))
        + τ ^ 7 * (2 * CR * s ^ 5 * Real.exp (-s)) + τ ^ 7 * (M * Real.exp (-s)))
        = τ ^ 7 * Real.exp (-s) * ((|g1| + |g2| + |g3| + |g4|) * (1 + s ^ 5)
          + 2 * CR * (1 + s + s ^ 2 + s ^ 3 + s ^ 4)
          + M * (s + s ^ 2 + s ^ 3 + s ^ 4 + s ^ 5)) := by ring
    rw [this]
    have := abs_nonneg g1; have := abs_nonneg g2; have := abs_nonneg g3; have := abs_nonneg g4
    positivity
  linarith


lemma int15 {f1 f2 f3 f4 f5 f6 f7 f8 f9 f10 f11 f12 f13 f14 f15 : ℝ → ℝ}
    (c1 c2 c3 c4 c5 c6 c7 c8 c9 c10 c11 c12 c13 c14 c15 : ℝ)
    (h1 : IntegrableOn f1 (Set.Ioi 0)) (h2 : IntegrableOn f2 (Set.Ioi 0)) (h3 : IntegrableOn f3 (Set.Ioi 0))
    (h4 : IntegrableOn f4 (Set.Ioi 0)) (h5 : IntegrableOn f5 (Set.Ioi 0)) (h6 : IntegrableOn f6 (Set.Ioi 0))
    (h7 : IntegrableOn f7 (Set.Ioi 0)) (h8 : IntegrableOn f8 (Set.Ioi 0)) (h9 : IntegrableOn f9 (Set.Ioi 0))
    (h10 : IntegrableOn f10 (Set.Ioi 0)) (h11 : IntegrableOn f11 (Set.Ioi 0))
    (h12 : IntegrableOn f12 (Set.Ioi 0)) (h13 : IntegrableOn f13 (Set.Ioi 0))
    (h14 : IntegrableOn f14 (Set.Ioi 0)) (h15 : IntegrableOn f15 (Set.Ioi 0)) :
    IntegrableOn (fun s => c1 * f1 s + c2 * f2 s + c3 * f3 s + c4 * f4 s + c5 * f5 s + c6 * f6 s
      + c7 * f7 s + c8 * f8 s + c9 * f9 s + c10 * f10 s + c11 * f11 s + c12 * f12 s + c13 * f13 s
      + c14 * f14 s + c15 * f15 s) (Set.Ioi 0) ∧
    ∫ s in Set.Ioi (0 : ℝ), (c1 * f1 s + c2 * f2 s + c3 * f3 s + c4 * f4 s + c5 * f5 s + c6 * f6 s
      + c7 * f7 s + c8 * f8 s + c9 * f9 s + c10 * f10 s + c11 * f11 s + c12 * f12 s + c13 * f13 s
      + c14 * f14 s + c15 * f15 s)
      = c1 * (∫ s in Set.Ioi (0 : ℝ), f1 s) + c2 * (∫ s in Set.Ioi (0 : ℝ), f2 s)
        + c3 * (∫ s in Set.Ioi (0 : ℝ), f3 s) + c4 * (∫ s in Set.Ioi (0 : ℝ), f4 s)
        + c5 * (∫ s in Set.Ioi (0 : ℝ), f5 s) + c6 * (∫ s in Set.Ioi (0 : ℝ), f6 s)
        + c7 * (∫ s in Set.Ioi (0 : ℝ), f7 s) + c8 * (∫ s in Set.Ioi (0 : ℝ), f8 s)
        + c9 * (∫ s in Set.Ioi (0 : ℝ), f9 s) + c10 * (∫ s in Set.Ioi (0 : ℝ), f10 s)
        + c11 * (∫ s in Set.Ioi (0 : ℝ), f11 s) + c12 * (∫ s in Set.Ioi (0 : ℝ), f12 s)
        + c13 * (∫ s in Set.Ioi (0 : ℝ), f13 s) + c14 * (∫ s in Set.Ioi (0 : ℝ), f14 s)
        + c15 * (∫ s in Set.Ioi (0 : ℝ), f15 s) := by
  obtain ⟨b10, e10⟩ := Z02.int10 c1 c2 c3 c4 c5 c6 c7 c8 c9 c10 h1 h2 h3 h4 h5 h6 h7 h8 h9 h10
  have a11 : Integrable (fun s => c11 * f11 s) (volume.restrict (Set.Ioi 0)) := h11.const_mul c11
  have a12 : Integrable (fun s => c12 * f12 s) (volume.restrict (Set.Ioi 0)) := h12.const_mul c12
  have a13 : Integrable (fun s => c13 * f13 s) (volume.restrict (Set.Ioi 0)) := h13.const_mul c13
  have a14 : Integrable (fun s => c14 * f14 s) (volume.restrict (Set.Ioi 0)) := h14.const_mul c14
  have a15 : Integrable (fun s => c15 * f15 s) (volume.restrict (Set.Ioi 0)) := h15.const_mul c15
  have b11 : Integrable (fun s => c1 * f1 s + c2 * f2 s + c3 * f3 s + c4 * f4 s + c5 * f5 s + c6 * f6 s
      + c7 * f7 s + c8 * f8 s + c9 * f9 s + c10 * f10 s + c11 * f11 s) (volume.restrict (Set.Ioi 0)) :=
    Integrable.add b10 a11
  have b12 : Integrable (fun s => c1 * f1 s + c2 * f2 s + c3 * f3 s + c4 * f4 s + c5 * f5 s + c6 * f6 s
      + c7 * f7 s + c8 * f8 s + c9 * f9 s + c10 * f10 s + c11 * f11 s + c12 * f12 s)
      (volume.restrict (Set.Ioi 0)) := b11.add a12
  have b13 : Integrable (fun s => c1 * f1 s + c2 * f2 s + c3 * f3 s + c4 * f4 s + c5 * f5 s + c6 * f6 s
      + c7 * f7 s + c8 * f8 s + c9 * f9 s + c10 * f10 s + c11 * f11 s + c12 * f12 s + c13 * f13 s)
      (volume.restrict (Set.Ioi 0)) := b12.add a13
  have b14 : Integrable (fun s => c1 * f1 s + c2 * f2 s + c3 * f3 s + c4 * f4 s + c5 * f5 s + c6 * f6 s
      + c7 * f7 s + c8 * f8 s + c9 * f9 s + c10 * f10 s + c11 * f11 s + c12 * f12 s + c13 * f13 s
      + c14 * f14 s) (volume.restrict (Set.Ioi 0)) := b13.add a14
  have b15 : Integrable (fun s => c1 * f1 s + c2 * f2 s + c3 * f3 s + c4 * f4 s + c5 * f5 s + c6 * f6 s
      + c7 * f7 s + c8 * f8 s + c9 * f9 s + c10 * f10 s + c11 * f11 s + c12 * f12 s + c13 * f13 s
      + c14 * f14 s + c15 * f15 s) (volume.restrict (Set.Ioi 0)) := b14.add a15
  refine ⟨b15, ?_⟩
  rw [integral_add b14 a15, integral_add b13 a14, integral_add b12 a13, integral_add b11 a12,
    integral_add b10 a11, e10, integral_const_mul, integral_const_mul, integral_const_mul,
    integral_const_mul, integral_const_mul]

lemma watson6 (Gf : ℝ → ℝ) (hG : ContDiff ℝ (⊤ : ℕ∞) Gf) (hb : ∃ M : ℝ, ∀ x : ℝ, 0 ≤ x → |Gf x| ≤ M) :
    ∃ C : ℝ, ∀ τ : ℝ, 0 < τ → τ < 1 →
      IntegrableOn (fun s : ℝ => ((1 + Real.exp (-s)) ^ τ - 1) * Gf (τ * s)) (Ioi 0) ∧
      |τ * (∫ s in Ioi (0 : ℝ), ((1 + Real.exp (-s)) ^ τ - 1) * Gf (τ * s)) -
          ∑ n ∈ Finset.Icc 2 6, τ ^ n * ∑ j ∈ Finset.Icc 1 (n - 1),
            iteratedDeriv (n - 1 - j) Gf 0 / ((n - 1 - j).factorial : ℝ) /
              (j.factorial : ℝ) * Y04.I (n - 1 - j) j| ≤ C * τ ^ 7 := by
  obtain ⟨M, hM⟩ := hb
  have hM0 : 0 ≤ M := (abs_nonneg _).trans (hM 0 le_rfl)
  obtain ⟨CR, hCR, hT⟩ := taylor4 hG hM
  set g0 := iteratedDeriv 0 Gf 0 / (Nat.factorial 0 : ℝ) with hg0
  set g1 := iteratedDeriv 1 Gf 0 / (Nat.factorial 1 : ℝ) with hg1
  set g2 := iteratedDeriv 2 Gf 0 / (Nat.factorial 2 : ℝ) with hg2
  set g3 := iteratedDeriv 3 Gf 0 / (Nat.factorial 3 : ℝ) with hg3
  set g4 := iteratedDeriv 4 Gf 0 / (Nat.factorial 4 : ℝ) with hg4
  refine ⟨(|g1| + |g2| + |g3| + |g4| + 2 * CR + M) * ∫ s in Set.Ioi (0 : ℝ), B7 s, fun τ hτ0 hτ1 => ?_⟩
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
  obtain ⟨hmi, hmv⟩ := int15 (τ ^ 2 * g0) (τ ^ 3 * g1) (τ ^ 3 * (g0 / 2)) (τ ^ 4 * g2) (τ ^ 4 * (g1 / 2))
    (τ ^ 4 * (g0 / 6)) (τ ^ 5 * g3) (τ ^ 5 * (g2 / 2)) (τ ^ 5 * (g1 / 6)) (τ ^ 5 * (g0 / 24))
    (τ ^ 6 * g4) (τ ^ 6 * (g3 / 2)) (τ ^ 6 * (g2 / 6)) (τ ^ 6 * (g1 / 24)) (τ ^ 6 * (g0 / 120))
    (Y04.I_int 0 1 le_rfl) (Y04.I_int 1 1 le_rfl) (Y04.I_int 0 2 (by norm_num))
    (Y04.I_int 2 1 le_rfl) (Y04.I_int 1 2 (by norm_num)) (Y04.I_int 0 3 (by norm_num))
    (Y04.I_int 3 1 le_rfl) (Y04.I_int 2 2 (by norm_num)) (Y04.I_int 1 3 (by norm_num))
    (Y04.I_int 0 4 (by norm_num))
    (Y04.I_int 4 1 le_rfl) (Y04.I_int 3 2 (by norm_num)) (Y04.I_int 2 3 (by norm_num))
    (Y04.I_int 1 4 (by norm_num)) (Y04.I_int 0 5 (by norm_num))
  have hsum : ∑ n ∈ Finset.Icc 2 6, τ ^ n * ∑ j ∈ Finset.Icc 1 (n - 1),
            iteratedDeriv (n - 1 - j) Gf 0 / ((n - 1 - j).factorial : ℝ) /
              (j.factorial : ℝ) * Y04.I (n - 1 - j) j
      = τ ^ 2 * g0 * Y04.I 0 1 + τ ^ 3 * g1 * Y04.I 1 1 + τ ^ 3 * (g0 / 2) * Y04.I 0 2
        + τ ^ 4 * g2 * Y04.I 2 1 + τ ^ 4 * (g1 / 2) * Y04.I 1 2 + τ ^ 4 * (g0 / 6) * Y04.I 0 3
        + τ ^ 5 * g3 * Y04.I 3 1 + τ ^ 5 * (g2 / 2) * Y04.I 2 2 + τ ^ 5 * (g1 / 6) * Y04.I 1 3
        + τ ^ 5 * (g0 / 24) * Y04.I 0 4
        + τ ^ 6 * g4 * Y04.I 4 1 + τ ^ 6 * (g3 / 2) * Y04.I 3 2 + τ ^ 6 * (g2 / 6) * Y04.I 2 3
        + τ ^ 6 * (g1 / 24) * Y04.I 1 4 + τ ^ 6 * (g0 / 120) * Y04.I 0 5 := by
    rw [show Finset.Icc 2 6 = {2, 3, 4, 5, 6} by decide]
    rw [Finset.sum_insert (by decide), Finset.sum_insert (by decide), Finset.sum_insert (by decide),
      Finset.sum_insert (by decide), Finset.sum_singleton]
    rw [show Finset.Icc 1 (2 - 1) = {1} by decide, show Finset.Icc 1 (3 - 1) = {1, 2} by decide,
      show Finset.Icc 1 (4 - 1) = {1, 2, 3} by decide, show Finset.Icc 1 (5 - 1) = {1, 2, 3, 4} by decide,
      show Finset.Icc 1 (6 - 1) = {1, 2, 3, 4, 5} by decide]
    simp [Finset.sum_insert, hg0, hg1, hg2, hg3, hg4, Nat.factorial]
    ring
  rw [hsum]
  unfold Y04.I
  rw [← hmv, ← integral_const_mul, ← integral_sub (hint.const_mul τ) hmi]
  have hB := (B7_int.const_mul ((|g1| + |g2| + |g3| + |g4| + 2 * CR + M))).const_mul (τ ^ 7)
  have := norm_integral_le_of_norm_le hB
    ((ae_restrict_iff' measurableSet_Ioi).2 (Filter.Eventually.of_forall
      (fun s (hs : 0 < s) => by
        rw [Real.norm_eq_abs]
        have hx : 0 ≤ τ * s := by positivity
        have := watson_pt6 (g0 := g0) (g1 := g1) (g2 := g2) (g3 := g3) (g4 := g4) hCR hM0 hτ0 hτ1 hs.le
          (hM (τ * s) hx)
          (by have := hT (τ * s) hx; simpa [hg0, hg1, hg2, hg3, hg4] using this)
        convert this using 2)))
  rw [Real.norm_eq_abs, integral_const_mul, integral_const_mul] at this
  calc _ ≤ _ := this
    _ = _ := by ring

end V02

namespace V02

lemma em5_bound {q : ℝ} (hq : 0 ≤ q) :
    |Real.exp (-q) - (1 - q + q ^ 2 / 2 - q ^ 3 / 6 + q ^ 4 / 24)| ≤ 4 * q ^ 5 := by
  rcases le_or_gt q 1 with h | h
  · have h1 : |-q| ≤ 1 := by rw [abs_neg, abs_of_nonneg hq]; exact h
    have he := Real.exp_bound h1 (n := 5) (by norm_num)
    simp only [Finset.sum_range_succ, Finset.sum_range_zero, Nat.factorial] at he
    norm_num at he
    rw [abs_of_nonneg hq] at he
    have e : Real.exp (-q) - (1 - q + q ^ 2 / 2 - q ^ 3 / 6 + q ^ 4 / 24)
        = Real.exp (-q) - (1 + -q + q ^ 2 / 2 + (-q) ^ 3 / 6 + q ^ 4 / 24) := by ring
    rw [e]
    have : 0 ≤ q ^ 5 := by positivity
    nlinarith
  · have he : Real.exp (-q) ≤ 1 := Real.exp_le_one_iff.2 (by linarith)
    have he0 := Real.exp_pos (-q)
    have h2 : q ≤ q ^ 2 := by nlinarith
    have h3 : q ^ 2 ≤ q ^ 3 := by nlinarith
    have h4 : q ^ 3 ≤ q ^ 4 := by nlinarith
    have h5 : q ^ 4 ≤ q ^ 5 := by nlinarith
    rw [abs_le]
    constructor <;> nlinarith

lemma ep6_bound {p : ℝ} (hp0 : 0 ≤ p) (hp1 : p ≤ 1) :
    |Real.exp (-p) - 1 - (-p + p ^ 2 / 2 - p ^ 3 / 6 + p ^ 4 / 24 - p ^ 5 / 120)| ≤ p ^ 6 := by
  have h1 : |-p| ≤ 1 := by rw [abs_neg, abs_of_nonneg hp0]; exact hp1
  have he := Real.exp_bound h1 (n := 6) (by norm_num)
  simp only [Finset.sum_range_succ, Finset.sum_range_zero, Nat.factorial] at he
  norm_num at he
  rw [abs_of_nonneg hp0] at he
  have e : Real.exp (-p) - 1 - (-p + p ^ 2 / 2 - p ^ 3 / 6 + p ^ 4 / 24 - p ^ 5 / 120)
      = Real.exp (-p) - (1 + -p + p ^ 2 / 2 + (-p) ^ 3 / 6 + p ^ 4 / 24 + (-p) ^ 5 / 120) := by ring
  rw [e]
  have : 0 ≤ p ^ 6 := by positivity
  nlinarith

lemma mono_bd6 {u v t L s : ℝ} (hu0 : 0 ≤ u) (hut : u ≤ t) (hv0 : 0 ≤ v) (hvt : v ≤ t) (ht1 : t ≤ 1)
    (hL0 : 0 ≤ L) (hL1 : L ≤ 1) (hLe : L ≤ Real.exp (-s)) (hs : 0 ≤ s) {i j : ℕ} (hi : 1 ≤ i)
    (hij : 6 ≤ i + j) :
    (u * L) ^ i * (v * s) ^ j ≤ t ^ 6 * (s ^ j * Real.exp (-s)) := by
  have ht0 : 0 ≤ t := hu0.trans hut
  rw [mul_pow, mul_pow]
  have a1 : u ^ i ≤ t ^ i := pow_le_pow_left₀ hu0 hut i
  have a2 : v ^ j ≤ t ^ j := pow_le_pow_left₀ hv0 hvt j
  have a3 : t ^ i * t ^ j ≤ t ^ 6 := by
    rw [← pow_add]; exact pow_le_pow_of_le_one ht0 ht1 hij
  have a4 : L ^ i ≤ Real.exp (-s) := (pow_le_of_le_one hL0 hL1 (by omega)).trans hLe
  have b1 : u ^ i * v ^ j ≤ t ^ 6 :=
    (mul_le_mul a1 a2 (by positivity) (by positivity)).trans a3
  calc u ^ i * L ^ i * (v ^ j * s ^ j) = (u ^ i * v ^ j) * L ^ i * s ^ j := by ring
    _ ≤ t ^ 6 * Real.exp (-s) * s ^ j := by
        apply mul_le_mul_of_nonneg_right _ (by positivity)
        exact mul_le_mul b1 a4 (by positivity) (by positivity)
    _ = _ := by ring

/-- the degree-`≤ 5` part of `(e^{-uℓ} - 1) e^{-vs}` -/
noncomputable def P5 (u v s : ℝ) : ℝ :=
  (-u) * (s ^ 0 * Y04.ell s ^ 1) + (u * v) * (s ^ 1 * Y04.ell s ^ 1)
    + (-(u * v ^ 2 / 2)) * (s ^ 2 * Y04.ell s ^ 1) + (u * v ^ 3 / 6) * (s ^ 3 * Y04.ell s ^ 1)
    + (-(u * v ^ 4 / 24)) * (s ^ 4 * Y04.ell s ^ 1)
    + (u ^ 2 / 2) * (s ^ 0 * Y04.ell s ^ 2) + (-(u ^ 2 * v / 2)) * (s ^ 1 * Y04.ell s ^ 2)
    + (u ^ 2 * v ^ 2 / 4) * (s ^ 2 * Y04.ell s ^ 2) + (-(u ^ 2 * v ^ 3 / 12)) * (s ^ 3 * Y04.ell s ^ 2)
    + (-(u ^ 3 / 6)) * (s ^ 0 * Y04.ell s ^ 3) + (u ^ 3 * v / 6) * (s ^ 1 * Y04.ell s ^ 3)
    + (-(u ^ 3 * v ^ 2 / 12)) * (s ^ 2 * Y04.ell s ^ 3)
    + (u ^ 4 / 24) * (s ^ 0 * Y04.ell s ^ 4) + (-(u ^ 4 * v / 24)) * (s ^ 1 * Y04.ell s ^ 4)
    + (-(u ^ 5 / 120)) * (s ^ 0 * Y04.ell s ^ 5)

lemma pt5 {u v t s : ℝ} (hu0 : 0 ≤ u) (hut : u ≤ t) (hv0 : 0 ≤ v) (hvt : v ≤ t) (ht1 : t ≤ 1)
    (hs : 0 ≤ s) :
    |(Real.exp (-u * Y04.ell s) - 1) * Real.exp (-v * s) - P5 u v s| ≤ 10 * t ^ 6 * B7 s := by
  have hl0 := Y04.ell_nonneg s
  have hl1 := Y04.ell_le_one s hs
  have hle := Y04.ell_le s
  have ht0 : 0 ≤ t := hu0.trans hut
  unfold P5
  set L := Y04.ell s
  set p := u * L with hp
  set q := v * s with hq
  have hp0 : 0 ≤ p := by positivity
  have hp1 : p ≤ 1 := by
    have := mul_le_mul (hut.trans ht1) hl1 hl0 zero_le_one; rw [hp]; linarith
  have hq0 : 0 ≤ q := by positivity
  set A := -p + p ^ 2 / 2 - p ^ 3 / 6 + p ^ 4 / 24 - p ^ 5 / 120 with hA
  set r1 := Real.exp (-p) - 1 - A with hr1
  set B := 1 - q + q ^ 2 / 2 - q ^ 3 / 6 + q ^ 4 / 24 with hB
  set r2 := Real.exp (-q) - B with hr2
  have hr1b : |r1| ≤ p ^ 6 := ep6_bound hp0 hp1
  have hr2b : |r2| ≤ 4 * q ^ 5 := em5_bound hq0
  set HT := p ^ 2 * q ^ 4 / 48 + p ^ 3 * q ^ 3 / 36 - p ^ 3 * q ^ 4 / 144 + p ^ 4 * q ^ 2 / 48
    - p ^ 4 * q ^ 3 / 144 + p ^ 4 * q ^ 4 / 576 + p ^ 5 * q / 120 - p ^ 5 * q ^ 2 / 240
    + p ^ 5 * q ^ 3 / 720 - p ^ 5 * q ^ 4 / 2880 with hHT
  have e1 : -u * L = -p := by rw [hp]; ring
  have e2 : -v * s = -q := by rw [hq]; ring
  rw [e1, e2]
  have key : (Real.exp (-p) - 1) * Real.exp (-q)
      - ((-u) * (s ^ 0 * L ^ 1) + (u * v) * (s ^ 1 * L ^ 1)
    + (-(u * v ^ 2 / 2)) * (s ^ 2 * L ^ 1) + (u * v ^ 3 / 6) * (s ^ 3 * L ^ 1)
    + (-(u * v ^ 4 / 24)) * (s ^ 4 * L ^ 1)
    + (u ^ 2 / 2) * (s ^ 0 * L ^ 2) + (-(u ^ 2 * v / 2)) * (s ^ 1 * L ^ 2)
    + (u ^ 2 * v ^ 2 / 4) * (s ^ 2 * L ^ 2) + (-(u ^ 2 * v ^ 3 / 12)) * (s ^ 3 * L ^ 2)
    + (-(u ^ 3 / 6)) * (s ^ 0 * L ^ 3) + (u ^ 3 * v / 6) * (s ^ 1 * L ^ 3)
    + (-(u ^ 3 * v ^ 2 / 12)) * (s ^ 2 * L ^ 3)
    + (u ^ 4 / 24) * (s ^ 0 * L ^ 4) + (-(u ^ 4 * v / 24)) * (s ^ 1 * L ^ 4)
    + (-(u ^ 5 / 120)) * (s ^ 0 * L ^ 5))
      = HT + A * r2 + r1 * Real.exp (-q) := by
    rw [hHT, hr1, hr2, hA, hB, hp, hq]; ring
  rw [key]
  have hes : 0 < Real.exp (-s) := Real.exp_pos _
  have m24 := mono_bd6 hu0 hut hv0 hvt ht1 hl0 hl1 hle hs (i := 2) (j := 4) (by norm_num) (by norm_num)
  have m33 := mono_bd6 hu0 hut hv0 hvt ht1 hl0 hl1 hle hs (i := 3) (j := 3) (by norm_num) (by norm_num)
  have m34 := mono_bd6 hu0 hut hv0 hvt ht1 hl0 hl1 hle hs (i := 3) (j := 4) (by norm_num) (by norm_num)
  have m42 := mono_bd6 hu0 hut hv0 hvt ht1 hl0 hl1 hle hs (i := 4) (j := 2) (by norm_num) (by norm_num)
  have m43 := mono_bd6 hu0 hut hv0 hvt ht1 hl0 hl1 hle hs (i := 4) (j := 3) (by norm_num) (by norm_num)
  have m44 := mono_bd6 hu0 hut hv0 hvt ht1 hl0 hl1 hle hs (i := 4) (j := 4) (by norm_num) (by norm_num)
  have m51 := mono_bd6 hu0 hut hv0 hvt ht1 hl0 hl1 hle hs (i := 5) (j := 1) (by norm_num) (by norm_num)
  have m52 := mono_bd6 hu0 hut hv0 hvt ht1 hl0 hl1 hle hs (i := 5) (j := 2) (by norm_num) (by norm_num)
  have m53 := mono_bd6 hu0 hut hv0 hvt ht1 hl0 hl1 hle hs (i := 5) (j := 3) (by norm_num) (by norm_num)
  have m54 := mono_bd6 hu0 hut hv0 hvt ht1 hl0 hl1 hle hs (i := 5) (j := 4) (by norm_num) (by norm_num)
  rw [← hp, ← hq] at m24 m33 m34 m42 m43 m44 m51 m52 m53 m54
  rw [pow_one] at m51
  have n24 : 0 ≤ p ^ 2 * q ^ 4 := by positivity
  have n33 : 0 ≤ p ^ 3 * q ^ 3 := by positivity
  have n34 : 0 ≤ p ^ 3 * q ^ 4 := by positivity
  have n42 : 0 ≤ p ^ 4 * q ^ 2 := by positivity
  have n43 : 0 ≤ p ^ 4 * q ^ 3 := by positivity
  have n44 : 0 ≤ p ^ 4 * q ^ 4 := by positivity
  have n51 : 0 ≤ p ^ 5 * q := by positivity
  have n52 : 0 ≤ p ^ 5 * q ^ 2 := by positivity
  have n53 : 0 ≤ p ^ 5 * q ^ 3 := by positivity
  have n54 : 0 ≤ p ^ 5 * q ^ 4 := by positivity
  have t6 : 0 ≤ t ^ 6 := by positivity
  have bHT : |HT| ≤ t ^ 6 * ((s + s ^ 2 + s ^ 3 + s ^ 4) * Real.exp (-s)) := by
    have w : t ^ 6 * ((s + s ^ 2 + s ^ 3 + s ^ 4) * Real.exp (-s))
        = t ^ 6 * (s ^ 1 * Real.exp (-s)) + t ^ 6 * (s ^ 2 * Real.exp (-s))
          + t ^ 6 * (s ^ 3 * Real.exp (-s)) + t ^ 6 * (s ^ 4 * Real.exp (-s)) := by ring
    have z1 : 0 ≤ t ^ 6 * (s ^ 1 * Real.exp (-s)) := by positivity
    have z2 : 0 ≤ t ^ 6 * (s ^ 2 * Real.exp (-s)) := by positivity
    have z3 : 0 ≤ t ^ 6 * (s ^ 3 * Real.exp (-s)) := by positivity
    have z4 : 0 ≤ t ^ 6 * (s ^ 4 * Real.exp (-s)) := by positivity
    rw [abs_le, hHT]; constructor <;> linarith
  -- A r2
  have hAb : |A| ≤ 2 * p := by
    have p2 : p ^ 2 ≤ p := pow_le_of_le_one hp0 hp1 (by norm_num)
    have p3 : p ^ 3 ≤ p := pow_le_of_le_one hp0 hp1 (by norm_num)
    have p4 : p ^ 4 ≤ p := pow_le_of_le_one hp0 hp1 (by norm_num)
    have p5 : p ^ 5 ≤ p := pow_le_of_le_one hp0 hp1 (by norm_num)
    have : 0 ≤ p ^ 2 := by positivity
    have : 0 ≤ p ^ 3 := by positivity
    have : 0 ≤ p ^ 4 := by positivity
    have : 0 ≤ p ^ 5 := by positivity
    rw [abs_le, hA]; constructor <;> linarith
  have bAr : |A * r2| ≤ 8 * t ^ 6 * (s ^ 5 * Real.exp (-s)) := by
    rw [abs_mul]
    have h1 : |A| * |r2| ≤ 2 * p * (4 * q ^ 5) :=
      mul_le_mul hAb hr2b (abs_nonneg _) (by positivity)
    have h2 : p * q ^ 5 ≤ t ^ 6 * (s ^ 5 * Real.exp (-s)) := by
      have := mono_bd6 hu0 hut hv0 hvt ht1 hl0 hl1 hle hs (i := 1) (j := 5) (by norm_num) (by norm_num)
      rw [← hp, ← hq, pow_one] at this; exact this
    calc |A| * |r2| ≤ 2 * p * (4 * q ^ 5) := h1
      _ = 8 * (p * q ^ 5) := by ring
      _ ≤ 8 * (t ^ 6 * (s ^ 5 * Real.exp (-s))) := by linarith
      _ = _ := by ring
  -- r1 e^{-q}
  have br1 : |r1 * Real.exp (-q)| ≤ t ^ 6 * Real.exp (-s) := by
    rw [abs_mul, abs_of_pos (Real.exp_pos _)]
    have h1 : Real.exp (-q) ≤ 1 := Real.exp_le_one_iff.2 (by linarith)
    have h2 : p ^ 6 ≤ t ^ 6 * Real.exp (-s) := by
      have := mono_bd6 hu0 hut hv0 hvt ht1 hl0 hl1 hle hs (i := 6) (j := 0) (by norm_num) (by norm_num)
      rw [← hp, ← hq, pow_zero, mul_one, pow_zero, one_mul] at this; exact this
    calc |r1| * Real.exp (-q) ≤ p ^ 6 * 1 :=
          mul_le_mul hr1b h1 (Real.exp_pos _).le (by positivity)
      _ ≤ _ := by linarith
  have t1 := abs_add_le (HT + A * r2) (r1 * Real.exp (-q))
  have t2 := abs_add_le HT (A * r2)
  have fin : t ^ 6 * ((s + s ^ 2 + s ^ 3 + s ^ 4) * Real.exp (-s)) + 8 * t ^ 6 * (s ^ 5 * Real.exp (-s))
      + t ^ 6 * Real.exp (-s) ≤ 10 * t ^ 6 * B7 s := by
    unfold B7
    have : 10 * t ^ 6 * ((1 + s + s ^ 2 + s ^ 3 + s ^ 4 + s ^ 5) * Real.exp (-s))
        - (t ^ 6 * ((s + s ^ 2 + s ^ 3 + s ^ 4) * Real.exp (-s)) + 8 * t ^ 6 * (s ^ 5 * Real.exp (-s))
          + t ^ 6 * Real.exp (-s))
        = t ^ 6 * Real.exp (-s) * (9 + 9 * s + 9 * s ^ 2 + 9 * s ^ 3 + 9 * s ^ 4 + 2 * s ^ 5) := by ring
    have : 0 ≤ t ^ 6 * Real.exp (-s) * (9 + 9 * s + 9 * s ^ 2 + 9 * s ^ 3 + 9 * s ^ 4 + 2 * s ^ 5) := by
      positivity
    linarith
  linarith

/-- five coefficients vanish if the quintic is `O(t⁶)` on `(0, d]` -/
lemma quintic_zero {d1 d2 d3 d4 d5 K d : ℝ} (hd : 0 < d) (hd1 : d ≤ 1) (hK : 0 ≤ K)
    (h : ∀ t : ℝ, 0 < t → t ≤ d → |d1 * t + d2 * t ^ 2 + d3 * t ^ 3 + d4 * t ^ 4 + d5 * t ^ 5| ≤ K * t ^ 6) :
    d1 = 0 ∧ d2 = 0 ∧ d3 = 0 ∧ d4 = 0 ∧ d5 = 0 := by
  have h1 : d1 = 0 := by
    refine Y04.eq_zero_of_le_lin hd (K := |d2| + |d3| + |d4| + |d5| + K)
      (by have := abs_nonneg d2; have := abs_nonneg d3; have := abs_nonneg d4; have := abs_nonneg d5
          linarith)
      (fun t ht0 htd => ?_)
    have ht1 : t ≤ 1 := htd.trans hd1
    have := h t ht0 htd
    set E := d1 * t + d2 * t ^ 2 + d3 * t ^ 3 + d4 * t ^ 4 + d5 * t ^ 5 with hE
    have e : d1 = E / t - d2 * t - d3 * t ^ 2 - d4 * t ^ 3 - d5 * t ^ 4 := by
      rw [hE]; field_simp; ring
    rw [e]
    have ht2 : t ^ 2 ≤ t := by nlinarith
    have ht3 : t ^ 3 ≤ t := by nlinarith [pow_le_one₀ ht0.le ht1 (n := 2)]
    have ht4 : t ^ 4 ≤ t := by nlinarith [pow_le_one₀ ht0.le ht1 (n := 3)]
    have ht5 : t ^ 5 ≤ t := by nlinarith [pow_le_one₀ ht0.le ht1 (n := 4)]
    have b1 : |E / t| ≤ K * t ^ 5 := by
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
    have t5 : K * t ^ 5 ≤ K * t := mul_le_mul_of_nonneg_left ht5 hK
    have := abs_sub (E / t - d2 * t - d3 * t ^ 2 - d4 * t ^ 3) (d5 * t ^ 4)
    have := abs_sub (E / t - d2 * t - d3 * t ^ 2) (d4 * t ^ 3)
    have := abs_sub (E / t - d2 * t) (d3 * t ^ 2)
    have := abs_sub (E / t) (d2 * t)
    have : (|d2| + |d3| + |d4| + |d5| + K) * t = |d2| * t + |d3| * t + |d4| * t + |d5| * t + K * t := by
      ring
    linarith
  subst h1
  obtain ⟨a, b, c, e⟩ := Z02.quartic_zero hd hd1 hK (d1 := d2) (d2 := d3) (d3 := d4) (d4 := d5)
    (fun t ht0 htd => by
    have := h t ht0 htd
    have e : 0 * t + d2 * t ^ 2 + d3 * t ^ 3 + d4 * t ^ 4 + d5 * t ^ 5
        = t * (d2 * t + d3 * t ^ 2 + d4 * t ^ 3 + d5 * t ^ 4) := by ring
    rw [e, abs_mul, abs_of_pos ht0, show K * t ^ 6 = t * (K * t ^ 5) by ring] at this
    exact le_of_mul_le_mul_left this ht0)
  exact ⟨rfl, a, b, c, e⟩


/-- the integrated expansion of `(e^{-uℓ} - 1) e^{-vs}` to order `5` -/
lemma int_expand5 {u v t : ℝ} (hu0 : 0 ≤ u) (hut : u ≤ t) (hv0 : 0 ≤ v) (hvt : v ≤ t) (ht1 : t ≤ 1) :
    IntegrableOn (fun s => (Real.exp (-u * Y04.ell s) - 1) * Real.exp (-v * s)) (Set.Ioi 0) ∧
    |(∫ s in Set.Ioi (0 : ℝ), (Real.exp (-u * Y04.ell s) - 1) * Real.exp (-v * s))
      - ((-u) * Y04.I 0 1 + (u * v) * Y04.I 1 1 + (-(u * v ^ 2 / 2)) * Y04.I 2 1
        + (u * v ^ 3 / 6) * Y04.I 3 1 + (-(u * v ^ 4 / 24)) * Y04.I 4 1
        + (u ^ 2 / 2) * Y04.I 0 2 + (-(u ^ 2 * v / 2)) * Y04.I 1 2
        + (u ^ 2 * v ^ 2 / 4) * Y04.I 2 2 + (-(u ^ 2 * v ^ 3 / 12)) * Y04.I 3 2
        + (-(u ^ 3 / 6)) * Y04.I 0 3 + (u ^ 3 * v / 6) * Y04.I 1 3 + (-(u ^ 3 * v ^ 2 / 12)) * Y04.I 2 3
        + (u ^ 4 / 24) * Y04.I 0 4 + (-(u ^ 4 * v / 24)) * Y04.I 1 4 + (-(u ^ 5 / 120)) * Y04.I 0 5)|
      ≤ 10 * t ^ 6 * ∫ s in Set.Ioi (0 : ℝ), B7 s := by
  have hint : IntegrableOn (fun s => (Real.exp (-u * Y04.ell s) - 1) * Real.exp (-v * s)) (Set.Ioi 0) :=
    (Z02.int_expand hu0 hut hv0 hvt ht1).1
  refine ⟨hint, ?_⟩
  obtain ⟨hmi, hmv⟩ := int15 (-u) (u * v) (-(u * v ^ 2 / 2)) (u * v ^ 3 / 6) (-(u * v ^ 4 / 24))
    (u ^ 2 / 2) (-(u ^ 2 * v / 2)) (u ^ 2 * v ^ 2 / 4) (-(u ^ 2 * v ^ 3 / 12))
    (-(u ^ 3 / 6)) (u ^ 3 * v / 6) (-(u ^ 3 * v ^ 2 / 12))
    (u ^ 4 / 24) (-(u ^ 4 * v / 24)) (-(u ^ 5 / 120))
    (Y04.I_int 0 1 le_rfl) (Y04.I_int 1 1 le_rfl) (Y04.I_int 2 1 le_rfl) (Y04.I_int 3 1 le_rfl)
    (Y04.I_int 4 1 le_rfl)
    (Y04.I_int 0 2 (by norm_num)) (Y04.I_int 1 2 (by norm_num)) (Y04.I_int 2 2 (by norm_num))
    (Y04.I_int 3 2 (by norm_num))
    (Y04.I_int 0 3 (by norm_num)) (Y04.I_int 1 3 (by norm_num)) (Y04.I_int 2 3 (by norm_num))
    (Y04.I_int 0 4 (by norm_num)) (Y04.I_int 1 4 (by norm_num)) (Y04.I_int 0 5 (by norm_num))
  unfold Y04.I
  rw [← hmv, ← integral_sub hint hmi]
  have hB := (B7_int.const_mul 10).const_mul (t ^ 6)
  have := norm_integral_le_of_norm_le hB
    ((ae_restrict_iff' measurableSet_Ioi).2 (Filter.Eventually.of_forall
      (fun s (hs : 0 < s) => by
        rw [Real.norm_eq_abs]
        have := pt5 hu0 hut hv0 hvt ht1 hs.le
        unfold P5 at this
        calc _ = _ := rfl
          _ ≤ _ := this
          _ = t ^ 6 * (10 * B7 s) := by ring)))
  rw [Real.norm_eq_abs, integral_const_mul, integral_const_mul] at this
  calc _ ≤ _ := this
    _ = _ := by ring


lemma series_trunc6 {α β γ x L : ℝ} (hα : |α| ≤ 3 / 2) (hβ : |β| ≤ 3 / 2) (hγ : |γ| ≤ 3 / 2)
    (hx0 : 0 ≤ x) (hx : x ≤ 1 / 3)
    (hL : HasSum (fun k : ℕ => (-1 : ℝ) ^ (k + 2) * Y04.zn (k + 2) / ((k + 2 : ℕ) : ℝ) *
      ((α * x) ^ (k + 2) + (β * x) ^ (k + 2) - (γ * x) ^ (k + 2))) L) :
    |L - (Y04.cf α β γ 0 * x ^ 2 + Y04.cf α β γ 1 * x ^ 3 + Y04.cf α β γ 2 * x ^ 4
      + Y04.cf α β γ 3 * x ^ 5 + Y04.cf α β γ 4 * x ^ 6)| ≤ 103 * x ^ 7 := by
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
  have hg := (hasSum_nat_add_iff' 5).mpr hL
  have hgeo : HasSum (fun k : ℕ => 3 * r ^ 7 * r ^ k) (3 * r ^ 7 * (1 - r)⁻¹) :=
    (hasSum_geometric_of_lt_one hr0 (by linarith)).mul_left _
  have hb1 : ∀ k : ℕ, Y04.cf α β γ (k + 5) * x ^ (k + 5 + 2) ≤ 3 * r ^ 7 * r ^ k := by
    intro k
    have := (abs_le.mp (hbd (k + 5))).2
    calc _ ≤ 3 * r ^ (k + 5 + 2) := this
      _ = _ := by ring
  have hb2 : ∀ k : ℕ, -(Y04.cf α β γ (k + 5) * x ^ (k + 5 + 2)) ≤ 3 * r ^ 7 * r ^ k := by
    intro k
    have := (abs_le.mp (hbd (k + 5))).1
    calc _ ≤ 3 * r ^ (k + 5 + 2) := by linarith
      _ = _ := by ring
  have u1 := hasSum_le hb1 hg hgeo
  have u2 := hasSum_le hb2 hg.neg hgeo
  have hinv : 3 * r ^ 7 * (1 - r)⁻¹ ≤ 6 * r ^ 7 := by
    rw [mul_inv_le_iff₀ (by linarith)]
    nlinarith [pow_nonneg hr0 7]
  have hr7 : 6 * r ^ 7 ≤ 103 * x ^ 7 := by
    rw [hr, mul_pow]
    nlinarith [pow_nonneg hx0 7]
  simp only [Finset.sum_range_succ, Finset.sum_range_zero, zero_add] at u1 u2
  rw [abs_le]
  constructor <;> norm_num at u1 u2 ⊢ <;> linarith

lemma abs_mul_le {a b A B : ℝ} (ha : |a| ≤ A) (hb : |b| ≤ B) : |a * b| ≤ A * B := by
  rw [abs_mul]; exact mul_le_mul ha hb (abs_nonneg _) ((abs_nonneg _).trans ha)

lemma exp_trunc6 {L a0 a1 a2 a3 a4 x B K : ℝ} (hB : 0 ≤ B) (_hK : 0 ≤ K)
    (ha0 : |a0| ≤ B) (ha1 : |a1| ≤ B) (ha2 : |a2| ≤ B) (ha3 : |a3| ≤ B) (ha4 : |a4| ≤ B)
    (hx0 : 0 ≤ x) (hx1 : x ≤ 1) (hs1 : 5 * B * x ^ 2 ≤ 1) (hs2 : K * x ^ 7 ≤ 1)
    (hL : |L - (a0 * x ^ 2 + a1 * x ^ 3 + a2 * x ^ 4 + a3 * x ^ 5 + a4 * x ^ 6)| ≤ K * x ^ 7) :
    |Real.exp L - (1 + a0 * x ^ 2 + a1 * x ^ 3 + (a2 + a0 ^ 2 / 2) * x ^ 4
        + (a3 + a0 * a1) * x ^ 5 + (a4 + a0 * a2 + a1 ^ 2 / 2 + a0 ^ 3 / 6) * x ^ 6)|
      ≤ (6 * K + (5 * B) ^ 4 + 10 * B ^ 2 + 21 * B ^ 3) * x ^ 7 := by
  set P := a0 * x ^ 2 + a1 * x ^ 3 + a2 * x ^ 4 + a3 * x ^ 5 + a4 * x ^ 6 with hPdef
  set Pd := a2 + a3 * x + a4 * x ^ 2 with hPd
  set Pb := a1 + x * Pd with hPb
  set Pa := a0 + x * Pb with hPa
  have hx2 : x ^ 2 ≤ 1 := pow_le_one₀ hx0 hx1
  have hP : P = x ^ 2 * Pa := by rw [hPdef, hPa, hPb, hPd]; ring
  have bPd : |Pd| ≤ 3 * B := by
    have t1 := abs_add_le (a2 + a3 * x) (a4 * x ^ 2)
    have t2 := abs_add_le a2 (a3 * x)
    have e1 : |a3 * x| ≤ B * 1 := abs_mul_le ha3 (by rw [abs_of_nonneg hx0]; exact hx1)
    have e2 : |a4 * x ^ 2| ≤ B * 1 := abs_mul_le ha4 (by rw [abs_of_nonneg (by positivity)]; exact hx2)
    rw [hPd]; linarith
  have bPb : |Pb| ≤ 4 * B := by
    have t1 := abs_add_le a1 (x * Pd)
    have e1 : |x * Pd| ≤ 1 * (3 * B) := abs_mul_le (by rw [abs_of_nonneg hx0]; exact hx1) bPd
    rw [hPb]; linarith
  have bPa : |Pa| ≤ 5 * B := by
    have t1 := abs_add_le a0 (x * Pb)
    have e1 : |x * Pb| ≤ 1 * (4 * B) := abs_mul_le (by rw [abs_of_nonneg hx0]; exact hx1) bPb
    rw [hPa]; linarith
  have bP : |P| ≤ 5 * B * x ^ 2 := by
    rw [hP, abs_mul, abs_of_nonneg (by positivity)]
    calc x ^ 2 * |Pa| ≤ x ^ 2 * (5 * B) := mul_le_mul_of_nonneg_left bPa (by positivity)
      _ = _ := by ring
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
  have hP4 : |P| ^ 4 ≤ (5 * B) ^ 4 * x ^ 7 := by
    calc |P| ^ 4 ≤ (5 * B * x ^ 2) ^ 4 := pow_le_pow_left₀ (abs_nonneg _) bP 4
      _ = (5 * B) ^ 4 * x ^ 7 * x := by ring
      _ ≤ (5 * B) ^ 4 * x ^ 7 * 1 := mul_le_mul_of_nonneg_left hx1 (by positivity)
      _ = _ := by ring
  set W := a0 * (a3 + a4 * x) + (2 * a1 * Pd + x * Pd ^ 2) / 2
    + (3 * a0 ^ 2 * Pb + 3 * a0 * x * Pb ^ 2 + x ^ 2 * Pb ^ 3) / 6 with hW
  have hWb : |W| ≤ 10 * B ^ 2 + 21 * B ^ 3 := by
    have ex : |x| ≤ 1 := by rw [abs_of_nonneg hx0]; exact hx1
    have ex2 : |x ^ 2| ≤ 1 := by rw [abs_of_nonneg (by positivity)]; exact hx2
    have q1 : |a3 + a4 * x| ≤ 2 * B := by
      have := abs_add_le a3 (a4 * x)
      have := abs_mul_le ha4 ex
      linarith
    have w1 := abs_mul_le ha0 q1
    have h2a1 : |(2:ℝ) * a1| ≤ 2 * B := by rw [abs_mul]; norm_num; linarith
    have w2 : |2 * a1 * Pd| ≤ (2 * B) * (3 * B) := abs_mul_le h2a1 bPd
    have hPd2 : |Pd ^ 2| ≤ (3 * B) ^ 2 := by rw [abs_pow]; exact pow_le_pow_left₀ (abs_nonneg _) bPd 2
    have hPb2 : |Pb ^ 2| ≤ (4 * B) ^ 2 := by rw [abs_pow]; exact pow_le_pow_left₀ (abs_nonneg _) bPb 2
    have hPb3 : |Pb ^ 3| ≤ (4 * B) ^ 3 := by rw [abs_pow]; exact pow_le_pow_left₀ (abs_nonneg _) bPb 3
    have ha02 : |3 * a0 ^ 2| ≤ 3 * B ^ 2 := by
      rw [abs_mul, abs_pow, abs_of_pos (by norm_num : (0:ℝ) < 3)]
      exact mul_le_mul_of_nonneg_left (pow_le_pow_left₀ (abs_nonneg a0) ha0 2) (by norm_num)
    have w3 := abs_mul_le ex hPd2
    have w4 := abs_mul_le ha02 bPb
    have w5 := abs_mul_le (abs_mul_le (show |(3:ℝ) * a0| ≤ 3 * B by
      rw [abs_mul]; norm_num; linarith) ex) hPb2
    have w6 := abs_mul_le ex2 hPb3
    have t1 := abs_add_le (a0 * (a3 + a4 * x) + (2 * a1 * Pd + x * Pd ^ 2) / 2)
      ((3 * a0 ^ 2 * Pb + 3 * a0 * x * Pb ^ 2 + x ^ 2 * Pb ^ 3) / 6)
    have t2 := abs_add_le (a0 * (a3 + a4 * x)) ((2 * a1 * Pd + x * Pd ^ 2) / 2)
    have t3 : |(2 * a1 * Pd + x * Pd ^ 2) / 2| ≤ (|2 * a1 * Pd| + |x * Pd ^ 2|) / 2 := by
      rw [abs_div, abs_of_pos (by norm_num : (0:ℝ) < 2)]
      exact div_le_div_of_nonneg_right (abs_add_le _ _) (by norm_num)
    have t4 : |(3 * a0 ^ 2 * Pb + 3 * a0 * x * Pb ^ 2 + x ^ 2 * Pb ^ 3) / 6|
        ≤ (|3 * a0 ^ 2 * Pb| + |3 * a0 * x * Pb ^ 2| + |x ^ 2 * Pb ^ 3|) / 6 := by
      rw [abs_div, abs_of_pos (by norm_num : (0:ℝ) < 6)]
      apply div_le_div_of_nonneg_right _ (by norm_num)
      have := abs_add_le (3 * a0 ^ 2 * Pb + 3 * a0 * x * Pb ^ 2) (x ^ 2 * Pb ^ 3)
      have := abs_add_le (3 * a0 ^ 2 * Pb) (3 * a0 * x * Pb ^ 2)
      linarith
    rw [hW]
    have hB2 : 0 ≤ B ^ 2 := by positivity
    have hB3 : 0 ≤ B ^ 3 := by positivity
    linarith
  have key : Real.exp L - (1 + a0 * x ^ 2 + a1 * x ^ 3 + (a2 + a0 ^ 2 / 2) * x ^ 4
        + (a3 + a0 * a1) * x ^ 5 + (a4 + a0 * a2 + a1 ^ 2 / 2 + a0 ^ 3 / 6) * x ^ 6)
      = Real.exp P * (Real.exp E - 1) + (Real.exp P - (1 + P + P ^ 2 / 2 + P ^ 3 / 6)) + x ^ 7 * W := by
    rw [hLPE, Real.exp_add, hW, hPdef, hPb, hPd]; ring
  rw [key]
  have b1 : |Real.exp P * (Real.exp E - 1)| ≤ 3 * (2 * (K * x ^ 7)) := by
    rw [abs_mul, abs_of_pos (Real.exp_pos _)]
    exact mul_le_mul hexpP (hexpE.trans (by linarith)) (abs_nonneg _) (by norm_num)
  have b2 : |Real.exp P - (1 + P + P ^ 2 / 2 + P ^ 3 / 6)| ≤ (5 * B) ^ 4 * x ^ 7 := by
    calc _ ≤ _ := heP
      _ ≤ |P| ^ 4 * 1 := by gcongr; norm_num [Nat.factorial]
      _ ≤ _ := by rw [mul_one]; exact hP4
  have b3 : |x ^ 7 * W| ≤ x ^ 7 * (10 * B ^ 2 + 21 * B ^ 3) := by
    rw [abs_mul, abs_of_nonneg (by positivity)]
    exact mul_le_mul_of_nonneg_left hWb (by positivity)
  have t1 := abs_add_le (Real.exp P * (Real.exp E - 1) + (Real.exp P - (1 + P + P ^ 2 / 2 + P ^ 3 / 6)))
    (x ^ 7 * W)
  have t2 := abs_add_le (Real.exp P * (Real.exp E - 1)) (Real.exp P - (1 + P + P ^ 2 / 2 + P ^ 3 / 6))
  linarith only [t1, t2, b1, b2, b3]


noncomputable def K6 : ℝ := 6 * 103 + (5 * 35) ^ 4 + 10 * 35 ^ 2 + 21 * 35 ^ 3

lemma A_expand6 {α β γ x L : ℝ} (hα : |α| ≤ 3 / 2) (hβ : |β| ≤ 3 / 2) (hγ : |γ| ≤ 3 / 2)
    (hx0 : 0 ≤ x) (hx : x ≤ 1 / 20)
    (hL : HasSum (fun k : ℕ => (-1 : ℝ) ^ (k + 2) * Y04.zn (k + 2) / ((k + 2 : ℕ) : ℝ) *
      ((α * x) ^ (k + 2) + (β * x) ^ (k + 2) - (γ * x) ^ (k + 2))) L) :
    |Real.exp L - (1 + Y04.cf α β γ 0 * x ^ 2 + Y04.cf α β γ 1 * x ^ 3
        + (Y04.cf α β γ 2 + Y04.cf α β γ 0 ^ 2 / 2) * x ^ 4
        + (Y04.cf α β γ 3 + Y04.cf α β γ 0 * Y04.cf α β γ 1) * x ^ 5
        + (Y04.cf α β γ 4 + Y04.cf α β γ 0 * Y04.cf α β γ 2 + Y04.cf α β γ 1 ^ 2 / 2
          + Y04.cf α β γ 0 ^ 3 / 6) * x ^ 6)| ≤ K6 * x ^ 7 := by
  have h := series_trunc6 hα hβ hγ hx0 (by linarith) hL
  have hB : ∀ k ≤ 4, |Y04.cf α β γ k| ≤ 35 := by
    intro k hk
    have := Y04.cf_bd hα hβ hγ k
    have : (3 / 2 : ℝ) ^ (k + 2) ≤ (3 / 2) ^ 6 := pow_le_pow_right₀ (by norm_num) (by omega)
    norm_num at this ⊢; linarith
  have hx2 : x ^ 2 ≤ 1 / 400 := by
    have := pow_le_pow_left₀ hx0 hx 2
    norm_num at this ⊢; linarith
  have hx7 : x ^ 7 ≤ 1 / 400 := by
    have : x ^ 7 ≤ x ^ 2 := pow_le_pow_of_le_one hx0 (by linarith) (by norm_num)
    linarith
  unfold K6
  exact exp_trunc6 (by norm_num) (by norm_num) (hB 0 (by norm_num)) (hB 1 (by norm_num))
    (hB 2 (by norm_num)) (hB 3 (by norm_num)) (hB 4 (by norm_num)) hx0 (by linarith)
    (by nlinarith) (by nlinarith) h

lemma K6_pos : 0 ≤ K6 := by unfold K6; norm_num

/-- the `x⁵` coefficient of the regularised Beta identity at `a = αx`, `b = βx` -/
theorem coef5 {α β : ℝ} (hα : 0 < α) (hβ : 0 < β) (h3 : α + β ≤ 3) :
    (α + β) * (α ^ 4 + β ^ 4) / 24 * Y04.I 4 1 + (α + β) ^ 2 * (α ^ 3 + β ^ 3) / 12 * Y04.I 3 2
      + (α + β) ^ 3 * (α ^ 2 + β ^ 2) / 12 * Y04.I 2 3 + (α + β) ^ 4 * (α + β) / 24 * Y04.I 1 4
      + (α + β) ^ 5 / 60 * Y04.I 0 5
      = -((α + β) / (α * β) * 64 * (Y04.cf (α / 2) (β / 2) ((α + β) / 2) 4
          + Y04.cf (α / 2) (β / 2) ((α + β) / 2) 0 * Y04.cf (α / 2) (β / 2) ((α + β) / 2) 2
          + Y04.cf (α / 2) (β / 2) ((α + β) / 2) 1 ^ 2 / 2
          + Y04.cf (α / 2) (β / 2) ((α + β) / 2) 0 ^ 3 / 6)) := by
  set γ := α + β with hγ
  set c0 := Y04.cf (α / 2) (β / 2) (γ / 2) 0 with hc0
  set c1 := Y04.cf (α / 2) (β / 2) (γ / 2) 1 with hc1
  set c2 := Y04.cf (α / 2) (β / 2) (γ / 2) 2 with hc2
  set c3 := Y04.cf (α / 2) (β / 2) (γ / 2) 3 with hc3
  set c4 := Y04.cf (α / 2) (β / 2) (γ / 2) 4 with hc4
  set KB := ∫ s in Set.Ioi (0 : ℝ), B7 s with hKB
  have hKB0 : 0 ≤ KB := setIntegral_nonneg measurableSet_Ioi (fun s hs => B7_nonneg (le_of_lt hs))
  have hK0p := K6_pos
  have hα3 : α ≤ 3 := by linarith
  have hβ3 : β ≤ 3 := by linarith
  have hγ0 : 0 < γ := by rw [hγ]; linarith
  set e := γ / (α * β) with he
  have he0 : 0 ≤ e := by positivity
  set d1 := -(2 * γ) * Y04.I 0 1 - e * 4 * c0 with hd1
  set d2 := γ * (α + β) * Y04.I 1 1 + γ ^ 2 * Y04.I 0 2 - e * 8 * c1 with hd2
  set d3 := -(γ * (α ^ 2 + β ^ 2) / 2) * Y04.I 2 1 - γ ^ 2 * (α + β) / 2 * Y04.I 1 2
    - γ ^ 3 / 3 * Y04.I 0 3 - e * 16 * (c2 + c0 ^ 2 / 2) with hd3
  set d4 := γ * (α ^ 3 + β ^ 3) / 6 * Y04.I 3 1 + γ ^ 2 * (α ^ 2 + β ^ 2) / 4 * Y04.I 2 2
      + γ ^ 3 * (α + β) / 6 * Y04.I 1 3 + γ ^ 4 / 12 * Y04.I 0 4 - e * 32 * (c3 + c0 * c1) with hd4
  set d5 := -(γ * (α ^ 4 + β ^ 4) / 24) * Y04.I 4 1 - γ ^ 2 * (α ^ 3 + β ^ 3) / 12 * Y04.I 3 2
      - γ ^ 3 * (α ^ 2 + β ^ 2) / 12 * Y04.I 2 3 - γ ^ 4 * (α + β) / 24 * Y04.I 1 4
      - γ ^ 5 / 60 * Y04.I 0 5 - e * 64 * (c4 + c0 * c2 + c1 ^ 2 / 2 + c0 ^ 3 / 6) with hd5
  have hmain : ∀ x : ℝ, 0 < x → x ≤ 1 / 40 →
      |d1 * x + d2 * x ^ 2 + d3 * x ^ 3 + d4 * x ^ 4 + d5 * x ^ 5|
        ≤ (20 * 729 * KB + e * 128 * K6) * x ^ 6 := by
    intro x hx0 hx1
    have hx3 : 3 * x ≤ 1 := by linarith
    have gx : γ * x ≤ 3 * x := mul_le_mul_of_nonneg_right (by rw [hγ]; exact h3) hx0.le
    have ax : α * x ≤ 3 * x := mul_le_mul_of_nonneg_right hα3 hx0.le
    have bx : β * x ≤ 3 * x := mul_le_mul_of_nonneg_right hβ3 hx0.le
    have dx1 : (α - β) * x ≤ 3 * x := mul_le_mul_of_nonneg_right (by linarith) hx0.le
    have dx2 : -3 * x ≤ (α - β) * x := mul_le_mul_of_nonneg_right (by linarith) hx0.le
    obtain ⟨_, E1⟩ := int_expand5 (u := γ * x) (v := α * x) (t := 3 * x) (by positivity)
      gx (by positivity) ax hx3
    obtain ⟨_, E2⟩ := int_expand5 (u := γ * x) (v := β * x) (t := 3 * x) (by positivity)
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
    have hA := A_expand6 (α := α / 2) (β := β / 2) (γ := γ / 2)
      (by rw [abs_of_pos (by positivity)]; linarith) (by rw [abs_of_pos (by positivity)]; linarith)
      (by rw [abs_of_pos (by positivity)]; linarith) (by positivity) (by linarith) h
    rw [Real.exp_log hpos] at hA
    rw [show α / 2 * (2 * x) = α * x by ring, show β / 2 * (2 * x) = β * x by ring,
      show γ / 2 * (2 * x) = γ * x by ring] at hA
    rw [← hc0, ← hc1, ← hc2, ← hc3, ← hc4] at hA
    set Rg := Real.Gamma (1 + α * x) * Real.Gamma (1 + β * x) / Real.Gamma (1 + γ * x) with hRg
    set Iα := ∫ s in Set.Ioi (0 : ℝ), (Real.exp (-(γ * x) * Y04.ell s) - 1) * Real.exp (-(α * x) * s)
    set Iβ := ∫ s in Set.Ioi (0 : ℝ), (Real.exp (-(γ * x) * Y04.ell s) - 1) * Real.exp (-(β * x) * s)
    set u := γ * x with hu
    set Sα := (-u) * Y04.I 0 1 + (u * (α * x)) * Y04.I 1 1 + (-(u * (α * x) ^ 2 / 2)) * Y04.I 2 1
        + (u * (α * x) ^ 3 / 6) * Y04.I 3 1 + (-(u * (α * x) ^ 4 / 24)) * Y04.I 4 1
        + (u ^ 2 / 2) * Y04.I 0 2 + (-(u ^ 2 * (α * x) / 2)) * Y04.I 1 2
        + (u ^ 2 * (α * x) ^ 2 / 4) * Y04.I 2 2 + (-(u ^ 2 * (α * x) ^ 3 / 12)) * Y04.I 3 2
        + (-(u ^ 3 / 6)) * Y04.I 0 3 + (u ^ 3 * (α * x) / 6) * Y04.I 1 3
        + (-(u ^ 3 * (α * x) ^ 2 / 12)) * Y04.I 2 3
        + (u ^ 4 / 24) * Y04.I 0 4 + (-(u ^ 4 * (α * x) / 24)) * Y04.I 1 4
        + (-(u ^ 5 / 120)) * Y04.I 0 5 with hSα
    set Sβ := (-u) * Y04.I 0 1 + (u * (β * x)) * Y04.I 1 1 + (-(u * (β * x) ^ 2 / 2)) * Y04.I 2 1
        + (u * (β * x) ^ 3 / 6) * Y04.I 3 1 + (-(u * (β * x) ^ 4 / 24)) * Y04.I 4 1
        + (u ^ 2 / 2) * Y04.I 0 2 + (-(u ^ 2 * (β * x) / 2)) * Y04.I 1 2
        + (u ^ 2 * (β * x) ^ 2 / 4) * Y04.I 2 2 + (-(u ^ 2 * (β * x) ^ 3 / 12)) * Y04.I 3 2
        + (-(u ^ 3 / 6)) * Y04.I 0 3 + (u ^ 3 * (β * x) / 6) * Y04.I 1 3
        + (-(u ^ 3 * (β * x) ^ 2 / 12)) * Y04.I 2 3
        + (u ^ 4 / 24) * Y04.I 0 4 + (-(u ^ 4 * (β * x) / 24)) * Y04.I 1 4
        + (-(u ^ 5 / 120)) * Y04.I 0 5 with hSβ
    set Q := 1 + c0 * (2 * x) ^ 2 + c1 * (2 * x) ^ 3 + (c2 + c0 ^ 2 / 2) * (2 * x) ^ 4
        + (c3 + c0 * c1) * (2 * x) ^ 5 + (c4 + c0 * c2 + c1 ^ 2 / 2 + c0 ^ 3 / 6) * (2 * x) ^ 6 with hQ
    have hpoly : d1 * x + d2 * x ^ 2 + d3 * x ^ 3 + d4 * x ^ 4 + d5 * x ^ 5
        = (Sα - Iα) + (Sβ - Iβ) + ((Iα + Iβ) - e * (4 * c0 * x + 8 * c1 * x ^ 2
          + 16 * (c2 + c0 ^ 2 / 2) * x ^ 3 + 32 * (c3 + c0 * c1) * x ^ 4
          + 64 * (c4 + c0 * c2 + c1 ^ 2 / 2 + c0 ^ 3 / 6) * x ^ 5)) := by
      rw [hd1, hd2, hd3, hd4, hd5, hSα, hSβ, hu]; ring
    have hlast : (Iα + Iβ) - e * (4 * c0 * x + 8 * c1 * x ^ 2
          + 16 * (c2 + c0 ^ 2 / 2) * x ^ 3 + 32 * (c3 + c0 * c1) * x ^ 4
          + 64 * (c4 + c0 * c2 + c1 ^ 2 / 2 + c0 ^ 3 / 6) * x ^ 5)
        = e / x * (Rg - Q) := by
      rw [hC, he, hQ, hu]
      field_simp
      ring
    rw [hpoly, hlast]
    have b3 : |e / x * (Rg - Q)| ≤ e * 128 * K6 * x ^ 6 := by
      rw [abs_mul, abs_of_nonneg (by positivity)]
      calc e / x * |Rg - Q| ≤ e / x * (K6 * (2 * x) ^ 7) :=
            mul_le_mul_of_nonneg_left hA (by positivity)
        _ = e * 128 * K6 * x ^ 6 := by field_simp; ring
    have b1 : |Sα - Iα| ≤ 10 * (3 * x) ^ 6 * KB := by rw [abs_sub_comm]; exact E1
    have b2 : |Sβ - Iβ| ≤ 10 * (3 * x) ^ 6 * KB := by rw [abs_sub_comm]; exact E2
    have t1 := abs_add_le ((Sα - Iα) + (Sβ - Iβ)) (e / x * (Rg - Q))
    have t2 := abs_add_le (Sα - Iα) (Sβ - Iβ)
    have : (20 * 729 * KB + e * 128 * K6) * x ^ 6
        = 10 * (3 * x) ^ 6 * KB + 10 * (3 * x) ^ 6 * KB + e * 128 * K6 * x ^ 6 := by ring
    linarith
  obtain ⟨_, _, _, _, h5⟩ := quintic_zero (by norm_num) (by norm_num) (by positivity) hmain
  rw [hd5] at h5
  linarith


lemma eqP1 : 1 / 6 * Y04.I 4 1 + 2 / 3 * Y04.I 3 2 + 4 / 3 * Y04.I 2 3 + 4 / 3 * Y04.I 1 4
    + 8 / 15 * Y04.I 0 5 = 79 / 8 * (π ^ 6 / 945) - 4 * Y04.z3 ^ 2 := by
  have h := coef5 (α := 1) (β := 1) one_pos one_pos (by norm_num)
  norm_num [Y04.cf, Y04.z_two, Y04.zn_three, Y04.z_four, Z02.zn_five, z_six] at h
  linear_combination h

lemma eqP2 : 17 / 8 * Y04.I 4 1 + 27 / 4 * Y04.I 3 2 + 45 / 4 * Y04.I 2 3 + 81 / 8 * Y04.I 1 4
    + 81 / 20 * Y04.I 0 5 = 363 / 4 * (π ^ 6 / 945) - 27 * Y04.z3 ^ 2 := by
  have h := coef5 (α := 1) (β := 2) one_pos two_pos (by norm_num)
  norm_num [Y04.cf, Y04.z_two, Y04.zn_three, Y04.z_four, Z02.zn_five, z_six] at h
  linear_combination h

lemma eqP3 : 41 / 729 * Y04.I 4 1 + 112 / 729 * Y04.I 3 2 + 160 / 729 * Y04.I 2 3
    + 128 / 729 * Y04.I 1 4 + 256 / 3645 * Y04.I 0 5
    = 1879 / 972 * (π ^ 6 / 945) - 32 / 81 * Y04.z3 ^ 2 := by
  have h := coef5 (α := 1) (β := 1 / 3) one_pos (by norm_num) (by norm_num)
  norm_num [Y04.cf, Y04.z_two, Y04.zn_three, Y04.z_four, Z02.zn_five, z_six] at h
  linear_combination h

lemma D61 : 2 * (-(23 / 24) * Y04.I 4 1 + (1 / 2) * Y04.I 3 2 / 2 + Y04.I 2 3 / 6) =
    -(733 / 16) * (π ^ 6 / 945) + Y04.z3 ^ 2 := by
  linear_combination (-1885 / 4) * eqP1 + (1672 / 15) * eqP2 + (-113967 / 40) * eqP3

lemma D60 : 2 * (121 / 384 * Y04.I 4 1 + (-(17 / 48)) * Y04.I 3 2 / 2 + (-(5 / 8)) * Y04.I 2 3 / 6
      + (1 / 2) * Y04.I 1 4 / 24 + Y04.I 0 5 / 120) =
    4039 / 256 * (π ^ 6 / 945) - 7 / 8 * Y04.z3 ^ 2 := by
  linear_combination (5695 / 32) * eqP1 + (-42) * eqP2 + (34263 / 32) * eqP3


end V02

namespace V02

lemma G_d4' (k : ℕ) : iteratedDeriv 4 (Y04.G k) 0 =
    24 * ((-1 : ℝ) ^ k * Real.sqrt 2 / 4 * (121 / 384 + 23 / 12 * (k : ℝ) ^ 2 + 2 / 3 * (k : ℝ) ^ 4)) := by
  have := G_d4 k; linarith

lemma theta6 (k : ℕ) :
    (-1 : ℝ) ^ k * (π / Real.sqrt 2) *
        (8 / π * ∑ j ∈ Finset.Icc 1 5,
          iteratedDeriv (5 - j) (Y04.G k) 0 / ((5 - j).factorial : ℝ) / (j.factorial : ℝ) *
            Y04.I (5 - j) j) =
      (4039 / 256 + 733 / 8 * (k : ℝ) ^ 2 + 31 * (k : ℝ) ^ 4) * (π ^ 6 / 945)
        - (7 / 8 + 2 * (k : ℝ) ^ 2) * Y04.z3 ^ 2 := by
  set E := (-1 : ℝ) ^ k with hE
  set Y := (121 / 384 + 23 / 12 * (k : ℝ) ^ 2 + 2 / 3 * (k : ℝ) ^ 4) * Y04.I 4 1
    + (-(17 / 48 + (k : ℝ) ^ 2)) / 2 * Y04.I 3 2
    + (-(2 * (k : ℝ) ^ 2 + 5 / 8)) / 6 * Y04.I 2 3
    + (1 / 2) / 24 * Y04.I 1 4 + 1 / 120 * Y04.I 0 5 with hY
  have hsum : ∑ j ∈ Finset.Icc 1 5,
          iteratedDeriv (5 - j) (Y04.G k) 0 / ((5 - j).factorial : ℝ) / (j.factorial : ℝ) *
            Y04.I (5 - j) j = E * Real.sqrt 2 / 4 * Y := by
    rw [show Finset.Icc 1 5 = {1, 2, 3, 4, 5} by decide]
    rw [Finset.sum_insert (by decide), Finset.sum_insert (by decide), Finset.sum_insert (by decide),
      Finset.sum_insert (by decide), Finset.sum_singleton]
    simp only [show 5 - 1 = 4 from rfl, show 5 - 2 = 3 from rfl, show 5 - 3 = 2 from rfl,
      show 5 - 4 = 1 from rfl, show 5 - 5 = 0 from rfl, iteratedDeriv_zero, iteratedDeriv_one]
    rw [G_d4', Z02.G_d3', Z02.G_d2', Y04.G_d1, Y04.G_zero, hY]
    simp [Nat.factorial]
    ring
  rw [hsum]
  have hE2 : E ^ 2 = 1 := by rw [hE, ← pow_mul, mul_comm, pow_mul, neg_one_sq, one_pow]
  have hs : Real.sqrt 2 ≠ 0 := by positivity
  have hpi : π ≠ 0 := Real.pi_pos.ne'
  have e : E * (π / Real.sqrt 2) * (8 / π * (E * Real.sqrt 2 / 4 * Y)) = 2 * E ^ 2 * Y := by
    field_simp; ring
  rw [e, hE2, hY]
  linear_combination D60 - 2 * (k : ℝ) ^ 2 * D61 + 4 / 3 * (k : ℝ) ^ 4 * I41

lemma sum_eval6 (Gf : ℝ → ℝ) (τ : ℝ) :
    ∑ n ∈ Finset.Icc 2 6, τ ^ n * ∑ j ∈ Finset.Icc 1 (n - 1),
        iteratedDeriv (n - 1 - j) Gf 0 / ((n - 1 - j).factorial : ℝ) / (j.factorial : ℝ) * Y04.I (n - 1 - j) j
      = τ ^ 2 * Gf 0 * Y04.I 0 1 + τ ^ 3 * deriv Gf 0 * Y04.I 1 1 + τ ^ 3 * (Gf 0 / 2) * Y04.I 0 2
        + τ ^ 4 * (iteratedDeriv 2 Gf 0 / 2) * Y04.I 2 1
        + τ ^ 4 * (deriv Gf 0 / 2) * Y04.I 1 2 + τ ^ 4 * (Gf 0 / 6) * Y04.I 0 3
        + τ ^ 5 * (iteratedDeriv 3 Gf 0 / 6) * Y04.I 3 1 + τ ^ 5 * (iteratedDeriv 2 Gf 0 / 4) * Y04.I 2 2
        + τ ^ 5 * (deriv Gf 0 / 6) * Y04.I 1 3 + τ ^ 5 * (Gf 0 / 24) * Y04.I 0 4
        + τ ^ 6 * ∑ j ∈ Finset.Icc 1 5,
          iteratedDeriv (5 - j) Gf 0 / ((5 - j).factorial : ℝ) / (j.factorial : ℝ) * Y04.I (5 - j) j := by
  rw [Finset.sum_Icc_succ_top (by norm_num : 2 ≤ 5 + 1), Z02.sum_eval5]

lemma fourier6 (k : ℕ) : ∃ C τ₀ : ℝ, 0 < τ₀ ∧ ∀ p : ℝ, 1 < p → Y04.tau p ≤ τ₀ →
    |Y04.a k (Y04.h p) - (-1 : ℝ) ^ k * (Real.sqrt 2 / π) *
        (-(4 / (16 * (k : ℝ) ^ 2 - 1)) + (π ^ 2 / 6) * Y04.tau p ^ 2 + Y04.z3 * Y04.tau p ^ 3
          - (33 / 16 + 7 * (k : ℝ) ^ 2) * (π ^ 4 / 90) * Y04.tau p ^ 4
          + (-(3 / 4 + 4 * (k : ℝ) ^ 2) * Z02.z5
              - (7 / 4 + 4 * (k : ℝ) ^ 2) * ((π ^ 2 / 6) * Y04.z3)) * Y04.tau p ^ 5
          + ((4039 / 256 + 733 / 8 * (k : ℝ) ^ 2 + 31 * (k : ℝ) ^ 4) * (π ^ 6 / 945)
              - (7 / 8 + 2 * (k : ℝ) ^ 2) * Y04.z3 ^ 2) * Y04.tau p ^ 6)|
      ≤ C * Y04.tau p ^ 7 := by
  obtain ⟨CW, hW⟩ := watson6 (Y04.G k) (Y04.G_contDiff k) ⟨1, fun x _ => Y04.G_bdd k x⟩
  refine ⟨8 / π * CW, 1 / 2, by norm_num, fun p hp _ => ?_⟩
  have hpi := Real.pi_pos
  have hτ0 := Y04.tau_pos hp
  have hτ1 := Y04.tau_lt_one hp
  obtain ⟨_, hb⟩ := hW (Y04.tau p) hτ0 hτ1
  rw [sum_eval6] at hb
  have h2 := Y04.layer_formula k hp
  have h1 := Y04.A_formula k
  set X := Y04.tau p * ∫ s in Set.Ioi (0 : ℝ), ((1 + Real.exp (-s)) ^ Y04.tau p - 1) * Y04.G k (Y04.tau p * s)
    with hX
  set τ := Y04.tau p with hτ
  set C6 := ∑ j ∈ Finset.Icc 1 5,
    iteratedDeriv (5 - j) (Y04.G k) 0 / ((5 - j).factorial : ℝ) / (j.factorial : ℝ) * Y04.I (5 - j) j
    with hC6
  set Θ := (4039 / 256 + 733 / 8 * (k : ℝ) ^ 2 + 31 * (k : ℝ) ^ 4) * (π ^ 6 / 945)
        - (7 / 8 + 2 * (k : ℝ) ^ 2) * Y04.z3 ^ 2 with hΘ
  have hth := theta6 k
  rw [← hC6, ← hΘ] at hth
  have hE2 : ((-1 : ℝ) ^ k) ^ 2 = 1 := by rw [← pow_mul, mul_comm, pow_mul, neg_one_sq, one_pow]
  have hsq : Real.sqrt 2 ^ 2 = 2 := Real.sq_sqrt (by norm_num)
  have hs0 : Real.sqrt 2 ≠ 0 := by positivity
  have hC6v : C6 = (-1 : ℝ) ^ k * Real.sqrt 2 / 8 * Θ := by
    rw [← hth]; field_simp
    linear_combination (-C6) * hE2
  set S := τ ^ 2 * Y04.G k 0 * Y04.I 0 1 + τ ^ 3 * deriv (Y04.G k) 0 * Y04.I 1 1
        + τ ^ 3 * (Y04.G k 0 / 2) * Y04.I 0 2
        + τ ^ 4 * (iteratedDeriv 2 (Y04.G k) 0 / 2) * Y04.I 2 1
        + τ ^ 4 * (deriv (Y04.G k) 0 / 2) * Y04.I 1 2 + τ ^ 4 * (Y04.G k 0 / 6) * Y04.I 0 3
        + τ ^ 5 * (iteratedDeriv 3 (Y04.G k) 0 / 6) * Y04.I 3 1
        + τ ^ 5 * (iteratedDeriv 2 (Y04.G k) 0 / 4) * Y04.I 2 2
        + τ ^ 5 * (deriv (Y04.G k) 0 / 6) * Y04.I 1 3 + τ ^ 5 * (Y04.G k 0 / 24) * Y04.I 0 4
        + τ ^ 6 * C6 with hS
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
              - (7 / 4 + 4 * (k : ℝ) ^ 2) * ((π ^ 2 / 6) * Y04.z3)) * τ ^ 5 + Θ * τ ^ 6)
        = 8 / π * (X - S) := by
    have hd2 : iteratedDeriv 2 (Y04.G k) 0 = 2 * -((-1 : ℝ) ^ k * Real.sqrt 2 * ((k : ℝ) ^ 2 / 2 + 5 / 32)) :=
      Z02.G_d2' k
    rw [ha, h1, hS, hC6v, hd2, Z02.G_d3', Y04.G_zero k, Y04.G_d1 k, Y04.I01, Y04.I11, Y04.I02, Y04.I21,
      hI03, Z02.I31, hI22, hI13]
    field_simp
    ring
  rw [key, abs_mul, abs_of_pos (by positivity)]
  calc 8 / π * |X - S| ≤ 8 / π * (CW * τ ^ 7) := mul_le_mul_of_nonneg_left hb (by positivity)
    _ = _ := by ring


end V02

theorem V02.proof_1 : V02.stmt_1 := V02.G_d4

theorem V02.proof_2 : V02.stmt_2 := fun Gf hG hb => V02.watson6 Gf hG hb

theorem V02.proof_3 : V02.stmt_3 := by unfold V02.stmt_3; rw [V02.I41]; ring

theorem V02.proof_4 : V02.stmt_4 := V02.D61

theorem V02.proof_5 : V02.stmt_5 := V02.D60

theorem V02.proof_6 : V02.stmt_6 := V02.theta6

theorem V02.proof_7 : V02.stmt_7 := by
  have hpi := Real.pi_pos
  have n6 : (31 * (π ^ 6 / 945) : ℝ) ≠ 0 := by positivity
  refine ⟨Polynomial.C (31 * (π ^ 6 / 945)) * Polynomial.X ^ 2
      + Polynomial.C (733 / 8 * (π ^ 6 / 945) - 2 * Y04.z3 ^ 2) * Polynomial.X
      + Polynomial.C (4039 / 256 * (π ^ 6 / 945) - 7 / 8 * Y04.z3 ^ 2), ?_, ?_, ?_⟩
  · rw [Polynomial.degree_quadratic n6]; rfl
  · rw [Polynomial.leadingCoeff_quadratic n6]; norm_num
  · intro k
    obtain ⟨C, τ₀, hτ₀, hC⟩ := V02.fourier6 k
    refine ⟨C, τ₀, hτ₀, fun p hp hτ => ?_⟩
    have := hC p hp hτ
    simp only [Polynomial.eval_add, Polynomial.eval_mul, Polynomial.eval_C, Polynomial.eval_X,
      Polynomial.eval_pow]
    convert this using 4
    ring

theorem V02.proof_8 : V02.stmt_8 := V02.fourier6
