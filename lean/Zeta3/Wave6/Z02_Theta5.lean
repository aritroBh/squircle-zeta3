import Mathlib
import Zeta3.Wave5.Y04_Fourier
/-!
# Z02 — the `τ⁵` coefficient of the Fourier coefficients `a_k(h_p)`: `Θ₅`, and `prop:degree` for `n ≤ 5`

Target id: Z02.
Paper labels: `prop:degree` and `prop:theta` (main.tex lines 1214–1225), extending
`thm:fourier` (main.tex 1204–1211, formalised in `Zeta3.Wave5.Y04_Fourier`) by one order.
Detailed proofs: paper/ieee/proofs/sec_g2_degree.tex (Lemmas `g2:lem:d:coef`, `g2:lem:d:e`,
`g2:lem:d:theta`, `g2:lem:d:atoms`, proof of Proposition 2) and
paper/ieee/proofs/sec_g2_theta.tex (Lemma `g2:lem:t:red`, Table `g2:tab:D`, the remark on
the case `n = 5`); the Watson expansion is Lemma `g2:lem:f:watson` of sec_g2_fourier.tex.

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
  $1,\frac12,-\frac{23}{24},-\frac{25}{48},\frac{731}{1152}$"
  (so `ρ₃(u) = -17/48 - u`, i.e. `g_{k,3} = -(-1)^k (√2/4)(17/48 + k²)`).

Supplement, Lemma `g2:lem:t:red`, Table `g2:tab:D` and the remark on `n = 5`:

  $\wD_{n,m}=2\sum_{j=1}^{n-1-2m}\dfrac{e_{m,\,n-1-j}}{j!}\,I_{n-1-j,\,j}$
  $\wD_{5,0}$ & $-\tfrac34,\ -\tfrac74$ \\   $\wD_{5,1}$ & $2,\ 2$ \\
  (on the products $\zeta(5),\ \zeta(2)\zeta(3)$)
  $\wD_{5,1}=I_{3,1}+I_{2,2}$ ... and $I_{3,1}=3!\,\eta(5)=\frac{45}8\zeta(5)$

Supplement, proof of Proposition 2, case `n = 2M`:

  the coefficient of $u^{M-1}$ is ... $(-1)^{M-1}(2^{2M-1}-1)\zeta(2M)$, which is not $0$.

Reused from `Zeta3.Wave5.Y04_Fourier` (not redefined): `Y04.tau` (`τ = 1 - 1/p`), `Y04.h`
(the support function `h_p`), `Y04.a` (the Fourier coefficient `a_k(f) = (4/π)∫_0^{π/2} f cos 4kθ`),
`Y04.G` (the layer function `G_k(x) = cos(4kθ_x) cos²θ_x sin θ_x`, `θ_x = arctan e^{-x}`),
`Y04.I` (`I_{m,j} = ∫_0^∞ s^m log^j(1+e^{-s}) ds`), `Y04.z3` (`ζ(3) = (riemannZeta 3).re`).
Y04 proves the reduction `a_k(h_p) - A_k = (8/π) τ ∫_0^∞ [(1+e^{-s})^τ - 1] G_k(τ s) ds`, the
values `g_{k,0}, g_{k,1}, g_{k,2}`, the atoms of weight `≤ 4`, Watson at `N = 4` and `thm:fourier`.

Conventions.
* `ζ(3) = Y04.z3 = (riemannZeta 3).re` and `ζ(5) = Z02.z5 = (riemannZeta 5).re`, the real parts
  of Mathlib's complex zeta values; by `zeta_nat_eq_tsum_of_gt_one` these are the real numbers
  `∑ 1/n^3`, `∑ 1/n^5`. `ζ(2) = π²/6` and `ζ(4) = π⁴/90` are written out.
* `g_{k,m} = G_k^{(m)}(0)/m!` is written `iteratedDeriv m (Y04.G k) 0 / m!`; `G_k` is `C^∞`
  on `ℝ` (Y04 `stmt_3`), so these are genuine derivatives.
* All half-line integrals have integrable integrands (`Y04.stmt_4` proves integrability of every
  `s^m ℓ(s)^j`, `j ≥ 1`; `stmt_2` asserts the Watson integrand integrable), so no statement holds
  by the "non-integrable ↦ 0" convention.
* As in Y04, the Watson lemma (`stmt_2`) asks `G` to be `C^∞` on all of `ℝ` and bounded on
  `[0,∞)` (the supplement: `C^∞` on `[0,∞)`); this covers `G = G_k`.
* "The coefficient of `τ^n`" in `prop:degree` is the coefficient of an asymptotic expansion
  (the supplement's "Setting"). `stmt_7` states it that way: polynomials `Θ₂, …, Θ₅ ∈ ℝ[X]`
  such that `(-1)^k (√2/π)[-4/(16k²-1) + Σ_{n=2}^5 Θ_n(k²) τ^n]` approximates `a_k(h_p)` to
  `O_k(τ⁶)` for every `k`. Since such coefficients are unique and a polynomial is determined by
  its values at `k²`, `k ∈ ℕ`, this pins down the paper's `Θ_n` for `n ≤ 5`. "Degree exactly
  `d`" is `Polynomial.degree Θ_n = d` (so `Θ_n ≠ 0`).
* The `O(τ⁶)` claims are explicit: `∃ C τ₀, 0 < τ₀ ∧ ∀ p > 1, τ ≤ τ₀ → |error| ≤ C τ⁶`, the
  same form as `Y04.stmt_8`.

NOT COVERED: `prop:degree` for `n ≥ 6` and in general, and its clause "coefficients are
alternating multiple zeta values of weight `n`" (no AMZV/iterated-integral library is set up;
for `n ≤ 5` the explicit closed forms are stated instead, which are products of single zeta
values); `Θ₆, …, Θ₉` and the first sentence of `prop:theta` for `6 ≤ n ≤ 9`; the iterated-
integral form of `I_{m,j}` (Lemma `g2:lem:d:atoms`), shuffle/stuffle/duality and the exact
certificates (we state the two `n = 5` identities `D̃_{5,0}`, `D̃_{5,1}` directly);
the general formula for `e_{r,m}`, `ρ_m` (only `g_{k,3}` is stated); Cauchy's formula
`P = ∫ h_p` and the `k = 0` comparison with `thm:pw`.
-/

-- BEGIN STATEMENTS
namespace Z02

open Real MeasureTheory Set Polynomial

/-- `ζ(5)` as a real number: the real part of Mathlib's `riemannZeta 5`. -/
noncomputable def z5 : ℝ := (riemannZeta 5).re

/-- Lemma `g2:lem:d:e`(b) at `m = 3` (with `ρ₃(u) = e_{0,3} - 2u e_{1,3} = -17/48 - u`):
the next Taylor coefficient of `G_k`, `g_{k,3} = G_k'''(0)/3! = -(-1)^k (√2/4)(17/48 + k²)`. -/
def stmt_1 : Prop :=
  ∀ k : ℕ, iteratedDeriv 3 (Y04.G k) 0 / 6 =
    -((-1 : ℝ) ^ k * Real.sqrt 2 / 4 * (17 / 48 + (k : ℝ) ^ 2))

/-- Lemma `g2:lem:f:watson` with `N = 5`: for `G` smooth and bounded on `[0,∞)`, with
`g_m = G^{(m)}(0)/m!`, there is `C` with, for `0 < τ < 1`, the integrand integrable and
`|τ ∫_0^∞ [(1+e^{-s})^τ - 1] G(τ s) ds - Σ_{n=2}^{5} τ^n Σ_{j=1}^{n-1} g_{n-1-j}/j! I_{n-1-j,j}|
  ≤ C τ⁶`. -/
def stmt_2 : Prop :=
  ∀ Gf : ℝ → ℝ, ContDiff ℝ (⊤ : ℕ∞) Gf → (∃ M : ℝ, ∀ x : ℝ, 0 ≤ x → |Gf x| ≤ M) →
    ∃ C : ℝ, ∀ τ : ℝ, 0 < τ → τ < 1 →
      IntegrableOn (fun s : ℝ => ((1 + Real.exp (-s)) ^ τ - 1) * Gf (τ * s)) (Ioi 0) ∧
      |τ * (∫ s in Ioi (0 : ℝ), ((1 + Real.exp (-s)) ^ τ - 1) * Gf (τ * s)) -
          ∑ n ∈ Finset.Icc 2 5, τ ^ n * ∑ j ∈ Finset.Icc 1 (n - 1),
            iteratedDeriv (n - 1 - j) Gf 0 / ((n - 1 - j).factorial : ℝ) /
              (j.factorial : ℝ) * Y04.I (n - 1 - j) j| ≤ C * τ ^ 6

/-- Remark on the case `n = 5` (sec_g2_theta): `I_{3,1} = 3! η(5) = 45 ζ(5)/8`. -/
def stmt_3 : Prop := Y04.I 3 1 = 45 / 8 * z5

/-- Table `g2:tab:D`, row `D̃_{5,1}`: `D̃_{5,1} = I_{3,1} + I_{2,2} = 2 ζ(5) + 2 ζ(2) ζ(3)`. -/
def stmt_4 : Prop := Y04.I 3 1 + Y04.I 2 2 = 2 * z5 + 2 * (π ^ 2 / 6) * Y04.z3

/-- Table `g2:tab:D`, row `D̃_{5,0}`: with `e_{0,0..3} = 1, 1/2, -5/8, -17/48`,
`D̃_{5,0} = 2(e_{0,3} I_{3,1} + e_{0,2} I_{2,2}/2! + e_{0,1} I_{1,3}/3! + e_{0,0} I_{0,4}/4!)
  = -(3/4) ζ(5) - (7/4) ζ(2) ζ(3)`. -/
def stmt_5 : Prop :=
  2 * (-(17 / 48) * Y04.I 3 1 + (-(5 / 8)) * Y04.I 2 2 / 2 + (1 / 2) * Y04.I 1 3 / 6
      + Y04.I 0 4 / 24) =
    -(3 / 4) * z5 - 7 / 4 * ((π ^ 2 / 6) * Y04.z3)

/-- Lemmas `g2:lem:d:coef`, `g2:lem:d:theta` and `g2:lem:t:red` at `n = 5`, together with
`prop:theta`'s closed form: the `τ⁵` coefficient `c_5(k) = (8/π) Σ_{j=1}^{4} g_{k,4-j}/j! I_{4-j,j}`
satisfies `(-1)^k (π/√2) c_5(k) = Θ₅(k²) = -(3/4 + 4k²) ζ(5) - (7/4 + 4k²) ζ(2) ζ(3)`. -/
def stmt_6 : Prop :=
  ∀ k : ℕ,
    (-1 : ℝ) ^ k * (π / Real.sqrt 2) *
        (8 / π * ∑ j ∈ Finset.Icc 1 4,
          iteratedDeriv (4 - j) (Y04.G k) 0 / ((4 - j).factorial : ℝ) / (j.factorial : ℝ) *
            Y04.I (4 - j) j) =
      -(3 / 4 + 4 * (k : ℝ) ^ 2) * z5 - (7 / 4 + 4 * (k : ℝ) ^ 2) * ((π ^ 2 / 6) * Y04.z3)

/-- `prop:degree` for `n = 2, 3, 4, 5`: there are polynomials `Θ₂, Θ₃, Θ₄, Θ₅` with
`deg Θ_n = ⌊n/2⌋ - 1` exactly (`0, 0, 1, 1`), whose values at `u = k²` are the `τ^n`
coefficients of `(-1)^k (π/√2) a_k(h_p)` in the asymptotic sense: for every `k` there are
`C_k`, `τ₀ > 0` with
`|a_k(h_p) - (-1)^k (√2/π)[-4/(16k²-1) + Σ_{n=2}^{5} Θ_n(k²) τ^n]| ≤ C_k τ⁶` for `τ ≤ τ₀`;
and for `n = 2m` (`m = 1, 2`) the leading coefficient is `(-1)^{m-1}(2^{2m-1}-1) ζ(2m)`. -/
def stmt_7 : Prop :=
  ∃ Θ₂ Θ₃ Θ₄ Θ₅ : ℝ[X],
    Θ₂.degree = ((2 / 2 - 1 : ℕ) : WithBot ℕ) ∧ Θ₃.degree = ((3 / 2 - 1 : ℕ) : WithBot ℕ) ∧
    Θ₄.degree = ((4 / 2 - 1 : ℕ) : WithBot ℕ) ∧ Θ₅.degree = ((5 / 2 - 1 : ℕ) : WithBot ℕ) ∧
    Θ₂.leadingCoeff = (-1 : ℝ) ^ (1 - 1) * (2 ^ (2 * 1 - 1) - 1) * (π ^ 2 / 6) ∧
    Θ₄.leadingCoeff = (-1 : ℝ) ^ (2 - 1) * (2 ^ (2 * 2 - 1) - 1) * (π ^ 4 / 90) ∧
    ∀ k : ℕ, ∃ C τ₀ : ℝ, 0 < τ₀ ∧ ∀ p : ℝ, 1 < p → Y04.tau p ≤ τ₀ →
      |Y04.a k (Y04.h p) - (-1 : ℝ) ^ k * (Real.sqrt 2 / π) *
          (-(4 / (16 * (k : ℝ) ^ 2 - 1))
            + Θ₂.eval ((k : ℝ) ^ 2) * Y04.tau p ^ 2 + Θ₃.eval ((k : ℝ) ^ 2) * Y04.tau p ^ 3
            + Θ₄.eval ((k : ℝ) ^ 2) * Y04.tau p ^ 4 + Θ₅.eval ((k : ℝ) ^ 2) * Y04.tau p ^ 5)|
        ≤ C * Y04.tau p ^ 6

/-- `thm:fourier` extended by one order with `prop:theta`'s `Θ₅`: for every `k ≥ 0` there are
`C_k` and `τ₀ > 0` such that for every `p > 1` with `τ = 1 - 1/p ≤ τ₀`, with `u = k²`,
`|a_k(h_p) - (-1)^k (√2/π)[-4/(16k²-1) + ζ(2)τ² + ζ(3)τ³ - (33/16 + 7k²)ζ(4)τ⁴
   + (-(3/4 + 4u) ζ(5) - (7/4 + 4u) ζ(2) ζ(3)) τ⁵]| ≤ C_k τ⁶`. -/
def stmt_8 : Prop :=
  ∀ k : ℕ, ∃ C τ₀ : ℝ, 0 < τ₀ ∧ ∀ p : ℝ, 1 < p → Y04.tau p ≤ τ₀ →
    |Y04.a k (Y04.h p) - (-1 : ℝ) ^ k * (Real.sqrt 2 / π) *
        (-(4 / (16 * (k : ℝ) ^ 2 - 1)) + (π ^ 2 / 6) * Y04.tau p ^ 2 + Y04.z3 * Y04.tau p ^ 3
          - (33 / 16 + 7 * (k : ℝ) ^ 2) * (π ^ 4 / 90) * Y04.tau p ^ 4
          + (-(3 / 4 + 4 * (k : ℝ) ^ 2) * z5
              - (7 / 4 + 4 * (k : ℝ) ^ 2) * ((π ^ 2 / 6) * Y04.z3)) * Y04.tau p ^ 5)|
      ≤ C * Y04.tau p ^ 6

end Z02
-- END STATEMENTS

open Real MeasureTheory Set

namespace Z02


lemma zn_five : Y04.zn 5 = z5 := by simp [Y04.zn, z5]

lemma I31 : Y04.I 3 1 = 45 / 8 * z5 := by
  rw [Y04.I_m1, Y04.eta_val 5 (by norm_num), zn_five]; norm_num [Nat.factorial]; ring

noncomputable def Fddd (k : ℕ) (t : ℝ) : ℝ :=
  (4 * k) ^ 3 * Real.sin (4 * k * t) * Real.cos t ^ 2 * Real.sin t
    - (4 * k) ^ 2 * Real.cos (4 * k * t) * (Real.cos t ^ 3 - 2 * Real.cos t * Real.sin t ^ 2)
    + 4 * (4 * k) ^ 2 * Real.cos (4 * k * t) * Real.cos t * Real.sin t ^ 2
    + 4 * (4 * k) * Real.sin (4 * k * t) * (2 * Real.cos t ^ 2 * Real.sin t - Real.sin t ^ 3)
    - 2 * (4 * k) ^ 2 * Real.cos (4 * k * t) * Real.cos t ^ 3
    + 6 * (4 * k) * Real.sin (4 * k * t) * Real.cos t ^ 2 * Real.sin t
    - (4 * k) * Real.sin (4 * k * t) * (2 * Real.sin t ^ 3 - 7 * Real.cos t ^ 2 * Real.sin t)
    + Real.cos (4 * k * t) * (20 * Real.cos t * Real.sin t ^ 2 - 7 * Real.cos t ^ 3)

lemma Fdd_hasDeriv (k : ℕ) (t : ℝ) : HasDerivAt (Y04.Fdd k) (Fddd k t) t := by
  have hS := ((hasDerivAt_id t).const_mul (4 * (k : ℝ))).sin
  have hC := ((hasDerivAt_id t).const_mul (4 * (k : ℝ))).cos
  have hc := Real.hasDerivAt_cos t
  have hs := Real.hasDerivAt_sin t
  have h := ((((hC.const_mul (-(4 * (k : ℝ)) ^ 2)).mul (hc.pow 2)).mul hs).add
      ((((hS.const_mul (4 * (4 * (k : ℝ)))).mul hc).mul (hs.pow 2)))).sub
      ((hS.const_mul (2 * (4 * (k : ℝ)))).mul (hc.pow 3)) |>.add
      (hC.mul (((hs.pow 3).const_mul 2).sub (((hc.pow 2).const_mul 7).mul hs)))
  convert h using 1
  · funext s; simp [Y04.Fdd, mul_assoc]
  · simp only [Fddd, id, Pi.mul_apply, Pi.pow_apply, Pi.sub_apply]; push_cast; ring

noncomputable def Hdd (k : ℕ) (t : ℝ) : ℝ :=
  -(Fddd k t * (Real.sin t * Real.cos t) + 2 * Y04.Fdd k t * (Real.cos t * Real.cos t - Real.sin t * Real.sin t)
    + Y04.Fd k t * (-(4 * Real.sin t * Real.cos t)))

lemma Hd_hasDeriv (k : ℕ) (t : ℝ) : HasDerivAt (Y04.Hd k) (Hdd k t) t := by
  have hc := Real.hasDerivAt_cos t
  have hs := Real.hasDerivAt_sin t
  have h : HasDerivAt (Y04.Hd k) _ t :=
    (((Fdd_hasDeriv k t).mul (hs.mul hc)).add ((Y04.Fd_hasDeriv k t).mul ((hc.mul hc).sub (hs.mul hs)))).neg
  refine h.congr_deriv ?_
  simp only [Hdd, Pi.mul_apply, Pi.sub_apply]; ring

noncomputable def Kf (k : ℕ) (t : ℝ) : ℝ := Y04.Hd k t * (-(Real.sin t * Real.cos t))

noncomputable def Kfd (k : ℕ) (t : ℝ) : ℝ :=
  Hdd k t * (-(Real.sin t * Real.cos t)) + Y04.Hd k t * (-(Real.cos t * Real.cos t - Real.sin t * Real.sin t))

lemma Kf_hasDeriv (k : ℕ) (t : ℝ) : HasDerivAt (Kf k) (Kfd k t) t := by
  have hc := Real.hasDerivAt_cos t
  have hs := Real.hasDerivAt_sin t
  have h : HasDerivAt (Kf k) _ t := (Hd_hasDeriv k t).mul ((hs.mul hc).neg)
  refine h.congr_deriv ?_
  simp only [Kfd]; ring

lemma G_dd_fun (k : ℕ) : iteratedDeriv 2 (Y04.G k) = fun x => Kf k (Y04.thetaX x) := by
  rw [iteratedDeriv_succ, iteratedDeriv_one, Y04.deriv_G]
  funext x
  exact ((Y04.H_hasDeriv k _).comp x (Y04.theta_hasDeriv x)).deriv

lemma G_d3 (k : ℕ) : iteratedDeriv 3 (Y04.G k) 0 / 6 =
    -((-1 : ℝ) ^ k * Real.sqrt 2 / 4 * (17 / 48 + (k : ℝ) ^ 2)) := by
  rw [iteratedDeriv_succ, G_dd_fun]
  have hc := (Kf_hasDeriv k (Y04.thetaX 0)).comp 0 (Y04.theta_hasDeriv 0)
  rw [show (fun x => Kf k (Y04.thetaX x)) = Kf k ∘ Y04.thetaX from rfl, hc.deriv]
  simp only [Y04.thetaX_zero, Kfd, Hdd, Fddd, Y04.Hd, Y04.Fdd, Y04.Fd, Y04.cos_k, Y04.sin_k,
    Real.cos_pi_div_four, Real.sin_pi_div_four]
  have : Real.sqrt 2 ^ 2 = 2 := Real.sq_sqrt (by norm_num)
  have h3 : Real.sqrt 2 ^ 3 = 2 * Real.sqrt 2 := by
    rw [show Real.sqrt 2 ^ 3 = Real.sqrt 2 ^ 2 * Real.sqrt 2 by ring, this]
  have h5 : Real.sqrt 2 ^ 5 = 4 * Real.sqrt 2 := by
    rw [show Real.sqrt 2 ^ 5 = (Real.sqrt 2 ^ 2) ^ 2 * Real.sqrt 2 by ring, this]; ring
  have h7 : Real.sqrt 2 ^ 7 = 8 * Real.sqrt 2 := by
    rw [show Real.sqrt 2 ^ 7 = (Real.sqrt 2 ^ 2) ^ 3 * Real.sqrt 2 by ring, this]; ring
  have h9 : Real.sqrt 2 ^ 9 = 16 * Real.sqrt 2 := by
    rw [show Real.sqrt 2 ^ 9 = (Real.sqrt 2 ^ 2) ^ 4 * Real.sqrt 2 by ring, this]; ring
  ring_nf
  rw [h9]; ring

end Z02

namespace Z02

/-- Taylor's formula at `0` to order 3 with a quartic remainder on `[0,∞)`. -/
lemma taylor3 {Gf : ℝ → ℝ} (hG : ContDiff ℝ (⊤ : ℕ∞) Gf) {M : ℝ}
    (hM : ∀ x : ℝ, 0 ≤ x → |Gf x| ≤ M) :
    ∃ C : ℝ, 0 ≤ C ∧ ∀ x : ℝ, 0 ≤ x →
      |Gf x - (iteratedDeriv 0 Gf 0 / (Nat.factorial 0 : ℝ) + iteratedDeriv 1 Gf 0 / (Nat.factorial 1 : ℝ) * x
        + iteratedDeriv 2 Gf 0 / (Nat.factorial 2 : ℝ) * x ^ 2
        + iteratedDeriv 3 Gf 0 / (Nat.factorial 3 : ℝ) * x ^ 3)| ≤ C * x ^ 4 := by
  set g0 := iteratedDeriv 0 Gf 0 / (Nat.factorial 0 : ℝ)
  set g1 := iteratedDeriv 1 Gf 0 / (Nat.factorial 1 : ℝ)
  set g2 := iteratedDeriv 2 Gf 0 / (Nat.factorial 2 : ℝ)
  set g3 := iteratedDeriv 3 Gf 0 / (Nat.factorial 3 : ℝ)
  have hc4 : Continuous (iteratedDeriv 4 Gf) :=
    hG.continuous_iteratedDeriv 4 (WithTop.coe_le_coe.2 le_top)
  obtain ⟨K4, hK4⟩ := (isCompact_Icc (a := (0 : ℝ)) (b := 1)).exists_bound_of_continuousOn
    hc4.continuousOn
  have hM0 : 0 ≤ M := (abs_nonneg _).trans (hM 0 le_rfl)
  refine ⟨max (|K4| / 24) (M + |g0| + |g1| + |g2| + |g3|), le_max_of_le_left (by positivity), ?_⟩
  intro x hx
  rcases eq_or_lt_of_le hx with h0 | hpos
  · subst h0; simp [g0]
  rcases le_or_gt x 1 with hx1 | hx1
  · have hcd : ContDiffOn ℝ (3 + 1) Gf (uIcc 0 x) := (hG.of_le (WithTop.coe_le_coe.2 le_top)).contDiffOn
    obtain ⟨x', hx', heq⟩ := taylor_mean_remainder_lagrange_iteratedDeriv (n := 3) hpos.ne hcd
    have hu : UniqueDiffOn ℝ (uIcc 0 x) := uniqueDiffOn_uIcc hpos.ne
    have hmem : (0 : ℝ) ∈ uIcc 0 x := left_mem_uIcc
    have hT : taylorWithinEval Gf 3 (uIcc 0 x) 0 x = g0 + g1 * x + g2 * x ^ 2 + g3 * x ^ 3 := by
      rw [taylor_within_apply]
      simp only [Finset.sum_range_succ, Finset.sum_range_zero, zero_add, sub_zero, smul_eq_mul]
      rw [iteratedDerivWithin_eq_iteratedDeriv hu (hG.contDiffAt.of_le (WithTop.coe_le_coe.2 le_top)) hmem,
        iteratedDerivWithin_eq_iteratedDeriv hu (hG.contDiffAt.of_le (WithTop.coe_le_coe.2 le_top)) hmem,
        iteratedDerivWithin_eq_iteratedDeriv hu (hG.contDiffAt.of_le (WithTop.coe_le_coe.2 le_top)) hmem,
        iteratedDerivWithin_eq_iteratedDeriv hu (hG.contDiffAt.of_le (WithTop.coe_le_coe.2 le_top)) hmem]
      simp only [g0, g1, g2, g3]; norm_num; ring
    rw [hT] at heq
    rw [heq]
    rw [uIoo_of_le hx] at hx'
    have hb := hK4 x' ⟨hx'.1.le, by linarith [hx'.2]⟩
    rw [Real.norm_eq_abs] at hb
    rw [sub_zero, abs_div, abs_mul, abs_pow, abs_of_pos hpos]
    norm_num [Nat.factorial]
    have : |iteratedDeriv 4 Gf x'| * x ^ 4 / 24 ≤ |K4| / 24 * x ^ 4 := by
      have := hb.trans (le_abs_self K4)
      have h3 : 0 ≤ x ^ 4 := by positivity
      nlinarith
    calc |iteratedDeriv 4 Gf x'| * x ^ 4 / 24 ≤ |K4| / 24 * x ^ 4 := this
      _ ≤ _ := by gcongr; exact le_max_left _ _
  · have h1 : 1 ≤ x := hx1.le
    have hx2' : 1 ≤ x ^ 2 := by nlinarith
    have hx3' : 1 ≤ x ^ 3 := by nlinarith
    have hx14 : x ≤ x ^ 4 := by
      rw [show x ^ 4 = x * x ^ 3 by ring]; nlinarith
    have hx24 : x ^ 2 ≤ x ^ 4 := by
      rw [show x ^ 4 = x ^ 2 * x ^ 2 by ring]; nlinarith
    have hx34 : x ^ 3 ≤ x ^ 4 := by
      rw [show x ^ 4 = x * x ^ 3 by ring]; nlinarith
    have h14 : 1 ≤ x ^ 4 := by linarith
    calc _ ≤ |Gf x| + (|g0| + |g1| * x + |g2| * x ^ 2 + |g3| * x ^ 3) := by
          have := abs_sub (Gf x) (g0 + g1 * x + g2 * x ^ 2 + g3 * x ^ 3)
          have a0 := abs_add_le (g0 + g1 * x + g2 * x ^ 2) (g3 * x ^ 3)
          have a1 := abs_add_le (g0 + g1 * x) (g2 * x ^ 2)
          have a2 := abs_add_le g0 (g1 * x)
          rw [abs_mul, abs_of_pos hpos] at a2
          rw [abs_mul, abs_of_nonneg (by positivity : (0:ℝ) ≤ x ^ 2)] at a1
          rw [abs_mul, abs_of_nonneg (by positivity : (0:ℝ) ≤ x ^ 3)] at a0
          linarith
      _ ≤ (M + |g0| + |g1| + |g2| + |g3|) * x ^ 4 := by
          have := hM x hx
          have e1 : |g1| * x ≤ |g1| * x ^ 4 := mul_le_mul_of_nonneg_left hx14 (abs_nonneg _)
          have e2 : |g2| * x ^ 2 ≤ |g2| * x ^ 4 := mul_le_mul_of_nonneg_left hx24 (abs_nonneg _)
          have e5 : |g3| * x ^ 3 ≤ |g3| * x ^ 4 := mul_le_mul_of_nonneg_left hx34 (abs_nonneg _)
          have e3 : M ≤ M * x ^ 4 := le_mul_of_one_le_right hM0 h14
          have e4 : |g0| ≤ |g0| * x ^ 4 := le_mul_of_one_le_right (abs_nonneg _) h14
          nlinarith
      _ ≤ _ := by gcongr; exact le_max_right _ _

/-- the dominating function of the order-5 Watson remainder -/
noncomputable def B6 (s : ℝ) : ℝ := (1 + s + s ^ 2 + s ^ 3 + s ^ 4) * Real.exp (-s)

lemma B6_int : IntegrableOn B6 (Set.Ioi 0) := by
  have h := ((((Y04.intOn_pow_exp1 0).add (Y04.intOn_pow_exp1 1)).add (Y04.intOn_pow_exp1 2)).add
    (Y04.intOn_pow_exp1 3)).add (Y04.intOn_pow_exp1 4)
  refine IntegrableOn.congr_fun h (fun s _ => ?_) measurableSet_Ioi
  simp only [B6, Pi.add_apply]; ring

lemma B6_nonneg {s : ℝ} (hs : 0 ≤ s) : 0 ≤ B6 s := by unfold B6; positivity

/-- the pointwise order-5 Watson remainder -/
lemma watson_pt5 {g0 g1 g2 g3 CR M τ s : ℝ} (hCR : 0 ≤ CR) (hM : 0 ≤ M) (hτ0 : 0 < τ) (hτ1 : τ < 1)
    (hs : 0 ≤ s) {Gv : ℝ} (hGv : |Gv| ≤ M)
    (hR : |Gv - (g0 + g1 * (τ * s) + g2 * (τ * s) ^ 2 + g3 * (τ * s) ^ 3)| ≤ CR * (τ * s) ^ 4) :
    |τ * ((Real.exp (τ * Y04.ell s) - 1) * Gv)
      - (τ ^ 2 * g0 * (s ^ 0 * Y04.ell s ^ 1) + τ ^ 3 * g1 * (s ^ 1 * Y04.ell s ^ 1)
        + τ ^ 3 * (g0 / 2) * (s ^ 0 * Y04.ell s ^ 2) + τ ^ 4 * g2 * (s ^ 2 * Y04.ell s ^ 1)
        + τ ^ 4 * (g1 / 2) * (s ^ 1 * Y04.ell s ^ 2) + τ ^ 4 * (g0 / 6) * (s ^ 0 * Y04.ell s ^ 3)
        + τ ^ 5 * g3 * (s ^ 3 * Y04.ell s ^ 1) + τ ^ 5 * (g2 / 2) * (s ^ 2 * Y04.ell s ^ 2)
        + τ ^ 5 * (g1 / 6) * (s ^ 1 * Y04.ell s ^ 3) + τ ^ 5 * (g0 / 24) * (s ^ 0 * Y04.ell s ^ 4))|
      ≤ τ ^ 6 * ((|g1| + |g2| + |g3| + 2 * CR + M) * B6 s) := by
  have hl0 := Y04.ell_nonneg s
  have hl1 := Y04.ell_le_one s hs
  have hle := Y04.ell_le s
  set L := Y04.ell s
  have hτL : |τ * L| ≤ 1 := by rw [abs_of_nonneg (by positivity)]; nlinarith
  have he := Real.exp_bound hτL (n := 5) (by norm_num)
  simp only [Finset.sum_range_succ, Finset.sum_range_zero, Nat.factorial] at he
  norm_num at he
  set a1 := τ * L + (τ * L) ^ 2 / 2 + (τ * L) ^ 3 / 6 + (τ * L) ^ 4 / 24 with ha1
  set r1 := Real.exp (τ * L) - 1 - a1 with hr1
  have hr1b : |r1| ≤ τ ^ 5 * L ^ 5 := by
    have : r1 = Real.exp (τ * L) - (1 + τ * L + (τ * L) ^ 2 / 2 + (τ * L) ^ 3 / 6 + (τ * L) ^ 4 / 24) := by
      rw [hr1, ha1]; ring
    rw [this]
    have h4 : (|τ| * |L|) ^ 5 = τ ^ 5 * L ^ 5 := by
      rw [abs_of_pos hτ0, abs_of_nonneg hl0]; ring
    rw [h4] at he
    have : 0 ≤ τ ^ 5 * L ^ 5 := by positivity
    linarith
  set q := g0 + g1 * (τ * s) + g2 * (τ * s) ^ 2 + g3 * (τ * s) ^ 3
  set R := Gv - q with hRdef
  set Q := g1 * s * (L ^ 4 / 24) + g2 * s ^ 2 * (L ^ 3 / 6 + τ * L ^ 4 / 24)
    + g3 * s ^ 3 * (L ^ 2 / 2 + τ * L ^ 3 / 6 + τ ^ 2 * L ^ 4 / 24) with hQ
  have key : τ * ((Real.exp (τ * L) - 1) * Gv)
      - (τ ^ 2 * g0 * (s ^ 0 * L ^ 1) + τ ^ 3 * g1 * (s ^ 1 * L ^ 1)
        + τ ^ 3 * (g0 / 2) * (s ^ 0 * L ^ 2) + τ ^ 4 * g2 * (s ^ 2 * L ^ 1)
        + τ ^ 4 * (g1 / 2) * (s ^ 1 * L ^ 2) + τ ^ 4 * (g0 / 6) * (s ^ 0 * L ^ 3)
        + τ ^ 5 * g3 * (s ^ 3 * L ^ 1) + τ ^ 5 * (g2 / 2) * (s ^ 2 * L ^ 2)
        + τ ^ 5 * (g1 / 6) * (s ^ 1 * L ^ 3) + τ ^ 5 * (g0 / 24) * (s ^ 0 * L ^ 4))
      = τ ^ 6 * Q + τ * a1 * R + τ * r1 * Gv := by
    rw [hRdef, hr1, ha1, hQ]; ring
  rw [key]
  have hes : 0 ≤ Real.exp (-s) := (Real.exp_pos _).le
  have hL2 : L ^ 2 ≤ Real.exp (-s) := (pow_le_of_le_one hl0 hl1 (by norm_num)).trans hle
  have hL3 : L ^ 3 ≤ Real.exp (-s) := (pow_le_of_le_one hl0 hl1 (by norm_num)).trans hle
  have hL4 : L ^ 4 ≤ Real.exp (-s) := (pow_le_of_le_one hl0 hl1 (by norm_num)).trans hle
  have hL5 : L ^ 5 ≤ Real.exp (-s) := (pow_le_of_le_one hl0 hl1 (by norm_num)).trans hle
  have hτ' : 0 ≤ τ := hτ0.le
  have hτL3 : τ * L ^ 3 ≤ L ^ 3 := mul_le_of_le_one_left (by positivity) hτ1.le
  have hτL4 : τ * L ^ 4 ≤ L ^ 4 := mul_le_of_le_one_left (by positivity) hτ1.le
  have hτ2 : τ ^ 2 ≤ 1 := pow_le_one₀ hτ' hτ1.le
  have hτ2L4 : τ ^ 2 * L ^ 4 ≤ L ^ 4 := mul_le_of_le_one_left (by positivity) hτ2
  -- term 1
  have b1 : |τ ^ 6 * Q| ≤ τ ^ 6 * ((|g1| + |g2| + |g3|) * (s + s ^ 2 + s ^ 3) * Real.exp (-s)) := by
    rw [abs_mul, abs_of_pos (by positivity)]
    apply mul_le_mul_of_nonneg_left _ (by positivity)
    have u1 : 0 ≤ L ^ 4 / 24 := by positivity
    have u1' : L ^ 4 / 24 ≤ Real.exp (-s) := by linarith
    have u2 : 0 ≤ L ^ 3 / 6 + τ * L ^ 4 / 24 := by positivity
    have u2' : L ^ 3 / 6 + τ * L ^ 4 / 24 ≤ Real.exp (-s) := by linarith
    have u3 : 0 ≤ L ^ 2 / 2 + τ * L ^ 3 / 6 + τ ^ 2 * L ^ 4 / 24 := by positivity
    have u3' : L ^ 2 / 2 + τ * L ^ 3 / 6 + τ ^ 2 * L ^ 4 / 24 ≤ Real.exp (-s) := by linarith
    have v1 : |g1 * s * (L ^ 4 / 24)| ≤ |g1| * (s * Real.exp (-s)) := by
      rw [abs_mul, abs_mul, abs_of_nonneg hs, abs_of_nonneg u1, mul_assoc]
      exact mul_le_mul_of_nonneg_left (mul_le_mul_of_nonneg_left u1' hs) (abs_nonneg _)
    have v2 : |g2 * s ^ 2 * (L ^ 3 / 6 + τ * L ^ 4 / 24)| ≤ |g2| * (s ^ 2 * Real.exp (-s)) := by
      rw [abs_mul, abs_mul, abs_of_nonneg (by positivity : (0:ℝ) ≤ s ^ 2), abs_of_nonneg u2, mul_assoc]
      exact mul_le_mul_of_nonneg_left (mul_le_mul_of_nonneg_left u2' (by positivity)) (abs_nonneg _)
    have v3 : |g3 * s ^ 3 * (L ^ 2 / 2 + τ * L ^ 3 / 6 + τ ^ 2 * L ^ 4 / 24)|
        ≤ |g3| * (s ^ 3 * Real.exp (-s)) := by
      rw [abs_mul, abs_mul, abs_of_nonneg (by positivity : (0:ℝ) ≤ s ^ 3), abs_of_nonneg u3, mul_assoc]
      exact mul_le_mul_of_nonneg_left (mul_le_mul_of_nonneg_left u3' (by positivity)) (abs_nonneg _)
    have t1 := abs_add_le (g1 * s * (L ^ 4 / 24) + g2 * s ^ 2 * (L ^ 3 / 6 + τ * L ^ 4 / 24))
      (g3 * s ^ 3 * (L ^ 2 / 2 + τ * L ^ 3 / 6 + τ ^ 2 * L ^ 4 / 24))
    have t2 := abs_add_le (g1 * s * (L ^ 4 / 24)) (g2 * s ^ 2 * (L ^ 3 / 6 + τ * L ^ 4 / 24))
    have ag1 := abs_nonneg g1; have ag2 := abs_nonneg g2; have ag3 := abs_nonneg g3
    have w : (|g1| + |g2| + |g3|) * (s + s ^ 2 + s ^ 3) * Real.exp (-s)
        - (|g1| * (s * Real.exp (-s)) + |g2| * (s ^ 2 * Real.exp (-s)) + |g3| * (s ^ 3 * Real.exp (-s)))
        = Real.exp (-s) * (|g1| * (s ^ 2 + s ^ 3) + |g2| * (s + s ^ 3) + |g3| * (s + s ^ 2)) := by ring
    have : 0 ≤ Real.exp (-s) * (|g1| * (s ^ 2 + s ^ 3) + |g2| * (s + s ^ 3) + |g3| * (s + s ^ 2)) := by
      positivity
    linarith
  -- term 2
  have b2 : |τ * a1 * R| ≤ τ ^ 6 * (2 * CR * s ^ 4 * Real.exp (-s)) := by
    have ha1b : |a1| ≤ 2 * τ * L := by
      have hu0 : 0 ≤ τ * L := mul_nonneg hτ0.le hl0
      have hu1 : τ * L ≤ 1 := by
        have := mul_le_mul hτ1.le hl1 hl0 zero_le_one; linarith
      have p2 : (τ * L) ^ 2 ≤ τ * L := pow_le_of_le_one hu0 hu1 (by norm_num)
      have p3 : (τ * L) ^ 3 ≤ τ * L := pow_le_of_le_one hu0 hu1 (by norm_num)
      have p4 : (τ * L) ^ 4 ≤ τ * L := pow_le_of_le_one hu0 hu1 (by norm_num)
      have q2 : 0 ≤ (τ * L) ^ 2 := sq_nonneg _
      have q3 : 0 ≤ (τ * L) ^ 3 := pow_nonneg hu0 3
      have q4 : 0 ≤ (τ * L) ^ 4 := pow_nonneg hu0 4
      have : 0 ≤ a1 := by rw [ha1]; linarith
      rw [abs_of_nonneg this, ha1]
      have : 2 * τ * L = 2 * (τ * L) := by ring
      linarith
    have hRb : |R| ≤ CR * (τ * s) ^ 4 := hR
    rw [abs_mul, abs_mul, abs_of_pos hτ0]
    calc τ * |a1| * |R| ≤ τ * (2 * τ * L) * (CR * (τ * s) ^ 4) :=
          mul_le_mul (mul_le_mul_of_nonneg_left ha1b hτ') hRb (abs_nonneg _) (by positivity)
      _ = τ ^ 6 * (2 * CR * s ^ 4 * L) := by ring
      _ ≤ τ ^ 6 * (2 * CR * s ^ 4 * Real.exp (-s)) :=
          mul_le_mul_of_nonneg_left (mul_le_mul_of_nonneg_left hle (by positivity)) (by positivity)
  -- term 3
  have b3 : |τ * r1 * Gv| ≤ τ ^ 6 * (M * Real.exp (-s)) := by
    rw [abs_mul, abs_mul, abs_of_pos hτ0]
    calc τ * |r1| * |Gv| ≤ τ * (τ ^ 5 * L ^ 5) * M :=
          mul_le_mul (mul_le_mul_of_nonneg_left hr1b hτ0.le) hGv (abs_nonneg _) (by positivity)
      _ = τ ^ 6 * (L ^ 5 * M) := by ring
      _ ≤ τ ^ 6 * (Real.exp (-s) * M) :=
          mul_le_mul_of_nonneg_left (mul_le_mul_of_nonneg_right hL5 hM) (by positivity)
      _ = _ := by ring
  have t1 := abs_add_le (τ ^ 6 * Q + τ * a1 * R) (τ * r1 * Gv)
  have t2 := abs_add_le (τ ^ 6 * Q) (τ * a1 * R)
  have hfin : τ ^ 6 * ((|g1| + |g2| + |g3|) * (s + s ^ 2 + s ^ 3) * Real.exp (-s))
      + τ ^ 6 * (2 * CR * s ^ 4 * Real.exp (-s))
      + τ ^ 6 * (M * Real.exp (-s)) ≤ τ ^ 6 * ((|g1| + |g2| + |g3| + 2 * CR + M) * B6 s) := by
    unfold B6
    apply le_of_sub_nonneg
    have : τ ^ 6 * ((|g1| + |g2| + |g3| + 2 * CR + M) * ((1 + s + s ^ 2 + s ^ 3 + s ^ 4) * Real.exp (-s)))
        - (τ ^ 6 * ((|g1| + |g2| + |g3|) * (s + s ^ 2 + s ^ 3) * Real.exp (-s))
        + τ ^ 6 * (2 * CR * s ^ 4 * Real.exp (-s)) + τ ^ 6 * (M * Real.exp (-s)))
        = τ ^ 6 * Real.exp (-s) * ((|g1| + |g2| + |g3|) * (1 + s ^ 4) + 2 * CR * (1 + s + s ^ 2 + s ^ 3)
          + M * (s + s ^ 2 + s ^ 3 + s ^ 4)) := by ring
    rw [this]
    have := abs_nonneg g1; have := abs_nonneg g2; have := abs_nonneg g3
    positivity
  linarith

lemma int10 {f1 f2 f3 f4 f5 f6 f7 f8 f9 f10 : ℝ → ℝ} (c1 c2 c3 c4 c5 c6 c7 c8 c9 c10 : ℝ)
    (h1 : IntegrableOn f1 (Set.Ioi 0)) (h2 : IntegrableOn f2 (Set.Ioi 0)) (h3 : IntegrableOn f3 (Set.Ioi 0))
    (h4 : IntegrableOn f4 (Set.Ioi 0)) (h5 : IntegrableOn f5 (Set.Ioi 0)) (h6 : IntegrableOn f6 (Set.Ioi 0))
    (h7 : IntegrableOn f7 (Set.Ioi 0)) (h8 : IntegrableOn f8 (Set.Ioi 0)) (h9 : IntegrableOn f9 (Set.Ioi 0))
    (h10 : IntegrableOn f10 (Set.Ioi 0)) :
    IntegrableOn (fun s => c1 * f1 s + c2 * f2 s + c3 * f3 s + c4 * f4 s + c5 * f5 s + c6 * f6 s
      + c7 * f7 s + c8 * f8 s + c9 * f9 s + c10 * f10 s) (Set.Ioi 0) ∧
    ∫ s in Set.Ioi (0 : ℝ), (c1 * f1 s + c2 * f2 s + c3 * f3 s + c4 * f4 s + c5 * f5 s + c6 * f6 s
      + c7 * f7 s + c8 * f8 s + c9 * f9 s + c10 * f10 s)
      = c1 * (∫ s in Set.Ioi (0 : ℝ), f1 s) + c2 * (∫ s in Set.Ioi (0 : ℝ), f2 s)
        + c3 * (∫ s in Set.Ioi (0 : ℝ), f3 s) + c4 * (∫ s in Set.Ioi (0 : ℝ), f4 s)
        + c5 * (∫ s in Set.Ioi (0 : ℝ), f5 s) + c6 * (∫ s in Set.Ioi (0 : ℝ), f6 s)
        + c7 * (∫ s in Set.Ioi (0 : ℝ), f7 s) + c8 * (∫ s in Set.Ioi (0 : ℝ), f8 s)
        + c9 * (∫ s in Set.Ioi (0 : ℝ), f9 s) + c10 * (∫ s in Set.Ioi (0 : ℝ), f10 s) := by
  obtain ⟨b6, e6⟩ := Y04.int6 c1 c2 c3 c4 c5 c6 h1 h2 h3 h4 h5 h6
  have a7 : Integrable (fun s => c7 * f7 s) (volume.restrict (Set.Ioi 0)) := h7.const_mul c7
  have a8 : Integrable (fun s => c8 * f8 s) (volume.restrict (Set.Ioi 0)) := h8.const_mul c8
  have a9 : Integrable (fun s => c9 * f9 s) (volume.restrict (Set.Ioi 0)) := h9.const_mul c9
  have a10 : Integrable (fun s => c10 * f10 s) (volume.restrict (Set.Ioi 0)) := h10.const_mul c10
  have b7 : Integrable (fun s => c1 * f1 s + c2 * f2 s + c3 * f3 s + c4 * f4 s + c5 * f5 s + c6 * f6 s
      + c7 * f7 s) (volume.restrict (Set.Ioi 0)) := Integrable.add b6 a7
  have b8 : Integrable (fun s => c1 * f1 s + c2 * f2 s + c3 * f3 s + c4 * f4 s + c5 * f5 s + c6 * f6 s
      + c7 * f7 s + c8 * f8 s) (volume.restrict (Set.Ioi 0)) := b7.add a8
  have b9 : Integrable (fun s => c1 * f1 s + c2 * f2 s + c3 * f3 s + c4 * f4 s + c5 * f5 s + c6 * f6 s
      + c7 * f7 s + c8 * f8 s + c9 * f9 s) (volume.restrict (Set.Ioi 0)) := b8.add a9
  have b10 : Integrable (fun s => c1 * f1 s + c2 * f2 s + c3 * f3 s + c4 * f4 s + c5 * f5 s + c6 * f6 s
      + c7 * f7 s + c8 * f8 s + c9 * f9 s + c10 * f10 s) (volume.restrict (Set.Ioi 0)) := b9.add a10
  refine ⟨b10, ?_⟩
  rw [integral_add b9 a10, integral_add b8 a9, integral_add b7 a8, integral_add b6 a7, e6,
    integral_const_mul, integral_const_mul, integral_const_mul, integral_const_mul]

lemma watson5 (Gf : ℝ → ℝ) (hG : ContDiff ℝ (⊤ : ℕ∞) Gf) (hb : ∃ M : ℝ, ∀ x : ℝ, 0 ≤ x → |Gf x| ≤ M) :
    ∃ C : ℝ, ∀ τ : ℝ, 0 < τ → τ < 1 →
      IntegrableOn (fun s : ℝ => ((1 + Real.exp (-s)) ^ τ - 1) * Gf (τ * s)) (Ioi 0) ∧
      |τ * (∫ s in Ioi (0 : ℝ), ((1 + Real.exp (-s)) ^ τ - 1) * Gf (τ * s)) -
          ∑ n ∈ Finset.Icc 2 5, τ ^ n * ∑ j ∈ Finset.Icc 1 (n - 1),
            iteratedDeriv (n - 1 - j) Gf 0 / ((n - 1 - j).factorial : ℝ) /
              (j.factorial : ℝ) * Y04.I (n - 1 - j) j| ≤ C * τ ^ 6 := by
  obtain ⟨M, hM⟩ := hb
  have hM0 : 0 ≤ M := (abs_nonneg _).trans (hM 0 le_rfl)
  obtain ⟨CR, hCR, hT⟩ := taylor3 hG hM
  set g0 := iteratedDeriv 0 Gf 0 / (Nat.factorial 0 : ℝ) with hg0
  set g1 := iteratedDeriv 1 Gf 0 / (Nat.factorial 1 : ℝ) with hg1
  set g2 := iteratedDeriv 2 Gf 0 / (Nat.factorial 2 : ℝ) with hg2
  set g3 := iteratedDeriv 3 Gf 0 / (Nat.factorial 3 : ℝ) with hg3
  refine ⟨(|g1| + |g2| + |g3| + 2 * CR + M) * ∫ s in Set.Ioi (0 : ℝ), B6 s, fun τ hτ0 hτ1 => ?_⟩
  have hGc : Continuous Gf := hG.continuous
  have e : (fun s : ℝ => ((1 + Real.exp (-s)) ^ τ - 1) * Gf (τ * s))
      = fun s => (Real.exp (τ * Y04.ell s) - 1) * Gf (τ * s) := by
    funext s; rw [Y04.rpow_eq_exp]
  rw [e]
  have hint : IntegrableOn (fun s => (Real.exp (τ * Y04.ell s) - 1) * Gf (τ * s)) (Set.Ioi 0) := by
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
  obtain ⟨hmi, hmv⟩ := int10 (τ ^ 2 * g0) (τ ^ 3 * g1) (τ ^ 3 * (g0 / 2)) (τ ^ 4 * g2) (τ ^ 4 * (g1 / 2))
    (τ ^ 4 * (g0 / 6)) (τ ^ 5 * g3) (τ ^ 5 * (g2 / 2)) (τ ^ 5 * (g1 / 6)) (τ ^ 5 * (g0 / 24))
    (Y04.I_int 0 1 le_rfl) (Y04.I_int 1 1 le_rfl) (Y04.I_int 0 2 (by norm_num))
    (Y04.I_int 2 1 le_rfl) (Y04.I_int 1 2 (by norm_num)) (Y04.I_int 0 3 (by norm_num))
    (Y04.I_int 3 1 le_rfl) (Y04.I_int 2 2 (by norm_num)) (Y04.I_int 1 3 (by norm_num))
    (Y04.I_int 0 4 (by norm_num))
  have hsum : ∑ n ∈ Finset.Icc 2 5, τ ^ n * ∑ j ∈ Finset.Icc 1 (n - 1),
            iteratedDeriv (n - 1 - j) Gf 0 / ((n - 1 - j).factorial : ℝ) /
              (j.factorial : ℝ) * Y04.I (n - 1 - j) j
      = τ ^ 2 * g0 * Y04.I 0 1 + τ ^ 3 * g1 * Y04.I 1 1 + τ ^ 3 * (g0 / 2) * Y04.I 0 2
        + τ ^ 4 * g2 * Y04.I 2 1 + τ ^ 4 * (g1 / 2) * Y04.I 1 2 + τ ^ 4 * (g0 / 6) * Y04.I 0 3
        + τ ^ 5 * g3 * Y04.I 3 1 + τ ^ 5 * (g2 / 2) * Y04.I 2 2 + τ ^ 5 * (g1 / 6) * Y04.I 1 3
        + τ ^ 5 * (g0 / 24) * Y04.I 0 4 := by
    rw [show Finset.Icc 2 5 = {2, 3, 4, 5} by decide]
    rw [Finset.sum_insert (by decide), Finset.sum_insert (by decide), Finset.sum_insert (by decide),
      Finset.sum_singleton]
    rw [show Finset.Icc 1 (2 - 1) = {1} by decide, show Finset.Icc 1 (3 - 1) = {1, 2} by decide,
      show Finset.Icc 1 (4 - 1) = {1, 2, 3} by decide, show Finset.Icc 1 (5 - 1) = {1, 2, 3, 4} by decide]
    simp [Finset.sum_insert, hg0, hg1, hg2, hg3, Nat.factorial]
    ring
  rw [hsum]
  unfold Y04.I
  rw [← hmv, ← integral_const_mul, ← integral_sub (hint.const_mul τ) hmi]
  have hB := (B6_int.const_mul ((|g1| + |g2| + |g3| + 2 * CR + M))).const_mul (τ ^ 6)
  have := norm_integral_le_of_norm_le hB
    ((ae_restrict_iff' measurableSet_Ioi).2 (Filter.Eventually.of_forall
      (fun s (hs : 0 < s) => by
        rw [Real.norm_eq_abs]
        have hx : 0 ≤ τ * s := by positivity
        have := watson_pt5 (g0 := g0) (g1 := g1) (g2 := g2) (g3 := g3) hCR hM0 hτ0 hτ1 hs.le
          (hM (τ * s) hx)
          (by have := hT (τ * s) hx; simpa [hg0, hg1, hg2, hg3] using this)
        convert this using 2)))
  rw [Real.norm_eq_abs, integral_const_mul, integral_const_mul] at this
  calc _ ≤ _ := this
    _ = _ := by ring

end Z02

namespace Z02

lemma logi_eq (s : ℝ) : Y04.logi s = Real.exp (-Y04.ell s) := by
  unfold Y04.logi Y04.ell
  rw [Real.exp_neg (Real.log (1 + Real.exp (-s))), Real.exp_log (by positivity), one_div]

lemma one_sub_logi (s : ℝ) : 1 - Y04.logi s = Real.exp (-(s + Y04.ell s)) := by
  rw [neg_add, Real.exp_add, ← logi_eq]
  unfold Y04.logi
  have : 0 < 1 + Real.exp (-s) := by positivity
  field_simp; ring

/-- one half of the Beta integral, on the half line -/
lemma beta_half {a b : ℝ} (ha : 0 < a) (hb : 0 < b) :
    IntegrableOn (fun s => Real.exp (-(a + b) * Y04.ell s - b * s)) (Set.Ioi 0) ∧
    IntegrableOn (fun t : ℝ => t ^ (a - 1) * (1 - t) ^ (b - 1)) (Set.Ioo (1 / 2) 1) ∧
    ∫ s in Set.Ioi (0 : ℝ), Real.exp (-(a + b) * Y04.ell s - b * s) =
      ∫ t in Set.Ioo (1 / 2 : ℝ) 1, t ^ (a - 1) * (1 - t) ^ (b - 1) := by
  set g : ℝ → ℝ := fun t => t ^ (a - 1) * (1 - t) ^ (b - 1) with hg
  have hder : ∀ s ∈ Set.Ioi (0 : ℝ), HasDerivWithinAt Y04.logi (Y04.logid s) (Set.Ioi 0) s :=
    fun s _ => (Y04.logi_hasDeriv s).hasDerivWithinAt
  have hinj : InjOn Y04.logi (Set.Ioi 0) := Y04.logi_strictMono.injective.injOn
  have hpt : ∀ s, |Y04.logid s| • g (Y04.logi s) = Real.exp (-(a + b) * Y04.ell s - b * s) := by
    intro s
    simp only [hg, smul_eq_mul]
    rw [Y04.logid_eq, one_sub_logi, logi_eq, ← Real.exp_mul, ← Real.exp_mul, ← Real.exp_add,
      abs_of_pos (Real.exp_pos _), ← Real.exp_add, ← Real.exp_add]
    congr 1; ring
  have hint : IntegrableOn (fun s => Real.exp (-(a + b) * Y04.ell s - b * s)) (Set.Ioi 0) := by
    have hb' : IntegrableOn (fun s : ℝ => Real.exp (-b * s)) (Set.Ioi 0) :=
      integrableOn_exp_mul_Ioi (by linarith) 0
    refine Integrable.mono' hb' ((Real.continuous_exp.comp ((continuous_const.mul
      Y04.continuous_ell).sub (continuous_const.mul continuous_id))).aestronglyMeasurable) ?_
    refine (ae_restrict_iff' measurableSet_Ioi).2 (Filter.Eventually.of_forall (fun s _ => ?_))
    rw [Real.norm_eq_abs, abs_of_pos (Real.exp_pos _)]
    apply Real.exp_le_exp.2
    have := Y04.ell_nonneg s
    nlinarith
  have himg := integral_image_eq_integral_abs_deriv_smul measurableSet_Ioi hder hinj g
  have himgI := (integrableOn_image_iff_integrableOn_abs_deriv_smul measurableSet_Ioi hder hinj g).2
    (by simp_rw [hpt]; exact hint)
  simp_rw [hpt] at himg
  rw [Y04.logi_image] at himg himgI
  exact ⟨hint, himgI, himg.symm⟩

/-- `∫_0^∞ e^{-(a+b)ℓ - bs} + ∫_0^∞ e^{-(a+b)ℓ - as} = Γ(a)Γ(b)/Γ(a+b)` -/
lemma beta_full {a b : ℝ} (ha : 0 < a) (hb : 0 < b) :
    (∫ s in Set.Ioi (0 : ℝ), Real.exp (-(a + b) * Y04.ell s - b * s)) +
      (∫ s in Set.Ioi (0 : ℝ), Real.exp (-(b + a) * Y04.ell s - a * s)) =
      Real.Gamma a * Real.Gamma b / Real.Gamma (a + b) := by
  obtain ⟨_, hI1, hv1⟩ := beta_half ha hb
  obtain ⟨_, hI2, hv2⟩ := beta_half hb ha
  set g : ℝ → ℝ := fun t => t ^ (a - 1) * (1 - t) ^ (b - 1) with hg
  set g' : ℝ → ℝ := fun t => t ^ (b - 1) * (1 - t) ^ (a - 1) with hg'
  have hsym : ∀ t, g (1 - t) = g' t := fun t => by simp only [hg, hg']; ring_nf
  have hII : IntervalIntegrable g volume (1 / 2) 1 := by
    rw [intervalIntegrable_iff_integrableOn_Ioc_of_le (by norm_num)]
    exact (integrableOn_Ioc_iff_integrableOn_Ioo (by simp)).2 hI1
  have hII' : IntervalIntegrable g' volume (1 / 2) 1 := by
    rw [intervalIntegrable_iff_integrableOn_Ioc_of_le (by norm_num)]
    exact (integrableOn_Ioc_iff_integrableOn_Ioo (by simp)).2 hI2
  have hII2 : IntervalIntegrable g volume 0 (1 / 2) := by
    have := hII'.comp_sub_left 1
    have e : (fun x => g' (1 - x)) = g := by
      funext x; simp only [hg, hg']; ring_nf
    rw [e] at this
    norm_num at this
    exact this.symm
  have hsplit := intervalIntegral.integral_add_adjacent_intervals hII2 hII
  have hrefl : ∫ t in (0:ℝ)..(1 / 2), g t = ∫ t in (1 / 2 : ℝ)..1, g' t := by
    have := intervalIntegral.integral_comp_sub_left (a := (1 / 2 : ℝ)) (b := 1) g 1
    simp only [hsym] at this
    norm_num at this
    exact this.symm
  have hbeta := Y04.beta_real ha hb
  have hG : 0 < Real.Gamma (a + b) := Real.Gamma_pos_of_pos (by linarith)
  have hset : ∀ f : ℝ → ℝ, ∫ t in Set.Ioo (1 / 2 : ℝ) 1, f t = ∫ t in (1 / 2 : ℝ)..1, f t := fun f => by
    rw [intervalIntegral.integral_of_le (by norm_num), integral_Ioc_eq_integral_Ioo]
  rw [hv1, hv2, hset, hset]
  rw [← hsplit, hrefl] at hbeta
  rw [eq_div_iff hG.ne', ← hbeta]
  simp only [hg, hg']
  ring

end Z02

namespace Z02


lemma em_bound {q : ℝ} (hq : 0 ≤ q) :
    |Real.exp (-q) - (1 - q + q ^ 2 / 2 - q ^ 3 / 6)| ≤ 4 * q ^ 4 := by
  rcases le_or_gt q 1 with h | h
  · have h1 : |-q| ≤ 1 := by rw [abs_neg, abs_of_nonneg hq]; exact h
    have he := Real.exp_bound h1 (n := 4) (by norm_num)
    simp only [Finset.sum_range_succ, Finset.sum_range_zero, Nat.factorial] at he
    norm_num at he
    rw [abs_of_nonneg hq] at he
    have e : Real.exp (-q) - (1 - q + q ^ 2 / 2 - q ^ 3 / 6)
        = Real.exp (-q) - (1 + -q + q ^ 2 / 2 + (-q) ^ 3 / 6) := by ring
    rw [e]
    have : 0 ≤ q ^ 4 := by positivity
    nlinarith
  · have he : Real.exp (-q) ≤ 1 := Real.exp_le_one_iff.2 (by linarith)
    have he0 := Real.exp_pos (-q)
    have h2 : q ≤ q ^ 2 := by nlinarith
    have h3 : q ^ 2 ≤ q ^ 3 := by nlinarith
    have h4 : q ^ 3 ≤ q ^ 4 := by nlinarith
    rw [abs_le]
    constructor <;> nlinarith

lemma ep_bound {p : ℝ} (hp0 : 0 ≤ p) (hp1 : p ≤ 1) :
    |Real.exp (-p) - 1 - (-p + p ^ 2 / 2 - p ^ 3 / 6 + p ^ 4 / 24)| ≤ p ^ 5 := by
  have h1 : |-p| ≤ 1 := by rw [abs_neg, abs_of_nonneg hp0]; exact hp1
  have he := Real.exp_bound h1 (n := 5) (by norm_num)
  simp only [Finset.sum_range_succ, Finset.sum_range_zero, Nat.factorial] at he
  norm_num at he
  rw [abs_of_nonneg hp0] at he
  have e : Real.exp (-p) - 1 - (-p + p ^ 2 / 2 - p ^ 3 / 6 + p ^ 4 / 24)
      = Real.exp (-p) - (1 + -p + p ^ 2 / 2 + (-p) ^ 3 / 6 + p ^ 4 / 24) := by ring
  rw [e]
  have : 0 ≤ p ^ 5 := by positivity
  nlinarith

lemma mono_bd {u v t L s : ℝ} (hu0 : 0 ≤ u) (hut : u ≤ t) (hv0 : 0 ≤ v) (hvt : v ≤ t) (ht1 : t ≤ 1)
    (hL0 : 0 ≤ L) (hL1 : L ≤ 1) (hLe : L ≤ Real.exp (-s)) (hs : 0 ≤ s) {i j : ℕ} (hi : 1 ≤ i)
    (hij : 5 ≤ i + j) :
    (u * L) ^ i * (v * s) ^ j ≤ t ^ 5 * (s ^ j * Real.exp (-s)) := by
  have ht0 : 0 ≤ t := hu0.trans hut
  rw [mul_pow, mul_pow]
  have a1 : u ^ i ≤ t ^ i := pow_le_pow_left₀ hu0 hut i
  have a2 : v ^ j ≤ t ^ j := pow_le_pow_left₀ hv0 hvt j
  have a3 : t ^ i * t ^ j ≤ t ^ 5 := by
    rw [← pow_add]; exact pow_le_pow_of_le_one ht0 ht1 hij
  have a4 : L ^ i ≤ Real.exp (-s) := (pow_le_of_le_one hL0 hL1 (by omega)).trans hLe
  have b1 : u ^ i * v ^ j ≤ t ^ 5 :=
    (mul_le_mul a1 a2 (by positivity) (by positivity)).trans a3
  calc u ^ i * L ^ i * (v ^ j * s ^ j) = (u ^ i * v ^ j) * L ^ i * s ^ j := by ring
    _ ≤ t ^ 5 * Real.exp (-s) * s ^ j := by
        apply mul_le_mul_of_nonneg_right _ (by positivity)
        exact mul_le_mul b1 a4 (by positivity) (by positivity)
    _ = _ := by ring

/-- the degree-`≤ 4` part of `(e^{-uℓ} - 1) e^{-vs}` -/
noncomputable def P4 (u v s : ℝ) : ℝ :=
  (-u) * (s ^ 0 * Y04.ell s ^ 1) + (u * v) * (s ^ 1 * Y04.ell s ^ 1)
    + (-(u * v ^ 2 / 2)) * (s ^ 2 * Y04.ell s ^ 1) + (u * v ^ 3 / 6) * (s ^ 3 * Y04.ell s ^ 1)
    + (u ^ 2 / 2) * (s ^ 0 * Y04.ell s ^ 2) + (-(u ^ 2 * v / 2)) * (s ^ 1 * Y04.ell s ^ 2)
    + (u ^ 2 * v ^ 2 / 4) * (s ^ 2 * Y04.ell s ^ 2) + (-(u ^ 3 / 6)) * (s ^ 0 * Y04.ell s ^ 3)
    + (u ^ 3 * v / 6) * (s ^ 1 * Y04.ell s ^ 3) + (u ^ 4 / 24) * (s ^ 0 * Y04.ell s ^ 4)

lemma pt4 {u v t s : ℝ} (hu0 : 0 ≤ u) (hut : u ≤ t) (hv0 : 0 ≤ v) (hvt : v ≤ t) (ht1 : t ≤ 1)
    (hs : 0 ≤ s) :
    |(Real.exp (-u * Y04.ell s) - 1) * Real.exp (-v * s) - P4 u v s| ≤ 10 * t ^ 5 * B6 s := by
  have hl0 := Y04.ell_nonneg s
  have hl1 := Y04.ell_le_one s hs
  have hle := Y04.ell_le s
  have ht0 : 0 ≤ t := hu0.trans hut
  unfold P4
  set L := Y04.ell s
  set p := u * L with hp
  set q := v * s with hq
  have hp0 : 0 ≤ p := by positivity
  have hp1 : p ≤ 1 := by
    have := mul_le_mul (hut.trans ht1) hl1 hl0 zero_le_one; rw [hp]; linarith
  have hq0 : 0 ≤ q := by positivity
  set A := -p + p ^ 2 / 2 - p ^ 3 / 6 + p ^ 4 / 24 with hA
  set r1 := Real.exp (-p) - 1 - A with hr1
  set B := 1 - q + q ^ 2 / 2 - q ^ 3 / 6 with hB
  set r2 := Real.exp (-q) - B with hr2
  have hr1b : |r1| ≤ p ^ 5 := ep_bound hp0 hp1
  have hr2b : |r2| ≤ 4 * q ^ 4 := em_bound hq0
  set HT := -(p ^ 4 * q) / 24 - p ^ 3 * q ^ 2 / 12 + p ^ 4 * q ^ 2 / 48 - p ^ 2 * q ^ 3 / 12
    + p ^ 3 * q ^ 3 / 36 - p ^ 4 * q ^ 3 / 144 with hHT
  have e1 : -u * L = -p := by rw [hp]; ring
  have e2 : -v * s = -q := by rw [hq]; ring
  rw [e1, e2]
  have key : (Real.exp (-p) - 1) * Real.exp (-q)
      - ((-u) * (s ^ 0 * L ^ 1) + (u * v) * (s ^ 1 * L ^ 1)
    + (-(u * v ^ 2 / 2)) * (s ^ 2 * L ^ 1) + (u * v ^ 3 / 6) * (s ^ 3 * L ^ 1)
    + (u ^ 2 / 2) * (s ^ 0 * L ^ 2) + (-(u ^ 2 * v / 2)) * (s ^ 1 * L ^ 2)
    + (u ^ 2 * v ^ 2 / 4) * (s ^ 2 * L ^ 2) + (-(u ^ 3 / 6)) * (s ^ 0 * L ^ 3)
    + (u ^ 3 * v / 6) * (s ^ 1 * L ^ 3) + (u ^ 4 / 24) * (s ^ 0 * L ^ 4))
      = HT + A * r2 + r1 * Real.exp (-q) := by
    rw [hHT, hr1, hr2, hA, hB, hp, hq]; ring
  rw [key]
  have hes : 0 < Real.exp (-s) := Real.exp_pos _
  -- monomial bounds
  have m41 := mono_bd hu0 hut hv0 hvt ht1 hl0 hl1 hle hs (i := 4) (j := 1) (by norm_num) (by norm_num)
  have m32 := mono_bd hu0 hut hv0 hvt ht1 hl0 hl1 hle hs (i := 3) (j := 2) (by norm_num) (by norm_num)
  have m42 := mono_bd hu0 hut hv0 hvt ht1 hl0 hl1 hle hs (i := 4) (j := 2) (by norm_num) (by norm_num)
  have m23 := mono_bd hu0 hut hv0 hvt ht1 hl0 hl1 hle hs (i := 2) (j := 3) (by norm_num) (by norm_num)
  have m33 := mono_bd hu0 hut hv0 hvt ht1 hl0 hl1 hle hs (i := 3) (j := 3) (by norm_num) (by norm_num)
  have m43 := mono_bd hu0 hut hv0 hvt ht1 hl0 hl1 hle hs (i := 4) (j := 3) (by norm_num) (by norm_num)
  rw [← hp, ← hq] at m41 m32 m42 m23 m33 m43
  have n41 : 0 ≤ p ^ 4 * q := by positivity
  have n32 : 0 ≤ p ^ 3 * q ^ 2 := by positivity
  have n42 : 0 ≤ p ^ 4 * q ^ 2 := by positivity
  have n23 : 0 ≤ p ^ 2 * q ^ 3 := by positivity
  have n33 : 0 ≤ p ^ 3 * q ^ 3 := by positivity
  have n43 : 0 ≤ p ^ 4 * q ^ 3 := by positivity
  have t5 : 0 ≤ t ^ 5 := by positivity
  have bHT : |HT| ≤ t ^ 5 * ((s + s ^ 2 + s ^ 3) * Real.exp (-s)) := by
    have w : t ^ 5 * ((s + s ^ 2 + s ^ 3) * Real.exp (-s))
        = t ^ 5 * (s ^ 1 * Real.exp (-s)) + t ^ 5 * (s ^ 2 * Real.exp (-s))
          + t ^ 5 * (s ^ 3 * Real.exp (-s)) := by ring
    have z1 : 0 ≤ t ^ 5 * (s ^ 1 * Real.exp (-s)) := by positivity
    have z2 : 0 ≤ t ^ 5 * (s ^ 2 * Real.exp (-s)) := by positivity
    have z3 : 0 ≤ t ^ 5 * (s ^ 3 * Real.exp (-s)) := by positivity
    rw [abs_le, hHT]; constructor <;> linarith
  -- A r2
  have hAb : |A| ≤ 2 * p := by
    have p2 : p ^ 2 ≤ p := pow_le_of_le_one hp0 hp1 (by norm_num)
    have p3 : p ^ 3 ≤ p := pow_le_of_le_one hp0 hp1 (by norm_num)
    have p4 : p ^ 4 ≤ p := pow_le_of_le_one hp0 hp1 (by norm_num)
    have : 0 ≤ p ^ 2 := by positivity
    have : 0 ≤ p ^ 3 := by positivity
    have : 0 ≤ p ^ 4 := by positivity
    rw [abs_le, hA]; constructor <;> linarith
  have bAr : |A * r2| ≤ 8 * t ^ 5 * (s ^ 4 * Real.exp (-s)) := by
    rw [abs_mul]
    have h1 : |A| * |r2| ≤ 2 * p * (4 * q ^ 4) :=
      mul_le_mul hAb hr2b (abs_nonneg _) (by positivity)
    have h2 : p * q ^ 4 ≤ t ^ 5 * (s ^ 4 * Real.exp (-s)) := by
      have := mono_bd hu0 hut hv0 hvt ht1 hl0 hl1 hle hs (i := 1) (j := 4) (by norm_num) (by norm_num)
      rw [← hp, ← hq, pow_one] at this; exact this
    calc |A| * |r2| ≤ 2 * p * (4 * q ^ 4) := h1
      _ = 8 * (p * q ^ 4) := by ring
      _ ≤ 8 * (t ^ 5 * (s ^ 4 * Real.exp (-s))) := by linarith
      _ = _ := by ring
  -- r1 e^{-q}
  have br1 : |r1 * Real.exp (-q)| ≤ t ^ 5 * Real.exp (-s) := by
    rw [abs_mul, abs_of_pos (Real.exp_pos _)]
    have h1 : Real.exp (-q) ≤ 1 := Real.exp_le_one_iff.2 (by linarith)
    have h2 : p ^ 5 ≤ t ^ 5 * Real.exp (-s) := by
      have := mono_bd hu0 hut hv0 hvt ht1 hl0 hl1 hle hs (i := 5) (j := 0) (by norm_num) (by norm_num)
      rw [← hp, ← hq, pow_zero, mul_one, pow_zero, one_mul] at this; exact this
    calc |r1| * Real.exp (-q) ≤ p ^ 5 * 1 :=
          mul_le_mul hr1b h1 (Real.exp_pos _).le (by positivity)
      _ ≤ _ := by linarith
  have t1 := abs_add_le (HT + A * r2) (r1 * Real.exp (-q))
  have t2 := abs_add_le HT (A * r2)
  have fin : t ^ 5 * ((s + s ^ 2 + s ^ 3) * Real.exp (-s)) + 8 * t ^ 5 * (s ^ 4 * Real.exp (-s))
      + t ^ 5 * Real.exp (-s) ≤ 10 * t ^ 5 * B6 s := by
    unfold B6
    have : 10 * t ^ 5 * ((1 + s + s ^ 2 + s ^ 3 + s ^ 4) * Real.exp (-s))
        - (t ^ 5 * ((s + s ^ 2 + s ^ 3) * Real.exp (-s)) + 8 * t ^ 5 * (s ^ 4 * Real.exp (-s))
          + t ^ 5 * Real.exp (-s))
        = t ^ 5 * Real.exp (-s) * (9 + 9 * s + 9 * s ^ 2 + 9 * s ^ 3 + 2 * s ^ 4) := by ring
    have : 0 ≤ t ^ 5 * Real.exp (-s) * (9 + 9 * s + 9 * s ^ 2 + 9 * s ^ 3 + 2 * s ^ 4) := by positivity
    linarith
  linarith

end Z02


namespace Z02

/-- the integrated expansion of `(e^{-uℓ} - 1) e^{-vs}` to order `4` -/
lemma int_expand {u v t : ℝ} (hu0 : 0 ≤ u) (hut : u ≤ t) (hv0 : 0 ≤ v) (hvt : v ≤ t) (ht1 : t ≤ 1) :
    IntegrableOn (fun s => (Real.exp (-u * Y04.ell s) - 1) * Real.exp (-v * s)) (Set.Ioi 0) ∧
    |(∫ s in Set.Ioi (0 : ℝ), (Real.exp (-u * Y04.ell s) - 1) * Real.exp (-v * s))
      - ((-u) * Y04.I 0 1 + (u * v) * Y04.I 1 1 + (-(u * v ^ 2 / 2)) * Y04.I 2 1
        + (u * v ^ 3 / 6) * Y04.I 3 1 + (u ^ 2 / 2) * Y04.I 0 2 + (-(u ^ 2 * v / 2)) * Y04.I 1 2
        + (u ^ 2 * v ^ 2 / 4) * Y04.I 2 2 + (-(u ^ 3 / 6)) * Y04.I 0 3
        + (u ^ 3 * v / 6) * Y04.I 1 3 + (u ^ 4 / 24) * Y04.I 0 4)|
      ≤ 10 * t ^ 5 * ∫ s in Set.Ioi (0 : ℝ), B6 s := by
  have hint : IntegrableOn (fun s => (Real.exp (-u * Y04.ell s) - 1) * Real.exp (-v * s)) (Set.Ioi 0) := by
    refine Y04.intOn_of_bound (by have := Y04.continuous_ell; fun_prop)
      ((Y04.intOn_pow_exp1 0).const_mul 2) (fun s hs => ?_)
    have hl0 := Y04.ell_nonneg s
    have hl1 := Y04.ell_le_one s hs.le
    have hle := Y04.ell_le s
    have hu1 : u ≤ 1 := hut.trans ht1
    have hx : |-u * Y04.ell s| ≤ 1 := by
      rw [abs_le]; constructor <;> nlinarith
    have h1 := Real.abs_exp_sub_one_le hx
    have h2 : |-u * Y04.ell s| ≤ Real.exp (-s) := by
      rw [abs_mul, abs_neg, abs_of_nonneg hu0, abs_of_nonneg hl0]
      have := mul_le_mul hu1 hle hl0 zero_le_one; linarith
    have h3 : |Real.exp (-v * s)| ≤ 1 := by
      rw [abs_of_pos (Real.exp_pos _)]; exact Real.exp_le_one_iff.2 (by nlinarith)
    rw [abs_mul]
    simp only [pow_zero, one_mul]
    calc |Real.exp (-u * Y04.ell s) - 1| * |Real.exp (-v * s)| ≤ (2 * Real.exp (-s)) * 1 :=
          mul_le_mul (by linarith) h3 (abs_nonneg _) (by positivity)
      _ = 2 * Real.exp (-s) := by ring
  refine ⟨hint, ?_⟩
  obtain ⟨hmi, hmv⟩ := int10 (-u) (u * v) (-(u * v ^ 2 / 2)) (u * v ^ 3 / 6) (u ^ 2 / 2)
    (-(u ^ 2 * v / 2)) (u ^ 2 * v ^ 2 / 4) (-(u ^ 3 / 6)) (u ^ 3 * v / 6) (u ^ 4 / 24)
    (Y04.I_int 0 1 le_rfl) (Y04.I_int 1 1 le_rfl) (Y04.I_int 2 1 le_rfl) (Y04.I_int 3 1 le_rfl)
    (Y04.I_int 0 2 (by norm_num)) (Y04.I_int 1 2 (by norm_num)) (Y04.I_int 2 2 (by norm_num))
    (Y04.I_int 0 3 (by norm_num)) (Y04.I_int 1 3 (by norm_num)) (Y04.I_int 0 4 (by norm_num))
  unfold Y04.I
  rw [← hmv, ← integral_sub hint hmi]
  have hB := (B6_int.const_mul 10).const_mul (t ^ 5)
  have := norm_integral_le_of_norm_le hB
    ((ae_restrict_iff' measurableSet_Ioi).2 (Filter.Eventually.of_forall
      (fun s (hs : 0 < s) => by
        rw [Real.norm_eq_abs]
        have := pt4 hu0 hut hv0 hvt ht1 hs.le
        unfold P4 at this
        calc _ = _ := rfl
          _ ≤ _ := this
          _ = t ^ 5 * (10 * B6 s) := by ring)))
  rw [Real.norm_eq_abs, integral_const_mul, integral_const_mul] at this
  calc _ ≤ _ := this
    _ = _ := by ring

/-- the two-parameter Beta identity, regularised: for `a, b > 0`,
`∫_0^∞ (e^{-(a+b)ℓ}-1)(e^{-as}+e^{-bs}) ds = (a+b)/(ab) (Γ(1+a)Γ(1+b)/Γ(1+a+b) - 1)` -/
lemma lhs_closed {a b : ℝ} (ha : 0 < a) (hb : 0 < b) :
    (∫ s in Set.Ioi (0 : ℝ), (Real.exp (-(a + b) * Y04.ell s) - 1) * Real.exp (-a * s))
      + (∫ s in Set.Ioi (0 : ℝ), (Real.exp (-(a + b) * Y04.ell s) - 1) * Real.exp (-b * s))
      = (a + b) / (a * b) * (Real.Gamma (1 + a) * Real.Gamma (1 + b) / Real.Gamma (1 + (a + b)) - 1) := by
  obtain ⟨i1, _, _⟩ := beta_half ha hb
  obtain ⟨i2, _, _⟩ := beta_half hb ha
  have hfull := beta_full ha hb
  have e1 : ∀ s, (Real.exp (-(a + b) * Y04.ell s) - 1) * Real.exp (-a * s)
      = Real.exp (-(b + a) * Y04.ell s - a * s) - Real.exp (-a * s) := fun s => by
    rw [sub_mul, one_mul, ← Real.exp_add]; congr 2; ring
  have e2 : ∀ s, (Real.exp (-(a + b) * Y04.ell s) - 1) * Real.exp (-b * s)
      = Real.exp (-(a + b) * Y04.ell s - b * s) - Real.exp (-b * s) := fun s => by
    rw [sub_mul, one_mul, ← Real.exp_add]; congr 2; ring
  simp_rw [e1, e2]
  have ha' : IntegrableOn (fun s : ℝ => Real.exp (-a * s)) (Set.Ioi 0) :=
    integrableOn_exp_mul_Ioi (by linarith) 0
  have hb' : IntegrableOn (fun s : ℝ => Real.exp (-b * s)) (Set.Ioi 0) :=
    integrableOn_exp_mul_Ioi (by linarith) 0
  rw [integral_sub i2 ha', integral_sub i1 hb', integral_exp_mul_Ioi (by linarith) 0,
    integral_exp_mul_Ioi (by linarith) 0]
  have hsum : (∫ s in Set.Ioi (0 : ℝ), Real.exp (-(b + a) * Y04.ell s - a * s))
      + (∫ s in Set.Ioi (0 : ℝ), Real.exp (-(a + b) * Y04.ell s - b * s))
      = Real.Gamma a * Real.Gamma b / Real.Gamma (a + b) := by rw [← hfull]; ring
  have hG : 0 < Real.Gamma (a + b) := Real.Gamma_pos_of_pos (by linarith)
  rw [add_comm 1 a, add_comm 1 b, add_comm 1 (a + b), Real.Gamma_add_one ha.ne',
    Real.Gamma_add_one hb.ne', Real.Gamma_add_one (by linarith : a + b ≠ 0)]
  simp only [mul_zero, Real.exp_zero]
  have hsub : ∀ X Y Z W : ℝ, X - Z + (Y - W) = (X + Y) - Z - W := fun _ _ _ _ => by ring
  rw [hsub, hsum]
  field_simp
  ring

/-- four coefficients vanish if the quartic is `O(t⁵)` on `(0, d]` -/
lemma quartic_zero {d1 d2 d3 d4 K d : ℝ} (hd : 0 < d) (hd1 : d ≤ 1) (hK : 0 ≤ K)
    (h : ∀ t : ℝ, 0 < t → t ≤ d → |d1 * t + d2 * t ^ 2 + d3 * t ^ 3 + d4 * t ^ 4| ≤ K * t ^ 5) :
    d1 = 0 ∧ d2 = 0 ∧ d3 = 0 ∧ d4 = 0 := by
  have h1 : d1 = 0 := by
    refine Y04.eq_zero_of_le_lin hd (K := |d2| + |d3| + |d4| + K)
      (by have := abs_nonneg d2; have := abs_nonneg d3; have := abs_nonneg d4; linarith)
      (fun t ht0 htd => ?_)
    have ht1 : t ≤ 1 := htd.trans hd1
    have := h t ht0 htd
    have e : d1 = (d1 * t + d2 * t ^ 2 + d3 * t ^ 3 + d4 * t ^ 4) / t - d2 * t - d3 * t ^ 2
        - d4 * t ^ 3 := by
      field_simp; ring
    rw [e]
    have ht2 : t ^ 2 ≤ t := by nlinarith
    have ht3 : t ^ 3 ≤ t := by nlinarith [pow_le_one₀ ht0.le ht1 (n := 2)]
    have ht4 : t ^ 4 ≤ t := by nlinarith [pow_le_one₀ ht0.le ht1 (n := 3)]
    have b1 : |(d1 * t + d2 * t ^ 2 + d3 * t ^ 3 + d4 * t ^ 4) / t| ≤ K * t ^ 4 := by
      rw [abs_div, abs_of_pos ht0, div_le_iff₀ ht0]; nlinarith
    have b2 : |d2 * t| ≤ |d2| * t := by rw [abs_mul, abs_of_pos ht0]
    have b3 : |d3 * t ^ 2| ≤ |d3| * t := by
      rw [abs_mul, abs_of_pos (pow_pos ht0 2)]
      exact mul_le_mul_of_nonneg_left ht2 (abs_nonneg _)
    have b4 : |d4 * t ^ 3| ≤ |d4| * t := by
      rw [abs_mul, abs_of_pos (pow_pos ht0 3)]
      exact mul_le_mul_of_nonneg_left ht3 (abs_nonneg _)
    have t4 : K * t ^ 4 ≤ K * t := mul_le_mul_of_nonneg_left ht4 hK
    have := abs_sub ((d1 * t + d2 * t ^ 2 + d3 * t ^ 3 + d4 * t ^ 4) / t - d2 * t - d3 * t ^ 2)
      (d4 * t ^ 3)
    have := abs_sub ((d1 * t + d2 * t ^ 2 + d3 * t ^ 3 + d4 * t ^ 4) / t - d2 * t) (d3 * t ^ 2)
    have := abs_sub ((d1 * t + d2 * t ^ 2 + d3 * t ^ 3 + d4 * t ^ 4) / t) (d2 * t)
    have : (|d2| + |d3| + |d4| + K) * t = |d2| * t + |d3| * t + |d4| * t + K * t := by ring
    linarith
  subst h1
  obtain ⟨a, b, c⟩ := Y04.cubic_zero hd hd1 hK (a1 := d2) (a2 := d3) (a3 := d4) (fun t ht0 htd => by
    have := h t ht0 htd
    have e : 0 * t + d2 * t ^ 2 + d3 * t ^ 3 + d4 * t ^ 4 = t * (d2 * t + d3 * t ^ 2 + d4 * t ^ 3) := by
      ring
    rw [e, abs_mul, abs_of_pos ht0, show K * t ^ 5 = t * (K * t ^ 4) by ring] at this
    exact le_of_mul_le_mul_left this ht0)
  exact ⟨rfl, a, b, c⟩

end Z02


namespace Z02

/-- the `x⁴` coefficient of the regularised Beta identity at `a = αx`, `b = βx` -/
theorem coef4 {α β : ℝ} (hα : 0 < α) (hβ : 0 < β) (h3 : α + β ≤ 3) :
    (α + β) * (α ^ 3 + β ^ 3) / 6 * Y04.I 3 1 + (α + β) ^ 2 * (α ^ 2 + β ^ 2) / 4 * Y04.I 2 2
      + (α + β) ^ 3 * (α + β) / 6 * Y04.I 1 3 + (α + β) ^ 4 / 12 * Y04.I 0 4
      = (α + β) / (α * β) * 32 * (Y04.cf (α / 2) (β / 2) ((α + β) / 2) 3
          + Y04.cf (α / 2) (β / 2) ((α + β) / 2) 0 * Y04.cf (α / 2) (β / 2) ((α + β) / 2) 1) := by
  set γ := α + β with hγ
  set c0 := Y04.cf (α / 2) (β / 2) (γ / 2) 0 with hc0
  set c1 := Y04.cf (α / 2) (β / 2) (γ / 2) 1 with hc1
  set c2 := Y04.cf (α / 2) (β / 2) (γ / 2) 2 with hc2
  set c3 := Y04.cf (α / 2) (β / 2) (γ / 2) 3 with hc3
  set KB := ∫ s in Set.Ioi (0 : ℝ), B6 s with hKB
  have hKB0 : 0 ≤ KB := setIntegral_nonneg measurableSet_Ioi (fun s hs => B6_nonneg (le_of_lt hs))
  set K0 : ℝ := 69 + 69 * (8 * 23 + 69) + 7 * 23 ^ 2 + (4 * 23 + 69) ^ 3 with hK0
  have hK0p : 0 ≤ K0 := by rw [hK0]; norm_num
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
  have hmain : ∀ x : ℝ, 0 < x → x ≤ 1 / 40 →
      |d1 * x + d2 * x ^ 2 + d3 * x ^ 3 + d4 * x ^ 4| ≤ (20 * 243 * KB + e * 64 * K0) * x ^ 5 := by
    intro x hx0 hx1
    have hx3 : 3 * x ≤ 1 := by linarith
    have gx : γ * x ≤ 3 * x := mul_le_mul_of_nonneg_right (by rw [hγ]; exact h3) hx0.le
    have ax : α * x ≤ 3 * x := mul_le_mul_of_nonneg_right hα3 hx0.le
    have bx : β * x ≤ 3 * x := mul_le_mul_of_nonneg_right hβ3 hx0.le
    have dx1 : (α - β) * x ≤ 3 * x := mul_le_mul_of_nonneg_right (by linarith) hx0.le
    have dx2 : -3 * x ≤ (α - β) * x := mul_le_mul_of_nonneg_right (by linarith) hx0.le
    obtain ⟨_, E1⟩ := int_expand (u := γ * x) (v := α * x) (t := 3 * x) (by positivity)
      gx (by positivity) ax hx3
    obtain ⟨_, E2⟩ := int_expand (u := γ * x) (v := β * x) (t := 3 * x) (by positivity)
      gx (by positivity) bx hx3
    have hC := lhs_closed (a := α * x) (b := β * x) (by positivity) (by positivity)
    rw [show α * x + β * x = γ * x by rw [hγ]; ring] at hC
    have hη : |γ * x| < 1 / 2 := by
      rw [abs_of_pos (by positivity)]; linarith
    have hs : |(α - β) * x / 2| < 1 / 2 := by
      rw [abs_lt]; constructor <;> linarith
    obtain ⟨hpos, h⟩ := Y04.gamma_log_series (γ * x) ((α - β) * x / 2) hη hs
    rw [show γ * x / 2 + (α - β) * x / 2 = α / 2 * (2 * x) by rw [hγ]; ring,
      show γ * x / 2 - (α - β) * x / 2 = β / 2 * (2 * x) by rw [hγ]; ring] at h hpos
    rw [show γ * x = γ / 2 * (2 * x) by ring] at h hpos
    have h2x : |2 * x| ≤ 1 / 20 := by rw [abs_of_pos (by positivity)]; linarith
    have hA := Y04.A_expand (α := α / 2) (β := β / 2) (γ := γ / 2)
      (by rw [abs_of_pos (by positivity)]; linarith) (by rw [abs_of_pos (by positivity)]; linarith)
      (by rw [abs_of_pos (by positivity)]; linarith) h2x h
    rw [Real.exp_log hpos] at hA
    rw [show α / 2 * (2 * x) = α * x by ring, show β / 2 * (2 * x) = β * x by ring,
      show γ / 2 * (2 * x) = γ * x by ring] at hA
    rw [← hc0, ← hc1, ← hc2, ← hc3, ← hK0] at hA
    set Rg := Real.Gamma (1 + α * x) * Real.Gamma (1 + β * x) / Real.Gamma (1 + γ * x) with hRg
    set Iα := ∫ s in Set.Ioi (0 : ℝ), (Real.exp (-(γ * x) * Y04.ell s) - 1) * Real.exp (-(α * x) * s)
    set Iβ := ∫ s in Set.Ioi (0 : ℝ), (Real.exp (-(γ * x) * Y04.ell s) - 1) * Real.exp (-(β * x) * s)
    set Sα := (-(γ * x)) * Y04.I 0 1 + (γ * x * (α * x)) * Y04.I 1 1
        + (-(γ * x * (α * x) ^ 2 / 2)) * Y04.I 2 1
        + (γ * x * (α * x) ^ 3 / 6) * Y04.I 3 1 + ((γ * x) ^ 2 / 2) * Y04.I 0 2
        + (-((γ * x) ^ 2 * (α * x) / 2)) * Y04.I 1 2
        + ((γ * x) ^ 2 * (α * x) ^ 2 / 4) * Y04.I 2 2 + (-((γ * x) ^ 3 / 6)) * Y04.I 0 3
        + ((γ * x) ^ 3 * (α * x) / 6) * Y04.I 1 3 + ((γ * x) ^ 4 / 24) * Y04.I 0 4 with hSα
    set Sβ := (-(γ * x)) * Y04.I 0 1 + (γ * x * (β * x)) * Y04.I 1 1
        + (-(γ * x * (β * x) ^ 2 / 2)) * Y04.I 2 1
        + (γ * x * (β * x) ^ 3 / 6) * Y04.I 3 1 + ((γ * x) ^ 2 / 2) * Y04.I 0 2
        + (-((γ * x) ^ 2 * (β * x) / 2)) * Y04.I 1 2
        + ((γ * x) ^ 2 * (β * x) ^ 2 / 4) * Y04.I 2 2 + (-((γ * x) ^ 3 / 6)) * Y04.I 0 3
        + ((γ * x) ^ 3 * (β * x) / 6) * Y04.I 1 3 + ((γ * x) ^ 4 / 24) * Y04.I 0 4 with hSβ
    set Q := 1 + c0 * (2 * x) ^ 2 + c1 * (2 * x) ^ 3 + (c2 + c0 ^ 2 / 2) * (2 * x) ^ 4
        + (c3 + c0 * c1) * (2 * x) ^ 5 with hQ
    have hpoly : d1 * x + d2 * x ^ 2 + d3 * x ^ 3 + d4 * x ^ 4
        = (Sα - Iα) + (Sβ - Iβ) + ((Iα + Iβ) - e * (4 * c0 * x + 8 * c1 * x ^ 2
          + 16 * (c2 + c0 ^ 2 / 2) * x ^ 3 + 32 * (c3 + c0 * c1) * x ^ 4)) := by
      rw [hd1, hd2, hd3, hd4, hSα, hSβ]; ring
    have hlast : (Iα + Iβ) - e * (4 * c0 * x + 8 * c1 * x ^ 2
          + 16 * (c2 + c0 ^ 2 / 2) * x ^ 3 + 32 * (c3 + c0 * c1) * x ^ 4)
        = e / x * (Rg - Q) := by
      rw [hC, he, hQ]
      field_simp
      ring
    rw [hpoly, hlast]
    have b3 : |e / x * (Rg - Q)| ≤ e * 64 * K0 * x ^ 5 := by
      rw [abs_mul, abs_of_nonneg (by positivity)]
      calc e / x * |Rg - Q| ≤ e / x * (K0 * (2 * x) ^ 6) :=
            mul_le_mul_of_nonneg_left hA (by positivity)
        _ = e * 64 * K0 * x ^ 5 := by field_simp; ring
    have b1 : |Sα - Iα| ≤ 10 * (3 * x) ^ 5 * KB := by rw [abs_sub_comm]; exact E1
    have b2 : |Sβ - Iβ| ≤ 10 * (3 * x) ^ 5 * KB := by rw [abs_sub_comm]; exact E2
    have t1 := abs_add_le ((Sα - Iα) + (Sβ - Iβ)) (e / x * (Rg - Q))
    have t2 := abs_add_le (Sα - Iα) (Sβ - Iβ)
    have : (20 * 243 * KB + e * 64 * K0) * x ^ 5
        = 10 * (3 * x) ^ 5 * KB + 10 * (3 * x) ^ 5 * KB + e * 64 * K0 * x ^ 5 := by ring
    linarith
  obtain ⟨_, _, _, h4⟩ := quartic_zero (by norm_num) (by norm_num) (by positivity) hmain
  rw [hd4] at h4
  linarith

end Z02


namespace Z02

lemma eqA : 2 / 3 * Y04.I 3 1 + 2 * Y04.I 2 2 + 8 / 3 * Y04.I 1 3 + 4 / 3 * Y04.I 0 4
    = 12 * z5 - 4 * ((π ^ 2 / 6) * Y04.z3) := by
  have h := coef4 (α := 1) (β := 1) one_pos one_pos (by norm_num)
  norm_num [Y04.cf, Y04.z_two, Y04.zn_three, zn_five] at h
  linear_combination h

lemma eqB : 9 / 2 * Y04.I 3 1 + 45 / 4 * Y04.I 2 2 + 27 / 2 * Y04.I 1 3 + 27 / 4 * Y04.I 0 4
    = 63 * z5 - 18 * ((π ^ 2 / 6) * Y04.z3) := by
  have h := coef4 (α := 1) (β := 2) one_pos two_pos (by norm_num)
  norm_num [Y04.cf, Y04.z_two, Y04.zn_three, zn_five] at h
  linear_combination h

lemma D51 : Y04.I 3 1 + Y04.I 2 2 = 2 * z5 + 2 * (π ^ 2 / 6) * Y04.z3 := by
  linear_combination (8 / 9) * eqB - (9 / 2) * eqA

lemma D50 : 2 * (-(17 / 48) * Y04.I 3 1 + (-(5 / 8)) * Y04.I 2 2 / 2 + (1 / 2) * Y04.I 1 3 / 6
      + Y04.I 0 4 / 24) =
    -(3 / 4) * z5 - 7 / 4 * ((π ^ 2 / 6) * Y04.z3) := by
  linear_combination (1 / 16) * eqA - (3 / 4) * D51

end Z02


namespace Z02

lemma G_d3' (k : ℕ) : iteratedDeriv 3 (Y04.G k) 0 =
    6 * -((-1 : ℝ) ^ k * Real.sqrt 2 / 4 * (17 / 48 + (k : ℝ) ^ 2)) := by
  have := G_d3 k; linarith

lemma G_d2' (k : ℕ) : iteratedDeriv 2 (Y04.G k) 0 =
    2 * -((-1 : ℝ) ^ k * Real.sqrt 2 * ((k : ℝ) ^ 2 / 2 + 5 / 32)) := by
  have := Y04.G_d2 k; linarith

lemma theta5 (k : ℕ) :
    (-1 : ℝ) ^ k * (π / Real.sqrt 2) *
        (8 / π * ∑ j ∈ Finset.Icc 1 4,
          iteratedDeriv (4 - j) (Y04.G k) 0 / ((4 - j).factorial : ℝ) / (j.factorial : ℝ) *
            Y04.I (4 - j) j) =
      -(3 / 4 + 4 * (k : ℝ) ^ 2) * z5 - (7 / 4 + 4 * (k : ℝ) ^ 2) * ((π ^ 2 / 6) * Y04.z3) := by
  set E := (-1 : ℝ) ^ k with hE
  set Y := -(17 / 48 + (k : ℝ) ^ 2) * Y04.I 3 1 + (-(2 * (k : ℝ) ^ 2 + 5 / 8)) / 2 * Y04.I 2 2
    + (1 / 2) / 6 * Y04.I 1 3 + 1 / 24 * Y04.I 0 4 with hY
  have hsum : ∑ j ∈ Finset.Icc 1 4,
          iteratedDeriv (4 - j) (Y04.G k) 0 / ((4 - j).factorial : ℝ) / (j.factorial : ℝ) *
            Y04.I (4 - j) j = E * Real.sqrt 2 / 4 * Y := by
    rw [show Finset.Icc 1 4 = {1, 2, 3, 4} by decide]
    rw [Finset.sum_insert (by decide), Finset.sum_insert (by decide), Finset.sum_insert (by decide),
      Finset.sum_singleton]
    simp only [show 4 - 1 = 3 from rfl, show 4 - 2 = 2 from rfl, show 4 - 3 = 1 from rfl,
      show 4 - 4 = 0 from rfl, iteratedDeriv_zero, iteratedDeriv_one]
    rw [G_d3', G_d2', Y04.G_d1, Y04.G_zero, hY]
    simp [Nat.factorial]
    ring
  rw [hsum]
  have hE2 : E ^ 2 = 1 := by rw [hE, ← pow_mul, mul_comm, pow_mul, neg_one_sq, one_pow]
  have hs : Real.sqrt 2 ≠ 0 := by positivity
  have hpi : π ≠ 0 := Real.pi_pos.ne'
  have e : E * (π / Real.sqrt 2) * (8 / π * (E * Real.sqrt 2 / 4 * Y)) = 2 * E ^ 2 * Y := by
    field_simp; ring
  rw [e, hE2, hY]
  linear_combination D50 - 2 * (k : ℝ) ^ 2 * D51

lemma sum_eval5 (Gf : ℝ → ℝ) (τ : ℝ) :
    ∑ n ∈ Finset.Icc 2 5, τ ^ n * ∑ j ∈ Finset.Icc 1 (n - 1),
        iteratedDeriv (n - 1 - j) Gf 0 / ((n - 1 - j).factorial : ℝ) / (j.factorial : ℝ) * Y04.I (n - 1 - j) j
      = τ ^ 2 * Gf 0 * Y04.I 0 1 + τ ^ 3 * deriv Gf 0 * Y04.I 1 1 + τ ^ 3 * (Gf 0 / 2) * Y04.I 0 2
        + τ ^ 4 * (iteratedDeriv 2 Gf 0 / 2) * Y04.I 2 1
        + τ ^ 4 * (deriv Gf 0 / 2) * Y04.I 1 2 + τ ^ 4 * (Gf 0 / 6) * Y04.I 0 3
        + τ ^ 5 * (iteratedDeriv 3 Gf 0 / 6) * Y04.I 3 1 + τ ^ 5 * (iteratedDeriv 2 Gf 0 / 4) * Y04.I 2 2
        + τ ^ 5 * (deriv Gf 0 / 6) * Y04.I 1 3 + τ ^ 5 * (Gf 0 / 24) * Y04.I 0 4 := by
  rw [show Finset.Icc 2 5 = {2, 3, 4, 5} by decide]
  rw [Finset.sum_insert (by decide), Finset.sum_insert (by decide), Finset.sum_insert (by decide),
    Finset.sum_singleton]
  rw [show Finset.Icc 1 (2 - 1) = {1} by decide, show Finset.Icc 1 (3 - 1) = {1, 2} by decide,
    show Finset.Icc 1 (4 - 1) = {1, 2, 3} by decide, show Finset.Icc 1 (5 - 1) = {1, 2, 3, 4} by decide]
  simp [Finset.sum_insert, Nat.factorial]
  ring

lemma fourier5 (k : ℕ) : ∃ C τ₀ : ℝ, 0 < τ₀ ∧ ∀ p : ℝ, 1 < p → Y04.tau p ≤ τ₀ →
    |Y04.a k (Y04.h p) - (-1 : ℝ) ^ k * (Real.sqrt 2 / π) *
        (-(4 / (16 * (k : ℝ) ^ 2 - 1)) + (π ^ 2 / 6) * Y04.tau p ^ 2 + Y04.z3 * Y04.tau p ^ 3
          - (33 / 16 + 7 * (k : ℝ) ^ 2) * (π ^ 4 / 90) * Y04.tau p ^ 4
          + (-(3 / 4 + 4 * (k : ℝ) ^ 2) * z5
              - (7 / 4 + 4 * (k : ℝ) ^ 2) * ((π ^ 2 / 6) * Y04.z3)) * Y04.tau p ^ 5)|
      ≤ C * Y04.tau p ^ 6 := by
  obtain ⟨CW, hW⟩ := watson5 (Y04.G k) (Y04.G_contDiff k) ⟨1, fun x _ => Y04.G_bdd k x⟩
  refine ⟨8 / π * CW, 1 / 2, by norm_num, fun p hp _ => ?_⟩
  have hpi := Real.pi_pos
  have hτ0 := Y04.tau_pos hp
  have hτ1 := Y04.tau_lt_one hp
  obtain ⟨_, hb⟩ := hW (Y04.tau p) hτ0 hτ1
  rw [sum_eval5] at hb
  have h2 := Y04.layer_formula k hp
  have h1 := Y04.A_formula k
  set X := Y04.tau p * ∫ s in Set.Ioi (0 : ℝ), ((1 + Real.exp (-s)) ^ Y04.tau p - 1) * Y04.G k (Y04.tau p * s)
    with hX
  set τ := Y04.tau p with hτ
  set S := τ ^ 2 * Y04.G k 0 * Y04.I 0 1 + τ ^ 3 * deriv (Y04.G k) 0 * Y04.I 1 1
        + τ ^ 3 * (Y04.G k 0 / 2) * Y04.I 0 2
        + τ ^ 4 * (iteratedDeriv 2 (Y04.G k) 0 / 2) * Y04.I 2 1
        + τ ^ 4 * (deriv (Y04.G k) 0 / 2) * Y04.I 1 2 + τ ^ 4 * (Y04.G k 0 / 6) * Y04.I 0 3
        + τ ^ 5 * (iteratedDeriv 3 (Y04.G k) 0 / 6) * Y04.I 3 1
        + τ ^ 5 * (iteratedDeriv 2 (Y04.G k) 0 / 4) * Y04.I 2 2
        + τ ^ 5 * (deriv (Y04.G k) 0 / 6) * Y04.I 1 3 + τ ^ 5 * (Y04.G k 0 / 24) * Y04.I 0 4 with hS
  have ha : Y04.a k (Y04.h p) = Y04.a k Y04.hDiamond + 8 / π * X := by rw [hX, ← mul_assoc, ← h2]; ring
  have a1 := Y04.four_k_ne k; have b1 := Y04.four_k_ne' k
  have c1 : (16 * (k : ℝ) ^ 2 - 1) ≠ 0 := by
    rw [show (16 * (k : ℝ) ^ 2 - 1) = (4 * k - 1) * (4 * k + 1) by ring]; exact mul_ne_zero a1 b1
  have hI03 : Y04.I 0 3 = (3 / 4 * (π ^ 4 / 90) - 3 * Y04.I 1 2) / 2 := by linarith [Y04.I0312]
  have hI22 : Y04.I 2 2 = 2 * z5 + 2 * (π ^ 2 / 6) * Y04.z3 - 45 / 8 * z5 := by
    have := D51; rw [I31] at this; linarith
  have hI13 : Y04.I 1 3 = 6 * (-(3 / 4) * z5 - 7 / 4 * ((π ^ 2 / 6) * Y04.z3))
      + 17 / 4 * (45 / 8 * z5) + 15 / 4 * (2 * z5 + 2 * (π ^ 2 / 6) * Y04.z3 - 45 / 8 * z5)
      - Y04.I 0 4 / 2 := by
    have := D50; rw [I31, hI22] at this; linarith
  have key : Y04.a k (Y04.h p) - (-1 : ℝ) ^ k * (Real.sqrt 2 / π) *
        (-(4 / (16 * (k : ℝ) ^ 2 - 1)) + (π ^ 2 / 6) * τ ^ 2 + Y04.z3 * τ ^ 3
          - (33 / 16 + 7 * (k : ℝ) ^ 2) * (π ^ 4 / 90) * τ ^ 4
          + (-(3 / 4 + 4 * (k : ℝ) ^ 2) * z5
              - (7 / 4 + 4 * (k : ℝ) ^ 2) * ((π ^ 2 / 6) * Y04.z3)) * τ ^ 5) = 8 / π * (X - S) := by
    have hd2 : iteratedDeriv 2 (Y04.G k) 0 = 2 * -((-1 : ℝ) ^ k * Real.sqrt 2 * ((k : ℝ) ^ 2 / 2 + 5 / 32)) :=
      G_d2' k
    rw [ha, h1, hS, hd2, G_d3', Y04.G_zero k, Y04.G_d1 k, Y04.I01, Y04.I11, Y04.I02, Y04.I21, hI03,
      I31, hI22, hI13]
    field_simp
    ring
  rw [key, abs_mul, abs_of_pos (by positivity)]
  calc 8 / π * |X - S| ≤ 8 / π * (CW * τ ^ 6) := mul_le_mul_of_nonneg_left hb (by positivity)
    _ = _ := by ring

lemma zn_pos {n : ℕ} (hn : 1 < n) : 0 < Y04.zn n := by
  have h := Y04.z_hasSum hn
  have := le_hasSum h 0 (fun j _ => by positivity)
  simp at this; linarith

lemma z3_pos : 0 < Y04.z3 := by rw [← Y04.zn_three]; exact zn_pos (by norm_num)

lemma z5_pos : 0 < z5 := by rw [← zn_five]; exact zn_pos (by norm_num)

end Z02



theorem Z02.proof_1 : Z02.stmt_1 := Z02.G_d3

theorem Z02.proof_2 : Z02.stmt_2 := fun Gf hG hb => Z02.watson5 Gf hG hb

theorem Z02.proof_3 : Z02.stmt_3 := Z02.I31

theorem Z02.proof_4 : Z02.stmt_4 := Z02.D51

theorem Z02.proof_5 : Z02.stmt_5 := Z02.D50

theorem Z02.proof_6 : Z02.stmt_6 := Z02.theta5

theorem Z02.proof_7 : Z02.stmt_7 := by
  have hpi := Real.pi_pos
  have h3 := Z02.z3_pos
  have h5 := Z02.z5_pos
  have n2 : (π ^ 2 / 6 : ℝ) ≠ 0 := by positivity
  have n4 : (-(7 * (π ^ 4 / 90)) : ℝ) ≠ 0 := by
    have : 0 < π ^ 4 / 90 := by positivity
    linarith
  have n5 : (-4 * Z02.z5 - 4 * ((π ^ 2 / 6) * Y04.z3) : ℝ) ≠ 0 := by
    have : 0 < (π ^ 2 / 6) * Y04.z3 := by positivity
    linarith
  refine ⟨Polynomial.C (π ^ 2 / 6), Polynomial.C Y04.z3,
    Polynomial.C (-(7 * (π ^ 4 / 90))) * Polynomial.X + Polynomial.C (-(33 / 16) * (π ^ 4 / 90)),
    Polynomial.C (-4 * Z02.z5 - 4 * ((π ^ 2 / 6) * Y04.z3)) * Polynomial.X
      + Polynomial.C (-(3 / 4) * Z02.z5 - 7 / 4 * ((π ^ 2 / 6) * Y04.z3)), ?_, ?_, ?_, ?_, ?_, ?_, ?_⟩
  · rw [Polynomial.degree_C n2]; rfl
  · rw [Polynomial.degree_C h3.ne']; rfl
  · rw [Polynomial.degree_linear n4]; rfl
  · rw [Polynomial.degree_linear n5]; rfl
  · rw [Polynomial.leadingCoeff_C]; norm_num
  · rw [Polynomial.leadingCoeff_linear n4]; norm_num
  · intro k
    obtain ⟨C, τ₀, hτ₀, hC⟩ := Z02.fourier5 k
    refine ⟨C, τ₀, hτ₀, fun p hp hτ => ?_⟩
    have := hC p hp hτ
    simp only [Polynomial.eval_add, Polynomial.eval_mul, Polynomial.eval_C, Polynomial.eval_X]
    convert this using 4
    ring

theorem Z02.proof_8 : Z02.stmt_8 := Z02.fourier5
