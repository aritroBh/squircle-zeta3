import Mathlib
/-!
# Y04 — Fourier coefficients of the squircle's support function at the diamond end

Target id: Y04.
Paper labels: `thm:fourier` (main.tex lines 1204–1211), with the definitions of `h_p` and
`a_k(h_p)` just before it (main.tex line 1202). Detailed proof:
paper/ieee/proofs/sec_g2_fourier.tex (Lemmas `g2:lem:f:red`, `g2:lem:f:g012`,
`g2:lem:f:watson`, `g2:lem:f:gen`, `g2:lem:f:atoms`, and Section `g2:sec:f:proof`).

Quoted verbatim.

main.tex 1202:

  $h_p(\theta)=\bigl(|\cos\theta|^{1/\tau}+|\sin\theta|^{1/\tau}\bigr)^{\tau}$, $\tau=1-\frac1p$.
  It is $\max(|\cos\theta|,|\sin\theta|)$ at the diamond ...
  $a_k(h_p)=\frac4\pi\int_0^{\pi/2}h_p\cos4k\theta\dd\theta$.

main.tex 1204–1211 (`thm:fourier`):

  For every $k\ge0$, with $p>1$ and $\tau=1-\frac1p\to0^+$,
  a_k(h_p)=(-1)^k\frac{\sqrt2}{\pi}\Bigl[-\frac{4}{16k^2-1}+\zeta(2)\tau^2+\zeta(3)\tau^3\\
  -\Bigl(\frac{33}{16}+7k^2\Bigr)\zeta(4)\tau^4+O_k(\tau^5)\Bigr].

Supplement, Lemma `g2:lem:f:red`:

  (a) $A_k=-(-1)^k\dfrac{4\sqrt2}{\pi(16k^2-1)}$.
  (b) $\displaystyle a_k(h_p)-A_k=\frac8\pi\,\tau\int_0^\infty\bigl[(1+\eu^{-s})^\tau-1\bigr]\,
      G_k(\tau s)\,ds .$
  with $\theta_x=\arctan\eu^{-x}$, $G_k(x)=\cos(4k\theta_x)\cos^2\theta_x\sin\theta_x$.

Supplement, Lemma `g2:lem:f:g012`:

  $g_{k,0}=(-1)^k\frac{\sqrt2}4$, $g_{k,1}=(-1)^k\frac{\sqrt2}8$,
  $g_{k,2}=-(-1)^k\sqrt2\bigl(\frac{k^2}2+\frac5{32}\bigr)$.   ($g_{k,m}=G_k^{(m)}(0)/m!$)

Supplement, Lemma `g2:lem:f:watson`:

  \Bigl|\,\tau\int_0^\infty\bigl[(1+\eu^{-s})^\tau-1\bigr]G(\tau s)\,ds
    -\sum_{n=2}^{N}\tau^n\sum_{j=1}^{n-1}\frac{g_{n-1-j}}{j!}\,I_{n-1-j,\,j}\Bigr|\le C\tau^{N+1}.
  with $I_{m,j}=\int_0^\infty s^m\,\ell(s)^j\,ds$, $\ell(s)=\log(1+\eu^{-s})$.

Supplement, Lemma `g2:lem:f:gen`:

  $J_j=\int_0^\infty\bigl(m(s)^j-(\tfrac s2)^j\bigr)ds$, $m(s)=\frac s2+\ell(s)$,
  $R(x)=\Gamma(1+x)^2/\Gamma(1+2x)$, ... $\Psi(\tau)=\int_0^\infty\eu^{\tau s/2}
  \bigl(\eu^{\tau\ell(s)}-1\bigr)ds$ ... $\Psi=\Phi$, $\Phi(\tau)=(R(-\tau/2)-1)/(-\tau/2)$ ...
  $J_1=\tfrac12\zeta(2),\qquad J_2=\zeta(3),\qquad J_3=\tfrac{27}{16}\zeta(4).$

Supplement, Lemma `g2:lem:f:atoms`:

  $I_{0,1}=\frac12\zeta(2)$, $I_{1,1}=\frac34\zeta(3)$, $I_{2,1}=\frac74\zeta(4)$,
  $I_{0,2}=\frac14\zeta(3)$ and $2I_{0,3}+3I_{1,2}=\frac34\zeta(4)$.

Conventions.
* `ζ(3)` is `(riemannZeta 3).re`, the real part of Mathlib's complex `riemannZeta 3`; by
  `zeta_nat_eq_tsum_of_gt_one` it is the real number `∑ 1/n^3`. `ζ(2) = π²/6`, `ζ(4) = π⁴/90`
  are written out.
* `h_p` is defined for `p > 1` through `τ = 1 - 1/p ∈ (0,1)`; all bases of `rpow` are `≥ 0`
  and the sum `|cos θ|^{1/τ} + |sin θ|^{1/τ}` is `> 0`, so `h_p` is continuous and the
  interval integral defining `a_k(h_p)` is a genuine integral. The diamond `h_1 = max(|cos|,|sin|)`
  is defined separately (the formula with `τ = 0` is not meaningful).
* Half-line integrals are Bochner integrals over `Ioi 0`; every statement using one also
  asserts integrability, so no statement holds by the "non-integrable ↦ 0" convention.
* In the expansion lemma (`stmt_5`) we ask `G` to be `C^∞` on all of `ℝ` and bounded on
  `[0,∞)` (the supplement asks `C^∞` on `[0,∞)`); this covers `G = G_k`, which is smooth
  on `ℝ` (Lemma `g2:lem:f:red`(c)). It is a restriction of the lemma, not of the theorem.
* In `stmt_6` the closed form of `Ψ` is stated for real `τ < 1`, `τ ≠ 0` (the supplement's
  holomorphic identity on `Re τ < 2`, restricted to a real interval where every `Γ` value is
  away from its poles); the values `J_1, J_2, J_3` are stated separately.

NOT COVERED: part (c) of Lemma `g2:lem:f:red` beyond what `stmt_3` records (the Gudermannian
form `G_k = (-1)^k cos(2k gd x) H(x)`); the holomorphy of `Φ`, `Ψ` in the complex variable;
the `log R` power series (it is the `log Γ(1+x)` series used in X05); Cauchy's formula
`P = ∫ h_p` and the `k = 0` consistency check with `thm:pw`; `prop:degree`, `prop:theta`.
-/

-- BEGIN STATEMENTS
namespace Y04

open Real MeasureTheory Set

/-- `ζ(3)` as a real number: the real part of Mathlib's `riemannZeta 3`. -/
noncomputable def z3 : ℝ := (riemannZeta 3).re

/-- `τ = 1 - 1/p`. -/
noncomputable def tau (p : ℝ) : ℝ := 1 - 1 / p

/-- The support function of the squircle, `h_p(θ) = (|cos θ|^{1/τ} + |sin θ|^{1/τ})^τ`,
`τ = 1 - 1/p` (main.tex 1202), for `p > 1`. -/
noncomputable def h (p θ : ℝ) : ℝ :=
  (|Real.cos θ| ^ (1 / tau p) + |Real.sin θ| ^ (1 / tau p)) ^ tau p

/-- The support function of the diamond, `h_1(θ) = max(|cos θ|, |sin θ|)`. -/
noncomputable def hDiamond (θ : ℝ) : ℝ := max |Real.cos θ| |Real.sin θ|

/-- The Fourier coefficient `a_k(f) = (4/π) ∫_0^{π/2} f(θ) cos(4kθ) dθ` (period `π/2`). -/
noncomputable def a (k : ℕ) (f : ℝ → ℝ) : ℝ :=
  4 / π * ∫ θ in (0 : ℝ)..(π / 2), f θ * Real.cos (4 * k * θ)

/-- `ℓ(s) = log(1 + e^{-s})`. -/
noncomputable def ell (s : ℝ) : ℝ := Real.log (1 + Real.exp (-s))

/-- `θ_x = arctan e^{-x}`. -/
noncomputable def thetaX (x : ℝ) : ℝ := Real.arctan (Real.exp (-x))

/-- `G_k(x) = cos(4k θ_x) cos² θ_x sin θ_x`. -/
noncomputable def G (k : ℕ) (x : ℝ) : ℝ :=
  Real.cos (4 * k * thetaX x) * Real.cos (thetaX x) ^ 2 * Real.sin (thetaX x)

/-- `I_{m,j} = ∫_0^∞ s^m ℓ(s)^j ds`. -/
noncomputable def I (m j : ℕ) : ℝ := ∫ s in Ioi (0 : ℝ), s ^ m * ell s ^ j

/-- `J_j = ∫_0^∞ (m(s)^j - (s/2)^j) ds`, `m(s) = s/2 + ℓ(s)`. -/
noncomputable def J (j : ℕ) : ℝ := ∫ s in Ioi (0 : ℝ), ((s / 2 + ell s) ^ j - (s / 2) ^ j)

/-- `R(x) = Γ(1+x)² / Γ(1+2x)`. -/
noncomputable def R (x : ℝ) : ℝ := Real.Gamma (1 + x) ^ 2 / Real.Gamma (1 + 2 * x)

/-- Lemma `g2:lem:f:red`(a): the diamond coefficients
`A_k = a_k(h_1) = -(-1)^k 4√2 / (π(16k²-1))`. -/
def stmt_1 : Prop :=
  ∀ k : ℕ, a k hDiamond = -((-1 : ℝ) ^ k * (4 * Real.sqrt 2) / (π * (16 * (k : ℝ) ^ 2 - 1)))

/-- Lemma `g2:lem:f:red`(b): for `p > 1`, the layer integrand is integrable on `(0,∞)` and
`a_k(h_p) - A_k = (8/π) τ ∫_0^∞ [(1+e^{-s})^τ - 1] G_k(τ s) ds`. -/
def stmt_2 : Prop :=
  ∀ k : ℕ, ∀ p : ℝ, 1 < p →
    IntegrableOn (fun s : ℝ => ((1 + Real.exp (-s)) ^ tau p - 1) * G k (tau p * s)) (Ioi 0) ∧
    a k (h p) - a k hDiamond =
      8 / π * tau p * ∫ s in Ioi (0 : ℝ), ((1 + Real.exp (-s)) ^ tau p - 1) * G k (tau p * s)

/-- Lemmas `g2:lem:f:red`(c) and `g2:lem:f:g012`: `G_k` is `C^∞` on `ℝ`, bounded on `[0,∞)`,
and `G_k(0) = (-1)^k √2/4`, `G_k'(0) = (-1)^k √2/8`,
`G_k''(0)/2 = -(-1)^k √2 (k²/2 + 5/32)`. -/
def stmt_3 : Prop :=
  ∀ k : ℕ, ContDiff ℝ (⊤ : ℕ∞) (G k) ∧ (∃ M : ℝ, ∀ x : ℝ, 0 ≤ x → |G k x| ≤ M) ∧
    G k 0 = (-1 : ℝ) ^ k * Real.sqrt 2 / 4 ∧
    deriv (G k) 0 = (-1 : ℝ) ^ k * Real.sqrt 2 / 8 ∧
    iteratedDeriv 2 (G k) 0 / 2 = -((-1 : ℝ) ^ k * Real.sqrt 2 * ((k : ℝ) ^ 2 / 2 + 5 / 32))

/-- Lemma `g2:lem:f:atoms` (with the finiteness of `I_{m,j}`): every `s^m ℓ(s)^j`, `j ≥ 1`, is
integrable on `(0,∞)`, and `I_{0,1} = ζ(2)/2`, `I_{1,1} = 3ζ(3)/4`, `I_{2,1} = 7ζ(4)/4`,
`I_{0,2} = ζ(3)/4`, `2 I_{0,3} + 3 I_{1,2} = 3ζ(4)/4`. -/
def stmt_4 : Prop :=
  (∀ m j : ℕ, 1 ≤ j → IntegrableOn (fun s : ℝ => s ^ m * ell s ^ j) (Ioi 0)) ∧
    I 0 1 = (π ^ 2 / 6) / 2 ∧ I 1 1 = 3 / 4 * z3 ∧ I 2 1 = 7 / 4 * (π ^ 4 / 90) ∧
    I 0 2 = z3 / 4 ∧ 2 * I 0 3 + 3 * I 1 2 = 3 / 4 * (π ^ 4 / 90)

/-- Lemma `g2:lem:f:watson` with `N = 4`: for `G` smooth and bounded on `[0,∞)`, with
`g_m = G^{(m)}(0)/m!`, there is `C` with, for `0 < τ < 1`,
`|τ ∫_0^∞ [(1+e^{-s})^τ - 1] G(τ s) ds - Σ_{n=2}^{4} τ^n Σ_{j=1}^{n-1} g_{n-1-j}/j! I_{n-1-j,j}|
  ≤ C τ^5`, and the integrand is integrable. -/
def stmt_5 : Prop :=
  ∀ Gf : ℝ → ℝ, ContDiff ℝ (⊤ : ℕ∞) Gf → (∃ M : ℝ, ∀ x : ℝ, 0 ≤ x → |Gf x| ≤ M) →
    ∃ C : ℝ, ∀ τ : ℝ, 0 < τ → τ < 1 →
      IntegrableOn (fun s : ℝ => ((1 + Real.exp (-s)) ^ τ - 1) * Gf (τ * s)) (Ioi 0) ∧
      |τ * (∫ s in Ioi (0 : ℝ), ((1 + Real.exp (-s)) ^ τ - 1) * Gf (τ * s)) -
          ∑ n ∈ Finset.Icc 2 4, τ ^ n * ∑ j ∈ Finset.Icc 1 (n - 1),
            iteratedDeriv (n - 1 - j) Gf 0 / ((n - 1 - j).factorial : ℝ) /
              (j.factorial : ℝ) * I (n - 1 - j) j| ≤ C * τ ^ 5

/-- Lemma `g2:lem:f:gen`, the closed form on a real interval: for real `τ < 1`, `τ ≠ 0`,
`e^{τ s/2}(e^{τ ℓ(s)} - 1)` is integrable on `(0,∞)` and
`Ψ(τ) = ∫_0^∞ e^{τ s/2}(e^{τ ℓ(s)} - 1) ds = (R(-τ/2) - 1)/(-τ/2)`. -/
def stmt_6 : Prop :=
  ∀ τ : ℝ, τ < 1 → τ ≠ 0 →
    IntegrableOn (fun s : ℝ => Real.exp (τ * s / 2) * (Real.exp (τ * ell s) - 1)) (Ioi 0) ∧
    ∫ s in Ioi (0 : ℝ), Real.exp (τ * s / 2) * (Real.exp (τ * ell s) - 1) =
      (R (-τ / 2) - 1) / (-τ / 2)

/-- Lemma `g2:lem:f:gen`, the values: `J_1 = ζ(2)/2`, `J_2 = ζ(3)`, `J_3 = 27ζ(4)/16`
(with the integrands of `J_1, J_2, J_3` integrable on `(0,∞)`). -/
def stmt_7 : Prop :=
  (∀ j : ℕ, 1 ≤ j → j ≤ 3 →
      IntegrableOn (fun s : ℝ => (s / 2 + ell s) ^ j - (s / 2) ^ j) (Ioi 0)) ∧
    J 1 = (π ^ 2 / 6) / 2 ∧ J 2 = z3 ∧ J 3 = 27 / 16 * (π ^ 4 / 90)

/-- `thm:fourier` (main.tex 1204–1211), with the `O_k(τ⁵)` made explicit: for every `k ≥ 0`
there are `C_k` and `τ_0 > 0` such that for every `p > 1` with `τ = 1 - 1/p ≤ τ_0`,
`|a_k(h_p) - (-1)^k (√2/π)[-4/(16k²-1) + ζ(2)τ² + ζ(3)τ³ - (33/16 + 7k²)ζ(4)τ⁴]| ≤ C_k τ⁵`. -/
def stmt_8 : Prop :=
  ∀ k : ℕ, ∃ C τ₀ : ℝ, 0 < τ₀ ∧ ∀ p : ℝ, 1 < p → tau p ≤ τ₀ →
    |a k (h p) - (-1 : ℝ) ^ k * (Real.sqrt 2 / π) *
        (-(4 / (16 * (k : ℝ) ^ 2 - 1)) + (π ^ 2 / 6) * tau p ^ 2 + z3 * tau p ^ 3
          - (33 / 16 + 7 * (k : ℝ) ^ 2) * (π ^ 4 / 90) * tau p ^ 4)|
      ≤ C * tau p ^ 5

end Y04
-- END STATEMENTS

open Real MeasureTheory Set

namespace Y04

lemma four_k_ne (k : ℕ) : (4 * (k : ℝ) - 1) ≠ 0 := by
  intro h
  have : (4 * k : ℝ) = 1 := by linarith
  have h2 : (4 * k : ℕ) = 1 := by exact_mod_cast this
  omega

lemma four_k_ne' (k : ℕ) : (4 * (k : ℝ) + 1) ≠ 0 := by positivity

lemma hd_cos (k : ℕ) (θ : ℝ) :
    HasDerivAt (fun θ => Real.sin ((4 * k + 1) * θ) / (2 * (4 * k + 1))
      + Real.sin ((4 * k - 1) * θ) / (2 * (4 * k - 1)))
      (Real.cos θ * Real.cos (4 * k * θ)) θ := by
  have h1 := ((hasDerivAt_id θ).const_mul (4 * (k : ℝ) + 1)).sin.div_const (2 * (4 * (k : ℝ) + 1))
  have h2 := ((hasDerivAt_id θ).const_mul (4 * (k : ℝ) - 1)).sin.div_const (2 * (4 * (k : ℝ) - 1))
  convert h1.add h2 using 1
  · funext x; simp
  have a := four_k_ne k; have b := four_k_ne' k
  simp only [id]
  rw [show (4 * (k : ℝ) + 1) * θ = 4 * k * θ + θ by ring,
    show (4 * (k : ℝ) - 1) * θ = 4 * k * θ - θ by ring, Real.cos_add, Real.cos_sub]
  field_simp; ring

lemma hd_sin (k : ℕ) (θ : ℝ) :
    HasDerivAt (fun θ => -Real.cos ((4 * k + 1) * θ) / (2 * (4 * k + 1))
      + Real.cos ((4 * k - 1) * θ) / (2 * (4 * k - 1)))
      (Real.sin θ * Real.cos (4 * k * θ)) θ := by
  have h1 := ((hasDerivAt_id θ).const_mul (4 * (k : ℝ) + 1)).cos.neg.div_const (2 * (4 * (k : ℝ) + 1))
  have h2 := ((hasDerivAt_id θ).const_mul (4 * (k : ℝ) - 1)).cos.div_const (2 * (4 * (k : ℝ) - 1))
  convert h1.add h2 using 1
  · funext x; simp
  have a := four_k_ne k; have b := four_k_ne' k
  simp only [id, mul_one]
  rw [show (4 * (k : ℝ) + 1) * θ = 4 * k * θ + θ by ring,
    show (4 * (k : ℝ) - 1) * θ = 4 * k * θ - θ by ring, Real.sin_add, Real.sin_sub]
  field_simp; ring

lemma sin_le_cos_of {θ : ℝ} (h0 : 0 ≤ θ) (h1 : θ ≤ π / 4) : Real.sin θ ≤ Real.cos θ := by
  rw [← Real.cos_pi_div_two_sub]
  exact Real.cos_le_cos_of_nonneg_of_le_pi h0 (by linarith [Real.pi_pos]) (by linarith)

lemma cos_le_sin_of {θ : ℝ} (h0 : π / 4 ≤ θ) (h1 : θ ≤ π / 2) : Real.cos θ ≤ Real.sin θ := by
  rw [← Real.cos_pi_div_two_sub]
  exact Real.cos_le_cos_of_nonneg_of_le_pi (by linarith) (by linarith [Real.pi_pos]) (by linarith)

lemma A_formula (k : ℕ) : a k hDiamond = -((-1 : ℝ) ^ k * (4 * Real.sqrt 2) / (π * (16 * (k : ℝ) ^ 2 - 1))) := by
  have hpi := Real.pi_pos
  have hc : ContinuousOn (fun θ => hDiamond θ * Real.cos (4 * k * θ)) (uIcc 0 (π / 2)) := by
    unfold hDiamond; fun_prop
  have hm : π / 4 ∈ uIcc 0 (π / 2) := by
    rw [uIcc_of_le (by linarith)]; constructor <;> linarith
  have hsplit := intervalIntegral.integral_add_adjacent_intervals
    (μ := volume) (a := 0) (b := π / 4) (c := π / 2)
    (hc.mono (uIcc_subset_uIcc left_mem_uIcc hm)).intervalIntegrable
    (hc.mono (uIcc_subset_uIcc hm right_mem_uIcc)).intervalIntegrable
  have e1 : ∫ θ in (0:ℝ)..(π / 4), hDiamond θ * Real.cos (4 * k * θ)
      = ∫ θ in (0:ℝ)..(π / 4), Real.cos θ * Real.cos (4 * k * θ) := by
    refine intervalIntegral.integral_congr (fun θ hθ => ?_)
    rw [uIcc_of_le (by linarith)] at hθ
    have hs : 0 ≤ Real.sin θ := Real.sin_nonneg_of_nonneg_of_le_pi hθ.1 (by linarith [hθ.2])
    have hcz : 0 ≤ Real.cos θ := Real.cos_nonneg_of_mem_Icc ⟨by linarith [hθ.1], by linarith [hθ.2]⟩
    simp only [hDiamond, abs_of_nonneg hs, abs_of_nonneg hcz, max_eq_left (sin_le_cos_of hθ.1 hθ.2)]
  have e2 : ∫ θ in (π / 4)..(π / 2), hDiamond θ * Real.cos (4 * k * θ)
      = ∫ θ in (π / 4)..(π / 2), Real.sin θ * Real.cos (4 * k * θ) := by
    refine intervalIntegral.integral_congr (fun θ hθ => ?_)
    rw [uIcc_of_le (by linarith)] at hθ
    have hs : 0 ≤ Real.sin θ := Real.sin_nonneg_of_nonneg_of_le_pi (by linarith [hθ.1]) (by linarith [hθ.2])
    have hcz : 0 ≤ Real.cos θ := Real.cos_nonneg_of_mem_Icc ⟨by linarith [hθ.1], by linarith [hθ.2]⟩
    simp only [hDiamond, abs_of_nonneg hs, abs_of_nonneg hcz, max_eq_right (cos_le_sin_of hθ.1 hθ.2)]
  have i1 := intervalIntegral.integral_eq_sub_of_hasDerivAt (a := 0) (b := π / 4)
    (fun θ _ => hd_cos k θ) (by apply Continuous.intervalIntegrable; fun_prop)
  have i2 := intervalIntegral.integral_eq_sub_of_hasDerivAt (a := π / 4) (b := π / 2)
    (fun θ _ => hd_sin k θ) (by apply Continuous.intervalIntegrable; fun_prop)
  unfold a
  rw [← hsplit, e1, e2, i1, i2]
  simp only [mul_zero, Real.sin_zero, zero_div, add_zero]
  have q1 : (4 * (k : ℝ) + 1) * (π / 4) = π / 4 + k * π := by ring
  have q2 : (4 * (k : ℝ) - 1) * (π / 4) = k * π - π / 4 := by ring
  have q3 : (4 * (k : ℝ) + 1) * (π / 2) = π / 2 + ((2 * k : ℕ) : ℝ) * π := by push_cast; ring
  have q4 : (4 * (k : ℝ) - 1) * (π / 2) = ((2 * k : ℕ) : ℝ) * π - π / 2 := by push_cast; ring
  rw [q1, q2, q3, q4, Real.sin_add_nat_mul_pi, Real.sin_nat_mul_pi_sub, Real.cos_add_nat_mul_pi,
    Real.cos_nat_mul_pi_sub, Real.cos_add_nat_mul_pi, Real.cos_nat_mul_pi_sub, Real.cos_pi_div_two,
    Real.sin_pi_div_four, Real.cos_pi_div_four]
  have a1 := four_k_ne k; have b1 := four_k_ne' k
  have c1 : (16 * (k : ℝ) ^ 2 - 1) ≠ 0 := by
    rw [show (16 * (k : ℝ) ^ 2 - 1) = (4 * k - 1) * (4 * k + 1) by ring]; exact mul_ne_zero a1 b1
  rw [show (16 * (k : ℝ) ^ 2 - 1) = (4 * k - 1) * (4 * k + 1) by ring]
  field_simp
  ring

end Y04

theorem Y04.proof_1 : Y04.stmt_1 := Y04.A_formula

open Real MeasureTheory Set

namespace Y04

lemma thetaX_zero : thetaX 0 = π / 4 := by
  simp [thetaX, Real.arctan_one]

lemma theta_hasDeriv (x : ℝ) :
    HasDerivAt thetaX (-(Real.sin (thetaX x) * Real.cos (thetaX x))) x := by
  have h : HasDerivAt thetaX _ x :=
    (Real.hasDerivAt_arctan (Real.exp (-x))).comp x ((hasDerivAt_neg x).exp)
  refine h.congr_deriv ?_
  simp only [thetaX, Real.sin_arctan, Real.cos_arctan]
  have hs : Real.sqrt (1 + Real.exp (-x) ^ 2) ^ 2 = 1 + Real.exp (-x) ^ 2 :=
    Real.sq_sqrt (by positivity)
  have hp : 0 < Real.sqrt (1 + Real.exp (-x) ^ 2) := Real.sqrt_pos.2 (by positivity)
  field_simp
  rw [hs]

/-- `F_k(t) = cos(4kt) cos² t sin t`, so that `G_k = F_k ∘ θ`. -/
noncomputable def F (k : ℕ) (t : ℝ) : ℝ := Real.cos (4 * k * t) * Real.cos t ^ 2 * Real.sin t

noncomputable def Fd (k : ℕ) (t : ℝ) : ℝ :=
  -(4 * k) * Real.sin (4 * k * t) * Real.cos t ^ 2 * Real.sin t
    + Real.cos (4 * k * t) * (Real.cos t ^ 3 - 2 * Real.cos t * Real.sin t ^ 2)

lemma G_eq (k : ℕ) : G k = fun x => F k (thetaX x) := rfl

lemma F_hasDeriv (k : ℕ) (t : ℝ) : HasDerivAt (F k) (Fd k t) t := by
  have h := ((((hasDerivAt_id t).const_mul (4 * (k : ℝ))).cos).mul
    ((Real.hasDerivAt_cos t).pow 2)).mul (Real.hasDerivAt_sin t)
  convert h using 1
  · funext s; simp [F, mul_assoc]
  · simp only [Fd, id, Pi.mul_apply, Pi.pow_apply]; push_cast; ring

noncomputable def H (k : ℕ) (t : ℝ) : ℝ := -(Fd k t * (Real.sin t * Real.cos t))

lemma G_hasDeriv (k : ℕ) (x : ℝ) : HasDerivAt (G k) (H k (thetaX x)) x := by
  have h : HasDerivAt (G k) _ x := (F_hasDeriv k (thetaX x)).comp x (theta_hasDeriv x)
  refine h.congr_deriv ?_
  simp only [H]; ring


lemma deriv_G (k : ℕ) : deriv (G k) = fun x => H k (thetaX x) :=
  funext fun x => (G_hasDeriv k x).deriv

lemma sin_k (k : ℕ) : Real.sin (4 * (k : ℝ) * (π / 4)) = 0 := by
  rw [show 4 * (k : ℝ) * (π / 4) = k * π by ring, Real.sin_nat_mul_pi]

lemma cos_k (k : ℕ) : Real.cos (4 * (k : ℝ) * (π / 4)) = (-1) ^ k := by
  rw [show 4 * (k : ℝ) * (π / 4) = k * π by ring, Real.cos_nat_mul_pi]

lemma G_contDiff (k : ℕ) : ContDiff ℝ (⊤ : ℕ∞) (G k) := by
  unfold G thetaX
  have ha : ContDiff ℝ (⊤ : ℕ∞) (fun x : ℝ => Real.arctan (Real.exp (-x))) :=
    Real.contDiff_arctan.comp (Real.contDiff_exp.comp contDiff_neg)
  exact ((Real.contDiff_cos.comp (contDiff_const.mul ha)).mul
    ((Real.contDiff_cos.comp ha).pow 2)).mul (Real.contDiff_sin.comp ha)

lemma G_bdd (k : ℕ) (x : ℝ) : |G k x| ≤ 1 := by
  unfold G
  rw [abs_mul, abs_mul, abs_pow]
  have h1 := Real.abs_cos_le_one (4 * k * thetaX x)
  have h2 := Real.abs_cos_le_one (thetaX x)
  have h3 := Real.abs_sin_le_one (thetaX x)
  have h2' : |Real.cos (thetaX x)| ^ 2 ≤ 1 := by
    have := abs_nonneg (Real.cos (thetaX x)); nlinarith
  calc |Real.cos (4 * ↑k * thetaX x)| * |Real.cos (thetaX x)| ^ 2 * |Real.sin (thetaX x)|
      ≤ 1 * 1 * 1 := by gcongr
    _ = 1 := by ring

lemma G_zero (k : ℕ) : G k 0 = (-1 : ℝ) ^ k * Real.sqrt 2 / 4 := by
  show F k (thetaX 0) = _
  simp only [thetaX_zero, F, cos_k, Real.cos_pi_div_four, Real.sin_pi_div_four]
  have : Real.sqrt 2 ^ 2 = 2 := Real.sq_sqrt (by norm_num)
  field_simp; rw [this]; ring

lemma G_d1 (k : ℕ) : deriv (G k) 0 = (-1 : ℝ) ^ k * Real.sqrt 2 / 8 := by
  rw [deriv_G]; simp only [thetaX_zero, H, Fd, cos_k, sin_k, Real.cos_pi_div_four,
    Real.sin_pi_div_four]
  have : Real.sqrt 2 ^ 2 = 2 := Real.sq_sqrt (by norm_num)
  have h5 : Real.sqrt 2 ^ 5 = 4 * Real.sqrt 2 := by
    rw [show Real.sqrt 2 ^ 5 = (Real.sqrt 2 ^ 2) ^ 2 * Real.sqrt 2 by ring, this]; ring
  ring_nf; simp only [h5]; ring

noncomputable def Fdd (k : ℕ) (t : ℝ) : ℝ :=
  -(4 * k) ^ 2 * Real.cos (4 * k * t) * Real.cos t ^ 2 * Real.sin t
    + 4 * (4 * k) * Real.sin (4 * k * t) * Real.cos t * Real.sin t ^ 2
    - 2 * (4 * k) * Real.sin (4 * k * t) * Real.cos t ^ 3
    + Real.cos (4 * k * t) * (2 * Real.sin t ^ 3 - 7 * Real.cos t ^ 2 * Real.sin t)

lemma Fd_hasDeriv (k : ℕ) (t : ℝ) : HasDerivAt (Fd k) (Fdd k t) t := by
  have h := ((((((hasDerivAt_id t).const_mul (4 * (k : ℝ))).sin).const_mul (-(4 * (k : ℝ)))).mul
      ((Real.hasDerivAt_cos t).pow 2)).mul (Real.hasDerivAt_sin t)).add
      ((((hasDerivAt_id t).const_mul (4 * (k : ℝ))).cos).mul
        (((Real.hasDerivAt_cos t).pow 3).sub
          (((Real.hasDerivAt_cos t).const_mul 2).mul ((Real.hasDerivAt_sin t).pow 2))))
  convert h using 1
  · funext s; simp [Fd, mul_assoc]
  · simp only [Fdd, id, Pi.mul_apply, Pi.pow_apply, Pi.sub_apply]; push_cast; ring

noncomputable def Hd (k : ℕ) (t : ℝ) : ℝ :=
  -(Fdd k t * (Real.sin t * Real.cos t) + Fd k t * (Real.cos t * Real.cos t - Real.sin t * Real.sin t))

lemma H_hasDeriv (k : ℕ) (t : ℝ) : HasDerivAt (H k) (Hd k t) t := by
  have h : HasDerivAt (H k) _ t :=
    ((Fd_hasDeriv k t).mul ((Real.hasDerivAt_sin t).mul (Real.hasDerivAt_cos t))).neg
  refine h.congr_deriv ?_
  simp only [Hd, Pi.mul_apply]; ring

lemma G_d2 (k : ℕ) :
    iteratedDeriv 2 (G k) 0 / 2 = -((-1 : ℝ) ^ k * Real.sqrt 2 * ((k : ℝ) ^ 2 / 2 + 5 / 32)) := by
  rw [iteratedDeriv_succ, iteratedDeriv_one, deriv_G]
  have hc := (H_hasDeriv k (thetaX 0)).comp 0 (theta_hasDeriv 0)
  rw [show (fun x => H k (thetaX x)) = H k ∘ thetaX from rfl, hc.deriv]
  simp only [thetaX_zero, Hd, Fdd, Fd, cos_k, sin_k, Real.cos_pi_div_four, Real.sin_pi_div_four]
  have : Real.sqrt 2 ^ 2 = 2 := Real.sq_sqrt (by norm_num)
  have h7 : Real.sqrt 2 ^ 7 = 8 * Real.sqrt 2 := by
    rw [show Real.sqrt 2 ^ 7 = (Real.sqrt 2 ^ 2) ^ 3 * Real.sqrt 2 by ring, this]; ring
  ring_nf
  rw [h7]; ring

end Y04

theorem Y04.proof_3 : Y04.stmt_3 := fun k =>
  ⟨Y04.G_contDiff k, ⟨1, fun x _ => Y04.G_bdd k x⟩, Y04.G_zero k, Y04.G_d1 k, Y04.G_d2 k⟩

open Real MeasureTheory Set
open Filter Topology Finset

namespace Y04

/-- `ζ(n)` as a real number (local helper). -/
noncomputable def zn (n : ℕ) : ℝ := (riemannZeta n).re

lemma z_hasSum {k : ℕ} (hk : 1 < k) :
    HasSum (fun n : ℕ => 1 / ((n : ℝ) + 1) ^ k) (zn k) := by
  have hs : Summable (fun n : ℕ => 1 / (n : ℝ) ^ k) := Real.summable_one_div_nat_pow.mpr hk
  have hz : zn k = ∑' n : ℕ, 1 / (n : ℝ) ^ k := by
    unfold zn
    rw [zeta_nat_eq_tsum_of_gt_one hk]
    rw [show (fun n : ℕ => 1 / (n : ℂ) ^ k) = fun n : ℕ => ((1 / (n : ℝ) ^ k : ℝ) : ℂ) by
      ext n; push_cast; rfl]
    rw [← Complex.ofReal_tsum, Complex.ofReal_re]
  rw [hz]
  have := hs.hasSum
  rw [← hasSum_nat_add_iff' 1] at this
  simpa [zero_pow (show k ≠ 0 by omega)] using this

lemma prod_shift (x : ℝ) (n : ℕ) :
    ∏ j ∈ range (n + 1), (1 + x + (j : ℝ)) =
      ((n + 1).factorial : ℝ) * ∏ k ∈ range (n + 1), (1 + x / ((k : ℝ) + 1)) := by
  rw [show (((n + 1).factorial : ℕ) : ℝ) = ∏ k ∈ range (n + 1), ((k : ℝ) + 1) by
    rw [← Finset.prod_range_add_one_eq_factorial]; push_cast; rfl, ← Finset.prod_mul_distrib]
  refine Finset.prod_congr rfl fun k _ => ?_
  have : (0 : ℝ) < (k : ℝ) + 1 := by positivity
  field_simp
  ring

lemma fac_pos {x : ℝ} (hx : |x| < 1) (k : ℕ) : 0 < 1 + x / ((k : ℝ) + 1) := by
  have h1 : (1 : ℝ) ≤ (k : ℝ) + 1 := by have := k.cast_nonneg (α := ℝ); linarith
  have : |x / ((k : ℝ) + 1)| < 1 := by
    rw [abs_div, abs_of_pos (by linarith : (0:ℝ) < (k : ℝ) + 1), div_lt_one (by linarith)]
    linarith
  linarith [neg_abs_le (x / ((k : ℝ) + 1))]

/-- the partial products `R_m` -/
noncomputable def Rp (η u v : ℝ) (m : ℕ) : ℝ :=
  ∏ k ∈ range m, ((1 + η / ((k : ℝ) + 1)) / ((1 + u / ((k : ℝ) + 1)) * (1 + v / ((k : ℝ) + 1))))

lemma gammaSeq_ratio {η u v : ℝ} (huv : u + v = η) (hu : |u| < 1) (hv : |v| < 1)
    (hη : |η| < 1) {n : ℕ} (hn : 1 ≤ n) :
    Real.GammaSeq (1 + u) n * Real.GammaSeq (1 + v) n / Real.GammaSeq (1 + η) n =
      (n : ℝ) / ((n : ℝ) + 1) * Rp η u v (n + 1) := by
  unfold Real.GammaSeq Rp
  rw [prod_shift, prod_shift, prod_shift, Finset.prod_div_distrib, Finset.prod_mul_distrib]
  have hn0 : (0 : ℝ) < n := by exact_mod_cast hn
  have pu : 0 < ∏ k ∈ range (n + 1), (1 + u / ((k : ℝ) + 1)) :=
    Finset.prod_pos fun k _ => fac_pos hu k
  have pv : 0 < ∏ k ∈ range (n + 1), (1 + v / ((k : ℝ) + 1)) :=
    Finset.prod_pos fun k _ => fac_pos hv k
  have pη : 0 < ∏ k ∈ range (n + 1), (1 + η / ((k : ℝ) + 1)) :=
    Finset.prod_pos fun k _ => fac_pos hη k
  have e : ∀ x : ℝ, (n : ℝ) ^ (1 + x) = n * (n : ℝ) ^ x := fun x => by
    rw [Real.rpow_add hn0, Real.rpow_one]
  have euv : (n : ℝ) ^ u * (n : ℝ) ^ v = (n : ℝ) ^ η := by rw [← Real.rpow_add hn0, huv]
  rw [e, e, e]
  have hpη : 0 < (n : ℝ) ^ η := Real.rpow_pos_of_pos hn0 _
  have hf : (0 : ℝ) < (n.factorial : ℝ) := by exact_mod_cast n.factorial_pos
  have hf1 : (0 : ℝ) < ((n + 1).factorial : ℝ) := by exact_mod_cast (n + 1).factorial_pos
  rw [Nat.factorial_succ] at hf1 ⊢
  push_cast at hf1 ⊢
  field_simp
  rw [← euv]

lemma log_series_k {η u v : ℝ} (huv : u + v = η) (hu : |u| < 1) (hv : |v| < 1)
    (hη : |η| < 1) (k : ℕ) :
    HasSum (fun j : ℕ => (-1 : ℝ) ^ (j + 2) * (1 / ((k : ℝ) + 1) ^ (j + 2)) / ((j + 2 : ℕ) : ℝ) *
        (u ^ (j + 2) + v ^ (j + 2) - η ^ (j + 2)))
      (Real.log (1 + η / ((k : ℝ) + 1)) - Real.log (1 + u / ((k : ℝ) + 1))
        - Real.log (1 + v / ((k : ℝ) + 1))) := by
  have hk : (0 : ℝ) < (k : ℝ) + 1 := by positivity
  have sm : ∀ x : ℝ, |x| < 1 → |-(x / ((k : ℝ) + 1))| < 1 := by
    intro x hx
    have h1 : (1 : ℝ) ≤ (k : ℝ) + 1 := by have := k.cast_nonneg (α := ℝ); linarith
    rw [abs_neg, abs_div, abs_of_pos hk, div_lt_one hk]; linarith
  have Lu := Real.hasSum_pow_div_log_of_abs_lt_one (sm u hu)
  have Lv := Real.hasSum_pow_div_log_of_abs_lt_one (sm v hv)
  have Lη := Real.hasSum_pow_div_log_of_abs_lt_one (sm η hη)
  have L := (Lu.add Lv).sub Lη
  simp only [sub_neg_eq_add] at L
  rw [← hasSum_nat_add_iff' 1] at L
  simp only [Finset.range_one, Finset.sum_singleton] at L
  convert L using 1
  · ext j
    simp only [neg_pow (_ / _), div_pow]
    push_cast
    ring
  · rw [← huv]
    ring_nf

lemma log_R {η u v : ℝ} (hu : |u| < 1) (hv : |v| < 1) (hη : |η| < 1) (m : ℕ) :
    Real.log (Rp η u v m) = ∑ k ∈ range m, (Real.log (1 + η / ((k : ℝ) + 1))
      - Real.log (1 + u / ((k : ℝ) + 1)) - Real.log (1 + v / ((k : ℝ) + 1))) := by
  unfold Rp
  rw [Real.log_prod]
  · refine Finset.sum_congr rfl fun k _ => ?_
    have a := fac_pos hu k; have b := fac_pos hv k; have c := fac_pos hη k
    rw [Real.log_div c.ne' (by positivity), Real.log_mul a.ne' b.ne']
    ring
  · intro k _
    have a := fac_pos hu k; have b := fac_pos hv k; have c := fac_pos hη k
    positivity

lemma R_tendsto {η u v : ℝ} (huv : u + v = η) (hu : |u| < 1) (hv : |v| < 1)
    (hη : |η| < 1) (hη1 : -1 < η) :
    Tendsto (fun n : ℕ => Rp η u v (n + 1)) atTop
      (𝓝 (Real.Gamma (1 + u) * Real.Gamma (1 + v) / Real.Gamma (1 + η))) := by
  have hG : Real.Gamma (1 + η) ≠ 0 := (Real.Gamma_pos_of_pos (by linarith)).ne'
  have h1 := ((Real.GammaSeq_tendsto_Gamma (1 + u)).mul
    (Real.GammaSeq_tendsto_Gamma (1 + v))).div (Real.GammaSeq_tendsto_Gamma (1 + η)) hG
  have h2 : Tendsto (fun n : ℕ => ((n : ℝ) + 1) / n) atTop (𝓝 1) := by
    have := (tendsto_one_div_atTop_nhds_zero_nat).const_add (1 : ℝ)
    rw [add_zero] at this
    refine this.congr' ?_
    filter_upwards [eventually_ge_atTop 1] with n hn
    have : (0 : ℝ) < n := by exact_mod_cast hn
    field_simp
  have h3 := h1.mul h2
  rw [mul_one] at h3
  refine h3.congr' ?_
  filter_upwards [eventually_ge_atTop 1] with n hn
  have : (0 : ℝ) < n := by exact_mod_cast hn
  simp only [Pi.div_apply]
  rw [gammaSeq_ratio huv hu hv hη hn]
  field_simp


theorem gamma_log_series (η s : ℝ) (hη : |η| < 1 / 2) (hs : |s| < 1 / 2) :
    0 < Real.Gamma (1 + (η / 2 + s)) * Real.Gamma (1 + (η / 2 - s)) / Real.Gamma (1 + η) ∧
    HasSum (fun k : ℕ => (-1 : ℝ) ^ (k + 2) * zn (k + 2) / ((k + 2 : ℕ) : ℝ) *
        ((η / 2 + s) ^ (k + 2) + (η / 2 - s) ^ (k + 2) - η ^ (k + 2)))
      (Real.log (Real.Gamma (1 + (η / 2 + s)) * Real.Gamma (1 + (η / 2 - s)) / Real.Gamma (1 + η))) := by
  set u := η / 2 + s with hudef
  set v := η / 2 - s with hvdef
  have huv : u + v = η := by rw [hudef, hvdef]; ring
  obtain ⟨hη1, hη2⟩ := abs_lt.mp hη
  obtain ⟨hs1, hs2⟩ := abs_lt.mp hs
  have hu : |u| < 3 / 4 := by rw [abs_lt]; constructor <;> linarith
  have hv : |v| < 3 / 4 := by rw [abs_lt]; constructor <;> linarith
  have hu' : |u| < 1 := by linarith
  have hv' : |v| < 1 := by linarith
  have hη' : |η| < 1 := by linarith
  have hApos : 0 < Real.Gamma (1 + u) * Real.Gamma (1 + v) / Real.Gamma (1 + η) := by
    have := Real.Gamma_pos_of_pos (show 0 < 1 + u by cases abs_lt.mp hu; linarith)
    have := Real.Gamma_pos_of_pos (show 0 < 1 + v by cases abs_lt.mp hv; linarith)
    have := Real.Gamma_pos_of_pos (show 0 < 1 + η by linarith)
    positivity
  refine ⟨hApos, ?_⟩
  -- the double series
  set f : ℕ × ℕ → ℝ := fun p => (-1 : ℝ) ^ (p.2 + 2) * (1 / ((p.1 : ℝ) + 1) ^ (p.2 + 2)) /
      ((p.2 + 2 : ℕ) : ℝ) * (u ^ (p.2 + 2) + v ^ (p.2 + 2) - η ^ (p.2 + 2)) with hfdef
  have hfs : Summable f := by
    have hg1 : Summable (fun k : ℕ => 1 / ((k : ℝ) + 1) ^ 2) := (z_hasSum (by norm_num)).summable
    have hg2 : Summable (fun j : ℕ => 3 * (3 / 4 : ℝ) ^ j) :=
      (summable_geometric_of_lt_one (by norm_num) (by norm_num)).mul_left 3
    refine Summable.of_norm_bounded (hg1.mul_of_nonneg hg2 (fun k => by positivity)
      (fun j => by positivity)) ?_
    rintro ⟨k, j⟩
    simp only [hfdef, Real.norm_eq_abs]
    have hk1 : (1 : ℝ) ≤ (k : ℝ) + 1 := by have := k.cast_nonneg (α := ℝ); linarith
    have hj : (1 : ℝ) ≤ ((j + 2 : ℕ) : ℝ) := by push_cast; have := j.cast_nonneg (α := ℝ); linarith
    rw [abs_mul, abs_div, abs_mul, abs_pow, abs_neg, abs_one, one_pow, one_mul,
      abs_of_pos (by positivity : (0 : ℝ) < 1 / ((k : ℝ) + 1) ^ (j + 2)),
      abs_of_pos (by positivity : (0 : ℝ) < ((j + 2 : ℕ) : ℝ))]
    have b1 : 1 / ((k : ℝ) + 1) ^ (j + 2) ≤ 1 / ((k : ℝ) + 1) ^ 2 := by
      apply one_div_le_one_div_of_le (by positivity)
      exact pow_le_pow_right₀ hk1 (by omega)
    have b2 : |u ^ (j + 2) + v ^ (j + 2) - η ^ (j + 2)| ≤ 3 * (3 / 4 : ℝ) ^ j := by
      have pu : |u| ^ (j + 2) ≤ (3 / 4 : ℝ) ^ j :=
        (pow_le_pow_left₀ (abs_nonneg u) hu.le _).trans
          (pow_le_pow_of_le_one (by norm_num) (by norm_num) (by omega))
      have pv : |v| ^ (j + 2) ≤ (3 / 4 : ℝ) ^ j :=
        (pow_le_pow_left₀ (abs_nonneg v) hv.le _).trans
          (pow_le_pow_of_le_one (by norm_num) (by norm_num) (by omega))
      have pη : |η| ^ (j + 2) ≤ (3 / 4 : ℝ) ^ j :=
        (pow_le_pow_left₀ (abs_nonneg η) (show |η| ≤ 3 / 4 by linarith) _).trans
          (pow_le_pow_of_le_one (by norm_num) (by norm_num) (by omega))
      calc |u ^ (j + 2) + v ^ (j + 2) - η ^ (j + 2)|
          ≤ |u ^ (j + 2)| + |v ^ (j + 2)| + |η ^ (j + 2)| := by
            have := abs_sub (u ^ (j + 2) + v ^ (j + 2)) (η ^ (j + 2))
            have := abs_add_le (u ^ (j + 2)) (v ^ (j + 2))
            linarith
        _ ≤ 3 * (3 / 4 : ℝ) ^ j := by rw [abs_pow, abs_pow, abs_pow]; linarith
    calc 1 / ((k : ℝ) + 1) ^ (j + 2) / ((j + 2 : ℕ) : ℝ) * |u ^ (j + 2) + v ^ (j + 2) - η ^ (j + 2)|
        ≤ 1 / ((k : ℝ) + 1) ^ 2 / 1 * (3 * (3 / 4 : ℝ) ^ j) := by
          gcongr
      _ = _ := by simp
  have hf := hfs.hasSum
  -- fibers in k
  have hT := hf.prod_fiberwise (fun k => log_series_k huv hu' hv' hη' k)
  -- fibers in j
  have hf' : HasSum (f ∘ Equiv.prodComm ℕ ℕ) (∑' p, f p) := (Equiv.prodComm ℕ ℕ).hasSum_iff.mpr hf
  have hZ := hf'.prod_fiberwise (fun j => by
    have := ((z_hasSum (show 1 < j + 2 by omega)).mul_left ((-1 : ℝ) ^ (j + 2))).div_const
      ((j + 2 : ℕ) : ℝ) |>.mul_right (u ^ (j + 2) + v ^ (j + 2) - η ^ (j + 2))
    exact this)
  -- identify the sum
  have hlim1 := (tendsto_add_atTop_iff_nat 1).mpr hT.tendsto_sum_nat
  have hlim2 : Tendsto (fun n : ℕ => ∑ k ∈ range (n + 1), (Real.log (1 + η / ((k : ℝ) + 1))
      - Real.log (1 + u / ((k : ℝ) + 1)) - Real.log (1 + v / ((k : ℝ) + 1)))) atTop
      (𝓝 (Real.log (Real.Gamma (1 + u) * Real.Gamma (1 + v) / Real.Gamma (1 + η)))) := by
    have := ((Real.continuousAt_log hApos.ne').tendsto).comp
      (R_tendsto huv hu' hv' hη' (by linarith))
    refine this.congr fun n => ?_
    simp only [Function.comp]
    exact log_R hu' hv' hη' (n + 1)
  have heq := tendsto_nhds_unique hlim1 hlim2
  rw [heq] at hZ
  convert hZ using 1


lemma abs_mono_bd {a b x B : ℝ} (ha : |a| ≤ B) (hb : |b| ≤ B) (hx : |x| ≤ 1) (m : ℕ) :
    |a * b * x ^ m| ≤ B ^ 2 := by
  rw [abs_mul, abs_mul, abs_pow]
  have h1 : |x| ^ m ≤ 1 := pow_le_one₀ (abs_nonneg _) hx
  have h2 : |a| * |b| ≤ B * B := mul_le_mul ha hb (abs_nonneg _) ((abs_nonneg _).trans ha)
  calc |a| * |b| * |x| ^ m ≤ (B * B) * 1 :=
        mul_le_mul h2 h1 (by positivity) (by nlinarith [abs_nonneg a])
    _ = B ^ 2 := by ring

lemma abs_coef_bd {a x B : ℝ} (ha : |a| ≤ B) (hx : |x| ≤ 1) (m : ℕ) :
    |a * x ^ (m + 2)| ≤ B * x ^ 2 := by
  rw [abs_mul, abs_pow]
  have : |x| ^ (m + 2) ≤ x ^ 2 := by
    rw [pow_add, ← sq_abs x, mul_comm]
    exact mul_le_of_le_one_right (by positivity) (pow_le_one₀ (abs_nonneg _) hx)
  exact mul_le_mul ha this (by positivity) ((abs_nonneg _).trans ha)

lemma exp_trunc {L a0 a1 a2 a3 x B K : ℝ} (hB : 0 ≤ B) (hK : 0 ≤ K)
    (ha0 : |a0| ≤ B) (ha1 : |a1| ≤ B) (ha2 : |a2| ≤ B) (ha3 : |a3| ≤ B)
    (hx : |x| ≤ 1) (hsmall : (4 * B + K) * x ^ 2 ≤ 1)
    (hL : |L - (a0 * x ^ 2 + a1 * x ^ 3 + a2 * x ^ 4 + a3 * x ^ 5)| ≤ K * x ^ 6) :
    |Real.exp L - (1 + a0 * x ^ 2 + a1 * x ^ 3 + (a2 + a0 ^ 2 / 2) * x ^ 4
        + (a3 + a0 * a1) * x ^ 5)| ≤ (K + K * (8 * B + K) + 7 * B ^ 2 + (4 * B + K) ^ 3) * x ^ 6 := by
  obtain ⟨P, hPdef⟩ : ∃ P, P = a0 * x ^ 2 + a1 * x ^ 3 + a2 * x ^ 4 + a3 * x ^ 5 := ⟨_, rfl⟩
  rw [← hPdef] at hL
  have hx2 : 0 ≤ x ^ 2 := sq_nonneg x
  have hx21 : x ^ 2 ≤ 1 := by rw [← sq_abs]; exact pow_le_one₀ (abs_nonneg _) hx
  have hx6 : x ^ 6 ≤ x ^ 2 := by
    rw [show x ^ 6 = x ^ 2 * (x ^ 2) ^ 2 by ring]
    exact mul_le_of_le_one_right hx2 (pow_le_one₀ hx2 hx21)
  have hP : |P| ≤ 4 * B * x ^ 2 := by
    have e0 := abs_coef_bd ha0 hx 0
    have e1 := abs_coef_bd ha1 hx 1
    have e2 := abs_coef_bd ha2 hx 2
    have e3 := abs_coef_bd ha3 hx 3
    simp only [zero_add] at e0 e1 e2 e3
    have t1 := abs_add_le (a0 * x ^ 2 + a1 * x ^ 3 + a2 * x ^ 4) (a3 * x ^ 5)
    have t2 := abs_add_le (a0 * x ^ 2 + a1 * x ^ 3) (a2 * x ^ 4)
    have t3 := abs_add_le (a0 * x ^ 2) (a1 * x ^ 3)
    rw [hPdef]
    linarith
  have hLb : |L| ≤ (4 * B + K) * x ^ 2 := by
    have : |L| ≤ |P| + |L - P| := by
      have := abs_add_le P (L - P); rwa [add_sub_cancel] at this
    nlinarith
  have hL1 : |L| ≤ 1 := hLb.trans hsmall
  have he := Real.exp_bound hL1 (n := 3) (by norm_num)
  simp only [Finset.sum_range_succ, Finset.sum_range_zero, Nat.factorial] at he
  norm_num at he
  have hL3 : |L| ^ 3 ≤ (4 * B + K) ^ 3 * x ^ 6 := by
    calc |L| ^ 3 ≤ ((4 * B + K) * x ^ 2) ^ 3 := pow_le_pow_left₀ (abs_nonneg _) hLb 3
      _ = (4 * B + K) ^ 3 * x ^ 6 := by ring
  have hE1 : |Real.exp L - (1 + L + L ^ 2 / 2)| ≤ (4 * B + K) ^ 3 * x ^ 6 := by
    have : |L| ^ 3 * (2 / 9) ≤ |L| ^ 3 := by nlinarith [pow_nonneg (abs_nonneg L) 3]
    linarith
  have hE2 : |(L ^ 2 - P ^ 2) / 2| ≤ K * (8 * B + K) * x ^ 6 := by
    have hLP : |L + P| ≤ (8 * B + K) * x ^ 2 := by
      calc |L + P| ≤ |L| + |P| := abs_add_le _ _
        _ ≤ (8 * B + K) * x ^ 2 := by linarith
    rw [show (L ^ 2 - P ^ 2) / 2 = (L - P) * (L + P) / 2 by ring, abs_div, abs_mul,
      abs_of_pos (by norm_num : (0:ℝ) < 2)]
    have := mul_le_mul hL hLP (abs_nonneg _) (by positivity)
    have h8 : 0 ≤ K * (8 * B + K) * x ^ 6 := by positivity
    have : K * x ^ 6 * ((8 * B + K) * x ^ 2) ≤ K * (8 * B + K) * x ^ 6 := by
      have := mul_le_mul_of_nonneg_left hx21 h8
      nlinarith
    linarith
  obtain ⟨W, hW⟩ : ∃ W : ℝ, W = a1 ^ 2 / 2 + a0 * a2 + (a0 * a3 + a1 * a2) * x + (a2 ^ 2 / 2 + a1 * a3) * x ^ 2
    + a2 * a3 * x ^ 3 + a3 ^ 2 / 2 * x ^ 4 := ⟨_, rfl⟩
  have hWb : |W| ≤ 7 * B ^ 2 := by
    have b1 := abs_mono_bd ha1 ha1 hx 0
    have b2 := abs_mono_bd ha0 ha2 hx 0
    have b3 := abs_mono_bd ha0 ha3 hx 1
    have b4 := abs_mono_bd ha1 ha2 hx 1
    have b5 := abs_mono_bd ha2 ha2 hx 2
    have b6 := abs_mono_bd ha1 ha3 hx 2
    have b7 := abs_mono_bd ha2 ha3 hx 3
    have b8 := abs_mono_bd ha3 ha3 hx 4
    have eW : W = (a1 * a1 * x ^ 0) / 2 + a0 * a2 * x ^ 0 + a0 * a3 * x ^ 1 + a1 * a2 * x ^ 1
        + (a2 * a2 * x ^ 2) / 2 + a1 * a3 * x ^ 2 + a2 * a3 * x ^ 3 + (a3 * a3 * x ^ 4) / 2 := by
      rw [hW]; ring
    rw [eW]
    rw [abs_le] at b1 b2 b3 b4 b5 b6 b7 b8 ⊢
    constructor <;> linarith [b1.1, b1.2, b2.1, b2.2, b3.1, b3.2, b4.1, b4.2, b5.1, b5.2,
      b6.1, b6.2, b7.1, b7.2, b8.1, b8.2]
  have hE3 : |x ^ 6 * W| ≤ 7 * B ^ 2 * x ^ 6 := by
    rw [abs_mul, abs_of_nonneg (by positivity : (0:ℝ) ≤ x ^ 6), mul_comm]
    exact mul_le_mul_of_nonneg_right hWb (by positivity)
  have key : Real.exp L - (1 + a0 * x ^ 2 + a1 * x ^ 3 + (a2 + a0 ^ 2 / 2) * x ^ 4
        + (a3 + a0 * a1) * x ^ 5) =
      (Real.exp L - (1 + L + L ^ 2 / 2)) + (L - P) + (L ^ 2 - P ^ 2) / 2 + x ^ 6 * W := by
    rw [hW, hPdef]; ring
  rw [key]
  have t1 := abs_add_le ((Real.exp L - (1 + L + L ^ 2 / 2)) + (L - P) + (L ^ 2 - P ^ 2) / 2)
    (x ^ 6 * W)
  have t2 := abs_add_le ((Real.exp L - (1 + L + L ^ 2 / 2)) + (L - P)) ((L ^ 2 - P ^ 2) / 2)
  have t3 := abs_add_le (Real.exp L - (1 + L + L ^ 2 / 2)) (L - P)
  have hk : 0 ≤ K * (8 * B + K) * x ^ 6 := by positivity
  linarith


lemma z_two : zn 2 = π ^ 2 / 6 := by
  unfold zn
  rw [show ((2 : ℕ) : ℂ) = 2 by norm_num, riemannZeta_two]
  rw [show (π : ℂ) ^ 2 / 6 = ((π ^ 2 / 6 : ℝ) : ℂ) by push_cast; ring, Complex.ofReal_re]

lemma z_four : zn 4 = π ^ 4 / 90 := by
  unfold zn
  rw [show ((4 : ℕ) : ℂ) = 4 by norm_num, riemannZeta_four]
  rw [show (π : ℂ) ^ 4 / 90 = ((π ^ 4 / 90 : ℝ) : ℂ) by push_cast; ring, Complex.ofReal_re]

lemma z_bd {n : ℕ} (hn : 2 ≤ n) : |zn n| ≤ 2 := by
  have h := z_hasSum (k := n) (by omega)
  have h2 := z_hasSum (k := 2) (by omega)
  have hnn : 0 ≤ zn n := h.nonneg fun m => by positivity
  have hle : zn n ≤ zn 2 := hasSum_le (fun m => by
    have h1 : (1 : ℝ) ≤ (m : ℝ) + 1 := by have := m.cast_nonneg (α := ℝ); linarith
    exact one_div_le_one_div_of_le (by positivity) (pow_le_pow_right₀ h1 hn)) h h2
  rw [z_two] at hle
  have := Real.pi_lt_d2
  have : π ^ 2 / 6 ≤ 2 := by nlinarith [Real.pi_pos]
  rw [abs_of_nonneg hnn]; linarith

/-- coefficient `a_k` of the one-variable series -/
noncomputable def cf (α β γ : ℝ) (k : ℕ) : ℝ :=
  (-1 : ℝ) ^ (k + 2) * zn (k + 2) / ((k + 2 : ℕ) : ℝ) * (α ^ (k + 2) + β ^ (k + 2) - γ ^ (k + 2))

lemma cf_bd {α β γ : ℝ} (hα : |α| ≤ 3 / 2) (hβ : |β| ≤ 3 / 2) (hγ : |γ| ≤ 3 / 2) (k : ℕ) :
    |cf α β γ k| ≤ 3 * (3 / 2) ^ (k + 2) := by
  unfold cf
  have hz := z_bd (n := k + 2) (by omega)
  have hn : (2 : ℝ) ≤ ((k + 2 : ℕ) : ℝ) := by push_cast; have := k.cast_nonneg (α := ℝ); linarith
  have hb : |α ^ (k + 2) + β ^ (k + 2) - γ ^ (k + 2)| ≤ 3 * (3 / 2) ^ (k + 2) := by
    have a := pow_le_pow_left₀ (abs_nonneg α) hα (k + 2)
    have b := pow_le_pow_left₀ (abs_nonneg β) hβ (k + 2)
    have c := pow_le_pow_left₀ (abs_nonneg γ) hγ (k + 2)
    rw [← abs_pow] at a b c
    have t1 := abs_sub (α ^ (k + 2) + β ^ (k + 2)) (γ ^ (k + 2))
    have t2 := abs_add_le (α ^ (k + 2)) (β ^ (k + 2))
    linarith
  rw [abs_mul, abs_div, abs_mul, abs_pow, abs_neg, abs_one, one_pow, one_mul,
    abs_of_pos (by linarith : (0:ℝ) < ((k + 2 : ℕ) : ℝ))]
  have : |zn (k + 2)| / ((k + 2 : ℕ) : ℝ) ≤ 1 := by
    rw [div_le_one (by linarith)]; linarith
  calc |zn (k + 2)| / ((k + 2 : ℕ) : ℝ) * |α ^ (k + 2) + β ^ (k + 2) - γ ^ (k + 2)|
      ≤ 1 * (3 * (3 / 2) ^ (k + 2)) := mul_le_mul this hb (abs_nonneg _) (by norm_num)
    _ = _ := one_mul _

lemma series_trunc {α β γ x L : ℝ} (hα : |α| ≤ 3 / 2) (hβ : |β| ≤ 3 / 2) (hγ : |γ| ≤ 3 / 2)
    (hx : |x| ≤ 1 / 3)
    (hL : HasSum (fun k : ℕ => (-1 : ℝ) ^ (k + 2) * zn (k + 2) / ((k + 2 : ℕ) : ℝ) *
      ((α * x) ^ (k + 2) + (β * x) ^ (k + 2) - (γ * x) ^ (k + 2))) L) :
    |L - (cf α β γ 0 * x ^ 2 + cf α β γ 1 * x ^ 3 + cf α β γ 2 * x ^ 4 + cf α β γ 3 * x ^ 5)|
      ≤ 69 * x ^ 6 := by
  have hf : (fun k : ℕ => (-1 : ℝ) ^ (k + 2) * zn (k + 2) / ((k + 2 : ℕ) : ℝ) *
      ((α * x) ^ (k + 2) + (β * x) ^ (k + 2) - (γ * x) ^ (k + 2))) =
      fun k => cf α β γ k * x ^ (k + 2) := by
    ext k; unfold cf; rw [mul_pow, mul_pow, mul_pow]; ring
  rw [hf] at hL
  set r := 3 / 2 * |x| with hr
  have hr0 : 0 ≤ r := by positivity
  have hr1 : r ≤ 1 / 2 := by rw [hr]; linarith
  have hbd : ∀ k, |cf α β γ k * x ^ (k + 2)| ≤ 3 * r ^ (k + 2) := by
    intro k
    rw [abs_mul, abs_pow, hr, mul_pow]
    have := cf_bd hα hβ hγ k
    calc |cf α β γ k| * |x| ^ (k + 2) ≤ 3 * (3 / 2) ^ (k + 2) * |x| ^ (k + 2) :=
          mul_le_mul_of_nonneg_right this (by positivity)
      _ = _ := by ring
  have hg := (hasSum_nat_add_iff' 4).mpr hL
  have hgeo : HasSum (fun k : ℕ => 3 * r ^ 6 * r ^ k) (3 * r ^ 6 * (1 - r)⁻¹) :=
    (hasSum_geometric_of_lt_one hr0 (by linarith)).mul_left _
  have hb1 : ∀ k : ℕ, cf α β γ (k + 4) * x ^ (k + 4 + 2) ≤ 3 * r ^ 6 * r ^ k := by
    intro k
    have := (abs_le.mp (hbd (k + 4))).2
    calc _ ≤ 3 * r ^ (k + 4 + 2) := this
      _ = _ := by ring
  have hb2 : ∀ k : ℕ, -(cf α β γ (k + 4) * x ^ (k + 4 + 2)) ≤ 3 * r ^ 6 * r ^ k := by
    intro k
    have := (abs_le.mp (hbd (k + 4))).1
    calc _ ≤ 3 * r ^ (k + 4 + 2) := by linarith
      _ = _ := by ring
  have u1 := hasSum_le hb1 hg hgeo
  have u2 := hasSum_le hb2 hg.neg hgeo
  have hinv : 3 * r ^ 6 * (1 - r)⁻¹ ≤ 6 * r ^ 6 := by
    rw [mul_inv_le_iff₀ (by linarith)]
    nlinarith [pow_nonneg hr0 6]
  have hr6 : 6 * r ^ 6 ≤ 69 * x ^ 6 := by
    rw [hr, mul_pow, pow_abs, abs_of_nonneg (by positivity : (0:ℝ) ≤ x ^ 6)]
    nlinarith [pow_nonneg (sq_nonneg x) 3]
  simp only [Finset.sum_range_succ, Finset.sum_range_zero, zero_add] at u1 u2
  rw [abs_le]
  constructor <;> norm_num at u1 u2 ⊢ <;> linarith

lemma A_expand {α β γ x L : ℝ} (hα : |α| ≤ 3 / 2) (hβ : |β| ≤ 3 / 2) (hγ : |γ| ≤ 3 / 2)
    (hx : |x| ≤ 1 / 20)
    (hL : HasSum (fun k : ℕ => (-1 : ℝ) ^ (k + 2) * zn (k + 2) / ((k + 2 : ℕ) : ℝ) *
      ((α * x) ^ (k + 2) + (β * x) ^ (k + 2) - (γ * x) ^ (k + 2))) L) :
    |Real.exp L - (1 + cf α β γ 0 * x ^ 2 + cf α β γ 1 * x ^ 3
        + (cf α β γ 2 + cf α β γ 0 ^ 2 / 2) * x ^ 4
        + (cf α β γ 3 + cf α β γ 0 * cf α β γ 1) * x ^ 5)|
      ≤ (69 + 69 * (8 * 23 + 69) + 7 * 23 ^ 2 + (4 * 23 + 69) ^ 3) * x ^ 6 := by
  have h := series_trunc hα hβ hγ (by linarith) hL
  have hB : ∀ k ≤ 3, |cf α β γ k| ≤ 23 := by
    intro k hk
    have := cf_bd hα hβ hγ k
    have : (3 / 2 : ℝ) ^ (k + 2) ≤ (3 / 2) ^ 5 := pow_le_pow_right₀ (by norm_num) (by omega)
    norm_num at this ⊢; linarith
  have hx2 : x ^ 2 ≤ 1 / 400 := by
    have := pow_le_pow_left₀ (abs_nonneg x) hx 2
    rw [sq_abs] at this; norm_num at this ⊢; linarith
  exact exp_trunc (by norm_num) (by norm_num) (hB 0 (by norm_num)) (hB 1 (by norm_num))
    (hB 2 (by norm_num)) (hB 3 (by norm_num)) (by linarith) (by nlinarith) h



lemma eq_zero_of_le {a K d : ℝ} (hd : 0 < d) (hK : 0 ≤ K)
    (h : ∀ s : ℝ, 0 < s → s ≤ d → |a| ≤ K * s ^ 2) : a = 0 := by
  by_contra ha
  have ha' : 0 < |a| := abs_pos.mpr ha
  set s := min d (min 1 (|a| / (K + 1))) with hs
  have hs0 : 0 < s := lt_min hd (lt_min one_pos (by positivity))
  have hs1 : s ≤ 1 := (min_le_right _ _).trans (min_le_left _ _)
  have hs2 : s ≤ |a| / (K + 1) := (min_le_right _ _).trans (min_le_right _ _)
  have := h s hs0 (min_le_left _ _)
  have hsq : s ^ 2 ≤ s := by nlinarith
  have : K * s ^ 2 ≤ K * (|a| / (K + 1)) := mul_le_mul_of_nonneg_left (hsq.trans hs2) hK
  have : K * (|a| / (K + 1)) < |a| := by
    rw [mul_div_assoc', div_lt_iff₀ (by linarith)]; nlinarith
  linarith

end Y04

namespace Y04

lemma zn_three : zn 3 = z3 := by simp [zn, z3]

noncomputable def CR : ℝ := 64 * (69 + 69 * (8 * 23 + 69) + 7 * 23 ^ 2 + (4 * 23 + 69) ^ 3) + 32 * (23 + 23 ^ 2)

lemma CR_pos : 0 < CR := by unfold CR; norm_num

lemma R_expand {x : ℝ} (hx : |x| ≤ 1 / 40) :
    |R x - (1 - (π ^ 2 / 6) * x ^ 2 + 2 * z3 * x ^ 3
      + ((π ^ 2 / 6) ^ 2 / 2 - 7 * (π ^ 4 / 90) / 2) * x ^ 4)| ≤ CR * |x| ^ 5 := by
  have h2x : |2 * x| ≤ 1 / 20 := by rw [abs_mul]; norm_num; linarith
  obtain ⟨hpos, h⟩ := gamma_log_series (2 * x) 0 (by linarith) (by norm_num)
  rw [show 2 * x / 2 + 0 = 1 / 2 * (2 * x) by ring, show 2 * x / 2 - 0 = 1 / 2 * (2 * x) by ring] at h hpos
  have hA := A_expand (α := 1 / 2) (β := 1 / 2) (γ := 1) (by norm_num) (by norm_num) (by norm_num) h2x
    (by simpa only [one_mul] using h)
  rw [Real.exp_log hpos] at hA
  have hRx : Real.Gamma (1 + 1 / 2 * (2 * x)) * Real.Gamma (1 + 1 / 2 * (2 * x)) / Real.Gamma (1 + 2 * x)
      = R x := by
    unfold R; rw [show 1 / 2 * (2 * x) = x by ring, sq]
  rw [hRx] at hA
  have hB : ∀ k ≤ 3, |cf (1 / 2) (1 / 2) 1 k| ≤ 23 := by
    intro k hk
    have := cf_bd (α := 1 / 2) (β := 1 / 2) (γ := 1) (by norm_num) (by norm_num) (by norm_num) k
    have : (3 / 2 : ℝ) ^ (k + 2) ≤ (3 / 2) ^ 5 := pow_le_pow_right₀ (by norm_num) (by omega)
    norm_num at this ⊢; linarith
  set c5 := cf (1 / 2) (1 / 2) 1 3 + cf (1 / 2) (1 / 2) 1 0 * cf (1 / 2) (1 / 2) 1 1 with hc5
  have hc5b : |c5| ≤ 23 + 23 ^ 2 := by
    have t := abs_add_le (cf (1 / 2) (1 / 2) 1 3) (cf (1 / 2) (1 / 2) 1 0 * cf (1 / 2) (1 / 2) 1 1)
    rw [abs_mul] at t
    have := hB 3 le_rfl; have h0 := hB 0 (by norm_num); have h1 := hB 1 (by norm_num)
    have : |cf (1 / 2) (1 / 2) 1 0| * |cf (1 / 2) (1 / 2) 1 1| ≤ 23 * 23 :=
      mul_le_mul h0 h1 (abs_nonneg _) (by norm_num)
    rw [hc5]; nlinarith
  have hpoly : (1 + cf (1 / 2) (1 / 2) 1 0 * (2 * x) ^ 2 + cf (1 / 2) (1 / 2) 1 1 * (2 * x) ^ 3
        + (cf (1 / 2) (1 / 2) 1 2 + cf (1 / 2) (1 / 2) 1 0 ^ 2 / 2) * (2 * x) ^ 4 + c5 * (2 * x) ^ 5)
      = (1 - (π ^ 2 / 6) * x ^ 2 + 2 * z3 * x ^ 3
      + ((π ^ 2 / 6) ^ 2 / 2 - 7 * (π ^ 4 / 90) / 2) * x ^ 4) + c5 * (2 * x) ^ 5 := by
    norm_num [cf, z_two, z_four, zn_three]
    ring
  rw [hpoly] at hA
  have hx1 : |x| ≤ 1 := by linarith
  have e6 : (2 * x) ^ 6 ≤ 64 * |x| ^ 5 := by
    have : x ^ 6 = |x| ^ 5 * |x| := by
      rw [← pow_succ, show x ^ 6 = (x ^ 2) ^ 3 by ring, ← sq_abs, ← pow_mul]
    rw [mul_pow, this]
    nlinarith [pow_nonneg (abs_nonneg x) 5, abs_nonneg x]
  have e5 : |c5 * (2 * x) ^ 5| ≤ (23 + 23 ^ 2) * 32 * |x| ^ 5 := by
    rw [abs_mul, abs_pow, abs_mul, mul_pow]; norm_num
    have := mul_le_mul_of_nonneg_right hc5b (by positivity : (0:ℝ) ≤ 32 * |x| ^ 5)
    nlinarith
  have key := abs_sub_le (R x) ((1 - (π ^ 2 / 6) * x ^ 2 + 2 * z3 * x ^ 3
      + ((π ^ 2 / 6) ^ 2 / 2 - 7 * (π ^ 4 / 90) / 2) * x ^ 4) + c5 * (2 * x) ^ 5)
      (1 - (π ^ 2 / 6) * x ^ 2 + 2 * z3 * x ^ 3
      + ((π ^ 2 / 6) ^ 2 / 2 - 7 * (π ^ 4 / 90) / 2) * x ^ 4)
  rw [add_sub_cancel_left] at key
  unfold CR
  nlinarith [pow_nonneg (abs_nonneg x) 5]

end Y04

namespace Y04

/-- the real Beta integral equals the Gamma quotient (as in X05) -/
lemma beta_real {a b : ℝ} (ha : 0 < a) (hb : 0 < b) :
    (∫ t in (0 : ℝ)..1, t ^ (a - 1) * (1 - t) ^ (b - 1)) * Real.Gamma (a + b) =
      Real.Gamma a * Real.Gamma b := by
  have key := Complex.Gamma_mul_Gamma_eq_betaIntegral (s := (a : ℂ)) (t := (b : ℂ))
    (by simpa using ha) (by simpa using hb)
  have hB : Complex.betaIntegral a b =
      ((∫ t in (0 : ℝ)..1, t ^ (a - 1) * (1 - t) ^ (b - 1) : ℝ) : ℂ) := by
    rw [Complex.betaIntegral, ← intervalIntegral.integral_ofReal]
    refine intervalIntegral.integral_congr_uIoo fun t ht => ?_
    rw [uIoo_of_le zero_le_one] at ht
    push_cast
    rw [Complex.ofReal_cpow ht.1.le, Complex.ofReal_cpow (by linarith [ht.2] : (0:ℝ) ≤ 1 - t)]
    push_cast
    ring_nf
  rw [hB, ← Complex.ofReal_add, Complex.Gamma_ofReal, Complex.Gamma_ofReal,
    Complex.Gamma_ofReal] at key
  have k2 := key.symm
  rw [mul_comm] at k2
  exact_mod_cast k2

/-- the logistic map `s ↦ 1/(1+e^{-s})` -/
noncomputable def logi (s : ℝ) : ℝ := 1 / (1 + Real.exp (-s))

noncomputable def logid (s : ℝ) : ℝ := Real.exp (-s) / (1 + Real.exp (-s)) ^ 2

lemma logi_hasDeriv (s : ℝ) : HasDerivAt logi (logid s) s := by
  have h1 : HasDerivAt (fun s => 1 + Real.exp (-s)) (Real.exp (-s) * -1) s :=
    ((hasDerivAt_neg s).exp).const_add 1
  have h := h1.inv (by positivity)
  have e : logi = fun s => (1 + Real.exp (-s))⁻¹ := by funext s; simp [logi]
  rw [e]; convert h using 1; simp [logid]

lemma logi_strictMono : StrictMono logi := by
  intro s t hst
  unfold logi
  have : Real.exp (-t) < Real.exp (-s) := Real.exp_lt_exp.2 (by linarith)
  apply one_div_lt_one_div_of_lt (by positivity); linarith

lemma logi_image : logi '' Set.Ioi 0 = Set.Ioo (1 / 2) 1 := by
  ext t; constructor
  · rintro ⟨s, hs, rfl⟩
    have h1 : Real.exp (-s) < 1 := Real.exp_lt_one_iff.2 (by simp at hs; linarith)
    have h2 := Real.exp_pos (-s)
    unfold logi
    constructor
    · rw [lt_div_iff₀ (by positivity)]; linarith
    · rw [div_lt_one (by positivity)]; linarith
  · rintro ⟨h1, h2⟩
    have ht0 : 0 < t := by linarith
    have h1t : 0 < 1 - t := by linarith
    refine ⟨Real.log (t / (1 - t)), ?_, ?_⟩
    · apply Real.log_pos; rw [lt_div_iff₀ h1t]; linarith
    · unfold logi
      rw [Real.exp_neg, Real.exp_log (by positivity), inv_div]
      field_simp; ring

lemma logid_eq (s : ℝ) : logid s = logi s * (1 - logi s) := by
  unfold logid logi
  have : 0 < 1 + Real.exp (-s) := by positivity
  field_simp; ring

lemma logid_pos (s : ℝ) : 0 < logid s := by unfold logid; positivity

/-- `log logid(s) = -2 m(s)`, `m(s) = s/2 + ℓ(s)` -/
lemma log_logid (s : ℝ) : Real.log (logid s) = -2 * (s / 2 + ell s) := by
  unfold logid ell
  rw [Real.log_div (by positivity) (by positivity), Real.log_exp, Real.log_pow]
  push_cast; ring

lemma ell_nonneg (s : ℝ) : 0 ≤ ell s := Real.log_nonneg (by linarith [Real.exp_pos (-s)])

lemma ell_le (s : ℝ) : ell s ≤ Real.exp (-s) := by
  have := Real.log_le_sub_one_of_pos (show 0 < 1 + Real.exp (-s) by positivity)
  unfold ell; linarith

lemma ell_le_one (s : ℝ) (hs : 0 ≤ s) : ell s ≤ 1 :=
  (ell_le s).trans (Real.exp_le_one_iff.2 (by linarith))

lemma continuous_ell : Continuous ell := by
  unfold ell
  exact (continuous_const.add (Real.continuous_exp.comp continuous_neg)).log
    (fun x => (show (0:ℝ) < 1 + Real.exp (-x) by positivity).ne')

/-- the Beta half-line integral: `∫_0^∞ e^{-2x m(s)} ds = Γ(x)²/(2Γ(2x))` for `x > 0` -/
lemma cosh_integral {x : ℝ} (hx : 0 < x) :
    IntegrableOn (fun s => Real.exp (-2 * x * (s / 2 + ell s))) (Set.Ioi 0) ∧
    ∫ s in Set.Ioi (0 : ℝ), Real.exp (-2 * x * (s / 2 + ell s)) =
      Real.Gamma x * Real.Gamma x / Real.Gamma (2 * x) / 2 := by
  set g : ℝ → ℝ := fun t => t ^ (x - 1) * (1 - t) ^ (x - 1) with hg
  have hder : ∀ s ∈ Set.Ioi (0 : ℝ), HasDerivWithinAt logi (logid s) (Set.Ioi 0) s :=
    fun s _ => (logi_hasDeriv s).hasDerivWithinAt
  have hinj : InjOn logi (Set.Ioi 0) := logi_strictMono.injective.injOn
  have hpt : ∀ s, |logid s| • g (logi s) = Real.exp (-2 * x * (s / 2 + ell s)) := by
    intro s
    have hl := logid_pos s
    have h0 : 0 < logi s := by unfold logi; positivity
    have h1 : 0 < 1 - logi s := by
      unfold logi; rw [sub_pos, div_lt_one (by positivity)]; linarith [Real.exp_pos (-s)]
    simp only [hg, smul_eq_mul, abs_of_pos hl]
    rw [← Real.mul_rpow h0.le h1.le, ← logid_eq]
    rw [show logid s * logid s ^ (x - 1) = logid s ^ (1 + (x - 1)) by
      rw [Real.rpow_add hl, Real.rpow_one]]
    rw [show 1 + (x - 1) = x by ring, Real.rpow_def_of_pos hl, log_logid]; ring_nf
  -- integrability of the left side
  have hint : IntegrableOn (fun s => Real.exp (-2 * x * (s / 2 + ell s))) (Set.Ioi 0) := by
    have hb : IntegrableOn (fun s : ℝ => Real.exp (-x * s)) (Set.Ioi 0) :=
      integrableOn_exp_mul_Ioi (by linarith) 0
    refine Integrable.mono' hb ((Real.continuous_exp.comp (continuous_const.mul
      ((continuous_id.div_const 2).add continuous_ell))).aestronglyMeasurable) ?_
    refine (ae_restrict_iff' measurableSet_Ioi).2 (Filter.Eventually.of_forall (fun s _ => ?_))
    rw [Real.norm_eq_abs, abs_of_pos (Real.exp_pos _)]
    apply Real.exp_le_exp.2
    have := ell_nonneg s
    nlinarith
  refine ⟨hint, ?_⟩
  have himg := integral_image_eq_integral_abs_deriv_smul measurableSet_Ioi hder hinj g
  have himgI := (integrableOn_image_iff_integrableOn_abs_deriv_smul measurableSet_Ioi hder hinj g).2
    (by simp_rw [hpt]; exact hint)
  simp_rw [hpt] at himg
  rw [logi_image] at himg himgI
  rw [← himg]
  -- the set integral over (1/2,1) is half the Beta integral
  have hII : IntervalIntegrable g volume (1 / 2) 1 := by
    rw [intervalIntegrable_iff_integrableOn_Ioc_of_le (by norm_num)]
    exact (integrableOn_Ioc_iff_integrableOn_Ioo (by simp [hg])).2 himgI
  have hsym : ∀ t, g (1 - t) = g t := fun t => by simp only [hg]; ring_nf
  have hII2 : IntervalIntegrable g volume 0 (1 / 2) := by
    have := hII.comp_sub_left 1
    simp only [hsym] at this
    norm_num at this
    exact this.symm
  have hsplit := intervalIntegral.integral_add_adjacent_intervals hII2 hII
  have hrefl : ∫ t in (0:ℝ)..(1 / 2), g t = ∫ t in (1 / 2 : ℝ)..1, g t := by
    have := intervalIntegral.integral_comp_sub_left (a := (1 / 2 : ℝ)) (b := 1) g 1
    simp only [hsym] at this
    norm_num at this
    exact this.symm
  have hbeta := beta_real hx hx
  have hG : 0 < Real.Gamma (x + x) := Real.Gamma_pos_of_pos (by linarith)
  have hset : ∫ t in Set.Ioo (1 / 2 : ℝ) 1, g t = ∫ t in (1 / 2 : ℝ)..1, g t := by
    rw [intervalIntegral.integral_of_le (by norm_num), integral_Ioc_eq_integral_Ioo]
  rw [hset]
  rw [← hsplit, hrefl] at hbeta
  rw [show 2 * x = x + x by ring, ← hbeta]
  field_simp
  ring

end Y04

namespace Y04

lemma intOn_pow_exp (j : ℕ) {b : ℝ} (hb : 0 < b) :
    IntegrableOn (fun s : ℝ => s ^ j * Real.exp (-(b * s))) (Set.Ioi 0) := by
  have := integrableOn_rpow_mul_exp_neg_mul_rpow (s := (j : ℝ)) (p := 1) (b := b)
    (by have : (0 : ℝ) ≤ j := Nat.cast_nonneg j; linarith) one_pos hb
  refine this.congr_fun (fun x _ => ?_) measurableSet_Ioi
  simp [Real.rpow_natCast, neg_mul]

lemma intOn_pow_exp1 (j : ℕ) : IntegrableOn (fun s : ℝ => s ^ j * Real.exp (-s)) (Set.Ioi 0) := by
  simpa using intOn_pow_exp j one_pos

/-- the dominating function `(s/2+1)^3 e^{-s}` -/
noncomputable def W (s : ℝ) : ℝ := (s / 2 + 1) ^ 3 * Real.exp (-s)

lemma W_int : IntegrableOn W (Set.Ioi 0) := by
  have h := ((((intOn_pow_exp1 3).const_mul (1 / 8)).add ((intOn_pow_exp1 2).const_mul (3 / 4))).add
    ((intOn_pow_exp1 1).const_mul (3 / 2))).add (intOn_pow_exp1 0)
  refine IntegrableOn.congr_fun h (fun s _ => ?_) measurableSet_Ioi
  simp only [W, Pi.add_apply]; ring_nf

lemma W_nonneg (s : ℝ) (hs : 0 ≤ s) : 0 ≤ W s := by unfold W; positivity

lemma continuous_W : Continuous W := by unfold W; fun_prop

lemma intOn_of_bound {f g : ℝ → ℝ} (hf : Continuous f) (hg : IntegrableOn g (Set.Ioi 0))
    (h : ∀ s : ℝ, 0 < s → |f s| ≤ g s) : IntegrableOn f (Set.Ioi 0) := by
  refine Integrable.mono' hg hf.aestronglyMeasurable ?_
  exact (ae_restrict_iff' measurableSet_Ioi).2 (Filter.Eventually.of_forall
    (fun s hs => by rw [Real.norm_eq_abs]; exact h s hs))

lemma Jint_bound {j : ℕ} (hj1 : 1 ≤ j) (hj3 : j ≤ 3) {s : ℝ} (hs : 0 ≤ s) :
    |(s / 2 + ell s) ^ j - (s / 2) ^ j| ≤ 3 * W s := by
  have hl0 := ell_nonneg s
  have hl1 := ell_le_one s hs
  have hle := ell_le s
  set l := ell s
  set a := s / 2 with ha
  have ha0 : 0 ≤ a := by rw [ha]; linarith
  have key : |(a + l) ^ j - a ^ j| ≤ 3 * (a + 1) ^ 2 * l := by
    interval_cases j
    · rw [show (a + l) ^ 1 - a ^ 1 = l by ring, abs_of_nonneg hl0]
      have : 1 ≤ 3 * (a + 1) ^ 2 := by nlinarith
      nlinarith
    · rw [show (a + l) ^ 2 - a ^ 2 = l * (2 * a + l) by ring, abs_of_nonneg (by positivity)]
      nlinarith [mul_nonneg hl0 ha0, mul_nonneg hl0 (mul_nonneg ha0 ha0)]
    · rw [show (a + l) ^ 3 - a ^ 3 = l * (3 * a ^ 2 + 3 * a * l + l ^ 2) by ring,
        abs_of_nonneg (by positivity)]
      have : 3 * a ^ 2 + 3 * a * l + l ^ 2 ≤ 3 * (a + 1) ^ 2 := by nlinarith
      nlinarith
  calc _ ≤ 3 * (a + 1) ^ 2 * l := key
    _ ≤ 3 * (a + 1) ^ 3 * Real.exp (-s) := by
      have : (a + 1) ^ 2 ≤ (a + 1) ^ 3 := by nlinarith [sq_nonneg (a + 1)]
      have h2 : 0 ≤ (a + 1) ^ 2 := by positivity
      calc 3 * (a + 1) ^ 2 * l ≤ 3 * (a + 1) ^ 2 * Real.exp (-s) := by gcongr
        _ ≤ _ := by gcongr
    _ = 3 * W s := by simp only [W, ha]; ring

lemma Jint_int {j : ℕ} (hj1 : 1 ≤ j) (hj3 : j ≤ 3) :
    IntegrableOn (fun s : ℝ => (s / 2 + ell s) ^ j - (s / 2) ^ j) (Set.Ioi 0) :=
  intOn_of_bound (by have := continuous_ell; fun_prop) (W_int.const_mul 3)
    (fun s hs => Jint_bound hj1 hj3 hs.le)

/-- `|e^c - 1 - c - c²/2| ≤ 4|c|³` for `c ≤ 0` -/
lemma r1_bound {c : ℝ} (hc : c ≤ 0) : |Real.exp c - 1 - c - c ^ 2 / 2| ≤ 4 * |c| ^ 3 := by
  rcases le_or_gt (|c|) 1 with h | h
  · have he := Real.exp_bound h (n := 3) (by norm_num)
    simp only [Finset.sum_range_succ, Finset.sum_range_zero, Nat.factorial] at he
    norm_num at he
    have : Real.exp c - 1 - c - c ^ 2 / 2 = Real.exp c - (1 + c + c ^ 2 / 2) := by ring
    rw [this]
    nlinarith [pow_nonneg (abs_nonneg c) 3]
  · have hc' : |c| = -c := abs_of_nonpos hc
    have he : Real.exp c ≤ 1 := Real.exp_le_one_iff.2 hc
    have he0 := Real.exp_pos c
    rw [abs_le]
    have h3 : 1 < |c| ^ 3 := one_lt_pow₀ h (by norm_num)
    have hc2 : c ^ 2 = |c| ^ 2 := (sq_abs c).symm
    have h2 : |c| ^ 2 ≤ |c| ^ 3 := by nlinarith
    constructor <;> nlinarith

lemma r_diff {a b M : ℝ} (ha : a ≤ 0) (hb : b ≤ 0) (haM : -M ≤ a) (hbM : -M ≤ b) :
    |(Real.exp a - 1 - a - a ^ 2 / 2 - a ^ 3 / 6) - (Real.exp b - 1 - b - b ^ 2 / 2 - b ^ 3 / 6)|
      ≤ 4 * M ^ 3 * |a - b| := by
  have hM : 0 ≤ M := by linarith
  have hd : ∀ u ∈ Set.Icc (-M) 0, HasDerivWithinAt
      (fun u => Real.exp u - 1 - u - u ^ 2 / 2 - u ^ 3 / 6)
      (Real.exp u - 1 - u - u ^ 2 / 2) (Set.Icc (-M) 0) u := by
    intro u _
    have h : HasDerivAt (fun u => Real.exp u - 1 - u - u ^ 2 / 2 - u ^ 3 / 6) _ u :=
      ((((Real.hasDerivAt_exp u).sub_const 1).sub (hasDerivAt_id u)).sub
      ((hasDerivAt_pow 2 u).div_const 2)).sub ((hasDerivAt_pow 3 u).div_const 6)
    refine (h.congr_deriv ?_).hasDerivWithinAt
    norm_num; ring
  have hbd : ∀ u ∈ Set.Icc (-M) 0, ‖Real.exp u - 1 - u - u ^ 2 / 2‖ ≤ 4 * M ^ 3 := by
    intro u hu
    rw [Real.norm_eq_abs]
    refine (r1_bound hu.2).trans ?_
    have : |u| ≤ M := by rw [abs_le]; constructor <;> linarith [hu.1, hu.2]
    gcongr
  have := Convex.norm_image_sub_le_of_norm_hasDerivWithin_le hd hbd (convex_Icc (-M) 0)
    ⟨hbM, hb⟩ ⟨haM, ha⟩
  simpa [Real.norm_eq_abs] using this

/-- the pointwise third-order expansion of the `Ψ` integrand -/
lemma pt_expand {τ s : ℝ} (hτ : τ ≤ 0) (_hτ1 : -1 ≤ τ) (hs : 0 ≤ s) :
    |Real.exp (τ * s / 2) * (Real.exp (τ * ell s) - 1)
      - (τ * ((s / 2 + ell s) ^ 1 - (s / 2) ^ 1) + τ ^ 2 / 2 * ((s / 2 + ell s) ^ 2 - (s / 2) ^ 2)
        + τ ^ 3 / 6 * ((s / 2 + ell s) ^ 3 - (s / 2) ^ 3))| ≤ 4 * τ ^ 4 * W s := by
  have hl0 := ell_nonneg s
  have hl1 := ell_le_one s hs
  have hle := ell_le s
  set l := ell s
  set m := s / 2 + l with hm
  have hm0 : 0 ≤ m := by rw [hm]; linarith
  have e1 : Real.exp (τ * s / 2) * (Real.exp (τ * l) - 1) = Real.exp (τ * m) - Real.exp (τ * (s / 2)) := by
    rw [mul_sub, ← Real.exp_add, hm]; ring_nf
  rw [e1]
  have hr := r_diff (a := τ * m) (b := τ * (s / 2)) (M := -τ * m) (by nlinarith) (by nlinarith)
    (by nlinarith) (by nlinarith)
  have eq : Real.exp (τ * m) - Real.exp (τ * (s / 2))
      - (τ * (m ^ 1 - (s / 2) ^ 1) + τ ^ 2 / 2 * (m ^ 2 - (s / 2) ^ 2) + τ ^ 3 / 6 * (m ^ 3 - (s / 2) ^ 3))
      = (Real.exp (τ * m) - 1 - τ * m - (τ * m) ^ 2 / 2 - (τ * m) ^ 3 / 6)
        - (Real.exp (τ * (s / 2)) - 1 - τ * (s / 2) - (τ * (s / 2)) ^ 2 / 2 - (τ * (s / 2)) ^ 3 / 6) := by
    ring
  rw [eq]
  refine hr.trans ?_
  have hab : |τ * m - τ * (s / 2)| = -τ * l := by
    rw [show τ * m - τ * (s / 2) = τ * l by rw [hm]; ring, abs_mul, abs_of_nonpos hτ,
      abs_of_nonneg hl0]
  rw [hab]
  have hmb : m ≤ s / 2 + 1 := by rw [hm]; linarith
  have : (-τ * m) ^ 3 * (-τ * l) = τ ^ 4 * (m ^ 3 * l) := by ring
  have h4 : 0 ≤ τ ^ 4 := by positivity
  have hm3 : m ^ 3 * l ≤ (s / 2 + 1) ^ 3 * Real.exp (-s) := by
    have : m ^ 3 ≤ (s / 2 + 1) ^ 3 := pow_le_pow_left₀ hm0 hmb 3
    exact mul_le_mul this hle hl0 (by positivity)
  calc 4 * (-τ * m) ^ 3 * (-τ * l) = 4 * τ ^ 4 * (m ^ 3 * l) := by ring
    _ ≤ 4 * τ ^ 4 * ((s / 2 + 1) ^ 3 * Real.exp (-s)) := by gcongr
    _ = 4 * τ ^ 4 * W s := rfl

lemma Psi_int {τ : ℝ} (hτ : τ ≤ 0) :
    IntegrableOn (fun s : ℝ => Real.exp (τ * s / 2) * (Real.exp (τ * ell s) - 1)) (Set.Ioi 0) := by
  refine intOn_of_bound (by have := continuous_ell; fun_prop)
    ((intOn_pow_exp1 0).const_mul (-τ)) (fun s hs => ?_)
  have hl0 := ell_nonneg s
  have hle := ell_le s
  have h1 : Real.exp (τ * s / 2) ≤ 1 := Real.exp_le_one_iff.2 (by nlinarith)
  have h2 : |Real.exp (τ * ell s) - 1| ≤ -τ * ell s := by
    have : τ * ell s ≤ 0 := by nlinarith
    rw [abs_of_nonpos (by linarith [Real.exp_le_one_iff.2 this])]
    linarith [Real.add_one_le_exp (τ * ell s)]
  rw [abs_mul, abs_of_pos (Real.exp_pos _)]
  simp only [pow_zero, one_mul]
  calc Real.exp (τ * s / 2) * |Real.exp (τ * ell s) - 1| ≤ 1 * (-τ * ell s) := by
        gcongr
    _ ≤ -τ * Real.exp (-s) := by nlinarith

lemma Psi_expand {τ : ℝ} (hτ : τ ≤ 0) (hτ1 : -1 ≤ τ) :
    |(∫ s in Set.Ioi (0 : ℝ), Real.exp (τ * s / 2) * (Real.exp (τ * ell s) - 1))
      - (τ * J 1 + τ ^ 2 / 2 * J 2 + τ ^ 3 / 6 * J 3)|
      ≤ 4 * τ ^ 4 * ∫ s in Set.Ioi (0 : ℝ), W s := by
  have i1 := (Jint_int (j := 1) le_rfl (by norm_num)).const_mul τ
  have i2 := (Jint_int (j := 2) (by norm_num) (by norm_num)).const_mul (τ ^ 2 / 2)
  have i3 := (Jint_int (j := 3) (by norm_num) le_rfl).const_mul (τ ^ 3 / 6)
  have hg12 : Integrable (fun s : ℝ => τ * ((s / 2 + ell s) ^ 1 - (s / 2) ^ 1)
      + τ ^ 2 / 2 * ((s / 2 + ell s) ^ 2 - (s / 2) ^ 2)) (volume.restrict (Set.Ioi 0)) := i1.add i2
  have hg : Integrable (fun s : ℝ => τ * ((s / 2 + ell s) ^ 1 - (s / 2) ^ 1)
      + τ ^ 2 / 2 * ((s / 2 + ell s) ^ 2 - (s / 2) ^ 2)
      + τ ^ 3 / 6 * ((s / 2 + ell s) ^ 3 - (s / 2) ^ 3)) (volume.restrict (Set.Ioi 0)) := hg12.add i3
  have hsplit : (∫ s in Set.Ioi (0 : ℝ), Real.exp (τ * s / 2) * (Real.exp (τ * ell s) - 1))
      - (τ * J 1 + τ ^ 2 / 2 * J 2 + τ ^ 3 / 6 * J 3)
      = ∫ s in Set.Ioi (0 : ℝ), (Real.exp (τ * s / 2) * (Real.exp (τ * ell s) - 1)
      - (τ * ((s / 2 + ell s) ^ 1 - (s / 2) ^ 1) + τ ^ 2 / 2 * ((s / 2 + ell s) ^ 2 - (s / 2) ^ 2)
        + τ ^ 3 / 6 * ((s / 2 + ell s) ^ 3 - (s / 2) ^ 3))) := by
    symm
    rw [integral_sub (Psi_int hτ) hg, integral_add hg12 i3, integral_add i1 i2,
      integral_const_mul, integral_const_mul, integral_const_mul]
    rfl
  rw [hsplit, ← integral_const_mul]
  have := norm_integral_le_of_norm_le (W_int.const_mul (4 * τ ^ 4))
    ((ae_restrict_iff' measurableSet_Ioi).2 (Filter.Eventually.of_forall
      (fun s (hs : 0 < s) => by rw [Real.norm_eq_abs]; exact pt_expand hτ hτ1 hs.le)))
  simpa [Real.norm_eq_abs] using this

/-- the closed form of `Ψ` for real `τ < 0` -/
lemma Psi_neg {τ : ℝ} (hτ : τ < 0) :
    ∫ s in Set.Ioi (0 : ℝ), Real.exp (τ * s / 2) * (Real.exp (τ * ell s) - 1) =
      (R (-τ / 2) - 1) / (-τ / 2) := by
  set x := -τ / 2 with hx
  have hx0 : 0 < x := by rw [hx]; linarith
  obtain ⟨hI, hval⟩ := cosh_integral hx0
  have e : ∀ s, Real.exp (τ * s / 2) * (Real.exp (τ * ell s) - 1)
      = Real.exp (-2 * x * (s / 2 + ell s)) - Real.exp (-x * s) := by
    intro s; rw [mul_sub, ← Real.exp_add, hx]; ring_nf
  simp_rw [e]
  have hb : IntegrableOn (fun s : ℝ => Real.exp (-x * s)) (Set.Ioi 0) :=
    integrableOn_exp_mul_Ioi (by linarith) 0
  rw [integral_sub hI hb, hval, integral_exp_mul_Ioi (by linarith) 0]
  unfold R
  rw [add_comm 1 x, add_comm 1 (2 * x), Real.Gamma_add_one hx0.ne',
    Real.Gamma_add_one (by positivity)]
  have hG : 0 < Real.Gamma (2 * x) := Real.Gamma_pos_of_pos (by linarith)
  simp only [mul_zero, Real.exp_zero]
  field_simp

lemma eq_zero_of_le_lin {a K d : ℝ} (hd : 0 < d) (hK : 0 ≤ K)
    (h : ∀ t : ℝ, 0 < t → t ≤ d → |a| ≤ K * t) : a = 0 := by
  by_contra ha
  have ha' : 0 < |a| := abs_pos.mpr ha
  set t := min d (|a| / (K + 1)) with ht
  have ht0 : 0 < t := lt_min hd (by positivity)
  have ht2 : t ≤ |a| / (K + 1) := min_le_right _ _
  have := h t ht0 (min_le_left _ _)
  have : K * t ≤ K * (|a| / (K + 1)) := mul_le_mul_of_nonneg_left ht2 hK
  have : K * (|a| / (K + 1)) < |a| := by
    rw [mul_div_assoc', div_lt_iff₀ (by linarith)]; nlinarith
  linarith

/-- three coefficients vanish if the cubic is `O(t⁴)` on `(0, d]` -/
lemma cubic_zero {a1 a2 a3 K d : ℝ} (hd : 0 < d) (hd1 : d ≤ 1) (hK : 0 ≤ K)
    (h : ∀ t : ℝ, 0 < t → t ≤ d → |a1 * t + a2 * t ^ 2 + a3 * t ^ 3| ≤ K * t ^ 4) :
    a1 = 0 ∧ a2 = 0 ∧ a3 = 0 := by
  have h1 : a1 = 0 := by
    refine eq_zero_of_le_lin hd (K := |a2| + |a3| + K)
      (by have := abs_nonneg a2; have := abs_nonneg a3; linarith) (fun t ht0 htd => ?_)
    have ht1 : t ≤ 1 := htd.trans hd1
    have := h t ht0 htd
    have e : a1 = (a1 * t + a2 * t ^ 2 + a3 * t ^ 3) / t - a2 * t - a3 * t ^ 2 := by
      field_simp; ring
    rw [e]
    have b1 : |(a1 * t + a2 * t ^ 2 + a3 * t ^ 3) / t| ≤ K * t ^ 3 := by
      rw [abs_div, abs_of_pos ht0, div_le_iff₀ ht0]; nlinarith
    have b2 : |a2 * t| ≤ |a2| * t := by rw [abs_mul, abs_of_pos ht0]
    have b3 : |a3 * t ^ 2| ≤ |a3| * t := by
      rw [abs_mul, abs_of_pos (pow_pos ht0 2)]
      exact mul_le_mul_of_nonneg_left (by nlinarith) (abs_nonneg _)
    have t3 : K * t ^ 3 ≤ K * t := mul_le_mul_of_nonneg_left (by
      have : t ^ 3 = t * t ^ 2 := by ring
      rw [this]; nlinarith [sq_nonneg t]) hK
    have := abs_sub (((a1 * t + a2 * t ^ 2 + a3 * t ^ 3) / t) - a2 * t) (a3 * t ^ 2)
    have := abs_sub ((a1 * t + a2 * t ^ 2 + a3 * t ^ 3) / t) (a2 * t)
    have : (|a2| + |a3| + K) * t = |a2| * t + |a3| * t + K * t := by ring
    linarith
  subst h1
  have h2 : a2 = 0 := by
    refine eq_zero_of_le_lin hd (K := |a3| + K)
      (by have := abs_nonneg a3; linarith) (fun t ht0 htd => ?_)
    have ht1 : t ≤ 1 := htd.trans hd1
    have := h t ht0 htd
    have e : a2 = (0 * t + a2 * t ^ 2 + a3 * t ^ 3) / t ^ 2 - a3 * t := by
      field_simp; ring
    rw [e]
    have b1 : |(0 * t + a2 * t ^ 2 + a3 * t ^ 3) / t ^ 2| ≤ K * t ^ 2 := by
      rw [abs_div, abs_of_pos (pow_pos ht0 2), div_le_iff₀ (pow_pos ht0 2)]; nlinarith
    have b3 : |a3 * t| ≤ |a3| * t := by rw [abs_mul, abs_of_pos ht0]
    have t3 : K * t ^ 2 ≤ K * t := mul_le_mul_of_nonneg_left (by nlinarith) hK
    have := abs_sub ((0 * t + a2 * t ^ 2 + a3 * t ^ 3) / t ^ 2) (a3 * t)
    have : (|a3| + K) * t = |a3| * t + K * t := by ring
    linarith
  subst h2
  refine ⟨rfl, rfl, eq_zero_of_le_lin hd (K := K) hK (fun t ht0 htd => ?_)⟩
  have := h t ht0 htd
  have e : a3 = (0 * t + 0 * t ^ 2 + a3 * t ^ 3) / t ^ 3 := by field_simp; ring
  rw [e, abs_div, abs_of_pos (pow_pos ht0 3), div_le_iff₀ (pow_pos ht0 3)]
  have ht1 : t ≤ 1 := htd.trans hd1
  have : K * t ^ 4 ≤ K * t * t ^ 3 := by
    rw [show K * t * t ^ 3 = K * t ^ 4 by ring]
  linarith

lemma J_values : J 1 = (π ^ 2 / 6) / 2 ∧ J 2 = z3 ∧ J 3 = 27 / 16 * (π ^ 4 / 90) := by
  set KW := ∫ s in Set.Ioi (0 : ℝ), W s
  have hKW : 0 ≤ KW := setIntegral_nonneg measurableSet_Ioi (fun s hs => W_nonneg s (le_of_lt hs))
  set c4 := (π ^ 2 / 6) ^ 2 / 2 - 7 * (π ^ 4 / 90) / 2 with hc4
  have hmain : ∀ t : ℝ, 0 < t → t ≤ 1 / 20 →
      |(-J 1 + (π ^ 2 / 6) / 2) * t + (J 2 / 2 - z3 / 2) * t ^ 2 + (-J 3 / 6 - c4 / 8) * t ^ 3|
        ≤ (4 * KW + CR / 16) * t ^ 4 := by
    intro t ht0 ht1
    have hP := Psi_expand (τ := -t) (by linarith) (by linarith)
    rw [Psi_neg (by linarith)] at hP
    have hx : |(- -t / 2)| ≤ 1 / 40 := by
      rw [neg_neg, abs_of_pos (by linarith)]; linarith
    have hR := R_expand hx
    set x := - -t / 2 with hxdef
    have hx0 : 0 < x := by rw [hxdef]; linarith
    have hxt : x = t / 2 := by rw [hxdef]; ring
    -- divide the R-expansion by x
    have hR' : |(R x - 1) / x - (-(π ^ 2 / 6) * x + 2 * z3 * x ^ 2 + c4 * x ^ 3)| ≤ CR * x ^ 4 := by
      have e : (R x - 1) / x - (-(π ^ 2 / 6) * x + 2 * z3 * x ^ 2 + c4 * x ^ 3)
          = (R x - (1 - (π ^ 2 / 6) * x ^ 2 + 2 * z3 * x ^ 3 + c4 * x ^ 4)) / x := by
        field_simp; ring
      rw [e, abs_div, abs_of_pos hx0, div_le_iff₀ hx0]
      rw [abs_of_pos hx0] at hR
      calc _ ≤ CR * x ^ 5 := hR
        _ = CR * x ^ 4 * x := by ring
    have e2 : (-(J 1) + (π ^ 2 / 6) / 2) * t + (J 2 / 2 - z3 / 2) * t ^ 2 + (-J 3 / 6 - c4 / 8) * t ^ 3
        = ((R x - 1) / x - (-(π ^ 2 / 6) * x + 2 * z3 * x ^ 2 + c4 * x ^ 3))
          - ((R x - 1) / x - (-t * J 1 + (-t) ^ 2 / 2 * J 2 + (-t) ^ 3 / 6 * J 3)) := by
      rw [hxt]; ring
    rw [e2]
    have t1 := abs_sub ((R x - 1) / x - (-(π ^ 2 / 6) * x + 2 * z3 * x ^ 2 + c4 * x ^ 3))
      ((R x - 1) / x - (-t * J 1 + (-t) ^ 2 / 2 * J 2 + (-t) ^ 3 / 6 * J 3))
    have e3 : CR * x ^ 4 = CR / 16 * t ^ 4 := by rw [hxt]; ring
    have e4 : 4 * (-t) ^ 4 * KW = 4 * KW * t ^ 4 := by ring
    rw [e3] at hR'; rw [e4] at hP
    linarith
  obtain ⟨h1, h2, h3⟩ := cubic_zero (by norm_num) (by norm_num)
    (by have := CR_pos; positivity) hmain
  refine ⟨by linarith, by linarith, ?_⟩
  have : J 3 = -(6 / 8) * c4 := by linarith
  rw [this, hc4]; ring

end Y04

namespace Y04

lemma summable_inv_pow (k : ℕ) (hk : 2 ≤ k) :
    Summable (fun n : ℕ => 1 / ((n : ℝ) + 1) ^ k) := by
  have := (Real.summable_one_div_nat_pow (p := k)).2 (by omega)
  exact (summable_nat_add_iff 1).2 this |>.congr (fun n => by push_cast; ring_nf)

lemma hasSum_alt (k : ℕ) {S : ℝ}
    (hS : HasSum (fun n : ℕ => 1 / ((n : ℝ) + 1) ^ k) S) :
    HasSum (fun n : ℕ => (-1) ^ n * (1 / ((n : ℝ) + 1) ^ k)) (S - 2 * (S / 2 ^ k)) := by
  let a : ℕ → ℝ := fun n => 1 / ((n : ℝ) + 1) ^ k
  have hodd : HasSum (fun m => a (2 * m + 1)) (S / 2 ^ k) := by
    have : (fun m => a (2 * m + 1)) = fun m => a m / 2 ^ k := by
      funext m; simp only [a]; push_cast
      rw [show (2 * (m : ℝ) + 1 + 1) = 2 * ((m : ℝ) + 1) by ring, mul_pow]
      ring
    rw [this]; exact hS.div_const _
  have heven_s : Summable (fun m => a (2 * m)) :=
    hS.summable.comp_injective (fun x y h => by simpa using h)
  have heven : HasSum (fun m => a (2 * m)) (S - S / 2 ^ k) := by
    have h := heven_s.hasSum.even_add_odd hodd
    have := h.unique hS
    convert heven_s.hasSum using 1; linarith
  have he' : HasSum (fun m => (fun n : ℕ => (-1) ^ n * a n) (2 * m)) (S - S / 2 ^ k) :=
    heven.congr_fun (fun m => by
      show (-1 : ℝ) ^ (2 * m) * a (2 * m) = a (2 * m)
      rw [pow_mul, neg_one_sq, one_pow, one_mul])
  have ho' : HasSum (fun m => (fun n : ℕ => (-1) ^ n * a n) (2 * m + 1)) (-(S / 2 ^ k)) :=
    hodd.neg.congr_fun (fun m => by
      show (-1 : ℝ) ^ (2 * m + 1) * a (2 * m + 1) = -a (2 * m + 1)
      rw [pow_succ, pow_mul, neg_one_sq, one_pow, one_mul, neg_one_mul])
  have := HasSum.even_add_odd (f := fun n : ℕ => (-1 : ℝ) ^ n * a n) he' ho'
  rw [show S - 2 * (S / 2 ^ k) = S - S / 2 ^ k + -(S / 2 ^ k) by ring]
  exact this

lemma integral_pow_exp (j : ℕ) {b : ℝ} (hb : 0 < b) :
    ∫ s in Set.Ioi (0 : ℝ), s ^ j * Real.exp (-(b * s)) = (j.factorial : ℝ) / b ^ (j + 1) := by
  have := Real.integral_rpow_mul_exp_neg_mul_Ioi (a := (j : ℝ) + 1) (r := b) (by positivity) hb
  rw [Real.Gamma_nat_eq_factorial] at this
  have e : ∀ s : ℝ, s ^ ((j : ℝ) + 1 - 1) = s ^ j := fun s => by
    rw [add_sub_cancel_right, Real.rpow_natCast]
  simp_rw [e] at this
  rw [this, show ((j : ℝ) + 1) = ((j + 1 : ℕ) : ℝ) by push_cast; ring, Real.rpow_natCast]
  rw [one_div, inv_pow]; field_simp

lemma integral_series (c : ℕ → ℝ) (hc : ∀ n, |c n| ≤ 1) (j : ℕ) :
    ∫ s in Set.Ioi (0 : ℝ), s ^ j * ∑' n : ℕ, c n * Real.exp (-(((n : ℝ) + 1) * s)) / ((n : ℝ) + 1)
      = (j.factorial : ℝ) * ∑' n : ℕ, c n / ((n : ℝ) + 1) ^ (j + 2) := by
  set F : ℕ → ℝ → ℝ := fun n s =>
    c n / ((n : ℝ) + 1) * (s ^ j * Real.exp (-(((n : ℝ) + 1) * s))) with hF
  have hpos : ∀ n : ℕ, (0 : ℝ) < (n : ℝ) + 1 := fun n => by positivity
  have hint : ∀ n, Integrable (F n) (volume.restrict (Set.Ioi 0)) := fun n =>
    (intOn_pow_exp j (hpos n)).const_mul _
  have hval : ∀ n, ∫ s in Set.Ioi (0 : ℝ), F n s
      = (j.factorial : ℝ) * (c n / ((n : ℝ) + 1) ^ (j + 2)) := fun n => by
    simp only [hF]
    rw [integral_const_mul, integral_pow_exp j (hpos n)]
    have := hpos n
    field_simp; ring
  have hnorm : ∀ n, ∫ s in Set.Ioi (0 : ℝ), ‖F n s‖
      = (j.factorial : ℝ) * (|c n| / ((n : ℝ) + 1) ^ (j + 2)) := fun n => by
    have : ∫ s in Set.Ioi (0 : ℝ), ‖F n s‖ = ∫ s in Set.Ioi (0 : ℝ), (|c n| / ((n : ℝ) + 1) *
        (s ^ j * Real.exp (-(((n : ℝ) + 1) * s)))) := by
      refine setIntegral_congr_fun measurableSet_Ioi (fun s hs => ?_)
      simp only [hF, Real.norm_eq_abs, abs_mul, abs_div]
      rw [abs_of_pos (hpos n), abs_of_nonneg (pow_nonneg (le_of_lt (mem_Ioi.mp hs)) j : (0:ℝ) ≤ s ^ j),
        abs_of_pos (Real.exp_pos _)]
    rw [this, integral_const_mul, integral_pow_exp j (hpos n)]
    have := hpos n
    field_simp; ring
  have hsum : Summable fun n => ∫ s in Set.Ioi (0 : ℝ), ‖F n s‖ := by
    simp_rw [hnorm]
    refine Summable.mul_left _ ?_
    refine Summable.of_nonneg_of_le (fun n => by positivity) (fun n => ?_)
      (summable_inv_pow (j + 2) (by omega))
    rw [div_le_div_iff_of_pos_right (by positivity)]; exact hc n
  have key := integral_tsum_of_summable_integral_norm hint hsum
  have e1 : ∀ s, s ^ j * ∑' n : ℕ, c n * Real.exp (-(((n : ℝ) + 1) * s)) / ((n : ℝ) + 1)
      = ∑' n, F n s := fun s => by
    rw [← tsum_mul_left]; congr 1; funext n; simp only [hF]; ring
  simp_rw [e1, ← key, hval]
  rw [tsum_mul_left]

lemma ell_series {s : ℝ} (hs : 0 < s) :
    ell s = ∑' n : ℕ, (-1) ^ n * Real.exp (-(((n : ℝ) + 1) * s)) / ((n : ℝ) + 1) := by
  have hx : |-Real.exp (-s)| < 1 := by
    rw [abs_neg, abs_of_pos (Real.exp_pos _)]; exact Real.exp_lt_one_iff.mpr (by linarith)
  have h := Real.hasSum_pow_div_log_of_abs_lt_one hx
  rw [sub_neg_eq_add] at h
  have e : ∀ n : ℕ, (-Real.exp (-s)) ^ (n + 1) / ((n : ℝ) + 1)
      = -((-1) ^ n * Real.exp (-(((n : ℝ) + 1) * s)) / ((n : ℝ) + 1)) := fun n => by
    rw [neg_pow, ← Real.exp_nat_mul]; push_cast; ring_nf
  simp_rw [e] at h
  have := h.neg.tsum_eq
  simp only [neg_neg] at this
  rw [ell, this]
lemma I_int (m j : ℕ) (hj : 1 ≤ j) : IntegrableOn (fun s : ℝ => s ^ m * ell s ^ j) (Set.Ioi 0) := by
  refine intOn_of_bound (by have := continuous_ell; fun_prop) (intOn_pow_exp1 m) (fun s hs => ?_)
  have hl0 := ell_nonneg s
  have hl1 := ell_le_one s hs.le
  have hle := ell_le s
  rw [abs_of_nonneg (by positivity)]
  have : ell s ^ j ≤ ell s := by
    calc ell s ^ j ≤ ell s ^ 1 := pow_le_pow_of_le_one hl0 hl1 hj
      _ = ell s := pow_one _
  have hm : 0 ≤ s ^ m := pow_nonneg hs.le m
  calc s ^ m * ell s ^ j ≤ s ^ m * ell s := mul_le_mul_of_nonneg_left this hm
    _ ≤ s ^ m * Real.exp (-s) := mul_le_mul_of_nonneg_left hle hm

lemma I_m1 (m : ℕ) : I m 1 = (m.factorial : ℝ) * ∑' n : ℕ, (-1) ^ n / ((n : ℝ) + 1) ^ (m + 2) := by
  have h := integral_series (fun n => (-1) ^ n) (fun n => by simp) m
  have e : I m 1 = ∫ s in Set.Ioi (0 : ℝ), s ^ m * ∑' n : ℕ,
      (-1) ^ n * Real.exp (-(((n : ℝ) + 1) * s)) / ((n : ℝ) + 1) :=
    setIntegral_congr_fun measurableSet_Ioi (fun s hs => by rw [pow_one, ell_series hs])
  rw [e, h]

lemma eta_val (k : ℕ) (hk : 1 < k) :
    ∑' n : ℕ, (-1 : ℝ) ^ n / ((n : ℝ) + 1) ^ k = zn k - 2 * (zn k / 2 ^ k) := by
  have := hasSum_alt k (z_hasSum hk)
  rw [← this.tsum_eq]; congr 1; funext n; ring

lemma I01 : I 0 1 = (π ^ 2 / 6) / 2 := by
  rw [← J_values.1]
  refine setIntegral_congr_fun measurableSet_Ioi (fun s _ => ?_)
  simp

lemma I11 : I 1 1 = 3 / 4 * z3 := by
  rw [I_m1, eta_val 3 (by norm_num), zn_three]; norm_num; ring

lemma I21 : I 2 1 = 7 / 4 * (π ^ 4 / 90) := by
  rw [I_m1, eta_val 4 (by norm_num), z_four]; norm_num; ring

lemma J2_split : J 2 = I 1 1 + I 0 2 := by
  unfold J I
  rw [← integral_add (I_int 1 1 le_rfl) (I_int 0 2 (by norm_num))]
  refine setIntegral_congr_fun measurableSet_Ioi (fun s _ => ?_)
  ring

lemma J3_split : J 3 = 3 / 4 * I 2 1 + 3 / 2 * I 1 2 + I 0 3 := by
  unfold J I
  have h12 : Integrable (fun s : ℝ => 3 / 4 * (s ^ 2 * ell s ^ 1) + 3 / 2 * (s ^ 1 * ell s ^ 2))
      (volume.restrict (Set.Ioi 0)) :=
    ((I_int 2 1 le_rfl).const_mul _).add ((I_int 1 2 (by norm_num)).const_mul _)
  rw [← integral_const_mul, ← integral_const_mul,
    ← integral_add ((I_int 2 1 le_rfl).const_mul _) ((I_int 1 2 (by norm_num)).const_mul _),
    ← integral_add h12 (I_int 0 3 (by norm_num))]
  refine setIntegral_congr_fun measurableSet_Ioi (fun s _ => ?_)
  ring

lemma I02 : I 0 2 = z3 / 4 := by
  have := J2_split; rw [J_values.2.1, I11] at this; linarith

lemma I0312 : 2 * I 0 3 + 3 * I 1 2 = 3 / 4 * (π ^ 4 / 90) := by
  have := J3_split; rw [J_values.2.2, I21] at this; linarith

end Y04

theorem Y04.proof_4 : Y04.stmt_4 :=
  ⟨fun m j hj => Y04.I_int m j hj, Y04.I01, Y04.I11, Y04.I21, Y04.I02, Y04.I0312⟩

theorem Y04.proof_7 : Y04.stmt_7 :=
  ⟨fun _ h1 h3 => Y04.Jint_int h1 h3, Y04.J_values⟩

namespace Y04

/-- Taylor's formula at `0` to order 2 with a cubic remainder on `[0,∞)`. -/
lemma taylor2 {Gf : ℝ → ℝ} (hG : ContDiff ℝ (⊤ : ℕ∞) Gf) {M : ℝ}
    (hM : ∀ x : ℝ, 0 ≤ x → |Gf x| ≤ M) :
    ∃ C : ℝ, 0 ≤ C ∧ ∀ x : ℝ, 0 ≤ x →
      |Gf x - (iteratedDeriv 0 Gf 0 / (Nat.factorial 0 : ℝ) + iteratedDeriv 1 Gf 0 / (Nat.factorial 1 : ℝ) * x
        + iteratedDeriv 2 Gf 0 / (Nat.factorial 2 : ℝ) * x ^ 2)| ≤ C * x ^ 3 := by
  set g0 := iteratedDeriv 0 Gf 0 / (Nat.factorial 0 : ℝ)
  set g1 := iteratedDeriv 1 Gf 0 / (Nat.factorial 1 : ℝ)
  set g2 := iteratedDeriv 2 Gf 0 / (Nat.factorial 2 : ℝ)
  have hc3 : Continuous (iteratedDeriv 3 Gf) :=
    hG.continuous_iteratedDeriv 3 (WithTop.coe_le_coe.2 le_top)
  obtain ⟨K3, hK3⟩ := (isCompact_Icc (a := (0 : ℝ)) (b := 1)).exists_bound_of_continuousOn
    hc3.continuousOn
  have hM0 : 0 ≤ M := (abs_nonneg _).trans (hM 0 le_rfl)
  refine ⟨max (|K3| / 6) (M + |g0| + |g1| + |g2|), le_max_of_le_left (by positivity), ?_⟩
  intro x hx
  rcases eq_or_lt_of_le hx with h0 | hpos
  · subst h0; simp [g0]
  rcases le_or_gt x 1 with hx1 | hx1
  · -- Lagrange remainder
    have hcd : ContDiffOn ℝ (2 + 1) Gf (uIcc 0 x) := (hG.of_le (WithTop.coe_le_coe.2 le_top)).contDiffOn
    obtain ⟨x', hx', heq⟩ := taylor_mean_remainder_lagrange_iteratedDeriv (n := 2) hpos.ne hcd
    have hu : UniqueDiffOn ℝ (uIcc 0 x) := uniqueDiffOn_uIcc hpos.ne
    have hmem : (0 : ℝ) ∈ uIcc 0 x := left_mem_uIcc
    have hT : taylorWithinEval Gf 2 (uIcc 0 x) 0 x = g0 + g1 * x + g2 * x ^ 2 := by
      rw [taylor_within_apply]
      simp only [Finset.sum_range_succ, Finset.sum_range_zero, zero_add, sub_zero, smul_eq_mul]
      rw [iteratedDerivWithin_eq_iteratedDeriv hu (hG.contDiffAt.of_le (WithTop.coe_le_coe.2 le_top)) hmem,
        iteratedDerivWithin_eq_iteratedDeriv hu (hG.contDiffAt.of_le (WithTop.coe_le_coe.2 le_top)) hmem,
        iteratedDerivWithin_eq_iteratedDeriv hu (hG.contDiffAt.of_le (WithTop.coe_le_coe.2 le_top)) hmem]
      simp only [g0, g1, g2]; norm_num; ring
    rw [hT] at heq
    rw [heq]
    rw [uIoo_of_le hx] at hx'
    have hb := hK3 x' ⟨hx'.1.le, by linarith [hx'.2]⟩
    rw [Real.norm_eq_abs] at hb
    rw [sub_zero, abs_div, abs_mul, abs_pow, abs_of_pos hpos]
    norm_num [Nat.factorial]
    have : |iteratedDeriv 3 Gf x'| * x ^ 3 / 6 ≤ |K3| / 6 * x ^ 3 := by
      have := hb.trans (le_abs_self K3)
      have h3 : 0 ≤ x ^ 3 := by positivity
      nlinarith
    calc |iteratedDeriv 3 Gf x'| * x ^ 3 / 6 ≤ |K3| / 6 * x ^ 3 := this
      _ ≤ _ := by gcongr; exact le_max_left _ _
  · -- crude bound for x > 1
    have h1 : 1 ≤ x := hx1.le
    have hx2' : 1 ≤ x ^ 2 := by nlinarith
    have hx22 : x ^ 2 ≤ x ^ 3 := by
      rw [show x ^ 3 = x * x ^ 2 by ring]; nlinarith
    have hx2 : x ≤ x ^ 3 := by
      rw [show x ^ 3 = x * x ^ 2 by ring]; nlinarith
    have h13 : 1 ≤ x ^ 3 := by linarith
    calc _ ≤ |Gf x| + (|g0| + |g1| * x + |g2| * x ^ 2) := by
          have := abs_sub (Gf x) (g0 + g1 * x + g2 * x ^ 2)
          have a1 := abs_add_le (g0 + g1 * x) (g2 * x ^ 2)
          have a2 := abs_add_le g0 (g1 * x)
          rw [abs_mul, abs_of_pos hpos] at a2
          rw [abs_mul, abs_of_nonneg (by positivity : (0:ℝ) ≤ x ^ 2)] at a1
          linarith
      _ ≤ (M + |g0| + |g1| + |g2|) * x ^ 3 := by
          have := hM x hx
          have e1 : |g1| * x ≤ |g1| * x ^ 3 := mul_le_mul_of_nonneg_left hx2 (abs_nonneg _)
          have e2 : |g2| * x ^ 2 ≤ |g2| * x ^ 3 := mul_le_mul_of_nonneg_left hx22 (abs_nonneg _)
          have e3 : M ≤ M * x ^ 3 := le_mul_of_one_le_right hM0 h13
          have e4 : |g0| ≤ |g0| * x ^ 3 := le_mul_of_one_le_right (abs_nonneg _) h13
          nlinarith
      _ ≤ _ := by gcongr; exact le_max_right _ _

/-- the dominating function of the Watson remainder -/
noncomputable def B5 (s : ℝ) : ℝ := (1 + s + s ^ 2 + s ^ 3) * Real.exp (-s)

lemma B5_int : IntegrableOn B5 (Set.Ioi 0) := by
  have h := (((intOn_pow_exp1 0).add (intOn_pow_exp1 1)).add (intOn_pow_exp1 2)).add (intOn_pow_exp1 3)
  refine IntegrableOn.congr_fun h (fun s _ => ?_) measurableSet_Ioi
  simp only [B5, Pi.add_apply]; ring

lemma B5_nonneg {s : ℝ} (hs : 0 ≤ s) : 0 ≤ B5 s := by unfold B5; positivity

lemma rpow_eq_exp (τ s : ℝ) : (1 + Real.exp (-s)) ^ τ = Real.exp (τ * ell s) := by
  rw [Real.rpow_def_of_pos (by positivity), ell, mul_comm]

/-- the pointwise Watson remainder -/
lemma watson_pt {g0 g1 g2 CR M τ s : ℝ} (hCR : 0 ≤ CR) (hM : 0 ≤ M) (hτ0 : 0 < τ) (hτ1 : τ < 1)
    (hs : 0 ≤ s) {Gv : ℝ} (hGv : |Gv| ≤ M)
    (hR : |Gv - (g0 + g1 * (τ * s) + g2 * (τ * s) ^ 2)| ≤ CR * (τ * s) ^ 3) :
    |τ * ((Real.exp (τ * ell s) - 1) * Gv)
      - (τ ^ 2 * g0 * (s ^ 0 * ell s ^ 1) + τ ^ 3 * g1 * (s ^ 1 * ell s ^ 1)
        + τ ^ 3 * (g0 / 2) * (s ^ 0 * ell s ^ 2) + τ ^ 4 * g2 * (s ^ 2 * ell s ^ 1)
        + τ ^ 4 * (g1 / 2) * (s ^ 1 * ell s ^ 2) + τ ^ 4 * (g0 / 6) * (s ^ 0 * ell s ^ 3))|
      ≤ τ ^ 5 * ((|g1| + |g2| + 2 * CR + M) * B5 s) := by
  have hl0 := ell_nonneg s
  have hl1 := ell_le_one s hs
  have hle := ell_le s
  set L := ell s
  have hτL : |τ * L| ≤ 1 := by rw [abs_of_nonneg (by positivity)]; nlinarith
  have he := Real.exp_bound hτL (n := 4) (by norm_num)
  simp only [Finset.sum_range_succ, Finset.sum_range_zero, Nat.factorial] at he
  norm_num at he
  set a1 := τ * L + (τ * L) ^ 2 / 2 + (τ * L) ^ 3 / 6 with ha1
  set r1 := Real.exp (τ * L) - 1 - a1 with hr1
  have hr1b : |r1| ≤ τ ^ 4 * L ^ 4 := by
    have : r1 = Real.exp (τ * L) - (1 + τ * L + (τ * L) ^ 2 / 2 + (τ * L) ^ 3 / 6) := by
      rw [hr1, ha1]; ring
    rw [this]
    have h4 : (|τ| * |L|) ^ 4 = τ ^ 4 * L ^ 4 := by
      rw [abs_of_pos hτ0, abs_of_nonneg hl0]; ring
    rw [h4] at he
    have : τ ^ 4 * L ^ 4 * (5 / 96) ≤ τ ^ 4 * L ^ 4 := by
      have : 0 ≤ τ ^ 4 * L ^ 4 := by positivity
      linarith
    linarith
  set q := g0 + g1 * (τ * s) + g2 * (τ * s) ^ 2
  set R := Gv - q with hRdef
  have key : τ * ((Real.exp (τ * L) - 1) * Gv)
      - (τ ^ 2 * g0 * (s ^ 0 * L ^ 1) + τ ^ 3 * g1 * (s ^ 1 * L ^ 1)
        + τ ^ 3 * (g0 / 2) * (s ^ 0 * L ^ 2) + τ ^ 4 * g2 * (s ^ 2 * L ^ 1)
        + τ ^ 4 * (g1 / 2) * (s ^ 1 * L ^ 2) + τ ^ 4 * (g0 / 6) * (s ^ 0 * L ^ 3))
      = τ ^ 5 * (g2 * s ^ 2 * L ^ 2 / 2 + g1 * s * L ^ 3 / 6 + τ * g2 * s ^ 2 * L ^ 3 / 6)
        + τ * a1 * R + τ * r1 * Gv := by
    rw [hRdef, hr1, ha1]; ring
  rw [key]
  have hes : 0 ≤ Real.exp (-s) := (Real.exp_pos _).le
  -- term 1
  have hL2 : L ^ 2 ≤ Real.exp (-s) := (pow_le_of_le_one hl0 hl1 (by norm_num)).trans hle
  have hL3 : L ^ 3 ≤ Real.exp (-s) := (pow_le_of_le_one hl0 hl1 (by norm_num)).trans hle
  have hL4 : L ^ 4 ≤ Real.exp (-s) := (pow_le_of_le_one hl0 hl1 (by norm_num)).trans hle
  have hτL3 : τ * L ^ 3 ≤ Real.exp (-s) := by
    have : τ * L ^ 3 ≤ 1 * L ^ 3 := mul_le_mul_of_nonneg_right hτ1.le (by positivity)
    linarith
  have b1 : |τ ^ 5 * (g2 * s ^ 2 * L ^ 2 / 2 + g1 * s * L ^ 3 / 6 + τ * g2 * s ^ 2 * L ^ 3 / 6)|
      ≤ τ ^ 5 * ((|g1| + |g2|) * (s + s ^ 2) * Real.exp (-s)) := by
    rw [abs_mul, abs_of_pos (by positivity)]
    apply mul_le_mul_of_nonneg_left _ (by positivity)
    have e : g2 * s ^ 2 * L ^ 2 / 2 + g1 * s * L ^ 3 / 6 + τ * g2 * s ^ 2 * L ^ 3 / 6
        = g2 * (s ^ 2 * (L ^ 2 / 2 + τ * L ^ 3 / 6)) + g1 * (s * (L ^ 3 / 6)) := by ring
    rw [e]
    have u1 : 0 ≤ s ^ 2 * (L ^ 2 / 2 + τ * L ^ 3 / 6) := by positivity
    have u1' : s ^ 2 * (L ^ 2 / 2 + τ * L ^ 3 / 6) ≤ s ^ 2 * Real.exp (-s) :=
      mul_le_mul_of_nonneg_left (by linarith) (sq_nonneg s)
    have u2 : 0 ≤ s * (L ^ 3 / 6) := by positivity
    have u2' : s * (L ^ 3 / 6) ≤ s * Real.exp (-s) := mul_le_mul_of_nonneg_left (by linarith) hs
    have v1 : |g2 * (s ^ 2 * (L ^ 2 / 2 + τ * L ^ 3 / 6))| ≤ |g2| * (s ^ 2 * Real.exp (-s)) := by
      rw [abs_mul, abs_of_nonneg u1]; exact mul_le_mul_of_nonneg_left u1' (abs_nonneg _)
    have v2 : |g1 * (s * (L ^ 3 / 6))| ≤ |g1| * (s * Real.exp (-s)) := by
      rw [abs_mul, abs_of_nonneg u2]; exact mul_le_mul_of_nonneg_left u2' (abs_nonneg _)
    have := abs_add_le (g2 * (s ^ 2 * (L ^ 2 / 2 + τ * L ^ 3 / 6))) (g1 * (s * (L ^ 3 / 6)))
    have w1 : |g2| * (s ^ 2 * Real.exp (-s)) ≤ |g2| * ((s + s ^ 2) * Real.exp (-s)) :=
      mul_le_mul_of_nonneg_left (mul_le_mul_of_nonneg_right (by linarith) hes) (abs_nonneg _)
    have w2 : |g1| * (s * Real.exp (-s)) ≤ |g1| * ((s + s ^ 2) * Real.exp (-s)) :=
      mul_le_mul_of_nonneg_left (mul_le_mul_of_nonneg_right (by nlinarith) hes) (abs_nonneg _)
    have : (|g1| + |g2|) * (s + s ^ 2) * Real.exp (-s)
        = |g1| * ((s + s ^ 2) * Real.exp (-s)) + |g2| * ((s + s ^ 2) * Real.exp (-s)) := by ring
    linarith
  -- term 2
  have b2 : |τ * a1 * R| ≤ τ ^ 5 * (2 * CR * s ^ 3 * Real.exp (-s)) := by
    have ha1b : |a1| ≤ 2 * τ * L := by
      have hu0 : 0 ≤ τ * L := mul_nonneg hτ0.le hl0
      have hu1 : τ * L ≤ 1 := by
        have := mul_le_mul hτ1.le hl1 hl0 zero_le_one; linarith
      have p2 : (τ * L) ^ 2 ≤ τ * L := pow_le_of_le_one hu0 hu1 (by norm_num)
      have p3 : (τ * L) ^ 3 ≤ τ * L := pow_le_of_le_one hu0 hu1 (by norm_num)
      have q2 : 0 ≤ (τ * L) ^ 2 := sq_nonneg _
      have q3 : 0 ≤ (τ * L) ^ 3 := pow_nonneg hu0 3
      have : 0 ≤ a1 := by rw [ha1]; linarith
      rw [abs_of_nonneg this, ha1]
      have : 2 * τ * L = 2 * (τ * L) := by ring
      linarith
    have hRb : |R| ≤ CR * (τ * s) ^ 3 := hR
    rw [abs_mul, abs_mul, abs_of_pos hτ0]
    have hτ' : 0 ≤ τ := hτ0.le
    calc τ * |a1| * |R| ≤ τ * (2 * τ * L) * (CR * (τ * s) ^ 3) :=
          mul_le_mul (mul_le_mul_of_nonneg_left ha1b hτ') hRb (abs_nonneg _) (by positivity)
      _ = τ ^ 5 * (2 * CR * s ^ 3 * L) := by ring
      _ ≤ τ ^ 5 * (2 * CR * s ^ 3 * Real.exp (-s)) :=
          mul_le_mul_of_nonneg_left (mul_le_mul_of_nonneg_left hle (by positivity)) (by positivity)
  -- term 3
  have b3 : |τ * r1 * Gv| ≤ τ ^ 5 * (M * Real.exp (-s)) := by
    rw [abs_mul, abs_mul, abs_of_pos hτ0]
    calc τ * |r1| * |Gv| ≤ τ * (τ ^ 4 * L ^ 4) * M :=
          mul_le_mul (mul_le_mul_of_nonneg_left hr1b hτ0.le) hGv (abs_nonneg _) (by positivity)
      _ = τ ^ 5 * (L ^ 4 * M) := by ring
      _ ≤ τ ^ 5 * (Real.exp (-s) * M) :=
          mul_le_mul_of_nonneg_left (mul_le_mul_of_nonneg_right hL4 hM) (by positivity)
      _ = _ := by ring
  have t1 := abs_add_le (τ ^ 5 * (g2 * s ^ 2 * L ^ 2 / 2 + g1 * s * L ^ 3 / 6 + τ * g2 * s ^ 2 * L ^ 3 / 6)
    + τ * a1 * R) (τ * r1 * Gv)
  have t2 := abs_add_le (τ ^ 5 * (g2 * s ^ 2 * L ^ 2 / 2 + g1 * s * L ^ 3 / 6 + τ * g2 * s ^ 2 * L ^ 3 / 6))
    (τ * a1 * R)
  have hfin : τ ^ 5 * ((|g1| + |g2|) * (s + s ^ 2) * Real.exp (-s)) + τ ^ 5 * (2 * CR * s ^ 3 * Real.exp (-s))
      + τ ^ 5 * (M * Real.exp (-s)) ≤ τ ^ 5 * ((|g1| + |g2| + 2 * CR + M) * B5 s) := by
    unfold B5
    have h5 : 0 ≤ τ ^ 5 := by positivity
    apply le_of_sub_nonneg
    have : τ ^ 5 * ((|g1| + |g2| + 2 * CR + M) * ((1 + s + s ^ 2 + s ^ 3) * Real.exp (-s)))
        - (τ ^ 5 * ((|g1| + |g2|) * (s + s ^ 2) * Real.exp (-s)) + τ ^ 5 * (2 * CR * s ^ 3 * Real.exp (-s))
        + τ ^ 5 * (M * Real.exp (-s)))
        = τ ^ 5 * Real.exp (-s) * ((|g1| + |g2|) * (1 + s ^ 3) + 2 * CR * (1 + s + s ^ 2)
          + M * (s + s ^ 2 + s ^ 3)) := by ring
    rw [this]
    have := abs_nonneg g1; have := abs_nonneg g2
    positivity
  linarith

end Y04

namespace Y04

lemma int6 {f1 f2 f3 f4 f5 f6 : ℝ → ℝ} (c1 c2 c3 c4 c5 c6 : ℝ)
    (h1 : IntegrableOn f1 (Set.Ioi 0)) (h2 : IntegrableOn f2 (Set.Ioi 0)) (h3 : IntegrableOn f3 (Set.Ioi 0))
    (h4 : IntegrableOn f4 (Set.Ioi 0)) (h5 : IntegrableOn f5 (Set.Ioi 0)) (h6 : IntegrableOn f6 (Set.Ioi 0)) :
    IntegrableOn (fun s => c1 * f1 s + c2 * f2 s + c3 * f3 s + c4 * f4 s + c5 * f5 s + c6 * f6 s) (Set.Ioi 0) ∧
    ∫ s in Set.Ioi (0 : ℝ), (c1 * f1 s + c2 * f2 s + c3 * f3 s + c4 * f4 s + c5 * f5 s + c6 * f6 s)
      = c1 * (∫ s in Set.Ioi (0 : ℝ), f1 s) + c2 * (∫ s in Set.Ioi (0 : ℝ), f2 s)
        + c3 * (∫ s in Set.Ioi (0 : ℝ), f3 s) + c4 * (∫ s in Set.Ioi (0 : ℝ), f4 s)
        + c5 * (∫ s in Set.Ioi (0 : ℝ), f5 s) + c6 * (∫ s in Set.Ioi (0 : ℝ), f6 s) := by
  have a1 : Integrable (fun s => c1 * f1 s) (volume.restrict (Set.Ioi 0)) := h1.const_mul c1
  have a2 : Integrable (fun s => c2 * f2 s) (volume.restrict (Set.Ioi 0)) := h2.const_mul c2
  have a3 : Integrable (fun s => c3 * f3 s) (volume.restrict (Set.Ioi 0)) := h3.const_mul c3
  have a4 : Integrable (fun s => c4 * f4 s) (volume.restrict (Set.Ioi 0)) := h4.const_mul c4
  have a5 : Integrable (fun s => c5 * f5 s) (volume.restrict (Set.Ioi 0)) := h5.const_mul c5
  have a6 : Integrable (fun s => c6 * f6 s) (volume.restrict (Set.Ioi 0)) := h6.const_mul c6
  have b2 : Integrable (fun s => c1 * f1 s + c2 * f2 s) (volume.restrict (Set.Ioi 0)) := a1.add a2
  have b3 : Integrable (fun s => c1 * f1 s + c2 * f2 s + c3 * f3 s) (volume.restrict (Set.Ioi 0)) :=
    b2.add a3
  have b4 : Integrable (fun s => c1 * f1 s + c2 * f2 s + c3 * f3 s + c4 * f4 s)
    (volume.restrict (Set.Ioi 0)) := b3.add a4
  have b5 : Integrable (fun s => c1 * f1 s + c2 * f2 s + c3 * f3 s + c4 * f4 s + c5 * f5 s)
    (volume.restrict (Set.Ioi 0)) := b4.add a5
  have b6 : Integrable (fun s => c1 * f1 s + c2 * f2 s + c3 * f3 s + c4 * f4 s + c5 * f5 s + c6 * f6 s)
    (volume.restrict (Set.Ioi 0)) := b5.add a6
  refine ⟨b6, ?_⟩
  rw [integral_add b5 a6, integral_add b4 a5, integral_add b3 a4, integral_add b2 a3,
    integral_add a1 a2, integral_const_mul, integral_const_mul, integral_const_mul,
    integral_const_mul, integral_const_mul, integral_const_mul]

lemma watson (Gf : ℝ → ℝ) (hG : ContDiff ℝ (⊤ : ℕ∞) Gf) (hb : ∃ M : ℝ, ∀ x : ℝ, 0 ≤ x → |Gf x| ≤ M) :
    ∃ C : ℝ, ∀ τ : ℝ, 0 < τ → τ < 1 →
      IntegrableOn (fun s : ℝ => ((1 + Real.exp (-s)) ^ τ - 1) * Gf (τ * s)) (Set.Ioi 0) ∧
      |τ * (∫ s in Set.Ioi (0 : ℝ), ((1 + Real.exp (-s)) ^ τ - 1) * Gf (τ * s)) -
          ∑ n ∈ Finset.Icc 2 4, τ ^ n * ∑ j ∈ Finset.Icc 1 (n - 1),
            iteratedDeriv (n - 1 - j) Gf 0 / ((n - 1 - j).factorial : ℝ) /
              (j.factorial : ℝ) * I (n - 1 - j) j| ≤ C * τ ^ 5 := by
  obtain ⟨M, hM⟩ := hb
  have hM0 : 0 ≤ M := (abs_nonneg _).trans (hM 0 le_rfl)
  obtain ⟨CR, hCR, hT⟩ := taylor2 hG hM
  set g0 := iteratedDeriv 0 Gf 0 / (Nat.factorial 0 : ℝ) with hg0
  set g1 := iteratedDeriv 1 Gf 0 / (Nat.factorial 1 : ℝ) with hg1
  set g2 := iteratedDeriv 2 Gf 0 / (Nat.factorial 2 : ℝ) with hg2
  refine ⟨(|g1| + |g2| + 2 * CR + M) * ∫ s in Set.Ioi (0 : ℝ), B5 s, fun τ hτ0 hτ1 => ?_⟩
  have hGc : Continuous Gf := hG.continuous
  have e : (fun s : ℝ => ((1 + Real.exp (-s)) ^ τ - 1) * Gf (τ * s))
      = fun s => (Real.exp (τ * ell s) - 1) * Gf (τ * s) := by
    funext s; rw [rpow_eq_exp]
  rw [e]
  have hint : IntegrableOn (fun s => (Real.exp (τ * ell s) - 1) * Gf (τ * s)) (Set.Ioi 0) := by
    refine intOn_of_bound (by have := continuous_ell; fun_prop) ((intOn_pow_exp1 0).const_mul (2 * M))
      (fun s hs => ?_)
    have hl0 := ell_nonneg s
    have hl1 := ell_le_one s hs.le
    have hτL : |τ * ell s| ≤ 1 := by
      rw [abs_of_nonneg (mul_nonneg hτ0.le hl0)]
      have := mul_le_mul hτ1.le hl1 hl0 zero_le_one; linarith
    have h1 := Real.abs_exp_sub_one_le hτL
    have h2 : |τ * ell s| ≤ Real.exp (-s) := by
      rw [abs_of_nonneg (mul_nonneg hτ0.le hl0)]
      have := mul_le_mul hτ1.le (ell_le s) hl0 zero_le_one; linarith
    have h3 := hM (τ * s) (by nlinarith)
    rw [abs_mul]
    simp only [pow_zero, one_mul]
    calc |Real.exp (τ * ell s) - 1| * |Gf (τ * s)| ≤ (2 * Real.exp (-s)) * M :=
          mul_le_mul (by linarith) h3 (abs_nonneg _) (by positivity)
      _ = 2 * M * Real.exp (-s) := by ring
  refine ⟨hint, ?_⟩
  -- the main terms
  obtain ⟨hmi, hmv⟩ := int6 (τ ^ 2 * g0) (τ ^ 3 * g1) (τ ^ 3 * (g0 / 2)) (τ ^ 4 * g2) (τ ^ 4 * (g1 / 2))
    (τ ^ 4 * (g0 / 6)) (I_int 0 1 le_rfl) (I_int 1 1 le_rfl) (I_int 0 2 (by norm_num))
    (I_int 2 1 le_rfl) (I_int 1 2 (by norm_num)) (I_int 0 3 (by norm_num))
  have hsum : ∑ n ∈ Finset.Icc 2 4, τ ^ n * ∑ j ∈ Finset.Icc 1 (n - 1),
            iteratedDeriv (n - 1 - j) Gf 0 / ((n - 1 - j).factorial : ℝ) /
              (j.factorial : ℝ) * I (n - 1 - j) j
      = τ ^ 2 * g0 * I 0 1 + τ ^ 3 * g1 * I 1 1 + τ ^ 3 * (g0 / 2) * I 0 2 + τ ^ 4 * g2 * I 2 1
        + τ ^ 4 * (g1 / 2) * I 1 2 + τ ^ 4 * (g0 / 6) * I 0 3 := by
    rw [show Finset.Icc 2 4 = {2, 3, 4} by decide]
    rw [Finset.sum_insert (by decide), Finset.sum_insert (by decide), Finset.sum_singleton]
    rw [show Finset.Icc 1 (2 - 1) = {1} by decide, show Finset.Icc 1 (3 - 1) = {1, 2} by decide,
      show Finset.Icc 1 (4 - 1) = {1, 2, 3} by decide]
    simp [Finset.sum_insert, hg0, hg1, hg2, Nat.factorial]
    ring
  rw [hsum]
  unfold I
  rw [← hmv, ← integral_const_mul, ← integral_sub (hint.const_mul τ) hmi]
  have hB := (B5_int.const_mul ((|g1| + |g2| + 2 * CR + M))).const_mul (τ ^ 5)
  have := norm_integral_le_of_norm_le hB
    ((ae_restrict_iff' measurableSet_Ioi).2 (Filter.Eventually.of_forall
      (fun s (hs : 0 < s) => by
        rw [Real.norm_eq_abs]
        have hx : 0 ≤ τ * s := by positivity
        have := watson_pt (g0 := g0) (g1 := g1) (g2 := g2) hCR hM0 hτ0 hτ1 hs.le (hM (τ * s) hx)
          (by have := hT (τ * s) hx; simpa [hg0, hg1, hg2] using this)
        convert this using 2)))
  rw [Real.norm_eq_abs, integral_const_mul, integral_const_mul] at this
  calc _ ≤ _ := this
    _ = _ := by ring

end Y04

theorem Y04.proof_5 : Y04.stmt_5 := fun Gf hG hb => Y04.watson Gf hG hb

namespace Y04

lemma tau_pos {p : ℝ} (hp : 1 < p) : 0 < tau p := by
  unfold tau; have : 1 / p < 1 := by rw [div_lt_one (by linarith)]; exact hp
  linarith

lemma tau_lt_one {p : ℝ} (hp : 1 < p) : tau p < 1 := by
  unfold tau; have : 0 < 1 / p := by positivity
  linarith

lemma h_cont {p : ℝ} (hp : 1 < p) : Continuous (h p) := by
  have ht := tau_pos hp
  unfold h
  refine Continuous.rpow_const ?_ (fun _ => Or.inr ht.le)
  exact ((Real.continuous_cos.abs).rpow_const (fun _ => Or.inr (by positivity))).add
    ((Real.continuous_sin.abs).rpow_const (fun _ => Or.inr (by positivity)))

lemma hDiamond_cont : Continuous hDiamond := by unfold hDiamond; fun_prop

lemma h_sym (p θ : ℝ) : h p (π / 2 - θ) = h p θ := by
  unfold h; rw [Real.cos_pi_div_two_sub, Real.sin_pi_div_two_sub, add_comm]

lemma hDiamond_sym (θ : ℝ) : hDiamond (π / 2 - θ) = hDiamond θ := by
  unfold hDiamond; rw [Real.cos_pi_div_two_sub, Real.sin_pi_div_two_sub, max_comm]

lemma cos_refl (k : ℕ) (θ : ℝ) : Real.cos (4 * (k : ℝ) * (π / 2 - θ)) = Real.cos (4 * k * θ) := by
  rw [show 4 * (k : ℝ) * (π / 2 - θ) = ((2 * k : ℕ) : ℝ) * π - 4 * k * θ by push_cast; ring,
    Real.cos_nat_mul_pi_sub, pow_mul, neg_one_sq, one_pow, one_mul]

/-- for a function symmetric about `π/4`, the coefficient is twice the integral over `[0, π/4]` -/
lemma a_sym (k : ℕ) {f : ℝ → ℝ} (hf : Continuous f) (hs : ∀ θ, f (π / 2 - θ) = f θ) :
    a k f = 8 / π * ∫ θ in (0 : ℝ)..(π / 4), f θ * Real.cos (4 * k * θ) := by
  have hpi := Real.pi_pos
  have hc : Continuous (fun θ => f θ * Real.cos (4 * k * θ)) := by fun_prop
  have hsplit := intervalIntegral.integral_add_adjacent_intervals (μ := volume)
    (a := 0) (b := π / 4) (c := π / 2) (hc.intervalIntegrable _ _) (hc.intervalIntegrable _ _)
  have hrefl : ∫ θ in (π / 4)..(π / 2), f θ * Real.cos (4 * k * θ)
      = ∫ θ in (0 : ℝ)..(π / 4), f θ * Real.cos (4 * k * θ) := by
    have := intervalIntegral.integral_comp_sub_left (a := (0 : ℝ)) (b := π / 4)
      (fun θ => f θ * Real.cos (4 * k * θ)) (π / 2)
    rw [show π / 2 - π / 4 = π / 4 by ring, sub_zero] at this
    rw [← this]
    refine intervalIntegral.integral_congr (fun θ _ => ?_)
    show f (π / 2 - θ) * Real.cos (4 * (k : ℝ) * (π / 2 - θ)) = f θ * Real.cos (4 * k * θ)
    rw [hs θ, cos_refl k θ]
  unfold a
  rw [← hsplit, hrefl]; ring

/-- `θ_x` on the half line: `thetaX (τ s) ∈ (0, π/4)` for `s > 0` -/
lemma thetaS_mem {τ s : ℝ} (hτ : 0 < τ) (hs : 0 < s) :
    0 < thetaX (τ * s) ∧ thetaX (τ * s) < π / 4 := by
  unfold thetaX
  constructor
  · exact Real.arctan_pos.2 (Real.exp_pos _)
  · rw [← Real.arctan_one]
    exact Real.arctan_strictMono (Real.exp_lt_one_iff.2 (by nlinarith))

lemma thetaS_image {τ : ℝ} (hτ : 0 < τ) :
    (fun s => thetaX (τ * s)) '' Set.Ioi 0 = Set.Ioo 0 (π / 4) := by
  ext θ; constructor
  · rintro ⟨s, hs, rfl⟩; exact thetaS_mem hτ hs
  · rintro ⟨h0, h1⟩
    have hpi := Real.pi_pos
    have ht0 : 0 < Real.tan θ := Real.tan_pos_of_pos_of_lt_pi_div_two h0 (by linarith)
    have ht1 : Real.tan θ < 1 := by
      rw [← Real.tan_pi_div_four]
      exact Real.tan_lt_tan_of_nonneg_of_lt_pi_div_two h0.le (by linarith) h1
    refine ⟨-Real.log (Real.tan θ) / τ, ?_, ?_⟩
    · have : Real.log (Real.tan θ) < 0 := Real.log_neg ht0 ht1
      simp only [Set.mem_Ioi]; apply div_pos (by linarith) hτ
    · simp only [thetaX]
      rw [show -(τ * (-Real.log (Real.tan θ) / τ)) = Real.log (Real.tan θ) by field_simp,
        Real.exp_log ht0, Real.arctan_tan (by linarith) (by linarith)]

lemma thetaS_hasDeriv (τ s : ℝ) :
    HasDerivAt (fun s => thetaX (τ * s))
      (-(Real.sin (thetaX (τ * s)) * Real.cos (thetaX (τ * s))) * τ) s := by
  have h : HasDerivAt (fun s => thetaX (τ * s)) _ s :=
    (theta_hasDeriv (τ * s)).comp s ((hasDerivAt_id s).const_mul τ)
  refine h.congr_deriv ?_
  simp

lemma thetaS_anti {τ : ℝ} (hτ : 0 < τ) : StrictAnti (fun s => thetaX (τ * s)) := by
  intro s t hst
  simp only [thetaX]
  exact Real.arctan_strictMono (Real.exp_lt_exp.2 (by nlinarith))

/-- the pointwise identity behind the substitution -/
lemma h_layer {p s : ℝ} (hp : 1 < p) (hs : 0 < s) :
    h p (thetaX (tau p * s)) - hDiamond (thetaX (tau p * s))
      = Real.cos (thetaX (tau p * s)) * ((1 + Real.exp (-s)) ^ tau p - 1) := by
  have hτ := tau_pos hp
  obtain ⟨h0, h1⟩ := thetaS_mem hτ hs
  set θ := thetaX (tau p * s) with hθ
  have hpi := Real.pi_pos
  have hc : 0 < Real.cos θ := Real.cos_pos_of_mem_Ioo ⟨by linarith, by linarith⟩
  have hsn : 0 < Real.sin θ := Real.sin_pos_of_pos_of_lt_pi h0 (by linarith)
  have hD : hDiamond θ = Real.cos θ := by
    unfold hDiamond
    rw [abs_of_pos hc, abs_of_pos hsn, max_eq_left (sin_le_cos_of h0.le h1.le)]
  rw [hD]
  -- tan θ = e^{-τ s}
  have htan : Real.sin θ = Real.exp (-(tau p * s)) * Real.cos θ := by
    have : Real.tan θ = Real.exp (-(tau p * s)) := by
      rw [hθ]; unfold thetaX; rw [Real.tan_arctan]
    rw [← this, Real.tan_eq_sin_div_cos]; field_simp
  unfold h
  rw [abs_of_pos hc, abs_of_pos hsn, htan, Real.mul_rpow (Real.exp_pos _).le hc.le]
  rw [show Real.cos θ ^ (1 / tau p) + Real.exp (-(tau p * s)) ^ (1 / tau p) * Real.cos θ ^ (1 / tau p)
      = Real.cos θ ^ (1 / tau p) * (1 + Real.exp (-(tau p * s)) ^ (1 / tau p)) by ring]
  rw [Real.mul_rpow (by positivity) (by positivity), ← Real.rpow_mul hc.le,
    show 1 / tau p * tau p = 1 by field_simp, Real.rpow_one, ← Real.exp_mul,
    show -(tau p * s) * (1 / tau p) = -s by field_simp]
  ring

lemma layer_formula (k : ℕ) {p : ℝ} (hp : 1 < p) :
    a k (h p) - a k hDiamond =
      8 / π * tau p * ∫ s in Set.Ioi (0 : ℝ), ((1 + Real.exp (-s)) ^ tau p - 1) * G k (tau p * s) := by
  have hτ := tau_pos hp
  have hpi := Real.pi_pos
  set τ := tau p with hτdef
  have hd : Continuous (fun θ => h p θ - hDiamond θ) := (h_cont hp).sub hDiamond_cont
  have hsub : a k (h p) - a k hDiamond = a k (fun θ => h p θ - hDiamond θ) := by
    unfold a
    rw [← mul_sub, ← intervalIntegral.integral_sub]
    · congr 1; refine intervalIntegral.integral_congr (fun θ _ => ?_); simp only; ring
    · exact ((h_cont hp).mul (by fun_prop)).intervalIntegrable _ _
    · exact (hDiamond_cont.mul (by fun_prop)).intervalIntegrable _ _
  rw [hsub, a_sym k hd (fun θ => by simp only [h_sym, hDiamond_sym])]
  -- the substitution
  set g : ℝ → ℝ := fun θ => (h p θ - hDiamond θ) * Real.cos (4 * k * θ) with hg
  have hder : ∀ s ∈ Set.Ioi (0 : ℝ), HasDerivWithinAt (fun s => thetaX (τ * s))
      (-(Real.sin (thetaX (τ * s)) * Real.cos (thetaX (τ * s))) * τ) (Set.Ioi 0) s :=
    fun s _ => (thetaS_hasDeriv τ s).hasDerivWithinAt
  have hinj : Set.InjOn (fun s => thetaX (τ * s)) (Set.Ioi 0) := (thetaS_anti hτ).injective.injOn
  have himg := integral_image_eq_integral_abs_deriv_smul measurableSet_Ioi hder hinj g
  rw [thetaS_image hτ] at himg
  have hset : ∫ θ in (0 : ℝ)..(π / 4), g θ = ∫ θ in Set.Ioo (0 : ℝ) (π / 4), g θ := by
    rw [intervalIntegral.integral_of_le (by linarith : (0 : ℝ) ≤ π / 4), integral_Ioc_eq_integral_Ioo]
  rw [show (∫ θ in (0 : ℝ)..(π / 4), (h p θ - hDiamond θ) * Real.cos (4 * k * θ))
      = ∫ θ in (0 : ℝ)..(π / 4), g θ from rfl, hset, himg, mul_assoc]
  congr 1
  rw [← integral_const_mul]
  refine setIntegral_congr_fun (measurableSet_Ioi : MeasurableSet (Set.Ioi (0 : ℝ))) (fun s hs => ?_)
  obtain ⟨h0, h1⟩ := thetaS_mem hτ hs
  have hc : 0 < Real.cos (thetaX (τ * s)) := Real.cos_pos_of_mem_Ioo ⟨by linarith, by linarith⟩
  have hsn : 0 < Real.sin (thetaX (τ * s)) := Real.sin_pos_of_pos_of_lt_pi h0 (by linarith)
  simp only [hg, smul_eq_mul]
  rw [hτdef, h_layer hp hs, ← hτdef]
  rw [abs_of_neg (by nlinarith [mul_pos hsn hc])]
  unfold G
  ring

end Y04

theorem Y04.proof_2 : Y04.stmt_2 := by
  intro k p hp
  refine ⟨?_, Y04.layer_formula k hp⟩
  obtain ⟨C, hC⟩ := Y04.watson (Y04.G k) (Y04.G_contDiff k) ⟨1, fun x _ => Y04.G_bdd k x⟩
  exact (hC (Y04.tau p) (Y04.tau_pos hp) (Y04.tau_lt_one hp)).1

namespace Y04

lemma sum_eval (Gf : ℝ → ℝ) (τ : ℝ) :
    ∑ n ∈ Finset.Icc 2 4, τ ^ n * ∑ j ∈ Finset.Icc 1 (n - 1),
        iteratedDeriv (n - 1 - j) Gf 0 / ((n - 1 - j).factorial : ℝ) / (j.factorial : ℝ) * I (n - 1 - j) j
      = τ ^ 2 * Gf 0 * I 0 1 + τ ^ 3 * deriv Gf 0 * I 1 1 + τ ^ 3 * (Gf 0 / 2) * I 0 2
        + τ ^ 4 * (iteratedDeriv 2 Gf 0 / 2) * I 2 1
        + τ ^ 4 * (deriv Gf 0 / 2) * I 1 2 + τ ^ 4 * (Gf 0 / 6) * I 0 3 := by
  rw [show Finset.Icc 2 4 = {2, 3, 4} by decide]
  rw [Finset.sum_insert (by decide), Finset.sum_insert (by decide), Finset.sum_singleton]
  rw [show Finset.Icc 1 (2 - 1) = {1} by decide, show Finset.Icc 1 (3 - 1) = {1, 2} by decide,
    show Finset.Icc 1 (4 - 1) = {1, 2, 3} by decide]
  simp [Finset.sum_insert, Nat.factorial]
  ring

lemma fourier_expansion (k : ℕ) : ∃ C τ₀ : ℝ, 0 < τ₀ ∧ ∀ p : ℝ, 1 < p → tau p ≤ τ₀ →
    |a k (h p) - (-1 : ℝ) ^ k * (Real.sqrt 2 / π) *
        (-(4 / (16 * (k : ℝ) ^ 2 - 1)) + (π ^ 2 / 6) * tau p ^ 2 + z3 * tau p ^ 3
          - (33 / 16 + 7 * (k : ℝ) ^ 2) * (π ^ 4 / 90) * tau p ^ 4)|
      ≤ C * tau p ^ 5 := by
  obtain ⟨CW, hW⟩ := watson (G k) (G_contDiff k) ⟨1, fun x _ => G_bdd k x⟩
  refine ⟨8 / π * CW, 1 / 2, by norm_num, fun p hp _ => ?_⟩
  have hpi := Real.pi_pos
  have hτ0 := tau_pos hp
  have hτ1 := tau_lt_one hp
  obtain ⟨_, hb⟩ := hW (tau p) hτ0 hτ1
  rw [sum_eval] at hb
  have h2 := layer_formula k hp
  have h1 := A_formula k
  set X := tau p * ∫ s in Set.Ioi (0 : ℝ), ((1 + Real.exp (-s)) ^ tau p - 1) * G k (tau p * s) with hX
  set τ := tau p with hτ
  set S := τ ^ 2 * G k 0 * I 0 1 + τ ^ 3 * deriv (G k) 0 * I 1 1 + τ ^ 3 * (G k 0 / 2) * I 0 2
        + τ ^ 4 * (iteratedDeriv 2 (G k) 0 / 2) * I 2 1
        + τ ^ 4 * (deriv (G k) 0 / 2) * I 1 2 + τ ^ 4 * (G k 0 / 6) * I 0 3 with hS
  have ha : a k (h p) = a k hDiamond + 8 / π * X := by rw [hX, ← mul_assoc, ← h2]; ring
  have a1 := four_k_ne k; have b1 := four_k_ne' k
  have c1 : (16 * (k : ℝ) ^ 2 - 1) ≠ 0 := by
    rw [show (16 * (k : ℝ) ^ 2 - 1) = (4 * k - 1) * (4 * k + 1) by ring]; exact mul_ne_zero a1 b1
  have hI03 : I 0 3 = (3 / 4 * (π ^ 4 / 90) - 3 * I 1 2) / 2 := by linarith [I0312]
  have key : a k (h p) - (-1 : ℝ) ^ k * (Real.sqrt 2 / π) *
        (-(4 / (16 * (k : ℝ) ^ 2 - 1)) + (π ^ 2 / 6) * τ ^ 2 + z3 * τ ^ 3
          - (33 / 16 + 7 * (k : ℝ) ^ 2) * (π ^ 4 / 90) * τ ^ 4) = 8 / π * (X - S) := by
    have hd2 : iteratedDeriv 2 (G k) 0 / 2 = -((-1 : ℝ) ^ k * Real.sqrt 2 * ((k : ℝ) ^ 2 / 2 + 5 / 32)) :=
      G_d2 k
    rw [ha, h1, hS, hd2, G_zero k, G_d1 k, I01, I11, I02, I21, hI03]
    field_simp
    ring
  rw [key, abs_mul, abs_of_pos (by positivity)]
  calc 8 / π * |X - S| ≤ 8 / π * (CW * τ ^ 5) := mul_le_mul_of_nonneg_left hb (by positivity)
    _ = _ := by ring

end Y04

theorem Y04.proof_8 : Y04.stmt_8 := Y04.fourier_expansion

namespace Y04

lemma ell_hasDeriv (s : ℝ) : HasDerivAt ell (-Real.exp (-s) / (1 + Real.exp (-s))) s := by
  have h1 : HasDerivAt (fun s => 1 + Real.exp (-s)) (Real.exp (-s) * -1) s :=
    ((hasDerivAt_neg s).exp).const_add 1
  have h := h1.log (by positivity)
  refine h.congr_deriv ?_
  ring

lemma exp_ell (s : ℝ) : Real.exp (ell s) = 1 + Real.exp (-s) := by
  unfold ell; exact Real.exp_log (by positivity)

lemma Psi_int' {τ : ℝ} (hτ : τ < 1) :
    IntegrableOn (fun s : ℝ => Real.exp (τ * s / 2) * (Real.exp (τ * ell s) - 1)) (Set.Ioi 0) := by
  rcases le_or_gt τ 0 with h0 | h0
  · exact Psi_int h0
  refine intOn_of_bound (by have := continuous_ell; fun_prop)
    ((integrableOn_exp_mul_Ioi (a := τ / 2 - 1) (by linarith) 0).const_mul 3) (fun s hs => ?_)
  have hl0 := ell_nonneg s
  have hl1 := ell_le_one s hs.le
  have hle := ell_le s
  have hu0 : 0 ≤ τ * ell s := by positivity
  have hu1 : τ * ell s ≤ 1 := by have := mul_le_mul hτ.le hl1 hl0 zero_le_one; linarith
  have hA : Real.exp (τ * ell s) - 1 ≤ τ * ell s * Real.exp (τ * ell s) := by
    have := Real.add_one_le_exp (-(τ * ell s))
    have h2 : Real.exp (-(τ * ell s)) * Real.exp (τ * ell s) = 1 := by
      rw [← Real.exp_add]; simp
    nlinarith [Real.exp_pos (τ * ell s)]
  have hE : Real.exp (τ * ell s) ≤ 3 := by
    have := Real.exp_le_exp.2 hu1
    have := Real.exp_one_lt_d9; linarith
  have hA0 : 0 ≤ Real.exp (τ * ell s) - 1 := by linarith [Real.add_one_le_exp (τ * ell s)]
  rw [abs_mul, abs_of_pos (Real.exp_pos _), abs_of_nonneg hA0]
  have hB : Real.exp (τ * ell s) - 1 ≤ 3 * Real.exp (-s) := by
    have : τ * ell s ≤ Real.exp (-s) := by
      have := mul_le_mul hτ.le hle hl0 zero_le_one; linarith
    nlinarith [Real.exp_pos (τ * ell s), Real.exp_pos (-s)]
  calc Real.exp (τ * s / 2) * (Real.exp (τ * ell s) - 1) ≤ Real.exp (τ * s / 2) * (3 * Real.exp (-s)) :=
        mul_le_mul_of_nonneg_left hB (Real.exp_pos _).le
    _ = 3 * Real.exp ((τ / 2 - 1) * s) := by
        rw [show (τ / 2 - 1) * s = τ * s / 2 + -s by ring, Real.exp_add]; ring

/-- the primitive used to evaluate `Ψ` for every real `τ < 1` -/
noncomputable def Fprim (τ s : ℝ) : ℝ :=
  1 / 2 * Real.exp (τ * s / 2) * ((1 - Real.exp (-s)) * Real.exp ((τ - 1) * ell s) - 1)

lemma Fprim_hasDeriv (τ s : ℝ) :
    HasDerivAt (Fprim τ)
      (τ / 4 * (Real.exp (τ * s / 2) * (Real.exp (τ * ell s) - 1))
        - (τ - 1) * Real.exp (-2 * (1 - τ / 2) * (s / 2 + ell s))) s := by
  have hP : HasDerivAt (fun s => Real.exp (τ * s / 2)) (Real.exp (τ * s / 2) * (τ / 2)) s := by
    have := ((hasDerivAt_id s).const_mul τ).div_const 2 |>.exp
    refine this.congr_deriv ?_; simp
  have hE : HasDerivAt (fun s => 1 - Real.exp (-s)) (-(Real.exp (-s) * -1)) s :=
    ((hasDerivAt_neg s).exp).const_sub 1
  have hQ : HasDerivAt (fun s => Real.exp ((τ - 1) * ell s))
      (Real.exp ((τ - 1) * ell s) * ((τ - 1) * (-Real.exp (-s) / (1 + Real.exp (-s))))) s :=
    ((ell_hasDeriv s).const_mul (τ - 1)).exp
  have h : HasDerivAt (Fprim τ) _ s := ((hP.const_mul (1 / 2)).mul ((hE.mul hQ).sub_const 1))
  refine h.congr_deriv ?_
  have e1 : Real.exp (τ * ell s) = Real.exp ((τ - 1) * ell s) * (1 + Real.exp (-s)) := by
    rw [← exp_ell, ← Real.exp_add]; ring_nf
  have e2 : Real.exp (-2 * (1 - τ / 2) * (s / 2 + ell s))
      = Real.exp (τ * s / 2) * Real.exp (-s) * Real.exp ((τ - 1) * ell s) / (1 + Real.exp (-s)) := by
    rw [← exp_ell, ← Real.exp_add, ← Real.exp_add, ← Real.exp_sub]; ring_nf
  rw [e1, e2]
  have : 0 < 1 + Real.exp (-s) := by positivity
  field_simp
  ring

lemma Fprim_zero (τ : ℝ) : Fprim τ 0 = -1 / 2 := by simp [Fprim]; norm_num

lemma Fprim_tendsto {τ : ℝ} (hτ : τ < 1) : Tendsto (Fprim τ) atTop (𝓝 0) := by
  have hlim : Tendsto (fun s : ℝ => (2 - τ) * Real.exp ((τ / 2 - 1) * s)) atTop (𝓝 0) := by
    have : Tendsto (fun s : ℝ => (τ / 2 - 1) * s) atTop atBot :=
      tendsto_id.const_mul_atTop_of_neg (by linarith)
    simpa using (Real.tendsto_exp_atBot.comp this).const_mul (2 - τ)
  refine squeeze_zero_norm' ?_ hlim
  filter_upwards [Filter.eventually_ge_atTop 0] with s hs
  have hl0 := ell_nonneg s
  have hle := ell_le s
  have hE0 := Real.exp_pos (-s)
  have hE1 : Real.exp (-s) ≤ 1 := Real.exp_le_one_iff.2 (by linarith)
  rw [Real.norm_eq_abs, Fprim, abs_mul, abs_mul, abs_of_pos (Real.exp_pos _),
    abs_of_pos (by norm_num : (0:ℝ) < 1 / 2)]
  set Q := Real.exp ((τ - 1) * ell s)
  have hQ0 : 0 < Q := Real.exp_pos _
  have hQ1 : Q ≤ 1 := Real.exp_le_one_iff.2 (by nlinarith)
  have hQ2 : 1 - Q ≤ (1 - τ) * Real.exp (-s) := by
    have := Real.add_one_le_exp ((τ - 1) * ell s)
    have : (1 - τ) * ell s ≤ (1 - τ) * Real.exp (-s) := mul_le_mul_of_nonneg_left hle (by linarith)
    linarith
  have hb : |(1 - Real.exp (-s)) * Q - 1| ≤ (2 - τ) * Real.exp (-s) := by
    rw [abs_le]; constructor <;> nlinarith
  calc 1 / 2 * Real.exp (τ * s / 2) * |(1 - Real.exp (-s)) * Q - 1|
      ≤ 1 * Real.exp (τ * s / 2) * ((2 - τ) * Real.exp (-s)) := by
        gcongr; norm_num
    _ = (2 - τ) * Real.exp ((τ / 2 - 1) * s) := by
        rw [show (τ / 2 - 1) * s = τ * s / 2 + -s by ring, Real.exp_add]; ring

/-- the closed form of `Ψ` for every real `τ < 1`, `τ ≠ 0` -/
lemma Psi_closed {τ : ℝ} (hτ : τ < 1) (hτ0 : τ ≠ 0) :
    ∫ s in Set.Ioi (0 : ℝ), Real.exp (τ * s / 2) * (Real.exp (τ * ell s) - 1) =
      (R (-τ / 2) - 1) / (-τ / 2) := by
  have hy : 0 < 1 - τ / 2 := by linarith
  obtain ⟨hIc, hvc⟩ := cosh_integral hy
  have hIP := Psi_int' hτ
  have hint : IntegrableOn (fun s => τ / 4 * (Real.exp (τ * s / 2) * (Real.exp (τ * ell s) - 1))
      - (τ - 1) * Real.exp (-2 * (1 - τ / 2) * (s / 2 + ell s))) (Set.Ioi 0) :=
    (hIP.const_mul _).sub (hIc.const_mul _)
  have hFTC := integral_Ioi_of_hasDerivAt_of_tendsto (f := Fprim τ)
    ((Fprim_hasDeriv τ 0).continuousAt.continuousWithinAt)
    (fun s _ => Fprim_hasDeriv τ s) hint (Fprim_tendsto hτ)
  rw [integral_sub (hIP.const_mul _) (hIc.const_mul _), integral_const_mul, integral_const_mul,
    hvc, Fprim_zero] at hFTC
  set Ψ := ∫ s in Set.Ioi (0 : ℝ), Real.exp (τ * s / 2) * (Real.exp (τ * ell s) - 1)
  -- Gamma bookkeeping
  have hx1 : 0 < 1 - τ := by linarith
  have hG2 : Real.Gamma (2 * (1 - τ / 2)) = (1 - τ) * Real.Gamma (1 - τ) := by
    rw [show 2 * (1 - τ / 2) = (1 - τ) + 1 by ring, Real.Gamma_add_one hx1.ne']
  have hA : Real.Gamma (1 - τ / 2) = Real.Gamma (1 + -τ / 2) := by congr 1; ring
  have hB : Real.Gamma (1 + 2 * (-τ / 2)) = Real.Gamma (1 - τ) := by congr 1; ring
  have hGp : 0 < Real.Gamma (1 - τ) := Real.Gamma_pos_of_pos hx1
  rw [hG2, hA] at hFTC
  unfold R
  rw [hB]
  generalize Real.Gamma (1 + -τ / 2) = A at hFTC ⊢
  generalize Real.Gamma (1 - τ) = B at hFTC hGp ⊢
  have hc : (τ - 1) * (A * A / ((1 - τ) * B) / 2) = -(A ^ 2 / B / 2) := by
    field_simp; ring
  rw [hc] at hFTC
  have key : τ * Ψ = 2 - 2 * (A ^ 2 / B) := by linarith
  have ht2 : -τ / 2 ≠ 0 := by intro h; apply hτ0; linarith
  rw [eq_div_iff ht2]
  linarith

end Y04

theorem Y04.proof_6 : Y04.stmt_6 := fun _ hτ hτ0 => ⟨Y04.Psi_int' hτ, Y04.Psi_closed hτ hτ0⟩
