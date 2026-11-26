import Mathlib
import Zeta3.Wave1.T13_Lambda
import Zeta3.Wave2.D00_Dilog
/-!
# Y02 -- exact length and asymptotics of Ahmed's power-law corners at every angle

Target id: Y02.
Paper labels: `thm:P2` and `thm:P2asym` (paper/ieee/main.tex, lines 1423-1460; setting and the
definition of `Λ_ν` at lines 1408-1421; the uniformity paragraph after the theorem at line 1462).
Detailed proofs: paper/ieee/proofs/sec_g3_P2.tex (Lemmas `g3:lem:Lam`, `g3:lem:Mellin`,
`g3:lem:defect`, `g3:lem:gegen`, Steps 1-5 of the proof of the length formula) and
paper/ieee/proofs/sec_g3_P2asym.tex (Lemmas `g3:lem:sums`, `g3:lem:ranges`, `g3:lem:taylor`,
`g3:prop:unif`, `g3:prop:n3`).

Quoted verbatim (main.tex, setting):

  His curve is $\gamma_n(x)=x\,u+x^n\,v$, $0\le x\le1$, where $u,v$ are unit vectors with
  $u\cdot v=c=\cos\theta$, $0\le\theta<\pi$.
  With $\varepsilon=1/(n-1)$, his starting point is
  \begin{gather*}
  2-L_\theta(n)=\varepsilon n^{-\varepsilon}\int_0^nh_c(s)\,s^{\varepsilon-1}ds,\\
  h_c(s)=1+s-\sqrt{1+2cs+s^2}.
  \end{gather*}
  \[
  \Lambda_\nu(\theta)=\int_{\cos\theta}^1P_\nu(t)\dd t=2\sin^2\tfrac\theta2\;{}_2F_1\bigl(-\nu,1+\nu;2;\sin^2\tfrac\theta2\bigr).
  \]

Quoted verbatim (main.tex, thm:P2):

  For $0<\varepsilon<1$,
  \[
  \int_0^\infty\Bigl(\sqrt{1+2cs+s^2}-s-c\Bigr)s^{\varepsilon-1}ds=\frac{\pi}{\sin\pi\varepsilon}\,\Lambda_\varepsilon(\theta),
  \]
  and for every $n>1$ (with removable singularities at $n=1+\frac1j$),
  \begin{multline*}
  L_\theta(n)=1+\cos\theta+\frac{\pi\varepsilon}{\sin\pi\varepsilon}\,\Lambda_\varepsilon(\theta)\,n^{-\varepsilon}\\
  +\sum_{m=1}^{\infty}(-1)^m\frac{\varepsilon}{m-\varepsilon}\,\Lambda_m(\theta)\,n^{-m},
  \end{multline*}
  where $\Lambda_m(\theta)=\frac{P_{m-1}(\cos\theta)-P_{m+1}(\cos\theta)}{2m+1}$ is a polynomial in
  $\cos\theta$. At $\theta=\pi/2$ this is Ahmed's exact formula.

Quoted verbatim (main.tex, before and in thm:P2asym):

  Expanding $\log\bigl(\frac{\pi\varepsilon}{\sin\pi\varepsilon}\Lambda_\varepsilon\bigr)=\log(2z)+\sum_ja_j\varepsilon^j$
  with $z=\sin^2\frac\theta2$, $w=\cot^2\frac\theta2$, $\ell=\log\cos^2\frac\theta2$ and $U_\theta=w\ell$
  gives the asymptotics.

  For $0<\theta<\pi$, with $\mu_\theta=e^{a_1}$,
  \begin{multline*}
  L_\theta(n)=1+\cos\theta+2\sin^2\tfrac\theta2\Bigl(\frac{\mu_\theta}{n}\Bigr)^{\frac1{n-1}}\\
  +\frac{c_2(\theta)}{n^2}+O\Bigl(\frac{\log n}{n^3}\Bigr),\qquad c_2(\theta)=2z(a_2-1+z),
  \end{multline*}
  where $a_1=-1-U_\theta$, $a_2=\zeta(2)+\frac12-\Li_2(z)-\frac12U_\theta^2$ and
  $a_3=2w\,S_{1,2}(z)-\frac13(1+U_\theta^3)$, with the Nielsen polylogarithm
  \[
  S_{1,2}(z)=\zeta(3)-\Li_3(1-z)+\ell\,\Li_2(1-z)+\tfrac12\ell^2\log z .
  \]

Quoted verbatim (main.tex line 1462, after the theorem):

  The remainder is uniform in the angle: for $n\ge9$ and every $0<\theta<\pi$ it is at most
  $C\sin^2\frac\theta2\,\log n/n^3$ with an absolute constant $C$, and its leading part is
  $-2z\,a_2\log n/n^3$.

## Conventions

* `L θ n` is the arc length the paper (and `g3:lem:defect`) gives:
  `L_θ(n) = ∫_0^1 √(1 + 2c y + y²) dx`, `y = n x^{n-1}`, which is `∫_0^1 |γ_n'(x)| dx` because
  `|u + y v|² = 1 + 2cy + y²` for unit `u, v` with `u·v = c`. For `n > 1` the integrand is
  continuous on `[0,1]` (`x^{n-1}` is `rpow` with positive exponent), so the interval integral is
  a genuine integral. Only `n > 1` is used.
* `Λ_ν(θ)` for real `ν` is the paper's hypergeometric expression `2z ₂F₁(-ν,1+ν;2;z)`,
  `z = sin²(θ/2)`, written as the `tsum` of `T13.lamTerm` (Mathlib has no Legendre function
  `P_ν` for non-integer `ν`). stmt_1 asserts the series is summable for `0 ≤ θ < π`, so the
  `tsum` is never a junk value where it is used.
* `Λ_m(θ)` for integer `m ≥ 1` is `T13.LamInt m θ = ∫_{cos θ}^1 P_m(t) dt` (Legendre polynomial);
  T13 proves it equals `(P_{m-1}(c) - P_{m+1}(c))/(2m+1)` and the hypergeometric series.
* Range of θ in thm:P2: the main text states `0 ≤ θ < π` in the setting (line 1408); we use it.
  thm:P2asym uses the paper's `0 < θ < π`.
* `ζ(2)` is written `π^2/6` (Mathlib `riemannZeta_two`). `ζ(3)` is written `z3 := (riemannZeta 3).re`,
  the real part of Mathlib's `riemannZeta 3`, which is real and equals `∑_{n≥1} n^{-3}`
  (`zeta_nat_eq_tsum_of_gt_one`).
* `Li_2 = D00.Li2` (`-∫_0^x log(1-t)/t dt`, the real dilogarithm). `Li_3` is defined here by
  `Li_3(x) = ∫_0^x Li_2(t)/t dt` (the standard integral definition; the integrand is continuous on
  `(0,1]` and bounded near `0`; the value `0/0 = 0` at the single point `t = 0` is irrelevant).
* Improper integrals over `(0,∞)` and `(0,n]` are Lebesgue set integrals, and each statement
  asserts integrability, so none holds by the "non-integrable ⇒ 0" convention.
* `ε = 1/(n-1)`. At `ε ∈ ℕ` (i.e. `n = 1 + 1/j`) the right side of the length formula has
  `sin(πε) = 0` and a zero denominator `m - ε`; stmt_6 excludes these points (where Lean's `x/0 = 0`
  would make the formula meaningless), and stmt_8 states the removability as a limit.
* The `O(log n/n³)` of thm:P2asym is stated as an explicit bound `∃ C N, ∀ n ≥ N` for real `n`
  (stmt_14), and the uniform version of line 1462 as `∃ C, ∀ θ, ∀ n ≥ 9` (stmt_13).
* The `a_j` are defined in the paper as Taylor coefficients of `log(R_θ(ε)/2z)`,
  `R_θ(ε) = (πε/sin πε) Λ_ε(θ)`. We state the paper's closed forms `a_1, a_2, a_3` and assert
  (stmt_11, stmt_12, from `g3:lem:taylor` and the proof of `g3:prop:n3`) that with these values
  `R_θ(ε) e^{-a_1 ε} = 2z(1 + a_2 ε² + a_3 ε³ + O(ε⁴))` with the supplement's explicit constants;
  this pins them down as the Taylor coefficients of `log(R_θ/2z)` up to order 3.

## Statements

* stmt_1: the series `Λ_ν(θ) = Σ_k lamTerm ν z k` converges (0 ≤ θ < π, every real ν).
* stmt_2: `|Λ_m(θ)| ≤ min(1 - cos θ, 2/(2m+1))` for `m ≥ 1` (g3:lem:Lam(c)).
* stmt_3: the defect formula `2 - L_θ(n) = ε n^{-ε} ∫_0^n h_c(s) s^{ε-1} ds` (g3:lem:defect).
* stmt_4: expansion at infinity `k_c(s) = Σ_{m≥2} (-1)^m Λ_{m-1} s^{1-m}`, `s > 1` (g3:lem:gegen).
* stmt_5: the Mellin identity of thm:P2, with integrability.
* stmt_6: the length formula of thm:P2 as a `HasSum`, for `n > 1`, `ε ∉ ℕ`.
* stmt_7: `n ↦ L_θ(n)` is continuous on `(1,∞)` (Step 5).
* stmt_8: removable singularities: the right side tends to `L_θ(n₀)` at `n₀ = 1 + 1/j`.
* stmt_9: ranges `-1 < U_θ < 0`, `-1 < a_1 < 0`, `0 < a_2 < ζ(2) + 1/2` (g3:lem:ranges).
* stmt_10: `S_{1,2}(z) = ½ ∫_0^z log²(1-t)/t dt` (g3:lem:sums(iv)), the ζ(3) identity.
* stmt_11: Taylor with uniform remainder: `|R e^{-a_1ε} - 2z(1 + a_2ε²)| ≤ 2z·3100 ε³`,
  `0 < ε ≤ 1/8` (g3:lem:taylor).
* stmt_12: one order further with `a_3`: `|R e^{-a_1ε} - 2z(1 + a_2ε² + a_3ε³)| ≤ 2z·31000 ε⁴`
  (proof of g3:prop:n3).
* stmt_13: uniform remainder: `∃ C, ∀ 0<θ<π, ∀ n ≥ 9, |E_n| ≤ C sin²(θ/2) log n / n³`
  (main.tex line 1462, g3:prop:unif first claim).
* stmt_14: thm:P2asym: for each `0 < θ < π`, `|E_n| ≤ C log n / n³` for `n ≥ N`.

## NOT COVERED

* The Legendre function `P_ν` for non-integer `ν` and the equality
  `∫_{cos θ}^1 P_ν = 2z ₂F₁(-ν,1+ν;2;z)` (no `P_ν` in Mathlib; the hypergeometric form is used as
  the definition). The ODE/uniqueness argument (g3:lem:Mellin Steps 2-5) is only a proof method.
* "At θ = π/2 this is Ahmed's exact formula": covered in T13 (Wave1), not repeated.
* `Λ_m` is a polynomial in `cos θ`: covered in T13 (stmt_2, stmt_3), not repeated.
* The second half of line 1462 (`E_n = -2z a_2 log n/n³ + Δ_n`, `|Δ_n| ≤ 9000 z/n³`) and the
  terms of order `n^{-3}`, `n^{-4}` and `a_4` ("explicit too (in our notes)"); `a_3` enters only
  through stmt_12, not through an expansion of `L_θ(n)`.
* Meromorphic continuation in complex `n` (only the real-variable statement of removability).
-/
-- BEGIN STATEMENTS
namespace Y02

open Real MeasureTheory Filter Topology

/-- `ζ(3)` as a real number: the real part of Mathlib's `riemannZeta 3` (which is real). -/
noncomputable def z3 : ℝ := (riemannZeta 3).re

/-- The arc length `L_θ(n) = ∫_0^1 √(1 + 2 cos θ · n x^{n-1} + (n x^{n-1})²) dx` of
`γ_n(x) = x u + x^n v`, `u·v = cos θ` (g3:lem:defect, first line of its proof). -/
noncomputable def L (θ n : ℝ) : ℝ :=
  ∫ x in (0 : ℝ)..1,
    Real.sqrt (1 + 2 * Real.cos θ * (n * x ^ (n - 1)) + (n * x ^ (n - 1)) ^ 2)

/-- `h_c(s) = 1 + s - √(1 + 2cs + s²)`. -/
noncomputable def h (c s : ℝ) : ℝ := 1 + s - Real.sqrt (1 + 2 * c * s + s ^ 2)

/-- `k_c(s) = √(1 + 2cs + s²) - s - c` (the integrand of the first identity). -/
noncomputable def k (c s : ℝ) : ℝ := Real.sqrt (1 + 2 * c * s + s ^ 2) - s - c

/-- `Λ_ν(θ) = 2z ₂F₁(-ν,1+ν;2;z)`, `z = sin²(θ/2)`, as the sum of the series `T13.lamTerm`
(summable by stmt_1). -/
noncomputable def Lam (ν θ : ℝ) : ℝ := ∑' j : ℕ, T13.lamTerm ν (Real.sin (θ / 2) ^ 2) j

/-- `ε = 1/(n-1)`. -/
noncomputable def eps (n : ℝ) : ℝ := 1 / (n - 1)

/-- The `m`-th term (`m = i+1`) of the series in thm:P2:
`(-1)^m ε/(m-ε) Λ_m(θ) n^{-m}`. -/
noncomputable def serTerm (θ n : ℝ) (i : ℕ) : ℝ :=
  (-1 : ℝ) ^ (i + 1) * (eps n / ((i + 1 : ℝ) - eps n)) * T13.LamInt (i + 1) θ / n ^ (i + 1)

/-- The right side of the length formula of thm:P2. -/
noncomputable def RHS (θ n : ℝ) : ℝ :=
  1 + Real.cos θ + π * eps n / Real.sin (π * eps n) * Lam (eps n) θ * n ^ (-eps n)
    + ∑' i : ℕ, serTerm θ n i

/-- `Li_3(x) = ∫_0^x Li_2(t)/t dt`. -/
noncomputable def Li3 (x : ℝ) : ℝ := ∫ t in (0 : ℝ)..x, D00.Li2 t / t

/-- `z = sin²(θ/2)`. -/
noncomputable def zz (θ : ℝ) : ℝ := Real.sin (θ / 2) ^ 2

/-- `w = cot²(θ/2) = (1-z)/z`. -/
noncomputable def ww (θ : ℝ) : ℝ := (1 - zz θ) / zz θ

/-- `ℓ = log cos²(θ/2) = log(1-z)`. -/
noncomputable def ell (θ : ℝ) : ℝ := Real.log (1 - zz θ)

/-- `U_θ = w ℓ`. -/
noncomputable def U (θ : ℝ) : ℝ := ww θ * ell θ

/-- `a_1 = -1 - U_θ`. -/
noncomputable def a1 (θ : ℝ) : ℝ := -1 - U θ

/-- `a_2 = ζ(2) + 1/2 - Li_2(z) - U_θ²/2`, with `ζ(2) = π²/6`. -/
noncomputable def a2 (θ : ℝ) : ℝ := π ^ 2 / 6 + 1 / 2 - D00.Li2 (zz θ) - U θ ^ 2 / 2

/-- Nielsen `S_{1,2}(z) = ζ(3) - Li_3(1-z) + ℓ Li_2(1-z) + ½ ℓ² log z`. -/
noncomputable def S12 (θ : ℝ) : ℝ :=
  z3 - Li3 (1 - zz θ) + ell θ * D00.Li2 (1 - zz θ) + ell θ ^ 2 * Real.log (zz θ) / 2

/-- `a_3 = 2w S_{1,2}(z) - (1 + U_θ³)/3`. -/
noncomputable def a3 (θ : ℝ) : ℝ := 2 * ww θ * S12 θ - (1 + U θ ^ 3) / 3

/-- `μ_θ = e^{a_1}`. -/
noncomputable def mu (θ : ℝ) : ℝ := Real.exp (a1 θ)

/-- `c_2(θ) = 2z(a_2 - 1 + z)`. -/
noncomputable def c2 (θ : ℝ) : ℝ := 2 * zz θ * (a2 θ - 1 + zz θ)

/-- `R_θ(ε) = (πε/sin πε) Λ_ε(θ)`. -/
noncomputable def R (θ ε : ℝ) : ℝ := π * ε / Real.sin (π * ε) * Lam ε θ

/-- The remainder `E_n = L_θ(n) - [1 + cos θ + 2z(μ_θ/n)^{1/(n-1)} + c_2(θ)/n²]`. -/
noncomputable def E (θ n : ℝ) : ℝ :=
  L θ n - (1 + Real.cos θ + 2 * zz θ * (mu θ / n) ^ (1 / (n - 1)) + c2 θ / n ^ 2)

/-- g3:lem:Lam(a): the hypergeometric series defining `Λ_ν(θ)` converges, `0 ≤ θ < π`. -/
def stmt_1 : Prop :=
  ∀ ν θ : ℝ, 0 ≤ θ → θ < π → Summable (fun j : ℕ => T13.lamTerm ν (Real.sin (θ / 2) ^ 2) j)

/-- g3:lem:Lam(c): `|Λ_m(θ)| ≤ min(1-c, 2/(2m+1))` for integers `m ≥ 1`. -/
def stmt_2 : Prop :=
  ∀ m : ℕ, 1 ≤ m → ∀ θ : ℝ, 0 ≤ θ → θ < π →
    |T13.LamInt m θ| ≤ min (1 - Real.cos θ) (2 / (2 * (m : ℝ) + 1))

/-- main.tex 1414 / g3:lem:defect: `2 - L_θ(n) = ε n^{-ε} ∫_0^n h_c(s) s^{ε-1} ds`, `n > 1`
(with integrability of the integrand on `(0,n]`). -/
def stmt_3 : Prop :=
  ∀ θ n : ℝ, 0 ≤ θ → θ < π → 1 < n →
    IntegrableOn (fun s => h (Real.cos θ) s * s ^ (eps n - 1)) (Set.Ioc 0 n) ∧
    2 - L θ n = eps n * n ^ (-eps n) *
      ∫ s in Set.Ioc 0 n, h (Real.cos θ) s * s ^ (eps n - 1)

/-- g3:lem:gegen: for `s > 1`, `k_c(s) = Σ_{m≥2} (-1)^m Λ_{m-1}(θ) s^{1-m}`
(index `i = m - 2`; `(-1)^m = (-1)^i`). -/
def stmt_4 : Prop :=
  ∀ θ s : ℝ, 0 ≤ θ → θ < π → 1 < s →
    HasSum (fun i : ℕ => (-1 : ℝ) ^ i * T13.LamInt (i + 1) θ / s ^ (i + 1)) (k (Real.cos θ) s)

/-- thm:P2, first identity: for `0 < ε < 1`,
`∫_0^∞ (√(1+2cs+s²) - s - c) s^{ε-1} ds = (π/sin πε) Λ_ε(θ)`, the integral converging. -/
def stmt_5 : Prop :=
  ∀ θ ε : ℝ, 0 ≤ θ → θ < π → 0 < ε → ε < 1 →
    IntegrableOn (fun s => k (Real.cos θ) s * s ^ (ε - 1)) (Set.Ioi 0) ∧
    ∫ s in Set.Ioi 0, k (Real.cos θ) s * s ^ (ε - 1) = π / Real.sin (π * ε) * Lam ε θ

/-- thm:P2, length formula, for `n > 1` with `ε = 1/(n-1) ∉ ℕ`:
`L_θ(n) = 1 + cos θ + (πε/sin πε) Λ_ε(θ) n^{-ε} + Σ_{m≥1} (-1)^m ε/(m-ε) Λ_m(θ) n^{-m}`. -/
def stmt_6 : Prop :=
  ∀ θ n : ℝ, 0 ≤ θ → θ < π → 1 < n → (∀ j : ℕ, eps n ≠ j) →
    HasSum (serTerm θ n)
      (L θ n - (1 + Real.cos θ + π * eps n / Real.sin (π * eps n) * Lam (eps n) θ * n ^ (-eps n)))

/-- Step 5 of the proof of thm:P2: `L_θ(n)` is continuous in `n > 1`. -/
def stmt_7 : Prop :=
  ∀ θ : ℝ, 0 ≤ θ → θ < π → ContinuousOn (L θ) (Set.Ioi 1)

/-- thm:P2, "with removable singularities at $n=1+\frac1j$": at `n₀ = 1 + 1/j` the right side of
the length formula tends to `L_θ(n₀)`. -/
def stmt_8 : Prop :=
  ∀ θ : ℝ, 0 ≤ θ → θ < π → ∀ j : ℕ, 1 ≤ j →
    Tendsto (RHS θ) (𝓝[≠] (1 + 1 / (j : ℝ))) (𝓝 (L θ (1 + 1 / (j : ℝ))))

/-- g3:lem:ranges: `-1 < U_θ < 0`, `-1 < a_1 < 0` (so `μ_θ ∈ (e^{-1},1)`), and
`0 < a_2 < ζ(2) + 1/2`, for `0 < θ < π`. -/
def stmt_9 : Prop :=
  ∀ θ : ℝ, 0 < θ → θ < π →
    -1 < U θ ∧ U θ < 0 ∧ -1 < a1 θ ∧ a1 θ < 0 ∧ 0 < a2 θ ∧ a2 θ < π ^ 2 / 6 + 1 / 2

/-- g3:lem:sums(iv), last claim: the closed form
`S_{1,2}(z) = ζ(3) - Li_3(1-z) + ℓ Li_2(1-z) + ½ℓ² log z` equals `½ ∫_0^z log²(1-t)/t dt`. -/
def stmt_10 : Prop :=
  ∀ θ : ℝ, 0 < θ → θ < π →
    S12 θ = (1 / 2) * ∫ t in (0 : ℝ)..zz θ, Real.log (1 - t) ^ 2 / t

/-- g3:lem:taylor (the `a_1, a_2` of the theorem are the Taylor coefficients of
`log(R_θ/2z)`): `R_θ(ε) e^{-a_1 ε} = 2z(1 + a_2 ε² + ρ)`, `|ρ| ≤ 3100 ε³` for `0 < ε ≤ 1/8`. -/
def stmt_11 : Prop :=
  ∀ θ ε : ℝ, 0 < θ → θ < π → 0 < ε → ε ≤ 1 / 8 →
    |R θ ε * Real.exp (-(a1 θ * ε)) - 2 * zz θ * (1 + a2 θ * ε ^ 2)| ≤ 2 * zz θ * (3100 * ε ^ 3)

/-- Proof of g3:prop:n3 (the `a_3` of the theorem is the next Taylor coefficient):
`R_θ(ε) e^{-a_1 ε} = 2z(1 + a_2 ε² + a_3 ε³ + ρ₄)`, `|ρ₄| ≤ 31000 ε⁴` for `0 < ε ≤ 1/8`. -/
def stmt_12 : Prop :=
  ∀ θ ε : ℝ, 0 < θ → θ < π → 0 < ε → ε ≤ 1 / 8 →
    |R θ ε * Real.exp (-(a1 θ * ε)) - 2 * zz θ * (1 + a2 θ * ε ^ 2 + a3 θ * ε ^ 3)|
      ≤ 2 * zz θ * (31000 * ε ^ 4)

/-- main.tex line 1462 / g3:prop:unif: the remainder is uniform in the angle,
`|E_n| ≤ C sin²(θ/2) log n / n³` for `n ≥ 9`, `0 < θ < π`, with an absolute constant `C`. -/
def stmt_13 : Prop :=
  ∃ C : ℝ, ∀ θ n : ℝ, 0 < θ → θ < π → 9 ≤ n →
    |E θ n| ≤ C * Real.sin (θ / 2) ^ 2 * Real.log n / n ^ 3

/-- thm:P2asym: for `0 < θ < π`,
`L_θ(n) = 1 + cos θ + 2 sin²(θ/2) (μ_θ/n)^{1/(n-1)} + c_2(θ)/n² + O(log n / n³)`. -/
def stmt_14 : Prop :=
  ∀ θ : ℝ, 0 < θ → θ < π → ∃ C N : ℝ, ∀ n : ℝ, N ≤ n →
    |L θ n - (1 + Real.cos θ + 2 * Real.sin (θ / 2) ^ 2 * (mu θ / n) ^ (1 / (n - 1))
        + c2 θ / n ^ 2)| ≤ C * Real.log n / n ^ 3

end Y02
-- END STATEMENTS

namespace Y02
open Real Filter Topology Polynomial MeasureTheory


lemma lamTerm_succ (ν z : ℝ) (k : ℕ) :
    T13.lamTerm ν z (k + 1) = T13.lamTerm ν z k *
      ((-ν + k) * (1 + ν + k) * z / (((k : ℝ) + 1) * (2 + k))) := by
  unfold T13.lamTerm
  rw [T13.pe_succ_right', T13.pe_succ_right', T13.pe_succ_right', Nat.factorial_succ]
  have h2 : (ascPochhammer ℝ k).eval 2 ≠ 0 := by
    rw [T13.pe_two]; positivity
  have hf : (k.factorial : ℝ) ≠ 0 := by positivity
  push_cast
  field_simp
  ring

lemma lamTerm_summable (ν z : ℝ) (hz0 : 0 ≤ z) (hz1 : z < 1) :
    Summable (fun j : ℕ => T13.lamTerm ν z j) := by
  set a := |ν|
  set d := (1 - z) / 2 with hd
  have hd0 : 0 < d := by rw [hd]; linarith
  refine summable_of_ratio_norm_eventually_le (r := (1 + z) / 2) (by linarith) ?_
  obtain ⟨N, hN⟩ := exists_nat_gt ((2 * a + a ^ 2 + 1) / d)
  rw [eventually_atTop]
  refine ⟨N, fun k hk => ?_⟩
  rw [lamTerm_succ, norm_mul]
  rw [mul_comm]
  apply mul_le_mul_of_nonneg_right _ (norm_nonneg _)
  have hk' : ((2 * a + a ^ 2 + 1) / d) < (k : ℝ) := lt_of_lt_of_le hN (by exact_mod_cast hk)
  have hx : (2 * a + a ^ 2 + 1) < d * ((k:ℝ) + 1) := by
    rw [div_lt_iff₀ hd0] at hk'; nlinarith
  have ha : 0 ≤ a := abs_nonneg ν
  have hpos : (0:ℝ) < ((k : ℝ) + 1) * (2 + k) := by positivity
  rw [Real.norm_eq_abs, abs_div, abs_of_pos hpos, div_le_iff₀ hpos, abs_mul, abs_mul,
    abs_of_nonneg hz0]
  have e1 : |-ν + k| ≤ k + a := by
    have := abs_add_le (-ν) (k:ℝ)
    rw [abs_neg, abs_of_nonneg (Nat.cast_nonneg (α := ℝ) k)] at this
    linarith
  have e2 : |1 + ν + k| ≤ k + 1 + a := by
    have := abs_add_le (ν) ((k:ℝ) + 1)
    rw [abs_of_nonneg (by positivity : (0:ℝ) ≤ k + 1)] at this
    have : 1 + ν + k = ν + (k + 1) := by ring
    rw [this]; linarith
  have hprod : |-ν + k| * |1 + ν + k| ≤ ((k:ℝ) + 1 + a) ^ 2 := by
    calc |-ν + k| * |1 + ν + k| ≤ (k + a) * (k + 1 + a) :=
          mul_le_mul e1 e2 (abs_nonneg _) (by positivity)
      _ ≤ _ := by nlinarith
  have hk0 : (0:ℝ) ≤ k := Nat.cast_nonneg k
  have key : z * ((k:ℝ) + 1 + a) ^ 2 ≤ (1 + z) / 2 * (((k : ℝ) + 1) * (2 + k)) := by
    have h3 : z * ((k:ℝ) + 1 + a) ^ 2 ≤ z * ((k:ℝ)+1)^2 + (2 * a * ((k:ℝ)+1) + a^2) := by
      have : z * (2 * a * ((k:ℝ)+1) + a^2) ≤ 2 * a * ((k:ℝ)+1) + a^2 := by
        have : 0 ≤ 2 * a * ((k:ℝ)+1) + a^2 := by positivity
        nlinarith
      nlinarith
    have h4 : (2 * a * ((k:ℝ)+1) + a^2) ≤ d * ((k:ℝ)+1)^2 := by
      have hk1 : (1:ℝ) ≤ (k:ℝ) + 1 := by linarith
      nlinarith
    rw [hd] at h4
    nlinarith
  calc |-ν + k| * |1 + ν + k| * z ≤ ((k:ℝ) + 1 + a) ^ 2 * z :=
        mul_le_mul_of_nonneg_right hprod hz0
    _ ≤ _ := by linarith


lemma zz_bounds (θ : ℝ) (h0 : 0 < θ) (h1 : θ < π) :
    0 < Real.sin (θ / 2) ^ 2 ∧ Real.sin (θ / 2) ^ 2 < 1 := by
  have hs : 0 < Real.sin (θ / 2) := Real.sin_pos_of_pos_of_lt_pi (by linarith) (by linarith)
  have hc : 0 < Real.cos (θ / 2) := Real.cos_pos_of_mem_Ioo ⟨by linarith, by linarith⟩
  refine ⟨by positivity, ?_⟩
  nlinarith [Real.sin_sq_add_cos_sq (θ / 2)]

lemma Li2_pos_lt {z : ℝ} (h0 : 0 < z) (h1 : z < 1) : 0 < D00.Li2 z ∧ D00.Li2 z < π ^ 2 / 6 := by
  have hz := D00.proof_6 z (by rw [abs_le]; constructor <;> linarith)
  have h1' := D00.proof_6 1 (by norm_num)
  rw [D00.proof_7] at h1'
  constructor
  · refine hasSum_lt (f := fun _ => (0:ℝ)) (i := 0) (fun n => by positivity) (by simp; linarith)
      hasSum_zero hz
  · refine hasSum_lt (i := 0) (fun n => ?_) (by simp; linarith) hz h1'
    apply div_le_div_of_nonneg_right _ (by positivity)
    rw [one_pow]; exact pow_le_one₀ h0.le h1.le

lemma U_bounds {z : ℝ} (h0 : 0 < z) (h1 : z < 1) :
    -1 < (1 - z) / z * Real.log (1 - z) ∧ (1 - z) / z * Real.log (1 - z) < 0 := by
  have hy : 0 < 1 - z := by linarith
  have hl : Real.log (1 - z) < 0 := Real.log_neg hy (by linarith)
  constructor
  · have h := Real.log_lt_sub_one_of_pos (x := 1 / (1 - z)) (by positivity)
      (by rw [ne_eq, div_eq_one_iff_eq hy.ne']; linarith)
    rw [one_div, Real.log_inv] at h
    rw [div_mul_eq_mul_div, lt_div_iff₀ h0]
    have : (1 - z)⁻¹ * (1 - z) = 1 := inv_mul_cancel₀ hy.ne'
    have h2 := mul_lt_mul_of_pos_right h hy
    nlinarith
  · exact mul_neg_of_pos_of_neg (by positivity) hl



lemma bc_ratio (n k : ℕ) :
    T13.bc n (k + 1) * ((k : ℝ) + 1) ^ 2 = T13.bc n k * ((-(n:ℝ) + k) * (1 + n + k)) := by
  unfold T13.bc
  rw [T13.pe_succ_right', T13.pe_succ_right', Nat.factorial_succ]
  have hf : (k.factorial : ℝ) ≠ 0 := by positivity
  push_cast
  field_simp

/-- Hypergeometric ODE for `B n`: `X(1-X) B'' + (1-2X) B' + n(n+1) B = 0`. -/
lemma B_ode (n : ℕ) :
    X * derivative (derivative (T13.B n)) - X * (X * derivative (derivative (T13.B n)))
      + derivative (T13.B n) - C 2 * (X * derivative (T13.B n))
      + C ((n:ℝ) * (n + 1)) * T13.B n = 0 := by
  ext k
  simp only [coeff_add, coeff_sub, coeff_C_mul, coeff_zero]
  rcases k with _ | j
  · simp [coeff_derivative, T13.coeff_B, T13.bc_zero_idx, T13.bc_one_idx]
  · rw [coeff_X_mul, coeff_X_mul, coeff_X_mul]
    simp only [coeff_derivative, T13.coeff_B]
    have hr := bc_ratio n (j + 1)
    rcases j with _ | i
    · simp [T13.bc_one_idx] at hr ⊢
      linear_combination hr
    · rw [coeff_X_mul, coeff_derivative, coeff_derivative, T13.coeff_B]
      push_cast at hr ⊢
      linear_combination hr

lemma B_eval_zero (n : ℕ) : (T13.B n).eval 0 = 1 := by
  rw [T13.eval_B, Finset.sum_range_succ']
  simp [T13.bc_zero_idx]

lemma B_eq_map (n : ℕ) : T13.B n = (shiftedLegendre n).map (Int.castRingHom ℝ) := by
  ext k
  rw [T13.coeff_B, coeff_map, coeff_shiftedLegendre, T13.bc_eq]
  simp

lemma B_eval_aeval (n : ℕ) (x : ℝ) : (T13.B n).eval x = aeval x (shiftedLegendre n) := by
  rw [B_eq_map, ← algebraMap_int_eq, eval_map_algebraMap]

lemma B_eval_one (n : ℕ) : (T13.B n).eval 1 = (-1) ^ n := by
  have h := shiftedLegendre_eval_symm n (R := ℝ) 0
  have h0 := B_eval_zero n
  rw [B_eval_aeval] at h0
  rw [B_eval_aeval]
  rw [sub_zero, h0] at h
  have hsq : ((-1:ℝ) ^ n) * ((-1) ^ n) = 1 := by rw [← pow_add, ← two_mul, pow_mul]; simp
  calc aeval (1:ℝ) (shiftedLegendre n) = ((-1:ℝ) ^ n) * ((-1) ^ n) * aeval (1:ℝ) (shiftedLegendre n) := by
        rw [hsq, one_mul]
    _ = (-1) ^ n := by rw [mul_assoc, ← h, mul_one]

/-- `|B_n(u)| ≤ 1` on `[0,1]`. -/
lemma B_abs_le (n : ℕ) (u : ℝ) (hu0 : 0 ≤ u) (hu1 : u ≤ 1) : |(T13.B n).eval u| ≤ 1 := by
  rcases Nat.eq_zero_or_pos n with rfl | hn
  · simp [T13.B, T13.bc_zero_idx]
  set N : ℝ := (n:ℝ) * (n + 1) with hN
  have hNpos : 0 < N := by positivity
  set B := T13.B n
  set B1 := derivative B
  set B2 := derivative B1
  let g : ℝ → ℝ := fun v => B.eval v ^ 2 + v * (1 - v) * B1.eval v ^ 2 / N
  have hode : ∀ v : ℝ, v * B2.eval v - v * (v * B2.eval v) + B1.eval v - 2 * (v * B1.eval v)
      + N * B.eval v = 0 := by
    intro v
    have := congrArg (eval v) (B_ode n)
    simpa using this
  have hg : ∀ v, HasDerivAt g (-(1 - 2 * v) * B1.eval v ^ 2 / N) v := by
    intro v
    have h0 := B.hasDerivAt v
    have h1 := B1.hasDerivAt v
    have hv := (hasDerivAt_id' v).mul ((hasDerivAt_id' v).const_sub 1)
    have := ((h0.pow 2).add ((hv.mul (h1.pow 2)).div_const N))
    convert this using 1
    simp only [Pi.pow_apply, Pi.mul_apply]
    have := hode v
    field_simp
    simp only [B2] at this
    push_cast
    linear_combination (-2 * B1.eval v) * this
  have hgc : Continuous g := by
    show Continuous fun v => B.eval v ^ 2 + v * (1 - v) * B1.eval v ^ 2 / N
    fun_prop
  have hg0 : g 0 = 1 := by simp [g, B, B_eval_zero]
  have hg1 : g 1 = 1 := by
    simp only [g, B, B_eval_one]; simp [← pow_mul, pow_mul']
  have hle : g u ≤ 1 := by
    rcases le_total u (1/2) with h | h
    · have hA : AntitoneOn g (Set.Icc 0 (1/2)) := by
        apply antitoneOn_of_hasDerivWithinAt_nonpos (convex_Icc _ _) hgc.continuousOn
          (fun x _ => (hg x).hasDerivWithinAt)
        intro x hx
        rw [interior_Icc] at hx
        have : 0 ≤ (1 - 2 * x) * B1.eval x ^ 2 := mul_nonneg (by linarith [hx.2]) (sq_nonneg _)
        rw [neg_mul, neg_div]
        exact neg_nonpos.mpr (div_nonneg this hNpos.le)
      have := hA ⟨le_refl 0, by norm_num⟩ ⟨hu0, h⟩ hu0
      linarith
    · have hM : MonotoneOn g (Set.Icc (1/2) 1) := by
        apply monotoneOn_of_hasDerivWithinAt_nonneg (convex_Icc _ _) hgc.continuousOn
          (fun x _ => (hg x).hasDerivWithinAt)
        intro x hx
        rw [interior_Icc] at hx
        have : 0 ≤ (2 * x - 1) * B1.eval x ^ 2 := mul_nonneg (by linarith [hx.1]) (sq_nonneg _)
        have e : -(1 - 2 * x) * B1.eval x ^ 2 = (2 * x - 1) * B1.eval x ^ 2 := by ring
        rw [e]
        exact div_nonneg this hNpos.le
      have := hM ⟨h, hu1⟩ ⟨by norm_num, le_refl 1⟩ hu1
      linarith
  have hsq : B.eval u ^ 2 ≤ 1 := by
    have : 0 ≤ u * (1 - u) * B1.eval u ^ 2 / N := by
      apply div_nonneg _ hNpos.le
      exact mul_nonneg (mul_nonneg hu0 (by linarith)) (sq_nonneg _)
    simp only [g] at hle
    linarith
  rw [← abs_one, ← sq_le_sq]
  simpa using hsq

lemma P_abs_le (n : ℕ) (t : ℝ) (h0 : -1 ≤ t) (h1 : t ≤ 1) : |(T13.P n).eval t| ≤ 1 := by
  rw [T13.P_eval']
  exact B_abs_le n _ (by linarith) (by linarith)

lemma stmt2_aux (m : ℕ) (hm : 1 ≤ m) (θ : ℝ) :
    |T13.LamInt m θ| ≤ min (1 - Real.cos θ) (2 / (2 * (m : ℝ) + 1)) := by
  have hc1 := Real.cos_le_one θ
  have hc0 := Real.neg_one_le_cos θ
  refine le_min ?_ ?_
  · unfold T13.LamInt
    have := intervalIntegral.norm_integral_le_of_norm_le_const (a := Real.cos θ) (b := 1)
      (C := 1) (f := fun t => (T13.P m).eval t) (by
        intro x hx
        rw [Set.uIoc_of_le hc1] at hx
        exact P_abs_le m x (by linarith [hx.1]) hx.2)
    rw [Real.norm_eq_abs, abs_of_nonneg (by linarith : (0:ℝ) ≤ 1 - Real.cos θ)] at this
    simpa using this
  · rw [T13.proof_2 m hm θ, abs_div, abs_of_pos (by positivity : (0:ℝ) < 2 * m + 1)]
    apply div_le_div_of_nonneg_right _ (by positivity)
    have a := P_abs_le (m - 1) (Real.cos θ) hc0 hc1
    have b := P_abs_le (m + 1) (Real.cos θ) hc0 hc1
    calc _ ≤ |(T13.P (m - 1)).eval (Real.cos θ)| + |(T13.P (m + 1)).eval (Real.cos θ)| :=
          abs_sub _ _
      _ ≤ 2 := by linarith


lemma continuous_h (c : ℝ) : Continuous (h c) := by unfold h; fun_prop

lemma continuous_rpow_pos {p : ℝ} (hp : 0 < p) : Continuous fun x : ℝ => x ^ p :=
  Real.continuous_rpow_const hp.le

lemma defect (c n : ℝ) (hn : 1 < n) :
    IntegrableOn (fun s => h c s * s ^ (eps n - 1)) (Set.Ioc 0 n) ∧
    (∫ x in (0:ℝ)..1, h c (n * x ^ (n - 1))) = eps n * n ^ (-eps n) *
      ∫ s in Set.Ioc 0 n, h c s * s ^ (eps n - 1) := by
  set ε := eps n with hεdef
  have hn0 : 0 < n := by linarith
  have hε : 0 < ε := by rw [hεdef, eps]; apply div_pos one_pos; linarith
  have hεn : ε * (n - 1) = 1 := by rw [hεdef, eps]; field_simp [show n - 1 ≠ 0 by linarith]
  have hII : IntervalIntegrable (fun s => h c s * s ^ (ε - 1)) volume 0 n := by
    have := (intervalIntegral.intervalIntegrable_rpow' (a := 0) (b := n)
      (show -1 < ε - 1 by linarith)).continuousOn_mul (continuous_h c).continuousOn
    exact this
  have hInt : IntegrableOn (fun s => h c s * s ^ (ε - 1)) (Set.Ioc 0 n) :=
    (intervalIntegrable_iff_integrableOn_Ioc_of_le hn0.le).mp hII
  refine ⟨hInt, ?_⟩
  set g : ℝ → ℝ := fun x => h c (n * x ^ (n - 1)) with hg
  have hgc : Continuous g := by
    rw [hg]
    exact (continuous_h c).comp (continuous_const.mul (continuous_rpow_pos (by linarith)))
  set f : ℝ → ℝ := fun s => (s / n) ^ ε with hf
  set f' : ℝ → ℝ := fun s => ε * n ^ (-ε) * s ^ (ε - 1) with hf'
  have hfc : Continuous f := (continuous_rpow_pos hε).comp (continuous_id.div_const n)
  have hderiv : ∀ s, 0 < s → HasDerivAt f (f' s) s := by
    intro s hs
    have h1 : HasDerivAt (fun s : ℝ => s / n) (1 / n) s := (hasDerivAt_id s).div_const n
    have h2 := h1.rpow_const (p := ε) (Or.inl (by positivity))
    convert h2 using 1
    simp only [hf']
    rw [Real.div_rpow hs.le hn0.le, Real.rpow_neg hn0.le, Real.rpow_sub_one hn0.ne']
    have : 0 < n ^ ε := Real.rpow_pos_of_pos hn0 ε
    field_simp
  have hcomp : ∀ s, 0 ≤ s → (g ∘ f) s * f' s = ε * n ^ (-ε) * (h c s * s ^ (ε - 1)) := by
    intro s hs
    have : n * (f s) ^ (n - 1) = s := by
      simp only [hf]
      rw [← Real.rpow_mul (by positivity), hεn, Real.rpow_one]
      field_simp
    simp only [Function.comp, hg, this, hf']
    ring
  have hsub := intervalIntegral.integral_comp_mul_deriv''' (a := 0) (b := n) (f := f) (f' := f')
    (g := g) hfc.continuousOn
    (by
      intro x hx
      rw [min_eq_left hn0.le] at hx
      exact (hderiv x hx.1).hasDerivWithinAt)
    hgc.continuousOn
    (hgc.continuousOn.integrableOn_compact (isCompact_uIcc.image hfc))
    (by
      rw [Set.uIcc_of_le hn0.le, integrableOn_Icc_iff_integrableOn_Ioc]
      have hI2 : IntegrableOn (fun s => ε * n ^ (-ε) * (h c s * s ^ (ε - 1))) (Set.Ioc 0 n) :=
        hInt.const_mul _
      exact hI2.congr_fun (fun s hs => (hcomp s hs.1.le).symm) measurableSet_Ioc)
  have hf0 : f 0 = 0 := by simp [hf, Real.zero_rpow hε.ne']
  have hfn : f n = 1 := by simp [hf, div_self hn0.ne']
  rw [hf0, hfn] at hsub
  rw [← hsub, intervalIntegral.integral_congr (g := fun s => ε * n ^ (-ε) * (h c s * s ^ (ε - 1)))
    (by
      intro s hs
      rw [Set.uIcc_of_le hn0.le] at hs
      exact hcomp s hs.1),
    intervalIntegral.integral_const_mul, intervalIntegral.integral_of_le hn0.le]

lemma L_eq (θ n : ℝ) (hn : 1 < n) :
    2 - L θ n = ∫ x in (0:ℝ)..1, h (Real.cos θ) (n * x ^ (n - 1)) := by
  have hp : Continuous fun x : ℝ => x ^ (n - 1) := continuous_rpow_pos (by linarith)
  unfold h L
  have i1 : IntervalIntegrable (fun x : ℝ => 1 + n * x ^ (n - 1)) volume 0 1 :=
    (continuous_const.add (continuous_const.mul hp)).intervalIntegrable _ _
  have i2 : IntervalIntegrable (fun x : ℝ =>
      Real.sqrt (1 + 2 * Real.cos θ * (n * x ^ (n - 1)) + (n * x ^ (n - 1)) ^ 2)) volume 0 1 := by
    apply Continuous.intervalIntegrable; fun_prop
  rw [intervalIntegral.integral_sub i1 i2]
  have e : (∫ x in (0:ℝ)..1, (1 + n * x ^ (n - 1))) = 2 := by
    have hm : IntervalIntegrable (fun x : ℝ => n * x ^ (n - 1)) volume 0 1 :=
      (show Continuous (fun x : ℝ => n * x ^ (n - 1)) from continuous_const.mul hp).intervalIntegrable _ _
    rw [intervalIntegral.integral_add intervalIntegrable_const hm, intervalIntegral.integral_const_mul,
      integral_rpow (Or.inl (by linarith))]
    have : n - 1 + 1 = n := by ring
    rw [this]
    have : n ≠ 0 := by linarith
    simp [Real.zero_rpow this]
    field_simp
    ring
  rw [e]


lemma L_cont (θ : ℝ) : ContinuousOn (L θ) (Set.Ioi 1) := by
  rw [continuousOn_iff_continuous_domRestrict]
  have hp : Continuous fun q : Set.Ioi (1:ℝ) × ℝ => q.2 ^ ((q.1 : ℝ) - 1) := by
    apply Continuous.rpow continuous_snd (continuous_subtype_val.comp continuous_fst |>.sub
      continuous_const)
    intro q; right; have := q.1.2; simp only [Set.mem_Ioi] at this; show 0 < (q.1:ℝ) - 1; linarith
  have hF : Continuous (Function.uncurry fun (n : Set.Ioi (1:ℝ)) (x : ℝ) =>
      Real.sqrt (1 + 2 * Real.cos θ * ((n:ℝ) * x ^ ((n:ℝ) - 1)) + ((n:ℝ) * x ^ ((n:ℝ) - 1)) ^ 2)) := by
    have hn : Continuous fun q : Set.Ioi (1:ℝ) × ℝ => (q.1 : ℝ) :=
      continuous_subtype_val.comp continuous_fst
    show Continuous fun q : Set.Ioi (1:ℝ) × ℝ =>
      Real.sqrt (1 + 2 * Real.cos θ * ((q.1:ℝ) * q.2 ^ ((q.1:ℝ) - 1)) + ((q.1:ℝ) * q.2 ^ ((q.1:ℝ) - 1)) ^ 2)
    fun_prop
  exact intervalIntegral.continuous_parametric_intervalIntegral_of_continuous' hF 0 1


lemma hasSum_z3 : HasSum (fun n : ℕ => 1 / ((n : ℝ) + 1) ^ 3) z3 := by
  have h1 : riemannZeta 3 = ((∑' n : ℕ, 1 / ((n : ℝ)) ^ 3 : ℝ) : ℂ) := by
    rw [show (3 : ℂ) = ((3 : ℕ) : ℂ) by norm_num, zeta_nat_eq_tsum_of_gt_one (by norm_num),
      Complex.ofReal_tsum]
    push_cast; rfl
  have h2 : z3 = ∑' n : ℕ, 1 / ((n : ℝ)) ^ 3 := by
    rw [z3, h1, Complex.ofReal_re]
  have hs : Summable (fun n : ℕ => 1 / ((n : ℝ)) ^ 3) :=
    (Real.summable_one_div_nat_pow (p := 3)).2 (by norm_num)
  rw [h2, hs.tsum_eq_zero_add]
  simp only [CharP.cast_eq_zero, ne_eq, OfNat.ofNat_ne_zero, not_false_eq_true, zero_pow,
    div_zero, zero_add, Nat.cast_add, Nat.cast_one]
  convert ((summable_nat_add_iff 1).2 hs).hasSum using 1 <;> simp

/-- `Q(t) = Σ t^n/(n+1)^2 = Li₂(t)/t`. -/
noncomputable def Qs (t : ℝ) : ℝ := ∑' n : ℕ, t ^ n / ((n : ℝ) + 1) ^ 2

lemma qterm_bound {x : ℝ} (hx : |x| ≤ 1) (n : ℕ) :
    ‖x ^ n / ((n : ℝ) + 1) ^ 2‖ ≤ 1 / ((n : ℝ) + 1) ^ 2 := by
  rw [Real.norm_eq_abs, abs_div, abs_pow, abs_of_pos (by positivity : (0:ℝ) < ((n:ℝ) + 1) ^ 2)]
  apply div_le_div_of_nonneg_right _ (by positivity)
  exact pow_le_one₀ (abs_nonneg x) hx

lemma continuousOn_Qs : ContinuousOn Qs (Set.Icc (-1) 1) := by
  apply continuousOn_tsum (fun n => (continuous_pow n).continuousOn.div_const _)
    D00.summable_inv_sq
  intro n x hx
  exact qterm_bound (abs_le.mpr hx) n

lemma Qs_eq {t : ℝ} (ht : |t| ≤ 1) (h0 : t ≠ 0) : Qs t = D00.Li2 t / t := by
  have h := (D00.proof_6 t ht).div_const t
  have e : (fun n : ℕ => t ^ (n + 1) / ((n : ℝ) + 1) ^ 2 / t) = fun n : ℕ => t ^ n / ((n:ℝ) + 1) ^ 2 := by
    funext n; rw [pow_succ]; field_simp
  rw [e] at h
  exact h.tsum_eq

lemma Li3_eq {x : ℝ} (hx : |x| ≤ 1) : Li3 x = ∫ t in (0 : ℝ)..x, Qs t := by
  rw [Li3, intervalIntegral.integral_congr_uIoo]
  intro t ht
  have := D00.mem_uIoo_zero ht
  exact (Qs_eq (by linarith [this.2]) this.1).symm

lemma Qs_intervalIntegrable {a b : ℝ} (ha : |a| ≤ 1) (hb : |b| ≤ 1) :
    IntervalIntegrable Qs volume a b := by
  apply ContinuousOn.intervalIntegrable
  apply continuousOn_Qs.mono
  intro t ht
  rw [Set.mem_uIcc] at ht
  rw [abs_le] at ha hb
  constructor <;> rcases ht with h | h <;> linarith [h.1, h.2]

lemma Li3_one : Li3 1 = z3 := by
  rw [Li3_eq (by norm_num)]
  have key := intervalIntegral.hasSum_integral_of_dominated_convergence (a := 0) (b := 1)
    (μ := volume) (F := fun (n : ℕ) (t : ℝ) => t ^ n / ((n : ℝ) + 1) ^ 2) (f := Qs)
    (fun n _ => 1 / ((n : ℝ) + 1) ^ 2)
    (fun n => (Continuous.aestronglyMeasurable (by fun_prop)))
    (by
      intro n
      filter_upwards with t ht
      rw [Set.uIoc_of_le zero_le_one] at ht
      exact qterm_bound (abs_le.mpr ⟨by linarith [ht.1], ht.2⟩) n)
    (by filter_upwards with t _ using D00.summable_inv_sq)
    intervalIntegrable_const
    (by
      filter_upwards with t ht
      rw [Set.uIoc_of_le zero_le_one] at ht
      exact (Summable.of_norm_bounded D00.summable_inv_sq
        (qterm_bound (abs_le.mpr ⟨by linarith [ht.1], ht.2⟩))).hasSum)
  have e : ∀ n : ℕ, (∫ t in (0:ℝ)..1, t ^ n / ((n : ℝ) + 1) ^ 2) = 1 / ((n:ℝ) + 1) ^ 3 := by
    intro n
    rw [intervalIntegral.integral_div, integral_pow]
    field_simp
    ring
  simp only [e] at key
  exact key.unique hasSum_z3

lemma hasDerivAt_Li3 {x : ℝ} (hx : |x| < 1) : HasDerivAt Li3 (Qs x) x := by
  have hx' := abs_lt.mp hx
  have hev : (fun y => ∫ t in (0:ℝ)..y, Qs t) =ᶠ[𝓝 x] Li3 := by
    filter_upwards [Ioo_mem_nhds hx'.1 hx'.2] with y hy
    exact (Li3_eq (abs_le.mpr ⟨hy.1.le, hy.2.le⟩)).symm
  have hco : ContinuousOn Qs (Set.Ioo (-1) 1) := continuousOn_Qs.mono Set.Ioo_subset_Icc_self
  refine HasDerivAt.congr_of_eventuallyEq ?_ hev.symm
  exact intervalIntegral.integral_hasDerivAt_right (Qs_intervalIntegrable (by norm_num) hx.le)
    (hco.stronglyMeasurableAtFilter isOpen_Ioo x ⟨hx'.1, hx'.2⟩)
    (hco.continuousAt (Ioo_mem_nhds hx'.1 hx'.2))

lemma continuousOn_Li3 : ContinuousOn Li3 (Set.Icc 0 1) := by
  have hI : IntegrableOn Qs (Set.uIcc 0 1) := by
    apply ContinuousOn.integrableOn_compact isCompact_uIcc
    apply continuousOn_Qs.mono
    rw [Set.uIcc_of_le zero_le_one]
    exact Set.Icc_subset_Icc (by norm_num) le_rfl
  have h := intervalIntegral.continuousOn_primitive_interval hI
  rw [Set.uIcc_of_le zero_le_one] at h
  refine h.congr (fun x hx => ?_)
  exact Li3_eq (abs_le.mpr ⟨by linarith [hx.1], hx.2⟩)

lemma integrand_eq (t : ℝ) : Real.log (1 - t) ^ 2 / t = D00.G t * Real.log (1 - t) := by
  rcases eq_or_ne t 0 with rfl | h
  · simp
  · rw [D00.G_of_ne h]; ring

lemma log_term_eq (w : ℝ) :
    Real.log (1 - w) ^ 2 * Real.log w = D00.G w ^ 2 * (w * (w * Real.log w)) := by
  rcases eq_or_ne w 0 with rfl | h
  · simp
  · rw [D00.G_of_ne h]; field_simp

lemma continuousAt_g {t : ℝ} (ht : t < 1) :
    ContinuousAt (fun t => D00.G t * Real.log (1 - t)) t := by
  apply (D00.continuousAt_G ht.ne).mul
  exact (continuousAt_const.sub continuousAt_id).log (sub_pos.mpr ht).ne'

lemma S12_identity {z : ℝ} (h0 : 0 < z) (h1 : z < 1) :
    z3 - Li3 (1 - z) + Real.log (1 - z) * D00.Li2 (1 - z) + Real.log (1 - z) ^ 2 * Real.log z / 2
      = (1 / 2) * ∫ t in (0 : ℝ)..z, Real.log (1 - t) ^ 2 / t := by
  simp_rw [integrand_eq]
  set g : ℝ → ℝ := fun t => D00.G t * Real.log (1 - t) with hg
  have gco : ContinuousOn g (Set.Iio 1) := fun t ht => (continuousAt_g ht).continuousWithinAt
  have gint : ∀ w, w < 1 → IntervalIntegrable g volume 0 w := by
    intro w hw
    apply ContinuousOn.intervalIntegrable
    apply gco.mono
    intro t ht
    rw [Set.mem_uIcc] at ht
    show t < 1
    rcases ht with h | h <;> linarith [h.1, h.2]
  set F : ℝ → ℝ := fun w => z3 - Li3 (1 - w) + Real.log (1 - w) * D00.Li2 (1 - w)
    + Real.log (1 - w) ^ 2 * Real.log w / 2 - (1 / 2) * ∫ t in (0 : ℝ)..w, g t with hF
  suffices F 0 = F z by
    have hF0 : F 0 = 0 := by simp [hF, Li3_one]
    simp only [hF] at this hF0 ⊢
    linarith
  apply D00.eq_of_deriv_zero h0.le
  · -- continuity
    have e : F = fun w => z3 - Li3 (1 - w) + Real.log (1 - w) * D00.Li2 (1 - w)
        + D00.G w ^ 2 * (w * (w * Real.log w)) / 2 - (1 / 2) * ∫ t in (0 : ℝ)..w, g t := by
      funext w; simp only [hF, log_term_eq]
    rw [e]
    have hIcc : ∀ w ∈ Set.Icc 0 z, w < 1 := fun w hw => lt_of_le_of_lt hw.2 h1
    have c1 : ContinuousOn (fun w : ℝ => Li3 (1 - w)) (Set.Icc 0 z) := by
      refine continuousOn_Li3.comp (continuous_const.sub continuous_id).continuousOn ?_
      intro w hw; exact ⟨by linarith [hIcc w hw], by linarith [hw.1]⟩
    have c2 : ContinuousOn (fun w : ℝ => Real.log (1 - w)) (Set.Icc 0 z) := by
      intro w hw
      exact ((continuousAt_const.sub continuousAt_id).log
        (sub_pos.mpr (hIcc w hw)).ne').continuousWithinAt
    have c3 : ContinuousOn (fun w : ℝ => D00.G w) (Set.Icc 0 z) :=
      fun w hw => (D00.continuousAt_G (hIcc w hw).ne).continuousWithinAt
    have c4 : ContinuousOn (fun w => ∫ t in (0 : ℝ)..w, g t) (Set.Icc 0 z) := by
      have := intervalIntegral.continuousOn_primitive_interval (a := (0:ℝ)) (b := z) (μ := volume)
        (f := g) (by
          apply ContinuousOn.integrableOn_compact isCompact_uIcc
          apply gco.mono
          rw [Set.uIcc_of_le h0.le]; exact fun w hw => hIcc w hw)
      rwa [Set.uIcc_of_le h0.le] at this
    have c5 : Continuous (fun w : ℝ => w * Real.log w) := Real.continuous_mul_log
    have c6 : ContinuousOn (fun w : ℝ => D00.Li2 (1 - w)) (Set.Icc 0 z) :=
      (D00.continuous_Li2.comp (continuous_const.sub continuous_id)).continuousOn
    exact (((continuousOn_const.sub c1).add (c2.mul c6)).add
      (((c3.pow 2).mul (continuous_id.mul c5).continuousOn).div_const 2)).sub
      (c4.const_smul (1/2 : ℝ))
  · intro w hw
    have hw1 : w < 1 := by linarith [hw.2]
    have hw0 : w ≠ 0 := hw.1.ne'
    have hsub : HasDerivAt (fun w : ℝ => 1 - w) (-1) w := by
      simpa using (hasDerivAt_id w).const_sub 1
    have d1 : HasDerivAt (fun w => Li3 (1 - w)) (Qs (1 - w) * (-1)) w :=
      (hasDerivAt_Li3 (by rw [abs_lt]; constructor <;> linarith [hw.1])).comp w hsub
    have d2 : HasDerivAt (fun w => Real.log (1 - w)) ((-1) / (1 - w)) w :=
      hsub.log (sub_pos.mpr hw1).ne'
    have d3 : HasDerivAt (fun w => D00.Li2 (1 - w)) (-Real.log (1 - (1 - w)) / (1 - w) * (-1)) w :=
      D00.Li2_comp hsub (by intro h; linarith [hw.1]) (sub_pos.mpr hw1).ne'
    have d4 : HasDerivAt Real.log w⁻¹ w := Real.hasDerivAt_log hw0
    have d5 : HasDerivAt (fun u => ∫ t in (0:ℝ)..u, g t) (g w) w :=
      intervalIntegral.integral_hasDerivAt_right (gint w hw1)
        (gco.stronglyMeasurableAtFilter isOpen_Iio w hw1) (continuousAt_g hw1)
    have hall := (((((hasDerivAt_const w z3).sub d1).add (d2.mul d3)).add
      (((d2.pow 2).mul d4).div_const 2)).sub (d5.const_mul (1/2 : ℝ)))
    convert hall using 1
    rw [Qs_eq (by rw [abs_le]; constructor <;> linarith [hw.1]) (by linarith)]
    simp only [hg, ← integrand_eq, Pi.pow_apply]
    have : (1:ℝ) - (1 - w) = w := by ring
    rw [this]
    have h1w : (1:ℝ) - w ≠ 0 := by linarith
    field_simp
    ring_nf


lemma bonnet_eval (n : ℕ) (c : ℝ) :
    ((n:ℝ) + 2) * (T13.P (n + 2)).eval c
      = (2 * (n:ℝ) + 3) * c * (T13.P (n + 1)).eval c - ((n:ℝ) + 1) * (T13.P n).eval c := by
  have h := congrArg (eval c) (T13.P_bonnet n)
  simp only [eval_smul, eval_mul, eval_sub, eval_X, smul_eq_mul] at h
  linear_combination h

lemma lam_mul (n : ℕ) (θ : ℝ) :
    T13.LamInt (n + 1) θ * (2 * (n:ℝ) + 3)
      = (T13.P n).eval (Real.cos θ) - (T13.P (n + 2)).eval (Real.cos θ) := by
  rw [T13.lam_succ, div_mul_eq_mul_div, div_eq_iff (by positivity)]
  push_cast; ring

lemma lam_rec (j : ℕ) (θ : ℝ) :
    ((j:ℝ) + 4) * T13.LamInt (j + 3) θ - (2 * (j:ℝ) + 5) * Real.cos θ * T13.LamInt (j + 2) θ
      + ((j:ℝ) + 1) * T13.LamInt (j + 1) θ = 0 := by
  set c := Real.cos θ
  have hL1 := lam_mul j θ
  have hL2 := lam_mul (j + 1) θ
  have hL3 := lam_mul (j + 2) θ
  have B1 := bonnet_eval j c
  have B2 := bonnet_eval (j + 2) c
  push_cast at hL2 hL3 B2
  have key : (2 * (j:ℝ) + 3) * (2 * (j:ℝ) + 7) *
      (((j:ℝ) + 4) * T13.LamInt (j + 3) θ - (2 * (j:ℝ) + 5) * c * T13.LamInt (j + 2) θ
      + ((j:ℝ) + 1) * T13.LamInt (j + 1) θ) = 0 := by
    have e1 : j + 1 + 2 = j + 3 := by omega
    have e2 : j + 2 + 2 = j + 4 := by omega
    have e3 : j + 1 + 1 = j + 2 := by omega
    have e4 : j + 2 + 1 = j + 3 := by omega
    rw [e1] at hL2; rw [e2] at hL3; rw [e2, e4] at B2
    rw [show j + 2 + 1 = j + 3 by omega] at hL3
    linear_combination (2 * (j:ℝ) + 3) * ((j:ℝ) + 4) * hL3
      - (2 * (j:ℝ) + 3) * (2 * (j:ℝ) + 7) * c * hL2 + (2 * (j:ℝ) + 7) * ((j:ℝ) + 1) * hL1
      + (2 * (j:ℝ) + 7) * B1 - (2 * (j:ℝ) + 3) * B2
  have hne : (2 * (j:ℝ) + 3) * (2 * (j:ℝ) + 7) ≠ 0 := by positivity
  exact (mul_eq_zero.mp key).resolve_left hne

/-- Taylor coefficients of `√(1+2ct+t²)`: `1, c, (-1)^m Λ_{m+1}`. -/
noncomputable def gc (θ : ℝ) : ℕ → ℝ
  | 0 => 1
  | 1 => Real.cos θ
  | (m + 2) => (-1) ^ m * T13.LamInt (m + 1) θ

lemma gc_rec (θ : ℝ) (m : ℕ) :
    ((m:ℝ) + 2) * gc θ (m + 2) + (2 * (m:ℝ) + 1) * Real.cos θ * gc θ (m + 1)
      + ((m:ℝ) - 1) * gc θ m = 0 := by
  set c := Real.cos θ
  rcases m with _ | _ | j
  · simp only [gc, pow_zero, one_mul, Nat.cast_zero]
    rw [T13.proof_9]
    nlinarith [Real.sin_sq_add_cos_sq θ]
  · simp only [gc, pow_one, pow_zero, one_mul]
    have hL1 := lam_mul 0 θ
    have hL2 := lam_mul 1 θ
    have B0 := bonnet_eval 0 c
    have B1 := bonnet_eval 1 c
    simp only [T13.P_zero, T13.P_one, eval_one, eval_X, Nat.cast_zero, Nat.cast_one] at hL1 hL2 B0 B1
    norm_num at hL1 hL2 B0 B1 ⊢
    linear_combination (-(3:ℝ)/5) * hL2 + c * hL1 + ((1:ℝ)/5) * B1
  · have hr := lam_rec j θ
    simp only [gc]
    push_cast
    have e : j + 1 + 1 = j + 2 := by omega
    rw [show j + 2 + 1 = j + 3 by omega, show j + 1 + 1 = j + 2 by omega] at *
    rw [pow_succ, pow_succ]
    linear_combination (-1:ℝ) ^ j * hr


section
variable (a : ℕ → ℝ) (c : ℝ)

lemma ps_summable (hb : ∀ m, |a m| ≤ 1) {t : ℝ} (ht : |t| < 1) :
    Summable (fun m => a m * t ^ m) := by
  refine Summable.of_norm_bounded (summable_geometric_of_lt_one (abs_nonneg t) ht) (fun m => ?_)
  rw [Real.norm_eq_abs, abs_mul, abs_pow]
  exact mul_le_of_le_one_left (by positivity) (hb m)

lemma ps_deriv (hb : ∀ m, |a m| ≤ 1) {t : ℝ} (ht : |t| < 1) :
    HasDerivAt (fun x => ∑' m, a m * x ^ m) (∑' m : ℕ, ((m:ℝ) + 1) * a (m + 1) * t ^ m) t := by
  set r := (|t| + 1) / 2 with hr
  have hr1 : r < 1 := by rw [hr]; linarith
  have hr0 : 0 < r := by rw [hr]; positivity
  have htr : t ∈ Set.Ioo (-r) r :=
    ⟨by linarith [neg_abs_le t, hr], by linarith [le_abs_self t, hr]⟩
  have hsum_u : Summable (fun n : ℕ => (n : ℝ) * r ^ n / r) := by
    have := (hasSum_coe_mul_geometric_of_norm_lt_one (𝕜 := ℝ)
      (by rw [Real.norm_eq_abs, abs_of_pos hr0]; exact hr1)).summable
    exact this.div_const r
  have hd := hasDerivAt_tsum_of_isPreconnected (t := Set.Ioo (-r) r) (y₀ := 0)
    (g := fun (n : ℕ) (z : ℝ) => a n * z ^ n)
    (g' := fun (n : ℕ) (z : ℝ) => a n * ((n : ℝ) * z ^ (n - 1)))
    (u := fun n : ℕ => (n : ℝ) * r ^ n / r)
    hsum_u isOpen_Ioo isPreconnected_Ioo
    (by
      intro n z _
      have := (hasDerivAt_pow n z).const_mul (a n)
      simpa using this)
    (by
      intro n z hz
      have hz' : |z| ≤ r := abs_le.mpr ⟨hz.1.le, hz.2.le⟩
      rw [Real.norm_eq_abs, abs_mul, abs_mul, abs_pow, Nat.abs_cast]
      rcases n with _ | n
      · simp
      · simp only [Nat.add_sub_cancel]
        have : (((n + 1 : ℕ) : ℝ)) * r ^ (n + 1) / r = ((n + 1 : ℕ) : ℝ) * r ^ n := by
          rw [pow_succ]; field_simp
        rw [this]
        calc |a (n + 1)| * (((n + 1 : ℕ) : ℝ) * |z| ^ n) ≤ 1 * (((n + 1 : ℕ) : ℝ) * r ^ n) := by
              apply mul_le_mul (hb _) _ (by positivity) zero_le_one
              exact mul_le_mul_of_nonneg_left (pow_le_pow_left₀ (abs_nonneg z) hz' n) (by positivity)
          _ = _ := one_mul _)
    ⟨by linarith, hr0⟩
    (by simpa using ps_summable a hb (t := 0) (by simp))
    htr
  convert hd using 1
  have hs : Summable (fun n : ℕ => a n * ((n : ℝ) * t ^ (n - 1))) := by
    refine Summable.of_norm_bounded hsum_u (fun n => ?_)
    have hz' : |t| ≤ r := abs_le.mpr ⟨htr.1.le, htr.2.le⟩
    rw [Real.norm_eq_abs, abs_mul, abs_mul, abs_pow, Nat.abs_cast]
    rcases n with _ | n
    · simp
    · simp only [Nat.add_sub_cancel]
      have : (((n + 1 : ℕ) : ℝ)) * r ^ (n + 1) / r = ((n + 1 : ℕ) : ℝ) * r ^ n := by
        rw [pow_succ]; field_simp
      rw [this]
      calc |a (n + 1)| * (((n + 1 : ℕ) : ℝ) * |t| ^ n) ≤ 1 * (((n + 1 : ℕ) : ℝ) * r ^ n) := by
            apply mul_le_mul (hb _) _ (by positivity) zero_le_one
            exact mul_le_mul_of_nonneg_left (pow_le_pow_left₀ (abs_nonneg t) hz' n) (by positivity)
        _ = _ := one_mul _
  rw [hs.tsum_eq_zero_add]
  simp only [CharP.cast_eq_zero, zero_mul, mul_zero, zero_add, Nat.add_sub_cancel]
  congr 1; funext n; push_cast; ring

lemma psd_summable (hb : ∀ m, |a m| ≤ 1) {t : ℝ} (ht : |t| < 1) :
    Summable (fun m : ℕ => ((m:ℝ) + 1) * a (m + 1) * t ^ m) := by
  have h1 := (hasSum_coe_mul_geometric_of_norm_lt_one (𝕜 := ℝ)
      (by rw [Real.norm_eq_abs, abs_abs]; exact ht)).summable
  have h2 := summable_geometric_of_lt_one (abs_nonneg t) ht
  refine Summable.of_norm_bounded (h1.add h2) (fun m => ?_)
  rw [Real.norm_eq_abs, abs_mul, abs_mul, abs_pow, abs_of_pos (by positivity : (0:ℝ) < m + 1)]
  have := hb (m + 1)
  calc ((m:ℝ) + 1) * |a (m + 1)| * |t| ^ m ≤ ((m:ℝ) + 1) * 1 * |t| ^ m := by gcongr
    _ = _ := by ring

lemma ps_ode (hb : ∀ m, |a m| ≤ 1) (h0 : a 0 = 1) (h1 : a 1 = c)
    (hrec : ∀ m : ℕ, ((m:ℝ) + 2) * a (m + 2) + (2 * (m:ℝ) + 1) * c * a (m + 1)
      + ((m:ℝ) - 1) * a m = 0) {t : ℝ} (ht : |t| < 1) :
    (1 + 2 * c * t + t ^ 2) * (∑' m : ℕ, ((m:ℝ) + 1) * a (m + 1) * t ^ m)
      = (c + t) * ∑' m, a m * t ^ m := by
  set D := ∑' m : ℕ, ((m:ℝ) + 1) * a (m + 1) * t ^ m
  set f := ∑' m, a m * t ^ m
  have hD : HasSum (fun m : ℕ => ((m:ℝ) + 1) * a (m + 1) * t ^ m) D :=
    (psd_summable a hb ht).hasSum
  have hf : HasSum (fun m => a m * t ^ m) f := (ps_summable a hb ht).hasSum
  have sD2 := (hasSum_nat_add_iff' 2).mpr hD
  have sD1 := ((hasSum_nat_add_iff' 1).mpr hD).mul_left (2 * c * t)
  have sD0 := hD.mul_left (t ^ 2)
  have sf2 := ((hasSum_nat_add_iff' 2).mpr hf).mul_left c
  have sf1 := ((hasSum_nat_add_iff' 1).mpr hf).mul_left t
  have hall := (((sD2.add sD1).add sD0).sub sf2).sub sf1
  have hz : (fun m : ℕ => ((m + 2 : ℕ) + 1 : ℝ) * a (m + 2 + 1) * t ^ (m + 2)
      + 2 * c * t * (((m + 1 : ℕ) + 1 : ℝ) * a (m + 1 + 1) * t ^ (m + 1))
      + t ^ 2 * (((m : ℕ) + 1 : ℝ) * a (m + 1) * t ^ m)
      - c * (a (m + 2) * t ^ (m + 2)) - t * (a (m + 1) * t ^ (m + 1))) = fun _ => 0 := by
    funext m
    have := hrec (m + 1)
    push_cast at this ⊢
    rw [show m + 1 + 2 = m + 3 by omega, show m + 1 + 1 = m + 2 by omega] at this
    rw [show m + 2 + 1 = m + 3 by omega, show m + 1 + 1 = m + 2 by omega]
    linear_combination t ^ (m + 2) * this
  simp only [Finset.sum_range_succ, Finset.sum_range_zero] at hall
  rw [hz] at hall
  have e := hall.unique hasSum_zero
  have r0 := hrec 0
  simp only [Nat.cast_zero, h0, h1, zero_add] at e r0
  push_cast at e
  linear_combination e + t * r0

lemma q_pos (hc : |c| ≤ 1) {x : ℝ} (hx : |x| < 1) : 0 < 1 + 2 * c * x + x ^ 2 := by
  have h1 : -(|x|) ≤ c * x := by
    have := abs_mul c x
    have h2 : |c * x| ≤ |x| := by rw [this]; exact mul_le_of_le_one_left (abs_nonneg x) hc
    linarith [neg_abs_le (c * x)]
  have h3 : 0 < (1 - |x|) ^ 2 := by have : 0 < 1 - |x| := by linarith
                                    positivity
  have h4 : x ^ 2 = |x| ^ 2 := (sq_abs x).symm
  nlinarith

lemma ps_sqrt (hc : |c| ≤ 1) (hb : ∀ m, |a m| ≤ 1) (h0 : a 0 = 1) (h1 : a 1 = c)
    (hrec : ∀ m : ℕ, ((m:ℝ) + 2) * a (m + 2) + (2 * (m:ℝ) + 1) * c * a (m + 1)
      + ((m:ℝ) - 1) * a m = 0) {t : ℝ} (ht : |t| < 1) :
    HasSum (fun m => a m * t ^ m) (Real.sqrt (1 + 2 * c * t + t ^ 2)) := by
  set F : ℝ → ℝ := fun x => ∑' m, a m * x ^ m with hF
  set φ : ℝ → ℝ := fun x => F x / Real.sqrt (1 + 2 * c * x + x ^ 2) with hφ
  have hφd : ∀ x, |x| < 1 → HasDerivAt φ 0 x := by
    intro x hx
    have hq := q_pos c hc hx
    have hdq : HasDerivAt (fun x => 1 + 2 * c * x + x ^ 2) (2 * c + 2 * x) x := by
      have h := ((hasDerivAt_id' x).const_mul (2 * c)).const_add 1 |>.add (hasDerivAt_pow 2 x)
      exact h.congr_deriv (by norm_num)
    have hg := hdq.sqrt hq.ne'
    have hFd := ps_deriv a hb hx
    have hdiv := hFd.div hg (Real.sqrt_pos.mpr hq).ne'
    convert hdiv using 1
    have hode := ps_ode a c hb h0 h1 hrec hx
    have hs := Real.sq_sqrt hq.le
    have hsp := (Real.sqrt_pos.mpr hq).ne'
    field_simp
    linear_combination (-1) * hode - (∑' m : ℕ, ((m:ℝ) + 1) * a (m + 1) * x ^ m) * hs
  have hφ0 : φ 0 = 1 := by
    simp only [hφ, hF]
    rw [tsum_eq_single 0 (fun m hm => by simp [zero_pow hm])]
    simp [h0]
  have hcont : ∀ x, |x| < 1 → ContinuousAt φ x := fun x hx => (hφd x hx).continuousAt
  have hmem : ∀ x ∈ Set.uIcc 0 t, |x| < 1 := by
    intro x hx
    rw [Set.mem_uIcc] at hx
    rw [abs_lt] at ht ⊢
    constructor <;> rcases hx with h | h <;> linarith [h.1, h.2]
  have key := D00.eq_of_deriv_zero' (f := φ) (a := 0) (b := t)
    (fun x hx => (hcont x (hmem x hx)).continuousWithinAt)
    (fun x hx => hφd x (hmem x (Set.Ioo_subset_Icc_self hx)))
  rw [hφ0] at key
  have hq := q_pos c hc ht
  have : F t = Real.sqrt (1 + 2 * c * t + t ^ 2) := by
    have := key.symm
    simp only [hφ] at this
    rw [div_eq_one_iff_eq (Real.sqrt_pos.mpr hq).ne'] at this
    exact this
  rw [← this]
  exact (ps_summable a hb ht).hasSum

end



/-- The removable-singularity claim follows from the length formula and continuity of `L`. -/
lemma stmt8_of_stmt6 (h6 : stmt_6) : stmt_8 := by
  intro θ h0 h1 j hj
  set n0 : ℝ := 1 + 1 / (j : ℝ) with hn0
  have hjpos : (0:ℝ) < j := by exact_mod_cast hj
  have hn01 : 1 < n0 := by rw [hn0]; have := one_div_pos.mpr hjpos; linarith
  have hepsn0 : eps n0 = j := by rw [eps, hn0]; field_simp; ring
  have hcont_eps : ContinuousAt eps n0 := by
    unfold eps
    exact continuousAt_const.div (continuousAt_id.sub continuousAt_const) (sub_pos.mpr hn01).ne'
  have hL : ContinuousAt (L θ) n0 :=
    (L_cont θ).continuousAt (Ioi_mem_nhds hn01)
  have hev1 : ∀ᶠ n in 𝓝 n0, |eps n - j| < 1 / 2 ∧ 1 < n := by
    have h2 : ∀ᶠ n in 𝓝 n0, eps n ∈ Set.Ioo ((j:ℝ) - 1 / 2) (j + 1 / 2) :=
      hcont_eps.preimage_mem_nhds (by rw [hepsn0]; exact Ioo_mem_nhds (by linarith) (by linarith))
    filter_upwards [h2, Ioi_mem_nhds hn01] with n hn hn'
    exact ⟨abs_lt.mpr ⟨by linarith [hn.1], by linarith [hn.2]⟩, hn'⟩
  have hev : ∀ᶠ n in 𝓝[≠] n0, RHS θ n = L θ n := by
    filter_upwards [nhdsWithin_le_nhds hev1, self_mem_nhdsWithin] with n hn hne
    obtain ⟨hclose, hn1⟩ := hn
    have hnat : ∀ i : ℕ, eps n ≠ i := by
      intro i hi
      by_cases hij : i = j
      · subst hij
        apply hne
        have : n - 1 = 1 / (i:ℝ) := by
          rw [eps] at hi
          have hn1' : n - 1 ≠ 0 := by linarith
          field_simp at hi ⊢; linarith
        simp only [Set.mem_singleton_iff, hn0]; linarith
      · rw [hi] at hclose
        have : (1:ℝ) ≤ |(i:ℝ) - j| := by
          rcases lt_or_gt_of_ne hij with h | h
          · have : (i:ℝ) + 1 ≤ j := by exact_mod_cast h
            rw [abs_of_neg (by linarith)]; linarith
          · have : (j:ℝ) + 1 ≤ i := by exact_mod_cast h
            rw [abs_of_pos (by linarith)]; linarith
        linarith
    have H := h6 θ n h0 h1 hn1 hnat
    rw [RHS, H.tsum_eq]
    ring
  exact (hL.tendsto.mono_left nhdsWithin_le_nhds).congr' (hev.mono fun n hn => hn.symm)

/-- thm:P2asym follows from the uniform remainder bound. -/
lemma stmt14_of_stmt13 (h13 : stmt_13) : stmt_14 := by
  obtain ⟨C, hC⟩ := h13
  intro θ h0 h1
  refine ⟨C * Real.sin (θ / 2) ^ 2, 9, fun n hn => ?_⟩
  have := hC θ n h0 h1 hn
  simpa [E, zz, mul_assoc] using this

section Mellin
open Set

lemma binom_half {x : ℝ} (hx : |x| < 1) :
    HasSum (fun n : ℕ => Ring.choose (1 / 2 : ℝ) n * x ^ n) (Real.sqrt (1 + x)) := by
  have H := (Real.one_add_rpow_hasFPowerSeriesOnBall_zero (a := (1 / 2 : ℝ))).hasSum
    (y := x) (by simpa [edist_dist, Real.dist_eq] using hx)
  have h1 : 0 ≤ 1 + x := by linarith [neg_abs_le x]
  rw [Real.sqrt_eq_rpow]
  simpa [binomialSeries, FormalMultilinearSeries.coeff_ofScalars, mul_comm] using H

lemma k_series (z ε s : ℝ) (hz0 : 0 ≤ z) (hz1 : z < 1) (hs : 0 < s) :
    HasSum (fun j : ℕ => (Ring.choose (1 / 2 : ℝ) (j + 2) * (-4) ^ (j + 2) * z ^ (j + 2)) *
        (s ^ (((j:ℝ) + 2 + ε) - 1) * (1 + s) ^ (-(((j:ℝ) + 2 + ε) + ((j:ℝ) + 1 - ε)))))
      (k (1 - 2 * z) s * s ^ (ε - 1) - 2 * z * (s ^ (ε - 1) * (1 + s) ^ (-(ε + (1 - ε))))) := by
  set x := -(z * (4 * s / (1 + s) ^ 2)) with hxdef
  have hs1 : 0 < 1 + s := by linarith
  have hu0 : 0 ≤ 4 * s / (1 + s) ^ 2 := by positivity
  have hu1 : 4 * s / (1 + s) ^ 2 ≤ 1 := by
    rw [div_le_one (by positivity)]; nlinarith [sq_nonneg (1 - s)]
  have hx : |x| < 1 := by
    rw [hxdef, abs_neg, abs_of_nonneg (mul_nonneg hz0 hu0)]
    calc z * (4 * s / (1 + s) ^ 2) ≤ z * 1 := mul_le_mul_of_nonneg_left hu1 hz0
      _ < 1 := by linarith
  have H := ((hasSum_nat_add_iff' 2).mpr (binom_half hx)).mul_left ((1 + s) * s ^ (ε - 1))
  convert H using 1
  · funext j
    have e1 : s ^ (((j:ℝ) + 2 + ε) - 1) = s ^ (ε - 1) * s ^ (j + 2) := by
      rw [← Real.rpow_natCast, ← Real.rpow_add hs]; push_cast; ring_nf
    have e2 : (1 + s) ^ (-(((j:ℝ) + 2 + ε) + ((j:ℝ) + 1 - ε))) = 1 / (1 + s) ^ (2 * j + 3) := by
      rw [Real.rpow_neg hs1.le, one_div, ← Real.rpow_natCast]; push_cast; ring_nf
    have e5 : x ^ (j + 2) = (-4) ^ (j + 2) * z ^ (j + 2) * s ^ (j + 2) / ((1 + s) ^ (2 * j + 3) * (1 + s)) := by
      rw [hxdef, show -(z * (4 * s / (1 + s) ^ 2)) = (-4) * z * s / (1 + s) ^ 2 by ring, div_pow,
        mul_pow, mul_pow, ← pow_mul, ← pow_succ]
      ring_nf
    rw [e1, e2, e5]
    field_simp
  · have e3 : (1 + s) ^ (-(ε + (1 - ε))) = 1 / (1 + s) := by
      rw [show -(ε + (1 - ε)) = (-1 : ℝ) by ring, Real.rpow_neg_one, one_div]
    have e4 : 1 + 2 * (1 - 2 * z) * s + s ^ 2 = (1 + s) ^ 2 * (1 + x) := by
      rw [hxdef]; field_simp; ring
    rw [e3, k, e4, Real.sqrt_mul (sq_nonneg _), Real.sqrt_sq hs1.le]
    simp only [Finset.sum_range_succ, Finset.sum_range_zero, Ring.choose_zero_right,
      Ring.choose_one_right, pow_zero, pow_one, hxdef]
    field_simp
    ring

lemma beta_real {α β : ℝ} (ha : 0 < α) (hb : 0 < β) :
    ∫ t in (0:ℝ)..1, t ^ (α - 1) * (1 - t) ^ (β - 1) = Gamma α * Gamma β / Gamma (α + β) := by
  have h := Complex.Gamma_mul_Gamma_eq_betaIntegral (s := (α:ℂ)) (t := (β:ℂ))
    (by simpa using ha) (by simpa using hb)
  have hB : Complex.betaIntegral α β = ((∫ t in (0:ℝ)..1, t ^ (α - 1) * (1 - t) ^ (β - 1) : ℝ) : ℂ) := by
    rw [Complex.betaIntegral, ← intervalIntegral.integral_ofReal]
    apply intervalIntegral.integral_congr
    intro t ht
    rw [uIcc_of_le zero_le_one] at ht
    simp only
    rw [Complex.ofReal_mul, Complex.ofReal_cpow ht.1, Complex.ofReal_cpow (by linarith [ht.2])]
    push_cast; ring_nf
  rw [hB, ← Complex.ofReal_add, Complex.Gamma_ofReal, Complex.Gamma_ofReal,
    Complex.Gamma_ofReal] at h
  have h' : Gamma α * Gamma β = Gamma (α + β) *
      ∫ t in (0:ℝ)..1, t ^ (α - 1) * (1 - t) ^ (β - 1) := by exact_mod_cast h
  have hg := Gamma_pos_of_pos (add_pos ha hb)
  field_simp
  linarith

/-- the map `t ↦ t/(1-t)` sends `(0,1)` onto `(0,∞)` -/
lemma image_logit : (fun t : ℝ => t / (1 - t)) '' Ioo 0 1 = Ioi 0 := by
  ext s; constructor
  · rintro ⟨t, ⟨h0, h1⟩, rfl⟩
    exact div_pos h0 (sub_pos.mpr h1)
  · intro hs
    have hs : 0 < s := hs
    refine ⟨s / (1 + s), ⟨by positivity, by rw [div_lt_one (by linarith)]; linarith⟩, ?_⟩
    simp only
    have : 1 - s / (1 + s) = 1 / (1 + s) := by field_simp; ring
    rw [this]; field_simp

lemma deriv_logit (t : ℝ) (ht : t ∈ Ioo (0:ℝ) 1) :
    HasDerivWithinAt (fun t : ℝ => t / (1 - t)) (1 / (1 - t) ^ 2) (Ioo 0 1) t := by
  have h1 : (1 - t) ≠ 0 := by linarith [ht.2]
  have := ((hasDerivAt_id' t).div ((hasDerivAt_id' t).const_sub 1) h1)
  convert this.hasDerivWithinAt using 1
  field_simp; ring

lemma injOn_logit : InjOn (fun t : ℝ => t / (1 - t)) (Ioo 0 1) := by
  intro x hx y hy hxy
  simp only at hxy
  have h1 : (1 - x) ≠ 0 := by linarith [hx.2]
  have h2 : (1 - y) ≠ 0 := by linarith [hy.2]
  rw [div_eq_div_iff h1 h2] at hxy
  linarith

lemma pull_logit (a b t : ℝ) (ht : t ∈ Ioo (0:ℝ) 1) :
    |1 / (1 - t) ^ 2| • ((t / (1 - t)) ^ (a - 1) * (1 + t / (1 - t)) ^ (-(a + b))) =
      t ^ (a - 1) * (1 - t) ^ (b - 1) := by
  have h0 : 0 < t := ht.1
  have h1 : 0 < 1 - t := by linarith [ht.2]
  have e : 1 + t / (1 - t) = (1 - t)⁻¹ := by field_simp; ring
  rw [e, smul_eq_mul, abs_of_pos (by positivity), Real.div_rpow h0.le h1.le,
    Real.inv_rpow h1.le, ← Real.rpow_neg h1.le, neg_neg]
  rw [show (1:ℝ) / (1 - t) ^ 2 = (1 - t) ^ (-2:ℝ) by
    rw [Real.rpow_neg h1.le, one_div]; norm_cast]
  rw [div_eq_mul_inv, ← Real.rpow_neg h1.le]
  have : (1 - t) ^ (-2:ℝ) * (t ^ (a - 1) * (1 - t) ^ (-(a - 1)) * (1 - t) ^ (a + b)) =
      t ^ (a - 1) * ((1 - t) ^ (-2:ℝ) * (1 - t) ^ (-(a - 1)) * (1 - t) ^ (a + b)) := by ring
  rw [this, ← Real.rpow_add h1, ← Real.rpow_add h1]
  congr 2; ring

/-- Beta integral on `(0,∞)`. -/
lemma beta_Ioi {a b : ℝ} (ha : 0 < a) (hb : 0 < b) :
    IntegrableOn (fun s : ℝ => s ^ (a - 1) * (1 + s) ^ (-(a + b))) (Ioi 0) ∧
    ∫ s in Ioi (0:ℝ), s ^ (a - 1) * (1 + s) ^ (-(a + b)) = Gamma a * Gamma b / Gamma (a + b) := by
  rw [← image_logit]
  constructor
  · rw [integrableOn_image_iff_integrableOn_abs_deriv_smul measurableSet_Ioo deriv_logit
      injOn_logit]
    refine (integrableOn_congr_fun (fun t ht => (pull_logit a b t ht).symm) measurableSet_Ioo).mp ?_
    have hc := Complex.betaIntegral_convergent (u := a) (v := b) (by simpa using ha)
      (by simpa using hb)
    rw [← integrableOn_Icc_iff_integrableOn_Ioo, ← intervalIntegrable_iff_integrableOn_Icc_of_le zero_le_one]
    refine (hc.norm).congr ?_
    rw [uIoc_of_le zero_le_one]
    intro t ht
    have h1 : 0 ≤ 1 - t := by linarith [ht.2]
    simp only
    have e : (t:ℂ) ^ ((a:ℂ) - 1) * (1 - (t:ℂ)) ^ ((b:ℂ) - 1) =
        ((t ^ (a - 1) * (1 - t) ^ (b - 1) : ℝ) : ℂ) := by
      rw [Complex.ofReal_mul, Complex.ofReal_cpow ht.1.le, Complex.ofReal_cpow h1]
      push_cast; ring_nf
    rw [e, Complex.norm_real, Real.norm_of_nonneg (mul_nonneg (Real.rpow_nonneg ht.1.le _) (Real.rpow_nonneg h1 _))]
  · rw [integral_image_eq_integral_abs_deriv_smul measurableSet_Ioo deriv_logit injOn_logit,
      setIntegral_congr_fun measurableSet_Ioo (fun t ht => pull_logit a b t ht),
      ← integral_Ioc_eq_integral_Ioo, ← intervalIntegral.integral_of_le zero_le_one]
    exact beta_real ha hb

lemma choose_succ' (a : ℝ) (n : ℕ) :
    Ring.choose a (n + 1) = Ring.choose a n * (a - n) / (n + 1) := by
  rw [Ring.choose_eq_smul, Ring.choose_eq_smul, descPochhammer_succ_right, smeval_mul,
    smeval_sub, smeval_X, smeval_natCast, Nat.factorial_succ]
  simp only [smul_eq_mul, npow_one, npow_zero, nsmul_eq_mul, mul_one]
  push_cast
  field_simp

lemma choose_half (j : ℕ) :
    Ring.choose (1 / 2 : ℝ) (j + 2) * (-4) ^ (j + 2) =
      -2 * ((2 * j + 2).factorial : ℝ) / ((j + 1).factorial * (j + 2).factorial) := by
  induction j with
  | zero =>
    rw [choose_succ', choose_succ', Ring.choose_zero_right]
    norm_num [Nat.factorial]
  | succ j ih =>
    rw [show j + 1 + 2 = (j + 2) + 1 by ring, choose_succ', pow_succ]
    have e : Ring.choose (1 / 2 : ℝ) (j + 2) * (1 / 2 - ((j + 2 : ℕ) : ℝ)) / (((j + 2 : ℕ) : ℝ) + 1)
        * ((-4) ^ (j + 2) * (-4)) =
        (Ring.choose (1 / 2 : ℝ) (j + 2) * (-4) ^ (j + 2)) *
          ((1 / 2 - ((j + 2 : ℕ) : ℝ)) * (-4) / (((j + 2 : ℕ) : ℝ) + 1)) := by ring
    rw [e, ih, show 2 * (j + 1) + 2 = (2 * j + 2) + 1 + 1 by ring, show j + 1 + 1 = j + 2 by ring,
      show j + 2 + 1 = (j + 2) + 1 by ring, Nat.factorial_succ (2 * j + 2 + 1),
      Nat.factorial_succ (2 * j + 2), Nat.factorial_succ (j + 2)]
    simp only [Nat.factorial_succ]
    push_cast
    field_simp
    ring

lemma Gamma_add_nat (x : ℝ) (hx : 0 < x) (n : ℕ) :
    Gamma (x + n) = Gamma x * (ascPochhammer ℝ n).eval x := by
  induction n with
  | zero => simp
  | succ n ih =>
    rw [T13.pe_succ_right', Nat.cast_succ, ← add_assoc, Real.Gamma_add_one (by positivity), ih]
    ring

lemma term_int (z ε : ℝ) (hε0 : 0 < ε) (hε1 : ε < 1) (j : ℕ) :
    (Ring.choose (1 / 2 : ℝ) (j + 2) * (-4) ^ (j + 2) * z ^ (j + 2)) *
      (Gamma ((j:ℝ) + 2 + ε) * Gamma ((j:ℝ) + 1 - ε) /
        Gamma (((j:ℝ) + 2 + ε) + ((j:ℝ) + 1 - ε))) =
      π / Real.sin (π * ε) * T13.lamTerm ε z (j + 1) := by
  have g1 : Gamma ((j:ℝ) + 2 + ε) = ε * Gamma ε * (ascPochhammer ℝ (j + 1)).eval (1 + ε) := by
    rw [show (j:ℝ) + 2 + ε = (1 + ε) + ((j + 1 : ℕ) : ℝ) by push_cast; ring,
      Gamma_add_nat _ (by linarith), add_comm 1 ε, Real.Gamma_add_one hε0.ne']
  have g2 : Gamma ((j:ℝ) + 1 - ε) = Gamma (1 - ε) * (ascPochhammer ℝ j).eval (1 - ε) := by
    rw [show (j:ℝ) + 1 - ε = (1 - ε) + (j : ℝ) by ring, Gamma_add_nat _ (by linarith)]
  have g3 : Gamma (((j:ℝ) + 2 + ε) + ((j:ℝ) + 1 - ε)) = ((2 * j + 2).factorial : ℝ) := by
    rw [show ((j:ℝ) + 2 + ε) + ((j:ℝ) + 1 - ε) = ((2 * j + 2 : ℕ) : ℝ) + 1 by push_cast; ring,
      Real.Gamma_nat_eq_factorial]
  have hrefl := Real.Gamma_mul_Gamma_one_sub ε
  have p1 : (ascPochhammer ℝ (j + 1)).eval (-ε) = (-ε) * (ascPochhammer ℝ j).eval (1 - ε) := by
    rw [ascPochhammer_succ_left, eval_mul, eval_X, eval_comp, eval_add, eval_X, eval_one]
    ring_nf
  rw [g1, g2, g3, choose_half]
  unfold T13.lamTerm
  rw [p1, T13.pe_two, ← hrefl]
  have hf1 : ((j + 1).factorial : ℝ) ≠ 0 := by positivity
  have hf2 : ((j + 1 + 1).factorial : ℝ) ≠ 0 := by positivity
  have hf3 : ((2 * j + 2).factorial : ℝ) ≠ 0 := by positivity
  rw [show j + 1 + 1 = j + 2 by ring] at hf2 ⊢
  field_simp
  ring

lemma lamTerm_zero (ν z : ℝ) : T13.lamTerm ν z 0 = 2 * z := by
  simp [T13.lamTerm]

lemma sqrt_rem_bound {x : ℝ} (h0 : -1 ≤ x) (h1 : x ≤ 0) :
    |Real.sqrt (1 + x) - 1 - x / 2| ≤ x ^ 2 / 2 := by
  set y := Real.sqrt (1 + x) with hy
  have hy0 : 0 ≤ y := Real.sqrt_nonneg _
  have hyy : y ^ 2 = 1 + x := Real.sq_sqrt (by linarith)
  have hy1 : y ≤ 1 := by nlinarith
  have e : y - 1 - x / 2 = -(y - 1) ^ 2 / 2 := by
    have : x = y ^ 2 - 1 := by linarith
    rw [this]; ring
  rw [e, abs_div, abs_neg, abs_of_nonneg (sq_nonneg _), abs_two]
  have : x = y ^ 2 - 1 := by linarith
  rw [this]
  have : (y ^ 2 - 1) ^ 2 = (1 - y) ^ 2 * (1 + y) ^ 2 := by ring
  rw [this]
  have : (1:ℝ) ≤ (1 + y) ^ 2 := by nlinarith
  nlinarith [sq_nonneg (1 - y)]

/-- the Mellin identity with `c = 1 - 2z`. -/
lemma mellin (z ε : ℝ) (hz0 : 0 ≤ z) (hz1 : z < 1) (hε0 : 0 < ε) (hε1 : ε < 1) :
    IntegrableOn (fun s => k (1 - 2 * z) s * s ^ (ε - 1)) (Ioi 0) ∧
    ∫ s in Ioi 0, k (1 - 2 * z) s * s ^ (ε - 1) =
      π / Real.sin (π * ε) * ∑' j : ℕ, T13.lamTerm ε z j := by
  set C : ℕ → ℝ := fun j => Ring.choose (1 / 2 : ℝ) (j + 2) * (-4) ^ (j + 2) * z ^ (j + 2) with hC
  set φ : ℕ → ℝ → ℝ := fun j s =>
    s ^ (((j:ℝ) + 2 + ε) - 1) * (1 + s) ^ (-(((j:ℝ) + 2 + ε) + ((j:ℝ) + 1 - ε))) with hφ
  set F0 : ℝ → ℝ := fun s => 2 * z * (s ^ (ε - 1) * (1 + s) ^ (-(ε + (1 - ε)))) with hF0
  set Rm : ℝ → ℝ := fun s => k (1 - 2 * z) s * s ^ (ε - 1) - F0 s with hRm
  have hB0 := beta_Ioi hε0 (by linarith : 0 < 1 - ε)
  have hBj : ∀ j : ℕ, IntegrableOn (φ j) (Ioi 0) ∧
      ∫ s in Ioi (0:ℝ), φ j s = Gamma ((j:ℝ) + 2 + ε) * Gamma ((j:ℝ) + 1 - ε) /
        Gamma (((j:ℝ) + 2 + ε) + ((j:ℝ) + 1 - ε)) := fun j =>
    beta_Ioi (by positivity) (by have := (Nat.cast_nonneg j : (0:ℝ) ≤ j); linarith)
  have hser : ∀ s ∈ Ioi (0:ℝ), HasSum (fun j => C j * φ j s) (Rm s) := fun s hs =>
    k_series z ε s hz0 hz1 hs
  have hφ0 : ∀ j, ∀ s ∈ Ioi (0:ℝ), 0 ≤ φ j s := fun j s hs => by
    have hs : (0:ℝ) < s := hs
    exact mul_nonneg (Real.rpow_nonneg hs.le _) (Real.rpow_nonneg (by linarith) _)
  -- integrals of the terms
  have hGint : ∀ j, ∫ s in Ioi (0:ℝ), C j * φ j s = π / Real.sin (π * ε) * T13.lamTerm ε z (j + 1) := by
    intro j
    rw [integral_const_mul, (hBj j).2]
    exact term_int z ε hε0 hε1 j
  have hGnorm : ∀ j, ∫ s in Ioi (0:ℝ), ‖C j * φ j s‖ = |π / Real.sin (π * ε) * T13.lamTerm ε z (j + 1)| := by
    intro j
    rw [← hGint j, setIntegral_congr_fun measurableSet_Ioi (fun s hs => by
      rw [norm_mul, Real.norm_eq_abs, Real.norm_eq_abs, abs_of_nonneg (hφ0 j s hs)]),
      integral_const_mul, integral_const_mul, abs_mul (C j), abs_of_nonneg
        (setIntegral_nonneg measurableSet_Ioi (hφ0 j))]
  have hsumL := Y02.lamTerm_summable ε z hz0 hz1
  have hsum1 : Summable (fun j => |π / Real.sin (π * ε) * T13.lamTerm ε z (j + 1)|) :=
    (((summable_nat_add_iff 1).mpr hsumL).mul_left _).abs
  have hHS := hasSum_integral_of_summable_integral_norm (μ := volume.restrict (Ioi (0:ℝ)))
    (F := fun j s => C j * φ j s) (fun j => ((hBj j).1.const_mul (C j)))
    (by simpa only [hGnorm] using hsum1)
  have hRm_eq : ∫ s in Ioi (0:ℝ), (∑' j, C j * φ j s) = ∫ s in Ioi (0:ℝ), Rm s :=
    setIntegral_congr_fun measurableSet_Ioi (fun s hs => (hser s hs).tsum_eq)
  rw [hRm_eq] at hHS
  simp only [hGint] at hHS
  -- integrability of the remainder
  have hF0int : IntegrableOn F0 (Ioi 0) := hB0.1.const_mul (2 * z)
  have hcont : ContinuousOn (fun s => k (1 - 2 * z) s * s ^ (ε - 1)) (Ioi 0) := by
    intro s hs
    have hs : (0:ℝ) < s := hs
    show ContinuousWithinAt (fun s => (Real.sqrt (1 + 2 * (1 - 2 * z) * s + s ^ 2) - s - (1 - 2 * z))
      * s ^ (ε - 1)) (Ioi 0) s
    have h1 : ContinuousAt (fun s => Real.sqrt (1 + 2 * (1 - 2 * z) * s + s ^ 2) - s - (1 - 2 * z)) s := by
      fun_prop
    have h2 : ContinuousAt (fun s : ℝ => s ^ (ε - 1)) s := (continuousAt_id).rpow_const (Or.inl hs.ne')
    exact (h1.mul h2).continuousWithinAt
  have hF0cont : ContinuousOn F0 (Ioi 0) := by
    intro s hs
    have hs : (0:ℝ) < s := hs
    have h2 : ContinuousAt (fun s : ℝ => s ^ (ε - 1)) s := (continuousAt_id).rpow_const (Or.inl hs.ne')
    have h3 : ContinuousAt (fun s : ℝ => (1 + s) ^ (-(ε + (1 - ε)))) s :=
      (continuousAt_const.add continuousAt_id).rpow_const (Or.inl (by show (1:ℝ) + s ≠ 0; linarith))
    exact ((continuousAt_const.mul (h2.mul h3))).continuousWithinAt
  have hRint : IntegrableOn Rm (Ioi 0) := by
    have hb := (beta_Ioi (a := ε + 2) (b := 1 - ε) (by linarith) (by linarith)).1.const_mul 8
    refine Integrable.mono' hb ((hcont.sub hF0cont).aestronglyMeasurable measurableSet_Ioi) ?_
    rw [ae_restrict_iff' measurableSet_Ioi]
    refine Filter.Eventually.of_forall (fun s hs => ?_)
    have hs : (0:ℝ) < s := hs
    have hs1 : 0 < 1 + s := by linarith
    set x := -(z * (4 * s / (1 + s) ^ 2)) with hxdef
    have hu0 : 0 ≤ 4 * s / (1 + s) ^ 2 := by positivity
    have hu1 : 4 * s / (1 + s) ^ 2 ≤ 1 := by
      rw [div_le_one (by positivity)]; nlinarith [sq_nonneg (1 - s)]
    have hx0 : x ≤ 0 := by rw [hxdef]; nlinarith [mul_nonneg hz0 hu0]
    have hx1 : -1 ≤ x := by rw [hxdef]; nlinarith [mul_le_mul_of_nonneg_left hu1 hz0]
    have e4 : 1 + 2 * (1 - 2 * z) * s + s ^ 2 = (1 + s) ^ 2 * (1 + x) := by
      rw [hxdef]; field_simp; ring
    have eR : Rm s = (1 + s) * s ^ (ε - 1) * (Real.sqrt (1 + x) - 1 - x / 2) := by
      simp only [hRm, hF0, k]
      rw [show -(ε + (1 - ε)) = (-1 : ℝ) by ring, Real.rpow_neg_one, e4,
        Real.sqrt_mul (sq_nonneg _), Real.sqrt_sq hs1.le, hxdef]
      field_simp
      ring
    have hbnd := sqrt_rem_bound hx1 hx0
    have hsp : 0 < s ^ (ε - 1) := Real.rpow_pos_of_pos hs _
    have hx2 : x ^ 2 ≤ (4 * s / (1 + s) ^ 2) ^ 2 := by
      rw [hxdef, neg_sq, mul_pow]
      have : z ^ 2 ≤ 1 := by nlinarith
      nlinarith [sq_nonneg (4 * s / (1 + s) ^ 2)]
    rw [Real.norm_eq_abs, eR, abs_mul, abs_of_pos (by positivity : 0 < (1 + s) * s ^ (ε - 1))]
    have e5 : 8 * (s ^ (ε + 2 - 1) * (1 + s) ^ (-(ε + 2 + (1 - ε)))) =
        (1 + s) * s ^ (ε - 1) * ((4 * s / (1 + s) ^ 2) ^ 2 / 2) := by
      rw [show ε + 2 - 1 = (ε - 1) + ((2:ℕ):ℝ) by push_cast; ring, Real.rpow_add_natCast hs.ne',
        show -(ε + 2 + (1 - ε)) = -((3:ℕ):ℝ) by push_cast; ring, Real.rpow_neg hs1.le,
        Real.rpow_natCast]
      field_simp
      ring
    rw [e5]
    gcongr
    linarith
  have hkint : IntegrableOn (fun s => k (1 - 2 * z) s * s ^ (ε - 1)) (Ioi 0) := by
    have := hRint.add hF0int
    refine this.congr_fun (fun s _ => ?_) measurableSet_Ioi
    simp [hRm]
  refine ⟨hkint, ?_⟩
  have hsplit : ∫ s in Ioi 0, k (1 - 2 * z) s * s ^ (ε - 1) =
      (∫ s in Ioi (0:ℝ), Rm s) + ∫ s in Ioi (0:ℝ), F0 s := by
    rw [← integral_add hRint hF0int]
    refine setIntegral_congr_fun measurableSet_Ioi (fun s _ => ?_)
    simp [hRm]
  have hF0v : ∫ s in Ioi (0:ℝ), F0 s = π / Real.sin (π * ε) * T13.lamTerm ε z 0 := by
    rw [hF0, integral_const_mul, hB0.2, lamTerm_zero, ← Real.Gamma_mul_Gamma_one_sub,
      show ε + (1 - ε) = 1 by ring, Real.Gamma_one]
    ring
  rw [hsplit, ← hHS.tsum_eq, hF0v, tsum_mul_left, hsumL.tsum_eq_zero_add]
  ring

end Mellin

lemma stmt4_aux : stmt_4 := by
  intro θ s h0 h1 hs
  set c := Real.cos θ with hcdef
  have hc : |c| ≤ 1 := Real.abs_cos_le_one θ
  have hb : ∀ m, |Y02.gc θ m| ≤ 1 := by
    intro m
    rcases m with _ | _ | m
    · simp [Y02.gc]
    · simpa [Y02.gc] using hc
    · simp only [Y02.gc, abs_mul, abs_pow, abs_neg, abs_one, one_pow, one_mul]
      have := (le_min_iff.mp (Y02.stmt2_aux (m + 1) (by omega) θ)).2
      refine this.trans ?_
      rw [div_le_one (by positivity)]
      push_cast
      linarith [(Nat.cast_nonneg m : (0:ℝ) ≤ m)]
  have hs0 : 0 < s := by linarith
  set t := 1 / s with htdef
  have ht : |t| < 1 := by
    rw [htdef, abs_of_pos (by positivity), div_lt_one hs0]; exact hs
  have H := Y02.ps_sqrt (Y02.gc θ) c hc hb rfl rfl (Y02.gc_rec θ) ht
  have H3 := ((hasSum_nat_add_iff' 2).mpr H).mul_left s
  convert H3 using 1
  · funext i
    simp only [Y02.gc, htdef, div_pow, one_pow]
    rw [pow_succ s i, pow_add s i 2]
    field_simp
  · have e1 : 1 + 2 * c * t + t ^ 2 = (1 + 2 * c * s + s ^ 2) / s ^ 2 := by
      rw [htdef]; field_simp; ring
    rw [e1, Real.sqrt_div' _ (sq_nonneg s), Real.sqrt_sq hs0.le]
    simp only [Finset.sum_range_succ, Finset.sum_range_zero, Y02.gc, htdef, Y02.k]
    field_simp
    ring


lemma stmt5_aux : stmt_5 := by
  intro θ ε h0 h1 hε0 hε1
  have hc : 0 < Real.cos (θ / 2) :=
    Real.cos_pos_of_mem_Ioo ⟨by linarith [Real.pi_pos], by linarith⟩
  have hz1 : Real.sin (θ / 2) ^ 2 < 1 := by nlinarith [Real.sin_sq_add_cos_sq (θ / 2)]
  have hcos : Real.cos θ = 1 - 2 * Real.sin (θ / 2) ^ 2 := by
    have := T13.sin_half_sq θ; linarith
  rw [hcos]
  exact mellin _ ε (sq_nonneg _) hz1 hε0 hε1

section Continuation
open Set

/-- `a_i(ε) = (-1)^i Λ_{i+1}/(ε-1-i)`: the coefficients of the expansion at infinity, integrated
against `s^{ε-1}` on `(1,∞)` (up to sign). -/
noncomputable def aCoef (θ ε : ℝ) (i : ℕ) : ℝ := (-1) ^ i * T13.LamInt (i + 1) θ / (ε - 1 - i)

/-- The Mellin transform of `k_c` regularized at infinity:
`A(ε) = ∫_0^1 k_c(s) s^{ε-1} ds - Σ_i a_i(ε)`. For `0<ε<1` it is `∫_0^∞ k_c s^{ε-1}`. -/
noncomputable def Areg (θ ε : ℝ) : ℝ :=
  (∫ s in Ioc (0:ℝ) 1, k (Real.cos θ) s * s ^ (ε - 1)) - ∑' i, aCoef θ ε i

lemma lam_abs_le (θ : ℝ) (i : ℕ) : |T13.LamInt (i + 1) θ| ≤ 2 / (2 * (i:ℝ) + 3) := by
  have := (le_min_iff.mp (stmt2_aux (i + 1) (by omega) θ)).2
  convert this using 2; push_cast; ring

lemma aCoef_abs_summable (θ ε : ℝ) : Summable (fun i => |aCoef θ ε i|) := by
  apply Summable.of_norm_bounded_eventually_nat (g := fun i : ℕ => 2 / ((i:ℝ) + 1) ^ 2)
  · have h := (summable_nat_add_iff 1).mpr (Real.summable_one_div_nat_pow.mpr one_lt_two)
    have h2 := h.mul_left 2
    refine h2.congr (fun i => ?_)
    push_cast; ring
  · filter_upwards [eventually_ge_atTop (Nat.ceil (2 * |ε|))] with i hi
    have hi' : 2 * |ε| ≤ (i:ℝ) := (Nat.ceil_le).mp hi
    have hab := abs_nonneg ε
    have hε := le_abs_self ε
    have hd : ((i:ℝ) + 1) / 2 ≤ |ε - 1 - i| := by
      rw [abs_sub_comm]; apply le_abs.mpr; left
      linarith
    have hd0 : 0 < ((i:ℝ) + 1) / 2 := by positivity
    rw [Real.norm_eq_abs, abs_abs, aCoef, abs_div, abs_mul, abs_pow, abs_neg, abs_one, one_pow,
      one_mul]
    calc |T13.LamInt (i + 1) θ| / |ε - 1 - i| ≤ (2 / (2 * (i:ℝ) + 3)) / (((i:ℝ) + 1) / 2) :=
          div_le_div₀ (by positivity) (lam_abs_le θ i) hd0 hd
      _ ≤ 2 / ((i:ℝ) + 1) ^ 2 := by
          have e : (2 / (2 * (i:ℝ) + 3)) / (((i:ℝ) + 1) / 2) = 4 / ((2 * i + 3) * (i + 1)) := by
            field_simp; ring
          rw [e, div_le_div_iff₀ (by positivity) (by positivity)]
          nlinarith [(Nat.cast_nonneg i : (0:ℝ) ≤ i)]

lemma aCoef_summable (θ ε : ℝ) : Summable (aCoef θ ε) :=
  (aCoef_abs_summable θ ε).of_abs

/-- termwise integration of the expansion at infinity over a set `S ⊆ (1,∞)` -/
lemma termwise (θ ε : ℝ) (h0 : 0 ≤ θ) (h1 : θ < π) {S : Set ℝ} (hS : MeasurableSet S)
    (hS1 : S ⊆ Ioi 1)
    (hint : ∀ i : ℕ, IntegrableOn (fun s : ℝ => s ^ (ε - 2 - i)) S)
    (hsum : Summable (fun i : ℕ => |T13.LamInt (i + 1) θ| * ∫ s in S, s ^ (ε - 2 - i))) :
    HasSum (fun i : ℕ => (-1) ^ i * T13.LamInt (i + 1) θ * ∫ s in S, s ^ (ε - 2 - i))
      (∫ s in S, k (Real.cos θ) s * s ^ (ε - 1)) := by
  have hpos : ∀ i : ℕ, ∀ s ∈ S, 0 ≤ s ^ (ε - 2 - i) := fun i s hs =>
    Real.rpow_nonneg (le_of_lt (lt_trans one_pos (hS1 hs))) _
  have H := hasSum_integral_of_summable_integral_norm (μ := volume.restrict S)
    (F := fun i s => (-1) ^ i * T13.LamInt (i + 1) θ * s ^ (ε - 2 - i))
    (fun i => (hint i).const_mul _) ?_
  · have e1 : ∀ i : ℕ, ∫ s in S, (-1) ^ i * T13.LamInt (i + 1) θ * s ^ (ε - 2 - i) =
        (-1) ^ i * T13.LamInt (i + 1) θ * ∫ s in S, s ^ (ε - 2 - i) :=
      fun i => integral_const_mul _ _
    simp only [e1] at H
    convert H using 1
    refine setIntegral_congr_fun hS (fun s hs => ?_)
    have hs1 : 1 < s := hS1 hs
    have hs0 : 0 < s := by linarith
    have H4 := (stmt4_aux θ s h0 h1 hs1).mul_right (s ^ (ε - 1))
    rw [← H4.tsum_eq]
    congr 1; funext i
    rw [show ε - 2 - (i:ℝ) = (ε - 1) - ((i + 1 : ℕ) : ℝ) by push_cast; ring,
      Real.rpow_sub_natCast hs0.ne']
    ring
  · convert hsum using 1
    funext i
    rw [← integral_const_mul]
    refine setIntegral_congr_fun hS (fun s hs => ?_)
    rw [norm_mul, norm_mul, norm_pow, norm_neg, norm_one, one_pow, one_mul, Real.norm_eq_abs,
      Real.norm_eq_abs, abs_of_nonneg (hpos i s hs)]

lemma int_Ioc_rpow {n p : ℝ} (hn : 1 < n) (hp : p ≠ -1) :
    ∫ s in Ioc 1 n, s ^ p = (n ^ (p + 1) - 1) / (p + 1) := by
  rw [← intervalIntegral.integral_of_le hn.le, integral_rpow (Or.inr ⟨hp, ?_⟩), Real.one_rpow]
  rw [Set.mem_uIcc]; push Not
  constructor <;> intro h <;> linarith

lemma intOn_Ioc_rpow {n : ℝ} (p : ℝ) : IntegrableOn (fun s : ℝ => s ^ p) (Ioc 1 n) := by
  refine (ContinuousOn.integrableOn_Icc ?_).mono_set Ioc_subset_Icc_self
  exact ContinuousOn.rpow_const continuousOn_id (fun x hx => Or.inl (by
    show x ≠ 0; linarith [hx.1]))

lemma nonint_ne (ε : ℝ) (hnat : ∀ j : ℕ, ε ≠ j) (i : ℕ) : ε - 1 - i ≠ 0 := by
  intro h; apply hnat (i + 1); push_cast; linarith

/-- `∫_1^n k_c s^{ε-1} = Σ_i a_i (n^{ε-1-i} - 1)`. -/
lemma int_one_n (θ ε n : ℝ) (h0 : 0 ≤ θ) (h1 : θ < π) (hn : 1 < n) (_hε : 0 < ε)
    (hnat : ∀ j : ℕ, ε ≠ j) :
    HasSum (fun i : ℕ => aCoef θ ε i * (n ^ (ε - 1 - i) - 1))
      (∫ s in Ioc 1 n, k (Real.cos θ) s * s ^ (ε - 1)) := by
  have hval : ∀ i : ℕ, ∫ s in Ioc 1 n, s ^ (ε - 2 - i) = (n ^ (ε - 1 - i) - 1) / (ε - 1 - i) := by
    intro i
    rw [int_Ioc_rpow hn (by intro h; exact nonint_ne ε hnat i (by linarith)),
      show ε - 2 - (i:ℝ) + 1 = ε - 1 - i by ring]
  have hnn : ∀ i : ℕ, 0 ≤ ∫ s in Ioc 1 n, s ^ (ε - 2 - i) := fun i =>
    setIntegral_nonneg measurableSet_Ioc (fun s hs => Real.rpow_nonneg (by linarith [hs.1]) _)
  have H := termwise θ ε h0 h1 (S := Ioc 1 n) measurableSet_Ioc (fun s hs => hs.1)
    (fun i => intOn_Ioc_rpow _) ?_
  · convert H using 1
    funext i
    rw [hval, aCoef]
    ring
  · have hb : Summable (fun i => |aCoef θ ε i| * (n ^ ε + 1)) :=
      (aCoef_abs_summable θ ε).mul_right _
    refine Summable.of_nonneg_of_le (fun i => mul_nonneg (abs_nonneg _) (hnn i)) (fun i => ?_) hb
    rw [← abs_of_nonneg (hnn i), ← abs_mul, hval, aCoef]
    have e : T13.LamInt (i + 1) θ * ((n ^ (ε - 1 - i) - 1) / (ε - 1 - i)) =
        ((-1) ^ i * T13.LamInt (i + 1) θ / (ε - 1 - i)) * ((-1) ^ i * (n ^ (ε - 1 - i) - 1)) := by
      have : ((-1 : ℝ) ^ i) ^ 2 = 1 := by rw [← pow_mul, mul_comm, pow_mul]; simp
      field_simp [nonint_ne ε hnat i]
      linear_combination (-(T13.LamInt (i + 1) θ * (n ^ (ε - 1 - i) - 1))) * this
    rw [e, abs_mul]
    gcongr
    rw [abs_mul, abs_pow, abs_neg, abs_one, one_pow, one_mul]
    have hq : n ^ (ε - 1 - i) ≤ n ^ ε :=
      Real.rpow_le_rpow_of_exponent_le hn.le (by linarith [(Nat.cast_nonneg i : (0:ℝ) ≤ i)])
    have hq0 : 0 < n ^ (ε - 1 - i) := Real.rpow_pos_of_pos (by linarith) _
    rw [abs_le]; constructor <;> linarith

/-- `∫_0^n k_c s^{ε-1} = A(ε) + Σ_i a_i n^{ε-1-i}` for `n > 1`, `ε > 0`, `ε ∉ ℕ`. -/
lemma int_zero_n (θ ε n : ℝ) (h0 : 0 ≤ θ) (h1 : θ < π) (hn : 1 < n) (hε : 0 < ε)
    (hnat : ∀ j : ℕ, ε ≠ j) :
    IntegrableOn (fun s => k (Real.cos θ) s * s ^ (ε - 1)) (Ioc 0 n) ∧
    ∫ s in Ioc 0 n, k (Real.cos θ) s * s ^ (ε - 1) =
      Areg θ ε + ∑' i, aCoef θ ε i * n ^ (ε - 1 - i) := by
  have hint : IntegrableOn (fun s => k (Real.cos θ) s * s ^ (ε - 1)) (Ioc 0 n) := by
    have hr : IntegrableOn (fun s : ℝ => s ^ (ε - 1)) (Ioc 0 n) := by
      have := intervalIntegral.intervalIntegrable_rpow' (a := 0) (b := n) (r := ε - 1) (by linarith)
      exact (intervalIntegrable_iff_integrableOn_Ioc_of_le (by linarith)).mp this
    have hkc : ContinuousOn (fun s => k (Real.cos θ) s) (Icc 0 n) := by
      apply Continuous.continuousOn; unfold k; fun_prop
    exact hr.continuousOn_mul_of_subset hkc isCompact_Icc measurableSet_Ioc Ioc_subset_Icc_self
  refine ⟨hint, ?_⟩
  have hS1 : Summable (fun i => aCoef θ ε i * n ^ (ε - 1 - i)) := by
    refine Summable.of_norm_bounded ((aCoef_abs_summable θ ε).mul_right (n ^ ε)) (fun i => ?_)
    rw [Real.norm_eq_abs, abs_mul, abs_of_pos (Real.rpow_pos_of_pos (by linarith) _)]
    exact mul_le_mul_of_nonneg_left (Real.rpow_le_rpow_of_exponent_le hn.le
      (by linarith [(Nat.cast_nonneg i : (0:ℝ) ≤ i)])) (abs_nonneg _)
  rw [← Ioc_union_Ioc_eq_Ioc zero_le_one hn.le,
    setIntegral_union (Set.Ioc_disjoint_Ioc_of_le le_rfl) measurableSet_Ioc
      (hint.mono_set (Ioc_subset_Ioc_right hn.le)) (hint.mono_set (Ioc_subset_Ioc_left zero_le_one)),
    (int_one_n θ ε n h0 h1 hn hε hnat).tsum_eq.symm, Areg]
  have : ∑' i, aCoef θ ε i * (n ^ (ε - 1 - i) - 1) =
      ∑' i, aCoef θ ε i * n ^ (ε - 1 - i) - ∑' i, aCoef θ ε i := by
    rw [← (hS1.tsum_sub (aCoef_summable θ ε))]
    congr 1; funext i; ring
  rw [this]; ring

/-- Reduction of stmt_6 to the identity `A(ε) = (π/sin πε) Λ_ε` for all `ε > 0`, `ε ∉ ℕ`. -/
lemma stmt6_of_Areg
    (hA : ∀ θ ε : ℝ, 0 ≤ θ → θ < π → 0 < ε → (∀ j : ℕ, ε ≠ j) →
      Areg θ ε = π / Real.sin (π * ε) * Lam ε θ) : stmt_6 := by
  intro θ n h0 h1 hn hnat
  set ε := eps n with hεdef
  have hε : 0 < ε := by rw [hεdef, eps]; exact div_pos one_pos (by linarith)
  have hn0 : 0 < n := by linarith
  obtain ⟨hint, hI⟩ := int_zero_n θ ε n h0 h1 hn hε hnat
  obtain ⟨hhint, hdef⟩ := defect (Real.cos θ) n hn
  rw [← L_eq θ n hn, ← hεdef] at hdef
  have hpow : ∫ s in Ioc 0 n, s ^ (ε - 1) = n ^ ε / ε := by
    rw [← intervalIntegral.integral_of_le hn0.le, integral_rpow (Or.inl (by linarith)),
      show ε - 1 + 1 = ε by ring, Real.zero_rpow hε.ne']
    ring
  have hr : IntegrableOn (fun s : ℝ => s ^ (ε - 1)) (Ioc 0 n) := by
    have := intervalIntegral.intervalIntegrable_rpow' (a := 0) (b := n) (r := ε - 1) (by linarith)
    exact (intervalIntegrable_iff_integrableOn_Ioc_of_le hn0.le).mp this
  have hhk : ∫ s in Ioc 0 n, h (Real.cos θ) s * s ^ (ε - 1) =
      (1 - Real.cos θ) * (n ^ ε / ε) - ∫ s in Ioc 0 n, k (Real.cos θ) s * s ^ (ε - 1) := by
    rw [← hpow, ← integral_const_mul, ← integral_sub (hr.const_mul _) hint]
    refine setIntegral_congr_fun measurableSet_Ioc (fun s _ => ?_)
    simp only [h, k]; ring
  have hnε : n ^ (-ε) * n ^ ε = 1 := by
    rw [← Real.rpow_add hn0, neg_add_cancel, Real.rpow_zero]
  have hL : L θ n = 1 + Real.cos θ + ε * n ^ (-ε) * ∫ s in Ioc 0 n, k (Real.cos θ) s * s ^ (ε - 1) := by
    have : 2 - L θ n = (1 - Real.cos θ) * (n ^ (-ε) * n ^ ε) -
        ε * n ^ (-ε) * ∫ s in Ioc 0 n, k (Real.cos θ) s * s ^ (ε - 1) := by
      rw [hdef, hhk]; field_simp
    rw [hnε] at this; linarith
  have hS1 : Summable (fun i => aCoef θ ε i * n ^ (ε - 1 - i)) := by
    refine Summable.of_norm_bounded ((aCoef_abs_summable θ ε).mul_right (n ^ ε)) (fun i => ?_)
    rw [Real.norm_eq_abs, abs_mul, abs_of_pos (Real.rpow_pos_of_pos (by linarith) _)]
    exact mul_le_mul_of_nonneg_left (Real.rpow_le_rpow_of_exponent_le hn.le
      (by linarith [(Nat.cast_nonneg i : (0:ℝ) ≤ i)])) (abs_nonneg _)
  have H := hS1.hasSum.mul_left (ε * n ^ (-ε))
  rw [hL, hI, hA θ ε h0 h1 hε hnat]
  convert H using 1
  · funext i
    have hne := nonint_ne ε hnat i
    have hne2 : ((i:ℝ) + 1) - ε ≠ 0 := by intro h; apply hne; linarith
    have hpw : n ^ (-ε) * n ^ (ε - 1 - i) = 1 / n ^ (i + 1) := by
      rw [← Real.rpow_add hn0, show -ε + (ε - 1 - i) = -((i + 1 : ℕ) : ℝ) by push_cast; ring,
        Real.rpow_neg hn0.le, Real.rpow_natCast, one_div]
    rw [serTerm, aCoef, ← hεdef]
    have e : ε * n ^ (-ε) * ((-1) ^ i * T13.LamInt (i + 1) θ / (ε - 1 - i) * n ^ (ε - 1 - i)) =
        ε * ((-1) ^ i * T13.LamInt (i + 1) θ / (ε - 1 - i)) * (n ^ (-ε) * n ^ (ε - 1 - i)) := by ring
    rw [e, hpw, pow_succ]
    field_simp
    ring
  · ring

/-- the identity `A(ε) = (π/sin πε) Λ_ε` for `0 < ε < 1` (from stmt_5). -/
lemma Areg_small (θ ε : ℝ) (h0 : 0 ≤ θ) (h1 : θ < π) (hε0 : 0 < ε) (hε1 : ε < 1) :
    Areg θ ε = π / Real.sin (π * ε) * Lam ε θ := by
  obtain ⟨hint, hval⟩ := stmt5_aux θ ε h0 h1 hε0 hε1
  have hnat : ∀ j : ℕ, ε ≠ j := by
    intro j hj
    rcases Nat.eq_zero_or_pos j with h | h
    · subst h; simp at hj; linarith
    · have : (1:ℝ) ≤ j := by exact_mod_cast h
      linarith
  have hlt : ∀ i : ℕ, ε - 2 - (i:ℝ) < -1 := fun i => by
    linarith [(Nat.cast_nonneg i : (0:ℝ) ≤ i)]
  have hval2 : ∀ i : ℕ, ∫ s in Ioi (1:ℝ), s ^ (ε - 2 - i) = -1 / (ε - 1 - i) := by
    intro i
    rw [integral_Ioi_rpow_of_lt (hlt i) one_pos, Real.one_rpow]
    congr 1; ring
  have hnn : ∀ i : ℕ, 0 ≤ ∫ s in Ioi (1:ℝ), s ^ (ε - 2 - i) := fun i =>
    setIntegral_nonneg measurableSet_Ioi (fun s hs => Real.rpow_nonneg (by
      have : (1:ℝ) < s := hs; linarith) _)
  have H := termwise θ ε h0 h1 measurableSet_Ioi (subset_refl _)
    (fun i => integrableOn_Ioi_rpow_of_lt (hlt i) one_pos) ?_
  · have hH : HasSum (fun i => -aCoef θ ε i) (∫ s in Ioi 1, k (Real.cos θ) s * s ^ (ε - 1)) := by
      convert H using 1
      funext i; rw [hval2, aCoef]; ring
    have hsplit : ∫ s in Ioi 0, k (Real.cos θ) s * s ^ (ε - 1) =
        (∫ s in Ioc 0 1, k (Real.cos θ) s * s ^ (ε - 1)) +
          ∫ s in Ioi 1, k (Real.cos θ) s * s ^ (ε - 1) := by
      rw [← Ioc_union_Ioi_eq_Ioi zero_le_one, setIntegral_union (Ioc_disjoint_Ioi le_rfl)
        measurableSet_Ioi (hint.mono_set Ioc_subset_Ioi_self)
        (hint.mono_set (Ioi_subset_Ioi zero_le_one))]
    rw [← hval, hsplit, ← hH.tsum_eq, tsum_neg, Areg]
    ring
  · refine (aCoef_abs_summable θ ε).congr (fun i => ?_)
    rw [hval2, aCoef, abs_div, abs_mul, abs_pow, abs_neg, abs_one, one_pow, one_mul,
      abs_of_neg (by linarith [hlt i] : ε - 1 - (i:ℝ) < 0)]
    field_simp

end Continuation

section Taylor
open Set

/-- `Π_{j<k} (1 - ε²/(j+1)²)` -/
noncomputable def pp (ε : ℝ) (k : ℕ) : ℝ :=
  ∏ j ∈ Finset.range k, (1 - ε ^ 2 / ((j:ℝ) + 1) ^ 2)

lemma pp_succ (ε : ℝ) (k : ℕ) : pp ε (k + 1) = pp ε k * (1 - ε ^ 2 / ((k:ℝ) + 1) ^ 2) := by
  unfold pp; rw [Finset.prod_range_succ]

/-- g3:lem:g(a) in the indexing of `T13.lamTerm`. -/
lemma lamTerm_closed (ε z : ℝ) (k : ℕ) :
    T13.lamTerm ε z (k + 1) =
      2 * z * z ^ (k + 1) * (-ε * ((k:ℝ) + 1 + ε) / (((k:ℝ) + 1) ^ 2 * ((k:ℝ) + 2))) * pp ε k := by
  induction k with
  | zero =>
    rw [lamTerm_succ, lamTerm_zero]; simp [pp]; ring
  | succ k ih =>
    rw [lamTerm_succ, ih, pp_succ]
    push_cast
    field_simp
    ring

lemma pp_bounds {ε : ℝ} (hε : ε ^ 2 ≤ 1) (k : ℕ) :
    0 ≤ pp ε k ∧ pp ε k ≤ 1 ∧ 1 - pp ε k ≤ ε ^ 2 * (2 - 2 / ((k:ℝ) + 1)) := by
  induction k with
  | zero => simp [pp]
  | succ k ih =>
    obtain ⟨h0, h1, h2⟩ := ih
    have hk : (1:ℝ) ≤ (k:ℝ) + 1 := by linarith [(Nat.cast_nonneg k : (0:ℝ) ≤ k)]
    have hq0 : 0 ≤ ε ^ 2 / ((k:ℝ) + 1) ^ 2 := by positivity
    have hq1 : ε ^ 2 / ((k:ℝ) + 1) ^ 2 ≤ ε ^ 2 := by
      apply div_le_self (sq_nonneg _); nlinarith
    rw [pp_succ]
    refine ⟨mul_nonneg h0 (by linarith), by nlinarith, ?_⟩
    have key : 1 / ((k:ℝ) + 1) ^ 2 ≤ 2 / ((k:ℝ) + 1) - 2 / ((k:ℝ) + 1 + 1) := by
      rw [div_sub_div _ _ (by positivity) (by positivity),
        div_le_div_iff₀ (by positivity) (by positivity)]
      nlinarith
    have e : ε ^ 2 / ((k:ℝ) + 1) ^ 2 = ε ^ 2 * (1 / ((k:ℝ) + 1) ^ 2) := by ring
    push_cast
    rw [e] at hq0 hq1 ⊢
    have : pp ε k * (ε ^ 2 * (1 / ((k:ℝ) + 1) ^ 2)) ≤ ε ^ 2 * (1 / ((k:ℝ) + 1) ^ 2) := by
      nlinarith
    nlinarith [sq_nonneg ε]

lemma sumA {z : ℝ} (h0 : 0 < z) (h1 : z < 1) :
    HasSum (fun k : ℕ => z ^ (k + 1) / (((k:ℝ) + 1) * ((k:ℝ) + 2)))
      (1 + (1 - z) / z * Real.log (1 - z)) := by
  have hz0 : z ≠ 0 := h0.ne'
  have hz : |z| < 1 := by rw [abs_lt]; constructor <;> linarith
  have hlog := Real.hasSum_pow_div_log_of_abs_lt_one hz
  have h2 : HasSum (fun n : ℕ => z ^ (n + 2) / ((n:ℝ) + 2)) (-Real.log (1 - z) - z) := by
    have := (hasSum_nat_add_iff' 1).mpr hlog
    convert this using 1
    · funext n; push_cast; ring_nf
    · simp
  have h3 := hlog.sub (h2.mul_left (1 / z))
  convert h3 using 1
  · funext k
    field_simp
    ring
  · field_simp
    ring

lemma sumB {z : ℝ} (h0 : 0 < z) (h1 : z < 1) :
    HasSum (fun k : ℕ => z ^ (k + 1) / (((k:ℝ) + 1) ^ 2 * ((k:ℝ) + 2)))
      (D00.Li2 z - 1 - (1 - z) / z * Real.log (1 - z)) := by
  have hz0 : z ≠ 0 := h0.ne'
  have hz : |z| < 1 := by rw [abs_lt]; constructor <;> linarith
  have hlog := Real.hasSum_pow_div_log_of_abs_lt_one hz
  have hL := D00.proof_6 z hz.le
  have h2 : HasSum (fun n : ℕ => z ^ (n + 2) / ((n:ℝ) + 2)) (-Real.log (1 - z) - z) := by
    have := (hasSum_nat_add_iff' 1).mpr hlog
    convert this using 1
    · funext n; push_cast; ring_nf
    · simp
  have h3 := (hL.sub hlog).add (h2.mul_left (1 / z))
  convert h3 using 1
  · funext k
    field_simp
    ring
  · field_simp
    ring

/-- `Λ_ε = 2z(1 - aε + (a - Li₂ z)ε² + r)`, `a = 1 + U`, `0 ≤ r ≤ (9/4)ε³`. -/
lemma Lam_expand (θ ε : ℝ) (h0 : 0 < θ) (h1 : θ < π) (hε0 : 0 < ε) (hε1 : ε ≤ 1 / 8) :
    ∃ r : ℝ, 0 ≤ r ∧ r ≤ 9 / 4 * ε ^ 3 ∧
      Lam ε θ = 2 * zz θ * (1 - (1 + U θ) * ε + ((1 + U θ) - D00.Li2 (zz θ)) * ε ^ 2 + r) := by
  obtain ⟨z0, z1⟩ := zz_bounds θ h0 h1
  have hzz : zz θ = Real.sin (θ / 2) ^ 2 := rfl
  set z := zz θ with hzdef
  have hU : U θ = (1 - z) / z * Real.log (1 - z) := rfl
  have hUb : -1 < U θ ∧ U θ < 0 := U_bounds z0 z1
  have hA : HasSum (fun k : ℕ => z ^ (k + 1) / (((k:ℝ) + 1) * ((k:ℝ) + 2))) (1 + U θ) := sumA z0 z1
  have hB : HasSum (fun k : ℕ => z ^ (k + 1) / (((k:ℝ) + 1) ^ 2 * ((k:ℝ) + 2)))
      (D00.Li2 z - 1 - U θ) := sumB z0 z1
  set fA := fun k : ℕ => z ^ (k + 1) / (((k:ℝ) + 1) * ((k:ℝ) + 2)) with hfA
  set fB := fun k : ℕ => z ^ (k + 1) / (((k:ℝ) + 1) ^ 2 * ((k:ℝ) + 2)) with hfB
  set rk := fun k : ℕ => z ^ (k + 1) * (ε * ((k:ℝ) + 1 + ε) / (((k:ℝ) + 1) ^ 2 * ((k:ℝ) + 2)))
      * (1 - pp ε k) with hrk
  have zpos : 0 < z := z0
  have hz0 : z ≠ 0 := z0.ne'
  have hsplit : ∀ k : ℕ, T13.lamTerm ε z (k + 1) = 2 * z * (-ε * fA k - ε ^ 2 * fB k + rk k) := by
    intro k
    rw [lamTerm_closed]
    simp only [hfA, hfB, hrk]
    field_simp
    ring
  have hsum : Summable (fun k : ℕ => T13.lamTerm ε z (k + 1)) :=
    (summable_nat_add_iff 1).mpr (lamTerm_summable ε z z0.le z1)
  have hrs : Summable rk := by
    have : Summable (fun k => T13.lamTerm ε z (k + 1) / (2 * z) + ε * fA k + ε ^ 2 * fB k) :=
      ((hsum.div_const _).add (hA.summable.mul_left ε)).add (hB.summable.mul_left (ε ^ 2))
    refine this.congr (fun k => ?_)
    rw [hsplit]; field_simp; ring
  have hpp : ∀ k, 0 ≤ pp ε k ∧ pp ε k ≤ 1 ∧ 1 - pp ε k ≤ ε ^ 2 * (2 - 2 / ((k:ℝ) + 1)) :=
    pp_bounds (by nlinarith)
  have hrk0 : ∀ k, 0 ≤ rk k := fun k => by
    simp only [hrk]
    have := (hpp k).2.1
    apply mul_nonneg (mul_nonneg (pow_nonneg zpos.le _)
      (div_nonneg (mul_nonneg hε0.le (by positivity)) (by positivity))); linarith
  have hrk1 : ∀ k, rk k ≤ 9 / 4 * ε ^ 3 * fA k := fun k => by
    simp only [hrk, hfA]
    obtain ⟨_, h2, h3⟩ := hpp k
    have hk : (1:ℝ) ≤ (k:ℝ) + 1 := by linarith [(Nat.cast_nonneg k : (0:ℝ) ≤ k)]
    have h3' : 1 - pp ε k ≤ 2 * ε ^ 2 := by
      have : 0 ≤ ε ^ 2 * (2 / ((k:ℝ) + 1)) := by positivity
      nlinarith
    have hfrac : ε * ((k:ℝ) + 1 + ε) / (((k:ℝ) + 1) ^ 2 * ((k:ℝ) + 2)) ≤
        9 / 8 * ε / (((k:ℝ) + 1) * ((k:ℝ) + 2)) := by
      rw [div_le_div_iff₀ (by positivity) (by positivity)]
      have : 0 ≤ ε * ((k:ℝ) + 1) * ((k:ℝ) + 2) := by positivity
      nlinarith
    have hzk : 0 ≤ z ^ (k + 1) := pow_nonneg zpos.le _
    calc z ^ (k + 1) * (ε * ((k:ℝ) + 1 + ε) / (((k:ℝ) + 1) ^ 2 * ((k:ℝ) + 2))) * (1 - pp ε k)
        ≤ z ^ (k + 1) * (9 / 8 * ε / (((k:ℝ) + 1) * ((k:ℝ) + 2))) * (2 * ε ^ 2) := by
          apply mul_le_mul (mul_le_mul_of_nonneg_left hfrac hzk) h3' (by linarith [h2])
          positivity
      _ = 9 / 4 * ε ^ 3 * (z ^ (k + 1) / (((k:ℝ) + 1) * ((k:ℝ) + 2))) := by
          field_simp; ring
  refine ⟨∑' k, rk k, tsum_nonneg hrk0, ?_, ?_⟩
  · have := hasSum_le hrk1 hrs.hasSum (hA.mul_left (9 / 4 * ε ^ 3))
    have hε3 : 0 ≤ ε ^ 3 := by positivity
    nlinarith
  · have H : HasSum (fun k : ℕ => T13.lamTerm ε z (k + 1))
        (2 * z * (-ε * (1 + U θ) - ε ^ 2 * (D00.Li2 z - 1 - U θ) + ∑' k, rk k)) := by
      have := (((hA.mul_left (-ε)).sub (hB.mul_left (ε ^ 2))).add hrs.hasSum).mul_left (2 * z)
      convert this using 1
      funext k; rw [hsplit]
    have H2 := (hasSum_nat_add_iff (f := fun k => T13.lamTerm ε z k) 1).mp H
    unfold Lam
    rw [← hzz, H2.tsum_eq]
    simp only [Finset.sum_range_one, lamTerm_zero]
    ring

/-- `|πε/sin(πε) - (1 + π²ε²/6)| ≤ 3ε³` for `0 < ε ≤ 1/8`. -/
lemma sinc_bound {ε : ℝ} (hε0 : 0 < ε) (hε1 : ε ≤ 1 / 8) :
    |π * ε / Real.sin (π * ε) - (1 + π ^ 2 / 6 * ε ^ 2)| ≤ 3 * ε ^ 3 := by
  set x := π * ε with hx
  have hpi := Real.pi_gt_three
  have hpi' := Real.pi_lt_d2
  have hx0 : 0 < x := by positivity
  have hx1 : x ≤ 2 / 5 := by rw [hx]; nlinarith
  have hb := Real.sin_bound (x := x) (by rw [abs_of_pos hx0]; linarith)
  rw [abs_of_pos hx0] at hb
  set sx := Real.sin x with hsx
  set δ := sx - (x - x ^ 3 / 6) with hδ
  have hδ' : |δ| ≤ x ^ 5 / 100 := hb
  have hδ2 := abs_le.mp hδ'
  have hx5 : x ^ 5 ≤ x ^ 3 * (4 / 25) := by
    have : x ^ 2 ≤ 4 / 25 := by nlinarith
    have : x ^ 5 = x ^ 3 * x ^ 2 := by ring
    rw [this]; exact mul_le_mul_of_nonneg_left (by assumption) (by positivity)
  have hsx0 : 9 / 10 * x ≤ sx := by
    have : x ^ 3 ≤ x * (4 / 25) := by nlinarith
    nlinarith
  have hsxp : 0 < sx := by linarith
  have e : π ^ 2 / 6 * ε ^ 2 = x ^ 2 / 6 := by rw [hx]; ring
  rw [e]
  have e2 : x / sx - (1 + x ^ 2 / 6) = (x ^ 5 / 36 - δ * (1 + x ^ 2 / 6)) / sx := by
    field_simp
    rw [hδ]; ring
  rw [e2, abs_div, abs_of_pos hsxp, div_le_iff₀ hsxp]
  have hnum : |x ^ 5 / 36 - δ * (1 + x ^ 2 / 6)| ≤ x ^ 5 / 36 + x ^ 5 / 100 * (1 + x ^ 2 / 6) := by
    have h1 : |δ * (1 + x ^ 2 / 6)| ≤ x ^ 5 / 100 * (1 + x ^ 2 / 6) := by
      rw [abs_mul, abs_of_pos (by positivity : (0:ℝ) < 1 + x ^ 2 / 6)]
      exact mul_le_mul_of_nonneg_right hδ' (by positivity)
    calc _ ≤ |x ^ 5 / 36| + |δ * (1 + x ^ 2 / 6)| := abs_sub _ _
      _ ≤ _ := by rw [abs_of_pos (by positivity)]; linarith
  refine hnum.trans ?_
  -- x^5 (1/36 + 1.03/100) ≤ 3 ε³ · 0.9 x, i.e. x^4 · 0.04 ≤ 2.7 ε³ with x = πε
  have hx4 : x ^ 4 ≤ x * (ε ^ 3 * 32) := by
    have : x ^ 3 ≤ ε ^ 3 * 32 := by
      rw [hx, mul_pow]
      have : π ^ 3 ≤ 32 := by nlinarith
      nlinarith [pow_pos hε0 3]
    have : x ^ 4 = x * x ^ 3 := by ring
    rw [this]; exact mul_le_mul_of_nonneg_left (by assumption) hx0.le
  have hx2 : x ^ 2 ≤ 4 / 25 := by nlinarith
  have : x ^ 5 / 36 + x ^ 5 / 100 * (1 + x ^ 2 / 6) ≤ x ^ 4 * (2 / 5) * (1 / 36 + 2 / 100) := by
    have h5 : x ^ 5 ≤ x ^ 4 * (2 / 5) := by
      have : x ^ 5 = x ^ 4 * x := by ring
      rw [this]; exact mul_le_mul_of_nonneg_left hx1 (by positivity)
    have : x ^ 5 / 100 * (1 + x ^ 2 / 6) ≤ x ^ 5 / 100 * 2 := by
      apply mul_le_mul_of_nonneg_left _ (by positivity); linarith
    nlinarith [pow_pos hx0 5, pow_pos hx0 4]
  nlinarith [pow_pos hε0 3, mul_pos hx0 (pow_pos hε0 3)]

/-- `|e^y - (1 + y + y²/2)| ≤ |y|³` for `|y| ≤ 1`. -/
lemma exp_bound3 {y : ℝ} (hy : |y| ≤ 1) : |Real.exp y - (1 + y + y ^ 2 / 2)| ≤ |y| ^ 3 := by
  have := Real.exp_bound hy (n := 3) (by norm_num)
  have e : ∑ m ∈ Finset.range 3, y ^ m / (m.factorial : ℝ) = 1 + y + y ^ 2 / 2 := by
    simp [Finset.sum_range_succ, Nat.factorial]
  rw [e] at this
  refine this.trans ?_
  have : ((Nat.succ 3 : ℕ) : ℝ) / ((Nat.factorial 3 : ℕ) * (3 : ℕ)) ≤ 1 := by
    norm_num [Nat.factorial]
  exact mul_le_of_le_one_right (by positivity) (by exact_mod_cast this)

/-- the algebra of g3:lem:taylor. -/
lemma taylor_alg (a L ε s F e : ℝ) (ha0 : 0 < a) (ha1 : a < 1) (hL0 : 0 < L) (hL1 : L < 2)
    (hε0 : 0 < ε) (hε1 : ε ≤ 1 / 8)
    (hs : |s - (1 + π ^ 2 / 6 * ε ^ 2)| ≤ 3 * ε ^ 3)
    (hF : |F - (1 - a * ε + (a - L) * ε ^ 2)| ≤ 3 * ε ^ 3)
    (he : |e - (1 + a * ε + a ^ 2 * ε ^ 2 / 2)| ≤ ε ^ 3) :
    |s * F * e - (1 + (π ^ 2 / 6 + a - a ^ 2 / 2 - L) * ε ^ 2)| ≤ 3100 * ε ^ 3 := by
  have hpi := Real.pi_gt_three
  have hpi' := Real.pi_lt_d2
  have hZ0 : 0 < π ^ 2 / 6 := by positivity
  have hZ1 : π ^ 2 / 6 < 2 := by nlinarith
  generalize π ^ 2 / 6 = Z at hZ0 hZ1 hs ⊢
  have hε2 : ε ^ 2 ≤ 1 / 64 := by nlinarith
  have hε3 : ε ^ 3 ≤ 1 / 512 := by nlinarith
  have hε2p : 0 ≤ ε ^ 2 := sq_nonneg ε
  have hε3p : 0 < ε ^ 3 := pow_pos hε0 3
  have ha2 : a ^ 2 ≤ 1 := pow_le_one₀ ha0.le ha1.le
  have ha2p : 0 ≤ a ^ 2 := sq_nonneg a
  have ha3 : a ^ 3 ≤ 1 := pow_le_one₀ ha0.le ha1.le
  have ha3p : 0 ≤ a ^ 3 := by positivity
  have haL : a * L ≤ 1 * 2 := mul_le_mul ha1.le hL1.le hL0.le (by norm_num)
  have haLp : 0 ≤ a * L := by positivity
  have hLa2 : L * a ^ 2 ≤ 2 * 1 := mul_le_mul hL1.le ha2 ha2p (by norm_num)
  have hLa2p : 0 ≤ L * a ^ 2 := by positivity
  have haε : a * ε ≤ 1 * (1 / 8) := mul_le_mul ha1.le hε1 hε0.le (by norm_num)
  have haεp : 0 ≤ a * ε := by positivity
  have hZε : Z * ε ^ 2 ≤ 2 * (1 / 64) := mul_le_mul hZ1.le hε2 hε2p (by norm_num)
  have hZεp : 0 ≤ Z * ε ^ 2 := by positivity
  have haε2 : a * ε ^ 2 ≤ 1 * (1 / 64) := mul_le_mul ha1.le hε2 hε2p (by norm_num)
  have hLε : L * ε ^ 2 ≤ 2 * (1 / 64) := mul_le_mul hL1.le hε2 hε2p (by norm_num)
  have haε2p : 0 ≤ a * ε ^ 2 := by positivity
  have hLεp : 0 ≤ L * ε ^ 2 := by positivity
  have ha2ε : a ^ 2 * ε ^ 2 ≤ 1 * (1 / 64) := mul_le_mul ha2 hε2 hε2p (by norm_num)
  have ha2εp : 0 ≤ a ^ 2 * ε ^ 2 := by positivity
  obtain ⟨c2, hc2d⟩ : ∃ c : ℝ, c = a - L - a ^ 2 / 2 := ⟨_, rfl⟩
  obtain ⟨c3, hc3d⟩ : ∃ c : ℝ, c = a ^ 2 - a ^ 3 / 2 - a * L := ⟨_, rfl⟩
  obtain ⟨c4, hc4d⟩ : ∃ c : ℝ, c = (a ^ 3 - L * a ^ 2) / 2 := ⟨_, rfl⟩
  have hpoly : (1 + Z * ε ^ 2) * (1 - a * ε + (a - L) * ε ^ 2) * (1 + a * ε + a ^ 2 * ε ^ 2 / 2)
      - (1 + (Z + a - a ^ 2 / 2 - L) * ε ^ 2) =
      ε ^ 3 * (c3 + c4 * ε + Z * c2 * ε + Z * c3 * ε ^ 2 + Z * c4 * ε ^ 3) := by
    rw [hc2d, hc3d, hc4d]; ring
  have hc2 : |c2| ≤ 4 := by rw [abs_le]; constructor <;> linarith
  have hc3 : |c3| ≤ 4 := by rw [abs_le]; constructor <;> linarith
  have hc4 : |c4| ≤ 1 := by rw [abs_le]; constructor <;> linarith
  have hZc2 : Z * |c2| ≤ 2 * 4 := mul_le_mul hZ1.le hc2 (abs_nonneg _) (by norm_num)
  have hZc3 : Z * |c3| ≤ 2 * 4 := mul_le_mul hZ1.le hc3 (abs_nonneg _) (by norm_num)
  have hZc4 : Z * |c4| ≤ 2 * 1 := mul_le_mul hZ1.le hc4 (abs_nonneg _) (by norm_num)
  have q1 : |c4 * ε| ≤ 1 / 8 := by
    rw [abs_mul, abs_of_pos hε0]
    have := mul_le_mul hc4 hε1 hε0.le (by norm_num : (0:ℝ) ≤ 1); linarith
  have q2 : |Z * c2 * ε| ≤ 1 := by
    rw [abs_mul, abs_mul, abs_of_pos hε0, abs_of_pos hZ0]
    have := mul_le_mul hZc2 hε1 hε0.le (by norm_num); linarith
  have q3 : |Z * c3 * ε ^ 2| ≤ 1 / 8 := by
    rw [abs_mul, abs_mul, abs_of_nonneg hε2p, abs_of_pos hZ0]
    have := mul_le_mul hZc3 hε2 hε2p (by norm_num); linarith
  have q4 : |Z * c4 * ε ^ 3| ≤ 1 / 8 := by
    rw [abs_mul, abs_mul, abs_of_pos hε3p, abs_of_pos hZ0]
    have := mul_le_mul hZc4 hε3 hε3p.le (by norm_num); linarith
  have hPb : |c3 + c4 * ε + Z * c2 * ε + Z * c3 * ε ^ 2 + Z * c4 * ε ^ 3| ≤ 30 := by
    have t1 := abs_add_le (c3 + c4 * ε + Z * c2 * ε + Z * c3 * ε ^ 2) (Z * c4 * ε ^ 3)
    have t2 := abs_add_le (c3 + c4 * ε + Z * c2 * ε) (Z * c3 * ε ^ 2)
    have t3 := abs_add_le (c3 + c4 * ε) (Z * c2 * ε)
    have t4 := abs_add_le c3 (c4 * ε)
    linarith
  set P := c3 + c4 * ε + Z * c2 * ε + Z * c3 * ε ^ 2 + Z * c4 * ε ^ 3 with hP
  set s0 := 1 + Z * ε ^ 2 with hs0d
  set F0 := 1 - a * ε + (a - L) * ε ^ 2 with hF0d
  set e0 := 1 + a * ε + a ^ 2 * ε ^ 2 / 2 with he0d
  have hs0 : |s0| ≤ 2 := by rw [hs0d, abs_le]; constructor <;> linarith
  have hF0 : |F0| ≤ 2 := by rw [hF0d, abs_le]; constructor <;> linarith
  have he0 : |e0| ≤ 2 := by rw [he0d, abs_le]; constructor <;> linarith
  set rs := s - s0 with hrs
  set rF := F - F0 with hrF
  set re := e - e0 with hre
  have hFb : |F| ≤ 3 := by
    have : F = F0 + rF := by rw [hrF]; ring
    rw [this]; calc |F0 + rF| ≤ |F0| + |rF| := abs_add_le _ _
      _ ≤ 3 := by linarith
  have heb : |e| ≤ 3 := by
    have : e = e0 + re := by rw [hre]; ring
    rw [this]; calc |e0 + re| ≤ |e0| + |re| := abs_add_le _ _
      _ ≤ 3 := by linarith
  have hdecomp : s * F * e - (1 + (Z + a - a ^ 2 / 2 - L) * ε ^ 2) =
      ε ^ 3 * P + rs * F * e + s0 * rF * e + s0 * F0 * re := by
    rw [← hpoly, hrs, hrF, hre]; ring
  rw [hdecomp]
  have t1 : |ε ^ 3 * P| ≤ ε ^ 3 * 30 := by
    rw [abs_mul, abs_of_pos hε3p]; exact mul_le_mul_of_nonneg_left hPb hε3p.le
  have t2 : |rs * F * e| ≤ 3 * ε ^ 3 * 9 := by
    rw [abs_mul, abs_mul]
    have := mul_le_mul hFb heb (abs_nonneg _) (by norm_num)
    calc |rs| * |F| * |e| = |rs| * (|F| * |e|) := by ring
      _ ≤ 3 * ε ^ 3 * (3 * 3) := mul_le_mul hs this (by positivity) (by positivity)
      _ = 3 * ε ^ 3 * 9 := by ring
  have t3 : |s0 * rF * e| ≤ 3 * ε ^ 3 * 6 := by
    rw [abs_mul, abs_mul]
    have := mul_le_mul hs0 heb (abs_nonneg _) (by norm_num)
    calc |s0| * |rF| * |e| = |rF| * (|s0| * |e|) := by ring
      _ ≤ 3 * ε ^ 3 * (2 * 3) := mul_le_mul hF this (by positivity) (by positivity)
      _ = 3 * ε ^ 3 * 6 := by ring
  have t4 : |s0 * F0 * re| ≤ ε ^ 3 * 4 := by
    rw [abs_mul, abs_mul]
    have := mul_le_mul hs0 hF0 (abs_nonneg _) (by norm_num)
    calc |s0| * |F0| * |re| = |re| * (|s0| * |F0|) := by ring
      _ ≤ ε ^ 3 * (2 * 2) := mul_le_mul he this (by positivity) (by positivity)
      _ = ε ^ 3 * 4 := by ring
  have u1 := abs_add_le (ε ^ 3 * P + rs * F * e + s0 * rF * e) (s0 * F0 * re)
  have u2 := abs_add_le (ε ^ 3 * P + rs * F * e) (s0 * rF * e)
  have u3 := abs_add_le (ε ^ 3 * P) (rs * F * e)
  linarith
lemma taylor_core (θ ε : ℝ) (h0 : 0 < θ) (h1 : θ < π) (hε0 : 0 < ε) (hε1 : ε ≤ 1 / 8) :
    |R θ ε * Real.exp (-(a1 θ * ε)) / (2 * zz θ) - (1 + a2 θ * ε ^ 2)| ≤ 3100 * ε ^ 3 := by
  obtain ⟨z0, z1⟩ := zz_bounds θ h0 h1
  have hzz : zz θ = Real.sin (θ / 2) ^ 2 := rfl
  have hU : U θ = (1 - zz θ) / zz θ * Real.log (1 - zz θ) := rfl
  have hUb : -1 < U θ ∧ U θ < 0 := U_bounds z0 z1
  have hLb : 0 < D00.Li2 (zz θ) ∧ D00.Li2 (zz θ) < π ^ 2 / 6 := Li2_pos_lt z0 z1
  have hpi := Real.pi_gt_three
  have hpi' := Real.pi_lt_d2
  obtain ⟨r, hr0, hr1, hLam⟩ := Lam_expand θ ε h0 h1 hε0 hε1
  have hz0 : zz θ ≠ 0 := by rw [hzz]; exact z0.ne'
  have e1 : R θ ε * Real.exp (-(a1 θ * ε)) / (2 * zz θ) =
      (π * ε / Real.sin (π * ε)) *
        (1 - (1 + U θ) * ε + ((1 + U θ) - D00.Li2 (zz θ)) * ε ^ 2 + r) *
        Real.exp ((1 + U θ) * ε) := by
    rw [R, hLam, a1]
    field_simp
    ring_nf
  have e2 : a2 θ = π ^ 2 / 6 + (1 + U θ) - (1 + U θ) ^ 2 / 2 - D00.Li2 (zz θ) := by
    rw [a2]; ring
  rw [e1, e2]
  refine taylor_alg (1 + U θ) (D00.Li2 (zz θ)) ε _ _ _ (by linarith) (by linarith) hLb.1
    (by nlinarith) hε0 hε1 (sinc_bound hε0 hε1) ?_ ?_
  · rw [show 1 - (1 + U θ) * ε + ((1 + U θ) - D00.Li2 (zz θ)) * ε ^ 2 + r -
      (1 - (1 + U θ) * ε + ((1 + U θ) - D00.Li2 (zz θ)) * ε ^ 2) = r by ring,
      abs_of_nonneg hr0]
    nlinarith [pow_pos hε0 3]
  · have hy : |(1 + U θ) * ε| ≤ 1 := by
      rw [abs_of_pos (by nlinarith)]; nlinarith
    have := exp_bound3 hy
    rw [show 1 + (1 + U θ) * ε + ((1 + U θ) * ε) ^ 2 / 2 =
      1 + (1 + U θ) * ε + (1 + U θ) ^ 2 * ε ^ 2 / 2 by ring] at this
    refine this.trans ?_
    rw [abs_of_pos (by nlinarith), mul_pow]
    have : (1 + U θ) ^ 3 ≤ 1 := pow_le_one₀ (by linarith) (by linarith)
    nlinarith [pow_pos hε0 3]

lemma stmt11_aux : stmt_11 := by
  intro θ ε h0 h1 hε0 hε1
  obtain ⟨z0, _⟩ := zz_bounds θ h0 h1
  have hz : 0 < 2 * zz θ := by unfold zz; linarith
  have hz0 : zz θ ≠ 0 := by unfold zz; exact z0.ne'
  have H := taylor_core θ ε h0 h1 hε0 hε1
  have e : R θ ε * Real.exp (-(a1 θ * ε)) - 2 * zz θ * (1 + a2 θ * ε ^ 2) =
      2 * zz θ * (R θ ε * Real.exp (-(a1 θ * ε)) / (2 * zz θ) - (1 + a2 θ * ε ^ 2)) := by
    field_simp
  rw [e, abs_mul, abs_of_pos hz]
  exact mul_le_mul_of_nonneg_left H hz.le

end Taylor

section Asym
open Set

/-- the length formula at one `n`, from the identity `A(ε) = (π/sin πε) Λ_ε` at `ε = eps n`. -/
lemma length_at (θ n : ℝ) (h0 : 0 ≤ θ) (h1 : θ < π) (hn : 1 < n) (hnat : ∀ j : ℕ, eps n ≠ j)
    (hA : Areg θ (eps n) = π / Real.sin (π * eps n) * Lam (eps n) θ) :
    HasSum (serTerm θ n)
      (L θ n - (1 + Real.cos θ + π * eps n / Real.sin (π * eps n) * Lam (eps n) θ * n ^ (-eps n))) := by
  set ε := eps n with hεdef
  have hε : 0 < ε := by rw [hεdef, eps]; exact div_pos one_pos (by linarith)
  have hn0 : 0 < n := by linarith
  obtain ⟨hint, hI⟩ := int_zero_n θ ε n h0 h1 hn hε hnat
  obtain ⟨hhint, hdef⟩ := defect (Real.cos θ) n hn
  rw [← L_eq θ n hn, ← hεdef] at hdef
  have hpow : ∫ s in Ioc 0 n, s ^ (ε - 1) = n ^ ε / ε := by
    rw [← intervalIntegral.integral_of_le hn0.le, integral_rpow (Or.inl (by linarith)),
      show ε - 1 + 1 = ε by ring, Real.zero_rpow hε.ne']
    ring
  have hr : IntegrableOn (fun s : ℝ => s ^ (ε - 1)) (Ioc 0 n) := by
    have := intervalIntegral.intervalIntegrable_rpow' (a := 0) (b := n) (r := ε - 1) (by linarith)
    exact (intervalIntegrable_iff_integrableOn_Ioc_of_le hn0.le).mp this
  have hhk : ∫ s in Ioc 0 n, h (Real.cos θ) s * s ^ (ε - 1) =
      (1 - Real.cos θ) * (n ^ ε / ε) - ∫ s in Ioc 0 n, k (Real.cos θ) s * s ^ (ε - 1) := by
    rw [← hpow, ← integral_const_mul, ← integral_sub (hr.const_mul _) hint]
    refine setIntegral_congr_fun measurableSet_Ioc (fun s _ => ?_)
    simp only [h, k]; ring
  have hnε : n ^ (-ε) * n ^ ε = 1 := by
    rw [← Real.rpow_add hn0, neg_add_cancel, Real.rpow_zero]
  have hL : L θ n = 1 + Real.cos θ + ε * n ^ (-ε) * ∫ s in Ioc 0 n, k (Real.cos θ) s * s ^ (ε - 1) := by
    have : 2 - L θ n = (1 - Real.cos θ) * (n ^ (-ε) * n ^ ε) -
        ε * n ^ (-ε) * ∫ s in Ioc 0 n, k (Real.cos θ) s * s ^ (ε - 1) := by
      rw [hdef, hhk]; field_simp
    rw [hnε] at this; linarith
  have hS1 : Summable (fun i => aCoef θ ε i * n ^ (ε - 1 - i)) := by
    refine Summable.of_norm_bounded ((aCoef_abs_summable θ ε).mul_right (n ^ ε)) (fun i => ?_)
    rw [Real.norm_eq_abs, abs_mul, abs_of_pos (Real.rpow_pos_of_pos (by linarith) _)]
    exact mul_le_mul_of_nonneg_left (Real.rpow_le_rpow_of_exponent_le hn.le
      (by linarith [(Nat.cast_nonneg i : (0:ℝ) ≤ i)])) (abs_nonneg _)
  have H := hS1.hasSum.mul_left (ε * n ^ (-ε))
  rw [hL, hI, hA]
  convert H using 1
  · funext i
    have hne := nonint_ne ε hnat i
    have hne2 : ((i:ℝ) + 1) - ε ≠ 0 := by intro h; apply hne; linarith
    have hpw : n ^ (-ε) * n ^ (ε - 1 - i) = 1 / n ^ (i + 1) := by
      rw [← Real.rpow_add hn0, show -ε + (ε - 1 - i) = -((i + 1 : ℕ) : ℝ) by push_cast; ring,
        Real.rpow_neg hn0.le, Real.rpow_natCast, one_div]
    rw [serTerm, aCoef, ← hεdef]
    have e : ε * n ^ (-ε) * ((-1) ^ i * T13.LamInt (i + 1) θ / (ε - 1 - i) * n ^ (ε - 1 - i)) =
        ε * ((-1) ^ i * T13.LamInt (i + 1) θ / (ε - 1 - i)) * (n ^ (-ε) * n ^ (ε - 1 - i)) := by ring
    rw [e, hpw, pow_succ]
    field_simp
    ring
  · ring

lemma nonint_small {ε : ℝ} (h0 : 0 < ε) (h1 : ε < 1) : ∀ j : ℕ, ε ≠ j := by
  intro j hj
  rcases Nat.eq_zero_or_pos j with h | h
  · subst h; simp at hj; linarith
  · have : (1:ℝ) ≤ j := by exact_mod_cast h
    linarith

lemma eps_facts {n : ℝ} (hn : 9 ≤ n) : 0 < eps n ∧ eps n ≤ 1 / 8 ∧ 1 / n ≤ eps n ∧
    eps n ≤ 2 * (1 / n) ∧ eps n ^ 2 - (1 / n) ^ 2 ≤ 3 * (1 / n) ^ 3 := by
  have hn1 : 0 < n - 1 := by linarith
  have hn0 : 0 < n := by linarith
  unfold eps
  refine ⟨by positivity, ?_, ?_, ?_, ?_⟩
  · rw [div_le_div_iff₀ hn1 (by norm_num)]; linarith
  · rw [div_le_div_iff₀ hn0 hn1]; linarith
  · rw [mul_one_div, div_le_div_iff₀ hn1 hn0]; linarith
  · simp only [div_pow, one_pow]
    rw [div_sub_div _ _ (by positivity) (by positivity), mul_one_div,
      div_le_div_iff₀ (by positivity) (by positivity)]
    have : 0 ≤ n ^ 2 * (n ^ 2 - 5 * n + 3) := mul_nonneg (by positivity) (by nlinarith)
    nlinarith

lemma two_le_log {n : ℝ} (hn : 9 ≤ n) : 2 ≤ Real.log n := by
  rw [Real.le_log_iff_exp_le (by linarith)]
  have h := Real.exp_one_lt_d9
  have e : Real.exp 2 = Real.exp 1 ^ 2 := by rw [← Real.exp_nat_mul]; norm_num
  rw [e]
  nlinarith [Real.exp_pos 1]

/-- `serTerm θ n 0 = -2z(1-z)/(n(n-2))` (`Λ_1 = 2z(1-z)`, `ε/(1-ε) = 1/(n-2)`). -/
lemma serTerm_zero (θ n : ℝ) (hn : 9 ≤ n) :
    serTerm θ n 0 = -(2 * zz θ * (1 - zz θ)) / (n * (n - 2)) := by
  have hzz : zz θ = Real.sin (θ / 2) ^ 2 := rfl
  have hsin : Real.sin θ ^ 2 = 4 * zz θ * (1 - zz θ) := by
    have h2 := Real.sin_two_mul (θ / 2)
    rw [show 2 * (θ / 2) = θ by ring] at h2
    rw [h2, hzz]
    nlinarith [Real.sin_sq_add_cos_sq (θ / 2)]
  have hn2 : n - 2 ≠ 0 := by intro h; linarith
  have hn1 : n - 1 ≠ 0 := by intro h; linarith
  have hn12 : n - 1 - 1 ≠ 0 := by intro h; linarith
  simp only [serTerm, eps, zero_add, pow_one, CharP.cast_eq_zero]
  rw [T13.proof_9 θ, hsin]
  have hn0 : n ≠ 0 := by intro h; linarith
  field_simp
  ring

lemma T3_bound (θ n : ℝ) (h0 : 0 < θ) (h1 : θ < π) (hn : 9 ≤ n) :
    |serTerm θ n 0 + 2 * zz θ * (1 - zz θ) * (1 / n) ^ 2| ≤ zz θ * (1 / n) ^ 3 * 8 := by
  obtain ⟨z0, z1⟩ := zz_bounds θ h0 h1
  have hz0 : 0 < zz θ := z0
  have hz1 : zz θ < 1 := z1
  rw [serTerm_zero θ n hn]
  have hn0 : 0 < n := by linarith
  have hn2 : 0 < n - 2 := by linarith
  have e : -(2 * zz θ * (1 - zz θ)) / (n * (n - 2)) + 2 * zz θ * (1 - zz θ) * (1 / n) ^ 2 =
      -(4 * zz θ * (1 - zz θ)) * ((1 / n) ^ 2 * (1 / (n - 2))) := by
    field_simp; ring
  rw [e, abs_mul, abs_neg, abs_of_pos (by nlinarith), abs_of_pos (by positivity)]
  have h1 : 1 / (n - 2) ≤ 2 * (1 / n) := by
    rw [mul_one_div, div_le_div_iff₀ hn2 hn0]; linarith
  have h2 : (1 / n) ^ 2 * (1 / (n - 2)) ≤ (1 / n) ^ 2 * (2 * (1 / n)) :=
    mul_le_mul_of_nonneg_left h1 (by positivity)
  have h3 : 4 * zz θ * (1 - zz θ) ≤ 4 * zz θ := by nlinarith
  have := mul_le_mul h3 h2 (by positivity) (by positivity)
  calc 4 * zz θ * (1 - zz θ) * ((1 / n) ^ 2 * (1 / (n - 2)))
      ≤ 4 * zz θ * ((1 / n) ^ 2 * (2 * (1 / n))) := this
    _ = zz θ * (1 / n) ^ 3 * 8 := by ring

lemma T4_bound (θ n : ℝ) (h0 : 0 < θ) (h1 : θ < π) (hn : 9 ≤ n) {T : ℝ}
    (hT : HasSum (fun i => serTerm θ n (i + 1)) T) : |T| ≤ zz θ * (1 / n) ^ 3 * 8 := by
  obtain ⟨z0, z1⟩ := zz_bounds θ h0 h1
  have hz0 : 0 < zz θ := z0
  obtain ⟨hε0, hε1, hεw, hε2, hε3⟩ := eps_facts hn
  have hn0 : 0 < n := by linarith
  set z := zz θ with hzdef
  set ε := eps n with hεdef
  set w := 1 / n with hw
  have hw0 : 0 < w := by positivity
  have hw19 : w ≤ 1 / 9 := by rw [hw, div_le_div_iff₀ hn0 (by norm_num)]; linarith
  have hw1 : w < 1 := by linarith
  have hg : HasSum (fun i : ℕ => 2 * z * ε * w ^ 2 * w ^ i) (2 * z * ε * w ^ 2 * (1 - w)⁻¹) :=
    (hasSum_geometric_of_lt_one hw0.le hw1).mul_left _
  have hc : 1 - Real.cos θ = 2 * z := by
    rw [hzdef]; have := T13.sin_half_sq θ; unfold zz; linarith
  have hb : ∀ i : ℕ, |serTerm θ n (i + 1)| ≤ 2 * z * ε * w ^ 2 * w ^ i := by
    intro i
    have hL2 := (le_min_iff.mp (stmt2_aux (i + 1 + 1) (by omega) θ)).1
    rw [hc] at hL2
    have hd : 1 ≤ ((i + 1 : ℕ) : ℝ) + 1 - ε := by
      push_cast; linarith [(Nat.cast_nonneg i : (0:ℝ) ≤ i)]
    have hd0 : (0:ℝ) < ((i + 1 : ℕ) : ℝ) + 1 - ε := by linarith
    simp only [serTerm, ← hεdef]
    rw [abs_div, abs_mul, abs_mul, abs_pow, abs_neg, abs_one, one_pow, one_mul,
      abs_div, abs_of_pos hε0, abs_of_pos hd0,
      abs_of_pos (by positivity : (0:ℝ) < n ^ (i + 1 + 1))]
    have hf : ε / (((i + 1 : ℕ) : ℝ) + 1 - ε) ≤ ε := div_le_self hε0.le hd
    have hp : 1 / n ^ (i + 1 + 1) = w ^ 2 * w ^ i := by
      rw [hw, ← pow_add, one_div_pow]; ring_nf
    rw [div_eq_mul_one_div, hp]
    have := mul_le_mul hf hL2 (abs_nonneg _) hε0.le
    calc ε / (((i + 1 : ℕ) : ℝ) + 1 - ε) * |T13.LamInt (i + 1 + 1) θ| * (w ^ 2 * w ^ i)
        ≤ ε * (2 * z) * (w ^ 2 * w ^ i) := mul_le_mul_of_nonneg_right this (by positivity)
      _ = 2 * z * ε * w ^ 2 * w ^ i := by ring
  have hup := hasSum_le (fun i => (abs_le.mp (hb i)).2) hT hg
  have hlo := hasSum_le (fun i => by have := (abs_le.mp (hb i)).1; linarith : ∀ i : ℕ,
    -serTerm θ n (i + 1) ≤ 2 * z * ε * w ^ 2 * w ^ i) hT.neg hg
  have hgeo : 2 * z * ε * w ^ 2 * (1 - w)⁻¹ ≤ z * w ^ 3 * 8 := by
    have h1 : (1 - w)⁻¹ ≤ 2 := by
      rw [inv_le_comm₀ (by linarith) (by norm_num)]
      linarith
    have h2 : 2 * z * ε * w ^ 2 ≤ 2 * z * (2 * w) * w ^ 2 := by
      have : 0 ≤ 2 * z * w ^ 2 := by positivity
      have := mul_le_mul_of_nonneg_left hε2 this
      linarith
    have hz' : 0 < z := hz0
    have := mul_le_mul h2 h1 (inv_nonneg.mpr (by linarith)) (by positivity)
    calc _ ≤ 2 * z * (2 * w) * w ^ 2 * 2 := this
      _ = z * w ^ 3 * 8 := by ring
  rw [abs_le]; constructor <;> linarith

lemma T1_bound {q ε w lg a z : ℝ} (hq0 : 0 < q) (hq1 : q ≤ 1) (hq2 : 1 - q ≤ (lg + 1) * ε)
    (hε0 : 0 < ε) (hw0 : 0 < w) (hεw : w ≤ ε) (hε3 : ε ^ 2 - w ^ 2 ≤ 3 * w ^ 3)
    (hεc : ε ^ 3 ≤ 8 * w ^ 3) (hlg : 0 < lg) (ha : 0 < a) (ha' : a < 3) (hz : 0 < z) :
    |2 * z * a * (q * ε ^ 2 - w ^ 2)| ≤ z * w ^ 3 * 66 + z * w ^ 3 * lg * 48 := by
  have hB : |q * ε ^ 2 - w ^ 2| ≤ w ^ 3 * 11 + w ^ 3 * lg * 8 := by
    have h1 : 0 ≤ (1 - q) * ε ^ 2 := mul_nonneg (by linarith) (sq_nonneg _)
    have h2 : (1 - q) * ε ^ 2 ≤ (lg + 1) * ε * ε ^ 2 :=
      mul_le_mul_of_nonneg_right hq2 (sq_nonneg _)
    have h3 : (lg + 1) * ε ^ 3 ≤ (lg + 1) * (8 * w ^ 3) :=
      mul_le_mul_of_nonneg_left hεc (by linarith)
    have h4 : w ^ 2 ≤ ε ^ 2 := by gcongr
    have e : (lg + 1) * ε * ε ^ 2 = (lg + 1) * ε ^ 3 := by ring
    rw [abs_le]; constructor <;> nlinarith
  rw [abs_mul, abs_of_pos (by positivity)]
  have h1 : 2 * z * a ≤ 2 * z * 3 := by nlinarith
  have := mul_le_mul h1 hB (abs_nonneg _) (by positivity)
  calc _ ≤ 2 * z * 3 * (w ^ 3 * 11 + w ^ 3 * lg * 8) := this
    _ = _ := by ring

lemma T2_bound {q ρ ε w z : ℝ} (hq0 : 0 < q) (hq1 : q ≤ 1) (hρ : |ρ| ≤ 3100 * ε ^ 3)
    (hεc : ε ^ 3 ≤ 8 * w ^ 3) (hz : 0 < z) : |2 * z * q * ρ| ≤ z * w ^ 3 * 49600 := by
  rw [abs_mul, abs_of_pos (by positivity)]
  have h1 : 2 * z * q ≤ 2 * z := by nlinarith
  have h2 : |ρ| ≤ 3100 * (8 * w ^ 3) := hρ.trans (by linarith)
  have := mul_le_mul h1 h2 (abs_nonneg _) (by positivity)
  calc _ ≤ 2 * z * (3100 * (8 * w ^ 3)) := this
    _ = _ := by ring

lemma final_comb {A B C D z w lg : ℝ} (hz : 0 < z) (hw : 0 < w) (hlg : 2 ≤ lg)
    (hA : |A| ≤ z * w ^ 3 * 66 + z * w ^ 3 * lg * 48) (hB : |B| ≤ z * w ^ 3 * 49600)
    (hC : |C| ≤ z * w ^ 3 * 8) (hD : |D| ≤ z * w ^ 3 * 8) :
    |A + B + C + D| ≤ 30000 * (z * w ^ 3 * lg) := by
  have u1 := abs_add_le (A + B + C) D
  have u2 := abs_add_le (A + B) C
  have u3 := abs_add_le A B
  have hp : 0 ≤ z * w ^ 3 := by positivity
  have hzw : z * w ^ 3 * 2 ≤ z * w ^ 3 * lg := mul_le_mul_of_nonneg_left hlg hp
  linarith

/-- g3:prop:unif, first claim, with the constant `30000` (the paper gets `4100`). -/
lemma E_bound (θ n : ℝ) (h0 : 0 < θ) (h1 : θ < π) (hn : 9 ≤ n) :
    |E θ n| ≤ 30000 * Real.sin (θ / 2) ^ 2 * Real.log n / n ^ 3 := by
  obtain ⟨z0, z1⟩ := zz_bounds θ h0 h1
  have hUb : -1 < U θ ∧ U θ < 0 := U_bounds z0 z1
  have hLb : 0 < D00.Li2 (zz θ) ∧ D00.Li2 (zz θ) < π ^ 2 / 6 := Li2_pos_lt z0 z1
  have hzz : zz θ = Real.sin (θ / 2) ^ 2 := rfl
  rw [← hzz] at z0 z1 ⊢
  obtain ⟨hε0, hε1, hεw, hε2, hε3⟩ := eps_facts hn
  have hρ := taylor_core θ (eps n) h0 h1 hε0 hε1
  have hlen := length_at θ n h0.le h1 (by linarith) (nonint_small hε0 (by linarith))
    (Areg_small θ (eps n) h0.le h1 hε0 (by linarith))
  have hT := (hasSum_nat_add_iff' 1).mpr hlen
  simp only [Finset.sum_range_one] at hT
  have t3 := T3_bound θ n h0 h1 hn
  have t4 := T4_bound θ n h0 h1 hn hT
  have hlog := two_le_log hn
  have hn0 : 0 < n := by linarith
  have ha2 : a2 θ = π ^ 2 / 6 + 1 / 2 - D00.Li2 (zz θ) - U θ ^ 2 / 2 := rfl
  have hpi := Real.pi_lt_d2
  have ha2b : 0 < a2 θ ∧ a2 θ < 3 := by
    rw [ha2]; constructor <;> nlinarith [hUb.1, hUb.2, hLb.1, hLb.2]
  have ha1 : a1 θ = -1 - U θ := rfl
  have hlogn0 : 0 < Real.log n := by linarith
  have hεc : eps n ^ 3 ≤ 8 * (1 / n) ^ 3 := by
    have := pow_le_pow_left₀ hε0.le hε2 3
    calc eps n ^ 3 ≤ (2 * (1 / n)) ^ 3 := this
      _ = 8 * (1 / n) ^ 3 := by ring
  set z := zz θ with hzdef
  set ε := eps n with hεdef
  set w := 1 / n with hw
  have hw0 : 0 < w := by positivity
  set T4 := ∑' i, serTerm θ n (i + 1) with hT4def
  have hT4 : T4 = _ := hT.tsum_eq
  set X := Real.exp (a1 θ * ε) with hX
  set Y := Real.exp (-(ε * Real.log n)) with hY
  have hYn : n ^ (-ε) = Y := by
    rw [Real.rpow_def_of_pos hn0, hY]; congr 1; ring
  have hmu : (mu θ / n) ^ (1 / (n - 1)) = X * Y := by
    rw [show (1:ℝ) / (n - 1) = ε from rfl, mu, Real.rpow_def_of_pos (div_pos (Real.exp_pos _) hn0),
      Real.log_div (Real.exp_pos _).ne' hn0.ne', Real.log_exp, hX, hY, ← Real.exp_add]
    congr 1; ring
  set ρ := R θ ε * Real.exp (-(a1 θ * ε)) / (2 * z) - (1 + a2 θ * ε ^ 2) with hρdef
  have hz0 : z ≠ 0 := z0.ne'
  have hR : π * ε / Real.sin (π * ε) * Lam ε θ = 2 * z * (1 + a2 θ * ε ^ 2 + ρ) * X := by
    have hXX : Real.exp (-(a1 θ * ε)) * X = 1 := by rw [hX, ← Real.exp_add]; simp
    have e1 : 1 + a2 θ * ε ^ 2 + ρ = R θ ε * Real.exp (-(a1 θ * ε)) / (2 * z) := by
      rw [hρdef]; ring
    rw [e1]
    calc π * ε / Real.sin (π * ε) * Lam ε θ = R θ ε * (Real.exp (-(a1 θ * ε)) * X) := by
          rw [hXX, mul_one]; rfl
      _ = _ := by field_simp
  have hE : E θ n = 2 * z * a2 θ * (X * Y * ε ^ 2 - w ^ 2) + 2 * z * (X * Y) * ρ
      + (serTerm θ n 0 + 2 * z * (1 - z) * w ^ 2) + T4 := by
    have hL : L θ n = 1 + Real.cos θ + π * ε / Real.sin (π * ε) * Lam ε θ * n ^ (-ε)
        + serTerm θ n 0 + T4 := by rw [hT4]; ring
    rw [E, hmu, hL, hR, hYn, c2, ← hzdef, hw]
    field_simp
    ring
  have hXY : X * Y = Real.exp ((a1 θ - Real.log n) * ε) := by
    rw [hX, hY, ← Real.exp_add]; congr 1; ring
  have hq0 : 0 < X * Y := by positivity
  have hq1 : X * Y ≤ 1 := by
    rw [hXY, Real.exp_le_one_iff]
    apply mul_nonpos_of_nonpos_of_nonneg _ hε0.le
    rw [ha1]; linarith [hUb.2]
  have hq2 : 1 - X * Y ≤ (Real.log n + 1) * ε := by
    have h := Real.add_one_le_exp ((a1 θ - Real.log n) * ε)
    rw [← hXY] at h
    have h' : -(a1 θ) ≤ 1 := by rw [ha1]; linarith [hUb.1]
    have : -(a1 θ) * ε ≤ 1 * ε := mul_le_mul_of_nonneg_right h' hε0.le
    linarith
  have t1 := T1_bound hq0 hq1 hq2 hε0 hw0 hεw hε3 hεc hlogn0 ha2b.1 ha2b.2 z0
  have t2 := T2_bound (ρ := ρ) hq0 hq1 hρ hεc z0
  have hfin : 30000 * z * Real.log n / n ^ 3 = 30000 * (z * w ^ 3 * Real.log n) := by
    rw [hw]; field_simp
  have t4' : |T4| ≤ z * w ^ 3 * 8 := by rw [hT4]; exact t4
  rw [hE, hfin]
  exact final_comb z0 hw0 hlog t1 t2 t3 t4'

lemma stmt13_aux : stmt_13 := ⟨30000, fun θ n h0 h1 hn => E_bound θ n h0 h1 hn⟩

end Asym

section Nielsen
open Set

/-- `h_{k+1} = Σ_{i<k} 1/(i+1)²` (`h_1 = 0`). -/
noncomputable def hk (k : ℕ) : ℝ := ∑ i ∈ Finset.range k, 1 / ((i:ℝ) + 1) ^ 2

lemma hk_succ (k : ℕ) : hk (k + 1) = hk k + 1 / ((k:ℝ) + 1) ^ 2 := by
  unfold hk; rw [Finset.sum_range_succ]

lemma hk_zero : hk 0 = 0 := by simp [hk]

lemma hk_bounds (k : ℕ) : 0 ≤ hk k ∧ hk k ≤ 2 - 2 / ((k:ℝ) + 1) := by
  induction k with
  | zero => norm_num [hk]
  | succ k ih =>
    rw [hk_succ]
    have key : 1 / ((k:ℝ) + 1) ^ 2 ≤ 2 / ((k:ℝ) + 1) - 2 / ((k:ℝ) + 1 + 1) := by
      rw [div_sub_div _ _ (by positivity) (by positivity),
        div_le_div_iff₀ (by positivity) (by positivity)]
      nlinarith
    push_cast
    have : 0 ≤ 1 / ((k:ℝ) + 1) ^ 2 := by positivity
    constructor
    · linarith [ih.1]
    · linarith [ih.2]

lemma hk_nonneg (k : ℕ) : 0 ≤ hk k := (hk_bounds k).1

lemma hk_le (k : ℕ) : hk k ≤ 2 := by
  have := (hk_bounds k).2
  have : 0 ≤ 2 / ((k:ℝ) + 1) := by positivity
  linarith

/-- `Σ_k h_{k+1} z^k = Li₂(z)/(1-z)`. -/
lemma hk_series {z : ℝ} (hz : |z| < 1) :
    HasSum (fun k => hk k * z ^ k) (D00.Li2 z / (1 - z)) := by
  have hs : Summable (fun k => hk k * z ^ k) := by
    refine Summable.of_norm_bounded
      ((summable_geometric_of_lt_one (abs_nonneg z) hz).mul_left 2) (fun k => ?_)
    rw [Real.norm_eq_abs, abs_mul, abs_of_nonneg (hk_nonneg k), abs_pow]
    exact mul_le_mul_of_nonneg_right (hk_le k) (by positivity)
  set S := ∑' k, hk k * z ^ k with hSdef
  have h1 : HasSum (fun k => hk (k + 1) * z ^ (k + 1)) S := by
    have := (hasSum_nat_add_iff' 1).mpr hs.hasSum
    simpa [hk_zero] using this
  have h2 : HasSum (fun k => hk (k + 1) * z ^ (k + 1)) (z * S + D00.Li2 z) := by
    have := (hs.hasSum.mul_left z).add (D00.proof_6 z hz.le)
    convert this using 1
    funext k; rw [hk_succ, pow_succ]; field_simp; ring
  have hS : S = z * S + D00.Li2 z := h1.unique h2
  have h1z : 1 - z ≠ 0 := by have := (abs_lt.mp hz).2; intro h; linarith
  convert hs.hasSum using 1
  field_simp
  linarith

/-- `Li_{1,2}(z) = Σ_{k≥2} h_k z^k/k`. -/
noncomputable def Mser (z : ℝ) : ℝ := ∑' k : ℕ, hk k * z ^ (k + 1) / ((k:ℝ) + 1)

lemma hasDerivAt_Mser {y : ℝ} (hy : |y| < 1) : HasDerivAt Mser (D00.Li2 y / (1 - y)) y := by
  set r := (|y| + 1) / 2 with hr
  have hr1 : r < 1 := by rw [hr]; linarith
  have hr0 : 0 ≤ r := by rw [hr]; positivity
  have ht : y ∈ Ioo (-r) r := ⟨by linarith [neg_abs_le y], by linarith [le_abs_self y]⟩
  have hd := hasDerivAt_tsum_of_isPreconnected (t := Ioo (-r) r) (y₀ := 0)
    (g := fun (n : ℕ) (z : ℝ) => hk n * z ^ (n + 1) / ((n : ℝ) + 1))
    (g' := fun (n : ℕ) (z : ℝ) => hk n * z ^ n) (u := fun n : ℕ => 2 * r ^ n)
    ((summable_geometric_of_lt_one hr0 hr1).mul_left 2) isOpen_Ioo isPreconnected_Ioo
    (by
      intro n z _
      have := ((hasDerivAt_pow (n + 1) z).const_mul (hk n)).div_const ((n : ℝ) + 1)
      refine this.congr_deriv ?_
      push_cast
      field_simp)
    (by
      intro n z hz
      rw [Real.norm_eq_abs, abs_mul, abs_of_nonneg (hk_nonneg n), abs_pow]
      have hz' : |z| ≤ r := abs_le.mpr ⟨hz.1.le, hz.2.le⟩
      exact mul_le_mul (hk_le n) (pow_le_pow_left₀ (abs_nonneg z) hz' n) (by positivity)
        (by norm_num))
    ⟨by linarith [abs_nonneg y, hr], by linarith [abs_nonneg y, hr]⟩
    (by simp)
    ht
  rw [(hk_series hy).tsum_eq] at hd
  exact hd

/-- `I(z) = ∫_0^z log²(1-t)/t dt`, written with `G(t) = log(1-t)/t`. -/
noncomputable def Iint (z : ℝ) : ℝ := ∫ t in (0:ℝ)..z, D00.G t * Real.log (1 - t)

lemma gfun_contOn : ContinuousOn (fun t => D00.G t * Real.log (1 - t)) (Iio 1) :=
  fun _ ht => (continuousAt_g ht).continuousWithinAt

lemma hasDerivAt_Iint {w : ℝ} (hw : w < 1) :
    HasDerivAt Iint (D00.G w * Real.log (1 - w)) w := by
  have gint : IntervalIntegrable (fun t => D00.G t * Real.log (1 - t)) volume 0 w := by
    apply ContinuousOn.intervalIntegrable
    apply gfun_contOn.mono
    intro t ht
    rw [Set.mem_uIcc] at ht
    show t < 1
    rcases ht with h | h <;> linarith [h.1, h.2]
  exact intervalIntegral.integral_hasDerivAt_right gint
    (gfun_contOn.stronglyMeasurableAtFilter isOpen_Iio w hw) (continuousAt_g hw)

lemma contAt_log1 {w : ℝ} (hw : w < 1) : ContinuousAt (fun w : ℝ => Real.log (1 - w)) w :=
  (continuousAt_const.sub continuousAt_id).log (sub_pos.mpr hw).ne'

/-- g3:lem:sums(iv): `Li_{1,2}(z) = -ℓ Li₂(z) - I(z)` (`I = 2 S_{1,2}`). -/
lemma Mser_eq {z : ℝ} (h0 : 0 ≤ z) (h1 : z < 1) :
    Mser z = -(Real.log (1 - z) * D00.Li2 z) - Iint z := by
  set F : ℝ → ℝ := fun w => Mser w + Real.log (1 - w) * D00.Li2 w + Iint w with hF
  have hF0 : F 0 = 0 := by simp [hF, Mser, Iint]
  suffices F 0 = F z by simp only [hF] at this hF0; linarith
  apply D00.eq_of_deriv_zero h0
  · intro w hw
    have hw1 : w < 1 := by linarith [hw.2]
    have hwa : |w| < 1 := by rw [abs_lt]; constructor <;> linarith [hw.1]
    apply ContinuousAt.continuousWithinAt
    exact (((hasDerivAt_Mser hwa).continuousAt).add
      ((contAt_log1 hw1).mul D00.continuous_Li2.continuousAt)).add
      (hasDerivAt_Iint hw1).continuousAt
  · intro w hw
    have hw1 : w < 1 := by linarith [hw.2]
    have hwa : |w| < 1 := by rw [abs_lt]; constructor <;> linarith [hw.1]
    have hsub : HasDerivAt (fun w : ℝ => 1 - w) (-1) w := by
      simpa using (hasDerivAt_id w).const_sub 1
    have d2 : HasDerivAt (fun w => Real.log (1 - w)) ((-1) / (1 - w)) w :=
      hsub.log (sub_pos.mpr hw1).ne'
    have d3 := D00.hasDerivAt_Li2_G (x := w) hw1.ne
    have hall := ((hasDerivAt_Mser hwa).add (d2.mul d3)).add (hasDerivAt_Iint hw1)
    convert hall using 1
    have h1w : (1:ℝ) - w ≠ 0 := by intro h; linarith
    field_simp
    ring

/-- `g(z) = Σ_k h_{k+1} z^{k+2}/((k+1)(k+2)) = z f_3(z)`. -/
noncomputable def gser (z : ℝ) : ℝ :=
  ∑' k : ℕ, hk k * z ^ (k + 2) / (((k:ℝ) + 1) * ((k:ℝ) + 2))

lemma hasDerivAt_gser {y : ℝ} (hy : |y| < 1) : HasDerivAt gser (Mser y) y := by
  set r := (|y| + 1) / 2 with hr
  have hr1 : r < 1 := by rw [hr]; linarith
  have hr0 : 0 ≤ r := by rw [hr]; positivity
  have ht : y ∈ Ioo (-r) r := ⟨by linarith [neg_abs_le y], by linarith [le_abs_self y]⟩
  have hd := hasDerivAt_tsum_of_isPreconnected (t := Ioo (-r) r) (y₀ := 0)
    (g := fun (n : ℕ) (z : ℝ) => hk n * z ^ (n + 2) / (((n : ℝ) + 1) * ((n : ℝ) + 2)))
    (g' := fun (n : ℕ) (z : ℝ) => hk n * z ^ (n + 1) / ((n : ℝ) + 1)) (u := fun n : ℕ => 2 * r ^ n)
    ((summable_geometric_of_lt_one hr0 hr1).mul_left 2) isOpen_Ioo isPreconnected_Ioo
    (by
      intro n z _
      have := ((hasDerivAt_pow (n + 2) z).const_mul (hk n)).div_const
        (((n : ℝ) + 1) * ((n : ℝ) + 2))
      refine this.congr_deriv ?_
      simp only [show n + 2 - 1 = n + 1 by omega]
      push_cast
      field_simp)
    (by
      intro n z hz
      rw [Real.norm_eq_abs, abs_div, abs_mul, abs_of_nonneg (hk_nonneg n), abs_pow,
        abs_of_pos (by positivity : (0:ℝ) < (n:ℝ) + 1)]
      have hz' : |z| ≤ r := abs_le.mpr ⟨hz.1.le, hz.2.le⟩
      have hr' : r ^ (n + 1) ≤ r ^ n := pow_le_pow_of_le_one hr0 hr1.le (by omega)
      have h1 : |z| ^ (n + 1) ≤ r ^ n := (pow_le_pow_left₀ (abs_nonneg z) hz' _).trans hr'
      have h2 : hk n * |z| ^ (n + 1) ≤ 2 * r ^ n :=
        mul_le_mul (hk_le n) h1 (by positivity) (by norm_num)
      exact (div_le_self (by have := hk_nonneg n; positivity)
        (by have : (0:ℝ) ≤ n := Nat.cast_nonneg n; linarith)).trans h2)
    ⟨by linarith [abs_nonneg y, hr], by linarith [abs_nonneg y, hr]⟩
    (by simp)
    ht
  exact hd

lemma gser_eq {z : ℝ} (h0 : 0 ≤ z) (h1 : z < 1) :
    gser z = (1 - z) * Iint z + (z + (1 - z) * Real.log (1 - z)) * (D00.Li2 z - 1) := by
  set F : ℝ → ℝ := fun w => gser w -
    ((1 - w) * Iint w + (w + (1 - w) * Real.log (1 - w)) * (D00.Li2 w - 1)) with hF
  have hF0 : F 0 = 0 := by simp [hF, gser, Iint]
  suffices F 0 = F z by simp only [hF] at this hF0; linarith
  apply D00.eq_of_deriv_zero h0
  · intro w hw
    have hw1 : w < 1 := by linarith [hw.2]
    have hwa : |w| < 1 := by rw [abs_lt]; constructor <;> linarith [hw.1]
    apply ContinuousAt.continuousWithinAt
    exact ((hasDerivAt_gser hwa).continuousAt).sub
      (((continuousAt_const.sub continuousAt_id).mul (hasDerivAt_Iint hw1).continuousAt).add
        ((continuousAt_id.add ((continuousAt_const.sub continuousAt_id).mul (contAt_log1 hw1))).mul
          (D00.continuous_Li2.continuousAt.sub continuousAt_const)))
  · intro w hw
    have hw1 : w < 1 := by linarith [hw.2]
    have hwa : |w| < 1 := by rw [abs_lt]; constructor <;> linarith [hw.1]
    have hw0 : w ≠ 0 := hw.1.ne'
    have hsub : HasDerivAt (fun w : ℝ => 1 - w) (-1) w := by
      simpa using (hasDerivAt_id w).const_sub 1
    have dlog : HasDerivAt (fun w => Real.log (1 - w)) ((-1) / (1 - w)) w :=
      hsub.log (sub_pos.mpr hw1).ne'
    have dLi := (D00.hasDerivAt_Li2_G (x := w) hw1.ne).sub_const 1
    have P1 := hsub.mul (hasDerivAt_Iint hw1)
    have Q := (hasDerivAt_id' w).add (hsub.mul dlog)
    have hall := (hasDerivAt_gser hwa).sub (P1.add (Q.mul dLi))
    convert hall using 1
    simp only [Pi.add_apply, Pi.mul_apply]
    rw [Mser_eq hw.1.le hw1, D00.G_of_ne hw0]
    have h1w : (1:ℝ) - w ≠ 0 := by intro h; linarith
    field_simp
    ring

/-- `f_3 = Σ_k h_{k+1} z^{k+1}/((k+1)(k+2)) = 2w S_{1,2}(z) + (1 + U)(Li₂ z - 1)`
(g3:lem:sums(iii),(iv) and the proof of g3:lem:a123). -/
lemma f3_hasSum (θ : ℝ) (h0 : 0 < θ) (h1 : θ < π) :
    HasSum (fun k => hk k * (zz θ ^ (k + 1) / (((k:ℝ) + 1) * ((k:ℝ) + 2))))
      (2 * ww θ * S12 θ + (1 + U θ) * (D00.Li2 (zz θ) - 1)) := by
  obtain ⟨z0, z1⟩ := zz_bounds θ h0 h1
  have hS := S12_identity z0 z1
  have hS' : S12 θ = (1 / 2) * Iint (zz θ) := by
    rw [S12, Iint, ell]; unfold zz; rw [hS]; simp_rw [integrand_eq]
  set z := zz θ with hzdef
  have z0' : 0 < z := z0
  have z1' : z < 1 := z1
  have hA := sumA z0' z1'
  have hs : Summable (fun k => hk k * (z ^ (k + 1) / (((k:ℝ) + 1) * ((k:ℝ) + 2)))) := by
    refine Summable.of_nonneg_of_le (fun k => mul_nonneg (hk_nonneg k) (by positivity))
      (fun k => mul_le_mul_of_nonneg_right (hk_le k) (by positivity)) (hA.summable.mul_left 2)
  have hg : gser z = z * ∑' k, hk k * (z ^ (k + 1) / (((k:ℝ) + 1) * ((k:ℝ) + 2))) := by
    rw [gser, ← tsum_mul_left]
    congr 1; funext k; ring
  have hz0 : z ≠ 0 := z0'.ne'
  set T := ∑' k, hk k * (z ^ (k + 1) / (((k:ℝ) + 1) * ((k:ℝ) + 2))) with hTdef
  have hval : T =
      2 * ww θ * S12 θ + (1 + U θ) * (D00.Li2 z - 1) := by
    have e := gser_eq z0'.le z1'
    rw [hg] at e
    rw [ww, U, ww, ell, hS', ← hzdef]
    field_simp
    linear_combination e
  rw [← hval]
  exact hs.hasSum

end Nielsen

section Taylor3
open Set

lemma pp_bounds2 {ε : ℝ} (hε : ε ^ 2 ≤ 1) (k : ℕ) :
    1 - ε ^ 2 * hk k ≤ pp ε k ∧ pp ε k ≤ 1 - ε ^ 2 * hk k + ε ^ 4 * hk k ^ 2 / 2 := by
  induction k with
  | zero => simp [pp, hk]
  | succ k ih =>
    obtain ⟨h1, h2⟩ := ih
    rw [pp_succ, hk_succ]
    have hk0 := hk_nonneg k
    have hy0 : 0 ≤ 1 / ((k:ℝ) + 1) ^ 2 := by positivity
    have hy1 : 1 / ((k:ℝ) + 1) ^ 2 ≤ 1 := by
      apply div_le_one_of_le₀ _ (by positivity)
      nlinarith [(Nat.cast_nonneg k : (0:ℝ) ≤ k)]
    have hx : ε ^ 2 / ((k:ℝ) + 1) ^ 2 = ε ^ 2 * (1 / ((k:ℝ) + 1) ^ 2) := by ring
    rw [hx]
    generalize 1 / ((k:ℝ) + 1) ^ 2 = y at hy0 hy1 ⊢
    have hε2 : 0 ≤ ε ^ 2 := sq_nonneg ε
    have hxy : ε ^ 2 * y ≤ 1 := by nlinarith
    have hf : 0 ≤ 1 - ε ^ 2 * y := by linarith
    have p1 : 0 ≤ ε ^ 4 * hk k * y := by positivity
    have p2 : 0 ≤ ε ^ 6 * hk k ^ 2 * y := by positivity
    have p3 : 0 ≤ ε ^ 4 * y ^ 2 := by positivity
    constructor
    · have := mul_le_mul_of_nonneg_right h1 hf
      nlinarith
    · have := mul_le_mul_of_nonneg_right h2 hf
      nlinarith

/-- `Λ_ε = 2z(1 - aε + (a - Li₂ z)ε² + f₃ε³ + r)`, `|r| ≤ 3ε⁴`, with the closed form of `f₃`. -/
lemma Lam_expand3 (θ ε : ℝ) (h0 : 0 < θ) (h1 : θ < π) (hε0 : 0 < ε) (hε1 : ε ≤ 1 / 8) :
    ∃ r : ℝ, |r| ≤ 3 * ε ^ 4 ∧
      Lam ε θ = 2 * zz θ * (1 - (1 + U θ) * ε + ((1 + U θ) - D00.Li2 (zz θ)) * ε ^ 2
        + (2 * ww θ * S12 θ + (1 + U θ) * (D00.Li2 (zz θ) - 1)) * ε ^ 3 + r) := by
  obtain ⟨z0, z1⟩ := zz_bounds θ h0 h1
  have hzz : zz θ = Real.sin (θ / 2) ^ 2 := rfl
  have hf3 := f3_hasSum θ h0 h1
  have hUb : -1 < U θ ∧ U θ < 0 := U_bounds z0 z1
  set V := 2 * ww θ * S12 θ + (1 + U θ) * (D00.Li2 (zz θ) - 1) with hV
  set z := zz θ with hzdef
  have zpos : 0 < z := z0
  have hA : HasSum (fun k : ℕ => z ^ (k + 1) / (((k:ℝ) + 1) * ((k:ℝ) + 2))) (1 + U θ) := sumA z0 z1
  have hB : HasSum (fun k : ℕ => z ^ (k + 1) / (((k:ℝ) + 1) ^ 2 * ((k:ℝ) + 2)))
      (D00.Li2 z - 1 - U θ) := sumB z0 z1
  set fA := fun k : ℕ => z ^ (k + 1) / (((k:ℝ) + 1) * ((k:ℝ) + 2)) with hfA
  set fB := fun k : ℕ => z ^ (k + 1) / (((k:ℝ) + 1) ^ 2 * ((k:ℝ) + 2)) with hfB
  set cc := fun k : ℕ => z ^ (k + 1) * (ε * ((k:ℝ) + 1 + ε) / (((k:ℝ) + 1) ^ 2 * ((k:ℝ) + 2)))
    with hcc
  set rk := fun k : ℕ => ε ^ 4 * hk k * fB k - cc k * (pp ε k - 1 + ε ^ 2 * hk k) with hrk
  have hz0 : z ≠ 0 := zpos.ne'
  have hsplit : ∀ k : ℕ, T13.lamTerm ε z (k + 1) =
      2 * z * (-ε * fA k - ε ^ 2 * fB k + ε ^ 3 * (hk k * fA k) + rk k) := by
    intro k
    rw [lamTerm_closed]
    simp only [hfA, hfB, hrk, hcc]
    field_simp
    ring
  have hsum : Summable (fun k : ℕ => T13.lamTerm ε z (k + 1)) :=
    (summable_nat_add_iff 1).mpr (lamTerm_summable ε z z0.le z1)
  have hrs : Summable rk := by
    have : Summable (fun k => T13.lamTerm ε z (k + 1) / (2 * z) + ε * fA k + ε ^ 2 * fB k
        - ε ^ 3 * (hk k * fA k)) :=
      (((hsum.div_const _).add (hA.summable.mul_left ε)).add (hB.summable.mul_left (ε ^ 2))).sub
        (hf3.summable.mul_left (ε ^ 3))
    refine this.congr (fun k => ?_)
    rw [hsplit]; field_simp; ring
  have hrk1 : ∀ k, |rk k| ≤ 3 * ε ^ 4 * fA k := fun k => by
    obtain ⟨hp1, hp2⟩ := pp_bounds2 (ε := ε) (by nlinarith) k
    have hk0 := hk_nonneg k
    have hk2 := hk_le k
    have hk' : (1:ℝ) ≤ (k:ℝ) + 1 := by linarith [(Nat.cast_nonneg k : (0:ℝ) ≤ k)]
    have hzk : 0 ≤ z ^ (k + 1) := pow_nonneg zpos.le _
    have hfA0 : 0 ≤ fA k := by simp only [hfA]; positivity
    have hBA : fB k ≤ fA k := by
      simp only [hfA, hfB]
      apply div_le_div_of_nonneg_left hzk (by positivity)
      nlinarith
    have hfB0 : 0 ≤ fB k := by simp only [hfB]; positivity
    have hfrac : ε * ((k:ℝ) + 1 + ε) / (((k:ℝ) + 1) ^ 2 * ((k:ℝ) + 2)) ≤
        9 / 8 * ε / (((k:ℝ) + 1) * ((k:ℝ) + 2)) := by
      rw [div_le_div_iff₀ (by positivity) (by positivity)]
      have : 0 ≤ ε * ((k:ℝ) + 1) * ((k:ℝ) + 2) := by positivity
      nlinarith
    have hc0 : 0 ≤ cc k := by
      simp only [hcc]; exact mul_nonneg hzk (div_nonneg (by positivity) (by positivity))
    have hc1 : cc k ≤ 9 / 8 * ε * fA k := by
      simp only [hcc, hfA]
      calc z ^ (k + 1) * (ε * ((k:ℝ) + 1 + ε) / (((k:ℝ) + 1) ^ 2 * ((k:ℝ) + 2)))
          ≤ z ^ (k + 1) * (9 / 8 * ε / (((k:ℝ) + 1) * ((k:ℝ) + 2))) :=
            mul_le_mul_of_nonneg_left hfrac hzk
        _ = 9 / 8 * ε * (z ^ (k + 1) / (((k:ℝ) + 1) * ((k:ℝ) + 2))) := by ring
    have hd0 : 0 ≤ pp ε k - 1 + ε ^ 2 * hk k := by linarith
    have hd1 : pp ε k - 1 + ε ^ 2 * hk k ≤ 2 * ε ^ 4 := by
      have : ε ^ 4 * hk k ^ 2 / 2 ≤ 2 * ε ^ 4 := by
        have : hk k ^ 2 ≤ 4 := by nlinarith
        have h4 : 0 ≤ ε ^ 4 := by positivity
        nlinarith
      linarith
    have t1 : 0 ≤ ε ^ 4 * hk k * fB k := by positivity
    have t2 : ε ^ 4 * hk k * fB k ≤ 2 * ε ^ 4 * fA k := by
      have h4 : 0 ≤ ε ^ 4 := by positivity
      have := mul_le_mul hk2 hBA hfB0 (by norm_num)
      calc ε ^ 4 * hk k * fB k = ε ^ 4 * (hk k * fB k) := by ring
        _ ≤ ε ^ 4 * (2 * fA k) := mul_le_mul_of_nonneg_left this h4
        _ = 2 * ε ^ 4 * fA k := by ring
    have t3 : 0 ≤ cc k * (pp ε k - 1 + ε ^ 2 * hk k) := mul_nonneg hc0 hd0
    have t4 : cc k * (pp ε k - 1 + ε ^ 2 * hk k) ≤ 9 / 8 * ε * fA k * (2 * ε ^ 4) :=
      mul_le_mul hc1 hd1 hd0 (by positivity)
    have t5 : 9 / 8 * ε * fA k * (2 * ε ^ 4) ≤ ε ^ 4 * fA k := by
      have : 0 ≤ ε ^ 4 * fA k := by positivity
      have e : 9 / 8 * ε * fA k * (2 * ε ^ 4) = (9 / 4 * ε) * (ε ^ 4 * fA k) := by ring
      rw [e]
      have : 9 / 4 * ε ≤ 1 := by linarith
      nlinarith
    simp only [hrk]
    rw [abs_le]; constructor <;> nlinarith
  refine ⟨∑' k, rk k, ?_, ?_⟩
  · have hup := hasSum_le (fun k => (abs_le.mp (hrk1 k)).2) hrs.hasSum (hA.mul_left (3 * ε ^ 4))
    have hlo := hasSum_le (fun k => by have := (abs_le.mp (hrk1 k)).1; linarith :
      ∀ k, -rk k ≤ 3 * ε ^ 4 * fA k) hrs.hasSum.neg (hA.mul_left (3 * ε ^ 4))
    have h4 : 0 < ε ^ 4 := by positivity
    have h5 : 3 * ε ^ 4 * (1 + U θ) ≤ 3 * ε ^ 4 * 1 :=
      mul_le_mul_of_nonneg_left (by linarith [hUb.2]) (by positivity)
    rw [abs_le]; constructor <;> linarith
  · have H : HasSum (fun k : ℕ => T13.lamTerm ε z (k + 1))
        (2 * z * (-ε * (1 + U θ) - ε ^ 2 * (D00.Li2 z - 1 - U θ) + ε ^ 3 * V + ∑' k, rk k)) := by
      have := ((((hA.mul_left (-ε)).sub (hB.mul_left (ε ^ 2))).add
        (hf3.mul_left (ε ^ 3))).add hrs.hasSum).mul_left (2 * z)
      convert this using 1
      funext k; rw [hsplit]
    have H2 := (hasSum_nat_add_iff (f := fun k => T13.lamTerm ε z k) 1).mp H
    unfold Lam
    rw [← hzz, H2.tsum_eq]
    simp only [Finset.sum_range_one, lamTerm_zero]
    ring

/-- `|πε/sin(πε) - (1 + π²ε²/6)| ≤ 5ε⁴` for `0 < ε ≤ 1/8`. -/
lemma sinc_bound4 {ε : ℝ} (hε0 : 0 < ε) (hε1 : ε ≤ 1 / 8) :
    |π * ε / Real.sin (π * ε) - (1 + π ^ 2 / 6 * ε ^ 2)| ≤ 5 * ε ^ 4 := by
  have hpi := Real.pi_gt_three
  have hpi' := Real.pi_lt_d2
  have hx0 : 0 < π * ε := by positivity
  have hx1 : π * ε ≤ 2 / 5 := by nlinarith
  have hb := Real.sin_bound (x := π * ε) (by rw [abs_of_pos hx0]; linarith)
  rw [abs_of_pos hx0] at hb
  have e : π ^ 2 / 6 * ε ^ 2 = (π * ε) ^ 2 / 6 := by ring
  have hx4 : (π * ε) ^ 4 ≤ 99 * ε ^ 4 := by
    rw [mul_pow]
    have h4 := pow_lt_pow_left₀ hpi' Real.pi_pos.le (four_ne_zero)
    have : π ^ 4 ≤ 99 := by linarith [h4, show (3.15:ℝ) ^ 4 ≤ 99 by norm_num]
    exact mul_le_mul_of_nonneg_right this (by positivity)
  rw [e]
  generalize π * ε = x at hx0 hx1 hb hx4 ⊢
  set sx := Real.sin x with hsx
  set δ := sx - (x - x ^ 3 / 6) with hδ
  have hδ' : |δ| ≤ x ^ 5 / 100 := hb
  have hδ2 := abs_le.mp hδ'
  have hx2 : x ^ 2 ≤ 4 / 25 := by nlinarith
  have hx5 : x ^ 5 ≤ x ^ 3 * (4 / 25) := by
    have : x ^ 5 = x ^ 3 * x ^ 2 := by ring
    rw [this]; exact mul_le_mul_of_nonneg_left hx2 (by positivity)
  have hsx0 : 9 / 10 * x ≤ sx := by
    have : x ^ 3 ≤ x * (4 / 25) := by nlinarith
    nlinarith
  have hsxp : 0 < sx := by linarith
  have e2 : x / sx - (1 + x ^ 2 / 6) = (x ^ 5 / 36 - δ * (1 + x ^ 2 / 6)) / sx := by
    field_simp
    rw [hδ]; ring
  rw [e2, abs_div, abs_of_pos hsxp, div_le_iff₀ hsxp]
  have hnum : |x ^ 5 / 36 - δ * (1 + x ^ 2 / 6)| ≤ x ^ 5 / 36 + x ^ 5 / 100 * (1 + x ^ 2 / 6) := by
    have h1 : |δ * (1 + x ^ 2 / 6)| ≤ x ^ 5 / 100 * (1 + x ^ 2 / 6) := by
      rw [abs_mul, abs_of_pos (by positivity : (0:ℝ) < 1 + x ^ 2 / 6)]
      exact mul_le_mul_of_nonneg_right hδ' (by positivity)
    calc _ ≤ |x ^ 5 / 36| + |δ * (1 + x ^ 2 / 6)| := abs_sub _ _
      _ ≤ _ := by rw [abs_of_pos (by positivity)]; linarith
  refine hnum.trans ?_
  have h1 : x ^ 5 / 36 + x ^ 5 / 100 * (1 + x ^ 2 / 6) ≤ x ^ 5 * (39 / 1000) := by
    have : x ^ 5 / 100 * (1 + x ^ 2 / 6) ≤ x ^ 5 / 100 * (1 + (4 / 25) / 6) :=
      mul_le_mul_of_nonneg_left (by linarith) (by positivity)
    nlinarith [pow_pos hx0 5]
  have h2 : x ^ 5 * (39 / 1000) ≤ 5 * ε ^ 4 * (9 / 10 * x) := by
    have : x ^ 5 = x * x ^ 4 := by ring
    rw [this]
    have := mul_le_mul_of_nonneg_left hx4 hx0.le
    nlinarith [pow_pos hε0 4]
  have h3 : 5 * ε ^ 4 * (9 / 10 * x) ≤ 5 * ε ^ 4 * sx :=
    mul_le_mul_of_nonneg_left hsx0 (by positivity)
  linarith

/-- `|e^y - (1 + y + y²/2 + y³/6)| ≤ |y|⁴` for `|y| ≤ 1`. -/
lemma exp_bound4 {y : ℝ} (hy : |y| ≤ 1) :
    |Real.exp y - (1 + y + y ^ 2 / 2 + y ^ 3 / 6)| ≤ |y| ^ 4 := by
  have := Real.exp_bound hy (n := 4) (by norm_num)
  have e : ∑ m ∈ Finset.range 4, y ^ m / (m.factorial : ℝ) = 1 + y + y ^ 2 / 2 + y ^ 3 / 6 := by
    simp [Finset.sum_range_succ, Nat.factorial]
  rw [e] at this
  refine this.trans ?_
  have : ((Nat.succ 4 : ℕ) : ℝ) / ((Nat.factorial 4 : ℕ) * (4 : ℕ)) ≤ 1 := by
    norm_num [Nat.factorial]
  exact mul_le_of_le_one_right (by positivity) (by exact_mod_cast this)

lemma abs_mul_le' {x y X Y : ℝ} (hx : |x| ≤ X) (hy : |y| ≤ Y) : |x * y| ≤ X * Y := by
  rw [abs_mul]; exact mul_le_mul hx hy (abs_nonneg _) ((abs_nonneg _).trans hx)

/-- the algebra of the proof of g3:prop:n3 (order 4). -/
lemma taylor_alg3 (a L Z f ε s F e : ℝ) (ha0 : 0 < a) (ha1 : a < 1) (hL0 : 0 < L) (hL1 : L < 2)
    (hZ0 : 0 < Z) (hZ1 : Z < 2) (hf : |f| ≤ 2) (hε0 : 0 < ε) (hε1 : ε ≤ 1 / 8)
    (hs : |s - (1 + Z * ε ^ 2)| ≤ 5 * ε ^ 4)
    (hF : |F - (1 - a * ε + (a - L) * ε ^ 2 + f * ε ^ 3)| ≤ 3 * ε ^ 4)
    (he : |e - (1 + a * ε + a ^ 2 * ε ^ 2 / 2 + a ^ 3 * ε ^ 3 / 6)| ≤ ε ^ 4) :
    |s * F * e - (1 + (Z + a - a ^ 2 / 2 - L) * ε ^ 2 + (f + a ^ 2 - a * L - a ^ 3 / 3) * ε ^ 3)|
      ≤ 31000 * ε ^ 4 := by
  have ha : |a| ≤ 1 := by rw [abs_le]; constructor <;> linarith
  have hZ : |Z| ≤ 2 := by rw [abs_le]; constructor <;> linarith
  have hε : |ε| ≤ 1 / 8 := by rw [abs_of_pos hε0]; exact hε1
  have ha2 : |a ^ 2| ≤ 1 := by rw [abs_pow]; exact pow_le_one₀ (abs_nonneg _) ha
  have ha3 : |a ^ 3| ≤ 1 := by rw [abs_pow]; exact pow_le_one₀ (abs_nonneg _) ha
  have ha4 : |a ^ 4| ≤ 1 := by rw [abs_pow]; exact pow_le_one₀ (abs_nonneg _) ha
  have hε2 : |ε ^ 2| ≤ 1 / 64 := by
    rw [abs_pow]; calc |ε| ^ 2 ≤ (1 / 8) ^ 2 := pow_le_pow_left₀ (abs_nonneg _) hε 2
      _ = 1 / 64 := by norm_num
  have hε3 : |ε ^ 3| ≤ 1 / 512 := by
    rw [abs_pow]; calc |ε| ^ 3 ≤ (1 / 8) ^ 3 := pow_le_pow_left₀ (abs_nonneg _) hε 3
      _ = 1 / 512 := by norm_num
  have hε4 : |ε ^ 4| ≤ 1 / 4096 := by
    rw [abs_pow]; calc |ε| ^ 4 ≤ (1 / 8) ^ 4 := pow_le_pow_left₀ (abs_nonneg _) hε 4
      _ = 1 / 4096 := by norm_num
  have hε4p : 0 < ε ^ 4 := by positivity
  have hε4' : ε ^ 4 ≤ 1 / 4096 := (le_abs_self _).trans hε4
  obtain ⟨b, hbd⟩ : ∃ b : ℝ, b = a - L := ⟨_, rfl⟩
  have hb : |b| ≤ 2 := by rw [hbd, abs_le]; constructor <;> linarith
  rw [← hbd] at hF
  have e2 : Z + a - a ^ 2 / 2 - L = Z + b - a ^ 2 / 2 := by rw [hbd]; ring
  have e3 : f + a ^ 2 - a * L - a ^ 3 / 3 = f + a * b - a ^ 3 / 3 := by rw [hbd]; ring
  rw [e2, e3]
  obtain ⟨c2, hc2⟩ : ∃ c : ℝ, c = b - a ^ 2 / 2 := ⟨_, rfl⟩
  obtain ⟨c3, hc3⟩ : ∃ c : ℝ, c = f + a * b - a ^ 3 / 3 := ⟨_, rfl⟩
  obtain ⟨c4, hc4⟩ : ∃ c : ℝ, c = f * a + b * a ^ 2 / 2 - a ^ 4 / 6 := ⟨_, rfl⟩
  obtain ⟨c5, hc5⟩ : ∃ c : ℝ, c = b * a ^ 3 / 6 + f * a ^ 2 / 2 := ⟨_, rfl⟩
  obtain ⟨c6, hc6⟩ : ∃ c : ℝ, c = f * a ^ 3 / 6 := ⟨_, rfl⟩
  have hpoly : (1 + Z * ε ^ 2) * (1 - a * ε + b * ε ^ 2 + f * ε ^ 3) *
      (1 + a * ε + a ^ 2 * ε ^ 2 / 2 + a ^ 3 * ε ^ 3 / 6)
      - (1 + (Z + b - a ^ 2 / 2) * ε ^ 2 + (f + a * b - a ^ 3 / 3) * ε ^ 3)
      = ε ^ 4 * ((c4 + Z * c2) + (c5 + Z * c3) * ε + (c6 + Z * c4) * ε ^ 2 + Z * c5 * ε ^ 3
        + Z * c6 * ε ^ 4) := by
    rw [hc2, hc3, hc4, hc5, hc6]; ring
  have pab := abs_le.mp (abs_mul_le' ha hb)
  have pfa := abs_le.mp (abs_mul_le' hf ha)
  have pba2 := abs_le.mp (abs_mul_le' hb ha2)
  have pba3 := abs_le.mp (abs_mul_le' hb ha3)
  have pfa2 := abs_le.mp (abs_mul_le' hf ha2)
  have pfa3 := abs_le.mp (abs_mul_le' hf ha3)
  have qa2 := abs_le.mp ha2
  have qa3 := abs_le.mp ha3
  have qa4 := abs_le.mp ha4
  have qb := abs_le.mp hb
  have qf := abs_le.mp hf
  have hC2 : |c2| ≤ 3 := by rw [hc2, abs_le]; constructor <;> linarith
  have hC3 : |c3| ≤ 5 := by rw [hc3, abs_le]; constructor <;> linarith
  have hC4 : |c4| ≤ 4 := by rw [hc4, abs_le]; constructor <;> linarith
  have hC5 : |c5| ≤ 2 := by rw [hc5, abs_le]; constructor <;> linarith
  have hC6 : |c6| ≤ 1 := by rw [hc6, abs_le]; constructor <;> linarith
  have z2 := abs_le.mp (abs_mul_le' hZ hC2)
  have z3' := abs_le.mp (abs_mul_le' hZ hC3)
  have z4 := abs_le.mp (abs_mul_le' hZ hC4)
  have qc4 := abs_le.mp hC4
  have qc5 := abs_le.mp hC5
  have qc6 := abs_le.mp hC6
  have hD4 : |c4 + Z * c2| ≤ 10 := by rw [abs_le]; constructor <;> linarith
  have hD5 : |c5 + Z * c3| ≤ 12 := by rw [abs_le]; constructor <;> linarith
  have hD6 : |c6 + Z * c4| ≤ 9 := by rw [abs_le]; constructor <;> linarith
  have hD7 : |Z * c5| ≤ 2 * 2 := abs_mul_le' hZ hC5
  have hD8 : |Z * c6| ≤ 2 * 1 := abs_mul_le' hZ hC6
  have r5 := abs_mul_le' hD5 hε
  have r6 := abs_mul_le' hD6 hε2
  have r7 := abs_mul_le' hD7 hε3
  have r8 := abs_mul_le' hD8 hε4
  have hPb : |(c4 + Z * c2) + (c5 + Z * c3) * ε + (c6 + Z * c4) * ε ^ 2 + Z * c5 * ε ^ 3
      + Z * c6 * ε ^ 4| ≤ 13 := by
    have u1 := abs_add_le ((c4 + Z * c2) + (c5 + Z * c3) * ε + (c6 + Z * c4) * ε ^ 2
      + Z * c5 * ε ^ 3) (Z * c6 * ε ^ 4)
    have u2 := abs_add_le ((c4 + Z * c2) + (c5 + Z * c3) * ε + (c6 + Z * c4) * ε ^ 2)
      (Z * c5 * ε ^ 3)
    have u3 := abs_add_le ((c4 + Z * c2) + (c5 + Z * c3) * ε) ((c6 + Z * c4) * ε ^ 2)
    have u4 := abs_add_le (c4 + Z * c2) ((c5 + Z * c3) * ε)
    linarith
  set P := (c4 + Z * c2) + (c5 + Z * c3) * ε + (c6 + Z * c4) * ε ^ 2 + Z * c5 * ε ^ 3
      + Z * c6 * ε ^ 4 with hP
  set s0 := 1 + Z * ε ^ 2 with hs0d
  set F0 := 1 - a * ε + b * ε ^ 2 + f * ε ^ 3 with hF0d
  set e0 := 1 + a * ε + a ^ 2 * ε ^ 2 / 2 + a ^ 3 * ε ^ 3 / 6 with he0d
  have m1 := abs_le.mp (abs_mul_le' hZ hε2)
  have m2 := abs_le.mp (abs_mul_le' ha hε)
  have m3 := abs_le.mp (abs_mul_le' hb hε2)
  have m4 := abs_le.mp (abs_mul_le' hf hε3)
  have m5 := abs_le.mp (abs_mul_le' ha2 hε2)
  have m6 := abs_le.mp (abs_mul_le' ha3 hε3)
  have hs0 : |s0| ≤ 2 := by rw [hs0d, abs_le]; constructor <;> linarith
  have hF0 : |F0| ≤ 2 := by rw [hF0d, abs_le]; constructor <;> linarith
  have he0 : |e0| ≤ 2 := by rw [he0d, abs_le]; constructor <;> linarith
  set rs := s - s0 with hrs
  set rF := F - F0 with hrF
  set re := e - e0 with hre
  have hFb : |F| ≤ 3 := by
    have : F = F0 + rF := by rw [hrF]; ring
    rw [this]; calc |F0 + rF| ≤ |F0| + |rF| := abs_add_le _ _
      _ ≤ 3 := by linarith
  have heb : |e| ≤ 3 := by
    have : e = e0 + re := by rw [hre]; ring
    rw [this]; calc |e0 + re| ≤ |e0| + |re| := abs_add_le _ _
      _ ≤ 3 := by linarith
  have hdecomp : s * F * e - (1 + (Z + b - a ^ 2 / 2) * ε ^ 2 + (f + a * b - a ^ 3 / 3) * ε ^ 3) =
      ε ^ 4 * P + rs * F * e + s0 * rF * e + s0 * F0 * re := by
    rw [← hpoly, hrs, hrF, hre]; ring
  rw [hdecomp]
  have t1 : |ε ^ 4 * P| ≤ ε ^ 4 * 13 := by
    rw [abs_mul, abs_of_pos hε4p]; exact mul_le_mul_of_nonneg_left hPb hε4p.le
  have t2 : |rs * F * e| ≤ 5 * ε ^ 4 * 9 := by
    rw [abs_mul, abs_mul]
    have := mul_le_mul hFb heb (abs_nonneg _) (by norm_num)
    calc |rs| * |F| * |e| = |rs| * (|F| * |e|) := by ring
      _ ≤ 5 * ε ^ 4 * (3 * 3) := mul_le_mul hs this (by positivity) (by positivity)
      _ = 5 * ε ^ 4 * 9 := by ring
  have t3 : |s0 * rF * e| ≤ 3 * ε ^ 4 * 6 := by
    rw [abs_mul, abs_mul]
    have := mul_le_mul hs0 heb (abs_nonneg _) (by norm_num)
    calc |s0| * |rF| * |e| = |rF| * (|s0| * |e|) := by ring
      _ ≤ 3 * ε ^ 4 * (2 * 3) := mul_le_mul hF this (by positivity) (by positivity)
      _ = 3 * ε ^ 4 * 6 := by ring
  have t4 : |s0 * F0 * re| ≤ ε ^ 4 * 4 := by
    rw [abs_mul, abs_mul]
    have := mul_le_mul hs0 hF0 (abs_nonneg _) (by norm_num)
    calc |s0| * |F0| * |re| = |re| * (|s0| * |F0|) := by ring
      _ ≤ ε ^ 4 * (2 * 2) := mul_le_mul he this (by positivity) (by positivity)
      _ = ε ^ 4 * 4 := by ring
  have u1 := abs_add_le (ε ^ 4 * P + rs * F * e + s0 * rF * e) (s0 * F0 * re)
  have u2 := abs_add_le (ε ^ 4 * P + rs * F * e) (s0 * rF * e)
  have u3 := abs_add_le (ε ^ 4 * P) (rs * F * e)
  linarith

lemma taylor_core3 (θ ε : ℝ) (h0 : 0 < θ) (h1 : θ < π) (hε0 : 0 < ε) (hε1 : ε ≤ 1 / 8) :
    |R θ ε * Real.exp (-(a1 θ * ε)) / (2 * zz θ) - (1 + a2 θ * ε ^ 2 + a3 θ * ε ^ 3)|
      ≤ 31000 * ε ^ 4 := by
  obtain ⟨z0, z1⟩ := zz_bounds θ h0 h1
  have hUb : -1 < U θ ∧ U θ < 0 := U_bounds z0 z1
  have hLb : 0 < D00.Li2 (zz θ) ∧ D00.Li2 (zz θ) < π ^ 2 / 6 := Li2_pos_lt z0 z1
  have hpi := Real.pi_gt_three
  have hpi' := Real.pi_lt_d2
  obtain ⟨r, hr, hLam⟩ := Lam_expand3 θ ε h0 h1 hε0 hε1
  -- bound on f₃
  have hf3 := f3_hasSum θ h0 h1
  have hA : HasSum (fun k : ℕ => zz θ ^ (k + 1) / (((k:ℝ) + 1) * ((k:ℝ) + 2))) (1 + U θ) :=
    sumA z0 z1
  set V := 2 * ww θ * S12 θ + (1 + U θ) * (D00.Li2 (zz θ) - 1) with hV
  have hz0' : 0 < zz θ := z0
  have hV0 : 0 ≤ V := hasSum_le (fun k => mul_nonneg (hk_nonneg k) (by positivity)) hasSum_zero hf3
  have hV1 : V ≤ 2 * (1 + U θ) := hasSum_le (fun k => mul_le_mul_of_nonneg_right (hk_le k)
    (by positivity)) hf3 (hA.mul_left 2)
  have hVb : |V| ≤ 2 := by rw [abs_le]; constructor <;> linarith
  have hz0 : zz θ ≠ 0 := hz0'.ne'
  have e1 : R θ ε * Real.exp (-(a1 θ * ε)) / (2 * zz θ) =
      (π * ε / Real.sin (π * ε)) *
        (1 - (1 + U θ) * ε + ((1 + U θ) - D00.Li2 (zz θ)) * ε ^ 2 + V * ε ^ 3 + r) *
        Real.exp ((1 + U θ) * ε) := by
    rw [R, hLam, a1]
    field_simp
    ring_nf
  have e2 : a2 θ = π ^ 2 / 6 + (1 + U θ) - (1 + U θ) ^ 2 / 2 - D00.Li2 (zz θ) := by
    rw [a2]; ring
  have e3 : a3 θ = V + (1 + U θ) ^ 2 - (1 + U θ) * D00.Li2 (zz θ) - (1 + U θ) ^ 3 / 3 := by
    rw [a3, hV]; ring
  rw [e1, e2, e3]
  refine taylor_alg3 (1 + U θ) (D00.Li2 (zz θ)) (π ^ 2 / 6) V ε _ _ _ (by linarith) (by linarith)
    hLb.1 (by nlinarith) (by positivity) (by nlinarith) hVb hε0 hε1 (sinc_bound4 hε0 hε1) ?_ ?_
  · rw [show 1 - (1 + U θ) * ε + ((1 + U θ) - D00.Li2 (zz θ)) * ε ^ 2 + V * ε ^ 3 + r -
      (1 - (1 + U θ) * ε + ((1 + U θ) - D00.Li2 (zz θ)) * ε ^ 2 + V * ε ^ 3) = r by ring]
    exact hr
  · have hy : |(1 + U θ) * ε| ≤ 1 := by
      rw [abs_of_pos (by nlinarith)]; nlinarith
    have := exp_bound4 hy
    rw [show 1 + (1 + U θ) * ε + ((1 + U θ) * ε) ^ 2 / 2 + ((1 + U θ) * ε) ^ 3 / 6 =
      1 + (1 + U θ) * ε + (1 + U θ) ^ 2 * ε ^ 2 / 2 + (1 + U θ) ^ 3 * ε ^ 3 / 6 by ring] at this
    refine this.trans ?_
    rw [abs_of_pos (by nlinarith), mul_pow]
    have : (1 + U θ) ^ 4 ≤ 1 := pow_le_one₀ (by linarith) (by linarith)
    exact mul_le_of_le_one_left (by positivity) this

lemma stmt12_aux : stmt_12 := by
  intro θ ε h0 h1 hε0 hε1
  obtain ⟨z0, _⟩ := zz_bounds θ h0 h1
  have hz : 0 < 2 * zz θ := by unfold zz; linarith
  have hz0 : zz θ ≠ 0 := by unfold zz; exact z0.ne'
  have H := taylor_core3 θ ε h0 h1 hε0 hε1
  have e : R θ ε * Real.exp (-(a1 θ * ε)) - 2 * zz θ * (1 + a2 θ * ε ^ 2 + a3 θ * ε ^ 3) =
      2 * zz θ * (R θ ε * Real.exp (-(a1 θ * ε)) / (2 * zz θ) -
        (1 + a2 θ * ε ^ 2 + a3 θ * ε ^ 3)) := by
    field_simp
  rw [e, abs_mul, abs_of_pos hz]
  exact mul_le_mul_of_nonneg_left H hz.le

end Taylor3

section ContA
open Set

def rect (a b c d : ℝ) : Set ℂ := {z | a < z.re ∧ z.re < b ∧ c < z.im ∧ z.im < d}

lemma rect_open (a b c d : ℝ) : IsOpen (rect a b c d) := by
  have hr := Complex.continuous_re; have hi := Complex.continuous_im
  exact (isOpen_lt continuous_const hr).and ((isOpen_lt hr continuous_const).and
    ((isOpen_lt continuous_const hi).and (isOpen_lt hi continuous_const)))

lemma rect_convex (a b c d : ℝ) : Convex ℝ (rect a b c d) := by
  have e : rect a b c d = Complex.reLm ⁻¹' Ioo a b ∩ Complex.imLm ⁻¹' Ioo c d := by
    ext z; simp [rect, and_assoc]
  rw [e]
  exact ((convex_Ioo a b).linear_preimage _).inter ((convex_Ioo c d).linear_preimage _)

lemma norm_sub_int_re (z : ℂ) (n : ℤ) : |z.re - n| ≤ ‖z - n‖ := by
  have := Complex.abs_re_le_norm (z - n); simpa using this

lemma norm_sub_int_im (z : ℂ) (n : ℤ) : |z.im| ≤ ‖z - n‖ := by
  have := Complex.abs_im_le_norm (z - n); simpa using this

lemma domain (ε0 : ℝ) (h1 : 1 < ε0) (hnat : ∀ j : ℕ, ε0 ≠ j) :
    ∃ U : Set ℂ, ∃ δ : ℝ, 0 < δ ∧ IsOpen U ∧ IsPreconnected U ∧ (1/2 : ℂ) ∈ U ∧ (ε0 : ℂ) ∈ U ∧
      ∀ z ∈ U, 0 < z.re ∧ ‖z‖ ≤ ε0 + 3 ∧ ∀ n : ℤ, δ ≤ ‖z - n‖ := by
  set f := Int.fract ε0 with hf
  set m := ⌊ε0⌋ with hm
  have hf0 : 0 < f := by
    rcases (Int.fract_nonneg ε0).lt_or_eq with h | h
    · exact h
    · exfalso
      have hm0 : 0 ≤ m := Int.floor_nonneg.mpr (by linarith)
      apply hnat m.toNat
      have : ε0 = m := by
        have := Int.floor_add_fract ε0; rw [← h] at this; linarith
      rw [this]; exact_mod_cast (Int.toNat_of_nonneg hm0).symm
  have hf1 : f < 1 := Int.fract_lt_one ε0
  have hεm : ε0 = m + f := (Int.floor_add_fract ε0).symm
  have hm0 : (0:ℝ) ≤ m := by exact_mod_cast Int.floor_nonneg.mpr (by linarith)
  set δ := min f (1 - f) / 2 with hδ
  have hδ0 : 0 < δ := by rw [hδ]; apply half_pos; exact lt_min hf0 (by linarith)
  have hδf : 2 * δ ≤ f := by rw [hδ]; linarith [min_le_left f (1 - f)]
  have hδf' : 2 * δ ≤ 1 - f := by rw [hδ]; linarith [min_le_right f (1 - f)]
  refine ⟨rect (1/4) (3/4) (-1) 1 ∪ rect (1/4) (ε0 + 1) (1/2) 1 ∪ rect (ε0 - δ) (ε0 + δ) (-1) 1,
    δ, hδ0, ((rect_open _ _ _ _).union (rect_open _ _ _ _)).union (rect_open _ _ _ _), ?_, ?_, ?_, ?_⟩
  · have h1 : IsPreconnected (rect (1/4) (3/4) (-1) 1 ∪ rect (1/4) (ε0 + 1) (1/2) 1) := by
      refine IsPreconnected.union ((1/2 : ℂ) + (3/4 : ℂ) * Complex.I) ?_ ?_
        (rect_convex _ _ _ _).isPreconnected (rect_convex _ _ _ _).isPreconnected
      · simp [rect]; norm_num
      · simp [rect]; norm_num; linarith
    refine IsPreconnected.union ((ε0 : ℂ) + (3/4 : ℂ) * Complex.I) ?_ ?_ h1
      (rect_convex _ _ _ _).isPreconnected
    · right; simp [rect]; norm_num; linarith
    · simp [rect]; norm_num; exact hδ0
  · left; left; simp [rect]; norm_num
  · right; simp [rect]; exact hδ0
  · intro z hz
    have hn : ‖z‖ ≤ |z.re| + |z.im| := Complex.norm_le_abs_re_add_abs_im z
    rcases hz with (hz | hz) | hz <;> obtain ⟨h1, h2, h3, h4⟩ := hz
    · refine ⟨by linarith, ?_, fun n => ?_⟩
      · rw [abs_of_pos (by linarith : 0 < z.re)] at hn
        have : |z.im| ≤ 1 := abs_le.mpr ⟨h3.le, h4.le⟩
        linarith
      · refine le_trans ?_ (norm_sub_int_re z n)
        rcases le_or_gt n 0 with hn0 | hn0
        · have : (n:ℝ) ≤ 0 := by exact_mod_cast hn0
          rw [abs_of_pos (by linarith)]; linarith
        · have : (1:ℝ) ≤ n := by exact_mod_cast hn0
          rw [abs_of_neg (by linarith)]; linarith
    · refine ⟨by linarith, ?_, fun n => ?_⟩
      · have : |z.im| ≤ 1 := abs_le.mpr ⟨by linarith, h4.le⟩
        rw [abs_of_pos (by linarith : 0 < z.re)] at hn
        linarith
      · refine le_trans ?_ (norm_sub_int_im z n)
        rw [abs_of_pos (by linarith)]; linarith
    · refine ⟨by linarith, ?_, fun n => ?_⟩
      · have : |z.im| ≤ 1 := abs_le.mpr ⟨h3.le, h4.le⟩
        rw [abs_of_pos (by linarith : 0 < z.re)] at hn
        linarith
      · refine le_trans ?_ (norm_sub_int_re z n)
        rcases le_or_gt n m with hn0 | hn0
        · have : (n:ℝ) ≤ m := by exact_mod_cast hn0
          rw [abs_of_pos (by linarith)]; linarith
        · have : (m:ℝ) + 1 ≤ n := by exact_mod_cast hn0
          rw [abs_of_neg (by linarith)]; linarith

/-- `‖ε - 1 - i‖ ≥ κ (i+1)` on a set bounded by `B` and `δ`-away from the integers. -/
lemma denom_lower {δ B : ℝ} (hδ : 0 < δ) (hB : 0 < B) {ε : ℂ} (hεB : ‖ε‖ ≤ B)
    (hd : ∀ n : ℤ, δ ≤ ‖ε - n‖) (i : ℕ) :
    min (1/2) (δ / (2 * B)) * ((i:ℝ) + 1) ≤ ‖ε - 1 - i‖ := by
  have hi : (0:ℝ) ≤ i := Nat.cast_nonneg i
  rcases le_or_gt (2 * B) ((i:ℝ) + 1) with h | h
  · have : ((i:ℝ) + 1) - ‖ε‖ ≤ ‖ε - 1 - i‖ := by
      have := norm_sub_norm_le ((1:ℂ) + i) ε
      rw [show ε - 1 - (i:ℂ) = -((1:ℂ) + i - ε) by ring, norm_neg]
      have e : ‖(1:ℂ) + (i:ℂ)‖ = (i:ℝ) + 1 := by
        rw [show (1:ℂ) + (i:ℂ) = ((i + 1 : ℕ) : ℂ) by push_cast; ring, Complex.norm_natCast]
        push_cast; ring
      linarith
    have : min (1/2) (δ / (2 * B)) * ((i:ℝ) + 1) ≤ (1/2) * ((i:ℝ) + 1) :=
      mul_le_mul_of_nonneg_right (min_le_left _ _) (by linarith)
    linarith
  · have h1 := hd ((i:ℤ) + 1)
    rw [show ε - (((i:ℤ) + 1 : ℤ) : ℂ) = ε - 1 - i by push_cast; ring] at h1
    have : min (1/2) (δ / (2 * B)) * ((i:ℝ) + 1) ≤ (δ / (2 * B)) * ((i:ℝ) + 1) :=
      mul_le_mul_of_nonneg_right (min_le_right _ _) (by linarith)
    have : (δ / (2 * B)) * ((i:ℝ) + 1) ≤ δ := by
      rw [div_mul_eq_mul_div, div_le_iff₀ (by positivity)]
      nlinarith
    linarith

/-- the series `Σ_i c_i/(ε - 1 - i)` is holomorphic on `U`. -/
lemma aSum_diff {δ B : ℝ} (hδ : 0 < δ) (hB : 0 < B) {U : Set ℂ} (hU : IsOpen U)
    (hUp : ∀ z ∈ U, ‖z‖ ≤ B ∧ ∀ n : ℤ, δ ≤ ‖z - n‖) (c : ℕ → ℝ)
    (hc : ∀ i, |c i| ≤ 2 / (2 * (i:ℝ) + 3)) :
    DifferentiableOn ℂ (fun w : ℂ => ∑' i : ℕ, (c i : ℂ) / (w - 1 - i)) U := by
  set κ := min (1/2) (δ / (2 * B)) with hκ
  have hκ0 : 0 < κ := lt_min (by norm_num) (by positivity)
  have hne : ∀ w ∈ U, ∀ i : ℕ, w - 1 - i ≠ 0 := by
    intro w hw i h
    have := denom_lower hδ hB (hUp w hw).1 (hUp w hw).2 i
    rw [h, norm_zero] at this
    have : 0 < κ * ((i:ℝ) + 1) := by positivity
    linarith
  refine Complex.differentiableOn_tsum_of_summable_norm
    (F := fun (i : ℕ) (w : ℂ) => (c i : ℂ) / (w - 1 - i)) (u := fun i : ℕ => (2 / κ) * (1 / ((i:ℝ) + 1) ^ 2))
    ?_ ?_ hU ?_
  · apply Summable.mul_left
    exact_mod_cast (summable_nat_add_iff 1).mpr (Real.summable_one_div_nat_pow.mpr one_lt_two)
  · intro i
    have hd : DifferentiableOn ℂ (fun w : ℂ => w - 1 - (i:ℂ)) U := by fun_prop
    exact (differentiableOn_const (c i : ℂ)).div hd (fun w hw => hne w hw i)
  · intro i w hw
    have hl := denom_lower hδ hB (hUp w hw).1 (hUp w hw).2 i
    have hi : (0:ℝ) ≤ i := Nat.cast_nonneg i
    show ‖(c i : ℂ) / (w - 1 - i)‖ ≤ (2 / κ) * (1 / ((i:ℝ) + 1) ^ 2)
    rw [norm_div, Complex.norm_real, Real.norm_eq_abs]
    calc |c i| / ‖w - 1 - i‖ ≤ (2 / (2 * (i:ℝ) + 3)) / (κ * ((i:ℝ) + 1)) :=
          div_le_div₀ (by positivity) (hc i) (by positivity) hl
      _ ≤ (2 / κ) * (1 / ((i:ℝ) + 1) ^ 2) := by
          rw [div_div, div_mul_div_comm, div_le_div_iff₀ (by positivity) (by positivity)]
          nlinarith [hκ0]

/-- complex continuation of `lamTerm (n+1)` (the closed form of `lamTerm_closed`). -/
noncomputable def tC (z : ℝ) (ε : ℂ) (n : ℕ) : ℂ :=
  ((2 * z * z ^ (n + 1) : ℝ) : ℂ) * (-ε * ((n:ℂ) + 1 + ε) / (((n:ℂ) + 1) ^ 2 * ((n:ℂ) + 2))) *
    ∏ j ∈ Finset.range n, (1 - ε ^ 2 / ((j:ℂ) + 1) ^ 2)

lemma prod_bound {B : ℝ} {ε : ℂ} (hεB : ‖ε‖ ≤ B) (n : ℕ) :
    ‖∏ j ∈ Finset.range n, (1 - ε ^ 2 / ((j:ℂ) + 1) ^ 2)‖ ≤ Real.exp (2 * B ^ 2) := by
  have hB : 0 ≤ B := (norm_nonneg _).trans hεB
  refine (Finset.norm_prod_le _ _).trans ?_
  have h1 : ∀ j ∈ Finset.range n, ‖1 - ε ^ 2 / ((j:ℂ) + 1) ^ 2‖ ≤
      Real.exp (B ^ 2 * (1 / ((j:ℝ) + 1) ^ 2)) := by
    intro j _
    refine (norm_sub_le _ _).trans ?_
    rw [norm_one, norm_div, norm_pow, norm_pow,
      show ((j:ℂ) + 1) = ((j + 1 : ℕ) : ℂ) by push_cast; ring, Complex.norm_natCast]
    push_cast
    have : ‖ε‖ ^ 2 / ((j:ℝ) + 1) ^ 2 ≤ B ^ 2 * (1 / ((j:ℝ) + 1) ^ 2) := by
      rw [mul_one_div]
      gcongr
    linarith [Real.add_one_le_exp (B ^ 2 * (1 / ((j:ℝ) + 1) ^ 2))]
  refine (Finset.prod_le_prod₀ (fun j _ => norm_nonneg _) h1).trans ?_
  rw [← Real.exp_sum, ← Finset.mul_sum, mul_comm 2]
  gcongr
  exact hk_le n

lemma tC_bound {B z : ℝ} (hz : 0 ≤ z) (hz1 : z ≤ 1) {ε : ℂ} (hεB : ‖ε‖ ≤ B) (n : ℕ) :
    ‖tC z ε n‖ ≤ (2 * B * (1 + B) * Real.exp (2 * B ^ 2)) * z ^ (n + 1) := by
  have hB : 0 ≤ B := (norm_nonneg _).trans hεB
  have hn : (0:ℝ) ≤ n := Nat.cast_nonneg n
  unfold tC
  rw [norm_mul, norm_mul, Complex.norm_real, Real.norm_eq_abs,
    abs_of_nonneg (by positivity)]
  have hq : ‖-ε * ((n:ℂ) + 1 + ε) / (((n:ℂ) + 1) ^ 2 * ((n:ℂ) + 2))‖ ≤ B * (1 + B) := by
    rw [norm_div, norm_mul, norm_neg, norm_mul, norm_pow,
      show ((n:ℂ) + 1) = ((n + 1 : ℕ) : ℂ) by push_cast; ring,
      show ((n:ℂ) + 2) = ((n + 2 : ℕ) : ℂ) by push_cast; ring,
      Complex.norm_natCast, Complex.norm_natCast]
    have h2 : ‖(n:ℂ) + 1 + ε‖ ≤ ((n:ℝ) + 1) * (1 + B) := by
      refine (norm_add_le _ _).trans ?_
      rw [show ((n:ℂ) + 1) = ((n + 1 : ℕ) : ℂ) by push_cast; ring, Complex.norm_natCast]
      push_cast; nlinarith
    push_cast
    rw [div_le_iff₀ (by positivity)]
    calc ‖ε‖ * ‖(n:ℂ) + 1 + ε‖ ≤ B * (((n:ℝ) + 1) * (1 + B)) :=
          mul_le_mul hεB h2 (norm_nonneg _) hB
      _ ≤ B * (1 + B) * (((n:ℝ) + 1) ^ 2 * ((n:ℝ) + 2)) := by
          have : ((n:ℝ) + 1) ≤ ((n:ℝ) + 1) ^ 2 * ((n:ℝ) + 2) := by nlinarith
          have : 0 ≤ B * (1 + B) := by positivity
          nlinarith
  have hp := prod_bound hεB n
  calc 2 * z * z ^ (n + 1) * ‖-ε * ((n:ℂ) + 1 + ε) / (((n:ℂ) + 1) ^ 2 * ((n:ℂ) + 2))‖ *
        ‖∏ j ∈ Finset.range n, (1 - ε ^ 2 / ((j:ℂ) + 1) ^ 2)‖
      ≤ 2 * z * z ^ (n + 1) * (B * (1 + B)) * Real.exp (2 * B ^ 2) := by
        gcongr
    _ ≤ (2 * B * (1 + B) * Real.exp (2 * B ^ 2)) * z ^ (n + 1) := by
        have h1 : z * z ^ (n + 1) ≤ z ^ (n + 1) := mul_le_of_le_one_left (by positivity) hz1
        have h2 : 0 ≤ B * (1 + B) * Real.exp (2 * B ^ 2) := by positivity
        nlinarith

end ContA

section ContC
open Set

lemma zz_lt_one (θ : ℝ) (h0 : 0 ≤ θ) (h1 : θ < π) : Real.sin (θ / 2) ^ 2 < 1 := by
  have hc : 0 < Real.cos (θ / 2) := Real.cos_pos_of_mem_Ioo ⟨by linarith [Real.pi_pos], by linarith⟩
  nlinarith [Real.sin_sq_add_cos_sq (θ / 2)]

/-- `k_c(s) 1_{(0,1]}(s)`, as a complex function (its Mellin transform is the first part of
`Areg`). -/
noncomputable def kInd (θ : ℝ) : ℝ → ℂ :=
  Set.indicator (Ioc 0 1) (fun s => ((k (Real.cos θ) s : ℝ) : ℂ))

lemma k_abs_le (θ s : ℝ) (hs0 : 0 ≤ s) (hs1 : s ≤ 1) : |k (Real.cos θ) s| ≤ 4 := by
  have hc1 := Real.cos_le_one θ
  have hc2 := Real.neg_one_le_cos θ
  have hq : Real.sqrt (1 + 2 * Real.cos θ * s + s ^ 2) ≤ 2 :=
    Real.sqrt_le_iff.mpr ⟨by norm_num, by nlinarith⟩
  have hq0 := Real.sqrt_nonneg (1 + 2 * Real.cos θ * s + s ^ 2)
  unfold k
  rw [abs_le]; constructor <;> linarith

lemma kInd_integrable (θ : ℝ) : Integrable (kInd θ) := by
  unfold kInd
  refine IntegrableOn.integrable_indicator ?_ measurableSet_Ioc
  refine (ContinuousOn.integrableOn_Icc ?_).mono_set Ioc_subset_Icc_self
  apply Continuous.continuousOn
  unfold k
  fun_prop

lemma mellin_kInd_diff (θ : ℝ) {s : ℂ} (hs : 0 < s.re) :
    DifferentiableAt ℂ (_root_.mellin (kInd θ)) s := by
  refine mellin_differentiableAt_of_isBigO_rpow (a := s.re + 1) (b := 0)
    ((kInd_integrable θ).locallyIntegrable.locallyIntegrableOn _) ?_ (by linarith) ?_ hs
  · apply Asymptotics.IsBigO.of_bound 0
    filter_upwards [eventually_gt_atTop 1] with t ht
    have : t ∉ Ioc (0:ℝ) 1 := fun h => by linarith [h.2]
    simp [kInd, Set.indicator_of_notMem this]
  · apply Asymptotics.IsBigO.of_bound 4
    filter_upwards [Ioo_mem_nhdsGT (show (0:ℝ) < 1 by norm_num)] with t ht
    rw [neg_zero, Real.rpow_zero, norm_one, mul_one, kInd,
      Set.indicator_of_mem (show t ∈ Ioc (0:ℝ) 1 from ⟨ht.1, ht.2.le⟩), Complex.norm_real,
      Real.norm_eq_abs]
    exact k_abs_le θ t ht.1.le ht.2.le

lemma mellin_kInd_real (θ ε : ℝ) : _root_.mellin (kInd θ) (ε : ℂ) =
    ((∫ s in Ioc (0:ℝ) 1, k (Real.cos θ) s * s ^ (ε - 1) : ℝ) : ℂ) := by
  unfold _root_.mellin
  rw [← integral_complex_ofReal]
  have e : ∫ t in Ioi (0:ℝ), (t : ℂ) ^ ((ε : ℂ) - 1) • kInd θ t =
      ∫ t in Ioi (0:ℝ), (Ioc (0:ℝ) 1).indicator
        (fun t => ((k (Real.cos θ) t * t ^ (ε - 1) : ℝ) : ℂ)) t := by
    refine setIntegral_congr_fun measurableSet_Ioi (fun t ht => ?_)
    have ht0 : (0:ℝ) < t := ht
    by_cases hm : t ∈ Ioc (0:ℝ) 1
    · rw [kInd, Set.indicator_of_mem hm, Set.indicator_of_mem hm, smul_eq_mul,
        Complex.ofReal_mul, Complex.ofReal_cpow ht0.le]
      push_cast; ring
    · rw [kInd, Set.indicator_of_notMem hm, Set.indicator_of_notMem hm, smul_zero]
  rw [e, setIntegral_indicator measurableSet_Ioc,
    Set.inter_eq_right.mpr Ioc_subset_Ioi_self]

/-- `(-1)^i Λ_{i+1}`. -/
noncomputable def aN (θ : ℝ) (i : ℕ) : ℝ := (-1) ^ i * T13.LamInt (i + 1) θ

/-- The complex continuation of `A(ε) - (π/sin πε) Λ_ε`. -/
noncomputable def Phi (θ : ℝ) (ε : ℂ) : ℂ :=
  _root_.mellin (kInd θ) ε - ∑' i : ℕ, (aN θ i : ℂ) / (ε - 1 - i)
    - π / Complex.sin (π * ε) * (((2 * zz θ : ℝ) : ℂ) + ∑' n, tC (zz θ) ε n)

lemma Phi_real (θ ε : ℝ) (h0 : 0 ≤ θ) (h1 : θ < π) :
    Phi θ (ε : ℂ) = ((Areg θ ε - π / Real.sin (π * ε) * Lam ε θ : ℝ) : ℂ) := by
  have hz1 := zz_lt_one θ h0 h1
  have hsum := lamTerm_summable ε (Real.sin (θ / 2) ^ 2) (sq_nonneg _) hz1
  have hL : Lam ε θ = 2 * zz θ + ∑' n : ℕ, T13.lamTerm ε (zz θ) (n + 1) := by
    rw [Lam, hsum.tsum_eq_zero_add, lamTerm_zero]; rfl
  have ht : ∀ n : ℕ, tC (zz θ) (ε : ℂ) n = ((T13.lamTerm ε (zz θ) (n + 1) : ℝ) : ℂ) := by
    intro n
    rw [lamTerm_closed, tC, pp]
    push_cast
    ring
  have ha : ∀ i : ℕ, (aN θ i : ℂ) / ((ε : ℂ) - 1 - i) = ((aCoef θ ε i : ℝ) : ℂ) := by
    intro i; rw [aCoef, aN]; push_cast; ring
  rw [Phi, mellin_kInd_real, hL, Areg]
  simp only [ht, ha]
  rw [← Complex.ofReal_tsum, ← Complex.ofReal_tsum]
  push_cast
  ring

lemma sin_ne_zero_of_far {δ : ℝ} (hδ : 0 < δ) {w : ℂ} (hd : ∀ n : ℤ, δ ≤ ‖w - n‖) :
    Complex.sin (π * w) ≠ 0 := by
  intro h
  obtain ⟨n, hn⟩ := Complex.sin_eq_zero_iff.mp h
  have hw : w = n := by
    have hπ : (π : ℂ) ≠ 0 := Complex.ofReal_ne_zero.mpr Real.pi_ne_zero
    apply mul_left_cancel₀ hπ; rw [hn]; ring
  have := hd n
  rw [hw, sub_self, norm_zero] at this
  linarith

lemma Phi_diff (θ : ℝ) (h0 : 0 ≤ θ) (h1 : θ < π) {δ B : ℝ} (hδ : 0 < δ) (hB : 0 < B)
    {U : Set ℂ} (hU : IsOpen U)
    (hUp : ∀ z ∈ U, 0 < z.re ∧ ‖z‖ ≤ B ∧ ∀ n : ℤ, δ ≤ ‖z - n‖) :
    DifferentiableOn ℂ (Phi θ) U := by
  have hz1 := zz_lt_one θ h0 h1
  have hz0 : 0 ≤ zz θ := sq_nonneg _
  have hM : DifferentiableOn ℂ (_root_.mellin (kInd θ)) U := fun w hw =>
    (mellin_kInd_diff θ (hUp w hw).1).differentiableWithinAt
  have hA : DifferentiableOn ℂ (fun w : ℂ => ∑' i : ℕ, (aN θ i : ℂ) / (w - 1 - i)) U := by
    refine aSum_diff hδ hB hU (fun z hz => ⟨(hUp z hz).2.1, (hUp z hz).2.2⟩) (aN θ) (fun i => ?_)
    rw [aN, abs_mul, abs_pow, abs_neg, abs_one, one_pow, one_mul]
    exact lam_abs_le θ i
  have hS : DifferentiableOn ℂ (fun w : ℂ => π / Complex.sin (π * w)) U := by
    have hs : DifferentiableOn ℂ (fun w : ℂ => Complex.sin (π * w)) U := by fun_prop
    refine (differentiableOn_const (π:ℂ)).div hs (fun w hw => ?_)
    exact sin_ne_zero_of_far hδ (hUp w hw).2.2
  have hT : DifferentiableOn ℂ (fun w : ℂ => ∑' n, tC (zz θ) w n) U := by
    refine Complex.differentiableOn_tsum_of_summable_norm
      (F := fun (n : ℕ) (w : ℂ) => tC (zz θ) w n)
      (u := fun n : ℕ => (2 * B * (1 + B) * Real.exp (2 * B ^ 2)) * zz θ ^ (n + 1)) ?_ ?_ hU ?_
    · apply Summable.mul_left
      exact (summable_geometric_of_lt_one hz0 hz1).comp_injective (add_left_injective 1)
    · intro n
      apply Differentiable.differentiableOn
      unfold tC
      fun_prop
    · intro n w hw
      exact tC_bound hz0 hz1.le (hUp w hw).2.1 n
  unfold Phi
  exact (hM.sub hA).sub (hS.mul ((differentiableOn_const _).add hT))

/-- the identity `A(ε) = (π/sin πε) Λ_ε` for non-integer `ε > 1`, by analytic continuation
from `0 < ε < 1` (Areg_small) along a chain of three rectangles avoiding the integers. -/
lemma Areg_large (θ ε : ℝ) (h0 : 0 ≤ θ) (h1 : θ < π) (hε1 : 1 < ε) (hnat : ∀ j : ℕ, ε ≠ j) :
    Areg θ ε = π / Real.sin (π * ε) * Lam ε θ := by
  obtain ⟨U, δ, hδ, hUo, hUc, h12, hεU, hUp⟩ := domain ε hε1 hnat
  have hAn : AnalyticOnNhd ℂ (Phi θ) U :=
    (Phi_diff θ h0 h1 hδ (by linarith : (0:ℝ) < ε + 3) hUo hUp).analyticOnNhd hUo
  have hfreq : ∃ᶠ w in 𝓝[≠] (1/2 : ℂ), Phi θ w = 0 := by
    set u : ℕ → ℝ := fun n => 1/2 + 1 / ((n:ℝ) + 1) / 4 with hu
    have hlim : Tendsto u atTop (𝓝 (1/2)) := by
      have := (tendsto_one_div_add_atTop_nhds_zero_nat (𝕜 := ℝ)).div_const 4
      have := (tendsto_const_nhds (x := (1/2 : ℝ))).add this
      simpa [hu] using this
    have hlimC : Tendsto (fun n => ((u n : ℝ) : ℂ)) atTop (𝓝[≠] (1/2 : ℂ)) := by
      refine tendsto_nhdsWithin_iff.mpr ⟨?_, Eventually.of_forall (fun n => ?_)⟩
      · have := (Complex.continuous_ofReal.tendsto _).comp hlim
        simpa [Function.comp_def] using this
      · show ((u n : ℝ) : ℂ) ≠ 1/2
        intro h
        have h' : u n = 1/2 := by have := congrArg Complex.re h; simpa using this
        rw [hu] at h'
        simp only at h'
        have : 0 < 1 / ((n:ℝ) + 1) / 4 := by positivity
        linarith
    refine hlimC.frequently (Frequently.of_forall (fun n => ?_))
    have hpos : 0 < 1 / ((n:ℝ) + 1) / 4 := by positivity
    have hle : 1 / ((n:ℝ) + 1) / 4 ≤ 1/4 := by
      rw [div_le_iff₀ (by norm_num : (0:ℝ) < 4)]
      rw [div_le_iff₀ (by positivity)]
      nlinarith [(Nat.cast_nonneg n : (0:ℝ) ≤ n)]
    rw [Phi_real θ (u n) h0 h1, Areg_small θ (u n) h0 h1 (by rw [hu]; linarith)
      (by rw [hu]; linarith), sub_self, Complex.ofReal_zero]
  have := hAn.eqOn_of_preconnected_of_frequently_eq analyticOnNhd_const hUc h12 hfreq hεU
  rw [Phi_real θ ε h0 h1] at this
  have h' : Areg θ ε - π / Real.sin (π * ε) * Lam ε θ = 0 := Complex.ofReal_eq_zero.mp this
  linarith

end ContC

end Y02

theorem Y02.proof_1 : Y02.stmt_1 := by
  intro ν θ h0 h1
  apply Y02.lamTerm_summable ν _ (sq_nonneg _)
  have hc : 0 < Real.cos (θ / 2) := Real.cos_pos_of_mem_Ioo ⟨by linarith [Real.pi_pos], by linarith⟩
  nlinarith [Real.sin_sq_add_cos_sq (θ / 2)]

theorem Y02.proof_2 : Y02.stmt_2 := fun m hm θ _ _ => Y02.stmt2_aux m hm θ

theorem Y02.proof_3 : Y02.stmt_3 := by
  intro θ n _ _ hn
  obtain ⟨h1, h2⟩ := Y02.defect (Real.cos θ) n hn
  exact ⟨h1, by rw [Y02.L_eq θ n hn, h2]⟩
theorem Y02.proof_4 : Y02.stmt_4 := Y02.stmt4_aux
theorem Y02.proof_5 : Y02.stmt_5 := Y02.stmt5_aux
theorem Y02.proof_6 : Y02.stmt_6 := by
  refine Y02.stmt6_of_Areg (fun θ ε h0 h1 hε hnat => ?_)
  rcases lt_or_gt_of_ne (show ε ≠ 1 by exact_mod_cast hnat 1) with h | h
  · exact Y02.Areg_small θ ε h0 h1 hε h
  · exact Y02.Areg_large θ ε h0 h1 h hnat
theorem Y02.proof_7 : Y02.stmt_7 := fun θ _ _ => Y02.L_cont θ
theorem Y02.proof_8 : Y02.stmt_8 := Y02.stmt8_of_stmt6 Y02.proof_6

theorem Y02.proof_9 : Y02.stmt_9 := by
  intro θ h0 h1
  obtain ⟨z0, z1⟩ := Y02.zz_bounds θ h0 h1
  have hU := Y02.U_bounds z0 z1
  have hL := Y02.Li2_pos_lt z0 z1
  have hU' : Y02.U θ = (1 - Real.sin (θ / 2) ^ 2) / Real.sin (θ / 2) ^ 2
      * Real.log (1 - Real.sin (θ / 2) ^ 2) := rfl
  have ha2 : Y02.a2 θ = Real.pi ^ 2 / 6 + 1 / 2 - D00.Li2 (Real.sin (θ / 2) ^ 2) - Y02.U θ ^ 2 / 2 := rfl
  rw [← hU'] at hU
  refine ⟨hU.1, hU.2, ?_, ?_, ?_, ?_⟩
  · unfold Y02.a1; linarith
  · unfold Y02.a1; linarith
  · rw [ha2]; nlinarith
  · rw [ha2]; nlinarith

theorem Y02.proof_10 : Y02.stmt_10 := by
  intro θ h0 h1
  obtain ⟨z0, z1⟩ := Y02.zz_bounds θ h0 h1
  exact Y02.S12_identity z0 z1
theorem Y02.proof_11 : Y02.stmt_11 := Y02.stmt11_aux
theorem Y02.proof_12 : Y02.stmt_12 := Y02.stmt12_aux
theorem Y02.proof_13 : Y02.stmt_13 := Y02.stmt13_aux
theorem Y02.proof_14 : Y02.stmt_14 := Y02.stmt14_of_stmt13 Y02.proof_13
