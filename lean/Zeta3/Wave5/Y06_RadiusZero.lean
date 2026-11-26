import Mathlib
import Zeta3.Wave4.X01_RhombusDerivs
/-!
# Y06 -- radius zero at the diamond end (`thm:zero`)

Target id: Y06.
Paper label: `thm:zero` (paper/ieee/main.tex, lines 1237-1240, proof sketch 1242-1247).
Supplement: paper/ieee/proofs/sec_g2_zero.tex (Theorem 10, Lemmas `g2:lem:z:laplace`,
`g2:lem:z:coef`, `g2:lem:z:H`, `g2:lem:z:rest`, and "Proof of Theorem~10").

Quoted verbatim (main.tex 1237-1240):

  \begin{theorem}[Radius zero$^\dagger$]\label{thm:zero}
  For every rhombus $|x/A|^p+|y/B|^p=1$, let $I(\eta)=P/(4r(1+\eta))$, with $r=\sqrt{A^2+B^2}$
  and $\eta=\frac1p-1$, and let $I_n=[\eta^n]I(\eta)=I^{(n)}(0)/n!$ be its Taylor coefficients.
  Then $\limsup_{n\to\infty}|I_n/n!|^{1/n}=1/\rho$ with $\rho=\sqrt{\log^2(B/A)+\pi^2/4}$.
  In particular the Taylor series of $P$ at $p=1$, in $\eta$ or in $d$, has radius of
  convergence $0$: $P$ is infinitely differentiable but not analytic at $p=1$.
  \end{theorem}

Quoted verbatim (sec_g2_zero.tex, the supplement's lemmas):

  R_1(z)=A^2+B^2\eu^{-2z},\quad R_2(z)=B^2+A^2\eu^{-2z},
  h_i=\frac{\sqrt{R_i}}{r},\quad H=h_1+h_2,\quad K(\eta,u)=(1+\eu^{-u})^{-2-\eta}
  ... the Taylor coefficients of $K$ in $\eta$ are $k_j(u)=(1+\eu^{-u})^{-2}(-\ell_u)^j/j!$.

  \begin{lemma}\label{g2:lem:z:laplace}
  For $-1<\eta<1$, $I(\eta)=\int_0^\infty\eu^{-u}H(\eta u)K(\eta,u)\,du$. Also
  $|\partial_\eta^jK(\eta,u)|=K\ell_u^j\le1$ for $\eta\ge-2$, $u\ge0$; $|k_j(u)|\le\eu^{-ju}/j!$
  for $j\ge1$; and $0\le1-k_0(u)\le2\eu^{-u}$.
  \end{lemma}

  Write $H(z)=\sum_mH_mz^m$ ... and $\eta_D(n)=\sum_{m\ge0}(-1)^m(m+1)^{-n}$. Put
  $M_{n,j}=\int_0^\infty\eu^{-u}u^{n-j}k_j(u)\,du$.
  \begin{lemma}\label{g2:lem:z:coef}
  $I\in C^\infty(-1,1)$ and $I_n=\sum_{j=0}^nH_{n-j}M_{n,j}$ for every $n\ge0$. Moreover
  $M_{n,0}=n!\,\eta_D(n)$ for $n\ge2$, and $|M_{n,j}|\le\dfrac{(n-j)!}{j!\,(j+1)^{n-j+1}}$ for
  $j\ge1$.
  \end{lemma}

  \begin{lemma}\label{g2:lem:z:H}
  (a) The zeros of $R_1$ are the points $L+i\pi(k+\frac12)$, $k\in\ZZ$; those of $R_2$ are
  $-L+i\pi(k+\frac12)$. All are simple. ...
  (b) ... In particular $H$ is holomorphic on $D(0,\rho)$, and for every $\sigma\in(0,\rho)$,
  $|H_m|\le C_\sigma\sigma^{-m}$ with $C_\sigma=\max_{|z|=\sigma}|H|$.
  (c) The Taylor series of $H$ at $0$ has radius of convergence exactly $\rho$.
  \end{lemma}

  For $n\ge1$ put $E_n=\dfrac{I_n}{n!}-H_n$.
  \begin{lemma}\label{g2:lem:z:rest}
  For every $\sigma\in(0,\rho)$ there is $C'_\sigma$ such that $|E_n|\le C'_\sigma(2\sigma)^{-n}$
  for all $n\ge2$. ...
  \end{lemma}

## Formalization choices

* The perimeter is `X01.P A B p = (4/p) ∫_0^1 √(A^2 t^(2/p-2) + B^2 (1-t)^(2/p-2)) dt`
  (imported from `Zeta3.Wave4.X01_RhombusDerivs`; `X01.stmt_3` there identifies it with the true
  length of the curve), and `r = X01.r A B = √(A^2+B^2)`.
* The paper's `I(η) = P/(4r(1+η))` is `Ieta A B η = X01.P A B (1/(1+η)) / (4 r (1+η))`. This is
  NOT `X01.I` (which is `P/(4r)` without the factor `1+η`); it is a new object, so it gets its
  own name. It is only used on `(-1,1)`.
* `I_n = I^{(n)}(0)/n!` is `In A B n = iteratedDeriv n (Ieta A B) 0 / n!`; it is genuine because
  stmt_3 asserts that `Ieta A B` is `C^∞` on `(-1,1)`. Smoothness is `ContDiffOn ℝ (⊤ : ℕ∞)`,
  i.e. `C^∞` (NOT the `ω`/analytic level `⊤ : WithTop ℕ∞`, which would be false at `η = 0`).
* `H` is used on the real line only, as the real function
  `H(x) = (√(A^2+B^2e^{-2x}) + √(B^2+A^2e^{-2x}))/r` (the arguments are positive, so these are the
  principal roots of the paper for real argument). Its Taylor coefficients are
  `Hc A B m = H^{(m)}(0)/m!`; `H` is smooth (indeed real analytic), so they are genuine. The
  holomorphic continuation of `H` to `D(0,ρ)` is not stated as such; its two consequences used in
  the proof are stated instead: the Cauchy bound `|H_m| ≤ C_σ σ^{-m}` (stmt_6) and radius
  exactly `ρ` in Cauchy--Hadamard form `limsup |H_m|^{1/m} = 1/ρ` (stmt_7).
* The zeros of `R_1`, `R_2` (part (a)) are stated over `ℂ` with `Complex.exp`.
* `M_{n,j}` is a Bochner integral over `(0,∞)`; its integrability is asserted (stmt_4), so
  the junk value `0` cannot be used.
* `η_D(n) = ∑ (-1)^m (m+1)^{-n}` is stated through `HasSum`.
* `limsup` is `Filter.limsup` along `atTop` of the real sequence `|I_n/n!|^{1/n}` (`Real.rpow`
  with exponent `1/n`; the base is `≥ 0`). Real `limsup` of a sequence that is not bounded above
  is a junk value; to rule this out the final statements also assert that the sequence is
  bounded above eventually (`IsBoundedUnder (· ≤ ·) atTop`). (The sequence is `≥ 0`, so it is
  automatically bounded below.)
* `ζ(3)` does not occur in this target.
* "Radius of convergence 0" is stated as: for every `η ≠ 0` (resp. `d ≠ 0`) the Taylor series
  `∑ c_n x^n` is not `Summable` (stmt_10, stmt_11). "Not analytic" is `¬ AnalyticAt ℝ`.

Statements:
* stmt_1: Laplace form `I(η) = ∫_0^∞ e^{-u} H(ηu) K(η,u) du` (integrable), `-1 < η < 1`.
* stmt_2: the bounds on `K`: `|∂_η^j K| = K ℓ_u^j ≤ 1` (as the closed form of the `j`-th
  derivative, for `η ≥ -2`, `u ≥ 0`), `|k_j(u)| ≤ e^{-ju}/j!` (`j ≥ 1`), `0 ≤ 1 - k_0 ≤ 2e^{-u}`.
* stmt_3: `I` is `C^∞` on `(-1,1)`.
* stmt_4: `I_n = ∑_{j=0}^n H_{n-j} M_{n,j}`, the `M_{n,j}` integrals exist,
  `M_{n,0} = n! η_D(n)` for `n ≥ 2`, and `|M_{n,j}| ≤ (n-j)!/(j!(j+1)^{n-j+1})` for `1 ≤ j ≤ n`.
* stmt_5: the zeros of `R_1` are `L + iπ(k+1/2)`, those of `R_2` are `-L + iπ(k+1/2)`, all simple.
* stmt_6: Cauchy bound: for every `σ ∈ (0,ρ)` there is `C` with `|H_m| ≤ C σ^{-m}` for all `m`.
* stmt_7: the Taylor series of `H` at `0` has radius exactly `ρ`: `limsup |H_m|^{1/m} = 1/ρ`.
* stmt_8: the rest: for every `σ ∈ (0,ρ)` there is `C'` with `|I_n/n! - H_n| ≤ C' (2σ)^{-n}`,
  `n ≥ 2`.
* stmt_9: `thm:zero`, main assertion: `limsup |I_n/n!|^{1/n} = 1/ρ` (with boundedness).
* stmt_10: corollary: `∑ I_n η^n` diverges (is not summable) for every `η ≠ 0`.
* stmt_11: `thm:zero` in full: the limsup identity, `I` and `P` are `C^∞` at `η = 0` resp.
  `p = 1`, the Taylor series of `P` at `p = 1` in `η` (of `η ↦ P(1/(1+η))` at `0`) and in `d`
  (of `P` at `1`) diverge at every nonzero point, and `P` is not analytic at `p = 1`.

## NOT COVERED

* The holomorphic continuation of `h_1`, `h_2`, `H` to the slit disk `Ω` and the description of
  the zeros inside `D(0,ρ_2)` (only the zero sets, simplicity, and the two consequences
  stmt_6/stmt_7 are stated).
* The explicit constants `C_σ = max_{|z|=σ}|H|` and
  `C'_σ = C_σ[1 + ½∑_{j≥1}(2σ)^j/(j!)^2]` (only existence of constants is stated).
* The bound `|H(z)| ≤ (A+B)(1+e^{-Re z})/r` on the strip and the Cauchy estimates on discs of
  radius `π/4` (used inside the proofs only).
* Everything after the theorem in the paper (the asymptotic law for `f_n` at `A = B`, Borel
  summability, Padé--Borel numerics) and Theorem `thm:one`.
-/
-- BEGIN STATEMENTS
namespace Y06

open MeasureTheory Set Filter

/-- `L = log(B/A)`. -/
noncomputable def L (A B : ℝ) : ℝ := Real.log (B / A)

/-- `ρ = √(log^2(B/A) + π^2/4)`. -/
noncomputable def ρ (A B : ℝ) : ℝ := Real.sqrt (L A B ^ 2 + Real.pi ^ 2 / 4)

/-- The paper's `I(η) = P/(4r(1+η))` with `p = 1/(1+η)` (not `X01.I`, which lacks `1+η`). -/
noncomputable def Ieta (A B η : ℝ) : ℝ :=
  X01.P A B (1 / (1 + η)) / (4 * X01.r A B * (1 + η))

/-- The Taylor coefficients `I_n = I^{(n)}(0)/n!`. -/
noncomputable def In (A B : ℝ) (n : ℕ) : ℝ :=
  iteratedDeriv n (Ieta A B) 0 / (n.factorial : ℝ)

/-- `H(x) = (√(A^2 + B^2 e^{-2x}) + √(B^2 + A^2 e^{-2x}))/r` on the real line. -/
noncomputable def H (A B x : ℝ) : ℝ :=
  (Real.sqrt (A ^ 2 + B ^ 2 * Real.exp (-2 * x)) +
    Real.sqrt (B ^ 2 + A ^ 2 * Real.exp (-2 * x))) / X01.r A B

/-- The Taylor coefficients `H_m = H^{(m)}(0)/m!` of `H` at `0`. -/
noncomputable def Hc (A B : ℝ) (m : ℕ) : ℝ :=
  iteratedDeriv m (H A B) 0 / (m.factorial : ℝ)

/-- `K(η,u) = (1+e^{-u})^{-2-η}`. -/
noncomputable def K (η u : ℝ) : ℝ := (1 + Real.exp (-u)) ^ (-2 - η)

/-- `ℓ_u = log(1+e^{-u})`. -/
noncomputable def ell (u : ℝ) : ℝ := Real.log (1 + Real.exp (-u))

/-- `k_j(u) = (1+e^{-u})^{-2} (-ℓ_u)^j / j!`, the Taylor coefficients of `K` in `η`. -/
noncomputable def kj (j : ℕ) (u : ℝ) : ℝ :=
  ((1 + Real.exp (-u)) ^ 2)⁻¹ * (-ell u) ^ j / (j.factorial : ℝ)

/-- `M_{n,j} = ∫_0^∞ e^{-u} u^{n-j} k_j(u) du` (used for `j ≤ n`). -/
noncomputable def Mnj (n j : ℕ) : ℝ :=
  ∫ u in Ioi (0 : ℝ), Real.exp (-u) * u ^ (n - j) * kj j u

/-- `R_1(z) = A^2 + B^2 e^{-2z}` over `ℂ`. -/
noncomputable def R1 (A B : ℝ) (z : ℂ) : ℂ := (A : ℂ) ^ 2 + (B : ℂ) ^ 2 * Complex.exp (-2 * z)

/-- `R_2(z) = B^2 + A^2 e^{-2z}` over `ℂ`. -/
noncomputable def R2 (A B : ℝ) (z : ℂ) : ℂ := (B : ℂ) ^ 2 + (A : ℂ) ^ 2 * Complex.exp (-2 * z)

/-- `g2:lem:z:laplace`: "For $-1<\eta<1$, $I(\eta)=\int_0^\infty\eu^{-u}H(\eta u)K(\eta,u)\,du$"
(the integrand is integrable on `(0,∞)`). -/
def stmt_1 : Prop :=
  ∀ A B : ℝ, 0 < A → 0 < B → ∀ η : ℝ, -1 < η → η < 1 →
    IntegrableOn (fun u : ℝ => Real.exp (-u) * H A B (η * u) * K η u) (Ioi 0) ∧
    Ieta A B η = ∫ u in Ioi (0 : ℝ), Real.exp (-u) * H A B (η * u) * K η u

/-- `g2:lem:z:laplace`, the bounds: "$|\partial_\eta^jK(\eta,u)|=K\ell_u^j\le1$ for $\eta\ge-2$,
$u\ge0$; $|k_j(u)|\le\eu^{-ju}/j!$ for $j\ge1$; and $0\le1-k_0(u)\le2\eu^{-u}$". -/
def stmt_2 : Prop :=
  (∀ u : ℝ, 0 ≤ u → ∀ j : ℕ,
    iteratedDeriv j (fun η => K η u) = (fun η => (-ell u) ^ j * K η u) ∧
    ∀ η : ℝ, -2 ≤ η → |(-ell u) ^ j * K η u| = K η u * ell u ^ j ∧ K η u * ell u ^ j ≤ 1) ∧
  (∀ u : ℝ, 0 ≤ u → ∀ j : ℕ, 1 ≤ j → |kj j u| ≤ Real.exp (-(j * u)) / (j.factorial : ℝ)) ∧
  (∀ u : ℝ, 0 ≤ u → 0 ≤ 1 - kj 0 u ∧ 1 - kj 0 u ≤ 2 * Real.exp (-u))

/-- `g2:lem:z:coef`: "$I\in C^\infty(-1,1)$". -/
def stmt_3 : Prop :=
  ∀ A B : ℝ, 0 < A → 0 < B → ContDiffOn ℝ (⊤ : ℕ∞) (Ieta A B) (Ioo (-1) 1)

/-- `g2:lem:z:coef`: "$I_n=\sum_{j=0}^nH_{n-j}M_{n,j}$ for every $n\ge0$. Moreover
$M_{n,0}=n!\,\eta_D(n)$ for $n\ge2$, and $|M_{n,j}|\le\dfrac{(n-j)!}{j!\,(j+1)^{n-j+1}}$ for
$j\ge1$" (with the integrals `M_{n,j}` existing; `η_D(n)` as a `HasSum`). -/
def stmt_4 : Prop :=
  (∀ n j : ℕ, j ≤ n →
    IntegrableOn (fun u : ℝ => Real.exp (-u) * u ^ (n - j) * kj j u) (Ioi 0)) ∧
  (∀ A B : ℝ, 0 < A → 0 < B → ∀ n : ℕ,
    In A B n = ∑ j ∈ Finset.range (n + 1), Hc A B (n - j) * Mnj n j) ∧
  (∀ n : ℕ, 2 ≤ n →
    HasSum (fun m : ℕ => (-1 : ℝ) ^ m / ((m : ℝ) + 1) ^ n) (Mnj n 0 / (n.factorial : ℝ))) ∧
  (∀ n j : ℕ, 1 ≤ j → j ≤ n →
    |Mnj n j| ≤ ((n - j).factorial : ℝ) / ((j.factorial : ℝ) * ((j : ℝ) + 1) ^ (n - j + 1)))

/-- `g2:lem:z:H`(a): "The zeros of $R_1$ are the points $L+i\pi(k+\frac12)$, $k\in\ZZ$; those of
$R_2$ are $-L+i\pi(k+\frac12)$. All are simple." -/
def stmt_5 : Prop :=
  ∀ A B : ℝ, 0 < A → 0 < B →
    (∀ z : ℂ, R1 A B z = 0 ↔
      ∃ k : ℤ, z = (L A B : ℂ) + Complex.I * Real.pi * ((k : ℂ) + 1 / 2)) ∧
    (∀ z : ℂ, R2 A B z = 0 ↔
      ∃ k : ℤ, z = -(L A B : ℂ) + Complex.I * Real.pi * ((k : ℂ) + 1 / 2)) ∧
    (∀ z : ℂ, R1 A B z = 0 → deriv (R1 A B) z ≠ 0) ∧
    (∀ z : ℂ, R2 A B z = 0 → deriv (R2 A B) z ≠ 0)

/-- `g2:lem:z:H`(b): "for every $\sigma\in(0,\rho)$, $|H_m|\le C_\sigma\sigma^{-m}$"
(existence of the constant). -/
def stmt_6 : Prop :=
  ∀ A B : ℝ, 0 < A → 0 < B → ∀ σ : ℝ, 0 < σ → σ < ρ A B →
    ∃ C : ℝ, ∀ m : ℕ, |Hc A B m| ≤ C * σ⁻¹ ^ m

/-- `g2:lem:z:H`(c): "The Taylor series of $H$ at $0$ has radius of convergence exactly $\rho$",
in Cauchy--Hadamard form `limsup |H_m|^{1/m} = 1/ρ` (with the sequence eventually bounded). -/
def stmt_7 : Prop :=
  ∀ A B : ℝ, 0 < A → 0 < B →
    IsBoundedUnder (· ≤ ·) atTop (fun m : ℕ => |Hc A B m| ^ (1 / (m : ℝ))) ∧
    limsup (fun m : ℕ => |Hc A B m| ^ (1 / (m : ℝ))) atTop = 1 / ρ A B

/-- `g2:lem:z:rest`: "For every $\sigma\in(0,\rho)$ there is $C'_\sigma$ such that
$|E_n|\le C'_\sigma(2\sigma)^{-n}$ for all $n\ge2$", `E_n = I_n/n! - H_n`. -/
def stmt_8 : Prop :=
  ∀ A B : ℝ, 0 < A → 0 < B → ∀ σ : ℝ, 0 < σ → σ < ρ A B →
    ∃ C : ℝ, ∀ n : ℕ, 2 ≤ n →
      |In A B n / (n.factorial : ℝ) - Hc A B n| ≤ C * (2 * σ)⁻¹ ^ n

/-- `thm:zero`, main assertion: "$\limsup_{n\to\infty}|I_n/n!|^{1/n}=1/\rho$ with
$\rho=\sqrt{\log^2(B/A)+\pi^2/4}$" (with the sequence eventually bounded above). -/
def stmt_9 : Prop :=
  ∀ A B : ℝ, 0 < A → 0 < B →
    IsBoundedUnder (· ≤ ·) atTop (fun n : ℕ => |In A B n / (n.factorial : ℝ)| ^ (1 / (n : ℝ))) ∧
    limsup (fun n : ℕ => |In A B n / (n.factorial : ℝ)| ^ (1 / (n : ℝ))) atTop = 1 / ρ A B

/-- "the Taylor series ... in $\eta$ ... has radius of convergence $0$", for `I`: the series
`∑ I_n η^n` diverges (is not summable) at every `η ≠ 0`. -/
def stmt_10 : Prop :=
  ∀ A B : ℝ, 0 < A → 0 < B → ∀ η : ℝ, η ≠ 0 →
    ¬ Summable (fun n : ℕ => In A B n * η ^ n)

/-- `thm:zero` in full: the limsup identity; "the Taylor series of $P$ at $p=1$, in $\eta$ or in
$d$, has radius of convergence $0$"; "$P$ is infinitely differentiable but not analytic at
$p=1$". -/
def stmt_11 : Prop :=
  ∀ A B : ℝ, 0 < A → 0 < B →
    (IsBoundedUnder (· ≤ ·) atTop
        (fun n : ℕ => |In A B n / (n.factorial : ℝ)| ^ (1 / (n : ℝ))) ∧
      limsup (fun n : ℕ => |In A B n / (n.factorial : ℝ)| ^ (1 / (n : ℝ))) atTop = 1 / ρ A B) ∧
    ContDiffAt ℝ (⊤ : ℕ∞) (Ieta A B) 0 ∧
    (∀ η : ℝ, η ≠ 0 → ¬ Summable (fun n : ℕ => In A B n * η ^ n)) ∧
    ContDiffAt ℝ (⊤ : ℕ∞) (fun η : ℝ => X01.P A B (1 / (1 + η))) 0 ∧
    (∀ η : ℝ, η ≠ 0 → ¬ Summable (fun n : ℕ =>
      iteratedDeriv n (fun η : ℝ => X01.P A B (1 / (1 + η))) 0 / (n.factorial : ℝ) * η ^ n)) ∧
    ContDiffAt ℝ (⊤ : ℕ∞) (X01.P A B) 1 ∧
    (∀ d : ℝ, d ≠ 0 → ¬ Summable (fun n : ℕ =>
      iteratedDeriv n (X01.P A B) 1 / (n.factorial : ℝ) * d ^ n)) ∧
    ¬ AnalyticAt ℝ (X01.P A B) 1

end Y06
-- END STATEMENTS

namespace Y06
open MeasureTheory Set Filter

lemma ell_nonneg (u : ℝ) : 0 ≤ ell u :=
  Real.log_nonneg (by linarith [Real.exp_pos (-u)])

lemma ell_le_exp (u : ℝ) : ell u ≤ Real.exp (-u) := by
  have := Real.log_le_sub_one_of_pos (by positivity : (0:ℝ) < 1 + Real.exp (-u))
  unfold ell; linarith

lemma exp_neg_le_one {u : ℝ} (hu : 0 ≤ u) : Real.exp (-u) ≤ 1 :=
  Real.exp_le_one_iff.mpr (by linarith)

lemma ell_le_one {u : ℝ} (hu : 0 ≤ u) : ell u ≤ 1 :=
  (ell_le_exp u).trans (exp_neg_le_one hu)

lemma K_eq (η u : ℝ) : K η u = Real.exp (ell u * (-2 - η)) := by
  unfold K ell; rw [Real.rpow_def_of_pos (by positivity)]

lemma K_pos (η u : ℝ) : 0 < K η u := by rw [K_eq]; exact Real.exp_pos _

lemma K_le_one {η u : ℝ} (hη : -2 ≤ η) : K η u ≤ 1 := by
  rw [K_eq, Real.exp_le_one_iff]
  exact mul_nonpos_of_nonneg_of_nonpos (ell_nonneg u) (by linarith)

lemma hasDerivAt_K (u η : ℝ) : HasDerivAt (fun η => K η u) (-ell u * K η u) η := by
  simp_rw [K_eq]
  have h : HasDerivAt (fun η => ell u * (-2 - η)) (ell u * (-1)) η :=
    ((hasDerivAt_id η).const_sub (-2)).const_mul (ell u)
  convert h.exp using 1; ring

lemma iter_K (u : ℝ) : ∀ j : ℕ,
    iteratedDeriv j (fun η => K η u) = (fun η => (-ell u) ^ j * K η u) := by
  intro j
  induction j with
  | zero => funext η; simp
  | succ n ih =>
    rw [iteratedDeriv_succ, ih]; funext η
    rw [((hasDerivAt_K u η).const_mul ((-ell u) ^ n)).deriv]; ring

lemma kj_abs_le {u : ℝ} (_hu : 0 ≤ u) (j : ℕ) : |kj j u| ≤ Real.exp (-(j * u)) / (j.factorial : ℝ) := by
  have hx := Real.exp_pos (-u)
  have h1 : ((1 + Real.exp (-u)) ^ 2)⁻¹ ≤ 1 := inv_le_one_of_one_le₀ (one_le_pow₀ (by linarith))
  have h0 : 0 ≤ ((1 + Real.exp (-u)) ^ 2)⁻¹ := by positivity
  unfold kj
  rw [abs_div, abs_mul, abs_pow, abs_neg, abs_of_nonneg (ell_nonneg u), abs_of_nonneg h0,
    Nat.abs_cast]
  gcongr
  calc ((1 + Real.exp (-u)) ^ 2)⁻¹ * ell u ^ j ≤ 1 * ell u ^ j := by
        gcongr; exact pow_nonneg (ell_nonneg u) _
    _ ≤ Real.exp (-u) ^ j := by rw [one_mul]; gcongr; exacts [ell_nonneg u, ell_le_exp u]
    _ = Real.exp (-(j * u)) := by rw [← Real.exp_nat_mul]; ring_nf

lemma kj_abs_le_one {u : ℝ} (hu : 0 ≤ u) (j : ℕ) : |kj j u| ≤ 1 := by
  refine (kj_abs_le hu j).trans ?_
  rw [div_le_one (by positivity)]
  calc Real.exp (-(j * u)) ≤ 1 := Real.exp_le_one_iff.mpr (by nlinarith [(j.cast_nonneg : (0:ℝ) ≤ j)])
    _ ≤ _ := by exact_mod_cast Nat.one_le_iff_ne_zero.mpr (Nat.factorial_ne_zero j)

end Y06

theorem Y06.proof_2 : Y06.stmt_2 := by
  refine ⟨fun u hu j => ⟨Y06.iter_K u j, fun η hη => ⟨?_, ?_⟩⟩, fun u hu j _ => Y06.kj_abs_le hu j,
    fun u hu => ⟨?_, ?_⟩⟩
  · rw [abs_mul, abs_pow, abs_neg, abs_of_nonneg (Y06.ell_nonneg u),
      abs_of_pos (Y06.K_pos η u)]; ring
  · exact (mul_le_mul (Y06.K_le_one hη) (pow_le_one₀ (Y06.ell_nonneg u) (Y06.ell_le_one hu))
      (pow_nonneg (Y06.ell_nonneg u) _) zero_le_one).trans (by norm_num)
  · simp only [Y06.kj, pow_zero, mul_one, Nat.factorial_zero, Nat.cast_one, div_one, sub_nonneg]
    exact inv_le_one_of_one_le₀ (one_le_pow₀ (by linarith [Real.exp_pos (-u)]))
  · simp only [Y06.kj, pow_zero, mul_one, Nat.factorial_zero, Nat.cast_one, div_one]
    have hx := Real.exp_pos (-u)
    set x := Real.exp (-u)
    have : 1 - 2 * x ≤ ((1 + x) ^ 2)⁻¹ := by
      rw [inv_eq_one_div, le_div_iff₀ (by positivity)]; nlinarith [sq_nonneg x, hx.le]
    linarith

namespace Y06
open MeasureTheory Set Filter

lemma exp_L (A B : ℝ) (hA : 0 < A) (hB : 0 < B) : Complex.exp (L A B : ℂ) = ((B / A : ℝ) : ℂ) := by
  rw [← Complex.ofReal_exp, L, Real.exp_log (div_pos hB hA)]

lemma R1_zero_iff (A B : ℝ) (hA : 0 < A) (hB : 0 < B) (z : ℂ) :
    R1 A B z = 0 ↔ ∃ k : ℤ, z = (L A B : ℂ) + Complex.I * Real.pi * ((k : ℂ) + 1 / 2) := by
  have hA' : (A : ℂ) ≠ 0 := by exact_mod_cast hA.ne'
  have hB' : (B : ℂ) ≠ 0 := by exact_mod_cast hB.ne'
  have hc : Complex.exp (-2 * (L A B : ℂ) + Real.pi * Complex.I) * (B : ℂ) ^ 2 = -(A : ℂ) ^ 2 := by
    rw [Complex.exp_add, Complex.exp_pi_mul_I,
      show -2 * (L A B : ℂ) = -(L A B : ℂ) + -(L A B : ℂ) by ring, Complex.exp_add,
      Complex.exp_neg, exp_L A B hA hB]
    push_cast; field_simp
  have key : R1 A B z = 0 ↔ Complex.exp (-2 * z) = Complex.exp (-2 * (L A B : ℂ) + Real.pi * Complex.I) := by
    unfold R1
    constructor
    · intro h
      have h2 : Complex.exp (-2 * z) * (B : ℂ) ^ 2 =
          Complex.exp (-2 * (L A B : ℂ) + Real.pi * Complex.I) * (B : ℂ) ^ 2 := by
        rw [hc]; linear_combination h
      exact mul_right_cancel₀ (pow_ne_zero 2 hB') h2
    · intro h; rw [h]; linear_combination hc
  rw [key, Complex.exp_eq_exp_iff_exists_int]
  constructor
  · rintro ⟨n, hn⟩
    refine ⟨-n - 1, ?_⟩; push_cast; linear_combination (-1/2 : ℂ) * hn
  · rintro ⟨k, hk⟩
    refine ⟨-k - 1, ?_⟩; push_cast; linear_combination (-2 : ℂ) * hk

lemma hasDerivAt_R1 (A B : ℝ) (z : ℂ) :
    HasDerivAt (R1 A B) ((B : ℂ) ^ 2 * (Complex.exp (-2 * z) * (-2))) z := by
  have h := (((hasDerivAt_id z).const_mul (-2 : ℂ)).cexp.const_mul ((B : ℂ) ^ 2)).const_add ((A : ℂ) ^ 2)
  convert h using 1
  · rfl
  · simp

lemma R1_simple (A B : ℝ) (hA : 0 < A) (z : ℂ) (hz : R1 A B z = 0) : deriv (R1 A B) z ≠ 0 := by
  rw [(hasDerivAt_R1 A B z).deriv]
  have hA' : (A : ℂ) ≠ 0 := by exact_mod_cast hA.ne'
  unfold R1 at hz
  have : (B : ℂ) ^ 2 * (Complex.exp (-2 * z) * (-2)) = 2 * (A : ℂ) ^ 2 := by
    linear_combination (-2 : ℂ) * hz
  rw [this]; exact mul_ne_zero two_ne_zero (pow_ne_zero 2 hA')

lemma R2_eq (A B : ℝ) : R2 A B = R1 B A := rfl

lemma L_swap (A B : ℝ) : L B A = -L A B := by
  unfold L; rw [← Real.log_inv, inv_div]

end Y06

theorem Y06.proof_5 : Y06.stmt_5 := by
  intro A B hA hB
  refine ⟨Y06.R1_zero_iff A B hA hB, fun z => ?_, Y06.R1_simple A B hA, ?_⟩
  · rw [Y06.R2_eq, Y06.R1_zero_iff B A hB hA, Y06.L_swap]; push_cast; rfl
  · rw [Y06.R2_eq]; exact Y06.R1_simple B A hB

namespace Y06
open MeasureTheory Set Filter

noncomputable def gg (u : ℝ) : ℝ := 1 / (1 + Real.exp (-u))
noncomputable def gd (u : ℝ) : ℝ := Real.exp (-u) / (1 + Real.exp (-u)) ^ 2
noncomputable def F (A B η t : ℝ) : ℝ :=
  Real.sqrt (A ^ 2 * t ^ (2 * η) + B ^ 2 * (1 - t) ^ (2 * η))

lemma hasDerivAt_gg (u : ℝ) : HasDerivAt gg (gd u) u := by
  have h1 : HasDerivAt (fun u => 1 + Real.exp (-u)) (Real.exp (-u) * -1) u :=
    ((Real.hasDerivAt_exp (-u)).comp u (hasDerivAt_neg u)).const_add 1
  have h2 := h1.inv (by positivity)
  unfold gg gd; simp only [one_div]
  convert h2 using 1; ring

lemma gg_strictMono : StrictMono gg := by
  intro u v huv
  unfold gg
  apply one_div_lt_one_div_of_lt (by positivity)
  have := Real.exp_lt_exp.mpr (neg_lt_neg huv); linarith

lemma gg_image : gg '' Ioi 0 = Ioo (1 / 2) 1 := by
  ext t; constructor
  · rintro ⟨u, hu, rfl⟩
    have hx := Real.exp_pos (-u)
    have hx1 : Real.exp (-u) < 1 := by
      have := Real.exp_lt_exp.mpr (show -u < 0 by simp at hu; linarith); simpa using this
    unfold gg
    refine ⟨one_div_lt_one_div_of_lt (by positivity) (by linarith), ?_⟩
    rw [div_lt_one (by positivity)]; linarith
  · rintro ⟨h1, h2⟩
    have ht0 : 0 < t := by linarith
    have ht1 : 0 < 1 - t := by linarith
    refine ⟨Real.log (t / (1 - t)), ?_, ?_⟩
    · rw [Set.mem_Ioi]
      apply Real.log_pos; rw [one_lt_div ht1]; linarith
    · unfold gg
      rw [Real.exp_neg, Real.exp_log (div_pos ht0 ht1)]
      field_simp; ring

lemma gd_pos (u : ℝ) : 0 < gd u := by unfold gd; positivity

lemma exp_ell (u : ℝ) : Real.exp (ell u) = 1 + Real.exp (-u) := by
  unfold ell; exact Real.exp_log (by positivity)

lemma F_gg (A B η u : ℝ) : F A B η (gg u) =
    Real.exp (-(ell u * η)) * Real.sqrt (A ^ 2 + B ^ 2 * Real.exp (-2 * (η * u))) := by
  have hg : gg u = Real.exp (-ell u) := by
    rw [Real.exp_neg, exp_ell]; unfold gg; rw [one_div]
  have hg1 : 1 - gg u = Real.exp (-u - ell u) := by
    rw [Real.exp_sub, exp_ell]; unfold gg; field_simp; ring
  have r1 : Real.exp (-ell u) ^ (2 * η) = Real.exp (-(ell u * η)) ^ 2 := by
    rw [← Real.exp_mul, ← Real.exp_nat_mul]; ring_nf
  have r2 : Real.exp (-u - ell u) ^ (2 * η) =
      Real.exp (-(ell u * η)) ^ 2 * Real.exp (-2 * (η * u)) := by
    rw [← Real.exp_mul, ← Real.exp_nat_mul, ← Real.exp_add]; ring_nf
  unfold F
  rw [hg1, hg, r1, r2]
  rw [show A ^ 2 * Real.exp (-(ell u * η)) ^ 2 +
      B ^ 2 * (Real.exp (-(ell u * η)) ^ 2 * Real.exp (-2 * (η * u))) =
      Real.exp (-(ell u * η)) ^ 2 * (A ^ 2 + B ^ 2 * Real.exp (-2 * (η * u))) by ring,
    Real.sqrt_mul (by positivity), Real.sqrt_sq (by positivity)]

lemma F_swap (A B η x : ℝ) : F A B η (1 - x) = F B A η x := by
  unfold F; rw [sub_sub_cancel, add_comm]

lemma integrand_eq (A B η t : ℝ) : X01.integrand A B (1 / (1 + η)) t = F A B η t := by
  unfold X01.integrand F
  rw [show 2 / (1 / (1 + η)) - 2 = 2 * η by rw [one_div, div_inv_eq_mul]; ring]

lemma iiF (A B : ℝ) (hA : 0 < A) (hB : 0 < B) (η : ℝ) (hη : -1 < η) :
    IntervalIntegrable (F A B η) volume 0 1 := by
  have := X01.proof_1 A B hA hB (1 / (1 + η)) (one_div_pos.mpr (by linarith))
  rwa [show X01.integrand A B (1 / (1 + η)) = F A B η from funext (integrand_eq A B η)] at this

lemma pt (A B : ℝ) (hA : 0 < A) (η u : ℝ) :
    |gd u| • (F A B η (gg u) + F B A η (gg u)) =
      X01.r A B * (Real.exp (-u) * H A B (η * u) * K η u) := by
  have hr : X01.r A B ≠ 0 := (Real.sqrt_pos.mpr (by positivity)).ne'
  rw [smul_eq_mul, abs_of_pos (gd_pos u), F_gg, F_gg, K_eq]
  unfold H gd
  rw [← exp_ell]
  have e1 : Real.exp (-u) / Real.exp (ell u) ^ 2 =
      Real.exp (-u) * Real.exp (ell u * (-2 - η)) / Real.exp (-(ell u * η)) := by
    rw [← Real.exp_nat_mul, ← Real.exp_sub, ← Real.exp_add, ← Real.exp_sub]; ring_nf
  rw [e1]
  field_simp

end Y06

theorem Y06.proof_1 : Y06.stmt_1 := by
  intro A B hA hB η hη1 _hη2
  open Y06 MeasureTheory Set in
  have hr : 0 < X01.r A B := Real.sqrt_pos.mpr (by positivity)
  have h1η : 0 < 1 + η := by linarith
  set φ : ℝ → ℝ := fun t => F A B η t + F B A η t with hφ
  have iA := iiF A B hA hB η hη1
  have iB := iiF B A hB hA η hη1
  have hsub : ∀ {a b : ℝ}, 0 ≤ a → b ≤ 1 → a ≤ b → Set.uIcc a b ⊆ Set.uIcc (0:ℝ) 1 := by
    intro a b ha hb hab
    rw [Set.uIcc_of_le hab, Set.uIcc_of_le zero_le_one]; exact Set.Icc_subset_Icc ha hb
  have iA1 := iA.mono_set (hsub (by norm_num) le_rfl (by norm_num : (1:ℝ)/2 ≤ 1))
  have iA0 := iA.mono_set (hsub le_rfl (by norm_num) (by norm_num : (0:ℝ) ≤ 1/2))
  have iB1 := iB.mono_set (hsub (by norm_num) le_rfl (by norm_num : (1:ℝ)/2 ≤ 1))
  -- the integral over (0,1) as one over (1/2,1)
  have hsplit : ∫ t in (0:ℝ)..1, F A B η t = ∫ t in Set.Ioo (1/2 : ℝ) 1, φ t := by
    rw [← intervalIntegral.integral_add_adjacent_intervals iA0 iA1]
    have hrefl : ∫ t in (0:ℝ)..1/2, F A B η t = ∫ t in (1/2:ℝ)..1, F B A η t := by
      rw [← funext (F_swap A B η), intervalIntegral.integral_comp_sub_left (F A B η) 1]
      norm_num
    rw [hrefl, add_comm, ← intervalIntegral.integral_add iA1 iB1,
      intervalIntegral.integral_of_le (by norm_num), integral_Ioc_eq_integral_Ioo]
  have hφint : IntegrableOn φ (Set.Ioo (1/2 : ℝ) 1) := by
    have := (iA1.add iB1)
    rw [intervalIntegrable_iff_integrableOn_Ioc_of_le (by norm_num)] at this
    exact this.mono_set Set.Ioo_subset_Ioc_self
  have hd : ∀ x ∈ Set.Ioi (0:ℝ), HasDerivWithinAt gg (gd x) (Set.Ioi 0) x :=
    fun x _ => (hasDerivAt_gg x).hasDerivWithinAt
  have hinj : Set.InjOn gg (Set.Ioi 0) := gg_strictMono.injective.injOn
  have hcov := integral_image_eq_integral_abs_deriv_smul measurableSet_Ioi hd hinj φ
  have hint := (integrableOn_image_iff_integrableOn_abs_deriv_smul measurableSet_Ioi hd hinj φ).mp
    (by rw [gg_image]; exact hφint)
  rw [gg_image] at hcov
  have hpt : ∀ u ∈ Set.Ioi (0:ℝ), Real.exp (-u) * H A B (η * u) * K η u =
      (X01.r A B)⁻¹ * (|gd u| • φ (gg u)) := by
    intro u _
    rw [hφ]; simp only; rw [pt A B hA η u]; field_simp
  refine ⟨IntegrableOn.congr_fun (hint.const_mul (X01.r A B)⁻¹) (fun u hu => (hpt u hu).symm) measurableSet_Ioi, ?_⟩
  rw [setIntegral_congr_fun measurableSet_Ioi hpt, integral_const_mul, ← hcov, ← hsplit]
  unfold Ieta X01.P
  simp_rw [integrand_eq]
  field_simp

namespace Y06
open MeasureTheory Set Filter

lemma intOn_pow_exp (m : ℕ) {b : ℝ} (hb : 0 < b) :
    IntegrableOn (fun u : ℝ => u ^ m * Real.exp (-(b * u))) (Ioi 0) := by
  have h := integrableOn_rpow_mul_exp_neg_mul_rpow (s := (m : ℝ)) (p := 1)
    (by have := m.cast_nonneg (α := ℝ); linarith) one_pos hb
  refine h.congr_fun (fun u _ => ?_) measurableSet_Ioi
  simp only [Real.rpow_natCast, Real.rpow_one, neg_mul]

lemma int_pow_exp (m : ℕ) {b : ℝ} (hb : 0 < b) :
    ∫ u in Ioi (0:ℝ), u ^ m * Real.exp (-(b * u)) = (m.factorial : ℝ) / b ^ (m + 1) := by
  have h := Real.integral_rpow_mul_exp_neg_mul_Ioi (a := (m : ℝ) + 1) (by positivity) hb
  have e : ∫ u in Ioi (0:ℝ), u ^ m * Real.exp (-(b * u)) =
      ∫ u in Ioi (0:ℝ), u ^ ((m : ℝ) + 1 - 1) * Real.exp (-(b * u)) := by
    refine setIntegral_congr_fun measurableSet_Ioi (fun u _ => ?_)
    simp only [add_sub_cancel_right, Real.rpow_natCast]
  rw [e, h, Real.Gamma_nat_eq_factorial, show (m : ℝ) + 1 = ((m + 1 : ℕ) : ℝ) by push_cast; ring,
    Real.rpow_natCast, one_div, inv_pow]
  ring

lemma kj_meas (j : ℕ) : Measurable (kj j) := by
  unfold kj ell; fun_prop

lemma Mint_meas (n j : ℕ) :
    AEStronglyMeasurable (fun u : ℝ => Real.exp (-u) * u ^ (n - j) * kj j u)
      (volume.restrict (Ioi 0)) :=
  (by have := kj_meas j; fun_prop : Measurable _).aestronglyMeasurable

lemma Mint_integrable (n j : ℕ) :
    IntegrableOn (fun u : ℝ => Real.exp (-u) * u ^ (n - j) * kj j u) (Ioi 0) := by
  have hg : IntegrableOn (fun u : ℝ => u ^ (n - j) * Real.exp (-u)) (Ioi 0) := by
    simpa using intOn_pow_exp (n - j) one_pos
  refine Integrable.mono' hg (Mint_meas n j)
    (ae_restrict_of_forall_mem measurableSet_Ioi (fun u hu => ?_))
  have hu : 0 ≤ u := le_of_lt hu
  rw [Real.norm_eq_abs, abs_mul, abs_mul, abs_of_pos (Real.exp_pos _),
    abs_of_nonneg (pow_nonneg hu _)]
  calc Real.exp (-u) * u ^ (n - j) * |kj j u| ≤ Real.exp (-u) * u ^ (n - j) * 1 := by
        gcongr; exact kj_abs_le_one hu j
    _ = _ := by ring

lemma Mnj_bound (n j : ℕ) :
    |Mnj n j| ≤ ((n - j).factorial : ℝ) / ((j.factorial : ℝ) * ((j : ℝ) + 1) ^ (n - j + 1)) := by
  have hb : (0:ℝ) < j + 1 := by positivity
  have hg := (intOn_pow_exp (n - j) hb).div_const (j.factorial : ℝ)
  have h := norm_integral_le_of_norm_le hg
    (ae_restrict_of_forall_mem measurableSet_Ioi (fun u (hu : 0 < u) => by
      have hu' : 0 ≤ u := le_of_lt hu
      show ‖Real.exp (-u) * u ^ (n - j) * kj j u‖ ≤ u ^ (n - j) * Real.exp (-((j + 1) * u)) / j.factorial
      rw [Real.norm_eq_abs, abs_mul, abs_mul, abs_of_pos (Real.exp_pos _),
        abs_of_nonneg (pow_nonneg hu' _)]
      calc Real.exp (-u) * u ^ (n - j) * |kj j u|
          ≤ Real.exp (-u) * u ^ (n - j) * (Real.exp (-(j * u)) / j.factorial) := by
            gcongr; exact kj_abs_le hu' j
        _ = u ^ (n - j) * (Real.exp (-u) * Real.exp (-(j * u))) / j.factorial := by ring
        _ = _ := by rw [← Real.exp_add]; ring_nf))
  rw [integral_div, int_pow_exp (n - j) hb] at h
  unfold Mnj; rw [← Real.norm_eq_abs]
  refine h.trans (le_of_eq ?_)
  field_simp

lemma inv_sq_hasSum {x : ℝ} (hx0 : 0 ≤ x) (hx1 : x < 1) :
    HasSum (fun m : ℕ => (-1 : ℝ) ^ m * ((m : ℝ) + 1) * x ^ m) (((1 + x) ^ 2)⁻¹) := by
  have hr : ‖-x‖ < 1 := by rw [norm_neg, Real.norm_eq_abs, abs_of_nonneg hx0]; exact hx1
  have h3 := (hasSum_coe_mul_geometric_of_norm_lt_one hr).add (hasSum_geometric_of_norm_lt_one hr)
  have hk : ((1 + x) ^ 2)⁻¹ = -x / (1 - -x) ^ 2 + (1 - -x)⁻¹ := by
    have : (1 + x) ≠ 0 := by positivity
    rw [sub_neg_eq_add]; field_simp; ring
  rw [hk]
  convert h3 using 1
  funext m; rw [neg_pow]; ring

lemma kj0_hasSum {u : ℝ} (hu : 0 < u) :
    HasSum (fun m : ℕ => (-1 : ℝ) ^ m * ((m : ℝ) + 1) * Real.exp (-u) ^ m) (kj 0 u) := by
  have hx1 : Real.exp (-u) < 1 := by
    rw [← Real.exp_zero]; exact Real.exp_lt_exp.mpr (by linarith)
  unfold kj; simp only [pow_zero, mul_one, Nat.factorial_zero, Nat.cast_one, div_one]
  exact inv_sq_hasSum (Real.exp_pos _).le hx1

lemma Mn0_hasSum (n : ℕ) (hn : 2 ≤ n) :
    HasSum (fun m : ℕ => (-1 : ℝ) ^ m * (n.factorial : ℝ) / ((m : ℝ) + 1) ^ n) (Mnj n 0) := by
  set Fm : ℕ → ℝ → ℝ := fun m u =>
    (-1 : ℝ) ^ m * ((m : ℝ) + 1) * (u ^ n * Real.exp (-(((m : ℝ) + 1) * u))) with hFm
  have hint : ∀ m, Integrable (Fm m) (volume.restrict (Ioi 0)) := fun m =>
    (intOn_pow_exp n (by positivity : (0:ℝ) < m + 1)).const_mul _
  have hval : ∀ m, ∫ u in Ioi (0:ℝ), Fm m u =
      (-1 : ℝ) ^ m * (n.factorial : ℝ) / ((m : ℝ) + 1) ^ n := by
    intro m
    rw [hFm]; simp only
    rw [integral_const_mul, int_pow_exp n (by positivity)]
    have : (m : ℝ) + 1 ≠ 0 := by positivity
    field_simp; ring
  have hnorm : ∀ m, ∫ u in Ioi (0:ℝ), ‖Fm m u‖ = (n.factorial : ℝ) / ((m : ℝ) + 1) ^ n := by
    intro m
    have : ∀ u ∈ Ioi (0:ℝ), ‖Fm m u‖ = ((m : ℝ) + 1) * (u ^ n * Real.exp (-(((m : ℝ) + 1) * u))) := by
      intro u hu
      rw [hFm]; simp only
      rw [Real.norm_eq_abs, abs_mul, abs_mul, abs_pow, abs_neg, abs_one, one_pow, one_mul,
        abs_of_pos (by positivity : (0:ℝ) < m + 1),
        abs_of_nonneg (by have : (0:ℝ) ≤ u := le_of_lt hu; positivity)]
    rw [setIntegral_congr_fun measurableSet_Ioi this, integral_const_mul,
      int_pow_exp n (by positivity)]
    have : (m : ℝ) + 1 ≠ 0 := by positivity
    field_simp; ring
  have hsum : Summable (fun m => ∫ u in Ioi (0:ℝ), ‖Fm m u‖) := by
    simp_rw [hnorm]
    have := (Real.summable_one_div_nat_pow.mpr (show 1 < n by omega))
    have := (summable_nat_add_iff 1).mpr this
    refine (this.mul_left (n.factorial : ℝ)).congr (fun m => ?_)
    push_cast; ring
  have H := hasSum_integral_of_summable_integral_norm hint hsum
  simp_rw [hval] at H
  convert H using 1
  unfold Mnj
  refine setIntegral_congr_fun measurableSet_Ioi (fun u hu => ?_)
  have hs := (kj0_hasSum hu).mul_left (Real.exp (-u) * u ^ (n - 0))
  rw [← hs.tsum_eq]
  congr 1; funext m
  rw [hFm]; simp only [Nat.sub_zero]
  rw [← Real.exp_nat_mul, show Real.exp (-(((m : ℝ) + 1) * u)) = Real.exp (-u) * Real.exp (m * -u) by
    rw [← Real.exp_add]; ring_nf]
  ring

lemma etaD_hasSum (n : ℕ) (hn : 2 ≤ n) :
    HasSum (fun m : ℕ => (-1 : ℝ) ^ m / ((m : ℝ) + 1) ^ n) (Mnj n 0 / (n.factorial : ℝ)) := by
  have h := (Mn0_hasSum n hn).div_const (n.factorial : ℝ)
  convert h using 1
  funext m
  have : (n.factorial : ℝ) ≠ 0 := by positivity
  field_simp

end Y06

namespace Y06
open MeasureTheory Set Filter Topology

/-- Real parts of iterated complex derivatives of a holomorphic extension. -/
lemma bridge {f : ℂ → ℂ} {g : ℝ → ℝ} {U : Set ℂ} (hU : IsOpen U) (hf : DifferentiableOn ℂ f U)
    {V : Set ℝ} (hV : IsOpen V) (hVU : ∀ x ∈ V, (x : ℂ) ∈ U) (hfg : ∀ x ∈ V, g x = (f x).re) :
    ∀ k : ℕ, ∀ x ∈ V, iteratedDeriv k g x = (iteratedDeriv k f x).re := by
  intro k
  induction k with
  | zero => intro x hx; simpa using hfg x hx
  | succ k ih =>
    intro x hx
    have hev : iteratedDeriv k g =ᶠ[𝓝 x] fun y : ℝ => (iteratedDeriv k f y).re :=
      Filter.eventually_of_mem (hV.mem_nhds hx) ih
    have han : AnalyticAt ℂ (iteratedDeriv k f) (x : ℂ) := by
      rw [iteratedDeriv_eq_iterate]; exact (hf.analyticOnNhd hU).iterated_deriv k _ (hVU x hx)
    have hd : HasDerivAt (iteratedDeriv k f) (iteratedDeriv (k + 1) f x) (x : ℂ) := by
      rw [iteratedDeriv_succ]; exact han.differentiableAt.hasDerivAt
    rw [iteratedDeriv_succ, hev.deriv_eq, hd.real_of_complex.deriv]

/-- The square root branch `w ↦ w^(1/2)`. -/
noncomputable def csq (w : ℂ) : ℂ := w ^ ((2 : ℂ)⁻¹)

lemma csq_ofReal {x : ℝ} (hx : 0 ≤ x) : csq (x : ℂ) = (Real.sqrt x : ℂ) := by
  unfold csq
  rw [Real.sqrt_eq_rpow, Complex.ofReal_cpow hx]; push_cast; ring_nf

lemma csq_sq (w : ℂ) : csq w ^ 2 = w := Complex.cpow_ofNat_inv_pow w 2

/-- The cut of `1 + c e^{2z}` avoids the strip and the disc of radius `√((log c/2)^2+π^2/4)`. -/
lemma slit_aux {c : ℝ} (hc0 : 0 < c) (hc1 : c ≤ 1) (z : ℂ)
    (hz : ‖z‖ ^ 2 < (Real.log c / 2) ^ 2 + Real.pi ^ 2 / 4 ∨ |z.im| < Real.pi / 2) :
    1 + (c : ℂ) * Complex.exp (2 * z) ∈ Complex.slitPlane := by
  rw [Complex.mem_slitPlane_iff]
  by_contra hcon
  push Not at hcon
  obtain ⟨hre, him⟩ := hcon
  set x := z.re
  set y := z.im
  have e2re : (2 * z).re = 2 * x := by simp [x]
  have e2im : (2 * z).im = 2 * y := by simp [y]
  simp only [Complex.add_re, Complex.add_im, Complex.one_re, Complex.one_im, Complex.mul_re,
    Complex.mul_im, Complex.ofReal_re, Complex.ofReal_im, Complex.exp_re, Complex.exp_im,
    e2re, e2im, zero_mul, sub_zero, zero_add] at hre him
  have hE := Real.exp_pos (2 * x)
  have hs : Real.sin (2 * y) = 0 := by
    have : c * (Real.exp (2 * x) * Real.sin (2 * y)) = 0 := by linarith
    rcases mul_eq_zero.mp this with h | h
    · linarith
    · rcases mul_eq_zero.mp h with h | h
      · linarith
      · exact h
  have hcneg : Real.cos (2 * y) < 0 := by
    by_contra h; push Not at h
    have : 0 ≤ c * (Real.exp (2 * x) * Real.cos (2 * y)) := by positivity
    linarith
  have hcos : Real.cos (2 * y) = -1 := by
    have := Real.sin_sq_add_cos_sq (2 * y)
    rw [hs] at this
    nlinarith
  have hge : 1 ≤ c * Real.exp (2 * x) := by
    rw [hcos] at hre; linarith
  -- |y| ≥ π/2
  have hy : Real.pi / 2 ≤ |y| := by
    by_contra h; push Not at h
    have hcy : 0 < Real.cos y := Real.cos_pos_of_mem_Ioo ⟨by linarith [neg_abs_le y],
      by linarith [le_abs_self y]⟩
    have : Real.cos (2 * y) = 2 * Real.cos y ^ 2 - 1 := Real.cos_two_mul y
    nlinarith
  -- x ≥ -log c / 2 ≥ 0
  have hx : -(Real.log c / 2) ≤ x := by
    by_contra h; push Not at h
    have : c * Real.exp (2 * x) < 1 := by
      rw [← Real.exp_log hc0, ← Real.exp_add, Real.exp_lt_one_iff]  
      linarith
    linarith
  have hlc : Real.log c ≤ 0 := Real.log_nonpos hc0.le hc1
  have hnorm : ‖z‖ ^ 2 = x ^ 2 + y ^ 2 := by
    rw [Complex.sq_norm, Complex.normSq_apply]; ring
  rcases hz with hz | hz
  · rw [hnorm] at hz
    have h1 : (Real.log c / 2) ^ 2 ≤ x ^ 2 := by nlinarith
    have h2 : Real.pi ^ 2 / 4 ≤ y ^ 2 := by
      have := mul_le_mul hy hy (by positivity) (abs_nonneg y)
      rw [abs_mul_abs_self] at this; nlinarith
    linarith
  · linarith

end Y06

namespace Y06
open MeasureTheory Set Filter Topology

noncomputable def cc (A B : ℝ) : ℝ := (A / B) ^ 2

/-- Holomorphic extension of `H` (for `A ≤ B`). -/
noncomputable def HC (A B : ℝ) (z : ℂ) : ℂ :=
  ((B / X01.r A B : ℝ) : ℂ) * (Complex.exp (-z) * csq (1 + (cc A B : ℂ) * Complex.exp (2 * z)) +
    csq (1 + (cc A B : ℂ) * Complex.exp (2 * (-z))))

/-- Its domain: the disc of radius `ρ` together with the strip `|Im z| < π/2`. -/
def Udom (A B : ℝ) : Set ℂ := {z : ℂ | ‖z‖ < ρ A B} ∪ {z : ℂ | |z.im| < Real.pi / 2}

lemma Udom_open (A B : ℝ) : IsOpen (Udom A B) := by
  apply IsOpen.union
  · exact isOpen_lt continuous_norm continuous_const
  · exact isOpen_lt (by fun_prop) continuous_const

lemma ρ_pos (A B : ℝ) : 0 < ρ A B := Real.sqrt_pos.mpr (by have := Real.pi_pos; positivity)

lemma ρ_sq (A B : ℝ) : ρ A B ^ 2 = L A B ^ 2 + Real.pi ^ 2 / 4 :=
  Real.sq_sqrt (by positivity)

lemma cc_pos {A B : ℝ} (hA : 0 < A) (hB : 0 < B) : 0 < cc A B := by unfold cc; positivity

lemma cc_le {A B : ℝ} (hA : 0 < A) (hAB : A ≤ B) : cc A B ≤ 1 := by
  unfold cc; exact pow_le_one₀ (div_nonneg hA.le (by linarith)) ((div_le_one (by linarith)).mpr hAB)

lemma log_cc {A B : ℝ} (_hA : 0 < A) (_hB : 0 < B) : (Real.log (cc A B) / 2) ^ 2 = L A B ^ 2 := by
  unfold cc L
  rw [Real.log_pow, show A / B = (B / A)⁻¹ by rw [inv_div], Real.log_inv]; push_cast; ring

lemma slit_U {A B : ℝ} (hA : 0 < A) (hB : 0 < B) (hAB : A ≤ B) {z : ℂ} (hz : z ∈ Udom A B) :
    1 + (cc A B : ℂ) * Complex.exp (2 * z) ∈ Complex.slitPlane := by
  apply slit_aux (cc_pos hA hB) (cc_le hA hAB)
  rcases hz with hz | hz
  · left
    rw [log_cc hA hB, ← ρ_sq]
    have : ‖z‖ < ρ A B := hz
    exact pow_lt_pow_left₀ this (norm_nonneg _) two_ne_zero
  · right; exact hz

lemma neg_U {A B : ℝ} {z : ℂ} (hz : z ∈ Udom A B) : -z ∈ Udom A B := by
  rcases hz with hz | hz
  · left; show ‖-z‖ < _; rw [norm_neg]; exact hz
  · right; show |(-z).im| < _; rw [Complex.neg_im, abs_neg]; exact hz

lemma HC_diff {A B : ℝ} (hA : 0 < A) (hB : 0 < B) (hAB : A ≤ B) :
    DifferentiableOn ℂ (HC A B) (Udom A B) := by
  intro z hz
  have h1 : DifferentiableAt ℂ (fun z : ℂ => csq (1 + (cc A B : ℂ) * Complex.exp (2 * z))) z := by
    unfold csq
    exact DifferentiableAt.cpow_const (by fun_prop) (slit_U hA hB hAB hz)
  have h2 : DifferentiableAt ℂ (fun z : ℂ => csq (1 + (cc A B : ℂ) * Complex.exp (2 * (-z)))) z := by
    unfold csq
    exact DifferentiableAt.cpow_const (by fun_prop) (slit_U hA hB hAB (neg_U hz))
  unfold HC
  exact (((by fun_prop : DifferentiableAt ℂ (fun z : ℂ => Complex.exp (-z)) z).mul h1).add h2
    |>.const_mul _).differentiableWithinAt

lemma HC_real {A B : ℝ} (hA : 0 < A) (hB : 0 < B) (x : ℝ) : HC A B x = (H A B x : ℂ) := by
  have hc := cc_pos hA hB
  have e1 : (1 + (cc A B : ℂ) * Complex.exp (2 * (x : ℂ))) = ((1 + cc A B * Real.exp (2 * x) : ℝ) : ℂ) := by
    push_cast; rfl
  have e2 : (1 + (cc A B : ℂ) * Complex.exp (2 * (-(x : ℂ)))) =
      ((1 + cc A B * Real.exp (2 * (-x)) : ℝ) : ℂ) := by
    push_cast; rfl
  unfold HC
  rw [e1, e2, csq_ofReal (by positivity), csq_ofReal (by positivity),
    show Complex.exp (-(x : ℂ)) = ((Real.exp (-x) : ℝ) : ℂ) by push_cast; rfl]
  norm_cast
  -- real identity
  have hr : X01.r A B ≠ 0 := (Real.sqrt_pos.mpr (by positivity)).ne'
  have hE := Real.exp_pos (-x)
  have h1 : Real.exp (2 * x) = (Real.exp (-x) ^ 2)⁻¹ := by
    rw [← Real.exp_nat_mul, ← Real.exp_neg]; ring_nf
  have h2 : Real.exp (-2 * x) = Real.exp (-x) ^ 2 := by rw [← Real.exp_nat_mul]; ring_nf
  have h3 : Real.exp (2 * (-x)) = Real.exp (-x) ^ 2 := by rw [← Real.exp_nat_mul]; ring_nf
  have s1 : Real.sqrt (A ^ 2 + B ^ 2 * Real.exp (-2 * x)) =
      B * Real.exp (-x) * Real.sqrt (1 + cc A B * Real.exp (2 * x)) := by
    rw [show B * Real.exp (-x) * Real.sqrt (1 + cc A B * Real.exp (2 * x)) =
        Real.sqrt ((B * Real.exp (-x)) ^ 2 * (1 + cc A B * Real.exp (2 * x))) by
      rw [Real.sqrt_mul (sq_nonneg _), Real.sqrt_sq (by positivity)]]
    congr 1
    rw [h1, h2]; unfold cc; field_simp; ring
  have s2 : Real.sqrt (B ^ 2 + A ^ 2 * Real.exp (-2 * x)) =
      B * Real.sqrt (1 + cc A B * Real.exp (2 * (-x))) := by
    rw [show B * Real.sqrt (1 + cc A B * Real.exp (2 * (-x))) =
        Real.sqrt (B ^ 2 * (1 + cc A B * Real.exp (2 * (-x)))) by
      rw [Real.sqrt_mul (sq_nonneg _), Real.sqrt_sq hB.le]]
    congr 1
    rw [h3, h2]; unfold cc; field_simp
  unfold H
  rw [s1, s2]; field_simp

lemma r_swap (A B : ℝ) : X01.r B A = X01.r A B := by unfold X01.r; rw [add_comm]

lemma H_swap (A B : ℝ) : H B A = H A B := by
  funext x; unfold H; rw [r_swap, add_comm]

lemma ρ_swap (A B : ℝ) : ρ B A = ρ A B := by
  unfold ρ; rw [L_swap, neg_sq]

lemma Hc_swap (A B : ℝ) : Hc B A = Hc A B := by
  funext m; unfold Hc; rw [H_swap]

/-- Real iterated derivatives of `H` are real parts of those of `HC`. -/
lemma H_iter {A B : ℝ} (hA : 0 < A) (hB : 0 < B) (hAB : A ≤ B) (k : ℕ) (x : ℝ) :
    iteratedDeriv k (H A B) x = (iteratedDeriv k (HC A B) x).re :=
  bridge (Udom_open A B) (HC_diff hA hB hAB) isOpen_univ
    (fun x _ => by right; simp [Real.pi_pos]) (fun x _ => by rw [HC_real hA hB]; simp) k x
    (mem_univ x)

lemma stmt6_aux {A B : ℝ} (hA : 0 < A) (hB : 0 < B) (hAB : A ≤ B) {σ : ℝ} (hσ : 0 < σ)
    (hσρ : σ < ρ A B) : ∃ C : ℝ, ∀ m : ℕ, |Hc A B m| ≤ C * σ⁻¹ ^ m := by
  have hdiff := HC_diff hA hB hAB
  have hsub : Metric.closedBall (0:ℂ) σ ⊆ Udom A B := by
    intro z hz; left; show ‖z‖ < _
    have : ‖z‖ ≤ σ := by simpa using hz
    exact lt_of_le_of_lt this hσρ
  have hDC : DiffContOnCl ℂ (HC A B) (Metric.ball 0 σ) :=
    (hdiff.mono (by rw [closure_ball 0 hσ.ne']; exact hsub)).diffContOnCl
  obtain ⟨C, hC⟩ := (isCompact_sphere (0:ℂ) σ).exists_bound_of_continuousOn
    (hdiff.continuousOn.mono (Metric.sphere_subset_closedBall.trans hsub))
  refine ⟨C, fun m => ?_⟩
  have hb := Complex.norm_iteratedDeriv_le_of_forall_mem_sphere_norm_le m hσ hDC hC
  have hbr := H_iter hA hB hAB m 0
  rw [Complex.ofReal_zero] at hbr
  have hf : (0:ℝ) < m.factorial := by positivity
  unfold Hc
  rw [hbr, abs_div, Nat.abs_cast, div_le_iff₀ hf]
  calc |(iteratedDeriv m (HC A B) 0).re| ≤ ‖iteratedDeriv m (HC A B) 0‖ := Complex.abs_re_le_norm _
    _ ≤ m.factorial * C / σ ^ m := hb
    _ = C * σ⁻¹ ^ m * m.factorial := by rw [inv_pow]; field_simp

end Y06

theorem Y06.proof_6 : Y06.stmt_6 := by
  intro A B hA hB σ hσ hσρ
  rcases le_total A B with hAB | hBA
  · exact Y06.stmt6_aux hA hB hAB hσ hσρ
  · rw [← Y06.Hc_swap]
    exact Y06.stmt6_aux hB hA hBA hσ (by rwa [Y06.ρ_swap])

namespace Y06
open MeasureTheory Set Filter Topology

lemma real_mem_U (A B x : ℝ) : (x : ℂ) ∈ Udom A B := by right; simp [Real.pi_pos]

lemma H_hasDeriv_aux {A B : ℝ} (hA : 0 < A) (hB : 0 < B) (hAB : A ≤ B) (m : ℕ) (x : ℝ) :
    HasDerivAt (iteratedDeriv m (H A B)) (iteratedDeriv (m + 1) (H A B) x) x := by
  have h := funext (H_iter hA hB hAB m)
  rw [h, H_iter hA hB hAB (m + 1)]
  have han : AnalyticAt ℂ (iteratedDeriv m (HC A B)) (x : ℂ) := by
    rw [iteratedDeriv_eq_iterate]
    exact ((HC_diff hA hB hAB).analyticOnNhd (Udom_open A B)).iterated_deriv m _ (real_mem_U A B x)
  have hd : HasDerivAt (iteratedDeriv m (HC A B)) (iteratedDeriv (m + 1) (HC A B) x) (x : ℂ) := by
    rw [iteratedDeriv_succ]; exact han.differentiableAt.hasDerivAt
  exact hd.real_of_complex

lemma H_hasDeriv {A B : ℝ} (hA : 0 < A) (hB : 0 < B) (m : ℕ) (x : ℝ) :
    HasDerivAt (iteratedDeriv m (H A B)) (iteratedDeriv (m + 1) (H A B) x) x := by
  rcases le_total A B with hAB | hBA
  · exact H_hasDeriv_aux hA hB hAB m x
  · rw [← H_swap A B]; exact H_hasDeriv_aux hB hA hBA m x

lemma H_contDiff {A B : ℝ} (hA : 0 < A) (hB : 0 < B) : ContDiff ℝ (⊤ : ℕ∞) (H A B) :=
  contDiff_of_differentiable_iteratedDeriv
    (fun m _ => fun x => (H_hasDeriv hA hB m x).differentiableAt)

lemma H_iter_cont {A B : ℝ} (hA : 0 < A) (hB : 0 < B) (m : ℕ) :
    Continuous (iteratedDeriv m (H A B)) :=
  continuous_iff_continuousAt.mpr (fun x => (H_hasDeriv hA hB m x).continuousAt)

lemma csq_norm_le {c : ℝ} (hc0 : 0 < c) (hc1 : c ≤ 1) (w : ℂ) :
    ‖csq (1 + (c : ℂ) * Complex.exp (2 * w))‖ ≤ 1 + Real.exp w.re := by
  have h2 : ‖csq (1 + (c : ℂ) * Complex.exp (2 * w))‖ ^ 2 = ‖1 + (c : ℂ) * Complex.exp (2 * w)‖ := by
    rw [← norm_pow, csq_sq]
  have h3 : ‖1 + (c : ℂ) * Complex.exp (2 * w)‖ ≤ 1 + c * Real.exp (2 * w.re) := by
    refine (norm_add_le _ _).trans ?_
    rw [norm_one, norm_mul, Complex.norm_exp, Complex.norm_real, Real.norm_of_nonneg hc0.le]
    simp
  have he := Real.exp_pos w.re
  have h4 : Real.exp (2 * w.re) = Real.exp w.re ^ 2 := by rw [← Real.exp_nat_mul]; ring_nf
  have h5 : c * Real.exp (2 * w.re) ≤ Real.exp w.re ^ 2 := by
    rw [h4]; nlinarith [sq_nonneg (Real.exp w.re)]
  nlinarith [norm_nonneg (csq (1 + (c : ℂ) * Complex.exp (2 * w)))]

lemma HC_norm_le {A B : ℝ} (hA : 0 < A) (hB : 0 < B) (hAB : A ≤ B) (z : ℂ) :
    ‖HC A B z‖ ≤ B / X01.r A B * (2 + 2 * Real.exp (-z.re)) := by
  have hr : 0 < X01.r A B := Real.sqrt_pos.mpr (by positivity)
  have hc0 := cc_pos hA hB
  have hc1 := cc_le hA hAB
  have h1 := csq_norm_le hc0 hc1 z
  have h2 := csq_norm_le hc0 hc1 (-z)
  rw [Complex.neg_re] at h2
  unfold HC
  rw [norm_mul, Complex.norm_real, Real.norm_of_nonneg (by positivity)]
  gcongr
  refine (norm_add_le _ _).trans ?_
  rw [norm_mul, Complex.norm_exp, Complex.neg_re]
  have hm : Real.exp (-z.re) * Real.exp z.re = 1 := by rw [← Real.exp_add]; simp
  have := mul_le_mul_of_nonneg_left h1 (Real.exp_pos (-z.re)).le
  nlinarith

lemma H_bound_aux {A B : ℝ} (hA : 0 < A) (hB : 0 < B) (hAB : A ≤ B) (m : ℕ) (x : ℝ) :
    |iteratedDeriv m (H A B) x| ≤
      m.factorial * (2 * Real.exp 1 * (B / X01.r A B) * (1 + Real.exp (-x))) * (4 / Real.pi) ^ m := by
  have hr : 0 < X01.r A B := Real.sqrt_pos.mpr (by positivity)
  have hpi := Real.pi_pos
  have hpi4 := Real.pi_lt_four
  have hR : 0 < Real.pi / 4 := by positivity
  have hdiff := HC_diff hA hB hAB
  have hsub : Metric.closedBall (x : ℂ) (Real.pi / 4) ⊆ Udom A B := by
    intro z hz; refine Or.inr ?_
    have h1 : ‖z - x‖ ≤ Real.pi / 4 := by simpa [dist_eq_norm] using hz
    have h2 := Complex.abs_im_le_norm (z - x)
    simp only [Complex.sub_im, Complex.ofReal_im, sub_zero] at h2
    exact (h2.trans h1).trans_lt (show Real.pi / 4 < Real.pi / 2 by linarith)
  have hDC : DiffContOnCl ℂ (HC A B) (Metric.ball (x : ℂ) (Real.pi / 4)) :=
    (hdiff.mono (by rw [closure_ball _ hR.ne']; exact hsub)).diffContOnCl
  have hC : ∀ z ∈ Metric.sphere (x : ℂ) (Real.pi / 4),
      ‖HC A B z‖ ≤ 2 * Real.exp 1 * (B / X01.r A B) * (1 + Real.exp (-x)) := by
    intro z hz
    have h1 : ‖z - x‖ = Real.pi / 4 := by simpa [dist_eq_norm] using hz
    have h2 := Complex.abs_re_le_norm (z - x)
    simp only [Complex.sub_re, Complex.ofReal_re] at h2
    rw [h1] at h2
    have h3 : Real.exp (-z.re) ≤ Real.exp 1 * Real.exp (-x) := by
      rw [← Real.exp_add]; apply Real.exp_le_exp.mpr
      have := (abs_le.mp h2).1; linarith
    have he1 : 1 ≤ Real.exp 1 := Real.one_le_exp (by norm_num)
    have hBr : 0 ≤ B / X01.r A B := by positivity
    refine (HC_norm_le hA hB hAB z).trans ?_
    have hex := Real.exp_pos (-x)
    have : 2 + 2 * Real.exp (-z.re) ≤ 2 * Real.exp 1 * (1 + Real.exp (-x)) := by nlinarith
    calc B / X01.r A B * (2 + 2 * Real.exp (-z.re)) ≤
        B / X01.r A B * (2 * Real.exp 1 * (1 + Real.exp (-x))) := by gcongr
      _ = _ := by ring
  have hb := Complex.norm_iteratedDeriv_le_of_forall_mem_sphere_norm_le m hR hDC hC
  rw [H_iter hA hB hAB]
  refine (Complex.abs_re_le_norm _).trans (hb.trans (le_of_eq ?_))
  rw [div_pow, div_pow]; field_simp

lemma H_bound {A B : ℝ} (hA : 0 < A) (hB : 0 < B) :
    ∃ C : ℝ, 0 ≤ C ∧ ∀ (m : ℕ) (x : ℝ),
      |iteratedDeriv m (H A B) x| ≤ m.factorial * (C * (1 + Real.exp (-x))) * (4 / Real.pi) ^ m := by
  rcases le_total A B with hAB | hBA
  · exact ⟨_, by have : 0 < X01.r A B := Real.sqrt_pos.mpr (by positivity)
                 positivity, H_bound_aux hA hB hAB⟩
  · refine ⟨2 * Real.exp 1 * (A / X01.r B A), by
      have : 0 < X01.r B A := Real.sqrt_pos.mpr (by positivity)
      positivity, fun m x => ?_⟩
    rw [← H_swap A B]; exact H_bound_aux hB hA hBA m x

end Y06

namespace Y06
open MeasureTheory Set Filter Topology

/-- The `n`-th `η`-derivative of `H(ηu)K(η,u)`, by the Leibniz rule. -/
noncomputable def Dn (A B : ℝ) (n : ℕ) (η u : ℝ) : ℝ :=
  ∑ i ∈ Finset.range (n + 1), (n.choose i : ℝ) * (u ^ i * iteratedDeriv i (H A B) (u * η)) *
    ((-ell u) ^ (n - i) * K η u)

lemma K_contDiff (u : ℝ) : ContDiff ℝ (⊤ : ℕ∞) (fun η => K η u) := by
  have : (fun η => K η u) = fun η => Real.exp (ell u * (-2 - η)) := funext (fun η => K_eq η u)
  rw [this]; fun_prop

lemma ψ_contDiff {A B : ℝ} (hA : 0 < A) (hB : 0 < B) (u : ℝ) :
    ContDiff ℝ (⊤ : ℕ∞) (fun η => H A B (η * u) * K η u) :=
  ((H_contDiff hA hB).comp (contDiff_id.mul contDiff_const)).mul (K_contDiff u)

lemma ψ_iter {A B : ℝ} (hA : 0 < A) (hB : 0 < B) (u : ℝ) (n : ℕ) (η : ℝ) :
    iteratedDeriv n (fun η => H A B (η * u) * K η u) η = Dn A B n η u := by
  have hH : ContDiff ℝ (⊤ : ℕ∞) (fun η => H A B (η * u)) :=
    (H_contDiff hA hB).comp (contDiff_id.mul contDiff_const)
  rw [iteratedDeriv_fun_mul (hH.contDiffAt.of_le (mod_cast le_top))
    ((K_contDiff u).contDiffAt.of_le (mod_cast le_top))]
  unfold Dn
  refine Finset.sum_congr rfl (fun i _ => ?_)
  have e : (fun η => H A B (η * u)) = fun η => H A B (u * η) := by funext η; rw [mul_comm]
  rw [e, iteratedDeriv_comp_const_mul ((H_contDiff hA hB).of_le (mod_cast le_top)) u, iter_K]

lemma Dn_hasDeriv {A B : ℝ} (hA : 0 < A) (hB : 0 < B) (n : ℕ) (u η : ℝ) :
    HasDerivAt (fun η => Dn A B n η u) (Dn A B (n + 1) η u) η := by
  have e : (fun η => Dn A B n η u) = iteratedDeriv n (fun η => H A B (η * u) * K η u) :=
    funext (fun η => (ψ_iter hA hB u n η).symm)
  rw [e, ← ψ_iter hA hB u (n + 1) η, iteratedDeriv_succ]
  exact (((ψ_contDiff hA hB u).of_le (mod_cast le_top : ((n + 1 : ℕ) : WithTop ℕ∞) ≤ _)).differentiable_iteratedDeriv' n
    η).hasDerivAt

lemma Dn_meas {A B : ℝ} (hA : 0 < A) (hB : 0 < B) (n : ℕ) (η : ℝ) :
    Measurable (fun u => Real.exp (-u) * Dn A B n η u) := by
  refine Measurable.mul (by fun_prop) ?_
  show Measurable (fun u => ∑ i ∈ Finset.range (n + 1), (n.choose i : ℝ) *
    (u ^ i * iteratedDeriv i (H A B) (u * η)) * ((-ell u) ^ (n - i) * K η u))
  refine Finset.measurable_sum _ (fun i _ => ?_)
  have := (H_iter_cont hA hB i).measurable
  unfold K ell
  fun_prop

lemma Kell_le {u : ℝ} (hu : 0 ≤ u) {η : ℝ} (hη : -2 ≤ η) (j : ℕ) : |(-ell u) ^ j * K η u| ≤ 1 := by
  rw [abs_mul, abs_pow, abs_neg, abs_of_nonneg (ell_nonneg u), abs_of_pos (K_pos η u)]
  calc ell u ^ j * K η u ≤ 1 * 1 :=
        mul_le_mul (pow_le_one₀ (ell_nonneg u) (ell_le_one hu)) (K_le_one hη) (K_pos η u).le zero_le_one
    _ = 1 := by ring

lemma pow_le_one_add {u : ℝ} (hu : 0 ≤ u) {i n : ℕ} (hi : i ≤ n) : u ^ i ≤ 1 + u ^ n := by
  rcases le_total u 1 with h | h
  · have := pow_le_one₀ hu h (n := i); have := pow_nonneg hu n; linarith
  · have := pow_le_pow_right₀ h hi; linarith

/-- The constant in the bound for `Dn`. -/
noncomputable def Cn (C1 : ℝ) (n : ℕ) : ℝ :=
  ∑ i ∈ Finset.range (n + 1), (n.choose i : ℝ) * ((i.factorial : ℝ) * C1 * (4 / Real.pi) ^ i)

lemma Dn_bound {A B C1 : ℝ} (hC1 : 0 ≤ C1)
    (hHb : ∀ (m : ℕ) (x : ℝ),
      |iteratedDeriv m (H A B) x| ≤ m.factorial * (C1 * (1 + Real.exp (-x))) * (4 / Real.pi) ^ m)
    (n : ℕ) {δ η u : ℝ} (hu : 0 ≤ u) (hη1 : -δ ≤ η) (hη2 : -2 ≤ η) :
    |Dn A B n η u| ≤ Cn C1 n * ((1 + u ^ n) * (1 + Real.exp (δ * u))) := by
  unfold Dn Cn
  rw [Finset.sum_mul]
  refine (Finset.abs_sum_le_sum_abs _ _).trans (Finset.sum_le_sum (fun i hi => ?_))
  have hin : i ≤ n := Nat.lt_succ_iff.mp (Finset.mem_range.mp hi)
  have hpi := Real.pi_pos
  have h1 := hHb i (u * η)
  have h2 := Kell_le hu hη2 (n - i)
  have h3 : 1 + Real.exp (-(u * η)) ≤ 1 + Real.exp (δ * u) := by
    have : -(u * η) ≤ δ * u := by nlinarith
    linarith [Real.exp_le_exp.mpr this]
  have h4 := pow_le_one_add hu hin
  rw [abs_mul, abs_mul, abs_mul, Nat.abs_cast, abs_of_nonneg (pow_nonneg hu i)]
  have hc : (0:ℝ) ≤ n.choose i := Nat.cast_nonneg _
  have hf : (0:ℝ) ≤ i.factorial := Nat.cast_nonneg _
  have hq : (0:ℝ) ≤ (4 / Real.pi) ^ i := by positivity
  have h5 : |iteratedDeriv i (H A B) (u * η)| ≤
      i.factorial * (C1 * (1 + Real.exp (δ * u))) * (4 / Real.pi) ^ i :=
    h1.trans (mul_le_mul_of_nonneg_right (mul_le_mul_of_nonneg_left
      (mul_le_mul_of_nonneg_left h3 hC1) hf) hq)
  calc (n.choose i : ℝ) * (u ^ i * |iteratedDeriv i (H A B) (u * η)|) * |(-ell u) ^ (n - i) * K η u|
      ≤ (n.choose i : ℝ) * (u ^ i * |iteratedDeriv i (H A B) (u * η)|) * 1 :=
        mul_le_mul_of_nonneg_left h2 (by positivity)
    _ ≤ (n.choose i : ℝ) * ((1 + u ^ n) *
          (i.factorial * (C1 * (1 + Real.exp (δ * u))) * (4 / Real.pi) ^ i)) * 1 := by
        refine mul_le_mul_of_nonneg_right (mul_le_mul_of_nonneg_left
          (mul_le_mul h4 h5 (abs_nonneg _) (by positivity)) hc) zero_le_one
    _ = _ := by ring

/-- The integrable dominating function. -/
lemma dom_integrable (n : ℕ) {δ : ℝ} (hδ : δ < 1) (c : ℝ) :
    IntegrableOn (fun u : ℝ => c * ((1 + u ^ n) * (Real.exp (-(1 * u)) + Real.exp (-((1 - δ) * u)))))
      (Ioi 0) := by
  have h1 := intOn_pow_exp 0 (one_pos)
  have h2 := intOn_pow_exp n (one_pos)
  have h3 := intOn_pow_exp 0 (by linarith : (0:ℝ) < 1 - δ)
  have h4 := intOn_pow_exp n (by linarith : (0:ℝ) < 1 - δ)
  refine IntegrableOn.congr_fun (((h1.add h2).add (h3.add h4)).const_mul c) (fun u _ => ?_)
    measurableSet_Ioi
  simp only [Pi.add_apply, pow_zero]; ring

lemma exp_dom (n : ℕ) {C δ u : ℝ} {D : ℝ}
    (hD : |D| ≤ C * ((1 + u ^ n) * (1 + Real.exp (δ * u)))) :
    ‖Real.exp (-u) * D‖ ≤ C * ((1 + u ^ n) * (Real.exp (-(1 * u)) + Real.exp (-((1 - δ) * u)))) := by
  rw [Real.norm_eq_abs, abs_mul, abs_of_pos (Real.exp_pos _)]
  have e : Real.exp (-u) * (1 + Real.exp (δ * u)) = Real.exp (-(1 * u)) + Real.exp (-((1 - δ) * u)) := by
    rw [mul_add, ← Real.exp_add]; ring_nf
  calc Real.exp (-u) * |D| ≤ Real.exp (-u) * (C * ((1 + u ^ n) * (1 + Real.exp (δ * u)))) := by
        gcongr
    _ = _ := by rw [← e]; ring

noncomputable def FF (A B : ℝ) (n : ℕ) (η : ℝ) : ℝ :=
  ∫ u in Ioi (0:ℝ), Real.exp (-u) * Dn A B n η u

lemma FF_hasDeriv {A B : ℝ} (hA : 0 < A) (hB : 0 < B) (n : ℕ) {η₀ : ℝ} (h : η₀ ∈ Ioo (-1 : ℝ) 1) :
    HasDerivAt (FF A B n) (FF A B (n + 1) η₀) η₀ := by
  obtain ⟨C1, hC1, hHb⟩ := H_bound hA hB
  have hab : |η₀| < 1 := abs_lt.mpr h
  set δ := (|η₀| + 1) / 2 with hδ
  have hδ1 : δ < 1 := by rw [hδ]; linarith
  have hδ0 : |η₀| < δ := by rw [hδ]; linarith
  have hs : Ioo (-δ) δ ∈ 𝓝 η₀ := Ioo_mem_nhds (by linarith [neg_abs_le η₀]) (by linarith [le_abs_self η₀])
  have hCn : ∀ m, 0 ≤ Cn C1 m := fun m => by unfold Cn; positivity
  have hbd : ∀ m, ∀ η ∈ Ioo (-δ) δ, ∀ u ∈ Ioi (0:ℝ),
      ‖Real.exp (-u) * Dn A B m η u‖ ≤
        Cn C1 m * ((1 + u ^ m) * (Real.exp (-(1 * u)) + Real.exp (-((1 - δ) * u)))) := by
    intro m η hη u hu
    exact exp_dom m
      (Dn_bound hC1 hHb m (le_of_lt hu) hη.1.le (by linarith [hη.1]))
  have hint0 : Integrable (fun u => Real.exp (-u) * Dn A B n η₀ u) (volume.restrict (Ioi 0)) :=
    Integrable.mono' (dom_integrable n hδ1 (Cn C1 n)) (Dn_meas hA hB n η₀).aestronglyMeasurable
      (ae_restrict_of_forall_mem measurableSet_Ioi
        (hbd n η₀ ⟨by linarith [neg_abs_le η₀], by linarith [le_abs_self η₀]⟩))
  have key := hasDerivAt_integral_of_dominated_loc_of_deriv_le (μ := volume.restrict (Ioi 0))
    (F := fun η u => Real.exp (-u) * Dn A B n η u)
    (F' := fun η u => Real.exp (-u) * Dn A B (n + 1) η u) hs
    (Filter.Eventually.of_forall (fun η => (Dn_meas hA hB n η).aestronglyMeasurable)) hint0
    (Dn_meas hA hB (n + 1) η₀).aestronglyMeasurable
    (ae_restrict_of_forall_mem measurableSet_Ioi (fun u hu η hη => hbd (n + 1) η hη u hu))
    (dom_integrable (n + 1) hδ1 (Cn C1 (n + 1)))
    (Filter.Eventually.of_forall (fun u η _ => (Dn_hasDeriv hA hB n u η).const_mul _))
  exact key.2

lemma FF_contDiffOn {A B : ℝ} (hA : 0 < A) (hB : 0 < B) :
    ∀ k : ℕ, ∀ n : ℕ, ContDiffOn ℝ k (FF A B n) (Ioo (-1) 1) := by
  intro k
  induction k with
  | zero =>
    intro n
    exact contDiffOn_zero.mpr (fun η h => (FF_hasDeriv hA hB n h).continuousAt.continuousWithinAt)
  | succ k ih =>
    intro n
    push_cast
    rw [contDiffOn_succ_iff_deriv_of_isOpen isOpen_Ioo]
    refine ⟨fun η h => (FF_hasDeriv hA hB n h).differentiableAt.differentiableWithinAt,
      fun h => absurd h (by simp), (ih (n + 1)).congr (fun η h => (FF_hasDeriv hA hB n h).deriv)⟩

lemma Ieta_eq_FF {A B : ℝ} (hA : 0 < A) (hB : 0 < B) {η : ℝ} (h : η ∈ Ioo (-1 : ℝ) 1) :
    Ieta A B η = FF A B 0 η := by
  rw [(Y06.proof_1 A B hA hB η h.1 h.2).2]
  unfold FF Dn
  congr 1; funext u
  simp only [zero_add, Finset.sum_range_one, Nat.choose_self, Nat.cast_one, one_mul, pow_zero,
    iteratedDeriv_zero, Nat.sub_self, mul_comm u η]
  ring

end Y06

theorem Y06.proof_3 : Y06.stmt_3 := by
  intro A B hA hB
  have : ContDiffOn ℝ (⊤ : ℕ∞) (Y06.FF A B 0) (Set.Ioo (-1) 1) :=
    contDiffOn_infty.mpr (fun k => Y06.FF_contDiffOn hA hB k 0)
  exact this.congr (fun η h => Y06.Ieta_eq_FF hA hB h)

namespace Y06
open MeasureTheory Set Filter Topology

lemma iter_FF {A B : ℝ} (hA : 0 < A) (hB : 0 < B) :
    ∀ n : ℕ, ∀ η ∈ Ioo (-1 : ℝ) 1, iteratedDeriv n (FF A B 0) η = FF A B n η := by
  intro n
  induction n with
  | zero => intro η _; simp
  | succ n ih =>
    intro η h
    rw [iteratedDeriv_succ]
    have hev : iteratedDeriv n (FF A B 0) =ᶠ[𝓝 η] FF A B n :=
      Filter.eventually_of_mem (isOpen_Ioo.mem_nhds h) ih
    rw [hev.deriv_eq, (FF_hasDeriv hA hB n h).deriv]

lemma In_eq_FF {A B : ℝ} (hA : 0 < A) (hB : 0 < B) (n : ℕ) :
    In A B n = FF A B n 0 / (n.factorial : ℝ) := by
  have h0 : (0:ℝ) ∈ Ioo (-1 : ℝ) 1 := ⟨by norm_num, by norm_num⟩
  have hev : Ieta A B =ᶠ[𝓝 0] FF A B 0 :=
    Filter.eventually_of_mem (isOpen_Ioo.mem_nhds h0) (fun η h => Ieta_eq_FF hA hB h)
  unfold In
  rw [hev.iteratedDeriv_eq n, iter_FF hA hB n 0 h0]

lemma K_zero (u : ℝ) : K 0 u = ((1 + Real.exp (-u)) ^ 2)⁻¹ := by
  unfold K
  rw [sub_zero, Real.rpow_neg (by positivity), show (2:ℝ) = ((2:ℕ):ℝ) by norm_num,
    Real.rpow_natCast]

lemma Dn_zero_eq {A B : ℝ} (n : ℕ) (u : ℝ) :
    Real.exp (-u) * Dn A B n 0 u = ∑ j ∈ Finset.range (n + 1),
      (n.factorial : ℝ) * Hc A B (n - j) * (Real.exp (-u) * u ^ (n - j) * kj j u) := by
  unfold Dn
  rw [Finset.mul_sum, ← Finset.sum_range_reflect]
  refine Finset.sum_congr rfl (fun j hj => ?_)
  have hjn : j ≤ n := Nat.lt_succ_iff.mp (Finset.mem_range.mp hj)
  rw [show n + 1 - 1 - j = n - j by omega, show n - (n - j) = j by omega, Nat.choose_symm hjn,
    mul_zero, K_zero]
  unfold Hc kj
  have hfac : (n.factorial : ℝ) = n.choose j * j.factorial * (n - j).factorial := by
    exact_mod_cast (Nat.choose_mul_factorial_mul_factorial hjn).symm
  rw [hfac]
  have h1 : (j.factorial : ℝ) ≠ 0 := by positivity
  have h2 : ((n - j).factorial : ℝ) ≠ 0 := by positivity
  field_simp

lemma In_sum {A B : ℝ} (hA : 0 < A) (hB : 0 < B) (n : ℕ) :
    In A B n = ∑ j ∈ Finset.range (n + 1), Hc A B (n - j) * Mnj n j := by
  rw [In_eq_FF hA hB, FF]
  simp_rw [Dn_zero_eq]
  rw [integral_finsetSum _ (fun j _ => (Mint_integrable n j).const_mul _), Finset.sum_div]
  refine Finset.sum_congr rfl (fun j _ => ?_)
  rw [integral_const_mul, Mnj]
  have : (n.factorial : ℝ) ≠ 0 := by positivity
  field_simp

end Y06

theorem Y06.proof_4 : Y06.stmt_4 := by
  refine ⟨fun n j _ => Y06.Mint_integrable n j, fun A B hA hB n => Y06.In_sum hA hB n,
    fun n hn => Y06.etaD_hasSum n hn, fun n j _ _ => Y06.Mnj_bound n j⟩

namespace Y06
open MeasureTheory Set Filter Topology

noncomputable def S2 : ℝ := ∑' m : ℕ, 1 / ((m : ℝ) + 1) ^ 2

lemma S2_summable : Summable (fun m : ℕ => 1 / ((m : ℝ) + 1) ^ 2) := by
  have := (summable_nat_add_iff 1).mpr (Real.summable_one_div_nat_pow.mpr one_lt_two)
  refine this.congr (fun m => ?_); push_cast; ring

lemma pow_bound (m n : ℕ) (hn : 2 ≤ n) :
    1 / ((m : ℝ) + 2) ^ n ≤ 4 * (2⁻¹ : ℝ) ^ n * (1 / ((m : ℝ) + 1) ^ 2) := by
  obtain ⟨k, rfl⟩ : ∃ k, n = k + 2 := ⟨n - 2, by omega⟩
  have hm : (0:ℝ) ≤ m := m.cast_nonneg
  have h1 : (2:ℝ) ^ k ≤ ((m : ℝ) + 2) ^ k := pow_le_pow_left₀ (by norm_num) (by linarith) k
  have h2 : ((m : ℝ) + 1) ^ 2 ≤ ((m : ℝ) + 2) ^ 2 := pow_le_pow_left₀ (by linarith) (by linarith) 2
  rw [inv_pow, div_le_iff₀ (by positivity)]
  have : (2:ℝ) ^ (k + 2) * ((m : ℝ) + 1) ^ 2 ≤ 4 * ((m : ℝ) + 2) ^ (k + 2) := by
    rw [pow_add, pow_add]
    have := mul_le_mul h1 h2 (by positivity) (by positivity)
    nlinarith
  field_simp
  linarith

lemma etaD_sub_bound (n : ℕ) (hn : 2 ≤ n) :
    |Mnj n 0 / (n.factorial : ℝ) - 1| ≤ 4 * (2⁻¹ : ℝ) ^ n * S2 := by
  have h := etaD_hasSum n hn
  have h1 : HasSum (fun m : ℕ => (-1 : ℝ) ^ (m + 1) / (((m + 1 : ℕ) : ℝ) + 1) ^ n)
      (Mnj n 0 / (n.factorial : ℝ) - 1) := by
    have := (hasSum_nat_add_iff' 1).mpr h
    simpa using this
  have hS := (S2_summable.hasSum).mul_left (4 * (2⁻¹ : ℝ) ^ n)
  have hb : ∀ m : ℕ, |(-1 : ℝ) ^ (m + 1) / (((m + 1 : ℕ) : ℝ) + 1) ^ n| ≤
      4 * (2⁻¹ : ℝ) ^ n * (1 / ((m : ℝ) + 1) ^ 2) := by
    intro m
    rw [abs_div, abs_pow, abs_neg, abs_one, one_pow, abs_of_pos (by positivity)]
    have := pow_bound m n hn
    push_cast
    calc 1 / ((m : ℝ) + 1 + 1) ^ n = 1 / ((m : ℝ) + 2) ^ n := by ring_nf
      _ ≤ _ := this
  rw [abs_le]
  constructor
  · have := hasSum_le (fun m => (abs_le.mp (hb m)).1) hS.neg h1
    unfold S2 at *; linarith
  · have := hasSum_le (fun m => (abs_le.mp (hb m)).2) h1 hS
    unfold S2 at *; linarith

lemma Mnj_term_bound {A B : ℝ} {σ C : ℝ} (hσ : 0 < σ) (hC : ∀ m : ℕ, |Hc A B m| ≤ C * σ⁻¹ ^ m)
    (n j : ℕ) (hj1 : 1 ≤ j) (hjn : j ≤ n) :
    |Hc A B (n - j) * Mnj n j| / (n.factorial : ℝ) ≤ C * (2 * σ)⁻¹ ^ n * ((2 * σ) ^ j / j.factorial) := by
  have hC0 : 0 ≤ C := by have := hC 0; simp at this; linarith [abs_nonneg (Hc A B 0)]
  have hM := Mnj_bound n j
  have hH := hC (n - j)
  rw [abs_mul]
  have hfac : (j.factorial : ℝ) * ((n - j).factorial : ℝ) ≤ n.factorial := by
    have := Nat.choose_mul_factorial_mul_factorial hjn
    have h2 : j.factorial * (n - j).factorial ≤ n.factorial := by
      rw [← this, mul_assoc]; exact Nat.le_mul_of_pos_left _ (Nat.choose_pos hjn)
    exact_mod_cast h2
  have hpow : ((j : ℝ) + 1) ^ (n - j + 1) ≥ 2 ^ (n - j) := by
    calc ((j : ℝ) + 1) ^ (n - j + 1) ≥ ((j : ℝ) + 1) ^ (n - j) :=
          pow_le_pow_right₀ (by have := j.cast_nonneg (α := ℝ); linarith) (by omega)
      _ ≥ 2 ^ (n - j) := pow_le_pow_left₀ (by norm_num) (by
          have : (1:ℝ) ≤ j := by exact_mod_cast hj1
          linarith) _
  have hjf : (1:ℝ) ≤ j.factorial := by exact_mod_cast Nat.one_le_iff_ne_zero.mpr (Nat.factorial_ne_zero j)
  have hnf : (0:ℝ) < n.factorial := by positivity
  have hjf0 : (0:ℝ) < j.factorial := by positivity
  have hnjf : (0:ℝ) < (n - j).factorial := by positivity
  -- |M| / n! ≤ (1/2)^(n-j) / j!
  have hM2 : |Mnj n j| / n.factorial ≤ (2⁻¹ : ℝ) ^ (n - j) / j.factorial := by
    rw [div_le_iff₀ hnf]
    refine hM.trans ?_
    rw [div_le_iff₀ (by positivity)]
    have e : (2⁻¹ : ℝ) ^ (n - j) / j.factorial * n.factorial *
        (j.factorial * ((j : ℝ) + 1) ^ (n - j + 1)) =
        (2⁻¹ : ℝ) ^ (n - j) * ((j : ℝ) + 1) ^ (n - j + 1) * n.factorial := by
      field_simp
    rw [e]
    have h2 : (1:ℝ) ≤ (2⁻¹ : ℝ) ^ (n - j) * ((j : ℝ) + 1) ^ (n - j + 1) := by
      rw [inv_pow]
      rw [le_inv_mul_iff₀ (by positivity)]; linarith
    calc ((n - j).factorial : ℝ) ≤ n.factorial / j.factorial := by
          rw [le_div_iff₀ hjf0]; linarith
      _ ≤ n.factorial := div_le_self hnf.le hjf
      _ ≤ _ := by nlinarith
  have hsplit : (2 * σ)⁻¹ ^ n * (2 * σ) ^ j = σ⁻¹ ^ (n - j) * (2⁻¹ : ℝ) ^ (n - j) := by
    rw [← mul_pow, ← mul_inv, mul_comm σ 2]
    rw [show n = (n - j) + j by omega, pow_add]
    rw [show n - j + j - j = n - j by omega]
    rw [mul_assoc, ← mul_pow, inv_mul_cancel₀ (by positivity), one_pow, mul_one]
  calc |Hc A B (n - j)| * |Mnj n j| / n.factorial = |Hc A B (n - j)| * (|Mnj n j| / n.factorial) := by
        ring
    _ ≤ (C * σ⁻¹ ^ (n - j)) * ((2⁻¹ : ℝ) ^ (n - j) / j.factorial) :=
        mul_le_mul hH hM2 (by positivity) (by positivity)
    _ = C * ((2 * σ)⁻¹ ^ n * (2 * σ) ^ j) / j.factorial := by rw [hsplit]; ring
    _ = _ := by ring

lemma stmt8_aux {A B : ℝ} (hA : 0 < A) (hB : 0 < B) {σ : ℝ} (hσ : 0 < σ) (hσρ : σ < ρ A B) :
    ∃ C : ℝ, ∀ n : ℕ, 2 ≤ n →
      |In A B n / (n.factorial : ℝ) - Hc A B n| ≤ C * (2 * σ)⁻¹ ^ n := by
  obtain ⟨C, hC⟩ := Y06.proof_6 A B hA hB σ hσ hσρ
  have hC0 : 0 ≤ C := by have := hC 0; simp at this; linarith [abs_nonneg (Hc A B 0)]
  refine ⟨C * (4 * S2 + Real.exp (2 * σ)), fun n hn => ?_⟩
  have hS2 : 0 ≤ S2 := tsum_nonneg (fun m => by positivity)
  have hnf : (0:ℝ) < n.factorial := by positivity
  rw [In_sum hA hB, Finset.sum_range_succ']
  simp only [Nat.sub_zero]
  have e : (∑ j ∈ Finset.range n, Hc A B (n - (j + 1)) * Mnj n (j + 1) + Hc A B n * Mnj n 0) /
      (n.factorial : ℝ) - Hc A B n =
      Hc A B n * (Mnj n 0 / n.factorial - 1) +
        ∑ j ∈ Finset.range n, Hc A B (n - (j + 1)) * Mnj n (j + 1) / n.factorial := by
    rw [← Finset.sum_div]; field_simp; ring
  rw [e]
  have t1 : |Hc A B n * (Mnj n 0 / n.factorial - 1)| ≤ C * (2 * σ)⁻¹ ^ n * (4 * S2) := by
    rw [abs_mul]
    calc |Hc A B n| * |Mnj n 0 / n.factorial - 1| ≤ (C * σ⁻¹ ^ n) * (4 * (2⁻¹ : ℝ) ^ n * S2) :=
          mul_le_mul (hC n) (etaD_sub_bound n hn) (abs_nonneg _) (by positivity)
      _ = _ := by rw [mul_inv, mul_pow]; ring
  have t2 : |∑ j ∈ Finset.range n, Hc A B (n - (j + 1)) * Mnj n (j + 1) / n.factorial| ≤
      C * (2 * σ)⁻¹ ^ n * Real.exp (2 * σ) := by
    refine (Finset.abs_sum_le_sum_abs _ _).trans ?_
    calc ∑ j ∈ Finset.range n, |Hc A B (n - (j + 1)) * Mnj n (j + 1) / n.factorial|
        ≤ ∑ j ∈ Finset.range n, C * (2 * σ)⁻¹ ^ n * ((2 * σ) ^ (j + 1) / (j + 1).factorial) := by
          refine Finset.sum_le_sum (fun j hj => ?_)
          rw [abs_div, abs_of_pos hnf]
          exact Mnj_term_bound hσ hC n (j + 1) (by omega)
            (by have := Finset.mem_range.mp hj; omega)
      _ = C * (2 * σ)⁻¹ ^ n * ∑ j ∈ Finset.range n, (2 * σ) ^ (j + 1) / (j + 1).factorial := by
          rw [Finset.mul_sum]
      _ ≤ C * (2 * σ)⁻¹ ^ n * Real.exp (2 * σ) := by
          gcongr
          have := Real.sum_le_exp_of_nonneg (by positivity : (0:ℝ) ≤ 2 * σ) (n + 1)
          rw [Finset.sum_range_succ'] at this
          simp only [pow_zero, Nat.factorial_zero, Nat.cast_one, div_one] at this
          linarith
  calc _ ≤ |Hc A B n * (Mnj n 0 / n.factorial - 1)| +
        |∑ j ∈ Finset.range n, Hc A B (n - (j + 1)) * Mnj n (j + 1) / n.factorial| := abs_add_le _ _
    _ ≤ C * (2 * σ)⁻¹ ^ n * (4 * S2) + C * (2 * σ)⁻¹ ^ n * Real.exp (2 * σ) := add_le_add t1 t2
    _ = _ := by ring

end Y06

theorem Y06.proof_8 : Y06.stmt_8 := fun _ _ hA hB _ hσ hσρ => Y06.stmt8_aux hA hB hσ hσρ

namespace Y06
open MeasureTheory Set Filter Topology

/-- Geometric upper bounds at every rate `σ < R`. -/
def UB (a : ℕ → ℝ) (R : ℝ) : Prop := ∀ σ : ℝ, 0 < σ → σ < R → ∃ C : ℝ, ∀ m : ℕ, |a m| ≤ C * σ⁻¹ ^ m

/-- No geometric bound at a rate `τ < 1/R`. -/
def NoBetter (a : ℕ → ℝ) (R : ℝ) : Prop :=
  ¬ ∃ K τ : ℝ, 0 < τ ∧ τ < 1 / R ∧ ∀ m : ℕ, |a m| ≤ K * τ ^ m

lemma root_le {x C s : ℝ} (hx : 0 ≤ x) (hs : 0 < s) {m : ℕ} (hm : m ≠ 0) (h : x ≤ C * s ^ m) :
    x ^ (1 / (m : ℝ)) ≤ C ^ (1 / (m : ℝ)) * s := by
  have hC : 0 ≤ C := by
    by_contra hC; push Not at hC
    have : C * s ^ m < 0 := mul_neg_of_neg_of_pos hC (by positivity)
    linarith
  calc x ^ (1 / (m : ℝ)) ≤ (C * s ^ m) ^ (1 / (m : ℝ)) := Real.rpow_le_rpow hx h (by positivity)
    _ = C ^ (1 / (m : ℝ)) * s := by
      rw [Real.mul_rpow hC (by positivity), one_div, Real.pow_rpow_inv_natCast hs.le hm]

lemma tendsto_root_const {C : ℝ} (hC : 0 < C) :
    Tendsto (fun m : ℕ => C ^ (1 / (m : ℝ))) atTop (𝓝 1) := by
  have h := (tendsto_const_nhds (x := C)).rpow (tendsto_one_div_atTop_nhds_zero_nat (𝕜 := ℝ))
    (Or.inl hC.ne')
  simpa using h

lemma root_nonneg (a : ℕ → ℝ) (m : ℕ) : 0 ≤ |a m| ^ (1 / (m : ℝ)) :=
  Real.rpow_nonneg (abs_nonneg _) _

/-- For every `y > 1/σ`... eventually `|a m|^{1/m} ≤ y` given a bound at rate `σ`. -/
lemma ev_root_le {a : ℕ → ℝ} {σ C y : ℝ} (hσ : 0 < σ) (hC : ∀ m : ℕ, |a m| ≤ C * σ⁻¹ ^ m)
    (hy : σ⁻¹ < y) : ∀ᶠ m in atTop, |a m| ^ (1 / (m : ℝ)) ≤ y := by
  set C' := max C 1 with hC'
  have hC'1 : 1 ≤ C' := le_max_right _ _
  have hb : ∀ m : ℕ, |a m| ≤ C' * σ⁻¹ ^ m := fun m =>
    (hC m).trans (mul_le_mul_of_nonneg_right (le_max_left _ _) (by positivity))
  have hlt : 1 < y * σ := by
    have := mul_lt_mul_of_pos_right hy hσ
    rwa [inv_mul_cancel₀ hσ.ne'] at this
  have hev := (tendsto_root_const (by linarith : (0:ℝ) < C')).eventually (gt_mem_nhds hlt)
  filter_upwards [hev, eventually_ge_atTop 1] with m hm hm1
  have := root_le (abs_nonneg (a m)) (inv_pos.mpr hσ) (by omega : m ≠ 0) (hb m)
  calc |a m| ^ (1 / (m : ℝ)) ≤ C' ^ (1 / (m : ℝ)) * σ⁻¹ := this
    _ ≤ (y * σ) * σ⁻¹ := by gcongr
    _ = y := by field_simp

lemma limsup_of {a : ℕ → ℝ} {R : ℝ} (hR : 0 < R) (h1 : UB a R) (h2 : NoBetter a R) :
    IsBoundedUnder (· ≤ ·) atTop (fun m : ℕ => |a m| ^ (1 / (m : ℝ))) ∧
    limsup (fun m : ℕ => |a m| ^ (1 / (m : ℝ))) atTop = 1 / R := by
  have hcob : IsCoboundedUnder (· ≤ ·) atTop (fun m : ℕ => |a m| ^ (1 / (m : ℝ))) :=
    isCoboundedUnder_le_of_le atTop (root_nonneg a)
  obtain ⟨C, hC⟩ := h1 (R / 2) (by positivity) (by linarith)
  have hbdd : IsBoundedUnder (· ≤ ·) atTop (fun m : ℕ => |a m| ^ (1 / (m : ℝ))) :=
    isBoundedUnder_of_eventually_le (ev_root_le (by positivity) hC (lt_add_one _))
  refine ⟨hbdd, le_antisymm ?_ ?_⟩
  · apply le_of_forall_gt_imp_ge_of_dense
    intro y hy
    have hy0 : 0 < y := lt_trans (by positivity) hy
    obtain ⟨σ, hσ1, hσ2⟩ : ∃ σ, y⁻¹ < σ ∧ σ < R := exists_between (by
      rw [one_div] at hy; exact (inv_lt_comm₀ hR hy0).mp hy)
    have hσ0 : 0 < σ := lt_trans (inv_pos.mpr hy0) hσ1
    obtain ⟨C, hC⟩ := h1 σ hσ0 hσ2
    exact limsup_le_of_le hcob (ev_root_le hσ0 hC ((inv_lt_comm₀ hy0 hσ0).mp hσ1))
  · by_contra hlt; push Not at hlt
    set ℓ := limsup (fun m : ℕ => |a m| ^ (1 / (m : ℝ))) atTop
    set τ := (max ℓ 0 + 1 / R) / 2 with hτ
    have hRi : 0 < 1 / R := by positivity
    have hτ0 : 0 < τ := by rw [hτ]; have := le_max_right ℓ 0; linarith
    have hτR : τ < 1 / R := by
      rw [hτ]; have : max ℓ 0 < 1 / R := max_lt hlt hRi; linarith
    have hℓτ : ℓ < τ := by rw [hτ]; have := le_max_left ℓ 0; have : max ℓ 0 < 1 / R := max_lt hlt hRi; linarith
    have hev := eventually_lt_of_limsup_lt hℓτ hbdd
    obtain ⟨N, hN⟩ := eventually_atTop.mp (hev.and (eventually_ge_atTop 1))
    apply h2
    refine ⟨1 + ∑ m ∈ Finset.range N, |a m| / τ ^ m, τ, hτ0, hτR, fun m => ?_⟩
    have hs : 0 ≤ ∑ m ∈ Finset.range N, |a m| / τ ^ m :=
      Finset.sum_nonneg (fun i _ => by positivity)
    by_cases hm : m < N
    · have : |a m| / τ ^ m ≤ ∑ m ∈ Finset.range N, |a m| / τ ^ m :=
        Finset.single_le_sum (f := fun i => |a i| / τ ^ i) (fun i _ => by positivity)
          (Finset.mem_range.mpr hm)
      have hp : 0 < τ ^ m := by positivity
      calc |a m| = |a m| / τ ^ m * τ ^ m := by field_simp
        _ ≤ (1 + ∑ m ∈ Finset.range N, |a m| / τ ^ m) * τ ^ m := by gcongr; linarith
    · push Not at hm
      obtain ⟨h1m, h2m⟩ := hN m hm
      have hm0 : m ≠ 0 := by omega
      have : |a m| ≤ τ ^ m := by
        rw [← Real.rpow_inv_natCast_pow (abs_nonneg (a m)) hm0]
        rw [one_div] at h1m
        exact pow_le_pow_left₀ (Real.rpow_nonneg (abs_nonneg _) _) h1m.le m
      calc |a m| ≤ 1 * τ ^ m := by rw [one_mul]; exact this
        _ ≤ _ := by gcongr; linarith

lemma transfer {a b : ℕ → ℝ} {R : ℝ} (hR : 0 < R) (ha1 : UB a R) (ha2 : NoBetter a R)
    (hab : ∀ σ : ℝ, 0 < σ → σ < R → ∃ C : ℝ, ∀ n : ℕ, 2 ≤ n → |b n - a n| ≤ C * (2 * σ)⁻¹ ^ n) :
    UB b R ∧ NoBetter b R := by
  constructor
  · intro σ hσ hσR
    obtain ⟨C, hC⟩ := ha1 σ hσ hσR
    obtain ⟨C', hC'⟩ := hab σ hσ hσR
    have hC0 : 0 ≤ C := by have := hC 0; simp at this; linarith [abs_nonneg (a 0)]
    refine ⟨C + |C'| + |b 0| + |b 1| * σ, fun n => ?_⟩
    have hsi : 0 < σ⁻¹ := inv_pos.mpr hσ
    rcases (show n = 0 ∨ n = 1 ∨ 2 ≤ n by omega) with h | h | h
    · subst h; simp; nlinarith [abs_nonneg C', abs_nonneg (b 1)]
    · subst h; simp
      have : |b 1| = |b 1| * σ * σ⁻¹ := by field_simp
      nlinarith [abs_nonneg C', abs_nonneg (b 0), abs_nonneg (b 1), mul_pos hσ hsi]
    · have h2 : (2 * σ)⁻¹ ^ n ≤ σ⁻¹ ^ n := by
        apply pow_le_pow_left₀ (by positivity)
        rw [inv_le_inv₀ (by positivity) hσ]; linarith
      have hb := hC' n h
      have e : |b n| ≤ |a n| + |b n - a n| := by
        have := abs_add_le (a n) (b n - a n); simp at this; linarith
      have hpn : 0 ≤ σ⁻¹ ^ n := by positivity
      calc |b n| ≤ C * σ⁻¹ ^ n + |C'| * σ⁻¹ ^ n := by
            have : C' * (2 * σ)⁻¹ ^ n ≤ |C'| * σ⁻¹ ^ n :=
              mul_le_mul (le_abs_self C') h2 (by positivity) (abs_nonneg _)
            linarith [hC n]
        _ ≤ _ := by
          have h0 : 0 ≤ (|b 0| + |b 1| * σ) * σ⁻¹ ^ n := by positivity
          nlinarith [h0]
  · rintro ⟨K, τ, hτ0, hτR, hK⟩
    obtain ⟨C', hC'⟩ := hab (3 * R / 4) (by positivity) (by linarith)
    set τ' := max τ (2 * (3 * R / 4))⁻¹ with hτ'
    have hτ'0 : 0 < τ' := lt_max_of_lt_left hτ0
    have hτ'R : τ' < 1 / R := by
      refine max_lt hτR ?_
      rw [one_div, inv_lt_inv₀ (by positivity) hR]; linarith
    apply ha2
    refine ⟨|K| + |C'| + |a 0| + |a 1| / τ', τ', hτ'0, hτ'R, fun n => ?_⟩
    rcases (show n = 0 ∨ n = 1 ∨ 2 ≤ n by omega) with h | h | h
    · subst h; simp; nlinarith [abs_nonneg C', abs_nonneg K, abs_nonneg (a 1), div_nonneg (abs_nonneg (a 1)) hτ'0.le]
    · subst h; simp
      have : |a 1| / τ' * τ' = |a 1| := by field_simp
      nlinarith [abs_nonneg C', abs_nonneg K, abs_nonneg (a 0), mul_pos hτ'0 hτ'0]
    · have h1 : τ ^ n ≤ τ' ^ n := pow_le_pow_left₀ hτ0.le (le_max_left _ _) n
      have h2 : (2 * (3 * R / 4))⁻¹ ^ n ≤ τ' ^ n := pow_le_pow_left₀ (by positivity) (le_max_right _ _) n
      have hb := hC' n h
      have e : |a n| ≤ |b n| + |b n - a n| := by
        have := abs_sub_le (a n) (b n) 0; rw [abs_sub_comm (a n) (b n)] at this; simp at this
        have := abs_sub (b n) (b n - a n); simp at this; linarith
      have hpn : 0 ≤ τ' ^ n := by positivity
      have t1 : |b n| ≤ |K| * τ' ^ n :=
        (hK n).trans (mul_le_mul (le_abs_self K) h1 (by positivity) (abs_nonneg _))
      have t2 : |b n - a n| ≤ |C'| * τ' ^ n :=
        hb.trans (mul_le_mul (le_abs_self C') h2 (by positivity) (abs_nonneg _))
      nlinarith [abs_nonneg (a 0), abs_nonneg (a 1), div_nonneg (abs_nonneg (a 1)) hτ'0.le]

/-- Coefficients with `c n / n!` not better than rate `1/R` give a divergent series. -/
lemma not_summable_of {c : ℕ → ℝ} {R : ℝ} (hR : 0 < R)
    (h : NoBetter (fun n => c n / (n.factorial : ℝ)) R) {η : ℝ} (hη : η ≠ 0) :
    ¬ Summable (fun n : ℕ => c n * η ^ n) := by
  intro hs
  obtain ⟨M, hM⟩ := (isBounded_iff_forall_norm_le.mp
    (Metric.isBounded_range_of_tendsto _ hs.tendsto_atTop_zero))
  have hM' : ∀ n, |c n| * |η| ^ n ≤ M := fun n => by
    have := hM _ ⟨n, rfl⟩; rwa [Real.norm_eq_abs, abs_mul, abs_pow] at this
  have hη0 : 0 < |η| := abs_pos.mpr hη
  set τ := 1 / (2 * R) with hτ
  have hτ0 : 0 < τ := by positivity
  set x := (|η| * τ)⁻¹
  apply h
  refine ⟨M * Real.exp x, τ, hτ0, by rw [hτ, one_div_lt_one_div (by positivity) hR]; linarith,
    fun n => ?_⟩
  have hf : (0:ℝ) < n.factorial := by positivity
  have hx := Real.pow_div_factorial_le_exp x (by positivity) n
  have hM0 : 0 ≤ M := (mul_nonneg (abs_nonneg _) (by positivity)).trans (hM' 0)
  rw [abs_div, abs_of_pos hf]
  have hcn : |c n| ≤ M * (|η| ^ n)⁻¹ := by
    rw [le_mul_inv_iff₀ (by positivity)]; exact hM' n
  calc |c n| / n.factorial ≤ M * (|η| ^ n)⁻¹ / n.factorial := by gcongr
    _ = M * (x ^ n / n.factorial) * τ ^ n := by
        simp only [x]; rw [inv_pow, mul_pow]; field_simp
    _ ≤ M * Real.exp x * τ ^ n := by gcongr

end Y06

namespace Y06
open MeasureTheory Set Filter Topology

lemma HC_analytic {A B : ℝ} (hA : 0 < A) (hB : 0 < B) (hAB : A ≤ B) :
    AnalyticOnNhd ℂ (HC A B) (Udom A B) :=
  (HC_diff hA hB hAB).analyticOnNhd (Udom_open A B)

lemma HC_im_iter {A B : ℝ} (hA : 0 < A) (hB : 0 < B) (hAB : A ≤ B) (k : ℕ) (x : ℝ) :
    (iteratedDeriv k (HC A B) x).im = 0 := by
  have hfd : DifferentiableOn ℂ (fun z => -Complex.I * HC A B z) (Udom A B) :=
    (HC_diff hA hB hAB).const_mul _
  have hb := bridge (Udom_open A B) hfd isOpen_univ (fun x _ => real_mem_U A B x)
    (g := fun _ => (0:ℝ)) (fun x _ => by simp [HC_real hA hB]) k x (mem_univ x)
  have hca : ContDiffAt ℂ k (HC A B) (x : ℂ) :=
    (HC_analytic hA hB hAB _ (real_mem_U A B x)).contDiffAt
  rw [iteratedDeriv_const_mul (-Complex.I) hca] at hb
  have h0 : iteratedDeriv k (fun _ : ℝ => (0:ℝ)) x = 0 := by
    rcases k with _ | k
    · simp
    · simp
  rw [h0] at hb
  simpa [Complex.mul_re] using hb.symm

lemma HC_iter0 {A B : ℝ} (hA : 0 < A) (hB : 0 < B) (hAB : A ≤ B) (m : ℕ) :
    iteratedDeriv m (HC A B) 0 = (((m.factorial : ℝ) * Hc A B m : ℝ) : ℂ) := by
  have h1 := H_iter hA hB hAB m 0
  have h2 := HC_im_iter hA hB hAB m 0
  rw [Complex.ofReal_zero] at h1 h2
  apply Complex.ext
  · simp only [Complex.ofReal_re]; rw [← h1]; unfold Hc; field_simp
  · simp only [Complex.ofReal_im]; exact h2

lemma HC_taylor {A B : ℝ} (hA : 0 < A) (hB : 0 < B) (hAB : A ≤ B) {z : ℂ} (hz : ‖z‖ < ρ A B) :
    HasSum (fun n : ℕ => (Hc A B n : ℂ) * z ^ n) (HC A B z) := by
  have hd : DifferentiableOn ℂ (HC A B) (Metric.ball 0 (ρ A B)) :=
    (HC_diff hA hB hAB).mono (fun w hw => Or.inl (by simpa using hw))
  have := Complex.hasSum_taylorSeries_on_ball hd (by simpa using hz)
  convert this using 1
  funext n
  rw [HC_iter0 hA hB hAB, sub_zero]
  simp only [smul_eq_mul]
  push_cast
  have : (n.factorial : ℂ) ≠ 0 := by exact_mod_cast Nat.factorial_ne_zero n
  field_simp

lemma ext_series {A B K τ r : ℝ} (hr0 : 0 < r) (hτ : 0 < τ) (hrτ : r * τ < 1)
    (hK : ∀ m : ℕ, |Hc A B m| ≤ K * τ ^ m) :
    ∃ G : ℂ → ℂ, AnalyticOnNhd ℂ G (Metric.ball 0 r) ∧
      ∀ z : ℂ, ‖z‖ < r → HasSum (fun n : ℕ => (Hc A B n : ℂ) * z ^ n) (G z) := by
  set p := FormalMultilinearSeries.ofScalars ℂ (fun n => (Hc A B n : ℂ))
  have hrad : ((r.toNNReal : NNReal) : ENNReal) ≤ p.radius := by
    refine p.le_radius_of_bound |K| (fun n => ?_)
    rw [FormalMultilinearSeries.ofScalars_norm, Complex.norm_real, Real.norm_eq_abs]
    rw [Real.coe_toNNReal r hr0.le]
    have h1 : (τ * r) ^ n ≤ 1 := pow_le_one₀ (by positivity) (by linarith)
    calc |Hc A B n| * r ^ n ≤ K * τ ^ n * r ^ n := by gcongr; exact hK n
      _ = K * (τ * r) ^ n := by rw [mul_pow]; ring
      _ ≤ |K| * (τ * r) ^ n := by gcongr; exact le_abs_self K
      _ ≤ |K| * 1 := by gcongr
      _ = |K| := mul_one _
  have hpos : 0 < p.radius := lt_of_lt_of_le (by simpa using hr0) hrad
  have hps := p.hasFPowerSeriesOnBall hpos
  have hmem : ∀ z : ℂ, ‖z‖ < r → z ∈ Metric.eball (0:ℂ) p.radius := by
    intro z hz
    have : z ∈ Metric.eball (0:ℂ) ((r.toNNReal : NNReal) : ENNReal) := by
      rw [Metric.eball_coe]; simpa [Real.coe_toNNReal r hr0.le] using hz
    exact Metric.eball_subset_eball hrad this
  refine ⟨p.sum, fun z hz => hps.analyticAt_of_mem (hmem z (by simpa using hz)), fun z hz => ?_⟩
  have := hps.hasSum (hmem z hz)
  simp only [zero_add] at this
  convert this using 1
  funext n
  rw [FormalMultilinearSeries.ofScalars_apply_eq, smul_eq_mul]

end Y06

namespace Y06
open MeasureTheory Set Filter Topology

lemma slit_re {c : ℝ} (hc0 : 0 < c) (w : ℂ) (hw : w.re < -(Real.log c / 2)) :
    1 + (c : ℂ) * Complex.exp (2 * w) ∈ Complex.slitPlane := by
  rw [Complex.mem_slitPlane_iff]
  by_contra hcon
  push Not at hcon
  obtain ⟨hre, him⟩ := hcon
  have e2re : (2 * w).re = 2 * w.re := by simp
  have e2im : (2 * w).im = 2 * w.im := by simp
  simp only [Complex.add_re, Complex.add_im, Complex.one_re, Complex.one_im, Complex.mul_re,
    Complex.mul_im, Complex.ofReal_re, Complex.ofReal_im, Complex.exp_re, Complex.exp_im,
    e2re, e2im, zero_mul, sub_zero, zero_add] at hre him
  have hcos : Real.cos (2 * w.im) ≤ 1 := Real.cos_le_one _
  have hge : 1 ≤ c * Real.exp (2 * w.re) := by
    have hE := Real.exp_pos (2 * w.re)
    have : c * (Real.exp (2 * w.re) * Real.cos (2 * w.im)) ≤ c * Real.exp (2 * w.re) := by
      have := mul_le_mul_of_nonneg_left hcos (le_of_lt (mul_pos hc0 hE))
      nlinarith
    -- cos(2y) = -1 is not needed: 1 + c e^{2x} cos ≤ 0 gives c e^{2x} cos ≤ -1
    have h1 : c * (Real.exp (2 * w.re) * Real.cos (2 * w.im)) ≤ -1 := by linarith
    have h2 : -1 ≤ Real.cos (2 * w.im) := Real.neg_one_le_cos _
    nlinarith [mul_pos hc0 hE]
  have : c * Real.exp (2 * w.re) < 1 := by
    rw [← Real.exp_log hc0, ← Real.exp_add, Real.exp_lt_one_iff]
    linarith
  linarith

lemma half_log_cc {A B : ℝ} (_hA : 0 < A) (_hB : 0 < B) : Real.log (cc A B) / 2 = -L A B := by
  unfold cc L
  rw [Real.log_pow, show A / B = (B / A)⁻¹ by rw [inv_div], Real.log_inv]; push_cast; ring

/-- The branch `U = (B/r) e^{-z} √(1 + c e^{2z})` of `√R₁/r`. -/
noncomputable def U (A B : ℝ) (z : ℂ) : ℂ :=
  ((B / X01.r A B : ℝ) : ℂ) * (Complex.exp (-z) * csq (1 + (cc A B : ℂ) * Complex.exp (2 * z)))

lemma U_sq {A B : ℝ} (hA : 0 < A) (hB : 0 < B) (z : ℂ) :
    U A B z ^ 2 = R1 A B z / (X01.r A B : ℂ) ^ 2 := by
  have hr : (X01.r A B : ℂ) ≠ 0 := by
    exact_mod_cast (Real.sqrt_pos.mpr (by positivity : 0 < A ^ 2 + B ^ 2)).ne'
  have hB' : (B : ℂ) ≠ 0 := by exact_mod_cast hB.ne'
  have hE : Complex.exp (-z) ≠ 0 := Complex.exp_ne_zero _
  have h1 : Complex.exp (2 * z) = (Complex.exp (-z) ^ 2)⁻¹ := by
    rw [← Complex.exp_nat_mul, ← Complex.exp_neg]; ring_nf
  have h2 : Complex.exp (-2 * z) = Complex.exp (-z) ^ 2 := by
    rw [← Complex.exp_nat_mul]; ring_nf
  unfold U R1 cc
  rw [mul_pow, mul_pow, csq_sq, h1, h2]
  push_cast
  field_simp
  ring

/-- The singular point `z₁ = L + iπ/2`. -/
noncomputable def z1 (A B : ℝ) : ℂ := (L A B : ℂ) + Complex.I * ((Real.pi / 2 : ℝ) : ℂ)

lemma R1_z1 {A B : ℝ} (hA : 0 < A) (hB : 0 < B) : R1 A B (z1 A B) = 0 := by
  rw [R1_zero_iff A B hA hB]
  exact ⟨0, by unfold z1; push_cast; ring⟩

lemma z1_re (A B : ℝ) : (z1 A B).re = L A B := by simp [z1]

lemma norm_z1 (A B : ℝ) : ‖z1 A B‖ = ρ A B := by
  have h : ‖z1 A B‖ ^ 2 = ρ A B ^ 2 := by
    rw [Complex.sq_norm, Complex.normSq_apply, ρ_sq]
    simp [z1]; ring
  have := ρ_pos A B
  nlinarith [norm_nonneg (z1 A B), sq_nonneg (‖z1 A B‖ - ρ A B), sq_nonneg (‖z1 A B‖ + ρ A B)]

lemma final_contra {A B : ℝ} (hA : 0 < A) (hB : 0 < B) {Φ : ℂ → ℂ} {D : Set ℂ}
    (hDo : IsOpen D) (hDc : Convex ℝ D) (hz1 : z1 A B ∈ D) (hΦ : AnalyticOnNhd ℂ Φ D)
    {z0 : ℂ} (hz0 : z0 ∈ D)
    (heq : ∀ᶠ z in 𝓝 z0, Φ z ^ 2 = R1 A B z / (X01.r A B : ℂ) ^ 2) : False := by
  have hdR : Differentiable ℂ (fun z => R1 A B z / (X01.r A B : ℂ) ^ 2) := by
    unfold R1; fun_prop
  have hRan : AnalyticOnNhd ℂ (fun z => R1 A B z / (X01.r A B : ℂ) ^ 2) D :=
    fun z _ => hdR.analyticAt z
  have hsq : AnalyticOnNhd ℂ (fun z => Φ z ^ 2) D := hΦ.pow 2
  have heqOn := hsq.eqOn_of_preconnected_of_eventuallyEq hRan hDc.isPreconnected hz0 heq
  have hv : Φ (z1 A B) = 0 := by
    have := heqOn hz1
    simp only [R1_z1 hA hB, zero_div] at this
    exact pow_eq_zero_iff (two_ne_zero) |>.mp this
  have hev1 : (fun z => Φ z ^ 2) =ᶠ[𝓝 (z1 A B)] (fun z => R1 A B z / (X01.r A B : ℂ) ^ 2) :=
    Filter.eventually_of_mem (hDo.mem_nhds hz1) heqOn
  have h1' : deriv (fun z => Φ z ^ 2) (z1 A B) = 0 := by
    have := ((hΦ _ hz1).differentiableAt.hasDerivAt).pow 2
    rw [show (fun z => Φ z ^ 2) = Φ ^ 2 from rfl, this.deriv, hv]; simp
  have hd2 := (hasDerivAt_R1 A B (z1 A B)).div_const ((X01.r A B : ℂ) ^ 2)
  have := hev1.deriv_eq
  rw [h1', hd2.deriv] at this
  have hz := R1_z1 hA hB
  unfold R1 at hz
  have hr : (X01.r A B : ℂ) ≠ 0 := by
    exact_mod_cast (Real.sqrt_pos.mpr (by positivity : 0 < A ^ 2 + B ^ 2)).ne'
  have hA' : (A : ℂ) ≠ 0 := by exact_mod_cast hA.ne'
  have e : (B : ℂ) ^ 2 * (Complex.exp (-2 * z1 A B) * (-2)) = 2 * (A : ℂ) ^ 2 := by
    linear_combination (-2 : ℂ) * hz
  rw [e] at this
  exact div_ne_zero (mul_ne_zero two_ne_zero (pow_ne_zero 2 hA')) (pow_ne_zero 2 hr) this.symm

/-- A point of the disc `D(0,ρ)` close to `z₁`. -/
lemma near_z1 (A B : ℝ) {δ : ℝ} (hδ : 0 < δ) (hδρ : δ ≤ ρ A B) :
    ∃ z0 : ℂ, ‖z0‖ < ρ A B ∧ ‖z0 - z1 A B‖ < δ := by
  have hρ := ρ_pos A B
  set t : ℝ := 1 - δ / (2 * ρ A B)
  have ht0 : 0 ≤ t := by
    have : δ / (2 * ρ A B) ≤ 1 / 2 := by rw [div_le_iff₀ (by positivity)]; linarith
    simp only [t]; linarith
  refine ⟨(t : ℂ) * z1 A B, ?_, ?_⟩
  · rw [norm_mul, Complex.norm_real, Real.norm_of_nonneg ht0, norm_z1]
    have : 0 < δ / (2 * ρ A B) := by positivity
    simp only [t]; nlinarith
  · rw [show (t : ℂ) * z1 A B - z1 A B = ((t - 1 : ℝ) : ℂ) * z1 A B by push_cast; ring, norm_mul,
      Complex.norm_real, norm_z1]
    simp only [t, sub_sub_cancel_left, Real.norm_eq_abs, abs_neg]
    rw [abs_of_pos (by positivity)]
    field_simp
    linarith

end Y06

namespace Y06
open MeasureTheory Set Filter Topology

lemma csq_diffAt {c : ℝ} {z : ℂ} (h : 1 + (c : ℂ) * Complex.exp (2 * z) ∈ Complex.slitPlane) :
    DifferentiableAt ℂ (fun z : ℂ => csq (1 + (c : ℂ) * Complex.exp (2 * z))) z := by
  unfold csq; exact DifferentiableAt.cpow_const (by fun_prop) h

lemma csq_neg_diffAt {c : ℝ} {z : ℂ} (h : 1 + (c : ℂ) * Complex.exp (2 * (-z)) ∈ Complex.slitPlane) :
    DifferentiableAt ℂ (fun z : ℂ => csq (1 + (c : ℂ) * Complex.exp (2 * (-z)))) z := by
  unfold csq; exact DifferentiableAt.cpow_const (by fun_prop) h

/-- Common setup: the series extension `G` of `H` beyond radius `ρ`. -/
lemma setup {A B : ℝ} (hA : 0 < A) (hB : 0 < B) (hAB : A ≤ B) {K τ : ℝ} (hτ0 : 0 < τ)
    (hτρ : τ < 1 / ρ A B) (hK : ∀ m : ℕ, |Hc A B m| ≤ K * τ ^ m) :
    ∃ r2 : ℝ, ρ A B < r2 ∧ ∃ G : ℂ → ℂ, AnalyticOnNhd ℂ G (Metric.ball 0 r2) ∧
      ∀ z : ℂ, ‖z‖ < ρ A B → G z = HC A B z := by
  have hρ := ρ_pos A B
  have hρτ : ρ A B < 1 / τ := (lt_one_div hτ0 hρ).mp hτρ
  set r2 := (ρ A B + 1 / τ) / 2 with hr2
  have hr2a : ρ A B < r2 := by rw [hr2]; linarith
  have hr2b : r2 * τ < 1 := by
    have : r2 < 1 / τ := by rw [hr2]; linarith
    rwa [lt_div_iff₀ hτ0] at this
  obtain ⟨G, hGan, hGsum⟩ := ext_series (by linarith : 0 < r2) hτ0 hr2b hK
  exact ⟨r2, hr2a, G, hGan, fun z hz =>
    (hGsum z (by linarith)).unique (HC_taylor hA hB hAB hz)⟩

lemma ball_sub {A B : ℝ} {r2 δ : ℝ} (hδ : δ ≤ r2 - ρ A B) :
    Metric.ball (z1 A B) δ ⊆ Metric.ball 0 r2 := by
  intro z hz
  rw [Metric.mem_ball, dist_eq_norm] at hz
  rw [Metric.mem_ball, dist_zero_right]
  have := norm_le_norm_add_norm_sub' z (z1 A B)
  rw [norm_z1] at this
  have : ‖z - z1 A B‖ = ‖z1 A B - z‖ := norm_sub_rev _ _
  have := norm_le_norm_add_norm_sub' z (z1 A B)
  rw [norm_z1] at this
  linarith [norm_sub_rev z (z1 A B)]

lemma sing_pos {A B : ℝ} (hA : 0 < A) (hB : 0 < B) (hAB : A < B) {K τ : ℝ} (hτ0 : 0 < τ)
    (hτρ : τ < 1 / ρ A B) (hK : ∀ m : ℕ, |Hc A B m| ≤ K * τ ^ m) : False := by
  obtain ⟨r2, hr2, G, hGan, hGHC⟩ := setup hA hB hAB.le hτ0 hτρ hK
  have hρ := ρ_pos A B
  have hL : 0 < L A B := Real.log_pos ((one_lt_div hA).mpr hAB)
  have hLρ : L A B ≤ ρ A B := by nlinarith [ρ_sq A B, Real.pi_pos]
  set v : ℂ → ℂ := fun z =>
    ((B / X01.r A B : ℝ) : ℂ) * csq (1 + (cc A B : ℂ) * Complex.exp (2 * (-z))) with hvdef
  set W : Set ℂ := {z : ℂ | -L A B < z.re}
  have hWo : IsOpen W := isOpen_lt continuous_const Complex.continuous_re
  have hv : DifferentiableOn ℂ v W := by
    intro z hz
    have hs : 1 + (cc A B : ℂ) * Complex.exp (2 * (-z)) ∈ Complex.slitPlane := by
      apply slit_re (cc_pos hA hB)
      rw [half_log_cc hA hB, Complex.neg_re]
      have : -L A B < z.re := hz
      linarith
    exact ((csq_neg_diffAt hs).const_mul _).differentiableWithinAt
  set δ := min (r2 - ρ A B) (L A B) with hδ
  have hδ0 : 0 < δ := lt_min (by linarith) hL
  have hδρ : δ ≤ ρ A B := (min_le_right _ _).trans hLρ
  set D := Metric.ball (z1 A B) δ
  have hD1 : D ⊆ Metric.ball 0 r2 := ball_sub (min_le_left _ _)
  have hDW : D ⊆ W := by
    intro z hz
    rw [Metric.mem_ball, dist_eq_norm] at hz
    have h1 := Complex.abs_re_le_norm (z - z1 A B)
    rw [Complex.sub_re, z1_re] at h1
    have h2 := (abs_le.mp (h1.trans hz.le)).1
    show -L A B < z.re
    have : δ ≤ L A B := min_le_right _ _
    linarith
  have hΦ : AnalyticOnNhd ℂ (fun z => G z - v z) D := fun z hz =>
    (hGan z (hD1 hz)).sub ((hv.analyticOnNhd hWo) z (hDW hz))
  obtain ⟨z0, hz0a, hz0b⟩ := near_z1 A B hδ0 hδρ
  have hz0D : z0 ∈ D := by rw [Metric.mem_ball, dist_eq_norm]; exact hz0b
  have hN : IsOpen (D ∩ {z : ℂ | ‖z‖ < ρ A B}) :=
    Metric.isOpen_ball.inter (isOpen_lt continuous_norm continuous_const)
  have heq : ∀ᶠ z in 𝓝 z0, (G z - v z) ^ 2 = R1 A B z / (X01.r A B : ℂ) ^ 2 := by
    filter_upwards [hN.mem_nhds ⟨hz0D, hz0a⟩] with z hz
    rw [hGHC z hz.2, ← U_sq hA hB]
    congr 1
    simp only [hvdef, HC, U]; ring
  exact final_contra hA hB Metric.isOpen_ball (convex_ball _ _) (Metric.mem_ball_self hδ0) hΦ
    hz0D heq

lemma cc_self {A : ℝ} (hA : 0 < A) : cc A A = 1 := by unfold cc; rw [div_self hA.ne', one_pow]

lemma sing_zero {A : ℝ} (hA : 0 < A) {K τ : ℝ} (hτ0 : 0 < τ)
    (hτρ : τ < 1 / ρ A A) (hK : ∀ m : ℕ, |Hc A A m| ≤ K * τ ^ m) : False := by
  obtain ⟨r2, hr2, G, hGan, hGHC⟩ := setup hA hA le_rfl hτ0 hτρ hK
  have hρ := ρ_pos A A
  set u : ℂ → ℂ := fun z => Complex.exp (-z) * csq (1 + (cc A A : ℂ) * Complex.exp (2 * z)) with hu
  set v : ℂ → ℂ := fun z => csq (1 + (cc A A : ℂ) * Complex.exp (2 * (-z))) with hv
  have hud : DifferentiableOn ℂ u (Udom A A) := fun z hz =>
    ((by fun_prop : DifferentiableAt ℂ (fun z : ℂ => Complex.exp (-z)) z).mul
      (csq_diffAt (slit_U hA hA le_rfl hz))).differentiableWithinAt
  have hvd : DifferentiableOn ℂ v (Udom A A) := fun z hz =>
    (csq_neg_diffAt (slit_U hA hA le_rfl (neg_U hz))).differentiableWithinAt
  have hball : Metric.ball (0:ℂ) (ρ A A) ⊆ Udom A A := fun w hw => Or.inl (by simpa using hw)
  have hfan : AnalyticOnNhd ℂ (fun z => u z - v z) (Metric.ball 0 (ρ A A)) :=
    ((hud.sub hvd).analyticOnNhd (Udom_open A A)).mono hball
  have hsq : ∀ z, u z ^ 2 = v z ^ 2 := by
    intro z
    simp only [hu, hv, mul_pow, csq_sq, cc_self hA]
    rw [← Complex.exp_nat_mul]
    push_cast
    have e : Complex.exp (2 * -z) * Complex.exp (2 * z) = 1 := by rw [← Complex.exp_add]; simp
    linear_combination e
  have hg0 : u 0 + v 0 ≠ 0 := by
    simp only [hu, hv, cc_self hA, neg_zero, mul_zero, Complex.exp_zero, one_mul, Complex.ofReal_one]
    rw [show (1 : ℂ) + 1 = ((2 : ℝ) : ℂ) by push_cast; ring, csq_ofReal (by norm_num)]
    have : (0:ℝ) < Real.sqrt 2 := Real.sqrt_pos.mpr (by norm_num)
    norm_cast; linarith
  have hgc : ContinuousAt (fun z => u z + v z) 0 :=
    ((hud.add hvd).differentiableAt ((Udom_open A A).mem_nhds (real_mem_U A A 0 |>.out
      |> fun h => by simpa using real_mem_U A A 0))).continuousAt
  have hev : (fun z => u z - v z) =ᶠ[𝓝 0] fun _ => 0 := by
    filter_upwards [hgc.eventually_ne hg0] with z hz
    have : (u z - v z) * (u z + v z) = 0 := by ring_nf; rw [hsq z]; ring
    exact (mul_eq_zero.mp this).resolve_right hz
  have huv := hfan.eqOn_of_preconnected_of_eventuallyEq analyticOnNhd_const
    (convex_ball _ _).isPreconnected (Metric.mem_ball_self hρ) hev
  set δ := min (r2 - ρ A A) (ρ A A) with hδ
  have hδ0 : 0 < δ := lt_min (by linarith) hρ
  set D := Metric.ball (z1 A A) δ
  have hD1 : D ⊆ Metric.ball 0 r2 := ball_sub (min_le_left _ _)
  have hΦ : AnalyticOnNhd ℂ (fun z => G z / 2) D := fun z hz =>
    (hGan z (hD1 hz)).div_const
  obtain ⟨z0, hz0a, hz0b⟩ := near_z1 A A hδ0 (min_le_right _ _)
  have hz0D : z0 ∈ D := by rw [Metric.mem_ball, dist_eq_norm]; exact hz0b
  have hN : IsOpen (D ∩ {z : ℂ | ‖z‖ < ρ A A}) :=
    Metric.isOpen_ball.inter (isOpen_lt continuous_norm continuous_const)
  have heq : ∀ᶠ z in 𝓝 z0, (G z / 2) ^ 2 = R1 A A z / (X01.r A A : ℂ) ^ 2 := by
    filter_upwards [hN.mem_nhds ⟨hz0D, hz0a⟩] with z hz
    have h1 : u z = v z := by
      have := huv (show z ∈ Metric.ball (0:ℂ) (ρ A A) by simpa using hz.2)
      simp only at this; exact sub_eq_zero.mp this
    rw [hGHC z hz.2, ← U_sq hA hA]
    congr 1
    have : HC A A z = ((A / X01.r A A : ℝ) : ℂ) * (u z + v z) := by simp only [HC, hu, hv]
    rw [this, ← h1]
    simp only [U, hu]; ring
  exact final_contra hA hA Metric.isOpen_ball (convex_ball _ _) (Metric.mem_ball_self hδ0) hΦ
    hz0D heq

lemma noBetter_aux {A B : ℝ} (hA : 0 < A) (hB : 0 < B) (hAB : A ≤ B) :
    NoBetter (Hc A B) (ρ A B) := by
  rintro ⟨K, τ, hτ0, hτρ, hK⟩
  rcases hAB.lt_or_eq with h | h
  · exact sing_pos hA hB h hτ0 hτρ hK
  · subst h; exact sing_zero hA hτ0 hτρ hK

lemma noBetter_Hc {A B : ℝ} (hA : 0 < A) (hB : 0 < B) : NoBetter (Hc A B) (ρ A B) := by
  rcases le_total A B with hAB | hBA
  · exact noBetter_aux hA hB hAB
  · rw [← Hc_swap, ← ρ_swap]; exact noBetter_aux hB hA hBA

lemma UB_Hc {A B : ℝ} (hA : 0 < A) (hB : 0 < B) : UB (Hc A B) (ρ A B) :=
  fun σ hσ hσρ => Y06.proof_6 A B hA hB σ hσ hσρ

end Y06

theorem Y06.proof_7 : Y06.stmt_7 := fun _ _ hA hB =>
  Y06.limsup_of (Y06.ρ_pos _ _) (Y06.UB_Hc hA hB) (Y06.noBetter_Hc hA hB)

namespace Y06

lemma In_growth {A B : ℝ} (hA : 0 < A) (hB : 0 < B) :
    UB (fun n => In A B n / (n.factorial : ℝ)) (ρ A B) ∧
    NoBetter (fun n => In A B n / (n.factorial : ℝ)) (ρ A B) :=
  transfer (ρ_pos A B) (UB_Hc hA hB) (noBetter_Hc hA hB)
    (fun σ hσ hσρ => Y06.proof_8 A B hA hB σ hσ hσρ)

end Y06

theorem Y06.proof_9 : Y06.stmt_9 := fun _ _ hA hB =>
  Y06.limsup_of (Y06.ρ_pos _ _) (Y06.In_growth hA hB).1 (Y06.In_growth hA hB).2

theorem Y06.proof_10 : Y06.stmt_10 := fun _ _ hA hB _ hη =>
  Y06.not_summable_of (Y06.ρ_pos _ _) (Y06.In_growth hA hB).2 hη


namespace Y06
open MeasureTheory Set Filter Topology

/-- Iterated derivatives of a power series at its centre. -/
lemma iter_of_hasFPS {f : ℝ → ℝ} {p : FormalMultilinearSeries ℝ ℝ ℝ} {x : ℝ} {r : ENNReal}
    (h : HasFPowerSeriesOnBall f p x r) (n : ℕ) :
    iteratedDeriv n f x = (n.factorial : ℝ) * p.coeff n := by
  have := h.factorial_smul 1 n
  rw [iteratedDeriv_eq_iteratedFDeriv, ← this, FormalMultilinearSeries.apply_eq_pow_smul_coeff]
  simp [nsmul_eq_mul]

/-- Taylor coefficients of a real analytic function decay geometrically. -/
lemma analytic_coeff {f : ℝ → ℝ} {x : ℝ} (hf : AnalyticAt ℝ f x) :
    ∃ C s : ℝ, 0 < s ∧ ∀ n : ℕ, |iteratedDeriv n f x / (n.factorial : ℝ)| ≤ C * s ^ n := by
  obtain ⟨p, r, hp⟩ := hf
  obtain ⟨r', hr'0, hr'r⟩ := ENNReal.lt_iff_exists_nnreal_btwn.mp hp.r_pos
  obtain ⟨C, _, hC⟩ := p.norm_mul_pow_le_of_lt_radius (lt_of_lt_of_le hr'r hp.r_le)
  have hr' : (0:ℝ) < r' := by exact_mod_cast hr'0
  refine ⟨C, (r' : ℝ)⁻¹, inv_pos.mpr hr', fun n => ?_⟩
  rw [iter_of_hasFPS hp n]
  have hf : (n.factorial : ℝ) ≠ 0 := by positivity
  rw [mul_div_cancel_left₀ _ hf, ← Real.norm_eq_abs, ← FormalMultilinearSeries.norm_apply_eq_norm_coef]
  have := hC n
  rw [inv_pow, le_mul_inv_iff₀ (by positivity)]
  exact this

/-- `I_n` cannot grow only like `(n+1) X^n`. -/
lemma no_geom_In {A B : ℝ} (hA : 0 < A) (hB : 0 < B) {C X : ℝ} (hX : 0 ≤ X)
    (h : ∀ n : ℕ, |In A B n| ≤ C * ((n : ℝ) + 1) * X ^ n) : False := by
  apply (In_growth hA hB).2
  have hρ := ρ_pos A B
  set τ := 1 / (2 * ρ A B) with hτ
  have hτ0 : 0 < τ := by positivity
  set y := 2 * X / τ
  have hy : 0 ≤ y := by positivity
  refine ⟨|C| * Real.exp y, τ, hτ0, by rw [hτ, one_div_lt_one_div (by positivity) hρ]; linarith,
    fun n => ?_⟩
  have hfn : (0:ℝ) < n.factorial := by positivity
  have h2 : ((n : ℝ) + 1) ≤ 2 ^ n := by
    have : n + 1 ≤ 2 ^ n := Nat.lt_two_pow_self
    exact_mod_cast this
  have hx := Real.pow_div_factorial_le_exp y hy n
  rw [abs_div, abs_of_pos hfn]
  calc |In A B n| / n.factorial ≤ |C| * ((n : ℝ) + 1) * X ^ n / n.factorial := by
        gcongr; exact (h n).trans (by gcongr; exact le_abs_self C)
    _ ≤ |C| * 2 ^ n * X ^ n / n.factorial := by gcongr
    _ = |C| * (y ^ n / n.factorial) * τ ^ n := by
        simp only [y]; rw [div_pow, mul_pow]; field_simp
    _ ≤ |C| * Real.exp y * τ ^ n := by gcongr

/-- `Q(η) = P(1/(1+η))`. -/
noncomputable def Qf (A B : ℝ) (η : ℝ) : ℝ := X01.P A B (1 / (1 + η))

lemma r_pos' {A B : ℝ} (hA : 0 < A) : 0 < X01.r A B := Real.sqrt_pos.mpr (by positivity)

lemma Qf_eq {A B : ℝ} (hA : 0 < A) {η : ℝ} (h : 1 + η ≠ 0) :
    Qf A B η = Ieta A B η * (4 * X01.r A B * (1 + η)) := by
  have := r_pos' (B := B) hA
  unfold Ieta Qf; field_simp

lemma Qf_ev {A B : ℝ} (hA : 0 < A) :
    Qf A B =ᶠ[𝓝 0] fun η => Ieta A B η * (4 * X01.r A B * (1 + η)) := by
  filter_upwards [Ioo_mem_nhds (by norm_num : (-1:ℝ) < 0) (by norm_num : (0:ℝ) < 1)] with η hη
  exact Qf_eq hA (by linarith [hη.1])

lemma Ieta_cd {A B : ℝ} (hA : 0 < A) (hB : 0 < B) : ContDiffAt ℝ (⊤ : ℕ∞) (Ieta A B) 0 :=
  (Y06.proof_3 A B hA hB).contDiffAt (Ioo_mem_nhds (by norm_num) (by norm_num))

lemma Qf_cd {A B : ℝ} (hA : 0 < A) (hB : 0 < B) : ContDiffAt ℝ (⊤ : ℕ∞) (Qf A B) 0 :=
  ((Ieta_cd hA hB).mul (by fun_prop)).congr_of_eventuallyEq (Qf_ev hA)

lemma P_eq_Qf (A B : ℝ) : X01.P A B = fun p => Qf A B (1 / p - 1) := by
  funext p; unfold Qf; congr 1; rw [add_sub_cancel, one_div_one_div]

lemma P_cd {A B : ℝ} (hA : 0 < A) (hB : 0 < B) : ContDiffAt ℝ (⊤ : ℕ∞) (X01.P A B) 1 := by
  have h1 : ContDiffAt ℝ (⊤ : ℕ∞) (fun p : ℝ => 1 / p - 1) 1 := by
    have : ContDiffAt ℝ (⊤ : ℕ∞) (fun p : ℝ => p⁻¹) 1 := contDiffAt_inv ℝ one_ne_zero
    simpa [one_div] using this.sub contDiffAt_const
  have h2 : ContDiffAt ℝ (⊤ : ℕ∞) (Qf A B) ((fun p : ℝ => 1 / p - 1) 1) := by
    norm_num; exact Qf_cd hA hB
  have h3 : ContDiffAt ℝ (⊤ : ℕ∞) (Qf A B ∘ (fun p : ℝ => 1 / p - 1)) 1 :=
    ContDiffAt.comp (g := Qf A B) (f := fun p : ℝ => 1 / p - 1) 1 h2 h1
  have e : X01.P A B = Qf A B ∘ (fun p : ℝ => 1 / p - 1) := P_eq_Qf A B
  rw [e]
  exact h3

lemma iter_affine (c : ℝ) (k : ℕ) :
    iteratedDeriv k (fun η : ℝ => c * (1 + η)) 0 = if k = 0 then c else if k = 1 then c else 0 := by
  rcases k with _ | k
  · simp
  · rw [iteratedDeriv_succ']
    have : deriv (fun η : ℝ => c * (1 + η)) = fun _ => c := by
      funext η
      have := ((hasDerivAt_id η).const_add 1).const_mul c
      simpa using this.deriv
    rw [this, iteratedDeriv_const]
    rcases k with _ | k <;> simp

lemma iter_Qf {A B : ℝ} (hA : 0 < A) (hB : 0 < B) (m : ℕ) :
    iteratedDeriv (m + 1) (Qf A B) 0 = 4 * X01.r A B *
      (iteratedDeriv (m + 1) (Ieta A B) 0 + (m + 1) * iteratedDeriv m (Ieta A B) 0) := by
  rw [(Qf_ev hA).iteratedDeriv_eq,
    iteratedDeriv_fun_mul ((Ieta_cd hA hB).of_le (mod_cast le_top)) (by fun_prop)]
  rw [Finset.sum_range_succ, Finset.sum_range_succ]
  have hz : ∑ i ∈ Finset.range m, ((m + 1).choose i : ℝ) * iteratedDeriv i (Ieta A B) 0 *
      iteratedDeriv (m + 1 - i) (fun η => 4 * X01.r A B * (1 + η)) 0 = 0 := by
    refine Finset.sum_eq_zero (fun i hi => ?_)
    have him := Finset.mem_range.mp hi
    rw [iter_affine, if_neg (by omega), if_neg (by omega), mul_zero]
  rw [hz, iter_affine, iter_affine, if_neg (by omega), if_pos (by omega), Nat.sub_self, if_pos rfl,
    Nat.choose_succ_self_right, Nat.choose_self]
  push_cast; ring

/-- Taylor coefficients of `Q` in terms of the `I_n`. -/
lemma qcoef {A B : ℝ} (hA : 0 < A) (hB : 0 < B) (m : ℕ) :
    iteratedDeriv (m + 1) (Qf A B) 0 / ((m + 1).factorial : ℝ) =
      4 * X01.r A B * (In A B (m + 1) + In A B m) := by
  rw [iter_Qf hA hB]
  unfold In
  rw [Nat.factorial_succ]; push_cast
  have : (m.factorial : ℝ) ≠ 0 := by positivity
  field_simp

lemma qcoef0 {A B : ℝ} (hA : 0 < A) (_hB : 0 < B) :
    iteratedDeriv 0 (Qf A B) 0 / ((0:ℕ).factorial : ℝ) = 4 * X01.r A B * In A B 0 := by
  simp only [iteratedDeriv_zero, Nat.factorial_zero, Nat.cast_one, div_one]
  unfold In; simp only [iteratedDeriv_zero, Nat.factorial_zero, Nat.cast_one, div_one]
  rw [Qf_eq hA (by norm_num)]; ring

/-- If the Taylor coefficients of `Q` are `O(X^n)`, contradiction. -/
lemma q_not_geom {A B : ℝ} (hA : 0 < A) (hB : 0 < B) {M X : ℝ} (hX : 1 ≤ X)
    (h : ∀ n : ℕ, |iteratedDeriv n (Qf A B) 0 / (n.factorial : ℝ)| ≤ M * X ^ n) : False := by
  have hr := r_pos' (B := B) hA
  have hc : 0 < 4 * X01.r A B := by positivity
  have hb : ∀ n : ℕ, |In A B n| ≤ M / (4 * X01.r A B) * ((n : ℝ) + 1) * X ^ n := by
    intro n
    induction n with
    | zero =>
      have := h 0
      rw [qcoef0 hA hB, abs_mul, abs_of_pos hc] at this
      simp only [pow_zero, mul_one] at this ⊢
      push_cast; rw [zero_add, mul_one, le_div_iff₀ hc]; linarith
    | succ m ih =>
      have h1 := h (m + 1)
      rw [qcoef hA hB, abs_mul, abs_of_pos hc] at h1
      have h2 : |In A B (m + 1)| ≤ |In A B (m + 1) + In A B m| + |In A B m| := by
        have := abs_sub (In A B (m + 1) + In A B m) (In A B m); simp at this; linarith
      have h3 : |In A B (m + 1) + In A B m| ≤ M / (4 * X01.r A B) * X ^ (m + 1) := by
        rw [div_mul_eq_mul_div, le_div_iff₀ hc]; linarith
      have hXm : X ^ m ≤ X ^ (m + 1) := pow_le_pow_right₀ hX (by omega)
      have hM : 0 ≤ M / (4 * X01.r A B) := by
        have := h 0; simp at this
        exact div_nonneg ((abs_nonneg _).trans this) hc.le
      push_cast
      calc |In A B (m + 1)| ≤ M / (4 * X01.r A B) * X ^ (m + 1) +
            M / (4 * X01.r A B) * ((m : ℝ) + 1) * X ^ m := by linarith
        _ ≤ M / (4 * X01.r A B) * X ^ (m + 1) +
            M / (4 * X01.r A B) * ((m : ℝ) + 1) * X ^ (m + 1) := by gcongr
        _ = _ := by ring
  exact no_geom_In hA hB (by linarith) hb

lemma q_not_summable {A B : ℝ} (hA : 0 < A) (hB : 0 < B) {η : ℝ} (hη : η ≠ 0) :
    ¬ Summable (fun n : ℕ => iteratedDeriv n (Qf A B) 0 / (n.factorial : ℝ) * η ^ n) := by
  intro hs
  obtain ⟨M, hM⟩ := (isBounded_iff_forall_norm_le.mp
    (Metric.isBounded_range_of_tendsto _ hs.tendsto_atTop_zero))
  have hη0 : 0 < |η| := abs_pos.mpr hη
  set X := max (|η|⁻¹) 1
  refine q_not_geom hA hB (le_max_right _ 1) (M := M) (X := X) (fun n => ?_)
  have := hM _ ⟨n, rfl⟩
  rw [Real.norm_eq_abs, abs_mul, abs_pow] at this
  have hp : 0 < |η| ^ n := by positivity
  have h1 : |iteratedDeriv n (Qf A B) 0 / (n.factorial : ℝ)| ≤ M * |η|⁻¹ ^ n := by
    rw [inv_pow, le_mul_inv_iff₀ hp]; exact this
  have hM0 : 0 ≤ M := (mul_nonneg (abs_nonneg _) hp.le).trans this
  exact h1.trans (by gcongr; exact le_max_left _ _)

end Y06

namespace Y06
open MeasureTheory Set Filter Topology

lemma hasDerivAt_psi {η : ℝ} (h : 1 + η ≠ 0) :
    HasDerivAt (fun η : ℝ => 1 / (1 + η)) (-(1 / (1 + η)) ^ 2) η := by
  have h2 : HasDerivAt (fun η : ℝ => 1 + η) 1 η := by simpa using (hasDerivAt_id η).const_add 1
  have h3 := h2.inv h
  simp only [one_div]
  convert h3 using 1
  rw [inv_pow]; ring

/-- A function flat at `1` stays flat at `0` after composing with `η ↦ 1/(1+η)`. -/
lemma flat_comp : ∀ n : ℕ, ∀ {h : ℝ → ℝ} {V : Set ℝ}, IsOpen V → (1:ℝ) ∈ V →
    ContDiffOn ℝ (⊤ : ℕ∞) h V → (∀ k : ℕ, iteratedDeriv k h 1 = 0) →
    iteratedDeriv n (fun η : ℝ => h (1 / (1 + η))) 0 = 0 := by
  intro n
  induction n with
  | zero =>
    intro h V _ _ _ hflat
    simpa using hflat 0
  | succ n ih =>
    intro h V hV h1 hcd hflat
    rw [iteratedDeriv_succ']
    set Dh : ℝ → ℝ := fun p => -(p ^ 2) * deriv h p with hDh
    have hψc : ContinuousAt (fun η : ℝ => 1 / (1 + η)) 0 :=
      (hasDerivAt_psi (by norm_num : (1:ℝ) + 0 ≠ 0)).continuousAt
    have hmem : ∀ᶠ η in 𝓝 (0:ℝ), 1 / (1 + η) ∈ V :=
      hψc.preimage_mem_nhds (by simpa using hV.mem_nhds h1)
    have hne : ∀ᶠ η in 𝓝 (0:ℝ), 1 + η ≠ 0 := by
      filter_upwards [Ioo_mem_nhds (show (-1:ℝ) < 0 by norm_num) (show (0:ℝ) < 1 by norm_num)]
        with η hη
      exact ne_of_gt (by linarith [hη.1])
    have hev : deriv (fun η : ℝ => h (1 / (1 + η))) =ᶠ[𝓝 0] fun η => Dh (1 / (1 + η)) := by
      filter_upwards [hmem, hne] with η hηV hη
      have hd : DifferentiableAt ℝ h (1 / (1 + η)) :=
        ((hcd _ hηV).contDiffAt (hV.mem_nhds hηV)).differentiableAt (by simp)
      have := hd.hasDerivAt.comp η (hasDerivAt_psi hη)
      rw [show (fun η : ℝ => h (1 / (1 + η))) = h ∘ (fun η : ℝ => 1 / (1 + η)) from rfl,
        this.deriv, hDh]; ring
    rw [hev.iteratedDeriv_eq]
    have hcdD : ContDiffOn ℝ (⊤ : ℕ∞) Dh V := by
      have := hcd.deriv_of_isOpen hV (m := (⊤ : ℕ∞)) (by simp)
      exact (contDiffOn_id.pow 2).neg.mul this
    refine ih hV h1 hcdD (fun k => ?_)
    have hc1 : ContDiffAt ℝ k (fun p : ℝ => -(p ^ 2)) 1 := by fun_prop
    have hc2 : ContDiffAt ℝ k (deriv h) 1 :=
      ((hcd.deriv_of_isOpen hV (m := (⊤ : ℕ∞)) (by simp)).contDiffAt (hV.mem_nhds h1)).of_le
        (mod_cast le_top)
    rw [hDh, iteratedDeriv_fun_mul hc1 hc2]
    refine Finset.sum_eq_zero (fun i _ => ?_)
    rw [← iteratedDeriv_succ', hflat, mul_zero]

lemma geom_le {C s : ℝ} (hs : 0 < s) (n : ℕ) : C * s ^ n ≤ |C| * (max s 1) ^ n :=
  mul_le_mul (le_abs_self C) (pow_le_pow_left₀ hs.le (le_max_left _ _) n) (by positivity)
    (abs_nonneg C)

lemma psi_analytic : AnalyticAt ℝ (fun η : ℝ => 1 / (1 + η)) 0 :=
  analyticAt_const.div (analyticAt_const.add analyticAt_id) (show (1:ℝ) + (0:ℝ) ≠ 0 by norm_num)

lemma Qf_on {A B : ℝ} (hA : 0 < A) (hB : 0 < B) :
    ContDiffOn ℝ (⊤ : ℕ∞) (Qf A B) (Ioo (-1) 1) :=
  ((Y06.proof_3 A B hA hB).mul (by fun_prop)).congr
    (fun η hη => Qf_eq hA (ne_of_gt (by linarith [hη.1])))

lemma P_on {A B : ℝ} (hA : 0 < A) (hB : 0 < B) :
    ContDiffOn ℝ (⊤ : ℕ∞) (X01.P A B) (Ioo (1 / 2) (3 / 2)) := by
  have h1 : ContDiffOn ℝ (⊤ : ℕ∞) (fun p : ℝ => 1 / p - 1) (Ioo (1 / 2) (3 / 2)) := by
    intro p hp
    have hp0 : p ≠ 0 := ne_of_gt (by linarith [hp.1])
    have : ContDiffAt ℝ (⊤ : ℕ∞) (fun p : ℝ => p⁻¹) p := contDiffAt_inv ℝ hp0
    exact (by simpa [one_div] using this.sub contDiffAt_const :
      ContDiffAt ℝ (⊤ : ℕ∞) (fun p : ℝ => 1 / p - 1) p).contDiffWithinAt
  have hmap : MapsTo (fun p : ℝ => 1 / p - 1) (Ioo (1 / 2) (3 / 2)) (Ioo (-1) 1) := by
    intro p hp
    have hp0 : 0 < p := by linarith [hp.1]
    constructor
    · have : 0 < 1 / p := by positivity
      linarith
    · have : 1 / p < 2 := by rw [div_lt_iff₀ hp0]; linarith [hp.1]
      linarith
  have h3 := ContDiffOn.comp (g := Qf A B) (f := fun p : ℝ => 1 / p - 1) (Qf_on hA hB) h1 hmap
  have e : X01.P A B = Qf A B ∘ (fun p : ℝ => 1 / p - 1) := P_eq_Qf A B
  rw [e]
  exact h3

/-- `P` is not real analytic at `p = 1`. -/
lemma P_not_analytic {A B : ℝ} (hA : 0 < A) (hB : 0 < B) : ¬ AnalyticAt ℝ (X01.P A B) 1 := by
  intro hP
  have h1 : AnalyticAt ℝ (X01.P A B) ((fun η : ℝ => 1 / (1 + η)) 0) := by norm_num; exact hP
  have hQ : AnalyticAt ℝ (Qf A B) 0 :=
    AnalyticAt.comp (g := X01.P A B) (f := fun η : ℝ => 1 / (1 + η)) h1 psi_analytic
  obtain ⟨C, s, hs, hC⟩ := analytic_coeff hQ
  exact q_not_geom hA hB (le_max_right s 1) (M := |C|) (fun n => (hC n).trans (geom_le hs n))

/-- The Taylor series of `P` at `p = 1` (in `d`) diverges at every `d ≠ 0`. -/
lemma P_not_summable {A B : ℝ} (hA : 0 < A) (hB : 0 < B) {d : ℝ} (hd : d ≠ 0) :
    ¬ Summable (fun n : ℕ => iteratedDeriv n (X01.P A B) 1 / (n.factorial : ℝ) * d ^ n) := by
  intro hs
  obtain ⟨M, hM⟩ := (isBounded_iff_forall_norm_le.mp
    (Metric.isBounded_range_of_tendsto _ hs.tendsto_atTop_zero))
  have hd0 : 0 < |d| := abs_pos.mpr hd
  set pc : ℕ → ℝ := fun n => iteratedDeriv n (X01.P A B) 1 / (n.factorial : ℝ) with hpc
  have hb : ∀ n, |pc n| * |d| ^ n ≤ M := fun n => by
    have := hM _ ⟨n, rfl⟩; rwa [Real.norm_eq_abs, abs_mul, abs_pow] at this
  have hM0 : 0 ≤ M := (mul_nonneg (abs_nonneg _) (by positivity)).trans (hb 0)
  set ps := FormalMultilinearSeries.ofScalars ℝ pc
  set r : ℝ := |d| / 2 with hr
  have hr0 : 0 < r := by positivity
  have hrad : ((r.toNNReal : NNReal) : ENNReal) ≤ ps.radius := by
    refine ps.le_radius_of_bound M (fun n => ?_)
    rw [FormalMultilinearSeries.ofScalars_norm, Real.norm_eq_abs, Real.coe_toNNReal r hr0.le]
    have h1 : r ^ n ≤ |d| ^ n := pow_le_pow_left₀ hr0.le (by rw [hr]; linarith) n
    calc |pc n| * r ^ n ≤ |pc n| * |d| ^ n := by gcongr
      _ ≤ M := hb n
  have hpos : 0 < ps.radius := lt_of_lt_of_le (by simpa using hr0) hrad
  have hps0 := ps.hasFPowerSeriesOnBall hpos
  have hps1 := hps0.comp_sub 1
  rw [zero_add] at hps1
  set Ph : ℝ → ℝ := fun z => ps.sum (z - 1) with hPh
  have hiter : ∀ n, iteratedDeriv n Ph 1 = iteratedDeriv n (X01.P A B) 1 := by
    intro n
    rw [iter_of_hasFPS hps1 n, FormalMultilinearSeries.coeff_ofScalars, hpc]
    simp only
    have : (n.factorial : ℝ) ≠ 0 := by positivity
    field_simp
  -- `Ph` is analytic on the ball of radius `r` about `1`
  have hPh_an : ∀ z ∈ Ioo (1 - r) (1 + r), AnalyticAt ℝ Ph z := by
    intro z hz
    apply hps1.analyticAt_of_mem
    have : z ∈ Metric.eball (1:ℝ) ((r.toNNReal : NNReal) : ENNReal) := by
      rw [Metric.eball_coe, Real.coe_toNNReal r hr0.le, Metric.mem_ball, Real.dist_eq, abs_lt]
      constructor <;> linarith [hz.1, hz.2]
    exact Metric.eball_subset_eball hrad this
  set V : Set ℝ := Ioo (max (1 / 2) (1 - r)) (min (3 / 2) (1 + r)) with hV
  have hVo : IsOpen V := isOpen_Ioo
  have h1V : (1:ℝ) ∈ V := ⟨max_lt (by norm_num) (by linarith), lt_min (by norm_num) (by linarith)⟩
  have hVP : V ⊆ Ioo (1 / 2) (3 / 2) := fun z hz =>
    ⟨lt_of_le_of_lt (le_max_left _ _) hz.1, lt_of_lt_of_le hz.2 (min_le_left _ _)⟩
  have hVr : V ⊆ Ioo (1 - r) (1 + r) := fun z hz =>
    ⟨lt_of_le_of_lt (le_max_right _ _) hz.1, lt_of_lt_of_le hz.2 (min_le_right _ _)⟩
  set hh : ℝ → ℝ := fun p => X01.P A B p - Ph p with hhh
  have hcd : ContDiffOn ℝ (⊤ : ℕ∞) hh V :=
    ((P_on hA hB).mono hVP).sub (fun z hz => (hPh_an z (hVr hz)).contDiffAt.contDiffWithinAt)
  have hflat : ∀ k : ℕ, iteratedDeriv k hh 1 = 0 := by
    intro k
    rw [hhh, iteratedDeriv_fun_sub ((P_cd hA hB).of_le (mod_cast le_top))
      ((hPh_an 1 ⟨by linarith, by linarith⟩).contDiffAt), hiter, sub_self]
  have hsplit : Qf A B = fun η => hh (1 / (1 + η)) + Ph (1 / (1 + η)) := by
    funext η; simp only [hhh, Qf]; ring
  have hcd1 : ContDiffAt ℝ (⊤ : ℕ∞) hh ((fun η : ℝ => 1 / (1 + η)) 0) := by
    norm_num; exact (hcd 1 h1V).contDiffAt (hVo.mem_nhds h1V)
  have hcomp1 : ContDiffAt ℝ (⊤ : ℕ∞) (fun η : ℝ => hh (1 / (1 + η))) 0 :=
    ContDiffAt.comp (g := hh) (f := fun η : ℝ => 1 / (1 + η)) 0 hcd1 psi_analytic.contDiffAt
  have hPh1 : AnalyticAt ℝ Ph ((fun η : ℝ => 1 / (1 + η)) 0) := by
    norm_num; exact hPh_an 1 ⟨by linarith, by linarith⟩
  have hcomp2 : AnalyticAt ℝ (fun η : ℝ => Ph (1 / (1 + η))) 0 :=
    AnalyticAt.comp (g := Ph) (f := fun η : ℝ => 1 / (1 + η)) hPh1 psi_analytic
  obtain ⟨C, s, hs, hC⟩ := analytic_coeff hcomp2
  refine q_not_geom hA hB (le_max_right s 1) (M := |C|) (fun n => ?_)
  rw [hsplit, iteratedDeriv_fun_add (hcomp1.of_le (mod_cast le_top)) hcomp2.contDiffAt,
    flat_comp n hVo h1V hcd hflat, zero_add]
  exact (hC n).trans (geom_le hs n)

end Y06

theorem Y06.proof_11 : Y06.stmt_11 := fun A B hA hB =>
  ⟨Y06.proof_9 A B hA hB, Y06.Ieta_cd hA hB, Y06.proof_10 A B hA hB, Y06.Qf_cd hA hB,
    fun _ hη => Y06.q_not_summable hA hB hη, Y06.P_cd hA hB,
    fun _ hd => Y06.P_not_summable hA hB hd, Y06.P_not_analytic hA hB⟩
