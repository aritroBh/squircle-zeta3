import Mathlib
/-!
# X01 -- derivatives of the rhombus perimeter at `p = 1`

Target id: X01.
Paper labels: `thm:rhombus` and its proof (paper/ieee/main.tex, lines 1136-1149), the derivation
"Steps 1-4" it reuses (main.tex lines 1101-1134), and `prop:an` for `a_0,...,a_3` only
(main.tex lines 1163-1164). Supplement: paper/ieee/proofs/sec_pureweight.tex, Section
"Rhombi (two dimensions)" (`eq:pw:rhombus`, "The Beta moments"), lines 54-80.

Quoted verbatim (main.tex 1102-1106, Step 1):

  \emph{Step 1: an exact formula.} Put $\eta=\frac1p-1$, so the diamond is $\eta=0$, and use the
  parametrization of Section~\ref{sec:squircle}. The speed is $(1+\eta)\sqrt{t^{2\eta}+(1-t)^{2\eta}}$.
  Let's factor out $[t(1-t)]^{\eta/2}$ and use the logit coordinate $l=\log\frac{t}{1-t}$. ...
  P(p)=4\sqrt2\,(1+\eta)\int_0^1\bigl(t(1-t)\bigr)^{\eta/2}\sqrt{\cosh(\eta l)}\;dt
  exactly, for every $p>0$. At $\eta=0$ the integrand is $1$ and $P(1)=4\sqrt2$.

Quoted verbatim (main.tex 1108-1112, Step 2):

  \emph{Step 2: Beta moments.} Every moment of $l$ comes from one Beta integral, because
  $e^{sl}=(t/(1-t))^s$:
  A(\eta,s)&=(1+\eta)\int_0^1\bigl(t(1-t)\bigr)^{\eta/2}e^{sl}\dd t\\
  &=\frac{\Gamma\bigl(1+\frac\eta2+s\bigr)\,\Gamma\bigl(1+\frac\eta2-s\bigr)}{\Gamma(1+\eta)} .

Quoted verbatim (main.tex 1130-1132, Step 4):

  \emph{Step 4: collecting terms.} From $A(\eta,0)=1-\frac{\pi^2}{24}\eta^2+\frac{\zeta(3)}{4}\eta^3+O(\eta^4)$
  and $\partial_s^2A(\eta,0)=\frac{\pi^2}{3}-2\zeta(3)\eta+O(\eta^2)$,
  \frac{P(p)}{4\sqrt2}=1+\frac{\pi^2}{24}\eta^2-\frac{\zeta(3)}{4}\eta^3+O(\eta^4).

Quoted verbatim (main.tex 1137-1148):

  The same four steps work for $|x/A|^p+|y/B|^p=1$ with semi-axes $A,B>0$. Put $r=\sqrt{A^2+B^2}$
  and $q=(A^2-B^2)/(A^2+B^2)$. Then $\sqrt{\cosh z}$ becomes
  $g(z)=\sqrt{\cosh z+q\sinh z}=1+\frac q2z+\frac{2-q^2}{8}z^2+\cdots$; its odd terms meet odd
  moments of $l$, which vanish by the symmetry $t\leftrightarrow1-t$.

  \begin{theorem}[Derivatives at the rhombus]\label{thm:rhombus}
  For $|x/A|^p+|y/B|^p=1$ with $r=\sqrt{A^2+B^2}$,
  \begin{gather*}
  P(1)=4r,\qquad P'(1)=0,\qquad P''(1)=\frac{4\pi^2A^2B^2}{3r^3},\\
  P'''(1)=\frac{24A^2B^2}{r^3}\Bigl(\zeta(3)-\frac{\pi^2}{3}\Bigr).
  \end{gather*}
  \end{theorem}

  \begin{proof}
  As above, $P/(4r)=1+\alpha\eta^2+\beta\eta^3+O(\eta^4)$ with $\alpha=\frac{(1-q^2)\pi^2}{24}$,
  $\beta=-\frac{(1-q^2)\zeta(3)}{4}$ and $1-q^2=4A^2B^2/r^4$. ... Differentiating under the
  integral is justified by an integrable bound (powers of $t$, $1-t$ with exponents $>-1$ times
  powers of logarithms).
  \end{proof}

Quoted verbatim (main.tex 1163-1164):

  \begin{proposition}\label{prop:an}
  Let $a_n=P^{(n)}(1)/(\sqrt2\,n!)$, the Taylor coefficients of $P(1+d)/\sqrt2$ at $d=0$ (the
  series diverges, by Theorem~\ref{thm:zero}). Then $a_0=4$, $a_1=0$, $a_2=\pi^2/6$,
  $a_3=\zeta(3)-\pi^2/3$,

## Formalization choices

* `ζ(3)` is `z3 := (riemannZeta 3).re`, the real part of Mathlib's complex `riemannZeta 3`
  (it is real; `zeta_nat_eq_tsum_of_gt_one` identifies `riemannZeta 3` with `∑ 1/n^3`).
* The perimeter function. The paper writes the perimeter integral (`eq:Pint`, main.tex 827)
  only for `A = B = 1` and says "the same four steps work" for semi-axes `A, B`. We use the
  natural integral from the parametrization `x = A t^{1/p}`, `y = B (1-t)^{1/p}` (as the
  supplement does, `eq:pw:rhombus`):
    `P A B p = (4/p) ∫_0^1 √(A^2 t^(2/p-2) + B^2 (1-t)^(2/p-2)) dt`.
  This makes sense for EVERY `p > 0` (the integrand is `O(t^{1/p-1})` at `0` and
  `O((1-t)^{1/p-1})` at `1`, integrable since `1/p - 1 > -1`); for `p < 1` the integrand is even
  continuous on `[0,1]`. So `P A B` is defined on the neighbourhood `(0,∞)` of `1`, as needed
  for derivatives at `1` (stmt_1). `Real.rpow` has `0^0 = 1` and `0^x = 0` for `x ≠ 0`; this
  only affects the null set `{0,1}`.
* True length (stmt_2, stmt_3): `loop A B p` on `[0,4]` is the quarter arc followed by its
  reflections (as `W04.loop`, scaled by `A`, `B`), valued in the Euclidean plane; the perimeter
  is its total variation `eVariationOn`. Stated for all `p > 0`.
* Derivatives are genuine: `ContDiffAt ℝ 3 (P A B) 1` is asserted (it gives an open
  neighbourhood of `1` on which `P A B` is `C^3`), and the values are `deriv` / `iteratedDeriv`
  at `1`.
* Log integrals (the intermediate steps that produce `π^2` and `ζ(3)`): with
  `L t = log(t(1-t))` and the paper's logit `l t = log(t/(1-t))`, differentiating
  `P/(4r) = (1+η) ∫ (t(1-t))^{η/2} g(η l) dt` three times at `η = 0` needs exactly
  `∫L = -2`, `∫L^2 = 8 - π^2/3`, `∫L^3 = 12ζ(3) - 48 + 2π^2`, `∫l^2 = π^2/3`,
  `∫L l^2 = -4ζ(3) - 2π^2/3`, and the vanishing of the odd moments of `l`. These are the
  derivatives at `η = 0` of the paper's `A(η,0)` and `∂_s^2 A(η,0)` (Step 4):
  `A(η,0) = (1+η)∫(t(1-t))^{η/2}`, `∂_s^2A(η,0) = (1+η)∫(t(1-t))^{η/2} l^2`.
* The paper's `O(η^4)` in the proof is recorded only through derivatives up to order 3 of
  `η ↦ P(1/(1+η))/(4r)` at `η = 0` (stmt_10), which is what the theorem uses.

Statements:
* stmt_1: the integrand is interval integrable on `[0,1]` for every `p > 0` (so `P` is a genuine
  integral near `p = 1`, including `p < 1`).
* stmt_2: `loop A B p` is a simple closed parametrization of `|x/A|^p+|y/B|^p=1`, `p > 0`.
* stmt_3: the perimeter (length of `loop A B p`) equals `P A B p`, `p > 0`.
* stmt_4: Step 1, exact formula `P = 4r(1+η)∫(t(1-t))^{η/2} √(cosh(ηl)+q sinh(ηl)) dt`, `p > 0`.
* stmt_5: Step 2, the Beta moments `A(η,s) = Γ(1+η/2+s)Γ(1+η/2-s)/Γ(1+η)` (with integrability),
  for `η > -1`, `|s| < 1 + η/2`.
* stmt_6: all mixed log moments `∫ L^j l^k` exist, and vanish for odd `k` (symmetry `t ↔ 1-t`).
* stmt_7: `∫L = -2`, `∫L^2 = 8 - π^2/3`, `∫L^3 = 12ζ(3) - 48 + 2π^2`.
* stmt_8: `∫ l^2 = π^2/3` (the source of `π^2`).
* stmt_9: `∫ L l^2 = -4ζ(3) - 2π^2/3` (the source of `ζ(3)`).
* stmt_10: in the variable `η`, `I(η) = P(1/(1+η))/(4r)` is `C^3` at `0` with `I(0)=1`, `I'(0)=0`,
  `I''(0) = 2α`, `I'''(0) = 6β`, `α = (1-q^2)π^2/24`, `β = -(1-q^2)ζ(3)/4`.
* stmt_11: `P A B` is `C^3` near `p = 1`.
* stmt_12: `P(1) = 4r`.
* stmt_13: `P'(1) = 0`.
* stmt_14: `P''(1) = 4π^2A^2B^2/(3r^3)`.
* stmt_15: `P'''(1) = (24A^2B^2/r^3)(ζ(3) - π^2/3)`.
* stmt_16: `prop:an` for the diamond `A = B = 1`: `a_0 = 4`, `a_1 = 0`, `a_2 = π^2/6`,
  `a_3 = ζ(3) - π^2/3`, with `a_n = P^{(n)}(1)/(√2 n!)`.
* stmt_17: `thm:rhombus` in full: for all `A, B > 0`, `P` is `C^3` at `1` and the four values hold.

## NOT COVERED

* The boxed expansion `eq:diamond` with its `O(d^4)` remainder and the `O(η^4)` remainders of
  Step 4 / the proof (only derivatives up to order 3 are stated; the fourth-order remainder
  needs `C^4`, which is beyond `thm:rhombus`).
* Step 3 (the Taylor series of `log Γ(1+x)` and of `log A(η,s)`); the log-integral route of
  stmt_7-stmt_9 replaces it.
* `prop:an` beyond `a_3` (`a_4,...,a_7`, `ω_n`, weights) and the divergence remark.
* The literature remarks (what was known, Rodríguez-Montilla, Abel summation).
* "Perimeter" is the total variation of the explicit loop; no parametrization-free notion
  (Hausdorff measure of the curve) is used.
-/
-- BEGIN STATEMENTS
namespace X01

open MeasureTheory Set Filter

/-- The Euclidean plane. -/
abbrev Plane := EuclideanSpace ℝ (Fin 2)

/-- `ζ(3)`, the real part of Mathlib's `riemannZeta 3`. -/
noncomputable def z3 : ℝ := (riemannZeta 3).re

/-- `r = √(A^2 + B^2)`. -/
noncomputable def r (A B : ℝ) : ℝ := Real.sqrt (A ^ 2 + B ^ 2)

/-- `q = (A^2 - B^2)/(A^2 + B^2)`. -/
noncomputable def q (A B : ℝ) : ℝ := (A ^ 2 - B ^ 2) / (A ^ 2 + B ^ 2)

/-- The integrand `√(A^2 t^(2/p-2) + B^2 (1-t)^(2/p-2))` of the perimeter integral. -/
noncomputable def integrand (A B p t : ℝ) : ℝ :=
  Real.sqrt (A ^ 2 * t ^ (2 / p - 2) + B ^ 2 * (1 - t) ^ (2 / p - 2))

/-- The perimeter integral `P(p) = (4/p) ∫_0^1 √(A^2 t^(2/p-2) + B^2 (1-t)^(2/p-2)) dt`
(the paper's `eq:Pint` with semi-axes `A, B`; see the module docstring). -/
noncomputable def P (A B p : ℝ) : ℝ := 4 / p * ∫ t in (0 : ℝ)..1, integrand A B p t

/-- The curve `|x/A|^p + |y/B|^p = 1`. -/
def curve (A B p : ℝ) : Set Plane := {v | |v 0 / A| ^ p + |v 1 / B| ^ p = 1}

/-- The whole curve traversed once on `[0,4]`: the quarter arc `x = A t^{1/p}`,
`y = B (1-t)^{1/p}` on `[0,1]`, then its reflections into the other three quadrants. -/
noncomputable def loop (A B p t : ℝ) : Plane :=
  if t ≤ 1 then !₂[A * t ^ (1 / p), B * (1 - t) ^ (1 / p)]
  else if t ≤ 2 then !₂[A * (2 - t) ^ (1 / p), -(B * (t - 1) ^ (1 / p))]
  else if t ≤ 3 then !₂[-(A * (t - 2) ^ (1 / p)), -(B * (3 - t) ^ (1 / p))]
  else !₂[-(A * (4 - t) ^ (1 / p)), B * (t - 3) ^ (1 / p)]

/-- The perimeter: the Euclidean length (total variation) of `loop A B p` over `[0,4]`. -/
noncomputable def perimeter (A B p : ℝ) : ENNReal := eVariationOn (loop A B p) (Icc (0 : ℝ) 4)

/-- `L(t) = log(t(1-t))`. -/
noncomputable def Lg (t : ℝ) : ℝ := Real.log (t * (1 - t))

/-- The paper's logit coordinate `l = log(t/(1-t))`. -/
noncomputable def lg (t : ℝ) : ℝ := Real.log (t / (1 - t))

/-- The paper's `P/(4r)` in the variable `η = 1/p - 1`, i.e. `p = 1/(1+η)`. -/
noncomputable def I (A B η : ℝ) : ℝ := P A B (1 / (1 + η)) / (4 * r A B)

/-- "exactly, for every $p>0$": the perimeter integral makes sense for every `p > 0`, in
particular on a neighbourhood of `p = 1` (its integrand is integrable). -/
def stmt_1 : Prop :=
  ∀ A B : ℝ, 0 < A → 0 < B → ∀ p : ℝ, 0 < p →
    IntervalIntegrable (integrand A B p) volume 0 1

/-- "For $|x/A|^p+|y/B|^p=1$": the quarter arc `x = A t^{1/p}`, `y = B(1-t)^{1/p}` and its three
reflections form a simple closed parametrization of the curve. -/
def stmt_2 : Prop :=
  ∀ A B : ℝ, 0 < A → 0 < B → ∀ p : ℝ, 0 < p →
    ContinuousOn (loop A B p) (Icc 0 4) ∧ loop A B p 4 = loop A B p 0 ∧
    InjOn (loop A B p) (Ico 0 4) ∧ loop A B p '' Icc 0 4 = curve A B p

/-- The perimeter integral is the true length of the curve `|x/A|^p+|y/B|^p=1`. -/
def stmt_3 : Prop :=
  ∀ A B : ℝ, 0 < A → 0 < B → ∀ p : ℝ, 0 < p →
    perimeter A B p = ENNReal.ofReal (P A B p)

/-- Step 1 for the rhombus (`eq:pw:rhombus`): with `η = 1/p - 1`,
`P = 4r(1+η) ∫_0^1 (t(1-t))^{η/2} √(cosh(ηl) + q sinh(ηl)) dt`, "exactly, for every $p>0$". -/
def stmt_4 : Prop :=
  ∀ A B : ℝ, 0 < A → 0 < B → ∀ p : ℝ, 0 < p →
    P A B p = 4 * r A B * (1 + (1 / p - 1)) *
      ∫ t in (0 : ℝ)..1, (t * (1 - t)) ^ ((1 / p - 1) / 2) *
        Real.sqrt (Real.cosh ((1 / p - 1) * lg t) + q A B * Real.sinh ((1 / p - 1) * lg t))

/-- Step 2, the Beta moments:
`A(η,s) = (1+η)∫_0^1 (t(1-t))^{η/2} e^{sl} dt = Γ(1+η/2+s)Γ(1+η/2-s)/Γ(1+η)`,
for `η > -1` (i.e. `p > 0`) and `|s| < 1 + η/2`, with the integrand integrable. -/
def stmt_5 : Prop :=
  ∀ η s : ℝ, -1 < η → |s| < 1 + η / 2 →
    IntervalIntegrable (fun t : ℝ => (t * (1 - t)) ^ (η / 2) * Real.exp (s * lg t)) volume 0 1 ∧
    (1 + η) * ∫ t in (0 : ℝ)..1, (t * (1 - t)) ^ (η / 2) * Real.exp (s * lg t) =
      Real.Gamma (1 + η / 2 + s) * Real.Gamma (1 + η / 2 - s) / Real.Gamma (1 + η)

/-- "odd moments of $l$, which vanish by the symmetry $t\leftrightarrow1-t$": every mixed log
moment `∫ L^j l^k` exists, and it vanishes for odd `k`. -/
def stmt_6 : Prop :=
  ∀ j k : ℕ, IntervalIntegrable (fun t : ℝ => Lg t ^ j * lg t ^ k) volume 0 1 ∧
    (Odd k → ∫ t in (0 : ℝ)..1, Lg t ^ j * lg t ^ k = 0)

/-- Step 4, the `η`-derivatives of `A(η,0) = 1 - π^2/24 η^2 + ζ(3)/4 η^3 + O(η^4)` written as
log integrals: `∫L = -2`, `∫L^2 = 8 - π^2/3`, `∫L^3 = 12ζ(3) - 48 + 2π^2`. -/
def stmt_7 : Prop :=
  ∫ t in (0 : ℝ)..1, Lg t = -2 ∧
  ∫ t in (0 : ℝ)..1, Lg t ^ 2 = 8 - Real.pi ^ 2 / 3 ∧
  ∫ t in (0 : ℝ)..1, Lg t ^ 3 = 12 * z3 - 48 + 2 * Real.pi ^ 2

/-- Step 4, `∂_s^2 A(0,0) = π^2/3`: `∫_0^1 log^2(t/(1-t)) dt = π^2/3` (where `π^2` enters). -/
def stmt_8 : Prop :=
  ∫ t in (0 : ℝ)..1, lg t ^ 2 = Real.pi ^ 2 / 3

/-- Step 4, `∂_s^2 A(η,0) = π^2/3 - 2ζ(3)η + O(η^2)` at first order in `η`:
`∫_0^1 log(t(1-t)) log^2(t/(1-t)) dt = -4ζ(3) - 2π^2/3` (where `ζ(3)` enters). -/
def stmt_9 : Prop :=
  ∫ t in (0 : ℝ)..1, Lg t * lg t ^ 2 = -4 * z3 - 2 * Real.pi ^ 2 / 3

/-- Proof of `thm:rhombus`: "$P/(4r)=1+\alpha\eta^2+\beta\eta^3+O(\eta^4)$ with
$\alpha=\frac{(1-q^2)\pi^2}{24}$, $\beta=-\frac{(1-q^2)\zeta(3)}{4}$", recorded through the
derivatives of order `≤ 3` of `η ↦ P(1/(1+η))/(4r)` at `η = 0`. -/
def stmt_10 : Prop :=
  ∀ A B : ℝ, 0 < A → 0 < B →
    ContDiffAt ℝ 3 (I A B) 0 ∧ I A B 0 = 1 ∧ deriv (I A B) 0 = 0 ∧
    iteratedDeriv 2 (I A B) 0 = 2 * ((1 - q A B ^ 2) * Real.pi ^ 2 / 24) ∧
    iteratedDeriv 3 (I A B) 0 = 6 * (-((1 - q A B ^ 2) * z3 / 4))

/-- "Differentiating under the integral is justified": `p ↦ P(p)` is `C^3` near `p = 1`. -/
def stmt_11 : Prop :=
  ∀ A B : ℝ, 0 < A → 0 < B → ContDiffAt ℝ 3 (P A B) 1

/-- `thm:rhombus`: `P(1) = 4r`. -/
def stmt_12 : Prop :=
  ∀ A B : ℝ, 0 < A → 0 < B → P A B 1 = 4 * r A B

/-- `thm:rhombus`: `P'(1) = 0`. -/
def stmt_13 : Prop :=
  ∀ A B : ℝ, 0 < A → 0 < B → HasDerivAt (P A B) 0 1

/-- `thm:rhombus`: `P''(1) = 4π^2A^2B^2/(3r^3)`. -/
def stmt_14 : Prop :=
  ∀ A B : ℝ, 0 < A → 0 < B →
    iteratedDeriv 2 (P A B) 1 = 4 * Real.pi ^ 2 * A ^ 2 * B ^ 2 / (3 * r A B ^ 3)

/-- `thm:rhombus`: `P'''(1) = (24A^2B^2/r^3)(ζ(3) - π^2/3)`. -/
def stmt_15 : Prop :=
  ∀ A B : ℝ, 0 < A → 0 < B →
    iteratedDeriv 3 (P A B) 1 = 24 * A ^ 2 * B ^ 2 / r A B ^ 3 * (z3 - Real.pi ^ 2 / 3)

/-- `prop:an` (first four coefficients), the diamond `A = B = 1`: with
`a_n = P^{(n)}(1)/(√2 n!)`, `a_0 = 4`, `a_1 = 0`, `a_2 = π^2/6`, `a_3 = ζ(3) - π^2/3`. -/
def stmt_16 : Prop :=
  ContDiffAt ℝ 3 (P 1 1) 1 ∧
  (let a : ℕ → ℝ := fun n => iteratedDeriv n (P 1 1) 1 / (Real.sqrt 2 * (n.factorial : ℝ))
   a 0 = 4 ∧ a 1 = 0 ∧ a 2 = Real.pi ^ 2 / 6 ∧ a 3 = z3 - Real.pi ^ 2 / 3)

/-- `thm:rhombus` in full: for every rhombus `|x/A|^p+|y/B|^p=1`, `P` is three times
continuously differentiable at `p = 1` and `P(1)=4r`, `P'(1)=0`,
`P''(1)=\frac{4\pi^2A^2B^2}{3r^3}`, `P'''(1)=\frac{24A^2B^2}{r^3}(\zeta(3)-\frac{\pi^2}{3})`. -/
def stmt_17 : Prop :=
  ∀ A B : ℝ, 0 < A → 0 < B →
    ContDiffAt ℝ 3 (P A B) 1 ∧
    P A B 1 = 4 * r A B ∧
    deriv (P A B) 1 = 0 ∧
    iteratedDeriv 2 (P A B) 1 = 4 * Real.pi ^ 2 * A ^ 2 * B ^ 2 / (3 * r A B ^ 3) ∧
    iteratedDeriv 3 (P A B) 1 = 24 * A ^ 2 * B ^ 2 / r A B ^ 3 * (z3 - Real.pi ^ 2 / 3)

end X01
-- END STATEMENTS
namespace X01
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

/-- `J n k = ∫_0^1 t^n log^k t` -/
noncomputable def J (n k : ℕ) : ℝ := ∫ t in (0:ℝ)..1, t ^ n * log t ^ k

lemma ii_J (n k : ℕ) : IntervalIntegrable (fun t : ℝ => t ^ n * log t ^ k) volume 0 1 := by
  refine ii_of_ls (by fun_prop : Measurable (fun t : ℝ => t ^ n * log t ^ k)).aestronglyMeasurable
    1 k fun t ht => ?_
  rw [abs_mul, abs_pow, abs_pow, one_mul]
  have h1 : |t| ^ n ≤ 1 := pow_le_one₀ (abs_nonneg t) (by rw [abs_of_pos ht.1]; exact ht.2.le)
  have a1 : |log t| ≤ ls t := by unfold ls; linarith [abs_nonneg (log (1 - t))]
  have := pow_le_pow_left₀ (abs_nonneg _) a1 k
  have := pow_nonneg (abs_nonneg (log t)) k
  nlinarith

lemma J_rec (n k : ℕ) : (n + 1 : ℝ) * J n (k+1) + (k + 1) * J n k = 0 := by
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
  simpa [J] using this

lemma J0 (n : ℕ) : J n 0 = 1 / (n + 1) := by
  simp [J, integral_pow]

lemma J1 (n : ℕ) : J n 1 = -1 / (n + 1) ^ 2 := by
  have := J_rec n 0; rw [J0] at this; push_cast at this
  have hn : (n:ℝ) + 1 ≠ 0 := by positivity
  field_simp; field_simp at this; linarith

lemma J2 (n : ℕ) : J n 2 = 2 / (n + 1) ^ 3 := by
  have := J_rec n 1; rw [J1] at this; push_cast at this
  have hn : (n:ℝ) + 1 ≠ 0 := by positivity
  field_simp; field_simp at this; linarith

lemma J3 (n : ℕ) : J n 3 = -6 / (n + 1) ^ 4 := by
  have := J_rec n 2; rw [J2] at this; push_cast at this
  have hn : (n:ℝ) + 1 ≠ 0 := by positivity
  field_simp; field_simp at this; linarith

lemma hasSum_logmul (k : ℕ) :
    HasSum (fun n : ℕ => -(J (n+1) k) / (n+1)) (∫ t in (0:ℝ)..1, log t ^ k * log (1 - t)) := by
  have key : ∀ t ∈ Ioo (0:ℝ) 1, HasSum (fun n : ℕ => -(log t ^ k * (t ^ (n+1) / (n+1))))
      (log t ^ k * log (1 - t)) := by
    intro t ht
    have h := (hasSum_pow_div_log_of_abs_lt_one
      (show |t| < 1 by rw [abs_of_pos ht.1]; exact ht.2)).mul_left (log t ^ k)
    convert h.neg using 1; ring
  have hint : ∀ n : ℕ, ∫ t in (0:ℝ)..1, -(log t ^ k * (t ^ (n+1) / (n+1))) = -(J (n+1) k) / (n+1) := by
    intro n
    rw [intervalIntegral.integral_neg]
    have : (fun t : ℝ => log t ^ k * (t ^ (n+1) / (n+1))) =
        fun t => (1 / (n+1 : ℝ)) * (t ^ (n+1) * log t ^ k) := by funext t; ring
    rw [this, intervalIntegral.integral_const_mul, J]; ring
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

lemma int_a (k : ℕ) : ∫ t in (0:ℝ)..1, log t ^ k = J 0 k := by simp [J]

lemma int_symm (f : ℝ → ℝ) : ∫ t in (0:ℝ)..1, f (1 - t) = ∫ t in (0:ℝ)..1, f t := by
  rw [intervalIntegral.integral_comp_sub_left]; norm_num

lemma int_b (k : ℕ) : ∫ t in (0:ℝ)..1, log (1 - t) ^ k = J 0 k := by
  rw [int_symm (fun t => log t ^ k), int_a]

lemma int_abb : ∫ t in (0:ℝ)..1, log t * log (1 - t) ^ 2 = 2 * z3 - 6 + π ^ 2 / 3 := by
  rw [← int_aab, ← int_symm (fun t => log t ^ 2 * log (1 - t))]
  congr 1; funext t; simp only [sub_sub_cancel]; ring

/-- monomial `log^i t log^j (1-t)` -/
noncomputable def M (i j : ℕ) (t : ℝ) : ℝ := log t ^ i * log (1 - t) ^ j

lemma iM (i j : ℕ) : IntervalIntegrable (M i j) volume 0 1 := ii_mono i j

lemma iMc (c : ℝ) (i j : ℕ) : IntervalIntegrable (fun t => c * M i j t) volume 0 1 :=
  (ii_mono i j).const_mul c

lemma M10 : ∫ t in (0:ℝ)..1, M 1 0 t = -1 := by
  have := int_a 1; simp only [pow_one] at this
  simp only [M, pow_one, pow_zero, mul_one]; rw [this, J1]; norm_num
lemma M01 : ∫ t in (0:ℝ)..1, M 0 1 t = -1 := by
  have := int_b 1; simp only [pow_one] at this
  simp only [M, pow_one, pow_zero, one_mul]; rw [this, J1]; norm_num
lemma M20 : ∫ t in (0:ℝ)..1, M 2 0 t = 2 := by
  simp only [M, pow_zero, mul_one]; rw [int_a, J2]; norm_num
lemma M02 : ∫ t in (0:ℝ)..1, M 0 2 t = 2 := by
  simp only [M, pow_zero, one_mul]; rw [int_b, J2]; norm_num
lemma M30 : ∫ t in (0:ℝ)..1, M 3 0 t = -6 := by
  simp only [M, pow_zero, mul_one]; rw [int_a, J3]; norm_num
lemma M03 : ∫ t in (0:ℝ)..1, M 0 3 t = -6 := by
  simp only [M, pow_zero, one_mul]; rw [int_b, J3]; norm_num
lemma M11 : ∫ t in (0:ℝ)..1, M 1 1 t = 2 - π ^ 2 / 6 := by
  simp only [M, pow_one]; exact int_ab
lemma M21 : ∫ t in (0:ℝ)..1, M 2 1 t = 2 * z3 - 6 + π ^ 2 / 3 := by
  simp only [M, pow_one]; exact int_aab
lemma M12 : ∫ t in (0:ℝ)..1, M 1 2 t = 2 * z3 - 6 + π ^ 2 / 3 := by
  simp only [M, pow_one]; exact int_abb

lemma Lg_eq (t : ℝ) : Lg t = log t + log (1 - t) := by
  unfold Lg
  by_cases h0 : t = 0
  · simp [h0]
  by_cases h1 : 1 - t = 0
  · obtain rfl : t = 1 := by linarith
    simp
  exact log_mul h0 h1

lemma lg_eq (t : ℝ) : lg t = log t - log (1 - t) := by
  unfold lg
  by_cases h0 : t = 0
  · simp [h0]
  by_cases h1 : 1 - t = 0
  · obtain rfl : t = 1 := by linarith
    simp
  exact log_div h0 h1


lemma int4 (c1 c2 c3 c4 : ℝ) (i1 j1 i2 j2 i3 j3 i4 j4 : ℕ) :
    ∫ t in (0:ℝ)..1, (c1 * M i1 j1 t + c2 * M i2 j2 t + c3 * M i3 j3 t + c4 * M i4 j4 t) =
    c1 * (∫ t in (0:ℝ)..1, M i1 j1 t) + c2 * (∫ t in (0:ℝ)..1, M i2 j2 t) +
    c3 * (∫ t in (0:ℝ)..1, M i3 j3 t) + c4 * (∫ t in (0:ℝ)..1, M i4 j4 t) := by
  rw [intervalIntegral.integral_add (((iMc c1 i1 j1).add (iMc c2 i2 j2)).add (iMc c3 i3 j3))
    (iMc c4 i4 j4), intervalIntegral.integral_add ((iMc c1 i1 j1).add (iMc c2 i2 j2))
    (iMc c3 i3 j3), intervalIntegral.integral_add (iMc c1 i1 j1) (iMc c2 i2 j2)]
  simp only [intervalIntegral.integral_const_mul]

lemma Lg_symm (t : ℝ) : Lg (1 - t) = Lg t := by
  rw [Lg_eq, Lg_eq, sub_sub_cancel]; ring

lemma lg_symm (t : ℝ) : lg (1 - t) = -lg t := by
  rw [lg_eq, lg_eq, sub_sub_cancel]; ring

lemma abs_Lg_le (t : ℝ) : |Lg t| ≤ ls t := by rw [Lg_eq]; exact abs_add_le _ _
lemma abs_lg_le (t : ℝ) : |lg t| ≤ ls t := by rw [lg_eq]; exact abs_sub _ _


/-- `S(η) = A² e^{2ηa} + B² e^{2ηb}` -/
noncomputable def Sx (A B a b η : ℝ) : ℝ := A ^ 2 * exp (2 * η * a) + B ^ 2 * exp (2 * η * b)
/-- `N_k(η) = A² (2a)^k e^{2ηa} + B² (2b)^k e^{2ηb}` -/
noncomputable def Nx (A B a b : ℝ) (k : ℕ) (η : ℝ) : ℝ :=
  A ^ 2 * (2 * a) ^ k * exp (2 * η * a) + B ^ 2 * (2 * b) ^ k * exp (2 * η * b)
noncomputable def mx (A B a b : ℝ) (k : ℕ) (η : ℝ) : ℝ := Nx A B a b k η / Sx A B a b η

/-- the `η`-derivatives of `√S` -/
noncomputable def DD (A B a b : ℝ) : ℕ → ℝ → ℝ
  | 0, η => sqrt (Sx A B a b η)
  | 1, η => sqrt (Sx A B a b η) * (mx A B a b 1 η / 2)
  | 2, η => sqrt (Sx A B a b η) * (mx A B a b 2 η / 2 - mx A B a b 1 η ^ 2 / 4)
  | _, η => sqrt (Sx A B a b η) *
      (mx A B a b 3 η / 2 - 3 * mx A B a b 1 η * mx A B a b 2 η / 4 + 3 * mx A B a b 1 η ^ 3 / 8)

section alg
variable {A B : ℝ} (a b : ℝ)

lemma Sx_pos (hA : 0 < A) (η : ℝ) : 0 < Sx A B a b η := by
  unfold Sx; have := exp_pos (2 * η * a); positivity

lemma Nx_zero (η : ℝ) : Nx A B a b 0 η = Sx A B a b η := by simp [Nx, Sx]

lemma hd_exp (c η : ℝ) : HasDerivAt (fun x => exp (2 * x * c)) (exp (2 * η * c) * (2 * c)) η := by
  have : HasDerivAt (fun x : ℝ => 2 * x * c) (2 * c) η := by
    simpa using ((hasDerivAt_id η).const_mul 2).mul_const c
  exact this.exp

lemma hd_Nx (k : ℕ) (η : ℝ) : HasDerivAt (Nx A B a b k) (Nx A B a b (k+1) η) η := by
  have h := ((hd_exp a η).const_mul (A ^ 2 * (2 * a) ^ k)).add
    ((hd_exp b η).const_mul (B ^ 2 * (2 * b) ^ k))
  unfold Nx
  convert h using 1
  ring

lemma hd_Sx (η : ℝ) : HasDerivAt (Sx A B a b) (Nx A B a b 1 η) η := by
  have : Sx A B a b = Nx A B a b 0 := by funext x; rw [Nx_zero]
  rw [this]; exact hd_Nx a b 0 η

lemma hd_mx (hA : 0 < A) (k : ℕ) (η : ℝ) :
    HasDerivAt (mx A B a b k) (mx A B a b (k+1) η - mx A B a b k η * mx A B a b 1 η) η := by
  have hS := (Sx_pos a b (B := B) hA η).ne'
  have h := (hd_Nx (A := A) (B := B) a b k η).div (hd_Sx (A := A) (B := B) a b η) hS
  unfold mx
  convert h using 1
  field_simp

lemma hd_sqrt (hA : 0 < A) (η : ℝ) :
    HasDerivAt (fun x => sqrt (Sx A B a b x)) (sqrt (Sx A B a b η) * (mx A B a b 1 η / 2)) η := by
  have hS := Sx_pos a b (B := B) hA η
  have h := (hd_Sx (A := A) (B := B) a b η).sqrt hS.ne'
  convert h using 1
  have hu := sq_sqrt hS.le
  have hu0 : 0 < sqrt (Sx A B a b η) := sqrt_pos.mpr hS
  unfold mx
  field_simp
  rw [hu]; ring

lemma hd_DD (hA : 0 < A) (k : ℕ) (hk : k < 3) (η : ℝ) :
    HasDerivAt (DD A B a b k) (DD A B a b (k+1) η) η := by
  have h0 := hd_sqrt a b (B := B) hA η
  have m1 := hd_mx a b (B := B) hA 1 η
  have m2 := hd_mx a b (B := B) hA 2 η
  have m3 := hd_mx a b (B := B) hA 3 η
  interval_cases k
  · exact h0
  · have h := h0.fun_mul (m1.div_const 2)
    have e1 : DD A B a b 1 = fun x => sqrt (Sx A B a b x) * (mx A B a b 1 x / 2) := rfl
    have e2 : DD A B a b (1+1) η =
        sqrt (Sx A B a b η) * (mx A B a b 2 η / 2 - mx A B a b 1 η ^ 2 / 4) := rfl
    rw [e1, e2]
    convert h using 1
    ring
  · have h := h0.fun_mul ((m2.div_const 2).fun_sub ((m1.fun_pow 2).div_const 4))
    have e1 : DD A B a b 2 = fun x => sqrt (Sx A B a b x) *
        (mx A B a b 2 x / 2 - mx A B a b 1 x ^ 2 / 4) := rfl
    have e2 : DD A B a b (2+1) η = sqrt (Sx A B a b η) *
      (mx A B a b 3 η / 2 - 3 * mx A B a b 1 η * mx A B a b 2 η / 4 +
        3 * mx A B a b 1 η ^ 3 / 8) := rfl
    rw [e1, e2]
    convert h using 1
    push_cast; ring

lemma abs_mx_le (hA : 0 < A) (k : ℕ) (η : ℝ) :
    |mx A B a b k η| ≤ (2 * (|a| + |b|)) ^ k := by
  have hS := Sx_pos a b (B := B) hA η
  unfold mx
  rw [abs_div, abs_of_pos hS, div_le_iff₀ hS]
  unfold Nx Sx
  have ea := exp_pos (2 * η * a)
  have eb := exp_pos (2 * η * b)
  have ha : |2 * a| ≤ 2 * (|a| + |b|) := by rw [abs_mul]; norm_num
  have hb : |2 * b| ≤ 2 * (|a| + |b|) := by rw [abs_mul]; norm_num
  calc _ ≤ |A ^ 2 * (2 * a) ^ k * exp (2 * η * a)| + |B ^ 2 * (2 * b) ^ k * exp (2 * η * b)| :=
        abs_add_le _ _
    _ = A ^ 2 * |2 * a| ^ k * exp (2 * η * a) + B ^ 2 * |2 * b| ^ k * exp (2 * η * b) := by
        simp only [abs_mul, abs_pow, sq_abs, abs_of_pos ea, abs_of_pos eb]
    _ ≤ A ^ 2 * (2 * (|a| + |b|)) ^ k * exp (2 * η * a) +
        B ^ 2 * (2 * (|a| + |b|)) ^ k * exp (2 * η * b) := by gcongr
    _ = _ := by ring

lemma abs_DD_le (hA : 0 < A) (k : ℕ) (hk : k ≤ 3) (η : ℝ) :
    |DD A B a b k η| ≤ 13 * (|a| + |b|) ^ k * sqrt (Sx A B a b η) := by
  have q0 := sqrt_nonneg (Sx A B a b η)
  set s := |a| + |b| with hs
  have s0 : 0 ≤ s := by positivity
  have b1 := abs_mx_le a b (B := B) hA 1 η
  have b2 := abs_mx_le a b (B := B) hA 2 η
  have b3 := abs_mx_le a b (B := B) hA 3 η
  rw [← hs] at b1 b2 b3
  set x1 := mx A B a b 1 η
  set x2 := mx A B a b 2 η
  set x3 := mx A B a b 3 η
  have a1 : |x1| ≤ 2 * s := by simpa using b1
  have a2 : |x2| ≤ 4 * s ^ 2 := by convert b2 using 1; ring
  have a3 : |x3| ≤ 8 * s ^ 3 := by convert b3 using 1; ring
  interval_cases k
  · simp only [DD, abs_of_nonneg q0]; nlinarith
  · simp only [DD]
    rw [abs_mul, abs_of_nonneg q0, abs_div, abs_two]
    have : |x1| / 2 ≤ 13 * s ^ 1 := by linarith
    calc _ ≤ sqrt (Sx A B a b η) * (13 * s ^ 1) := by gcongr
      _ = _ := by ring
  · simp only [DD]
    rw [abs_mul, abs_of_nonneg q0]
    have e : |x2 / 2 - x1 ^ 2 / 4| ≤ 13 * s ^ 2 := by
      have h1 : |x2 / 2 - x1 ^ 2 / 4| ≤ |x2| / 2 + |x1| ^ 2 / 4 := by
        calc _ ≤ |x2 / 2| + |x1 ^ 2 / 4| := abs_sub _ _
          _ = _ := by rw [abs_div, abs_div, abs_pow]; norm_num
      have h2 : |x1| ^ 2 ≤ (2 * s) ^ 2 := pow_le_pow_left₀ (abs_nonneg _) a1 2
      nlinarith
    calc _ ≤ sqrt (Sx A B a b η) * (13 * s ^ 2) := by gcongr
      _ = _ := by ring
  · simp only [DD]
    rw [abs_mul, abs_of_nonneg q0]
    have e : |x3 / 2 - 3 * x1 * x2 / 4 + 3 * x1 ^ 3 / 8| ≤ 13 * s ^ 3 := by
      have h1 : |x3 / 2 - 3 * x1 * x2 / 4 + 3 * x1 ^ 3 / 8| ≤
          |x3| / 2 + 3 * (|x1| * |x2|) / 4 + 3 * |x1| ^ 3 / 8 := by
        calc _ ≤ |x3 / 2 - 3 * x1 * x2 / 4| + |3 * x1 ^ 3 / 8| := abs_add_le _ _
          _ ≤ |x3 / 2| + |3 * x1 * x2 / 4| + |3 * x1 ^ 3 / 8| := by
              gcongr; exact abs_sub _ _
          _ = _ := by
              rw [abs_div, abs_div, abs_div, abs_mul, abs_mul, abs_mul, abs_pow]; norm_num; ring
      have h2 : |x1| ^ 3 ≤ (2 * s) ^ 3 := pow_le_pow_left₀ (abs_nonneg _) a1 3
      have h3 : |x1| * |x2| ≤ (2 * s) * (4 * s ^ 2) :=
        mul_le_mul a1 a2 (abs_nonneg _) (by positivity)
      nlinarith
    calc _ ≤ sqrt (Sx A B a b η) * (13 * s ^ 3) := by gcongr
      _ = _ := by ring

end alg



section integ
variable {A B : ℝ}

/-- `F_k(η) = ∫_0^1 D_k(η, log t, log(1-t)) dt` -/
noncomputable def FF (A B : ℝ) (k : ℕ) (η : ℝ) : ℝ :=
  ∫ t in (0:ℝ)..1, DD A B (log t) (log (1 - t)) k η

/-- the dominating function -/
noncomputable def bnd (A B : ℝ) (k : ℕ) (t : ℝ) : ℝ :=
  13 * (A + B) * (4 * ((k:ℝ) + 1)) ^ k * (t * (1 - t)) ^ (-(3/4:ℝ))

lemma ii_bnd (k : ℕ) : IntervalIntegrable (bnd A B k) volume 0 1 :=
  (ii_tt (-(3/4)) (by norm_num)).const_mul _

lemma sqrtS_le (hA : 0 < A) (hB : 0 < B) {t η : ℝ} (ht : t ∈ Ioo (0:ℝ) 1) (hη : |η| ≤ 1/4) :
    sqrt (Sx A B (log t) (log (1 - t)) η) ≤ (A + B) * (t * (1 - t)) ^ (-(1/4:ℝ)) := by
  have hp : 0 < t * (1 - t) := mul_pos ht.1 (by linarith [ht.2])
  have hη' := abs_le.mp hη
  have key : ∀ u : ℝ, 0 < u → u ≤ 1 → t * (1 - t) ≤ u →
      exp (η * log u) ≤ (t * (1 - t)) ^ (-(1/4:ℝ)) := by
    intro u hu hu1 hle
    have hlog : log u ≤ 0 := log_nonpos hu.le hu1
    have h1 : η * log u ≤ (-(1/4)) * log u := by nlinarith
    calc exp (η * log u) ≤ exp ((-(1/4)) * log u) := exp_le_exp.mpr h1
      _ = u ^ (-(1/4:ℝ)) := by rw [rpow_def_of_pos hu]; ring_nf
      _ ≤ _ := rpow_le_rpow_of_nonpos hp hle (by norm_num)
  have e1 := key t ht.1 ht.2.le (by nlinarith [ht.1, ht.2])
  have e2 := key (1 - t) (by linarith [ht.2]) (by linarith [ht.1]) (by nlinarith [ht.1, ht.2])
  have x1 := exp_pos (η * log t)
  have x2 := exp_pos (η * log (1 - t))
  have hS : Sx A B (log t) (log (1 - t)) η ≤
      (A * exp (η * log t) + B * exp (η * log (1 - t))) ^ 2 := by
    unfold Sx
    rw [show 2 * η * log t = η * log t + η * log t by ring,
      show 2 * η * log (1 - t) = η * log (1 - t) + η * log (1 - t) by ring, exp_add, exp_add]
    nlinarith [mul_pos (mul_pos hA hB) (mul_pos x1 x2)]
  calc sqrt (Sx A B (log t) (log (1 - t)) η)
      ≤ sqrt ((A * exp (η * log t) + B * exp (η * log (1 - t))) ^ 2) := sqrt_le_sqrt hS
    _ = A * exp (η * log t) + B * exp (η * log (1 - t)) := sqrt_sq (by positivity)
    _ ≤ A * (t * (1 - t)) ^ (-(1/4:ℝ)) + B * (t * (1 - t)) ^ (-(1/4:ℝ)) := by gcongr
    _ = _ := by ring

lemma DD_bound (hA : 0 < A) (hB : 0 < B) (k : ℕ) (hk : k ≤ 3) {t η : ℝ} (ht : t ∈ Ioo (0:ℝ) 1)
    (hη : |η| ≤ 1/4) : |DD A B (log t) (log (1 - t)) k η| ≤ bnd A B k t := by
  have h1 := abs_DD_le (log t) (log (1 - t)) (B := B) hA k hk η
  have h2 := sqrtS_le hA hB ht hη
  have h3 := ls_pow_le ht k
  have hp : 0 < t * (1 - t) := mul_pos ht.1 (by linarith [ht.2])
  have e : (t * (1 - t)) ^ (-(3/4:ℝ)) =
      (t * (1 - t)) ^ (-(1/2:ℝ)) * (t * (1 - t)) ^ (-(1/4:ℝ)) := by
    rw [← rpow_add hp]; norm_num
  unfold ls at h3
  have q0 := sqrt_nonneg (Sx A B (log t) (log (1 - t)) η)
  have q1 : 0 ≤ (|log t| + |log (1 - t)|) ^ k := by positivity
  unfold bnd
  calc _ ≤ 13 * (|log t| + |log (1 - t)|) ^ k * sqrt (Sx A B (log t) (log (1 - t)) η) := h1
    _ ≤ 13 * ((4 * ((k:ℝ) + 1)) ^ k * (t * (1 - t)) ^ (-(1/2:ℝ))) *
        ((A + B) * (t * (1 - t)) ^ (-(1/4:ℝ))) := by gcongr
    _ = _ := by rw [e]; ring

lemma meas_DD (k : ℕ) (η : ℝ) :
    AEStronglyMeasurable (fun t => DD A B (log t) (log (1 - t)) k η)
      (volume.restrict (Ι (0:ℝ) 1)) := by
  apply Measurable.aestronglyMeasurable
  rcases k with _ | _ | _ | k <;> simp only [DD, mx, Nx, Sx] <;> fun_prop

lemma ii_DD (hA : 0 < A) (hB : 0 < B) (k : ℕ) (hk : k ≤ 3) {η : ℝ} (hη : |η| ≤ 1/4) :
    IntervalIntegrable (fun t => DD A B (log t) (log (1 - t)) k η) volume 0 1 := by
  refine IntervalIntegrable.mono_fun' (ii_bnd (A := A) (B := B) k) (meas_DD k η) ?_
  filter_upwards [ae_Ioo01] with t ht
  exact DD_bound hA hB k hk ht hη

/-- the parameter window `|η| < 1/4` -/
noncomputable def U : Set ℝ := Ioo (-(1/4)) (1/4)

lemma U_open : IsOpen U := isOpen_Ioo

lemma abs_of_U {η : ℝ} (h : η ∈ U) : |η| ≤ 1/4 := abs_le.mpr ⟨h.1.le, h.2.le⟩

lemma hd_FF (hA : 0 < A) (hB : 0 < B) (k : ℕ) (hk : k < 3) {η₀ : ℝ} (h₀ : η₀ ∈ U) :
    HasDerivAt (FF A B k) (FF A B (k+1) η₀) η₀ := by
  have := intervalIntegral.hasDerivAt_integral_of_dominated_loc_of_deriv_le
    (F := fun η t => DD A B (log t) (log (1 - t)) k η)
    (F' := fun η t => DD A B (log t) (log (1 - t)) (k+1) η) (bound := bnd A B (k+1))
    (U_open.mem_nhds h₀) (Eventually.of_forall fun η => meas_DD k η)
    (ii_DD hA hB k hk.le (abs_of_U h₀)) (meas_DD (k+1) η₀) ?_ (ii_bnd (k+1)) ?_
  · exact this.2
  · filter_upwards [ae_Ioo01'] with t ht ht0 η hη
    exact DD_bound hA hB (k+1) (by omega) (ht ht0) (abs_of_U hη)
  · exact Eventually.of_forall fun t _ η _ => hd_DD (log t) (log (1 - t)) (B := B) hA k hk η

lemma cont_mx (hA : 0 < A) (a b : ℝ) (k : ℕ) : Continuous (mx A B a b k) :=
  Continuous.div (show Continuous (Nx A B a b k) by unfold Nx; fun_prop)
    (show Continuous (Sx A B a b) by unfold Sx; fun_prop)
    (fun x => (Sx_pos a b (B := B) hA x).ne')

lemma cont_DD3 (hA : 0 < A) (a b : ℝ) : Continuous (DD A B a b 3) := by
  have h1 := cont_mx (B := B) hA a b 1
  have h2 := cont_mx (B := B) hA a b 2
  have h3 := cont_mx (B := B) hA a b 3
  have hs : Continuous (fun x => sqrt (Sx A B a b x)) := by unfold Sx; fun_prop
  show Continuous (fun η => sqrt (Sx A B a b η) *
      (mx A B a b 3 η / 2 - 3 * mx A B a b 1 η * mx A B a b 2 η / 4 + 3 * mx A B a b 1 η ^ 3 / 8))
  fun_prop

lemma cont_FF3 (hA : 0 < A) (hB : 0 < B) : ContinuousOn (FF A B 3) U := by
  intro η₀ h₀
  apply ContinuousAt.continuousWithinAt
  apply intervalIntegral.continuousAt_of_dominated_interval (bound := bnd A B 3)
    (Eventually.of_forall fun η => meas_DD 3 η)
  · filter_upwards [U_open.mem_nhds h₀] with η hη
    filter_upwards [ae_Ioo01'] with t ht ht0
    exact DD_bound hA hB 3 le_rfl (ht ht0) (abs_of_U hη)
  · exact ii_bnd 3
  · exact Eventually.of_forall fun t _ => (cont_DD3 (B := B) hA _ _).continuousAt

end integ


section chain

lemma cd_step {f g : ℝ → ℝ} {s : Set ℝ} {n : ℕ} (hs : IsOpen s)
    (hd : ∀ x ∈ s, HasDerivAt f (g x) x) (hg : ContDiffOn ℝ n g s) :
    ContDiffOn ℝ ((n : WithTop ℕ∞) + 1) f s := by
  rw [contDiffOn_succ_iff_deriv_of_isOpen hs]
  refine ⟨fun x hx => (hd x hx).differentiableAt.differentiableWithinAt, fun h => ?_,
    hg.congr fun x hx => (hd x hx).deriv⟩
  exact absurd h (by simp)

lemma chain3 {f0 f1 f2 f3 : ℝ → ℝ} {s : Set ℝ} (hs : IsOpen s)
    (h0 : ∀ x ∈ s, HasDerivAt f0 (f1 x) x) (h1 : ∀ x ∈ s, HasDerivAt f1 (f2 x) x)
    (h2 : ∀ x ∈ s, HasDerivAt f2 (f3 x) x) (h3 : ContinuousOn f3 s) :
    ContDiffOn ℝ 3 f0 s ∧ ∀ x ∈ s, deriv f0 x = f1 x ∧ iteratedDeriv 2 f0 x = f2 x ∧
      iteratedDeriv 3 f0 x = f3 x := by
  have c2 : ContDiffOn ℝ ((0:ℕ) + 1 : WithTop ℕ∞) f2 s :=
    cd_step hs h2 (by simpa using h3)
  have c1 : ContDiffOn ℝ ((1:ℕ) + 1 : WithTop ℕ∞) f1 s := cd_step hs h1 (by simpa using c2)
  have e2 : ((1:ℕ) + 1 : WithTop ℕ∞) = ((2:ℕ) : WithTop ℕ∞) := by norm_num
  have e3 : ((2:ℕ) + 1 : WithTop ℕ∞) = 3 := by norm_num
  rw [e2] at c1
  have c0 : ContDiffOn ℝ ((2:ℕ) + 1 : WithTop ℕ∞) f0 s := cd_step hs h0 c1
  rw [e3] at c0
  have d1 : ∀ x ∈ s, deriv f0 x = f1 x := fun x hx => (h0 x hx).deriv
  have d2 : ∀ x ∈ s, iteratedDeriv 2 f0 x = f2 x := by
    intro x hx
    rw [iteratedDeriv_succ, iteratedDeriv_one]
    rw [Filter.EventuallyEq.deriv_eq (f := f1)
      (by filter_upwards [hs.mem_nhds hx] with y hy using d1 y hy)]
    exact (h1 x hx).deriv
  have d3 : ∀ x ∈ s, iteratedDeriv 3 f0 x = f3 x := by
    intro x hx
    rw [iteratedDeriv_succ]
    rw [Filter.EventuallyEq.deriv_eq (f := f2)
      (by filter_upwards [hs.mem_nhds hx] with y hy using d2 y hy)]
    exact (h2 x hx).deriv
  exact ⟨c0, fun x hx => ⟨d1 x hx, d2 x hx, d3 x hx⟩⟩

end chain

section Ifun
variable {A B : ℝ}

lemma r_pos (hA : 0 < A) : 0 < r A B := sqrt_pos.mpr (by positivity)

lemma r_sq (A B : ℝ) : r A B ^ 2 = A ^ 2 + B ^ 2 := sq_sqrt (by positivity)

noncomputable def I0 (A B η : ℝ) : ℝ := (1 + η) * FF A B 0 η / r A B
noncomputable def I1 (A B η : ℝ) : ℝ := (FF A B 0 η + (1 + η) * FF A B 1 η) / r A B
noncomputable def I2 (A B η : ℝ) : ℝ := (2 * FF A B 1 η + (1 + η) * FF A B 2 η) / r A B
noncomputable def I3 (A B η : ℝ) : ℝ := (3 * FF A B 2 η + (1 + η) * FF A B 3 η) / r A B

lemma hd_I (hA : 0 < A) (hB : 0 < B) {x : ℝ} (hx : x ∈ U) :
    HasDerivAt (I0 A B) (I1 A B x) x ∧ HasDerivAt (I1 A B) (I2 A B x) x ∧
      HasDerivAt (I2 A B) (I3 A B x) x := by
  have f0 := hd_FF hA hB 0 (by norm_num) hx
  have f1 := hd_FF hA hB 1 (by norm_num) hx
  have f2 := hd_FF hA hB 2 (by norm_num) hx
  have e := (hasDerivAt_id' x).const_add (1:ℝ)
  refine ⟨?_, ?_, ?_⟩
  · have h := (e.fun_mul f0).div_const (r A B)
    unfold I0 I1; convert h using 1; ring
  · have h := (f0.add (e.fun_mul f1)).div_const (r A B)
    unfold I1 I2; convert h using 1; ring
  · have h := ((f1.const_mul 2).add (e.fun_mul f2)).div_const (r A B)
    unfold I2 I3; convert h using 1; ring

lemma cont_I3 (hA : 0 < A) (hB : 0 < B) : ContinuousOn (I3 A B) U := by
  have h2 : ContinuousOn (FF A B 2) U := fun x hx =>
    (hd_FF hA hB 2 (by norm_num) hx).continuousAt.continuousWithinAt
  have h3 := cont_FF3 hA hB
  have : ContinuousOn (fun η => (3 * FF A B 2 η + (1 + η) * FF A B 3 η) / r A B) U :=
    ((continuousOn_const.mul h2).add ((continuousOn_const.add continuousOn_id).mul h3)).div_const _
  exact this

lemma P_eq (hη : -1 < η) : P A B (1 / (1 + η)) = 4 * (1 + η) * FF A B 0 η := by
  have hp : 0 < 1 + η := by linarith
  unfold P FF
  rw [show 4 / (1 / (1 + η)) = 4 * (1 + η) by field_simp]
  congr 1
  apply intervalIntegral.integral_congr_uIoo
  rw [uIoo01]; intro t ht
  have h1 : 0 < 1 - t := by linarith [ht.2]
  simp only [integrand, DD, Sx]
  rw [show 2 / (1 / (1 + η)) - 2 = 2 * η by field_simp; ring,
    rpow_def_of_pos ht.1, rpow_def_of_pos h1]
  ring_nf

lemma I_eq (hA : 0 < A) {η : ℝ} (hη : η ∈ U) : I A B η = I0 A B η := by
  have h : -1 < η := by have := hη.1; linarith
  have hr := r_pos (B := B) hA
  unfold I I0
  rw [P_eq h]; field_simp

lemma hd_I_real (hA : 0 < A) (hB : 0 < B) {x : ℝ} (hx : x ∈ U) :
    HasDerivAt (I A B) (I1 A B x) x :=
  (hd_I hA hB hx).1.congr_of_eventuallyEq
    (by filter_upwards [U_open.mem_nhds hx] with y hy using I_eq hA hy)

lemma chainI (hA : 0 < A) (hB : 0 < B) :
    ContDiffOn ℝ 3 (I A B) U ∧ ∀ x ∈ U, deriv (I A B) x = I1 A B x ∧
      iteratedDeriv 2 (I A B) x = I2 A B x ∧ iteratedDeriv 3 (I A B) x = I3 A B x :=
  chain3 U_open (fun x hx => hd_I_real hA hB hx) (fun x hx => (hd_I hA hB hx).2.1)
    (fun x hx => (hd_I hA hB hx).2.2) (cont_I3 hA hB)

end Ifun


section values
variable {A B : ℝ}

lemma FF0_0 : FF A B 0 0 = r A B := by
  simp [FF, DD, Sx, r]

lemma FF1_0 (hA : 0 < A) :
    FF A B 1 0 = r A B * A ^ 2 / (A ^ 2 + B ^ 2) * (-1) + r A B * B ^ 2 / (A ^ 2 + B ^ 2) * (-1) := by
  have hR : A ^ 2 + B ^ 2 ≠ 0 := by positivity
  have e : FF A B 1 0 = ∫ t in (0:ℝ)..1, (r A B * A ^ 2 / (A ^ 2 + B ^ 2) * M 1 0 t +
      r A B * B ^ 2 / (A ^ 2 + B ^ 2) * M 0 1 t + 0 * M 0 0 t + 0 * M 0 0 t) := by
    unfold FF
    refine intervalIntegral.integral_congr fun t _ => ?_
    simp only [DD, mx, Nx, Sx, M, mul_zero, zero_mul, exp_zero, mul_one, pow_zero, pow_one]
    rw [show sqrt (A ^ 2 + B ^ 2) = r A B from rfl]
    field_simp
    ring
  rw [e, int4, M10, M01]; ring

lemma FF2_0 (hA : 0 < A) :
    FF A B 2 0 = r A B * (2 * A ^ 2 / (A ^ 2 + B ^ 2) - (A ^ 2) ^ 2 / (A ^ 2 + B ^ 2) ^ 2) * 2 +
      r A B * (2 * B ^ 2 / (A ^ 2 + B ^ 2) - (B ^ 2) ^ 2 / (A ^ 2 + B ^ 2) ^ 2) * 2 +
      (-2 * r A B * A ^ 2 * B ^ 2 / (A ^ 2 + B ^ 2) ^ 2) * (2 - π ^ 2 / 6) := by
  have hR : A ^ 2 + B ^ 2 ≠ 0 := by positivity
  have e : FF A B 2 0 = ∫ t in (0:ℝ)..1,
      (r A B * (2 * A ^ 2 / (A ^ 2 + B ^ 2) - (A ^ 2) ^ 2 / (A ^ 2 + B ^ 2) ^ 2) * M 2 0 t +
      r A B * (2 * B ^ 2 / (A ^ 2 + B ^ 2) - (B ^ 2) ^ 2 / (A ^ 2 + B ^ 2) ^ 2) * M 0 2 t +
      (-2 * r A B * A ^ 2 * B ^ 2 / (A ^ 2 + B ^ 2) ^ 2) * M 1 1 t + 0 * M 0 0 t) := by
    unfold FF
    refine intervalIntegral.integral_congr fun t _ => ?_
    simp only [DD, mx, Nx, Sx, M, mul_zero, zero_mul, exp_zero, mul_one, pow_zero, pow_one]
    rw [show sqrt (A ^ 2 + B ^ 2) = r A B from rfl]
    field_simp
    ring
  rw [e, int4, M20, M02, M11]; ring

lemma FF3_0 (hA : 0 < A) :
    FF A B 3 0 =
      r A B * (4 * A ^ 2 / (A ^ 2 + B ^ 2) - 6 * (A ^ 2) ^ 2 / (A ^ 2 + B ^ 2) ^ 2 +
        3 * (A ^ 2) ^ 3 / (A ^ 2 + B ^ 2) ^ 3) * (-6) +
      r A B * (4 * B ^ 2 / (A ^ 2 + B ^ 2) - 6 * (B ^ 2) ^ 2 / (A ^ 2 + B ^ 2) ^ 2 +
        3 * (B ^ 2) ^ 3 / (A ^ 2 + B ^ 2) ^ 3) * (-6) +
      r A B * (-6 * A ^ 2 * B ^ 2 / (A ^ 2 + B ^ 2) ^ 2 +
        9 * (A ^ 2) ^ 2 * B ^ 2 / (A ^ 2 + B ^ 2) ^ 3) * (2 * z3 - 6 + π ^ 2 / 3) +
      r A B * (-6 * A ^ 2 * B ^ 2 / (A ^ 2 + B ^ 2) ^ 2 +
        9 * A ^ 2 * (B ^ 2) ^ 2 / (A ^ 2 + B ^ 2) ^ 3) * (2 * z3 - 6 + π ^ 2 / 3) := by
  have hR : A ^ 2 + B ^ 2 ≠ 0 := by positivity
  have e : FF A B 3 0 = ∫ t in (0:ℝ)..1,
      (r A B * (4 * A ^ 2 / (A ^ 2 + B ^ 2) - 6 * (A ^ 2) ^ 2 / (A ^ 2 + B ^ 2) ^ 2 +
        3 * (A ^ 2) ^ 3 / (A ^ 2 + B ^ 2) ^ 3) * M 3 0 t +
      r A B * (4 * B ^ 2 / (A ^ 2 + B ^ 2) - 6 * (B ^ 2) ^ 2 / (A ^ 2 + B ^ 2) ^ 2 +
        3 * (B ^ 2) ^ 3 / (A ^ 2 + B ^ 2) ^ 3) * M 0 3 t +
      r A B * (-6 * A ^ 2 * B ^ 2 / (A ^ 2 + B ^ 2) ^ 2 +
        9 * (A ^ 2) ^ 2 * B ^ 2 / (A ^ 2 + B ^ 2) ^ 3) * M 2 1 t +
      r A B * (-6 * A ^ 2 * B ^ 2 / (A ^ 2 + B ^ 2) ^ 2 +
        9 * A ^ 2 * (B ^ 2) ^ 2 / (A ^ 2 + B ^ 2) ^ 3) * M 1 2 t) := by
    unfold FF
    refine intervalIntegral.integral_congr fun t _ => ?_
    simp only [DD, mx, Nx, Sx, M, mul_zero, zero_mul, exp_zero, mul_one, pow_zero, pow_one]
    rw [show sqrt (A ^ 2 + B ^ 2) = r A B from rfl]
    field_simp
    ring
  rw [e, int4, M30, M03, M21, M12]

end values

end X01

theorem X01.proof_6 : X01.stmt_6 := by
  intro j k
  have hm : Measurable (fun t : ℝ => X01.Lg t ^ j * X01.lg t ^ k) := by
    unfold X01.Lg X01.lg; fun_prop
  constructor
  · refine X01.ii_of_ls hm.aestronglyMeasurable 1 (j + k) fun t _ => ?_
    rw [abs_mul, abs_pow, abs_pow, one_mul, pow_add]
    have := X01.ls_nonneg t
    gcongr
    · exact X01.abs_Lg_le t
    · exact X01.abs_lg_le t
  · intro hk
    have h := X01.int_symm (fun t => X01.Lg t ^ j * X01.lg t ^ k)
    simp only [X01.Lg_symm, X01.lg_symm, hk.neg_pow, mul_neg, intervalIntegral.integral_neg] at h
    linarith

theorem X01.proof_7 : X01.stmt_7 := by
  open X01 in
  refine ⟨?_, ?_, ?_⟩
  · have e : ∫ t in (0:ℝ)..1, Lg t =
        ∫ t in (0:ℝ)..1, (1 * M 1 0 t + 1 * M 0 1 t + 0 * M 0 0 t + 0 * M 0 0 t) :=
      intervalIntegral.integral_congr fun t _ => by simp only [M, Lg_eq]; ring
    rw [e, int4, M10, M01]; ring
  · have e : ∫ t in (0:ℝ)..1, Lg t ^ 2 =
        ∫ t in (0:ℝ)..1, (1 * M 2 0 t + 2 * M 1 1 t + 1 * M 0 2 t + 0 * M 0 0 t) :=
      intervalIntegral.integral_congr fun t _ => by simp only [M, Lg_eq]; ring
    rw [e, int4, M20, M11, M02]; ring
  · have e : ∫ t in (0:ℝ)..1, Lg t ^ 3 =
        ∫ t in (0:ℝ)..1, (1 * M 3 0 t + 3 * M 2 1 t + 3 * M 1 2 t + 1 * M 0 3 t) :=
      intervalIntegral.integral_congr fun t _ => by simp only [M, Lg_eq]; ring
    rw [e, int4, M30, M21, M12, M03]; ring

theorem X01.proof_8 : X01.stmt_8 := by
  open X01 in
  have e : ∫ t in (0:ℝ)..1, lg t ^ 2 =
      ∫ t in (0:ℝ)..1, (1 * M 2 0 t + (-2) * M 1 1 t + 1 * M 0 2 t + 0 * M 0 0 t) :=
    intervalIntegral.integral_congr fun t _ => by simp only [M, lg_eq]; ring
  unfold stmt_8; rw [e, int4, M20, M11, M02]; ring

theorem X01.proof_9 : X01.stmt_9 := by
  open X01 in
  have e : ∫ t in (0:ℝ)..1, Lg t * lg t ^ 2 =
      ∫ t in (0:ℝ)..1, (1 * M 3 0 t + (-1) * M 2 1 t + (-1) * M 1 2 t + 1 * M 0 3 t) :=
    intervalIntegral.integral_congr fun t _ => by simp only [M, Lg_eq, lg_eq]; ring
  unfold stmt_9; rw [e, int4, M30, M21, M12, M03]; ring

theorem X01.proof_10 : X01.stmt_10 := by
  intro A B hA hB
  have h0U : (0:ℝ) ∈ X01.U := by unfold X01.U; constructor <;> norm_num
  obtain ⟨hc, hv⟩ := X01.chainI hA hB
  obtain ⟨d1, d2, d3⟩ := hv 0 h0U
  have hr := X01.r_pos (B := B) hA
  have hr2 := X01.r_sq A B
  have hB2 : B ^ 2 = X01.r A B ^ 2 - A ^ 2 := by linarith
  refine ⟨hc.contDiffAt (X01.U_open.mem_nhds h0U), ?_, ?_, ?_, ?_⟩
  · rw [X01.I_eq hA h0U]; unfold X01.I0; rw [X01.FF0_0]; field_simp; norm_num
  · rw [d1]; unfold X01.I1; rw [X01.FF0_0, X01.FF1_0 hA, ← hr2, hB2]
    field_simp; ring
  · rw [d2]; unfold X01.I2; rw [X01.FF1_0 hA, X01.FF2_0 hA]; unfold X01.q
    rw [← hr2, hB2]
    field_simp; ring
  · rw [d3]; unfold X01.I3; rw [X01.FF2_0 hA, X01.FF3_0 hA]; unfold X01.q
    rw [← hr2, hB2]
    field_simp; ring

namespace X01
open MeasureTheory Set Filter Real Interval Topology

section Pfun
variable {A B : ℝ}

/-- the window `4/5 < p < 4/3` around `p = 1` (it maps into `U` under `p ↦ 1/p - 1`) -/
noncomputable def V : Set ℝ := Ioo (4/5) (4/3)

lemma V_open : IsOpen V := isOpen_Ioo
lemma one_mem_V : (1:ℝ) ∈ V := by unfold V; constructor <;> norm_num
lemma V_pos {p : ℝ} (h : p ∈ V) : 0 < p := by have := h.1; linarith

lemma phi_mem {p : ℝ} (h : p ∈ V) : 1 / p - 1 ∈ U := by
  have hp := V_pos h
  unfold U; constructor
  · have : 3/4 < 1/p := by rw [lt_one_div (by norm_num) hp]; norm_num; linarith [h.2]
    linarith
  · have : 1/p < 5/4 := by rw [one_div_lt hp (by norm_num)]; norm_num; linarith [h.1]
    linarith

lemma hd_phi {p : ℝ} (hp : p ≠ 0) : HasDerivAt (fun x : ℝ => 1 / x - 1) (-1 / p ^ 2) p := by
  have := (hasDerivAt_inv hp).sub_const 1
  convert this using 1
  · funext x; rw [one_div]
  · rw [neg_div, one_div]

lemma hd_invpow (n : ℕ) {p : ℝ} (hp : p ≠ 0) :
    HasDerivAt (fun x : ℝ => 1 / x ^ (n+1)) (-((n:ℝ) + 1) / p ^ (n+2)) p := by
  have h := (hasDerivAt_pow (n+1) p).inv (pow_ne_zero _ hp)
  simp only [Nat.add_sub_cancel] at h
  convert h using 1
  · funext x; simp only [one_div, Pi.inv_apply]
  · field_simp; push_cast; ring

noncomputable def P0 (A B p : ℝ) : ℝ := 4 * r A B * I0 A B (1 / p - 1)
noncomputable def P1 (A B p : ℝ) : ℝ := 4 * r A B * (I1 A B (1 / p - 1) * (-1 / p ^ 2))
noncomputable def P2 (A B p : ℝ) : ℝ :=
  4 * r A B * (I2 A B (1 / p - 1) / p ^ 4 + I1 A B (1 / p - 1) * (2 / p ^ 3))
noncomputable def P3 (A B p : ℝ) : ℝ :=
  4 * r A B * (-(I3 A B (1 / p - 1)) / p ^ 6 - 6 * I2 A B (1 / p - 1) / p ^ 5 -
    6 * I1 A B (1 / p - 1) / p ^ 4)

lemma P_eqV (hA : 0 < A) {p : ℝ} (h : p ∈ V) : P A B p = P0 A B p := by
  have hp := V_pos h
  have hm := phi_mem h
  have h1 : -1 < 1 / p - 1 := by have := hm.1; linarith
  have := P_eq (A := A) (B := B) h1
  rw [show 1 / (1 + (1 / p - 1)) = p by field_simp; ring] at this
  have hr := r_pos (B := B) hA
  rw [this]; unfold P0 I0; field_simp

lemma hd_P (hA : 0 < A) (hB : 0 < B) {p : ℝ} (h : p ∈ V) :
    HasDerivAt (P0 A B) (P1 A B p) p ∧ HasDerivAt (P1 A B) (P2 A B p) p ∧
      HasDerivAt (P2 A B) (P3 A B p) p := by
  have hp := (V_pos h).ne'
  have hφ := hd_phi hp
  obtain ⟨g0, g1, g2⟩ := hd_I hA hB (phi_mem h)
  have k0 := g0.comp p hφ
  have k1 := g1.comp p hφ
  have k2 := g2.comp p hφ
  have q1 := hd_invpow 1 hp
  have q2 := hd_invpow 2 hp
  have q3 := hd_invpow 3 hp
  refine ⟨?_, ?_, ?_⟩
  · have hh := k0.const_mul (4 * r A B)
    unfold P0 P1; convert hh using 1
    rfl
  · have hh := (k1.fun_mul (q1.const_mul (-1))).const_mul (4 * r A B)
    unfold P1 P2; convert hh using 1
    · funext x; simp only [Function.comp]; ring
    · simp only [Function.comp]; field_simp; ring
  · have hh := ((k2.fun_mul q3).fun_add (k1.fun_mul (q2.const_mul 2))).const_mul (4 * r A B)
    unfold P2 P3; convert hh using 1
    · funext x; simp only [Function.comp]; ring
    · simp only [Function.comp]; field_simp; ring

lemma cont_P3 (hA : 0 < A) (hB : 0 < B) : ContinuousOn (P3 A B) V := by
  have hφc : ContinuousOn (fun p : ℝ => 1 / p - 1) V :=
    (continuousOn_const.div continuousOn_id (fun x hx => (V_pos hx).ne')).sub continuousOn_const
  have hmap : MapsTo (fun p : ℝ => 1 / p - 1) V U := fun p hp => phi_mem hp
  have hI1 : ContinuousOn (I1 A B) U := fun x hx =>
    (hd_I hA hB hx).2.1.continuousAt.continuousWithinAt
  have hI2 : ContinuousOn (I2 A B) U := fun x hx =>
    (hd_I hA hB hx).2.2.continuousAt.continuousWithinAt
  have c1 := hI1.comp hφc hmap
  have c2 := hI2.comp hφc hmap
  have c3 := (cont_I3 hA hB).comp hφc hmap
  have hpow : ∀ n : ℕ, ContinuousOn (fun p : ℝ => p ^ n) V := fun n => (continuous_pow n).continuousOn
  have hne : ∀ n : ℕ, ∀ p ∈ V, p ^ n ≠ 0 := fun n p hp => pow_ne_zero n (V_pos hp).ne'
  have := continuousOn_const (c := 4 * r A B) |>.mul
    ((((c3.neg).div (hpow 6) (hne 6)).sub ((continuousOn_const (c := (6:ℝ)).mul c2).div (hpow 5) (hne 5))).sub
      ((continuousOn_const (c := (6:ℝ)).mul c1).div (hpow 4) (hne 4)))
  exact this

lemma chainP (hA : 0 < A) (hB : 0 < B) :
    ContDiffOn ℝ 3 (P A B) V ∧ ∀ x ∈ V, deriv (P A B) x = P1 A B x ∧
      iteratedDeriv 2 (P A B) x = P2 A B x ∧ iteratedDeriv 3 (P A B) x = P3 A B x :=
  chain3 V_open
    (fun x hx => (hd_P hA hB hx).1.congr_of_eventuallyEq
      (by filter_upwards [V_open.mem_nhds hx] with y hy using P_eqV hA hy))
    (fun x hx => (hd_P hA hB hx).2.1) (fun x hx => (hd_P hA hB hx).2.2) (cont_P3 hA hB)

lemma hd_P_real (hA : 0 < A) (hB : 0 < B) : HasDerivAt (P A B) (P1 A B 1) 1 :=
  (hd_P hA hB one_mem_V).1.congr_of_eventuallyEq
    (by filter_upwards [V_open.mem_nhds one_mem_V] with y hy using P_eqV hA hy)

lemma h0U : (0:ℝ) ∈ U := by unfold U; constructor <;> norm_num

end Pfun

end X01

theorem X01.Ivals {A B : ℝ} (hA : 0 < A) (hB : 0 < B) :
    X01.I1 A B 0 = 0 ∧ X01.I2 A B 0 = 2 * ((1 - X01.q A B ^ 2) * Real.pi ^ 2 / 24) ∧
      X01.I3 A B 0 = 6 * (-((1 - X01.q A B ^ 2) * X01.z3 / 4)) := by
  obtain ⟨_, _, h1, h2, h3⟩ := X01.proof_10 A B hA hB
  obtain ⟨d1, d2, d3⟩ := (X01.chainI hA hB).2 0 X01.h0U
  exact ⟨d1.symm.trans h1, d2.symm.trans h2, d3.symm.trans h3⟩

theorem X01.proof_11 : X01.stmt_11 := fun _ _ hA hB =>
  (X01.chainP hA hB).1.contDiffAt (X01.V_open.mem_nhds X01.one_mem_V)

theorem X01.proof_12 : X01.stmt_12 := by
  intro A B _ _
  simp [X01.P, X01.integrand, X01.r]

theorem X01.proof_13 : X01.stmt_13 := by
  intro A B hA hB
  have h := X01.hd_P_real hA hB
  have e : X01.P1 A B 1 = 0 := by
    unfold X01.P1; norm_num; right; exact (X01.Ivals hA hB).1
  rwa [e] at h

theorem X01.proof_14 : X01.stmt_14 := by
  intro A B hA hB
  rw [((X01.chainP hA hB).2 1 X01.one_mem_V).2.1]
  obtain ⟨v1, v2, _⟩ := X01.Ivals hA hB
  unfold X01.P2; norm_num; rw [v1, v2]; unfold X01.q
  have hr := X01.r_pos (B := B) hA
  have hr2 := X01.r_sq A B
  have hB2 : B ^ 2 = X01.r A B ^ 2 - A ^ 2 := by linarith
  rw [← hr2, hB2]
  field_simp; ring

theorem X01.proof_15 : X01.stmt_15 := by
  intro A B hA hB
  rw [((X01.chainP hA hB).2 1 X01.one_mem_V).2.2]
  obtain ⟨v1, v2, v3⟩ := X01.Ivals hA hB
  unfold X01.P3; norm_num; rw [v1, v2, v3]; unfold X01.q
  have hr := X01.r_pos (B := B) hA
  have hr2 := X01.r_sq A B
  have hB2 : B ^ 2 = X01.r A B ^ 2 - A ^ 2 := by linarith
  rw [← hr2, hB2]
  field_simp; ring

theorem X01.proof_17 : X01.stmt_17 := fun A B hA hB =>
  ⟨X01.proof_11 A B hA hB, X01.proof_12 A B hA hB, (X01.proof_13 A B hA hB).deriv,
    X01.proof_14 A B hA hB, X01.proof_15 A B hA hB⟩

theorem X01.proof_16 : X01.stmt_16 := by
  obtain ⟨h11, h12, h13, h14, h15⟩ := X01.proof_17 1 1 one_pos one_pos
  refine ⟨h11, ?_⟩
  have hr : X01.r 1 1 = Real.sqrt 2 := by unfold X01.r; norm_num
  have h2 : Real.sqrt 2 ^ 2 = 2 := Real.sq_sqrt (by norm_num)
  have hs : 0 < Real.sqrt 2 := by positivity
  rw [hr] at h12 h14 h15
  simp only [iteratedDeriv_zero, iteratedDeriv_one, h12, h13, h14, h15, Nat.factorial]
  refine ⟨?_, ?_, ?_, ?_⟩
  · field_simp; norm_num
  · simp
  · field_simp; ring_nf; rw [show Real.sqrt 2 ^ 4 = (Real.sqrt 2 ^ 2) ^ 2 by ring, h2]; ring
  · field_simp; ring_nf; rw [show Real.sqrt 2 ^ 4 = (Real.sqrt 2 ^ 2) ^ 2 by ring, h2]; ring

theorem X01.proof_4 : X01.stmt_4 := by
  intro A B hA hB p hp
  set η := 1 / p - 1 with hη
  have hr := X01.r_pos (B := B) hA
  have hr2 := X01.r_sq A B
  have hR : A ^ 2 + B ^ 2 ≠ 0 := by positivity
  have key : ∀ t ∈ Set.uIoo (0:ℝ) 1, X01.integrand A B p t = X01.r A B *
      ((t * (1 - t)) ^ (η / 2) *
        Real.sqrt (Real.cosh (η * X01.lg t) + X01.q A B * Real.sinh (η * X01.lg t))) := by
    rw [X01.uIoo01]; intro t ht
    have h0 := ht.1
    have h1 : 0 < 1 - t := by linarith [ht.2]
    have hlog : Real.log (t * (1 - t)) = Real.log t + Real.log (1 - t) := X01.Lg_eq t
    unfold X01.integrand
    rw [Real.rpow_def_of_pos h0, Real.rpow_def_of_pos h1, Real.rpow_def_of_pos (mul_pos h0 h1),
      hlog, X01.lg_eq]
    set a := Real.log t
    set b := Real.log (1 - t)
    have E : 0 ≤ X01.r A B * Real.exp ((a + b) * (η / 2)) := by positivity
    rw [← mul_assoc, ← Real.sqrt_sq E, ← Real.sqrt_mul (sq_nonneg _)]
    congr 1
    rw [Real.cosh_eq, Real.sinh_eq, mul_pow, hr2]
    unfold X01.q
    have e1 : Real.exp (a * (2 / p - 2)) = Real.exp (η * a) ^ 2 := by
      rw [← Real.exp_nat_mul]; congr 1; rw [hη]; push_cast; ring
    have e2 : Real.exp (b * (2 / p - 2)) = Real.exp (η * b) ^ 2 := by
      rw [← Real.exp_nat_mul]; congr 1; rw [hη]; push_cast; ring
    have e3 : Real.exp ((a + b) * (η / 2)) ^ 2 = Real.exp (η * a) * Real.exp (η * b) := by
      rw [← Real.exp_nat_mul, ← Real.exp_add]; congr 1; push_cast; ring
    have e4 : Real.exp (η * (a - b)) = Real.exp (η * a) / Real.exp (η * b) := by
      rw [← Real.exp_sub]; congr 1; ring
    have e5 : Real.exp (-(η * (a - b))) = Real.exp (η * b) / Real.exp (η * a) := by
      rw [← Real.exp_sub]; congr 1; ring
    rw [e1, e2, e3, e4, e5]
    have pa := Real.exp_pos (η * a)
    have pb := Real.exp_pos (η * b)
    field_simp
    ring
  unfold X01.P
  rw [intervalIntegral.integral_congr_uIoo key, intervalIntegral.integral_const_mul]
  rw [hη]; field_simp; ring

namespace X01
open MeasureTheory Set Filter Real Interval Topology

lemma ii_beta {c1 c2 : ℝ} (h1 : -1 < c1) (h2 : -1 < c2) :
    IntervalIntegrable (fun t : ℝ => t ^ c1 * (1 - t) ^ c2) volume 0 1 := by
  have hl : IntervalIntegrable (fun t : ℝ => t ^ c1 * (1 - t) ^ c2) volume 0 (1/2) := by
    apply IntervalIntegrable.mul_continuousOn (intervalIntegral.intervalIntegrable_rpow' h1)
    apply ContinuousOn.rpow_const (by fun_prop)
    intro x hx; left; rw [uIcc_of_le (by norm_num)] at hx; linarith [hx.2]
  have hr : IntervalIntegrable (fun t : ℝ => t ^ c1 * (1 - t) ^ c2) volume (1/2) 1 := by
    have i2 : IntervalIntegrable (fun t : ℝ => (1 - t) ^ c2) volume (1/2) 1 := by
      have := (intervalIntegral.intervalIntegrable_rpow' h2 (a := 0) (b := 1/2)).comp_sub_left 1
      norm_num at this; exact this.symm
    apply IntervalIntegrable.continuousOn_mul i2
    apply ContinuousOn.rpow_const continuousOn_id
    intro x hx; left; rw [uIcc_of_le (by norm_num)] at hx; simp only [id]; linarith [hx.1]
  exact hl.trans hr

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

end X01

theorem X01.proof_5 : X01.stmt_5 := by
  intro η s hη hs
  have hs' := abs_lt.mp hs
  have ha : 0 < 1 + η / 2 + s := by linarith
  have hb : 0 < 1 + η / 2 - s := by linarith
  have key : Set.EqOn (fun t : ℝ => (t * (1 - t)) ^ (η / 2) * Real.exp (s * X01.lg t))
      (fun t : ℝ => t ^ (1 + η / 2 + s - 1) * (1 - t) ^ (1 + η / 2 - s - 1)) (Set.uIoo 0 1) := by
    rw [X01.uIoo01]; intro t ht
    have h0 := ht.1
    have h1 : 0 < 1 - t := by linarith [ht.2]
    simp only
    rw [Real.rpow_def_of_pos h0, Real.rpow_def_of_pos h1, Real.rpow_def_of_pos (mul_pos h0 h1),
      show Real.log (t * (1 - t)) = Real.log t + Real.log (1 - t) from X01.Lg_eq t, X01.lg_eq,
      ← Real.exp_add, ← Real.exp_add]
    congr 1
    ring
  constructor
  · exact (X01.ii_beta (by linarith) (by linarith)).congr_uIoo key.symm
  · rw [intervalIntegral.integral_congr_uIoo key, X01.beta_real ha hb]
    have h1 : (0:ℝ) < 1 + η := by linarith
    rw [show 1 + η / 2 + s + (1 + η / 2 - s) = (1 + η) + 1 by ring, Real.Gamma_add_one h1.ne']
    have hg := Real.Gamma_pos_of_pos h1
    field_simp

theorem X01.proof_1 : X01.stmt_1 := by
  intro A B hA hB p hp
  have hc : -1 < 1 / p - 1 := by have : 0 < 1 / p := by positivity
                                 linarith
  have i1 : IntervalIntegrable (fun t : ℝ => t ^ (1 / p - 1)) MeasureTheory.volume 0 1 :=
    intervalIntegral.intervalIntegrable_rpow' hc
  have i2 : IntervalIntegrable (fun t : ℝ => (1 - t) ^ (1 / p - 1)) MeasureTheory.volume 0 1 := by
    refine (X01.ii_beta (c1 := 0) (by norm_num) hc).congr_uIoo ?_
    intro t _; simp
  refine IntervalIntegrable.mono_fun' ((i1.const_mul A).add (i2.const_mul B))
    (by unfold X01.integrand; fun_prop : Measurable (X01.integrand A B p)).aestronglyMeasurable ?_
  filter_upwards [X01.ae_Ioo01] with t ht
  have h0 := ht.1
  have h1 : 0 < 1 - t := by linarith [ht.2]
  simp only [Real.norm_eq_abs, X01.integrand]
  rw [abs_of_nonneg (Real.sqrt_nonneg _)]
  have ea : t ^ (2 / p - 2) = (t ^ (1 / p - 1)) ^ 2 := by
    rw [← Real.rpow_natCast, ← Real.rpow_mul h0.le]; congr 1; push_cast; ring
  have eb : (1 - t) ^ (2 / p - 2) = ((1 - t) ^ (1 / p - 1)) ^ 2 := by
    rw [← Real.rpow_natCast, ← Real.rpow_mul h1.le]; congr 1; push_cast; ring
  rw [ea, eb]
  have u := Real.rpow_pos_of_pos h0 (1 / p - 1)
  have v := Real.rpow_pos_of_pos h1 (1 / p - 1)
  calc Real.sqrt (A ^ 2 * (t ^ (1 / p - 1)) ^ 2 + B ^ 2 * ((1 - t) ^ (1 / p - 1)) ^ 2)
      ≤ Real.sqrt ((A * t ^ (1 / p - 1) + B * (1 - t) ^ (1 / p - 1)) ^ 2) :=
        Real.sqrt_le_sqrt (by nlinarith [mul_pos (mul_pos hA hB) (mul_pos u v)])
    _ = A * t ^ (1 / p - 1) + B * (1 - t) ^ (1 / p - 1) := Real.sqrt_sq (by positivity)

namespace X01
open MeasureTheory Set Filter Real Interval Topology

section loopsec
variable {A B p : ℝ}

lemma pw1 (hp : 0 < p) {x : ℝ} (hx : 0 ≤ x) : (x ^ (1 / p)) ^ p = x := by
  rw [← Real.rpow_mul hx, one_div_mul_cancel hp.ne', Real.rpow_one]

lemma pw2 (hp : 0 < p) {x : ℝ} (hx : 0 ≤ x) : (x ^ p) ^ (1 / p) = x := by
  rw [← Real.rpow_mul hx, mul_one_div_cancel hp.ne', Real.rpow_one]

@[simp] lemma vec0 (a b : ℝ) : (!₂[a, b] : Plane) 0 = a := rfl
@[simp] lemma vec1 (a b : ℝ) : (!₂[a, b] : Plane) 1 = b := rfl

lemma plane_ext {v w : Plane} (h0 : v 0 = w 0) (h1 : v 1 = w 1) : v = w := by
  ext i; fin_cases i <;> simp [h0, h1]

lemma plane_eta (v : Plane) : (!₂[v 0, v 1] : Plane) = v := plane_ext rfl rfl

lemma norm_vec (a b : ℝ) : ‖(!₂[a, b] : Plane)‖ = Real.sqrt (a ^ 2 + b ^ 2) := by
  rw [EuclideanSpace.norm_eq, Fin.sum_univ_two]; simp

lemma dist_vec (a b c d : ℝ) :
    dist (!₂[a, b] : Plane) !₂[c, d] = Real.sqrt ((a - c) ^ 2 + (b - d) ^ 2) := by
  rw [EuclideanSpace.dist_eq, Fin.sum_univ_two]; simp [Real.dist_eq, sq_abs]

/-- Embedding of `ℝ × ℝ` into the plane. -/
noncomputable def emb (z : ℝ × ℝ) : Plane := !₂[z.1, z.2]

lemma emb_inj : Function.Injective emb := by
  intro z w h
  have h0 := congrArg (fun v : Plane => v 0) h
  have h1 := congrArg (fun v : Plane => v 1) h
  simp [emb] at h0 h1
  exact Prod.ext h0 h1

lemma emb_cont : Continuous emb := by
  unfold emb
  exact (PiLp.continuousLinearEquiv 2 ℝ (fun _ : Fin 2 => ℝ)).symm.continuous.comp
    (continuous_pi fun i => by fin_cases i <;> simp <;> fun_prop)

lemma edist_emb (z w : ℝ × ℝ) :
    edist (emb z) (emb w) = ENNReal.ofReal (Real.sqrt ((z.1 - w.1) ^ 2 + (z.2 - w.2) ^ 2)) := by
  rw [edist_dist, emb, emb, dist_vec]

noncomputable def loopP (A B p t : ℝ) : ℝ × ℝ :=
  if t ≤ 1 then (A * t ^ (1 / p), B * (1 - t) ^ (1 / p))
  else if t ≤ 2 then (A * (2 - t) ^ (1 / p), -(B * (t - 1) ^ (1 / p)))
  else if t ≤ 3 then (-(A * (t - 2) ^ (1 / p)), -(B * (3 - t) ^ (1 / p)))
  else (-(A * (4 - t) ^ (1 / p)), B * (t - 3) ^ (1 / p))

lemma loop_eq (A B p t : ℝ) : loop A B p t = emb (loopP A B p t) := by
  unfold loop loopP emb; split_ifs <;> rfl

/-- the quarter arc as a map into `ℝ × ℝ` -/
noncomputable def arc (A B p s : ℝ) : ℝ × ℝ := (A * s ^ (1 / p), B * (1 - s) ^ (1 / p))

lemma loop_Q1 (t : ℝ) (ht : t ∈ Icc (0:ℝ) 1) : loopP A B p t = arc A B p t := by
  simp [loopP, arc, ht.2]

lemma loop_Q2 (hp : 0 < p) (t : ℝ) (ht : t ∈ Icc (1:ℝ) 2) :
    loopP A B p t = ((arc A B p (2 - t)).1, -(arc A B p (2 - t)).2) := by
  have hr : (1 / p) ≠ 0 := by positivity
  have h0 : (0:ℝ) ^ p⁻¹ = 0 := by simpa using Real.zero_rpow hr
  rcases ht.1.eq_or_lt with rfl | h1
  · norm_num [loopP, arc, h0]
  · have h : ¬ t ≤ 1 := not_le.2 h1
    simp only [loopP, arc, h, ite_false, ht.2, ite_true]
    rw [show 1 - (2 - t) = t - 1 by ring]

lemma loop_Q3 (hp : 0 < p) (t : ℝ) (ht : t ∈ Icc (2:ℝ) 3) :
    loopP A B p t = (-(arc A B p (t - 2)).1, -(arc A B p (t - 2)).2) := by
  have hr : (1 / p) ≠ 0 := by positivity
  have h0 : (0:ℝ) ^ p⁻¹ = 0 := by simpa using Real.zero_rpow hr
  have h : ¬ t ≤ 1 := by have := ht.1; intro; linarith
  rcases ht.1.eq_or_lt with rfl | h1
  · norm_num [loopP, arc, h0]
  · have h' : ¬ t ≤ 2 := not_le.2 h1
    simp only [loopP, arc, h, h', ite_false, ht.2, ite_true]
    rw [show 1 - (t - 2) = 3 - t by ring]

lemma loop_Q4 (hp : 0 < p) (t : ℝ) (ht : t ∈ Icc (3:ℝ) 4) :
    loopP A B p t = (-(arc A B p (4 - t)).1, (arc A B p (4 - t)).2) := by
  have hr : (1 / p) ≠ 0 := by positivity
  have h0 : (0:ℝ) ^ p⁻¹ = 0 := by simpa using Real.zero_rpow hr
  have h : ¬ t ≤ 1 := by have := ht.1; intro; linarith
  have h' : ¬ t ≤ 2 := by have := ht.1; intro; linarith
  rcases ht.1.eq_or_lt with rfl | h1
  · norm_num [loopP, arc, h0]
  · have h'' : ¬ t ≤ 3 := not_le.2 h1
    simp only [loopP, arc, h, h', h'', ite_false]
    rw [show 1 - (4 - t) = t - 3 by ring]

lemma arc_cont (hp : 0 < p) : Continuous (arc A B p) := by
  have hr : 0 ≤ 1 / p := by positivity
  exact (continuous_const.mul ((Real.continuous_rpow_const hr).comp continuous_id)).prodMk
    (continuous_const.mul ((Real.continuous_rpow_const hr).comp (continuous_const.sub continuous_id)))

lemma loopP_cont (hp : 0 < p) : ContinuousOn (loopP A B p) (Icc 0 4) := by
  have hc := arc_cont (A := A) (B := B) hp
  have c1 : ContinuousOn (loopP A B p) (Icc 0 1) :=
    hc.continuousOn.congr (fun t ht => loop_Q1 t ht)
  have c2 : ContinuousOn (loopP A B p) (Icc 1 2) :=
    ((hc.comp (continuous_const.sub continuous_id)).fst.prodMk
      (hc.comp (continuous_const.sub continuous_id)).snd.neg |>.continuousOn).congr
      (fun t ht => loop_Q2 hp t ht)
  have c3 : ContinuousOn (loopP A B p) (Icc 2 3) :=
    ((hc.comp (continuous_id.sub continuous_const)).fst.neg.prodMk
      (hc.comp (continuous_id.sub continuous_const)).snd.neg |>.continuousOn).congr
      (fun t ht => loop_Q3 hp t ht)
  have c4 : ContinuousOn (loopP A B p) (Icc 3 4) :=
    ((hc.comp (continuous_const.sub continuous_id)).fst.neg.prodMk
      (hc.comp (continuous_const.sub continuous_id)).snd |>.continuousOn).congr
      (fun t ht => loop_Q4 hp t ht)
  rw [← Icc_union_Icc_eq_Icc (a := 0) (b := 1) (c := 4) (by norm_num) (by norm_num),
    ← Icc_union_Icc_eq_Icc (a := 1) (b := 2) (c := 4) (by norm_num) (by norm_num),
    ← Icc_union_Icc_eq_Icc (a := 2) (b := 3) (c := 4) (by norm_num) (by norm_num)]
  exact c1.union_of_isClosed (c2.union_of_isClosed (c3.union_of_isClosed c4 isClosed_Icc isClosed_Icc)
    isClosed_Icc (isClosed_Icc.union isClosed_Icc)) isClosed_Icc
    (isClosed_Icc.union (isClosed_Icc.union isClosed_Icc))

noncomputable def ginv (A B p : ℝ) (z : ℝ × ℝ) : ℝ :=
  if 0 ≤ z.1 ∧ 0 < z.2 then |z.1 / A| ^ p
  else if 0 < z.1 ∧ z.2 ≤ 0 then 1 + |z.2 / B| ^ p
  else if z.1 ≤ 0 ∧ z.2 < 0 then 2 + |z.1 / A| ^ p
  else 3 + |z.2 / B| ^ p

lemma ginv_loop (hA : 0 < A) (hB : 0 < B) (hp : 0 < p) (t : ℝ) (ht : t ∈ Ico (0:ℝ) 4) :
    ginv A B p (loopP A B p t) = t := by
  have nn : ∀ x : ℝ, 0 ≤ x → 0 ≤ x ^ (1 / p) := fun x hx => Real.rpow_nonneg hx _
  have ps : ∀ x : ℝ, 0 < x → 0 < x ^ (1 / p) := fun x hx => Real.rpow_pos_of_pos hx _
  have dA : ∀ x : ℝ, A * x / A = x := fun x => by field_simp
  have dB : ∀ x : ℝ, B * x / B = x := fun x => by field_simp
  have ht2 := ht.2
  by_cases h1 : t < 1
  · rw [loop_Q1 (A := A) (B := B) (p := p) t ⟨ht.1, h1.le⟩]
    simp only [arc]
    have a1 := nn t ht.1
    have a2 := ps (1 - t) (by linarith)
    rw [ginv, if_pos ⟨by positivity, by positivity⟩, dA, abs_of_nonneg a1, pw1 hp ht.1]
  by_cases h2 : t < 2
  · rw [loop_Q2 hp t ⟨not_lt.1 h1, h2.le⟩]
    simp only [arc]
    have a1 := nn (1 - (2 - t)) (by linarith)
    have a2 := ps (2 - t) (by linarith)
    rw [ginv, if_neg (by intro h; have := h.2; nlinarith [mul_pos hA a2]),
      if_pos ⟨by positivity, by simp only [neg_nonpos]; positivity⟩, neg_div, dB, abs_neg,
      abs_of_nonneg a1, pw1 hp (by linarith)]
    ring
  by_cases h3 : t < 3
  · rw [loop_Q3 hp t ⟨not_lt.1 h2, h3.le⟩]
    simp only [arc]
    have a1 := nn (t - 2) (by linarith)
    have a2 := ps (1 - (t - 2)) (by linarith)
    have b1 : 0 < B * (1 - (t - 2)) ^ (1 / p) := by positivity
    have b0 : 0 ≤ A * (t - 2) ^ (1 / p) := by positivity
    rw [ginv, if_neg (by intro h; linarith [h.2]), if_neg (by intro h; linarith [h.1]),
      if_pos ⟨by linarith, by linarith⟩, neg_div, dA, abs_neg, abs_of_nonneg a1,
      pw1 hp (by linarith)]
    ring
  · rw [loop_Q4 hp t ⟨not_lt.1 h3, ht.2.le⟩]
    simp only [arc]
    have a1 := nn (1 - (4 - t)) (by linarith)
    have a2 := ps (4 - t) (by linarith)
    have b1 : 0 < A * (4 - t) ^ (1 / p) := by positivity
    have b0 : 0 ≤ B * (1 - (4 - t)) ^ (1 / p) := by positivity
    rw [ginv, if_neg (by intro h; linarith [h.1]), if_neg (by intro h; linarith [h.1]),
      if_neg (by intro h; linarith [h.2]), dB, abs_of_nonneg a1, pw1 hp (by linarith)]
    ring

lemma norm_loop (hA : 0 < A) (hB : 0 < B) (hp : 0 < p) (t : ℝ) (ht : t ∈ Icc (0:ℝ) 4) :
    |(loopP A B p t).1 / A| ^ p + |(loopP A B p t).2 / B| ^ p = 1 := by
  have nn : ∀ x : ℝ, 0 ≤ x → 0 ≤ x ^ (1 / p) := fun x hx => Real.rpow_nonneg hx _
  have dA : ∀ x : ℝ, A * x / A = x := fun x => by field_simp
  have dB : ∀ x : ℝ, B * x / B = x := fun x => by field_simp
  have key : ∀ s : ℝ, 0 ≤ s → s ≤ 1 →
      |(arc A B p s).1 / A| ^ p + |(arc A B p s).2 / B| ^ p = 1 := by
    intro s h0 h1
    simp only [arc, dA, dB]
    rw [abs_of_nonneg (nn _ h0), abs_of_nonneg (nn _ (by linarith)), pw1 hp h0,
      pw1 hp (by linarith)]
    ring
  by_cases h1 : t ≤ 1
  · rw [loop_Q1 t ⟨ht.1, h1⟩]; exact key _ (by linarith [ht.1]) (by linarith)
  by_cases h2 : t ≤ 2
  · rw [loop_Q2 hp t ⟨(not_le.1 h1).le, h2⟩]
    have := key (2 - t) (by linarith) (by linarith [not_le.1 h1])
    simpa [neg_div, abs_neg] using this
  by_cases h3 : t ≤ 3
  · rw [loop_Q3 hp t ⟨(not_le.1 h2).le, h3⟩]
    have := key (t - 2) (by linarith [not_le.1 h2]) (by linarith)
    simpa [neg_div, abs_neg] using this
  · rw [loop_Q4 hp t ⟨(not_le.1 h3).le, ht.2⟩]
    have := key (4 - t) (by linarith [ht.2]) (by linarith [not_le.1 h3])
    simpa [neg_div, abs_neg] using this

lemma exists_t (hA : 0 < A) (hB : 0 < B) (hp : 0 < p) (z : ℝ × ℝ)
    (hz : |z.1 / A| ^ p + |z.2 / B| ^ p = 1) : ∃ t ∈ Icc (0:ℝ) 4, loopP A B p t = z := by
  set u := |z.1 / A| ^ p with hu
  set v := |z.2 / B| ^ p with hv
  have u0 : 0 ≤ u := Real.rpow_nonneg (abs_nonneg _) _
  have v0 : 0 ≤ v := Real.rpow_nonneg (abs_nonneg _) _
  have ur : A * u ^ (1 / p) = |z.1| := by
    rw [hu, pw2 hp (abs_nonneg _), abs_div, abs_of_pos hA]; field_simp
  have vr : B * v ^ (1 / p) = |z.2| := by
    rw [hv, pw2 hp (abs_nonneg _), abs_div, abs_of_pos hB]; field_simp
  have u1 : (1 - u) = v := by linarith
  have v1 : (1 - v) = u := by linarith
  by_cases c1 : 0 ≤ z.1 ∧ 0 < z.2
  · refine ⟨u, ⟨u0, by linarith⟩, ?_⟩
    rw [loop_Q1 u ⟨u0, by linarith⟩]
    simp only [arc, u1, ur, vr]
    ext <;> simp [abs_of_nonneg c1.1, abs_of_pos c1.2]
  by_cases c2 : 0 < z.1 ∧ z.2 ≤ 0
  · refine ⟨1 + v, ⟨by linarith, by linarith⟩, ?_⟩
    rw [loop_Q2 hp (1 + v) ⟨by linarith, by linarith⟩]
    simp only [arc, show 2 - (1 + v) = u by linarith, u1, ur, vr]
    ext <;> simp [abs_of_pos c2.1, abs_of_nonpos c2.2]
  by_cases c3 : z.1 ≤ 0 ∧ z.2 < 0
  · refine ⟨2 + u, ⟨by linarith, by linarith⟩, ?_⟩
    rw [loop_Q3 hp (2 + u) ⟨by linarith, by linarith⟩]
    simp only [arc, add_sub_cancel_left, u1, ur, vr]
    ext <;> simp [abs_of_nonpos c3.1, abs_of_neg c3.2]
  · have h1 : z.1 ≤ 0 := by
      by_contra h
      have hx := not_le.1 h
      have : 0 < z.2 := by
        by_contra h'
        exact c2 ⟨hx, not_lt.1 h'⟩
      exact c1 ⟨hx.le, this⟩
    have h2 : 0 ≤ z.2 := by
      by_contra h
      exact c3 ⟨h1, not_le.1 h⟩
    refine ⟨3 + v, ⟨by linarith, by linarith⟩, ?_⟩
    rw [loop_Q4 hp (3 + v) ⟨by linarith, by linarith⟩]
    simp only [arc, show 4 - (3 + v) = u by linarith, u1, ur, vr]
    ext <;> simp [abs_of_nonpos h1, abs_of_nonneg h2]

end loopsec

end X01

theorem X01.proof_2 : X01.stmt_2 := by
  intro A B hA hB p hp
  have hr : (1 / p) ≠ 0 := by positivity
  have h0 : (0:ℝ) ^ p⁻¹ = 0 := by simpa using Real.zero_rpow hr
  have hl : X01.loop A B p = X01.emb ∘ X01.loopP A B p := funext fun t => X01.loop_eq A B p t
  refine ⟨?_, ?_, ?_, ?_⟩
  · rw [hl]; exact X01.emb_cont.comp_continuousOn (X01.loopP_cont hp)
  · rw [X01.loop_eq, X01.loop_eq]; norm_num [X01.loopP, h0]
  · rw [hl]
    exact X01.emb_inj.comp_injOn
      (Set.LeftInvOn.injOn (f₁' := X01.ginv A B p) (fun t ht => X01.ginv_loop hA hB hp t ht))
  · ext v
    constructor
    · rintro ⟨t, ht, rfl⟩
      have := X01.norm_loop hA hB hp t ht
      rw [X01.loop_eq]
      exact this
    · intro hv
      have hv' : |v 0 / A| ^ p + |v 1 / B| ^ p = 1 := hv
      obtain ⟨t, ht, he⟩ := X01.exists_t hA hB hp (v 0, v 1) hv'
      refine ⟨t, ht, ?_⟩
      rw [X01.loop_eq, he]
      exact X01.plane_eta v

namespace X01
open MeasureTheory Set Filter Real Interval Topology

section lengthsec
variable {A B p : ℝ}

/-- the quarter arc in the plane -/
noncomputable def qa (A B p : ℝ) (s : ℝ) : Plane := emb (arc A B p s)

lemma qa_cont (hp : 0 < p) : Continuous (qa A B p) := emb_cont.comp (arc_cont hp)

/-- the speed of the quarter arc -/
noncomputable def spd (A B p t : ℝ) : ℝ := 1 / p * integrand A B p t

lemma qa_deriv (hp : 0 < p) {t : ℝ} (ht : t ∈ Ioo (0 : ℝ) 1) :
    ∃ v : Plane, HasDerivAt (qa A B p) v t ∧ ‖v‖ = spd A B p t := by
  have ht1 : 0 < 1 - t := by linarith [ht.2]
  set a := A * (1 / p * t ^ (1 / p - 1))
  set b := -(B * (1 / p * (1 - t) ^ (1 / p - 1)))
  refine ⟨!₂[a, b], ?_, ?_⟩
  · have hA : HasDerivAt (fun s : ℝ => A * s ^ (1 / p)) a t :=
      (Real.hasDerivAt_rpow_const (Or.inl ht.1.ne')).const_mul A
    have hB : HasDerivAt (fun s : ℝ => B * (1 - s) ^ (1 / p)) b t := by
      have := (((hasDerivAt_id' t).const_sub 1).rpow_const (p := 1 / p) (Or.inl ht1.ne')).const_mul B
      convert this using 1
      simp only [b]
      ring
    have hpi : HasDerivAt (fun s : ℝ => (![A * s ^ (1 / p), B * (1 - s) ^ (1 / p)] : Fin 2 → ℝ))
        ![a, b] t := by
      rw [hasDerivAt_pi]
      intro i; fin_cases i
      · simpa using hA
      · simpa using hB
    exact (PiLp.continuousLinearEquiv 2 ℝ (fun _ : Fin 2 => ℝ)).symm.hasFDerivAt.comp_hasDerivAt t hpi
  · rw [norm_vec, spd, integrand]
    have e1 : a ^ 2 = (1 / p) ^ 2 * (A ^ 2 * t ^ (2 / p - 2)) := by
      simp only [a]; rw [mul_pow, mul_pow, ← Real.rpow_mul_natCast ht.1.le]
      rw [show (1 / p - 1) * ((2:ℕ):ℝ) = 2 / p - 2 by push_cast; ring]; ring
    have e2 : b ^ 2 = (1 / p) ^ 2 * (B ^ 2 * (1 - t) ^ (2 / p - 2)) := by
      simp only [b]; rw [neg_sq, mul_pow, mul_pow, ← Real.rpow_mul_natCast ht1.le]
      rw [show (1 / p - 1) * ((2:ℕ):ℝ) = 2 / p - 2 by push_cast; ring]; ring
    rw [e1, e2, ← mul_add, Real.sqrt_mul (by positivity), Real.sqrt_sq (by positivity)]

/-- Arc length from `0` to `x`, as an integral of the speed. -/
noncomputable def Gp (A B p x : ℝ) : ℝ := ∫ t in (0:ℝ)..x, spd A B p t

lemma spd_iInt (hA : 0 < A) (hB : 0 < B) (hp : 0 < p) :
    IntervalIntegrable (spd A B p) volume 0 1 :=
  (proof_1 A B hA hB p hp).const_mul _

lemma spd_nonneg (hp : 0 < p) (t : ℝ) : 0 ≤ spd A B p t := by
  unfold spd integrand
  have : 0 < 1 / p := by positivity
  positivity

lemma spd_sub (hA : 0 < A) (hB : 0 < B) (hp : 0 < p) {a b : ℝ} (ha : 0 ≤ a) (hab : a ≤ b)
    (hb : b ≤ 1) : IntervalIntegrable (spd A B p) volume a b := by
  refine (spd_iInt hA hB hp).mono_set ?_
  rw [uIcc_of_le zero_le_one, uIcc_of_le hab]
  exact Icc_subset_Icc ha hb

lemma Gp_cont (hA : 0 < A) (hB : 0 < B) (hp : 0 < p) : ContinuousOn (Gp A B p) (Icc 0 1) := by
  have := intervalIntegral.continuousOn_primitive_interval' (spd_iInt hA hB hp) left_mem_uIcc
  rw [uIcc_of_le zero_le_one] at this
  exact this

lemma spd_contAt {t : ℝ} (ht : t ∈ Ioo (0:ℝ) 1) : ContinuousAt (spd A B p) t := by
  have h1 : 0 < 1 - t := by linarith [ht.2]
  unfold spd integrand
  have c1 : ContinuousAt (fun s : ℝ => s ^ (2 / p - 2)) t :=
    Real.continuousAt_rpow_const _ _ (Or.inl ht.1.ne')
  have c2 : ContinuousAt (fun s : ℝ => (1 - s) ^ (2 / p - 2)) t :=
    (Real.continuousAt_rpow_const _ _ (Or.inl h1.ne')).comp (continuousAt_const.sub continuousAt_id)
  exact continuousAt_const.mul ((continuousAt_const.mul c1).add (continuousAt_const.mul c2)).sqrt

lemma Gp_deriv (hA : 0 < A) (hB : 0 < B) (hp : 0 < p) {t : ℝ} (ht : t ∈ Ioo (0:ℝ) 1) :
    HasDerivAt (Gp A B p) (spd A B p t) t := by
  have hm : Measurable (spd A B p) := by unfold spd integrand; fun_prop
  exact intervalIntegral.integral_hasDerivAt_right (spd_sub hA hB hp le_rfl ht.1.le ht.2.le)
    hm.stronglyMeasurable.stronglyMeasurableAtFilter (spd_contAt ht)

lemma Gp_sub (hA : 0 < A) (hB : 0 < B) (hp : 0 < p) {a b : ℝ} (ha : 0 ≤ a) (hab : a ≤ b)
    (hb : b ≤ 1) : Gp A B p b - Gp A B p a = ∫ t in a..b, spd A B p t := by
  unfold Gp
  rw [intervalIntegral.integral_interval_sub_left (spd_sub hA hB hp le_rfl (ha.trans hab) hb)
    (spd_sub hA hB hp le_rfl ha (hab.trans hb))]

lemma Gp_mono (hA : 0 < A) (hB : 0 < B) (hp : 0 < p) {a b : ℝ} (ha : 0 ≤ a) (hab : a ≤ b)
    (hb : b ≤ 1) : Gp A B p a ≤ Gp A B p b := by
  have := Gp_sub hA hB hp ha hab hb
  have h2 : 0 ≤ ∫ t in a..b, spd A B p t :=
    intervalIntegral.integral_nonneg hab fun u _ => spd_nonneg hp u
  linarith

lemma qa_hasDeriv (hp : 0 < p) {t : ℝ} (ht : t ∈ Ioo (0:ℝ) 1) :
    HasDerivAt (qa A B p) (deriv (qa A B p) t) t ∧ ‖deriv (qa A B p) t‖ = spd A B p t := by
  obtain ⟨v, hv, hn⟩ := qa_deriv (A := A) (B := B) hp ht
  rw [hv.deriv]; exact ⟨hv, hn⟩

lemma up_open (hA : 0 < A) (hB : 0 < B) (hp : 0 < p) {x y : ℝ} (hx : 0 < x) (hxy : x ≤ y)
    (hy : y < 1) : ‖qa A B p y - qa A B p x‖ ≤ Gp A B p y - Gp A B p x := by
  have hI : ∀ z ∈ Ico x y, z ∈ Ioo (0:ℝ) 1 := fun z hz => ⟨by linarith [hz.1], by linarith [hz.2]⟩
  refine image_norm_le_of_norm_deriv_right_le_deriv_boundary' (a := x) (b := y)
    (f := fun z => qa A B p z - qa A B p x) (f' := deriv (qa A B p))
    ((qa_cont hp).continuousOn.sub continuousOn_const)
    (fun z hz => ((qa_hasDeriv hp (hI z hz)).1.sub_const _).hasDerivWithinAt)
    (B := fun z => Gp A B p z - Gp A B p x) (B' := spd A B p) (by simp)
    (((Gp_cont hA hB hp).mono (Icc_subset_Icc hx.le hy.le)).sub continuousOn_const)
    (fun z hz => ((Gp_deriv hA hB hp (hI z hz)).sub_const _).hasDerivWithinAt)
    (fun z hz => (qa_hasDeriv hp (hI z hz)).2.le) ⟨hxy, le_rfl⟩

lemma eps_tend {x : ℝ} (hx : x ∈ Icc (0:ℝ) 1) :
    Tendsto (fun e : ℝ => e + (1 - 2 * e) * x) (𝓝[>] 0) (𝓝[Icc 0 1] x) := by
  apply tendsto_nhdsWithin_iff.2 ⟨?_, ?_⟩
  · have : Tendsto (fun e : ℝ => e + (1 - 2 * e) * x) (𝓝 0) (𝓝 (0 + (1 - 2 * 0) * x)) :=
      ((continuous_id.add ((continuous_const.sub (continuous_const.mul continuous_id)).mul
        continuous_const)).tendsto 0)
    simpa using this.mono_left nhdsWithin_le_nhds
  · filter_upwards [Ioo_mem_nhdsGT (show (0:ℝ) < 1/2 by norm_num)] with e he
    constructor <;> nlinarith [he.1, he.2, hx.1, hx.2]

lemma up_closed (hA : 0 < A) (hB : 0 < B) (hp : 0 < p) {x y : ℝ} (hx : 0 ≤ x) (hxy : x ≤ y)
    (hy : y ≤ 1) : ‖qa A B p y - qa A B p x‖ ≤ Gp A B p y - Gp A B p x := by
  have tx := eps_tend ⟨hx, hxy.trans hy⟩
  have ty := eps_tend ⟨hx.trans hxy, hy⟩
  have L : Tendsto (fun e : ℝ => ‖qa A B p (e + (1 - 2 * e) * y) - qa A B p (e + (1 - 2 * e) * x)‖)
      (𝓝[>] 0) (𝓝 ‖qa A B p y - qa A B p x‖) :=
    (((qa_cont hp).tendsto y).comp (ty.mono_right nhdsWithin_le_nhds)).sub
      (((qa_cont hp).tendsto x).comp (tx.mono_right nhdsWithin_le_nhds)) |>.norm
  have R : Tendsto (fun e : ℝ => Gp A B p (e + (1 - 2 * e) * y) - Gp A B p (e + (1 - 2 * e) * x))
      (𝓝[>] 0) (𝓝 (Gp A B p y - Gp A B p x)) :=
    ((Gp_cont hA hB hp y ⟨hx.trans hxy, hy⟩).tendsto.comp ty).sub
      ((Gp_cont hA hB hp x ⟨hx, hxy.trans hy⟩).tendsto.comp tx)
  refine le_of_tendsto_of_tendsto L R ?_
  filter_upwards [Ioo_mem_nhdsGT (show (0:ℝ) < 1/2 by norm_num)] with e he
  exact up_open hA hB hp (by nlinarith [he.1, he.2]) (by nlinarith [he.1, he.2])
    (by nlinarith [he.1, he.2])

lemma var_le {E F : Type*} [PseudoEMetricSpace E] [PseudoEMetricSpace F]
    (f : ℝ → E) (h : ℝ → F) (s : Set ℝ)
    (hd : ∀ a ∈ s, ∀ b ∈ s, edist (f a) (f b) ≤ edist (h a) (h b)) :
    eVariationOn f s ≤ eVariationOn h s := by
  unfold eVariationOn
  exact iSup_mono fun p => Finset.sum_le_sum fun i _ => hd _ (p.2.2.2 _) _ (p.2.2.2 _)

lemma var_congr {E F : Type*} [PseudoEMetricSpace E] [PseudoEMetricSpace F]
    (f : ℝ → E) (h : ℝ → F) (s : Set ℝ)
    (hd : ∀ a ∈ s, ∀ b ∈ s, edist (f a) (f b) = edist (h a) (h b)) :
    eVariationOn f s = eVariationOn h s :=
  le_antisymm (var_le f h s fun a ha b hb => (hd a ha b hb).le)
    (var_le h f s fun a ha b hb => (hd a ha b hb).ge)

lemma Gp_zero : Gp A B p 0 = 0 := intervalIntegral.integral_same

lemma var_up (hA : 0 < A) (hB : 0 < B) (hp : 0 < p) :
    eVariationOn (qa A B p) (Icc 0 1) ≤ ENNReal.ofReal (Gp A B p 1) := by
  have key : ∀ a ∈ Icc (0:ℝ) 1, ∀ b ∈ Icc (0:ℝ) 1, a ≤ b →
      edist (qa A B p a) (qa A B p b) ≤ edist (Gp A B p a) (Gp A B p b) := by
    intro a ha b hb hab
    rw [edist_dist, edist_dist, dist_eq_norm, norm_sub_rev, Real.dist_eq, abs_sub_comm,
      abs_of_nonneg (by linarith [Gp_mono hA hB hp ha.1 hab hb.2])]
    exact ENNReal.ofReal_le_ofReal (up_closed hA hB hp ha.1 hab hb.2)
  calc eVariationOn (qa A B p) (Icc 0 1) ≤ eVariationOn (Gp A B p) (Icc 0 1) := by
        refine var_le _ _ _ fun a ha b hb => ?_
        rcases le_total a b with h | h
        · exact key a ha b hb h
        · rw [edist_comm, edist_comm (Gp A B p a)]; exact key b hb a ha h
    _ = ENNReal.ofReal (Gp A B p 1) := by
        have hm : MonotoneOn (Gp A B p) (Icc 0 1) := fun a ha b hb hab =>
          Gp_mono hA hB hp ha.1 hab hb.2
        have := hm.eVariationOn_eq (a := 0) (b := 1) ⟨le_rfl, zero_le_one⟩ ⟨zero_le_one, le_rfl⟩
        simpa [Gp_zero] using this

lemma V_inc (hbv : BoundedVariationOn (qa A B p) (Icc 0 1)) {x z : ℝ}
    (hx : 0 ≤ x) (hxz : x ≤ z) (hz : z ≤ 1) :
    ‖qa A B p z - qa A B p x‖ ≤
      variationOnFromTo (qa A B p) (Icc 0 1) 0 z - variationOnFromTo (qa A B p) (Icc 0 1) 0 x := by
  rw [variationOnFromTo.sub_right hbv.locallyBoundedVariationOn ⟨le_rfl, zero_le_one⟩
    ⟨hx.trans hxz, hz⟩ ⟨hx, hxz.trans hz⟩, variationOnFromTo.eq_of_le _ _ hxz, ← dist_eq_norm,
    dist_comm, dist_edist]
  refine ENNReal.toReal_mono (hbv.mono inter_subset_left) ?_
  exact eVariationOn.edist_le _ ⟨⟨hx, hxz.trans hz⟩, le_rfl, hxz⟩ ⟨⟨hx.trans hxz, hz⟩, hxz, le_rfl⟩

lemma V_cont (hp : 0 < p) (hbv : BoundedVariationOn (qa A B p) (Icc 0 1)) :
    ContinuousOn (variationOnFromTo (qa A B p) (Icc 0 1) 0) (Icc 0 1) := fun x hx =>
  (hbv.continuousWithinAt_variationOnFromTo_iff ⟨le_rfl, zero_le_one⟩ hx).2
    (qa_cont hp).continuousWithinAt

lemma low_open (hA : 0 < A) (hB : 0 < B) (hp : 0 < p)
    (hbv : BoundedVariationOn (qa A B p) (Icc 0 1))
    {a b : ℝ} (ha : 0 < a) (hab : a ≤ b) (hb : b < 1) :
    Gp A B p b - Gp A B p a ≤
      variationOnFromTo (qa A B p) (Icc 0 1) 0 b - variationOnFromTo (qa A B p) (Icc 0 1) 0 a := by
  set V := variationOnFromTo (qa A B p) (Icc 0 1) 0 with hV
  set f : ℝ → ℝ := fun z => (Gp A B p z - Gp A B p a) - (V z - V a) with hf
  have hsub : Icc a b ⊆ Icc 0 1 := Icc_subset_Icc ha.le hb.le
  have hfc : ContinuousOn f (Icc a b) :=
    (((Gp_cont hA hB hp).mono hsub).sub continuousOn_const).sub
      (((V_cont hp hbv).mono hsub).sub continuousOn_const)
  have key := image_le_of_liminf_slope_right_le_deriv_boundary (a := a) (b := b) hfc
    (B := fun _ => 0) (B' := fun _ => 0) (by simp [f]) continuousOn_const
    (fun _ _ => hasDerivWithinAt_const _ _ _) ?_ ⟨hab, le_rfl⟩
  · simp only [f] at key; linarith
  intro x hx r hr
  have hx01 : x ∈ Ioo (0:ℝ) 1 := ⟨by linarith [hx.1], by linarith [hx.2]⟩
  obtain ⟨hd, hn⟩ := qa_hasDeriv (A := A) (B := B) hp hx01
  have t1 : Tendsto (slope (qa A B p) x) (𝓝[>] x) (𝓝 (deriv (qa A B p) x)) :=
    (hasDerivAt_iff_tendsto_slope.1 hd).mono_left (nhdsWithin_mono _ fun z hz => ne_of_gt hz)
  have t2 : Tendsto (slope (Gp A B p) x) (𝓝[>] x) (𝓝 (spd A B p x)) :=
    (hasDerivAt_iff_tendsto_slope.1 (Gp_deriv hA hB hp hx01)).mono_left
      (nhdsWithin_mono _ fun z hz => ne_of_gt hz)
  have t3 : Tendsto (fun z => slope (Gp A B p) x z - ‖slope (qa A B p) x z‖) (𝓝[>] x)
      (𝓝 (spd A B p x - ‖deriv (qa A B p) x‖)) := t2.sub t1.norm
  rw [hn, sub_self] at t3
  refine Filter.Eventually.frequently ?_
  filter_upwards [t3 (Iio_mem_nhds hr), self_mem_nhdsWithin,
    nhdsWithin_le_nhds (Iio_mem_nhds hx01.2)] with z hz hzx hz1
  simp only [mem_preimage, mem_Iio] at hz
  have hzx' : x < z := hzx
  have hpos : 0 < z - x := by linarith
  have hVi := V_inc hbv hx01.1.le hzx'.le (le_of_lt hz1)
  refine lt_of_le_of_lt ?_ hz
  rw [slope_def_field, slope_def_field, slope_def_module, norm_smul, Real.norm_eq_abs,
    abs_of_pos (inv_pos.2 hpos)]
  simp only [f]
  rw [div_le_iff₀ hpos]
  have : (z - x)⁻¹ * ‖qa A B p z - qa A B p x‖ * (z - x) = ‖qa A B p z - qa A B p x‖ := by
    field_simp
  rw [sub_mul, div_mul_cancel₀ _ hpos.ne', this]
  linarith

lemma var_eq (hA : 0 < A) (hB : 0 < B) (hp : 0 < p) :
    eVariationOn (qa A B p) (Icc 0 1) = ENNReal.ofReal (Gp A B p 1) := by
  have hup := var_up hA hB hp
  have hne : eVariationOn (qa A B p) (Icc 0 1) ≠ ⊤ := ne_top_of_le_ne_top ENNReal.ofReal_ne_top hup
  have hbv : BoundedVariationOn (qa A B p) (Icc 0 1) := hne
  set W := (eVariationOn (qa A B p) (Icc 0 1)).toReal with hW
  have hab : ∀ a b : ℝ, 0 < a → a ≤ b → b < 1 → Gp A B p b - Gp A B p a ≤ W := by
    intro a b ha hab hb
    have h1 := low_open hA hB hp hbv ha hab hb
    rw [variationOnFromTo.sub_right hbv.locallyBoundedVariationOn ⟨le_rfl, zero_le_one⟩
      ⟨by linarith, hb.le⟩ ⟨ha.le, by linarith⟩, variationOnFromTo.eq_of_le _ _ hab] at h1
    exact h1.trans (ENNReal.toReal_mono hne (eVariationOn.mono _ inter_subset_left))
  have t0 := eps_tend (x := 0) ⟨le_rfl, zero_le_one⟩
  have t1 := eps_tend (x := 1) ⟨zero_le_one, le_rfl⟩
  have R : Tendsto (fun e : ℝ => Gp A B p (e + (1 - 2 * e) * 1) - Gp A B p (e + (1 - 2 * e) * 0))
      (𝓝[>] 0) (𝓝 (Gp A B p 1 - Gp A B p 0)) :=
    ((Gp_cont hA hB hp 1 ⟨zero_le_one, le_rfl⟩).tendsto.comp t1).sub
      ((Gp_cont hA hB hp 0 ⟨le_rfl, zero_le_one⟩).tendsto.comp t0)
  have hG : Gp A B p 1 ≤ W := by
    have := le_of_tendsto (b := W) R ?_
    · simpa [Gp_zero] using this
    filter_upwards [Ioo_mem_nhdsGT (show (0:ℝ) < 1/2 by norm_num)] with e he
    exact hab _ _ (by nlinarith [he.1, he.2]) (by nlinarith [he.1, he.2]) (by nlinarith [he.1, he.2])
  refine le_antisymm hup ?_
  rw [← ENNReal.ofReal_toReal hne]
  exact ENNReal.ofReal_le_ofReal hG

lemma var_split {E : Type*} [PseudoEMetricSpace E] (f : ℝ → E) {a b c : ℝ}
    (hab : a ≤ b) (hbc : b ≤ c) :
    eVariationOn f (Icc a b) + eVariationOn f (Icc b c) = eVariationOn f (Icc a c) := by
  have := eVariationOn.Icc_add_Icc f (s := Icc a c) hab hbc ⟨hab, hbc⟩
  rwa [Set.inter_eq_right.2 (Icc_subset_Icc le_rfl hbc),
    Set.inter_eq_right.2 (Icc_subset_Icc hab le_rfl), Set.inter_self] at this

lemma loop_split {E : Type*} [PseudoEMetricSpace E] (f : ℝ → E) :
    eVariationOn f (Icc 0 4) = eVariationOn f (Icc 0 1) + eVariationOn f (Icc 1 2) +
      eVariationOn f (Icc 2 3) + eVariationOn f (Icc 3 4) := by
  rw [← var_split f (b := 1) (by norm_num) (by norm_num),
    ← var_split f (a := 1) (b := 2) (by norm_num) (by norm_num),
    ← var_split f (a := 2) (b := 3) (by norm_num) (by norm_num)]
  ring

lemma piece_var (R : ℝ × ℝ → ℝ × ℝ)
    (hR : ∀ z w, edist (emb (R z)) (emb (R w)) = edist (emb z) (emb w))
    (φ : ℝ → ℝ) (s : Set ℝ) (hφ : MonotoneOn φ s ∨ AntitoneOn φ s) (himg : φ '' s = Icc 0 1)
    (hf : ∀ t ∈ s, loopP A B p t = R (arc A B p (φ t))) :
    eVariationOn (loop A B p) s = eVariationOn (qa A B p) (Icc 0 1) := by
  rw [eVariationOn.eq_of_eqOn (f' := (fun z => emb (R (arc A B p z))) ∘ φ)
    (fun t ht => by simp only [Function.comp, loop_eq, hf t ht])]
  rcases hφ with h | h
  · rw [eVariationOn.comp_eq_of_monotoneOn _ _ h, himg]
    exact var_congr _ _ _ fun a _ b _ => hR _ _
  · rw [eVariationOn.comp_eq_of_antitoneOn _ _ h, himg]
    exact var_congr _ _ _ fun a _ b _ => hR _ _

end lengthsec

end X01

theorem X01.proof_3 : X01.stmt_3 := by
  intro A B hA hB p hp
  unfold X01.perimeter
  have hR1 : ∀ z w : ℝ × ℝ, edist (X01.emb ((fun z : ℝ × ℝ => (z.1, -z.2)) z))
      (X01.emb ((fun z : ℝ × ℝ => (z.1, -z.2)) w)) = edist (X01.emb z) (X01.emb w) := by
    intro z w; rw [X01.edist_emb, X01.edist_emb]; congr 2; ring
  have hR2 : ∀ z w : ℝ × ℝ, edist (X01.emb ((fun z : ℝ × ℝ => (-z.1, -z.2)) z))
      (X01.emb ((fun z : ℝ × ℝ => (-z.1, -z.2)) w)) = edist (X01.emb z) (X01.emb w) := by
    intro z w; rw [X01.edist_emb, X01.edist_emb]; congr 2; ring
  have hR3 : ∀ z w : ℝ × ℝ, edist (X01.emb ((fun z : ℝ × ℝ => (-z.1, z.2)) z))
      (X01.emb ((fun z : ℝ × ℝ => (-z.1, z.2)) w)) = edist (X01.emb z) (X01.emb w) := by
    intro z w; rw [X01.edist_emb, X01.edist_emb]; congr 2; ring
  rw [X01.loop_split,
    X01.piece_var (fun z => z) (fun z w => rfl) (fun t => t) (Set.Icc 0 1)
      (Or.inl monotoneOn_id) (by simp) (fun t ht => X01.loop_Q1 t ht),
    X01.piece_var _ hR1 (fun t => 2 - t) (Set.Icc 1 2)
      (Or.inr fun a _ b _ h => by simp only; linarith) (by rw [Set.image_const_sub_Icc]; norm_num)
      (fun t ht => X01.loop_Q2 hp t ht),
    X01.piece_var _ hR2 (fun t => t - 2) (Set.Icc 2 3)
      (Or.inl fun a _ b _ h => by simp only; linarith) (by rw [Set.image_sub_const_Icc]; norm_num)
      (fun t ht => X01.loop_Q3 hp t ht),
    X01.piece_var _ hR3 (fun t => 4 - t) (Set.Icc 3 4)
      (Or.inr fun a _ b _ h => by simp only; linarith) (by rw [Set.image_const_sub_Icc]; norm_num)
      (fun t ht => X01.loop_Q4 hp t ht),
    X01.var_eq hA hB hp, X01.Gp]
  unfold X01.spd X01.P
  rw [intervalIntegral.integral_const_mul]
  rw [show ENNReal.ofReal (1 / p * ∫ t in (0:ℝ)..1, X01.integrand A B p t) +
      ENNReal.ofReal (1 / p * ∫ t in (0:ℝ)..1, X01.integrand A B p t) +
      ENNReal.ofReal (1 / p * ∫ t in (0:ℝ)..1, X01.integrand A B p t) +
      ENNReal.ofReal (1 / p * ∫ t in (0:ℝ)..1, X01.integrand A B p t) =
      4 * ENNReal.ofReal (1 / p * ∫ t in (0:ℝ)..1, X01.integrand A B p t) by ring]
  rw [show (4:ENNReal) = ENNReal.ofReal 4 by norm_num, ← ENNReal.ofReal_mul (by norm_num)]
  congr 1; ring
