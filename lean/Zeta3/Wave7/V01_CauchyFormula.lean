import Mathlib
import Zeta3.Wave3.W02_Holo
import Zeta3.Wave3.W04_PerimeterIntegral
import Zeta3.Wave5.Y04_Fourier
import Zeta3.Wave6.Z02_Theta5
/-!
# V01 — Cauchy's formula for the squircle, `P = ∫_0^{2π} h_p = π a_0(h_p)`, and the
`τ⁵` expansion of `thm:pw` for the true perimeter

Target id: V01.
Paper labels: the paragraph before `thm:fourier` (main.tex line 1202: support function,
Cauchy's formula, `P(p) = π a_0(h_p)`); the displayed expansion of `thm:pw` (main.tex
lines 1183–1190); `thm:fourier` (main.tex 1204–1211) and `prop:theta` (main.tex 1220–1225)
at `k = 0`. The supplement does not prove Cauchy's formula: sec_g2_fourier.tex, remark (iii)
after the proof, says "$P=\pi a_0(h_p)$ (Cauchy's formula, used in the paper before the
theorem, not proved here)". The intermediate statements below follow the classical
argument (support function, normal-angle parametrisation, change of variables, symmetry).

Quoted verbatim.

main.tex 1202:

  What should we expand? Not the curve, but its \emph{support function}, the distance from the
  origin to the tangent line with normal $(\cos\theta,\sin\theta)$:
  $h_p(\theta)=\bigl(|\cos\theta|^{1/\tau}+|\sin\theta|^{1/\tau}\bigr)^{\tau}$, $\tau=1-\frac1p$.
  It is $\max(|\cos\theta|,|\sin\theta|)$ at the diamond and $|\cos\theta|+|\sin\theta|$ at the
  square, and Cauchy's classical formula gives $P=\int_0^{2\pi}h_p\dd\theta$. Since $h_p$ has
  period $\pi/2$, let's expand it as in Section~\ref{sec:zeta2}, with $T=\pi/2$:
  $h_p=\frac{a_0(h_p)}{2}+\sum_{k\ge1}a_k(h_p)\cos4k\theta$,
  $a_k(h_p)=\frac4\pi\int_0^{\pi/2}h_p\cos4k\theta\dd\theta$. Then $P(p)=\pi a_0(h_p)$: the
  perimeter is a single Fourier coefficient.

main.tex 1183–1190 (`thm:pw`, the displayed expansion):

  \frac{P}{\sqrt2}=4+\zeta(2)\tau^2+\zeta(3)\tau^3-\frac{33}{16}\zeta(4)\tau^4\\
  -\Bigl(\frac34\zeta(5)+\frac74\zeta(2)\zeta(3)\Bigr)\tau^5+O(\tau^6).

main.tex 1211 (after `thm:fourier`): "For $k=0$ this is Theorem~\ref{thm:pw}."
main.tex 1225 (after `prop:theta`): "At $k=0$, $\Theta_5$ is the $\tau^5$ coefficient of
Theorem~\ref{thm:pw}".

Reused, not redefined:
* `Y04.tau p = 1 - 1/p`, `Y04.h p θ = (|cos θ|^{1/τ} + |sin θ|^{1/τ})^τ` (the paper's `h_p`,
  meant for `p > 1`), `Y04.hDiamond θ = max |cos θ| |sin θ|` (the diamond `h_1`),
  `Y04.a k f = (4/π) ∫_0^{π/2} f θ cos(4kθ) dθ`, `Y04.z3 = (riemannZeta 3).re`
  (Zeta3.Wave5.Y04_Fourier); `Z02.z5 = (riemannZeta 5).re` (Zeta3.Wave6.Z02_Theta5).
* `W02.pIntegrand`, `W02.Pint p = (4/p) ∫_0^1 √(t^{2/p-2} + (1-t)^{2/p-2}) dt` (eq:Pint,
  Zeta3.Wave3.W02_Holo).
* `W04.Plane` (Euclidean plane), `W04.squircle p = {v | |v 0|^p + |v 1|^p = 1}`, and
  `W04.perimeter p = eVariationOn (W04.loop p) [0,4]`, the true Euclidean length of the curve
  (Zeta3.Wave3.W04_PerimeterIntegral; there `W04.stmt_7` proves `perimeter p = ofReal (Pint p)`,
  `W04.Pint` being literally `W02.Pint`).

Conventions.
* `ζ(3) = Y04.z3 = (riemannZeta 3).re`, `ζ(5) = Z02.z5 = (riemannZeta 5).re` (real parts of
  Mathlib's complex values; by `zeta_nat_eq_tsum_of_gt_one` these are `∑ 1/n^3`, `∑ 1/n^5`);
  `ζ(2) = π²/6`, `ζ(4) = π⁴/90` are written out.
* `Y04.h p` is used only for `p > 1` (`τ ∈ (0,1)`); for `p = 1` the formula has `τ = 0` and
  Lean's `1/0 = 0` would give the junk value `1`, so the diamond uses `Y04.hDiamond`, as the
  paper's sentence "It is max(...) at the diamond" does.
* `h_p` is continuous, so its interval integrals are genuine; `stmt_6`, `stmt_8` still assert
  the integrability explicitly. The quarter-length integral in `stmt_4` comes with the
  integrability of the (endpoint-singular) integrand.
* The support function is stated as the maximum of `⟨v, (cos θ, sin θ)⟩` over the curve
  (`IsGreatest`), the standard definition, which is what "the distance from the origin to the
  tangent line with normal (cos θ, sin θ)" means for a convex curve around the origin.
* Normal-angle parametrisation of the open first-quadrant arc `(t^{1/p}, (1-t)^{1/p})`,
  `0 < t < 1`: `θ(t) = arctan(((1-t)/t)^τ)` (`normalAngle`); its normal `(cos θ, sin θ)` is
  proportional to the gradient `(x^{p-1}, y^{p-1}) = (t^τ, (1-t)^τ)`.
* `O(τ⁶)` is explicit: `∃ C τ₀, 0 < τ₀ ∧ ∀ p > 1, τ ≤ τ₀ → |error| ≤ C τ⁶`, the form of
  `Y04.stmt_8` and `Z02.stmt_8`. The last statement is for the true perimeter
  `(W04.perimeter p).toReal` (finite by W04, so `toReal` is not a junk value).

Statements.
* stmt_1: for `p > 1`, `h_p(θ)` is the maximum of `x cos θ + y sin θ` over the squircle.
* stmt_2: for `p = 1`, the same with `max(|cos θ|, |sin θ|)`.
* stmt_3: normal-angle parametrisation of the quarter arc (`p > 1`): `θ(t)` is a strictly
  decreasing bijection of `(0,1)` onto `(0, π/2)`, the normal `(cos θ(t), sin θ(t))` is
  orthogonal to the velocity of `(t^{1/p}, (1-t)^{1/p})`, and the point attains the support
  value: `h_p(θ(t)) = t^{1/p} cos θ(t) + (1-t)^{1/p} sin θ(t)`.
* stmt_4: change of variables (`p > 1`): the quarter-arc length integral equals `∫_0^{π/2} h_p`:
  `(1/p) ∫_0^1 √(t^{2/p-2} + (1-t)^{2/p-2}) dt = ∫_0^{π/2} h_p(θ) dθ`.
* stmt_5: `h_p` has period `π/2` and `∫_0^{2π} h_p = 4 ∫_0^{π/2} h_p` (`p > 1`).
* stmt_6: Cauchy's formula for eq:Pint (`p > 1`): `P(p) = ∫_0^{2π} h_p = π a_0(h_p)`.
* stmt_7: Cauchy's formula for the true perimeter (`p > 1`).
* stmt_8: the case `p = 1` (diamond): `P(1) = ∫_0^{2π} h_1 = π a_0(h_1) = 4√2`, also for the
  true perimeter.
* stmt_9: `thm:pw`'s displayed expansion for the true perimeter, through `τ⁵` with `O(τ⁶)`.

NOT COVERED: the all-order claim and the weight/homogeneity claims of `thm:pw` (`log(P/4√2) =
Σ w_k τ^k` with `w_k` homogeneous of weight `k`), the rhombus and higher-dimensional
statements of `thm:pw`, and its proof sketch; the square end `h_∞ = |cos θ| + |sin θ|`
(a limit `p → ∞` statement, not needed for Cauchy's formula); the general convex-curve form
`length = ∫ (h + h'') dθ` (we state the direct substitution instead); the Fourier series
`h_p = a_0/2 + Σ a_k cos 4kθ` itself (its convergence is not used; only `a_0` enters);
`thm:fourier` and `prop:theta` for `k ≥ 1` (these are Y04, Z02).
-/

-- BEGIN STATEMENTS
namespace V01

open Real MeasureTheory Set

/-- The normal angle at the point `(t^{1/p}, (1-t)^{1/p})` of the first-quadrant arc,
`0 < t < 1`: `θ(t) = arctan(((1-t)/t)^τ)`, `τ = 1 - 1/p`. -/
noncomputable def normalAngle (p t : ℝ) : ℝ := Real.arctan (((1 - t) / t) ^ Y04.tau p)

/-- main.tex 1202, "its support function, the distance from the origin to the tangent line
with normal $(\cos\theta,\sin\theta)$: $h_p(\theta)=\bigl(|\cos\theta|^{1/\tau}+
|\sin\theta|^{1/\tau}\bigr)^{\tau}$": for `p > 1`, `h_p(θ)` is the maximum of
`x cos θ + y sin θ` over the squircle `|x|^p + |y|^p = 1`. -/
def stmt_1 : Prop :=
  ∀ p : ℝ, 1 < p → ∀ θ : ℝ,
    IsGreatest ((fun v : W04.Plane => v 0 * Real.cos θ + v 1 * Real.sin θ) '' W04.squircle p)
      (Y04.h p θ)

/-- main.tex 1202, "It is $\max(|\cos\theta|,|\sin\theta|)$ at the diamond": the support
function of `|x| + |y| = 1` is `max(|cos θ|, |sin θ|)`. -/
def stmt_2 : Prop :=
  ∀ θ : ℝ,
    IsGreatest ((fun v : W04.Plane => v 0 * Real.cos θ + v 1 * Real.sin θ) '' W04.squircle 1)
      (Y04.hDiamond θ)

/-- The normal-angle parametrisation used to prove Cauchy's formula: for `p > 1`,
`θ(t) = arctan(((1-t)/t)^τ)` maps `(0,1)` strictly decreasingly onto `(0, π/2)`; at
`t ∈ (0,1)` the normal `(cos θ(t), sin θ(t))` is orthogonal to the velocity
`((1/p) t^{1/p-1}, -(1/p)(1-t)^{1/p-1})` of the arc `(t^{1/p}, (1-t)^{1/p})`, and the point of
the arc attains the support value: `h_p(θ(t)) = t^{1/p} cos θ(t) + (1-t)^{1/p} sin θ(t)`. -/
def stmt_3 : Prop :=
  ∀ p : ℝ, 1 < p →
    StrictAntiOn (normalAngle p) (Ioo 0 1) ∧
    normalAngle p '' Ioo 0 1 = Ioo 0 (π / 2) ∧
    ∀ t ∈ Ioo (0 : ℝ) 1,
      1 / p * t ^ (1 / p - 1) * Real.cos (normalAngle p t)
          - 1 / p * (1 - t) ^ (1 / p - 1) * Real.sin (normalAngle p t) = 0 ∧
      Y04.h p (normalAngle p t) =
        t ^ (1 / p) * Real.cos (normalAngle p t) + (1 - t) ^ (1 / p) * Real.sin (normalAngle p t)

/-- The change of variables from the arc parameter `t` to the normal angle `θ` (the heart of
Cauchy's formula on one quarter): for `p > 1` the integrand of eq:Pint is integrable on `[0,1]`
and the quarter-arc length `(1/p) ∫_0^1 √(t^{2/p-2} + (1-t)^{2/p-2}) dt` equals
`∫_0^{π/2} h_p(θ) dθ`. -/
def stmt_4 : Prop :=
  ∀ p : ℝ, 1 < p →
    IntervalIntegrable (W02.pIntegrand p) volume 0 1 ∧
    1 / p * ∫ t in (0 : ℝ)..1, W02.pIntegrand p t = ∫ θ in (0 : ℝ)..(π / 2), Y04.h p θ

/-- main.tex 1202, "Since $h_p$ has period $\pi/2$": `h_p(θ + π/2) = h_p(θ)`, hence
`∫_0^{2π} h_p = 4 ∫_0^{π/2} h_p` (`p > 1`). -/
def stmt_5 : Prop :=
  ∀ p : ℝ, 1 < p →
    (∀ θ : ℝ, Y04.h p (θ + π / 2) = Y04.h p θ) ∧
    ∫ θ in (0 : ℝ)..(2 * π), Y04.h p θ = 4 * ∫ θ in (0 : ℝ)..(π / 2), Y04.h p θ

/-- main.tex 1202, "Cauchy's classical formula gives $P=\int_0^{2\pi}h_p\dd\theta$ ... Then
$P(p)=\pi a_0(h_p)$", for `P` given by eq:Pint and `p > 1`: `h_p` is integrable on
`[0, 2π]`, `P(p) = ∫_0^{2π} h_p(θ) dθ` and `P(p) = π a_0(h_p)`. -/
def stmt_6 : Prop :=
  ∀ p : ℝ, 1 < p →
    IntervalIntegrable (Y04.h p) volume 0 (2 * π) ∧
    W02.Pint p = ∫ θ in (0 : ℝ)..(2 * π), Y04.h p θ ∧
    W02.Pint p = π * Y04.a 0 (Y04.h p)

/-- Cauchy's formula for the true perimeter (Euclidean length of the closed curve
`|x|^p + |y|^p = 1`), `p > 1`: `perimeter = ∫_0^{2π} h_p` and `P(p) = π a_0(h_p)`. -/
def stmt_7 : Prop :=
  ∀ p : ℝ, 1 < p →
    W04.perimeter p = ENNReal.ofReal (∫ θ in (0 : ℝ)..(2 * π), Y04.h p θ) ∧
    (W04.perimeter p).toReal = π * Y04.a 0 (Y04.h p)

/-- Cauchy's formula at the diamond `p = 1` (support function `h_1 = max(|cos θ|, |sin θ|)`):
`P(1) = ∫_0^{2π} h_1 = π a_0(h_1) = 4√2`, for eq:Pint and for the true perimeter. -/
def stmt_8 : Prop :=
  IntervalIntegrable Y04.hDiamond volume 0 (2 * π) ∧
  W02.Pint 1 = ∫ θ in (0 : ℝ)..(2 * π), Y04.hDiamond θ ∧
  W02.Pint 1 = π * Y04.a 0 Y04.hDiamond ∧
  π * Y04.a 0 Y04.hDiamond = 4 * Real.sqrt 2 ∧
  W04.perimeter 1 = ENNReal.ofReal (∫ θ in (0 : ℝ)..(2 * π), Y04.hDiamond θ)

/-- `thm:pw`, the displayed expansion, for the true perimeter `P(p)`:
"$\frac{P}{\sqrt2}=4+\zeta(2)\tau^2+\zeta(3)\tau^3-\frac{33}{16}\zeta(4)\tau^4
-\Bigl(\frac34\zeta(5)+\frac74\zeta(2)\zeta(3)\Bigr)\tau^5+O(\tau^6)$":
there are `C` and `τ₀ > 0` such that for every `p > 1` with `τ = 1 - 1/p ≤ τ₀`,
`|P(p)/√2 - (4 + ζ(2)τ² + ζ(3)τ³ - (33/16)ζ(4)τ⁴ - ((3/4)ζ(5) + (7/4)ζ(2)ζ(3))τ⁵)| ≤ C τ⁶`. -/
def stmt_9 : Prop :=
  ∃ C τ₀ : ℝ, 0 < τ₀ ∧ ∀ p : ℝ, 1 < p → Y04.tau p ≤ τ₀ →
    |(W04.perimeter p).toReal / Real.sqrt 2 -
        (4 + (π ^ 2 / 6) * Y04.tau p ^ 2 + Y04.z3 * Y04.tau p ^ 3
          - 33 / 16 * (π ^ 4 / 90) * Y04.tau p ^ 4
          - (3 / 4 * Z02.z5 + 7 / 4 * ((π ^ 2 / 6) * Y04.z3)) * Y04.tau p ^ 5)|
      ≤ C * Y04.tau p ^ 6

end V01
-- END STATEMENTS

open Real MeasureTheory Set

namespace V01

lemma signed_pt {m c : ℝ} (hm : 0 ≤ m) : ∃ x : ℝ, |x| = m ∧ x * c = m * |c| := by
  rcases le_total 0 c with hc | hc
  · exact ⟨m, abs_of_nonneg hm, by rw [abs_of_nonneg hc]⟩
  · exact ⟨-m, by rw [abs_neg, abs_of_nonneg hm], by rw [abs_of_nonpos hc]; ring⟩

lemma q_facts {p : ℝ} (hp : 1 < p) :
    1 < 1 / Y04.tau p ∧ (1 / Y04.tau p - 1) * p = 1 / Y04.tau p ∧
      p.HolderConjugate (1 / Y04.tau p) := by
  have ht0 := Y04.tau_pos hp
  have ht1 := Y04.tau_lt_one hp
  have hp0 : (0 : ℝ) < p := by linarith
  refine ⟨by rw [lt_div_iff₀ ht0]; linarith, ?_, ?_⟩
  · unfold Y04.tau
    have : p - 1 ≠ 0 := by linarith
    field_simp
    ring
  · rw [Real.holderConjugate_iff]
    refine ⟨hp, ?_⟩
    unfold Y04.tau
    field_simp
    ring

lemma h_pos {p : ℝ} (hp : 1 < p) (θ : ℝ) : 0 < Y04.h p θ := by
  have ht0 := Y04.tau_pos hp
  unfold Y04.h
  apply Real.rpow_pos_of_pos
  have hq : 0 < 1 / Y04.tau p := by positivity
  rcases eq_or_ne (Real.cos θ) 0 with hc | hc
  · have : Real.sin θ ≠ 0 := by
      intro hs; have := Real.sin_sq_add_cos_sq θ; rw [hc, hs] at this; norm_num at this
    have := Real.rpow_pos_of_pos (abs_pos.2 this) (1 / Y04.tau p)
    have := Real.rpow_nonneg (abs_nonneg (Real.cos θ)) (1 / Y04.tau p)
    linarith
  · have := Real.rpow_pos_of_pos (abs_pos.2 hc) (1 / Y04.tau p)
    have := Real.rpow_nonneg (abs_nonneg (Real.sin θ)) (1 / Y04.tau p)
    linarith

lemma rpow_succ' {x y : ℝ} (hx : 0 ≤ x) (hy : y + 1 ≠ 0) : x ^ y * x = x ^ (y + 1) :=
  (Real.rpow_add_one' hx hy).symm

lemma support_p {p : ℝ} (hp : 1 < p) (θ : ℝ) :
    IsGreatest ((fun v : W04.Plane => v 0 * Real.cos θ + v 1 * Real.sin θ) '' W04.squircle p)
      (Y04.h p θ) := by
  obtain ⟨hq1, hqp, hpq⟩ := q_facts hp
  have ht0 := Y04.tau_pos hp
  have hp0 : (0 : ℝ) < p := by linarith
  set q := 1 / Y04.tau p with hqdef
  have hτq : Y04.tau p = 1 / q := by rw [hqdef]; field_simp
  set c := Real.cos θ
  set s := Real.sin θ
  set H := Y04.h p θ with hH
  have hH0 : 0 < H := h_pos hp θ
  have hHq : H ^ q = |c| ^ q + |s| ^ q := by
    rw [hH]; unfold Y04.h
    rw [← Real.rpow_mul (by positivity), ← hqdef, hτq, one_div_mul_cancel (by positivity),
      Real.rpow_one]
  constructor
  · -- the maximiser
    obtain ⟨x, hx, hxc⟩ := signed_pt (c := c) (by positivity : 0 ≤ |c| ^ (q - 1) / H ^ (q - 1))
    obtain ⟨y, hy, hys⟩ := signed_pt (c := s) (by positivity : 0 ≤ |s| ^ (q - 1) / H ^ (q - 1))
    refine ⟨!₂[x, y], ?_, ?_⟩
    · show |x| ^ p + |y| ^ p = 1
      rw [hx, hy, Real.div_rpow (by positivity) (by positivity), Real.div_rpow (by positivity) (by positivity),
        ← Real.rpow_mul (abs_nonneg _), ← Real.rpow_mul (abs_nonneg _), ← Real.rpow_mul hH0.le, hqp,
        ← add_div, hHq, div_self]
      rw [← hHq]; exact (Real.rpow_pos_of_pos hH0 q).ne'
    · show x * c + y * s = H
      rw [hxc, hys, div_mul_eq_mul_div, div_mul_eq_mul_div, rpow_succ' (abs_nonneg _) (by linarith),
        rpow_succ' (abs_nonneg _) (by linarith), ← add_div, sub_add_cancel, ← hHq,
        Real.rpow_sub hH0, Real.rpow_one]
      field_simp
  · rintro _ ⟨v, hv, rfl⟩
    have hv' : |v 0| ^ p + |v 1| ^ p = 1 := hv
    have := Real.inner_le_Lp_mul_Lq (s := Finset.univ) ![v 0, v 1] ![c, s] hpq
    simp only [Fin.sum_univ_two, Matrix.cons_val_zero, Matrix.cons_val_one] at this
    rw [hv', Real.one_rpow, one_mul, hqdef, one_div_one_div] at this
    simpa [hH, Y04.h, ← hqdef] using this

end V01

namespace V01

lemma support_one (θ : ℝ) :
    IsGreatest ((fun v : W04.Plane => v 0 * Real.cos θ + v 1 * Real.sin θ) '' W04.squircle 1)
      (Y04.hDiamond θ) := by
  set c := Real.cos θ
  set s := Real.sin θ
  constructor
  · rcases le_total |s| |c| with hcs | hcs
    · obtain ⟨x, hx, hxc⟩ := signed_pt (m := 1) (c := c) zero_le_one
      refine ⟨!₂[x, 0], ?_, ?_⟩
      · show |x| ^ (1 : ℝ) + |(0 : ℝ)| ^ (1 : ℝ) = 1
        rw [hx]; simp
      · show x * c + 0 * s = Y04.hDiamond θ
        rw [hxc]; unfold Y04.hDiamond; rw [max_eq_left hcs]; ring
    · obtain ⟨y, hy, hys⟩ := signed_pt (m := 1) (c := s) zero_le_one
      refine ⟨!₂[0, y], ?_, ?_⟩
      · show |(0 : ℝ)| ^ (1 : ℝ) + |y| ^ (1 : ℝ) = 1
        rw [hy]; simp
      · show 0 * c + y * s = Y04.hDiamond θ
        rw [hys]; unfold Y04.hDiamond; rw [max_eq_right hcs]; ring
  · rintro _ ⟨v, hv, rfl⟩
    have hv' : |v 0| ^ (1 : ℝ) + |v 1| ^ (1 : ℝ) = 1 := hv
    rw [Real.rpow_one, Real.rpow_one] at hv'
    show v 0 * c + v 1 * s ≤ max |c| |s|
    have h1 : v 0 * c ≤ |v 0| * max |c| |s| := by
      calc v 0 * c ≤ |v 0 * c| := le_abs_self _
        _ = |v 0| * |c| := abs_mul _ _
        _ ≤ |v 0| * max |c| |s| := mul_le_mul_of_nonneg_left (le_max_left _ _) (abs_nonneg _)
    have h2 : v 1 * s ≤ |v 1| * max |c| |s| := by
      calc v 1 * s ≤ |v 1 * s| := le_abs_self _
        _ = |v 1| * |s| := abs_mul _ _
        _ ≤ |v 1| * max |c| |s| := mul_le_mul_of_nonneg_left (le_max_right _ _) (abs_nonneg _)
    have h3 : |v 0| * max |c| |s| + |v 1| * max |c| |s| = max |c| |s| := by
      rw [← add_mul, hv', one_mul]
    linarith

lemma cos_sin_arctan_div {A B : ℝ} (hA : 0 < A) (_hB : 0 ≤ B) :
    Real.cos (Real.arctan (B / A)) = A / Real.sqrt (A ^ 2 + B ^ 2) ∧
    Real.sin (Real.arctan (B / A)) = B / Real.sqrt (A ^ 2 + B ^ 2) := by
  have hs : Real.sqrt (1 + (B / A) ^ 2) = Real.sqrt (A ^ 2 + B ^ 2) / A := by
    rw [show 1 + (B / A) ^ 2 = (A ^ 2 + B ^ 2) / A ^ 2 by field_simp, Real.sqrt_div' _ (by positivity),
      Real.sqrt_sq hA.le]
  have hD : 0 < Real.sqrt (A ^ 2 + B ^ 2) := Real.sqrt_pos.2 (by positivity)
  rw [Real.cos_arctan, Real.sin_arctan, hs]
  constructor <;> field_simp

lemma exps {p : ℝ} (hp : 1 < p) {x : ℝ} (hx : 0 < x) :
    x ^ (1 / p) * x ^ Y04.tau p = x ∧ x ^ (1 / p - 1) * x ^ Y04.tau p = 1 ∧
      (x ^ Y04.tau p) ^ (1 / Y04.tau p) = x := by
  have ht0 := Y04.tau_pos hp
  refine ⟨?_, ?_, ?_⟩
  · rw [← Real.rpow_add hx, show 1 / p + Y04.tau p = 1 by unfold Y04.tau; ring, Real.rpow_one]
  · rw [← Real.rpow_add hx, show 1 / p - 1 + Y04.tau p = 0 by unfold Y04.tau; ring, Real.rpow_zero]
  · rw [← Real.rpow_mul hx.le, mul_one_div_cancel ht0.ne', Real.rpow_one]

/-- `normalAngle` as `arctan (B/A)` with `A = t^τ`, `B = (1-t)^τ` -/
lemma normalAngle_eq {p t : ℝ} (ht0 : 0 < t) (ht1 : t < 1) :
    normalAngle p t = Real.arctan ((1 - t) ^ Y04.tau p / t ^ Y04.tau p) := by
  unfold normalAngle; rw [Real.div_rpow (by linarith) ht0.le]

lemma cos_sin_nA {p t : ℝ} (ht0 : 0 < t) (ht1 : t < 1) :
    Real.cos (normalAngle p t) = t ^ Y04.tau p /
        Real.sqrt ((t ^ Y04.tau p) ^ 2 + ((1 - t) ^ Y04.tau p) ^ 2) ∧
    Real.sin (normalAngle p t) = (1 - t) ^ Y04.tau p /
        Real.sqrt ((t ^ Y04.tau p) ^ 2 + ((1 - t) ^ Y04.tau p) ^ 2) := by
  rw [normalAngle_eq ht0 ht1]
  exact cos_sin_arctan_div (Real.rpow_pos_of_pos ht0 _) (Real.rpow_nonneg (by linarith) _)

lemma nA_strictAnti {p : ℝ} (hp : 1 < p) : StrictAntiOn (normalAngle p) (Ioo 0 1) := by
  have ht0 := Y04.tau_pos hp
  intro s hs t ht hst
  unfold normalAngle
  apply Real.arctan_strictMono
  apply Real.rpow_lt_rpow (div_nonneg (by linarith [ht.2]) ht.1.le) _ ht0
  rw [div_lt_div_iff₀ ht.1 hs.1]; nlinarith [hs.1]

lemma nA_image {p : ℝ} (hp : 1 < p) : normalAngle p '' Ioo 0 1 = Ioo 0 (π / 2) := by
  have ht0 := Y04.tau_pos hp
  ext θ; constructor
  · rintro ⟨t, ⟨h0, h1⟩, rfl⟩
    have hu : 0 < ((1 - t) / t) ^ Y04.tau p := Real.rpow_pos_of_pos (div_pos (by linarith) h0) _
    exact ⟨Real.arctan_pos.2 hu, Real.arctan_lt_pi_div_two _⟩
  · rintro ⟨h0, h1⟩
    have hpi := Real.pi_pos
    have htan : 0 < Real.tan θ := Real.tan_pos_of_pos_of_lt_pi_div_two h0 h1
    set w := Real.tan θ ^ (1 / Y04.tau p) with hw
    have hw0 : 0 < w := Real.rpow_pos_of_pos htan _
    refine ⟨1 / (1 + w), ⟨by positivity, ?_⟩, ?_⟩
    · rw [div_lt_one (by positivity)]; linarith
    · unfold normalAngle
      rw [show (1 - 1 / (1 + w)) / (1 / (1 + w)) = w by field_simp; ring, hw,
        ← Real.rpow_mul htan.le, one_div_mul_cancel ht0.ne', Real.rpow_one,
        Real.arctan_tan (by linarith) h1]

lemma nA_props {p : ℝ} (hp : 1 < p) {t : ℝ} (ht : t ∈ Ioo (0 : ℝ) 1) :
    1 / p * t ^ (1 / p - 1) * Real.cos (normalAngle p t)
          - 1 / p * (1 - t) ^ (1 / p - 1) * Real.sin (normalAngle p t) = 0 ∧
      Y04.h p (normalAngle p t) =
        t ^ (1 / p) * Real.cos (normalAngle p t) + (1 - t) ^ (1 / p) * Real.sin (normalAngle p t) := by
  obtain ⟨h0, h1⟩ := ht
  have h1' : 0 < 1 - t := by linarith
  have ht0 := Y04.tau_pos hp
  obtain ⟨hc, hs⟩ := cos_sin_nA (p := p) h0 h1
  obtain ⟨a1, a2, a3⟩ := exps hp h0
  obtain ⟨b1, b2, b3⟩ := exps hp h1'
  set A := t ^ Y04.tau p
  set B := (1 - t) ^ Y04.tau p
  have hA : 0 < A := Real.rpow_pos_of_pos h0 _
  have hB : 0 < B := Real.rpow_pos_of_pos h1' _
  set D := Real.sqrt (A ^ 2 + B ^ 2)
  have hD : 0 < D := Real.sqrt_pos.2 (by positivity)
  have eA : t ^ (1 / p - 1) = 1 / A := by rw [eq_div_iff hA.ne']; exact a2
  have eB : (1 - t) ^ (1 / p - 1) = 1 / B := by rw [eq_div_iff hB.ne']; exact b2
  have fA : t ^ (1 / p) = t / A := by rw [eq_div_iff hA.ne']; exact a1
  have fB : (1 - t) ^ (1 / p) = (1 - t) / B := by rw [eq_div_iff hB.ne']; exact b1
  refine ⟨?_, ?_⟩
  · rw [hc, hs, eA, eB]; field_simp; ring
  · unfold Y04.h
    rw [hc, hs, fA, fB]
    rw [abs_of_pos (by positivity), abs_of_pos (by positivity),
      Real.div_rpow hA.le hD.le, Real.div_rpow hB.le hD.le, a3, b3, ← add_div,
      show t + (1 - t) = 1 by ring, Real.div_rpow zero_le_one (by positivity), Real.one_rpow,
      ← Real.rpow_mul hD.le, one_div_mul_cancel ht0.ne', Real.rpow_one]
    field_simp
    ring

end V01

namespace V01

open Filter Topology

/-- the derivative of the normal angle -/
noncomputable def dTheta (p t : ℝ) : ℝ :=
  -(1 / (1 + (((1 - t) / t) ^ Y04.tau p) ^ 2) *
      (1 / t ^ 2 * Y04.tau p * ((1 - t) / t) ^ (Y04.tau p - 1)))

lemma dTheta_neg {p t : ℝ} (hp : 1 < p) (ht : t ∈ Ioo (0 : ℝ) 1) : dTheta p t < 0 := by
  have ht0 := Y04.tau_pos hp
  have hr : 0 < (1 - t) / t := div_pos (by linarith [ht.2]) ht.1
  have := ht.1
  unfold dTheta
  have : 0 < ((1 - t) / t) ^ (Y04.tau p - 1) := Real.rpow_pos_of_pos hr _
  have : 0 < 1 / (1 + (((1 - t) / t) ^ Y04.tau p) ^ 2) *
      (1 / t ^ 2 * Y04.tau p * ((1 - t) / t) ^ (Y04.tau p - 1)) := by positivity
  linarith

lemma nA_hasDeriv {p t : ℝ} (ht : t ∈ Ioo (0 : ℝ) 1) :
    HasDerivAt (normalAngle p) (dTheta p t) t := by
  obtain ⟨h0, h1⟩ := ht
  have hr : HasDerivAt (fun s => (1 - s) / s) (((-1) * t - (1 - t) * 1) / t ^ 2) t :=
    ((hasDerivAt_id t).const_sub 1).div (hasDerivAt_id t) h0.ne'
  have hu := hr.rpow_const (p := Y04.tau p) (Or.inl (div_pos (by linarith) h0).ne')
  have hθ := hu.arctan
  unfold normalAngle
  refine hθ.congr_deriv ?_
  unfold dTheta
  rw [show ((-1) * t - (1 - t) * 1) / t ^ 2 = -(1 / t ^ 2) by field_simp; ring]
  ring

lemma pInt_eq {p : ℝ} (hp : 1 < p) {t : ℝ} (ht : t ∈ Ioo (0 : ℝ) 1) :
    t ^ (1 / p - 1) * Real.sin (normalAngle p t) + (1 - t) ^ (1 / p - 1) * Real.cos (normalAngle p t)
      = W02.pIntegrand p t := by
  obtain ⟨h0, h1⟩ := ht
  have h1' : 0 < 1 - t := by linarith
  obtain ⟨hc, hs⟩ := cos_sin_nA (p := p) h0 h1
  obtain ⟨_, a2, _⟩ := exps hp h0
  obtain ⟨_, b2, _⟩ := exps hp h1'
  set A := t ^ Y04.tau p
  set B := (1 - t) ^ Y04.tau p
  have hA : 0 < A := Real.rpow_pos_of_pos h0 _
  have hB : 0 < B := Real.rpow_pos_of_pos h1' _
  have hD2 : Real.sqrt (A ^ 2 + B ^ 2) ^ 2 = A ^ 2 + B ^ 2 := Real.sq_sqrt (by positivity)
  set D := Real.sqrt (A ^ 2 + B ^ 2)
  have hD : 0 < D := Real.sqrt_pos.2 (by positivity)
  have eA : t ^ (1 / p - 1) = 1 / A := by rw [eq_div_iff hA.ne']; exact a2
  have eB : (1 - t) ^ (1 / p - 1) = 1 / B := by rw [eq_div_iff hB.ne']; exact b2
  have hP : W02.pIntegrand p t = D / (A * B) := by
    unfold W02.pIntegrand
    rw [W04.rpow_sq' h0.le, W04.rpow_sq' h1'.le, eA, eB]
    have e : (D / (A * B)) ^ 2 = (A ^ 2 + B ^ 2) / (A * B) ^ 2 := by rw [div_pow, hD2]
    rw [show (1 / A) ^ 2 + (1 / B) ^ 2 = (D / (A * B)) ^ 2 by rw [e]; field_simp; ring]
    exact Real.sqrt_sq (by positivity)
  rw [hP, hc, hs, eA, eB]
  field_simp
  rw [hD2]; ring

/-- the boundary term `F = ⟨x, n'⟩`, extended continuously to `[0,1]` -/
noncomputable def Fb (p t : ℝ) : ℝ :=
  ((1 - t) ^ (1 / p) * t ^ Y04.tau p - t ^ (1 / p) * (1 - t) ^ Y04.tau p) /
    Real.sqrt ((t ^ Y04.tau p) ^ 2 + ((1 - t) ^ Y04.tau p) ^ 2)

lemma Fb_cont {p : ℝ} (hp : 1 < p) : ContinuousOn (Fb p) (Icc 0 1) := by
  have ht0 := Y04.tau_pos hp
  have hp0 : 0 < 1 / p := one_div_pos.2 (by linarith)
  have c1 : Continuous (fun t : ℝ => t ^ (1 / p)) := Real.continuous_rpow_const hp0.le
  have c2 : Continuous (fun t : ℝ => t ^ Y04.tau p) := Real.continuous_rpow_const ht0.le
  have c3 : Continuous (fun t : ℝ => (1 - t) ^ (1 / p)) := c1.comp (continuous_const.sub continuous_id)
  have c4 : Continuous (fun t : ℝ => (1 - t) ^ Y04.tau p) :=
    c2.comp (continuous_const.sub continuous_id)
  have hn : Continuous (fun t : ℝ => (1 - t) ^ (1 / p) * t ^ Y04.tau p - t ^ (1 / p) * (1 - t) ^ Y04.tau p) :=
    (c3.mul c2).sub (c1.mul c4)
  have hd : Continuous (fun t : ℝ => Real.sqrt ((t ^ Y04.tau p) ^ 2 + ((1 - t) ^ Y04.tau p) ^ 2)) :=
    ((c2.pow 2).add (c4.pow 2)).sqrt
  refine hn.continuousOn.div hd.continuousOn (fun t ht => ?_)
  apply (Real.sqrt_pos.2 _).ne'
  rcases eq_or_lt_of_le ht.1 with h | h
  · rw [← h, sub_zero, Real.one_rpow]
    have : 0 ≤ ((0:ℝ) ^ Y04.tau p) ^ 2 := sq_nonneg _
    linarith
  · have := Real.rpow_pos_of_pos h (Y04.tau p)
    have : 0 ≤ ((1 - t) ^ Y04.tau p) ^ 2 := sq_nonneg _
    nlinarith

lemma Fb_zero {p : ℝ} (hp : 1 < p) : Fb p 0 = 0 := by
  have ht0 := Y04.tau_pos hp
  have : p⁻¹ ≠ 0 := inv_ne_zero (by linarith)
  simp [Fb, Real.zero_rpow ht0.ne', Real.zero_rpow this]

lemma Fb_one {p : ℝ} (hp : 1 < p) : Fb p 1 = 0 := by
  have ht0 := Y04.tau_pos hp
  have : p⁻¹ ≠ 0 := inv_ne_zero (by linarith)
  simp [Fb, Real.zero_rpow ht0.ne', Real.zero_rpow this]

/-- the integrand on the `θ` side, pulled back -/
noncomputable def kθ (p t : ℝ) : ℝ := |dTheta p t| * Y04.h p (normalAngle p t)

lemma Fb_hasDeriv {p : ℝ} (hp : 1 < p) {t : ℝ} (ht : t ∈ Ioo (0 : ℝ) 1) :
    HasDerivAt (Fb p) (kθ p t - 1 / p * W02.pIntegrand p t) t := by
  obtain ⟨h0, h1⟩ := ht
  have h1' : 0 < 1 - t := by linarith
  have hθ := nA_hasDeriv (p := p) ⟨h0, h1⟩
  have hx : HasDerivAt (fun s : ℝ => s ^ (1 / p)) (1 / p * t ^ (1 / p - 1)) t :=
    Real.hasDerivAt_rpow_const (Or.inl h0.ne')
  have hy : HasDerivAt (fun s : ℝ => (1 - s) ^ (1 / p)) ((-1) * (1 / p) * (1 - t) ^ (1 / p - 1)) t :=
    ((hasDerivAt_id t).const_sub 1).rpow_const (Or.inl h1'.ne')
  have hF := (hx.mul hθ.sin).neg.add (hy.mul hθ.cos)
  have heq : (fun s => -(s ^ (1 / p) * Real.sin (normalAngle p s)) +
      (1 - s) ^ (1 / p) * Real.cos (normalAngle p s)) =ᶠ[𝓝 t] Fb p := by
    filter_upwards [Ioo_mem_nhds h0 h1] with s hs
    obtain ⟨hc, hs'⟩ := cos_sin_nA (p := p) hs.1 hs.2
    rw [hc, hs']; unfold Fb; ring
  refine (hF.congr_of_eventuallyEq heq.symm).congr_deriv ?_
  have hL := pInt_eq hp ⟨h0, h1⟩
  have hH := (nA_props hp ⟨h0, h1⟩).2
  unfold kθ
  rw [abs_of_neg (dTheta_neg hp ⟨h0, h1⟩), hH, ← hL]
  ring

lemma quarter_cov {p : ℝ} (hp : 1 < p) :
    IntervalIntegrable (W02.pIntegrand p) volume 0 1 ∧
    1 / p * ∫ t in (0 : ℝ)..1, W02.pIntegrand p t = ∫ θ in (0 : ℝ)..(π / 2), Y04.h p θ := by
  have hpi := Real.pi_pos
  have hP : IntervalIntegrable (W02.pIntegrand p) volume 0 1 := W04.pInt_iInt hp.le
  refine ⟨hP, ?_⟩
  have hder : ∀ t ∈ Ioo (0 : ℝ) 1, HasDerivWithinAt (normalAngle p) (dTheta p t) (Ioo 0 1) t :=
    fun t ht => (nA_hasDeriv ht).hasDerivWithinAt
  have hinj : InjOn (normalAngle p) (Ioo 0 1) := (nA_strictAnti hp).injOn
  have himg := integral_image_eq_integral_abs_deriv_smul measurableSet_Ioo hder hinj (Y04.h p)
  have hint := (integrableOn_image_iff_integrableOn_abs_deriv_smul measurableSet_Ioo hder hinj
    (Y04.h p)).1
  rw [nA_image hp] at himg hint
  have hhI : IntegrableOn (Y04.h p) (Ioo 0 (π / 2)) :=
    ((Y04.h_cont hp).integrableOn_Icc).mono_set Ioo_subset_Icc_self
  have hk : IntervalIntegrable (kθ p) volume 0 1 := by
    rw [intervalIntegrable_iff_integrableOn_Ioo_of_le zero_le_one]
    have := hint hhI
    simp only [smul_eq_mul] at this
    exact this
  have hf' : IntervalIntegrable (fun t => kθ p t - 1 / p * W02.pIntegrand p t) volume 0 1 :=
    hk.sub (hP.const_mul _)
  have hFTC := intervalIntegral.integral_eq_sub_of_hasDerivAt_of_le zero_le_one (Fb_cont hp)
    (fun t ht => Fb_hasDeriv hp ht) hf'
  rw [Fb_zero hp, Fb_one hp, sub_zero, intervalIntegral.integral_sub hk (hP.const_mul _),
    intervalIntegral.integral_const_mul] at hFTC
  have e1 : ∫ θ in (0 : ℝ)..(π / 2), Y04.h p θ = ∫ θ in Ioo (0 : ℝ) (π / 2), Y04.h p θ := by
    rw [intervalIntegral.integral_of_le (by linarith), integral_Ioc_eq_integral_Ioo]
  have e2 : ∫ t in (0 : ℝ)..1, kθ p t = ∫ t in Ioo (0 : ℝ) 1, kθ p t := by
    rw [intervalIntegral.integral_of_le zero_le_one, integral_Ioc_eq_integral_Ioo]
  rw [e1, himg]
  simp only [smul_eq_mul]
  rw [show (∫ x in Ioo (0 : ℝ) 1, |dTheta p x| * Y04.h p (normalAngle p x))
      = ∫ t in Ioo (0 : ℝ) 1, kθ p t from rfl, ← e2]
  linarith

end V01

namespace V01

lemma int_period {f : ℝ → ℝ} (hf : Continuous f) (hper : ∀ θ, f (θ + π / 2) = f θ) :
    ∫ θ in (0 : ℝ)..(2 * π), f θ = 4 * ∫ θ in (0 : ℝ)..(π / 2), f θ := by
  have hP : Function.Periodic f (π / 2) := hper
  have := hP.intervalIntegral_add_zsmul_eq 4 0 (fun a b => hf.intervalIntegrable a b)
  rw [show (0 : ℝ) + (4 : ℤ) • (π / 2) = 2 * π by simp; ring, zero_add] at this
  rw [this]; simp

lemma h_period {p : ℝ} (θ : ℝ) : Y04.h p (θ + π / 2) = Y04.h p θ := by
  unfold Y04.h
  rw [Real.cos_add_pi_div_two, Real.sin_add_pi_div_two, abs_neg, add_comm]

lemma hDiamond_period (θ : ℝ) : Y04.hDiamond (θ + π / 2) = Y04.hDiamond θ := by
  unfold Y04.hDiamond
  rw [Real.cos_add_pi_div_two, Real.sin_add_pi_div_two, abs_neg, max_comm]

lemma a_zero (f : ℝ → ℝ) : π * Y04.a 0 f = 4 * ∫ θ in (0 : ℝ)..(π / 2), f θ := by
  have hpi := Real.pi_pos
  unfold Y04.a
  simp only [Nat.cast_zero, mul_zero, zero_mul, Real.cos_zero, mul_one]
  field_simp

lemma cauchy_p {p : ℝ} (hp : 1 < p) :
    IntervalIntegrable (Y04.h p) volume 0 (2 * π) ∧
    W02.Pint p = ∫ θ in (0 : ℝ)..(2 * π), Y04.h p θ ∧
    W02.Pint p = π * Y04.a 0 (Y04.h p) := by
  have hq := (quarter_cov hp).2
  have hper := int_period (Y04.h_cont hp) (fun θ => h_period θ)
  have hP : W02.Pint p = 4 * ∫ θ in (0 : ℝ)..(π / 2), Y04.h p θ := by
    unfold W02.Pint; rw [← hq]; ring
  refine ⟨(Y04.h_cont hp).intervalIntegrable _ _, by rw [hP, hper], by rw [hP, a_zero]⟩

lemma Pint_nonneg {p : ℝ} (hp : 1 < p) : 0 ≤ W02.Pint p := by
  have hpi := Real.pi_pos
  rw [(cauchy_p hp).2.1]
  apply intervalIntegral.integral_nonneg (by positivity)
  intro θ _; unfold Y04.h; positivity

lemma Pint_W04 (p : ℝ) : W04.Pint p = W02.Pint p := rfl

lemma cauchy_true {p : ℝ} (hp : 1 < p) :
    W04.perimeter p = ENNReal.ofReal (∫ θ in (0 : ℝ)..(2 * π), Y04.h p θ) ∧
    (W04.perimeter p).toReal = π * Y04.a 0 (Y04.h p) := by
  have h7 := W04.proof_7 p hp.le
  rw [Pint_W04] at h7
  obtain ⟨_, h2, h3⟩ := cauchy_p hp
  refine ⟨by rw [h7, h2], ?_⟩
  rw [h7, ENNReal.toReal_ofReal (Pint_nonneg hp), h3]

lemma diamond :
    IntervalIntegrable Y04.hDiamond volume 0 (2 * π) ∧
    W02.Pint 1 = ∫ θ in (0 : ℝ)..(2 * π), Y04.hDiamond θ ∧
    W02.Pint 1 = π * Y04.a 0 Y04.hDiamond ∧
    π * Y04.a 0 Y04.hDiamond = 4 * Real.sqrt 2 ∧
    W04.perimeter 1 = ENNReal.ofReal (∫ θ in (0 : ℝ)..(2 * π), Y04.hDiamond θ) := by
  have hpi := Real.pi_pos
  have hA : π * Y04.a 0 Y04.hDiamond = 4 * Real.sqrt 2 := by
    rw [Y04.proof_1 0]; field_simp; norm_num
  have hper := int_period Y04.hDiamond_cont hDiamond_period
  have hI : ∫ θ in (0 : ℝ)..(2 * π), Y04.hDiamond θ = 4 * Real.sqrt 2 := by
    rw [hper, ← a_zero, hA]
  have hP : W02.Pint 1 = 4 * Real.sqrt 2 := by rw [← Pint_W04]; exact W04.proof_8
  refine ⟨Y04.hDiamond_cont.intervalIntegrable _ _, by rw [hP, hI], by rw [hP, hA], hA, ?_⟩
  rw [W04.proof_7 1 le_rfl, Pint_W04, hP, hI]

lemma expansion : ∃ C τ₀ : ℝ, 0 < τ₀ ∧ ∀ p : ℝ, 1 < p → Y04.tau p ≤ τ₀ →
    |(W04.perimeter p).toReal / Real.sqrt 2 -
        (4 + (π ^ 2 / 6) * Y04.tau p ^ 2 + Y04.z3 * Y04.tau p ^ 3
          - 33 / 16 * (π ^ 4 / 90) * Y04.tau p ^ 4
          - (3 / 4 * Z02.z5 + 7 / 4 * ((π ^ 2 / 6) * Y04.z3)) * Y04.tau p ^ 5)|
      ≤ C * Y04.tau p ^ 6 := by
  obtain ⟨C, τ₀, hτ₀, hC⟩ := Z02.proof_8 0
  have hpi := Real.pi_pos
  have hs : 0 < Real.sqrt 2 := by positivity
  refine ⟨π / Real.sqrt 2 * C, τ₀, hτ₀, fun p hp hτ => ?_⟩
  have hb := hC p hp hτ
  rw [(cauchy_true hp).2]
  set E := Y04.a 0 (Y04.h p) - (-1 : ℝ) ^ (0 : ℕ) * (Real.sqrt 2 / π) *
        (-(4 / (16 * ((0 : ℕ) : ℝ) ^ 2 - 1)) + (π ^ 2 / 6) * Y04.tau p ^ 2 + Y04.z3 * Y04.tau p ^ 3
          - (33 / 16 + 7 * ((0 : ℕ) : ℝ) ^ 2) * (π ^ 4 / 90) * Y04.tau p ^ 4
          + (-(3 / 4 + 4 * ((0 : ℕ) : ℝ) ^ 2) * Z02.z5
              - (7 / 4 + 4 * ((0 : ℕ) : ℝ) ^ 2) * ((π ^ 2 / 6) * Y04.z3)) * Y04.tau p ^ 5) with hE
  have key : π * Y04.a 0 (Y04.h p) / Real.sqrt 2 -
        (4 + (π ^ 2 / 6) * Y04.tau p ^ 2 + Y04.z3 * Y04.tau p ^ 3
          - 33 / 16 * (π ^ 4 / 90) * Y04.tau p ^ 4
          - (3 / 4 * Z02.z5 + 7 / 4 * ((π ^ 2 / 6) * Y04.z3)) * Y04.tau p ^ 5) = π / Real.sqrt 2 * E := by
    rw [hE]; simp only [Nat.cast_zero, pow_zero, one_mul]
    field_simp
    ring
  rw [key, abs_mul, abs_of_pos (by positivity), mul_assoc]
  exact mul_le_mul_of_nonneg_left hb (by positivity)

end V01

theorem V01.proof_1 : V01.stmt_1 := fun _ hp θ => V01.support_p hp θ

theorem V01.proof_2 : V01.stmt_2 := V01.support_one

theorem V01.proof_3 : V01.stmt_3 := fun _ hp =>
  ⟨V01.nA_strictAnti hp, V01.nA_image hp, fun _ ht => V01.nA_props hp ht⟩

theorem V01.proof_4 : V01.stmt_4 := fun _ hp => V01.quarter_cov hp

theorem V01.proof_5 : V01.stmt_5 := fun _ hp =>
  ⟨fun θ => V01.h_period θ, V01.int_period (Y04.h_cont hp) (fun θ => V01.h_period θ)⟩

theorem V01.proof_6 : V01.stmt_6 := fun _ hp => V01.cauchy_p hp

theorem V01.proof_7 : V01.stmt_7 := fun _ hp => V01.cauchy_true hp

theorem V01.proof_8 : V01.stmt_8 := V01.diamond

theorem V01.proof_9 : V01.stmt_9 := V01.expansion
