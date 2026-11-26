import Mathlib
import Zeta3.Wave2.D00_Dilog
/-!
# Z03 -- the special values used by Proposition `prop:nonlocal`

Target id: Z03.
Paper label: `prop:nonlocal` (main.tex 1697-1708), its last sentence (the four special values).
Detailed proof: paper/ieee/proofs/sec_g4_nonlocal.tex (Lemma `g4:lem:values`, and the proof of
Proposition 8), which uses Lemma `g4:lem:NK` of paper/ieee/proofs/sec_g4_hexagon.tex (116-124)
and the Nielsen function of the proof sketch of `thm:hexagon` (main.tex 1720-1730).

Quoted verbatim (main.tex 1697-1708, Proposition prop:nonlocal):

  Let $A(x)=\sum_{n\ge1}H_nx^n/n^2$ and let $\zeta_3[\,\cdot\,]$ be the coefficient of
  $\zeta(3)$ in the values listed below. For $r\in\{2,\frac12\}$ the $\zeta(3)$ of $P'''(1)$
  beyond $Z_e$ is $N_e=\Phi_A\,\zeta_3[A(1-r)]+\Phi_L\,\zeta_3[\Li_3(1-r)]$, and $N_e=0$ for
  $r=1$, where, with $\gamma_e=Bh^2/L^3$, ...
  Here $A(-1)=-\frac58\zeta(3)$, $\Li_3(-1)=-\frac34\zeta(3)$,
  $A(\frac12)=\zeta(3)-\frac{\pi^2\log2}{12}$,
  $\Li_3(\frac12)=\frac78\zeta(3)-\frac{\pi^2\log2}{12}+\frac{\log^32}6$.

Quoted verbatim (main.tex, proof sketch of thm:hexagon):

  $A(x)$ means the analytic continuation $\Li_3(x)+S_{1,2}(x)$, with
  $S_{1,2}(x)=\frac12\int_0^x\log^2(1-t)\,\frac{dt}t$, which equals the series for $|x|\le1$
  ... the classical Nielsen--K\"olbig identity~\cite{Kolbig86}, valid for $x<1$
  (every term is real for $x<1$),
  6A(x)=12\Li_3(x)+6\Li_3\bigl(\tfrac x{x-1}\bigr)-6\Li_2(x)\log(1-x)-\log^3(1-x),

Quoted verbatim (supplement, Lemma g4:lem:values):

  (a) $\Li_3(-1)=-\frac34\zeta(3)$, $\Li_2(-1)=-\frac12\zeta(2)$,
  $\Li_2(\frac12)=\frac12\zeta(2)-\frac12\log^22$, and
  $\Li_3(\frac12)=\frac78\zeta(3)-\frac{\pi^2\log2}{12}+\frac{\log^32}6$.
  (b) $A(-1)=-\frac58\zeta(3)$ and $A(\frac12)=\zeta(3)-\frac{\pi^2\log2}{12}$.

  and from its proof: "The fourth is Landen's trilogarithm identity, for $0<x<1$,
  $\Li_3(x)+\Li_3(1-x)+\Li_3(1-\frac1x)=\zeta(3)+\frac16\log^3x+\zeta(2)\log x
  -\frac12\log^2x\log(1-x)$"; "The series for $A$ converges absolutely on $|x|\le1$".

Quoted verbatim (supplement, proof of Lemma g4:lem:NK): the derivative of
  $2\mathcal L_3+\Li_3(u)+\mathcal L_1\mathcal L_2+\frac16\mathcal L_1^3$ is
  $\frac{\mathcal L_2}x+\frac{\mathcal L_1^2}{2x}$, which is $\mathcal A'$
  (with $\mathcal L_1=-\log(1-x)$, $\mathcal L_2=\Li_2(x)$, $\mathcal L_3=\Li_3(x)$).

Quoted verbatim (supplement, proof of Proposition 8): "The coefficients are
  $c_1/L^5=2\Phi_A+\Phi_L$ for $\Li_3(x)$ and $c_2/L^5=\Phi_A$ for $\Li_3(\frac x{x-1})$".

## Conventions

* `ζ(3)` is `z3 := (riemannZeta 3).re`; stmt_1 also says `riemannZeta 3 = ↑z3` and
  `HasSum (fun n => 1/(n+1)^3) z3` (proof via `zeta_nat_eq_tsum_of_gt_one` /
  `hasSum_zeta_nat`-style lemmas). `ζ(2)` is written `π²/6`.
* `H_n` is Mathlib's `harmonic n : ℚ` cast to `ℝ`; the series are indexed from `n = 0`
  with `n+1` in place of the paper's `n ≥ 1`.
* `Li₂ = D00.Li2` (imported, `-∫_0^x log(1-t)/t dt`). Its values `Li₂(-1)`, `Li₂(½)` of
  Lemma g4:lem:values (a) and Landen's dilogarithm identity (the first identity of g4:lem:NK)
  are already stmt_8, stmt_14, stmt_10 of D00 and are not restated.
* `Li3 x := ∫_0^x Li₂(t)/t dt` (the same definition as the unfinished Y01, re-made here since
  Y01 cannot be imported); stmt_2, stmt_3 record that it is integrable for `x ≤ 1` and equals
  the series `Σ x^k/k³` for `|x| ≤ 1`. All arguments `x/(x-1)`, `1-1/x` used below are `< 1`.
* `S12 x := ½ ∫_0^x log²(1-t)/t dt` (the Nielsen function `S_{1,2}`); stmt_9 asserts the
  integrability of its integrand on `[0,x]` for every `x ≤ 1`.
* `Aterm x n = H_{n+1} x^{n+1}/(n+1)²` and `A x := ∑' n, Aterm x n`; `A` is only used for
  `|x| ≤ 1`, where stmt_7 gives absolute summability, so the `tsum` is never the junk `0`.
  The four final values are stated as `HasSum` of the series themselves, not through `A`.

## Statements

* stmt_1: `riemannZeta 3 = ↑z3` and `Σ_{n≥1} 1/n³ = z3`.
* stmt_2: `t ↦ Li₂(t)/t` is integrable on `[0,x]` for every `x ≤ 1`.
* stmt_3: `Σ_{n≥1} x^n/n³ = Li3 x` for `|x| ≤ 1`.
* stmt_4: `Li3 1 = ζ(3)` and `Li3 (-1) = -¾ ζ(3)` (Lemma g4:lem:values (a), first value).
* stmt_5: Landen's trilogarithm identity on `0 < x < 1`.
* stmt_6: `Li3 (½) = ⅞ζ(3) - π² log 2/12 + log³2/6` (function value).
* stmt_7: the series of `A` converges absolutely for `|x| ≤ 1`.
* stmt_8: `A'(x) = (Li₂(x) + ½ log²(1-x))/x` for `|x| < 1`, `x ≠ 0` (proof of g4:lem:NK).
* stmt_9: the Nielsen function: integrability of `log²(1-t)/t` on `[0,x]` for `x ≤ 1`, and
  `Σ H_n x^n/n² = Li₃(x) + S_{1,2}(x)` for `|x| ≤ 1`.
* stmt_10: Nielsen--Kölbig on all of `x < 1`, in the form
  `S_{1,2}(x) = Li₃(x) + Li₃(x/(x-1)) - Li₂(x)log(1-x) - log³(1-x)/6`
  (equivalent to the displayed `6A = ...` with `A = Li₃ + S_{1,2}`).
* stmt_11: Lemma g4:lem:NK for the series on `-1 ≤ x < 1`:
  `Σ H_n x^n/n² = (12Li₃(x) + 6Li₃(x/(x-1)) - 6Li₂(x)log(1-x) - log³(1-x))/6`.
* stmt_12: the regrouping used in the proof of Proposition 8: for all reals `Φ_A, Φ_L` and
  `-1 ≤ x < 1`, `(2Φ_A+Φ_L)Li₃(x) + Φ_A Li₃(x/(x-1))
  = Φ_A A(x) + Φ_L Li₃(x) + Φ_A (Li₂(x)log(1-x) + log³(1-x)/6)`.
* stmt_13: `A(-1) = -⅝ ζ(3)` (HasSum).
* stmt_14: `Li₃(-1) = -¾ ζ(3)` (HasSum).
* stmt_15: `A(½) = ζ(3) - π² log 2/12` (HasSum).
* stmt_16: `Li₃(½) = ⅞ζ(3) - π² log 2/12 + log³2/6` (HasSum).
* stmt_17: the last sentence of prop:nonlocal: the conjunction of stmt_13 to stmt_16.

## NOT COVERED

* The formula `N_e = Φ_A ζ₃[A(1-r)] + Φ_L ζ₃[Li₃(1-r)]` for `r ∈ {2, ½}` and `N_e = 0` for
  `r = 1`: `N_e` is the coefficient of `ζ(3)` in the computer-generated expression `E_e` of
  Theorem thm:hexagon (machine steps M4-M6), with `ζ(3), ζ(2), log 2` treated as formal
  symbols; it is not a concrete identity between real numbers and is not formalized. The
  concrete analytic input of that proof (the regrouping of the `Li₃` coefficients
  `c₁ = 2Φ_A+Φ_L`, `c₂ = Φ_A` into `A` and `Li₃`) is stmt_12; the `r = 1` case only uses
  `Li₃(0) = Li₂(0) = log 1 = 0`, which is not worth a separate statement.
* The coefficients `Φ_A`, `Φ_L`, `γ_e` and the edge data, the claim about `P'''(1)`, and the
  numerical example (`N_e = 159√5/500`): not part of this target.
* The values `Li₂(-1) = -ζ(2)/2`, `Li₂(½) = ζ(2)/2 - ½log²2` (in D00).
-/
-- BEGIN STATEMENTS
namespace Z03

open MeasureTheory Real

/-- `ζ(3)` as a real number: the real part of Mathlib's `riemannZeta 3`. -/
noncomputable def z3 : ℝ := (riemannZeta 3).re

/-- The real trilogarithm `Li₃(x) = ∫_0^x Li₂(t)/t dt` (genuine on `(-∞,1]`). -/
noncomputable def Li3 (x : ℝ) : ℝ := ∫ t in (0 : ℝ)..x, D00.Li2 t / t

/-- The Nielsen function `S_{1,2}(x) = ½ ∫_0^x log²(1-t)/t dt`. -/
noncomputable def S12 (x : ℝ) : ℝ := (1 / 2) * ∫ t in (0 : ℝ)..x, Real.log (1 - t) ^ 2 / t

/-- The `n`-th term (from `n = 0`) of `A(x) = Σ_{n≥1} H_n x^n/n²`. -/
noncomputable def Aterm (x : ℝ) (n : ℕ) : ℝ :=
  (harmonic (n + 1) : ℝ) * x ^ (n + 1) / ((n : ℝ) + 1) ^ 2

/-- `A(x) = Σ_{n≥1} H_n x^n/n²` (used only for `|x| ≤ 1`, where the series converges). -/
noncomputable def A (x : ℝ) : ℝ := ∑' n : ℕ, Aterm x n

/-- The `n`-th term (from `n = 0`) of `Li₃(x) = Σ_{n≥1} x^n/n³`. -/
noncomputable def Li3term (x : ℝ) (n : ℕ) : ℝ := x ^ (n + 1) / ((n : ℝ) + 1) ^ 3

/-- `ζ(3)`: `riemannZeta 3` is the real number `z3`, and `z3 = Σ_{n≥1} 1/n³`. -/
def stmt_1 : Prop :=
  riemannZeta 3 = ((z3 : ℝ) : ℂ) ∧ HasSum (fun n : ℕ => 1 / ((n : ℝ) + 1) ^ 3) z3

/-- The integrand `Li₂(t)/t` of `Li3` is integrable on `[0,x]` for every `x ≤ 1`. -/
def stmt_2 : Prop :=
  ∀ x : ℝ, x ≤ 1 → IntervalIntegrable (fun t : ℝ => D00.Li2 t / t) volume 0 x

/-- "$\Li_3(x)=\sum_{n\ge1}x^n/n^3$": the integral `Li3` is the series for `|x| ≤ 1`. -/
def stmt_3 : Prop :=
  ∀ x : ℝ, |x| ≤ 1 → HasSum (Li3term x) (Li3 x)

/-- Lemma g4:lem:values (a), first value, with `ζ(3) = Li₃(1)`:
`Li₃(1) = ζ(3)` and `Li₃(-1) = -¾ζ(3)`. -/
def stmt_4 : Prop := Li3 1 = z3 ∧ Li3 (-1) = -(3 / 4) * z3

/-- Landen's trilogarithm identity (proof of g4:lem:values), for `0 < x < 1`. -/
def stmt_5 : Prop :=
  ∀ x : ℝ, 0 < x → x < 1 →
    Li3 x + Li3 (1 - x) + Li3 (1 - 1 / x)
      = z3 + (1 / 6) * Real.log x ^ 3 + π ^ 2 / 6 * Real.log x
        - (1 / 2) * Real.log x ^ 2 * Real.log (1 - x)

/-- Lemma g4:lem:values (a), fourth value: `Li₃(½) = ⅞ζ(3) - π² log 2/12 + log³2/6`. -/
def stmt_6 : Prop :=
  Li3 (1 / 2) = (7 / 8) * z3 - π ^ 2 * Real.log 2 / 12 + Real.log 2 ^ 3 / 6

/-- "The series for $A$ converges absolutely on $|x|\le1$". -/
def stmt_7 : Prop :=
  ∀ x : ℝ, |x| ≤ 1 → Summable (fun n : ℕ => |Aterm x n|)

/-- Proof of g4:lem:NK: `A'(x) = (Li₂(x) + ½ log²(1-x))/x` for `|x| < 1`, `x ≠ 0`. -/
def stmt_8 : Prop :=
  ∀ x : ℝ, |x| < 1 → x ≠ 0 →
    HasDerivAt A ((D00.Li2 x + (1 / 2) * Real.log (1 - x) ^ 2) / x) x

/-- Proof sketch of thm:hexagon: the Nielsen integrand is integrable on `[0,x]` for `x ≤ 1`,
and `S_{1,2}` "equals the series for $|x|\le1$", i.e. `Σ H_n x^n/n² = Li₃(x) + S_{1,2}(x)`. -/
def stmt_9 : Prop :=
  (∀ x : ℝ, x ≤ 1 →
      IntervalIntegrable (fun t : ℝ => Real.log (1 - t) ^ 2 / t) volume 0 x) ∧
    ∀ x : ℝ, |x| ≤ 1 → HasSum (Aterm x) (Li3 x + S12 x)

/-- The Nielsen--Kölbig identity "valid for $x<1$" for the continuation
`A = Li₃ + S_{1,2}`: `S_{1,2}(x) = Li₃(x) + Li₃(x/(x-1)) - Li₂(x)log(1-x) - log³(1-x)/6`. -/
def stmt_10 : Prop :=
  ∀ x : ℝ, x < 1 →
    S12 x = Li3 x + Li3 (x / (x - 1)) - D00.Li2 x * Real.log (1 - x)
      - Real.log (1 - x) ^ 3 / 6

/-- Lemma g4:lem:NK for the series, on `-1 ≤ x < 1`:
`6A(x) = 12Li₃(x) + 6Li₃(x/(x-1)) - 6Li₂(x)log(1-x) - log³(1-x)`. -/
def stmt_11 : Prop :=
  ∀ x : ℝ, -1 ≤ x → x < 1 →
    HasSum (Aterm x)
      ((12 * Li3 x + 6 * Li3 (x / (x - 1)) - 6 * D00.Li2 x * Real.log (1 - x)
        - Real.log (1 - x) ^ 3) / 6)

/-- Proof of Proposition 8: with the coefficients `c₁ = 2Φ_A+Φ_L` of `Li₃(x)` and
`c₂ = Φ_A` of `Li₃(x/(x-1))`, the trilogarithm part regroups as `Φ_A A(x) + Φ_L Li₃(x)` plus
terms in `Li₂(x)`, `log(1-x)` only, for `-1 ≤ x < 1`. -/
def stmt_12 : Prop :=
  ∀ PhiA PhiL x : ℝ, -1 ≤ x → x < 1 →
    (2 * PhiA + PhiL) * Li3 x + PhiA * Li3 (x / (x - 1))
      = PhiA * A x + PhiL * Li3 x
        + PhiA * (D00.Li2 x * Real.log (1 - x) + Real.log (1 - x) ^ 3 / 6)

/-- prop:nonlocal: "$A(-1)=-\frac58\zeta(3)$". -/
def stmt_13 : Prop := HasSum (Aterm (-1)) (-(5 / 8) * z3)

/-- prop:nonlocal: "$\Li_3(-1)=-\frac34\zeta(3)$". -/
def stmt_14 : Prop := HasSum (Li3term (-1)) (-(3 / 4) * z3)

/-- prop:nonlocal: "$A(\frac12)=\zeta(3)-\frac{\pi^2\log2}{12}$". -/
def stmt_15 : Prop := HasSum (Aterm (1 / 2)) (z3 - π ^ 2 * Real.log 2 / 12)

/-- prop:nonlocal: "$\Li_3(\frac12)=\frac78\zeta(3)-\frac{\pi^2\log2}{12}+\frac{\log^32}6$". -/
def stmt_16 : Prop :=
  HasSum (Li3term (1 / 2)) ((7 / 8) * z3 - π ^ 2 * Real.log 2 / 12 + Real.log 2 ^ 3 / 6)

/-- The last sentence of prop:nonlocal: all four special values. -/
def stmt_17 : Prop :=
  HasSum (Aterm (-1)) (-(5 / 8) * z3) ∧
    HasSum (Li3term (-1)) (-(3 / 4) * z3) ∧
    HasSum (Aterm (1 / 2)) (z3 - π ^ 2 * Real.log 2 / 12) ∧
    HasSum (Li3term (1 / 2)) ((7 / 8) * z3 - π ^ 2 * Real.log 2 / 12 + Real.log 2 ^ 3 / 6)

end Z03
-- END STATEMENTS

namespace Z03

open MeasureTheory Real Filter Topology Set

/-! ### `ζ(3)` -/

lemma zeta3_eq_tsum : riemannZeta 3 = ((∑' n : ℕ, 1 / ((n : ℝ)) ^ 3 : ℝ) : ℂ) := by
  rw [show (3 : ℂ) = ((3 : ℕ) : ℂ) by norm_num, zeta_nat_eq_tsum_of_gt_one (by norm_num),
    Complex.ofReal_tsum]
  push_cast; rfl

lemma zeta3_eq : riemannZeta 3 = ((z3 : ℝ) : ℂ) := by
  rw [z3, zeta3_eq_tsum, Complex.ofReal_re]

lemma hasSum_z3 : HasSum (fun n : ℕ => 1 / ((n : ℝ) + 1) ^ 3) z3 := by
  have h2 : z3 = ∑' n : ℕ, 1 / ((n : ℝ)) ^ 3 := by
    rw [z3, zeta3_eq_tsum, Complex.ofReal_re]
  have hs : Summable (fun n : ℕ => 1 / ((n : ℝ)) ^ 3) :=
    (Real.summable_one_div_nat_pow (p := 3)).2 (by norm_num)
  rw [h2, hs.tsum_eq_zero_add]
  simp only [CharP.cast_eq_zero, ne_eq, OfNat.ofNat_ne_zero, not_false_eq_true, zero_pow,
    div_zero, zero_add, Nat.cast_add, Nat.cast_one]
  convert ((summable_nat_add_iff 1).2 hs).hasSum using 1 <;> simp

/-! ### The trilogarithm as an integral -/

/-- The continuous extension of `Li₂(t)/t` (value `1` at `0`). -/
noncomputable def H (t : ℝ) : ℝ := if t = 0 then 1 else D00.Li2 t / t

lemma H_of_ne {t : ℝ} (h : t ≠ 0) : H t = D00.Li2 t / t := by simp [H, h]

lemma continuousAt_H_zero : ContinuousAt H 0 := by
  rw [← continuousWithinAt_compl_self, ContinuousWithinAt, show H 0 = 1 by simp [H] ]
  have hd := D00.hasDerivAt_Li2_zero
  rw [hasDerivAt_iff_tendsto_slope] at hd
  apply hd.congr'
  filter_upwards [self_mem_nhdsWithin] with t ht
  rw [H_of_ne ht, slope_def_field]
  simp

lemma continuous_H : Continuous H := by
  rw [continuous_iff_continuousAt]
  intro t
  rcases eq_or_ne t 0 with rfl | h0
  · exact continuousAt_H_zero
  · have hc : ContinuousAt (fun s => D00.Li2 s / s) t :=
      D00.continuous_Li2.continuousAt.div continuousAt_id h0
    apply hc.congr
    filter_upwards [isOpen_ne.mem_nhds h0] with s hs
    exact (H_of_ne hs).symm

lemma ae_eq_H (a b : ℝ) :
    (fun t : ℝ => D00.Li2 t / t) =ᵐ[volume.restrict (Set.uIoc a b)] H := by
  apply MeasureTheory.ae_restrict_of_ae
  filter_upwards [(Set.countable_singleton (0:ℝ)).ae_notMem volume] with t ht
  exact (H_of_ne ht).symm

lemma Li3_eq_H (x : ℝ) : Li3 x = ∫ t in (0 : ℝ)..x, H t := by
  unfold Li3
  apply intervalIntegral.integral_congr_ae
  filter_upwards [(Set.countable_singleton (0:ℝ)).ae_notMem volume] with t ht _
  exact (H_of_ne ht).symm

lemma intervalIntegrable_Li2_div (a b : ℝ) :
    IntervalIntegrable (fun t : ℝ => D00.Li2 t / t) volume a b :=
  (continuous_H.intervalIntegrable a b).congr_ae (ae_eq_H a b).symm

lemma hasDerivAt_Li3 (x : ℝ) : HasDerivAt Li3 (H x) x := by
  have := (continuous_H.integral_hasStrictDerivAt 0 x).hasDerivAt
  convert this using 1
  funext y
  exact Li3_eq_H y

lemma continuous_Li3 : Continuous Li3 :=
  continuous_iff_continuousAt.2 fun x => (hasDerivAt_Li3 x).continuousAt

@[simp] lemma Li3_zero : Li3 0 = 0 := by simp [Li3]

lemma hasDerivAt_Li3_comp {f : ℝ → ℝ} {f' x : ℝ} (hf : HasDerivAt f f' x) :
    HasDerivAt (fun y => Li3 (f y)) (H (f x) * f') x :=
  (hasDerivAt_Li3 (f x)).comp x hf

/-! ### Power series with derivative -/

lemma summable_geom_mul {r : ℝ} (hr0 : 0 ≤ r) (hr1 : r < 1) :
    Summable (fun n : ℕ => ((n : ℝ) + 1) * r ^ n) := by
  have h1 := summable_pow_mul_geometric_of_norm_lt_one 1
    (by rwa [Real.norm_eq_abs, abs_of_nonneg hr0] : ‖r‖ < 1)
  have h2 := summable_geometric_of_lt_one hr0 hr1
  simpa [add_mul] using h1.add h2

lemma hasDerivAt_powser (a : ℕ → ℝ) (ha : ∀ n, |a n| ≤ (n : ℝ) + 1) {y : ℝ} (hy : |y| < 1) :
    HasDerivAt (fun z => ∑' n : ℕ, a n * z ^ (n + 1) / ((n : ℝ) + 1))
      (∑' n : ℕ, a n * y ^ n) y := by
  set r := (|y| + 1) / 2 with hr
  have hr1 : r < 1 := by rw [hr]; linarith
  have hr0 : 0 ≤ r := by rw [hr]; positivity
  have hrp : 0 < r := by rw [hr]; positivity
  have ht : y ∈ Ioo (-r) r := ⟨by linarith [neg_abs_le y], by linarith [le_abs_self y]⟩
  exact hasDerivAt_tsum_of_isPreconnected (t := Ioo (-r) r) (y₀ := 0)
    (g := fun (n : ℕ) (z : ℝ) => a n * z ^ (n + 1) / ((n : ℝ) + 1))
    (g' := fun (n : ℕ) (z : ℝ) => a n * z ^ n) (summable_geom_mul hr0 hr1) isOpen_Ioo
    isPreconnected_Ioo
    (by
      intro n z _
      have := ((hasDerivAt_pow (n + 1) z).const_mul (a n)).div_const ((n : ℝ) + 1)
      refine this.congr_deriv ?_
      push_cast
      field_simp)
    (by
      intro n z hz
      rw [Real.norm_eq_abs, abs_mul, abs_pow]
      have hz' : |z| ≤ r := abs_le.mpr ⟨hz.1.le, hz.2.le⟩
      exact mul_le_mul (ha n) (pow_le_pow_left₀ (abs_nonneg z) hz' n) (by positivity)
        (by positivity))
    ⟨by linarith, by linarith⟩
    (by simp)
    ht

/-! ### `Li3` equals its series on `[-1,1]` -/

noncomputable def T (x : ℝ) : ℝ := ∑' n : ℕ, Li3term x n

lemma summable_inv_cube : Summable (fun n : ℕ => 1 / ((n : ℝ) + 1) ^ 3) := hasSum_z3.summable

lemma Li3term_bound {x : ℝ} (hx : |x| ≤ 1) (n : ℕ) :
    ‖Li3term x n‖ ≤ 1 / ((n : ℝ) + 1) ^ 3 := by
  rw [Li3term, Real.norm_eq_abs, abs_div, abs_pow,
    abs_of_pos (by positivity : (0:ℝ) < ((n:ℝ) + 1) ^ 3)]
  apply div_le_div_of_nonneg_right _ (by positivity)
  exact pow_le_one₀ (abs_nonneg x) hx

lemma summable_Li3term {x : ℝ} (hx : |x| ≤ 1) : Summable (Li3term x) :=
  Summable.of_norm_bounded summable_inv_cube (Li3term_bound hx)

lemma continuousOn_T : ContinuousOn T (Icc (-1) 1) := by
  apply continuousOn_tsum (fun n => ((continuous_pow (n + 1)).continuousOn.div_const _))
    summable_inv_cube
  intro n x hx
  exact Li3term_bound (abs_le.mpr hx) n

lemma hasDerivAt_T {y : ℝ} (hy : |y| < 1) (h0 : y ≠ 0) : HasDerivAt T (H y) y := by
  have hT : T = fun z => ∑' n : ℕ, (1 / ((n : ℝ) + 1) ^ 2) * z ^ (n + 1) / ((n : ℝ) + 1) := by
    funext z
    unfold T Li3term
    congr 1
    funext n
    field_simp
  rw [hT]
  have hd := hasDerivAt_powser (fun n : ℕ => 1 / ((n : ℝ) + 1) ^ 2) (by
    intro n
    rw [abs_of_pos (by positivity)]
    have h1 : (1:ℝ) ≤ ((n:ℝ) + 1) ^ 2 := one_le_pow₀ (by linarith [(Nat.cast_nonneg n : (0:ℝ) ≤ n)])
    calc 1 / ((n:ℝ) + 1) ^ 2 ≤ 1 := by rw [div_le_one (by positivity)]; exact h1
      _ ≤ (n:ℝ) + 1 := by linarith [(Nat.cast_nonneg n : (0:ℝ) ≤ n)]) hy
  have hsum : HasSum (fun n : ℕ => (1 / ((n : ℝ) + 1) ^ 2) * y ^ n) (H y) := by
    have := (D00.hasSum_Li2 hy.le).div_const y
    rw [H_of_ne h0]
    convert this using 1
    funext n
    field_simp
    ring
  rw [hsum.tsum_eq] at hd
  exact hd

lemma Li3_eq_T {x : ℝ} (hx : |x| ≤ 1) : Li3 x = T x := by
  have hx' := abs_le.mp hx
  have key : (fun z => Li3 z - T z) 0 = (fun z => Li3 z - T z) x := by
    apply D00.eq_of_deriv_zero' (f := fun z => Li3 z - T z)
    · apply continuous_Li3.continuousOn.sub
      apply continuousOn_T.mono
      intro z hz
      rw [mem_uIcc] at hz
      constructor <;> rcases hz with h | h <;> linarith [h.1, h.2]
    · intro z hz
      obtain ⟨hz0, hzx⟩ := D00.mem_uIoo_zero hz
      have hz1 : |z| < 1 := lt_of_lt_of_le hzx hx
      have := (hasDerivAt_Li3 z).sub (hasDerivAt_T hz1 hz0)
      rw [sub_self] at this
      exact this
  have hT0 : T 0 = 0 := by simp [T, Li3term]
  simp only [hT0, Li3_zero, sub_zero] at key
  linarith

lemma hasSum_Li3 {x : ℝ} (hx : |x| ≤ 1) : HasSum (Li3term x) (Li3 x) := by
  rw [Li3_eq_T hx, T]
  exact (summable_Li3term hx).hasSum

lemma Li3_one : Li3 1 = z3 := by
  have h := hasSum_Li3 (x := 1) (by norm_num)
  have h2 : Li3term 1 = fun n : ℕ => 1 / ((n : ℝ) + 1) ^ 3 := by
    funext n; simp [Li3term]
  rw [h2] at h
  exact h.unique hasSum_z3

lemma Li3_duplication {x : ℝ} (hx : |x| ≤ 1) :
    Li3 x + Li3 (-x) = (1 / 4) * Li3 (x ^ 2) := by
  have key := D00.eq_of_deriv_zero' (f := fun z => Li3 z + Li3 (-z) - (1/4) * Li3 (z^2))
    (a := 0) (b := x)
    ((continuous_Li3.add (continuous_Li3.comp continuous_neg)).sub
      (continuous_const.mul (continuous_Li3.comp (continuous_pow 2)))).continuousOn
    (by
      intro z hz
      obtain ⟨hz0, hzx⟩ := D00.mem_uIoo_zero hz
      have hz1 : |z| ≤ 1 := (lt_of_lt_of_le hzx hx).le
      have hsq0 : z ^ 2 ≠ 0 := pow_ne_zero 2 hz0
      have d1 := hasDerivAt_Li3 z
      have d2 := hasDerivAt_Li3_comp (hasDerivAt_neg z)
      have d3 := hasDerivAt_Li3_comp (hasDerivAt_pow 2 z)
      have := (d1.add d2).sub (d3.const_mul (1/4))
      refine this.congr_deriv ?_
      have e := D00.duplication hz1
      have e' : D00.Li2 (-z) = (1 / 2) * D00.Li2 (z ^ 2) - D00.Li2 z := by linarith
      rw [H_of_ne hz0, H_of_ne (neg_ne_zero.mpr hz0), H_of_ne hsq0, e']
      field_simp
      ring)
  simp at key
  linarith

lemma Li3_neg_one : Li3 (-1) = -(3 / 4) * z3 := by
  have := Li3_duplication (x := 1) (by norm_num)
  rw [one_pow, Li3_one] at this
  linarith

/-! ### Landen's trilogarithm identity -/

/-- `log z * log (1-z)` is continuous away from `0` (including at `z = 1`). -/
lemma continuousAt_log_mul_log {z : ℝ} (hz : z ≠ 0) :
    ContinuousAt (fun w => Real.log w * Real.log (1 - w)) z := by
  have hK : ContinuousAt (fun w : ℝ => (1 - w) * Real.log (1 - w) * D00.G (1 - w)) z := by
    have h1 : ContinuousAt (fun w : ℝ => (1 - w) * Real.log (1 - w)) z :=
      continuous_mul_log.continuousAt.comp (continuousAt_const.sub continuousAt_id)
    have h2 : ContinuousAt (fun w : ℝ => D00.G (1 - w)) z :=
      (D00.continuousAt_G (by intro h; apply hz; linarith)).comp
        (continuousAt_const.sub continuousAt_id)
    exact h1.mul h2
  apply hK.congr
  apply Filter.Eventually.of_forall
  intro w
  rcases eq_or_ne (1 - w) 0 with h | h
  · have hw : w = 1 := by linarith
    subst hw
    simp
  · simp only
    rw [D00.G_of_ne h, sub_sub_cancel]
    field_simp

lemma landen3_deriv {z : ℝ} (hzp : 0 < z) (hz1 : z < 1) :
    HasDerivAt (fun z => Li3 z + Li3 (1 - z) + Li3 (1 - z⁻¹)
      - ((1 / 6) * Real.log z ^ 3 + π ^ 2 / 6 * Real.log z
        - (1 / 2) * (Real.log z ^ 2 * Real.log (1 - z)))) 0 z := by
  have hz0 : z ≠ 0 := hzp.ne'
  have hz2 : 1 - z ≠ 0 := by linarith
  have hz3 : 1 - z⁻¹ ≠ 0 := by
    rw [sub_ne_zero]; intro h
    have : z = 1 := by
      have := congrArg (·⁻¹) h; simpa using this.symm
    linarith
  have d1 := hasDerivAt_Li3 z
  have d2 := hasDerivAt_Li3_comp ((hasDerivAt_id' z).const_sub 1)
  have d3 := hasDerivAt_Li3_comp ((hasDerivAt_inv hz0).const_sub 1)
  have dl := Real.hasDerivAt_log hz0
  have dm := ((hasDerivAt_id' z).const_sub 1).log hz2
  have := ((d1.fun_add d2).fun_add d3).fun_sub
    ((((dl.fun_pow 3).const_mul (1 / 6)).fun_add (dl.const_mul (π ^ 2 / 6))).fun_sub
      (((dl.fun_pow 2).fun_mul dm).const_mul (1 / 2)))
  refine this.congr_deriv ?_
  have eL := D00.landen (x := 1 - z) (by linarith)
  have ew : (1 - z) / (1 - z - 1) = 1 - z⁻¹ := by
    rw [show 1 - z - 1 = -z by ring, div_neg, sub_div, div_self hz0, one_div]; ring
  rw [ew, sub_sub_cancel] at eL
  have eR := D00.reflection hzp hz1
  have e1 : D00.Li2 (1 - z⁻¹) = -(1 / 2) * Real.log z ^ 2 - D00.Li2 (1 - z) := by linarith
  have e2 : D00.Li2 (1 - z) = π ^ 2 / 6 - Real.log z * Real.log (1 - z) - D00.Li2 z := by
    linarith
  have hz3' : 1 - z⁻¹ = (z - 1) / z := by field_simp
  rw [H_of_ne hz0, H_of_ne hz2, H_of_ne hz3, e1, e2, hz3']
  have hz4 : z - 1 ≠ 0 := by linarith
  field_simp
  ring

lemma landen3 {x : ℝ} (h0 : 0 < x) (h1 : x < 1) :
    Li3 x + Li3 (1 - x) + Li3 (1 - x⁻¹)
      = z3 + (1 / 6) * Real.log x ^ 3 + π ^ 2 / 6 * Real.log x
        - (1 / 2) * Real.log x ^ 2 * Real.log (1 - x) := by
  have key := D00.eq_of_deriv_zero
    (f := fun z => Li3 z + Li3 (1 - z) + Li3 (1 - z⁻¹)
      - ((1 / 6) * Real.log z ^ 3 + π ^ 2 / 6 * Real.log z
        - (1 / 2) * (Real.log z ^ 2 * Real.log (1 - z)))) (a := x) (b := 1) h1.le
    (by
      intro z hz
      have hz0 : z ≠ 0 := by linarith [hz.1]
      apply ContinuousAt.continuousWithinAt
      have c1 : ContinuousAt (fun z : ℝ => Li3 (1 - z)) z :=
        continuous_Li3.continuousAt.comp (continuousAt_const.sub continuousAt_id)
      have c2 : ContinuousAt (fun z : ℝ => Li3 (1 - z⁻¹)) z :=
        continuous_Li3.continuousAt.comp (continuousAt_const.sub (continuousAt_inv₀ hz0))
      have c3 : ContinuousAt (fun z : ℝ => Real.log z) z := Real.continuousAt_log hz0
      have c4 : ContinuousAt (fun z : ℝ => Real.log z ^ 2 * Real.log (1 - z)) z := by
        apply (c3.mul (continuousAt_log_mul_log hz0)).congr
        apply Filter.Eventually.of_forall
        intro w; simp only [Pi.mul_apply]; ring
      exact ((continuous_Li3.continuousAt.add c1).add c2).sub
        (((continuousAt_const.mul (c3.pow 3)).add (continuousAt_const.mul c3)).sub
          (continuousAt_const.mul c4)))
    (fun z hz => landen3_deriv (by linarith [hz.1]) hz.2)
  simp only [inv_one, sub_self, Li3_zero, Real.log_one, Li3_one] at key
  norm_num at key
  linarith

/-! ### Harmonic numbers -/

lemma harmonic_nonneg' (n : ℕ) : (0:ℝ) ≤ (harmonic n : ℝ) := by
  rcases Nat.eq_zero_or_pos n with rfl | h
  · simp
  · exact_mod_cast (harmonic_pos h.ne').le

lemma harmonic_le_succ (n : ℕ) : (harmonic (n + 1) : ℝ) ≤ (n : ℝ) + 1 := by
  have h := harmonic_le_one_add_log (n + 1)
  have h2 : Real.log (((n + 1 : ℕ) : ℝ)) ≤ ((n + 1 : ℕ) : ℝ) - 1 :=
    Real.log_le_sub_one_of_pos (by positivity)
  push_cast at h h2
  linarith

lemma summable_HC {y : ℝ} (hy : |y| < 1) :
    Summable (fun n : ℕ => (harmonic (n + 1) : ℝ) * y ^ (n + 1)) := by
  apply Summable.of_norm_bounded (summable_geom_mul (abs_nonneg y) hy)
  intro n
  rw [Real.norm_eq_abs, abs_mul, abs_pow, abs_of_nonneg (harmonic_nonneg' _)]
  apply mul_le_mul (harmonic_le_succ n) _ (by positivity) (by positivity)
  rw [pow_succ]
  exact mul_le_of_le_one_right (by positivity) hy.le

lemma hasSum_C {y : ℝ} (hy : |y| < 1) :
    HasSum (fun n : ℕ => (harmonic (n + 1) : ℝ) * y ^ (n + 1)) (-Real.log (1 - y) / (1 - y)) := by
  have hs := (summable_HC hy).hasSum
  set c := ∑' n : ℕ, (harmonic (n + 1) : ℝ) * y ^ (n + 1)
  have h1 : HasSum (fun n : ℕ => (harmonic n : ℝ) * y ^ (n + 1)) (y * c) := by
    rw [← hasSum_nat_add_iff' 1]
    simp only [Finset.range_one, Finset.sum_singleton, harmonic_zero, Rat.cast_zero, zero_mul,
      sub_zero]
    convert hs.mul_left y using 1
    funext n; ring
  have h2 := hs.sub h1
  have h3 : HasSum (fun n : ℕ => y ^ (n + 1) / ((n : ℝ) + 1)) (c - y * c) := by
    convert h2 using 1
    funext n
    rw [harmonic_succ]; push_cast; ring
  have h4 := h3.unique (Real.hasSum_pow_div_log_of_abs_lt_one hy)
  have hy1 : 1 - y ≠ 0 := by linarith [le_abs_self y]
  rw [show -Real.log (1 - y) / (1 - y) = c by rw [← h4]; field_simp]
  exact hs

/-! ### `B(x) = Σ H_n x^n/n = Li₂(x) + ½ log²(1-x)` -/

noncomputable def B (x : ℝ) : ℝ := ∑' n : ℕ, (harmonic (n + 1) : ℝ) * x ^ (n + 1) / ((n : ℝ) + 1)

lemma abs_harmonic_le (n : ℕ) : |(harmonic (n + 1) : ℝ)| ≤ (n : ℝ) + 1 := by
  rw [abs_of_nonneg (harmonic_nonneg' _)]; exact harmonic_le_succ n

lemma hasDerivAt_B {y : ℝ} (hy : |y| < 1) (h0 : y ≠ 0) :
    HasDerivAt B (-Real.log (1 - y) / (y * (1 - y))) y := by
  have hd := hasDerivAt_powser (fun n => (harmonic (n + 1) : ℝ)) abs_harmonic_le hy
  have hy1 : 1 - y ≠ 0 := by linarith [le_abs_self y]
  have hsum : HasSum (fun n : ℕ => (harmonic (n + 1) : ℝ) * y ^ n)
      (-Real.log (1 - y) / (y * (1 - y))) := by
    have := (hasSum_C hy).div_const y
    convert this using 1
    · funext n; rw [pow_succ]; field_simp
    · field_simp
  rw [hsum.tsum_eq] at hd
  exact hd

lemma continuousAt_B {y : ℝ} (hy : |y| < 1) : ContinuousAt B y :=
  (hasDerivAt_powser (fun n => (harmonic (n + 1) : ℝ)) abs_harmonic_le hy).continuousAt

lemma abs_le_of_mem_uIcc {x z : ℝ} (hz : z ∈ uIcc 0 x) : |z| ≤ |x| := by
  rw [mem_uIcc] at hz
  rcases hz with ⟨h1, h2⟩ | ⟨h1, h2⟩
  · rw [abs_of_nonneg h1, abs_of_nonneg (h1.trans h2)]; exact h2
  · rw [abs_of_nonpos h2, abs_of_nonpos (h1.trans h2)]; linarith

lemma B_eq {x : ℝ} (hx : |x| < 1) : B x = D00.Li2 x + (1 / 2) * Real.log (1 - x) ^ 2 := by
  have key := D00.eq_of_deriv_zero'
    (f := fun z => B z - D00.Li2 z - (1 / 2) * Real.log (1 - z) ^ 2) (a := 0) (b := x)
    (by
      intro z hz
      have hz1 : |z| < 1 := lt_of_le_of_lt (abs_le_of_mem_uIcc hz) hx
      have hz2 : 1 - z ≠ 0 := by linarith [le_abs_self z]
      apply ContinuousAt.continuousWithinAt
      exact ((continuousAt_B hz1).sub D00.continuous_Li2.continuousAt).sub
        (continuousAt_const.mul (((continuousAt_const.sub continuousAt_id).log hz2).pow 2)))
    (by
      intro z hz
      obtain ⟨hz0, hzx⟩ := D00.mem_uIoo_zero hz
      have hz1 : |z| < 1 := lt_trans hzx hx
      have hz1' : z < 1 := lt_of_le_of_lt (le_abs_self z) hz1
      have hz2 : 1 - z ≠ 0 := by linarith
      have := ((hasDerivAt_B hz1 hz0).fun_sub (D00.hasDerivAt_Li2 hz1' hz0)).fun_sub
        (((((hasDerivAt_id' z).const_sub 1).log hz2).fun_pow 2).const_mul (1 / 2))
      refine this.congr_deriv ?_
      field_simp
      ring)
  have hB0 : B 0 = 0 := by simp [B]
  simp only [hB0, D00.Li2_zero, sub_zero, Real.log_one] at key
  norm_num at key
  linarith

lemma summable_B {y : ℝ} (hy : |y| < 1) :
    Summable (fun n : ℕ => (harmonic (n + 1) : ℝ) * y ^ (n + 1) / ((n : ℝ) + 1)) := by
  apply Summable.of_norm_bounded (summable_geometric_of_lt_one (abs_nonneg y) hy)
  intro n
  have hp : (0:ℝ) < (n:ℝ) + 1 := by positivity
  rw [Real.norm_eq_abs, abs_div, abs_mul, abs_pow, abs_of_nonneg (harmonic_nonneg' _),
    abs_of_pos hp, div_le_iff₀ hp]
  have : |y| ^ (n + 1) ≤ |y| ^ n := by
    rw [pow_succ]; exact mul_le_of_le_one_right (by positivity) hy.le
  calc (harmonic (n + 1) : ℝ) * |y| ^ (n + 1) ≤ ((n:ℝ) + 1) * |y| ^ n :=
        mul_le_mul (harmonic_le_succ n) this (by positivity) (by positivity)
    _ = |y| ^ n * ((n:ℝ) + 1) := by ring

/-! ### The series `A` -/

lemma A_eq_powser : A = fun z => ∑' n : ℕ,
    ((harmonic (n + 1) : ℝ) / ((n : ℝ) + 1)) * z ^ (n + 1) / ((n : ℝ) + 1) := by
  funext z; unfold A Aterm; congr 1; funext n; field_simp

lemma hasDerivAt_A {y : ℝ} (hy : |y| < 1) (h0 : y ≠ 0) :
    HasDerivAt A ((D00.Li2 y + (1 / 2) * Real.log (1 - y) ^ 2) / y) y := by
  rw [A_eq_powser]
  have hd := hasDerivAt_powser (fun n => (harmonic (n + 1) : ℝ) / ((n : ℝ) + 1)) (fun n => by
    have hp : (0:ℝ) < (n:ℝ) + 1 := by positivity
    rw [abs_of_nonneg (div_nonneg (harmonic_nonneg' _) hp.le)]
    have : (harmonic (n + 1) : ℝ) / ((n : ℝ) + 1) ≤ 1 := by
      rw [div_le_one hp]; exact harmonic_le_succ n
    linarith [(Nat.cast_nonneg n : (0:ℝ) ≤ n)]) hy
  have hsum : HasSum (fun n : ℕ => (harmonic (n + 1) : ℝ) / ((n : ℝ) + 1) * y ^ n)
      ((D00.Li2 y + (1 / 2) * Real.log (1 - y) ^ 2) / y) := by
    rw [← B_eq hy]
    have := (summable_B hy).hasSum.div_const y
    convert this using 1
    swap
    · rfl
    funext n
    rw [pow_succ]
    field_simp
  rw [hsum.tsum_eq] at hd
  exact hd

lemma harmonic_le_sqrt (n : ℕ) : (harmonic (n + 1) : ℝ) ≤ 3 * ((n : ℝ) + 1) ^ (1 / 2 : ℝ) := by
  have h := harmonic_le_one_add_log (n + 1)
  have hp : (0:ℝ) ≤ (n:ℝ) + 1 := by positivity
  have h2 := Real.log_le_rpow_div (x := ((n:ℝ) + 1)) hp (by norm_num : (0:ℝ) < 1 / 2)
  have h3 : (1:ℝ) ≤ ((n:ℝ) + 1) ^ (1 / 2 : ℝ) :=
    Real.one_le_rpow (by linarith [(Nat.cast_nonneg n : (0:ℝ) ≤ n)]) (by norm_num)
  push_cast at h
  linarith

lemma Aterm_bound {x : ℝ} (hx : |x| ≤ 1) (n : ℕ) :
    |Aterm x n| ≤ 3 * (((n : ℝ) + 1) ^ (3 / 2 : ℝ))⁻¹ := by
  have hp : (0:ℝ) < (n:ℝ) + 1 := by positivity
  set s := ((n : ℝ) + 1) ^ (1 / 2 : ℝ) with hs_def
  have hs : 0 < s := Real.rpow_pos_of_pos hp _
  have hss : s * s = (n : ℝ) + 1 := by
    rw [hs_def, ← Real.rpow_add hp]; norm_num
  have e : ((n : ℝ) + 1) ^ (3 / 2 : ℝ) = s * ((n : ℝ) + 1) := by
    rw [show (3 / 2 : ℝ) = 1 / 2 + 1 by norm_num, Real.rpow_add hp, Real.rpow_one]
  rw [Aterm, abs_div, abs_mul, abs_pow, abs_of_nonneg (harmonic_nonneg' _),
    abs_of_pos (by positivity : (0:ℝ) < ((n:ℝ) + 1) ^ 2), e,
    show ((n:ℝ) + 1) ^ 2 = s * s * ((n:ℝ) + 1) by rw [hss]; ring,
    div_le_iff₀ (by positivity)]
  have e2 : 3 * (s * ((n:ℝ) + 1))⁻¹ * (s * s * ((n:ℝ) + 1)) = 3 * s := by field_simp
  rw [e2]
  have hxp : |x| ^ (n + 1) ≤ 1 := pow_le_one₀ (abs_nonneg x) hx
  calc (harmonic (n + 1) : ℝ) * |x| ^ (n + 1) ≤ (harmonic (n + 1) : ℝ) * 1 :=
        mul_le_mul_of_nonneg_left hxp (harmonic_nonneg' _)
    _ ≤ 3 * s := by rw [mul_one]; exact harmonic_le_sqrt n

lemma summable_M : Summable (fun n : ℕ => 3 * (((n : ℝ) + 1) ^ (3 / 2 : ℝ))⁻¹) := by
  have h := (summable_nat_add_iff 1).mpr
    (Real.summable_nat_rpow_inv.mpr (by norm_num : (1:ℝ) < 3 / 2))
  exact (h.mul_left 3).congr (fun n => by push_cast; ring)

lemma summable_abs_Aterm {x : ℝ} (hx : |x| ≤ 1) : Summable (fun n : ℕ => |Aterm x n|) :=
  Summable.of_nonneg_of_le (fun _ => abs_nonneg _) (Aterm_bound hx) summable_M

lemma continuousOn_A : ContinuousOn A (Icc (-1) 1) := by
  apply continuousOn_tsum (fun n => by unfold Aterm; fun_prop) summable_M
  intro n x hx
  rw [Real.norm_eq_abs]
  exact Aterm_bound (abs_le.mpr hx) n

/-! ### The Nielsen function -/

noncomputable def Q (t : ℝ) : ℝ := Real.log (1 - t) ^ 2 / t

lemma Q_eq (t : ℝ) : Q t = Real.log (1 - t) * D00.G t := by
  rcases eq_or_ne t 0 with rfl | h
  · simp [Q, D00.G_zero]
  · rw [Q, D00.G_of_ne h]; ring

lemma continuousAt_Q {t : ℝ} (ht : t ≠ 1) : ContinuousAt Q t := by
  have : ContinuousAt (fun t => Real.log (1 - t) * D00.G t) t :=
    ((continuousAt_const.sub continuousAt_id).log (sub_ne_zero.mpr (Ne.symm ht))).mul
      (D00.continuousAt_G ht)
  exact this.congr (Eventually.of_forall fun s => (Q_eq s).symm)

lemma log_sq_bound {u : ℝ} (hu0 : 0 < u) (hu1 : u ≤ 1) :
    Real.log u ^ 2 ≤ 16 * u ^ (-1 / 2 : ℝ) := by
  have h1 := Real.neg_rpow_div_le_log hu0.le (by norm_num : (0:ℝ) < 1 / 4)
  have h2 : Real.log u ≤ 0 := Real.log_nonpos hu0.le hu1
  have h3 : u ^ (-(1 / 4 : ℝ)) * u ^ (-(1 / 4 : ℝ)) = u ^ (-1 / 2 : ℝ) := by
    rw [← Real.rpow_add hu0]; norm_num
  have h4 : 0 ≤ u ^ (-(1 / 4 : ℝ)) := Real.rpow_nonneg hu0.le _
  rw [← h3]
  set a := u ^ (-(1 / 4 : ℝ))
  have h5 : -(4 * a) ≤ Real.log u := by
    have : -a / (1 / 4) = -(4 * a) := by ring
    linarith
  nlinarith [mul_nonneg (by linarith : 0 ≤ 4 * a + Real.log u)
    (by linarith : 0 ≤ 4 * a - Real.log u)]

lemma intervalIntegrable_log_sq :
    IntervalIntegrable (fun u : ℝ => Real.log u ^ 2) volume 0 1 := by
  apply IntervalIntegrable.mono_fun' (g := fun u : ℝ => 16 * u ^ (-1 / 2 : ℝ))
  · exact (intervalIntegral.intervalIntegrable_rpow' (by norm_num)).const_mul 16
  · exact (Real.measurable_log.pow_const 2).aestronglyMeasurable
  · rw [Filter.EventuallyLE, ae_restrict_iff' measurableSet_uIoc]
    apply Eventually.of_forall
    intro u hu
    rw [uIoc_of_le zero_le_one] at hu
    rw [Real.norm_eq_abs, abs_of_nonneg (sq_nonneg _)]
    exact log_sq_bound hu.1 hu.2

lemma intervalIntegrable_Q_half_one : IntervalIntegrable Q volume (1 / 2) 1 := by
  have h0 : IntervalIntegrable (fun u : ℝ => Real.log u ^ 2) volume (1 - 1 / 2) (1 - 1) :=
    intervalIntegrable_log_sq.mono_set (uIcc_subset_uIcc
      (by rw [mem_uIcc]; left; constructor <;> norm_num)
      (by rw [mem_uIcc]; left; constructor <;> norm_num))
  have h1 : IntervalIntegrable (fun t : ℝ => Real.log (1 - t) ^ 2) volume (1 / 2) 1 := by
    have := h0.comp_sub_left 1
    convert this using 1 <;> norm_num
  have h2 := h1.mul_continuousOn (g := fun t : ℝ => t⁻¹) (by
    intro t ht
    rw [uIcc_of_le (by norm_num : (1 / 2 : ℝ) ≤ 1)] at ht
    exact (continuousAt_inv₀ (by linarith [ht.1] : t ≠ 0)).continuousWithinAt)
  convert h2 using 1
  funext t
  rw [Q, div_eq_mul_inv]

lemma intervalIntegrable_Q {a b : ℝ} (ha : a ≤ 1) (hb : b ≤ 1) :
    IntervalIntegrable Q volume a b := by
  have key : ∀ c : ℝ, c ≤ 1 → IntervalIntegrable Q volume c (1 / 2) := by
    intro c hc
    rcases le_total c (1 / 2) with h | h
    · apply ContinuousOn.intervalIntegrable
      intro t ht
      rw [uIcc_of_le h] at ht
      exact (continuousAt_Q (by linarith [ht.2])).continuousWithinAt
    · exact (intervalIntegrable_Q_half_one.mono_set (uIcc_subset_uIcc
        (by rw [mem_uIcc]; left; constructor <;> norm_num)
        (by rw [mem_uIcc]; left; exact ⟨h, hc⟩))).symm
  exact (key a ha).trans (key b hb).symm

lemma S12_def : S12 = fun x => (1 / 2) * ∫ t in (0 : ℝ)..x, Q t := rfl

@[simp] lemma S12_zero : S12 0 = 0 := by simp [S12]

lemma hasDerivAt_S12 {z : ℝ} (hz : z < 1) : HasDerivAt S12 ((1 / 2) * Q z) z := by
  have hm : StronglyMeasurableAtFilter Q (𝓝 z) volume := by
    have hco : ContinuousOn Q (Iio 1) := fun t ht =>
      (continuousAt_Q (ne_of_lt ht)).continuousWithinAt
    exact hco.stronglyMeasurableAtFilter isOpen_Iio z hz
  have := (intervalIntegral.integral_hasDerivAt_right
    (intervalIntegrable_Q (a := 0) (by norm_num) hz.le) hm (continuousAt_Q hz.ne)).const_mul (1 / 2 : ℝ)
  rw [S12_def]
  exact this

lemma continuousOn_S12 {x : ℝ} (hx : x ≤ 1) : ContinuousOn S12 (uIcc 0 x) := by
  rw [S12_def]
  exact continuousOn_const.mul (intervalIntegral.continuousOn_primitive_interval'
    (intervalIntegrable_Q (by norm_num) hx) left_mem_uIcc)

lemma A_eq {x : ℝ} (hx : |x| ≤ 1) : A x = Li3 x + S12 x := by
  have hx1 : x ≤ 1 := le_trans (le_abs_self x) hx
  have key := D00.eq_of_deriv_zero' (f := fun z => A z - Li3 z - S12 z) (a := 0) (b := x)
    (by
      apply ContinuousOn.sub (ContinuousOn.sub (continuousOn_A.mono ?_)
        continuous_Li3.continuousOn) (continuousOn_S12 hx1)
      intro z hz
      exact abs_le.mp ((abs_le_of_mem_uIcc hz).trans hx))
    (by
      intro z hz
      obtain ⟨hz0, hzx⟩ := D00.mem_uIoo_zero hz
      have hz1 : |z| < 1 := lt_of_lt_of_le hzx hx
      have hz1' : z < 1 := lt_of_le_of_lt (le_abs_self z) hz1
      have := ((hasDerivAt_A hz1 hz0).fun_sub (hasDerivAt_Li3 z)).fun_sub (hasDerivAt_S12 hz1')
      refine this.congr_deriv ?_
      rw [H_of_ne hz0, Q]
      field_simp
      ring)
  have hA0 : A 0 = 0 := by simp [A, Aterm]
  simp only [hA0, Li3_zero, S12_zero, sub_zero] at key
  linarith

/-! ### Nielsen--Kölbig -/

lemma NK {x : ℝ} (hx : x < 1) :
    S12 x = Li3 x + Li3 (x / (x - 1)) - D00.Li2 x * Real.log (1 - x)
      - Real.log (1 - x) ^ 3 / 6 := by
  have key := D00.eq_of_deriv_zero'
    (f := fun z => S12 z - Li3 z - Li3 (z / (z - 1)) + D00.Li2 z * Real.log (1 - z)
      + Real.log (1 - z) ^ 3 / 6) (a := 0) (b := x)
    (by
      intro z hz
      have hz' : z < 1 := by
        rw [mem_uIcc] at hz
        rcases hz with h | h <;> linarith [h.1, h.2]
      have hz1 : z - 1 ≠ 0 := by linarith
      have hz2 : 1 - z ≠ 0 := by linarith
      apply ContinuousAt.continuousWithinAt
      have c0 := (hasDerivAt_S12 hz').continuousAt
      have c1 : ContinuousAt (fun z : ℝ => Li3 (z / (z - 1))) z :=
        continuous_Li3.continuousAt.comp
          (continuousAt_id.div (continuousAt_id.sub continuousAt_const) hz1)
      have c2 : ContinuousAt (fun z : ℝ => Real.log (1 - z)) z :=
        (continuousAt_const.sub continuousAt_id).log hz2
      exact (((c0.sub continuous_Li3.continuousAt).sub c1).add
        (D00.continuous_Li2.continuousAt.mul c2)).add ((c2.pow 3).div_const 6))
    (by
      intro z hz
      obtain ⟨hz0, -⟩ := D00.mem_uIoo_zero hz
      have hz' : z < 1 := by
        simp only [uIoo, mem_Ioo] at hz
        rcases le_total 0 x with h | h
        · rw [min_eq_left h, max_eq_right h] at hz; linarith [hz.2]
        · rw [min_eq_right h, max_eq_left h] at hz; linarith [hz.2]
      have hz1 : z - 1 ≠ 0 := by linarith
      have hz2 : 1 - z ≠ 0 := by linarith
      have hq0 : z / (z - 1) ≠ 0 := div_ne_zero hz0 hz1
      have dS := hasDerivAt_S12 hz'
      have d1 := hasDerivAt_Li3 z
      have d2 := hasDerivAt_Li3_comp
        ((hasDerivAt_id' z).fun_div ((hasDerivAt_id' z).sub_const 1) hz1)
      have dL := ((hasDerivAt_id' z).const_sub 1).log hz2
      have d3 := (D00.hasDerivAt_Li2 hz' hz0).fun_mul dL
      have d4 := (dL.fun_pow 3).div_const 6
      have := (((dS.fun_sub d1).fun_sub d2).fun_add d3).fun_add d4
      refine this.congr_deriv ?_
      have eL := D00.landen hz'
      have e1 : D00.Li2 (z / (z - 1)) = -(1 / 2) * Real.log (1 - z) ^ 2 - D00.Li2 z := by
        linarith
      rw [H_of_ne hz0, H_of_ne hq0, e1, Q]
      field_simp
      ring)
  norm_num at key
  linarith

end Z03

theorem Z03.proof_1 : Z03.stmt_1 := ⟨Z03.zeta3_eq, Z03.hasSum_z3⟩
theorem Z03.proof_2 : Z03.stmt_2 := fun x _ => Z03.intervalIntegrable_Li2_div 0 x
theorem Z03.proof_3 : Z03.stmt_3 := fun _ hx => Z03.hasSum_Li3 hx
theorem Z03.proof_4 : Z03.stmt_4 := ⟨Z03.Li3_one, Z03.Li3_neg_one⟩
theorem Z03.proof_5 : Z03.stmt_5 := by
  intro x h0 h1
  rw [one_div]
  exact Z03.landen3 h0 h1
theorem Z03.proof_6 : Z03.stmt_6 := by
  have h := Z03.landen3 (x := 1 / 2) (by norm_num) (by norm_num)
  rw [show (1:ℝ) - 1 / 2 = 1 / 2 by norm_num, show (1:ℝ) - (1 / 2)⁻¹ = -1 by norm_num,
    Z03.Li3_neg_one, show Real.log (1 / 2) = -Real.log 2 by rw [one_div, Real.log_inv]] at h
  show Z03.Li3 (1 / 2) = _
  linear_combination h / 2
theorem Z03.proof_7 : Z03.stmt_7 := fun _ hx => Z03.summable_abs_Aterm hx
theorem Z03.proof_8 : Z03.stmt_8 := fun _ hx h0 => Z03.hasDerivAt_A hx h0
theorem Z03.proof_9 : Z03.stmt_9 := by
  refine ⟨fun x hx => Z03.intervalIntegrable_Q (by norm_num) hx, fun x hx => ?_⟩
  rw [← Z03.A_eq hx]
  have hs : Summable (Z03.Aterm x) :=
    Summable.of_norm (by simpa [Real.norm_eq_abs] using Z03.summable_abs_Aterm hx)
  exact hs.hasSum
theorem Z03.proof_10 : Z03.stmt_10 := fun _ hx => Z03.NK hx
theorem Z03.proof_11 : Z03.stmt_11 := by
  intro x h1 h2
  have hx : |x| ≤ 1 := abs_le.mpr ⟨h1, h2.le⟩
  have := Z03.proof_9.2 x hx
  convert this using 1
  rw [Z03.NK h2]; ring
theorem Z03.proof_12 : Z03.stmt_12 := by
  intro PhiA PhiL x h1 h2
  have h : Z03.A x = _ := (Z03.proof_11 x h1 h2).tsum_eq
  rw [h]; ring
theorem Z03.proof_13 : Z03.stmt_13 := by
  have h6 : Z03.Li3 (1 / 2) = (7 / 8) * Z03.z3 - Real.pi ^ 2 * Real.log 2 / 12 + Real.log 2 ^ 3 / 6 :=
    Z03.proof_6
  have h := Z03.proof_11 (-1) le_rfl (by norm_num)
  rw [show (-1:ℝ) / (-1 - 1) = 1 / 2 by norm_num, show (1:ℝ) - -1 = 2 by norm_num,
    Z03.Li3_neg_one, h6, D00.Li2_neg_one] at h
  show HasSum _ _
  convert h using 1
  ring
theorem Z03.proof_14 : Z03.stmt_14 := by
  have := Z03.hasSum_Li3 (x := -1) (by norm_num)
  rwa [Z03.Li3_neg_one] at this
theorem Z03.proof_15 : Z03.stmt_15 := by
  have h6 : Z03.Li3 (1 / 2) = (7 / 8) * Z03.z3 - Real.pi ^ 2 * Real.log 2 / 12 + Real.log 2 ^ 3 / 6 :=
    Z03.proof_6
  have h := Z03.proof_11 (1 / 2) (by norm_num) (by norm_num)
  rw [show (1 / 2 : ℝ) / (1 / 2 - 1) = -1 by norm_num, show (1:ℝ) - 1 / 2 = 1 / 2 by norm_num,
    Z03.Li3_neg_one, h6, D00.Li2_half,
    show Real.log (1 / 2) = -Real.log 2 by rw [one_div, Real.log_inv]] at h
  show HasSum _ _
  convert h using 1
  ring
theorem Z03.proof_16 : Z03.stmt_16 := by
  have h6 : Z03.Li3 (1 / 2) = (7 / 8) * Z03.z3 - Real.pi ^ 2 * Real.log 2 / 12 + Real.log 2 ^ 3 / 6 :=
    Z03.proof_6
  have := Z03.hasSum_Li3 (x := 1 / 2) (by norm_num [abs_of_pos])
  rwa [h6] at this
theorem Z03.proof_17 : Z03.stmt_17 :=
  ⟨Z03.proof_13, Z03.proof_14, Z03.proof_15, Z03.proof_16⟩
