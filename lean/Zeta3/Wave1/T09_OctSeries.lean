import Mathlib
/-!
# T09 — consistency of the two series in the octahedron-end theorem (paper label `thm:oct`)

Source: `paper/ieee/main.tex`, Section "The Octahedron End" (`sec:octahedron`), line 1973
(setup) and Theorem `thm:oct`, lines 1975–1990.

Setup (line 1973):

> Let's write $d=p-1$, as at the diamond, and $\tau=1-\frac1p=\frac d{1+d}$, and $[d^k]S$ for
> the coefficient of $d^k$ in $S(1+d)$.

Claim 1 of `thm:oct` (lines 1977–1980):

> As $d\to0$ from either side,
> \frac{S(1+d)}{4\sqrt3}=1+d+\Bigl(1-\frac{\pi^2}{9}\Bigr)d^2+\Bigl(1+\frac{\pi^2}{9}-\frac{16\zeta(3)}{9}\Bigr)d^3
> +\Bigl(1-\frac{2\pi^2}{9}+\frac{32\zeta(3)}{9}-\frac{16\pi^4}{405}\Bigr)d^4+O(d^5),

Claim 2 of `thm:oct` (lines 1982–1986):

> \log\frac{(2-p)\,S(p)}{4\sqrt3}=-\frac{\pi^2}{9}\tau^2-\frac{16\zeta(3)}{9}\tau^3
> -\frac{37\pi^4}{810}\tau^4-\frac{160\zeta(5)}{27}\tau^5
> +\Bigl(\frac{16}{81}\zeta(3)^2-\frac{572}{76545}\pi^6\Bigr)\tau^6+O(\tau^7),

What is formalized: the internal consistency of the two displays through order 4. Writing
`F(d) = S(1+d)/(4√3)` (so `2-p = 1-d`), `stmt_1` says: for ANY real function `F` that is
positive near `d = 0` and whose logarithm `log((1-d) F(d))` agrees with the paper's `τ`-series
through `τ^4` (error `O(τ^5)`, `τ = d/(1+d)`) as `d → 0` (two-sided), `F` agrees with the paper's
`d`-series through `d^4` (error `O(d^5)`). So exponentiating the `τ`-series and dividing by
`(2-p)` yields exactly the stated `d`-coefficients.

Conventions:
* `ζ(3)` is a universally quantified real `z3` (the claim is an identity of formal
  coefficients; it holds for every value, in particular for `ζ(3)`). `π` is `Real.pi`.
* The `τ^5`, `τ^6` terms (with `ζ(5)`, `ζ(3)^2`) are absorbed into the `O(τ^5)` hypothesis; they
  do not affect the `d`-series through `d^4`.
* Positivity of `F` near `0` is assumed so that `Real.log` is the true logarithm (the paper's
  `S` is a surface area, hence positive); no other hypothesis is added.
* "From either side" is the two-sided filter `𝓝 0`.

NOT COVERED: the actual analytic statements about the superellipsoid surface area `S(p)`
(the definition of `S` and the truth of either expansion), the weight/rationality claim for
every `[d^k]S`, the homogeneity of the `τ^k` coefficients, and the `τ^5`, `τ^6` coefficients.
-/
-- BEGIN STATEMENTS
namespace T09

/-- The paper's `τ`-series of `log((2-p)S(p)/(4√3))`, truncated after `τ^4`
(`z3` stands for `ζ(3)`). -/
noncomputable def logSeries (z3 τ : ℝ) : ℝ :=
  -(Real.pi ^ 2 / 9) * τ ^ 2 - (16 * z3 / 9) * τ ^ 3 - (37 * Real.pi ^ 4 / 810) * τ ^ 4

/-- The paper's `d`-series of `S(1+d)/(4√3)` through `d^4` (`z3` stands for `ζ(3)`). -/
noncomputable def dSeries (z3 d : ℝ) : ℝ :=
  1 + d + (1 - Real.pi ^ 2 / 9) * d ^ 2 + (1 + Real.pi ^ 2 / 9 - 16 * z3 / 9) * d ^ 3
    + (1 - 2 * Real.pi ^ 2 / 9 + 32 * z3 / 9 - 16 * Real.pi ^ 4 / 405) * d ^ 4

/-- thm:oct, lines 1973–1986: with `p = 1+d`, `τ = d/(1+d)` and `F(d) = S(1+d)/(4√3)`,
the `τ`-series of `log((2-p) F)` through `τ^4` implies the `d`-series of `F` through `d^4`. -/
def stmt_1 : Prop :=
  ∀ (z3 : ℝ) (F : ℝ → ℝ),
    (∀ᶠ d in nhds (0 : ℝ), 0 < F d) →
    (fun d : ℝ => Real.log ((1 - d) * F d) - logSeries z3 (d / (1 + d)))
        =O[nhds 0] (fun d : ℝ => (d / (1 + d)) ^ 5) →
    (fun d : ℝ => F d - dSeries z3 d) =O[nhds 0] (fun d : ℝ => d ^ 5)

end T09
-- END STATEMENTS
open Filter Asymptotics

noncomputable def T09.poly8 (c0 c1 c2 c3 c4 c5 c6 c7 c8 d : ℝ) : ℝ :=
  c0 + c1 * d + c2 * d ^ 2 + c3 * d ^ 3 + c4 * d ^ 4 + c5 * d ^ 5 + c6 * d ^ 6 + c7 * d ^ 7
    + c8 * d ^ 8

noncomputable def T09.Rp (p z d : ℝ) : ℝ :=
  T09.poly8
    ((1/1 : ℝ) * p^0 * z^0 + (-64/9 : ℝ) * p^0 * z^1 + (2/9 : ℝ) * p^2 * z^0 + (16/81 : ℝ) * p^2 * z^1 + (16/135 : ℝ) * p^4 * z^0)
    ((8/1 : ℝ) * p^0 * z^0 + (-352/9 : ℝ) * p^0 * z^1 + (128/81 : ℝ) * p^0 * z^2 + (11/9 : ℝ) * p^2 * z^0 + (16/27 : ℝ) * p^2 * z^1 + (224/405 : ℝ) * p^4 * z^0 + (37/7290 : ℝ) * p^6 * z^0)
    ((28/1 : ℝ) * p^0 * z^0 + (-752/9 : ℝ) * p^0 * z^1 + (256/81 : ℝ) * p^0 * z^2 + (22/9 : ℝ) * p^2 * z^0 + (16/27 : ℝ) * p^2 * z^1 + (128/135 : ℝ) * p^4 * z^0 + (296/3645 : ℝ) * p^4 * z^1 + (37/3645 : ℝ) * p^6 * z^0)
    ((56/1 : ℝ) * p^0 * z^0 + (-688/9 : ℝ) * p^0 * z^1 + (128/81 : ℝ) * p^0 * z^2 + (13/9 : ℝ) * p^2 * z^0 + (16/81 : ℝ) * p^2 * z^1 + (208/405 : ℝ) * p^4 * z^0 + (296/3645 : ℝ) * p^4 * z^1 + (37/7290 : ℝ) * p^6 * z^0 + (1369/1312200 : ℝ) * p^8 * z^0)
    ((70/1 : ℝ) * p^0 * z^0 + (-20/9 : ℝ) * p^2 * z^0 + (-224/405 : ℝ) * p^4 * z^0)
    ((56/1 : ℝ) * p^0 * z^0 + (64/1 : ℝ) * p^0 * z^1 + (-43/9 : ℝ) * p^2 * z^0 + (-448/405 : ℝ) * p^4 * z^0)
    ((28/1 : ℝ) * p^0 * z^0 + (176/3 : ℝ) * p^0 * z^1 + (-34/9 : ℝ) * p^2 * z^0 + (-64/81 : ℝ) * p^4 * z^0)
    ((8/1 : ℝ) * p^0 * z^0 + (208/9 : ℝ) * p^0 * z^1 + (-13/9 : ℝ) * p^2 * z^0 + (-112/405 : ℝ) * p^4 * z^0)
    ((1/1 : ℝ) * p^0 * z^0 + (32/9 : ℝ) * p^0 * z^1 + (-2/9 : ℝ) * p^2 * z^0 + (-16/405 : ℝ) * p^4 * z^0)
    d / ((1 - d) * (1 + d) ^ 8)

lemma T09.alg (z d : ℝ) (h1 : 1 - d ≠ 0) (h2 : 1 + d ≠ 0) :
    (1 + T09.logSeries z (d / (1 + d)) + (T09.logSeries z (d / (1 + d))) ^ 2 / 2) / (1 - d)
      - T09.dSeries z d = d ^ 5 * T09.Rp Real.pi z d := by
  unfold T09.logSeries T09.dSeries T09.Rp T09.poly8
  field_simp
  ring

lemma T09.poly8_cont (c0 c1 c2 c3 c4 c5 c6 c7 c8 : ℝ) :
    Continuous (fun d : ℝ => T09.poly8 c0 c1 c2 c3 c4 c5 c6 c7 c8 d) := by
  unfold T09.poly8
  fun_prop

lemma T09.Rp_cont (p z : ℝ) : ContinuousAt (fun d : ℝ => T09.Rp p z d) 0 := by
  unfold T09.Rp
  exact ContinuousAt.div (T09.poly8_cont _ _ _ _ _ _ _ _ _).continuousAt
    (by fun_prop : Continuous fun d : ℝ => (1 - d) * (1 + d) ^ 8).continuousAt
    (by show (1 - (0:ℝ)) * (1 + 0) ^ 8 ≠ 0; norm_num)

lemma T09.bdd {w : ℝ → ℝ} (hw : ContinuousAt w 0) (g : ℝ → ℝ) :
    (fun d => g d * w d) =O[nhds 0] g := by
  simpa using (isBigO_refl g (nhds 0)).mul (hw.tendsto.isBigO_one ℝ)

lemma T09.hpow (m k : ℕ) : (fun d : ℝ => d ^ (m + k)) =O[nhds 0] fun d => d ^ m := by
  simpa [pow_add] using T09.bdd (w := fun d : ℝ => d ^ k) (by fun_prop) (fun d => d ^ m)

lemma T09.cinv (a : ℝ) (ha : a ≠ 0) : ContinuousAt (fun d : ℝ => (a + d)⁻¹) 0 :=
  ContinuousAt.inv₀ (continuousAt_const.add continuousAt_id) (by show a + (0:ℝ) ≠ 0; rwa [add_zero])

lemma T09.cinv' : ContinuousAt (fun d : ℝ => (1 - d)⁻¹) 0 :=
  ContinuousAt.inv₀ (continuousAt_const.sub continuousAt_id) (by show (1:ℝ) - 0 ≠ 0; rw [sub_zero]; exact one_ne_zero)

lemma T09.dpow_one (k : ℕ) (hk : 0 < k) : (fun d : ℝ => d ^ k) =O[nhds 0] (1 : ℝ → ℝ) := by
  apply Filter.Tendsto.isBigO_one ℝ (c := (0:ℝ))
  simpa [zero_pow hk.ne'] using ((continuous_pow k).tendsto (0:ℝ))

lemma T09.expstep (x G y : ℝ → ℝ) (hy : y =O[nhds 0] fun d : ℝ => d ^ 2)
    (hx : x =O[nhds 0] fun d : ℝ => d ^ 2) (hG : G =O[nhds 0] fun d : ℝ => d ^ 5)
    (hyx : ∀ d, y d = x d + G d) :
    (fun d => Real.exp (y d) - (1 + x d + x d ^ 2 / 2)) =O[nhds 0] fun d : ℝ => d ^ 5 := by
  have hy0 : Tendsto y (nhds 0) (nhds 0) :=
    hy.trans_tendsto (by simpa using ((continuous_pow 2).tendsto (0:ℝ)))
  have hy1 : ∀ᶠ d in nhds (0:ℝ), |y d| ≤ 1 := by
    filter_upwards [hy0.eventually (Metric.closedBall_mem_nhds (0:ℝ) one_pos)] with d hd
    simpa [Real.norm_eq_abs] using hd
  have h3 : (fun d => Real.exp (y d) - (1 + y d + y d ^ 2 / 2)) =O[nhds 0] fun d => y d ^ 3 := by
    refine IsBigO.of_bound 1 ?_
    filter_upwards [hy1] with d hd
    have := Real.exp_bound hd (n := 3) (by norm_num)
    simp [Finset.sum_range_succ, Nat.factorial] at this
    simp only [Real.norm_eq_abs, abs_pow]
    have h' : |Real.exp (y d) - (1 + y d + y d ^ 2 / 2)| ≤ |y d| ^ 3 * (4 / 18) := by
      convert this using 2 <;> norm_num <;> ring_nf
    nlinarith [pow_nonneg (abs_nonneg (y d)) 3]
  have h6 : (fun d => y d ^ 3) =O[nhds 0] fun d : ℝ => d ^ 6 := by
    have := hy.pow 3
    refine this.trans (isBigO_of_le _ (fun d => ?_))
    simp only [Real.norm_eq_abs, ← pow_mul]
    simp
  have h65 : (fun d : ℝ => d ^ 6) =O[nhds 0] fun d : ℝ => d ^ 5 := T09.hpow 5 1
  have h52 : (fun d : ℝ => d ^ 5) =O[nhds 0] fun d : ℝ => d ^ 2 := T09.hpow 2 3
  have hv : (fun d => y d + x d) =O[nhds 0] (1 : ℝ → ℝ) :=
    (hy.add hx).trans (T09.dpow_one 2 (by norm_num))
  have hGv : (fun d => G d * (1 + (y d + x d) / 2)) =O[nhds 0] fun d : ℝ => d ^ 5 := by
    have h1 : (fun d => 1 + (y d + x d) / 2) =O[nhds 0] (1 : ℝ → ℝ) := by
      refine (isBigO_const_const (1:ℝ) one_ne_zero _).add ?_
      exact (hv.const_mul_left (1/2)).congr_left (fun d => by ring)
    simpa using (hG.mul h1)
  have := (h3.trans (h6.trans h65)).add hGv
  refine this.congr_left (fun d => ?_)
  simp only [hyx d]
  ring

theorem T09.proof_1 : T09.stmt_1 := by
  intro z3 F hpos hG
  have hI : ∀ᶠ d in nhds (0:ℝ), d ∈ Set.Ioo (-1:ℝ) 1 :=
    Ioo_mem_nhds (by norm_num) (by norm_num)
  have hτ5 : (fun d : ℝ => (d / (1 + d)) ^ 5) =O[nhds 0] fun d : ℝ => d ^ 5 := by
    have hc : ContinuousAt (fun d : ℝ => ((1 + d)⁻¹) ^ 5) 0 := by
      have := T09.cinv 1 one_ne_zero
      fun_prop
    have := T09.bdd hc (fun d => d ^ 5)
    refine this.congr_left (fun d => ?_)
    simp [div_eq_mul_inv, mul_pow]
  have hG5 := hG.trans hτ5
  have hx : (fun d : ℝ => T09.logSeries z3 (d / (1 + d))) =O[nhds 0] fun d : ℝ => d ^ 2 := by
    have hw : ContinuousAt (fun d : ℝ => -(Real.pi ^ 2 / 9) * ((1 + d)⁻¹) ^ 2
        - (16 * z3 / 9) * d * ((1 + d)⁻¹) ^ 3 - (37 * Real.pi ^ 4 / 810) * d ^ 2 * ((1 + d)⁻¹) ^ 4) 0 := by
      have := T09.cinv 1 one_ne_zero
      fun_prop
    have := T09.bdd hw (fun d => d ^ 2)
    refine this.congr_left (fun d => ?_)
    simp only [T09.logSeries, div_eq_mul_inv]
    ring
  have hy : (fun d : ℝ => Real.log ((1 - d) * F d)) =O[nhds 0] fun d : ℝ => d ^ 2 := by
    have := hx.add (hG5.trans (T09.hpow 2 3))
    refine this.congr_left (fun d => ?_)
    simp
  have hE := T09.expstep _ _ _ hy hx hG5 (fun d => by simp)
  have hT1 : (fun d : ℝ => (Real.exp (Real.log ((1 - d) * F d))
      - (1 + T09.logSeries z3 (d / (1 + d)) + (T09.logSeries z3 (d / (1 + d))) ^ 2 / 2))
      * (1 - d)⁻¹) =O[nhds 0] fun d : ℝ => d ^ 5 := by
    simpa using hE.mul (T09.cinv'.tendsto.isBigO_one ℝ)
  have hRq : ContinuousAt (fun d : ℝ => T09.Rp Real.pi z3 d) 0 := T09.Rp_cont _ _
  have hT2 : (fun d : ℝ => d ^ 5 * T09.Rp Real.pi z3 d) =O[nhds 0] fun d : ℝ => d ^ 5 :=
    T09.bdd hRq _
  refine (hT1.add hT2).congr' ?_ EventuallyEq.rfl
  filter_upwards [hI, hpos] with d hd hF
  have h1 : 1 - d ≠ 0 := by linarith [hd.2]
  have h2 : 1 + d ≠ 0 := by linarith [hd.1]
  have hexp : Real.exp (Real.log ((1 - d) * F d)) = (1 - d) * F d :=
    Real.exp_log (mul_pos (by linarith [hd.2]) hF)
  have hA := T09.alg z3 d h1 h2
  show _ = F d - T09.dSeries z3 d
  rw [hexp, ← hA]
  field_simp
  ring
