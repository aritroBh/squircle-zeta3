import Mathlib
/-!
# D00 -- the real dilogarithm library

Target id: D00.
Paper labels served: `lem:dilog` (its proof), `thm:M1`, `prop:M1short`, `thm:M2` eq. `(T)`,
`thm:Malpha`. This file supplies the real dilogarithm and the classical identities that the
paper's proofs use; it does not itself state any of the paper's theorems.

Quoted verbatim (main.tex 2325, Appendix "Proof of Theorem thm:M1"):

  We use the dilogarithm on $(-\infty,1]$, its real part for arguments above $1$, and the
  classical identities~\cite{Lewin81,Zagier}: reflection
  $\Li_2(z)+\Li_2(1-z)=\frac{\pi^2}{6}-\log z\log(1-z)$ ($0<z<1$); inversion
  $\Li_2(z)+\Li_2(1/z)=-\frac{\pi^2}{6}-\frac12\log^2(-z)$ ($z<0$) and
  $\operatorname{Re}\Li_2(z)+\Li_2(1/z)=\frac{\pi^2}3-\frac12\log^2z$ ($z>1$); Landen
  $\Li_2(z)+\Li_2(\frac{z}{z-1})=-\frac12\log^2(1-z)$ ($z<1$); duplication
  $\Li_2(w^2)=2\Li_2(w)+2\Li_2(-w)$; and Abel's five-term relation
  \begin{multline*}
  \Li_2\Bigl(\frac{x}{1-y}\Bigr)+\Li_2\Bigl(\frac{y}{1-x}\Bigr)-\Li_2\Bigl(\frac{xy}{(1-x)(1-y)}\Bigr)\\
  =\Li_2(x)+\Li_2(y)+\log(1-x)\log(1-y)
  \end{multline*}
  on $\{x<1,\ y<1,\ x+y<1\}$.

Quoted verbatim (main.tex 2346, proof of `lem:dilog`):

  For the first line, $\chi_2(x)=\frac12(\Li_2(x)-\Li_2(-x))$ satisfies Landen's identity
  $\chi_2(\frac{1-x}{1+x})+\chi_2(x)=\frac{\pi^2}8+\frac12\log x\log\frac{1+x}{1-x}$, and $x=b$
  is its fixed point, with $\frac{1+b}{1-b}=a$. The second line is duplication at $w=b$. The
  third is reflection at $z=1-b$. The fourth is Landen at $z=-\sqrt2$. The fifth is the
  five-term relation at $(x,y)=(-\sqrt2,1-\sqrt2)$ ...

Quoted verbatim (main.tex 877): We use $\Li_n(z)=\sum_{k\ge1}z^k/k^n$, ...
Quoted verbatim (main.tex 904, eq. (T)):
  T=\frac{\pi^2+6\lambda^2}{24}=\frac{\zeta(2)-\Li_2(\frac12)}{2}

## Definition

`Li2 x = -∫_0^x log(1-t)/t dt` (interval integral, `Real.log`).
* For `x ≤ 1` this is the genuine real dilogarithm (the integrand extends continuously to
  `t = 0` with value `-1`; at `x = 1` it has an integrable log singularity at `t = 1`).
  For `x < 0` the oriented interval integral gives the standard analytic continuation.
* For `x > 1`, because `Real.log` is `log |·|`, the integrand is `log|1-t|/t`, still integrable,
  and `Li2 x` equals `Re Li₂(x)` -- exactly the paper's convention "its real part for arguments
  above 1". So `Li2` is meaningful on all of `ℝ`; only `stmt_12` uses it above 1.
* `stmt_1` asserts the integrand is interval integrable on `[0,x]` for every real `x`, so the
  definition never falls back on Lean's "non-integrable ⇒ 0" convention.

`chi2 x = (Li2 x - Li2 (-x))/2` (Legendre chi, the paper's `χ_2`).

## Statements

* stmt_1: integrability of `t ↦ log(1-t)/t` on `[0,x]` for every real `x`.
* stmt_2: `Li2 0 = 0`.
* stmt_3: `Li2` is continuous on `(-∞,1]`.
* stmt_4: `HasDerivAt Li2 (-log(1-x)/x) x` for `x < 1`, `x ≠ 0`.
* stmt_5: `HasDerivAt Li2 1 0`.
* stmt_6: power series `Σ_{n≥0} x^(n+1)/(n+1)^2 = Li2 x` (`HasSum`) for `|x| ≤ 1`.
* stmt_7: `Li2 1 = π²/6`.
* stmt_8: `Li2 (-1) = -π²/12`.
* stmt_9: reflection, `0 < x < 1`.
* stmt_10: Landen, `x < 1`.
* stmt_11: inversion, `x < 0`.
* stmt_12: inversion above 1 (real part), `x > 1`: `Li2 x + Li2 (1/x) = π²/3 - ½ log² x`.
* stmt_13: duplication `Li2 x + Li2 (-x) = ½ Li2 (x²)` for `|x| ≤ 1` (the paper's
  `Li2(w²) = 2Li2(w) + 2Li2(-w)`; the paper gives no domain, we use `|w| ≤ 1`, where all three
  arguments are in `[-1,1]`).
* stmt_14: `Li2 (1/2) = π²/12 - ½ log² 2`.
* stmt_15: the Legendre-chi ("Landen") identity, `0 < x < 1`.
* stmt_16: Abel's five-term relation on `{x < 1, y < 1, x + y < 1}`.

NOT COVERED: the complex dilogarithm; `Li_n` for `n ≥ 3`; duplication for `|w| > 1`
(not needed: the paper applies it at `w = b = √2-1`); the specific evaluations of `lem:dilog`
(those belong to the files of `thm:M1`, which import this one).
-/
-- BEGIN STATEMENTS
namespace D00

open Real

/-- The real dilogarithm `Li₂(x) = -∫_0^x log(1-t)/t dt`; genuine on `(-∞,1]`, equal to
`Re Li₂(x)` for `x > 1` (since `Real.log` is `log |·|`). -/
noncomputable def Li2 (x : ℝ) : ℝ := -∫ t in (0 : ℝ)..x, Real.log (1 - t) / t

/-- Legendre's chi function `χ₂(x) = ½(Li₂(x) - Li₂(-x))` (main.tex 2346). -/
noncomputable def chi2 (x : ℝ) : ℝ := (Li2 x - Li2 (-x)) / 2

/-- The defining integrand is integrable on every `[0,x]`, so `Li2` never uses the
junk value `0` of a non-integrable interval integral. -/
def stmt_1 : Prop :=
  ∀ x : ℝ, IntervalIntegrable (fun t : ℝ => Real.log (1 - t) / t) MeasureTheory.volume 0 x

/-- `Li₂(0) = 0`. -/
def stmt_2 : Prop := Li2 0 = 0

/-- "the dilogarithm on `(-∞,1]`": continuity there. -/
def stmt_3 : Prop := ContinuousOn Li2 (Set.Iic 1)

/-- `Li₂'(x) = -log(1-x)/x` for `x < 1`, `x ≠ 0`. -/
def stmt_4 : Prop :=
  ∀ x : ℝ, x < 1 → x ≠ 0 → HasDerivAt Li2 (-Real.log (1 - x) / x) x

/-- `Li₂'(0) = 1`. -/
def stmt_5 : Prop := HasDerivAt Li2 1 0

/-- `Li₂(x) = Σ_{k≥1} x^k/k²` for `|x| ≤ 1` (main.tex 877, `n = 2`). -/
def stmt_6 : Prop :=
  ∀ x : ℝ, |x| ≤ 1 →
    HasSum (fun n : ℕ => x ^ (n + 1) / ((n : ℝ) + 1) ^ 2) (Li2 x)

/-- `Li₂(1) = π²/6`. -/
def stmt_7 : Prop := Li2 1 = π ^ 2 / 6

/-- `Li₂(-1) = -π²/12`. -/
def stmt_8 : Prop := Li2 (-1) = -(π ^ 2 / 12)

/-- Reflection: `Li₂(z)+Li₂(1-z) = π²/6 - log z log(1-z)` for `0<z<1`. -/
def stmt_9 : Prop :=
  ∀ x : ℝ, 0 < x → x < 1 →
    Li2 x + Li2 (1 - x) = π ^ 2 / 6 - Real.log x * Real.log (1 - x)

/-- Landen: `Li₂(z)+Li₂(z/(z-1)) = -½ log²(1-z)` for `z<1`. -/
def stmt_10 : Prop :=
  ∀ x : ℝ, x < 1 →
    Li2 x + Li2 (x / (x - 1)) = -(1 / 2) * Real.log (1 - x) ^ 2

/-- Inversion: `Li₂(z)+Li₂(1/z) = -π²/6 - ½ log²(-z)` for `z<0`. -/
def stmt_11 : Prop :=
  ∀ x : ℝ, x < 0 →
    Li2 x + Li2 (1 / x) = -(π ^ 2 / 6) - (1 / 2) * Real.log (-x) ^ 2

/-- Inversion above 1: `Re Li₂(z)+Li₂(1/z) = π²/3 - ½ log² z` for `z>1`
(here `Li2 z` is `Re Li₂(z)`). -/
def stmt_12 : Prop :=
  ∀ x : ℝ, 1 < x →
    Li2 x + Li2 (1 / x) = π ^ 2 / 3 - (1 / 2) * Real.log x ^ 2

/-- Duplication: `Li₂(w²) = 2Li₂(w) + 2Li₂(-w)`, stated for `|w| ≤ 1`. -/
def stmt_13 : Prop :=
  ∀ x : ℝ, |x| ≤ 1 → Li2 x + Li2 (-x) = (1 / 2) * Li2 (x ^ 2)

/-- `Li₂(½) = π²/12 - ½ log² 2` (used in eq. (T)). -/
def stmt_14 : Prop := Li2 (1 / 2) = π ^ 2 / 12 - (1 / 2) * Real.log 2 ^ 2

/-- Landen's identity for `χ₂`:
`χ₂((1-x)/(1+x)) + χ₂(x) = π²/8 + ½ log x log((1+x)/(1-x))` for `0<x<1`. -/
def stmt_15 : Prop :=
  ∀ x : ℝ, 0 < x → x < 1 →
    chi2 ((1 - x) / (1 + x)) + chi2 x
      = π ^ 2 / 8 + (1 / 2) * Real.log x * Real.log ((1 + x) / (1 - x))

/-- Abel's five-term relation on `{x<1, y<1, x+y<1}`. -/
def stmt_16 : Prop :=
  ∀ x y : ℝ, x < 1 → y < 1 → x + y < 1 →
    Li2 (x / (1 - y)) + Li2 (y / (1 - x)) - Li2 (x * y / ((1 - x) * (1 - y)))
      = Li2 x + Li2 y + Real.log (1 - x) * Real.log (1 - y)

end D00
-- END STATEMENTS

namespace D00

open Real Filter Topology Set

/-- The continuous extension of the integrand `log(1-t)/t` (value `-1` at `0`). -/
noncomputable def G (t : ℝ) : ℝ := if t = 0 then -1 else Real.log (1 - t) / t

lemma G_of_ne {t : ℝ} (h : t ≠ 0) : G t = Real.log (1 - t) / t := by simp [G, h]

lemma G_zero : G 0 = -1 := by simp [G]

lemma continuousAt_G_zero : ContinuousAt G 0 := by
  rw [← continuousWithinAt_compl_self, ContinuousWithinAt, G_zero]
  have hd : HasDerivAt (fun t : ℝ => Real.log (1 - t)) (-1) 0 := by
    have := ((hasDerivAt_id (0:ℝ)).const_sub 1).log (by norm_num)
    simpa using this
  rw [hasDerivAt_iff_tendsto_slope] at hd
  apply hd.congr'
  filter_upwards [self_mem_nhdsWithin] with t ht
  rw [G_of_ne ht, slope_def_field]
  simp

lemma continuousAt_G {t : ℝ} (ht : t ≠ 1) : ContinuousAt G t := by
  rcases eq_or_ne t 0 with rfl | h0
  · exact continuousAt_G_zero
  · have hc : ContinuousAt (fun s : ℝ => Real.log (1 - s) / s) t :=
      ((continuousAt_const.sub continuousAt_id).log (sub_ne_zero.mpr (Ne.symm ht))).div
        continuousAt_id h0
    apply hc.congr
    filter_upwards [isOpen_ne.mem_nhds h0] with s hs
    exact (G_of_ne hs).symm

lemma intervalIntegrable_G_le (c : ℝ) (hc : c ≤ 1 / 2) :
    IntervalIntegrable G MeasureTheory.volume c (1 / 2) := by
  apply ContinuousOn.intervalIntegrable
  intro t ht
  rw [uIcc_of_le hc] at ht
  exact (continuousAt_G (by linarith [ht.2])).continuousWithinAt

lemma intervalIntegrable_g_le (c : ℝ) (hc : c ≤ 1 / 2) :
    IntervalIntegrable (fun t : ℝ => Real.log (1 - t) / t) MeasureTheory.volume c (1 / 2) := by
  apply (intervalIntegrable_G_le c hc).congr_ae
  apply MeasureTheory.ae_restrict_of_ae
  filter_upwards [(Set.countable_singleton (0:ℝ)).ae_notMem MeasureTheory.volume] with t ht
  exact G_of_ne ht

lemma intervalIntegrable_g_ge (c : ℝ) (hc : 1 / 2 ≤ c) :
    IntervalIntegrable (fun t : ℝ => Real.log (1 - t) / t) MeasureTheory.volume (1 / 2) c := by
  have h1 : IntervalIntegrable (fun t : ℝ => Real.log (1 - t)) MeasureTheory.volume (1 / 2) c := by
    have := (intervalIntegral.intervalIntegrable_log' (a := 1 - 1 / 2) (b := 1 - c)).comp_sub_left 1
    simpa using this
  have h2 := h1.mul_continuousOn (g := fun t : ℝ => t⁻¹) (by
    intro t ht
    rw [uIcc_of_le hc] at ht
    exact (continuousAt_inv₀ (by linarith [ht.1] : t ≠ 0)).continuousWithinAt)
  simpa [div_eq_mul_inv] using h2

lemma intervalIntegrable_g (a b : ℝ) :
    IntervalIntegrable (fun t : ℝ => Real.log (1 - t) / t) MeasureTheory.volume a b := by
  have key : ∀ c, IntervalIntegrable (fun t : ℝ => Real.log (1 - t) / t)
      MeasureTheory.volume c (1 / 2) := by
    intro c
    rcases le_total c (1 / 2) with h | h
    · exact intervalIntegrable_g_le c h
    · exact (intervalIntegrable_g_ge c h).symm
  exact (key a).trans (key b).symm

lemma Li2_zero' : Li2 0 = 0 := by simp [Li2]

lemma Li2_eq_G (x : ℝ) : Li2 x = -∫ t in (0 : ℝ)..x, G t := by
  rw [Li2, intervalIntegral.integral_congr_uIoo]
  intro t ht
  have h0 : t ≠ 0 := by
    rintro rfl
    simp only [uIoo, mem_Ioo] at ht
    rcases le_total 0 x with h | h
    · simp [h] at ht
    · simp [h] at ht
  exact (G_of_ne h0).symm

lemma intervalIntegrable_G (x : ℝ) : IntervalIntegrable G MeasureTheory.volume 0 x := by
  apply (intervalIntegrable_g 0 x).congr_ae
  apply MeasureTheory.ae_restrict_of_ae
  filter_upwards [(Set.countable_singleton (0:ℝ)).ae_notMem MeasureTheory.volume] with t ht
  exact (G_of_ne ht).symm

lemma hasDerivAt_Li2_G {x : ℝ} (hx : x ≠ 1) : HasDerivAt Li2 (-G x) x := by
  have hm : StronglyMeasurableAtFilter G (𝓝 x) MeasureTheory.volume := by
    have hco : ContinuousOn G {1}ᶜ := fun t ht => (continuousAt_G ht).continuousWithinAt
    exact hco.stronglyMeasurableAtFilter isOpen_compl_singleton x hx
  have := (intervalIntegral.integral_hasDerivAt_right (intervalIntegrable_G x) hm
    (continuousAt_G hx)).neg
  convert this using 1
  funext y
  exact Li2_eq_G y

lemma hasDerivAt_Li2_ne {x : ℝ} (hx : x ≠ 1) (h0 : x ≠ 0) :
    HasDerivAt Li2 (-Real.log (1 - x) / x) x := by
  have := hasDerivAt_Li2_G hx
  rwa [G_of_ne h0, ← neg_div] at this

lemma continuous_Li2 : Continuous Li2 := by
  unfold Li2
  exact (intervalIntegral.continuous_primitive (fun a b => intervalIntegrable_g a b) 0).neg

/-- Mean value consequence: zero derivative inside, continuity on the closed segment. -/
lemma eq_of_deriv_zero {f : ℝ → ℝ} {a b : ℝ} (hab : a ≤ b) (hc : ContinuousOn f (Icc a b))
    (hd : ∀ x ∈ Ioo a b, HasDerivAt f 0 x) : f a = f b := by
  rcases hab.eq_or_lt with rfl | hlt
  · rfl
  obtain ⟨c, _, hc'⟩ := exists_hasDerivAt_eq_slope f (fun _ => (0:ℝ)) hlt hc hd
  have hne : b - a ≠ 0 := sub_ne_zero.mpr hlt.ne'
  field_simp at hc'
  linarith

/-- Symmetric version: for any two points. -/
lemma eq_of_deriv_zero' {f : ℝ → ℝ} {a b : ℝ} (hc : ContinuousOn f (uIcc a b))
    (hd : ∀ x ∈ uIoo a b, HasDerivAt f 0 x) : f a = f b := by
  rcases le_total a b with h | h
  · apply eq_of_deriv_zero h
    · simpa [uIcc_of_le h] using hc
    · intro x hx; apply hd; simpa [uIoo, h] using hx
  · symm
    apply eq_of_deriv_zero h
    · simpa [uIcc_of_ge h] using hc
    · intro x hx; apply hd; simpa [uIoo, h] using hx

/-- The power series `Σ x^(n+1)/(n+1)²`. -/
noncomputable def S (x : ℝ) : ℝ := ∑' n : ℕ, x ^ (n + 1) / ((n : ℝ) + 1) ^ 2

lemma summable_inv_sq : Summable (fun n : ℕ => 1 / ((n : ℝ) + 1) ^ 2) := by
  have := (summable_nat_add_iff 1).mpr (Real.summable_one_div_nat_pow.mpr (by norm_num : 1 < 2))
  simpa using this

lemma term_bound {x : ℝ} (hx : |x| ≤ 1) (n : ℕ) :
    ‖x ^ (n + 1) / ((n : ℝ) + 1) ^ 2‖ ≤ 1 / ((n : ℝ) + 1) ^ 2 := by
  rw [Real.norm_eq_abs, abs_div, abs_pow, abs_of_pos (by positivity : (0:ℝ) < ((n:ℝ) + 1) ^ 2)]
  apply div_le_div_of_nonneg_right _ (by positivity)
  exact pow_le_one₀ (abs_nonneg x) hx

lemma summable_S {x : ℝ} (hx : |x| ≤ 1) :
    Summable (fun n : ℕ => x ^ (n + 1) / ((n : ℝ) + 1) ^ 2) :=
  Summable.of_norm_bounded summable_inv_sq (term_bound hx)

lemma continuousOn_S : ContinuousOn S (Icc (-1) 1) := by
  apply continuousOn_tsum (fun n => (continuous_pow (n + 1)).continuousOn.div_const _)
    summable_inv_sq
  intro n x hx
  exact term_bound (abs_le.mpr hx) n

lemma S_zero : S 0 = 0 := by simp [S]

lemma hasDerivAt_S {y : ℝ} (hy : |y| < 1) (h0 : y ≠ 0) : HasDerivAt S (-G y) y := by
  set r := (|y| + 1) / 2 with hr
  have hr1 : r < 1 := by rw [hr]; linarith
  have hyr : |y| < r := by rw [hr]; linarith
  have hr0 : 0 ≤ r := by rw [hr]; positivity
  have ht : y ∈ Ioo (-r) r := ⟨by linarith [neg_abs_le y], by linarith [le_abs_self y]⟩
  have hd := hasDerivAt_tsum_of_isPreconnected (t := Ioo (-r) r) (y₀ := 0)
    (g := fun (n : ℕ) (z : ℝ) => z ^ (n + 1) / ((n : ℝ) + 1) ^ 2)
    (g' := fun (n : ℕ) (z : ℝ) => z ^ n / ((n : ℝ) + 1)) (u := fun n : ℕ => r ^ n)
    (summable_geometric_of_lt_one hr0 hr1) isOpen_Ioo isPreconnected_Ioo
    (by
      intro n z _
      have := (hasDerivAt_pow (n + 1) z).div_const (((n : ℝ) + 1) ^ 2)
      refine this.congr_deriv ?_
      push_cast
      field_simp)
    (by
      intro n z hz
      rw [Real.norm_eq_abs, abs_div, abs_pow, abs_of_pos (by positivity : (0:ℝ) < (n:ℝ) + 1)]
      have hz' : |z| ≤ r := abs_le.mpr ⟨hz.1.le, hz.2.le⟩
      calc |z| ^ n / ((n:ℝ) + 1) ≤ |z| ^ n := by
            apply div_le_self (by positivity)
            have : (0:ℝ) ≤ n := Nat.cast_nonneg n
            linarith
        _ ≤ r ^ n := pow_le_pow_left₀ (abs_nonneg z) hz' n)
    ⟨by linarith [abs_nonneg y, hr], by linarith [abs_nonneg y, hr]⟩
    (by simp)
    ht
  have hsum : HasSum (fun n : ℕ => y ^ n / ((n : ℝ) + 1)) (-G y) := by
    have := (Real.hasSum_pow_div_log_of_abs_lt_one hy).div_const y
    rw [G_of_ne h0, ← neg_div]
    convert this using 1
    funext n
    field_simp
    ring
  rw [hsum.tsum_eq] at hd
  exact hd

lemma Li2_eq_S {x : ℝ} (hx : |x| ≤ 1) : Li2 x = S x := by
  have hx' := abs_le.mp hx
  have key : (fun z => Li2 z - S z) 0 = (fun z => Li2 z - S z) x := by
    apply eq_of_deriv_zero' (f := fun z => Li2 z - S z)
    · apply continuous_Li2.continuousOn.sub
      apply continuousOn_S.mono
      intro z hz
      rw [mem_uIcc] at hz
      constructor <;> rcases hz with h | h <;> linarith [h.1, h.2]
    · intro z hz
      simp only [uIoo, mem_Ioo] at hz
      have hz0 : z ≠ 0 := by
        rintro rfl
        rcases le_total 0 x with h | h <;> simp [h] at hz
      have hz1 : |z| < 1 := by
        rw [abs_lt]
        rcases le_total 0 x with h | h <;> simp [h] at hz <;> constructor <;> linarith
      have := (hasDerivAt_Li2_G (by linarith [le_abs_self z] : z ≠ 1)).sub (hasDerivAt_S hz1 hz0)
      rw [sub_self] at this
      exact this
  simp [S_zero, Li2_zero'] at key
  linarith

lemma Li2_comp {f : ℝ → ℝ} {f' x : ℝ} (hf : HasDerivAt f f' x) (h1 : f x ≠ 1) (h0 : f x ≠ 0) :
    HasDerivAt (fun y => Li2 (f y)) (-Real.log (1 - f x) / f x * f') x :=
  (hasDerivAt_Li2_ne h1 h0).comp x hf

lemma mem_uIoo_zero {x z : ℝ} (hz : z ∈ uIoo 0 x) : z ≠ 0 ∧ |z| < |x| := by
  simp only [uIoo, mem_Ioo] at hz
  rcases le_total 0 x with h | h
  · rw [min_eq_left h, max_eq_right h] at hz
    exact ⟨hz.1.ne', by rw [abs_of_pos hz.1, abs_of_nonneg h]; exact hz.2⟩
  · rw [min_eq_right h, max_eq_left h] at hz
    exact ⟨hz.2.ne, by rw [abs_of_neg hz.2, abs_of_nonpos h]; linarith [hz.1]⟩

lemma duplication' {x : ℝ} (hx : |x| ≤ 1) : Li2 x + Li2 (-x) = (1 / 2) * Li2 (x ^ 2) := by
  have key := eq_of_deriv_zero' (f := fun z => Li2 z + Li2 (-z) - (1/2) * Li2 (z^2))
    (a := 0) (b := x)
    ((continuous_Li2.add (continuous_Li2.comp continuous_neg)).sub
      (continuous_const.mul (continuous_Li2.comp (continuous_pow 2)))).continuousOn
    (by
      intro z hz
      obtain ⟨hz0, hzx⟩ := mem_uIoo_zero hz
      have hz1 : |z| < 1 := lt_of_lt_of_le hzx hx
      have hzl := abs_lt.mp hz1
      have hsq0 : z ^ 2 ≠ 0 := pow_ne_zero 2 hz0
      have hsq1 : z ^ 2 < 1 := by nlinarith
      have d1 := Li2_comp (hasDerivAt_id z) (by simp; linarith) (by simpa using hz0)
      have d2 := Li2_comp (hasDerivAt_neg z) (by simp; linarith) (by simpa using hz0)
      have d3 := Li2_comp (hasDerivAt_pow 2 z) hsq1.ne hsq0
      have := (d1.add d2).sub (d3.const_mul (1/2))
      refine this.congr_deriv ?_
      have h1 : (0:ℝ) < 1 - z := by linarith
      have h2 : (0:ℝ) < 1 + z := by linarith
      have e : Real.log (1 - z ^ 2) = Real.log (1 - z) + Real.log (1 + z) := by
        rw [show 1 - z ^ 2 = (1 - z) * (1 + z) by ring, Real.log_mul h1.ne' h2.ne']
      simp only [id, e, sub_neg_eq_add]
      field_simp
      ring)
  simp [Li2_zero'] at key
  linarith

lemma Li2_one' : Li2 1 = π ^ 2 / 6 := by
  have h1 := Li2_eq_S (x := 1) (by norm_num)
  have h2 : HasSum (fun n : ℕ => (1 : ℝ) / ((n : ℝ) + 1) ^ 2) (Real.pi ^ 2 / 6) := by
    have := (hasSum_nat_add_iff' 1).mpr hasSum_zeta_two
    simpa using this
  rw [h1, S]
  simp only [one_pow]
  exact h2.tsum_eq

lemma Li2_neg_one' : Li2 (-1) = -(π ^ 2 / 12) := by
  have := duplication' (x := 1) (by norm_num)
  rw [one_pow, Li2_one'] at this
  linarith

lemma inversion_neg' {x : ℝ} (hx : x < 0) :
    Li2 x + Li2 (1 / x) = -(π ^ 2 / 6) - (1 / 2) * Real.log (-x) ^ 2 := by
  have key := eq_of_deriv_zero'
    (f := fun z => Li2 z + Li2 z⁻¹ + (1/2) * Real.log (-z) ^ 2) (a := -1) (b := x)
    (by
      intro z hz
      have hz' : z < 0 := by
        rw [mem_uIcc] at hz
        rcases hz with h | h <;> linarith [h.1, h.2]
      have hz0 : z ≠ 0 := hz'.ne
      have hn : -z ≠ 0 := by linarith
      apply ContinuousAt.continuousWithinAt
      have c1 : ContinuousAt (fun z : ℝ => Li2 z⁻¹) z :=
        continuous_Li2.continuousAt.comp (continuousAt_inv₀ hz0)
      have c2 : ContinuousAt (fun z : ℝ => Real.log (-z)) z :=
        continuousAt_id.neg.log hn
      exact (continuous_Li2.continuousAt.add c1).add (continuousAt_const.mul (c2.pow 2)))
    (by
      intro z hz
      have hz' : z < 0 := by
        simp only [uIoo, mem_Ioo] at hz
        rcases le_total (-1) x with h | h
        · rw [min_eq_left h, max_eq_right h] at hz; linarith [hz.2]
        · rw [min_eq_right h, max_eq_left h] at hz; linarith [hz.2]
      have hz0 : z ≠ 0 := hz'.ne
      have hn : -z ≠ 0 := by linarith
      have d1 := Li2_comp (hasDerivAt_id z) (by simp; linarith) (by simpa using hz0)
      have d2 := Li2_comp (hasDerivAt_inv hz0) (by
        intro h; have : z = 1 := by simpa using h
        linarith) (inv_ne_zero hz0)
      have d3 := (((hasDerivAt_neg z).log hn).pow 2).const_mul (1/2 : ℝ)
      have := (d1.add d2).add d3
      refine this.congr_deriv ?_
      have h1 : (0:ℝ) < 1 - z := by linarith
      have e : Real.log (1 - z⁻¹) = Real.log (1 - z) - Real.log (-z) := by
        rw [show 1 - z⁻¹ = (1 - z) / (-z) by field_simp; ring, Real.log_div h1.ne' hn]
      simp only [id, e]
      field_simp
      ring)
  simp only [inv_neg, inv_one, neg_neg, Real.log_one] at key
  rw [Li2_neg_one'] at key
  rw [one_div, Real.log_neg_eq_log]
  norm_num at key
  linarith

lemma inversion_gt_one' {x : ℝ} (hx : 1 < x) :
    Li2 x + Li2 (1 / x) = π ^ 2 / 3 - (1 / 2) * Real.log x ^ 2 := by
  have key := eq_of_deriv_zero'
    (f := fun z => Li2 z + Li2 z⁻¹ + (1/2) * Real.log z ^ 2) (a := 1) (b := x)
    (by
      intro z hz
      have hz' : 1 ≤ z := by
        rw [mem_uIcc] at hz
        rcases hz with h | h <;> linarith [h.1, h.2]
      have hz0 : z ≠ 0 := by linarith
      apply ContinuousAt.continuousWithinAt
      have c1 : ContinuousAt (fun z : ℝ => Li2 z⁻¹) z :=
        continuous_Li2.continuousAt.comp (continuousAt_inv₀ hz0)
      have c2 : ContinuousAt (fun z : ℝ => Real.log z) z := Real.continuousAt_log hz0
      exact (continuous_Li2.continuousAt.add c1).add (continuousAt_const.mul (c2.pow 2)))
    (by
      intro z hz
      have hz' : 1 < z := by
        simp only [uIoo, mem_Ioo] at hz
        rw [min_eq_left hx.le, max_eq_right hx.le] at hz; exact hz.1
      have hz0 : z ≠ 0 := by linarith
      have d1 := Li2_comp (hasDerivAt_id z) (by simp; linarith) (by simpa using hz0)
      have d2 := Li2_comp (hasDerivAt_inv hz0) (by
        intro h; have : z = 1 := by simpa using h
        linarith) (inv_ne_zero hz0)
      have d3 := ((Real.hasDerivAt_log hz0).pow 2).const_mul (1/2 : ℝ)
      have := (d1.add d2).add d3
      refine this.congr_deriv ?_
      have h1 : (0:ℝ) < z - 1 := by linarith
      have e : Real.log (1 - z⁻¹) = Real.log (z - 1) - Real.log z := by
        rw [show 1 - z⁻¹ = (z - 1) / z by field_simp, Real.log_div h1.ne' hz0]
      have e2 : Real.log (1 - z) = Real.log (z - 1) := by
        rw [show 1 - z = -(z - 1) by ring, Real.log_neg_eq_log]
      simp only [id, e, e2]
      field_simp
      ring)
  simp only [inv_one, Real.log_one] at key
  rw [Li2_one'] at key
  rw [one_div]
  norm_num at key
  linarith

lemma landen' {x : ℝ} (hx : x < 1) :
    Li2 x + Li2 (x / (x - 1)) = -(1 / 2) * Real.log (1 - x) ^ 2 := by
  have key := eq_of_deriv_zero'
    (f := fun z => Li2 z + Li2 (z / (z - 1)) + (1/2) * Real.log (1 - z) ^ 2) (a := 0) (b := x)
    (by
      intro z hz
      have hz' : z < 1 := by
        rw [mem_uIcc] at hz
        rcases hz with h | h <;> linarith [h.1, h.2]
      have hz1 : z - 1 ≠ 0 := by linarith
      have hz2 : 1 - z ≠ 0 := by linarith
      apply ContinuousAt.continuousWithinAt
      have c1 : ContinuousAt (fun z : ℝ => Li2 (z / (z - 1))) z :=
        continuous_Li2.continuousAt.comp (continuousAt_id.div (continuousAt_id.sub continuousAt_const) hz1)
      have c2 : ContinuousAt (fun z : ℝ => Real.log (1 - z)) z :=
        (continuousAt_const.sub continuousAt_id).log hz2
      exact (continuous_Li2.continuousAt.add c1).add (continuousAt_const.mul (c2.pow 2)))
    (by
      intro z hz
      obtain ⟨hz0, -⟩ := mem_uIoo_zero hz
      have hz' : z < 1 := by
        simp only [uIoo, mem_Ioo] at hz
        rcases le_total 0 x with h | h
        · rw [min_eq_left h, max_eq_right h] at hz; linarith [hz.2]
        · rw [min_eq_right h, max_eq_left h] at hz; linarith [hz.2]
      have hz1 : z - 1 ≠ 0 := by linarith
      have hz2 : 1 - z ≠ 0 := by linarith
      have hq0 : z / (z - 1) ≠ 0 := div_ne_zero hz0 hz1
      have hq1 : z / (z - 1) ≠ 1 := by
        intro h; field_simp at h; linarith
      have d1 := Li2_comp (hasDerivAt_id z) (by simp; linarith) (by simpa using hz0)
      have d2 := Li2_comp ((hasDerivAt_id z).div ((hasDerivAt_id z).sub_const 1) hz1) hq1 hq0
      have d3 := ((((hasDerivAt_id z).const_sub 1).log hz2).pow 2).const_mul (1/2 : ℝ)
      have := (d1.add d2).add d3
      refine this.congr_deriv ?_
      have e : Real.log (1 - z / (z - 1)) = - Real.log (1 - z) := by
        rw [show 1 - z / (z - 1) = (1 - z)⁻¹ by field_simp; ring, Real.log_inv]
      simp only [id, Pi.div_apply, e]
      field_simp
      ring)
  simp [Li2_zero'] at key
  linarith

lemma reflection' {x : ℝ} (h0 : 0 < x) (h1 : x < 1) :
    Li2 x + Li2 (1 - x) = π ^ 2 / 6 - Real.log x * Real.log (1 - x) := by
  have hH : ∀ w : ℝ, w ≠ 0 → w * Real.log w * G w = Real.log w * Real.log (1 - w) := by
    intro w hw; rw [G_of_ne hw]; field_simp
  have key := eq_of_deriv_zero (f := fun z => Li2 z + Li2 (1 - z) + z * Real.log z * G z)
    (a := 0) (b := x) h0.le
    (by
      intro z hz
      have hz1 : z ≠ 1 := by linarith [hz.2]
      apply ContinuousAt.continuousWithinAt
      have c1 : ContinuousAt (fun z : ℝ => Li2 (1 - z)) z :=
        continuous_Li2.continuousAt.comp (continuousAt_const.sub continuousAt_id)
      exact (continuous_Li2.continuousAt.add c1).add
        (continuous_mul_log.continuousAt.mul (continuousAt_G hz1)))
    (by
      intro z hz
      have hz0 : z ≠ 0 := hz.1.ne'
      have hz1' : z < 1 := by linarith [hz.2]
      have hz2 : 1 - z ≠ 0 := by linarith
      have d1 := Li2_comp (hasDerivAt_id z) (by simp; linarith) (by simpa using hz0)
      have d2 := Li2_comp ((hasDerivAt_id z).const_sub 1) (by simpa using hz0) (by simpa using hz2)
      have d3 := (Real.hasDerivAt_log hz0).mul (((hasDerivAt_id z).const_sub 1).log hz2)
      have hd : HasDerivAt (fun z => Li2 z + Li2 (1 - z) + Real.log z * Real.log (1 - z)) 0 z := by
        refine ((d1.add d2).add d3).congr_deriv ?_
        simp only [id, sub_sub_cancel]
        field_simp
        ring
      apply hd.congr_of_eventuallyEq
      filter_upwards [isOpen_ne.mem_nhds hz0] with w hw
      simp only [hH w hw])
  beta_reduce at key
  rw [hH x h0.ne'] at key
  simp [Li2_zero', Li2_one'] at key
  linarith

lemma chi2_landen' {x : ℝ} (h0 : 0 < x) (h1 : x < 1) :
    chi2 ((1 - x) / (1 + x)) + chi2 x
      = π ^ 2 / 8 + (1 / 2) * Real.log x * Real.log ((1 + x) / (1 - x)) := by
  have hK : ∀ w : ℝ, w ≠ 0 → (1/2) * (w * Real.log w * (-G (-w) - G w))
      = (1/2) * (Real.log w * (Real.log (1 + w) - Real.log (1 - w))) := by
    intro w hw
    rw [G_of_ne hw, G_of_ne (neg_ne_zero.mpr hw), sub_neg_eq_add]
    field_simp
  have key := eq_of_deriv_zero
    (f := fun z => (Li2 ((1 - z) / (1 + z)) - Li2 (-((1 - z) / (1 + z)))) / 2
      + (Li2 z - Li2 (-z)) / 2 - (1/2) * (z * Real.log z * (-G (-z) - G z)))
    (a := 0) (b := x) h0.le
    (by
      intro z hz
      have hz1 : z ≠ 1 := by linarith [hz.2]
      have hz1' : -z ≠ 1 := by linarith [hz.1]
      have hz2 : 1 + z ≠ 0 := by linarith [hz.1]
      apply ContinuousAt.continuousWithinAt
      have cq : ContinuousAt (fun z : ℝ => (1 - z) / (1 + z)) z :=
        (continuousAt_const.sub continuousAt_id).div (continuousAt_const.add continuousAt_id) hz2
      have c1 : ContinuousAt (fun z : ℝ => Li2 ((1 - z) / (1 + z))) z :=
        continuous_Li2.continuousAt.comp cq
      have c2 : ContinuousAt (fun z : ℝ => Li2 (-((1 - z) / (1 + z)))) z :=
        continuous_Li2.continuousAt.comp cq.neg
      have c3 : ContinuousAt (fun z : ℝ => Li2 (-z)) z :=
        continuous_Li2.continuousAt.comp continuousAt_id.neg
      have c4 : ContinuousAt (fun z : ℝ => G (-z)) z :=
        (continuousAt_G hz1').comp continuousAt_id.neg
      exact (((c1.sub c2).div_const 2).add ((continuous_Li2.continuousAt.sub c3).div_const 2)).sub
        (continuousAt_const.mul (continuous_mul_log.continuousAt.mul (c4.neg.sub (continuousAt_G hz1)))))
    (by
      intro z hz
      have hz0 : z ≠ 0 := hz.1.ne'
      have hzp : 0 < z := hz.1
      have hz1' : z < 1 := by linarith [hz.2]
      have hA : 0 < 1 - z := by linarith
      have hB : 0 < 1 + z := by linarith
      have hq0 : (1 - z) / (1 + z) ≠ 0 := div_ne_zero hA.ne' hB.ne'
      have hq1 : (1 - z) / (1 + z) ≠ 1 := by
        rw [Ne, div_eq_one_iff_eq hB.ne']; linarith
      have hq : HasDerivAt (fun y : ℝ => (1 - id y) / (1 + id y))
          ((-1 * (1 + id z) - (1 - id z) * 1) / (1 + id z) ^ 2) z :=
        ((hasDerivAt_id z).const_sub 1).div ((hasDerivAt_id z).const_add 1) hB.ne'
      have dA := Li2_comp hq (by simpa using hq1) (by simpa using hq0)
      have dB := Li2_comp hq.neg (by
        simp only [id, Pi.neg_apply]; intro h
        have : (1 - z) / (1 + z) = -1 := by linarith
        rw [div_eq_iff hB.ne'] at this; linarith) (by simp only [id, Pi.neg_apply]; exact neg_ne_zero.mpr hq0)
      have dC := Li2_comp (hasDerivAt_id z) (by simp; linarith) (by simpa using hz0)
      have dD := Li2_comp (hasDerivAt_neg z) (by simp; linarith) (by simpa using hz0)
      have dE := ((Real.hasDerivAt_log hz0).mul ((((hasDerivAt_id z).const_add 1).log hB.ne').sub
        (((hasDerivAt_id z).const_sub 1).log hA.ne'))).const_mul (1/2 : ℝ)
      have e1 : Real.log (1 - (1 - z) / (1 + z)) = Real.log 2 + Real.log z - Real.log (1 + z) := by
        rw [show 1 - (1 - z) / (1 + z) = 2 * z / (1 + z) by (field_simp; ring),
          Real.log_div (by positivity) hB.ne', Real.log_mul (by norm_num) hz0]
      have e2 : Real.log (1 + (1 - z) / (1 + z)) = Real.log 2 - Real.log (1 + z) := by
        rw [show 1 + (1 - z) / (1 + z) = 2 / (1 + z) by (field_simp; ring),
          Real.log_div (by norm_num) hB.ne']
      have hd : HasDerivAt (fun z => (Li2 ((1 - z) / (1 + z)) - Li2 (-((1 - z) / (1 + z)))) / 2
          + (Li2 z - Li2 (-z)) / 2
          - (1/2) * (Real.log z * (Real.log (1 + z) - Real.log (1 - z)))) 0 z := by
        refine ((((dA.sub dB).div_const 2).add ((dC.sub dD).div_const 2)).sub dE).congr_deriv ?_
        simp only [id, Pi.neg_apply, Pi.sub_apply, sub_neg_eq_add, e1, e2]
        field_simp
        ring
      apply hd.congr_of_eventuallyEq
      filter_upwards [isOpen_ne.mem_nhds hz0] with w hw
      simp only [hK w hw])
  beta_reduce at key
  rw [hK x h0.ne'] at key
  simp [Li2_zero', Li2_one', Li2_neg_one'] at key
  have hA : (1 - x) ≠ 0 := by linarith
  have hB : (1 + x) ≠ 0 := by linarith
  unfold chi2
  rw [Real.log_div hB hA]
  ring_nf at key ⊢
  linarith

lemma five_term' {x y : ℝ} (hx : x < 1) (hy : y < 1) (hxy : x + y < 1) :
    Li2 (x / (1 - y)) + Li2 (y / (1 - x)) - Li2 (x * y / ((1 - x) * (1 - y)))
      = Li2 x + Li2 y + Real.log (1 - x) * Real.log (1 - y) := by
  rcases eq_or_ne y 0 with rfl | hy0
  · simp [Li2_zero']
  have hB : 0 < 1 - y := by linarith
  have key := eq_of_deriv_zero'
    (f := fun t => Li2 (t / (1 - y)) + Li2 (y / (1 - t)) - Li2 (t * y / ((1 - t) * (1 - y)))
      - Li2 t - Li2 y - Real.log (1 - t) * Real.log (1 - y)) (a := 0) (b := x)
    (by
      intro t ht
      have ht' : t < 1 := by
        rw [mem_uIcc] at ht
        rcases ht with h | h <;> linarith [h.1, h.2]
      have hA : 1 - t ≠ 0 := by linarith
      apply ContinuousAt.continuousWithinAt
      have c1 : ContinuousAt (fun t : ℝ => Li2 (t / (1 - y))) t :=
        continuous_Li2.continuousAt.comp (continuousAt_id.div_const _)
      have c2 : ContinuousAt (fun t : ℝ => Li2 (y / (1 - t))) t :=
        continuous_Li2.continuousAt.comp
          (continuousAt_const.div (continuousAt_const.sub continuousAt_id) hA)
      have c3 : ContinuousAt (fun t : ℝ => Li2 (t * y / ((1 - t) * (1 - y)))) t :=
        continuous_Li2.continuousAt.comp
          ((continuousAt_id.mul continuousAt_const).div
            ((continuousAt_const.sub continuousAt_id).mul continuousAt_const)
            (mul_ne_zero hA hB.ne'))
      have c4 : ContinuousAt (fun t : ℝ => Real.log (1 - t)) t :=
        (continuousAt_const.sub continuousAt_id).log hA
      exact ((((c1.add c2).sub c3).sub continuous_Li2.continuousAt).sub continuousAt_const).sub
        (c4.mul continuousAt_const))
    (by
      intro t ht
      obtain ⟨ht0, -⟩ := mem_uIoo_zero ht
      have hlt : t < 1 ∧ t + y < 1 := by
        simp only [uIoo, mem_Ioo] at ht
        rcases le_total 0 x with h | h
        · rw [min_eq_left h, max_eq_right h] at ht; constructor <;> linarith [ht.2]
        · rw [min_eq_right h, max_eq_left h] at ht; constructor <;> linarith [ht.2]
      have hA : 0 < 1 - t := by linarith [hlt.1]
      have hL : 0 < 1 - t - y := by linarith [hlt.2]
      have hu0 : t / (1 - y) ≠ 0 := div_ne_zero ht0 hB.ne'
      have hu1 : t / (1 - y) ≠ 1 := by rw [Ne, div_eq_one_iff_eq hB.ne']; linarith
      have hv0 : y / (1 - t) ≠ 0 := div_ne_zero hy0 hA.ne'
      have hv1 : y / (1 - t) ≠ 1 := by rw [Ne, div_eq_one_iff_eq hA.ne']; linarith
      have hw0 : t * y / ((1 - t) * (1 - y)) ≠ 0 :=
        div_ne_zero (mul_ne_zero ht0 hy0) (mul_ne_zero hA.ne' hB.ne')
      have hw1 : t * y / ((1 - t) * (1 - y)) ≠ 1 := by
        rw [Ne, div_eq_one_iff_eq (mul_ne_zero hA.ne' hB.ne')]
        intro h
        have : (1 - t) * (1 - y) - t * y = 1 - t - y := by ring
        linarith
      have du := Li2_comp ((hasDerivAt_id t).div_const (1 - y)) (by simp only [id]; exact hu1)
        (by simp only [id]; exact hu0)
      have dv := Li2_comp ((hasDerivAt_const t y).div ((hasDerivAt_id t).const_sub 1) hA.ne')
        (by simp only [id]; exact hv1) (by simp only [id]; exact hv0)
      have dw := Li2_comp (((hasDerivAt_id t).mul_const y).div
        (((hasDerivAt_id t).const_sub 1).mul_const (1 - y)) (mul_ne_zero hA.ne' hB.ne'))
        (by simp only [id]; exact hw1) (by simp only [id]; exact hw0)
      have dx := Li2_comp (hasDerivAt_id t) (by simp only [id]; linarith [hlt.1])
        (by simp only [id]; exact ht0)
      have dl := (((hasDerivAt_id t).const_sub 1).log hA.ne').mul_const (Real.log (1 - y))
      have e1 : Real.log (1 - t / (1 - y)) = Real.log (1 - t - y) - Real.log (1 - y) := by
        rw [show 1 - t / (1 - y) = (1 - t - y) / (1 - y) by (field_simp; ring),
          Real.log_div hL.ne' hB.ne']
      have e2 : Real.log (1 - y / (1 - t)) = Real.log (1 - t - y) - Real.log (1 - t) := by
        rw [show 1 - y / (1 - t) = (1 - t - y) / (1 - t) by field_simp,
          Real.log_div hL.ne' hA.ne']
      have e3 : Real.log (1 - t * y / ((1 - t) * (1 - y)))
          = Real.log (1 - t - y) - (Real.log (1 - t) + Real.log (1 - y)) := by
        rw [show 1 - t * y / ((1 - t) * (1 - y)) = (1 - t - y) / ((1 - t) * (1 - y)) by
            (field_simp; ring),
          Real.log_div hL.ne' (mul_ne_zero hA.ne' hB.ne'), Real.log_mul hA.ne' hB.ne']
      refine (((((du.add dv).sub dw).sub dx).sub (hasDerivAt_const t (Li2 y))).sub dl).congr_deriv ?_
      simp only [id, Pi.div_apply, e1, e2, e3]
      field_simp
      ring)
  simp [Li2_zero'] at key
  linarith

end D00


theorem D00.proof_1 : D00.stmt_1 := fun x => D00.intervalIntegrable_g 0 x
theorem D00.proof_2 : D00.stmt_2 := by simp [D00.stmt_2, D00.Li2]
theorem D00.proof_3 : D00.stmt_3 := D00.continuous_Li2.continuousOn
theorem D00.proof_4 : D00.stmt_4 := fun _ hx h0 => D00.hasDerivAt_Li2_ne hx.ne h0
theorem D00.proof_5 : D00.stmt_5 := by
  have := D00.hasDerivAt_Li2_G (x := 0) (by norm_num)
  simpa [D00.stmt_5, D00.G_zero] using this
theorem D00.proof_6 : D00.stmt_6 := by
  intro x hx
  rw [D00.Li2_eq_S hx, D00.S]
  exact (D00.summable_S hx).hasSum
theorem D00.proof_7 : D00.stmt_7 := by
  have h1 := D00.proof_6 1 (by norm_num)
  have h2 : HasSum (fun n : ℕ => (1 : ℝ) / ((n : ℝ) + 1) ^ 2) (Real.pi ^ 2 / 6) := by
    have := (hasSum_nat_add_iff' 1).mpr hasSum_zeta_two
    simpa using this
  simp only [one_pow] at h1
  exact h1.unique h2
theorem D00.proof_8 : D00.stmt_8 := D00.Li2_neg_one'
theorem D00.proof_9 : D00.stmt_9 := fun _ h0 h1 => D00.reflection' h0 h1
theorem D00.proof_10 : D00.stmt_10 := fun _ hx => D00.landen' hx
theorem D00.proof_11 : D00.stmt_11 := fun _ hx => D00.inversion_neg' hx
theorem D00.proof_12 : D00.stmt_12 := fun _ hx => D00.inversion_gt_one' hx
theorem D00.proof_13 : D00.stmt_13 := fun _ hx => D00.duplication' hx
theorem D00.proof_14 : D00.stmt_14 := by
  have h := D00.reflection' (x := 1 / 2) (by norm_num) (by norm_num)
  rw [show (1:ℝ) - 1 / 2 = 1 / 2 by norm_num, one_div, Real.log_inv] at h
  show D00.Li2 (1 / 2) = _
  rw [one_div]
  linear_combination h / 2
theorem D00.proof_15 : D00.stmt_15 := fun _ h0 h1 => D00.chi2_landen' h0 h1
theorem D00.proof_16 : D00.stmt_16 := fun _ _ hx hy hxy => D00.five_term' hx hy hxy

/-! ## Convenient restatements for importing files -/

namespace D00

open Real

theorem intervalIntegrable_integrand (x : ℝ) :
    IntervalIntegrable (fun t : ℝ => Real.log (1 - t) / t) MeasureTheory.volume 0 x :=
  proof_1 x

@[simp] theorem Li2_zero : Li2 0 = 0 := proof_2

theorem continuousOn_Li2 : ContinuousOn Li2 (Set.Iic 1) := proof_3

theorem hasDerivAt_Li2 {x : ℝ} (hx : x < 1) (h0 : x ≠ 0) :
    HasDerivAt Li2 (-Real.log (1 - x) / x) x := proof_4 x hx h0

theorem hasDerivAt_Li2_zero : HasDerivAt Li2 1 0 := proof_5

theorem hasSum_Li2 {x : ℝ} (hx : |x| ≤ 1) :
    HasSum (fun n : ℕ => x ^ (n + 1) / ((n : ℝ) + 1) ^ 2) (Li2 x) := proof_6 x hx

@[simp] theorem Li2_one : Li2 1 = π ^ 2 / 6 := proof_7

@[simp] theorem Li2_neg_one : Li2 (-1) = -(π ^ 2 / 12) := proof_8

theorem reflection {x : ℝ} (h0 : 0 < x) (h1 : x < 1) :
    Li2 x + Li2 (1 - x) = π ^ 2 / 6 - Real.log x * Real.log (1 - x) := proof_9 x h0 h1

theorem landen {x : ℝ} (hx : x < 1) :
    Li2 x + Li2 (x / (x - 1)) = -(1 / 2) * Real.log (1 - x) ^ 2 := proof_10 x hx

theorem inversion_neg {x : ℝ} (hx : x < 0) :
    Li2 x + Li2 (1 / x) = -(π ^ 2 / 6) - (1 / 2) * Real.log (-x) ^ 2 := proof_11 x hx

theorem inversion_gt_one {x : ℝ} (hx : 1 < x) :
    Li2 x + Li2 (1 / x) = π ^ 2 / 3 - (1 / 2) * Real.log x ^ 2 := proof_12 x hx

theorem duplication {x : ℝ} (hx : |x| ≤ 1) :
    Li2 x + Li2 (-x) = (1 / 2) * Li2 (x ^ 2) := proof_13 x hx

theorem Li2_half : Li2 (1 / 2) = π ^ 2 / 12 - (1 / 2) * Real.log 2 ^ 2 := proof_14

theorem chi2_landen {x : ℝ} (h0 : 0 < x) (h1 : x < 1) :
    chi2 ((1 - x) / (1 + x)) + chi2 x
      = π ^ 2 / 8 + (1 / 2) * Real.log x * Real.log ((1 + x) / (1 - x)) := proof_15 x h0 h1

theorem five_term {x y : ℝ} (hx : x < 1) (hy : y < 1) (hxy : x + y < 1) :
    Li2 (x / (1 - y)) + Li2 (y / (1 - x)) - Li2 (x * y / ((1 - x) * (1 - y)))
      = Li2 x + Li2 y + Real.log (1 - x) * Real.log (1 - y) := proof_16 x y hx hy hxy

end D00
