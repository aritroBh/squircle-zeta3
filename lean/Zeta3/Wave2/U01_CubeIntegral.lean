import Mathlib
/-!
# U01 — the cube-end integral `I` behind `prop:cuberadius`

Source: `paper/ieee/main.tex`, Proposition `prop:cuberadius` (lines 1960-1966) and the proof
sketch right after it (line 1968).

**Proposition, lines 1961-1964:**

> In $\varepsilon$, $S$ has a simple pole at $-1$ ($p=-1$) with residue $48R$, where
> \[
> R=-1-\frac{2\pi^{3/2}}{\Gamma(\frac14)^2}=-1.8472130847939790866\ldots
> \]
> is the residue of $V$ there, and $E$ is regular there.

**Proof sketch, line 1968:**

> Its residue is $I-2$, where $I=\int_0^1(\sqrt{1+a^4}-1)\,da/a^2=1-2\pi^{3/2}/\Gamma(\frac14)^2$
> by integration by parts and a Beta integral: the integral itself gives $I-1$, and a boundary
> layer at $a=0$ adds a point mass $-1$.

Handling of `a = 0`: the integrand is `f a = (√(1+a⁴) - 1) / a²`. For `a ≠ 0` this is the paper's
integrand; as `a → 0` it behaves like `a²/2 → 0`, so its continuous extension takes the value `0`
at `a = 0`. In Lean `f 0 = 0 / 0 = 0`, which is exactly that extension value, so `f` *is* the
continuous extension on `[0,1]` (and `stmt_1` says so). The value at a single point does not
affect the integral anyway. Integrability is stated explicitly in `stmt_2`, so the integral
identity cannot hold by the `intervalIntegral`-of-non-integrable-is-`0` convention.

`π^{3/2}` is `Real.pi ^ (3/2 : ℝ)` (rpow of a positive base) and `Γ(1/4)` is `Real.Gamma (1/4)`
(positive argument, so no convention issue).

* `stmt_1`: `f` is continuous on `[0,1]` (the integrand extends continuously to `a = 0`).
* `stmt_2`: `f` is interval integrable on `[0,1]` and `I = ∫₀¹ f = 1 - 2π^{3/2}/Γ(1/4)²`.
* `stmt_3`: `R := I - 2` equals the printed closed form `-1 - 2π^{3/2}/Γ(1/4)²`.
* `stmt_4`: the closed form of `R` lies strictly between `-1.8472130847939790867` and
  `-1.8472130847939790866` (all 19 printed decimals; the true value is
  `-1.84721308479397908660649…`, checked via `2π^{3/2}/Γ(1/4)² = AGM(√2,1)/√2`). A proof needs
  `Γ(1/4)` to about 20 significant digits, which Mathlib does not provide directly.

NOT COVERED:
* The holomorphy of `E` and `V` on `|ε| < 1`, the radius of convergence one, the simple pole of
  `S` at `ε = -1` with residue `48R`, the claim that `I - 2` *is* the residue of `V`, the
  boundary-layer point mass `-1`, and the regularity of `E` at `ε = -1` (the paper's `E`, `V` are
  only described as "explicit in our notes"; the task excludes these).
* The intermediate claim "the integral itself gives `I - 1`" (it refers to the unstated
  residue integral of `V`).
-/
-- BEGIN STATEMENTS
namespace U01

/-- The integrand `(√(1+a⁴) - 1)/a²` of `I` (proof sketch of `prop:cuberadius`, line 1968).
At `a = 0` Lean gives `0/0 = 0`, the value of the continuous extension. -/
noncomputable def f (a : ℝ) : ℝ := (Real.sqrt (1 + a ^ 4) - 1) / a ^ 2

/-- `I = ∫_0^1 (√(1+a⁴) - 1) da/a²` (line 1968). -/
noncomputable def I : ℝ := ∫ a in (0 : ℝ)..1, f a

/-- `R = I - 2`, "Its residue is $I-2$" (line 1968). -/
noncomputable def R : ℝ := I - 2

/-- The integrand extends continuously to `a = 0` (with value `0`, which is `f 0` in Lean),
so it is continuous on the closed interval `[0,1]`. -/
def stmt_1 : Prop := ContinuousOn f (Set.Icc 0 1)

/-- "$I=\int_0^1(\sqrt{1+a^4}-1)\,da/a^2=1-2\pi^{3/2}/\Gamma(\frac14)^2$" (line 1968),
with integrability stated. -/
def stmt_2 : Prop :=
  IntervalIntegrable f MeasureTheory.volume 0 1 ∧
    I = 1 - 2 * Real.pi ^ (3 / 2 : ℝ) / Real.Gamma (1 / 4) ^ 2

/-- "Its residue is $I-2$" together with "$R=-1-\frac{2\pi^{3/2}}{\Gamma(\frac14)^2}$"
(lines 1963 and 1968): `I - 2` equals the printed closed form of `R`. -/
def stmt_3 : Prop := R = -1 - 2 * Real.pi ^ (3 / 2 : ℝ) / Real.Gamma (1 / 4) ^ 2

/-- "$R=-1-\frac{2\pi^{3/2}}{\Gamma(\frac14)^2}=-1.8472130847939790866\ldots$" (line 1963):
all 19 printed decimals of the closed form. -/
def stmt_4 : Prop :=
  (-1.8472130847939790867 : ℝ) < -1 - 2 * Real.pi ^ (3 / 2 : ℝ) / Real.Gamma (1 / 4) ^ 2 ∧
    -1 - 2 * Real.pi ^ (3 / 2 : ℝ) / Real.Gamma (1 / 4) ^ 2 < (-1.8472130847939790866 : ℝ)

end U01
-- END STATEMENTS
namespace U01Aux

open Set intervalIntegral MeasureTheory

lemma f_eq (a : ℝ) : U01.f a = a ^ 2 / (Real.sqrt (1 + a ^ 4) + 1) := by
  have hpos : 0 < 1 + a^4 := by positivity
  have hs : Real.sqrt (1 + a^4) ^ 2 = 1 + a^4 := Real.sq_sqrt hpos.le
  have hs1 : 0 < Real.sqrt (1 + a^4) + 1 := by positivity
  unfold U01.f
  by_cases ha : a = 0
  · subst ha; simp
  · field_simp
    linear_combination hs

lemma f_cont : Continuous U01.f := by
  have : U01.f = fun a => a ^ 2 / (Real.sqrt (1 + a ^ 4) + 1) := funext f_eq
  rw [this]
  have h1 : Continuous fun a : ℝ => a ^ 2 := by fun_prop
  have h2 : Continuous fun a : ℝ => Real.sqrt (1 + a ^ 4) + 1 := by fun_prop
  exact h1.div h2 (fun a => by positivity)

noncomputable def H (a : ℝ) : ℝ :=
  -(a ^ 3 / (Real.sqrt (1 + a ^ 4) + 1)) + 2 * (a ^ 3 / Real.sqrt (1 + a ^ 4))

lemma H_deriv (a : ℝ) : HasDerivAt H (a ^ 2 / (Real.sqrt (1 + a ^ 4) + 1)
    + 4 * (a ^ 2 / Real.sqrt (1 + a ^ 4) ^ 3)) a := by
  have hpos : 0 < 1 + a^4 := by positivity
  have hs : Real.sqrt (1 + a^4) ^ 2 = 1 + a^4 := Real.sq_sqrt hpos.le
  have hs0 : 0 < Real.sqrt (1 + a^4) := Real.sqrt_pos.2 hpos
  have hsd : HasDerivAt (fun a : ℝ => Real.sqrt (1 + a ^ 4)) (4 * a ^ 3 / (2 * Real.sqrt (1 + a ^ 4))) a := by
    have h1 : HasDerivAt (fun a : ℝ => 1 + a ^ 4) (4 * a ^ 3) a := by
      simpa using ((hasDerivAt_pow 4 a).const_add 1)
    exact h1.sqrt hpos.ne'
  have h3 : HasDerivAt (fun a : ℝ => a ^ 3) (3 * a ^ 2) a := by
    simpa using hasDerivAt_pow 3 a
  have hA := h3.div (hsd.add_const 1) (by positivity)
  have hB := h3.div hsd hs0.ne'
  have hC := hA.neg.add (hB.const_mul 2)
  convert hC using 1
  · funext x; simp [H, Pi.div_apply]
  generalize Real.sqrt (1 + a ^ 4) = s at *
  have hs1 : 0 < s + 1 := by positivity
  field_simp
  linear_combination (-(a^2)*(8*(s+1)^2-4*s^2)) * hs

lemma g3_cont : Continuous fun a : ℝ => a ^ 2 / Real.sqrt (1 + a ^ 4) ^ 3 := by
  have h1 : Continuous fun a : ℝ => a ^ 2 := by fun_prop
  have h2 : Continuous fun a : ℝ => Real.sqrt (1 + a ^ 4) ^ 3 := by fun_prop
  exact h1.div h2 (fun a => by have : 0 < Real.sqrt (1 + a^4) := Real.sqrt_pos.2 (by positivity); positivity)

lemma I_eq : U01.I = 1 - 4 * ∫ a in (0:ℝ)..1, a ^ 2 / Real.sqrt (1 + a ^ 4) ^ 3 := by
  have h1 : Continuous fun a : ℝ => a ^ 2 / (Real.sqrt (1 + a ^ 4) + 1) := by
    have : U01.f = fun a => a ^ 2 / (Real.sqrt (1 + a ^ 4) + 1) := funext f_eq
    rw [← this]; exact f_cont
  have hftc := intervalIntegral.integral_eq_sub_of_hasDerivAt (a := 0) (b := 1)
    (f := H) (f' := fun a => a ^ 2 / (Real.sqrt (1 + a ^ 4) + 1)
    + 4 * (a ^ 2 / Real.sqrt (1 + a ^ 4) ^ 3)) (fun x _ => H_deriv x)
    ((h1.add (g3_cont.const_mul 4)).intervalIntegrable _ _)
  have hH1 : H 1 = 1 := by
    have h2 : Real.sqrt 2 ^ 2 = 2 := Real.sq_sqrt (by norm_num)
    have h2p : 0 < Real.sqrt 2 := by positivity
    unfold H
    norm_num
    field_simp
    nlinarith
  have hH0 : H 0 = 0 := by simp [H]
  rw [hH1, hH0, intervalIntegral.integral_add (h1.intervalIntegrable _ _)
    ((g3_cont.const_mul 4).intervalIntegrable _ _), intervalIntegral.integral_const_mul] at hftc
  have : U01.I = ∫ a in (0:ℝ)..1, a ^ 2 / (Real.sqrt (1 + a ^ 4) + 1) := by
    unfold U01.I; simp_rw [f_eq]
  rw [this]; linarith
end U01Aux

namespace U01Aux
open Set intervalIntegral MeasureTheory

noncomputable def gb (u : ℝ) : ℝ := u ^ (-(1/4 : ℝ)) * (1 - u) ^ (-(1/4 : ℝ))

noncomputable def phi (a : ℝ) : ℝ := a ^ 4 / (1 + a ^ 4)

noncomputable def phi' (a : ℝ) : ℝ := 4 * a ^ 3 / (1 + a ^ 4) ^ 2

lemma phi_deriv (a : ℝ) : HasDerivAt phi (phi' a) a := by
  have hpos : 0 < 1 + a^4 := by positivity
  have h1 : HasDerivAt (fun a : ℝ => a ^ 4) (4 * a ^ 3) a := by
    simpa using hasDerivAt_pow 4 a
  have h2 : HasDerivAt (fun a : ℝ => 1 + a ^ 4) (4 * a ^ 3) a := by
    simpa using h1.const_add 1
  have h3 := h1.div h2 hpos.ne'
  convert h3 using 1
  · funext x; simp [phi, Pi.div_apply]
  unfold phi'
  field_simp
  ring

lemma gb_comp (a : ℝ) (ha : 0 ≤ a) :
    (1/4 : ℝ) * gb (phi a) * phi' a = a ^ 2 / Real.sqrt (1 + a ^ 4) ^ 3 := by
  have hpos : 0 < 1 + a^4 := by positivity
  rcases ha.eq_or_lt with h0 | hapos
  · subst h0
    simp [gb, phi, phi']
  · set w : ℝ := (1 + a ^ 4) ^ (1/4 : ℝ) with hw
    have hwpos : 0 < w := Real.rpow_pos_of_pos hpos _
    have hw4 : w ^ 4 = 1 + a ^ 4 := by
      rw [hw, ← Real.rpow_natCast, ← Real.rpow_mul hpos.le]; norm_num
    have hs : Real.sqrt (1 + a ^ 4) = w ^ 2 := by
      rw [← hw4, show w ^ 4 = (w ^ 2) ^ 2 by ring, Real.sqrt_sq (by positivity)]
    have hphi : phi a = (a / w) ^ 4 := by
      unfold phi; rw [div_pow, hw4]
    have h1phi : 1 - phi a = (1 / w) ^ 4 := by
      rw [hphi, div_pow, one_div, inv_pow, hw4]
      field_simp; ring
    have e1 : (phi a) ^ (-(1/4 : ℝ)) = w / a := by
      rw [hphi, Real.rpow_neg (by positivity), show (1/4:ℝ) = ((4:ℕ):ℝ)⁻¹ by norm_num,
        Real.pow_rpow_inv_natCast (by positivity) (by norm_num)]
      field_simp
    have e2 : (1 - phi a) ^ (-(1/4 : ℝ)) = w := by
      rw [h1phi, Real.rpow_neg (by positivity), show (1/4:ℝ) = ((4:ℕ):ℝ)⁻¹ by norm_num,
        Real.pow_rpow_inv_natCast (by positivity) (by norm_num)]
      field_simp
    unfold gb
    rw [e1, e2, hs]
    unfold phi'
    rw [← hw4]
    field_simp

end U01Aux

namespace U01Aux
open Set intervalIntegral MeasureTheory

lemma gb_cplx (x : ℝ) (hx : x ∈ Icc (0:ℝ) 1) :
    ((gb x : ℝ) : ℂ) = (x:ℂ) ^ ((3/4:ℂ) - 1) * (1 - (x:ℂ)) ^ ((3/4:ℂ) - 1) := by
  have e : ((3/4:ℂ) - 1) = (((-(1/4:ℝ)) : ℝ) : ℂ) := by push_cast; ring
  rw [e]
  unfold gb
  have h1 : (0:ℝ) ≤ 1 - x := by linarith [hx.2]
  rw [Complex.ofReal_mul, Complex.ofReal_cpow hx.1, Complex.ofReal_cpow h1]
  push_cast
  rfl

lemma gb_int : IntervalIntegrable gb volume 0 1 := by
  have hc := Complex.betaIntegral_convergent (u := 3/4) (v := 3/4) (by norm_num) (by norm_num)
  rw [intervalIntegrable_iff_integrableOn_Icc_of_le (by norm_num)] at hc ⊢
  have hre := hc.re
  refine hre.congr_fun (fun x hx => ?_) measurableSet_Icc
  have := gb_cplx x hx
  show RCLike.re _ = _
  rw [← this]; simp

lemma beta_eq : Complex.betaIntegral (3/4) (3/4) = (((∫ u in (0:ℝ)..1, gb u) : ℝ) : ℂ) := by
  rw [Complex.betaIntegral, ← intervalIntegral.integral_ofReal]
  refine intervalIntegral.integral_congr (fun x hx => ?_)
  rw [uIcc_of_le (by norm_num)] at hx
  exact (gb_cplx x hx).symm

lemma gb_total : (∫ u in (0:ℝ)..1, gb u) =
    2 * Real.Gamma (3/4) ^ 2 / Real.sqrt Real.pi := by
  have h := Complex.Gamma_mul_Gamma_eq_betaIntegral (s := 3/4) (t := 3/4) (by norm_num) (by norm_num)
  rw [beta_eq] at h
  have e1 : (3/4 : ℂ) = (((3/4:ℝ)) : ℂ) := by push_cast; ring
  have e2 : (3/4 + 3/4 : ℂ) = (((3/2:ℝ)) : ℂ) := by push_cast; ring
  rw [e2, e1, Complex.Gamma_ofReal, Complex.Gamma_ofReal] at h
  rw [← Complex.ofReal_mul, ← Complex.ofReal_mul] at h
  have h' := Complex.ofReal_injective h
  have hg : Real.Gamma (3/2) = Real.sqrt Real.pi / 2 := by
    have : (3/2:ℝ) = 1/2 + 1 := by norm_num
    rw [this, Real.Gamma_add_one (by norm_num), Real.Gamma_one_half_eq]; ring
  rw [hg] at h'
  have hp : 0 < Real.sqrt Real.pi := Real.sqrt_pos.2 Real.pi_pos
  field_simp
  nlinarith [h']

end U01Aux

namespace U01Aux
open Set intervalIntegral MeasureTheory

lemma phi_cont : Continuous phi := continuous_iff_continuousAt.2 fun x => (phi_deriv x).continuousAt

lemma phi_mem (a : ℝ) (ha : a ∈ Icc (0:ℝ) 1) : phi a ∈ Icc (0:ℝ) (1/2) := by
  have h4 : 0 ≤ a ^ 4 := by positivity
  have h41 : a ^ 4 ≤ 1 := pow_le_one₀ ha.1 ha.2
  unfold phi
  constructor
  · positivity
  · rw [div_le_iff₀ (by positivity)]; linarith

lemma phi_mem' (a : ℝ) (ha : a ∈ Ioo (0:ℝ) 1) : phi a ∈ Ioo (0:ℝ) 1 := by
  have h4 : 0 < a ^ 4 := by have := ha.1; positivity
  unfold phi
  constructor
  · positivity
  · rw [div_lt_one (by positivity)]; linarith

lemma sub_eq : (∫ a in (0:ℝ)..1, a ^ 2 / Real.sqrt (1 + a ^ 4) ^ 3) =
    (1/4) * ∫ u in (0:ℝ)..(1/2), gb u := by
  have hgi : IntegrableOn gb (Icc (0:ℝ) 1) := by
    have := gb_int
    rwa [intervalIntegrable_iff_integrableOn_Icc_of_le (by norm_num)] at this
  have key := intervalIntegral.integral_comp_mul_deriv''' (a := 0) (b := 1) (f := phi) (f' := phi')
    (g := fun u => (1/4 : ℝ) * gb u)
    phi_cont.continuousOn
    (fun x _ => (phi_deriv x).hasDerivWithinAt)
    (by
      have hsub : phi '' Ioo (min (0:ℝ) 1) (max 0 1) ⊆ Ioo 0 1 := by
        rintro _ ⟨a, ha, rfl⟩
        simp only [min_eq_left zero_le_one, max_eq_right zero_le_one] at ha
        exact phi_mem' a ha
      have hc : ContinuousOn gb (Ioo (0:ℝ) 1) := by
        unfold gb
        refine ContinuousOn.mul ?_ ?_
        · exact ContinuousOn.rpow_const continuousOn_id (fun x hx => Or.inl hx.1.ne')
        · exact ContinuousOn.rpow_const (continuousOn_const.sub continuousOn_id)
            (fun x hx => Or.inl (by have := hx.2; linarith))
      exact (continuousOn_const.mul hc).mono hsub)
    (by
      have hsub : phi '' uIcc (0:ℝ) 1 ⊆ Icc 0 1 := by
        rintro _ ⟨a, ha, rfl⟩
        rw [uIcc_of_le zero_le_one] at ha
        have := phi_mem a ha
        exact ⟨this.1, by linarith [this.2]⟩
      exact IntegrableOn.mono_set (Integrable.const_mul hgi (1/4)) hsub)
    (by
      have : IntegrableOn (fun a : ℝ => a ^ 2 / Real.sqrt (1 + a ^ 4) ^ 3) (uIcc (0:ℝ) 1) :=
        g3_cont.continuousOn.integrableOn_uIcc
      refine this.congr_fun (fun a ha => ?_) measurableSet_uIcc
      rw [uIcc_of_le zero_le_one] at ha
      simp only [Function.comp]
      rw [← gb_comp a ha.1]) 
  have e0 : phi 0 = 0 := by simp [phi]
  have e1 : phi 1 = 1/2 := by norm_num [phi]
  rw [e0, e1] at key
  rw [intervalIntegral.integral_const_mul] at key
  rw [← key]
  refine intervalIntegral.integral_congr (fun a ha => ?_)
  rw [uIcc_of_le zero_le_one] at ha
  simp only [Function.comp]
  exact (gb_comp a ha.1).symm

lemma half_eq : (∫ u in (0:ℝ)..(1/2), gb u) = (1/2) * ∫ u in (0:ℝ)..1, gb u := by
  have hi := gb_int
  have h1 : IntervalIntegrable gb volume 0 (1/2) := hi.mono_set (by
    rw [uIcc_of_le (by norm_num), uIcc_of_le (by norm_num)]; exact Icc_subset_Icc le_rfl (by norm_num))
  have h2 : IntervalIntegrable gb volume (1/2) 1 := hi.mono_set (by
    rw [uIcc_of_le (by norm_num), uIcc_of_le (by norm_num)]; exact Icc_subset_Icc (by norm_num) le_rfl)
  have hs : (∫ u in (1/2:ℝ)..1, gb u) = ∫ u in (0:ℝ)..(1/2), gb u := by
    have this : (∫ x in (0:ℝ)..(1/2), gb (1 - x)) = ∫ x in (1 - 1/2 : ℝ)..(1 - 0 : ℝ), gb x :=
      intervalIntegral.integral_comp_sub_left (fun u => gb u) (1:ℝ)
    have e1 : (1:ℝ) - 1/2 = 1/2 := by norm_num
    rw [e1, sub_zero] at this
    rw [← this]
    refine intervalIntegral.integral_congr (fun x _ => ?_)
    simp only [gb, sub_sub_cancel]
    ring
  have := intervalIntegral.integral_add_adjacent_intervals h1 h2
  rw [hs] at this
  linarith

lemma gamma34 : Real.Gamma (3/4) ^ 2 = 2 * Real.pi ^ 2 / Real.Gamma (1/4) ^ 2 := by
  have h := Real.Gamma_mul_Gamma_one_sub (1/4)
  have e : (1:ℝ) - 1/4 = 3/4 := by norm_num
  rw [e] at h
  have hsin : Real.sin (Real.pi * (1/4)) = Real.sqrt 2 / 2 := by
    rw [show Real.pi * (1/4) = Real.pi / 4 by ring]; exact Real.sin_pi_div_four
  rw [hsin] at h
  have hg : 0 < Real.Gamma (1/4) := Real.Gamma_pos_of_pos (by norm_num)
  have h2 : Real.sqrt 2 ^ 2 = 2 := Real.sq_sqrt (by norm_num)
  have h2p : 0 < Real.sqrt 2 := by positivity
  have hG : Real.Gamma (3/4) = Real.pi * Real.sqrt 2 / Real.Gamma (1/4) := by
    field_simp at h ⊢
    have h5 : (Real.Gamma (1/4) * Real.Gamma (3/4) - Real.pi * Real.sqrt 2) * Real.sqrt 2 = 0 := by
      linear_combination h - Real.pi * h2
    rcases mul_eq_zero.1 h5 with h6 | h6
    · linarith
    · exact absurd h6 h2p.ne'
  rw [hG, div_pow, mul_pow, h2]; ring

end U01Aux

namespace U01Aux
open Set intervalIntegral MeasureTheory

noncomputable def cc : ℕ → ℝ
  | 0 => 1
  | (n+1) => cc n * ((n:ℝ) + 1/4) / ((n:ℝ) + 1)

lemma cc_pos (n : ℕ) : 0 < cc n := by
  induction n with
  | zero => simp [cc]
  | succ n ih => simp only [cc]; positivity

noncomputable def PP (N : ℕ) (u : ℝ) : ℝ := ∑ n ∈ Finset.range N, cc n * u ^ n

noncomputable def QQ (N : ℕ) (u : ℝ) : ℝ := (1 - u) ^ (1/4 : ℝ) * PP N u

lemma w_id (u : ℝ) (hu : u < 1) : (1 - u) ^ (1/4 : ℝ) = (1 - u) * (1 - u) ^ (-(3/4) : ℝ) := by
  have h : 0 < 1 - u := by linarith
  conv_lhs => rw [show (1/4:ℝ) = 1 + (-(3/4)) by norm_num, Real.rpow_add h, Real.rpow_one]

lemma QQ_deriv (N : ℕ) (u : ℝ) (hu : u < 1) :
    HasDerivAt (QQ (N+1)) (-(((N:ℝ)+1) * cc (N+1)) * u ^ N * (1 - u) ^ (-(3/4) : ℝ)) u := by
  have h1 : 0 < 1 - u := by linarith
  have hr : HasDerivAt (fun u : ℝ => (1 - u) ^ (1/4 : ℝ))
      (-(1/4) * (1 - u) ^ (-(3/4) : ℝ)) u := by
    have := ((hasDerivAt_id u).const_sub 1).rpow_const (p := 1/4) (Or.inl h1.ne')
    refine this.congr_deriv ?_
    rw [show (1/4:ℝ) - 1 = -(3/4) by norm_num]
    simp only [id_eq]
    ring
  induction N with
  | zero =>
    have : QQ 1 = fun u : ℝ => (1 - u) ^ (1/4 : ℝ) := by
      funext v; simp [QQ, PP, cc]
    rw [this]
    convert hr using 1
    simp [cc]
  | succ N ih =>
    have hq : QQ (N + 1 + 1) = fun v => QQ (N+1) v + cc (N+1) * v ^ (N+1) * (1 - v) ^ (1/4 : ℝ) := by
      funext v; simp only [QQ, PP, Finset.sum_range_succ]; ring
    rw [hq]
    have hp : HasDerivAt (fun v : ℝ => v ^ (N+1)) (((N+1:ℕ):ℝ) * u ^ N) u := by
      simpa using hasDerivAt_pow (N+1) u
    have hs := ((hp.const_mul (cc (N+1))).mul hr)
    have := ih.add hs
    convert this using 1
    rw [w_id u hu]
    simp only [cc]
    push_cast
    field_simp
    ring
end U01Aux

namespace U01Aux
open Set intervalIntegral MeasureTheory

lemma rp_le_two (u p : ℝ) (hu0 : 0 ≤ u) (hu : u ≤ 1/2) (hp1 : p ≤ 1) :
    (1 - u) ^ (-p) ≤ 2 := by
  have h0 : 0 < 1 - u := by linarith
  have h1 : 1 - u ≤ 1 := by linarith
  calc (1 - u) ^ (-p) ≤ (1 - u) ^ (-1 : ℝ) :=
        Real.rpow_le_rpow_of_exponent_ge h0 h1 (by linarith)
    _ = (1 - u)⁻¹ := by rw [Real.rpow_neg_one]
    _ ≤ 2 := by rw [inv_le_comm₀ h0 (by norm_num)]; linarith

lemma PP_zero (N : ℕ) : PP (N+1) 0 = 1 := by
  simp [PP, Finset.sum_range_succ', cc]

lemma E_bounds (M : ℕ) (u : ℝ) (hu : u ∈ Icc (0:ℝ) (1/2)) :
    0 ≤ 1 - QQ (M+1) u ∧ 1 - QQ (M+1) u ≤ 2 * cc (M+1) * u ^ (M+1) := by
  have hD : Convex ℝ (Icc (0:ℝ) (1/2)) := convex_Icc (0:ℝ) (1/2)
  have hd : ∀ v ∈ Icc (0:ℝ) (1/2), HasDerivAt (fun v => 1 - QQ (M+1) v)
      (((M:ℝ)+1) * cc (M+1) * v ^ M * (1 - v) ^ (-(3/4) : ℝ)) v := by
    intro v hv
    have := (QQ_deriv M v (by linarith [hv.2])).const_sub 1
    refine this.congr_deriv ?_
    ring
  have hQ0 : QQ (M+1) 0 = 1 := by simp [QQ, PP_zero]
  constructor
  · have hm : MonotoneOn (fun v => 1 - QQ (M+1) v) (Icc (0:ℝ) (1/2)) := by
      refine monotoneOn_of_hasDerivWithinAt_nonneg hD
        (fun y hy => (hd y hy).continuousAt.continuousWithinAt)
        (fun y hy => (hd y (interior_subset hy)).hasDerivWithinAt) ?_
      intro y hy
      simp only [interior_Icc, Set.mem_Ioo] at hy
      have := cc_pos (M+1)
      have h4 : 0 ≤ y ^ M := pow_nonneg hy.1.le _
      have hw : 0 ≤ (1 - y) ^ (-(3/4) : ℝ) := Real.rpow_nonneg (by linarith [hy.2]) _
      positivity
    have := hm ⟨le_rfl, by norm_num⟩ hu hu.1
    simpa [hQ0] using this
  · have hm : MonotoneOn (fun v => 2 * cc (M+1) * v ^ (M+1) - (1 - QQ (M+1) v)) (Icc (0:ℝ) (1/2)) := by
      refine monotoneOn_of_hasDerivWithinAt_nonneg hD
        (f' := fun v => ((M:ℝ)+1) * cc (M+1) * v ^ M * (2 - (1 - v) ^ (-(3/4) : ℝ)))
        (fun y hy => ?_) (fun y hy => ?_) ?_
      · have h1 : HasDerivAt (fun v : ℝ => 2 * cc (M+1) * v ^ (M+1) - (1 - QQ (M+1) v))
            (((M:ℝ)+1) * cc (M+1) * y ^ M * (2 - (1 - y) ^ (-(3/4) : ℝ))) y := by
          have hp : HasDerivAt (fun v : ℝ => v ^ (M+1)) (((M+1:ℕ):ℝ) * y ^ M) y := by
            simpa using hasDerivAt_pow (M+1) y
          have := (hp.const_mul (2 * cc (M+1))).sub (hd y hy)
          refine this.congr_deriv ?_
          push_cast; ring
        exact h1.continuousAt.continuousWithinAt
      · have hy' := interior_subset hy
        have h1 : HasDerivAt (fun v : ℝ => 2 * cc (M+1) * v ^ (M+1) - (1 - QQ (M+1) v))
            (((M:ℝ)+1) * cc (M+1) * y ^ M * (2 - (1 - y) ^ (-(3/4) : ℝ))) y := by
          have hp : HasDerivAt (fun v : ℝ => v ^ (M+1)) (((M+1:ℕ):ℝ) * y ^ M) y := by
            simpa using hasDerivAt_pow (M+1) y
          have := (hp.const_mul (2 * cc (M+1))).sub (hd y hy')
          refine this.congr_deriv ?_
          push_cast; ring
        exact h1.hasDerivWithinAt
      · intro y hy
        simp only [interior_Icc, Set.mem_Ioo] at hy
        have := cc_pos (M+1)
        have h2 := rp_le_two y (3/4) hy.1.le hy.2.le (by norm_num)
        have h3 : 0 ≤ 2 - (1 - y) ^ (-(3/4) : ℝ) := by linarith
        have h4 : 0 ≤ y ^ M := pow_nonneg hy.1.le _
        exact mul_nonneg (mul_nonneg (by positivity) h4) h3
    have := hm ⟨le_rfl, by norm_num⟩ hu hu.1
    simp [hQ0] at this
    linarith

lemma y_bounds (M : ℕ) (u : ℝ) (hu : u ∈ Icc (0:ℝ) (1/2)) :
    PP (M+1) u ≤ (1 - u) ^ (-(1/4) : ℝ) ∧
    (1 - u) ^ (-(1/4) : ℝ) ≤ PP (M+1) u + 4 * cc (M+1) * u ^ (M+1) := by
  obtain ⟨hE0, hE1⟩ := E_bounds M u hu
  have h0 : 0 < 1 - u := by linarith [hu.2]
  set y := (1 - u) ^ (-(1/4) : ℝ) with hy
  set r := (1 - u) ^ (1/4 : ℝ) with hr
  have hyr : y * r = 1 := by
    rw [hy, hr, ← Real.rpow_add h0]; norm_num
  have hr0 : 0 < r := Real.rpow_pos_of_pos h0 _
  have hy0 : 0 < y := Real.rpow_pos_of_pos h0 _
  have hy2 : y ≤ 2 := rp_le_two u (1/4) hu.1 hu.2 (by norm_num)
  have hQ : QQ (M+1) u = r * PP (M+1) u := rfl
  have hdiff : y - PP (M+1) u = y * (1 - QQ (M+1) u) := by
    rw [hQ]; linear_combination (PP (M+1) u) * hyr
  have hc := cc_pos (M+1)
  have hup : 0 ≤ u ^ (M+1) := pow_nonneg hu.1 _
  constructor
  · nlinarith
  · have : y * (1 - QQ (M+1) u) ≤ 2 * (2 * cc (M+1) * u ^ (M+1)) := by
      apply mul_le_mul hy2 hE1 hE0 (by norm_num)
    nlinarith
end U01Aux

namespace U01Aux
open Set intervalIntegral MeasureTheory

noncomputable def ss : ℝ := (1/2 : ℝ) ^ (3/4 : ℝ)

lemma tm (n : ℕ) (u : ℝ) (hu : 0 ≤ u) :
    u ^ (-(1/4) : ℝ) * u ^ n = u ^ ((n:ℝ) - 1/4) := by
  have hne : (-(1/4) : ℝ) + n ≠ 0 := by
    rcases Nat.eq_zero_or_pos n with h | h
    · subst h; norm_num
    · have : (1:ℝ) ≤ n := by exact_mod_cast h
      intro h'; linarith
  rw [← Real.rpow_natCast, ← Real.rpow_add' hu hne]
  ring_nf

lemma ii_term (n : ℕ) : IntervalIntegrable (fun u : ℝ => u ^ (-(1/4) : ℝ) * u ^ n) volume 0 (1/2) := by
  have h : IntervalIntegrable (fun u : ℝ => u ^ ((n:ℝ) - 1/4)) volume 0 (1/2) := by
    apply intervalIntegral.intervalIntegrable_rpow'
    have : (0:ℝ) ≤ n := Nat.cast_nonneg n
    linarith
  rw [intervalIntegrable_iff_integrableOn_Icc_of_le (by norm_num)] at h ⊢
  exact h.congr_fun (fun u hu => (tm n u hu.1).symm) measurableSet_Icc

lemma moment (n : ℕ) : (∫ u in (0:ℝ)..(1/2), u ^ (-(1/4) : ℝ) * u ^ n) = ss * (1/2) ^ n / ((n:ℝ) + 3/4) := by
  have h1 : (∫ u in (0:ℝ)..(1/2), u ^ (-(1/4) : ℝ) * u ^ n) = ∫ u in (0:ℝ)..(1/2), u ^ ((n:ℝ) - 1/4) := by
    refine intervalIntegral.integral_congr (fun u hu => ?_)
    rw [uIcc_of_le (by norm_num)] at hu
    exact tm n u hu.1
  have hn : (-1:ℝ) < (n:ℝ) - 1/4 := by
    have : (0:ℝ) ≤ n := Nat.cast_nonneg n
    linarith
  rw [h1, integral_rpow (Or.inl hn)]
  have e : (n:ℝ) - 1/4 + 1 = n + 3/4 := by ring
  rw [e, Real.zero_rpow (by positivity), sub_zero, Real.rpow_add (by norm_num), Real.rpow_natCast]
  unfold ss
  ring

end U01Aux

namespace U01Aux
open Set intervalIntegral MeasureTheory

noncomputable def Lo (M : ℕ) (u : ℝ) : ℝ := ∑ n ∈ Finset.range (M+1), cc n * (u ^ (-(1/4) : ℝ) * u ^ n)

lemma Lo_int (M : ℕ) : IntervalIntegrable (Lo M) volume 0 (1/2) := by
  have := IntervalIntegrable.sum (Finset.range (M+1)) (fun n _ => (ii_term n).const_mul (cc n))
  convert this using 1
  funext u
  simp [Lo]

lemma Lo_integral (M : ℕ) : (∫ u in (0:ℝ)..(1/2), Lo M u) =
    ∑ n ∈ Finset.range (M+1), cc n * (ss * (1/2) ^ n / ((n:ℝ) + 3/4)) := by
  unfold Lo
  rw [intervalIntegral.integral_of_le (by norm_num),
    MeasureTheory.integral_finsetSum _ (fun n _ =>
      ((ii_term n).const_mul (cc n)).1)]
  refine Finset.sum_congr rfl (fun n _ => ?_)
  rw [← intervalIntegral.integral_of_le (by norm_num), intervalIntegral.integral_const_mul, moment]

lemma Lo_eq (M : ℕ) (u : ℝ) : Lo M u = u ^ (-(1/4) : ℝ) * PP (M+1) u := by
  unfold Lo PP
  rw [Finset.mul_sum]
  exact Finset.sum_congr rfl (fun n _ => by ring)

lemma half_int : IntervalIntegrable gb volume 0 (1/2) :=
  gb_int.mono_set (by
    rw [uIcc_of_le (by norm_num), uIcc_of_le (by norm_num)]; exact Icc_subset_Icc le_rfl (by norm_num))

lemma J_lo (M : ℕ) : ∑ n ∈ Finset.range (M+1), cc n * (ss * (1/2) ^ n / ((n:ℝ) + 3/4))
    ≤ ∫ u in (0:ℝ)..(1/2), gb u := by
  rw [← Lo_integral]
  apply intervalIntegral.integral_mono_on (by norm_num) (Lo_int M) half_int
  intro u hu
  rw [Lo_eq]
  unfold gb
  have h0 : 0 ≤ u ^ (-(1/4) : ℝ) := Real.rpow_nonneg hu.1 _
  exact mul_le_mul_of_nonneg_left (y_bounds M u hu).1 h0

lemma J_hi (M : ℕ) : ∫ u in (0:ℝ)..(1/2), gb u ≤
    ∑ n ∈ Finset.range (M+1), cc n * (ss * (1/2) ^ n / ((n:ℝ) + 3/4))
    + 4 * cc (M+1) * (ss * (1/2) ^ (M+1) / (((M+1:ℕ):ℝ) + 3/4)) := by
  have hup : IntervalIntegrable (fun u : ℝ => Lo M u + 4 * cc (M+1) * (u ^ (-(1/4) : ℝ) * u ^ (M+1)))
      volume 0 (1/2) := (Lo_int M).add ((ii_term (M+1)).const_mul _)
  have hi := intervalIntegral.integral_mono_on (by norm_num) half_int hup (fun u hu => by
    rw [Lo_eq]
    unfold gb
    have h0 : 0 ≤ u ^ (-(1/4) : ℝ) := Real.rpow_nonneg hu.1 _
    have := mul_le_mul_of_nonneg_left (y_bounds M u hu).2 h0
    linarith)
  rw [intervalIntegral.integral_add (Lo_int M) ((ii_term (M+1)).const_mul _),
    intervalIntegral.integral_const_mul, moment, Lo_integral] at hi
  exact hi

end U01Aux

namespace U01Aux
open Set intervalIntegral MeasureTheory

def ccq : ℕ → ℚ
  | 0 => 1
  | (n+1) => ccq n * ((n:ℚ) + 1/4) / ((n:ℚ) + 1)

def Sq (m : ℕ) : ℚ := ∑ n ∈ Finset.range m, ccq n * ((1/2)^n / ((n:ℚ) + 3/4))

def Tq (m : ℕ) : ℚ := 4 * ccq m * ((1/2)^m / ((m:ℚ) + 3/4))

lemma cc_eq (n : ℕ) : cc n = (ccq n : ℝ) := by
  induction n with
  | zero => simp [cc, ccq]
  | succ n ih => simp only [cc, ccq, ih]; push_cast; ring

lemma lo_cast (m : ℕ) : ∑ n ∈ Finset.range m, cc n * (ss * (1/2) ^ n / ((n:ℝ) + 3/4)) = ss * (Sq m : ℝ) := by
  unfold Sq
  push_cast
  rw [Finset.mul_sum]
  refine Finset.sum_congr rfl (fun n _ => ?_)
  rw [cc_eq]; ring

lemma hi_cast (m : ℕ) : 4 * cc m * (ss * (1/2) ^ m / ((m:ℝ) + 3/4)) = ss * (Tq m : ℝ) := by
  unfold Tq
  push_cast
  rw [cc_eq]; ring

lemma ss_pos : 0 < ss := Real.rpow_pos_of_pos (by norm_num) _

lemma ss_pow4 : ss ^ 4 = 1/8 := by
  unfold ss
  rw [← Real.rpow_natCast, ← Real.rpow_mul (by norm_num)]
  norm_num

lemma ss_lo : ((594603557501360533358749 : ℚ) / 10^24 : ℚ) < ss := by
  by_contra h
  push Not at h
  have h4 := pow_le_pow_left₀ ss_pos.le h 4
  rw [ss_pow4] at h4
  have : (((594603557501360533358749 : ℚ) / 10^24)^4 : ℚ) < 1/8 := by norm_num
  have h5 : ((((594603557501360533358749 : ℚ) / 10^24)^4 : ℚ) : ℝ) < ((1/8 : ℚ) : ℝ) := by
    exact_mod_cast this
  push_cast at h5
  linarith

lemma ss_hi : ss < ((594603557501360533358751 : ℚ) / 10^24 : ℚ) := by
  by_contra h
  push Not at h
  have h4 := pow_le_pow_left₀ (by positivity) h 4
  rw [ss_pow4] at h4
  have : (1/8 : ℚ) < (((594603557501360533358751 : ℚ) / 10^24)^4 : ℚ) := by norm_num
  have h5 : ((1/8 : ℚ) : ℝ) < ((((594603557501360533358751 : ℚ) / 10^24)^4 : ℚ) : ℝ) := by
    exact_mod_cast this
  push_cast at h5
  linarith

theorem num_Sq : (8472130847939790866 : ℚ) / 10^19 < (594603557501360533358749 : ℚ) / 10^24 * Sq 64 := by
  decide +kernel
theorem num_ST : (594603557501360533358751 : ℚ) / 10^24 * (Sq 64 + Tq 64) < (8472130847939790867 : ℚ) / 10^19 := by
  decide +kernel
theorem num_pos : 0 < Sq 64 ∧ 0 < Sq 64 + Tq 64 := by
  decide +kernel

lemma J_bounds : (8472130847939790866 : ℝ) / 10^19 < ∫ u in (0:ℝ)..(1/2), gb u ∧
    ∫ u in (0:ℝ)..(1/2), gb u < (8472130847939790867 : ℝ) / 10^19 := by
  have hlo := J_lo 63
  have hhi := J_hi 63
  rw [lo_cast 64] at hlo
  rw [lo_cast 64] at hhi
  have h63 : ((63 + 1 : ℕ) : ℝ) = (64 : ℝ) := by norm_num
  have hh := hi_cast 64
  have e1 : (4:ℝ) * cc (63+1) * (ss * (1/2) ^ (63+1) / (((63+1:ℕ):ℝ) + 3/4)) = ss * (Tq 64 : ℝ) := hh
  rw [e1] at hhi
  obtain ⟨hSp, hSTp⟩ := num_pos
  have hSp' : (0:ℝ) < (Sq 64 : ℝ) := by exact_mod_cast hSp
  have hSTp' : (0:ℝ) < (Sq 64 : ℝ) + (Tq 64 : ℝ) := by exact_mod_cast hSTp
  have h1 : ((8472130847939790866 : ℚ) / 10^19 : ℝ) < (((594603557501360533358749 : ℚ) / 10^24 * Sq 64 : ℚ) : ℝ) := by
    exact_mod_cast num_Sq
  have h2 : (((594603557501360533358751 : ℚ) / 10^24 * (Sq 64 + Tq 64) : ℚ) : ℝ) < ((8472130847939790867 : ℚ) / 10^19 : ℝ) := by
    exact_mod_cast num_ST
  push_cast at h1 h2
  have s1 := ss_lo
  have s2 := ss_hi
  have s1' : ((594603557501360533358749 : ℚ) / 10^24 : ℝ) < ss := by exact_mod_cast s1
  have s2' : ss < ((594603557501360533358751 : ℚ) / 10^24 : ℝ) := by exact_mod_cast s2
  push_cast at s1' s2'
  constructor
  · nlinarith
  · nlinarith
end U01Aux

namespace U01Aux
lemma I_closed : U01.I = 1 - 2 * Real.pi ^ (3 / 2 : ℝ) / Real.Gamma (1 / 4) ^ 2 := by
  rw [I_eq, sub_eq, half_eq, gb_total, gamma34]
  have h32 : Real.pi ^ (3/2 : ℝ) = Real.pi * Real.sqrt Real.pi := by
    rw [show (3/2:ℝ) = 1 + 1/2 by norm_num, Real.rpow_add Real.pi_pos, Real.rpow_one,
      ← Real.sqrt_eq_rpow]
  rw [h32]
  have hp : 0 < Real.sqrt Real.pi := Real.sqrt_pos.2 Real.pi_pos
  have hsq : Real.sqrt Real.pi ^ 2 = Real.pi := Real.sq_sqrt Real.pi_pos.le
  have hg : 0 < Real.Gamma (1/4) := Real.Gamma_pos_of_pos (by norm_num)
  field_simp
  nlinarith [hsq]
end U01Aux

theorem U01.proof_1 : U01.stmt_1 := U01Aux.f_cont.continuousOn

theorem U01.proof_2 : U01.stmt_2 :=
  ⟨U01Aux.f_cont.intervalIntegrable _ _, U01Aux.I_closed⟩

theorem U01.proof_3 : U01.stmt_3 := by
  unfold U01.stmt_3 U01.R
  rw [U01Aux.I_closed]; ring

namespace U01Aux
lemma X_eq : 2 * Real.pi ^ (3 / 2 : ℝ) / Real.Gamma (1 / 4) ^ 2 = ∫ u in (0:ℝ)..(1/2), gb u := by
  have h1 := I_eq
  rw [sub_eq, I_closed] at h1
  linarith
end U01Aux

theorem U01.proof_4 : U01.stmt_4 := by
  have h := U01Aux.J_bounds
  rw [← U01Aux.X_eq] at h
  unfold U01.stmt_4
  constructor <;> norm_num at h ⊢ <;> linarith [h.1, h.2]
