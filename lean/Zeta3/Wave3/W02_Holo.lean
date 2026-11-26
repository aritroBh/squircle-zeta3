import Mathlib
/-!
# W02 -- the perimeter integral, the function `F`, and its holomorphy (Lemma `lem:holo`)

Target id: W02.
Paper labels: `eq:Pint`, `eq:F`, `lem:holo`, `eq:expansion`, `eq:M1t`
(paper/ieee/main.tex, lines 820-860). The proofs supplement has no separate proof of
`lem:holo`; the intermediate statements follow the short proof printed in main.tex.

Quoted verbatim (main.tex, eq:Pint):

  P(p)=\frac4p\int_0^1\sqrt{t^{2/p-2}+(1-t)^{2/p-2}}\;dt .

Quoted verbatim (main.tex, before and in eq:F):

  To expand near the square we put $\varepsilon=1/p$, $a_\varepsilon=t^{\varepsilon-1}$ and
  $b_\varepsilon=(1-t)^{\varepsilon-1}$. Since $\int_0^1\varepsilon(a_\varepsilon+b_\varepsilon)\dd t=2$,
  \eqref{eq:Pint} can be rewritten as
  \frac{P(p)}4=2-\varepsilon F(\varepsilon),\qquad
  F(\varepsilon)=\int_0^1\Bigl(a_\varepsilon+b_\varepsilon-\sqrt{a_\varepsilon^2+b_\varepsilon^2}\Bigr)dt .

Quoted verbatim (main.tex, Lemma lem:holo with eq:expansion):

  $F$ extends to a holomorphic function on $\{\operatorname{Re}\varepsilon<1\}$. Hence
  $P(p)/8=1+\varepsilon M(\varepsilon)$ with $M=-F/2$ holomorphic in $|\varepsilon|<1$, the expansion
  P(p)=8\Bigl(1+\frac{M_0}{p}+\frac{M_1}{p^2}+\frac{M_2}{p^3}+\frac{M_3}{p^4}+\cdots\Bigr)
  converges for every $p>1$, and $M_k=-F^{(k)}(0)/(2\,k!)$.

Quoted verbatim (main.tex, proof of lem:holo):

  By the symmetry $t\leftrightarrow1-t$, $F=2\int_0^{1/2}$. For $0<t\le\frac12$ put
  $z=(b_\varepsilon/a_\varepsilon)^2=(t/(1-t))^{2-2\varepsilon}$. For $\operatorname{Re}\varepsilon<1$
  we have $|z|\le1$, with equality only at $t=\frac12$, so
  $\sqrt{a_\varepsilon^2+b_\varepsilon^2}:=a_\varepsilon\sqrt{1+z}$, with the principal branch of the
  root ($1+z\notin(-\infty,0]$), is holomorphic in $\varepsilon$ and is the positive root for real
  $\varepsilon$. Writing the integrand as $b_\varepsilon-b_\varepsilon^2/(a_\varepsilon(1+\sqrt{1+z}))$
  and using $|b_\varepsilon|\le2^{1-\operatorname{Re}\varepsilon}\le|a_\varepsilon|$ bounds it by
  $2^{2-\operatorname{Re}\varepsilon}$, uniformly on compact sets, so we may differentiate under the
  integral sign.

Quoted verbatim (main.tex, eq:M1t):

  M_1=\int_0^1\frac{\log t}{t}\Bigl(\frac{1-t}{q}-1\Bigr)dt,\qquad q=\sqrt{t^2+(1-t)^2},

## Conventions

* Real side: `aR ε t = t^(ε-1)`, `bR ε t = (1-t)^(ε-1)` (real `rpow`; only `t ∈ (0,1)` matters,
  the endpoints are a null set), `FR ε = ∫_0^1 (aR + bR - √(aR² + bR²))`, and
  `Pint p = (4/p) ∫_0^1 √(t^(2/p-2) + (1-t)^(2/p-2))` exactly as in eq:Pint.
  Every interval integral whose value is asserted is accompanied by an explicit
  `IntervalIntegrable` claim, so no statement holds through Lean's "non-integrable ↦ 0" junk value.
* Complex side (the proof of lem:holo): for `ε : ℂ` and real `t ∈ (0, 1/2]`,
  `aC = t^(ε-1)`, `bC = (1-t)^(ε-1)` (complex `cpow` of a positive real base, i.e.
  `exp((ε-1) log t)`), `zC = (t/(1-t))^(2-2ε)`, and the integrand is
  `gC = aC + bC - aC · (1+zC)^(1/2)` with the principal branch `w^(1/2) = exp(log w / 2)`.
  The holomorphic extension is `FC ε = 2 ∫_0^{1/2} gC ε t dt` (the paper's `F = 2∫_0^{1/2}`).
* `M_k := -F^{(k)}(0)/(2 k!)` is `-(iteratedDeriv k FC 0) / (2 k!)` (complex derivatives of the
  extension `FC`). Indexing checked against eq:expansion: `M(ε) = Σ_{k≥0} M_k ε^k`, so
  `P/8 = 1 + ε M(ε) = 1 + Σ_{k≥0} M_k p^{-(k+1)}`, i.e. `M_0` goes with `1/p`, `M_1` with `1/p²`.

## Statements

* stmt_1: for real `0 < ε < 1`, `ε(aR + bR)` is integrable on `[0,1]` and `∫_0^1 ε(aR+bR) = 2`.
* stmt_2: for real `p > 1`, `ε = 1/p`: both integrands are integrable on `[0,1]` and
  `P(p)/4 = 2 - ε F(ε)` (eq:F).
* stmt_3: the pointwise bound of the proof: `‖gC ε t‖ ≤ 2^(2 - Re ε)` for `0 < t ≤ 1/2`,
  `Re ε < 1`.
* stmt_4: for real `ε < 1`, the real integrand of `F` is integrable on `[0,1]` and the complex
  extension agrees with it: `FC ε = FR ε`.
* stmt_5: lem:holo, first sentence: `gC ε` is integrable on `[0,1/2]` for `Re ε < 1`, and `FC` is
  holomorphic (`DifferentiableOn ℂ`) on `{Re ε < 1}`.
* stmt_6: eq:expansion: for every real `p > 1`, `Σ_{k≥0} M_k / p^(k+1)` converges (`HasSum`) to
  `P(p)/8 - 1`, with `M_k = -F^{(k)}(0)/(2 k!)`.
* stmt_7: eq:M1t: `(log t / t)((1-t)/q - 1)` is integrable on `[0,1]` and
  `M_1 = ∫_0^1 (log t / t)((1-t)/q - 1) dt`, `q = √(t² + (1-t)²)`.

## NOT COVERED

* "with equality only at t = 1/2" (`|z| = 1` iff `t = 1/2`) and "is the positive root for real ε"
  are not stated separately; the latter is contained in stmt_4 (agreement with the real `F`).
* "uniformly on compact sets, so we may differentiate under the integral sign" is the proof
  method of stmt_5, not a separate claim.
* "$M=-F/2$ holomorphic in $|\varepsilon|<1$" is not stated separately: it follows from stmt_5
  since the unit disk lies in `{Re ε < 1}`.
* The identification `κ = 8 M_1` with Ahmed's constant (a citation) is not formalized.
-/
-- BEGIN STATEMENTS
namespace W02

open MeasureTheory Real

/-- `a_ε(t) = t^(ε-1)` (real). -/
noncomputable def aR (ε t : ℝ) : ℝ := t ^ (ε - 1)

/-- `b_ε(t) = (1-t)^(ε-1)` (real). -/
noncomputable def bR (ε t : ℝ) : ℝ := (1 - t) ^ (ε - 1)

/-- The integrand of `F` for real `ε`: `a + b - √(a² + b²)`. -/
noncomputable def gR (ε t : ℝ) : ℝ := aR ε t + bR ε t - Real.sqrt (aR ε t ^ 2 + bR ε t ^ 2)

/-- eq:F, real `ε`: `F(ε) = ∫_0^1 (a + b - √(a² + b²)) dt`. -/
noncomputable def FR (ε : ℝ) : ℝ := ∫ t in (0 : ℝ)..1, gR ε t

/-- The integrand of eq:Pint: `√(t^(2/p-2) + (1-t)^(2/p-2))`. -/
noncomputable def pIntegrand (p t : ℝ) : ℝ :=
  Real.sqrt (t ^ (2 / p - 2) + (1 - t) ^ (2 / p - 2))

/-- eq:Pint: `P(p) = (4/p) ∫_0^1 √(t^(2/p-2) + (1-t)^(2/p-2)) dt`. -/
noncomputable def Pint (p : ℝ) : ℝ := 4 / p * ∫ t in (0 : ℝ)..1, pIntegrand p t

/-- `a_ε(t) = t^(ε-1)` for complex `ε` (cpow of the positive real `t`). -/
noncomputable def aC (ε : ℂ) (t : ℝ) : ℂ := (t : ℂ) ^ (ε - 1)

/-- `b_ε(t) = (1-t)^(ε-1)` for complex `ε`. -/
noncomputable def bC (ε : ℂ) (t : ℝ) : ℂ := ((1 - t : ℝ) : ℂ) ^ (ε - 1)

/-- `z = (t/(1-t))^(2-2ε)`. -/
noncomputable def zC (ε : ℂ) (t : ℝ) : ℂ := ((t / (1 - t) : ℝ) : ℂ) ^ (2 - 2 * ε)

/-- The complex integrand on `(0, 1/2]`: `a + b - a √(1+z)`, principal branch of the root. -/
noncomputable def gC (ε : ℂ) (t : ℝ) : ℂ :=
  aC ε t + bC ε t - aC ε t * (1 + zC ε t) ^ (1 / 2 : ℂ)

/-- The extension of `F` from the proof of lem:holo: `F(ε) = 2 ∫_0^{1/2} gC ε t dt`. -/
noncomputable def FC (ε : ℂ) : ℂ := 2 * ∫ t in (0 : ℝ)..(1 / 2), gC ε t

/-- `M_k = -F^{(k)}(0) / (2 k!)`. -/
noncomputable def Mk (k : ℕ) : ℂ := -(iteratedDeriv k FC 0) / (2 * (Nat.factorial k : ℂ))

/-- `q = √(t² + (1-t)²)` of eq:M1t. -/
noncomputable def q (t : ℝ) : ℝ := Real.sqrt (t ^ 2 + (1 - t) ^ 2)

/-- Text before eq:F: "$\int_0^1\varepsilon(a_\varepsilon+b_\varepsilon)\dd t=2$", for real
`0 < ε < 1` (that is, `p > 1`), with the integrability. -/
def stmt_1 : Prop :=
  ∀ ε : ℝ, 0 < ε → ε < 1 →
    IntervalIntegrable (fun t => ε * (aR ε t + bR ε t)) volume 0 1 ∧
      ∫ t in (0 : ℝ)..1, ε * (aR ε t + bR ε t) = 2

/-- eq:F: "$\frac{P(p)}4=2-\varepsilon F(\varepsilon)$" for real `p > 1`, `ε = 1/p`, with the
integrability of both integrands. -/
def stmt_2 : Prop :=
  ∀ p : ℝ, 1 < p →
    IntervalIntegrable (pIntegrand p) volume 0 1 ∧
      IntervalIntegrable (gR (1 / p)) volume 0 1 ∧
      Pint p / 4 = 2 - (1 / p) * FR (1 / p)

/-- Proof of lem:holo: for `0 < t ≤ 1/2` and `Re ε < 1` the integrand (with the branch
`a √(1+z)`) is bounded by `2^(2 - Re ε)`. -/
def stmt_3 : Prop :=
  ∀ ε : ℂ, ε.re < 1 → ∀ t : ℝ, 0 < t → t ≤ 1 / 2 →
    ‖gC ε t‖ ≤ (2 : ℝ) ^ (2 - ε.re)

/-- Proof of lem:holo: the branch "is the positive root for real ε", and "$F=2\int_0^{1/2}$":
for real `ε < 1` the real integrand is integrable on `[0,1]` and the extension `FC` agrees with
the real `F`. -/
def stmt_4 : Prop :=
  ∀ ε : ℝ, ε < 1 →
    IntervalIntegrable (gR ε) volume 0 1 ∧ FC (ε : ℂ) = ((FR ε : ℝ) : ℂ)

/-- lem:holo: "$F$ extends to a holomorphic function on $\{\operatorname{Re}\varepsilon<1\}$",
with the integrability of the complex integrand on `[0, 1/2]`. -/
def stmt_5 : Prop :=
  (∀ ε : ℂ, ε.re < 1 → IntervalIntegrable (gC ε) volume 0 (1 / 2)) ∧
    DifferentiableOn ℂ FC {ε : ℂ | ε.re < 1}

/-- lem:holo, eq:expansion: "$P(p)=8(1+\frac{M_0}{p}+\frac{M_1}{p^2}+\cdots)$ converges for every
$p>1$", with `M_k = -F^{(k)}(0)/(2 k!)`. Stated as `Σ_{k≥0} M_k / p^(k+1) = P(p)/8 - 1`. -/
def stmt_6 : Prop :=
  ∀ p : ℝ, 1 < p →
    HasSum (fun k : ℕ => Mk k / (p : ℂ) ^ (k + 1)) ((Pint p / 8 - 1 : ℝ) : ℂ)

/-- eq:M1t: "$M_1=\int_0^1\frac{\log t}{t}\Bigl(\frac{1-t}{q}-1\Bigr)dt$", with the integrability. -/
def stmt_7 : Prop :=
  IntervalIntegrable (fun t => Real.log t / t * ((1 - t) / q t - 1)) volume 0 1 ∧
    Mk 1 = (((∫ t in (0 : ℝ)..1, Real.log t / t * ((1 - t) / q t - 1)) : ℝ) : ℂ)

end W02
-- END STATEMENTS

open MeasureTheory Real Set

namespace W02

lemma rpow_iInt {ε : ℝ} (h0 : 0 < ε) :
    IntervalIntegrable (fun t : ℝ => t ^ (ε - 1)) volume 0 1 :=
  intervalIntegral.intervalIntegrable_rpow' (by linarith)

lemma rpow_iInt' {ε : ℝ} (h0 : 0 < ε) :
    IntervalIntegrable (fun t : ℝ => (1 - t) ^ (ε - 1)) volume 0 1 := by
  have := (rpow_iInt h0).comp_sub_left 1
  simpa using this.symm

lemma rpow_int {ε : ℝ} (h0 : 0 < ε) : ∫ t in (0:ℝ)..1, t ^ (ε - 1) = 1 / ε := by
  rw [integral_rpow (Or.inl (by linarith))]
  simp [Real.zero_rpow h0.ne']

lemma rpow_int' {ε : ℝ} (h0 : 0 < ε) : ∫ t in (0:ℝ)..1, (1 - t) ^ (ε - 1) = 1 / ε := by
  rw [intervalIntegral.integral_comp_sub_left (fun t : ℝ => t ^ (ε - 1)) 1]
  simpa using rpow_int h0

end W02

theorem W02.proof_1 : W02.stmt_1 := by
  intro ε h0 h1
  have ha := W02.rpow_iInt h0
  have hb := W02.rpow_iInt' h0
  refine ⟨((ha.add hb).const_mul ε), ?_⟩
  simp only [W02.aR, W02.bR]
  rw [intervalIntegral.integral_const_mul, intervalIntegral.integral_add ha hb,
    W02.rpow_int h0, W02.rpow_int' h0]
  field_simp
  ring

namespace W02

lemma half_rpow (ε : ℝ) : (1 / 2 : ℝ) ^ (ε - 1) = (2 : ℝ) ^ (1 - ε) := by
  rw [one_div, Real.inv_rpow (by norm_num), ← Real.rpow_neg (by norm_num)]
  ring_nf

lemma gR_bound {ε t : ℝ} (hε : ε < 1) (ht0 : 0 < t) (ht1 : t < 1) :
    0 ≤ gR ε t ∧ gR ε t ≤ (2 : ℝ) ^ (1 - ε) := by
  have ha : 0 < aR ε t := Real.rpow_pos_of_pos ht0 _
  have hb : 0 < bR ε t := Real.rpow_pos_of_pos (by linarith) _
  have hs1 : Real.sqrt (aR ε t ^ 2 + bR ε t ^ 2) ≤ aR ε t + bR ε t := by
    rw [Real.sqrt_le_iff]; constructor <;> nlinarith
  have hs2 : aR ε t ≤ Real.sqrt (aR ε t ^ 2 + bR ε t ^ 2) := by
    rw [Real.le_sqrt ha.le (by positivity)]; nlinarith
  have hs3 : bR ε t ≤ Real.sqrt (aR ε t ^ 2 + bR ε t ^ 2) := by
    rw [Real.le_sqrt hb.le (by positivity)]; nlinarith
  refine ⟨by unfold gR; linarith, ?_⟩
  rcases le_total t (1 / 2) with h | h
  · have : bR ε t ≤ (1 / 2 : ℝ) ^ (ε - 1) :=
      Real.rpow_le_rpow_of_nonpos (by norm_num) (by linarith) (by linarith)
    rw [half_rpow] at this
    unfold gR; linarith
  · have : aR ε t ≤ (1 / 2 : ℝ) ^ (ε - 1) :=
      Real.rpow_le_rpow_of_nonpos (by norm_num) h (by linarith)
    rw [half_rpow] at this
    unfold gR; linarith

lemma gR_meas (ε : ℝ) : Measurable (gR ε) := by
  unfold gR aR bR; fun_prop

lemma gR_iInt {ε : ℝ} (hε : ε < 1) : IntervalIntegrable (gR ε) volume 0 1 := by
  refine IntervalIntegrable.mono_fun' (g := fun _ => (2 : ℝ) ^ (1 - ε))
    intervalIntegrable_const (gR_meas ε).aestronglyMeasurable ?_
  rw [uIoc_of_le zero_le_one, ← restrict_Ioo_eq_restrict_Ioc]
  refine ae_restrict_of_forall_mem measurableSet_Ioo fun t ht => ?_
  have := gR_bound hε ht.1 ht.2
  simp only [Real.norm_eq_abs, abs_of_nonneg this.1]
  exact this.2

lemma gR_symm (ε t : ℝ) : gR ε (1 - t) = gR ε t := by
  unfold gR aR bR
  simp only [sub_sub_cancel]
  ring_nf

lemma FR_half {ε : ℝ} (hε : ε < 1) : FR ε = 2 * ∫ t in (0:ℝ)..(1/2), gR ε t := by
  have hi := gR_iInt hε
  have h1 : IntervalIntegrable (gR ε) volume 0 (1/2) :=
    hi.mono_set (by rw [uIcc_of_le (by norm_num), uIcc_of_le (by norm_num)]
                    exact Icc_subset_Icc le_rfl (by norm_num))
  have h2 : IntervalIntegrable (gR ε) volume (1/2) 1 :=
    hi.mono_set (by rw [uIcc_of_le (by norm_num), uIcc_of_le (by norm_num)]
                    exact Icc_subset_Icc (by norm_num) le_rfl)
  unfold FR
  rw [← intervalIntegral.integral_add_adjacent_intervals h1 h2]
  have : ∫ t in (1/2:ℝ)..1, gR ε t = ∫ t in (0:ℝ)..(1/2), gR ε t := by
    have := intervalIntegral.integral_comp_sub_left (gR ε) (a := 0) (b := 1/2) 1
    simp only [gR_symm] at this
    rw [this]; norm_num
  rw [this]; ring

end W02

theorem W02.proof_2 : W02.stmt_2 := by
  intro p hp
  set ε := 1 / p with hεdef
  have h0 : 0 < ε := by positivity
  have h1 : ε < 1 := by rw [hεdef, div_lt_one (by linarith)]; exact hp
  have ha := W02.rpow_iInt h0
  have hb := W02.rpow_iInt' h0
  have hab : IntervalIntegrable (fun t => W02.aR ε t + W02.bR ε t) volume 0 1 := ha.add hb
  have hg := W02.gR_iInt h1
  have heq : EqOn (fun t => W02.aR ε t + W02.bR ε t - W02.gR ε t) (W02.pIntegrand p)
      (uIoo 0 1) := by
    intro t ht
    rw [uIoo_of_le zero_le_one] at ht
    have hsq : ∀ x : ℝ, 0 < x → x ^ (2 / p - 2) = (x ^ (ε - 1)) ^ 2 := by
      intro x hx
      rw [sq, ← Real.rpow_add hx]; congr 1; rw [hεdef]; ring
    simp only [W02.gR, W02.pIntegrand, W02.aR, W02.bR]
    rw [hsq t ht.1, hsq (1 - t) (by linarith [ht.2])]
    ring
  have hP : IntervalIntegrable (W02.pIntegrand p) volume 0 1 :=
    (intervalIntegrable_congr_uIoo heq).1 (hab.sub hg)
  refine ⟨hP, hg, ?_⟩
  have h2 := (W02.proof_1 ε h0 h1).2
  rw [intervalIntegral.integral_const_mul] at h2
  unfold W02.Pint W02.FR
  rw [← intervalIntegral.integral_congr_uIoo heq, intervalIntegral.integral_sub hab hg]
  have : (4 / p) = 4 * ε := by rw [hεdef]; ring
  rw [this]
  have hpe : 1 / p = ε := hεdef.symm
  linear_combination h2 + (∫ t in (0:ℝ)..1, W02.gR ε t) * hpe

namespace W02

lemma cpow_pos {x : ℝ} (hx : 0 < x) (w : ℂ) :
    ((x : ℝ) : ℂ) ^ w = Complex.exp ((Real.log x : ℂ) * w) := by
  rw [Complex.cpow_def_of_ne_zero (by exact_mod_cast hx.ne'), Complex.ofReal_log hx.le]

lemma aC_ne (ε : ℂ) {t : ℝ} (ht : 0 < t) : aC ε t ≠ 0 := by
  unfold aC; rw [cpow_pos ht]; exact Complex.exp_ne_zero _

lemma azC {ε : ℂ} {t : ℝ} (ht0 : 0 < t) (ht1 : t < 1) :
    aC ε t * zC ε t * aC ε t = bC ε t ^ 2 := by
  have h1 : 0 < 1 - t := by linarith
  unfold aC bC zC
  rw [cpow_pos ht0, cpow_pos h1, cpow_pos (div_pos ht0 h1),
    Real.log_div ht0.ne' h1.ne', sq, ← Complex.exp_add, ← Complex.exp_add, ← Complex.exp_add]
  congr 1
  push_cast
  ring

lemma s_sq (w : ℂ) : (w ^ (1 / 2 : ℂ)) ^ 2 = w := by
  rw [one_div]; exact Complex.cpow_ofNat_inv_pow w 2

lemma s_re (w : ℂ) : 0 ≤ (w ^ (1 / 2 : ℂ)).re := by
  rw [one_div, Complex.cpow_inv_two_re]; positivity

lemma one_add_s_ne (w : ℂ) : 1 + w ^ (1 / 2 : ℂ) ≠ 0 := by
  intro h
  have := congrArg Complex.re h
  simp only [Complex.add_re, Complex.one_re, Complex.zero_re] at this
  linarith [s_re w]

lemma one_le_norm_one_add_s (w : ℂ) : 1 ≤ ‖1 + w ^ (1 / 2 : ℂ)‖ := by
  have h := Complex.re_le_norm (1 + w ^ (1 / 2 : ℂ))
  simp only [Complex.add_re, Complex.one_re] at h
  linarith [s_re w]

lemma gC_alt {ε : ℂ} {t : ℝ} :
    gC ε t = bC ε t - aC ε t * zC ε t / (1 + (1 + zC ε t) ^ (1 / 2 : ℂ)) := by
  have h1s := one_add_s_ne (1 + zC ε t)
  have hs2 := s_sq (1 + zC ε t)
  have : aC ε t * (1 - (1 + zC ε t) ^ (1 / 2 : ℂ)) =
      -(aC ε t * zC ε t) / (1 + (1 + zC ε t) ^ (1 / 2 : ℂ)) := by
    rw [eq_div_iff h1s]; linear_combination (-aC ε t) * hs2
  unfold gC
  rw [neg_div] at this
  linear_combination this

lemma gC_bound {ε : ℂ} (hε : ε.re < 1) {t : ℝ} (ht0 : 0 < t) (ht : t ≤ 1 / 2) :
    ‖gC ε t‖ ≤ 2 * ‖bC ε t‖ ∧ ‖bC ε t‖ ≤ (2 : ℝ) ^ (1 - ε.re) := by
  have ht1 : t < 1 := by linarith
  have h1 : 0 < 1 - t := by linarith
  have hna : ‖aC ε t‖ = t ^ (ε.re - 1) := by
    unfold aC; rw [Complex.norm_cpow_eq_rpow_re_of_pos ht0]; simp
  have hnb : ‖bC ε t‖ = (1 - t) ^ (ε.re - 1) := by
    unfold bC; rw [Complex.norm_cpow_eq_rpow_re_of_pos h1]; simp
  have hba : ‖bC ε t‖ ≤ ‖aC ε t‖ := by
    rw [hna, hnb]; exact Real.rpow_le_rpow_of_nonpos ht0 (by linarith) (by linarith)
  have hb2 : ‖bC ε t‖ ≤ (2 : ℝ) ^ (1 - ε.re) := by
    rw [hnb, ← half_rpow]
    exact Real.rpow_le_rpow_of_nonpos (by norm_num) (by linarith) (by linarith)
  have hapos : 0 < ‖aC ε t‖ := norm_pos_iff.2 (aC_ne ε ht0)
  have haz : ‖aC ε t * zC ε t‖ * ‖aC ε t‖ = ‖bC ε t‖ ^ 2 := by
    rw [← norm_mul, azC ht0 ht1, norm_pow]
  have haz2 : ‖aC ε t * zC ε t‖ ≤ ‖bC ε t‖ := by
    by_contra hc
    replace hc := not_le.1 hc
    have : ‖bC ε t‖ * ‖aC ε t‖ < ‖aC ε t * zC ε t‖ * ‖aC ε t‖ :=
      mul_lt_mul_of_pos_right hc hapos
    nlinarith [norm_nonneg (bC ε t)]
  refine ⟨?_, hb2⟩
  rw [gC_alt]
  calc ‖bC ε t - aC ε t * zC ε t / (1 + (1 + zC ε t) ^ (1 / 2 : ℂ))‖
      ≤ ‖bC ε t‖ + ‖aC ε t * zC ε t / (1 + (1 + zC ε t) ^ (1 / 2 : ℂ))‖ := norm_sub_le _ _
    _ ≤ ‖bC ε t‖ + ‖aC ε t * zC ε t‖ := by
        gcongr
        rw [norm_div]
        exact div_le_self (norm_nonneg _) (one_le_norm_one_add_s _)
    _ ≤ 2 * ‖bC ε t‖ := by linarith

end W02

theorem W02.proof_3 : W02.stmt_3 := by
  intro ε hε t ht0 ht
  obtain ⟨h1, h2⟩ := W02.gC_bound hε ht0 ht
  have : (2 : ℝ) ^ (2 - ε.re) = 2 * (2 : ℝ) ^ (1 - ε.re) := by
    rw [show (2 : ℝ) - ε.re = 1 + (1 - ε.re) by ring, Real.rpow_add (by norm_num), Real.rpow_one]
  rw [this]; linarith

namespace W02

lemma aC_ofReal (ε : ℝ) {t : ℝ} (ht : 0 ≤ t) : aC (ε : ℂ) t = ((aR ε t : ℝ) : ℂ) := by
  unfold aC aR; rw [Complex.ofReal_cpow ht]; push_cast; rfl

lemma bC_ofReal (ε : ℝ) {t : ℝ} (ht : t ≤ 1) : bC (ε : ℂ) t = ((bR ε t : ℝ) : ℂ) := by
  unfold bC bR; rw [Complex.ofReal_cpow (by linarith)]; push_cast; rfl

lemma zC_ofReal (ε : ℝ) {t : ℝ} (ht0 : 0 ≤ t) (ht : t < 1) :
    zC (ε : ℂ) t = (((t / (1 - t)) ^ (2 - 2 * ε) : ℝ) : ℂ) := by
  unfold zC
  rw [Complex.ofReal_cpow (div_nonneg ht0 (by linarith))]; push_cast; rfl

lemma gC_ofReal (ε : ℝ) {t : ℝ} (ht0 : 0 < t) (ht : t < 1) :
    gC (ε : ℂ) t = ((gR ε t : ℝ) : ℂ) := by
  have hz0 : 0 ≤ (t / (1 - t)) ^ (2 - 2 * ε) :=
    Real.rpow_nonneg (div_nonneg ht0.le (by linarith)) _
  have haz := azC (ε := (ε : ℂ)) ht0 ht
  rw [aC_ofReal ε ht0.le, bC_ofReal ε ht.le, zC_ofReal ε ht0.le ht] at haz
  have haz' : aR ε t * (t / (1 - t)) ^ (2 - 2 * ε) * aR ε t = bR ε t ^ 2 := by
    exact_mod_cast haz
  have hs' : ∀ z : ℝ, 0 ≤ z → (1 + (z : ℂ)) ^ (1 / 2 : ℂ) = ((Real.sqrt (1 + z) : ℝ) : ℂ) := by
    intro z hz
    rw [Real.sqrt_eq_rpow, Complex.ofReal_cpow (by linarith)]
    push_cast; rfl
  have hs := hs' _ hz0
  have ha0 : 0 ≤ aR ε t := Real.rpow_nonneg ht0.le _
  have hsq : aR ε t * Real.sqrt (1 + (t / (1 - t)) ^ (2 - 2 * ε)) =
      Real.sqrt (aR ε t ^ 2 + bR ε t ^ 2) := by
    rw [← Real.sqrt_sq ha0, ← Real.sqrt_mul (sq_nonneg _), Real.sqrt_sq ha0]
    congr 1
    rw [← haz']; ring
  unfold gC gR
  rw [aC_ofReal ε ht0.le, bC_ofReal ε ht.le, zC_ofReal ε ht0.le ht, hs, ← hsq]
  push_cast; ring

end W02

theorem W02.proof_4 : W02.stmt_4 := by
  intro ε hε
  refine ⟨W02.gR_iInt hε, ?_⟩
  unfold W02.FC
  rw [W02.FR_half hε]
  have : ∫ t in (0:ℝ)..(1/2), W02.gC (ε : ℂ) t = ∫ t in (0:ℝ)..(1/2), ((W02.gR ε t : ℝ) : ℂ) := by
    refine intervalIntegral.integral_congr_ae (Filter.Eventually.of_forall fun t ht => ?_)
    rw [uIoc_of_le (by norm_num)] at ht
    exact W02.gC_ofReal ε ht.1 (by linarith [ht.2])
  rw [this, intervalIntegral.integral_ofReal]
  push_cast; ring

namespace W02

lemma gC_meas (ε : ℂ) : Measurable (gC ε) := by
  unfold gC aC bC zC
  fun_prop

/-- explicit `ε`-derivative of `gC`. -/
noncomputable def gCd (e : ℂ) (t : ℝ) : ℂ :=
  aC e t * Complex.log (t : ℂ) + bC e t * Complex.log ((1 - t : ℝ) : ℂ) -
    (aC e t * Complex.log (t : ℂ) * (1 + zC e t) ^ (1 / 2 : ℂ) +
      aC e t * ((1 / 2 : ℂ) * (1 + zC e t) ^ ((1 / 2 : ℂ) - 1) *
        (zC e t * Complex.log ((t / (1 - t) : ℝ) : ℂ) * (-2))))

lemma gCd_meas (ε : ℂ) : Measurable (gCd ε) := by
  unfold gCd aC bC zC
  have := Complex.measurable_log
  fun_prop

lemma zC_slit {e : ℂ} (he : e.re < 1) {t : ℝ} (ht0 : 0 < t) (ht : t ≤ 1 / 2) :
    1 + zC e t ∈ Complex.slitPlane := by
  have h1 : 0 < 1 - t := by linarith
  rcases eq_or_lt_of_le ht with h | h
  · subst h
    unfold zC
    norm_num
  · have hu0 : 0 < t / (1 - t) := div_pos ht0 h1
    have hu1 : t / (1 - t) < 1 := by rw [div_lt_one h1]; linarith
    have hn : ‖zC e t‖ < 1 := by
      unfold zC
      rw [Complex.norm_cpow_eq_rpow_re_of_pos hu0]
      apply Real.rpow_lt_one hu0.le hu1
      simp; linarith
    refine Complex.mem_slitPlane_iff.2 (Or.inl ?_)
    have := Complex.abs_re_le_norm (zC e t)
    simp only [Complex.add_re, Complex.one_re]
    linarith [neg_abs_le (zC e t).re]

lemma gC_hasDerivAt {e : ℂ} (he : e.re < 1) {t : ℝ} (ht0 : 0 < t) (ht : t ≤ 1 / 2) :
    HasDerivAt (fun x => gC x t) (gCd e t) e := by
  have h1 : 0 < 1 - t := by linarith
  have ha : HasDerivAt (fun x => aC x t) (aC e t * Complex.log (t : ℂ) * 1) e :=
    ((hasDerivAt_id e).sub_const 1).const_cpow (Or.inl (by exact_mod_cast ht0.ne'))
  have hb : HasDerivAt (fun x => bC x t) (bC e t * Complex.log ((1 - t : ℝ) : ℂ) * 1) e :=
    ((hasDerivAt_id e).sub_const 1).const_cpow (Or.inl (by exact_mod_cast h1.ne'))
  have hz : HasDerivAt (fun x => zC x t)
      (zC e t * Complex.log ((t / (1 - t) : ℝ) : ℂ) * (0 - 2 * 1)) e :=
    ((hasDerivAt_const e (2 : ℂ)).sub ((hasDerivAt_id e).const_mul 2)).const_cpow
      (Or.inl (by exact_mod_cast (div_pos ht0 h1).ne'))
  have hS := ((hasDerivAt_const e (1 : ℂ)).add hz).cpow_const (c := (1 / 2 : ℂ))
    (zC_slit he ht0 ht)
  have := (ha.add hb).sub (ha.mul hS)
  unfold gC gCd
  convert this using 1
  simp only [Pi.add_apply]
  ring

lemma gC_iInt {ε : ℂ} (hε : ε.re < 1) : IntervalIntegrable (gC ε) volume 0 (1 / 2) := by
  refine IntervalIntegrable.mono_fun' (g := fun _ => (2 : ℝ) ^ (2 - ε.re))
    intervalIntegrable_const (gC_meas ε).aestronglyMeasurable ?_
  rw [uIoc_of_le (by norm_num)]
  exact ae_restrict_of_forall_mem measurableSet_Ioc fun t ht => proof_3 ε hε t ht.1 ht.2

lemma gCd_bound {x₀ : ℂ} (h : x₀.re < 1) {x : ℂ}
    (hx : x ∈ Metric.ball x₀ ((1 - x₀.re) / 3)) {t : ℝ} (ht0 : 0 < t) (ht : t ≤ 1 / 2) :
    ‖gCd x t‖ ≤ (2 : ℝ) ^ (2 - x₀.re + 2 * ((1 - x₀.re) / 3)) / ((1 - x₀.re) / 3) := by
  set δ := (1 - x₀.re) / 3 with hδ
  have hδ0 : 0 < δ := by rw [hδ]; linarith
  have hxre : |x.re - x₀.re| < δ := by
    have := Complex.abs_re_le_norm (x - x₀)
    rw [Metric.mem_ball, Complex.dist_eq] at hx
    simp only [Complex.sub_re] at this
    linarith
  have hsub : Metric.closedBall x δ ⊆ {e : ℂ | e.re < 1} := by
    intro y hy
    rw [Metric.mem_closedBall, Complex.dist_eq] at hy
    have := Complex.abs_re_le_norm (y - x)
    simp only [Complex.sub_re] at this
    simp only [Set.mem_ofPred_eq]
    linarith [(abs_le.1 (this.trans hy)).2, (abs_lt.1 hxre).2]
  have hdiff : DifferentiableOn ℂ (fun e => gC e t) {e : ℂ | e.re < 1} := fun e he =>
    (gC_hasDerivAt he ht0 ht).differentiableAt.differentiableWithinAt
  have hxU : x.re < 1 := hsub (Metric.mem_closedBall_self hδ0.le)
  rw [← (gC_hasDerivAt hxU ht0 ht).deriv]
  apply Complex.norm_deriv_le_of_forall_mem_sphere_norm_le hδ0 (hdiff.diffContOnCl_ball hsub)
  intro y hy
  have hyU : y.re < 1 := hsub (Metric.sphere_subset_closedBall hy)
  rw [mem_sphere_iff_norm] at hy
  have hyre : x₀.re - 2 * δ ≤ y.re := by
    have := Complex.abs_re_le_norm (y - x)
    simp only [Complex.sub_re] at this
    linarith [(abs_le.1 (this.trans hy.le)).1, (abs_lt.1 hxre).1]
  refine (proof_3 y hyU t ht0 ht).trans ?_
  exact Real.rpow_le_rpow_of_exponent_le (by norm_num) (by linarith)

lemma FC_hasDerivAt {x₀ : ℂ} (h : x₀.re < 1) :
    HasDerivAt FC (2 * ∫ t in (0:ℝ)..(1 / 2), gCd x₀ t) x₀ := by
  have hδ0 : 0 < (1 - x₀.re) / 3 := by linarith
  have key := intervalIntegral.hasDerivAt_integral_of_dominated_loc_of_deriv_le
    (μ := volume) (a := 0) (b := 1 / 2) (F := fun e t => gC e t) (F' := gCd)
    (x₀ := x₀) (s := Metric.ball x₀ ((1 - x₀.re) / 3))
    (bound := fun _ => (2 : ℝ) ^ (2 - x₀.re + 2 * ((1 - x₀.re) / 3)) / ((1 - x₀.re) / 3))
    (Metric.ball_mem_nhds x₀ hδ0)
    (Filter.Eventually.of_forall fun x => (gC_meas x).aestronglyMeasurable)
    (gC_iInt h) ((gCd_meas x₀).aestronglyMeasurable)
    (Filter.Eventually.of_forall fun t ht x hx => by
      rw [uIoc_of_le (by norm_num)] at ht
      exact gCd_bound h hx ht.1 ht.2)
    intervalIntegrable_const
    (Filter.Eventually.of_forall fun t ht x hx => by
      rw [uIoc_of_le (by norm_num)] at ht
      have hxre : x.re < 1 := by
        have := Complex.abs_re_le_norm (x - x₀)
        rw [Metric.mem_ball, Complex.dist_eq] at hx
        simp only [Complex.sub_re] at this
        linarith [(abs_le.1 (this.trans hx.le)).2]
      exact gC_hasDerivAt hxre ht.1 ht.2)
  exact key.2.const_mul (2 : ℂ)

end W02

theorem W02.proof_5 : W02.stmt_5 :=
  ⟨fun _ hε => W02.gC_iInt hε, fun _ hx₀ =>
    (W02.FC_hasDerivAt hx₀).differentiableAt.differentiableWithinAt⟩

theorem W02.proof_6 : W02.stmt_6 := by
  intro p hp
  have hp0 : (0 : ℝ) < p := by linarith
  have h1 : (1 / p : ℝ) < 1 := by rw [div_lt_one hp0]; exact hp
  have hz : ((1 / p : ℝ) : ℂ) ∈ Metric.ball (0 : ℂ) 1 := by
    rw [mem_ball_zero_iff, Complex.norm_real, Real.norm_eq_abs, abs_of_pos (by positivity)]
    exact h1
  have hf : DifferentiableOn ℂ W02.FC (Metric.ball 0 1) := by
    refine W02.proof_5.2.mono fun e he => ?_
    rw [mem_ball_zero_iff] at he
    have := Complex.re_le_norm e
    simp only [Set.mem_ofPred_eq]; linarith
  have H := Complex.hasSum_taylorSeries_on_ball hf hz
  have H2 := H.mul_left (-(1 / (2 * (p : ℂ))))
  have h4 := (W02.proof_4 (1 / p) h1).2
  have h2 := (W02.proof_2 p hp).2.2
  have hpc : (p : ℂ) ≠ 0 := by exact_mod_cast hp0.ne'
  convert H2 using 1
  · funext k
    simp only [smul_eq_mul, sub_zero, W02.Mk]
    push_cast
    rw [div_pow, one_pow]
    have : (p : ℂ) ^ k ≠ 0 := pow_ne_zero _ hpc
    field_simp
    ring
  · rw [h4]
    have : W02.Pint p / 8 - 1 = -(1 / (2 * p)) * W02.FR (1 / p) := by
      rw [show W02.Pint p / 8 - 1 = (W02.Pint p / 4) / 2 - 1 by ring, h2]; ring
    rw [this]; push_cast; ring

namespace W02

/-- the integrand of eq:M1t. -/
noncomputable def hM (t : ℝ) : ℝ := Real.log t / t * ((1 - t) / q t - 1)

lemma q_pos (t : ℝ) : 0 < q t := by
  unfold q; apply Real.sqrt_pos.2; nlinarith [sq_nonneg (t - 1/2)]

lemma q_sq (t : ℝ) : q t ^ 2 = t ^ 2 + (1 - t) ^ 2 := by
  unfold q; rw [Real.sq_sqrt (by positivity)]

lemma q_symm (t : ℝ) : q (1 - t) = q t := by
  unfold q; ring_nf

lemma hM_alt {t : ℝ} (ht0 : 0 < t) (ht1 : t < 1) :
    hM t = -Real.log t * t / (q t * (q t + 1 - t)) := by
  have hq := q_pos t
  have hq2 := q_sq t
  have hd : 0 < q t + 1 - t := by linarith
  unfold hM
  field_simp
  linear_combination (-(Real.log t)) * hq2

lemma hM_bound {t : ℝ} (ht0 : 0 < t) (ht1 : t < 1) : 0 ≤ hM t ∧ hM t ≤ 2 := by
  have hq := q_pos t
  have hq2 := q_sq t
  have hd : 0 < q t + 1 - t := by linarith
  have hl : Real.log t < 0 := Real.log_neg ht0 ht1
  have hlt := Real.abs_log_mul_self_lt t ht0 ht1.le
  rw [abs_lt] at hlt
  have hden : 1 / 2 ≤ q t * (q t + 1 - t) := by nlinarith
  rw [hM_alt ht0 ht1]
  constructor
  · apply div_nonneg _ (by positivity); nlinarith
  · rw [div_le_iff₀ (by positivity)]; nlinarith

lemma hM_meas : Measurable hM := by
  unfold hM q; fun_prop

lemma hM_iInt : IntervalIntegrable hM volume 0 1 := by
  refine IntervalIntegrable.mono_fun' (g := fun _ => (2 : ℝ))
    intervalIntegrable_const hM_meas.aestronglyMeasurable ?_
  rw [uIoc_of_le zero_le_one, ← restrict_Ioo_eq_restrict_Ioc]
  refine ae_restrict_of_forall_mem measurableSet_Ioo fun t ht => ?_
  have := hM_bound ht.1 ht.2
  simp only [Real.norm_eq_abs, abs_of_nonneg this.1]
  exact this.2

lemma hM_half : ∫ t in (0:ℝ)..(1/2), (hM t + hM (1 - t)) = ∫ t in (0:ℝ)..1, hM t := by
  have hi := hM_iInt
  have h1 : IntervalIntegrable hM volume 0 (1/2) :=
    hi.mono_set (by rw [uIcc_of_le (by norm_num), uIcc_of_le (by norm_num)]
                    exact Icc_subset_Icc le_rfl (by norm_num))
  have h2 : IntervalIntegrable hM volume (1/2) 1 :=
    hi.mono_set (by rw [uIcc_of_le (by norm_num), uIcc_of_le (by norm_num)]
                    exact Icc_subset_Icc (by norm_num) le_rfl)
  have h3 : IntervalIntegrable (fun t => hM (1 - t)) volume 0 (1/2) := by
    have := h2.comp_sub_left 1
    norm_num at this
    exact this.symm
  rw [intervalIntegral.integral_add h1 h3, ← intervalIntegral.integral_add_adjacent_intervals h1 h2,
    intervalIntegral.integral_comp_sub_left hM 1]
  norm_num

end W02

namespace W02

lemma sqrt_cpow (z : ℝ) (hz : 0 ≤ z) :
    (1 + (z : ℂ)) ^ (1 / 2 : ℂ) = ((Real.sqrt (1 + z) : ℝ) : ℂ) := by
  rw [Real.sqrt_eq_rpow, Complex.ofReal_cpow (by linarith)]
  push_cast; rfl

lemma sqrt_inv_cpow (z : ℝ) (hz : 0 ≤ z) :
    (1 + (z : ℂ)) ^ ((1 / 2 : ℂ) - 1) = (((Real.sqrt (1 + z))⁻¹ : ℝ) : ℂ) := by
  rw [show (1 / 2 : ℂ) - 1 = ((-(1 / 2) : ℝ) : ℂ) by push_cast; ring,
    Real.sqrt_eq_rpow, ← Real.rpow_neg (by linarith), Complex.ofReal_cpow (by linarith)]
  push_cast; rfl

lemma real_core {t : ℝ} (ht0 : 0 < t) (ht1 : t < 1) :
    aR 0 t * Real.log t + bR 0 t * Real.log (1 - t) -
      (aR 0 t * Real.log t * Real.sqrt (1 + (t / (1 - t)) ^ (2 - 2 * (0:ℝ))) +
        aR 0 t * ((1 / 2) * (Real.sqrt (1 + (t / (1 - t)) ^ (2 - 2 * (0:ℝ))))⁻¹ *
          ((t / (1 - t)) ^ (2 - 2 * (0:ℝ)) * Real.log (t / (1 - t)) * (-2)))) =
      -(hM t + hM (1 - t)) := by
  have h1 : 0 < 1 - t := by linarith
  have hq := q_pos t
  have hq2 := q_sq t
  have ha : aR 0 t = t⁻¹ := by
    unfold aR; rw [show (0:ℝ) - 1 = -1 by ring, Real.rpow_neg_one]
  have hb : bR 0 t = (1 - t)⁻¹ := by
    unfold bR; rw [show (0:ℝ) - 1 = -1 by ring, Real.rpow_neg_one]
  have hz : (t / (1 - t)) ^ (2 - 2 * (0:ℝ)) = (t / (1 - t)) ^ 2 := by
    rw [show (2:ℝ) - 2 * 0 = ((2:ℕ):ℝ) by norm_num, Real.rpow_natCast]
  have hs : Real.sqrt (1 + (t / (1 - t)) ^ 2) = q t / (1 - t) := by
    rw [show 1 + (t / (1 - t)) ^ 2 = (q t / (1 - t)) ^ 2 by
      rw [div_pow, div_pow, hq2]; field_simp; ring]
    exact Real.sqrt_sq (by positivity)
  have hl : Real.log (t / (1 - t)) = Real.log t - Real.log (1 - t) :=
    Real.log_div ht0.ne' h1.ne'
  have hM1 : hM (1 - t) = Real.log (1 - t) / (1 - t) * (t / q t - 1) := by
    unfold hM; rw [q_symm, sub_sub_cancel]
  rw [ha, hb, hz, hs, hl, hM1]
  unfold hM
  field_simp
  linear_combination (-(Real.log t)) * hq2

lemma gCd_zero {t : ℝ} (ht0 : 0 < t) (ht : t ≤ 1 / 2) :
    gCd 0 t = -(((hM t + hM (1 - t) : ℝ)) : ℂ) := by
  have h1 : t < 1 := by linarith
  have h1' : 0 < 1 - t := by linarith
  have hz0 : 0 ≤ (t / (1 - t)) ^ (2 - 2 * (0:ℝ)) :=
    Real.rpow_nonneg (div_nonneg ht0.le h1'.le) _
  have hc := real_core ht0 h1
  rw [← Complex.ofReal_zero]
  unfold gCd
  rw [aC_ofReal 0 ht0.le, bC_ofReal 0 h1.le, zC_ofReal 0 ht0.le h1, sqrt_cpow _ hz0,
    sqrt_inv_cpow _ hz0, ← Complex.ofReal_log ht0.le, ← Complex.ofReal_log h1'.le,
    ← Complex.ofReal_log (div_nonneg ht0.le h1'.le)]
  have hc' := congrArg Complex.ofReal hc
  push_cast at hc' ⊢
  linear_combination hc'

end W02

theorem W02.proof_7 : W02.stmt_7 := by
  refine ⟨W02.hM_iInt, ?_⟩
  have hD := (W02.FC_hasDerivAt (x₀ := 0) (by simp)).deriv
  unfold W02.Mk
  rw [iteratedDeriv_one, hD]
  have : ∫ t in (0:ℝ)..(1 / 2), W02.gCd 0 t =
      ∫ t in (0:ℝ)..(1 / 2), -(((W02.hM t + W02.hM (1 - t) : ℝ)) : ℂ) := by
    refine intervalIntegral.integral_congr_ae (Filter.Eventually.of_forall fun t ht => ?_)
    rw [uIoc_of_le (by norm_num)] at ht
    exact W02.gCd_zero ht.1 ht.2
  rw [this, intervalIntegral.integral_neg, intervalIntegral.integral_ofReal, W02.hM_half]
  rw [show (∫ t in (0:ℝ)..1, W02.hM t) =
      ∫ t in (0:ℝ)..1, Real.log t / t * ((1 - t) / W02.q t - 1) from rfl]
  simp only [Nat.factorial_one, Nat.cast_one]
  ring
