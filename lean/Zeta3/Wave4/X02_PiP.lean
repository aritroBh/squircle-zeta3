import Mathlib
/-!
# X02 — the expansion of `π_p` to order `p^{-4}` (paper label `thm:pip`)

Source: `paper/ieee/main.tex`, definition of `π_p` (lines 1330–1334) and Theorem `thm:pip`
("Expansion of $\pi_p$"), lines 1336–1350; detailed proof in `paper/ieee/proofs/sec_pip.tex`:
Lemma `lem:pip:split` (lines 52–66), eq. `eq:pip:Ehat` (line 71), the function `W` (line 79),
Lemma `lem:pip:W` (lines 83–102), Lemma `lem:pip:expansion` (lines 105–135),
Lemma `lem:pip:J` (lines 139–159), Lemma `prop:pip:first` (lines 162–189).

Definition (main.tex line 1332):

> \pi_p=\frac2p\int_0^1\bigl(t^{1-p}+(1-t)^{1-p}\bigr)^{1/p}dt

The claims formalized (main.tex lines 1337–1344):

> $\pi_p=4\cdot2^{-1/p}+E(p)$ with an explicit integral $E(p)$ (given in the proof), and
> \[ \pi_p=4\Bigl(1+\frac{N_0}{p}+\frac{N_1}{p^2}+\frac{N_2}{p^3}+\frac{N_3}{p^4}+O(p^{-5})\Bigr), \]
> where $N_0=-\log2$, $N_1=\frac12\log^22$, $N_2=\frac{\pi^2}{24}-\frac{\log^32}{6}$ and
> \[ N_3=\frac{\log^42}{24}+\frac{\pi^2}{24}(1-\log2)+\frac{\zeta(3)}{4} . \]

The explicit integral (sec_pip.tex, Lemma `lem:pip:split`, lines 52–57):

> For every $p>1$,
> \[ \pi_p=4\cdot2^{-1/p}+E(p),\qquad
> E(p)=\frac{4}{p(p-1)}\int_0^\infty \eu^{-\frac{s}{p(p-1)}}\bigl(1+\eu^{-\frac{s}{p-1}}\bigr)^{-1-\frac1p}
> \Bigl[\bigl(1+\eu^{-s}\bigr)^{1/p}-1\Bigr]ds . \]

Intermediate steps used (sec_pip.tex):

> (line 79) W(x,y)=\exp\Bigl(-\frac{xy}{1-x}\Bigr)\Bigl(1+\exp\Bigl(-\frac y{1-x}\Bigr)\Bigr)^{-1-x}.
> (line 118) Write $\hat E=\frac{h^2}{1-h}Q(h)$ with
>   $Q(h)=\int_0^\infty W(h,sh)\bigl[(1+\eu^{-s})^h-1\bigr]ds$.
> (lines 89–91, Lemma lem:pip:W with K = 1, and line 178)
>   \Bigl|W(h,sh)-\sum_{m\le K}h^mW_m(s)\Bigr|\le C_K\bigl(h(1+s)\bigr)^{K+1}\qquad(0<h\le h_K,\ 0\le s\le h^{-1/2}),
>   W_0(s)=\tfrac12,\qquad W_1(s)=\tfrac s4-\tfrac{\log2}2 .
> (line 46) J_{r,j}=\int_0^\infty s^j\mathcal L(s)^r\,ds,  $\mathcal L(s)=\log(1+\eu^{-s})$
> (lines 141–143) J_{1,0}=\frac{\pi^2}{12},\quad J_{1,1}=\frac{3\zeta(3)}4, ...  (ii) $J_{2,0}=\frac{\zeta(3)}4$.
> (lines 180–186) $q_1=\int W_0\mathcal L=\frac12J_{1,0}$ ... $q_1=\frac{\pi^2}{24}$ and
>   q_2=\frac{3\zeta(3)}{16}-\frac{\pi^2\log2}{24}+\frac{\zeta(3)}{16}=\frac{\zeta(3)}4-\frac{\pi^2\log2}{24}.
> (lines 165–167) E(p)=\frac{\pi^2}6p^{-3}+\Bigl(\frac{\pi^2}6+\zeta(3)-\frac{\pi^2\log2}6\Bigr)p^{-4}+O(p^{-5}).

Conventions.
* `ζ(3)` is written `(riemannZeta 3).re`, the real part of Mathlib's `riemannZeta 3`. Since
  `riemannZeta 3 = ∑' n, 1/(n:ℂ)^3` (`zeta_nat_eq_tsum_of_gt_one`) is real, this is the real
  number `∑_{n≥1} n^{-3}`.
* `π_p` is the interval integral `∫ t in 0..1`. Its integrand is `t^{-1+1/p}`-singular at `0`
  and `1` but integrable for `p > 1`; `stmt_1` asserts this integrability, so no statement
  holds only by the convention "non-integrable ⇒ integral = 0". The value of the integrand at
  the two endpoints (Lean's `0^{1-p} = 0`) is irrelevant (measure zero).
* `E(p)`, `Q(h)` and the `J_{r,j}` are integrals over `Ioi 0`; `stmt_1` and `stmt_2`–`stmt_4`
  assert integrability of the respective integrands. All bases of `rpow` are positive.
* The `O(·)` claims are stated as explicit bounds `∃ C, ∃ p₀ (or h₀ > 0), ∀ … , |error| ≤ C·…`;
  the paper prints no constants.
* `h = 1/p`. The bound for `W` (stmt_5) is Lemma `lem:pip:W` with `K = 1` and the paper's
  `W_0, W_1`, on the paper's range `0 ≤ s ≤ h^{-1/2}`.

NOT COVERED:
* Holomorphy of `W` on the bidisc and the all-orders bound of Lemma `lem:pip:W` (only `K = 1`).
* The all-orders expansion (Lemma `lem:pip:expansion` for general `K`, formula `eq:pip:Nk`),
  `N_4`, `N_5`, the divergence of the series, and the large-`k` asymptotics of `N_k`.
* `J_{1,2}`, `J_{1,3}` and the general formula `J_{1,j} = j! η(j+2)`.
* The main term `4·2^{-1/p}` expansion is done separately in `Zeta3.Wave1.T10_PiP`; here it is
  only used implicitly inside the final statement.
-/
-- BEGIN STATEMENTS
namespace X02

/-- `ζ(3)` as a real number: the real part of Mathlib's `riemannZeta 3` (which is real). -/
noncomputable def z3 : ℝ := (riemannZeta 3).re

/-- The integrand of `π_p`: `(t^{1-p} + (1-t)^{1-p})^{1/p}` (main.tex line 1332). -/
noncomputable def piIntegrand (p t : ℝ) : ℝ :=
  (t ^ (1 - p) + (1 - t) ^ (1 - p)) ^ (1 / p)

/-- `π_p = (2/p) ∫_0^1 (t^{1-p} + (1-t)^{1-p})^{1/p} dt` (main.tex line 1332). -/
noncomputable def piP (p : ℝ) : ℝ :=
  (2 / p) * ∫ t in (0 : ℝ)..1, piIntegrand p t

/-- The integrand of `E(p)` (sec_pip.tex Lemma `lem:pip:split`):
`e^{-s/(p(p-1))} (1 + e^{-s/(p-1)})^{-1-1/p} [(1+e^{-s})^{1/p} - 1]`. -/
noncomputable def EIntegrand (p s : ℝ) : ℝ :=
  Real.exp (-(s / (p * (p - 1)))) * (1 + Real.exp (-(s / (p - 1)))) ^ (-1 - 1 / p) *
    ((1 + Real.exp (-s)) ^ (1 / p) - 1)

/-- `E(p) = 4/(p(p-1)) ∫_0^∞ EIntegrand p s ds` (sec_pip.tex Lemma `lem:pip:split`). -/
noncomputable def E (p : ℝ) : ℝ :=
  4 / (p * (p - 1)) * ∫ s in Set.Ioi (0 : ℝ), EIntegrand p s

/-- `W(x,y) = exp(-xy/(1-x)) (1 + exp(-y/(1-x)))^{-1-x}` (sec_pip.tex line 79). -/
noncomputable def W (x y : ℝ) : ℝ :=
  Real.exp (-(x * y / (1 - x))) * (1 + Real.exp (-(y / (1 - x)))) ^ (-1 - x)

/-- `Q(h) = ∫_0^∞ W(h,sh) [(1+e^{-s})^h - 1] ds` (sec_pip.tex line 118). -/
noncomputable def Q (h : ℝ) : ℝ :=
  ∫ s in Set.Ioi (0 : ℝ), W h (s * h) * ((1 + Real.exp (-s)) ^ h - 1)

/-- `𝓛(s) = log(1 + e^{-s})` (sec_pip.tex line 46). -/
noncomputable def calL (s : ℝ) : ℝ := Real.log (1 + Real.exp (-s))

/-- `N_0 = -log 2` (main.tex line 1341). -/
noncomputable def N0 : ℝ := -Real.log 2
/-- `N_1 = ½ log² 2` (main.tex line 1341). -/
noncomputable def N1 : ℝ := (1 / 2) * Real.log 2 ^ 2
/-- `N_2 = π²/24 - log³2/6` (main.tex line 1341). -/
noncomputable def N2 : ℝ := Real.pi ^ 2 / 24 - Real.log 2 ^ 3 / 6
/-- `N_3 = log⁴2/24 + (π²/24)(1 - log 2) + ζ(3)/4` (main.tex line 1343). -/
noncomputable def N3 : ℝ :=
  Real.log 2 ^ 4 / 24 + Real.pi ^ 2 / 24 * (1 - Real.log 2) + z3 / 4

/-- Lemma `lem:pip:split` (first sentence of thm:pip): for every `p > 1` the integrands of
`π_p` (on `[0,1]`) and of `E(p)` (on `(0,∞)`) are integrable, and `π_p = 4·2^{-1/p} + E(p)`. -/
def stmt_1 : Prop :=
  ∀ p : ℝ, 1 < p →
    IntervalIntegrable (piIntegrand p) MeasureTheory.volume 0 1 ∧
    MeasureTheory.IntegrableOn (EIntegrand p) (Set.Ioi 0) ∧
    piP p = 4 * (2 : ℝ) ^ (-(1 / p)) + E p

/-- Lemma `lem:pip:J`(i), `j = 0`: `J_{1,0} = ∫_0^∞ log(1+e^{-s}) ds = π²/12`. -/
def stmt_2 : Prop :=
  MeasureTheory.IntegrableOn calL (Set.Ioi 0) ∧
    ∫ s in Set.Ioi (0 : ℝ), calL s = Real.pi ^ 2 / 12

/-- Lemma `lem:pip:J`(i), `j = 1`: `J_{1,1} = ∫_0^∞ s log(1+e^{-s}) ds = 3ζ(3)/4`. -/
def stmt_3 : Prop :=
  MeasureTheory.IntegrableOn (fun s => s * calL s) (Set.Ioi 0) ∧
    ∫ s in Set.Ioi (0 : ℝ), s * calL s = 3 * z3 / 4

/-- Lemma `lem:pip:J`(ii): `J_{2,0} = ∫_0^∞ log²(1+e^{-s}) ds = ζ(3)/4`. -/
def stmt_4 : Prop :=
  MeasureTheory.IntegrableOn (fun s => calL s ^ 2) (Set.Ioi 0) ∧
    ∫ s in Set.Ioi (0 : ℝ), calL s ^ 2 = z3 / 4

/-- Lemma `lem:pip:W` with `K = 1` and `W_0 = 1/2`, `W_1(s) = s/4 - log2/2`
(sec_pip.tex lines 89–91, 178): `|W(h,sh) - 1/2 - h(s/4 - log 2/2)| ≤ C (h(1+s))²`
for `0 < h ≤ h₀` and `0 ≤ s ≤ h^{-1/2}`. -/
def stmt_5 : Prop :=
  ∃ C h₀ : ℝ, 0 < h₀ ∧ ∀ h s : ℝ, 0 < h → h ≤ h₀ → 0 ≤ s → s ≤ 1 / Real.sqrt h →
    |W h (s * h) - (1 / 2 + h * (s / 4 - Real.log 2 / 2))| ≤ C * (h * (1 + s)) ^ 2

/-- Lemma `lem:pip:expansion` with `K = 2` and the values of `q_1, q_2` from Lemma
`prop:pip:first`: `Q(h) = (π²/24) h + (ζ(3)/4 - π² log 2/24) h² + O(h³)` as `h → 0⁺`. -/
def stmt_6 : Prop :=
  ∃ C h₀ : ℝ, 0 < h₀ ∧ ∀ h : ℝ, 0 < h → h ≤ h₀ →
    |Q h - (Real.pi ^ 2 / 24 * h + (z3 / 4 - Real.pi ^ 2 * Real.log 2 / 24) * h ^ 2)|
      ≤ C * h ^ 3

/-- Lemma `prop:pip:first`, expansion of `E`:
`E(p) = (π²/6) p^{-3} + (π²/6 + ζ(3) - π² log 2/6) p^{-4} + O(p^{-5})`. -/
def stmt_7 : Prop :=
  ∃ C p₀ : ℝ, 1 < p₀ ∧ ∀ p : ℝ, p₀ ≤ p →
    |E p - (Real.pi ^ 2 / 6 / p ^ 3
        + (Real.pi ^ 2 / 6 + z3 - Real.pi ^ 2 * Real.log 2 / 6) / p ^ 4)| ≤ C / p ^ 5

/-- Theorem `thm:pip` (main.tex lines 1338–1344):
`π_p = 4(1 + N_0/p + N_1/p² + N_2/p³ + N_3/p⁴ + O(p^{-5}))` as `p → ∞`. -/
def stmt_8 : Prop :=
  ∃ C p₀ : ℝ, 1 < p₀ ∧ ∀ p : ℝ, p₀ ≤ p →
    |piP p - 4 * (1 + N0 / p + N1 / p ^ 2 + N2 / p ^ 3 + N3 / p ^ 4)| ≤ C / p ^ 5

end X02
-- END STATEMENTS

open MeasureTheory Set Real

namespace X02


lemma summable_inv_pow (k : ℕ) (hk : 2 ≤ k) :
    Summable (fun n : ℕ => 1 / ((n : ℝ) + 1) ^ k) := by
  have := (Real.summable_one_div_nat_pow (p := k)).2 (by omega)
  exact (summable_nat_add_iff 1).2 this |>.congr (fun n => by push_cast; ring_nf)

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

lemma hasSum_zeta2 : HasSum (fun n : ℕ => 1 / ((n : ℝ) + 1) ^ 2) (π ^ 2 / 6) := by
  have := (hasSum_nat_add_iff' 1).2 hasSum_zeta_two
  simpa using this

/-- Alternating version: `Σ (-1)^n a_n = (1 - 2/2^k) Σ a_n` for `a_n = 1/(n+1)^k`. -/
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

lemma hasSum_eta2 : HasSum (fun n : ℕ => (-1) ^ n * (1 / ((n : ℝ) + 1) ^ 2)) (π ^ 2 / 12) := by
  have := hasSum_alt 2 hasSum_zeta2
  convert this using 1; ring

lemma hasSum_eta3 : HasSum (fun n : ℕ => (-1) ^ n * (1 / ((n : ℝ) + 1) ^ 3)) (3 * z3 / 4) := by
  have := hasSum_alt 3 hasSum_z3
  convert this using 1; ring


/-- `M(t) = -log(1 - e^{-t})`. -/
noncomputable def Mf (t : ℝ) : ℝ := -Real.log (1 - Real.exp (-t))

lemma intOn_pow_exp (j : ℕ) {b : ℝ} (hb : 0 < b) :
    IntegrableOn (fun s : ℝ => s ^ j * Real.exp (-(b * s))) (Ioi 0) := by
  have := integrableOn_rpow_mul_exp_neg_mul_rpow (s := (j : ℝ)) (p := 1) (b := b)
    (by have : (0 : ℝ) ≤ j := Nat.cast_nonneg j; linarith) one_pos hb
  refine this.congr_fun (fun x _ => ?_) measurableSet_Ioi
  simp [Real.rpow_natCast, neg_mul]

lemma integral_pow_exp (j : ℕ) {b : ℝ} (hb : 0 < b) :
    ∫ s in Ioi (0 : ℝ), s ^ j * Real.exp (-(b * s)) = (j.factorial : ℝ) / b ^ (j + 1) := by
  have := Real.integral_rpow_mul_exp_neg_mul_Ioi (a := (j : ℝ) + 1) (r := b) (by positivity) hb
  rw [Real.Gamma_nat_eq_factorial] at this
  have e : ∀ s : ℝ, s ^ ((j : ℝ) + 1 - 1) = s ^ j := fun s => by
    rw [add_sub_cancel_right, Real.rpow_natCast]
  simp_rw [e] at this
  rw [this, show ((j : ℝ) + 1) = ((j + 1 : ℕ) : ℝ) by push_cast; ring, Real.rpow_natCast]
  rw [one_div, inv_pow]; field_simp

lemma integral_series (c : ℕ → ℝ) (hc : ∀ n, |c n| ≤ 1) (j : ℕ) :
    ∫ s in Ioi (0 : ℝ), s ^ j * ∑' n : ℕ, c n * Real.exp (-(((n : ℝ) + 1) * s)) / ((n : ℝ) + 1)
      = (j.factorial : ℝ) * ∑' n : ℕ, c n / ((n : ℝ) + 1) ^ (j + 2) := by
  set F : ℕ → ℝ → ℝ := fun n s =>
    c n / ((n : ℝ) + 1) * (s ^ j * Real.exp (-(((n : ℝ) + 1) * s))) with hF
  have hpos : ∀ n : ℕ, (0 : ℝ) < (n : ℝ) + 1 := fun n => by positivity
  have hint : ∀ n, Integrable (F n) (volume.restrict (Ioi 0)) := fun n =>
    (intOn_pow_exp j (hpos n)).const_mul _
  have hval : ∀ n, ∫ s in Ioi (0 : ℝ), F n s
      = (j.factorial : ℝ) * (c n / ((n : ℝ) + 1) ^ (j + 2)) := fun n => by
    simp only [hF]
    rw [integral_const_mul, integral_pow_exp j (hpos n)]
    have := hpos n
    field_simp; ring
  have hnorm : ∀ n, ∫ s in Ioi (0 : ℝ), ‖F n s‖
      = (j.factorial : ℝ) * (|c n| / ((n : ℝ) + 1) ^ (j + 2)) := fun n => by
    have : ∫ s in Ioi (0 : ℝ), ‖F n s‖ = ∫ s in Ioi (0 : ℝ), (|c n| / ((n : ℝ) + 1) *
        (s ^ j * Real.exp (-(((n : ℝ) + 1) * s)))) := by
      refine setIntegral_congr_fun measurableSet_Ioi (fun s hs => ?_)
      simp only [hF, Real.norm_eq_abs, abs_mul, abs_div]
      rw [abs_of_pos (hpos n), abs_of_nonneg (pow_nonneg (le_of_lt (mem_Ioi.mp hs)) j : (0:ℝ) ≤ s ^ j),
        abs_of_pos (Real.exp_pos _)]
    rw [this, integral_const_mul, integral_pow_exp j (hpos n)]
    have := hpos n
    field_simp; ring
  have hsum : Summable fun n => ∫ s in Ioi (0 : ℝ), ‖F n s‖ := by
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

lemma calL_series {s : ℝ} (hs : 0 < s) :
    calL s = ∑' n : ℕ, (-1) ^ n * Real.exp (-(((n : ℝ) + 1) * s)) / ((n : ℝ) + 1) := by
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
  rw [calL, this]

lemma Mf_series {t : ℝ} (ht : 0 < t) :
    Mf t = ∑' n : ℕ, (1 : ℝ) * Real.exp (-(((n : ℝ) + 1) * t)) / ((n : ℝ) + 1) := by
  have hx : |Real.exp (-t)| < 1 := by
    rw [abs_of_pos (Real.exp_pos _)]; exact Real.exp_lt_one_iff.mpr (by linarith)
  have h := Real.hasSum_pow_div_log_of_abs_lt_one hx
  have e : ∀ n : ℕ, (Real.exp (-t)) ^ (n + 1) / ((n : ℝ) + 1)
      = (1 : ℝ) * Real.exp (-(((n : ℝ) + 1) * t)) / ((n : ℝ) + 1) := fun n => by
    rw [← Real.exp_nat_mul]; push_cast; ring_nf
  simp_rw [e] at h
  rw [Mf, h.tsum_eq]


lemma continuous_calL : Continuous calL := by
  unfold calL
  exact (continuous_const.add (Real.continuous_exp.comp continuous_neg)).log
    (fun x => (show (0:ℝ) < 1 + Real.exp (-x) by positivity).ne')

lemma calL_nonneg (s : ℝ) : 0 ≤ calL s :=
  Real.log_nonneg (by linarith [Real.exp_pos (-s)])

lemma calL_le (s : ℝ) : calL s ≤ Real.exp (-s) := by
  have := Real.log_le_sub_one_of_pos (show 0 < 1 + Real.exp (-s) by positivity)
  unfold calL; linarith

lemma intOn_of_bound {f g : ℝ → ℝ} (hf : Continuous f) (hg : IntegrableOn g (Ioi 0))
    (h : ∀ s : ℝ, 0 < s → |f s| ≤ g s) : IntegrableOn f (Ioi 0) := by
  refine Integrable.mono' hg hf.aestronglyMeasurable ?_
  exact (ae_restrict_iff' measurableSet_Ioi).2 (Filter.Eventually.of_forall
    (fun s hs => by rw [Real.norm_eq_abs]; exact h s hs))

lemma intOn_exp : IntegrableOn (fun s : ℝ => Real.exp (-s)) (Ioi 0) := by
  have := intOn_pow_exp 0 one_pos
  simpa using this

lemma intOn_calL : IntegrableOn calL (Ioi 0) :=
  intOn_of_bound continuous_calL intOn_exp (fun s _ => by
    rw [abs_of_nonneg (calL_nonneg s)]; exact calL_le s)

lemma intOn_s_calL : IntegrableOn (fun s => s * calL s) (Ioi 0) := by
  have := intOn_pow_exp 1 one_pos
  refine intOn_of_bound (continuous_id.mul continuous_calL) (g := fun s => s ^ 1 * Real.exp (-(1 * s)))
    this (fun s hs => ?_)
  rw [abs_of_nonneg (mul_nonneg hs.le (calL_nonneg s))]
  simp only [pow_one, one_mul]
  exact mul_le_mul_of_nonneg_left (calL_le s) hs.le

lemma intOn_calL_sq : IntegrableOn (fun s => calL s ^ 2) (Ioi 0) :=
  intOn_of_bound (continuous_calL.pow 2) intOn_exp (fun s hs => by
    rw [abs_of_nonneg (by positivity)]
    have h1 := calL_le s
    have h2 : Real.exp (-s) ≤ 1 := Real.exp_le_one_iff.mpr (by linarith)
    have h0 := calL_nonneg s
    nlinarith)

lemma J10 : ∫ s in Ioi (0 : ℝ), calL s = π ^ 2 / 12 := by
  have h := integral_series (fun n => (-1) ^ n) (fun n => by simp) 0
  have e : ∫ s in Ioi (0 : ℝ), calL s = ∫ s in Ioi (0 : ℝ), s ^ 0 * ∑' n : ℕ,
      (-1) ^ n * Real.exp (-(((n : ℝ) + 1) * s)) / ((n : ℝ) + 1) :=
    setIntegral_congr_fun measurableSet_Ioi (fun s hs => by rw [pow_zero, one_mul, calL_series hs])
  rw [e, h, ← hasSum_eta2.tsum_eq]
  simp [div_eq_mul_inv]

lemma J11 : ∫ s in Ioi (0 : ℝ), s * calL s = 3 * z3 / 4 := by
  have h := integral_series (fun n => (-1) ^ n) (fun n => by simp) 1
  have e : ∫ s in Ioi (0 : ℝ), s * calL s = ∫ s in Ioi (0 : ℝ), s ^ 1 * ∑' n : ℕ,
      (-1) ^ n * Real.exp (-(((n : ℝ) + 1) * s)) / ((n : ℝ) + 1) :=
    setIntegral_congr_fun measurableSet_Ioi (fun s hs => by rw [pow_one, calL_series hs])
  rw [e, h, ← hasSum_eta3.tsum_eq]
  simp [div_eq_mul_inv]

lemma JM : ∫ t in Ioi (0 : ℝ), t * Mf t = z3 := by
  have h := integral_series (fun _ => 1) (fun n => by simp) 1
  have e : ∫ t in Ioi (0 : ℝ), t * Mf t = ∫ t in Ioi (0 : ℝ), t ^ 1 * ∑' n : ℕ,
      (1 : ℝ) * Real.exp (-(((n : ℝ) + 1) * t)) / ((n : ℝ) + 1) :=
    setIntegral_congr_fun measurableSet_Ioi (fun s hs => by rw [pow_one, Mf_series hs])
  rw [e, h, ← hasSum_z3.tsum_eq]
  simp


lemma intOn_of_boundOn {f g : ℝ → ℝ} (hf : ContinuousOn f (Ioi 0)) (hg : IntegrableOn g (Ioi 0))
    (h : ∀ s : ℝ, 0 < s → |f s| ≤ g s) : IntegrableOn f (Ioi 0) := by
  refine Integrable.mono' hg (hf.aestronglyMeasurable measurableSet_Ioi) ?_
  exact (ae_restrict_iff' measurableSet_Ioi).2 (Filter.Eventually.of_forall
    (fun s hs => by rw [Real.norm_eq_abs]; exact h s hs))

lemma phi_eq (s : ℝ) : Real.log (1 + Real.exp s) = s + calL s := by
  have : 1 + Real.exp s = Real.exp s * (1 + Real.exp (-s)) := by
    rw [mul_add, ← Real.exp_add]; simp; ring
  rw [this, Real.log_mul (Real.exp_pos _).ne' (by positivity), Real.log_exp, calL]

lemma Mf_phi (s : ℝ) : Mf (Real.log (1 + Real.exp s)) = calL s := by
  have h1 : Real.exp (-Real.log (1 + Real.exp s)) = 1 / (1 + Real.exp s) := by
    rw [Real.exp_neg, Real.exp_log (by positivity), one_div]
  have h2 : 1 - 1 / (1 + Real.exp s) = Real.exp s / (1 + Real.exp s) := by
    field_simp; ring
  rw [Mf, h1, h2, Real.log_div (Real.exp_pos _).ne' (show (0:ℝ) < 1 + Real.exp s by positivity).ne', Real.log_exp, phi_eq]
  ring

lemma Mf_nonneg {t : ℝ} (ht : 0 < t) : 0 ≤ Mf t := by
  have h1 : Real.exp (-t) < 1 := Real.exp_lt_one_iff.mpr (by linarith)
  have h2 := Real.exp_pos (-t)
  unfold Mf
  have := Real.log_nonpos (by linarith : (0:ℝ) ≤ 1 - Real.exp (-t)) (by linarith)
  linarith

lemma tMf_le {t : ℝ} (ht : 0 < t) : t * Mf t ≤ (1 + t) * Real.exp (-t) := by
  set x := Real.exp (-t) with hx
  have h1 : x < 1 := Real.exp_lt_one_iff.mpr (by linarith)
  have h2 : 0 < x := Real.exp_pos _
  have hl := Real.one_sub_inv_le_log_of_pos (show 0 < 1 - x by linarith)
  have hM : Mf t ≤ x / (1 - x) := by
    unfold Mf; rw [← hx]
    have : 1 - (1 - x)⁻¹ = -(x / (1 - x)) := by field_simp; ring
    linarith
  have he : (1 + t) * x ≤ 1 := by
    have := Real.add_one_le_exp t
    rw [hx, Real.exp_neg]
    rw [mul_inv_le_iff₀ (Real.exp_pos t)]; linarith
  calc t * Mf t ≤ t * (x / (1 - x)) := mul_le_mul_of_nonneg_left hM ht.le
    _ ≤ (1 + t) * x := by
      rw [mul_div_assoc', div_le_iff₀ (by linarith)]
      nlinarith

lemma intOn_exp_lin : IntegrableOn (fun t : ℝ => (1 + t) * Real.exp (-t)) (Ioi 0) := by
  have h0 := intOn_pow_exp 0 one_pos
  have h1 := intOn_pow_exp 1 one_pos
  refine (h0.add h1).congr_fun (fun t _ => ?_) measurableSet_Ioi
  simp; ring

lemma continuousOn_Mf : ContinuousOn Mf (Ioi 0) := by
  intro t ht
  have h1 : Real.exp (-t) < 1 := Real.exp_lt_one_iff.mpr (by linarith [mem_Ioi.mp ht])
  unfold Mf
  have hc : ContinuousAt (fun t : ℝ => 1 - Real.exp (-t)) t := by fun_prop
  exact (hc.log (by linarith : (0:ℝ) < 1 - Real.exp (-t)).ne').neg.continuousWithinAt

lemma intOn_tMf : IntegrableOn (fun t => t * Mf t) (Ioi 0) :=
  intOn_of_boundOn (continuousOn_id.mul continuousOn_Mf) intOn_exp_lin (fun t ht => by
    rw [abs_of_nonneg (mul_nonneg ht.le (Mf_nonneg ht))]; exact tMf_le ht)

/-- `J_{1,1} + J_{2,0} = ζ(3)` (the Beta-function route of Remark rem:pip:C45x, done by the
substitution `t = log(1+e^s)`). -/
lemma JS : ∫ s in Ioi (0 : ℝ), (s * calL s + calL s ^ 2) = z3 := by
  set φ : ℝ → ℝ := fun s => Real.log (1 + Real.exp s) with hφ
  set σ : ℝ → ℝ := fun s => Real.exp s / (1 + Real.exp s) with hσ
  set f : ℝ → ℝ := fun s => φ s * calL s with hf
  have hder : ∀ x, HasDerivAt φ (σ x) x := fun x => by
    exact ((Real.hasDerivAt_exp x).const_add 1).log (by positivity)
  have hinj : InjOn φ univ := by
    intro a _ b _ hab
    by_contra hne
    rcases lt_or_gt_of_ne hne with h | h
    · have : φ a < φ b := Real.log_lt_log (by positivity)
        (by linarith [Real.exp_lt_exp.2 h])
      linarith
    · have : φ b < φ a := Real.log_lt_log (by positivity)
        (by linarith [Real.exp_lt_exp.2 h])
      linarith
  have himg : φ '' univ = Ioi 0 := by
    ext t; constructor
    · rintro ⟨s, _, rfl⟩
      exact Real.log_pos (by linarith [Real.exp_pos s])
    · intro ht
      have ht' : 1 < Real.exp t := Real.one_lt_exp_iff.mpr ht
      refine ⟨Real.log (Real.exp t - 1), mem_univ _, ?_⟩
      simp only [hφ]
      rw [Real.exp_log (by linarith)]; simp
  have hσpos : ∀ x, 0 < σ x := fun x => by simp only [hσ]; positivity
  have hσlt : ∀ x, σ x < 1 := fun x => by
    simp only [hσ]; rw [div_lt_one (by positivity)]; linarith
  have hGf : ∀ x, |σ x| • ((fun t => t * Mf t) (φ x)) = σ x * f x := fun x => by
    simp only [smul_eq_mul, hf, hφ, Mf_phi, abs_of_pos (hσpos x)]
  have hI := integral_image_eq_integral_abs_deriv_smul MeasurableSet.univ
    (fun x _ => (hder x).hasDerivWithinAt) hinj (fun t => t * Mf t)
  have hII := (integrableOn_image_iff_integrableOn_abs_deriv_smul MeasurableSet.univ
    (fun x _ => (hder x).hasDerivWithinAt) hinj (fun t => t * Mf t)).1
    (by rw [himg]; exact intOn_tMf)
  simp_rw [hGf] at hI hII
  rw [himg, JM, setIntegral_univ] at hI
  rw [integrableOn_univ] at hII
  -- `f` is even and equals `s 𝓛 + 𝓛²`
  have hfeq : ∀ x, f x = x * calL x + calL x ^ 2 := fun x => by
    simp only [hf, hφ, phi_eq]; ring
  have hfeven : ∀ x, f (-x) = f x := fun x => by
    simp only [hf, hφ, calL, neg_neg]; ring
  have hσneg : ∀ x, σ (-x) = 1 - σ x := fun x => by
    simp only [hσ, Real.exp_neg]
    field_simp; ring
  have hfint : IntegrableOn f (Ioi 0) := by
    have := intOn_s_calL.add intOn_calL_sq
    refine this.congr_fun (fun x _ => (hfeq x).symm) measurableSet_Ioi
  have hcont_f : Continuous f := by
    simp only [hf, hφ]
    exact ((continuous_const.add Real.continuous_exp).log
      (fun x => (show (0:ℝ) < 1 + Real.exp x by positivity).ne')).mul continuous_calL
  have hcont_σ : Continuous σ := by
    simp only [hσ]
    exact Real.continuous_exp.div (continuous_const.add Real.continuous_exp)
      (fun x => (show (0:ℝ) < 1 + Real.exp x by positivity).ne')
  have hfnn : ∀ x, 0 < x → 0 ≤ f x := fun x hx => by
    rw [hfeq]; have := calL_nonneg x; positivity
  have hA : IntegrableOn (fun x => σ x * f x) (Ioi 0) :=
    intOn_of_bound (hcont_σ.mul hcont_f) hfint (fun x hx => by
      rw [abs_of_nonneg (mul_nonneg (hσpos x).le (hfnn x hx))]
      nlinarith [hσlt x, hfnn x hx])
  have hB : IntegrableOn (fun x => (1 - σ x) * f x) (Ioi 0) :=
    intOn_of_bound ((continuous_const.sub hcont_σ).mul hcont_f) hfint (fun x hx => by
      rw [abs_of_nonneg (mul_nonneg (by linarith [hσlt x]) (hfnn x hx))]
      nlinarith [hσpos x, hfnn x hx])
  have hsplit : ∫ x, σ x * f x = (∫ x in Iic 0, σ x * f x) + ∫ x in Ioi 0, σ x * f x := by
    exact (intervalIntegral.integral_Iic_add_Ioi (b := 0) hII.integrableOn hII.integrableOn).symm
  have hneg : ∫ x in Iic 0, σ x * f x = ∫ x in Ioi 0, (1 - σ x) * f x := by
    have := integral_comp_neg_Ioi (0 : ℝ) (fun x => σ x * f x)
    rw [show -(0:ℝ) = 0 from neg_zero] at this
    rw [← this]
    refine setIntegral_congr_fun measurableSet_Ioi (fun x _ => ?_)
    simp only [hσneg, hfeven]
  have htot : ∫ x in Ioi 0, f x = ∫ x in Ioi 0, ((1 - σ x) * f x + σ x * f x) :=
    setIntegral_congr_fun measurableSet_Ioi (fun x _ => by ring)
  rw [integral_add hB hA] at htot
  have : ∫ x in Ioi (0:ℝ), (x * calL x + calL x ^ 2) = ∫ x in Ioi 0, f x :=
    setIntegral_congr_fun measurableSet_Ioi (fun x _ => (hfeq x).symm)
  linarith

lemma J20 : ∫ s in Ioi (0 : ℝ), calL s ^ 2 = z3 / 4 := by
  have h := JS
  rw [integral_add intOn_s_calL intOn_calL_sq, J11] at h
  linarith



lemma log_one_add_exp_neg (y : ℝ) :
    Real.log (1 + Real.exp (-y)) = Real.log 2 - y / 2 + Real.log (Real.cosh (y / 2)) := by
  have hc : 0 < Real.cosh (y / 2) := Real.cosh_pos _
  have : 1 + Real.exp (-y) = 2 * Real.exp (-(y / 2)) * Real.cosh (y / 2) := by
    rw [Real.cosh_eq, mul_assoc, mul_div_assoc', mul_add, ← Real.exp_add, ← Real.exp_add]
    ring_nf; rw [Real.exp_zero]; ring
  rw [this, Real.log_mul (by positivity) hc.ne', Real.log_mul (by norm_num) (by positivity),
    Real.log_exp]
  ring

lemma logcosh_bounds (z : ℝ) : 0 ≤ Real.log (Real.cosh z) ∧ Real.log (Real.cosh z) ≤ z ^ 2 / 2 := by
  constructor
  · exact Real.log_nonneg (Real.one_le_cosh z)
  · have := Real.cosh_le_exp_half_sq z
    calc Real.log (Real.cosh z) ≤ Real.log (Real.exp (z ^ 2 / 2)) :=
          Real.log_le_log (Real.cosh_pos z) this
      _ = z ^ 2 / 2 := Real.log_exp _

set_option maxHeartbeats 1000000 in
/-- The local bound: for `u = h(1+s) ≤ 1/4`, `|W(h,sh) - 1/2 - h(s/4 - log2/2)| ≤ 6u²`. -/
lemma W_local {h s : ℝ} (hh : 0 < h) (hs : 0 ≤ s) (hu : h * (1 + s) ≤ 1 / 4) :
    |W h (s * h) - (1 / 2 + h * (s / 4 - Real.log 2 / 2))| ≤ 6 * (h * (1 + s)) ^ 2 := by
  have ha0 : 0 ≤ s * h := mul_nonneg hs hh.le
  have hua : h * (1 + s) = h + s * h := by ring
  have hh4 : h ≤ 1 / 4 := by linarith
  have h1h : 0 < 1 - h := by linarith
  have hy0 : 0 ≤ s * h / (1 - h) := div_nonneg ha0 h1h.le
  have hya : s * h / (1 - h) - s * h = s * h * h / (1 - h) := by field_simp; ring
  have hL := log_one_add_exp_neg (s * h / (1 - h))
  obtain ⟨hc0, hc1⟩ := logcosh_bounds (s * h / (1 - h) / 2)
  have hW : W h (s * h) = Real.exp (-(h * (s * h / (1 - h))) + (-1 - h) *
      (Real.log 2 - s * h / (1 - h) / 2 + Real.log (Real.cosh (s * h / (1 - h) / 2)))
      + Real.log 2) / 2 := by
    unfold W
    rw [Real.rpow_def_of_pos (by positivity), hL, ← Real.exp_add, Real.exp_add _ (Real.log 2),
      Real.exp_log (by norm_num)]
    field_simp
  rw [hW]
  generalize hu' : h * (1 + s) = u at *
  generalize ha' : s * h = a at *
  generalize hy' : a / (1 - h) = y at *
  generalize hc' : Real.log (Real.cosh (y / 2)) = c at *
  generalize hX : -(h * y) + (-1 - h) * (Real.log 2 - y / 2 + c) + Real.log 2 = X
  have hya0 : 0 ≤ y - a := by rw [hya]; exact div_nonneg (mul_nonneg ha0 hh.le) h1h.le
  have hya1 : y - a ≤ 4 / 3 * (a * h) := by
    rw [hya, div_le_iff₀ h1h]; nlinarith [mul_nonneg ha0 hh.le]
  have hah4 : a * h ≤ a / 4 := by nlinarith
  have hu0 : 0 ≤ u := by linarith
  have hyu : y ≤ 4 / 3 * u := by linarith
  have hhu : h ≤ u := by linarith
  have hau : a ≤ u := by linarith
  have hlog2 : Real.log 2 < 1 := by
    have := Real.log_two_lt_d9; linarith
  have hlog2' : 0 < Real.log 2 := Real.log_pos (by norm_num)
  have hXX0 : X - (-(h * Real.log 2) + a / 2) = -(h * y) / 2 + (y - a) / 2 - (1 + h) * c := by
    rw [← hX]; ring
  have hhy : h * y ≤ 4 / 3 * u ^ 2 := by
    have := mul_le_mul hhu hyu hy0 hu0; nlinarith
  have hhy0 : 0 ≤ h * y := mul_nonneg hh.le hy0
  have hysq : y ^ 2 ≤ (4 / 3 * u) ^ 2 := pow_le_pow_left₀ hy0 hyu 2
  have hc1' : c ≤ 2 / 9 * u ^ 2 := by nlinarith
  have hc2 : (1 + h) * c ≤ u ^ 2 := by nlinarith
  have hc20 : 0 ≤ (1 + h) * c := mul_nonneg (by linarith) hc0
  have hah : a * h ≤ u ^ 2 := by
    have := mul_le_mul hau hhu hh.le hu0; nlinarith
  have hdiff : |X - (-(h * Real.log 2) + a / 2)| ≤ 3 * u ^ 2 := by
    rw [hXX0, abs_le]; constructor <;> linarith
  have hhl : h * Real.log 2 ≤ h := by nlinarith
  have hhl0 : 0 ≤ h * Real.log 2 := by positivity
  have hX0 : |-(h * Real.log 2) + a / 2| ≤ 3 / 2 * u := by
    rw [abs_le]; constructor <;> linarith
  have hXb : |X| ≤ 3 / 2 * u + 3 * u ^ 2 := by
    have := abs_sub_abs_le_abs_sub X (-(h * Real.log 2) + a / 2)
    linarith
  have hu2 : u ^ 2 ≤ u / 4 := by
    have := mul_le_mul_of_nonneg_left hu hu0; rw [sq]; linarith
  have hq : 3 / 2 * u + 3 * u ^ 2 ≤ 9 / 4 * u := by linarith
  have hX1 : |X| ≤ 1 := by linarith
  have hE := Real.abs_exp_sub_one_sub_id_le hX1
  have hXsq : X ^ 2 ≤ (9 / 4 * u) ^ 2 := by
    rw [← sq_abs]; exact pow_le_pow_left₀ (abs_nonneg _) (le_trans hXb hq) 2
  have key : Real.exp X / 2 - (1 / 2 + h * (s / 4 - Real.log 2 / 2))
      = (Real.exp X - 1 - X) / 2 + (X - (-(h * Real.log 2) + a / 2)) / 2 := by
    rw [← ha']; ring
  rw [key]
  have t1 : |(Real.exp X - 1 - X) / 2| ≤ (9 / 4 * u) ^ 2 / 2 := by
    rw [abs_div, abs_two]; linarith
  have t2 : |(X - (-(h * Real.log 2) + a / 2)) / 2| ≤ 3 * u ^ 2 / 2 := by
    rw [abs_div, abs_two]; linarith
  calc _ ≤ |(Real.exp X - 1 - X) / 2| + |(X - (-(h * Real.log 2) + a / 2)) / 2| := abs_add_le _ _
    _ ≤ (9 / 4 * u) ^ 2 / 2 + 3 * u ^ 2 / 2 := add_le_add t1 t2
    _ ≤ 6 * u ^ 2 := by nlinarith

lemma W_pos_le {h s : ℝ} (hh : 0 ≤ h) (h1 : h < 1) (hs : 0 ≤ s) :
    0 < W h (s * h) ∧ W h (s * h) ≤ 1 := by
  unfold W
  have h1h : 0 < 1 - h := by linarith
  have hy0 : 0 ≤ s * h / (1 - h) := div_nonneg (mul_nonneg hs hh) h1h.le
  have hb : 1 ≤ 1 + Real.exp (-(s * h / (1 - h))) := by linarith [Real.exp_pos (-(s * h / (1 - h)))]
  refine ⟨by positivity, ?_⟩
  have e1 : Real.exp (-(h * (s * h) / (1 - h))) ≤ 1 := by
    rw [Real.exp_le_one_iff, neg_nonpos]
    exact div_nonneg (mul_nonneg hh (mul_nonneg hs hh)) h1h.le
  have e2 : (1 + Real.exp (-(s * h / (1 - h)))) ^ (-1 - h) ≤ 1 :=
    Real.rpow_le_one_of_one_le_of_nonpos hb (by linarith)
  have e3 : 0 ≤ (1 + Real.exp (-(s * h / (1 - h)))) ^ (-1 - h) := by positivity
  calc _ ≤ 1 * 1 := mul_le_mul e1 e2 e3 zero_le_one
    _ = 1 := one_mul 1

/-- The global bound: `|W(h,sh) - 1/2 - h(s/4 - log2/2)| ≤ 32 h²(1+s)³` for `0 < h ≤ 1/4`. -/
lemma W_global {h s : ℝ} (hh : 0 < h) (hh4 : h ≤ 1 / 4) (hs : 0 ≤ s) :
    |W h (s * h) - (1 / 2 + h * (s / 4 - Real.log 2 / 2))| ≤ 32 * h ^ 2 * (1 + s) ^ 3 := by
  have hs1 : 1 ≤ 1 + s := by linarith
  rcases le_or_gt (h * (1 + s)) (1 / 4) with hu | hu
  · have := W_local hh hs hu
    have : 0 ≤ h ^ 2 * (1 + s) ^ 3 := by positivity
    have : (h * (1 + s)) ^ 2 ≤ h ^ 2 * (1 + s) ^ 3 := by
      rw [mul_pow]
      exact mul_le_mul_of_nonneg_left (by nlinarith) (by positivity)
    linarith
  · obtain ⟨w0, w1⟩ := W_pos_le hh.le (by linarith) hs
    have hlog2 : Real.log 2 < 1 := by
      have := Real.log_two_lt_d9; linarith
    have hlog2' : 0 < Real.log 2 := Real.log_pos (by norm_num)
    have hb : |W h (s * h) - (1 / 2 + h * (s / 4 - Real.log 2 / 2))| ≤ 2 * (1 + s) := by
      rw [abs_le]; constructor <;> nlinarith
    have h16 : 1 ≤ 16 * (h * (1 + s)) ^ 2 := by nlinarith
    have : 2 * (1 + s) ≤ 32 * h ^ 2 * (1 + s) ^ 3 := by
      have : 2 * (1 + s) * 1 ≤ 2 * (1 + s) * (16 * (h * (1 + s)) ^ 2) :=
        mul_le_mul_of_nonneg_left h16 (by linarith)
      nlinarith
    linarith


/-- Lemma `lem:pip:W` with `K = 1` on the paper's range `s ≤ h^{-1/2}`. -/
theorem W_stmt5 : ∃ C h₀ : ℝ, 0 < h₀ ∧ ∀ h s : ℝ, 0 < h → h ≤ h₀ → 0 ≤ s →
    s ≤ 1 / Real.sqrt h →
    |W h (s * h) - (1 / 2 + h * (s / 4 - Real.log 2 / 2))| ≤ C * (h * (1 + s)) ^ 2 := by
  refine ⟨6, 1 / 100, by norm_num, fun h s hh hh0 hs hs1 => ?_⟩
  have hsq : 0 < Real.sqrt h := Real.sqrt_pos.2 hh
  have hsq1 : Real.sqrt h ≤ 1 / 10 := by
    have : Real.sqrt (1 / 100) = 1 / 10 := by
      rw [show (1 / 100 : ℝ) = (1 / 10) ^ 2 by norm_num, Real.sqrt_sq (by norm_num)]
    rw [← this]; exact Real.sqrt_le_sqrt hh0
  have hhs : h * s ≤ Real.sqrt h := by
    calc h * s ≤ h * (1 / Real.sqrt h) := mul_le_mul_of_nonneg_left hs1 hh.le
      _ = Real.sqrt h := by
        have := Real.sq_sqrt hh.le
        field_simp
        linarith
  exact W_local hh hs (by nlinarith)


lemma piIntegrand_symm (p t : ℝ) : piIntegrand p (1 - t) = piIntegrand p t := by
  unfold piIntegrand; rw [sub_sub_cancel, add_comm]

lemma piIntegrand_le {p t : ℝ} (hp : 1 < p) (ht0 : 0 ≤ t) (ht1 : t ≤ 1) :
    piIntegrand p t ≤ t ^ (-1 + 1 / p) + (1 - t) ^ (-1 + 1 / p) := by
  have hp0 : 0 < p := by linarith
  unfold piIntegrand
  have e : (1 - p) * (1 / p) = -1 + 1 / p := by field_simp; ring
  have h1 := Real.rpow_add_le_add_rpow (Real.rpow_nonneg ht0 (1 - p))
    (Real.rpow_nonneg (by linarith : (0:ℝ) ≤ 1 - t) (1 - p)) (by positivity : (0:ℝ) ≤ 1 / p)
    (by rw [div_le_one hp0]; linarith)
  rw [← Real.rpow_mul ht0, ← Real.rpow_mul (by linarith), e] at h1
  exact h1

lemma measurable_piIntegrand (p : ℝ) : Measurable (piIntegrand p) := by
  unfold piIntegrand; fun_prop

lemma intervalIntegrable_bound {p : ℝ} (hp : 1 < p) (a b : ℝ) :
    IntervalIntegrable (fun t : ℝ => t ^ (-1 + 1 / p) + (1 - t) ^ (-1 + 1 / p)) volume a b := by
  have hp0 : 0 < p := by linarith
  have hr : -1 < -1 + 1 / p := by have : 0 < 1 / p := by positivity
                                  linarith
  have h1 : IntervalIntegrable (fun t : ℝ => t ^ (-1 + 1 / p)) volume a b :=
    intervalIntegral.intervalIntegrable_rpow' hr
  have h2 : IntervalIntegrable (fun t : ℝ => (1 - t) ^ (-1 + 1 / p)) volume a b := by
    have := (intervalIntegral.intervalIntegrable_rpow' (a := 1 - a) (b := 1 - b) hr).comp_sub_left 1
    simpa using this
  exact h1.add h2

lemma intervalIntegrable_piIntegrand {p : ℝ} (hp : 1 < p) :
    IntervalIntegrable (piIntegrand p) volume 0 1 := by
  refine (intervalIntegrable_bound hp 0 1).mono_fun' (measurable_piIntegrand p).aestronglyMeasurable ?_
  rw [uIoc_of_le zero_le_one]
  refine (ae_restrict_iff' measurableSet_Ioc).2 (Filter.Eventually.of_forall (fun t ht => ?_))
  have hnn : 0 ≤ piIntegrand p t := by
    unfold piIntegrand
    exact Real.rpow_nonneg (add_nonneg (Real.rpow_nonneg ht.1.le _)
      (Real.rpow_nonneg (by linarith [ht.2]) _)) _
  show ‖piIntegrand p t‖ ≤ _
  rw [Real.norm_eq_abs, abs_of_nonneg hnn]
  exact piIntegrand_le hp ht.1.le ht.2


/-- The part of the integrand of `π_p` beyond the main term `t^{-1+1/p}`. -/
noncomputable def Dp (p t : ℝ) : ℝ := piIntegrand p t - t ^ (-1 + 1 / p)

/-- The substitution `t = 1/(1+e^{s/(p-1)})`, pointwise. -/
lemma subst_pointwise {p : ℝ} (hp : 1 < p) (s : ℝ) :
    |-(Real.exp (s / (p - 1)) * (1 / (p - 1))) / (1 + Real.exp (s / (p - 1))) ^ 2| *
      Dp p (1 + Real.exp (s / (p - 1)))⁻¹ = 1 / (p - 1) * EIntegrand p s := by
  have hp0 : 0 < p := by linarith
  have hp1 : 0 < p - 1 := by linarith
  set E := Real.exp (s / (p - 1)) with hE_def
  have hE : 0 < E := Real.exp_pos _
  have h1E : 0 < 1 + E := by linarith
  set ℓ := Real.log (1 + E) with hℓ_def
  have hℓ : Real.exp ℓ = 1 + E := Real.exp_log h1E
  set lam := Real.log (1 + Real.exp (-s)) with hlam_def
  have hlam : Real.exp lam = 1 + Real.exp (-s) := Real.exp_log (by positivity)
  have ht : 0 < (1 + E)⁻¹ := inv_pos.2 h1E
  have hlt : Real.log (1 + E)⁻¹ = -ℓ := Real.log_inv _
  have h1t : 1 - (1 + E)⁻¹ = E * (1 + E)⁻¹ := by field_simp; ring
  have h1t0 : 0 < E * (1 + E)⁻¹ := mul_pos hE ht
  have hl1t : Real.log (E * (1 + E)⁻¹) = s / (p - 1) - ℓ := by
    rw [Real.log_mul hE.ne' ht.ne', hlt, hE_def, Real.log_exp]; ring
  have f1 : ((1 + E)⁻¹) ^ (1 - p) = Real.exp (ℓ * (p - 1)) := by
    rw [Real.rpow_def_of_pos ht, hlt]; congr 1; ring
  have f2 : (1 - (1 + E)⁻¹) ^ (1 - p) = Real.exp (ℓ * (p - 1) - s) := by
    rw [h1t, Real.rpow_def_of_pos h1t0, hl1t]; congr 1; field_simp; ring
  have f3 : ((1 + E)⁻¹) ^ (1 - p) + (1 - (1 + E)⁻¹) ^ (1 - p) = Real.exp (ℓ * (p - 1) + lam) := by
    rw [f1, f2, Real.exp_add, hlam, Real.exp_sub, Real.exp_neg]; field_simp
  have f4 : piIntegrand p (1 + E)⁻¹ = Real.exp ((ℓ * (p - 1) + lam) * (1 / p)) := by
    unfold piIntegrand; rw [f3, Real.rpow_def_of_pos (Real.exp_pos _), Real.log_exp]
  have f5 : ((1 + E)⁻¹) ^ (-1 + 1 / p) = Real.exp (-ℓ * (-1 + 1 / p)) := by
    rw [Real.rpow_def_of_pos ht, hlt]
  have f6 : (1 + Real.exp (-s)) ^ (1 / p) = Real.exp (lam * (1 / p)) := by
    rw [Real.rpow_def_of_pos (by positivity)]
  have hb : 1 + Real.exp (-(s / (p - 1))) = (1 + E) * E⁻¹ := by
    rw [Real.exp_neg, ← hE_def]; field_simp; ring
  have f7 : (1 + Real.exp (-(s / (p - 1)))) ^ (-1 - 1 / p)
      = Real.exp ((ℓ - s / (p - 1)) * (-1 - 1 / p)) := by
    rw [hb, Real.rpow_def_of_pos (mul_pos h1E (inv_pos.2 hE)),
      Real.log_mul h1E.ne' (inv_pos.2 hE).ne', Real.log_inv, hE_def, Real.log_exp]
    rfl
  have f8 : |-(E * (1 / (p - 1))) / (1 + E) ^ 2| = 1 / (p - 1) * Real.exp (s / (p - 1) - 2 * ℓ) := by
    rw [abs_div, abs_neg, abs_of_pos (by positivity), abs_of_pos (by positivity),
      Real.exp_sub, ← hE_def, show 2 * ℓ = ℓ + ℓ by ring, Real.exp_add, hℓ]
    field_simp
  unfold Dp EIntegrand
  rw [f8, f4, f5, f6, f7]
  have lhs : 1 / (p - 1) * Real.exp (s / (p - 1) - 2 * ℓ) *
      (Real.exp ((ℓ * (p - 1) + lam) * (1 / p)) - Real.exp (-ℓ * (-1 + 1 / p)))
      = 1 / (p - 1) * (Real.exp (s / (p - 1) - 2 * ℓ + (ℓ * (p - 1) + lam) * (1 / p))
        - Real.exp (s / (p - 1) - 2 * ℓ + -ℓ * (-1 + 1 / p))) := by
    rw [Real.exp_add, Real.exp_add]; ring
  have rhs : 1 / (p - 1) * (Real.exp (-(s / (p * (p - 1)))) *
      Real.exp ((ℓ - s / (p - 1)) * (-1 - 1 / p)) * (Real.exp (lam * (1 / p)) - 1))
      = 1 / (p - 1) * (Real.exp (-(s / (p * (p - 1))) + (ℓ - s / (p - 1)) * (-1 - 1 / p)
        + lam * (1 / p)) - Real.exp (-(s / (p * (p - 1))) + (ℓ - s / (p - 1)) * (-1 - 1 / p))) := by
    rw [Real.exp_add, Real.exp_add]; ring
  rw [lhs, rhs]
  have e1 : s / (p - 1) - 2 * ℓ + (ℓ * (p - 1) + lam) * (1 / p)
      = -(s / (p * (p - 1))) + (ℓ - s / (p - 1)) * (-1 - 1 / p) + lam * (1 / p) := by
    field_simp; ring
  have e2 : s / (p - 1) - 2 * ℓ + -ℓ * (-1 + 1 / p)
      = -(s / (p * (p - 1))) + (ℓ - s / (p - 1)) * (-1 - 1 / p) := by
    field_simp; ring
  rw [e1, e2]


lemma intervalIntegrable_mainterm {p : ℝ} (hp : 1 < p) (a b : ℝ) :
    IntervalIntegrable (fun t : ℝ => t ^ (-1 + 1 / p)) volume a b := by
  have hp0 : 0 < p := by linarith
  have : 0 < 1 / p := by positivity
  exact intervalIntegral.intervalIntegrable_rpow' (by linarith)

lemma subst_main {p : ℝ} (hp : 1 < p) :
    IntegrableOn (EIntegrand p) (Ioi 0) ∧
      ∫ t in (0 : ℝ)..(1 / 2), Dp p t = 1 / (p - 1) * ∫ s in Ioi (0 : ℝ), EIntegrand p s := by
  have hp1 : 0 < p - 1 := by linarith
  set σ : ℝ → ℝ := fun s => (1 + Real.exp (s / (p - 1)))⁻¹ with hσ
  set σ' : ℝ → ℝ := fun s =>
    -(Real.exp (s / (p - 1)) * (1 / (p - 1))) / (1 + Real.exp (s / (p - 1))) ^ 2 with hσ'
  have hder : ∀ s, HasDerivAt σ (σ' s) s := fun s => by
    have h1 := (((hasDerivAt_id' s).div_const (p - 1)).exp).const_add 1
    exact h1.inv (by positivity)
  have hinj : InjOn σ (Ioi 0) := by
    intro a _ b _ hab
    simp only [hσ] at hab
    have := inv_injective hab
    have h2 : Real.exp (a / (p - 1)) = Real.exp (b / (p - 1)) := by linarith
    have h3 := Real.exp_injective h2
    field_simp at h3; exact h3
  have himg : σ '' Ioi 0 = Ioo 0 (1 / 2) := by
    ext t; constructor
    · rintro ⟨s, hs, rfl⟩
      have h1 : 1 < Real.exp (s / (p - 1)) :=
        Real.one_lt_exp_iff.mpr (div_pos hs hp1)
      simp only [hσ]
      refine ⟨by positivity, ?_⟩
      rw [inv_lt_comm₀ (by positivity) (by norm_num)]; linarith
    · rintro ⟨h0, h1⟩
      have hq : 1 < 1 / t - 1 := by
        rw [lt_sub_iff_add_lt, lt_div_iff₀ h0]; linarith
      refine ⟨(p - 1) * Real.log (1 / t - 1), ?_, ?_⟩
      · exact mul_pos hp1 (Real.log_pos hq)
      · simp only [hσ]
        rw [mul_div_cancel_left₀ _ hp1.ne', Real.exp_log (by linarith)]
        field_simp; ring
  -- integrability of `Dp` on `(0, 1/2)`
  have hDi : IntervalIntegrable (Dp p) volume 0 (1 / 2) := by
    have h1 : IntervalIntegrable (piIntegrand p) volume 0 (1 / 2) :=
      (intervalIntegrable_piIntegrand hp).mono_set (by
        rw [uIcc_of_le (by norm_num), uIcc_of_le (by norm_num)]
        exact Icc_subset_Icc le_rfl (by norm_num))
    exact h1.sub (intervalIntegrable_mainterm hp 0 (1 / 2))
  have hDo : IntegrableOn (Dp p) (Ioo 0 (1 / 2)) :=
    ((intervalIntegrable_iff_integrableOn_Ioc_of_le (by norm_num)).1 hDi).mono_set
      Ioo_subset_Ioc_self
  have hpt : ∀ s, |σ' s| • Dp p (σ s) = 1 / (p - 1) * EIntegrand p s := fun s => by
    simp only [hσ, hσ', smul_eq_mul]; exact subst_pointwise hp s
  have hII := (integrableOn_image_iff_integrableOn_abs_deriv_smul measurableSet_Ioi
    (fun x _ => (hder x).hasDerivWithinAt) hinj (Dp p)).1 (by rw [himg]; exact hDo)
  have hI := integral_image_eq_integral_abs_deriv_smul measurableSet_Ioi
    (fun x _ => (hder x).hasDerivWithinAt) hinj (Dp p)
  simp_rw [hpt] at hII hI
  rw [himg] at hI
  refine ⟨?_, ?_⟩
  · have : IntegrableOn (fun s => (p - 1) * (1 / (p - 1) * EIntegrand p s)) (Ioi 0) :=
      hII.const_mul (p - 1)
    refine this.congr_fun (fun s _ => ?_) measurableSet_Ioi
    field_simp
  · rw [intervalIntegral.integral_of_le (by norm_num), integral_Ioc_eq_integral_Ioo, hI,
      integral_const_mul]

theorem split {p : ℝ} (hp : 1 < p) :
    IntervalIntegrable (piIntegrand p) MeasureTheory.volume 0 1 ∧
    MeasureTheory.IntegrableOn (EIntegrand p) (Set.Ioi 0) ∧
    piP p = 4 * (2 : ℝ) ^ (-(1 / p)) + E p := by
  have hp0 : 0 < p := by linarith
  have hp1 : 0 < p - 1 := by linarith
  obtain ⟨hEi, hEq⟩ := subst_main hp
  refine ⟨intervalIntegrable_piIntegrand hp, hEi, ?_⟩
  have hI := intervalIntegrable_piIntegrand hp
  have hI1 : IntervalIntegrable (piIntegrand p) volume 0 (1 / 2) :=
    hI.mono_set (by
      rw [uIcc_of_le (by norm_num), uIcc_of_le (by norm_num)]
      exact Icc_subset_Icc le_rfl (by norm_num))
  have hI2 : IntervalIntegrable (piIntegrand p) volume (1 / 2) 1 :=
    hI.mono_set (by
      rw [uIcc_of_le (by norm_num), uIcc_of_le (by norm_num)]
      exact Icc_subset_Icc (by norm_num) le_rfl)
  have hsym : ∫ t in (1 / 2 : ℝ)..1, piIntegrand p t = ∫ t in (0 : ℝ)..(1 / 2), piIntegrand p t := by
    have := intervalIntegral.integral_comp_sub_left (a := (1 / 2 : ℝ)) (b := 1)
      (piIntegrand p) 1
    simp only [piIntegrand_symm] at this
    rw [this]; norm_num
  have hsplit : ∫ t in (0 : ℝ)..1, piIntegrand p t
      = 2 * ∫ t in (0 : ℝ)..(1 / 2), piIntegrand p t := by
    rw [← intervalIntegral.integral_add_adjacent_intervals hI1 hI2, hsym]; ring
  have hmain : ∫ t in (0 : ℝ)..(1 / 2), t ^ (-1 + 1 / p) = p * (2 : ℝ) ^ (-(1 / p)) := by
    have hr : (-1 + 1 / p) + 1 = 1 / p := by ring
    rw [integral_rpow (Or.inl (by have : 0 < 1 / p := by positivity
                                  linarith))]
    rw [hr, Real.zero_rpow (by positivity)]
    have : ((1:ℝ) / 2) ^ (1 / p) = (2:ℝ) ^ (-(1 / p)) := by
      rw [Real.rpow_neg (by norm_num), one_div, Real.inv_rpow (by norm_num)]
    rw [this]
    field_simp
    simp
  have hD : ∫ t in (0 : ℝ)..(1 / 2), piIntegrand p t
      = (∫ t in (0 : ℝ)..(1 / 2), t ^ (-1 + 1 / p)) + ∫ t in (0 : ℝ)..(1 / 2), Dp p t := by
    unfold Dp
    rw [← intervalIntegral.integral_add (intervalIntegrable_mainterm hp 0 (1 / 2))
      (hI1.sub (intervalIntegrable_mainterm hp 0 (1 / 2)))]
    congr 1; funext t; ring
  unfold piP E
  rw [hsplit, hD, hmain, hEq]
  field_simp
  ring



lemma exp_taylor3 {x : ℝ} (h0 : 0 ≤ x) (h1 : x ≤ 1) :
    |Real.exp x - 1 - x - x ^ 2 / 2| ≤ x ^ 3 := by
  have := Real.exp_bound (x := x) (by rw [abs_of_nonneg h0]; exact h1) (n := 3) (by norm_num)
  simp only [Finset.sum_range_succ, Finset.sum_range_zero, Nat.factorial] at this
  rw [abs_of_nonneg h0] at this
  have e : Real.exp x - 1 - x - x ^ 2 / 2 = Real.exp x - (0 + x ^ 0 / ((1 : ℕ) : ℝ) +
      x ^ 1 / ((1 * 1 : ℕ) : ℝ) + x ^ 2 / ((2 * (1 * 1) : ℕ) : ℝ)) := by
    push_cast; ring
  rw [e]
  refine le_trans this ?_
  have : 0 ≤ x ^ 3 := by positivity
  norm_num
  nlinarith

lemma intOn_cube : IntegrableOn (fun s : ℝ => (1 + s) ^ 3 * Real.exp (-s)) (Ioi 0) := by
  have h0 := intOn_pow_exp 0 one_pos
  have h1 := intOn_pow_exp 1 one_pos
  have h2 := intOn_pow_exp 2 one_pos
  have h3 := intOn_pow_exp 3 one_pos
  have := ((h0.add (h1.const_mul 3)).add (h2.const_mul 3)).add h3
  refine this.congr_fun (fun s _ => ?_) measurableSet_Ioi
  simp only [Pi.add_apply, one_mul, pow_zero, pow_one]; ring

/-- The main part `G` of the integrand of `Q`. -/
noncomputable def Gq (h s : ℝ) : ℝ :=
  (h / 2 - h ^ 2 * Real.log 2 / 2) * calL s + h ^ 2 / 4 * (s * calL s) + h ^ 2 / 4 * calL s ^ 2

lemma Q_pointwise {h s : ℝ} (hh : 0 < h) (hh4 : h ≤ 1 / 4) (hs : 0 ≤ s) :
    |W h (s * h) * ((1 + Real.exp (-s)) ^ h - 1) - Gq h s|
      ≤ 66 * h ^ 3 * ((1 + s) ^ 3 * Real.exp (-s)) := by
  have hL0 := calL_nonneg s
  have hLe := calL_le s
  have he1 : Real.exp (-s) ≤ 1 := Real.exp_le_one_iff.mpr (by linarith)
  have he0 := Real.exp_pos (-s)
  have hL1 : calL s ≤ 1 := le_trans hLe he1
  have hg : (1 + Real.exp (-s)) ^ h - 1 = Real.exp (h * calL s) - 1 := by
    rw [Real.rpow_def_of_pos (by positivity), calL, mul_comm]
  rw [hg]
  obtain ⟨w0, w1⟩ := W_pos_le hh.le (by linarith) hs
  have hρ := W_global hh hh4 hs
  generalize hx : h * calL s = x
  have hx0 : 0 ≤ x := by rw [← hx]; positivity
  have hxL : x ≤ h * calL s := le_of_eq hx.symm
  have hx1 : x ≤ 1 / 4 := by
    rw [← hx]
    have := mul_le_mul hh4 hL1 hL0 (by norm_num : (0:ℝ) ≤ 1 / 4)
    linarith
  have hT := exp_taylor3 hx0 (by linarith)
  set ρ := W h (s * h) - (1 / 2 + h * (s / 4 - Real.log 2 / 2)) with hρ_def
  have key : W h (s * h) * (Real.exp x - 1) - Gq h s
      = W h (s * h) * (Real.exp x - 1 - x - x ^ 2 / 2) + ρ * (x + x ^ 2 / 2)
        + h ^ 3 * (s / 4 - Real.log 2 / 2) * calL s ^ 2 / 2 := by
    rw [hρ_def, Gq, ← hx]; ring
  rw [key]
  have hs1 : 1 ≤ 1 + s := by linarith
  have hc1 : 1 + s ≤ (1 + s) ^ 3 := le_self_pow₀ hs1 (by norm_num)
  have hc0 : 1 ≤ (1 + s) ^ 3 := le_trans hs1 hc1
  have hlog2 : Real.log 2 < 1 := by
    have := Real.log_two_lt_d9; linarith
  have hlog2' : 0 < Real.log 2 := Real.log_pos (by norm_num)
  -- term 1
  have t1 : |W h (s * h) * (Real.exp x - 1 - x - x ^ 2 / 2)| ≤ h ^ 3 * Real.exp (-s) := by
    rw [abs_mul, abs_of_pos w0]
    have hx3 : x ^ 3 ≤ h ^ 3 * Real.exp (-s) := by
      rw [← hx, mul_pow]
      apply mul_le_mul_of_nonneg_left _ (by positivity)
      have : calL s ^ 3 ≤ calL s := by
        have : calL s ^ 3 = calL s * calL s ^ 2 := by ring
        rw [this]; exact mul_le_of_le_one_right hL0 (by nlinarith)
      linarith
    calc _ ≤ 1 * x ^ 3 := mul_le_mul w1 hT (abs_nonneg _) zero_le_one
      _ ≤ _ := by linarith
  -- term 2
  have t2 : |ρ * (x + x ^ 2 / 2)| ≤ 64 * h ^ 3 * ((1 + s) ^ 3 * Real.exp (-s)) := by
    rw [abs_mul]
    have hxx : |x + x ^ 2 / 2| ≤ 2 * h * Real.exp (-s) := by
      rw [abs_of_nonneg (by positivity)]
      have : x ≤ h * Real.exp (-s) := by rw [← hx]; exact mul_le_mul_of_nonneg_left hLe hh.le
      nlinarith
    calc _ ≤ (32 * h ^ 2 * (1 + s) ^ 3) * (2 * h * Real.exp (-s)) :=
          mul_le_mul hρ hxx (abs_nonneg _) (by positivity)
      _ = _ := by ring
  -- term 3
  have t3 : |h ^ 3 * (s / 4 - Real.log 2 / 2) * calL s ^ 2 / 2|
      ≤ h ^ 3 * ((1 + s) * Real.exp (-s)) := by
    have hW1 : |s / 4 - Real.log 2 / 2| ≤ 1 + s := by
      rw [abs_le]; constructor <;> linarith
    have hL2 : calL s ^ 2 ≤ Real.exp (-s) := by nlinarith
    rw [abs_div, abs_two, abs_mul, abs_mul, abs_of_pos (by positivity : 0 < h ^ 3),
      abs_of_nonneg (by positivity : 0 ≤ calL s ^ 2)]
    have : h ^ 3 * |s / 4 - Real.log 2 / 2| * calL s ^ 2 ≤ h ^ 3 * ((1 + s) * Real.exp (-s)) := by
      rw [mul_assoc]
      apply mul_le_mul_of_nonneg_left _ (by positivity)
      exact mul_le_mul hW1 hL2 (by positivity) (by linarith)
    have : 0 ≤ h ^ 3 * |s / 4 - Real.log 2 / 2| * calL s ^ 2 := by positivity
    linarith
  have hb1 : h ^ 3 * Real.exp (-s) ≤ h ^ 3 * ((1 + s) ^ 3 * Real.exp (-s)) := by
    apply mul_le_mul_of_nonneg_left _ (by positivity); exact le_mul_of_one_le_left he0.le hc0
  have hb3 : h ^ 3 * ((1 + s) * Real.exp (-s)) ≤ h ^ 3 * ((1 + s) ^ 3 * Real.exp (-s)) := by
    apply mul_le_mul_of_nonneg_left _ (by positivity); exact mul_le_mul_of_nonneg_right hc1 he0.le
  calc _ ≤ |W h (s * h) * (Real.exp x - 1 - x - x ^ 2 / 2) + ρ * (x + x ^ 2 / 2)|
        + |h ^ 3 * (s / 4 - Real.log 2 / 2) * calL s ^ 2 / 2| := abs_add_le _ _
    _ ≤ |W h (s * h) * (Real.exp x - 1 - x - x ^ 2 / 2)| + |ρ * (x + x ^ 2 / 2)|
        + |h ^ 3 * (s / 4 - Real.log 2 / 2) * calL s ^ 2 / 2| := by
          gcongr; exact abs_add_le _ _
    _ ≤ _ := by linarith

lemma continuous_Qint {h : ℝ} (_h1 : h < 1) :
    Continuous (fun s => W h (s * h) * ((1 + Real.exp (-s)) ^ h - 1)) := by
  unfold W
  have hb : ∀ s : ℝ, 1 + Real.exp (-(s * h / (1 - h))) ≠ 0 := fun s => by positivity
  have hb2 : ∀ s : ℝ, 1 + Real.exp (-s) ≠ 0 := fun s => by positivity
  refine ((by fun_prop : Continuous fun s : ℝ => Real.exp (-(h * (s * h) / (1 - h)))).mul
    ((by fun_prop : Continuous fun s : ℝ => 1 + Real.exp (-(s * h / (1 - h)))).rpow_const
      (fun s => Or.inl (hb s)))).mul
    (((by fun_prop : Continuous fun s : ℝ => 1 + Real.exp (-s)).rpow_const
      (fun s => Or.inl (hb2 s))).sub continuous_const)

lemma intOn_Gq (h : ℝ) : IntegrableOn (Gq h) (Ioi 0) := by
  have := ((intOn_calL.const_mul (h / 2 - h ^ 2 * Real.log 2 / 2)).add
    (intOn_s_calL.const_mul (h ^ 2 / 4))).add (intOn_calL_sq.const_mul (h ^ 2 / 4))
  exact this

lemma integral_Gq (h : ℝ) : ∫ s in Ioi (0 : ℝ), Gq h s
    = π ^ 2 / 24 * h + (z3 / 4 - π ^ 2 * Real.log 2 / 24) * h ^ 2 := by
  unfold Gq
  have i1 : IntegrableOn (fun s => (h / 2 - h ^ 2 * Real.log 2 / 2) * calL s) (Ioi 0) :=
    intOn_calL.const_mul _
  have i2 : IntegrableOn (fun s => h ^ 2 / 4 * (s * calL s)) (Ioi 0) := intOn_s_calL.const_mul _
  have i3 : IntegrableOn (fun s => h ^ 2 / 4 * calL s ^ 2) (Ioi 0) := intOn_calL_sq.const_mul _
  have i12 : IntegrableOn (fun s => (h / 2 - h ^ 2 * Real.log 2 / 2) * calL s
      + h ^ 2 / 4 * (s * calL s)) (Ioi 0) := i1.add i2
  rw [integral_add i12 i3, integral_add i1 i2,
    integral_const_mul, integral_const_mul, integral_const_mul, J10, J11, J20]
  ring

theorem Q_expansion : ∃ C h₀ : ℝ, 0 < h₀ ∧ ∀ h : ℝ, 0 < h → h ≤ h₀ →
    |Q h - (Real.pi ^ 2 / 24 * h + (z3 / 4 - Real.pi ^ 2 * Real.log 2 / 24) * h ^ 2)|
      ≤ C * h ^ 3 := by
  refine ⟨66 * ∫ s in Ioi (0 : ℝ), (1 + s) ^ 3 * Real.exp (-s), 1 / 4, by norm_num,
    fun h hh hh4 => ?_⟩
  set F := fun s => W h (s * h) * ((1 + Real.exp (-s)) ^ h - 1) with hF
  have hFG : IntegrableOn (fun s => F s - Gq h s) (Ioi 0) := by
    refine intOn_of_bound ((continuous_Qint (by linarith)).sub ?_)
      (intOn_cube.const_mul (66 * h ^ 3)) (fun s hs => Q_pointwise hh hh4 hs.le)
    unfold Gq
    exact ((continuous_const.mul continuous_calL).add
      (continuous_const.mul (continuous_id.mul continuous_calL))).add
      (continuous_const.mul (continuous_calL.pow 2))
  have hQ : Q h = (∫ s in Ioi (0 : ℝ), (F s - Gq h s)) + ∫ s in Ioi (0 : ℝ), Gq h s := by
    rw [← integral_add hFG (intOn_Gq h)]
    unfold Q
    congr 1; funext s; simp only [hF]; ring
  rw [hQ, integral_Gq, add_sub_cancel_right]
  have hb : ‖∫ s in Ioi (0 : ℝ), (F s - Gq h s)‖
      ≤ ∫ s in Ioi (0 : ℝ), 66 * h ^ 3 * ((1 + s) ^ 3 * Real.exp (-s)) := by
    refine norm_integral_le_of_norm_le (intOn_cube.const_mul (66 * h ^ 3)) ?_
    exact (ae_restrict_iff' measurableSet_Ioi).2 (Filter.Eventually.of_forall
      (fun s hs => by rw [Real.norm_eq_abs]; exact Q_pointwise hh hh4 (le_of_lt hs)))
  rw [integral_const_mul, Real.norm_eq_abs] at hb
  calc _ ≤ _ := hb
    _ = _ := by ring



lemma E_eq_Q {p : ℝ} (hp : 1 < p) : E p = 4 / (p * (p - 1)) * Q (1 / p) := by
  have hp0 : 0 < p := by linarith
  have hp1 : 0 < p - 1 := by linarith
  have hpt : ∀ s, EIntegrand p s = W (1 / p) (s * (1 / p)) * ((1 + Real.exp (-s)) ^ (1 / p) - 1) :=
    fun s => by
      unfold EIntegrand W
      have a1 : 1 / p * (s * (1 / p)) / (1 - 1 / p) = s / (p * (p - 1)) := by
        field_simp
      have a2 : s * (1 / p) / (1 - 1 / p) = s / (p - 1) := by
        field_simp
      rw [a1, a2]
  unfold E Q
  rw [setIntegral_congr_fun measurableSet_Ioi (fun s _ => hpt s)]

theorem E_expansion : ∃ C p₀ : ℝ, 1 < p₀ ∧ ∀ p : ℝ, p₀ ≤ p →
    |E p - (Real.pi ^ 2 / 6 / p ^ 3
        + (Real.pi ^ 2 / 6 + z3 - Real.pi ^ 2 * Real.log 2 / 6) / p ^ 4)| ≤ C / p ^ 5 := by
  obtain ⟨C, h₀, hh₀, hQ⟩ := Q_expansion
  set q1 : ℝ := Real.pi ^ 2 / 24
  set q2 : ℝ := z3 / 4 - Real.pi ^ 2 * Real.log 2 / 24
  refine ⟨8 * (|C| + |q1 + q2|), max 4 (1 / h₀), lt_of_lt_of_le (by norm_num) (le_max_left _ _),
    fun p hp => ?_⟩
  have hp4 : 4 ≤ p := le_trans (le_max_left _ _) hp
  have hp0 : 0 < p := by linarith
  have hp1 : 0 < p - 1 := by linarith
  have hph : 1 / p ≤ h₀ := by
    rw [div_le_iff₀ hp0]
    have := le_trans (le_max_right _ _) hp
    rw [div_le_iff₀ hh₀] at this; linarith
  have hR := hQ (1 / p) (by positivity) hph
  rw [E_eq_Q (by linarith)]
  set R := Q (1 / p) - (q1 * (1 / p) + q2 * (1 / p) ^ 2) with hR_def
  have key : 4 / (p * (p - 1)) * Q (1 / p) - (Real.pi ^ 2 / 6 / p ^ 3
        + (Real.pi ^ 2 / 6 + z3 - Real.pi ^ 2 * Real.log 2 / 6) / p ^ 4)
      = 4 / (p * (p - 1)) * R + 4 * (q1 + q2) / (p ^ 4 * (p - 1)) := by
    rw [hR_def]
    field_simp
    ring
  rw [key]
  have hpp : 4 / (p * (p - 1)) ≤ 8 / p ^ 2 := by
    rw [div_le_div_iff₀ (by positivity) (by positivity)]; nlinarith
  have hpp' : 4 / (p ^ 4 * (p - 1)) ≤ 8 / p ^ 5 := by
    rw [div_le_div_iff₀ (by positivity) (by positivity)]
    have : 0 < p ^ 4 := by positivity
    nlinarith
  have t1 : |4 / (p * (p - 1)) * R| ≤ 8 * |C| / p ^ 5 := by
    rw [abs_mul, abs_of_pos (by positivity)]
    have hR' : |R| ≤ |C| * (1 / p) ^ 3 :=
      le_trans hR (mul_le_mul_of_nonneg_right (le_abs_self C) (by positivity))
    calc _ ≤ 8 / p ^ 2 * (|C| * (1 / p) ^ 3) :=
          mul_le_mul hpp hR' (abs_nonneg _) (by positivity)
      _ = _ := by field_simp
  have t2 : |4 * (q1 + q2) / (p ^ 4 * (p - 1))| ≤ 8 * |q1 + q2| / p ^ 5 := by
    rw [abs_div, abs_mul, abs_of_pos (by positivity : (0:ℝ) < p ^ 4 * (p - 1)), abs_of_pos
      (by norm_num : (0:ℝ) < 4)]
    calc 4 * |q1 + q2| / (p ^ 4 * (p - 1)) = |q1 + q2| * (4 / (p ^ 4 * (p - 1))) := by ring
      _ ≤ |q1 + q2| * (8 / p ^ 5) := mul_le_mul_of_nonneg_left hpp' (abs_nonneg _)
      _ = _ := by ring
  calc _ ≤ _ := abs_add_le _ _
    _ ≤ 8 * |C| / p ^ 5 + 8 * |q1 + q2| / p ^ 5 := add_le_add t1 t2
    _ = _ := by ring

theorem piP_expansion : ∃ C p₀ : ℝ, 1 < p₀ ∧ ∀ p : ℝ, p₀ ≤ p →
    |piP p - 4 * (1 + N0 / p + N1 / p ^ 2 + N2 / p ^ 3 + N3 / p ^ 4)| ≤ C / p ^ 5 := by
  obtain ⟨C, p₁, hp₁, hE⟩ := E_expansion
  refine ⟨C + 1, p₁, hp₁, fun p hp => ?_⟩
  have hp1 : 1 < p := lt_of_lt_of_le hp₁ hp
  have hp0 : 0 < p := by linarith
  rw [(split hp1).2.2]
  have hlog2 : Real.log 2 < 1 := by
    have := Real.log_two_lt_d9; linarith
  have hlog2' : 0 < Real.log 2 := Real.log_pos (by norm_num)
  set x := Real.log 2 * (-(1 / p)) with hx
  have h2 : (2 : ℝ) ^ (-(1 / p)) = Real.exp x := Real.rpow_def_of_pos (by norm_num) _
  have hxa : |x| = Real.log 2 / p := by
    rw [hx, abs_mul, abs_neg, abs_of_pos hlog2', abs_of_pos (by positivity)]; ring
  have hx1 : |x| ≤ 1 / p := by
    rw [hxa]; exact div_le_div_of_nonneg_right hlog2.le hp0.le
  have hx1' : |x| ≤ 1 := le_trans hx1 (by rw [div_le_one hp0]; linarith)
  have hT := Real.exp_bound hx1' (n := 5) (by norm_num)
  set S := ∑ m ∈ Finset.range 5, x ^ m / (m.factorial : ℝ) with hS
  have key : 4 * Real.exp x + E p - 4 * (1 + N0 / p + N1 / p ^ 2 + N2 / p ^ 3 + N3 / p ^ 4)
      = 4 * (Real.exp x - S) + (E p - (Real.pi ^ 2 / 6 / p ^ 3
        + (Real.pi ^ 2 / 6 + z3 - Real.pi ^ 2 * Real.log 2 / 6) / p ^ 4)) := by
    rw [hS, hx]
    simp only [Finset.sum_range_succ, Finset.sum_range_zero, Nat.factorial]
    unfold N0 N1 N2 N3
    push_cast
    field_simp
    ring
  rw [h2, key]
  have hx5 : |x| ^ 5 ≤ (1 / p) ^ 5 := pow_le_pow_left₀ (abs_nonneg _) hx1 5
  have t1 : |4 * (Real.exp x - S)| ≤ 1 / p ^ 5 := by
    rw [abs_mul, abs_of_pos (by norm_num : (0:ℝ) < 4)]
    have hc : ((Nat.succ 5 : ℕ) / ((Nat.factorial 5 : ℕ) * (5 : ℕ)) : ℝ) = 1 / 100 := by
      norm_num [Nat.factorial]
    rw [hc] at hT
    rw [div_pow, one_pow] at hx5
    have : 0 < 1 / p ^ 5 := by positivity
    have : |Real.exp x - S| ≤ 1 / p ^ 5 * (1 / 100) :=
      le_trans hT (mul_le_mul_of_nonneg_right hx5 (by norm_num))
    linarith
  have t2 := hE p hp
  calc _ ≤ _ := abs_add_le _ _
    _ ≤ 1 / p ^ 5 + C / p ^ 5 := add_le_add t1 t2
    _ = _ := by ring

end X02

theorem X02.proof_1 : X02.stmt_1 := fun _ hp => X02.split hp
theorem X02.proof_2 : X02.stmt_2 := ⟨X02.intOn_calL, X02.J10⟩
theorem X02.proof_3 : X02.stmt_3 := ⟨X02.intOn_s_calL, X02.J11⟩
theorem X02.proof_4 : X02.stmt_4 := ⟨X02.intOn_calL_sq, X02.J20⟩
theorem X02.proof_5 : X02.stmt_5 := X02.W_stmt5
theorem X02.proof_6 : X02.stmt_6 := X02.Q_expansion
theorem X02.proof_7 : X02.stmt_7 := X02.E_expansion
theorem X02.proof_8 : X02.stmt_8 := X02.piP_expansion
