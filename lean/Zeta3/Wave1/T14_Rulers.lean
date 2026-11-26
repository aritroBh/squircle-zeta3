import Mathlib
/-!
# T14 -- two exact rulers: `L_1(p) = 8` and `L_∞(p) = 8 * 2^(-1/p)`

Target id: T14.
Paper label: `prop:ruler` (paper/ieee/main.tex, lines 1384-1390), only the two exact values
`L_1(p)=8` and `L_\infty(p)=8\cdot2^{-1/p}`.
Supporting setup: paper/ieee/proofs/sec_ruler.tex, lines 5-12 (setting and statement),
19-25 (Lemma `lem:ruler:formula`: the length formula), 28-41 (Lemma `lem:ruler:ends`).

Quoted verbatim (main.tex 1384):

  Let's give ourselves a free ruler: let
  $L_q(p)=\frac4p\int_0^1\bigl(t^{q/p-q}+(1-t)^{q/p-q}\bigr)^{1/q}\dd t$ be the length of the
  $\ell^p$ circle measured with the $\ell^q$ norm, so that $L_2=P$ and $L_p=2\pi_p$.

Quoted verbatim (main.tex 1386-1388):

  \begin{proposition}[Two exact rulers$^\dagger$]\label{prop:ruler}
  For every $p\ge1$, $L_1(p)=8$ and $L_\infty(p)=8\cdot2^{-1/p}$. For every fixed $q\ge1$, the
  expansion of $L_q(p)$ in powers of $1/p$ converges for every $p>1$.
  \end{proposition}

Quoted verbatim (sec_ruler.tex 19-21, Lemma `lem:ruler:formula`):

  Let $\gamma(t)=(t^{1/p},(1-t)^{1/p})$, $0\le t\le1$, parametrize the arc of
  $|x|^p+|y|^p=1$ in the first quadrant. The length of the whole curve, measured with the norm
  $\|(u,v)\|_q=(|u|^q+|v|^q)^{1/q}$ ($1\le q<\infty$) or $\|(u,v)\|_\infty=\max(|u|,|v|)$, is
  $4\int_0^1\|\gamma'(t)\|\,dt$. ... For $q=\infty$ it equals
  $\frac4p\int_0^1\max\bigl(t^{1/p-1},(1-t)^{1/p-1}\bigr)dt$.

Quoted verbatim (sec_ruler.tex 31, Lemma `lem:ruler:ends`):

  For every $p\ge1$: $L_1(p)=8$, and $L_\infty(p)=\lim_{q\to\infty}L_q(p)=8\cdot2^{-1/p}$.

## Formalization

`loop p : ℝ → ℝ × ℝ` is a closed parametrization of the whole curve `|x|^p+|y|^p=1` on `[0,4]`,
made of the paper's quarter arc `γ(t) = (t^{1/p}, (1-t)^{1/p})` on `[0,1]` and its three images
under the sign changes of the coordinates (one per quadrant, glued end to end, clockwise).
Length is measured honestly as total variation (`eVariationOn`, the supremum over partitions of
`∑ ‖loop t_i - loop t_{i-1}‖`):
* in the `ℓ¹` norm, via `WithLp 1 (ℝ × ℝ)` (whose norm is `|u| + |v|`);
* in the `ℓ^∞` norm, via `ℝ × ℝ` itself (Mathlib's product norm/distance is the max).

* stmt_1: `loop p` is a continuous closed curve on `[0,4]`, injective on `[0,4)`, whose image
  is exactly the curve `{(x,y) | |x|^p + |y|^p = 1}` (so its variation is the length of the curve).
* stmt_2: the `ℓ¹` length (total variation in `WithLp 1 (ℝ × ℝ)`) of the curve is `8`.
* stmt_3: the `ℓ^∞` length (total variation in `ℝ × ℝ`, sup norm) is `8 * 2^(-1/p)`.
* stmt_4: the paper's formula `L_q(p)` (defined literally as `Lq` below) gives `L_1(p) = 8`.
* stmt_5: the supplement's integral for the sup ruler,
  `(4/p) ∫_0^1 max(t^{1/p-1}, (1-t)^{1/p-1}) dt = 8 * 2^(-1/p)`.
* stmt_6: `L_q(p) → 8 * 2^(-1/p)` as `q → ∞` (the supplement's identification of `L_∞` as the
  limit of the `L_q`, which removes any ambiguity in reading `L_∞` off the `L_q` formula).

All statements are for real `p ≥ 1`, as in the paper. Real powers are `Real.rpow`; on `[0,1]` the
bases are nonnegative, and the integrands are integrable for `p ≥ 1, q ≥ 1` (they blow up only like
`t^{1/p-1}`), so no statement relies on the "non-integrable integral = 0" convention. The value
`0 ^ (negative) = 0` of `Real.rpow` only affects the endpoints, a null set.

## NOT COVERED

* The third sentence of `prop:ruler` ("For every fixed $q\ge1$, the expansion of $L_q(p)$ in
  powers of $1/p$ converges for every $p>1$") -- outside this target.
* The general identity "`ℓ^q` length of the curve = `L_q(p)`" for `1 < q < ∞`
  (Lemma `lem:ruler:formula` for general `q`); only the `q = 1` and `q = ∞` rulers are stated,
  both as geometric lengths (stmt_2, stmt_3) and as the paper's integrals (stmt_4, stmt_5).
-/
-- BEGIN STATEMENTS
namespace T14

open MeasureTheory Set Filter

/-- The whole curve `|x|^p + |y|^p = 1`, traversed once on `[0,4]`: the paper's quarter arc
`(t^{1/p}, (1-t)^{1/p})` on `[0,1]`, followed by its reflections in the other three quadrants. -/
noncomputable def loop (p : ℝ) (t : ℝ) : ℝ × ℝ :=
  if t ≤ 1 then (t ^ (1 / p), (1 - t) ^ (1 / p))
  else if t ≤ 2 then ((2 - t) ^ (1 / p), -((t - 1) ^ (1 / p)))
  else if t ≤ 3 then (-((t - 2) ^ (1 / p)), -((3 - t) ^ (1 / p)))
  else (-((4 - t) ^ (1 / p)), (t - 3) ^ (1 / p))

/-- Length of the curve measured with the `ℓ¹` norm `|u| + |v|`: total variation of `loop p`
viewed in `WithLp 1 (ℝ × ℝ)`. -/
noncomputable def l1Length (p : ℝ) : ENNReal :=
  eVariationOn (fun t => WithLp.toLp 1 (loop p t)) (Icc (0 : ℝ) 4)

/-- Length of the curve measured with the `ℓ^∞` norm `max(|u|, |v|)`: total variation of
`loop p` in `ℝ × ℝ`, whose Mathlib distance is the max of the coordinate distances. -/
noncomputable def supLength (p : ℝ) : ENNReal :=
  eVariationOn (loop p) (Icc (0 : ℝ) 4)

/-- The paper's `L_q(p)=\frac4p\int_0^1(t^{q/p-q}+(1-t)^{q/p-q})^{1/q}\,dt` (main.tex 1384). -/
noncomputable def Lq (p q : ℝ) : ℝ :=
  (4 / p) * ∫ t in (0 : ℝ)..1, (t ^ (q / p - q) + (1 - t) ^ (q / p - q)) ^ (1 / q)

/-- `loop p` is a simple closed parametrization of the curve `|x|^p+|y|^p=1`
(sec_ruler.tex 19-23: the four quadrant arcs are images of each other under sign changes). -/
def stmt_1 : Prop :=
  ∀ p : ℝ, 1 ≤ p →
    ContinuousOn (loop p) (Icc 0 4) ∧ loop p 4 = loop p 0 ∧ InjOn (loop p) (Ico 0 4) ∧
    loop p '' Icc 0 4 = {z : ℝ × ℝ | |z.1| ^ p + |z.2| ^ p = 1}

/-- "For every $p\ge1$, $L_1(p)=8$", with `L_1` the `ℓ¹` length of the curve. -/
def stmt_2 : Prop :=
  ∀ p : ℝ, 1 ≤ p → l1Length p = ENNReal.ofReal 8

/-- "For every $p\ge1$, ... $L_\infty(p)=8\cdot2^{-1/p}$", with `L_∞` the sup-norm length. -/
def stmt_3 : Prop :=
  ∀ p : ℝ, 1 ≤ p → supLength p = ENNReal.ofReal (8 * (2 : ℝ) ^ (-1 / p))

/-- "$L_1(p)=8$" for the paper's integral formula `L_q(p)` at `q = 1`. -/
def stmt_4 : Prop :=
  ∀ p : ℝ, 1 ≤ p → Lq p 1 = 8

/-- sec_ruler.tex 21 and 37-39: `\frac4p\int_0^1\max(t^{1/p-1},(1-t)^{1/p-1})\,dt = 8\cdot2^{-1/p}`. -/
def stmt_5 : Prop :=
  ∀ p : ℝ, 1 ≤ p →
    (4 / p) * ∫ t in (0 : ℝ)..1, max (t ^ (1 / p - 1)) ((1 - t) ^ (1 / p - 1))
      = 8 * (2 : ℝ) ^ (-1 / p)

/-- sec_ruler.tex 31: "$L_\infty(p)=\lim_{q\to\infty}L_q(p)=8\cdot2^{-1/p}$". -/
def stmt_6 : Prop :=
  ∀ p : ℝ, 1 ≤ p → Tendsto (fun q => Lq p q) atTop (nhds (8 * (2 : ℝ) ^ (-1 / p)))

end T14
-- END STATEMENTS
open MeasureTheory Set Filter

theorem T14.var_congr {E F : Type*} [PseudoEMetricSpace E] [PseudoEMetricSpace F]
    (f : ℝ → E) (h : ℝ → F) (s : Set ℝ)
    (hd : ∀ a ∈ s, ∀ b ∈ s, edist (f a) (f b) = edist (h a) (h b)) :
    eVariationOn f s = eVariationOn h s := by
  unfold eVariationOn
  refine iSup_congr fun p => Finset.sum_congr rfl fun i _ => ?_
  exact hd _ (p.2.2.2 _) _ (p.2.2.2 _)

theorem T14.var_shift {E F : Type*} [PseudoEMetricSpace E] [PseudoEMetricSpace F]
    (f : ℝ → E) (g : ℝ → F) (k : ℝ)
    (hd : ∀ a ∈ Icc k (k+1), ∀ b ∈ Icc k (k+1), edist (f a) (f b) = edist (g (a-k)) (g (b-k))) :
    eVariationOn f (Icc k (k+1)) = eVariationOn g (Icc 0 1) := by
  rw [T14.var_congr f (g ∘ fun t => t - k) _ hd,
    eVariationOn.comp_eq_of_monotoneOn g (fun t => t - k) (fun a _ b _ h => by simpa using h)]
  congr 1
  ext x
  simp only [mem_image, mem_Icc]
  constructor
  · rintro ⟨y, ⟨h1, h2⟩, rfl⟩; constructor <;> linarith
  · rintro ⟨h1, h2⟩; exact ⟨x + k, ⟨by linarith, by linarith⟩, by ring⟩

theorem T14.var_mono (h : ℝ → ℝ) {a b : ℝ} (hab : a ≤ b) (hm : MonotoneOn h (Icc a b)) :
    eVariationOn h (Icc a b) = ENNReal.ofReal (h b - h a) := by
  have := hm.eVariationOn_eq (a := a) (b := b) ⟨le_rfl, hab⟩ ⟨hab, le_rfl⟩
  simpa using this

theorem T14.l1_edist (z w : ℝ × ℝ) :
    edist (WithLp.toLp 1 z) (WithLp.toLp 1 w) = ENNReal.ofReal (|z.1 - w.1| + |z.2 - w.2|) := by
  rw [WithLp.prod_edist_eq_add (by simp)]
  simp [edist_dist, Real.dist_eq, ENNReal.ofReal_add, abs_nonneg]

theorem T14.sup_edist (z w : ℝ × ℝ) :
    edist z w = ENNReal.ofReal (max |z.1 - w.1| |z.2 - w.2|) := by
  rw [Prod.edist_eq]
  simp [edist_dist, Real.dist_eq, ENNReal.ofReal_max]

theorem T14.core_conc {r a b : ℝ} (hr0 : 0 ≤ r) (hr1 : r ≤ 1) (ha : 0 ≤ a) (hab : a ≤ b)
    (hb : b ≤ 1/2) : (1 - a) ^ r - (1 - b) ^ r ≤ b ^ r - a ^ r := by
  rcases eq_or_lt_of_le (hab.trans hb |>.trans_eq rfl : a ≤ 1/2) with h | h
  · have : b = a := by linarith
    subst this; linarith
  · have hc := Real.concaveOn_rpow hr0 hr1
    set lam : ℝ := (1 - a - b) / (1 - 2 * a) with hlam
    have hd : 0 < 1 - 2 * a := by linarith
    have hd' : 1 - 2 * a ≠ 0 := hd.ne'
    have hl0 : 0 ≤ lam := div_nonneg (by linarith) hd.le
    have hl1 : lam ≤ 1 := by rw [hlam, div_le_one hd]; linarith
    have hm : lam * (1 - 2 * a) = 1 - a - b := by rw [hlam]; exact div_mul_cancel₀ _ hd'
    have e1 : lam * a + (1 - lam) * (1 - a) = b := by linear_combination (-1) * hm
    have e2 : (1 - lam) * a + lam * (1 - a) = 1 - b := by linear_combination hm
    have h1 := hc.2 (x := a) (y := 1 - a) (by simpa using ha) (by simp; linarith) hl0
      (by linarith : 0 ≤ 1 - lam) (by ring)
    have h2 := hc.2 (x := a) (y := 1 - a) (by simpa using ha) (by simp; linarith)
      (by linarith : 0 ≤ 1 - lam) hl0 (by ring)
    simp only [smul_eq_mul] at h1 h2
    rw [e1] at h1
    rw [e2] at h2
    nlinarith

theorem T14.var_split {E : Type*} [PseudoEMetricSpace E] (f : ℝ → E) {a b c : ℝ}
    (hab : a ≤ b) (hbc : b ≤ c) :
    eVariationOn f (Icc a b) + eVariationOn f (Icc b c) = eVariationOn f (Icc a c) := by
  have := eVariationOn.Icc_add_Icc f (s := Icc a c) hab hbc ⟨hab, hbc⟩
  rwa [Set.inter_eq_right.2 (Icc_subset_Icc le_rfl hbc),
    Set.inter_eq_right.2 (Icc_subset_Icc hab le_rfl), Set.inter_self] at this

noncomputable def T14.arc (p s : ℝ) : ℝ × ℝ := (s ^ (1 / p), (1 - s) ^ (1 / p))

theorem T14.edist_real (u v : ℝ) : edist u v = ENNReal.ofReal |u - v| := by
  rw [edist_dist, Real.dist_eq]

theorem T14.arc_l1 {p : ℝ} (hp : 1 ≤ p) :
    eVariationOn (fun t => WithLp.toLp 1 (T14.arc p t)) (Icc 0 1) = ENNReal.ofReal 2 := by
  have hr : 0 < 1 / p := by positivity
  rw [T14.var_congr _ (fun t => t ^ (1 / p) - (1 - t) ^ (1 / p)) _ ?_]
  · rw [T14.var_mono _ zero_le_one ?_]
    · show ENNReal.ofReal ((1:ℝ) ^ (1 / p) - (1 - 1) ^ (1 / p) - ((0:ℝ) ^ (1 / p) - (1 - 0) ^ (1 / p))) = _
      rw [sub_self, Real.zero_rpow hr.ne', Real.one_rpow, sub_zero]
      norm_num
    · intro a ha b hb hab
      have h1 : a ^ (1 / p) ≤ b ^ (1 / p) := Real.rpow_le_rpow ha.1 hab hr.le
      have h2 : (1 - b) ^ (1 / p) ≤ (1 - a) ^ (1 / p) :=
        Real.rpow_le_rpow (by linarith [hb.2]) (by linarith) hr.le
      simp only; linarith
  · intro a ha b hb
    rw [T14.l1_edist, T14.edist_real]
    congr 1
    simp only [T14.arc]
    rcases le_total a b with hab | hab
    · have h1 : a ^ (1 / p) ≤ b ^ (1 / p) := Real.rpow_le_rpow ha.1 hab hr.le
      have h2 : (1 - b) ^ (1 / p) ≤ (1 - a) ^ (1 / p) :=
        Real.rpow_le_rpow (by linarith [hb.2]) (by linarith) hr.le
      rw [abs_of_nonpos (by linarith : a ^ (1 / p) - b ^ (1 / p) ≤ 0),
        abs_of_nonneg (by linarith : (1 - a) ^ (1 / p) - (1 - b) ^ (1 / p) ≥ 0),
        abs_of_nonpos (by linarith)]
      ring
    · have h1 : b ^ (1 / p) ≤ a ^ (1 / p) := Real.rpow_le_rpow hb.1 hab hr.le
      have h2 : (1 - a) ^ (1 / p) ≤ (1 - b) ^ (1 / p) :=
        Real.rpow_le_rpow (by linarith [ha.2]) (by linarith) hr.le
      rw [abs_of_nonneg (by linarith : a ^ (1 / p) - b ^ (1 / p) ≥ 0),
        abs_of_nonpos (by linarith : (1 - a) ^ (1 / p) - (1 - b) ^ (1 / p) ≤ 0),
        abs_of_nonneg (by linarith)]
      ring

theorem T14.arc_sup {p : ℝ} (hp : 1 ≤ p) :
    eVariationOn (T14.arc p) (Icc 0 1) = ENNReal.ofReal (2 * (2 : ℝ) ^ (-1 / p)) := by
  have hr : 0 < 1 / p := by positivity
  have hr1 : 1 / p ≤ 1 := by rw [div_le_one (by linarith)]; exact hp
  have hc : (1 / 2 : ℝ) ^ (1 / p) = (2 : ℝ) ^ (-1 / p) := by
    rw [one_div, Real.inv_rpow (by norm_num), ← Real.rpow_neg (by norm_num)]
    congr 1; ring
  rw [← T14.var_split (T14.arc p) (b := 1 / 2) (by norm_num) (by norm_num)]
  have L : eVariationOn (T14.arc p) (Icc 0 (1 / 2)) = ENNReal.ofReal ((2 : ℝ) ^ (-1 / p)) := by
    rw [T14.var_congr _ (fun t => t ^ (1 / p)) _ ?_]
    · rw [T14.var_mono _ (by norm_num) ?_]
      · rw [Real.zero_rpow hr.ne', sub_zero, hc]
      · intro a ha b hb hab
        exact Real.rpow_le_rpow ha.1 hab hr.le
    · intro a ha b hb
      rw [T14.sup_edist, T14.edist_real]
      congr 1
      simp only [T14.arc]
      have key : ∀ a b : ℝ, 0 ≤ a → a ≤ b → b ≤ 1 / 2 →
          |(1 - a) ^ (1 / p) - (1 - b) ^ (1 / p)| ≤ |a ^ (1 / p) - b ^ (1 / p)| := by
        intro a b ha hab hb
        have := T14.core_conc hr.le hr1 ha hab hb
        have h1 : a ^ (1 / p) ≤ b ^ (1 / p) := Real.rpow_le_rpow ha hab hr.le
        have h2 : (1 - b) ^ (1 / p) ≤ (1 - a) ^ (1 / p) :=
          Real.rpow_le_rpow (by linarith) (by linarith) hr.le
        rw [abs_of_nonneg (by linarith), abs_of_nonpos (by linarith)]
        linarith
      rcases le_total a b with hab | hab
      · exact max_eq_left (key a b ha.1 hab hb.2)
      · rw [abs_sub_comm ((1 - a) ^ (1 / p)), abs_sub_comm (a ^ (1 / p))]
        exact max_eq_left (key b a hb.1 hab ha.2)
  have R : eVariationOn (T14.arc p) (Icc (1 / 2) 1) = ENNReal.ofReal ((2 : ℝ) ^ (-1 / p)) := by
    rw [T14.var_congr _ (fun t => -((1 - t) ^ (1 / p))) _ ?_]
    · rw [T14.var_mono _ (by norm_num) ?_]
      · show ENNReal.ofReal (-((1 - 1 : ℝ) ^ (1 / p)) - -((1 - 1 / 2 : ℝ) ^ (1 / p))) = _
        rw [sub_self, Real.zero_rpow hr.ne', show (1:ℝ) - 1/2 = 1/2 by norm_num, hc]
        congr 1; ring
      · intro a ha b hb hab
        have : (1 - b) ^ (1 / p) ≤ (1 - a) ^ (1 / p) :=
          Real.rpow_le_rpow (by linarith [hb.2]) (by linarith) hr.le
        simp only; linarith
    · intro a ha b hb
      rw [T14.sup_edist, T14.edist_real]
      congr 1
      simp only [T14.arc]
      have key : ∀ a b : ℝ, 1 / 2 ≤ a → a ≤ b → b ≤ 1 →
          |a ^ (1 / p) - b ^ (1 / p)| ≤ |(1 - a) ^ (1 / p) - (1 - b) ^ (1 / p)| := by
        intro a b ha hab hb
        have := T14.core_conc (a := 1 - b) (b := 1 - a) hr.le hr1 (by linarith) (by linarith)
          (by linarith)
        simp only [sub_sub_cancel] at this
        have h1 : a ^ (1 / p) ≤ b ^ (1 / p) := Real.rpow_le_rpow (by linarith) hab hr.le
        have h2 : (1 - b) ^ (1 / p) ≤ (1 - a) ^ (1 / p) :=
          Real.rpow_le_rpow (by linarith) (by linarith) hr.le
        rw [abs_of_nonpos (by linarith), abs_of_nonneg (by linarith)]
        linarith
      rw [show |-(1 - a) ^ (1 / p) - -(1 - b) ^ (1 / p)| = |(1 - a) ^ (1 / p) - (1 - b) ^ (1 / p)|
        by rw [← abs_neg]; congr 1; ring]
      rcases le_total a b with hab | hab
      · exact (max_eq_right (key a b ha.1 hab hb.2))
      · rw [abs_sub_comm (a ^ (1 / p)), abs_sub_comm ((1 - a) ^ (1 / p))]
        exact (max_eq_right (key b a hb.1 hab ha.2))
  rw [L, R, ← ENNReal.ofReal_add (by positivity) (by positivity)]
  congr 1; ring

theorem T14.loop_Q1 {p : ℝ} (t : ℝ) (ht : t ∈ Icc (0:ℝ) 1) :
    T14.loop p t = T14.arc p (t - 0) := by
  simp [T14.loop, T14.arc, ht.2]

theorem T14.loop_Q2 {p : ℝ} (hp : 1 ≤ p) (t : ℝ) (ht : t ∈ Icc (1:ℝ) 2) :
    T14.loop p t = ((T14.arc p (t - 1)).2, -(T14.arc p (t - 1)).1) := by
  have hr : (1 / p) ≠ 0 := by positivity
  have h0 : (0:ℝ) ^ p⁻¹ = 0 := by simpa using Real.zero_rpow hr
  rcases ht.1.eq_or_lt with rfl | h1
  · simp [T14.loop, T14.arc, h0]
  · have h : ¬ t ≤ 1 := not_le.2 h1
    simp only [T14.loop, T14.arc, h, ite_false, ht.2, ite_true]
    rw [show 1 - (t - 1) = 2 - t by ring]

theorem T14.loop_Q3 {p : ℝ} (hp : 1 ≤ p) (t : ℝ) (ht : t ∈ Icc (2:ℝ) 3) :
    T14.loop p t = (-(T14.arc p (t - 2)).1, -(T14.arc p (t - 2)).2) := by
  have hr : (1 / p) ≠ 0 := by positivity
  have h0 : (0:ℝ) ^ p⁻¹ = 0 := by simpa using Real.zero_rpow hr
  have h : ¬ t ≤ 1 := by have := ht.1; intro; linarith
  rcases ht.1.eq_or_lt with rfl | h1
  · norm_num [T14.loop, T14.arc, h0]
  · have h' : ¬ t ≤ 2 := not_le.2 h1
    simp only [T14.loop, T14.arc, h, h', ite_false, ht.2, ite_true]
    rw [show 1 - (t - 2) = 3 - t by ring]

theorem T14.loop_Q4 {p : ℝ} (hp : 1 ≤ p) (t : ℝ) (ht : t ∈ Icc (3:ℝ) 4) :
    T14.loop p t = (-(T14.arc p (t - 3)).2, (T14.arc p (t - 3)).1) := by
  have hr : (1 / p) ≠ 0 := by positivity
  have h0 : (0:ℝ) ^ p⁻¹ = 0 := by simpa using Real.zero_rpow hr
  have h : ¬ t ≤ 1 := by have := ht.1; intro; linarith
  have h' : ¬ t ≤ 2 := by have := ht.1; intro; linarith
  rcases ht.1.eq_or_lt with rfl | h1
  · norm_num [T14.loop, T14.arc, h0]
  · have h'' : ¬ t ≤ 3 := not_le.2 h1
    simp only [T14.loop, T14.arc, h, h', h'', ite_false]
    rw [show 1 - (t - 3) = 4 - t by ring]

theorem T14.abs_neg_sub_neg (u v : ℝ) : |-u - -v| = |u - v| := by
  rw [neg_sub_neg, abs_sub_comm]

theorem T14.piece_sup {p : ℝ} (k m : ℝ) (hm : m = k + 1) (A : ℝ × ℝ → ℝ × ℝ)
    (hA : ∀ z w, edist (A z) (A w) = edist z w)
    (hf : ∀ t ∈ Icc k m, T14.loop p t = A (T14.arc p (t - k))) :
    eVariationOn (T14.loop p) (Icc k m) = eVariationOn (T14.arc p) (Icc 0 1) := by
  subst hm
  refine T14.var_shift _ _ k fun a ha b hb => ?_
  rw [hf a ha, hf b hb, hA]

theorem T14.piece_l1 {p : ℝ} (k m : ℝ) (hm : m = k + 1) (A : ℝ × ℝ → ℝ × ℝ)
    (hA : ∀ z w, edist (WithLp.toLp 1 (A z)) (WithLp.toLp 1 (A w)) =
      edist (WithLp.toLp 1 z) (WithLp.toLp 1 w))
    (hf : ∀ t ∈ Icc k m, T14.loop p t = A (T14.arc p (t - k))) :
    eVariationOn (fun t => WithLp.toLp 1 (T14.loop p t)) (Icc k m) =
      eVariationOn (fun t => WithLp.toLp 1 (T14.arc p t)) (Icc 0 1) := by
  subst hm
  refine T14.var_shift _ _ k fun a ha b hb => ?_
  simp only [hf a ha, hf b hb, hA]

theorem T14.loop_split {E : Type*} [PseudoEMetricSpace E] (f : ℝ → E) :
    eVariationOn f (Icc 0 4) = eVariationOn f (Icc 0 1) + eVariationOn f (Icc 1 2) +
      eVariationOn f (Icc 2 3) + eVariationOn f (Icc 3 4) := by
  rw [← T14.var_split f (b := 1) (by norm_num) (by norm_num),
    ← T14.var_split f (a := 1) (b := 2) (by norm_num) (by norm_num),
    ← T14.var_split f (a := 2) (b := 3) (by norm_num) (by norm_num)]
  ring

theorem T14.proof_2' : T14.stmt_2 := by
  intro p hp
  unfold T14.l1Length
  rw [T14.loop_split]
  rw [T14.piece_l1 0 1 (by norm_num) (fun z => z) (fun z w => rfl) (fun t ht => T14.loop_Q1 t ht),
    T14.piece_l1 1 2 (by norm_num) (fun z => (z.2, -z.1))
      (fun z w => by
        rw [T14.l1_edist, T14.l1_edist]
        simp only [T14.abs_neg_sub_neg]
        congr 1; ring)
      (fun t ht => T14.loop_Q2 hp t ht),
    T14.piece_l1 2 3 (by norm_num) (fun z => (-z.1, -z.2))
      (fun z w => by
        rw [T14.l1_edist, T14.l1_edist]
        simp only [T14.abs_neg_sub_neg])
      (fun t ht => T14.loop_Q3 hp t ht),
    T14.piece_l1 3 4 (by norm_num) (fun z => (-z.2, z.1))
      (fun z w => by
        rw [T14.l1_edist, T14.l1_edist]
        simp only [T14.abs_neg_sub_neg]
        congr 1; ring)
      (fun t ht => T14.loop_Q4 hp t ht),
    T14.arc_l1 hp]
  rw [← ENNReal.ofReal_add (by norm_num) (by norm_num), ← ENNReal.ofReal_add (by norm_num) (by norm_num),
    ← ENNReal.ofReal_add (by norm_num) (by norm_num)]
  norm_num

theorem T14.proof_3' : T14.stmt_3 := by
  intro p hp
  unfold T14.supLength
  rw [T14.loop_split]
  rw [T14.piece_sup 0 1 (by norm_num) (fun z => z) (fun z w => rfl) (fun t ht => T14.loop_Q1 t ht),
    T14.piece_sup 1 2 (by norm_num) (fun z => (z.2, -z.1))
      (fun z w => by
        rw [T14.sup_edist, T14.sup_edist]
        simp only [T14.abs_neg_sub_neg]
        rw [max_comm])
      (fun t ht => T14.loop_Q2 hp t ht),
    T14.piece_sup 2 3 (by norm_num) (fun z => (-z.1, -z.2))
      (fun z w => by
        rw [T14.sup_edist, T14.sup_edist]
        simp only [T14.abs_neg_sub_neg])
      (fun t ht => T14.loop_Q3 hp t ht),
    T14.piece_sup 3 4 (by norm_num) (fun z => (-z.2, z.1))
      (fun z w => by
        rw [T14.sup_edist, T14.sup_edist]
        simp only [T14.abs_neg_sub_neg]
        rw [max_comm])
      (fun t ht => T14.loop_Q4 hp t ht),
    T14.arc_sup hp]
  have hc : (0:ℝ) ≤ 2 * (2 : ℝ) ^ (-1 / p) := by positivity
  rw [← ENNReal.ofReal_add hc hc, ← ENNReal.ofReal_add (by positivity) hc,
    ← ENNReal.ofReal_add (by positivity) hc]
  congr 1; ring

theorem T14.pw1 {p : ℝ} (hp : 1 ≤ p) {x : ℝ} (hx : 0 ≤ x) : (x ^ (1 / p)) ^ p = x := by
  rw [← Real.rpow_mul hx, one_div_mul_cancel (by linarith), Real.rpow_one]

theorem T14.pw2 {p : ℝ} (hp : 1 ≤ p) {x : ℝ} (hx : 0 ≤ x) : (x ^ p) ^ (1 / p) = x := by
  rw [← Real.rpow_mul hx, mul_one_div_cancel (by linarith), Real.rpow_one]

theorem T14.arc_cont {p : ℝ} (hp : 1 ≤ p) : Continuous (T14.arc p) := by
  have hr : 0 ≤ 1 / p := by positivity
  exact ((Real.continuous_rpow_const hr).comp continuous_id).prodMk
    ((Real.continuous_rpow_const hr).comp (continuous_const.sub continuous_id))

theorem T14.loop_cont {p : ℝ} (hp : 1 ≤ p) : ContinuousOn (T14.loop p) (Icc 0 4) := by
  have hc := T14.arc_cont hp
  have c1 : ContinuousOn (T14.loop p) (Icc 0 1) :=
    (hc.comp (continuous_id.sub continuous_const) |>.continuousOn).congr
      (fun t ht => T14.loop_Q1 t ht)
  have c2 : ContinuousOn (T14.loop p) (Icc 1 2) :=
    ((hc.comp (continuous_id.sub continuous_const)).snd.prodMk
      (hc.comp (continuous_id.sub continuous_const)).fst.neg |>.continuousOn).congr
      (fun t ht => T14.loop_Q2 hp t ht)
  have c3 : ContinuousOn (T14.loop p) (Icc 2 3) :=
    ((hc.comp (continuous_id.sub continuous_const)).fst.neg.prodMk
      (hc.comp (continuous_id.sub continuous_const)).snd.neg |>.continuousOn).congr
      (fun t ht => T14.loop_Q3 hp t ht)
  have c4 : ContinuousOn (T14.loop p) (Icc 3 4) :=
    ((hc.comp (continuous_id.sub continuous_const)).snd.neg.prodMk
      (hc.comp (continuous_id.sub continuous_const)).fst |>.continuousOn).congr
      (fun t ht => T14.loop_Q4 hp t ht)
  rw [← Icc_union_Icc_eq_Icc (a := 0) (b := 1) (c := 4) (by norm_num) (by norm_num),
    ← Icc_union_Icc_eq_Icc (a := 1) (b := 2) (c := 4) (by norm_num) (by norm_num),
    ← Icc_union_Icc_eq_Icc (a := 2) (b := 3) (c := 4) (by norm_num) (by norm_num)]
  exact c1.union_of_isClosed (c2.union_of_isClosed (c3.union_of_isClosed c4 isClosed_Icc isClosed_Icc)
    isClosed_Icc (isClosed_Icc.union isClosed_Icc)) isClosed_Icc
    (isClosed_Icc.union (isClosed_Icc.union isClosed_Icc))

noncomputable def T14.ginv (p : ℝ) (z : ℝ × ℝ) : ℝ :=
  if 0 ≤ z.1 ∧ 0 < z.2 then |z.1| ^ p
  else if 0 < z.1 ∧ z.2 ≤ 0 then 1 + |z.2| ^ p
  else if z.1 ≤ 0 ∧ z.2 < 0 then 2 + |z.1| ^ p
  else 3 + |z.2| ^ p

theorem T14.ginv_loop {p : ℝ} (hp : 1 ≤ p) (t : ℝ) (ht : t ∈ Ico (0:ℝ) 4) :
    T14.ginv p (T14.loop p t) = t := by
  have hr : 0 < 1 / p := by positivity
  have nn : ∀ x : ℝ, 0 ≤ x → 0 ≤ x ^ (1 / p) := fun x hx => Real.rpow_nonneg hx _
  have ps : ∀ x : ℝ, 0 < x → 0 < x ^ (1 / p) := fun x hx => Real.rpow_pos_of_pos hx _
  have ht2 := ht.2
  by_cases h1 : t < 1
  · have hz := T14.loop_Q1 (p := p) t ⟨ht.1, h1.le⟩
    rw [hz]
    simp only [T14.arc, sub_zero]
    have a1 := nn t ht.1
    have a2 := ps (1 - t) (by linarith)
    rw [T14.ginv, if_pos ⟨a1, a2⟩, abs_of_nonneg a1, T14.pw1 hp ht.1]
  by_cases h2 : t < 2
  · have hz := T14.loop_Q2 hp t ⟨not_lt.1 h1, h2.le⟩
    rw [hz]
    simp only [T14.arc]
    have a1 := nn (t - 1) (by linarith)
    have a2 := ps (1 - (t - 1)) (by linarith)
    rw [T14.ginv, if_neg (by intro h; linarith [h.2]), if_pos ⟨a2, by linarith⟩, abs_neg,
      abs_of_nonneg a1, T14.pw1 hp (by linarith)]
    ring
  by_cases h3 : t < 3
  · have hz := T14.loop_Q3 hp t ⟨not_lt.1 h2, h3.le⟩
    rw [hz]
    simp only [T14.arc]
    have a1 := nn (t - 2) (by linarith)
    have a2 := ps (1 - (t - 2)) (by linarith)
    rw [T14.ginv, if_neg (by intro h; linarith [h.2]), if_neg (by intro h; linarith [h.1]),
      if_pos ⟨by linarith, by linarith⟩, abs_neg, abs_of_nonneg a1, T14.pw1 hp (by linarith)]
    ring
  · have hz := T14.loop_Q4 hp t ⟨not_lt.1 h3, ht.2.le⟩
    rw [hz]
    simp only [T14.arc]
    have a1 := nn (t - 3) (by linarith)
    have a2 := ps (1 - (t - 3)) (by linarith)
    rw [T14.ginv, if_neg (by intro h; linarith [h.1]), if_neg (by intro h; linarith [h.1]),
      if_neg (by intro h; linarith [h.2]), abs_of_nonneg a1, T14.pw1 hp (by linarith)]
    ring

theorem T14.norm_loop {p : ℝ} (hp : 1 ≤ p) (t : ℝ) (ht : t ∈ Icc (0:ℝ) 4) :
    |(T14.loop p t).1| ^ p + |(T14.loop p t).2| ^ p = 1 := by
  have nn : ∀ x : ℝ, 0 ≤ x → 0 ≤ x ^ (1 / p) := fun x hx => Real.rpow_nonneg hx _
  have key : ∀ s : ℝ, 0 ≤ s → s ≤ 1 →
      |(T14.arc p s).1| ^ p + |(T14.arc p s).2| ^ p = 1 := by
    intro s h0 h1
    simp only [T14.arc]
    rw [abs_of_nonneg (nn _ h0), abs_of_nonneg (nn _ (by linarith)), T14.pw1 hp h0,
      T14.pw1 hp (by linarith)]
    ring
  by_cases h1 : t ≤ 1
  · rw [T14.loop_Q1 t ⟨ht.1, h1⟩]; exact key _ (by linarith [ht.1]) (by linarith)
  by_cases h2 : t ≤ 2
  · rw [T14.loop_Q2 hp t ⟨(not_le.1 h1).le, h2⟩]
    have := key (t - 1) (by linarith [not_le.1 h1]) (by linarith)
    simpa [abs_neg, add_comm] using this
  by_cases h3 : t ≤ 3
  · rw [T14.loop_Q3 hp t ⟨(not_le.1 h2).le, h3⟩]
    have := key (t - 2) (by linarith [not_le.1 h2]) (by linarith)
    simpa [abs_neg] using this
  · rw [T14.loop_Q4 hp t ⟨(not_le.1 h3).le, ht.2⟩]
    have := key (t - 3) (by linarith [not_le.1 h3]) (by linarith [ht.2])
    simpa [abs_neg, add_comm] using this

theorem T14.exists_t {p : ℝ} (hp : 1 ≤ p) (z : ℝ × ℝ) (hz : |z.1| ^ p + |z.2| ^ p = 1) :
    ∃ t ∈ Icc (0:ℝ) 4, T14.loop p t = z := by
  have hp0 : 0 < p := by linarith
  set u := |z.1| ^ p with hu
  set v := |z.2| ^ p with hv
  have u0 : 0 ≤ u := Real.rpow_nonneg (abs_nonneg _) _
  have v0 : 0 ≤ v := Real.rpow_nonneg (abs_nonneg _) _
  have ur : u ^ (1 / p) = |z.1| := T14.pw2 hp (abs_nonneg _)
  have vr : v ^ (1 / p) = |z.2| := T14.pw2 hp (abs_nonneg _)
  have u1 : (1 - u) = v := by linarith
  have v1 : (1 - v) = u := by linarith
  by_cases c1 : 0 ≤ z.1 ∧ 0 < z.2
  · refine ⟨u, ⟨u0, by linarith⟩, ?_⟩
    rw [T14.loop_Q1 u ⟨u0, by linarith⟩]
    simp only [T14.arc, sub_zero, u1, ur, vr]
    ext <;> simp [abs_of_nonneg c1.1, abs_of_pos c1.2]
  by_cases c2 : 0 < z.1 ∧ z.2 ≤ 0
  · refine ⟨1 + v, ⟨by linarith, by linarith⟩, ?_⟩
    rw [T14.loop_Q2 hp (1 + v) ⟨by linarith, by linarith⟩]
    simp only [T14.arc, add_sub_cancel_left, v1, ur, vr]
    ext <;> simp [abs_of_pos c2.1, abs_of_nonpos c2.2]
  by_cases c3 : z.1 ≤ 0 ∧ z.2 < 0
  · refine ⟨2 + u, ⟨by linarith, by linarith⟩, ?_⟩
    rw [T14.loop_Q3 hp (2 + u) ⟨by linarith, by linarith⟩]
    simp only [T14.arc, add_sub_cancel_left, u1, ur, vr]
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
    rw [T14.loop_Q4 hp (3 + v) ⟨by linarith, by linarith⟩]
    simp only [T14.arc, add_sub_cancel_left, v1, ur, vr]
    ext <;> simp [abs_of_nonpos h1, abs_of_nonneg h2]

theorem T14.proof_1' : T14.stmt_1 := by
  intro p hp
  have hr : (1 / p) ≠ 0 := by positivity
  have h0 : (0:ℝ) ^ p⁻¹ = 0 := by simpa using Real.zero_rpow hr
  refine ⟨T14.loop_cont hp, ?_, ?_, ?_⟩
  · norm_num [T14.loop, h0]
  · exact (Set.LeftInvOn.injOn (f₁' := T14.ginv p) (fun t ht => T14.ginv_loop hp t ht))
  · ext z
    constructor
    · rintro ⟨t, ht, rfl⟩
      exact T14.norm_loop hp t ht
    · intro hz
      exact T14.exists_t hp z hz

theorem T14.int_rpow {s : ℝ} (hs : -1 < s) :
    IntervalIntegrable (fun x : ℝ => x ^ s) volume 0 1 :=
  intervalIntegral.intervalIntegrable_rpow' hs

theorem T14.int_rpow' {s : ℝ} (hs : -1 < s) :
    IntervalIntegrable (fun x : ℝ => (1 - x) ^ s) volume 0 1 := by
  have := (T14.int_rpow hs).comp_sub_left 1
  simpa using this.symm

theorem T14.int_rpow_val {p : ℝ} (hp : 1 ≤ p) :
    ∫ t in (0:ℝ)..1, t ^ (1 / p - 1) = p := by
  have hr : 0 < 1 / p := by positivity
  rw [integral_rpow (Or.inl (by linarith))]
  rw [show 1 / p - 1 + 1 = 1 / p by ring, Real.zero_rpow hr.ne', Real.one_rpow, sub_zero]
  field_simp

theorem T14.int_rpow_half {p : ℝ} (hp : 1 ≤ p) :
    ∫ t in (0:ℝ)..(1 / 2), t ^ (1 / p - 1) = p * (2 : ℝ) ^ (-1 / p) := by
  have hr : 0 < 1 / p := by positivity
  have hc : (1 / 2 : ℝ) ^ (1 / p) = (2 : ℝ) ^ (-1 / p) := by
    rw [one_div, Real.inv_rpow (by norm_num), ← Real.rpow_neg (by norm_num)]
    congr 1; ring
  rw [integral_rpow (Or.inl (by linarith))]
  rw [show 1 / p - 1 + 1 = 1 / p by ring, Real.zero_rpow hr.ne', hc, sub_zero]
  field_simp

theorem T14.proof_4' : T14.stmt_4 := by
  intro p hp
  have hs : -1 < 1 / p - 1 := by have : 0 < 1 / p := by positivity
                                 linarith
  unfold T14.Lq
  simp only [div_one, Real.rpow_one]
  rw [intervalIntegral.integral_add (T14.int_rpow hs) (T14.int_rpow' hs)]
  have h2 : ∫ t in (0:ℝ)..1, (1 - t) ^ (1 / p - 1) = p := by
    rw [intervalIntegral.integral_comp_sub_left (fun x : ℝ => x ^ (1 / p - 1)) 1]
    simpa using T14.int_rpow_val hp
  rw [h2, T14.int_rpow_val hp]
  field_simp
  ring

theorem T14.proof_5' : T14.stmt_5 := by
  intro p hp
  have hr : 0 < 1 / p := by positivity
  have hr1 : 1 / p ≤ 1 := by rw [div_le_one (by linarith)]; exact hp
  set s : ℝ := 1 / p - 1 with hsdef
  have hs : -1 < s := by rw [hsdef]; linarith
  have hs0 : s ≤ 0 := by rw [hsdef]; linarith
  set F : ℝ → ℝ := fun t => max (t ^ s) ((1 - t) ^ s) with hF
  have hFi : IntervalIntegrable F volume 0 1 := by
    refine ((T14.int_rpow hs).add (T14.int_rpow' hs)).mono_fun' ?_ ?_
    · exact (by fun_prop : Measurable F).aestronglyMeasurable
    · refine ae_restrict_of_forall_mem measurableSet_uIoc fun x hx => ?_
      rw [Set.uIoc_of_le zero_le_one] at hx
      have a1 : 0 ≤ x ^ s := Real.rpow_nonneg hx.1.le _
      have a2 : 0 ≤ (1 - x) ^ s := Real.rpow_nonneg (by linarith [hx.2]) _
      simp only [hF, Real.norm_eq_abs, Pi.add_apply]
      rw [abs_of_nonneg (le_max_of_le_left a1)]
      rcases max_choice (x ^ s) ((1 - x) ^ s) with h | h <;> rw [h] <;> linarith
  have sub : ∀ a b : ℝ, 0 ≤ a → a ≤ b → b ≤ 1 → IntervalIntegrable F volume a b := by
    intro a b ha hab hb
    refine hFi.mono_set ?_
    rw [Set.uIcc_of_le zero_le_one, Set.uIcc_of_le hab]
    exact Icc_subset_Icc ha hb
  have hc : (1 / 2 : ℝ) ^ (1 / p) = (2 : ℝ) ^ (-1 / p) := by
    rw [one_div, Real.inv_rpow (by norm_num), ← Real.rpow_neg (by norm_num)]
    congr 1; ring
  show 4 / p * ∫ t in (0:ℝ)..1, F t = _
  rw [← intervalIntegral.integral_add_adjacent_intervals (b := 1 / 2)
    (sub 0 (1 / 2) le_rfl (by norm_num) (by norm_num)) (sub (1 / 2) 1 (by norm_num) (by norm_num) le_rfl)]
  have L : ∫ t in (0:ℝ)..(1 / 2), F t = p * (2 : ℝ) ^ (-1 / p) := by
    rw [← T14.int_rpow_half hp]
    refine intervalIntegral.integral_congr_Ioo_of_le (by norm_num) fun t ht => ?_
    simp only [hF]
    apply max_eq_left
    exact Real.rpow_le_rpow_of_nonpos ht.1 (by linarith [ht.2]) hs0
  have R : ∫ t in (1 / 2 : ℝ)..1, F t = p * (2 : ℝ) ^ (-1 / p) := by
    rw [← T14.int_rpow_half hp]
    have : ∫ t in (1 / 2 : ℝ)..1, F t = ∫ t in (1 / 2 : ℝ)..1, (1 - t) ^ s := by
      refine intervalIntegral.integral_congr_Ioo_of_le (by norm_num) fun t ht => ?_
      simp only [hF]
      apply max_eq_right
      exact Real.rpow_le_rpow_of_nonpos (by linarith [ht.2]) (by linarith [ht.1]) hs0
    rw [this, intervalIntegral.integral_comp_sub_left (fun x : ℝ => x ^ s) 1,
      show (1:ℝ) - 1 = 0 by norm_num, show (1:ℝ) - 1 / 2 = 1 / 2 by norm_num]
  rw [L, R]
  field_simp
  ring

theorem T14.pow_bounds {A B q : ℝ} (hA : 0 ≤ A) (hB : 0 ≤ B) (hq : 1 ≤ q) :
    max A B ≤ (A ^ q + B ^ q) ^ (1 / q) ∧ (A ^ q + B ^ q) ^ (1 / q) ≤ A + B := by
  have hq0 : 0 < q := by linarith
  have hM : 0 ≤ max A B := le_max_of_le_left hA
  have Aq : 0 ≤ A ^ q := Real.rpow_nonneg hA _
  have Bq : 0 ≤ B ^ q := Real.rpow_nonneg hB _
  constructor
  · have h1 : (max A B) ^ q ≤ A ^ q + B ^ q := by
      rcases max_choice A B with h | h <;> rw [h] <;> linarith
    calc max A B = ((max A B) ^ q) ^ (1 / q) := (T14.pw2 hq hM).symm
      _ ≤ _ := Real.rpow_le_rpow (Real.rpow_nonneg hM _) h1 (by positivity)
  · have h1 := Real.add_rpow_le_rpow_add hA hB hq
    calc (A ^ q + B ^ q) ^ (1 / q) ≤ ((A + B) ^ q) ^ (1 / q) :=
          Real.rpow_le_rpow (by positivity) h1 (by positivity)
      _ = A + B := T14.pw2 hq (by positivity)

theorem T14.lim_pow {A B : ℝ} (hA : 0 ≤ A) (hB : 0 ≤ B) :
    Tendsto (fun q : ℝ => (A ^ q + B ^ q) ^ (1 / q)) atTop (nhds (max A B)) := by
  have hM : 0 ≤ max A B := le_max_of_le_left hA
  have h2 : Tendsto (fun q : ℝ => (2 : ℝ) ^ (1 / q)) atTop (nhds 1) := by
    have h0 : Tendsto (fun q : ℝ => 1 / q) atTop (nhds 0) := by
      simpa [one_div] using tendsto_inv_atTop_zero
    simpa using (tendsto_const_nhds (x := (2 : ℝ))).rpow h0 (Or.inl two_ne_zero)
  have h3 : Tendsto (fun q : ℝ => (2 : ℝ) ^ (1 / q) * max A B) atTop (nhds (max A B)) := by
    simpa using h2.mul_const (max A B)
  refine tendsto_of_tendsto_of_tendsto_of_le_of_le' tendsto_const_nhds h3 ?_ ?_
  · filter_upwards [eventually_ge_atTop (1 : ℝ)] with q hq
    exact (T14.pow_bounds hA hB hq).1
  · filter_upwards [eventually_ge_atTop (1 : ℝ)] with q hq
    have hq0 : 0 < q := by linarith
    have Aq : 0 ≤ A ^ q := Real.rpow_nonneg hA _
    have Bq : 0 ≤ B ^ q := Real.rpow_nonneg hB _
    have h1 : A ^ q + B ^ q ≤ 2 * (max A B) ^ q := by
      have a := Real.rpow_le_rpow hA (le_max_left A B) hq0.le
      have b := Real.rpow_le_rpow hB (le_max_right A B) hq0.le
      linarith
    calc (A ^ q + B ^ q) ^ (1 / q) ≤ (2 * (max A B) ^ q) ^ (1 / q) :=
          Real.rpow_le_rpow (by positivity) h1 (by positivity)
      _ = 2 ^ (1 / q) * max A B := by
          rw [Real.mul_rpow (by norm_num) (Real.rpow_nonneg hM _), T14.pw2 hq hM]

theorem T14.rpow_q {p q t : ℝ} (ht : 0 ≤ t) :
    t ^ (q / p - q) = (t ^ (1 / p - 1)) ^ q := by
  rw [← Real.rpow_mul ht]; congr 1; ring

theorem T14.proof_6' : T14.stmt_6 := by
  intro p hp
  have hr : 0 < 1 / p := by positivity
  have hs : -1 < 1 / p - 1 := by linarith
  have key : Tendsto (fun q : ℝ => ∫ t in (0:ℝ)..1,
      (t ^ (q / p - q) + (1 - t) ^ (q / p - q)) ^ (1 / q)) atTop
      (nhds (∫ t in (0:ℝ)..1, max (t ^ (1 / p - 1)) ((1 - t) ^ (1 / p - 1)))) := by
    refine intervalIntegral.tendsto_integral_filter_of_dominated_convergence
      (fun t => t ^ (1 / p - 1) + (1 - t) ^ (1 / p - 1)) ?_ ?_
      ((T14.int_rpow hs).add (T14.int_rpow' hs)) ?_
    · exact Filter.Eventually.of_forall fun q =>
        (by fun_prop : Measurable fun t : ℝ =>
          (t ^ (q / p - q) + (1 - t) ^ (q / p - q)) ^ (1 / q)).aestronglyMeasurable
    · filter_upwards [eventually_ge_atTop (1 : ℝ)] with q hq
      refine ae_of_all _ fun x hx => ?_
      rw [Set.uIoc_of_le zero_le_one] at hx
      have a1 : 0 ≤ x := hx.1.le
      have a2 : 0 ≤ 1 - x := by linarith [hx.2]
      rw [T14.rpow_q a1, T14.rpow_q a2, Real.norm_eq_abs, abs_of_nonneg (by positivity)]
      exact (T14.pow_bounds (Real.rpow_nonneg a1 _) (Real.rpow_nonneg a2 _) hq).2
    · refine ae_of_all _ fun x hx => ?_
      rw [Set.uIoc_of_le zero_le_one] at hx
      have a1 : 0 ≤ x := hx.1.le
      have a2 : 0 ≤ 1 - x := by linarith [hx.2]
      have e : (fun q : ℝ => (x ^ (q / p - q) + (1 - x) ^ (q / p - q)) ^ (1 / q)) =
          fun q : ℝ => ((x ^ (1 / p - 1)) ^ q + ((1 - x) ^ (1 / p - 1)) ^ q) ^ (1 / q) := by
        funext q
        rw [T14.rpow_q (q := q) a1, T14.rpow_q (q := q) a2]
      rw [e]
      exact T14.lim_pow (Real.rpow_nonneg a1 _) (Real.rpow_nonneg a2 _)
  have := key.const_mul (4 / p)
  have h5 := T14.proof_5' p hp
  rw [h5] at this
  exact this
theorem T14.proof_1 : T14.stmt_1 := T14.proof_1'
theorem T14.proof_2 : T14.stmt_2 := T14.proof_2'
theorem T14.proof_3 : T14.stmt_3 := T14.proof_3'
theorem T14.proof_4 : T14.stmt_4 := T14.proof_4'
theorem T14.proof_5 : T14.stmt_5 := T14.proof_5'
theorem T14.proof_6 : T14.stmt_6 := T14.proof_6'
