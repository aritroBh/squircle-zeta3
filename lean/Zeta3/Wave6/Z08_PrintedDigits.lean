import Mathlib
import Zeta3.Wave3.W03_M1Full
import Zeta3.Wave5.Y01_M2
/-!
# Z08 -- the printed decimals of `M_0`, `κ = 8 M_1` and `M_2` are truncations

Target id: Z08.
Paper labels: the unlabelled display of `M_0` (main.tex 882), Theorem `thm:M1` (main.tex 889-895,
the decimal on line 894) and Theorem `thm:M2` (main.tex 901-907, the decimal on line 906), as
printed in draft 59 of the paper.

Quoted verbatim (main.tex 882):

  M_0=\lambda-\sqrt2\,L=-0.55330329972051571737\ldots,

Quoted verbatim (main.tex 894, Theorem thm:M1; draft 59 changed the last printed digit from
`5` to `4`):

  Numerically $\kappa=8M_1=2.02095276965635835694\ldots$

Quoted verbatim (main.tex 906, Theorem thm:M2):

  Numerically $M_2=-0.02216191108726939942\ldots$

## Reading of the dots

The decimals followed by `\ldots` are read as truncations to the 20 printed decimals: for a
positive value `x`, `printed ≤ x < printed + 10⁻²⁰`; for a negative value `x`,
`printed - 10⁻²⁰ < x ≤ printed` (the printed string is the true value with all later digits
dropped, sign kept).

## Objects used (imported, not redefined)

* `W02.Mk k = -F^{(k)}(0)/(2 k!)`, the coefficient `M_k` of eq:expansion (complex numbers).
* `W03.M0closed = D03.lam - √2 · D03.L = λ - √2 L` (main.tex 882); `W03.proof_2` gives
  `W02.Mk 0 = M0closed`.
* `D03.M1closed`, the closed form of thm:M1; `W03.proof_4` gives `W02.Mk 1 = M1closed`.
* `Y01.M2closed`, the closed form eq:M2 = (T1) of Table tab:closed, where `ζ(3)` enters as
  `Y01.z3 = (riemannZeta 3).re` (and `Y01.proof_4` shows `riemannZeta 3 = z3` with
  `HasSum (fun n => 1/(n+1)^3) z3`); `Y01.proof_12` gives `W02.Mk 2 = M2closed`.

## Statements

* stmt_1: the closed form `λ - √2 L` truncates to `-0.55330329972051571737`.
* stmt_2: `8 · M1closed` truncates to `2.02095276965635835694`.
* stmt_3: `M2closed` truncates to `-0.02216191108726939942`.
* stmt_4 (final): the coefficients themselves: `M_0 = W02.Mk 0`, `8 M_1 = 8 · W02.Mk 1` and
  `M_2 = W02.Mk 2` have imaginary part `0`, and their real parts satisfy the three truncation
  inequalities.

## Numerical check (Python `decimal`, 70 digits; mpmath was not available offline)

* `λ - √2 L   = -0.553303299720515717370808039042959437…`, minus printed `= -8.08·10⁻²²`.
* `8 M1closed =  2.020952769656358356946908167725238884…`, minus printed `= +6.91·10⁻²¹`.
* `M2closed   = -0.022161911087269399423754937739791514…`, minus printed `= -3.75·10⁻²¹`.
All three truncation claims hold. (The earlier draft's `…695` for `κ` was a rounding, not a
truncation; draft 59's `…694` is the truncation.)

## NOT COVERED

* Nothing of the three printed decimals is left out. The closed forms themselves
  (`M_0 = λ - √2 L`, thm:M1, thm:M2) are not restated here: they are W03 and Y01.
-/
-- BEGIN STATEMENTS
namespace Z08

/-- main.tex 882: "$M_0=\lambda-\sqrt2\,L=-0.55330329972051571737\ldots$", for the closed form
`λ - √2 L`: the printed decimal is its truncation (negative value). -/
def stmt_1 : Prop :=
  (-0.55330329972051571737 : ℝ) - (1 : ℝ) / 10 ^ 20 < W03.M0closed ∧
    W03.M0closed ≤ (-0.55330329972051571737 : ℝ)

/-- main.tex 894 (thm:M1): "Numerically $\kappa=8M_1=2.02095276965635835694\ldots$", for the
closed form of thm:M1: the printed decimal is the truncation of `8 · M1closed` (positive). -/
def stmt_2 : Prop :=
  (2.02095276965635835694 : ℝ) ≤ 8 * D03.M1closed ∧
    8 * D03.M1closed < (2.02095276965635835694 : ℝ) + (1 : ℝ) / 10 ^ 20

/-- main.tex 906 (thm:M2): "Numerically $M_2=-0.02216191108726939942\ldots$", for the closed
form eq:M2: the printed decimal is its truncation (negative value). -/
def stmt_3 : Prop :=
  (-0.02216191108726939942 : ℝ) - (1 : ℝ) / 10 ^ 20 < Y01.M2closed ∧
    Y01.M2closed ≤ (-0.02216191108726939942 : ℝ)

/-- main.tex 882, 894, 906: the three printed decimals are truncations of the true coefficients
`M_0 = W02.Mk 0`, `κ = 8 M_1 = 8 · W02.Mk 1`, `M_2 = W02.Mk 2` of eq:expansion, which are real
(imaginary part `0`). -/
def stmt_4 : Prop :=
  ((W02.Mk 0).im = 0 ∧
      (-0.55330329972051571737 : ℝ) - (1 : ℝ) / 10 ^ 20 < (W02.Mk 0).re ∧
      (W02.Mk 0).re ≤ (-0.55330329972051571737 : ℝ)) ∧
    ((8 * W02.Mk 1).im = 0 ∧
      (2.02095276965635835694 : ℝ) ≤ (8 * W02.Mk 1).re ∧
      (8 * W02.Mk 1).re < (2.02095276965635835694 : ℝ) + (1 : ℝ) / 10 ^ 20) ∧
    ((W02.Mk 2).im = 0 ∧
      (-0.02216191108726939942 : ℝ) - (1 : ℝ) / 10 ^ 20 < (W02.Mk 2).re ∧
      (W02.Mk 2).re ≤ (-0.02216191108726939942 : ℝ))

end Z08
-- END STATEMENTS

namespace Z08

open Real

lemma M0_trunc (s Lg : ℝ)
    (hs : (1.41421356237309504880168872 : ℝ) ≤ s ∧ s ≤ (1.41421356237309504880168873 : ℝ))
    (l : ℝ) (hl : (0.69314718055994530941723212 : ℝ) ≤ l ∧ l ≤ (0.69314718055994530941723213 : ℝ))
    (hL : (0.88137358701954302523260931 : ℝ) ≤ Lg ∧ Lg ≤ (0.88137358701954302523260934 : ℝ)) :
    (-0.55330329972051571737 : ℝ) - (1 : ℝ) / 10 ^ 20 < l - s * Lg ∧
      l - s * Lg ≤ (-0.55330329972051571737 : ℝ) := by
  obtain ⟨hs1, hs2⟩ := hs
  obtain ⟨hl1, hl2⟩ := hl
  obtain ⟨hL1, hL2⟩ := hL
  have p : s * Lg ≤ (1.41421356237309504880168873 : ℝ) * (0.88137358701954302523260934 : ℝ) :=
    mul_le_mul hs2 hL2 (by linarith) (by norm_num)
  have p' : (1.41421356237309504880168872 : ℝ) * (0.88137358701954302523260931 : ℝ) ≤ s * Lg :=
    mul_le_mul hs1 hL1 (by norm_num) (by linarith)
  constructor <;> norm_num at p p' ⊢ <;> linarith

lemma M1_trunc (s u w l Lg : ℝ) (hs : (1.41421356237309504880168872 : ℝ) ≤ s ∧ s ≤ (1.41421356237309504880168873 : ℝ))
    (hu : (0.46753399702300461762696479 : ℝ) ≤ u ∧ u ≤ (0.46753399702300461762696481 : ℝ)) (hw : (0.58224052646501250590265631 : ℝ) ≤ w ∧ w ≤ (0.58224052646501250590265633 : ℝ))
    (hl : (0.69314718055994530941723212 : ℝ) ≤ l ∧ l ≤ (0.69314718055994530941723213 : ℝ)) (hL : (0.88137358701954302523260931 : ℝ) ≤ Lg ∧ Lg ≤ (0.88137358701954302523260934 : ℝ)) :
    (2.02095276965635835694 : ℝ) ≤ 8 * ((19/4) * s * w - w / 2 + (19/8 * s + 1/4) * l ^ 2 + 2 * s * l * Lg - Lg ^ 2
      - (5/2) * s * Lg ^ 2 - 5 * s * u) ∧
    8 * ((19/4) * s * w - w / 2 + (19/8 * s + 1/4) * l ^ 2 + 2 * s * l * Lg - Lg ^ 2
      - (5/2) * s * Lg ^ 2 - 5 * s * u) < (2.02095276965635835694 : ℝ) + (1 : ℝ) / 10 ^ 20 := by
  obtain ⟨hs1, hs2⟩ := hs
  obtain ⟨hu1, hu2⟩ := hu
  obtain ⟨hw1, hw2⟩ := hw
  obtain ⟨hl1, hl2⟩ := hl
  obtain ⟨hL1, hL2⟩ := hL
  have hs0 : 0 ≤ s := by linarith
  have hu0 : 0 ≤ u := by linarith
  have hw0 : 0 ≤ w := by linarith
  have hl0 : 0 ≤ l := by linarith
  have hL0 : 0 ≤ Lg := by linarith
  have p1 : s * w ≤ (1.41421356237309504880168873 : ℝ) * (0.58224052646501250590265633 : ℝ) := mul_le_mul hs2 hw2 hw0 (by norm_num)
  have p1' : (1.41421356237309504880168872 : ℝ) * (0.58224052646501250590265631 : ℝ) ≤ s * w := mul_le_mul hs1 hw1 (by norm_num) hs0
  have q : l ^ 2 ≤ (0.69314718055994530941723213 : ℝ) ^ 2 := pow_le_pow_left₀ hl0 hl2 2
  have q' : (0.69314718055994530941723212 : ℝ) ^ 2 ≤ l ^ 2 := pow_le_pow_left₀ (by norm_num) hl1 2
  have p2 : s * l ^ 2 ≤ (1.41421356237309504880168873 : ℝ) * (0.69314718055994530941723213 : ℝ) ^ 2 := mul_le_mul hs2 q (by positivity) (by norm_num)
  have p2' : (1.41421356237309504880168872 : ℝ) * (0.69314718055994530941723212 : ℝ) ^ 2 ≤ s * l ^ 2 := mul_le_mul hs1 q' (by norm_num) hs0
  have r1 : s * l ≤ (1.41421356237309504880168873 : ℝ) * (0.69314718055994530941723213 : ℝ) := mul_le_mul hs2 hl2 hl0 (by norm_num)
  have r1' : (1.41421356237309504880168872 : ℝ) * (0.69314718055994530941723212 : ℝ) ≤ s * l := mul_le_mul hs1 hl1 (by norm_num) hs0
  have p3 : s * l * Lg ≤ (1.41421356237309504880168873 : ℝ) * (0.69314718055994530941723213 : ℝ) * (0.88137358701954302523260934 : ℝ) := mul_le_mul r1 hL2 hL0 (by norm_num)
  have p3' : (1.41421356237309504880168872 : ℝ) * (0.69314718055994530941723212 : ℝ) * (0.88137358701954302523260931 : ℝ) ≤ s * l * Lg := mul_le_mul r1' hL1 (by norm_num) (by positivity)
  have t : Lg ^ 2 ≤ (0.88137358701954302523260934 : ℝ) ^ 2 := pow_le_pow_left₀ hL0 hL2 2
  have t' : (0.88137358701954302523260931 : ℝ) ^ 2 ≤ Lg ^ 2 := pow_le_pow_left₀ (by norm_num) hL1 2
  have p4 : s * Lg ^ 2 ≤ (1.41421356237309504880168873 : ℝ) * (0.88137358701954302523260934 : ℝ) ^ 2 := mul_le_mul hs2 t (by positivity) (by norm_num)
  have p4' : (1.41421356237309504880168872 : ℝ) * (0.88137358701954302523260931 : ℝ) ^ 2 ≤ s * Lg ^ 2 := mul_le_mul hs1 t' (by norm_num) hs0
  have p5 : s * u ≤ (1.41421356237309504880168873 : ℝ) * (0.46753399702300461762696481 : ℝ) := mul_le_mul hs2 hu2 hu0 (by norm_num)
  have p5' : (1.41421356237309504880168872 : ℝ) * (0.46753399702300461762696479 : ℝ) ≤ s * u := mul_le_mul hs1 hu1 (by norm_num) hs0
  constructor
  · nlinarith
  · nlinarith

open Y01.PF in
set_option maxHeartbeats 4000000 in
lemma M2_trunc (s l L w u T P Q O : ℝ) (hs2 : s ^ 2 = 2)
    (hs : (1.41421356237309504880168872 : ℝ) ≤ s ∧ s ≤ (1.41421356237309504880168873 : ℝ))
    (hu : (0.46753399702300461762696479 : ℝ) ≤ u ∧ u ≤ (0.46753399702300461762696481 : ℝ))
    (hw : (0.58224052646501250590265631 : ℝ) ≤ w ∧ w ≤ (0.58224052646501250590265633 : ℝ))
    (hl : (0.69314718055994530941723212 : ℝ) ≤ l ∧ l ≤ (0.69314718055994530941723213 : ℝ))
    (hL : (0.88137358701954302523260931 : ℝ) ≤ L ∧ L ≤ (0.88137358701954302523260934 : ℝ))
    (hT : (0.537213193608040200940623182388 : ℝ) ≤ T ∧ T ≤ (0.537213193608040200940623227450 : ℝ))
    (hP : (0.438881700078672395657262704418 : ℝ) ≤ P ∧ P ≤ (0.438881700078672395657263102955 : ℝ))
    (hQ : (0.175454454019511266772526112476 : ℝ) ≤ Q ∧ Q ≤ (0.175454454019511266772526153583 : ℝ))
    (hO : (1.021002840767543829163857617042 : ℝ) ≤ O ∧ O ≤ (1.021002840767543829163857839087 : ℝ)) :
    let X := ((12 * w + 6 * l ^ 2) + 6 * l ^ 2) / 24
      + (189 * s / 64 - 19 / 16) * (8 / 7 * T + 8 / 7 * w * l + 8 / 21 * l ^ 3)
      + (11 + 15 * s) * P + (11 - 15 * s) * (Q / 4 - P)
      - 27 * s * (T / 8 + s * O / 2)
      + (12 * L - 7 * s / 2 * l) * u
      + (15 * s / 8 - 7 / 12) * (12 * w + 6 * l ^ 2) * L
      - (1 / 24 + 65 * s / 48) * (12 * w + 6 * l ^ 2) * l
      + (1 / 6 + 9 * s / 16) * l ^ 3
      + 11 * s / 8 * l ^ 2 * L
      - (1 + 7 * s / 4) * l * L ^ 2
      + (25 / 6 - 7 * s / 3) * L ^ 3
    (-0.02216191108726939942 : ℝ) - (1 : ℝ) / 10 ^ 20 < X ∧ X ≤ (-0.02216191108726939942 : ℝ) := by
  intro X
  simp only [X]
  have hsO : s ^ 2 * O = 2 * O := by rw [hs2]
  have m1 := (mb hL hu (by norm_num) (by norm_num))
  have m2 := (mb hL hw (by norm_num) (by norm_num))
  have m3 := (mb (mb hL hL (by norm_num) (by norm_num)) hL (by norm_num) (by norm_num))
  have m4 := (mb hl hw (by norm_num) (by norm_num))
  have m5 := (mb (mb hl hL (by norm_num) (by norm_num)) hL (by norm_num) (by norm_num))
  have m6 := (mb hl hl (by norm_num) (by norm_num))
  have m7 := (mb (mb hl hl (by norm_num) (by norm_num)) hL (by norm_num) (by norm_num))
  have m8 := (mb (mb hl hl (by norm_num) (by norm_num)) hl (by norm_num) (by norm_num))
  have m9 := (mb hs hQ (by norm_num) (by norm_num))
  have m10 := (mb hs hP (by norm_num) (by norm_num))
  have m11 := (mb (mb hs hL (by norm_num) (by norm_num)) hw (by norm_num) (by norm_num))
  have m12 := (mb (mb (mb hs hL (by norm_num) (by norm_num)) hL (by norm_num) (by norm_num)) hL (by norm_num) (by norm_num))
  have m13 := (mb (mb hs hl (by norm_num) (by norm_num)) hu (by norm_num) (by norm_num))
  have m14 := (mb (mb hs hl (by norm_num) (by norm_num)) hw (by norm_num) (by norm_num))
  have m15 := (mb (mb (mb hs hl (by norm_num) (by norm_num)) hL (by norm_num) (by norm_num)) hL (by norm_num) (by norm_num))
  have m16 := (mb (mb (mb hs hl (by norm_num) (by norm_num)) hl (by norm_num) (by norm_num)) hL (by norm_num) (by norm_num))
  have m17 := (mb (mb (mb hs hl (by norm_num) (by norm_num)) hl (by norm_num) (by norm_num)) hl (by norm_num) (by norm_num))
  constructor
  · linarith
  · linarith

end Z08

theorem Z08.proof_1 : Z08.stmt_1 := by
  unfold Z08.stmt_1 W03.M0closed
  exact Z08.M0_trunc _ _ D03.s_bounds _ D03.lam_bounds D03.L_bounds

theorem Z08.proof_2 : Z08.stmt_2 := by
  unfold Z08.stmt_2
  rw [D03.M1closed_expand]
  exact Z08.M1_trunc _ _ _ _ _ D03.s_bounds D03.u_bounds D03.w_bounds D03.lam_bounds D03.L_bounds

theorem Z08.proof_3 : Z08.stmt_3 := by
  unfold Z08.stmt_3
  rw [Y01.PF.M2closed_eq]
  exact Z08.M2_trunc _ _ _ _ _ _ _ _ _ Y01.PF.s_sq Y01.PF.s_bounds Y01.PF.u_bounds
    Y01.PF.w_bounds Y01.PF.lam_bounds Y01.PF.L_bounds Y01.PF.T_bounds Y01.PF.P_bounds
    Y01.PF.Q_bounds Y01.PF.O_bounds

theorem Z08.proof_4 : Z08.stmt_4 := by
  unfold Z08.stmt_4
  have h1 := Z08.proof_1
  have h2 := Z08.proof_2
  have h3 := Z08.proof_3
  unfold Z08.stmt_1 at h1
  unfold Z08.stmt_2 at h2
  unfold Z08.stmt_3 at h3
  rw [W03.proof_2, W03.proof_4, Y01.proof_12]
  have e : (8 : ℂ) * ((D03.M1closed : ℝ) : ℂ) = ((8 * D03.M1closed : ℝ) : ℂ) := by push_cast; ring
  rw [e]
  simp only [Complex.ofReal_im, Complex.ofReal_re, true_and]
  exact ⟨h1, h2, h3⟩
