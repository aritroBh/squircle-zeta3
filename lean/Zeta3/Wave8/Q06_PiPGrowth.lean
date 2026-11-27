import Mathlib
import Zeta3.Wave8.Q02_PiPAllOrders
/-!
# Q06 — continuation of Q02: `lem:pip:Tk`, `lem:pip:TkT1`, `lem:pip:twoterm`, the growth law and the divergence
(paper label `thm:pip`)

Target Q06 is the continuation of `Zeta3.Wave8.Q02_PiPAllOrders` (planner split, 2026-11-26).
The notation (`Q02.N`, `Q02.T`, `Q02.T1`, `Q02.psi`, …) is that of Q02; see its header for sources
and conventions. Statements, with the old Q02 number and the original docstring text:

* `stmt_1` (was `Q02.stmt_13`):
/-- Lemma `lem:pip:Tk`: `|N_k - (-log2)^{k+1}/(k+1)! - T_k| ≤ C (2^k + k (k-2)! 6^{-k})`. -/

* `stmt_2` (was `Q02.stmt_15`):
/-- Lemma `lem:pip:TkT1`: `T_k = T_k^{(1)} + O((k-2)! π^{-(k-1)} k^{-2} log k)`. -/

* `stmt_3` (was `Q02.stmt_16`):
/-- Lemma `lem:pip:twoterm`: with `y_0 = iπ`, `X = ψ(k-2) + γ - log π`,
`N_k = 2 Re[(k-2)! y_0^{1-k} + (k-3)! y_0^{2-k}(X - iπ/2)] + O((k-2)! π^{-(k-1)} log²k/k²)`. -/

* `stmt_4` (was `Q02.stmt_17`):
/-- thm:pip / Lemma `prop:pip:growth`, odd `k`:
`N_k = (-1)^{(k-1)/2} (2(k-2)!/π^{k-1}) (1 + π²/(2(k-2)) + O(log²k/k²))`. -/

* `stmt_5` (was `Q02.stmt_18`):
/-- thm:pip / Lemma `prop:pip:growth`, even `k`:
`N_k = (-1)^{(k-2)/2} (2(k-3)!/π^{k-2}) (ψ(k-2) + γ - log π + O(log²k/k))`. -/

* `stmt_6` (was `Q02.stmt_19`):
/-- thm:pip / Lemma `prop:pip:div`: for every real `p ≠ 0` the terms `N_k p^{-k-1}` do not tend
to `0` (so `Σ N_k p^{-k-1}` diverges), and for every `x ≠ 0` the terms `N_k x^{k+1}` are unbounded
(radius of convergence `0`). -/

* `stmt_7` (was `Q02.stmt_20`):
/-- thm:pip, the part not proved in `X02` (main.tex lines 1344–1349), in full: the expansion to
every order, the divergence of the series for every `p`, and the growth law for odd and even `k`. -/
-/
-- BEGIN STATEMENTS
namespace Q06

/-- Lemma `lem:pip:Tk`: `|N_k - (-log2)^{k+1}/(k+1)! - T_k| ≤ C (2^k + k (k-2)! 6^{-k})`. (was Q02.stmt_13) -/
def stmt_1 : Prop :=
  ∃ C : ℝ, ∀ k : ℕ, 2 ≤ k →
    |Q02.N k - (-Real.log 2) ^ (k + 1) / ((k + 1).factorial : ℝ) - Q02.T k|
      ≤ C * (2 ^ k + (k : ℝ) * ((k - 2).factorial : ℝ) / 6 ^ k)

/-- Lemma `lem:pip:TkT1`: `T_k = T_k^{(1)} + O((k-2)! π^{-(k-1)} k^{-2} log k)`. (was Q02.stmt_15) -/
def stmt_2 : Prop :=
  ∃ C : ℝ, ∃ k₀ : ℕ, ∀ k : ℕ, k₀ ≤ k →
    |Q02.T k - Q02.T1 k| ≤ C * (((k - 2).factorial : ℝ) / Real.pi ^ (k - 1) * (Real.log k / (k : ℝ) ^ 2))

/-- Lemma `lem:pip:twoterm`: with `y_0 = iπ`, `X = ψ(k-2) + γ - log π`,
`N_k = 2 Re[(k-2)! y_0^{1-k} + (k-3)! y_0^{2-k}(X - iπ/2)] + O((k-2)! π^{-(k-1)} log²k/k²)`.
(was Q02.stmt_16) -/
def stmt_3 : Prop :=
  ∃ C : ℝ, ∃ k₀ : ℕ, ∀ k : ℕ, k₀ ≤ k →
    |Q02.N k - 2 * ((((k - 2).factorial : ℂ) * ((Complex.I * Real.pi) ^ (k - 1))⁻¹ +
        ((k - 3).factorial : ℂ) * ((Complex.I * Real.pi) ^ (k - 2))⁻¹ *
          (((Q02.psi ((k : ℝ) - 2) + Real.eulerMascheroniConstant - Real.log Real.pi : ℝ) : ℂ) -
            Complex.I * Real.pi / 2)).re)|
      ≤ C * (((k - 2).factorial : ℝ) / Real.pi ^ (k - 1) * (Real.log k ^ 2 / (k : ℝ) ^ 2))

/-- thm:pip / Lemma `prop:pip:growth`, odd `k`:
`N_k = (-1)^{(k-1)/2} (2(k-2)!/π^{k-1}) (1 + π²/(2(k-2)) + O(log²k/k²))`.
(was Q02.stmt_17) -/
def stmt_4 : Prop :=
  ∃ C : ℝ, ∃ k₀ : ℕ, ∀ k : ℕ, k₀ ≤ k → Odd k →
    |Q02.N k - (-1 : ℝ) ^ ((k - 1) / 2) * (2 * ((k - 2).factorial : ℝ) / Real.pi ^ (k - 1)) *
        (1 + Real.pi ^ 2 / (2 * ((k : ℝ) - 2)))|
      ≤ C * (2 * ((k - 2).factorial : ℝ) / Real.pi ^ (k - 1)) * (Real.log k ^ 2 / (k : ℝ) ^ 2)

/-- thm:pip / Lemma `prop:pip:growth`, even `k`:
`N_k = (-1)^{(k-2)/2} (2(k-3)!/π^{k-2}) (ψ(k-2) + γ - log π + O(log²k/k))`.
(was Q02.stmt_18) -/
def stmt_5 : Prop :=
  ∃ C : ℝ, ∃ k₀ : ℕ, ∀ k : ℕ, k₀ ≤ k → Even k →
    |Q02.N k - (-1 : ℝ) ^ ((k - 2) / 2) * (2 * ((k - 3).factorial : ℝ) / Real.pi ^ (k - 2)) *
        (Q02.psi ((k : ℝ) - 2) + Real.eulerMascheroniConstant - Real.log Real.pi)|
      ≤ C * (2 * ((k - 3).factorial : ℝ) / Real.pi ^ (k - 2)) * (Real.log k ^ 2 / (k : ℝ))

/-- thm:pip / Lemma `prop:pip:div`: for every real `p ≠ 0` the terms `N_k p^{-k-1}` do not tend
to `0` (so `Σ N_k p^{-k-1}` diverges), and for every `x ≠ 0` the terms `N_k x^{k+1}` are unbounded
(radius of convergence `0`).
(was Q02.stmt_19) -/
def stmt_6 : Prop :=
  (∀ p : ℝ, p ≠ 0 →
    ¬ Filter.Tendsto (fun k : ℕ => Q02.N k * (1 / p) ^ (k + 1)) Filter.atTop (nhds 0)) ∧
  (∀ x : ℝ, x ≠ 0 → ¬ BddAbove (Set.range fun k : ℕ => |Q02.N k * x ^ (k + 1)|))

/-- thm:pip, the part not proved in `X02` (main.tex lines 1344–1349), in full: the expansion to
every order, the divergence of the series for every `p`, and the growth law for odd and even `k`.
(was Q02.stmt_20) -/
def stmt_7 : Prop :=
  Q02.stmt_8 ∧ stmt_6 ∧ stmt_4 ∧ stmt_5

end Q06
-- END STATEMENTS

section

open MeasureTheory Set Real
open Filter Topology Finset

namespace Q06
open Q02

/-! ### `prop:pip:div` from the growth law for odd `k` -/

lemma log_sq_div_le {k : ℝ} (hk : 1 ≤ k) : Real.log k ^ 2 / k ^ 2 ≤ 4 / k := by
  have hk0 : 0 < k := by linarith
  have hl0 : 0 ≤ Real.log k := Real.log_nonneg hk
  have h1 : Real.log k ≤ k ^ (1 / 2 : ℝ) / (1 / 2) :=
    Real.log_le_rpow_div hk0.le (by norm_num)
  have h2 : (k ^ (1 / 2 : ℝ)) ^ 2 = k := by
    rw [← Real.rpow_natCast, ← Real.rpow_mul hk0.le]; norm_num
  have h3 : Real.log k ^ 2 ≤ 4 * k := by
    have := pow_le_pow_left₀ hl0 h1 2
    calc Real.log k ^ 2 ≤ (k ^ (1 / 2 : ℝ) / (1 / 2)) ^ 2 := this
      _ = 4 * (k ^ (1 / 2 : ℝ)) ^ 2 := by ring
      _ = 4 * k := by rw [h2]
  rw [div_le_div_iff₀ (by positivity) hk0]
  nlinarith

/-- a sequence with the growth law of `stmt_4` has `|a_k x^{k+1}|` unbounded for `x ≠ 0` -/
lemma unbounded_of_growth (a : ℕ → ℝ)
    (h : ∃ C : ℝ, ∃ k₀ : ℕ, ∀ k : ℕ, k₀ ≤ k → Odd k →
      |a k - (-1 : ℝ) ^ ((k - 1) / 2) * (2 * ((k - 2).factorial : ℝ) / Real.pi ^ (k - 1)) *
          (1 + Real.pi ^ 2 / (2 * ((k : ℝ) - 2)))|
        ≤ C * (2 * ((k - 2).factorial : ℝ) / Real.pi ^ (k - 1)) * (Real.log k ^ 2 / (k : ℝ) ^ 2))
    {x : ℝ} (hx : x ≠ 0) : ¬ BddAbove (Set.range fun k : ℕ => |a k * x ^ (k + 1)|) := by
  obtain ⟨C, k₀, hC⟩ := h
  rintro ⟨B, hB⟩
  have hB' : ∀ k, |a k * x ^ (k + 1)| ≤ B := fun k => hB ⟨k, rfl⟩
  have hpi := Real.pi_pos
  have hxa : 0 < |x| := abs_pos.mpr hx
  -- `r^n / n! → 0` for `r = π / |x|`
  set r := Real.pi / |x| with hr
  have ht := FloorSemiring.tendsto_pow_div_factorial_atTop r
  set B' := |B| * Real.pi / |x| ^ 3 + 1 with hB'def
  have hB'0 : 0 < B' := by positivity
  obtain ⟨n₀, hn₀⟩ := (ht.eventually (gt_mem_nhds (show (0 : ℝ) < 1 / B' by positivity))).exists_forall_of_atTop
  -- choose an odd `k` large enough
  set k := 2 * (n₀ + k₀ + 8 * Nat.ceil |C| + 10) + 1 with hk
  have hkodd : Odd k := ⟨_, rfl⟩
  have hk0 : k₀ ≤ k := by omega
  have hn : n₀ ≤ k - 2 := by omega
  have hk3 : 3 ≤ k := by omega
  have hkC : 8 * |C| ≤ (k : ℝ) := by
    have h1 : |C| ≤ (Nat.ceil |C| : ℝ) := Nat.le_ceil _
    have h2 : (8 * Nat.ceil |C| : ℕ) ≤ k := by omega
    have h3 : ((8 * Nat.ceil |C| : ℕ) : ℝ) ≤ (k : ℝ) := by exact_mod_cast h2
    push_cast at h3; linarith
  have hkR : (3 : ℝ) ≤ (k : ℝ) := by exact_mod_cast hk3
  set A := 2 * ((k - 2).factorial : ℝ) / Real.pi ^ (k - 1) with hA
  have hA0 : 0 < A := by positivity
  have hgr := hC k hk0 hkodd
  -- lower bound `|a k| ≥ A / 2`
  have hlog := log_sq_div_le (k := (k : ℝ)) (by linarith)
  have herr : C * A * (Real.log k ^ 2 / (k : ℝ) ^ 2) ≤ A / 2 := by
    have h1 : C * A * (Real.log k ^ 2 / (k : ℝ) ^ 2) ≤ |C| * A * (4 / k) := by
      have hl : 0 ≤ Real.log k ^ 2 / (k : ℝ) ^ 2 := by positivity
      calc C * A * (Real.log k ^ 2 / (k : ℝ) ^ 2) ≤ |C| * A * (Real.log k ^ 2 / (k : ℝ) ^ 2) :=
            mul_le_mul_of_nonneg_right (mul_le_mul_of_nonneg_right (le_abs_self C) hA0.le) hl
        _ ≤ |C| * A * (4 / k) := mul_le_mul_of_nonneg_left hlog (by positivity)
    have h2 : |C| * A * (4 / k) ≤ A / 2 := by
      rw [show |C| * A * (4 / (k : ℝ)) = A * (4 * |C|) / k by ring, div_le_iff₀ (by linarith)]
      nlinarith
    linarith
  have hmain : A ≤ |(-1 : ℝ) ^ ((k - 1) / 2) * A * (1 + Real.pi ^ 2 / (2 * ((k : ℝ) - 2)))| := by
    rw [abs_mul, abs_mul, abs_pow, abs_neg, abs_one, one_pow, one_mul, abs_of_pos hA0,
      abs_of_pos (by have : (0 : ℝ) < 2 * ((k : ℝ) - 2) := by linarith
                     positivity)]
    have : 0 ≤ Real.pi ^ 2 / (2 * ((k : ℝ) - 2)) := by
      have : (0 : ℝ) < 2 * ((k : ℝ) - 2) := by linarith
      positivity
    nlinarith
  have hak : A / 2 ≤ |a k| := by
    have := abs_sub_abs_le_abs_sub ((-1 : ℝ) ^ ((k - 1) / 2) * A *
      (1 + Real.pi ^ 2 / (2 * ((k : ℝ) - 2)))) (a k)
    rw [abs_sub_comm] at this
    linarith
  -- upper bound from boundedness
  have hbd := hB' k
  rw [abs_mul, abs_pow] at hbd
  have hfac : ((k - 2).factorial : ℝ) * |x| ^ (k + 1) / Real.pi ^ (k - 1) ≤ |B| := by
    have : A / 2 * |x| ^ (k + 1) ≤ |a k| * |x| ^ (k + 1) :=
      mul_le_mul_of_nonneg_right hak (by positivity)
    have e : A / 2 * |x| ^ (k + 1) = ((k - 2).factorial : ℝ) * |x| ^ (k + 1) / Real.pi ^ (k - 1) := by
      rw [hA]; ring
    linarith [le_abs_self B]
  -- contradiction with `r^{k-2}/(k-2)! < 1/B'`
  have hsmall := hn₀ (k - 2) hn
  have hpos : (0 : ℝ) < ((k - 2).factorial : ℝ) := by exact_mod_cast Nat.factorial_pos _
  have e1 : k - 1 = (k - 2) + 1 := by omega
  have e2 : k + 1 = (k - 2) + 3 := by omega
  rw [e1, e2] at hfac
  have hrk : r ^ (k - 2) = Real.pi ^ (k - 2) / |x| ^ (k - 2) := by rw [hr, div_pow]
  rw [div_lt_div_iff₀ hpos hB'0, one_mul, hrk] at hsmall
  have hxn : 0 < |x| ^ (k - 2) := by positivity
  have hs2 : Real.pi ^ (k - 2) * B' < ((k - 2).factorial : ℝ) * |x| ^ (k - 2) := by
    have := mul_lt_mul_of_pos_right hsmall hxn
    rw [div_mul_eq_mul_div, div_mul_cancel₀ _ hxn.ne'] at this
    linarith
  have hfac' : ((k - 2).factorial : ℝ) * |x| ^ (k - 2) * |x| ^ 3 ≤
      |B| * Real.pi * Real.pi ^ (k - 2) := by
    rw [pow_add, pow_add, pow_one] at hfac
    rw [div_le_iff₀ (by positivity)] at hfac
    nlinarith
  have hB'e : B' * (|x| ^ 3) = |B| * Real.pi + |x| ^ 3 := by
    rw [hB'def]; field_simp
  have hx3 : 0 < |x| ^ 3 := by positivity
  have h5 : Real.pi ^ (k - 2) * B' * |x| ^ 3 < ((k - 2).factorial : ℝ) * |x| ^ (k - 2) * |x| ^ 3 :=
    mul_lt_mul_of_pos_right hs2 hx3
  have h6 : Real.pi ^ (k - 2) * B' * |x| ^ 3 = Real.pi ^ (k - 2) * (|B| * Real.pi + |x| ^ 3) := by
    rw [mul_assoc, hB'e]
  have h7 : 0 < Real.pi ^ (k - 2) * |x| ^ 3 := by positivity
  nlinarith

lemma div_of_growth (a : ℕ → ℝ)
    (h : ∃ C : ℝ, ∃ k₀ : ℕ, ∀ k : ℕ, k₀ ≤ k → Odd k →
      |a k - (-1 : ℝ) ^ ((k - 1) / 2) * (2 * ((k - 2).factorial : ℝ) / Real.pi ^ (k - 1)) *
          (1 + Real.pi ^ 2 / (2 * ((k : ℝ) - 2)))|
        ≤ C * (2 * ((k - 2).factorial : ℝ) / Real.pi ^ (k - 1)) * (Real.log k ^ 2 / (k : ℝ) ^ 2)) :
    (∀ p : ℝ, p ≠ 0 →
      ¬ Filter.Tendsto (fun k : ℕ => a k * (1 / p) ^ (k + 1)) Filter.atTop (nhds 0)) ∧
    (∀ x : ℝ, x ≠ 0 → ¬ BddAbove (Set.range fun k : ℕ => |a k * x ^ (k + 1)|)) := by
  refine ⟨fun p hp ht => ?_, fun x hx => unbounded_of_growth a h hx⟩
  have hx : (1 / p) ≠ 0 := one_div_ne_zero hp
  apply unbounded_of_growth a h hx
  exact ht.abs.bddAbove_range

end Q06

namespace Q06
open Q02

/-! ### `prop:pip:growth` from `lem:pip:twoterm`: the real parts of the two terms -/

lemma ipow_even (n : ℕ) :
    (Complex.I * (Real.pi : ℂ)) ^ (2 * n) = (((-1 : ℝ) ^ n * Real.pi ^ (2 * n) : ℝ) : ℂ) := by
  rw [mul_pow, pow_mul, Complex.I_sq]; push_cast; ring

lemma ipow_odd (n : ℕ) :
    (Complex.I * (Real.pi : ℂ)) ^ (2 * n + 1) =
      (((-1 : ℝ) ^ n * Real.pi ^ (2 * n + 1) : ℝ) : ℂ) * Complex.I := by
  rw [pow_succ, ipow_even]; push_cast; ring

lemma inv_ipow_even (n : ℕ) :
    ((Complex.I * (Real.pi : ℂ)) ^ (2 * n))⁻¹ = (((-1 : ℝ) ^ n * Real.pi ^ (2 * n))⁻¹ : ℝ) := by
  rw [ipow_even, Complex.ofReal_inv]

lemma inv_ipow_odd (n : ℕ) :
    ((Complex.I * (Real.pi : ℂ)) ^ (2 * n + 1))⁻¹ =
      (((-1 : ℝ) ^ n * Real.pi ^ (2 * n + 1))⁻¹ : ℝ) * (-Complex.I) := by
  rw [ipow_odd, mul_inv, Complex.ofReal_inv, Complex.inv_I]

lemma neg_one_pow_inv (n : ℕ) : ((-1 : ℝ) ^ n)⁻¹ = (-1) ^ n := by
  rcases neg_one_pow_eq_or ℝ n with h | h <;> rw [h] <;> norm_num

/-- odd `k = 2i+3`: the two-term main part equals the main part of `stmt_4` -/
lemma main_odd (i : ℕ) (X : ℝ) :
    2 * (((2 * i + 1).factorial : ℂ) * ((Complex.I * Real.pi) ^ (2 * i + 2))⁻¹ +
        ((2 * i).factorial : ℂ) * ((Complex.I * Real.pi) ^ (2 * i + 1))⁻¹ *
          ((X : ℂ) - Complex.I * Real.pi / 2)).re =
      (-1 : ℝ) ^ (i + 1) * (2 * ((2 * i + 1).factorial : ℝ) / Real.pi ^ (2 * i + 2)) *
        (1 + Real.pi ^ 2 / (2 * (2 * (i : ℝ) + 1))) := by
  rw [show 2 * i + 2 = 2 * (i + 1) by ring, inv_ipow_even, inv_ipow_odd]
  have hre : (((2 * i + 1).factorial : ℂ) *
        ((((-1 : ℝ) ^ (i + 1) * Real.pi ^ (2 * (i + 1)))⁻¹ : ℝ) : ℂ) +
      ((2 * i).factorial : ℂ) * (((((-1 : ℝ) ^ i * Real.pi ^ (2 * i + 1))⁻¹ : ℝ) : ℂ) *
        (-Complex.I)) * ((X : ℂ) - Complex.I * Real.pi / 2)).re =
      ((2 * i + 1).factorial : ℝ) * ((-1 : ℝ) ^ (i + 1) * Real.pi ^ (2 * (i + 1)))⁻¹ -
        ((2 * i).factorial : ℝ) * ((-1 : ℝ) ^ i * Real.pi ^ (2 * i + 1))⁻¹ * (Real.pi / 2) := by
    generalize ((-1 : ℝ) ^ (i + 1) * Real.pi ^ (2 * (i + 1)))⁻¹ = a
    generalize ((-1 : ℝ) ^ i * Real.pi ^ (2 * i + 1))⁻¹ = b
    generalize ((2 * i + 1).factorial : ℕ) = m
    generalize ((2 * i).factorial : ℕ) = n
    simp [Complex.mul_re, Complex.mul_im, Complex.sub_re, Complex.sub_im]
    ring
  rw [hre]
  have hpi := Real.pi_pos
  have hf : ((2 * i + 1).factorial : ℝ) = (2 * (i : ℝ) + 1) * ((2 * i).factorial : ℝ) := by
    rw [Nat.factorial_succ]; push_cast; ring
  rw [hf, mul_inv, mul_inv, neg_one_pow_inv, neg_one_pow_inv]
  have h2i : (0 : ℝ) < 2 * (i : ℝ) + 1 := by positivity
  field_simp
  ring

/-- even `k = 2i+4`: the two-term main part equals the main part of `stmt_5` -/
lemma main_even (i : ℕ) (X : ℝ) :
    2 * (((2 * i + 2).factorial : ℂ) * ((Complex.I * Real.pi) ^ (2 * (i + 1) + 1))⁻¹ +
        ((2 * i + 1).factorial : ℂ) * ((Complex.I * Real.pi) ^ (2 * (i + 1)))⁻¹ *
          ((X : ℂ) - Complex.I * Real.pi / 2)).re =
      (-1 : ℝ) ^ (i + 1) * (2 * ((2 * i + 1).factorial : ℝ) / Real.pi ^ (2 * (i + 1))) * X := by
  rw [inv_ipow_even, inv_ipow_odd]
  have hre : (((2 * i + 2).factorial : ℂ) *
        (((((-1 : ℝ) ^ (i + 1) * Real.pi ^ (2 * (i + 1) + 1))⁻¹ : ℝ) : ℂ) * (-Complex.I)) +
      ((2 * i + 1).factorial : ℂ) * ((((-1 : ℝ) ^ (i + 1) * Real.pi ^ (2 * (i + 1)))⁻¹ : ℝ) : ℂ) *
        ((X : ℂ) - Complex.I * Real.pi / 2)).re =
      ((2 * i + 1).factorial : ℝ) * ((-1 : ℝ) ^ (i + 1) * Real.pi ^ (2 * (i + 1)))⁻¹ * X := by
    generalize ((-1 : ℝ) ^ (i + 1) * Real.pi ^ (2 * (i + 1) + 1))⁻¹ = a
    generalize ((-1 : ℝ) ^ (i + 1) * Real.pi ^ (2 * (i + 1)))⁻¹ = b
    generalize ((2 * i + 2).factorial : ℕ) = m
    generalize ((2 * i + 1).factorial : ℕ) = n
    simp [Complex.mul_re, Complex.mul_im, Complex.sub_re, Complex.sub_im]
  rw [hre, mul_inv, neg_one_pow_inv]
  have hpi := Real.pi_pos
  field_simp


/-- the hypothesis shape of `stmt_3` for a sequence `a` and a function `ps` (`= ψ`) -/
def TwoTerm (a : ℕ → ℝ) (ps : ℝ → ℝ) : Prop :=
  ∃ C : ℝ, ∃ k₀ : ℕ, ∀ k : ℕ, k₀ ≤ k →
    |a k - 2 * ((((k - 2).factorial : ℂ) * ((Complex.I * Real.pi) ^ (k - 1))⁻¹ +
        ((k - 3).factorial : ℂ) * ((Complex.I * Real.pi) ^ (k - 2))⁻¹ *
          (((ps ((k : ℝ) - 2) + Real.eulerMascheroniConstant - Real.log Real.pi : ℝ) : ℂ) -
            Complex.I * Real.pi / 2)).re)|
      ≤ C * (((k - 2).factorial : ℝ) / Real.pi ^ (k - 1) * (Real.log k ^ 2 / (k : ℝ) ^ 2))

lemma growth_odd (a : ℕ → ℝ) (ps : ℝ → ℝ) (h : TwoTerm a ps) :
    ∃ C : ℝ, ∃ k₀ : ℕ, ∀ k : ℕ, k₀ ≤ k → Odd k →
      |a k - (-1 : ℝ) ^ ((k - 1) / 2) * (2 * ((k - 2).factorial : ℝ) / Real.pi ^ (k - 1)) *
          (1 + Real.pi ^ 2 / (2 * ((k : ℝ) - 2)))|
        ≤ C * (2 * ((k - 2).factorial : ℝ) / Real.pi ^ (k - 1)) * (Real.log k ^ 2 / (k : ℝ) ^ 2) := by
  obtain ⟨C, k₀, hC⟩ := h
  refine ⟨C / 2, max k₀ 3, fun k hk hodd => ?_⟩
  have hk0 : k₀ ≤ k := le_of_max_le_left hk
  have hk3 : 3 ≤ k := le_of_max_le_right hk
  have H := hC k hk0
  obtain ⟨m, rfl⟩ := hodd
  obtain ⟨i, rfl⟩ : ∃ i, m = i + 1 := ⟨m - 1, by omega⟩
  have e1 : 2 * (i + 1) + 1 - 2 = 2 * i + 1 := by omega
  have e2 : 2 * (i + 1) + 1 - 1 = 2 * i + 2 := by omega
  have e3 : 2 * (i + 1) + 1 - 3 = 2 * i := by omega
  have e4 : (2 * (i + 1) + 1 - 1) / 2 = i + 1 := by omega
  have e5 : (((2 * (i + 1) + 1 : ℕ) : ℝ) - 2) = 2 * (i : ℝ) + 1 := by push_cast; ring
  rw [e1, e2, e3, e5] at H
  rw [e4, e1, e2, e5]
  rw [main_odd] at H
  calc _ ≤ _ := H
    _ = _ := by ring

lemma growth_even (a : ℕ → ℝ) (ps : ℝ → ℝ) (h : TwoTerm a ps) :
    ∃ C : ℝ, ∃ k₀ : ℕ, ∀ k : ℕ, k₀ ≤ k → Even k →
      |a k - (-1 : ℝ) ^ ((k - 2) / 2) * (2 * ((k - 3).factorial : ℝ) / Real.pi ^ (k - 2)) *
          (ps ((k : ℝ) - 2) + Real.eulerMascheroniConstant - Real.log Real.pi)|
        ≤ C * (2 * ((k - 3).factorial : ℝ) / Real.pi ^ (k - 2)) * (Real.log k ^ 2 / (k : ℝ)) := by
  obtain ⟨C, k₀, hC⟩ := h
  refine ⟨|C| / (2 * Real.pi), max k₀ 4, fun k hk hev => ?_⟩
  have hk0 : k₀ ≤ k := le_of_max_le_left hk
  have hk4 : 4 ≤ k := le_of_max_le_right hk
  have H := hC k hk0
  obtain ⟨m, rfl⟩ := hev
  obtain ⟨i, rfl⟩ : ∃ i, m = i + 2 := ⟨m - 2, by omega⟩
  have e1 : i + 2 + (i + 2) - 2 = 2 * (i + 1) := by omega
  have e2 : i + 2 + (i + 2) - 1 = 2 * (i + 1) + 1 := by omega
  have e3 : i + 2 + (i + 2) - 3 = 2 * i + 1 := by omega
  have e4 : (i + 2 + (i + 2) - 2) / 2 = i + 1 := by omega
  have e6 : (2 * (i + 1)).factorial = (2 * i + 2).factorial := by ring_nf
  rw [e1, e2, e3] at H
  rw [e4, e1, e3]
  rw [e6, main_even] at H
  refine H.trans ?_
  -- compare the error terms
  have hpi := Real.pi_pos
  set k := i + 2 + (i + 2) with hk
  have hkR : (k : ℝ) = 2 * (i : ℝ) + 4 := by rw [hk]; push_cast; ring
  set L := Real.log (k : ℝ) ^ 2 with hLdef
  have hL : 0 ≤ L := sq_nonneg _
  have hf : ((2 * i + 2).factorial : ℝ) = (2 * (i : ℝ) + 2) * ((2 * i + 1).factorial : ℝ) := by
    rw [show 2 * i + 2 = (2 * i + 1) + 1 by ring, Nat.factorial_succ]; push_cast; ring
  have hp : Real.pi ^ (2 * (i + 1) + 1) = Real.pi * Real.pi ^ (2 * (i + 1)) := by ring
  rw [hf, hp, hkR]
  have hF : (0 : ℝ) < ((2 * i + 1).factorial : ℝ) := by exact_mod_cast Nat.factorial_pos _
  have hP : (0 : ℝ) < Real.pi ^ (2 * (i + 1)) := by positivity
  have hi : (0 : ℝ) ≤ (i : ℝ) := Nat.cast_nonneg i
  calc C * ((2 * (i : ℝ) + 2) * ((2 * i + 1).factorial : ℝ) / (Real.pi * Real.pi ^ (2 * (i + 1))) *
        (L / (2 * (i : ℝ) + 4) ^ 2))
      ≤ |C| * ((2 * (i : ℝ) + 2) * ((2 * i + 1).factorial : ℝ) / (Real.pi * Real.pi ^ (2 * (i + 1))) *
        (L / (2 * (i : ℝ) + 4) ^ 2)) :=
        mul_le_mul_of_nonneg_right (le_abs_self C) (by positivity)
    _ ≤ |C| / (2 * Real.pi) * (2 * ((2 * i + 1).factorial : ℝ) / Real.pi ^ (2 * (i + 1))) *
        (L / (2 * (i : ℝ) + 4)) := by
        have key : (2 * (i : ℝ) + 2) / (2 * (i : ℝ) + 4) ^ 2 ≤ 1 / (2 * (i : ℝ) + 4) := by
          rw [div_le_div_iff₀ (by positivity) (by positivity)]; nlinarith
        have e : |C| * ((2 * (i : ℝ) + 2) * ((2 * i + 1).factorial : ℝ) /
            (Real.pi * Real.pi ^ (2 * (i + 1))) * (L / (2 * (i : ℝ) + 4) ^ 2)) =
            (|C| * ((2 * i + 1).factorial : ℝ) / (Real.pi * Real.pi ^ (2 * (i + 1))) *
              L) * ((2 * (i : ℝ) + 2) / (2 * (i : ℝ) + 4) ^ 2) := by ring
        have e' : |C| / (2 * Real.pi) * (2 * ((2 * i + 1).factorial : ℝ) / Real.pi ^ (2 * (i + 1))) *
            (L / (2 * (i : ℝ) + 4)) =
            (|C| * ((2 * i + 1).factorial : ℝ) / (Real.pi * Real.pi ^ (2 * (i + 1))) *
              L) * (1 / (2 * (i : ℝ) + 4)) := by field_simp
        rw [e, e']
        exact mul_le_mul_of_nonneg_left key (by positivity)

end Q06

end

section
-- merged from work/Q06/PartA.lean
/-!
# Q06, part A: Lemma `lem:pip:Tk`

`|N_k - (-log 2)^(k+1)/(k+1)! - T_k| ≤ C (2^k + k (k-2)! 6^{-k})` for `k ≥ 2`.

Plan (sec_pip.tex, `lem:pip:laplace`, `lem:pip:watson`, `lem:pip:Tk`):
* `Ehat h = h² I_1(h) + h ∫_0^∞ e^{-t/h} F(t,h) dt` (Laplace form `Q02.proof_12`, the term `n = 1`
  split off, the substitution `t = n y` in the others), `F(t,h) = Σ_{n≥2} binom(h,n)/n H_n(t/n,h)`;
* `I_1`: Taylor's theorem twice for `𝓐` (as `Q02.W_taylor`), then Watson's lemma (`watson`);
* `F`: `F(t,h) = Σ_r h^r φ_r(t)` with `φ_r` holomorphic on `|t| < 6.2`, `|φ_r| ≤ C 2^r`, Cauchy's
  estimate on `|t| = 6`, then Watson's lemma;
* uniqueness of asymptotic expansions against `Q02.proof_6`/`Q02.proof_7`.
-/

open MeasureTheory Set Real Filter Topology Finset
open scoped ContDiff

namespace Q06.PA

/-! ### Uniqueness of asymptotic expansions -/

lemma poly_coeff_zero_aux (p : Polynomial ℝ) (D : ℕ) (C δ : ℝ) (hδ : 0 < δ)
    (h : ∀ x : ℝ, 0 < x → x ≤ δ → |p.eval x| ≤ C * x ^ (D + 1)) : p.coeff 0 = 0 := by
  have h1 : Tendsto (fun x => |p.eval x|) (𝓝[>] 0) (𝓝 |p.eval 0|) :=
    ((p.continuous.abs).tendsto 0).mono_left nhdsWithin_le_nhds
  have h2 : Tendsto (fun x : ℝ => C * x ^ (D + 1)) (𝓝[>] 0) (𝓝 (C * 0 ^ (D + 1))) :=
    (((continuous_const.mul (continuous_pow (D + 1)))).tendsto 0).mono_left nhdsWithin_le_nhds
  have hev : ∀ᶠ x in 𝓝[>] (0 : ℝ), |p.eval x| ≤ C * x ^ (D + 1) := by
    have : Ioo (0 : ℝ) δ ∈ 𝓝[>] (0 : ℝ) := Ioo_mem_nhdsGT hδ
    filter_upwards [this] with x hx using h x hx.1 hx.2.le
  have := le_of_tendsto_of_tendsto h1 h2 hev
  rw [zero_pow (by omega), mul_zero] at this
  rw [Polynomial.coeff_zero_eq_eval_zero]
  exact abs_nonpos_iff.mp this

lemma poly_coeff_zero (D : ℕ) : ∀ (p : Polynomial ℝ) (C δ : ℝ), 0 < δ →
    (∀ x : ℝ, 0 < x → x ≤ δ → |p.eval x| ≤ C * x ^ (D + 1)) → ∀ j ≤ D, p.coeff j = 0 := by
  induction D with
  | zero =>
    intro p C δ hδ h j hj
    rw [Nat.le_zero.mp hj]; exact poly_coeff_zero_aux p 0 C δ hδ h
  | succ D ih =>
    intro p C δ hδ h j hj
    have h0 := poly_coeff_zero_aux p (D + 1) C δ hδ h
    rcases j with _ | j
    · exact h0
    · rw [← Polynomial.coeff_divX]
      refine ih p.divX C δ hδ (fun x hx hxδ => ?_) j (by omega)
      have e := Polynomial.X_mul_divX_add p
      rw [h0, map_zero, add_zero] at e
      have := h x hx hxδ
      rw [← e, Polynomial.eval_mul, Polynomial.eval_X, abs_mul, abs_of_pos hx, pow_succ] at this
      have : x * |p.divX.eval x| ≤ x * (C * x ^ (D + 1)) := by linarith
      exact le_of_mul_le_mul_left this hx


/-! ### Watson's lemma, generic form -/

/-- the polynomial part `Σ_{p ∈ S} c_p t^{p.1} h^{p.2}` -/
noncomputable def Pp (S : Finset (ℕ × ℕ)) (c : ℕ × ℕ → ℝ) (t h : ℝ) : ℝ :=
  ∑ p ∈ S, c p * t ^ p.1 * h ^ p.2

lemma intOn_mono_exp (j : ℕ) {h : ℝ} (hh : 0 < h) :
    IntegrableOn (fun t => Real.exp (-(t / h)) * t ^ j) (Ioi 0) := by
  have := X02.intOn_pow_exp j (inv_pos.mpr hh)
  refine this.congr_fun (fun t _ => ?_) measurableSet_Ioi
  rw [mul_comm, div_eq_inv_mul]

lemma integral_mono_exp (j : ℕ) {h : ℝ} (hh : 0 < h) :
    ∫ t in Ioi (0 : ℝ), Real.exp (-(t / h)) * t ^ j = (j.factorial : ℝ) * h ^ (j + 1) := by
  have := X02.integral_pow_exp j (inv_pos.mpr hh)
  rw [show (fun t : ℝ => Real.exp (-(t / h)) * t ^ j) = fun t => t ^ j * Real.exp (-(h⁻¹ * t)) by
    funext t; rw [mul_comm, div_eq_inv_mul], this, inv_pow, div_inv_eq_mul]

lemma intOn_Pp (S : Finset (ℕ × ℕ)) (c : ℕ × ℕ → ℝ) {h : ℝ} (hh : 0 < h) :
    IntegrableOn (fun t => Real.exp (-(t / h)) * Pp S c t h) (Ioi 0) := by
  have : (fun t => Real.exp (-(t / h)) * Pp S c t h) =
      fun t => ∑ p ∈ S, (c p * h ^ p.2) * (Real.exp (-(t / h)) * t ^ p.1) := by
    funext t; rw [Pp, Finset.mul_sum]; congr 1; funext p; ring
  rw [this]
  exact integrable_finsetSum _ (fun p _ => (intOn_mono_exp p.1 hh).const_mul _)

lemma integral_Pp (S : Finset (ℕ × ℕ)) (c : ℕ × ℕ → ℝ) {h : ℝ} (hh : 0 < h) :
    ∫ t in Ioi (0 : ℝ), Real.exp (-(t / h)) * Pp S c t h =
      ∑ p ∈ S, c p * (p.1.factorial : ℝ) * h ^ (p.1 + p.2 + 1) := by
  have : (fun t => Real.exp (-(t / h)) * Pp S c t h) =
      fun t => ∑ p ∈ S, (c p * h ^ p.2) * (Real.exp (-(t / h)) * t ^ p.1) := by
    funext t; rw [Pp, Finset.mul_sum]; congr 1; funext p; ring
  rw [this, integral_finsetSum _ (fun p _ => (intOn_mono_exp p.1 hh).const_mul _)]
  refine Finset.sum_congr rfl (fun p _ => ?_)
  rw [integral_const_mul, integral_mono_exp p.1 hh]
  rw [show p.1 + p.2 + 1 = p.2 + (p.1 + 1) by ring, pow_add]; ring

lemma Pp_bd (S : Finset (ℕ × ℕ)) (c : ℕ × ℕ → ℝ) {t h : ℝ} (ht : 0 ≤ t) (hh : 0 ≤ h)
    (h1 : h ≤ 1) : |Pp S c t h| ≤ (∑ p ∈ S, |c p|) * (1 + t) ^ (S.sup Prod.fst) := by
  rw [Pp, Finset.sum_mul]
  refine (Finset.abs_sum_le_sum_abs _ _).trans (Finset.sum_le_sum (fun p hp => ?_))
  rw [abs_mul, abs_mul, abs_of_nonneg (pow_nonneg ht _), abs_of_nonneg (pow_nonneg hh _)]
  have e1 : t ^ p.1 ≤ (1 + t) ^ (S.sup Prod.fst) :=
    (pow_le_pow_left₀ ht (by linarith) _).trans
      (pow_le_pow_right₀ (by linarith) (Finset.le_sup (f := Prod.fst) hp))
  have e2 : h ^ p.2 ≤ 1 := pow_le_one₀ hh h1
  calc |c p| * t ^ p.1 * h ^ p.2 ≤ |c p| * (1 + t) ^ (S.sup Prod.fst) * 1 :=
        mul_le_mul (mul_le_mul_of_nonneg_left e1 (abs_nonneg _)) e2 (by positivity) (by positivity)
    _ = _ := mul_one _

/-- the majorant after the substitution `t = h u` -/
noncomputable def Psi (K D : ℕ) (a C B CP u : ℝ) : ℝ :=
  C * (u + 1) ^ (K + 1) + (u / a) ^ (K + 1) * (B * Real.exp (u / 2) + CP * (1 + u) ^ D)

lemma intOn_Psi (K D : ℕ) {a : ℝ} (ha : 0 < a) (C B CP : ℝ) :
    IntegrableOn (fun u => Real.exp (-u) * Psi K D a C B CP u) (Ioi 0) := by
  have i1 := (Q02.intOn_onepow (K + 1)).const_mul C
  have i2 := (X02.intOn_pow_exp (K + 1) (by norm_num : (0 : ℝ) < 1 / 2)).const_mul
    (B / a ^ (K + 1))
  have i3 : IntegrableOn (fun u : ℝ => u ^ (K + 1) * (1 + u) ^ D * Real.exp (-u)) (Ioi 0) := by
    refine X02.intOn_of_bound (by fun_prop) (Q02.intOn_onepow (K + 1 + D)) (fun u hu => ?_)
    rw [abs_of_nonneg (by positivity), pow_add (1 + u) (K + 1) D]
    have : u ^ (K + 1) ≤ (1 + u) ^ (K + 1) := pow_le_pow_left₀ hu.le (by linarith) _
    exact mul_le_mul_of_nonneg_right (mul_le_mul_of_nonneg_right this (by positivity))
      (Real.exp_pos _).le
  have i4 := i3.const_mul (CP / a ^ (K + 1))
  have : IntegrableOn (fun u => C * ((1 + u) ^ (K + 1) * Real.exp (-u)) +
      B / a ^ (K + 1) * (u ^ (K + 1) * Real.exp (-(1 / 2 * u))) +
      CP / a ^ (K + 1) * (u ^ (K + 1) * (1 + u) ^ D * Real.exp (-u))) (Ioi 0) :=
    (i1.add i2).add i4
  refine this.congr_fun (fun u _ => ?_) measurableSet_Ioi
  simp only [Psi]
  have e : Real.exp (-(1 / 2 * u)) = Real.exp (-u) * Real.exp (u / 2) := by
    rw [← Real.exp_add]; congr 1; ring
  rw [e, div_pow]
  field_simp
  ring

lemma pt_bd (S : Finset (ℕ × ℕ)) (c : ℕ × ℕ → ℝ) (Φ : ℝ → ℝ → ℝ) (K : ℕ) {a h0 C B : ℝ}
    (ha : 0 < a) (hh0' : h0 ≤ 1 / 2) (hC : 0 ≤ C) (hB : 0 ≤ B)
    (hT : ∀ h t, 0 < h → h ≤ h0 → 0 < t → t ≤ a →
      |Φ t h - Pp S c t h| ≤ C * (t + h) ^ (K + 1))
    (hE : ∀ h t, 0 < h → h ≤ h0 → 0 < t → |Φ t h| ≤ B * Real.exp t)
    {h t : ℝ} (hh : 0 < h) (hhh : h ≤ h0) (ht : 0 < t) :
    |Φ t h - Pp S c t h| ≤
      h ^ (K + 1) * Psi K (S.sup Prod.fst) a C B (∑ p ∈ S, |c p|) (t / h) := by
  set D := S.sup Prod.fst
  set CP := ∑ p ∈ S, |c p|
  have hCP : 0 ≤ CP := Finset.sum_nonneg (fun p _ => abs_nonneg _)
  set u := t / h with hu
  have hu0 : 0 < u := div_pos ht hh
  have htu : t = h * u := by rw [hu]; field_simp
  have hpsi1 : h ^ (K + 1) * (C * (u + 1) ^ (K + 1)) = C * (t + h) ^ (K + 1) := by
    rw [htu, show h * u + h = h * (u + 1) by ring, mul_pow]; ring
  have hsec : 0 ≤ (u / a) ^ (K + 1) * (B * Real.exp (u / 2) + CP * (1 + u) ^ D) := by positivity
  by_cases hta : t ≤ a
  · refine (hT h t hh hhh ht hta).trans ?_
    rw [Psi, mul_add, hpsi1]; nlinarith [pow_nonneg hh.le (K + 1)]
  · push Not at hta
    have h1 : |Φ t h - Pp S c t h| ≤ B * Real.exp t + CP * (1 + t) ^ D :=
      (abs_sub _ _).trans (add_le_add (hE h t hh hhh ht)
        (Pp_bd S c ht.le hh.le (by linarith)))
    have hh1 : h ≤ 1 / 2 := hhh.trans hh0'
    have htu2 : t ≤ u / 2 := by rw [htu]; nlinarith
    have htu1 : t ≤ u := by linarith
    have h2 : B * Real.exp t + CP * (1 + t) ^ D ≤ B * Real.exp (u / 2) + CP * (1 + u) ^ D := by
      gcongr
    have h3 : (1 : ℝ) ≤ h ^ (K + 1) * (u / a) ^ (K + 1) := by
      rw [← mul_pow, show h * (u / a) = t / a by rw [htu]; ring]
      exact one_le_pow₀ (by rw [le_div_iff₀ ha]; linarith)
    calc |Φ t h - Pp S c t h| ≤ 1 * (B * Real.exp (u / 2) + CP * (1 + u) ^ D) := by
          rw [one_mul]; exact h1.trans h2
      _ ≤ (h ^ (K + 1) * (u / a) ^ (K + 1)) * (B * Real.exp (u / 2) + CP * (1 + u) ^ D) :=
          mul_le_mul_of_nonneg_right h3 (by positivity)
      _ ≤ h ^ (K + 1) * Psi K D a C B CP u := by
          rw [Psi, mul_add]
          have : 0 ≤ h ^ (K + 1) * (C * (u + 1) ^ (K + 1)) := by positivity
          nlinarith

/-- **Watson's lemma** (generic): a two-variable expansion of `Φ` near `t = 0` and exponential
growth give the expansion of `∫_0^∞ e^{-t/h} Φ(t,h) dt` with error `O(h^{K+2})`. -/
theorem watson (S : Finset (ℕ × ℕ)) (c : ℕ × ℕ → ℝ) (Φ : ℝ → ℝ → ℝ) (K : ℕ) {a h0 C B : ℝ}
    (ha : 0 < a) (hh0' : h0 ≤ 1 / 2) (hC : 0 ≤ C) (hB : 0 ≤ B)
    (hT : ∀ h t, 0 < h → h ≤ h0 → 0 < t → t ≤ a →
      |Φ t h - Pp S c t h| ≤ C * (t + h) ^ (K + 1))
    (hE : ∀ h t, 0 < h → h ≤ h0 → 0 < t → |Φ t h| ≤ B * Real.exp t)
    (hcont : ∀ h, 0 < h → h ≤ h0 → ContinuousOn (fun t => Φ t h) (Ioi 0)) :
    ∃ C' : ℝ, ∀ h, 0 < h → h ≤ h0 →
      |(∫ t in Ioi (0 : ℝ), Real.exp (-(t / h)) * Φ t h) -
        ∑ p ∈ S, c p * (p.1.factorial : ℝ) * h ^ (p.1 + p.2 + 1)| ≤ C' * h ^ (K + 2) := by
  set D := S.sup Prod.fst
  set CP := ∑ p ∈ S, |c p|
  set g : ℝ → ℝ := fun u => Real.exp (-u) * Psi K D a C B CP u with hg
  refine ⟨∫ u in Ioi (0 : ℝ), g u, fun h hh hhh => ?_⟩
  have hgi : IntegrableOn g (Ioi 0) := intOn_Psi K D ha C B CP
  have hsc : IntegrableOn (fun t => g (h⁻¹ * t)) (Ioi 0) := by
    rw [integrableOn_Ioi_comp_mul_left_iff g 0 (inv_pos.mpr hh), mul_zero]; exact hgi
  have hbd : ∀ t, 0 < t → |Real.exp (-(t / h)) * (Φ t h - Pp S c t h)| ≤
      h ^ (K + 1) * g (h⁻¹ * t) := by
    intro t ht
    rw [abs_mul, abs_of_pos (Real.exp_pos _)]
    have := pt_bd S c Φ K ha hh0' hC hB hT hE hh hhh ht
    simp only [hg]
    rw [← div_eq_inv_mul]
    calc Real.exp (-(t / h)) * |Φ t h - Pp S c t h|
        ≤ Real.exp (-(t / h)) * (h ^ (K + 1) * Psi K D a C B CP (t / h)) :=
          mul_le_mul_of_nonneg_left this (Real.exp_pos _).le
      _ = _ := by ring
  have hdi : IntegrableOn (fun t => Real.exp (-(t / h)) * (Φ t h - Pp S c t h)) (Ioi 0) := by
    refine X02.intOn_of_boundOn ?_ (hsc.const_mul (h ^ (K + 1))) hbd
    have hPc : Continuous (fun t => Pp S c t h) := by unfold Pp; fun_prop
    exact (Continuous.continuousOn (by fun_prop)).mul ((hcont h hh hhh).sub hPc.continuousOn)
  have hsplit : (∫ t in Ioi (0 : ℝ), Real.exp (-(t / h)) * Φ t h) =
      (∫ t in Ioi (0 : ℝ), Real.exp (-(t / h)) * Pp S c t h) +
      ∫ t in Ioi (0 : ℝ), Real.exp (-(t / h)) * (Φ t h - Pp S c t h) := by
    rw [← integral_add (intOn_Pp S c hh) hdi]
    congr 1; funext t; ring
  rw [hsplit, integral_Pp S c hh, add_sub_cancel_left]
  have hn := norm_integral_le_of_norm_le (hsc.const_mul (h ^ (K + 1)))
    ((ae_restrict_iff' measurableSet_Ioi).2 (Filter.Eventually.of_forall
      (fun t ht => by rw [Real.norm_eq_abs]; exact hbd t ht)))
  rw [Real.norm_eq_abs] at hn
  refine hn.trans (le_of_eq ?_)
  rw [integral_const_mul, integral_comp_mul_left_Ioi g 0 (inv_pos.mpr hh), mul_zero, inv_inv,
    smul_eq_mul]
  ring


/-! ### The triangle `m + i ≤ K` -/

/-- `{(m,i) : m + i ≤ K}` -/
def Stri (K : ℕ) : Finset (ℕ × ℕ) :=
  (Finset.range (K + 1) ×ˢ Finset.range (K + 1)).filter (fun p => p.1 + p.2 ≤ K)

lemma tri_sum {M : Type*} [AddCommMonoid M] (K : ℕ) (g : ℕ → ℕ → M) :
    ∑ m ∈ Finset.range (K + 1), ∑ i ∈ Finset.range (K - m + 1), g m i =
      ∑ p ∈ Stri K, g p.1 p.2 := by
  rw [Stri, Finset.sum_filter, Finset.sum_product]
  refine Finset.sum_congr rfl (fun m hm => ?_)
  rw [← Finset.sum_filter]
  congr 1
  ext i; simp only [Finset.mem_range, Finset.mem_filter] at hm ⊢; omega

lemma tri_sum' {M : Type*} [AddCommMonoid M] (K : ℕ) (g : ℕ → ℕ → M) :
    ∑ r ∈ Finset.range (K + 1), ∑ m ∈ Finset.range (K - r + 1), g m r =
      ∑ p ∈ Stri K, g p.1 p.2 := by
  rw [tri_sum K (fun r m => g m r)]
  refine Finset.sum_nbij' (fun p => (p.2, p.1)) (fun p => (p.2, p.1)) ?_ ?_ ?_ ?_ ?_
  · intro p hp; simp only [Stri, Finset.mem_filter, Finset.mem_product,
      Finset.mem_range] at hp ⊢; omega
  · intro p hp; simp only [Stri, Finset.mem_filter, Finset.mem_product,
      Finset.mem_range] at hp ⊢; omega
  · intro p _; rfl
  · intro p _; rfl
  · intro p _; rfl

/-! ### `n = 1`: Taylor's theorem twice for `𝓐` (as `Q02.W_taylor`) -/

lemma Acal_contDiff : ContDiff ℝ ∞ (fun p : ℝ × ℝ => Q02.Acal p.2 p.1) := by
  have e : (fun p : ℝ × ℝ => Q02.Acal p.2 p.1) = fun p : ℝ × ℝ =>
      Real.exp (2 * p.2) * Real.exp (Real.log (1 + Real.exp p.2) * (-1 - p.1)) := by
    funext p; simp only [Q02.Acal]; rw [Real.rpow_def_of_pos (by positivity)]
  rw [e]
  have hl : ContDiff ℝ ∞ (fun p : ℝ × ℝ => Real.log (1 + Real.exp p.2)) :=
    (contDiff_const.add (contDiff_snd.exp)).log (fun p => by positivity)
  exact (contDiff_const.mul contDiff_snd).exp.mul ((hl.mul (contDiff_const.sub contDiff_fst)).exp)

/-- `GA m (h, y) = ∂_y^m 𝓐(y, h)` -/
noncomputable def GA (m : ℕ) (p : ℝ × ℝ) : ℝ := iteratedDeriv m (fun y => Q02.Acal y p.1) p.2

lemma GA_contDiff (m : ℕ) : ContDiff ℝ ∞ (GA m) := by
  induction m with
  | zero =>
    have : GA 0 = fun p : ℝ × ℝ => Q02.Acal p.2 p.1 := by funext q; simp [GA]
    rw [this]; exact Acal_contDiff
  | succ m ih =>
    have hf : ContDiff ℝ ∞ (Function.uncurry fun (p : ℝ × ℝ) (y : ℝ) => GA m (p.1, y)) :=
      ih.comp (contDiff_fst.fst.prodMk contDiff_snd)
    have := ContDiff.fderiv_apply (m := ∞) (n := ∞) hf contDiff_snd (contDiff_const (c := (1 : ℝ))) (by simp)
    convert this using 1
    funext p
    simp only [GA, iteratedDeriv_succ]
    rfl

lemma Acal_y_contDiffAt (h y : ℝ) : ContDiffAt ℝ ∞ (fun y => Q02.Acal y h) y :=
  (Acal_contDiff.comp (contDiff_const.prodMk contDiff_id)).contDiffAt

/-- Taylor's theorem in `y`, uniformly for `h ∈ [0, 1/2]`. -/
lemma Ataylor_y (K : ℕ) : ∃ C, ∀ h ∈ Icc (0 : ℝ) (1 / 2), ∀ y ∈ Icc (0 : ℝ) 1,
    |Q02.Acal y h - ∑ j ∈ Finset.range (K + 1), GA j (h, 0) / (j.factorial : ℝ) * y ^ j| ≤
      C * y ^ (K + 1) := by
  have hcont : ContinuousOn (GA (K + 1)) (Icc (0 : ℝ) (1 / 2) ×ˢ Icc (0 : ℝ) 1) :=
    (GA_contDiff (K + 1)).continuous.continuousOn
  obtain ⟨M, hM⟩ := (isCompact_Icc.prod isCompact_Icc).exists_bound_of_continuousOn hcont
  refine ⟨M / (K.factorial : ℝ), fun h hh y hy => ?_⟩
  have hcd : ContDiffOn ℝ ((K + 1 : ℕ) : ℕ∞ω) (fun y => Q02.Acal y h) (Icc 0 1) :=
    fun y' _ => (Q02.contDiffAt_le_top (Acal_y_contDiffAt h y') (K + 1)).contDiffWithinAt
  have hb : ∀ y' ∈ Icc (0 : ℝ) 1,
      ‖iteratedDerivWithin (K + 1) (fun y => Q02.Acal y h) (Icc 0 1) y'‖ ≤ M := by
    intro y' hy'
    rw [iteratedDerivWithin_eq_iteratedDeriv (uniqueDiffOn_Icc (by norm_num))
      (Q02.contDiffAt_le_top (Acal_y_contDiffAt h y') (K + 1)) hy']
    exact hM (h, y') ⟨hh, hy'⟩
  have T := taylor_mean_remainder_bound (by exact_mod_cast hcd) hy hb
  rw [taylor_within_apply, Real.norm_eq_abs, sub_zero] at T
  have e : ∑ k ∈ Finset.range (K + 1), ((k.factorial : ℝ)⁻¹ * y ^ k) •
        iteratedDerivWithin k (fun y => Q02.Acal y h) (Icc 0 1) 0 =
      ∑ j ∈ Finset.range (K + 1), GA j (h, 0) / (j.factorial : ℝ) * y ^ j := by
    refine Finset.sum_congr rfl (fun k _ => ?_)
    rw [iteratedDerivWithin_eq_iteratedDeriv (uniqueDiffOn_Icc (by norm_num))
      (Q02.contDiffAt_le_top (Acal_y_contDiffAt h 0) k) (left_mem_Icc.2 (by norm_num)),
      smul_eq_mul]
    simp only [GA]; ring
  rw [e] at T
  calc _ ≤ M * y ^ (K + 1) / (K.factorial : ℝ) := T
    _ = _ := by ring

lemma GAj_contDiffOn (j : ℕ) : ContDiffOn ℝ ∞ (fun x => GA j (x, 0)) (Iio 1) :=
  ((GA_contDiff j).comp (contDiff_id.prodMk contDiff_const)).contDiffOn

theorem A_taylor (K : ℕ) : ∃ C : ℝ, 0 ≤ C ∧ ∀ h y : ℝ, 0 < h → h ≤ 1 / 2 → 0 < y → y ≤ 1 →
    |Q02.Acal y h - Pp (Stri K) (fun p => Q02.Acoef p.1 p.2) y h| ≤ C * (y + h) ^ (K + 1) := by
  obtain ⟨C1, hC1⟩ := Ataylor_y K
  have hx := fun j : ℕ => Q02.taylor_x (GAj_contDiffOn j) (K - j)
  choose Cx hCx using hx
  refine ⟨|C1| + ∑ j ∈ Finset.range (K + 1), |Cx j|, by positivity,
    fun h y hh hh2 hy0 hy1 => ?_⟩
  have hhI : h ∈ Icc (0 : ℝ) (1 / 2) := ⟨hh.le, hh2⟩
  have hyI : y ∈ Icc (0 : ℝ) 1 := ⟨hy0.le, hy1⟩
  have hre : Pp (Stri K) (fun p => Q02.Acoef p.1 p.2) y h =
      ∑ j ∈ Finset.range (K + 1), y ^ j / (j.factorial : ℝ) *
        ∑ i ∈ Finset.range (K - j + 1),
          iteratedDeriv i (fun x => GA j (x, 0)) 0 / (i.factorial : ℝ) * h ^ i := by
    rw [Pp, ← tri_sum K (fun m i => Q02.Acoef m i * y ^ m * h ^ i)]
    refine Finset.sum_congr rfl (fun j _ => ?_)
    rw [Finset.mul_sum]
    refine Finset.sum_congr rfl (fun i _ => ?_)
    simp only [Q02.Acoef, GA]
    field_simp
  rw [hre, add_comm y h]
  set A := ∑ j ∈ Finset.range (K + 1), GA j (h, 0) / (j.factorial : ℝ) * y ^ j
  have key : Q02.Acal y h - ∑ j ∈ Finset.range (K + 1), y ^ j / (j.factorial : ℝ) *
        ∑ i ∈ Finset.range (K - j + 1),
          iteratedDeriv i (fun x => GA j (x, 0)) 0 / (i.factorial : ℝ) * h ^ i =
      (Q02.Acal y h - A) + ∑ j ∈ Finset.range (K + 1), y ^ j / (j.factorial : ℝ) *
        (GA j (h, 0) - ∑ i ∈ Finset.range (K - j + 1),
          iteratedDeriv i (fun x => GA j (x, 0)) 0 / (i.factorial : ℝ) * h ^ i) := by
    simp only [A, mul_sub, Finset.sum_sub_distrib]
    have : ∑ j ∈ Finset.range (K + 1), GA j (h, 0) / (j.factorial : ℝ) * y ^ j =
        ∑ j ∈ Finset.range (K + 1), y ^ j / (j.factorial : ℝ) * GA j (h, 0) := by
      refine Finset.sum_congr rfl (fun j _ => ?_); ring
    rw [this]; ring
  rw [key]
  have b1 : |Q02.Acal y h - A| ≤ |C1| * (h + y) ^ (K + 1) := by
    have := hC1 h hhI y hyI
    calc _ ≤ C1 * y ^ (K + 1) := this
      _ ≤ |C1| * y ^ (K + 1) := mul_le_mul_of_nonneg_right (le_abs_self _) (by positivity)
      _ ≤ |C1| * (h + y) ^ (K + 1) :=
          mul_le_mul_of_nonneg_left (pow_le_pow_left₀ hy0.le (by linarith) _) (abs_nonneg _)
  have b2 : ∀ j ∈ Finset.range (K + 1), |y ^ j / (j.factorial : ℝ) *
        (GA j (h, 0) - ∑ i ∈ Finset.range (K - j + 1),
          iteratedDeriv i (fun x => GA j (x, 0)) 0 / (i.factorial : ℝ) * h ^ i)|
        ≤ |Cx j| * (h + y) ^ (K + 1) := by
    intro j hj
    have hjK : j ≤ K := by simp at hj; omega
    have hf : (1 : ℝ) ≤ (j.factorial : ℝ) := by exact_mod_cast Nat.factorial_pos j
    have hyj : y ^ j / (j.factorial : ℝ) ≤ y ^ j := div_le_self (by positivity) hf
    have hT := hCx j h hhI
    rw [abs_mul, abs_of_nonneg (by positivity : (0 : ℝ) ≤ y ^ j / (j.factorial : ℝ))]
    calc y ^ j / (j.factorial : ℝ) * |_| ≤ y ^ j * (|Cx j| * h ^ (K - j + 1)) := by
          refine mul_le_mul hyj (hT.trans ?_) (abs_nonneg _) (by positivity)
          exact mul_le_mul_of_nonneg_right (le_abs_self _) (by positivity)
      _ = |Cx j| * (y ^ j * h ^ (K - j + 1)) := by ring
      _ ≤ |Cx j| * (h + y) ^ (K + 1) :=
          mul_le_mul_of_nonneg_left (Q02.mono_bd hh.le hy0.le hjK) (abs_nonneg _)
  calc _ ≤ |Q02.Acal y h - A| + |∑ j ∈ Finset.range (K + 1), y ^ j / (j.factorial : ℝ) *
        (GA j (h, 0) - ∑ i ∈ Finset.range (K - j + 1),
          iteratedDeriv i (fun x => GA j (x, 0)) 0 / (i.factorial : ℝ) * h ^ i)| := abs_add_le _ _
    _ ≤ |C1| * (h + y) ^ (K + 1) + ∑ j ∈ Finset.range (K + 1), |Cx j| * (h + y) ^ (K + 1) :=
        add_le_add b1 ((Finset.abs_sum_le_sum_abs _ _).trans (Finset.sum_le_sum b2))
    _ = _ := by rw [← Finset.sum_mul]; ring

lemma Acal_bd {t h : ℝ} (ht : 0 < t) (hh : 0 ≤ h) : |Q02.Acal t h| ≤ 1 * Real.exp t := by
  unfold Q02.Acal
  rw [abs_of_pos (by positivity), one_mul]
  have hle : (1 + Real.exp t) ^ (-1 - h) ≤ (Real.exp t) ^ (-1 - h) :=
    Real.rpow_le_rpow_of_nonpos (Real.exp_pos t) (by linarith [Real.exp_pos t]) (by linarith)
  rw [← Real.exp_mul] at hle
  calc Real.exp (2 * t) * (1 + Real.exp t) ^ (-1 - h) ≤ Real.exp (2 * t) * Real.exp (t * (-1 - h)) :=
        mul_le_mul_of_nonneg_left hle (Real.exp_pos _).le
    _ = Real.exp (t - t * h) := by rw [← Real.exp_add]; congr 1; ring
    _ ≤ Real.exp t := Real.exp_le_exp.mpr (by nlinarith)

lemma Acal_cont (h : ℝ) : Continuous (fun t => Q02.Acal t h) := by
  unfold Q02.Acal
  exact (by fun_prop : Continuous fun t : ℝ => Real.exp (2 * t)).mul
    ((by fun_prop : Continuous fun t : ℝ => 1 + Real.exp t).rpow_const
      (fun t => Or.inl (by positivity)))

/-- the term `n = 1`: `I_1(h) = ∫_0^∞ e^{-t/h} 𝓐(t,h) dt = Σ_{m+i ≤ K} m! [y^m h^i]𝓐 h^{m+i+1} + O(h^{K+2})` -/
theorem I1_expand (K : ℕ) : ∃ C' : ℝ, ∀ h, 0 < h → h ≤ 1 / 2 →
    |(∫ t in Ioi (0 : ℝ), Real.exp (-(t / h)) * Q02.Acal t h) -
      ∑ p ∈ Stri K, Q02.Acoef p.1 p.2 * (p.1.factorial : ℝ) * h ^ (p.1 + p.2 + 1)| ≤
        C' * h ^ (K + 2) := by
  obtain ⟨C, hC, hT⟩ := A_taylor K
  exact watson (Stri K) (fun p => Q02.Acoef p.1 p.2) (fun t h => Q02.Acal t h) K
    (a := 1) (h0 := 1 / 2) (C := C) (B := 1) one_pos le_rfl hC zero_le_one
    (fun h t hh hh2 ht ht1 => hT h t hh hh2 ht ht1)
    (fun h t hh _ ht => Acal_bd ht hh.le) (fun h _ _ => (Acal_cont h).continuousOn)

/-! ### The coefficients of `binom h n` as a polynomial in `h` -/

/-- `γ_{n,j} = [h^j] binom(h, n)` -/
noncomputable def gam (n j : ℕ) : ℝ := (descPochhammer ℝ n).coeff j / (n.factorial : ℝ)

/-- `P_n(ρ) = ρ(ρ+1)⋯(ρ+n-1)/n!` -/
noncomputable def Pn (n : ℕ) (ρ : ℝ) : ℝ := (∏ i ∈ Finset.range n, (ρ + i)) / (n.factorial : ℝ)

lemma desc_coeff_bd (n : ℕ) : ∀ (N : ℕ) {ρ : ℝ}, 0 ≤ ρ →
    ∑ j ∈ Finset.range N, |(descPochhammer ℝ n).coeff j| * ρ ^ j ≤
      ∏ i ∈ Finset.range n, (ρ + i) := by
  induction n with
  | zero =>
    intro N ρ hρ
    rw [descPochhammer_zero, Finset.prod_range_zero]
    rcases N with _ | N
    · simp
    · rw [Finset.sum_range_succ']
      have : ∀ j ∈ Finset.range N, |(1 : Polynomial ℝ).coeff (j + 1)| * ρ ^ (j + 1) = 0 := by
        intro j _; rw [Polynomial.coeff_one]; simp
      rw [Finset.sum_eq_zero this]; simp
  | succ n ih =>
    intro N ρ hρ
    have hP : 0 ≤ ∏ i ∈ Finset.range n, (ρ + i) :=
      Finset.prod_nonneg (fun i _ => by positivity)
    rw [Finset.prod_range_succ]
    set c := (descPochhammer ℝ n).coeff with hc
    have hsucc : descPochhammer ℝ (n + 1) = descPochhammer ℝ n * (Polynomial.X - Polynomial.C (n : ℝ)) := by
      rw [descPochhammer_succ_right, map_natCast]
    have e0 : (descPochhammer ℝ (n + 1)).coeff 0 = -(c 0 * n) := by
      rw [hsucc, Polynomial.mul_coeff_zero]; simp [hc]
    have e1 : ∀ j, (descPochhammer ℝ (n + 1)).coeff (j + 1) = c j - c (j + 1) * n := by
      intro j; rw [hsucc, Polynomial.coeff_mul_X_sub_C]
    rcases N with _ | M
    · simp only [Finset.sum_range_zero]; positivity
    · rw [Finset.sum_range_succ', e0]
      simp only [e1]
      have h1 : ∀ j ∈ Finset.range M, |c j - c (j + 1) * n| * ρ ^ (j + 1) ≤
          ρ * (|c j| * ρ ^ j) + n * (|c (j + 1)| * ρ ^ (j + 1)) := by
        intro j _
        have : |c j - c (j + 1) * n| ≤ |c j| + |c (j + 1)| * n := by
          refine (abs_sub _ _).trans ?_
          rw [abs_mul, Nat.abs_cast]
        calc |c j - c (j + 1) * n| * ρ ^ (j + 1) ≤ (|c j| + |c (j + 1)| * n) * ρ ^ (j + 1) :=
              mul_le_mul_of_nonneg_right this (by positivity)
          _ = _ := by rw [pow_succ]; ring
      have h2 := Finset.sum_le_sum h1
      rw [Finset.sum_add_distrib, ← Finset.mul_sum, ← Finset.mul_sum] at h2
      have hA := ih M hρ
      have hB := ih (M + 1) hρ
      rw [Finset.sum_range_succ'] at hB
      simp only [pow_zero, mul_one] at hB
      have hn : (0 : ℝ) ≤ n := Nat.cast_nonneg n
      have e2 : |-(c 0 * n)| * ρ ^ 0 = n * |c 0| := by
        rw [pow_zero, mul_one, abs_neg, abs_mul, Nat.abs_cast]; ring
      rw [e2]
      calc ∑ j ∈ Finset.range M, |c j - c (j + 1) * n| * ρ ^ (j + 1) + n * |c 0|
          ≤ ρ * ∑ j ∈ Finset.range M, |c j| * ρ ^ j +
              n * (∑ j ∈ Finset.range M, |c (j + 1)| * ρ ^ (j + 1) + |c 0|) := by nlinarith
        _ ≤ ρ * ∏ i ∈ Finset.range n, (ρ + i) + n * ∏ i ∈ Finset.range n, (ρ + i) := by
            gcongr
        _ = _ := by ring

lemma gam_sum_bd (n N : ℕ) {ρ : ℝ} (hρ : 0 ≤ ρ) :
    ∑ j ∈ Finset.range N, |gam n j| * ρ ^ j ≤ Pn n ρ := by
  have hf : (0 : ℝ) < (n.factorial : ℝ) := by exact_mod_cast Nat.factorial_pos n
  have := desc_coeff_bd n N hρ
  rw [Pn, le_div_iff₀ hf, Finset.sum_mul]
  refine le_of_eq_of_le (Finset.sum_congr rfl (fun j _ => ?_)) this
  rw [gam, abs_div, abs_of_pos hf]; field_simp

lemma gam_eq_zero {n j : ℕ} (hj : n < j) : gam n j = 0 := by
  rw [gam, Polynomial.coeff_eq_zero_of_natDegree_lt (by rw [descPochhammer_natDegree]; exact hj),
    zero_div]

lemma binom_eq_sum (h : ℝ) (n : ℕ) :
    Q02.binom h n = ∑ j ∈ Finset.range (n + 1), gam n j * h ^ j := by
  rw [Q02.binom, Polynomial.eval_eq_sum_range' (n := n + 1) (by rw [descPochhammer_natDegree]; omega),
    Finset.sum_div]
  refine Finset.sum_congr rfl (fun j _ => ?_)
  rw [gam]; ring

lemma Pn_succ (n : ℕ) (ρ : ℝ) : Pn (n + 1) ρ = Pn n ρ * ((ρ + n) / (n + 1)) := by
  rw [Pn, Pn, Finset.prod_range_succ, Nat.factorial_succ]
  push_cast
  have : (0 : ℝ) < (n.factorial : ℝ) := by exact_mod_cast Nat.factorial_pos n
  field_simp

lemma Pn_nonneg (n : ℕ) {ρ : ℝ} (hρ : 0 ≤ ρ) : 0 ≤ Pn n ρ := by
  unfold Pn; exact div_nonneg (Finset.prod_nonneg (fun i _ => by positivity)) (by positivity)

lemma Pn_half_sq (n : ℕ) : Pn n (1 / 2) ^ 2 ≤ 1 / ((n : ℝ) + 1) := by
  induction n with
  | zero => simp [Pn]
  | succ n ih =>
    rw [Pn_succ, mul_pow]
    have hn : (0 : ℝ) ≤ n := Nat.cast_nonneg n
    have h1 : ((1 / 2 + (n : ℝ)) / ((n : ℝ) + 1)) ^ 2 ≤ ((n : ℝ) + 1) / ((n : ℝ) + 2) := by
      rw [div_pow, div_le_div_iff₀ (by positivity) (by positivity)]
      nlinarith [sq_nonneg (n : ℝ)]
    calc Pn n (1 / 2) ^ 2 * ((1 / 2 + (n : ℝ)) / ((n : ℝ) + 1)) ^ 2
        ≤ 1 / ((n : ℝ) + 1) * (((n : ℝ) + 1) / ((n : ℝ) + 2)) :=
          mul_le_mul ih h1 (by positivity) (by positivity)
      _ = _ := by push_cast; field_simp; ring

lemma Pn_mono (n : ℕ) {ρ σ : ℝ} (hρ : 0 ≤ ρ) (hσ : ρ ≤ σ) : Pn n ρ ≤ Pn n σ := by
  unfold Pn
  gcongr

/-- `u_k = P_{k+2}(1/2)/(k+2)`, summable -/
noncomputable def uk (k : ℕ) : ℝ := Pn (k + 2) (1 / 2) / ((k : ℝ) + 2)

lemma uk_nonneg (k : ℕ) : 0 ≤ uk k := div_nonneg (Pn_nonneg _ (by norm_num)) (by positivity)

lemma uk_summable : Summable uk := by
  have hs : Summable (fun k : ℕ => 1 / ((k : ℝ) + 1) ^ (3 / 2 : ℝ)) := by
    have := (Real.summable_one_div_nat_rpow (p := 3 / 2)).2 (by norm_num)
    exact (summable_nat_add_iff 1).2 this |>.congr (fun n => by push_cast; ring_nf)
  refine Summable.of_nonneg_of_le uk_nonneg (fun k => ?_) hs
  have hk : (0 : ℝ) < (k : ℝ) + 1 := by positivity
  have hP := Pn_half_sq (k + 2)
  have hP0 := Pn_nonneg (k + 2) (by norm_num : (0 : ℝ) ≤ 1 / 2)
  have hsq : Pn (k + 2) (1 / 2) ≤ 1 / Real.sqrt ((k : ℝ) + 1) := by
    rw [le_div_iff₀ (Real.sqrt_pos.2 hk)]
    have h3 : (Pn (k + 2) (1 / 2) * Real.sqrt ((k : ℝ) + 1)) ^ 2 ≤ 1 := by
      rw [mul_pow, Real.sq_sqrt hk.le]
      calc Pn (k + 2) (1 / 2) ^ 2 * ((k : ℝ) + 1) ≤ 1 / (((k + 2 : ℕ) : ℝ) + 1) * ((k : ℝ) + 1) :=
            mul_le_mul_of_nonneg_right hP hk.le
        _ ≤ 1 := by
            rw [div_mul_eq_mul_div, one_mul, div_le_one (by positivity)]; push_cast; linarith
    nlinarith [Real.sqrt_nonneg ((k : ℝ) + 1), sq_nonneg (Pn (k + 2) (1 / 2) * Real.sqrt ((k : ℝ) + 1) - 1)]
  have hr : ((k : ℝ) + 1) ^ (3 / 2 : ℝ) = ((k : ℝ) + 1) * Real.sqrt ((k : ℝ) + 1) := by
    rw [show (3 / 2 : ℝ) = 1 + 1 / 2 by norm_num, Real.rpow_add hk, Real.rpow_one,
      Real.sqrt_eq_rpow]
  rw [hr, uk]
  calc Pn (k + 2) (1 / 2) / ((k : ℝ) + 2) ≤ (1 / Real.sqrt ((k : ℝ) + 1)) / ((k : ℝ) + 1) :=
        div_le_div₀ (by positivity) hsq hk (by linarith)
    _ = _ := by rw [div_div, mul_comm]

/-! ### `n ≥ 2`: the complex functions `φ_r` -/

/-- `λ(w) = Log(1 + e^w)` -/
noncomputable def lam (w : ℂ) : ℂ := Complex.log (1 + Complex.exp w)

lemma slit {w : ℂ} (hw : ‖w‖ < Real.pi) : 1 + Complex.exp w ∈ Complex.slitPlane := by
  rw [Complex.mem_slitPlane_iff]
  have him : |w.im| < Real.pi := (Complex.abs_im_le_norm w).trans_lt hw
  have him' := abs_lt.mp him
  by_cases h0 : w.im = 0
  · left
    simp only [Complex.add_re, Complex.one_re, Complex.exp_re, h0, Real.cos_zero, mul_one]
    positivity
  · right
    simp only [Complex.add_im, Complex.one_im, zero_add, Complex.exp_im]
    apply mul_ne_zero (Real.exp_pos _).ne'
    intro hs
    rcases lt_or_gt_of_ne h0 with hneg | hpos
    · have := Real.sin_pos_of_pos_of_lt_pi (x := -w.im) (by linarith) (by linarith)
      rw [Real.sin_neg, hs] at this; linarith
    · have := Real.sin_pos_of_pos_of_lt_pi hpos (by linarith)
      linarith

lemma lam_diffAt {w : ℂ} (hw : ‖w‖ < Real.pi) : DifferentiableAt ℂ lam w :=
  ((differentiableAt_const (1 : ℂ)).add Complex.differentiableAt_exp).clog (slit hw)

lemma r31 : (31 / 10 : ℝ) < Real.pi := by linarith [Real.pi_gt_d2]

lemma lam_bounds : ∃ L B0 : ℝ, 0 ≤ L ∧ 0 ≤ B0 ∧ ∀ w : ℂ, ‖w‖ ≤ 31 / 10 →
    ‖lam w‖ ≤ L ∧ ‖Complex.exp (-lam w)‖ ≤ B0 := by
  have hm : ∀ w ∈ Metric.closedBall (0 : ℂ) (31 / 10), ‖w‖ < Real.pi := fun w hw => by
    rw [Metric.mem_closedBall, dist_zero_right] at hw; linarith [r31]
  have hc : ContinuousOn lam (Metric.closedBall 0 (31 / 10)) := fun w hw =>
    (lam_diffAt (hm w hw)).continuousAt.continuousWithinAt
  have hc2 : ContinuousOn (fun w => Complex.exp (-lam w)) (Metric.closedBall 0 (31 / 10)) :=
    fun w hw => ((lam_diffAt (hm w hw)).neg.cexp).continuousAt.continuousWithinAt
  obtain ⟨L, hL⟩ := (isCompact_closedBall (0 : ℂ) (31 / 10)).exists_bound_of_continuousOn hc
  obtain ⟨B0, hB0⟩ := (isCompact_closedBall (0 : ℂ) (31 / 10)).exists_bound_of_continuousOn hc2
  refine ⟨max L 0, max B0 0, le_max_right _ _, le_max_right _ _, fun w hw => ?_⟩
  have hw' : w ∈ Metric.closedBall (0 : ℂ) (31 / 10) := by
    rw [Metric.mem_closedBall, dist_zero_right]; exact hw
  exact ⟨(hL w hw').trans (le_max_left _ _), (hB0 w hw').trans (le_max_left _ _)⟩

/-- `G_n(t) = e^{(n+1)t/n} e^{-λ(t/n)}` -/
noncomputable def Gc (n : ℕ) (t : ℂ) : ℂ :=
  Complex.exp ((((n : ℂ) + 1) / n) * t) * Complex.exp (-lam (t / n))

/-- `x_{n,i}(t) = G_n(t) (-λ(t/n))^i / i!` -/
noncomputable def xc (n i : ℕ) (t : ℂ) : ℂ := Gc n t * (-lam (t / n)) ^ i / (i.factorial : ℂ)

/-- `c_{n,r}(t) = (1/n) Σ_{j ≤ r} γ_{n,j} x_{n,r-j}(t)` -/
noncomputable def cc (n r : ℕ) (t : ℂ) : ℂ :=
  1 / (n : ℂ) * ∑ j ∈ Finset.range (r + 1), (gam n j : ℂ) * xc n (r - j) t

/-- `φ_r(t) = Σ_{n ≥ 2} c_{n,r}(t)` -/
noncomputable def phi (r : ℕ) (t : ℂ) : ℂ := ∑' k : ℕ, cc (k + 2) r t

lemma tn_bd {n : ℕ} (hn : 2 ≤ n) {t : ℂ} (ht : ‖t‖ < 62 / 10) : ‖t / (n : ℂ)‖ < 31 / 10 := by
  have hn' : (2 : ℝ) ≤ n := by exact_mod_cast hn
  rw [norm_div, Complex.norm_natCast, div_lt_iff₀ (by linarith)]
  nlinarith [norm_nonneg t]

lemma Gc_bd {L B0 : ℝ} (hLB : ∀ w : ℂ, ‖w‖ ≤ 31 / 10 → ‖lam w‖ ≤ L ∧ ‖Complex.exp (-lam w)‖ ≤ B0)
    {n : ℕ} (hn : 2 ≤ n) {t : ℂ} (ht : ‖t‖ < 62 / 10) :
    ‖Gc n t‖ ≤ Real.exp (93 / 10) * B0 := by
  have hn' : (2 : ℝ) ≤ n := by exact_mod_cast hn
  rw [Gc, norm_mul]
  refine mul_le_mul ?_ (hLB _ (tn_bd hn ht).le).2 (norm_nonneg _) (Real.exp_pos _).le
  refine (Complex.norm_exp_le_exp_norm _).trans (Real.exp_le_exp.mpr ?_)
  have e : (((n : ℂ) + 1) / n) = ((((n : ℝ) + 1) / n : ℝ) : ℂ) := by push_cast; rfl
  rw [norm_mul, e, Complex.norm_real, Real.norm_eq_abs, abs_of_pos (by positivity)]
  have h1 : ((n : ℝ) + 1) / n ≤ 3 / 2 := by rw [div_le_iff₀ (by linarith)]; linarith
  nlinarith [norm_nonneg t]

lemma xc_bd {L B0 : ℝ} (hLB : ∀ w : ℂ, ‖w‖ ≤ 31 / 10 → ‖lam w‖ ≤ L ∧ ‖Complex.exp (-lam w)‖ ≤ B0)
    (hL : 0 ≤ L) {n : ℕ} (hn : 2 ≤ n) (i : ℕ) {t : ℂ} (ht : ‖t‖ < 62 / 10) :
    ‖xc n i t‖ ≤ Real.exp (93 / 10) * B0 * L ^ i / (i.factorial : ℝ) := by
  rw [xc, norm_div, norm_mul, norm_pow, norm_neg, Complex.norm_natCast]
  have h1 := Gc_bd hLB hn ht
  have h2 := (hLB _ (tn_bd hn ht).le).1
  have h3 : 0 ≤ Real.exp (93 / 10) * B0 := (norm_nonneg _).trans h1
  gcongr

lemma cc_bd {L B0 : ℝ} (hLB : ∀ w : ℂ, ‖w‖ ≤ 31 / 10 → ‖lam w‖ ≤ L ∧ ‖Complex.exp (-lam w)‖ ≤ B0)
    (hL : 0 ≤ L) (hB0 : 0 ≤ B0) (k r : ℕ) {t : ℂ} (ht : ‖t‖ < 62 / 10) :
    ‖cc (k + 2) r t‖ ≤ Real.exp (93 / 10) * B0 * Real.exp (L / 2) * 2 ^ r * uk k := by
  set BG := Real.exp (93 / 10) * B0
  have hBG : 0 ≤ BG := by positivity
  have hterm : ∀ j ∈ Finset.range (r + 1), ‖(gam (k + 2) j : ℂ) * xc (k + 2) (r - j) t‖ ≤
      BG * Real.exp (L / 2) * 2 ^ r * (|gam (k + 2) j| * (1 / 2) ^ j) := by
    intro j hj
    have hjr : j ≤ r := by simp at hj; omega
    rw [norm_mul, Complex.norm_real, Real.norm_eq_abs]
    have hx := xc_bd hLB hL (by omega : 2 ≤ k + 2) (r - j) ht
    have hf : (L / 2) ^ (r - j) / ((r - j).factorial : ℝ) ≤ Real.exp (L / 2) :=
      Real.pow_div_factorial_le_exp (L / 2) (by positivity) _
    have e2 : (2 : ℝ) ^ r = 2 ^ (r - j) * 2 ^ j := by rw [← pow_add, Nat.sub_add_cancel hjr]
    have hx2 : BG * L ^ (r - j) / ((r - j).factorial : ℝ) ≤ BG * Real.exp (L / 2) * 2 ^ (r - j) := by
      have e3 : BG * L ^ (r - j) / ((r - j).factorial : ℝ) =
          BG * 2 ^ (r - j) * ((L / 2) ^ (r - j) / ((r - j).factorial : ℝ)) := by
        rw [div_pow]; field_simp
      rw [e3]
      calc BG * 2 ^ (r - j) * ((L / 2) ^ (r - j) / ((r - j).factorial : ℝ)) ≤
          BG * 2 ^ (r - j) * Real.exp (L / 2) := mul_le_mul_of_nonneg_left hf (by positivity)
        _ = _ := by ring
    calc |gam (k + 2) j| * ‖xc (k + 2) (r - j) t‖ ≤
        |gam (k + 2) j| * (BG * Real.exp (L / 2) * 2 ^ (r - j)) :=
          mul_le_mul_of_nonneg_left (hx.trans hx2) (abs_nonneg _)
      _ = _ := by rw [e2, one_div, inv_pow]; field_simp
  rw [cc, norm_mul]
  have hs := (norm_sum_le _ _).trans (Finset.sum_le_sum hterm)
  rw [← Finset.mul_sum] at hs
  have hg := gam_sum_bd (k + 2) (r + 1) (by norm_num : (0 : ℝ) ≤ 1 / 2)
  have e1 : ‖1 / ((k + 2 : ℕ) : ℂ)‖ = 1 / ((k : ℝ) + 2) := by
    rw [norm_div, norm_one, Complex.norm_natCast]; push_cast; ring
  rw [e1, uk]
  calc 1 / ((k : ℝ) + 2) * ‖∑ j ∈ Finset.range (r + 1), (gam (k + 2) j : ℂ) * xc (k + 2) (r - j) t‖
      ≤ 1 / ((k : ℝ) + 2) * (BG * Real.exp (L / 2) * 2 ^ r * Pn (k + 2) (1 / 2)) := by
        gcongr
        exact hs.trans (mul_le_mul_of_nonneg_left hg (by positivity))
    _ = _ := by ring

lemma cc_diffAt (k r : ℕ) {t : ℂ} (ht : ‖t‖ < 62 / 10) : DifferentiableAt ℂ (cc (k + 2) r) t := by
  have hl : DifferentiableAt ℂ (fun t : ℂ => lam (t / ((k + 2 : ℕ) : ℂ))) t :=
    (lam_diffAt (by linarith [tn_bd (by omega : 2 ≤ k + 2) ht, r31])).comp t
      (differentiableAt_id.div_const _)
  have hG : DifferentiableAt ℂ (Gc (k + 2)) t :=
    (((differentiableAt_const _).mul differentiableAt_id).cexp).mul hl.neg.cexp
  have hx : ∀ i, DifferentiableAt ℂ (xc (k + 2) i) t := fun i =>
    (hG.mul (hl.neg.pow i)).div_const _
  unfold cc
  exact (differentiableAt_const _).mul (DifferentiableAt.fun_sum (fun j _ =>
    (differentiableAt_const _).mul (hx (r - j))))

lemma phi_diff (r : ℕ) : DifferentiableOn ℂ (phi r) (Metric.ball 0 (62 / 10)) := by
  obtain ⟨L, B0, hL, hB0, hLB⟩ := lam_bounds
  refine Complex.differentiableOn_tsum_of_summable_norm
    (uk_summable.mul_left (Real.exp (93 / 10) * B0 * Real.exp (L / 2) * 2 ^ r))
    (fun k => fun t ht => (cc_diffAt k r ?_).differentiableWithinAt) Metric.isOpen_ball
    (fun k t ht => cc_bd hLB hL hB0 k r ?_)
  · rwa [Metric.mem_ball, dist_zero_right] at ht
  · rwa [Metric.mem_ball, dist_zero_right] at ht

lemma phi_bd : ∃ Cφ : ℝ, 0 ≤ Cφ ∧ ∀ (r : ℕ) (t : ℂ), ‖t‖ < 62 / 10 → ‖phi r t‖ ≤ Cφ * 2 ^ r := by
  obtain ⟨L, B0, hL, hB0, hLB⟩ := lam_bounds
  set A := Real.exp (93 / 10) * B0 * Real.exp (L / 2)
  have hA : 0 ≤ A := by positivity
  refine ⟨A * ∑' k, uk k, mul_nonneg hA (tsum_nonneg uk_nonneg), fun r t ht => ?_⟩
  have hb := fun k => cc_bd hLB hL hB0 k r ht
  have hs : Summable (fun k => A * 2 ^ r * uk k) := uk_summable.mul_left _
  have hn : Summable (fun k => ‖cc (k + 2) r t‖) :=
    Summable.of_nonneg_of_le (fun k => norm_nonneg _) hb hs
  calc ‖phi r t‖ ≤ ∑' k, ‖cc (k + 2) r t‖ := norm_tsum_le_tsum_norm hn
    _ ≤ ∑' k, A * 2 ^ r * uk k := hn.tsum_le_tsum hb hs
    _ = _ := by rw [tsum_mul_left]; ring

/-- the Taylor coefficient `[t^m] φ_r` -/
noncomputable def co (m r : ℕ) : ℂ := iteratedDeriv m (phi r) 0 / (m.factorial : ℂ)

lemma co_bd : ∃ Cφ : ℝ, 0 ≤ Cφ ∧
    (∀ (r : ℕ) (t : ℂ), ‖t‖ < 62 / 10 → ‖phi r t‖ ≤ Cφ * 2 ^ r) ∧
    ∀ m r : ℕ, ‖co m r‖ ≤ Cφ * 2 ^ r / 6 ^ m := by
  obtain ⟨Cφ, hC0, hC⟩ := phi_bd
  refine ⟨Cφ, hC0, hC, fun m r => ?_⟩
  have hsub : Metric.closedBall (0 : ℂ) 6 ⊆ Metric.ball 0 (62 / 10) :=
    Metric.closedBall_subset_ball (by norm_num)
  have hdc : DiffContOnCl ℂ (phi r) (Metric.ball 0 6) := (phi_diff r).diffContOnCl_ball hsub
  have hcau := Complex.norm_iteratedDeriv_le_of_forall_mem_sphere_norm_le m (by norm_num) hdc
    (fun z hz => hC r z (by rw [mem_sphere_zero_iff_norm] at hz; rw [hz]; norm_num))
  have hf : (0 : ℝ) < (m.factorial : ℝ) := by exact_mod_cast Nat.factorial_pos m
  rw [co, norm_div, Complex.norm_natCast, div_le_iff₀ hf]
  calc ‖iteratedDeriv m (phi r) 0‖ ≤ m.factorial * (Cφ * 2 ^ r) / 6 ^ m := hcau
    _ = _ := by ring

lemma phi_taylor (r : ℕ) {t : ℂ} (ht : ‖t‖ < 62 / 10) :
    HasSum (fun m => co m r * t ^ m) (phi r t) := by
  have := Complex.hasSum_taylorSeries_on_ball (phi_diff r)
    (by rw [Metric.mem_ball, dist_zero_right]; exact ht)
  convert this using 1
  funext m
  rw [co, sub_zero, smul_eq_mul, smul_eq_mul]
  field_simp

/-! ### `F(t,h) = Σ_r h^r φ_r(t)` -/

lemma cexp_hasSum (z : ℂ) : HasSum (fun i : ℕ => z ^ i / (i.factorial : ℂ)) (Complex.exp z) := by
  rw [Complex.exp_eq_exp_ℂ]; exact NormedSpace.expSeries_div_hasSum_exp z

lemma Hn_complex (n : ℕ) (hn : 1 ≤ n) (t h : ℝ) :
    ((Q02.Hn n (t / n) h : ℝ) : ℂ) = Gc n (t : ℂ) * Complex.exp (-((h : ℂ) * lam ((t : ℂ) / n))) := by
  have hn' : (n : ℝ) ≠ 0 := Nat.cast_ne_zero.mpr (by omega)
  have hl : lam ((t : ℂ) / n) = ((Real.log (1 + Real.exp (t / n)) : ℝ) : ℂ) := by
    rw [lam, Complex.ofReal_log (by positivity)]
    push_cast; rfl
  rw [Gc, hl, Q02.Hn, Real.rpow_def_of_pos (by positivity)]
  rw [mul_assoc, ← Complex.exp_add]
  have e1 : (((n : ℂ) + 1) / n) * (t : ℂ) = (((((n : ℝ) + 1) * (t / n)) : ℝ) : ℂ) := by
    push_cast; field_simp
  rw [e1]
  have e2 : -(((Real.log (1 + Real.exp (t / n)) : ℝ) : ℂ)) +
      -((h : ℂ) * ((Real.log (1 + Real.exp (t / n)) : ℝ) : ℂ)) =
      ((Real.log (1 + Real.exp (t / n)) * (-1 - h) : ℝ) : ℂ) := by push_cast; ring
  rw [e2, ← Complex.ofReal_exp, ← Complex.ofReal_exp]
  push_cast; rfl

lemma Gexp_hasSum {L B0 : ℝ} (hLB : ∀ w : ℂ, ‖w‖ ≤ 31 / 10 → ‖lam w‖ ≤ L ∧ ‖Complex.exp (-lam w)‖ ≤ B0)
    {n : ℕ} (hn : 2 ≤ n) {t : ℂ} (ht : ‖t‖ < 62 / 10) (h : ℂ) :
    HasSum (fun i => h ^ i * xc n i t) (Gc n t * Complex.exp (-(h * lam (t / n)))) := by
  have := (cexp_hasSum (-(h * lam (t / n)))).mul_left (Gc n t)
  convert this using 1
  funext i
  rw [xc, show -(h * lam (t / n)) = h * (-lam (t / n)) by ring, mul_pow]
  ring

lemma g_norm_summable {L B0 : ℝ}
    (hLB : ∀ w : ℂ, ‖w‖ ≤ 31 / 10 → ‖lam w‖ ≤ L ∧ ‖Complex.exp (-lam w)‖ ≤ B0) (hL : 0 ≤ L)
    {n : ℕ} (hn : 2 ≤ n) {t : ℂ} (ht : ‖t‖ < 62 / 10) (h : ℂ) :
    Summable (fun i => ‖h ^ i * xc n i t‖) := by
  refine Summable.of_nonneg_of_le (fun i => norm_nonneg _) (fun i => ?_)
    ((Real.summable_pow_div_factorial (‖h‖ * L)).mul_left (Real.exp (93 / 10) * B0))
  rw [norm_mul, norm_pow]
  have := xc_bd hLB hL hn i ht
  calc ‖h‖ ^ i * ‖xc n i t‖ ≤ ‖h‖ ^ i * (Real.exp (93 / 10) * B0 * L ^ i / (i.factorial : ℝ)) :=
        mul_le_mul_of_nonneg_left this (by positivity)
    _ = _ := by rw [mul_pow]; ring

/-- one term `n` of `F`: `binom(h,n)/n · H_n(t/n,h) = Σ_r h^r c_{n,r}(t)` -/
lemma term_hasSum {L B0 : ℝ}
    (hLB : ∀ w : ℂ, ‖w‖ ≤ 31 / 10 → ‖lam w‖ ≤ L ∧ ‖Complex.exp (-lam w)‖ ≤ B0) (hL : 0 ≤ L)
    (k : ℕ) {t h : ℝ} (ht : |t| < 62 / 10) :
    HasSum (fun r => (h : ℂ) ^ r * cc (k + 2) r (t : ℂ))
      (((Q02.binom h (k + 2) / ((k + 2 : ℕ) : ℝ) *
        Q02.Hn (k + 2) (t / ((k + 2 : ℕ) : ℝ)) h : ℝ) : ℂ)) := by
  have ht' : ‖(t : ℂ)‖ < 62 / 10 := by rw [Complex.norm_real, Real.norm_eq_abs]; exact ht
  set n := k + 2 with hn
  have hn2 : 2 ≤ n := by omega
  set f : ℕ → ℂ := fun j => (gam n j : ℂ) * (h : ℂ) ^ j with hf
  set g : ℕ → ℂ := fun i => (h : ℂ) ^ i * xc n i (t : ℂ) with hg
  have hf0 : ∀ j ∉ Finset.range (n + 1), f j = 0 := by
    intro j hj
    simp only [Finset.mem_range, not_lt] at hj
    simp [hf, gam_eq_zero (by omega : n < j)]
  have hfs : Summable (fun j => ‖f j‖) :=
    summable_of_ne_finset_zero (s := Finset.range (n + 1)) (fun j hj => by rw [hf0 j hj, norm_zero])
  have hgs := g_norm_summable hLB hL hn2 ht' (h : ℂ)
  have hft : ∑' j, f j = ((Q02.binom h n : ℝ) : ℂ) := by
    rw [tsum_eq_sum (s := Finset.range (n + 1)) hf0, binom_eq_sum]; push_cast; rfl
  have hgt : ∑' i, g i = Gc n (t : ℂ) * Complex.exp (-((h : ℂ) * lam ((t : ℂ) / n))) :=
    (Gexp_hasSum hLB hn2 ht' (h : ℂ)).tsum_eq
  have H := (hasSum_sum_range_mul_of_summable_norm hfs hgs).mul_left (1 / (n : ℂ))
  rw [hft, hgt] at H
  convert H using 1
  · funext r
    rw [cc, Finset.mul_sum, Finset.mul_sum, Finset.mul_sum]
    refine Finset.sum_congr rfl (fun j hj => ?_)
    have hjr : j ≤ r := by simp at hj; omega
    simp only [hf, hg]
    rw [show r = j + (r - j) by omega, pow_add]
    simp only [show j + (r - j) - j = r - j by omega]
    ring
  · rw [show (((k + 2 : ℕ) : ℝ)) = (n : ℝ) from rfl]
    have := Hn_complex n (by omega) t h
    push_cast at this ⊢
    rw [this]; ring

/-- `F(t,h) = Σ_{n≥2} binom(h,n)/n · H_n(t/n,h)` -/
noncomputable def Fr (t h : ℝ) : ℝ :=
  ∑' k : ℕ, Q02.binom h (k + 2) / ((k + 2 : ℕ) : ℝ) * Q02.Hn (k + 2) (t / ((k + 2 : ℕ) : ℝ)) h

lemma Fr_eq {t h : ℝ} (ht : |t| < 62 / 10) (hh : 0 ≤ h) (hh4 : h ≤ 1 / 4) :
    ((Fr t h : ℝ) : ℂ) = ∑' r, (h : ℂ) ^ r * phi r (t : ℂ) := by
  obtain ⟨L, B0, hL, hB0, hLB⟩ := lam_bounds
  have ht' : ‖(t : ℂ)‖ < 62 / 10 := by rw [Complex.norm_real, Real.norm_eq_abs]; exact ht
  set A := Real.exp (93 / 10) * B0 * Real.exp (L / 2)
  set F2 : ℕ → ℕ → ℂ := fun k r => (h : ℂ) ^ r * cc (k + 2) r (t : ℂ) with hF2
  have hb : ∀ p : ℕ × ℕ, ‖Function.uncurry F2 p‖ ≤ A * (uk p.1 * (1 / 2) ^ p.2) := by
    intro p
    simp only [Function.uncurry, hF2, norm_mul, norm_pow, Complex.norm_real, Real.norm_eq_abs,
      abs_of_nonneg hh]
    have h1 := cc_bd hLB hL hB0 p.1 p.2 ht'
    have h2 : h ^ p.2 ≤ (1 / 4) ^ p.2 := pow_le_pow_left₀ hh hh4 _
    calc h ^ p.2 * ‖cc (p.1 + 2) p.2 (t : ℂ)‖ ≤ (1 / 4) ^ p.2 * (A * 2 ^ p.2 * uk p.1) :=
          mul_le_mul h2 h1 (norm_nonneg _) (by positivity)
      _ = A * (uk p.1 * ((1 / 4) ^ p.2 * 2 ^ p.2)) := by ring
      _ = _ := by rw [← mul_pow]; norm_num
  have hgs : Summable (fun p : ℕ × ℕ => A * (uk p.1 * (1 / 2) ^ p.2)) :=
    (Summable.mul_of_nonneg uk_summable (summable_geometric_of_lt_one (by norm_num) (by norm_num))
      (fun k => uk_nonneg k) (fun r => by positivity)).mul_left A
  have hsum : Summable (Function.uncurry F2) := Summable.of_norm_bounded hgs hb
  rw [Fr, Complex.ofReal_tsum]
  have e1 : ∀ k, (((Q02.binom h (k + 2) / ((k + 2 : ℕ) : ℝ) *
      Q02.Hn (k + 2) (t / ((k + 2 : ℕ) : ℝ)) h : ℝ) : ℂ)) = ∑' r, F2 k r :=
    fun k => (term_hasSum hLB hL k ht).tsum_eq.symm
  simp_rw [e1]
  rw [← hsum.tsum_comm]
  refine tsum_congr (fun r => ?_)
  simp only [hF2]
  rw [tsum_mul_left, phi]

lemma Fr_re {t h : ℝ} (ht : |t| < 62 / 10) (hh : 0 ≤ h) (hh4 : h ≤ 1 / 4) :
    Fr t h = (∑' r, (h : ℂ) ^ r * phi r (t : ℂ)).re := by
  rw [← Fr_eq ht hh hh4, Complex.ofReal_re]

lemma geom_tail_bd {x : ℝ} (hx0 : 0 ≤ x) (hx : x ≤ 1 / 2) (M : ℕ) {Cst : ℝ} (hC : 0 ≤ Cst)
    {a : ℕ → ℂ} (ha : ∀ m, ‖a m‖ ≤ Cst * x ^ m) :
    Summable a ∧ ‖∑' m, a (m + M)‖ ≤ 2 * Cst * x ^ M := by
  have hg : Summable (fun m : ℕ => Cst * x ^ m) :=
    (summable_geometric_of_lt_one hx0 (by linarith)).mul_left Cst
  have hs : Summable a := Summable.of_norm_bounded hg ha
  refine ⟨hs, ?_⟩
  have hg' : Summable (fun m : ℕ => Cst * x ^ M * x ^ m) :=
    (summable_geometric_of_lt_one hx0 (by linarith)).mul_left _
  have hb : ∀ m, ‖a (m + M)‖ ≤ Cst * x ^ M * x ^ m := fun m => by
    have := ha (m + M); rw [pow_add] at this; linarith [show Cst * (x ^ m * x ^ M) = Cst * x ^ M * x ^ m by ring]
  have hn : Summable (fun m => ‖a (m + M)‖) :=
    Summable.of_nonneg_of_le (fun m => norm_nonneg _) hb hg'
  calc ‖∑' m, a (m + M)‖ ≤ ∑' m, ‖a (m + M)‖ := norm_tsum_le_tsum_norm hn
    _ ≤ ∑' m, Cst * x ^ M * x ^ m := hn.tsum_le_tsum hb hg'
    _ = Cst * x ^ M * (1 - x)⁻¹ := by rw [tsum_mul_left, tsum_geometric_of_lt_one hx0 (by linarith)]
    _ ≤ Cst * x ^ M * 2 := by
        gcongr
        rw [inv_le_comm₀ (by linarith) (by norm_num)]; linarith
    _ = _ := by ring

lemma mono2 {t h : ℝ} (ht : 0 ≤ t) (hh : 0 ≤ h) {r K : ℕ} (hr : r ≤ K) :
    (2 * h) ^ r * t ^ (K - r + 1) ≤ (2 * (t + h)) ^ (K + 1) := by
  have e : (2 * (t + h)) ^ (K + 1) = (2 * (t + h)) ^ r * (2 * (t + h)) ^ (K - r + 1) := by
    rw [← pow_add]; congr 1; omega
  rw [e]
  exact mul_le_mul (pow_le_pow_left₀ (by positivity) (by linarith) _)
    (pow_le_pow_left₀ ht (by linarith) _) (by positivity) (by positivity)

/-- the two-variable expansion of `F` near `t = 0` -/
theorem F_taylor (K : ℕ) : ∃ C : ℝ, 0 ≤ C ∧ ∀ h t : ℝ, 0 < h → h ≤ 1 / 4 → 0 < t → t ≤ 1 →
    |Fr t h - Pp (Stri K) (fun p => (co p.1 p.2).re) t h| ≤ C * (t + h) ^ (K + 1) := by
  obtain ⟨Cφ, hφ0, hφ, hco⟩ := co_bd
  refine ⟨2 * Cφ * ((K : ℝ) + 2) * 2 ^ (K + 1), by positivity, fun h t hh hh4 ht ht1 => ?_⟩
  have htab : |t| < 62 / 10 := by rw [abs_of_pos ht]; linarith
  have hz : ‖(t : ℂ)‖ < 62 / 10 := by rw [Complex.norm_real, Real.norm_eq_abs]; exact htab
  -- the polynomial part, complex form
  have hP : Pp (Stri K) (fun p => (co p.1 p.2).re) t h =
      (∑ r ∈ Finset.range (K + 1), (h : ℂ) ^ r *
        ∑ m ∈ Finset.range (K - r + 1), co m r * (t : ℂ) ^ m).re := by
    rw [Pp, ← tri_sum' K (fun m r => (co m r).re * t ^ m * h ^ r), Complex.re_sum]
    refine Finset.sum_congr rfl (fun r _ => ?_)
    rw [Finset.mul_sum, Complex.re_sum]
    refine Finset.sum_congr rfl (fun m _ => ?_)
    rw [show (h : ℂ) ^ r * (co m r * (t : ℂ) ^ m) = co m r * ((t ^ m * h ^ r : ℝ) : ℂ) by
      push_cast; ring, Complex.re_mul_ofReal]
    ring
  rw [Fr_re htab hh.le hh4, hP, ← Complex.sub_re]
  refine (Complex.abs_re_le_norm _).trans ?_
  -- the tail in `h`
  have hrb : ∀ r, ‖(h : ℂ) ^ r * phi r (t : ℂ)‖ ≤ Cφ * (2 * h) ^ r := fun r => by
    rw [norm_mul, norm_pow, Complex.norm_real, Real.norm_eq_abs, abs_of_pos hh, mul_pow]
    calc h ^ r * ‖phi r (t : ℂ)‖ ≤ h ^ r * (Cφ * 2 ^ r) :=
          mul_le_mul_of_nonneg_left (hφ r _ hz) (by positivity)
      _ = _ := by ring
  obtain ⟨hsr, htail⟩ := geom_tail_bd (by positivity : 0 ≤ 2 * h) (by linarith) (K + 1) hφ0 hrb
  rw [← hsr.sum_add_tsum_nat_add (K + 1)]
  -- the Taylor tails in `t`
  have hTt : ∀ r ∈ Finset.range (K + 1), ‖(h : ℂ) ^ r * (phi r (t : ℂ) -
      ∑ m ∈ Finset.range (K - r + 1), co m r * (t : ℂ) ^ m)‖ ≤ 2 * Cφ * (2 * (t + h)) ^ (K + 1) := by
    intro r hr
    have hrK : r ≤ K := by simp at hr; omega
    have hmb : ∀ m, ‖co m r * (t : ℂ) ^ m‖ ≤ (Cφ * 2 ^ r) * (t / 6) ^ m := fun m => by
      rw [norm_mul, norm_pow, Complex.norm_real, Real.norm_eq_abs, abs_of_pos ht, div_pow]
      calc ‖co m r‖ * t ^ m ≤ Cφ * 2 ^ r / 6 ^ m * t ^ m :=
            mul_le_mul_of_nonneg_right (hco m r) (by positivity)
        _ = _ := by ring
    obtain ⟨hsm, hmt⟩ := geom_tail_bd (by positivity : 0 ≤ t / 6) (by linarith) (K - r + 1)
      (by positivity : 0 ≤ Cφ * 2 ^ r) hmb
    have e : phi r (t : ℂ) - ∑ m ∈ Finset.range (K - r + 1), co m r * (t : ℂ) ^ m =
        ∑' m, co (m + (K - r + 1)) r * (t : ℂ) ^ (m + (K - r + 1)) := by
      rw [← (phi_taylor r hz).tsum_eq, ← hsm.sum_add_tsum_nat_add (K - r + 1)]; ring
    rw [e, norm_mul, norm_pow, Complex.norm_real, Real.norm_eq_abs, abs_of_pos hh]
    have h6 : (t / 6) ^ (K - r + 1) ≤ t ^ (K - r + 1) :=
      pow_le_pow_left₀ (by positivity) (by linarith) _
    calc h ^ r * ‖∑' m, co (m + (K - r + 1)) r * (t : ℂ) ^ (m + (K - r + 1))‖
        ≤ h ^ r * (2 * (Cφ * 2 ^ r) * t ^ (K - r + 1)) :=
          mul_le_mul_of_nonneg_left (hmt.trans (by gcongr)) (by positivity)
      _ = 2 * Cφ * ((2 * h) ^ r * t ^ (K - r + 1)) := by rw [mul_pow]; ring
      _ ≤ 2 * Cφ * (2 * (t + h)) ^ (K + 1) :=
          mul_le_mul_of_nonneg_left (mono2 ht.le hh.le hrK) (by positivity)
  have e2 : ∑ r ∈ Finset.range (K + 1), (h : ℂ) ^ r * phi r (t : ℂ) +
      ∑' r, (h : ℂ) ^ (r + (K + 1)) * phi (r + (K + 1)) (t : ℂ) -
      ∑ r ∈ Finset.range (K + 1), (h : ℂ) ^ r *
        ∑ m ∈ Finset.range (K - r + 1), co m r * (t : ℂ) ^ m =
      ∑ r ∈ Finset.range (K + 1), (h : ℂ) ^ r * (phi r (t : ℂ) -
        ∑ m ∈ Finset.range (K - r + 1), co m r * (t : ℂ) ^ m) +
      ∑' r, (h : ℂ) ^ (r + (K + 1)) * phi (r + (K + 1)) (t : ℂ) := by
    simp only [mul_sub, Finset.sum_sub_distrib]; ring
  rw [e2]
  have hs1 := (norm_sum_le _ _).trans (Finset.sum_le_sum hTt)
  rw [Finset.sum_const, Finset.card_range, nsmul_eq_mul] at hs1
  have htail' : 2 * Cφ * (2 * h) ^ (K + 1) ≤ 2 * Cφ * (2 * (t + h)) ^ (K + 1) := by
    gcongr; linarith
  calc _ ≤ _ := norm_add_le _ _
    _ ≤ ((K + 1 : ℕ) : ℝ) * (2 * Cφ * (2 * (t + h)) ^ (K + 1)) + 2 * Cφ * (2 * (t + h)) ^ (K + 1) :=
        add_le_add hs1 (htail.trans htail')
    _ = _ := by rw [mul_pow]; push_cast; ring

lemma Hn_le {n : ℕ} {y h : ℝ} (hy : 0 ≤ y) (hh : 0 ≤ h) :
    0 < Q02.Hn n y h ∧ Q02.Hn n y h ≤ Real.exp ((n : ℝ) * y) := by
  unfold Q02.Hn
  refine ⟨by positivity, ?_⟩
  have hle : (1 + Real.exp y) ^ (-1 - h) ≤ (Real.exp y) ^ (-1 - h) :=
    Real.rpow_le_rpow_of_nonpos (Real.exp_pos y) (by linarith [Real.exp_pos y]) (by linarith)
  rw [← Real.exp_mul] at hle
  calc Real.exp (((n : ℝ) + 1) * y) * (1 + Real.exp y) ^ (-1 - h)
      ≤ Real.exp (((n : ℝ) + 1) * y) * Real.exp (y * (-1 - h)) :=
        mul_le_mul_of_nonneg_left hle (Real.exp_pos _).le
    _ = Real.exp ((n : ℝ) * y - h * y) := by rw [← Real.exp_add]; congr 1; ring
    _ ≤ _ := Real.exp_le_exp.mpr (by nlinarith)

lemma Fterm_bd {k : ℕ} {t h : ℝ} (ht : 0 ≤ t) (hh : 0 < h) (h1 : h < 1) :
    |Q02.binom h (k + 2) / ((k + 2 : ℕ) : ℝ) * Q02.Hn (k + 2) (t / ((k + 2 : ℕ) : ℝ)) h| ≤
      1 / ((k : ℝ) + 1) ^ 2 * Real.exp t := by
  have hk : (0 : ℝ) < ((k + 2 : ℕ) : ℝ) := by positivity
  obtain ⟨hH0, hH⟩ := Hn_le (n := k + 2) (y := t / ((k + 2 : ℕ) : ℝ)) (by positivity) hh.le
  rw [mul_div_cancel₀ _ hk.ne'] at hH
  have hb := Q02.binom_abs hh h1 (k + 1)
  rw [abs_mul, abs_div, abs_of_pos hk, abs_of_pos hH0]
  have e : ((k + 1 : ℕ) : ℝ) + 1 = ((k + 2 : ℕ) : ℝ) := by push_cast; ring
  rw [e] at hb
  calc |Q02.binom h (k + 2)| / ((k + 2 : ℕ) : ℝ) * Q02.Hn (k + 2) (t / ((k + 2 : ℕ) : ℝ)) h
      ≤ h / ((k + 2 : ℕ) : ℝ) / ((k + 2 : ℕ) : ℝ) * Real.exp t := by gcongr
    _ ≤ 1 / ((k : ℝ) + 1) ^ 2 * Real.exp t := by
        have key : h / ((k + 2 : ℕ) : ℝ) / ((k + 2 : ℕ) : ℝ) ≤ 1 / ((k : ℝ) + 1) ^ 2 := by
          rw [div_div, div_le_div_iff₀ (by positivity) (by positivity)]
          push_cast; nlinarith
        exact mul_le_mul_of_nonneg_right key (Real.exp_pos _).le

lemma Fr_bd {t h : ℝ} (ht : 0 ≤ t) (hh : 0 < h) (h1 : h < 1) :
    |Fr t h| ≤ (∑' k : ℕ, 1 / ((k : ℝ) + 1) ^ 2) * Real.exp t := by
  have hs : Summable (fun k : ℕ => 1 / ((k : ℝ) + 1) ^ 2 * Real.exp t) :=
    (X02.summable_inv_pow 2 le_rfl).mul_right _
  have hn : Summable (fun k : ℕ => ‖Q02.binom h (k + 2) / ((k + 2 : ℕ) : ℝ) *
      Q02.Hn (k + 2) (t / ((k + 2 : ℕ) : ℝ)) h‖) :=
    Summable.of_nonneg_of_le (fun k => norm_nonneg _)
      (fun k => by rw [Real.norm_eq_abs]; exact Fterm_bd ht hh h1) hs
  rw [Fr, ← Real.norm_eq_abs, ← tsum_mul_right]
  refine (norm_tsum_le_tsum_norm hn).trans (hn.tsum_le_tsum (fun k => ?_) hs)
  rw [Real.norm_eq_abs]; exact Fterm_bd ht hh h1

lemma Fr_cont {h : ℝ} (hh : 0 < h) (h1 : h < 1) : ContinuousOn (fun t => Fr t h) (Ioi 0) := by
  intro t0 ht0
  have hc : ContinuousOn (fun t => Fr t h) (Ioo 0 (t0 + 1)) := by
    refine continuousOn_tsum (u := fun k : ℕ => 1 / ((k : ℝ) + 1) ^ 2 * Real.exp (t0 + 1))
      (fun k => ?_) ((X02.summable_inv_pow 2 le_rfl).mul_right _) (fun k t ht => ?_)
    · apply Continuous.continuousOn
      unfold Q02.Hn
      have c1 : Continuous fun t : ℝ =>
          Real.exp ((((k + 2 : ℕ) : ℝ) + 1) * (t / ((k + 2 : ℕ) : ℝ))) := by fun_prop
      have c2 : Continuous fun t : ℝ => (1 + Real.exp (t / ((k + 2 : ℕ) : ℝ))) ^ (-1 - h) :=
        (by fun_prop : Continuous fun t : ℝ => 1 + Real.exp (t / ((k + 2 : ℕ) : ℝ))).rpow_const
          (fun t => Or.inl (by positivity))
      exact continuous_const.mul (c1.mul c2)
    · rw [Real.norm_eq_abs]
      refine (Fterm_bd ht.1.le hh h1).trans ?_
      gcongr; exact ht.2.le
  exact (hc.continuousAt (Ioo_mem_nhds ht0 (by linarith))).continuousWithinAt

/-- the terms `n ≥ 2`: `∫_0^∞ e^{-t/h} F(t,h) dt = Σ_{m+r ≤ K} m! [t^m]φ_r h^{m+r+1} + O(h^{K+2})` -/
theorem F_expand (K : ℕ) : ∃ C' : ℝ, ∀ h, 0 < h → h ≤ 1 / 4 →
    |(∫ t in Ioi (0 : ℝ), Real.exp (-(t / h)) * Fr t h) -
      ∑ p ∈ Stri K, (co p.1 p.2).re * (p.1.factorial : ℝ) * h ^ (p.1 + p.2 + 1)| ≤
        C' * h ^ (K + 2) := by
  obtain ⟨C, hC, hT⟩ := F_taylor K
  exact watson (Stri K) (fun p => (co p.1 p.2).re) Fr K
    (a := 1) (h0 := 1 / 4) (C := C) (B := ∑' k : ℕ, 1 / ((k : ℝ) + 1) ^ 2) one_pos (by norm_num) hC
    (tsum_nonneg (fun k => by positivity))
    (fun h t hh hh4 ht ht1 => hT h t hh hh4 ht ht1)
    (fun h t hh hh4 ht => Fr_bd ht.le hh (by linarith))
    (fun h hh hh4 => Fr_cont hh (by linarith))

/-! ### The splitting `Ehat = h² I_1 + h ∫ e^{-t/h} F` -/

/-- the term `n = k + 2` of `F` -/
noncomputable def Fk (k : ℕ) (h t : ℝ) : ℝ :=
  Q02.binom h (k + 2) / ((k + 2 : ℕ) : ℝ) * Q02.Hn (k + 2) (t / ((k + 2 : ℕ) : ℝ)) h

lemma Fk_cont (k : ℕ) (h : ℝ) : Continuous (fun t => Fk k h t) := by
  unfold Fk Q02.Hn
  have c1 : Continuous fun t : ℝ =>
      Real.exp ((((k + 2 : ℕ) : ℝ) + 1) * (t / ((k + 2 : ℕ) : ℝ))) := by fun_prop
  have c2 : Continuous fun t : ℝ => (1 + Real.exp (t / ((k + 2 : ℕ) : ℝ))) ^ (-1 - h) :=
    (by fun_prop : Continuous fun t : ℝ => 1 + Real.exp (t / ((k + 2 : ℕ) : ℝ))).rpow_const
      (fun t => Or.inl (by positivity))
  exact continuous_const.mul (c1.mul c2)

lemma Fkint_bd (k : ℕ) {h t : ℝ} (hh : 0 < h) (h1 : h < 1) (ht : 0 < t) :
    |Real.exp (-(t / h)) * Fk k h t| ≤ 1 / ((k : ℝ) + 1) ^ 2 * Real.exp (-((1 / h - 1) * t)) := by
  rw [abs_mul, abs_of_pos (Real.exp_pos _)]
  have := Fterm_bd (k := k) ht.le hh h1
  rw [← Fk] at this
  calc Real.exp (-(t / h)) * |Fk k h t| ≤ Real.exp (-(t / h)) * (1 / ((k : ℝ) + 1) ^ 2 * Real.exp t) :=
        mul_le_mul_of_nonneg_left this (Real.exp_pos _).le
    _ = 1 / ((k : ℝ) + 1) ^ 2 * (Real.exp (-(t / h)) * Real.exp t) := by ring
    _ = _ := by rw [← Real.exp_add]; congr 2; ring

theorem Ehat_split {h : ℝ} (hh : 0 < h) (h1 : h < 1) :
    Q02.Ehat h = h * h * (∫ t in Ioi (0 : ℝ), Real.exp (-(t / h)) * Q02.Acal t h) +
      h * ∫ t in Ioi (0 : ℝ), Real.exp (-(t / h)) * Fr t h := by
  have H := (Q02.proof_12 h hh h1).2
  rw [← H.tsum_eq, H.summable.tsum_eq_zero_add]
  have hc : 0 < 1 / h - 1 := by rw [sub_pos, lt_one_div one_pos hh]; simpa using h1
  congr 1
  · -- the term `n = 1`
    have hb : Q02.binom h (0 + 1) = h := by simp [Q02.binom]
    rw [hb]
    have e : ∀ y : ℝ, Real.exp (-(((0 + 1 : ℕ) : ℝ) * y / h)) * Q02.Hn (0 + 1) y h =
        Real.exp (-(y / h)) * Q02.Acal y h := fun y => by
      unfold Q02.Hn Q02.Acal; norm_num
    simp_rw [e]
  · -- the terms `n ≥ 2`
    have hint : ∀ k : ℕ, Integrable (fun t => Real.exp (-(t / h)) * Fk k h t)
        (volume.restrict (Ioi 0)) := fun k =>
      X02.intOn_of_bound ((by fun_prop : Continuous fun t : ℝ => Real.exp (-(t / h))).mul
        (Fk_cont k h)) ((Q02.intOn_exp_mul hc).const_mul _) (fun t ht => Fkint_bd k hh h1 ht)
    have hnorm : ∀ k : ℕ, ∫ t in Ioi (0 : ℝ), ‖Real.exp (-(t / h)) * Fk k h t‖ ≤
        1 / ((k : ℝ) + 1) ^ 2 * (1 / (1 / h - 1)) := fun k => by
      calc ∫ t in Ioi (0 : ℝ), ‖Real.exp (-(t / h)) * Fk k h t‖
          ≤ ∫ t in Ioi (0 : ℝ), 1 / ((k : ℝ) + 1) ^ 2 * Real.exp (-((1 / h - 1) * t)) := by
            refine setIntegral_mono_on (hint k).norm ((Q02.intOn_exp_mul hc).const_mul _)
              measurableSet_Ioi (fun t ht => ?_)
            rw [Real.norm_eq_abs]; exact Fkint_bd k hh h1 ht
        _ = _ := by rw [integral_const_mul, Q02.integral_exp_mul hc]
    have hsumN : Summable fun k : ℕ => ∫ t in Ioi (0 : ℝ), ‖Real.exp (-(t / h)) * Fk k h t‖ :=
      Summable.of_nonneg_of_le (fun k => integral_nonneg (fun t => norm_nonneg _)) hnorm
        ((X02.summable_inv_pow 2 le_rfl).mul_right _)
    have key := integral_tsum_of_summable_integral_norm hint hsumN
    have hF : ∀ t, ∑' k, Real.exp (-(t / h)) * Fk k h t = Real.exp (-(t / h)) * Fr t h := by
      intro t; rw [tsum_mul_left]; rfl
    simp_rw [hF] at key
    rw [← key, ← tsum_mul_left]
    refine tsum_congr (fun k => ?_)
    -- the substitution `t = (k+2) y`
    have hk : (0 : ℝ) < ((k + 2 : ℕ) : ℝ) := by positivity
    have hsub := integral_comp_mul_left_Ioi
      (fun t => Real.exp (-(t / h)) * Q02.Hn (k + 2) (t / ((k + 2 : ℕ) : ℝ)) h) 0 hk
    rw [mul_zero, smul_eq_mul] at hsub
    have e1 : ∀ y : ℝ, (fun t => Real.exp (-(t / h)) * Q02.Hn (k + 2) (t / ((k + 2 : ℕ) : ℝ)) h)
        (((k + 2 : ℕ) : ℝ) * y) = Real.exp (-(((k + 1 + 1 : ℕ) : ℝ) * y / h)) * Q02.Hn (k + 1 + 1) y h := by
      intro y; simp only; rw [mul_div_cancel_left₀ _ hk.ne', mul_div_assoc]
    simp_rw [e1] at hsub
    rw [hsub]
    have e2 : ∫ t in Ioi (0 : ℝ), Real.exp (-(t / h)) * Fk k h t =
        Q02.binom h (k + 2) / ((k + 2 : ℕ) : ℝ) *
          ∫ t in Ioi (0 : ℝ), Real.exp (-(t / h)) * Q02.Hn (k + 2) (t / ((k + 2 : ℕ) : ℝ)) h := by
      rw [← integral_const_mul]; congr 1; funext t; unfold Fk; ring
    rw [e2]
    field_simp

/-! ### The expansion of `Ehat` with the coefficients `N_k - (-log 2)^{k+1}/(k+1)!` -/

theorem Ehat_expand (K : ℕ) : ∃ C h₀ : ℝ, 0 < h₀ ∧ ∀ h : ℝ, 0 < h → h ≤ h₀ →
    |Q02.Ehat h - ∑ k ∈ Finset.range (K + 1), h ^ (k + 1) * ∑ n ∈ Finset.Ico 1 k, Q02.q n| ≤
      C * h ^ (K + 2) := by
  obtain ⟨C, h0, h0pos, hQ⟩ := Q02.proof_7 (K + 1) (by omega)
  set S := ∑ n ∈ Finset.Icc 1 K, Q02.q n
  refine ⟨2 * (|C| + |Q02.q (K + 1)| + |S|), min h0 (1 / 2), lt_min h0pos (by norm_num),
    fun h hh hhh => ?_⟩
  have hh0 : h ≤ h0 := hhh.trans (min_le_left _ _)
  have hh2 : h ≤ 1 / 2 := hhh.trans (min_le_right _ _)
  have h1h : 0 < 1 - h := by linarith
  have hE : Q02.Ehat h = h ^ 2 / (1 - h) * X02.Q h := (Q02.proof_6 h hh (by linarith)).2
  have hG : ∑ k ∈ Finset.range (K + 1), h ^ (k + 1) * ∑ n ∈ Finset.Ico 1 k, Q02.q n =
      h ^ 2 / (1 - h) * (∑ n ∈ Finset.Icc 1 K, Q02.q n * h ^ n - h ^ K * S) := by
    rw [div_mul_eq_mul_div, eq_div_iff h1h.ne', mul_comm]
    exact Q02.geom_id K h
  rw [hE, hG]
  set Tt := ∑ n ∈ Finset.Icc 1 K, Q02.q n * h ^ n with hT
  have hQK := hQ h hh hh0
  rw [Finset.sum_Icc_succ_top (by omega), ← hT] at hQK
  have hh1 : h ≤ 1 := by linarith
  have hhK : h ^ (K + 1) ≤ h ^ K := pow_le_pow_of_le_one hh.le hh1 (by omega)
  have hhK2 : h ^ (K + 2) ≤ h ^ K := pow_le_pow_of_le_one hh.le hh1 (by omega)
  have hB : |X02.Q h - Tt + h ^ K * S| ≤ (|C| + |Q02.q (K + 1)| + |S|) * h ^ K := by
    have e : X02.Q h - Tt + h ^ K * S = (X02.Q h - (Tt + Q02.q (K + 1) * h ^ (K + 1)))
        + Q02.q (K + 1) * h ^ (K + 1) + h ^ K * S := by ring
    rw [e]
    have t1 := abs_add_le (X02.Q h - (Tt + Q02.q (K + 1) * h ^ (K + 1)) + Q02.q (K + 1) * h ^ (K + 1))
      (h ^ K * S)
    have t2 := abs_add_le (X02.Q h - (Tt + Q02.q (K + 1) * h ^ (K + 1))) (Q02.q (K + 1) * h ^ (K + 1))
    have a1 : C * h ^ (K + 1 + 1) ≤ |C| * h ^ K :=
      mul_le_mul (le_abs_self C) hhK2 (by positivity) (abs_nonneg _)
    have a2 : |Q02.q (K + 1) * h ^ (K + 1)| ≤ |Q02.q (K + 1)| * h ^ K := by
      rw [abs_mul, abs_of_nonneg (pow_nonneg hh.le (K + 1))]
      exact mul_le_mul_of_nonneg_left hhK (abs_nonneg _)
    have a3 : |h ^ K * S| = |S| * h ^ K := by
      rw [abs_mul, abs_of_nonneg (pow_nonneg hh.le K), mul_comm]
    linarith
  have hc : h ^ 2 / (1 - h) ≤ 2 * h ^ 2 := by
    rw [div_le_iff₀ h1h]; nlinarith
  rw [← mul_sub, show X02.Q h - (Tt - h ^ K * S) = X02.Q h - Tt + h ^ K * S by ring, abs_mul,
    abs_of_nonneg (by positivity)]
  calc h ^ 2 / (1 - h) * |X02.Q h - Tt + h ^ K * S| ≤
      (2 * h ^ 2) * ((|C| + |Q02.q (K + 1)| + |S|) * h ^ K) :=
        mul_le_mul hc hB (abs_nonneg _) (by positivity)
    _ = _ := by rw [show K + 2 = 2 + K by ring, pow_add]; ring

/-! ### The arithmetic of `Σ_{m<k} m! 2^{k-1-m} 6^{-m}` -/

lemma fac12_le_one {m : ℕ} (hm : m ≤ 12) : (m.factorial : ℝ) / 12 ^ m ≤ 1 := by
  rw [div_le_one (by positivity)]
  have : m.factorial ≤ 12 ^ m := by interval_cases m <;> decide
  exact_mod_cast this

lemma fac12_mono {m n : ℕ} (hm : 12 ≤ m) (hmn : m ≤ n) :
    (m.factorial : ℝ) / 12 ^ m ≤ (n.factorial : ℝ) / 12 ^ n := by
  induction n with
  | zero => omega
  | succ n ih =>
    rcases Nat.lt_or_ge m (n + 1) with h | h
    · refine (ih (by omega)).trans ?_
      rw [Nat.factorial_succ, pow_succ]
      push_cast
      rw [div_le_div_iff₀ (by positivity) (by positivity)]
      have : (12 : ℝ) ≤ (n : ℝ) + 1 := by
        have : (12 : ℕ) ≤ n + 1 := by omega
        exact_mod_cast this
      have h0 : (0 : ℝ) ≤ (n.factorial : ℝ) * 12 ^ n := by positivity
      nlinarith
    · rw [show m = n + 1 by omega]

lemma fac12_bd {m n : ℕ} (hmn : m ≤ n) :
    (m.factorial : ℝ) / 12 ^ m ≤ (if m ≤ 12 then 1 else 0) + (n.factorial : ℝ) / 12 ^ n := by
  have h0 : (0 : ℝ) ≤ (n.factorial : ℝ) / 12 ^ n := by positivity
  split_ifs with h
  · linarith [fac12_le_one h]
  · rw [zero_add]; exact fac12_mono (by omega) hmn

lemma sum_fac12 (n : ℕ) :
    ∑ m ∈ Finset.range (n + 1), (m.factorial : ℝ) / 12 ^ m ≤
      13 + ((n : ℝ) + 1) * ((n.factorial : ℝ) / 12 ^ n) := by
  have h1 := Finset.sum_le_sum (fun m (hm : m ∈ Finset.range (n + 1)) =>
    fac12_bd (n := n) (m := m) (by simp at hm; omega))
  rw [Finset.sum_add_distrib, Finset.sum_const, Finset.card_range, nsmul_eq_mul] at h1
  have h2 : ∑ m ∈ Finset.range (n + 1), (if m ≤ 12 then (1 : ℝ) else 0) ≤ 13 := by
    rw [Finset.sum_boole]
    have : ((Finset.range (n + 1)).filter (fun m => m ≤ 12)) ⊆ Finset.range 13 := by
      intro m hm; simp at hm ⊢; omega
    have := Finset.card_le_card this
    rw [Finset.card_range] at this
    exact_mod_cast this
  push_cast at h1
  linarith

lemma U_arith (n : ℕ) :
    ∑ m ∈ Finset.range (n + 2), (m.factorial : ℝ) * 2 ^ (n + 1 - m) / 6 ^ m ≤
      78 * (2 ^ (n + 2) + ((n + 2 : ℕ) : ℝ) * (n.factorial : ℝ) / 6 ^ (n + 2)) := by
  rw [Finset.sum_range_succ]
  have hterm : ∀ m ∈ Finset.range (n + 1),
      (m.factorial : ℝ) * 2 ^ (n + 1 - m) / 6 ^ m = 2 ^ (n + 1) * ((m.factorial : ℝ) / 12 ^ m) := by
    intro m hm
    have hmn : m ≤ n + 1 := by simp at hm; omega
    have e : (2 : ℝ) ^ (n + 1) = 2 ^ (n + 1 - m) * 2 ^ m := by
      rw [← pow_add, Nat.sub_add_cancel hmn]
    rw [e, show (12 : ℝ) ^ m = 2 ^ m * 6 ^ m by rw [← mul_pow]; norm_num]
    field_simp
  rw [Finset.sum_congr rfl hterm, ← Finset.mul_sum]
  have hs := sum_fac12 n
  have hlast : ((n + 1).factorial : ℝ) * 2 ^ (n + 1 - (n + 1)) / 6 ^ (n + 1) =
      6 * ((n : ℝ) + 1) * (n.factorial : ℝ) / 6 ^ (n + 2) := by
    rw [Nat.sub_self, pow_zero, mul_one, Nat.factorial_succ, pow_succ 6 (n + 1)]
    push_cast; field_simp
  rw [hlast]
  have hmid : (2 : ℝ) ^ (n + 1) * (((n : ℝ) + 1) * ((n.factorial : ℝ) / 12 ^ n)) =
      72 * ((n : ℝ) + 1) * (n.factorial : ℝ) / 6 ^ (n + 2) := by
    rw [show (12 : ℝ) ^ n = 2 ^ n * 6 ^ n by rw [← mul_pow]; norm_num]
    rw [pow_succ, pow_succ 6 (n + 1), pow_succ 6 n]
    field_simp
    ring
  have hF : (0 : ℝ) ≤ (n.factorial : ℝ) / 6 ^ (n + 2) := by positivity
  have h2 : (0 : ℝ) ≤ 2 ^ (n + 1) := by positivity
  calc (2 : ℝ) ^ (n + 1) * ∑ m ∈ Finset.range (n + 1), (m.factorial : ℝ) / 12 ^ m +
        6 * ((n : ℝ) + 1) * (n.factorial : ℝ) / 6 ^ (n + 2)
      ≤ 2 ^ (n + 1) * (13 + ((n : ℝ) + 1) * ((n.factorial : ℝ) / 12 ^ n)) +
        6 * ((n : ℝ) + 1) * (n.factorial : ℝ) / 6 ^ (n + 2) := by gcongr
    _ = 13 * 2 ^ (n + 1) + 78 * (((n : ℝ) + 1) * ((n.factorial : ℝ) / 6 ^ (n + 2))) := by
        rw [mul_add, hmid]; ring
    _ ≤ 78 * 2 ^ (n + 2) + 78 * (((n + 2 : ℕ) : ℝ) * ((n.factorial : ℝ) / 6 ^ (n + 2))) := by
        have : (0 : ℝ) ≤ (n : ℝ) := Nat.cast_nonneg n
        have a : 13 * (2 : ℝ) ^ (n + 1) ≤ 78 * 2 ^ (n + 2) := by
          rw [pow_succ 2 (n + 1)]; nlinarith
        have b : ((n : ℝ) + 1) * ((n.factorial : ℝ) / 6 ^ (n + 2)) ≤
            (((n + 2 : ℕ) : ℝ)) * ((n.factorial : ℝ) / 6 ^ (n + 2)) := by
          push_cast; nlinarith
        linarith
    _ = _ := by ring

/-! ### Comparison of the two expansions -/

lemma filter_antidiag (K d s : ℕ) (hd : d ≤ K) (h : K + 1 = d + s) (f : ℕ × ℕ → ℝ) :
    ∑ p ∈ Stri K, (if K + 1 = p.1 + p.2 + s then f p else 0) =
      ∑ m ∈ Finset.range (d + 1), f (m, d - m) := by
  rw [← Finset.sum_filter, ← Finset.Nat.sum_antidiagonal_eq_sum_range_succ_mk]
  congr 1
  ext p
  simp only [Stri, Finset.mem_filter, Finset.mem_product, Finset.mem_range,
    Finset.HasAntidiagonal.mem_antidiagonal]
  omega

/-- the coefficient of `h^{k+1}`: `N_k - (-log 2)^{k+1}/(k+1)! = T_k + U_k` -/
theorem coeff_identity (n : ℕ) :
    ∑ i ∈ Finset.Ico 1 (n + 2), Q02.q i =
      ∑ m ∈ Finset.range (n + 1), Q02.Acoef m (n - m) * (m.factorial : ℝ) +
      ∑ m ∈ Finset.range (n + 2), (co m (n + 1 - m)).re * (m.factorial : ℝ) := by
  set K := n + 2 with hK
  obtain ⟨C1, hI1⟩ := I1_expand K
  obtain ⟨C2, hF⟩ := F_expand K
  obtain ⟨C3, h3, h3pos, hE⟩ := Ehat_expand K
  set S : ℕ → ℝ := fun j => ∑ i ∈ Finset.Ico 1 j, Q02.q i with hS
  set a : ℕ × ℕ → ℝ := fun p => Q02.Acoef p.1 p.2 * (p.1.factorial : ℝ) with ha
  set b : ℕ × ℕ → ℝ := fun p => (co p.1 p.2).re * (p.1.factorial : ℝ) with hb
  set P : Polynomial ℝ :=
    (∑ j ∈ Finset.range (K + 1), Polynomial.C (S j) * Polynomial.X ^ (j + 1)) -
    (∑ p ∈ Stri K, Polynomial.C (a p) * Polynomial.X ^ (p.1 + p.2 + 3)) -
    (∑ p ∈ Stri K, Polynomial.C (b p) * Polynomial.X ^ (p.1 + p.2 + 2)) with hP
  have hev : ∀ x : ℝ, P.eval x =
      (∑ j ∈ Finset.range (K + 1), x ^ (j + 1) * S j) -
      x * x * (∑ p ∈ Stri K, Q02.Acoef p.1 p.2 * (p.1.factorial : ℝ) * x ^ (p.1 + p.2 + 1)) -
      x * (∑ p ∈ Stri K, (co p.1 p.2).re * (p.1.factorial : ℝ) * x ^ (p.1 + p.2 + 1)) := by
    intro x
    simp only [hP, Polynomial.eval_sub, Polynomial.eval_finsetSum, Polynomial.eval_mul,
      Polynomial.eval_C, Polynomial.eval_pow, Polynomial.eval_X]
    rw [Finset.mul_sum, Finset.mul_sum]
    congr 1; congr 1
    · refine Finset.sum_congr rfl (fun j _ => ?_); ring
    · refine Finset.sum_congr rfl (fun p _ => ?_); simp only [ha]; ring
    · refine Finset.sum_congr rfl (fun p _ => ?_); simp only [hb]; ring
  have hδ : 0 < min (1 / 4) h3 := lt_min (by norm_num) h3pos
  have hbd : ∀ x : ℝ, 0 < x → x ≤ min (1 / 4) h3 →
      |P.eval x| ≤ (|C1| + |C2| + |C3|) * x ^ (K + 1 + 1) := by
    intro x hx hxδ
    have hx4 : x ≤ 1 / 4 := hxδ.trans (min_le_left _ _)
    have hx3 : x ≤ h3 := hxδ.trans (min_le_right _ _)
    have hsplit := Ehat_split hx (by linarith)
    have e1 := hI1 x hx (by linarith)
    have e2 := hF x hx hx4
    have e3 := hE x hx hx3
    rw [hev]
    set I1 := ∫ t in Ioi (0 : ℝ), Real.exp (-(t / x)) * Q02.Acal t x
    set IF := ∫ t in Ioi (0 : ℝ), Real.exp (-(t / x)) * Fr t x
    set A := ∑ p ∈ Stri K, Q02.Acoef p.1 p.2 * (p.1.factorial : ℝ) * x ^ (p.1 + p.2 + 1)
    set B := ∑ p ∈ Stri K, (co p.1 p.2).re * (p.1.factorial : ℝ) * x ^ (p.1 + p.2 + 1)
    set E := ∑ j ∈ Finset.range (K + 1), x ^ (j + 1) * S j
    have e : E - x * x * A - x * B = -(Q02.Ehat x - E) + x * x * (I1 - A) + x * (IF - B) := by
      rw [hsplit]; ring
    rw [e]
    have hx1 : x ≤ 1 := by linarith
    have hxx : x * x ≤ 1 := by nlinarith
    have b1 : |x * x * (I1 - A)| ≤ |C1| * x ^ (K + 1 + 1) := by
      rw [abs_mul, abs_of_nonneg (by positivity)]
      calc x * x * |I1 - A| ≤ 1 * (C1 * x ^ (K + 2)) :=
            mul_le_mul hxx e1 (abs_nonneg _) zero_le_one
        _ ≤ |C1| * x ^ (K + 1 + 1) := by
            rw [one_mul]; exact mul_le_mul_of_nonneg_right (le_abs_self _) (by positivity)
    have b2 : |x * (IF - B)| ≤ |C2| * x ^ (K + 1 + 1) := by
      rw [abs_mul, abs_of_pos hx]
      calc x * |IF - B| ≤ 1 * (C2 * x ^ (K + 2)) := mul_le_mul hx1 e2 (abs_nonneg _) zero_le_one
        _ ≤ |C2| * x ^ (K + 1 + 1) := by
            rw [one_mul]; exact mul_le_mul_of_nonneg_right (le_abs_self _) (by positivity)
    have b3 : |-(Q02.Ehat x - E)| ≤ |C3| * x ^ (K + 1 + 1) := by
      rw [abs_neg]
      exact e3.trans (mul_le_mul_of_nonneg_right (le_abs_self _) (by positivity))
    calc _ ≤ |-(Q02.Ehat x - E) + x * x * (I1 - A)| + |x * (IF - B)| := abs_add_le _ _
      _ ≤ |-(Q02.Ehat x - E)| + |x * x * (I1 - A)| + |x * (IF - B)| := by
          gcongr; exact abs_add_le _ _
      _ ≤ _ := by linarith
  have hc := poly_coeff_zero (K + 1) P _ _ hδ hbd (K + 1) le_rfl
  rw [hP, Polynomial.coeff_sub, Polynomial.coeff_sub, Polynomial.finsetSum_coeff,
    Polynomial.finsetSum_coeff, Polynomial.finsetSum_coeff] at hc
  simp_rw [Polynomial.coeff_C_mul_X_pow] at hc
  rw [filter_antidiag K n 3 (by omega) (by omega),
    filter_antidiag K (n + 1) 2 (by omega) (by omega)] at hc
  have hS : ∑ j ∈ Finset.range (K + 1), (if K + 1 = j + 1 then S j else 0) = S K := by
    rw [Finset.sum_eq_single K]
    · simp
    · intro j _ hj; rw [if_neg (by omega)]
    · intro h; simp at h
  rw [hS] at hc
  simp only [ha, hb, hS] at hc
  linarith

/-- **Lemma `lem:pip:Tk`**: `|N_k - (-log 2)^{k+1}/(k+1)! - T_k| ≤ C (2^k + k (k-2)! 6^{-k})`. -/
theorem main : Q06.stmt_1 := by
  obtain ⟨Cφ, hφ0, -, hco⟩ := co_bd
  refine ⟨78 * Cφ, fun k hk => ?_⟩
  obtain ⟨n, rfl⟩ : ∃ n, k = n + 2 := ⟨k - 2, by omega⟩
  have hN : Q02.N (n + 2) - (-Real.log 2) ^ (n + 2 + 1) / ((n + 2 + 1).factorial : ℝ) =
      ∑ i ∈ Finset.Ico 1 (n + 2), Q02.q i := by unfold Q02.N; ring
  have hT : Q02.T (n + 2) = ∑ m ∈ Finset.range (n + 1), Q02.Acoef m (n - m) * (m.factorial : ℝ) := by
    unfold Q02.T
    rw [show n + 2 - 1 = n + 1 by omega]
    refine Finset.sum_congr rfl (fun m _ => ?_)
    rw [show n + 2 - 2 - m = n - m by omega]; ring
  rw [hN, coeff_identity, hT, add_sub_cancel_left, show n + 2 - 2 = n by omega]
  have hb : ∀ m ∈ Finset.range (n + 2), |(co m (n + 1 - m)).re * (m.factorial : ℝ)| ≤
      Cφ * ((m.factorial : ℝ) * 2 ^ (n + 1 - m) / 6 ^ m) := by
    intro m _
    rw [abs_mul, Nat.abs_cast]
    calc |(co m (n + 1 - m)).re| * (m.factorial : ℝ) ≤
        (Cφ * 2 ^ (n + 1 - m) / 6 ^ m) * (m.factorial : ℝ) :=
          mul_le_mul_of_nonneg_right ((Complex.abs_re_le_norm _).trans (hco m _)) (by positivity)
      _ = _ := by ring
  calc _ ≤ ∑ m ∈ Finset.range (n + 2), |(co m (n + 1 - m)).re * (m.factorial : ℝ)| :=
        Finset.abs_sum_le_sum_abs _ _
    _ ≤ ∑ m ∈ Finset.range (n + 2), Cφ * ((m.factorial : ℝ) * 2 ^ (n + 1 - m) / 6 ^ m) :=
        Finset.sum_le_sum hb
    _ = Cφ * ∑ m ∈ Finset.range (n + 2), (m.factorial : ℝ) * 2 ^ (n + 1 - m) / 6 ^ m := by
        rw [Finset.mul_sum]
    _ ≤ Cφ * (78 * (2 ^ (n + 2) + ((n + 2 : ℕ) : ℝ) * (n.factorial : ℝ) / 6 ^ (n + 2))) :=
        mul_le_mul_of_nonneg_left (U_arith n) hφ0
    _ = _ := by ring

end Q06.PA
end

section
-- merged from work/Q06/PartB.lean
/-!
# Q06, part B: `lem:pip:TkT1`

`Q06.PB.main : Q06.stmt_2`: `|T_k - T_k^{(1)}| ≤ C (k-2)! π^{-(k-1)} log k / k²` for `k ≥ k₀`.

Route (no Hankel contour; only Cauchy estimates on circles `|z| = π(m+1)/(m+2)`):
* `reduction`: `T_k - T_k^{(1)} = Σ_m Re ∂^m E_{k-2-m}(0) / (k-2-m)!` with
  `E_i = B(-λ)^i - S_+(-L_+)^i - S_-(-L_-)^i`, `B = e^{2z}/(1+e^z)`, `λ = z/2 + Log(2cosh(z/2))`,
  `S_± = (±iπ - z)^{-1}`, `L_± = Log(±iπ - z)`. This uses the Schwarz swap of the `y`/`h`
  derivatives in `Acoef` (`Acoef_eq`), the real/complex coefficients (`bpart_iter_real`), the
  recursion identifying the terms of `T_k^{(1)}` (`gpart_eq_Qf`) and conjugation (`gpart_conj`).
* `bound0` (`E_0 = Q02.g2`, Cauchy on `|z| = 2π`), `bound3` (all `i`, `|E_i| ≤ C(i+1)Λ^i` with
  `Λ = 11 + log(2/t)`), `bound2` (Cauchy for `E_2'`, `|E_2'| ≤ C(Λ² + Λ/t)`), `bound1` (Cauchy for
  `E_1' + (3/2)(S_+ + S_-)`, bounded by `CΛ`; the poles `S_±` have explicit coefficients).
  The local facts at `c = ±iπ`: `Bc_local`, `lcd_local`, `d_local` (via `λ = L_c + u/2 + Log σ(u)`).
* `assemble`: the sum over `i` (terms `i = 0,1,2` and a ratio test for `i ≥ 3`).
-/

open Filter Topology
open scoped ContDiff

namespace Q06.PB

/-! ### Part 1: the summation (pure real analysis) -/

lemma geom_tail (s : ℕ → ℝ) (p : ℕ) (h : ∀ j, j < p → s (j + 1) ≤ s j / 2) (h0 : 0 ≤ s 0) :
    ∑ j ∈ Finset.range (p + 1), s j ≤ 2 * s 0 := by
  have hb : ∀ j, j ≤ p → s j ≤ s 0 * (1 / 2) ^ j := by
    intro j
    induction j with
    | zero => intro _; simp
    | succ j ih =>
      intro hj
      have h1 := h j (by omega)
      have h2 := ih (by omega)
      calc s (j + 1) ≤ s j / 2 := h1
        _ ≤ s 0 * (1 / 2) ^ j / 2 := by linarith
        _ = s 0 * (1 / 2) ^ (j + 1) := by ring
  calc ∑ j ∈ Finset.range (p + 1), s j
      ≤ ∑ j ∈ Finset.range (p + 1), s 0 * (1 / 2 : ℝ) ^ j :=
        Finset.sum_le_sum (fun j hj => hb j (by simp at hj; omega))
    _ = s 0 * ∑ j ∈ Finset.range (p + 1), (1 / 2 : ℝ) ^ j := by rw [Finset.mul_sum]
    _ ≤ s 0 * 2 := mul_le_mul_of_nonneg_left (sum_geometric_two_le _) h0
    _ = 2 * s 0 := by ring

/-- eventually `(log(n+2) + K)^3 ≤ ε n` -/
lemma log_cube_small (K ε : ℝ) (hK : 0 ≤ K) (hε : 0 < ε) :
    ∃ N : ℕ, ∀ n : ℕ, N ≤ n → (Real.log ((n : ℝ) + 2) + K) ^ 3 ≤ ε * n := by
  have ht := Real.tendsto_pow_log_div_mul_add_atTop 1 0 3 one_ne_zero
  have hev := ht.eventually (gt_mem_nhds (show (0 : ℝ) < ε / 16 by positivity))
  obtain ⟨X0, hX0⟩ := Filter.eventually_atTop.mp hev
  refine ⟨Nat.ceil (max X0 (Real.exp K)) + 2, fun n hn => ?_⟩
  have hn' : ((Nat.ceil (max X0 (Real.exp K)) + 2 : ℕ) : ℝ) ≤ n := by exact_mod_cast hn
  push_cast at hn'
  have hc := Nat.le_ceil (max X0 (Real.exp K))
  set x : ℝ := (n : ℝ) + 2 with hx
  have hx0 : X0 ≤ x := by linarith [le_max_left X0 (Real.exp K)]
  have hxK : Real.exp K ≤ x := by linarith [le_max_right X0 (Real.exp K)]
  have hxpos : 0 < x := by linarith [Real.exp_pos K]
  have hlogK : K ≤ Real.log x := by
    rw [← Real.log_exp K]; exact Real.log_le_log (Real.exp_pos K) hxK
  have h1 := hX0 x hx0
  simp only [one_mul, add_zero] at h1
  rw [div_lt_iff₀ hxpos] at h1
  have hn2 : (2 : ℝ) ≤ n := by
    have := Real.exp_pos K
    have : (0 : ℝ) ≤ (Nat.ceil (max X0 (Real.exp K)) : ℝ) := Nat.cast_nonneg _
    linarith
  have hl0 : 0 ≤ Real.log x := le_trans hK hlogK
  have h2 : (Real.log x + K) ^ 3 ≤ (2 * Real.log x) ^ 3 :=
    pow_le_pow_left₀ (by linarith) (by linarith) 3
  calc (Real.log x + K) ^ 3 ≤ (2 * Real.log x) ^ 3 := h2
    _ = 8 * Real.log x ^ 3 := by ring
    _ ≤ 8 * (ε / 16 * x) := by linarith
    _ ≤ ε * n := by rw [hx]; nlinarith

/-- the scale `n!/π^{n+1} · log(n+2)/(n+2)^2` (with `n = k - 2`) -/
noncomputable def X (n : ℕ) : ℝ :=
  (n.factorial : ℝ) / Real.pi ^ (n + 1) * (Real.log ((n : ℝ) + 2) / ((n : ℝ) + 2) ^ 2)

lemma X_nonneg (n : ℕ) : 0 ≤ X n := by
  unfold X
  have : 0 ≤ Real.log ((n : ℝ) + 2) := Real.log_nonneg (by linarith [(Nat.cast_nonneg n : (0:ℝ) ≤ n)])
  have := Real.pi_pos
  positivity

lemma sq_le_two_pow (n : ℕ) (hn : 6 ≤ n) : ((n : ℝ) + 2) ^ 2 ≤ 2 ^ n := by
  induction n, hn using Nat.le_induction with
  | base => norm_num
  | succ n hn ih =>
    have h6 : (6 : ℝ) ≤ n := by exact_mod_cast hn
    have : ((n : ℝ) + 1 + 2) ^ 2 ≤ 2 * ((n : ℝ) + 2) ^ 2 := by nlinarith
    have e : (2 : ℝ) ^ (n + 1) = 2 * 2 ^ n := by ring
    push_cast; rw [e]; linarith

lemma one_le_log (n : ℕ) (hn : 1 ≤ n) : 1 ≤ Real.log ((n : ℝ) + 2) := by
  rw [← Real.log_exp 1]
  apply Real.log_le_log (Real.exp_pos 1)
  have := Real.exp_one_lt_d9
  have : (1 : ℝ) ≤ n := by exact_mod_cast hn
  linarith

lemma term0_le (n : ℕ) (hn : 6 ≤ n) :
    (n.factorial : ℝ) / (2 * Real.pi) ^ n ≤ Real.pi * X n := by
  have hpi := Real.pi_pos
  have hl := one_le_log n (by omega)
  have hs := sq_le_two_pow n hn
  have h1 : 1 / (2 : ℝ) ^ n ≤ Real.log ((n : ℝ) + 2) / ((n : ℝ) + 2) ^ 2 := by
    rw [div_le_div_iff₀ (by positivity) (by positivity)]; nlinarith
  calc (n.factorial : ℝ) / (2 * Real.pi) ^ n = (n.factorial : ℝ) / Real.pi ^ n * (1 / 2 ^ n) := by
        rw [mul_pow]; field_simp
    _ ≤ (n.factorial : ℝ) / Real.pi ^ n * (Real.log ((n : ℝ) + 2) / ((n : ℝ) + 2) ^ 2) :=
        mul_le_mul_of_nonneg_left h1 (by positivity)
    _ = Real.pi * X n := by unfold X; field_simp; ring

lemma term12_le (p : ℕ) (hp : 1 ≤ p) {L : ℝ} (hL : L ≤ Real.log ((p : ℝ) + 5)) :
    ((p + 1).factorial : ℝ) * L / Real.pi ^ (p + 2) ≤ 4 * Real.pi ^ 2 * X (p + 3) := by
  have hpi := Real.pi_pos
  have hp1 : (1 : ℝ) ≤ p := by exact_mod_cast hp
  have hf : ((p + 3).factorial : ℝ) = ((p : ℝ) + 3) * ((p : ℝ) + 2) * ((p + 1).factorial : ℝ) := by
    rw [show p + 3 = (p + 1) + 1 + 1 by ring, Nat.factorial_succ, Nat.factorial_succ]; push_cast; ring
  have e : 4 * Real.pi ^ 2 * X (p + 3) =
      4 * ((p + 3).factorial : ℝ) * Real.log ((p : ℝ) + 5) / ((p : ℝ) + 5) ^ 2 / Real.pi ^ (p + 2) := by
    unfold X; push_cast
    rw [show (p : ℝ) + 3 + 2 = (p : ℝ) + 5 by ring, show p + 3 + 1 = (p + 2) + 2 by ring, pow_add]
    field_simp
  rw [e, hf]
  have hF : (0 : ℝ) < ((p + 1).factorial : ℝ) := by exact_mod_cast Nat.factorial_pos _
  apply div_le_div_of_nonneg_right _ (by positivity)
  rw [le_div_iff₀ (by positivity)]
  have hq : ((p : ℝ) + 5) ^ 2 ≤ 4 * (((p : ℝ) + 3) * ((p : ℝ) + 2)) := by nlinarith
  have hl0 : 0 ≤ Real.log ((p : ℝ) + 5) := Real.log_nonneg (by linarith)
  calc ((p + 1).factorial : ℝ) * L * ((p : ℝ) + 5) ^ 2
      ≤ ((p + 1).factorial : ℝ) * Real.log ((p : ℝ) + 5) * (4 * (((p : ℝ) + 3) * ((p : ℝ) + 2))) := by
        apply mul_le_mul (mul_le_mul_of_nonneg_left hL hF.le) hq (by positivity) (by positivity)
    _ = _ := by ring

lemma tail_ratio (p j : ℕ) (hj : j < p) (x : ℝ) (hx : 0 ≤ x) (h2x : 2 * x ≤ (p : ℝ) + 3) :
    ((p - (j + 1)).factorial : ℝ) * x ^ (j + 1 + 3) / ((j + 1 + 3).factorial : ℝ) ≤
      ((p - j).factorial : ℝ) * x ^ (j + 3) / ((j + 3).factorial : ℝ) / 2 := by
  obtain ⟨q, rfl⟩ : ∃ q, p = j + 1 + q := ⟨p - j - 1, by omega⟩
  rw [show j + 1 + q - (j + 1) = q by omega, show j + 1 + q - j = q + 1 by omega,
    show j + 1 + 3 = (j + 3) + 1 by ring, Nat.factorial_succ (j + 3), Nat.factorial_succ q]
  set A := (q.factorial : ℝ) * x ^ (j + 3) / ((j + 3).factorial : ℝ) with hA
  have hA0 : 0 ≤ A := by positivity
  have hF : (0 : ℝ) < ((j + 3).factorial : ℝ) := by exact_mod_cast Nat.factorial_pos _
  have e1 : (q.factorial : ℝ) * x ^ (j + 3 + 1) / ((((j + 3 + 1 : ℕ) : ℝ)) * ((j + 3).factorial : ℝ)) =
      A * (x / ((j : ℝ) + 4)) := by
    rw [hA]; push_cast; field_simp; ring
  have e2 : (((q + 1 : ℕ) : ℝ) * (q.factorial : ℝ)) * x ^ (j + 3) / ((j + 3).factorial : ℝ) / 2 =
      A * (((q : ℝ) + 1) / 2) := by
    rw [hA]; push_cast; field_simp
  push_cast at e1 e2 ⊢
  rw [e1, e2]
  apply mul_le_mul_of_nonneg_left _ hA0
  rw [div_le_div_iff₀ (by positivity) (by norm_num)]
  push_cast at h2x
  have hq0 : (0 : ℝ) ≤ q := Nat.cast_nonneg q
  have hj0 : (0 : ℝ) ≤ j := Nat.cast_nonneg j
  nlinarith

lemma tail_sum (p : ℕ) (x : ℝ) (hx : 0 ≤ x) (h2x : 2 * x ≤ (p : ℝ) + 3) :
    ∑ j ∈ Finset.range (p + 1), ((p - j).factorial : ℝ) * x ^ (j + 3) / ((j + 3).factorial : ℝ) ≤
      (p.factorial : ℝ) * x ^ 3 / 3 := by
  have := geom_tail (fun j => ((p - j).factorial : ℝ) * x ^ (j + 3) / ((j + 3).factorial : ℝ)) p
    (fun j hj => tail_ratio p j hj x hx h2x) (by positivity)
  refine this.trans (le_of_eq ?_)
  simp only [Nat.sub_zero, zero_add]
  rw [show (3 : ℕ).factorial = 6 by rfl]
  push_cast; ring

lemma tail_final (p : ℕ) (hp : 1 ≤ p) {x : ℝ} (hx0 : 0 ≤ x)
    (hx3 : x ^ 3 ≤ 2 * Real.pi ^ 2 * ((p : ℝ) + 3)) :
    (p.factorial : ℝ) * x ^ 3 / 3 / Real.pi ^ (p + 3) ≤ 6 * Real.pi ^ 3 * X (p + 3) := by
  have hpi := Real.pi_pos
  have hp1 : (1 : ℝ) ≤ p := by exact_mod_cast hp
  have hf : ((p + 3).factorial : ℝ) =
      ((p : ℝ) + 3) * ((p : ℝ) + 2) * ((p : ℝ) + 1) * (p.factorial : ℝ) := by
    rw [show p + 3 = p + 1 + 1 + 1 by ring, Nat.factorial_succ, Nat.factorial_succ,
      Nat.factorial_succ]; push_cast; ring
  have hl : 1 ≤ Real.log ((p : ℝ) + 5) := by
    have := one_le_log (p + 3) (by omega); push_cast at this
    rw [show (p : ℝ) + 3 + 2 = p + 5 by ring] at this; exact this
  have hF : (0 : ℝ) < (p.factorial : ℝ) := by exact_mod_cast Nat.factorial_pos _
  have e : 6 * Real.pi ^ 3 * X (p + 3) =
      6 * ((p + 3).factorial : ℝ) * Real.log ((p : ℝ) + 5) / ((p : ℝ) + 5) ^ 2 / Real.pi ^ (p + 1) := by
    unfold X; push_cast
    rw [show (p : ℝ) + 3 + 2 = (p : ℝ) + 5 by ring, show p + 3 + 1 = (p + 1) + 3 by ring, pow_add]
    field_simp
  rw [e, hf]
  have e2 : (p.factorial : ℝ) * x ^ 3 / 3 / Real.pi ^ (p + 3) =
      (p.factorial : ℝ) * (x ^ 3 / Real.pi ^ 2) / 3 / Real.pi ^ (p + 1) := by
    rw [show p + 3 = (p + 1) + 2 by ring, pow_add]; field_simp
  rw [e2]
  apply div_le_div_of_nonneg_right _ (by positivity)
  have hx3' : x ^ 3 / Real.pi ^ 2 ≤ 2 * ((p : ℝ) + 3) := by
    rw [div_le_iff₀ (by positivity)]; linarith
  rw [div_le_div_iff₀ (by norm_num) (by positivity)]
  have hq : ((p : ℝ) + 5) ^ 2 ≤ 9 * (((p : ℝ) + 2) * ((p : ℝ) + 1)) := by nlinarith
  have h1 : (p.factorial : ℝ) * (x ^ 3 / Real.pi ^ 2) ≤ (p.factorial : ℝ) * (2 * ((p : ℝ) + 3)) :=
    mul_le_mul_of_nonneg_left hx3' hF.le
  have h2 : 0 ≤ (p.factorial : ℝ) * (x ^ 3 / Real.pi ^ 2) := by positivity
  calc (p.factorial : ℝ) * (x ^ 3 / Real.pi ^ 2) * ((p : ℝ) + 5) ^ 2
      ≤ (p.factorial : ℝ) * (2 * ((p : ℝ) + 3)) * (9 * (((p : ℝ) + 2) * ((p : ℝ) + 1))) :=
        mul_le_mul h1 hq (by positivity) (by positivity)
    _ ≤ 6 * (((p : ℝ) + 3) * ((p : ℝ) + 2) * ((p : ℝ) + 1) * (p.factorial : ℝ)) *
        Real.log ((p : ℝ) + 5) * 3 := by
        have : 0 ≤ (p.factorial : ℝ) * (((p : ℝ) + 3) * ((p : ℝ) + 2) * ((p : ℝ) + 1)) := by positivity
        nlinarith

lemma nonneg_of_le_mul_div {a C A D : ℝ} (ha : 0 ≤ a) (hA : 0 < A) (hD : 0 < D)
    (h : a ≤ C * A / D) : 0 ≤ C := by
  by_contra hc
  push Not at hc
  have : C * A / D < 0 := div_neg_of_neg_of_pos (mul_neg_of_neg_of_pos hc hA) hD
  linarith

theorem assemble (a : ℕ → ℕ → ℝ) (C0 C1 C2 C3 K : ℝ) (hK : 0 ≤ K)
    (ha : ∀ i m, 0 ≤ a i m)
    (h0 : ∀ m, a 0 m ≤ C0 * (m.factorial : ℝ) / (2 * Real.pi) ^ m)
    (h1 : ∀ m, 1 ≤ m → a 1 m ≤ C1 * (((m - 1).factorial : ℝ) * Real.log ((m : ℝ) + 2)) / Real.pi ^ m)
    (h2 : ∀ m, 1 ≤ m → a 2 m ≤ C2 * ((m.factorial : ℝ) * Real.log ((m : ℝ) + 2)) / Real.pi ^ m)
    (h3 : ∀ i m, a i m ≤ C3 * (((i : ℝ) + 1) * (m.factorial : ℝ) *
      (Real.log ((m : ℝ) + 2) + K) ^ i) / Real.pi ^ m) :
    ∃ C : ℝ, ∃ N : ℕ, ∀ n, N ≤ n →
      ∑ i ∈ Finset.range (n + 1), a i (n - i) / (i.factorial : ℝ) ≤ C * X n := by
  have hpi := Real.pi_pos
  have hC0 : 0 ≤ C0 := nonneg_of_le_mul_div (ha 0 0) (by positivity) (by positivity) (h0 0)
  have hC1 : 0 ≤ C1 := nonneg_of_le_mul_div (ha 1 1)
    (by have := Real.log_pos (show (1 : ℝ) < (1 : ℕ) + 2 by norm_num); positivity)
    (by positivity) (h1 1 le_rfl)
  have hC2 : 0 ≤ C2 := nonneg_of_le_mul_div (ha 2 1)
    (by have := Real.log_pos (show (1 : ℝ) < (1 : ℕ) + 2 by norm_num); positivity)
    (by positivity) (h2 1 le_rfl)
  have hC3 : 0 ≤ C3 := nonneg_of_le_mul_div (ha 0 0)
    (by positivity) (by positivity) (h3 0 0)
  obtain ⟨N1, hN1⟩ := log_cube_small K (1 / (4 * Real.pi)) hK (by positivity)
  refine ⟨Real.pi * C0 + 4 * Real.pi ^ 2 * C1 + 2 * Real.pi ^ 3 * C2 + 6 * Real.pi ^ 3 * C3,
    max N1 8, fun n hn => ?_⟩
  obtain ⟨p, rfl⟩ : ∃ p, n = p + 3 := ⟨n - 3, by omega⟩
  have hp : 5 ≤ p := by omega
  have hXn := X_nonneg (p + 3)
  rw [show p + 3 + 1 = p + 1 + 1 + 1 + 1 by ring, Finset.sum_range_succ', Finset.sum_range_succ',
    Finset.sum_range_succ']
  simp only [Nat.sub_zero, Nat.factorial_zero, Nat.cast_one, div_one, zero_add]
  -- the term `i = 0`
  have b0 : a 0 (p + 3) ≤ Real.pi * C0 * X (p + 3) := by
    calc a 0 (p + 3) ≤ C0 * ((p + 3).factorial : ℝ) / (2 * Real.pi) ^ (p + 3) := h0 _
      _ = C0 * (((p + 3).factorial : ℝ) / (2 * Real.pi) ^ (p + 3)) := by ring
      _ ≤ C0 * (Real.pi * X (p + 3)) := mul_le_mul_of_nonneg_left (term0_le _ (by omega)) hC0
      _ = _ := by ring
  -- the term `i = 1`
  have b1 : a 1 (p + 3 - 1) / ((1 : ℕ).factorial : ℝ) ≤ 4 * Real.pi ^ 2 * C1 * X (p + 3) := by
    rw [show p + 3 - 1 = p + 2 by omega, Nat.factorial_one, Nat.cast_one, div_one]
    have := h1 (p + 2) (by omega)
    rw [show p + 2 - 1 = p + 1 by omega] at this
    have hl : Real.log (((p + 2 : ℕ) : ℝ) + 2) ≤ Real.log ((p : ℝ) + 5) :=
      Real.log_le_log (by positivity) (by push_cast; linarith)
    have hl0 : 0 ≤ Real.log (((p + 2 : ℕ) : ℝ) + 2) := Real.log_nonneg (by push_cast; linarith)
    have t := term12_le p (by omega) hl
    calc a 1 (p + 2) ≤ C1 * (((p + 1).factorial : ℝ) * Real.log (((p + 2 : ℕ) : ℝ) + 2)) /
          Real.pi ^ (p + 2) := this
      _ = C1 * (((p + 1).factorial : ℝ) * Real.log (((p + 2 : ℕ) : ℝ) + 2) / Real.pi ^ (p + 2)) := by
          ring
      _ ≤ C1 * (4 * Real.pi ^ 2 * X (p + 3)) := mul_le_mul_of_nonneg_left t hC1
      _ = _ := by ring
  -- the term `i = 2`
  have b2 : a 2 (p + 3 - (1 + 1)) / ((1 + 1 : ℕ).factorial : ℝ) ≤
      2 * Real.pi ^ 3 * C2 * X (p + 3) := by
    rw [show p + 3 - (1 + 1) = p + 1 by omega, show (1 + 1 : ℕ).factorial = 2 by rfl]
    have := h2 (p + 1) (by omega)
    have hl : Real.log (((p + 1 : ℕ) : ℝ) + 2) ≤ Real.log ((p : ℝ) + 5) :=
      Real.log_le_log (by positivity) (by push_cast; linarith)
    have t := term12_le p (by omega) hl
    have e : C2 * (((p + 1).factorial : ℝ) * Real.log (((p + 1 : ℕ) : ℝ) + 2)) / Real.pi ^ (p + 1) =
        C2 * Real.pi * (((p + 1).factorial : ℝ) * Real.log (((p + 1 : ℕ) : ℝ) + 2) /
          Real.pi ^ (p + 2)) := by
      field_simp; ring
    rw [e] at this
    have h4 : C2 * Real.pi * (((p + 1).factorial : ℝ) * Real.log (((p + 1 : ℕ) : ℝ) + 2) /
          Real.pi ^ (p + 2)) ≤ C2 * Real.pi * (4 * Real.pi ^ 2 * X (p + 3)) :=
      mul_le_mul_of_nonneg_left t (by positivity)
    push_cast
    linarith
  -- the terms `i ≥ 3`
  set Y := Real.log ((p : ℝ) + 5) + K with hY
  have hY1 : 1 ≤ Y := by
    have := one_le_log (p + 3) (by omega); push_cast at this
    rw [show (p : ℝ) + 3 + 2 = p + 5 by ring] at this; rw [hY]; linarith
  have hY3 : Y ^ 3 ≤ 1 / (4 * Real.pi) * ((p : ℝ) + 3) := by
    have := hN1 (p + 3) (le_trans (le_max_left _ _) hn)
    push_cast at this
    rw [show (p : ℝ) + 3 + 2 = p + 5 by ring] at this
    exact this
  set x := 2 * Real.pi * Y with hx
  have hx0 : 0 ≤ x := by positivity
  have hYY : Y ≤ Y ^ 3 := by
    have h1 : 1 ≤ Y ^ 2 := by nlinarith
    nlinarith
  have h2x : 2 * x ≤ (p : ℝ) + 3 := by
    have : 4 * Real.pi * Y ≤ 4 * Real.pi * (1 / (4 * Real.pi) * ((p : ℝ) + 3)) :=
      mul_le_mul_of_nonneg_left (hYY.trans hY3) (by positivity)
    rw [hx]; field_simp at this; linarith
  have hx3 : x ^ 3 ≤ 2 * Real.pi ^ 2 * ((p : ℝ) + 3) := by
    have : 8 * Real.pi ^ 3 * Y ^ 3 ≤ 8 * Real.pi ^ 3 * (1 / (4 * Real.pi) * ((p : ℝ) + 3)) :=
      mul_le_mul_of_nonneg_left hY3 (by positivity)
    rw [hx, mul_pow, mul_pow]; field_simp at this; nlinarith
  have bt : ∑ j ∈ Finset.range (p + 1), a (j + 1 + 1 + 1) (p + 3 - (j + 1 + 1 + 1)) /
      ((j + 1 + 1 + 1).factorial : ℝ) ≤ 6 * Real.pi ^ 3 * C3 * X (p + 3) := by
    have hterm : ∀ j ∈ Finset.range (p + 1), a (j + 1 + 1 + 1) (p + 3 - (j + 1 + 1 + 1)) /
        ((j + 1 + 1 + 1).factorial : ℝ) ≤ C3 / Real.pi ^ (p + 3) *
          (((p - j).factorial : ℝ) * x ^ (j + 3) / ((j + 3).factorial : ℝ)) := by
      intro j hj
      have hjp : j ≤ p := by simp at hj; omega
      rw [show j + 1 + 1 + 1 = j + 3 by ring, show p + 3 - (j + 3) = p - j by omega]
      have hF : (0 : ℝ) < ((j + 3).factorial : ℝ) := by exact_mod_cast Nat.factorial_pos _
      have e : C3 / Real.pi ^ (p + 3) *
          (((p - j).factorial : ℝ) * x ^ (j + 3) / ((j + 3).factorial : ℝ)) =
          C3 * (2 ^ (j + 3) * ((p - j).factorial : ℝ) * Y ^ (j + 3)) / Real.pi ^ (p - j) /
            ((j + 3).factorial : ℝ) := by
        rw [hx, mul_pow, mul_pow, show p + 3 = (p - j) + (j + 3) by omega, pow_add]
        field_simp
      rw [e]
      apply div_le_div_of_nonneg_right _ hF.le
      refine (h3 (j + 3) (p - j)).trans ?_
      apply div_le_div_of_nonneg_right _ (by positivity)
      apply mul_le_mul_of_nonneg_left _ hC3
      have h4 : ((j + 3 : ℕ) : ℝ) + 1 ≤ 2 ^ (j + 3) := by
        have : j + 3 + 1 ≤ 2 ^ (j + 3) := Nat.lt_two_pow_self
        exact_mod_cast this
      have hL0 : 0 ≤ Real.log (((p - j : ℕ) : ℝ) + 2) + K :=
        add_nonneg (Real.log_nonneg (by linarith [(Nat.cast_nonneg (p - j) : (0 : ℝ) ≤ _)])) hK
      have hLY : Real.log (((p - j : ℕ) : ℝ) + 2) + K ≤ Y := by
        rw [hY]
        have : ((p - j : ℕ) : ℝ) ≤ p := by exact_mod_cast Nat.sub_le p j
        have := Real.log_le_log (by positivity) (show ((p - j : ℕ) : ℝ) + 2 ≤ (p : ℝ) + 5 by linarith)
        linarith
      have hpow := pow_le_pow_left₀ hL0 hLY (j + 3)
      apply mul_le_mul (mul_le_mul_of_nonneg_right h4 (by positivity)) hpow (by positivity)
        (by positivity)
    calc _ ≤ ∑ j ∈ Finset.range (p + 1), C3 / Real.pi ^ (p + 3) *
          (((p - j).factorial : ℝ) * x ^ (j + 3) / ((j + 3).factorial : ℝ)) :=
          Finset.sum_le_sum hterm
      _ = C3 / Real.pi ^ (p + 3) * ∑ j ∈ Finset.range (p + 1),
          (((p - j).factorial : ℝ) * x ^ (j + 3) / ((j + 3).factorial : ℝ)) := by
          rw [Finset.mul_sum]
      _ ≤ C3 / Real.pi ^ (p + 3) * ((p.factorial : ℝ) * x ^ 3 / 3) :=
          mul_le_mul_of_nonneg_left (tail_sum p x hx0 h2x) (by positivity)
      _ = C3 * ((p.factorial : ℝ) * x ^ 3 / 3 / Real.pi ^ (p + 3)) := by ring
      _ ≤ C3 * (6 * Real.pi ^ 3 * X (p + 3)) :=
          mul_le_mul_of_nonneg_left (tail_final p (by omega) hx0 hx3) hC3
      _ = _ := by ring
  linarith

/-! ### Part 2: the terms of `T_k^{(1)}` as Taylor coefficients of `(c-z)^{-1} (-Log(c-z))^i` -/

/-- `y_0 = iπ` -/
noncomputable def ip : ℂ := Complex.I * Real.pi

/-- `Q_m(h) = c^{-(m+1)} e^{-h Log c} (1+h)_m` -/
noncomputable def Qf (c : ℂ) (m : ℕ) (h : ℂ) : ℂ :=
  (c ^ (m + 1))⁻¹ * Complex.exp (-(h * Complex.log c)) * (ascPochhammer ℂ m).eval (1 + h)

/-- `g_i(z) = (c-z)^{-1} (-Log(c-z))^i` -/
noncomputable def gpart (c : ℂ) (i : ℕ) (z : ℂ) : ℂ := (c - z)⁻¹ * (-Complex.log (c - z)) ^ i

lemma Qf_contDiff (c : ℂ) (m : ℕ) : ContDiff ℂ ⊤ (Qf c m) := by
  have : Differentiable ℂ (Qf c m) := by
    unfold Qf
    have hp0 : Differentiable ℂ (fun x : ℂ => (ascPochhammer ℂ m).eval x) :=
      Polynomial.differentiable _
    have hp : Differentiable ℂ (fun h : ℂ => (ascPochhammer ℂ m).eval (1 + h)) :=
      hp0.comp (differentiable_const _ |>.add differentiable_id)
    fun_prop
  exact this.contDiff

lemma Qf_succ (c : ℂ) (m : ℕ) (h : ℂ) :
    Qf c (m + 1) h = c⁻¹ * (((m : ℂ) + 1) * Qf c m h + h * Qf c m h) := by
  unfold Qf
  rw [ascPochhammer_succ_eval, pow_succ, mul_inv]
  ring

lemma iter_mul_id (f : ℂ → ℂ) (n : ℕ) (hf : ContDiffAt ℂ n f 0) :
    iteratedDeriv n (fun h => h * f h) 0 = n * iteratedDeriv (n - 1) f 0 := by
  rw [iteratedDeriv_fun_mul (f := fun h : ℂ => h) (g := f) contDiffAt_id hf]
  simp only [iteratedDeriv_fun_id_zero]
  rcases n with _ | n
  · simp
  · rw [Finset.sum_eq_single 1]
    · simp
    · intro b _ hb; simp [hb]
    · intro h; simp at h

lemma Qf_iter_succ (c : ℂ) (m i : ℕ) :
    iteratedDeriv i (Qf c (m + 1)) 0 = c⁻¹ * (((m : ℂ) + 1) * iteratedDeriv i (Qf c m) 0 +
      (i : ℂ) * iteratedDeriv (i - 1) (Qf c m) 0) := by
  have e : Qf c (m + 1) = fun h => c⁻¹ * (((m : ℂ) + 1) * Qf c m h + h * Qf c m h) := by
    funext h; exact Qf_succ c m h
  have hQ := Qf_contDiff c m
  rw [e, iteratedDeriv_const_mul_field]
  have h1 : ContDiffAt ℂ i (fun h => ((m : ℂ) + 1) * Qf c m h) 0 :=
    (contDiffAt_const.mul (hQ.contDiffAt)).of_le le_top
  have h2 : ContDiffAt ℂ i (fun h => h * Qf c m h) 0 :=
    (contDiffAt_id.mul (hQ.contDiffAt)).of_le le_top
  rw [iteratedDeriv_fun_add h1 h2, iteratedDeriv_const_mul_field,
    iter_mul_id _ _ (hQ.contDiffAt.of_le le_top)]

lemma Qf_iter_zero (c : ℂ) (i : ℕ) :
    iteratedDeriv i (Qf c 0) 0 = c⁻¹ * (-Complex.log c) ^ i := by
  have e : Qf c 0 = fun h => c⁻¹ * Complex.exp ((-Complex.log c) * h) := by
    funext h
    simp only [Qf, zero_add, pow_one, ascPochhammer_zero, Polynomial.eval_one, mul_one]
    ring_nf
  rw [e, iteratedDeriv_const_mul_field, iteratedDeriv_cexp_const_mul]
  simp

lemma isOpen_slitU (c : ℂ) : IsOpen {z : ℂ | c - z ∈ Complex.slitPlane} :=
  Complex.isOpen_slitPlane.preimage (continuous_const.sub continuous_id)

lemma gpart_analyticAt (c z : ℂ) (hz : c - z ∈ Complex.slitPlane) (i : ℕ) :
    AnalyticAt ℂ (gpart c i) z := by
  have hne : c - z ≠ 0 := Complex.slitPlane_ne_zero hz
  have h1 : AnalyticAt ℂ (fun z : ℂ => c - z) z := analyticAt_const.sub analyticAt_id
  unfold gpart
  exact (h1.inv hne).mul ((h1.clog hz).neg.pow i)

lemma gpart_hasDerivAt (c z : ℂ) (hz : c - z ∈ Complex.slitPlane) (i : ℕ) :
    HasDerivAt (gpart c i) ((c - z)⁻¹ * (gpart c i z + (i : ℂ) * gpart c (i - 1) z)) z := by
  have hne : c - z ≠ 0 := Complex.slitPlane_ne_zero hz
  have h1 : HasDerivAt (fun z : ℂ => c - z) (-1) z := by
    simpa using (hasDerivAt_id z).const_sub c
  have h2 : HasDerivAt (fun z => (c - z)⁻¹) (-(-1) / (c - z) ^ 2) z := h1.inv hne
  have h3 : HasDerivAt (fun z => (-Complex.log (c - z)) ^ i)
      ((i : ℂ) * (-Complex.log (c - z)) ^ (i - 1) * (-((-1) / (c - z)))) z :=
    ((h1.clog hz).neg).pow i
  have h4 : HasDerivAt (gpart c i) (-(-1) / (c - z) ^ 2 * (-Complex.log (c - z)) ^ i +
      (c - z)⁻¹ * ((i : ℂ) * (-Complex.log (c - z)) ^ (i - 1) * (-((-1) / (c - z))))) z :=
    h2.mul h3
  convert h4 using 1
  rcases i with _ | i
  · simp [gpart]; field_simp
  · simp only [gpart, Nat.add_sub_cancel]
    field_simp

lemma gpart_iter_succ (c : ℂ) (hc : c ∈ Complex.slitPlane) (m i : ℕ) :
    c * iteratedDeriv (m + 1) (gpart c i) 0 =
      ((m : ℂ) + 1) * iteratedDeriv m (gpart c i) 0 +
        (i : ℂ) * iteratedDeriv m (gpart c (i - 1)) 0 := by
  have hU := isOpen_slitU c
  have h0 : (0 : ℂ) ∈ {z : ℂ | c - z ∈ Complex.slitPlane} := by simpa using hc
  have han : ∀ j, AnalyticOnNhd ℂ (gpart c j) {z : ℂ | c - z ∈ Complex.slitPlane} :=
    fun j z hz => gpart_analyticAt c z hz j
  have hd : AnalyticAt ℂ (deriv (gpart c i)) 0 := ((han i).deriv) 0 h0
  have heq : (fun z => c * deriv (gpart c i) z - z * deriv (gpart c i) z) =ᶠ[𝓝 0]
      (fun z => gpart c i z + (i : ℂ) * gpart c (i - 1) z) := by
    filter_upwards [hU.mem_nhds h0] with z hz
    rw [(gpart_hasDerivAt c z hz i).deriv]
    have hne : c - z ≠ 0 := Complex.slitPlane_ne_zero hz
    field_simp
  have H := heq.iteratedDeriv_eq m
  have hA : ContDiffAt ℂ m (fun z => c * deriv (gpart c i) z) 0 :=
    (contDiffAt_const.mul hd.contDiffAt).of_le le_top
  have hB : ContDiffAt ℂ m (fun z => z * deriv (gpart c i) z) 0 :=
    (contDiffAt_id.mul hd.contDiffAt).of_le le_top
  have hC : ContDiffAt ℂ m (gpart c i) 0 := ((han i) 0 h0).contDiffAt.of_le le_top
  have hD : ContDiffAt ℂ m (fun z => (i : ℂ) * gpart c (i - 1) z) 0 :=
    (contDiffAt_const.mul ((han (i - 1)) 0 h0).contDiffAt).of_le le_top
  rw [iteratedDeriv_fun_sub hA hB, iteratedDeriv_fun_add hC hD, iteratedDeriv_const_mul_field,
    iteratedDeriv_const_mul_field, iter_mul_id _ _ (hd.contDiffAt.of_le le_top)] at H
  have e1 : ∀ n, iteratedDeriv n (deriv (gpart c i)) = iteratedDeriv (n + 1) (gpart c i) :=
    fun n => iteratedDeriv_succ'.symm
  rw [e1, e1] at H
  rcases m with _ | m
  · simp at H ⊢; linear_combination H
  · rw [show m + 1 - 1 = m by omega] at H
    push_cast at H ⊢
    linear_combination H

/-- the terms of `T_k^{(1)}`: `[h^i] Q_m = m! [z^m] g_i` -/
theorem gpart_eq_Qf (c : ℂ) (hc : c ∈ Complex.slitPlane) (m : ℕ) :
    ∀ i, iteratedDeriv m (gpart c i) 0 = iteratedDeriv i (Qf c m) 0 := by
  have hc0 : c ≠ 0 := Complex.slitPlane_ne_zero hc
  induction m with
  | zero =>
    intro i
    rw [Qf_iter_zero, iteratedDeriv_zero]
    simp [gpart]
  | succ m ih =>
    intro i
    rw [Qf_iter_succ, ← ih i, ← ih (i - 1)]
    have H := gpart_iter_succ c hc m i
    field_simp
    linear_combination H

lemma Qf_conj (c : ℂ) (hl : Complex.log ((starRingEnd ℂ) c) = (starRingEnd ℂ) (Complex.log c))
    (m : ℕ) : ∀ i, iteratedDeriv i (Qf ((starRingEnd ℂ) c) m) 0 =
      (starRingEnd ℂ) (iteratedDeriv i (Qf c m) 0) := by
  induction m with
  | zero => intro i; rw [Qf_iter_zero, Qf_iter_zero, hl]; simp
  | succ m ih => intro i; rw [Qf_iter_succ, Qf_iter_succ, ih i, ih (i - 1)]; simp

lemma ip_slit : ip ∈ Complex.slitPlane := by
  rw [Complex.mem_slitPlane_iff]; right; simp [ip, Real.pi_ne_zero]

lemma nip_slit : -ip ∈ Complex.slitPlane := by
  rw [Complex.mem_slitPlane_iff]; right; simp [ip, Real.pi_ne_zero]

lemma conj_ip : (starRingEnd ℂ) ip = -ip := by simp [ip]

lemma log_conj_ip : Complex.log ((starRingEnd ℂ) ip) = (starRingEnd ℂ) (Complex.log ip) := by
  apply Complex.log_conj
  rw [ip, show Complex.I * (Real.pi : ℂ) = (Real.pi : ℂ) * Complex.I by ring,
    Complex.arg_real_mul _ Real.pi_pos, Complex.arg_I]
  linarith [Real.pi_pos]

/-- the `-iπ` terms are the conjugates of the `iπ` terms -/
lemma gpart_conj (m i : ℕ) :
    iteratedDeriv m (gpart (-ip) i) 0 = (starRingEnd ℂ) (iteratedDeriv m (gpart ip i) 0) := by
  rw [gpart_eq_Qf _ nip_slit, gpart_eq_Qf _ ip_slit, ← conj_ip, Qf_conj _ log_conj_ip]

/-! ### Part 3: swapping the `y`- and `h`-derivatives in `[y^m h^i] 𝓐` -/

/-- `∂_y F` -/
noncomputable def pY (F : ℝ × ℝ → ℝ) (q : ℝ × ℝ) : ℝ := deriv (fun y => F (y, q.2)) q.1
/-- `∂_h F` -/
noncomputable def pH (F : ℝ × ℝ → ℝ) (q : ℝ × ℝ) : ℝ := deriv (fun h => F (q.1, h)) q.2

lemma hasFDerivAt_of_smooth {F : ℝ × ℝ → ℝ} (hF : ContDiff ℝ ∞ F) (q : ℝ × ℝ) :
    HasFDerivAt F (fderiv ℝ F q) q :=
  ((hF.differentiable (by simp)) q).hasFDerivAt

lemma pY_eq {F : ℝ × ℝ → ℝ} (hF : ContDiff ℝ ∞ F) : pY F = fun q => fderiv ℝ F q (1, 0) := by
  funext q
  have h1 : HasDerivAt (fun y : ℝ => (y, q.2)) ((1 : ℝ), (0 : ℝ)) q.1 :=
    (hasDerivAt_id q.1).prodMk (hasDerivAt_const q.1 q.2)
  have h2 := (hasFDerivAt_of_smooth hF (q.1, q.2)).comp_hasDerivAt q.1 h1
  exact h2.deriv

lemma pH_eq {F : ℝ × ℝ → ℝ} (hF : ContDiff ℝ ∞ F) : pH F = fun q => fderiv ℝ F q (0, 1) := by
  funext q
  have h1 : HasDerivAt (fun h : ℝ => (q.1, h)) ((0 : ℝ), (1 : ℝ)) q.2 :=
    (hasDerivAt_const q.2 q.1).prodMk (hasDerivAt_id q.2)
  have h2 := (hasFDerivAt_of_smooth hF (q.1, q.2)).comp_hasDerivAt q.2 h1
  exact h2.deriv

lemma apply_smooth {F : ℝ × ℝ → ℝ} (hF : ContDiff ℝ ∞ F) (v : ℝ × ℝ) :
    ContDiff ℝ ∞ (fun q => fderiv ℝ F q v) :=
  (hF.fderiv_right (m := ∞) (by simp)).clm_apply contDiff_const

lemma pY_smooth {F : ℝ × ℝ → ℝ} (hF : ContDiff ℝ ∞ F) : ContDiff ℝ ∞ (pY F) := by
  rw [pY_eq hF]; exact apply_smooth hF _

lemma pH_smooth {F : ℝ × ℝ → ℝ} (hF : ContDiff ℝ ∞ F) : ContDiff ℝ ∞ (pH F) := by
  rw [pH_eq hF]; exact apply_smooth hF _

/-- the derivative of `q ↦ F'(q) v` along a path -/
lemma deriv_path {F : ℝ × ℝ → ℝ} (hF : ContDiff ℝ ∞ F) (v w : ℝ × ℝ) (γ : ℝ → ℝ × ℝ) (t0 : ℝ)
    (hγ : HasDerivAt γ w t0) :
    deriv (fun t : ℝ => fderiv ℝ F (γ t) v) t0 = fderiv ℝ (fderiv ℝ F) (γ t0) w v := by
  have hd : HasFDerivAt (fderiv ℝ F) (fderiv ℝ (fderiv ℝ F) (γ t0)) (γ t0) :=
    (((hF.fderiv_right (m := ∞) (by simp)).differentiable (by simp)) (γ t0)).hasFDerivAt
  have h1 : HasFDerivAt (fun p => fderiv ℝ F p v) ((fderiv ℝ F (γ t0)).comp 0 +
      (fderiv ℝ (fderiv ℝ F) (γ t0)).flip v) (γ t0) := hd.clm_apply (hasFDerivAt_const v _)
  have h3 := h1.comp_hasDerivAt t0 hγ
  rw [show (fun t : ℝ => fderiv ℝ F (γ t) v) = (fun p => fderiv ℝ F p v) ∘ γ from rfl, h3.deriv]
  simp

/-- Schwarz: `∂_h ∂_y F = ∂_y ∂_h F` -/
lemma swap1 {F : ℝ × ℝ → ℝ} (hF : ContDiff ℝ ∞ F) (y h : ℝ) :
    deriv (fun h' => pY F (y, h')) h = deriv (fun y' => pH F (y', h)) y := by
  have hs : IsSymmSndFDerivAt ℝ F (y, h) :=
    hF.contDiffAt.isSymmSndFDerivAt (by simp)
  rw [pY_eq hF, pH_eq hF]
  have h1 : HasDerivAt (fun h' : ℝ => (y, h')) ((0 : ℝ), (1 : ℝ)) h :=
    (hasDerivAt_const h y).prodMk (hasDerivAt_id h)
  have h2 : HasDerivAt (fun y' : ℝ => (y', h)) ((1 : ℝ), (0 : ℝ)) y :=
    (hasDerivAt_id y).prodMk (hasDerivAt_const y h)
  have e1 := deriv_path hF ((1 : ℝ), (0 : ℝ)) _ _ h h1
  have e2 := deriv_path hF ((0 : ℝ), (1 : ℝ)) _ _ y h2
  rw [e1, e2]
  exact hs _ _

lemma swapP (m : ℕ) : ∀ {F : ℝ × ℝ → ℝ}, ContDiff ℝ ∞ F → ∀ y h : ℝ,
    deriv (fun h' => iteratedDeriv m (fun y' => F (y', h')) y) h =
      iteratedDeriv m (fun y' => pH F (y', h)) y := by
  induction m with
  | zero => intro F hF y h; simp [pH]
  | succ m ih =>
    intro F hF y h
    have e1 : ∀ h', iteratedDeriv (m + 1) (fun y' => F (y', h')) y =
        iteratedDeriv m (fun y' => pY F (y', h')) y := by
      intro h'; rw [iteratedDeriv_succ']; rfl
    simp_rw [e1]
    rw [ih (pY_smooth hF) y h, iteratedDeriv_succ']
    congr 1
    funext y'
    simp only [pH]
    exact swap1 hF y' h

lemma swapS (i : ℕ) : ∀ {F : ℝ × ℝ → ℝ}, ContDiff ℝ ∞ F → ∀ (m : ℕ) (y h : ℝ),
    iteratedDeriv i (fun h' => iteratedDeriv m (fun y' => F (y', h')) y) h =
      iteratedDeriv m (fun y' => iteratedDeriv i (fun h' => F (y', h')) h) y := by
  induction i with
  | zero => intro F hF m y h; simp
  | succ i ih =>
    intro F hF m y h
    rw [iteratedDeriv_succ']
    have e1 : deriv (fun h' => iteratedDeriv m (fun y' => F (y', h')) y) =
        fun h' => iteratedDeriv m (fun y' => pH F (y', h')) y := by
      funext h'; exact swapP m hF y h'
    rw [e1, ih (pH_smooth hF) m y h]
    congr 1
    funext y'
    rw [iteratedDeriv_succ']
    rfl

/-- `e^{2y}/(1+e^y)` -/
noncomputable def Bre (y : ℝ) : ℝ := Real.exp (2 * y) / (1 + Real.exp y)
/-- `log(1+e^y)` -/
noncomputable def lre (y : ℝ) : ℝ := Real.log (1 + Real.exp y)

lemma Acal_eq (y h : ℝ) : Q02.Acal y h = Bre y * Real.exp ((-lre y) * h) := by
  unfold Q02.Acal Bre lre
  have hp : 0 < 1 + Real.exp y := by positivity
  rw [Real.rpow_def_of_pos hp]
  have : Real.exp (Real.log (1 + Real.exp y) * (-1 - h)) =
      Real.exp (-Real.log (1 + Real.exp y)) * Real.exp (-Real.log (1 + Real.exp y) * h) := by
    rw [← Real.exp_add]; ring_nf
  rw [this, Real.exp_neg, Real.exp_log hp]
  field_simp

lemma Acal_smooth : ContDiff ℝ ∞ (fun q : ℝ × ℝ => Q02.Acal q.1 q.2) := by
  have e : (fun q : ℝ × ℝ => Q02.Acal q.1 q.2) = fun q : ℝ × ℝ =>
      Real.exp (2 * q.1) * Real.exp (Real.log (1 + Real.exp q.1) * (-1 - q.2)) := by
    funext q
    unfold Q02.Acal
    rw [Real.rpow_def_of_pos (by positivity)]
  rw [e]
  have hl : ContDiff ℝ ∞ (fun q : ℝ × ℝ => Real.log (1 + Real.exp q.1)) :=
    (contDiff_const.add (contDiff_fst.exp)).log (fun q => by positivity)
  exact (contDiff_const.mul contDiff_fst).exp.mul ((hl.mul (contDiff_const.sub contDiff_snd)).exp)

lemma iter_h_Acal (y : ℝ) (i : ℕ) :
    iteratedDeriv i (fun h => Q02.Acal y h) 0 = Bre y * (-lre y) ^ i := by
  simp_rw [Acal_eq]
  rw [iteratedDeriv_const_mul_field, iteratedDeriv_exp_const_mul]
  simp

/-- `m! [y^m h^i] 𝓐 = ∂_y^m (B(y) (-λ(y))^i)(0) / i!` -/
theorem Acoef_eq (m i : ℕ) :
    (m.factorial : ℝ) * Q02.Acoef m i =
      iteratedDeriv m (fun y => Bre y * (-lre y) ^ i) 0 / (i.factorial : ℝ) := by
  unfold Q02.Acoef
  have := swapS i Acal_smooth m 0 0
  simp only at this
  rw [this]
  simp_rw [iter_h_Acal]
  have hm : (m.factorial : ℝ) ≠ 0 := by exact_mod_cast Nat.factorial_ne_zero m
  field_simp

/-! ### Part 4: the holomorphic extension to `|z| < π` -/

/-- `B(z) = e^{2z}/(1+e^z)` -/
noncomputable def Bc (z : ℂ) : ℂ := Complex.exp (2 * z) / (1 + Complex.exp z)
/-- `λ(z) = z/2 + Log(2 cosh(z/2))`, a branch of `log(1+e^z)` on `|Im z| < π` -/
noncomputable def lc (z : ℂ) : ℂ := z / 2 + Complex.log (2 * Complex.cosh (z / 2))
/-- `B(z) (-λ(z))^i` -/
noncomputable def bpart (i : ℕ) (z : ℂ) : ℂ := Bc z * (-lc z) ^ i
/-- `E_i = B(-λ)^i - S_+ (-L_+)^i - S_- (-L_-)^i` -/
noncomputable def E (i : ℕ) (z : ℂ) : ℂ := bpart i z - gpart ip i z - gpart (-ip) i z

lemma cosh_decomp (w : ℂ) : Complex.cosh w =
    ((Real.cosh w.re * Real.cos w.im : ℝ) : ℂ) + ((Real.sinh w.re * Real.sin w.im : ℝ) : ℂ) * Complex.I := by
  conv_lhs => rw [← Complex.re_add_im w]
  rw [Complex.cosh_add, Complex.cosh_mul_I, Complex.sinh_mul_I, ← Complex.ofReal_cosh,
    ← Complex.ofReal_sinh, ← Complex.ofReal_cos, ← Complex.ofReal_sin]
  push_cast; ring

lemma cosh_re_eq (w : ℂ) : (Complex.cosh w).re = Real.cosh w.re * Real.cos w.im := by
  rw [cosh_decomp]
  simp only [Complex.add_re, Complex.ofReal_re, Complex.mul_re, Complex.I_re, Complex.I_im,
    Complex.ofReal_im]
  ring

lemma norm_ip : ‖ip‖ = Real.pi := by
  rw [ip, norm_mul, Complex.norm_I, one_mul, Complex.norm_real, Real.norm_eq_abs,
    abs_of_pos Real.pi_pos]

lemma abs_im_lt {z : ℂ} (hz : ‖z‖ < Real.pi) : |z.im| < Real.pi :=
  lt_of_le_of_lt (Complex.abs_im_le_norm z) hz

lemma two_cosh_slit {z : ℂ} (hz : |z.im| < Real.pi) :
    2 * Complex.cosh (z / 2) ∈ Complex.slitPlane := by
  rw [Complex.mem_slitPlane_iff]
  left
  have h1 : (2 * Complex.cosh (z / 2)).re = 2 * (Complex.cosh (z / 2)).re := by simp
  rw [h1, cosh_re_eq]
  have hre : (z / 2).im = z.im / 2 := by simp
  rw [hre]
  have hc : 0 < Real.cos (z.im / 2) := by
    apply Real.cos_pos_of_mem_Ioo
    constructor <;> [skip; skip] <;> rw [abs_lt] at hz <;> linarith
  have := Real.cosh_pos (z / 2).re
  positivity

lemma one_add_exp_ne' {z : ℂ} (hz : ‖z‖ < Real.pi) : 1 + Complex.exp z ≠ 0 := by
  have hpi := Real.pi_pos
  apply Q02.one_add_exp_ne (by linarith)
  · intro h; rw [h, Q02.norm_ipi] at hz; linarith
  · intro h; rw [h, norm_neg, Q02.norm_ipi] at hz; linarith

lemma lc_analyticAt {z : ℂ} (hz : |z.im| < Real.pi) : AnalyticAt ℂ lc z := by
  unfold lc
  have h1 : AnalyticAt ℂ (fun z : ℂ => 2 * Complex.cosh (z / 2)) z := by
    have : AnalyticAt ℂ (fun z : ℂ => Complex.cosh (z / 2)) z :=
      (Complex.analyticAt_cosh).comp (analyticAt_id.div_const)
    exact analyticAt_const.mul this
  exact (analyticAt_id.div_const).add (h1.clog (two_cosh_slit hz))

lemma Bc_analyticAt {z : ℂ} (hz : ‖z‖ < Real.pi) : AnalyticAt ℂ Bc z := by
  unfold Bc
  have h1 : AnalyticAt ℂ (fun z : ℂ => Complex.exp (2 * z)) z :=
    (by fun_prop : Differentiable ℂ fun z : ℂ => Complex.exp (2 * z)).analyticAt z
  have h2 : AnalyticAt ℂ (fun z : ℂ => 1 + Complex.exp z) z :=
    (by fun_prop : Differentiable ℂ fun z : ℂ => 1 + Complex.exp z).analyticAt z
  exact h1.div h2 (one_add_exp_ne' hz)

lemma bpart_analyticAt {z : ℂ} (hz : ‖z‖ < Real.pi) (i : ℕ) : AnalyticAt ℂ (bpart i) z :=
  (Bc_analyticAt hz).mul ((lc_analyticAt (abs_im_lt hz)).neg.pow i)

lemma lre_eq (y : ℝ) : lre y = y / 2 + Real.log (2 * Real.cosh (y / 2)) := by
  unfold lre
  have h : 1 + Real.exp y = Real.exp (y / 2) * (2 * Real.cosh (y / 2)) := by
    rw [Real.cosh_eq]
    have e1 : Real.exp y = Real.exp (y / 2) * Real.exp (y / 2) := by rw [← Real.exp_add]; ring_nf
    have e2 : Real.exp (y / 2) * Real.exp (-(y / 2)) = 1 := by rw [← Real.exp_add]; simp
    rw [e1]; linear_combination -e2
  rw [h, Real.log_mul (by positivity) (by positivity), Real.log_exp]

lemma lc_real (y : ℝ) : lc (y : ℂ) = ((lre y : ℝ) : ℂ) := by
  rw [lre_eq, lc]
  have hc : 0 ≤ 2 * Real.cosh (y / 2) := by positivity
  push_cast
  rw [← Complex.ofReal_ofNat, ← Complex.ofReal_div, ← Complex.ofReal_cosh, ← Complex.ofReal_mul,
    ← Complex.ofReal_log hc]

lemma bpart_real (i : ℕ) (y : ℝ) : bpart i (y : ℂ) = ((Bre y * (-lre y) ^ i : ℝ) : ℂ) := by
  rw [bpart, lc_real, Bc, Bre]
  push_cast
  ring

lemma bpart_analyticOnNhd (i : ℕ) : AnalyticOnNhd ℂ (bpart i) (Metric.ball 0 Real.pi) :=
  fun z hz => bpart_analyticAt (by simpa using hz) i

lemma bpart_iter_real (i k : ℕ) : ∀ y : ℝ, |y| < Real.pi →
    iteratedDeriv k (fun y : ℝ => Bre y * (-lre y) ^ i) y = (iteratedDeriv k (bpart i) (y : ℂ)).re := by
  have han := bpart_analyticOnNhd i
  induction k with
  | zero =>
    intro y _
    simp only [iteratedDeriv_zero]
    rw [bpart_real, Complex.ofReal_re]
  | succ k ih =>
    intro y hy
    have hopen : IsOpen {t : ℝ | |t| < Real.pi} := isOpen_lt continuous_abs continuous_const
    have heq : iteratedDeriv k (fun y : ℝ => Bre y * (-lre y) ^ i) =ᶠ[𝓝 y]
        fun t : ℝ => (iteratedDeriv k (bpart i) (t : ℂ)).re := by
      filter_upwards [hopen.mem_nhds hy] with t ht using ih t ht
    rw [iteratedDeriv_succ, heq.deriv_eq]
    have hmem : (y : ℂ) ∈ Metric.ball (0 : ℂ) Real.pi := by
      rw [Metric.mem_ball, dist_zero_right, Complex.norm_real, Real.norm_eq_abs]; exact hy
    have hdk : DifferentiableAt ℂ (iteratedDeriv k (bpart i)) (y : ℂ) := by
      rw [iteratedDeriv_eq_iterate]
      exact ((han.iterated_deriv k) (y : ℂ) hmem).differentiableAt
    have := hdk.hasDerivAt.real_of_complex
    rw [this.deriv, iteratedDeriv_succ]

/-- the reduction: `T_k - T_k^{(1)} = Σ_m Re(∂^m E_{k-2-m}(0)) / (k-2-m)!` -/
theorem reduction (k : ℕ) :
    Q02.T k - Q02.T1 k = ∑ m ∈ Finset.range (k - 1),
      (iteratedDeriv m (E (k - 2 - m)) 0).re / ((k - 2 - m).factorial : ℝ) := by
  unfold Q02.T Q02.T1
  have h0 : ((0 : ℝ) : ℂ) = 0 := by simp
  have hT : ∀ m ∈ Finset.range (k - 1), (m.factorial : ℝ) * Q02.Acoef m (k - 2 - m) =
      (iteratedDeriv m (bpart (k - 2 - m)) 0).re / ((k - 2 - m).factorial : ℝ) := by
    intro m _
    rw [Acoef_eq, bpart_iter_real _ _ 0 (by simp [Real.pi_pos]), h0]
  rw [Finset.sum_congr rfl hT]
  have hQ : ∀ m ∈ Finset.range (k - 1),
      iteratedDeriv (k - 2 - m)
        (fun h : ℂ => ((Complex.I * Real.pi) ^ (m + 1))⁻¹ *
          Complex.exp (-(h * Complex.log (Complex.I * Real.pi))) *
          (ascPochhammer ℂ m).eval (1 + h)) 0 / ((k - 2 - m).factorial : ℂ) =
      iteratedDeriv m (gpart ip (k - 2 - m)) 0 / ((k - 2 - m).factorial : ℂ) := by
    intro m _
    rw [gpart_eq_Qf ip ip_slit]
    rfl
  rw [Finset.sum_congr rfl hQ, Complex.re_sum, Finset.mul_sum, ← Finset.sum_sub_distrib]
  refine Finset.sum_congr rfl (fun m _ => ?_)
  set i := k - 2 - m
  have h0' : ‖(0 : ℂ)‖ < Real.pi := by simp [Real.pi_pos]
  have hb : ContDiffAt ℂ m (bpart i) 0 := (bpart_analyticAt h0' i).contDiffAt.of_le le_top
  have hg1 : ContDiffAt ℂ m (gpart ip i) 0 :=
    (gpart_analyticAt ip 0 (by simpa using ip_slit) i).contDiffAt.of_le le_top
  have hg2 : ContDiffAt ℂ m (gpart (-ip) i) 0 :=
    (gpart_analyticAt (-ip) 0 (by simpa using nip_slit) i).contDiffAt.of_le le_top
  have hE : iteratedDeriv m (E i) 0 = iteratedDeriv m (bpart i) 0 -
      iteratedDeriv m (gpart ip i) 0 - iteratedDeriv m (gpart (-ip) i) 0 := by
    have e : E i = (bpart i - gpart ip i) - gpart (-ip) i := by funext z; simp [E]
    have h1 : ContDiffAt ℂ m (bpart i - gpart ip i) 0 := hb.sub hg1
    rw [e, iteratedDeriv_sub h1 hg2, iteratedDeriv_sub hb hg1]
  rw [hE, gpart_conj]
  have hf : ((i.factorial : ℕ) : ℂ) = ((i.factorial : ℝ) : ℂ) := by simp
  rw [hf, Complex.div_ofReal_re]
  simp only [Complex.sub_re, Complex.conj_re]
  ring

lemma ne_ip {z : ℂ} (hz : ‖z‖ < Real.pi) : z ≠ ip := by
  intro h; rw [h, norm_ip] at hz; linarith

lemma ne_nip {z : ℂ} (hz : ‖z‖ < Real.pi) : z ≠ -ip := by
  intro h; rw [h, norm_neg, norm_ip] at hz; linarith

lemma ip_sq : ip ^ 2 = -(Real.pi : ℂ) ^ 2 := by
  rw [ip, mul_pow, Complex.I_sq]; ring

lemma E0_eq_gR {z : ℂ} (hz : ‖z‖ < Real.pi) : E 0 z = Q02.gR z := by
  have h1 : ip - z ≠ 0 := sub_ne_zero.mpr (ne_ip hz).symm
  have h2 : -ip - z ≠ 0 := sub_ne_zero.mpr (ne_nip hz).symm
  have h3 : z ^ 2 + (Real.pi : ℂ) ^ 2 ≠ 0 := by
    have : z ^ 2 + (Real.pi : ℂ) ^ 2 = (ip - z) * (-ip - z) := by
      ring_nf; rw [ip_sq]; ring
    rw [this]; exact mul_ne_zero h1 h2
  simp only [E, bpart, gpart, pow_zero, mul_one, Q02.gR, Bc]
  have e : (ip - z)⁻¹ + (-ip - z)⁻¹ = -(2 * z / (z ^ 2 + (Real.pi : ℂ) ^ 2)) := by
    have : z ^ 2 + (Real.pi : ℂ) ^ 2 = (ip - z) * (-ip - z) := by
      ring_nf; rw [ip_sq]; ring
    rw [this]; field_simp; ring
  linear_combination -e

lemma E0_eq_g2 {z : ℂ} (hz : ‖z‖ < Real.pi) : E 0 z = Q02.g2 z := by
  rw [E0_eq_gR hz, Q02.g2,
    Function.update_of_ne (show z ≠ -(Complex.I * (Real.pi : ℂ)) from ne_nip hz), Q02.g1,
    Function.update_of_ne (show z ≠ Complex.I * (Real.pi : ℂ) from ne_ip hz)]

/-! ### Part 5a: local estimates at `c = ±iπ` -/

lemma exp_c {c : ℂ} (hc : c = ip ∨ c = -ip) : Complex.exp c = -1 := by
  rcases hc with rfl | rfl
  · exact Q02.exp_ipi
  · exact Q02.exp_neg_ipi

lemma norm_c {c : ℂ} (hc : c = ip ∨ c = -ip) : ‖c‖ = Real.pi := by
  rcases hc with rfl | rfl
  · exact norm_ip
  · rw [norm_neg, norm_ip]

lemma exp3_bound {w : ℂ} (hw : ‖w‖ ≤ 1) : ‖Complex.exp w - 1 - w - w ^ 2 / 2‖ ≤ ‖w‖ ^ 3 := by
  have h := Complex.exp_bound hw (n := 3) (by norm_num)
  have e : ∑ m ∈ Finset.range 3, w ^ m / (m.factorial : ℂ) = 1 + w + w ^ 2 / 2 := by
    simp [Finset.sum_range_succ, Nat.factorial]
  rw [e] at h
  have e2 : Complex.exp w - 1 - w - w ^ 2 / 2 = Complex.exp w - (1 + w + w ^ 2 / 2) := by ring
  rw [e2]
  refine h.trans ?_
  have : (0 : ℝ) ≤ ‖w‖ ^ 3 := by positivity
  norm_num
  nlinarith

lemma exp2_bound {w : ℂ} (hw : ‖w‖ ≤ 1) : ‖Complex.exp w - 1 - w‖ ≤ ‖w‖ ^ 2 :=
  Complex.norm_exp_sub_one_sub_id_le hw

lemma one_sub_exp_lower {u : ℂ} (hu : ‖u‖ ≤ 1 / 2) : ‖u‖ / 2 ≤ ‖1 - Complex.exp u‖ := by
  have h := exp2_bound (w := u) (by linarith)
  have e : 1 - Complex.exp u = -u - (Complex.exp u - 1 - u) := by ring
  rw [e]
  have h2 : ‖u‖ - ‖Complex.exp u - 1 - u‖ ≤ ‖-u - (Complex.exp u - 1 - u)‖ := by
    have := norm_sub_norm_le (-u) (Complex.exp u - 1 - u)
    rwa [norm_neg] at this
  have h3 : ‖u‖ ^ 2 ≤ ‖u‖ / 2 := by nlinarith [norm_nonneg u]
  linarith

/-- `B(z) - 1/(c-z) + 3/2 = O(z - c)` -/
lemma Bc_local {c : ℂ} (hc : c = ip ∨ c = -ip) {z : ℂ} (hu : ‖z - c‖ ≤ 1 / 2) (hz0 : z ≠ c) :
    ‖Bc z - (c - z)⁻¹ + 3 / 2‖ ≤ 13 * ‖z - c‖ := by
  set u := z - c with hudef
  have hz : z = c + u := by rw [hudef]; ring
  have hu0 : u ≠ 0 := sub_ne_zero.mpr hz0
  have hn0 : 0 < ‖u‖ := norm_pos_iff.mpr hu0
  have hec := exp_c hc
  have hlow := one_sub_exp_lower hu
  have hne : 1 - Complex.exp u ≠ 0 := by
    intro h; rw [h, norm_zero] at hlow; linarith
  have hB : Bc z = Complex.exp (2 * u) / (1 - Complex.exp u) := by
    rw [Bc, hz, mul_add, Complex.exp_add, Complex.exp_add, show (2 : ℂ) * c = c + c by ring,
      Complex.exp_add, hec]
    ring_nf
  set R1 := Complex.exp (2 * u) - 1 - 2 * u
  set R2 := Complex.exp u - 1 - u - u ^ 2 / 2
  set R3 := Complex.exp u - 1 - u
  have hR1 : ‖R1‖ ≤ 4 * ‖u‖ ^ 2 := by
    have := exp2_bound (w := 2 * u) (by rw [norm_mul]; norm_num; linarith)
    rw [norm_mul] at this; norm_num at this; nlinarith
  have hR2 : ‖R2‖ ≤ ‖u‖ ^ 3 := exp3_bound (by linarith)
  have hR3 : ‖R3‖ ≤ ‖u‖ ^ 2 := exp2_bound (by linarith)
  have key : Bc z - (c - z)⁻¹ + 3 / 2 = (u * R1 - R2 - 3 / 2 * u * R3) / (u * (1 - Complex.exp u)) := by
    rw [hB, show c - z = -u by rw [hudef]; ring]
    field_simp
    simp only [R1, R2, R3]
    ring
  rw [key, norm_div, norm_mul, div_le_iff₀ (by positivity)]
  have hN : ‖u * R1 - R2 - 3 / 2 * u * R3‖ ≤ 13 / 2 * ‖u‖ ^ 3 := by
    calc ‖u * R1 - R2 - 3 / 2 * u * R3‖ ≤ ‖u * R1‖ + ‖R2‖ + ‖3 / 2 * u * R3‖ := by
          refine (norm_sub_le _ _).trans ?_
          gcongr
          exact norm_sub_le _ _
      _ = ‖u‖ * ‖R1‖ + ‖R2‖ + 3 / 2 * ‖u‖ * ‖R3‖ := by
          rw [norm_mul, norm_mul, norm_mul]; norm_num
      _ ≤ ‖u‖ * (4 * ‖u‖ ^ 2) + ‖u‖ ^ 3 + 3 / 2 * ‖u‖ * ‖u‖ ^ 2 := by gcongr
      _ = 13 / 2 * ‖u‖ ^ 3 := by ring
  have : 13 / 2 * ‖u‖ ^ 3 ≤ 13 * ‖u‖ * (‖u‖ * ‖1 - Complex.exp u‖) := by
    have := mul_le_mul_of_nonneg_left hlow (by positivity : (0 : ℝ) ≤ 13 * ‖u‖ * ‖u‖)
    nlinarith
  linarith

/-- `λ'(z) = e^z/(1+e^z)` -/
noncomputable def lcd (z : ℂ) : ℂ := Complex.exp z / (1 + Complex.exp z)

lemma one_add_exp_eq (z : ℂ) :
    1 + Complex.exp z = Complex.exp (z / 2) * (2 * Complex.cosh (z / 2)) := by
  rw [Complex.cosh, mul_div_cancel₀ _ (two_ne_zero), mul_add, ← Complex.exp_add,
    ← Complex.exp_add]
  ring_nf
  simp [add_comm]

lemma lc_hasDerivAt {z : ℂ} (hz : |z.im| < Real.pi) : HasDerivAt lc (lcd z) z := by
  have hs := two_cosh_slit hz
  have hc0 : 2 * Complex.cosh (z / 2) ≠ 0 := Complex.slitPlane_ne_zero hs
  have h1 : HasDerivAt (fun z : ℂ => z / 2) (1 / 2) z := by
    simpa using (hasDerivAt_id z).div_const 2
  have h2a : HasDerivAt (fun z : ℂ => Complex.cosh (z / 2)) (Complex.sinh (z / 2) * (1 / 2)) z :=
    h1.ccosh
  have h2 : HasDerivAt (fun z : ℂ => 2 * Complex.cosh (z / 2)) (2 * (Complex.sinh (z / 2) * (1 / 2))) z :=
    h2a.const_mul 2
  have h3 : HasDerivAt (fun z : ℂ => z / 2 + Complex.log (2 * Complex.cosh (z / 2)))
      (1 / 2 + 2 * (Complex.sinh (z / 2) * (1 / 2)) / (2 * Complex.cosh (z / 2))) z :=
    h1.add (h2.clog hs)
  have hq : Complex.exp (z / 2) ≠ 0 := Complex.exp_ne_zero _
  have he : 1 + Complex.exp z ≠ 0 := by rw [one_add_exp_eq]; exact mul_ne_zero hq hc0
  have e1 : Complex.exp z = Complex.exp (z / 2) * Complex.exp (z / 2) := by
    rw [← Complex.exp_add]; ring_nf
  have hsinh : Complex.sinh (z / 2) = (Complex.exp (z / 2) - (Complex.exp (z / 2))⁻¹) / 2 := by
    rw [Complex.sinh, Complex.exp_neg]
  have hcosh : 2 * Complex.cosh (z / 2) = Complex.exp (z / 2) + (Complex.exp (z / 2))⁻¹ := by
    rw [Complex.cosh, Complex.exp_neg]; ring
  have hv : 1 / 2 + 2 * (Complex.sinh (z / 2) * (1 / 2)) / (2 * Complex.cosh (z / 2)) = lcd z := by
    rw [lcd, hcosh, hsinh, e1]
    rw [hcosh] at hc0
    rw [e1] at he
    generalize Complex.exp (z / 2) = q at hq hc0 he ⊢
    rw [show q + q⁻¹ = (1 + q * q) * q⁻¹ by field_simp; ring,
      show q - q⁻¹ = (q * q - 1) * q⁻¹ by field_simp]
    field_simp
    have he2 : 1 + q ^ 2 ≠ 0 := by rwa [sq]
    have : (1 + q ^ 2) * (1 + q ^ 2)⁻¹ = 1 := mul_inv_cancel₀ he2
    linear_combination (-1 : ℂ) * this
  rw [← hv]
  exact h3

/-- `λ'(z) - 1/(z - c) = O(1)` -/
lemma lcd_local {c : ℂ} (hc : c = ip ∨ c = -ip) {z : ℂ} (hu : ‖z - c‖ ≤ 1 / 2) (hz0 : z ≠ c) :
    ‖lcd z - (z - c)⁻¹‖ ≤ 3 := by
  set u := z - c with hudef
  have hz : z = c + u := by rw [hudef]; ring
  have hu0 : u ≠ 0 := sub_ne_zero.mpr hz0
  have hn0 : 0 < ‖u‖ := norm_pos_iff.mpr hu0
  have hec := exp_c hc
  have hlow := one_sub_exp_lower hu
  have hne : 1 - Complex.exp u ≠ 0 := by
    intro h; rw [h, norm_zero] at hlow; linarith
  have hL : lcd z = -Complex.exp u / (1 - Complex.exp u) := by
    rw [lcd, hz, Complex.exp_add, hec]; ring_nf
  set R2 := Complex.exp u - 1 - u - u ^ 2 / 2
  set R3 := Complex.exp u - 1 - u
  have hR2 : ‖R2‖ ≤ ‖u‖ ^ 3 := exp3_bound (by linarith)
  have hR3 : ‖R3‖ ≤ ‖u‖ ^ 2 := exp2_bound (by linarith)
  have key : lcd z - (z - c)⁻¹ = -(u ^ 2 / 2 + u * R3 - R2) / (u * (1 - Complex.exp u)) := by
    rw [hL]
    field_simp
    simp only [R2, R3]
    ring
  rw [key, norm_div, norm_neg, norm_mul, div_le_iff₀ (by positivity)]
  have hN : ‖u ^ 2 / 2 + u * R3 - R2‖ ≤ 3 / 2 * ‖u‖ ^ 2 := by
    calc ‖u ^ 2 / 2 + u * R3 - R2‖ ≤ ‖u ^ 2 / 2‖ + ‖u * R3‖ + ‖R2‖ := by
          refine (norm_sub_le _ _).trans ?_
          gcongr
          exact norm_add_le _ _
      _ = ‖u‖ ^ 2 / 2 + ‖u‖ * ‖R3‖ + ‖R2‖ := by
          rw [norm_mul, norm_div, norm_pow]; norm_num
      _ ≤ ‖u‖ ^ 2 / 2 + ‖u‖ * ‖u‖ ^ 2 + ‖u‖ ^ 3 := by gcongr
      _ ≤ 3 / 2 * ‖u‖ ^ 2 := by nlinarith [norm_nonneg u]
  have : 3 / 2 * ‖u‖ ^ 2 ≤ 3 * (‖u‖ * ‖1 - Complex.exp u‖) := by
    have := mul_le_mul_of_nonneg_left hlow (by positivity : (0 : ℝ) ≤ 3 * ‖u‖)
    nlinarith
  linarith

/-- `σ(u) = 2 sinh(u/2)/u` -/
noncomputable def sig (u : ℂ) : ℂ := 2 * Complex.sinh (u / 2) / u

lemma sinh_local {u : ℂ} (hu : ‖u‖ ≤ 1 / 2) : ‖2 * Complex.sinh (u / 2) - u‖ ≤ ‖u‖ ^ 3 / 4 := by
  have hv : ‖u / 2‖ ≤ 1 := by rw [norm_div]; norm_num; linarith
  have hv' : ‖-(u / 2)‖ ≤ 1 := by rwa [norm_neg]
  have h1 := exp3_bound hv
  have h2 := exp3_bound hv'
  have e : 2 * Complex.sinh (u / 2) - u = (Complex.exp (u / 2) - 1 - u / 2 - (u / 2) ^ 2 / 2) -
      (Complex.exp (-(u / 2)) - 1 - -(u / 2) - (-(u / 2)) ^ 2 / 2) := by
    rw [Complex.sinh]; ring
  rw [e]
  refine (norm_sub_le _ _).trans ?_
  rw [norm_neg] at h2
  have : ‖u / 2‖ ^ 3 = ‖u‖ ^ 3 / 8 := by rw [norm_div]; norm_num; ring
  linarith

lemma cosh_local {u : ℂ} (hu : ‖u‖ ≤ 1 / 2) : ‖Complex.cosh (u / 2) - 1‖ ≤ ‖u‖ ^ 2 / 4 := by
  have hv : ‖u / 2‖ ≤ 1 := by rw [norm_div]; norm_num; linarith
  have hv' : ‖-(u / 2)‖ ≤ 1 := by rwa [norm_neg]
  have h1 := exp2_bound hv
  have h2 := exp2_bound hv'
  have e : Complex.cosh (u / 2) - 1 = ((Complex.exp (u / 2) - 1 - u / 2) +
      (Complex.exp (-(u / 2)) - 1 - -(u / 2))) / 2 := by
    rw [Complex.cosh]; ring
  rw [e, norm_div]
  have h3 := norm_add_le (Complex.exp (u / 2) - 1 - u / 2) (Complex.exp (-(u / 2)) - 1 - -(u / 2))
  rw [norm_neg] at h2
  have : ‖u / 2‖ ^ 2 = ‖u‖ ^ 2 / 4 := by rw [norm_div]; norm_num; ring
  have h2' : ‖(2 : ℂ)‖ = 2 := by norm_num
  rw [h2', div_le_iff₀ (by norm_num)]
  linarith

lemma sig_local {u : ℂ} (hu : ‖u‖ ≤ 1 / 2) (hu0 : u ≠ 0) : ‖sig u - 1‖ ≤ ‖u‖ ^ 2 / 4 := by
  have hn : 0 < ‖u‖ := norm_pos_iff.mpr hu0
  have e : sig u - 1 = (2 * Complex.sinh (u / 2) - u) / u := by
    rw [sig]; field_simp
  rw [e, norm_div, div_le_iff₀ hn]
  have := sinh_local hu
  nlinarith

/-- the data of `c = ±iπ`: `c = sπ`, `cosh(c/2) = 0`, `sinh(c/2) = s`, `Log s = c/2` -/
lemma c_data {c : ℂ} (hc : c = ip ∨ c = -ip) : ∃ s : ℂ, (s = Complex.I ∨ s = -Complex.I) ∧
    c = s * Real.pi ∧ Complex.cosh (c / 2) = 0 ∧ Complex.sinh (c / 2) = s ∧
    Complex.log s = c / 2 := by
  rcases hc with rfl | rfl
  · refine ⟨Complex.I, Or.inl rfl, by rw [ip, mul_comm], ?_, ?_, ?_⟩
    · rw [show ip / 2 = ((Real.pi / 2 : ℝ) : ℂ) * Complex.I by rw [ip]; push_cast; ring,
        Complex.cosh_mul_I, ← Complex.ofReal_cos, Real.cos_pi_div_two]; simp
    · rw [show ip / 2 = ((Real.pi / 2 : ℝ) : ℂ) * Complex.I by rw [ip]; push_cast; ring,
        Complex.sinh_mul_I, ← Complex.ofReal_sin, Real.sin_pi_div_two]; simp
    · rw [Complex.log_I, ip]; ring
  · refine ⟨-Complex.I, Or.inr rfl, by rw [ip]; ring, ?_, ?_, ?_⟩
    · rw [show -ip / 2 = ((-(Real.pi / 2) : ℝ) : ℂ) * Complex.I by rw [ip]; push_cast; ring,
        Complex.cosh_mul_I, ← Complex.ofReal_cos, Real.cos_neg, Real.cos_pi_div_two]; simp
    · rw [show -ip / 2 = ((-(Real.pi / 2) : ℝ) : ℂ) * Complex.I by rw [ip]; push_cast; ring,
        Complex.sinh_mul_I, ← Complex.ofReal_sin, Real.sin_neg, Real.sin_pi_div_two]; simp
    · rw [Complex.log_neg_I, ip]; ring

lemma arg_lt_of_re_pos {w : ℂ} (hw : 0 < w.re) : |Complex.arg w| < Real.pi / 2 :=
  Complex.abs_arg_lt_pi_div_two_iff.mpr (Or.inl hw)

/-- `λ(z) = Log(c - z) + (z - c)/2 + Log σ(z - c)` near `c` -/
lemma lc_local_eq {c : ℂ} (hc : c = ip ∨ c = -ip) {z : ℂ} (hz : ‖z‖ < Real.pi)
    (hu : ‖z - c‖ ≤ 1 / 2) :
    lc z = Complex.log (c - z) + (z - c) / 2 + Complex.log (sig (z - c)) := by
  obtain ⟨s, hs, hcs, hch, hsh, hls⟩ := c_data hc
  have hpi := Real.pi_pos
  set u := z - c with hudef
  have hzc : z = c + u := by rw [hudef]; ring
  have hu0 : u ≠ 0 := by
    intro h; rw [h, norm_zero] at hu
    have : z = c := by rw [hzc, h, add_zero]
    rw [this, norm_c hc] at hz; linarith
  have hs2 : s * s = -1 := by rcases hs with rfl | rfl <;> simp
  have hs0 : s ≠ 0 := by rcases hs with rfl | rfl <;> simp
  -- `2 cosh(z/2) = s u σ(u)`
  have hcosh : 2 * Complex.cosh (z / 2) = (s * u) * sig u := by
    rw [hzc, show (c + u) / 2 = c / 2 + u / 2 by ring, Complex.cosh_add, hch, hsh, sig]
    field_simp
    ring
  -- `Re(s u) > 0`
  have hre : 0 < (s * u).re := by
    have hsu : s * u = s * z + Real.pi := by
      rw [hudef, hcs]; ring_nf; rw [show s ^ 2 = s * s by ring, hs2]; ring
    have him := abs_im_lt hz
    rw [abs_lt] at him
    rw [hsu]
    rcases hs with rfl | rfl <;> simp <;> linarith
  have hsig := sig_local hu hu0
  have hsigre : 0 < (sig u).re := by
    have h1 : |(sig u - 1).re| ≤ ‖sig u - 1‖ := Complex.abs_re_le_norm _
    have h2 : ‖u‖ ^ 2 / 4 ≤ 1 / 16 := by nlinarith [norm_nonneg u]
    have : (sig u - 1).re = (sig u).re - 1 := by simp
    rw [this] at h1
    rw [abs_le] at h1
    linarith
  have hsu0 : s * u ≠ 0 := mul_ne_zero hs0 hu0
  have hsig0 : sig u ≠ 0 := by intro h; rw [h] at hsigre; simp at hsigre
  have ha1 := arg_lt_of_re_pos hre
  have ha2 := arg_lt_of_re_pos hsigre
  have hl1 : Complex.log ((s * u) * sig u) = Complex.log (s * u) + Complex.log (sig u) := by
    apply Complex.log_mul hsu0 hsig0
    rw [abs_lt] at ha1 ha2
    constructor <;> linarith
  have hargs : |Complex.arg s| = Real.pi / 2 := by
    rcases hs with rfl | rfl
    · rw [Complex.arg_I, abs_of_pos (by positivity)]
    · rw [Complex.arg_neg_I, abs_neg, abs_of_pos (by positivity)]
  have hl2 : Complex.log (c - z) = Complex.log s + Complex.log (s * u) := by
    have e : c - z = s * (s * u) := by rw [← mul_assoc, hs2, hudef]; ring
    rw [e]
    apply Complex.log_mul hs0 hsu0
    rw [abs_lt] at ha1
    have := abs_le.mp hargs.le
    constructor <;> linarith
  rw [lc, hcosh, hl1, hl2, hls]
  ring

/-- `d_c(z) = λ(z) - Log(c - z)` -/
noncomputable def dfun (c z : ℂ) : ℂ := lc z - Complex.log (c - z)

lemma c_sub_slit {c : ℂ} (hc : c = ip ∨ c = -ip) {z : ℂ} (hz : ‖z‖ < Real.pi) :
    c - z ∈ Complex.slitPlane := by
  rw [Complex.mem_slitPlane_iff]; right
  have him := abs_im_lt hz
  rw [abs_lt] at him
  rcases hc with rfl | rfl <;> simp [ip] <;> intro h <;> linarith

lemma dfun_hasDerivAt {c : ℂ} (hc : c = ip ∨ c = -ip) {z : ℂ} (hz : ‖z‖ < Real.pi) :
    HasDerivAt (dfun c) (lcd z + (c - z)⁻¹) z := by
  have h1 := lc_hasDerivAt (abs_im_lt hz)
  have h2 : HasDerivAt (fun z : ℂ => c - z) (-1) z := by
    simpa using (hasDerivAt_id z).const_sub c
  have h3 : HasDerivAt (fun z : ℂ => Complex.log (c - z)) (-1 / (c - z)) z :=
    h2.clog (c_sub_slit hc hz)
  have h4 : HasDerivAt (fun z : ℂ => lc z - Complex.log (c - z)) (lcd z - -1 / (c - z)) z :=
    h1.sub h3
  have e : lcd z - -1 / (c - z) = lcd z + (c - z)⁻¹ := by ring
  rw [e] at h4
  exact h4

lemma d_local {c : ℂ} (hc : c = ip ∨ c = -ip) {z : ℂ} (hz : ‖z‖ < Real.pi)
    (hu : ‖z - c‖ < 1 / 2) :
    ‖dfun c z‖ ≤ ‖z - c‖ ∧ ‖dfun c z - (z - c) * (lcd z + (c - z)⁻¹)‖ ≤ ‖z - c‖ ^ 2 := by
  have hpi := Real.pi_pos
  set u := z - c with hudef
  have hu0 : u ≠ 0 := by
    intro h
    have : z = c := by rw [← sub_eq_zero]; exact h
    rw [this, norm_c hc] at hz; linarith
  have hn : 0 < ‖u‖ := norm_pos_iff.mpr hu0
  -- the local form on an open set
  set V := Metric.ball (0 : ℂ) Real.pi ∩ Metric.ball c (1 / 2) with hV
  have hVo : IsOpen V := Metric.isOpen_ball.inter Metric.isOpen_ball
  have hzV : z ∈ V := ⟨by simpa using hz, by rw [Metric.mem_ball, dist_eq_norm]; exact hu⟩
  set φ : ℂ → ℂ := fun w => (w - c) / 2 + Complex.log (2 * Complex.sinh ((w - c) / 2) / (w - c))
    with hφ
  have hEq : dfun c =ᶠ[𝓝 z] φ := by
    filter_upwards [hVo.mem_nhds hzV] with w hw
    have hw1 : ‖w‖ < Real.pi := by simpa using hw.1
    have hw2 : ‖w - c‖ ≤ 1 / 2 := by
      have := hw.2; rw [Metric.mem_ball, dist_eq_norm] at this; linarith
    rw [dfun, lc_local_eq hc hw1 hw2, hφ, sig]; ring
  -- the derivative of `φ`
  have hsig := sig_local hu.le hu0
  have hsigre : 0 < (sig u).re := by
    have h1 : |(sig u - 1).re| ≤ ‖sig u - 1‖ := Complex.abs_re_le_norm _
    have h2 : ‖u‖ ^ 2 / 4 ≤ 1 / 16 := by nlinarith [norm_nonneg u]
    have : (sig u - 1).re = (sig u).re - 1 := by simp
    rw [this, abs_le] at h1
    linarith
  have hsig0 : sig u ≠ 0 := by intro h; rw [h] at hsigre; simp at hsigre
  have hslit : 2 * Complex.sinh (u / 2) / u ∈ Complex.slitPlane :=
    Complex.mem_slitPlane_iff.mpr (Or.inl hsigre)
  have hw : HasDerivAt (fun w : ℂ => w - c) 1 z := (hasDerivAt_id' z).sub_const c
  have hw2 : HasDerivAt (fun w : ℂ => (w - c) / 2) (1 / 2) z := hw.div_const 2
  have hsh0 : HasDerivAt (fun w : ℂ => Complex.sinh ((w - c) / 2))
      (Complex.cosh (u / 2) * (1 / 2)) z := hw2.csinh
  have hsh : HasDerivAt (fun w : ℂ => 2 * Complex.sinh ((w - c) / 2))
      (2 * (Complex.cosh (u / 2) * (1 / 2))) z := hsh0.const_mul 2
  have hsg : HasDerivAt (fun w : ℂ => 2 * Complex.sinh ((w - c) / 2) / (w - c))
      ((2 * (Complex.cosh (u / 2) * (1 / 2)) * u - 2 * Complex.sinh (u / 2) * 1) / u ^ 2) z :=
    hsh.div hw hu0
  have hφd : HasDerivAt φ (1 / 2 + ((2 * (Complex.cosh (u / 2) * (1 / 2)) * u -
      2 * Complex.sinh (u / 2) * 1) / u ^ 2) / (2 * Complex.sinh (u / 2) / u)) z :=
    hw2.add (hsg.clog hslit)
  have hd1 := dfun_hasDerivAt hc hz
  have hd2 := hφd.congr_of_eventuallyEq hEq
  have huniq : lcd z + (c - z)⁻¹ = 1 / 2 + ((2 * (Complex.cosh (u / 2) * (1 / 2)) * u -
      2 * Complex.sinh (u / 2) * 1) / u ^ 2) / sig u := hd1.unique hd2
  have hdz : dfun c z = u / 2 + Complex.log (sig u) := by
    rw [hEq.eq_of_nhds, hφ, sig]
  -- the bounds
  have hlog : ‖Complex.log (sig u)‖ ≤ 3 / 2 * ‖sig u - 1‖ := by
    have := Complex.norm_log_one_add_half_le_self (z := sig u - 1) (by nlinarith [norm_nonneg u])
    simpa using this
  have hsign : 15 / 16 ≤ ‖sig u‖ := by
    have := norm_sub_norm_le (1 : ℂ) (1 - sig u)
    have e : ‖1 - sig u‖ = ‖sig u - 1‖ := norm_sub_rev _ _
    simp only [sub_sub_cancel, norm_one] at this
    nlinarith [norm_nonneg u]
  have hcl := cosh_local hu.le
  constructor
  · rw [hdz]
    calc ‖u / 2 + Complex.log (sig u)‖ ≤ ‖u / 2‖ + ‖Complex.log (sig u)‖ := norm_add_le _ _
      _ ≤ ‖u‖ / 2 + 3 / 2 * (‖u‖ ^ 2 / 4) := by
          rw [norm_div]; norm_num; linarith
      _ ≤ ‖u‖ := by nlinarith
  · rw [huniq, hdz]
    have e : u / 2 + Complex.log (sig u) - u * (1 / 2 + ((2 * (Complex.cosh (u / 2) * (1 / 2)) * u -
        2 * Complex.sinh (u / 2) * 1) / u ^ 2) / sig u) =
        Complex.log (sig u) - (Complex.cosh (u / 2) - sig u) / sig u := by
      rw [sig]; field_simp; ring
    rw [e]
    have h1 : ‖Complex.cosh (u / 2) - sig u‖ ≤ ‖u‖ ^ 2 / 2 := by
      have : Complex.cosh (u / 2) - sig u = (Complex.cosh (u / 2) - 1) - (sig u - 1) := by ring
      rw [this]
      refine (norm_sub_le _ _).trans ?_
      linarith
    have h2 : ‖(Complex.cosh (u / 2) - sig u) / sig u‖ ≤ ‖u‖ ^ 2 / 2 / (15 / 16) := by
      rw [norm_div]
      exact div_le_div₀ (by positivity) h1 (by norm_num) hsign
    calc ‖Complex.log (sig u) - (Complex.cosh (u / 2) - sig u) / sig u‖
        ≤ ‖Complex.log (sig u)‖ + ‖(Complex.cosh (u / 2) - sig u) / sig u‖ := norm_sub_le _ _
      _ ≤ 3 / 2 * (‖u‖ ^ 2 / 4) + ‖u‖ ^ 2 / 2 / (15 / 16) := by linarith
      _ ≤ ‖u‖ ^ 2 := by nlinarith [norm_nonneg u]

/-! ### Part 5b: global estimates on `|z| < π` -/

lemma norm_cosh_sq (w : ℂ) :
    ‖Complex.cosh w‖ ^ 2 = Real.cos w.im ^ 2 + Real.sinh w.re ^ 2 := by
  rw [cosh_decomp, ← Complex.normSq_eq_norm_sq]
  rw [Complex.normSq_apply]
  simp only [Complex.add_re, Complex.ofReal_re, Complex.mul_re, Complex.I_re, Complex.I_im,
    Complex.ofReal_im, Complex.add_im, Complex.mul_im]
  have h1 := Real.cosh_sq w.re
  have h2 := Real.sin_sq_add_cos_sq w.im
  nlinarith [h1, h2]

lemma norm_exp_half_le {z : ℂ} (hz : ‖z‖ < Real.pi) : ‖Complex.exp (z / 2)‖ ≤ 8 := by
  rw [Complex.norm_exp]
  have h1 : (z / 2).re ≤ Real.pi / 2 := by
    have := Complex.re_le_norm z
    simp; linarith
  have h2 : Real.exp (Real.pi / 2) ≤ 8 := by
    have := Real.pi_lt_d2
    have h3 : Real.exp (Real.pi / 2) ≤ Real.exp 2 := Real.exp_le_exp.mpr (by linarith)
    have h4 := Real.exp_one_lt_d9
    have h5 : Real.exp 2 = Real.exp 1 * Real.exp 1 := by rw [← Real.exp_add]; norm_num
    nlinarith [Real.exp_pos 1]
  exact (Real.exp_le_exp.mpr h1).trans h2

lemma norm_cosh_half_le {z : ℂ} (hz : ‖z‖ < Real.pi) : ‖2 * Complex.cosh (z / 2)‖ ≤ 16 := by
  have h1 : ‖Complex.exp (-(z / 2))‖ ≤ 8 := by
    have := norm_exp_half_le (z := -z) (by rwa [norm_neg])
    rwa [neg_div] at this
  have h2 := norm_exp_half_le hz
  rw [Complex.cosh, mul_div_cancel₀ _ two_ne_zero]
  exact (norm_add_le _ _).trans (by linarith)

/-- `|2 cosh(z/2)| ≥ 2t/π` when `z` is at distance `≥ t` from `±iπ` -/
lemma norm_cosh_half_ge {z : ℂ} (hz : ‖z‖ < Real.pi) {t : ℝ} (ht : 0 < t)
    (ha : t ≤ ‖z - ip‖) (hb : t ≤ ‖z + ip‖) : 2 * t / Real.pi ≤ ‖2 * Complex.cosh (z / 2)‖ := by
  have hpi := Real.pi_pos
  set x := z.re
  set y := z.im
  have him := abs_im_lt hz
  have hsq := norm_cosh_sq (z / 2)
  have hre : (z / 2).re = x / 2 := by simp [x]
  have hime : (z / 2).im = y / 2 := by simp [y]
  rw [hre, hime] at hsq
  -- `cos(y/2) ≥ (π - |y|)/π`
  have hcos : (Real.pi - |y|) / Real.pi ≤ Real.cos (y / 2) := by
    have e : Real.cos (y / 2) = Real.sin ((Real.pi - |y|) / 2) := by
      rw [show (Real.pi - |y|) / 2 = Real.pi / 2 - |y| / 2 by ring, Real.sin_pi_div_two_sub]
      rcases abs_choice y with h | h <;> rw [h]
      rw [neg_div, Real.cos_neg]
    rw [e]
    have := Real.mul_le_sin (x := (Real.pi - |y|) / 2) (by linarith [abs_nonneg y]) (by linarith [abs_nonneg y])
    calc (Real.pi - |y|) / Real.pi = 2 / Real.pi * ((Real.pi - |y|) / 2) := by field_simp
      _ ≤ _ := this
  have hcos0 : 0 ≤ (Real.pi - |y|) / Real.pi := by
    apply div_nonneg _ hpi.le; linarith
  have hsinh : (x / 2) ^ 2 ≤ Real.sinh (x / 2) ^ 2 := by
    have h1 : |x / 2| ≤ |Real.sinh (x / 2)| := by
      rw [Real.abs_sinh]; exact Real.self_le_sinh_iff.mpr (abs_nonneg _)
    nlinarith [abs_nonneg (x / 2), sq_abs (x / 2), sq_abs (Real.sinh (x / 2))]
  -- `t² ≤ x² + (π - |y|)²`
  have ht2 : t ^ 2 ≤ x ^ 2 + (Real.pi - |y|) ^ 2 := by
    rcases le_or_gt 0 y with hy | hy
    · have e : ‖z - ip‖ ^ 2 = x ^ 2 + (Real.pi - |y|) ^ 2 := by
        rw [← Complex.normSq_eq_norm_sq, Complex.normSq_apply, abs_of_nonneg hy]
        simp [ip, x, y]; ring
      nlinarith
    · have e : ‖z + ip‖ ^ 2 = x ^ 2 + (Real.pi - |y|) ^ 2 := by
        rw [← Complex.normSq_eq_norm_sq, Complex.normSq_apply, abs_of_neg hy]
        simp [ip, x, y]; ring
      nlinarith
  have hpi2 : Real.pi ^ 2 ≥ 4 := by nlinarith [Real.pi_gt_three]
  have hc2 : ((Real.pi - |y|) / Real.pi) ^ 2 ≤ Real.cos (y / 2) ^ 2 :=
    pow_le_pow_left₀ hcos0 hcos 2
  have key : (t / Real.pi) ^ 2 ≤ ‖Complex.cosh (z / 2)‖ ^ 2 := by
    rw [hsq]
    have e1 : (t / Real.pi) ^ 2 = t ^ 2 / Real.pi ^ 2 := by ring
    have e2 : ((Real.pi - |y|) / Real.pi) ^ 2 = (Real.pi - |y|) ^ 2 / Real.pi ^ 2 := by ring
    rw [e1]; rw [e2] at hc2
    have h3 : x ^ 2 / Real.pi ^ 2 ≤ (x / 2) ^ 2 := by
      rw [div_le_iff₀ (by positivity)]; nlinarith [sq_nonneg x]
    have h4 : t ^ 2 / Real.pi ^ 2 ≤ x ^ 2 / Real.pi ^ 2 + (Real.pi - |y|) ^ 2 / Real.pi ^ 2 := by
      rw [← add_div]; exact div_le_div_of_nonneg_right ht2 (by positivity)
    linarith
  have key2 : t / Real.pi ≤ ‖Complex.cosh (z / 2)‖ :=
    (pow_le_pow_iff_left₀ (by positivity) (norm_nonneg _) (by norm_num)).mp key
  rw [norm_mul]; norm_num
  have : 2 * t / Real.pi = 2 * (t / Real.pi) := by ring
  linarith

lemma abs_log_le_of {x lo hi : ℝ} (hlo : 0 < lo) (h1 : lo ≤ x) (h2 : x ≤ hi) :
    |Real.log x| ≤ |Real.log lo| + |Real.log hi| := by
  have a1 := Real.log_le_log hlo h1
  have a2 := Real.log_le_log (by linarith) h2
  rw [abs_le]; constructor
  · have := neg_abs_le (Real.log lo); linarith [abs_nonneg (Real.log hi)]
  · have := le_abs_self (Real.log hi); linarith [abs_nonneg (Real.log lo)]

lemma norm_log_le_of {W : ℂ} {lo hi : ℝ} (hlo : 0 < lo) (h1 : lo ≤ ‖W‖) (h2 : ‖W‖ ≤ hi) :
    ‖Complex.log W‖ ≤ |Real.log lo| + |Real.log hi| + Real.pi := by
  refine (Complex.norm_le_abs_re_add_abs_im _).trans ?_
  rw [Complex.log_re, Complex.log_im]
  have := abs_log_le_of hlo h1 h2
  have := Complex.abs_arg_le_pi W
  linarith

lemma abs_log_t {t : ℝ} (ht : 0 < t) (ht2 : t ≤ 2) : |Real.log t| ≤ Real.log 2 + Real.log (2 / t) := by
  have hl2 : Real.log (2 / t) = Real.log 2 - Real.log t := Real.log_div (by norm_num) ht.ne'
  have h0 : 0 ≤ Real.log (2 / t) := Real.log_nonneg (by rw [le_div_iff₀ ht]; linarith)
  have hl2p : 0 < Real.log 2 := Real.log_pos (by norm_num)
  rcases le_or_gt t 1 with h | h
  · rw [abs_of_nonpos (Real.log_nonpos ht.le h)]; linarith
  · rw [abs_of_pos (Real.log_pos h)]
    have := Real.log_le_log ht ht2
    linarith

lemma log_consts : Real.log 2 < 0.7 ∧ Real.log Real.pi < 1.4 ∧ Real.log 16 < 2.8 ∧
    Real.log (2 * Real.pi) < 2.1 := by
  have h2 := Real.log_two_lt_d9
  have hpi := Real.pi_lt_d2
  have hpi0 := Real.pi_pos
  have hpi4 : Real.log Real.pi < Real.log 4 := Real.log_lt_log hpi0 (by linarith)
  have h4 : Real.log 4 = 2 * Real.log 2 := by
    rw [show (4 : ℝ) = 2 ^ 2 by norm_num, Real.log_pow]; norm_num
  have h16 : Real.log 16 = 4 * Real.log 2 := by
    rw [show (16 : ℝ) = 2 ^ 4 by norm_num, Real.log_pow]; norm_num
  have h8 : Real.log (2 * Real.pi) < Real.log 8 := Real.log_lt_log (by positivity) (by linarith)
  have h8' : Real.log 8 = 3 * Real.log 2 := by
    rw [show (8 : ℝ) = 2 ^ 3 by norm_num, Real.log_pow]; norm_num
  refine ⟨by linarith, by linarith, by linarith, by linarith⟩

/-- `|λ(z)| ≤ 11 + log(2/t)` -/
lemma lc_bound {z : ℂ} (hz : ‖z‖ < Real.pi) {t : ℝ} (ht : 0 < t) (ht2 : t ≤ 2)
    (ha : t ≤ ‖z - ip‖) (hb : t ≤ ‖z + ip‖) : ‖lc z‖ ≤ 11 + Real.log (2 / t) := by
  have hpi := Real.pi_pos
  have hlo := norm_cosh_half_ge hz ht ha hb
  have hhi := norm_cosh_half_le hz
  have hW := norm_log_le_of (by positivity) hlo hhi
  have hl : Real.log (2 * t / Real.pi) = Real.log 2 + Real.log t - Real.log Real.pi := by
    rw [Real.log_div (by positivity) hpi.ne', Real.log_mul (by norm_num) ht.ne']
  have habs : |Real.log (2 * t / Real.pi)| ≤ Real.log 2 + |Real.log t| + Real.log Real.pi := by
    rw [hl]
    have h2 : 0 ≤ Real.log 2 := Real.log_nonneg (by norm_num)
    have h3 : 0 ≤ Real.log Real.pi := Real.log_nonneg (by linarith [Real.pi_gt_three])
    calc |Real.log 2 + Real.log t - Real.log Real.pi| ≤ |Real.log 2 + Real.log t| + |Real.log Real.pi| :=
          abs_sub _ _
      _ ≤ |Real.log 2| + |Real.log t| + |Real.log Real.pi| := by gcongr; exact abs_add_le _ _
      _ = _ := by rw [abs_of_nonneg h2, abs_of_nonneg h3]
  have ht' := abs_log_t ht ht2
  obtain ⟨c1, c2, c3, _⟩ := log_consts
  have h16 : |Real.log 16| = Real.log 16 := abs_of_nonneg (Real.log_nonneg (by norm_num))
  have hz2 : ‖z / 2‖ ≤ Real.pi / 2 := by rw [norm_div]; norm_num; linarith
  have hpi' := Real.pi_lt_d2
  calc ‖lc z‖ ≤ ‖z / 2‖ + ‖Complex.log (2 * Complex.cosh (z / 2))‖ := norm_add_le _ _
    _ ≤ Real.pi / 2 + (|Real.log (2 * t / Real.pi)| + |Real.log 16| + Real.pi) := by gcongr
    _ ≤ 11 + Real.log (2 / t) := by rw [h16]; linarith

/-- `|Log(c - z)| ≤ 6 + log(2/t)` -/
lemma Lc_bound {c : ℂ} (hc : c = ip ∨ c = -ip) {z : ℂ} (hz : ‖z‖ < Real.pi) {t : ℝ} (ht : 0 < t)
    (ht2 : t ≤ 2) (hd : t ≤ ‖z - c‖) : ‖Complex.log (c - z)‖ ≤ 6 + Real.log (2 / t) := by
  have hpi := Real.pi_pos
  have h1 : t ≤ ‖c - z‖ := by rwa [norm_sub_rev]
  have h2 : ‖c - z‖ ≤ 2 * Real.pi := by
    have := norm_sub_le c z; rw [norm_c hc] at this; linarith
  have hW := norm_log_le_of ht h1 h2
  have ht' := abs_log_t ht ht2
  obtain ⟨c1, _, _, c4⟩ := log_consts
  have h2p : |Real.log (2 * Real.pi)| = Real.log (2 * Real.pi) :=
    abs_of_nonneg (Real.log_nonneg (by linarith [Real.pi_gt_three]))
  have hpi' := Real.pi_lt_d2
  rw [h2p] at hW
  linarith

/-- `|λ'(z)| ≤ 13/t` -/
lemma lcd_bound {z : ℂ} (hz : ‖z‖ < Real.pi) {t : ℝ} (ht : 0 < t)
    (ha : t ≤ ‖z - ip‖) (hb : t ≤ ‖z + ip‖) : ‖lcd z‖ ≤ 13 / t := by
  have hpi := Real.pi_pos
  have hlo := norm_cosh_half_ge hz ht ha hb
  have he := norm_exp_half_le hz
  have hq : Complex.exp (z / 2) ≠ 0 := Complex.exp_ne_zero _
  have hc0 : 2 * Complex.cosh (z / 2) ≠ 0 :=
    Complex.slitPlane_ne_zero (two_cosh_slit (abs_im_lt hz))
  have e : lcd z = Complex.exp (z / 2) / (2 * Complex.cosh (z / 2)) := by
    rw [lcd, one_add_exp_eq]
    have e1 : Complex.exp z = Complex.exp (z / 2) * Complex.exp (z / 2) := by
      rw [← Complex.exp_add]; ring_nf
    rw [e1]; field_simp
  rw [e, norm_div]
  have hpos : 0 < 2 * t / Real.pi := by positivity
  calc ‖Complex.exp (z / 2)‖ / ‖2 * Complex.cosh (z / 2)‖ ≤ 8 / (2 * t / Real.pi) :=
        div_le_div₀ (by norm_num) he hpos hlo
    _ = 4 * Real.pi / t := by field_simp; ring
    _ ≤ 13 / t := by
        apply div_le_div_of_nonneg_right _ ht.le
        linarith [Real.pi_lt_d2]

lemma E0_bounds : ∃ K : ℝ, 0 ≤ K ∧ ∀ z : ℂ, ‖z‖ < Real.pi →
    ‖E 0 z‖ ≤ K ∧ HasDerivAt (E 0) (deriv Q02.g2 z) z ∧ ‖deriv Q02.g2 z‖ ≤ K := by
  have hpi := Real.pi_pos
  have han : AnalyticOnNhd ℂ Q02.g2 (Metric.ball 0 (3 * Real.pi)) :=
    Q02.g2_diffOn.analyticOnNhd Metric.isOpen_ball
  have hsub : Metric.closedBall (0 : ℂ) Real.pi ⊆ Metric.ball 0 (3 * Real.pi) :=
    Metric.closedBall_subset_ball (by linarith)
  have hc1 : ContinuousOn Q02.g2 (Metric.closedBall 0 Real.pi) :=
    (Q02.g2_diffOn.continuousOn).mono hsub
  have hc2 : ContinuousOn (deriv Q02.g2) (Metric.closedBall 0 Real.pi) :=
    (han.deriv.continuousOn).mono hsub
  obtain ⟨K1, hK1⟩ := (isCompact_closedBall (0 : ℂ) Real.pi).exists_bound_of_continuousOn hc1
  obtain ⟨K2, hK2⟩ := (isCompact_closedBall (0 : ℂ) Real.pi).exists_bound_of_continuousOn hc2
  refine ⟨max (max K1 K2) 0, le_max_right _ _, fun z hz => ?_⟩
  have hzb : z ∈ Metric.closedBall (0 : ℂ) Real.pi := by simpa using hz.le
  have hev : E 0 =ᶠ[𝓝 z] Q02.g2 := by
    filter_upwards [Metric.isOpen_ball.mem_nhds (show z ∈ Metric.ball (0 : ℂ) Real.pi by simpa using hz)]
      with w hw
    exact E0_eq_g2 (by simpa using hw)
  refine ⟨?_, ?_, ?_⟩
  · rw [hev.eq_of_nhds]
    exact (hK1 z hzb).trans ((le_max_left _ _).trans (le_max_left _ _))
  · have hd : DifferentiableAt ℂ Q02.g2 z := (han z (hsub hzb)).differentiableAt
    exact hd.hasDerivAt.congr_of_eventuallyEq hev
  · exact (hK2 z hzb).trans ((le_max_right _ _).trans (le_max_left _ _))

lemma norm_pow_sub_pow_le' (x y : ℂ) (M : ℝ) (hx : ‖x‖ ≤ M) (hy : ‖y‖ ≤ M) (i : ℕ) :
    ‖x ^ i - y ^ i‖ ≤ i * ‖x - y‖ * M ^ (i - 1) := by
  rw [← geom_sum₂_mul, norm_mul, mul_comm]
  have hM : 0 ≤ M := (norm_nonneg x).trans hx
  have h : ‖∑ j ∈ Finset.range i, x ^ j * y ^ (i - 1 - j)‖ ≤ i * M ^ (i - 1) := by
    refine (norm_sum_le _ _).trans ?_
    have : ∀ j ∈ Finset.range i, ‖x ^ j * y ^ (i - 1 - j)‖ ≤ M ^ (i - 1) := by
      intro j hj
      simp only [Finset.mem_range] at hj
      rw [norm_mul, norm_pow, norm_pow]
      calc ‖x‖ ^ j * ‖y‖ ^ (i - 1 - j) ≤ M ^ j * M ^ (i - 1 - j) := by gcongr
        _ = M ^ (i - 1) := by rw [← pow_add]; congr 1; omega
    refine (Finset.sum_le_sum this).trans ?_
    simp
  calc ‖x - y‖ * ‖∑ j ∈ Finset.range i, x ^ j * y ^ (i - 1 - j)‖ ≤ ‖x - y‖ * (i * M ^ (i - 1)) :=
        mul_le_mul_of_nonneg_left h (norm_nonneg _)
    _ = _ := by ring

/-! ### Part 5c: pointwise bounds for `E_i` -/

lemma E_decomp (i : ℕ) (z : ℂ) : E i z = E 0 z * (-lc z) ^ i +
    (ip - z)⁻¹ * ((-lc z) ^ i - (-Complex.log (ip - z)) ^ i) +
    (-ip - z)⁻¹ * ((-lc z) ^ i - (-Complex.log (-ip - z)) ^ i) := by
  simp only [E, bpart, gpart, pow_zero, mul_one, Bc]
  ring

lemma Lam_ge {t : ℝ} (ht : 0 < t) (ht2 : t ≤ 2) : 11 ≤ 11 + Real.log (2 / t) := by
  have : 0 ≤ Real.log (2 / t) := Real.log_nonneg (by rw [le_div_iff₀ ht]; linarith)
  linarith

lemma Sterm_bound {c : ℂ} (hc : c = ip ∨ c = -ip) {z : ℂ} (hz : ‖z‖ < Real.pi) {t : ℝ}
    (ht : 0 < t) (ht2 : t ≤ 2) (ha : t ≤ ‖z - ip‖) (hb : t ≤ ‖z + ip‖) (hd : t ≤ ‖z - c‖)
    (i : ℕ) :
    ‖(c - z)⁻¹ * ((-lc z) ^ i - (-Complex.log (c - z)) ^ i)‖ ≤
      ((i : ℝ) + 4) * (11 + Real.log (2 / t)) ^ i := by
  set Λ := 11 + Real.log (2 / t) with hΛ
  have hΛ1 : 11 ≤ Λ := Lam_ge ht ht2
  have h1 : ‖-lc z‖ ≤ Λ := by rw [norm_neg]; exact lc_bound hz ht ht2 ha hb
  have h2 : ‖-Complex.log (c - z)‖ ≤ Λ := by
    rw [norm_neg]; have := Lc_bound hc hz ht ht2 hd; linarith
  have hΛp : 0 ≤ Λ ^ i := by positivity
  have hpos : 0 < ‖z - c‖ := ht.trans_le hd
  rw [norm_mul, norm_inv, norm_sub_rev c z]
  rcases lt_or_ge ‖z - c‖ (1 / 2) with hn | hn
  · have hdl := (d_local hc hz hn).1
    have hdiff : ‖-lc z - -Complex.log (c - z)‖ = ‖dfun c z‖ := by
      rw [dfun, ← norm_neg]; congr 1; ring
    have hp := norm_pow_sub_pow_le' _ _ Λ h1 h2 i
    rw [hdiff] at hp
    have hp2 : ‖(-lc z) ^ i - (-Complex.log (c - z)) ^ i‖ ≤ i * ‖z - c‖ * Λ ^ (i - 1) :=
      hp.trans (by gcongr)
    have hΛi : Λ ^ (i - 1) ≤ Λ ^ i := pow_le_pow_right₀ (by linarith) (Nat.sub_le i 1)
    calc ‖z - c‖⁻¹ * ‖(-lc z) ^ i - (-Complex.log (c - z)) ^ i‖
        ≤ ‖z - c‖⁻¹ * (i * ‖z - c‖ * Λ ^ (i - 1)) := by gcongr
      _ = i * Λ ^ (i - 1) := by field_simp
      _ ≤ ((i : ℝ) + 4) * Λ ^ i := by
          have : (i : ℝ) * Λ ^ (i - 1) ≤ i * Λ ^ i := by gcongr
          nlinarith
  · have hinv : ‖z - c‖⁻¹ ≤ 2 := by
      rw [inv_le_comm₀ hpos (by norm_num)]; linarith
    have hp : ‖(-lc z) ^ i - (-Complex.log (c - z)) ^ i‖ ≤ 2 * Λ ^ i := by
      refine (norm_sub_le _ _).trans ?_
      rw [norm_pow, norm_pow]
      have := pow_le_pow_left₀ (norm_nonneg _) h1 i
      have := pow_le_pow_left₀ (norm_nonneg _) h2 i
      linarith
    calc ‖z - c‖⁻¹ * ‖(-lc z) ^ i - (-Complex.log (c - z)) ^ i‖ ≤ 2 * (2 * Λ ^ i) := by
          gcongr
      _ ≤ ((i : ℝ) + 4) * Λ ^ i := by
          have : (0 : ℝ) ≤ i := Nat.cast_nonneg i
          nlinarith

lemma E_pointwise : ∃ K : ℝ, 0 ≤ K ∧ ∀ (i : ℕ) (z : ℂ), ‖z‖ < Real.pi → ∀ t : ℝ, 0 < t → t ≤ 2 →
    t ≤ ‖z - ip‖ → t ≤ ‖z + ip‖ →
    ‖E i z‖ ≤ K * (((i : ℝ) + 1) * (11 + Real.log (2 / t)) ^ i) := by
  obtain ⟨K4, hK4, hE0⟩ := E0_bounds
  refine ⟨K4 + 8, by linarith, fun i z hz t ht ht2 ha hb => ?_⟩
  set Λ := 11 + Real.log (2 / t) with hΛ
  have hΛ1 : 11 ≤ Λ := Lam_ge ht ht2
  have h1 : ‖-lc z‖ ≤ Λ := by rw [norm_neg]; exact lc_bound hz ht ht2 ha hb
  have hb' : t ≤ ‖z - -ip‖ := by rwa [sub_neg_eq_add]
  have s1 := Sterm_bound (Or.inl rfl) hz ht ht2 ha hb ha i
  have s2 := Sterm_bound (Or.inr rfl) hz ht ht2 ha hb hb' i
  rw [E_decomp]
  have e0 : ‖E 0 z * (-lc z) ^ i‖ ≤ K4 * Λ ^ i := by
    rw [norm_mul, norm_pow]
    exact mul_le_mul (hE0 z hz).1 (pow_le_pow_left₀ (norm_nonneg _) h1 i) (by positivity) hK4
  have hΛp : 0 ≤ Λ ^ i := by positivity
  have hi : (0 : ℝ) ≤ i := Nat.cast_nonneg i
  calc _ ≤ ‖E 0 z * (-lc z) ^ i‖ + ‖(ip - z)⁻¹ * ((-lc z) ^ i - (-Complex.log (ip - z)) ^ i)‖ +
        ‖(-ip - z)⁻¹ * ((-lc z) ^ i - (-Complex.log (-ip - z)) ^ i)‖ := norm_add₃_le
    _ ≤ K4 * Λ ^ i + ((i : ℝ) + 4) * Λ ^ i + ((i : ℝ) + 4) * Λ ^ i := by gcongr
    _ ≤ (K4 + 8) * (((i : ℝ) + 1) * Λ ^ i) := by nlinarith [mul_nonneg hK4 hΛp, mul_nonneg hi hΛp]

/-! ### Part 5d: Cauchy's estimate on `|z| = π(m+1)/(m+2)` -/

lemma pow_ratio_le (m : ℕ) :
    Real.pi ^ m ≤ 3 * (Real.pi * ((m : ℝ) + 1) / ((m : ℝ) + 2)) ^ m := by
  have hpi := Real.pi_pos
  have hm : (0 : ℝ) ≤ m := Nat.cast_nonneg m
  have e : Real.pi ^ m = (Real.pi * ((m : ℝ) + 1) / ((m : ℝ) + 2)) ^ m *
      (((m : ℝ) + 2) / ((m : ℝ) + 1)) ^ m := by
    rw [← mul_pow]; congr 1; field_simp
  rw [e]
  have hq : (((m : ℝ) + 2) / ((m : ℝ) + 1)) ^ m ≤ 3 := by
    have h1 : ((m : ℝ) + 2) / ((m : ℝ) + 1) ≤ Real.exp (1 / ((m : ℝ) + 1)) := by
      have := Real.add_one_le_exp (1 / ((m : ℝ) + 1))
      have e2 : ((m : ℝ) + 2) / ((m : ℝ) + 1) = 1 / ((m : ℝ) + 1) + 1 := by field_simp; ring
      linarith
    have h2 := pow_le_pow_left₀ (by positivity) h1 m
    rw [← Real.exp_nat_mul] at h2
    have h3 : (m : ℝ) * (1 / ((m : ℝ) + 1)) ≤ 1 := by
      rw [mul_one_div, div_le_one (by positivity)]; linarith
    have h4 := Real.exp_le_exp.mpr h3
    have h5 := Real.exp_one_lt_d9
    linarith
  have : 0 ≤ (Real.pi * ((m : ℝ) + 1) / ((m : ℝ) + 2)) ^ m := by positivity
  nlinarith

lemma cauchy_circle {f : ℂ → ℂ} (hf : AnalyticOnNhd ℂ f (Metric.ball 0 Real.pi)) (m : ℕ) {M : ℝ}
    (hM0 : 0 ≤ M)
    (hM : ∀ z : ℂ, ‖z‖ = Real.pi * ((m : ℝ) + 1) / ((m : ℝ) + 2) → ‖f z‖ ≤ M) :
    ‖iteratedDeriv m f 0‖ ≤ 3 * (m.factorial : ℝ) * M / Real.pi ^ m := by
  have hpi := Real.pi_pos
  have hm : (0 : ℝ) ≤ m := Nat.cast_nonneg m
  set r := Real.pi * ((m : ℝ) + 1) / ((m : ℝ) + 2) with hr
  have hr0 : 0 < r := by positivity
  have hrpi : r < Real.pi := by
    rw [hr, div_lt_iff₀ (by positivity)]; nlinarith
  have hsub : Metric.closedBall (0 : ℂ) r ⊆ Metric.ball 0 Real.pi :=
    Metric.closedBall_subset_ball hrpi
  have hdc : DiffContOnCl ℂ f (Metric.ball 0 r) :=
    (hf.differentiableOn).diffContOnCl_ball hsub
  have hC : ∀ z ∈ Metric.sphere (0 : ℂ) r, ‖f z‖ ≤ M := by
    intro z hz; exact hM z (by simpa using hz)
  have h := Complex.norm_iteratedDeriv_le_of_forall_mem_sphere_norm_le m hr0 hdc hC
  refine h.trans ?_
  have hpr := pow_ratio_le m
  rw [← hr] at hpr
  have hF : (0 : ℝ) ≤ (m.factorial : ℝ) * M := by positivity
  rw [div_le_div_iff₀ (by positivity) (by positivity)]
  calc (m.factorial : ℝ) * M * Real.pi ^ m ≤ (m.factorial : ℝ) * M * (3 * r ^ m) :=
        mul_le_mul_of_nonneg_left hpr hF
    _ = 3 * (m.factorial : ℝ) * M * r ^ m := by ring

lemma sphere_dist {m : ℕ} {z : ℂ} (hz : ‖z‖ = Real.pi * ((m : ℝ) + 1) / ((m : ℝ) + 2)) :
    ‖z‖ < Real.pi ∧ Real.pi / ((m : ℝ) + 2) ≤ ‖z - ip‖ ∧ Real.pi / ((m : ℝ) + 2) ≤ ‖z + ip‖ := by
  have hpi := Real.pi_pos
  have hm : (0 : ℝ) ≤ m := Nat.cast_nonneg m
  have e : Real.pi / ((m : ℝ) + 2) = Real.pi - ‖z‖ := by rw [hz]; field_simp; ring
  refine ⟨?_, ?_, ?_⟩
  · rw [hz, div_lt_iff₀ (by positivity)]; nlinarith
  · rw [e]; have := norm_sub_norm_le ip z; rw [norm_ip, norm_sub_rev] at this; linarith
  · rw [e]
    have := norm_sub_norm_le (-ip) z
    rw [norm_neg, norm_ip, norm_sub_rev, sub_neg_eq_add] at this; linarith

lemma log_t_le (m : ℕ) : Real.log (2 / (Real.pi / ((m : ℝ) + 2))) ≤ Real.log ((m : ℝ) + 2) := by
  have hpi := Real.pi_pos
  have hm : (0 : ℝ) ≤ m := Nat.cast_nonneg m
  apply Real.log_le_log (by positivity)
  rw [div_div_eq_mul_div, div_le_iff₀ hpi]
  nlinarith [Real.pi_gt_three]

lemma t_le_two (m : ℕ) : Real.pi / ((m : ℝ) + 2) ≤ 2 := by
  have hm : (0 : ℝ) ≤ m := Nat.cast_nonneg m
  rw [div_le_iff₀ (by positivity)]; nlinarith [Real.pi_lt_d2]

lemma E_analytic (i : ℕ) : AnalyticOnNhd ℂ (E i) (Metric.ball 0 Real.pi) := by
  intro z hz
  have hz' : ‖z‖ < Real.pi := by simpa using hz
  unfold E
  exact ((bpart_analyticAt hz' i).sub (gpart_analyticAt ip z (c_sub_slit (Or.inl rfl) hz') i)).sub
    (gpart_analyticAt (-ip) z (c_sub_slit (Or.inr rfl) hz') i)

/-! ### Part 5e: the derivatives of `E_1`, `E_2` -/

/-- `(S_c d_c)' = S_c^2 d_c + S_c (λ' + S_c)` -/
noncomputable def Aterm (c z : ℂ) : ℂ :=
  ((c - z)⁻¹) ^ 2 * dfun c z + (c - z)⁻¹ * (lcd z + (c - z)⁻¹)

lemma S_hasDerivAt (c z : ℂ) (h : c - z ≠ 0) :
    HasDerivAt (fun w : ℂ => (c - w)⁻¹) (((c - z)⁻¹) ^ 2) z := by
  have h1 : HasDerivAt (fun w : ℂ => c - w) (-1) z := by
    simpa using (hasDerivAt_id z).const_sub c
  have h2 := h1.inv h
  convert h2 using 1
  field_simp

lemma Sd_hasDerivAt {c : ℂ} (hc : c = ip ∨ c = -ip) {z : ℂ} (hz : ‖z‖ < Real.pi) :
    HasDerivAt (fun w => (c - w)⁻¹ * dfun c w) (Aterm c z) z := by
  have hne : c - z ≠ 0 := Complex.slitPlane_ne_zero (c_sub_slit hc hz)
  have := (S_hasDerivAt c z hne).mul (dfun_hasDerivAt hc hz)
  convert this using 1
  rw [Aterm]

lemma Aterm_bound {c : ℂ} (hc : c = ip ∨ c = -ip) {z : ℂ} (hz : ‖z‖ < Real.pi) {t : ℝ}
    (ht : 0 < t) (ht2 : t ≤ 2) (ha : t ≤ ‖z - ip‖) (hb : t ≤ ‖z + ip‖) (hd : t ≤ ‖z - c‖) :
    ‖(c - z)⁻¹ * dfun c z‖ ≤ 4 * (11 + Real.log (2 / t)) ∧
      ‖Aterm c z‖ ≤ 8 * (11 + Real.log (2 / t)) + 34 / t ∧
      (‖z - c‖ < 1 / 2 → ‖Aterm c z‖ ≤ 1) := by
  set Λ := 11 + Real.log (2 / t) with hΛ
  have hΛ1 : 11 ≤ Λ := Lam_ge ht ht2
  have h1 : ‖lc z‖ ≤ Λ := lc_bound hz ht ht2 ha hb
  have h2 : ‖Complex.log (c - z)‖ ≤ Λ := by have := Lc_bound hc hz ht ht2 hd; linarith
  have hdd : ‖dfun c z‖ ≤ 2 * Λ := by
    rw [dfun]; refine (norm_sub_le _ _).trans ?_; linarith
  have hl := lcd_bound hz ht ha hb
  have hpos : 0 < ‖z - c‖ := ht.trans_le hd
  have hu0 : z - c ≠ 0 := norm_pos_iff.mp hpos
  have hS : (c - z)⁻¹ = -(z - c)⁻¹ := by rw [← inv_neg, neg_sub]
  have hnear : ‖z - c‖ < 1 / 2 → ‖(c - z)⁻¹ * dfun c z‖ ≤ 1 ∧ ‖Aterm c z‖ ≤ 1 := by
    intro hn
    obtain ⟨d1, d2⟩ := d_local hc hz hn
    constructor
    · rw [norm_mul, norm_inv, norm_sub_rev c z]
      calc ‖z - c‖⁻¹ * ‖dfun c z‖ ≤ ‖z - c‖⁻¹ * ‖z - c‖ := by gcongr
        _ = 1 := inv_mul_cancel₀ hpos.ne'
    · have e : Aterm c z = (dfun c z - (z - c) * (lcd z + (c - z)⁻¹)) / (z - c) ^ 2 := by
        rw [Aterm, hS]; field_simp; ring
      rw [e, norm_div, norm_pow, div_le_one (by positivity)]
      exact d2
  have hinv2 : ‖z - c‖ < 1 / 2 ∨ ‖(c - z)⁻¹‖ ≤ 2 := by
    rcases lt_or_ge ‖z - c‖ (1 / 2) with hn | hn
    · exact Or.inl hn
    · right; rw [norm_inv, norm_sub_rev, inv_le_comm₀ hpos (by norm_num)]; linarith
  have ht1 : 1 / t ≥ 1 / 2 := by rw [ge_iff_le, div_le_div_iff₀ (by norm_num) ht]; linarith
  refine ⟨?_, ?_, fun hn => (hnear hn).2⟩
  · rcases hinv2 with hn | hS2
    · linarith [(hnear hn).1]
    · rw [norm_mul]
      calc ‖(c - z)⁻¹‖ * ‖dfun c z‖ ≤ 2 * (2 * Λ) := by gcongr
        _ = 4 * Λ := by ring
  · rcases hinv2 with hn | hS2
    · have := (hnear hn).2
      have : 0 ≤ 34 / t := by positivity
      linarith
    · have hS0 : 0 ≤ ‖(c - z)⁻¹‖ := norm_nonneg _
      have e1 : ‖Aterm c z‖ ≤ ‖(c - z)⁻¹‖ ^ 2 * ‖dfun c z‖ + ‖(c - z)⁻¹‖ * (‖lcd z‖ + ‖(c - z)⁻¹‖) := by
        rw [Aterm]
        refine (norm_add_le _ _).trans ?_
        rw [norm_mul, norm_mul, norm_pow]
        gcongr
        exact norm_add_le _ _
      have e2 : ‖(c - z)⁻¹‖ ^ 2 * ‖dfun c z‖ ≤ 4 * (2 * Λ) := by
        have : ‖(c - z)⁻¹‖ ^ 2 ≤ 4 := by nlinarith
        exact mul_le_mul this hdd (norm_nonneg _) (by norm_num)
      have e3 : ‖(c - z)⁻¹‖ * (‖lcd z‖ + ‖(c - z)⁻¹‖) ≤ 2 * (13 / t + 2) := by
        apply mul_le_mul hS2 (by linarith) (by positivity) (by norm_num)
      have e4 : (2 : ℝ) * (13 / t + 2) ≤ 34 / t := by
        have : (4 : ℝ) ≤ 8 / t := by rw [le_div_iff₀ ht]; linarith
        have : 2 * (13 / t + 2) = 26 / t + 4 := by ring
        have : 34 / t = 26 / t + 8 / t := by ring
        linarith
      linarith

lemma E2_eq : E 2 = fun w => E 0 w * (-lc w) ^ 2 -
    (ip - w)⁻¹ * dfun ip w * (-lc w - Complex.log (ip - w)) -
    (-ip - w)⁻¹ * dfun (-ip) w * (-lc w - Complex.log (-ip - w)) := by
  funext w
  simp only [E, bpart, gpart, pow_zero, mul_one, dfun, Bc]
  ring

lemma nuL_hasDerivAt {c : ℂ} (hc : c = ip ∨ c = -ip) {z : ℂ} (hz : ‖z‖ < Real.pi) :
    HasDerivAt (fun w => -lc w - Complex.log (c - w)) (-lcd z + (c - z)⁻¹) z := by
  have h1 := lc_hasDerivAt (abs_im_lt hz)
  have h2 : HasDerivAt (fun z : ℂ => c - z) (-1) z := by
    simpa using (hasDerivAt_id z).const_sub c
  have h3 : HasDerivAt (fun z : ℂ => Complex.log (c - z)) (-1 / (c - z)) z :=
    h2.clog (c_sub_slit hc hz)
  have h4 := h1.neg.sub h3
  convert h4 using 1
  ring

lemma E2_deriv_bound : ∃ C : ℝ, 0 ≤ C ∧ ∀ z : ℂ, ‖z‖ < Real.pi → ∀ t : ℝ, 0 < t → t ≤ 2 →
    t ≤ ‖z - ip‖ → t ≤ ‖z + ip‖ →
    ‖deriv (E 2) z‖ ≤ C * ((11 + Real.log (2 / t)) ^ 2 + (11 + Real.log (2 / t)) / t) := by
  obtain ⟨K4, hK4, hE0⟩ := E0_bounds
  refine ⟨26 * K4 + 280, by linarith, fun z hz t ht ht2 ha hb => ?_⟩
  obtain ⟨hE0n, hE0d, hE0dn⟩ := hE0 z hz
  set e0' := deriv Q02.g2 z
  have hb' : t ≤ ‖z - -ip‖ := by rwa [sub_neg_eq_add]
  have hlc := lc_hasDerivAt (abs_im_lt hz)
  have hP : HasDerivAt (E 2) (e0' * (-lc z) ^ 2 + E 0 z * (2 * (-lc z) * (-lcd z)) -
      (Aterm ip z * (-lc z - Complex.log (ip - z)) + (ip - z)⁻¹ * dfun ip z * (-lcd z + (ip - z)⁻¹)) -
      (Aterm (-ip) z * (-lc z - Complex.log (-ip - z)) +
        (-ip - z)⁻¹ * dfun (-ip) z * (-lcd z + (-ip - z)⁻¹))) z := by
    rw [E2_eq]
    have h1 : HasDerivAt (fun w => (-lc w) ^ 2) (2 * (-lc z) * (-lcd z)) z := by
      have h1' := (hlc.neg).pow 2
      convert h1' using 1
      push_cast; simp only [Pi.neg_apply]; ring
    have h2 := hE0d.mul h1
    have h3 := (Sd_hasDerivAt (Or.inl rfl) hz).mul (nuL_hasDerivAt (Or.inl rfl) hz)
    have h4 := (Sd_hasDerivAt (Or.inr rfl) hz).mul (nuL_hasDerivAt (Or.inr rfl) hz)
    exact (h2.sub h3).sub h4
  rw [hP.deriv]
  set Λ := 11 + Real.log (2 / t) with hΛ
  have hΛ1 : 11 ≤ Λ := Lam_ge ht ht2
  have hν : ‖-lc z‖ ≤ Λ := by rw [norm_neg]; exact lc_bound hz ht ht2 ha hb
  have hν' : ‖lc z‖ ≤ Λ := lc_bound hz ht ht2 ha hb
  have hl := lcd_bound hz ht ha hb
  obtain ⟨p1, p2, -⟩ := Aterm_bound (Or.inl rfl) hz ht ht2 ha hb ha
  obtain ⟨q1, q2, -⟩ := Aterm_bound (Or.inr rfl) hz ht ht2 ha hb hb'
  have hL1 : ‖-lc z - Complex.log (ip - z)‖ ≤ 2 * Λ := by
    have := Lc_bound (Or.inl rfl) hz ht ht2 ha
    refine (norm_sub_le _ _).trans ?_; rw [norm_neg]; linarith
  have hL2 : ‖-lc z - Complex.log (-ip - z)‖ ≤ 2 * Λ := by
    have := Lc_bound (Or.inr rfl) hz ht ht2 hb'
    refine (norm_sub_le _ _).trans ?_; rw [norm_neg]; linarith
  have hS1 : ‖(ip - z)⁻¹‖ ≤ 1 / t := by
    rw [norm_inv, norm_sub_rev, one_div]; exact inv_anti₀ ht ha
  have hS2 : ‖(-ip - z)⁻¹‖ ≤ 1 / t := by
    rw [norm_inv, norm_sub_rev, one_div]; exact inv_anti₀ ht hb'
  have hD1 : ‖-lcd z + (ip - z)⁻¹‖ ≤ 14 / t := by
    refine (norm_add_le _ _).trans ?_; rw [norm_neg]
    have : 13 / t + 1 / t = 14 / t := by ring
    linarith
  have hD2 : ‖-lcd z + (-ip - z)⁻¹‖ ≤ 14 / t := by
    refine (norm_add_le _ _).trans ?_; rw [norm_neg]
    have : 13 / t + 1 / t = 14 / t := by ring
    linarith
  have hΛt : 0 ≤ Λ / t := by positivity
  have t1 : ‖e0' * (-lc z) ^ 2‖ ≤ K4 * Λ ^ 2 := by
    rw [norm_mul, norm_pow]
    exact mul_le_mul hE0dn (pow_le_pow_left₀ (norm_nonneg _) hν 2) (by positivity) hK4
  have t2 : ‖E 0 z * (2 * (-lc z) * (-lcd z))‖ ≤ K4 * (2 * Λ * (13 / t)) := by
    rw [norm_mul, norm_mul, norm_mul, norm_neg]
    norm_num
    apply mul_le_mul hE0n _ (by positivity) hK4
    have : ‖lc z‖ ≤ Λ := by rw [← norm_neg]; exact hν
    exact mul_le_mul (by linarith) hl (norm_nonneg _) (by positivity)
  have t3 : ‖Aterm ip z * (-lc z - Complex.log (ip - z)) +
      (ip - z)⁻¹ * dfun ip z * (-lcd z + (ip - z)⁻¹)‖ ≤
      (8 * Λ + 34 / t) * (2 * Λ) + 4 * Λ * (14 / t) := by
    refine (norm_add_le _ _).trans ?_
    rw [norm_mul, norm_mul (_ * _)]
    gcongr
  have t4 : ‖Aterm (-ip) z * (-lc z - Complex.log (-ip - z)) +
      (-ip - z)⁻¹ * dfun (-ip) z * (-lcd z + (-ip - z)⁻¹)‖ ≤
      (8 * Λ + 34 / t) * (2 * Λ) + 4 * Λ * (14 / t) := by
    refine (norm_add_le _ _).trans ?_
    rw [norm_mul, norm_mul (_ * _)]
    gcongr
  have hsum := norm_sub_le (e0' * (-lc z) ^ 2 + E 0 z * (2 * (-lc z) * (-lcd z)) -
      (Aterm ip z * (-lc z - Complex.log (ip - z)) + (ip - z)⁻¹ * dfun ip z * (-lcd z + (ip - z)⁻¹)))
    (Aterm (-ip) z * (-lc z - Complex.log (-ip - z)) +
        (-ip - z)⁻¹ * dfun (-ip) z * (-lcd z + (-ip - z)⁻¹))
  have hsum2 := norm_sub_le (e0' * (-lc z) ^ 2 + E 0 z * (2 * (-lc z) * (-lcd z)))
    (Aterm ip z * (-lc z - Complex.log (ip - z)) + (ip - z)⁻¹ * dfun ip z * (-lcd z + (ip - z)⁻¹))
  have hsum3 := norm_add_le (e0' * (-lc z) ^ 2) (E 0 z * (2 * (-lc z) * (-lcd z)))
  have e1 : (8 * Λ + 34 / t) * (2 * Λ) + 4 * Λ * (14 / t) = 16 * Λ ^ 2 + 124 * (Λ / t) := by ring
  have e2 : K4 * (2 * Λ * (13 / t)) = 26 * K4 * (Λ / t) := by ring
  have hK : 0 ≤ K4 * Λ ^ 2 := by positivity
  have hK' : 0 ≤ K4 * (Λ / t) := by positivity
  nlinarith

lemma E1_eq : E 1 = fun w => E 0 w * (-lc w) - (ip - w)⁻¹ * dfun ip w -
    (-ip - w)⁻¹ * dfun (-ip) w := by
  funext w
  simp only [E, bpart, gpart, pow_zero, pow_one, mul_one, dfun, Bc]
  ring

lemma E0_eq_c {c : ℂ} (hc : c = ip ∨ c = -ip) (z : ℂ) :
    E 0 z = Bc z - (c - z)⁻¹ - (-c - z)⁻¹ := by
  rcases hc with rfl | rfl
  · simp only [E, bpart, gpart, pow_zero, mul_one]
  · simp only [E, bpart, gpart, pow_zero, mul_one, neg_neg]; ring

/-- the pole-corrected derivative of `E_1` near `c` -/
lemma Q_near {c : ℂ} (hc : c = ip ∨ c = -ip) {z : ℂ} (hz : ‖z‖ < Real.pi) (hn : ‖z - c‖ < 1 / 2)
    {t : ℝ} (ht : 0 < t) (ht2 : t ≤ 2) (ha : t ≤ ‖z - ip‖) (hb : t ≤ ‖z + ip‖)
    {e0' : ℂ} {K4 : ℝ} (hK4 : 0 ≤ K4) (he0 : ‖e0'‖ ≤ K4) :
    ‖e0' * (-lc z) + E 0 z * (-lcd z) - Aterm c z - Aterm (-c) z +
        3 / 2 * ((c - z)⁻¹ + (-c - z)⁻¹)‖ ≤ (K4 + 5) * (11 + Real.log (2 / t)) := by
  set Λ := 11 + Real.log (2 / t) with hΛ
  have hΛ1 : 11 ≤ Λ := Lam_ge ht ht2
  have hpi := Real.pi_pos
  have hcn := norm_c hc
  have hc' : -c = ip ∨ -c = -ip := by rcases hc with rfl | rfl; right; rfl; left; simp
  have hzc : z ≠ c := by intro h; rw [h, hcn] at hz; linarith
  have hu0 : z - c ≠ 0 := sub_ne_zero.mpr hzc
  have hpos : 0 < ‖z - c‖ := norm_pos_iff.mpr hu0
  have hdc : t ≤ ‖z - c‖ := by rcases hc with rfl | rfl; exact ha; rwa [sub_neg_eq_add]
  -- the other point is far
  have hfar : 5 ≤ ‖z - -c‖ := by
    have e : z - -c = 2 * c + (z - c) := by ring
    have h1 := norm_sub_norm_le (2 * c) (-(z - c))
    rw [norm_neg, sub_neg_eq_add, ← e, norm_mul, hcn] at h1
    norm_num at h1
    have := Real.pi_gt_three
    rw [sub_neg_eq_add]
    linarith
  have hdc' : t ≤ ‖z - -c‖ := by linarith
  have hS' : ‖(-c - z)⁻¹‖ ≤ 1 / 5 := by
    rw [norm_inv, norm_sub_rev, one_div]; exact inv_anti₀ (by norm_num) hfar
  have hν : ‖-lc z‖ ≤ Λ := by rw [norm_neg]; exact lc_bound hz ht ht2 ha hb
  have hd' : ‖dfun (-c) z‖ ≤ 2 * Λ := by
    rw [dfun]; refine (norm_sub_le _ _).trans ?_
    have := lc_bound hz ht ht2 ha hb
    have := Lc_bound hc' hz ht ht2 hdc'
    linarith
  have hB := Bc_local hc hn.le hzc
  have hL := lcd_local hc hn.le hzc
  have hA := (Aterm_bound hc hz ht ht2 ha hb hdc).2.2 hn
  have hS : (c - z)⁻¹ = -(z - c)⁻¹ := by rw [← inv_neg, neg_sub]
  have key : e0' * (-lc z) + E 0 z * (-lcd z) - Aterm c z - Aterm (-c) z +
        3 / 2 * ((c - z)⁻¹ + (-c - z)⁻¹) =
      e0' * (-lc z) - (Bc z - (c - z)⁻¹ + 3 / 2) * lcd z + 3 / 2 * (lcd z - (z - c)⁻¹) -
        Aterm c z - (((-c - z)⁻¹) ^ 2 * dfun (-c) z + ((-c - z)⁻¹) ^ 2) +
        3 / 2 * (-c - z)⁻¹ := by
    rw [E0_eq_c hc z]
    simp only [Aterm]
    rw [show lcd z - (z - c)⁻¹ = lcd z + (c - z)⁻¹ by rw [hS]; ring]
    ring
  rw [key]
  have t1 : ‖e0' * (-lc z)‖ ≤ K4 * Λ := by
    rw [norm_mul]; exact mul_le_mul he0 hν (norm_nonneg _) hK4
  have t2 : ‖(Bc z - (c - z)⁻¹ + 3 / 2) * lcd z‖ ≤ 33 := by
    rw [norm_mul]
    have hl : ‖lcd z‖ ≤ 3 + ‖z - c‖⁻¹ := by
      have := norm_add_le (lcd z - (z - c)⁻¹) ((z - c)⁻¹)
      rw [sub_add_cancel, norm_inv] at this; linarith
    calc ‖Bc z - (c - z)⁻¹ + 3 / 2‖ * ‖lcd z‖ ≤ 13 * ‖z - c‖ * (3 + ‖z - c‖⁻¹) :=
          mul_le_mul hB hl (norm_nonneg _) (by positivity)
      _ = 39 * ‖z - c‖ + 13 := by field_simp; ring
      _ ≤ 33 := by linarith
  have t3 : ‖3 / 2 * (lcd z - (z - c)⁻¹)‖ ≤ 9 / 2 := by
    rw [norm_mul]; norm_num; linarith
  have t4 : ‖((-c - z)⁻¹) ^ 2 * dfun (-c) z + ((-c - z)⁻¹) ^ 2‖ ≤ 2 * Λ / 25 + 1 / 25 := by
    refine (norm_add_le _ _).trans ?_
    rw [norm_mul, norm_pow]
    have h2 : ‖(-c - z)⁻¹‖ ^ 2 ≤ 1 / 25 := by
      have := pow_le_pow_left₀ (norm_nonneg _) hS' 2; norm_num at this ⊢; linarith
    have : ‖(-c - z)⁻¹‖ ^ 2 * ‖dfun (-c) z‖ ≤ 1 / 25 * (2 * Λ) :=
      mul_le_mul h2 hd' (norm_nonneg _) (by norm_num)
    linarith
  have t5 : ‖3 / 2 * (-c - z)⁻¹‖ ≤ 3 / 10 := by
    have h32 : ‖(3 / 2 : ℂ)‖ = 3 / 2 := by norm_num
    rw [norm_mul, h32]; linarith
  have hsum : ‖e0' * (-lc z) - (Bc z - (c - z)⁻¹ + 3 / 2) * lcd z + 3 / 2 * (lcd z - (z - c)⁻¹) -
        Aterm c z - (((-c - z)⁻¹) ^ 2 * dfun (-c) z + ((-c - z)⁻¹) ^ 2) +
        3 / 2 * (-c - z)⁻¹‖ ≤ K4 * Λ + 33 + 9 / 2 + 1 + (2 * Λ / 25 + 1 / 25) + 3 / 10 := by
    refine (norm_add_le _ _).trans ?_
    gcongr
    refine (norm_sub_le _ _).trans ?_
    gcongr
    refine (norm_sub_le _ _).trans ?_
    gcongr
    refine (norm_add_le _ _).trans ?_
    gcongr
    refine (norm_sub_le _ _).trans ?_
    gcongr
  have : 0 ≤ K4 * Λ := by positivity
  nlinarith

lemma log4_lt : Real.log (2 / (1 / 2)) < 1.4 := by
  have h := Real.log_two_lt_d9
  rw [show (2 : ℝ) / (1 / 2) = 2 ^ 2 by norm_num, Real.log_pow]
  push_cast; linarith

lemma Q_far {z : ℂ} (hz : ‖z‖ < Real.pi) (ha : 1 / 2 ≤ ‖z - ip‖) (hb : 1 / 2 ≤ ‖z + ip‖)
    {e0' : ℂ} {K4 : ℝ} (hK4 : 0 ≤ K4) (he0 : ‖e0'‖ ≤ K4) (hE : ‖E 0 z‖ ≤ K4) :
    ‖e0' * (-lc z) + E 0 z * (-lcd z) - Aterm ip z - Aterm (-ip) z +
        3 / 2 * ((ip - z)⁻¹ + (-ip - z)⁻¹)‖ ≤ 39 * K4 + 350 := by
  have hl4 := log4_lt
  have ht : (0 : ℝ) < 1 / 2 := by norm_num
  have ht2 : (1 / 2 : ℝ) ≤ 2 := by norm_num
  have hb' : 1 / 2 ≤ ‖z - -ip‖ := by rwa [sub_neg_eq_add]
  have hν : ‖-lc z‖ ≤ 13 := by
    rw [norm_neg]; have := lc_bound hz ht ht2 ha hb; linarith
  have hl : ‖lcd z‖ ≤ 26 := by have := lcd_bound hz ht ha hb; norm_num at this; linarith
  have hA1 : ‖Aterm ip z‖ ≤ 172 := by
    have := (Aterm_bound (Or.inl rfl) hz ht ht2 ha hb ha).2.1; norm_num at this ⊢; linarith
  have hA2 : ‖Aterm (-ip) z‖ ≤ 172 := by
    have := (Aterm_bound (Or.inr rfl) hz ht ht2 ha hb hb').2.1; norm_num at this ⊢; linarith
  have hS1 : ‖(ip - z)⁻¹‖ ≤ 2 := by
    rw [norm_inv, norm_sub_rev, inv_le_comm₀ (by linarith) (by norm_num)]; linarith
  have hS2 : ‖(-ip - z)⁻¹‖ ≤ 2 := by
    rw [norm_inv, norm_sub_rev, inv_le_comm₀ (by linarith) (by norm_num)]; linarith
  have t1 : ‖e0' * (-lc z)‖ ≤ K4 * 13 := by
    rw [norm_mul]; exact mul_le_mul he0 hν (norm_nonneg _) hK4
  have t2 : ‖E 0 z * (-lcd z)‖ ≤ K4 * 26 := by
    rw [norm_mul, norm_neg]; exact mul_le_mul hE hl (norm_nonneg _) hK4
  have t3 : ‖3 / 2 * ((ip - z)⁻¹ + (-ip - z)⁻¹)‖ ≤ 6 := by
    have h32 : ‖(3 / 2 : ℂ)‖ = 3 / 2 := by norm_num
    rw [norm_mul, h32]
    have := norm_add_le (ip - z)⁻¹ (-ip - z)⁻¹
    linarith
  calc _ ≤ ‖e0' * (-lc z)‖ + ‖E 0 z * (-lcd z)‖ + ‖Aterm ip z‖ + ‖Aterm (-ip) z‖ +
        ‖3 / 2 * ((ip - z)⁻¹ + (-ip - z)⁻¹)‖ := by
        refine (norm_add_le _ _).trans ?_
        gcongr
        refine (norm_sub_le _ _).trans ?_
        gcongr
        refine (norm_sub_le _ _).trans ?_
        gcongr
        exact norm_add_le _ _
    _ ≤ K4 * 13 + K4 * 26 + 172 + 172 + 6 := by gcongr
    _ = 39 * K4 + 350 := by ring

lemma E1_hasDerivAt {z : ℂ} (hz : ‖z‖ < Real.pi) {e0' : ℂ} (he0 : HasDerivAt (E 0) e0' z) :
    HasDerivAt (E 1) (e0' * (-lc z) + E 0 z * (-lcd z) - Aterm ip z - Aterm (-ip) z) z := by
  rw [E1_eq]
  have h1 := he0.mul (lc_hasDerivAt (abs_im_lt hz)).neg
  have h2 := Sd_hasDerivAt (Or.inl rfl) hz
  have h3 := Sd_hasDerivAt (Or.inr rfl) hz
  exact (h1.sub h2).sub h3

lemma Q_bound : ∃ C : ℝ, 0 ≤ C ∧ ∀ z : ℂ, ‖z‖ < Real.pi → ∀ t : ℝ, 0 < t → t ≤ 2 →
    t ≤ ‖z - ip‖ → t ≤ ‖z + ip‖ →
    ‖deriv (E 1) z + 3 / 2 * ((ip - z)⁻¹ + (-ip - z)⁻¹)‖ ≤ C * (11 + Real.log (2 / t)) := by
  obtain ⟨K4, hK4, hE0⟩ := E0_bounds
  refine ⟨39 * K4 + 350, by linarith, fun z hz t ht ht2 ha hb => ?_⟩
  obtain ⟨hE0n, hE0d, hE0dn⟩ := hE0 z hz
  rw [(E1_hasDerivAt hz hE0d).deriv]
  have hΛ1 := Lam_ge ht ht2
  rcases lt_or_ge ‖z - ip‖ (1 / 2) with h1 | h1
  · have := Q_near (Or.inl rfl) hz h1 ht ht2 ha hb hK4 hE0dn
    nlinarith
  rcases lt_or_ge ‖z - -ip‖ (1 / 2) with h2 | h2
  · have := Q_near (Or.inr rfl) hz h2 ht ht2 ha hb hK4 hE0dn
    rw [neg_neg] at this
    have e : deriv Q02.g2 z * -lc z + E 0 z * -lcd z - Aterm (-ip) z - Aterm ip z +
        3 / 2 * ((-ip - z)⁻¹ + (ip - z)⁻¹) = deriv Q02.g2 z * -lc z + E 0 z * -lcd z -
        Aterm ip z - Aterm (-ip) z + 3 / 2 * ((ip - z)⁻¹ + (-ip - z)⁻¹) := by ring
    rw [e] at this
    nlinarith
  · have h2' : 1 / 2 ≤ ‖z + ip‖ := by rwa [sub_neg_eq_add] at h2
    have := Q_far hz h1 h2' hK4 hE0dn hE0n
    nlinarith

/-! ### Part 5: the coefficient bounds (TODO) -/

theorem bound0 : ∃ C0 : ℝ, ∀ m : ℕ,
    ‖iteratedDeriv m (E 0) 0‖ ≤ C0 * (m.factorial : ℝ) / (2 * Real.pi) ^ m := by
  have hpi := Real.pi_pos
  have hsub : Metric.closedBall (0 : ℂ) (2 * Real.pi) ⊆ Metric.ball 0 (3 * Real.pi) :=
    Metric.closedBall_subset_ball (by linarith)
  have hdc : DiffContOnCl ℂ Q02.g2 (Metric.ball 0 (2 * Real.pi)) :=
    Q02.g2_diffOn.diffContOnCl_ball hsub
  have hcont : ContinuousOn Q02.g2 (Metric.sphere (0 : ℂ) (2 * Real.pi)) :=
    (Q02.g2_diffOn.continuousOn).mono (Metric.sphere_subset_closedBall.trans hsub)
  obtain ⟨C, hC⟩ := (isCompact_sphere (0 : ℂ) (2 * Real.pi)).exists_bound_of_continuousOn hcont
  refine ⟨C, fun m => ?_⟩
  have hev : E 0 =ᶠ[𝓝 0] Q02.g2 := by
    filter_upwards [Metric.isOpen_ball.mem_nhds (Metric.mem_ball_self hpi)] with z hz
    exact E0_eq_g2 (by simpa using hz)
  rw [hev.iteratedDeriv_eq m]
  have := Complex.norm_iteratedDeriv_le_of_forall_mem_sphere_norm_le m (by positivity) hdc hC
  calc _ ≤ m.factorial * C / (2 * Real.pi) ^ m := this
    _ = _ := by ring

lemma Lam_n_bounds (n : ℕ) :
    11 + Real.log ((n : ℝ) + 2) ≤ 12 * Real.log ((n : ℝ) + 3) ∧
    11 + Real.log ((n : ℝ) + 2) ≤ 12 * ((n : ℝ) + 1) ∧ 1 ≤ Real.log ((n : ℝ) + 3) := by
  have hn : (0 : ℝ) ≤ n := Nat.cast_nonneg n
  have h3 : 1 ≤ Real.log ((n : ℝ) + 3) := by
    have := one_le_log (n + 1) (by omega); push_cast at this
    rw [show (n : ℝ) + 1 + 2 = n + 3 by ring] at this; exact this
  have h23 : Real.log ((n : ℝ) + 2) ≤ Real.log ((n : ℝ) + 3) :=
    Real.log_le_log (by positivity) (by linarith)
  have hx : Real.log ((n : ℝ) + 2) ≤ (n : ℝ) + 2 - 1 := Real.log_le_sub_one_of_pos (by positivity)
  refine ⟨by linarith, by linarith, h3⟩

lemma S_iter (c : ℂ) (hc : c ∈ Complex.slitPlane) (n : ℕ) :
    iteratedDeriv n (fun w : ℂ => (c - w)⁻¹) 0 = (c ^ (n + 1))⁻¹ * (n.factorial : ℂ) := by
  have e : (fun w : ℂ => (c - w)⁻¹) = gpart c 0 := by funext w; simp [gpart]
  rw [e, gpart_eq_Qf c hc n 0, iteratedDeriv_zero, Qf]
  simp [ascPochhammer_eval_one]

lemma norm_S_iter {c : ℂ} (hc : c = ip ∨ c = -ip) (n : ℕ) :
    ‖iteratedDeriv n (fun w : ℂ => (c - w)⁻¹) 0‖ = (n.factorial : ℝ) / Real.pi ^ (n + 1) := by
  have hs : c ∈ Complex.slitPlane := by rcases hc with rfl | rfl; exact ip_slit; exact nip_slit
  rw [S_iter c hs, norm_mul, norm_inv, norm_pow, norm_c hc]
  simp [div_eq_inv_mul]

theorem bound1 : ∃ C1 : ℝ, ∀ m : ℕ, 1 ≤ m →
    ‖iteratedDeriv m (E 1) 0‖ ≤
      C1 * (((m - 1).factorial : ℝ) * Real.log ((m : ℝ) + 2)) / Real.pi ^ m := by
  obtain ⟨C, hC0, hC⟩ := Q_bound
  refine ⟨36 * Real.pi * C + 3, fun m hm => ?_⟩
  obtain ⟨n, rfl⟩ : ∃ n, m = n + 1 := ⟨m - 1, by omega⟩
  have hpi := Real.pi_pos
  have hn : (0 : ℝ) ≤ n := Nat.cast_nonneg n
  set Qfun : ℂ → ℂ := fun w => deriv (E 1) w + 3 / 2 * ((ip - w)⁻¹ + (-ip - w)⁻¹) with hQ
  set P : ℂ → ℂ := fun w => 3 / 2 * ((ip - w)⁻¹ + (-ip - w)⁻¹) with hP
  have hd : deriv (E 1) = fun w => Qfun w - P w := by funext w; simp only [hQ, hP]; ring
  -- analyticity
  have hS1 : AnalyticAt ℂ (fun w : ℂ => (ip - w)⁻¹) 0 :=
    (analyticAt_const.sub analyticAt_id).inv (by simpa using (Complex.slitPlane_ne_zero ip_slit))
  have hS2 : AnalyticAt ℂ (fun w : ℂ => (-ip - w)⁻¹) 0 :=
    (analyticAt_const.sub analyticAt_id).inv (by simpa using (Complex.slitPlane_ne_zero nip_slit))
  have hPa : AnalyticAt ℂ P 0 := analyticAt_const.mul (hS1.add hS2)
  have hQan : AnalyticOnNhd ℂ Qfun (Metric.ball 0 Real.pi) := by
    intro w hw
    have hw' : ‖w‖ < Real.pi := by simpa using hw
    have h1 := ((E_analytic 1).deriv) w hw
    have h2 : AnalyticAt ℂ (fun w : ℂ => (ip - w)⁻¹) w :=
      (analyticAt_const.sub analyticAt_id).inv
        (Complex.slitPlane_ne_zero (c_sub_slit (Or.inl rfl) hw'))
    have h3 : AnalyticAt ℂ (fun w : ℂ => (-ip - w)⁻¹) w :=
      (analyticAt_const.sub analyticAt_id).inv
        (Complex.slitPlane_ne_zero (c_sub_slit (Or.inr rfl) hw'))
    exact h1.add (analyticAt_const.mul (h2.add h3))
  have hQ0 : AnalyticAt ℂ Qfun 0 := hQan 0 (by simp [hpi])
  -- split the `n`-th derivative
  have hsplit : iteratedDeriv (n + 1) (E 1) 0 = iteratedDeriv n Qfun 0 -
      3 / 2 * (iteratedDeriv n (fun w : ℂ => (ip - w)⁻¹) 0 +
        iteratedDeriv n (fun w : ℂ => (-ip - w)⁻¹) 0) := by
    rw [iteratedDeriv_succ', hd, iteratedDeriv_fun_sub (hQ0.contDiffAt.of_le le_top)
      (hPa.contDiffAt.of_le le_top), hP, iteratedDeriv_const_mul_field,
      iteratedDeriv_fun_add (hS1.contDiffAt.of_le le_top) (hS2.contDiffAt.of_le le_top)]
  -- the Cauchy bound for `Qfun`
  set t := Real.pi / ((n : ℝ) + 2) with ht
  have ht0 : 0 < t := by positivity
  have hlt := log_t_le n
  rw [← ht] at hlt
  obtain ⟨b1, b2, b3⟩ := Lam_n_bounds n
  have hMz : ∀ z : ℂ, ‖z‖ = Real.pi * ((n : ℝ) + 1) / ((n : ℝ) + 2) →
      ‖Qfun z‖ ≤ C * (11 + Real.log ((n : ℝ) + 2)) := by
    intro z hz
    obtain ⟨h1, h2, h3⟩ := sphere_dist hz
    refine (hC z h1 t ht0 (t_le_two n) h2 h3).trans ?_
    apply mul_le_mul_of_nonneg_left _ hC0
    linarith
  have hL0 : 0 ≤ 11 + Real.log ((n : ℝ) + 2) := by
    have := Real.log_nonneg (show (1 : ℝ) ≤ n + 2 by linarith); linarith
  have hcau := cauchy_circle hQan n (by positivity) hMz
  rw [hsplit]
  have hS1n := norm_S_iter (Or.inl rfl) n
  have hS2n := norm_S_iter (Or.inr rfl) n
  have hF : (0 : ℝ) < (n.factorial : ℝ) := by exact_mod_cast Nat.factorial_pos n
  have hnorm : ‖iteratedDeriv n Qfun 0 - 3 / 2 * (iteratedDeriv n (fun w : ℂ => (ip - w)⁻¹) 0 +
        iteratedDeriv n (fun w : ℂ => (-ip - w)⁻¹) 0)‖ ≤
      3 * (n.factorial : ℝ) * (C * (11 + Real.log ((n : ℝ) + 2))) / Real.pi ^ n +
        3 * ((n.factorial : ℝ) / Real.pi ^ (n + 1)) := by
    refine (norm_sub_le _ _).trans ?_
    have h32 : ‖(3 / 2 : ℂ)‖ = 3 / 2 := by norm_num
    rw [norm_mul, h32]
    have := norm_add_le (iteratedDeriv n (fun w : ℂ => (ip - w)⁻¹) 0)
      (iteratedDeriv n (fun w : ℂ => (-ip - w)⁻¹) 0)
    rw [hS1n, hS2n] at this
    linarith
  refine hnorm.trans ?_
  simp only [Nat.add_sub_cancel]
  push_cast
  rw [show (n : ℝ) + 1 + 2 = n + 3 by ring]
  have e1 : 3 * (n.factorial : ℝ) * (C * (11 + Real.log ((n : ℝ) + 2))) / Real.pi ^ n =
      3 * Real.pi * C * (n.factorial : ℝ) * (11 + Real.log ((n : ℝ) + 2)) / Real.pi ^ (n + 1) := by
    rw [pow_succ]; field_simp
  rw [e1, show 3 * ((n.factorial : ℝ) / Real.pi ^ (n + 1)) = 3 * (n.factorial : ℝ) / Real.pi ^ (n + 1)
    by ring, ← add_div, div_le_div_iff_of_pos_right (by positivity)]
  have hpC : 0 ≤ 3 * Real.pi * C * (n.factorial : ℝ) := by positivity
  have k1 : 3 * Real.pi * C * (n.factorial : ℝ) * (11 + Real.log ((n : ℝ) + 2)) ≤
      3 * Real.pi * C * (n.factorial : ℝ) * (12 * Real.log ((n : ℝ) + 3)) :=
    mul_le_mul_of_nonneg_left b1 hpC
  have k2 : 3 * ((n.factorial : ℝ)) ≤ 3 * (n.factorial : ℝ) * Real.log ((n : ℝ) + 3) := by
    nlinarith
  nlinarith

theorem bound2 : ∃ C2 : ℝ, ∀ m : ℕ, 1 ≤ m →
    ‖iteratedDeriv m (E 2) 0‖ ≤ C2 * ((m.factorial : ℝ) * Real.log ((m : ℝ) + 2)) / Real.pi ^ m := by
  obtain ⟨C, hC0, hC⟩ := E2_deriv_bound
  refine ⟨456 * Real.pi * C, fun m hm => ?_⟩
  obtain ⟨n, rfl⟩ : ∃ n, m = n + 1 := ⟨m - 1, by omega⟩
  have hpi := Real.pi_pos
  have hn : (0 : ℝ) ≤ n := Nat.cast_nonneg n
  set t := Real.pi / ((n : ℝ) + 2) with ht
  have ht0 : 0 < t := by positivity
  have hlt := log_t_le n
  rw [← ht] at hlt
  have hl0 : 0 ≤ Real.log (2 / t) := Real.log_nonneg (by rw [le_div_iff₀ ht0]; linarith [t_le_two n])
  set Λ := 11 + Real.log ((n : ℝ) + 2) with hΛ
  obtain ⟨b1, b2, b3⟩ := Lam_n_bounds n
  have hΛ0 : 0 ≤ Λ := by have := Real.log_nonneg (show (1 : ℝ) ≤ n + 2 by linarith); linarith
  have hit : 1 / t = ((n : ℝ) + 2) / Real.pi := by rw [ht]; field_simp
  set M := C * (Λ ^ 2 + Λ * (((n : ℝ) + 2) / Real.pi)) with hM
  have hMz : ∀ z : ℂ, ‖z‖ = Real.pi * ((n : ℝ) + 1) / ((n : ℝ) + 2) → ‖deriv (E 2) z‖ ≤ M := by
    intro z hz
    obtain ⟨h1, h2, h3⟩ := sphere_dist hz
    refine (hC z h1 t ht0 (t_le_two n) h2 h3).trans ?_
    apply mul_le_mul_of_nonneg_left _ hC0
    have hΛt : 11 + Real.log (2 / t) ≤ Λ := by linarith
    have hpos : 0 ≤ 11 + Real.log (2 / t) := by linarith
    have e : (11 + Real.log (2 / t)) / t = (11 + Real.log (2 / t)) * (1 / t) := by ring
    rw [e, hit]
    gcongr
  have hcau := cauchy_circle ((E_analytic 2).deriv) n (by positivity) hMz
  rw [iteratedDeriv_succ']
  refine hcau.trans ?_
  -- compare
  have hM2 : M ≤ C * (152 * (((n : ℝ) + 1) * Real.log ((n : ℝ) + 3))) := by
    apply mul_le_mul_of_nonneg_left _ hC0
    have hpi3 := Real.pi_gt_three
    have c1 : Λ ^ 2 ≤ 144 * (((n : ℝ) + 1) * Real.log ((n : ℝ) + 3)) := by
      have := mul_le_mul b2 b1 hΛ0 (by positivity)
      nlinarith
    have c2 : Λ * (((n : ℝ) + 2) / Real.pi) ≤ 8 * (((n : ℝ) + 1) * Real.log ((n : ℝ) + 3)) := by
      have h1 : ((n : ℝ) + 2) / Real.pi ≤ 2 * ((n : ℝ) + 1) / 3 := by
        rw [div_le_div_iff₀ hpi (by norm_num)]; nlinarith
      have := mul_le_mul b1 h1 (by positivity) (by positivity)
      nlinarith
    linarith
  push_cast
  rw [show (n : ℝ) + 1 + 2 = n + 3 by ring, Nat.factorial_succ]
  push_cast
  have hF : (0 : ℝ) < (n.factorial : ℝ) := by exact_mod_cast Nat.factorial_pos n
  rw [div_le_div_iff₀ (by positivity) (by positivity), pow_succ]
  have hL3 : 0 ≤ Real.log ((n : ℝ) + 3) := by linarith
  have key : 3 * (n.factorial : ℝ) * M ≤
      3 * (n.factorial : ℝ) * (C * (152 * (((n : ℝ) + 1) * Real.log ((n : ℝ) + 3)))) :=
    mul_le_mul_of_nonneg_left hM2 (by positivity)
  have hpn : 0 < Real.pi ^ n := by positivity
  calc 3 * (n.factorial : ℝ) * M * (Real.pi ^ n * Real.pi)
      ≤ 3 * (n.factorial : ℝ) * (C * (152 * (((n : ℝ) + 1) * Real.log ((n : ℝ) + 3)))) *
          (Real.pi ^ n * Real.pi) := mul_le_mul_of_nonneg_right key (by positivity)
    _ = 456 * Real.pi * C * ((((n : ℝ) + 1) * (n.factorial : ℝ)) * Real.log ((n : ℝ) + 3)) *
          Real.pi ^ n := by ring

theorem bound3 : ∃ C3 K : ℝ, 0 ≤ K ∧ ∀ i m : ℕ,
    ‖iteratedDeriv m (E i) 0‖ ≤ C3 * (((i : ℝ) + 1) * (m.factorial : ℝ) *
      (Real.log ((m : ℝ) + 2) + K) ^ i) / Real.pi ^ m := by
  obtain ⟨K, hK, hE⟩ := E_pointwise
  refine ⟨3 * K, 11, by norm_num, fun i m => ?_⟩
  have hpi := Real.pi_pos
  have hm : (0 : ℝ) ≤ m := Nat.cast_nonneg m
  set t := Real.pi / ((m : ℝ) + 2) with ht
  have ht0 : 0 < t := by positivity
  have hlt := log_t_le m
  rw [← ht] at hlt
  have hl0 : 0 ≤ Real.log (2 / t) := Real.log_nonneg (by rw [le_div_iff₀ ht0]; linarith [t_le_two m])
  have hL : 0 ≤ Real.log ((m : ℝ) + 2) := Real.log_nonneg (by linarith)
  have hM : ∀ z : ℂ, ‖z‖ = Real.pi * ((m : ℝ) + 1) / ((m : ℝ) + 2) →
      ‖E i z‖ ≤ K * (((i : ℝ) + 1) * (Real.log ((m : ℝ) + 2) + 11) ^ i) := by
    intro z hz
    obtain ⟨h1, h2, h3⟩ := sphere_dist hz
    refine (hE i z h1 t ht0 (t_le_two m) h2 h3).trans ?_
    apply mul_le_mul_of_nonneg_left _ hK
    apply mul_le_mul_of_nonneg_left _ (by positivity)
    apply pow_le_pow_left₀ (by linarith)
    linarith
  have := cauchy_circle (E_analytic i) m (by positivity) hM
  calc _ ≤ 3 * (m.factorial : ℝ) * (K * (((i : ℝ) + 1) * (Real.log ((m : ℝ) + 2) + 11) ^ i)) /
        Real.pi ^ m := this
    _ = _ := by ring

/-! ### The main theorem -/

theorem main : Q06.stmt_2 := by
  obtain ⟨C0, h0⟩ := bound0
  obtain ⟨C1, h1⟩ := bound1
  obtain ⟨C2, h2⟩ := bound2
  obtain ⟨C3, K, hK, h3⟩ := bound3
  obtain ⟨C, N, hC⟩ := assemble (fun i m => ‖iteratedDeriv m (E i) 0‖) C0 C1 C2 C3 K hK
    (fun i m => norm_nonneg _) h0 h1 h2 h3
  refine ⟨C, N + 2, fun k hk => ?_⟩
  obtain ⟨n, rfl⟩ : ∃ n, k = n + 2 := ⟨k - 2, by omega⟩
  have hn : N ≤ n := by omega
  rw [reduction]
  have e1 : n + 2 - 1 = n + 1 := by omega
  have e2 : n + 2 - 2 = n := by omega
  rw [e1, e2]
  have hX : X n = (n.factorial : ℝ) / Real.pi ^ (n + 1) *
      (Real.log ((n + 2 : ℕ) : ℝ) / ((n + 2 : ℕ) : ℝ) ^ 2) := by
    rw [X]; push_cast; ring_nf
  rw [← hX]
  refine (Finset.abs_sum_le_sum_abs _ _).trans ?_
  have hre : ∀ m ∈ Finset.range (n + 1),
      |(iteratedDeriv m (E (n - m)) 0).re / ((n - m).factorial : ℝ)| ≤
        ‖iteratedDeriv m (E (n - m)) 0‖ / ((n - m).factorial : ℝ) := by
    intro m _
    have hf : (0 : ℝ) < ((n - m).factorial : ℝ) := by exact_mod_cast Nat.factorial_pos _
    rw [abs_div, abs_of_pos hf]
    exact div_le_div_of_nonneg_right (Complex.abs_re_le_norm _) hf.le
  refine (Finset.sum_le_sum hre).trans ?_
  have hrefl := Finset.sum_range_reflect
    (fun i => ‖iteratedDeriv (n - i) (E i) 0‖ / (i.factorial : ℝ)) (n + 1)
  have e3 : ∀ j ∈ Finset.range (n + 1), (fun i => ‖iteratedDeriv (n - i) (E i) 0‖ /
      (i.factorial : ℝ)) (n + 1 - 1 - j) = ‖iteratedDeriv j (E (n - j)) 0‖ /
        ((n - j).factorial : ℝ) := by
    intro j hj
    simp only [Finset.mem_range] at hj
    simp only [show n + 1 - 1 - j = n - j by omega, show n - (n - j) = j by omega]
  rw [Finset.sum_congr rfl e3] at hrefl
  rw [hrefl]
  exact hC n hn

end Q06.PB
end

section
-- merged from work/Q06/PartC.lean
/-!
# Q06, part C: `lem:pip:twoterm` from `lem:pip:Tk` and `lem:pip:TkT1`

`Q06.PC.main : Q06.stmt_1 → Q06.stmt_2 → Q06.stmt_3`.

Plan. `T_k^{(1)} = 2 Re Σ_{m ≤ k-2} (iπ)^{-(m+1)} d_{m,k-2-m}` with
`d_{m,i} = [h^i] e^{-hL} (1+h)_m`, `L = Log(iπ) = log π + iπ/2`. With `P_m(h) = (1+h)_m`
(a polynomial with coefficients `p_{m,j}`, `p_{m,0} = m!`, `p_{m,1} = m! H_m`,
`|p_{m,j}| ≤ m! H_m^j / j!`), Leibniz gives `d_{m,i} = Σ_j (-L)^j/j! p_{m,i-j}`, hence
`d_{m,0} = m!`, `d_{m,1} = m!(H_m - L)` and `|d_{m,i}| ≤ m!(|L| + H_m)^i/i!`.
The terms `i = 0, 1` are the two displayed terms (`ψ(n) + γ = H_{n-1}`), and the terms `i ≥ 2`
are summed by a ratio test.
-/

open Filter Topology Polynomial

namespace Q06.PC

/-! ### The polynomial `(1+h)_m` and its coefficients -/

/-- `P_m(h) = (1+h)_m` as a polynomial -/
noncomputable def P (m : ℕ) : ℂ[X] := (ascPochhammer ℂ m).comp (X + 1)

lemma P_eval (m : ℕ) (h : ℂ) : (ascPochhammer ℂ m).eval (1 + h) = (P m).eval h := by
  simp [P, eval_comp, add_comm]

lemma P_succ (m : ℕ) : P (m + 1) = P m * (X + C ((m : ℂ) + 1)) := by
  rw [P, ascPochhammer_succ_right, mul_comp, ← P]
  congr 1
  simp only [add_comp, X_comp, natCast_comp, C_add, C_1, map_natCast]
  ring

lemma P_coeff_succ (m j : ℕ) :
    (P (m + 1)).coeff (j + 1) = (P m).coeff j + (P m).coeff (j + 1) * ((m : ℂ) + 1) := by
  rw [P_succ, mul_add, coeff_add, coeff_mul_X, coeff_mul_C]

lemma P_coeff_succ_zero (m : ℕ) :
    (P (m + 1)).coeff 0 = (P m).coeff 0 * ((m : ℂ) + 1) := by
  rw [P_succ, mul_add, coeff_add, coeff_mul_X_zero, coeff_mul_C, zero_add]

lemma P_zero : P 0 = 1 := by simp [P]

lemma P_coeff0 (m : ℕ) : (P m).coeff 0 = (m.factorial : ℂ) := by
  induction m with
  | zero => simp [P_zero]
  | succ m ih => rw [P_coeff_succ_zero, ih, Nat.factorial_succ]; push_cast; ring

lemma P_coeff1 (m : ℕ) : (P m).coeff 1 = (m.factorial : ℂ) * ((harmonic m : ℝ) : ℂ) := by
  induction m with
  | zero => simp [P_zero, coeff_one]
  | succ m ih =>
    rw [P_coeff_succ, ih, P_coeff0, harmonic_succ, Nat.factorial_succ]
    push_cast
    have : ((m : ℂ) + 1) ≠ 0 := by
      have : (0 : ℝ) < (m : ℝ) + 1 := by positivity
      exact_mod_cast this.ne'
    field_simp
    ring

lemma bern {H δ : ℝ} (hH : 0 ≤ H) (hδ : 0 ≤ δ) (j : ℕ) :
    H ^ (j + 1) + ((j : ℝ) + 1) * H ^ j * δ ≤ (H + δ) ^ (j + 1) := by
  induction j with
  | zero => simp
  | succ j ih =>
    have h1 : 0 ≤ H ^ j := pow_nonneg hH j
    have h2 : 0 ≤ H ^ (j + 1) := pow_nonneg hH _
    calc H ^ (j + 1 + 1) + ((((j + 1 : ℕ) : ℝ)) + 1) * H ^ (j + 1) * δ
        ≤ (H + δ) * (H ^ (j + 1) + ((j : ℝ) + 1) * H ^ j * δ) := by
          push_cast
          have : 0 ≤ ((j : ℝ) + 1) * H ^ j * δ * δ := by positivity
          have e : (H + δ) * (H ^ (j + 1) + ((j : ℝ) + 1) * H ^ j * δ) -
              (H ^ (j + 1 + 1) + ((j : ℝ) + 1 + 1) * H ^ (j + 1) * δ) =
              ((j : ℝ) + 1) * H ^ j * δ * δ := by ring
          linarith
      _ ≤ (H + δ) * (H + δ) ^ (j + 1) := mul_le_mul_of_nonneg_left ih (by positivity)
      _ = (H + δ) ^ (j + 1 + 1) := by ring

lemma harm_nonneg (m : ℕ) : (0 : ℝ) ≤ (harmonic m : ℝ) := by
  rcases m with _ | m
  · simp
  · exact_mod_cast (harmonic_pos (Nat.succ_ne_zero m)).le

lemma P_coeff_bd (m j : ℕ) :
    ‖(P m).coeff j‖ ≤ (m.factorial : ℝ) * (harmonic m : ℝ) ^ j / (j.factorial : ℝ) := by
  induction m generalizing j with
  | zero =>
    rcases j with _ | j
    · simp [P_zero]
    · simp [P_zero, coeff_one]
  | succ m ih =>
    rcases j with _ | j
    · rw [P_coeff0]; simp
    · rw [P_coeff_succ]
      have hm1 : (0 : ℝ) < (m : ℝ) + 1 := by positivity
      have hH := harm_nonneg m
      have hfj : (0 : ℝ) < (j.factorial : ℝ) := by exact_mod_cast Nat.factorial_pos j
      have hfm : (0 : ℝ) < (m.factorial : ℝ) := by exact_mod_cast Nat.factorial_pos m
      have hnorm : ‖((m : ℂ) + 1)‖ = (m : ℝ) + 1 := by
        rw [show ((m : ℂ) + 1) = (((m : ℝ) + 1 : ℝ) : ℂ) by push_cast; ring, Complex.norm_real,
          Real.norm_of_nonneg hm1.le]
      calc ‖(P m).coeff j + (P m).coeff (j + 1) * ((m : ℂ) + 1)‖
          ≤ ‖(P m).coeff j‖ + ‖(P m).coeff (j + 1)‖ * ((m : ℝ) + 1) := by
            refine (norm_add_le _ _).trans ?_
            rw [norm_mul, hnorm]
        _ ≤ (m.factorial : ℝ) * (harmonic m : ℝ) ^ j / (j.factorial : ℝ) +
            (m.factorial : ℝ) * (harmonic m : ℝ) ^ (j + 1) / ((j + 1).factorial : ℝ) *
              ((m : ℝ) + 1) := by
            gcongr
            · exact ih j
            · exact ih (j + 1)
        _ = ((m + 1).factorial : ℝ) / ((j + 1).factorial : ℝ) *
            ((harmonic m : ℝ) ^ (j + 1) + ((j : ℝ) + 1) * (harmonic m : ℝ) ^ j *
              ((m : ℝ) + 1)⁻¹) := by
            rw [Nat.factorial_succ, Nat.factorial_succ]; push_cast
            field_simp
            ring
        _ ≤ ((m + 1).factorial : ℝ) / ((j + 1).factorial : ℝ) *
            ((harmonic m : ℝ) + ((m : ℝ) + 1)⁻¹) ^ (j + 1) := by
            gcongr
            exact bern hH (by positivity) j
        _ = _ := by
            rw [harmonic_succ]; push_cast; ring

/-! ### Taylor coefficients of `e^{-hL} P_m(h)` -/

lemma iteratedDeriv_eval (p : ℂ[X]) (j : ℕ) :
    iteratedDeriv j (fun h : ℂ => p.eval h) = fun h => (derivative^[j] p).eval h := by
  induction j with
  | zero => simp
  | succ j ih =>
    rw [iteratedDeriv_succ, ih]
    funext h
    rw [Polynomial.deriv, Function.iterate_succ_apply']

lemma iteratedDeriv_eval_zero (p : ℂ[X]) (j : ℕ) :
    iteratedDeriv j (fun h : ℂ => p.eval h) 0 = (j.factorial : ℂ) * p.coeff j := by
  rw [iteratedDeriv_eval]
  dsimp only
  rw [← coeff_zero_eq_eval_zero, coeff_iterate_derivative, zero_add,
    Nat.descFactorial_self, nsmul_eq_mul]

/-- `L = Log(iπ)` -/
noncomputable def L : ℂ := Complex.log (Complex.I * (Real.pi : ℂ))

lemma L_eq : L = (Real.log Real.pi : ℂ) + Complex.I * Real.pi / 2 := by
  rw [L, Complex.log_mul_ofReal Real.pi Real.pi_pos Complex.I Complex.I_ne_zero,
    Complex.log_I]
  ring

lemma norm_L_le : ‖L‖ ≤ 4 := by
  rw [L_eq]
  have hpi := Real.pi_pos
  have hp3 := Real.pi_lt_d2
  have hl : Real.log Real.pi ≤ Real.pi - 1 := Real.log_le_sub_one_of_pos hpi
  have hl0 : 0 ≤ Real.log Real.pi := Real.log_nonneg (by linarith [Real.pi_gt_three])
  calc ‖(Real.log Real.pi : ℂ) + Complex.I * Real.pi / 2‖
      ≤ ‖(Real.log Real.pi : ℂ)‖ + ‖Complex.I * Real.pi / 2‖ := norm_add_le _ _
    _ = Real.log Real.pi + Real.pi / 2 := by
      rw [Complex.norm_real, Real.norm_of_nonneg hl0, norm_div, norm_mul, Complex.norm_I,
        Complex.norm_real, Real.norm_of_nonneg hpi.le]
      simp
    _ ≤ 4 := by linarith

/-- `d_{m,i} = [h^i] (e^{-hL} P_m(h))` -/
noncomputable def dcoef (m i : ℕ) : ℂ :=
  iteratedDeriv i (fun h : ℂ => Complex.exp (-L * h) * (P m).eval h) 0 / (i.factorial : ℂ)

lemma dcoef_eq (m i : ℕ) :
    dcoef m i = ∑ j ∈ Finset.range (i + 1),
      (-L) ^ j / (j.factorial : ℂ) * (P m).coeff (i - j) := by
  rw [dcoef]
  have hf : ContDiffAt ℂ i (fun h : ℂ => Complex.exp (-L * h)) 0 := by fun_prop
  have hg : ContDiffAt ℂ i (fun h : ℂ => (P m).eval h) 0 :=
    ((Polynomial.differentiable (P m)).contDiff).contDiffAt
  rw [iteratedDeriv_fun_mul hf hg, Finset.sum_div]
  refine Finset.sum_congr rfl (fun j hj => ?_)
  have hj' : j ≤ i := Nat.lt_succ_iff.mp (Finset.mem_range.mp hj)
  rw [iteratedDeriv_cexp_const_mul, iteratedDeriv_eval_zero]
  simp only [mul_zero, Complex.exp_zero, mul_one]
  have hc := Nat.choose_mul_factorial_mul_factorial hj'
  have hc' : ((i.choose j : ℕ) : ℂ) * (j.factorial : ℂ) * ((i - j).factorial : ℂ) =
      (i.factorial : ℂ) := by exact_mod_cast hc
  have h1 : (j.factorial : ℂ) ≠ 0 := by exact_mod_cast (Nat.factorial_pos j).ne'
  have h2 : (i.factorial : ℂ) ≠ 0 := by exact_mod_cast (Nat.factorial_pos i).ne'
  field_simp
  linear_combination (-L) ^ j * (P m).coeff (i - j) * hc'

lemma dcoef_zero (m : ℕ) : dcoef m 0 = (m.factorial : ℂ) := by
  rw [dcoef_eq]; simp [P_coeff0]

lemma dcoef_one (m : ℕ) :
    dcoef m 1 = (m.factorial : ℂ) * (((harmonic m : ℝ) : ℂ) - L) := by
  rw [dcoef_eq, Finset.sum_range_succ, Finset.sum_range_one]
  simp [P_coeff0, P_coeff1]
  ring

lemma dcoef_bd (m i : ℕ) :
    ‖dcoef m i‖ ≤ (m.factorial : ℝ) * (‖L‖ + (harmonic m : ℝ)) ^ i / (i.factorial : ℝ) := by
  rw [dcoef_eq]
  have hfi : (0 : ℝ) < (i.factorial : ℝ) := by exact_mod_cast Nat.factorial_pos i
  calc ‖∑ j ∈ Finset.range (i + 1), (-L) ^ j / (j.factorial : ℂ) * (P m).coeff (i - j)‖
      ≤ ∑ j ∈ Finset.range (i + 1), ‖(-L) ^ j / (j.factorial : ℂ) * (P m).coeff (i - j)‖ :=
        norm_sum_le _ _
    _ ≤ ∑ j ∈ Finset.range (i + 1), ‖L‖ ^ j / (j.factorial : ℝ) *
          ((m.factorial : ℝ) * (harmonic m : ℝ) ^ (i - j) / ((i - j).factorial : ℝ)) := by
        refine Finset.sum_le_sum (fun j _ => ?_)
        rw [norm_mul, norm_div, norm_pow, norm_neg, Complex.norm_natCast]
        exact mul_le_mul_of_nonneg_left (P_coeff_bd m (i - j)) (by positivity)
    _ = (m.factorial : ℝ) / (i.factorial : ℝ) * ∑ j ∈ Finset.range (i + 1),
          ‖L‖ ^ j * (harmonic m : ℝ) ^ (i - j) * (i.choose j : ℝ) := by
        rw [Finset.mul_sum]
        refine Finset.sum_congr rfl (fun j hj => ?_)
        have hj' : j ≤ i := Nat.lt_succ_iff.mp (Finset.mem_range.mp hj)
        have hc := Nat.choose_mul_factorial_mul_factorial hj'
        have hc' : ((i.choose j : ℕ) : ℝ) * (j.factorial : ℝ) * ((i - j).factorial : ℝ) =
            (i.factorial : ℝ) := by exact_mod_cast hc
        have h1 : (0 : ℝ) < (j.factorial : ℝ) := by exact_mod_cast Nat.factorial_pos j
        have h2 : (0 : ℝ) < ((i - j).factorial : ℝ) := by exact_mod_cast Nat.factorial_pos _
        have h3 : (0 : ℝ) < ((i.choose j : ℕ) : ℝ) := by exact_mod_cast Nat.choose_pos hj'
        rw [← hc']
        field_simp
    _ = _ := by rw [← add_pow]; ring

/-! ### `T_k^{(1)}` through the coefficients `d_{m,i}` -/

lemma T1_eq (k : ℕ) : Q02.T1 k = 2 * (∑ m ∈ Finset.range (k - 1),
    ((Complex.I * Real.pi) ^ (m + 1))⁻¹ * dcoef m (k - 2 - m)).re := by
  unfold Q02.T1
  congr 2
  refine Finset.sum_congr rfl (fun m _ => ?_)
  have hfun : (fun h : ℂ => ((Complex.I * Real.pi) ^ (m + 1))⁻¹ *
      Complex.exp (-(h * Complex.log (Complex.I * Real.pi))) *
      (ascPochhammer ℂ m).eval (1 + h)) =
      fun h => ((Complex.I * Real.pi) ^ (m + 1))⁻¹ * (Complex.exp (-L * h) * (P m).eval h) := by
    funext h
    rw [P_eval, show Complex.log (Complex.I * (Real.pi : ℂ)) = L from rfl,
      show -(h * L) = -L * h by ring, mul_assoc]
  rw [hfun, iteratedDeriv_const_mul_field, dcoef, mul_div_assoc]

lemma psi_nat (m : ℕ) :
    Q02.psi ((m : ℝ) + 1) = (harmonic m : ℝ) - Real.eulerMascheroniConstant := by
  rw [Q02.psi, logDeriv_apply, Real.deriv_Gamma_nat, Real.Gamma_nat_eq_factorial]
  have : (m.factorial : ℝ) ≠ 0 := by positivity
  field_simp
  ring

/-! ### The ratio test for the terms `i ≥ 2` -/

lemma tail_step (n j : ℕ) (hj : j + 1 ≤ n - 2) {x : ℝ} (hx : 0 ≤ x) (hxn : 2 * x ≤ n) :
    ((n - 2 - (j + 1)).factorial : ℝ) * x ^ (2 + (j + 1)) / ((2 + (j + 1)).factorial : ℝ) ≤
      ((n - 2 - j).factorial : ℝ) * x ^ (2 + j) / ((2 + j).factorial : ℝ) / 2 := by
  obtain ⟨t, ht⟩ : ∃ t, n - 2 - (j + 1) = t := ⟨_, rfl⟩
  have ht' : n - 2 - j = t + 1 := by omega
  have hn : (n : ℝ) = (t : ℝ) + 1 + 2 + j := by
    have : n = t + 1 + 2 + j := by omega
    rw [this]; push_cast; ring
  rw [ht, ht', show 2 + (j + 1) = (2 + j) + 1 by ring, Nat.factorial_succ (2 + j),
    Nat.factorial_succ t, pow_succ]
  push_cast
  have hf1 : (0 : ℝ) < (t.factorial : ℝ) := by exact_mod_cast Nat.factorial_pos t
  have hf2 : (0 : ℝ) < ((2 + j).factorial : ℝ) := by exact_mod_cast Nat.factorial_pos _
  have hxp : 0 ≤ x ^ (2 + j) := pow_nonneg hx _
  have hkey : 2 * x ≤ ((t : ℝ) + 1) * ((2 : ℝ) + j + 1) := by
    have ht0 : (0 : ℝ) ≤ t := Nat.cast_nonneg t
    have hj0 : (0 : ℝ) ≤ j := Nat.cast_nonneg j
    nlinarith
  rw [div_div, div_le_div_iff₀ (by positivity) (by positivity)]
  have hA : 0 ≤ (t.factorial : ℝ) * x ^ (2 + j) * ((2 + j).factorial : ℝ) := by positivity
  nlinarith [mul_le_mul_of_nonneg_left hkey hA]

lemma geom_tail (n : ℕ) {x : ℝ} (hx : 0 ≤ x) (hxn : 2 * x ≤ n) :
    ∑ j ∈ Finset.range (n - 1),
      ((n - 2 - j).factorial : ℝ) * x ^ (2 + j) / ((2 + j).factorial : ℝ) ≤
      ((n - 2).factorial : ℝ) * x ^ 2 := by
  set a : ℕ → ℝ := fun j => ((n - 2 - j).factorial : ℝ) * x ^ (2 + j) / ((2 + j).factorial : ℝ)
    with ha
  have hgeo : ∀ j, j ≤ n - 2 → a j ≤ a 0 * (1 / 2) ^ j := by
    intro j
    induction j with
    | zero => intro _; simp
    | succ j ih =>
      intro hj
      have h1 := tail_step n j hj hx hxn
      have h2 := ih (by omega)
      calc a (j + 1) ≤ a j / 2 := h1
        _ ≤ a 0 * (1 / 2) ^ j / 2 := by linarith
        _ = a 0 * (1 / 2) ^ (j + 1) := by ring
  have ha0 : 0 ≤ a 0 := by simp only [ha]; positivity
  calc ∑ j ∈ Finset.range (n - 1), a j ≤ ∑ j ∈ Finset.range (n - 1), a 0 * (1 / 2) ^ j :=
        Finset.sum_le_sum (fun j hj => hgeo j (by have := Finset.mem_range.mp hj; omega))
    _ = a 0 * ∑ j ∈ Finset.range (n - 1), (1 / 2 : ℝ) ^ j := by rw [Finset.mul_sum]
    _ ≤ a 0 * 2 := mul_le_mul_of_nonneg_left (sum_geometric_two_le _) ha0
    _ = ((n - 2).factorial : ℝ) * x ^ 2 := by simp only [ha]; norm_num

/-! ### Splitting `T_{n+2}^{(1)}` into the two main terms and the rest -/

/-- the term `i = n - m` of `T_{n+2}^{(1)}` -/
noncomputable def g (n i : ℕ) : ℂ := ((Complex.I * Real.pi) ^ (n - i + 1))⁻¹ * dcoef (n - i) i

/-- the two-term main part, `k = n + 2` -/
noncomputable def M (n : ℕ) : ℂ :=
  (n.factorial : ℂ) * ((Complex.I * Real.pi) ^ (n + 1))⁻¹ +
    ((n - 1).factorial : ℂ) * ((Complex.I * Real.pi) ^ n)⁻¹ *
      (((Q02.psi (n : ℝ) + Real.eulerMascheroniConstant - Real.log Real.pi : ℝ) : ℂ) -
        Complex.I * Real.pi / 2)

lemma sum_eq (n : ℕ) (hn : 1 ≤ n) :
    ∑ m ∈ Finset.range (n + 2 - 1),
      ((Complex.I * Real.pi) ^ (m + 1))⁻¹ * dcoef m (n + 2 - 2 - m) =
      M n + ∑ j ∈ Finset.range (n - 1), g n (2 + j) := by
  have e1 : n + 2 - 1 = n + 1 := by omega
  have e2 : n + 2 - 2 = n := by omega
  rw [e1, e2, ← Finset.sum_range_reflect]
  have h3 : ∀ i ∈ Finset.range (n + 1),
      ((Complex.I * Real.pi) ^ (n + 1 - 1 - i + 1))⁻¹ * dcoef (n + 1 - 1 - i) (n - (n + 1 - 1 - i))
        = g n i := by
    intro i hi
    have hi' := Finset.mem_range.mp hi
    rw [g, show n + 1 - 1 - i = n - i by omega, show n - (n - i) = i by omega]
  rw [Finset.sum_congr rfl h3, show n + 1 = 2 + (n - 1) by omega, Finset.sum_range_add]
  congr 1
  rw [Finset.sum_range_succ, Finset.sum_range_one]
  obtain ⟨p, rfl⟩ : ∃ p, n = p + 1 := ⟨n - 1, by omega⟩
  simp only [g, dcoef_zero, dcoef_one, M, Nat.sub_zero, Nat.add_sub_cancel]
  rw [show ((p + 1 : ℕ) : ℝ) = (p : ℝ) + 1 by push_cast; ring, psi_nat, L_eq]
  push_cast
  ring

lemma norm_c (p : ℕ) : ‖((Complex.I * (Real.pi : ℂ)) ^ p)⁻¹‖ = (Real.pi ^ p)⁻¹ := by
  rw [norm_inv, norm_pow, Q02.norm_ipi]

lemma rem_bd (n : ℕ) (b : ℝ) (hb : ∀ m ≤ n, ‖L‖ + (harmonic m : ℝ) ≤ b) (hb0 : 0 ≤ b)
    (hbn : 2 * (Real.pi * b) ≤ n) :
    ‖∑ j ∈ Finset.range (n - 1), g n (2 + j)‖ ≤
      ((n - 2).factorial : ℝ) * (Real.pi * b) ^ 2 / Real.pi ^ (n + 1) := by
  have hpi := Real.pi_pos
  calc ‖∑ j ∈ Finset.range (n - 1), g n (2 + j)‖
      ≤ ∑ j ∈ Finset.range (n - 1), ‖g n (2 + j)‖ := norm_sum_le _ _
    _ ≤ ∑ j ∈ Finset.range (n - 1), ((n - 2 - j).factorial : ℝ) * (Real.pi * b) ^ (2 + j) /
          ((2 + j).factorial : ℝ) / Real.pi ^ (n + 1) := by
        refine Finset.sum_le_sum (fun j hj => ?_)
        have hj' := Finset.mem_range.mp hj
        have e : n - (2 + j) = n - 2 - j := by omega
        rw [g, norm_mul, norm_c, e]
        have hd := dcoef_bd (n - 2 - j) (2 + j)
        have hH0 := harm_nonneg (n - 2 - j)
        have hH : (‖L‖ + (harmonic (n - 2 - j) : ℝ)) ^ (2 + j) ≤ b ^ (2 + j) :=
          pow_le_pow_left₀ (by positivity) (hb _ (by omega)) _
        have hpow : Real.pi ^ (n + 1) = Real.pi ^ (n - 2 - j + 1) * Real.pi ^ (2 + j) := by
          rw [← pow_add]; congr 1; omega
        have hf1 : (0 : ℝ) < ((n - 2 - j).factorial : ℝ) := by exact_mod_cast Nat.factorial_pos _
        have hf2 : (0 : ℝ) < ((2 + j).factorial : ℝ) := by exact_mod_cast Nat.factorial_pos _
        calc (Real.pi ^ (n - 2 - j + 1))⁻¹ * ‖dcoef (n - 2 - j) (2 + j)‖
            ≤ (Real.pi ^ (n - 2 - j + 1))⁻¹ *
              (((n - 2 - j).factorial : ℝ) * b ^ (2 + j) / ((2 + j).factorial : ℝ)) := by
              refine mul_le_mul_of_nonneg_left (hd.trans ?_) (by positivity)
              gcongr
          _ = _ := by
              rw [hpow, mul_pow]
              field_simp
    _ = (∑ j ∈ Finset.range (n - 1), ((n - 2 - j).factorial : ℝ) * (Real.pi * b) ^ (2 + j) /
          ((2 + j).factorial : ℝ)) / Real.pi ^ (n + 1) := by rw [Finset.sum_div]
    _ ≤ _ := by
        gcongr
        exact geom_tail n (by positivity) hbn

/-! ### The scale `(k-2)! π^{-(k-1)} log²k / k²` and the error terms -/

/-- the scale, `k = n + 2` -/
noncomputable def S (n : ℕ) : ℝ :=
  (n.factorial : ℝ) / Real.pi ^ (n + 1) *
    (Real.log ((n + 2 : ℕ) : ℝ) ^ 2 / ((n + 2 : ℕ) : ℝ) ^ 2)

lemma log_ge_one {n : ℕ} (hn : 1 ≤ n) : 1 ≤ Real.log ((n + 2 : ℕ) : ℝ) := by
  rw [Real.le_log_iff_exp_le (by positivity)]
  have := Real.exp_one_lt_d9
  have h3 : (3 : ℝ) ≤ ((n + 2 : ℕ) : ℝ) := by exact_mod_cast (by omega : 3 ≤ n + 2)
  linarith

lemma S_eq (n : ℕ) : S n = (n.factorial : ℝ) * Real.log ((n + 2 : ℕ) : ℝ) ^ 2 /
    (Real.pi ^ (n + 1) * ((n + 2 : ℕ) : ℝ) ^ 2) := by
  rw [S, div_mul_div_comm]

lemma S_nonneg (n : ℕ) : 0 ≤ S n := by
  rw [S]; have := Real.pi_pos; positivity

/-- `f ≤ S n` as soon as `f · π^{n+1} (n+2)² ≤ n!` -/
lemma le_S {n : ℕ} (hn : 1 ≤ n) {f : ℝ}
    (hf : f * (Real.pi ^ (n + 1) * ((n + 2 : ℕ) : ℝ) ^ 2) ≤ (n.factorial : ℝ)) : f ≤ S n := by
  have hpi := Real.pi_pos
  rw [S_eq, le_div_iff₀ (by positivity)]
  have hl := log_ge_one hn
  have hl2 : 1 ≤ Real.log ((n + 2 : ℕ) : ℝ) ^ 2 := one_le_pow₀ hl
  have hf0 : (0 : ℝ) ≤ (n.factorial : ℝ) := by positivity
  nlinarith

lemma ev_a : ∀ᶠ n : ℕ in atTop, (2 : ℝ) ^ (n + 2) ≤ S n := by
  have hpi := Real.pi_pos
  have ht := FloorSemiring.tendsto_pow_div_factorial_atTop (8 * Real.pi)
  have hev := ht.eventually (gt_mem_nhds (show (0 : ℝ) < 1 / (64 * Real.pi) by positivity))
  filter_upwards [hev, eventually_ge_atTop 1] with n hn hn1
  apply le_S hn1
  have hf : (0 : ℝ) < (n.factorial : ℝ) := by exact_mod_cast Nat.factorial_pos n
  rw [div_lt_div_iff₀ hf (by positivity), one_mul] at hn
  have hsq : ((n + 2 : ℕ) : ℝ) ^ 2 ≤ (4 : ℝ) ^ (n + 2) := by
    have h1 : n + 2 < 2 ^ (n + 2) := Nat.lt_two_pow_self
    have h2 : (n + 2) ^ 2 ≤ 4 ^ (n + 2) := by
      calc (n + 2) ^ 2 ≤ (2 ^ (n + 2)) ^ 2 := Nat.pow_le_pow_left h1.le 2
        _ = 4 ^ (n + 2) := by rw [← pow_mul, mul_comm, pow_mul]; norm_num
    exact_mod_cast h2
  calc (2 : ℝ) ^ (n + 2) * (Real.pi ^ (n + 1) * ((n + 2 : ℕ) : ℝ) ^ 2)
      ≤ (2 : ℝ) ^ (n + 2) * (Real.pi ^ (n + 1) * (4 : ℝ) ^ (n + 2)) := by gcongr
    _ = 64 * Real.pi * (8 * Real.pi) ^ n := by
        rw [show (8 : ℝ) * Real.pi = 2 * 4 * Real.pi by norm_num, mul_pow, mul_pow]; ring
    _ ≤ (n.factorial : ℝ) := by linarith

lemma ev_b : ∀ᶠ n : ℕ in atTop,
    ((n + 2 : ℕ) : ℝ) * ((n + 2 - 2).factorial : ℝ) / 6 ^ (n + 2) ≤ S n := by
  have hpi := Real.pi_pos
  have hr : |Real.pi / 6| < 1 := by
    rw [abs_of_pos (by positivity)]; have := Real.pi_lt_d2; linarith
  have ht := (tendsto_pow_const_mul_const_pow_of_abs_lt_one 3 hr).comp (tendsto_add_atTop_nat 2)
  have hev := ht.eventually (gt_mem_nhds (show (0 : ℝ) < 1 / Real.pi by positivity))
  filter_upwards [hev, eventually_ge_atTop 1] with n hn hn1
  simp only [Function.comp] at hn
  rw [Nat.add_sub_cancel]
  apply le_S hn1
  have hf : (0 : ℝ) ≤ (n.factorial : ℝ) := by positivity
  have h6 : (0 : ℝ) < 6 ^ (n + 2) := by positivity
  have hk : (((n + 2 : ℕ) : ℝ)) = ((n + 2 : ℕ) : ℝ) := rfl
  have key : ((n + 2 : ℕ) : ℝ) ^ 3 * Real.pi ^ (n + 1) ≤ 6 ^ (n + 2) := by
    have e : ((n + 2 : ℕ) : ℝ) ^ 3 * (Real.pi / 6) ^ (n + 2) * 6 ^ (n + 2) =
        ((n + 2 : ℕ) : ℝ) ^ 3 * Real.pi ^ (n + 1) * Real.pi := by
      rw [div_pow]; field_simp; ring
    have h1 : ((n + 2 : ℕ) : ℝ) ^ 3 * (Real.pi / 6) ^ (n + 2) * 6 ^ (n + 2) ≤
        1 / Real.pi * 6 ^ (n + 2) := by
      push_cast at hn ⊢
      exact mul_le_mul_of_nonneg_right hn.le h6.le
    rw [e] at h1
    have h3 : (3 : ℝ) ≤ Real.pi := Real.pi_gt_three.le
    have hA : 0 ≤ ((n + 2 : ℕ) : ℝ) ^ 3 * Real.pi ^ (n + 1) := by positivity
    rw [div_mul_eq_mul_div, one_mul, le_div_iff₀ hpi] at h1
    have hpp : 1 ≤ Real.pi * Real.pi := by nlinarith
    have := mul_le_mul_of_nonneg_left hpp hA
    linarith [show ((n + 2 : ℕ) : ℝ) ^ 3 * Real.pi ^ (n + 1) * (Real.pi * Real.pi) =
      ((n + 2 : ℕ) : ℝ) ^ 3 * Real.pi ^ (n + 1) * Real.pi * Real.pi by ring]
  rw [div_mul_eq_mul_div, div_le_iff₀ h6]
  calc ((n + 2 : ℕ) : ℝ) * (n.factorial : ℝ) * (Real.pi ^ (n + 1) * ((n + 2 : ℕ) : ℝ) ^ 2)
      = (n.factorial : ℝ) * (((n + 2 : ℕ) : ℝ) ^ 3 * Real.pi ^ (n + 1)) := by ring
    _ ≤ (n.factorial : ℝ) * 6 ^ (n + 2) := mul_le_mul_of_nonneg_left key hf

lemma ev_c : ∀ᶠ n : ℕ in atTop,
    2 * (Real.pi * (5 + Real.log ((n + 2 : ℕ) : ℝ))) ≤ (n : ℝ) := by
  have hpi := Real.pi_pos
  have hp := Real.pi_lt_d2
  have hb := Real.isLittleO_log_id_atTop.bound (show (0 : ℝ) < 1 / 14 by norm_num)
  have ht : Tendsto (fun n : ℕ => ((n + 2 : ℕ) : ℝ)) atTop atTop :=
    tendsto_natCast_atTop_atTop.comp (tendsto_add_atTop_nat 2)
  filter_upwards [ht.eventually hb, eventually_ge_atTop 60] with n hn hn60
  simp only [id, Real.norm_eq_abs] at hn
  have hk0 : (0 : ℝ) ≤ ((n + 2 : ℕ) : ℝ) := by positivity
  rw [abs_of_nonneg hk0] at hn
  have hl := (le_abs_self _).trans hn
  have h60 : (60 : ℝ) ≤ (n : ℝ) := by exact_mod_cast hn60
  have hk : ((n + 2 : ℕ) : ℝ) = (n : ℝ) + 2 := by push_cast; ring
  rw [hk] at hl ⊢
  have h1 : Real.pi * ((n : ℝ) + 2) ≤ 3.15 * ((n : ℝ) + 2) :=
    mul_le_mul_of_nonneg_right hp.le (by positivity)
  have h2 := mul_le_mul_of_nonneg_left hl hpi.le
  linarith

lemma e_bd (k : ℕ) : |(-Real.log 2) ^ (k + 1) / ((k + 1).factorial : ℝ)| ≤ 1 := by
  have hl0 : 0 ≤ Real.log 2 := Real.log_nonneg (by norm_num)
  have hl1 : Real.log 2 ≤ 1 := by have := Real.log_two_lt_d9; linarith
  have hf : (1 : ℝ) ≤ ((k + 1).factorial : ℝ) := Nat.one_le_cast.mpr (Nat.factorial_pos _)
  rw [abs_div, abs_pow, abs_neg, abs_of_nonneg hl0, abs_of_pos (by linarith)]
  rw [div_le_one (by linarith)]
  exact (pow_le_one₀ hl0 hl1).trans hf

lemma rem_scale {n : ℕ} (hn : 4 ≤ n) {b : ℝ} (hb0 : 0 ≤ b)
    (hb : b ≤ 6 * Real.log ((n + 2 : ℕ) : ℝ)) :
    2 * (((n - 2).factorial : ℝ) * (Real.pi * b) ^ 2 / Real.pi ^ (n + 1)) ≤
      288 * Real.pi ^ 2 * S n := by
  have hpi := Real.pi_pos
  obtain ⟨q, rfl⟩ : ∃ q, n = q + 2 := ⟨n - 2, by omega⟩
  rw [show q + 2 - 2 = q by omega, S_eq]
  have hq : (2 : ℝ) ≤ q := by exact_mod_cast (by omega : 2 ≤ q)
  have hfac : ((q + 2).factorial : ℝ) = ((q : ℝ) + 2) * ((q : ℝ) + 1) * (q.factorial : ℝ) := by
    rw [Nat.factorial_succ, Nat.factorial_succ]; push_cast; ring
  have hk : ((q + 2 + 2 : ℕ) : ℝ) = (q : ℝ) + 4 := by push_cast; ring
  generalize hl : Real.log ((q + 2 + 2 : ℕ) : ℝ) = l at hb ⊢
  rw [hfac, hk]
  have hl0 : 0 ≤ l := by
    rw [← hl]; exact Real.log_nonneg (by exact_mod_cast (by omega : 1 ≤ q + 2 + 2))
  have hb2 : b ^ 2 ≤ 36 * l ^ 2 := by nlinarith
  have hf : (0 : ℝ) < (q.factorial : ℝ) := by exact_mod_cast Nat.factorial_pos q
  have hP : (0 : ℝ) < Real.pi ^ (q + 2 + 1) := by positivity
  have h1 : 2 * b ^ 2 * ((q : ℝ) + 4) ^ 2 ≤ 288 * (((q : ℝ) + 2) * ((q : ℝ) + 1)) * l ^ 2 := by
    have hr : ((q : ℝ) + 4) ^ 2 ≤ 4 * (((q : ℝ) + 2) * ((q : ℝ) + 1)) := by nlinarith
    have hq4 : 0 ≤ ((q : ℝ) + 4) ^ 2 := by positivity
    calc 2 * b ^ 2 * ((q : ℝ) + 4) ^ 2 ≤ 2 * (36 * l ^ 2) * ((q : ℝ) + 4) ^ 2 := by gcongr
      _ = 72 * l ^ 2 * ((q : ℝ) + 4) ^ 2 := by ring
      _ ≤ 72 * l ^ 2 * (4 * (((q : ℝ) + 2) * ((q : ℝ) + 1))) := by gcongr
      _ = _ := by ring
  rw [show 2 * ((q.factorial : ℝ) * (Real.pi * b) ^ 2 / Real.pi ^ (q + 2 + 1)) =
      (Real.pi ^ 2 * (q.factorial : ℝ)) * (2 * b ^ 2) / Real.pi ^ (q + 2 + 1) by ring,
    show 288 * Real.pi ^ 2 * (((q : ℝ) + 2) * ((q : ℝ) + 1) * (q.factorial : ℝ) * l ^ 2 /
        (Real.pi ^ (q + 2 + 1) * ((q : ℝ) + 4) ^ 2)) =
      (Real.pi ^ 2 * (q.factorial : ℝ)) * (288 * (((q : ℝ) + 2) * ((q : ℝ) + 1)) * l ^ 2) /
        (Real.pi ^ (q + 2 + 1) * ((q : ℝ) + 4) ^ 2) by ring,
    div_le_div_iff₀ hP (by positivity)]
  have hA : 0 ≤ Real.pi ^ 2 * (q.factorial : ℝ) * Real.pi ^ (q + 2 + 1) := by positivity
  have := mul_le_mul_of_nonneg_left h1 hA
  nlinarith

/-! ### Assembly -/

lemma T1_split (n : ℕ) (hn : 1 ≤ n) :
    Q02.T1 (n + 2) - 2 * (M n).re = 2 * (∑ j ∈ Finset.range (n - 1), g n (2 + j)).re := by
  rw [T1_eq, sum_eq n hn, Complex.add_re]
  ring

lemma harm_le {n m : ℕ} (hm : m ≤ n) :
    ‖L‖ + (harmonic m : ℝ) ≤ 5 + Real.log ((n + 2 : ℕ) : ℝ) := by
  have h1 := harmonic_le_one_add_log m
  have h2 : Real.log (m : ℝ) ≤ Real.log ((n + 2 : ℕ) : ℝ) := by
    rcases Nat.eq_zero_or_pos m with h | h
    · rw [h]; simp only [Nat.cast_zero, Real.log_zero]
      exact Real.log_nonneg (by exact_mod_cast (by omega : 1 ≤ n + 2))
    · exact Real.log_le_log (by exact_mod_cast h) (by exact_mod_cast (by omega : m ≤ n + 2))
  have h3 := norm_L_le
  linarith

lemma main_n (h1 : Q06.stmt_1) (h2 : Q06.stmt_2) : ∃ C : ℝ, ∃ n₀ : ℕ, ∀ n : ℕ, n₀ ≤ n →
    |Q02.N (n + 2) - 2 * (M n).re| ≤ C * S n := by
  obtain ⟨C1, hC1⟩ := h1
  obtain ⟨C2, k2, hC2⟩ := h2
  obtain ⟨n₁, hn₁⟩ := Filter.eventually_atTop.mp ((ev_a.and ev_b).and ev_c)
  refine ⟨2 * |C1| + 1 + |C2| + 288 * Real.pi ^ 2, max n₁ (max k2 60), fun n hn => ?_⟩
  have hn1 : n₁ ≤ n := le_of_max_le_left hn
  have hk2 : k2 ≤ n := (le_of_max_le_left (le_of_max_le_right hn))
  have hn60 : 60 ≤ n := (le_of_max_le_right (le_of_max_le_right hn))
  obtain ⟨⟨ha, hb⟩, hc⟩ := hn₁ n hn1
  have hpi := Real.pi_pos
  have hS := S_nonneg n
  -- lem:pip:Tk
  have H1 := hC1 (n + 2) (by omega)
  rw [Nat.add_sub_cancel] at H1
  have hE0 : 0 ≤ ((n + 2 : ℕ) : ℝ) * (n.factorial : ℝ) / 6 ^ (n + 2) := by positivity
  rw [Nat.add_sub_cancel] at hb
  have H1' : C1 * (2 ^ (n + 2) + ((n + 2 : ℕ) : ℝ) * (n.factorial : ℝ) / 6 ^ (n + 2)) ≤
      |C1| * (2 * S n) := by
    calc C1 * (2 ^ (n + 2) + ((n + 2 : ℕ) : ℝ) * (n.factorial : ℝ) / 6 ^ (n + 2))
        ≤ |C1| * (2 ^ (n + 2) + ((n + 2 : ℕ) : ℝ) * (n.factorial : ℝ) / 6 ^ (n + 2)) :=
          mul_le_mul_of_nonneg_right (le_abs_self _) (by positivity)
      _ ≤ |C1| * (2 * S n) := by
          apply mul_le_mul_of_nonneg_left _ (abs_nonneg _); linarith
  -- the term `(-log 2)^{k+1}/(k+1)!`
  have he := e_bd (n + 2)
  have h2k : (1 : ℝ) ≤ 2 ^ (n + 2) := one_le_pow₀ (by norm_num)
  -- lem:pip:TkT1
  have H2 := hC2 (n + 2) (by omega)
  rw [show n + 2 - 2 = n by omega, show n + 2 - 1 = n + 1 by omega] at H2
  have hl1 := log_ge_one (n := n) (by omega)
  have H2' : C2 * ((n.factorial : ℝ) / Real.pi ^ (n + 1) *
      (Real.log ((n + 2 : ℕ) : ℝ) / ((n + 2 : ℕ) : ℝ) ^ 2)) ≤ |C2| * S n := by
    rw [S]
    apply (mul_le_mul_of_nonneg_right (le_abs_self C2) (by positivity)).trans
    apply mul_le_mul_of_nonneg_left _ (abs_nonneg _)
    apply mul_le_mul_of_nonneg_left _ (by positivity)
    apply div_le_div_of_nonneg_right _ (by positivity)
    nlinarith
  -- the terms `i ≥ 2`
  set b := 5 + Real.log ((n + 2 : ℕ) : ℝ) with hbdef
  have hb0 : 0 ≤ b := by rw [hbdef]; linarith
  have hR := rem_bd n b (fun m hm => harm_le hm) hb0 hc
  have hR' := rem_scale (n := n) (by omega) hb0 (by rw [hbdef]; linarith)
  have hsplit := T1_split n (by omega)
  have hre : |2 * (∑ j ∈ Finset.range (n - 1), g n (2 + j)).re| ≤
      2 * ‖∑ j ∈ Finset.range (n - 1), g n (2 + j)‖ := by
    rw [abs_mul, abs_two]
    exact mul_le_mul_of_nonneg_left (Complex.abs_re_le_norm _) (by norm_num)
  -- combine
  have e : Q02.N (n + 2) - 2 * (M n).re =
      (Q02.N (n + 2) - (-Real.log 2) ^ (n + 2 + 1) / ((n + 2 + 1).factorial : ℝ) - Q02.T (n + 2)) +
      (-Real.log 2) ^ (n + 2 + 1) / ((n + 2 + 1).factorial : ℝ) +
      (Q02.T (n + 2) - Q02.T1 (n + 2)) + (Q02.T1 (n + 2) - 2 * (M n).re) := by ring
  rw [e, hsplit]
  set A := Q02.N (n + 2) - (-Real.log 2) ^ (n + 2 + 1) / ((n + 2 + 1).factorial : ℝ) -
    Q02.T (n + 2)
  set B := (-Real.log 2) ^ (n + 2 + 1) / ((n + 2 + 1).factorial : ℝ)
  set D := Q02.T (n + 2) - Q02.T1 (n + 2)
  set E := 2 * (∑ j ∈ Finset.range (n - 1), g n (2 + j)).re
  have t1 := abs_add_le (A + B + D) E
  have t2 := abs_add_le (A + B) D
  have t3 := abs_add_le A B
  have := abs_nonneg C1
  have hfin : |A| + |B| + |D| + |E| ≤
      |C1| * (2 * S n) + S n + |C2| * S n + 288 * Real.pi ^ 2 * S n := by linarith
  nlinarith

end Q06.PC

/-- Part C: Lemma `lem:pip:twoterm` from Lemmas `lem:pip:Tk` (`stmt_1`) and `lem:pip:TkT1`
(`stmt_2`). -/
theorem Q06.PC.main : Q06.stmt_1 → Q06.stmt_2 → Q06.stmt_3 := by
  intro h1 h2
  obtain ⟨C, n₀, hC⟩ := Q06.PC.main_n h1 h2
  refine ⟨C, n₀ + 2, fun k hk => ?_⟩
  obtain ⟨n, rfl⟩ : ∃ n, k = n + 2 := ⟨k - 2, by omega⟩
  have H := hC n (by omega)
  have e1 : n + 2 - 2 = n := by omega
  have e2 : n + 2 - 1 = n + 1 := by omega
  have e3 : n + 2 - 3 = n - 1 := by omega
  have e4 : ((n + 2 : ℕ) : ℝ) - 2 = (n : ℝ) := by push_cast; ring
  rw [e1, e2, e3, e4]
  simpa only [Q06.PC.M, Q06.PC.S] using H
end

theorem Q06.proof_1 : Q06.stmt_1 := Q06.PA.main
theorem Q06.proof_2 : Q06.stmt_2 := Q06.PB.main
theorem Q06.proof_3 : Q06.stmt_3 := Q06.PC.main Q06.proof_1 Q06.proof_2
/-- Lemma `prop:pip:growth`, odd `k`: the real part of the two-term formula of `stmt_3`. -/
theorem Q06.proof_4 : Q06.stmt_4 := Q06.growth_odd Q02.N Q02.psi Q06.proof_3
/-- Lemma `prop:pip:growth`, even `k`: the real part of the two-term formula of `stmt_3`. -/
theorem Q06.proof_5 : Q06.stmt_5 := Q06.growth_even Q02.N Q02.psi Q06.proof_3
/-- Lemma `prop:pip:div` from the growth law for odd `k` (`stmt_4`): `|N_k| ≥ (k-2)!/π^{k-1}`
for large odd `k`, which beats every geometric sequence. -/
theorem Q06.proof_6 : Q06.stmt_6 := Q06.div_of_growth Q02.N Q06.proof_4
theorem Q06.proof_7 : Q06.stmt_7 := ⟨Q02.proof_8, Q06.proof_6, Q06.proof_4, Q06.proof_5⟩
