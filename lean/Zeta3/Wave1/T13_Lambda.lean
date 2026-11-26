import Mathlib
/-!
# T13 -- the Legendre integrals `Λ_m(θ)` of Theorem `thm:P2`

Target id: T13.
Paper label: `thm:P2` (paper/ieee/main.tex, lines 1423-1434), the clause about `Λ_m` and the
clause "At θ = π/2 this is Ahmed's exact formula". Setup: main.tex lines 1418-1421 (definition of
`Λ_ν`). Supporting lemmas in the proofs supplement paper/ieee/proofs/sec_g3_P2.tex:
`g3:lem:Lam` (lines 36-43) and `g3:lem:pi2` (lines 150-161). Also main.tex line 1458
(`Λ_1 = ½ sin²θ`).

Quoted verbatim (main.tex 1418-1421):

  The key object is the Legendre integral, with $P_\nu$ the Legendre function~\cite[Ch.~14]{DLMF}:
  \[
  \Lambda_\nu(\theta)=\int_{\cos\theta}^1P_\nu(t)\dd t=2\sin^2\tfrac\theta2\;{}_2F_1\bigl(-\nu,1+\nu;2;\sin^2\tfrac\theta2\bigr).
  \]

Quoted verbatim (main.tex 1433, last sentences of thm:P2):

  where $\Lambda_m(\theta)=\frac{P_{m-1}(\cos\theta)-P_{m+1}(\cos\theta)}{2m+1}$ is a polynomial in
  $\cos\theta$. At $\theta=\pi/2$ this is Ahmed's exact formula.

Quoted verbatim (main.tex 1458, proof of thm:P2asym):

  because $\Lambda_1=\frac12\sin^2\theta=2z(1-z)$.

Quoted verbatim (proofs/sec_g3_P2.tex 150-154, Lemma g3:lem:pi2):

  At $\theta=\pi/2$ ($c=0$, $z=\frac12$): $\Lambda_m(\frac\pi2)=0$ for even $m$, and
  \[
  \frac{\pi\varepsilon}{\sin\pi\varepsilon}\Lambda_\varepsilon\bigl(\tfrac\pi2\bigr)=\frac{\Gamma(1+\frac\varepsilon2)\,\Gamma(\frac{1-\varepsilon}2)}{\sqrt\pi\,(1+\varepsilon)} .
  \]

and (line 161) "with $\Lambda_m(\frac\pi2)=\frac{P_{m-1}(0)-P_{m+1}(0)}{2m+1}$".

The main paper does not print Ahmed's exact formula; the supplement reduces the θ = π/2 case to
the Gamma closed form above plus the values `Λ_m(π/2)`. So the "Ahmed" part is formalized as
(i) the Gamma closed form of the `ε`-term (stmt_8) and (ii) explicit values of `Λ_m(π/2)` through
the Legendre values at `0` (stmt_5, stmt_6, stmt_7).

## Conventions

* Legendre polynomials: Mathlib has only `Polynomial.shiftedLegendre n : ℤ[X]`, which is
  `Σ (-1)^k C(n,k) C(n+k,n) x^k = P_n(1-2x)`. We define `P n := shiftedLegendre n ∘ ((1-X)/2)`
  over `ℝ`, and stmt_1 checks that this is the classical Legendre polynomial (`P_0 = 1`,
  `P_1 = X`, Bonnet's recurrence), so the definition is not taken on trust.
* For integer `m`, `Λ_m(θ)` is the paper's integral `∫_{cos θ}^1 P_m(t) dt` (interval integral).
* For non-integer `ν` there is no Legendre function in Mathlib, so `Λ_ν(θ)` is taken from the
  paper's second (hypergeometric) expression `2z ₂F₁(-ν,1+ν;2;z)`, `z = sin²(θ/2)`, written as an
  explicit series with `HasSum` (no `tsum` junk values). stmt_4 checks that this series agrees
  with the integral for integer `ν = m`.
* Hypotheses: stmt_2, stmt_3, stmt_9 are stated for all real `θ` (the paper assumes
  `0 ≤ θ < π`; the polynomial identities hold for all `θ`, so this is a weakening of the
  hypotheses, not a strengthening). stmt_4 keeps the paper's `0 ≤ θ < π`.
* stmt_8 is stated for real `ε > 0` not a natural number (the paper's `ε = 1/(n-1)` with the
  removable points `ε ∈ ℕ` excluded; at those points `sin(πε) = 0` and Lean's `x/0 = 0` would
  make the left side meaningless).

## Statements

* stmt_1: `P` is the Legendre family: `P_0 = 1`, `P_1 = X`,
  `(n+2) P_{n+2} = (2n+3) X P_{n+1} - (n+1) P_n`.
* stmt_2: for `m ≥ 1`, `∫_{cos θ}^1 P_m = (P_{m-1}(cos θ) - P_{m+1}(cos θ))/(2m+1)`.
* stmt_3: for `m ≥ 1`, `Λ_m(θ)` is a polynomial with rational coefficients in `cos θ`.
* stmt_4: for integer `m`, `0 ≤ θ < π`, the hypergeometric series `2z ₂F₁(-m,1+m;2;z)` sums to
  `∫_{cos θ}^1 P_m`.
* stmt_5: Legendre values at `0`: `P_{2k+1}(0) = 0`, `P_{2k}(0) = (-1)^k C(2k,k)/4^k`.
* stmt_6: `Λ_m(π/2) = 0` for even `m ≥ 2`.
* stmt_7: `Λ_{2k+1}(π/2) = (-1)^k C(2k,k) / (4^k (2k+2))` (explicit; derived from stmt_2 and
  stmt_5, not printed in the paper).
* stmt_8: the `ε`-term at `θ = π/2`: `(πε/sin πε) Λ_ε(π/2) = Γ(1+ε/2)Γ((1-ε)/2)/(√π (1+ε))`.
* stmt_9: `Λ_1(θ) = ½ sin²θ`.

## NOT COVERED

* The Mellin identity `∫_0^∞ (√(1+2cs+s²) - s - c) s^{ε-1} ds = (π/sin πε) Λ_ε(θ)` and the exact
  length series for `L_θ(n)` (the first two displays of thm:P2): statable, but outside the scope
  assigned to this target (Λ_m and the θ = π/2 reduction).
* Ahmed's exact formula itself as printed in Ahmed's preprint (not reproduced in the paper), and
  the full length formula `L_{π/2}(n)` at `θ = π/2`; only its ingredients (stmt_5..stmt_8) are
  stated.
* `Λ_ν` for complex `ν`, entireness in `ν`, the Legendre function `P_ν` for non-integer `ν`
  (not in Mathlib), the bound `|Λ_m| ≤ min(1-c, 2/(2m+1))`, and `Λ_2`.
-/
-- BEGIN STATEMENTS
namespace T13

open Polynomial Real

/-- Legendre polynomial `P_n` over `ℝ`, built from Mathlib's shifted Legendre polynomial:
`P_n(t) = shiftedLegendre n ((1 - t)/2)`. -/
noncomputable def P (n : ℕ) : ℝ[X] :=
  ((shiftedLegendre n).map (Int.castRingHom ℝ)).comp (C (1 / 2 : ℝ) * (1 - X))

/-- The paper's Legendre integral at integer order: `Λ_m(θ) = ∫_{cos θ}^1 P_m(t) dt`. -/
noncomputable def LamInt (m : ℕ) (θ : ℝ) : ℝ :=
  ∫ t in Real.cos θ..1, (P m).eval t

/-- `k`-th term of `2z ₂F₁(-ν, 1+ν; 2; z)`:
`2z · (-ν)_k (1+ν)_k / (k! (2)_k) · z^k`. -/
noncomputable def lamTerm (ν z : ℝ) (k : ℕ) : ℝ :=
  2 * z * ((ascPochhammer ℝ k).eval (-ν) * (ascPochhammer ℝ k).eval (1 + ν) /
    ((Nat.factorial k : ℝ) * (ascPochhammer ℝ k).eval 2)) * z ^ k

/-- `P` is the classical Legendre family (initial values and Bonnet's recurrence).
Sanity check of the definition, not a sentence of the paper. -/
def stmt_1 : Prop :=
  P 0 = 1 ∧ P 1 = X ∧
  ∀ n : ℕ, ((n : ℝ) + 2) • P (n + 2) = (2 * (n : ℝ) + 3) • (X * P (n + 1)) - ((n : ℝ) + 1) • P n

/-- thm:P2: "$\Lambda_m(\theta)=\frac{P_{m-1}(\cos\theta)-P_{m+1}(\cos\theta)}{2m+1}$",
with `Λ_m` the Legendre integral `∫_{cos θ}^1 P_m`. -/
def stmt_2 : Prop :=
  ∀ m : ℕ, 1 ≤ m → ∀ θ : ℝ,
    LamInt m θ = ((P (m - 1)).eval (Real.cos θ) - (P (m + 1)).eval (Real.cos θ)) / (2 * m + 1)

/-- thm:P2: "[$\Lambda_m(\theta)$] is a polynomial in $\cos\theta$" (with rational coefficients). -/
def stmt_3 : Prop :=
  ∀ m : ℕ, 1 ≤ m → ∃ Q : ℚ[X], ∀ θ : ℝ,
    LamInt m θ = (Q.map (Rat.castHom ℝ)).eval (Real.cos θ)

/-- Display before thm:P2 (and g3:lem:Lam(a)) at integer order:
`∫_{cos θ}^1 P_m = 2 sin²(θ/2) ₂F₁(-m, 1+m; 2; sin²(θ/2))`, `0 ≤ θ < π`. -/
def stmt_4 : Prop :=
  ∀ m : ℕ, ∀ θ : ℝ, 0 ≤ θ → θ < π →
    HasSum (fun k : ℕ => lamTerm (m : ℝ) (Real.sin (θ / 2) ^ 2) k) (LamInt m θ)

/-- Legendre values at `0` (used in g3:lem:pi2 and its last line): `P_{2k+1}(0) = 0` and
`P_{2k}(0) = (-1)^k C(2k,k) / 4^k`. -/
def stmt_5 : Prop :=
  ∀ k : ℕ, (P (2 * k + 1)).eval 0 = 0 ∧
    (P (2 * k)).eval 0 = (-1 : ℝ) ^ k * ((2 * k).choose k : ℝ) / 4 ^ k

/-- g3:lem:pi2: "$\Lambda_m(\frac\pi2)=0$ for even $m$" (for `m ≥ 2`; `Λ_m` is defined for `m ≥ 1`
in thm:P2). -/
def stmt_6 : Prop :=
  ∀ m : ℕ, 2 ≤ m → Even m → LamInt m (π / 2) = 0

/-- Explicit odd values at `θ = π/2`, from "$\Lambda_m(\frac\pi2)=\frac{P_{m-1}(0)-P_{m+1}(0)}{2m+1}$":
`Λ_{2k+1}(π/2) = (-1)^k C(2k,k) / (4^k (2k+2))`. -/
def stmt_7 : Prop :=
  ∀ k : ℕ, LamInt (2 * k + 1) (π / 2) =
    (-1 : ℝ) ^ k * ((2 * k).choose k : ℝ) / (4 ^ k * (2 * k + 2))

/-- g3:lem:pi2 (the `ε`-term of thm:P2 at `θ = π/2`, Ahmed's Gamma factor):
`(πε/sin πε) Λ_ε(π/2) = Γ(1+ε/2) Γ((1-ε)/2) / (√π (1+ε))`, with `Λ_ε` the hypergeometric series,
for real `ε > 0`, `ε ∉ ℕ`. -/
def stmt_8 : Prop :=
  ∀ ε : ℝ, 0 < ε → (∀ j : ℕ, ε ≠ j) →
    HasSum (fun k : ℕ => π * ε / Real.sin (π * ε) * lamTerm ε (Real.sin (π / 2 / 2) ^ 2) k)
      (Real.Gamma (1 + ε / 2) * Real.Gamma ((1 - ε) / 2) / (Real.sqrt π * (1 + ε)))

/-- main.tex 1458: "$\Lambda_1=\frac12\sin^2\theta$". -/
def stmt_9 : Prop :=
  ∀ θ : ℝ, LamInt 1 θ = Real.sin θ ^ 2 / 2

end T13
-- END STATEMENTS

namespace T13
open Polynomial Real

noncomputable def bc (n k : ℕ) : ℝ :=
  (ascPochhammer ℝ k).eval (-(n:ℝ)) * (ascPochhammer ℝ k).eval (1 + (n:ℝ)) / ((Nat.factorial k : ℝ)^2)

lemma bc_eq (n k : ℕ) :
    bc n k = (-1:ℝ)^k * (n.choose k : ℝ) * ((n+k).choose n : ℝ) := by
  have h1 : (ascPochhammer ℝ k).eval (-(n:ℝ)) = (-1)^k * (n.descFactorial k : ℝ) := by
    rw [ascPochhammer_eval_neg_eq_descPochhammer, descPochhammer_eval_eq_descFactorial]
  have h2 : (ascPochhammer ℝ k).eval (1 + (n:ℝ)) = ((n+1).ascFactorial k : ℝ) := by
    have := ascPochhammer_eval_cast (S := ℝ) k (n+1)
    rw [ascPochhammer_nat_eq_ascFactorial] at this
    push_cast at this
    rw [add_comm]; exact this.symm
  rw [bc, h1, h2, Nat.descFactorial_eq_factorial_mul_choose, Nat.ascFactorial_eq_factorial_mul_choose,
    ← Nat.choose_symm_add]
  have : (k.factorial : ℝ) ≠ 0 := by positivity
  push_cast
  field_simp

lemma bc_zero (n k : ℕ) (h : n < k) : bc n k = 0 := by
  rw [bc, ascPochhammer_eval_neg_natCast_of_lt h]; simp

lemma P_eval (n : ℕ) (t : ℝ) :
    (P n).eval t = ∑ k ∈ Finset.range (n+1), bc n k * ((1-t)/2)^k := by
  simp only [P, shiftedLegendre, eval_comp, eval_map]
  rw [eval₂_finsetSum]
  refine Finset.sum_congr rfl fun k _ => ?_
  rw [bc_eq]
  rw [eval₂_mul, eval₂_C, eval₂_X_pow]
  simp only [eval_mul, eval_C, eval_sub, eval_one, eval_X, eq_intCast]
  push_cast
  ring

lemma pe_succ_left' (x : ℝ) (m : ℕ) :
    (ascPochhammer ℝ (m+1)).eval x = x * (ascPochhammer ℝ m).eval (x+1) := by
  rw [ascPochhammer_succ_left]; simp [eval_comp]

lemma pe_succ_right' (x : ℝ) (m : ℕ) :
    (ascPochhammer ℝ (m+1)).eval x = (ascPochhammer ℝ m).eval x * (x + m) :=
  ascPochhammer_succ_eval m x

lemma bc_zero_idx (n : ℕ) : bc n 0 = 1 := by simp [bc]

lemma bc_one_idx (n : ℕ) : bc n 1 = -((n:ℝ) * ((n:ℝ) + 1)) := by
  simp only [bc, ascPochhammer_one, eval_X, Nat.factorial_one]; push_cast; ring

lemma bc_f2 (n i : ℕ) : bc (n+2) (i+2) =
    ((n:ℝ)+2)*((n:ℝ)+1) * (ascPochhammer ℝ i).eval (-(n:ℝ)) *
      ((ascPochhammer ℝ i).eval ((n:ℝ)+3) * (((n:ℝ)+3+i)*((n:ℝ)+4+i))) /
      (((i+2).factorial : ℝ)^2) := by
  have e1 : (ascPochhammer ℝ (i+2)).eval (-((n+2:ℕ):ℝ)) =
      ((n:ℝ)+2)*((n:ℝ)+1) * (ascPochhammer ℝ i).eval (-(n:ℝ)) := by
    rw [pe_succ_left', pe_succ_left']
    have h1 : -((n+2:ℕ):ℝ) + 1 = -((n:ℝ)+1) := by push_cast; ring
    have h2 : -((n:ℝ)+1) + 1 = -(n:ℝ) := by ring
    rw [h1, h2]; push_cast; ring
  have e2 : (ascPochhammer ℝ (i+2)).eval (1 + ((n+2:ℕ):ℝ)) =
      (ascPochhammer ℝ i).eval ((n:ℝ)+3) * (((n:ℝ)+3+i)*((n:ℝ)+4+i)) := by
    have h : 1 + ((n+2:ℕ):ℝ) = (n:ℝ)+3 := by push_cast; ring
    rw [h, pe_succ_right', pe_succ_right']
    push_cast; ring
  rw [bc, e1, e2]

lemma bc_f1 (n i : ℕ) : bc (n+1) (i+2) =
    -((n:ℝ)+1) * (ascPochhammer ℝ i).eval (-(n:ℝ)) * (-(n:ℝ)+i) *
      (((n:ℝ)+2) * (ascPochhammer ℝ i).eval ((n:ℝ)+3) * ((n:ℝ)+3+i)) /
      (((i+2).factorial : ℝ)^2) := by
  have e1 : (ascPochhammer ℝ (i+2)).eval (-((n+1:ℕ):ℝ)) =
      -((n:ℝ)+1) * (ascPochhammer ℝ i).eval (-(n:ℝ)) * (-(n:ℝ)+i) := by
    rw [pe_succ_left']
    have h1 : -((n+1:ℕ):ℝ) + 1 = -(n:ℝ) := by push_cast; ring
    rw [h1, pe_succ_right']
    push_cast; ring
  have e2 : (ascPochhammer ℝ (i+2)).eval (1 + ((n+1:ℕ):ℝ)) =
      ((n:ℝ)+2) * (ascPochhammer ℝ i).eval ((n:ℝ)+3) * ((n:ℝ)+3+i) := by
    rw [pe_succ_left']
    have h1 : 1 + ((n+1:ℕ):ℝ) = (n:ℝ)+2 := by push_cast; ring
    have h2 : (n:ℝ)+2+1 = (n:ℝ)+3 := by ring
    rw [h1, h2, pe_succ_right']
    push_cast; ring
  rw [bc, e1, e2]

lemma bc_f1' (n i : ℕ) : bc (n+1) (i+1) =
    -((n:ℝ)+1) * (ascPochhammer ℝ i).eval (-(n:ℝ)) *
      (((n:ℝ)+2) * (ascPochhammer ℝ i).eval ((n:ℝ)+3)) /
      (((i+1).factorial : ℝ)^2) := by
  have e1 : (ascPochhammer ℝ (i+1)).eval (-((n+1:ℕ):ℝ)) =
      -((n:ℝ)+1) * (ascPochhammer ℝ i).eval (-(n:ℝ)) := by
    rw [pe_succ_left']
    have h1 : -((n+1:ℕ):ℝ) + 1 = -(n:ℝ) := by push_cast; ring
    rw [h1]; push_cast; ring
  have e2 : (ascPochhammer ℝ (i+1)).eval (1 + ((n+1:ℕ):ℝ)) =
      ((n:ℝ)+2) * (ascPochhammer ℝ i).eval ((n:ℝ)+3) := by
    rw [pe_succ_left']
    have h1 : 1 + ((n+1:ℕ):ℝ) = (n:ℝ)+2 := by push_cast; ring
    have h2 : (n:ℝ)+2+1 = (n:ℝ)+3 := by ring
    rw [h1, h2]
  rw [bc, e1, e2]

lemma bc_f0 (n i : ℕ) : bc n (i+2) =
    (ascPochhammer ℝ i).eval (-(n:ℝ)) * ((-(n:ℝ)+i) * (-(n:ℝ)+i+1)) *
      (((n:ℝ)+1)*((n:ℝ)+2) * (ascPochhammer ℝ i).eval ((n:ℝ)+3)) /
      (((i+2).factorial : ℝ)^2) := by
  have e1 : (ascPochhammer ℝ (i+2)).eval (-(n:ℝ)) =
      (ascPochhammer ℝ i).eval (-(n:ℝ)) * ((-(n:ℝ)+i) * (-(n:ℝ)+i+1)) := by
    rw [pe_succ_right', pe_succ_right']
    push_cast; ring
  have e2 : (ascPochhammer ℝ (i+2)).eval (1 + (n:ℝ)) =
      ((n:ℝ)+1)*((n:ℝ)+2) * (ascPochhammer ℝ i).eval ((n:ℝ)+3) := by
    rw [pe_succ_left', pe_succ_left']
    have h1 : 1 + (n:ℝ) + 1 = (n:ℝ)+2 := by ring
    have h2 : (n:ℝ)+2+1 = (n:ℝ)+3 := by ring
    rw [h1, h2]; ring
  rw [bc, e1, e2]

lemma I1_main (n i : ℕ) :
    ((n:ℝ)+2) * bc (n+2) (i+2) =
      (2*(n:ℝ)+3) * (bc (n+1) (i+2) - 2 * bc (n+1) (i+1)) - ((n:ℝ)+1) * bc n (i+2) := by
  rw [bc_f2, bc_f1, bc_f1', bc_f0]
  have hf : ((i.factorial : ℕ) : ℝ) ≠ 0 := by positivity
  simp only [Nat.factorial_succ]
  push_cast
  field_simp
  ring

lemma I1_zero (n : ℕ) :
    ((n:ℝ)+2) * bc (n+2) 1 =
      (2*(n:ℝ)+3) * (bc (n+1) 1 - 2 * bc (n+1) 0) - ((n:ℝ)+1) * bc n 1 := by
  simp only [bc_one_idx, bc_zero_idx]; push_cast; ring

lemma I2_main (n i : ℕ) :
    2 * (2*(n:ℝ)+3) * bc (n+1) (i+1) / (((i+1 : ℕ):ℝ)+1) =
      bc n (i+2) - bc (n+2) (i+2) := by
  rw [bc_f2, bc_f1', bc_f0]
  have hf : ((i.factorial : ℕ) : ℝ) ≠ 0 := by positivity
  simp only [Nat.factorial_succ]
  push_cast
  field_simp
  ring

lemma I2_zero (n : ℕ) :
    2 * (2*(n:ℝ)+3) * bc (n+1) 0 / (((0 : ℕ):ℝ)+1) = bc n 1 - bc (n+2) 1 := by
  simp only [bc_one_idx, bc_zero_idx]; push_cast; ring



noncomputable def B (n : ℕ) : ℝ[X] :=
  ∑ k ∈ Finset.range (n+1), C (bc n k) * X ^ k

lemma coeff_B (n k : ℕ) : (B n).coeff k = bc n k := by
  simp only [B, finsetSum_coeff, coeff_C_mul_X_pow]
  by_cases h : k ≤ n
  · rw [Finset.sum_eq_single k]
    · simp
    · intro b _ hb; simp [Ne.symm hb]
    · intro hk; exact absurd (Finset.mem_range.2 (Nat.lt_succ_of_le h)) hk
  · replace h := Nat.lt_of_not_ge h
    rw [bc_zero n k h]
    apply Finset.sum_eq_zero
    intro b hb
    have : b ≠ k := by
      intro e; subst e; exact absurd (Finset.mem_range.1 hb) (by omega)
    simp [this.symm]

lemma eval_B (n : ℕ) (u : ℝ) : (B n).eval u = ∑ k ∈ Finset.range (n+1), bc n k * u ^ k := by
  simp [B, eval_finsetSum]

lemma P_eval' (n : ℕ) (t : ℝ) : (P n).eval t = (B n).eval ((1 - t)/2) := by
  rw [P_eval, eval_B]

lemma P_eval_range (n N : ℕ) (h : n + 1 ≤ N) (t : ℝ) :
    (P n).eval t = ∑ k ∈ Finset.range N, bc n k * ((1 - t)/2) ^ k := by
  rw [P_eval]
  apply Finset.sum_subset (Finset.range_subset_range.2 h)
  intro k hk hk'
  have : n < k := by
    simp only [Finset.mem_range] at hk hk'; omega
  rw [bc_zero n k this]; simp

lemma bonnet_B (n : ℕ) :
    C ((n:ℝ)+2) * B (n+2) =
      C (2*(n:ℝ)+3) * (B (n+1) - C 2 * (X * B (n+1))) - C ((n:ℝ)+1) * B n := by
  ext k
  simp only [coeff_C_mul, coeff_sub, coeff_B]
  rcases k with _ | j
  · simp [bc_zero_idx]; ring
  · rw [coeff_X_mul]
    simp only [coeff_B]
    rcases j with _ | i
    · simpa using I1_zero n
    · exact I1_main n i

lemma I2_gen (n k : ℕ) :
    2 * (2*(n:ℝ)+3) * bc (n+1) k / ((k:ℝ)+1) = bc n (k+1) - bc (n+2) (k+1) := by
  rcases k with _ | i
  · simpa using I2_zero n
  · exact I2_main n i

lemma int_pow_term (k : ℕ) (c : ℝ) :
    ∫ t in c..1, ((1-t)/2)^k = 2/((k:ℝ)+1) * ((1-c)/2)^(k+1) := by
  have hd : ∀ x ∈ Set.uIcc c 1, HasDerivAt (fun t : ℝ => -(2/((k:ℝ)+1)) * ((1-t)/2)^(k+1))
      (((1-x)/2)^k) x := by
    intro x _
    have h1 : HasDerivAt (fun t : ℝ => (1-t)/2) (-1/2) x := by
      have := ((hasDerivAt_id x).const_sub 1).div_const 2
      simpa using this
    have h2 := (h1.pow (k+1)).const_mul (-(2/((k:ℝ)+1)))
    convert h2 using 1
    have : ((k:ℝ)+1) ≠ 0 := by positivity
    push_cast
    field_simp
  rw [intervalIntegral.integral_eq_sub_of_hasDerivAt hd
    (by apply Continuous.intervalIntegrable; fun_prop)]
  simp

lemma LamInt_sum (m : ℕ) (θ : ℝ) :
    LamInt m θ = ∑ k ∈ Finset.range (m+1),
      bc m k * (2/((k:ℝ)+1) * ((1 - Real.cos θ)/2)^(k+1)) := by
  unfold LamInt
  have : ∀ t : ℝ, (P m).eval t = ∑ k ∈ Finset.range (m+1), bc m k * ((1-t)/2)^k :=
    fun t => P_eval m t
  simp_rw [this]
  rw [intervalIntegral.integral_finsetSum]
  · refine Finset.sum_congr rfl fun k _ => ?_
    rw [intervalIntegral.integral_const_mul, int_pow_term]
  · intro k _
    apply Continuous.intervalIntegrable; fun_prop


lemma P_zero : P 0 = 1 := by
  apply Polynomial.funext; intro t
  rw [P_eval, Finset.sum_range_one, bc_zero_idx]; simp

lemma P_one : P 1 = X := by
  apply Polynomial.funext; intro t
  rw [P_eval, Finset.sum_range_succ, Finset.sum_range_one, bc_zero_idx, bc_one_idx]
  simp; ring

lemma P_bonnet (n : ℕ) :
    ((n : ℝ) + 2) • P (n + 2) = (2 * (n : ℝ) + 3) • (X * P (n + 1)) - ((n : ℝ) + 1) • P n := by
  apply Polynomial.funext; intro r
  have h := congrArg (eval ((1 - r)/2)) (bonnet_B n)
  simp only [eval_mul, eval_C, eval_sub, eval_X] at h
  simp only [eval_smul, eval_mul, eval_sub, eval_X, smul_eq_mul, P_eval']
  linear_combination h

lemma Q_sum (n : ℕ) (θ : ℝ) :
    (P n).eval (Real.cos θ) = ∑ k ∈ Finset.range (n+3), bc n k * ((1 - Real.cos θ)/2) ^ k :=
  P_eval_range n (n+3) (by omega) _

lemma lam_succ (n : ℕ) (θ : ℝ) :
    LamInt (n+1) θ = ((P n).eval (Real.cos θ) - (P (n+2)).eval (Real.cos θ)) /
      (2 * ((n+1 : ℕ) : ℝ) + 1) := by
  set z := (1 - Real.cos θ)/2 with hz
  have h1 : (P n).eval (Real.cos θ) = ∑ k ∈ Finset.range (n+3), bc n k * z ^ k :=
    P_eval_range n (n+3) (by omega) _
  have h2 : (P (n+2)).eval (Real.cos θ) = ∑ k ∈ Finset.range (n+3), bc (n+2) k * z ^ k :=
    P_eval_range (n+2) (n+3) (by omega) _
  rw [h1, h2, ← Finset.sum_sub_distrib, Finset.sum_range_succ']
  rw [LamInt_sum, eq_div_iff (by positivity)]
  simp only [bc_zero_idx, pow_zero, sub_self, add_zero]
  rw [Finset.sum_mul]
  refine Finset.sum_congr rfl fun k _ => ?_
  have := I2_gen n k
  have hk : ((k:ℝ)+1) ≠ 0 := by positivity
  rw [← sub_mul, ← this]
  push_cast
  field_simp
  ring
lemma sin_half_sq (θ : ℝ) : Real.sin (θ/2)^2 = (1 - Real.cos θ)/2 := by
  have h := Real.cos_two_mul (θ/2)
  rw [show 2*(θ/2) = θ by ring] at h
  nlinarith [Real.sin_sq_add_cos_sq (θ/2)]

lemma pe_two (k : ℕ) : (ascPochhammer ℝ k).eval 2 = (((k+1).factorial : ℕ) : ℝ) := by
  induction k with
  | zero => simp
  | succ k ih =>
    rw [pe_succ_right', ih, Nat.factorial_succ (k+1)]
    push_cast; ring

lemma lamTerm_eq (m k : ℕ) (z : ℝ) :
    lamTerm (m:ℝ) z k = bc m k * (2/((k:ℝ)+1) * z^(k+1)) := by
  unfold lamTerm bc
  rw [pe_two, Nat.factorial_succ k]
  have hf : (k.factorial : ℝ) ≠ 0 := by positivity
  push_cast
  field_simp
  ring

lemma P_rec0 (n : ℕ) :
    ((n:ℝ)+2) * (P (n+2)).eval 0 = -(((n:ℝ)+1) * (P n).eval 0) := by
  have h := congrArg (eval 0) (P_bonnet n)
  simp only [eval_smul, eval_mul, eval_sub, eval_X, smul_eq_mul] at h
  linear_combination h

lemma P_at0 (k : ℕ) : (P (2 * k + 1)).eval 0 = 0 ∧
    (P (2 * k)).eval 0 = (-1 : ℝ) ^ k * ((2 * k).choose k : ℝ) / 4 ^ k := by
  induction k with
  | zero => simp [P_zero, P_one]
  | succ k ih =>
    obtain ⟨ih1, ih2⟩ := ih
    have r1 := P_rec0 (2*k+1)
    have r2 := P_rec0 (2*k)
    have e1 : 2 * (k+1) + 1 = (2*k+1) + 2 := by ring
    have e2 : 2 * (k+1) = (2*k) + 2 := by ring
    rw [e1, e2]
    constructor
    · rw [ih1] at r1
      have : ((2*k+1 : ℕ):ℝ) + 2 ≠ 0 := by positivity
      have h0 : (((2*k+1 : ℕ):ℝ) + 2) * (P (2*k+1+2)).eval 0 = 0 := by rw [r1]; ring
      exact (mul_eq_zero.mp h0).resolve_left this
    · have hcb := Nat.succ_mul_centralBinom_succ k
      rw [Nat.centralBinom_eq_two_mul_choose, Nat.centralBinom_eq_two_mul_choose] at hcb
      have hcb' : ((k:ℝ)+1) * ((2*(k+1)).choose (k+1) : ℝ) = 2 * (2*k+1) * ((2*k).choose k : ℝ) := by
        exact_mod_cast hcb
      push_cast at r2
      rw [ih2] at r2
      have hk : ((2*k:ℝ)+2) ≠ 0 := by positivity
      have e3 : (2*k+2 : ℕ) = 2*(k+1) := by ring
      have e4 : ((2*k+2).choose (k+1) : ℝ) = ((2*(k+1)).choose (k+1) : ℝ) := by rw [e3]
      rw [e4]
      have hk1 : ((k:ℝ)+1) ≠ 0 := by positivity
      have hE : (P (2*k+2)).eval 0 =
          (-((2*k+1:ℝ)) * ((-1)^k * ((2*k).choose k : ℝ) / 4^k)) / (2*k+2) := by
        rw [eq_div_iff hk]; linear_combination r2
      have hC : ((2*(k+1)).choose (k+1) : ℝ) = 2*(2*k+1) * ((2*k).choose k : ℝ) / (k+1) := by
        field_simp; linarith [hcb']
      rw [hE, hC]
      field_simp
      ring

noncomputable def Pq (n : ℕ) : ℚ[X] :=
  ((shiftedLegendre n).map (Int.castRingHom ℚ)).comp (C (1 / 2 : ℚ) * (1 - X))

lemma Pq_map (n : ℕ) : (Pq n).map (Rat.castHom ℝ) = P n := by
  have h : (Rat.castHom ℝ).comp (Int.castRingHom ℚ) = Int.castRingHom ℝ := RingHom.ext_int _ _
  simp only [Pq, P, Polynomial.map_comp, Polynomial.map_map, h, Polynomial.map_mul,
    Polynomial.map_sub, Polynomial.map_one, Polynomial.map_X, map_C]
  simp

end T13

theorem T13.proof_1 : T13.stmt_1 := by
  refine ⟨T13.P_zero, T13.P_one, fun n => T13.P_bonnet n⟩

theorem T13.proof_2 : T13.stmt_2 := by
  intro m hm θ
  obtain ⟨n, rfl⟩ : ∃ n, m = n + 1 := ⟨m - 1, by omega⟩
  simpa using T13.lam_succ n θ

theorem T13.proof_3 : T13.stmt_3 := by
  intro m hm
  refine ⟨Polynomial.C (1 / (2 * (m:ℚ) + 1)) * (T13.Pq (m-1) - T13.Pq (m+1)), fun θ => ?_⟩
  rw [T13.proof_2 m hm θ]
  have h1 := congrArg (Polynomial.eval (Real.cos θ)) (T13.Pq_map (m-1))
  have h2 := congrArg (Polynomial.eval (Real.cos θ)) (T13.Pq_map (m+1))
  simp only [Polynomial.map_mul, Polynomial.map_sub, Polynomial.map_C, Polynomial.eval_mul,
    Polynomial.eval_C, Polynomial.eval_sub, h1, h2]
  simp
  ring

theorem T13.proof_4 : T13.stmt_4 := by
  intro m θ _ _
  have hz : Real.sin (θ/2)^2 = (1 - Real.cos θ)/2 := T13.sin_half_sq θ
  rw [T13.LamInt_sum, hz]
  have : HasSum (fun k : ℕ => T13.lamTerm (m:ℝ) ((1 - Real.cos θ)/2) k)
      (∑ k ∈ Finset.range (m+1), T13.lamTerm (m:ℝ) ((1 - Real.cos θ)/2) k) :=
    hasSum_sum_of_ne_finset_zero (by
      intro k hk
      rw [T13.lamTerm_eq, T13.bc_zero m k (by simpa using hk)]; simp)
  simpa [T13.lamTerm_eq] using this

theorem T13.proof_5 : T13.stmt_5 := fun k => T13.P_at0 k

theorem T13.proof_6 : T13.stmt_6 := by
  intro m hm hev
  obtain ⟨j, rfl⟩ := hev
  obtain ⟨i, rfl⟩ : ∃ i, j = i + 1 := ⟨j - 1, by omega⟩
  rw [T13.proof_2 (i + 1 + (i + 1)) (by omega), Real.cos_pi_div_two]
  have e1 : i + 1 + (i + 1) - 1 = 2 * i + 1 := by omega
  have e2 : i + 1 + (i + 1) + 1 = 2 * (i + 1) + 1 := by omega
  rw [e1, e2, (T13.P_at0 i).1, (T13.P_at0 (i+1)).1]
  simp

theorem T13.proof_7 : T13.stmt_7 := by
  intro k
  rw [T13.proof_2 (2 * k + 1) (by omega), Real.cos_pi_div_two]
  have e1 : 2 * k + 1 - 1 = 2 * k := by omega
  rw [e1, (T13.P_at0 k).2]
  have r := T13.P_rec0 (2 * k)
  rw [(T13.P_at0 k).2] at r
  have hk : ((2 * k : ℕ) : ℝ) + 2 ≠ 0 := by positivity
  have e2 : 2 * k + 1 + 1 = 2 * k + 2 := by omega
  rw [e2]
  have hE : (T13.P (2*k+2)).eval 0 =
      (-(((2*k : ℕ):ℝ) + 1) * ((-1:ℝ)^k * ((2*k).choose k : ℝ) / 4^k)) / (((2*k:ℕ):ℝ) + 2) := by
    rw [eq_div_iff hk]; linear_combination r
  rw [hE]
  push_cast
  have : (2 * (k:ℝ) + 2) ≠ 0 := by positivity
  field_simp
  ring

theorem T13.proof_9 : T13.stmt_9 := by
  intro θ
  unfold T13.LamInt
  rw [T13.P_one]
  simp only [Polynomial.eval_X]
  rw [integral_id]
  nlinarith [Real.sin_sq_add_cos_sq θ]

namespace T13
open Polynomial Real Filter Topology

/-! Helpers for `proof_8`: Bailey's sum `₂F₁(-ε,1+ε;2;1/2)` via the contiguous relation
`c(c+1)F(c) = (c+1+ε)(c-ε)F(c+2)` at `z = 1/2`, iterated, with Euler's limit formula. -/
/-- `k`-th term of `₂F₁(-ε, 1+ε; c; 1/2)`. -/
noncomputable def hT (ε c : ℝ) (k : ℕ) : ℝ :=
  (ascPochhammer ℝ k).eval (-ε) * (ascPochhammer ℝ k).eval (1 + ε) /
    ((Nat.factorial k : ℝ) * (ascPochhammer ℝ k).eval c) * (1/2 : ℝ) ^ k

lemma hT_zero (ε c : ℝ) : hT ε c 0 = 1 := by simp [hT]

lemma hT_succ (ε c : ℝ) (hc : 0 < c) (k : ℕ) :
    hT ε c (k+1) = hT ε c k * ((-ε + k) * (1 + ε + k) / (((k:ℝ) + 1) * (c + k)) * (1/2)) := by
  have hp := ascPochhammer_pos k c hc
  have hf : (k.factorial : ℝ) ≠ 0 := by positivity
  have hck : c + k ≠ 0 := by positivity
  unfold hT
  rw [pe_succ_right', pe_succ_right', pe_succ_right', Nat.factorial_succ]
  push_cast
  field_simp
  ring

lemma hT_shift2 (ε c : ℝ) (hc : 0 < c) (k : ℕ) :
    hT ε (c+2) k = hT ε c k * (c * (c+1) / ((c + k) * (c + k + 1))) := by
  have hp := ascPochhammer_pos k c hc
  have hp2 := ascPochhammer_pos k (c+2) (by linarith)
  have key : (ascPochhammer ℝ k).eval (c+2) * (c * (c+1)) =
      (ascPochhammer ℝ k).eval c * ((c + k) * (c + k + 1)) := by
    have e1 := pe_succ_left' c (k+1)
    rw [pe_succ_left' (c+1) k] at e1
    have e2 := pe_succ_right' c (k+1)
    rw [pe_succ_right' c k] at e2
    have h3 : c + 1 + 1 = c + 2 := by ring
    rw [h3] at e1
    rw [e1] at e2
    push_cast at e2
    linear_combination e2
  have hf : (k.factorial : ℝ) ≠ 0 := by positivity
  have hck : c + k ≠ 0 := by positivity
  have hck1 : c + k + 1 ≠ 0 := by positivity
  have e : (ascPochhammer ℝ k).eval (c+2) =
      (ascPochhammer ℝ k).eval c * ((c + k) * (c + k + 1)) / (c * (c+1)) := by
    rw [eq_div_iff (by positivity)]; exact key
  unfold hT
  rw [e]
  field_simp

lemma hT_summable (ε c : ℝ) (hε : 0 < ε) (hc : 1 ≤ c) : Summable (hT ε c) := by
  refine summable_of_ratio_norm_eventually_le (r := 1/2) (by norm_num) ?_
  refine eventually_atTop.2 ⟨⌈ε⌉₊, fun k hk => ?_⟩
  have hk' : ε ≤ (k:ℝ) := Nat.ceil_le.1 hk
  rw [hT_succ ε c (by linarith) k, norm_mul, mul_comm]
  apply mul_le_mul_of_nonneg_right _ (norm_nonneg _)
  have h1 : 0 ≤ (-ε + k) * (1 + ε + k) := by apply mul_nonneg <;> linarith
  have h2 : 0 < ((k:ℝ) + 1) * (c + k) := by positivity
  rw [Real.norm_eq_abs, abs_of_nonneg (by positivity)]
  have : (-ε + k) * (1 + ε + k) / (((k:ℝ) + 1) * (c + k)) ≤ 1 := by
    rw [div_le_one h2]; nlinarith
  linarith

/-- The contiguous relation `c(c+1) F(c) = (c+1+ε)(c-ε) F(c+2)` at `z = 1/2`. -/
lemma hT_contig (ε c : ℝ) (hε : 0 < ε) (hc : 1 ≤ c) :
    c * (c+1) * ∑' k, hT ε c k = (c + 1 + ε) * (c - ε) * ∑' k, hT ε (c+2) k := by
  have hc0 : 0 < c := by linarith
  set R : ℕ → ℝ := fun k => -2 * c * (c+1) * k / (c + k) * hT ε c k with hR
  have tele : ∀ k : ℕ, c * (c+1) * hT ε c k - (c + 1 + ε) * (c - ε) * hT ε (c+2) k =
      R (k+1) - R k := by
    intro k
    simp only [hR]
    rw [hT_succ ε c hc0 k, hT_shift2 ε c hc0 k]
    have hck : c + k ≠ 0 := by positivity
    have hck1 : c + k + 1 ≠ 0 := by positivity
    have hk1 : (k:ℝ) + 1 ≠ 0 := by positivity
    push_cast
    have : c + ((k:ℝ) + 1) = c + k + 1 := by ring
    rw [this]
    field_simp
    ring
  have hs1 := (hT_summable ε c hε hc).hasSum
  have hs2 := (hT_summable ε (c+2) hε (by linarith)).hasSum
  have hs := (hs1.mul_left (c * (c+1))).sub (hs2.mul_left ((c + 1 + ε) * (c - ε)))
  have hlim := hs.tendsto_sum_nat
  simp_rw [tele, Finset.sum_range_sub] at hlim
  have hR0 : Tendsto R atTop (𝓝 0) := by
    have ht := (hT_summable ε c hε hc).tendsto_atTop_zero
    have hb : Tendsto (fun k : ℕ => 2 * c * (c+1) * |hT ε c k|) atTop (𝓝 0) := by
      simpa using (ht.abs.const_mul (2 * c * (c+1)))
    refine squeeze_zero_norm (fun k => ?_) hb
    simp only [hR, Real.norm_eq_abs]
    rw [abs_mul]
    have hck : 0 < c + k := by positivity
    have : |-2 * c * (c+1) * (k:ℝ) / (c + k)| ≤ 2 * c * (c+1) := by
      rw [abs_div, abs_of_pos hck, div_le_iff₀ hck]
      rw [show -2 * c * (c+1) * (k:ℝ) = -(2 * c * (c+1) * k) by ring, abs_neg,
        abs_of_nonneg (by positivity)]
      nlinarith [sq_nonneg c]
    exact mul_le_mul_of_nonneg_right this (abs_nonneg _)
  have h0 : R 0 = 0 := by simp [hR]
  rw [h0] at hlim
  simp only [sub_zero] at hlim
  have := tendsto_nhds_unique hlim hR0
  linarith

lemma poch_cmp (c : ℝ) (hc : 2 ≤ c) (k : ℕ) :
    c * (ascPochhammer ℝ (k+1)).eval 2 ≤ 2 * (ascPochhammer ℝ (k+1)).eval c := by
  induction k with
  | zero => simp; linarith
  | succ k ih =>
    rw [pe_succ_right' 2 (k+1), pe_succ_right' c (k+1)]
    have h2 := ascPochhammer_pos (k+1) (2:ℝ) (by norm_num)
    have hc' := ascPochhammer_pos (k+1) c (by linarith)
    push_cast
    nlinarith

lemma hT_bound (ε c : ℝ) (hc : 2 ≤ c) (k : ℕ) :
    |hT ε c (k+1)| ≤ 2 / c * |hT ε 2 (k+1)| := by
  have h2 := ascPochhammer_pos (k+1) (2:ℝ) (by norm_num)
  have hc' := ascPochhammer_pos (k+1) c (by linarith)
  have hf : ((k+1).factorial : ℝ) ≠ 0 := by positivity
  have e : hT ε c (k+1) = hT ε 2 (k+1) *
      ((ascPochhammer ℝ (k+1)).eval 2 / (ascPochhammer ℝ (k+1)).eval c) := by
    unfold hT; field_simp
  rw [e, abs_mul, abs_of_pos (div_pos h2 hc'), mul_comm]
  apply mul_le_mul_of_nonneg_right _ (abs_nonneg _)
  rw [div_le_div_iff₀ hc' (by linarith)]
  linarith [poch_cmp c hc k]

lemma S_near_one (ε c : ℝ) (hε : 0 < ε) (hc : 2 ≤ c) :
    |∑' k, hT ε c k - 1| ≤ 2 / c * ∑' k, |hT ε 2 (k+1)| := by
  rw [(hT_summable ε c hε (by linarith)).tsum_eq_zero_add, hT_zero, add_sub_cancel_left]
  have hs : Summable (fun k => |hT ε 2 (k+1)|) :=
    (summable_nat_add_iff (f := fun k => |hT ε 2 k|) 1).mpr (hT_summable ε 2 hε (by norm_num)).abs
  rw [← Real.norm_eq_abs]
  refine tsum_of_norm_bounded (f := fun k => hT ε c (k+1)) (hs.hasSum.mul_left (2 / c)) (fun k => ?_)
  rw [Real.norm_eq_abs]
  exact hT_bound ε c hc k

lemma S_tendsto (ε : ℝ) (hε : 0 < ε) :
    Tendsto (fun N : ℕ => ∑' k, hT ε (2 + 2 * (N:ℝ)) k) atTop (𝓝 1) := by
  set A := ∑' k, |hT ε 2 (k+1)|
  rw [tendsto_iff_norm_sub_tendsto_zero]
  have hb : Tendsto (fun N : ℕ => 1 / ((N:ℝ) + 1) * A) atTop (𝓝 0) := by
    simpa using (tendsto_one_div_add_atTop_nhds_zero_nat (𝕜 := ℝ)).mul_const A
  refine squeeze_zero_norm (fun N => ?_) hb
  rw [norm_norm, Real.norm_eq_abs]
  have := S_near_one ε (2 + 2 * (N:ℝ)) hε (by have := N.cast_nonneg (α := ℝ); linarith)
  have e : 2 / (2 + 2 * (N:ℝ)) = 1 / ((N:ℝ) + 1) := by
    field_simp; ring
  rwa [e] at this

/-- Factor of the iterated contiguous relation. -/
noncomputable def fct (ε : ℝ) (j : ℕ) : ℝ :=
  ((3 + ε) / 2 + j) * (1 - ε / 2 + j) / ((1 + j) * (3 / 2 + j))

lemma S_iter (ε : ℝ) (hε : 0 < ε) (N : ℕ) :
    ∑' k, hT ε 2 k = (∏ j ∈ Finset.range N, fct ε j) * ∑' k, hT ε (2 + 2 * (N:ℝ)) k := by
  induction N with
  | zero => simp
  | succ N ih =>
    rw [ih, Finset.prod_range_succ]
    have hN := N.cast_nonneg (α := ℝ)
    have hcon := hT_contig ε (2 + 2 * (N:ℝ)) hε (by linarith)
    have e : (2:ℝ) + 2 * ((N+1 : ℕ):ℝ) = 2 + 2 * (N:ℝ) + 2 := by push_cast; ring
    rw [e]
    have h0 : (2 + 2 * (N:ℝ)) * (2 + 2 * (N:ℝ) + 1) ≠ 0 := by positivity
    have hS : ∑' k, hT ε (2 + 2 * (N:ℝ)) k =
        (2 + 2 * (N:ℝ) + 1 + ε) * (2 + 2 * (N:ℝ) - ε) / ((2 + 2 * (N:ℝ)) * (2 + 2 * (N:ℝ) + 1)) *
          ∑' k, hT ε (2 + 2 * (N:ℝ) + 2) k := by
      rw [div_mul_eq_mul_div, eq_div_iff h0]; linear_combination hcon
    rw [hS]
    unfold fct
    have h1 : (1 + (N:ℝ)) ≠ 0 := by positivity
    have h2 : (3 / 2 + (N:ℝ)) ≠ 0 := by positivity
    field_simp
    ring

lemma prod_fct_eq (ε : ℝ) (hε : 0 < ε) (hnat : ∀ j : ℕ, ε ≠ j) (N : ℕ) (hN : 1 ≤ N) :
    ∏ j ∈ Finset.range (N+1), fct ε j =
      Real.GammaSeq 1 N * Real.GammaSeq (3/2) N /
        (Real.GammaSeq ((3 + ε) / 2) N * Real.GammaSeq (1 - ε / 2) N) := by
  have hN0 : (0:ℝ) < N := by exact_mod_cast hN
  have hPa : ∏ j ∈ Finset.range (N+1), ((3 + ε) / 2 + (j:ℝ)) ≠ 0 :=
    Finset.prod_ne_zero_iff.2 fun j _ => by positivity
  have hPb : ∏ j ∈ Finset.range (N+1), (1 - ε / 2 + (j:ℝ)) ≠ 0 :=
    Finset.prod_ne_zero_iff.2 fun j _ => by
      intro h; apply hnat (2 + 2 * j); push_cast; linarith
  have hP1 : ∏ j ∈ Finset.range (N+1), (1 + (j:ℝ)) ≠ 0 :=
    Finset.prod_ne_zero_iff.2 fun j _ => by positivity
  have hP32 : ∏ j ∈ Finset.range (N+1), (3 / 2 + (j:ℝ)) ≠ 0 :=
    Finset.prod_ne_zero_iff.2 fun j _ => by positivity
  have hpow : (N:ℝ) ^ (1:ℝ) * (N:ℝ) ^ (3/2:ℝ) = (N:ℝ) ^ ((3 + ε) / 2) * (N:ℝ) ^ (1 - ε / 2) := by
    rw [← Real.rpow_add hN0, ← Real.rpow_add hN0]; congr 1; ring
  have h1 : (N:ℝ) ^ (1:ℝ) = (N:ℝ) ^ ((3 + ε) / 2) * (N:ℝ) ^ (1 - ε / 2) / (N:ℝ) ^ (3/2:ℝ) := by
    rw [← hpow]; field_simp
  have hxa : (N:ℝ) ^ ((3 + ε) / 2) ≠ 0 := (Real.rpow_pos_of_pos hN0 _).ne'
  have hxb : (N:ℝ) ^ (1 - ε / 2) ≠ 0 := (Real.rpow_pos_of_pos hN0 _).ne'
  have hf : (N.factorial : ℝ) ≠ 0 := by positivity
  unfold fct Real.GammaSeq
  rw [Finset.prod_div_distrib, Finset.prod_mul_distrib, Finset.prod_mul_distrib, h1]
  field_simp

lemma S2_value (ε : ℝ) (hε : 0 < ε) (hnat : ∀ j : ℕ, ε ≠ j) :
    ∑' k, hT ε 2 k = Real.Gamma 1 * Real.Gamma (3/2) /
      (Real.Gamma ((3 + ε) / 2) * Real.Gamma (1 - ε / 2)) := by
  have hGa : Real.Gamma ((3 + ε) / 2) ≠ 0 := (Real.Gamma_pos_of_pos (by linarith)).ne'
  have hGb : Real.Gamma (1 - ε / 2) ≠ 0 := Real.Gamma_ne_zero fun m h => by
    apply hnat (2 + 2 * m); push_cast; linarith
  have hlim : Tendsto (fun N : ℕ => ∏ j ∈ Finset.range (N+1), fct ε j) atTop
      (𝓝 (Real.Gamma 1 * Real.Gamma (3/2) /
        (Real.Gamma ((3 + ε) / 2) * Real.Gamma (1 - ε / 2)))) := by
    have h := ((Real.GammaSeq_tendsto_Gamma 1).mul (Real.GammaSeq_tendsto_Gamma (3/2))).div
      ((Real.GammaSeq_tendsto_Gamma ((3 + ε) / 2)).mul (Real.GammaSeq_tendsto_Gamma (1 - ε / 2)))
      (mul_ne_zero hGa hGb)
    refine h.congr' (eventually_atTop.2 ⟨1, fun N hN => ?_⟩)
    exact (prod_fct_eq ε hε hnat N hN).symm
  have hS := (S_tendsto ε hε).comp (tendsto_add_atTop_nat 1)
  have hprod := hlim.mul hS
  rw [mul_one] at hprod
  refine tendsto_nhds_unique (f := fun _ : ℕ => ∑' k, hT ε 2 k) tendsto_const_nhds
    (hprod.congr fun N => ?_)
  simp only [Function.comp]
  exact (S_iter ε hε (N+1)).symm

lemma final_alg (p q ε s co X Y U V : ℝ) (hq : q ^ 2 = p) (hq0 : q ≠ 0) (hs : s ≠ 0)
    (hco : co ≠ 0) (hY : Y ≠ 0) (hU : U ≠ 0) (he : 1 + ε ≠ 0)
    (hXY : X * Y = p / s) (hUV : U * V = p / co) :
    p * ε / (2 * s * co) * (1 * (q / 2) / ((1 + ε) / 2 * U * Y)) =
      ε / 2 * X * V / (q * (1 + ε)) := by
  have hX : X = p / s / Y := by rw [eq_div_iff hY]; exact hXY
  have hV : V = p / co / U := by rw [eq_div_iff hU, mul_comm]; exact hUV
  rw [hX, hV, ← hq]
  field_simp

theorem stmt8_proof (ε : ℝ) (hε : 0 < ε) (hnat : ∀ j : ℕ, ε ≠ j) :
    HasSum (fun k : ℕ => π * ε / Real.sin (π * ε) * lamTerm ε (Real.sin (π / 2 / 2) ^ 2) k)
      (Real.Gamma (1 + ε / 2) * Real.Gamma ((1 - ε) / 2) / (Real.sqrt π * (1 + ε))) := by
  have hz : Real.sin (π / 2 / 2) ^ 2 = 1 / 2 := by
    rw [show π / 2 / 2 = π / 4 by ring, Real.sin_pi_div_four, div_pow,
      Real.sq_sqrt (by norm_num)]; norm_num
  have hterm : ∀ k, lamTerm ε (1 / 2) k = hT ε 2 k := by
    intro k; unfold lamTerm hT; ring
  have hsum := (hT_summable ε 2 hε (by norm_num)).hasSum.mul_left (π * ε / Real.sin (π * ε))
  rw [S2_value ε hε hnat] at hsum
  convert hsum using 1
  · funext k; rw [hz, hterm]
  · -- Gamma algebra
    have hsin : Real.sin (π * ε) ≠ 0 := by
      intro h
      obtain ⟨n, hn⟩ := Real.sin_eq_zero_iff.1 h
      have hn' : (n:ℝ) = ε := by
        have := Real.pi_pos; field_simp at hn; linarith
      have hn0 : 0 ≤ n := by
        have : (0:ℝ) < n := by linarith
        exact_mod_cast this.le
      apply hnat n.toNat
      rw [← hn']; exact_mod_cast (Int.toNat_of_nonneg hn0).symm
    have h2 : Real.sin (π * ε) = 2 * Real.sin (π * (ε / 2)) * Real.cos (π * (ε / 2)) := by
      rw [← Real.sin_two_mul]; ring_nf
    rw [h2] at hsin ⊢
    have hs : Real.sin (π * (ε / 2)) ≠ 0 := by
      intro h; apply hsin; rw [h]; ring
    have hco : Real.cos (π * (ε / 2)) ≠ 0 := by
      intro h; apply hsin; rw [h]; ring
    have hG32 : Real.Gamma (3 / 2) = Real.sqrt π / 2 := by
      rw [show (3 / 2 : ℝ) = 1 / 2 + 1 by norm_num, Real.Gamma_add_one (by norm_num),
        Real.Gamma_one_half_eq]; ring
    have hGa : Real.Gamma ((3 + ε) / 2) = (1 + ε) / 2 * Real.Gamma ((1 + ε) / 2) := by
      rw [show (3 + ε) / 2 = (1 + ε) / 2 + 1 by ring, Real.Gamma_add_one (by positivity)]
    have hG1 : Real.Gamma (1 + ε / 2) = ε / 2 * Real.Gamma (ε / 2) := by
      rw [show 1 + ε / 2 = ε / 2 + 1 by ring, Real.Gamma_add_one (by positivity)]
    have r1 := Real.Gamma_mul_Gamma_one_sub (ε / 2)
    have r2 := Real.Gamma_mul_Gamma_one_sub ((1 + ε) / 2)
    rw [show 1 - (1 + ε) / 2 = (1 - ε) / 2 by ring,
      show π * ((1 + ε) / 2) = π * (ε / 2) + π / 2 by ring, Real.sin_add_pi_div_two] at r2
    have hY : Real.Gamma (1 - ε / 2) ≠ 0 := Real.Gamma_ne_zero fun m h => by
      apply hnat (2 + 2 * m); push_cast; linarith
    have hU : Real.Gamma ((1 + ε) / 2) ≠ 0 := (Real.Gamma_pos_of_pos (by linarith)).ne'
    rw [hG32, hGa, hG1, Real.Gamma_one]
    exact final_alg π (Real.sqrt π) ε _ _ _ _ _ _ (Real.sq_sqrt Real.pi_pos.le)
      (Real.sqrt_pos.2 Real.pi_pos).ne' hs hco hY hU (by positivity) r1 r2 |>.symm

end T13

theorem T13.proof_8 : T13.stmt_8 := fun ε hε hnat => T13.stmt8_proof ε hε hnat

