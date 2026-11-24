# W1 / Route A: function-level identity in c (working notes, updated as the work proceeds)

Folder: `research/10_frontier/W1_corner_family/A_function_level/`. Python: `.venv/bin/python`.
Labels: PROVED (exact algebra / complete argument), NUMERICAL, CONJECTURE.

(Interim file; the final summary is at the end once the work is complete.)

## Established so far

1. **Exact constant-free HPL form of the squircle pieces (PROVED, `a04_rfamily.py`).**
   Scale invariance v = a w turns t = (v−b)(v+a)/(4v) into (w−r)(w+1)/(2(1−r)w) with r = b/a = b².
   So U = U_r, V = V_r at r = b², where U_r, V_r are the same integrals with letters {0, ±1, ±r}.
   In c = (1−r)/(1+r) (r = b² ⟺ c = 1/√2) the GPL→HPL reduction (Goncharov differential, all letters
   → −∞ as c → 0) gives
   - U(c) = 16 weight-4 HPLs H_w(c), letters {0,±1}, rational coefficients, **no constants**;
   - V(c) = 24 weight-4 HPLs + ζ(2)·(4 weight-2 HPLs) + ζ(3)·(2 weight-1 HPLs).
   Numerical check at c = 1/√2 against direct quadrature of U, V: 1.9e-50, 2.4e-50 (`a04_rfamily.out`).
2. **These are identical, term by term, to the weight-4 part of G1's fitted corner formula B3**
   (U ↔ c⁰ part, V ↔ c⁻¹ part + ½ c¹ part) (`a05_compare_B3.out`: 0 difference terms). So the weight-4
   data of the corner family at the relevant combinations are now derived exactly, not fitted.
3. **Route A as stated is impossible (PROVED, exact symbol algebra).** As functions of c, U(c) and V(c)
   have nonzero Λ²B₂ (δ₂,₂) obstruction, so no classical R(c) (Li₄ + products) can equal them up to a
   constant (`a01_obstruction_check.out`). With L = (1−c)/(1+c) = τ², τ = tan(α/4):
   ob(U) = {L}₂∧{−L}₂,  ob(V) = ({c}₂+{−c}₂)∧{−L}₂ = ½{c²}₂∧{−L}₂.
   The obstruction survives pull-back to τ; the complementary-angle involution τ ↦ (1−τ)/(1+τ) (c ↔ s,
   fixing α = π/2) cannot cancel it, and no rational-prefactor combination can: the four vectors
   ob(F₀∘c), ob(F₀∘s), ob(F₋₁∘c), ob(F₋₁∘s) have exact rank 4 (`a02`, `a03`).
   The same obstruction is the one Kimi found for the s-family (s-family = this r-family with r = 1/s²).
