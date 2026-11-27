import Mathlib
import Zeta3.Wave9.R01_ZetaLadder
/-!
# S04 -- the zeta parts of `M_4(α)` at every angle

Target id: S04.
Paper label: `thm:M4alpha` (paper/ieee/main.tex, lines 1572-1583; proof idea 1585-1587; the three
consequences 1589, first two sentences; `F_2` printed at line 1615).
Supplement: paper/ieee/proofs/sec_g8_M4alpha.tex (setting 7-8; `g8:thm:M4alpha` 10-20; the
precise form `g8:eq:M4formal` with `F_1`, `F_2` 24-31; gap box 33; "Why a finite computation
proves the identity for every n" 35-50; the `ζ(3)` part by the ladder 54; M8, M9 56-58; the
`π⁴` part by hand 60-66; proof of the theorem 68-70; "What rests on what" 78-89).

Earlier Lean, imported and not restated (`Zeta3.Wave9.R01_ZetaLadder`, which imports
`Zeta3.Wave8.Q04_FlatAllK` and `Zeta3.Wave7.V05_FlatAngle`), all from frozen statement blocks:
`M_k(α) = V05.Mk k α` (the Taylor coefficient `M^{(k)}(0;α)/k!` of the corner integral),
`c = D04.c α = cos(α/2)`, the real zeta values `ζ(w) = X05.z w = (riemannZeta w).re`
(`zeta_nat_eq_tsum_of_gt_one` identifies `riemannZeta w` with `Σ_{n≥1} n^{-w}` for `w ≥ 2`, and
it is real there), the formal coefficients `P_{k,n} = Q04.Pf k n ∈ ℚ[𝔷_2, 𝔷_3, …]`
(Definition `g6:def:formal`; the variable `X w` is the formal zeta value `𝔷_w`), and the
supplement's parts `[z]M_k(α) = R01.zM s k α = Σ_{n≥1} ([𝔷_z]P_{k,n}) c^{2n}` (a `tsum`, always
used together with `R01.AllSummable α`, which says every such series converges absolutely).
`Q04.stmt_5` (proved) gives `M_k(α) = Σ_n P_{k,n}(ζ) c^{2n}`; `R01.stmt_7` (proved) is the
first rung `[ζ(3)]M_{k+1} = -[ζ(2)]M_k`, whose case `k = 3` is the `ζ(3)` part of this theorem
(`[ζ(3)]M_4 = -[ζ(2)]M_3 = -Ψ`); it is used, not restated.

Quoted verbatim (main.tex 1572-1583):

  \begin{theorem}[$M_4$ at every angle$^\dagger$]\label{thm:M4alpha}
  For $0<\alpha<\pi$, with $c=\cos\frac\alpha2$ and $\Psi$ the $\zeta(2)$ part of $M_3$
  (Theorem~\ref{thm:flat}),
  \[
  \begin{split}
  M_4(\alpha)={}&R_4(c)+\zeta(2)\Psi_2(c)-\zeta(3)\Psi(c)\\
  &+\frac{\pi^4}{2880}\Bigl[\Bigl(9c-\frac{30}c+\frac{21}{c^3}\Bigr)\operatorname{artanh}c+23-\frac{21}{c^2}\Bigr],
  \end{split}
  \]
  where $R_4$ is a rational combination of $109$ terms $c^eH_w(\pm c)$ with harmonic
  polylogarithms of weights two to five and $-3\le e\le1$, and $\Psi_2$ is a rational combination
  of $21$ such terms of weights one to three and the constant $\frac12$ (in our notes).
  \end{theorem}

Quoted verbatim (main.tex 1589): "The $\zeta(3)$ part of $M_4$ is $-\zeta(3)$ times the
$\zeta(2)$ part of $M_3$, ... The $\pi^4$ part is elementary: only $\operatorname{artanh}c$ and
rational functions of $c$."

Quoted verbatim (main.tex 1615): "Here $F_2=\bigl(-\frac c{16}-\frac1{8c}+\frac3{16c^3}\bigr)
\operatorname{artanh}c+\frac1{16}-\frac3{16c^2}$ is the $\zeta(2)^2$ part of $M_4$."

Quoted verbatim (supplement, lines 24-31, the precise form that is proved):

  M_4(\alpha)=R_4(c)+\zeta(2)\Psi_2(c)-\zeta(3)\Psi(c)+\zeta(4)F_1(c)+\zeta(2)^2F_2(c),
  F_1=\Bigl(\tfrac{7c}{16}-\tfrac5{8c}+\tfrac3{16c^3}\Bigr)\operatorname{artanh}c+\tfrac9{16}-\tfrac3{16c^2},\qquad
  F_2=\Bigl(-\tfrac c{16}-\tfrac1{8c}+\tfrac3{16c^3}\Bigr)\operatorname{artanh}c+\tfrac1{16}-\tfrac3{16c^2},
  "where, more precisely, $R_4$, $\Psi_2$, $\Psi$, $F_1$, $F_2$ are the sums
  $\sum_n[\gvs z]P_{4,n}c^{2n}$ for $\gvs z=1,\gvs2,-\gvs3,\gvs4,\gvs2^2$ ... Then
  $\frac{F_1}{90}+\frac{F_2}{36}=\frac1{2880}\bigl[(9c-\frac{30}c+\frac{21}{c^3})\operatorname{artanh}c+23-\frac{21}{c^2}\bigr]$".

Quoted verbatim (supplement, M9, line 58): "The same run prints
$[\gvs4]P_{4,n}=\frac{n+3}{(2n-1)(2n+1)(2n+3)}$ and
$[\gvs2^2]P_{4,n}=-\frac n{(2n-1)(2n+1)(2n+3)}$, each as an element of $\mathfrak D$ with only the
empty word, i.e.\ as a rational function of $n$".

Quoted verbatim (supplement, line 62, the three sums): for $0<c<1$,
  \sum_{n\ge1}\frac{c^{2n}}{2n-1}=c\operatorname{artanh}c,\quad
  \sum_{n\ge1}\frac{c^{2n}}{2n+1}=\frac{\operatorname{artanh}c}c-1,\quad
  \sum_{n\ge1}\frac{c^{2n}}{2n+3}=\frac{\operatorname{artanh}c}{c^3}-\frac1{c^2}-\frac13

Quoted verbatim (supplement, proof, line 69): "$P_{4,n}$ a polynomial in $\gvs2,\gvs3,\gvs4$ of
weight $\le4$: the monomials are $1,\gvs2,\gvs3,\gvs4,\gvs2^2$ (the only monomials of weight
$\le4$)."

## Formalization choices

* Everything about zeta parts is in the supplement's formal sense (gap box (iii): with real
  numbers the split is not unique). A part is `R01.zM s 4 α` with `s` the monomial: `1` is
  `s = 0`, `𝔷2` is `single 2 1`, `𝔷3` is `single 3 1`, `𝔷4` is `single 4 1`, `𝔷2²` is
  `single 2 2`. The real values `ζ(2), ζ(3), ζ(4)` are `X05.z 2, X05.z 3, X05.z 4`.
* `Ψ` ("the `ζ(2)` part of `M_3`", Theorem `thm:flat`; supplement line 7: `Ψ = Σ q'_n c^{2n}`,
  `q'_n = [𝔷2]P_{3,n}`) is `R01.zM (single 2 1) 3 α`. V05 has no closed form for `Ψ` (it is a
  series there too), so none is stated.
* `R_4 := [1]M_4(α) = R01.zM 0 4 α` and `Ψ_2 := [ζ(2)]M_4(α) = R01.zM (single 2 1) 4 α` are the
  SERIES parts, as in the supplement's "more precisely" sentence. They are not the 109-term and
  21+1-term closed forms in harmonic polylogarithms (see NOT COVERED).
* The `ζ(3)` part `[ζ(3)]M_4 = -Ψ` is `R01.stmt_7` at `k = 3` (proved there); it is not
  restated as its own statement, but the final statement contains it in the form `-ζ(3)Ψ`.
* `F_1`, `F_2` are the supplement's closed forms (`F_2` is also printed in main.tex 1615), as
  functions of a real `c`; they are only evaluated at `0 < c < 1`, so `1/c` is not a junk value.
  The `π⁴` bracket is written out literally in stmt_6 and stmt_7.
* stmt_5 is stated for every `0 < c < 1` (the formal coefficients do not depend on the angle);
  stmt_6 and stmt_7 at `c = cos(α/2)`, `0 < α < π`, as in the paper.
* Sums over `n ≥ 1` are over `n : ℕ` with index `n + 1`; `HasSum` is used except inside
  `R01.zM`, which is always guarded by `R01.AllSummable`.

## Statements

* stmt_1: `P_{4,n}` has only the monomials `1, 𝔷2, 𝔷3, 𝔷4, 𝔷2²`, for every `n ≥ 1`.
* stmt_2: for `0 < α < π`, `M_4(α) = [1]M_4 + ζ(2)[ζ(2)]M_4 + ζ(3)[ζ(3)]M_4 + ζ(4)[ζ(4)]M_4
  + ζ(2)²[ζ(2)²]M_4`, all parts absolutely convergent.
* stmt_3 (M9): `[𝔷4]P_{4,n} = (n+3)/((2n-1)(2n+1)(2n+3))`,
  `[𝔷2²]P_{4,n} = -n/((2n-1)(2n+1)(2n+3))` for every `n ≥ 1`.
* stmt_4: the three `artanh` series for `0 < c < 1`.
* stmt_5: `Σ_n [𝔷4]P_{4,n} c^{2n} = F_1(c)` and `Σ_n [𝔷2²]P_{4,n} c^{2n} = F_2(c)`, `0 < c < 1`.
* stmt_6: the `π⁴` part: `ζ(4)[ζ(4)]M_4(α) + ζ(2)²[ζ(2)²]M_4(α) = (π⁴/2880)[…]`.
* stmt_7: `thm:M4alpha` with `R_4`, `Ψ_2` the series parts:
  `M_4(α) = R_4 + ζ(2)Ψ_2 - ζ(3)Ψ + (π⁴/2880)[(9c - 30/c + 21/c³) artanh c + 23 - 21/c²]`.

## Numerical check before freezing

A SymPy script mirroring `Q04.Pf` exactly (the same truncations of `E_n` and `exp Δ_⋆`, exact
rationals) gave for `n = 1, …, 10`: the support of `P_{4,n}` lies in `{1, 𝔷2, 𝔷3, 𝔷4, 𝔷2²}`,
and `[𝔷4]P_{4,n}`, `[𝔷2²]P_{4,n}` equal the two rational functions of stmt_3 exactly (e.g.
`n = 1`: `4/15`, `-1/15`, as in `P_{4,1}` of Q04). The Taylor series of `F_1`, `F_2` and of the
bracket agree exactly with `Σ f(n) c^{2n}` for the two rational functions `f` (and with
`32 f_4 + 80 f_{22}`) through `c^{25}` (odd coefficients `0`), and `F_1/90 + F_2/36 - bracket/2880`
simplifies to `0`. With mpmath (40 digits), `Σ_n (ζ(4) f_4(n) + ζ(2)² f_{22}(n)) c^{2n}` agrees
with `(π⁴/2880)[…]` to `2·10⁻⁴¹` at `α = π/3` and `4·10⁻⁴¹` at `α = 2.5`.

## NOT COVERED

* The closed forms of `R_4` (109 terms `c^e H_w(±c)`, weights 2 to 5, `-3 ≤ e ≤ 1`) and of `Ψ_2`
  (21 terms of weights 1 to 3 and the constant `1/2`): their coefficients are only in the
  authors' data files (`g8_terms_k4_part1.json`, `g8_terms_k4_partz2.json`), not printed in the
  paper or the supplement, so they cannot be stated from the text. `R_4` and `Ψ_2` appear only
  as the series parts.
* The harmonic polylogarithms `H_w(±c)`, the algebra `𝔇` of Z-sums and Lemma `g8:prop:crit`.
* The depth-two remarks after the theorem (the `(3,2)` symbol pairings at the generic angle, at
  the hexagon corner, at `c_0 = 1/3`, at the golden corners, `N_5^{hex}`), and the numerical
  checks of the proof idea.
* Statements about the zeta parts of the real number `M_4(α)` (not well defined without the
  independence of `1, π², ζ(3)`, gap box (iii)).
-/
-- BEGIN STATEMENTS
namespace S04

open Real

/-- The supplement's `F_1(c) = (7c/16 - 5/(8c) + 3/(16c³)) artanh c + 9/16 - 3/(16c²)`. -/
noncomputable def F1 (c : ℝ) : ℝ :=
  (7 * c / 16 - 5 / (8 * c) + 3 / (16 * c ^ 3)) * Real.artanh c + 9 / 16 - 3 / (16 * c ^ 2)

/-- `F_2(c) = (-c/16 - 1/(8c) + 3/(16c³)) artanh c + 1/16 - 3/(16c²)` (main.tex 1615, supplement). -/
noncomputable def F2 (c : ℝ) : ℝ :=
  (-c / 16 - 1 / (8 * c) + 3 / (16 * c ^ 3)) * Real.artanh c + 1 / 16 - 3 / (16 * c ^ 2)

/-- Supplement, proof (line 69): `P_{4,n}` is a polynomial in `𝔷2, 𝔷3, 𝔷4` whose monomials are
among `1, 𝔷2, 𝔷3, 𝔷4, 𝔷2²`, for every `n ≥ 1`. -/
def stmt_1 : Prop :=
  ∀ n : ℕ, 1 ≤ n →
    Q04.Pf 4 n =
      MvPolynomial.C ((Q04.Pf 4 n).coeff 0)
      + MvPolynomial.C ((Q04.Pf 4 n).coeff (Finsupp.single 2 1)) * MvPolynomial.X 2
      + MvPolynomial.C ((Q04.Pf 4 n).coeff (Finsupp.single 3 1)) * MvPolynomial.X 3
      + MvPolynomial.C ((Q04.Pf 4 n).coeff (Finsupp.single 4 1)) * MvPolynomial.X 4
      + MvPolynomial.C ((Q04.Pf 4 n).coeff (Finsupp.single 2 2)) * MvPolynomial.X 2 ^ 2

/-- The decomposition behind `g8:eq:M4formal`: for `0 < α < π` all parts converge absolutely and
`M_4(α) = [1]M_4 + ζ(2)[ζ(2)]M_4 + ζ(3)[ζ(3)]M_4 + ζ(4)[ζ(4)]M_4 + ζ(2)²[ζ(2)²]M_4`. -/
def stmt_2 : Prop :=
  ∀ α : ℝ, 0 < α → α < π →
    R01.AllSummable α ∧
    V05.Mk 4 α =
      R01.zM 0 4 α + X05.z 2 * R01.zM (Finsupp.single 2 1) 4 α
        + X05.z 3 * R01.zM (Finsupp.single 3 1) 4 α
        + X05.z 4 * R01.zM (Finsupp.single 4 1) 4 α
        + X05.z 2 ^ 2 * R01.zM (Finsupp.single 2 2) 4 α

/-- M9 (supplement line 58), for every `n ≥ 1`:
`[𝔷4]P_{4,n} = (n+3)/((2n-1)(2n+1)(2n+3))` and `[𝔷2²]P_{4,n} = -n/((2n-1)(2n+1)(2n+3))`. -/
def stmt_3 : Prop :=
  ∀ n : ℕ, 1 ≤ n →
    (Q04.Pf 4 n).coeff (Finsupp.single 4 1) =
      ((n : ℚ) + 3) / ((2 * (n : ℚ) - 1) * (2 * (n : ℚ) + 1) * (2 * (n : ℚ) + 3)) ∧
    (Q04.Pf 4 n).coeff (Finsupp.single 2 2) =
      -(n : ℚ) / ((2 * (n : ℚ) - 1) * (2 * (n : ℚ) + 1) * (2 * (n : ℚ) + 3))

/-- The three sums of the supplement (line 62), for `0 < c < 1`:
`Σ_{n≥1} c^{2n}/(2n-1) = c artanh c`, `Σ_{n≥1} c^{2n}/(2n+1) = artanh c / c - 1`,
`Σ_{n≥1} c^{2n}/(2n+3) = artanh c / c³ - 1/c² - 1/3`. -/
def stmt_4 : Prop :=
  ∀ c : ℝ, 0 < c → c < 1 →
    HasSum (fun n : ℕ => c ^ (2 * (n + 1)) / (2 * ((n : ℝ) + 1) - 1)) (c * Real.artanh c) ∧
    HasSum (fun n : ℕ => c ^ (2 * (n + 1)) / (2 * ((n : ℝ) + 1) + 1)) (Real.artanh c / c - 1) ∧
    HasSum (fun n : ℕ => c ^ (2 * (n + 1)) / (2 * ((n : ℝ) + 1) + 3))
      (Real.artanh c / c ^ 3 - 1 / c ^ 2 - 1 / 3)

/-- The supplement's `F_1`, `F_2` (line 66): for `0 < c < 1`,
`Σ_{n≥1} [𝔷4]P_{4,n} c^{2n} = F_1(c)` and `Σ_{n≥1} [𝔷2²]P_{4,n} c^{2n} = F_2(c)`. -/
def stmt_5 : Prop :=
  ∀ c : ℝ, 0 < c → c < 1 →
    HasSum (fun n : ℕ => (((Q04.Pf 4 (n + 1)).coeff (Finsupp.single 4 1) : ℚ) : ℝ)
      * c ^ (2 * (n + 1))) (F1 c) ∧
    HasSum (fun n : ℕ => (((Q04.Pf 4 (n + 1)).coeff (Finsupp.single 2 2) : ℚ) : ℝ)
      * c ^ (2 * (n + 1))) (F2 c)

/-- The `π⁴` part of `thm:M4alpha` ("The `π⁴` part is elementary"): for `0 < α < π`, `c = cos(α/2)`,
`ζ(4)[ζ(4)]M_4(α) + ζ(2)²[ζ(2)²]M_4(α)
  = (π⁴/2880)[(9c - 30/c + 21/c³) artanh c + 23 - 21/c²]`. -/
def stmt_6 : Prop :=
  ∀ α : ℝ, 0 < α → α < π →
    R01.AllSummable α ∧
    X05.z 4 * R01.zM (Finsupp.single 4 1) 4 α + X05.z 2 ^ 2 * R01.zM (Finsupp.single 2 2) 4 α =
      π ^ 4 / 2880 *
        ((9 * D04.c α - 30 / D04.c α + 21 / D04.c α ^ 3) * Real.artanh (D04.c α)
          + 23 - 21 / D04.c α ^ 2)

/-- `thm:M4alpha`, with `R_4 = [1]M_4(α)` and `Ψ_2 = [ζ(2)]M_4(α)` the series parts (not their
109-term and 21+1-term closed forms) and `Ψ = [ζ(2)]M_3(α)`: for `0 < α < π`, `c = cos(α/2)`, all
parts converge absolutely and
`M_4(α) = R_4 + ζ(2)Ψ_2 - ζ(3)Ψ + (π⁴/2880)[(9c - 30/c + 21/c³) artanh c + 23 - 21/c²]`. -/
def stmt_7 : Prop :=
  ∀ α : ℝ, 0 < α → α < π →
    R01.AllSummable α ∧
    V05.Mk 4 α =
      R01.zM 0 4 α + X05.z 2 * R01.zM (Finsupp.single 2 1) 4 α
        - X05.z 3 * R01.zM (Finsupp.single 2 1) 3 α
        + π ^ 4 / 2880 *
          ((9 * D04.c α - 30 / D04.c α + 21 / D04.c α ^ 3) * Real.artanh (D04.c α)
            + 23 - 21 / D04.c α ^ 2)

end S04
-- END STATEMENTS

namespace S04
open Polynomial

abbrev RR := MvPolynomial (Fin 3) ℚ

lemma cmul (p q : RR[X]) (k : ℕ) :
    (p * q).coeff k = ∑ i ∈ Finset.range (k + 1), p.coeff i * q.coeff (k - i) := by
  rw [Polynomial.coeff_mul, Finset.Nat.sum_antidiagonal_eq_sum_range_succ (fun i j => p.coeff i * q.coeff j)]

lemma csmul (q : ℚ) (p : RR[X]) (k : ℕ) : (q • p).coeff k = q • p.coeff k :=
  AddMonoidAlgebra.coeff_smul_apply ..

lemma pow2c (p : RR[X]) (h0 : p.coeff 0 = 0) :
    (p ^ 2).coeff 0 = 0 ∧ (p ^ 2).coeff 1 = 0 ∧ (p ^ 2).coeff 2 = p.coeff 1 ^ 2 ∧
    (p ^ 2).coeff 3 = 2 * p.coeff 1 * p.coeff 2 ∧
    (p ^ 2).coeff 4 = 2 * p.coeff 1 * p.coeff 3 + p.coeff 2 ^ 2 := by
  simp only [pow_two, cmul, Finset.sum_range_succ, Finset.sum_range_zero, h0]
  refine ⟨by ring, by ring, by ring, by ring, by ring⟩

lemma pow3c (p : RR[X]) (h0 : p.coeff 0 = 0) :
    (p ^ 3).coeff 0 = 0 ∧ (p ^ 3).coeff 1 = 0 ∧ (p ^ 3).coeff 2 = 0 ∧
    (p ^ 3).coeff 3 = p.coeff 1 ^ 3 ∧ (p ^ 3).coeff 4 = 3 * p.coeff 1 ^ 2 * p.coeff 2 := by
  obtain ⟨a0, a1, a2, a3, a4⟩ := pow2c p h0
  rw [pow_succ]
  simp only [cmul, Finset.sum_range_succ, Finset.sum_range_zero, h0, a0, a1, a2, a3, a4]
  refine ⟨by ring, by ring, by ring, by ring, by ring⟩

lemma pow4c (p : RR[X]) (h0 : p.coeff 0 = 0) :
    (p ^ 4).coeff 0 = 0 ∧ (p ^ 4).coeff 1 = 0 ∧ (p ^ 4).coeff 2 = 0 ∧
    (p ^ 4).coeff 3 = 0 ∧ (p ^ 4).coeff 4 = p.coeff 1 ^ 4 := by
  obtain ⟨a0, a1, a2, a3, a4⟩ := pow3c p h0
  rw [pow_succ]
  simp only [cmul, Finset.sum_range_succ, Finset.sum_range_zero, h0, a0, a1, a2, a3, a4]
  refine ⟨by ring, by ring, by ring, by ring, by ring⟩

lemma Icc14 : Finset.Icc 1 4 = {1, 2, 3, 4} := by decide
lemma Icc12 : Finset.Icc 1 2 = {1, 2} := by decide

lemma expT4c (L : RR[X]) (h0 : L.coeff 0 = 0) :
    (Q04.expT 4 L).coeff 4 = L.coeff 4 + ((1:ℚ)/2) • (2 * L.coeff 1 * L.coeff 3 + L.coeff 2 ^ 2)
      + ((1:ℚ)/6) • (3 * L.coeff 1 ^ 2 * L.coeff 2) + ((1:ℚ)/24) • L.coeff 1 ^ 4 := by
  obtain ⟨-, -, -, -, b4⟩ := pow2c L h0
  obtain ⟨-, -, -, -, c4⟩ := pow3c L h0
  obtain ⟨-, -, -, -, d4⟩ := pow4c L h0
  have e1 : (1 : RR[X]).coeff 4 = 0 := Polynomial.coeff_one.trans (if_neg (by norm_num))
  simp only [Q04.expT, Polynomial.finsetSum_coeff, csmul, Finset.sum_range_succ,
    Finset.sum_range_zero, pow_zero, pow_one, Polynomial.coeff_add, e1, b4, c4, d4]
  norm_num [Nat.factorial]

lemma logT4_eq (w : RR[X]) : Q04.logT 4 w =
    w + (-1/2:ℚ) • w ^ 2 + (1/3:ℚ) • w ^ 3 + (-1/4:ℚ) • w ^ 4 := by
  simp only [Q04.logT, Icc14]
  rw [Finset.sum_insert (by decide), Finset.sum_insert (by decide), Finset.sum_insert (by decide),
    Finset.sum_singleton]
  norm_num
  abel

lemma logT4c (w : RR[X]) (h0 : w.coeff 0 = 0) :
    (Q04.logT 4 w).coeff 0 = 0 ∧
    (Q04.logT 4 w).coeff 1 = w.coeff 1 ∧
    (Q04.logT 4 w).coeff 2 = w.coeff 2 - ((1:ℚ)/2) • w.coeff 1 ^ 2 ∧
    (Q04.logT 4 w).coeff 3 = w.coeff 3 - w.coeff 1 * w.coeff 2 + ((1:ℚ)/3) • w.coeff 1 ^ 3 ∧
    (Q04.logT 4 w).coeff 4 = w.coeff 4 - ((1:ℚ)/2) • (2 * w.coeff 1 * w.coeff 3 + w.coeff 2 ^ 2)
      + w.coeff 1 ^ 2 * w.coeff 2 - ((1:ℚ)/4) • w.coeff 1 ^ 4 := by
  obtain ⟨b0, b1, b2, b3, b4⟩ := pow2c w h0
  obtain ⟨c0, c1, c2, c3, c4⟩ := pow3c w h0
  obtain ⟨d0, d1, d2, d3, d4⟩ := pow4c w h0
  simp only [logT4_eq, Polynomial.coeff_add, csmul,
    b0, b1, b2, b3, b4, c0, c1, c2, c3, c4, d0, d1, d2, d3, d4, h0]
  refine ⟨by simp, by simp, by algebra, by algebra, by algebra⟩

/-- `ℓ = b - a` -/
noncomputable def ell : RR := MvPolynomial.X 2 - MvPolynomial.X 1

lemma wK4c : (Q04.wK 4).coeff 0 = 0 ∧ (Q04.wK 4).coeff 1 = MvPolynomial.X 0 * ell ∧
    (Q04.wK 4).coeff 2 = ((1:ℚ)/2) • (MvPolynomial.X 0 * ell ^ 2) ∧
    (Q04.wK 4).coeff 3 = ((1:ℚ)/6) • (MvPolynomial.X 0 * ell ^ 3) ∧
    (Q04.wK 4).coeff 4 = ((1:ℚ)/24) • (MvPolynomial.X 0 * ell ^ 4) := by
  simp only [Q04.wK, Icc14]
  rw [Finset.sum_insert (by decide), Finset.sum_insert (by decide), Finset.sum_insert (by decide),
    Finset.sum_singleton]
  simp only [pow_one, mul_pow, ← Polynomial.C_pow, Polynomial.coeff_C_mul, Polynomial.coeff_add, csmul,
    Polynomial.coeff_C_mul_X_pow, Polynomial.coeff_X, ell]
  refine ⟨by simp, by simp, ?_, ?_, ?_⟩ <;> norm_num [Nat.factorial, Polynomial.coeff_X]

lemma monoX (j i m : ℕ) (c : ℚ) : MvPolynomial.monomial (R01.mono j i m) c =
    algebraMap ℚ RR c * MvPolynomial.X 0 ^ j * MvPolynomial.X 1 ^ i * MvPolynomial.X 2 ^ m := by
  simp only [R01.mono, MvPolynomial.X_pow_eq_monomial, MvPolynomial.algebraMap_eq,
    MvPolynomial.C_mul_monomial, MvPolynomial.monomial_mul, mul_one]

/-- the coefficients `e^{(4)}_{j,i,4-i}(n)` of `y^j a^i b^{4-i}` in `E_{n,4}` -/
def cf4 (n : ℚ) : ℕ → ℕ → ℚ
  | 0, 0 => n^4/24
  | 0, 1 => -n^4/6 + n^3/6
  | 0, 2 => n^4/4 - n^3/2 + n^2/4
  | 0, 3 => -n^4/6 + n^3/2 - n^2/2 + n/6
  | 0, 4 => n^4/24 - n^3/6 + n^2/4 - n/6 + 1/24
  | 1, 0 => -n^4/3 - n^3/3 - n^2/12 + n/12 + 1/24
  | 1, 1 => 4*n^4/3 + n^3/3 - n^2/6 - n/6
  | 1, 2 => -2*n^4 + n^3
  | 1, 3 => 4*n^4/3 - 5*n^3/3 + 5*n^2/6 - n/6
  | 1, 4 => -n^4/3 + 2*n^3/3 - 7*n^2/12 + n/4 - 1/24
  | 2, 0 => n^4 + 3*n^3/2 + n^2/6 - 7*n/12
  | 2, 1 => -4*n^4 - 4*n^3 + n^2/3 + 4*n/3
  | 2, 2 => 6*n^4 + 3*n^3 - n^2 - n
  | 2, 3 => -4*n^4 + n^2/3 + n/3
  | 2, 4 => n^4 - n^3/2 + n^2/6 - n/12
  | 3, 0 => -4*n^4/3 - 2*n^3 + n^2/3 + n/2
  | 3, 1 => 16*n^4/3 + 20*n^3/3 - 4*n^2/3 - 5*n/3
  | 3, 2 => -8*n^4 - 8*n^3 + 2*n^2 + 2*n
  | 3, 3 => 16*n^4/3 + 4*n^3 - 4*n^2/3 - n
  | 3, 4 => -4*n^4/3 - 2*n^3/3 + n^2/3 + n/6
  | 4, 0 => 2*n^4/3 + 2*n^3/3 - n^2/6 - n/6
  | 4, 1 => -8*n^4/3 - 8*n^3/3 + 2*n^2/3 + 2*n/3
  | 4, 2 => 4*n^4 + 4*n^3 - n^2 - n
  | 4, 3 => -8*n^4/3 - 8*n^3/3 + 2*n^2/3 + 2*n/3
  | 4, 4 => 2*n^4/3 + 2*n^3/3 - n^2/6 - n/6
  | _, _ => 0

lemma LK4c (n : ℕ) :
    (Q04.LK n 4).coeff 0 = 0 ∧
    (Q04.LK n 4).coeff 1 = ((-((n : ℚ) - 1)) • (MvPolynomial.X 1 : RR) + (n : ℚ) • MvPolynomial.X 2)
      + (1 - 2 * (n : ℚ)) • (Q04.wK 4).coeff 1 ∧
    (Q04.LK n 4).coeff 2 = (1 - 2 * (n : ℚ)) • (Q04.logT 4 (Q04.wK 4)).coeff 2 ∧
    (Q04.LK n 4).coeff 3 = (1 - 2 * (n : ℚ)) • (Q04.logT 4 (Q04.wK 4)).coeff 3 ∧
    (Q04.LK n 4).coeff 4 = (1 - 2 * (n : ℚ)) • (Q04.logT 4 (Q04.wK 4)).coeff 4 := by
  obtain ⟨w0, -⟩ := wK4c
  obtain ⟨l0, l1, -⟩ := logT4c (Q04.wK 4) w0
  simp only [Q04.LK_eq, Polynomial.coeff_add, csmul, Polynomial.coeff_C_mul_X, l0, l1]
  norm_num

lemma Ecoef4 (n : ℕ) : Q04.Ecoef 4 n = ∑ j ∈ Finset.range 5, ∑ i ∈ Finset.range 5,
    MvPolynomial.monomial (R01.mono j i (4 - i)) (cf4 n j i) := by
  obtain ⟨w0, w1, w2, w3, w4⟩ := wK4c
  obtain ⟨-, -, l2, l3, l4⟩ := logT4c (Q04.wK 4) w0
  obtain ⟨L0, L1, L2, L3, L4⟩ := LK4c n
  rw [Q04.Ecoef, Q04.Eser_eq, expT4c _ L0, L1, L2, L3, L4, l2, l3, l4, w1, w2, w3, w4]
  simp only [Finset.sum_range_succ, Finset.sum_range_zero, monoX, cf4, ell, zero_add]
  algebra

lemma expT2c (L : RR[X]) (h0 : L.coeff 0 = 0) :
    (Q04.expT 2 L).coeff 2 = L.coeff 2 + ((1:ℚ)/2) • L.coeff 1 ^ 2 := by
  obtain ⟨-, -, b2, -, -⟩ := pow2c L h0
  have e1 : (1 : RR[X]).coeff 2 = 0 := Polynomial.coeff_one.trans (if_neg (by norm_num))
  simp only [Q04.expT, Polynomial.finsetSum_coeff, csmul, Finset.sum_range_succ,
    Finset.sum_range_zero, pow_zero, pow_one, Polynomial.coeff_add, e1, b2]
  norm_num [Nat.factorial]

lemma logT2_eq (w : RR[X]) : Q04.logT 2 w = w + (-1/2:ℚ) • w ^ 2 := by
  simp only [Q04.logT, Icc12]
  rw [Finset.sum_insert (by decide), Finset.sum_singleton]
  norm_num

lemma wK2c : (Q04.wK 2).coeff 0 = 0 ∧ (Q04.wK 2).coeff 1 = MvPolynomial.X 0 * ell ∧
    (Q04.wK 2).coeff 2 = ((1:ℚ)/2) • (MvPolynomial.X 0 * ell ^ 2) := by
  simp only [Q04.wK, Icc12]
  rw [Finset.sum_insert (by decide), Finset.sum_singleton]
  simp only [pow_one, mul_pow, ← Polynomial.C_pow, Polynomial.coeff_C_mul, Polynomial.coeff_add, csmul,
    Polynomial.coeff_C_mul_X_pow, Polynomial.coeff_X, ell]
  refine ⟨by simp, by simp, ?_⟩
  norm_num [Nat.factorial, Polynomial.coeff_X]

/-- the coefficients `e^{(2)}_{j,i,2-i}(n)` of `y^j a^i b^{2-i}` in `E_{n,2}` -/
def cf2 (n : ℚ) : ℕ → ℕ → ℚ
  | 0, 0 => n^2/2
  | 0, 1 => -n^2 + n
  | 0, 2 => n^2/2 - n + 1/2
  | 1, 0 => 1/2 - 2*n^2
  | 1, 1 => 4*n^2 - 2*n
  | 1, 2 => -2*n^2 + 2*n - 1/2
  | 2, 0 => 2*n^2 - n
  | 2, 1 => -4*n^2 + 2*n
  | 2, 2 => 2*n^2 - n
  | _, _ => 0

lemma Ecoef2 (n : ℕ) : Q04.Ecoef 2 n = ∑ j ∈ Finset.range 3, ∑ i ∈ Finset.range 3,
    MvPolynomial.monomial (R01.mono j i (2 - i)) (cf2 n j i) := by
  obtain ⟨w0, w1, w2⟩ := wK2c
  have hlg : (Q04.logT 2 (Q04.wK 2)).coeff 0 = 0 ∧ (Q04.logT 2 (Q04.wK 2)).coeff 1 = (Q04.wK 2).coeff 1 ∧
      (Q04.logT 2 (Q04.wK 2)).coeff 2 = (Q04.wK 2).coeff 2 - ((1:ℚ)/2) • (Q04.wK 2).coeff 1 ^ 2 := by
    obtain ⟨b0, b1, b2, -, -⟩ := pow2c (Q04.wK 2) w0
    simp only [logT2_eq, Polynomial.coeff_add, csmul, b0, b1, b2, w0]
    refine ⟨by simp, by simp, by algebra⟩
  obtain ⟨l0, l1, l2⟩ := hlg
  have L0 : (Q04.LK n 2).coeff 0 = 0 := by
    simp only [Q04.LK_eq, Polynomial.coeff_add, csmul, Polynomial.coeff_C_mul_X, l0]; simp
  have L1 : (Q04.LK n 2).coeff 1 = ((-((n : ℚ) - 1)) • (MvPolynomial.X 1 : RR) + (n : ℚ) • MvPolynomial.X 2)
      + (1 - 2 * (n : ℚ)) • (Q04.wK 2).coeff 1 := by
    simp only [Q04.LK_eq, Polynomial.coeff_add, csmul, Polynomial.coeff_C_mul_X, l1]; simp
  have L2 : (Q04.LK n 2).coeff 2 = (1 - 2 * (n : ℚ)) • (Q04.logT 2 (Q04.wK 2)).coeff 2 := by
    simp only [Q04.LK_eq, Polynomial.coeff_add, csmul, Polynomial.coeff_C_mul_X]; simp
  rw [Q04.Ecoef, Q04.Eser_eq, expT2c _ L0, L1, L2, l2, w1, w2]
  simp only [Finset.sum_range_succ, Finset.sum_range_zero, monoX, cf2, ell, zero_add]
  algebra

lemma mono_inj {j i m j' i' m' : ℕ} (h : R01.mono j i m = R01.mono j' i' m') :
    j = j' ∧ i = i' ∧ m = m' := by
  obtain ⟨a0, a1, a2⟩ := R01.mono_apply j i m
  obtain ⟨b0, b1, b2⟩ := R01.mono_apply j' i' m'
  rw [h] at a0 a1 a2
  exact ⟨a0.symm.trans b0, a1.symm.trans b1, a2.symm.trans b2⟩

lemma coeff_table (K : ℕ) (f : ℕ → ℕ → ℚ) (j i : ℕ) (hj : j ≤ K) (hi : i ≤ K) :
    (∑ j' ∈ Finset.range (K + 1), ∑ i' ∈ Finset.range (K + 1),
      MvPolynomial.monomial (R01.mono j' i' (K - i')) (f j' i')).coeff (R01.mono j i (K - i)) =
      f j i := by
  classical
  simp only [MvPolynomial.coeff_sum, MvPolynomial.coeff_monomial]
  rw [Finset.sum_eq_single j, Finset.sum_eq_single i]
  · simp
  · intro b _ hb; rw [if_neg]; intro h; exact hb (mono_inj h).2.1
  · intro h; exact absurd (Finset.mem_range.mpr (by omega)) h
  · intro b _ hb
    refine Finset.sum_eq_zero fun c _ => ?_
    rw [if_neg]; intro h; exact hb (mono_inj h).1
  · intro h; exact absurd (Finset.mem_range.mpr (by omega)) h

lemma eF4 (n j i : ℕ) (hj : j ≤ 4) (hi : i ≤ 4) : Q04.eF 4 n j i (4 - i) = cf4 n j i := by
  rw [Q04.eF, ← R01.mono, Ecoef4]; exact coeff_table 4 _ j i hj hi

lemma eF2 (n j i : ℕ) (hj : j ≤ 2) (hi : i ≤ 2) : Q04.eF 2 n j i (2 - i) = cf2 n j i := by
  rw [Q04.eF, ← R01.mono, Ecoef2]; exact coeff_table 2 _ j i hj hi

lemma D00 (p q : ℕ) : (Q04.DsharpF 0 0 p q).coeff 0 = 1 := by
  simp [Q04.DsharpF, Q04.expT]

lemma Dz2e (e i m p q : ℕ) (he : e ≠ 0) :
    (Q04.DsharpF i m p q).coeff (Finsupp.single 2 e) =
      (-1 : ℚ) ^ e * ((e.factorial : ℚ)⁻¹) * ((i.descFactorial e * m.descFactorial e : ℕ) : ℚ) *
        (Q04.DsharpF (i - e) (m - e) p q).coeff 0 := by
  classical
  have hQ : R01.Qz (Finsupp.single 2 e) =
      MvPolynomial.monomial (Finsupp.single 0 e + Finsupp.single 1 e) ((-1 : ℚ) ^ e * ((e.factorial : ℚ)⁻¹)) := by
    rw [R01.Qz_single2 e he, mul_pow, MvPolynomial.X_pow_eq_monomial, MvPolynomial.X_pow_eq_monomial,
      MvPolynomial.monomial_mul, MvPolynomial.C_mul_monomial, mul_one, mul_one]
  have hc : ((-1 : ℚ) ^ e * ((e.factorial : ℚ)⁻¹)) ≠ 0 := by
    have : (e.factorial : ℚ) ≠ 0 := by positivity
    simp [this]
  rw [R01.core i m p q _ (R01.zeta2_mono e), hQ, MvPolynomial.support_monomial, if_neg hc,
    Finset.sum_singleton, MvPolynomial.coeff_monomial, if_pos rfl]
  simp only [Finsupp.coe_add, Pi.add_apply, Finsupp.single_eq_same, Finsupp.single_apply]
  norm_num
  ring

lemma D22 (p q : ℕ) : (Q04.DsharpF 0 4 p q).coeff (Finsupp.single 2 2) = 0 ∧
    (Q04.DsharpF 1 3 p q).coeff (Finsupp.single 2 2) = 0 ∧
    (Q04.DsharpF 2 2 p q).coeff (Finsupp.single 2 2) = 2 ∧
    (Q04.DsharpF 3 1 p q).coeff (Finsupp.single 2 2) = 0 ∧
    (Q04.DsharpF 4 0 p q).coeff (Finsupp.single 2 2) = 0 := by
  simp only [Dz2e 2 _ _ p q two_ne_zero]
  refine ⟨by simp, by simp, ?_, by simp, by simp⟩
  rw [D00]; norm_num [Nat.descFactorial, Nat.factorial]

lemma D21 (p q : ℕ) : (Q04.DsharpF 0 2 p q).coeff (Finsupp.single 2 1) = 0 ∧
    (Q04.DsharpF 1 1 p q).coeff (Finsupp.single 2 1) = -1 ∧
    (Q04.DsharpF 2 0 p q).coeff (Finsupp.single 2 1) = 0 := by
  simp only [Dz2e 1 _ _ p q one_ne_zero]
  refine ⟨by simp, ?_, by simp⟩
  rw [D00]; norm_num

lemma z22 (n : ℕ) (hn : 1 ≤ n) : (Q04.Pf 4 n).coeff (Finsupp.single 2 2) =
    -(n : ℚ) / ((2 * (n : ℚ) - 1) * (2 * (n : ℚ) + 1) * (2 * (n : ℚ) + 3)) := by
  have d0 := fun p q => (D22 p q).1
  have d1 := fun p q => (D22 p q).2.1
  have d2 := fun p q => (D22 p q).2.2.1
  have d3 := fun p q => (D22 p q).2.2.2.1
  have d4 := fun p q => (D22 p q).2.2.2.2
  rw [R01.proof_2 4 n hn, MvPolynomial.coeff_smul, smul_eq_mul, R01.PhiN_Ecoef]
  simp only [MvPolynomial.coeff_sum, MvPolynomial.coeff_add, MvPolynomial.coeff_smul, smul_eq_mul, Finset.sum_range_succ,
    Finset.sum_range_zero, zero_add, d0, d1, d2, d3, d4, mul_zero, add_zero, zero_add]
  have e : ∀ j, j ≤ 4 → Q04.eF 4 n j 2 2 = cf4 n j 2 := fun j hj => eF4 n j 2 hj (by norm_num)
  rw [e 0 (by norm_num), e 1 (by norm_num), e 2 (by norm_num), e 3 (by norm_num), e 4 (by norm_num)]
  simp only [cf4, V05.rhoQ, Finset.prod_range_succ, Finset.prod_range_zero]
  have h1 : (1 : ℚ) ≤ n := by exact_mod_cast hn
  have a : (n : ℚ) ≠ 0 := by positivity
  have b : 2 * (n : ℚ) - 1 ≠ 0 := by linarith
  have c : 2 * (n : ℚ) + 1 ≠ 0 := by linarith
  have d : 2 * (n : ℚ) + 2 ≠ 0 := by linarith
  have f : 2 * (n : ℚ) + 3 ≠ 0 := by linarith
  have g : 2 * (n : ℚ) ≠ 0 := by linarith
  push_cast
  field_simp
  ring

lemma z21 (n : ℕ) (hn : 1 ≤ n) : (Q04.Pf 2 n).coeff (Finsupp.single 2 1) =
    1 / ((2 * (n : ℚ) - 1) * (2 * (n : ℚ) + 1)) := by
  have d0 := fun p q => (D21 p q).1
  have d1 := fun p q => (D21 p q).2.1
  have d2 := fun p q => (D21 p q).2.2
  rw [R01.proof_2 2 n hn, MvPolynomial.coeff_smul, smul_eq_mul, R01.PhiN_Ecoef]
  simp only [MvPolynomial.coeff_sum, MvPolynomial.coeff_add, MvPolynomial.coeff_smul, smul_eq_mul, Finset.sum_range_succ,
    Finset.sum_range_zero, zero_add, d0, d1, d2, mul_zero, add_zero, zero_add]
  have e : ∀ j, j ≤ 2 → Q04.eF 2 n j 1 1 = cf2 n j 1 := fun j hj => eF2 n j 1 hj (by norm_num)
  rw [e 0 (by norm_num), e 1 (by norm_num), e 2 (by norm_num)]
  simp only [cf2, V05.rhoQ, Finset.prod_range_succ, Finset.prod_range_zero]
  have h1 : (1 : ℚ) ≤ n := by exact_mod_cast hn
  have a : (n : ℚ) ≠ 0 := by positivity
  have b : 2 * (n : ℚ) - 1 ≠ 0 := by linarith
  have c : 2 * (n : ℚ) + 1 ≠ 0 := by linarith
  have g : 2 * (n : ℚ) ≠ 0 := by linarith
  push_cast
  field_simp
  ring

lemma z4 (n : ℕ) (hn : 1 ≤ n) : (Q04.Pf 4 n).coeff (Finsupp.single 4 1) =
    ((n : ℚ) + 3) / ((2 * (n : ℚ) - 1) * (2 * (n : ℚ) + 1) * (2 * (n : ℚ) + 3)) := by
  rw [R01.proof_5 (Finsupp.single 4 1) (R01.mono_single 4 1 (by norm_num))
    (R01.single_ne 4 1 one_ne_zero) 4 n hn, R01.wt_single,
    show 1 * 4 / 2 = 2 from rfl, R01.Icc12, Finset.sum_insert (by decide), Finset.sum_singleton,
    R01.gam4.1, R01.gam4.2, if_pos (by norm_num), if_pos (by norm_num),
    show 4 + 2 * 1 - 1 * 4 = 2 from rfl, show 4 + 2 * 2 - 1 * 4 = 4 from rfl, z21 n hn, z22 n hn]
  have h1 : (1 : ℚ) ≤ n := by exact_mod_cast hn
  have b : 2 * (n : ℚ) - 1 ≠ 0 := by linarith
  have c : 2 * (n : ℚ) + 1 ≠ 0 := by linarith
  have f : 2 * (n : ℚ) + 3 ≠ 0 := by linarith
  field_simp
  norm_num [Nat.factorial]
  ring

lemma wz5 (n : ℕ) : (Q04.Pf 4 n).coeff (Finsupp.single 5 1) = 0 := by
  by_contra hne
  have h := (Q04.Pf_weight 4 n _ (MvPolynomial.mem_support_iff.mpr hne)).2
  rw [Q04.wsum_single] at h; omega

lemma wz23 (n : ℕ) : (Q04.Pf 4 n).coeff (Finsupp.single 2 1 + Finsupp.single 3 1) = 0 := by
  by_contra hne
  have h := (Q04.Pf_weight 4 n _ (MvPolynomial.mem_support_iff.mpr hne)).2
  rw [Finsupp.sum_add_index' (fun _ => by simp) (fun _ _ _ => by ring), Q04.wsum_single,
    Q04.wsum_single] at h
  omega

lemma zMnat (s : ℕ →₀ ℕ) (k : ℕ) (α : ℝ) : R01.zM s (k : ℤ) α = ∑' n, R01.cser s k α n :=
  R01.zM_nat s k α

end S04

theorem S04.proof_1 : S04.stmt_1 := by
  intro n _
  have h := Q04.Pf_weight 4 n
  conv_lhs => rw [(Q04.Pf 4 n).as_sum]
  rw [Finset.sum_subset (fun s hs => Q04.cases5 s (h s hs).1 ((h s hs).2.trans (by norm_num)))
    (fun s _ hs => by rw [MvPolynomial.notMem_support_iff.mp hs]; simp)]
  repeat rw [Finset.sum_insert (by
    simp only [Finset.mem_insert, Finset.mem_singleton, not_or]
    repeat' apply And.intro
    all_goals exact Q04.fs_ne (by simp [Finsupp.single_apply]))]
  rw [Finset.sum_singleton, S04.wz5, S04.wz23]
  simp only [MvPolynomial.monomial_zero, add_zero, ← MvPolynomial.C_mul_X_pow_eq_monomial,
    ← MvPolynomial.C_apply, pow_one]
  ring

theorem S04.proof_2 : S04.stmt_2 := by
  intro α hα0 hα1
  have hA := R01.allSummable hα0 hα1
  refine ⟨hA, ?_⟩
  have hM := (Q04.proof_5.2 α hα0 hα1 4).2
  have hS : ∀ s : ℕ →₀ ℕ, HasSum (R01.cser s 4 α) (R01.zM s (4 : ℕ) α) := fun s => by
    rw [S04.zMnat]; exact (R01.summable_cser hα0 hα1 s 4).hasSum
  have hT := ((((hS 0).add ((hS (Finsupp.single 2 1)).mul_left (X05.z 2))).add
    ((hS (Finsupp.single 3 1)).mul_left (X05.z 3))).add
    ((hS (Finsupp.single 4 1)).mul_left (X05.z 4))).add
    ((hS (Finsupp.single 2 2)).mul_left (X05.z 2 ^ 2))
  refine hM.unique (hT.congr_fun fun n => ?_)
  rw [Q04.aeval_S5 _ (fun s hs => ⟨(Q04.Pf_weight 4 (n + 1) s hs).1,
    (Q04.Pf_weight 4 (n + 1) s hs).2.trans (by norm_num)⟩), S04.wz5, S04.wz23]
  simp only [R01.cser, Pi.add_apply]
  push_cast; ring

theorem S04.proof_3 : S04.stmt_3 := fun n hn => ⟨S04.z4 n hn, S04.z22 n hn⟩

theorem S04.proof_4 : S04.stmt_4 := by
  intro c hc0 hc1
  have hc : |c| < 1 := by rw [abs_of_pos hc0]; exact hc1
  have h := V05.artanh_hasSum hc
  have hcne : c ≠ 0 := hc0.ne'
  refine ⟨?_, ?_, ?_⟩
  · refine (h.mul_left c).congr_fun fun n => ?_
    rw [show 2 * (n + 1) = 2 * n + 1 + 1 by ring, pow_succ]
    push_cast; ring
  · have h1 := ((hasSum_nat_add_iff' 1).mpr h).mul_left (1 / c)
    convert h1 using 1
    · funext n
      push_cast
      field_simp
      ring
    · simp only [Finset.sum_range_one]
      field_simp
      ring
  · have h1 := ((hasSum_nat_add_iff' 2).mpr h).mul_left (1 / c ^ 3)
    convert h1 using 1
    · funext n
      push_cast
      field_simp
      ring
    · simp only [Finset.sum_range_succ, Finset.sum_range_zero]
      field_simp
      ring

theorem S04.proof_5 : S04.stmt_5 := by
  intro c hc0 hc1
  obtain ⟨h1, h2, h3⟩ := S04.proof_4 c hc0 hc1
  have hc : c ≠ 0 := hc0.ne'
  refine ⟨?_, ?_⟩
  · have h := ((h1.mul_left (7 / 16)).sub (h2.mul_left (5 / 8))).add (h3.mul_left (3 / 16))
    convert h using 1
    · funext n
      rw [S04.z4 (n + 1) (by omega)]
      push_cast
      have := n.cast_nonneg (α := ℝ)
      rw [show 2 * ((n : ℝ) + 1) - 1 = 2 * n + 1 by ring, show 2 * ((n : ℝ) + 1) + 1 = 2 * n + 3 by ring,
        show 2 * ((n : ℝ) + 1) + 3 = 2 * n + 5 by ring]
      have e1 : 2 * (n : ℝ) + 1 ≠ 0 := by linarith
      have e2 : 2 * (n : ℝ) + 3 ≠ 0 := by linarith
      have e3 : 2 * (n : ℝ) + 5 ≠ 0 := by linarith
      field_simp
      ring
    · unfold S04.F1
      field_simp
      ring
  · have h := ((h1.mul_left (-1 / 16)).sub (h2.mul_left (1 / 8))).add (h3.mul_left (3 / 16))
    convert h using 1
    · funext n
      rw [S04.z22 (n + 1) (by omega)]
      push_cast
      have := n.cast_nonneg (α := ℝ)
      rw [show 2 * ((n : ℝ) + 1) - 1 = 2 * n + 1 by ring, show 2 * ((n : ℝ) + 1) + 1 = 2 * n + 3 by ring,
        show 2 * ((n : ℝ) + 1) + 3 = 2 * n + 5 by ring]
      have e1 : 2 * (n : ℝ) + 1 ≠ 0 := by linarith
      have e2 : 2 * (n : ℝ) + 3 ≠ 0 := by linarith
      have e3 : 2 * (n : ℝ) + 5 ≠ 0 := by linarith
      field_simp
      ring
    · unfold S04.F2
      field_simp
      ring

theorem S04.proof_6 : S04.stmt_6 := by
  intro α hα0 hα1
  refine ⟨R01.allSummable hα0 hα1, ?_⟩
  obtain ⟨hc0, hc1, -⟩ := D04.cs_facts hα0 hα1
  obtain ⟨g1, g2⟩ := S04.proof_5 _ hc0 hc1
  have e1 : R01.zM (Finsupp.single 4 1) 4 α = S04.F1 (D04.c α) := (S04.zMnat _ 4 α).trans g1.tsum_eq
  have e2 : R01.zM (Finsupp.single 2 2) 4 α = S04.F2 (D04.c α) := (S04.zMnat _ 4 α).trans g2.tsum_eq
  have hc : D04.c α ≠ 0 := hc0.ne'
  rw [e1, e2, X05.z_four, X05.z_two]
  unfold S04.F1 S04.F2
  field_simp
  ring

theorem S04.proof_7 : S04.stmt_7 := by
  intro α hα0 hα1
  obtain ⟨hA, h2⟩ := S04.proof_2 α hα0 hα1
  obtain ⟨-, h6⟩ := S04.proof_6 α hα0 hα1
  have h3 := (R01.proof_7 α hα0 hα1).2 3
  have h3' : R01.zM (Finsupp.single 3 1) 4 α = -R01.zM (Finsupp.single 2 1) 3 α := h3
  refine ⟨hA, ?_⟩
  rw [h2, h3']
  linear_combination h6
