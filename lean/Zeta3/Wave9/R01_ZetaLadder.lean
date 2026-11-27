import Mathlib
import Zeta3.Wave8.Q04_FlatAllK
/-!
# R01 -- the zeta ladder: every zeta part of `M_k(α)` is a fixed rational combination of the
# `ζ(2)^e` parts

Target id: R01.
Paper label: `thm:zladder` (paper/ieee/main.tex, lines 1591-1593; context 1589).
Supplement: paper/ieee/proofs/sec_g6_ladder.tex (meaning of `[z]M` lines 8-12; `g6:thm:zladder`
15-17; `P_w`, `Q_z`, `γ_e(z)` and `g6:eq:Qz` 19-27; Lemma `g6:thm:ladderprecise` 29-35; `Φ_n`,
`Φ^0_n` 44-49; Lemma `g6:lem:op` 51-53; Lemma `g6:lem:DE` 70-72; proof of the precise form 80-95;
the examples and further rungs 99-119).

Earlier Lean, imported and not restated: `Zeta3.Wave8.Q04_FlatAllK` (which imports
`Zeta3.Wave7.V05_FlatAngle`). From its frozen statement block we use the kernel coefficients
`E_{n,k} = Q04.Ecoef k n ∈ ℚ[y,a,b]` (`y = X 0`, `a = X 1`, `b = X 2` of `MvPolynomial (Fin 3) ℚ`),
`D^♯_{i,m}(p,q) = Q04.DsharpF i m p q ∈ ℚ[𝔷_2, 𝔷_3, …]` and the formal coefficients
`P_{k,n} = Q04.Pf k n ∈ MvPolynomial ℕ ℚ` of Definition `g6:def:formal` (the variable `X w` is the
formal zeta value `𝔷_w`), and `c = D04.c α = cos(α/2)`. `Q04.stmt_5` (proved) says
`M_k(α) = Σ_{n≥1} P_{k,n}(ζ) c^{2n}`, and `Q04.stmt_9` (proved) is the first rung in the formal
sense; stmt_7 below is the same rung in the supplement's notation `[z]M_k`.

Quoted verbatim (main.tex 1589, 1591-1593):

  As in Theorem~\ref{thm:flat}, we treat the zeta values as symbols, and we write $[z]M$ for the
  coefficient of the product $z$ of zeta values in $M$.

  \begin{theorem}[The zeta ladder$^\dagger$]\label{thm:zladder}
  For every $k\ge0$ and $0<\alpha<\pi$, $[\zeta(3)]M_{k+1}(\alpha)=-[\zeta(2)]M_k(\alpha)$. More
  generally, for every product $z$ of zeta values, $[z]M_k(\alpha)$ is a fixed combination, with
  rational coefficients, of the parts $[\zeta(2)^e]M_{k'}(\alpha)$ with $k'\le k$. For example,
  $[\zeta(5)]M_k=-[\zeta(2)]M_{k-3}-2[\zeta(2)^2]M_{k-1}$.
  \end{theorem}

Quoted verbatim (supplement, line 10, the meaning of `[z]M`):

  [z]M_k(\alpha):=\sum_{n\ge1}\bigl([\gvs z]P_{k,n}\bigr)c^{2n}\ ,\qquad [z]M_k:=0\ \text{ for }k<0 .

Quoted verbatim (supplement, line 19): "``Fixed'' means: the rational coefficients depend on $z$
only, not on $k$, $\alpha$ or $n$."

Quoted verbatim (supplement, lines 21-27):

  \mathsf P_w(x,x')=\frac{(-1)^{w+1}}w\Bigl[(x+x')^w-x^w-x'^w\Bigr],\qquad
  \mathsf Q_z=\prod_{w}\frac{\mathsf P_w^{\,e_w}}{e_w!}\quad(z\leftrightarrow\gvs z=\textstyle\prod\gvs w^{e_w}).
  ...
  \mathsf Q_z=\sum_{e=1}^{\lfloor\gwt(z)/2\rfloor}\gamma_e(z)\,\sigma_1^{\gwt(z)-2e}\sigma_2^{\,e},\qquad\gamma_e(z)\in\QQ .

Quoted verbatim (supplement, Lemma `g6:thm:ladderprecise`, lines 30-34):

  For every product $z\ne1$ of zeta values, every $k\ge0$ and $0<\alpha<\pi$,
  [z]M_k(\alpha)=\sum_{e=1}^{\lfloor\gwt(z)/2\rfloor}(-1)^e\,e!\,\gamma_e(z)\,[\zeta(2)^e]M_{k-\gwt(z)+2e}(\alpha),
  and in fact $[\gvs z]P_{k,n}=\sum_e(-1)^ee!\gamma_e(z)[\gvs2^e]P_{k-\gwt(z)+2e,n}$ for every $n\ge1$.

Quoted verbatim (supplement, lines 44-47): "Let $\Phi_n:\QQ[y,a,b]\to\QQ[\gvs2,\gvs3,\dots]$ be the
$\QQ$-linear map with $\Phi_n[y^ja^ib^m]=\rho_j(n)\,D^\sharp_{i,m}(n+j,n)$, and let $\Phi^0_n$ be the
composition of $\Phi_n$ with $\gvs w\mapsto0$ (the part free of zeta values). ... By
Definition~\ref{g6:def:formal}, $P_{k,n}=-\frac1{n(2n-1)}\Phi_n[E_{n,k}]$".

Quoted verbatim (supplement, Lemma `g6:lem:op`, line 52): "For every $F\in\QQ[y,a,b]$ and every
product $z$ of zeta values (including $z=1$, with $\mathsf Q_1=1$):
$[\gvs z]\Phi_n[F]=\Phi^0_n\bigl[\mathsf Q_z(\partial_a,\partial_b)F\bigr]$."

Quoted verbatim (supplement, Lemma `g6:lem:DE`, line 71): "In $\QQ[y,a,b][[u]]$,
$(\partial_a+\partial_b)E_n=u\,E_n$."

Quoted verbatim (supplement, lines 109, 112-118, the examples and further rungs):

  For $z=\zeta(5)$: $W=5$, $\gamma_1=1$, $\gamma_2=-1$, so
  $[\zeta(5)]M_k=(-1)\cdot1\cdot[\zeta(2)]M_{k-3}+(+1)\cdot2\cdot(-1)[\zeta(2)^2]M_{k-1}=-[\zeta(2)]M_{k-3}-2[\zeta(2)^2]M_{k-1}$.
  ...
  [\zeta(4)]M_k&=M^{(1)}_{k-2}+M^{(2)}_k, &\qquad [\zeta(6)]M_k&=M^{(1)}_{k-4}+3M^{(2)}_{k-2}+2M^{(3)}_k,\\
  [\zeta(2)\zeta(3)]M_k&=-2M^{(2)}_{k-1}, &\qquad [\zeta(3)^2]M_k&=M^{(2)}_{k-2}
  (with $M_{k}^{(e)}=[\zeta(2)^e]M_k$).

## Formalization choices

* Everything is in the supplement's formal sense (coefficient extraction on the explicit
  polynomials `P_{k,n} = Q04.Pf k n`, never on real numbers). A product
  `z = ∏ ζ(w)^{e_w}` is a finitely supported `s : ℕ →₀ ℕ` (`s w = e_w`) with all `w ≥ 2`
  (`IsZetaMono s`); `𝔷_z` is the monomial `s`, `[𝔷_z]P = P.coeff s`, and `wt s = Σ e_w w` (the
  same weight as in `Q04.WeightLE`). `z = 1` is `s = 0`; `ζ(2)^e` is `Finsupp.single 2 e`.
* `[z]M_k(α)` is `zM s k α` with an integer index `k`: `Σ_{n≥1} [𝔷_z]P_{k,n} c^{2n}` (written
  with `tsum`) for `k ≥ 0`, and `0` for `k < 0`, exactly as on supplement line 10. Because the
  series is written with `tsum`, every statement that uses `zM` also asserts
  `AllSummable α`: every series `Σ_n |[𝔷_z]P_{k,n} c^{2n}|` converges for `0 < α < π` (the
  supplement's remark after line 10). So no identity below can hold by the junk value `0`.
* `ρ_j(n) = ∏_{t<j} (n+t)/(2n+t)` (Section `g6:sec:flat`) is defined again here (`rhoR`) as an
  audited definition, since `V05.rhoQ` (used inside `Q04.Pf`) sits below V05's END marker; the two
  are the same formula.
* `Φ_n` is `PhiN n`, the `ℚ`-linear map on `ℚ[y,a,b]` given on monomials by
  `y^j a^i b^m ↦ ρ_j(n) D^♯_{i,m}(n+j,n)` (a sum over the support of `F`). `Φ^0_n[F]` (`𝔷_w ↦ 0`) is the
  constant coefficient `(PhiN n F).coeff 0`.
* `P_w`, `Q_z` are polynomials in `x = X 0`, `x' = X 1` of `MvPolynomial (Fin 2) ℚ`;
  `Q_z = ∏_{w ∈ supp s} P_w^{e_w}/e_w!`. The operator `Q(∂_a,∂_b)` (`QOp`) sends `F` to
  `Σ_t q_t ∂_a^{t_0} ∂_b^{t_1} F` with `∂_a = pderiv 1`, `∂_b = pderiv 2` on `ℚ[y,a,b]`.
* `γ_e(z)`: the supplement defines it by `g6:eq:Qz`, the expansion of the symmetric polynomial
  `Q_z` in `σ_1 = x+x'`, `σ_2 = xx'` (unique, since `σ_1, σ_2` are algebraically independent),
  and computes it with the Newton recurrence `p_w = σ_1 p_{w-1} - σ_2 p_{w-2}`, `p_0 = 2`,
  `p_1 = σ_1` (lines 100-108). Here `γ_e(z) = gam e s` is defined by that computation: `QsymN s` is
  `∏ (P_w written in σ_1, σ_2)^{e_w}/e_w!` with `P_w = (-1)^{w+1}(σ_1^w - p_w)/w` (variables
  `σ_1 = X 0`, `σ_2 = X 1`), and `gam e s` is its coefficient of `σ_1^{W-2e} σ_2^e`
  (only used for `1 ≤ e ≤ W/2`, so the truncated subtraction `W - 2e` is exact). stmt_4 states
  that with these numbers `Q_z = Σ_{e=1}^{⌊W/2⌋} γ_e(z) σ_1^{W-2e} σ_2^e` as polynomials in
  `x, x'`, which is `g6:eq:Qz` and identifies `gam` with the supplement's `γ_e(z)`.
* Lemma `g6:lem:DE` lives in `ℚ[y,a,b][[u]]`; it is stated coefficientwise in `u` on the
  coefficients `E_{n,k}` (which do not depend on the truncation used to compute them):
  `(∂_a+∂_b)E_{n,0} = 0` and `(∂_a+∂_b)E_{n,k+1} = E_{n,k}`. It is stated for every `n ∈ ℕ` (the
  identity `(∂_a+∂_b)A = 1` holds for every `n`).
* Terms with `k' = k - W + 2e < 0` are absent (supplement line 94, "with the convention that
  $P_{k',n}=0$ for $k'<0$"): in the formal statements they are written
  `if W ≤ k + 2e then [𝔷_2^e]P_{k+2e-W,n} else 0`; in the series statements `zM` is `0` at negative
  indices.
* The final statement (stmt_10) is the theorem as the paper states it: `AllSummable`; the first
  rung for every `k ≥ 0`; for every product `z` (`z = 1` included) a finite family of rationals
  `r_{e,d}` chosen before `k`, `n`, `α` with `[z]M_k = Σ_{e,d} r_{e,d} [ζ(2)^e]M_{k-d}` (so
  `k' = k - d ≤ k`), both for the formal coefficients `[𝔷_z]P_{k,n}` and for the series; and the
  `ζ(5)` example for every `k ≥ 0`.

## Statements

* stmt_1: Lemma `g6:lem:DE`: `(∂_a+∂_b)E_{n,0} = 0`, `(∂_a+∂_b)E_{n,k+1} = E_{n,k}`.
* stmt_2: `P_{k,n} = -Φ_n[E_{n,k}]/(n(2n-1))` (supplement line 46).
* stmt_3: Lemma `g6:lem:op`: `[𝔷_z]Φ_n[F] = Φ^0_n[Q_z(∂_a,∂_b)F]` for every `F`, every `z`, `n ≥ 1`.
* stmt_4: `g6:eq:Qz`: `Q_z = Σ_{e=1}^{⌊W/2⌋} γ_e(z) (x+x')^{W-2e} (xx')^e` for `z ≠ 1`.
* stmt_5: Lemma `g6:thm:ladderprecise`, formal: `[𝔷_z]P_{k,n} = Σ_e (-1)^e e! γ_e(z) [𝔷_2^e]P_{k-W+2e,n}`.
* stmt_6: Lemma `g6:thm:ladderprecise`, series: the same for `[z]M_k(α)`, with convergence.
* stmt_7: the first rung `[ζ(3)]M_{k+1}(α) = -[ζ(2)]M_k(α)`, `k ≥ 0` (agrees with `Q04.stmt_9`).
* stmt_8: the `ζ(5)` example, formal and series, every `k ≥ 0`.
* stmt_9: the further printed rungs for `ζ(4)`, `ζ(6)`, `ζ(2)ζ(3)`, `ζ(3)²`, every `k ≥ 0`.
* stmt_10: `thm:zladder` as stated in the paper (see the last formalization bullet).

## Numerical check before freezing

A SymPy script mirroring these definitions exactly (the same truncations as `Q04.Eser`,
`Q04.DeltaT`, `Q04.Pf`; exact rationals) gave: stmt_1 for `k ≤ 4`, `n ≤ 3` (`n = 0` included;
`E_{n,0} = 1`); stmt_2 for `k ≤ 4`, `n ≤ 3`; stmt_3 for `n ≤ 2`, all `F = y^j a^i b^m` with
`j ≤ 2`, `i + m ≤ 3`, and all monomials `z` of weight `≤ i+m+1` (204 identities); stmt_4 for all
14 monomials `z ≠ 1` of weight `≤ 7`; stmt_5 for those 14 monomials, `k ≤ 6`, `n ≤ 3` (294
identities, 0 mismatches). The computed constants agree with the supplement: `γ(ζ(3)) = (1)`,
`γ(ζ(5)) = (1,-1)`, `γ(ζ(4)) = (-1,1/2)`, `γ(ζ(6)) = (-1,3/2,-1/3)`, `γ(ζ(2)ζ(3)) = (0,-1)`,
`γ(ζ(3)²) = (0,1/2,0)`. Non-vacuity: `P_{3,2} = 737/51840 + 43𝔷2/900 - 𝔷3/15`,
`P_{2,2} = -103/864 + 𝔷2/15`; `[𝔷5]P_{5,1} = -1/5 = -[𝔷2]P_{2,1} - 2[𝔷2²]P_{4,1}
= -1/3 + 2/15`.

## NOT COVERED

* Remark `g6:rem:kernel` and the proof-idea formula of the paper
  `M(u;α) = Σ_e T_e (-∂_s)^e [Γ(1+ξ)Γ(1+η)/Γ(1+u)]_{s=0}` (a formal series in two variables
  with a real-analytic kernel); only its coefficientwise content (stmt_5, stmt_6) is stated.
* The intermediate objects `T_e(u)`, `Θ_e` of the proof (they are coefficient series of
  `Σ_k P_{k,n} u^k`; their content is in stmt_5).
* The checks paragraph (3150 single-zeta rungs for `k ≤ 8`, `n ≤ 150`; eight weight-nine rungs) and
  the `M_5(α)` display after the theorem (it needs `Ψ_3`, `G`, `F_2`, `R_5`).
* Real-number zeta parts: as the supplement says (line 12), `[z]` of the real number `M_k(α)` is
  not well defined; nothing is stated about it.
-/
-- BEGIN STATEMENTS
namespace R01

open Real Filter MvPolynomial

/-- `ρ_j(n) = ∏_{t<j} (n+t)/(2n+t)` (Section `g6:sec:flat`; the same formula as `V05.rhoQ`). -/
noncomputable def rhoR (j n : ℕ) : ℚ :=
  ∏ t ∈ Finset.range j, ((n : ℚ) + t) / (2 * (n : ℚ) + t)

/-- `Φ_n : ℚ[y,a,b] → ℚ[𝔷_2, 𝔷_3, …]`, linear, `y^j a^i b^m ↦ ρ_j(n) D^♯_{i,m}(n+j,n)`. -/
noncomputable def PhiN (n : ℕ) (F : MvPolynomial (Fin 3) ℚ) : MvPolynomial ℕ ℚ :=
  ∑ t ∈ F.support, (F.coeff t * rhoR (t 0) n) • Q04.DsharpF (t 1) (t 2) (n + t 0) n

/-- The weight `Σ_w e_w w` of the monomial `∏ 𝔷_w^{e_w}`. -/
def wt (s : ℕ →₀ ℕ) : ℕ := s.sum fun w e => e * w

/-- `s` is a product of zeta values: only `𝔷_w` with `w ≥ 2` occur. -/
def IsZetaMono (s : ℕ →₀ ℕ) : Prop := ∀ w ∈ s.support, 2 ≤ w

/-- `P_w(x,x') = (-1)^{w+1}/w [(x+x')^w - x^w - x'^w]`, `x = X 0`, `x' = X 1`. -/
noncomputable def Pw (w : ℕ) : MvPolynomial (Fin 2) ℚ :=
  ((-1 : ℚ) ^ (w + 1) / (w : ℚ)) • ((X 0 + X 1) ^ w - X 0 ^ w - X 1 ^ w)

/-- `Q_z = ∏_w P_w^{e_w}/e_w!` for `𝔷_z = ∏ 𝔷_w^{e_w}` (`Q_1 = 1`). -/
noncomputable def Qz (s : ℕ →₀ ℕ) : MvPolynomial (Fin 2) ℚ :=
  ∏ w ∈ s.support, (((s w).factorial : ℚ)⁻¹) • Pw w ^ s w

/-- The operator `Q(∂_a,∂_b)` on `ℚ[y,a,b]` (`a = X 1`, `b = X 2`): `Σ_t q_t ∂_a^{t_0}∂_b^{t_1}`. -/
noncomputable def QOp (Q : MvPolynomial (Fin 2) ℚ) (F : MvPolynomial (Fin 3) ℚ) :
    MvPolynomial (Fin 3) ℚ :=
  ∑ t ∈ Q.support, Q.coeff t • ((⇑(pderiv (1 : Fin 3)))^[t 0] ((⇑(pderiv (2 : Fin 3)))^[t 1] F))

/-- Newton's power sums `p_w = x^w + x'^w` written in `σ_1 = X 0`, `σ_2 = X 1`:
`p_0 = 2`, `p_1 = σ_1`, `p_{w+2} = σ_1 p_{w+1} - σ_2 p_w`. -/
noncomputable def pN : ℕ → MvPolynomial (Fin 2) ℚ
  | 0 => 2
  | 1 => X 0
  | (w + 2) => X 0 * pN (w + 1) - X 1 * pN w

/-- `Q_z` written in `σ_1 = X 0`, `σ_2 = X 1`, via `P_w = (-1)^{w+1}(σ_1^w - p_w)/w`. -/
noncomputable def QsymN (s : ℕ →₀ ℕ) : MvPolynomial (Fin 2) ℚ :=
  ∏ w ∈ s.support, (((s w).factorial : ℚ)⁻¹) •
    (((-1 : ℚ) ^ (w + 1) / (w : ℚ)) • (X 0 ^ w - pN w)) ^ s w

/-- `γ_e(z)`: the coefficient of `σ_1^{W-2e} σ_2^e` in `Q_z` (`W = wt z`, used for `1 ≤ e ≤ W/2`). -/
noncomputable def gam (e : ℕ) (s : ℕ →₀ ℕ) : ℚ :=
  (QsymN s).coeff (Finsupp.single 0 (wt s - 2 * e) + Finsupp.single 1 e)

/-- The `n`-th term (`n ≥ 0` standing for `n + 1`) of `[z]M_k(α) = Σ_{n≥1} [𝔷_z]P_{k,n} c^{2n}`. -/
noncomputable def cser (s : ℕ →₀ ℕ) (k : ℕ) (α : ℝ) (n : ℕ) : ℝ :=
  (((Q04.Pf k (n + 1)).coeff s : ℚ) : ℝ) * D04.c α ^ (2 * (n + 1))

/-- `[z]M_k(α)` (supplement line 10): the series for `k ≥ 0`, and `0` for `k < 0`. -/
noncomputable def zM (s : ℕ →₀ ℕ) (k : ℤ) (α : ℝ) : ℝ :=
  if 0 ≤ k then ∑' n : ℕ, cser s k.toNat α n else 0

/-- Every series `[z]M_k(α)` converges absolutely. -/
def AllSummable (α : ℝ) : Prop :=
  ∀ (s : ℕ →₀ ℕ) (k : ℕ), Summable (fun n : ℕ => |cser s k α n|)

/-- Lemma `g6:lem:DE`, coefficientwise in `u`: `(∂_a+∂_b)E_n = uE_n` in `ℚ[y,a,b][[u]]`, i.e.
`(∂_a+∂_b)E_{n,0} = 0` and `(∂_a+∂_b)E_{n,k+1} = E_{n,k}`. -/
def stmt_1 : Prop :=
  ∀ n : ℕ,
    (pderiv (1 : Fin 3) (Q04.Ecoef 0 n) + pderiv (2 : Fin 3) (Q04.Ecoef 0 n) = 0) ∧
    ∀ k : ℕ, pderiv (1 : Fin 3) (Q04.Ecoef (k + 1) n) + pderiv (2 : Fin 3) (Q04.Ecoef (k + 1) n)
      = Q04.Ecoef k n

/-- Supplement line 46: "By Definition g6:def:formal, `P_{k,n} = -Φ_n[E_{n,k}]/(n(2n-1))`". -/
def stmt_2 : Prop :=
  ∀ k n : ℕ, 1 ≤ n →
    Q04.Pf k n = (-(1 / ((n : ℚ) * (2 * (n : ℚ) - 1)))) • PhiN n (Q04.Ecoef k n)

/-- Lemma `g6:lem:op`: for every `F ∈ ℚ[y,a,b]` and every product `z` of zeta values
(`z = 1` included), `[𝔷_z]Φ_n[F] = Φ^0_n[Q_z(∂_a,∂_b)F]`. -/
def stmt_3 : Prop :=
  ∀ n : ℕ, 1 ≤ n → ∀ (F : MvPolynomial (Fin 3) ℚ) (s : ℕ →₀ ℕ), IsZetaMono s →
    (PhiN n F).coeff s = (PhiN n (QOp (Qz s) F)).coeff 0

/-- `g6:eq:Qz`: for a product `z ≠ 1` of zeta values of weight `W`,
`Q_z = Σ_{e=1}^{⌊W/2⌋} γ_e(z) (x+x')^{W-2e} (xx')^e`. -/
def stmt_4 : Prop :=
  ∀ s : ℕ →₀ ℕ, IsZetaMono s → s ≠ 0 →
    Qz s = ∑ e ∈ Finset.Icc 1 (wt s / 2),
      C (gam e s) * (X 0 + X 1) ^ (wt s - 2 * e) * (X 0 * X 1) ^ e

/-- Lemma `g6:thm:ladderprecise`, formal form: for every product `z ≠ 1` of zeta values of weight
`W`, every `k ≥ 0`, `n ≥ 1`,
`[𝔷_z]P_{k,n} = Σ_{e=1}^{⌊W/2⌋} (-1)^e e! γ_e(z) [𝔷_2^e]P_{k-W+2e,n}` (terms with
`k - W + 2e < 0` absent). -/
def stmt_5 : Prop :=
  ∀ s : ℕ →₀ ℕ, IsZetaMono s → s ≠ 0 → ∀ k n : ℕ, 1 ≤ n →
    (Q04.Pf k n).coeff s = ∑ e ∈ Finset.Icc 1 (wt s / 2),
      (-1 : ℚ) ^ e * (e.factorial : ℚ) * gam e s *
        (if wt s ≤ k + 2 * e then (Q04.Pf (k + 2 * e - wt s) n).coeff (Finsupp.single 2 e) else 0)

/-- Lemma `g6:thm:ladderprecise`, series form: for `0 < α < π` all the series `[z]M_k(α)` converge
absolutely, and for every product `z ≠ 1` of zeta values of weight `W` and every `k ≥ 0`,
`[z]M_k(α) = Σ_{e=1}^{⌊W/2⌋} (-1)^e e! γ_e(z) [ζ(2)^e]M_{k-W+2e}(α)`. -/
def stmt_6 : Prop :=
  ∀ α : ℝ, 0 < α → α < π →
    AllSummable α ∧
    ∀ s : ℕ →₀ ℕ, IsZetaMono s → s ≠ 0 → ∀ k : ℕ,
      zM s k α = ∑ e ∈ Finset.Icc 1 (wt s / 2),
        ((-1 : ℝ) ^ e * (e.factorial : ℝ) * (gam e s : ℝ)) *
          zM (Finsupp.single 2 e) ((k : ℤ) - wt s + 2 * e) α

/-- The first rung (first sentence of `thm:zladder`): for every `k ≥ 0` and `0 < α < π`,
`[ζ(3)]M_{k+1}(α) = -[ζ(2)]M_k(α)`, the series converging absolutely. -/
def stmt_7 : Prop :=
  ∀ α : ℝ, 0 < α → α < π →
    AllSummable α ∧
    ∀ k : ℕ, zM (Finsupp.single 3 1) ((k : ℤ) + 1) α = -zM (Finsupp.single 2 1) k α

/-- The example of `thm:zladder`: `[ζ(5)]M_k = -[ζ(2)]M_{k-3} - 2[ζ(2)²]M_{k-1}` for every
`k ≥ 0` (parts with negative index are `0`), for the formal coefficients `[𝔷_5]P_{k,n}` and for the
series at every `0 < α < π`. -/
def stmt_8 : Prop :=
  (∀ k n : ℕ, 1 ≤ n →
    (Q04.Pf k n).coeff (Finsupp.single 5 1) =
      -(if 3 ≤ k then (Q04.Pf (k - 3) n).coeff (Finsupp.single 2 1) else 0)
      - 2 * (if 1 ≤ k then (Q04.Pf (k - 1) n).coeff (Finsupp.single 2 2) else 0)) ∧
  ∀ α : ℝ, 0 < α → α < π →
    AllSummable α ∧
    ∀ k : ℕ, zM (Finsupp.single 5 1) k α =
      -zM (Finsupp.single 2 1) ((k : ℤ) - 3) α - 2 * zM (Finsupp.single 2 2) ((k : ℤ) - 1) α

/-- The further rungs printed in the supplement (lines 114-118), with `M^{(e)}_k = [ζ(2)^e]M_k`:
`[ζ(4)]M_k = M^{(1)}_{k-2} + M^{(2)}_k`, `[ζ(6)]M_k = M^{(1)}_{k-4} + 3M^{(2)}_{k-2} + 2M^{(3)}_k`,
`[ζ(2)ζ(3)]M_k = -2M^{(2)}_{k-1}`, `[ζ(3)²]M_k = M^{(2)}_{k-2}`, for every `k ≥ 0` and
`0 < α < π`. -/
def stmt_9 : Prop :=
  ∀ α : ℝ, 0 < α → α < π →
    AllSummable α ∧
    ∀ k : ℕ,
      zM (Finsupp.single 4 1) k α =
        zM (Finsupp.single 2 1) ((k : ℤ) - 2) α + zM (Finsupp.single 2 2) k α ∧
      zM (Finsupp.single 6 1) k α =
        zM (Finsupp.single 2 1) ((k : ℤ) - 4) α + 3 * zM (Finsupp.single 2 2) ((k : ℤ) - 2) α
          + 2 * zM (Finsupp.single 2 3) k α ∧
      zM (Finsupp.single 2 1 + Finsupp.single 3 1) k α =
        -2 * zM (Finsupp.single 2 2) ((k : ℤ) - 1) α ∧
      zM (Finsupp.single 3 2) k α = zM (Finsupp.single 2 2) ((k : ℤ) - 2) α

/-- `thm:zladder` as stated: for `0 < α < π` all `[z]M_k(α)` converge absolutely; for every
`k ≥ 0`, `[ζ(3)]M_{k+1}(α) = -[ζ(2)]M_k(α)`; for every product `z` of zeta values there are
rationals `r_{e,d}` (finitely many, depending on `z` only, not on `k`, `n`, `α`) with
`[𝔷_z]P_{k,n} = Σ_{e,d} r_{e,d} [𝔷_2^e]P_{k-d,n}` and `[z]M_k(α) = Σ_{e,d} r_{e,d} [ζ(2)^e]M_{k-d}(α)`
(so every `k' = k - d ≤ k`; parts with `k' < 0` are `0`); and
`[ζ(5)]M_k = -[ζ(2)]M_{k-3} - 2[ζ(2)²]M_{k-1}` for every `k ≥ 0`. -/
def stmt_10 : Prop :=
  (∀ α : ℝ, 0 < α → α < π → AllSummable α) ∧
  (∀ k : ℕ, ∀ α : ℝ, 0 < α → α < π →
    zM (Finsupp.single 3 1) ((k : ℤ) + 1) α = -zM (Finsupp.single 2 1) k α) ∧
  (∀ s : ℕ →₀ ℕ, IsZetaMono s → ∃ (N : ℕ) (r : ℕ → ℕ → ℚ), ∀ k : ℕ,
    (∀ n : ℕ, 1 ≤ n → (Q04.Pf k n).coeff s =
      ∑ e ∈ Finset.range N, ∑ d ∈ Finset.range N,
        r e d * (if d ≤ k then (Q04.Pf (k - d) n).coeff (Finsupp.single 2 e) else 0)) ∧
    ∀ α : ℝ, 0 < α → α < π →
      zM s k α = ∑ e ∈ Finset.range N, ∑ d ∈ Finset.range N,
        (r e d : ℝ) * zM (Finsupp.single 2 e) ((k : ℤ) - d) α) ∧
  (∀ k : ℕ, ∀ α : ℝ, 0 < α → α < π →
    zM (Finsupp.single 5 1) k α =
      -zM (Finsupp.single 2 1) ((k : ℤ) - 3) α - 2 * zM (Finsupp.single 2 2) ((k : ℤ) - 1) α)

end R01
-- END STATEMENTS

namespace R01
open MvPolynomial

lemma Ecoef_zero (n : ℕ) : Q04.Ecoef 0 n = 1 := by
  simp [Q04.Ecoef, Q04.Eser, Q04.expT]

lemma rhoR_eq (j n : ℕ) : rhoR j n = V05.rhoQ j n := rfl

/-- the monomial `y^j a^i b^m` -/
noncomputable def mono (j i m : ℕ) : Fin 3 →₀ ℕ := Finsupp.single 0 j + Finsupp.single 1 i + Finsupp.single 2 m

lemma mono_apply (j i m : ℕ) : (mono j i m) 0 = j ∧ (mono j i m) 1 = i ∧ (mono j i m) 2 = m := by
  simp [mono]

lemma mono_eta (t : Fin 3 →₀ ℕ) : mono (t 0) (t 1) (t 2) = t := by
  ext x; fin_cases x <;> simp [mono]

lemma PhiN_Ecoef (k n : ℕ) : PhiN n (Q04.Ecoef k n) =
    ∑ j ∈ Finset.range (k + 1), ∑ i ∈ Finset.range (k + 1),
      (Q04.eF k n j i (k - i) * V05.rhoQ j n) • Q04.DsharpF i (k - i) (n + j) n := by
  classical
  unfold PhiN
  set g : (Fin 3 →₀ ℕ) → MvPolynomial ℕ ℚ := fun t =>
    ((Q04.Ecoef k n).coeff t * rhoR (t 0) n) • Q04.DsharpF (t 1) (t 2) (n + t 0) n with hg
  set f : ℕ × ℕ → (Fin 3 →₀ ℕ) := fun p => mono p.1 p.2 (k - p.2) with hf
  have hinj : Set.InjOn f ↑(Finset.range (k + 1) ×ˢ Finset.range (k + 1)) := by
    intro p _ q _ h
    have h0 := DFunLike.congr_fun h 0
    have h1 := DFunLike.congr_fun h 1
    simp [hf, mono] at h0 h1
    exact Prod.ext h0 h1
  have hsub : (Q04.Ecoef k n).support ⊆ (Finset.range (k + 1) ×ˢ Finset.range (k + 1)).image f := by
    intro t ht
    rw [MvPolynomial.mem_support_iff] at ht
    have ht' : Q04.eF k n (t 0) (t 1) (t 2) ≠ 0 := by
      unfold Q04.eF; rw [← mono, mono_eta]; exact ht
    obtain ⟨h1, h2⟩ := Q04.eF_supp _ _ _ _ _ ht'
    refine Finset.mem_image.mpr ⟨(t 0, t 1), ?_, ?_⟩
    · simp only [Finset.mem_product, Finset.mem_range]; omega
    · rw [hf]; dsimp only
      rw [show k - t 1 = t 2 by omega, mono_eta]
  have e1 : ∑ t ∈ (Q04.Ecoef k n).support, g t =
      ∑ t ∈ (Finset.range (k + 1) ×ˢ Finset.range (k + 1)).image f, g t := by
    refine Finset.sum_subset hsub fun t _ ht => ?_
    simp only [hg]; rw [MvPolynomial.notMem_support_iff.mp ht]; simp
  rw [e1, Finset.sum_image hinj, Finset.sum_product]
  refine Finset.sum_congr rfl fun j _ => Finset.sum_congr rfl fun i _ => ?_
  obtain ⟨a0, a1, a2⟩ := mono_apply j i (k - i)
  simp only [hg, hf, a0, a1, a2, Q04.eF, rhoR_eq]
  rfl

end R01

theorem R01.proof_1 : R01.stmt_1 := by
  intro n
  refine ⟨?_, fun k => ?_⟩
  · rw [R01.Ecoef_zero]; simp
  · have := Q04.Dd_Ecoef n k
    simpa [Q04.Dd] using this

theorem R01.proof_2 : R01.stmt_2 := by
  intro k n _
  rw [R01.PhiN_Ecoef]; rfl

namespace R01

lemma allSummable {α : ℝ} (hα0 : 0 < α) (hα1 : α < Real.pi) : AllSummable α := by
  intro s k
  obtain ⟨hc0, hc1, -, -, -⟩ := D04.cs_facts hα0 hα1
  obtain ⟨C, hC⟩ := Q04.proof_6 k
  have hr0 : 0 ≤ D04.c α ^ 2 := by positivity
  have hr1 : D04.c α ^ 2 < 1 := by nlinarith
  refine (V05.maj_summable k hr0 hr1 (|C| * ((k : ℝ) + 2) ^ k)).of_norm_bounded fun n => ?_
  have hm : (1 : ℝ) ≤ (n : ℝ) + 1 := by have := n.cast_nonneg (α := ℝ); linarith
  have hb := hC (n + 1) (by omega) s
  push_cast at hb
  have hL0 : 0 ≤ Real.log (2 * ((n : ℝ) + 1) + k) := Real.log_nonneg (by
    have := k.cast_nonneg (α := ℝ); linarith)
  have hL1 : 1 + Real.log (2 * ((n : ℝ) + 1) + k) ≤ ((k : ℝ) + 2) * ((n : ℝ) + 1) := by
    have := Real.log_le_sub_one_of_pos (show (0 : ℝ) < 2 * ((n : ℝ) + 1) + k by positivity)
    have := k.cast_nonneg (α := ℝ)
    nlinarith
  have hP : ((n : ℝ) + 1) ^ k * (1 + Real.log (2 * ((n : ℝ) + 1) + k)) ^ k ≤
      ((k : ℝ) + 2) ^ k * ((n : ℝ) + 1) ^ (2 * k) := by
    calc ((n : ℝ) + 1) ^ k * (1 + Real.log (2 * ((n : ℝ) + 1) + k)) ^ k
        ≤ ((n : ℝ) + 1) ^ k * (((k : ℝ) + 2) * ((n : ℝ) + 1)) ^ k :=
          mul_le_mul_of_nonneg_left (pow_le_pow_left₀ (by linarith) hL1 k) (by positivity)
      _ = _ := by rw [mul_pow, two_mul, pow_add]; ring
  have hb2 : |(((Q04.Pf k (n + 1)).coeff s : ℚ) : ℝ)| ≤
      |C| * (((k : ℝ) + 2) ^ k * ((n : ℝ) + 1) ^ (2 * k)) := by
    refine hb.trans ?_
    rw [mul_assoc]
    exact mul_le_mul (le_abs_self C) hP (by positivity) (abs_nonneg C)
  unfold cser
  rw [Real.norm_eq_abs, abs_abs, abs_mul, abs_pow, pow_mul, abs_of_pos hc0]
  calc _ ≤ |C| * (((k : ℝ) + 2) ^ k * ((n : ℝ) + 1) ^ (2 * k)) * (D04.c α ^ 2) ^ (n + 1) :=
        mul_le_mul_of_nonneg_right hb2 (by positivity)
    _ = _ := by ring

lemma summable_cser {α : ℝ} (hα0 : 0 < α) (hα1 : α < Real.pi) (s : ℕ →₀ ℕ) (k : ℕ) :
    Summable (cser s k α) :=
  (allSummable hα0 hα1 s k).of_abs

lemma zM_nat (s : ℕ →₀ ℕ) (k : ℕ) (α : ℝ) : zM s (k : ℤ) α = ∑' n, cser s k α n := by
  simp [zM]

end R01

theorem R01.proof_7 : R01.stmt_7 := by
  intro α hα0 hα1
  refine ⟨R01.allSummable hα0 hα1, fun k => ?_⟩
  rw [show ((k : ℤ) + 1) = ((k + 1 : ℕ) : ℤ) by push_cast; ring, R01.zM_nat, R01.zM_nat,
    ← tsum_neg]
  refine tsum_congr fun n => ?_
  unfold R01.cser
  rw [Q04.stmt9_formal]
  push_cast; ring

namespace R01
open MvPolynomial

/-- `σ_1 ↦ x + x'`, `σ_2 ↦ x x'` -/
noncomputable def sig : MvPolynomial (Fin 2) ℚ →ₐ[ℚ] MvPolynomial (Fin 2) ℚ :=
  aeval ![X 0 + X 1, X 0 * X 1]

lemma sig_X0 : sig (X 0) = X 0 + X 1 := by simp [sig]
lemma sig_X1 : sig (X 1) = X 0 * X 1 := by simp [sig]

lemma sig_pN (w : ℕ) : sig (pN w) = X 0 ^ w + X 1 ^ w ∧ sig (pN (w + 1)) = X 0 ^ (w + 1) + X 1 ^ (w + 1) := by
  induction w with
  | zero =>
    refine ⟨?_, ?_⟩
    · simp only [pN, pow_zero]; rw [map_ofNat]; norm_num
    · simp [pN, sig_X0]
  | succ w ih =>
    refine ⟨ih.2, ?_⟩
    rw [show w + 1 + 1 = w + 2 by ring, pN, map_sub, map_mul, map_mul, ih.2, ih.1, sig_X0, sig_X1]
    ring

/-- weights `σ_1 ↦ 1`, `σ_2 ↦ 2` -/
def wv : Fin 2 → ℕ := ![1, 2]

lemma pN_hom (w : ℕ) : IsWeightedHomogeneous wv (pN w) w ∧ IsWeightedHomogeneous wv (pN (w + 1)) (w + 1) := by
  induction w with
  | zero =>
    refine ⟨?_, ?_⟩
    · show IsWeightedHomogeneous wv (2 : MvPolynomial (Fin 2) ℚ) 0
      rw [show (2 : MvPolynomial (Fin 2) ℚ) = C 2 from (map_ofNat C 2).symm]
      exact isWeightedHomogeneous_C _ _
    · rw [show pN (0 + 1) = X 0 from rfl]; exact isWeightedHomogeneous_X ℚ wv 0
  | succ w ih =>
    refine ⟨ih.2, ?_⟩
    rw [show w + 1 + 1 = w + 2 by ring, pN]
    refine IsWeightedHomogeneous.sub ?_ ?_
    · have := (isWeightedHomogeneous_X ℚ wv 0).mul ih.2
      convert this using 1; simp [wv]; ring
    · have := (isWeightedHomogeneous_X ℚ wv 1).mul ih.1
      convert this using 1; simp [wv]; ring

lemma pN_dvd (w : ℕ) : (X 1 : MvPolynomial (Fin 2) ℚ) ∣ X 0 ^ (w + 1) - pN (w + 1) := by
  induction w with
  | zero => simp [pN]
  | succ w ih =>
    rw [show w + 1 + 1 = w + 2 by ring, pN]
    have : (X 0 : MvPolynomial (Fin 2) ℚ) ^ (w + 2) - (X 0 * pN (w + 1) - X 1 * pN w) =
        X 0 * (X 0 ^ (w + 1) - pN (w + 1)) + X 1 * pN w := by ring
    rw [this]
    exact dvd_add (dvd_mul_of_dvd_right ih _) (dvd_mul_right _ _)

lemma wt_eq (s : ℕ →₀ ℕ) : wt s = ∑ w ∈ s.support, s w * w := rfl

lemma QsymN_hom (s : ℕ →₀ ℕ) : IsWeightedHomogeneous wv (QsymN s) (wt s) := by
  unfold QsymN
  rw [wt_eq]
  refine IsWeightedHomogeneous.prod _ _ _ fun w _ => ?_
  rw [smul_eq_C_mul]
  have h1 : IsWeightedHomogeneous wv (((-1 : ℚ) ^ (w + 1) / (w : ℚ)) •
      ((X 0 : MvPolynomial (Fin 2) ℚ) ^ w - pN w)) w := by
    rw [smul_eq_C_mul]
    refine IsWeightedHomogeneous.C_mul (IsWeightedHomogeneous.sub ?_ (pN_hom w).1) _
    have := (isWeightedHomogeneous_X ℚ wv 0).pow w
    simpa [wv] using this
  have := (h1.pow (s w)).C_mul (((s w).factorial : ℚ)⁻¹)
  simpa [mul_comm] using this

lemma QsymN_dvd (s : ℕ →₀ ℕ) (hs : IsZetaMono s) (h0 : s ≠ 0) :
    (X 1 : MvPolynomial (Fin 2) ℚ) ∣ QsymN s := by
  obtain ⟨w, hw⟩ := Finsupp.ne_iff.mp h0
  have hws : w ∈ s.support := by simpa using hw
  have h2 := hs w hws
  unfold QsymN
  refine dvd_trans ?_ (Finset.dvd_prod_of_mem _ hws)
  rw [smul_eq_C_mul]
  refine dvd_mul_of_dvd_right ?_ _
  have hsw : s w ≠ 0 := by simpa using hw
  obtain ⟨e, he⟩ : ∃ e, s w = e + 1 := ⟨s w - 1, by omega⟩
  rw [he, pow_succ]
  refine dvd_mul_of_dvd_right ?_ _
  rw [smul_eq_C_mul]
  refine dvd_mul_of_dvd_right ?_ _
  obtain ⟨v, hv⟩ : ∃ v, w = v + 1 := ⟨w - 1, by omega⟩
  rw [hv]; exact pN_dvd v

lemma wv_weight (u : Fin 2 →₀ ℕ) : Finsupp.weight wv u = u 0 + 2 * u 1 := by
  rw [Finsupp.weight_apply, Finsupp.sum_fintype _ _ (by simp)]
  simp [Fin.sum_univ_two, wv]; ring

lemma mon2_eq (a b : ℕ) : (C (1 : ℚ) * X 0 ^ a * X 1 ^ b : MvPolynomial (Fin 2) ℚ) =
    monomial (Finsupp.single 0 a + Finsupp.single 1 b) 1 := by
  rw [C_1, one_mul, ← Q04.mon_eq]

lemma mon_C (a b : ℕ) (c : ℚ) : (monomial (Finsupp.single 0 a + Finsupp.single 1 b) c :
    MvPolynomial (Fin 2) ℚ) = C c * (X 0 ^ a * X 1 ^ b) := by
  rw [Q04.mon_eq, C_mul_monomial, mul_one]

lemma expand_sym (F : MvPolynomial (Fin 2) ℚ) (W : ℕ) (hF : IsWeightedHomogeneous wv F W)
    (hd : (X 1 : MvPolynomial (Fin 2) ℚ) ∣ F) :
    F = ∑ e ∈ Finset.Icc 1 (W / 2),
      monomial (Finsupp.single 0 (W - 2 * e) + Finsupp.single 1 e) (F.coeff (Finsupp.single 0 (W - 2 * e) + Finsupp.single 1 e)) := by
  classical
  ext u
  rw [coeff_sum]
  simp only [coeff_monomial]
  by_cases hu : F.coeff u = 0
  · rw [hu]; refine (Finset.sum_eq_zero fun e _ => ?_).symm
    split_ifs with h
    · rw [h, hu]
    · rfl
  · have hw := hF hu
    rw [wv_weight] at hw
    have h1 : 1 ≤ u 1 := by
      obtain ⟨G, hG⟩ := hd
      by_contra hc
      apply hu
      rw [hG, coeff_X_mul']
      rw [if_neg]; simp; omega
    have hu' : Finsupp.single 0 (W - 2 * u 1) + Finsupp.single 1 (u 1) = u := by
      ext x; fin_cases x <;> simp <;> omega
    rw [Finset.sum_eq_single (u 1)]
    · rw [if_pos hu', hu']
    · intro e _ he
      rw [if_neg]
      intro h; apply he
      have := DFunLike.congr_fun h 1; simp at this; exact this
    · intro h; exfalso; apply h; rw [Finset.mem_Icc]; omega

lemma sig_QsymN (s : ℕ →₀ ℕ) : sig (QsymN s) = Qz s := by
  unfold QsymN Qz
  rw [map_prod]
  refine Finset.prod_congr rfl fun w _ => ?_
  rw [map_smul, map_pow, map_smul, map_sub, map_pow, sig_X0, (sig_pN w).1, Pw]
  congr 2; ring

end R01

theorem R01.proof_4 : R01.stmt_4 := by
  intro s hs h0
  show _ = ∑ e ∈ Finset.Icc 1 (R01.wt s / 2), _
  rw [← R01.sig_QsymN, R01.expand_sym _ _ (R01.QsymN_hom s) (R01.QsymN_dvd s hs h0), map_sum]
  refine Finset.sum_congr rfl fun e _ => ?_
  rw [show (R01.QsymN s).coeff (Finsupp.single 0 (R01.wt s - 2 * e) + Finsupp.single 1 e) =
    R01.gam e s from rfl]
  rw [R01.mon_C]
  rw [map_mul, map_mul, map_pow, map_pow, R01.sig_X0, R01.sig_X1, R01.sig, MvPolynomial.aeval_C,
    MvPolynomial.algebraMap_eq]
  ring

namespace R01
open MvPolynomial

abbrev S2 := MvPolynomial (Fin 2) (MvPolynomial ℕ ℚ)
abbrev Q2 := MvPolynomial (Fin 2) ℚ

/-- `∂/∂𝔷_w` applied to each coefficient -/
noncomputable def Dz (w : ℕ) (F : S2) : S2 :=
  ∑ t ∈ F.support, monomial t (pderiv w (F.coeff t))

lemma Dz_coeff (w : ℕ) (F : S2) (t : Fin 2 →₀ ℕ) : (Dz w F).coeff t = pderiv w (F.coeff t) := by
  classical
  unfold Dz
  rw [coeff_sum]
  simp only [coeff_monomial]
  rw [Finset.sum_ite_eq']
  split_ifs with h
  · rfl
  · rw [notMem_support_iff.mp h, map_zero]

lemma Dz_add (w : ℕ) (F G : S2) : Dz w (F + G) = Dz w F + Dz w G := by
  ext t; simp [Dz_coeff]

lemma Dz_smul (w : ℕ) (c : ℚ) (F : S2) : Dz w (c • F) = c • Dz w F := by
  ext t; simp [Dz_coeff, coeff_smul]

lemma Dz_sum (w : ℕ) {ι : Type*} (S : Finset ι) (f : ι → S2) :
    Dz w (∑ i ∈ S, f i) = ∑ i ∈ S, Dz w (f i) := by
  ext t; simp [Dz_coeff, coeff_sum]

lemma Dz_mul (w : ℕ) (F G : S2) : Dz w (F * G) = Dz w F * G + F * Dz w G := by
  refine MvPolynomial.ext _ _ fun t => ?_
  rw [Dz_coeff, MvPolynomial.coeff_add, coeff_mul, coeff_mul, coeff_mul, map_sum, ← Finset.sum_add_distrib]
  refine Finset.sum_congr rfl (fun x _ => ?_)
  rw [Derivation.leibniz, Dz_coeff, Dz_coeff, smul_eq_mul, smul_eq_mul]
  ring

lemma Dz_one (w : ℕ) : Dz w 1 = 0 := by
  ext t; rw [Dz_coeff, coeff_one]; split_ifs <;> simp

lemma Dz_pow (w : ℕ) (F : S2) (r : ℕ) :
    Dz w (F ^ (r + 1)) = ((r : ℚ) + 1) • (F ^ r * Dz w F) := by
  induction r with
  | zero => simp
  | succ r ih =>
    rw [pow_succ, Dz_mul, ih]
    rw [smul_mul_assoc, Algebra.smul_def, Algebra.smul_def]
    push_cast
    try simp only [map_add, map_one, map_natCast]
    ring

lemma Dz_expT (w K : ℕ) (F : S2) : Dz w (Q04.expT (K + 1) F) = Q04.expT K F * Dz w F := by
  unfold Q04.expT
  rw [Finset.sum_range_succ', Dz_add, Dz_sum, pow_zero, Dz_smul, Dz_one, smul_zero, add_zero,
    Finset.sum_mul]
  refine Finset.sum_congr rfl fun r _ => ?_
  rw [Dz_smul, Dz_pow, smul_smul, smul_mul_assoc]
  congr 1
  rw [Nat.factorial_succ]
  push_cast
  have : (r.factorial : ℚ) ≠ 0 := by positivity
  field_simp

/-- `Q2 → S2`, rational coefficients as constants in `𝔷` -/
noncomputable def Lq : Q2 →+* S2 := MvPolynomial.map (C : ℚ →+* MvPolynomial ℕ ℚ)

lemma Lq_coeff (R : Q2) (t : Fin 2 →₀ ℕ) : (Lq R).coeff t = C (R.coeff t) := by
  unfold Lq; rw [coeff_map]

lemma Dz_Lq (w : ℕ) (R : Q2) : Dz w (Lq R) = 0 := by
  ext t; rw [Dz_coeff, Lq_coeff, pderiv_C]; simp

lemma Dz_iter_coeff (w e : ℕ) (F : S2) (t : Fin 2 →₀ ℕ) :
    ((Dz w)^[e] F).coeff t = (⇑(pderiv w))^[e] (F.coeff t) := by
  induction e generalizing F with
  | zero => rfl
  | succ e ih =>
    rw [Function.iterate_succ_apply, Function.iterate_succ_apply, ih, Dz_coeff]

lemma Dz_iter_expT (w : ℕ) (Δ : S2) (R : Q2) (hΔ : Dz w Δ = Lq R) (M e : ℕ) (he : e ≤ M) :
    (Dz w)^[e] (Q04.expT M Δ) = Q04.expT (M - e) Δ * Lq (R ^ e) := by
  induction e with
  | zero => simp
  | succ e ih =>
    rw [Function.iterate_succ_apply', ih (by omega), Dz_mul, Dz_Lq, mul_zero, add_zero,
      show M - e = (M - (e + 1)) + 1 by omega, Dz_expT, hΔ, mul_assoc, ← map_mul, pow_succ,
      mul_comm (R ^ e) R]

/-- `[𝔷^s]` applied to each coefficient -/
noncomputable def cz (s : ℕ →₀ ℕ) (F : S2) : Q2 :=
  ∑ t ∈ F.support, monomial t ((F.coeff t).coeff s)

lemma cz_coeff (s : ℕ →₀ ℕ) (F : S2) (t : Fin 2 →₀ ℕ) : (cz s F).coeff t = (F.coeff t).coeff s := by
  classical
  unfold cz
  rw [coeff_sum]
  simp only [coeff_monomial]
  rw [Finset.sum_ite_eq']
  split_ifs with h
  · rfl
  · rw [notMem_support_iff.mp h]; rfl

lemma cz_zero (F : S2) : cz 0 F = Q04.K0 F := by
  ext t; rw [cz_coeff, Q04.K0_coeff, constantCoeff_eq]

lemma cz_mulL (s : ℕ →₀ ℕ) (F : S2) (R : Q2) : cz s (F * Lq R) = cz s F * R := by
  refine MvPolynomial.ext _ _ fun t => ?_
  rw [cz_coeff, coeff_mul, coeff_mul, coeff_sum]
  refine Finset.sum_congr rfl fun x _ => ?_
  rw [Lq_coeff, mul_comm, coeff_C_mul, cz_coeff, mul_comm]

lemma coeff_iter_pderiv (w e : ℕ) (P : MvPolynomial ℕ ℚ) (f : ℕ →₀ ℕ) (hf : f w = 0) :
    ((⇑(pderiv w))^[e] P).coeff f = (e.factorial : ℚ) * P.coeff (f + Finsupp.single w e) := by
  induction e generalizing P with
  | zero => simp
  | succ e ih =>
    rw [Function.iterate_succ_apply, ih, coeff_pderiv]
    have : (f + Finsupp.single w e) w = e := by simp [hf]
    rw [this, add_assoc, ← Finsupp.single_add, Nat.factorial_succ]
    push_cast; ring

lemma cz_split (w e : ℕ) (f : ℕ →₀ ℕ) (hf : f w = 0) (F : S2) :
    cz (f + Finsupp.single w e) F = ((e.factorial : ℚ)⁻¹) • cz f ((Dz w)^[e] F) := by
  ext t
  rw [coeff_smul, cz_coeff, cz_coeff, Dz_iter_coeff, coeff_iter_pderiv _ _ _ _ hf, smul_eq_mul]
  have : (e.factorial : ℚ) ≠ 0 := by positivity
  field_simp

/-- number of zeta factors -/
def deg (s : ℕ →₀ ℕ) : ℕ := s.sum fun _ e => e

lemma deg_single_add (w e : ℕ) (f : ℕ →₀ ℕ) : deg (Finsupp.single w e + f) = e + deg f := by
  unfold deg
  rw [Finsupp.sum_add_index' (by simp) (by intros; rfl), Finsupp.sum_single_index rfl]

lemma Pw_eq (w : ℕ) : Pw w = Q04.Qb w := rfl

lemma Lq_X (i : Fin 2) : Lq (X i) = X i := by unfold Lq; rw [map_X]

lemma Dz_C (w : ℕ) (a : MvPolynomial ℕ ℚ) : Dz w (C a : S2) = C (pderiv w a) := by
  refine MvPolynomial.ext _ _ fun t => ?_
  rw [Dz_coeff, coeff_C, coeff_C]; split_ifs <;> simp

lemma Delta_split (K p q : ℕ) : Q04.DeltaT K p q = Lq (Q04.Hb K p q) +
    ∑ v ∈ Finset.Icc 2 K, ((-1 : ℚ) ^ (v + 1) / (v : ℚ)) •
      (C (X v) * Lq ((X 0 + X 1) ^ v - X 0 ^ v - X 1 ^ v)) := by
  unfold Q04.DeltaT Q04.Hb
  simp only [map_sum, map_rat_smul, map_add, map_sub, map_pow, Lq_X]

lemma Dz_Delta (K p q w : ℕ) (hw : 2 ≤ w) (hwK : w ≤ K) :
    Dz w (Q04.DeltaT K p q) = Lq (Q04.Qb w) := by
  rw [Delta_split, Dz_add, Dz_Lq, zero_add, Dz_sum]
  simp only [Dz_smul, Dz_mul, Dz_C, Dz_Lq, mul_zero, add_zero, pderiv_X]
  rw [Finset.sum_eq_single w]
  · simp [Q04.Qb, map_rat_smul]
  · intro v _ hv
    simp [Pi.single_apply, hv.symm]
  · intro h; exact absurd (Finset.mem_Icc.mpr ⟨hw, hwK⟩) h

/-- the main formula: `[𝔷^s] exp(Δ) = Q_s · exp(Δ^H)` -/
lemma cz_expT (K p q : ℕ) (s : ℕ →₀ ℕ) :
    ∀ M, deg s ≤ M → (∀ w ∈ s.support, 2 ≤ w ∧ w ≤ K) →
      cz s (Q04.expT M (Q04.DeltaT K p q)) =
        Qz s * Q04.expT (M - deg s) (Q04.Hb K p q) := by
  induction s using Finsupp.induction with
  | zero =>
    intro M _ _
    simp [cz_zero, Qz, deg, Q04.expT_map, Q04.K0_Delta]
  | single_add w e f hw he ih =>
    intro M hM hs
    have hfw : f w = 0 := by simpa using hw
    have hsw := hs w (by rw [Finsupp.support_single_add hw he]; simp)
    rw [deg_single_add] at hM
    have hΔ : Dz w (Q04.DeltaT K p q) = Lq (Q04.Qb w) := Dz_Delta K p q w hsw.1 hsw.2
    rw [add_comm, cz_split w e f hfw, Dz_iter_expT w _ _ hΔ M e (by omega), cz_mulL,
      ih (M - e) (by omega) (fun v hv => hs v (by
        rw [Finsupp.support_single_add hw he]; exact Finset.mem_cons_of_mem hv))]
    rw [add_comm, deg_single_add]
    unfold Qz
    rw [Finsupp.support_single_add hw he, Finset.prod_cons]
    have e1 : ∀ v ∈ f.support, (((Finsupp.single w e + f) v).factorial : ℚ)⁻¹ •
        Pw v ^ ((Finsupp.single w e + f) v) = ((f v).factorial : ℚ)⁻¹ • Pw v ^ f v := by
      intro v hv
      have : v ≠ w := fun h => hw (h ▸ hv)
      simp [Finsupp.single_apply, this.symm]
    rw [Finset.prod_congr rfl e1]
    simp only [Finsupp.add_apply, Finsupp.single_eq_same, hfw, add_zero, Pw_eq]
    rw [show M - (e + deg f) = M - e - deg f by omega]
    rw [Algebra.smul_def, Algebra.smul_def]
    ring

end R01

namespace R01
open MvPolynomial

lemma deg_le_wt (s : ℕ →₀ ℕ) (hs : IsZetaMono s) : deg s ≤ wt s := by
  unfold deg; rw [wt_eq, Finsupp.sum]
  refine Finset.sum_le_sum fun w hw => ?_
  have := hs w hw; nlinarith

lemma supp_le_wt (s : ℕ →₀ ℕ) (w : ℕ) (hw : w ∈ s.support) : w ≤ wt s := by
  have h1 := Q04.mono_bound s w
  have h2 : 1 ≤ s w := Nat.one_le_iff_ne_zero.mpr (Finsupp.mem_support_iff.mp hw)
  unfold wt; nlinarith

lemma Dsharp_coeff (i m p q : ℕ) (s : ℕ →₀ ℕ) (hs : IsZetaMono s) (hw : wt s ≤ i + m) :
    (Q04.DsharpF i m p q).coeff s = ((i.factorial * m.factorial : ℕ) : ℚ) *
      (Qz s * Q04.expT (i + m - deg s) (Q04.Hb (i + m) p q)).coeff
        (Finsupp.single 0 i + Finsupp.single 1 m) := by
  unfold Q04.DsharpF
  rw [coeff_smul, smul_eq_mul, ← cz_coeff, cz_expT (i + m) p q s (i + m)
    ((deg_le_wt s hs).trans hw) (fun w hws => ⟨hs w hws, (supp_le_wt s w hws).trans hw⟩)]

lemma wt_zero : wt 0 = 0 := by simp [wt]
lemma deg_zero : deg 0 = 0 := by simp [deg]
lemma Qz_zero : Qz 0 = 1 := by simp [Qz]

lemma Dsharp_coeff0 (i m p q : ℕ) :
    (Q04.DsharpF i m p q).coeff 0 = ((i.factorial * m.factorial : ℕ) : ℚ) * Q04.gH p q i m := by
  rw [Dsharp_coeff i m p q 0 (by intro w hw; simp at hw) (by rw [wt_zero]; omega), Qz_zero,
    one_mul, deg_zero, Nat.sub_zero]
  rfl

lemma Qz_hom (s : ℕ →₀ ℕ) : IsWeightedHomogeneous (fun _ : Fin 2 => 1) (Qz s) (wt s) := by
  unfold Qz
  rw [wt_eq]
  refine IsWeightedHomogeneous.prod _ _ _ fun w _ => ?_
  rw [smul_eq_C_mul]
  have hX : ∀ i : Fin 2, IsWeightedHomogeneous (fun _ : Fin 2 => 1) (X i : Q2) 1 :=
    fun i => isWeightedHomogeneous_X ℚ _ i
  have h1 : IsWeightedHomogeneous (fun _ : Fin 2 => 1) (Pw w) w := by
    unfold Pw
    rw [smul_eq_C_mul]
    refine IsWeightedHomogeneous.C_mul ?_ _
    have a := ((hX 0).add (hX 1)).pow w
    have b := (hX 0).pow w
    have c := (hX 1).pow w
    simp only [smul_eq_mul, mul_one] at a b c
    exact (a.sub b).sub c
  have := (h1.pow (s w)).C_mul (((s w).factorial : ℚ)⁻¹)
  simpa [mul_comm] using this

lemma one_weight (u : Fin 2 →₀ ℕ) : Finsupp.weight (fun _ : Fin 2 => 1) u = u 0 + u 1 := by
  rw [Finsupp.weight_apply, Finsupp.sum_fintype _ _ (by simp)]
  simp [Fin.sum_univ_two]

lemma sub_mon (i m : ℕ) (u : Fin 2 →₀ ℕ) :
    Finsupp.single 0 i + Finsupp.single 1 m - u = Finsupp.single 0 (i - u 0) + Finsupp.single 1 (m - u 1) := by
  ext x; fin_cases x <;> simp

lemma le_mon (i m : ℕ) (u : Fin 2 →₀ ℕ) :
    u ≤ Finsupp.single 0 i + Finsupp.single 1 m ↔ u 0 ≤ i ∧ u 1 ≤ m := by
  rw [Finsupp.le_def]
  constructor
  · intro h; exact ⟨by simpa using h 0, by simpa using h 1⟩
  · intro h j; fin_cases j <;> simp [h.1, h.2]

/-- the core identity `(★)` -/
lemma core (i m p q : ℕ) (s : ℕ →₀ ℕ) (hs : IsZetaMono s) :
    (Q04.DsharpF i m p q).coeff s = ∑ u ∈ (Qz s).support, (Qz s).coeff u *
      (((i.descFactorial (u 0) * m.descFactorial (u 1) : ℕ) : ℚ) *
        (Q04.DsharpF (i - u 0) (m - u 1) p q).coeff 0) := by
  have hhom := Qz_hom s
  by_cases hw : wt s ≤ i + m
  · rw [Dsharp_coeff i m p q s hs hw]
    conv_lhs => rw [(Qz s).as_sum, Finset.sum_mul, coeff_sum, Finset.mul_sum]
    refine Finset.sum_congr rfl fun u hu => ?_
    have hu0 : u 0 + u 1 = wt s := by
      have := hhom (mem_support_iff.mp hu); rwa [one_weight] at this
    rw [coeff_monomial_mul', Dsharp_coeff0]
    by_cases hle : u 0 ≤ i ∧ u 1 ≤ m
    · rw [if_pos ((le_mon i m u).mpr hle), sub_mon,
        Q04.Hb_trunc _ _ _ _ _ _ (by have := deg_le_wt s hs; omega) (by omega)]
      have e1 := Nat.factorial_mul_descFactorial hle.1
      have e2 := Nat.factorial_mul_descFactorial hle.2
      have e1' : ((i - u 0).factorial : ℚ) * (i.descFactorial (u 0)) = i.factorial := by
        exact_mod_cast e1
      have e2' : ((m - u 1).factorial : ℚ) * (m.descFactorial (u 1)) = m.factorial := by
        exact_mod_cast e2
      push_cast
      rw [← e1', ← e2']; ring
    · rw [if_neg (fun h => hle ((le_mon i m u).mp h))]
      have : i.descFactorial (u 0) * m.descFactorial (u 1) = 0 := by
        rcases not_and_or.mp hle with h | h
        · rw [show i.descFactorial (u 0) = 0 from Nat.descFactorial_eq_zero_iff_lt.mpr (by omega),
            zero_mul]
        · rw [show m.descFactorial (u 1) = 0 from Nat.descFactorial_eq_zero_iff_lt.mpr (by omega),
            mul_zero]
      rw [this]; simp
  · have h0 : (Q04.DsharpF i m p q).coeff s = 0 := by
      by_contra hne
      have := (Q04.DsharpF_weight i m p q s (mem_support_iff.mpr hne)).2
      exact hw this
    rw [h0]
    refine (Finset.sum_eq_zero fun u hu => ?_).symm
    have hu0 : u 0 + u 1 = wt s := by
      have := hhom (mem_support_iff.mp hu); rwa [one_weight] at this
    have : i.descFactorial (u 0) * m.descFactorial (u 1) = 0 := by
      by_cases h : u 0 ≤ i
      · rw [show m.descFactorial (u 1) = 0 from Nat.descFactorial_eq_zero_iff_lt.mpr (by omega),
          mul_zero]
      · rw [show i.descFactorial (u 0) = 0 from Nat.descFactorial_eq_zero_iff_lt.mpr (by omega),
          zero_mul]
    rw [this]; simp

end R01

namespace R01
open MvPolynomial

lemma PhiN_sub (n : ℕ) (F : MvPolynomial (Fin 3) ℚ) (S : Finset (Fin 3 →₀ ℕ)) (h : F.support ⊆ S) :
    PhiN n F = ∑ t ∈ S, (F.coeff t * rhoR (t 0) n) • Q04.DsharpF (t 1) (t 2) (n + t 0) n := by
  unfold PhiN
  refine Finset.sum_subset h fun t _ ht => ?_
  rw [notMem_support_iff.mp ht]; simp

lemma PhiN_add (n : ℕ) (F G : MvPolynomial (Fin 3) ℚ) : PhiN n (F + G) = PhiN n F + PhiN n G := by
  classical
  rw [PhiN_sub n (F + G) (F.support ∪ G.support) support_add,
    PhiN_sub n F (F.support ∪ G.support) Finset.subset_union_left,
    PhiN_sub n G (F.support ∪ G.support) Finset.subset_union_right, ← Finset.sum_add_distrib]
  refine Finset.sum_congr rfl fun t _ => ?_
  rw [MvPolynomial.coeff_add, add_mul, add_smul]

lemma PhiN_zero (n : ℕ) : PhiN n 0 = 0 := by simp [PhiN]

lemma PhiN_smul (n : ℕ) (c : ℚ) (F : MvPolynomial (Fin 3) ℚ) : PhiN n (c • F) = c • PhiN n F := by
  rw [PhiN_sub n (c • F) F.support support_smul, PhiN, Finset.smul_sum]
  refine Finset.sum_congr rfl fun t _ => ?_
  rw [coeff_smul, smul_smul, smul_eq_mul, mul_assoc]

lemma PhiN_sum (n : ℕ) {ι : Type*} (S : Finset ι) (f : ι → MvPolynomial (Fin 3) ℚ) :
    PhiN n (∑ i ∈ S, f i) = ∑ i ∈ S, PhiN n (f i) := by
  classical
  induction S using Finset.induction_on with
  | empty => simp [PhiN_zero]
  | insert a S ha ih => rw [Finset.sum_insert ha, Finset.sum_insert ha, PhiN_add, ih]

lemma PhiN_monomial (n : ℕ) (t : Fin 3 →₀ ℕ) (a : ℚ) :
    PhiN n (monomial t a) = (a * rhoR (t 0) n) • Q04.DsharpF (t 1) (t 2) (n + t 0) n := by
  classical
  rw [PhiN_sub n _ {t} (by rw [support_monomial]; split_ifs <;> simp), Finset.sum_singleton,
    coeff_monomial, if_pos rfl]

lemma iter_add (i : Fin 3) (r : ℕ) (F G : MvPolynomial (Fin 3) ℚ) :
    (⇑(pderiv i))^[r] (F + G) = (⇑(pderiv i))^[r] F + (⇑(pderiv i))^[r] G := by
  induction r with
  | zero => rfl
  | succ r ih => rw [Function.iterate_succ_apply', Function.iterate_succ_apply',
      Function.iterate_succ_apply', ih, map_add]

lemma iter_mono (i : Fin 3) (r : ℕ) (t : Fin 3 →₀ ℕ) (c : ℚ) :
    (⇑(pderiv i))^[r] (monomial t c) = monomial (t - Finsupp.single i r) (c * ((t i).descFactorial r : ℕ)) := by
  induction r with
  | zero => simp
  | succ r ih =>
    rw [Function.iterate_succ_apply', ih, pderiv_monomial]
    have e1 : t - Finsupp.single i r - Finsupp.single i 1 = t - Finsupp.single i (r + 1) := by
      ext x
      simp only [Finsupp.tsub_apply, Finsupp.single_apply]
      split_ifs <;> omega
    have e2 : c * (((t i).descFactorial r : ℕ) : ℚ) * (((t - Finsupp.single i r) i : ℕ) : ℚ) =
        c * (((t i).descFactorial (r + 1) : ℕ) : ℚ) := by
      simp only [Finsupp.tsub_apply, Finsupp.single_eq_same]
      rw [Nat.descFactorial_succ]; push_cast; ring
    rw [e1, e2]

lemma QOp_add (Q : Q2) (F G : MvPolynomial (Fin 3) ℚ) : QOp Q (F + G) = QOp Q F + QOp Q G := by
  unfold QOp
  rw [← Finset.sum_add_distrib]
  refine Finset.sum_congr rfl fun u _ => ?_
  rw [iter_add, iter_add, smul_add]

lemma QOp_monomial (Q : Q2) (t : Fin 3 →₀ ℕ) (a : ℚ) : QOp Q (monomial t a) =
    ∑ u ∈ Q.support, Q.coeff u • monomial (t - Finsupp.single 2 (u 1) - Finsupp.single 1 (u 0))
      (a * ((t 2).descFactorial (u 1) : ℕ) * ((t 1).descFactorial (u 0) : ℕ)) := by
  unfold QOp
  refine Finset.sum_congr rfl fun u _ => ?_
  rw [iter_mono, iter_mono]
  congr 3
  simp

end R01

theorem R01.proof_3 : R01.stmt_3 := by
  intro n _ F s hs
  induction F using MvPolynomial.induction_on' with
  | monomial t a =>
    rw [R01.QOp_monomial, R01.PhiN_sum, R01.PhiN_monomial, MvPolynomial.coeff_smul, smul_eq_mul,
      R01.core _ _ _ _ s hs, MvPolynomial.coeff_sum, Finset.mul_sum]
    refine Finset.sum_congr rfl fun u _ => ?_
    rw [R01.PhiN_smul, R01.PhiN_monomial, MvPolynomial.coeff_smul, MvPolynomial.coeff_smul,
      smul_eq_mul, smul_eq_mul]
    simp only [Finsupp.coe_tsub, Pi.sub_apply, Finsupp.single_apply]
    norm_num
    push_cast
    ring
  | add F G hF hG =>
    rw [R01.PhiN_add, R01.QOp_add, R01.PhiN_add, MvPolynomial.coeff_add, MvPolynomial.coeff_add,
      hF, hG]

namespace R01
open MvPolynomial

/-- `∂_a^{u_0} ∂_b^{u_1}` -/
noncomputable def opU (u : Fin 2 →₀ ℕ) (F : MvPolynomial (Fin 3) ℚ) : MvPolynomial (Fin 3) ℚ :=
  (⇑(pderiv (1 : Fin 3)))^[u 0] ((⇑(pderiv (2 : Fin 3)))^[u 1] F)

lemma QOp_sub (Q : Q2) (F : MvPolynomial (Fin 3) ℚ) (S : Finset (Fin 2 →₀ ℕ)) (h : Q.support ⊆ S) :
    QOp Q F = ∑ u ∈ S, Q.coeff u • opU u F := by
  unfold QOp
  refine Finset.sum_subset h fun u _ hu => ?_
  rw [notMem_support_iff.mp hu, zero_smul]

lemma QOp_addQ (Q R : Q2) (F : MvPolynomial (Fin 3) ℚ) : QOp (Q + R) F = QOp Q F + QOp R F := by
  classical
  rw [QOp_sub (Q + R) F (Q.support ∪ R.support) support_add,
    QOp_sub Q F (Q.support ∪ R.support) Finset.subset_union_left,
    QOp_sub R F (Q.support ∪ R.support) Finset.subset_union_right, ← Finset.sum_add_distrib]
  refine Finset.sum_congr rfl fun t _ => ?_
  rw [MvPolynomial.coeff_add, add_smul]

lemma QOp_monoQ (v : Fin 2 →₀ ℕ) (c : ℚ) (F : MvPolynomial (Fin 3) ℚ) :
    QOp (monomial v c) F = c • opU v F := by
  classical
  rw [QOp_sub _ F {v} (by rw [support_monomial]; split_ifs <;> simp), Finset.sum_singleton,
    coeff_monomial, if_pos rfl]

lemma iter_smul (i : Fin 3) (r : ℕ) (c : ℚ) (F : MvPolynomial (Fin 3) ℚ) :
    (⇑(pderiv i))^[r] (c • F) = c • (⇑(pderiv i))^[r] F := by
  induction r with
  | zero => rfl
  | succ r ih => rw [Function.iterate_succ_apply', Function.iterate_succ_apply', ih, (pderiv i).map_smul]

lemma iter_zero (i : Fin 3) (r : ℕ) : (⇑(pderiv i))^[r] (0 : MvPolynomial (Fin 3) ℚ) = 0 := by
  induction r with
  | zero => rfl
  | succ r ih => rw [Function.iterate_succ_apply', ih, map_zero]

lemma QOp_smulF (Q : Q2) (c : ℚ) (F : MvPolynomial (Fin 3) ℚ) : QOp Q (c • F) = c • QOp Q F := by
  unfold QOp; rw [Finset.smul_sum]
  refine Finset.sum_congr rfl fun u _ => ?_
  rw [iter_smul, iter_smul, smul_comm]

lemma QOp_zeroF (Q : Q2) : QOp Q 0 = 0 := by
  unfold QOp; simp [iter_zero]

lemma QOp_mulX (Q : Q2) (i : Fin 2) (F : MvPolynomial (Fin 3) ℚ) :
    QOp (Q * X i) F = QOp Q (pderiv (i.succ) F) := by
  induction Q using MvPolynomial.induction_on' with
  | monomial v c =>
    rw [X, monomial_mul, mul_one, QOp_monoQ, QOp_monoQ]
    congr 1
    unfold opU
    rcases (show i = 0 ∨ i = 1 by fin_cases i <;> simp) with rfl | rfl
    · rw [show (v + Finsupp.single (0 : Fin 2) 1 : Fin 2 →₀ ℕ) 0 = v 0 + 1 by simp,
        show (v + Finsupp.single (0 : Fin 2) 1 : Fin 2 →₀ ℕ) 1 = v 1 by simp]
      rw [Function.iterate_succ_apply]
      congr 1
      exact ((commute_pderiv (R := ℚ) (1 : Fin 3) 2).iterate_right (v 1)) F
    · rw [show (v + Finsupp.single (1 : Fin 2) 1 : Fin 2 →₀ ℕ) 0 = v 0 by simp,
        show (v + Finsupp.single (1 : Fin 2) 1 : Fin 2 →₀ ℕ) 1 = v 1 + 1 by simp]
      rw [Function.iterate_succ_apply]
      rfl
  | add P R hP hR => rw [add_mul, QOp_addQ, QOp_addQ, hP, hR]

lemma QOp_C (a : ℚ) (F : MvPolynomial (Fin 3) ℚ) : QOp (C a) F = a • F := by
  rw [C_apply, QOp_monoQ]; rfl

lemma QOp_smulQ (a : ℚ) (Q : Q2) (F : MvPolynomial (Fin 3) ℚ) : QOp (a • Q) F = a • QOp Q F := by
  rw [QOp_sub _ _ Q.support support_smul, QOp_sub Q F Q.support le_rfl, Finset.smul_sum]
  refine Finset.sum_congr rfl fun u _ => ?_
  rw [coeff_smul, smul_eq_mul, mul_smul]

lemma QOp_mul (Q R : Q2) (F : MvPolynomial (Fin 3) ℚ) : QOp (Q * R) F = QOp Q (QOp R F) := by
  induction R using MvPolynomial.induction_on generalizing F with
  | C a => rw [QOp_C, QOp_smulF, mul_comm, ← smul_eq_C_mul, QOp_smulQ]
  | add P R hP hR => rw [mul_add, QOp_addQ, QOp_addQ, hP, hR, QOp_add]
  | mul_X P i hP => rw [← mul_assoc, QOp_mulX, QOp_mulX, hP]

end R01

namespace R01
open MvPolynomial

lemma QOp_X (i : Fin 2) (F : MvPolynomial (Fin 3) ℚ) : QOp (X i) F = pderiv i.succ F := by
  rw [← one_mul (X i), QOp_mulX, ← C_1, QOp_C, one_smul]

lemma QOp_sigma1 : QOp (X 0 + X 1) = ⇑Q04.Dd := by
  funext F
  rw [QOp_addQ, QOp_X, QOp_X]
  simp [Q04.Dd]

lemma QOp_pow (Q : Q2) (r : ℕ) (F : MvPolynomial (Fin 3) ℚ) : QOp (Q ^ r) F = (QOp Q)^[r] F := by
  induction r generalizing F with
  | zero => rw [pow_zero, ← C_1, QOp_C, one_smul]; rfl
  | succ r ih => rw [pow_succ, QOp_mul, ih, Function.iterate_succ_apply]

lemma QOp_zeroQ (F : MvPolynomial (Fin 3) ℚ) : QOp 0 F = 0 := by
  simp [QOp]

lemma QOp_sumQ {ι : Type*} (S : Finset ι) (f : ι → Q2) (F : MvPolynomial (Fin 3) ℚ) :
    QOp (∑ i ∈ S, f i) F = ∑ i ∈ S, QOp (f i) F := by
  classical
  induction S using Finset.induction_on with
  | empty => simp [QOp_zeroQ]
  | insert a S ha ih => rw [Finset.sum_insert ha, Finset.sum_insert ha, QOp_addQ, ih]

lemma Dd_iter (r k n : ℕ) :
    (⇑Q04.Dd)^[r] (Q04.Ecoef k n) = if r ≤ k then Q04.Ecoef (k - r) n else 0 := by
  induction r with
  | zero => simp
  | succ r ih =>
    rw [Function.iterate_succ_apply', ih]
    by_cases h : r + 1 ≤ k
    · rw [if_pos (by omega), if_pos h, show k - r = (k - (r + 1)) + 1 by omega, Q04.Dd_Ecoef]
    · rw [if_neg h]
      by_cases h' : r ≤ k
      · rw [if_pos h', show k - r = 0 by omega, Ecoef_zero, Derivation.map_one_eq_zero]
      · rw [if_neg h', map_zero]

lemma Pw2 : Pw 2 = -(X 0 * X 1) := by
  calc Pw 2 = C ((-1 : ℚ) ^ (2 + 1) / ((2 : ℕ) : ℚ)) * ((X 0 + X 1) ^ 2 - X 0 ^ 2 - X 1 ^ 2) := by
          rw [Pw, smul_eq_C_mul]
    _ = C ((-1 : ℚ) ^ (2 + 1) / ((2 : ℕ) : ℚ)) * C 2 * (X 0 * X 1) := by
          rw [map_ofNat C 2]; ring
    _ = -(X 0 * X 1) := by rw [← C_mul]; norm_num

lemma Qz_single2 (e : ℕ) (he : e ≠ 0) :
    Qz (Finsupp.single 2 e) = C ((-1 : ℚ) ^ e * ((e.factorial : ℚ)⁻¹)) * (X 0 * X 1) ^ e := by
  unfold Qz
  rw [Finsupp.support_single_ne_zero _ he, Finset.prod_singleton, Finsupp.single_eq_same, Pw2,
    neg_pow, smul_eq_C_mul, C_mul, map_pow, map_neg, map_one]
  ring

lemma zeta2_mono (e : ℕ) : IsZetaMono (Finsupp.single 2 e) := by
  intro w hw
  have := Finsupp.support_single_subset hw
  simp at this; omega

lemma phi_sigma2 (n : ℕ) (hn : 1 ≤ n) (e : ℕ) (he : e ≠ 0) (G : MvPolynomial (Fin 3) ℚ) :
    (PhiN n (QOp ((X 0 * X 1) ^ e) G)).coeff 0 =
      (-1 : ℚ) ^ e * (e.factorial : ℚ) * (PhiN n G).coeff (Finsupp.single 2 e) := by
  rw [R01.proof_3 n hn G _ (zeta2_mono e), Qz_single2 e he, QOp_mul, QOp_C, PhiN_smul, coeff_smul,
    smul_eq_mul, ← mul_assoc, ← mul_assoc]
  have h1 : (e.factorial : ℚ) ≠ 0 := by positivity
  have h2 : ((-1 : ℚ) ^ e) ^ 2 = 1 := by rw [← pow_mul, mul_comm, pow_mul]; norm_num
  field_simp
  rw [h2, mul_one]

end R01

theorem R01.proof_5 : R01.stmt_5 := by
  intro s hs h0 k n hn
  rw [R01.proof_2 k n hn, MvPolynomial.coeff_smul, smul_eq_mul, R01.proof_3 n hn _ s hs,
    R01.proof_4 s hs h0, R01.QOp_sumQ, R01.PhiN_sum, MvPolynomial.coeff_sum, Finset.mul_sum]
  refine Finset.sum_congr rfl fun e he => ?_
  have he' := Finset.mem_Icc.mp he
  have h2e : 2 * e ≤ R01.wt s := by omega
  rw [mul_assoc, R01.QOp_mul, R01.QOp_C, R01.PhiN_smul, MvPolynomial.coeff_smul, smul_eq_mul,
    mul_comm ((MvPolynomial.X 0 + MvPolynomial.X 1) ^ (R01.wt s - 2 * e)), R01.QOp_mul,
    R01.QOp_pow (MvPolynomial.X 0 + MvPolynomial.X 1), R01.QOp_sigma1, R01.Dd_iter]
  by_cases hk : R01.wt s ≤ k + 2 * e
  · rw [if_pos (by omega), if_pos hk, R01.phi_sigma2 n hn e (by omega), R01.proof_2 _ n hn,
      MvPolynomial.coeff_smul, smul_eq_mul, show k - (R01.wt s - 2 * e) = k + 2 * e - R01.wt s by omega]
    ring
  · rw [if_neg (by omega), if_neg hk, R01.QOp_zeroF, R01.PhiN_zero]
    simp

namespace R01

/-- `[ζ(2)^e]M_{k'}` for an integer `k' = k - W + 2e`, written as a series of the `n`-th terms -/
lemma zM_shift (e k W : ℕ) (α : ℝ) :
    zM (Finsupp.single 2 e) ((k : ℤ) - W + 2 * e) α =
      ∑' n, (if W ≤ k + 2 * e then cser (Finsupp.single 2 e) (k + 2 * e - W) α n else 0) := by
  by_cases h : W ≤ k + 2 * e
  · simp only [if_pos h]
    rw [show (k : ℤ) - W + 2 * e = ((k + 2 * e - W : ℕ) : ℤ) by push_cast [h]; ring, zM_nat]
  · simp only [if_neg h, tsum_zero]
    unfold zM; rw [if_neg (by omega)]

lemma ladder_series {α : ℝ} (hα0 : 0 < α) (hα1 : α < Real.pi) (s : ℕ →₀ ℕ) (hs : IsZetaMono s)
    (h0 : s ≠ 0) (k : ℕ) :
    zM s k α = ∑ e ∈ Finset.Icc 1 (wt s / 2),
        ((-1 : ℝ) ^ e * (e.factorial : ℝ) * (gam e s : ℝ)) *
          zM (Finsupp.single 2 e) ((k : ℤ) - wt s + 2 * e) α := by
  simp only [zM_shift]
  rw [zM_nat]
  have hterm : ∀ n, cser s k α n = ∑ e ∈ Finset.Icc 1 (wt s / 2),
      ((-1 : ℝ) ^ e * (e.factorial : ℝ) * (gam e s : ℝ)) *
        (if wt s ≤ k + 2 * e then cser (Finsupp.single 2 e) (k + 2 * e - wt s) α n else 0) := by
    intro n
    unfold cser
    rw [R01.proof_5 s hs h0 k (n + 1) (by omega)]
    push_cast
    rw [Finset.sum_mul]
    refine Finset.sum_congr rfl fun e _ => ?_
    split_ifs <;> push_cast <;> ring
  rw [tsum_congr hterm, Summable.tsum_finsetSum]
  · refine Finset.sum_congr rfl fun e _ => ?_
    rw [tsum_mul_left]
  · intro e _
    refine Summable.mul_left _ ?_
    by_cases h : wt s ≤ k + 2 * e
    · simp only [if_pos h]; exact summable_cser hα0 hα1 _ _
    · simp only [if_neg h]; exact summable_zero

end R01

theorem R01.proof_6 : R01.stmt_6 := fun α hα0 hα1 =>
  ⟨R01.allSummable hα0 hα1, fun s hs h0 k => R01.ladder_series hα0 hα1 s hs h0 k⟩

namespace R01
open MvPolynomial

lemma pN0 : pN 0 = 2 := rfl
lemma pN1 : pN 1 = X 0 := rfl
lemma pN_ss (w : ℕ) : pN (w + 2) = X 0 * pN (w + 1) - X 1 * pN w := rfl
lemma pN2 : pN 2 = X 0 ^ 2 - 2 * X 1 := by rw [pN_ss, pN1, pN0]; ring
lemma pN3 : pN 3 = X 0 ^ 3 - 3 * X 0 * X 1 := by rw [pN_ss, pN2, pN1]; ring
lemma pN4 : pN 4 = X 0 ^ 4 - 4 * X 0 ^ 2 * X 1 + 2 * X 1 ^ 2 := by rw [pN_ss, pN3, pN2]; ring
lemma pN5 : pN 5 = X 0 ^ 5 - 5 * X 0 ^ 3 * X 1 + 5 * X 0 * X 1 ^ 2 := by
  rw [pN_ss, pN4, pN3]; ring
lemma pN6 : pN 6 = X 0 ^ 6 - 6 * X 0 ^ 4 * X 1 + 9 * X 0 ^ 2 * X 1 ^ 2 - 2 * X 1 ^ 3 := by
  rw [pN_ss, pN5, pN4]; ring

lemma Qsym_single (w : ℕ) :
    QsymN (Finsupp.single w 1) = C ((-1 : ℚ) ^ (w + 1) / (w : ℚ)) * (X 0 ^ w - pN w) := by
  unfold QsymN
  rw [Finsupp.support_single_ne_zero _ one_ne_zero, Finset.prod_singleton, Finsupp.single_eq_same]
  simp [smul_eq_C_mul]

lemma mon2_inj (a b c d : ℕ) :
    (Finsupp.single 0 a + Finsupp.single 1 b : Fin 2 →₀ ℕ) = Finsupp.single 0 c + Finsupp.single 1 d ↔
      a = c ∧ b = d := by
  constructor
  · intro h
    exact ⟨by simpa using DFunLike.congr_fun h 0, by simpa using DFunLike.congr_fun h 1⟩
  · rintro ⟨rfl, rfl⟩; rfl

lemma Cnum (c : ℚ) (m : ℕ) [m.AtLeastTwo] : C c * (OfNat.ofNat m : Q2) = C (c * OfNat.ofNat m) := by
  rw [C_mul, map_ofNat]

lemma mC (a b : ℕ) (c : ℚ) : C c * X 0 ^ a * X 1 ^ b =
    (monomial (Finsupp.single 0 a + Finsupp.single 1 b) c : Q2) := by
  rw [mon_C]; ring

lemma QsymN5 : QsymN (Finsupp.single 5 1) = monomial (Finsupp.single 0 3 + Finsupp.single 1 1) 1 +
    monomial (Finsupp.single 0 1 + Finsupp.single 1 2) (-1) := by
  rw [Qsym_single, pN5, ← mC, ← mC]
  have h : C ((-1 : ℚ) ^ (5 + 1) / ((5 : ℕ) : ℚ)) * (5 : Q2) = 1 := by
    rw [Cnum]; norm_num
  rw [map_neg, map_one]
  linear_combination (X 0 ^ 3 * X 1 - X 0 * X 1 ^ 2) * h

lemma QsymN4 : QsymN (Finsupp.single 4 1) = monomial (Finsupp.single 0 2 + Finsupp.single 1 1) (-1) +
    monomial (Finsupp.single 0 0 + Finsupp.single 1 2) (1 / 2) := by
  rw [Qsym_single, pN4, ← mC, ← mC]
  have h : C ((-1 : ℚ) ^ (4 + 1) / ((4 : ℕ) : ℚ)) * (4 : Q2) = -1 := by
    rw [Cnum]; norm_num
  have h2 : C ((-1 : ℚ) ^ (4 + 1) / ((4 : ℕ) : ℚ)) * (2 : Q2) = -C (1 / 2) := by
    rw [Cnum, ← map_neg]; norm_num
  rw [map_neg, map_one]
  linear_combination (X 0 ^ 2 * X 1) * h - X 1 ^ 2 * h2

lemma QsymN6 : QsymN (Finsupp.single 6 1) = monomial (Finsupp.single 0 4 + Finsupp.single 1 1) (-1) +
    monomial (Finsupp.single 0 2 + Finsupp.single 1 2) (3 / 2) +
    monomial (Finsupp.single 0 0 + Finsupp.single 1 3) (-1 / 3) := by
  rw [Qsym_single, pN6, ← mC, ← mC, ← mC]
  have h : C ((-1 : ℚ) ^ (6 + 1) / ((6 : ℕ) : ℚ)) * (6 : Q2) = -1 := by
    rw [Cnum]; norm_num
  have h2 : C ((-1 : ℚ) ^ (6 + 1) / ((6 : ℕ) : ℚ)) * (9 : Q2) = -C (3 / 2) := by
    rw [Cnum, ← map_neg]; norm_num
  have h3 : C ((-1 : ℚ) ^ (6 + 1) / ((6 : ℕ) : ℚ)) * (2 : Q2) = C (-1 / 3) := by
    rw [Cnum]; norm_num
  rw [map_neg, map_one]
  linear_combination (X 0 ^ 4 * X 1) * h - X 0 ^ 2 * X 1 ^ 2 * h2 + X 1 ^ 3 * h3

end R01

namespace R01
open MvPolynomial

lemma QsymN23 : QsymN (Finsupp.single 2 1 + Finsupp.single 3 1) =
    monomial (Finsupp.single 0 1 + Finsupp.single 1 2) (-1) := by
  unfold QsymN
  have h2 : (2 : ℕ) ∉ (Finsupp.single 3 1 : ℕ →₀ ℕ).support := by simp
  rw [Finsupp.support_single_add h2 one_ne_zero, Finset.prod_cons,
    Finsupp.support_single_ne_zero _ one_ne_zero, Finset.prod_singleton]
  simp only [Finsupp.add_apply, Finsupp.single_eq_same, Finsupp.single_apply]
  norm_num
  rw [pN2, pN3, ← mC]
  simp only [smul_eq_C_mul]
  simp only [one_div, map_one]
  have h : C (2⁻¹ : ℚ) * (2 : Q2) = 1 := by rw [Cnum]; norm_num
  have h' : C (3⁻¹ : ℚ) * (3 : Q2) = 1 := by rw [Cnum]; norm_num
  linear_combination (X 0 * X 1 ^ 2) * (C (3⁻¹ : ℚ) * (3 : Q2)) * h + (X 0 * X 1 ^ 2) * h'

lemma QsymN33 : QsymN (Finsupp.single 3 2) =
    monomial (Finsupp.single 0 2 + Finsupp.single 1 2) (1 / 2) := by
  unfold QsymN
  rw [Finsupp.support_single_ne_zero _ two_ne_zero, Finset.prod_singleton, Finsupp.single_eq_same]
  norm_num
  rw [pN3, ← mC]
  simp only [smul_eq_C_mul]
  simp only [one_div]
  have h' : C (3⁻¹ : ℚ) * (3 : Q2) = 1 := by rw [Cnum]; norm_num
  linear_combination (C (2⁻¹ : ℚ) * X 0 ^ 2 * X 1 ^ 2) * (C (3⁻¹ : ℚ) * (3 : Q2) + 1) * h'

end R01

namespace R01
open MvPolynomial

lemma wt_single (w e : ℕ) : wt (Finsupp.single w e) = e * w := Q04.wsum_single w e

lemma wt23 : wt (Finsupp.single 2 1 + Finsupp.single 3 1) = 5 := by
  unfold wt
  rw [Finsupp.sum_add_index' (by simp) (by intros; ring), Finsupp.sum_single_index (by simp),
    Finsupp.sum_single_index (by simp)]

lemma single_ne (w e : ℕ) (he : e ≠ 0) : (Finsupp.single w e : ℕ →₀ ℕ) ≠ 0 :=
  Finsupp.single_ne_zero.mpr he

lemma mono_single (w e : ℕ) (hw : 2 ≤ w) : IsZetaMono (Finsupp.single w e) := by
  intro v hv
  have := Finsupp.support_single_subset hv
  simp at this; omega

lemma mono23 : IsZetaMono (Finsupp.single 2 1 + Finsupp.single 3 1) := by
  intro v hv
  have := Finsupp.support_add hv
  rcases Finset.mem_union.mp this with h | h
  · have := Finsupp.support_single_subset h; simp at this; omega
  · have := Finsupp.support_single_subset h; simp at this; omega

lemma ne23 : (Finsupp.single 2 1 + Finsupp.single 3 1 : ℕ →₀ ℕ) ≠ 0 := by
  intro h; have := DFunLike.congr_fun h 2; simp at this

lemma mon2_inj' (b c d : ℕ) :
    (Finsupp.single 1 b : Fin 2 →₀ ℕ) = Finsupp.single 0 c + Finsupp.single 1 d ↔ c = 0 ∧ b = d := by
  rw [show (Finsupp.single 1 b : Fin 2 →₀ ℕ) = Finsupp.single 0 0 + Finsupp.single 1 b by simp,
    mon2_inj]; omega

lemma mon2_inj'' (a b d : ℕ) :
    (Finsupp.single 0 a + Finsupp.single 1 b : Fin 2 →₀ ℕ) = Finsupp.single 1 d ↔ a = 0 ∧ b = d := by
  rw [show (Finsupp.single 1 d : Fin 2 →₀ ℕ) = Finsupp.single 0 0 + Finsupp.single 1 d by simp,
    mon2_inj]

lemma gam5 : gam 1 (Finsupp.single 5 1) = 1 ∧ gam 2 (Finsupp.single 5 1) = -1 := by
  unfold gam; rw [wt_single, QsymN5]; simp [coeff_monomial, mon2_inj, mon2_inj', mon2_inj'']

lemma gam4 : gam 1 (Finsupp.single 4 1) = -1 ∧ gam 2 (Finsupp.single 4 1) = 1 / 2 := by
  unfold gam; rw [wt_single, QsymN4]; simp [coeff_monomial, mon2_inj, mon2_inj', mon2_inj'']

lemma gam6 : gam 1 (Finsupp.single 6 1) = -1 ∧ gam 2 (Finsupp.single 6 1) = 3 / 2 ∧
    gam 3 (Finsupp.single 6 1) = -1 / 3 := by
  unfold gam; rw [wt_single, QsymN6]; simp [coeff_monomial, mon2_inj, mon2_inj', mon2_inj'']

lemma gam23 : gam 1 (Finsupp.single 2 1 + Finsupp.single 3 1) = 0 ∧
    gam 2 (Finsupp.single 2 1 + Finsupp.single 3 1) = -1 := by
  unfold gam; rw [wt23, QsymN23]; simp [coeff_monomial, mon2_inj, mon2_inj', mon2_inj'']

lemma gam33 : gam 1 (Finsupp.single 3 2) = 0 ∧ gam 2 (Finsupp.single 3 2) = 1 / 2 ∧
    gam 3 (Finsupp.single 3 2) = 0 := by
  unfold gam; rw [wt_single, QsymN33]; simp [coeff_monomial, mon2_inj, mon2_inj', mon2_inj'']

end R01

namespace R01

lemma ite_Pf (W d e k n : ℕ) (h : W = d + 2 * e) :
    (if W ≤ k + 2 * e then (Q04.Pf (k + 2 * e - W) n).coeff (Finsupp.single 2 e) else 0) =
      (if d ≤ k then (Q04.Pf (k - d) n).coeff (Finsupp.single 2 e) else 0) := by
  subst h
  by_cases hd : d ≤ k
  · rw [if_pos (by omega), if_pos hd, show k + 2 * e - (d + 2 * e) = k - d by omega]
  · rw [if_neg (by omega), if_neg hd]

lemma Icc12 : Finset.Icc 1 2 = {1, 2} := by decide
lemma Icc13 : Finset.Icc 1 3 = {1, 2, 3} := by decide

end R01

theorem R01.proof_8 : R01.stmt_8 := by
  have hs := R01.mono_single 5 1 (by norm_num)
  have h0 := R01.single_ne 5 1 one_ne_zero
  refine ⟨fun k n hn => ?_, fun α hα0 hα1 => ⟨R01.allSummable hα0 hα1, fun k => ?_⟩⟩
  · rw [R01.proof_5 _ hs h0 k n hn, R01.wt_single, show 1 * 5 / 2 = 2 by norm_num, R01.Icc12,
      Finset.sum_pair (by norm_num), R01.gam5.1, R01.gam5.2,
      R01.ite_Pf (1 * 5) 3 1 k n (by norm_num), R01.ite_Pf (1 * 5) 1 2 k n (by norm_num)]
    norm_num
    split_ifs <;> ring
  · rw [R01.ladder_series hα0 hα1 _ hs h0 k, R01.wt_single, show 1 * 5 / 2 = 2 by norm_num, R01.Icc12,
      Finset.sum_pair (by norm_num), R01.gam5.1, R01.gam5.2]
    norm_num
    ring_nf

theorem R01.proof_9 : R01.stmt_9 := by
  intro α hα0 hα1
  refine ⟨R01.allSummable hα0 hα1, fun k => ⟨?_, ?_, ?_, ?_⟩⟩
  · rw [R01.ladder_series hα0 hα1 _ (R01.mono_single 4 1 (by norm_num)) (R01.single_ne 4 1 one_ne_zero) k,
      R01.wt_single, show 1 * 4 / 2 = 2 by norm_num, R01.Icc12,
      Finset.sum_pair (by norm_num), R01.gam4.1, R01.gam4.2]
    norm_num
    ring_nf
  · rw [R01.ladder_series hα0 hα1 _ (R01.mono_single 6 1 (by norm_num)) (R01.single_ne 6 1 one_ne_zero) k,
      R01.wt_single, show 1 * 6 / 2 = 3 by norm_num, R01.Icc13,
      Finset.sum_insert (by decide), Finset.sum_pair (by norm_num), R01.gam6.1, R01.gam6.2.1,
      R01.gam6.2.2]
    norm_num
    ring_nf
  · rw [R01.ladder_series hα0 hα1 _ R01.mono23 R01.ne23 k,
      R01.wt23, show 5 / 2 = 2 by norm_num, R01.Icc12,
      Finset.sum_pair (by norm_num), R01.gam23.1, R01.gam23.2]
    norm_num
    ring_nf
  · rw [R01.ladder_series hα0 hα1 _ (R01.mono_single 3 2 (by norm_num)) (R01.single_ne 3 2 two_ne_zero) k,
      R01.wt_single, show 2 * 3 / 2 = 3 by norm_num, R01.Icc13,
      Finset.sum_insert (by decide), Finset.sum_pair (by norm_num), R01.gam33.1, R01.gam33.2.1,
      R01.gam33.2.2]
    norm_num
    ring_nf

namespace R01

lemma collapse {M : Type*} [CommRing M] (W N : ℕ) (hN : W < N) (c : ℕ → M) (g : ℕ → ℕ → M) :
    ∑ e ∈ Finset.range N, ∑ d ∈ Finset.range N,
      (if 1 ≤ e ∧ e ≤ W / 2 ∧ d = W - 2 * e then c e else 0) * g e d =
    ∑ e ∈ Finset.Icc 1 (W / 2), c e * g e (W - 2 * e) := by
  classical
  have inner : ∀ e, ∑ d ∈ Finset.range N,
      (if 1 ≤ e ∧ e ≤ W / 2 ∧ d = W - 2 * e then c e else 0) * g e d =
      if 1 ≤ e ∧ e ≤ W / 2 then c e * g e (W - 2 * e) else 0 := by
    intro e
    by_cases hc : 1 ≤ e ∧ e ≤ W / 2
    · rw [if_pos hc, Finset.sum_eq_single (W - 2 * e)]
      · rw [if_pos ⟨hc.1, hc.2, rfl⟩]
      · intro d _ hd; rw [if_neg (by omega), zero_mul]
      · intro h; exfalso; apply h; rw [Finset.mem_range]; omega
    · rw [if_neg hc]
      exact Finset.sum_eq_zero fun d _ => by rw [if_neg (by tauto), zero_mul]
  rw [Finset.sum_congr rfl fun e _ => inner e, ← Finset.sum_filter]
  refine Finset.sum_congr ?_ fun _ _ => rfl
  ext e; simp only [Finset.mem_filter, Finset.mem_range, Finset.mem_Icc]; omega

lemma zladder_general (s : ℕ →₀ ℕ) (hs : IsZetaMono s) : ∃ (N : ℕ) (r : ℕ → ℕ → ℚ), ∀ k : ℕ,
    (∀ n : ℕ, 1 ≤ n → (Q04.Pf k n).coeff s =
      ∑ e ∈ Finset.range N, ∑ d ∈ Finset.range N,
        r e d * (if d ≤ k then (Q04.Pf (k - d) n).coeff (Finsupp.single 2 e) else 0)) ∧
    ∀ α : ℝ, 0 < α → α < Real.pi →
      zM s k α = ∑ e ∈ Finset.range N, ∑ d ∈ Finset.range N,
        (r e d : ℝ) * zM (Finsupp.single 2 e) ((k : ℤ) - d) α := by
  by_cases h0 : s = 0
  · subst h0
    refine ⟨1, fun _ _ => 1, fun k => ⟨fun n _ => ?_, fun _ _ _ => ?_⟩⟩
    · simp
    · simp
  · refine ⟨wt s + 1, fun e d => if 1 ≤ e ∧ e ≤ wt s / 2 ∧ d = wt s - 2 * e then
        (-1 : ℚ) ^ e * (e.factorial : ℚ) * gam e s else 0, fun k => ⟨fun n hn => ?_, fun α hα0 hα1 => ?_⟩⟩
    · rw [collapse (wt s) (wt s + 1) (by omega), R01.proof_5 s hs h0 k n hn]
      refine Finset.sum_congr rfl fun e he => ?_
      have he' := Finset.mem_Icc.mp he
      rw [ite_Pf (wt s) (wt s - 2 * e) e k n (by omega)]
    · simp only [apply_ite (Rat.cast : ℚ → ℝ), Rat.cast_zero]
      rw [collapse (wt s) (wt s + 1) (by omega), ladder_series hα0 hα1 s hs h0 k]
      refine Finset.sum_congr rfl fun e he => ?_
      have he' := Finset.mem_Icc.mp he
      push_cast
      congr 2
      rw [Nat.cast_sub (by omega)]; push_cast; ring

end R01

theorem R01.proof_10 : R01.stmt_10 :=
  ⟨fun α hα0 hα1 => R01.allSummable hα0 hα1,
   fun k α hα0 hα1 => (R01.proof_7 α hα0 hα1).2 k,
   R01.zladder_general,
   fun k α hα0 hα1 => (R01.proof_8.2 α hα0 hα1).2 k⟩
