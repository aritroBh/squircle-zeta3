import Mathlib
import Zeta3.Wave10.S07_HexagonFormula
import Zeta3.Wave1.T08_PolygonSums
import Zeta3.Wave10.S03_FieldRule
import Zeta3.Wave4.X05_DiamondA4A5
/-! Statements stmt_2, stmt_10, stmt_11, stmt_12, stmt_13 of this audited block are in `P07_Hexagon5Rest.lean` (target P07); the header below still describes all of them. -/
/-!
# P03 -- `P^{(5)}(1)` of every hexagon (thm:hexagon5)

Target id: P03.
Paper labels: `thm:hexagon5` (main.tex 1783-1790, Theorem "$P^{(5)}(1)$ of every hexagon") and
`prop:order5` (main.tex 1755-1773, Proposition "Order five", for `Z_e^{(5)}`, `Z_e^{(2,3)}`).
Detailed proof: paper/ieee/proofs/sec_g9_hex5.tex (Theorem 26 = `g9:thm:hexagon5`, lines 12-19;
"What we prove" 24; "What the words mean" 27; first sentence 41-45; Lemma `g9:lem:shape5` 50-56;
Lemma `g9:lem:Jseries5` 65-75; the formal elements 77; Lemma `g9:lem:hpl` 79-87; Lemma
`g9:lem:nocst` 100-106; Lemma `g9:lem:precise5` 110-124; census 127; Lemma `g9:lem:z5` 129-135;
M6 137; "Proof of Theorem 26" 143). Order three: sec_g4_hexagon.tex, Lean S05/S07.

## Quoted verbatim

main.tex 1783-1790 (Theorem thm:hexagon5):
  For every non-degenerate centrally symmetric convex hexagon and all weights $w_j>0$, with the
  notation and the logarithms of Theorem~\ref{thm:hexagon}, $P^{(5)}(1)=2\sum_eE^{(5)}_e$, summed
  over three edges, one from each pair of opposite edges. With $x=1-r$ and $\ell=-\log r$, each
  $E^{(5)}_e$ is an explicit polynomial, of weight at most five, in $\zeta(2)$, $\zeta(3)$,
  $\zeta(5)$, the logarithms of the corner term values $A_0$, $A_1$, $B$ (each divided by its
  weight) and the thirteen harmonic polylogarithms $H_w(x)$
    \ell,\ \Li_2,\ \Li_3,\ H_{011},\ \Li_4,\ H_{0011},\ H_{0111},
    \Li_5,\ H_{00011},\ H_{00101},\ H_{00111},\ H_{01011},\ H_{01111},
  with coefficients that are polynomials in $A_1$, $B$, $r$, $q_0$, $|V_0|^2$, $L^2$ and $1/x$,
  divided by $L^9$ (the poles at $x=0$ cancel; at $r=1$, $E^{(5)}_e$ is the limit as $x\to0$).
  In this form the coefficients of $\zeta(5)$ and $\zeta(2)\zeta(3)$ are exactly $Z_e^{(5)}$ and
  $Z_e^{(2,3)}$ of Proposition~\ref{prop:order5}, for every $r$ (as coefficients of the number,
  this assumes that the constants are linearly independent).

main.tex 1757-1765 (prop:order5, the two laws; `\bar\xi=\frac{\xi_0L_0+\xi_1L_1}L`,
`\xi_+=\xi_0+\xi_1`, `\xi_\times=\xi_0\xi_1`):
  Z_e^{(5)}&=360\,\xi_\times\frac{h^2}L\Bigl[20\bar\xi^3-32\bar\xi^2+15\bar\xi-2
   -4(5\bar\xi-4)\,\xi_\times\frac{L_0L_1}{L^2}+\frac{h^2}{L^2}\Bigl(5\xi_\times\xi_++8\xi_+^2
   -28\xi_\times-15(\xi_+^2-3\xi_\times)\bar\xi\Bigr)\Bigr],
  Z_e^{(2,3)}&=120\,\xi_\times\frac{h^2}L\Bigl[12\bar\xi^2-12\bar\xi+2
   -12(5\bar\xi-4)\,\xi_\times\frac{L_0L_1}{L^2}+\frac{h^2}{L^2}\Bigl(15\xi_\times\xi_+
   -3\xi_+^2+24\xi_\times-45\xi_\times\bar\xi\Bigr)\Bigr].
  (Lean: `T08.Z5Data`, `T08.Z23Data`, imported; for the hexagon edge `ξ_0 = A_0`, `ξ_1 = A_1`,
  supplement line 134: "with $\xi_0=A_0$, $\xi_1=A_1$".)

Supplement 9: "The \emph{harmonic polylogarithms} are defined for words $w$ in the letters $0,1$
that end in $1$ by $H_1(x)=-\log(1-x)$, $H_{0w}(x)=\int_0^xH_w(s)\frac{\dd s}s$,
$H_{1w}(x)=\int_0^xH_w(s)\frac{\dd s}{1-s}$ ... For $x=1-r$, $H_1(x)=-\log r=\ell$. The
\emph{weight} of $\zeta(i)$ is $i$, of $H_w$ is the length of $w$, of
$\mathrm l_0,\mathrm l_1,\mathrm l_B$ is $1$."

Supplement 21: "$E_e^{(5)}=5!\int_0^1E_5\dd t$ is the contribution of the edge $e$ to
$P^{(5)}(1)$, with $E_5$ the $d^5$ coefficient of the length element of $e$".

Supplement 27: "there are $307$ distinct monomials
$M_k=x^{e_k}\prod_i(H_{w_i}(x))^{\epsilon_{ki}}\,\fz_2^{a_k}\fz_3^{b_k}\fz_5^{c_k}$ in the $13$
functions and in $\zeta(2),\zeta(3),\zeta(5)$, and $307$ coefficients $C_k$ with
$L^9C_k\in\QQ[A_1,B,\kappa,q_0,|V_0|^2,L^2,\mathrm l_0,\mathrm l_1,\mathrm l_B]$
($A_0=1-B-B\kappa$, $\kappa=-x$), such that $E^{(5)}_e=120\sum_kC_kM_k$ for every $r>0$. ...
\emph{``Weight at most five''} is the weight defined above, with the coefficients $C_k$ of weight
$0$ except for the powers of $\mathrm l_0,\mathrm l_1,\mathrm l_B$."

Supplement 43 (first sentence):
  P^{(5)}(1)=5!\sum_{\text{6 edges}}\int_0^1E_5\dd t=2\sum_{e}E_e^{(5)},\qquad
  E_e^{(5)}:=5!\int_0^1E_5\dd t .

Supplement 51 (Lemma g9:lem:shape5): "$E_5=\sum c_{n,a,b,c}\,t^nX^aY^bW^c$, a finite sum with
$a+b+c\le5$, where $c_{n,a,b,c}\in L^{-9}\QQ[B,\kappa,q_0,|V_0|^2,L^2,\mathrm l_0,\mathrm l_1,
\mathrm l_B]$ (with $A_1=1-B$, $A_0=1-B-B\kappa$), and the degree of $c_{n,a,b,c}$ in
$\mathrm l_0,\mathrm l_1,\mathrm l_B$ is at most $5-a-b-c$." (`X=\log t`, `Y=\log(1-t)`,
`W=\log(1+\kappa t)`, line 48.)

Supplement 66-70 (Lemma g9:lem:Jseries5): "For $|x|<1$, $c\ge1$ and $a+b+c\le5$ (so $a+b\le4$),
  J(n;a,b,c)=\sum_{q\ge c}x^q\lambda_c(q)\,I_{n+q}(a,b),\qquad \lambda_c(q)=\frac{(-1)^cc!}q\,
  e_{c-1}\bigl(1,\tfrac12,\dots,\tfrac1{q-1}\bigr),
with $I_N(a,b)=\int_0^1s^N\log^as\log^b(1-s)\dd s$ ... and $e_j$ the elementary symmetric
polynomial. The function $J(n;a,b,c)$ is analytic in $x\in\CC\setminus[1,\infty)$, and
$J=O(x^c)$ as $x\to0$." Here `J(n;a,b,c)=\int_0^1t^nX^aY^b\log^c(1-xt)\dd t` (line 63).
In the proof (line 74): "$-\frac vq\exp(-\sum_{i\ge1}\frac{v^i}iH_{q-1}^{(i)})$ ... whose
coefficient of $v^c$ is ... the printed $\lambda_c(q)/c!$", i.e. `log^c(1-z) = Σ_q λ_c(q) z^q`.

Supplement 77: "$\Ihat(N,0;a,b)$ is, for $a+b\le5$, a polynomial in $\fz_2,\fz_3,\fz_5$, in
$\frac1{N+1}$ and in $H^{(i)}_{N+1}$, $i\le5$, with rational coefficients"; Lemma
g9:lem:beta5(a) (used in the proof of g9:lem:nocst): its monomials "have weight $\le a+b$".

Supplement 80-83 (Lemma g9:lem:hpl):
  (a) For $s_1,\dots,s_k\ge1$ and $|x|<1$, $\Li_{s_1,\dots,s_k}(x):=\sum_{n_1>\dots>n_k\ge1}
  \frac{x^{n_1}}{n_1^{s_1}\cdots n_k^{s_k}}=H_w(x)$ with $w=0^{s_1-1}1\cdots0^{s_k-1}1$.
  (b) $H_uH_v=\sum_{w\in u\sqcup\!\sqcup\,v}H_w$ (shuffle product) for words $u,v$ ending in $1$.
  (c) (Stuffle.) $\sum_{n\ge1}\frac{x^n}{n^a}\prod_j(H_n^{(j)})^{e_j}$ is a rational combination
  of $\Li_{s_1,\dots,s_k}(x)$ of weights $\le a+\sum je_j$ ...
  (d) Every $H_w$ with $|w|\le5$ is a polynomial with rational coefficients in the $13$ functions
  of the theorem. The Lyndon words ending in $1$ of length $\le5$ are $1$; $01$; $001$, $011$;
  $0001$, $0011$, $0111$; $00001$, $00011$, $00101$, $00111$, $01011$, $01111$.

Supplement 101 (Lemma g9:lem:nocst): "Let $\widehat J(n;a,b,c)$ be the closed form of S4, an
element of $\QQ(x)[\text{13 generators}][\fz_2,\fz_3]$. Every monomial of $\widehat J$ in
$\fz_2,\fz_3,\fz_5$ has weight $\le a+b\le4$. In particular $\widehat J$ contains no monomial
$\fz_5$ and no monomial $\fz_2\fz_3$." (S4, line 96: "$210$ closed forms $J(n;a,b,c)$,
polynomials in $x^{\pm1}$, the $13$ generators and $\fz_2,\fz_3$", keys `n = 0,…,5`.)

Supplement 111-119 (Lemma g9:lem:precise5): "... collecting terms gives
$E_e^{(5)}=120\sum_{k=1}^{307}C_kM_k$ with distinct monomials $M_k$ in $x$, the $13$ functions
and $\fz_2,\fz_3,\fz_5$, with $L^9C_k\in\QQ\bigl[A_1,B,\kappa,q_0,|V_0|^2,L^2,\mathrm l_0,
\mathrm l_1,\mathrm l_B\bigr]$, and every term of $C_kM_k$ of weight $\le5$. The identity holds
for every $r>0$ (at $r=1$ as a limit)."  Proof, line 123: "keep
$B,q_0,L,|V_0|^2,\mathrm l_0,\mathrm l_1,\mathrm l_B$ fixed and let $\kappa$ be the variable".

Supplement 130 (Lemma g9:lem:z5): "In the sum of Lemma~\ref{g9:lem:precise5}, the coefficient of
the monomial $\fz_5$ is $Z_e^{(5)}/120$ and the coefficient of $\fz_2\fz_3$ is
$Z_e^{(2,3)}/120$, for every $r>0$."  M6, line 137: "the coefficient of $\fz_5$ is the monomial
with empty generator list and $x^0$, the coefficient of $\fz_2\fz_3$ likewise ...; these are the
only two monomials of the table that contain $\fz_5$ or $\fz_2\fz_3$."

## Conventions (how the paper's objects are modelled)

* The hexagon, weights, edges, edge data, `P(p)` and the length element are those of S05
  (imported, audited block): `S05.Hex`, `S05.IsHex`, `S05.Perim`, `S05.Gl`, `S05.A0`, `S05.A1`,
  `S05.Bv`, `S05.rr`, `S05.q0`, `S05.V02`, `S05.Lsq`, `S05.Len`, `S05.l0`, `S05.l1`, `S05.lB`,
  `S05.evars` (the eight edge variables `B, r, q_0, V_{02}, L², l_0, l_1, l_B`, in this order),
  `S05.Jint`, `S05.Iint`. Since `A_1 = 1 - B` and `κ = r - 1 = -x`, polynomials in
  `A_1, B, κ, …` are the same as polynomials in `B, r, …` (rational coefficients here).
* `E^{(5)}_e` (`Ee5`) is `∫_0^1 ∂_p^5 |x_p'(t)| |_{p=1} dt`, which is `5!∫_0^1 E_5 dt` since
  `E_5` is the `d⁵` coefficient of `G(t,d) = |x_p'(t)|`, `d = p - 1` (exactly as `S05.Ee` at
  order three). Integrability of the integrand is stated (stmt_1).
* The harmonic polylogarithms are `S03.H` (imported, audited): `H_∅ = 1`,
  `H_{a w}(x) = ∫_0^x f_a(t) H_w(t) dt` with `f_0 = 1/t`, `f_1 = 1/(1-t)`; so `H_{[1]}(x) =
  -log(1-x) = ℓ`, exactly the supplement's definition. On `x < 1` all defining integrals are
  genuine. The thirteen functions are `S03.H (gen i)`, `i : Fin 13`, in the order of the theorem.
* `ζ(n)` is `X05.z n = (riemannZeta n).re` (imported, audited), for `n = 2, 3, 5`.
* Nested sums `Σ_{M ≥ i_1 > ⋯ > i_d ≥ 1} Π i_r^{-s_r}` are `S03.Z` with positive entries;
  `H^{(j)}_n` is `S03.hn j n`.
* The formal representative of the theorem (the "table") is modelled as a finitely supported map
  `T : Key →₀ ℚ[8 edge variables]`, `Key = ℕ × (Fin 13 →₀ ℕ) × (ℕ × ℕ × ℕ)`: the key
  `(e, ε, (a,b,c))` is the monomial `M = x^{-e} Π_i H_{w_i}(x)^{ε_i} ζ(2)^a ζ(3)^b ζ(5)^c`
  (distinct keys = distinct monomials; the theorem says "polynomials in ... and $1/x$", so only
  non-positive powers of `x`), and the coefficient is `C = T(key)(edge data)/L⁹`. The weight of
  a term is the `l`-degree of the coefficient monomial plus `Σ ε_i|w_i| + 2a + 3b + 5c`.
  The table itself (`scripts/g9_hex5_table.txt`, 307 lines, machine output) is quantified
  existentially; the two lines the paper prints (through prop:order5) are fixed (stmt_12).
* "In this form the coefficient of ζ(5)" is formalised, as the supplement says (lines 24, 27,
  137), as the coefficient of the monomial `ζ(5)` (key `(0, 0, (0,0,1))`) of the formal
  representative, together with M6's statement that no other monomial contains `ζ(5)` or
  `ζ(2)ζ(3)`. No linear independence is claimed.
* "at $r=1$, $E^{(5)}_e$ is the limit as $x\to0$": `Tendsto` of the formula along `x → 0`,
  `x ≠ 0`, with the other seven edge variables fixed (supplement line 123).

## Already proved elsewhere (imported, not restated)

* `ΣZ_e^{(5)} = C_m^{(5)}`, `ΣZ_e^{(2,3)} = C_m^{(2,3)}` for the regular `2m`-gon: T08 stmt_5,
  stmt_6. The setting of the hexagon (every hexagon with every weight is a `K_1`; the three
  terms on an edge; `0 < B < 1`, `A_0, A_1, r > 0`): S05 stmt_5 to stmt_8.

## NOT COVERED

* The 307 coefficients of the table and the 336 coefficients `c_{n,a,b,c}` of `E_5` (machine
  steps S1, S5-S6; only machine output, ~5 MB): quantified existentially (stmt_2, stmt_11);
  only the two printed blocks (`ζ(5)`, `ζ(2)ζ(3)`) are fixed. The census (307 monomials,
  x-exponents `-6..0`, which generators occur) is not stated.
* The 210 closed forms `Ĵ(n;a,b,c)` of S4 (machine output): only their existence in the stated
  shape, with the weight bound of g9:lem:nocst, is stated (stmt_10). That `Ĵ` is the specific
  output of S4 cannot be expressed.
* The complex analyticity in Lemma g9:lem:Jseries5 (on `ℂ ∖ [1,∞)`): stated only as real
  analyticity on `(-∞,1)` (stmt_4); the complex continuation is a proof device.
* Lemma g9:lem:hpl(d)'s determinants `1,-2,-6,-96,5760` (machine step S2'): only the conclusion
  (every word is a polynomial in the 13 generators) is stated (stmt_9). Lemma (c) is stated for
  `a ≥ 1` (for `a = 0` the sum is not a combination of `Li_s`; the lemma tacitly has `a ≥ 1`).
* The remark "eleven of the thirteen functions are Nielsen polylogarithms $S_{n,p}=H_{0^n1^p}$"
  is a statement about the shape of words (and the supplement says the irreducibility of
  `H_{00101}`, `H_{01011}` is neither cited nor used); not stated.
* Corollary g9:cor:r1 ("the local value of that edge", M2 of Section g9:sec:order5) and the
  regular-hexagon closed form: the local value is not defined here; only the `r = 1` limit
  (stmt_11) is stated.
* The `r = 2` corollary (the thirteen values at `x = -1` as alternating MZVs, the coefficients
  `-20187/2000 √5`, `-858/125 √5`, the "far kink lines add"): basis-dependent and relying on
  the alternating-MZV relations; not stated. The numerical checks (golden hexagon, six
  hexagons, 76 digits) are not statements.
* The intermediate formula of g9:lem:precise5 in terms of `c_{n,a,b,c}`, `Î` and `Ĵ` is not
  stated separately; it is the composition of stmt_2, stmt_5, stmt_10 and is subsumed by
  stmt_11.
* The g9:lem:z5 proof's identification of the coefficient as `24∫Σ g_{ab}` (local `g_{ab}` of
  Section g9:sec:order5): not stated; only its conclusion (stmt_12).

## Numerical check before freezing

Own mpmath code (40 digits). Two random non-regular hexagons with random weights (corners
`(1.19437,-0.136084), (0.947446,1.63997), (-0.255319,0.974098)`, weights `2.09321, 0.463807,
2.19403`, ratios `r ≈ 1.1887, 0.5406, 1.5561`; and a second one with `r ≈ 0.832, 0.799, 1.504`):
for each edge `Ee5` computed from the definition here (`∂_p^5|x_p'(t)|` by numerical
differentiation, integrated in `t`) agrees with the 307-line table (`scripts/g9_hex5_table.txt`,
the thirteen `H_w` summed as nested series) to `2·10⁻³¹`; the `ζ(5)` and `ζ(2)ζ(3)` lines of the
table times `120` equal `T08.Z5Data`, `T08.Z23Data` at `ξ_0 = A_0`, `ξ_1 = A_1` to `10⁻³⁷`.
Sentence 1: for the supplement's hexagon `(1,1), (-27/10,33/10), (-2,1)`, weights `5,1,2`,
`2ΣEe5` from the definition is `-1474.678480049841208073026`, the referee's `P^{(5)}(1)` to 25
digits. In the table the largest weight (generators + `ζ` + `l`-degree) is `5` and only odd
powers of `L` occur. Also: the `λ_c` series for `c = 1..5` (`10⁻³¹`), the `J` series for
`(n;a,b,c) = (2;1,1,2)`, `x = 0.4` (`10⁻³²`), (a) for `0101`, `011` against iterated integrals,
(b) for `01 ш 1` at `x = -1.5` (`2·10⁻³¹`), (c) for `Σx^nH_n/n² = Li_{2,1} + Li_3`, and
`I_N(0,1) = -H_{N+1}/(N+1)`.
-/

-- BEGIN STATEMENTS
namespace P03

open MeasureTheory Real Set Filter

/-- `E^{(5)}_e = 5!∫_0^1 E_5 dt = ∫_0^1 ∂_p^5 G(t, p-1)|_{p=1} dt` for edge `k`
(the order-five analogue of `S05.Ee`). -/
noncomputable def Ee5 (H : S05.Hex) (k : Fin 3) : ℝ :=
  ∫ t in (0 : ℝ)..1, iteratedDeriv 5 (fun p => S05.Gl H k p t) 1

/-- Evaluation of a rational polynomial in the eight edge variables
`(B, r, q_0, V_{02}, L², l_0, l_1, l_B)` at the data of edge `k`. -/
noncomputable def evq (H : S05.Hex) (k : Fin 3) (P : MvPolynomial (Fin 8) ℚ) : ℝ :=
  MvPolynomial.eval₂ (algebraMap ℚ ℝ) (S05.evars H k) P

/-- Degree of a monomial of the edge variables in the logarithms `l_0, l_1, l_B`
(variables `5, 6, 7`). -/
def ldeg (m : Fin 8 →₀ ℕ) : ℕ := m 5 + m 6 + m 7

/-- The thirteen words of the theorem (the Lyndon words ending in `1` of length `≤ 5`), in the
order `ℓ, Li₂, Li₃, H₀₁₁, Li₄, H₀₀₁₁, H₀₁₁₁, Li₅, H₀₀₀₁₁, H₀₀₁₀₁, H₀₀₁₁₁, H₀₁₀₁₁, H₀₁₁₁₁`. -/
def gen : Fin 13 → List ℤ :=
  ![[1], [0, 1], [0, 0, 1], [0, 1, 1], [0, 0, 0, 1], [0, 0, 1, 1], [0, 1, 1, 1],
    [0, 0, 0, 0, 1], [0, 0, 0, 1, 1], [0, 0, 1, 0, 1], [0, 0, 1, 1, 1], [0, 1, 0, 1, 1],
    [0, 1, 1, 1, 1]]

/-- A word in the letters `0, 1`. -/
def IsWord01 (w : List ℤ) : Prop := ∀ a ∈ w, a = 0 ∨ a = 1

/-- The word `0^{s_1-1}1 0^{s_2-1}1 ⋯ 0^{s_k-1}1` of the index `(s_1, …, s_k)`. -/
def wordOf (s : List ℕ) : List ℤ := s.flatMap (fun n => List.replicate (n - 1) 0 ++ [1])

/-- The shuffle product of two words, as the list of all shuffles (with multiplicity). -/
def shuffle : List ℤ → List ℤ → List (List ℤ)
  | [], v => [v]
  | u, [] => [u]
  | a :: u, b :: v =>
      (shuffle u (b :: v)).map (a :: ·) ++ (shuffle (a :: u) v).map (b :: ·)

/-- `λ_c(q) = (-1)^c c!/q · e_{c-1}(1, 1/2, …, 1/(q-1))` (`e_j` the elementary symmetric
polynomial). -/
noncomputable def lamc (c q : ℕ) : ℝ :=
  (-1 : ℝ) ^ c * (c.factorial : ℝ) / (q : ℝ) *
    ((Finset.Icc 1 (q - 1)).val.map (fun i : ℕ => (1 : ℝ) / (i : ℝ))).esymm (c - 1)

/-- A key `(e, ε, (a, b, c))` of the formal representative: the monomial
`x^{-e} Π_i H_{w_i}(x)^{ε_i} ζ(2)^a ζ(3)^b ζ(5)^c`. -/
abbrev Key := ℕ × (Fin 13 →₀ ℕ) × (ℕ × ℕ × ℕ)

/-- Weight of a key: `Σ_i ε_i |w_i| + 2a + 3b + 5c` (`x` has weight `0`). -/
def kwt (κ : Key) : ℕ :=
  (κ.2.1.sum fun i e => e * (gen i).length) + 2 * κ.2.2.1 + 3 * κ.2.2.2.1 + 5 * κ.2.2.2.2

/-- Value of the monomial of key `κ` at `x`. -/
noncomputable def mval (κ : Key) (x : ℝ) : ℝ :=
  (x ^ κ.1)⁻¹ * (κ.2.1.prod fun i e => S03.H (gen i) x ^ e) *
    X05.z 2 ^ κ.2.2.1 * X05.z 3 ^ κ.2.2.2.1 * X05.z 5 ^ κ.2.2.2.2

/-- The formula `120 Σ_k C_k M_k` of the theorem at the edge variables `v`
(`v 1 = r`, `x = 1 - r`, `L = √(v 4)`), with `C_k = T(κ_k)(v)/L⁹`. -/
noncomputable def fval (T : Key →₀ MvPolynomial (Fin 8) ℚ) (v : Fin 8 → ℝ) : ℝ :=
  120 * T.sum fun κ P =>
    MvPolynomial.eval₂ (algebraMap ℚ ℝ) v P / Real.sqrt (v 4) ^ 9 * mval κ (1 - v 1)

/-- The key of the monomial `ζ(5)`. -/
def key5 : Key := (0, 0, (0, 0, 1))

/-- The key of the monomial `ζ(2)ζ(3)`. -/
def key23 : Key := (0, 0, (1, 1, 0))

/-- `Z_e^{(5)}` of prop:order5 for edge `k` of the hexagon, with `ξ_0 = A_0`, `ξ_1 = A_1`
(supplement line 134) and `h, L, L_0` from the corners (T08's edge data). -/
noncomputable def Ze5 (H : S05.Hex) (k : Fin 3) : ℝ :=
  let V0 := S05.toPair (S05.eV0 H k)
  let V1 := S05.toPair (S05.eV1 H k)
  T08.Z5Data (S05.A0 H k) (S05.A1 H k) (T08.edgeH V0 V1) (T08.edgeL V0 V1) (T08.edgeL0 V0 V1)

/-- `Z_e^{(2,3)}` of prop:order5 for edge `k`, with `ξ_0 = A_0`, `ξ_1 = A_1`. -/
noncomputable def Ze23 (H : S05.Hex) (k : Fin 3) : ℝ :=
  let V0 := S05.toPair (S05.eV0 H k)
  let V1 := S05.toPair (S05.eV1 H k)
  T08.Z23Data (S05.A0 H k) (S05.A1 H k) (T08.edgeH V0 V1) (T08.edgeL V0 V1) (T08.edgeL0 V0 V1)

/-- First sentence of thm:hexagon5 (supplement 41-45): `P` is five times continuously
differentiable at `p = 1`, the edge integrands are integrable, and
`P^{(5)}(1) = 2 Σ_e E^{(5)}_e` over the three edges `C_0C_1`, `C_1C_2`, `C_2(-C_0)`. -/
def stmt_1 : Prop :=
  ∀ H : S05.Hex, S05.IsHex H →
    ContDiffAt ℝ 5 (S05.Perim H) 1 ∧
    (∀ k : Fin 3,
      IntervalIntegrable (fun t => iteratedDeriv 5 (fun p => S05.Gl H k p t) 1) volume 0 1) ∧
    iteratedDeriv 5 (S05.Perim H) 1 = 2 * ∑ k : Fin 3, Ee5 H k

/-- Proof of g9:lem:Jseries5: `log^c(1-z) = Σ_{q≥1} λ_c(q) z^q` for `c ≥ 1`, `|z| < 1`
(the terms with `q < c` vanish). -/
def stmt_3 : Prop :=
  ∀ c : ℕ, 1 ≤ c → ∀ z : ℝ, |z| < 1 →
    HasSum (fun q : ℕ => lamc c (q + 1) * z ^ (q + 1)) (Real.log (1 - z) ^ c)

/-- Lemma g9:lem:Jseries5: for `c ≥ 1`, `a + b + c ≤ 5`, the integrand of `J(n;a,b,c)` is
integrable for every `x < 1`, `J` is (real-)analytic on `x < 1` and `O(x^c)` at `0`, and for
`|x| < 1`, `J(n;a,b,c) = Σ_{q≥1} x^q λ_c(q) I_{n+q}(a,b)`. -/
def stmt_4 : Prop :=
  ∀ n a b c : ℕ, 1 ≤ c → a + b + c ≤ 5 →
    (∀ x : ℝ, x < 1 →
      IntervalIntegrable
        (fun t : ℝ => t ^ n * Real.log t ^ a * Real.log (1 - t) ^ b * Real.log (1 - x * t) ^ c)
        volume 0 1) ∧
    AnalyticOnNhd ℝ (S05.Jint n a b c) (Iio 1) ∧
    (S05.Jint n a b c) =O[nhds 0] (fun x : ℝ => x ^ c) ∧
    ∀ x : ℝ, |x| < 1 →
      HasSum (fun q : ℕ => x ^ (q + 1) * lamc c (q + 1) * S05.Iint (n + q + 1) a b)
        (S05.Jint n a b c x)

/-- Supplement line 77 with Lemma g9:lem:beta5(a): for `a + b ≤ 5`, `I_N(a,b)` is a polynomial
with rational coefficients in `1/(N+1)`, `H^{(i)}_{N+1}` (`i = 1,…,5`) and `ζ(2), ζ(3), ζ(5)`,
whose monomials in the `ζ`'s have weight `≤ a + b`. Variables: `0 ↦ 1/(N+1)`,
`i ↦ H^{(i)}_{N+1}` (`i = 1..5`), `6, 7, 8 ↦ ζ(2), ζ(3), ζ(5)`. -/
def stmt_5 : Prop :=
  ∀ a b : ℕ, a + b ≤ 5 →
    ∃ P : MvPolynomial (Fin 9) ℚ,
      (∀ m ∈ P.support, 2 * m 6 + 3 * m 7 + 5 * m 8 ≤ a + b) ∧
      ∀ N : ℕ, S05.Iint N a b =
        MvPolynomial.eval₂ (algebraMap ℚ ℝ)
          (![1 / ((N : ℝ) + 1), (S03.hn 1 (N + 1) : ℝ), (S03.hn 2 (N + 1) : ℝ),
            (S03.hn 3 (N + 1) : ℝ), (S03.hn 4 (N + 1) : ℝ), (S03.hn 5 (N + 1) : ℝ),
            X05.z 2, X05.z 3, X05.z 5]) P

/-- Lemma g9:lem:hpl(a): for `s_1, …, s_k ≥ 1` and `|x| < 1`,
`Li_{s_1,…,s_k}(x) = Σ_{n_1 > ⋯ > n_k ≥ 1} x^{n_1}/(n_1^{s_1}⋯n_k^{s_k}) = H_w(x)`,
`w = 0^{s_1-1}1⋯0^{s_k-1}1` (the inner sum is `S03.Z (s_2,…,s_k) (n_1 - 1)`). -/
def stmt_6 : Prop :=
  ∀ (s1 : ℕ) (s : List ℕ), 1 ≤ s1 → (∀ si ∈ s, 1 ≤ si) → ∀ x : ℝ, |x| < 1 →
    HasSum
      (fun n : ℕ => x ^ (n + 1) / ((n : ℝ) + 1) ^ s1 *
        (S03.Z (s.map (fun si : ℕ => (si : ℤ))) n : ℝ))
      (S03.H (wordOf (s1 :: s)) x)

/-- Lemma g9:lem:hpl(b): the shuffle product `H_u H_v = Σ_{w ∈ u ш v} H_w` for words `u, v` in
the letters `0, 1` ending in `1`, at every `x < 1`. -/
def stmt_7 : Prop :=
  ∀ u v : List ℤ, IsWord01 u → IsWord01 v → S03.Regular u → S03.Regular v →
    ∀ x : ℝ, x < 1 →
      S03.H u x * S03.H v x = ((shuffle u v).map (fun w => S03.H w x)).sum

/-- Lemma g9:lem:hpl(c) (stuffle), for `a ≥ 1`: `Σ_{n≥1} x^n/n^a Π_{j=1}^J (H^{(j)}_n)^{e_j}` is
a rational combination of `Li_{s}(x) = H_{w(s)}(x)` (by (a)) of weights
`|s| ≤ a + Σ_j j e_j`, for `|x| < 1`. -/
def stmt_8 : Prop :=
  ∀ (a : ℕ), 1 ≤ a → ∀ (J : ℕ) (e : ℕ → ℕ),
    ∃ terms : List (ℚ × List ℕ),
      (∀ τ ∈ terms, τ.2 ≠ [] ∧ (∀ si ∈ τ.2, 1 ≤ si) ∧
        τ.2.sum ≤ a + ∑ j ∈ Finset.Icc 1 J, j * e j) ∧
      ∀ x : ℝ, |x| < 1 →
        HasSum
          (fun n : ℕ => x ^ (n + 1) / ((n : ℝ) + 1) ^ a *
            ∏ j ∈ Finset.Icc 1 J, (S03.hn j (n + 1) : ℝ) ^ e j)
          ((terms.map fun τ => (τ.1 : ℝ) * S03.H (wordOf τ.2) x).sum)

/-- Lemma g9:lem:hpl(d): every `H_w`, `w` a word in `0, 1` ending in `1` with `|w| ≤ 5`, is a
polynomial with rational coefficients in the thirteen functions, on all of `x < 1`; the
polynomial is homogeneous of weight `|w|` (the grading of the proof). -/
def stmt_9 : Prop :=
  ∀ w : List ℤ, IsWord01 w → S03.Regular w → w.length ≤ 5 →
    ∃ P : MvPolynomial (Fin 13) ℚ,
      (∀ m ∈ P.support, (m.sum fun i e => e * (gen i).length) = w.length) ∧
      ∀ x : ℝ, x < 1 →
        S03.H w x = MvPolynomial.eval₂ (algebraMap ℚ ℝ) (fun i => S03.H (gen i) x) P

end P03
-- END STATEMENTS


namespace P03.H5
open Real MeasureTheory Set
open scoped Interval

lemma intInt_bdd (g : ℝ → ℝ) (R r : ℝ) (hc : ContinuousOn g (Icc (-R) r \ {0})) (K : ℝ)
    (hK : ∀ t ∈ Icc (-R) r, |g t| ≤ K) (x : ℝ) (hx : x ∈ Icc (-R) r) (h0 : (0:ℝ) ∈ Icc (-R) r) :
    IntervalIntegrable g volume 0 x := by
  have hsub : Ι (0 : ℝ) x ⊆ Icc (-R) r := by
    intro t ht
    have := uIoc_subset_uIcc ht
    rw [mem_uIcc] at this
    constructor <;> rcases this with h | h <;> linarith [h.1, h.2, hx.1, hx.2, h0.1, h0.2]
  rw [intervalIntegrable_iff]
  refine IntegrableOn.of_bound measure_Ioc_lt_top ?_ K ?_
  · have := (hc.mono (sdiff_subset_sdiff_left hsub)).aestronglyMeasurable
      (measurableSet_uIoc.diff (measurableSet_singleton 0)) (μ := volume)
    rwa [Measure.restrict_congr_set (sdiff_null_ae_eq_self (measure_singleton 0))] at this
  · exact ae_restrict_of_forall_mem measurableSet_uIoc (fun t ht => by
      rw [Real.norm_eq_abs]; exact hK t (hsub ht))

lemma f0e (t : ℝ) : S03.f 0 t = 1 / t := by simp [S03.f]
lemma f1e (t : ℝ) : S03.f 1 t = 1 / (1 - t) := by simp [S03.f]

lemma f_contOn' (a : ℤ) (ha : a = 0 ∨ a = 1) : ContinuousOn (S03.f a) (Iio 1 \ {0}) := by
  intro t ht
  obtain ⟨h1, h0⟩ := ht
  have h0 : t ≠ 0 := h0
  have h1 : t < 1 := h1
  apply ContinuousAt.continuousWithinAt
  rcases ha with rfl | rfl
  · rw [funext f0e]
    exact continuousAt_const.div continuousAt_id h0
  · rw [funext f1e]
    exact continuousAt_const.div (continuousAt_const.sub continuousAt_id)
      (ne_of_gt (show (1:ℝ) - t > 0 by linarith))

lemma f1_bound (r : ℝ) (hr : r < 1) (t : ℝ) (ht : t ≤ r) : |S03.f 1 t| ≤ 1 / (1 - r) := by
  rw [funext f1e]
  have hr1 : 0 < 1 - r := by linarith
  rw [abs_of_pos (by apply div_pos one_pos; linarith)]
  apply one_div_le_one_div_of_le hr1; linarith

lemma IsWord01_tail {a : ℤ} {w : List ℤ} (h : IsWord01 (a :: w)) : IsWord01 w := fun b hb => h b (by simp [hb])

/-- Continuity of `H_w` on `(-∞,1)`, `|H_w(t)| ≤ C|t|` and boundedness on `[-R,r]`. -/
lemma HG (w : List ℤ) (hw : IsWord01 w) (hr : S03.Regular w) :
    ContinuousOn (S03.H w) (Iio 1) ∧
    (w ≠ [] → ∀ R r : ℝ, 0 ≤ R → 0 ≤ r → r < 1 → ∃ C : ℝ, 0 ≤ C ∧
      ∀ t ∈ Icc (-R) r, |S03.H w t| ≤ C * |t|) ∧
    (∀ R r : ℝ, 0 ≤ R → 0 ≤ r → r < 1 → ∃ B : ℝ, 0 ≤ B ∧ ∀ t ∈ Icc (-R) r, |S03.H w t| ≤ B) := by
  induction w with
  | nil =>
    refine ⟨continuousOn_const, fun h => absurd rfl h, fun R r _ _ _ => ⟨1, by norm_num, fun t _ => ?_⟩⟩
    simp [S03.H]
  | cons a w ih =>
    have hw' := IsWord01_tail hw
    have hr' := S03.reg_tail a w hr
    obtain ⟨ihc, ihC, ihB⟩ := ih hw' hr'
    have ha := hw a (by simp)
    set g := fun t => S03.f a t * S03.H w t with hg
    have gc : ContinuousOn g (Iio 1 \ {0}) := (f_contOn' a ha).mul (ihc.mono sdiff_subset)
    have gK : ∀ R r : ℝ, 0 ≤ R → 0 ≤ r → r < 1 → ∃ K : ℝ, 0 ≤ K ∧
        ∀ t ∈ Icc (-R) r, |g t| ≤ K := by
      intro R r hR0 hr0 hr1
      rcases ha with rfl | rfl
      · obtain ⟨C, hC0, hC⟩ := ihC (S03.reg_zero_ne w hr) R r hR0 hr0 hr1
        refine ⟨C, hC0, fun t ht => ?_⟩
        simp only [hg, f0e]
        by_cases ht0 : t = 0
        · subst ht0; simp [hC0]
        · rw [abs_mul, abs_div, abs_one]
          have := hC t ht
          have hpos : 0 < |t| := abs_pos.mpr ht0
          rw [div_mul_eq_mul_div, one_mul, div_le_iff₀ hpos]; exact this
      · obtain ⟨B, hB0, hB⟩ := ihB R r hR0 hr0 hr1
        refine ⟨1 / (1 - r) * B, by apply mul_nonneg; apply div_nonneg; norm_num; linarith; exact hB0,
          fun t ht => ?_⟩
        rw [abs_mul]
        exact mul_le_mul (f1_bound r hr1 t ht.2) (hB t ht) (abs_nonneg _)
          (by apply div_nonneg; norm_num; linarith)
    have gI : ∀ R r : ℝ, 0 ≤ R → 0 ≤ r → r < 1 → ∀ x ∈ Icc (-R) r,
        IntervalIntegrable g volume 0 x := by
      intro R r hR0 hr0 hr1 x hx
      obtain ⟨K, _, hK⟩ := gK R r hR0 hr0 hr1
      refine intInt_bdd g R r (gc.mono ?_) K hK x hx ⟨by linarith, hr0⟩
      intro t ⟨⟨h1, h2⟩, h3⟩
      exact ⟨by simp only [mem_Iio]; linarith, h3⟩
    have Hlin : ∀ R r : ℝ, 0 ≤ R → 0 ≤ r → r < 1 → ∃ C : ℝ, 0 ≤ C ∧
        ∀ t ∈ Icc (-R) r, |S03.H (a :: w) t| ≤ C * |t| := by
      intro R r hR0 hr0 hr1
      obtain ⟨K, hK0, hK⟩ := gK R r hR0 hr0 hr1
      refine ⟨K, hK0, fun t ht => ?_⟩
      rw [S03.H_cons]
      have := intervalIntegral.norm_integral_le_of_norm_le_const (a := 0) (b := t) (C := K)
        (f := g) (fun s hs => by
          rw [Real.norm_eq_abs]
          apply hK
          have := uIoc_subset_uIcc hs
          rw [mem_uIcc] at this
          constructor <;> rcases this with h | h <;> linarith [h.1, h.2, ht.1, ht.2])
      rw [Real.norm_eq_abs, sub_zero] at this
      exact this
    refine ⟨?_, fun _ => Hlin, ?_⟩
    · intro x0 hx0
      have hx0 : x0 < 1 := hx0
      set r := (max x0 0 + 1) / 2 with hrdef
      set R := |x0| + 1 with hRdef
      have hm : max x0 0 < 1 := max_lt hx0 one_pos
      have hm0 : 0 ≤ max x0 0 := le_max_right _ _
      have hmx : x0 ≤ max x0 0 := le_max_left _ _
      have hr0 : 0 ≤ r := by positivity
      have hr1 : r < 1 := by rw [hrdef]; linarith
      have hxr : x0 < r := by rw [hrdef]; linarith
      have hR0 : 0 ≤ R := by positivity
      have hxR : -R < x0 := by rw [hRdef]; linarith [neg_abs_le x0]
      have hint : IntervalIntegrable g volume (-R) r := by
        apply IntervalIntegrable.trans (b := 0)
        · exact (gI R r hR0 hr0 hr1 (-R) ⟨le_rfl, by linarith⟩).symm
        · exact gI R r hR0 hr0 hr1 r ⟨by linarith, le_rfl⟩
      have hcont := intervalIntegral.continuousOn_primitive_interval' hint
        (a := 0) (by rw [mem_uIcc]; left; constructor <;> linarith)
      have hmem : uIcc (-R) r ∈ nhds x0 := by
        rw [uIcc_of_le (by linarith)]
        exact Icc_mem_nhds hxR hxr
      exact (hcont.continuousAt hmem).continuousWithinAt
    · intro R r hR0 hr0 hr1
      obtain ⟨C, hC0, hC⟩ := Hlin R r hR0 hr0 hr1
      refine ⟨C * (R + r), mul_nonneg hC0 (by linarith), fun t ht => ?_⟩
      refine (hC t ht).trans (mul_le_mul_of_nonneg_left ?_ hC0)
      rw [abs_le]; constructor <;> linarith [ht.1, ht.2]

lemma integrand_ok' (a : ℤ) (w : List ℤ) (hw : IsWord01 (a :: w)) (hr : S03.Regular (a :: w)) :
    (∀ x : ℝ, x < 1 → IntervalIntegrable (fun t => S03.f a t * S03.H w t) volume 0 x) ∧
    ContinuousOn (fun t => S03.f a t * S03.H w t) (Iio 1 \ {0}) := by
  have hw' := IsWord01_tail hw
  have hr' := S03.reg_tail a w hr
  obtain ⟨ihc, ihC, ihB⟩ := HG w hw' hr'
  have ha := hw a (by simp)
  have gc : ContinuousOn (fun t => S03.f a t * S03.H w t) (Iio 1 \ {0}) :=
    (f_contOn' a ha).mul (ihc.mono sdiff_subset)
  refine ⟨fun x hx => ?_, gc⟩
  set r := max x 0
  set R := |x|
  have hr0 : 0 ≤ r := le_max_right _ _
  have hr1 : r < 1 := max_lt hx one_pos
  have hR0 : 0 ≤ R := abs_nonneg x
  have hxI : x ∈ Icc (-R) r := ⟨neg_abs_le x, le_max_left _ _⟩
  have h0I : (0:ℝ) ∈ Icc (-R) r := ⟨by linarith, hr0⟩
  have gc' : ContinuousOn (fun t => S03.f a t * S03.H w t) (Icc (-R) r \ {0}) := by
    refine gc.mono ?_
    intro t ⟨⟨h1, h2⟩, h3⟩
    exact ⟨by simp only [mem_Iio]; linarith, h3⟩
  rcases ha with rfl | rfl
  · obtain ⟨C, hC0, hC⟩ := ihC (S03.reg_zero_ne w hr) R r hR0 hr0 hr1
    refine intInt_bdd _ R r gc' C (fun t ht => ?_) x hxI h0I
    simp only [f0e]
    by_cases ht0 : t = 0
    · subst ht0; simp [hC0]
    · rw [abs_mul, abs_div, abs_one]
      have := hC t ht
      have hpos : 0 < |t| := abs_pos.mpr ht0
      rw [div_mul_eq_mul_div, one_mul, div_le_iff₀ hpos]; exact this
  · obtain ⟨B, hB0, hB⟩ := ihB R r hR0 hr0 hr1
    refine intInt_bdd _ R r gc' (1 / (1 - r) * B) (fun t ht => ?_) x hxI h0I
    rw [abs_mul]
    exact mul_le_mul (f1_bound r hr1 t ht.2) (hB t ht) (abs_nonneg _)
      (by apply div_nonneg; norm_num; linarith)

lemma H_deriv' (a : ℤ) (w : List ℤ) (hw : IsWord01 (a :: w)) (hr : S03.Regular (a :: w)) (x : ℝ)
    (hx0 : x ≠ 0) (hx1 : x < 1) : HasDerivAt (S03.H (a :: w)) (S03.f a x * S03.H w x) x := by
  obtain ⟨hI, hc⟩ := integrand_ok' a w hw hr
  have hxs : x ∈ Iio (1 : ℝ) \ {0} := ⟨hx1, hx0⟩
  have hopen : IsOpen (Iio (1 : ℝ) \ {0}) := isOpen_Iio.sdiff isClosed_singleton
  exact intervalIntegral.integral_hasDerivAt_right (hI x hx1)
    (hc.stronglyMeasurableAtFilter hopen x hxs) (hc.continuousAt (hopen.mem_nhds hxs))

lemma H_zero (w : List ℤ) (hw : w ≠ []) : S03.H w 0 = 0 := by
  cases w with
  | nil => exact absurd rfl hw
  | cons a w => simp [S03.H_cons]


lemma reg_cons (a : ℤ) (w : List ℤ) (hw : w ≠ []) (hr : S03.Regular w) : S03.Regular (a :: w) := by
  cases w with
  | nil => exact absurd rfl hw
  | cons b w => intro c hc; exact hr c (by simpa [List.getLast?_cons_cons] using hc)

lemma shuffle_len : ∀ u v : List ℤ, ∀ w ∈ shuffle u v, w.length = u.length + v.length
  | [], v, w, hw => by simp [shuffle] at hw; simp [hw]
  | a :: u, [], w, hw => by simp [shuffle] at hw; simp [hw]
  | a :: u, b :: v, w, hw => by
    simp only [shuffle, List.mem_append, List.mem_map] at hw
    rcases hw with ⟨w', hw', rfl⟩ | ⟨w', hw', rfl⟩
    · have := shuffle_len u (b :: v) w' hw'; simp [this]; omega
    · have := shuffle_len (a :: u) v w' hw'; simp [this]; omega

lemma shuffle_w01 : ∀ u v : List ℤ, IsWord01 u → IsWord01 v → ∀ w ∈ shuffle u v, IsWord01 w
  | [], v, _, hv, w, hw => by simp [shuffle] at hw; rwa [hw]
  | a :: u, [], hu, _, w, hw => by simp [shuffle] at hw; rwa [hw]
  | a :: u, b :: v, hu, hv, w, hw => by
    simp only [shuffle, List.mem_append, List.mem_map] at hw
    rcases hw with ⟨w', hw', rfl⟩ | ⟨w', hw', rfl⟩
    · have := shuffle_w01 u (b :: v) (IsWord01_tail hu) hv w' hw'
      intro c hc; simp at hc; rcases hc with rfl | hc
      · exact hu c (by simp)
      · exact this c hc
    · have := shuffle_w01 (a :: u) v hu (IsWord01_tail hv) w' hw'
      intro c hc; simp at hc; rcases hc with rfl | hc
      · exact hv c (by simp)
      · exact this c hc

lemma shuffle_reg : ∀ u v : List ℤ, S03.Regular u → S03.Regular v → ∀ w ∈ shuffle u v,
    S03.Regular w
  | [], v, _, hv, w, hw => by simp [shuffle] at hw; rwa [hw]
  | a :: u, [], hu, _, w, hw => by simp [shuffle] at hw; rwa [hw]
  | a :: u, b :: v, hu, hv, w, hw => by
    simp only [shuffle, List.mem_append, List.mem_map] at hw
    rcases hw with ⟨w', hw', rfl⟩ | ⟨w', hw', rfl⟩
    · have hne : w' ≠ [] := by
        intro h; have := shuffle_len u (b :: v) w' hw'; rw [h] at this; simp at this
      cases u with
      | nil =>
        simp [shuffle] at hw'; subst hw'
        exact reg_cons a _ (by simp) hv
      | cons c u =>
        exact reg_cons a w' hne (shuffle_reg (c :: u) (b :: v) (S03.reg_tail a _ hu) hv w' hw')
    · have hne : w' ≠ [] := by
        intro h; have := shuffle_len (a :: u) v w' hw'; rw [h] at this; simp at this; omega
      cases v with
      | nil =>
        simp [shuffle] at hw'; subst hw'
        exact reg_cons b _ (by simp) hu
      | cons c v =>
        exact reg_cons b w' hne (shuffle_reg (a :: u) (c :: v) hu (S03.reg_tail b _ hv) w' hw')

lemma hasDerivAt_lsum {ι : Type*} (L : List ι) (φ : ι → ℝ → ℝ) (φ' : ι → ℝ) (x : ℝ)
    (h : ∀ i ∈ L, HasDerivAt (φ i) (φ' i) x) :
    HasDerivAt (fun y => (L.map fun i => φ i y).sum) ((L.map φ').sum) x := by
  induction L with
  | nil => simpa using hasDerivAt_const x (0:ℝ)
  | cons i L ih =>
    simp only [List.map_cons, List.sum_cons]
    exact (h i (by simp)).add (ih fun j hj => h j (by simp [hj]))

lemma contOn_lsum {ι : Type*} (L : List ι) (φ : ι → ℝ → ℝ) (s : Set ℝ)
    (h : ∀ i ∈ L, ContinuousOn (φ i) s) :
    ContinuousOn (fun y => (L.map fun i => φ i y).sum) s := by
  induction L with
  | nil => simpa using continuousOn_const
  | cons i L ih =>
    simp only [List.map_cons, List.sum_cons]
    exact (h i (by simp)).add (ih fun j hj => h j (by simp [hj]))

lemma eq_of_deriv (F G : ℝ → ℝ) (d : ℝ → ℝ) (hF : ContinuousOn F (Iio 1))
    (hG : ContinuousOn G (Iio 1)) (hdF : ∀ x, x ≠ 0 → x < 1 → HasDerivAt F (d x) x)
    (hdG : ∀ x, x ≠ 0 → x < 1 → HasDerivAt G (d x) x) (h0 : F 0 = G 0) :
    ∀ x, x < 1 → F x = G x := by
  set D := fun x => F x - G x with hD
  have hDc : ContinuousOn D (Iio 1) := hF.sub hG
  have hdD : ∀ x, x ≠ 0 → x < 1 → HasDerivAt D 0 x := fun x h0 h1 => by
    exact ((hdF x h0 h1).sub (hdG x h0 h1)).congr_deriv (by simp)
  intro x hx
  rcases lt_trichotomy x 0 with h | h | h
  · obtain ⟨c, hc, hc'⟩ := exists_hasDerivAt_eq_slope D (fun _ => 0) h
      (hDc.mono fun t ht => by simp only [mem_Iio]; linarith [ht.2])
      (fun t ht => hdD t (ne_of_lt ht.2) (by linarith [ht.2]))
    have : D 0 - D x = 0 := by
      have hne : (0:ℝ) - x ≠ 0 := by linarith
      field_simp at hc'; linarith
    simp only [hD] at this; linarith
  · subst h; exact h0
  · obtain ⟨c, hc, hc'⟩ := exists_hasDerivAt_eq_slope D (fun _ => 0) h
      (hDc.mono fun t ht => by simp only [mem_Iio]; linarith [ht.2])
      (fun t ht => hdD t (ne_of_gt ht.1) (by linarith [ht.2]))
    have : D x - D 0 = 0 := by
      have hne : x - 0 ≠ 0 := by linarith
      field_simp at hc'; linarith
    simp only [hD] at this; linarith

lemma mem1 {a b : ℤ} {u v w : List ℤ} (hw : w ∈ shuffle u (b :: v)) :
    a :: w ∈ shuffle (a :: u) (b :: v) := by
  simp only [shuffle, List.mem_append, List.mem_map]; exact Or.inl ⟨w, hw, rfl⟩

lemma mem2 {a b : ℤ} {u v w : List ℤ} (hw : w ∈ shuffle (a :: u) v) :
    b :: w ∈ shuffle (a :: u) (b :: v) := by
  simp only [shuffle, List.mem_append, List.mem_map]; exact Or.inr ⟨w, hw, rfl⟩

lemma shuffle_thm : ∀ n : ℕ, ∀ u v : List ℤ, u.length + v.length = n → IsWord01 u →
    IsWord01 v → S03.Regular u → S03.Regular v →
    ∀ x : ℝ, x < 1 → S03.H u x * S03.H v x = ((shuffle u v).map (fun w => S03.H w x)).sum := by
  intro n
  induction n using Nat.strong_induction_on with
  | _ n ih =>
  intro u v hn hu hv hru hrv
  match u, v with
  | [], v => intro x _; simp [shuffle, S03.H]
  | a :: u, [] => intro x _; simp [shuffle, S03.H]
  | a :: u, b :: v =>
    have hu' := IsWord01_tail hu
    have hv' := IsWord01_tail hv
    have hru' := S03.reg_tail a u hru
    have hrv' := S03.reg_tail b v hrv
    have ih1 := ih (u.length + (b :: v).length) (by simp at hn ⊢; omega) u (b :: v) rfl hu' hv hru' hrv
    have ih2 := ih ((a :: u).length + v.length) (by simp at hn ⊢; omega) (a :: u) v rfl hu hv' hru hrv'
    have cu := (HG (a :: u) hu hru).1
    have cv := (HG (b :: v) hv hrv).1
    refine eq_of_deriv _ _ (fun x => S03.f a x * S03.H u x * S03.H (b :: v) x +
      S03.H (a :: u) x * (S03.f b x * S03.H v x)) (cu.mul cv) ?_ ?_ ?_ ?_
    · simp only [shuffle, List.map_append, List.sum_append, List.map_map]
      refine (contOn_lsum _ _ _ fun w hw => ?_).add (contOn_lsum _ _ _ fun w hw => ?_)
      · exact (HG (a :: w) (shuffle_w01 _ _ (by simpa using hu) hv _
          (mem1 hw))
          (shuffle_reg _ _ hru hrv _ (mem1 hw))).1
      · exact (HG (b :: w) (shuffle_w01 _ _ hu (by simpa using hv) _
          (mem2 hw))
          (shuffle_reg _ _ hru hrv _ (mem2 hw))).1
    · intro x hx0 hx1
      exact ((H_deriv' a u hu hru x hx0 hx1).mul (H_deriv' b v hv hrv x hx0 hx1)).congr_deriv
        (by ring)
    · intro x hx0 hx1
      simp only [shuffle, List.map_append, List.sum_append, List.map_map]
      have e1 := ih1 x hx1
      have e2 := ih2 x hx1
      have d1 := hasDerivAt_lsum (shuffle u (b :: v)) (fun w y => S03.H (a :: w) y)
        (fun w => S03.f a x * S03.H w x) x (fun w hw => H_deriv' a w
          (shuffle_w01 _ _ (by simpa using hu) hv _ (mem1 hw))
          (shuffle_reg _ _ hru hrv _ (mem1 hw)) x hx0 hx1)
      have d2 := hasDerivAt_lsum (shuffle (a :: u) v) (fun w y => S03.H (b :: w) y)
        (fun w => S03.f b x * S03.H w x) x (fun w hw => H_deriv' b w
          (shuffle_w01 _ _ hu (by simpa using hv) _ (mem2 hw))
          (shuffle_reg _ _ hru hrv _ (mem2 hw)) x hx0 hx1)
      refine (d1.add d2).congr_deriv ?_
      rw [List.sum_map_mul_left, List.sum_map_mul_left, ← e1, ← e2]; ring
    · rw [H_zero _ (by simp)]
      simp only [zero_mul, shuffle, List.map_append, List.sum_append, List.map_map]
      simp [Function.comp_def, H_zero]

end P03.H5

theorem P03.proof_7 : P03.stmt_7 := fun u v hu hv hru hrv x hx =>
  P03.H5.shuffle_thm _ u v rfl hu hv hru hrv x hx

namespace P03.H5
open MvPolynomial

/-- The weight of a generator. -/
def gw (i : Fin 13) : ℕ := (gen i).length

/-- `H_w` is a weighted-homogeneous rational polynomial of weight `|w|` in the 13 generators. -/
def Good (w : List ℤ) : Prop :=
  IsWord01 w ∧ S03.Regular w ∧ ∃ P : MvPolynomial (Fin 13) ℚ, IsWeightedHomogeneous gw P w.length ∧
    ∀ x : ℝ, x < 1 → S03.H w x = MvPolynomial.eval₂ (algebraMap ℚ ℝ) (fun i => S03.H (gen i) x) P

lemma good_nil : Good [] :=
  ⟨by simp [IsWord01], by simp [S03.Regular], 1, isWeightedHomogeneous_one ℚ gw, fun x _ => by
    simp [S03.H]⟩

lemma good_gen (i : Fin 13) : Good (gen i) :=
  ⟨by fin_cases i <;> simp [IsWord01, gen], by fin_cases i <;> simp [S03.Regular, gen],
    X i, isWeightedHomogeneous_X ℚ gw i, fun x _ => by simp⟩

lemma list_poly (n : ℕ) : ∀ L : List (List ℤ), (∀ w ∈ L, Good w ∧ w.length = n) →
    ∃ P : MvPolynomial (Fin 13) ℚ, IsWeightedHomogeneous gw P n ∧ ∀ x : ℝ, x < 1 →
      (L.map fun w => S03.H w x).sum =
        MvPolynomial.eval₂ (algebraMap ℚ ℝ) (fun i => S03.H (gen i) x) P
  | [], _ => ⟨0, isWeightedHomogeneous_zero ℚ gw n, fun x _ => by simp⟩
  | w :: L, h => by
    obtain ⟨⟨_, _, Pw, hPw, ew⟩, hl⟩ := h w (by simp)
    obtain ⟨PL, hPL, eL⟩ := list_poly n L (fun w' hw' => h w' (by simp [hw']))
    refine ⟨Pw + PL, (hl ▸ hPw).add hPL, fun x hx => ?_⟩
    simp only [List.map_cons, List.sum_cons, eval₂_add, ew x hx, eL x hx]

lemma step (u v w : List ℤ) (L : List (List ℤ)) (k : ℕ) (hk : k ≠ 0) (hu : Good u) (hv : Good v)
    (hp : (shuffle u v).Perm (List.replicate k w ++ L)) (hL : ∀ w' ∈ L, Good w')
    (hw : IsWord01 w) (hr : S03.Regular w) (hlen : w.length = u.length + v.length) : Good w := by
  obtain ⟨hu1, hu2, Pu, hPu, eu⟩ := hu
  obtain ⟨hv1, hv2, Pv, hPv, ev⟩ := hv
  have hLlen : ∀ w' ∈ L, Good w' ∧ w'.length = w.length := fun w' hw' =>
    ⟨hL w' hw', by rw [hlen]; exact shuffle_len u v w' (hp.symm.subset (by simp [hw']))⟩
  obtain ⟨PL, hPL, eL⟩ := list_poly w.length L hLlen
  refine ⟨hw, hr, C (1 / (k : ℚ)) * (Pu * Pv - PL), ?_, fun x hx => ?_⟩
  · refine IsWeightedHomogeneous.C_mul ?_ _
    exact (hlen ▸ hPu.mul hPv).sub hPL
  · have hs := shuffle_thm _ u v rfl hu1 hv1 hu2 hv2 x hx
    rw [(hp.map (fun w => S03.H w x)).sum_eq] at hs
    simp only [List.map_append, List.sum_append, List.map_replicate, List.sum_replicate,
      nsmul_eq_mul] at hs
    rw [eu x hx, ev x hx, eL x hx] at hs
    simp only [eval₂_mul, eval₂_sub, eval₂_C]
    have hk' : (k : ℝ) ≠ 0 := by exact_mod_cast hk
    rw [map_div₀, map_one, map_natCast]
    field_simp
    linarith

lemma g1 : Good [1] := good_gen 0

lemma g01 : Good [0, 1] := good_gen 1

lemma g001 : Good [0, 0, 1] := good_gen 2

lemma g011 : Good [0, 1, 1] := good_gen 3

lemma g0001 : Good [0, 0, 0, 1] := good_gen 4

lemma g0011 : Good [0, 0, 1, 1] := good_gen 5

lemma g0111 : Good [0, 1, 1, 1] := good_gen 6

lemma g00001 : Good [0, 0, 0, 0, 1] := good_gen 7

lemma g00011 : Good [0, 0, 0, 1, 1] := good_gen 8

lemma g00101 : Good [0, 0, 1, 0, 1] := good_gen 9

lemma g00111 : Good [0, 0, 1, 1, 1] := good_gen 10

lemma g01011 : Good [0, 1, 0, 1, 1] := good_gen 11

lemma g01111 : Good [0, 1, 1, 1, 1] := good_gen 12

lemma g11 : Good [1, 1] :=
  step [1] [1] [1, 1] [] 2 (by norm_num) g1 g1
    (by rw [show shuffle [1] [1] = [[1, 1], [1, 1]] by simp [shuffle]]; decide) (by simp) (by simp [IsWord01]) (by simp [S03.Regular]) rfl

lemma g101 : Good [1, 0, 1] :=
  step [1] [0, 1] [1, 0, 1] [[0, 1, 1], [0, 1, 1]] 1 (by norm_num) g1 g01
    (by rw [show shuffle [1] [0, 1] = [[1, 0, 1], [0, 1, 1], [0, 1, 1]] by simp [shuffle]]; decide) (by intro w' hw'; simp only [List.mem_cons, List.not_mem_nil, or_false] at hw'; rcases hw' with rfl | rfl <;> first | exact g011) (by simp [IsWord01]) (by simp [S03.Regular]) rfl

lemma g111 : Good [1, 1, 1] :=
  step [1] [1, 1] [1, 1, 1] [] 3 (by norm_num) g1 g11
    (by rw [show shuffle [1] [1, 1] = [[1, 1, 1], [1, 1, 1], [1, 1, 1]] by simp [shuffle]]; decide) (by simp) (by simp [IsWord01]) (by simp [S03.Regular]) rfl

lemma g0101 : Good [0, 1, 0, 1] :=
  step [0, 1] [0, 1] [0, 1, 0, 1] [[0, 0, 1, 1], [0, 0, 1, 1], [0, 0, 1, 1], [0, 0, 1, 1]] 2 (by norm_num) g01 g01
    (by rw [show shuffle [0, 1] [0, 1] = [[0, 1, 0, 1], [0, 0, 1, 1], [0, 0, 1, 1], [0, 0, 1, 1], [0, 0, 1, 1], [0, 1, 0, 1]] by simp [shuffle]]; decide) (by intro w' hw'; simp only [List.mem_cons, List.not_mem_nil, or_false] at hw'; rcases hw' with rfl | rfl | rfl | rfl <;> first | exact g0011) (by simp [IsWord01]) (by simp [S03.Regular]) rfl

lemma g1001 : Good [1, 0, 0, 1] :=
  step [1] [0, 0, 1] [1, 0, 0, 1] [[0, 1, 0, 1], [0, 0, 1, 1], [0, 0, 1, 1]] 1 (by norm_num) g1 g001
    (by rw [show shuffle [1] [0, 0, 1] = [[1, 0, 0, 1], [0, 1, 0, 1], [0, 0, 1, 1], [0, 0, 1, 1]] by simp [shuffle]]; decide) (by intro w' hw'; simp only [List.mem_cons, List.not_mem_nil, or_false] at hw'; rcases hw' with rfl | rfl | rfl <;> first | exact g0011 | exact g0101) (by simp [IsWord01]) (by simp [S03.Regular]) rfl

lemma g1011 : Good [1, 0, 1, 1] :=
  step [1] [0, 1, 1] [1, 0, 1, 1] [[0, 1, 1, 1], [0, 1, 1, 1], [0, 1, 1, 1]] 1 (by norm_num) g1 g011
    (by rw [show shuffle [1] [0, 1, 1] = [[1, 0, 1, 1], [0, 1, 1, 1], [0, 1, 1, 1], [0, 1, 1, 1]] by simp [shuffle]]; decide) (by intro w' hw'; simp only [List.mem_cons, List.not_mem_nil, or_false] at hw'; rcases hw' with rfl | rfl | rfl <;> first | exact g0111) (by simp [IsWord01]) (by simp [S03.Regular]) rfl

lemma g1101 : Good [1, 1, 0, 1] :=
  step [1] [1, 0, 1] [1, 1, 0, 1] [[1, 0, 1, 1], [1, 0, 1, 1]] 2 (by norm_num) g1 g101
    (by rw [show shuffle [1] [1, 0, 1] = [[1, 1, 0, 1], [1, 1, 0, 1], [1, 0, 1, 1], [1, 0, 1, 1]] by simp [shuffle]]; decide) (by intro w' hw'; simp only [List.mem_cons, List.not_mem_nil, or_false] at hw'; rcases hw' with rfl | rfl <;> first | exact g1011) (by simp [IsWord01]) (by simp [S03.Regular]) rfl

lemma g1111 : Good [1, 1, 1, 1] :=
  step [1] [1, 1, 1] [1, 1, 1, 1] [] 4 (by norm_num) g1 g111
    (by rw [show shuffle [1] [1, 1, 1] = [[1, 1, 1, 1], [1, 1, 1, 1], [1, 1, 1, 1], [1, 1, 1, 1]] by simp [shuffle]]; decide) (by simp) (by simp [IsWord01]) (by simp [S03.Regular]) rfl

lemma g01001 : Good [0, 1, 0, 0, 1] :=
  step [0, 1] [0, 0, 1] [0, 1, 0, 0, 1] [[0, 0, 1, 0, 1], [0, 0, 0, 1, 1], [0, 0, 0, 1, 1], [0, 0, 1, 0, 1], [0, 0, 0, 1, 1], [0, 0, 0, 1, 1], [0, 0, 0, 1, 1], [0, 0, 0, 1, 1], [0, 0, 1, 0, 1]] 1 (by norm_num) g01 g001
    (by rw [show shuffle [0, 1] [0, 0, 1] = [[0, 1, 0, 0, 1], [0, 0, 1, 0, 1], [0, 0, 0, 1, 1], [0, 0, 0, 1, 1], [0, 0, 1, 0, 1], [0, 0, 0, 1, 1], [0, 0, 0, 1, 1], [0, 0, 0, 1, 1], [0, 0, 0, 1, 1], [0, 0, 1, 0, 1]] by simp [shuffle]]; decide) (by intro w' hw'; simp only [List.mem_cons, List.not_mem_nil, or_false] at hw'; rcases hw' with rfl | rfl | rfl | rfl | rfl | rfl | rfl | rfl | rfl <;> first | exact g00011 | exact g00101) (by simp [IsWord01]) (by simp [S03.Regular]) rfl

lemma g01101 : Good [0, 1, 1, 0, 1] :=
  step [0, 1] [0, 1, 1] [0, 1, 1, 0, 1] [[0, 1, 0, 1, 1], [0, 0, 1, 1, 1], [0, 0, 1, 1, 1], [0, 0, 1, 1, 1], [0, 0, 1, 1, 1], [0, 0, 1, 1, 1], [0, 0, 1, 1, 1], [0, 1, 0, 1, 1], [0, 1, 0, 1, 1]] 1 (by norm_num) g01 g011
    (by rw [show shuffle [0, 1] [0, 1, 1] = [[0, 1, 0, 1, 1], [0, 0, 1, 1, 1], [0, 0, 1, 1, 1], [0, 0, 1, 1, 1], [0, 0, 1, 1, 1], [0, 0, 1, 1, 1], [0, 0, 1, 1, 1], [0, 1, 0, 1, 1], [0, 1, 0, 1, 1], [0, 1, 1, 0, 1]] by simp [shuffle]]; decide) (by intro w' hw'; simp only [List.mem_cons, List.not_mem_nil, or_false] at hw'; rcases hw' with rfl | rfl | rfl | rfl | rfl | rfl | rfl | rfl | rfl <;> first | exact g00111 | exact g01011) (by simp [IsWord01]) (by simp [S03.Regular]) rfl

lemma g10001 : Good [1, 0, 0, 0, 1] :=
  step [1] [0, 0, 0, 1] [1, 0, 0, 0, 1] [[0, 1, 0, 0, 1], [0, 0, 1, 0, 1], [0, 0, 0, 1, 1], [0, 0, 0, 1, 1]] 1 (by norm_num) g1 g0001
    (by rw [show shuffle [1] [0, 0, 0, 1] = [[1, 0, 0, 0, 1], [0, 1, 0, 0, 1], [0, 0, 1, 0, 1], [0, 0, 0, 1, 1], [0, 0, 0, 1, 1]] by simp [shuffle]]; decide) (by intro w' hw'; simp only [List.mem_cons, List.not_mem_nil, or_false] at hw'; rcases hw' with rfl | rfl | rfl | rfl <;> first | exact g00011 | exact g00101 | exact g01001) (by simp [IsWord01]) (by simp [S03.Regular]) rfl

lemma g10011 : Good [1, 0, 0, 1, 1] :=
  step [1] [0, 0, 1, 1] [1, 0, 0, 1, 1] [[0, 1, 0, 1, 1], [0, 0, 1, 1, 1], [0, 0, 1, 1, 1], [0, 0, 1, 1, 1]] 1 (by norm_num) g1 g0011
    (by rw [show shuffle [1] [0, 0, 1, 1] = [[1, 0, 0, 1, 1], [0, 1, 0, 1, 1], [0, 0, 1, 1, 1], [0, 0, 1, 1, 1], [0, 0, 1, 1, 1]] by simp [shuffle]]; decide) (by intro w' hw'; simp only [List.mem_cons, List.not_mem_nil, or_false] at hw'; rcases hw' with rfl | rfl | rfl | rfl <;> first | exact g00111 | exact g01011) (by simp [IsWord01]) (by simp [S03.Regular]) rfl

lemma g10101 : Good [1, 0, 1, 0, 1] :=
  step [1] [0, 1, 0, 1] [1, 0, 1, 0, 1] [[0, 1, 1, 0, 1], [0, 1, 1, 0, 1], [0, 1, 0, 1, 1], [0, 1, 0, 1, 1]] 1 (by norm_num) g1 g0101
    (by rw [show shuffle [1] [0, 1, 0, 1] = [[1, 0, 1, 0, 1], [0, 1, 1, 0, 1], [0, 1, 1, 0, 1], [0, 1, 0, 1, 1], [0, 1, 0, 1, 1]] by simp [shuffle]]; decide) (by intro w' hw'; simp only [List.mem_cons, List.not_mem_nil, or_false] at hw'; rcases hw' with rfl | rfl | rfl | rfl <;> first | exact g01011 | exact g01101) (by simp [IsWord01]) (by simp [S03.Regular]) rfl

lemma g10111 : Good [1, 0, 1, 1, 1] :=
  step [1] [0, 1, 1, 1] [1, 0, 1, 1, 1] [[0, 1, 1, 1, 1], [0, 1, 1, 1, 1], [0, 1, 1, 1, 1], [0, 1, 1, 1, 1]] 1 (by norm_num) g1 g0111
    (by rw [show shuffle [1] [0, 1, 1, 1] = [[1, 0, 1, 1, 1], [0, 1, 1, 1, 1], [0, 1, 1, 1, 1], [0, 1, 1, 1, 1], [0, 1, 1, 1, 1]] by simp [shuffle]]; decide) (by intro w' hw'; simp only [List.mem_cons, List.not_mem_nil, or_false] at hw'; rcases hw' with rfl | rfl | rfl | rfl <;> first | exact g01111) (by simp [IsWord01]) (by simp [S03.Regular]) rfl

lemma g11001 : Good [1, 1, 0, 0, 1] :=
  step [1] [1, 0, 0, 1] [1, 1, 0, 0, 1] [[1, 0, 1, 0, 1], [1, 0, 0, 1, 1], [1, 0, 0, 1, 1]] 2 (by norm_num) g1 g1001
    (by rw [show shuffle [1] [1, 0, 0, 1] = [[1, 1, 0, 0, 1], [1, 1, 0, 0, 1], [1, 0, 1, 0, 1], [1, 0, 0, 1, 1], [1, 0, 0, 1, 1]] by simp [shuffle]]; decide) (by intro w' hw'; simp only [List.mem_cons, List.not_mem_nil, or_false] at hw'; rcases hw' with rfl | rfl | rfl <;> first | exact g10011 | exact g10101) (by simp [IsWord01]) (by simp [S03.Regular]) rfl

lemma g11011 : Good [1, 1, 0, 1, 1] :=
  step [1] [1, 0, 1, 1] [1, 1, 0, 1, 1] [[1, 0, 1, 1, 1], [1, 0, 1, 1, 1], [1, 0, 1, 1, 1]] 2 (by norm_num) g1 g1011
    (by rw [show shuffle [1] [1, 0, 1, 1] = [[1, 1, 0, 1, 1], [1, 1, 0, 1, 1], [1, 0, 1, 1, 1], [1, 0, 1, 1, 1], [1, 0, 1, 1, 1]] by simp [shuffle]]; decide) (by intro w' hw'; simp only [List.mem_cons, List.not_mem_nil, or_false] at hw'; rcases hw' with rfl | rfl | rfl <;> first | exact g10111) (by simp [IsWord01]) (by simp [S03.Regular]) rfl

lemma g11101 : Good [1, 1, 1, 0, 1] :=
  step [1] [1, 1, 0, 1] [1, 1, 1, 0, 1] [[1, 1, 0, 1, 1], [1, 1, 0, 1, 1]] 3 (by norm_num) g1 g1101
    (by rw [show shuffle [1] [1, 1, 0, 1] = [[1, 1, 1, 0, 1], [1, 1, 1, 0, 1], [1, 1, 1, 0, 1], [1, 1, 0, 1, 1], [1, 1, 0, 1, 1]] by simp [shuffle]]; decide) (by intro w' hw'; simp only [List.mem_cons, List.not_mem_nil, or_false] at hw'; rcases hw' with rfl | rfl <;> first | exact g11011) (by simp [IsWord01]) (by simp [S03.Regular]) rfl

lemma g11111 : Good [1, 1, 1, 1, 1] :=
  step [1] [1, 1, 1, 1] [1, 1, 1, 1, 1] [] 5 (by norm_num) g1 g1111
    (by rw [show shuffle [1] [1, 1, 1, 1] = [[1, 1, 1, 1, 1], [1, 1, 1, 1, 1], [1, 1, 1, 1, 1], [1, 1, 1, 1, 1], [1, 1, 1, 1, 1]] by simp [shuffle]]; decide) (by simp) (by simp [IsWord01]) (by simp [S03.Regular]) rfl

lemma good_all : ∀ w ∈ ([[], [1], [0, 1], [1, 1], [0, 0, 1], [0, 1, 1], [1, 0, 1], [1, 1, 1], [0, 0, 0, 1], [0, 0, 1, 1], [0, 1, 0, 1], [0, 1, 1, 1], [1, 0, 0, 1], [1, 0, 1, 1], [1, 1, 0, 1], [1, 1, 1, 1], [0, 0, 0, 0, 1], [0, 0, 0, 1, 1], [0, 0, 1, 0, 1], [0, 0, 1, 1, 1], [0, 1, 0, 0, 1], [0, 1, 0, 1, 1], [0, 1, 1, 0, 1], [0, 1, 1, 1, 1], [1, 0, 0, 0, 1], [1, 0, 0, 1, 1], [1, 0, 1, 0, 1], [1, 0, 1, 1, 1], [1, 1, 0, 0, 1], [1, 1, 0, 1, 1], [1, 1, 1, 0, 1], [1, 1, 1, 1, 1]] : List (List ℤ)), Good w := by
  simp only [List.forall_mem_cons, List.not_mem_nil, IsEmpty.forall_iff, implies_true, and_true]
  exact ⟨good_nil, g1, g01, g11, g001, g011, g101, g111, g0001, g0011, g0101, g0111, g1001, g1011, g1101, g1111, g00001, g00011, g00101, g00111, g01001, g01011, g01101, g01111, g10001, g10011, g10101, g10111, g11001, g11011, g11101, g11111⟩

end P03.H5

theorem P03.proof_9 : P03.stmt_9 := by
  intro w hw hr hl
  have hmem : w ∈ ([[], [1], [0, 1], [1, 1], [0, 0, 1], [0, 1, 1], [1, 0, 1], [1, 1, 1], [0, 0, 0, 1], [0, 0, 1, 1], [0, 1, 0, 1], [0, 1, 1, 1], [1, 0, 0, 1], [1, 0, 1, 1], [1, 1, 0, 1], [1, 1, 1, 1], [0, 0, 0, 0, 1], [0, 0, 0, 1, 1], [0, 0, 1, 0, 1], [0, 0, 1, 1, 1], [0, 1, 0, 0, 1], [0, 1, 0, 1, 1], [0, 1, 1, 0, 1], [0, 1, 1, 1, 1], [1, 0, 0, 0, 1], [1, 0, 0, 1, 1], [1, 0, 1, 0, 1], [1, 0, 1, 1, 1], [1, 1, 0, 0, 1], [1, 1, 0, 1, 1], [1, 1, 1, 0, 1], [1, 1, 1, 1, 1]] : List (List ℤ)) := by
    have l : ∀ a ∈ w, a = 0 ∨ a = 1 := hw
    rcases w with _ | ⟨a, _ | ⟨b, _ | ⟨c, _ | ⟨d, _ | ⟨e, _ | ⟨f, w⟩⟩⟩⟩⟩⟩
    · decide
    all_goals try (simp at hl; omega)
    all_goals
      rcases l a (by simp) with rfl | rfl
    all_goals first
      | decide
      | (exfalso; simp [S03.Regular] at hr; done)
      | (rcases l b (by simp) with rfl | rfl)
    all_goals first
      | decide
      | (exfalso; simp [S03.Regular] at hr; done)
      | (rcases l c (by simp) with rfl | rfl)
    all_goals first
      | decide
      | (exfalso; simp [S03.Regular] at hr; done)
      | (rcases l d (by simp) with rfl | rfl)
    all_goals first
      | decide
      | (exfalso; simp [S03.Regular] at hr; done)
      | (rcases l e (by simp) with rfl | rfl)
    all_goals first
      | decide
      | (exfalso; simp [S03.Regular] at hr)
  obtain ⟨_, _, P, hP, e⟩ := P03.H5.good_all w hmem
  refine ⟨P, fun m hm => ?_, e⟩
  have := hP (MvPolynomial.mem_support_iff.mp hm)
  simpa [Finsupp.weight_apply, P03.H5.gw, smul_eq_mul] using this

namespace P03.H5

lemma wordOf_cons (s1 : ℕ) (s : List ℕ) :
    wordOf (s1 :: s) = List.replicate (s1 - 1) 0 ++ 1 :: wordOf s := by
  simp [wordOf]

lemma dd_one_cons (u : List ℤ) (hu : u ≠ []) :
    S03.dd (1 :: u) = (1, 1, (S03.dd u).2.2.1,
      ((S03.dd u).1 * ((S03.dd u).2.1 : ℤ)) :: (S03.dd u).2.2.2) := by
  cases u with
  | nil => exact absurd rfl hu
  | cons b w => simp [S03.dd]

lemma dd_zero_cons (u : List ℤ) (hu : u ≠ []) :
    S03.dd (0 :: u) = ((S03.dd u).1, (S03.dd u).2.1 + 1, (S03.dd u).2.2.1, (S03.dd u).2.2.2) := by
  cases u with
  | nil => exact absurd rfl hu
  | cons b w => simp [S03.dd]

lemma dd_rep (u : List ℤ) (hu : u ≠ []) : ∀ m : ℕ, S03.dd (List.replicate m 0 ++ u) =
    ((S03.dd u).1, (S03.dd u).2.1 + m, (S03.dd u).2.2.1, (S03.dd u).2.2.2)
  | 0 => by simp
  | m + 1 => by
    rw [List.replicate_succ, List.cons_append, dd_zero_cons _ (by simp [hu]), dd_rep u hu m]
    simp only [Prod.mk.injEq, true_and, and_true]; omega

lemma dd_word : ∀ s : List ℕ, (∀ si ∈ s, 1 ≤ si) → ∀ s1 : ℕ, 1 ≤ s1 →
    S03.dd (wordOf (s1 :: s)) = (1, s1, 1, s.map (fun si : ℕ => (si : ℤ)))
  | [], _, s1, h1 => by
    rw [wordOf_cons, dd_rep _ (by simp)]
    simp [wordOf, S03.dd]; omega
  | s2 :: s, hs, s1, h1 => by
    have ih := dd_word s (fun si hsi => hs si (by simp [hsi])) s2 (hs s2 (by simp))
    rw [wordOf_cons, dd_rep _ (by simp), dd_one_cons _ (by rw [wordOf_cons]; simp), ih]
    simp; omega

lemma wordOf_letters (s : List ℕ) : S03.IsLetters (wordOf s) := by
  intro a ha
  simp only [wordOf, List.mem_flatMap, List.mem_append, List.mem_replicate, List.mem_singleton] at ha
  obtain ⟨_, _, h | h⟩ := ha
  · left; exact h.2
  · right; left; exact h

lemma wordOf_last : ∀ s : List ℕ, s ≠ [] → (wordOf s).getLast? = some 1
  | [], h => absurd rfl h
  | [s1], _ => by simp [wordOf]
  | s1 :: s2 :: s, _ => by
    rw [wordOf_cons, List.getLast?_append_of_ne_nil _ (by rw [wordOf_cons]; simp)]
    rw [List.getLast?_cons, wordOf_last (s2 :: s) (by simp)]
    simp

end P03.H5

theorem P03.proof_6 : P03.stmt_6 := by
  intro s1 s h1 hs x hx
  have hne : wordOf (s1 :: s) ≠ [] := by rw [P03.H5.wordOf_cons]; simp
  have hreg : S03.Regular (wordOf (s1 :: s)) := by
    intro a ha
    rw [P03.H5.wordOf_last _ (by simp)] at ha
    simp at ha; rw [← ha]; norm_num
  obtain ⟨_, _, _, _, _, hS⟩ := (S03.proof_5).2.1 _ hne (P03.H5.wordOf_letters _) hreg
  have := hS x hx
  rw [P03.H5.dd_word s hs s1 h1] at this
  convert this using 1
  funext n
  simp only [Int.cast_one, one_pow, one_mul]
  ring

namespace P03.H5

/-- Cast of a positive index list to an index word. -/
def ci (u : List ℕ) : List ℤ := u.map (fun a : ℕ => (a : ℤ))

def Pos (u : List ℕ) : Prop := ∀ a ∈ u, 1 ≤ a

lemma gq_pos (a : ℕ) (ha : 1 ≤ a) (i : ℕ) : S03.gq (a : ℤ) i = 1 / (i : ℚ) ^ a := by
  unfold S03.gq
  rw [Int.sign_eq_one_of_pos (by omega)]
  simp

lemma pstuffle : ∀ u v : List ℕ, Pos u → Pos v →
    ∃ L : List (List ℕ), (∀ w ∈ L, Pos w ∧ w.sum = u.sum + v.sum) ∧
      ∀ M : ℕ, S03.Z (ci u) M * S03.Z (ci v) M = (L.map fun w => S03.Z (ci w) M).sum := by
  intro u
  induction u with
  | nil =>
    intro v _ hv
    exact ⟨[v], by simp [hv], fun M => by simp [S03.Z, ci]⟩
  | cons a u ihu =>
    intro v hu hv
    have ha : 1 ≤ a := hu a (by simp)
    have hu' : Pos u := fun x hx => hu x (by simp [hx])
    induction v with
    | nil => exact ⟨[a :: u], by simp [hu], fun M => by simp [S03.Z, ci]⟩
    | cons b v ihv =>
      have hb : 1 ≤ b := hv b (by simp)
      have hv' : Pos v := fun x hx => hv x (by simp [hx])
      obtain ⟨L1, hL1, hZ1⟩ := ihu (b :: v) hu' hv
      obtain ⟨L2, hL2, hZ2⟩ := ihv hv'
      obtain ⟨L3, hL3, hZ3⟩ := ihu v hu' hv'
      refine ⟨L1.map (a :: ·) ++ L2.map (b :: ·) ++ L3.map ((a + b) :: ·), ?_, ?_⟩
      · intro w hw
        simp only [List.mem_append, List.mem_map] at hw
        rcases hw with (⟨w', hw', rfl⟩ | ⟨w', hw', rfl⟩) | ⟨w', hw', rfl⟩
        · obtain ⟨h1, h2⟩ := hL1 w' hw'
          refine ⟨fun x hx => ?_, ?_⟩
          · simp at hx; rcases hx with rfl | hx; exact ha; exact h1 x hx
          · simp only [List.sum_cons] at h2 ⊢; omega
        · obtain ⟨h1, h2⟩ := hL2 w' hw'
          refine ⟨fun x hx => ?_, ?_⟩
          · simp at hx; rcases hx with rfl | hx; exact hb; exact h1 x hx
          · simp only [List.sum_cons] at h2 ⊢; omega
        · obtain ⟨h1, h2⟩ := hL3 w' hw'
          refine ⟨fun x hx => ?_, ?_⟩
          · simp at hx; rcases hx with rfl | hx; omega; exact h1 x hx
          · simp only [List.sum_cons] at h2 ⊢; omega
      · intro M
        rw [List.map_append, List.map_append, List.sum_append, List.sum_append,
          List.map_map, List.map_map, List.map_map]
        have e : ∀ (c : ℕ) (w : List ℕ), ci (c :: w) = (c : ℤ) :: ci w := fun c w => rfl
        rw [show ((fun w => S03.Z (ci w) M) ∘ (a :: ·)) = fun w => S03.Z ((a : ℤ) :: ci w) M from rfl,
          show ((fun w => S03.Z (ci w) M) ∘ (b :: ·)) = fun w => S03.Z ((b : ℤ) :: ci w) M from rfl,
          show ((fun w => S03.Z (ci w) M) ∘ ((a + b) :: ·)) =
            fun w => S03.Z (((a + b : ℕ) : ℤ) :: ci w) M from rfl]
        have s1 := S03.sum_map_cons (a : ℤ) (L1.map ci) M
        have s2 := S03.sum_map_cons (b : ℤ) (L2.map ci) M
        have s3 := S03.sum_map_cons ((a + b : ℕ) : ℤ) (L3.map ci) M
        simp only [List.map_map, Function.comp_def] at s1 s2 s3
        rw [s1, s2, s3, e, e, S03.stuffle_core]
        simp_rw [← hZ1, ← hZ2, ← hZ3, gq_pos a ha, gq_pos b hb, gq_pos (a + b) (by omega)]
        congr 1
        refine Finset.sum_congr rfl fun i _ => ?_
        rw [pow_add]; field_simp

/-- `F(n)` is a rational combination of `(n+1)^{-k} Z_v(n)`, `v` positive, `k + |v| ≤ d`. -/
def Rep (d : ℕ) (F : ℕ → ℚ) : Prop :=
  ∃ T : List (ℚ × ℕ × List ℕ), (∀ τ ∈ T, Pos τ.2.2 ∧ τ.2.1 + τ.2.2.sum ≤ d) ∧
    ∀ n : ℕ, F n = (T.map fun τ => τ.1 / ((n : ℚ) + 1) ^ τ.2.1 * S03.Z (ci τ.2.2) n).sum

lemma rep_zero (d : ℕ) : Rep d (fun _ => 0) := ⟨[], by simp, fun n => by simp⟩

lemma rep_mono {d d' : ℕ} {F : ℕ → ℚ} (h : Rep d F) (hd : d ≤ d') : Rep d' F := by
  obtain ⟨T, hT, e⟩ := h
  exact ⟨T, fun τ hτ => ⟨(hT τ hτ).1, (hT τ hτ).2.trans hd⟩, e⟩

lemma rep_add {d : ℕ} {F G : ℕ → ℚ} (hF : Rep d F) (hG : Rep d G) : Rep d (fun n => F n + G n) := by
  obtain ⟨T, hT, e⟩ := hF
  obtain ⟨T', hT', e'⟩ := hG
  refine ⟨T ++ T', fun τ hτ => ?_, fun n => by beta_reduce; rw [e n, e' n, List.map_append, List.sum_append]⟩
  rcases List.mem_append.mp hτ with h | h
  · exact hT τ h
  · exact hT' τ h

lemma rep_tt (c c' : ℚ) (k k' : ℕ) (u v : List ℕ) (hu : Pos u) (hv : Pos v) :
    Rep (k + u.sum + (k' + v.sum)) (fun n => c / ((n : ℚ) + 1) ^ k * S03.Z (ci u) n *
      (c' / ((n : ℚ) + 1) ^ k' * S03.Z (ci v) n)) := by
  obtain ⟨L, hL, hZ⟩ := pstuffle u v hu hv
  refine ⟨L.map fun w => (c * c', k + k', w), fun τ hτ => ?_, fun n => ?_⟩
  · simp only [List.mem_map] at hτ
    obtain ⟨w, hw, rfl⟩ := hτ
    refine ⟨(hL w hw).1, ?_⟩
    simp only [(hL w hw).2]; omega
  · simp only [List.map_map, Function.comp_def]
    rw [List.sum_map_mul_left, ← hZ n, pow_add]
    field_simp

lemma rep_tG (c : ℚ) (k : ℕ) (u : List ℕ) (hu : Pos u) (d : ℕ) :
    ∀ T : List (ℚ × ℕ × List ℕ), (∀ τ ∈ T, Pos τ.2.2 ∧ τ.2.1 + τ.2.2.sum ≤ d) →
    Rep (k + u.sum + d) (fun n => c / ((n : ℚ) + 1) ^ k * S03.Z (ci u) n *
      (T.map fun τ => τ.1 / ((n : ℚ) + 1) ^ τ.2.1 * S03.Z (ci τ.2.2) n).sum)
  | [], _ => by simpa using rep_zero _
  | σ :: T, h => by
    have h1 := rep_mono (rep_tt c σ.1 k σ.2.1 u σ.2.2 hu (h σ (by simp)).1)
      (show k + u.sum + (σ.2.1 + σ.2.2.sum) ≤ k + u.sum + d by have := (h σ (by simp)).2; omega)
    have h2 := rep_tG c k u hu d T (fun τ hτ => h τ (by simp [hτ]))
    have := rep_add h1 h2
    obtain ⟨T', hT', e'⟩ := this
    exact ⟨T', hT', fun n => by rw [← e' n]; beta_reduce; simp only [List.map_cons, List.sum_cons]; ring⟩

lemma rep_mul {d1 d2 : ℕ} {F G : ℕ → ℚ} (hF : Rep d1 F) (hG : Rep d2 G) :
    Rep (d1 + d2) (fun n => F n * G n) := by
  obtain ⟨T, hT, e⟩ := hF
  obtain ⟨T', hT', e'⟩ := hG
  have key : ∀ T : List (ℚ × ℕ × List ℕ), (∀ τ ∈ T, Pos τ.2.2 ∧ τ.2.1 + τ.2.2.sum ≤ d1) →
      Rep (d1 + d2) (fun n => (T.map fun τ => τ.1 / ((n : ℚ) + 1) ^ τ.2.1 *
        S03.Z (ci τ.2.2) n).sum * G n) := by
    intro T
    induction T with
    | nil => intro _; simpa using rep_zero _
    | cons τ T ih =>
      intro h
      have h1 := rep_mono (rep_tG τ.1 τ.2.1 τ.2.2 (h τ (by simp)).1 d2 T' hT')
        (show τ.2.1 + τ.2.2.sum + d2 ≤ d1 + d2 by have := (h τ (by simp)).2; omega)
      have h2 := ih (fun σ hσ => h σ (by simp [hσ]))
      obtain ⟨T'', hT'', e''⟩ := rep_add h1 h2
      exact ⟨T'', hT'', fun n => by rw [← e'' n]; beta_reduce; rw [e' n]; simp only [List.map_cons, List.sum_cons]; ring⟩
  obtain ⟨T'', hT'', e''⟩ := key T hT
  exact ⟨T'', hT'', fun n => by rw [← e'' n]; beta_reduce; rw [e n]⟩

lemma rep_hn (j : ℕ) (hj : 1 ≤ j) : Rep j (fun n => S03.hn j (n + 1)) := by
  refine ⟨[(1, 0, [j]), (1, j, [])], fun τ hτ => ?_, fun n => ?_⟩
  · simp only [List.mem_cons, List.not_mem_nil, or_false] at hτ
    rcases hτ with rfl | rfl
    · exact ⟨fun a ha => by simp at ha; omega, by simp⟩
    · exact ⟨fun a ha => by simp at ha, by simp⟩
  · beta_reduce; rw [S03.hn_succ]
    simp only [List.map_cons, List.map_nil, List.sum_cons, List.sum_nil, ci, S03.Z,
      pow_zero, div_one, one_mul, add_zero, mul_one]
    congr 1
    · unfold S03.hn
      refine Finset.sum_congr rfl fun i _ => ?_
      rw [Int.sign_eq_one_of_pos (by omega)]; simp

lemma rep_pow (j : ℕ) (hj : 1 ≤ j) : ∀ e : ℕ, Rep (j * e) (fun n => S03.hn j (n + 1) ^ e)
  | 0 => ⟨[(1, 0, [])], by simp [Pos], fun n => by simp [S03.Z, ci]⟩
  | e + 1 => by
    have := rep_mul (rep_pow j hj e) (rep_hn j hj)
    obtain ⟨T, hT, h⟩ := this
    exact ⟨T, fun τ hτ => ⟨(hT τ hτ).1, by have := (hT τ hτ).2; rw [Nat.mul_succ]; exact this⟩,
      fun n => by rw [← h n]; beta_reduce; rw [pow_succ]⟩

lemma rep_prod (e : ℕ → ℕ) (s : Finset ℕ) (hs : ∀ j ∈ s, 1 ≤ j) :
    Rep (∑ j ∈ s, j * e j) (fun n => ∏ j ∈ s, S03.hn j (n + 1) ^ e j) := by
  induction s using Finset.induction_on with
  | empty => exact ⟨[(1, 0, [])], by simp [Pos], fun n => by simp [S03.Z, ci]⟩
  | insert j s hj ih =>
    have h1 := rep_pow j (hs j (by simp)) (e j)
    have h2 := ih (fun i hi => hs i (by simp [hi]))
    obtain ⟨T, hT, h⟩ := rep_mul h1 h2
    refine ⟨T, fun τ hτ => ⟨(hT τ hτ).1, ?_⟩, fun n => ?_⟩
    · rw [Finset.sum_insert hj]; exact (hT τ hτ).2
    · rw [← h n]; beta_reduce; rw [Finset.prod_insert hj]

lemma hasSum_terms (a : ℕ) (ha : 1 ≤ a) (x : ℝ) (hx : |x| < 1) :
    ∀ T : List (ℚ × ℕ × List ℕ), (∀ τ ∈ T, Pos τ.2.2) →
    HasSum (fun n : ℕ => x ^ (n + 1) / ((n : ℝ) + 1) ^ a *
      (((T.map fun τ => τ.1 / ((n : ℚ) + 1) ^ τ.2.1 * S03.Z (ci τ.2.2) n).sum : ℚ) : ℝ))
      ((T.map fun τ => (τ.1 : ℝ) * S03.H (wordOf ((a + τ.2.1) :: τ.2.2)) x).sum)
  | [], _ => by simp
  | τ :: T, h => by
    have h1 := (P03.proof_6 (a + τ.2.1) τ.2.2 (by omega) (h τ (by simp)) x hx).mul_left (τ.1 : ℝ)
    have h2 := hasSum_terms a ha x hx T (fun σ hσ => h σ (by simp [hσ]))
    convert h1.add h2 using 1
    · funext n
      simp only [List.map_cons, List.sum_cons, Rat.cast_add, mul_add]
      congr 1
      simp only [ci, Rat.cast_mul, Rat.cast_div, Rat.cast_pow, Rat.cast_add, Rat.cast_natCast,
        Rat.cast_one]
      rw [pow_add]; ring
    · simp

end P03.H5

theorem P03.proof_8 : P03.stmt_8 := by
  intro a ha J e
  obtain ⟨T, hT, hF⟩ := P03.H5.rep_prod e (Finset.Icc 1 J) (fun j hj => (Finset.mem_Icc.mp hj).1)
  refine ⟨T.map fun τ => (τ.1, (a + τ.2.1) :: τ.2.2), fun τ hτ => ?_, fun x hx => ?_⟩
  · simp only [List.mem_map] at hτ
    obtain ⟨σ, hσ, rfl⟩ := hτ
    refine ⟨by simp, fun si hsi => ?_, ?_⟩
    · simp at hsi; rcases hsi with rfl | hsi; omega; exact (hT σ hσ).1 si hsi
    · simp only [List.sum_cons]; have := (hT σ hσ).2; omega
  · have := P03.H5.hasSum_terms a ha x hx T (fun τ hτ => (hT τ hτ).1)
    rw [List.map_map]
    convert this using 1
    funext n
    have hn := hF n
    beta_reduce at hn
    rw [← hn]
    push_cast
    rfl
    rfl

namespace P03.H5
open Real MeasureTheory Set
open scoped Interval

/-- `e_j(1, 1/2, …, 1/n)`. -/
noncomputable def E (j n : ℕ) : ℝ :=
  ((Finset.Icc 1 n).val.map (fun i : ℕ => (1 : ℝ) / (i : ℝ))).esymm j

lemma esymm_cons (a : ℝ) (s : Multiset ℝ) (j : ℕ) :
    (a ::ₘ s).esymm (j + 1) = s.esymm (j + 1) + a * s.esymm j := by
  simp only [Multiset.esymm, Multiset.powersetCard_cons, Multiset.map_add, Multiset.sum_add,
    Multiset.map_map, Function.comp_def, Multiset.prod_cons, Multiset.sum_map_mul_left]

lemma E_zero_left (n : ℕ) : E 0 n = 1 := by simp [E]

lemma E_zero_right (j : ℕ) : E (j + 1) 0 = 0 := by
  unfold E; apply Multiset.esymm_of_card_lt; simp

lemma E_succ (j n : ℕ) : E (j + 1) (n + 1) = E (j + 1) n + E j n / ((n : ℝ) + 1) := by
  unfold E
  rw [← Finset.insert_Icc_right_eq_Icc_add_one (by omega),
    Finset.insert_val_of_notMem (by simp), Multiset.map_cons, esymm_cons]
  push_cast; ring

lemma E_bound : ∀ n j : ℕ, 0 ≤ E j n ∧ E j n ≤ ((n : ℝ) + 1) ^ j
  | n, 0 => by simp [E_zero_left]
  | 0, j + 1 => by simp [E_zero_right]
  | n + 1, j + 1 => by
    obtain ⟨a0, a1⟩ := E_bound n (j + 1)
    obtain ⟨b0, b1⟩ := E_bound n j
    rw [E_succ]
    have hn : (0:ℝ) < (n:ℝ) + 1 := by positivity
    have hq : E j n / ((n:ℝ) + 1) ≤ ((n:ℝ) + 1) ^ j := by
      rw [div_le_iff₀ hn]
      have : (1:ℝ) ≤ (n:ℝ) + 1 := by linarith [(Nat.cast_nonneg n : (0:ℝ) ≤ n)]
      nlinarith [pow_nonneg hn.le j]
    have hp : ((n:ℝ) + 1) ^ j ≤ ((n:ℝ) + 1 + 1) ^ j := pow_le_pow_left₀ hn.le (by linarith) j
    refine ⟨by positivity, ?_⟩
    push_cast
    calc E (j + 1) n + E j n / ((n:ℝ) + 1) ≤ ((n:ℝ) + 1) ^ (j + 1) + ((n:ℝ) + 1) ^ j := by
          linarith
      _ ≤ ((n:ℝ) + 1 + 1) ^ (j + 1) := by
          rw [pow_succ, pow_succ]; nlinarith [pow_nonneg hn.le j]

lemma lamc_eq (c n : ℕ) : lamc c (n + 1) = (-1 : ℝ) ^ c * (c.factorial : ℝ) / ((n : ℝ) + 1) *
    E (c - 1) n := by
  simp [lamc, E]

lemma lamc_bound (c n : ℕ) : |lamc c (n + 1)| ≤ (c.factorial : ℝ) * ((n : ℝ) + 1) ^ c := by
  rw [lamc_eq]
  obtain ⟨h0, h1⟩ := E_bound n (c - 1)
  have hn : (0:ℝ) < (n:ℝ) + 1 := by positivity
  have h1' : (1:ℝ) ≤ (n:ℝ) + 1 := by linarith [(Nat.cast_nonneg n : (0:ℝ) ≤ n)]
  rw [abs_mul, abs_div, abs_mul, abs_pow, abs_neg, abs_one, one_pow, one_mul, abs_of_pos hn,
    abs_of_nonneg h0, Nat.abs_cast]
  have hf : (0:ℝ) ≤ c.factorial := Nat.cast_nonneg _
  have : E (c - 1) n ≤ ((n:ℝ) + 1) ^ c := h1.trans (pow_le_pow_right₀ h1' (by omega))
  calc (c.factorial : ℝ) / ((n:ℝ) + 1) * E (c - 1) n ≤ (c.factorial : ℝ) / 1 * E (c - 1) n := by
        gcongr
    _ ≤ (c.factorial : ℝ) * ((n:ℝ) + 1) ^ c := by rw [div_one]; gcongr

lemma log_series (c : ℕ) (hc : 1 ≤ c) (z : ℝ) (hz : |z| < 1) :
    HasSum (fun q : ℕ => lamc c (q + 1) * z ^ (q + 1)) (Real.log (1 - z) ^ c) := by
  induction c, hc using Nat.le_induction generalizing z with
  | base =>
    have := (Real.hasSum_pow_div_log_of_abs_lt_one hz).neg
    simp only [neg_neg] at this
    rw [pow_one]
    convert this using 1
    funext q
    rw [lamc_eq]; simp [E_zero_left]; ring
  | succ c hc ih =>
    -- scaled coefficients
    set K : ℝ := ((c + 1).factorial : ℝ) with hK
    have hK0 : 0 < K := by positivity
    set cc : ℕ → ℝ := fun n => lamc (c + 1) (n + 1) / K with hcc
    have hcc_eq : ∀ n, cc n * ((n : ℝ) + 1) = (-1 : ℝ) ^ (c + 1) * E c n := by
      intro n; simp only [hcc, lamc_eq, Nat.add_sub_cancel]; rw [hK]
      field_simp
    have hbd : ∀ n, |cc n| ≤ ((n : ℝ) + 1) ^ c := by
      intro n
      have hn : (0:ℝ) < (n:ℝ) + 1 := by positivity
      have h1' : (1:ℝ) ≤ (n:ℝ) + 1 := by linarith [(Nat.cast_nonneg n : (0:ℝ) ≤ n)]
      have e1 : cc n = (-1 : ℝ) ^ (c + 1) * E c n / ((n:ℝ) + 1) := by
        rw [← hcc_eq n]; field_simp
      obtain ⟨h0, h1⟩ := E_bound n c
      rw [e1, abs_div, abs_mul, abs_pow, abs_neg, abs_one, one_pow, one_mul, abs_of_pos hn,
        abs_of_nonneg h0, div_le_iff₀ hn]
      nlinarith [pow_nonneg hn.le c]
    -- the derivative series
    set g : ℝ → ℝ := fun t => -((c : ℝ) + 1) * Real.log (1 - t) ^ c / (1 - t) / K with hg
    have hsum : ∀ t : ℝ, |t| < 1 → ∑' n, cc n * ((n : ℝ) + 1) * t ^ n = g t := by
      intro t ht
      have hA := ((S03.psD cc c hbd t ht).2.1).hasSum
      set A := ∑' n, cc n * ((n : ℝ) + 1) * t ^ n
      have h1 := (hasSum_nat_add_iff' 1).mpr hA
      have h2 := hA.mul_right t
      have h3 := h1.sub h2
      have h4 := (ih t ht).mul_left (-((c : ℝ) + 1) / K)
      have hEq : (fun n => cc (n + 1) * (((n + 1 : ℕ) : ℝ) + 1) * t ^ (n + 1) -
          cc n * ((n : ℝ) + 1) * t ^ n * t) =
          fun n => -((c : ℝ) + 1) / K * (lamc c (n + 1) * t ^ (n + 1)) := by
        funext n
        obtain ⟨d, rfl⟩ : ∃ d, c = d + 1 := ⟨c - 1, by omega⟩
        rw [hcc_eq, hcc_eq, mul_assoc _ (t ^ n), ← pow_succ, lamc_eq, E_succ d n]
        simp only [Nat.add_sub_cancel, hK, Nat.factorial_succ]
        push_cast
        field_simp
        ring
      rw [hEq] at h3
      have hv := h3.unique h4
      simp only [Finset.range_one, Finset.sum_singleton, CharP.cast_eq_zero, zero_add, pow_zero,
        mul_one] at hv
      have h0 : cc 0 = 0 := by
        simp only [hcc, lamc_eq]
        obtain ⟨d, rfl⟩ : ∃ d, c = d + 1 := ⟨c - 1, by omega⟩
        simp [E_zero_right]
      rw [h0] at hv
      have ht1 : 1 - t ≠ 0 := by
        have := (abs_lt.mp ht).2; linarith
      have h5 : A * (1 - t) = -((c : ℝ) + 1) / K * Real.log (1 - t) ^ c := by linarith
      simp only [hg]
      rw [eq_div_iff hK0.ne', eq_div_iff ht1, mul_right_comm, h5]
      field_simp
    -- integrate
    have hz1 : ∀ t ∈ uIcc 0 z, |t| < 1 := by
      intro t ht
      rw [mem_uIcc] at ht
      rw [abs_lt] at hz ⊢
      rcases ht with h | h <;> constructor <;> linarith [h.1, h.2]
    have hF : ∀ t ∈ uIcc 0 z, HasDerivAt (fun y => Real.log (1 - y) ^ (c + 1) / K) (g t) t := by
      intro t ht
      have ht1 : 0 < 1 - t := by have := (abs_lt.mp (hz1 t ht)).2; linarith
      have := ((((hasDerivAt_id t).const_sub 1).log ht1.ne').pow (c + 1)).div_const K
      convert this using 1
      · funext y; simp
      · simp only [hg, id, Nat.add_sub_cancel]; push_cast; field_simp
    have hgc : ContinuousOn g (uIcc 0 z) := by
      intro t ht
      have ht1 : 0 < 1 - t := by have := (abs_lt.mp (hz1 t ht)).2; linarith
      apply ContinuousAt.continuousWithinAt
      simp only [hg]
      have : ContinuousAt (fun y : ℝ => Real.log (1 - y)) t :=
        (continuousAt_const.sub continuousAt_id).log ht1.ne'
      exact (((continuousAt_const.mul (this.pow c)).div (continuousAt_const.sub continuousAt_id)
        ht1.ne').div_const K)
    have hint : IntervalIntegrable g volume 0 z := hgc.intervalIntegrable
    have key := S03.ftc_series cc c hbd g z hz hint (fun t ht _ => hsum t (by linarith))
    rw [intervalIntegral.integral_eq_sub_of_hasDerivAt hF hint] at key
    simp only [sub_zero, Real.log_one, ne_eq, Nat.add_eq_zero_iff, one_ne_zero, and_false,
      not_false_eq_true, zero_pow, zero_div] at key
    have := key.mul_left K
    convert this using 1
    · funext n; simp only [hcc]; field_simp
    · field_simp

end P03.H5

theorem P03.proof_3 : P03.stmt_3 := fun c hc z hz => P03.H5.log_series c hc z hz

namespace P03.H5
open Real MeasureTheory Set Filter Topology
open scoped Interval

lemma re_lower {z : ℂ} {δ : ℝ} (hδ1 : δ ≤ 1) (hz : z.re ≤ 1 - δ) {t : ℝ} (ht0 : 0 ≤ t)
    (ht1 : t ≤ 1) : δ ≤ (1 - z * (t : ℂ)).re := by
  simp only [Complex.sub_re, Complex.one_re, Complex.mul_re, Complex.ofReal_re, Complex.ofReal_im,
    mul_zero, sub_zero]
  rcases le_total z.re 0 with h | h
  · nlinarith
  · nlinarith

lemma norm_clog_le {w : ℂ} {δ M : ℝ} (hδ : 0 < δ) (hw : δ ≤ w.re) (hM : ‖w‖ ≤ M) :
    ‖Complex.log w‖ ≤ M + 1 / δ + π := by
  have hwn : δ ≤ ‖w‖ := hw.trans (Complex.re_le_norm w)
  have hwp : 0 < ‖w‖ := lt_of_lt_of_le hδ hwn
  refine (Complex.norm_le_abs_re_add_abs_im _).trans ?_
  rw [Complex.log_re, Complex.log_im]
  have h1 := S05.abs_log_le_sum hwp
  have h2 := Complex.abs_arg_le_pi w
  have h3 : 1 / ‖w‖ ≤ 1 / δ := one_div_le_one_div_of_le hδ hwn
  linarith

lemma J_analytic (n a b c : ℕ) : AnalyticOnNhd ℝ (S05.Jint n a b c) (Iio 1) := by
  set g : ℝ → ℝ := fun t => t ^ n * Real.log t ^ a * Real.log (1 - t) ^ b with hg
  have gi : IntervalIntegrable g volume 0 1 := S05.ii_g n a b
  have gm : Measurable g := by rw [hg]; fun_prop
  set F : ℂ → ℝ → ℂ := fun z t => (g t : ℂ) * Complex.log (1 - z * t) ^ c with hF
  set F' : ℂ → ℝ → ℂ := fun z t => (g t : ℂ) *
    ((c : ℂ) * Complex.log (1 - z * t) ^ (c - 1) * (-(t : ℂ) / (1 - z * t))) with hF'
  set Jc : ℂ → ℂ := fun z => ∫ t in (0:ℝ)..1, F z t with hJc
  have Fm : ∀ z, Measurable (F z) := fun z => by rw [hF]; fun_prop
  have F'm : ∀ z, Measurable (F' z) := fun z => by rw [hF']; fun_prop
  have hd : DifferentiableOn ℂ Jc {z : ℂ | z.re < 1} := by
    intro z0 hz0
    have hz0 : z0.re < 1 := hz0
    set ε : ℝ := (1 - z0.re) / 2 with hε
    have hε0 : 0 < ε := by rw [hε]; linarith
    set δ : ℝ := min ε 1 with hδ
    have hδ0 : 0 < δ := lt_min hε0 one_pos
    have hδ1 : δ ≤ 1 := min_le_right _ _
    set M : ℝ := 1 + ‖z0‖ + ε with hM
    set Λ : ℝ := M + 1 / δ + π with hΛ
    have hΛ0 : 0 ≤ Λ := by rw [hΛ, hM]; positivity
    have key : ∀ z ∈ Metric.ball z0 ε, ∀ t ∈ Icc (0:ℝ) 1,
        δ ≤ (1 - z * (t : ℂ)).re ∧ ‖1 - z * (t : ℂ)‖ ≤ M := by
      intro z hz t ht
      have hzz : ‖z - z0‖ < ε := by simpa [dist_eq_norm] using hz
      have hre : z.re ≤ 1 - δ := by
        have := Complex.abs_re_le_norm (z - z0)
        rw [Complex.sub_re] at this
        have := (abs_le.mp this).2
        have : δ ≤ ε := min_le_left _ _
        linarith
      refine ⟨re_lower hδ1 hre ht.1 ht.2, ?_⟩
      calc ‖1 - z * (t : ℂ)‖ ≤ ‖(1 : ℂ)‖ + ‖z‖ * ‖(t : ℂ)‖ := by
            refine (norm_sub_le _ _).trans ?_; rw [norm_mul]
        _ ≤ 1 + (‖z0‖ + ε) * 1 := by
            have : ‖z‖ ≤ ‖z0‖ + ε := by
              have := norm_le_norm_add_norm_sub' z z0
              linarith
            have ht' : ‖(t : ℂ)‖ ≤ 1 := by
              rw [Complex.norm_real, Real.norm_eq_abs, abs_of_nonneg ht.1]; exact ht.2
            gcongr
            · simp
        _ = M := by rw [hM]; ring
    have hslit : ∀ z ∈ Metric.ball z0 ε, ∀ t ∈ Icc (0:ℝ) 1, (1 - z * (t : ℂ)) ∈ Complex.slitPlane :=
      fun z hz t ht => Complex.mem_slitPlane_iff.mpr (Or.inl (lt_of_lt_of_le hδ0 (key z hz t ht).1))
    have hIoc : ∀ t ∈ Ι (0:ℝ) 1, t ∈ Icc (0:ℝ) 1 := fun t ht => by
      rw [uIoc_of_le zero_le_one] at ht; exact ⟨ht.1.le, ht.2⟩
    have hres := intervalIntegral.hasDerivAt_integral_of_dominated_loc_of_deriv_le
      (μ := volume) (a := 0) (b := 1) (F := F) (F' := F') (x₀ := z0)
      (bound := fun t => ‖g t‖ * ((c : ℝ) * Λ ^ (c - 1) * (1 / δ)))
      (Metric.ball_mem_nhds z0 hε0)
      (Eventually.of_forall fun z => (Fm z).aestronglyMeasurable)
      (by
        refine IntervalIntegrable.mono_fun' (g := fun t => ‖g t‖ * Λ ^ c)
          (gi.norm.mul_const _) (Fm z0).aestronglyMeasurable ?_
        refine (ae_restrict_iff' measurableSet_uIoc).mpr (Eventually.of_forall fun t ht => ?_)
        obtain ⟨k1, k2⟩ := key z0 (Metric.mem_ball_self hε0) t (hIoc t ht)
        have hl := norm_clog_le hδ0 k1 k2
        simp only [hF, norm_mul, norm_pow, Complex.norm_real]
        gcongr)
      (F'm z0).aestronglyMeasurable
      (Eventually.of_forall fun t ht z hz => by
        obtain ⟨k1, k2⟩ := key z hz t (hIoc t ht)
        have hl := norm_clog_le hδ0 k1 k2
        have hwn : δ ≤ ‖1 - z * (t : ℂ)‖ := k1.trans (Complex.re_le_norm _)
        have hwp : 0 < ‖1 - z * (t : ℂ)‖ := lt_of_lt_of_le hδ0 hwn
        have ht1 : ‖(t : ℂ)‖ ≤ 1 := by
          rw [Complex.norm_real, Real.norm_eq_abs, abs_of_nonneg (hIoc t ht).1]; exact (hIoc t ht).2
        simp only [hF', norm_mul, norm_pow, Complex.norm_real, norm_div, norm_neg, Complex.norm_natCast]
        gcongr
        rw [Real.norm_eq_abs, abs_of_nonneg (hIoc t ht).1]; exact (hIoc t ht).2)
      (gi.norm.mul_const _)
      (Eventually.of_forall fun t ht z hz => by
        have h1 : HasDerivAt (fun z : ℂ => 1 - z * (t : ℂ)) (-(t : ℂ)) z := by
          simpa using ((hasDerivAt_id z).mul_const (t : ℂ)).const_sub 1
        have h2 := (h1.clog (hslit z hz t (hIoc t ht))).pow c
        exact h2.const_mul (g t : ℂ))
    exact hres.2.differentiableAt.differentiableWithinAt
  have hU : IsOpen {z : ℂ | z.re < 1} := isOpen_lt Complex.continuous_re continuous_const
  have han : AnalyticOnNhd ℂ Jc {z : ℂ | z.re < 1} := hd.analyticOnNhd hU
  intro x0 hx0
  have hx0 : x0 < 1 := hx0
  have h1 : AnalyticAt ℂ Jc (x0 : ℂ) := han _ (by simpa using hx0)
  have h2 : AnalyticAt ℝ Jc (x0 : ℂ) := h1.restrictScalars
  have h3 : AnalyticAt ℝ (fun x : ℝ => (Jc (x : ℂ)).re) x0 :=
    Complex.reCLM.analyticAt _ |>.comp (h2.comp (Complex.ofRealCLM.analyticAt x0))
  refine h3.congr ?_
  filter_upwards [isOpen_Iio.mem_nhds hx0] with x hx
  have hx : x < 1 := hx
  have e : Jc (x : ℂ) = ((S05.Jint n a b c x : ℝ) : ℂ) := by
    rw [hJc, S05.Jint, ← intervalIntegral.integral_ofReal]
    refine intervalIntegral.integral_congr fun t ht => ?_
    rw [uIcc_of_le zero_le_one] at ht
    have hpos : 0 < 1 - x * t := by
      rcases le_total x 0 with h | h
      · nlinarith [ht.1]
      · nlinarith [ht.2]
    simp only [hF, hg]
    rw [show (1 : ℂ) - (x : ℂ) * (t : ℂ) = ((1 - x * t : ℝ) : ℂ) by push_cast; ring,
      ← Complex.ofReal_log hpos.le]
    push_cast; ring
  simp only [e, Complex.ofReal_re]

end P03.H5

namespace P03.H5
open Real MeasureTheory Set Filter Topology
open scoped Interval

lemma abs_log_lin {x t : ℝ} (hx : x < 1) (ht0 : 0 ≤ t) (ht1 : t ≤ 1) :
    |Real.log (1 - x * t)| ≤ |Real.log (1 - x)| := by
  rcases le_total x 0 with h | h
  · have h1 : 1 ≤ 1 - x * t := by nlinarith
    have h2 : 1 - x * t ≤ 1 - x := by nlinarith
    rw [abs_of_nonneg (Real.log_nonneg h1), abs_of_nonneg (Real.log_nonneg (by linarith))]
    exact Real.log_le_log (by linarith) h2
  · have h1 : 1 - x * t ≤ 1 := by nlinarith
    have h2 : 1 - x ≤ 1 - x * t := by nlinarith
    have hp : 0 < 1 - x := by linarith
    rw [abs_of_nonpos (Real.log_nonpos (by linarith) h1),
      abs_of_nonpos (Real.log_nonpos hp.le (by linarith))]
    linarith [Real.log_le_log hp h2]

lemma ii_J' (x : ℝ) (hx : x < 1) (n a b c : ℕ) :
    IntervalIntegrable
        (fun t : ℝ => t ^ n * Real.log t ^ a * Real.log (1 - t) ^ b * Real.log (1 - x * t) ^ c)
        volume 0 1 := by
  refine S05.ii_of_ls (by fun_prop : Measurable (fun t : ℝ =>
      t ^ n * Real.log t ^ a * Real.log (1 - t) ^ b * Real.log (1 - x * t) ^ c)).aestronglyMeasurable
    (|Real.log (1 - x)| ^ c) (a + b) (fun t ht => ?_)
  rw [abs_mul, abs_pow]
  have h1 := S05.g_bound (n := n) (a := a) (b := b) ht
  have h2 := pow_le_pow_left₀ (abs_nonneg _) (abs_log_lin hx ht.1.le ht.2.le) c
  have := pow_nonneg (S05.ls_nonneg t) (a + b)
  calc |t ^ n * Real.log t ^ a * Real.log (1 - t) ^ b| * |Real.log (1 - x * t)| ^ c
      ≤ (1 * S05.ls t ^ (a + b)) * |Real.log (1 - x)| ^ c := by gcongr
    _ = _ := by ring

lemma J_bigO (n a b c : ℕ) : (S05.Jint n a b c) =O[𝓝 0] (fun x : ℝ => x ^ c) := by
  set g : ℝ → ℝ := fun t => t ^ n * Real.log t ^ a * Real.log (1 - t) ^ b with hg
  have gi : IntervalIntegrable g volume 0 1 := S05.ii_g n a b
  set C : ℝ := (∫ t in (0:ℝ)..1, ‖g t‖) * 2 ^ c
  refine Asymptotics.IsBigO.of_bound C ?_
  have hball : Ioo (-(1/2 : ℝ)) (1/2) ∈ 𝓝 (0 : ℝ) := Ioo_mem_nhds (by norm_num) (by norm_num)
  filter_upwards [hball] with x hx
  have hxa : |x| ≤ 1 / 2 := abs_le.mpr ⟨hx.1.le, hx.2.le⟩
  have hbd : ∀ t ∈ Ι (0:ℝ) 1, ‖g t * Real.log (1 - x * t) ^ c‖ ≤ ‖g t‖ * (2 * |x|) ^ c := by
    intro t ht
    rw [uIoc_of_le zero_le_one] at ht
    have hxt : |x * t| ≤ |x| := by
      rw [abs_mul, abs_of_pos ht.1]; exact mul_le_of_le_one_right (abs_nonneg _) ht.2
    have hxt1 : |x * t| < 1 := by linarith
    have hl := Real.abs_log_sub_add_sum_range_le hxt1 0
    simp only [Finset.range_zero, Finset.sum_empty, zero_add, pow_one] at hl
    have hl2 : |Real.log (1 - x * t)| ≤ 2 * |x| := by
      have h1 : 1 / 2 ≤ 1 - |x * t| := by linarith
      calc |Real.log (1 - x * t)| ≤ |x * t| / (1 - |x * t|) := hl
        _ ≤ |x| / (1 / 2) := by
            gcongr
        _ = 2 * |x| := by ring
    rw [norm_mul, norm_pow, Real.norm_eq_abs (Real.log _)]
    gcongr
  have := intervalIntegral.norm_integral_le_of_norm_le zero_le_one
    (Eventually.of_forall fun t ht => hbd t ht) (gi.norm.mul_const _)
  rw [intervalIntegral.integral_mul_const] at this
  calc ‖S05.Jint n a b c x‖ = ‖∫ t in (0:ℝ)..1, g t * Real.log (1 - x * t) ^ c‖ := rfl
    _ ≤ (∫ t in (0:ℝ)..1, ‖g t‖) * (2 * |x|) ^ c := this
    _ = C * ‖x ^ c‖ := by
        rw [mul_pow, norm_pow, Real.norm_eq_abs]; ring

lemma J_series (x : ℝ) (hx : |x| < 1) (n a b c : ℕ) (hc1 : 1 ≤ c) :
    HasSum (fun q : ℕ => x ^ (q + 1) * lamc c (q + 1) * S05.Iint (n + q + 1) a b)
      (S05.Jint n a b c x) := by
  set g := fun t : ℝ => t ^ n * Real.log t ^ a * Real.log (1 - t) ^ b with hg
  have hgi : IntervalIntegrable g MeasureTheory.volume 0 1 := S05.ii_g n a b
  have hS : Summable (fun q : ℕ => |lamc c (q + 1)| * |x| ^ (q + 1)) := by
    have hx' : ‖|x|‖ < 1 := by rwa [Real.norm_eq_abs, abs_abs]
    have h0 : Summable (fun n : ℕ => (n : ℝ) ^ c * |x| ^ n) :=
      summable_pow_mul_geometric_of_norm_lt_one c hx'
    have h1 : Summable (fun n : ℕ => ((n + 1 : ℕ) : ℝ) ^ c * |x| ^ (n + 1)) :=
      (summable_nat_add_iff (f := fun n : ℕ => (n : ℝ) ^ c * |x| ^ n) 1).mpr h0
    have hs := h1.mul_left (c.factorial : ℝ)
    apply Summable.of_nonneg_of_le (fun q => by positivity) (fun q => ?_) hs
    have := lamc_bound c q
    push_cast
    calc |lamc c (q + 1)| * |x| ^ (q + 1) ≤ ((c.factorial : ℝ) * ((q : ℝ) + 1) ^ c) * |x| ^ (q + 1) :=
          mul_le_mul_of_nonneg_right this (by positivity)
      _ = _ := by ring
  have key := intervalIntegral.hasSum_integral_of_dominated_convergence
    (μ := MeasureTheory.volume) (a := 0) (b := 1)
    (F := fun (q : ℕ) (t : ℝ) => g t * (lamc c (q + 1) * (x * t) ^ (q + 1)))
    (f := fun t => g t * Real.log (1 - x * t) ^ c)
    (fun q t => |lamc c (q + 1)| * |x| ^ (q + 1) * |g t|)
    (fun q => (by fun_prop : Measurable (fun t : ℝ =>
      g t * (lamc c (q + 1) * (x * t) ^ (q + 1)))).aestronglyMeasurable)
    (fun q => by
      refine Filter.Eventually.of_forall (fun t ht => ?_)
      rw [Set.uIoc_of_le zero_le_one] at ht
      have hgt : ∀ G L : ℝ, ‖G * (L * (x * t) ^ (q + 1))‖ =
          |L| * |x| ^ (q + 1) * |G| * t ^ (q + 1) := by
        intro G L
        rw [Real.norm_eq_abs, abs_mul, abs_mul, abs_pow, abs_mul, abs_of_pos ht.1]; ring
      rw [hgt]
      have : t ^ (q + 1) ≤ 1 := pow_le_one₀ ht.1.le ht.2
      have h0 : 0 ≤ |lamc c (q + 1)| * |x| ^ (q + 1) * |g t| := by positivity
      calc |lamc c (q + 1)| * |x| ^ (q + 1) * |g t| * t ^ (q + 1)
          ≤ |lamc c (q + 1)| * |x| ^ (q + 1) * |g t| * 1 := by gcongr
        _ = _ := by ring)
    (Filter.Eventually.of_forall (fun t _ => hS.mul_right _))
    (by
      simp_rw [tsum_mul_right]
      exact hgi.norm.const_mul _)
    (Filter.Eventually.of_forall (fun t ht => by
      rw [Set.uIoc_of_le zero_le_one] at ht
      have hxt : |x * t| < 1 := by
        rw [abs_mul, abs_of_pos ht.1]
        calc |x| * t ≤ |x| * 1 := by gcongr; exact ht.2
          _ < 1 := by linarith
      exact (log_series c hc1 (x * t) hxt).mul_left (g t)))
  have e1 : (∫ t in (0:ℝ)..1, g t * Real.log (1 - x * t) ^ c) = S05.Jint n a b c x := by
    unfold S05.Jint; rfl
  rw [e1] at key
  convert key using 1
  funext q
  unfold S05.Iint
  rw [← intervalIntegral.integral_const_mul]
  congr 1; funext t
  simp only [hg]
  rw [mul_pow, show n + q + 1 = n + (q + 1) by ring, pow_add]
  ring

end P03.H5

theorem P03.proof_4 : P03.stmt_4 := by
  intro n a b c hc1 _
  exact ⟨fun x hx => P03.H5.ii_J' x hx n a b c, P03.H5.J_analytic n a b c,
    P03.H5.J_bigO n a b c, fun x hx => P03.H5.J_series x hx n a b c hc1⟩

namespace P03.P5
open S05 Set Real Filter Topology MeasureTheory
open scoped Interval

lemma nd_exp5 {h : ℝ → ℝ} (hh : ContDiffOn ℝ 5 h (Ioi 0)) {p : ℝ} (hp : 0 < p) {n : ℕ}
    (hn : n ≤ 5) {D : ℝ} (hD : ∀ i, 1 ≤ i → i ≤ n → nd i h p ≤ D ^ i) :
    nd n (fun q => Real.exp (h q)) p ≤ (n.factorial : ℝ) * Real.exp (h p) * D ^ n := by
  have := norm_iteratedFDerivWithin_comp_le (g := Real.exp) (f := h) (t := univ) (s := Ioi 0)
    (x := p) (n := n) (N := 5) (contDiff_exp.contDiffOn) hh (by exact_mod_cast hn)
    uniqueDiffOn_univ uS (mapsTo_univ _ _) hp (C := Real.exp (h p)) (D := D) ?_
    (fun i h1 h2 => by rw [nd_F]; exact hD i h1 h2)
  · rw [← nd_F]; exact this
  · intro i _
    rw [iteratedFDerivWithin_univ, norm_iteratedFDeriv_eq_norm_iteratedDeriv,
      iteratedDeriv_eq_iterate, Real.iter_deriv_exp, Real.norm_eq_abs, abs_of_pos (exp_pos _)]

lemma log_derivs5 (i : ℕ) (hi : i ≤ 5) : |iteratedDeriv i Real.log 1| ≤ 24 := by
  rcases i with _ | i
  · simp
  · rw [iteratedDeriv_succ', Real.deriv_log', iteratedDeriv_eq_iterate, iter_deriv_inv]
    simp only [one_zpow, mul_one, abs_mul, abs_pow, abs_neg, abs_one, one_pow, one_mul]
    have : i ≤ 4 := by omega
    interval_cases i <;> simp [Nat.factorial] <;> norm_num

lemma nd_log5 {f : ℝ → ℝ} (hf : ContDiffOn ℝ 5 f (Ioi 0)) (hpos : ∀ q ∈ Ioi (0:ℝ), 0 < f q)
    {p : ℝ} (hp : 0 < p) {n : ℕ} (hn1 : 1 ≤ n) (hn : n ≤ 5) {D : ℝ}
    (hD : ∀ i, 1 ≤ i → i ≤ n → nd i f p ≤ f p * D ^ i) :
    nd n (fun q => Real.log (f q)) p ≤ (n.factorial : ℝ) * 24 * D ^ n := by
  have hc : 0 < f p := hpos p hp
  set φ : ℝ → ℝ := fun q => (f p)⁻¹ * f q with hφ
  have hEq : EqOn (fun q => Real.log (f q)) (fun q => Real.log (f p) + Real.log (φ q)) (Ioi 0) := by
    intro q hq
    simp only [hφ]
    rw [Real.log_mul (inv_ne_zero hc.ne') (hpos q hq).ne', Real.log_inv]; ring
  have h1 : nd n (fun q => Real.log (f q)) p = nd n (fun q => Real.log (φ q)) p := by
    unfold nd
    rw [iteratedDerivWithin_congr hEq hp, iteratedDerivWithin_const_add (by omega)]
  rw [h1]
  have hφc : ContDiffOn ℝ 5 φ (Ioi 0) := contDiffOn_const.mul hf
  have hφD : ∀ i, 1 ≤ i → i ≤ n → nd i φ p ≤ D ^ i := by
    intro i hi1 hi2
    have : nd i φ p = (f p)⁻¹ * nd i f p := by
      unfold nd
      rw [hφ, iteratedDerivWithin_const_mul_field, abs_mul, abs_of_pos (inv_pos.2 hc)]
    rw [this, inv_mul_le_iff₀ hc]; exact hD i hi1 hi2
  have := norm_iteratedFDerivWithin_comp_le (g := Real.log) (f := φ) (t := Ioi 0) (s := Ioi 0)
    (x := p) (n := n) (N := 5) (Real.contDiffOn_log.mono (fun x hx => ne_of_gt hx)) hφc
    (by exact_mod_cast hn) uS uS (fun q hq => by
      simp only [hφ, mem_Ioi]; exact mul_pos (inv_pos.2 hc) (hpos q hq)) hp (C := 24) (D := D) ?_
    (fun i h1 h2 => by rw [nd_F]; exact hφD i h1 h2)
  · rw [← nd_F]; simpa [Function.comp_def, mul_assoc] using this
  · intro i hi
    have h1 : φ p = 1 := by simp [hφ, inv_mul_cancel₀ hc.ne']
    rw [h1, iteratedFDerivWithin_of_isOpen i isOpen_Ioi (show (1:ℝ) ∈ Ioi 0 by norm_num),
      norm_iteratedFDeriv_eq_norm_iteratedDeriv, Real.norm_eq_abs]
    exact log_derivs5 i (hi.trans hn)

lemma nd_mul5 {f g : ℝ → ℝ} (hf : ContDiffOn ℝ 5 f (Ioi 0)) (hg : ContDiffOn ℝ 5 g (Ioi 0))
    {p : ℝ} (hp : 0 < p) {n : ℕ} (hn : n ≤ 5) {a b D : ℝ}
    (hfD : ∀ i ≤ n, nd i f p ≤ a * D ^ i) (hgD : ∀ i ≤ n, nd i g p ≤ b * D ^ i) :
    nd n (fun q => f q * g q) p ≤ 2 ^ n * (a * b) * D ^ n := by
  have := norm_iteratedFDerivWithin_mul_le (N := 5) hf hg uS hp (n := n) (by exact_mod_cast hn)
  rw [nd_F] at this
  simp only [nd_F] at this
  refine this.trans ?_
  calc ∑ i ∈ Finset.range (n + 1), (n.choose i : ℝ) * nd i f p * nd (n - i) g p
      ≤ ∑ i ∈ Finset.range (n + 1), (n.choose i : ℝ) * (a * b * D ^ n) := by
        refine Finset.sum_le_sum fun i hi => ?_
        have hi' : i ≤ n := Nat.lt_succ_iff.1 (Finset.mem_range.1 hi)
        have e : a * b * D ^ n = (a * D ^ i) * (b * D ^ (n - i)) := by
          rw [show D ^ n = D ^ i * D ^ (n - i) by rw [← pow_add, Nat.add_sub_cancel' hi'] ]; ring
        rw [mul_assoc, e]
        exact mul_le_mul_of_nonneg_left (mul_le_mul (hfD i hi') (hgD _ (Nat.sub_le _ _))
          (nd_nonneg _ _ _) ((nd_nonneg _ _ _).trans (hfD i hi'))) (Nat.cast_nonneg _)
    _ = 2 ^ n * (a * b) * D ^ n := by rw [← Finset.sum_mul, choose_sum]; ring

lemma nd_add5 {f g : ℝ → ℝ} (hf : ContDiffOn ℝ 5 f (Ioi 0)) (hg : ContDiffOn ℝ 5 g (Ioi 0))
    {p : ℝ} (hp : 0 < p) {n : ℕ} (hn : n ≤ 5) :
    nd n (fun q => f q + g q) p ≤ nd n f p + nd n g p := by
  unfold nd
  rw [iteratedDerivWithin_fun_add hp uS ((hf p hp).of_le (by exact_mod_cast hn))
    ((hg p hp).of_le (by exact_mod_cast hn))]
  exact abs_add_le _ _

lemma nd_inv_part5 (i : ℕ) (hi : i ≤ 5) {p : ℝ} (hp : 1 / 2 ≤ p) :
    nd i (fun q => -(1 + q⁻¹)) p ≤ 7680 := by
  have hp0 : 0 < p := by linarith
  have e : (fun q : ℝ => -(1 + q⁻¹)) = fun q => (-1) * (1 + q⁻¹) := by ext q; ring
  rw [e, nd_cmul]
  simp only [abs_neg, abs_one, one_mul]
  unfold nd
  rw [iteratedDerivWithin_of_isOpen isOpen_Ioi hp0]
  have hinv : p⁻¹ ≤ 2 := by rw [inv_le_comm₀ hp0 (by norm_num)]; linarith
  rcases Nat.eq_zero_or_pos i with rfl | hi0
  · simp only [iteratedDeriv_zero]
    rw [abs_le]; constructor <;> nlinarith [inv_pos.2 hp0]
  · rw [iteratedDeriv_const_add hi0, iteratedDeriv_eq_iterate, iter_deriv_inv, abs_mul, abs_mul,
      abs_pow, abs_neg, abs_one, one_pow, one_mul, Nat.abs_cast]
    have hz : |p ^ (-1 - (i : ℤ))| = (p⁻¹) ^ (i + 1) := by
      rw [abs_of_pos (zpow_pos hp0 _), show (-1 - (i : ℤ)) = -((i + 1 : ℕ) : ℤ) by push_cast; ring,
        zpow_neg, zpow_natCast, inv_pow]
    rw [hz]
    have h1 : (p⁻¹) ^ (i + 1) ≤ 2 ^ (i + 1) := pow_le_pow_left₀ (inv_pos.2 hp0).le hinv _
    have h2 : (2:ℝ) ^ (i + 1) ≤ 64 := by
      calc (2:ℝ) ^ (i + 1) ≤ 2 ^ 6 := pow_le_pow_right₀ (by norm_num) (by omega)
        _ = 64 := by norm_num
    have h3 : (i.factorial : ℝ) ≤ 120 := by
      have : i.factorial ≤ Nat.factorial 5 := Nat.factorial_le hi
      exact_mod_cast this.trans (by decide)
    calc (i.factorial : ℝ) * (p⁻¹) ^ (i + 1) ≤ 120 * 64 :=
          mul_le_mul h3 (h1.trans h2) (by positivity) (by norm_num)
      _ = 7680 := by norm_num

lemma fact_le5 {n : ℕ} (hn : n ≤ 5) : (n.factorial : ℝ) ≤ 120 := by
  have : n.factorial ≤ Nat.factorial 5 := Nat.factorial_le hn
  exact_mod_cast this.trans (by decide)

section Core
variable (w α β : Fin 3 → ℝ) (E0 E1 Y0 Y1 : ℝ)

lemma cF_cd5 : ContDiff ℝ 5 (cF w α) := by rw [cF_eq]; fun_prop
lemma cS_cd5 : ContDiff ℝ 5 (cS w α β) := by rw [cS_eq]; fun_prop

variable {w α β E0 E1 Y0 Y1}

lemma nd_three5 {f0 f1 f2 : ℝ → ℝ} (h0 : ContDiff ℝ 5 f0) (h1 : ContDiff ℝ 5 f1)
    (h2 : ContDiff ℝ 5 f2) {p : ℝ} (hp : 0 < p) {n : ℕ} (hn : n ≤ 5) :
    nd n (fun q => f0 q + f1 q + f2 q) p ≤ nd n f0 p + nd n f1 p + nd n f2 p := by
  refine (nd_add5 (h0.add h1).contDiffOn h2.contDiffOn hp hn).trans ?_
  have := nd_add5 h0.contDiffOn h1.contDiffOn hp hn
  linarith

lemma nd_cF5 (hw : ∀ j, 0 < w j) {D : ℝ} (hD : ∀ j, |Real.log (α j)| ≤ D) {p : ℝ} (hp : 0 < p)
    {n : ℕ} (hn : n ≤ 5) : nd n (cF w α) p ≤ cF w α p * D ^ n := by
  rw [cF_eq]
  refine (nd_three5 (by fun_prop) (by fun_prop) (by fun_prop) hp hn).trans ?_
  simp only [nd_base _ _ _ hp]
  have k : ∀ j, |w j| * |Real.log (α j)| ^ n * Real.exp (Real.log (α j) * p) ≤
      w j * Real.exp (Real.log (α j) * p) * D ^ n := by
    intro j
    rw [abs_of_pos (hw j)]
    have := pow_le_pow_left₀ (abs_nonneg _) (hD j) n
    have := hw j
    calc w j * |Real.log (α j)| ^ n * Real.exp (Real.log (α j) * p)
        ≤ w j * D ^ n * Real.exp (Real.log (α j) * p) := by gcongr
      _ = _ := by ring
  linarith [k 0, k 1, k 2]

lemma nd_cS5 {D : ℝ} (hD : ∀ j, |Real.log (α j)| ≤ D) {p : ℝ} (hp : 0 < p)
    {n : ℕ} (hn : n ≤ 5) : nd n (cS w α β) p ≤ cSt w α β p * D ^ n := by
  rw [cS_eq]
  refine (nd_three5 (by fun_prop) (by fun_prop) (by fun_prop) hp hn).trans ?_
  simp only [nd_base _ _ _ hp]
  have k : ∀ j, |w j * β j / α j| * |Real.log (α j)| ^ n * Real.exp (Real.log (α j) * p) ≤
      |w j * β j / α j| * Real.exp (Real.log (α j) * p) * D ^ n := by
    intro j
    have := pow_le_pow_left₀ (abs_nonneg _) (hD j) n
    calc |w j * β j / α j| * |Real.log (α j)| ^ n * Real.exp (Real.log (α j) * p)
        ≤ |w j * β j / α j| * D ^ n * Real.exp (Real.log (α j) * p) := by gcongr
      _ = _ := by ring
  have e : cSt w α β p = |w 0 * β 0 / α 0| * Real.exp (Real.log (α 0) * p) +
      |w 1 * β 1 / α 1| * Real.exp (Real.log (α 1) * p) +
      |w 2 * β 2 / α 2| * Real.exp (Real.log (α 2) * p) := by simp [cSt, Fin.sum_univ_three]
  rw [e]
  linarith [k 0, k 1, k 2]

lemma nd_cQ5 (hw : ∀ j, 0 < w j) {D : ℝ} (hD : ∀ j, |Real.log (α j)| ≤ D) {p : ℝ} (hp : 0 < p)
    {n : ℕ} (hn : n ≤ 5) :
    nd n (cQ w α β E0 E1 Y0 Y1) p ≤ 2 ^ n * ((|E0 ^ 2 + E1 ^ 2| + 2 * |E0 * Y0 + E1 * Y1| +
      |Y0 ^ 2 + Y1 ^ 2|) * (cF w α p + cSt w α β p) ^ 2) * D ^ n := by
  set F := cF w α
  set S := cS w α β
  have e : cQ w α β E0 E1 Y0 Y1 = fun q => (E0 ^ 2 + E1 ^ 2) * (F q * F q) +
      (-2 * (E0 * Y0 + E1 * Y1)) * (F q * S q) + (Y0 ^ 2 + Y1 ^ 2) * (S q * S q) := by
    ext q; simp only [cQ, F, S]; ring
  have hFc : ContDiff ℝ 5 F := cF_cd5 w α
  have hSc : ContDiff ℝ 5 S := cS_cd5 w α β
  rw [e]
  refine (nd_three5 (by fun_prop) (by fun_prop) (by fun_prop) hp hn).trans ?_
  simp only [nd_cmul]
  have hF : ∀ i ≤ n, nd i F p ≤ F p * D ^ i := fun i hi => nd_cF5 hw hD hp (hi.trans hn)
  have hS : ∀ i ≤ n, nd i S p ≤ cSt w α β p * D ^ i := fun i hi => nd_cS5 hD hp (hi.trans hn)
  have a1 := nd_mul5 hFc.contDiffOn hFc.contDiffOn hp hn hF hF
  have a2 := nd_mul5 hFc.contDiffOn hSc.contDiffOn hp hn hF hS
  have a3 := nd_mul5 hSc.contDiffOn hSc.contDiffOn hp hn hS hS
  have hFp : 0 < F p := cF_pos hw p
  have hSt := cSt_nonneg (w := w) (α := α) (β := β) p
  have hDn : 0 ≤ (2:ℝ) ^ n * D ^ n := by
    have : 0 ≤ nd 0 F p := nd_nonneg _ _ _
    have h0 := hF 0 (Nat.zero_le _)
    simp only [pow_zero, mul_one, nd_zero] at h0
    have hD0 : 0 ≤ D := (abs_nonneg _).trans (hD 0)
    positivity
  have := abs_nonneg (E0 ^ 2 + E1 ^ 2)
  have := abs_nonneg (E0 * Y0 + E1 * Y1)
  have := abs_nonneg (Y0 ^ 2 + Y1 ^ 2)
  rw [abs_mul, abs_neg, abs_two] at *
  calc |E0 ^ 2 + E1 ^ 2| * nd n (fun q => F q * F q) p +
        2 * |E0 * Y0 + E1 * Y1| * nd n (fun q => F q * S q) p +
        |Y0 ^ 2 + Y1 ^ 2| * nd n (fun q => S q * S q) p
      ≤ |E0 ^ 2 + E1 ^ 2| * (2 ^ n * (F p * F p) * D ^ n) +
        2 * |E0 * Y0 + E1 * Y1| * (2 ^ n * (F p * cSt w α β p) * D ^ n) +
        |Y0 ^ 2 + Y1 ^ 2| * (2 ^ n * (cSt w α β p * cSt w α β p) * D ^ n) := by gcongr
    _ ≤ _ := by
        have e2 : ∀ X : ℝ, 2 ^ n * X * D ^ n = X * (2 ^ n * D ^ n) := fun X => by ring
        rw [e2, e2, e2, e2]
        have : |E0 ^ 2 + E1 ^ 2| * (F p * F p) + 2 * |E0 * Y0 + E1 * Y1| * (F p * cSt w α β p) +
            |Y0 ^ 2 + Y1 ^ 2| * (cSt w α β p * cSt w α β p) ≤
            (|E0 ^ 2 + E1 ^ 2| + 2 * |E0 * Y0 + E1 * Y1| + |Y0 ^ 2 + Y1 ^ 2|) *
              (F p + cSt w α β p) ^ 2 := by
          nlinarith [mul_nonneg hFp.le hSt, mul_nonneg hSt hSt, mul_nonneg hFp.le hFp.le]
        nlinarith

lemma cQ_cd5 (hw : ∀ j, 0 < w j) : ContDiff ℝ 5 (cQ w α β E0 E1 Y0 Y1) := by
  have := cF_cd5 w α; have := cS_cd5 w α β
  unfold cQ; fun_prop

lemma nd_cQ_rel5 (hw : ∀ j, 0 < w j) {D : ℝ} (hD1 : 1 ≤ D) (hD : ∀ j, |Real.log (α j)| ≤ D)
    {p : ℝ} (hp : 0 < p) {K0 : ℝ} (hK0 : 0 ≤ K0)
    (hQ : (cF w α p + cSt w α β p) ^ 2 ≤ K0 * cQ w α β E0 E1 Y0 Y1 p)
    {i : ℕ} (hi1 : 1 ≤ i) (hi : i ≤ 5) :
    nd i (cQ w α β E0 E1 Y0 Y1) p ≤
      cQ w α β E0 E1 Y0 Y1 p * ((2 * (cM E0 E1 Y0 Y1 * K0 + 1)) * D) ^ i := by
  refine (nd_cQ5 hw hD hp hi).trans ?_
  have hM := cM_nonneg (E0 := E0) (E1 := E1) (Y0 := Y0) (Y1 := Y1)
  set Q := cQ w α β E0 E1 Y0 Y1 p
  set X := (cF w α p + cSt w α β p) ^ 2
  set M := cM E0 E1 Y0 Y1
  have h1 : M * K0 ≤ (M * K0 + 1) ^ i := by
    calc M * K0 ≤ M * K0 + 1 := by linarith
      _ = (M * K0 + 1) ^ 1 := (pow_one _).symm
      _ ≤ (M * K0 + 1) ^ i := pow_le_pow_right₀ (by nlinarith) hi1
  have hD0 : 0 ≤ D ^ i := by positivity
  rw [mul_pow, mul_pow]
  change 2 ^ i * (M * X) * D ^ i ≤ _
  have hQ0 : 0 ≤ Q := by simp only [Q, cQ]; positivity
  have : M * X ≤ Q * (M * K0 + 1) ^ i := by
    calc M * X ≤ M * (K0 * Q) := mul_le_mul_of_nonneg_left hQ hM
      _ = Q * (M * K0) := by ring
      _ ≤ Q * (M * K0 + 1) ^ i := mul_le_mul_of_nonneg_left h1 hQ0
  have h2 : (0:ℝ) ≤ 2 ^ i := by positivity
  calc 2 ^ i * (M * X) * D ^ i ≤ 2 ^ i * (Q * (M * K0 + 1) ^ i) * D ^ i := by gcongr
    _ = _ := by ring

lemma ch_cd5 (hw : ∀ j, 0 < w j) (hdet : E0 * Y1 - E1 * Y0 ≠ 0) :
    ContDiffOn ℝ 5 (ch w α β E0 E1 Y0 Y1) (Ioi 0) := by
  unfold ch
  have hF := cF_cd5 w α
  have hQ := cQ_cd5 (α := α) (β := β) (E0 := E0) (E1 := E1) (Y0 := Y0) (Y1 := Y1) hw
  refine ContDiffOn.add ?_ ?_
  · refine ContDiffOn.mul ?_ ?_
    · exact (contDiffOn_const.add (contDiffOn_inv ℝ |>.mono fun x hx => ne_of_gt hx)).neg
    · exact (hF.log fun q => (cF_pos hw q).ne').contDiffOn
  · exact contDiffOn_const.mul (hQ.log fun q => (cQ_pos hw hdet q).ne').contDiffOn

/-- the constant `K₅` in the bound for `log G`. -/
noncomputable def cK5 (b K2 : ℝ) : ℝ := 32 * 7680 * (|b| + 2880) + 1440 * K2 ^ 5 + 1

lemma nd_ch5 (hw : ∀ j, 0 < w j) (hdet : E0 * Y1 - E1 * Y0 ≠ 0) {D : ℝ} (hD1 : 1 ≤ D)
    (hD : ∀ j, |Real.log (α j)| ≤ D) {p : ℝ} (hp : 1 / 2 ≤ p) {K0 : ℝ} (hK0 : 0 ≤ K0)
    (hQ : (cF w α p + cSt w α β p) ^ 2 ≤ K0 * cQ w α β E0 E1 Y0 Y1 p)
    {n : ℕ} (hn1 : 1 ≤ n) (hn : n ≤ 5) :
    nd n (ch w α β E0 E1 Y0 Y1) p ≤
      cK5 (Real.log (cF w α p)) (2 * (cM E0 E1 Y0 Y1 * K0 + 1)) * D ^ n := by
  have hp0 : 0 < p := by linarith
  set K2 := 2 * (cM E0 E1 Y0 Y1 * K0 + 1)
  have hK2 : 1 ≤ K2 := by have := cM_nonneg (E0 := E0) (E1 := E1) (Y0 := Y0) (Y1 := Y1); nlinarith
  have hF := cF_cd5 w α
  have hQc := cQ_cd5 (α := α) (β := β) (E0 := E0) (E1 := E1) (Y0 := Y0) (Y1 := Y1) hw
  have hlF : ContDiffOn ℝ 5 (fun q => Real.log (cF w α q)) (Ioi 0) :=
    (hF.log fun q => (cF_pos hw q).ne').contDiffOn
  have hlQ : ContDiffOn ℝ 5 (fun q => Real.log (cQ w α β E0 E1 Y0 Y1 q)) (Ioi 0) :=
    (hQc.log fun q => (cQ_pos hw hdet q).ne').contDiffOn
  have hinv : ContDiffOn ℝ 5 (fun q : ℝ => -(1 + q⁻¹)) (Ioi 0) :=
    (contDiffOn_const.add (contDiffOn_inv ℝ |>.mono fun x hx => ne_of_gt hx)).neg
  have hDpow : ∀ i : ℕ, 1 ≤ D ^ i := fun i => one_le_pow₀ hD1
  unfold ch
  refine (nd_add5 (hinv.mul hlF) (contDiffOn_const.mul hlQ) hp0 hn).trans ?_
  rw [nd_cmul]
  have hA : ∀ i ≤ n, nd i (fun q : ℝ => -(1 + q⁻¹)) p ≤ 7680 * D ^ i := fun i hi => by
    have := nd_inv_part5 i (hi.trans hn) hp
    nlinarith [hDpow i]
  have hB : ∀ i ≤ n, nd i (fun q => Real.log (cF w α q)) p ≤
      (|Real.log (cF w α p)| + 2880) * D ^ i := by
    intro i hi
    rcases Nat.eq_zero_or_pos i with rfl | hi0
    · rw [nd_zero]; simp
    · have := nd_log5 hF.contDiffOn (fun q _ => cF_pos hw q) hp0 hi0 (hi.trans hn)
        (fun j _ hj => nd_cF5 hw hD hp0 (hj.trans (hi.trans hn)))
      have hf := fact_le5 (hi.trans hn)
      have : (i.factorial : ℝ) * 24 * D ^ i ≤ 2880 * D ^ i := by nlinarith [hDpow i]
      nlinarith [abs_nonneg (Real.log (cF w α p)), hDpow i]
  have h1 := nd_mul5 hinv hlF hp0 hn hA hB
  have h2 := nd_log5 hQc.contDiffOn (fun q _ => cQ_pos hw hdet q) hp0 hn1 hn
    (fun i hi1 hi => nd_cQ_rel5 hw hD1 hD hp0 hK0 hQ hi1 (hi.trans hn))
  have hf := fact_le5 hn
  have h2' : |(1:ℝ) / 2| * nd n (fun q => Real.log (cQ w α β E0 E1 Y0 Y1 q)) p ≤
      1440 * K2 ^ 5 * D ^ n := by
    rw [mul_pow] at h2
    have hk : K2 ^ n ≤ K2 ^ 5 := pow_le_pow_right₀ hK2 hn
    have : (n.factorial : ℝ) * 24 * (K2 ^ n * D ^ n) ≤ 2880 * K2 ^ 5 * D ^ n := by
      have := hDpow n
      have : 0 ≤ K2 ^ n := by positivity
      nlinarith [mul_le_mul hf hk (by positivity) (by norm_num)]
    rw [abs_of_pos (by norm_num : (0:ℝ) < 1 / 2)]
    linarith
  have h1' : 2 ^ n * (7680 * (|Real.log (cF w α p)| + 2880)) * D ^ n ≤
      32 * 7680 * (|Real.log (cF w α p)| + 2880) * D ^ n := by
    have : (2:ℝ) ^ n ≤ 32 := by
      calc (2:ℝ) ^ n ≤ 2 ^ 5 := pow_le_pow_right₀ (by norm_num) hn
        _ = 32 := by norm_num
    have hX : 0 ≤ 7680 * (|Real.log (cF w α p)| + 2880) * D ^ n := by
      have := hDpow n; positivity
    calc 2 ^ n * (7680 * (|Real.log (cF w α p)| + 2880)) * D ^ n
        = 2 ^ n * (7680 * (|Real.log (cF w α p)| + 2880) * D ^ n) := by ring
      _ ≤ 32 * (7680 * (|Real.log (cF w α p)| + 2880) * D ^ n) := mul_le_mul_of_nonneg_right this hX
      _ = _ := by ring
  unfold cK5
  nlinarith [hDpow n]

lemma cK5_one (b : ℝ) {K2 : ℝ} (hK : 0 ≤ K2) : 1 ≤ cK5 b K2 := by
  unfold cK5; have := pow_nonneg hK 5; nlinarith [abs_nonneg b]

lemma cK5_mono {b b' K2 K2' : ℝ} (hb : |b| ≤ b') (hK : 0 ≤ K2) (hK' : K2 ≤ K2') :
    cK5 b K2 ≤ cK5 b' K2' := by
  unfold cK5
  have : K2 ^ 5 ≤ K2' ^ 5 := pow_le_pow_left₀ hK hK' 5
  have : |b| ≤ |b'| := hb.trans (le_abs_self _)
  nlinarith

lemma nd_cG5 (hw : ∀ j, 0 < w j) (hdet : E0 * Y1 - E1 * Y0 ≠ 0) {D : ℝ} (hD1 : 1 ≤ D)
    (hD : ∀ j, |Real.log (α j)| ≤ D) {p : ℝ} (hp : 1 / 2 ≤ p) {K0 : ℝ} (hK0 : 0 ≤ K0)
    (hQ : (cF w α p + cSt w α β p) ^ 2 ≤ K0 * cQ w α β E0 E1 Y0 Y1 p)
    {n : ℕ} (hn : n ≤ 5) :
    nd n (cG w α β E0 E1 Y0 Y1) p ≤ 120 * cG w α β E0 E1 Y0 Y1 p *
      (cK5 (Real.log (cF w α p)) (2 * (cM E0 E1 Y0 Y1 * K0 + 1)) * D) ^ n := by
  have hp0 : 0 < p := by linarith
  set K3 := cK5 (Real.log (cF w α p)) (2 * (cM E0 E1 Y0 Y1 * K0 + 1))
  have hK3 : 1 ≤ K3 := by
    have : 0 ≤ (2 * (cM E0 E1 Y0 Y1 * K0 + 1)) := by
      have := cM_nonneg (E0 := E0) (E1 := E1) (Y0 := Y0) (Y1 := Y1); nlinarith
    exact cK5_one _ this
  have := nd_exp5 (ch_cd5 hw hdet) hp0 hn (D := K3 * D) (fun i hi1 hi => by
    refine (nd_ch5 hw hdet hD1 hD hp hK0 hQ hi1 (hi.trans hn)).trans ?_
    rw [mul_pow]
    have : K3 ≤ K3 ^ i := by
      calc K3 = K3 ^ 1 := (pow_one _).symm
        _ ≤ K3 ^ i := pow_le_pow_right₀ hK3 hi1
    have : 0 ≤ D ^ i := by positivity
    nlinarith)
  have hf := fact_le5 hn
  unfold cG
  refine this.trans ?_
  have h0 : 0 ≤ (K3 * D) ^ n := pow_nonneg (mul_nonneg (by linarith) (by linarith)) n
  have := exp_pos (ch w α β E0 E1 Y0 Y1 p)
  exact mul_le_mul_of_nonneg_right (mul_le_mul_of_nonneg_right hf this.le) h0

end Core

end P03.P5

namespace P03.P5
open S05 Set Real Filter Topology MeasureTheory
open scoped Interval

lemma pow5_le {x : ℝ} (hx : 0 ≤ x) : (1 + x) ^ 5 ≤ 32 * (1 + x ^ 5) := by
  rcases le_total x 1 with h | h
  · have : (1 + x) ^ 5 ≤ 2 ^ 5 := pow_le_pow_left₀ (by linarith) (by linarith) 5
    have : 0 ≤ x ^ 5 := pow_nonneg hx 5
    nlinarith
  · have : (1 + x) ^ 5 ≤ (2 * x) ^ 5 := pow_le_pow_left₀ (by linarith) (by linarith) 5
    rw [mul_pow] at this
    nlinarith

lemma gen_bound5 (w : Fin 3 → ℝ) (hw : ∀ j, 0 < w j) {a b m2 M2 c2 d2 CY CE : ℝ}
    (ha : 0 < a) (hb : 0 < b) (hm2 : 0 < m2) (hc2 : 0 ≤ c2) (hd2 : 0 < d2) (hCY : 0 ≤ CY)
    (hCE : 0 ≤ CE) :
    ∃ K, ∀ (t p : ℝ) (α β : Fin 3 → ℝ) (E0 E1 Y0 Y1 : ℝ), t ∈ Ioo (0:ℝ) 1 →
      p ∈ Icc (3/4:ℝ) (5/4) → α 0 = a * t → α 1 = b * (1 - t) → m2 ≤ α 2 → α 2 ≤ M2 →
      β 0 = a → β 1 = -b → |β 2| ≤ c2 * α 2 → d2 ≤ (E0 * Y1 - E1 * Y0) ^ 2 →
      Y0 ^ 2 + Y1 ^ 2 ≤ CY → E0 ^ 2 + E1 ^ 2 ≤ CE →
      ∀ n ≤ 5, nd n (cG w α β E0 E1 Y0 Y1) p ≤ K * (t * (1 - t)) ^ (-(3/4:ℝ)) := by
  have hw0 := hw 0; have hw1 := hw 1; have hw2 := hw 2
  obtain ⟨U, hU⟩ : ∃ U, U = 1 + a + b + |M2| := ⟨_, rfl⟩
  have hU1 : 1 ≤ U := by have := abs_nonneg M2; linarith
  obtain ⟨L0, hL0d⟩ : ∃ L0, L0 = w 2 * (min 1 m2) ^ 2 := ⟨_, rfl⟩
  have hL0 : 0 < L0 := by rw [hL0d]; have : 0 < min 1 m2 := lt_min one_pos hm2; positivity
  obtain ⟨U0, hU0⟩ : ∃ U0, U0 = (w 0 + w 1 + w 2) * U ^ 2 := ⟨_, rfl⟩
  have hU00 : 0 ≤ U0 := by rw [hU0]; positivity
  obtain ⟨CF, hCF⟩ : ∃ CF, CF = |Real.log L0| + |Real.log U0| := ⟨_, rfl⟩
  have hCF0 : 0 ≤ CF := by rw [hCF]; positivity
  obtain ⟨Kd, hKdd⟩ : ∃ Kd, Kd = 1 + |Real.log a| + |Real.log b| + |Real.log m2| +
    |Real.log M2| := ⟨_, rfl⟩
  have hKd : 1 ≤ Kd := by
    rw [hKdd]; have := abs_nonneg (Real.log a); have := abs_nonneg (Real.log b)
    have := abs_nonneg (Real.log m2); have := abs_nonneg (Real.log M2); linarith
  obtain ⟨cc, hccd⟩ : ∃ cc, cc = 4 + 2 * c2 := ⟨_, rfl⟩
  have hcc : 1 ≤ cc + 1 := by rw [hccd]; linarith
  obtain ⟨K0, hK0d⟩ : ∃ K0, K0 = 2 * (cc + 1) ^ 2 * (CE + CY) / d2 := ⟨_, rfl⟩
  have hK0 : 0 ≤ K0 := by rw [hK0d]; positivity
  obtain ⟨CM, hCMd⟩ : ∃ CM, CM = 2 * (CE + CY) := ⟨_, rfl⟩
  have hCM : 0 ≤ CM := by rw [hCMd]; positivity
  obtain ⟨K3, hK3d⟩ : ∃ K3, K3 = cK5 CF (2 * (CM * K0 + 1)) := ⟨_, rfl⟩
  have hK3 : 1 ≤ K3 := by rw [hK3d]; exact cK5_one _ (by positivity)
  obtain ⟨CS, hCSd⟩ : ∃ CS, CS = (w 0 + w 1 + w 2 * c2) * U ^ 2 := ⟨_, rfl⟩
  have hCS : 0 ≤ CS := by rw [hCSd]; positivity
  obtain ⟨CG, hCGd⟩ : ∃ CG, CG = Real.exp (3 * CF) * Real.sqrt CM * (U0 + CS) := ⟨_, rfl⟩
  have hCG : 0 ≤ CG := by rw [hCGd]; positivity
  refine ⟨120 * CG * (K3 * Kd) ^ 5 * (32 * 7962625), ?_⟩
  intro t p α β E0 E1 Y0 Y1 ht hp hα0 hα1 hm hM hβ0 hβ1 hβ2 hdd hY hE n hn
  have hp0 : 0 < p := by linarith [hp.1]
  have ht0 := ht.1
  have ht1 : 0 < 1 - t := by linarith [ht.2]
  have hpos : ∀ j, 0 < α j := by
    intro j; fin_cases j
    · simp only [Fin.zero_eta]; rw [hα0]; positivity
    · simp only [Fin.mk_one]; rw [hα1]; positivity
    · exact hm2.trans_le hm
  have hαU : ∀ j, α j ≤ U := by
    have := abs_nonneg M2
    intro j; fin_cases j
    · simp only [Fin.zero_eta]; rw [hα0]
      calc a * t ≤ a * 1 := mul_le_mul_of_nonneg_left ht.2.le ha.le
        _ ≤ U := by rw [hU]; linarith
    · simp only [Fin.mk_one]; rw [hα1]
      calc b * (1 - t) ≤ b * 1 := mul_le_mul_of_nonneg_left (by linarith) hb.le
        _ ≤ U := by rw [hU]; linarith
    · exact hM.trans ((le_abs_self M2).trans (by linarith))
  have hex : ∀ j, Real.exp (Real.log (α j) * p) = α j ^ p := fun j =>
    (Real.rpow_def_of_pos (hpos j) p).symm
  obtain ⟨F, hFdef⟩ : ∃ F, F = cF w α p := ⟨_, rfl⟩
  have hF3 : F = w 0 * α 0 ^ p + w 1 * α 1 ^ p + w 2 * α 2 ^ p := by
    simp only [hFdef, cF, Fin.sum_univ_three, hex]
  have hpow : ∀ j, α j ^ p ≤ U ^ 2 := fun j =>
    rpow_upper (hpos j) (hαU j) hU1 hp0.le (by linarith [hp.2])
  have hpow0 : ∀ j, 0 ≤ α j ^ p := fun j => (Real.rpow_pos_of_pos (hpos j) p).le
  have hFU : F ≤ U0 := by
    rw [hF3, hU0]
    have := mul_le_mul_of_nonneg_left (hpow 0) hw0.le
    have := mul_le_mul_of_nonneg_left (hpow 1) hw1.le
    have := mul_le_mul_of_nonneg_left (hpow 2) hw2.le
    have : (w 0 + w 1 + w 2) * U ^ 2 = w 0 * U ^ 2 + w 1 * U ^ 2 + w 2 * U ^ 2 := by ring
    linarith
  have hFL : L0 ≤ F := by
    rw [hF3, hL0d]
    have := mul_le_mul_of_nonneg_left (rpow_lower hm2 hm hp0.le (by linarith [hp.2])) hw2.le
    have := mul_nonneg hw0.le (hpow0 0); have := mul_nonneg hw1.le (hpow0 1)
    linarith
  have hlogF : |Real.log F| ≤ CF := hCF ▸ abs_log_le_of hL0 hFL hFU
  -- `D`
  obtain ⟨D, hDd⟩ : ∃ D, D = Kd * (1 + ls t) := ⟨_, rfl⟩
  have hls := ls_nonneg t
  have hD1 : 1 ≤ D := by rw [hDd]; exact one_le_mul_of_one_le_of_one_le hKd (by linarith)
  have hKdls : Kd + ls t ≤ D := by
    rw [hDd]; have := mul_le_mul_of_nonneg_right hKd hls; linarith
  have hD : ∀ j, |Real.log (α j)| ≤ D := by
    have := abs_nonneg (Real.log a); have := abs_nonneg (Real.log b)
    have := abs_nonneg (Real.log m2); have := abs_nonneg (Real.log M2)
    intro j; fin_cases j
    · simp only [Fin.zero_eta]; rw [hα0, Real.log_mul ha.ne' ht0.ne']
      refine (abs_add_le _ _).trans ?_
      have : |Real.log t| ≤ ls t := abs_log_le_ls t
      rw [hKdd] at hKdls; linarith
    · simp only [Fin.mk_one]; rw [hα1, Real.log_mul hb.ne' ht1.ne']
      refine (abs_add_le _ _).trans ?_
      have : |Real.log (1 - t)| ≤ ls t := abs_log1_le_ls t
      rw [hKdd] at hKdls; linarith
    · show |Real.log (α 2)| ≤ D
      have := abs_log_le_of hm2 hm hM
      rw [hKdd] at hKdls; linarith
  -- `S`, `St`
  obtain ⟨S, hSd⟩ : ∃ S, S = cS w α β p := ⟨_, rfl⟩
  obtain ⟨St, hStd⟩ : ∃ St, St = cSt w α β p := ⟨_, rfl⟩
  have hx0 : 0 ≤ w 0 * α 0 ^ p := mul_nonneg hw0.le (hpow0 0)
  have hx1 : 0 ≤ w 1 * α 1 ^ p := mul_nonneg hw1.le (hpow0 1)
  have hx2 : 0 ≤ w 2 * α 2 ^ p := mul_nonneg hw2.le (hpow0 2)
  have hS3 : S = (w 0 * α 0 ^ p) / t - (w 1 * α 1 ^ p) / (1 - t) +
      (β 2 / α 2) * (w 2 * α 2 ^ p) := by
    rw [hSd]; simp only [cS, Fin.sum_univ_three, hex, hβ0, hβ1]
    rw [hα0, hα1]
    field_simp
    ring
  have hSt3 : St = (w 0 * α 0 ^ p) / t + (w 1 * α 1 ^ p) / (1 - t) +
      |β 2 / α 2| * (w 2 * α 2 ^ p) := by
    rw [hStd]; simp only [cSt, Fin.sum_univ_three, hex, hβ0, hβ1]
    rw [hα0, hα1]
    have q0 : w 0 * a / (a * t) = w 0 / t := by field_simp
    have q1 : w 1 * -b / (b * (1 - t)) = -(w 1 / (1 - t)) := by field_simp
    have q2 : |w 2 * β 2 / α 2| = w 2 * |β 2 / α 2| := by rw [mul_div_assoc, abs_mul, abs_of_pos hw2]
    rw [q0, q1, q2, abs_neg, abs_of_pos (div_pos hw0 ht0), abs_of_pos (div_pos hw1 ht1)]; ring
  have hq2 : |β 2 / α 2| ≤ c2 := by
    rw [abs_div, abs_of_pos (hpos 2), div_le_iff₀ (hpos 2)]; exact hβ2
  have hStle : St ≤ |S| + cc * F := by
    rw [hSt3, hS3, hF3, hccd]
    exact St_le hx0 hx1 hx2 ht hq2
  have hSt0 : 0 ≤ St := hStd ▸ cSt_nonneg p
  have hFpos : 0 < F := hL0.trans_le hFL
  obtain ⟨Q, hQd⟩ : ∃ Q, Q = cQ w α β E0 E1 Y0 Y1 p := ⟨_, rfl⟩
  have hdd0 : E0 * Y1 - E1 * Y0 ≠ 0 := by
    intro h; rw [h] at hdd; norm_num at hdd; linarith
  have hQpos : 0 < Q := hQd ▸ cQ_pos hw hdd0 p
  have hQlow := Q_lower F S E0 E1 Y0 Y1
  have hQe : Q = (F * E0 - S * Y0) ^ 2 + (F * E1 - S * Y1) ^ 2 := by rw [hQd, hFdef, hSd]; rfl
  rw [← hQe] at hQlow
  have hFS : F ^ 2 + S ^ 2 ≤ Q * (CE + CY) / d2 := by
    rw [le_div_iff₀ hd2]
    have h1 : d2 * (F ^ 2 + S ^ 2) ≤ (E0 * Y1 - E1 * Y0) ^ 2 * (F ^ 2 + S ^ 2) :=
      mul_le_mul_of_nonneg_right hdd (by positivity)
    have h2 : Q * ((E0 ^ 2 + E1 ^ 2) + (Y0 ^ 2 + Y1 ^ 2)) ≤ Q * (CE + CY) :=
      mul_le_mul_of_nonneg_left (by linarith) hQpos.le
    linarith
  have hQ : (F + St) ^ 2 ≤ K0 * Q := by
    have h1 : F + St ≤ |S| + (cc + 1) * F := by linarith
    have h2 : (F + St) ^ 2 ≤ (|S| + (cc + 1) * F) ^ 2 :=
      pow_le_pow_left₀ (by linarith) h1 2
    have h3 : (|S| + (cc + 1) * F) ^ 2 ≤ 2 * (cc + 1) ^ 2 * (F ^ 2 + S ^ 2) :=
      sq_sum_le S (cc + 1) F hcc
    have h4 : 2 * (cc + 1) ^ 2 * (F ^ 2 + S ^ 2) ≤ 2 * (cc + 1) ^ 2 * (Q * (CE + CY) / d2) :=
      mul_le_mul_of_nonneg_left hFS (by positivity)
    have : 2 * (cc + 1) ^ 2 * (Q * (CE + CY) / d2) = K0 * Q := by rw [hK0d]; ring
    linarith
  have hcM : cM E0 E1 Y0 Y1 ≤ CM := by
    unfold cM
    rw [abs_of_nonneg (by positivity : (0:ℝ) ≤ E0 ^ 2 + E1 ^ 2),
      abs_of_nonneg (by positivity : (0:ℝ) ≤ Y0 ^ 2 + Y1 ^ 2)]
    have := dot_abs_le E0 E1 Y0 Y1
    rw [hCMd]; linarith
  have hcM0 : 0 ≤ cM E0 E1 Y0 Y1 := cM_nonneg
  -- the derivative bound
  rw [hFdef, hStd, hQd] at hQ
  have hnd := nd_cG5 hw hdd0 hD1 hD (by linarith [hp.1]) hK0 hQ hn
  rw [← hFdef] at hnd
  have hk3 : cK5 (Real.log F) (2 * (cM E0 E1 Y0 Y1 * K0 + 1)) ≤ K3 := by
    rw [hK3d]; refine cK5_mono hlogF (by positivity) ?_
    have := mul_le_mul_of_nonneg_right hcM hK0; linarith
  have hk3' : 1 ≤ cK5 (Real.log F) (2 * (cM E0 E1 Y0 Y1 * K0 + 1)) := cK5_one _ (by positivity)
  have hpowK : (cK5 (Real.log F) (2 * (cM E0 E1 Y0 Y1 * K0 + 1)) * D) ^ n ≤ (K3 * D) ^ 5 := by
    calc (cK5 (Real.log F) (2 * (cM E0 E1 Y0 Y1 * K0 + 1)) * D) ^ n ≤ (K3 * D) ^ n :=
          pow_le_pow_left₀ (by positivity) (mul_le_mul_of_nonneg_right hk3 (by linarith)) n
      _ ≤ (K3 * D) ^ 5 := pow_le_pow_right₀ (one_le_mul_of_one_le_of_one_le hK3 hD1) hn
  -- the size of `G`
  set T := (t * (1 - t)) ^ (-(1/4:ℝ)) with hT
  have hT1 : 1 ≤ T := one_le_tt ht (by norm_num)
  have hStT : St ≤ CS * T := by
    rw [hSt3, hα0, hα1, Real.mul_rpow ha.le ht0.le, Real.mul_rpow hb.le ht1.le, hCSd]
    have e0 : w 0 * (a ^ p * t ^ p) / t = w 0 * a ^ p * t ^ (p - 1) := by
      rw [Real.rpow_sub_one ht0.ne']; ring
    have e1 : w 1 * (b ^ p * (1 - t) ^ p) / (1 - t) = w 1 * b ^ p * (1 - t) ^ (p - 1) := by
      rw [Real.rpow_sub_one ht1.ne']; ring
    rw [e0, e1]
    have ha' : a ^ p ≤ U ^ 2 := rpow_upper ha (by rw [hU]; have := abs_nonneg M2; linarith) hU1
      hp0.le (by linarith [hp.2])
    have hb' : b ^ p ≤ U ^ 2 := rpow_upper hb (by rw [hU]; have := abs_nonneg M2; linarith) hU1
      hp0.le (by linarith [hp.2])
    have t0 : t ^ (p - 1) ≤ T := tpow_le ht hp.1
    have t1 : (1 - t) ^ (p - 1) ≤ T := by
      have := tpow_le (t := 1 - t) ⟨ht1, by linarith⟩ hp.1
      rwa [show (1 - t) * (1 - (1 - t)) = t * (1 - t) by ring] at this
    have e2 : |β 2 / α 2| * (w 2 * α 2 ^ p) ≤ c2 * (w 2 * U ^ 2) * T := by
      have := hpow 2
      calc |β 2 / α 2| * (w 2 * α 2 ^ p) ≤ c2 * (w 2 * U ^ 2) :=
            mul_le_mul hq2 (mul_le_mul_of_nonneg_left this hw2.le) (by positivity) hc2
        _ ≤ c2 * (w 2 * U ^ 2) * T := le_mul_of_one_le_right (by positivity) hT1
    have := Real.rpow_pos_of_pos ha p
    have := Real.rpow_pos_of_pos hb p
    have := Real.rpow_pos_of_pos ht0 (p - 1)
    have := Real.rpow_pos_of_pos ht1 (p - 1)
    have f0 : w 0 * a ^ p * t ^ (p - 1) ≤ w 0 * U ^ 2 * T := by gcongr
    have f1 : w 1 * b ^ p * (1 - t) ^ (p - 1) ≤ w 1 * U ^ 2 * T := by gcongr
    have : (w 0 + w 1 + w 2 * c2) * U ^ 2 * T =
      w 0 * U ^ 2 * T + w 1 * U ^ 2 * T + c2 * (w 2 * U ^ 2) * T := by ring
    linarith
  have hQup : Q ≤ CM * (F + St) ^ 2 := by
    have := nd_cQ5 (α := α) (β := β) (E0 := E0) (E1 := E1) (Y0 := Y0) (Y1 := Y1) hw hD hp0
      (n := 0) (by norm_num)
    rw [nd_zero, pow_zero, pow_zero, one_mul, mul_one, ← hQd, ← hFdef, ← hStd,
      abs_of_pos hQpos] at this
    refine this.trans (mul_le_mul_of_nonneg_right ?_ (sq_nonneg _))
    exact hcM
  have hGle : cG w α β E0 E1 Y0 Y1 p ≤ CG * T := by
    unfold cG ch
    rw [← hFdef, ← hQd, Real.exp_add]
    have hinv : p⁻¹ ≤ 2 := by rw [inv_le_comm₀ hp0 (by norm_num)]; linarith [hp.1]
    have h1 : Real.exp (-(1 + p⁻¹) * Real.log F) ≤ Real.exp (3 * CF) := by
      apply Real.exp_le_exp.2
      have hpi : 0 < p⁻¹ := inv_pos.2 hp0
      have h1 := mul_le_mul_of_nonneg_left (neg_abs_le (Real.log F)) (by linarith : (0:ℝ) ≤ 1 + p⁻¹)
      have h2 := mul_le_mul (by linarith : 1 + p⁻¹ ≤ 3) hlogF (abs_nonneg _) (by norm_num)
      have : -(1 + p⁻¹) * Real.log F = -((1 + p⁻¹) * Real.log F) := by ring
      linarith
    have h2 : Real.exp (1 / 2 * Real.log Q) = Real.sqrt Q := by
      rw [Real.sqrt_eq_rpow, Real.rpow_def_of_pos hQpos, mul_comm]
    rw [h2]
    have h3 : Real.sqrt Q ≤ Real.sqrt CM * (F + St) := by
      rw [← Real.sqrt_sq (by linarith : 0 ≤ F + St), ← Real.sqrt_mul hCM]
      exact Real.sqrt_le_sqrt hQup
    have h4 : F + St ≤ (U0 + CS) * T := by
      have := le_mul_of_one_le_right hU00 hT1
      have : (U0 + CS) * T = U0 * T + CS * T := by ring
      linarith
    rw [hCGd]
    calc Real.exp (-(1 + p⁻¹) * Real.log F) * Real.sqrt Q
        ≤ Real.exp (3 * CF) * (Real.sqrt CM * ((U0 + CS) * T)) :=
          mul_le_mul h1 (h3.trans (mul_le_mul_of_nonneg_left h4 (Real.sqrt_nonneg _)))
            (Real.sqrt_nonneg _) (exp_pos _).le
      _ = _ := by ring
  -- the log factor
  have hls3 : (1 + ls t) ^ 5 ≤ 32 * 7962625 * (t * (1 - t)) ^ (-(1/2:ℝ)) := by
    have h1 := ls_pow_le ht 5
    norm_num at h1
    have h2 : 1 ≤ (t * (1 - t)) ^ (-(1/2:ℝ)) := one_le_tt ht (by norm_num)
    have h3 : (1 + ls t) ^ 5 ≤ 32 * (1 + ls t ^ 5) := pow5_le hls
    linarith
  have hTT : T * (t * (1 - t)) ^ (-(1/2:ℝ)) = (t * (1 - t)) ^ (-(3/4:ℝ)) := by
    rw [hT, ← Real.rpow_add (mul_pos ht0 ht1)]; norm_num
  have hG0 : 0 ≤ cG w α β E0 E1 Y0 Y1 p := (exp_pos _).le
  calc nd n (cG w α β E0 E1 Y0 Y1) p
      ≤ 120 * cG w α β E0 E1 Y0 Y1 p *
          (cK5 (Real.log F) (2 * (cM E0 E1 Y0 Y1 * K0 + 1)) * D) ^ n := hnd
    _ ≤ 120 * (CG * T) * ((K3 * Kd) ^ 5 * (1 + ls t) ^ 5) := by
        rw [show (K3 * Kd) ^ 5 * (1 + ls t) ^ 5 = (K3 * D) ^ 5 by rw [hDd]; ring]
        gcongr
    _ ≤ 120 * (CG * T) * ((K3 * Kd) ^ 5 * (32 * 7962625 * (t * (1 - t)) ^ (-(1/2:ℝ)))) := by
        gcongr
    _ = 120 * CG * (K3 * Kd) ^ 5 * (32 * 7962625) * (T * (t * (1 - t)) ^ (-(1/2:ℝ))) := by ring
    _ = _ := by rw [hTT]

lemma edge_bound5 (H : Hex) (hH : IsHex H) (k : Fin 3) :
    ∃ K, ∀ p ∈ Icc (3/4:ℝ) (5/4), ∀ t ∈ Ioo (0:ℝ) 1, ∀ n ≤ 5,
      nd n (fun q => Gf H k q t) p ≤ K * (t * (1 - t)) ^ (-(3/4:ℝ)) := by
  obtain ⟨-, -, -, hB, hB1, hA0, hA1, hr⟩ := S05.proof_8 H hH k
  have hw := hH.1
  have hmr : 0 < min 1 (rr H k) := lt_min one_pos hr
  set V0 := eV0 H k
  set E := eE H k
  have hdet := edge_det_pos H hH k
  obtain ⟨K, hK⟩ := gen_bound5 (ew H k) (ew_pos H hH k) (a := A0 H k / H.w k)
    (b := A1 H k / H.w (k + 1)) (m2 := Bv H k * min 1 (rr H k) / H.w (k + 2))
    (M2 := Bv H k * max 1 (rr H k) / H.w (k + 2)) (c2 := |rr H k - 1| / min 1 (rr H k))
    (d2 := (det2 (eV0 H k) (eV1 H k)) ^ 2) (CY := 2 * (V0 0 ^ 2 + V0 1 ^ 2) + 2 * (E 0 ^ 2 + E 1 ^ 2))
    (CE := E 0 ^ 2 + E 1 ^ 2) (div_pos hA0 (hw _)) (div_pos hA1 (hw _))
    (div_pos (mul_pos hB hmr) (hw _)) (by positivity) (by positivity) (by positivity)
    (by positivity)
  refine ⟨K, fun p hp t ht n hn => ?_⟩
  have ht0 := ht.1.le
  have ht1 := ht.2.le
  have hml := min_le_aff (r := rr H k) ht0 ht1
  have hw2 := hw (k + 2)
  refine hK t p (eα H k t) (eβ H k) (E 0) (E 1) (ept H k t 0) (ept H k t 1) ht hp
    ?_ ?_ ?_ ?_ ?_ ?_ ?_ ?_ ?_ ?_ n hn
  · simp [eα]; ring
  · simp [eα]; ring
  · simp only [eα, Matrix.cons_val_two, Matrix.tail_cons, Matrix.head_cons]
    gcongr
  · simp only [eα, Matrix.cons_val_two, Matrix.tail_cons, Matrix.head_cons]
    gcongr
    exact aff_le_max ht0 ht1
  · simp [eα, eβ]
  · simp [eα, eβ]
  · simp only [eα, eβ, Matrix.cons_val_two, Matrix.tail_cons, Matrix.head_cons]
    have : 0 ≤ |rr H k - 1| := abs_nonneg _
    calc |Bv H k * (rr H k - 1) / H.w (k + 2)|
        = |rr H k - 1| * Bv H k * min 1 (rr H k) / (min 1 (rr H k) * H.w (k + 2)) := by
          rw [abs_div, abs_mul, abs_of_pos hB, abs_of_pos hw2]; field_simp
      _ ≤ |rr H k - 1| * Bv H k * (1 + (rr H k - 1) * t) / (min 1 (rr H k) * H.w (k + 2)) := by
          gcongr
      _ = _ := by field_simp
  · have e : E 0 * ept H k t 1 - E 1 * ept H k t 0 = -det2 (eV0 H k) (eV1 H k) := by
      simp only [E, eE, ept, det2, PiLp.add_apply, PiLp.sub_apply, PiLp.smul_apply, smul_eq_mul]
      ring
    rw [e, neg_sq]
  · have e0 : ept H k t 0 = V0 0 + t * E 0 := by
      simp [ept, V0, E, eE]
    have e1 : ept H k t 1 = V0 1 + t * E 1 := by
      simp [ept, V0, E, eE]
    rw [e0, e1]
    have : t ^ 2 ≤ 1 := by nlinarith
    nlinarith [sq_nonneg (V0 0 - t * E 0), sq_nonneg (V0 1 - t * E 1), sq_nonneg (E 0),
      sq_nonneg (E 1), mul_le_mul_of_nonneg_right this (sq_nonneg (E 0)),
      mul_le_mul_of_nonneg_right this (sq_nonneg (E 1))]
  · exact le_rfl

/-- `C^N` from a chain of derivatives on an open set. -/
lemma chainN (N : ℕ) (f : ℕ → ℝ → ℝ) {U : Set ℝ} (hU : IsOpen U)
    (hd : ∀ j < N, ∀ y ∈ U, HasDerivAt (f j) (f (j + 1) y) y) (hc : ContinuousOn (f N) U) :
    (∀ m j : ℕ, j + m ≤ N → ContDiffOn ℝ m (f j) U) ∧
    (∀ j ≤ N, ∀ y ∈ U, iteratedDeriv j (f 0) y = f j y) := by
  have hdv : ∀ j < N, EqOn (deriv (f j)) (f (j + 1)) U := fun j hj y hy => (hd j hj y hy).deriv
  refine ⟨?_, ?_⟩
  · intro m
    induction m with
    | zero =>
      intro j hj
      rw [Nat.cast_zero, contDiffOn_zero]
      rcases lt_or_eq_of_le (show j ≤ N by omega) with h | h
      · exact fun y hy => (hd j h y hy).continuousAt.continuousWithinAt
      · subst h; exact hc
    | succ m ih =>
      intro j hj
      rw [show ((m + 1 : ℕ) : WithTop ℕ∞) = (m : WithTop ℕ∞) + 1 by push_cast; rfl,
        contDiffOn_succ_iff_deriv_of_isOpen hU]
      refine ⟨fun y hy => (hd j (by omega) y hy).differentiableAt.differentiableWithinAt,
        by simp, (ih (j + 1) (by omega)).congr (hdv j (by omega))⟩
  · intro j
    induction j with
    | zero => intro _ y _; simp
    | succ j ih =>
      intro hj y hy
      have hev : iteratedDeriv j (f 0) =ᶠ[𝓝 y] f j := by
        filter_upwards [hU.mem_nhds hy] with z hz
        exact ih (by omega) z hz
      rw [iteratedDeriv_succ, hev.deriv_eq]; exact hdv j (by omega) hy

lemma Gf_cd5 (H : Hex) (hH : IsHex H) (k : Fin 3) (t : ℝ) :
    ContDiffOn ℝ 5 (fun q => Gf H k q t) (Ioi 0) := by
  have := ch_cd5 (α := eα H k t) (β := eβ H k) (ew_pos H hH k) (edge_det_ne H hH k t)
  exact this.exp

lemma Gw_hasDeriv5 (H : Hex) (hH : IsHex H) (k : Fin 3) {j : ℕ} (hj : j < 5) {p : ℝ}
    (hp : 0 < p) (t : ℝ) : HasDerivAt (fun q => Gw H k j q t) (Gw H k (j + 1) p t) p := by
  have hd := (Gf_cd5 H hH k t).differentiableOn_iteratedDerivWithin (m := j)
    (by exact_mod_cast hj) uS
  have hda := (hd p hp).differentiableAt (isOpen_Ioi.mem_nhds hp)
  have := hda.hasDerivAt
  unfold Gw
  rw [iteratedDerivWithin_succ, derivWithin_of_isOpen isOpen_Ioi hp]
  exact this

lemma Gw_contOn5 (H : Hex) (hH : IsHex H) (k : Fin 3) (t : ℝ) :
    ContinuousOn (fun q => Gw H k 5 q t) (Ioi 0) :=
  (Gf_cd5 H hH k t).continuousOn_iteratedDerivWithin le_rfl uS

lemma Gw_bound5 (H : Hex) (hH : IsHex H) (k : Fin 3) :
    ∃ K, ∀ p ∈ Icc (3/4:ℝ) (5/4), ∀ t ∈ Ioo (0:ℝ) 1, ∀ j ≤ 5,
      |Gw H k j p t| ≤ K * (t * (1 - t)) ^ (-(3/4:ℝ)) := by
  obtain ⟨K, hK⟩ := edge_bound5 H hH k
  exact ⟨K, fun p hp t ht j hj => hK p hp t ht j hj⟩

lemma Gw_meas5 (H : Hex) (hH : IsHex H) (k : Fin 3) :
    ∀ j ≤ 5, ∀ q > 0, AEStronglyMeasurable (fun t => Gw H k j q t)
      (volume.restrict (Ι (0:ℝ) 1)) := by
  intro j
  induction j with
  | zero =>
    intro _ q hq
    refine (Gl_meas H k q).aestronglyMeasurable.congr ?_
    filter_upwards [ae_Ioo01] with t ht
    rw [Gw_zero]; exact (arc_deriv H hH k hq ht).2
  | succ j ih =>
    intro hj q hq
    have hq' : ∀ m : ℕ, 0 < q + 1 / ((m : ℝ) + 1) := fun m => by positivity
    refine aestronglyMeasurable_of_tendsto_ae atTop
      (f := fun (m : ℕ) t => ((m : ℝ) + 1) * (Gw H k j (q + 1 / ((m : ℝ) + 1)) t - Gw H k j q t))
      (fun m => ((ih (by omega) _ (hq' m)).sub (ih (by omega) q hq)).const_mul _) ?_
    refine Filter.Eventually.of_forall fun t => ?_
    have hD := Gw_hasDeriv5 H hH k (j := j) (by omega) hq t
    have hs := hD.tendsto_slope_zero_right
    have hseq : Tendsto (fun m : ℕ => 1 / ((m : ℝ) + 1)) atTop (𝓝[>] 0) := by
      refine tendsto_nhdsWithin_iff.2 ⟨tendsto_one_div_add_atTop_nhds_zero_nat, ?_⟩
      exact Filter.Eventually.of_forall fun m => mem_Ioi.2 (by positivity)
    refine (hs.comp hseq).congr fun m => ?_
    simp only [Function.comp_apply, smul_eq_mul, one_div, inv_inv]

lemma Gw_int5 (H : Hex) (hH : IsHex H) (k : Fin 3) {j : ℕ} (hj : j ≤ 5) {p : ℝ}
    (hp : p ∈ Icc (3/4:ℝ) (5/4)) : IntervalIntegrable (fun t => Gw H k j p t) volume 0 1 := by
  obtain ⟨K, hK⟩ := Gw_bound5 H hH k
  refine (bnd_int K).mono_fun' (Gw_meas5 H hH k j hj p (by linarith [hp.1])) ?_
  filter_upwards [ae_Ioo01] with t ht
  rw [Real.norm_eq_abs]; exact hK p hp t ht j hj

lemma Phi_deriv5 (H : Hex) (hH : IsHex H) (k : Fin 3) {j : ℕ} (hj : j < 5) {p0 : ℝ}
    (hp0 : p0 ∈ Ioo (3/4:ℝ) (5/4)) : HasDerivAt (Phi H k j) (Phi H k (j + 1) p0) p0 := by
  obtain ⟨K, hK⟩ := Gw_bound5 H hH k
  have := intervalIntegral.hasDerivAt_integral_of_dominated_loc_of_deriv_le (μ := volume)
    (a := 0) (b := 1) (F := fun p t => Gw H k j p t) (F' := fun p t => Gw H k (j + 1) p t)
    (x₀ := p0) (s := Ioo (3/4) (5/4)) (bound := fun t => K * (t * (1 - t)) ^ (-(3/4:ℝ)))
    (isOpen_Ioo.mem_nhds hp0)
    (by filter_upwards [isOpen_Ioo.mem_nhds hp0] with q hq
        exact Gw_meas5 H hH k j hj.le q (by linarith [hq.1]))
    (Gw_int5 H hH k hj.le (Ioo_subset_Icc_self hp0))
    (Gw_meas5 H hH k (j + 1) hj p0 (by linarith [hp0.1]))
    (ae_io fun t ht x hx => by
      rw [Real.norm_eq_abs]; exact hK x (Ioo_subset_Icc_self hx) t ht (j + 1) hj)
    (bnd_int K)
    (ae_io fun t _ x hx => Gw_hasDeriv5 H hH k hj (by linarith [hx.1]) t)
  exact this.2

lemma Phi_cont5 (H : Hex) (hH : IsHex H) (k : Fin 3) :
    ContinuousOn (Phi H k 5) (Ioo (3/4:ℝ) (5/4)) := by
  obtain ⟨K, hK⟩ := Gw_bound5 H hH k
  intro p0 hp0
  refine (intervalIntegral.continuousAt_of_dominated_interval (μ := volume)
    (F := fun p t => Gw H k 5 p t) (bound := fun t => K * (t * (1 - t)) ^ (-(3/4:ℝ)))
    ?_ ?_ (bnd_int K) ?_).continuousWithinAt
  · filter_upwards [isOpen_Ioo.mem_nhds hp0] with q hq
    exact Gw_meas5 H hH k 5 le_rfl q (by linarith [hq.1])
  · filter_upwards [isOpen_Ioo.mem_nhds hp0] with q hq
    exact ae_io fun t ht => by
      rw [Real.norm_eq_abs]; exact hK q (Ioo_subset_Icc_self hq) t ht 5 le_rfl
  · exact ae_io fun t _ => ((Gw_contOn5 H hH k t) p0 (by simp; linarith [hp0.1])).continuousAt
      (isOpen_Ioi.mem_nhds (by simp; linarith [hp0.1]))

lemma iter_Gl5 (H : Hex) (hH : IsHex H) (k : Fin 3) {t : ℝ} (ht : t ∈ Ioo (0:ℝ) 1) :
    iteratedDeriv 5 (fun p => Gl H k p t) 1 = Gw H k 5 1 t := by
  have hev : (fun p => Gl H k p t) =ᶠ[𝓝 1] (fun p => Gf H k p t) := by
    filter_upwards [isOpen_Ioi.mem_nhds (show (1:ℝ) ∈ Ioi 0 by norm_num)] with p hp
    exact (arc_deriv H hH k hp ht).2
  rw [hev.iteratedDeriv_eq 5, Gw, iteratedDerivWithin_of_isOpen isOpen_Ioi (by norm_num : (1:ℝ) ∈ Ioi 0)]

end P03.P5

theorem P03.proof_1 : P03.stmt_1 := by
  intro H hH
  set f : ℕ → ℝ → ℝ := fun j p => 2 * ∑ k, S05.Phi H k j p with hf
  have hd : ∀ j < 5, ∀ y ∈ Set.Ioo (3/4:ℝ) (5/4), HasDerivAt (f j) (f (j + 1) y) y :=
    fun j hj y hy => (HasDerivAt.fun_sum fun k _ => P03.P5.Phi_deriv5 H hH k hj hy).const_mul 2
  have hc : ContinuousOn (f 5) (Set.Ioo (3/4:ℝ) (5/4)) :=
    continuousOn_const.mul (continuousOn_finsetSum (f := fun k p => S05.Phi H k 5 p) Finset.univ
      fun k _ => P03.P5.Phi_cont5 H hH k)
  obtain ⟨hcd, hit⟩ := P03.P5.chainN 5 f isOpen_Ioo hd hc
  have hev : S05.Perim H =ᶠ[nhds 1] f 0 := by
    filter_upwards [isOpen_Ioo.mem_nhds S05.one_mem_U] with p hp
    exact S05.perim_eq H hH (Set.Ioo_subset_Icc_self hp)
  have hEe : ∀ k, P03.Ee5 H k = S05.Phi H k 5 1 := fun k =>
    intervalIntegral.integral_congr_ae (S05.ae_io fun t ht => P03.P5.iter_Gl5 H hH k ht)
  refine ⟨((hcd 5 0 (by norm_num)).contDiffAt (isOpen_Ioo.mem_nhds S05.one_mem_U)).congr_of_eventuallyEq hev,
    fun k => ?_, ?_⟩
  · exact (P03.P5.Gw_int5 H hH k (j := 5) le_rfl (Set.Ioo_subset_Icc_self S05.one_mem_U)).congr_ae
      (by filter_upwards [S05.ae_Ioo01] with t ht; exact (P03.P5.iter_Gl5 H hH k ht).symm)
  · rw [hev.iteratedDeriv_eq 5, hit 5 le_rfl 1 S05.one_mem_U]
    simp only [hf, hEe]

namespace P03.B5
open Real MeasureTheory Set Filter Topology

/-- power series on a ball give the Taylor coefficients -/
lemma ps_coef {F : ℝ → ℝ} {c : ℕ → ℝ} {r : ℝ} (hr : 0 < r)
    (h : ∀ t : ℝ, |t| < r → HasSum (fun k => c k * t ^ k) (F t)) :
    AnalyticAt ℝ F 0 ∧ ∀ k : ℕ, iteratedDeriv k F 0 = (k.factorial : ℝ) * c k := by
  set p : FormalMultilinearSeries ℝ ℝ ℝ := FormalMultilinearSeries.ofScalars ℝ c with hp
  have hrad : ENNReal.ofReal r ≤ p.radius := by
    refine ENNReal.le_of_forall_nnreal_lt fun r' hr' => ?_
    have hr'1 : (r' : ℝ) < r := by
      have := (ENNReal.lt_ofReal_iff_toReal_lt (by simp)).1 hr'
      simpa using this
    have hs := (h r' (by rw [NNReal.abs_eq]; exact hr'1)).summable
    refine p.le_radius_of_tendsto (l := 0) ?_
    have := hs.tendsto_atTop_zero.norm
    rw [norm_zero] at this
    refine this.congr fun n => ?_
    rw [hp, FormalMultilinearSeries.ofScalars_norm, norm_mul, norm_pow, Real.norm_eq_abs,
      Real.norm_eq_abs, NNReal.abs_eq]
  have hF : HasFPowerSeriesOnBall F p 0 (ENNReal.ofReal r) := by
    refine ⟨hrad, by simpa using hr, fun {y} hy => ?_⟩
    have hy' : |y| < r := by
      rw [Metric.mem_eball, edist_zero_right, Real.eNorm_eq_ofReal_abs] at hy
      exact (ENNReal.ofReal_lt_ofReal_iff hr).1 hy
    rw [zero_add]
    have := h y hy'
    simpa [hp, FormalMultilinearSeries.ofScalars_apply_eq, mul_comm] using this
  have hA : AnalyticAt ℝ F 0 := hF.analyticAt
  refine ⟨hA, fun k => ?_⟩
  have h1 := hA.hasFPowerSeriesAt
  have h2 := h1.eq_formalMultilinearSeries hF.hasFPowerSeriesAt
  have := congrFun (FormalMultilinearSeries.ofScalars_series_injective ℝ ℝ h2) k
  rw [← this]
  have : (k.factorial : ℝ) ≠ 0 := by positivity
  field_simp

/-- the evaluation vector of stmt_5 -/
noncomputable def vI (N : ℕ) : Fin 9 → ℝ :=
  ![1 / ((N : ℝ) + 1), (S03.hn 1 (N + 1) : ℝ), (S03.hn 2 (N + 1) : ℝ),
    (S03.hn 3 (N + 1) : ℝ), (S03.hn 4 (N + 1) : ℝ), (S03.hn 5 (N + 1) : ℝ),
    X05.z 2, X05.z 3, X05.z 5]

def wt9 (m : Fin 9 →₀ ℕ) : ℕ := 2 * m 6 + 3 * m 7 + 5 * m 8

lemma wt9_add (a b : Fin 9 →₀ ℕ) : wt9 (a + b) = wt9 a + wt9 b := by
  simp [wt9]; ring

def RepI (w : ℕ) (f : ℕ → ℝ) : Prop :=
  ∃ P : MvPolynomial (Fin 9) ℚ, (∀ m ∈ P.support, wt9 m ≤ w) ∧
    ∀ N, f N = MvPolynomial.eval₂ (algebraMap ℚ ℝ) (vI N) P

lemma rep_congr {w : ℕ} {f g : ℕ → ℝ} (h : ∀ N, f N = g N) (hf : RepI w f) : RepI w g := by
  obtain ⟨P, h1, h2⟩ := hf
  exact ⟨P, h1, fun N => by rw [← h N, h2]⟩

lemma rep_mono {w w' : ℕ} {f : ℕ → ℝ} (hf : RepI w f) (h : w ≤ w') : RepI w' f := by
  obtain ⟨P, h1, h2⟩ := hf
  exact ⟨P, fun m hm => (h1 m hm).trans h, h2⟩

lemma rep_add {w : ℕ} {f g : ℕ → ℝ} (hf : RepI w f) (hg : RepI w g) :
    RepI w (fun N => f N + g N) := by
  classical
  obtain ⟨P, h1, h2⟩ := hf
  obtain ⟨Q, h3, h4⟩ := hg
  refine ⟨P + Q, fun m hm => ?_, fun N => by rw [MvPolynomial.eval₂_add, ← h2, ← h4]⟩
  rcases Finset.mem_union.1 (MvPolynomial.support_add hm) with h | h
  · exact h1 m h
  · exact h3 m h

lemma rep_mul {w w' : ℕ} {f g : ℕ → ℝ} (hf : RepI w f) (hg : RepI w' g) :
    RepI (w + w') (fun N => f N * g N) := by
  classical
  obtain ⟨P, h1, h2⟩ := hf
  obtain ⟨Q, h3, h4⟩ := hg
  refine ⟨P * Q, fun m hm => ?_, fun N => by rw [MvPolynomial.eval₂_mul, ← h2, ← h4]⟩
  obtain ⟨a, ha, b, hb, rfl⟩ := Finset.mem_add.1 (MvPolynomial.support_mul P Q hm)
  rw [wt9_add]
  exact Nat.add_le_add (h1 a ha) (h3 b hb)

lemma rep_const (c : ℚ) : RepI 0 (fun _ => (c : ℝ)) := by
  classical
  refine ⟨MvPolynomial.C c, fun m hm => ?_, fun N => by simp⟩
  rw [MvPolynomial.support_C] at hm
  split_ifs at hm
  · simp at hm
  · rw [Finset.mem_singleton] at hm; subst hm; simp [wt9]

lemma rep_X (i : Fin 9) : RepI (wt9 (Finsupp.single i 1)) (fun N => vI N i) := by
  classical
  refine ⟨MvPolynomial.X i, fun m hm => ?_, fun N => by simp⟩
  rw [MvPolynomial.support_X, Finset.mem_singleton] at hm
  rw [hm]

lemma rep_smul {w : ℕ} {f : ℕ → ℝ} (c : ℚ) (hf : RepI w f) : RepI w (fun N => (c : ℝ) * f N) := by
  have := rep_mul (rep_const c) hf
  simpa using this

lemma rep_zero (w : ℕ) : RepI w (fun _ => 0) := by
  have := rep_mono (rep_const 0) (Nat.zero_le w)
  simpa using this

lemma rep_sum {w : ℕ} {ι : Type*} (s : Finset ι) (f : ι → ℕ → ℝ) (hf : ∀ i ∈ s, RepI w (f i)) :
    RepI w (fun N => ∑ i ∈ s, f i N) := by
  classical
  induction s using Finset.induction_on with
  | empty => simpa using rep_zero w
  | insert i s hi ih =>
    simp only [Finset.sum_insert hi]
    exact rep_add (hf i (Finset.mem_insert_self _ _))
      (ih fun j hj => hf j (Finset.mem_insert_of_mem hj))

lemma rep_pow {w : ℕ} {f : ℕ → ℝ} (hf : RepI 0 f) (k : ℕ) : RepI w (fun N => f N ^ k) := by
  induction k with
  | zero => simpa using rep_mono (rep_const 1) (Nat.zero_le w)
  | succ k ih =>
    have := rep_mul ih hf
    simpa [pow_succ] using this

lemma rep_inv (k : ℕ) (w : ℕ) : RepI w (fun N => 1 / ((N : ℝ) + 1) ^ k) := by
  have h0 : RepI 0 (fun N => vI N 0) := by simpa [wt9] using rep_X 0
  have := rep_pow (w := w) h0 k
  refine rep_congr (fun N => ?_) this
  simp [vI]

lemma rep_hn (k : ℕ) (hk1 : 1 ≤ k) (hk5 : k ≤ 5) (w : ℕ) :
    RepI w (fun N => (S03.hn k (N + 1) : ℝ)) := by
  have hX : ∀ i : Fin 9, i.1 < 6 → RepI w (fun N => vI N i) := fun i hi =>
    rep_mono (by
      have h := rep_X i
      have : wt9 (Finsupp.single i 1) = 0 := by
        simp only [wt9, Finsupp.single_apply]
        have h6 : i ≠ 6 := fun h => by rw [h] at hi; simp at hi
        have h7 : i ≠ 7 := fun h => by rw [h] at hi; simp at hi
        have h8 : i ≠ 8 := fun h => by rw [h] at hi; simp at hi
        simp [h6, h7, h8]
      rwa [this] at h) (Nat.zero_le w)
  interval_cases k
  · exact rep_congr (fun N => by simp [vI]) (hX 1 (by decide))
  · exact rep_congr (fun N => by simp [vI]) (hX 2 (by decide))
  · exact rep_congr (fun N => by simp [vI]) (hX 3 (by decide))
  · exact rep_congr (fun N => by simp [vI]) (hX 4 (by decide))
  · exact rep_congr (fun N => by simp [vI]) (hX 5 (by decide))

lemma rep_z2 : RepI 2 (fun _ => X05.z 2) := by
  have h : RepI 2 (fun N => vI N 6) := by simpa [wt9] using rep_X 6
  exact rep_congr (fun N => by simp [vI]) h

lemma rep_z3 : RepI 3 (fun _ => X05.z 3) := by
  have h : RepI 3 (fun N => vI N 7) := by simpa [wt9] using rep_X 7
  exact rep_congr (fun N => by simp [vI]) h

lemma rep_z5 : RepI 5 (fun _ => X05.z 5) := by
  have h : RepI 5 (fun N => vI N 8) := by simpa [wt9] using rep_X 8
  exact rep_congr (fun N => by simp [vI]) h

lemma z4_eq : X05.z 4 = 2 / 5 * X05.z 2 ^ 2 := by
  rw [X05.z_four, X05.z_two]; ring

lemma rep_z4 : RepI 4 (fun _ => X05.z 4) := by
  have := rep_smul (2 / 5) (rep_mul rep_z2 rep_z2)
  refine rep_congr (fun N => ?_) this
  rw [z4_eq]; push_cast; ring

/-- `ℓ(s) = α log s + log(1-s)` -/
noncomputable def lf (α : ℕ) (s : ℝ) : ℝ := (α : ℝ) * Real.log s + Real.log (1 - s)

/-- `g(t) = ∫ s^N e^{t ℓ(s)} ds = B(N+1+αt, 1+t)` -/
noncomputable def gfun (α N : ℕ) (t : ℝ) : ℝ := ∫ s in (0:ℝ)..1, s ^ N * Real.exp (t * lf α s)

/-- `K_m = ∫ s^N ℓ(s)^m ds` -/
noncomputable def Kint (α N m : ℕ) : ℝ := ∫ s in (0:ℝ)..1, s ^ N * lf α s ^ m

lemma lf_bd {α : ℕ} (hα : α ≤ 5) (s : ℝ) : |lf α s| ≤ 6 * S05.ls s := by
  unfold lf S05.ls
  have h1 : (α : ℝ) ≤ 5 := by exact_mod_cast hα
  have h0 : (0:ℝ) ≤ α := Nat.cast_nonneg α
  calc |(α : ℝ) * Real.log s + Real.log (1 - s)| ≤ |(α : ℝ) * Real.log s| + |Real.log (1 - s)| :=
        abs_add_le _ _
    _ = α * |Real.log s| + |Real.log (1 - s)| := by rw [abs_mul, abs_of_nonneg h0]
    _ ≤ 6 * (|Real.log s| + |Real.log (1 - s)|) := by
        nlinarith [abs_nonneg (Real.log s), abs_nonneg (Real.log (1 - s))]

lemma ls_eq {s : ℝ} (hs : s ∈ Ioo (0:ℝ) 1) : S05.ls s = -Real.log (s * (1 - s)) := by
  unfold S05.ls
  have h1 : 0 < 1 - s := by linarith [hs.2]
  rw [abs_of_nonpos (Real.log_nonpos hs.1.le hs.2.le),
    abs_of_nonpos (Real.log_nonpos h1.le (by linarith [hs.1])), Real.log_mul hs.1.ne' h1.ne']
  ring

lemma exp_ls_le {s θ : ℝ} (hs : s ∈ Ioo (0:ℝ) 1) (hθ : θ ≤ 1 / 2) :
    Real.exp (θ * S05.ls s) ≤ (s * (1 - s)) ^ (-(1/2 : ℝ)) := by
  have h1 : 0 < 1 - s := by linarith [hs.2]
  have hp : 0 < s * (1 - s) := mul_pos hs.1 h1
  have hp1 : s * (1 - s) ≤ 1 := by nlinarith [hs.1, hs.2]
  rw [ls_eq hs]
  have : Real.exp (θ * -Real.log (s * (1 - s))) = (s * (1 - s)) ^ (-θ) := by
    rw [Real.rpow_def_of_pos hp]; ring_nf
  rw [this]
  exact Real.rpow_le_rpow_of_exponent_ge hp hp1 (by linarith)

lemma g_series {α : ℕ} (hα : α ≤ 5) (N : ℕ) {t : ℝ} (ht : |t| < 1 / 12) :
    HasSum (fun m : ℕ => Kint α N m / (m.factorial : ℝ) * t ^ m) (gfun α N t) := by
  have key := intervalIntegral.hasSum_integral_of_dominated_convergence
    (μ := MeasureTheory.volume) (a := 0) (b := 1)
    (F := fun (m : ℕ) (s : ℝ) => s ^ N * (t ^ m / (m.factorial : ℝ) * lf α s ^ m))
    (f := fun s => s ^ N * Real.exp (t * lf α s))
    (fun m s => (|t| * (6 * S05.ls s)) ^ m / (m.factorial : ℝ))
    (fun m => (by unfold lf; fun_prop : Measurable (fun s : ℝ =>
      s ^ N * (t ^ m / (m.factorial : ℝ) * lf α s ^ m))).aestronglyMeasurable)
    (fun m => by
      refine Filter.Eventually.of_forall (fun s hs => ?_)
      rw [Set.uIoc_of_le zero_le_one] at hs
      rw [Real.norm_eq_abs, abs_mul, abs_mul, abs_div, abs_pow, abs_pow, abs_pow,
        Nat.abs_cast]
      have h1 : |s| ^ N ≤ 1 := pow_le_one₀ (abs_nonneg s) (by rw [abs_of_pos hs.1]; exact hs.2)
      have h2 : |lf α s| ^ m ≤ (6 * S05.ls s) ^ m := pow_le_pow_left₀ (abs_nonneg _) (lf_bd hα s) m
      have h3 : 0 ≤ |t| ^ m / (m.factorial : ℝ) * |lf α s| ^ m := by positivity
      calc |s| ^ N * (|t| ^ m / (m.factorial : ℝ) * |lf α s| ^ m)
          ≤ 1 * (|t| ^ m / (m.factorial : ℝ) * (6 * S05.ls s) ^ m) := by gcongr
        _ = _ := by rw [mul_pow]; ring)
    (Filter.Eventually.of_forall (fun s _ =>
      (NormedSpace.expSeries_div_hasSum_exp (|t| * (6 * S05.ls s))).summable))
    (by
      have he : (fun s => ∑' m : ℕ, (|t| * (6 * S05.ls s)) ^ m / (m.factorial : ℝ)) =
          fun s => Real.exp (6 * |t| * S05.ls s) := by
        funext s
        rw [(NormedSpace.expSeries_div_hasSum_exp (|t| * (6 * S05.ls s))).tsum_eq,
          ← Real.exp_eq_exp_ℝ]; ring_nf
      rw [he]
      refine IntervalIntegrable.mono_fun' (S05.ii_tt (-(1/2)) (by norm_num))
        (by unfold S05.ls; fun_prop : Measurable (fun s : ℝ =>
          Real.exp (6 * |t| * S05.ls s))).aestronglyMeasurable ?_
      filter_upwards [S05.ae_Ioo01] with s hs
      rw [Real.norm_eq_abs, abs_of_pos (Real.exp_pos _)]
      exact exp_ls_le hs (by linarith))
    (Filter.Eventually.of_forall (fun s _ => by
      have := (NormedSpace.expSeries_div_hasSum_exp (t * lf α s)).mul_left (s ^ N)
      rw [← Real.exp_eq_exp_ℝ] at this
      convert this using 1
      funext m; rw [mul_pow]; ring))
  convert key using 1
  funext m
  unfold Kint
  have : (fun s : ℝ => s ^ N * (t ^ m / (m.factorial : ℝ) * lf α s ^ m)) =
      fun s => (t ^ m / (m.factorial : ℝ)) * (s ^ N * lf α s ^ m) := by funext s; ring
  rw [this, intervalIntegral.integral_const_mul]; ring
  rfl

lemma Gamma_shift {x : ℝ} (hx : 0 < x) (n : ℕ) :
    Real.Gamma (x + n) = Real.Gamma x * ∏ j ∈ Finset.range n, (x + j) := by
  induction n with
  | zero => simp
  | succ n ih =>
    rw [Nat.cast_succ, Finset.prod_range_succ,
      show x + ((n : ℝ) + 1) = (x + n) + 1 by ring,
      Real.Gamma_add_one (by positivity), ih]
    ring

lemma abs_small {α : ℕ} (_hα : α ≤ 5) {t : ℝ} (ht : |t| < 1 / 12) (c : ℝ) (hc0 : 0 ≤ c)
    (hc : c ≤ 6) : |c * t| < 1 / 2 := by
  rw [abs_mul, abs_of_nonneg hc0]
  calc c * |t| ≤ 6 * |t| := by gcongr
    _ < 1 / 2 := by linarith

lemma g_beta {α : ℕ} (hα : α ≤ 5) (N : ℕ) {t : ℝ} (ht : |t| < 1 / 12) :
    gfun α N t = (∏ j ∈ Finset.range N, (1 + (α : ℝ) * t + j)) /
      (∏ j ∈ Finset.range (N + 1), (1 + ((α : ℝ) + 1) * t + j)) *
      X05.A (((α : ℝ) + 1) * t) (((α : ℝ) - 1) * t / 2) := by
  have hα' : (α : ℝ) ≤ 5 := by exact_mod_cast hα
  have hα0 : (0 : ℝ) ≤ α := Nat.cast_nonneg α
  have h1 := abs_small hα ht α hα0 (by linarith)
  have h2 := abs_small hα ht (α + 1) (by linarith) (by linarith)
  have h3 := abs_small hα ht 1 (by norm_num) (by norm_num)
  rw [one_mul] at h3
  obtain ⟨h1a, h1b⟩ := abs_lt.1 h1
  obtain ⟨h2a, h2b⟩ := abs_lt.1 h2
  obtain ⟨h3a, h3b⟩ := abs_lt.1 h3
  set a : ℝ := (1 + α * t) + N with ha
  set b : ℝ := 1 + t with hb
  have hN : (0 : ℝ) ≤ N := Nat.cast_nonneg N
  have ha0 : 0 < a := by rw [ha]; linarith
  have hb0 : 0 < b := by rw [hb]; linarith
  have hint : gfun α N t = ∫ s in (0:ℝ)..1, s ^ (a - 1) * (1 - s) ^ (b - 1) := by
    unfold gfun
    refine intervalIntegral.integral_congr_uIoo fun s hs => ?_
    rw [uIoo_of_le zero_le_one] at hs
    have hs1 : 0 < 1 - s := by linarith [hs.2]
    rw [Real.rpow_def_of_pos hs.1, Real.rpow_def_of_pos hs1, ← Real.rpow_natCast,
      Real.rpow_def_of_pos hs.1, ← Real.exp_add, ← Real.exp_add, lf, ha, hb]
    congr 1; ring
  have hB := X05.beta_real ha0 hb0
  rw [← hint] at hB
  have hab : a + b = (1 + ((α : ℝ) + 1) * t) + ((N + 1 : ℕ) : ℝ) := by
    rw [ha, hb]; push_cast; ring
  rw [hab, Gamma_shift (by linarith) (N + 1), ha, Gamma_shift (by linarith) N] at hB
  have hA := X05.proof_2 (((α : ℝ) + 1) * t) (((α : ℝ) - 1) * t / 2) (by linarith)
    (by rw [abs_lt]; constructor <;> nlinarith)
  rw [show 1 + ((α : ℝ) + 1) * t / 2 + ((α : ℝ) - 1) * t / 2 = 1 + α * t by ring,
    show 1 + ((α : ℝ) + 1) * t / 2 - ((α : ℝ) - 1) * t / 2 = 1 + t by ring] at hA
  rw [hA]
  have gp : 0 < Real.Gamma (1 + ((α : ℝ) + 1) * t) := Real.Gamma_pos_of_pos (by linarith)
  have pp : 0 < ∏ j ∈ Finset.range (N + 1), (1 + ((α : ℝ) + 1) * t + j) :=
    Finset.prod_pos fun j _ => by have := (Nat.cast_nonneg j : (0:ℝ) ≤ j); linarith
  rw [hb] at hB
  have hne : Real.Gamma (1 + ((α : ℝ) + 1) * t) *
      ∏ j ∈ Finset.range (N + 1), (1 + ((α : ℝ) + 1) * t + j) ≠ 0 := by positivity
  rw [eq_div_of_mul_eq hne hB]
  field_simp

/-- `F = log g` near `0` -/
noncomputable def Ff (α N : ℕ) (t : ℝ) : ℝ :=
  ∑ j ∈ Finset.range N, Real.log (1 + (α : ℝ) * t + j) -
    ∑ j ∈ Finset.range (N + 1), Real.log (1 + ((α : ℝ) + 1) * t + j) +
    Real.log (X05.A (((α : ℝ) + 1) * t) (((α : ℝ) - 1) * t / 2))

lemma A_pos {α : ℕ} (hα : α ≤ 5) {t : ℝ} (ht : |t| < 1 / 12) :
    0 < X05.A (((α : ℝ) + 1) * t) (((α : ℝ) - 1) * t / 2) := by
  have hα' : (α : ℝ) ≤ 5 := by exact_mod_cast hα
  have hα0 : (0 : ℝ) ≤ α := Nat.cast_nonneg α
  obtain ⟨h3a, h3b⟩ := abs_lt.1 ht
  refine (X05.proof_3 _ _ ?_ ?_).1
  · rw [abs_lt]; constructor <;> nlinarith
  · rw [abs_lt]; constructor <;> nlinarith

lemma g_exp {α : ℕ} (hα : α ≤ 5) (N : ℕ) {t : ℝ} (ht : |t| < 1 / 12) :
    gfun α N t = Real.exp (Ff α N t) := by
  have hα' : (α : ℝ) ≤ 5 := by exact_mod_cast hα
  have hα0 : (0 : ℝ) ≤ α := Nat.cast_nonneg α
  have h1 := abs_small hα ht α hα0 (by linarith)
  have h2 := abs_small hα ht (α + 1) (by linarith) (by linarith)
  obtain ⟨h1a, h1b⟩ := abs_lt.1 h1
  obtain ⟨h2a, h2b⟩ := abs_lt.1 h2
  rw [g_beta hα N ht, Ff, Real.exp_add, Real.exp_sub, Real.exp_sum, Real.exp_sum,
    Real.exp_log (A_pos hα ht)]
  congr 2
  · exact Finset.prod_congr rfl fun j _ => by
      rw [Real.exp_log]; have := (Nat.cast_nonneg j : (0:ℝ) ≤ j); linarith
  · exact Finset.prod_congr rfl fun j _ => by
      rw [Real.exp_log]; have := (Nat.cast_nonneg j : (0:ℝ) ≤ j); linarith

noncomputable def lc (c : ℝ) (j k : ℕ) : ℝ :=
  if k = 0 then Real.log ((j : ℝ) + 1) else (-1) ^ (k + 1) * c ^ k / ((k : ℝ) * ((j : ℝ) + 1) ^ k)

noncomputable def zc (α k : ℕ) : ℝ :=
  if 2 ≤ k then (-1) ^ k * X05.z k / (k : ℝ) * ((α : ℝ) ^ k + 1 - ((α : ℝ) + 1) ^ k) else 0

noncomputable def phi (α N k : ℕ) : ℝ :=
  ∑ j ∈ Finset.range N, lc α j k - ∑ j ∈ Finset.range (N + 1), lc ((α : ℝ) + 1) j k + zc α k

lemma log_lin_series (c : ℝ) (j : ℕ) {t : ℝ} (h : |c * t| < 1) :
    HasSum (fun k => lc c j k * t ^ k) (Real.log (1 + c * t + j)) := by
  have hj : (0 : ℝ) < (j : ℝ) + 1 := by positivity
  have hj1 : (1 : ℝ) ≤ (j : ℝ) + 1 := by have := (Nat.cast_nonneg j : (0:ℝ) ≤ j); linarith
  set x : ℝ := -(c * t / ((j : ℝ) + 1)) with hx
  have hxa : |x| < 1 := by
    rw [hx, abs_neg, abs_div, abs_of_pos hj, div_lt_one hj]; linarith
  have H := (Real.hasSum_pow_div_log_of_abs_lt_one hxa).neg
  rw [neg_neg] at H
  have hpos : 0 < 1 + c * t + j := by
    have := neg_abs_le (c * t); linarith
  rw [← hasSum_nat_add_iff' 1]
  simp only [Finset.range_one, Finset.sum_singleton, lc, pow_zero, mul_one]
  convert H using 1
  · funext n
    simp only [Nat.add_eq_zero_iff, one_ne_zero, and_false, ite_false]
    rw [hx, neg_pow (c * t / ((j : ℝ) + 1)), div_pow, mul_pow]
    push_cast
    field_simp
    ring
  · rw [if_pos trivial, ← Real.log_div hpos.ne' hj.ne']
    congr 1
    rw [hx]; field_simp; ring

lemma F_series {α : ℕ} (hα : α ≤ 5) (N : ℕ) {t : ℝ} (ht : |t| < 1 / 12) :
    HasSum (fun k => phi α N k * t ^ k) (Ff α N t) := by
  have hα' : (α : ℝ) ≤ 5 := by exact_mod_cast hα
  have hα0 : (0 : ℝ) ≤ α := Nat.cast_nonneg α
  have h1 := abs_small hα ht α hα0 (by linarith)
  have h2 := abs_small hα ht (α + 1) (by linarith) (by linarith)
  have S1 := hasSum_sum (s := Finset.range N) fun j _ =>
    log_lin_series (α : ℝ) j (t := t) (by linarith)
  have S2 := hasSum_sum (s := Finset.range (N + 1)) fun j _ =>
    log_lin_series ((α : ℝ) + 1) j (t := t) (by linarith)
  obtain ⟨h3a, h3b⟩ := abs_lt.1 ht
  have HA0 := (X05.proof_3 (((α : ℝ) + 1) * t) (((α : ℝ) - 1) * t / 2)
    (by rw [abs_lt]; constructor <;> nlinarith) (by rw [abs_lt]; constructor <;> nlinarith)).2
  have HA : HasSum (fun k => zc α k * t ^ k)
      (Real.log (X05.A (((α : ℝ) + 1) * t) (((α : ℝ) - 1) * t / 2))) := by
    rw [← hasSum_nat_add_iff' 2]
    simp only [Finset.sum_range_succ, Finset.sum_range_zero, zc]
    norm_num
    convert HA0 using 1
    funext n
    rw [show ((α : ℝ) + 1) * t / 2 + ((α : ℝ) - 1) * t / 2 = α * t by ring,
      show ((α : ℝ) + 1) * t / 2 - ((α : ℝ) - 1) * t / 2 = t by ring, mul_pow, mul_pow]
    push_cast
    ring
  have := (S1.sub S2).add HA
  convert this using 1
  · funext k
    rw [phi, add_mul, sub_mul, Finset.sum_mul, Finset.sum_mul]
  · rfl

lemma F_coef {α : ℕ} (hα : α ≤ 5) (N : ℕ) :
    AnalyticAt ℝ (Ff α N) 0 ∧
      ∀ k : ℕ, iteratedDeriv k (Ff α N) 0 = (k.factorial : ℝ) * phi α N k :=
  ps_coef (by norm_num : (0:ℝ) < 1 / 12) fun t ht => F_series hα N ht

lemma g_coef {α : ℕ} (hα : α ≤ 5) (N m : ℕ) : iteratedDeriv m (gfun α N) 0 = Kint α N m := by
  rw [(ps_coef (by norm_num : (0:ℝ) < 1 / 12) fun t ht => g_series hα N ht).2 m]
  have : (m.factorial : ℝ) ≠ 0 := by positivity
  field_simp

lemma K_rec {α : ℕ} (hα : α ≤ 5) (N : ℕ) (n : ℕ) (hn : 1 ≤ n) :
    (n : ℝ) * (Kint α N n / (n.factorial : ℝ)) = ∑ k ∈ Finset.Icc 1 n,
      (k : ℝ) * phi α N k * (Kint α N (n - k) / ((n - k).factorial : ℝ)) := by
  obtain ⟨hF, hFc⟩ := F_coef hα N
  have hE : AnalyticAt ℝ (fun x => Real.exp (Ff α N x)) 0 := hF.rexp'
  have hgE : gfun α N =ᶠ[𝓝 0] fun x => Real.exp (Ff α N x) := by
    have : {t : ℝ | |t| < 1 / 12} ∈ 𝓝 (0 : ℝ) := by
      have := Metric.ball_mem_nhds (0 : ℝ) (by norm_num : (0:ℝ) < 1 / 12)
      simpa [Metric.ball, Real.dist_eq] using this
    filter_upwards [this] with t ht
    exact g_exp hα N ht
  have hK : ∀ m, iteratedDeriv m (fun x => Real.exp (Ff α N x)) 0 = Kint α N m := fun m => by
    rw [← hgE.iteratedDeriv_eq, g_coef hα]
  have hw : ∀ k, phi α N k = iteratedDeriv k (Ff α N) 0 / (k.factorial : ℝ) := fun k => by
    rw [hFc k]; have : (k.factorial : ℝ) ≠ 0 := by positivity
    field_simp
  simp only [← hK, hw]
  have hder : deriv (fun x => Real.exp (Ff α N x)) =ᶠ[𝓝 0]
      fun x => deriv (Ff α N) x * Real.exp (Ff α N x) := by
    filter_upwards [hF.eventually_analyticAt] with x hx
    rw [(hx.differentiableAt.hasDerivAt.exp).deriv]; ring
  obtain ⟨m, rfl⟩ : ∃ m, n = m + 1 := ⟨n - 1, by omega⟩
  have h1 : iteratedDeriv (m + 1) (fun x => Real.exp (Ff α N x)) 0 =
      ∑ i ∈ Finset.range (m + 1), (m.choose i : ℝ) * iteratedDeriv (i + 1) (Ff α N) 0 *
        iteratedDeriv (m - i) (fun x => Real.exp (Ff α N x)) 0 := by
    rw [iteratedDeriv_succ', hder.iteratedDeriv_eq]
    have := iteratedDeriv_mul (n := m) hF.deriv.contDiffAt hE.contDiffAt
    rw [show (fun x => deriv (Ff α N) x * Real.exp (Ff α N x)) =
      deriv (Ff α N) * (fun x => Real.exp (Ff α N x)) from rfl, this]
    refine Finset.sum_congr rfl fun i _ => ?_
    rw [← iteratedDeriv_succ']
  have hIcc : Finset.Icc 1 (m + 1) = Finset.Ico 1 (m + 2) := rfl
  rw [hIcc, Finset.sum_Ico_eq_sum_range, show m + 2 - 1 = m + 1 by omega]
  rw [h1, Finset.sum_div, Finset.mul_sum]
  refine Finset.sum_congr rfl fun i hi => ?_
  have hi' : i ≤ m := Nat.lt_succ_iff.1 (Finset.mem_range.1 hi)
  rw [show m + 1 - (1 + i) = m - i by omega, add_comm 1 i, Nat.cast_choose ℝ hi']
  rw [Nat.factorial_succ, Nat.factorial_succ]
  push_cast
  have h1 : (i.factorial : ℝ) ≠ 0 := by positivity
  have h2 : ((m - i).factorial : ℝ) ≠ 0 := by positivity
  have h3 : (m.factorial : ℝ) ≠ 0 := by positivity
  field_simp

lemma hn_cast (k n : ℕ) :
    (S03.hn k n : ℝ) = ∑ j ∈ Finset.range n, 1 / ((j : ℝ) + 1) ^ k := by
  induction n with
  | zero => simp [S03.hn]
  | succ n ih =>
    rw [Finset.sum_range_succ, ← ih]
    simp only [S03.hn]
    rw [Finset.sum_Icc_succ_top (by omega)]
    push_cast
    ring

lemma phi_eq {α N k : ℕ} (hk : 1 ≤ k) :
    phi α N k = (((-1) ^ (k + 1) * (α : ℚ) ^ k / k : ℚ) : ℝ) *
        ((S03.hn k (N + 1) : ℝ) - 1 / ((N : ℝ) + 1) ^ k) -
      (((-1) ^ (k + 1) * ((α : ℚ) + 1) ^ k / k : ℚ) : ℝ) * (S03.hn k (N + 1) : ℝ) + zc α k := by
  have hk0 : k ≠ 0 := by omega
  have e1 : (S03.hn k (N + 1) : ℝ) - 1 / ((N : ℝ) + 1) ^ k =
      ∑ j ∈ Finset.range N, 1 / ((j : ℝ) + 1) ^ k := by
    rw [hn_cast, Finset.sum_range_succ]; ring
  rw [e1, hn_cast, phi, Finset.mul_sum, Finset.mul_sum]
  congr 2
  · refine Finset.sum_congr rfl fun j _ => ?_
    simp only [lc, hk0, if_false]; push_cast
    have : ((j : ℝ) + 1) ≠ 0 := by positivity
    have : (k : ℝ) ≠ 0 := by exact_mod_cast hk0
    field_simp
  · refine Finset.sum_congr rfl fun j _ => ?_
    simp only [lc, hk0, if_false]; push_cast
    have : ((j : ℝ) + 1) ≠ 0 := by positivity
    have : (k : ℝ) ≠ 0 := by exact_mod_cast hk0
    field_simp

lemma zc_rep (α k : ℕ) (hk5 : k ≤ 5) : RepI k (fun _ => zc α k) := by
  by_cases h2 : 2 ≤ k
  · have hz : RepI k (fun _ => X05.z k) := by
      interval_cases k
      · exact rep_z2
      · exact rep_z3
      · exact rep_z4
      · exact rep_z5
    have := rep_smul ((-1) ^ k / k * ((α : ℚ) ^ k + 1 - ((α : ℚ) + 1) ^ k)) hz
    refine rep_congr (fun N => ?_) this
    simp only [zc, if_pos h2]; push_cast; ring
  · refine rep_congr (fun N => ?_) (rep_zero k)
    simp [zc, h2]

lemma phi_rep (α k : ℕ) (hk1 : 1 ≤ k) (hk5 : k ≤ 5) : RepI k (fun N => phi α N k) := by
  have hH := rep_hn k hk1 hk5 0
  have hI := rep_inv k 0
  have h1 := rep_smul (((-1) ^ (k + 1) * (α : ℚ) ^ k / k : ℚ))
    (rep_add hH (rep_smul (-1) hI))
  have h2 := rep_smul (-(((-1) ^ (k + 1) * ((α : ℚ) + 1) ^ k / k : ℚ))) hH
  have h3 := rep_add (rep_mono (rep_add h1 h2) (Nat.zero_le k)) (zc_rep α k hk5)
  refine rep_congr (fun N => ?_) h3
  rw [phi_eq hk1]; push_cast; ring

lemma K0 (α N : ℕ) : Kint α N 0 = 1 / ((N : ℝ) + 1) := by
  simp [Kint, integral_pow]

lemma K_rep {α : ℕ} (hα : α ≤ 5) :
    ∀ n : ℕ, n ≤ 5 → RepI n (fun N => Kint α N n / (n.factorial : ℝ)) := by
  intro n
  induction n using Nat.strong_induction_on with
  | _ n ih =>
    intro hn5
    rcases Nat.eq_zero_or_pos n with h0 | hpos
    · subst h0
      refine rep_congr (fun N => ?_) (rep_inv 1 0)
      rw [K0]; simp
    · have hrec := fun N => K_rec hα N n hpos
      have hsum : RepI n (fun N => ∑ k ∈ Finset.Icc 1 n,
          (k : ℝ) * phi α N k * (Kint α N (n - k) / ((n - k).factorial : ℝ))) := by
        refine rep_sum _ _ fun k hk => ?_
        obtain ⟨hk1, hkn⟩ := Finset.mem_Icc.1 hk
        have := rep_smul (k : ℚ) (rep_mul (phi_rep α k hk1 (by omega))
          (ih (n - k) (by omega) (by omega)))
        refine rep_congr (fun N => ?_) (rep_mono this (by omega))
        push_cast; ring
      have := rep_smul (1 / (n : ℚ)) hsum
      refine rep_congr (fun N => ?_) this
      rw [← hrec N]
      have : (n : ℝ) ≠ 0 := by exact_mod_cast hpos.ne'
      push_cast
      field_simp

lemma K_expand (α N m : ℕ) :
    Kint α N m = ∑ a ∈ Finset.range (m + 1),
      ((α : ℝ) ^ a * (m.choose a : ℝ)) * S05.Iint N a (m - a) := by
  unfold Kint S05.Iint
  have : (fun s : ℝ => s ^ N * lf α s ^ m) = fun s => ∑ a ∈ Finset.range (m + 1),
      ((α : ℝ) ^ a * (m.choose a : ℝ)) * (s ^ N * Real.log s ^ a * Real.log (1 - s) ^ (m - a)) := by
    funext s
    rw [lf, add_pow, Finset.mul_sum]
    refine Finset.sum_congr rfl fun a _ => ?_
    rw [mul_pow]; ring
  rw [this, intervalIntegral.integral_finsetSum]
  · refine Finset.sum_congr rfl fun a _ => ?_
    rw [intervalIntegral.integral_const_mul]
  · intro a _
    exact (S05.ii_g N a (m - a)).const_mul _

lemma I_rep (m : ℕ) (hm : m ≤ 5) (a : ℕ) (ha : a ≤ m) :
    RepI m (fun N => S05.Iint N a (m - a)) := by
  classical
  set V : Matrix (Fin (m + 1)) (Fin (m + 1)) ℚ := Matrix.vandermonde (fun i => ((i : ℕ) : ℚ))
    with hV
  have hdet : V.det ≠ 0 := Matrix.det_vandermonde_ne_zero_iff.2 fun i j h => by
    apply Fin.ext; exact_mod_cast h
  have hWV : V⁻¹ * V = 1 := Matrix.nonsing_inv_mul V (isUnit_iff_ne_zero.2 hdet)
  set x : Fin (m + 1) → ℕ → ℝ := fun j N => (m.choose j : ℝ) * S05.Iint N j (m - j) with hx
  have hK : ∀ i : Fin (m + 1), RepI m (fun N => ∑ j : Fin (m + 1), (V i j : ℝ) * x j N) := by
    intro i
    have h1 := rep_smul ((m.factorial : ℚ)) (K_rep (α := i) (by omega) m hm)
    refine rep_congr (fun N => ?_) h1
    push_cast
    have : (m.factorial : ℝ) ≠ 0 := by positivity
    rw [mul_div_cancel₀ _ this, K_expand, ← Fin.sum_univ_eq_sum_range
      (fun a => ((i : ℕ) : ℝ) ^ a * (m.choose a : ℝ) * S05.Iint N a (m - a)) (m + 1)]
    refine Finset.sum_congr rfl fun j _ => ?_
    rw [hV, Matrix.vandermonde_apply, hx]; push_cast; ring
  have hxa : RepI m (x ⟨a, by omega⟩) := by
    have := rep_sum Finset.univ _ fun i (_ : i ∈ Finset.univ) =>
      rep_smul (V⁻¹ ⟨a, by omega⟩ i) (hK i)
    refine rep_congr (fun N => ?_) this
    simp only [Finset.mul_sum]
    rw [Finset.sum_comm]
    have hcol : ∀ j : Fin (m + 1), ∑ i : Fin (m + 1), ((V⁻¹ ⟨a, by omega⟩ i : ℚ) : ℝ) *
        ((V i j : ℝ) * x j N) = (((V⁻¹ * V) ⟨a, by omega⟩ j : ℚ) : ℝ) * x j N := by
      intro j
      rw [Matrix.mul_apply]; push_cast
      rw [Finset.sum_mul]
      refine Finset.sum_congr rfl fun i _ => by ring
    simp only [hcol, hWV, Matrix.one_apply]
    simp only [apply_ite (fun q : ℚ => (q : ℝ)), Rat.cast_one, Rat.cast_zero, ite_mul, one_mul,
      zero_mul, Finset.sum_ite_eq, Finset.mem_univ, if_true]
  have hc : (m.choose a : ℝ) ≠ 0 := by exact_mod_cast (Nat.choose_pos ha).ne'
  have := rep_smul (1 / (m.choose a : ℚ)) hxa
  refine rep_congr (fun N => ?_) this
  simp only [hx]
  push_cast
  field_simp

lemma stmt11 : P03.stmt_5 := by
  intro a b hab
  have := I_rep (a + b) hab a (by omega)
  rw [show a + b - a = b by omega] at this
  obtain ⟨P, h1, h2⟩ := this
  exact ⟨P, fun m hm => h1 m hm, fun N => h2 N⟩

end P03.B5

theorem P03.proof_5 : P03.stmt_5 := P03.B5.stmt11
