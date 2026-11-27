import Mathlib
import Zeta3.Wave10.S08_FieldRuleMain
/-!
# P02 -- `M_3` at other angles: the symbol of the weight-four part and its `(2,2)` cobracket

Target id: P02.
Paper labels: `prop:M3angle` (paper/ieee/main.tex, Proposition "$M_3$ at other angles", lines
1546-1553) and the sentence before it (line 1544); the sentence after it on one depth-two term
(line 1555).
Supplement: paper/ieee/proofs/sec_g8_M3alpha.tex (statement and table of the 56 terms, lines 5-41;
"What is proved, precisely", 43-46; "What the words mean", 48-72; Lemma `g8:prop:R3`, 76-82 with
machine step M5, 84; the table of 24 pure tensors, 90-124; Lemma `g8:lem:vanish`, 132-138; Lemma
`g8:prop:one`, 145-153 and the ratio remark, 155); sec_g8_M2alpha.tex (Lemma `g8:prop:crit`,
76-91, with `g8:lem:hcoef`, 62-68). Data: paper/ieee/proofs/scripts/g8_terms_k3_part1.json (the 74
terms of `R_3`) and g8_terms_k3_partz2.json (`Ψ`: 14 terms and the constant `1/2`).

Quoted verbatim (main.tex 1544): "The same method gives $M_3(\alpha)$ at every angle in the form
$-\zeta(3)\Phi(c)+\zeta(2)\Psi(c)+R_3(c)$, with $R_3$ a rational combination of 74 terms
$c^eH_w(c)$ (34 distinct harmonic polylogarithms, of weights two to four) (proved for every Taylor
coefficient in our notes)" ... "The $(2,2)$ part of a weight-four symbol lies in $\Lambda^2B_2$,
where $B_2$ is the space of weight-two symbols modulo products; for $x_1\otimes x_2\otimes
x_3\otimes x_4$ it is $(x_1\wedge x_2)\wedge(x_3\wedge x_4)$. It vanishes on products and on
$\Li_4$, so it vanishes on every combination of classical polylogarithms."

Quoted verbatim (main.tex 1546-1553):

  \begin{proposition}[$M_3$ at other angles$^\dagger$]\label{prop:M3angle}
  Let $\Omega(c)$ be the weight-four part of $R_3(c)$ (56 of the 74 terms), whose letters are
  $[c]$, $[1-c]$ and $[1+c]$. Its $(2,2)$ part (of the function, modulo constants), with each
  symbol written with the outermost integration letter last, is
  \begin{multline*}
  \delta_{2,2}(\Omega)=\Bigl(1-\frac1c\Bigr)([c]\wedge[1-c])\wedge([1-c]\wedge[1+c])\\
  -\Bigl(1+\frac1c\Bigr)([c]\wedge[1+c])\wedge([1-c]\wedge[1+c]),
  \end{multline*}
  and it is not zero.
  \end{proposition}

Quoted verbatim (main.tex 1555): "One depth-two term is not enough: for $r(c)\,G$, with $G$ a
function whose symbol has constant coefficients, the two components of $\delta_{2,2}$ above are in a
constant ratio, but for $\Omega$ their ratio is $(1-c)/(1+c)$ (proved at the level of symbols;
refereed)."

Quoted verbatim (supplement 44, "What is proved, precisely"): "Take the displayed expression
$\Omega$ of $56$ terms. Its constant-free symbol (letters modulo constants, Section~\ref{...}), with
the coefficients in $\QQ(c)$ treated as scalars and not differentiated, has the $(2,2)$ cobracket
displayed, and this is not zero."

Quoted verbatim (supplement 55): "The three elements $\mathsf x=[c]$, $\mathsf y=[1-c]$,
$\mathsf z=[1+c]$ form a basis of $V$ (a relation $c^a(1-c)^b(1+c)^d=$ const would force
$a=b=d=0$ by the zeros at $0,1,-1$)."
(supplement 60): "\mathcal S(H_w)=g(a_n)\otimes\dots\otimes g(a_1)\in T^n(V),\qquad g(0)=\mathsf x,\
g(1)=-\mathsf y,\ g(-1)=\mathsf z ." (62) "we put $\mathcal S(F)=\sum_ir_i(c)\,\mathcal
S(H_{w_i})\in T^n(V)\otimes\QQ(c)$: the coefficients $r_i(c)$ are scalars and are not
differentiated." (67) "\delta_{2,2}:T^4(V)\otimes\QQ(c)\to\Lambda^2(\Lambda^2V)\otimes\QQ(c),\qquad
x_1\otimes x_2\otimes x_3\otimes x_4\mapsto(x_1\wedge x_2)\wedge(x_3\wedge x_4)." (69) "The space
$\Lambda^2(\Lambda^2V)$ has the basis $p\wedge q$, $p\wedge r$, $q\wedge r$. ... ``Not zero''
means that these two coefficients are nonzero elements of $\QQ(c)$."
(120-122) "Adding the $24$ rows with their coefficients, the components along $p\wedge q$ cancel,
and \delta_{2,2}(\mathcal S(\Omega))=\Bigl(1-\frac1c\Bigr)p\wedge r-\Bigl(1+\frac1c\Bigr)q\wedge r".

Quoted verbatim (`g8:lem:vanish`, 133): "Let $V'$ be any $\QQ$-vector space. The map
$\delta_{2,2}:T^4(V')\to\Lambda^2(\Lambda^2V')$, $x_1\otimes x_2\otimes x_3\otimes x_4\mapsto(x_1\wedge
x_2)\wedge(x_3\wedge x_4)$, vanishes on (a) every shuffle product $s\sqcup\!\sqcup t$ of tensors of
degrees $(1,3)$, $(2,2)$ or $(3,1)$, and (b) on $-a\otimes b\otimes b\otimes b$ for every $a,b\in
V'$ (the symbol of $\Li_4(r)$ with $a=[1-r]$, $b=[r]$)."

Quoted verbatim (`g8:prop:one`, 146): "Let $U_0$ be any $\QQ$-vector space, $\pi:V\oplus U_0\to V$
the projection along $U_0$, $r\in\QQ(c)^\times$, and $\tau\in T^4(V\oplus U_0)$ any tensor with
\emph{rational} coefficients. Then $\delta_{2,2}\bigl(r(c)\,(\pi^{\otimes4}\tau)\bigr)\ne
\delta_{2,2}(\mathcal S(\Omega))$." (proof, 152) "$\delta_{2,2}(\pi^{\otimes4}\tau)\in
\Lambda^2(\Lambda^2V)$ has rational coefficients $a,b,d$ on $p\wedge q$, $p\wedge r$, $q\wedge r$
... $b/d=\frac{1-1/c}{-(1+1/c)}=\frac{1-c}{1+c}$, which is not constant."

Quoted verbatim (`g8:prop:R3`, 77): "For $0<c<1$, $M_3(\alpha)=R_3(c)-\zeta(3)\Phi(c)+
\zeta(2)\Psi(c)$, where $R_3$ is the combination of $74$ terms of \script{g8_terms_k3_part1.json}
and $\Psi$ the combination of \script{g8_terms_k3_partz2.json} ($14$ terms $q\,c^eH_w$ and the
constant $\frac12$)." (setting, line 8) "$R_3(c)=\sum_{n\ge1}r'_nc^{2n}$, $r'_n=[1]P_{3,n}$", "
$\Psi(c)=\sum_{n\ge1}q'_nc^{2n}$ ... $q'_n=[\gvs2]P_{3,n}$".

Quoted verbatim (`g8:prop:crit`, sec_g8_M2alpha 82-86): "(i) $L_0=T_z$ and $L_1=0$ in
$\mathfrak D$; (ii) the coefficients of $c^N$, for all $N\le1$, in the Laurent series of
$L-\sum_{n\ge1}\rho_nc^{2n}$ ... vanish ... Then $L(c)=\sum_{n\ge1}p_nc^{2n}$ for $0<c<1$."
(proof, 90) "$[c^{2n}]c^eH_w=[x^{2n-e}]H_w=\gsub{h_{w,-e}}{n}$ and
$[c^{2n+1}]c^eH_w=\gsub{h_{w,1-e}}{n}$".

## Formalization choices

* Data. A term `q c^e H_w` is a triple `(q, e, w) : ℚ × ℤ × List ℤ`, the word `w = (a_1,…,a_n)`
  written left to right as in `H_{a_1,…,a_n}` (outermost integration letter `a_1`); `R3terms` is
  the 74 terms of g8_terms_k3_part1.json in file order, `PsiTerms` the 14 terms of
  g8_terms_k3_partz2.json (its constant `1/2` is added separately). `Omega` is the sublist of the
  terms with `|w| = 4`. The 56 entries agree with the table of the supplement (lines 10-29).
* Coefficients. `K = RatFunc ℚ = ℚ(c)`, `c = RatFunc.X`, `c^e` the integer power in the field.
* The letter space `V` (basis `x = [c]`, `y = [1-c]`, `z = [1+c]`) is coordinatized by `Fin 3`
  (`0 ↦ x`, `1 ↦ y`, `2 ↦ z`); `V ⊗ ℚ(c)` is `Fin 3 → K`, `ev i` the basis vector. `stmt_2` is
  the supplement's justification that these three classes are linearly independent modulo
  constants. `T⁴(V) ⊗ ℚ(c)` is `T4 = (Fin 4 → Fin 3) → K` (coordinates on the basis tensors
  `ev (t 0) ⊗ ⋯ ⊗ ev (t 3)`), and `ptens v` is the pure tensor `v 0 ⊗ ⋯ ⊗ v 3`.
* `𝒮(H_w) = g(a_4) ⊗ g(a_3) ⊗ g(a_2) ⊗ g(a_1)` for `|w| = 4` is `symbH w = ptens (i ↦ g(w.reverse[i]))`
  (slot `0` carries the innermost letter `a_4`, slot `3` the outermost `a_1`: "outermost integration
  letter last"); `g 0 = x`, `g 1 = -y`, `g (-1) = z`. `SOmega = Σ q c^e 𝒮(H_w)` over `Omega`, the
  coefficients as scalars.
* Exterior powers, honestly and in coordinates: for a free module `ι → K`, `Λ²` is modelled by
  antisymmetric matrices, `u ∧ v ↦ (i,j) ↦ u_i v_j - u_j v_i` (`wedge`; this is the injective map
  `u ∧ v ↦ u⊗v - v⊗u` into `V ⊗ V`), and `Λ²(Λ²)` the same way one level up (`wedge22`, indices in
  `ι × ι`). A coordinate model of this kind is faithful (the image of `Λ²F` in `F ⊗ F` is a direct
  summand for a free `F`), so an element is zero iff its model is zero. `dpure [x₁,x₂,x₃,x₄] =
  (x₁ ∧ x₂) ∧ (x₃ ∧ x₄)`, and `delta22` is its `K`-linear extension from the basis tensors of `T4`
  (`δ_{2,2} = Σ_t S(t) (e_{t0} ∧ e_{t1}) ∧ (e_{t2} ∧ e_{t3})`). The basis of `Λ²V` is
  `p = x∧y`, `q = x∧z`, `r = y∧z` (`wp`, `wq`, `wr`); the coordinate of an element of
  `Λ²(Λ²V)` on `p∧q`, `p∧r`, `q∧r` is its value at `((0,1),(0,2))`, `((0,1),(1,2))`,
  `((0,2),(1,2))`.
* `g8:lem:vanish` is stated for every index type `ι` and every commutative ring `K` (the paper's
  `ℚ`-vector space `V'` is `ι →₀ ℚ` for a basis `ι`, a submodule of `ι → ℚ`; nothing is lost). The
  shuffle product of two pure tensors `s = s₁⊗⋯`, `t = t₁⊗⋯` is the sum of the pure tensors of the
  shuffled lists (`shuf`), so `δ_{2,2}(s ⧢ t) = Σ_{shuffles} dpure`. By bilinearity this is the
  lemma for all tensors (the reduction the supplement makes in the first line of its proof).
* `g8:prop:one`: `U₀` has a basis `ι`, so `V ⊕ U₀` has the basis `Fin 3 ⊕ ι`; a tensor with
  rational coefficients is `τ : (Fin 4 → Fin 3 ⊕ ι) → ℚ` (its coordinates), and `π^{⊗4}τ` has the
  coordinate `τ (inl ∘ t)` at the basis tensor `t` (`π` kills the basis letters of `U₀`).
* PART 2. `R3fun c = Σ q c^e H_w(c)` over the 74 terms and `PsiFun c = 1/2 + Σ q c^e H_w(c)` over
  the 14 terms, with `S03.H` the harmonic polylogarithm of S03. `[x^N] H_w` is given by the
  dictionary `g6:lem:dict` (`S03.stmt_5`, proved): `hcoef w N = 𝔠_w ε^N N^{-κ} Z_{v_w}(N-1)` for
  `N ≥ 1` and `0` for `N ≤ 0`, with `(ε, κ, 𝔠_w, v_w) = S03.dd w`. Then the coefficient of `c^N`
  in the Laurent series of `Σ q c^e H_w` is `coefTerms L N = Σ q · hcoef w (N - e)`. The finite
  identities of `g8:prop:crit` (i) are identities in `𝔇` whose values at every `n` are these
  coefficients; they are stated here through their values, for every `n ≥ 1` ("proved for every
  Taylor coefficient"), and (ii) as the vanishing for every `N ≤ 1`.
* `r'_n = [1]P_{3,n}` is `V05.P 3 n 0 0` and `q'_n = [𝔷2]P_{3,n}` is `V05.P 3 n 1 0 - V05.P 3 n 0 0`
  (`V05.stmt_11`, proved: `P_{3,n}` is affine in `𝔷2, 𝔷3` with these coefficients, and
  `[𝔷3]P_{3,n} = -1/(4n²-1)`; that part is not restated). `M_3(α)` is `V05.Mk 3 α`,
  `c = D04.c α = cos(α/2)`, `Φ = V05.Phi` (`Φ(c) = Σ c^{2n}/(4n²-1)` by `V05.stmt_1`, proved),
  `ζ(2) = V05.z2 = π²/6`, `ζ(3) = V05.z3 = (riemannZeta 3).re` (real since `riemannZeta 3 =
  Σ n^{-3}` by `zeta_nat_eq_tsum_of_gt_one`).

## Statements

* stmt_1: the data: 74 terms, 2, 16, 56 of weights 2, 3, 4, 34 distinct words, all regular over
  `{0,1,-1}`, exponents in `[-2, 1]`; `Ψ`: 14 terms of weights 1, 2, exponents in `[-2, 1]`.
* stmt_2: `x, y, z` are a basis modulo constants: `c^a(1-c)^b(1+c)^d` constant forces `a=b=d=0`.
* stmt_3: Lemma `g8:lem:vanish` (a) and (b).
* stmt_4: `𝒮(Ω)` is the sum of the 24 pure tensors of the table with the printed coefficients.
* stmt_5: the `δ_{2,2}` column of the table: the image of each of the 24 basis tensors.
* stmt_6: the three coordinates of `δ_{2,2}(𝒮(Ω))`: `0` on `p∧q`, `1 - 1/c` on `p∧r`,
  `-(1 + 1/c)` on `q∧r`.
* stmt_7: the ratio for `Ω` is `(1-c)/(1+c)` (`A(1+c) = B(1-c)`), not a rational constant.
* stmt_8: proof of `g8:prop:one`: for `r ∈ ℚ(c)` and `τ` with rational coefficients,
  `δ_{2,2}(r π^{⊗4}τ) = r (a p∧q + b p∧r + d q∧r)` with `a, b, d ∈ ℚ` (constant ratio).
* stmt_9: Lemma `g8:prop:one`.
* stmt_10: Proposition `prop:M3angle` as proved ("What is proved, precisely"): the displayed
  formula for `δ_{2,2}(𝒮(Ω))`, both coefficients nonzero in `ℚ(c)`, and `δ_{2,2}(𝒮(Ω)) ≠ 0`.
* stmt_11: M5 for `R_3`: `[c^{2n}] = r'_n`, `[c^{2n+1}] = 0` for `n ≥ 1`, `[c^N] = 0` for `N ≤ 1`.
* stmt_12: M5 for `Ψ` (with its constant `1/2`): `[c^{2n}] = q'_n`, odd and low orders `0`.
* stmt_13: `g8:prop:crit` applied: `R_3(c) = Σ r'_n c^{2n}` and `Ψ(c) = Σ q'_n c^{2n}`, `0 < c < 1`.
* stmt_14: Lemma `g8:prop:R3`: `M_3(α) = R_3(c) - ζ(3)Φ(c) + ζ(2)Ψ(c)` for `0 < α < π`.

## Numerical and exact checks before freezing

A Python script (exact rationals, SymPy, mpmath; written from these definitions) checked: the
data counts of `stmt_1`; that `𝒮(Ω)` computed from the 56 terms has exactly the 24 tensors of the
table, typed from the TeX, with equal coefficients (`stmt_4`); the 24 images in the antisymmetric
model (`stmt_5`); `δ_{2,2}(𝒮(Ω))` (`stmt_6`, `stmt_10`), in agreement with the output of the
supplement's g8_M3_delta22.py (run read-only); `hcoef` against an independent series recursion of
`H_w` (`H_{(0,w)} : h_N/N`, `H_{(1,w)} : Σ_{m<N} h_m / N`, `H_{(-1,w)} : Σ_{m<N} (-1)^{N-1-m}h_m/N`)
for all 34 + 7 words, `N ≤ 30`; `stmt_11` and `stmt_12` exactly for `n ≤ 20` and all `N ≤ 1`, with
`P_{3,n}` computed from the definition `V05.P` (`ecoef`, `rho`, `Dsharp`), and
`[𝔷3]P_{3,n} = -1/(4n²-1)`; `stmt_14` at `α = 2π/3` and `α = 1.8` (45 digits, `H_w` summed to
700 terms) against `Σ_{n<120} P_{3,n}(ζ(2),ζ(3)) c^{2n}` (`V05.stmt_8`): differences
`2.9·10⁻⁴⁷` and `1.7·10⁻⁴⁶`; values `0.010717362932289465418409576291224896` and
`0.017173835093232526245372316031957722`, which agree with the supplement's table (line 248-249).

## NOT COVERED

* Lemma `g8:prop:two` and the sentence "Two terms match it" (the two-variable symbol of
  `Li_{2,2}`, Lemma `g8:lem:diff`, the specialisation to `x = 1/2`), and the four residual
  components with the constant letter `[2]` kept.
* "of the function": that the symbol depends only on the function (Section `g8b:sec:M3func`,
  Theorem G8B(b)), and the consequence "`M_3(α)` is not classical as a function of the angle"; the
  identification of `Λ²B₂` with `Λ²(Λ²V)` and Goncharov's criterion (quoted in the supplement,
  not proved there); the passage from shuffle products of symbols to products of functions.
* The square corner specialisation (`Λ²(Λ²)` of a plane is `0`), the hexagon corner, the corner
  `c₀ = 1/3`, the lists of corners and the period-conjecture remarks.
* The construction of the elements `L_0, L_1, T_z ∈ 𝔇` of `g8:prop:crit` and the identities
  "in `𝔇`" themselves: only their values at every `n` (which is what the conclusion uses) are
  stated, in `stmt_11` and `stmt_12`.
* The remark "at `c = 1/√2` it gives `W_4`".
-/
-- BEGIN STATEMENTS
namespace P02

open Real

/-- The 74 terms `(q, e, w)` of `R_3 = Σ q c^e H_w` (g8_terms_k3_part1.json, in file order). -/
def R3terms : List (ℚ × ℤ × List ℤ) := [
  ((-1 : ℚ), (0 : ℤ), [0, 1]),
  ((1 : ℚ), (0 : ℤ), [0, -1]),
  ((1 : ℚ), (0 : ℤ), [0, 0, 1]),
  ((-1 : ℚ), (0 : ℤ), [0, 0, -1]),
  ((1 : ℚ), (0 : ℤ), [0, 1, -1]),
  ((1 : ℚ), (0 : ℤ), [0, -1, 1]),
  ((1 : ℚ), (-2 : ℤ), [1, 0, 1]),
  ((1/2 : ℚ), (-1 : ℤ), [1, 0, 1]),
  ((-3/2 : ℚ), (0 : ℤ), [1, 0, 1]),
  ((-1 : ℚ), (-2 : ℤ), [1, 0, -1]),
  ((-1/2 : ℚ), (-1 : ℤ), [1, 0, -1]),
  ((3/2 : ℚ), (0 : ℤ), [1, 0, -1]),
  ((-1 : ℚ), (-2 : ℤ), [-1, 0, 1]),
  ((1/2 : ℚ), (-1 : ℤ), [-1, 0, 1]),
  ((3/2 : ℚ), (0 : ℤ), [-1, 0, 1]),
  ((1 : ℚ), (-2 : ℤ), [-1, 0, -1]),
  ((-1/2 : ℚ), (-1 : ℤ), [-1, 0, -1]),
  ((-3/2 : ℚ), (0 : ℤ), [-1, 0, -1]),
  ((-1 : ℚ), (0 : ℤ), [0, 0, 0, 1]),
  ((1 : ℚ), (1 : ℤ), [0, 0, 0, 1]),
  ((1 : ℚ), (0 : ℤ), [0, 0, 0, -1]),
  ((1 : ℚ), (1 : ℤ), [0, 0, 0, -1]),
  ((-1 : ℚ), (0 : ℤ), [0, 0, 1, -1]),
  ((1 : ℚ), (1 : ℤ), [0, 0, 1, -1]),
  ((-1 : ℚ), (0 : ℤ), [0, 0, -1, 1]),
  ((-1 : ℚ), (1 : ℤ), [0, 0, -1, 1]),
  ((-1/2 : ℚ), (-1 : ℤ), [0, 1, 0, 1]),
  ((1/2 : ℚ), (1 : ℤ), [0, 1, 0, 1]),
  ((1/2 : ℚ), (-1 : ℤ), [0, 1, 0, -1]),
  ((-1 : ℚ), (0 : ℤ), [0, 1, 0, -1]),
  ((1/2 : ℚ), (1 : ℤ), [0, 1, 0, -1]),
  ((1 : ℚ), (0 : ℤ), [0, 1, -1, 1]),
  ((-1 : ℚ), (1 : ℤ), [0, 1, -1, 1]),
  ((-1/2 : ℚ), (-1 : ℤ), [0, -1, 0, 1]),
  ((-1 : ℚ), (0 : ℤ), [0, -1, 0, 1]),
  ((-1/2 : ℚ), (1 : ℤ), [0, -1, 0, 1]),
  ((1/2 : ℚ), (-1 : ℤ), [0, -1, 0, -1]),
  ((-1/2 : ℚ), (1 : ℤ), [0, -1, 0, -1]),
  ((-1 : ℚ), (0 : ℤ), [0, -1, 1, -1]),
  ((-1 : ℚ), (1 : ℤ), [0, -1, 1, -1]),
  ((-1/2 : ℚ), (-1 : ℤ), [1, 0, 0, 1]),
  ((1/2 : ℚ), (1 : ℤ), [1, 0, 0, 1]),
  ((1/2 : ℚ), (-1 : ℤ), [1, 0, 0, -1]),
  ((-1 : ℚ), (0 : ℤ), [1, 0, 0, -1]),
  ((1/2 : ℚ), (1 : ℤ), [1, 0, 0, -1]),
  ((-1/2 : ℚ), (-1 : ℤ), [1, 0, 1, -1]),
  ((1/2 : ℚ), (1 : ℤ), [1, 0, 1, -1]),
  ((-1/2 : ℚ), (-1 : ℤ), [1, 0, -1, 1]),
  ((1 : ℚ), (0 : ℤ), [1, 0, -1, 1]),
  ((-1/2 : ℚ), (1 : ℤ), [1, 0, -1, 1]),
  ((-1/2 : ℚ), (-1 : ℤ), [1, -1, 0, 1]),
  ((1 : ℚ), (0 : ℤ), [1, -1, 0, 1]),
  ((-1/2 : ℚ), (1 : ℤ), [1, -1, 0, 1]),
  ((1/2 : ℚ), (-1 : ℤ), [1, -1, 0, -1]),
  ((-1/2 : ℚ), (1 : ℤ), [1, -1, 0, -1]),
  ((1 : ℚ), (0 : ℤ), [1, -1, 1, -1]),
  ((-1 : ℚ), (1 : ℤ), [1, -1, 1, -1]),
  ((-1/2 : ℚ), (-1 : ℤ), [-1, 0, 0, 1]),
  ((-1 : ℚ), (0 : ℤ), [-1, 0, 0, 1]),
  ((-1/2 : ℚ), (1 : ℤ), [-1, 0, 0, 1]),
  ((1/2 : ℚ), (-1 : ℤ), [-1, 0, 0, -1]),
  ((-1/2 : ℚ), (1 : ℤ), [-1, 0, 0, -1]),
  ((-1/2 : ℚ), (-1 : ℤ), [-1, 0, 1, -1]),
  ((-1 : ℚ), (0 : ℤ), [-1, 0, 1, -1]),
  ((-1/2 : ℚ), (1 : ℤ), [-1, 0, 1, -1]),
  ((-1/2 : ℚ), (-1 : ℤ), [-1, 0, -1, 1]),
  ((1/2 : ℚ), (1 : ℤ), [-1, 0, -1, 1]),
  ((1/2 : ℚ), (-1 : ℤ), [-1, 1, 0, 1]),
  ((-1/2 : ℚ), (1 : ℤ), [-1, 1, 0, 1]),
  ((-1/2 : ℚ), (-1 : ℤ), [-1, 1, 0, -1]),
  ((-1 : ℚ), (0 : ℤ), [-1, 1, 0, -1]),
  ((-1/2 : ℚ), (1 : ℤ), [-1, 1, 0, -1]),
  ((1 : ℚ), (0 : ℤ), [-1, 1, -1, 1]),
  ((1 : ℚ), (1 : ℤ), [-1, 1, -1, 1])]

/-- The 14 terms `(q, e, w)` of `Ψ - 1/2 = Σ q c^e H_w` (g8_terms_k3_partz2.json; its constant
`1/2` is the first entry of that file and is added in `PsiFun`). -/
def PsiTerms : List (ℚ × ℤ × List ℤ) := [
  ((-1/2 : ℚ), (-2 : ℤ), [1]),
  ((-1/4 : ℚ), (-1 : ℤ), [1]),
  ((3/4 : ℚ), (0 : ℤ), [1]),
  ((1/2 : ℚ), (-2 : ℤ), [-1]),
  ((-1/4 : ℚ), (-1 : ℤ), [-1]),
  ((-3/4 : ℚ), (0 : ℤ), [-1]),
  ((1/4 : ℚ), (-1 : ℤ), [0, 1]),
  ((-1/4 : ℚ), (1 : ℤ), [0, 1]),
  ((1/4 : ℚ), (-1 : ℤ), [0, -1]),
  ((-1/4 : ℚ), (1 : ℤ), [0, -1]),
  ((1/4 : ℚ), (-1 : ℤ), [1, -1]),
  ((-1/4 : ℚ), (1 : ℤ), [1, -1]),
  ((-1/4 : ℚ), (-1 : ℤ), [-1, 1]),
  ((1/4 : ℚ), (1 : ℤ), [-1, 1])]

/-- `Ω`: the terms of `R_3` with `|w| = 4` (56 of the 74). -/
def Omega : List (ℚ × ℤ × List ℤ) := R3terms.filter fun τ => τ.2.2.length = 4

/-- The field of coefficients `ℚ(c)`. -/
abbrev K := RatFunc ℚ

/-- The variable `c ∈ ℚ(c)`. -/
noncomputable def cK : K := RatFunc.X

/-- Basis vectors of `V ⊗ ℚ(c) = Fin 3 → K`: `ev 0 = x = [c]`, `ev 1 = y = [1-c]`,
`ev 2 = z = [1+c]`. -/
noncomputable def ev (i : Fin 3) : Fin 3 → K := Pi.single i 1

/-- `g(0) = x`, `g(1) = -y`, `g(-1) = z` (any other integer gets `z`; only letters `0, ±1`
occur). -/
noncomputable def g (a : ℤ) : Fin 3 → K :=
  if a = 0 then ev 0 else if a = 1 then -ev 1 else ev 2

/-- `T⁴(V) ⊗ ℚ(c)`: coordinates on the basis tensors `ev (t 0) ⊗ ⋯ ⊗ ev (t 3)`. -/
abbrev T4 := (Fin 4 → Fin 3) → K

/-- The pure tensor `v 0 ⊗ v 1 ⊗ v 2 ⊗ v 3`. -/
noncomputable def ptens (v : Fin 4 → Fin 3 → K) : T4 := fun t => ∏ i, v i (t i)

/-- The basis tensor at `l`. -/
noncomputable def btens (l : Fin 4 → Fin 3) : T4 := fun t => if t = l then 1 else 0

/-- `𝒮(H_w) = g(a_4) ⊗ g(a_3) ⊗ g(a_2) ⊗ g(a_1)` for `w = (a_1, a_2, a_3, a_4)`. -/
noncomputable def symbH (w : List ℤ) : T4 := ptens fun i => g (w.reverse.getD i 0)

/-- `𝒮(Ω) = Σ_{(q,e,w) ∈ Ω} q c^e 𝒮(H_w)`, the coefficients as scalars. -/
noncomputable def SOmega : T4 :=
  fun t => (Omega.map fun τ => (τ.1 : K) * cK ^ τ.2.1 * symbH τ.2.2 t).sum

/-- `u ∧ v` in `Λ²`, as the antisymmetric matrix `u_i v_j - u_j v_i`. -/
def wedge {ι R : Type} [CommRing R] (u v : ι → R) : ι → ι → R :=
  fun i j => u i * v j - u j * v i

/-- `P ∧ Q` in `Λ²(Λ²)`, the same one level up (indices in `ι × ι`). -/
def wedge22 {ι R : Type} [CommRing R] (P Q : ι → ι → R) : ι × ι → ι × ι → R :=
  fun a b => P a.1 a.2 * Q b.1 b.2 - P b.1 b.2 * Q a.1 a.2

/-- `δ_{2,2}` of a pure tensor `x₁ ⊗ x₂ ⊗ x₃ ⊗ x₄`: `(x₁ ∧ x₂) ∧ (x₃ ∧ x₄)` (`0` for lists of
another length, which do not occur). -/
def dpure {ι R : Type} [CommRing R] : List (ι → R) → (ι × ι → ι × ι → R)
  | [x₁, x₂, x₃, x₄] => wedge22 (wedge x₁ x₂) (wedge x₃ x₄)
  | _ => 0

/-- `δ_{2,2} : T⁴(V) ⊗ ℚ(c) → Λ²(Λ²V) ⊗ ℚ(c)`, the `ℚ(c)`-linear extension from the basis
tensors. -/
noncomputable def delta22 (S : T4) : Fin 3 × Fin 3 → Fin 3 × Fin 3 → K :=
  ∑ t : Fin 4 → Fin 3, S t • dpure [ev (t 0), ev (t 1), ev (t 2), ev (t 3)]

/-- `p = x ∧ y`. -/
noncomputable def wp : Fin 3 → Fin 3 → K := wedge (ev 0) (ev 1)
/-- `q = x ∧ z`. -/
noncomputable def wq : Fin 3 → Fin 3 → K := wedge (ev 0) (ev 2)
/-- `r = y ∧ z`. -/
noncomputable def wr : Fin 3 → Fin 3 → K := wedge (ev 1) (ev 2)

/-- The coordinates of an element of `Λ²(Λ²V)` on `p∧q`, `p∧r`, `q∧r`. -/
def compPQ (D : Fin 3 × Fin 3 → Fin 3 × Fin 3 → K) : K := D (0, 1) (0, 2)
/-- Coordinate on `p ∧ r`. -/
def compPR (D : Fin 3 × Fin 3 → Fin 3 × Fin 3 → K) : K := D (0, 1) (1, 2)
/-- Coordinate on `q ∧ r`. -/
def compQR (D : Fin 3 × Fin 3 → Fin 3 × Fin 3 → K) : K := D (0, 2) (1, 2)

/-- The shuffles of two lists. -/
def shuf {α : Type} : List α → List α → List (List α)
  | [], v => [v]
  | a :: u, [] => [a :: u]
  | a :: u, b :: v => (shuf u (b :: v)).map (a :: ·) ++ (shuf (a :: u) v).map (b :: ·)

/-- The table of the supplement (lines 93-117): basis tensor `x₁x₂x₃x₄` (`0 = x`, `1 = y`,
`2 = z`), its coefficient `a c^{-1} + b + d c` as `(a, b, d)`, and its image under `δ_{2,2}` as
integer coordinates on `p∧q`, `p∧r`, `q∧r`. -/
def Table24 : List ((Fin 4 → Fin 3) × ℚ × ℚ × ℚ × ℤ × ℤ × ℤ) := [
  (![1, 0, 0, 0], (0 : ℚ), (1 : ℚ), (-1 : ℚ), (0 : ℤ), (0 : ℤ), (0 : ℤ)),
  (![1, 0, 0, 1], (-1/2 : ℚ), (0 : ℚ), (1/2 : ℚ), (0 : ℤ), (0 : ℤ), (0 : ℤ)),
  (![1, 0, 0, 2], (1/2 : ℚ), (1 : ℚ), (1/2 : ℚ), (-1 : ℤ), (0 : ℤ), (0 : ℤ)),
  (![1, 0, 1, 0], (-1/2 : ℚ), (0 : ℚ), (1/2 : ℚ), (0 : ℤ), (0 : ℤ), (0 : ℤ)),
  (![1, 0, 1, 2], (1/2 : ℚ), (0 : ℚ), (-1/2 : ℚ), (0 : ℤ), (-1 : ℤ), (0 : ℤ)),
  (![1, 0, 2, 0], (1/2 : ℚ), (1 : ℚ), (1/2 : ℚ), (1 : ℤ), (0 : ℤ), (0 : ℤ)),
  (![1, 0, 2, 1], (-1/2 : ℚ), (1 : ℚ), (-1/2 : ℚ), (0 : ℤ), (1 : ℤ), (0 : ℤ)),
  (![1, 2, 0, 0], (0 : ℚ), (1 : ℚ), (1 : ℚ), (0 : ℤ), (0 : ℤ), (0 : ℤ)),
  (![1, 2, 0, 1], (-1/2 : ℚ), (1 : ℚ), (-1/2 : ℚ), (0 : ℤ), (-1 : ℤ), (0 : ℤ)),
  (![1, 2, 0, 2], (1/2 : ℚ), (0 : ℚ), (-1/2 : ℚ), (0 : ℤ), (0 : ℤ), (-1 : ℤ)),
  (![1, 2, 1, 0], (0 : ℚ), (1 : ℚ), (-1 : ℚ), (0 : ℤ), (1 : ℤ), (0 : ℤ)),
  (![1, 2, 1, 2], (0 : ℚ), (1 : ℚ), (1 : ℚ), (0 : ℤ), (0 : ℤ), (0 : ℤ)),
  (![2, 0, 0, 0], (0 : ℚ), (1 : ℚ), (1 : ℚ), (0 : ℤ), (0 : ℤ), (0 : ℤ)),
  (![2, 0, 0, 1], (-1/2 : ℚ), (1 : ℚ), (-1/2 : ℚ), (1 : ℤ), (0 : ℤ), (0 : ℤ)),
  (![2, 0, 0, 2], (1/2 : ℚ), (0 : ℚ), (-1/2 : ℚ), (0 : ℤ), (0 : ℤ), (0 : ℤ)),
  (![2, 0, 1, 0], (-1/2 : ℚ), (1 : ℚ), (-1/2 : ℚ), (-1 : ℤ), (0 : ℤ), (0 : ℤ)),
  (![2, 0, 1, 2], (1/2 : ℚ), (1 : ℚ), (1/2 : ℚ), (0 : ℤ), (0 : ℤ), (-1 : ℤ)),
  (![2, 0, 2, 0], (1/2 : ℚ), (0 : ℚ), (-1/2 : ℚ), (0 : ℤ), (0 : ℤ), (0 : ℤ)),
  (![2, 0, 2, 1], (-1/2 : ℚ), (0 : ℚ), (1/2 : ℚ), (0 : ℤ), (0 : ℤ), (1 : ℤ)),
  (![2, 1, 0, 0], (0 : ℚ), (1 : ℚ), (-1 : ℚ), (0 : ℤ), (0 : ℤ), (0 : ℤ)),
  (![2, 1, 0, 1], (-1/2 : ℚ), (0 : ℚ), (1/2 : ℚ), (0 : ℤ), (1 : ℤ), (0 : ℤ)),
  (![2, 1, 0, 2], (1/2 : ℚ), (1 : ℚ), (1/2 : ℚ), (0 : ℤ), (0 : ℤ), (1 : ℤ)),
  (![2, 1, 2, 0], (0 : ℚ), (1 : ℚ), (1 : ℚ), (0 : ℤ), (0 : ℤ), (-1 : ℤ)),
  (![2, 1, 2, 1], (0 : ℚ), (1 : ℚ), (-1 : ℚ), (0 : ℤ), (0 : ℤ), (0 : ℤ))]

/-- The tensor `Σ_rows (a c^{-1} + b + d c) · (basis tensor of the row)`. -/
noncomputable def tab24 : T4 := fun t =>
  (Table24.map fun row => if row.1 = t then
    ((row.2.1 : K) * cK⁻¹ + (row.2.2.1 : K) + (row.2.2.2.1 : K) * cK) else 0).sum

/-- `R_3(c) = Σ_{74 terms} q c^e H_w(c)` (real). -/
noncomputable def R3fun (c : ℝ) : ℝ :=
  (R3terms.map fun τ => (τ.1 : ℝ) * c ^ τ.2.1 * S03.H τ.2.2 c).sum

/-- `Ψ(c) = 1/2 + Σ_{14 terms} q c^e H_w(c)` (real). -/
noncomputable def PsiFun (c : ℝ) : ℝ :=
  1 / 2 + (PsiTerms.map fun τ => (τ.1 : ℝ) * c ^ τ.2.1 * S03.H τ.2.2 c).sum

/-- `[x^N] H_w = 𝔠_w ε^N N^{-κ} Z_{v_w}(N-1)` for `N ≥ 1` (Lemma `g6:lem:dict`, `S03.stmt_5`),
`0` for `N ≤ 0`. -/
def hcoef (w : List ℤ) (N : ℤ) : ℚ :=
  if 1 ≤ N then
    ((S03.dd w).2.2.1 : ℚ) * ((S03.dd w).1 : ℚ) ^ N.toNat / (N.toNat : ℚ) ^ (S03.dd w).2.1
      * S03.Z (S03.dd w).2.2.2 (N.toNat - 1)
  else 0

/-- The coefficient of `c^N` in the Laurent series of `Σ q c^e H_w(c)`: `Σ q [x^{N-e}] H_w`. -/
def coefTerms (L : List (ℚ × ℤ × List ℤ)) (N : ℤ) : ℚ :=
  (L.map fun τ => τ.1 * hcoef τ.2.2 (N - τ.2.1)).sum

/-- The data (main.tex 1544, supplement line 8): `R_3` has 74 terms, 2, 16 and 56 of weights
2, 3, 4 (no other weights), 34 distinct harmonic polylogarithms, all words regular over
`{0, 1, -1}`, exponents `e ∈ {-2, …, 1}`; `Ψ - 1/2` has 14 terms, of weights 1 and 2, regular words
over `{0,1,-1}`, exponents `e ≤ 1`. -/
def stmt_1 : Prop :=
  R3terms.length = 74 ∧
  (R3terms.filter fun τ => τ.2.2.length = 2).length = 2 ∧
  (R3terms.filter fun τ => τ.2.2.length = 3).length = 16 ∧
  Omega.length = 56 ∧
  (R3terms.map fun τ => τ.2.2).dedup.length = 34 ∧
  (∀ τ ∈ R3terms, S03.IsLetters τ.2.2 ∧ S03.Regular τ.2.2 ∧ -2 ≤ τ.2.1 ∧ τ.2.1 ≤ 1) ∧
  PsiTerms.length = 14 ∧
  (∀ τ ∈ PsiTerms, S03.IsLetters τ.2.2 ∧ S03.Regular τ.2.2 ∧
    (τ.2.2.length = 1 ∨ τ.2.2.length = 2) ∧ τ.2.1 ≤ 1)

/-- "What the words mean": `x = [c]`, `y = [1-c]`, `z = [1+c]` form a basis of `V`: a relation
`c^a (1-c)^b (1+c)^d = const` in `ℚ(c)` forces `a = b = d = 0`. -/
def stmt_2 : Prop :=
  ∀ a b d : ℤ, ∀ k : ℚ,
    cK ^ a * (1 - cK) ^ b * (1 + cK) ^ d = (k : K) → a = 0 ∧ b = 0 ∧ d = 0

/-- Lemma `g8:lem:vanish`: for every free module `ι → R` over a commutative ring, (a) `δ_{2,2}`
vanishes on the shuffle product of two pure tensors of degrees `(1,3)`, `(2,2)`, `(3,1)`, and (b)
on `-a ⊗ b ⊗ b ⊗ b` (the symbol of `Li_4`). -/
def stmt_3 : Prop :=
  (∀ (ι R : Type) [CommRing R] (s t : List (ι → R)), s ≠ [] → t ≠ [] →
      s.length + t.length = 4 → ((shuf s t).map dpure).sum = 0) ∧
  ∀ (ι R : Type) [CommRing R] (a b : ι → R), dpure [-a, b, b, b] = 0

/-- The table of 24 pure tensors (supplement lines 90-119): `𝒮(Ω)` is the sum of the 24 basis
tensors of the table with the printed coefficients (the 56 terms collapse to 24 tensors). -/
def stmt_4 : Prop := SOmega = tab24

/-- The `δ_{2,2}` column of the table: for each of the 24 rows, `δ_{2,2}` of its basis tensor is
the printed element `a p∧q + b p∧r + d q∧r` (`a, b, d ∈ {0, ±1}`). -/
def stmt_5 : Prop :=
  ∀ row ∈ Table24, delta22 (btens row.1) =
    ((row.2.2.2.2.1 : K) • wedge22 wp wq + (row.2.2.2.2.2.1 : K) • wedge22 wp wr
      + (row.2.2.2.2.2.2 : K) • wedge22 wq wr)

/-- Supplement lines 120-122: adding the rows, the components of `δ_{2,2}(𝒮(Ω))` along `p∧q`
cancel, the component along `p∧r` is `1 - 1/c` and the one along `q∧r` is `-(1 + 1/c)`. -/
def stmt_6 : Prop :=
  compPQ (delta22 SOmega) = 0 ∧ compPR (delta22 SOmega) = 1 - cK⁻¹ ∧
    compQR (delta22 SOmega) = -(1 + cK⁻¹)

/-- main.tex 1555 and supplement 152, 155: for `Ω` the two components `A` (on `p∧r`) and `B` (on
`q∧r`) satisfy `A (1+c) = B (1-c)`, i.e. `A/B = (1-c)/(1+c)`, and their ratio is not a rational
constant. -/
def stmt_7 : Prop :=
  compPR (delta22 SOmega) * (1 + cK) = compQR (delta22 SOmega) * (1 - cK) ∧
    compPR (delta22 SOmega) = (1 - cK) / (1 + cK) * compQR (delta22 SOmega) ∧
    ∀ k : ℚ, compPR (delta22 SOmega) ≠ (k : K) * compQR (delta22 SOmega)

/-- Proof of `g8:prop:one` (main.tex 1555, "constant ratio"): for any `U₀` with basis `ι`, any
`r ∈ ℚ(c)` and any tensor `τ ∈ T⁴(V ⊕ U₀)` with rational coefficients, there are rationals
`a, b, d` with `δ_{2,2}(r π^{⊗4}τ) = r (a p∧q + b p∧r + d q∧r)`. -/
def stmt_8 : Prop :=
  ∀ (ι : Type) (τ : (Fin 4 → Fin 3 ⊕ ι) → ℚ) (r : K), ∃ a b d : ℚ,
    delta22 (fun t => r * (τ (fun i => Sum.inl (t i)) : K)) =
      r • ((a : K) • wedge22 wp wq + (b : K) • wedge22 wp wr + (d : K) • wedge22 wq wr)

/-- Lemma `g8:prop:one`: for any `U₀` (basis `ι`), `r ∈ ℚ(c)^×` and `τ ∈ T⁴(V ⊕ U₀)` with rational
coefficients, `δ_{2,2}(r π^{⊗4}τ) ≠ δ_{2,2}(𝒮(Ω))`. -/
def stmt_9 : Prop :=
  ∀ (ι : Type) (τ : (Fin 4 → Fin 3 ⊕ ι) → ℚ) (r : K), r ≠ 0 →
    delta22 (fun t => r * (τ (fun i => Sum.inl (t i)) : K)) ≠ delta22 SOmega

/-- Proposition `prop:M3angle`, in the form of "What is proved, precisely": the constant-free
symbol of the 56-term expression `Ω`, coefficients in `ℚ(c)` as scalars, has
`δ_{2,2}(𝒮(Ω)) = (1 - 1/c) ([c]∧[1-c])∧([1-c]∧[1+c]) - (1 + 1/c) ([c]∧[1+c])∧([1-c]∧[1+c])`,
and it is not zero: both coefficients are nonzero elements of `ℚ(c)`, and `δ_{2,2}(𝒮(Ω)) ≠ 0`. -/
def stmt_10 : Prop :=
  delta22 SOmega = (1 - cK⁻¹) • wedge22 wp wr - (1 + cK⁻¹) • wedge22 wq wr ∧
    (1 - cK⁻¹) ≠ 0 ∧ (1 + cK⁻¹) ≠ 0 ∧ delta22 SOmega ≠ 0

/-- Machine step M5 for `R_3` (`g8:prop:crit` (i), (ii) with `𝔷 = 1`, through their values): for
every `n ≥ 1` the coefficient of `c^{2n}` in `Σ q c^e H_w` (74 terms) is `r'_n = [1]P_{3,n}` and
that of `c^{2n+1}` is `0`; the coefficients of `c^N`, `N ≤ 1`, vanish. -/
def stmt_11 : Prop :=
  (∀ n : ℕ, 1 ≤ n → ((coefTerms R3terms (2 * (n : ℤ)) : ℚ) : ℝ) = V05.P 3 n 0 0) ∧
  (∀ n : ℕ, 1 ≤ n → coefTerms R3terms (2 * (n : ℤ) + 1) = 0) ∧
  ∀ N : ℤ, N ≤ 1 → coefTerms R3terms N = 0

/-- Machine step M5 for `Ψ` (`𝔷 = 𝔷2`): for every `n ≥ 1` the coefficient of `c^{2n}` in the
14 terms is `q'_n = [𝔷2]P_{3,n}`, that of `c^{2n+1}` is `0`, and the coefficients of `c^N`,
`N ≤ 1`, of the 14 terms plus the constant `1/2` vanish. -/
def stmt_12 : Prop :=
  (∀ n : ℕ, 1 ≤ n →
    ((coefTerms PsiTerms (2 * (n : ℤ)) : ℚ) : ℝ) = V05.P 3 n 1 0 - V05.P 3 n 0 0) ∧
  (∀ n : ℕ, 1 ≤ n → coefTerms PsiTerms (2 * (n : ℤ) + 1) = 0) ∧
  ∀ N : ℤ, N ≤ 1 → coefTerms PsiTerms N + (if N = 0 then (1 / 2 : ℚ) else 0) = 0

/-- `g8:prop:crit` applied to `R_3` and `Ψ` (supplement line 8): for `0 < c < 1`,
`R_3(c) = Σ_{n≥1} r'_n c^{2n}` and `Ψ(c) = Σ_{n≥1} q'_n c^{2n}`. -/
def stmt_13 : Prop :=
  ∀ c : ℝ, 0 < c → c < 1 →
    HasSum (fun n : ℕ => V05.P 3 (n + 1) 0 0 * c ^ (2 * (n + 1))) (R3fun c) ∧
    HasSum (fun n : ℕ => (V05.P 3 (n + 1) 1 0 - V05.P 3 (n + 1) 0 0) * c ^ (2 * (n + 1)))
      (PsiFun c)

/-- Lemma `g8:prop:R3` (main.tex 1544): for `0 < α < π`, `c = cos(α/2)`,
`M_3(α) = R_3(c) - ζ(3) Φ(c) + ζ(2) Ψ(c)` with the 74 terms of `R_3` and the 14 terms and the
constant `1/2` of `Ψ`. -/
def stmt_14 : Prop :=
  ∀ α : ℝ, 0 < α → α < π →
    V05.Mk 3 α = R3fun (D04.c α) - V05.z3 * V05.Phi (D04.c α) + V05.z2 * PsiFun (D04.c α)

end P02
-- END STATEMENTS

theorem P02.proof_1 : P02.stmt_1 := by
  unfold P02.stmt_1 P02.Omega
  refine ⟨by decide +kernel, by decide +kernel, by decide +kernel, by decide +kernel,
    by decide +kernel, ?_, by decide +kernel, ?_⟩
  · intro τ hτ
    simp only [P02.R3terms, List.mem_cons, List.not_mem_nil, or_false] at hτ
    rcases hτ with h|h|h|h|h|h|h|h|h|h|h|h|h|h|h|h|h|h|h|h|h|h|h|h|h|h|h|h|h|h|h|h|h|h|h|h|h|h|h|h|h|h|h|h|h|h|h|h|h|h|h|h|h|h|h|h|h|h|h|h|h|h|h|h|h|h|h|h|h|h|h|h|h|h <;>
      subst h <;> simp [S03.IsLetters, S03.Regular]
  · intro τ hτ
    simp only [P02.PsiTerms, List.mem_cons, List.not_mem_nil, or_false] at hτ
    rcases hτ with h|h|h|h|h|h|h|h|h|h|h|h|h|h <;>
      subst h <;> simp [S03.IsLetters, S03.Regular]

namespace P02.H2

open Polynomial

lemma rm (A : ℚ[X]) (t : ℚ) (n : ℕ) (hA : A.eval t ≠ 0) :
    rootMultiplicity t (A * (X - C t) ^ n) = n := by
  have hA0 : A ≠ 0 := by rintro rfl; simp at hA
  rw [rootMultiplicity_mul_X_sub_C_pow hA0, rootMultiplicity_eq_zero (by simpa [IsRoot] using hA)]
  simp

/-- `C k X^m0 (1-X)^m1 (1+X)^m2`. -/
noncomputable def PP (k : ℚ) (m0 m1 m2 : ℕ) : ℚ[X] := C k * X ^ m0 * (1 - X) ^ m1 * (1 + X) ^ m2

lemma rm0 (k : ℚ) (hk : k ≠ 0) (m0 m1 m2 : ℕ) : rootMultiplicity 0 (PP k m0 m1 m2) = m0 := by
  have : PP k m0 m1 m2 = (C k * (1 - X) ^ m1 * (1 + X) ^ m2) * (X - C 0) ^ m0 := by
    simp [PP]; ring
  rw [this, rm]; simp [hk]

lemma rm1 (k : ℚ) (hk : k ≠ 0) (m0 m1 m2 : ℕ) : rootMultiplicity 1 (PP k m0 m1 m2) = m1 := by
  have h1 : (1 - X : ℚ[X]) = (-1) * (X - C 1) := by simp
  have : PP k m0 m1 m2 = (C k * (-1) ^ m1 * X ^ m0 * (1 + X) ^ m2) * (X - C 1) ^ m1 := by
    simp only [PP, h1, mul_pow]; ring
  rw [this, rm]; simp [hk]

lemma rm2 (k : ℚ) (hk : k ≠ 0) (m0 m1 m2 : ℕ) : rootMultiplicity (-1) (PP k m0 m1 m2) = m2 := by
  have h1 : (1 + X : ℚ[X]) = X - C (-1) := by simp [add_comm]
  have : PP k m0 m1 m2 = (C k * X ^ m0 * (1 - X) ^ m1) * (X - C (-1)) ^ m2 := by
    rw [PP, h1]
  rw [this, rm]; simp [hk]

lemma zp (y : P02.K) (hy : y ≠ 0) (a : ℤ) : y ^ a * y ^ ((-a).toNat) = y ^ (a.toNat) := by
  rw [← zpow_natCast, ← zpow_natCast, ← zpow_add₀ hy]; congr 1; omega

end P02.H2

theorem P02.proof_2 : P02.stmt_2 := by
  intro a b d k h
  open Polynomial in
  have hX : (P02.cK : P02.K) = algebraMap ℚ[X] P02.K X := (RatFunc.algebraMap_X).symm
  set u := algebraMap ℚ[X] P02.K with hu
  have inj : Function.Injective u := RatFunc.algebraMap_injective ℚ
  have h1 : 1 - P02.cK = u (1 - X) := by simp [hX]
  have h2 : 1 + P02.cK = u (1 + X) := by simp [hX]
  have nz : ∀ p : ℚ[X], p ≠ 0 → u p ≠ 0 := fun p hp h => hp (inj (by simpa using h))
  have y0 : P02.cK ≠ 0 := by rw [hX]; exact nz _ X_ne_zero
  have y1 : 1 - P02.cK ≠ 0 := by
    rw [h1]; refine nz _ (fun h => ?_); have := congrArg (eval 0) h; simp at this
  have y2 : 1 + P02.cK ≠ 0 := by
    rw [h2]; refine nz _ (fun h => ?_); have := congrArg (eval 0) h; simp at this
  have hk : k ≠ 0 := by
    rintro rfl
    simp at h
    rcases h with (h | h) | h <;> simp_all [zpow_ne_zero]
  have key := congrArg (· * (P02.cK ^ (-a).toNat * (1 - P02.cK) ^ (-b).toNat *
    (1 + P02.cK) ^ (-d).toNat)) h
  have e : P02.cK ^ a * (1 - P02.cK) ^ b * (1 + P02.cK) ^ d *
      (P02.cK ^ (-a).toNat * (1 - P02.cK) ^ (-b).toNat * (1 + P02.cK) ^ (-d).toNat) =
      (P02.cK ^ a * P02.cK ^ (-a).toNat) * ((1 - P02.cK) ^ b * (1 - P02.cK) ^ (-b).toNat) *
      ((1 + P02.cK) ^ d * (1 + P02.cK) ^ (-d).toNat) := by ring
  rw [e, P02.H2.zp _ y0, P02.H2.zp _ y1, P02.H2.zp _ y2] at key
  have key2 : u (P02.H2.PP 1 a.toNat b.toNat d.toNat) =
      u (P02.H2.PP k (-a).toNat (-b).toNat (-d).toNat) := by
    simp only [P02.H2.PP, map_mul, map_pow, ← h1, ← h2, ← hX]
    have hC : ∀ q : ℚ, u (C q) = (q : P02.K) := fun q => by
      rw [hu, RatFunc.algebraMap_C, ← RatFunc.algebraMap_eq_C]; exact (eq_ratCast _ q)
    rw [hC, hC, mul_assoc, mul_assoc, mul_assoc, mul_assoc, ← mul_assoc (_ ^ _), key]
    push_cast; ring
  have key3 := inj key2
  have r0 := P02.H2.rm0 1 one_ne_zero a.toNat b.toNat d.toNat
  have r1 := P02.H2.rm1 1 one_ne_zero a.toNat b.toNat d.toNat
  have r2 := P02.H2.rm2 1 one_ne_zero a.toNat b.toNat d.toNat
  rw [key3, P02.H2.rm0 k hk] at r0
  rw [key3, P02.H2.rm1 k hk] at r1
  rw [key3, P02.H2.rm2 k hk] at r2
  omega
theorem P02.proof_3 : P02.stmt_3 := by
  refine ⟨?_, ?_⟩
  · intro ι R _ s t hs ht hl
    rcases s with _ | ⟨a1, _ | ⟨a2, _ | ⟨a3, _ | ⟨a4, s'⟩⟩⟩⟩ <;>
    rcases t with _ | ⟨b1, _ | ⟨b2, _ | ⟨b3, _ | ⟨b4, t'⟩⟩⟩⟩ <;>
    simp at hs ht hl <;> try omega
    all_goals
      funext x y
      simp [P02.shuf, P02.dpure, P02.wedge22, P02.wedge]
      ring
  · intro ι R _ a b
    funext x y
    simp [P02.dpure, P02.wedge22, P02.wedge]
    ring

namespace P02.H4

open P02

def gi (a : ℤ) : Fin 3 := if a = 0 then 0 else if a = 1 then 1 else 2
def gs (a : ℤ) : ℤ := if a = 1 then -1 else 1

lemma g_eq (a : ℤ) : g a = ((gs a : ℤ) : K) • ev (gi a) := by
  unfold g gs gi
  by_cases h0 : a = 0
  · subst h0; simp
  · by_cases h1 : a = 1
    · subst h1; simp
    · simp [h0, h1]

lemma ev_apply (i j : Fin 3) : ev i j = if j = i then 1 else 0 := by
  unfold ev; rw [Pi.single_apply]

lemma prod_ind (l t : Fin 4 → Fin 3) :
    (∏ i, (if t i = l i then (1 : K) else 0)) = if t = l then 1 else 0 := by
  by_cases h : t = l
  · subst h; simp
  · rw [if_neg h]
    obtain ⟨i, hi⟩ : ∃ i, t i ≠ l i := by
      by_contra hc; exact h (funext fun i => not_not.mp fun hi => hc ⟨i, hi⟩)
    exact Finset.prod_eq_zero (Finset.mem_univ i) (by simp [hi])

lemma ptens_smul (s : Fin 4 → ℤ) (l : Fin 4 → Fin 3) (t : Fin 4 → Fin 3) :
    ptens (fun i => ((s i : ℤ) : K) • ev (l i)) t = ((∏ i, s i : ℤ) : K) * (if l = t then 1 else 0) := by
  unfold ptens
  simp only [Pi.smul_apply, smul_eq_mul, ev_apply, Finset.prod_mul_distrib, prod_ind, Int.cast_prod]
  by_cases h : t = l
  · subst h; simp
  · rw [if_neg h, if_neg (Ne.symm h)]

/-- The basis tensor index of `𝒮(H_w)`. -/
def lW (w : List ℤ) : Fin 4 → Fin 3 := fun i => gi (w.reverse.getD i 0)
/-- Its sign. -/
def sW (w : List ℤ) : ℤ := ∏ i : Fin 4, gs (w.reverse.getD i 0)

lemma symbH_eq (w : List ℤ) (t : Fin 4 → Fin 3) :
    symbH w t = ((sW w : ℤ) : K) * (if lW w = t then 1 else 0) := by
  unfold symbH
  simp_rw [g_eq]
  exact ptens_smul _ _ t

/-- A tensor given by a list of `(q, e, l)`: `Σ q c^e [l]`. -/
noncomputable def FL (L : List (ℚ × ℤ × (Fin 4 → Fin 3))) : T4 := fun t =>
  (L.map fun x => (x.1 : K) * cK ^ x.2.1 * (if x.2.2 = t then 1 else 0)).sum

def coefL (L : List (ℚ × ℤ × (Fin 4 → Fin 3))) (e : ℤ) (t : Fin 4 → Fin 3) : ℚ :=
  ((L.filter (fun x => decide (x.2.1 = e ∧ x.2.2 = t))).map Prod.fst).sum

lemma FL_eq (L : List (ℚ × ℤ × (Fin 4 → Fin 3))) (hL : ∀ x ∈ L, x.2.1 = -1 ∨ x.2.1 = 0 ∨ x.2.1 = 1)
    (t : Fin 4 → Fin 3) :
    FL L t = cK ^ (-1 : ℤ) * (coefL L (-1) t : K) + cK ^ (0 : ℤ) * (coefL L 0 t : K)
      + cK ^ (1 : ℤ) * (coefL L 1 t : K) := by
  induction L with
  | nil => simp [FL, coefL]
  | cons x L ih =>
    have ih' := ih (fun y hy => hL y (List.mem_cons_of_mem _ hy))
    unfold FL at ih' ⊢
    unfold coefL at ih' ⊢
    rw [List.map_cons, List.sum_cons, ih']
    obtain ⟨q, e, l⟩ := x
    have he := hL _ List.mem_cons_self
    simp only at he
    by_cases hl : l = t
    · subst hl
      rcases he with rfl | rfl | rfl <;> simp <;> ring
    · simp [hl]

def LOm : List (ℚ × ℤ × (Fin 4 → Fin 3)) :=
  Omega.map fun τ => (τ.1 * (sW τ.2.2 : ℚ), τ.2.1, lW τ.2.2)

def LTab : List (ℚ × ℤ × (Fin 4 → Fin 3)) :=
  Table24.flatMap fun row => [(row.2.1, -1, row.1), (row.2.2.1, 0, row.1), (row.2.2.2.1, 1, row.1)]

lemma SOmega_eq : SOmega = FL LOm := by
  funext t
  unfold SOmega FL LOm
  rw [List.map_map]
  congr 1
  apply List.map_congr_left
  intro τ _
  simp only [Function.comp, symbH_eq]
  push_cast
  ring

lemma tab24_eq : tab24 = FL LTab := by
  funext t
  unfold tab24 FL LTab
  induction Table24 with
  | nil => simp
  | cons row rows ih =>
    rw [List.map_cons, List.sum_cons, ih, List.flatMap_cons, List.map_append, List.sum_append]
    congr 1
    by_cases h : row.1 = t
    · simp [h, zpow_neg]; ring
    · simp [h]

lemma coef_check : ∀ t : Fin 4 → Fin 3, coefL LOm (-1) t = coefL LTab (-1) t ∧
    coefL LOm 0 t = coefL LTab 0 t ∧ coefL LOm 1 t = coefL LTab 1 t := by
  decide +kernel

end P02.H4

theorem P02.proof_4 : P02.stmt_4 := by
  unfold P02.stmt_4
  rw [P02.H4.SOmega_eq, P02.H4.tab24_eq]
  funext t
  rw [P02.H4.FL_eq _ (by decide +kernel), P02.H4.FL_eq _ (by decide +kernel)]
  obtain ⟨h1, h2, h3⟩ := P02.H4.coef_check t
  rw [h1, h2, h3]

namespace P02.H5

open P02

def evZ (i : Fin 3) : Fin 3 → ℤ := fun j => if j = i then 1 else 0

lemma ev_cast (i : Fin 3) : ev i = fun j => ((evZ i j : ℤ) : K) := by
  funext j; rw [H4.ev_apply]; unfold evZ; split_ifs <;> simp

/-- Integer model of `(e_{l0} ∧ e_{l1}) ∧ (e_{l2} ∧ e_{l3})`. -/
def dZ (l : Fin 4 → Fin 3) : Fin 3 × Fin 3 → Fin 3 × Fin 3 → ℤ :=
  wedge22 (wedge (evZ (l 0)) (evZ (l 1))) (wedge (evZ (l 2)) (evZ (l 3)))

noncomputable def castM (M : Fin 3 × Fin 3 → Fin 3 × Fin 3 → ℤ) : Fin 3 × Fin 3 → Fin 3 × Fin 3 → K :=
  fun a b => ((M a b : ℤ) : K)

lemma dpure_ev (l : Fin 4 → Fin 3) :
    dpure [ev (l 0), ev (l 1), ev (l 2), ev (l 3)] = castM (dZ l) := by
  funext a b
  simp only [dpure, dZ, castM, wedge22, wedge, ev_cast]
  push_cast; ring

lemma delta22_btens (l : Fin 4 → Fin 3) : delta22 (btens l) = castM (dZ l) := by
  rw [← dpure_ev]
  unfold delta22 btens
  simp

def wpZ : Fin 3 → Fin 3 → ℤ := wedge (evZ 0) (evZ 1)
def wqZ : Fin 3 → Fin 3 → ℤ := wedge (evZ 0) (evZ 2)
def wrZ : Fin 3 → Fin 3 → ℤ := wedge (evZ 1) (evZ 2)

lemma w22_cast (P Q : Fin 3 → Fin 3 → ℤ) :
    wedge22 (fun i j => ((P i j : ℤ) : K)) (fun i j => ((Q i j : ℤ) : K)) = castM (wedge22 P Q) := by
  funext a b; simp [wedge22, castM]

lemma wp_cast : wp = fun i j => ((wpZ i j : ℤ) : K) := by
  funext i j; simp [wp, wpZ, wedge, ev_cast]
lemma wq_cast : wq = fun i j => ((wqZ i j : ℤ) : K) := by
  funext i j; simp [wq, wqZ, wedge, ev_cast]
lemma wr_cast : wr = fun i j => ((wrZ i j : ℤ) : K) := by
  funext i j; simp [wr, wrZ, wedge, ev_cast]

lemma comb_cast (a b d : ℤ) :
    (a : K) • wedge22 wp wq + (b : K) • wedge22 wp wr + (d : K) • wedge22 wq wr =
      castM (a • wedge22 wpZ wqZ + b • wedge22 wpZ wrZ + d • wedge22 wqZ wrZ) := by
  rw [wp_cast, wq_cast, wr_cast, w22_cast, w22_cast, w22_cast]
  funext x y; simp [castM]

lemma rows_check : ∀ row ∈ Table24, dZ row.1 =
    row.2.2.2.2.1 • wedge22 wpZ wqZ + row.2.2.2.2.2.1 • wedge22 wpZ wrZ
      + row.2.2.2.2.2.2 • wedge22 wqZ wrZ := by
  decide +kernel

end P02.H5

theorem P02.proof_5 : P02.stmt_5 := by
  intro row hrow
  rw [P02.H5.delta22_btens, P02.H5.comb_cast, P02.H5.rows_check row hrow]

namespace P02.H6

open P02

lemma delta22_add (S1 S2 : T4) : delta22 (S1 + S2) = delta22 S1 + delta22 S2 := by
  unfold delta22; rw [← Finset.sum_add_distrib]; congr 1; funext t; rw [Pi.add_apply, add_smul]

lemma delta22_zero : delta22 0 = 0 := by unfold delta22; simp

lemma delta22_smul (c : K) (S : T4) : delta22 (c • S) = c • delta22 S := by
  unfold delta22; rw [Finset.smul_sum]; congr 1; funext t; rw [Pi.smul_apply, smul_eq_mul, mul_smul]

lemma delta22_list {α : Type} (rows : List α) (l : α → (Fin 4 → Fin 3)) (f : α → K) :
    delta22 (fun t => (rows.map fun row => if l row = t then f row else 0).sum) =
      (rows.map fun row => f row • delta22 (btens (l row))).sum := by
  induction rows with
  | nil => simp only [List.map_nil, List.sum_nil]; exact delta22_zero
  | cons row rows ih =>
    have : (fun t => ((row :: rows).map fun row => if l row = t then f row else 0).sum) =
        f row • btens (l row) + (fun t => (rows.map fun row => if l row = t then f row else 0).sum) := by
      funext t; simp only [List.map_cons, List.sum_cons, Pi.add_apply, Pi.smul_apply, btens,
        smul_eq_mul]
      congr 1; split_ifs with h1 h2 h2 <;> simp_all [eq_comm]
    rw [this, delta22_add, delta22_smul, ih, List.map_cons, List.sum_cons]

lemma list_comb {α : Type} (rows : List α) (f a b d : α → K) (X Y Z : Fin 3 × Fin 3 → Fin 3 × Fin 3 → K) :
    (rows.map fun row => f row • (a row • X + b row • Y + d row • Z)).sum =
      (rows.map fun row => f row * a row).sum • X + (rows.map fun row => f row * b row).sum • Y
        + (rows.map fun row => f row * d row).sum • Z := by
  induction rows with
  | nil => simp
  | cons row rows ih =>
    rw [List.map_cons, List.sum_cons, ih]
    simp only [List.map_cons, List.sum_cons, add_smul, smul_add, mul_smul]
    abel

lemma cK_ne_const (q : ℚ) : cK ≠ (q : K) := by
  intro h
  have h' : algebraMap (Polynomial ℚ) K Polynomial.X = algebraMap (Polynomial ℚ) K (Polynomial.C q) := by
    rw [RatFunc.algebraMap_X, RatFunc.algebraMap_C, ← RatFunc.algebraMap_eq_C]
    exact h.trans (eq_ratCast _ q).symm
  exact Polynomial.X_ne_C q (RatFunc.algebraMap_injective ℚ h')

lemma cK_ne_zero : cK ≠ 0 := by simpa using cK_ne_const 0

lemma one_sub_inv_ne : (1 : K) - cK⁻¹ ≠ 0 := by
  intro h
  have : cK⁻¹ = 1 := by linear_combination -h
  have : cK = ((1 : ℚ) : K) := by rw [inv_eq_one] at this; simpa using this
  exact cK_ne_const 1 this

lemma one_add_inv_ne : (1 : K) + cK⁻¹ ≠ 0 := by
  intro h
  have : cK⁻¹ = -1 := by linear_combination h
  have : cK = ((-1 : ℚ) : K) := by
    have := congrArg (·⁻¹) this; simp only [inv_inv] at this; rw [this]; simp
  exact cK_ne_const (-1) this

lemma delta22_SOmega : delta22 SOmega = (1 - cK⁻¹) • wedge22 wp wr - (1 + cK⁻¹) • wedge22 wq wr := by
  rw [P02.proof_4]
  unfold tab24
  rw [delta22_list Table24 (fun row => row.1)]
  have : ∀ row ∈ Table24, delta22 (btens row.1) =
      ((row.2.2.2.2.1 : K) • wedge22 wp wq + (row.2.2.2.2.2.1 : K) • wedge22 wp wr
        + (row.2.2.2.2.2.2 : K) • wedge22 wq wr) := P02.proof_5
  rw [List.map_congr_left (fun row hrow => by rw [this row hrow]), list_comb]
  simp only [Table24, List.map_cons, List.map_nil, List.sum_cons, List.sum_nil]
  push_cast
  have e1 : ∀ X Y Z : Fin 3 × Fin 3 → Fin 3 × Fin 3 → K, ∀ A B C A' B' C' : K, A = A' → B = B' → C = C' →
      A • X + B • Y + C • Z = A' • X + B' • Y + C' • Z := by intros; subst_vars; rfl
  rw [e1 _ _ _ _ _ _ 0 (1 - cK⁻¹) (-(1 + cK⁻¹)) (by ring) (by ring) (by ring)]
  rw [zero_smul, zero_add, neg_smul]; abel

lemma comps (A B : K) :
    compPQ (A • wedge22 wp wr - B • wedge22 wq wr) = 0 ∧
    compPR (A • wedge22 wp wr - B • wedge22 wq wr) = A ∧
    compQR (A • wedge22 wp wr - B • wedge22 wq wr) = -B := by
  simp [compPQ, compPR, compQR, wedge22, wp, wq, wr, wedge, H4.ev_apply]

end P02.H6

theorem P02.proof_6 : P02.stmt_6 := by
  unfold P02.stmt_6
  rw [P02.H6.delta22_SOmega]
  obtain ⟨h1, h2, h3⟩ := P02.H6.comps (1 - P02.cK⁻¹) (1 + P02.cK⁻¹)
  exact ⟨h1, h2, h3⟩
theorem P02.proof_7 : P02.stmt_7 := by
  unfold P02.stmt_7
  obtain ⟨_, h2, h3⟩ := P02.proof_6
  rw [h2, h3]
  have hc := P02.H6.cK_ne_zero
  have hc1 : (1 : P02.K) + P02.cK ≠ 0 := by
    intro h
    exact P02.H6.cK_ne_const (-1) (by push_cast; linear_combination h)
  refine ⟨?_, ?_, ?_⟩
  · field_simp; ring
  · field_simp; ring
  · intro k hk
    have e : (1 + (k : P02.K)) * P02.cK = 1 - k := by
      field_simp at hk; linear_combination hk
    by_cases hk1 : (1 + (k : P02.K)) = 0
    · rw [hk1, zero_mul] at e
      have : (k : P02.K) = 1 := by linear_combination e
      rw [this] at hk1
      norm_num at hk1
    · apply P02.H6.cK_ne_const ((1 - k) / (1 + k))
      push_cast
      field_simp
      linear_combination e
namespace P02.H8

open P02 H5

lemma dZ_decomp : ∀ t : Fin 4 → Fin 3, dZ t = dZ t (0, 1) (0, 2) • wedge22 wpZ wqZ +
    dZ t (0, 1) (1, 2) • wedge22 wpZ wrZ + dZ t (0, 2) (1, 2) • wedge22 wqZ wrZ := by
  decide +kernel

lemma delta22_expand (S : T4) : delta22 S =
    (∑ t, S t * ((dZ t (0, 1) (0, 2) : ℤ) : K)) • wedge22 wp wq +
    (∑ t, S t * ((dZ t (0, 1) (1, 2) : ℤ) : K)) • wedge22 wp wr +
    (∑ t, S t * ((dZ t (0, 2) (1, 2) : ℤ) : K)) • wedge22 wq wr := by
  unfold delta22
  have : ∀ t : Fin 4 → Fin 3, dpure [ev (t 0), ev (t 1), ev (t 2), ev (t 3)] =
      ((dZ t (0, 1) (0, 2) : ℤ) : K) • wedge22 wp wq + ((dZ t (0, 1) (1, 2) : ℤ) : K) • wedge22 wp wr
        + ((dZ t (0, 2) (1, 2) : ℤ) : K) • wedge22 wq wr := by
    intro t; rw [dpure_ev, comb_cast, ← dZ_decomp t]
  simp only [this, smul_add, Finset.sum_add_distrib, Finset.sum_smul, mul_smul]

lemma comps3 (A B C : K) :
    compPQ (A • wedge22 wp wq + B • wedge22 wp wr + C • wedge22 wq wr) = A ∧
    compPR (A • wedge22 wp wq + B • wedge22 wp wr + C • wedge22 wq wr) = B ∧
    compQR (A • wedge22 wp wq + B • wedge22 wp wr + C • wedge22 wq wr) = C := by
  simp [compPQ, compPR, compQR, wedge22, wp, wq, wr, wedge, H4.ev_apply]

end P02.H8

theorem P02.proof_8 : P02.stmt_8 := by
  intro ι τ r
  refine ⟨∑ t : Fin 4 → Fin 3, τ (fun i => Sum.inl (t i)) * ((P02.H5.dZ t (0, 1) (0, 2) : ℤ) : ℚ),
    ∑ t : Fin 4 → Fin 3, τ (fun i => Sum.inl (t i)) * ((P02.H5.dZ t (0, 1) (1, 2) : ℤ) : ℚ),
    ∑ t : Fin 4 → Fin 3, τ (fun i => Sum.inl (t i)) * ((P02.H5.dZ t (0, 2) (1, 2) : ℤ) : ℚ), ?_⟩
  rw [P02.H8.delta22_expand]
  simp only [smul_add, smul_smul]
  push_cast
  simp only [Finset.mul_sum, mul_assoc]
theorem P02.proof_9 : P02.stmt_9 := by
  intro ι τ r _ h
  obtain ⟨a, b, d, hab⟩ := P02.proof_8 ι τ r
  rw [hab] at h
  have hPR := P02.proof_6.2.1
  have hQR := P02.proof_6.2.2
  rw [← h, smul_add, smul_add, smul_smul, smul_smul, smul_smul] at hPR hQR
  rw [(P02.H8.comps3 _ _ _).2.1] at hPR
  rw [(P02.H8.comps3 _ _ _).2.2] at hQR
  by_cases hd : (d : P02.K) = 0
  · rw [hd, mul_zero] at hQR
    exact P02.H6.one_add_inv_ne (by linear_combination hQR)
  · apply (P02.proof_7.2.2 (b / d))
    rw [← h, smul_add, smul_add, smul_smul, smul_smul, smul_smul, (P02.H8.comps3 _ _ _).2.1,
      (P02.H8.comps3 _ _ _).2.2]
    push_cast
    field_simp
theorem P02.proof_10 : P02.stmt_10 := by
  refine ⟨P02.H6.delta22_SOmega, P02.H6.one_sub_inv_ne, P02.H6.one_add_inv_ne, ?_⟩
  intro h
  have := P02.proof_6.2.1
  rw [h] at this
  exact P02.H6.one_sub_inv_ne (by simpa [P02.compPR] using this.symm)

namespace P02.H12

open P02

lemma hcoef_succ (w : List ℤ) (k : ℕ) : hcoef w ((k : ℤ) + 1) =
    ((S03.dd w).2.2.1 : ℚ) * ((S03.dd w).1 : ℚ) ^ (k + 1) / ((k : ℚ) + 1) ^ (S03.dd w).2.1
      * S03.Z (S03.dd w).2.2.2 k := by
  unfold hcoef
  rw [if_pos (by omega), show ((k : ℤ) + 1).toNat = k + 1 by omega]
  push_cast
  simp

lemma hcoef_nonpos (w : List ℤ) (N : ℤ) (h : N ≤ 0) : hcoef w N = 0 := by
  unfold hcoef; rw [if_neg (by omega)]

lemma Zm1_succ (k : ℕ) : S03.Z [-1] (k + 1) = S03.Z [-1] k + (-1) ^ (k + 1) / ((k : ℚ) + 1) := by
  rw [S03.Z_cons_succ]; simp [S03.gq, S03.Z]

lemma Znil (k : ℕ) : S03.Z [] k = 1 := rfl

lemma Zm1_even (m : ℕ) : ((S03.Z [-1] (2 * m) : ℚ) : ℝ) = V05.Hr 1 m - V05.Hr 1 (2 * m) := by
  induction m with
  | zero => simp [S03.Z, V05.Hr]
  | succ m ih =>
    rw [show 2 * (m + 1) = 2 * m + 1 + 1 by ring, Zm1_succ, Zm1_succ, V05.Hr_succ, V05.Hr_succ,
      V05.Hr_succ]
    push_cast
    rw [ih]
    have : (0 : ℝ) ≤ m := by positivity
    simp only [pow_succ, pow_mul, pow_zero, neg_one_sq, one_pow]
    field_simp
    ring

lemma dd1 : S03.dd [1] = (1, 1, 1, []) := rfl
lemma ddm1 : S03.dd [-1] = (-1, 1, -1, []) := rfl
lemma dd01 : S03.dd [0, 1] = (1, 2, 1, []) := rfl
lemma dd0m1 : S03.dd [0, -1] = (-1, 2, -1, []) := rfl
lemma dd1m1 : S03.dd [1, -1] = (1, 1, -1, [-1]) := rfl
lemma ddm11 : S03.dd [-1, 1] = (-1, 1, -1, [-1]) := rfl

lemma e321_eq (n : ℕ) : V05.ecoef 3 n 2 1 =
    Polynomial.C (1/2 * (n : ℝ) - n ^ 2 + 1/2 * n ^ 3)
    + Polynomial.C (-1/2 * (n : ℝ) + 5/2 * n ^ 2 - 3 * n ^ 3) * Polynomial.X
    + Polynomial.C (-(n : ℝ) - n ^ 2 + 6 * n ^ 3) * Polynomial.X ^ 2
    + Polynomial.C ((n : ℝ) - 4 * n ^ 3) * Polynomial.X ^ 3 := by
  apply Polynomial.funext; intro y
  simp [V05.ecoef, V05.al1, V05.be1, V05.ga, V05.ga3]; ring

lemma e312_eq (n : ℕ) : V05.ecoef 3 n 1 2 =
    Polynomial.C (1/2 * (n : ℝ) ^ 2 - 1/2 * n ^ 3)
    + Polynomial.C (-1/2 * (n : ℝ) - 1/2 * n ^ 2 + 3 * n ^ 3) * Polynomial.X
    + Polynomial.C (2 * (n : ℝ) - n ^ 2 - 6 * n ^ 3) * Polynomial.X ^ 2
    + Polynomial.C (-(n : ℝ) + 4 * n ^ 3) * Polynomial.X ^ 3 := by
  apply Polynomial.funext; intro y
  simp [V05.ecoef, V05.al1, V05.be1, V05.ga, V05.ga3]; ring

lemma part1 (m : ℕ) : ((coefTerms PsiTerms (2 * ((m + 1 : ℕ) : ℤ)) : ℚ) : ℝ) =
    V05.P 3 (m + 1) 1 0 - V05.P 3 (m + 1) 0 0 := by
  rw [V05.P3_x2]
  unfold coefTerms PsiTerms
  simp only [List.map_cons, List.map_nil, List.sum_cons, List.sum_nil]
  rw [show 2 * ((m + 1 : ℕ) : ℤ) - -2 = ((2 * m + 1 + 1 + 1 : ℕ) : ℤ) + 1 by push_cast; ring,
    show 2 * ((m + 1 : ℕ) : ℤ) - -1 = ((2 * m + 1 + 1 : ℕ) : ℤ) + 1 by push_cast; ring,
    show 2 * ((m + 1 : ℕ) : ℤ) - 0 = ((2 * m + 1 : ℕ) : ℤ) + 1 by push_cast; ring,
    show 2 * ((m + 1 : ℕ) : ℤ) - 1 = ((2 * m : ℕ) : ℤ) + 1 by push_cast; ring]
  simp only [hcoef_succ, dd1, ddm1, dd01, dd0m1, dd1m1, ddm11, Znil, Zm1_succ]
  push_cast
  rw [Zm1_even]
  simp only [V05.L2, Finset.sum_range_succ, Finset.sum_range_zero, V05.D2, V05.rho,
    Finset.prod_range_succ, Finset.prod_range_zero]
  norm_num
  simp only [e321_eq, e312_eq, V05.coeff3]
  norm_num
  rw [show m + 1 + m = 2 * m + 1 by ring, show m + 1 + 1 + m = 2 * m + 1 + 1 by ring,
    show m + 1 + 2 + m = 2 * m + 1 + 1 + 1 by ring, show m + 1 + 3 + m = 2 * m + 1 + 1 + 1 + 1 by ring,
    show m + 2 = m + 1 + 1 by ring, show m + 3 = m + 1 + 1 + 1 by ring]
  simp only [V05.Hr_succ, pow_succ, pow_mul, neg_one_sq, one_pow]
  push_cast
  have : (0 : ℝ) ≤ m := by positivity
  have h1 : (2 * (m : ℝ) + 1) ≠ 0 := by positivity
  have h2 : (2 * ((m : ℝ) + 1) - 1) ≠ 0 := by linarith
  field_simp
  ring

lemma part2 (m : ℕ) : coefTerms PsiTerms (2 * ((m + 1 : ℕ) : ℤ) + 1) = 0 := by
  unfold coefTerms PsiTerms
  simp only [List.map_cons, List.map_nil, List.sum_cons, List.sum_nil]
  rw [show 2 * ((m + 1 : ℕ) : ℤ) + 1 - -2 = ((2 * m + 1 + 1 + 1 + 1 : ℕ) : ℤ) + 1 by push_cast; ring,
    show 2 * ((m + 1 : ℕ) : ℤ) + 1 - -1 = ((2 * m + 1 + 1 + 1 : ℕ) : ℤ) + 1 by push_cast; ring,
    show 2 * ((m + 1 : ℕ) : ℤ) + 1 - 0 = ((2 * m + 1 + 1 : ℕ) : ℤ) + 1 by push_cast; ring,
    show 2 * ((m + 1 : ℕ) : ℤ) + 1 - 1 = ((2 * m + 1 : ℕ) : ℤ) + 1 by push_cast; ring]
  simp only [hcoef_succ, dd1, ddm1, dd01, dd0m1, dd1m1, ddm11, Znil, Zm1_succ]
  push_cast
  simp only [pow_succ, pow_mul, neg_one_sq, one_pow]
  have : (0 : ℚ) ≤ m := by positivity
  have h1 : (2 * (m : ℚ) + 1) ≠ 0 := by positivity
  field_simp
  ring

lemma coefTerms_low (L : List (ℚ × ℤ × List ℤ)) (N : ℤ) (h : ∀ τ ∈ L, N ≤ τ.2.1) :
    coefTerms L N = 0 := by
  unfold coefTerms
  induction L with
  | nil => simp
  | cons τ L ih =>
    rw [List.map_cons, List.sum_cons, ih (fun σ hσ => h σ (List.mem_cons_of_mem _ hσ)),
      hcoef_nonpos _ _ (by have := h τ List.mem_cons_self; omega)]
    simp

lemma part3 : ∀ N : ℤ, N ≤ 1 → coefTerms PsiTerms N + (if N = 0 then (1 / 2 : ℚ) else 0) = 0 := by
  intro N hN
  rcases (show N ≤ -2 ∨ N = -1 ∨ N = 0 ∨ N = 1 by omega) with h | rfl | rfl | rfl
  · rw [coefTerms_low _ _ (fun τ hτ => by
      simp only [PsiTerms, List.mem_cons, List.not_mem_nil, or_false] at hτ
      rcases hτ with h' | h' | h' | h' | h' | h' | h' | h' | h' | h' | h' | h' | h' | h' <;>
        subst h' <;> simp only <;> omega)]
    rw [if_neg (by omega)]; simp
  · decide +kernel
  · decide +kernel
  · decide +kernel

end P02.H12

namespace P02.H11

open P02 H12

lemma ddR0 : S03.dd [0, 1] = (1, 2, 1, []) := rfl
lemma ddR1 : S03.dd [0, -1] = (-1, 2, -1, []) := rfl
lemma ddR2 : S03.dd [0, 0, 1] = (1, 3, 1, []) := rfl
lemma ddR3 : S03.dd [0, 0, -1] = (-1, 3, -1, []) := rfl
lemma ddR4 : S03.dd [0, 1, -1] = (1, 2, -1, [-1]) := rfl
lemma ddR5 : S03.dd [0, -1, 1] = (-1, 2, -1, [-1]) := rfl
lemma ddR6 : S03.dd [1, 0, 1] = (1, 1, 1, [2]) := rfl
lemma ddR7 : S03.dd [1, 0, -1] = (1, 1, -1, [-2]) := rfl
lemma ddR8 : S03.dd [-1, 0, 1] = (-1, 1, -1, [-2]) := rfl
lemma ddR9 : S03.dd [-1, 0, -1] = (-1, 1, 1, [2]) := rfl
lemma ddR10 : S03.dd [0, 0, 0, 1] = (1, 4, 1, []) := rfl
lemma ddR11 : S03.dd [0, 0, 0, -1] = (-1, 4, -1, []) := rfl
lemma ddR12 : S03.dd [0, 0, 1, -1] = (1, 3, -1, [-1]) := rfl
lemma ddR13 : S03.dd [0, 0, -1, 1] = (-1, 3, -1, [-1]) := rfl
lemma ddR14 : S03.dd [0, 1, 0, 1] = (1, 2, 1, [2]) := rfl
lemma ddR15 : S03.dd [0, 1, 0, -1] = (1, 2, -1, [-2]) := rfl
lemma ddR16 : S03.dd [0, 1, -1, 1] = (1, 2, -1, [-1, -1]) := rfl
lemma ddR17 : S03.dd [0, -1, 0, 1] = (-1, 2, -1, [-2]) := rfl
lemma ddR18 : S03.dd [0, -1, 0, -1] = (-1, 2, 1, [2]) := rfl
lemma ddR19 : S03.dd [0, -1, 1, -1] = (-1, 2, 1, [-1, -1]) := rfl
lemma ddR20 : S03.dd [1, 0, 0, 1] = (1, 1, 1, [3]) := rfl
lemma ddR21 : S03.dd [1, 0, 0, -1] = (1, 1, -1, [-3]) := rfl
lemma ddR22 : S03.dd [1, 0, 1, -1] = (1, 1, -1, [2, -1]) := rfl
lemma ddR23 : S03.dd [1, 0, -1, 1] = (1, 1, -1, [-2, -1]) := rfl
lemma ddR24 : S03.dd [1, -1, 0, 1] = (1, 1, -1, [-1, -2]) := rfl
lemma ddR25 : S03.dd [1, -1, 0, -1] = (1, 1, 1, [-1, 2]) := rfl
lemma ddR26 : S03.dd [1, -1, 1, -1] = (1, 1, 1, [-1, -1, -1]) := rfl
lemma ddR27 : S03.dd [-1, 0, 0, 1] = (-1, 1, -1, [-3]) := rfl
lemma ddR28 : S03.dd [-1, 0, 0, -1] = (-1, 1, 1, [3]) := rfl
lemma ddR29 : S03.dd [-1, 0, 1, -1] = (-1, 1, 1, [-2, -1]) := rfl
lemma ddR30 : S03.dd [-1, 0, -1, 1] = (-1, 1, 1, [2, -1]) := rfl
lemma ddR31 : S03.dd [-1, 1, 0, 1] = (-1, 1, -1, [-1, 2]) := rfl
lemma ddR32 : S03.dd [-1, 1, 0, -1] = (-1, 1, 1, [-1, -2]) := rfl
lemma ddR33 : S03.dd [-1, 1, -1, 1] = (-1, 1, 1, [-1, -1, -1]) := rfl

lemma sg2 : Int.sign 2 = 1 := rfl
lemma sg3 : Int.sign 3 = 1 := rfl
lemma gqm1 (k : ℕ) : S03.gq (-1) k = (-1) ^ k / (k : ℚ) := by simp [S03.gq]
lemma gq2 (k : ℕ) : S03.gq 2 k = 1 / (k : ℚ) ^ 2 := by simp [S03.gq, sg2]
lemma gqm2 (k : ℕ) : S03.gq (-2) k = (-1) ^ k / (k : ℚ) ^ 2 := by simp [S03.gq, sg2]
lemma gq3 (k : ℕ) : S03.gq 3 k = 1 / (k : ℚ) ^ 3 := by simp [S03.gq, sg3]
lemma gqm3 (k : ℕ) : S03.gq (-3) k = (-1) ^ k / (k : ℚ) ^ 3 := by simp [S03.gq, sg3]

lemma st2 (a b c : ℤ) (hc : ∀ i, S03.gq a i * S03.gq b i = S03.gq c i) (M : ℕ) :
    S03.Z [a] M * S03.Z [b] M = S03.Z [a, b] M + S03.Z [b, a] M + S03.Z [c] M := by
  induction M with
  | zero => simp [S03.Z]
  | succ M ih =>
    simp only [S03.Z_cons_succ, Znil]
    rw [← hc]
    linear_combination ih

lemma st3 (a b c ab ac : ℤ) (hab : ∀ i, S03.gq a i * S03.gq b i = S03.gq ab i)
    (hac : ∀ i, S03.gq a i * S03.gq c i = S03.gq ac i) (M : ℕ) :
    S03.Z [a] M * S03.Z [b, c] M = S03.Z [a, b, c] M + S03.Z [b, a, c] M + S03.Z [b, c, a] M
      + S03.Z [ab, c] M + S03.Z [b, ac] M := by
  induction M with
  | zero => simp [S03.Z]
  | succ M ih =>
    simp only [S03.Z_cons_succ, Znil]
    rw [← hab]
    linear_combination ih + S03.gq b (M + 1) * st2 a c ac hac M

lemma Hr_Z (r : ℕ) (hr : 1 ≤ r) (N : ℕ) : V05.Hr r N = ((S03.Z [(r : ℤ)] N : ℚ) : ℝ) := by
  induction N with
  | zero => simp [V05.Hr, S03.Z]
  | succ N ih =>
    rw [V05.Hr_succ, ih, S03.Z_cons_succ]
    have : Int.sign (r : ℤ) = 1 := Int.sign_eq_one_of_pos (by omega)
    simp [S03.gq, this, Znil]

lemma Hr_half (r : ℕ) (hr : 1 ≤ r) (m : ℕ) : V05.Hr r m =
    (2 : ℝ) ^ (r - 1) * (((S03.Z [(r : ℤ)] (2 * m) : ℚ) : ℝ) + ((S03.Z [-(r : ℤ)] (2 * m) : ℚ) : ℝ)) := by
  have := S03.stmt3_aux r hr m
  rw [S08.PA.hn_eq_HrQ] at this
  rw [← V05.HrQ_cast, ← this]
  push_cast; ring

lemma e330_eq (n : ℕ) : V05.ecoef 3 n 3 0 =
    Polynomial.C (1/6 - 1/2 * (n : ℝ) + 1/2 * n ^ 2 - 1/6 * n ^ 3)
    + Polynomial.C (-1/6 + 5/6 * (n : ℝ) - 3/2 * n ^ 2 + n ^ 3) * Polynomial.X
    + Polynomial.C ((n : ℝ) ^ 2 - 2 * n ^ 3) * Polynomial.X ^ 2
    + Polynomial.C (-1/3 * (n : ℝ) + 4/3 * n ^ 3) * Polynomial.X ^ 3 := by
  apply Polynomial.funext; intro y
  simp [V05.ecoef, V05.al1, V05.be1, V05.ga, V05.ga3]; ring

lemma e303_eq (n : ℕ) : V05.ecoef 3 n 0 3 =
    Polynomial.C (1/6 * (n : ℝ) ^ 3)
    + Polynomial.C (1/6 + 1/6 * (n : ℝ) - 1/2 * n ^ 2 - n ^ 3) * Polynomial.X
    + Polynomial.C (-(n : ℝ) + n ^ 2 + 2 * n ^ 3) * Polynomial.X ^ 2
    + Polynomial.C (1/3 * (n : ℝ) - 4/3 * n ^ 3) * Polynomial.X ^ 3 := by
  apply Polynomial.funext; intro y
  simp [V05.ecoef, V05.al1, V05.be1, V05.ga, V05.ga3]; ring


lemma Hr1m (m : ℕ) : V05.Hr 1 m = ((S03.Z [1] (2 * m) : ℚ) : ℝ) + ((S03.Z [-1] (2 * m) : ℚ) : ℝ) := by
  have := Hr_half 1 le_rfl m; simpa using this
lemma Hr2m (m : ℕ) : V05.Hr 2 m = 2 * (((S03.Z [2] (2 * m) : ℚ) : ℝ) + ((S03.Z [-2] (2 * m) : ℚ) : ℝ)) := by
  have := Hr_half 2 (by norm_num) m; simpa using this
lemma Hr3m (m : ℕ) : V05.Hr 3 m = 4 * (((S03.Z [3] (2 * m) : ℚ) : ℝ) + ((S03.Z [-3] (2 * m) : ℚ) : ℝ)) := by
  have := Hr_half 3 (by norm_num) m; norm_num at this; linear_combination this
lemma Hr12m (m : ℕ) : V05.Hr 1 (2 * m) = ((S03.Z [1] (2 * m) : ℚ) : ℝ) := by
  have := Hr_Z 1 le_rfl (2 * m); simpa using this
lemma Hr22m (m : ℕ) : V05.Hr 2 (2 * m) = ((S03.Z [2] (2 * m) : ℚ) : ℝ) := by
  have := Hr_Z 2 (by norm_num) (2 * m); simpa using this
lemma Hr32m (m : ℕ) : V05.Hr 3 (2 * m) = ((S03.Z [3] (2 * m) : ℚ) : ℝ) := by
  have := Hr_Z 3 (by norm_num) (2 * m); simpa using this

lemma g11 : ∀ i, S03.gq (-1) i * S03.gq (-1) i = S03.gq 2 i := fun i => by
  rw [gqm1, gq2, div_mul_div_comm, ← mul_pow]; norm_num; try ring
lemma g12 : ∀ i, S03.gq (-1) i * S03.gq 2 i = S03.gq (-3) i := fun i => by
  rw [gqm1, gq2, gqm3]; ring
lemma g1m2 : ∀ i, S03.gq (-1) i * S03.gq (-2) i = S03.gq 3 i := fun i => by
  rw [gqm1, gqm2, gq3, div_mul_div_comm, ← mul_pow]; norm_num; try ring

lemma t1Q (M : ℕ) : S03.Z [-1, -1] M = (S03.Z [-1] M ^ 2 - S03.Z [2] M) / 2 := by
  have := st2 (-1) (-1) 2 g11 M; linear_combination -this / 2
lemma t2Q (M : ℕ) : S03.Z [2, -1] M = S03.Z [-1] M * S03.Z [2] M - S03.Z [-3] M - S03.Z [-1, 2] M := by
  have := st2 (-1) 2 (-3) g12 M; linear_combination -this
lemma t3Q (M : ℕ) : S03.Z [-2, -1] M = S03.Z [-1] M * S03.Z [-2] M - S03.Z [3] M - S03.Z [-1, -2] M := by
  have := st2 (-1) (-2) 3 g1m2 M; linear_combination -this
lemma t4Q (M : ℕ) : S03.Z [-1, -1, -1] M =
    (S03.Z [-1] M * S03.Z [-1, -1] M - S03.Z [2, -1] M - S03.Z [-1, 2] M) / 3 := by
  have := st3 (-1) (-1) (-1) 2 2 g11 g11 M; linear_combination -this / 3

/-- The `(j, i)` summand of `P_{3,n}` at `𝔷 = 0`, `n = m + 1`. -/
noncomputable def Tji (j i m : ℕ) : ℝ :=
  -(1 / (((m + 1 : ℕ) : ℝ) * (2 * ((m + 1 : ℕ) : ℝ) - 1))) * (V05.ecoef 3 (m + 1) i (3 - i)).coeff j
    * V05.rho j (m + 1) * V05.Dsharp i (3 - i) (m + 1 + j) (m + 1) 0 0

lemma P_split (m : ℕ) : V05.P 3 (m + 1) 0 0 =
    ∑ j ∈ Finset.range 4, ∑ i ∈ Finset.range 4, Tji j i m := by
  unfold V05.P Tji
  rw [Finset.mul_sum]
  refine Finset.sum_congr rfl fun j _ => ?_
  rw [Finset.mul_sum]
  refine Finset.sum_congr rfl fun i _ => ?_
  ring

lemma hE0 (m : ℕ) : hcoef [0, 1] (2 * ((m + 1 : ℕ) : ℤ) - 0) =
      (1 / 4 : ℚ) * (1 / ((m : ℚ) + 1) ^ 2) := by
  rw [show 2 * ((m + 1 : ℕ) : ℤ) - 0 = ((2 * m + 1 : ℕ) : ℤ) + 1 by push_cast; ring, hcoef_succ]
  simp only [ddR0, Znil, S03.Z_cons_succ, gqm1, gq2, gqm2, gq3, gqm3]
  push_cast
  simp only [pow_succ, pow_mul, neg_one_sq, one_pow]
  field_simp
  ring

lemma hO0 (m : ℕ) : hcoef [0, 1] (2 * ((m + 1 : ℕ) : ℤ) + 1 - 0) =
      (1 : ℚ) * (1 / (2 * (m : ℚ) + 3) ^ 2) := by
  rw [show 2 * ((m + 1 : ℕ) : ℤ) + 1 - 0 = ((2 * m + 1 + 1 : ℕ) : ℤ) + 1 by push_cast; ring, hcoef_succ]
  simp only [ddR0, Znil, S03.Z_cons_succ, gqm1, gq2, gqm2, gq3, gqm3]
  push_cast
  simp only [pow_succ, pow_mul, neg_one_sq, one_pow]
  field_simp
  ring

lemma hE1 (m : ℕ) : hcoef [0, -1] (2 * ((m + 1 : ℕ) : ℤ) - 0) =
      (-1 / 4 : ℚ) * (1 / ((m : ℚ) + 1) ^ 2) := by
  rw [show 2 * ((m + 1 : ℕ) : ℤ) - 0 = ((2 * m + 1 : ℕ) : ℤ) + 1 by push_cast; ring, hcoef_succ]
  simp only [ddR1, Znil, S03.Z_cons_succ, gqm1, gq2, gqm2, gq3, gqm3]
  push_cast
  simp only [pow_succ, pow_mul, neg_one_sq, one_pow]
  field_simp
  ring

lemma hO1 (m : ℕ) : hcoef [0, -1] (2 * ((m + 1 : ℕ) : ℤ) + 1 - 0) =
      (1 : ℚ) * (1 / (2 * (m : ℚ) + 3) ^ 2) := by
  rw [show 2 * ((m + 1 : ℕ) : ℤ) + 1 - 0 = ((2 * m + 1 + 1 : ℕ) : ℤ) + 1 by push_cast; ring, hcoef_succ]
  simp only [ddR1, Znil, S03.Z_cons_succ, gqm1, gq2, gqm2, gq3, gqm3]
  push_cast
  simp only [pow_succ, pow_mul, neg_one_sq, one_pow]
  field_simp
  ring

lemma hE2 (m : ℕ) : hcoef [0, 0, 1] (2 * ((m + 1 : ℕ) : ℤ) - 0) =
      (1 / 8 : ℚ) * (1 / ((m : ℚ) + 1) ^ 3) := by
  rw [show 2 * ((m + 1 : ℕ) : ℤ) - 0 = ((2 * m + 1 : ℕ) : ℤ) + 1 by push_cast; ring, hcoef_succ]
  simp only [ddR2, Znil, S03.Z_cons_succ, gqm1, gq2, gqm2, gq3, gqm3]
  push_cast
  simp only [pow_succ, pow_mul, neg_one_sq, one_pow]
  field_simp
  ring

lemma hO2 (m : ℕ) : hcoef [0, 0, 1] (2 * ((m + 1 : ℕ) : ℤ) + 1 - 0) =
      (1 : ℚ) * (1 / (2 * (m : ℚ) + 3) ^ 3) := by
  rw [show 2 * ((m + 1 : ℕ) : ℤ) + 1 - 0 = ((2 * m + 1 + 1 : ℕ) : ℤ) + 1 by push_cast; ring, hcoef_succ]
  simp only [ddR2, Znil, S03.Z_cons_succ, gqm1, gq2, gqm2, gq3, gqm3]
  push_cast
  simp only [pow_succ, pow_mul, neg_one_sq, one_pow]
  field_simp
  ring

lemma hE3 (m : ℕ) : hcoef [0, 0, -1] (2 * ((m + 1 : ℕ) : ℤ) - 0) =
      (-1 / 8 : ℚ) * (1 / ((m : ℚ) + 1) ^ 3) := by
  rw [show 2 * ((m + 1 : ℕ) : ℤ) - 0 = ((2 * m + 1 : ℕ) : ℤ) + 1 by push_cast; ring, hcoef_succ]
  simp only [ddR3, Znil, S03.Z_cons_succ, gqm1, gq2, gqm2, gq3, gqm3]
  push_cast
  simp only [pow_succ, pow_mul, neg_one_sq, one_pow]
  field_simp
  ring

lemma hO3 (m : ℕ) : hcoef [0, 0, -1] (2 * ((m + 1 : ℕ) : ℤ) + 1 - 0) =
      (1 : ℚ) * (1 / (2 * (m : ℚ) + 3) ^ 3) := by
  rw [show 2 * ((m + 1 : ℕ) : ℤ) + 1 - 0 = ((2 * m + 1 + 1 : ℕ) : ℤ) + 1 by push_cast; ring, hcoef_succ]
  simp only [ddR3, Znil, S03.Z_cons_succ, gqm1, gq2, gqm2, gq3, gqm3]
  push_cast
  simp only [pow_succ, pow_mul, neg_one_sq, one_pow]
  field_simp
  ring

lemma hE4 (m : ℕ) : hcoef [0, 1, -1] (2 * ((m + 1 : ℕ) : ℤ) - 0) =
      (-1 / 2 : ℚ) * (1 / ((m : ℚ) + 1) ^ 1) +
      (-1 / 4 : ℚ) * (1 / ((m : ℚ) + 1) ^ 2) +
      (1 : ℚ) * (1 / (2 * (m : ℚ) + 1) ^ 1) +
      (-1 / 4 : ℚ) * (1 / ((m : ℚ) + 1) ^ 2) * S03.Z [-1] (2 * m) := by
  rw [show 2 * ((m + 1 : ℕ) : ℤ) - 0 = ((2 * m + 1 : ℕ) : ℤ) + 1 by push_cast; ring, hcoef_succ]
  simp only [ddR4, Znil, S03.Z_cons_succ, gqm1, gq2, gqm2, gq3, gqm3]
  push_cast
  simp only [pow_succ, pow_mul, neg_one_sq, one_pow]
  field_simp
  ring

lemma hO4 (m : ℕ) : hcoef [0, 1, -1] (2 * ((m + 1 : ℕ) : ℤ) + 1 - 0) =
      (-1 / 2 : ℚ) * (1 / ((m : ℚ) + 1) ^ 1) +
      (1 / 4 : ℚ) * (1 / (2 * (m : ℚ) + 1) ^ 1) +
      (3 / 4 : ℚ) * (1 / (2 * (m : ℚ) + 3) ^ 1) +
      (1 / 2 : ℚ) * (1 / (2 * (m : ℚ) + 3) ^ 2) +
      (-1 : ℚ) * (1 / (2 * (m : ℚ) + 3) ^ 2) * S03.Z [-1] (2 * m) := by
  rw [show 2 * ((m + 1 : ℕ) : ℤ) + 1 - 0 = ((2 * m + 1 + 1 : ℕ) : ℤ) + 1 by push_cast; ring, hcoef_succ]
  simp only [ddR4, Znil, S03.Z_cons_succ, gqm1, gq2, gqm2, gq3, gqm3]
  push_cast
  simp only [pow_succ, pow_mul, neg_one_sq, one_pow]
  field_simp
  ring

lemma hE5 (m : ℕ) : hcoef [0, -1, 1] (2 * ((m + 1 : ℕ) : ℤ) - 0) =
      (-1 / 2 : ℚ) * (1 / ((m : ℚ) + 1) ^ 1) +
      (-1 / 4 : ℚ) * (1 / ((m : ℚ) + 1) ^ 2) +
      (1 : ℚ) * (1 / (2 * (m : ℚ) + 1) ^ 1) +
      (-1 / 4 : ℚ) * (1 / ((m : ℚ) + 1) ^ 2) * S03.Z [-1] (2 * m) := by
  rw [show 2 * ((m + 1 : ℕ) : ℤ) - 0 = ((2 * m + 1 : ℕ) : ℤ) + 1 by push_cast; ring, hcoef_succ]
  simp only [ddR5, Znil, S03.Z_cons_succ, gqm1, gq2, gqm2, gq3, gqm3]
  push_cast
  simp only [pow_succ, pow_mul, neg_one_sq, one_pow]
  field_simp
  ring

lemma hO5 (m : ℕ) : hcoef [0, -1, 1] (2 * ((m + 1 : ℕ) : ℤ) + 1 - 0) =
      (1 / 2 : ℚ) * (1 / ((m : ℚ) + 1) ^ 1) +
      (-1 / 4 : ℚ) * (1 / (2 * (m : ℚ) + 1) ^ 1) +
      (-3 / 4 : ℚ) * (1 / (2 * (m : ℚ) + 3) ^ 1) +
      (-1 / 2 : ℚ) * (1 / (2 * (m : ℚ) + 3) ^ 2) +
      (1 : ℚ) * (1 / (2 * (m : ℚ) + 3) ^ 2) * S03.Z [-1] (2 * m) := by
  rw [show 2 * ((m + 1 : ℕ) : ℤ) + 1 - 0 = ((2 * m + 1 + 1 : ℕ) : ℤ) + 1 by push_cast; ring, hcoef_succ]
  simp only [ddR5, Znil, S03.Z_cons_succ, gqm1, gq2, gqm2, gq3, gqm3]
  push_cast
  simp only [pow_succ, pow_mul, neg_one_sq, one_pow]
  field_simp
  ring

lemma hE6 (m : ℕ) : hcoef [1, 0, 1] (2 * ((m + 1 : ℕ) : ℤ) - -2) =
      (-1 / 8 : ℚ) * (1 / ((m : ℚ) + 1) ^ 1) +
      (1 / 8 : ℚ) * (1 / ((m : ℚ) + 1) ^ 2) +
      (-1 / 9 : ℚ) * (1 / (2 * (m : ℚ) + 1) ^ 1) +
      (1 / 3 : ℚ) * (1 / (2 * (m : ℚ) + 1) ^ 2) +
      (49 / 72 : ℚ) * (1 / ((m : ℚ) + 2) ^ 1) +
      (-1 : ℚ) * (1 / (2 * (m : ℚ) + 3) ^ 1) +
      (1 : ℚ) * (1 / (2 * (m : ℚ) + 3) ^ 2) +
      (1 / 2 : ℚ) * (1 / ((m : ℚ) + 2) ^ 1) * S03.Z [2] (2 * m) := by
  rw [show 2 * ((m + 1 : ℕ) : ℤ) - -2 = ((2 * m + 1 + 1 + 1 : ℕ) : ℤ) + 1 by push_cast; ring, hcoef_succ]
  simp only [ddR6, Znil, S03.Z_cons_succ, gqm1, gq2, gqm2, gq3, gqm3]
  push_cast
  simp only [pow_succ, pow_mul, neg_one_sq, one_pow]
  field_simp
  ring

lemma hO6 (m : ℕ) : hcoef [1, 0, 1] (2 * ((m + 1 : ℕ) : ℤ) + 1 - -2) =
      (-1 / 18 : ℚ) * (1 / ((m : ℚ) + 1) ^ 1) +
      (1 / 12 : ℚ) * (1 / ((m : ℚ) + 1) ^ 2) +
      (-1 / 16 : ℚ) * (1 / (2 * (m : ℚ) + 1) ^ 1) +
      (1 / 4 : ℚ) * (1 / (2 * (m : ℚ) + 1) ^ 2) +
      (-1 / 2 : ℚ) * (1 / ((m : ℚ) + 2) ^ 1) +
      (1 / 4 : ℚ) * (1 / ((m : ℚ) + 2) ^ 2) +
      (-1 / 4 : ℚ) * (1 / (2 * (m : ℚ) + 3) ^ 1) +
      (1 / 2 : ℚ) * (1 / (2 * (m : ℚ) + 3) ^ 2) +
      (205 / 144 : ℚ) * (1 / (2 * (m : ℚ) + 5) ^ 1) +
      (1 : ℚ) * (1 / (2 * (m : ℚ) + 5) ^ 1) * S03.Z [2] (2 * m) := by
  rw [show 2 * ((m + 1 : ℕ) : ℤ) + 1 - -2 = ((2 * m + 1 + 1 + 1 + 1 : ℕ) : ℤ) + 1 by push_cast; ring, hcoef_succ]
  simp only [ddR6, Znil, S03.Z_cons_succ, gqm1, gq2, gqm2, gq3, gqm3]
  push_cast
  simp only [pow_succ, pow_mul, neg_one_sq, one_pow]
  field_simp
  ring

lemma hE7 (m : ℕ) : hcoef [1, 0, 1] (2 * ((m + 1 : ℕ) : ℤ) - -1) =
      (-1 / 2 : ℚ) * (1 / ((m : ℚ) + 1) ^ 1) +
      (1 / 4 : ℚ) * (1 / ((m : ℚ) + 1) ^ 2) +
      (-1 / 4 : ℚ) * (1 / (2 * (m : ℚ) + 1) ^ 1) +
      (1 / 2 : ℚ) * (1 / (2 * (m : ℚ) + 1) ^ 2) +
      (5 / 4 : ℚ) * (1 / (2 * (m : ℚ) + 3) ^ 1) +
      (1 : ℚ) * (1 / (2 * (m : ℚ) + 3) ^ 1) * S03.Z [2] (2 * m) := by
  rw [show 2 * ((m + 1 : ℕ) : ℤ) - -1 = ((2 * m + 1 + 1 : ℕ) : ℤ) + 1 by push_cast; ring, hcoef_succ]
  simp only [ddR6, Znil, S03.Z_cons_succ, gqm1, gq2, gqm2, gq3, gqm3]
  push_cast
  simp only [pow_succ, pow_mul, neg_one_sq, one_pow]
  field_simp
  ring

lemma hO7 (m : ℕ) : hcoef [1, 0, 1] (2 * ((m + 1 : ℕ) : ℤ) + 1 - -1) =
      (-1 / 8 : ℚ) * (1 / ((m : ℚ) + 1) ^ 1) +
      (1 / 8 : ℚ) * (1 / ((m : ℚ) + 1) ^ 2) +
      (-1 / 9 : ℚ) * (1 / (2 * (m : ℚ) + 1) ^ 1) +
      (1 / 3 : ℚ) * (1 / (2 * (m : ℚ) + 1) ^ 2) +
      (49 / 72 : ℚ) * (1 / ((m : ℚ) + 2) ^ 1) +
      (-1 : ℚ) * (1 / (2 * (m : ℚ) + 3) ^ 1) +
      (1 : ℚ) * (1 / (2 * (m : ℚ) + 3) ^ 2) +
      (1 / 2 : ℚ) * (1 / ((m : ℚ) + 2) ^ 1) * S03.Z [2] (2 * m) := by
  rw [show 2 * ((m + 1 : ℕ) : ℤ) + 1 - -1 = ((2 * m + 1 + 1 + 1 : ℕ) : ℤ) + 1 by push_cast; ring, hcoef_succ]
  simp only [ddR6, Znil, S03.Z_cons_succ, gqm1, gq2, gqm2, gq3, gqm3]
  push_cast
  simp only [pow_succ, pow_mul, neg_one_sq, one_pow]
  field_simp
  ring

lemma hE8 (m : ℕ) : hcoef [1, 0, 1] (2 * ((m + 1 : ℕ) : ℤ) - 0) =
      (1 / 2 : ℚ) * (1 / ((m : ℚ) + 1) ^ 1) +
      (-1 : ℚ) * (1 / (2 * (m : ℚ) + 1) ^ 1) +
      (1 : ℚ) * (1 / (2 * (m : ℚ) + 1) ^ 2) +
      (1 / 2 : ℚ) * (1 / ((m : ℚ) + 1) ^ 1) * S03.Z [2] (2 * m) := by
  rw [show 2 * ((m + 1 : ℕ) : ℤ) - 0 = ((2 * m + 1 : ℕ) : ℤ) + 1 by push_cast; ring, hcoef_succ]
  simp only [ddR6, Znil, S03.Z_cons_succ, gqm1, gq2, gqm2, gq3, gqm3]
  push_cast
  simp only [pow_succ, pow_mul, neg_one_sq, one_pow]
  field_simp
  ring

lemma hO8 (m : ℕ) : hcoef [1, 0, 1] (2 * ((m + 1 : ℕ) : ℤ) + 1 - 0) =
      (-1 / 2 : ℚ) * (1 / ((m : ℚ) + 1) ^ 1) +
      (1 / 4 : ℚ) * (1 / ((m : ℚ) + 1) ^ 2) +
      (-1 / 4 : ℚ) * (1 / (2 * (m : ℚ) + 1) ^ 1) +
      (1 / 2 : ℚ) * (1 / (2 * (m : ℚ) + 1) ^ 2) +
      (5 / 4 : ℚ) * (1 / (2 * (m : ℚ) + 3) ^ 1) +
      (1 : ℚ) * (1 / (2 * (m : ℚ) + 3) ^ 1) * S03.Z [2] (2 * m) := by
  rw [show 2 * ((m + 1 : ℕ) : ℤ) + 1 - 0 = ((2 * m + 1 + 1 : ℕ) : ℤ) + 1 by push_cast; ring, hcoef_succ]
  simp only [ddR6, Znil, S03.Z_cons_succ, gqm1, gq2, gqm2, gq3, gqm3]
  push_cast
  simp only [pow_succ, pow_mul, neg_one_sq, one_pow]
  field_simp
  ring

lemma hE9 (m : ℕ) : hcoef [1, 0, -1] (2 * ((m + 1 : ℕ) : ℤ) - -2) =
      (1 / 8 : ℚ) * (1 / ((m : ℚ) + 1) ^ 1) +
      (-1 / 8 : ℚ) * (1 / ((m : ℚ) + 1) ^ 2) +
      (-1 / 9 : ℚ) * (1 / (2 * (m : ℚ) + 1) ^ 1) +
      (1 / 3 : ℚ) * (1 / (2 * (m : ℚ) + 1) ^ 2) +
      (31 / 72 : ℚ) * (1 / ((m : ℚ) + 2) ^ 1) +
      (-1 : ℚ) * (1 / (2 * (m : ℚ) + 3) ^ 1) +
      (1 : ℚ) * (1 / (2 * (m : ℚ) + 3) ^ 2) +
      (-1 / 2 : ℚ) * (1 / ((m : ℚ) + 2) ^ 1) * S03.Z [-2] (2 * m) := by
  rw [show 2 * ((m + 1 : ℕ) : ℤ) - -2 = ((2 * m + 1 + 1 + 1 : ℕ) : ℤ) + 1 by push_cast; ring, hcoef_succ]
  simp only [ddR7, Znil, S03.Z_cons_succ, gqm1, gq2, gqm2, gq3, gqm3]
  push_cast
  simp only [pow_succ, pow_mul, neg_one_sq, one_pow]
  field_simp
  ring

lemma hO9 (m : ℕ) : hcoef [1, 0, -1] (2 * ((m + 1 : ℕ) : ℤ) + 1 - -2) =
      (1 / 18 : ℚ) * (1 / ((m : ℚ) + 1) ^ 1) +
      (-1 / 12 : ℚ) * (1 / ((m : ℚ) + 1) ^ 2) +
      (-1 / 16 : ℚ) * (1 / (2 * (m : ℚ) + 1) ^ 1) +
      (1 / 4 : ℚ) * (1 / (2 * (m : ℚ) + 1) ^ 2) +
      (1 / 2 : ℚ) * (1 / ((m : ℚ) + 2) ^ 1) +
      (-1 / 4 : ℚ) * (1 / ((m : ℚ) + 2) ^ 2) +
      (-1 / 4 : ℚ) * (1 / (2 * (m : ℚ) + 3) ^ 1) +
      (1 / 2 : ℚ) * (1 / (2 * (m : ℚ) + 3) ^ 2) +
      (-115 / 144 : ℚ) * (1 / (2 * (m : ℚ) + 5) ^ 1) +
      (-1 : ℚ) * (1 / (2 * (m : ℚ) + 5) ^ 1) * S03.Z [-2] (2 * m) := by
  rw [show 2 * ((m + 1 : ℕ) : ℤ) + 1 - -2 = ((2 * m + 1 + 1 + 1 + 1 : ℕ) : ℤ) + 1 by push_cast; ring, hcoef_succ]
  simp only [ddR7, Znil, S03.Z_cons_succ, gqm1, gq2, gqm2, gq3, gqm3]
  push_cast
  simp only [pow_succ, pow_mul, neg_one_sq, one_pow]
  field_simp
  ring

lemma hE10 (m : ℕ) : hcoef [1, 0, -1] (2 * ((m + 1 : ℕ) : ℤ) - -1) =
      (1 / 2 : ℚ) * (1 / ((m : ℚ) + 1) ^ 1) +
      (-1 / 4 : ℚ) * (1 / ((m : ℚ) + 1) ^ 2) +
      (-1 / 4 : ℚ) * (1 / (2 * (m : ℚ) + 1) ^ 1) +
      (1 / 2 : ℚ) * (1 / (2 * (m : ℚ) + 1) ^ 2) +
      (-3 / 4 : ℚ) * (1 / (2 * (m : ℚ) + 3) ^ 1) +
      (-1 : ℚ) * (1 / (2 * (m : ℚ) + 3) ^ 1) * S03.Z [-2] (2 * m) := by
  rw [show 2 * ((m + 1 : ℕ) : ℤ) - -1 = ((2 * m + 1 + 1 : ℕ) : ℤ) + 1 by push_cast; ring, hcoef_succ]
  simp only [ddR7, Znil, S03.Z_cons_succ, gqm1, gq2, gqm2, gq3, gqm3]
  push_cast
  simp only [pow_succ, pow_mul, neg_one_sq, one_pow]
  field_simp
  ring

lemma hO10 (m : ℕ) : hcoef [1, 0, -1] (2 * ((m + 1 : ℕ) : ℤ) + 1 - -1) =
      (1 / 8 : ℚ) * (1 / ((m : ℚ) + 1) ^ 1) +
      (-1 / 8 : ℚ) * (1 / ((m : ℚ) + 1) ^ 2) +
      (-1 / 9 : ℚ) * (1 / (2 * (m : ℚ) + 1) ^ 1) +
      (1 / 3 : ℚ) * (1 / (2 * (m : ℚ) + 1) ^ 2) +
      (31 / 72 : ℚ) * (1 / ((m : ℚ) + 2) ^ 1) +
      (-1 : ℚ) * (1 / (2 * (m : ℚ) + 3) ^ 1) +
      (1 : ℚ) * (1 / (2 * (m : ℚ) + 3) ^ 2) +
      (-1 / 2 : ℚ) * (1 / ((m : ℚ) + 2) ^ 1) * S03.Z [-2] (2 * m) := by
  rw [show 2 * ((m + 1 : ℕ) : ℤ) + 1 - -1 = ((2 * m + 1 + 1 + 1 : ℕ) : ℤ) + 1 by push_cast; ring, hcoef_succ]
  simp only [ddR7, Znil, S03.Z_cons_succ, gqm1, gq2, gqm2, gq3, gqm3]
  push_cast
  simp only [pow_succ, pow_mul, neg_one_sq, one_pow]
  field_simp
  ring

lemma hE11 (m : ℕ) : hcoef [1, 0, -1] (2 * ((m + 1 : ℕ) : ℤ) - 0) =
      (1 / 2 : ℚ) * (1 / ((m : ℚ) + 1) ^ 1) +
      (-1 : ℚ) * (1 / (2 * (m : ℚ) + 1) ^ 1) +
      (1 : ℚ) * (1 / (2 * (m : ℚ) + 1) ^ 2) +
      (-1 / 2 : ℚ) * (1 / ((m : ℚ) + 1) ^ 1) * S03.Z [-2] (2 * m) := by
  rw [show 2 * ((m + 1 : ℕ) : ℤ) - 0 = ((2 * m + 1 : ℕ) : ℤ) + 1 by push_cast; ring, hcoef_succ]
  simp only [ddR7, Znil, S03.Z_cons_succ, gqm1, gq2, gqm2, gq3, gqm3]
  push_cast
  simp only [pow_succ, pow_mul, neg_one_sq, one_pow]
  field_simp
  ring

lemma hO11 (m : ℕ) : hcoef [1, 0, -1] (2 * ((m + 1 : ℕ) : ℤ) + 1 - 0) =
      (1 / 2 : ℚ) * (1 / ((m : ℚ) + 1) ^ 1) +
      (-1 / 4 : ℚ) * (1 / ((m : ℚ) + 1) ^ 2) +
      (-1 / 4 : ℚ) * (1 / (2 * (m : ℚ) + 1) ^ 1) +
      (1 / 2 : ℚ) * (1 / (2 * (m : ℚ) + 1) ^ 2) +
      (-3 / 4 : ℚ) * (1 / (2 * (m : ℚ) + 3) ^ 1) +
      (-1 : ℚ) * (1 / (2 * (m : ℚ) + 3) ^ 1) * S03.Z [-2] (2 * m) := by
  rw [show 2 * ((m + 1 : ℕ) : ℤ) + 1 - 0 = ((2 * m + 1 + 1 : ℕ) : ℤ) + 1 by push_cast; ring, hcoef_succ]
  simp only [ddR7, Znil, S03.Z_cons_succ, gqm1, gq2, gqm2, gq3, gqm3]
  push_cast
  simp only [pow_succ, pow_mul, neg_one_sq, one_pow]
  field_simp
  ring

lemma hE12 (m : ℕ) : hcoef [-1, 0, 1] (2 * ((m + 1 : ℕ) : ℤ) - -2) =
      (1 / 8 : ℚ) * (1 / ((m : ℚ) + 1) ^ 1) +
      (-1 / 8 : ℚ) * (1 / ((m : ℚ) + 1) ^ 2) +
      (-1 / 9 : ℚ) * (1 / (2 * (m : ℚ) + 1) ^ 1) +
      (1 / 3 : ℚ) * (1 / (2 * (m : ℚ) + 1) ^ 2) +
      (31 / 72 : ℚ) * (1 / ((m : ℚ) + 2) ^ 1) +
      (-1 : ℚ) * (1 / (2 * (m : ℚ) + 3) ^ 1) +
      (1 : ℚ) * (1 / (2 * (m : ℚ) + 3) ^ 2) +
      (-1 / 2 : ℚ) * (1 / ((m : ℚ) + 2) ^ 1) * S03.Z [-2] (2 * m) := by
  rw [show 2 * ((m + 1 : ℕ) : ℤ) - -2 = ((2 * m + 1 + 1 + 1 : ℕ) : ℤ) + 1 by push_cast; ring, hcoef_succ]
  simp only [ddR8, Znil, S03.Z_cons_succ, gqm1, gq2, gqm2, gq3, gqm3]
  push_cast
  simp only [pow_succ, pow_mul, neg_one_sq, one_pow]
  field_simp
  ring

lemma hO12 (m : ℕ) : hcoef [-1, 0, 1] (2 * ((m + 1 : ℕ) : ℤ) + 1 - -2) =
      (-1 / 18 : ℚ) * (1 / ((m : ℚ) + 1) ^ 1) +
      (1 / 12 : ℚ) * (1 / ((m : ℚ) + 1) ^ 2) +
      (1 / 16 : ℚ) * (1 / (2 * (m : ℚ) + 1) ^ 1) +
      (-1 / 4 : ℚ) * (1 / (2 * (m : ℚ) + 1) ^ 2) +
      (-1 / 2 : ℚ) * (1 / ((m : ℚ) + 2) ^ 1) +
      (1 / 4 : ℚ) * (1 / ((m : ℚ) + 2) ^ 2) +
      (1 / 4 : ℚ) * (1 / (2 * (m : ℚ) + 3) ^ 1) +
      (-1 / 2 : ℚ) * (1 / (2 * (m : ℚ) + 3) ^ 2) +
      (115 / 144 : ℚ) * (1 / (2 * (m : ℚ) + 5) ^ 1) +
      (1 : ℚ) * (1 / (2 * (m : ℚ) + 5) ^ 1) * S03.Z [-2] (2 * m) := by
  rw [show 2 * ((m + 1 : ℕ) : ℤ) + 1 - -2 = ((2 * m + 1 + 1 + 1 + 1 : ℕ) : ℤ) + 1 by push_cast; ring, hcoef_succ]
  simp only [ddR8, Znil, S03.Z_cons_succ, gqm1, gq2, gqm2, gq3, gqm3]
  push_cast
  simp only [pow_succ, pow_mul, neg_one_sq, one_pow]
  field_simp
  ring

lemma hE13 (m : ℕ) : hcoef [-1, 0, 1] (2 * ((m + 1 : ℕ) : ℤ) - -1) =
      (-1 / 2 : ℚ) * (1 / ((m : ℚ) + 1) ^ 1) +
      (1 / 4 : ℚ) * (1 / ((m : ℚ) + 1) ^ 2) +
      (1 / 4 : ℚ) * (1 / (2 * (m : ℚ) + 1) ^ 1) +
      (-1 / 2 : ℚ) * (1 / (2 * (m : ℚ) + 1) ^ 2) +
      (3 / 4 : ℚ) * (1 / (2 * (m : ℚ) + 3) ^ 1) +
      (1 : ℚ) * (1 / (2 * (m : ℚ) + 3) ^ 1) * S03.Z [-2] (2 * m) := by
  rw [show 2 * ((m + 1 : ℕ) : ℤ) - -1 = ((2 * m + 1 + 1 : ℕ) : ℤ) + 1 by push_cast; ring, hcoef_succ]
  simp only [ddR8, Znil, S03.Z_cons_succ, gqm1, gq2, gqm2, gq3, gqm3]
  push_cast
  simp only [pow_succ, pow_mul, neg_one_sq, one_pow]
  field_simp
  ring

lemma hO13 (m : ℕ) : hcoef [-1, 0, 1] (2 * ((m + 1 : ℕ) : ℤ) + 1 - -1) =
      (1 / 8 : ℚ) * (1 / ((m : ℚ) + 1) ^ 1) +
      (-1 / 8 : ℚ) * (1 / ((m : ℚ) + 1) ^ 2) +
      (-1 / 9 : ℚ) * (1 / (2 * (m : ℚ) + 1) ^ 1) +
      (1 / 3 : ℚ) * (1 / (2 * (m : ℚ) + 1) ^ 2) +
      (31 / 72 : ℚ) * (1 / ((m : ℚ) + 2) ^ 1) +
      (-1 : ℚ) * (1 / (2 * (m : ℚ) + 3) ^ 1) +
      (1 : ℚ) * (1 / (2 * (m : ℚ) + 3) ^ 2) +
      (-1 / 2 : ℚ) * (1 / ((m : ℚ) + 2) ^ 1) * S03.Z [-2] (2 * m) := by
  rw [show 2 * ((m + 1 : ℕ) : ℤ) + 1 - -1 = ((2 * m + 1 + 1 + 1 : ℕ) : ℤ) + 1 by push_cast; ring, hcoef_succ]
  simp only [ddR8, Znil, S03.Z_cons_succ, gqm1, gq2, gqm2, gq3, gqm3]
  push_cast
  simp only [pow_succ, pow_mul, neg_one_sq, one_pow]
  field_simp
  ring

lemma hE14 (m : ℕ) : hcoef [-1, 0, 1] (2 * ((m + 1 : ℕ) : ℤ) - 0) =
      (1 / 2 : ℚ) * (1 / ((m : ℚ) + 1) ^ 1) +
      (-1 : ℚ) * (1 / (2 * (m : ℚ) + 1) ^ 1) +
      (1 : ℚ) * (1 / (2 * (m : ℚ) + 1) ^ 2) +
      (-1 / 2 : ℚ) * (1 / ((m : ℚ) + 1) ^ 1) * S03.Z [-2] (2 * m) := by
  rw [show 2 * ((m + 1 : ℕ) : ℤ) - 0 = ((2 * m + 1 : ℕ) : ℤ) + 1 by push_cast; ring, hcoef_succ]
  simp only [ddR8, Znil, S03.Z_cons_succ, gqm1, gq2, gqm2, gq3, gqm3]
  push_cast
  simp only [pow_succ, pow_mul, neg_one_sq, one_pow]
  field_simp
  ring

lemma hO14 (m : ℕ) : hcoef [-1, 0, 1] (2 * ((m + 1 : ℕ) : ℤ) + 1 - 0) =
      (-1 / 2 : ℚ) * (1 / ((m : ℚ) + 1) ^ 1) +
      (1 / 4 : ℚ) * (1 / ((m : ℚ) + 1) ^ 2) +
      (1 / 4 : ℚ) * (1 / (2 * (m : ℚ) + 1) ^ 1) +
      (-1 / 2 : ℚ) * (1 / (2 * (m : ℚ) + 1) ^ 2) +
      (3 / 4 : ℚ) * (1 / (2 * (m : ℚ) + 3) ^ 1) +
      (1 : ℚ) * (1 / (2 * (m : ℚ) + 3) ^ 1) * S03.Z [-2] (2 * m) := by
  rw [show 2 * ((m + 1 : ℕ) : ℤ) + 1 - 0 = ((2 * m + 1 + 1 : ℕ) : ℤ) + 1 by push_cast; ring, hcoef_succ]
  simp only [ddR8, Znil, S03.Z_cons_succ, gqm1, gq2, gqm2, gq3, gqm3]
  push_cast
  simp only [pow_succ, pow_mul, neg_one_sq, one_pow]
  field_simp
  ring

lemma hE15 (m : ℕ) : hcoef [-1, 0, -1] (2 * ((m + 1 : ℕ) : ℤ) - -2) =
      (-1 / 8 : ℚ) * (1 / ((m : ℚ) + 1) ^ 1) +
      (1 / 8 : ℚ) * (1 / ((m : ℚ) + 1) ^ 2) +
      (-1 / 9 : ℚ) * (1 / (2 * (m : ℚ) + 1) ^ 1) +
      (1 / 3 : ℚ) * (1 / (2 * (m : ℚ) + 1) ^ 2) +
      (49 / 72 : ℚ) * (1 / ((m : ℚ) + 2) ^ 1) +
      (-1 : ℚ) * (1 / (2 * (m : ℚ) + 3) ^ 1) +
      (1 : ℚ) * (1 / (2 * (m : ℚ) + 3) ^ 2) +
      (1 / 2 : ℚ) * (1 / ((m : ℚ) + 2) ^ 1) * S03.Z [2] (2 * m) := by
  rw [show 2 * ((m + 1 : ℕ) : ℤ) - -2 = ((2 * m + 1 + 1 + 1 : ℕ) : ℤ) + 1 by push_cast; ring, hcoef_succ]
  simp only [ddR9, Znil, S03.Z_cons_succ, gqm1, gq2, gqm2, gq3, gqm3]
  push_cast
  simp only [pow_succ, pow_mul, neg_one_sq, one_pow]
  field_simp
  ring

lemma hO15 (m : ℕ) : hcoef [-1, 0, -1] (2 * ((m + 1 : ℕ) : ℤ) + 1 - -2) =
      (1 / 18 : ℚ) * (1 / ((m : ℚ) + 1) ^ 1) +
      (-1 / 12 : ℚ) * (1 / ((m : ℚ) + 1) ^ 2) +
      (1 / 16 : ℚ) * (1 / (2 * (m : ℚ) + 1) ^ 1) +
      (-1 / 4 : ℚ) * (1 / (2 * (m : ℚ) + 1) ^ 2) +
      (1 / 2 : ℚ) * (1 / ((m : ℚ) + 2) ^ 1) +
      (-1 / 4 : ℚ) * (1 / ((m : ℚ) + 2) ^ 2) +
      (1 / 4 : ℚ) * (1 / (2 * (m : ℚ) + 3) ^ 1) +
      (-1 / 2 : ℚ) * (1 / (2 * (m : ℚ) + 3) ^ 2) +
      (-205 / 144 : ℚ) * (1 / (2 * (m : ℚ) + 5) ^ 1) +
      (-1 : ℚ) * (1 / (2 * (m : ℚ) + 5) ^ 1) * S03.Z [2] (2 * m) := by
  rw [show 2 * ((m + 1 : ℕ) : ℤ) + 1 - -2 = ((2 * m + 1 + 1 + 1 + 1 : ℕ) : ℤ) + 1 by push_cast; ring, hcoef_succ]
  simp only [ddR9, Znil, S03.Z_cons_succ, gqm1, gq2, gqm2, gq3, gqm3]
  push_cast
  simp only [pow_succ, pow_mul, neg_one_sq, one_pow]
  field_simp
  ring

lemma hE16 (m : ℕ) : hcoef [-1, 0, -1] (2 * ((m + 1 : ℕ) : ℤ) - -1) =
      (1 / 2 : ℚ) * (1 / ((m : ℚ) + 1) ^ 1) +
      (-1 / 4 : ℚ) * (1 / ((m : ℚ) + 1) ^ 2) +
      (1 / 4 : ℚ) * (1 / (2 * (m : ℚ) + 1) ^ 1) +
      (-1 / 2 : ℚ) * (1 / (2 * (m : ℚ) + 1) ^ 2) +
      (-5 / 4 : ℚ) * (1 / (2 * (m : ℚ) + 3) ^ 1) +
      (-1 : ℚ) * (1 / (2 * (m : ℚ) + 3) ^ 1) * S03.Z [2] (2 * m) := by
  rw [show 2 * ((m + 1 : ℕ) : ℤ) - -1 = ((2 * m + 1 + 1 : ℕ) : ℤ) + 1 by push_cast; ring, hcoef_succ]
  simp only [ddR9, Znil, S03.Z_cons_succ, gqm1, gq2, gqm2, gq3, gqm3]
  push_cast
  simp only [pow_succ, pow_mul, neg_one_sq, one_pow]
  field_simp
  ring

lemma hO16 (m : ℕ) : hcoef [-1, 0, -1] (2 * ((m + 1 : ℕ) : ℤ) + 1 - -1) =
      (-1 / 8 : ℚ) * (1 / ((m : ℚ) + 1) ^ 1) +
      (1 / 8 : ℚ) * (1 / ((m : ℚ) + 1) ^ 2) +
      (-1 / 9 : ℚ) * (1 / (2 * (m : ℚ) + 1) ^ 1) +
      (1 / 3 : ℚ) * (1 / (2 * (m : ℚ) + 1) ^ 2) +
      (49 / 72 : ℚ) * (1 / ((m : ℚ) + 2) ^ 1) +
      (-1 : ℚ) * (1 / (2 * (m : ℚ) + 3) ^ 1) +
      (1 : ℚ) * (1 / (2 * (m : ℚ) + 3) ^ 2) +
      (1 / 2 : ℚ) * (1 / ((m : ℚ) + 2) ^ 1) * S03.Z [2] (2 * m) := by
  rw [show 2 * ((m + 1 : ℕ) : ℤ) + 1 - -1 = ((2 * m + 1 + 1 + 1 : ℕ) : ℤ) + 1 by push_cast; ring, hcoef_succ]
  simp only [ddR9, Znil, S03.Z_cons_succ, gqm1, gq2, gqm2, gq3, gqm3]
  push_cast
  simp only [pow_succ, pow_mul, neg_one_sq, one_pow]
  field_simp
  ring

lemma hE17 (m : ℕ) : hcoef [-1, 0, -1] (2 * ((m + 1 : ℕ) : ℤ) - 0) =
      (1 / 2 : ℚ) * (1 / ((m : ℚ) + 1) ^ 1) +
      (-1 : ℚ) * (1 / (2 * (m : ℚ) + 1) ^ 1) +
      (1 : ℚ) * (1 / (2 * (m : ℚ) + 1) ^ 2) +
      (1 / 2 : ℚ) * (1 / ((m : ℚ) + 1) ^ 1) * S03.Z [2] (2 * m) := by
  rw [show 2 * ((m + 1 : ℕ) : ℤ) - 0 = ((2 * m + 1 : ℕ) : ℤ) + 1 by push_cast; ring, hcoef_succ]
  simp only [ddR9, Znil, S03.Z_cons_succ, gqm1, gq2, gqm2, gq3, gqm3]
  push_cast
  simp only [pow_succ, pow_mul, neg_one_sq, one_pow]
  field_simp
  ring

lemma hO17 (m : ℕ) : hcoef [-1, 0, -1] (2 * ((m + 1 : ℕ) : ℤ) + 1 - 0) =
      (1 / 2 : ℚ) * (1 / ((m : ℚ) + 1) ^ 1) +
      (-1 / 4 : ℚ) * (1 / ((m : ℚ) + 1) ^ 2) +
      (1 / 4 : ℚ) * (1 / (2 * (m : ℚ) + 1) ^ 1) +
      (-1 / 2 : ℚ) * (1 / (2 * (m : ℚ) + 1) ^ 2) +
      (-5 / 4 : ℚ) * (1 / (2 * (m : ℚ) + 3) ^ 1) +
      (-1 : ℚ) * (1 / (2 * (m : ℚ) + 3) ^ 1) * S03.Z [2] (2 * m) := by
  rw [show 2 * ((m + 1 : ℕ) : ℤ) + 1 - 0 = ((2 * m + 1 + 1 : ℕ) : ℤ) + 1 by push_cast; ring, hcoef_succ]
  simp only [ddR9, Znil, S03.Z_cons_succ, gqm1, gq2, gqm2, gq3, gqm3]
  push_cast
  simp only [pow_succ, pow_mul, neg_one_sq, one_pow]
  field_simp
  ring

lemma hE18 (m : ℕ) : hcoef [0, 0, 0, 1] (2 * ((m + 1 : ℕ) : ℤ) - 0) =
      (1 / 16 : ℚ) * (1 / ((m : ℚ) + 1) ^ 4) := by
  rw [show 2 * ((m + 1 : ℕ) : ℤ) - 0 = ((2 * m + 1 : ℕ) : ℤ) + 1 by push_cast; ring, hcoef_succ]
  simp only [ddR10, Znil, S03.Z_cons_succ, gqm1, gq2, gqm2, gq3, gqm3]
  push_cast
  simp only [pow_succ, pow_mul, neg_one_sq, one_pow]
  field_simp
  ring

lemma hO18 (m : ℕ) : hcoef [0, 0, 0, 1] (2 * ((m + 1 : ℕ) : ℤ) + 1 - 0) =
      (1 : ℚ) * (1 / (2 * (m : ℚ) + 3) ^ 4) := by
  rw [show 2 * ((m + 1 : ℕ) : ℤ) + 1 - 0 = ((2 * m + 1 + 1 : ℕ) : ℤ) + 1 by push_cast; ring, hcoef_succ]
  simp only [ddR10, Znil, S03.Z_cons_succ, gqm1, gq2, gqm2, gq3, gqm3]
  push_cast
  simp only [pow_succ, pow_mul, neg_one_sq, one_pow]
  field_simp
  ring

lemma hE19 (m : ℕ) : hcoef [0, 0, 0, 1] (2 * ((m + 1 : ℕ) : ℤ) - 1) =
      (1 : ℚ) * (1 / (2 * (m : ℚ) + 1) ^ 4) := by
  rw [show 2 * ((m + 1 : ℕ) : ℤ) - 1 = ((2 * m : ℕ) : ℤ) + 1 by push_cast; ring, hcoef_succ]
  simp only [ddR10, Znil, S03.Z_cons_succ, gqm1, gq2, gqm2, gq3, gqm3]
  push_cast
  simp only [pow_succ, pow_mul, neg_one_sq, one_pow]
  field_simp
  ring

lemma hO19 (m : ℕ) : hcoef [0, 0, 0, 1] (2 * ((m + 1 : ℕ) : ℤ) + 1 - 1) =
      (1 / 16 : ℚ) * (1 / ((m : ℚ) + 1) ^ 4) := by
  rw [show 2 * ((m + 1 : ℕ) : ℤ) + 1 - 1 = ((2 * m + 1 : ℕ) : ℤ) + 1 by push_cast; ring, hcoef_succ]
  simp only [ddR10, Znil, S03.Z_cons_succ, gqm1, gq2, gqm2, gq3, gqm3]
  push_cast
  simp only [pow_succ, pow_mul, neg_one_sq, one_pow]
  field_simp
  ring

lemma hE20 (m : ℕ) : hcoef [0, 0, 0, -1] (2 * ((m + 1 : ℕ) : ℤ) - 0) =
      (-1 / 16 : ℚ) * (1 / ((m : ℚ) + 1) ^ 4) := by
  rw [show 2 * ((m + 1 : ℕ) : ℤ) - 0 = ((2 * m + 1 : ℕ) : ℤ) + 1 by push_cast; ring, hcoef_succ]
  simp only [ddR11, Znil, S03.Z_cons_succ, gqm1, gq2, gqm2, gq3, gqm3]
  push_cast
  simp only [pow_succ, pow_mul, neg_one_sq, one_pow]
  field_simp
  ring

lemma hO20 (m : ℕ) : hcoef [0, 0, 0, -1] (2 * ((m + 1 : ℕ) : ℤ) + 1 - 0) =
      (1 : ℚ) * (1 / (2 * (m : ℚ) + 3) ^ 4) := by
  rw [show 2 * ((m + 1 : ℕ) : ℤ) + 1 - 0 = ((2 * m + 1 + 1 : ℕ) : ℤ) + 1 by push_cast; ring, hcoef_succ]
  simp only [ddR11, Znil, S03.Z_cons_succ, gqm1, gq2, gqm2, gq3, gqm3]
  push_cast
  simp only [pow_succ, pow_mul, neg_one_sq, one_pow]
  field_simp
  ring

lemma hE21 (m : ℕ) : hcoef [0, 0, 0, -1] (2 * ((m + 1 : ℕ) : ℤ) - 1) =
      (1 : ℚ) * (1 / (2 * (m : ℚ) + 1) ^ 4) := by
  rw [show 2 * ((m + 1 : ℕ) : ℤ) - 1 = ((2 * m : ℕ) : ℤ) + 1 by push_cast; ring, hcoef_succ]
  simp only [ddR11, Znil, S03.Z_cons_succ, gqm1, gq2, gqm2, gq3, gqm3]
  push_cast
  simp only [pow_succ, pow_mul, neg_one_sq, one_pow]
  field_simp
  ring

lemma hO21 (m : ℕ) : hcoef [0, 0, 0, -1] (2 * ((m + 1 : ℕ) : ℤ) + 1 - 1) =
      (-1 / 16 : ℚ) * (1 / ((m : ℚ) + 1) ^ 4) := by
  rw [show 2 * ((m + 1 : ℕ) : ℤ) + 1 - 1 = ((2 * m + 1 : ℕ) : ℤ) + 1 by push_cast; ring, hcoef_succ]
  simp only [ddR11, Znil, S03.Z_cons_succ, gqm1, gq2, gqm2, gq3, gqm3]
  push_cast
  simp only [pow_succ, pow_mul, neg_one_sq, one_pow]
  field_simp
  ring

lemma hE22 (m : ℕ) : hcoef [0, 0, 1, -1] (2 * ((m + 1 : ℕ) : ℤ) - 0) =
      (-1 / 2 : ℚ) * (1 / ((m : ℚ) + 1) ^ 1) +
      (-1 / 4 : ℚ) * (1 / ((m : ℚ) + 1) ^ 2) +
      (-1 / 8 : ℚ) * (1 / ((m : ℚ) + 1) ^ 3) +
      (1 : ℚ) * (1 / (2 * (m : ℚ) + 1) ^ 1) +
      (-1 / 8 : ℚ) * (1 / ((m : ℚ) + 1) ^ 3) * S03.Z [-1] (2 * m) := by
  rw [show 2 * ((m + 1 : ℕ) : ℤ) - 0 = ((2 * m + 1 : ℕ) : ℤ) + 1 by push_cast; ring, hcoef_succ]
  simp only [ddR12, Znil, S03.Z_cons_succ, gqm1, gq2, gqm2, gq3, gqm3]
  push_cast
  simp only [pow_succ, pow_mul, neg_one_sq, one_pow]
  field_simp
  ring

lemma hO22 (m : ℕ) : hcoef [0, 0, 1, -1] (2 * ((m + 1 : ℕ) : ℤ) + 1 - 0) =
      (-1 / 2 : ℚ) * (1 / ((m : ℚ) + 1) ^ 1) +
      (1 / 8 : ℚ) * (1 / (2 * (m : ℚ) + 1) ^ 1) +
      (7 / 8 : ℚ) * (1 / (2 * (m : ℚ) + 3) ^ 1) +
      (3 / 4 : ℚ) * (1 / (2 * (m : ℚ) + 3) ^ 2) +
      (1 / 2 : ℚ) * (1 / (2 * (m : ℚ) + 3) ^ 3) +
      (-1 : ℚ) * (1 / (2 * (m : ℚ) + 3) ^ 3) * S03.Z [-1] (2 * m) := by
  rw [show 2 * ((m + 1 : ℕ) : ℤ) + 1 - 0 = ((2 * m + 1 + 1 : ℕ) : ℤ) + 1 by push_cast; ring, hcoef_succ]
  simp only [ddR12, Znil, S03.Z_cons_succ, gqm1, gq2, gqm2, gq3, gqm3]
  push_cast
  simp only [pow_succ, pow_mul, neg_one_sq, one_pow]
  field_simp
  ring

lemma hE23 (m : ℕ) : hcoef [0, 0, 1, -1] (2 * ((m + 1 : ℕ) : ℤ) - 1) =
      (-1 : ℚ) * (1 / (2 * (m : ℚ) + 1) ^ 3) * S03.Z [-1] (2 * m) := by
  rw [show 2 * ((m + 1 : ℕ) : ℤ) - 1 = ((2 * m : ℕ) : ℤ) + 1 by push_cast; ring, hcoef_succ]
  simp only [ddR12, Znil, S03.Z_cons_succ, gqm1, gq2, gqm2, gq3, gqm3]
  push_cast
  simp only [pow_succ, pow_mul, neg_one_sq, one_pow]
  field_simp
  ring

lemma hO23 (m : ℕ) : hcoef [0, 0, 1, -1] (2 * ((m + 1 : ℕ) : ℤ) + 1 - 1) =
      (-1 / 2 : ℚ) * (1 / ((m : ℚ) + 1) ^ 1) +
      (-1 / 4 : ℚ) * (1 / ((m : ℚ) + 1) ^ 2) +
      (-1 / 8 : ℚ) * (1 / ((m : ℚ) + 1) ^ 3) +
      (1 : ℚ) * (1 / (2 * (m : ℚ) + 1) ^ 1) +
      (-1 / 8 : ℚ) * (1 / ((m : ℚ) + 1) ^ 3) * S03.Z [-1] (2 * m) := by
  rw [show 2 * ((m + 1 : ℕ) : ℤ) + 1 - 1 = ((2 * m + 1 : ℕ) : ℤ) + 1 by push_cast; ring, hcoef_succ]
  simp only [ddR12, Znil, S03.Z_cons_succ, gqm1, gq2, gqm2, gq3, gqm3]
  push_cast
  simp only [pow_succ, pow_mul, neg_one_sq, one_pow]
  field_simp
  ring

lemma hE24 (m : ℕ) : hcoef [0, 0, -1, 1] (2 * ((m + 1 : ℕ) : ℤ) - 0) =
      (-1 / 2 : ℚ) * (1 / ((m : ℚ) + 1) ^ 1) +
      (-1 / 4 : ℚ) * (1 / ((m : ℚ) + 1) ^ 2) +
      (-1 / 8 : ℚ) * (1 / ((m : ℚ) + 1) ^ 3) +
      (1 : ℚ) * (1 / (2 * (m : ℚ) + 1) ^ 1) +
      (-1 / 8 : ℚ) * (1 / ((m : ℚ) + 1) ^ 3) * S03.Z [-1] (2 * m) := by
  rw [show 2 * ((m + 1 : ℕ) : ℤ) - 0 = ((2 * m + 1 : ℕ) : ℤ) + 1 by push_cast; ring, hcoef_succ]
  simp only [ddR13, Znil, S03.Z_cons_succ, gqm1, gq2, gqm2, gq3, gqm3]
  push_cast
  simp only [pow_succ, pow_mul, neg_one_sq, one_pow]
  field_simp
  ring

lemma hO24 (m : ℕ) : hcoef [0, 0, -1, 1] (2 * ((m + 1 : ℕ) : ℤ) + 1 - 0) =
      (1 / 2 : ℚ) * (1 / ((m : ℚ) + 1) ^ 1) +
      (-1 / 8 : ℚ) * (1 / (2 * (m : ℚ) + 1) ^ 1) +
      (-7 / 8 : ℚ) * (1 / (2 * (m : ℚ) + 3) ^ 1) +
      (-3 / 4 : ℚ) * (1 / (2 * (m : ℚ) + 3) ^ 2) +
      (-1 / 2 : ℚ) * (1 / (2 * (m : ℚ) + 3) ^ 3) +
      (1 : ℚ) * (1 / (2 * (m : ℚ) + 3) ^ 3) * S03.Z [-1] (2 * m) := by
  rw [show 2 * ((m + 1 : ℕ) : ℤ) + 1 - 0 = ((2 * m + 1 + 1 : ℕ) : ℤ) + 1 by push_cast; ring, hcoef_succ]
  simp only [ddR13, Znil, S03.Z_cons_succ, gqm1, gq2, gqm2, gq3, gqm3]
  push_cast
  simp only [pow_succ, pow_mul, neg_one_sq, one_pow]
  field_simp
  ring

lemma hE25 (m : ℕ) : hcoef [0, 0, -1, 1] (2 * ((m + 1 : ℕ) : ℤ) - 1) =
      (1 : ℚ) * (1 / (2 * (m : ℚ) + 1) ^ 3) * S03.Z [-1] (2 * m) := by
  rw [show 2 * ((m + 1 : ℕ) : ℤ) - 1 = ((2 * m : ℕ) : ℤ) + 1 by push_cast; ring, hcoef_succ]
  simp only [ddR13, Znil, S03.Z_cons_succ, gqm1, gq2, gqm2, gq3, gqm3]
  push_cast
  simp only [pow_succ, pow_mul, neg_one_sq, one_pow]
  field_simp
  ring

lemma hO25 (m : ℕ) : hcoef [0, 0, -1, 1] (2 * ((m + 1 : ℕ) : ℤ) + 1 - 1) =
      (-1 / 2 : ℚ) * (1 / ((m : ℚ) + 1) ^ 1) +
      (-1 / 4 : ℚ) * (1 / ((m : ℚ) + 1) ^ 2) +
      (-1 / 8 : ℚ) * (1 / ((m : ℚ) + 1) ^ 3) +
      (1 : ℚ) * (1 / (2 * (m : ℚ) + 1) ^ 1) +
      (-1 / 8 : ℚ) * (1 / ((m : ℚ) + 1) ^ 3) * S03.Z [-1] (2 * m) := by
  rw [show 2 * ((m + 1 : ℕ) : ℤ) + 1 - 1 = ((2 * m + 1 : ℕ) : ℤ) + 1 by push_cast; ring, hcoef_succ]
  simp only [ddR13, Znil, S03.Z_cons_succ, gqm1, gq2, gqm2, gq3, gqm3]
  push_cast
  simp only [pow_succ, pow_mul, neg_one_sq, one_pow]
  field_simp
  ring

lemma hE26 (m : ℕ) : hcoef [0, 1, 0, 1] (2 * ((m + 1 : ℕ) : ℤ) - -1) =
      (-1 : ℚ) * (1 / ((m : ℚ) + 1) ^ 1) +
      (1 / 4 : ℚ) * (1 / ((m : ℚ) + 1) ^ 2) +
      (-1 / 4 : ℚ) * (1 / (2 * (m : ℚ) + 1) ^ 1) +
      (1 / 4 : ℚ) * (1 / (2 * (m : ℚ) + 1) ^ 2) +
      (9 / 4 : ℚ) * (1 / (2 * (m : ℚ) + 3) ^ 1) +
      (5 / 4 : ℚ) * (1 / (2 * (m : ℚ) + 3) ^ 2) +
      (1 : ℚ) * (1 / (2 * (m : ℚ) + 3) ^ 2) * S03.Z [2] (2 * m) := by
  rw [show 2 * ((m + 1 : ℕ) : ℤ) - -1 = ((2 * m + 1 + 1 : ℕ) : ℤ) + 1 by push_cast; ring, hcoef_succ]
  simp only [ddR14, Znil, S03.Z_cons_succ, gqm1, gq2, gqm2, gq3, gqm3]
  push_cast
  simp only [pow_succ, pow_mul, neg_one_sq, one_pow]
  field_simp
  ring

lemma hO26 (m : ℕ) : hcoef [0, 1, 0, 1] (2 * ((m + 1 : ℕ) : ℤ) + 1 - -1) =
      (-1 / 8 : ℚ) * (1 / ((m : ℚ) + 1) ^ 1) +
      (1 / 16 : ℚ) * (1 / ((m : ℚ) + 1) ^ 2) +
      (-2 / 27 : ℚ) * (1 / (2 * (m : ℚ) + 1) ^ 1) +
      (1 / 9 : ℚ) * (1 / (2 * (m : ℚ) + 1) ^ 2) +
      (251 / 216 : ℚ) * (1 / ((m : ℚ) + 2) ^ 1) +
      (49 / 144 : ℚ) * (1 / ((m : ℚ) + 2) ^ 2) +
      (-2 : ℚ) * (1 / (2 * (m : ℚ) + 3) ^ 1) +
      (1 : ℚ) * (1 / (2 * (m : ℚ) + 3) ^ 2) +
      (1 / 4 : ℚ) * (1 / ((m : ℚ) + 2) ^ 2) * S03.Z [2] (2 * m) := by
  rw [show 2 * ((m + 1 : ℕ) : ℤ) + 1 - -1 = ((2 * m + 1 + 1 + 1 : ℕ) : ℤ) + 1 by push_cast; ring, hcoef_succ]
  simp only [ddR14, Znil, S03.Z_cons_succ, gqm1, gq2, gqm2, gq3, gqm3]
  push_cast
  simp only [pow_succ, pow_mul, neg_one_sq, one_pow]
  field_simp
  ring

lemma hE27 (m : ℕ) : hcoef [0, 1, 0, 1] (2 * ((m + 1 : ℕ) : ℤ) - 1) =
      (1 : ℚ) * (1 / (2 * (m : ℚ) + 1) ^ 2) * S03.Z [2] (2 * m) := by
  rw [show 2 * ((m + 1 : ℕ) : ℤ) - 1 = ((2 * m : ℕ) : ℤ) + 1 by push_cast; ring, hcoef_succ]
  simp only [ddR14, Znil, S03.Z_cons_succ, gqm1, gq2, gqm2, gq3, gqm3]
  push_cast
  simp only [pow_succ, pow_mul, neg_one_sq, one_pow]
  field_simp
  ring

lemma hO27 (m : ℕ) : hcoef [0, 1, 0, 1] (2 * ((m + 1 : ℕ) : ℤ) + 1 - 1) =
      (1 : ℚ) * (1 / ((m : ℚ) + 1) ^ 1) +
      (1 / 4 : ℚ) * (1 / ((m : ℚ) + 1) ^ 2) +
      (-2 : ℚ) * (1 / (2 * (m : ℚ) + 1) ^ 1) +
      (1 : ℚ) * (1 / (2 * (m : ℚ) + 1) ^ 2) +
      (1 / 4 : ℚ) * (1 / ((m : ℚ) + 1) ^ 2) * S03.Z [2] (2 * m) := by
  rw [show 2 * ((m + 1 : ℕ) : ℤ) + 1 - 1 = ((2 * m + 1 : ℕ) : ℤ) + 1 by push_cast; ring, hcoef_succ]
  simp only [ddR14, Znil, S03.Z_cons_succ, gqm1, gq2, gqm2, gq3, gqm3]
  push_cast
  simp only [pow_succ, pow_mul, neg_one_sq, one_pow]
  field_simp
  ring

lemma hE28 (m : ℕ) : hcoef [0, 1, 0, -1] (2 * ((m + 1 : ℕ) : ℤ) - -1) =
      (1 : ℚ) * (1 / ((m : ℚ) + 1) ^ 1) +
      (-1 / 4 : ℚ) * (1 / ((m : ℚ) + 1) ^ 2) +
      (-1 / 4 : ℚ) * (1 / (2 * (m : ℚ) + 1) ^ 1) +
      (1 / 4 : ℚ) * (1 / (2 * (m : ℚ) + 1) ^ 2) +
      (-7 / 4 : ℚ) * (1 / (2 * (m : ℚ) + 3) ^ 1) +
      (-3 / 4 : ℚ) * (1 / (2 * (m : ℚ) + 3) ^ 2) +
      (-1 : ℚ) * (1 / (2 * (m : ℚ) + 3) ^ 2) * S03.Z [-2] (2 * m) := by
  rw [show 2 * ((m + 1 : ℕ) : ℤ) - -1 = ((2 * m + 1 + 1 : ℕ) : ℤ) + 1 by push_cast; ring, hcoef_succ]
  simp only [ddR15, Znil, S03.Z_cons_succ, gqm1, gq2, gqm2, gq3, gqm3]
  push_cast
  simp only [pow_succ, pow_mul, neg_one_sq, one_pow]
  field_simp
  ring

lemma hO28 (m : ℕ) : hcoef [0, 1, 0, -1] (2 * ((m + 1 : ℕ) : ℤ) + 1 - -1) =
      (1 / 8 : ℚ) * (1 / ((m : ℚ) + 1) ^ 1) +
      (-1 / 16 : ℚ) * (1 / ((m : ℚ) + 1) ^ 2) +
      (-2 / 27 : ℚ) * (1 / (2 * (m : ℚ) + 1) ^ 1) +
      (1 / 9 : ℚ) * (1 / (2 * (m : ℚ) + 1) ^ 2) +
      (197 / 216 : ℚ) * (1 / ((m : ℚ) + 2) ^ 1) +
      (31 / 144 : ℚ) * (1 / ((m : ℚ) + 2) ^ 2) +
      (-2 : ℚ) * (1 / (2 * (m : ℚ) + 3) ^ 1) +
      (1 : ℚ) * (1 / (2 * (m : ℚ) + 3) ^ 2) +
      (-1 / 4 : ℚ) * (1 / ((m : ℚ) + 2) ^ 2) * S03.Z [-2] (2 * m) := by
  rw [show 2 * ((m + 1 : ℕ) : ℤ) + 1 - -1 = ((2 * m + 1 + 1 + 1 : ℕ) : ℤ) + 1 by push_cast; ring, hcoef_succ]
  simp only [ddR15, Znil, S03.Z_cons_succ, gqm1, gq2, gqm2, gq3, gqm3]
  push_cast
  simp only [pow_succ, pow_mul, neg_one_sq, one_pow]
  field_simp
  ring

lemma hE29 (m : ℕ) : hcoef [0, 1, 0, -1] (2 * ((m + 1 : ℕ) : ℤ) - 0) =
      (1 : ℚ) * (1 / ((m : ℚ) + 1) ^ 1) +
      (1 / 4 : ℚ) * (1 / ((m : ℚ) + 1) ^ 2) +
      (-2 : ℚ) * (1 / (2 * (m : ℚ) + 1) ^ 1) +
      (1 : ℚ) * (1 / (2 * (m : ℚ) + 1) ^ 2) +
      (-1 / 4 : ℚ) * (1 / ((m : ℚ) + 1) ^ 2) * S03.Z [-2] (2 * m) := by
  rw [show 2 * ((m + 1 : ℕ) : ℤ) - 0 = ((2 * m + 1 : ℕ) : ℤ) + 1 by push_cast; ring, hcoef_succ]
  simp only [ddR15, Znil, S03.Z_cons_succ, gqm1, gq2, gqm2, gq3, gqm3]
  push_cast
  simp only [pow_succ, pow_mul, neg_one_sq, one_pow]
  field_simp
  ring

lemma hO29 (m : ℕ) : hcoef [0, 1, 0, -1] (2 * ((m + 1 : ℕ) : ℤ) + 1 - 0) =
      (1 : ℚ) * (1 / ((m : ℚ) + 1) ^ 1) +
      (-1 / 4 : ℚ) * (1 / ((m : ℚ) + 1) ^ 2) +
      (-1 / 4 : ℚ) * (1 / (2 * (m : ℚ) + 1) ^ 1) +
      (1 / 4 : ℚ) * (1 / (2 * (m : ℚ) + 1) ^ 2) +
      (-7 / 4 : ℚ) * (1 / (2 * (m : ℚ) + 3) ^ 1) +
      (-3 / 4 : ℚ) * (1 / (2 * (m : ℚ) + 3) ^ 2) +
      (-1 : ℚ) * (1 / (2 * (m : ℚ) + 3) ^ 2) * S03.Z [-2] (2 * m) := by
  rw [show 2 * ((m + 1 : ℕ) : ℤ) + 1 - 0 = ((2 * m + 1 + 1 : ℕ) : ℤ) + 1 by push_cast; ring, hcoef_succ]
  simp only [ddR15, Znil, S03.Z_cons_succ, gqm1, gq2, gqm2, gq3, gqm3]
  push_cast
  simp only [pow_succ, pow_mul, neg_one_sq, one_pow]
  field_simp
  ring

lemma hE30 (m : ℕ) : hcoef [0, 1, 0, -1] (2 * ((m + 1 : ℕ) : ℤ) - 1) =
      (-1 : ℚ) * (1 / (2 * (m : ℚ) + 1) ^ 2) * S03.Z [-2] (2 * m) := by
  rw [show 2 * ((m + 1 : ℕ) : ℤ) - 1 = ((2 * m : ℕ) : ℤ) + 1 by push_cast; ring, hcoef_succ]
  simp only [ddR15, Znil, S03.Z_cons_succ, gqm1, gq2, gqm2, gq3, gqm3]
  push_cast
  simp only [pow_succ, pow_mul, neg_one_sq, one_pow]
  field_simp
  ring

lemma hO30 (m : ℕ) : hcoef [0, 1, 0, -1] (2 * ((m + 1 : ℕ) : ℤ) + 1 - 1) =
      (1 : ℚ) * (1 / ((m : ℚ) + 1) ^ 1) +
      (1 / 4 : ℚ) * (1 / ((m : ℚ) + 1) ^ 2) +
      (-2 : ℚ) * (1 / (2 * (m : ℚ) + 1) ^ 1) +
      (1 : ℚ) * (1 / (2 * (m : ℚ) + 1) ^ 2) +
      (-1 / 4 : ℚ) * (1 / ((m : ℚ) + 1) ^ 2) * S03.Z [-2] (2 * m) := by
  rw [show 2 * ((m + 1 : ℕ) : ℤ) + 1 - 1 = ((2 * m + 1 : ℕ) : ℤ) + 1 by push_cast; ring, hcoef_succ]
  simp only [ddR15, Znil, S03.Z_cons_succ, gqm1, gq2, gqm2, gq3, gqm3]
  push_cast
  simp only [pow_succ, pow_mul, neg_one_sq, one_pow]
  field_simp
  ring

lemma hE31 (m : ℕ) : hcoef [0, 1, -1, 1] (2 * ((m + 1 : ℕ) : ℤ) - 0) =
      (-1 / 2 : ℚ) * (1 / ((m : ℚ) + 1) ^ 1) * S03.Z [-1] (2 * m) +
      (-1 / 4 : ℚ) * (1 / ((m : ℚ) + 1) ^ 2) * S03.Z [-1] (2 * m) +
      (1 : ℚ) * (1 / (2 * (m : ℚ) + 1) ^ 1) * S03.Z [-1] (2 * m) +
      (-1 / 4 : ℚ) * (1 / ((m : ℚ) + 1) ^ 2) * S03.Z [-1, -1] (2 * m) := by
  rw [show 2 * ((m + 1 : ℕ) : ℤ) - 0 = ((2 * m + 1 : ℕ) : ℤ) + 1 by push_cast; ring, hcoef_succ]
  simp only [ddR16, Znil, S03.Z_cons_succ, gqm1, gq2, gqm2, gq3, gqm3]
  push_cast
  simp only [pow_succ, pow_mul, neg_one_sq, one_pow]
  field_simp
  ring

lemma hO31 (m : ℕ) : hcoef [0, 1, -1, 1] (2 * ((m + 1 : ℕ) : ℤ) + 1 - 0) =
      (-1 / 2 : ℚ) * (1 / ((m : ℚ) + 1) ^ 1) +
      (1 / 4 : ℚ) * (1 / (2 * (m : ℚ) + 1) ^ 1) +
      (3 / 4 : ℚ) * (1 / (2 * (m : ℚ) + 3) ^ 1) +
      (1 / 2 : ℚ) * (1 / (2 * (m : ℚ) + 3) ^ 2) +
      (-1 / 2 : ℚ) * (1 / ((m : ℚ) + 1) ^ 1) * S03.Z [-1] (2 * m) +
      (1 / 4 : ℚ) * (1 / (2 * (m : ℚ) + 1) ^ 1) * S03.Z [-1] (2 * m) +
      (3 / 4 : ℚ) * (1 / (2 * (m : ℚ) + 3) ^ 1) * S03.Z [-1] (2 * m) +
      (1 / 2 : ℚ) * (1 / (2 * (m : ℚ) + 3) ^ 2) * S03.Z [-1] (2 * m) +
      (-1 : ℚ) * (1 / (2 * (m : ℚ) + 3) ^ 2) * S03.Z [-1, -1] (2 * m) := by
  rw [show 2 * ((m + 1 : ℕ) : ℤ) + 1 - 0 = ((2 * m + 1 + 1 : ℕ) : ℤ) + 1 by push_cast; ring, hcoef_succ]
  simp only [ddR16, Znil, S03.Z_cons_succ, gqm1, gq2, gqm2, gq3, gqm3]
  push_cast
  simp only [pow_succ, pow_mul, neg_one_sq, one_pow]
  field_simp
  ring

lemma hE32 (m : ℕ) : hcoef [0, 1, -1, 1] (2 * ((m + 1 : ℕ) : ℤ) - 1) =
      (-1 : ℚ) * (1 / (2 * (m : ℚ) + 1) ^ 2) * S03.Z [-1, -1] (2 * m) := by
  rw [show 2 * ((m + 1 : ℕ) : ℤ) - 1 = ((2 * m : ℕ) : ℤ) + 1 by push_cast; ring, hcoef_succ]
  simp only [ddR16, Znil, S03.Z_cons_succ, gqm1, gq2, gqm2, gq3, gqm3]
  push_cast
  simp only [pow_succ, pow_mul, neg_one_sq, one_pow]
  field_simp
  ring

lemma hO32 (m : ℕ) : hcoef [0, 1, -1, 1] (2 * ((m + 1 : ℕ) : ℤ) + 1 - 1) =
      (-1 / 2 : ℚ) * (1 / ((m : ℚ) + 1) ^ 1) * S03.Z [-1] (2 * m) +
      (-1 / 4 : ℚ) * (1 / ((m : ℚ) + 1) ^ 2) * S03.Z [-1] (2 * m) +
      (1 : ℚ) * (1 / (2 * (m : ℚ) + 1) ^ 1) * S03.Z [-1] (2 * m) +
      (-1 / 4 : ℚ) * (1 / ((m : ℚ) + 1) ^ 2) * S03.Z [-1, -1] (2 * m) := by
  rw [show 2 * ((m + 1 : ℕ) : ℤ) + 1 - 1 = ((2 * m + 1 : ℕ) : ℤ) + 1 by push_cast; ring, hcoef_succ]
  simp only [ddR16, Znil, S03.Z_cons_succ, gqm1, gq2, gqm2, gq3, gqm3]
  push_cast
  simp only [pow_succ, pow_mul, neg_one_sq, one_pow]
  field_simp
  ring

lemma hE33 (m : ℕ) : hcoef [0, -1, 0, 1] (2 * ((m + 1 : ℕ) : ℤ) - -1) =
      (-1 : ℚ) * (1 / ((m : ℚ) + 1) ^ 1) +
      (1 / 4 : ℚ) * (1 / ((m : ℚ) + 1) ^ 2) +
      (1 / 4 : ℚ) * (1 / (2 * (m : ℚ) + 1) ^ 1) +
      (-1 / 4 : ℚ) * (1 / (2 * (m : ℚ) + 1) ^ 2) +
      (7 / 4 : ℚ) * (1 / (2 * (m : ℚ) + 3) ^ 1) +
      (3 / 4 : ℚ) * (1 / (2 * (m : ℚ) + 3) ^ 2) +
      (1 : ℚ) * (1 / (2 * (m : ℚ) + 3) ^ 2) * S03.Z [-2] (2 * m) := by
  rw [show 2 * ((m + 1 : ℕ) : ℤ) - -1 = ((2 * m + 1 + 1 : ℕ) : ℤ) + 1 by push_cast; ring, hcoef_succ]
  simp only [ddR17, Znil, S03.Z_cons_succ, gqm1, gq2, gqm2, gq3, gqm3]
  push_cast
  simp only [pow_succ, pow_mul, neg_one_sq, one_pow]
  field_simp
  ring

lemma hO33 (m : ℕ) : hcoef [0, -1, 0, 1] (2 * ((m + 1 : ℕ) : ℤ) + 1 - -1) =
      (1 / 8 : ℚ) * (1 / ((m : ℚ) + 1) ^ 1) +
      (-1 / 16 : ℚ) * (1 / ((m : ℚ) + 1) ^ 2) +
      (-2 / 27 : ℚ) * (1 / (2 * (m : ℚ) + 1) ^ 1) +
      (1 / 9 : ℚ) * (1 / (2 * (m : ℚ) + 1) ^ 2) +
      (197 / 216 : ℚ) * (1 / ((m : ℚ) + 2) ^ 1) +
      (31 / 144 : ℚ) * (1 / ((m : ℚ) + 2) ^ 2) +
      (-2 : ℚ) * (1 / (2 * (m : ℚ) + 3) ^ 1) +
      (1 : ℚ) * (1 / (2 * (m : ℚ) + 3) ^ 2) +
      (-1 / 4 : ℚ) * (1 / ((m : ℚ) + 2) ^ 2) * S03.Z [-2] (2 * m) := by
  rw [show 2 * ((m + 1 : ℕ) : ℤ) + 1 - -1 = ((2 * m + 1 + 1 + 1 : ℕ) : ℤ) + 1 by push_cast; ring, hcoef_succ]
  simp only [ddR17, Znil, S03.Z_cons_succ, gqm1, gq2, gqm2, gq3, gqm3]
  push_cast
  simp only [pow_succ, pow_mul, neg_one_sq, one_pow]
  field_simp
  ring

lemma hE34 (m : ℕ) : hcoef [0, -1, 0, 1] (2 * ((m + 1 : ℕ) : ℤ) - 0) =
      (1 : ℚ) * (1 / ((m : ℚ) + 1) ^ 1) +
      (1 / 4 : ℚ) * (1 / ((m : ℚ) + 1) ^ 2) +
      (-2 : ℚ) * (1 / (2 * (m : ℚ) + 1) ^ 1) +
      (1 : ℚ) * (1 / (2 * (m : ℚ) + 1) ^ 2) +
      (-1 / 4 : ℚ) * (1 / ((m : ℚ) + 1) ^ 2) * S03.Z [-2] (2 * m) := by
  rw [show 2 * ((m + 1 : ℕ) : ℤ) - 0 = ((2 * m + 1 : ℕ) : ℤ) + 1 by push_cast; ring, hcoef_succ]
  simp only [ddR17, Znil, S03.Z_cons_succ, gqm1, gq2, gqm2, gq3, gqm3]
  push_cast
  simp only [pow_succ, pow_mul, neg_one_sq, one_pow]
  field_simp
  ring

lemma hO34 (m : ℕ) : hcoef [0, -1, 0, 1] (2 * ((m + 1 : ℕ) : ℤ) + 1 - 0) =
      (-1 : ℚ) * (1 / ((m : ℚ) + 1) ^ 1) +
      (1 / 4 : ℚ) * (1 / ((m : ℚ) + 1) ^ 2) +
      (1 / 4 : ℚ) * (1 / (2 * (m : ℚ) + 1) ^ 1) +
      (-1 / 4 : ℚ) * (1 / (2 * (m : ℚ) + 1) ^ 2) +
      (7 / 4 : ℚ) * (1 / (2 * (m : ℚ) + 3) ^ 1) +
      (3 / 4 : ℚ) * (1 / (2 * (m : ℚ) + 3) ^ 2) +
      (1 : ℚ) * (1 / (2 * (m : ℚ) + 3) ^ 2) * S03.Z [-2] (2 * m) := by
  rw [show 2 * ((m + 1 : ℕ) : ℤ) + 1 - 0 = ((2 * m + 1 + 1 : ℕ) : ℤ) + 1 by push_cast; ring, hcoef_succ]
  simp only [ddR17, Znil, S03.Z_cons_succ, gqm1, gq2, gqm2, gq3, gqm3]
  push_cast
  simp only [pow_succ, pow_mul, neg_one_sq, one_pow]
  field_simp
  ring

lemma hE35 (m : ℕ) : hcoef [0, -1, 0, 1] (2 * ((m + 1 : ℕ) : ℤ) - 1) =
      (1 : ℚ) * (1 / (2 * (m : ℚ) + 1) ^ 2) * S03.Z [-2] (2 * m) := by
  rw [show 2 * ((m + 1 : ℕ) : ℤ) - 1 = ((2 * m : ℕ) : ℤ) + 1 by push_cast; ring, hcoef_succ]
  simp only [ddR17, Znil, S03.Z_cons_succ, gqm1, gq2, gqm2, gq3, gqm3]
  push_cast
  simp only [pow_succ, pow_mul, neg_one_sq, one_pow]
  field_simp
  ring

lemma hO35 (m : ℕ) : hcoef [0, -1, 0, 1] (2 * ((m + 1 : ℕ) : ℤ) + 1 - 1) =
      (1 : ℚ) * (1 / ((m : ℚ) + 1) ^ 1) +
      (1 / 4 : ℚ) * (1 / ((m : ℚ) + 1) ^ 2) +
      (-2 : ℚ) * (1 / (2 * (m : ℚ) + 1) ^ 1) +
      (1 : ℚ) * (1 / (2 * (m : ℚ) + 1) ^ 2) +
      (-1 / 4 : ℚ) * (1 / ((m : ℚ) + 1) ^ 2) * S03.Z [-2] (2 * m) := by
  rw [show 2 * ((m + 1 : ℕ) : ℤ) + 1 - 1 = ((2 * m + 1 : ℕ) : ℤ) + 1 by push_cast; ring, hcoef_succ]
  simp only [ddR17, Znil, S03.Z_cons_succ, gqm1, gq2, gqm2, gq3, gqm3]
  push_cast
  simp only [pow_succ, pow_mul, neg_one_sq, one_pow]
  field_simp
  ring

lemma hE36 (m : ℕ) : hcoef [0, -1, 0, -1] (2 * ((m + 1 : ℕ) : ℤ) - -1) =
      (1 : ℚ) * (1 / ((m : ℚ) + 1) ^ 1) +
      (-1 / 4 : ℚ) * (1 / ((m : ℚ) + 1) ^ 2) +
      (1 / 4 : ℚ) * (1 / (2 * (m : ℚ) + 1) ^ 1) +
      (-1 / 4 : ℚ) * (1 / (2 * (m : ℚ) + 1) ^ 2) +
      (-9 / 4 : ℚ) * (1 / (2 * (m : ℚ) + 3) ^ 1) +
      (-5 / 4 : ℚ) * (1 / (2 * (m : ℚ) + 3) ^ 2) +
      (-1 : ℚ) * (1 / (2 * (m : ℚ) + 3) ^ 2) * S03.Z [2] (2 * m) := by
  rw [show 2 * ((m + 1 : ℕ) : ℤ) - -1 = ((2 * m + 1 + 1 : ℕ) : ℤ) + 1 by push_cast; ring, hcoef_succ]
  simp only [ddR18, Znil, S03.Z_cons_succ, gqm1, gq2, gqm2, gq3, gqm3]
  push_cast
  simp only [pow_succ, pow_mul, neg_one_sq, one_pow]
  field_simp
  ring

lemma hO36 (m : ℕ) : hcoef [0, -1, 0, -1] (2 * ((m + 1 : ℕ) : ℤ) + 1 - -1) =
      (-1 / 8 : ℚ) * (1 / ((m : ℚ) + 1) ^ 1) +
      (1 / 16 : ℚ) * (1 / ((m : ℚ) + 1) ^ 2) +
      (-2 / 27 : ℚ) * (1 / (2 * (m : ℚ) + 1) ^ 1) +
      (1 / 9 : ℚ) * (1 / (2 * (m : ℚ) + 1) ^ 2) +
      (251 / 216 : ℚ) * (1 / ((m : ℚ) + 2) ^ 1) +
      (49 / 144 : ℚ) * (1 / ((m : ℚ) + 2) ^ 2) +
      (-2 : ℚ) * (1 / (2 * (m : ℚ) + 3) ^ 1) +
      (1 : ℚ) * (1 / (2 * (m : ℚ) + 3) ^ 2) +
      (1 / 4 : ℚ) * (1 / ((m : ℚ) + 2) ^ 2) * S03.Z [2] (2 * m) := by
  rw [show 2 * ((m + 1 : ℕ) : ℤ) + 1 - -1 = ((2 * m + 1 + 1 + 1 : ℕ) : ℤ) + 1 by push_cast; ring, hcoef_succ]
  simp only [ddR18, Znil, S03.Z_cons_succ, gqm1, gq2, gqm2, gq3, gqm3]
  push_cast
  simp only [pow_succ, pow_mul, neg_one_sq, one_pow]
  field_simp
  ring

lemma hE37 (m : ℕ) : hcoef [0, -1, 0, -1] (2 * ((m + 1 : ℕ) : ℤ) - 1) =
      (-1 : ℚ) * (1 / (2 * (m : ℚ) + 1) ^ 2) * S03.Z [2] (2 * m) := by
  rw [show 2 * ((m + 1 : ℕ) : ℤ) - 1 = ((2 * m : ℕ) : ℤ) + 1 by push_cast; ring, hcoef_succ]
  simp only [ddR18, Znil, S03.Z_cons_succ, gqm1, gq2, gqm2, gq3, gqm3]
  push_cast
  simp only [pow_succ, pow_mul, neg_one_sq, one_pow]
  field_simp
  ring

lemma hO37 (m : ℕ) : hcoef [0, -1, 0, -1] (2 * ((m + 1 : ℕ) : ℤ) + 1 - 1) =
      (1 : ℚ) * (1 / ((m : ℚ) + 1) ^ 1) +
      (1 / 4 : ℚ) * (1 / ((m : ℚ) + 1) ^ 2) +
      (-2 : ℚ) * (1 / (2 * (m : ℚ) + 1) ^ 1) +
      (1 : ℚ) * (1 / (2 * (m : ℚ) + 1) ^ 2) +
      (1 / 4 : ℚ) * (1 / ((m : ℚ) + 1) ^ 2) * S03.Z [2] (2 * m) := by
  rw [show 2 * ((m + 1 : ℕ) : ℤ) + 1 - 1 = ((2 * m + 1 : ℕ) : ℤ) + 1 by push_cast; ring, hcoef_succ]
  simp only [ddR18, Znil, S03.Z_cons_succ, gqm1, gq2, gqm2, gq3, gqm3]
  push_cast
  simp only [pow_succ, pow_mul, neg_one_sq, one_pow]
  field_simp
  ring

lemma hE38 (m : ℕ) : hcoef [0, -1, 1, -1] (2 * ((m + 1 : ℕ) : ℤ) - 0) =
      (1 / 2 : ℚ) * (1 / ((m : ℚ) + 1) ^ 1) * S03.Z [-1] (2 * m) +
      (1 / 4 : ℚ) * (1 / ((m : ℚ) + 1) ^ 2) * S03.Z [-1] (2 * m) +
      (-1 : ℚ) * (1 / (2 * (m : ℚ) + 1) ^ 1) * S03.Z [-1] (2 * m) +
      (1 / 4 : ℚ) * (1 / ((m : ℚ) + 1) ^ 2) * S03.Z [-1, -1] (2 * m) := by
  rw [show 2 * ((m + 1 : ℕ) : ℤ) - 0 = ((2 * m + 1 : ℕ) : ℤ) + 1 by push_cast; ring, hcoef_succ]
  simp only [ddR19, Znil, S03.Z_cons_succ, gqm1, gq2, gqm2, gq3, gqm3]
  push_cast
  simp only [pow_succ, pow_mul, neg_one_sq, one_pow]
  field_simp
  ring

lemma hO38 (m : ℕ) : hcoef [0, -1, 1, -1] (2 * ((m + 1 : ℕ) : ℤ) + 1 - 0) =
      (-1 / 2 : ℚ) * (1 / ((m : ℚ) + 1) ^ 1) +
      (1 / 4 : ℚ) * (1 / (2 * (m : ℚ) + 1) ^ 1) +
      (3 / 4 : ℚ) * (1 / (2 * (m : ℚ) + 3) ^ 1) +
      (1 / 2 : ℚ) * (1 / (2 * (m : ℚ) + 3) ^ 2) +
      (-1 / 2 : ℚ) * (1 / ((m : ℚ) + 1) ^ 1) * S03.Z [-1] (2 * m) +
      (1 / 4 : ℚ) * (1 / (2 * (m : ℚ) + 1) ^ 1) * S03.Z [-1] (2 * m) +
      (3 / 4 : ℚ) * (1 / (2 * (m : ℚ) + 3) ^ 1) * S03.Z [-1] (2 * m) +
      (1 / 2 : ℚ) * (1 / (2 * (m : ℚ) + 3) ^ 2) * S03.Z [-1] (2 * m) +
      (-1 : ℚ) * (1 / (2 * (m : ℚ) + 3) ^ 2) * S03.Z [-1, -1] (2 * m) := by
  rw [show 2 * ((m + 1 : ℕ) : ℤ) + 1 - 0 = ((2 * m + 1 + 1 : ℕ) : ℤ) + 1 by push_cast; ring, hcoef_succ]
  simp only [ddR19, Znil, S03.Z_cons_succ, gqm1, gq2, gqm2, gq3, gqm3]
  push_cast
  simp only [pow_succ, pow_mul, neg_one_sq, one_pow]
  field_simp
  ring

lemma hE39 (m : ℕ) : hcoef [0, -1, 1, -1] (2 * ((m + 1 : ℕ) : ℤ) - 1) =
      (-1 : ℚ) * (1 / (2 * (m : ℚ) + 1) ^ 2) * S03.Z [-1, -1] (2 * m) := by
  rw [show 2 * ((m + 1 : ℕ) : ℤ) - 1 = ((2 * m : ℕ) : ℤ) + 1 by push_cast; ring, hcoef_succ]
  simp only [ddR19, Znil, S03.Z_cons_succ, gqm1, gq2, gqm2, gq3, gqm3]
  push_cast
  simp only [pow_succ, pow_mul, neg_one_sq, one_pow]
  field_simp
  ring

lemma hO39 (m : ℕ) : hcoef [0, -1, 1, -1] (2 * ((m + 1 : ℕ) : ℤ) + 1 - 1) =
      (1 / 2 : ℚ) * (1 / ((m : ℚ) + 1) ^ 1) * S03.Z [-1] (2 * m) +
      (1 / 4 : ℚ) * (1 / ((m : ℚ) + 1) ^ 2) * S03.Z [-1] (2 * m) +
      (-1 : ℚ) * (1 / (2 * (m : ℚ) + 1) ^ 1) * S03.Z [-1] (2 * m) +
      (1 / 4 : ℚ) * (1 / ((m : ℚ) + 1) ^ 2) * S03.Z [-1, -1] (2 * m) := by
  rw [show 2 * ((m + 1 : ℕ) : ℤ) + 1 - 1 = ((2 * m + 1 : ℕ) : ℤ) + 1 by push_cast; ring, hcoef_succ]
  simp only [ddR19, Znil, S03.Z_cons_succ, gqm1, gq2, gqm2, gq3, gqm3]
  push_cast
  simp only [pow_succ, pow_mul, neg_one_sq, one_pow]
  field_simp
  ring

lemma hE40 (m : ℕ) : hcoef [1, 0, 0, 1] (2 * ((m + 1 : ℕ) : ℤ) - -1) =
      (1 / 2 : ℚ) * (1 / ((m : ℚ) + 1) ^ 1) +
      (-1 / 4 : ℚ) * (1 / ((m : ℚ) + 1) ^ 2) +
      (1 / 8 : ℚ) * (1 / ((m : ℚ) + 1) ^ 3) +
      (1 / 8 : ℚ) * (1 / (2 * (m : ℚ) + 1) ^ 1) +
      (-1 / 4 : ℚ) * (1 / (2 * (m : ℚ) + 1) ^ 2) +
      (1 / 2 : ℚ) * (1 / (2 * (m : ℚ) + 1) ^ 3) +
      (-9 / 8 : ℚ) * (1 / (2 * (m : ℚ) + 3) ^ 1) +
      (1 : ℚ) * (1 / (2 * (m : ℚ) + 3) ^ 1) * S03.Z [3] (2 * m) := by
  rw [show 2 * ((m + 1 : ℕ) : ℤ) - -1 = ((2 * m + 1 + 1 : ℕ) : ℤ) + 1 by push_cast; ring, hcoef_succ]
  simp only [ddR20, Znil, S03.Z_cons_succ, gqm1, gq2, gqm2, gq3, gqm3]
  push_cast
  simp only [pow_succ, pow_mul, neg_one_sq, one_pow]
  field_simp
  ring

lemma hO40 (m : ℕ) : hcoef [1, 0, 0, 1] (2 * ((m + 1 : ℕ) : ℤ) + 1 - -1) =
      (1 / 16 : ℚ) * (1 / ((m : ℚ) + 1) ^ 1) +
      (-1 / 16 : ℚ) * (1 / ((m : ℚ) + 1) ^ 2) +
      (1 / 16 : ℚ) * (1 / ((m : ℚ) + 1) ^ 3) +
      (1 / 27 : ℚ) * (1 / (2 * (m : ℚ) + 1) ^ 1) +
      (-1 / 9 : ℚ) * (1 / (2 * (m : ℚ) + 1) ^ 2) +
      (1 / 3 : ℚ) * (1 / (2 * (m : ℚ) + 1) ^ 3) +
      (-251 / 432 : ℚ) * (1 / ((m : ℚ) + 2) ^ 1) +
      (1 : ℚ) * (1 / (2 * (m : ℚ) + 3) ^ 1) +
      (-1 : ℚ) * (1 / (2 * (m : ℚ) + 3) ^ 2) +
      (1 : ℚ) * (1 / (2 * (m : ℚ) + 3) ^ 3) +
      (1 / 2 : ℚ) * (1 / ((m : ℚ) + 2) ^ 1) * S03.Z [3] (2 * m) := by
  rw [show 2 * ((m + 1 : ℕ) : ℤ) + 1 - -1 = ((2 * m + 1 + 1 + 1 : ℕ) : ℤ) + 1 by push_cast; ring, hcoef_succ]
  simp only [ddR20, Znil, S03.Z_cons_succ, gqm1, gq2, gqm2, gq3, gqm3]
  push_cast
  simp only [pow_succ, pow_mul, neg_one_sq, one_pow]
  field_simp
  ring

lemma hE41 (m : ℕ) : hcoef [1, 0, 0, 1] (2 * ((m + 1 : ℕ) : ℤ) - 1) =
      (1 : ℚ) * (1 / (2 * (m : ℚ) + 1) ^ 1) * S03.Z [3] (2 * m) := by
  rw [show 2 * ((m + 1 : ℕ) : ℤ) - 1 = ((2 * m : ℕ) : ℤ) + 1 by push_cast; ring, hcoef_succ]
  simp only [ddR20, Znil, S03.Z_cons_succ, gqm1, gq2, gqm2, gq3, gqm3]
  push_cast
  simp only [pow_succ, pow_mul, neg_one_sq, one_pow]
  field_simp
  ring

lemma hO41 (m : ℕ) : hcoef [1, 0, 0, 1] (2 * ((m + 1 : ℕ) : ℤ) + 1 - 1) =
      (-1 / 2 : ℚ) * (1 / ((m : ℚ) + 1) ^ 1) +
      (1 : ℚ) * (1 / (2 * (m : ℚ) + 1) ^ 1) +
      (-1 : ℚ) * (1 / (2 * (m : ℚ) + 1) ^ 2) +
      (1 : ℚ) * (1 / (2 * (m : ℚ) + 1) ^ 3) +
      (1 / 2 : ℚ) * (1 / ((m : ℚ) + 1) ^ 1) * S03.Z [3] (2 * m) := by
  rw [show 2 * ((m + 1 : ℕ) : ℤ) + 1 - 1 = ((2 * m + 1 : ℕ) : ℤ) + 1 by push_cast; ring, hcoef_succ]
  simp only [ddR20, Znil, S03.Z_cons_succ, gqm1, gq2, gqm2, gq3, gqm3]
  push_cast
  simp only [pow_succ, pow_mul, neg_one_sq, one_pow]
  field_simp
  ring

lemma hE42 (m : ℕ) : hcoef [1, 0, 0, -1] (2 * ((m + 1 : ℕ) : ℤ) - -1) =
      (-1 / 2 : ℚ) * (1 / ((m : ℚ) + 1) ^ 1) +
      (1 / 4 : ℚ) * (1 / ((m : ℚ) + 1) ^ 2) +
      (-1 / 8 : ℚ) * (1 / ((m : ℚ) + 1) ^ 3) +
      (1 / 8 : ℚ) * (1 / (2 * (m : ℚ) + 1) ^ 1) +
      (-1 / 4 : ℚ) * (1 / (2 * (m : ℚ) + 1) ^ 2) +
      (1 / 2 : ℚ) * (1 / (2 * (m : ℚ) + 1) ^ 3) +
      (7 / 8 : ℚ) * (1 / (2 * (m : ℚ) + 3) ^ 1) +
      (-1 : ℚ) * (1 / (2 * (m : ℚ) + 3) ^ 1) * S03.Z [-3] (2 * m) := by
  rw [show 2 * ((m + 1 : ℕ) : ℤ) - -1 = ((2 * m + 1 + 1 : ℕ) : ℤ) + 1 by push_cast; ring, hcoef_succ]
  simp only [ddR21, Znil, S03.Z_cons_succ, gqm1, gq2, gqm2, gq3, gqm3]
  push_cast
  simp only [pow_succ, pow_mul, neg_one_sq, one_pow]
  field_simp
  ring

lemma hO42 (m : ℕ) : hcoef [1, 0, 0, -1] (2 * ((m + 1 : ℕ) : ℤ) + 1 - -1) =
      (-1 / 16 : ℚ) * (1 / ((m : ℚ) + 1) ^ 1) +
      (1 / 16 : ℚ) * (1 / ((m : ℚ) + 1) ^ 2) +
      (-1 / 16 : ℚ) * (1 / ((m : ℚ) + 1) ^ 3) +
      (1 / 27 : ℚ) * (1 / (2 * (m : ℚ) + 1) ^ 1) +
      (-1 / 9 : ℚ) * (1 / (2 * (m : ℚ) + 1) ^ 2) +
      (1 / 3 : ℚ) * (1 / (2 * (m : ℚ) + 1) ^ 3) +
      (-197 / 432 : ℚ) * (1 / ((m : ℚ) + 2) ^ 1) +
      (1 : ℚ) * (1 / (2 * (m : ℚ) + 3) ^ 1) +
      (-1 : ℚ) * (1 / (2 * (m : ℚ) + 3) ^ 2) +
      (1 : ℚ) * (1 / (2 * (m : ℚ) + 3) ^ 3) +
      (-1 / 2 : ℚ) * (1 / ((m : ℚ) + 2) ^ 1) * S03.Z [-3] (2 * m) := by
  rw [show 2 * ((m + 1 : ℕ) : ℤ) + 1 - -1 = ((2 * m + 1 + 1 + 1 : ℕ) : ℤ) + 1 by push_cast; ring, hcoef_succ]
  simp only [ddR21, Znil, S03.Z_cons_succ, gqm1, gq2, gqm2, gq3, gqm3]
  push_cast
  simp only [pow_succ, pow_mul, neg_one_sq, one_pow]
  field_simp
  ring

lemma hE43 (m : ℕ) : hcoef [1, 0, 0, -1] (2 * ((m + 1 : ℕ) : ℤ) - 0) =
      (-1 / 2 : ℚ) * (1 / ((m : ℚ) + 1) ^ 1) +
      (1 : ℚ) * (1 / (2 * (m : ℚ) + 1) ^ 1) +
      (-1 : ℚ) * (1 / (2 * (m : ℚ) + 1) ^ 2) +
      (1 : ℚ) * (1 / (2 * (m : ℚ) + 1) ^ 3) +
      (-1 / 2 : ℚ) * (1 / ((m : ℚ) + 1) ^ 1) * S03.Z [-3] (2 * m) := by
  rw [show 2 * ((m + 1 : ℕ) : ℤ) - 0 = ((2 * m + 1 : ℕ) : ℤ) + 1 by push_cast; ring, hcoef_succ]
  simp only [ddR21, Znil, S03.Z_cons_succ, gqm1, gq2, gqm2, gq3, gqm3]
  push_cast
  simp only [pow_succ, pow_mul, neg_one_sq, one_pow]
  field_simp
  ring

lemma hO43 (m : ℕ) : hcoef [1, 0, 0, -1] (2 * ((m + 1 : ℕ) : ℤ) + 1 - 0) =
      (-1 / 2 : ℚ) * (1 / ((m : ℚ) + 1) ^ 1) +
      (1 / 4 : ℚ) * (1 / ((m : ℚ) + 1) ^ 2) +
      (-1 / 8 : ℚ) * (1 / ((m : ℚ) + 1) ^ 3) +
      (1 / 8 : ℚ) * (1 / (2 * (m : ℚ) + 1) ^ 1) +
      (-1 / 4 : ℚ) * (1 / (2 * (m : ℚ) + 1) ^ 2) +
      (1 / 2 : ℚ) * (1 / (2 * (m : ℚ) + 1) ^ 3) +
      (7 / 8 : ℚ) * (1 / (2 * (m : ℚ) + 3) ^ 1) +
      (-1 : ℚ) * (1 / (2 * (m : ℚ) + 3) ^ 1) * S03.Z [-3] (2 * m) := by
  rw [show 2 * ((m + 1 : ℕ) : ℤ) + 1 - 0 = ((2 * m + 1 + 1 : ℕ) : ℤ) + 1 by push_cast; ring, hcoef_succ]
  simp only [ddR21, Znil, S03.Z_cons_succ, gqm1, gq2, gqm2, gq3, gqm3]
  push_cast
  simp only [pow_succ, pow_mul, neg_one_sq, one_pow]
  field_simp
  ring

lemma hE44 (m : ℕ) : hcoef [1, 0, 0, -1] (2 * ((m + 1 : ℕ) : ℤ) - 1) =
      (-1 : ℚ) * (1 / (2 * (m : ℚ) + 1) ^ 1) * S03.Z [-3] (2 * m) := by
  rw [show 2 * ((m + 1 : ℕ) : ℤ) - 1 = ((2 * m : ℕ) : ℤ) + 1 by push_cast; ring, hcoef_succ]
  simp only [ddR21, Znil, S03.Z_cons_succ, gqm1, gq2, gqm2, gq3, gqm3]
  push_cast
  simp only [pow_succ, pow_mul, neg_one_sq, one_pow]
  field_simp
  ring

lemma hO44 (m : ℕ) : hcoef [1, 0, 0, -1] (2 * ((m + 1 : ℕ) : ℤ) + 1 - 1) =
      (-1 / 2 : ℚ) * (1 / ((m : ℚ) + 1) ^ 1) +
      (1 : ℚ) * (1 / (2 * (m : ℚ) + 1) ^ 1) +
      (-1 : ℚ) * (1 / (2 * (m : ℚ) + 1) ^ 2) +
      (1 : ℚ) * (1 / (2 * (m : ℚ) + 1) ^ 3) +
      (-1 / 2 : ℚ) * (1 / ((m : ℚ) + 1) ^ 1) * S03.Z [-3] (2 * m) := by
  rw [show 2 * ((m + 1 : ℕ) : ℤ) + 1 - 1 = ((2 * m + 1 : ℕ) : ℤ) + 1 by push_cast; ring, hcoef_succ]
  simp only [ddR21, Znil, S03.Z_cons_succ, gqm1, gq2, gqm2, gq3, gqm3]
  push_cast
  simp only [pow_succ, pow_mul, neg_one_sq, one_pow]
  field_simp
  ring

lemma hE45 (m : ℕ) : hcoef [1, 0, 1, -1] (2 * ((m + 1 : ℕ) : ℤ) - -1) =
      (-1 / 4 : ℚ) * (1 / ((m : ℚ) + 1) ^ 2) +
      (1 / 2 : ℚ) * (1 / (2 * (m : ℚ) + 1) ^ 1) +
      (-1 / 2 : ℚ) * (1 / (2 * (m : ℚ) + 3) ^ 1) +
      (-1 : ℚ) * (1 / (2 * (m : ℚ) + 3) ^ 1) * S03.Z [2, -1] (2 * m) +
      (1 / 2 : ℚ) * (1 / ((m : ℚ) + 1) ^ 1) * S03.Z [-1] (2 * m) +
      (-1 / 4 : ℚ) * (1 / ((m : ℚ) + 1) ^ 2) * S03.Z [-1] (2 * m) +
      (1 / 4 : ℚ) * (1 / (2 * (m : ℚ) + 1) ^ 1) * S03.Z [-1] (2 * m) +
      (-1 / 2 : ℚ) * (1 / (2 * (m : ℚ) + 1) ^ 2) * S03.Z [-1] (2 * m) +
      (-5 / 4 : ℚ) * (1 / (2 * (m : ℚ) + 3) ^ 1) * S03.Z [-1] (2 * m) := by
  rw [show 2 * ((m + 1 : ℕ) : ℤ) - -1 = ((2 * m + 1 + 1 : ℕ) : ℤ) + 1 by push_cast; ring, hcoef_succ]
  simp only [ddR22, Znil, S03.Z_cons_succ, gqm1, gq2, gqm2, gq3, gqm3]
  push_cast
  simp only [pow_succ, pow_mul, neg_one_sq, one_pow]
  field_simp
  ring

lemma hO45 (m : ℕ) : hcoef [1, 0, 1, -1] (2 * ((m + 1 : ℕ) : ℤ) + 1 - -1) =
      (-3 / 8 : ℚ) * (1 / ((m : ℚ) + 1) ^ 1) +
      (-1 / 8 : ℚ) * (1 / ((m : ℚ) + 1) ^ 2) +
      (5 / 12 : ℚ) * (1 / (2 * (m : ℚ) + 1) ^ 1) +
      (1 / 24 : ℚ) * (1 / ((m : ℚ) + 2) ^ 1) +
      (1 / 4 : ℚ) * (1 / (2 * (m : ℚ) + 3) ^ 1) +
      (1 / 2 : ℚ) * (1 / (2 * (m : ℚ) + 3) ^ 2) +
      (-1 / 2 : ℚ) * (1 / ((m : ℚ) + 2) ^ 1) * S03.Z [2, -1] (2 * m) +
      (1 / 8 : ℚ) * (1 / ((m : ℚ) + 1) ^ 1) * S03.Z [-1] (2 * m) +
      (-1 / 8 : ℚ) * (1 / ((m : ℚ) + 1) ^ 2) * S03.Z [-1] (2 * m) +
      (1 / 9 : ℚ) * (1 / (2 * (m : ℚ) + 1) ^ 1) * S03.Z [-1] (2 * m) +
      (-1 / 3 : ℚ) * (1 / (2 * (m : ℚ) + 1) ^ 2) * S03.Z [-1] (2 * m) +
      (-49 / 72 : ℚ) * (1 / ((m : ℚ) + 2) ^ 1) * S03.Z [-1] (2 * m) +
      (1 : ℚ) * (1 / (2 * (m : ℚ) + 3) ^ 1) * S03.Z [-1] (2 * m) +
      (-1 : ℚ) * (1 / (2 * (m : ℚ) + 3) ^ 2) * S03.Z [-1] (2 * m) := by
  rw [show 2 * ((m + 1 : ℕ) : ℤ) + 1 - -1 = ((2 * m + 1 + 1 + 1 : ℕ) : ℤ) + 1 by push_cast; ring, hcoef_succ]
  simp only [ddR22, Znil, S03.Z_cons_succ, gqm1, gq2, gqm2, gq3, gqm3]
  push_cast
  simp only [pow_succ, pow_mul, neg_one_sq, one_pow]
  field_simp
  ring

lemma hE46 (m : ℕ) : hcoef [1, 0, 1, -1] (2 * ((m + 1 : ℕ) : ℤ) - 1) =
      (-1 : ℚ) * (1 / (2 * (m : ℚ) + 1) ^ 1) * S03.Z [2, -1] (2 * m) := by
  rw [show 2 * ((m + 1 : ℕ) : ℤ) - 1 = ((2 * m : ℕ) : ℤ) + 1 by push_cast; ring, hcoef_succ]
  simp only [ddR22, Znil, S03.Z_cons_succ, gqm1, gq2, gqm2, gq3, gqm3]
  push_cast
  simp only [pow_succ, pow_mul, neg_one_sq, one_pow]
  field_simp
  ring

lemma hO46 (m : ℕ) : hcoef [1, 0, 1, -1] (2 * ((m + 1 : ℕ) : ℤ) + 1 - 1) =
      (-1 / 2 : ℚ) * (1 / ((m : ℚ) + 1) ^ 1) * S03.Z [2, -1] (2 * m) +
      (-1 / 2 : ℚ) * (1 / ((m : ℚ) + 1) ^ 1) * S03.Z [-1] (2 * m) +
      (1 : ℚ) * (1 / (2 * (m : ℚ) + 1) ^ 1) * S03.Z [-1] (2 * m) +
      (-1 : ℚ) * (1 / (2 * (m : ℚ) + 1) ^ 2) * S03.Z [-1] (2 * m) := by
  rw [show 2 * ((m + 1 : ℕ) : ℤ) + 1 - 1 = ((2 * m + 1 : ℕ) : ℤ) + 1 by push_cast; ring, hcoef_succ]
  simp only [ddR22, Znil, S03.Z_cons_succ, gqm1, gq2, gqm2, gq3, gqm3]
  push_cast
  simp only [pow_succ, pow_mul, neg_one_sq, one_pow]
  field_simp
  ring

lemma hE47 (m : ℕ) : hcoef [1, 0, -1, 1] (2 * ((m + 1 : ℕ) : ℤ) - -1) =
      (-1 / 4 : ℚ) * (1 / ((m : ℚ) + 1) ^ 2) +
      (1 / 2 : ℚ) * (1 / (2 * (m : ℚ) + 1) ^ 1) +
      (-1 / 2 : ℚ) * (1 / (2 * (m : ℚ) + 3) ^ 1) +
      (1 / 2 : ℚ) * (1 / ((m : ℚ) + 1) ^ 1) * S03.Z [-1] (2 * m) +
      (-1 / 4 : ℚ) * (1 / ((m : ℚ) + 1) ^ 2) * S03.Z [-1] (2 * m) +
      (-1 / 4 : ℚ) * (1 / (2 * (m : ℚ) + 1) ^ 1) * S03.Z [-1] (2 * m) +
      (1 / 2 : ℚ) * (1 / (2 * (m : ℚ) + 1) ^ 2) * S03.Z [-1] (2 * m) +
      (-3 / 4 : ℚ) * (1 / (2 * (m : ℚ) + 3) ^ 1) * S03.Z [-1] (2 * m) +
      (-1 : ℚ) * (1 / (2 * (m : ℚ) + 3) ^ 1) * S03.Z [-2, -1] (2 * m) := by
  rw [show 2 * ((m + 1 : ℕ) : ℤ) - -1 = ((2 * m + 1 + 1 : ℕ) : ℤ) + 1 by push_cast; ring, hcoef_succ]
  simp only [ddR23, Znil, S03.Z_cons_succ, gqm1, gq2, gqm2, gq3, gqm3]
  push_cast
  simp only [pow_succ, pow_mul, neg_one_sq, one_pow]
  field_simp
  ring

lemma hO47 (m : ℕ) : hcoef [1, 0, -1, 1] (2 * ((m + 1 : ℕ) : ℤ) + 1 - -1) =
      (1 / 8 : ℚ) * (1 / ((m : ℚ) + 1) ^ 1) +
      (-1 / 8 : ℚ) * (1 / ((m : ℚ) + 1) ^ 2) +
      (1 / 4 : ℚ) * (1 / (2 * (m : ℚ) + 1) ^ 1) +
      (-1 / 8 : ℚ) * (1 / ((m : ℚ) + 2) ^ 1) +
      (-1 / 4 : ℚ) * (1 / (2 * (m : ℚ) + 3) ^ 1) +
      (-1 / 2 : ℚ) * (1 / (2 * (m : ℚ) + 3) ^ 2) +
      (1 / 8 : ℚ) * (1 / ((m : ℚ) + 1) ^ 1) * S03.Z [-1] (2 * m) +
      (-1 / 8 : ℚ) * (1 / ((m : ℚ) + 1) ^ 2) * S03.Z [-1] (2 * m) +
      (-1 / 9 : ℚ) * (1 / (2 * (m : ℚ) + 1) ^ 1) * S03.Z [-1] (2 * m) +
      (1 / 3 : ℚ) * (1 / (2 * (m : ℚ) + 1) ^ 2) * S03.Z [-1] (2 * m) +
      (31 / 72 : ℚ) * (1 / ((m : ℚ) + 2) ^ 1) * S03.Z [-1] (2 * m) +
      (-1 : ℚ) * (1 / (2 * (m : ℚ) + 3) ^ 1) * S03.Z [-1] (2 * m) +
      (1 : ℚ) * (1 / (2 * (m : ℚ) + 3) ^ 2) * S03.Z [-1] (2 * m) +
      (-1 / 2 : ℚ) * (1 / ((m : ℚ) + 2) ^ 1) * S03.Z [-2, -1] (2 * m) := by
  rw [show 2 * ((m + 1 : ℕ) : ℤ) + 1 - -1 = ((2 * m + 1 + 1 + 1 : ℕ) : ℤ) + 1 by push_cast; ring, hcoef_succ]
  simp only [ddR23, Znil, S03.Z_cons_succ, gqm1, gq2, gqm2, gq3, gqm3]
  push_cast
  simp only [pow_succ, pow_mul, neg_one_sq, one_pow]
  field_simp
  ring

lemma hE48 (m : ℕ) : hcoef [1, 0, -1, 1] (2 * ((m + 1 : ℕ) : ℤ) - 0) =
      (1 / 2 : ℚ) * (1 / ((m : ℚ) + 1) ^ 1) * S03.Z [-1] (2 * m) +
      (-1 : ℚ) * (1 / (2 * (m : ℚ) + 1) ^ 1) * S03.Z [-1] (2 * m) +
      (1 : ℚ) * (1 / (2 * (m : ℚ) + 1) ^ 2) * S03.Z [-1] (2 * m) +
      (-1 / 2 : ℚ) * (1 / ((m : ℚ) + 1) ^ 1) * S03.Z [-2, -1] (2 * m) := by
  rw [show 2 * ((m + 1 : ℕ) : ℤ) - 0 = ((2 * m + 1 : ℕ) : ℤ) + 1 by push_cast; ring, hcoef_succ]
  simp only [ddR23, Znil, S03.Z_cons_succ, gqm1, gq2, gqm2, gq3, gqm3]
  push_cast
  simp only [pow_succ, pow_mul, neg_one_sq, one_pow]
  field_simp
  ring

lemma hO48 (m : ℕ) : hcoef [1, 0, -1, 1] (2 * ((m + 1 : ℕ) : ℤ) + 1 - 0) =
      (-1 / 4 : ℚ) * (1 / ((m : ℚ) + 1) ^ 2) +
      (1 / 2 : ℚ) * (1 / (2 * (m : ℚ) + 1) ^ 1) +
      (-1 / 2 : ℚ) * (1 / (2 * (m : ℚ) + 3) ^ 1) +
      (1 / 2 : ℚ) * (1 / ((m : ℚ) + 1) ^ 1) * S03.Z [-1] (2 * m) +
      (-1 / 4 : ℚ) * (1 / ((m : ℚ) + 1) ^ 2) * S03.Z [-1] (2 * m) +
      (-1 / 4 : ℚ) * (1 / (2 * (m : ℚ) + 1) ^ 1) * S03.Z [-1] (2 * m) +
      (1 / 2 : ℚ) * (1 / (2 * (m : ℚ) + 1) ^ 2) * S03.Z [-1] (2 * m) +
      (-3 / 4 : ℚ) * (1 / (2 * (m : ℚ) + 3) ^ 1) * S03.Z [-1] (2 * m) +
      (-1 : ℚ) * (1 / (2 * (m : ℚ) + 3) ^ 1) * S03.Z [-2, -1] (2 * m) := by
  rw [show 2 * ((m + 1 : ℕ) : ℤ) + 1 - 0 = ((2 * m + 1 + 1 : ℕ) : ℤ) + 1 by push_cast; ring, hcoef_succ]
  simp only [ddR23, Znil, S03.Z_cons_succ, gqm1, gq2, gqm2, gq3, gqm3]
  push_cast
  simp only [pow_succ, pow_mul, neg_one_sq, one_pow]
  field_simp
  ring

lemma hE49 (m : ℕ) : hcoef [1, 0, -1, 1] (2 * ((m + 1 : ℕ) : ℤ) - 1) =
      (-1 : ℚ) * (1 / (2 * (m : ℚ) + 1) ^ 1) * S03.Z [-2, -1] (2 * m) := by
  rw [show 2 * ((m + 1 : ℕ) : ℤ) - 1 = ((2 * m : ℕ) : ℤ) + 1 by push_cast; ring, hcoef_succ]
  simp only [ddR23, Znil, S03.Z_cons_succ, gqm1, gq2, gqm2, gq3, gqm3]
  push_cast
  simp only [pow_succ, pow_mul, neg_one_sq, one_pow]
  field_simp
  ring

lemma hO49 (m : ℕ) : hcoef [1, 0, -1, 1] (2 * ((m + 1 : ℕ) : ℤ) + 1 - 1) =
      (1 / 2 : ℚ) * (1 / ((m : ℚ) + 1) ^ 1) * S03.Z [-1] (2 * m) +
      (-1 : ℚ) * (1 / (2 * (m : ℚ) + 1) ^ 1) * S03.Z [-1] (2 * m) +
      (1 : ℚ) * (1 / (2 * (m : ℚ) + 1) ^ 2) * S03.Z [-1] (2 * m) +
      (-1 / 2 : ℚ) * (1 / ((m : ℚ) + 1) ^ 1) * S03.Z [-2, -1] (2 * m) := by
  rw [show 2 * ((m + 1 : ℕ) : ℤ) + 1 - 1 = ((2 * m + 1 : ℕ) : ℤ) + 1 by push_cast; ring, hcoef_succ]
  simp only [ddR23, Znil, S03.Z_cons_succ, gqm1, gq2, gqm2, gq3, gqm3]
  push_cast
  simp only [pow_succ, pow_mul, neg_one_sq, one_pow]
  field_simp
  ring

lemma hE50 (m : ℕ) : hcoef [1, -1, 0, 1] (2 * ((m + 1 : ℕ) : ℤ) - -1) =
      (1 / 2 : ℚ) * (1 / ((m : ℚ) + 1) ^ 1) +
      (-3 / 4 : ℚ) * (1 / (2 * (m : ℚ) + 1) ^ 1) +
      (1 / 2 : ℚ) * (1 / (2 * (m : ℚ) + 1) ^ 2) +
      (-1 / 4 : ℚ) * (1 / (2 * (m : ℚ) + 3) ^ 1) +
      (-1 : ℚ) * (1 / (2 * (m : ℚ) + 3) ^ 1) * S03.Z [-1, -2] (2 * m) +
      (-1 / 2 : ℚ) * (1 / ((m : ℚ) + 1) ^ 1) * S03.Z [-2] (2 * m) +
      (1 / 2 : ℚ) * (1 / (2 * (m : ℚ) + 1) ^ 1) * S03.Z [-2] (2 * m) +
      (1 / 2 : ℚ) * (1 / (2 * (m : ℚ) + 3) ^ 1) * S03.Z [-2] (2 * m) := by
  rw [show 2 * ((m + 1 : ℕ) : ℤ) - -1 = ((2 * m + 1 + 1 : ℕ) : ℤ) + 1 by push_cast; ring, hcoef_succ]
  simp only [ddR24, Znil, S03.Z_cons_succ, gqm1, gq2, gqm2, gq3, gqm3]
  push_cast
  simp only [pow_succ, pow_mul, neg_one_sq, one_pow]
  field_simp
  ring

lemma hO50 (m : ℕ) : hcoef [1, -1, 0, 1] (2 * ((m + 1 : ℕ) : ℤ) + 1 - -1) =
      (-1 / 8 : ℚ) * (1 / ((m : ℚ) + 1) ^ 1) +
      (1 / 8 : ℚ) * (1 / ((m : ℚ) + 1) ^ 2) +
      (-11 / 36 : ℚ) * (1 / (2 * (m : ℚ) + 1) ^ 1) +
      (1 / 6 : ℚ) * (1 / (2 * (m : ℚ) + 1) ^ 2) +
      (-7 / 72 : ℚ) * (1 / ((m : ℚ) + 2) ^ 1) +
      (3 / 4 : ℚ) * (1 / (2 * (m : ℚ) + 3) ^ 1) +
      (-1 / 2 : ℚ) * (1 / ((m : ℚ) + 2) ^ 1) * S03.Z [-1, -2] (2 * m) +
      (-1 / 4 : ℚ) * (1 / ((m : ℚ) + 1) ^ 1) * S03.Z [-2] (2 * m) +
      (1 / 3 : ℚ) * (1 / (2 * (m : ℚ) + 1) ^ 1) * S03.Z [-2] (2 * m) +
      (-5 / 12 : ℚ) * (1 / ((m : ℚ) + 2) ^ 1) * S03.Z [-2] (2 * m) +
      (1 : ℚ) * (1 / (2 * (m : ℚ) + 3) ^ 1) * S03.Z [-2] (2 * m) := by
  rw [show 2 * ((m + 1 : ℕ) : ℤ) + 1 - -1 = ((2 * m + 1 + 1 + 1 : ℕ) : ℤ) + 1 by push_cast; ring, hcoef_succ]
  simp only [ddR24, Znil, S03.Z_cons_succ, gqm1, gq2, gqm2, gq3, gqm3]
  push_cast
  simp only [pow_succ, pow_mul, neg_one_sq, one_pow]
  field_simp
  ring

lemma hE51 (m : ℕ) : hcoef [1, -1, 0, 1] (2 * ((m + 1 : ℕ) : ℤ) - 0) =
      (-1 / 2 : ℚ) * (1 / ((m : ℚ) + 1) ^ 1) * S03.Z [-1, -2] (2 * m) +
      (-1 / 2 : ℚ) * (1 / ((m : ℚ) + 1) ^ 1) * S03.Z [-2] (2 * m) +
      (1 : ℚ) * (1 / (2 * (m : ℚ) + 1) ^ 1) * S03.Z [-2] (2 * m) := by
  rw [show 2 * ((m + 1 : ℕ) : ℤ) - 0 = ((2 * m + 1 : ℕ) : ℤ) + 1 by push_cast; ring, hcoef_succ]
  simp only [ddR24, Znil, S03.Z_cons_succ, gqm1, gq2, gqm2, gq3, gqm3]
  push_cast
  simp only [pow_succ, pow_mul, neg_one_sq, one_pow]
  field_simp
  ring

lemma hO51 (m : ℕ) : hcoef [1, -1, 0, 1] (2 * ((m + 1 : ℕ) : ℤ) + 1 - 0) =
      (1 / 2 : ℚ) * (1 / ((m : ℚ) + 1) ^ 1) +
      (-3 / 4 : ℚ) * (1 / (2 * (m : ℚ) + 1) ^ 1) +
      (1 / 2 : ℚ) * (1 / (2 * (m : ℚ) + 1) ^ 2) +
      (-1 / 4 : ℚ) * (1 / (2 * (m : ℚ) + 3) ^ 1) +
      (-1 : ℚ) * (1 / (2 * (m : ℚ) + 3) ^ 1) * S03.Z [-1, -2] (2 * m) +
      (-1 / 2 : ℚ) * (1 / ((m : ℚ) + 1) ^ 1) * S03.Z [-2] (2 * m) +
      (1 / 2 : ℚ) * (1 / (2 * (m : ℚ) + 1) ^ 1) * S03.Z [-2] (2 * m) +
      (1 / 2 : ℚ) * (1 / (2 * (m : ℚ) + 3) ^ 1) * S03.Z [-2] (2 * m) := by
  rw [show 2 * ((m + 1 : ℕ) : ℤ) + 1 - 0 = ((2 * m + 1 + 1 : ℕ) : ℤ) + 1 by push_cast; ring, hcoef_succ]
  simp only [ddR24, Znil, S03.Z_cons_succ, gqm1, gq2, gqm2, gq3, gqm3]
  push_cast
  simp only [pow_succ, pow_mul, neg_one_sq, one_pow]
  field_simp
  ring

lemma hE52 (m : ℕ) : hcoef [1, -1, 0, 1] (2 * ((m + 1 : ℕ) : ℤ) - 1) =
      (-1 : ℚ) * (1 / (2 * (m : ℚ) + 1) ^ 1) * S03.Z [-1, -2] (2 * m) := by
  rw [show 2 * ((m + 1 : ℕ) : ℤ) - 1 = ((2 * m : ℕ) : ℤ) + 1 by push_cast; ring, hcoef_succ]
  simp only [ddR24, Znil, S03.Z_cons_succ, gqm1, gq2, gqm2, gq3, gqm3]
  push_cast
  simp only [pow_succ, pow_mul, neg_one_sq, one_pow]
  field_simp
  ring

lemma hO52 (m : ℕ) : hcoef [1, -1, 0, 1] (2 * ((m + 1 : ℕ) : ℤ) + 1 - 1) =
      (-1 / 2 : ℚ) * (1 / ((m : ℚ) + 1) ^ 1) * S03.Z [-1, -2] (2 * m) +
      (-1 / 2 : ℚ) * (1 / ((m : ℚ) + 1) ^ 1) * S03.Z [-2] (2 * m) +
      (1 : ℚ) * (1 / (2 * (m : ℚ) + 1) ^ 1) * S03.Z [-2] (2 * m) := by
  rw [show 2 * ((m + 1 : ℕ) : ℤ) + 1 - 1 = ((2 * m + 1 : ℕ) : ℤ) + 1 by push_cast; ring, hcoef_succ]
  simp only [ddR24, Znil, S03.Z_cons_succ, gqm1, gq2, gqm2, gq3, gqm3]
  push_cast
  simp only [pow_succ, pow_mul, neg_one_sq, one_pow]
  field_simp
  ring

lemma hE53 (m : ℕ) : hcoef [1, -1, 0, -1] (2 * ((m + 1 : ℕ) : ℤ) - -1) =
      (1 / 2 : ℚ) * (1 / ((m : ℚ) + 1) ^ 1) +
      (-3 / 4 : ℚ) * (1 / (2 * (m : ℚ) + 1) ^ 1) +
      (1 / 2 : ℚ) * (1 / (2 * (m : ℚ) + 1) ^ 2) +
      (-1 / 4 : ℚ) * (1 / (2 * (m : ℚ) + 3) ^ 1) +
      (1 / 2 : ℚ) * (1 / ((m : ℚ) + 1) ^ 1) * S03.Z [2] (2 * m) +
      (-1 / 2 : ℚ) * (1 / (2 * (m : ℚ) + 1) ^ 1) * S03.Z [2] (2 * m) +
      (-1 / 2 : ℚ) * (1 / (2 * (m : ℚ) + 3) ^ 1) * S03.Z [2] (2 * m) +
      (1 : ℚ) * (1 / (2 * (m : ℚ) + 3) ^ 1) * S03.Z [-1, 2] (2 * m) := by
  rw [show 2 * ((m + 1 : ℕ) : ℤ) - -1 = ((2 * m + 1 + 1 : ℕ) : ℤ) + 1 by push_cast; ring, hcoef_succ]
  simp only [ddR25, Znil, S03.Z_cons_succ, gqm1, gq2, gqm2, gq3, gqm3]
  push_cast
  simp only [pow_succ, pow_mul, neg_one_sq, one_pow]
  field_simp
  ring

lemma hO53 (m : ℕ) : hcoef [1, -1, 0, -1] (2 * ((m + 1 : ℕ) : ℤ) + 1 - -1) =
      (5 / 8 : ℚ) * (1 / ((m : ℚ) + 1) ^ 1) +
      (-1 / 8 : ℚ) * (1 / ((m : ℚ) + 1) ^ 2) +
      (-11 / 36 : ℚ) * (1 / (2 * (m : ℚ) + 1) ^ 1) +
      (1 / 6 : ℚ) * (1 / (2 * (m : ℚ) + 1) ^ 2) +
      (11 / 72 : ℚ) * (1 / ((m : ℚ) + 2) ^ 1) +
      (-5 / 4 : ℚ) * (1 / (2 * (m : ℚ) + 3) ^ 1) +
      (1 / 4 : ℚ) * (1 / ((m : ℚ) + 1) ^ 1) * S03.Z [2] (2 * m) +
      (-1 / 3 : ℚ) * (1 / (2 * (m : ℚ) + 1) ^ 1) * S03.Z [2] (2 * m) +
      (5 / 12 : ℚ) * (1 / ((m : ℚ) + 2) ^ 1) * S03.Z [2] (2 * m) +
      (-1 : ℚ) * (1 / (2 * (m : ℚ) + 3) ^ 1) * S03.Z [2] (2 * m) +
      (1 / 2 : ℚ) * (1 / ((m : ℚ) + 2) ^ 1) * S03.Z [-1, 2] (2 * m) := by
  rw [show 2 * ((m + 1 : ℕ) : ℤ) + 1 - -1 = ((2 * m + 1 + 1 + 1 : ℕ) : ℤ) + 1 by push_cast; ring, hcoef_succ]
  simp only [ddR25, Znil, S03.Z_cons_succ, gqm1, gq2, gqm2, gq3, gqm3]
  push_cast
  simp only [pow_succ, pow_mul, neg_one_sq, one_pow]
  field_simp
  ring

lemma hE54 (m : ℕ) : hcoef [1, -1, 0, -1] (2 * ((m + 1 : ℕ) : ℤ) - 1) =
      (1 : ℚ) * (1 / (2 * (m : ℚ) + 1) ^ 1) * S03.Z [-1, 2] (2 * m) := by
  rw [show 2 * ((m + 1 : ℕ) : ℤ) - 1 = ((2 * m : ℕ) : ℤ) + 1 by push_cast; ring, hcoef_succ]
  simp only [ddR25, Znil, S03.Z_cons_succ, gqm1, gq2, gqm2, gq3, gqm3]
  push_cast
  simp only [pow_succ, pow_mul, neg_one_sq, one_pow]
  field_simp
  ring

lemma hO54 (m : ℕ) : hcoef [1, -1, 0, -1] (2 * ((m + 1 : ℕ) : ℤ) + 1 - 1) =
      (1 / 2 : ℚ) * (1 / ((m : ℚ) + 1) ^ 1) * S03.Z [2] (2 * m) +
      (-1 : ℚ) * (1 / (2 * (m : ℚ) + 1) ^ 1) * S03.Z [2] (2 * m) +
      (1 / 2 : ℚ) * (1 / ((m : ℚ) + 1) ^ 1) * S03.Z [-1, 2] (2 * m) := by
  rw [show 2 * ((m + 1 : ℕ) : ℤ) + 1 - 1 = ((2 * m + 1 : ℕ) : ℤ) + 1 by push_cast; ring, hcoef_succ]
  simp only [ddR25, Znil, S03.Z_cons_succ, gqm1, gq2, gqm2, gq3, gqm3]
  push_cast
  simp only [pow_succ, pow_mul, neg_one_sq, one_pow]
  field_simp
  ring

lemma hE55 (m : ℕ) : hcoef [1, -1, 1, -1] (2 * ((m + 1 : ℕ) : ℤ) - 0) =
      (1 / 2 : ℚ) * (1 / ((m : ℚ) + 1) ^ 1) * S03.Z [-1, -1] (2 * m) +
      (-1 : ℚ) * (1 / (2 * (m : ℚ) + 1) ^ 1) * S03.Z [-1, -1] (2 * m) +
      (1 / 2 : ℚ) * (1 / ((m : ℚ) + 1) ^ 1) * S03.Z [-1, -1, -1] (2 * m) := by
  rw [show 2 * ((m + 1 : ℕ) : ℤ) - 0 = ((2 * m + 1 : ℕ) : ℤ) + 1 by push_cast; ring, hcoef_succ]
  simp only [ddR26, Znil, S03.Z_cons_succ, gqm1, gq2, gqm2, gq3, gqm3]
  push_cast
  simp only [pow_succ, pow_mul, neg_one_sq, one_pow]
  field_simp
  ring

lemma hO55 (m : ℕ) : hcoef [1, -1, 1, -1] (2 * ((m + 1 : ℕ) : ℤ) + 1 - 0) =
      (1 / 2 : ℚ) * (1 / ((m : ℚ) + 1) ^ 1) * S03.Z [-1] (2 * m) +
      (-1 / 2 : ℚ) * (1 / (2 * (m : ℚ) + 1) ^ 1) * S03.Z [-1] (2 * m) +
      (-1 / 2 : ℚ) * (1 / (2 * (m : ℚ) + 3) ^ 1) * S03.Z [-1] (2 * m) +
      (1 / 2 : ℚ) * (1 / ((m : ℚ) + 1) ^ 1) * S03.Z [-1, -1] (2 * m) +
      (-1 / 2 : ℚ) * (1 / (2 * (m : ℚ) + 1) ^ 1) * S03.Z [-1, -1] (2 * m) +
      (-1 / 2 : ℚ) * (1 / (2 * (m : ℚ) + 3) ^ 1) * S03.Z [-1, -1] (2 * m) +
      (1 : ℚ) * (1 / (2 * (m : ℚ) + 3) ^ 1) * S03.Z [-1, -1, -1] (2 * m) := by
  rw [show 2 * ((m + 1 : ℕ) : ℤ) + 1 - 0 = ((2 * m + 1 + 1 : ℕ) : ℤ) + 1 by push_cast; ring, hcoef_succ]
  simp only [ddR26, Znil, S03.Z_cons_succ, gqm1, gq2, gqm2, gq3, gqm3]
  push_cast
  simp only [pow_succ, pow_mul, neg_one_sq, one_pow]
  field_simp
  ring

lemma hE56 (m : ℕ) : hcoef [1, -1, 1, -1] (2 * ((m + 1 : ℕ) : ℤ) - 1) =
      (1 : ℚ) * (1 / (2 * (m : ℚ) + 1) ^ 1) * S03.Z [-1, -1, -1] (2 * m) := by
  rw [show 2 * ((m + 1 : ℕ) : ℤ) - 1 = ((2 * m : ℕ) : ℤ) + 1 by push_cast; ring, hcoef_succ]
  simp only [ddR26, Znil, S03.Z_cons_succ, gqm1, gq2, gqm2, gq3, gqm3]
  push_cast
  simp only [pow_succ, pow_mul, neg_one_sq, one_pow]
  field_simp
  ring

lemma hO56 (m : ℕ) : hcoef [1, -1, 1, -1] (2 * ((m + 1 : ℕ) : ℤ) + 1 - 1) =
      (1 / 2 : ℚ) * (1 / ((m : ℚ) + 1) ^ 1) * S03.Z [-1, -1] (2 * m) +
      (-1 : ℚ) * (1 / (2 * (m : ℚ) + 1) ^ 1) * S03.Z [-1, -1] (2 * m) +
      (1 / 2 : ℚ) * (1 / ((m : ℚ) + 1) ^ 1) * S03.Z [-1, -1, -1] (2 * m) := by
  rw [show 2 * ((m + 1 : ℕ) : ℤ) + 1 - 1 = ((2 * m + 1 : ℕ) : ℤ) + 1 by push_cast; ring, hcoef_succ]
  simp only [ddR26, Znil, S03.Z_cons_succ, gqm1, gq2, gqm2, gq3, gqm3]
  push_cast
  simp only [pow_succ, pow_mul, neg_one_sq, one_pow]
  field_simp
  ring

lemma hE57 (m : ℕ) : hcoef [-1, 0, 0, 1] (2 * ((m + 1 : ℕ) : ℤ) - -1) =
      (1 / 2 : ℚ) * (1 / ((m : ℚ) + 1) ^ 1) +
      (-1 / 4 : ℚ) * (1 / ((m : ℚ) + 1) ^ 2) +
      (1 / 8 : ℚ) * (1 / ((m : ℚ) + 1) ^ 3) +
      (-1 / 8 : ℚ) * (1 / (2 * (m : ℚ) + 1) ^ 1) +
      (1 / 4 : ℚ) * (1 / (2 * (m : ℚ) + 1) ^ 2) +
      (-1 / 2 : ℚ) * (1 / (2 * (m : ℚ) + 1) ^ 3) +
      (-7 / 8 : ℚ) * (1 / (2 * (m : ℚ) + 3) ^ 1) +
      (1 : ℚ) * (1 / (2 * (m : ℚ) + 3) ^ 1) * S03.Z [-3] (2 * m) := by
  rw [show 2 * ((m + 1 : ℕ) : ℤ) - -1 = ((2 * m + 1 + 1 : ℕ) : ℤ) + 1 by push_cast; ring, hcoef_succ]
  simp only [ddR27, Znil, S03.Z_cons_succ, gqm1, gq2, gqm2, gq3, gqm3]
  push_cast
  simp only [pow_succ, pow_mul, neg_one_sq, one_pow]
  field_simp
  ring

lemma hO57 (m : ℕ) : hcoef [-1, 0, 0, 1] (2 * ((m + 1 : ℕ) : ℤ) + 1 - -1) =
      (-1 / 16 : ℚ) * (1 / ((m : ℚ) + 1) ^ 1) +
      (1 / 16 : ℚ) * (1 / ((m : ℚ) + 1) ^ 2) +
      (-1 / 16 : ℚ) * (1 / ((m : ℚ) + 1) ^ 3) +
      (1 / 27 : ℚ) * (1 / (2 * (m : ℚ) + 1) ^ 1) +
      (-1 / 9 : ℚ) * (1 / (2 * (m : ℚ) + 1) ^ 2) +
      (1 / 3 : ℚ) * (1 / (2 * (m : ℚ) + 1) ^ 3) +
      (-197 / 432 : ℚ) * (1 / ((m : ℚ) + 2) ^ 1) +
      (1 : ℚ) * (1 / (2 * (m : ℚ) + 3) ^ 1) +
      (-1 : ℚ) * (1 / (2 * (m : ℚ) + 3) ^ 2) +
      (1 : ℚ) * (1 / (2 * (m : ℚ) + 3) ^ 3) +
      (-1 / 2 : ℚ) * (1 / ((m : ℚ) + 2) ^ 1) * S03.Z [-3] (2 * m) := by
  rw [show 2 * ((m + 1 : ℕ) : ℤ) + 1 - -1 = ((2 * m + 1 + 1 + 1 : ℕ) : ℤ) + 1 by push_cast; ring, hcoef_succ]
  simp only [ddR27, Znil, S03.Z_cons_succ, gqm1, gq2, gqm2, gq3, gqm3]
  push_cast
  simp only [pow_succ, pow_mul, neg_one_sq, one_pow]
  field_simp
  ring

lemma hE58 (m : ℕ) : hcoef [-1, 0, 0, 1] (2 * ((m + 1 : ℕ) : ℤ) - 0) =
      (-1 / 2 : ℚ) * (1 / ((m : ℚ) + 1) ^ 1) +
      (1 : ℚ) * (1 / (2 * (m : ℚ) + 1) ^ 1) +
      (-1 : ℚ) * (1 / (2 * (m : ℚ) + 1) ^ 2) +
      (1 : ℚ) * (1 / (2 * (m : ℚ) + 1) ^ 3) +
      (-1 / 2 : ℚ) * (1 / ((m : ℚ) + 1) ^ 1) * S03.Z [-3] (2 * m) := by
  rw [show 2 * ((m + 1 : ℕ) : ℤ) - 0 = ((2 * m + 1 : ℕ) : ℤ) + 1 by push_cast; ring, hcoef_succ]
  simp only [ddR27, Znil, S03.Z_cons_succ, gqm1, gq2, gqm2, gq3, gqm3]
  push_cast
  simp only [pow_succ, pow_mul, neg_one_sq, one_pow]
  field_simp
  ring

lemma hO58 (m : ℕ) : hcoef [-1, 0, 0, 1] (2 * ((m + 1 : ℕ) : ℤ) + 1 - 0) =
      (1 / 2 : ℚ) * (1 / ((m : ℚ) + 1) ^ 1) +
      (-1 / 4 : ℚ) * (1 / ((m : ℚ) + 1) ^ 2) +
      (1 / 8 : ℚ) * (1 / ((m : ℚ) + 1) ^ 3) +
      (-1 / 8 : ℚ) * (1 / (2 * (m : ℚ) + 1) ^ 1) +
      (1 / 4 : ℚ) * (1 / (2 * (m : ℚ) + 1) ^ 2) +
      (-1 / 2 : ℚ) * (1 / (2 * (m : ℚ) + 1) ^ 3) +
      (-7 / 8 : ℚ) * (1 / (2 * (m : ℚ) + 3) ^ 1) +
      (1 : ℚ) * (1 / (2 * (m : ℚ) + 3) ^ 1) * S03.Z [-3] (2 * m) := by
  rw [show 2 * ((m + 1 : ℕ) : ℤ) + 1 - 0 = ((2 * m + 1 + 1 : ℕ) : ℤ) + 1 by push_cast; ring, hcoef_succ]
  simp only [ddR27, Znil, S03.Z_cons_succ, gqm1, gq2, gqm2, gq3, gqm3]
  push_cast
  simp only [pow_succ, pow_mul, neg_one_sq, one_pow]
  field_simp
  ring

lemma hE59 (m : ℕ) : hcoef [-1, 0, 0, 1] (2 * ((m + 1 : ℕ) : ℤ) - 1) =
      (1 : ℚ) * (1 / (2 * (m : ℚ) + 1) ^ 1) * S03.Z [-3] (2 * m) := by
  rw [show 2 * ((m + 1 : ℕ) : ℤ) - 1 = ((2 * m : ℕ) : ℤ) + 1 by push_cast; ring, hcoef_succ]
  simp only [ddR27, Znil, S03.Z_cons_succ, gqm1, gq2, gqm2, gq3, gqm3]
  push_cast
  simp only [pow_succ, pow_mul, neg_one_sq, one_pow]
  field_simp
  ring

lemma hO59 (m : ℕ) : hcoef [-1, 0, 0, 1] (2 * ((m + 1 : ℕ) : ℤ) + 1 - 1) =
      (-1 / 2 : ℚ) * (1 / ((m : ℚ) + 1) ^ 1) +
      (1 : ℚ) * (1 / (2 * (m : ℚ) + 1) ^ 1) +
      (-1 : ℚ) * (1 / (2 * (m : ℚ) + 1) ^ 2) +
      (1 : ℚ) * (1 / (2 * (m : ℚ) + 1) ^ 3) +
      (-1 / 2 : ℚ) * (1 / ((m : ℚ) + 1) ^ 1) * S03.Z [-3] (2 * m) := by
  rw [show 2 * ((m + 1 : ℕ) : ℤ) + 1 - 1 = ((2 * m + 1 : ℕ) : ℤ) + 1 by push_cast; ring, hcoef_succ]
  simp only [ddR27, Znil, S03.Z_cons_succ, gqm1, gq2, gqm2, gq3, gqm3]
  push_cast
  simp only [pow_succ, pow_mul, neg_one_sq, one_pow]
  field_simp
  ring

lemma hE60 (m : ℕ) : hcoef [-1, 0, 0, -1] (2 * ((m + 1 : ℕ) : ℤ) - -1) =
      (-1 / 2 : ℚ) * (1 / ((m : ℚ) + 1) ^ 1) +
      (1 / 4 : ℚ) * (1 / ((m : ℚ) + 1) ^ 2) +
      (-1 / 8 : ℚ) * (1 / ((m : ℚ) + 1) ^ 3) +
      (-1 / 8 : ℚ) * (1 / (2 * (m : ℚ) + 1) ^ 1) +
      (1 / 4 : ℚ) * (1 / (2 * (m : ℚ) + 1) ^ 2) +
      (-1 / 2 : ℚ) * (1 / (2 * (m : ℚ) + 1) ^ 3) +
      (9 / 8 : ℚ) * (1 / (2 * (m : ℚ) + 3) ^ 1) +
      (-1 : ℚ) * (1 / (2 * (m : ℚ) + 3) ^ 1) * S03.Z [3] (2 * m) := by
  rw [show 2 * ((m + 1 : ℕ) : ℤ) - -1 = ((2 * m + 1 + 1 : ℕ) : ℤ) + 1 by push_cast; ring, hcoef_succ]
  simp only [ddR28, Znil, S03.Z_cons_succ, gqm1, gq2, gqm2, gq3, gqm3]
  push_cast
  simp only [pow_succ, pow_mul, neg_one_sq, one_pow]
  field_simp
  ring

lemma hO60 (m : ℕ) : hcoef [-1, 0, 0, -1] (2 * ((m + 1 : ℕ) : ℤ) + 1 - -1) =
      (1 / 16 : ℚ) * (1 / ((m : ℚ) + 1) ^ 1) +
      (-1 / 16 : ℚ) * (1 / ((m : ℚ) + 1) ^ 2) +
      (1 / 16 : ℚ) * (1 / ((m : ℚ) + 1) ^ 3) +
      (1 / 27 : ℚ) * (1 / (2 * (m : ℚ) + 1) ^ 1) +
      (-1 / 9 : ℚ) * (1 / (2 * (m : ℚ) + 1) ^ 2) +
      (1 / 3 : ℚ) * (1 / (2 * (m : ℚ) + 1) ^ 3) +
      (-251 / 432 : ℚ) * (1 / ((m : ℚ) + 2) ^ 1) +
      (1 : ℚ) * (1 / (2 * (m : ℚ) + 3) ^ 1) +
      (-1 : ℚ) * (1 / (2 * (m : ℚ) + 3) ^ 2) +
      (1 : ℚ) * (1 / (2 * (m : ℚ) + 3) ^ 3) +
      (1 / 2 : ℚ) * (1 / ((m : ℚ) + 2) ^ 1) * S03.Z [3] (2 * m) := by
  rw [show 2 * ((m + 1 : ℕ) : ℤ) + 1 - -1 = ((2 * m + 1 + 1 + 1 : ℕ) : ℤ) + 1 by push_cast; ring, hcoef_succ]
  simp only [ddR28, Znil, S03.Z_cons_succ, gqm1, gq2, gqm2, gq3, gqm3]
  push_cast
  simp only [pow_succ, pow_mul, neg_one_sq, one_pow]
  field_simp
  ring

lemma hE61 (m : ℕ) : hcoef [-1, 0, 0, -1] (2 * ((m + 1 : ℕ) : ℤ) - 1) =
      (-1 : ℚ) * (1 / (2 * (m : ℚ) + 1) ^ 1) * S03.Z [3] (2 * m) := by
  rw [show 2 * ((m + 1 : ℕ) : ℤ) - 1 = ((2 * m : ℕ) : ℤ) + 1 by push_cast; ring, hcoef_succ]
  simp only [ddR28, Znil, S03.Z_cons_succ, gqm1, gq2, gqm2, gq3, gqm3]
  push_cast
  simp only [pow_succ, pow_mul, neg_one_sq, one_pow]
  field_simp
  ring

lemma hO61 (m : ℕ) : hcoef [-1, 0, 0, -1] (2 * ((m + 1 : ℕ) : ℤ) + 1 - 1) =
      (-1 / 2 : ℚ) * (1 / ((m : ℚ) + 1) ^ 1) +
      (1 : ℚ) * (1 / (2 * (m : ℚ) + 1) ^ 1) +
      (-1 : ℚ) * (1 / (2 * (m : ℚ) + 1) ^ 2) +
      (1 : ℚ) * (1 / (2 * (m : ℚ) + 1) ^ 3) +
      (1 / 2 : ℚ) * (1 / ((m : ℚ) + 1) ^ 1) * S03.Z [3] (2 * m) := by
  rw [show 2 * ((m + 1 : ℕ) : ℤ) + 1 - 1 = ((2 * m + 1 : ℕ) : ℤ) + 1 by push_cast; ring, hcoef_succ]
  simp only [ddR28, Znil, S03.Z_cons_succ, gqm1, gq2, gqm2, gq3, gqm3]
  push_cast
  simp only [pow_succ, pow_mul, neg_one_sq, one_pow]
  field_simp
  ring

lemma hE62 (m : ℕ) : hcoef [-1, 0, 1, -1] (2 * ((m + 1 : ℕ) : ℤ) - -1) =
      (-1 / 4 : ℚ) * (1 / ((m : ℚ) + 1) ^ 2) +
      (1 / 2 : ℚ) * (1 / (2 * (m : ℚ) + 1) ^ 1) +
      (-1 / 2 : ℚ) * (1 / (2 * (m : ℚ) + 3) ^ 1) +
      (1 / 2 : ℚ) * (1 / ((m : ℚ) + 1) ^ 1) * S03.Z [-1] (2 * m) +
      (-1 / 4 : ℚ) * (1 / ((m : ℚ) + 1) ^ 2) * S03.Z [-1] (2 * m) +
      (-1 / 4 : ℚ) * (1 / (2 * (m : ℚ) + 1) ^ 1) * S03.Z [-1] (2 * m) +
      (1 / 2 : ℚ) * (1 / (2 * (m : ℚ) + 1) ^ 2) * S03.Z [-1] (2 * m) +
      (-3 / 4 : ℚ) * (1 / (2 * (m : ℚ) + 3) ^ 1) * S03.Z [-1] (2 * m) +
      (-1 : ℚ) * (1 / (2 * (m : ℚ) + 3) ^ 1) * S03.Z [-2, -1] (2 * m) := by
  rw [show 2 * ((m + 1 : ℕ) : ℤ) - -1 = ((2 * m + 1 + 1 : ℕ) : ℤ) + 1 by push_cast; ring, hcoef_succ]
  simp only [ddR29, Znil, S03.Z_cons_succ, gqm1, gq2, gqm2, gq3, gqm3]
  push_cast
  simp only [pow_succ, pow_mul, neg_one_sq, one_pow]
  field_simp
  ring

lemma hO62 (m : ℕ) : hcoef [-1, 0, 1, -1] (2 * ((m + 1 : ℕ) : ℤ) + 1 - -1) =
      (-1 / 8 : ℚ) * (1 / ((m : ℚ) + 1) ^ 1) +
      (1 / 8 : ℚ) * (1 / ((m : ℚ) + 1) ^ 2) +
      (-1 / 4 : ℚ) * (1 / (2 * (m : ℚ) + 1) ^ 1) +
      (1 / 8 : ℚ) * (1 / ((m : ℚ) + 2) ^ 1) +
      (1 / 4 : ℚ) * (1 / (2 * (m : ℚ) + 3) ^ 1) +
      (1 / 2 : ℚ) * (1 / (2 * (m : ℚ) + 3) ^ 2) +
      (-1 / 8 : ℚ) * (1 / ((m : ℚ) + 1) ^ 1) * S03.Z [-1] (2 * m) +
      (1 / 8 : ℚ) * (1 / ((m : ℚ) + 1) ^ 2) * S03.Z [-1] (2 * m) +
      (1 / 9 : ℚ) * (1 / (2 * (m : ℚ) + 1) ^ 1) * S03.Z [-1] (2 * m) +
      (-1 / 3 : ℚ) * (1 / (2 * (m : ℚ) + 1) ^ 2) * S03.Z [-1] (2 * m) +
      (-31 / 72 : ℚ) * (1 / ((m : ℚ) + 2) ^ 1) * S03.Z [-1] (2 * m) +
      (1 : ℚ) * (1 / (2 * (m : ℚ) + 3) ^ 1) * S03.Z [-1] (2 * m) +
      (-1 : ℚ) * (1 / (2 * (m : ℚ) + 3) ^ 2) * S03.Z [-1] (2 * m) +
      (1 / 2 : ℚ) * (1 / ((m : ℚ) + 2) ^ 1) * S03.Z [-2, -1] (2 * m) := by
  rw [show 2 * ((m + 1 : ℕ) : ℤ) + 1 - -1 = ((2 * m + 1 + 1 + 1 : ℕ) : ℤ) + 1 by push_cast; ring, hcoef_succ]
  simp only [ddR29, Znil, S03.Z_cons_succ, gqm1, gq2, gqm2, gq3, gqm3]
  push_cast
  simp only [pow_succ, pow_mul, neg_one_sq, one_pow]
  field_simp
  ring

lemma hE63 (m : ℕ) : hcoef [-1, 0, 1, -1] (2 * ((m + 1 : ℕ) : ℤ) - 0) =
      (-1 / 2 : ℚ) * (1 / ((m : ℚ) + 1) ^ 1) * S03.Z [-1] (2 * m) +
      (1 : ℚ) * (1 / (2 * (m : ℚ) + 1) ^ 1) * S03.Z [-1] (2 * m) +
      (-1 : ℚ) * (1 / (2 * (m : ℚ) + 1) ^ 2) * S03.Z [-1] (2 * m) +
      (1 / 2 : ℚ) * (1 / ((m : ℚ) + 1) ^ 1) * S03.Z [-2, -1] (2 * m) := by
  rw [show 2 * ((m + 1 : ℕ) : ℤ) - 0 = ((2 * m + 1 : ℕ) : ℤ) + 1 by push_cast; ring, hcoef_succ]
  simp only [ddR29, Znil, S03.Z_cons_succ, gqm1, gq2, gqm2, gq3, gqm3]
  push_cast
  simp only [pow_succ, pow_mul, neg_one_sq, one_pow]
  field_simp
  ring

lemma hO63 (m : ℕ) : hcoef [-1, 0, 1, -1] (2 * ((m + 1 : ℕ) : ℤ) + 1 - 0) =
      (-1 / 4 : ℚ) * (1 / ((m : ℚ) + 1) ^ 2) +
      (1 / 2 : ℚ) * (1 / (2 * (m : ℚ) + 1) ^ 1) +
      (-1 / 2 : ℚ) * (1 / (2 * (m : ℚ) + 3) ^ 1) +
      (1 / 2 : ℚ) * (1 / ((m : ℚ) + 1) ^ 1) * S03.Z [-1] (2 * m) +
      (-1 / 4 : ℚ) * (1 / ((m : ℚ) + 1) ^ 2) * S03.Z [-1] (2 * m) +
      (-1 / 4 : ℚ) * (1 / (2 * (m : ℚ) + 1) ^ 1) * S03.Z [-1] (2 * m) +
      (1 / 2 : ℚ) * (1 / (2 * (m : ℚ) + 1) ^ 2) * S03.Z [-1] (2 * m) +
      (-3 / 4 : ℚ) * (1 / (2 * (m : ℚ) + 3) ^ 1) * S03.Z [-1] (2 * m) +
      (-1 : ℚ) * (1 / (2 * (m : ℚ) + 3) ^ 1) * S03.Z [-2, -1] (2 * m) := by
  rw [show 2 * ((m + 1 : ℕ) : ℤ) + 1 - 0 = ((2 * m + 1 + 1 : ℕ) : ℤ) + 1 by push_cast; ring, hcoef_succ]
  simp only [ddR29, Znil, S03.Z_cons_succ, gqm1, gq2, gqm2, gq3, gqm3]
  push_cast
  simp only [pow_succ, pow_mul, neg_one_sq, one_pow]
  field_simp
  ring

lemma hE64 (m : ℕ) : hcoef [-1, 0, 1, -1] (2 * ((m + 1 : ℕ) : ℤ) - 1) =
      (-1 : ℚ) * (1 / (2 * (m : ℚ) + 1) ^ 1) * S03.Z [-2, -1] (2 * m) := by
  rw [show 2 * ((m + 1 : ℕ) : ℤ) - 1 = ((2 * m : ℕ) : ℤ) + 1 by push_cast; ring, hcoef_succ]
  simp only [ddR29, Znil, S03.Z_cons_succ, gqm1, gq2, gqm2, gq3, gqm3]
  push_cast
  simp only [pow_succ, pow_mul, neg_one_sq, one_pow]
  field_simp
  ring

lemma hO64 (m : ℕ) : hcoef [-1, 0, 1, -1] (2 * ((m + 1 : ℕ) : ℤ) + 1 - 1) =
      (-1 / 2 : ℚ) * (1 / ((m : ℚ) + 1) ^ 1) * S03.Z [-1] (2 * m) +
      (1 : ℚ) * (1 / (2 * (m : ℚ) + 1) ^ 1) * S03.Z [-1] (2 * m) +
      (-1 : ℚ) * (1 / (2 * (m : ℚ) + 1) ^ 2) * S03.Z [-1] (2 * m) +
      (1 / 2 : ℚ) * (1 / ((m : ℚ) + 1) ^ 1) * S03.Z [-2, -1] (2 * m) := by
  rw [show 2 * ((m + 1 : ℕ) : ℤ) + 1 - 1 = ((2 * m + 1 : ℕ) : ℤ) + 1 by push_cast; ring, hcoef_succ]
  simp only [ddR29, Znil, S03.Z_cons_succ, gqm1, gq2, gqm2, gq3, gqm3]
  push_cast
  simp only [pow_succ, pow_mul, neg_one_sq, one_pow]
  field_simp
  ring

lemma hE65 (m : ℕ) : hcoef [-1, 0, -1, 1] (2 * ((m + 1 : ℕ) : ℤ) - -1) =
      (-1 / 4 : ℚ) * (1 / ((m : ℚ) + 1) ^ 2) +
      (1 / 2 : ℚ) * (1 / (2 * (m : ℚ) + 1) ^ 1) +
      (-1 / 2 : ℚ) * (1 / (2 * (m : ℚ) + 3) ^ 1) +
      (-1 : ℚ) * (1 / (2 * (m : ℚ) + 3) ^ 1) * S03.Z [2, -1] (2 * m) +
      (1 / 2 : ℚ) * (1 / ((m : ℚ) + 1) ^ 1) * S03.Z [-1] (2 * m) +
      (-1 / 4 : ℚ) * (1 / ((m : ℚ) + 1) ^ 2) * S03.Z [-1] (2 * m) +
      (1 / 4 : ℚ) * (1 / (2 * (m : ℚ) + 1) ^ 1) * S03.Z [-1] (2 * m) +
      (-1 / 2 : ℚ) * (1 / (2 * (m : ℚ) + 1) ^ 2) * S03.Z [-1] (2 * m) +
      (-5 / 4 : ℚ) * (1 / (2 * (m : ℚ) + 3) ^ 1) * S03.Z [-1] (2 * m) := by
  rw [show 2 * ((m + 1 : ℕ) : ℤ) - -1 = ((2 * m + 1 + 1 : ℕ) : ℤ) + 1 by push_cast; ring, hcoef_succ]
  simp only [ddR30, Znil, S03.Z_cons_succ, gqm1, gq2, gqm2, gq3, gqm3]
  push_cast
  simp only [pow_succ, pow_mul, neg_one_sq, one_pow]
  field_simp
  ring

lemma hO65 (m : ℕ) : hcoef [-1, 0, -1, 1] (2 * ((m + 1 : ℕ) : ℤ) + 1 - -1) =
      (3 / 8 : ℚ) * (1 / ((m : ℚ) + 1) ^ 1) +
      (1 / 8 : ℚ) * (1 / ((m : ℚ) + 1) ^ 2) +
      (-5 / 12 : ℚ) * (1 / (2 * (m : ℚ) + 1) ^ 1) +
      (-1 / 24 : ℚ) * (1 / ((m : ℚ) + 2) ^ 1) +
      (-1 / 4 : ℚ) * (1 / (2 * (m : ℚ) + 3) ^ 1) +
      (-1 / 2 : ℚ) * (1 / (2 * (m : ℚ) + 3) ^ 2) +
      (1 / 2 : ℚ) * (1 / ((m : ℚ) + 2) ^ 1) * S03.Z [2, -1] (2 * m) +
      (-1 / 8 : ℚ) * (1 / ((m : ℚ) + 1) ^ 1) * S03.Z [-1] (2 * m) +
      (1 / 8 : ℚ) * (1 / ((m : ℚ) + 1) ^ 2) * S03.Z [-1] (2 * m) +
      (-1 / 9 : ℚ) * (1 / (2 * (m : ℚ) + 1) ^ 1) * S03.Z [-1] (2 * m) +
      (1 / 3 : ℚ) * (1 / (2 * (m : ℚ) + 1) ^ 2) * S03.Z [-1] (2 * m) +
      (49 / 72 : ℚ) * (1 / ((m : ℚ) + 2) ^ 1) * S03.Z [-1] (2 * m) +
      (-1 : ℚ) * (1 / (2 * (m : ℚ) + 3) ^ 1) * S03.Z [-1] (2 * m) +
      (1 : ℚ) * (1 / (2 * (m : ℚ) + 3) ^ 2) * S03.Z [-1] (2 * m) := by
  rw [show 2 * ((m + 1 : ℕ) : ℤ) + 1 - -1 = ((2 * m + 1 + 1 + 1 : ℕ) : ℤ) + 1 by push_cast; ring, hcoef_succ]
  simp only [ddR30, Znil, S03.Z_cons_succ, gqm1, gq2, gqm2, gq3, gqm3]
  push_cast
  simp only [pow_succ, pow_mul, neg_one_sq, one_pow]
  field_simp
  ring

lemma hE66 (m : ℕ) : hcoef [-1, 0, -1, 1] (2 * ((m + 1 : ℕ) : ℤ) - 1) =
      (-1 : ℚ) * (1 / (2 * (m : ℚ) + 1) ^ 1) * S03.Z [2, -1] (2 * m) := by
  rw [show 2 * ((m + 1 : ℕ) : ℤ) - 1 = ((2 * m : ℕ) : ℤ) + 1 by push_cast; ring, hcoef_succ]
  simp only [ddR30, Znil, S03.Z_cons_succ, gqm1, gq2, gqm2, gq3, gqm3]
  push_cast
  simp only [pow_succ, pow_mul, neg_one_sq, one_pow]
  field_simp
  ring

lemma hO66 (m : ℕ) : hcoef [-1, 0, -1, 1] (2 * ((m + 1 : ℕ) : ℤ) + 1 - 1) =
      (1 / 2 : ℚ) * (1 / ((m : ℚ) + 1) ^ 1) * S03.Z [2, -1] (2 * m) +
      (1 / 2 : ℚ) * (1 / ((m : ℚ) + 1) ^ 1) * S03.Z [-1] (2 * m) +
      (-1 : ℚ) * (1 / (2 * (m : ℚ) + 1) ^ 1) * S03.Z [-1] (2 * m) +
      (1 : ℚ) * (1 / (2 * (m : ℚ) + 1) ^ 2) * S03.Z [-1] (2 * m) := by
  rw [show 2 * ((m + 1 : ℕ) : ℤ) + 1 - 1 = ((2 * m + 1 : ℕ) : ℤ) + 1 by push_cast; ring, hcoef_succ]
  simp only [ddR30, Znil, S03.Z_cons_succ, gqm1, gq2, gqm2, gq3, gqm3]
  push_cast
  simp only [pow_succ, pow_mul, neg_one_sq, one_pow]
  field_simp
  ring

lemma hE67 (m : ℕ) : hcoef [-1, 1, 0, 1] (2 * ((m + 1 : ℕ) : ℤ) - -1) =
      (1 / 2 : ℚ) * (1 / ((m : ℚ) + 1) ^ 1) +
      (-3 / 4 : ℚ) * (1 / (2 * (m : ℚ) + 1) ^ 1) +
      (1 / 2 : ℚ) * (1 / (2 * (m : ℚ) + 1) ^ 2) +
      (-1 / 4 : ℚ) * (1 / (2 * (m : ℚ) + 3) ^ 1) +
      (1 / 2 : ℚ) * (1 / ((m : ℚ) + 1) ^ 1) * S03.Z [2] (2 * m) +
      (-1 / 2 : ℚ) * (1 / (2 * (m : ℚ) + 1) ^ 1) * S03.Z [2] (2 * m) +
      (-1 / 2 : ℚ) * (1 / (2 * (m : ℚ) + 3) ^ 1) * S03.Z [2] (2 * m) +
      (1 : ℚ) * (1 / (2 * (m : ℚ) + 3) ^ 1) * S03.Z [-1, 2] (2 * m) := by
  rw [show 2 * ((m + 1 : ℕ) : ℤ) - -1 = ((2 * m + 1 + 1 : ℕ) : ℤ) + 1 by push_cast; ring, hcoef_succ]
  simp only [ddR31, Znil, S03.Z_cons_succ, gqm1, gq2, gqm2, gq3, gqm3]
  push_cast
  simp only [pow_succ, pow_mul, neg_one_sq, one_pow]
  field_simp
  ring

lemma hO67 (m : ℕ) : hcoef [-1, 1, 0, 1] (2 * ((m + 1 : ℕ) : ℤ) + 1 - -1) =
      (-5 / 8 : ℚ) * (1 / ((m : ℚ) + 1) ^ 1) +
      (1 / 8 : ℚ) * (1 / ((m : ℚ) + 1) ^ 2) +
      (11 / 36 : ℚ) * (1 / (2 * (m : ℚ) + 1) ^ 1) +
      (-1 / 6 : ℚ) * (1 / (2 * (m : ℚ) + 1) ^ 2) +
      (-11 / 72 : ℚ) * (1 / ((m : ℚ) + 2) ^ 1) +
      (5 / 4 : ℚ) * (1 / (2 * (m : ℚ) + 3) ^ 1) +
      (-1 / 4 : ℚ) * (1 / ((m : ℚ) + 1) ^ 1) * S03.Z [2] (2 * m) +
      (1 / 3 : ℚ) * (1 / (2 * (m : ℚ) + 1) ^ 1) * S03.Z [2] (2 * m) +
      (-5 / 12 : ℚ) * (1 / ((m : ℚ) + 2) ^ 1) * S03.Z [2] (2 * m) +
      (1 : ℚ) * (1 / (2 * (m : ℚ) + 3) ^ 1) * S03.Z [2] (2 * m) +
      (-1 / 2 : ℚ) * (1 / ((m : ℚ) + 2) ^ 1) * S03.Z [-1, 2] (2 * m) := by
  rw [show 2 * ((m + 1 : ℕ) : ℤ) + 1 - -1 = ((2 * m + 1 + 1 + 1 : ℕ) : ℤ) + 1 by push_cast; ring, hcoef_succ]
  simp only [ddR31, Znil, S03.Z_cons_succ, gqm1, gq2, gqm2, gq3, gqm3]
  push_cast
  simp only [pow_succ, pow_mul, neg_one_sq, one_pow]
  field_simp
  ring

lemma hE68 (m : ℕ) : hcoef [-1, 1, 0, 1] (2 * ((m + 1 : ℕ) : ℤ) - 1) =
      (1 : ℚ) * (1 / (2 * (m : ℚ) + 1) ^ 1) * S03.Z [-1, 2] (2 * m) := by
  rw [show 2 * ((m + 1 : ℕ) : ℤ) - 1 = ((2 * m : ℕ) : ℤ) + 1 by push_cast; ring, hcoef_succ]
  simp only [ddR31, Znil, S03.Z_cons_succ, gqm1, gq2, gqm2, gq3, gqm3]
  push_cast
  simp only [pow_succ, pow_mul, neg_one_sq, one_pow]
  field_simp
  ring

lemma hO68 (m : ℕ) : hcoef [-1, 1, 0, 1] (2 * ((m + 1 : ℕ) : ℤ) + 1 - 1) =
      (-1 / 2 : ℚ) * (1 / ((m : ℚ) + 1) ^ 1) * S03.Z [2] (2 * m) +
      (1 : ℚ) * (1 / (2 * (m : ℚ) + 1) ^ 1) * S03.Z [2] (2 * m) +
      (-1 / 2 : ℚ) * (1 / ((m : ℚ) + 1) ^ 1) * S03.Z [-1, 2] (2 * m) := by
  rw [show 2 * ((m + 1 : ℕ) : ℤ) + 1 - 1 = ((2 * m + 1 : ℕ) : ℤ) + 1 by push_cast; ring, hcoef_succ]
  simp only [ddR31, Znil, S03.Z_cons_succ, gqm1, gq2, gqm2, gq3, gqm3]
  push_cast
  simp only [pow_succ, pow_mul, neg_one_sq, one_pow]
  field_simp
  ring

lemma hE69 (m : ℕ) : hcoef [-1, 1, 0, -1] (2 * ((m + 1 : ℕ) : ℤ) - -1) =
      (1 / 2 : ℚ) * (1 / ((m : ℚ) + 1) ^ 1) +
      (-3 / 4 : ℚ) * (1 / (2 * (m : ℚ) + 1) ^ 1) +
      (1 / 2 : ℚ) * (1 / (2 * (m : ℚ) + 1) ^ 2) +
      (-1 / 4 : ℚ) * (1 / (2 * (m : ℚ) + 3) ^ 1) +
      (-1 : ℚ) * (1 / (2 * (m : ℚ) + 3) ^ 1) * S03.Z [-1, -2] (2 * m) +
      (-1 / 2 : ℚ) * (1 / ((m : ℚ) + 1) ^ 1) * S03.Z [-2] (2 * m) +
      (1 / 2 : ℚ) * (1 / (2 * (m : ℚ) + 1) ^ 1) * S03.Z [-2] (2 * m) +
      (1 / 2 : ℚ) * (1 / (2 * (m : ℚ) + 3) ^ 1) * S03.Z [-2] (2 * m) := by
  rw [show 2 * ((m + 1 : ℕ) : ℤ) - -1 = ((2 * m + 1 + 1 : ℕ) : ℤ) + 1 by push_cast; ring, hcoef_succ]
  simp only [ddR32, Znil, S03.Z_cons_succ, gqm1, gq2, gqm2, gq3, gqm3]
  push_cast
  simp only [pow_succ, pow_mul, neg_one_sq, one_pow]
  field_simp
  ring

lemma hO69 (m : ℕ) : hcoef [-1, 1, 0, -1] (2 * ((m + 1 : ℕ) : ℤ) + 1 - -1) =
      (1 / 8 : ℚ) * (1 / ((m : ℚ) + 1) ^ 1) +
      (-1 / 8 : ℚ) * (1 / ((m : ℚ) + 1) ^ 2) +
      (11 / 36 : ℚ) * (1 / (2 * (m : ℚ) + 1) ^ 1) +
      (-1 / 6 : ℚ) * (1 / (2 * (m : ℚ) + 1) ^ 2) +
      (7 / 72 : ℚ) * (1 / ((m : ℚ) + 2) ^ 1) +
      (-3 / 4 : ℚ) * (1 / (2 * (m : ℚ) + 3) ^ 1) +
      (1 / 2 : ℚ) * (1 / ((m : ℚ) + 2) ^ 1) * S03.Z [-1, -2] (2 * m) +
      (1 / 4 : ℚ) * (1 / ((m : ℚ) + 1) ^ 1) * S03.Z [-2] (2 * m) +
      (-1 / 3 : ℚ) * (1 / (2 * (m : ℚ) + 1) ^ 1) * S03.Z [-2] (2 * m) +
      (5 / 12 : ℚ) * (1 / ((m : ℚ) + 2) ^ 1) * S03.Z [-2] (2 * m) +
      (-1 : ℚ) * (1 / (2 * (m : ℚ) + 3) ^ 1) * S03.Z [-2] (2 * m) := by
  rw [show 2 * ((m + 1 : ℕ) : ℤ) + 1 - -1 = ((2 * m + 1 + 1 + 1 : ℕ) : ℤ) + 1 by push_cast; ring, hcoef_succ]
  simp only [ddR32, Znil, S03.Z_cons_succ, gqm1, gq2, gqm2, gq3, gqm3]
  push_cast
  simp only [pow_succ, pow_mul, neg_one_sq, one_pow]
  field_simp
  ring

lemma hE70 (m : ℕ) : hcoef [-1, 1, 0, -1] (2 * ((m + 1 : ℕ) : ℤ) - 0) =
      (1 / 2 : ℚ) * (1 / ((m : ℚ) + 1) ^ 1) * S03.Z [-1, -2] (2 * m) +
      (1 / 2 : ℚ) * (1 / ((m : ℚ) + 1) ^ 1) * S03.Z [-2] (2 * m) +
      (-1 : ℚ) * (1 / (2 * (m : ℚ) + 1) ^ 1) * S03.Z [-2] (2 * m) := by
  rw [show 2 * ((m + 1 : ℕ) : ℤ) - 0 = ((2 * m + 1 : ℕ) : ℤ) + 1 by push_cast; ring, hcoef_succ]
  simp only [ddR32, Znil, S03.Z_cons_succ, gqm1, gq2, gqm2, gq3, gqm3]
  push_cast
  simp only [pow_succ, pow_mul, neg_one_sq, one_pow]
  field_simp
  ring

lemma hO70 (m : ℕ) : hcoef [-1, 1, 0, -1] (2 * ((m + 1 : ℕ) : ℤ) + 1 - 0) =
      (1 / 2 : ℚ) * (1 / ((m : ℚ) + 1) ^ 1) +
      (-3 / 4 : ℚ) * (1 / (2 * (m : ℚ) + 1) ^ 1) +
      (1 / 2 : ℚ) * (1 / (2 * (m : ℚ) + 1) ^ 2) +
      (-1 / 4 : ℚ) * (1 / (2 * (m : ℚ) + 3) ^ 1) +
      (-1 : ℚ) * (1 / (2 * (m : ℚ) + 3) ^ 1) * S03.Z [-1, -2] (2 * m) +
      (-1 / 2 : ℚ) * (1 / ((m : ℚ) + 1) ^ 1) * S03.Z [-2] (2 * m) +
      (1 / 2 : ℚ) * (1 / (2 * (m : ℚ) + 1) ^ 1) * S03.Z [-2] (2 * m) +
      (1 / 2 : ℚ) * (1 / (2 * (m : ℚ) + 3) ^ 1) * S03.Z [-2] (2 * m) := by
  rw [show 2 * ((m + 1 : ℕ) : ℤ) + 1 - 0 = ((2 * m + 1 + 1 : ℕ) : ℤ) + 1 by push_cast; ring, hcoef_succ]
  simp only [ddR32, Znil, S03.Z_cons_succ, gqm1, gq2, gqm2, gq3, gqm3]
  push_cast
  simp only [pow_succ, pow_mul, neg_one_sq, one_pow]
  field_simp
  ring

lemma hE71 (m : ℕ) : hcoef [-1, 1, 0, -1] (2 * ((m + 1 : ℕ) : ℤ) - 1) =
      (-1 : ℚ) * (1 / (2 * (m : ℚ) + 1) ^ 1) * S03.Z [-1, -2] (2 * m) := by
  rw [show 2 * ((m + 1 : ℕ) : ℤ) - 1 = ((2 * m : ℕ) : ℤ) + 1 by push_cast; ring, hcoef_succ]
  simp only [ddR32, Znil, S03.Z_cons_succ, gqm1, gq2, gqm2, gq3, gqm3]
  push_cast
  simp only [pow_succ, pow_mul, neg_one_sq, one_pow]
  field_simp
  ring

lemma hO71 (m : ℕ) : hcoef [-1, 1, 0, -1] (2 * ((m + 1 : ℕ) : ℤ) + 1 - 1) =
      (1 / 2 : ℚ) * (1 / ((m : ℚ) + 1) ^ 1) * S03.Z [-1, -2] (2 * m) +
      (1 / 2 : ℚ) * (1 / ((m : ℚ) + 1) ^ 1) * S03.Z [-2] (2 * m) +
      (-1 : ℚ) * (1 / (2 * (m : ℚ) + 1) ^ 1) * S03.Z [-2] (2 * m) := by
  rw [show 2 * ((m + 1 : ℕ) : ℤ) + 1 - 1 = ((2 * m + 1 : ℕ) : ℤ) + 1 by push_cast; ring, hcoef_succ]
  simp only [ddR32, Znil, S03.Z_cons_succ, gqm1, gq2, gqm2, gq3, gqm3]
  push_cast
  simp only [pow_succ, pow_mul, neg_one_sq, one_pow]
  field_simp
  ring

lemma hE72 (m : ℕ) : hcoef [-1, 1, -1, 1] (2 * ((m + 1 : ℕ) : ℤ) - 0) =
      (1 / 2 : ℚ) * (1 / ((m : ℚ) + 1) ^ 1) * S03.Z [-1, -1] (2 * m) +
      (-1 : ℚ) * (1 / (2 * (m : ℚ) + 1) ^ 1) * S03.Z [-1, -1] (2 * m) +
      (1 / 2 : ℚ) * (1 / ((m : ℚ) + 1) ^ 1) * S03.Z [-1, -1, -1] (2 * m) := by
  rw [show 2 * ((m + 1 : ℕ) : ℤ) - 0 = ((2 * m + 1 : ℕ) : ℤ) + 1 by push_cast; ring, hcoef_succ]
  simp only [ddR33, Znil, S03.Z_cons_succ, gqm1, gq2, gqm2, gq3, gqm3]
  push_cast
  simp only [pow_succ, pow_mul, neg_one_sq, one_pow]
  field_simp
  ring

lemma hO72 (m : ℕ) : hcoef [-1, 1, -1, 1] (2 * ((m + 1 : ℕ) : ℤ) + 1 - 0) =
      (-1 / 2 : ℚ) * (1 / ((m : ℚ) + 1) ^ 1) * S03.Z [-1] (2 * m) +
      (1 / 2 : ℚ) * (1 / (2 * (m : ℚ) + 1) ^ 1) * S03.Z [-1] (2 * m) +
      (1 / 2 : ℚ) * (1 / (2 * (m : ℚ) + 3) ^ 1) * S03.Z [-1] (2 * m) +
      (-1 / 2 : ℚ) * (1 / ((m : ℚ) + 1) ^ 1) * S03.Z [-1, -1] (2 * m) +
      (1 / 2 : ℚ) * (1 / (2 * (m : ℚ) + 1) ^ 1) * S03.Z [-1, -1] (2 * m) +
      (1 / 2 : ℚ) * (1 / (2 * (m : ℚ) + 3) ^ 1) * S03.Z [-1, -1] (2 * m) +
      (-1 : ℚ) * (1 / (2 * (m : ℚ) + 3) ^ 1) * S03.Z [-1, -1, -1] (2 * m) := by
  rw [show 2 * ((m + 1 : ℕ) : ℤ) + 1 - 0 = ((2 * m + 1 + 1 : ℕ) : ℤ) + 1 by push_cast; ring, hcoef_succ]
  simp only [ddR33, Znil, S03.Z_cons_succ, gqm1, gq2, gqm2, gq3, gqm3]
  push_cast
  simp only [pow_succ, pow_mul, neg_one_sq, one_pow]
  field_simp
  ring

lemma hE73 (m : ℕ) : hcoef [-1, 1, -1, 1] (2 * ((m + 1 : ℕ) : ℤ) - 1) =
      (-1 : ℚ) * (1 / (2 * (m : ℚ) + 1) ^ 1) * S03.Z [-1, -1, -1] (2 * m) := by
  rw [show 2 * ((m + 1 : ℕ) : ℤ) - 1 = ((2 * m : ℕ) : ℤ) + 1 by push_cast; ring, hcoef_succ]
  simp only [ddR33, Znil, S03.Z_cons_succ, gqm1, gq2, gqm2, gq3, gqm3]
  push_cast
  simp only [pow_succ, pow_mul, neg_one_sq, one_pow]
  field_simp
  ring

lemma hO73 (m : ℕ) : hcoef [-1, 1, -1, 1] (2 * ((m + 1 : ℕ) : ℤ) + 1 - 1) =
      (1 / 2 : ℚ) * (1 / ((m : ℚ) + 1) ^ 1) * S03.Z [-1, -1] (2 * m) +
      (-1 : ℚ) * (1 / (2 * (m : ℚ) + 1) ^ 1) * S03.Z [-1, -1] (2 * m) +
      (1 / 2 : ℚ) * (1 / ((m : ℚ) + 1) ^ 1) * S03.Z [-1, -1, -1] (2 * m) := by
  rw [show 2 * ((m + 1 : ℕ) : ℤ) + 1 - 1 = ((2 * m + 1 : ℕ) : ℤ) + 1 by push_cast; ring, hcoef_succ]
  simp only [ddR33, Znil, S03.Z_cons_succ, gqm1, gq2, gqm2, gq3, gqm3]
  push_cast
  simp only [pow_succ, pow_mul, neg_one_sq, one_pow]
  field_simp
  ring

set_option maxRecDepth 20000 in
lemma T_0_0 (m : ℕ) : Tji 0 0 m =
      (1 / 4 : ℝ) * (1 / (2 * (m : ℝ) + 1) ^ 2) +
      (1 / 2 : ℝ) * (1 / (2 * (m : ℝ) + 1) ^ 3) +
      (1 / 4 : ℝ) * (1 / (2 * (m : ℝ) + 1) ^ 4) +
      (-1 / 8 : ℝ) * ((S03.Z [2] (2 * m) : ℚ) : ℝ) +
      (-1 / 4 : ℝ) * (1 / (2 * (m : ℝ) + 1) ^ 1) * ((S03.Z [2] (2 * m) : ℚ) : ℝ) +
      (-1 / 8 : ℝ) * (1 / (2 * (m : ℝ) + 1) ^ 2) * ((S03.Z [2] (2 * m) : ℚ) : ℝ) +
      (3 / 8 : ℝ) * ((S03.Z [2] (2 * m) : ℚ) : ℝ) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) +
      (1 / 4 : ℝ) * ((m : ℝ) ^ 1) * ((S03.Z [2] (2 * m) : ℚ) : ℝ) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) +
      (1 / 8 : ℝ) * (1 / (2 * (m : ℝ) + 1) ^ 1) * ((S03.Z [2] (2 * m) : ℚ) : ℝ) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) +
      (-3 / 4 : ℝ) * ((S03.Z [3] (2 * m) : ℚ) : ℝ) +
      (-1 / 2 : ℝ) * ((m : ℝ) ^ 1) * ((S03.Z [3] (2 * m) : ℚ) : ℝ) +
      (-1 / 4 : ℝ) * (1 / (2 * (m : ℝ) + 1) ^ 1) * ((S03.Z [3] (2 * m) : ℚ) : ℝ) +
      (-1 / 4 : ℝ) * (1 / (2 * (m : ℝ) + 1) ^ 1) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) +
      (-1 / 2 : ℝ) * (1 / (2 * (m : ℝ) + 1) ^ 2) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) +
      (-1 / 4 : ℝ) * (1 / (2 * (m : ℝ) + 1) ^ 3) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) +
      (1 / 8 : ℝ) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) +
      (1 / 4 : ℝ) * (1 / (2 * (m : ℝ) + 1) ^ 1) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) +
      (1 / 8 : ℝ) * (1 / (2 * (m : ℝ) + 1) ^ 2) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) +
      (-1 / 8 : ℝ) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) +
      (-1 / 12 : ℝ) * ((m : ℝ) ^ 1) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) +
      (-1 / 24 : ℝ) * (1 / (2 * (m : ℝ) + 1) ^ 1) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) +
      (3 / 4 : ℝ) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) * ((S03.Z [-2] (2 * m) : ℚ) : ℝ) +
      (1 / 2 : ℝ) * ((m : ℝ) ^ 1) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) * ((S03.Z [-2] (2 * m) : ℚ) : ℝ) +
      (1 / 4 : ℝ) * (1 / (2 * (m : ℝ) + 1) ^ 1) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) * ((S03.Z [-2] (2 * m) : ℚ) : ℝ) +
      (-1 / 4 : ℝ) * ((S03.Z [-2] (2 * m) : ℚ) : ℝ) +
      (-1 / 2 : ℝ) * (1 / (2 * (m : ℝ) + 1) ^ 1) * ((S03.Z [-2] (2 * m) : ℚ) : ℝ) +
      (-1 / 4 : ℝ) * (1 / (2 * (m : ℝ) + 1) ^ 2) * ((S03.Z [-2] (2 * m) : ℚ) : ℝ) +
      (-1 : ℝ) * ((S03.Z [-3] (2 * m) : ℚ) : ℝ) +
      (-2 / 3 : ℝ) * ((m : ℝ) ^ 1) * ((S03.Z [-3] (2 * m) : ℚ) : ℝ) +
      (-1 / 3 : ℝ) * (1 / (2 * (m : ℝ) + 1) ^ 1) * ((S03.Z [-3] (2 * m) : ℚ) : ℝ) := by
  have q1 : (0 : ℝ) ≤ m := by positivity
  simp only [Tji, V05.rho, Finset.prod_range_succ, Finset.prod_range_zero]
  norm_num
  simp only [e303_eq, e312_eq, e321_eq, e330_eq, V05.coeff3, V05.Dsharp]
  norm_num
  try rw [show m + 1 + m = 2 * m + 1 by ring]
  try rw [show m + 1 + 1 + m = 2 * m + 1 + 1 by ring]
  try rw [show m + 1 + 2 + m = 2 * m + 1 + 1 + 1 by ring]
  try rw [show m + 1 + 3 + m = 2 * m + 1 + 1 + 1 + 1 by ring]
  try rw [show m + 2 = m + 1 + 1 by ring]
  try rw [show m + 3 = m + 1 + 1 + 1 by ring]
  simp only [V05.Hr_succ]
  simp only [Hr1m m, Hr2m m, Hr3m m, Hr12m m, Hr22m m, Hr32m m]
  push_cast
  have h1 : (2 * ((m : ℝ) + 1) - 1) ≠ 0 := by linarith
  have h2 : (2 * (m : ℝ) + 1) ≠ 0 := by positivity
  field_simp
  ring

set_option maxRecDepth 20000 in
lemma T_0_1 (m : ℕ) : Tji 0 1 m =
      (-3 / 4 : ℝ) * (1 / (2 * (m : ℝ) + 1) ^ 2) +
      (3 / 4 : ℝ) * (1 / (2 * (m : ℝ) + 1) ^ 4) +
      (-1 / 8 : ℝ) * ((S03.Z [2] (2 * m) : ℚ) : ℝ) +
      (1 / 8 : ℝ) * (1 / (2 * (m : ℝ) + 1) ^ 2) * ((S03.Z [2] (2 * m) : ℚ) : ℝ) +
      (1 / 8 : ℝ) * ((S03.Z [2] (2 * m) : ℚ) : ℝ) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) +
      (1 / 4 : ℝ) * ((m : ℝ) ^ 1) * ((S03.Z [2] (2 * m) : ℚ) : ℝ) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) +
      (-1 / 8 : ℝ) * (1 / (2 * (m : ℝ) + 1) ^ 1) * ((S03.Z [2] (2 * m) : ℚ) : ℝ) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) +
      (-1 / 4 : ℝ) * ((S03.Z [3] (2 * m) : ℚ) : ℝ) +
      (-1 / 2 : ℝ) * ((m : ℝ) ^ 1) * ((S03.Z [3] (2 * m) : ℚ) : ℝ) +
      (1 / 4 : ℝ) * (1 / (2 * (m : ℝ) + 1) ^ 1) * ((S03.Z [3] (2 * m) : ℚ) : ℝ) +
      (3 / 4 : ℝ) * (1 / (2 * (m : ℝ) + 1) ^ 1) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) +
      (-3 / 4 : ℝ) * (1 / (2 * (m : ℝ) + 1) ^ 3) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) +
      (-3 / 8 : ℝ) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) +
      (3 / 8 : ℝ) * (1 / (2 * (m : ℝ) + 1) ^ 2) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) +
      (1 / 8 : ℝ) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) +
      (1 / 4 : ℝ) * ((m : ℝ) ^ 1) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) +
      (-1 / 8 : ℝ) * (1 / (2 * (m : ℝ) + 1) ^ 1) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) +
      (-1 / 4 : ℝ) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) * ((S03.Z [-2] (2 * m) : ℚ) : ℝ) +
      (-1 / 2 : ℝ) * ((m : ℝ) ^ 1) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) * ((S03.Z [-2] (2 * m) : ℚ) : ℝ) +
      (1 / 4 : ℝ) * (1 / (2 * (m : ℝ) + 1) ^ 1) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) * ((S03.Z [-2] (2 * m) : ℚ) : ℝ) +
      (1 / 4 : ℝ) * ((S03.Z [-2] (2 * m) : ℚ) : ℝ) +
      (-1 / 4 : ℝ) * (1 / (2 * (m : ℝ) + 1) ^ 2) * ((S03.Z [-2] (2 * m) : ℚ) : ℝ) := by
  have q1 : (0 : ℝ) ≤ m := by positivity
  simp only [Tji, V05.rho, Finset.prod_range_succ, Finset.prod_range_zero]
  norm_num
  simp only [e303_eq, e312_eq, e321_eq, e330_eq, V05.coeff3, V05.Dsharp]
  norm_num
  try rw [show m + 1 + m = 2 * m + 1 by ring]
  try rw [show m + 1 + 1 + m = 2 * m + 1 + 1 by ring]
  try rw [show m + 1 + 2 + m = 2 * m + 1 + 1 + 1 by ring]
  try rw [show m + 1 + 3 + m = 2 * m + 1 + 1 + 1 + 1 by ring]
  try rw [show m + 2 = m + 1 + 1 by ring]
  try rw [show m + 3 = m + 1 + 1 + 1 by ring]
  simp only [V05.Hr_succ]
  simp only [Hr1m m, Hr2m m, Hr3m m, Hr12m m, Hr22m m, Hr32m m]
  push_cast
  have h1 : (2 * ((m : ℝ) + 1) - 1) ≠ 0 := by linarith
  have h2 : (2 * (m : ℝ) + 1) ≠ 0 := by positivity
  field_simp
  ring

set_option maxRecDepth 20000 in
lemma T_0_2 (m : ℕ) : Tji 0 2 m =
      (3 / 4 : ℝ) * (1 / (2 * (m : ℝ) + 1) ^ 2) +
      (-3 / 2 : ℝ) * (1 / (2 * (m : ℝ) + 1) ^ 3) +
      (3 / 4 : ℝ) * (1 / (2 * (m : ℝ) + 1) ^ 4) +
      (1 / 8 : ℝ) * ((S03.Z [2] (2 * m) : ℚ) : ℝ) +
      (-1 / 4 : ℝ) * (1 / (2 * (m : ℝ) + 1) ^ 1) * ((S03.Z [2] (2 * m) : ℚ) : ℝ) +
      (1 / 8 : ℝ) * (1 / (2 * (m : ℝ) + 1) ^ 2) * ((S03.Z [2] (2 * m) : ℚ) : ℝ) +
      (1 / 8 : ℝ) * ((S03.Z [2] (2 * m) : ℚ) : ℝ) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) +
      (-1 / 4 : ℝ) * ((m : ℝ) ^ 1) * ((S03.Z [2] (2 * m) : ℚ) : ℝ) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) +
      (-1 / 8 : ℝ) * (1 / (2 * (m : ℝ) + 1) ^ 1) * ((S03.Z [2] (2 * m) : ℚ) : ℝ) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) +
      (-1 / 4 : ℝ) * ((S03.Z [3] (2 * m) : ℚ) : ℝ) +
      (1 / 2 : ℝ) * ((m : ℝ) ^ 1) * ((S03.Z [3] (2 * m) : ℚ) : ℝ) +
      (1 / 4 : ℝ) * (1 / (2 * (m : ℝ) + 1) ^ 1) * ((S03.Z [3] (2 * m) : ℚ) : ℝ) +
      (-3 / 4 : ℝ) * (1 / (2 * (m : ℝ) + 1) ^ 1) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) +
      (3 / 2 : ℝ) * (1 / (2 * (m : ℝ) + 1) ^ 2) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) +
      (-3 / 4 : ℝ) * (1 / (2 * (m : ℝ) + 1) ^ 3) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) +
      (3 / 8 : ℝ) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) +
      (-3 / 4 : ℝ) * (1 / (2 * (m : ℝ) + 1) ^ 1) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) +
      (3 / 8 : ℝ) * (1 / (2 * (m : ℝ) + 1) ^ 2) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) +
      (1 / 8 : ℝ) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) +
      (-1 / 4 : ℝ) * ((m : ℝ) ^ 1) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) +
      (-1 / 8 : ℝ) * (1 / (2 * (m : ℝ) + 1) ^ 1) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) +
      (-1 / 4 : ℝ) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) * ((S03.Z [-2] (2 * m) : ℚ) : ℝ) +
      (1 / 2 : ℝ) * ((m : ℝ) ^ 1) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) * ((S03.Z [-2] (2 * m) : ℚ) : ℝ) +
      (1 / 4 : ℝ) * (1 / (2 * (m : ℝ) + 1) ^ 1) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) * ((S03.Z [-2] (2 * m) : ℚ) : ℝ) +
      (-1 / 4 : ℝ) * ((S03.Z [-2] (2 * m) : ℚ) : ℝ) +
      (1 / 2 : ℝ) * (1 / (2 * (m : ℝ) + 1) ^ 1) * ((S03.Z [-2] (2 * m) : ℚ) : ℝ) +
      (-1 / 4 : ℝ) * (1 / (2 * (m : ℝ) + 1) ^ 2) * ((S03.Z [-2] (2 * m) : ℚ) : ℝ) := by
  have q1 : (0 : ℝ) ≤ m := by positivity
  simp only [Tji, V05.rho, Finset.prod_range_succ, Finset.prod_range_zero]
  norm_num
  simp only [e303_eq, e312_eq, e321_eq, e330_eq, V05.coeff3, V05.Dsharp]
  norm_num
  try rw [show m + 1 + m = 2 * m + 1 by ring]
  try rw [show m + 1 + 1 + m = 2 * m + 1 + 1 by ring]
  try rw [show m + 1 + 2 + m = 2 * m + 1 + 1 + 1 by ring]
  try rw [show m + 1 + 3 + m = 2 * m + 1 + 1 + 1 + 1 by ring]
  try rw [show m + 2 = m + 1 + 1 by ring]
  try rw [show m + 3 = m + 1 + 1 + 1 by ring]
  simp only [V05.Hr_succ]
  simp only [Hr1m m, Hr2m m, Hr3m m, Hr12m m, Hr22m m, Hr32m m]
  push_cast
  have h1 : (2 * ((m : ℝ) + 1) - 1) ≠ 0 := by linarith
  have h2 : (2 * (m : ℝ) + 1) ≠ 0 := by positivity
  field_simp
  ring

set_option maxRecDepth 20000 in
lemma T_0_3 (m : ℕ) : Tji 0 3 m =
      (1 : ℝ) * (1 / ((m : ℝ) + 1) ^ 1) +
      (-2 : ℝ) * (1 / (2 * (m : ℝ) + 1) ^ 1) +
      (7 / 4 : ℝ) * (1 / (2 * (m : ℝ) + 1) ^ 2) +
      (-1 : ℝ) * (1 / (2 * (m : ℝ) + 1) ^ 3) +
      (1 / 4 : ℝ) * (1 / (2 * (m : ℝ) + 1) ^ 4) +
      (1 / 8 : ℝ) * ((S03.Z [2] (2 * m) : ℚ) : ℝ) +
      (-1 / 2 : ℝ) * (1 / ((m : ℝ) + 1) ^ 1) * ((S03.Z [2] (2 * m) : ℚ) : ℝ) +
      (1 / 2 : ℝ) * (1 / (2 * (m : ℝ) + 1) ^ 1) * ((S03.Z [2] (2 * m) : ℚ) : ℝ) +
      (-1 / 8 : ℝ) * (1 / (2 * (m : ℝ) + 1) ^ 2) * ((S03.Z [2] (2 * m) : ℚ) : ℝ) +
      (3 / 8 : ℝ) * ((S03.Z [2] (2 * m) : ℚ) : ℝ) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) +
      (-1 / 4 : ℝ) * ((m : ℝ) ^ 1) * ((S03.Z [2] (2 * m) : ℚ) : ℝ) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) +
      (-1 / 2 : ℝ) * (1 / ((m : ℝ) + 1) ^ 1) * ((S03.Z [2] (2 * m) : ℚ) : ℝ) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) +
      (1 / 8 : ℝ) * (1 / (2 * (m : ℝ) + 1) ^ 1) * ((S03.Z [2] (2 * m) : ℚ) : ℝ) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) +
      (-3 / 4 : ℝ) * ((S03.Z [3] (2 * m) : ℚ) : ℝ) +
      (1 / 2 : ℝ) * ((m : ℝ) ^ 1) * ((S03.Z [3] (2 * m) : ℚ) : ℝ) +
      (1 : ℝ) * (1 / ((m : ℝ) + 1) ^ 1) * ((S03.Z [3] (2 * m) : ℚ) : ℝ) +
      (-1 / 4 : ℝ) * (1 / (2 * (m : ℝ) + 1) ^ 1) * ((S03.Z [3] (2 * m) : ℚ) : ℝ) +
      (1 : ℝ) * (1 / ((m : ℝ) + 1) ^ 1) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) +
      (-7 / 4 : ℝ) * (1 / (2 * (m : ℝ) + 1) ^ 1) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) +
      (1 : ℝ) * (1 / (2 * (m : ℝ) + 1) ^ 2) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) +
      (-1 / 4 : ℝ) * (1 / (2 * (m : ℝ) + 1) ^ 3) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) +
      (-1 / 8 : ℝ) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) +
      (1 / 2 : ℝ) * (1 / ((m : ℝ) + 1) ^ 1) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) +
      (-1 / 2 : ℝ) * (1 / (2 * (m : ℝ) + 1) ^ 1) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) +
      (1 / 8 : ℝ) * (1 / (2 * (m : ℝ) + 1) ^ 2) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) +
      (-1 / 8 : ℝ) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) +
      (1 / 12 : ℝ) * ((m : ℝ) ^ 1) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) +
      (1 / 6 : ℝ) * (1 / ((m : ℝ) + 1) ^ 1) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) +
      (-1 / 24 : ℝ) * (1 / (2 * (m : ℝ) + 1) ^ 1) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) +
      (3 / 4 : ℝ) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) * ((S03.Z [-2] (2 * m) : ℚ) : ℝ) +
      (-1 / 2 : ℝ) * ((m : ℝ) ^ 1) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) * ((S03.Z [-2] (2 * m) : ℚ) : ℝ) +
      (-1 : ℝ) * (1 / ((m : ℝ) + 1) ^ 1) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) * ((S03.Z [-2] (2 * m) : ℚ) : ℝ) +
      (1 / 4 : ℝ) * (1 / (2 * (m : ℝ) + 1) ^ 1) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) * ((S03.Z [-2] (2 * m) : ℚ) : ℝ) +
      (1 / 4 : ℝ) * ((S03.Z [-2] (2 * m) : ℚ) : ℝ) +
      (-1 : ℝ) * (1 / ((m : ℝ) + 1) ^ 1) * ((S03.Z [-2] (2 * m) : ℚ) : ℝ) +
      (1 : ℝ) * (1 / (2 * (m : ℝ) + 1) ^ 1) * ((S03.Z [-2] (2 * m) : ℚ) : ℝ) +
      (-1 / 4 : ℝ) * (1 / (2 * (m : ℝ) + 1) ^ 2) * ((S03.Z [-2] (2 * m) : ℚ) : ℝ) +
      (-1 : ℝ) * ((S03.Z [-3] (2 * m) : ℚ) : ℝ) +
      (2 / 3 : ℝ) * ((m : ℝ) ^ 1) * ((S03.Z [-3] (2 * m) : ℚ) : ℝ) +
      (4 / 3 : ℝ) * (1 / ((m : ℝ) + 1) ^ 1) * ((S03.Z [-3] (2 * m) : ℚ) : ℝ) +
      (-1 / 3 : ℝ) * (1 / (2 * (m : ℝ) + 1) ^ 1) * ((S03.Z [-3] (2 * m) : ℚ) : ℝ) := by
  have q1 : (0 : ℝ) ≤ m := by positivity
  simp only [Tji, V05.rho, Finset.prod_range_succ, Finset.prod_range_zero]
  norm_num
  simp only [e303_eq, e312_eq, e321_eq, e330_eq, V05.coeff3, V05.Dsharp]
  norm_num
  try rw [show m + 1 + m = 2 * m + 1 by ring]
  try rw [show m + 1 + 1 + m = 2 * m + 1 + 1 by ring]
  try rw [show m + 1 + 2 + m = 2 * m + 1 + 1 + 1 by ring]
  try rw [show m + 1 + 3 + m = 2 * m + 1 + 1 + 1 + 1 by ring]
  try rw [show m + 2 = m + 1 + 1 by ring]
  try rw [show m + 3 = m + 1 + 1 + 1 by ring]
  simp only [V05.Hr_succ]
  simp only [Hr1m m, Hr2m m, Hr3m m, Hr12m m, Hr22m m, Hr32m m]
  push_cast
  have h1 : (2 * ((m : ℝ) + 1) - 1) ≠ 0 := by linarith
  have h2 : (2 * (m : ℝ) + 1) ≠ 0 := by positivity
  field_simp
  ring

set_option maxRecDepth 20000 in
lemma T_1_0 (m : ℕ) : Tji 1 0 m =
      (3 / 8 : ℝ) * (1 / ((m : ℝ) + 1) ^ 1) +
      (3 / 16 : ℝ) * (1 / ((m : ℝ) + 1) ^ 2) +
      (-1 / 16 : ℝ) * (1 / ((m : ℝ) + 1) ^ 3) +
      (-1 / 16 : ℝ) * (1 / ((m : ℝ) + 1) ^ 4) +
      (-3 / 4 : ℝ) * (1 / (2 * (m : ℝ) + 1) ^ 1) +
      (-3 : ℝ) * (1 / (2 * (m : ℝ) + 1) ^ 2) +
      (-13 / 4 : ℝ) * (1 / (2 * (m : ℝ) + 1) ^ 3) +
      (3 / 4 : ℝ) * ((S03.Z [2] (2 * m) : ℚ) : ℝ) +
      (1 / 8 : ℝ) * (1 / ((m : ℝ) + 1) ^ 1) * ((S03.Z [2] (2 * m) : ℚ) : ℝ) +
      (1 / 8 : ℝ) * (1 / ((m : ℝ) + 1) ^ 2) * ((S03.Z [2] (2 * m) : ℚ) : ℝ) +
      (13 / 8 : ℝ) * (1 / (2 * (m : ℝ) + 1) ^ 1) * ((S03.Z [2] (2 * m) : ℚ) : ℝ) +
      (-3 / 2 : ℝ) * ((S03.Z [2] (2 * m) : ℚ) : ℝ) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) +
      (-3 / 4 : ℝ) * ((m : ℝ) ^ 1) * ((S03.Z [2] (2 * m) : ℚ) : ℝ) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) +
      (-1 / 4 : ℝ) * (1 / ((m : ℝ) + 1) ^ 1) * ((S03.Z [2] (2 * m) : ℚ) : ℝ) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) +
      (3 : ℝ) * ((S03.Z [3] (2 * m) : ℚ) : ℝ) +
      (3 / 2 : ℝ) * ((m : ℝ) ^ 1) * ((S03.Z [3] (2 * m) : ℚ) : ℝ) +
      (1 / 2 : ℝ) * (1 / ((m : ℝ) + 1) ^ 1) * ((S03.Z [3] (2 * m) : ℚ) : ℝ) +
      (-3 / 8 : ℝ) * (1 / ((m : ℝ) + 1) ^ 1) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) +
      (1 / 8 : ℝ) * (1 / ((m : ℝ) + 1) ^ 2) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) +
      (1 / 8 : ℝ) * (1 / ((m : ℝ) + 1) ^ 3) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) +
      (3 : ℝ) * (1 / (2 * (m : ℝ) + 1) ^ 1) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) +
      (13 / 4 : ℝ) * (1 / (2 * (m : ℝ) + 1) ^ 2) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) +
      (-3 / 4 : ℝ) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) +
      (-1 / 8 : ℝ) * (1 / ((m : ℝ) + 1) ^ 1) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) +
      (-1 / 8 : ℝ) * (1 / ((m : ℝ) + 1) ^ 2) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) +
      (-13 / 8 : ℝ) * (1 / (2 * (m : ℝ) + 1) ^ 1) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) +
      (1 / 2 : ℝ) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) +
      (1 / 4 : ℝ) * ((m : ℝ) ^ 1) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) +
      (1 / 12 : ℝ) * (1 / ((m : ℝ) + 1) ^ 1) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) +
      (-3 : ℝ) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) * ((S03.Z [-2] (2 * m) : ℚ) : ℝ) +
      (-3 / 2 : ℝ) * ((m : ℝ) ^ 1) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) * ((S03.Z [-2] (2 * m) : ℚ) : ℝ) +
      (-1 / 2 : ℝ) * (1 / ((m : ℝ) + 1) ^ 1) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) * ((S03.Z [-2] (2 * m) : ℚ) : ℝ) +
      (3 / 2 : ℝ) * ((S03.Z [-2] (2 * m) : ℚ) : ℝ) +
      (1 / 4 : ℝ) * (1 / ((m : ℝ) + 1) ^ 1) * ((S03.Z [-2] (2 * m) : ℚ) : ℝ) +
      (1 / 4 : ℝ) * (1 / ((m : ℝ) + 1) ^ 2) * ((S03.Z [-2] (2 * m) : ℚ) : ℝ) +
      (13 / 4 : ℝ) * (1 / (2 * (m : ℝ) + 1) ^ 1) * ((S03.Z [-2] (2 * m) : ℚ) : ℝ) +
      (4 : ℝ) * ((S03.Z [-3] (2 * m) : ℚ) : ℝ) +
      (2 : ℝ) * ((m : ℝ) ^ 1) * ((S03.Z [-3] (2 * m) : ℚ) : ℝ) +
      (2 / 3 : ℝ) * (1 / ((m : ℝ) + 1) ^ 1) * ((S03.Z [-3] (2 * m) : ℚ) : ℝ) := by
  have q1 : (0 : ℝ) ≤ m := by positivity
  simp only [Tji, V05.rho, Finset.prod_range_succ, Finset.prod_range_zero]
  norm_num
  simp only [e303_eq, e312_eq, e321_eq, e330_eq, V05.coeff3, V05.Dsharp]
  norm_num
  try rw [show m + 1 + m = 2 * m + 1 by ring]
  try rw [show m + 1 + 1 + m = 2 * m + 1 + 1 by ring]
  try rw [show m + 1 + 2 + m = 2 * m + 1 + 1 + 1 by ring]
  try rw [show m + 1 + 3 + m = 2 * m + 1 + 1 + 1 + 1 by ring]
  try rw [show m + 2 = m + 1 + 1 by ring]
  try rw [show m + 3 = m + 1 + 1 + 1 by ring]
  simp only [V05.Hr_succ]
  simp only [Hr1m m, Hr2m m, Hr3m m, Hr12m m, Hr22m m, Hr32m m]
  push_cast
  have h1 : (2 * ((m : ℝ) + 1) - 1) ≠ 0 := by linarith
  have h2 : (2 * (m : ℝ) + 1) ≠ 0 := by positivity
  field_simp
  ring

set_option maxRecDepth 20000 in
lemma T_1_1 (m : ℕ) : Tji 1 1 m =
      (-3 / 8 : ℝ) * (1 / ((m : ℝ) + 1) ^ 1) +
      (1 / 16 : ℝ) * (1 / ((m : ℝ) + 1) ^ 2) +
      (1 / 16 : ℝ) * (1 / ((m : ℝ) + 1) ^ 3) +
      (3 / 4 : ℝ) * (1 / (2 * (m : ℝ) + 1) ^ 1) +
      (7 / 2 : ℝ) * (1 / (2 * (m : ℝ) + 1) ^ 2) +
      (15 / 4 : ℝ) * (1 / (2 * (m : ℝ) + 1) ^ 3) +
      (3 / 2 : ℝ) * ((S03.Z [2] (2 * m) : ℚ) : ℝ) +
      (3 / 8 : ℝ) * (1 / ((m : ℝ) + 1) ^ 1) * ((S03.Z [2] (2 * m) : ℚ) : ℝ) +
      (5 / 8 : ℝ) * (1 / (2 * (m : ℝ) + 1) ^ 1) * ((S03.Z [2] (2 * m) : ℚ) : ℝ) +
      (-1 : ℝ) * ((S03.Z [2] (2 * m) : ℚ) : ℝ) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) +
      (-3 / 4 : ℝ) * ((m : ℝ) ^ 1) * ((S03.Z [2] (2 * m) : ℚ) : ℝ) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) +
      (2 : ℝ) * ((S03.Z [3] (2 * m) : ℚ) : ℝ) +
      (3 / 2 : ℝ) * ((m : ℝ) ^ 1) * ((S03.Z [3] (2 * m) : ℚ) : ℝ) +
      (-1 / 8 : ℝ) * (1 / ((m : ℝ) + 1) ^ 1) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) +
      (-1 / 8 : ℝ) * (1 / ((m : ℝ) + 1) ^ 2) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) +
      (-7 / 2 : ℝ) * (1 / (2 * (m : ℝ) + 1) ^ 1) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) +
      (-15 / 4 : ℝ) * (1 / (2 * (m : ℝ) + 1) ^ 2) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) +
      (3 / 2 : ℝ) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) +
      (1 / 8 : ℝ) * (1 / ((m : ℝ) + 1) ^ 1) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) +
      (15 / 8 : ℝ) * (1 / (2 * (m : ℝ) + 1) ^ 1) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) +
      (-1 : ℝ) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) +
      (-3 / 4 : ℝ) * ((m : ℝ) ^ 1) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) +
      (2 : ℝ) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) * ((S03.Z [-2] (2 * m) : ℚ) : ℝ) +
      (3 / 2 : ℝ) * ((m : ℝ) ^ 1) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) * ((S03.Z [-2] (2 * m) : ℚ) : ℝ) +
      (1 / 4 : ℝ) * (1 / ((m : ℝ) + 1) ^ 1) * ((S03.Z [-2] (2 * m) : ℚ) : ℝ) +
      (-5 / 4 : ℝ) * (1 / (2 * (m : ℝ) + 1) ^ 1) * ((S03.Z [-2] (2 * m) : ℚ) : ℝ) := by
  have q1 : (0 : ℝ) ≤ m := by positivity
  simp only [Tji, V05.rho, Finset.prod_range_succ, Finset.prod_range_zero]
  norm_num
  simp only [e303_eq, e312_eq, e321_eq, e330_eq, V05.coeff3, V05.Dsharp]
  norm_num
  try rw [show m + 1 + m = 2 * m + 1 by ring]
  try rw [show m + 1 + 1 + m = 2 * m + 1 + 1 by ring]
  try rw [show m + 1 + 2 + m = 2 * m + 1 + 1 + 1 by ring]
  try rw [show m + 1 + 3 + m = 2 * m + 1 + 1 + 1 + 1 by ring]
  try rw [show m + 2 = m + 1 + 1 by ring]
  try rw [show m + 3 = m + 1 + 1 + 1 by ring]
  simp only [V05.Hr_succ]
  simp only [Hr1m m, Hr2m m, Hr3m m, Hr12m m, Hr22m m, Hr32m m]
  push_cast
  have h1 : (2 * ((m : ℝ) + 1) - 1) ≠ 0 := by linarith
  have h2 : (2 * (m : ℝ) + 1) ≠ 0 := by positivity
  field_simp
  ring

set_option maxRecDepth 20000 in
lemma T_1_2 (m : ℕ) : Tji 1 2 m =
      (-3 / 8 : ℝ) * (1 / ((m : ℝ) + 1) ^ 1) +
      (5 / 16 : ℝ) * (1 / ((m : ℝ) + 1) ^ 2) +
      (-1 / 16 : ℝ) * (1 / ((m : ℝ) + 1) ^ 3) +
      (3 / 4 : ℝ) * (1 / (2 * (m : ℝ) + 1) ^ 1) +
      (-2 : ℝ) * (1 / (2 * (m : ℝ) + 1) ^ 2) +
      (-3 / 4 : ℝ) * (1 / (2 * (m : ℝ) + 1) ^ 3) +
      (3 / 4 : ℝ) * ((S03.Z [2] (2 * m) : ℚ) : ℝ) +
      (-3 / 8 : ℝ) * (1 / ((m : ℝ) + 1) ^ 1) * ((S03.Z [2] (2 * m) : ℚ) : ℝ) +
      (-1 / 8 : ℝ) * (1 / (2 * (m : ℝ) + 1) ^ 1) * ((S03.Z [2] (2 * m) : ℚ) : ℝ) +
      (1 / 2 : ℝ) * ((S03.Z [2] (2 * m) : ℚ) : ℝ) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) +
      (3 / 4 : ℝ) * ((m : ℝ) ^ 1) * ((S03.Z [2] (2 * m) : ℚ) : ℝ) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) +
      (-1 : ℝ) * ((S03.Z [3] (2 * m) : ℚ) : ℝ) +
      (-3 / 2 : ℝ) * ((m : ℝ) ^ 1) * ((S03.Z [3] (2 * m) : ℚ) : ℝ) +
      (-5 / 8 : ℝ) * (1 / ((m : ℝ) + 1) ^ 1) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) +
      (1 / 8 : ℝ) * (1 / ((m : ℝ) + 1) ^ 2) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) +
      (2 : ℝ) * (1 / (2 * (m : ℝ) + 1) ^ 1) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) +
      (3 / 4 : ℝ) * (1 / (2 * (m : ℝ) + 1) ^ 2) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) +
      (-3 / 4 : ℝ) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) +
      (-1 / 8 : ℝ) * (1 / ((m : ℝ) + 1) ^ 1) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) +
      (-3 / 8 : ℝ) * (1 / (2 * (m : ℝ) + 1) ^ 1) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) +
      (1 / 2 : ℝ) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) +
      (3 / 4 : ℝ) * ((m : ℝ) ^ 1) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) +
      (-1 : ℝ) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) * ((S03.Z [-2] (2 * m) : ℚ) : ℝ) +
      (-3 / 2 : ℝ) * ((m : ℝ) ^ 1) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) * ((S03.Z [-2] (2 * m) : ℚ) : ℝ) +
      (3 / 2 : ℝ) * ((S03.Z [-2] (2 * m) : ℚ) : ℝ) +
      (-1 / 4 : ℝ) * (1 / ((m : ℝ) + 1) ^ 1) * ((S03.Z [-2] (2 * m) : ℚ) : ℝ) +
      (1 / 4 : ℝ) * (1 / (2 * (m : ℝ) + 1) ^ 1) * ((S03.Z [-2] (2 * m) : ℚ) : ℝ) := by
  have q1 : (0 : ℝ) ≤ m := by positivity
  simp only [Tji, V05.rho, Finset.prod_range_succ, Finset.prod_range_zero]
  norm_num
  simp only [e303_eq, e312_eq, e321_eq, e330_eq, V05.coeff3, V05.Dsharp]
  norm_num
  try rw [show m + 1 + m = 2 * m + 1 by ring]
  try rw [show m + 1 + 1 + m = 2 * m + 1 + 1 by ring]
  try rw [show m + 1 + 2 + m = 2 * m + 1 + 1 + 1 by ring]
  try rw [show m + 1 + 3 + m = 2 * m + 1 + 1 + 1 + 1 by ring]
  try rw [show m + 2 = m + 1 + 1 by ring]
  try rw [show m + 3 = m + 1 + 1 + 1 by ring]
  simp only [V05.Hr_succ]
  simp only [Hr1m m, Hr2m m, Hr3m m, Hr12m m, Hr22m m, Hr32m m]
  push_cast
  have h1 : (2 * ((m : ℝ) + 1) - 1) ≠ 0 := by linarith
  have h2 : (2 * (m : ℝ) + 1) ≠ 0 := by positivity
  field_simp
  ring

set_option maxRecDepth 20000 in
lemma T_1_3 (m : ℕ) : Tji 1 3 m =
      (-5 / 8 : ℝ) * (1 / ((m : ℝ) + 1) ^ 1) +
      (-9 / 16 : ℝ) * (1 / ((m : ℝ) + 1) ^ 2) +
      (5 / 16 : ℝ) * (1 / ((m : ℝ) + 1) ^ 3) +
      (-1 / 16 : ℝ) * (1 / ((m : ℝ) + 1) ^ 4) +
      (5 / 4 : ℝ) * (1 / (2 * (m : ℝ) + 1) ^ 1) +
      (-1 / 2 : ℝ) * (1 / (2 * (m : ℝ) + 1) ^ 2) +
      (1 / 4 : ℝ) * (1 / (2 * (m : ℝ) + 1) ^ 3) +
      (-1 / 8 : ℝ) * (1 / ((m : ℝ) + 1) ^ 1) * ((S03.Z [2] (2 * m) : ℚ) : ℝ) +
      (1 / 8 : ℝ) * (1 / ((m : ℝ) + 1) ^ 2) * ((S03.Z [2] (2 * m) : ℚ) : ℝ) +
      (-1 / 8 : ℝ) * (1 / (2 * (m : ℝ) + 1) ^ 1) * ((S03.Z [2] (2 * m) : ℚ) : ℝ) +
      (3 / 4 : ℝ) * ((m : ℝ) ^ 1) * ((S03.Z [2] (2 * m) : ℚ) : ℝ) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) +
      (1 / 4 : ℝ) * (1 / ((m : ℝ) + 1) ^ 1) * ((S03.Z [2] (2 * m) : ℚ) : ℝ) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) +
      (-3 / 2 : ℝ) * ((m : ℝ) ^ 1) * ((S03.Z [3] (2 * m) : ℚ) : ℝ) +
      (-1 / 2 : ℝ) * (1 / ((m : ℝ) + 1) ^ 1) * ((S03.Z [3] (2 * m) : ℚ) : ℝ) +
      (1 / 8 : ℝ) * (1 / ((m : ℝ) + 1) ^ 1) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) +
      (-5 / 8 : ℝ) * (1 / ((m : ℝ) + 1) ^ 2) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) +
      (1 / 8 : ℝ) * (1 / ((m : ℝ) + 1) ^ 3) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) +
      (1 / 2 : ℝ) * (1 / (2 * (m : ℝ) + 1) ^ 1) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) +
      (-1 / 4 : ℝ) * (1 / (2 * (m : ℝ) + 1) ^ 2) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) +
      (1 / 8 : ℝ) * (1 / ((m : ℝ) + 1) ^ 1) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) +
      (-1 / 8 : ℝ) * (1 / ((m : ℝ) + 1) ^ 2) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) +
      (1 / 8 : ℝ) * (1 / (2 * (m : ℝ) + 1) ^ 1) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) +
      (-1 / 4 : ℝ) * ((m : ℝ) ^ 1) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) +
      (-1 / 12 : ℝ) * (1 / ((m : ℝ) + 1) ^ 1) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) +
      (3 / 2 : ℝ) * ((m : ℝ) ^ 1) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) * ((S03.Z [-2] (2 * m) : ℚ) : ℝ) +
      (1 / 2 : ℝ) * (1 / ((m : ℝ) + 1) ^ 1) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) * ((S03.Z [-2] (2 * m) : ℚ) : ℝ) +
      (-1 / 4 : ℝ) * (1 / ((m : ℝ) + 1) ^ 1) * ((S03.Z [-2] (2 * m) : ℚ) : ℝ) +
      (1 / 4 : ℝ) * (1 / ((m : ℝ) + 1) ^ 2) * ((S03.Z [-2] (2 * m) : ℚ) : ℝ) +
      (-1 / 4 : ℝ) * (1 / (2 * (m : ℝ) + 1) ^ 1) * ((S03.Z [-2] (2 * m) : ℚ) : ℝ) +
      (-2 : ℝ) * ((m : ℝ) ^ 1) * ((S03.Z [-3] (2 * m) : ℚ) : ℝ) +
      (-2 / 3 : ℝ) * (1 / ((m : ℝ) + 1) ^ 1) * ((S03.Z [-3] (2 * m) : ℚ) : ℝ) := by
  have q1 : (0 : ℝ) ≤ m := by positivity
  simp only [Tji, V05.rho, Finset.prod_range_succ, Finset.prod_range_zero]
  norm_num
  simp only [e303_eq, e312_eq, e321_eq, e330_eq, V05.coeff3, V05.Dsharp]
  norm_num
  try rw [show m + 1 + m = 2 * m + 1 by ring]
  try rw [show m + 1 + 1 + m = 2 * m + 1 + 1 by ring]
  try rw [show m + 1 + 2 + m = 2 * m + 1 + 1 + 1 by ring]
  try rw [show m + 1 + 3 + m = 2 * m + 1 + 1 + 1 + 1 by ring]
  try rw [show m + 2 = m + 1 + 1 by ring]
  try rw [show m + 3 = m + 1 + 1 + 1 by ring]
  simp only [V05.Hr_succ]
  simp only [Hr1m m, Hr2m m, Hr3m m, Hr12m m, Hr22m m, Hr32m m]
  push_cast
  have h1 : (2 * ((m : ℝ) + 1) - 1) ≠ 0 := by linarith
  have h2 : (2 * (m : ℝ) + 1) ≠ 0 := by positivity
  field_simp
  ring

set_option maxRecDepth 20000 in
lemma T_2_0 (m : ℕ) : Tji 2 0 m =
      (-9 / 8 : ℝ) * (1 / ((m : ℝ) + 1) ^ 1) +
      (3 / 8 : ℝ) * (1 / ((m : ℝ) + 1) ^ 3) +
      (21 / 8 : ℝ) * (1 / (2 * (m : ℝ) + 1) ^ 1) +
      (45 / 8 : ℝ) * (1 / (2 * (m : ℝ) + 1) ^ 2) +
      (27 / 8 : ℝ) * (1 / (2 * (m : ℝ) + 1) ^ 3) +
      (-3 / 8 : ℝ) * (1 / (2 * (m : ℝ) + 3) ^ 1) +
      (-9 / 8 : ℝ) * (1 / (2 * (m : ℝ) + 3) ^ 2) +
      (3 / 8 : ℝ) * (1 / (2 * (m : ℝ) + 3) ^ 3) +
      (3 / 4 : ℝ) * (1 / (2 * (m : ℝ) + 3) ^ 4) +
      (-9 / 8 : ℝ) * ((S03.Z [2] (2 * m) : ℚ) : ℝ) +
      (-3 / 4 : ℝ) * (1 / ((m : ℝ) + 1) ^ 1) * ((S03.Z [2] (2 * m) : ℚ) : ℝ) +
      (-27 / 16 : ℝ) * (1 / (2 * (m : ℝ) + 1) ^ 1) * ((S03.Z [2] (2 * m) : ℚ) : ℝ) +
      (-3 / 16 : ℝ) * (1 / (2 * (m : ℝ) + 3) ^ 1) * ((S03.Z [2] (2 * m) : ℚ) : ℝ) +
      (-3 / 8 : ℝ) * (1 / (2 * (m : ℝ) + 3) ^ 2) * ((S03.Z [2] (2 * m) : ℚ) : ℝ) +
      (15 / 8 : ℝ) * ((S03.Z [2] (2 * m) : ℚ) : ℝ) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) +
      (3 / 4 : ℝ) * ((m : ℝ) ^ 1) * ((S03.Z [2] (2 * m) : ℚ) : ℝ) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) +
      (3 / 8 : ℝ) * (1 / (2 * (m : ℝ) + 3) ^ 1) * ((S03.Z [2] (2 * m) : ℚ) : ℝ) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) +
      (-15 / 4 : ℝ) * ((S03.Z [3] (2 * m) : ℚ) : ℝ) +
      (-3 / 2 : ℝ) * ((m : ℝ) ^ 1) * ((S03.Z [3] (2 * m) : ℚ) : ℝ) +
      (-3 / 4 : ℝ) * (1 / (2 * (m : ℝ) + 3) ^ 1) * ((S03.Z [3] (2 * m) : ℚ) : ℝ) +
      (-3 / 4 : ℝ) * (1 / ((m : ℝ) + 1) ^ 2) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) +
      (-45 / 8 : ℝ) * (1 / (2 * (m : ℝ) + 1) ^ 1) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) +
      (-27 / 8 : ℝ) * (1 / (2 * (m : ℝ) + 1) ^ 2) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) +
      (9 / 8 : ℝ) * (1 / (2 * (m : ℝ) + 3) ^ 1) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) +
      (-3 / 8 : ℝ) * (1 / (2 * (m : ℝ) + 3) ^ 2) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) +
      (-3 / 4 : ℝ) * (1 / (2 * (m : ℝ) + 3) ^ 3) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) +
      (9 / 8 : ℝ) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) +
      (3 / 4 : ℝ) * (1 / ((m : ℝ) + 1) ^ 1) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) +
      (27 / 16 : ℝ) * (1 / (2 * (m : ℝ) + 1) ^ 1) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) +
      (3 / 16 : ℝ) * (1 / (2 * (m : ℝ) + 3) ^ 1) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) +
      (3 / 8 : ℝ) * (1 / (2 * (m : ℝ) + 3) ^ 2) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) +
      (-5 / 8 : ℝ) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) +
      (-1 / 4 : ℝ) * ((m : ℝ) ^ 1) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) +
      (-1 / 8 : ℝ) * (1 / (2 * (m : ℝ) + 3) ^ 1) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) +
      (15 / 4 : ℝ) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) * ((S03.Z [-2] (2 * m) : ℚ) : ℝ) +
      (3 / 2 : ℝ) * ((m : ℝ) ^ 1) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) * ((S03.Z [-2] (2 * m) : ℚ) : ℝ) +
      (3 / 4 : ℝ) * (1 / (2 * (m : ℝ) + 3) ^ 1) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) * ((S03.Z [-2] (2 * m) : ℚ) : ℝ) +
      (-9 / 4 : ℝ) * ((S03.Z [-2] (2 * m) : ℚ) : ℝ) +
      (-3 / 2 : ℝ) * (1 / ((m : ℝ) + 1) ^ 1) * ((S03.Z [-2] (2 * m) : ℚ) : ℝ) +
      (-27 / 8 : ℝ) * (1 / (2 * (m : ℝ) + 1) ^ 1) * ((S03.Z [-2] (2 * m) : ℚ) : ℝ) +
      (-3 / 8 : ℝ) * (1 / (2 * (m : ℝ) + 3) ^ 1) * ((S03.Z [-2] (2 * m) : ℚ) : ℝ) +
      (-3 / 4 : ℝ) * (1 / (2 * (m : ℝ) + 3) ^ 2) * ((S03.Z [-2] (2 * m) : ℚ) : ℝ) +
      (-5 : ℝ) * ((S03.Z [-3] (2 * m) : ℚ) : ℝ) +
      (-2 : ℝ) * ((m : ℝ) ^ 1) * ((S03.Z [-3] (2 * m) : ℚ) : ℝ) +
      (-1 : ℝ) * (1 / (2 * (m : ℝ) + 3) ^ 1) * ((S03.Z [-3] (2 * m) : ℚ) : ℝ) := by
  have q1 : (0 : ℝ) ≤ m := by positivity
  simp only [Tji, V05.rho, Finset.prod_range_succ, Finset.prod_range_zero]
  norm_num
  simp only [e303_eq, e312_eq, e321_eq, e330_eq, V05.coeff3, V05.Dsharp]
  norm_num
  try rw [show m + 1 + m = 2 * m + 1 by ring]
  try rw [show m + 1 + 1 + m = 2 * m + 1 + 1 by ring]
  try rw [show m + 1 + 2 + m = 2 * m + 1 + 1 + 1 by ring]
  try rw [show m + 1 + 3 + m = 2 * m + 1 + 1 + 1 + 1 by ring]
  try rw [show m + 2 = m + 1 + 1 by ring]
  try rw [show m + 3 = m + 1 + 1 + 1 by ring]
  simp only [V05.Hr_succ]
  simp only [Hr1m m, Hr2m m, Hr3m m, Hr12m m, Hr22m m, Hr32m m]
  push_cast
  have h1 : (2 * ((m : ℝ) + 1) - 1) ≠ 0 := by linarith
  have h2 : (2 * (m : ℝ) + 1) ≠ 0 := by positivity
  field_simp
  ring

set_option maxRecDepth 20000 in
lemma T_2_1 (m : ℕ) : Tji 2 1 m =
      (5 / 8 : ℝ) * (1 / ((m : ℝ) + 1) ^ 1) +
      (3 / 8 : ℝ) * (1 / ((m : ℝ) + 1) ^ 2) +
      (-1 / 4 : ℝ) * (1 / ((m : ℝ) + 1) ^ 3) +
      (-5 / 8 : ℝ) * (1 / (2 * (m : ℝ) + 1) ^ 1) +
      (-55 / 8 : ℝ) * (1 / (2 * (m : ℝ) + 1) ^ 2) +
      (-63 / 8 : ℝ) * (1 / (2 * (m : ℝ) + 1) ^ 3) +
      (-5 / 8 : ℝ) * (1 / (2 * (m : ℝ) + 3) ^ 1) +
      (7 / 8 : ℝ) * (1 / (2 * (m : ℝ) + 3) ^ 2) +
      (-15 / 8 : ℝ) * (1 / (2 * (m : ℝ) + 3) ^ 3) +
      (-3 / 4 : ℝ) * (1 / (2 * (m : ℝ) + 3) ^ 4) +
      (-21 / 8 : ℝ) * ((S03.Z [2] (2 * m) : ℚ) : ℝ) +
      (-3 / 2 : ℝ) * (1 / ((m : ℝ) + 1) ^ 1) * ((S03.Z [2] (2 * m) : ℚ) : ℝ) +
      (-21 / 16 : ℝ) * (1 / (2 * (m : ℝ) + 1) ^ 1) * ((S03.Z [2] (2 * m) : ℚ) : ℝ) +
      (-5 / 16 : ℝ) * (1 / (2 * (m : ℝ) + 3) ^ 1) * ((S03.Z [2] (2 * m) : ℚ) : ℝ) +
      (-1 / 8 : ℝ) * (1 / (2 * (m : ℝ) + 3) ^ 2) * ((S03.Z [2] (2 * m) : ℚ) : ℝ) +
      (13 / 8 : ℝ) * ((S03.Z [2] (2 * m) : ℚ) : ℝ) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) +
      (3 / 4 : ℝ) * ((m : ℝ) ^ 1) * ((S03.Z [2] (2 * m) : ℚ) : ℝ) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) +
      (1 / 8 : ℝ) * (1 / (2 * (m : ℝ) + 3) ^ 1) * ((S03.Z [2] (2 * m) : ℚ) : ℝ) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) +
      (-13 / 4 : ℝ) * ((S03.Z [3] (2 * m) : ℚ) : ℝ) +
      (-3 / 2 : ℝ) * ((m : ℝ) ^ 1) * ((S03.Z [3] (2 * m) : ℚ) : ℝ) +
      (-1 / 4 : ℝ) * (1 / (2 * (m : ℝ) + 3) ^ 1) * ((S03.Z [3] (2 * m) : ℚ) : ℝ) +
      (-3 / 4 : ℝ) * (1 / ((m : ℝ) + 1) ^ 1) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) +
      (1 / 2 : ℝ) * (1 / ((m : ℝ) + 1) ^ 2) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) +
      (55 / 8 : ℝ) * (1 / (2 * (m : ℝ) + 1) ^ 1) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) +
      (63 / 8 : ℝ) * (1 / (2 * (m : ℝ) + 1) ^ 2) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) +
      (-7 / 8 : ℝ) * (1 / (2 * (m : ℝ) + 3) ^ 1) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) +
      (15 / 8 : ℝ) * (1 / (2 * (m : ℝ) + 3) ^ 2) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) +
      (3 / 4 : ℝ) * (1 / (2 * (m : ℝ) + 3) ^ 3) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) +
      (-15 / 8 : ℝ) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) +
      (-1 / 2 : ℝ) * (1 / ((m : ℝ) + 1) ^ 1) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) +
      (-63 / 16 : ℝ) * (1 / (2 * (m : ℝ) + 1) ^ 1) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) +
      (-15 / 16 : ℝ) * (1 / (2 * (m : ℝ) + 3) ^ 1) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) +
      (-3 / 8 : ℝ) * (1 / (2 * (m : ℝ) + 3) ^ 2) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) +
      (13 / 8 : ℝ) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) +
      (3 / 4 : ℝ) * ((m : ℝ) ^ 1) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) +
      (1 / 8 : ℝ) * (1 / (2 * (m : ℝ) + 3) ^ 1) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) +
      (-13 / 4 : ℝ) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) * ((S03.Z [-2] (2 * m) : ℚ) : ℝ) +
      (-3 / 2 : ℝ) * ((m : ℝ) ^ 1) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) * ((S03.Z [-2] (2 * m) : ℚ) : ℝ) +
      (-1 / 4 : ℝ) * (1 / (2 * (m : ℝ) + 3) ^ 1) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) * ((S03.Z [-2] (2 * m) : ℚ) : ℝ) +
      (-3 / 4 : ℝ) * ((S03.Z [-2] (2 * m) : ℚ) : ℝ) +
      (-1 : ℝ) * (1 / ((m : ℝ) + 1) ^ 1) * ((S03.Z [-2] (2 * m) : ℚ) : ℝ) +
      (21 / 8 : ℝ) * (1 / (2 * (m : ℝ) + 1) ^ 1) * ((S03.Z [-2] (2 * m) : ℚ) : ℝ) +
      (5 / 8 : ℝ) * (1 / (2 * (m : ℝ) + 3) ^ 1) * ((S03.Z [-2] (2 * m) : ℚ) : ℝ) +
      (1 / 4 : ℝ) * (1 / (2 * (m : ℝ) + 3) ^ 2) * ((S03.Z [-2] (2 * m) : ℚ) : ℝ) := by
  have q1 : (0 : ℝ) ≤ m := by positivity
  simp only [Tji, V05.rho, Finset.prod_range_succ, Finset.prod_range_zero]
  norm_num
  simp only [e303_eq, e312_eq, e321_eq, e330_eq, V05.coeff3, V05.Dsharp]
  norm_num
  try rw [show m + 1 + m = 2 * m + 1 by ring]
  try rw [show m + 1 + 1 + m = 2 * m + 1 + 1 by ring]
  try rw [show m + 1 + 2 + m = 2 * m + 1 + 1 + 1 by ring]
  try rw [show m + 1 + 3 + m = 2 * m + 1 + 1 + 1 + 1 by ring]
  try rw [show m + 2 = m + 1 + 1 by ring]
  try rw [show m + 3 = m + 1 + 1 + 1 by ring]
  simp only [V05.Hr_succ]
  simp only [Hr1m m, Hr2m m, Hr3m m, Hr12m m, Hr22m m, Hr32m m]
  push_cast
  have h1 : (2 * ((m : ℝ) + 1) - 1) ≠ 0 := by linarith
  have h2 : (2 * (m : ℝ) + 1) ≠ 0 := by positivity
  field_simp
  ring

set_option maxRecDepth 20000 in
lemma T_2_2 (m : ℕ) : Tji 2 2 m =
      (5 / 8 : ℝ) * (1 / ((m : ℝ) + 1) ^ 1) +
      (-1 / 4 : ℝ) * (1 / ((m : ℝ) + 1) ^ 2) +
      (-1 / 8 : ℝ) * (1 / ((m : ℝ) + 1) ^ 3) +
      (-17 / 8 : ℝ) * (1 / (2 * (m : ℝ) + 1) ^ 1) +
      (7 / 8 : ℝ) * (1 / (2 * (m : ℝ) + 1) ^ 2) +
      (45 / 8 : ℝ) * (1 / (2 * (m : ℝ) + 1) ^ 3) +
      (7 / 8 : ℝ) * (1 / (2 * (m : ℝ) + 3) ^ 1) +
      (-11 / 8 : ℝ) * (1 / (2 * (m : ℝ) + 3) ^ 2) +
      (21 / 8 : ℝ) * (1 / (2 * (m : ℝ) + 3) ^ 3) +
      (-3 / 4 : ℝ) * (1 / (2 * (m : ℝ) + 3) ^ 4) +
      (-15 / 8 : ℝ) * ((S03.Z [2] (2 * m) : ℚ) : ℝ) +
      (-3 / 4 : ℝ) * (1 / ((m : ℝ) + 1) ^ 1) * ((S03.Z [2] (2 * m) : ℚ) : ℝ) +
      (15 / 16 : ℝ) * (1 / (2 * (m : ℝ) + 1) ^ 1) * ((S03.Z [2] (2 * m) : ℚ) : ℝ) +
      (7 / 16 : ℝ) * (1 / (2 * (m : ℝ) + 3) ^ 1) * ((S03.Z [2] (2 * m) : ℚ) : ℝ) +
      (-1 / 8 : ℝ) * (1 / (2 * (m : ℝ) + 3) ^ 2) * ((S03.Z [2] (2 * m) : ℚ) : ℝ) +
      (-11 / 8 : ℝ) * ((S03.Z [2] (2 * m) : ℚ) : ℝ) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) +
      (-3 / 4 : ℝ) * ((m : ℝ) ^ 1) * ((S03.Z [2] (2 * m) : ℚ) : ℝ) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) +
      (1 / 8 : ℝ) * (1 / (2 * (m : ℝ) + 3) ^ 1) * ((S03.Z [2] (2 * m) : ℚ) : ℝ) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) +
      (11 / 4 : ℝ) * ((S03.Z [3] (2 * m) : ℚ) : ℝ) +
      (3 / 2 : ℝ) * ((m : ℝ) ^ 1) * ((S03.Z [3] (2 * m) : ℚ) : ℝ) +
      (-1 / 4 : ℝ) * (1 / (2 * (m : ℝ) + 3) ^ 1) * ((S03.Z [3] (2 * m) : ℚ) : ℝ) +
      (1 / 2 : ℝ) * (1 / ((m : ℝ) + 1) ^ 1) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) +
      (1 / 4 : ℝ) * (1 / ((m : ℝ) + 1) ^ 2) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) +
      (-7 / 8 : ℝ) * (1 / (2 * (m : ℝ) + 1) ^ 1) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) +
      (-45 / 8 : ℝ) * (1 / (2 * (m : ℝ) + 1) ^ 2) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) +
      (11 / 8 : ℝ) * (1 / (2 * (m : ℝ) + 3) ^ 1) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) +
      (-21 / 8 : ℝ) * (1 / (2 * (m : ℝ) + 3) ^ 2) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) +
      (3 / 4 : ℝ) * (1 / (2 * (m : ℝ) + 3) ^ 3) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) +
      (3 / 8 : ℝ) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) +
      (-1 / 4 : ℝ) * (1 / ((m : ℝ) + 1) ^ 1) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) +
      (45 / 16 : ℝ) * (1 / (2 * (m : ℝ) + 1) ^ 1) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) +
      (21 / 16 : ℝ) * (1 / (2 * (m : ℝ) + 3) ^ 1) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) +
      (-3 / 8 : ℝ) * (1 / (2 * (m : ℝ) + 3) ^ 2) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) +
      (-11 / 8 : ℝ) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) +
      (-3 / 4 : ℝ) * ((m : ℝ) ^ 1) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) +
      (1 / 8 : ℝ) * (1 / (2 * (m : ℝ) + 3) ^ 1) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) +
      (11 / 4 : ℝ) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) * ((S03.Z [-2] (2 * m) : ℚ) : ℝ) +
      (3 / 2 : ℝ) * ((m : ℝ) ^ 1) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) * ((S03.Z [-2] (2 * m) : ℚ) : ℝ) +
      (-1 / 4 : ℝ) * (1 / (2 * (m : ℝ) + 3) ^ 1) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) * ((S03.Z [-2] (2 * m) : ℚ) : ℝ) +
      (-9 / 4 : ℝ) * ((S03.Z [-2] (2 * m) : ℚ) : ℝ) +
      (-1 / 2 : ℝ) * (1 / ((m : ℝ) + 1) ^ 1) * ((S03.Z [-2] (2 * m) : ℚ) : ℝ) +
      (-15 / 8 : ℝ) * (1 / (2 * (m : ℝ) + 1) ^ 1) * ((S03.Z [-2] (2 * m) : ℚ) : ℝ) +
      (-7 / 8 : ℝ) * (1 / (2 * (m : ℝ) + 3) ^ 1) * ((S03.Z [-2] (2 * m) : ℚ) : ℝ) +
      (1 / 4 : ℝ) * (1 / (2 * (m : ℝ) + 3) ^ 2) * ((S03.Z [-2] (2 * m) : ℚ) : ℝ) := by
  have q1 : (0 : ℝ) ≤ m := by positivity
  simp only [Tji, V05.rho, Finset.prod_range_succ, Finset.prod_range_zero]
  norm_num
  simp only [e303_eq, e312_eq, e321_eq, e330_eq, V05.coeff3, V05.Dsharp]
  norm_num
  try rw [show m + 1 + m = 2 * m + 1 by ring]
  try rw [show m + 1 + 1 + m = 2 * m + 1 + 1 by ring]
  try rw [show m + 1 + 2 + m = 2 * m + 1 + 1 + 1 by ring]
  try rw [show m + 1 + 3 + m = 2 * m + 1 + 1 + 1 + 1 by ring]
  try rw [show m + 2 = m + 1 + 1 by ring]
  try rw [show m + 3 = m + 1 + 1 + 1 by ring]
  simp only [V05.Hr_succ]
  simp only [Hr1m m, Hr2m m, Hr3m m, Hr12m m, Hr22m m, Hr32m m]
  push_cast
  have h1 : (2 * ((m : ℝ) + 1) - 1) ≠ 0 := by linarith
  have h2 : (2 * (m : ℝ) + 1) ≠ 0 := by positivity
  field_simp
  ring

set_option maxRecDepth 20000 in
lemma T_2_3 (m : ℕ) : Tji 2 3 m =
      (-9 / 8 : ℝ) * (1 / ((m : ℝ) + 1) ^ 1) +
      (3 / 8 : ℝ) * (1 / ((m : ℝ) + 1) ^ 2) +
      (9 / 8 : ℝ) * (1 / (2 * (m : ℝ) + 1) ^ 1) +
      (3 / 8 : ℝ) * (1 / (2 * (m : ℝ) + 1) ^ 2) +
      (-9 / 8 : ℝ) * (1 / (2 * (m : ℝ) + 1) ^ 3) +
      (9 / 8 : ℝ) * (1 / (2 * (m : ℝ) + 3) ^ 1) +
      (-3 / 8 : ℝ) * (1 / (2 * (m : ℝ) + 3) ^ 2) +
      (-9 / 8 : ℝ) * (1 / (2 * (m : ℝ) + 3) ^ 3) +
      (3 / 4 : ℝ) * (1 / (2 * (m : ℝ) + 3) ^ 4) +
      (-3 / 8 : ℝ) * ((S03.Z [2] (2 * m) : ℚ) : ℝ) +
      (9 / 16 : ℝ) * (1 / (2 * (m : ℝ) + 1) ^ 1) * ((S03.Z [2] (2 * m) : ℚ) : ℝ) +
      (9 / 16 : ℝ) * (1 / (2 * (m : ℝ) + 3) ^ 1) * ((S03.Z [2] (2 * m) : ℚ) : ℝ) +
      (-3 / 8 : ℝ) * (1 / (2 * (m : ℝ) + 3) ^ 2) * ((S03.Z [2] (2 * m) : ℚ) : ℝ) +
      (-9 / 8 : ℝ) * ((S03.Z [2] (2 * m) : ℚ) : ℝ) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) +
      (-3 / 4 : ℝ) * ((m : ℝ) ^ 1) * ((S03.Z [2] (2 * m) : ℚ) : ℝ) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) +
      (3 / 8 : ℝ) * (1 / (2 * (m : ℝ) + 3) ^ 1) * ((S03.Z [2] (2 * m) : ℚ) : ℝ) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) +
      (9 / 4 : ℝ) * ((S03.Z [3] (2 * m) : ℚ) : ℝ) +
      (3 / 2 : ℝ) * ((m : ℝ) ^ 1) * ((S03.Z [3] (2 * m) : ℚ) : ℝ) +
      (-3 / 4 : ℝ) * (1 / (2 * (m : ℝ) + 3) ^ 1) * ((S03.Z [3] (2 * m) : ℚ) : ℝ) +
      (-3 / 4 : ℝ) * (1 / ((m : ℝ) + 1) ^ 1) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) +
      (-3 / 8 : ℝ) * (1 / (2 * (m : ℝ) + 1) ^ 1) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) +
      (9 / 8 : ℝ) * (1 / (2 * (m : ℝ) + 1) ^ 2) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) +
      (3 / 8 : ℝ) * (1 / (2 * (m : ℝ) + 3) ^ 1) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) +
      (9 / 8 : ℝ) * (1 / (2 * (m : ℝ) + 3) ^ 2) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) +
      (-3 / 4 : ℝ) * (1 / (2 * (m : ℝ) + 3) ^ 3) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) +
      (3 / 8 : ℝ) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) +
      (-9 / 16 : ℝ) * (1 / (2 * (m : ℝ) + 1) ^ 1) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) +
      (-9 / 16 : ℝ) * (1 / (2 * (m : ℝ) + 3) ^ 1) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) +
      (3 / 8 : ℝ) * (1 / (2 * (m : ℝ) + 3) ^ 2) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) +
      (3 / 8 : ℝ) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) +
      (1 / 4 : ℝ) * ((m : ℝ) ^ 1) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) +
      (-1 / 8 : ℝ) * (1 / (2 * (m : ℝ) + 3) ^ 1) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) +
      (-9 / 4 : ℝ) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) * ((S03.Z [-2] (2 * m) : ℚ) : ℝ) +
      (-3 / 2 : ℝ) * ((m : ℝ) ^ 1) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) * ((S03.Z [-2] (2 * m) : ℚ) : ℝ) +
      (3 / 4 : ℝ) * (1 / (2 * (m : ℝ) + 3) ^ 1) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) * ((S03.Z [-2] (2 * m) : ℚ) : ℝ) +
      (-3 / 4 : ℝ) * ((S03.Z [-2] (2 * m) : ℚ) : ℝ) +
      (9 / 8 : ℝ) * (1 / (2 * (m : ℝ) + 1) ^ 1) * ((S03.Z [-2] (2 * m) : ℚ) : ℝ) +
      (9 / 8 : ℝ) * (1 / (2 * (m : ℝ) + 3) ^ 1) * ((S03.Z [-2] (2 * m) : ℚ) : ℝ) +
      (-3 / 4 : ℝ) * (1 / (2 * (m : ℝ) + 3) ^ 2) * ((S03.Z [-2] (2 * m) : ℚ) : ℝ) +
      (3 : ℝ) * ((S03.Z [-3] (2 * m) : ℚ) : ℝ) +
      (2 : ℝ) * ((m : ℝ) ^ 1) * ((S03.Z [-3] (2 * m) : ℚ) : ℝ) +
      (-1 : ℝ) * (1 / (2 * (m : ℝ) + 3) ^ 1) * ((S03.Z [-3] (2 * m) : ℚ) : ℝ) := by
  have q1 : (0 : ℝ) ≤ m := by positivity
  simp only [Tji, V05.rho, Finset.prod_range_succ, Finset.prod_range_zero]
  norm_num
  simp only [e303_eq, e312_eq, e321_eq, e330_eq, V05.coeff3, V05.Dsharp]
  norm_num
  try rw [show m + 1 + m = 2 * m + 1 by ring]
  try rw [show m + 1 + 1 + m = 2 * m + 1 + 1 by ring]
  try rw [show m + 1 + 2 + m = 2 * m + 1 + 1 + 1 by ring]
  try rw [show m + 1 + 3 + m = 2 * m + 1 + 1 + 1 + 1 by ring]
  try rw [show m + 2 = m + 1 + 1 by ring]
  try rw [show m + 3 = m + 1 + 1 + 1 by ring]
  simp only [V05.Hr_succ]
  simp only [Hr1m m, Hr2m m, Hr3m m, Hr12m m, Hr22m m, Hr32m m]
  push_cast
  have h1 : (2 * ((m : ℝ) + 1) - 1) ≠ 0 := by linarith
  have h2 : (2 * (m : ℝ) + 1) ≠ 0 := by positivity
  field_simp
  ring

set_option maxRecDepth 20000 in
lemma T_3_0 (m : ℕ) : Tji 3 0 m =
      (7 / 16 : ℝ) * (1 / ((m : ℝ) + 1) ^ 1) +
      (-3 / 16 : ℝ) * (1 / ((m : ℝ) + 1) ^ 2) +
      (-1 / 8 : ℝ) * (1 / ((m : ℝ) + 1) ^ 3) +
      (-41 / 24 : ℝ) * (1 / (2 * (m : ℝ) + 1) ^ 1) +
      (-61 / 24 : ℝ) * (1 / (2 * (m : ℝ) + 1) ^ 2) +
      (-5 / 4 : ℝ) * (1 / (2 * (m : ℝ) + 1) ^ 3) +
      (-1 / 48 : ℝ) * (1 / ((m : ℝ) + 2) ^ 1) +
      (1 / 6 : ℝ) * (1 / ((m : ℝ) + 2) ^ 2) +
      (-1 / 16 : ℝ) * (1 / ((m : ℝ) + 2) ^ 3) +
      (7 / 8 : ℝ) * (1 / (2 * (m : ℝ) + 3) ^ 1) +
      (1 / 8 : ℝ) * (1 / (2 * (m : ℝ) + 3) ^ 2) +
      (-3 / 4 : ℝ) * (1 / (2 * (m : ℝ) + 3) ^ 3) +
      (1 / 2 : ℝ) * ((S03.Z [2] (2 * m) : ℚ) : ℝ) +
      (1 / 4 : ℝ) * (1 / ((m : ℝ) + 1) ^ 1) * ((S03.Z [2] (2 * m) : ℚ) : ℝ) +
      (5 / 8 : ℝ) * (1 / (2 * (m : ℝ) + 1) ^ 1) * ((S03.Z [2] (2 * m) : ℚ) : ℝ) +
      (1 / 8 : ℝ) * (1 / ((m : ℝ) + 2) ^ 1) * ((S03.Z [2] (2 * m) : ℚ) : ℝ) +
      (3 / 8 : ℝ) * (1 / (2 * (m : ℝ) + 3) ^ 1) * ((S03.Z [2] (2 * m) : ℚ) : ℝ) +
      (-3 / 4 : ℝ) * ((S03.Z [2] (2 * m) : ℚ) : ℝ) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) +
      (-1 / 4 : ℝ) * ((m : ℝ) ^ 1) * ((S03.Z [2] (2 * m) : ℚ) : ℝ) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) +
      (3 / 2 : ℝ) * ((S03.Z [3] (2 * m) : ℚ) : ℝ) +
      (1 / 2 : ℝ) * ((m : ℝ) ^ 1) * ((S03.Z [3] (2 * m) : ℚ) : ℝ) +
      (3 / 8 : ℝ) * (1 / ((m : ℝ) + 1) ^ 1) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) +
      (1 / 4 : ℝ) * (1 / ((m : ℝ) + 1) ^ 2) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) +
      (61 / 24 : ℝ) * (1 / (2 * (m : ℝ) + 1) ^ 1) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) +
      (5 / 4 : ℝ) * (1 / (2 * (m : ℝ) + 1) ^ 2) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) +
      (-1 / 3 : ℝ) * (1 / ((m : ℝ) + 2) ^ 1) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) +
      (1 / 8 : ℝ) * (1 / ((m : ℝ) + 2) ^ 2) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) +
      (-1 / 8 : ℝ) * (1 / (2 * (m : ℝ) + 3) ^ 1) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) +
      (3 / 4 : ℝ) * (1 / (2 * (m : ℝ) + 3) ^ 2) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) +
      (-1 / 2 : ℝ) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) +
      (-1 / 4 : ℝ) * (1 / ((m : ℝ) + 1) ^ 1) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) +
      (-5 / 8 : ℝ) * (1 / (2 * (m : ℝ) + 1) ^ 1) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) +
      (-1 / 8 : ℝ) * (1 / ((m : ℝ) + 2) ^ 1) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) +
      (-3 / 8 : ℝ) * (1 / (2 * (m : ℝ) + 3) ^ 1) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) +
      (1 / 4 : ℝ) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) +
      (1 / 12 : ℝ) * ((m : ℝ) ^ 1) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) +
      (-3 / 2 : ℝ) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) * ((S03.Z [-2] (2 * m) : ℚ) : ℝ) +
      (-1 / 2 : ℝ) * ((m : ℝ) ^ 1) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) * ((S03.Z [-2] (2 * m) : ℚ) : ℝ) +
      (1 : ℝ) * ((S03.Z [-2] (2 * m) : ℚ) : ℝ) +
      (1 / 2 : ℝ) * (1 / ((m : ℝ) + 1) ^ 1) * ((S03.Z [-2] (2 * m) : ℚ) : ℝ) +
      (5 / 4 : ℝ) * (1 / (2 * (m : ℝ) + 1) ^ 1) * ((S03.Z [-2] (2 * m) : ℚ) : ℝ) +
      (1 / 4 : ℝ) * (1 / ((m : ℝ) + 2) ^ 1) * ((S03.Z [-2] (2 * m) : ℚ) : ℝ) +
      (3 / 4 : ℝ) * (1 / (2 * (m : ℝ) + 3) ^ 1) * ((S03.Z [-2] (2 * m) : ℚ) : ℝ) +
      (2 : ℝ) * ((S03.Z [-3] (2 * m) : ℚ) : ℝ) +
      (2 / 3 : ℝ) * ((m : ℝ) ^ 1) * ((S03.Z [-3] (2 * m) : ℚ) : ℝ) := by
  have q1 : (0 : ℝ) ≤ m := by positivity
  simp only [Tji, V05.rho, Finset.prod_range_succ, Finset.prod_range_zero]
  norm_num
  simp only [e303_eq, e312_eq, e321_eq, e330_eq, V05.coeff3, V05.Dsharp]
  norm_num
  try rw [show m + 1 + m = 2 * m + 1 by ring]
  try rw [show m + 1 + 1 + m = 2 * m + 1 + 1 by ring]
  try rw [show m + 1 + 2 + m = 2 * m + 1 + 1 + 1 by ring]
  try rw [show m + 1 + 3 + m = 2 * m + 1 + 1 + 1 + 1 by ring]
  try rw [show m + 2 = m + 1 + 1 by ring]
  try rw [show m + 3 = m + 1 + 1 + 1 by ring]
  simp only [V05.Hr_succ]
  simp only [Hr1m m, Hr2m m, Hr3m m, Hr12m m, Hr22m m, Hr32m m]
  push_cast
  have h1 : (2 * ((m : ℝ) + 1) - 1) ≠ 0 := by linarith
  have h2 : (2 * (m : ℝ) + 1) ≠ 0 := by positivity
  field_simp
  ring

set_option maxRecDepth 20000 in
lemma T_3_1 (m : ℕ) : Tji 3 1 m =
      (-11 / 16 : ℝ) * (1 / ((m : ℝ) + 1) ^ 1) +
      (-3 / 16 : ℝ) * (1 / ((m : ℝ) + 1) ^ 2) +
      (1 / 8 : ℝ) * (1 / ((m : ℝ) + 1) ^ 3) +
      (5 / 24 : ℝ) * (1 / (2 * (m : ℝ) + 1) ^ 1) +
      (91 / 24 : ℝ) * (1 / (2 * (m : ℝ) + 1) ^ 2) +
      (15 / 4 : ℝ) * (1 / (2 * (m : ℝ) + 1) ^ 3) +
      (13 / 48 : ℝ) * (1 / ((m : ℝ) + 2) ^ 1) +
      (-1 / 6 : ℝ) * (1 / ((m : ℝ) + 2) ^ 2) +
      (1 / 16 : ℝ) * (1 / ((m : ℝ) + 2) ^ 3) +
      (5 / 8 : ℝ) * (1 / (2 * (m : ℝ) + 3) ^ 1) +
      (-7 / 8 : ℝ) * (1 / (2 * (m : ℝ) + 3) ^ 2) +
      (9 / 4 : ℝ) * (1 / (2 * (m : ℝ) + 3) ^ 3) +
      (5 / 4 : ℝ) * ((S03.Z [2] (2 * m) : ℚ) : ℝ) +
      (3 / 4 : ℝ) * (1 / ((m : ℝ) + 1) ^ 1) * ((S03.Z [2] (2 * m) : ℚ) : ℝ) +
      (5 / 8 : ℝ) * (1 / (2 * (m : ℝ) + 1) ^ 1) * ((S03.Z [2] (2 * m) : ℚ) : ℝ) +
      (3 / 8 : ℝ) * (1 / ((m : ℝ) + 2) ^ 1) * ((S03.Z [2] (2 * m) : ℚ) : ℝ) +
      (3 / 8 : ℝ) * (1 / (2 * (m : ℝ) + 3) ^ 1) * ((S03.Z [2] (2 * m) : ℚ) : ℝ) +
      (-3 / 4 : ℝ) * ((S03.Z [2] (2 * m) : ℚ) : ℝ) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) +
      (-1 / 4 : ℝ) * ((m : ℝ) ^ 1) * ((S03.Z [2] (2 * m) : ℚ) : ℝ) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) +
      (3 / 2 : ℝ) * ((S03.Z [3] (2 * m) : ℚ) : ℝ) +
      (1 / 2 : ℝ) * ((m : ℝ) ^ 1) * ((S03.Z [3] (2 * m) : ℚ) : ℝ) +
      (3 / 8 : ℝ) * (1 / ((m : ℝ) + 1) ^ 1) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) +
      (-1 / 4 : ℝ) * (1 / ((m : ℝ) + 1) ^ 2) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) +
      (-91 / 24 : ℝ) * (1 / (2 * (m : ℝ) + 1) ^ 1) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) +
      (-15 / 4 : ℝ) * (1 / (2 * (m : ℝ) + 1) ^ 2) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) +
      (1 / 3 : ℝ) * (1 / ((m : ℝ) + 2) ^ 1) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) +
      (-1 / 8 : ℝ) * (1 / ((m : ℝ) + 2) ^ 2) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) +
      (7 / 8 : ℝ) * (1 / (2 * (m : ℝ) + 3) ^ 1) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) +
      (-9 / 4 : ℝ) * (1 / (2 * (m : ℝ) + 3) ^ 2) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) +
      (3 / 4 : ℝ) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) +
      (1 / 4 : ℝ) * (1 / ((m : ℝ) + 1) ^ 1) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) +
      (15 / 8 : ℝ) * (1 / (2 * (m : ℝ) + 1) ^ 1) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) +
      (1 / 8 : ℝ) * (1 / ((m : ℝ) + 2) ^ 1) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) +
      (9 / 8 : ℝ) * (1 / (2 * (m : ℝ) + 3) ^ 1) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) +
      (-3 / 4 : ℝ) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) +
      (-1 / 4 : ℝ) * ((m : ℝ) ^ 1) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) +
      (3 / 2 : ℝ) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) * ((S03.Z [-2] (2 * m) : ℚ) : ℝ) +
      (1 / 2 : ℝ) * ((m : ℝ) ^ 1) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) * ((S03.Z [-2] (2 * m) : ℚ) : ℝ) +
      (1 / 2 : ℝ) * ((S03.Z [-2] (2 * m) : ℚ) : ℝ) +
      (1 / 2 : ℝ) * (1 / ((m : ℝ) + 1) ^ 1) * ((S03.Z [-2] (2 * m) : ℚ) : ℝ) +
      (-5 / 4 : ℝ) * (1 / (2 * (m : ℝ) + 1) ^ 1) * ((S03.Z [-2] (2 * m) : ℚ) : ℝ) +
      (1 / 4 : ℝ) * (1 / ((m : ℝ) + 2) ^ 1) * ((S03.Z [-2] (2 * m) : ℚ) : ℝ) +
      (-3 / 4 : ℝ) * (1 / (2 * (m : ℝ) + 3) ^ 1) * ((S03.Z [-2] (2 * m) : ℚ) : ℝ) := by
  have q1 : (0 : ℝ) ≤ m := by positivity
  simp only [Tji, V05.rho, Finset.prod_range_succ, Finset.prod_range_zero]
  norm_num
  simp only [e303_eq, e312_eq, e321_eq, e330_eq, V05.coeff3, V05.Dsharp]
  norm_num
  try rw [show m + 1 + m = 2 * m + 1 by ring]
  try rw [show m + 1 + 1 + m = 2 * m + 1 + 1 by ring]
  try rw [show m + 1 + 2 + m = 2 * m + 1 + 1 + 1 by ring]
  try rw [show m + 1 + 3 + m = 2 * m + 1 + 1 + 1 + 1 by ring]
  try rw [show m + 2 = m + 1 + 1 by ring]
  try rw [show m + 3 = m + 1 + 1 + 1 by ring]
  simp only [V05.Hr_succ]
  simp only [Hr1m m, Hr2m m, Hr3m m, Hr12m m, Hr22m m, Hr32m m]
  push_cast
  have h1 : (2 * ((m : ℝ) + 1) - 1) ≠ 0 := by linarith
  have h2 : (2 * (m : ℝ) + 1) ≠ 0 := by positivity
  field_simp
  ring

set_option maxRecDepth 20000 in
lemma T_3_2 (m : ℕ) : Tji 3 2 m =
      (-11 / 16 : ℝ) * (1 / ((m : ℝ) + 1) ^ 1) +
      (3 / 16 : ℝ) * (1 / ((m : ℝ) + 1) ^ 2) +
      (1 / 8 : ℝ) * (1 / ((m : ℝ) + 1) ^ 3) +
      (41 / 24 : ℝ) * (1 / (2 * (m : ℝ) + 1) ^ 1) +
      (1 / 24 : ℝ) * (1 / (2 * (m : ℝ) + 1) ^ 2) +
      (-15 / 4 : ℝ) * (1 / (2 * (m : ℝ) + 1) ^ 3) +
      (13 / 48 : ℝ) * (1 / ((m : ℝ) + 2) ^ 1) +
      (-1 / 6 : ℝ) * (1 / ((m : ℝ) + 2) ^ 2) +
      (1 / 16 : ℝ) * (1 / ((m : ℝ) + 2) ^ 3) +
      (-7 / 8 : ℝ) * (1 / (2 * (m : ℝ) + 3) ^ 1) +
      (11 / 8 : ℝ) * (1 / (2 * (m : ℝ) + 3) ^ 2) +
      (-9 / 4 : ℝ) * (1 / (2 * (m : ℝ) + 3) ^ 3) +
      (1 : ℝ) * ((S03.Z [2] (2 * m) : ℚ) : ℝ) +
      (3 / 4 : ℝ) * (1 / ((m : ℝ) + 1) ^ 1) * ((S03.Z [2] (2 * m) : ℚ) : ℝ) +
      (-5 / 8 : ℝ) * (1 / (2 * (m : ℝ) + 1) ^ 1) * ((S03.Z [2] (2 * m) : ℚ) : ℝ) +
      (3 / 8 : ℝ) * (1 / ((m : ℝ) + 2) ^ 1) * ((S03.Z [2] (2 * m) : ℚ) : ℝ) +
      (-3 / 8 : ℝ) * (1 / (2 * (m : ℝ) + 3) ^ 1) * ((S03.Z [2] (2 * m) : ℚ) : ℝ) +
      (3 / 4 : ℝ) * ((S03.Z [2] (2 * m) : ℚ) : ℝ) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) +
      (1 / 4 : ℝ) * ((m : ℝ) ^ 1) * ((S03.Z [2] (2 * m) : ℚ) : ℝ) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) +
      (-3 / 2 : ℝ) * ((S03.Z [3] (2 * m) : ℚ) : ℝ) +
      (-1 / 2 : ℝ) * ((m : ℝ) ^ 1) * ((S03.Z [3] (2 * m) : ℚ) : ℝ) +
      (-3 / 8 : ℝ) * (1 / ((m : ℝ) + 1) ^ 1) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) +
      (-1 / 4 : ℝ) * (1 / ((m : ℝ) + 1) ^ 2) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) +
      (-1 / 24 : ℝ) * (1 / (2 * (m : ℝ) + 1) ^ 1) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) +
      (15 / 4 : ℝ) * (1 / (2 * (m : ℝ) + 1) ^ 2) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) +
      (1 / 3 : ℝ) * (1 / ((m : ℝ) + 2) ^ 1) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) +
      (-1 / 8 : ℝ) * (1 / ((m : ℝ) + 2) ^ 2) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) +
      (-11 / 8 : ℝ) * (1 / (2 * (m : ℝ) + 3) ^ 1) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) +
      (9 / 4 : ℝ) * (1 / (2 * (m : ℝ) + 3) ^ 2) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) +
      (1 / 4 : ℝ) * (1 / ((m : ℝ) + 1) ^ 1) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) +
      (-15 / 8 : ℝ) * (1 / (2 * (m : ℝ) + 1) ^ 1) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) +
      (1 / 8 : ℝ) * (1 / ((m : ℝ) + 2) ^ 1) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) +
      (-9 / 8 : ℝ) * (1 / (2 * (m : ℝ) + 3) ^ 1) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) +
      (3 / 4 : ℝ) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) +
      (1 / 4 : ℝ) * ((m : ℝ) ^ 1) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) +
      (-3 / 2 : ℝ) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) * ((S03.Z [-2] (2 * m) : ℚ) : ℝ) +
      (-1 / 2 : ℝ) * ((m : ℝ) ^ 1) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) * ((S03.Z [-2] (2 * m) : ℚ) : ℝ) +
      (1 : ℝ) * ((S03.Z [-2] (2 * m) : ℚ) : ℝ) +
      (1 / 2 : ℝ) * (1 / ((m : ℝ) + 1) ^ 1) * ((S03.Z [-2] (2 * m) : ℚ) : ℝ) +
      (5 / 4 : ℝ) * (1 / (2 * (m : ℝ) + 1) ^ 1) * ((S03.Z [-2] (2 * m) : ℚ) : ℝ) +
      (1 / 4 : ℝ) * (1 / ((m : ℝ) + 2) ^ 1) * ((S03.Z [-2] (2 * m) : ℚ) : ℝ) +
      (3 / 4 : ℝ) * (1 / (2 * (m : ℝ) + 3) ^ 1) * ((S03.Z [-2] (2 * m) : ℚ) : ℝ) := by
  have q1 : (0 : ℝ) ≤ m := by positivity
  simp only [Tji, V05.rho, Finset.prod_range_succ, Finset.prod_range_zero]
  norm_num
  simp only [e303_eq, e312_eq, e321_eq, e330_eq, V05.coeff3, V05.Dsharp]
  norm_num
  try rw [show m + 1 + m = 2 * m + 1 by ring]
  try rw [show m + 1 + 1 + m = 2 * m + 1 + 1 by ring]
  try rw [show m + 1 + 2 + m = 2 * m + 1 + 1 + 1 by ring]
  try rw [show m + 1 + 3 + m = 2 * m + 1 + 1 + 1 + 1 by ring]
  try rw [show m + 2 = m + 1 + 1 by ring]
  try rw [show m + 3 = m + 1 + 1 + 1 by ring]
  simp only [V05.Hr_succ]
  simp only [Hr1m m, Hr2m m, Hr3m m, Hr12m m, Hr22m m, Hr32m m]
  push_cast
  have h1 : (2 * ((m : ℝ) + 1) - 1) ≠ 0 := by linarith
  have h2 : (2 * (m : ℝ) + 1) ≠ 0 := by positivity
  field_simp
  ring

set_option maxRecDepth 20000 in
lemma T_3_3 (m : ℕ) : Tji 3 3 m =
      (7 / 16 : ℝ) * (1 / ((m : ℝ) + 1) ^ 1) +
      (3 / 16 : ℝ) * (1 / ((m : ℝ) + 1) ^ 2) +
      (-1 / 8 : ℝ) * (1 / ((m : ℝ) + 1) ^ 3) +
      (-5 / 24 : ℝ) * (1 / (2 * (m : ℝ) + 1) ^ 1) +
      (-31 / 24 : ℝ) * (1 / (2 * (m : ℝ) + 1) ^ 2) +
      (5 / 4 : ℝ) * (1 / (2 * (m : ℝ) + 1) ^ 3) +
      (-1 / 48 : ℝ) * (1 / ((m : ℝ) + 2) ^ 1) +
      (1 / 6 : ℝ) * (1 / ((m : ℝ) + 2) ^ 2) +
      (-1 / 16 : ℝ) * (1 / ((m : ℝ) + 2) ^ 3) +
      (-5 / 8 : ℝ) * (1 / (2 * (m : ℝ) + 3) ^ 1) +
      (-5 / 8 : ℝ) * (1 / (2 * (m : ℝ) + 3) ^ 2) +
      (3 / 4 : ℝ) * (1 / (2 * (m : ℝ) + 3) ^ 3) +
      (1 / 4 : ℝ) * ((S03.Z [2] (2 * m) : ℚ) : ℝ) +
      (1 / 4 : ℝ) * (1 / ((m : ℝ) + 1) ^ 1) * ((S03.Z [2] (2 * m) : ℚ) : ℝ) +
      (-5 / 8 : ℝ) * (1 / (2 * (m : ℝ) + 1) ^ 1) * ((S03.Z [2] (2 * m) : ℚ) : ℝ) +
      (1 / 8 : ℝ) * (1 / ((m : ℝ) + 2) ^ 1) * ((S03.Z [2] (2 * m) : ℚ) : ℝ) +
      (-3 / 8 : ℝ) * (1 / (2 * (m : ℝ) + 3) ^ 1) * ((S03.Z [2] (2 * m) : ℚ) : ℝ) +
      (3 / 4 : ℝ) * ((S03.Z [2] (2 * m) : ℚ) : ℝ) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) +
      (1 / 4 : ℝ) * ((m : ℝ) ^ 1) * ((S03.Z [2] (2 * m) : ℚ) : ℝ) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) +
      (-3 / 2 : ℝ) * ((S03.Z [3] (2 * m) : ℚ) : ℝ) +
      (-1 / 2 : ℝ) * ((m : ℝ) ^ 1) * ((S03.Z [3] (2 * m) : ℚ) : ℝ) +
      (-3 / 8 : ℝ) * (1 / ((m : ℝ) + 1) ^ 1) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) +
      (1 / 4 : ℝ) * (1 / ((m : ℝ) + 1) ^ 2) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) +
      (31 / 24 : ℝ) * (1 / (2 * (m : ℝ) + 1) ^ 1) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) +
      (-5 / 4 : ℝ) * (1 / (2 * (m : ℝ) + 1) ^ 2) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) +
      (-1 / 3 : ℝ) * (1 / ((m : ℝ) + 2) ^ 1) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) +
      (1 / 8 : ℝ) * (1 / ((m : ℝ) + 2) ^ 2) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) +
      (5 / 8 : ℝ) * (1 / (2 * (m : ℝ) + 3) ^ 1) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) +
      (-3 / 4 : ℝ) * (1 / (2 * (m : ℝ) + 3) ^ 2) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) +
      (-1 / 4 : ℝ) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) +
      (-1 / 4 : ℝ) * (1 / ((m : ℝ) + 1) ^ 1) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) +
      (5 / 8 : ℝ) * (1 / (2 * (m : ℝ) + 1) ^ 1) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) +
      (-1 / 8 : ℝ) * (1 / ((m : ℝ) + 2) ^ 1) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) +
      (3 / 8 : ℝ) * (1 / (2 * (m : ℝ) + 3) ^ 1) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) +
      (-1 / 4 : ℝ) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) +
      (-1 / 12 : ℝ) * ((m : ℝ) ^ 1) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) +
      (3 / 2 : ℝ) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) * ((S03.Z [-2] (2 * m) : ℚ) : ℝ) +
      (1 / 2 : ℝ) * ((m : ℝ) ^ 1) * ((S03.Z [-1] (2 * m) : ℚ) : ℝ) * ((S03.Z [-2] (2 * m) : ℚ) : ℝ) +
      (1 / 2 : ℝ) * ((S03.Z [-2] (2 * m) : ℚ) : ℝ) +
      (1 / 2 : ℝ) * (1 / ((m : ℝ) + 1) ^ 1) * ((S03.Z [-2] (2 * m) : ℚ) : ℝ) +
      (-5 / 4 : ℝ) * (1 / (2 * (m : ℝ) + 1) ^ 1) * ((S03.Z [-2] (2 * m) : ℚ) : ℝ) +
      (1 / 4 : ℝ) * (1 / ((m : ℝ) + 2) ^ 1) * ((S03.Z [-2] (2 * m) : ℚ) : ℝ) +
      (-3 / 4 : ℝ) * (1 / (2 * (m : ℝ) + 3) ^ 1) * ((S03.Z [-2] (2 * m) : ℚ) : ℝ) +
      (-2 : ℝ) * ((S03.Z [-3] (2 * m) : ℚ) : ℝ) +
      (-2 / 3 : ℝ) * ((m : ℝ) ^ 1) * ((S03.Z [-3] (2 * m) : ℚ) : ℝ) := by
  have q1 : (0 : ℝ) ≤ m := by positivity
  simp only [Tji, V05.rho, Finset.prod_range_succ, Finset.prod_range_zero]
  norm_num
  simp only [e303_eq, e312_eq, e321_eq, e330_eq, V05.coeff3, V05.Dsharp]
  norm_num
  try rw [show m + 1 + m = 2 * m + 1 by ring]
  try rw [show m + 1 + 1 + m = 2 * m + 1 + 1 by ring]
  try rw [show m + 1 + 2 + m = 2 * m + 1 + 1 + 1 by ring]
  try rw [show m + 1 + 3 + m = 2 * m + 1 + 1 + 1 + 1 by ring]
  try rw [show m + 2 = m + 1 + 1 by ring]
  try rw [show m + 3 = m + 1 + 1 + 1 by ring]
  simp only [V05.Hr_succ]
  simp only [Hr1m m, Hr2m m, Hr3m m, Hr12m m, Hr22m m, Hr32m m]
  push_cast
  have h1 : (2 * ((m : ℝ) + 1) - 1) ≠ 0 := by linarith
  have h2 : (2 * (m : ℝ) + 1) ≠ 0 := by positivity
  field_simp
  ring


lemma part1 (m : ℕ) : ((coefTerms R3terms (2 * ((m + 1 : ℕ) : ℤ)) : ℚ) : ℝ) = V05.P 3 (m + 1) 0 0 := by
  rw [P_split]
  simp only [Finset.sum_range_succ, Finset.sum_range_zero, zero_add]
  simp only [T_0_0, T_0_1, T_0_2, T_0_3, T_1_0, T_1_1, T_1_2, T_1_3, T_2_0, T_2_1, T_2_2, T_2_3,
    T_3_0, T_3_1, T_3_2, T_3_3]
  unfold coefTerms R3terms
  simp only [List.map_cons, List.map_nil, List.sum_cons, List.sum_nil]
  simp only [hE0, hE1, hE2, hE3, hE4, hE5, hE6, hE7, hE8, hE9, hE10, hE11, hE12, hE13, hE14, hE15, hE16, hE17, hE18, hE19, hE20, hE21, hE22, hE23, hE24, hE25, hE26, hE27, hE28, hE29, hE30, hE31, hE32, hE33, hE34, hE35, hE36, hE37, hE38, hE39, hE40, hE41, hE42, hE43, hE44, hE45, hE46, hE47, hE48, hE49, hE50, hE51, hE52, hE53, hE54, hE55, hE56, hE57, hE58, hE59, hE60, hE61, hE62, hE63, hE64, hE65, hE66, hE67, hE68, hE69, hE70, hE71, hE72, hE73]
  rw [t4Q, t2Q, t3Q, t1Q]
  push_cast
  ring

lemma part2 (m : ℕ) : coefTerms R3terms (2 * ((m + 1 : ℕ) : ℤ) + 1) = 0 := by
  unfold coefTerms R3terms
  simp only [List.map_cons, List.map_nil, List.sum_cons, List.sum_nil]
  simp only [hO0, hO1, hO2, hO3, hO4, hO5, hO6, hO7, hO8, hO9, hO10, hO11, hO12, hO13, hO14, hO15, hO16, hO17, hO18, hO19, hO20, hO21, hO22, hO23, hO24, hO25, hO26, hO27, hO28, hO29, hO30, hO31, hO32, hO33, hO34, hO35, hO36, hO37, hO38, hO39, hO40, hO41, hO42, hO43, hO44, hO45, hO46, hO47, hO48, hO49, hO50, hO51, hO52, hO53, hO54, hO55, hO56, hO57, hO58, hO59, hO60, hO61, hO62, hO63, hO64, hO65, hO66, hO67, hO68, hO69, hO70, hO71, hO72, hO73]
  rw [t4Q, t2Q, t3Q, t1Q]
  ring

end P02.H11

theorem P02.proof_11 : P02.stmt_11 := by
  refine ⟨fun n hn => ?_, fun n hn => ?_, fun N hN => ?_⟩
  · obtain ⟨m, rfl⟩ : ∃ m, n = m + 1 := ⟨n - 1, by omega⟩
    exact P02.H11.part1 m
  · obtain ⟨m, rfl⟩ : ∃ m, n = m + 1 := ⟨n - 1, by omega⟩
    exact P02.H11.part2 m
  · rcases (show N ≤ -2 ∨ N = -1 ∨ N = 0 ∨ N = 1 by omega) with h | rfl | rfl | rfl
    · exact P02.H12.coefTerms_low _ _ (fun τ hτ => by
        have := (P02.proof_1.2.2.2.2.2.1 τ hτ).2.2.1; omega)
    · decide +kernel
    · decide +kernel
    · decide +kernel

theorem P02.proof_12 : P02.stmt_12 := by
  refine ⟨fun n hn => ?_, fun n hn => ?_, P02.H12.part3⟩
  · obtain ⟨m, rfl⟩ : ∃ m, n = m + 1 := ⟨n - 1, by omega⟩
    exact P02.H12.part1 m
  · obtain ⟨m, rfl⟩ : ∃ m, n = m + 1 := ⟨n - 1, by omega⟩
    exact P02.H12.part2 m

namespace P02.H13

open P02 H12

lemma termHasSum (q : ℚ) (e : ℤ) (w : List ℤ) (hw : w ≠ []) (hl : S03.IsLetters w)
    (hr : S03.Regular w) (he : -2 ≤ e) (c : ℝ) (hc0 : 0 < c) (hc1 : c < 1) :
    HasSum (fun k : ℕ => (q : ℝ) * (hcoef w ((k : ℤ) - 2 - e) : ℝ) * c ^ ((k : ℤ) - 2))
      ((q : ℝ) * c ^ e * S03.H w c) := by
  obtain ⟨j, hj⟩ : ∃ j : ℕ, (j : ℤ) = e + 3 := ⟨(e + 3).toNat, by omega⟩
  have hH := (S03.proof_5.2.1 w hw hl hr).2.2.2.2.2 c (by rw [abs_lt]; constructor <;> linarith)
  rw [← hasSum_nat_add_iff' j]
  have hz : ∑ i ∈ Finset.range j,
      (q : ℝ) * (hcoef w ((i : ℤ) - 2 - e) : ℝ) * c ^ ((i : ℤ) - 2) = 0 := by
    refine Finset.sum_eq_zero fun i hi => ?_
    have : (i : ℤ) < j := by exact_mod_cast Finset.mem_range.mp hi
    rw [hcoef_nonpos _ _ (by omega)]; simp
  rw [hz, sub_zero]
  convert hH.mul_left ((q : ℝ) * c ^ e) using 1
  funext N
  rw [show ((N + j : ℕ) : ℤ) - 2 - e = (N : ℤ) + 1 by push_cast; omega, hcoef_succ,
    show ((N + j : ℕ) : ℤ) - 2 = e + ((N + 1 : ℕ) : ℤ) by push_cast; omega,
    zpow_add₀ hc0.ne', zpow_natCast]
  push_cast
  ring

lemma listHasSum (L : List (ℚ × ℤ × List ℤ))
    (hL : ∀ τ ∈ L, τ.2.2 ≠ [] ∧ S03.IsLetters τ.2.2 ∧ S03.Regular τ.2.2 ∧ -2 ≤ τ.2.1)
    (c : ℝ) (hc0 : 0 < c) (hc1 : c < 1) :
    HasSum (fun k : ℕ => ((coefTerms L ((k : ℤ) - 2) : ℚ) : ℝ) * c ^ ((k : ℤ) - 2))
      ((L.map fun τ => (τ.1 : ℝ) * c ^ τ.2.1 * S03.H τ.2.2 c).sum) := by
  induction L with
  | nil => simp [coefTerms]
  | cons τ L ih =>
    have h1 := hL τ List.mem_cons_self
    have h2 := ih (fun σ hσ => hL σ (List.mem_cons_of_mem _ hσ))
    have h3 := termHasSum τ.1 τ.2.1 τ.2.2 h1.1 h1.2.1 h1.2.2.1 h1.2.2.2 c hc0 hc1
    rw [List.map_cons, List.sum_cons]
    convert h3.add h2 using 1
    funext k
    simp only [coefTerms, List.map_cons, List.sum_cons]
    push_cast
    ring

lemma R3_ok : ∀ τ ∈ R3terms, τ.2.2 ≠ [] ∧ S03.IsLetters τ.2.2 ∧ S03.Regular τ.2.2 ∧ -2 ≤ τ.2.1 := by
  intro τ hτ
  have := P02.proof_1.2.2.2.2.2.1 τ hτ
  refine ⟨?_, this.1, this.2.1, this.2.2.1⟩
  have hne : ∀ σ ∈ R3terms, σ.2.2 ≠ [] := by decide +kernel
  exact hne τ hτ

lemma Psi_ok : ∀ τ ∈ PsiTerms, τ.2.2 ≠ [] ∧ S03.IsLetters τ.2.2 ∧ S03.Regular τ.2.2 ∧ -2 ≤ τ.2.1 := by
  intro τ hτ
  have := P02.proof_1.2.2.2.2.2.2.2 τ hτ
  have hne : ∀ σ ∈ PsiTerms, σ.2.2 ≠ [] ∧ -2 ≤ σ.2.1 := by decide +kernel
  exact ⟨(hne τ hτ).1, this.1, this.2.1, (hne τ hτ).2⟩

/-- Reindexing: a sequence on `ℕ` supported on `k = 2n + 4` with values `b n`. -/
lemma reindex (f : ℕ → ℝ) (b : ℕ → ℝ) (S : ℝ) (hf : HasSum f S)
    (h1 : ∀ n : ℕ, f (2 * n + 4) = b n) (h0 : ∀ k, k ∉ Set.range (fun n : ℕ => 2 * n + 4) → f k = 0) :
    HasSum b S := by
  have hinj : Function.Injective (fun n : ℕ => 2 * n + 4) := fun a b h => by simp at h; omega
  have := (hinj.hasSum_iff h0).mpr hf
  convert this using 1
  funext n; exact (h1 n).symm

lemma not_range (k : ℕ) (hk : k ∉ Set.range (fun n : ℕ => 2 * n + 4)) :
    k < 4 ∨ ∃ n : ℕ, 1 ≤ n ∧ k = 2 * n + 3 := by
  by_contra h
  push_neg at h
  apply hk
  obtain ⟨h4, h5⟩ := h
  rcases Nat.even_or_odd k with ⟨r, hr⟩ | ⟨r, hr⟩
  · exact ⟨r - 2, show 2 * (r - 2) + 4 = k by omega⟩
  · exact absurd (h5 (r - 1) (by omega)) (by omega)

end P02.H13

theorem P02.proof_13 : P02.stmt_13 := by
  intro c hc0 hc1
  obtain ⟨p1, p2, p3⟩ := P02.proof_11
  obtain ⟨s1, s2, s3⟩ := P02.proof_12
  constructor
  · have hf := P02.H13.listHasSum P02.R3terms P02.H13.R3_ok c hc0 hc1
    refine P02.H13.reindex _ _ _ hf (fun n => ?_) (fun k hk => ?_)
    · rw [show ((2 * n + 4 : ℕ) : ℤ) - 2 = 2 * ((n + 1 : ℕ) : ℤ) by push_cast; ring, p1 (n + 1) (by omega),
        show (2 * ((n + 1 : ℕ) : ℤ)) = ((2 * (n + 1) : ℕ) : ℤ) by push_cast; ring, zpow_natCast]
    · rcases P02.H13.not_range k hk with h | ⟨n, hn, rfl⟩
      · rw [p3 _ (by omega)]; simp
      · rw [show ((2 * n + 3 : ℕ) : ℤ) - 2 = 2 * (n : ℤ) + 1 by push_cast; ring, p2 n hn]; simp
  · have hf := (P02.H13.listHasSum P02.PsiTerms P02.H13.Psi_ok c hc0 hc1).add
      (hasSum_ite_eq (2 : ℕ) (1 / 2 : ℝ))
    have e : P02.PsiFun c = (P02.PsiTerms.map fun τ => (τ.1 : ℝ) * c ^ τ.2.1 * S03.H τ.2.2 c).sum
        + 1 / 2 := by unfold P02.PsiFun; ring
    rw [e]
    refine P02.H13.reindex _ _ _ hf (fun n => ?_) (fun k hk => ?_)
    · rw [ite_eq_right_iff.mpr (fun h => by omega),
        show ((2 * n + 4 : ℕ) : ℤ) - 2 = 2 * ((n + 1 : ℕ) : ℤ) by push_cast; ring, s1 (n + 1) (by omega),
        show (2 * ((n + 1 : ℕ) : ℤ)) = ((2 * (n + 1) : ℕ) : ℤ) by push_cast; ring, zpow_natCast]
      ring
    · rcases P02.H13.not_range k hk with h | ⟨n, hn, rfl⟩
      · have := s3 ((k : ℤ) - 2) (by omega)
        by_cases hk2 : k = 2
        · subst hk2
          rw [if_pos rfl]
          have h' : ((P02.coefTerms P02.PsiTerms (((2 : ℕ) : ℤ) - 2)) : ℝ) = -(1 / 2) := by
            rw [if_pos (by norm_num)] at this
            have : P02.coefTerms P02.PsiTerms (((2 : ℕ) : ℤ) - 2) = -(1 / 2) := by linarith
            rw [this]; norm_num
          rw [h']; norm_num
        · rw [if_neg hk2]
          rw [if_neg (by omega)] at this
          rw [show P02.coefTerms P02.PsiTerms ((k : ℤ) - 2) = 0 by linarith]; simp
      · rw [if_neg (by omega),
          show ((2 * n + 3 : ℕ) : ℤ) - 2 = 2 * (n : ℤ) + 1 by push_cast; ring, s2 n hn]; simp
theorem P02.proof_14 : P02.stmt_14 := by
  intro α hα0 hα1
  obtain ⟨hc0, hc1, -, -, -⟩ := D04.cs_facts hα0 hα1
  have h3 := V05.proof_8 α hα0 hα1 3 (by norm_num)
  have hPhi := V05.proof_1 (D04.c α) (by rw [abs_of_pos hc0]; exact hc0)
    (by rw [abs_of_pos hc0]; exact hc1)
  obtain ⟨hR, hP⟩ := P02.proof_13 (D04.c α) hc0 hc1
  have hsum := (hR.sub (hPhi.mul_left V05.z3)).add (hP.mul_left V05.z2)
  refine h3.unique ?_
  convert hsum using 1
  funext n
  rw [V05.P_lin 3 (n + 1) V05.z2 V05.z3, V05.L3_three (by omega), ← V05.P3_x2]
  push_cast
  ring
