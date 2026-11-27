import Mathlib
import Zeta3.Wave8.Q04_FlatAllK
/-!
# S03 -- the field rule: `M_k(α)` as a finite combination of harmonic polylogarithms in `c`

Target id: S03.
Paper label: `conj:field` (paper/ieee/main.tex, Theorem "Field rule", lines 1486-1494; the
"More precisely" sentence, line 1493).
Supplement: paper/ieee/proofs/sec_g6_field.tex (harmonic polylogarithms, line 9; theorem
`g6:thm:field` 19-26 and the reading of its last sentence, line 28; nested sums 35-39;
Lemma `g6:lem:stuffle` 41-51; Lemma `g6:lem:window` 53-63; the space `𝔇` and its filtration
67-71; Lemma `g6:lem:A` 73-102; Lemma `g6:lem:dict` 106-116 and `g6:eq:block` 118-121;
Lemma `g6:lem:ops` 123-129; Lemma `g6:prop:R` 133-137 with the window functions 139-142;
Lemma `g6:lem:conv` 144-151; proof of the last sentence 184-188).

Imported, not restated (all from frozen statement blocks): `M_k(α)` is `V05.Mk k α`
(`= iteratedDeriv k (D04.Mint · α) 0 / k!`), `c = D04.c α = cos(α/2)`, the real zeta value
`ζ(n) = X05.z n = (riemannZeta n).re` (`riemannZeta n` is real for real `n > 1`, and
`zeta_nat_eq_tsum_of_gt_one` identifies it with `Σ_{m≥1} m^{-n}`), and the formal coefficient
`P_{k,n} ∈ ℚ[𝔷_2, 𝔷_3, …]` of Definition `g6:def:formal`, `Q04.Pf k n` (Q04 proved
`M_k(α) = Σ_{n≥1} P_{k,n}(ζ) c^{2n}`, Lemma `g6:prop:formula`, as `Q04.stmt_5`). The iterated
integrals of R05 (`R05.iterInt`, Goncharov integrals over `1 > t_1 > ⋯ > t_w > 0`) do not fit:
the harmonic polylogarithms here are functions of the upper end point `x`, so they are defined
again below.

Quoted verbatim (main.tex 1486-1494, `conj:field`):

  \begin{theorem}[Field rule$^\dagger$]\label{conj:field}
  Let $0<\alpha<\pi$, $c=\cos\frac\alpha2$, $s=\sin\frac\alpha2$ and $K=\mathbb Q(c,s)$.
  \begin{enumerate}
  \item[(a)] If $c$ is algebraic, $K$ is a number field. Let $S$ be the set of its infinite places
  and of the primes that divide $2$, $c$, $1-c$ or $1+c$. Then every $M_k(\alpha)$ is a
  combination of polylogarithms at $S$-units of $K$.
  \item[(b)] If $c$ is transcendental, $K$ is a function field, not a number field. Let $S$ be the
  set of places of $K$ over $c=0,1,-1,\infty$. Then, as a function of $c$, every $M_k(\alpha)$ is
  a combination of polylogarithms at $S$-units of $K$, up to nonzero rational constants.
  \end{enumerate}
  More precisely, for every $k$, $M_k(\alpha)=\sum_jz_j\,r_j(c)\,H_{w_j}(c)$, a finite sum, where
  the $z_j$ are monomials in $\zeta(2),\dots,\zeta(k)$, the $r_j$ are rational functions of $c$
  whose denominators divide products of powers of $c$, $1-c$ and $1+c$, the $H_{w_j}$ are harmonic
  polylogarithms with letters $0,1,-1$, and $\operatorname{wt}(z_j)+|w_j|\le k+1$.
  \end{theorem}

Quoted verbatim (supplement, line 9): "Let $f_0(t)=1/t$, $f_1(t)=1/(1-t)$, $f_{-1}(t)=1/(1+t)$. A
\emph{word} is a finite sequence of letters in $\{0,1,-1\}$; its length is $|w|$ (the weight). It
is \emph{regular} if it is empty or its last letter is nonzero. For regular words put
$H_\emptyset=1$ and $H_{(a,w')}(x)=\int_0^xf_a(t)H_{w'}(t)\,dt$ (the integrand is integrable at
$0$: $H_{w'}(0)=0$ for $w'\ne\emptyset$)."

Quoted verbatim (supplement, line 28): "We prove it with
$r_j\in\QQ[c,c^{-1},(1-c)^{-1},(1+c)^{-1}]$ and regular words $w_j$."

Quoted verbatim (supplement, lines 35-37, nested sums): "An \emph{index word} is a finite
sequence $v=(a_1,\dots,a_d)$ of nonzero integers, with weight $|v|=\sum_r|a_r|$ ... Put
$\sigma_r=\gsgn(a_r)$. For an integer $M\ge0$,
  Z_v(M)=\sum_{M\ge i_1>\dots>i_d\ge1}\ \prod_{r=1}^d\frac{\sigma_r^{\,i_r}}{i_r^{|a_r|}}"

Quoted verbatim (`g6:lem:stuffle`, line 42): "For index words $u,v$ there is a finite multiset
$u*v$ of index words, each of weight $|u|+|v|$, such that $Z_u(M)Z_v(M)=\sum_{w\in u*v}Z_w(M)$ for
all $M\ge0$."

Quoted verbatim (`g6:lem:window`, lines 54-58): "Let $v=(a_1,\dots,a_d)$ be an index word and
$0\le M_0\le M$ integers. Then
  Z_v(M)=\sum_{t=0}^dW_t\,Z_{v_{>t}}(M_0),\qquad W_t=\sum_{M\ge i_1>\dots>i_t>M_0}\ \prod_{r\le t}
  \frac{\sigma_r^{\,i_r}}{i_r^{|a_r|}},\quad W_0=1,
with $v_{>t}=(a_{t+1},\dots,a_d)$."

Quoted verbatim (line 67): "Let $\mathfrak D$ be the $\QQ$-vector space of finite formal sums
$T=\sum_wR_w[w]$, where $w$ runs over index words and $R_w\in\QQ(N)$ has poles only at integers
$p\le1$. ... We say that $T\in\mathfrak D^{(B)}$ if $|w|+\gord_p(R_w)\le B$ for every term
$R_w[w]$ of $T$ and every $p\in\ZZ$." (line 71) "put $\gsub{T}{n}:=\sum_wR_w(2n)Z_w(2n)$".

Quoted verbatim (`g6:lem:A`, line 74): "Let $k\ge0$ and let $\gvs z$ be a monomial in
$\gvs2,\dots,\gvs k$ of weight $W\le k$. There is $T_z\in\mathfrak D^{(k+1-W)}$ such that
$[\gvs z]P_{k,n}=\gsub{T_z}{n}$ for every $n\ge1$."
(proof, line 88) "\gsub{\chi^{(r)}_1}{n}=H^{(r)}_{n+j-1},\qquad\gsub{\chi^{(r)}_2}{n}=H^{(r)}_{n-1},
\qquad\gsub{\chi^{(r)}_3}{n}=H^{(r)}_{2n+j-1}" with (line 90) "$Z_r(2n)+Z_{-r}(2n)=\dots=
2^{1-r}H^{(r)}_n$ ... Next $Z_r(2n)-(2n)^{-r}=H^{(r)}_{2n-1}$."

Quoted verbatim (`g6:lem:dict`, lines 107-111): "Let $w$ be a nonempty regular word with nonzero
letters $b_1,\dots,b_\nu$ (left to right), $m_0$ zeros before $b_1$ and $m_t$ zeros between $b_t$
and $b_{t+1}$ $(1\le t<\nu)$. Put $\varepsilon=b_1$, $\kappa=1+m_0$, $\gcf_w=b_1\cdots b_\nu$, and
let $v_w$ be the index word $\bigl(b_1b_2(1+m_1),\dots,b_{\nu-1}b_\nu(1+m_{\nu-1})\bigr)$. Then for
every $N\ge1$ [x^N]H_w(x)=\gcf_w\,\varepsilon^N\,N^{-\kappa}\,Z_{v_w}(N-1),\qquad |w|=\kappa+|v_w| .
The map $w\mapsto(\varepsilon,\kappa,v_w)$ is a bijection from the nonempty regular words onto
$\{\pm1\}\times\ZZ_{\ge1}\times\{\text{index words}\}$." (line 9) "By induction each $H_w$,
$w\ne\emptyset$, is a power series in $x$ without constant term, convergent for $|x|<1$."

Quoted verbatim (`g6:lem:ops`, line 124): "(a) For a nonempty regular word $w$,
$H_w(-x)=(-1)^{\nu}H_{\bar w}(x)$, where $\nu$ is the number of nonzero letters of $w$ and $\bar w$
is $w$ with every nonzero letter replaced by its negative. (b) Let $\theta=x\frac d{dx}$. Then
$\theta H_{(0,w')}=H_{w'}$ and $\theta H_{(\pm1,w')}=\frac x{1\mp x}H_{w'}$. (c) Let
$\mathcal K=\QQ[x,x^{-1},(1-x)^{-1},(1+x)^{-1}]$ and $\mathcal G_B$ the $\QQ$-span of the products
$r(x)H_w(x)$, $r\in\mathcal K$, $w$ regular with $|w|\le B$ (with $H_\emptyset=1$). Then
$\mathcal G_B\supset\mathcal G_{B'}$ for $B'\le B$, $\mathcal G_B$ is stable under $x\mapsto-x$ and
under $\theta$."

Quoted verbatim (`g6:prop:R`, line 136; line 133): "Let $B\ge0$ and $T\in\mathfrak D^{(B)}$. There is
$F\in\mathcal G_B$ such that $[x^N]F=\mathcal S_T(N)$ for every even $N\ge2$." "All series in $x$
below are Laurent series at $x=0$. For $T=\sum_wR_w[w]\in\mathfrak D$ and $N\ge2$ put
$\mathcal S_T(N)=\sum_wR_w(N)Z_w(N)$".

Quoted verbatim (window functions, line 141, and `g6:lem:conv`, lines 145-146):
  W^{(\tau)}_t(v;N)=\sum_{\tau-2\ge\delta_1>\dots>\delta_t\ge-1}\ \prod_{r=1}^t\sigma_r^{\delta_r}\,
  (N+\delta_r)^{-|a_r|}\in\QQ(N),\qquad W^{(\tau)}_0=1 .
"(a) For every even $N\ge2$, every index word $v$ and $\tau\ge0$:
$Z_v(N-2+\tau)=\sum_{t=0}^dW^{(\tau)}_t(v;N)\,Z_{v_{>t}}(N-2)$. (b) The poles of $W^{(\tau)}_t$ are
at integers $\le1$ (at $N=-\delta$, $-1\le\delta\le\tau-2$, if $\tau\ge1$; for $\tau=0$ only $t=0$
occurs), and $\gord_p(W^{(\tau)}_t(v;N))\le|v|-|v_{>t}|$ for every $p\in\ZZ$."

Quoted verbatim (proof of the theorem, line 185): "By Lemma~\ref{g6:prop:R} there is
$F_0\in\mathcal G_B$ with $[x^N]F_0=E_z(N/2)$ for even $N\ge2$. ... Put $F_z=F_1-L\in\mathcal G_B$.
Then $F_z=\sum_{n\ge1}E_z(n)x^{2n}$ as Laurent series. ... $F_z(x)=\sum_nE_z(n)x^{2n}$ for
$0<|x|<1$." (line 187) "This is $M_k(\alpha)=\sum_jz_jr_j(c)H_{w_j}(c)$ with $z_j=\zeta^z$ and
$\gwt(z_j)+|w_j|\le W+B=k+1$."

Quoted verbatim (main.tex 1472): "M_0(\alpha)=-(1+c)\log(1+c)-(1-c)\log(1-c)."

## Formalization choices

* Letters and index letters are integers (`List ℤ`); a word over `{0,1,-1}` satisfies
  `IsLetters`, an index word satisfies `IsIndex` (nonzero letters). Regularity (`Regular`) is
  "empty or last letter nonzero". `|w|` is `w.length` for a word and `wt v = Σ |a_r|` for an
  index word.
* `H w x` is the iterated integral of the supplement, by recursion on the word, with the
  interval integral `∫_0^x` (so it is defined for `|x| < 1`, negative `x` included, which
  `g6:lem:ops`(a) needs). It is only ever used for regular words over `{0,1,-1}`; for those the
  integrands are integrable (that is stated, in `stmt_5`), so the interval integral is never a
  junk value. Trailing zeros (non-regular words, powers of `log x`) never occur, as in the
  supplement (line 28: regular words suffice).
* `Z v M` is the nested sum, by the recursion of line 39; `Zw v M₀ M` is the window sum `W_t`
  (indices in `(M₀, M]`), so `W_t = Zw (v.take t) M₀ M`, `v_{>t} = v.drop t`.
* An element `Σ R_w [w]` of `𝔇` is a list of triples `(Q, e, w)` standing for
  `R_w(N) = Q(N) / Π_p (N - p)^{e_p}` with `Q ∈ ℚ[N]` and a finitely supported exponent
  `e : ℤ →₀ ℕ` (`ratEval`). "Poles only at integers `p ≤ 1`, `|w| + ord_p R_w ≤ B`" is `InDB B`:
  every `p` in the support of `e` has `p ≤ 1` and `|w| + e_p ≤ B`, and `|w| ≤ B`. Every rational
  function with poles at integers `≤ 1` of orders `o_p` has such a form with `e_p = o_p`, and
  conversely the form bounds the pole orders by `e_p`, so this is the paper's filtration.
  `ST T N = Σ R_w(N) Z_w(N)` is `𝒮_T(N)`, and `⟨T⟩_n = ST T (2n)`.
* `𝒦 = ℚ[x, x⁻¹, (1-x)⁻¹, (1+x)⁻¹]` is written as `p(x) / (x^a (1-x)^b (1+x)^d)`,
  `p ∈ ℚ[x]` (`kval`), and `𝒢_B` (`InG B F`) is the set of real functions `F` that agree on
  `0 < |x| < 1` with a finite sum `Σ r_i(x) H_{w_i}(x)`, `r_i ∈ 𝒦`, `w_i` regular over
  `{0,1,-1}`, `|w_i| ≤ B`. The statement "`[x^N]F = …`" for `F ∈ 𝒢_B` (a Laurent series with
  finitely many negative powers) is written as: there are a polynomial `q` and rationals `a_N`
  with `F(x) = q(1/x) + Σ_{N≥0} a_N x^N` on `0 < |x| < 1`, and `a_N = …`. (Laurent coefficients
  are unique except that the constant term may be shared between `q` and `a_0`; only `N ≥ 2` is
  constrained.)
* The window functions `W^{(τ)}_t(v;N)` are `Wc (v.take t) (τ - 1) N`, a nested sum over
  `-1 ≤ δ_t < ⋯ < δ_1 < τ - 1` with `σ^{δ}` an integer power (`σ = ±1`, so `σ^{-1} = σ`). The
  rational-function claim of `g6:lem:conv`(b) is stated as: `W` agrees with some `Q/Π(N-p)^{e_p}`
  (poles at `p ≤ 1`, `e_p ≤ |v| - |v_{>t}| = wt (v.take t)`) at every rational `N ≥ 2`.
* The dictionary data `(ε, κ, 𝔠_w, v_w)` of `g6:lem:dict` are computed by `dd` with the
  recursion of the lemma's proof (prepending `0` raises `κ`; prepending `a = ±1` to `w'` with data
  `(ε', κ', 𝔠', v')` gives `(a, 1, a𝔠', (aε'κ', v'))`), which is the lemma's explicit
  description `ε = b_1`, `κ = 1 + m_0`, `𝔠_w = b_1⋯b_ν`, `v_w = (b_1b_2(1+m_1), …)`; this was
  checked by machine on all regular words of length `≤ 8` against the explicit description.
* Zeta monomials are `s : ℕ →₀ ℕ` (`𝔷^s = Π 𝔷_w^{s_w}`) with all `w ∈ supp s` in `[2, k]`; the
  weight is `mwt s = Σ s_w w`; the real value is `zmon s = Π ζ(w)^{s_w}`.
* The coefficient `[𝔷^s] P_{k,n}` is `(Q04.Pf k n).coeff s`.
* The final statement (`stmt_11`) is the "More precisely" sentence: for every `k` one finite list
  of `(s_j, p_j, a_j, b_j, d_j, w_j)`, chosen before `α`, with
  `M_k(α) = Σ_j ζ^{s_j} · p_j(c)/(c^{a_j}(1-c)^{b_j}(1+c)^{d_j}) · H_{w_j}(c)` for every
  `0 < α < π`, `ζ^{s_j}` a monomial in `ζ(2),…,ζ(k)`, `w_j` regular over `{0,1,-1}`, and
  `wt(ζ^{s_j}) + |w_j| ≤ k + 1`.

## Statements

* stmt_1: Lemma `g6:lem:stuffle` (stuffle product of nested sums).
* stmt_2: Lemma `g6:lem:window` (windows).
* stmt_3: the harmonic-number identities in the proof of `g6:lem:A`:
  `2^{r-1}(Z_r(2n) + Z_{-r}(2n)) - n^{-r} = H^{(r)}_{n-1}` and `Z_r(2n) - (2n)^{-r} = H^{(r)}_{2n-1}`.
* stmt_4: Lemma `g6:lem:A`, for every `k` and every zeta monomial of weight `W ≤ k`.
* stmt_5: the harmonic polylogarithms are well defined (integrable integrands) and Lemma
  `g6:lem:dict`: power series on `|x| < 1` with `[x^N]H_w = 𝔠_w ε^N N^{-κ} Z_{v_w}(N-1)`,
  `|w| = κ + |v_w|`; and the map `w ↦ (ε, κ, v_w)` is a bijection onto `{±1} × ℤ_{≥1} × {index
  words}`.
* stmt_6: Lemma `g6:lem:ops` (a), (b), and (c) (stability of `𝒢_B` under `x ↦ -x` and `θ`;
  `θF = G` is stated as `F' = G/x`).
* stmt_7: Lemma `g6:lem:conv` (a) and (b).
* stmt_8: Lemma `g6:prop:R`.
* stmt_9: the first paragraph of the proof of the theorem: for every `k` and zeta monomial `𝔷^s`
  of weight `W ≤ k` there is `F_s ∈ 𝒢_{k+1-W}` with `F_s(x) = Σ_{n≥1} [𝔷^s]P_{k,n} x^{2n}` on
  `0 < |x| < 1`.
* stmt_10: the case `k = 0` explicitly: `M_0(α) = -(1+c) H_{-1}(c) + (1-c) H_1(c)` (main.tex 1472,
  with `H_{-1} = log(1+x)`, `H_1 = -log(1-x)`); it is independent of the general construction.
* stmt_11: the "More precisely" sentence of `conj:field` (final claim).

## Numerical check before freezing

A Python script (exact rationals; floats for quadrature) mirroring these definitions checked:
`stmt_1` with the explicit quasi-shuffle of the proof, all pairs of index words of total weight
`≤ 5`, `M ≤ 8`; `stmt_2` for all index words of weight `≤ 5`, `M₀ ≤ M ≤ 9`; `stmt_7`(a) for the
same words, `τ ≤ 4`, even `N ≤ 12`; `stmt_3` for `r ≤ 5`, `n ≤ 10`; `dd` against the lemma's
explicit data `(b_1, 1+m_0, Π b, (b_tb_{t+1}(1+m_t)))` and `|w| = κ + |v_w|` for all 2186 regular
words of length `≤ 7`, injectivity of `w ↦ (ε, κ, v_w)` on them and that their images exhaust
the 2186 triples with `κ + |v| ≤ 7`; the coefficient formula of `stmt_5` against the recursion of
the coefficients of `H_w` (proof, line 115), regular words of length `≤ 6`, `N ≤ 30`;
`stmt_6`(a) on the coefficients, length `≤ 6`; the series of `stmt_5` against the iterated
integral (Simpson) at `x = 0.37` and `x = -0.6`, length `≤ 3`.
Final statement: for `k = 0` and `k = 1` an exact fit of `P_{k,n}` (`P_{0,n} = -1/(n(2n-1))`,
`P_{1,n} = -(H_{n-1} - H_{2n-1} + (2n-1)/(2n))/(n(2n-1))`, worked out from Definition
`g6:def:formal`) by `Σ x_{e,w} c^e H_w(c)`, `e ∈ {0,1}`, regular `|w| ≤ k+1`, has full rank and is
exact on all coefficients `N ≤ 140`. It gives `M_0 = (1-c)H_1 - (1+c)H_{-1}` (this is `stmt_10`,
the closed form of main.tex 1472) and
`M_1 = -H_{0,1} + H_{0,-1} - H_{1,-1} - H_{-1,1} + c(H_{0,1} + H_{0,-1} + H_{1,-1} - H_{-1,1})`;
at `c = 1/2` this is `0.125795678937973`, equal (to `3·10⁻¹⁶`) to the closed form
`24M_1(2π/3) = π² + 48λ² - 6log²3 - 24λlog3 - 12Li₂(1/3)` of main.tex, and the two series agree
at `c = 0.1, 0.5, 0.9`. For `k = 2, 3` (with zeta monomials `1, 𝔷_2, 𝔷_3`) the supplement's own exact
fit `g6_field_check.out` (read, not re-run: SymPy and python-flint are not installed here) reports
exact fits with `e ∈ [-(k-1), 1]`, `|w| ≤ k + 1 - W`, i.e. the shape of `stmt_11` with
`a_j = k - 1`, `b_j = d_j = 0`, and agreement with `M_k(2π/3)` from the definition to `4·10⁻²⁸`
(`k = 2`) and `4·10⁻²⁷` (`k = 3`).

## NOT COVERED

* Parts (a) and (b) of `conj:field` (combinations of multiple polylogarithms at `S`-units of the
  number field / function field `K = ℚ(c,s)`): `S`-units, places of `K` and multiple
  polylogarithms `Li_{s_1,…,s_d}` are not formalized. (Their proof in the supplement is the
  identity `H_w(x) = 𝔠_w Li_{κ,|a_1|,…}(εx, sgn a_1, …)`, which is `stmt_5` written out.)
* The substitution step of the proof of `g6:lem:A` with `Λ_j` and the `χ^{(r)}_•` as elements of
  `𝔇` (only the numerical identities behind `⟨χ⟩_n = H^{(r)}_•` are stated, in `stmt_3`); the
  multiplicativity of `T ↦ ⟨T⟩_n` is the content of `stmt_1`.
* The inclusion `𝒢_{B'} ⊂ 𝒢_B` of `g6:lem:ops`(c) (trivial).
* Remark `g6:rem:k7` (denominators are powers of `c` alone for `k ≤ 7`; pole set in `n`): not
  claimed by the theorem, not proved in the supplement.

Statements 4, 8, 9, 11 moved to `Zeta3.Wave10.S08_FieldRuleMain` (continuation of S03); this file keeps 1, 2, 3, 5, 6, 7, 10.
-/
-- BEGIN STATEMENTS
namespace S03

open Real MeasureTheory Set Filter

/-- The kernels `f_0(t) = 1/t`, `f_1(t) = 1/(1-t)`, `f_{-1}(t) = 1/(1+t)` (only used for
letters `0, 1, -1`; any other integer gets `f_{-1}`). -/
noncomputable def f (a : ℤ) (t : ℝ) : ℝ :=
  if a = 0 then 1 / t else if a = 1 then 1 / (1 - t) else 1 / (1 + t)

/-- A word over the letters `0, 1, -1`. -/
def IsLetters (w : List ℤ) : Prop := ∀ a ∈ w, a = 0 ∨ a = 1 ∨ a = -1

/-- A regular word: empty, or its last letter is nonzero. -/
def Regular (w : List ℤ) : Prop := ∀ a ∈ w.getLast?, a ≠ 0

/-- Harmonic polylogarithms: `H_∅ = 1`, `H_{(a,w')}(x) = ∫_0^x f_a(t) H_{w'}(t) dt`. -/
noncomputable def H : List ℤ → ℝ → ℝ
  | [], _ => 1
  | a :: w, x => ∫ t in (0 : ℝ)..x, f a t * H w t

/-- An index word: a finite sequence of nonzero integers. -/
def IsIndex (v : List ℤ) : Prop := ∀ a ∈ v, a ≠ 0

/-- Weight `|v| = Σ_r |a_r|` of an index word. -/
def wt (v : List ℤ) : ℕ := (v.map Int.natAbs).sum

/-- Nested sums `Z_v(M) = Σ_{M ≥ i_1 > ⋯ > i_d ≥ 1} Π σ_r^{i_r} / i_r^{|a_r|}`, by
`Z_{(a,v')}(M) = Σ_{i=1}^M σ^i i^{-|a|} Z_{v'}(i-1)`, `Z_∅ = 1`. -/
def Z : List ℤ → ℕ → ℚ
  | [], _ => 1
  | a :: v, M => ∑ i ∈ Finset.Icc 1 M,
      (Int.sign a : ℚ) ^ i / (i : ℚ) ^ a.natAbs * Z v (i - 1)

/-- Window sums `Σ_{M ≥ i_1 > ⋯ > i_t > M₀} Π σ_r^{i_r} / i_r^{|a_r|}` (the `W_t` of
`g6:lem:window` for the word `(a_1, …, a_t)`). -/
def Zw : List ℤ → ℕ → ℕ → ℚ
  | [], _, _ => 1
  | a :: v, M₀, M => ∑ i ∈ Finset.Ioc M₀ M,
      (Int.sign a : ℚ) ^ i / (i : ℚ) ^ a.natAbs * Zw v M₀ (i - 1)

/-- Generalized harmonic number `H^{(r)}_N = Σ_{l=1}^N l^{-r}` (rational). -/
def hn (r N : ℕ) : ℚ := ∑ l ∈ Finset.Icc 1 N, 1 / (l : ℚ) ^ r

/-- The rational function `Q(N) / Π_p (N - p)^{e_p}`. -/
noncomputable def ratEval (Q : Polynomial ℚ) (e : ℤ →₀ ℕ) (N : ℚ) : ℚ :=
  Q.eval N / e.prod (fun p k => (N - (p : ℚ)) ^ k)

/-- A term `R_w [w]` of `𝔇`: `(Q, e, w)` with `R_w = Q / Π (N-p)^{e_p}`. -/
abbrev DTerm := Polynomial ℚ × (ℤ →₀ ℕ) × List ℤ

/-- `T ∈ 𝔇^{(B)}`: index words, poles only at integers `p ≤ 1`, `|w| + ord_p R_w ≤ B`. -/
def InDB (B : ℕ) (T : List DTerm) : Prop :=
  ∀ τ ∈ T, IsIndex τ.2.2 ∧ wt τ.2.2 ≤ B ∧
    ∀ p ∈ τ.2.1.support, p ≤ 1 ∧ wt τ.2.2 + τ.2.1 p ≤ B

/-- `𝒮_T(N) = Σ_w R_w(N) Z_w(N)` (so `⟨T⟩_n = 𝒮_T(2n)`). -/
noncomputable def ST (T : List DTerm) (N : ℕ) : ℚ :=
  (T.map fun τ => ratEval τ.1 τ.2.1 (N : ℚ) * Z τ.2.2 N).sum

/-- The data `(ε, κ, 𝔠_w, v_w)` of `g6:lem:dict`, by the recursion of its proof. -/
def dd : List ℤ → ℤ × ℕ × ℤ × List ℤ
  | [] => (1, 0, 1, [])
  | [b] => (b, 1, b, [])
  | a :: b :: w =>
      let d := dd (b :: w)
      if a = 0 then (d.1, d.2.1 + 1, d.2.2.1, d.2.2.2)
      else (a, 1, a * d.2.2.1, (a * d.1 * (d.2.1 : ℤ)) :: d.2.2.2)

/-- An element `r(x) = p(x) / (x^a (1-x)^b (1+x)^d)` of `𝒦`, `p ∈ ℚ[x]`. -/
noncomputable def kval (p : Polynomial ℚ) (a b d : ℕ) (x : ℝ) : ℝ :=
  Polynomial.aeval x p / (x ^ a * (1 - x) ^ b * (1 + x) ^ d)

/-- A term `r(x) H_w(x)` of `𝒢_B`: `(p, a, b, d, w)`. -/
abbrev GTerm := Polynomial ℚ × ℕ × ℕ × ℕ × List ℤ

/-- `F ∈ 𝒢_B`: on `0 < |x| < 1`, `F = Σ_i r_i H_{w_i}` with `r_i ∈ 𝒦` and `w_i` regular over
`{0,1,-1}` with `|w_i| ≤ B`. -/
def InG (B : ℕ) (F : ℝ → ℝ) : Prop :=
  ∃ L : List GTerm,
    (∀ g ∈ L, IsLetters g.2.2.2.2 ∧ Regular g.2.2.2.2 ∧ g.2.2.2.2.length ≤ B) ∧
    ∀ x : ℝ, 0 < |x| → |x| < 1 →
      F x = (L.map fun g => kval g.1 g.2.1 g.2.2.1 g.2.2.2.1 x * H g.2.2.2.2 x).sum

/-- The window functions: `Wc v U N = Σ_{U > δ_1 > ⋯ > δ_t ≥ -1} Π σ_r^{δ_r} (N+δ_r)^{-|a_r|}`
for `v = (a_1, …, a_t)`; `W^{(τ)}_t(v;N) = Wc (v.take t) (τ - 1) N`. -/
def Wc : List ℤ → ℤ → ℚ → ℚ
  | [], _, _ => 1
  | a :: v, U, N => ∑ δ ∈ Finset.Ico (-1 : ℤ) U,
      (Int.sign a : ℚ) ^ δ / (N + (δ : ℚ)) ^ a.natAbs * Wc v δ N

/-- Weight of a zeta monomial `𝔷^s = Π 𝔷_w^{s_w}`: `Σ s_w w`. -/
def mwt (s : ℕ →₀ ℕ) : ℕ := s.sum fun w e => e * w

/-- The real value `ζ^s = Π ζ(w)^{s_w}`, `ζ(w) = (riemannZeta w).re`. -/
noncomputable def zmon (s : ℕ →₀ ℕ) : ℝ := s.prod fun w e => X05.z w ^ e

/-- A term `ζ^{s} r(c) H_w(c)` of the final formula: `(s, p, a, b, d, w)`. -/
abbrev FTerm := (ℕ →₀ ℕ) × Polynomial ℚ × ℕ × ℕ × ℕ × List ℤ

/-- Lemma `g6:lem:stuffle`: for index words `u, v` there is a finite multiset of index words, each
of weight `|u| + |v|`, with `Z_u(M) Z_v(M) = Σ_{w ∈ u*v} Z_w(M)` for all `M ≥ 0`. -/
def stmt_1 : Prop :=
  ∀ u v : List ℤ, IsIndex u → IsIndex v →
    ∃ S : Multiset (List ℤ), (∀ w ∈ S, IsIndex w ∧ wt w = wt u + wt v) ∧
      ∀ M : ℕ, Z u M * Z v M = (S.map fun w => Z w M).sum

/-- Lemma `g6:lem:window`: for an index word `v = (a_1,…,a_d)` and `0 ≤ M₀ ≤ M`,
`Z_v(M) = Σ_{t=0}^d W_t Z_{v_{>t}}(M₀)`. -/
def stmt_2 : Prop :=
  ∀ v : List ℤ, IsIndex v → ∀ M₀ M : ℕ, M₀ ≤ M →
    Z v M = ∑ t ∈ Finset.range (v.length + 1), Zw (v.take t) M₀ M * Z (v.drop t) M₀

/-- Proof of `g6:lem:A`, the harmonic-number identities: for `r, n ≥ 1`,
`2^{r-1}(Z_r(2n) + Z_{-r}(2n)) - n^{-r} = H^{(r)}_{n-1}` and `Z_r(2n) - (2n)^{-r} = H^{(r)}_{2n-1}`. -/
def stmt_3 : Prop :=
  ∀ r n : ℕ, 1 ≤ r → 1 ≤ n →
    (2 : ℚ) ^ (r - 1) * (Z [(r : ℤ)] (2 * n) + Z [-(r : ℤ)] (2 * n)) - 1 / (n : ℚ) ^ r
        = hn r (n - 1) ∧
      Z [(r : ℤ)] (2 * n) - 1 / (2 * (n : ℚ)) ^ r = hn r (2 * n - 1)

/-- Harmonic polylogarithms are well defined, and Lemma `g6:lem:dict`. For a word `(a, w')` over
`{0,1,-1}` that is regular, the integrand `f_a H_{w'}` is interval integrable on `[0, x]`,
`|x| < 1`. For a nonempty regular word `w` with data `(ε, κ, 𝔠_w, v_w) = dd w`: `ε = ±1`,
`κ ≥ 1`, `𝔠_w = ±1`, `v_w` is an index word, `|w| = κ + |v_w|`, and for `|x| < 1`,
`H_w(x) = Σ_{N≥1} 𝔠_w ε^N N^{-κ} Z_{v_w}(N-1) x^N`. The map `w ↦ (ε, κ, v_w)` is a bijection from
the nonempty regular words onto `{±1} × ℤ_{≥1} × {index words}`. -/
def stmt_5 : Prop :=
  (∀ a : ℤ, ∀ w : List ℤ, IsLetters (a :: w) → Regular (a :: w) → ∀ x : ℝ, |x| < 1 →
    IntervalIntegrable (fun t => f a t * H w t) volume 0 x) ∧
  (∀ w : List ℤ, w ≠ [] → IsLetters w → Regular w →
    ((dd w).1 = 1 ∨ (dd w).1 = -1) ∧ 1 ≤ (dd w).2.1 ∧
    ((dd w).2.2.1 = 1 ∨ (dd w).2.2.1 = -1) ∧ IsIndex (dd w).2.2.2 ∧
    w.length = (dd w).2.1 + wt (dd w).2.2.2 ∧
    ∀ x : ℝ, |x| < 1 →
      HasSum (fun N : ℕ => ((dd w).2.2.1 : ℝ) * ((dd w).1 : ℝ) ^ (N + 1) /
          ((N : ℝ) + 1) ^ (dd w).2.1 * (Z (dd w).2.2.2 N : ℝ) * x ^ (N + 1))
        (H w x)) ∧
  (∀ ε : ℤ, ∀ κ : ℕ, ∀ v : List ℤ, (ε = 1 ∨ ε = -1) → 1 ≤ κ → IsIndex v →
    ∃! w : List ℤ, w ≠ [] ∧ IsLetters w ∧ Regular w ∧
      (dd w).1 = ε ∧ (dd w).2.1 = κ ∧ (dd w).2.2.2 = v)

/-- Lemma `g6:lem:ops`. (a) For a nonempty regular word `w` over `{0,1,-1}` and `|x| < 1`,
`H_w(-x) = (-1)^ν H_{w̄}(x)`, `ν` the number of nonzero letters, `w̄ = -w`. (b) For `0 < |x| < 1`
and a regular word `(a, w')`: `H_{(0,w')}' = H_{w'}/x`, `H_{(1,w')}' = H_{w'}/(1-x)`,
`H_{(-1,w')}' = H_{w'}/(1+x)` (that is `θH_{(0,w')} = H_{w'}`, `θH_{(±1,w')} = x/(1∓x) H_{w'}`).
(c) `𝒢_B` is stable under `x ↦ -x` and under `θ = x d/dx`. -/
def stmt_6 : Prop :=
  (∀ w : List ℤ, w ≠ [] → IsLetters w → Regular w → ∀ x : ℝ, |x| < 1 →
    H w (-x) = (-1 : ℝ) ^ (w.filter fun a => a ≠ 0).length * H (w.map fun a => -a) x) ∧
  (∀ a : ℤ, ∀ w : List ℤ, IsLetters (a :: w) → Regular (a :: w) → ∀ x : ℝ, 0 < |x| → |x| < 1 →
    HasDerivAt (H (a :: w)) (f a x * H w x) x) ∧
  (∀ B : ℕ, ∀ F : ℝ → ℝ, InG B F →
    InG B (fun x => F (-x)) ∧
    ∃ G : ℝ → ℝ, InG B G ∧ ∀ x : ℝ, 0 < |x| → |x| < 1 → HasDerivAt F (G x / x) x)

/-- Lemma `g6:lem:conv`. (a) For even `N ≥ 2`, an index word `v` and `τ ≥ 0`:
`Z_v(N-2+τ) = Σ_{t=0}^d W^{(τ)}_t(v;N) Z_{v_{>t}}(N-2)`. (b) `W^{(τ)}_t(v;·)` is a rational
function with poles only at integers `p ≤ 1`, of order `≤ |v| - |v_{>t}|`: it agrees with
`Q/Π(N-p)^{e_p}`, `p ≤ 1`, `e_p ≤ wt(a_1,…,a_t)`, at every rational `N ≥ 2`. -/
def stmt_7 : Prop :=
  (∀ v : List ℤ, IsIndex v → ∀ τ N : ℕ, 2 ≤ N → Even N →
    Z v (N - 2 + τ) = ∑ t ∈ Finset.range (v.length + 1),
      Wc (v.take t) ((τ : ℤ) - 1) (N : ℚ) * Z (v.drop t) (N - 2)) ∧
  (∀ v : List ℤ, IsIndex v → ∀ τ t : ℕ, ∃ Q : Polynomial ℚ, ∃ e : ℤ →₀ ℕ,
    (∀ p ∈ e.support, p ≤ 1 ∧ e p ≤ wt (v.take t)) ∧
    ∀ N : ℚ, 2 ≤ N → Wc (v.take t) ((τ : ℤ) - 1) N = ratEval Q e N)

/-- The case `k = 0` written out: `M_0(α) = -(1+c) H_{-1}(c) + (1-c) H_1(c)` for `0 < α < π`
(main.tex 1472, `H_{-1}(x) = log(1+x)`, `H_1(x) = -log(1-x)`). -/
def stmt_10 : Prop :=
  ∀ α : ℝ, 0 < α → α < π →
    V05.Mk 0 α = -(1 + D04.c α) * H [-1] (D04.c α) + (1 - D04.c α) * H [1] (D04.c α)

end S03
-- END STATEMENTS


namespace S03

lemma Z_single (a : ℤ) (M : ℕ) :
    Z [a] M = ∑ i ∈ Finset.Icc 1 M, (Int.sign a : ℚ) ^ i / (i : ℚ) ^ a.natAbs := by
  simp [Z]

lemma hn_succ (r N : ℕ) : hn r (N + 1) = hn r N + 1 / ((N : ℚ) + 1) ^ r := by
  unfold hn
  rw [Finset.sum_Icc_succ_top (by omega)]
  push_cast; ring

lemma stmt3_aux (r : ℕ) (hr : 1 ≤ r) (n : ℕ) :
    (2 : ℚ) ^ (r - 1) * (Z [(r : ℤ)] (2 * n) + Z [-(r : ℤ)] (2 * n)) = hn r n := by
  have h1 : Int.sign (r : ℤ) = 1 := Int.sign_eq_one_of_pos (by omega)
  have h2 : Int.sign (-(r : ℤ)) = -1 := by
    rw [Int.sign_neg, h1]
  have h3 : (-(r : ℤ)).natAbs = r := by simp
  rw [Z_single, Z_single, h1, h2, h3]
  simp only [Int.natAbs_natCast, Int.cast_one, one_pow, Int.cast_neg]
  induction n with
  | zero => simp [hn]
  | succ n ih =>
    rw [hn_succ, ← ih, show 2 * (n + 1) = 2 * n + 1 + 1 by ring]
    rw [Finset.sum_Icc_succ_top (by omega), Finset.sum_Icc_succ_top (by omega),
      Finset.sum_Icc_succ_top (by omega), Finset.sum_Icc_succ_top (by omega)]
    have hp : (-1 : ℚ) ^ (2 * n + 1) = -1 := by
      rw [pow_succ, pow_mul]; norm_num
    have hq : (-1 : ℚ) ^ (2 * n + 1 + 1) = 1 := by
      rw [pow_succ, hp]; norm_num
    rw [hp, hq]
    have h2r : (2 : ℚ) ^ r = 2 * 2 ^ (r - 1) := by
      rw [← pow_succ']; congr 1; omega
    have hn1 : ((2 * n + 1 + 1 : ℕ) : ℚ) ^ r = 2 ^ r * ((n : ℚ) + 1) ^ r := by
      rw [← mul_pow]; push_cast; ring
    rw [hn1, h2r]
    have : ((n : ℚ) + 1) ^ r ≠ 0 := by positivity
    have : (2 : ℚ) ^ (r - 1) ≠ 0 := by positivity
    field_simp
    ring

theorem proof3 : stmt_3 := by
  intro r n hr hn1
  refine ⟨?_, ?_⟩
  · rw [stmt3_aux r hr n]
    obtain ⟨m, rfl⟩ : ∃ m, n = m + 1 := ⟨n - 1, by omega⟩
    rw [hn_succ]; simp
  · rw [Z_single]
    have h1 : Int.sign (r : ℤ) = 1 := Int.sign_eq_one_of_pos (by omega)
    rw [h1]
    simp only [Int.natAbs_natCast, Int.cast_one, one_pow]
    obtain ⟨m, hm⟩ : ∃ m, 2 * n = m + 1 := ⟨2 * n - 1, by omega⟩
    rw [hm, Finset.sum_Icc_succ_top (by omega), show m + 1 - 1 = m by omega]
    have : (2 * (n : ℚ)) = ((m + 1 : ℕ) : ℚ) := by exact_mod_cast hm
    rw [this]; unfold hn; ring

end S03

theorem S03.proof_3 : S03.stmt_3 := S03.proof3

namespace S03
open Real

lemma H_m1 (x : ℝ) (hx : -1 < x) : H [-1] x = Real.log (1 + x) := by
  have : H [-1] x = ∫ t in (0 : ℝ)..x, (1 + t)⁻¹ := by
    simp [H, f]
  rw [this, intervalIntegral.integral_comp_add_left (fun t => t⁻¹) (1 : ℝ),
    integral_inv_of_pos (by norm_num) (by linarith)]
  simp

lemma H_p1 (x : ℝ) (hx : x < 1) : H [1] x = -Real.log (1 - x) := by
  have : H [1] x = ∫ t in (0 : ℝ)..x, (1 - t)⁻¹ := by
    simp [H, f]
  rw [this, intervalIntegral.integral_comp_sub_left (fun t => t⁻¹) (1 : ℝ),
    integral_inv_of_pos (by linarith) (by norm_num)]
  simp [Real.log_inv]

theorem proof10 : stmt_10 := by
  intro α h0 hπ
  set c := D04.c α with hc
  have hc0 : 0 < c := by
    rw [hc, D04.c]; apply Real.cos_pos_of_mem_Ioo; constructor <;> linarith
  have hc1 : c < 1 := by
    rw [hc, D04.c]
    have := Real.cos_lt_cos_of_nonneg_of_le_pi (le_refl 0) (by linarith : α / 2 ≤ π)
      (by linarith : 0 < α / 2)
    simpa using this
  have habs : |c| < 1 := by rw [abs_lt]; constructor <;> linarith
  have habs2 : |c ^ 2| < 1 := by
    rw [abs_pow]; exact pow_lt_one₀ (abs_nonneg _) habs (by norm_num)
  have hS := V05.proof_8 α h0 hπ 0 (by norm_num)
  have hP : ∀ n : ℕ, V05.P 0 (n + 1) V05.z2 V05.z3 * c ^ (2 * (n + 1)) =
      -(c * ((2 : ℝ) * (1 / (2 * n + 1)) * c ^ (2 * n + 1))) + (c ^ 2) ^ (n + 1) / (n + 1) := by
    intro n
    rw [(V05.proof_9 (n + 1) (by omega) V05.z2 V05.z3).1]
    rw [← pow_mul, show 2 * (n + 1) = 2 * n + 1 + 1 by ring, pow_succ]
    push_cast
    have h1 : (2 * ((n : ℝ) + 1) - 1) = 2 * n + 1 := by ring
    rw [h1]
    field_simp
    ring
  rw [← hc] at hS
  simp_rw [hP] at hS
  have hA := ((Real.hasSum_log_sub_log_of_abs_lt_one habs).mul_left c).neg
  have hB := Real.hasSum_pow_div_log_of_abs_lt_one habs2
  have hT := hA.add hB
  have hval := hS.unique hT
  rw [H_m1 c (by linarith), H_p1 c hc1]
  rw [hval]
  have : 1 - c ^ 2 = (1 + c) * (1 - c) := by ring
  rw [this, Real.log_mul (by linarith) (by linarith)]
  ring

end S03

theorem S03.proof_10 : S03.stmt_10 := S03.proof10

namespace S03

lemma Z_cons (a : ℤ) (v : List ℤ) (M : ℕ) :
    Z (a :: v) M = ∑ i ∈ Finset.Icc 1 M, (Int.sign a : ℚ) ^ i / (i : ℚ) ^ a.natAbs * Z v (i - 1) :=
  rfl

lemma Zw_cons (a : ℤ) (v : List ℤ) (M₀ M : ℕ) :
    Zw (a :: v) M₀ M = ∑ i ∈ Finset.Ioc M₀ M,
      (Int.sign a : ℚ) ^ i / (i : ℚ) ^ a.natAbs * Zw v M₀ (i - 1) := rfl

lemma window (v : List ℤ) : ∀ M₀ M : ℕ, M₀ ≤ M →
    Z v M = ∑ t ∈ Finset.range (v.length + 1), Zw (v.take t) M₀ M * Z (v.drop t) M₀ := by
  induction v with
  | nil => intro M₀ M _; simp [Z, Zw]
  | cons a v ih =>
    intro M₀ M h
    rw [List.length_cons, Finset.sum_range_succ']
    simp only [List.take_succ_cons, List.drop_succ_cons, List.take_zero, List.drop_zero]
    rw [show Zw [] M₀ M = 1 from rfl, one_mul, Z_cons, Z_cons, show Finset.Icc 1 M = Finset.Ioc 0 M from Finset.Icc_add_one_left_eq_Ioc 0 M,
      show Finset.Icc 1 M₀ = Finset.Ioc 0 M₀ from Finset.Icc_add_one_left_eq_Ioc 0 M₀, ← Finset.sum_Ioc_consecutive _ (Nat.zero_le M₀) h]
    rw [add_comm]
    congr 1
    rw [Finset.sum_congr rfl (fun i hi => by
      rw [ih M₀ (i - 1) (by simp at hi; omega)])]
    simp_rw [Finset.mul_sum, Zw_cons, Finset.sum_mul]
    rw [Finset.sum_comm]
    refine Finset.sum_congr rfl (fun t _ => Finset.sum_congr rfl (fun i _ => by ring))

end S03

theorem S03.proof_2 : S03.stmt_2 := fun v _ M₀ M h => S03.window v M₀ M h

namespace S03

lemma Wc_cons (a : ℤ) (v : List ℤ) (U : ℤ) (N : ℚ) :
    Wc (a :: v) U N = ∑ δ ∈ Finset.Ico (-1 : ℤ) U,
      (Int.sign a : ℚ) ^ δ / (N + (δ : ℚ)) ^ a.natAbs * Wc v δ N := rfl

lemma sign_pow_shift (a : ℤ) (ha : a ≠ 0) (N i : ℕ) (hN : Even N) :
    (Int.sign a : ℚ) ^ i = (Int.sign a : ℚ) ^ ((i : ℤ) - N) := by
  have h : Int.sign a = 1 ∨ Int.sign a = -1 := by
    rcases lt_or_gt_of_ne ha with h | h
    · exact Or.inr (Int.sign_eq_neg_one_of_neg h)
    · exact Or.inl (Int.sign_eq_one_of_pos h)
  rcases h with h | h <;> rw [h]
  · simp
  · rw [zpow_sub₀ (by norm_num), zpow_natCast, zpow_natCast]
    push_cast
    rw [hN.neg_one_pow]; simp

lemma Zw_eq_Wc (N : ℕ) (hN2 : 2 ≤ N) (hN : Even N) (v : List ℤ) (hv : IsIndex v) :
    ∀ M : ℕ, N - 2 ≤ M → Zw v (N - 2) M = Wc v ((M : ℤ) - N + 1) (N : ℚ) := by
  induction v with
  | nil => intro M _; rfl
  | cons a v ih =>
    intro M hM
    have ha : a ≠ 0 := hv a (by simp)
    have hv' : IsIndex v := fun b hb => hv b (by simp [hb])
    rw [Zw_cons, Wc_cons]
    refine Finset.sum_nbij' (fun i => (i : ℤ) - N) (fun δ => (δ + N).toNat) ?_ ?_ ?_ ?_ ?_
    · intro i hi; simp only [Finset.mem_Ioc] at hi
      simp only [Finset.mem_Ico]; omega
    · intro δ hδ; simp only [Finset.mem_Ico] at hδ
      simp only [Finset.mem_Ioc]; omega
    · intro i hi; simp only [Finset.mem_Ioc] at hi; omega
    · intro δ hδ; simp only [Finset.mem_Ico] at hδ; omega
    · intro i hi
      simp only [Finset.mem_Ioc] at hi
      rw [ih hv' (i - 1) (by omega), sign_pow_shift a ha N i hN]
      have h1 : ((i : ℚ)) = (N : ℚ) + (((i : ℤ) - N : ℤ) : ℚ) := by push_cast; ring
      have h2 : ((i - 1 : ℕ) : ℤ) - N + 1 = (i : ℤ) - N := by omega
      rw [← h1, h2]

theorem proof7a : (∀ v : List ℤ, IsIndex v → ∀ τ N : ℕ, 2 ≤ N → Even N →
    Z v (N - 2 + τ) = ∑ t ∈ Finset.range (v.length + 1),
      Wc (v.take t) ((τ : ℤ) - 1) (N : ℚ) * Z (v.drop t) (N - 2)) := by
  intro v hv τ N hN2 hN
  rw [window v (N - 2) (N - 2 + τ) (by omega)]
  refine Finset.sum_congr rfl (fun t _ => ?_)
  rw [Zw_eq_Wc N hN2 hN (v.take t) (fun b hb => hv b (List.mem_of_mem_take hb)) _ (by omega)]
  congr 2
  omega

end S03

namespace S03

lemma wt_cons (a : ℤ) (v : List ℤ) : wt (a :: v) = a.natAbs + wt v := by simp [wt]

/-- `D_U^m = Π_{δ ∈ [-1, U)} (N+δ)^m`. -/
lemma Wc_poly (u : List ℤ) : ∀ U : ℤ, ∃ Q : Polynomial ℚ, ∀ N : ℚ, 2 ≤ N →
    Wc u U N * ∏ δ ∈ Finset.Ico (-1 : ℤ) U, (N + (δ : ℚ)) ^ wt u = Q.eval N := by
  induction u with
  | nil => intro U; exact ⟨1, fun N _ => by simp [Wc, wt]⟩
  | cons a v ih =>
    intro U
    choose Qf hQf using ih
    refine ⟨∑ δ ∈ Finset.Ico (-1 : ℤ) U, Polynomial.C ((Int.sign a : ℚ) ^ δ) *
      (∏ δ' ∈ (Finset.Ico (-1 : ℤ) U).erase δ, (Polynomial.X + Polynomial.C (δ' : ℚ)) ^ a.natAbs) *
      (∏ δ' ∈ Finset.Ico δ U, (Polynomial.X + Polynomial.C (δ' : ℚ)) ^ wt v) * Qf δ, ?_⟩
    intro N hN
    rw [Wc_cons, Finset.sum_mul, Polynomial.eval_finsetSum]
    refine Finset.sum_congr rfl (fun δ hδ => ?_)
    have hδ' := Finset.mem_Ico.mp hδ
    simp only [Polynomial.eval_mul, Polynomial.eval_C, Polynomial.eval_prod, Polynomial.eval_pow,
      Polynomial.eval_add, Polynomial.eval_X]
    rw [← hQf δ N hN, wt_cons]
    simp_rw [pow_add]
    rw [Finset.prod_mul_distrib, ← Finset.mul_prod_erase _ _ hδ,
      ← Finset.Ico_union_Ico_eq_Ico (show (-1 : ℤ) ≤ δ by omega) (show δ ≤ U by omega),
      Finset.prod_union (Finset.Ico_disjoint_Ico_consecutive _ _ _)]
    have hne : N + (δ : ℚ) ≠ 0 := by
      have : (-1 : ℚ) ≤ δ := by exact_mod_cast hδ'.1
      linarith
    field_simp

theorem proof7b : (∀ v : List ℤ, IsIndex v → ∀ τ t : ℕ, ∃ Q : Polynomial ℚ, ∃ e : ℤ →₀ ℕ,
    (∀ p ∈ e.support, p ≤ 1 ∧ e p ≤ wt (v.take t)) ∧
    ∀ N : ℚ, 2 ≤ N → Wc (v.take t) ((τ : ℤ) - 1) N = ratEval Q e N) := by
  intro v _ τ t
  set u := v.take t
  set U : ℤ := (τ : ℤ) - 1
  obtain ⟨Q, hQ⟩ := Wc_poly u U
  set s := Finset.Ico (-1 : ℤ) U
  refine ⟨Q, ∑ δ ∈ s, Finsupp.single (-δ) (wt u), ?_, ?_⟩
  · intro p hp
    have hval : (∑ δ ∈ s, Finsupp.single (-δ) (wt u)) p = if -p ∈ s then wt u else 0 := by
      rw [Finsupp.finsetSum_apply]
      rw [Finset.sum_congr rfl (fun δ _ => show Finsupp.single (-δ) (wt u) p =
        if δ = -p then wt u else 0 by
          rw [Finsupp.single_apply]; congr 1; apply propext; omega)]
      rw [Finset.sum_ite_eq']
    rw [Finsupp.mem_support_iff, hval] at hp
    rw [hval]
    split_ifs at hp ⊢ with h
    · refine ⟨?_, le_refl _⟩
      have := (Finset.mem_Ico.mp h).1; omega
    · exact absurd rfl hp
  · intro N hN
    have hprod : (∑ δ ∈ s, Finsupp.single (-δ) (wt u)).prod (fun p k => (N - (p : ℚ)) ^ k) =
        ∏ δ ∈ s, (N + (δ : ℚ)) ^ wt u := by
      rw [← Finsupp.prod_finsetSum_index (h := fun (p : ℤ) (k : ℕ) => (N - (p : ℚ)) ^ k) (fun _ => by simp)
        (fun _ _ _ => by rw [pow_add])]
      refine Finset.prod_congr rfl (fun δ _ => ?_)
      rw [Finsupp.prod_single_index (by simp)]
      push_cast; ring_nf
    unfold ratEval
    rw [hprod, ← hQ N hN]
    have : ∏ δ ∈ s, (N + (δ : ℚ)) ^ wt u ≠ 0 := by
      rw [Finset.prod_ne_zero_iff]
      intro δ hδ
      have : (-1 : ℚ) ≤ δ := by exact_mod_cast (Finset.mem_Ico.mp hδ).1
      apply pow_ne_zero; linarith
    field_simp

end S03

theorem S03.proof_7 : S03.stmt_7 := ⟨S03.proof7a, S03.proof7b⟩

namespace S03

/-- `g_a(i) = σ_a^i / i^{|a|}`. -/
def gq (a : ℤ) (i : ℕ) : ℚ := (Int.sign a : ℚ) ^ i / (i : ℚ) ^ a.natAbs

lemma Z_cons' (a : ℤ) (v : List ℤ) (M : ℕ) :
    Z (a :: v) M = ∑ i ∈ Finset.Icc 1 M, gq a i * Z v (i - 1) := rfl

lemma Z_cons_succ (a : ℤ) (v : List ℤ) (M : ℕ) :
    Z (a :: v) (M + 1) = Z (a :: v) M + gq a (M + 1) * Z v M := by
  rw [Z_cons', Z_cons', Finset.sum_Icc_succ_top (by omega)]; simp

lemma stuffle_core (a b : ℤ) (u v : List ℤ) (M : ℕ) :
    Z (a :: u) M * Z (b :: v) M =
      ∑ i ∈ Finset.Icc 1 M, gq a i * (Z u (i - 1) * Z (b :: v) (i - 1)) +
      ∑ i ∈ Finset.Icc 1 M, gq b i * (Z (a :: u) (i - 1) * Z v (i - 1)) +
      ∑ i ∈ Finset.Icc 1 M, gq a i * gq b i * (Z u (i - 1) * Z v (i - 1)) := by
  induction M with
  | zero => simp [Z]
  | succ M ih =>
    rw [Z_cons_succ, Z_cons_succ, Finset.sum_Icc_succ_top (by omega),
      Finset.sum_Icc_succ_top (by omega), Finset.sum_Icc_succ_top (by omega)]
    simp only [Nat.add_sub_cancel]
    have := ih
    linear_combination this

lemma sum_map_cons (a : ℤ) (L : List (List ℤ)) (M : ℕ) :
    (L.map fun w => Z (a :: w) M).sum =
      ∑ i ∈ Finset.Icc 1 M, gq a i * (L.map fun w => Z w (i - 1)).sum := by
  induction L with
  | nil => simp
  | cons w L ih =>
    simp only [List.map_cons, List.sum_cons]
    rw [ih, Z_cons', ← Finset.sum_add_distrib]
    refine Finset.sum_congr rfl (fun i _ => by ring)

lemma sign_ne (a : ℤ) (ha : a ≠ 0) : Int.sign a = 1 ∨ Int.sign a = -1 := by
  rcases lt_or_gt_of_ne ha with h | h
  · exact Or.inr (Int.sign_eq_neg_one_of_neg h)
  · exact Or.inl (Int.sign_eq_one_of_pos h)

lemma gq_mul (a b : ℤ) (ha : a ≠ 0) (hb : b ≠ 0) (i : ℕ) :
    gq a i * gq b i = gq (Int.sign a * Int.sign b * ((a.natAbs + b.natAbs : ℕ) : ℤ)) i ∧
    (Int.sign a * Int.sign b * ((a.natAbs + b.natAbs : ℕ) : ℤ)).natAbs = a.natAbs + b.natAbs ∧
    Int.sign a * Int.sign b * ((a.natAbs + b.natAbs : ℕ) : ℤ) ≠ 0 := by
  have hpos : 0 < a.natAbs + b.natAbs := by
    have := Int.natAbs_pos.mpr ha; omega
  have hs : Int.sign (((a.natAbs + b.natAbs : ℕ) : ℤ)) = 1 := Int.sign_eq_one_of_pos (by omega)
  have hna : (Int.sign a * Int.sign b * ((a.natAbs + b.natAbs : ℕ) : ℤ)).natAbs =
      a.natAbs + b.natAbs := by
    rw [Int.natAbs_mul, Int.natAbs_mul, Int.natAbs_natCast]
    rcases sign_ne a ha with h1 | h1 <;> rcases sign_ne b hb with h2 | h2 <;>
      rw [h1, h2] <;> simp
  refine ⟨?_, hna, ?_⟩
  · unfold gq
    rw [hna, Int.sign_mul, Int.sign_mul, hs, Int.sign_sign, Int.sign_sign]
    push_cast
    rw [mul_one, mul_pow, pow_add]
    ring
  · intro h
    have := congrArg Int.natAbs h
    rw [hna] at this; simp at this; omega

lemma stuffle_list : ∀ u v : List ℤ, IsIndex u → IsIndex v →
    ∃ L : List (List ℤ), (∀ w ∈ L, IsIndex w ∧ wt w = wt u + wt v) ∧
      ∀ M : ℕ, Z u M * Z v M = (L.map fun w => Z w M).sum := by
  intro u
  induction u with
  | nil =>
    intro v _ hv
    exact ⟨[v], by simp [hv, wt], fun M => by simp [Z]⟩
  | cons a u ihu =>
    intro v hu hv
    have ha : a ≠ 0 := hu a (by simp)
    have hu' : IsIndex u := fun x hx => hu x (by simp [hx])
    induction v with
    | nil => exact ⟨[a :: u], by simp [hu, wt], fun M => by simp [Z]⟩
    | cons b v ihv =>
      have hb : b ≠ 0 := hv b (by simp)
      have hv' : IsIndex v := fun x hx => hv x (by simp [hx])
      obtain ⟨L1, hL1, hZ1⟩ := ihu (b :: v) hu' hv
      obtain ⟨L2, hL2, hZ2⟩ := ihv hv'
      obtain ⟨L3, hL3, hZ3⟩ := ihu v hu' hv'
      obtain ⟨_, hc, hc0⟩ := gq_mul a b ha hb 0
      have hg := fun i => (gq_mul a b ha hb i).1
      set c := Int.sign a * Int.sign b * ((a.natAbs + b.natAbs : ℕ) : ℤ)
      refine ⟨L1.map (a :: ·) ++ L2.map (b :: ·) ++ L3.map (c :: ·), ?_, ?_⟩
      · intro w hw
        simp only [List.mem_append, List.mem_map] at hw
        rcases hw with (⟨w', hw', rfl⟩ | ⟨w', hw', rfl⟩) | ⟨w', hw', rfl⟩
        · obtain ⟨h1, h2⟩ := hL1 w' hw'
          refine ⟨fun x hx => ?_, ?_⟩
          · simp at hx; rcases hx with rfl | hx; exact ha; exact h1 x hx
          · simp only [wt_cons] at h2 ⊢; omega
        · obtain ⟨h1, h2⟩ := hL2 w' hw'
          refine ⟨fun x hx => ?_, ?_⟩
          · simp at hx; rcases hx with rfl | hx; exact hb; exact h1 x hx
          · simp only [wt_cons] at h2 ⊢; omega
        · obtain ⟨h1, h2⟩ := hL3 w' hw'
          refine ⟨fun x hx => ?_, ?_⟩
          · simp at hx; rcases hx with rfl | hx; exact hc0; exact h1 x hx
          · simp only [wt_cons, hc] at h2 ⊢; omega
      · intro M
        rw [List.map_append, List.map_append, List.sum_append, List.sum_append,
          List.map_map, List.map_map, List.map_map]
        rw [show ((fun w => Z w M) ∘ (a :: ·)) = fun w => Z (a :: w) M from rfl,
          show ((fun w => Z w M) ∘ (b :: ·)) = fun w => Z (b :: w) M from rfl,
          show ((fun w => Z w M) ∘ (c :: ·)) = fun w => Z (c :: w) M from rfl,
          sum_map_cons, sum_map_cons, sum_map_cons, stuffle_core]
        simp_rw [← hZ1, ← hZ2, ← hZ3, hg]

theorem proof1 : stmt_1 := by
  intro u v hu hv
  obtain ⟨L, hL, hZ⟩ := stuffle_list u v hu hv
  refine ⟨(L : Multiset (List ℤ)), fun w hw => hL w (Multiset.mem_coe.mp hw), fun M => ?_⟩
  rw [hZ M, Multiset.map_coe, Multiset.sum_coe]

end S03

theorem S03.proof_1 : S03.stmt_1 := S03.proof1

namespace S03
open Real MeasureTheory Set

lemma H_cons (a : ℤ) (w : List ℤ) (x : ℝ) :
    H (a :: w) x = ∫ t in (0 : ℝ)..x, f a t * H w t := rfl

lemma H_neg (w : List ℤ) (hw : IsLetters w) : ∀ x : ℝ,
    H w (-x) = (-1 : ℝ) ^ (w.filter fun a => a ≠ 0).length * H (w.map fun a => -a) x := by
  induction w with
  | nil => intro x; simp [H]
  | cons a w ih =>
    intro x
    have hw' : IsLetters w := fun b hb => hw b (by simp [hb])
    have ha := hw a (by simp)
    have key : H (a :: w) (-x) = -∫ s in (0 : ℝ)..x, f a (-s) * H w (-s) := by
      rw [intervalIntegral.integral_comp_neg (fun t => f a t * H w t), H_cons]
      rw [neg_zero, intervalIntegral.integral_symm 0 (-x), neg_neg]
    rw [key]
    simp_rw [ih hw']
    simp only [List.map_cons, H_cons, List.filter_cons]
    rcases ha with rfl | rfl | rfl
    · simp only [ne_eq, not_true_eq_false, decide_false, Bool.false_eq_true, ite_false, neg_zero]
      rw [← intervalIntegral.integral_neg, ← intervalIntegral.integral_const_mul]
      congr 1; ext s; simp [f]; ring
    · simp only [ne_eq, one_ne_zero, not_false_eq_true, decide_true, ite_true, List.length_cons]
      rw [← intervalIntegral.integral_const_mul, ← intervalIntegral.integral_neg]
      congr 1; ext s; simp [f]; ring
    · simp only [ne_eq, neg_eq_zero, one_ne_zero, not_false_eq_true, decide_true, ite_true,
        List.length_cons, neg_neg]
      rw [← intervalIntegral.integral_const_mul, ← intervalIntegral.integral_neg]
      congr 1; ext s; simp [f]; ring

end S03

namespace S03
open Real MeasureTheory Set
open scoped Interval

lemma reg_tail (a : ℤ) (w : List ℤ) (h : Regular (a :: w)) : Regular w := by
  cases w with
  | nil => intro b hb; simp at hb
  | cons b w => intro c hc; exact h c (by simpa [List.getLast?_cons_cons] using hc)

lemma reg_zero_ne (w : List ℤ) (h : Regular (0 :: w)) : w ≠ [] := by
  rintro rfl; exact h 0 (by simp) rfl

lemma intInt_of_bdd (g : ℝ → ℝ) (r : ℝ) (hc : ContinuousOn g (Icc (-r) r \ {0})) (K : ℝ)
    (hK : ∀ t ∈ Icc (-r) r, |g t| ≤ K) (x : ℝ) (hx : |x| ≤ r) :
    IntervalIntegrable g volume 0 x := by
  have hsub : Ι (0 : ℝ) x ⊆ Icc (-r) r := by
    intro t ht
    have := uIoc_subset_uIcc ht
    rw [mem_uIcc] at this
    rw [abs_le] at hx
    constructor <;> rcases this with h | h <;> linarith [h.1, h.2]
  rw [intervalIntegrable_iff]
  refine IntegrableOn.of_bound measure_Ioc_lt_top ?_ K ?_
  · have := (hc.mono (sdiff_subset_sdiff_left hsub)).aestronglyMeasurable
      (measurableSet_uIoc.diff (measurableSet_singleton 0)) (μ := volume)
    rwa [Measure.restrict_congr_set (sdiff_null_ae_eq_self (measure_singleton 0))] at this
  · exact ae_restrict_of_forall_mem measurableSet_uIoc (fun t ht => by
      rw [Real.norm_eq_abs]; exact hK t (hsub ht))

lemma f_contOn (a : ℤ) : ContinuousOn (f a) (Ioo (-1) 1 \ {0}) := by
  intro t ht
  obtain ⟨⟨h1, h2⟩, h0⟩ := ht
  have h0 : t ≠ 0 := h0
  apply ContinuousAt.continuousWithinAt
  unfold f
  split_ifs
  · exact continuousAt_const.div continuousAt_id h0
  · exact continuousAt_const.div (continuousAt_const.sub continuousAt_id)
      (ne_of_gt (show (1:ℝ) - t > 0 by linarith))
  · exact continuousAt_const.div (continuousAt_const.add continuousAt_id)
      (ne_of_gt (show (1:ℝ) + t > 0 by linarith))

lemma f_bound (a : ℤ) (ha : a ≠ 0) (r : ℝ) (hr : r < 1) (t : ℝ) (ht : |t| ≤ r) :
    |f a t| ≤ 1 / (1 - r) := by
  rw [abs_le] at ht
  have hr1 : 0 < 1 - r := by linarith
  unfold f
  rw [if_neg ha]
  split_ifs
  · rw [abs_of_pos (by apply div_pos one_pos; linarith)]
    apply one_div_le_one_div_of_le hr1; linarith
  · rw [abs_of_pos (by apply div_pos one_pos; linarith)]
    apply one_div_le_one_div_of_le hr1; linarith

/-- Continuity of `H_w` on `(-1,1)` and `|H_w(t)| ≤ C|t|` on `[-r,r]` for nonempty `w`. -/
lemma H_good (w : List ℤ) (hw : IsLetters w) (hr : Regular w) :
    ContinuousOn (H w) (Ioo (-1) 1) ∧
    (w ≠ [] → ∀ r : ℝ, 0 ≤ r → r < 1 → ∃ C : ℝ, 0 ≤ C ∧ ∀ t : ℝ, |t| ≤ r → |H w t| ≤ C * |t|) ∧
    (∀ r : ℝ, 0 ≤ r → r < 1 → ∃ B : ℝ, 0 ≤ B ∧ ∀ t : ℝ, |t| ≤ r → |H w t| ≤ B) := by
  induction w with
  | nil =>
    refine ⟨continuousOn_const, fun h => absurd rfl h, fun r _ _ => ⟨1, by norm_num, fun t _ => ?_⟩⟩
    simp [H]
  | cons a w ih =>
    have hw' : IsLetters w := fun b hb => hw b (by simp [hb])
    have hr' := reg_tail a w hr
    obtain ⟨ihc, ihC, ihB⟩ := ih hw' hr'
    -- the integrand
    set g := fun t => f a t * H w t with hg
    have gc : ContinuousOn g (Ioo (-1) 1 \ {0}) :=
      (f_contOn a).mul (ihc.mono sdiff_subset)
    have gK : ∀ r : ℝ, 0 ≤ r → r < 1 → ∃ K : ℝ, 0 ≤ K ∧ ∀ t ∈ Icc (-r) r, |g t| ≤ K := by
      intro r hr0 hr1
      by_cases ha : a = 0
      · subst ha
        obtain ⟨C, hC0, hC⟩ := ihC (reg_zero_ne w hr) r hr0 hr1
        refine ⟨C, hC0, fun t ht => ?_⟩
        have ht' : |t| ≤ r := abs_le.mpr ht
        simp only [hg, f, if_pos]
        by_cases ht0 : t = 0
        · subst ht0; simp [hC0]
        · rw [abs_mul, abs_div, abs_one]
          have := hC t ht'
          have hpos : 0 < |t| := abs_pos.mpr ht0
          rw [div_mul_eq_mul_div, one_mul, div_le_iff₀ hpos]; exact this
      · obtain ⟨B, hB0, hB⟩ := ihB r hr0 hr1
        refine ⟨1 / (1 - r) * B, by apply mul_nonneg; apply div_nonneg; norm_num; linarith; exact hB0,
          fun t ht => ?_⟩
        have ht' : |t| ≤ r := abs_le.mpr ht
        rw [abs_mul]
        exact mul_le_mul (f_bound a ha r hr1 t ht') (hB t ht') (abs_nonneg _)
          (by apply div_nonneg; norm_num; linarith)
    have gI : ∀ r : ℝ, 0 ≤ r → r < 1 → ∀ x : ℝ, |x| ≤ r → IntervalIntegrable g volume 0 x := by
      intro r hr0 hr1 x hx
      obtain ⟨K, _, hK⟩ := gK r hr0 hr1
      refine intInt_of_bdd g r (gc.mono ?_) K hK x hx
      intro t ⟨⟨h1, h2⟩, h3⟩
      exact ⟨⟨by linarith, by linarith⟩, h3⟩
    have Hlin : ∀ r : ℝ, 0 ≤ r → r < 1 → ∃ C : ℝ, 0 ≤ C ∧ ∀ t : ℝ, |t| ≤ r → |H (a :: w) t| ≤ C * |t| := by
      intro r hr0 hr1
      obtain ⟨K, hK0, hK⟩ := gK r hr0 hr1
      refine ⟨K, hK0, fun t ht => ?_⟩
      rw [H_cons]
      have := intervalIntegral.norm_integral_le_of_norm_le_const (a := 0) (b := t) (C := K)
        (f := g) (fun s hs => by
          rw [Real.norm_eq_abs]
          apply hK
          have := uIoc_subset_uIcc hs
          rw [mem_uIcc] at this
          rw [abs_le] at ht
          constructor <;> rcases this with h | h <;> linarith [h.1, h.2])
      rw [Real.norm_eq_abs, sub_zero] at this
      exact this
    refine ⟨?_, fun _ => Hlin, ?_⟩
    · intro x0 hx0
      obtain ⟨h1, h2⟩ := hx0
      set r := (|x0| + 1) / 2 with hrdef
      have hax : |x0| < 1 := abs_lt.mpr ⟨h1, h2⟩
      have hr0 : 0 ≤ r := by positivity
      have hr1 : r < 1 := by rw [hrdef]; linarith
      have hxr : |x0| < r := by rw [hrdef]; linarith
      have hint : IntervalIntegrable g volume (-r) r := by
        apply IntervalIntegrable.trans (b := 0)
        · exact (gI r hr0 hr1 (-r) (by rw [abs_neg, abs_of_nonneg hr0])).symm
        · exact gI r hr0 hr1 r (by rw [abs_of_nonneg hr0])
      have hcont := intervalIntegral.continuousOn_primitive_interval' hint
        (a := 0) (by rw [mem_uIcc]; left; constructor <;> linarith)
      have hmem : uIcc (-r) r ∈ nhds x0 := by
        rw [uIcc_of_le (by linarith)]
        apply Icc_mem_nhds <;> [skip; skip] <;> (rw [abs_lt] at hxr; linarith [hxr.1, hxr.2])
      exact (hcont.continuousAt hmem).continuousWithinAt
    · intro r hr0 hr1
      obtain ⟨C, hC0, hC⟩ := Hlin r hr0 hr1
      refine ⟨C * r, mul_nonneg hC0 hr0, fun t ht => ?_⟩
      exact (hC t ht).trans (mul_le_mul_of_nonneg_left ht hC0)

lemma integrand_ok (a : ℤ) (w : List ℤ) (hw : IsLetters (a :: w)) (hr : Regular (a :: w)) :
    (∀ x : ℝ, |x| < 1 → IntervalIntegrable (fun t => f a t * H w t) volume 0 x) ∧
    ContinuousOn (fun t => f a t * H w t) (Ioo (-1) 1 \ {0}) := by
  have hw' : IsLetters w := fun b hb => hw b (by simp [hb])
  have hr' := reg_tail a w hr
  obtain ⟨ihc, ihC, ihB⟩ := H_good w hw' hr'
  have gc : ContinuousOn (fun t => f a t * H w t) (Ioo (-1) 1 \ {0}) :=
    (f_contOn a).mul (ihc.mono sdiff_subset)
  refine ⟨fun x hx => ?_, gc⟩
  set r := |x|
  have hr0 : 0 ≤ r := abs_nonneg x
  by_cases ha : a = 0
  · subst ha
    obtain ⟨C, hC0, hC⟩ := ihC (reg_zero_ne w hr) r hr0 hx
    refine intInt_of_bdd _ r (gc.mono ?_) C (fun t ht => ?_) x le_rfl
    · intro t ⟨⟨h1, h2⟩, h3⟩
      exact ⟨⟨by linarith, by linarith⟩, h3⟩
    · have ht' : |t| ≤ r := abs_le.mpr ht
      simp only [f, if_pos]
      by_cases ht0 : t = 0
      · subst ht0; simp [hC0]
      · rw [abs_mul, abs_div, abs_one]
        have := hC t ht'
        have hpos : 0 < |t| := abs_pos.mpr ht0
        rw [div_mul_eq_mul_div, one_mul, div_le_iff₀ hpos]; exact this
  · obtain ⟨B, hB0, hB⟩ := ihB r hr0 hx
    refine intInt_of_bdd _ r (gc.mono ?_) (1 / (1 - r) * B) (fun t ht => ?_) x le_rfl
    · intro t ⟨⟨h1, h2⟩, h3⟩
      exact ⟨⟨by linarith, by linarith⟩, h3⟩
    · have ht' : |t| ≤ r := abs_le.mpr ht
      rw [abs_mul]
      exact mul_le_mul (f_bound a ha r hx t ht') (hB t ht') (abs_nonneg _)
        (by apply div_nonneg; norm_num; linarith)

lemma H_deriv (a : ℤ) (w : List ℤ) (hw : IsLetters (a :: w)) (hr : Regular (a :: w)) (x : ℝ)
    (hx0 : 0 < |x|) (hx1 : |x| < 1) : HasDerivAt (H (a :: w)) (f a x * H w x) x := by
  obtain ⟨hI, hc⟩ := integrand_ok a w hw hr
  have hxs : x ∈ Ioo (-1 : ℝ) 1 \ {0} :=
    ⟨abs_lt.mp hx1, by simpa [abs_pos] using hx0⟩
  have hopen : IsOpen (Ioo (-1 : ℝ) 1 \ {0}) := isOpen_Ioo.sdiff isClosed_singleton
  have := intervalIntegral.integral_hasDerivAt_right (hI x hx1)
    (hc.stronglyMeasurableAtFilter hopen x hxs) (hc.continuousAt (hopen.mem_nhds hxs))
  exact this

end S03

namespace S03
open Real MeasureTheory Set
open scoped Interval

lemma geo_poly_summable (D : ℕ) (r : ℝ) (hr0 : 0 < r) (hr : r < 1) :
    Summable (fun n : ℕ => ((n : ℝ) + 1) ^ D * r ^ n) := by
  have h := (summable_nat_add_iff 1).mpr (summable_pow_mul_geometric_of_norm_lt_one D
    (r := r) (by rw [Real.norm_eq_abs, abs_of_pos hr0]; exact hr))
  refine (h.mul_left (1 / r)).congr (fun n => ?_)
  push_cast; field_simp; ring

lemma psD (c : ℕ → ℝ) (D : ℕ) (hc : ∀ n, |c n| ≤ ((n : ℝ) + 1) ^ D) (x : ℝ) (hx : |x| < 1) :
    Summable (fun n => c n * x ^ (n + 1)) ∧ Summable (fun n => c n * ((n : ℝ) + 1) * x ^ n) ∧
    HasDerivAt (fun y => ∑' n, c n * y ^ (n + 1)) (∑' n, c n * ((n : ℝ) + 1) * x ^ n) x := by
  set r := (|x| + 1) / 2 with hrdef
  have hr0 : 0 < r := by positivity
  have hr1 : r < 1 := by rw [hrdef]; linarith
  have hxr : |x| < r := by rw [hrdef]; linarith
  have hu := geo_poly_summable (D + 1) r hr0 hr1
  have hu' := geo_poly_summable D r hr0 hr1
  have bnd : ∀ n (y : ℝ), |y| < r → ‖c n * ((n : ℝ) + 1) * y ^ n‖ ≤ ((n : ℝ) + 1) ^ (D + 1) * r ^ n := by
    intro n y hy
    rw [Real.norm_eq_abs, abs_mul, abs_mul, abs_pow, abs_of_pos (by positivity : (0:ℝ) < n + 1),
      pow_succ]
    have h1 : |y| ^ n ≤ r ^ n := pow_le_pow_left₀ (abs_nonneg _) hy.le n
    calc |c n| * ((n : ℝ) + 1) * |y| ^ n ≤ ((n : ℝ) + 1) ^ D * ((n : ℝ) + 1) * r ^ n := by
          gcongr; exact hc n
      _ = _ := by ring
  refine ⟨?_, ?_, ?_⟩
  · refine hu'.of_norm_bounded (fun n => ?_)
    rw [Real.norm_eq_abs, abs_mul, abs_pow]
    have h1 : |x| ^ (n + 1) ≤ r ^ n := by
      rw [pow_succ]
      calc |x| ^ n * |x| ≤ r ^ n * 1 := by
            gcongr
        _ = r ^ n := mul_one _
    exact mul_le_mul (hc n) h1 (by positivity) (by positivity)
  · exact hu.of_norm_bounded (fun n => bnd n x hxr)
  · refine hasDerivAt_tsum_of_isPreconnected (t := Ioo (-r) r) (y₀ := 0)
      (g := fun (n : ℕ) (y : ℝ) => c n * y ^ (n + 1))
      (g' := fun (n : ℕ) (y : ℝ) => c n * ((n : ℝ) + 1) * y ^ n)
      (u := fun (n : ℕ) => ((n : ℝ) + 1) ^ (D + 1) * r ^ n) hu isOpen_Ioo
      isPreconnected_Ioo (fun n y _ => ?_) (fun n y hy => bnd n y ?_) ?_ ?_ ?_
    · have := (hasDerivAt_pow (n + 1) y).const_mul (c n)
      convert this using 1; push_cast; ring
    · rw [abs_lt]; exact hy
    · constructor <;> linarith
    · simp
    · rw [mem_Ioo]; rw [abs_lt] at hxr; exact hxr

lemma ftc_series (c : ℕ → ℝ) (D : ℕ) (hc : ∀ n, |c n| ≤ ((n : ℝ) + 1) ^ D) (g : ℝ → ℝ) (x : ℝ)
    (hx : |x| < 1) (hgi : IntervalIntegrable g volume 0 x)
    (hg : ∀ t : ℝ, |t| < |x| → t ≠ 0 → ∑' n, c n * ((n : ℝ) + 1) * t ^ n = g t) :
    HasSum (fun n => c n * x ^ (n + 1)) (∫ t in (0 : ℝ)..x, g t) := by
  have hS : ∀ t ∈ uIcc 0 x, HasDerivAt (fun y => ∑' n, c n * y ^ (n + 1))
      (∑' n, c n * ((n : ℝ) + 1) * t ^ n) t := by
    intro t ht
    have : |t| ≤ |x| := by
      rw [mem_uIcc] at ht
      rcases ht with h | h <;> rw [abs_le] <;> constructor <;>
        cases abs_cases x <;> linarith [h.1, h.2]
    exact (psD c D hc t (by linarith)).2.2
  have hEq : EqOn (fun t => ∑' n, c n * ((n : ℝ) + 1) * t ^ n) g (uIoo 0 x) := by
    intro t ht
    rcases le_total 0 x with h | h
    · simp only [uIoo, mem_Ioo, inf_of_le_left h, sup_of_le_right h] at ht
      exact hg t (by rw [abs_of_pos ht.1, abs_of_nonneg h]; exact ht.2) (ne_of_gt ht.1)
    · simp only [uIoo, mem_Ioo, inf_of_le_right h, sup_of_le_left h] at ht
      exact hg t (by rw [abs_of_neg ht.2, abs_of_nonpos h]; linarith [ht.1]) (ne_of_lt ht.2)
  have hint := hgi.congr_uIoo hEq.symm
  rw [← intervalIntegral.integral_congr_uIoo hEq,
    intervalIntegral.integral_eq_sub_of_hasDerivAt hS hint]
  simp only [ne_eq, Nat.add_eq_zero_iff, one_ne_zero, and_false, not_false_eq_true, zero_pow,
    mul_zero, tsum_zero, sub_zero]
  exact (psD c D hc x hx).1.hasSum

end S03

namespace S03
open Real MeasureTheory Set
open scoped Interval

lemma Z_bound (v : List ℤ) : ∀ M : ℕ, |(Z v M : ℝ)| ≤ ((M : ℝ) + 1) ^ v.length := by
  induction v with
  | nil => intro M; simp [Z]
  | cons a v ih =>
    intro M
    rw [Z_cons]
    push_cast
    calc |∑ i ∈ Finset.Icc 1 M, ((Int.sign a : ℝ)) ^ i / (i : ℝ) ^ a.natAbs * (Z v (i - 1) : ℝ)|
        ≤ ∑ i ∈ Finset.Icc 1 M, |((Int.sign a : ℝ)) ^ i / (i : ℝ) ^ a.natAbs * (Z v (i - 1) : ℝ)| :=
          Finset.abs_sum_le_sum_abs _ _
      _ ≤ ∑ i ∈ Finset.Icc 1 M, ((M : ℝ) + 1) ^ v.length := by
          apply Finset.sum_le_sum
          intro i hi
          rw [Finset.mem_Icc] at hi
          rw [abs_mul]
          have h1 : |((Int.sign a : ℝ)) ^ i / (i : ℝ) ^ a.natAbs| ≤ 1 := by
            rw [abs_div, abs_pow, abs_pow]
            have hs : |((Int.sign a : ℤ) : ℝ)| ≤ 1 := by
              rcases Int.sign_trichotomy a with h | h | h <;> rw [h] <;> simp
            have hi1 : (1 : ℝ) ≤ |(i : ℝ)| := by
              rw [abs_of_nonneg (by positivity)]; exact_mod_cast hi.1
            rw [div_le_one (by positivity)]
            calc |((Int.sign a : ℤ) : ℝ)| ^ i ≤ 1 := pow_le_one₀ (abs_nonneg _) hs
              _ ≤ |(i : ℝ)| ^ a.natAbs := one_le_pow₀ hi1
          have h2 := ih (i - 1)
          have h3 : ((i - 1 : ℕ) : ℝ) + 1 ≤ (M : ℝ) + 1 := by
            have : i - 1 ≤ M := by omega
            exact_mod_cast (by omega : i - 1 + 1 ≤ M + 1)
          calc |((Int.sign a : ℝ)) ^ i / (i : ℝ) ^ a.natAbs| * |(Z v (i - 1) : ℝ)|
              ≤ 1 * ((M : ℝ) + 1) ^ v.length := by
                apply mul_le_mul h1 (h2.trans (pow_le_pow_left₀ (by positivity) h3 _))
                  (abs_nonneg _) (by norm_num)
            _ = _ := one_mul _
      _ ≤ ((M : ℝ) + 1) ^ (v.length + 1) := by
          rw [Finset.sum_const, Nat.card_Icc, nsmul_eq_mul, pow_succ]
          simp only [Nat.add_sub_cancel]
          rw [mul_comm]
          gcongr
          linarith

/-- The coefficients of `g6:lem:dict`. -/
noncomputable def cf (w : List ℤ) (N : ℕ) : ℝ :=
  ((dd w).2.2.1 : ℝ) * ((dd w).1 : ℝ) ^ (N + 1) / ((N : ℝ) + 1) ^ (dd w).2.1 *
    (Z (dd w).2.2.2 N : ℝ)

/-- The combinatorial part of `g6:lem:dict`. -/
def Facts (w : List ℤ) : Prop :=
  ((dd w).1 = 1 ∨ (dd w).1 = -1) ∧ 1 ≤ (dd w).2.1 ∧
    ((dd w).2.2.1 = 1 ∨ (dd w).2.2.1 = -1) ∧ IsIndex (dd w).2.2.2 ∧
    w.length = (dd w).2.1 + wt (dd w).2.2.2

lemma pm_abs (e : ℤ) (h : e = 1 ∨ e = -1) (n : ℕ) : |((e : ℝ)) ^ n| = 1 := by
  rcases h with rfl | rfl <;> simp

lemma cf_bound (w : List ℤ) (hF : Facts w) (N : ℕ) :
    |cf w N| ≤ ((N : ℝ) + 1) ^ (dd w).2.2.2.length := by
  obtain ⟨h1, h2, h3, _, _⟩ := hF
  unfold cf
  rw [abs_mul, abs_div, abs_mul, pm_abs _ h1, show |((dd w).2.2.1 : ℝ)| = 1 by
    have := pm_abs _ h3 1; simpa using this, one_mul]
  have hk : 1 ≤ |((N : ℝ) + 1) ^ (dd w).2.1| := by
    rw [abs_of_pos (by positivity)]; exact one_le_pow₀ (by linarith [(Nat.cast_nonneg N : (0:ℝ) ≤ N)])
  calc 1 / |((N : ℝ) + 1) ^ (dd w).2.1| * |(Z (dd w).2.2.2 N : ℝ)| ≤ 1 * |(Z (dd w).2.2.2 N : ℝ)| := by
        gcongr
        rw [div_le_one (by linarith)]; exact hk
    _ ≤ _ := by rw [one_mul]; exact Z_bound _ N

lemma dd_zero (b : ℤ) (w : List ℤ) : dd (0 :: b :: w) =
    ((dd (b :: w)).1, (dd (b :: w)).2.1 + 1, (dd (b :: w)).2.2.1, (dd (b :: w)).2.2.2) := by
  simp [dd]

lemma dd_nz (a b : ℤ) (w : List ℤ) (ha : a ≠ 0) : dd (a :: b :: w) =
    (a, 1, a * (dd (b :: w)).2.2.1,
      (a * (dd (b :: w)).1 * (((dd (b :: w)).2.1 : ℕ) : ℤ)) :: (dd (b :: w)).2.2.2) := by
  simp [dd, ha]

lemma log_series_m1 (x : ℝ) (hx : |x| < 1) :
    HasSum (fun N : ℕ => cf [-1] N * x ^ (N + 1)) (H [-1] x) := by
  have h := (Real.hasSum_pow_div_log_of_abs_lt_one (x := -x) (by rwa [abs_neg])).neg
  rw [sub_neg_eq_add, neg_neg] at h
  rw [H_m1 x (by rw [abs_lt] at hx; linarith)]
  refine h.congr_fun (fun N => ?_)
  simp only [cf, dd, Z]
  push_cast
  rw [neg_pow x]
  ring

lemma log_series_p1 (x : ℝ) (hx : |x| < 1) :
    HasSum (fun N : ℕ => cf [1] N * x ^ (N + 1)) (H [1] x) := by
  have h := Real.hasSum_pow_div_log_of_abs_lt_one hx
  rw [← H_p1 x (by rw [abs_lt] at hx; linarith)] at h
  refine h.congr_fun (fun N => ?_)
  simp only [cf, dd, Z]
  push_cast
  ring

end S03

namespace S03
open Real MeasureTheory Set
open scoped Interval

lemma sign_u (a e : ℤ) (k : ℕ) (ha : a = 1 ∨ a = -1) (he : e = 1 ∨ e = -1) (hk : 1 ≤ k) :
    Int.sign (a * e * (k : ℤ)) = a * e ∧ (a * e * (k : ℤ)).natAbs = k := by
  have hs : Int.sign (k : ℤ) = 1 := Int.sign_eq_one_of_pos (by omega)
  rw [Int.sign_mul, Int.sign_mul, hs, Int.natAbs_mul, Int.natAbs_mul, Int.natAbs_natCast]
  rcases ha with rfl | rfl <;> rcases he with rfl | rfl <;> simp

lemma cf_zero (b : ℤ) (w : List ℤ) (n : ℕ) :
    cf (0 :: b :: w) n * ((n : ℝ) + 1) = cf (b :: w) n := by
  unfold cf
  rw [dd_zero]
  simp only
  rw [pow_succ]
  have : ((n : ℝ) + 1) ≠ 0 := by positivity
  field_simp
  ring

lemma key_nz (a b : ℤ) (w : List ℤ) (ha : a = 1 ∨ a = -1) (hF : Facts (b :: w)) (t : ℝ)
    (n : ℕ) :
    cf (a :: b :: w) (n + 1) * (((n + 1 : ℕ) : ℝ) + 1) * t ^ (n + 1) -
      a * t * (cf (a :: b :: w) n * ((n : ℝ) + 1) * t ^ n) = cf (b :: w) n * t ^ (n + 1) := by
  have ha0 : a ≠ 0 := by rcases ha with rfl | rfl <;> norm_num
  obtain ⟨he, hk, _, _, _⟩ := hF
  obtain ⟨hsg, hab⟩ := sign_u a (dd (b :: w)).1 (dd (b :: w)).2.1 ha he hk
  unfold cf
  rw [dd_nz a b w ha0]
  simp only
  rw [Z_cons_succ, gq, hsg, hab]
  set e := (dd (b :: w)).1
  set k := (dd (b :: w)).2.1
  set C := (dd (b :: w)).2.2.1
  set v := (dd (b :: w)).2.2.2
  set u := a * e * (k : ℤ)
  push_cast
  have h1 : ((n : ℝ) + 1) ≠ 0 := by positivity
  have h2 : ((n : ℝ) + 1 + 1) ≠ 0 := by positivity
  have hsq : ((-1 : ℝ) ^ (n * 2)) = 1 := by
    rw [mul_comm, pow_mul]; norm_num
  field_simp
  rcases ha with rfl | rfl
  · push_cast; ring
  · push_cast
    simp only [mul_pow, pow_succ]
    linear_combination ((C : ℝ) * e * (e : ℝ) ^ n * (Z v n : ℝ) * t * t ^ n) * hsq

lemma dict (w : List ℤ) : w ≠ [] → IsLetters w → Regular w →
    Facts w ∧ ∀ x : ℝ, |x| < 1 → HasSum (fun N : ℕ => cf w N * x ^ (N + 1)) (H w x) := by
  induction w with
  | nil => intro h; exact absurd rfl h
  | cons a w ih =>
    intro _ hl hr
    have hw' : IsLetters w := fun b hb => hl b (by simp [hb])
    have hr' := reg_tail a w hr
    cases w with
    | nil =>
      have ha : a ≠ 0 := hr a (by simp)
      rcases hl a (by simp) with h | h | h
      · exact absurd h ha
      · subst h
        exact ⟨by simp [Facts, dd, IsIndex, wt], log_series_p1⟩
      · subst h
        exact ⟨by simp [Facts, dd, IsIndex, wt], log_series_m1⟩
    | cons b w' =>
      obtain ⟨hF, hS⟩ := ih (by simp) hw' hr'
      by_cases ha : a = 0
      · subst ha
        have hF0 : Facts (0 :: b :: w') := by
          obtain ⟨h1, h2, h3, h4, h5⟩ := hF
          refine ⟨?_, ?_, ?_, ?_, ?_⟩ <;> rw [dd_zero] <;> simp only
          · exact h1
          · omega
          · exact h3
          · exact h4
          · rw [List.length_cons, h5]; ring
        refine ⟨hF0, fun x hx => ?_⟩
        rw [H_cons]
        refine ftc_series (cf (0 :: b :: w')) _ (cf_bound _ hF0) _ x hx
          ((integrand_ok 0 (b :: w') hl hr).1 x hx) (fun t ht ht0 => ?_)
        have ht1 : |t| < 1 := by linarith
        have h := (hS t ht1).div_const t
        have h' : HasSum (fun n => cf (0 :: b :: w') n * ((n : ℝ) + 1) * t ^ n)
            (H (b :: w') t / t) := by
          refine h.congr_fun (fun n => ?_)
          rw [cf_zero, pow_succ]; field_simp
        rw [h'.tsum_eq]
        simp [f]; ring
      · have ha' : a = 1 ∨ a = -1 := by
          rcases hl a (by simp) with h | h | h
          · exact absurd h ha
          · exact Or.inl h
          · exact Or.inr h
        obtain ⟨he, hk, hc, hv, hlen⟩ := hF
        obtain ⟨hsg, hab⟩ := sign_u a (dd (b :: w')).1 (dd (b :: w')).2.1 ha' he hk
        have hF1 : Facts (a :: b :: w') := by
          refine ⟨?_, ?_, ?_, ?_, ?_⟩ <;> rw [dd_nz a b w' ha] <;> simp only
          · exact ha'
          · rcases ha' with rfl | rfl <;> rcases hc with h | h <;> rw [h] <;> simp
          · intro y hy
            simp only [List.mem_cons] at hy
            rcases hy with rfl | hy
            · intro h0; rw [h0] at hab; simp at hab; omega
            · exact hv y hy
          · rw [wt_cons, hab, List.length_cons, hlen]; ring
        refine ⟨hF1, fun x hx => ?_⟩
        rw [H_cons]
        refine ftc_series (cf (a :: b :: w')) _ (cf_bound _ hF1) _ x hx
          ((integrand_ok a (b :: w') hl hr).1 x hx) (fun t ht _ => ?_)
        have ht1 : |t| < 1 := by linarith
        set s := fun n : ℕ => cf (a :: b :: w') n * ((n : ℝ) + 1) * t ^ n with hsdef
        have hsum : Summable s := (psD _ _ (cf_bound _ hF1) t ht1).2.1
        have hs0 : s 0 = 0 := by
          simp only [hsdef, cf]
          rw [dd_nz a b w' ha]
          simp [Z_cons]
        have h1 : HasSum (fun n => s (n + 1)) (∑' n, s n) := by
          have := (hasSum_nat_add_iff' 1).mpr hsum.hasSum
          simpa [hs0] using this
        have h2 := hsum.hasSum.mul_left (a * t)
        have h3 := h1.sub h2
        have h4 : (fun n => s (n + 1) - a * t * s n) = fun n => cf (b :: w') n * t ^ (n + 1) :=
          funext (fun n => key_nz a b w' ha' ⟨he, hk, hc, hv, hlen⟩ t n)
        rw [h4] at h3
        have h5 := h3.unique (hS t ht1)
        show ∑' n, s n = _
        rw [abs_lt] at ht1
        rcases ha' with rfl | rfl
        · simp only [f, one_ne_zero, if_false, if_true]
          have : 1 - t ≠ 0 := ne_of_gt (by linarith)
          field_simp
          push_cast at h5
          linarith
        · simp only [f, show (-1 : ℤ) ≠ 0 by norm_num, show (-1 : ℤ) ≠ 1 by norm_num, if_false]
          have : 1 + t ≠ 0 := ne_of_gt (by linarith)
          field_simp
          push_cast at h5
          linarith

end S03

namespace S03

/-- The inverse of `w ↦ (ε, κ, v_w)`. -/
def build : ℤ → ℕ → List ℤ → List ℤ
  | e, k, [] => List.replicate (k - 1) 0 ++ [e]
  | e, k, u :: v => List.replicate (k - 1) 0 ++ e :: build (e * Int.sign u) u.natAbs v

lemma dd_zeros (m : ℕ) (w : List ℤ) (hw : w ≠ []) :
    dd (List.replicate m 0 ++ w) = ((dd w).1, (dd w).2.1 + m, (dd w).2.2.1, (dd w).2.2.2) := by
  induction m with
  | zero => simp
  | succ m ih =>
    rw [List.replicate_succ, List.cons_append]
    obtain ⟨b, r, hbr⟩ : ∃ b r, List.replicate m 0 ++ w = b :: r := by
      cases h : List.replicate m 0 ++ w with
      | nil => simp at h; exact absurd h.2 hw
      | cons b r => exact ⟨b, r, rfl⟩
    rw [hbr, dd_zero, ← hbr, ih]
    simp only; ring_nf

lemma reg_append (l w : List ℤ) (hw : w ≠ []) (hr : Regular w) : Regular (l ++ w) := by
  intro a ha
  apply hr a
  rw [List.getLast?_append] at ha
  cases h : w.getLast? with
  | none => simp [List.getLast?_eq_none_iff] at h; exact absurd h hw
  | some b => rw [h] at ha; simpa using ha

lemma build_props (v : List ℤ) (hv : IsIndex v) : ∀ (e : ℤ) (k : ℕ), (e = 1 ∨ e = -1) → 1 ≤ k →
    build e k v ≠ [] ∧ IsLetters (build e k v) ∧ Regular (build e k v) ∧
    (dd (build e k v)).1 = e ∧ (dd (build e k v)).2.1 = k ∧ (dd (build e k v)).2.2.2 = v := by
  induction v with
  | nil =>
    intro e k he hk
    have hne : [e] ≠ [] := by simp
    refine ⟨by simp [build], ?_, ?_, ?_⟩
    · intro a ha
      simp only [build, List.mem_append, List.mem_replicate, List.mem_singleton] at ha
      rcases ha with ⟨_, rfl⟩ | rfl
      · left; rfl
      · rcases he with rfl | rfl <;> simp
    · apply reg_append _ _ hne
      intro a ha; simp at ha; subst ha; rcases he with rfl | rfl <;> norm_num
    · simp only [build]
      rw [dd_zeros _ _ hne]
      simp only [dd]
      refine ⟨by simp, by omega, by simp⟩
  | cons u v ih =>
    intro e k he hk
    have hu : u ≠ 0 := hv u (by simp)
    have hv' : IsIndex v := fun b hb => hv b (by simp [hb])
    have he' : e * Int.sign u = 1 ∨ e * Int.sign u = -1 := by
      rcases lt_or_gt_of_ne hu with h | h
      · rw [Int.sign_eq_neg_one_of_neg h]; rcases he with rfl | rfl <;> simp
      · rw [Int.sign_eq_one_of_pos h]; rcases he with rfl | rfl <;> simp
    obtain ⟨h1, h2, h3, h4, h5, h6⟩ := ih hv' (e * Int.sign u) u.natAbs he'
      (Int.natAbs_pos.mpr hu)
    set w' := build (e * Int.sign u) u.natAbs v
    obtain ⟨b, r, hbr⟩ : ∃ b r, w' = b :: r := by
      cases h : w' with
      | nil => exact absurd h h1
      | cons b r => exact ⟨b, r, rfl⟩
    have hne : e :: w' ≠ [] := by simp
    have he0 : e ≠ 0 := by rcases he with rfl | rfl <;> norm_num
    refine ⟨by simp [build], ?_, ?_, ?_⟩
    · intro a ha
      simp only [build, List.mem_append, List.mem_replicate, List.mem_cons] at ha
      rcases ha with ⟨_, rfl⟩ | rfl | ha
      · left; rfl
      · rcases he with rfl | rfl <;> simp
      · exact h2 a ha
    · simp only [build]
      apply reg_append _ _ hne
      rw [hbr]; rw [hbr] at h3
      intro a ha; apply h3 a; simpa [List.getLast?_cons_cons] using ha
    · simp only [build]
      rw [dd_zeros _ _ hne]
      simp only
      rw [hbr, dd_nz e b r he0, ← hbr, h4, h5, h6]
      simp only
      refine ⟨by simp, by omega, ?_⟩
      congr 1
      rw [← mul_assoc, show e * e = 1 by rcases he with rfl | rfl <;> norm_num, one_mul,
        Int.sign_mul_natAbs]

lemma build_dd (w : List ℤ) : w ≠ [] → IsLetters w → Regular w →
    build (dd w).1 (dd w).2.1 (dd w).2.2.2 = w := by
  induction w with
  | nil => intro h; exact absurd rfl h
  | cons a w ih =>
    intro _ hl hr
    have hw' : IsLetters w := fun b hb => hl b (by simp [hb])
    have hr' := reg_tail a w hr
    cases w with
    | nil => simp [dd, build]
    | cons b w' =>
      have ihw := ih (by simp) hw' hr'
      obtain ⟨_, hk, _, _, _⟩ := (dict (b :: w') (by simp) hw' hr').1
      by_cases ha : a = 0
      · subst ha
        rw [dd_zero]
        simp only
        conv_rhs => rw [← ihw]
        obtain ⟨k, hk'⟩ : ∃ k, (dd (b :: w')).2.1 = k + 1 := ⟨_, (Nat.succ_pred_eq_of_pos hk).symm⟩
        rw [hk']
        cases h : (dd (b :: w')).2.2.2 <;> simp [build, List.replicate_succ]
      · have ha' : a = 1 ∨ a = -1 := by
          rcases hl a (by simp) with h | h | h
          · exact absurd h ha
          · exact Or.inl h
          · exact Or.inr h
        obtain ⟨he, _, _, _, _⟩ := (dict (b :: w') (by simp) hw' hr').1
        obtain ⟨hsg, hab⟩ := sign_u a (dd (b :: w')).1 (dd (b :: w')).2.1 ha' he hk
        rw [dd_nz a b w' ha]
        simp only [build, Nat.sub_self, List.replicate_zero, List.nil_append]
        rw [hsg, hab, ← mul_assoc, show a * a = 1 by rcases ha' with rfl | rfl <;> norm_num,
          one_mul, ihw]

theorem proof5 : stmt_5 := by
  refine ⟨fun a w hl hr x hx => (integrand_ok a w hl hr).1 x hx, fun w hne hl hr => ?_, ?_⟩
  · obtain ⟨⟨h1, h2, h3, h4, h5⟩, hS⟩ := dict w hne hl hr
    exact ⟨h1, h2, h3, h4, h5, hS⟩
  · intro e k v he hk hv
    obtain ⟨h1, h2, h3, h4, h5, h6⟩ := build_props v hv e k he hk
    refine ⟨build e k v, ⟨h1, h2, h3, h4, h5, h6⟩, ?_⟩
    rintro w ⟨hne, hl, hr, he', hk', hv'⟩
    rw [← build_dd w hne hl hr, he', hk', hv']

end S03

theorem S03.proof_5 : S03.stmt_5 := S03.proof5

namespace S03
open Real MeasureTheory Set Polynomial

lemma hd_pow (n : ℕ) (x : ℝ) (hx : x ≠ 0) : HasDerivAt (fun y => y ^ n) (n * x ^ n / x) x := by
  convert hasDerivAt_pow n x using 1
  cases n with
  | zero => simp
  | succ n => field_simp; push_cast; ring

lemma hd_sub (n : ℕ) (x : ℝ) (hx : 1 - x ≠ 0) :
    HasDerivAt (fun y => (1 - y) ^ n) (-(n * (1 - x) ^ n / (1 - x))) x := by
  have h := (hasDerivAt_pow n (1 - x)).comp x ((hasDerivAt_id' x).const_sub 1)
  refine h.congr_deriv ?_
  cases n with
  | zero => simp
  | succ n => field_simp; push_cast; ring

lemma hd_add (n : ℕ) (x : ℝ) (hx : 1 + x ≠ 0) :
    HasDerivAt (fun y => (1 + y) ^ n) (n * (1 + x) ^ n / (1 + x)) x := by
  have h := (hasDerivAt_pow n (1 + x)).comp x ((hasDerivAt_id' x).const_add 1)
  refine h.congr_deriv ?_
  cases n with
  | zero => simp
  | succ n => field_simp; push_cast; ring

/-- `x d/dx` of `p/(x^a(1-x)^b(1+x)^d)` has numerator `θP`, denominators `a, b+1, d+1`. -/
noncomputable def thP (p : Polynomial ℚ) (a b d : ℕ) : Polynomial ℚ :=
  X * derivative p * (1 - X) * (1 + X) - C (a : ℚ) * p * (1 - X) * (1 + X)
    + C (b : ℚ) * X * p * (1 + X) - C (d : ℚ) * X * p * (1 - X)

lemma kval_deriv (p : Polynomial ℚ) (a b d : ℕ) (x : ℝ) (hx : x ≠ 0) (h1 : 1 - x ≠ 0)
    (h2 : 1 + x ≠ 0) :
    HasDerivAt (kval p a b d) (kval (thP p a b d) a (b + 1) (d + 1) x / x) x := by
  have hD : HasDerivAt (fun y => y ^ a * (1 - y) ^ b * (1 + y) ^ d) _ x :=
    ((hd_pow a x hx).mul (hd_sub b x h1)).mul (hd_add d x h2)
  have hN := Polynomial.hasDerivAt_aeval p x
  have hD0 : x ^ a * (1 - x) ^ b * (1 + x) ^ d ≠ 0 := by
    apply mul_ne_zero (mul_ne_zero (pow_ne_zero _ hx) (pow_ne_zero _ h1)) (pow_ne_zero _ h2)
  have := hN.div hD hD0
  show HasDerivAt (fun y => aeval y p / (y ^ a * (1 - y) ^ b * (1 + y) ^ d)) _ x
  convert this using 1
  unfold kval
  simp only [thP, map_sub, map_mul, map_add, aeval_X, aeval_C, map_one, eq_ratCast,
    Rat.cast_natCast, Pi.mul_apply]
  field_simp
  ring

/-- The terms of `θ(r H_w)`. -/
noncomputable def thT (g : GTerm) : List GTerm :=
  (thP g.1 g.2.1 g.2.2.1 g.2.2.2.1, g.2.1, g.2.2.1 + 1, g.2.2.2.1 + 1, g.2.2.2.2) ::
    (match g.2.2.2.2 with
      | [] => []
      | c :: w' => if c = 0 then [(g.1, g.2.1, g.2.2.1, g.2.2.2.1, w')]
          else if c = 1 then [(X * g.1, g.2.1, g.2.2.1 + 1, g.2.2.2.1, w')]
          else [(X * g.1, g.2.1, g.2.2.1, g.2.2.2.1 + 1, w')])

noncomputable def gv (x : ℝ) (g : GTerm) : ℝ :=
  kval g.1 g.2.1 g.2.2.1 g.2.2.2.1 x * H g.2.2.2.2 x

lemma thT_deriv (g : GTerm) (hw : IsLetters g.2.2.2.2) (hr : Regular g.2.2.2.2) (x : ℝ)
    (hx0 : 0 < |x|) (hx1 : |x| < 1) :
    HasDerivAt (fun y => gv y g) (((thT g).map (gv x)).sum / x) x := by
  obtain ⟨p, a, b, d, w⟩ := g
  simp only at hw hr
  have hx : x ≠ 0 := abs_pos.mp hx0
  have h1 : 1 - x ≠ 0 := by rw [abs_lt] at hx1; exact ne_of_gt (by linarith)
  have h2 : 1 + x ≠ 0 := by rw [abs_lt] at hx1; exact ne_of_gt (by linarith)
  have hk := kval_deriv p a b d x hx h1 h2
  cases w with
  | nil =>
    simp only [gv, thT, H, mul_one, List.map_cons, List.map_nil, List.sum_cons, List.sum_nil,
      add_zero]
    exact hk
  | cons c w' =>
    have hH := H_deriv c w' hw hr x hx0 hx1
    have := hk.mul hH
    simp only [gv] at this ⊢
    convert this using 1
    have hc := hw c (by simp)
    simp only [thT, List.map_cons, List.map_nil, List.sum_cons, List.sum_nil, add_zero]
    rcases hc with rfl | rfl | rfl
    · simp only [if_true, List.map_cons, List.map_nil, List.sum_cons, List.sum_nil, add_zero, gv]
      simp only [f, if_true]
      field_simp
    · simp only [one_ne_zero, if_false, if_true, List.map_cons, List.map_nil, List.sum_cons,
        List.sum_nil, add_zero, gv, kval, map_mul, aeval_X]
      simp only [f, one_ne_zero, if_false, if_true]
      field_simp
      ring
    · simp only [show (-1 : ℤ) ≠ 0 by norm_num, show (-1 : ℤ) ≠ 1 by norm_num, if_false,
        List.map_cons, List.map_nil, List.sum_cons, List.sum_nil, add_zero, gv, kval, map_mul,
        aeval_X]
      simp only [f, show (-1 : ℤ) ≠ 0 by norm_num, show (-1 : ℤ) ≠ 1 by norm_num, if_false]
      field_simp
      ring

end S03

namespace S03
open Real MeasureTheory Set Polynomial

noncomputable def negT (g : GTerm) : GTerm :=
  (C ((-1 : ℚ) ^ (g.2.1 + (g.2.2.2.2.filter fun a => a ≠ 0).length)) * g.1.comp (-X),
    g.2.1, g.2.2.2.1, g.2.2.1, g.2.2.2.2.map fun a => -a)

lemma gv_neg (g : GTerm) (hw : IsLetters g.2.2.2.2) (x : ℝ) :
    gv (-x) g = gv x (negT g) := by
  obtain ⟨p, a, b, d, w⟩ := g
  simp only at hw
  simp only [gv, negT, kval, H_neg w hw x, map_mul, aeval_C, aeval_comp, map_neg, aeval_X,
    eq_ratCast, Rat.cast_pow, Rat.cast_neg, Rat.cast_one, sub_neg_eq_add, neg_pow x a, pow_add]
  rcases neg_one_pow_eq_or ℝ a with h | h <;> rw [h] <;> ring

lemma InG_neg (B : ℕ) (F : ℝ → ℝ) (hF : InG B F) : InG B (fun x => F (-x)) := by
  obtain ⟨L, hL, hFL⟩ := hF
  refine ⟨L.map negT, ?_, ?_⟩
  · intro g hg
    rw [List.mem_map] at hg
    obtain ⟨g0, hg0, rfl⟩ := hg
    obtain ⟨h1, h2, h3⟩ := hL g0 hg0
    refine ⟨?_, ?_, ?_⟩
    · intro c hc
      simp only [negT, List.mem_map] at hc
      obtain ⟨c0, hc0, rfl⟩ := hc
      rcases h1 c0 hc0 with h | h | h <;> rw [h] <;> norm_num
    · intro c hc
      simp only [negT, List.getLast?_map, Option.mem_def, Option.map_eq_some_iff] at hc
      obtain ⟨c0, hc0, rfl⟩ := hc
      exact neg_ne_zero.mpr (h2 c0 hc0)
    · simpa [negT] using h3
  · intro x hx0 hx1
    have := hFL (-x) (by rwa [abs_neg]) (by rwa [abs_neg])
    show F (-x) = _
    rw [this, List.map_map]
    show (L.map (gv (-x))).sum = _
    congr 1
    refine List.map_congr_left (fun g hg => ?_)
    exact gv_neg g (hL g hg).1 x

lemma sum_deriv (x : ℝ) (hx0 : 0 < |x|) (hx1 : |x| < 1) : ∀ L : List GTerm,
    (∀ g ∈ L, IsLetters g.2.2.2.2 ∧ Regular g.2.2.2.2) →
    HasDerivAt (fun y => (L.map (gv y)).sum) (((L.flatMap thT).map (gv x)).sum / x) x := by
  intro L
  induction L with
  | nil => intro _; simpa using hasDerivAt_const x (0 : ℝ)
  | cons g L ih =>
    intro hL
    have h1 := thT_deriv g (hL g (by simp)).1 (hL g (by simp)).2 x hx0 hx1
    have h2 := ih (fun g' hg' => hL g' (by simp [hg']))
    have := h1.add h2
    simp only [List.map_cons, List.sum_cons, List.flatMap_cons, List.map_append,
      List.sum_append, add_div]
    exact this

lemma InG_theta (B : ℕ) (F : ℝ → ℝ) (hF : InG B F) :
    ∃ G : ℝ → ℝ, InG B G ∧ ∀ x : ℝ, 0 < |x| → |x| < 1 → HasDerivAt F (G x / x) x := by
  obtain ⟨L, hL, hFL⟩ := hF
  refine ⟨fun x => ((L.flatMap thT).map (gv x)).sum, ⟨L.flatMap thT, ?_, fun x _ _ => rfl⟩, ?_⟩
  · intro g' hg'
    rw [List.mem_flatMap] at hg'
    obtain ⟨g, hg, hg'⟩ := hg'
    obtain ⟨h1, h2, h3⟩ := hL g hg
    obtain ⟨p, a, b, d, w⟩ := g
    simp only at h1 h2 h3
    simp only [thT, List.mem_cons] at hg'
    rcases hg' with rfl | hg'
    · exact ⟨h1, h2, h3⟩
    · cases w with
      | nil => simp at hg'
      | cons c w' =>
        have h1' : IsLetters w' := fun b hb => h1 b (by simp [hb])
        have h2' := reg_tail c w' h2
        have h3' : w'.length ≤ B := by simp at h3; omega
        simp only at hg'
        split_ifs at hg' <;> (simp only [List.mem_singleton] at hg'; subst hg'; exact ⟨h1', h2', h3'⟩)
  · intro x hx0 hx1
    have hd := sum_deriv x hx0 hx1 L (fun g hg => ⟨(hL g hg).1, (hL g hg).2.1⟩)
    refine hd.congr_of_eventuallyEq ?_
    have hopen : IsOpen (Ioo (-1 : ℝ) 1 \ {0}) := isOpen_Ioo.sdiff isClosed_singleton
    have hxs : x ∈ Ioo (-1 : ℝ) 1 \ {0} := ⟨abs_lt.mp hx1, by simpa [abs_pos] using hx0⟩
    filter_upwards [hopen.mem_nhds hxs] with y hy
    obtain ⟨hy1, hy0⟩ := hy
    have hy0' : y ≠ 0 := hy0
    exact hFL y (abs_pos.mpr hy0') (abs_lt.mpr hy1)

theorem proof6 : stmt_6 :=
  ⟨fun w _ hw _ x _ => H_neg w hw x,
   fun a w hw hr x hx0 hx1 => H_deriv a w hw hr x hx0 hx1,
   fun B F hF => ⟨InG_neg B F hF, InG_theta B F hF⟩⟩

end S03

theorem S03.proof_6 : S03.stmt_6 := S03.proof6
