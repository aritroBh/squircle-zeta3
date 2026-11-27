import Mathlib
import Zeta3.Wave3.W02_Holo
import Zeta3.Wave5.Y01_M2
/-! Statements stmt_3, stmt_4, stmt_16, stmt_17 of this audited block are in `P06_WeightTwoRest.lean` (target P06); the header below still describes all of them. -/
/-!
# P01 -- the weight-two part of every `M_k` (Theorem `thm:wt2`)

Target id: P01.
Paper labels: `thm:wt2` (main.tex 1022-1061: the theorem 1022-1024, its proof sketch 1026-1061);
supplement paper/ieee/proofs/sec_g1_wt2.tex, all of it: Theorem `g1:thm:wt2` (lines 10-14),
the reduction (lines 16-42, Definition `g1:def:red`), Lemmas `g1:lem:init` (46-55), `g1:lem:A`
(67-69), `g1:lem:B` (79-83), `g1:lem:I` (98-105), `g1:lem:D` (115-117), `g1:lem:kernel`
(129-131), the proof of the theorem (139-141), and the second half: `g1:lem:inv` (156-168),
`g1:lem:par` (177-179), `g1:lem:res` (197-203), `g1:lem:real` (209-215), `g1:lem:split`
(229-235). Notation from sec_g1_M2.tex (`g1:lem:Mk`, `g1:lem:euler`, folding, `g1:lem:k2`,
`g1:lem:T`).

Quoted verbatim (main.tex 1022-1024):

  \begin{theorem}[The weight-two part of every $M_k$$^\dagger$]\label{thm:wt2}
  For every $k\ge2$, the weight-two part of $M_k$ (in the sense of the reduction of
  Section~\ref{sec:method}) is exactly $T=(\pi^2+6\lambda^2)/24$.
  \end{theorem}

Quoted verbatim (supplement, after g1:thm:wt2):

  the theorem then says: \emph{the degree-one kernel produced by the reduction at level $k$ is
  $\Rwt$, and $\int_b^a\Rwt\,\ell_1\,dv=T$}.

Quoted verbatim (supplement, the list of symbols, line 7):

  $a=1+\sqrt2$, $b=\sqrt2-1$, $L$, $\lambda$, $t(v)=\frac{(v-b)(v+a)}{4v}$, $\ell_1=\log t$,
  $\ell_2=\log(1-t)$, $D_j=\ell_j'$ (Lemma~\ref{g1:lem:euler}), $\Fold(R)(v)=R(1/v)/v^2$,
  $\Rwt=\frac1{2v}+\frac1{2(v-a)}-\frac1{v+a}-\frac1{2(v+b)}$, and $M_k$ the Taylor coefficients
  of Lemma~\ref{g1:lem:Mk}.

Quoted verbatim (supplement, "The reduction, made precise"):

  $D_1$, $D_2$ have simple zeros at $\pm i$ (the numerator $v^2+1$), and
  D_1'(i)=\frac{2}{(i-b)(i+a)}=\tau_1,\quad D_2'(i)=\frac2{(i-a)(i+b)}=\tau_2,\qquad
  D_1'(-i)=\tau_2,\quad D_2'(-i)=\tau_1 ,

  \paragraph{Initial data.} By Lemma~\ref{g1:lem:Mk},
  $-\frac12\varphi_k\frac{dt}{dv}=\sum_{r+s=k}\tilde R^{(k)}_{rs}\,\ell_1^r\ell_2^s$ with rational
  $\tilde R^{(k)}_{rs}\in\QQ(\sqrt2)(v)$. Fold to canonical keys:
  $R^{(k)}_{rs}=\tilde R^{(k)}_{rs}+\Fold(\tilde R^{(k)}_{sr})$ for $r>s$,
  $R^{(k)}_{rr}=\tilde R^{(k)}_{rr}$. Then
  $M_k=\int_b^a\sum_{r\ge s,\ r+s=k}R^{(k)}_{rs}\ell_1^r\ell_2^s\,dv$ (folding rule of
  Section~\ref{g1:sub:M2:int}). Put $R^{(k)}_{rs}=0$ if $r+s\ne k$.

  \paragraph{The operators $\Hp$ and $\Jp$.} For a rational function $W$ let $\Hp(W)$ be the sum
  of the principal parts of $W$ at $v=i$ and $v=-i$ of order $\ge2$ (the terms $c_m(v\mp i)^{-m}$,
  $m\ge2$). It has real coefficients and no residue, hence has a rational primitive; let
  $\Jp(W)$ be the primitive of $\Hp(W)$ that vanishes at $v=b$.

  \begin{definition}[the reduction]\label{g1:def:red}
  For $k\ge1$ define rational functions $W^{(k)}_{rs}$ and $G^{(k)}_{rs}$ for all keys, by
  descending induction on $d=r+s$ (starting with $W^{(k)}_{rs}=G^{(k)}_{rs}=0$ for $d>k$):
  W^{(k)}_{rs}&=R^{(k)}_{rs}-(r+1)G^{(k)}_{r+1,s}D_1-[r>s]\,(s+1)\,G^{(k)}_{r,s+1}D_2
    -[s=r-1]\,r\,\Fold\bigl(G^{(k)}_{rr}D_1\bigr),\\
  G^{(k)}_{rs}&=\Jp\bigl(W^{(k)}_{rs}\bigr)\ \ (d\ge2),\qquad G^{(k)}_{rs}=0\ \ (d=1),
  where $[\cdot]$ is an indicator and a term whose index is not a key is $0$.
  \end{definition}

  The \emph{weight-two part of $M_k$} is
  W_2(k)=\int_b^aK^{(k)}_{10}(v)\,\ell_1\,dv,
  where $K^{(k)}_{10}$ is the simple-pole part of $W^{(k)}_{10}$ at the real letters $0,\pm a,\pm b$

Quoted verbatim (Lemma g1:lem:init):

  Let $d=r+s=k\ge2$. Put $V=(1-t)^2(e^{2x}-1)+t^2(e^{2y}-1)$ and $U_j^{(r,s)}=[x^ry^s]V^j$, a
  polynomial in $t$. Then
  \tilde R^{(k)}_{rs}=\frac1{2\sqrt2\,v\,t(1-t)}\sum_{j=1}^d\binom{1/2}{j}q^{2-2j}\,U_j^{(r,s)}+E_{rs},
  where $E_{k0}=-\frac{D_1}{2\,k!}$, $E_{0k}=\frac{D_2}{2\,k!}$ and $E_{rs}=0$ otherwise; $E_{rs}$
  has no pole at $\pm i$. At $v=\pm i$ the functions $\tilde R_{rs}$ and $R_{rs}=R^{(k)}_{rs}$ have
  poles of order at most $2d-2$, and the coefficient of $(v-i)^{-(2d-2)}$ in $R_{rs}$ is
  \rho_{rs}=K_d\,(-1)^r\bigl(1+[r>s]\bigr)\binom dr,\qquad
  K_d=\frac{(-1)^d(2d-3)!!\;i^{d+1}}{2\sqrt2\;d!}.

Quoted verbatim (Lemma g1:lem:A):

  For every $k\ge1$ and every key $(r,s)$ of degree $d$, $W^{(k)}_{rs}$ has poles at $\pm i$ of
  order at most $2d-2$ (so for $d=1$ it has none there), and for $d\ge2$, $G^{(k)}_{rs}$ has poles
  at $\pm i$ of order at most $2d-3$ and no other poles.

Quoted verbatim (before and in Lemma g1:lem:B):

  For $d\ge2$ let $h^{(k)}_{rs}$ be the coefficient of $(v-i)^{-(2d-2)}$ in $W^{(k)}_{rs}$
  ... $h^{(k)}_{rs}=0$ if $d>k$.
  For $k\ge d\ge2$,
  (2d-1)\bigl(h^{(k)}_{rs}-\delta_{kd}\rho_{rs}\bigr)=(r+1)\tau_1h^{(k)}_{r+1,s}
   +[r\ge s+2]\,(s+1)\tau_2h^{(k)}_{r,s+1}+[s=r-1]\,r\,\tau_2\bigl(h^{(k)}_{rr}-\overline{h^{(k)}_{rr}}\bigr).

Quoted verbatim (Lemma g1:lem:I):

  Write $\rho_{rs}=K_d(-1)^r(1+[r>s])\binom dr$ and $\hat h_p=(p+1)\tau_1\rho_{p+1,p}/(4p-1)$.
  For all $d\ge2$, $p\ge1$:
  \text{(I1)}\quad&(d+1)K_{d+1}=-i(2d-1)K_d,\\
  \text{(I2)}\quad&(2d-1)\rho_{rs}=(r+1)\tau_1\rho_{r+1,s}+(s+1)\tau_2\rho_{r,s+1}\qquad(r+s=d,\ r\ge s+2),\\
  \text{(I3)}\quad&(4p-3)\rho_{p,p-1}=(p+1)\tau_1\rho_{p+1,p-1}+p\,\tau_2\bigl(\rho_{pp}-\overline{\rho_{pp}}\bigr)\qquad(p\ge2),\\
  \text{(I4)}\quad&\hat h_p=(1-i)\rho_{pp},\quad \rho_{pp}\in i\RR,\quad \hat h_p-\overline{\hat h_p}=\rho_{pp}-\overline{\rho_{pp}}.

Quoted verbatim (Lemma g1:lem:D):

  For every $k\ge r+s=d\ge2$: if $r>s$ then $h^{(k)}_{rs}=\rho_{rs}$; on the diagonal
  $h^{(2p)}_{pp}=\rho_{pp}$ and $h^{(k)}_{pp}=\hat h_p$ for $k\ge2p+1$.

Quoted verbatim (Lemma g1:lem:kernel):

  For every $k\ge2$: $W^{(k)}_{10}=\Rwt$. More precisely, $G^{(k)}_{20}=g_{20}$ for all $k\ge2$,
  $G^{(2)}_{11}=-g_{20}$ and $G^{(k)}_{11}=g_{11}$ for $k\ge3$, with $g_{20}$ as in
  Lemma~\ref{g1:lem:k2} and $g_{11}=-\frac{\sqrt2(v^2+2v-1)}{4(v^2+1)}$.
  (g1:lem:k2: g_{20}(v)=\frac a4-\frac{\sqrt2}{2(v^2+1)}.)

Quoted verbatim (proof of g1:thm:wt2):

  By Lemma~\ref{g1:lem:kernel}, for every $k\ge2$ the degree-one kernel is
  $W^{(k)}_{10}=\Rwt=\frac1{2v}+\frac1{2(v-a)}-\frac1{v+a}-\frac1{2(v+b)}$. It has only simple
  poles, all at real letters, so it is its own simple-pole part $K^{(k)}_{10}$
  (Section~\ref{g1:sub:wt2:def}). Therefore $W_2(k)=\int_b^a\Rwt\,\ell_1\,dv=T=\frac{\pi^2+6\lambda^2}{24}$
  by Lemma~\ref{g1:lem:T}. \qed

Quoted verbatim (second half): $\iota(v)=-1/v$, $\mathcal I(W)(v)=\frac{W(-1/v)}{v^2}$;

  Lemma g1:lem:inv. For rational functions $W$, $G$:
  (a) $\mathcal I$ is linear and $\mathcal I(\mathcal I(W))=W$; $\mathcal I\circ\Fold=\Fold\circ\mathcal I$
  (b) $\mathcal I(D_1)=D_1$, $\mathcal I(D_2)=D_2$, and $\frac{d}{dv}\bigl[G(-1/v)\bigr]=\mathcal I(G')$.
  (c) If $W$ has no pole at $i$ (resp.\ $-i$), neither has $\mathcal I(W)$. For $m\ge1$,
  \mathcal I\bigl((v-i)^{-m}\bigr)=i^m\,\frac{v^{m-2}}{(v-i)^m},\qquad
  \mathcal I\bigl((v+i)^{-m}\bigr)=(-i)^m\,\frac{v^{m-2}}{(v+i)^m};
  (d) $\Hp(\mathcal I(W))=\mathcal I(\Hp(W))$, and $\Res_{v=\pm i}\mathcal I(W)=\Res_{v=\pm i}W$.

  Lemma g1:lem:par. Let $k\ge2$. For every key $(r,s)$ the function
  $W^{(k)}_{rs}+\mathcal I(W^{(k)}_{rs})$ has no pole at $\pm i$. For every key of degree $d\ge2$
  the function $G^{(k)}_{rs}(v)+G^{(k)}_{rs}(-1/v)$ is constant.

  Lemma g1:lem:res. For every $k\ge2$ and every key $(r,s)$, $\Res_{v=\pm i}W^{(k)}_{rs}=0$.
  Hence the \emph{simple kernel}
  S^{(k)}_{rs}:=W^{(k)}_{rs}-\Hp\bigl(W^{(k)}_{rs}\bigr)=W^{(k)}_{rs}-\bigl(G^{(k)}_{rs}\bigr)'\quad(d\ge2),
  \qquad S^{(k)}_{10}:=W^{(k)}_{10},
  has no pole at $\pm i$.

  Lemma g1:lem:real. For every $k\ge2$ and every key $(r,s)$ of degree $d\ge1$ there are
  $c_x=c^{(k)}_{rs,x}\in\Qs$ with
  S^{(k)}_{rs}=\sum_{x\in\{0,\,a,\,-a,\,-b\}}\frac{c_x}{v-x}.
  In particular $S^{(k)}_{rs}$ has no pole at $b$.

  Lemma g1:lem:split. For every $k\ge2$,
  M_k=\sum_{d=1}^k\ \sum_{\substack{r\ge s\ge0\\ r+s=d}}\int_b^aS^{(k)}_{rs}(v)\,\ell_1^r\ell_2^s\,dv,
  where every integral converges absolutely and
  $S^{(k)}_{rs}=\sum_{x\in\{0,a,-a,-b\}}c^{(k)}_{rs,x}/(v-x)$ (Lemma~\ref{g1:lem:real}). The term
  $d=1$ is $\int_b^a\Rwt\,\ell_1\,dv=T$.

## Conventions

* Rational functions in `v` are elements of `RatFunc ℂ` (complex coefficients; the paper's
  functions have coefficients in `ℚ(√2) ⊂ ℝ ⊂ ℂ`, and `ℚ(√2)` membership is stated separately in
  stmt_16). `vX` is the variable `v`, `cst c` the constant `c`. `√2 = s2`, `a = aC`, `b = bC`
  are the casts of `Real.sqrt 2`, `Y01.a = 1+√2`, `Y01.b = √2-1` (imported, statement block).
* Laurent coefficients: `lcoeff c f n` is the coefficient of `(v-c)^n` in `f`, defined as the
  `n`-th coefficient of the Laurent series (at `X = 0`) of `RatFunc.laurent c f = f(X + c)`.
  `poleLE c f m` ("pole of order at most `m` at `c`") says `lcoeff c f n = 0` for all `n < -m`;
  `poleLE c f 0` is "no pole at `c`". The residue at `c` is `lcoeff c f (-1)`.
* `ev f z` is `RatFunc.eval` at `z` (numerator and denominator in lowest terms); at a pole it is
  Lean's junk value `x/0 = 0`. It is used only at real `v ∈ (b,a)` (inside the integrals the
  endpoints are a null set), where none of the functions below has a pole, except in stmt_1
  where `z ≠ ±i` is assumed and `Jp`, `Hp` have poles only at `±i`.
* Composition: `subst g f = f.num(g) / f.denom(g)` is `f ∘ g` for `g ∈ RatFunc ℂ` not constant;
  `Fold f = f(1/v)/v²` and `Inv f = 𝓘(f) = f(-1/v)/v²`.
* `Hp f` is written out as `Σ_{c=±i} Σ_{m=2}^{N} lcoeff c f (-m) (v-c)^{-m}` with
  `N = deg(denominator of f)`, which bounds every pole order of `f`. `Jp f` is written out as the
  termwise primitive `Σ_{c=±i} Σ_{m=2}^{N} lcoeff c f (-m) · (-(1/(m-1))) ((v-c)^{-(m-1)} - (b-c)^{-(m-1)})`.
  stmt_1 states that these are the paper's objects (the principal parts of order `≥2` at `±i`,
  no residue, and `Jp f` the primitive of `Hp f` vanishing at `b`).
* `t(v)`, `q`, `D₁`, `D₂`, `R_wt`, `g₂₀`, `g₁₁` are the rational functions of the supplement.
* `U_j^{(r,s)} = [x^r y^s] V^j` is taken literally: `V` is a power series in `x` whose coefficients
  are power series in `y` over `RatFunc ℂ`, with `t` already replaced by `t(v)` (so `Ucoef j r s`
  is `U_j^{(r,s)}(t(v))`); `e^{2x} - 1` is `rescale 2 (exp) - 1`. `binom(1/2, j)` is
  `binHalf j = (Π_{m<j} (1/2 - m)) / j!`. `q^{2-2j}` is an integer power in the field `RatFunc ℂ`.
* `Rt r s` is `R̃_{rs}` given by the formula of g1:lem:init (with `d = r+s`); stmt_3 states that it
  is the coefficient of `ℓ₁^r ℓ₂^s` in `-½ φ_k dt/dv` (Lemma g1:lem:Mk), with
  `φ_k(t) = [ε^k] Φ_ε(t)` the `k`-th Taylor coefficient at `ε = 0` of `W02.gC ε t` (imported,
  statement block: `t^{ε-1} + (1-t)^{ε-1} - t^{ε-1}(1+z)^{1/2}`, which for real `ε` near `0` and any
  `t ∈ (0,1)` is the positive root, so it is `Φ_ε(t)` near `ε = 0`). `Rinit k r s` is `R^{(k)}_{rs}`
  (folded, `0` unless `r+s = k`, `r ≥ s`).
* The reduction: `W k r s` is `W^{(k)}_{rs}`, defined by recursion on `k+1-(r+s)` (descending in
  the degree); inside it the parent kernels `G^{(k)}_{parent}` are written `Jp (W k parent)`,
  which is the paper's `G` because a parent has degree `d+1 ≥ 2`; for `d > k`, `W = 0` and
  `Jp 0 = 0`. `G k r s` is `Jp (W k r s)` for `d ≥ 2` and `0` for `d ≤ 1`. Only keys
  `r ≥ s` are ever used (the `D₂` feed is guarded by `s < r`).
* `h k r s = lcoeff i (W k r s) (-(2d-2))`, `d = r+s`. `K_d`, `ρ_{rs}`, `τ₁`, `τ₂`, `ĥ_p` as in the
  supplement; `(2d-3)!!` is `Nat.doubleFactorial (2d-3)` (used only for `d ≥ 2`).
* `K^{(k)}_{10}` is `simplePart (W k 1 0)` `= Σ_{x∈{0,a,-a,b,-b}} Res_x(W)/(v-x)`, and
  `W₂(k) = ∫_b^a K^{(k)}_{10}(v) ℓ₁(v) dv` (`W2 k`), a complex interval integral of a function that is
  real on `(b,a)`; `ℓ₁ = Y01.l1`, `ℓ₂ = Y01.l2` (imported). `λ = Y01.lam = log 2`.
* `M_k` is `W02.Mk k = -F^{(k)}(0)/(2 k!)` (imported). No `ζ` value occurs in this file.
* Every integral whose value is asserted comes with an `IntervalIntegrable` claim.

## Checked before freezing

A sympy transcription of exactly these definitions (`Rt` from the formula with the literal
`[x^r y^s]V^j`, `Rinit`, `Fold`, `lcoeff` by Laurent expansion at `±i`, `Hp`/`Jp` with
`N = deg denom`, the recursion `W`) was run for `k = 2, 3, 4`: every key satisfies the pole bound
of stmt_6, the recursion of stmt_7, the values of stmt_9 (in particular `h^{(4)}_{21} = -3√2/4`),
residue `0` at `±i`; `G^{(k)}_{20} = g₂₀`, `G^{(2)}_{11} = -g₂₀`, `G^{(k)}_{11} = g₁₁` (`k = 3, 4`),
and `W^{(k)}_{10} - R_wt = 0` for `k = 2, 3, 4`. This agrees with the supplement's script
g1_wt2_check.py.

## Statements

* stmt_1: `Hp`, `Jp` are the paper's operators.
* stmt_2: `Fold(D₁) = -D₂` (g1:lem:euler); `D₁`, `D₂` vanish at `±i` with `D₁'(i) = τ₁`, ...
* stmt_3: g1:lem:init, first claim: `-½ φ_k dt/dv = Σ_{r+s=k} R̃_{rs} ℓ₁^r ℓ₂^s` on `(b,a)`, `k ≥ 2`.
* stmt_4: "Then `M_k = ∫_b^a Σ_{r≥s, r+s=k} R^{(k)}_{rs} ℓ₁^r ℓ₂^s dv`", `k ≥ 1`, with integrability.
* stmt_5: g1:lem:init, pole bound `2d-2` for `R̃` and `R`, leading coefficient `ρ_{rs}`.
* stmt_6: g1:lem:A.
* stmt_7: g1:lem:B (the corrected recursion with `[r ≥ s+2]`).
* stmt_8: g1:lem:I, (I1)-(I4).
* stmt_9: g1:lem:D.
* stmt_10: g1:lem:kernel, "more precisely": `G^{(k)}_{20}`, `G^{(k)}_{11}`.
* stmt_11: g1:lem:kernel: `W^{(k)}_{10} = R_wt` for every `k ≥ 2`.
* stmt_12: Theorem g1:thm:wt2 / thm:wt2 in the supplement's precise form, for every `k ≥ 2`.
* stmt_13: g1:lem:inv (a), (b) first part, (c) formulas, (d).
* stmt_14: g1:lem:par.
* stmt_15: g1:lem:res.
* stmt_16: g1:lem:real.
* stmt_17: g1:lem:split.

## NOT COVERED

* The value `∫_b^a R_wt ℓ₁ dv = (π²+6λ²)/24` (g1:lem:T) is not restated: it is `Y01.stmt_8`
  (`Y01.Tint`, proved in Zeta3.Wave5.Y01_M2). stmt_12 states the theorem's integral for the
  rational function `K^{(k)}_{10}` of this file.
* "It has real coefficients" (for `Hp(W)`), "$\tilde R^{(k)}_{rs}\in\QQ(\sqrt2)(v)$" and
  "$E_{rs}$ has no pole at $\pm i$" are not stated separately (`ℚ(√2)` membership is stated for
  the final coefficients `c_x` in stmt_16).
* The derivative identity of g1:lem:inv (b), `d/dv[G(-1/v)] = 𝓘(G')`, and the first sentence of
  (c) (no pole is created at `±i`) are not stated: `RatFunc` has no derivative in Mathlib, and both
  are proof steps of stmt_14.
* The `S = W - (G)'` form in g1:lem:res (only `S = W - Hp(W)` is stated; stmt_1 gives `G' = Hp(W)`
  pointwise).
* The paper's sentences on the sum over `k` (the generating-series identity
  `(1-ε)(W₁₀-εR₁₀) = ε²R` of the main-text sketch, and `G₂₀ = ε²g₂₀/(1-ε)` etc.) are the summed
  form of stmt_10, stmt_11 and are not stated separately. The remark on `k ≤ 15` computations, the
  remark that `k = 1` gives `M₁`, and the remark after g1:lem:split (branch points, the 1-form
  argument) are not formalized.
* "weight" is not formalized beyond its definition in the supplement (the integral `W₂(k)` of the
  simple-pole part of the degree-one kernel).
-/
-- BEGIN STATEMENTS
namespace P01

open MeasureTheory Real

/-- `√2` as a complex number. -/
noncomputable def s2 : ℂ := ((Real.sqrt 2 : ℝ) : ℂ)

/-- `a = 1 + √2` as a complex number (`Y01.a`). -/
noncomputable def aC : ℂ := ((Y01.a : ℝ) : ℂ)

/-- `b = √2 - 1` as a complex number (`Y01.b`). -/
noncomputable def bC : ℂ := ((Y01.b : ℝ) : ℂ)

/-- The variable `v` of the rational functions. -/
noncomputable def vX : RatFunc ℂ := RatFunc.X

/-- The constant rational function `c`. -/
noncomputable def cst (c : ℂ) : RatFunc ℂ := RatFunc.C c

/-- Evaluation of a rational function at `z` (numerator and denominator in lowest terms). -/
noncomputable def ev (f : RatFunc ℂ) (z : ℂ) : ℂ := RatFunc.eval (RingHom.id ℂ) z f

/-- Composition `f ∘ g` of rational functions: `f.num(g) / f.denom(g)`. -/
noncomputable def subst (g f : RatFunc ℂ) : RatFunc ℂ :=
  Polynomial.aeval g f.num / Polynomial.aeval g f.denom

/-- The fold `Fold(R)(v) = R(1/v)/v²`. -/
noncomputable def Fold (f : RatFunc ℂ) : RatFunc ℂ := subst vX⁻¹ f / vX ^ 2

/-- The involution operator `𝓘(W)(v) = W(-1/v)/v²` (g1:sub:wt2:split). -/
noncomputable def Inv (f : RatFunc ℂ) : RatFunc ℂ := subst (-vX⁻¹) f / vX ^ 2

/-- `lcoeff c f n`: the coefficient of `(v-c)^n` in the Laurent expansion of `f` at `v = c`. -/
noncomputable def lcoeff (c : ℂ) (f : RatFunc ℂ) (n : ℤ) : ℂ :=
  ((RatFunc.laurent c f : RatFunc ℂ) : LaurentSeries ℂ).coeff n

/-- `f` has at `c` a pole of order at most `m` (`m = 0`: no pole at `c`). -/
def poleLE (c : ℂ) (f : RatFunc ℂ) (m : ℕ) : Prop :=
  ∀ n : ℤ, n < -(m : ℤ) → lcoeff c f n = 0

/-- The principal part of `f` at `c` of order `≥ 2`: `Σ_{m=2}^{N} lcoeff c f (-m) (v-c)^{-m}`,
`N` the degree of the denominator of `f` (a bound for every pole order of `f`). -/
noncomputable def ppart2 (c : ℂ) (f : RatFunc ℂ) : RatFunc ℂ :=
  ∑ m ∈ Finset.Icc 2 f.denom.natDegree,
    cst (lcoeff c f (-(m : ℤ))) * (vX - cst c) ^ (-(m : ℤ))

/-- `Hp(W)`: the sum of the principal parts of `W` at `v = i` and `v = -i` of order `≥ 2`. -/
noncomputable def Hp (f : RatFunc ℂ) : RatFunc ℂ :=
  ppart2 Complex.I f + ppart2 (-Complex.I) f

/-- The termwise primitive of `ppart2 c f` that vanishes at `v = b`. -/
noncomputable def jpart (c : ℂ) (f : RatFunc ℂ) : RatFunc ℂ :=
  ∑ m ∈ Finset.Icc 2 f.denom.natDegree,
    cst (lcoeff c f (-(m : ℤ)) * (-(1 / ((m : ℂ) - 1)))) *
      ((vX - cst c) ^ (-((m : ℤ) - 1)) - cst ((bC - c) ^ (-((m : ℤ) - 1))))

/-- `Jp(W)`: the primitive of `Hp(W)` that vanishes at `v = b`. -/
noncomputable def Jp (f : RatFunc ℂ) : RatFunc ℂ :=
  jpart Complex.I f + jpart (-Complex.I) f

/-- `t(v) = (v² + 2v - 1)/(4v) = (v-b)(v+a)/(4v)`. -/
noncomputable def tR : RatFunc ℂ := (vX ^ 2 + 2 * vX - 1) / (4 * vX)

/-- `q = (v² + 1)/(2√2 v)`. -/
noncomputable def qR : RatFunc ℂ := (vX ^ 2 + 1) / (cst (2 * s2) * vX)

/-- `D₁ = ℓ₁' = (v²+1)/(v(v-b)(v+a))`. -/
noncomputable def D1 : RatFunc ℂ := (vX ^ 2 + 1) / (vX * (vX - cst bC) * (vX + cst aC))

/-- `D₂ = ℓ₂' = (v²+1)/(v(v-a)(v+b))`. -/
noncomputable def D2 : RatFunc ℂ := (vX ^ 2 + 1) / (vX * (vX - cst aC) * (vX + cst bC))

/-- `R_wt = 1/(2v) + 1/(2(v-a)) - 1/(v+a) - 1/(2(v+b))`. -/
noncomputable def Rwt : RatFunc ℂ :=
  1 / (2 * vX) + 1 / (2 * (vX - cst aC)) - 1 / (vX + cst aC) - 1 / (2 * (vX + cst bC))

/-- `g₂₀ = a/4 - √2/(2(v²+1))` (g1:lem:k2). -/
noncomputable def g20 : RatFunc ℂ := cst (aC / 4) - cst (s2 / 2) / (vX ^ 2 + 1)

/-- `g₁₁ = -√2(v²+2v-1)/(4(v²+1))` (g1:lem:kernel). -/
noncomputable def g11 : RatFunc ℂ := -(cst s2 * (vX ^ 2 + 2 * vX - 1)) / (4 * (vX ^ 2 + 1))

/-- `τ₁ = (-1-i)/2`. -/
noncomputable def tau1 : ℂ := (-1 - Complex.I) / 2

/-- `τ₂ = (-1+i)/2`. -/
noncomputable def tau2 : ℂ := (-1 + Complex.I) / 2

/-- `binom(1/2, j) = (Π_{m<j} (1/2 - m)) / j!`. -/
noncomputable def binHalf (j : ℕ) : ℂ :=
  (∏ m ∈ Finset.range j, ((1 / 2 : ℂ) - m)) / (j.factorial : ℂ)

/-- `e^{2x} - 1` as a power series. -/
noncomputable def E2 : PowerSeries ℂ := PowerSeries.rescale (2 : ℂ) (PowerSeries.exp ℂ) - 1

/-- `V = (1-t)²(e^{2x}-1) + t²(e^{2y}-1)`, with `t = t(v)`, as a power series in `x` whose
coefficients are power series in `y` over `RatFunc ℂ`. -/
noncomputable def Vser : PowerSeries (PowerSeries (RatFunc ℂ)) :=
  PowerSeries.C (PowerSeries.C ((1 - tR) ^ 2)) *
      PowerSeries.map ((PowerSeries.C).comp (algebraMap ℂ (RatFunc ℂ))) E2 +
    PowerSeries.C (PowerSeries.C (tR ^ 2)) *
      PowerSeries.C (PowerSeries.map (algebraMap ℂ (RatFunc ℂ)) E2)

/-- `U_j^{(r,s)} = [x^r y^s] V^j` (at `t = t(v)`). -/
noncomputable def Ucoef (j r s : ℕ) : RatFunc ℂ :=
  PowerSeries.coeff s (PowerSeries.coeff r (Vser ^ j))

/-- `R̃_{rs}` (degree `d = r+s`) by the formula of g1:lem:init:
`1/(2√2 v t(1-t)) Σ_{j=1}^d binom(1/2,j) q^{2-2j} U_j^{(r,s)} + E_{rs}`, with
`E_{d0} = -D₁/(2 d!)`, `E_{0d} = D₂/(2 d!)`. -/
noncomputable def Rt (r s : ℕ) : RatFunc ℂ :=
  1 / (cst (2 * s2) * vX * tR * (1 - tR)) *
      ∑ j ∈ Finset.Icc 1 (r + s), cst (binHalf j) * qR ^ (2 - 2 * (j : ℤ)) * Ucoef j r s
    + (if s = 0 then -D1 / cst (2 * ((r + s).factorial : ℂ)) else 0)
    + (if r = 0 then D2 / cst (2 * ((r + s).factorial : ℂ)) else 0)

/-- The folded initial data `R^{(k)}_{rs}`: `R̃_{rs} + Fold(R̃_{sr})` for `r > s`, `R̃_{rr}` on the
diagonal, `0` unless `r + s = k` and `r ≥ s`. -/
noncomputable def Rinit (k r s : ℕ) : RatFunc ℂ :=
  if r + s = k ∧ s ≤ r then (if s < r then Rt r s + Fold (Rt s r) else Rt r s) else 0

/-- Definition g1:def:red: `W^{(k)}_{rs}`, by descending induction on `d = r+s`
(`W = 0` for `d > k`). The parent kernels `G^{(k)}_{r+1,s}`, `G^{(k)}_{r,s+1}`, `G^{(k)}_{rr}` have
degree `d+1 ≥ 2` and are `Jp` of the parent `W`. -/
noncomputable def W (k : ℕ) : ℕ → ℕ → RatFunc ℂ
  | r, s =>
    if _h : r + s ≤ k then
      Rinit k r s
        - cst ((r : ℂ) + 1) * Jp (W k (r + 1) s) * D1
        - (if s < r then cst ((s : ℂ) + 1) * Jp (W k r (s + 1)) * D2 else 0)
        - (if _h2 : s + 1 = r then cst (r : ℂ) * Fold (Jp (W k r r) * D1) else 0)
    else 0
termination_by r s => k + 1 - (r + s)
decreasing_by all_goals omega

/-- Definition g1:def:red: `G^{(k)}_{rs} = Jp(W^{(k)}_{rs})` for `d ≥ 2`, `0` for `d ≤ 1`. -/
noncomputable def G (k r s : ℕ) : RatFunc ℂ := if 2 ≤ r + s then Jp (W k r s) else 0

/-- `h^{(k)}_{rs}`: the coefficient of `(v-i)^{-(2d-2)}` in `W^{(k)}_{rs}`, `d = r+s`. -/
noncomputable def h (k r s : ℕ) : ℂ :=
  lcoeff Complex.I (W k r s) (-(2 * ((r + s : ℕ) : ℤ) - 2))

/-- `K_d = (-1)^d (2d-3)!! i^{d+1} / (2√2 d!)`. -/
noncomputable def Kd (d : ℕ) : ℂ :=
  (-1) ^ d * ((Nat.doubleFactorial (2 * d - 3) : ℕ) : ℂ) * Complex.I ^ (d + 1) /
    (2 * s2 * (d.factorial : ℂ))

/-- `ρ_{rs} = K_d (-1)^r (1 + [r>s]) binom(d, r)`, `d = r+s`. -/
noncomputable def rho (r s : ℕ) : ℂ :=
  Kd (r + s) * (-1) ^ r * (if s < r then 2 else 1) * (((r + s).choose r : ℕ) : ℂ)

/-- `ĥ_p = (p+1) τ₁ ρ_{p+1,p} / (4p-1)`. -/
noncomputable def hhat (p : ℕ) : ℂ := ((p : ℂ) + 1) * tau1 * rho (p + 1) p / (4 * (p : ℂ) - 1)

/-- The residue of `f` at `c`. -/
noncomputable def res (c : ℂ) (f : RatFunc ℂ) : ℂ := lcoeff c f (-1)

/-- The simple-pole part of `f` at the real letters `0, a, -a, b, -b`. -/
noncomputable def simplePart (f : RatFunc ℂ) : RatFunc ℂ :=
  cst (res 0 f) / vX + cst (res aC f) / (vX - cst aC) + cst (res (-aC) f) / (vX + cst aC)
    + cst (res bC f) / (vX - cst bC) + cst (res (-bC) f) / (vX + cst bC)

/-- The weight-two part `W₂(k) = ∫_b^a K^{(k)}_{10}(v) ℓ₁(v) dv`. -/
noncomputable def W2 (k : ℕ) : ℂ :=
  ∫ v in Y01.b..Y01.a, ev (simplePart (W k 1 0)) (v : ℂ) * ((Y01.l1 v : ℝ) : ℂ)

/-- The simple kernel `S^{(k)}_{rs} = W - Hp(W)` (`d ≥ 2`), `S^{(k)}_{10} = W^{(k)}_{10}`. -/
noncomputable def S (k r s : ℕ) : RatFunc ℂ :=
  if 2 ≤ r + s then W k r s - Hp (W k r s) else W k r s

/-- `c ∈ ℚ(√2)`. -/
def inQs2 (c : ℂ) : Prop := ∃ p q : ℚ, c = (p : ℂ) + (q : ℂ) * s2

/-- "The operators Hp and Jp": `Hp(W)` collects the terms of order `≥ 2` of `W` at `±i`
(`W - Hp(W)` has at most simple poles there), has poles only at `±i`, no residue, and is
`O(v^{-2})` at infinity; `Jp(W)` is a primitive of `Hp(W)` (off `±i`) that vanishes at `v = b`. -/
def stmt_1 : Prop :=
  ∀ f : RatFunc ℂ,
    poleLE Complex.I (f - Hp f) 1 ∧ poleLE (-Complex.I) (f - Hp f) 1 ∧
    (∀ c : ℂ, c ≠ Complex.I → c ≠ -Complex.I → poleLE c (Hp f) 0) ∧
    res Complex.I (Hp f) = 0 ∧ res (-Complex.I) (Hp f) = 0 ∧
    (Hp f = 0 ∨ (Hp f).intDegree ≤ -2) ∧
    ev (Jp f) bC = 0 ∧
    ∀ z : ℂ, z ≠ Complex.I → z ≠ -Complex.I →
      HasDerivAt (fun w => ev (Jp f) w) (ev (Hp f) z) z

/-- g1:lem:euler "$\Fold(D_1)=-D_2$", and "$D_1$, $D_2$ have simple zeros at $\pm i$ ...
$D_1'(i)=\tau_1$, $D_2'(i)=\tau_2$, $D_1'(-i)=\tau_2$, $D_2'(-i)=\tau_1$". -/
def stmt_2 : Prop :=
  Fold D1 = -D2 ∧
    poleLE Complex.I D1 0 ∧ poleLE Complex.I D2 0 ∧
    poleLE (-Complex.I) D1 0 ∧ poleLE (-Complex.I) D2 0 ∧
    lcoeff Complex.I D1 0 = 0 ∧ lcoeff Complex.I D2 0 = 0 ∧
    lcoeff (-Complex.I) D1 0 = 0 ∧ lcoeff (-Complex.I) D2 0 = 0 ∧
    lcoeff Complex.I D1 1 = tau1 ∧ lcoeff Complex.I D2 1 = tau2 ∧
    lcoeff (-Complex.I) D1 1 = tau2 ∧ lcoeff (-Complex.I) D2 1 = tau1

/-- g1:lem:init, second claim: for `d = r+s = k ≥ 2`, `R̃_{rs}` and (for keys `r ≥ s`)
`R^{(k)}_{rs}` have poles of order at most `2d-2` at `±i`, and the coefficient of
`(v-i)^{-(2d-2)}` in `R^{(k)}_{rs}` is `ρ_{rs}`. -/
def stmt_5 : Prop :=
  ∀ k : ℕ, 2 ≤ k → ∀ r s : ℕ, r + s = k →
    (poleLE Complex.I (Rt r s) (2 * k - 2) ∧ poleLE (-Complex.I) (Rt r s) (2 * k - 2)) ∧
    (s ≤ r →
      poleLE Complex.I (Rinit k r s) (2 * k - 2) ∧ poleLE (-Complex.I) (Rinit k r s) (2 * k - 2) ∧
      lcoeff Complex.I (Rinit k r s) (-(2 * (k : ℤ) - 2)) = rho r s)

/-- g1:lem:A: for every `k ≥ 1` and every key `(r,s)` of degree `d ≥ 1`, `W^{(k)}_{rs}` has poles at
`±i` of order at most `2d-2`, and for `d ≥ 2`, `G^{(k)}_{rs}` has poles at `±i` of order at most
`2d-3` and no other poles. -/
def stmt_6 : Prop :=
  ∀ k : ℕ, 1 ≤ k → ∀ r s : ℕ, s ≤ r → 1 ≤ r + s →
    poleLE Complex.I (W k r s) (2 * (r + s) - 2) ∧
    poleLE (-Complex.I) (W k r s) (2 * (r + s) - 2) ∧
    (2 ≤ r + s →
      poleLE Complex.I (G k r s) (2 * (r + s) - 3) ∧
      poleLE (-Complex.I) (G k r s) (2 * (r + s) - 3) ∧
      ∀ c : ℂ, c ≠ Complex.I → c ≠ -Complex.I → poleLE c (G k r s) 0)

/-- g1:lem:B, the corrected recursion: for `k ≥ d = r+s ≥ 2`, `r ≥ s`,
`(2d-1)(h_{rs} - δ_{kd} ρ_{rs}) = (r+1)τ₁h_{r+1,s} + [r ≥ s+2](s+1)τ₂h_{r,s+1}
 + [s = r-1] r τ₂ (h_{rr} - conj h_{rr})`. -/
def stmt_7 : Prop :=
  ∀ k r s : ℕ, s ≤ r → 2 ≤ r + s → r + s ≤ k →
    (2 * ((r + s : ℕ) : ℂ) - 1) * (h k r s - (if k = r + s then rho r s else 0))
      = ((r : ℂ) + 1) * tau1 * h k (r + 1) s
        + (if s + 2 ≤ r then ((s : ℂ) + 1) * tau2 * h k r (s + 1) else 0)
        + (if s + 1 = r then (r : ℂ) * tau2 * (h k r r - (starRingEnd ℂ) (h k r r)) else 0)

/-- g1:lem:I, the four identities (I1)-(I4). -/
def stmt_8 : Prop :=
  (∀ d : ℕ, 2 ≤ d → ((d : ℂ) + 1) * Kd (d + 1) = -Complex.I * (2 * (d : ℂ) - 1) * Kd d) ∧
  (∀ r s : ℕ, 2 ≤ r + s → s + 2 ≤ r →
    (2 * ((r + s : ℕ) : ℂ) - 1) * rho r s
      = ((r : ℂ) + 1) * tau1 * rho (r + 1) s + ((s : ℂ) + 1) * tau2 * rho r (s + 1)) ∧
  (∀ p : ℕ, 2 ≤ p →
    (4 * (p : ℂ) - 3) * rho p (p - 1)
      = ((p : ℂ) + 1) * tau1 * rho (p + 1) (p - 1)
        + (p : ℂ) * tau2 * (rho p p - (starRingEnd ℂ) (rho p p))) ∧
  (∀ p : ℕ, 1 ≤ p →
    hhat p = (1 - Complex.I) * rho p p ∧ (rho p p).re = 0 ∧
    hhat p - (starRingEnd ℂ) (hhat p) = rho p p - (starRingEnd ℂ) (rho p p))

/-- g1:lem:D: for every `k ≥ r+s = d ≥ 2`: if `r > s` then `h^{(k)}_{rs} = ρ_{rs}`; on the diagonal
`h^{(2p)}_{pp} = ρ_{pp}` and `h^{(k)}_{pp} = ĥ_p` for `k ≥ 2p+1`. -/
def stmt_9 : Prop :=
  ∀ k r s : ℕ, s ≤ r → 2 ≤ r + s → r + s ≤ k →
    (s < r → h k r s = rho r s) ∧
    (s = r → (k = 2 * r → h k r r = rho r r) ∧ (2 * r + 1 ≤ k → h k r r = hhat r))

/-- g1:lem:kernel, "More precisely": `G^{(k)}_{20} = g₂₀` for all `k ≥ 2`, `G^{(2)}_{11} = -g₂₀` and
`G^{(k)}_{11} = g₁₁` for `k ≥ 3`. -/
def stmt_10 : Prop :=
  (∀ k : ℕ, 2 ≤ k → G k 2 0 = g20) ∧ G 2 1 1 = -g20 ∧ (∀ k : ℕ, 3 ≤ k → G k 1 1 = g11)

/-- g1:lem:kernel: "For every $k\ge2$: $W^{(k)}_{10}=\Rwt$." -/
def stmt_11 : Prop := ∀ k : ℕ, 2 ≤ k → W k 1 0 = Rwt

/-- Theorem g1:thm:wt2 (thm:wt2) in the supplement's precise form: for every `k ≥ 2` the
degree-one kernel produced by the reduction is `W^{(k)}_{10} = R_wt`, it is its own simple-pole
part `K^{(k)}_{10}` at the real letters, and the weight-two part
`W₂(k) = ∫_b^a K^{(k)}_{10} ℓ₁ dv` (an integrable integrand) is exactly `T = (π²+6λ²)/24`. -/
def stmt_12 : Prop :=
  ∀ k : ℕ, 2 ≤ k →
    W k 1 0 = Rwt ∧ simplePart (W k 1 0) = W k 1 0 ∧
    IntervalIntegrable
      (fun v : ℝ => ev (simplePart (W k 1 0)) (v : ℂ) * ((Y01.l1 v : ℝ) : ℂ)) volume Y01.b Y01.a ∧
    W2 k = (((π ^ 2 + 6 * Y01.lam ^ 2) / 24 : ℝ) : ℂ)

/-- g1:lem:inv (a), (b) first part, (c) the two formulas, (d). -/
def stmt_13 : Prop :=
  (∀ f g : RatFunc ℂ, Inv (f + g) = Inv f + Inv g) ∧
  (∀ (c : ℂ) (f : RatFunc ℂ), Inv (cst c * f) = cst c * Inv f) ∧
  (∀ f : RatFunc ℂ, Inv (Inv f) = f) ∧
  (∀ f : RatFunc ℂ, Inv (Fold f) = Fold (Inv f)) ∧
  Inv D1 = D1 ∧ Inv D2 = D2 ∧
  (∀ m : ℕ, 1 ≤ m →
    Inv ((vX - cst Complex.I) ^ (-(m : ℤ)))
        = cst (Complex.I ^ m) * vX ^ ((m : ℤ) - 2) / (vX - cst Complex.I) ^ m ∧
    Inv ((vX + cst Complex.I) ^ (-(m : ℤ)))
        = cst ((-Complex.I) ^ m) * vX ^ ((m : ℤ) - 2) / (vX + cst Complex.I) ^ m) ∧
  (∀ f : RatFunc ℂ,
    Hp (Inv f) = Inv (Hp f) ∧ res Complex.I (Inv f) = res Complex.I f ∧
      res (-Complex.I) (Inv f) = res (-Complex.I) f)

/-- g1:lem:par: for `k ≥ 2` and every key, `W + 𝓘(W)` has no pole at `±i`; for every key of
degree `d ≥ 2`, `G(v) + G(-1/v)` is constant. -/
def stmt_14 : Prop :=
  ∀ k : ℕ, 2 ≤ k → ∀ r s : ℕ, s ≤ r → 1 ≤ r + s →
    poleLE Complex.I (W k r s + Inv (W k r s)) 0 ∧
    poleLE (-Complex.I) (W k r s + Inv (W k r s)) 0 ∧
    (2 ≤ r + s → ∃ c : ℂ, G k r s + subst (-vX⁻¹) (G k r s) = cst c)

/-- g1:lem:res: for `k ≥ 2` and every key, `Res_{±i} W^{(k)}_{rs} = 0`, and the simple kernel
`S^{(k)}_{rs}` has no pole at `±i`. -/
def stmt_15 : Prop :=
  ∀ k : ℕ, 2 ≤ k → ∀ r s : ℕ, s ≤ r → 1 ≤ r + s →
    res Complex.I (W k r s) = 0 ∧ res (-Complex.I) (W k r s) = 0 ∧
    poleLE Complex.I (S k r s) 0 ∧ poleLE (-Complex.I) (S k r s) 0

end P01
-- END STATEMENTS

namespace P01

lemma s2_ne : s2 ≠ 0 := by
  unfold s2; exact_mod_cast (Real.sqrt_pos.2 (by norm_num : (0:ℝ) < 2)).ne'

lemma Kd_succ (d : ℕ) (hd : 2 ≤ d) :
    ((d : ℂ) + 1) * Kd (d + 1) = -Complex.I * (2 * (d : ℂ) - 1) * Kd d := by
  obtain ⟨n, rfl⟩ : ∃ n, d = n + 2 := ⟨d - 2, by omega⟩
  unfold Kd
  have e1 : 2 * (n + 2 + 1) - 3 = (2 * n + 1) + 2 := by omega
  have e2 : 2 * (n + 2) - 3 = 2 * n + 1 := by omega
  rw [e1, e2, Nat.doubleFactorial_add_two, Nat.factorial_succ (n + 2)]
  have := s2_ne
  have hf : ((n + 2).factorial : ℂ) ≠ 0 := by exact_mod_cast (Nat.factorial_pos _).ne'
  have h3 : (n : ℂ) + 2 + 1 ≠ 0 := by norm_cast
  push_cast
  field_simp
  ring_nf

lemma rho_pp (p : ℕ) : ∃ x : ℝ, rho p p = (x : ℂ) * Complex.I := by
  unfold rho Kd
  refine ⟨(-1) ^ (p + p) * ((Nat.doubleFactorial (2 * (p + p) - 3) : ℕ) : ℝ) * (-1) ^ p /
      (2 * Real.sqrt 2 * ((p + p).factorial : ℝ)) * (-1) ^ p * 1 * (((p + p).choose p : ℕ) : ℝ), ?_⟩
  have h : Complex.I ^ (p + p + 1) = (-1) ^ p * Complex.I := by
    rw [pow_succ, ← two_mul, pow_mul, Complex.I_sq]
  simp only [lt_irrefl, ite_false, h, s2]
  push_cast
  ring

lemma conj_rho_pp (p : ℕ) : (starRingEnd ℂ) (rho p p) = - rho p p := by
  obtain ⟨x, hx⟩ := rho_pp p
  rw [hx, map_mul, Complex.conj_ofReal, Complex.conj_I]; ring

end P01


namespace P01

open Polynomial in
lemma X_sub_C_ne (c : ℂ) : (RatFunc.X : RatFunc ℂ) - RatFunc.C c ≠ 0 := by
  rw [← RatFunc.algebraMap_X, ← RatFunc.algebraMap_C, ← map_sub]
  exact (map_ne_zero_iff _ (RatFunc.algebraMap_injective ℂ)).2 (X_sub_C_ne_zero c)

lemma X_ne_C (c : ℂ) : (RatFunc.X : RatFunc ℂ) ≠ RatFunc.C c := fun h => X_sub_C_ne c (by rw [h, sub_self])

lemma X_ne_zero' : (RatFunc.X : RatFunc ℂ) ≠ 0 := by
  have := X_ne_C 0; rwa [map_zero] at this

lemma Xinv_sub_C_ne (c : ℂ) : (RatFunc.X : RatFunc ℂ)⁻¹ - RatFunc.C c ≠ 0 := by
  intro h
  rw [sub_eq_zero] at h
  have h2 : (RatFunc.X : RatFunc ℂ) = RatFunc.C c⁻¹ := by
    rw [map_inv₀, ← h, inv_inv]
  exact X_ne_C _ h2

lemma not_const_of (g : RatFunc ℂ) (hg : ∀ c, g - RatFunc.C c ≠ 0) : ¬ ∃ c, g = RatFunc.C c := by
  rintro ⟨c, hc⟩; exact hg c (by rw [hc, sub_self])

lemma inj_of (g : RatFunc ℂ) (hg : ∀ c, g - RatFunc.C c ≠ 0) :
    Function.Injective (Polynomial.aeval (R := ℂ) g) :=
  (transcendental_iff_injective).1 (RatFunc.transcendental_of_ne_C _ (not_const_of g hg))

/-- substitution as an algebra hom. -/
noncomputable def sig (g : RatFunc ℂ) (hg : ∀ c, g - RatFunc.C c ≠ 0) : RatFunc ℂ →ₐ[ℂ] RatFunc ℂ :=
  RatFunc.liftAlgHom (Polynomial.aeval g)
    (nonZeroDivisors_le_comap_nonZeroDivisors_of_injective _ (inj_of g hg))

lemma subst_eq (g : RatFunc ℂ) (hg : ∀ c, g - RatFunc.C c ≠ 0) (f : RatFunc ℂ) :
    subst g f = sig g hg f := by
  unfold subst sig; rw [RatFunc.liftAlgHom_apply]

lemma sig_X (g : RatFunc ℂ) (hg : ∀ c, g - RatFunc.C c ≠ 0) : sig g hg RatFunc.X = g := by
  have := RatFunc.liftAlgHom_apply_div (Polynomial.aeval g)
    (nonZeroDivisors_le_comap_nonZeroDivisors_of_injective _ (inj_of g hg)) Polynomial.X 1
  rw [map_one, div_one, RatFunc.algebraMap_X] at this
  unfold sig; rw [this]; simp

lemma sig_C (g : RatFunc ℂ) (hg : ∀ c, g - RatFunc.C c ≠ 0) (c : ℂ) : sig g hg (RatFunc.C c) = RatFunc.C c := by
  rw [← RatFunc.algebraMap_eq_C, AlgHom.commutes]

lemma hXinv : ∀ c, (RatFunc.X : RatFunc ℂ)⁻¹ - RatFunc.C c ≠ 0 := Xinv_sub_C_ne

lemma hnXinv : ∀ c, -(RatFunc.X : RatFunc ℂ)⁻¹ - RatFunc.C c ≠ 0 := by
  intro c h
  have : (RatFunc.X : RatFunc ℂ)⁻¹ - RatFunc.C (-c) = 0 := by
    rw [map_neg]; linear_combination -h
  exact Xinv_sub_C_ne _ this

lemma s2_sq : s2 ^ 2 = 2 := by
  unfold s2; rw [← Complex.ofReal_pow, Real.sq_sqrt (by norm_num)]; norm_num

lemma aC_eq : aC = 1 + s2 := by unfold aC Y01.a s2; push_cast; ring
lemma bC_eq : bC = s2 - 1 := by unfold bC Y01.b s2; push_cast; ring

lemma fold_D1 : Fold D1 = -D2 := by
  unfold Fold D1 D2
  simp only [vX, cst]
  rw [subst_eq _ hXinv]
  simp only [vX, cst, map_div₀, map_mul, map_add, map_sub, map_pow, map_one, sig_X, sig_C]
  have hX := X_ne_zero'
  have h1 := Xinv_sub_C_ne bC
  have h2 : (RatFunc.X : RatFunc ℂ)⁻¹ + RatFunc.C aC ≠ 0 := by
    have := Xinv_sub_C_ne (-aC); rwa [map_neg, sub_neg_eq_add] at this
  have h3 := X_sub_C_ne aC
  have h4 : (RatFunc.X : RatFunc ℂ) + RatFunc.C bC ≠ 0 := by
    have := X_sub_C_ne (-bC); rwa [map_neg, sub_neg_eq_add] at this
  have ha : RatFunc.C aC = RatFunc.C bC + 2 := by
    rw [aC_eq, bC_eq, show (1 : ℂ) + s2 = (s2 - 1) + 2 by ring, map_add, map_ofNat]
  have hb : RatFunc.C bC * RatFunc.C bC + 2 * RatFunc.C bC = (1 : RatFunc ℂ) := by
    rw [← map_mul, show (2 : RatFunc ℂ) = RatFunc.C 2 from (map_ofNat _ 2).symm, ← map_mul, ← map_add, bC_eq]
    rw [show (s2 - 1) * (s2 - 1) + 2 * (s2 - 1) = s2 ^ 2 - 1 by ring, s2_sq]; norm_num
  have e1 : (RatFunc.X : RatFunc ℂ)⁻¹ ^ 2 + 1 = (1 + RatFunc.X ^ 2) / RatFunc.X ^ 2 := by
    field_simp
  have e2 : (RatFunc.X : RatFunc ℂ)⁻¹ - RatFunc.C bC = (1 - RatFunc.C bC * RatFunc.X) / RatFunc.X := by
    field_simp
  have e3 : (RatFunc.X : RatFunc ℂ)⁻¹ + RatFunc.C aC = (1 + RatFunc.C aC * RatFunc.X) / RatFunc.X := by
    field_simp
  have h1' : 1 - RatFunc.C bC * RatFunc.X ≠ 0 := by
    intro h; apply h1; rw [e2, h, zero_div]
  have h2' : 1 + RatFunc.C aC * RatFunc.X ≠ 0 := by
    intro h; apply h2; rw [e3, h, zero_div]
  rw [e1, e2, e3, neg_div', div_div, div_eq_div_iff (by
    refine mul_ne_zero (mul_ne_zero (mul_ne_zero (inv_ne_zero hX) (div_ne_zero h1' hX))
      (div_ne_zero h2' hX)) (pow_ne_zero _ hX)) (mul_ne_zero (mul_ne_zero hX h3) h4)]
  field_simp
  rw [ha]
  linear_combination (-(1 + (RatFunc.X : RatFunc ℂ) ^ 2) ^ 2) * hb

open Polynomial in
lemma lcoeff_reg (c : ℂ) (N Q : ℂ[X]) (hQ : Q.eval c ≠ 0) (n : ℤ) :
    lcoeff c (algebraMap _ (RatFunc ℂ) N / algebraMap _ _ Q) n
      = if n < 0 then 0 else PowerSeries.coeff n.natAbs
          (((taylor c N : ℂ[X]) : PowerSeries ℂ) * ((taylor c Q : ℂ[X]) : PowerSeries ℂ)⁻¹) := by
  unfold lcoeff
  rw [RatFunc.laurent_div]
  set P := taylor c N
  set R := taylor c Q
  have hR0 : PowerSeries.constantCoeff (R : PowerSeries ℂ) ≠ 0 := by
    rw [← PowerSeries.coeff_zero_eq_constantCoeff_apply, Polynomial.coeff_coe, taylor_coeff_zero]
    exact hQ
  have key : ((algebraMap _ (RatFunc ℂ) P / algebraMap _ _ R : RatFunc ℂ) : LaurentSeries ℂ)
      = (((P : PowerSeries ℂ) * (R : PowerSeries ℂ)⁻¹ : PowerSeries ℂ) : LaurentSeries ℂ) := by
    rw [map_div₀]
    have hP : ((algebraMap ℂ[X] (RatFunc ℂ) P : RatFunc ℂ) : LaurentSeries ℂ)
        = ((P : PowerSeries ℂ) : LaurentSeries ℂ) := (RatFunc.coe_coe P).symm
    have hR : ((algebraMap ℂ[X] (RatFunc ℂ) R : RatFunc ℂ) : LaurentSeries ℂ)
        = ((R : PowerSeries ℂ) : LaurentSeries ℂ) := (RatFunc.coe_coe R).symm
    rw [hP, hR, PowerSeries.coe_mul, div_eq_mul_inv]
    congr 1
    symm
    apply eq_inv_of_mul_eq_one_right
    rw [← PowerSeries.coe_mul, PowerSeries.mul_inv_cancel _ hR0, PowerSeries.coe_one]
  rw [key, PowerSeries.coeff_coe]

open Polynomial in
lemma ps_coeff0 (c : ℂ) (N Q : ℂ[X]) :
    PowerSeries.coeff 0 (((taylor c N : ℂ[X]) : PowerSeries ℂ) * ((taylor c Q : ℂ[X]) : PowerSeries ℂ)⁻¹)
      = N.eval c * (Q.eval c)⁻¹ := by
  rw [PowerSeries.coeff_zero_eq_constantCoeff_apply, map_mul, PowerSeries.constantCoeff_inv,
    ← PowerSeries.coeff_zero_eq_constantCoeff_apply, ← PowerSeries.coeff_zero_eq_constantCoeff_apply,
    Polynomial.coeff_coe, Polynomial.coeff_coe, taylor_coeff_zero, taylor_coeff_zero]

open Polynomial in
lemma ps_coeff1 (c : ℂ) (N Q : ℂ[X]) (hN : N.eval c = 0) :
    PowerSeries.coeff 1 (((taylor c N : ℂ[X]) : PowerSeries ℂ) * ((taylor c Q : ℂ[X]) : PowerSeries ℂ)⁻¹)
      = N.derivative.eval c * (Q.eval c)⁻¹ := by
  rw [PowerSeries.coeff_mul, Finset.Nat.antidiagonal_succ]
  simp only [Finset.Nat.antidiagonal_zero, Finset.sum_cons, Finset.sum_singleton, Prod.fst,
    Prod.snd, Polynomial.coeff_coe, taylor_coeff_zero, hN, zero_mul, zero_add, taylor_coeff_one,
    Finset.map_singleton]
  simp [PowerSeries.coeff_inv, Polynomial.coeff_coe, taylor_coeff_zero]

open Polynomial in
lemma D_alg (p q : ℂ) :
    (vX ^ 2 + 1) / (vX * (vX - cst p) * (vX + cst q))
      = algebraMap ℂ[X] (RatFunc ℂ) (X ^ 2 + 1) / algebraMap ℂ[X] (RatFunc ℂ) (X * (X - C p) * (X + C q)) := by
  simp only [map_add, map_mul, map_pow, map_sub, map_one, RatFunc.algebraMap_X, RatFunc.algebraMap_C, vX, cst]

open Polynomial in
lemma D_facts (p q c : ℂ) (hc : c ^ 2 + 1 = 0) (hQ : c * (c - p) * (c + q) ≠ 0) :
    poleLE c ((vX ^ 2 + 1) / (vX * (vX - cst p) * (vX + cst q))) 0 ∧
    lcoeff c ((vX ^ 2 + 1) / (vX * (vX - cst p) * (vX + cst q))) 0 = 0 ∧
    lcoeff c ((vX ^ 2 + 1) / (vX * (vX - cst p) * (vX + cst q))) 1 = 2 * c * (c * (c - p) * (c + q))⁻¹ := by
  have hQ' : (X * (X - C p) * (X + C q) : ℂ[X]).eval c ≠ 0 := by simpa using hQ
  have hN : (X ^ 2 + 1 : ℂ[X]).eval c = 0 := by simpa using hc
  rw [D_alg]
  refine ⟨?_, ?_, ?_⟩
  · intro n hn
    rw [lcoeff_reg c _ _ hQ']
    simp only [Nat.cast_zero, neg_zero] at hn
    simp [hn]
  · rw [lcoeff_reg c _ _ hQ']
    simp only [lt_irrefl, ite_false]
    rw [show (0 : ℤ).natAbs = 0 from rfl, ps_coeff0, hN, zero_mul]
  · rw [lcoeff_reg c _ _ hQ']
    simp only [show ¬ ((1 : ℤ) < 0) by norm_num, ite_false]
    rw [show (1 : ℤ).natAbs = 1 from rfl, ps_coeff1 c _ _ hN]
    simp <;> norm_num

lemma hI3 : Complex.I ^ 3 = -Complex.I := by rw [pow_succ, Complex.I_sq]; ring

lemma num1 : Complex.I * (Complex.I - bC) * (Complex.I + aC) = -2 - 2 * Complex.I := by
  rw [aC_eq, bC_eq]; ring_nf; simp only [s2_sq, Complex.I_sq, hI3]; ring
lemma num2 : Complex.I * (Complex.I - aC) * (Complex.I + bC) = 2 - 2 * Complex.I := by
  rw [aC_eq, bC_eq]; ring_nf; simp only [s2_sq, Complex.I_sq, hI3]; ring
lemma num3 : -Complex.I * (-Complex.I - bC) * (-Complex.I + aC) = -2 + 2 * Complex.I := by
  rw [aC_eq, bC_eq]; ring_nf; simp only [s2_sq, Complex.I_sq, hI3]; ring
lemma num4 : -Complex.I * (-Complex.I - aC) * (-Complex.I + bC) = 2 + 2 * Complex.I := by
  rw [aC_eq, bC_eq]; ring_nf; simp only [s2_sq, Complex.I_sq, hI3]; ring

end P01




namespace P01

lemma lcoeff_add (c : ℂ) (f g : RatFunc ℂ) (n : ℤ) :
    lcoeff c (f + g) n = lcoeff c f n + lcoeff c g n := by
  unfold lcoeff; rw [map_add, map_add, HahnSeries.coeff_add]

lemma lcoeff_neg (c : ℂ) (f : RatFunc ℂ) (n : ℤ) :
    lcoeff c (-f) n = -lcoeff c f n := by
  unfold lcoeff; rw [map_neg, map_neg, HahnSeries.coeff_neg]

lemma lcoeff_sub (c : ℂ) (f g : RatFunc ℂ) (n : ℤ) :
    lcoeff c (f - g) n = lcoeff c f n - lcoeff c g n := by
  rw [sub_eq_add_neg, lcoeff_add, lcoeff_neg, ← sub_eq_add_neg]

lemma lcoeff_zero (c : ℂ) (n : ℤ) : lcoeff c 0 n = 0 := by
  unfold lcoeff; simp

lemma lcoeff_sum {ι : Type*} (s : Finset ι) (c : ℂ) (F : ι → RatFunc ℂ) (n : ℤ) :
    lcoeff c (∑ i ∈ s, F i) n = ∑ i ∈ s, lcoeff c (F i) n := by
  classical
  induction s using Finset.induction_on with
  | empty => simp [lcoeff_zero]
  | insert a s ha ih => rw [Finset.sum_insert ha, Finset.sum_insert ha, lcoeff_add, ih]

lemma lcoeff_cst_mul (c a : ℂ) (f : RatFunc ℂ) (n : ℤ) :
    lcoeff c (cst a * f) n = a * lcoeff c f n := by
  unfold lcoeff cst
  rw [map_mul, RatFunc.laurent_C, map_mul]
  have : ((RatFunc.C a : RatFunc ℂ) : LaurentSeries ℂ) = HahnSeries.C a := by
    rw [← RatFunc.algebraMap_C]
    exact (RatFunc.coe_coe _).symm.trans (by rw [Polynomial.coe_C, PowerSeries.coe_C])
  rw [this, HahnSeries.C_mul_eq_smul, HahnSeries.coeff_smul, smul_eq_mul]

lemma lcoeff_one (c : ℂ) (n : ℤ) : lcoeff c 1 n = if n = 0 then 1 else 0 := by
  unfold lcoeff; rw [map_one, map_one, HahnSeries.coeff_one]; split_ifs <;> simp_all

lemma lcoeff_shift (c : ℂ) (k : ℕ) (g : RatFunc ℂ) (n : ℤ) :
    lcoeff c ((vX - cst c) ^ (-(k : ℤ)) * g) n = lcoeff c g (n + k) := by
  unfold lcoeff vX cst
  rw [map_mul, map_zpow₀, map_sub, RatFunc.laurent_X, RatFunc.laurent_C, add_sub_cancel_right,
    map_mul, map_zpow₀, RatFunc.coe_X, ← RatFunc.single_zpow]
  have : n = (n + k) + (-(k : ℤ)) := by ring
  conv_lhs => rw [this]
  rw [HahnSeries.coeff_single_mul_add, one_mul]

lemma lcoeff_pole (c : ℂ) (m : ℕ) (n : ℤ) :
    lcoeff c ((vX - cst c) ^ (-(m : ℤ))) n = if n = -(m : ℤ) then 1 else 0 := by
  have := lcoeff_shift c m 1 n
  rw [mul_one] at this
  rw [this, lcoeff_one]
  congr 1
  apply propext; constructor <;> intro h <;> omega

open Polynomial in
lemma pole_alg (d : ℂ) (m : ℕ) :
    (vX - cst d) ^ (-(m : ℤ)) = algebraMap ℂ[X] (RatFunc ℂ) 1 / algebraMap ℂ[X] (RatFunc ℂ) ((X - C d) ^ m) := by
  simp only [map_one, map_pow, map_sub, RatFunc.algebraMap_X, RatFunc.algebraMap_C, vX, cst,
    zpow_neg, zpow_natCast, one_div]

lemma lcoeff_other (c d : ℂ) (hcd : c ≠ d) (m : ℕ) (n : ℤ) (hn : n < 0) :
    lcoeff c ((vX - cst d) ^ (-(m : ℤ))) n = 0 := by
  rw [pole_alg, lcoeff_reg]
  · simp [hn]
  · simp [sub_ne_zero.2 hcd]

open Polynomial in
lemma pole_bound (c : ℂ) (f : RatFunc ℂ) (n : ℤ) (hn : n < -(f.denom.natDegree : ℤ)) :
    lcoeff c f n = 0 := by
  obtain ⟨q, hq, hnd⟩ := exists_eq_pow_rootMultiplicity_mul_and_not_dvd f.denom f.denom_ne_zero c
  set k := f.denom.rootMultiplicity c
  have hq0 : q ≠ 0 := by
    rintro rfl; rw [mul_zero] at hq; exact f.denom_ne_zero hq
  have hk : k ≤ f.denom.natDegree := by
    have := congrArg natDegree hq
    rw [natDegree_mul (pow_ne_zero _ (X_sub_C_ne_zero c)) hq0, natDegree_pow, natDegree_X_sub_C,
      mul_one] at this
    omega
  have hqc : q.eval c ≠ 0 := by
    intro h; exact hnd (dvd_iff_isRoot.2 h)
  have hf : f = (vX - cst c) ^ (-(k : ℤ)) * (algebraMap ℂ[X] (RatFunc ℂ) f.num / algebraMap _ _ q) := by
    conv_lhs => rw [← RatFunc.num_div_denom f, hq]
    rw [pole_alg]
    simp only [map_mul, map_one]
    have h1 : algebraMap ℂ[X] (RatFunc ℂ) ((X - C c) ^ k) ≠ 0 := by
      rw [map_ne_zero_iff _ (RatFunc.algebraMap_injective ℂ)]; exact pow_ne_zero _ (X_sub_C_ne_zero c)
    have h2 : algebraMap ℂ[X] (RatFunc ℂ) q ≠ 0 := by
      rw [map_ne_zero_iff _ (RatFunc.algebraMap_injective ℂ)]; exact hq0
    field_simp
  rw [hf, lcoeff_shift, lcoeff_reg c _ _ hqc]
  have : n + k < 0 := by omega
  simp [this]

end P01

namespace P01

lemma lcoeff_ppart_self (c : ℂ) (f : RatFunc ℂ) (j : ℕ) :
    lcoeff c (ppart2 c f) (-(j : ℤ)) =
      if j ∈ Finset.Icc 2 f.denom.natDegree then lcoeff c f (-(j : ℤ)) else 0 := by
  unfold ppart2
  rw [lcoeff_sum]
  simp only [lcoeff_cst_mul, lcoeff_pole, neg_inj, Nat.cast_inj, mul_ite, mul_one, mul_zero]
  rw [Finset.sum_ite_eq]

lemma lcoeff_ppart_other (c d : ℂ) (hcd : c ≠ d) (f : RatFunc ℂ) (n : ℤ) (hn : n < 0) :
    lcoeff c (ppart2 d f) n = 0 := by
  unfold ppart2
  rw [lcoeff_sum]
  simp only [lcoeff_cst_mul, lcoeff_other c d hcd _ n hn, mul_zero, Finset.sum_const_zero]

lemma I_ne_negI : Complex.I ≠ -Complex.I := by
  intro h; have := congrArg Complex.im h; simp at this; norm_num at this

lemma Hp_neg_coeff (f : RatFunc ℂ) (c : ℂ) (hc : c = Complex.I ∨ c = -Complex.I) (j : ℕ) :
    lcoeff c (Hp f) (-(j : ℤ)) = (if j ∈ Finset.Icc 2 f.denom.natDegree then lcoeff c f (-(j : ℤ)) else 0) ∨ j = 0 := by
  rcases Nat.eq_zero_or_pos j with h | h
  · right; exact h
  left
  have hn : (-(j : ℤ)) < 0 := by omega
  unfold Hp
  rw [lcoeff_add]
  rcases hc with rfl | rfl
  · rw [lcoeff_ppart_self, lcoeff_ppart_other _ _ I_ne_negI _ _ hn, add_zero]
  · rw [lcoeff_ppart_self, lcoeff_ppart_other _ _ (Ne.symm I_ne_negI) _ _ hn, zero_add]

lemma sub_Hp_pole (f : RatFunc ℂ) (c : ℂ) (hc : c = Complex.I ∨ c = -Complex.I) :
    poleLE c (f - Hp f) 1 := by
  intro n hn
  obtain ⟨j, rfl⟩ : ∃ j : ℕ, n = -(j : ℤ) := ⟨n.natAbs, by omega⟩
  rw [lcoeff_sub]
  rcases Hp_neg_coeff f c hc j with h | h
  · rw [h]
    split_ifs with hj
    · ring
    · rw [Finset.mem_Icc] at hj
      rw [pole_bound c f _ (by omega), sub_zero]
  · omega

lemma Hp_regular (f : RatFunc ℂ) (c : ℂ) (h1 : c ≠ Complex.I) (h2 : c ≠ -Complex.I) :
    poleLE c (Hp f) 0 := by
  intro n hn
  have hn' : n < 0 := by simpa using hn
  unfold Hp
  rw [lcoeff_add, lcoeff_ppart_other _ _ h1 _ _ hn', lcoeff_ppart_other _ _ h2 _ _ hn', add_zero]

lemma Hp_res (f : RatFunc ℂ) (c : ℂ) (hc : c = Complex.I ∨ c = -Complex.I) :
    res c (Hp f) = 0 := by
  unfold res
  have := Hp_neg_coeff f c hc 1
  simp only [Nat.cast_one, Finset.mem_Icc, one_ne_zero, or_false] at this
  rw [this]
  simp

end P01

namespace P01

def Small (g : RatFunc ℂ) : Prop := g = 0 ∨ g.intDegree ≤ -2

lemma small_add (x y : RatFunc ℂ) (hx : Small x) (hy : Small y) : Small (x + y) := by
  by_cases hxy : x + y = 0
  · left; exact hxy
  right
  rcases hx with rfl | hx
  · rcases hy with rfl | hy
    · simp at hxy
    · simpa using hy
  rcases hy with rfl | hy
  · simpa using hx
  exact (RatFunc.intDegree_add_le hxy).trans (max_le hx hy)

lemma small_sum {ι : Type*} (s : Finset ι) (F : ι → RatFunc ℂ) (h : ∀ i ∈ s, Small (F i)) :
    Small (∑ i ∈ s, F i) := by
  classical
  induction s using Finset.induction_on with
  | empty => left; simp
  | insert a s ha ih =>
    rw [Finset.sum_insert ha]
    exact small_add _ _ (h a (Finset.mem_insert_self a s))
      (ih (fun i hi => h i (Finset.mem_insert_of_mem hi)))

open Polynomial in
lemma small_term (a d : ℂ) (m : ℕ) (hm : 2 ≤ m) : Small (cst a * (vX - cst d) ^ (-(m : ℤ))) := by
  by_cases ha : a = 0
  · left; simp [ha, cst]
  right
  have hc : cst a ≠ 0 := by unfold cst; simpa using ha
  have hp : (vX - cst d) ^ (-(m : ℤ)) ≠ 0 := zpow_ne_zero _ (X_sub_C_ne d)
  rw [RatFunc.intDegree_mul hc hp]
  have h1 : (cst a).intDegree = 0 := RatFunc.intDegree_C a
  have h2 : ((vX - cst d) ^ (-(m : ℤ))).intDegree = -(m : ℤ) := by
    rw [pole_alg, map_one, one_div, RatFunc.intDegree_inv, RatFunc.intDegree_polynomial,
      natDegree_pow, natDegree_X_sub_C, mul_one]
  rw [h1, h2]; omega

lemma Hp_small (f : RatFunc ℂ) : Hp f = 0 ∨ (Hp f).intDegree ≤ -2 := by
  have h : ∀ c, Small (ppart2 c f) := fun c =>
    small_sum _ _ (fun m hm => small_term _ _ m (Finset.mem_Icc.1 hm).1)
  exact small_add _ _ (h _) (h _)

end P01

namespace P01

open Polynomial

/-- `f` is `p/q` with `q(z) ≠ 0` and value `v = p(z)/q(z)` at `z`. -/
def NV (z : ℂ) (f : RatFunc ℂ) (v : ℂ) : Prop :=
  ∃ p q : ℂ[X], q.eval z ≠ 0 ∧ f = algebraMap ℂ[X] (RatFunc ℂ) p / algebraMap _ _ q ∧
    v = p.eval z / q.eval z

lemma alg_ne (q : ℂ[X]) (hq : q ≠ 0) : algebraMap ℂ[X] (RatFunc ℂ) q ≠ 0 :=
  (map_ne_zero_iff _ (RatFunc.algebraMap_injective ℂ)).2 hq

lemma ne_of_eval (q : ℂ[X]) (z : ℂ) (h : q.eval z ≠ 0) : q ≠ 0 := by
  rintro rfl; simp at h

lemma ev_of_NV {z : ℂ} {f : RatFunc ℂ} {v : ℂ} (h : NV z f v) : ev f z = v := by
  obtain ⟨p, q, hq, rfl, rfl⟩ := h
  have hq0 := ne_of_eval q z hq
  set f := algebraMap ℂ[X] (RatFunc ℂ) p / algebraMap _ _ q with hf
  obtain ⟨t, ht⟩ : f.denom ∣ q := (RatFunc.denom_dvd hq0).2 ⟨p, rfl⟩
  have hd : f.denom.eval z ≠ 0 := by
    intro h0; apply hq; rw [ht, eval_mul, h0, zero_mul]
  have hcross : p * f.denom = f.num * q := by
    apply RatFunc.algebraMap_injective ℂ
    have h1 : algebraMap ℂ[X] (RatFunc ℂ) f.num / algebraMap ℂ[X] (RatFunc ℂ) f.denom
        = algebraMap ℂ[X] (RatFunc ℂ) p / algebraMap _ _ q := RatFunc.num_div_denom f
    rw [div_eq_div_iff (alg_ne _ f.denom_ne_zero) (alg_ne _ hq0)] at h1
    rw [map_mul, map_mul]
    linear_combination -h1
  unfold ev RatFunc.eval
  rw [eval₂_id, eval₂_id]
  rw [div_eq_div_iff hd hq]
  have := congrArg (eval z) hcross
  simp only [eval_mul] at this
  linear_combination -this

lemma NV_add {z : ℂ} {f g : RatFunc ℂ} {v w : ℂ} (hf : NV z f v) (hg : NV z g w) :
    NV z (f + g) (v + w) := by
  obtain ⟨p1, q1, h1, rfl, rfl⟩ := hf
  obtain ⟨p2, q2, h2, rfl, rfl⟩ := hg
  refine ⟨p1 * q2 + p2 * q1, q1 * q2, by simp [h1, h2], ?_, ?_⟩
  · have a1 := alg_ne _ (ne_of_eval q1 z h1)
    have a2 := alg_ne _ (ne_of_eval q2 z h2)
    rw [map_add, map_mul, map_mul, map_mul]
    field_simp
  · simp only [eval_add, eval_mul]; field_simp

lemma NV_mul {z : ℂ} {f g : RatFunc ℂ} {v w : ℂ} (hf : NV z f v) (hg : NV z g w) :
    NV z (f * g) (v * w) := by
  obtain ⟨p1, q1, h1, rfl, rfl⟩ := hf
  obtain ⟨p2, q2, h2, rfl, rfl⟩ := hg
  refine ⟨p1 * p2, q1 * q2, by simp [h1, h2], ?_, ?_⟩
  · rw [map_mul, map_mul]; rw [div_mul_div_comm]
  · simp only [eval_mul]; field_simp

lemma NV_cst (z a : ℂ) : NV z (cst a) a :=
  ⟨C a, 1, by simp, by simp [cst, RatFunc.algebraMap_C], by simp⟩

lemma NV_neg {z : ℂ} {f : RatFunc ℂ} {v : ℂ} (hf : NV z f v) : NV z (-f) (-v) := by
  have := NV_mul (NV_cst z (-1)) hf
  simp only [cst, map_neg, map_one, neg_one_mul] at this
  exact this

lemma NV_sub {z : ℂ} {f g : RatFunc ℂ} {v w : ℂ} (hf : NV z f v) (hg : NV z g w) :
    NV z (f - g) (v - w) := by
  rw [sub_eq_add_neg, sub_eq_add_neg]; exact NV_add hf (NV_neg hg)

lemma NV_sum {ι : Type*} (s : Finset ι) (z : ℂ) (F : ι → RatFunc ℂ) (V : ι → ℂ)
    (h : ∀ i ∈ s, NV z (F i) (V i)) : NV z (∑ i ∈ s, F i) (∑ i ∈ s, V i) := by
  classical
  induction s using Finset.induction_on with
  | empty => exact ⟨0, 1, by simp, by simp, by simp⟩
  | insert a s ha ih =>
    rw [Finset.sum_insert ha, Finset.sum_insert ha]
    exact NV_add (h a (Finset.mem_insert_self a s)) (ih (fun i hi => h i (Finset.mem_insert_of_mem hi)))

lemma NV_zpow (z d : ℂ) (hz : z ≠ d) (e : ℤ) : NV z ((vX - cst d) ^ e) ((z - d) ^ e) := by
  have hzd : z - d ≠ 0 := sub_ne_zero.2 hz
  rcases Int.eq_nat_or_neg e with ⟨k, rfl | rfl⟩
  · refine ⟨(X - C d) ^ k, 1, by simp, ?_, by simp⟩
    simp [vX, cst, RatFunc.algebraMap_C]
  · refine ⟨1, (X - C d) ^ k, by simp [hzd], pole_alg d k, by simp⟩

end P01

namespace P01

noncomputable def JV (c : ℂ) (f : RatFunc ℂ) (w : ℂ) : ℂ :=
  ∑ m ∈ Finset.Icc 2 f.denom.natDegree,
    (lcoeff c f (-(m : ℤ)) * (-(1 / ((m : ℂ) - 1)))) *
      ((w - c) ^ (-((m : ℤ) - 1)) - (bC - c) ^ (-((m : ℤ) - 1)))

noncomputable def HV (c : ℂ) (f : RatFunc ℂ) (w : ℂ) : ℂ :=
  ∑ m ∈ Finset.Icc 2 f.denom.natDegree, lcoeff c f (-(m : ℤ)) * (w - c) ^ (-(m : ℤ))

lemma NV_jpart (c : ℂ) (f : RatFunc ℂ) (w : ℂ) (hw : w ≠ c) : NV w (jpart c f) (JV c f w) := by
  unfold jpart JV
  exact NV_sum _ _ _ _ (fun m _ => NV_mul (NV_cst _ _) (NV_sub (NV_zpow w c hw _) (NV_cst _ _)))

lemma NV_ppart (c : ℂ) (f : RatFunc ℂ) (w : ℂ) (hw : w ≠ c) : NV w (ppart2 c f) (HV c f w) := by
  unfold ppart2 HV
  exact NV_sum _ _ _ _ (fun m _ => NV_mul (NV_cst _ _) (NV_zpow w c hw _))

lemma term_deriv (a c z : ℂ) (m : ℕ) (hm : 2 ≤ m) (hz : z ≠ c) :
    HasDerivAt (fun w => (a * (-(1 / ((m : ℂ) - 1)))) *
      ((w - c) ^ (-((m : ℤ) - 1)) - (bC - c) ^ (-((m : ℤ) - 1)))) (a * (z - c) ^ (-(m : ℤ))) z := by
  have hzc : z - c ≠ 0 := sub_ne_zero.2 hz
  have h1 := hasDerivAt_zpow (-((m : ℤ) - 1)) (z - c) (Or.inl hzc)
  have h2 := h1.comp z ((hasDerivAt_id z).sub_const c)
  have h3 := (h2.sub_const ((bC - c) ^ (-((m : ℤ) - 1)))).const_mul (a * (-(1 / ((m : ℂ) - 1))))
  have hm1 : (m : ℂ) - 1 ≠ 0 := by
    intro h
    have : ((m : ℝ) : ℂ) = 1 := by push_cast; linear_combination h
    have : (m : ℝ) = 1 := by exact_mod_cast this
    have : (2 : ℝ) ≤ m := by exact_mod_cast hm
    linarith
  convert h3 using 1
  · rfl
  rw [show -((m : ℤ) - 1) - 1 = -(m : ℤ) by ring]
  push_cast
  field_simp

lemma b_ne_I : bC ≠ Complex.I := by
  intro h; have := congrArg Complex.im h; simp [bC] at this

lemma b_ne_negI : bC ≠ -Complex.I := by
  intro h; have := congrArg Complex.im h; simp [bC] at this

lemma Jp_b (f : RatFunc ℂ) : ev (Jp f) bC = 0 := by
  unfold Jp
  rw [ev_of_NV (NV_add (NV_jpart _ f _ b_ne_I) (NV_jpart _ f _ b_ne_negI))]
  simp [JV]

lemma Jp_deriv (f : RatFunc ℂ) (z : ℂ) (h1 : z ≠ Complex.I) (h2 : z ≠ -Complex.I) :
    HasDerivAt (fun w => ev (Jp f) w) (ev (Hp f) z) z := by
  have hH : ev (Hp f) z = HV Complex.I f z + HV (-Complex.I) f z := by
    unfold Hp; exact ev_of_NV (NV_add (NV_ppart _ f _ h1) (NV_ppart _ f _ h2))
  rw [hH]
  have hG : HasDerivAt (fun w => JV Complex.I f w + JV (-Complex.I) f w)
      (HV Complex.I f z + HV (-Complex.I) f z) z := by
    unfold JV HV
    apply HasDerivAt.add
    · exact HasDerivAt.fun_sum (fun m hm => term_deriv _ _ z m (Finset.mem_Icc.1 hm).1 h1)
    · exact HasDerivAt.fun_sum (fun m hm => term_deriv _ _ z m (Finset.mem_Icc.1 hm).1 h2)
  apply hG.congr_of_eventuallyEq
  filter_upwards [(eventually_ne_nhds h1).and (eventually_ne_nhds h2)] with w hw
  unfold Jp
  exact ev_of_NV (NV_add (NV_jpart _ f _ hw.1) (NV_jpart _ f _ hw.2))

end P01



namespace P01

open Polynomial

lemma alg_apply (φ : RatFunc ℂ →ₐ[ℂ] RatFunc ℂ) (p : ℂ[X]) :
    φ (algebraMap ℂ[X] (RatFunc ℂ) p) = aeval (φ RatFunc.X) p := by
  rw [← RatFunc.aeval_X_left_eq_algebraMap, Polynomial.aeval_algHom_apply]

lemma algHom_eq (φ ψ : RatFunc ℂ →ₐ[ℂ] RatFunc ℂ) (h : φ RatFunc.X = ψ RatFunc.X) (f : RatFunc ℂ) :
    φ f = ψ f := by
  rw [← RatFunc.num_div_denom f, map_div₀, map_div₀, alg_apply, alg_apply, alg_apply, alg_apply, h]

/-- `S f = f(-1/v)`. -/
noncomputable abbrev Sg : RatFunc ℂ →ₐ[ℂ] RatFunc ℂ := sig (-(RatFunc.X : RatFunc ℂ)⁻¹) hnXinv
/-- `T f = f(1/v)`. -/
noncomputable abbrev Tg : RatFunc ℂ →ₐ[ℂ] RatFunc ℂ := sig ((RatFunc.X : RatFunc ℂ)⁻¹) hXinv

lemma Inv_eq (f : RatFunc ℂ) : Inv f = Sg f / RatFunc.X ^ 2 := by
  unfold Inv; simp only [vX]; rw [subst_eq _ hnXinv]

lemma Fold_eq (f : RatFunc ℂ) : Fold f = Tg f / RatFunc.X ^ 2 := by
  unfold Fold; simp only [vX]; rw [subst_eq _ hXinv]

lemma Sg_X : Sg RatFunc.X = -(RatFunc.X : RatFunc ℂ)⁻¹ := sig_X _ _
lemma Tg_X : Tg RatFunc.X = (RatFunc.X : RatFunc ℂ)⁻¹ := sig_X _ _

lemma SS (f : RatFunc ℂ) : Sg (Sg f) = f := by
  have := algHom_eq (Sg.comp Sg) (AlgHom.id ℂ _) (by
    simp only [AlgHom.comp_apply, Sg_X, map_neg, map_inv₀, AlgHom.id_apply]
    rw [neg_inv, neg_neg, inv_inv])
    f
  simpa using this

lemma ST (f : RatFunc ℂ) : Sg (Tg f) = Tg (Sg f) := by
  have := algHom_eq (Sg.comp Tg) (Tg.comp Sg) (by
    simp only [AlgHom.comp_apply, Sg_X, Tg_X, map_neg, map_inv₀]
    exact neg_inv.symm) f
  simpa using this

lemma Inv_add (f g : RatFunc ℂ) : Inv (f + g) = Inv f + Inv g := by
  rw [Inv_eq, Inv_eq, Inv_eq, map_add, add_div]

lemma Inv_cst_mul (c : ℂ) (f : RatFunc ℂ) : Inv (cst c * f) = cst c * Inv f := by
  rw [Inv_eq, Inv_eq, map_mul, cst, sig_C, mul_div_assoc]

lemma Inv_Inv (f : RatFunc ℂ) : Inv (Inv f) = f := by
  rw [Inv_eq, Inv_eq, map_div₀, SS, map_pow, Sg_X]
  have := X_ne_zero'
  field_simp

lemma Inv_Fold (f : RatFunc ℂ) : Inv (Fold f) = Fold (Inv f) := by
  rw [Inv_eq, Fold_eq, Fold_eq, Inv_eq, map_div₀, map_div₀, map_pow, map_pow, Sg_X, Tg_X, ST]
  have := X_ne_zero'
  field_simp

lemma Inv_D (p q : ℂ) (hpq : p * q = 1) :
    Inv ((vX ^ 2 + 1) / (vX * (vX - cst p) * (vX + cst q))) = (vX ^ 2 + 1) / (vX * (vX - cst p) * (vX + cst q)) := by
  rw [Inv_eq]
  simp only [vX, cst, map_div₀, map_mul, map_add, map_sub, map_pow, map_one, Sg_X, sig_C]
  have hX := X_ne_zero'
  have e1 : (-(RatFunc.X : RatFunc ℂ)⁻¹) ^ 2 + 1 = (1 + RatFunc.X ^ 2) / RatFunc.X ^ 2 := by
    field_simp
    try ring
  have e2 : -(RatFunc.X : RatFunc ℂ)⁻¹ - RatFunc.C p = -(1 + RatFunc.C p * RatFunc.X) / RatFunc.X := by
    field_simp
    try ring
  have e3 : -(RatFunc.X : RatFunc ℂ)⁻¹ + RatFunc.C q = -(1 - RatFunc.C q * RatFunc.X) / RatFunc.X := by
    field_simp
    try ring
  have hpq' : RatFunc.C p * RatFunc.C q = (1 : RatFunc ℂ) := by rw [← map_mul, hpq, map_one]
  have hp : RatFunc.C p ≠ (0 : RatFunc ℂ) := by
    intro h; rw [h, zero_mul] at hpq'; exact zero_ne_one hpq'
  have hq : RatFunc.C q ≠ (0 : RatFunc ℂ) := by
    intro h; rw [h, mul_zero] at hpq'; exact zero_ne_one hpq'
  have f1 : 1 + RatFunc.C p * RatFunc.X = RatFunc.C p * (RatFunc.X + RatFunc.C q) := by
    linear_combination (-1 : RatFunc ℂ) * hpq'
  have f2 : 1 - RatFunc.C q * RatFunc.X = -RatFunc.C q * (RatFunc.X - RatFunc.C p) := by
    linear_combination (-1 : RatFunc ℂ) * hpq'
  have h3 := X_sub_C_ne p
  have h4 : (RatFunc.X : RatFunc ℂ) + RatFunc.C q ≠ 0 := by
    have := X_sub_C_ne (-q); rwa [map_neg, sub_neg_eq_add] at this
  have h5 := X_sub_C_ne q
  have h6 : (RatFunc.X : RatFunc ℂ) + RatFunc.C p ≠ 0 := by
    have := X_sub_C_ne (-p); rwa [map_neg, sub_neg_eq_add] at this
  rw [e1, e2, e3, f1, f2]
  rw [div_div, div_eq_div_iff (by
    refine mul_ne_zero (mul_ne_zero (mul_ne_zero (by simpa using hX) ?_) ?_) (pow_ne_zero _ hX)
    · exact div_ne_zero (neg_ne_zero.2 (mul_ne_zero hp h4)) hX
    · exact div_ne_zero (neg_ne_zero.2 (mul_ne_zero (neg_ne_zero.2 hq) h3)) hX)
    (mul_ne_zero (mul_ne_zero hX h3) h4)]
  field_simp
  linear_combination (-(1 + (RatFunc.X : RatFunc ℂ) ^ 2)) * hpq'

end P01

namespace P01

open Polynomial

lemma Inv_pole (c : ℂ) (hc : c ^ 2 = -1) (m : ℕ) :
    Inv ((vX - cst c) ^ (-(m : ℤ))) = cst (c ^ m) * vX ^ ((m : ℤ) - 2) / (vX - cst c) ^ m := by
  have hX := X_ne_zero'
  have key : -(RatFunc.X : RatFunc ℂ)⁻¹ - RatFunc.C c = RatFunc.C (-c) * (RatFunc.X - RatFunc.C c) * RatFunc.X⁻¹ := by
    have hcc : RatFunc.C c * RatFunc.C c = (-1 : RatFunc ℂ) := by
      rw [← map_mul, ← sq, hc, map_neg, map_one]
    rw [map_neg]
    field_simp
    linear_combination (-1 : RatFunc ℂ) * hcc
  have hcm : (-c) ^ (-(m : ℤ)) = c ^ m := by
    rw [zpow_neg, zpow_natCast, inv_eq_iff_eq_inv, ← inv_pow]
    congr 1
    have hc0 : c ≠ 0 := by rintro rfl; norm_num at hc
    field_simp
    linear_combination -hc
  have hXc := X_sub_C_ne c
  rw [Inv_eq]
  simp only [vX, cst, map_zpow₀, map_sub, Sg_X, sig_C]
  rw [key, mul_zpow, mul_zpow, ← map_zpow₀, hcm, zpow_sub₀ hX, inv_zpow', neg_neg]
  simp only [zpow_neg, zpow_natCast, zpow_ofNat]
  field_simp

lemma Inv_sum {ι : Type*} (s : Finset ι) (F : ι → RatFunc ℂ) :
    Inv (∑ i ∈ s, F i) = ∑ i ∈ s, Inv (F i) := by
  classical
  induction s using Finset.induction_on with
  | empty => simp [Inv_eq]
  | insert a s ha ih => rw [Finset.sum_insert ha, Finset.sum_insert ha, Inv_add, ih]

/-- reduced local form at `c`. -/
lemma red_form (c : ℂ) (f : RatFunc ℂ) : ∃ k : ℕ, ∃ p q : ℂ[X], q.eval c ≠ 0 ∧
    f = (vX - cst c) ^ (-(k : ℤ)) * (algebraMap ℂ[X] (RatFunc ℂ) p / algebraMap _ _ q) ∧
    (1 ≤ k → p.eval c ≠ 0) ∧ (f.denom.eval c = 0 → 1 ≤ k) := by
  obtain ⟨q, hq, hnd⟩ := exists_eq_pow_rootMultiplicity_mul_and_not_dvd f.denom f.denom_ne_zero c
  set k := f.denom.rootMultiplicity c
  have hq0 : q ≠ 0 := by
    rintro rfl; rw [mul_zero] at hq; exact f.denom_ne_zero hq
  have hqc : q.eval c ≠ 0 := by
    intro h; exact hnd (dvd_iff_isRoot.2 h)
  refine ⟨k, f.num, q, hqc, ?_, ?_, ?_⟩
  · conv_lhs => rw [← RatFunc.num_div_denom f, hq]
    rw [pole_alg]
    simp only [map_mul, map_one]
    have h1 : algebraMap ℂ[X] (RatFunc ℂ) ((X - C c) ^ k) ≠ 0 := by
      rw [map_ne_zero_iff _ (RatFunc.algebraMap_injective ℂ)]; exact pow_ne_zero _ (X_sub_C_ne_zero c)
    have h2 : algebraMap ℂ[X] (RatFunc ℂ) q ≠ 0 := by
      rw [map_ne_zero_iff _ (RatFunc.algebraMap_injective ℂ)]; exact hq0
    field_simp
  · intro hk hnum
    have hden : f.denom.eval c = 0 := by
      rw [hq, eval_mul, eval_pow, eval_sub, eval_X, eval_C, sub_self, zero_pow (by omega), zero_mul]
    obtain ⟨u, v, huv⟩ := RatFunc.isCoprime_num_denom f
    have := congrArg (eval c) huv
    simp [hnum, hden] at this
  · intro h0
    exact (rootMultiplicity_pos f.denom_ne_zero).2 h0

lemma NV_reg {c : ℂ} {h : RatFunc ℂ} {v : ℂ} (hh : NV c h v) :
    (∀ n : ℤ, n < 0 → lcoeff c h n = 0) ∧ lcoeff c h 0 = v := by
  obtain ⟨p, q, hq, rfl, rfl⟩ := hh
  refine ⟨fun n hn => ?_, ?_⟩
  · rw [lcoeff_reg c _ _ hq]; simp [hn]
  · rw [lcoeff_reg c _ _ hq]
    simp only [lt_irrefl, ite_false]
    rw [show (0 : ℤ).natAbs = 0 from rfl, ps_coeff0, div_eq_mul_inv]

lemma NV_alg (c : ℂ) (p q : ℂ[X]) (hq : q.eval c ≠ 0) :
    NV c (algebraMap ℂ[X] (RatFunc ℂ) p / algebraMap _ _ q) (p.eval c / q.eval c) :=
  ⟨p, q, hq, rfl, rfl⟩

lemma NV_aeval {c : ℂ} {g : RatFunc ℂ} {w : ℂ} (hg : NV c g w) (p : ℂ[X]) :
    NV c (aeval g p) (p.eval w) := by
  induction p using Polynomial.induction_on with
  | C a => simpa [cst, RatFunc.algebraMap_eq_C] using NV_cst c a
  | add p q hp hq => simpa using NV_add hp hq
  | monomial n a ih =>
    have := NV_mul ih hg
    simpa [pow_succ, mul_assoc] using this

lemma NV_div {z : ℂ} {f g : RatFunc ℂ} {v w : ℂ} (hf : NV z f v) (hg : NV z g w) (hw : w ≠ 0) :
    NV z (f / g) (v / w) := by
  obtain ⟨p1, q1, h1, rfl, rfl⟩ := hf
  obtain ⟨p2, q2, h2, rfl, rfl⟩ := hg
  have hp2 : p2.eval z ≠ 0 := by intro h; apply hw; rw [h, zero_div]
  refine ⟨p1 * q2, q1 * p2, by simp [h1, hp2], ?_, ?_⟩
  · have a1 := alg_ne _ (ne_of_eval q1 z h1)
    have a2 := alg_ne _ (ne_of_eval q2 z h2)
    have a3 := alg_ne _ (ne_of_eval p2 z hp2)
    rw [map_mul, map_mul]
    field_simp
  · simp only [eval_mul]; field_simp

lemma NV_X_zpow (z : ℂ) (hz : z ≠ 0) (e : ℤ) : NV z ((RatFunc.X : RatFunc ℂ) ^ e) (z ^ e) := by
  have := NV_zpow z 0 hz e
  simpa [vX, cst] using this

lemma NV_Sg {c : ℂ} (hc : c ^ 2 = -1) {h : RatFunc ℂ} {v : ℂ} (hh : NV c h v) : NV c (Sg h) v := by
  have hc0 : c ≠ 0 := by rintro rfl; norm_num at hc
  have hSX : NV c (Sg RatFunc.X) c := by
    rw [Sg_X]
    have := NV_neg (NV_X_zpow c hc0 (-1))
    have e : -c ^ (-1 : ℤ) = c := by
      rw [zpow_neg_one]; field_simp; linear_combination -hc
    rw [e] at this
    simpa using this
  obtain ⟨p, q, hq, rfl, rfl⟩ := hh
  rw [map_div₀, alg_apply, alg_apply]
  exact NV_div (NV_aeval hSX p) (NV_aeval hSX q) hq

end P01

namespace P01

open Polynomial

lemma Inv_local (c : ℂ) (hc : c ^ 2 = -1) (k : ℕ) (h : RatFunc ℂ) (v : ℂ) (hh : NV c h v) :
    ∃ H : RatFunc ℂ, Inv ((vX - cst c) ^ (-(k : ℤ)) * h) = (vX - cst c) ^ (-(k : ℤ)) * H ∧
      NV c H ((-c) ^ (-(k : ℤ)) * c ^ k * v / c ^ 2) := by
  have hc0 : c ≠ 0 := by rintro rfl; norm_num at hc
  have hX := X_ne_zero'
  have key : -(RatFunc.X : RatFunc ℂ)⁻¹ - RatFunc.C c = RatFunc.C (-c) * (RatFunc.X - RatFunc.C c) * RatFunc.X⁻¹ := by
    have hcc : RatFunc.C c * RatFunc.C c = (-1 : RatFunc ℂ) := by
      rw [← map_mul, ← sq, hc, map_neg, map_one]
    rw [map_neg]
    field_simp
    linear_combination (-1 : RatFunc ℂ) * hcc
  refine ⟨cst ((-c) ^ (-(k : ℤ))) * RatFunc.X ^ k * Sg h / RatFunc.X ^ 2, ?_, ?_⟩
  · rw [Inv_eq]
    simp only [vX, cst, map_mul, map_zpow₀, map_sub, Sg_X, sig_C]
    rw [key, mul_zpow, mul_zpow, ← map_zpow₀, inv_zpow', neg_neg, zpow_natCast]
    ring
  · refine NV_div (NV_mul (NV_mul (NV_cst _ _) ?_) (NV_Sg hc hh)) ?_ (pow_ne_zero _ hc0)
    · simpa using NV_X_zpow c hc0 k
    · simpa using NV_X_zpow c hc0 2

lemma lcoeff_red (c : ℂ) (k : ℕ) (h : RatFunc ℂ) (v : ℂ) (hh : NV c h v) (n : ℤ) :
    lcoeff c ((vX - cst c) ^ (-(k : ℤ)) * h) n = if n + k < 0 then 0 else lcoeff c h (n + k) := by
  rw [lcoeff_shift]
  split_ifs with hn
  · exact (NV_reg hh).1 _ hn
  · rfl

lemma Inv_pole1 (c : ℂ) (hc : c ^ 2 = -1) (F : RatFunc ℂ) (hF : poleLE c F 1) :
    poleLE c (Inv F) 1 ∧ res c (Inv F) = res c F := by
  obtain ⟨k, p, q, hq, rfl, hk, -⟩ := red_form c F
  have hh := NV_alg c p q hq
  have hk1 : k ≤ 1 := by
    by_contra hk2
    have := hF (-(k : ℤ)) (by omega)
    rw [lcoeff_red c k _ _ hh] at this
    simp only [neg_add_cancel, lt_irrefl, ite_false] at this
    rw [(NV_reg hh).2] at this
    exact (div_ne_zero (hk (by omega)) hq) this
  obtain ⟨H, hH, hHv⟩ := Inv_local c hc k _ _ hh
  rw [hH]
  refine ⟨fun n hn => ?_, ?_⟩
  · rw [lcoeff_red c k _ _ hHv]
    rw [if_pos (by omega)]
  · unfold res
    rw [lcoeff_red c k _ _ hHv, lcoeff_red c k _ _ hh]
    rcases Nat.eq_zero_or_pos k with h0 | h0
    · subst h0; simp
    · have : k = 1 := by omega
      subst this
      simp only [Nat.cast_one, neg_add_cancel, lt_irrefl, ite_false]
      rw [(NV_reg hHv).2, (NV_reg hh).2]
      have hc0 : c ≠ 0 := by rintro rfl; norm_num at hc
      rw [Nat.cast_one, zpow_neg_one, pow_one]
      have e : (-c)⁻¹ * c = -1 := by rw [← neg_inv, neg_mul, inv_mul_cancel₀ hc0]
      rw [e, hc]; ring

end P01

namespace P01

open Polynomial

lemma NV_Xpow (z : ℂ) (n : ℕ) : NV z ((RatFunc.X : RatFunc ℂ) ^ n) (z ^ n) :=
  ⟨X ^ n, 1, by simp, by simp, by simp⟩

lemma vX_zpow (m : ℕ) (hm : 2 ≤ m) : vX ^ ((m : ℤ) - 2) = algebraMap ℂ[X] (RatFunc ℂ) (X ^ (m - 2)) := by
  rw [show (m : ℤ) - 2 = ((m - 2 : ℕ) : ℤ) by omega, zpow_natCast, map_pow, RatFunc.algebraMap_X]; rfl

/-- a term of `𝓘(Hp f)`. -/
noncomputable def tm (c a : ℂ) (m : ℕ) : RatFunc ℂ :=
  cst a * (cst (c ^ m) * vX ^ ((m : ℤ) - 2) / (vX - cst c) ^ m)

lemma Inv_ppart (c : ℂ) (hc : c ^ 2 = -1) (f : RatFunc ℂ) :
    Inv (ppart2 c f) = ∑ m ∈ Finset.Icc 2 f.denom.natDegree, tm c (lcoeff c f (-(m : ℤ))) m := by
  unfold ppart2 tm
  rw [Inv_sum]
  refine Finset.sum_congr rfl (fun m _ => ?_)
  rw [Inv_cst_mul, Inv_pole c hc]

lemma tm_NV (c c' a : ℂ) (h : c' ≠ c) (m : ℕ) (hm : 2 ≤ m) :
    NV c' (tm c a m) (a * (c ^ m * c' ^ (m - 2) / (c' - c) ^ m)) := by
  unfold tm
  rw [vX_zpow m hm, map_pow, RatFunc.algebraMap_X]
  refine NV_mul (NV_cst _ _) (NV_div (NV_mul (NV_cst _ _) (NV_Xpow _ _)) ?_ (pow_ne_zero _ (sub_ne_zero.2 h)))
  have := NV_zpow c' c h m
  simpa using this

lemma tm_small (c a : ℂ) (hc : c ≠ 0) (m : ℕ) (hm : 2 ≤ m) : Small (tm c a m) := by
  by_cases ha : a = 0
  · left; simp [tm, ha, cst]
  right
  unfold tm
  rw [vX_zpow m hm]
  have h1 : cst a ≠ 0 := by unfold cst; simpa using ha
  have h2 : cst (c ^ m) ≠ 0 := by unfold cst; simpa using pow_ne_zero m hc
  have h3 : algebraMap ℂ[X] (RatFunc ℂ) (X ^ (m - 2)) ≠ 0 := alg_ne _ (pow_ne_zero _ X_ne_zero)
  have h4 : (vX - cst c) ^ m ≠ 0 := pow_ne_zero _ (X_sub_C_ne c)
  rw [div_eq_mul_inv, RatFunc.intDegree_mul h1 (mul_ne_zero (mul_ne_zero h2 h3) (inv_ne_zero h4)),
    RatFunc.intDegree_mul (mul_ne_zero h2 h3) (inv_ne_zero h4), RatFunc.intDegree_mul h2 h3,
    RatFunc.intDegree_inv]
  have e4 : (vX - cst c) ^ m = algebraMap ℂ[X] (RatFunc ℂ) ((X - C c) ^ m) := by
    simp [vX, cst, RatFunc.algebraMap_C]
  rw [e4, RatFunc.intDegree_polynomial, RatFunc.intDegree_polynomial]
  simp only [natDegree_pow, natDegree_X_sub_C, natDegree_X, mul_one]
  unfold cst
  rw [RatFunc.intDegree_C, RatFunc.intDegree_C]
  omega

lemma lcoeff_poly (c : ℂ) (P : ℂ[X]) (j : ℕ) :
    lcoeff c (algebraMap ℂ[X] (RatFunc ℂ) P) j = (taylor c P).coeff j := by
  have := lcoeff_reg c P 1 (by simp) j
  rw [map_one, div_one] at this
  rw [this]
  simp only [show ¬ ((j : ℤ) < 0) by omega, ite_false, Int.natAbs_natCast]
  rw [taylor_one, map_one, Polynomial.coe_one, inv_one, mul_one, Polynomial.coeff_coe]

lemma tm_res (c a : ℂ) (m : ℕ) (hm : 2 ≤ m) : res c (tm c a m) = 0 := by
  have e : tm c a m = (vX - cst c) ^ (-(m : ℤ)) * (cst (a * c ^ m) * algebraMap ℂ[X] (RatFunc ℂ) (X ^ (m - 2))) := by
    unfold tm
    rw [vX_zpow m hm, zpow_neg, zpow_natCast]
    unfold cst; rw [map_mul]; ring
  unfold res
  rw [e, lcoeff_shift, lcoeff_cst_mul]
  rw [show (-1 : ℤ) + m = ((m - 1 : ℕ) : ℤ) by omega, lcoeff_poly]
  rw [coeff_eq_zero_of_natDegree_lt, mul_zero]
  rw [natDegree_taylor, natDegree_pow, natDegree_X]
  omega

lemma no_poles_zero (E : RatFunc ℂ) (h : ∀ c, poleLE c E 0) (hs : Small E) : E = 0 := by
  have hroot : ∀ c, E.denom.eval c ≠ 0 := by
    intro c h0
    obtain ⟨k, p, q, hq, hE, hk, hk2⟩ := red_form c E
    have hk1 := hk2 h0
    obtain ⟨-, hv⟩ := NV_reg (NV_alg c p q hq)
    have := h c (-(k : ℤ)) (by omega)
    rw [hE, lcoeff_shift, neg_add_cancel, hv] at this
    exact div_ne_zero (hk hk1) hq this
  have hdeg : E.denom.natDegree = 0 := by
    by_contra hd
    obtain ⟨c, hc⟩ := IsAlgClosed.exists_root E.denom (by
      rw [Polynomial.degree_eq_natDegree E.denom_ne_zero]; exact_mod_cast hd)
    exact hroot c hc
  have hd1 : E.denom = 1 := Polynomial.eq_one_of_monic_natDegree_zero E.monic_denom hdeg
  rcases hs with hs | hs
  · exact hs
  exfalso
  have : E = algebraMap ℂ[X] (RatFunc ℂ) E.num := by
    conv_lhs => rw [← RatFunc.num_div_denom E, hd1, map_one, div_one]
  rw [this, RatFunc.intDegree_polynomial] at hs
  omega

end P01

namespace P01

open Polynomial

lemma sqI (c : ℂ) (hc : c = Complex.I ∨ c = -Complex.I) : c ^ 2 = -1 := by
  rcases hc with rfl | rfl
  · exact Complex.I_sq
  · rw [neg_sq, Complex.I_sq]

lemma I_ne_zero' : Complex.I ≠ 0 := Complex.I_ne_zero
lemma negI_ne_zero : -Complex.I ≠ 0 := neg_ne_zero.2 Complex.I_ne_zero

lemma Inv_Hp_eq (f : RatFunc ℂ) :
    Inv (Hp f) = (∑ m ∈ Finset.Icc 2 f.denom.natDegree, tm Complex.I (lcoeff Complex.I f (-(m : ℤ))) m)
      + ∑ m ∈ Finset.Icc 2 f.denom.natDegree, tm (-Complex.I) (lcoeff (-Complex.I) f (-(m : ℤ))) m := by
  unfold Hp
  rw [Inv_add, Inv_ppart _ (sqI _ (Or.inl rfl)), Inv_ppart _ (sqI _ (Or.inr rfl))]

lemma sum_reg (c d : ℂ) (hcd : c ≠ d) (s : Finset ℕ) (A : ℕ → ℂ) (hs : ∀ m ∈ s, 2 ≤ m) (n : ℤ) (hn : n < 0) :
    lcoeff c (∑ m ∈ s, tm d (A m) m) n = 0 := by
  rw [lcoeff_sum]
  exact Finset.sum_eq_zero (fun m hm => (NV_reg (tm_NV d c (A m) hcd m (hs m hm))).1 n hn)

lemma sum_res_self (c : ℂ) (s : Finset ℕ) (A : ℕ → ℂ) (hs : ∀ m ∈ s, 2 ≤ m) :
    lcoeff c (∑ m ∈ s, tm c (A m) m) (-1) = 0 := by
  rw [lcoeff_sum]
  exact Finset.sum_eq_zero (fun m hm => tm_res c (A m) m (hs m hm))

lemma Icc_two (N : ℕ) : ∀ m ∈ Finset.Icc 2 N, 2 ≤ m := fun m hm => (Finset.mem_Icc.1 hm).1

lemma P_reg (f : RatFunc ℂ) (c : ℂ) (h1 : c ≠ Complex.I) (h2 : c ≠ -Complex.I) :
    poleLE c (Inv (Hp f)) 0 := by
  intro n hn
  have hn' : n < 0 := by simpa using hn
  rw [Inv_Hp_eq, lcoeff_add, sum_reg _ _ h1 _ _ (Icc_two _) _ hn', sum_reg _ _ h2 _ _ (Icc_two _) _ hn',
    add_zero]

lemma P_small (f : RatFunc ℂ) : Small (Inv (Hp f)) := by
  rw [Inv_Hp_eq]
  exact small_add _ _ (small_sum _ _ (fun m hm => tm_small _ _ I_ne_zero' m (Icc_two _ m hm)))
    (small_sum _ _ (fun m hm => tm_small _ _ negI_ne_zero m (Icc_two _ m hm)))

lemma P_res (f : RatFunc ℂ) (c : ℂ) (hc : c = Complex.I ∨ c = -Complex.I) : res c (Inv (Hp f)) = 0 := by
  unfold res
  rw [Inv_Hp_eq, lcoeff_add]
  rcases hc with rfl | rfl
  · rw [sum_res_self _ _ _ (Icc_two _), sum_reg _ _ I_ne_negI _ _ (Icc_two _) _ (by norm_num), add_zero]
  · rw [sum_res_self _ _ _ (Icc_two _), sum_reg _ _ (Ne.symm I_ne_negI) _ _ (Icc_two _) _ (by norm_num),
      zero_add]

lemma small_neg (x : RatFunc ℂ) (hx : Small x) : Small (-x) := by
  rcases hx with rfl | hx
  · left; simp
  · right; rwa [RatFunc.intDegree_neg]

lemma Hp_Inv (f : RatFunc ℂ) :
    Hp (Inv f) = Inv (Hp f) ∧ res Complex.I (Inv f) = res Complex.I f ∧
      res (-Complex.I) (Inv f) = res (-Complex.I) f := by
  set F := f - Hp f with hFdef
  have hg : Inv f = Inv (Hp f) + Inv F := by
    rw [← Inv_add, hFdef, add_sub_cancel]
  have loc : ∀ c, (c = Complex.I ∨ c = -Complex.I) →
      poleLE c (Inv F) 1 ∧ res c (Inv F) = res c F := fun c hc =>
    Inv_pole1 c (sqI c hc) F (sub_Hp_pole f c hc)
  have hE : Inv (Hp f) - Hp (Inv f) = 0 := by
    apply no_poles_zero
    · intro c n hn
      have hn' : n < 0 := by simpa using hn
      rw [lcoeff_sub]
      by_cases hc : c = Complex.I ∨ c = -Complex.I
      · rcases lt_or_eq_of_le (show n ≤ -1 by omega) with h1 | h1
        · have a1 := (loc c hc).1 n h1
          have a2 := sub_Hp_pole (Inv f) c hc n h1
          rw [lcoeff_sub] at a2
          have a3 : lcoeff c (Inv f) n = lcoeff c (Inv (Hp f)) n := by
            rw [hg, lcoeff_add, a1, add_zero]
          linear_combination a2 - a3
        · subst h1
          have := P_res f c hc
          have := Hp_res (Inv f) c hc
          unfold res at *
          simp_all
      · push Not at hc
        rw [P_reg f c hc.1 hc.2 n hn, Hp_regular (Inv f) c hc.1 hc.2 n hn, sub_zero]
    · rw [sub_eq_add_neg]; exact small_add _ _ (P_small f) (small_neg _ (Hp_small _))
  have hEq : Hp (Inv f) = Inv (Hp f) := (sub_eq_zero.1 hE).symm
  have resf : ∀ c, (c = Complex.I ∨ c = -Complex.I) → res c (Inv f) = res c f := by
    intro c hc
    unfold res at *
    rw [hg, lcoeff_add]
    have h1 := P_res f c hc
    have h2 := (loc c hc).2
    have h3 := Hp_res f c hc
    unfold res at h1 h2 h3
    rw [h1, h2, hFdef, lcoeff_sub, h3]; ring
  exact ⟨hEq, resf _ (Or.inl rfl), resf _ (Or.inr rfl)⟩

end P01



namespace P01

open Polynomial

/-- regular at `c`. -/
def Reg (c : ℂ) (f : RatFunc ℂ) : Prop := ∃ v, NV c f v

lemma poleLE_mono {c : ℂ} {f : RatFunc ℂ} {n m : ℕ} (h : poleLE c f n) (hnm : n ≤ m) : poleLE c f m :=
  fun k hk => h k (by omega)

lemma poleLE_add {c : ℂ} {f g : RatFunc ℂ} {n : ℕ} (hf : poleLE c f n) (hg : poleLE c g n) :
    poleLE c (f + g) n := fun k hk => by rw [lcoeff_add, hf k hk, hg k hk, add_zero]

lemma poleLE_neg {c : ℂ} {f : RatFunc ℂ} {n : ℕ} (hf : poleLE c f n) : poleLE c (-f) n :=
  fun k hk => by rw [lcoeff_neg, hf k hk, neg_zero]

lemma poleLE_sub {c : ℂ} {f g : RatFunc ℂ} {n : ℕ} (hf : poleLE c f n) (hg : poleLE c g n) :
    poleLE c (f - g) n := by rw [sub_eq_add_neg]; exact poleLE_add hf (poleLE_neg hg)

lemma poleLE_zero' (c : ℂ) (n : ℕ) : poleLE c 0 n := fun k _ => lcoeff_zero c k

lemma poleLE_cst_mul {c : ℂ} {f : RatFunc ℂ} {n : ℕ} (a : ℂ) (hf : poleLE c f n) :
    poleLE c (cst a * f) n := fun k hk => by rw [lcoeff_cst_mul, hf k hk, mul_zero]

lemma poleLE_sum {ι : Type*} (s : Finset ι) {c : ℂ} (F : ι → RatFunc ℂ) {n : ℕ}
    (h : ∀ i ∈ s, poleLE c (F i) n) : poleLE c (∑ i ∈ s, F i) n := fun k hk => by
  rw [lcoeff_sum]; exact Finset.sum_eq_zero (fun i hi => h i hi k hk)

lemma poleLE_of_reg {c : ℂ} {f : RatFunc ℂ} (h : Reg c f) (n : ℕ) : poleLE c f n := by
  obtain ⟨v, hv⟩ := h
  intro k hk; exact (NV_reg hv).1 k (by omega)

lemma poleLE_form_of {c : ℂ} {h : RatFunc ℂ} (hh : Reg c h) (k : ℕ) :
    poleLE c ((vX - cst c) ^ (-(k : ℤ)) * h) k := by
  obtain ⟨v, hv⟩ := hh
  intro n hn
  rw [lcoeff_red c k h v hv, if_pos (by omega)]

lemma form_of_poleLE {c : ℂ} {f : RatFunc ℂ} {n : ℕ} (hf : poleLE c f n) :
    ∃ k ≤ n, ∃ h, Reg c h ∧ f = (vX - cst c) ^ (-(k : ℤ)) * h := by
  obtain ⟨k, p, q, hq, hfe, hk, -⟩ := red_form c f
  refine ⟨k, ?_, _, ⟨_, NV_alg c p q hq⟩, hfe⟩
  by_contra hkn
  have := hf (-(k : ℤ)) (by omega)
  rw [hfe, lcoeff_red c k _ _ (NV_alg c p q hq)] at this
  simp only [neg_add_cancel, lt_irrefl, ite_false] at this
  rw [(NV_reg (NV_alg c p q hq)).2] at this
  exact div_ne_zero (hk (by omega)) hq this

lemma Reg_mul {c : ℂ} {f g : RatFunc ℂ} (hf : Reg c f) (hg : Reg c g) : Reg c (f * g) := by
  obtain ⟨v, hv⟩ := hf; obtain ⟨w, hw⟩ := hg; exact ⟨_, NV_mul hv hw⟩

lemma Reg_add {c : ℂ} {f g : RatFunc ℂ} (hf : Reg c f) (hg : Reg c g) : Reg c (f + g) := by
  obtain ⟨v, hv⟩ := hf; obtain ⟨w, hw⟩ := hg; exact ⟨_, NV_add hv hw⟩

lemma Reg_cst (c a : ℂ) : Reg c (cst a) := ⟨_, NV_cst c a⟩

lemma poleLE_mul_reg {c : ℂ} {f g : RatFunc ℂ} {n : ℕ} (hf : poleLE c f n) (hg : Reg c g) :
    poleLE c (f * g) n := by
  obtain ⟨k, hk, h, hh, rfl⟩ := form_of_poleLE hf
  rw [mul_assoc]
  exact poleLE_mono (poleLE_form_of (Reg_mul hh hg) k) hk

lemma zero_factor {c : ℂ} {D : RatFunc ℂ} (hD : NV c D 0) :
    ∃ D', Reg c D' ∧ D = (vX - cst c) * D' := by
  obtain ⟨p, q, hq, rfl, hv⟩ := hD
  have hp : p.IsRoot c := by
    have : p.eval c = 0 := by
      rcases div_eq_zero_iff.1 hv.symm with h | h
      · exact h
      · exact absurd h hq
    exact this
  refine ⟨algebraMap ℂ[X] (RatFunc ℂ) (p /ₘ (X - C c)) / algebraMap _ _ q, ⟨_, NV_alg c _ q hq⟩, ?_⟩
  conv_lhs => rw [← (mul_divByMonic_eq_iff_isRoot).2 hp]
  rw [map_mul, mul_div_assoc]
  simp [vX, cst, RatFunc.algebraMap_C]

lemma poleLE_mul_zero {c : ℂ} {f D : RatFunc ℂ} {n : ℕ} (hf : poleLE c f n) (hD : NV c D 0) :
    poleLE c (f * D) (n - 1) := by
  obtain ⟨D', hD', rfl⟩ := zero_factor hD
  obtain ⟨k, hk, h, hh, rfl⟩ := form_of_poleLE hf
  rcases Nat.eq_zero_or_pos k with h0 | h0
  · subst h0
    apply poleLE_of_reg
    simp only [Nat.cast_zero, neg_zero, zpow_zero, one_mul]
    refine Reg_mul hh (Reg_mul ?_ hD')
    exact ⟨_, NV_sub (by simpa [vX] using NV_Xpow c 1) (NV_cst c c)⟩
  · have hXc : vX - cst c ≠ 0 := X_sub_C_ne c
    have e : (vX - cst c) ^ (-(k : ℤ)) * h * ((vX - cst c) * D')
        = (vX - cst c) ^ (-((k - 1 : ℕ) : ℤ)) * (h * D') := by
      have : (vX - cst c) ^ (-(k : ℤ)) * (vX - cst c) = (vX - cst c) ^ (-((k - 1 : ℕ) : ℤ)) := by
        rw [show (-((k - 1 : ℕ) : ℤ)) = -(k : ℤ) + 1 by omega, zpow_add_one₀ hXc]
      rw [← this]; ring
    rw [e]
    exact poleLE_mono (poleLE_form_of (Reg_mul hh hD') _) (by omega)

end P01

namespace P01

open Polynomial

lemma NV_sig (g : RatFunc ℂ) (hg : ∀ c, g - RatFunc.C c ≠ 0) {c w : ℂ} (hgX : NV c g w)
    {h : RatFunc ℂ} {v : ℂ} (hh : NV w h v) : NV c (sig g hg h) v := by
  obtain ⟨p, q, hq, rfl, rfl⟩ := hh
  rw [map_div₀, alg_apply, alg_apply, sig_X]
  exact NV_div (NV_aeval hgX p) (NV_aeval hgX q) hq

lemma NV_Tg {c : ℂ} (hc : c ^ 2 = -1) {h : RatFunc ℂ} {v : ℂ} (hh : NV (-c) h v) : NV c (Tg h) v := by
  have hc0 : c ≠ 0 := by rintro rfl; norm_num at hc
  have hX : NV c ((RatFunc.X : RatFunc ℂ)⁻¹) (-c) := by
    have := NV_X_zpow c hc0 (-1)
    have e : c ^ (-1 : ℤ) = -c := by
      rw [zpow_neg_one]; field_simp; linear_combination hc
    rw [e] at this
    simpa using this
  exact NV_sig _ _ hX hh

lemma Fold_pole (c : ℂ) (hc : c ^ 2 = -1) (F : RatFunc ℂ) (n : ℕ) (hF : poleLE (-c) F n) :
    poleLE c (Fold F) n := by
  have hc0 : c ≠ 0 := by rintro rfl; norm_num at hc
  obtain ⟨k, hk, h, ⟨v, hv⟩, rfl⟩ := form_of_poleLE hF
  have hX := X_ne_zero'
  have key : (RatFunc.X : RatFunc ℂ)⁻¹ - RatFunc.C (-c) = RatFunc.C c * (RatFunc.X - RatFunc.C c) * RatFunc.X⁻¹ := by
    have hcc : RatFunc.C c * RatFunc.C c = (-1 : RatFunc ℂ) := by
      rw [← map_mul, ← sq, hc, map_neg, map_one]
    rw [map_neg]
    field_simp
    linear_combination hcc
  have e : Fold ((vX - cst (-c)) ^ (-(k : ℤ)) * h)
      = (vX - cst c) ^ (-(k : ℤ)) * (cst (c ^ (-(k : ℤ))) * RatFunc.X ^ k * Tg h / RatFunc.X ^ 2) := by
    rw [Fold_eq]
    simp only [vX, cst, map_mul, map_zpow₀, map_sub, Tg_X, sig_C]
    rw [key, mul_zpow, mul_zpow, ← map_zpow₀, inv_zpow', neg_neg, zpow_natCast]
    ring
  rw [e]
  refine poleLE_mono (poleLE_form_of ⟨_, NV_div (NV_mul (NV_mul (NV_cst _ _) (NV_Xpow c k))
    (NV_Tg hc hv)) (NV_Xpow c 2) (pow_ne_zero _ hc0)⟩ k) hk

lemma pole_pow_bound (c : ℂ) (j : ℕ) : poleLE c ((vX - cst c) ^ (-(j : ℤ))) j := by
  have := poleLE_form_of (c := c) (h := 1) ⟨1, by simpa [cst] using NV_cst c 1⟩ j
  simpa using this

lemma Jp_pole (c : ℂ) (hc : c = Complex.I ∨ c = -Complex.I) (W : RatFunc ℂ) (M : ℕ)
    (hW : poleLE c W M) : poleLE c (Jp W) (M - 1) := by
  have self : poleLE c (jpart c W) (M - 1) := by
    unfold jpart
    apply poleLE_sum
    intro m hm
    have hm2 := (Finset.mem_Icc.1 hm).1
    by_cases hmM : m ≤ M
    · apply poleLE_cst_mul
      apply poleLE_sub
      · rw [show -((m : ℤ) - 1) = -((m - 1 : ℕ) : ℤ) by omega]
        exact poleLE_mono (pole_pow_bound c (m - 1)) (by omega)
      · exact poleLE_of_reg (Reg_cst c _) _
    · rw [hW (-(m : ℤ)) (by omega), zero_mul]
      unfold cst; rw [map_zero, zero_mul]
      exact poleLE_zero' c _
  unfold Jp
  rcases hc with rfl | rfl
  · exact poleLE_add self (poleLE_of_reg ⟨_, NV_jpart _ W _ I_ne_negI⟩ _)
  · exact poleLE_add (poleLE_of_reg ⟨_, NV_jpart _ W _ (Ne.symm I_ne_negI)⟩ _) self

lemma Jp_reg (W : RatFunc ℂ) (c : ℂ) (h1 : c ≠ Complex.I) (h2 : c ≠ -Complex.I) : Reg c (Jp W) :=
  ⟨_, NV_add (NV_jpart _ W _ h1) (NV_jpart _ W _ h2)⟩

lemma D1_zero (c : ℂ) (hc : c = Complex.I ∨ c = -Complex.I) : NV c D1 0 := by
  have hQ : (X * (X - C bC) * (X + C aC) : ℂ[X]).eval c ≠ 0 := by
    rcases hc with rfl | rfl
    · have := num1; simp only [eval_mul, eval_X, eval_sub, eval_C, eval_add]; rw [this]
      intro h; have := congrArg Complex.re h; simp at this
    · have := num3; simp only [eval_mul, eval_X, eval_sub, eval_C, eval_add]; rw [this]
      intro h; have := congrArg Complex.re h; simp at this
  have := NV_alg c (X ^ 2 + 1) _ hQ
  rw [← D_alg] at this
  have h0 : (X ^ 2 + 1 : ℂ[X]).eval c = 0 := by simp [sqI c hc]
  rw [h0, zero_div] at this
  exact this

lemma D2_zero (c : ℂ) (hc : c = Complex.I ∨ c = -Complex.I) : NV c D2 0 := by
  have hQ : (X * (X - C aC) * (X + C bC) : ℂ[X]).eval c ≠ 0 := by
    rcases hc with rfl | rfl
    · have := num2; simp only [eval_mul, eval_X, eval_sub, eval_C, eval_add]; rw [this]
      intro h; have := congrArg Complex.re h; simp at this
    · have := num4; simp only [eval_mul, eval_X, eval_sub, eval_C, eval_add]; rw [this]
      intro h; have := congrArg Complex.re h; simp at this
  have := NV_alg c (X ^ 2 + 1) _ hQ
  rw [← D_alg] at this
  have h0 : (X ^ 2 + 1 : ℂ[X]).eval c = 0 := by simp [sqI c hc]
  rw [h0, zero_div] at this
  exact this

end P01

namespace P01

open Polynomial

lemma Reg_zero (c : ℂ) : Reg c 0 := ⟨0, by simpa [cst] using NV_cst c 0⟩
lemma Reg_one (c : ℂ) : Reg c 1 := ⟨1, by simpa [cst] using NV_cst c 1⟩

lemma Reg_ite (c : ℂ) (P : Prop) [Decidable P] {x : RatFunc ℂ} (hx : Reg c x) :
    Reg c (if P then x else 0) := by
  split_ifs
  · exact hx
  · exact Reg_zero c

lemma Reg_sum {ι : Type*} (s : Finset ι) (c : ℂ) (F : ι → RatFunc ℂ) (h : ∀ i ∈ s, Reg c (F i)) :
    Reg c (∑ i ∈ s, F i) := by
  choose! v hv using h
  exact ⟨_, NV_sum s c F v hv⟩

/-- all coefficients of a double power series are regular at `c`. -/
def DReg (c : ℂ) (F : PowerSeries (PowerSeries (RatFunc ℂ))) : Prop :=
  ∀ r s : ℕ, Reg c (PowerSeries.coeff s (PowerSeries.coeff r F))

lemma DReg_add {c : ℂ} {F G : PowerSeries (PowerSeries (RatFunc ℂ))} (hF : DReg c F) (hG : DReg c G) :
    DReg c (F + G) := fun r s => by
  simp only [map_add]; exact Reg_add (hF r s) (hG r s)

lemma DReg_mul {c : ℂ} {F G : PowerSeries (PowerSeries (RatFunc ℂ))} (hF : DReg c F) (hG : DReg c G) :
    DReg c (F * G) := fun r s => by
  rw [PowerSeries.coeff_mul, map_sum]
  apply Reg_sum; intro p _
  rw [PowerSeries.coeff_mul]
  apply Reg_sum; intro q _
  exact Reg_mul (hF _ _) (hG _ _)

lemma DReg_one (c : ℂ) : DReg c 1 := fun r s => by
  rw [PowerSeries.coeff_one]
  split_ifs
  · rw [PowerSeries.coeff_one]; exact Reg_ite c _ (Reg_one c)
  · rw [map_zero]; exact Reg_zero c

lemma DReg_pow {c : ℂ} {F : PowerSeries (PowerSeries (RatFunc ℂ))} (hF : DReg c F) (j : ℕ) :
    DReg c (F ^ j) := by
  induction j with
  | zero => simpa using DReg_one c
  | succ j ih => rw [pow_succ]; exact DReg_mul ih hF

lemma DReg_CC {c : ℂ} {A : RatFunc ℂ} (hA : Reg c A) : DReg c (PowerSeries.C (PowerSeries.C A)) :=
  fun r s => by
  rw [PowerSeries.coeff_C]
  split_ifs
  · rw [PowerSeries.coeff_C]; exact Reg_ite c _ hA
  · rw [map_zero]; exact Reg_zero c

lemma Reg_alg (c x : ℂ) : Reg c (algebraMap ℂ (RatFunc ℂ) x) := by
  rw [RatFunc.algebraMap_eq_C]; exact Reg_cst c x

lemma DReg_E2x (c : ℂ) :
    DReg c (PowerSeries.map ((PowerSeries.C).comp (algebraMap ℂ (RatFunc ℂ))) E2) := fun r s => by
  rw [PowerSeries.coeff_map, RingHom.comp_apply, PowerSeries.coeff_C]
  exact Reg_ite c _ (Reg_alg c _)

lemma DReg_E2y (c : ℂ) :
    DReg c (PowerSeries.C (PowerSeries.map (algebraMap ℂ (RatFunc ℂ)) E2)) := fun r s => by
  rw [PowerSeries.coeff_C]
  split_ifs
  · rw [PowerSeries.coeff_map]; exact Reg_alg c _
  · rw [map_zero]; exact Reg_zero c

lemma tR_alg : tR = algebraMap ℂ[X] (RatFunc ℂ) (X ^ 2 + C 2 * X - 1) / algebraMap ℂ[X] (RatFunc ℂ) (C 4 * X) := by
  unfold tR vX
  simp only [map_add, map_mul, map_pow, map_sub, map_one, RatFunc.algebraMap_X, RatFunc.algebraMap_C]
  rw [map_ofNat, map_ofNat]

lemma NV_tR (c : ℂ) (hc : c ≠ 0) : NV c tR ((c ^ 2 + 2 * c - 1) / (4 * c)) := by
  rw [tR_alg]
  have := NV_alg c (X ^ 2 + C 2 * X - 1) (C 4 * X) (by simp [hc])
  simpa using this

lemma NV_one (c : ℂ) : NV c 1 1 := by simpa [cst] using NV_cst c 1

lemma Reg_Ucoef (c : ℂ) (hc : c ≠ 0) (j r s : ℕ) : Reg c (Ucoef j r s) := by
  have ht : Reg c tR := ⟨_, NV_tR c hc⟩
  have hV : DReg c Vser := by
    unfold Vser
    refine DReg_add (DReg_mul (DReg_CC ?_) (DReg_E2x c)) (DReg_mul (DReg_CC ?_) (DReg_E2y c))
    · rw [sq]; have h1 : Reg c (1 - tR) := ⟨_, NV_sub (NV_one c) (NV_tR c hc)⟩
      exact Reg_mul h1 h1
    · rw [sq]; exact Reg_mul ht ht
  exact DReg_pow hV j r s

end P01

namespace P01

open Polynomial

lemma NV_zpow_gen {c : ℂ} {f : RatFunc ℂ} {v : ℂ} (hf : NV c f v) (hv : v ≠ 0) (e : ℤ) :
    NV c (f ^ e) (v ^ e) := by
  have hpow : ∀ n : ℕ, NV c (f ^ n) (v ^ n) := by
    intro n; induction n with
    | zero => simpa using NV_one c
    | succ n ih => rw [pow_succ, pow_succ]; exact NV_mul ih hf
  rcases Int.eq_nat_or_neg e with ⟨n, rfl | rfl⟩
  · simpa using hpow n
  · have := NV_div (NV_one c) (hpow n) (pow_ne_zero _ hv)
    simpa [zpow_neg] using this

lemma c_ne_zero {c : ℂ} (hc : c ^ 2 = -1) : c ≠ 0 := by rintro rfl; norm_num at hc

lemma Reg_pre (c : ℂ) (hc : c ^ 2 = -1) :
    Reg c (1 / (cst (2 * s2) * vX * tR * (1 - tR))) := by
  have hc0 := c_ne_zero hc
  have hX : NV c vX c := by simpa [vX] using NV_Xpow c 1
  have ht := NV_tR c hc0
  have hnum : c ^ 2 + 2 * c - 1 ≠ 0 := by
    intro h
    have : c = 1 := by linear_combination (h - hc) / 2
    rw [this] at hc; norm_num at hc
  have htc : (c ^ 2 + 2 * c - 1) / (4 * c) ≠ 0 := div_ne_zero hnum (by simp [hc0])
  have h1t : 1 - (c ^ 2 + 2 * c - 1) / (4 * c) ≠ 0 := by
    intro h
    have h' : c ^ 2 + 2 * c - 1 = 4 * c := by
      field_simp at h; linear_combination -h
    have : c = -1 := by linear_combination (h' - hc) / (-2)
    rw [this] at hc; norm_num at hc
  refine ⟨_, NV_div (NV_one c) (NV_mul (NV_mul (NV_mul (NV_cst c _) hX) ht) (NV_sub (NV_one c) ht)) ?_⟩
  exact mul_ne_zero (mul_ne_zero (mul_ne_zero (mul_ne_zero two_ne_zero s2_ne) hc0) htc) h1t

/-- `q = (v - c) Q_c`. -/
noncomputable def Qc (c : ℂ) : RatFunc ℂ := (vX + cst c) / (cst (2 * s2) * vX)

lemma qR_eq (c : ℂ) (hc : c ^ 2 = -1) : qR = (vX - cst c) * Qc c := by
  unfold qR Qc
  rw [← mul_div_assoc]
  congr 1
  have : cst c * cst c = -1 := by unfold cst; rw [← map_mul, ← sq, hc, map_neg, map_one]
  linear_combination this

lemma NV_Qc (c : ℂ) (hc : c ^ 2 = -1) : NV c (Qc c) ((c + c) / (2 * s2 * c)) := by
  have hc0 := c_ne_zero hc
  have hX : NV c vX c := by simpa [vX] using NV_Xpow c 1
  exact NV_div (NV_add hX (NV_cst c c)) (NV_mul (NV_cst c _) hX)
    (mul_ne_zero (mul_ne_zero two_ne_zero s2_ne) hc0)

lemma Qc_val_ne (c : ℂ) (hc : c ^ 2 = -1) : (c + c) / (2 * s2 * c) ≠ 0 := by
  have hc0 := c_ne_zero hc
  exact div_ne_zero (by rw [← two_mul]; exact mul_ne_zero two_ne_zero hc0)
    (mul_ne_zero (mul_ne_zero two_ne_zero s2_ne) hc0)

lemma Rt_pole (c : ℂ) (hc : c ^ 2 = -1) (r s : ℕ) : poleLE c (Rt r s) (2 * (r + s) - 2) := by
  have hc0 := c_ne_zero hc
  unfold Rt
  have hcc : c = Complex.I ∨ c = -Complex.I := by
    have : (c - Complex.I) * (c + Complex.I) = 0 := by
      linear_combination hc - Complex.I_sq
    rcases mul_eq_zero.1 this with h | h
    · left; linear_combination h
    · right; linear_combination h
  have hE1 : Reg c (if s = 0 then -D1 / cst (2 * ((r + s).factorial : ℂ)) else 0) := by
    refine Reg_ite c _ ⟨_, NV_div (NV_neg (D1_zero c hcc)) (NV_cst c _) ?_⟩
    exact mul_ne_zero two_ne_zero (by exact_mod_cast (Nat.factorial_pos _).ne')
  have hE2 : Reg c (if r = 0 then D2 / cst (2 * ((r + s).factorial : ℂ)) else 0) := by
    refine Reg_ite c _ ⟨_, NV_div (D2_zero c hcc) (NV_cst c _) ?_⟩
    exact mul_ne_zero two_ne_zero (by exact_mod_cast (Nat.factorial_pos _).ne')
  refine poleLE_add (poleLE_add ?_ (poleLE_of_reg hE1 _)) (poleLE_of_reg hE2 _)
  rw [mul_comm]
  refine poleLE_mul_reg ?_ (Reg_pre c hc)
  apply poleLE_sum
  intro j hj
  have hj1 := (Finset.mem_Icc.1 hj).1
  have hj2 := (Finset.mem_Icc.1 hj).2
  have e : cst (binHalf j) * qR ^ (2 - 2 * (j : ℤ)) * Ucoef j r s
      = (vX - cst c) ^ (-((2 * j - 2 : ℕ) : ℤ)) *
          (cst (binHalf j) * Qc c ^ (2 - 2 * (j : ℤ)) * Ucoef j r s) := by
    rw [qR_eq c hc, mul_zpow, show (-((2 * j - 2 : ℕ) : ℤ)) = 2 - 2 * (j : ℤ) by omega]
    ring
  rw [e]
  refine poleLE_mono (poleLE_form_of ?_ _) (by omega)
  obtain ⟨u, hu⟩ := Reg_Ucoef c hc0 j r s
  exact ⟨_, NV_mul (NV_mul (NV_cst c _) (NV_zpow_gen (NV_Qc c hc) (Qc_val_ne c hc) _)) hu⟩

lemma Rinit_pole (c : ℂ) (hc : c ^ 2 = -1) (k r s : ℕ) :
    poleLE c (Rinit k r s) (2 * (r + s) - 2) := by
  unfold Rinit
  split_ifs with h1 h2
  · refine poleLE_add (Rt_pole c hc r s) (Fold_pole c hc _ _ ?_)
    have := Rt_pole (-c) (by rw [neg_sq, hc]) s r
    rwa [add_comm s r] at this
  · exact Rt_pole c hc r s
  · exact poleLE_zero' c _

end P01

namespace P01

lemma W_eq (k r s : ℕ) (h : r + s ≤ k) :
    W k r s = Rinit k r s
        - cst ((r : ℂ) + 1) * Jp (W k (r + 1) s) * D1
        - (if s < r then cst ((s : ℂ) + 1) * Jp (W k r (s + 1)) * D2 else 0)
        - (if _h2 : s + 1 = r then cst (r : ℂ) * Fold (Jp (W k r r) * D1) else 0) := by
  rw [W, dif_pos h]

lemma W_gt (k r s : ℕ) (h : k < r + s) : W k r s = 0 := by
  rw [W, dif_neg (by omega)]

lemma sq_of (c : ℂ) (hc : c = Complex.I ∨ c = -Complex.I) : c ^ 2 = -1 := sqI c hc

lemma neg_mem (c : ℂ) (hc : c = Complex.I ∨ c = -Complex.I) : -c = Complex.I ∨ -c = -Complex.I := by
  rcases hc with rfl | rfl
  · right; rfl
  · left; rw [neg_neg]

lemma W_bound (k : ℕ) : ∀ r s : ℕ, s ≤ r → 1 ≤ r + s → ∀ c : ℂ, (c = Complex.I ∨ c = -Complex.I) →
    poleLE c (W k r s) (2 * (r + s) - 2) := by
  suffices H : ∀ n : ℕ, ∀ r s : ℕ, k + 1 - (r + s) = n → s ≤ r → 1 ≤ r + s →
      ∀ c : ℂ, (c = Complex.I ∨ c = -Complex.I) → poleLE c (W k r s) (2 * (r + s) - 2) by
    intro r s h1 h2 c hc; exact H _ r s rfl h1 h2 c hc
  intro n
  induction n using Nat.strong_induction_on with
  | _ n ih =>
  intro r s hn hsr hd c hc
  by_cases hk : r + s ≤ k
  swap
  · rw [W_gt k r s (by omega)]; exact poleLE_zero' c _
  rw [W_eq k r s hk]
  have IH : ∀ r' s', s' ≤ r' → r' + s' = r + s + 1 → ∀ c' : ℂ, (c' = Complex.I ∨ c' = -Complex.I) →
      poleLE c' (W k r' s') (2 * (r + s + 1) - 2) := by
    intro r' s' h1 h2 c' hc'
    have := ih (k + 1 - (r' + s')) (by omega) r' s' rfl h1 (by omega) c' hc'
    rwa [h2] at this
  have par : ∀ r' s', s' ≤ r' → r' + s' = r + s + 1 → ∀ c' : ℂ, (c' = Complex.I ∨ c' = -Complex.I) →
      ∀ D : RatFunc ℂ, NV c' D 0 → ∀ a : ℂ,
      poleLE c' (cst a * Jp (W k r' s') * D) (2 * (r + s) - 2) := by
    intro r' s' h1 h2 c' hc' D hD a
    have h3 := poleLE_mul_zero (poleLE_cst_mul a (Jp_pole c' hc' _ _ (IH r' s' h1 h2 c' hc'))) hD
    exact poleLE_mono h3 (by omega)
  refine poleLE_sub (poleLE_sub (poleLE_sub (Rinit_pole c (sq_of c hc) k r s)
    (par (r + 1) s (by omega) (by omega) c hc D1 (D1_zero c hc) _)) ?_) ?_
  · split_ifs with h
    · exact par r (s + 1) (by omega) (by omega) c hc D2 (D2_zero c hc) _
    · exact poleLE_zero' c _
  · split_ifs with h
    · apply poleLE_cst_mul
      apply Fold_pole c (sq_of c hc)
      have h3 := poleLE_mul_zero (Jp_pole (-c) (neg_mem c hc) _ _
        (IH r r le_rfl (by omega) (-c) (neg_mem c hc))) (D1_zero (-c) (neg_mem c hc))
      exact poleLE_mono h3 (by omega)
    · exact poleLE_zero' c _

end P01



namespace P01

lemma E2_coeff0 : PowerSeries.coeff 0 E2 = 0 := by
  unfold E2
  rw [map_sub, PowerSeries.coeff_rescale, PowerSeries.coeff_exp, PowerSeries.coeff_one]
  simp

lemma E2_coeff1 : PowerSeries.coeff 1 E2 = 2 := by
  unfold E2
  rw [map_sub, PowerSeries.coeff_rescale, PowerSeries.coeff_exp, PowerSeries.coeff_one]
  simp

lemma pow_coeff_X_dvd {R : Type*} [CommRing R] (φ : PowerSeries R) (h0 : PowerSeries.coeff 0 φ = 0)
    (i n : ℕ) :
    PowerSeries.coeff n (φ ^ i) = if n < i then 0 else if n = i then (PowerSeries.coeff 1 φ) ^ i
      else PowerSeries.coeff n (φ ^ i) := by
  obtain ⟨ψ, hψ⟩ : (PowerSeries.X : PowerSeries R) ∣ φ := by
    rw [PowerSeries.X_dvd_iff, ← PowerSeries.coeff_zero_eq_constantCoeff_apply]; exact h0
  split_ifs with h1 h2
  · rw [hψ, mul_pow, PowerSeries.coeff_X_pow_mul', if_neg (by omega)]
  · subst h2
    rw [hψ, mul_pow, PowerSeries.coeff_X_pow_mul', if_pos le_rfl, Nat.sub_self,
      PowerSeries.coeff_zero_eq_constantCoeff_apply, map_pow, ← PowerSeries.coeff_zero_eq_constantCoeff_apply,
      show (1 : ℕ) = 0 + 1 from rfl, PowerSeries.coeff_succ_X_mul]
  · rfl

lemma E2_pow_lt (i n : ℕ) (h : n < i) : PowerSeries.coeff n (E2 ^ i) = 0 := by
  rw [pow_coeff_X_dvd E2 E2_coeff0, if_pos h]

lemma E2_pow_eq (i : ℕ) : PowerSeries.coeff i (E2 ^ i) = 2 ^ i := by
  rw [pow_coeff_X_dvd E2 E2_coeff0, if_neg (lt_irrefl _), if_pos rfl, E2_coeff1]

lemma term_eq (A B : RatFunc ℂ) (i j n : ℕ) :
    (PowerSeries.C (PowerSeries.C A) * PowerSeries.map ((PowerSeries.C).comp (algebraMap ℂ (RatFunc ℂ))) E2) ^ i *
      (PowerSeries.C (PowerSeries.C B) * PowerSeries.C (PowerSeries.map (algebraMap ℂ (RatFunc ℂ)) E2)) ^ j *
      (n : PowerSeries (PowerSeries (RatFunc ℂ)))
    = PowerSeries.C (PowerSeries.C (A ^ i * B ^ j * n)) *
      (PowerSeries.map ((PowerSeries.C).comp (algebraMap ℂ (RatFunc ℂ))) (E2 ^ i) *
        PowerSeries.C (PowerSeries.map (algebraMap ℂ (RatFunc ℂ)) (E2 ^ j))) := by
  simp only [map_mul, map_pow, map_natCast]; ring

lemma term_coeff (α : RatFunc ℂ) (G H : PowerSeries ℂ) (r s : ℕ) :
    PowerSeries.coeff s (PowerSeries.coeff r (PowerSeries.C (PowerSeries.C α) *
      (PowerSeries.map ((PowerSeries.C).comp (algebraMap ℂ (RatFunc ℂ))) G *
        PowerSeries.C (PowerSeries.map (algebraMap ℂ (RatFunc ℂ)) H))))
    = α * (algebraMap ℂ (RatFunc ℂ) (PowerSeries.coeff r G) *
        algebraMap ℂ (RatFunc ℂ) (PowerSeries.coeff s H)) := by
  rw [PowerSeries.coeff_C_mul, PowerSeries.coeff_mul_C, PowerSeries.coeff_map, RingHom.comp_apply,
    PowerSeries.coeff_C_mul, PowerSeries.coeff_C_mul, PowerSeries.coeff_map]

lemma Ucoef_top (r s : ℕ) :
    Ucoef (r + s) r s = cst (((r + s).choose r : ℕ) * 2 ^ (r + s)) * ((1 - tR) ^ 2) ^ r * (tR ^ 2) ^ s := by
  unfold Ucoef Vser
  rw [add_pow, map_sum, map_sum]
  simp only [term_eq, term_coeff]
  rw [Finset.sum_eq_single r]
  · rw [show r + s - r = s by omega, E2_pow_eq, E2_pow_eq]
    unfold cst
    rw [← RatFunc.algebraMap_eq_C]
    simp only [map_mul, map_pow, map_natCast, map_ofNat]
    ring
  · intro i hi hir
    have hi' := Finset.mem_range.1 hi
    rcases lt_or_gt_of_ne hir with h | h
    · rw [E2_pow_lt (r + s - i) s (by omega), map_zero, mul_zero, mul_zero]
    · rw [E2_pow_lt i r h, map_zero, zero_mul, mul_zero]
  · intro h; exact absurd (Finset.mem_range.2 (by omega)) h

end P01

namespace P01

lemma binHalf_succ (j : ℕ) : binHalf (j + 1) = binHalf j * ((1 / 2 : ℂ) - j) / ((j : ℂ) + 1) := by
  unfold binHalf
  rw [Finset.prod_range_succ, Nat.factorial_succ]
  have : (j.factorial : ℂ) ≠ 0 := by exact_mod_cast (Nat.factorial_pos _).ne'
  have h2 : (j : ℂ) + 1 ≠ 0 := by exact_mod_cast Nat.succ_ne_zero j
  push_cast
  field_simp

lemma binHalf_formula (n : ℕ) :
    binHalf (n + 2) = (-1) ^ (n + 1) * ((Nat.doubleFactorial (2 * n + 1) : ℕ) : ℂ) /
      (2 ^ (n + 2) * ((n + 2).factorial : ℂ)) := by
  induction n with
  | zero =>
    rw [binHalf_succ, binHalf_succ]
    unfold binHalf
    simp
    norm_num
  | succ n ih =>
    rw [binHalf_succ, ih]
    rw [show 2 * (n + 1) + 1 = (2 * n + 1) + 2 by ring, Nat.doubleFactorial_add_two,
      Nat.factorial_succ (n + 2)]
    have : ((n + 2).factorial : ℂ) ≠ 0 := by exact_mod_cast (Nat.factorial_pos _).ne'
    have h3 : (n : ℂ) + 2 + 1 ≠ 0 := by norm_cast
    push_cast
    field_simp
    ring

lemma NV_npow {c : ℂ} {f : RatFunc ℂ} {v : ℂ} (hf : NV c f v) (n : ℕ) : NV c (f ^ n) (v ^ n) := by
  induction n with
  | zero => simpa using NV_one c
  | succ n ih => rw [pow_succ, pow_succ]; exact NV_mul ih hf

/-- value of `t` at `c`. -/
noncomputable def tv (c : ℂ) : ℂ := (c ^ 2 + 2 * c - 1) / (4 * c)

lemma NV_pre (c : ℂ) (hc : c ^ 2 = -1) :
    NV c (1 / (cst (2 * s2) * vX * tR * (1 - tR))) (1 / (2 * s2 * c * tv c * (1 - tv c))) := by
  have hc0 := c_ne_zero hc
  have hX : NV c vX c := by simpa [vX] using NV_Xpow c 1
  have ht := NV_tR c hc0
  obtain ⟨v, hv⟩ := Reg_pre c hc
  -- recompute with explicit value
  have hnum : c ^ 2 + 2 * c - 1 ≠ 0 := by
    intro h
    have : c = 1 := by linear_combination (h - hc) / 2
    rw [this] at hc; norm_num at hc
  have htc : (c ^ 2 + 2 * c - 1) / (4 * c) ≠ 0 := div_ne_zero hnum (by simp [hc0])
  have h1t : 1 - (c ^ 2 + 2 * c - 1) / (4 * c) ≠ 0 := by
    intro h
    have h' : c ^ 2 + 2 * c - 1 = 4 * c := by
      field_simp at h; linear_combination -h
    have : c = -1 := by linear_combination (h' - hc) / (-2)
    rw [this] at hc; norm_num at hc
  exact NV_div (NV_one c) (NV_mul (NV_mul (NV_mul (NV_cst c _) hX) ht) (NV_sub (NV_one c) ht))
    (mul_ne_zero (mul_ne_zero (mul_ne_zero (mul_ne_zero two_ne_zero s2_ne) hc0) htc) h1t)

lemma Rt_lead (c : ℂ) (hc : c ^ 2 = -1) (r s : ℕ) (hd : 2 ≤ r + s) :
    lcoeff c (Rt r s) (-(2 * ((r + s : ℕ) : ℤ) - 2)) =
      binHalf (r + s) * ((c + c) / (2 * s2 * c)) ^ (2 - 2 * ((r + s : ℕ) : ℤ)) *
        ((((r + s).choose r : ℕ) * 2 ^ (r + s) : ℂ) * ((1 - tv c) ^ 2) ^ r * (tv c ^ 2) ^ s) *
        (1 / (2 * s2 * c * tv c * (1 - tv c))) := by
  have hc0 := c_ne_zero hc
  have hcc : c = Complex.I ∨ c = -Complex.I := by
    have : (c - Complex.I) * (c + Complex.I) = 0 := by
      linear_combination hc - Complex.I_sq
    rcases mul_eq_zero.1 this with h | h
    · left; linear_combination h
    · right; linear_combination h
  have hneg : (-(2 * ((r + s : ℕ) : ℤ) - 2)) < 0 := by omega
  have ht := NV_tR c hc0
  have hA : NV c ((1 - tR) ^ 2) ((1 - tv c) ^ 2) := NV_npow (NV_sub (NV_one c) ht) 2
  have hB : NV c (tR ^ 2) (tv c ^ 2) := NV_npow ht 2
  unfold Rt
  rw [lcoeff_add, lcoeff_add]
  have hE1 : Reg c (if s = 0 then -D1 / cst (2 * ((r + s).factorial : ℂ)) else 0) := by
    refine Reg_ite c _ ⟨_, NV_div (NV_neg (D1_zero c hcc)) (NV_cst c _) ?_⟩
    exact mul_ne_zero two_ne_zero (by exact_mod_cast (Nat.factorial_pos _).ne')
  have hE2 : Reg c (if r = 0 then D2 / cst (2 * ((r + s).factorial : ℂ)) else 0) := by
    refine Reg_ite c _ ⟨_, NV_div (D2_zero c hcc) (NV_cst c _) ?_⟩
    exact mul_ne_zero two_ne_zero (by exact_mod_cast (Nat.factorial_pos _).ne')
  obtain ⟨_, h1⟩ := hE1
  obtain ⟨_, h2⟩ := hE2
  rw [(NV_reg h1).1 _ hneg, (NV_reg h2).1 _ hneg, add_zero, add_zero]
  rw [Finset.mul_sum, lcoeff_sum]
  have key : ∀ j ∈ Finset.Icc 1 (r + s),
      lcoeff c (1 / (cst (2 * s2) * vX * tR * (1 - tR)) *
        (cst (binHalf j) * qR ^ (2 - 2 * (j : ℤ)) * Ucoef j r s)) (-(2 * ((r + s : ℕ) : ℤ) - 2))
      = if j = r + s then binHalf (r + s) * ((c + c) / (2 * s2 * c)) ^ (2 - 2 * ((r + s : ℕ) : ℤ)) *
        ((((r + s).choose r : ℕ) * 2 ^ (r + s) : ℂ) * ((1 - tv c) ^ 2) ^ r * (tv c ^ 2) ^ s) *
        (1 / (2 * s2 * c * tv c * (1 - tv c))) else 0 := by
    intro j hj
    have hj1 := (Finset.mem_Icc.1 hj).1
    have hj2 := (Finset.mem_Icc.1 hj).2
    have e : 1 / (cst (2 * s2) * vX * tR * (1 - tR)) *
        (cst (binHalf j) * qR ^ (2 - 2 * (j : ℤ)) * Ucoef j r s)
        = (vX - cst c) ^ (-((2 * j - 2 : ℕ) : ℤ)) *
          (cst (binHalf j) * Qc c ^ (2 - 2 * (j : ℤ)) * Ucoef j r s *
            (1 / (cst (2 * s2) * vX * tR * (1 - tR)))) := by
      rw [qR_eq c hc, mul_zpow, show (-((2 * j - 2 : ℕ) : ℤ)) = 2 - 2 * (j : ℤ) by omega]
      ring
    rw [e]
    split_ifs with hjd
    · subst hjd
      have hU : NV c (Ucoef (r + s) r s)
          ((((r + s).choose r : ℕ) * 2 ^ (r + s) : ℂ) * ((1 - tv c) ^ 2) ^ r * (tv c ^ 2) ^ s) := by
        rw [Ucoef_top]
        exact NV_mul (NV_mul (NV_cst c _) (NV_npow hA r)) (NV_npow hB s)
      have hbr := NV_mul (NV_mul (NV_mul (NV_cst c (binHalf (r + s)))
        (NV_zpow_gen (NV_Qc c hc) (Qc_val_ne c hc) (2 - 2 * ((r + s : ℕ) : ℤ)))) hU) (NV_pre c hc)
      rw [lcoeff_red c _ _ _ hbr, if_neg (by omega)]
      rw [show (-(2 * ((r + s : ℕ) : ℤ) - 2)) + ((2 * (r + s) - 2 : ℕ) : ℤ) = 0 by omega]
      rw [(NV_reg hbr).2]
    · obtain ⟨u, hu⟩ := Reg_Ucoef c hc0 j r s
      have hbr := NV_mul (NV_mul (NV_mul (NV_cst c (binHalf j))
        (NV_zpow_gen (NV_Qc c hc) (Qc_val_ne c hc) (2 - 2 * (j : ℤ)))) hu) (NV_pre c hc)
      rw [lcoeff_red c _ _ _ hbr, if_pos (by omega)]
  rw [Finset.sum_congr rfl key, Finset.sum_ite_eq']
  rw [if_pos (Finset.mem_Icc.2 ⟨by omega, le_rfl⟩)]

end P01

namespace P01

lemma Fold_form (c : ℂ) (hc : c ^ 2 = -1) (k : ℕ) (h : RatFunc ℂ) (v : ℂ) (hv : NV (-c) h v) :
    ∃ H, Fold ((vX - cst (-c)) ^ (-(k : ℤ)) * h) = (vX - cst c) ^ (-(k : ℤ)) * H ∧ NV c H (-v) := by
  have hc0 : c ≠ 0 := c_ne_zero hc
  have hX := X_ne_zero'
  have key : (RatFunc.X : RatFunc ℂ)⁻¹ - RatFunc.C (-c) = RatFunc.C c * (RatFunc.X - RatFunc.C c) * RatFunc.X⁻¹ := by
    have hcc : RatFunc.C c * RatFunc.C c = (-1 : RatFunc ℂ) := by
      rw [← map_mul, ← sq, hc, map_neg, map_one]
    rw [map_neg]
    field_simp
    linear_combination hcc
  refine ⟨cst (c ^ (-(k : ℤ))) * RatFunc.X ^ k * Tg h / RatFunc.X ^ 2, ?_, ?_⟩
  · rw [Fold_eq]
    simp only [vX, cst, map_mul, map_zpow₀, map_sub, Tg_X, sig_C]
    rw [key, mul_zpow, mul_zpow, ← map_zpow₀, inv_zpow', neg_neg, zpow_natCast]
    ring
  · have := NV_div (NV_mul (NV_mul (NV_cst _ (c ^ (-(k : ℤ)))) (NV_Xpow c k))
      (NV_Tg hc hv)) (NV_Xpow c 2) (pow_ne_zero _ hc0)
    have e : c ^ (-(k : ℤ)) * c ^ k * v / c ^ 2 = -v := by
      rw [zpow_neg, zpow_natCast, inv_mul_cancel₀ (pow_ne_zero _ hc0), hc]; ring
    rwa [e] at this

lemma Fold_lead (c : ℂ) (hc : c ^ 2 = -1) (F : RatFunc ℂ) (K : ℕ) (hF : poleLE (-c) F K) :
    lcoeff c (Fold F) (-(K : ℤ)) = - lcoeff (-c) F (-(K : ℤ)) := by
  obtain ⟨k, hk, h, ⟨v, hv⟩, rfl⟩ := form_of_poleLE hF
  obtain ⟨H, hH, hHv⟩ := Fold_form c hc k h v hv
  rw [hH, lcoeff_red c k H _ hHv, lcoeff_red (-c) k h v hv]
  split_ifs with h1
  · simp
  · have : -(K : ℤ) + k = 0 := by omega
    rw [this, (NV_reg hHv).2, (NV_reg hv).2]

lemma I_pow_three' : Complex.I ^ 3 = -Complex.I := by rw [pow_succ, Complex.I_sq]; ring

lemma tvI : tv Complex.I = (Complex.I - 1) / (2 * Complex.I) := by
  unfold tv; rw [Complex.I_sq]; field_simp; ring
lemma tvnI : tv (-Complex.I) = (Complex.I + 1) / (2 * Complex.I) := by
  unfold tv; rw [neg_sq, Complex.I_sq]; field_simp; ring

lemma valsI : (1 - tv Complex.I) ^ 2 = -Complex.I / 2 ∧ tv Complex.I ^ 2 = Complex.I / 2 ∧
    1 / (2 * s2 * Complex.I * tv Complex.I * (1 - tv Complex.I)) = -Complex.I / s2 := by
  rw [tvI]
  have := s2_ne
  refine ⟨?_, ?_, ?_⟩
  · field_simp; ring_nf; simp only [Complex.I_sq, I_pow_three']; ring
  · field_simp; ring_nf; simp only [Complex.I_sq, I_pow_three']; ring
  · field_simp; ring_nf; simp only [Complex.I_sq, I_pow_three']; ring

lemma valsnI : (1 - tv (-Complex.I)) ^ 2 = Complex.I / 2 ∧ tv (-Complex.I) ^ 2 = -Complex.I / 2 ∧
    1 / (2 * s2 * (-Complex.I) * tv (-Complex.I) * (1 - tv (-Complex.I))) = Complex.I / s2 := by
  rw [tvnI]
  have := s2_ne
  refine ⟨?_, ?_, ?_⟩
  · field_simp; ring_nf; simp only [Complex.I_sq, I_pow_three']; ring
  · field_simp; ring_nf; simp only [Complex.I_sq, I_pow_three']; ring
  · field_simp; ring_nf; simp only [Complex.I_sq, I_pow_three']; ring

end P01

namespace P01

lemma Qpow (c : ℂ) (hc : c ≠ 0) (d : ℕ) (hd : 1 ≤ d) :
    ((c + c) / (2 * s2 * c)) ^ (2 - 2 * (d : ℤ)) = 2 ^ (d - 1) := by
  have := s2_ne
  have e : (c + c) / (2 * s2 * c) = s2⁻¹ := by field_simp; ring
  rw [e, inv_zpow', show -(2 - 2 * (d : ℤ)) = ((2 * (d - 1) : ℕ) : ℤ) by omega, zpow_natCast, pow_mul,
    s2_sq]

lemma lead_core (n : ℕ) (C : ℂ) (σ : ℂ) :
    binHalf (n + 2) * 2 ^ (n + 1) * (C * 2 ^ (n + 2) * (σ * (Complex.I / 2) ^ (n + 2))) * (-Complex.I / s2)
      = Kd (n + 2) * σ * C := by
  rw [binHalf_formula]
  unfold Kd
  rw [show 2 * (n + 2) - 3 = 2 * n + 1 by omega]
  have := s2_ne
  have hf : ((n + 2).factorial : ℂ) ≠ 0 := by exact_mod_cast (Nat.factorial_pos _).ne'
  rw [div_pow]
  field_simp
  ring

lemma lead_I (r s : ℕ) (h : 2 ≤ r + s) :
    binHalf (r + s) * ((Complex.I + Complex.I) / (2 * s2 * Complex.I)) ^ (2 - 2 * ((r + s : ℕ) : ℤ)) *
        ((((r + s).choose r : ℕ) * 2 ^ (r + s) : ℂ) * ((1 - tv Complex.I) ^ 2) ^ r * (tv Complex.I ^ 2) ^ s) *
        (1 / (2 * s2 * Complex.I * tv Complex.I * (1 - tv Complex.I)))
      = Kd (r + s) * (-1) ^ r * (((r + s).choose r : ℕ) : ℂ) := by
  obtain ⟨hA, hB, hP⟩ := valsI
  rw [hA, hB, hP, Qpow _ Complex.I_ne_zero _ (by omega)]
  have e : (-Complex.I / 2) ^ r * (Complex.I / 2) ^ s = (-1) ^ r * (Complex.I / 2) ^ (r + s) := by
    rw [show -Complex.I / 2 = (-1) * (Complex.I / 2) by ring, mul_pow, mul_assoc, ← pow_add]
  obtain ⟨n, hn⟩ : ∃ n, r + s = n + 2 := ⟨r + s - 2, by omega⟩
  rw [mul_assoc (((r + s).choose r : ℕ) * 2 ^ (r + s) : ℂ), e]
  rw [hn, show n + 2 - 1 = n + 1 by omega]
  exact lead_core n _ _

lemma lead_negI (r s : ℕ) (h : 2 ≤ r + s) :
    binHalf (s + r) * ((-Complex.I + -Complex.I) / (2 * s2 * (-Complex.I))) ^ (2 - 2 * ((s + r : ℕ) : ℤ)) *
        ((((s + r).choose s : ℕ) * 2 ^ (s + r) : ℂ) * ((1 - tv (-Complex.I)) ^ 2) ^ s * (tv (-Complex.I) ^ 2) ^ r) *
        (1 / (2 * s2 * (-Complex.I) * tv (-Complex.I) * (1 - tv (-Complex.I))))
      = -(Kd (r + s) * (-1) ^ r * (((r + s).choose r : ℕ) : ℂ)) := by
  obtain ⟨hA, hB, hP⟩ := valsnI
  rw [hA, hB, hP, Qpow _ (neg_ne_zero.2 Complex.I_ne_zero) _ (by omega)]
  have e : (Complex.I / 2) ^ s * (-Complex.I / 2) ^ r = (-1) ^ r * (Complex.I / 2) ^ (r + s) := by
    rw [show -Complex.I / 2 = (-1) * (Complex.I / 2) by ring, mul_pow, add_comm r s, pow_add]; ring
  obtain ⟨n, hn⟩ : ∃ n, r + s = n + 2 := ⟨r + s - 2, by omega⟩
  rw [mul_assoc (((s + r).choose s : ℕ) * 2 ^ (s + r) : ℂ), e]
  rw [show s + r = r + s by ring, show (r + s).choose s = (r + s).choose r by
    rw [Nat.choose_symm_add]]
  rw [hn, show n + 2 - 1 = n + 1 by omega]
  rw [← lead_core n _ _]
  ring

end P01


theorem P01.proof_1 : P01.stmt_1 := by
  intro f
  refine ⟨P01.sub_Hp_pole f _ (Or.inl rfl), P01.sub_Hp_pole f _ (Or.inr rfl),
    fun c h1 h2 => P01.Hp_regular f c h1 h2, P01.Hp_res f _ (Or.inl rfl), P01.Hp_res f _ (Or.inr rfl),
    P01.Hp_small f, P01.Jp_b f, fun z h1 h2 => P01.Jp_deriv f z h1 h2⟩

theorem P01.proof_2 : P01.stmt_2 := by
  have hcI : Complex.I ^ 2 + 1 = 0 := by rw [Complex.I_sq]; ring
  have hcm : (-Complex.I) ^ 2 + 1 = 0 := by rw [neg_sq, Complex.I_sq]; ring
  have z1 : (-2 - 2 * Complex.I : ℂ) ≠ 0 := by
    intro h; have := congrArg Complex.re h; simp at this
  have z2 : (2 - 2 * Complex.I : ℂ) ≠ 0 := by
    intro h; have := congrArg Complex.re h; simp at this
  have z3 : (-2 + 2 * Complex.I : ℂ) ≠ 0 := by
    intro h; have := congrArg Complex.re h; simp at this
  have z4 : (2 + 2 * Complex.I : ℂ) ≠ 0 := by
    intro h; have := congrArg Complex.re h; simp at this
  obtain ⟨a1, a2, a3⟩ := P01.D_facts P01.bC P01.aC Complex.I hcI (by rw [P01.num1]; exact z1)
  obtain ⟨b1, b2, b3⟩ := P01.D_facts P01.aC P01.bC Complex.I hcI (by rw [P01.num2]; exact z2)
  obtain ⟨c1, c2, c3⟩ := P01.D_facts P01.bC P01.aC (-Complex.I) hcm (by rw [P01.num3]; exact z3)
  obtain ⟨d1, d2, d3⟩ := P01.D_facts P01.aC P01.bC (-Complex.I) hcm (by rw [P01.num4]; exact z4)
  rw [P01.num1] at a3
  rw [P01.num2] at b3
  rw [P01.num3] at c3
  rw [P01.num4] at d3
  refine ⟨P01.fold_D1, a1, b1, c1, d1, a2, b2, c2, d2, ?_, ?_, ?_, ?_⟩
  · rw [P01.D1, a3, mul_inv_eq_iff_eq_mul₀ z1]; unfold P01.tau1; ring_nf; rw [Complex.I_sq]; ring
  · rw [P01.D2, b3, mul_inv_eq_iff_eq_mul₀ z2]; unfold P01.tau2; ring_nf; rw [Complex.I_sq]; ring
  · rw [P01.D1, c3, mul_inv_eq_iff_eq_mul₀ z3]; unfold P01.tau2; ring_nf; rw [Complex.I_sq]; ring
  · rw [P01.D2, d3, mul_inv_eq_iff_eq_mul₀ z4]; unfold P01.tau1; ring_nf; rw [Complex.I_sq]; ring

theorem P01.proof_5 : P01.stmt_5 := by
  intro k hk r s hrs
  have hI : (Complex.I) ^ 2 = -1 := Complex.I_sq
  have hnI : (-Complex.I) ^ 2 = -1 := by rw [neg_sq, Complex.I_sq]
  refine ⟨⟨?_, ?_⟩, fun hsr => ⟨?_, ?_, ?_⟩⟩
  · have := P01.Rt_pole _ hI r s; rwa [hrs] at this
  · have := P01.Rt_pole _ hnI r s; rwa [hrs] at this
  · have := P01.Rinit_pole _ hI k r s; rwa [hrs] at this
  · have := P01.Rinit_pole _ hnI k r s; rwa [hrs] at this
  · subst hrs
    have hd : 2 ≤ r + s := hk
    unfold P01.Rinit
    rw [if_pos ⟨rfl, hsr⟩]
    have hlead := P01.Rt_lead _ hI r s hd
    rw [P01.lead_I r s hd] at hlead
    have hidx : (-(2 * ((r + s : ℕ) : ℤ) - 2)) = -(((2 * (r + s) - 2 : ℕ)) : ℤ) := by omega
    split_ifs with hlt
    · rw [P01.lcoeff_add]
      rw [hlead]
      have hpole : P01.poleLE (-Complex.I) (P01.Rt s r) (2 * (r + s) - 2) := by
        have := P01.Rt_pole _ hnI s r; rwa [add_comm s r] at this
      have hF := P01.Fold_lead _ hI (P01.Rt s r) _ hpole
      rw [hidx, hF, ← hidx]
      have hl2 := P01.Rt_lead _ hnI s r (by omega)
      rw [P01.lead_negI r s hd] at hl2
      rw [show (-(2 * ((r + s : ℕ) : ℤ) - 2)) = (-(2 * ((s + r : ℕ) : ℤ) - 2)) by push_cast; ring, hl2]
      unfold P01.rho
      rw [if_pos hlt]
      ring
    · rw [hlead]
      unfold P01.rho
      rw [if_neg hlt]
      ring

theorem P01.proof_6 : P01.stmt_6 := by
  intro k _ r s hsr hd
  have hI := P01.W_bound k r s hsr hd _ (Or.inl rfl)
  have hJ := P01.W_bound k r s hsr hd _ (Or.inr rfl)
  refine ⟨hI, hJ, fun h2 => ?_⟩
  have e : P01.G k r s = P01.Jp (P01.W k r s) := by unfold P01.G; rw [if_pos h2]
  rw [e]
  refine ⟨?_, ?_, fun c h1 h2 => P01.poleLE_of_reg (P01.Jp_reg _ c h1 h2) 0⟩
  · exact P01.poleLE_mono (P01.Jp_pole _ (Or.inl rfl) _ _ hI) (by omega)
  · exact P01.poleLE_mono (P01.Jp_pole _ (Or.inr rfl) _ _ hJ) (by omega)


namespace P01

lemma lcoeff_mulX (c : ℂ) (g : RatFunc ℂ) (n : ℤ) :
    lcoeff c ((vX - cst c) * g) (n + 1) = lcoeff c g n := by
  have h := lcoeff_shift c 1 ((vX - cst c) * g) n
  have hX : vX - cst c ≠ 0 := X_sub_C_ne c
  rw [← mul_assoc, show (vX - cst c) ^ (-((1 : ℕ) : ℤ)) * (vX - cst c) = 1 by
    rw [Nat.cast_one, zpow_neg_one, inv_mul_cancel₀ hX], one_mul] at h
  rw [h]; push_cast; ring_nf

lemma zero_factor' {c : ℂ} {D : RatFunc ℂ} (hD : NV c D 0) :
    ∃ D', D = (vX - cst c) * D' ∧ NV c D' (lcoeff c D 1) := by
  obtain ⟨D', ⟨v, hv⟩, rfl⟩ := zero_factor hD
  refine ⟨D', rfl, ?_⟩
  have := lcoeff_mulX c D' 0
  rw [zero_add] at this
  rw [this, (NV_reg hv).2]; exact hv

lemma pole_bound_of {c : ℂ} {W : RatFunc ℂ} {M : ℕ} (hW : poleLE c W M) (m : ℕ) (hm : M < m) :
    lcoeff c W (-(m : ℤ)) = 0 := hW _ (by omega)

lemma lead_JpD (c : ℂ) (hc : c = Complex.I ∨ c = -Complex.I) (W' : RatFunc ℂ) (K : ℕ) (hK : 1 ≤ K)
    (hW : poleLE c W' (K + 2)) (D : RatFunc ℂ) (hD : NV c D 0) :
    lcoeff c (Jp W' * D) (-(K : ℤ)) = lcoeff c W' (-((K : ℤ) + 2)) * (-(1 / ((K : ℂ) + 1))) * lcoeff c D 1 := by
  obtain ⟨D', rfl, hD'⟩ := zero_factor' hD
  set τ := lcoeff c ((vX - cst c) * D') 1
  have hneg : (-(K : ℤ)) < 0 := by omega
  have hother : ∀ d : ℂ, d ≠ c → lcoeff c (jpart d W' * ((vX - cst c) * D')) (-(K : ℤ)) = 0 := by
    intro d hdc
    obtain ⟨u, hu⟩ : Reg c (jpart d W') := ⟨_, NV_jpart d W' c (Ne.symm hdc)⟩
    exact (NV_reg (NV_mul hu (NV_mul (NV_sub (by simpa [vX] using NV_Xpow c 1) (NV_cst c c)) hD'))).1 _ hneg
  have hself : lcoeff c (jpart c W' * ((vX - cst c) * D')) (-(K : ℤ))
      = lcoeff c W' (-((K : ℤ) + 2)) * (-(1 / ((K : ℂ) + 1))) * τ := by
    unfold jpart
    rw [Finset.sum_mul, lcoeff_sum]
    have term : ∀ m ∈ Finset.Icc 2 W'.denom.natDegree,
        lcoeff c (cst (lcoeff c W' (-(m : ℤ)) * (-(1 / ((m : ℂ) - 1)))) *
          ((vX - cst c) ^ (-((m : ℤ) - 1)) - cst ((bC - c) ^ (-((m : ℤ) - 1)))) * ((vX - cst c) * D'))
          (-(K : ℤ))
        = if m = K + 2 then lcoeff c W' (-((K : ℤ) + 2)) * (-(1 / ((K : ℂ) + 1))) * τ else 0 := by
      intro m hm
      have hm2 := (Finset.mem_Icc.1 hm).1
      have hX : vX - cst c ≠ 0 := X_sub_C_ne c
      have e : cst (lcoeff c W' (-(m : ℤ)) * (-(1 / ((m : ℂ) - 1)))) *
          ((vX - cst c) ^ (-((m : ℤ) - 1)) - cst ((bC - c) ^ (-((m : ℤ) - 1)))) * ((vX - cst c) * D')
          = cst (lcoeff c W' (-(m : ℤ)) * (-(1 / ((m : ℂ) - 1)))) *
            ((vX - cst c) ^ (-((m - 2 : ℕ) : ℤ)) * D'
              - cst ((bC - c) ^ (-((m : ℤ) - 1))) * ((vX - cst c) * D')) := by
        have : (vX - cst c) ^ (-((m : ℤ) - 1)) * (vX - cst c) = (vX - cst c) ^ (-((m - 2 : ℕ) : ℤ)) := by
          rw [← zpow_add_one₀ hX]; congr 1; omega
        rw [← this]; ring
      rw [e, lcoeff_cst_mul, lcoeff_sub, lcoeff_cst_mul]
      have hreg : lcoeff c ((vX - cst c) * D') (-(K : ℤ)) = 0 :=
        (NV_reg (NV_mul (NV_sub (by simpa [vX] using NV_Xpow c 1) (NV_cst c c)) hD')).1 _ hneg
      rw [hreg, mul_zero, sub_zero, lcoeff_red c _ D' _ hD']
      split_ifs with h1 h2 h2
      · omega
      · rw [mul_zero]
      · subst h2
        rw [show (-(K : ℤ)) + ((K + 2 - 2 : ℕ) : ℤ) = 0 by omega, (NV_reg hD').2]
        push_cast
        ring_nf
      · by_cases hmK : m ≤ K + 2
        · omega
        · rw [pole_bound_of hW m (by omega), zero_mul, zero_mul]
    rw [Finset.sum_congr rfl term, Finset.sum_ite_eq']
    split_ifs with hmem
    · rfl
    · rw [Finset.mem_Icc] at hmem
      have : W'.denom.natDegree < K + 2 := by omega
      rw [pole_bound c W' _ (by omega), zero_mul, zero_mul]
  unfold Jp
  rw [add_mul, lcoeff_add]
  rcases hc with rfl | rfl
  · rw [hself, hother _ (Ne.symm I_ne_negI), add_zero]
  · rw [hself, hother _ I_ne_negI, zero_add]

end P01

namespace P01

/-- `h` at `c`. -/
noncomputable def hg (c : ℂ) (k r s : ℕ) : ℂ := lcoeff c (W k r s) (-(2 * ((r + s : ℕ) : ℤ) - 2))

lemma W_rec (c : ℂ) (hc : c = Complex.I ∨ c = -Complex.I) (k r s : ℕ) (hsr : s ≤ r) (hd : 2 ≤ r + s)
    (hk : r + s ≤ k) :
    (2 * ((r + s : ℕ) : ℂ) - 1) * hg c k r s
      = (2 * ((r + s : ℕ) : ℂ) - 1) * lcoeff c (Rinit k r s) (-(2 * ((r + s : ℕ) : ℤ) - 2))
        + ((r : ℂ) + 1) * lcoeff c D1 1 * hg c k (r + 1) s
        + (if s < r then ((s : ℂ) + 1) * lcoeff c D2 1 * hg c k r (s + 1) else 0)
        - (if s + 1 = r then (r : ℂ) * lcoeff (-c) D1 1 * hg (-c) k r r else 0) := by
  set K := 2 * (r + s) - 2 with hKdef
  have hK1 : 1 ≤ K := by omega
  have hidx : (-(2 * ((r + s : ℕ) : ℤ) - 2)) = -(K : ℤ) := by omega
  have hKc : ((K : ℂ) + 1) = 2 * ((r + s : ℕ) : ℂ) - 1 := by
    rw [hKdef, Nat.cast_sub (by omega)]; push_cast; ring
  have hKne : ((K : ℂ) + 1) ≠ 0 := by
    rw [hKc]; intro h
    have : (2 * ((r + s : ℕ) : ℝ) - 1) = 0 := by exact_mod_cast h
    have : (2 : ℝ) ≤ ((r + s : ℕ) : ℝ) := by exact_mod_cast hd
    linarith
  have hnc := neg_mem c hc
  have bnd : ∀ c' : ℂ, (c' = Complex.I ∨ c' = -Complex.I) → ∀ r' s', s' ≤ r' → r' + s' = r + s + 1 →
      poleLE c' (W k r' s') (K + 2) := by
    intro c' hc' r' s' h1 h2
    have := W_bound k r' s' h1 (by omega) c' hc'
    exact poleLE_mono this (by omega)
  have idx2 : ∀ r' s', r' + s' = r + s + 1 →
      (-(2 * ((r' + s' : ℕ) : ℤ) - 2)) = -((K : ℤ) + 2) := by
    intro r' s' h; rw [h]; omega
  unfold hg
  rw [hidx, W_eq k r s hk]
  rw [lcoeff_sub, lcoeff_sub, lcoeff_sub, mul_assoc, lcoeff_cst_mul,
    lead_JpD c hc _ K hK1 (bnd c hc (r + 1) s (by omega) (by omega)) D1 (D1_zero c hc)]
  rw [idx2 (r + 1) s (by omega)]
  have t3 : lcoeff c (if s < r then cst ((s : ℂ) + 1) * Jp (W k r (s + 1)) * D2 else 0) (-(K : ℤ))
      = if s < r then ((s : ℂ) + 1) * (lcoeff c (W k r (s + 1)) (-((K : ℤ) + 2)) * (-(1 / ((K : ℂ) + 1)))
          * lcoeff c D2 1) else 0 := by
    split_ifs with h
    · rw [mul_assoc, lcoeff_cst_mul,
        lead_JpD c hc _ K hK1 (bnd c hc r (s + 1) (by omega) (by omega)) D2 (D2_zero c hc)]
    · exact lcoeff_zero c _
  have t4 : lcoeff c (if _h2 : s + 1 = r then cst (r : ℂ) * Fold (Jp (W k r r) * D1) else 0) (-(K : ℤ))
      = if s + 1 = r then (r : ℂ) * (-(lcoeff (-c) (W k r r) (-((K : ℤ) + 2)) * (-(1 / ((K : ℂ) + 1)))
          * lcoeff (-c) D1 1)) else 0 := by
    split_ifs with h
    · rw [lcoeff_cst_mul]
      have hp : poleLE (-c) (Jp (W k r r) * D1) K := by
        have := poleLE_mul_zero (Jp_pole (-c) hnc _ _ (bnd (-c) hnc r r le_rfl (by omega))) (D1_zero (-c) hnc)
        exact poleLE_mono this (by omega)
      rw [Fold_lead c (sqI c hc) _ K hp,
        lead_JpD (-c) hnc _ K hK1 (bnd (-c) hnc r r le_rfl (by omega)) D1 (D1_zero (-c) hnc)]
    · exact lcoeff_zero c _
  rw [t3, t4]
  have e3 : (-(2 * ((r + (s + 1) : ℕ) : ℤ) - 2)) = -((K : ℤ) + 2) := idx2 r (s + 1) (by omega)
  rw [e3, ← hKc]
  split_ifs with h1 h2
  · have e4 : (-(2 * ((r + r : ℕ) : ℤ) - 2)) = -((K : ℤ) + 2) := idx2 r r (by omega)
    rw [e4]; field_simp; ring
  · field_simp; ring
  · omega
  · field_simp; ring

end P01

namespace P01

lemma conj_s2 : (starRingEnd ℂ) s2 = s2 := by unfold s2; exact Complex.conj_ofReal _

lemma conj_Kd (d : ℕ) : (starRingEnd ℂ) (Kd d) = (-1) ^ (d + 1) * Kd d := by
  unfold Kd
  simp only [map_div₀, map_mul, map_pow, map_neg, map_one, Complex.conj_I, map_natCast, conj_s2,
    map_ofNat]
  rw [neg_pow Complex.I]; ring

lemma conj_rho (r s : ℕ) : (starRingEnd ℂ) (rho r s) =
    (-1) ^ (r + s + 1) * Kd (r + s) * (-1) ^ r * (if s < r then 2 else 1) * (((r + s).choose r : ℕ) : ℂ) := by
  unfold rho
  simp only [map_mul, map_pow, map_neg, map_one, map_natCast, conj_Kd]
  split_ifs
  · rw [map_ofNat]
  · rw [map_one]

lemma lead_negI_self (r s : ℕ) (h : 2 ≤ r + s) :
    binHalf (r + s) * ((-Complex.I + -Complex.I) / (2 * s2 * (-Complex.I))) ^ (2 - 2 * ((r + s : ℕ) : ℤ)) *
        ((((r + s).choose r : ℕ) * 2 ^ (r + s) : ℂ) * ((1 - tv (-Complex.I)) ^ 2) ^ r * (tv (-Complex.I) ^ 2) ^ s) *
        (1 / (2 * s2 * (-Complex.I) * tv (-Complex.I) * (1 - tv (-Complex.I))))
      = Kd (r + s) * (-(-1) ^ s) * (((r + s).choose r : ℕ) : ℂ) := by
  obtain ⟨hA, hB, hP⟩ := valsnI
  rw [hA, hB, hP, Qpow _ (neg_ne_zero.2 Complex.I_ne_zero) _ (by omega)]
  have e : (Complex.I / 2) ^ r * (-Complex.I / 2) ^ s = (-1) ^ s * (Complex.I / 2) ^ (r + s) := by
    rw [show -Complex.I / 2 = (-1) * (Complex.I / 2) by ring, mul_pow, pow_add]; ring
  obtain ⟨n, hn⟩ : ∃ n, r + s = n + 2 := ⟨r + s - 2, by omega⟩
  rw [mul_assoc (((r + s).choose r : ℕ) * 2 ^ (r + s) : ℂ), e]
  rw [hn, show n + 2 - 1 = n + 1 by omega]
  rw [← lead_core n _ (-(-1) ^ s)]
  ring
lemma proof_5' : stmt_5 := P01.proof_5

lemma Rinit_lead_I (k r s : ℕ) (hsr : s ≤ r) (hd : 2 ≤ r + s) :
    lcoeff Complex.I (Rinit k r s) (-(2 * ((r + s : ℕ) : ℤ) - 2)) = if k = r + s then rho r s else 0 := by
  split_ifs with hk
  · subst hk
    have := (proof_5' (r + s) hd r s rfl).2 hsr
    exact this.2.2
  · unfold Rinit; rw [if_neg (by omega)]; exact lcoeff_zero _ _

end P01

namespace P01

lemma neg_one_aux (r s : ℕ) : (-1 : ℂ) ^ (r + s + 1) * (-1) ^ r = -(-1) ^ s := by
  rw [← pow_add, show r + s + 1 + r = 2 * r + (s + 1) by ring, pow_add, pow_mul]; norm_num; ring

lemma Rinit_lead_negI (k r s : ℕ) (hsr : s ≤ r) (hd : 2 ≤ r + s) :
    lcoeff (-Complex.I) (Rinit k r s) (-(2 * ((r + s : ℕ) : ℤ) - 2))
      = if k = r + s then (starRingEnd ℂ) (rho r s) else 0 := by
  have hnI : (-Complex.I) ^ 2 = -1 := by rw [neg_sq, Complex.I_sq]
  split_ifs with hk
  · subst hk
    rw [conj_rho]
    unfold Rinit
    rw [if_pos ⟨rfl, hsr⟩]
    have hl := Rt_lead _ hnI r s hd
    rw [lead_negI_self r s hd] at hl
    split_ifs with hlt
    · rw [lcoeff_add, hl]
      have hpole : poleLE (-(-Complex.I)) (Rt s r) (2 * (r + s) - 2) := by
        rw [neg_neg]; have := Rt_pole _ Complex.I_sq s r; rwa [add_comm s r] at this
      have hidx : (-(2 * ((r + s : ℕ) : ℤ) - 2)) = -(((2 * (r + s) - 2 : ℕ)) : ℤ) := by omega
      rw [hidx, Fold_lead _ hnI _ _ hpole, ← hidx, neg_neg]
      have hl2 := Rt_lead _ Complex.I_sq s r (by omega)
      rw [lead_I s r (by omega)] at hl2
      rw [show (-(2 * ((r + s : ℕ) : ℤ) - 2)) = (-(2 * ((s + r : ℕ) : ℤ) - 2)) by push_cast; ring, hl2]
      rw [show s + r = r + s by ring, show (r + s).choose s = (r + s).choose r by rw [Nat.choose_symm_add]]
      have := neg_one_aux r s
      linear_combination (-(2 : ℂ) * Kd (r + s) * (((r + s).choose r : ℕ) : ℂ)) * this
    · rw [hl]
      have := neg_one_aux r s
      have hrs : s = r := by omega
      subst hrs
      linear_combination (-(1 : ℂ) * Kd (s + s) * (((s + s).choose s : ℕ) : ℂ)) * this
  · unfold Rinit; rw [if_neg (by omega)]; exact lcoeff_zero _ _

lemma proof_2' : stmt_2 := P01.proof_2

lemma D_vals : lcoeff Complex.I D1 1 = tau1 ∧ lcoeff Complex.I D2 1 = tau2 ∧
    lcoeff (-Complex.I) D1 1 = tau2 ∧ lcoeff (-Complex.I) D2 1 = tau1 := by
  have := proof_2'
  exact ⟨this.2.2.2.2.2.2.2.2.2.1, this.2.2.2.2.2.2.2.2.2.2.1, this.2.2.2.2.2.2.2.2.2.2.2.1,
    this.2.2.2.2.2.2.2.2.2.2.2.2⟩

lemma conj_tau1 : (starRingEnd ℂ) tau1 = tau2 := by
  unfold tau1 tau2; rw [map_div₀, map_sub, map_neg, map_one, Complex.conj_I, map_ofNat]; ring
lemma conj_tau2 : (starRingEnd ℂ) tau2 = tau1 := by
  unfold tau1 tau2; rw [map_div₀, map_add, map_neg, map_one, Complex.conj_I, map_ofNat]; ring

lemma hg_conj (k : ℕ) : ∀ r s : ℕ, s ≤ r → 2 ≤ r + s →
    hg (-Complex.I) k r s = (starRingEnd ℂ) (hg Complex.I k r s) := by
  obtain ⟨d1I, d2I, d1n, d2n⟩ := D_vals
  suffices H : ∀ n : ℕ, ∀ r s : ℕ, k + 1 - (r + s) = n → s ≤ r → 2 ≤ r + s →
      hg (-Complex.I) k r s = (starRingEnd ℂ) (hg Complex.I k r s) by
    intro r s h1 h2; exact H _ r s rfl h1 h2
  intro n
  induction n using Nat.strong_induction_on with
  | _ n ih =>
  intro r s hn hsr hd
  by_cases hk : r + s ≤ k
  swap
  · unfold hg; rw [W_gt k r s (by omega), lcoeff_zero, lcoeff_zero, map_zero]
  have IH : ∀ r' s', s' ≤ r' → r' + s' = r + s + 1 →
      hg (-Complex.I) k r' s' = (starRingEnd ℂ) (hg Complex.I k r' s') := fun r' s' h1 h2 =>
    ih (k + 1 - (r' + s')) (by omega) r' s' rfl h1 (by omega)
  have e1 := W_rec _ (Or.inl rfl) k r s hsr hd hk
  have e2 := W_rec _ (Or.inr rfl) k r s hsr hd hk
  rw [neg_neg] at e2
  rw [Rinit_lead_I k r s hsr hd, d1I, d2I, d1n] at e1
  rw [Rinit_lead_negI k r s hsr hd, d1n, d2n, d1I] at e2
  have hne : (2 * ((r + s : ℕ) : ℂ) - 1) ≠ 0 := by
    intro h
    have : (2 * ((r + s : ℕ) : ℝ) - 1) = 0 := by exact_mod_cast h
    have : (2 : ℝ) ≤ ((r + s : ℕ) : ℝ) := by exact_mod_cast hd
    linarith
  have IHrr : hg Complex.I k r r = (starRingEnd ℂ) (hg (-Complex.I) k r r) → True := fun _ => trivial
  clear IHrr
  split_ifs at e1 e2 with h1 h2 h3 <;>
  · apply mul_left_cancel₀ hne
    rw [e2]
    have c1 := congrArg (starRingEnd ℂ) e1
    simp only [map_mul, map_add, map_sub, map_natCast, map_one, map_ofNat, conj_tau1, conj_tau2,
      map_zero] at c1
    try rw [IH r (s + 1) (by omega) (by omega)]
    try rw [IH r r le_rfl (by omega), Complex.conj_conj] at c1
    rw [IH (r + 1) s (by omega) (by omega)]
    linear_combination (-1 : ℂ) * c1

end P01

theorem P01.proof_7 : P01.stmt_7 := by
  intro k r s hsr hd hk
  obtain ⟨d1I, d2I, d1n, d2n⟩ := P01.D_vals
  have e1 := P01.W_rec _ (Or.inl rfl) k r s hsr hd hk
  rw [P01.Rinit_lead_I k r s hsr hd, d1I, d2I, d1n] at e1
  have hh : ∀ a b, P01.h k a b = P01.hg Complex.I k a b := fun _ _ => rfl
  simp only [hh]
  rcases (show s + 2 ≤ r ∨ s + 1 = r ∨ s = r by omega) with h1 | h1 | h1
  · rw [if_pos (by omega : s < r), if_neg (by omega : ¬ s + 1 = r)] at e1
    rw [if_pos h1, if_neg (by omega : ¬ s + 1 = r)]
    linear_combination e1
  · obtain rfl : r = s + 1 := h1.symm
    rw [if_pos (by omega : s < s + 1), if_pos rfl] at e1
    rw [if_neg (by omega : ¬ s + 2 ≤ s + 1), if_pos rfl]
    rw [P01.hg_conj k (s + 1) (s + 1) le_rfl (by omega)] at e1
    push_cast at e1 ⊢
    linear_combination e1
  · subst h1
    rw [if_neg (lt_irrefl s), if_neg (by omega : ¬ s + 1 = s)] at e1
    rw [if_neg (by omega : ¬ s + 2 ≤ s), if_neg (by omega : ¬ s + 1 = s)]
    linear_combination e1

theorem P01.proof_8 : P01.stmt_8 := by
  refine ⟨P01.Kd_succ, ?_, ?_, ?_⟩
  · intro r s hd hrs
    have hK := P01.Kd_succ (r + s) hd
    unfold P01.rho
    have i1 : r + 1 + s = r + s + 1 := by ring
    have i2 : r + (s + 1) = r + s + 1 := by ring
    rw [i1, i2]
    have c1 : ((r : ℂ) + 1) * (((r + s + 1).choose (r + 1) : ℕ) : ℂ)
        = ((r : ℂ) + s + 1) * (((r + s).choose r : ℕ) : ℂ) := by
      have h' := Nat.add_one_mul_choose_eq (r + s) r
      have h'' : (((r + s + 1) * (r + s).choose r : ℕ) : ℂ) = (((r + s + 1).choose (r + 1) * (r + 1) : ℕ) : ℂ) := by rw [h']
      push_cast at h''
      linear_combination -h''
    have c2 : ((s : ℂ) + 1) * (((r + s + 1).choose r : ℕ) : ℂ)
        = ((r : ℂ) + s + 1) * (((r + s).choose r : ℕ) : ℂ) := by
      have h' := Nat.choose_mul_succ_eq (r + s) r
      have h'' : (((r + s).choose r * (r + s + 1) : ℕ) : ℂ) = (((r + s + 1).choose r * (r + s + 1 - r) : ℕ) : ℂ) := by rw [h']
      rw [show r + s + 1 - r = s + 1 by omega] at h''
      push_cast at h''
      linear_combination -h''
    have l1 : s < r + 1 := by omega
    have l2 : s + 1 < r := by omega
    have l3 : s < r := by omega
    simp only [l1, l2, l3, ite_true]
    push_cast at hK ⊢
    unfold P01.tau1 P01.tau2
    rw [pow_succ]
    linear_combination (-(-1:ℂ)^r * P01.Kd (r+s+1) * (1 + Complex.I)) * c1
      + (-(-1:ℂ)^r * P01.Kd (r+s+1) * (-1 + Complex.I)) * c2
      + (-2 * Complex.I * (-1:ℂ)^r * (((r+s).choose r : ℕ) : ℂ)) * hK
      + (2 * (-1:ℂ)^r * (((r+s).choose r : ℕ) : ℂ) * (2 * ((r:ℂ)+s) - 1) * P01.Kd (r+s)) * Complex.I_mul_I
  · intro p hp
    obtain ⟨n, rfl⟩ : ∃ n, p = n + 1 := ⟨p - 1, by omega⟩
    have hK := P01.Kd_succ (2 * n + 1) (by omega)
    rw [show n + 1 - 1 = n by omega, P01.conj_rho_pp]
    unfold P01.rho
    rw [show n + 1 + n = 2 * n + 1 by ring, show n + 1 + 1 + n = 2 * n + 1 + 1 by ring,
      show n + 1 + (n + 1) = 2 * n + 1 + 1 by ring]
    have c1 : ((n : ℂ) + 2) * (((2 * n + 1 + 1).choose (n + 1 + 1) : ℕ) : ℂ)
        = (2 * (n : ℂ) + 2) * (((2 * n + 1).choose (n + 1) : ℕ) : ℂ) := by
      have h' := Nat.add_one_mul_choose_eq (2 * n + 1) (n + 1)
      have h'' : (((2 * n + 1 + 1) * (2 * n + 1).choose (n + 1) : ℕ) : ℂ)
          = (((2 * n + 1 + 1).choose (n + 1 + 1) * (n + 1 + 1) : ℕ) : ℂ) := by rw [h']
      push_cast at h''
      linear_combination -h''
    have c2 : ((n : ℂ) + 1) * (((2 * n + 1 + 1).choose (n + 1) : ℕ) : ℂ)
        = (2 * (n : ℂ) + 2) * (((2 * n + 1).choose (n + 1) : ℕ) : ℂ) := by
      have h' := Nat.add_one_mul_choose_eq (2 * n + 1) n
      rw [← Nat.choose_symm_half n] at h'
      have h'' : (((2 * n + 1 + 1) * (2 * n + 1).choose (n + 1) : ℕ) : ℂ)
          = (((2 * n + 1 + 1).choose (n + 1) * (n + 1) : ℕ) : ℂ) := by rw [h']
      push_cast at h''
      linear_combination -h''
    have l1 : n < n + 1 := by omega
    have l2 : n < n + 1 + 1 := by omega
    simp only [l1, l2, lt_irrefl, ite_true, ite_false]
    push_cast at hK ⊢
    unfold P01.tau1 P01.tau2
    rw [pow_succ (-1 : ℂ) (n + 1)]
    linear_combination (-(-1:ℂ)^(n+1) * P01.Kd (2*n+1+1) * (1 + Complex.I)) * c1
      + (-(-1:ℂ)^(n+1) * P01.Kd (2*n+1+1) * (-1 + Complex.I)) * c2
      + (-2 * Complex.I * (-1:ℂ)^(n+1) * (((2*n+1).choose (n+1) : ℕ) : ℂ)) * hK
      + (2 * (-1:ℂ)^(n+1) * (((2*n+1).choose (n+1) : ℕ) : ℂ) * (4 * (n:ℂ) + 1)
          * P01.Kd (2*n+1)) * Complex.I_mul_I
  · intro p hp
    have e1 : P01.hhat p = (1 - Complex.I) * P01.rho p p := by
      have hK := P01.Kd_succ (p + p) (by omega)
      have h4 : (4 * (p : ℂ) - 1) ≠ 0 := by
        intro h
        have : (4 * (p : ℝ) - 1) = 0 := by exact_mod_cast h
        have : (1 : ℝ) ≤ p := by exact_mod_cast hp
        linarith
      unfold P01.hhat
      rw [div_eq_iff h4]
      unfold P01.rho
      rw [show p + 1 + p = p + p + 1 by ring]
      have c : ((p : ℂ) + 1) * (((p + p + 1).choose (p + 1) : ℕ) : ℂ)
          = ((p : ℂ) + p + 1) * (((p + p).choose p : ℕ) : ℂ) := by
        have h' := Nat.add_one_mul_choose_eq (p + p) p
        have h'' : (((p + p + 1) * (p + p).choose p : ℕ) : ℂ)
            = (((p + p + 1).choose (p + 1) * (p + 1) : ℕ) : ℂ) := by rw [h']
        push_cast at h''
        linear_combination -h''
      have l1 : p < p + 1 := by omega
      simp only [l1, lt_irrefl, ite_true, ite_false]
      push_cast at hK ⊢
      unfold P01.tau1
      rw [pow_succ (-1 : ℂ) p]
      linear_combination (-2 * (-1:ℂ)^p * ((-1 - Complex.I) / 2) * P01.Kd (p+p+1)) * c
        + (-2 * (-1:ℂ)^p * ((-1 - Complex.I) / 2) * (((p + p).choose p : ℕ) : ℂ)) * hK
        + (-(-1:ℂ)^p * (((p + p).choose p : ℕ) : ℂ) * (4 * (p:ℂ) - 1) * P01.Kd (p+p))
          * Complex.I_mul_I
    refine ⟨e1, ?_, ?_⟩
    · obtain ⟨x, hx⟩ := P01.rho_pp p
      rw [hx]; simp
    · rw [e1, map_mul, P01.conj_rho_pp]
      simp only [map_sub, map_one, Complex.conj_I]
      ring

namespace P01

lemma h_gt (k r s : ℕ) (h : k < r + s) : P01.h k r s = 0 := by
  unfold P01.h; rw [W_gt k r s h]; exact lcoeff_zero _ _

lemma two_d_ne (d : ℕ) (hd : 2 ≤ d) : (2 * (d : ℂ) - 1) ≠ 0 := by
  intro h
  have : (2 * (d : ℝ) - 1) = 0 := by exact_mod_cast h
  have : (2 : ℝ) ≤ d := by exact_mod_cast hd
  linarith

lemma proof_7' : stmt_7 := P01.proof_7
lemma proof_8' : stmt_8 := P01.proof_8

lemma lemD (k : ℕ) : ∀ r s : ℕ, s ≤ r → 2 ≤ r + s → r + s ≤ k →
    (s < r → P01.h k r s = rho r s) ∧ (s = r → (k = 2 * r → P01.h k r r = rho r r) ∧
      (2 * r + 1 ≤ k → P01.h k r r = hhat r)) := by
  obtain ⟨I1, I2, I3, I4⟩ := proof_8'
  suffices H : ∀ n : ℕ, ∀ r s : ℕ, k + 1 - (r + s) = n → s ≤ r → 2 ≤ r + s → r + s ≤ k →
      (s < r → P01.h k r s = rho r s) ∧ (s = r → (k = 2 * r → P01.h k r r = rho r r) ∧
        (2 * r + 1 ≤ k → P01.h k r r = hhat r)) by
    intro r s h1 h2 h3; exact H _ r s rfl h1 h2 h3
  intro n
  induction n using Nat.strong_induction_on with
  | _ n ih =>
  intro r s hn hsr hd hk
  have IH : ∀ r' s', s' ≤ r' → r' + s' = r + s + 1 → r' + s' ≤ k →
      (s' < r' → P01.h k r' s' = rho r' s') ∧ (s' = r' → (k = 2 * r' → P01.h k r' r' = rho r' r') ∧
        (2 * r' + 1 ≤ k → P01.h k r' r' = hhat r')) := fun r' s' h1 h2 h3 =>
    ih (k + 1 - (r' + s')) (by omega) r' s' rfl h1 (by omega) h3
  have e := proof_7' k r s hsr hd hk
  have hne := two_d_ne (r + s) hd
  rcases eq_or_lt_of_le hk with hkd | hkd
  · -- top degree: all parents vanish
    rw [if_pos hkd.symm, h_gt k (r + 1) s (by omega)] at e
    have hpar : (if s + 2 ≤ r then ((s : ℂ) + 1) * tau2 * P01.h k r (s + 1) else 0) = 0 := by
      split_ifs
      · rw [h_gt k r (s + 1) (by omega), mul_zero]
      · rfl
    have hpar2 : (if s + 1 = r then (r : ℂ) * tau2 * (P01.h k r r - (starRingEnd ℂ) (P01.h k r r)) else 0) = 0 := by
      split_ifs
      · rw [h_gt k r r (by omega), map_zero, sub_zero, mul_zero]
      · rfl
    rw [hpar, hpar2, mul_zero, zero_add, add_zero] at e
    have hh : P01.h k r s = rho r s := by
      have := mul_eq_zero.1 e
      rcases this with h0 | h0
      · exact absurd h0 hne
      · exact sub_eq_zero.1 h0
    refine ⟨fun _ => hh, fun hrs => ⟨fun _ => by subst hrs; exact hh, fun h2 => by omega⟩⟩
  · rw [if_neg (by omega)] at e
    rw [sub_zero] at e
    refine ⟨fun hlt => ?_, fun hrs => ⟨fun h2 => by omega, fun _ => ?_⟩⟩
    · have h1 := (IH (r + 1) s (by omega) (by omega) (by omega)).1 (by omega)
      rw [h1] at e
      rcases (show s + 2 ≤ r ∨ s + 1 = r by omega) with h2 | h2
      · rw [if_pos h2, if_neg (by omega)] at e
        have h3 := (IH r (s + 1) (by omega) (by omega) (by omega)).1 (by omega)
        rw [h3, add_zero, ← I2 r s hd h2] at e
        exact mul_left_cancel₀ hne e
      · rw [if_neg (by omega), if_pos h2, add_zero] at e
        obtain rfl : r = s + 1 := h2.symm
        have hdiag := (IH (s + 1) (s + 1) le_rfl (by omega) (by omega)).2 rfl
        have hc : P01.h k (s + 1) (s + 1) - (starRingEnd ℂ) (P01.h k (s + 1) (s + 1))
            = rho (s + 1) (s + 1) - (starRingEnd ℂ) (rho (s + 1) (s + 1)) := by
          rcases eq_or_lt_of_le (show 2 * (s + 1) ≤ k by omega) with h4 | h4
          · rw [hdiag.1 h4.symm]
          · rw [hdiag.2 (by omega)]; exact (I4 (s + 1) (by omega)).2.2
        rw [hc] at e
        have i3 := I3 (s + 1) (by omega)
        rw [show s + 1 - 1 = s by omega] at i3
        rw [← i3] at e
        apply mul_left_cancel₀ hne
        rw [e]; push_cast; ring
    · subst hrs
      have h1 := (IH (s + 1) s (by omega) (by omega) (by omega)).1 (by omega)
      rw [h1, if_neg (by omega), if_neg (by omega), add_zero, add_zero] at e
      unfold hhat
      have h4 : (4 * (s : ℂ) - 1) ≠ 0 := by
        have := two_d_ne (s + s) hd; push_cast at this; convert this using 1; ring
      rw [eq_div_iff h4]
      push_cast at e
      linear_combination e

end P01

theorem P01.proof_9 : P01.stmt_9 := fun k r s h1 h2 h3 => P01.lemD k r s h1 h2 h3


namespace P01

lemma jpart_two (c : ℂ) (W : RatFunc ℂ) (hW : poleLE c W 2) :
    jpart c W = cst (-(lcoeff c W (-2))) * ((vX - cst c) ^ (-1 : ℤ) - cst ((bC - c) ^ (-1 : ℤ))) := by
  unfold jpart
  by_cases hN : 2 ≤ W.denom.natDegree
  · rw [Finset.sum_eq_single 2]
    · norm_num
    · intro m hm hm2
      have h3 : 2 < m := by have := (Finset.mem_Icc.1 hm).1; omega
      rw [hW (-(m : ℤ)) (by omega), zero_mul]; unfold cst; rw [map_zero, zero_mul]
    · intro h; exact absurd (Finset.mem_Icc.2 ⟨le_rfl, hN⟩) h
  · rw [Finset.Icc_eq_empty (by omega), Finset.sum_empty]
    rw [pole_bound c W (-2) (by omega), neg_zero]; unfold cst; rw [map_zero, zero_mul]

lemma Jp_two (W : RatFunc ℂ) (hI : poleLE Complex.I W 2) (hJ : poleLE (-Complex.I) W 2) :
    Jp W = cst (-(lcoeff Complex.I W (-2))) * (vX - cst Complex.I)⁻¹
      + cst (-(lcoeff (-Complex.I) W (-2))) * (vX + cst Complex.I)⁻¹
      + cst (lcoeff Complex.I W (-2) * (bC - Complex.I)⁻¹ + lcoeff (-Complex.I) W (-2) * (bC + Complex.I)⁻¹) := by
  unfold Jp
  rw [jpart_two _ W hI, jpart_two _ W hJ]
  have e : vX - cst (-Complex.I) = vX + cst Complex.I := by unfold cst; rw [map_neg, sub_neg_eq_add]
  rw [e, sub_neg_eq_add]
  simp only [zpow_neg_one]
  unfold cst
  simp only [map_add, map_mul, map_neg]
  ring

lemma pf_X2 : (vX ^ 2 + 1)⁻¹ = cst (-Complex.I / 2) * (vX - cst Complex.I)⁻¹ + cst (Complex.I / 2) * (vX + cst Complex.I)⁻¹ := by
  have h1 : vX - cst Complex.I ≠ 0 := X_sub_C_ne _
  have h2 : vX + cst Complex.I ≠ 0 := by
    have := X_sub_C_ne (-Complex.I); unfold cst at *; rwa [map_neg, sub_neg_eq_add] at this
  have hI2 : cst Complex.I * cst Complex.I = -1 := by
    unfold cst; rw [← map_mul, Complex.I_mul_I, map_neg, map_one]
  have hX2 : vX ^ 2 + 1 = (vX - cst Complex.I) * (vX + cst Complex.I) := by
    linear_combination hI2
  have hc2 : cst (-Complex.I / 2) = -cst Complex.I / 2 := by
    unfold cst; rw [map_div₀, map_neg, map_ofNat]
  have hc3 : cst (Complex.I / 2) = cst Complex.I / 2 := by
    unfold cst; rw [map_div₀, map_ofNat]
  rw [hX2, hc2, hc3]
  field_simp
  linear_combination (2 : RatFunc ℂ) * hI2

lemma rho20 : rho 2 0 = -Complex.I / (2 * s2) := by
  unfold rho Kd
  have := s2_ne
  norm_num [Nat.doubleFactorial]
  field_simp
  try ring

lemma rho11 : rho 1 1 = Complex.I / (2 * s2) := by
  unfold rho Kd
  have := s2_ne
  norm_num [Nat.doubleFactorial]
  field_simp
  try ring

lemma hhat1 : hhat 1 = (1 + Complex.I) / (2 * s2) := by
  rw [(proof_8'.2.2.2 1 le_rfl).1, rho11]
  have := s2_ne
  field_simp; ring_nf; rw [Complex.I_sq]; ring

end P01
namespace P01

lemma G_eq_form (k r s : ℕ) (hsr : s ≤ r) (hd : r + s = 2) (a : ℂ) (ha : P01.h k r s = a) :
    G k r s = cst (-a) * (vX - cst Complex.I)⁻¹ + cst (-(starRingEnd ℂ) a) * (vX + cst Complex.I)⁻¹
      + cst (a * (bC - Complex.I)⁻¹ + (starRingEnd ℂ) a * (bC + Complex.I)⁻¹) := by
  unfold G
  rw [if_pos (by omega)]
  have hI := W_bound k r s hsr (by omega) _ (Or.inl rfl)
  have hJ := W_bound k r s hsr (by omega) _ (Or.inr rfl)
  rw [hd] at hI hJ
  rw [Jp_two _ hI hJ]
  have e1 : lcoeff Complex.I (W k r s) (-2) = a := by
    rw [← ha]; unfold P01.h; congr 1; rw [hd]; norm_num
  have e2 : lcoeff (-Complex.I) (W k r s) (-2) = (starRingEnd ℂ) a := by
    have := hg_conj k r s hsr (by omega)
    unfold hg at this
    rw [hd] at this
    norm_num at this
    rw [this, ← e1]
  rw [e1, e2]

lemma g20_form : g20 = cst (s2 * Complex.I / 4) * (vX - cst Complex.I)⁻¹
    + cst (-(s2 * Complex.I / 4)) * (vX + cst Complex.I)⁻¹ + cst (aC / 4) := by
  unfold g20
  rw [div_eq_mul_inv (cst (s2 / 2)), pf_X2]
  unfold cst
  simp only [map_mul, map_div₀, map_neg, map_ofNat]
  ring

lemma pf2 : (2 * vX - 2) * (vX ^ 2 + 1)⁻¹
    = cst (1 + Complex.I) * (vX - cst Complex.I)⁻¹ + cst (1 - Complex.I) * (vX + cst Complex.I)⁻¹ := by
  have h1 : vX - cst Complex.I ≠ 0 := X_sub_C_ne _
  have h2 : vX + cst Complex.I ≠ 0 := by
    have := X_sub_C_ne (-Complex.I); unfold cst at *; rwa [map_neg, sub_neg_eq_add] at this
  have hI2 : cst Complex.I * cst Complex.I = -1 := by
    unfold cst; rw [← map_mul, Complex.I_mul_I, map_neg, map_one]
  have hX2 : vX ^ 2 + 1 = (vX - cst Complex.I) * (vX + cst Complex.I) := by
    linear_combination hI2
  have hc2 : cst (1 + Complex.I) = 1 + cst Complex.I := by unfold cst; rw [map_add, map_one]
  have hc3 : cst (1 - Complex.I) = 1 - cst Complex.I := by unfold cst; rw [map_sub, map_one]
  rw [hX2, hc2, hc3]
  field_simp
  linear_combination (-2 : RatFunc ℂ) * hI2

lemma g11_form : g11 = cst (-(s2 * (1 + Complex.I) / 4)) * (vX - cst Complex.I)⁻¹
    + cst (-(s2 * (1 - Complex.I) / 4)) * (vX + cst Complex.I)⁻¹ + cst (-(s2 / 4)) := by
  have hX : vX ^ 2 + 1 ≠ 0 := by
    have h1 : vX - cst Complex.I ≠ 0 := X_sub_C_ne _
    have h2 : vX + cst Complex.I ≠ 0 := by
      have := X_sub_C_ne (-Complex.I); unfold cst at *; rwa [map_neg, sub_neg_eq_add] at this
    have hI2 : cst Complex.I * cst Complex.I = -1 := by
      unfold cst; rw [← map_mul, Complex.I_mul_I, map_neg, map_one]
    have hX2 : vX ^ 2 + 1 = (vX - cst Complex.I) * (vX + cst Complex.I) := by
      linear_combination hI2
    rw [hX2]; exact mul_ne_zero h1 h2
  have e : g11 = cst (-(s2 / 4)) * (1 + (2 * vX - 2) * (vX ^ 2 + 1)⁻¹) := by
    unfold g11
    unfold cst; rw [map_neg, map_div₀, map_ofNat]
    field_simp
    ring
  rw [e, pf2]
  unfold cst
  simp only [map_mul, map_div₀, map_neg, map_ofNat, map_add, map_sub, map_one]
  ring

lemma bI_ne : bC - Complex.I ≠ 0 := sub_ne_zero.2 b_ne_I
lemma bnI_ne : bC + Complex.I ≠ 0 := by
  intro h; exact b_ne_negI (by linear_combination h)

lemma s2_pow3 : s2 ^ 3 = 2 * s2 := by rw [pow_succ, s2_sq]
lemma s2_pow4 : s2 ^ 4 = 4 := by rw [show (4 : ℕ) = 2 * 2 from rfl, pow_mul, s2_sq]; norm_num

end P01

theorem P01.proof_10 : P01.stmt_10 := by
  have hs := P01.s2_ne
  have hb1 := P01.bI_ne
  have hb2 := P01.bnI_ne
  have cs := P01.conj_s2
  refine ⟨fun k hk => ?_, ?_, fun k hk => ?_⟩
  · have ha : P01.h k 2 0 = P01.rho 2 0 := (P01.lemD k 2 0 (by omega) le_rfl hk).1 (by omega)
    rw [P01.G_eq_form k 2 0 (by omega) rfl _ ha, P01.g20_form, P01.rho20]
    have c1 : (starRingEnd ℂ) (-Complex.I / (2 * P01.s2)) = Complex.I / (2 * P01.s2) := by
      rw [map_div₀, map_neg, Complex.conj_I, map_mul, cs, map_ofNat]; ring
    rw [c1]
    have p1 : -(-Complex.I / (2 * P01.s2)) = P01.s2 * Complex.I / 4 := by
      field_simp; linear_combination (-2 : ℂ) * P01.s2_sq
    have p2 : -(Complex.I / (2 * P01.s2)) = -(P01.s2 * Complex.I / 4) := by
      field_simp; linear_combination (2 : ℂ) * P01.s2_sq
    have p3 : -Complex.I / (2 * P01.s2) * (P01.bC - Complex.I)⁻¹ + Complex.I / (2 * P01.s2) * (P01.bC + Complex.I)⁻¹
        = P01.aC / 4 := by
      rw [P01.aC_eq]
      rw [P01.bC_eq] at hb1 hb2 ⊢
      field_simp
      ring_nf
      simp only [P01.s2_sq, P01.s2_pow3, P01.s2_pow4, Complex.I_sq, P01.I_pow_three']
      ring
    rw [p1, p2, p3]
  · have ha : P01.h 2 1 1 = P01.rho 1 1 := ((P01.lemD 2 1 1 le_rfl (by omega) le_rfl).2 rfl).1 rfl
    rw [P01.G_eq_form 2 1 1 le_rfl rfl _ ha, P01.g20_form, P01.rho11]
    have c1 : (starRingEnd ℂ) (Complex.I / (2 * P01.s2)) = -Complex.I / (2 * P01.s2) := by
      rw [map_div₀, Complex.conj_I, map_mul, cs, map_ofNat]
    rw [c1]
    have p1 : -(Complex.I / (2 * P01.s2)) = -(P01.s2 * Complex.I / 4) := by
      field_simp
      ring_nf
      simp only [P01.s2_sq, P01.s2_pow3, P01.s2_pow4, Complex.I_sq, P01.I_pow_three']
      ring
    have p2 : -(-Complex.I / (2 * P01.s2)) = P01.s2 * Complex.I / 4 := by
      field_simp
      ring_nf
      simp only [P01.s2_sq, P01.s2_pow3, P01.s2_pow4, Complex.I_sq, P01.I_pow_three']
      ring
    have p3 : Complex.I / (2 * P01.s2) * (P01.bC - Complex.I)⁻¹ + -Complex.I / (2 * P01.s2) * (P01.bC + Complex.I)⁻¹
        = -(P01.aC / 4) := by
      rw [P01.aC_eq]
      rw [P01.bC_eq] at hb1 hb2 ⊢
      field_simp
      ring_nf
      simp only [P01.s2_sq, P01.s2_pow3, P01.s2_pow4, Complex.I_sq, P01.I_pow_three']
      ring
    rw [p1, p2, p3]
    unfold P01.cst; simp only [map_neg]; ring
  · have ha : P01.h k 1 1 = P01.hhat 1 := ((P01.lemD k 1 1 le_rfl (by omega) (by omega)).2 rfl).2 (by omega)
    rw [P01.G_eq_form k 1 1 le_rfl rfl _ ha, P01.g11_form, P01.hhat1]
    have c1 : (starRingEnd ℂ) ((1 + Complex.I) / (2 * P01.s2)) = (1 - Complex.I) / (2 * P01.s2) := by
      rw [map_div₀, map_add, map_one, Complex.conj_I, map_mul, cs, map_ofNat]; ring
    rw [c1]
    have p1 : -((1 + Complex.I) / (2 * P01.s2)) = -(P01.s2 * (1 + Complex.I) / 4) := by
      field_simp
      ring_nf
      simp only [P01.s2_sq, P01.s2_pow3, P01.s2_pow4, Complex.I_sq, P01.I_pow_three']
      ring
    have p2 : -((1 - Complex.I) / (2 * P01.s2)) = -(P01.s2 * (1 - Complex.I) / 4) := by
      field_simp
      ring_nf
      simp only [P01.s2_sq, P01.s2_pow3, P01.s2_pow4, Complex.I_sq, P01.I_pow_three']
      ring
    have p3 : (1 + Complex.I) / (2 * P01.s2) * (P01.bC - Complex.I)⁻¹
        + (1 - Complex.I) / (2 * P01.s2) * (P01.bC + Complex.I)⁻¹ = -(P01.s2 / 4) := by
      rw [P01.bC_eq] at hb1 hb2 ⊢
      field_simp
      ring_nf
      simp only [P01.s2_sq, P01.s2_pow3, P01.s2_pow4, Complex.I_sq, P01.I_pow_three']
      ring
    rw [p1, p2, p3]


namespace P01

lemma Fold_mul_D1 (g : RatFunc ℂ) : Fold (g * D1) = -(Tg g * D2) := by
  have h := fold_D1
  rw [Fold_eq] at h ⊢
  rw [map_mul, mul_div_assoc, h]; ring

lemma ne_aux (c : ℂ) : (RatFunc.X : RatFunc ℂ) + RatFunc.C c ≠ 0 := by
  have := X_sub_C_ne (-c); rwa [map_neg, sub_neg_eq_add] at this

lemma X2_ne : (RatFunc.X : RatFunc ℂ) ^ 2 + 1 ≠ 0 := by
  have hI2 : RatFunc.C Complex.I * RatFunc.C Complex.I = (-1 : RatFunc ℂ) := by
    rw [← map_mul, Complex.I_mul_I, map_neg, map_one]
  have : (RatFunc.X : RatFunc ℂ) ^ 2 + 1 = (RatFunc.X - RatFunc.C Complex.I) * (RatFunc.X + RatFunc.C Complex.I) := by
    linear_combination hI2
  rw [this]; exact mul_ne_zero (X_sub_C_ne _) (ne_aux _)

lemma Delta_g11 : Tg g11 - g11 = -(cst s2 * (1 - vX ^ 2)) / (2 * (vX ^ 2 + 1)) := by
  unfold g11
  simp only [vX, cst, map_div₀, map_neg, map_mul, map_add, map_sub, map_pow, map_one, map_ofNat, Tg_X, sig_C]
  have hX := X_ne_zero'
  have h2 := X2_ne
  have h3 : (RatFunc.X : RatFunc ℂ)⁻¹ ^ 2 + 1 ≠ 0 := by
    intro h; apply h2
    have : (RatFunc.X : RatFunc ℂ) ^ 2 + 1 = RatFunc.X ^ 2 * ((RatFunc.X : RatFunc ℂ)⁻¹ ^ 2 + 1) := by
      rw [mul_add, ← mul_pow, mul_inv_cancel₀ hX]; ring
    rw [this, h, mul_zero]
  have h2' : 1 + (RatFunc.X : RatFunc ℂ) ^ 2 ≠ 0 := by rwa [add_comm]
  field_simp
  ring

lemma Delta_g20 : Tg (-g20) - (-g20) = -(cst s2 * (1 - vX ^ 2)) / (2 * (vX ^ 2 + 1)) := by
  unfold g20
  simp only [vX, cst, map_div₀, map_neg, map_mul, map_add, map_sub, map_pow, map_one, map_ofNat, Tg_X, sig_C]
  have hX := X_ne_zero'
  have h2 := X2_ne
  have h3 : (RatFunc.X : RatFunc ℂ)⁻¹ ^ 2 + 1 ≠ 0 := by
    intro h; apply h2
    have : (RatFunc.X : RatFunc ℂ) ^ 2 + 1 = RatFunc.X ^ 2 * ((RatFunc.X : RatFunc ℂ)⁻¹ ^ 2 + 1) := by
      rw [mul_add, ← mul_pow, mul_inv_cancel₀ hX]; ring
    rw [this, h, mul_zero]
  have h2' : 1 + (RatFunc.X : RatFunc ℂ) ^ 2 ≠ 0 := by rwa [add_comm]
  field_simp
  ring

lemma main_id : -(2 * g20 * D1) + (-(cst s2 * (1 - vX ^ 2)) / (2 * (vX ^ 2 + 1))) * D2 = Rwt := by
  unfold g20 D1 D2 Rwt
  have hX := X_ne_zero'
  have h2 := X2_ne
  have ha := X_sub_C_ne aC
  have hb := X_sub_C_ne bC
  have ha' := ne_aux aC
  have hb' := ne_aux bC
  have h2' : 1 + (RatFunc.X : RatFunc ℂ) ^ 2 ≠ 0 := by rwa [add_comm]
  have hS : RatFunc.C s2 ^ 2 = (2 : RatFunc ℂ) := by rw [← map_pow, s2_sq, map_ofNat]
  have hA : RatFunc.C aC = RatFunc.C s2 + 1 := by rw [aC_eq, map_add, map_one, add_comm]
  have hB : RatFunc.C bC = RatFunc.C s2 - 1 := by rw [bC_eq, map_sub, map_one]
  simp only [vX, cst, map_div₀, map_ofNat]
  field_simp
  rw [hA, hB]
  ring_nf
  have h3 : RatFunc.C s2 ^ 3 = 2 * RatFunc.C s2 := by rw [pow_succ, hS]
  have h4 : RatFunc.C s2 ^ 4 = (4 : RatFunc ℂ) := by rw [show 4 = 2 * 2 from rfl, pow_mul, hS]; norm_num
  have h5 : RatFunc.C s2 ^ 5 = 4 * RatFunc.C s2 := by rw [pow_succ, h4]
  have h6 : RatFunc.C s2 ^ 6 = (8 : RatFunc ℂ) := by rw [pow_succ, h5, mul_assoc, ← sq, hS]; norm_num
  have h7 : RatFunc.C s2 ^ 7 = 8 * RatFunc.C s2 := by rw [pow_succ, h6]
  have h8 : RatFunc.C s2 ^ 8 = (16 : RatFunc ℂ) := by rw [pow_succ, h7, mul_assoc, ← sq, hS]; norm_num
  simp only [hS, h3, h4, h5, h6, h7, h8]
  ring

end P01
namespace P01

lemma proof_10' : stmt_10 := P01.proof_10

lemma W10_eq (k : ℕ) (hk : 2 ≤ k) :
    W k 1 0 = -(2 * g20 * D1) + (Tg (G k 1 1) - G k 1 1) * D2 := by
  rw [W_eq k 1 0 (by omega)]
  have hR : Rinit k 1 0 = 0 := by unfold Rinit; rw [if_neg (by omega)]
  have hG1 : Jp (W k (1 + 1) 0) = g20 := by
    rw [← proof_10'.1 k hk]; unfold G; rw [if_pos (by norm_num)]
  have hG2 : Jp (W k 1 (0 + 1)) = G k 1 1 := by unfold G; rw [if_pos (by norm_num)]
  have hG3 : Jp (W k 1 1) = G k 1 1 := by unfold G; rw [if_pos (by norm_num)]
  rw [hR, hG1, hG2, if_pos (by norm_num), dif_pos rfl, Fold_mul_D1]
  unfold cst
  simp only [Nat.cast_one, Nat.cast_zero, zero_add, map_one, show (1 : ℂ) + 1 = 2 by norm_num, map_ofNat]
  ring

end P01

theorem P01.proof_11 : P01.stmt_11 := by
  intro k hk
  rw [P01.W10_eq k hk, ← P01.main_id]
  congr 2
  rcases eq_or_lt_of_le hk with h | h
  · subst h; rw [P01.proof_10'.2.1, P01.Delta_g20]
  · rw [P01.proof_10'.2.2 k (by omega), P01.Delta_g11]


namespace P01

lemma res_term (c d α : ℂ) : res c (cst α * (vX - cst d) ^ (-1 : ℤ)) = if c = d then α else 0 := by
  unfold res
  split_ifs with h
  · subst h
    rw [lcoeff_cst_mul]
    have := lcoeff_pole c 1 (-1)
    simp only [Nat.cast_one] at this
    rw [this]; simp
  · have := lcoeff_other c d h 1 (-1) (by norm_num)
    simp only [Nat.cast_one] at this
    rw [lcoeff_cst_mul, this, mul_zero]

lemma Rwt_sum : Rwt = cst (1 / 2) * (vX - cst 0) ^ (-1 : ℤ) + cst (1 / 2) * (vX - cst aC) ^ (-1 : ℤ)
    + cst (-1) * (vX - cst (-aC)) ^ (-1 : ℤ) + cst (-1 / 2) * (vX - cst (-bC)) ^ (-1 : ℤ) := by
  have c1 : RatFunc.C (1 / 2 : ℂ) = (1 / 2 : RatFunc ℂ) := by rw [map_div₀, map_one, map_ofNat]
  have c2 : RatFunc.C (-1 / 2 : ℂ) = (-1 / 2 : RatFunc ℂ) := by rw [map_div₀, map_neg, map_one, map_ofNat]
  unfold Rwt cst
  rw [c1, c2]
  simp only [zpow_neg_one, map_zero, sub_zero, map_neg, sub_neg_eq_add, map_one, one_div, mul_inv]
  ring

lemma Y01_b_pos : 0 < Y01.b := by
  unfold Y01.b
  have : (1 : ℝ) < Real.sqrt 2 := by
    rw [show (1 : ℝ) = Real.sqrt 1 by simp]; exact Real.sqrt_lt_sqrt (by norm_num) (by norm_num)
  linarith

lemma Y01_a_eq : Y01.a = Y01.b + 2 := by unfold Y01.a Y01.b; ring

lemma ofReal_ne {x y : ℝ} (h : x ≠ y) : (x : ℂ) ≠ (y : ℂ) := fun e => h (Complex.ofReal_injective e)

lemma simplePart_Rwt : simplePart Rwt = Rwt := by
  have hb := Y01_b_pos
  have ha := Y01_a_eq
  have d1 : (0 : ℂ) ≠ aC := by unfold aC; exact_mod_cast (show (0 : ℝ) ≠ Y01.a by linarith)
  have d2 : (0 : ℂ) ≠ -aC := by
    unfold aC; rw [← Complex.ofReal_neg]; exact_mod_cast (show (0 : ℝ) ≠ -Y01.a by linarith)
  have d3 : (0 : ℂ) ≠ -bC := by
    unfold bC; rw [← Complex.ofReal_neg]; exact_mod_cast (show (0 : ℝ) ≠ -Y01.b by linarith)
  have d4 : aC ≠ -aC := by
    unfold aC; rw [← Complex.ofReal_neg]; exact ofReal_ne (by linarith)
  have d5 : aC ≠ -bC := by
    unfold aC bC; rw [← Complex.ofReal_neg]; exact ofReal_ne (by linarith)
  have d6 : -aC ≠ -bC := by
    unfold aC bC; rw [← Complex.ofReal_neg, ← Complex.ofReal_neg]; exact ofReal_ne (by linarith)
  have d7 : bC ≠ 0 := by unfold bC; exact_mod_cast (show Y01.b ≠ 0 by linarith)
  have d8 : bC ≠ aC := by unfold aC bC; exact ofReal_ne (by linarith)
  have d9 : bC ≠ -aC := by unfold aC bC; rw [← Complex.ofReal_neg]; exact ofReal_ne (by linarith)
  have d10 : bC ≠ -bC := by unfold bC; rw [← Complex.ofReal_neg]; exact ofReal_ne (by linarith)
  have hres : ∀ c : ℂ, res c Rwt = (if c = 0 then 1 / 2 else 0) + (if c = aC then 1 / 2 else 0)
      + (if c = -aC then -1 else 0) + (if c = -bC then -1 / 2 else 0) := by
    intro c
    rw [Rwt_sum]
    unfold res
    rw [lcoeff_add, lcoeff_add, lcoeff_add]
    have e := fun d α => res_term c d α
    unfold res at e
    simp only [e]
  unfold simplePart
  rw [hres, hres, hres, hres, hres]
  simp only [d1, d2, d3, d4, d5, d6, d7, d8, d9, d10, d1.symm, d2.symm, d3.symm, d4.symm, d5.symm,
    d6.symm, if_true, if_false, eq_self_iff_true, add_zero, zero_add, neg_eq_zero, neg_inj]
  conv_rhs => rw [Rwt_sum]
  unfold cst
  simp only [zpow_neg_one, map_zero, sub_zero, map_neg, sub_neg_eq_add, map_div₀, map_one, map_ofNat,
    zero_div, add_zero]
  ring

end P01
namespace P01

lemma ev_Rwt (v : ℝ) (h0 : v ≠ 0) (ha : v ≠ Y01.a) (ha' : v ≠ -Y01.a) (hb' : v ≠ -Y01.b) :
    ev Rwt (v : ℂ) = ((Y01.Rwt v : ℝ) : ℂ) := by
  have n0 : (v : ℂ) ≠ 0 := by exact_mod_cast h0
  have n1 : (v : ℂ) ≠ aC := by unfold aC; exact ofReal_ne ha
  have n2 : (v : ℂ) ≠ -aC := by unfold aC; rw [← Complex.ofReal_neg]; exact ofReal_ne ha'
  have n3 : (v : ℂ) ≠ -bC := by unfold bC; rw [← Complex.ofReal_neg]; exact ofReal_ne hb'
  rw [Rwt_sum]
  rw [ev_of_NV (NV_add (NV_add (NV_add (NV_mul (NV_cst _ _) (NV_zpow _ _ (by simpa using n0) _))
    (NV_mul (NV_cst _ _) (NV_zpow _ _ n1 _))) (NV_mul (NV_cst _ _) (NV_zpow _ _ n2 _)))
    (NV_mul (NV_cst _ _) (NV_zpow _ _ n3 _)))]
  unfold Y01.Rwt
  have e1 : ((Y01.a : ℝ) : ℂ) = aC := rfl
  have e2 : ((Y01.b : ℝ) : ℂ) = bC := rfl
  push_cast
  rw [e1, e2]
  simp only [zpow_neg_one, sub_zero, sub_neg_eq_add, one_div, mul_inv]
  ring

lemma proof_11' : stmt_11 := P01.proof_11

lemma eqOn_Rwt : Set.EqOn (fun v : ℝ => ev Rwt (v : ℂ) * ((Y01.l1 v : ℝ) : ℂ))
    (fun v : ℝ => (((Y01.l1 v * Y01.Rwt v) : ℝ) : ℂ)) (Set.uIoo Y01.b Y01.a) := by
  have hb := Y01_b_pos
  have ha := Y01_a_eq
  intro v hv
  rw [Set.uIoo_of_lt (by linarith)] at hv
  obtain ⟨h1, h2⟩ := hv
  simp only
  rw [ev_Rwt v (by linarith) (by linarith) (by linarith) (by linarith)]
  push_cast; ring

end P01

theorem P01.proof_12 : P01.stmt_12 := by
  intro k hk
  have h11 := P01.proof_11' k hk
  have hY := Y01.proof_7.2.1
  have hT := Y01.proof_8.1
  refine ⟨h11, by rw [h11, P01.simplePart_Rwt], ?_, ?_⟩
  · rw [h11, P01.simplePart_Rwt]
    have hc : IntervalIntegrable (fun v : ℝ => (((Y01.l1 v * Y01.Rwt v) : ℝ) : ℂ)) MeasureTheory.volume Y01.b Y01.a := by
      rw [intervalIntegrable_iff] at hY ⊢
      exact hY.ofReal
    exact hc.congr_uIoo (P01.eqOn_Rwt.symm)
  · unfold P01.W2
    rw [h11, P01.simplePart_Rwt, intervalIntegral.integral_congr_uIoo P01.eqOn_Rwt,
      intervalIntegral.integral_ofReal]
    rw [← hT]
    rfl

theorem P01.proof_13 : P01.stmt_13 := by
  refine ⟨P01.Inv_add, P01.Inv_cst_mul, P01.Inv_Inv, P01.Inv_Fold, ?_, ?_, ?_, P01.Hp_Inv⟩
  · have hab : P01.bC * P01.aC = 1 := by
      rw [P01.aC_eq, P01.bC_eq]; linear_combination P01.s2_sq
    exact P01.Inv_D _ _ hab
  · have hab : P01.aC * P01.bC = 1 := by
      rw [P01.aC_eq, P01.bC_eq]; linear_combination P01.s2_sq
    exact P01.Inv_D _ _ hab
  · intro m _
    refine ⟨P01.Inv_pole _ Complex.I_sq m, ?_⟩
    have e : P01.vX + P01.cst Complex.I = P01.vX - P01.cst (-Complex.I) := by
      unfold P01.cst; rw [map_neg, sub_neg_eq_add]
    rw [e]
    exact P01.Inv_pole _ (by rw [neg_sq, Complex.I_sq]) m


namespace P01

open Polynomial

/-! ### The involution and the pairing `W + 𝓘(W)` (g1:lem:par, g1:lem:res) -/

lemma Inv_neg (f : RatFunc ℂ) : Inv (-f) = -Inv f := by
  rw [Inv_eq, Inv_eq, map_neg]; ring

lemma Inv_sub (f g : RatFunc ℂ) : Inv (f - g) = Inv f - Inv g := by
  rw [sub_eq_add_neg, Inv_add, Inv_neg, ← sub_eq_add_neg]

lemma Inv_zero : Inv 0 = 0 := by rw [Inv_eq, map_zero, zero_div]

lemma Inv_mul (A B : RatFunc ℂ) : Inv (A * B) = Sg A * Inv B := by
  rw [Inv_eq, Inv_eq, map_mul, mul_div_assoc]

lemma Inv_mul' (A B : RatFunc ℂ) : Inv (A * B) = Inv A * Sg B := by
  rw [Inv_eq, Inv_eq, map_mul, mul_div_right_comm]

lemma Inv_div_cst (f : RatFunc ℂ) (x : ℂ) : Inv (f / cst x) = Inv f / cst x := by
  rw [Inv_eq, Inv_eq, map_div₀, cst, sig_C, div_right_comm]

lemma Fold_add (f g : RatFunc ℂ) : Fold (f + g) = Fold f + Fold g := by
  rw [Fold_eq, Fold_eq, Fold_eq, map_add, add_div]

lemma Inv_D1 : Inv D1 = D1 := by
  have hab : bC * aC = 1 := by rw [aC_eq, bC_eq]; linear_combination s2_sq
  exact Inv_D _ _ hab

lemma Inv_D2 : Inv D2 = D2 := by
  have hab : aC * bC = 1 := by rw [aC_eq, bC_eq]; linear_combination s2_sq
  exact Inv_D _ _ hab

lemma ev_zero (z : ℂ) : ev 0 z = 0 := by simp [ev]

lemma ne_zero_of_NV {z : ℂ} {f : RatFunc ℂ} {v : ℂ} (h : NV z f v) (hv : v ≠ 0) : f ≠ 0 := by
  rintro rfl; exact hv ((ev_of_NV h).symm.trans (ev_zero z))

/-- if `W + 𝓘(W)` has no pole at `±i`, then `Hp(W) + 𝓘(Hp(W)) = 0`. -/
lemma Hp_par (f : RatFunc ℂ) (hI : poleLE Complex.I (f + Inv f) 0)
    (hJ : poleLE (-Complex.I) (f + Inv f) 0) : Hp f + Inv (Hp f) = 0 := by
  rw [← (Hp_Inv f).1]
  apply no_poles_zero
  · intro c n hn
    have hn' : n < 0 := by simpa using hn
    by_cases hc : c = Complex.I ∨ c = -Complex.I
    · rw [lcoeff_add]
      rcases lt_or_eq_of_le (show n ≤ -1 by omega) with h1 | h1
      · have a1 := sub_Hp_pole f c hc n h1
        have a2 := sub_Hp_pole (Inv f) c hc n h1
        rw [lcoeff_sub] at a1 a2
        have a3 : lcoeff c (f + Inv f) n = 0 := by
          rcases hc with rfl | rfl
          · exact hI n hn
          · exact hJ n hn
        rw [lcoeff_add] at a3
        linear_combination -a1 - a2 + a3
      · subst h1
        have e1 := Hp_res f c hc
        have e2 := Hp_res (Inv f) c hc
        unfold res at e1 e2
        rw [e1, e2, add_zero]
    · push Not at hc
      rw [lcoeff_add, Hp_regular f c hc.1 hc.2 n hn, Hp_regular (Inv f) c hc.1 hc.2 n hn, add_zero]
  · exact small_add _ _ (Hp_small _) (Hp_small _)

lemma NV_Sg_pt {w : ℂ} (hw : w ≠ 0) {h : RatFunc ℂ} {v : ℂ} (hh : NV (-w⁻¹) h v) : NV w (Sg h) v := by
  have hX : NV w (-(RatFunc.X : RatFunc ℂ)⁻¹) (-w⁻¹) := by
    have := NV_neg (NV_X_zpow w hw (-1))
    simpa [zpow_neg_one] using this
  exact NV_sig _ _ hX hh

noncomputable def JVt (f : RatFunc ℂ) (w : ℂ) : ℂ := JV Complex.I f w + JV (-Complex.I) f w
noncomputable def HVt (f : RatFunc ℂ) (w : ℂ) : ℂ := HV Complex.I f w + HV (-Complex.I) f w

lemma NV_Jp (f : RatFunc ℂ) {z : ℂ} (h1 : z ≠ Complex.I) (h2 : z ≠ -Complex.I) :
    NV z (Jp f) (JVt f z) := NV_add (NV_jpart _ f _ h1) (NV_jpart _ f _ h2)

lemma NV_Hp (f : RatFunc ℂ) {z : ℂ} (h1 : z ≠ Complex.I) (h2 : z ≠ -Complex.I) :
    NV z (Hp f) (HVt f z) := NV_add (NV_ppart _ f _ h1) (NV_ppart _ f _ h2)

lemma JVt_deriv (f : RatFunc ℂ) (z : ℂ) (h1 : z ≠ Complex.I) (h2 : z ≠ -Complex.I) :
    HasDerivAt (JVt f) (HVt f z) z := by
  unfold JVt HVt JV HV
  apply HasDerivAt.add
  · exact HasDerivAt.fun_sum (fun m hm => term_deriv _ _ z m (Finset.mem_Icc.1 hm).1 h1)
  · exact HasDerivAt.fun_sum (fun m hm => term_deriv _ _ z m (Finset.mem_Icc.1 hm).1 h2)

lemma ninv_ne_I {x : ℂ} (hx0 : x ≠ 0) (hx : x ≠ Complex.I) : -x⁻¹ ≠ Complex.I := by
  intro h; apply hx
  have e : Complex.I * x + 1 = 0 := by
    rw [← h]; field_simp; ring
  linear_combination (-Complex.I) * e + x * Complex.I_sq

lemma ninv_ne_negI {x : ℂ} (hx0 : x ≠ 0) (hx : x ≠ -Complex.I) : -x⁻¹ ≠ -Complex.I := by
  intro h; apply hx
  have e : Complex.I * x - 1 = 0 := by
    have h' : x⁻¹ = Complex.I := by linear_combination -h
    rw [← h']; field_simp; ring
  linear_combination (-Complex.I) * e + x * Complex.I_sq

lemma eq_zero_of_NV (f : RatFunc ℂ) (U : Set ℂ) (hU : U.Infinite) (h : ∀ x ∈ U, NV x f 0) : f = 0 := by
  rw [← RatFunc.num_eq_zero_iff]
  apply Polynomial.eq_zero_of_infinite_isRoot
  apply hU.mono
  intro x hx
  obtain ⟨p, q, hq, hf, hv⟩ := h x hx
  have hq0 := ne_of_eval q x hq
  have hp : p.eval x = 0 := by
    rcases div_eq_zero_iff.1 hv.symm with h | h
    · exact h
    · exact absurd h hq
  have hcross : f.num * q = p * f.denom := by
    apply RatFunc.algebraMap_injective ℂ
    have h1 := (RatFunc.num_div_denom f).trans hf
    rw [div_eq_div_iff (alg_ne _ f.denom_ne_zero) (alg_ne _ hq0)] at h1
    rw [map_mul, map_mul]; exact h1
  have := congrArg (eval x) hcross
  simp only [eval_mul, hp, zero_mul] at this
  exact (mul_eq_zero.1 this).resolve_right hq

/-- "the function `G(v) + G(-1/v)` is constant": if `W + 𝓘(W)` has no pole at `±i`, then
`Jp(W) + Jp(W)(-1/v)` is a constant (its derivative `Hp(W) + 𝓘(Hp(W))` vanishes). -/
lemma par_const (f : RatFunc ℂ) (hI : poleLE Complex.I (f + Inv f) 0)
    (hJ : poleLE (-Complex.I) (f + Inv f) 0) : ∃ c : ℂ, Jp f + Sg (Jp f) = cst c := by
  have hH := Hp_par f hI hJ
  set U : Set ℂ := ({0, Complex.I, -Complex.I} : Set ℂ)ᶜ with hUdef
  have memU : ∀ x, x ∈ U ↔ x ≠ 0 ∧ x ≠ Complex.I ∧ x ≠ -Complex.I := by
    intro x; simp [hUdef]
  have hfin : ({0, Complex.I, -Complex.I} : Set ℂ).Finite := by
    exact (Set.toFinite _)
  have hUo : IsOpen U := hfin.isClosed.isOpen_compl
  have hUc : IsPreconnected U :=
    (Set.Countable.isConnected_compl_of_one_lt_rank (by simp) hfin.countable).isPreconnected
  have hUi : U.Infinite := hfin.infinite_compl
  have hder : ∀ x ∈ U, HasDerivAt (fun w => JVt f w + JVt f (-w⁻¹)) 0 x := by
    intro x hx
    obtain ⟨hx0, hx1, hx2⟩ := (memU x).1 hx
    have hi : HasDerivAt (fun w : ℂ => -w⁻¹) ((x ^ 2)⁻¹) x := by
      have := ((hasDerivAt_id x).inv hx0).fun_neg
      exact this.congr_deriv (by rw [neg_div, neg_neg, one_div]; rfl)
    have hJ : HasDerivAt (JVt f) (HVt f (-x⁻¹)) ((fun w : ℂ => -w⁻¹) x) :=
      JVt_deriv f (-x⁻¹) (ninv_ne_I hx0 hx1) (ninv_ne_negI hx0 hx2)
    have h2 : HasDerivAt (JVt f ∘ fun w : ℂ => -w⁻¹) (HVt f (-x⁻¹) * (x ^ 2)⁻¹) x :=
      HasDerivAt.comp x hJ hi
    have h3 := (JVt_deriv f x hx1 hx2).add h2
    -- the derivative vanishes by `Hp(W) + 𝓘(Hp(W)) = 0`
    have nv1 := NV_Hp f hx1 hx2
    have nv2 : NV x (Inv (Hp f)) (HVt f (-x⁻¹) / x ^ 2) := by
      rw [Inv_eq]
      exact NV_div (NV_Sg_pt hx0 (NV_Hp f (ninv_ne_I hx0 hx1) (ninv_ne_negI hx0 hx2)))
        (NV_Xpow x 2) (pow_ne_zero _ hx0)
    have := ev_of_NV (NV_add nv1 nv2)
    rw [hH, ev_zero] at this
    have hval : HVt f x + HVt f (-x⁻¹) * (x ^ 2)⁻¹ = 0 := by
      rw [this]; field_simp
    rw [hval] at h3; exact h3
  obtain ⟨a, ha⟩ := hUo.exists_is_const_of_deriv_eq_zero hUc
    (fun x hx => (hder x hx).differentiableAt.differentiableWithinAt)
    (fun x hx => (hder x hx).deriv)
  refine ⟨a, sub_eq_zero.1 (eq_zero_of_NV _ U hUi (fun x hx => ?_))⟩
  obtain ⟨hx0, hx1, hx2⟩ := (memU x).1 hx
  have := NV_sub (NV_add (NV_Jp f hx1 hx2)
    (NV_Sg_pt hx0 (NV_Jp f (ninv_ne_I hx0 hx1) (ninv_ne_negI hx0 hx2)))) (NV_cst x a)
  have e : JVt f x + JVt f (-x⁻¹) - a = 0 := by
    rw [ha x hx, sub_self]
  rwa [e] at this

/-! ### Invariance of the initial data -/

lemma X_ne : (RatFunc.X : RatFunc ℂ) ≠ 0 := X_ne_zero'

lemma Sg_tR : Sg tR = tR := by
  unfold tR vX
  simp only [map_div₀, map_add, map_mul, map_pow, map_sub, map_one, map_ofNat, Sg_X]
  have := X_ne
  field_simp
  ring

lemma cst_ne {x : ℂ} (hx : x ≠ 0) : cst x ≠ 0 := by unfold cst; simpa using hx

lemma Sg_qR : Sg qR = -qR := by
  unfold qR vX
  have hC := cst_ne (mul_ne_zero two_ne_zero s2_ne)
  unfold cst at hC ⊢
  simp only [map_div₀, map_add, map_mul, map_pow, map_one, Sg_X, sig_C]
  have := X_ne
  field_simp
  ring

def DInv (F : PowerSeries (PowerSeries (RatFunc ℂ))) : Prop :=
  ∀ r s : ℕ, Sg (PowerSeries.coeff s (PowerSeries.coeff r F)) = PowerSeries.coeff s (PowerSeries.coeff r F)

lemma DInv_add {F G : PowerSeries (PowerSeries (RatFunc ℂ))} (hF : DInv F) (hG : DInv G) :
    DInv (F + G) := fun r s => by
  simp only [map_add]; rw [hF r s, hG r s]

lemma DInv_mul {F G : PowerSeries (PowerSeries (RatFunc ℂ))} (hF : DInv F) (hG : DInv G) :
    DInv (F * G) := fun r s => by
  rw [PowerSeries.coeff_mul, map_sum, map_sum]
  refine Finset.sum_congr rfl (fun p _ => ?_)
  rw [PowerSeries.coeff_mul, map_sum]
  refine Finset.sum_congr rfl (fun q _ => ?_)
  rw [map_mul, hF, hG]

lemma DInv_one : DInv 1 := fun r s => by
  rw [PowerSeries.coeff_one]
  split_ifs
  · rw [PowerSeries.coeff_one]; split_ifs <;> simp
  · simp

lemma DInv_pow {F : PowerSeries (PowerSeries (RatFunc ℂ))} (hF : DInv F) (j : ℕ) : DInv (F ^ j) := by
  induction j with
  | zero => simpa using DInv_one
  | succ j ih => rw [pow_succ]; exact DInv_mul ih hF

lemma DInv_CC {A : RatFunc ℂ} (hA : Sg A = A) : DInv (PowerSeries.C (PowerSeries.C A)) := fun r s => by
  rw [PowerSeries.coeff_C]
  split_ifs
  · rw [PowerSeries.coeff_C]; split_ifs
    · exact hA
    · simp
  · simp

lemma DInv_E2x : DInv (PowerSeries.map ((PowerSeries.C).comp (algebraMap ℂ (RatFunc ℂ))) E2) :=
  fun r s => by
  rw [PowerSeries.coeff_map, RingHom.comp_apply, PowerSeries.coeff_C]
  split_ifs
  · exact AlgHom.commutes Sg _
  · simp

lemma DInv_E2y : DInv (PowerSeries.C (PowerSeries.map (algebraMap ℂ (RatFunc ℂ)) E2)) := fun r s => by
  rw [PowerSeries.coeff_C]
  split_ifs
  · rw [PowerSeries.coeff_map]; exact AlgHom.commutes Sg _
  · simp

lemma Sg_Ucoef (j r s : ℕ) : Sg (Ucoef j r s) = Ucoef j r s := by
  have hV : DInv Vser := by
    unfold Vser
    refine DInv_add (DInv_mul (DInv_CC ?_) DInv_E2x) (DInv_mul (DInv_CC ?_) DInv_E2y)
    · rw [map_pow, map_sub, map_one, Sg_tR]
    · rw [map_pow, Sg_tR]
  exact DInv_pow hV j r s

lemma tR_ne : tR ≠ 0 := by
  have hc : Complex.I ≠ 0 := Complex.I_ne_zero
  refine ne_zero_of_NV (NV_tR Complex.I hc) ?_
  refine div_ne_zero ?_ (by simp [hc])
  intro h
  have : Complex.I = 1 := by linear_combination (h - Complex.I_sq) / 2
  have := congrArg Complex.re this; simp at this

lemma one_sub_tR_ne : 1 - tR ≠ 0 := by
  have hc : Complex.I ≠ 0 := Complex.I_ne_zero
  refine ne_zero_of_NV (NV_sub (NV_one _) (NV_tR Complex.I hc)) ?_
  intro h
  have h' : Complex.I ^ 2 + 2 * Complex.I - 1 = 4 * Complex.I := by
    field_simp at h; linear_combination -h
  have : Complex.I = -1 := by linear_combination (h' - Complex.I_sq) / (-2)
  have := congrArg Complex.im this; simp at this

lemma Inv_pre : Inv (1 / (cst (2 * s2) * vX * tR * (1 - tR))) = -(1 / (cst (2 * s2) * vX * tR * (1 - tR))) := by
  have hC := cst_ne (mul_ne_zero two_ne_zero s2_ne)
  have h1 := tR_ne
  have h2 := one_sub_tR_ne
  rw [Inv_eq]
  unfold vX cst at *
  simp only [map_div₀, map_mul, map_sub, map_one, sig_C, Sg_X, Sg_tR]
  have := X_ne
  field_simp

lemma Sg_sum_Rt (r s : ℕ) :
    Sg (∑ j ∈ Finset.Icc 1 (r + s), cst (binHalf j) * qR ^ (2 - 2 * (j : ℤ)) * Ucoef j r s)
      = ∑ j ∈ Finset.Icc 1 (r + s), cst (binHalf j) * qR ^ (2 - 2 * (j : ℤ)) * Ucoef j r s := by
  rw [map_sum]
  refine Finset.sum_congr rfl (fun j _ => ?_)
  rw [map_mul, map_mul, map_zpow₀, Sg_qR, Sg_Ucoef, cst, sig_C]
  rw [show (2 - 2 * (j : ℤ)) = 2 * (1 - j) by ring, zpow_mul, zpow_mul, zpow_two, zpow_two, neg_mul_neg]

lemma Rt_par (r s : ℕ) (c : ℂ) (hc : c = Complex.I ∨ c = -Complex.I) :
    poleLE c (Rt r s + Inv (Rt r s)) 0 := by
  have hE1 : Reg c (if s = 0 then -D1 / cst (2 * ((r + s).factorial : ℂ)) else 0) := by
    refine Reg_ite c _ ⟨_, NV_div (NV_neg (D1_zero c hc)) (NV_cst c _) ?_⟩
    exact mul_ne_zero two_ne_zero (by exact_mod_cast (Nat.factorial_pos _).ne')
  have hE2 : Reg c (if r = 0 then D2 / cst (2 * ((r + s).factorial : ℂ)) else 0) := by
    refine Reg_ite c _ ⟨_, NV_div (D2_zero c hc) (NV_cst c _) ?_⟩
    exact mul_ne_zero two_ne_zero (by exact_mod_cast (Nat.factorial_pos _).ne')
  have iE1 : Inv (if s = 0 then -D1 / cst (2 * ((r + s).factorial : ℂ)) else 0)
      = (if s = 0 then -D1 / cst (2 * ((r + s).factorial : ℂ)) else 0) := by
    split_ifs
    · rw [Inv_div_cst, Inv_neg, Inv_D1]
    · exact Inv_zero
  have iE2 : Inv (if r = 0 then D2 / cst (2 * ((r + s).factorial : ℂ)) else 0)
      = (if r = 0 then D2 / cst (2 * ((r + s).factorial : ℂ)) else 0) := by
    split_ifs
    · rw [Inv_div_cst, Inv_D2]
    · exact Inv_zero
  have iM : Inv (1 / (cst (2 * s2) * vX * tR * (1 - tR)) *
      ∑ j ∈ Finset.Icc 1 (r + s), cst (binHalf j) * qR ^ (2 - 2 * (j : ℤ)) * Ucoef j r s)
      = -(1 / (cst (2 * s2) * vX * tR * (1 - tR)) *
      ∑ j ∈ Finset.Icc 1 (r + s), cst (binHalf j) * qR ^ (2 - 2 * (j : ℤ)) * Ucoef j r s) := by
    rw [Inv_mul', Inv_pre, Sg_sum_Rt, neg_mul]
  have e : Rt r s + Inv (Rt r s)
      = cst 2 * (if s = 0 then -D1 / cst (2 * ((r + s).factorial : ℂ)) else 0)
        + cst 2 * (if r = 0 then D2 / cst (2 * ((r + s).factorial : ℂ)) else 0) := by
    unfold Rt
    rw [Inv_add, Inv_add, iE1, iE2, iM]
    unfold cst; rw [map_ofNat]; ring
  rw [e]
  exact poleLE_of_reg (Reg_add (Reg_mul (Reg_cst c 2) hE1) (Reg_mul (Reg_cst c 2) hE2)) 0

lemma Rinit_par (k r s : ℕ) (c : ℂ) (hc : c = Complex.I ∨ c = -Complex.I) :
    poleLE c (Rinit k r s + Inv (Rinit k r s)) 0 := by
  unfold Rinit
  split_ifs with h1 h2
  · have e : Rt r s + Fold (Rt s r) + Inv (Rt r s + Fold (Rt s r))
        = (Rt r s + Inv (Rt r s)) + Fold (Rt s r + Inv (Rt s r)) := by
      rw [Inv_add, Inv_Fold, Fold_add]; ring
    rw [e]
    exact poleLE_add (Rt_par r s c hc) (Fold_pole c (sqI c hc) _ 0 (Rt_par s r (-c) (neg_mem c hc)))
  · exact Rt_par r s c hc
  · rw [Inv_zero, add_zero]; exact poleLE_zero' c 0

lemma pair_term (a : ℂ) (J D : RatFunc ℂ) (hD : Inv D = D) (c0 : ℂ) (hJ : J + Sg J = cst c0) :
    cst a * J * D + Inv (cst a * J * D) = cst (a * c0) * D := by
  rw [mul_assoc, Inv_cst_mul, Inv_mul, hD]
  rw [show cst a * (J * D) + cst a * (Sg J * D) = cst a * ((J + Sg J) * D) by ring, hJ]
  unfold cst; rw [map_mul]; ring

lemma pair_fold (a : ℂ) (J : RatFunc ℂ) (c0 : ℂ) (hJ : J + Sg J = cst c0) :
    cst a * Fold (J * D1) + Inv (cst a * Fold (J * D1)) = cst a * Fold (cst c0 * D1) := by
  rw [Inv_cst_mul, Inv_Fold, Inv_mul, Inv_D1, ← mul_add, ← Fold_add, ← add_mul, hJ]

lemma pair_expand (R T1 T2 T3 : RatFunc ℂ) :
    R - T1 - T2 - T3 + Inv (R - T1 - T2 - T3)
      = (R + Inv R) - (T1 + Inv T1) - (T2 + Inv T2) - (T3 + Inv T3) := by
  rw [Inv_sub, Inv_sub, Inv_sub]; ring

/-- g1:lem:par, first claim, for every level and key. -/
lemma W_par (k : ℕ) : ∀ r s : ℕ, s ≤ r → 1 ≤ r + s → ∀ c : ℂ, (c = Complex.I ∨ c = -Complex.I) →
    poleLE c (W k r s + Inv (W k r s)) 0 := by
  suffices H : ∀ n : ℕ, ∀ r s : ℕ, k + 1 - (r + s) = n → s ≤ r → 1 ≤ r + s →
      ∀ c : ℂ, (c = Complex.I ∨ c = -Complex.I) → poleLE c (W k r s + Inv (W k r s)) 0 by
    intro r s h1 h2 c hc; exact H _ r s rfl h1 h2 c hc
  intro n
  induction n using Nat.strong_induction_on with
  | _ n ih =>
  intro r s hn hsr hd c hc
  by_cases hk : r + s ≤ k
  swap
  · rw [W_gt k r s (by omega), Inv_zero, add_zero]; exact poleLE_zero' c _
  have PC : ∀ r' s', s' ≤ r' → r' + s' = r + s + 1 →
      ∃ c0 : ℂ, Jp (W k r' s') + Sg (Jp (W k r' s')) = cst c0 := by
    intro r' s' h1 h2
    exact par_const _
      (ih (k + 1 - (r' + s')) (by omega) r' s' rfl h1 (by omega) _ (Or.inl rfl))
      (ih (k + 1 - (r' + s')) (by omega) r' s' rfl h1 (by omega) _ (Or.inr rfl))
  have regD : ∀ (c' x : ℂ) (D : RatFunc ℂ), NV c' D 0 → poleLE c' (cst x * D) 0 := fun c' x D hD =>
    poleLE_of_reg ⟨_, NV_mul (NV_cst _ _) hD⟩ 0
  rw [W_eq k r s hk, pair_expand]
  refine poleLE_sub (poleLE_sub (poleLE_sub (Rinit_par k r s c hc) ?_) ?_) ?_
  · obtain ⟨c0, hc0⟩ := PC (r + 1) s (by omega) (by omega)
    rw [pair_term _ _ _ Inv_D1 c0 hc0]
    exact regD _ _ _ (D1_zero c hc)
  · split_ifs with h
    · obtain ⟨c0, hc0⟩ := PC r (s + 1) (by omega) (by omega)
      rw [pair_term _ _ _ Inv_D2 c0 hc0]
      exact regD _ _ _ (D2_zero c hc)
    · rw [Inv_zero, add_zero]; exact poleLE_zero' c _
  · split_ifs with h
    · obtain ⟨c0, hc0⟩ := PC r r le_rfl (by omega)
      rw [pair_fold _ _ c0 hc0]
      apply poleLE_cst_mul
      exact Fold_pole c (sqI c hc) _ 0 (regD _ _ _ (D1_zero (-c) (neg_mem c hc)))
    · rw [Inv_zero, add_zero]; exact poleLE_zero' c _

lemma G_sub (k r s : ℕ) : subst (-vX⁻¹) (G k r s) = Sg (G k r s) := by
  unfold vX; rw [subst_eq _ hnXinv]

end P01

theorem P01.proof_14 : P01.stmt_14 := by
  intro k _ r s hsr hd
  refine ⟨P01.W_par k r s hsr hd _ (Or.inl rfl), P01.W_par k r s hsr hd _ (Or.inr rfl), fun h2 => ?_⟩
  rw [P01.G_sub]
  unfold P01.G
  simp only [h2, ↓reduceIte]
  exact P01.par_const _ (P01.W_par k r s hsr hd _ (Or.inl rfl)) (P01.W_par k r s hsr hd _ (Or.inr rfl))

namespace P01

lemma res_W (k r s : ℕ) (hsr : s ≤ r) (hd : 1 ≤ r + s) (c : ℂ) (hc : c = Complex.I ∨ c = -Complex.I) :
    res c (W k r s) = 0 := by
  have h0 := W_par k r s hsr hd c hc (-1) (by norm_num)
  have h1 : res c (Inv (W k r s)) = res c (W k r s) := by
    rcases hc with rfl | rfl
    · exact (Hp_Inv _).2.1
    · exact (Hp_Inv _).2.2
  unfold res at h1 ⊢
  rw [lcoeff_add, h1] at h0
  linear_combination h0 / 2

lemma pole1_res0 {c : ℂ} {f : RatFunc ℂ} (h1 : poleLE c f 1) (h2 : res c f = 0) : poleLE c f 0 := by
  intro n hn
  rcases lt_or_eq_of_le (show n ≤ -1 by omega) with h | h
  · exact h1 n (by omega)
  · subst h; exact h2

lemma S_reg (k r s : ℕ) (hsr : s ≤ r) (hd : 1 ≤ r + s) (c : ℂ) (hc : c = Complex.I ∨ c = -Complex.I) :
    poleLE c (S k r s) 0 := by
  unfold S
  split_ifs with h
  · refine pole1_res0 (sub_Hp_pole _ c hc) ?_
    unfold res; rw [lcoeff_sub]
    have a := res_W k r s hsr hd c hc
    have b := Hp_res (W k r s) c hc
    unfold res at a b; rw [a, b, sub_zero]
  · have := W_bound k r s hsr hd c hc
    rwa [show 2 * (r + s) - 2 = 0 by omega] at this

end P01

theorem P01.proof_15 : P01.stmt_15 := by
  intro k _ r s hsr hd
  exact ⟨P01.res_W k r s hsr hd _ (Or.inl rfl), P01.res_W k r s hsr hd _ (Or.inr rfl),
    P01.S_reg k r s hsr hd _ (Or.inl rfl), P01.S_reg k r s hsr hd _ (Or.inr rfl)⟩

