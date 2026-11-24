# G9 targeted second referee check

2026-11-23. Supplied sources and scripts were not edited. Line numbers below refer to this supplied revision. New checks are in `referee_g9b_checks.py`; outputs are `referee_g9b_algebra.out`, `referee_g9b_coincident.out`, and `referee_g9b_oct_rerun.out`.

The arbitrary-argument obstruction is correct. I found no counterexample to it, to the two obstruction coordinates, or to the repaired factor-120 calculation. The revision still needs a complete definition of the empty-word contribution to W_5, the short regularisation justification supplied below, and corrections to several explanatory sentences. One new/existing assertion in the revised hexagon explanation is literally false: the proposed independence of all thirteen specialised HPL values together with the zeta constants is already disproved at r=2. This does not invalidate the formal coefficient theorem.

## 1. Arbitrary arguments, complex coefficients, reversed order, and products

**Verdict: CORRECT.** `sec_g9_oct5.tex:85–100,132`.

The new Lemma `g9:lem:arb` supplies the missing argument. In detail:

- Tensoring the abelian group C^× with Q kills exactly its torsion. A rational dependence between [2] and [1+sqrt(2)] would, after clearing denominators, give a positive root of unity 2^a(1+sqrt(2))^b=1. Taking the quadratic norm gives 4^a(-1)^b=1, hence a=0 and then b=0. This proves the required Q-linear independence.
- A Q-linear retraction pi:U→V therefore exists. Its tensor power is well defined over Q. Extending scalars gives `id_C tensor pi^{tensor 5}` on `C tensor_Q U^{tensor 5}`. This fixes `C tensor_Q V^{tensor 5}`. The identification with `(C tensor_Q V)^{tensor_C 5}` used in Lemma zero is legitimate. In particular, extending scalars does not identify the formal classes X,Y with their numerical logarithms inside the one-dimensional C-vector space C.
- Applying pi in every slot commutes with every shuffle permutation, including the multiplication of external complex scalar coefficients.
- A classical symbol projects to `-a tensor b^{tensor 4}` or its reverse, regardless of the argument. Neither a nor b needs to be the logarithm of an S-unit, or to satisfy any relation after projection.
- In the forward shape, the first contraction contains B evaluated on b^2 and the second contains A evaluated on b^3; both vanish. In the reverse shape, these two vanishing factors occur in the opposite contractions. This works over C, not just over Q(sqrt(2)).
- The pairing on a word is exactly its coefficient in the degree-five Lie polynomial `[A,B]=AB-BA`. Lie polynomials annihilate proper shuffles. This also gives a shorter proof of Lemma zero(b), independently of the detailed deconcatenation explanation. It covers arbitrary tensors by bilinearity and arbitrary products of positive weights by grouping factors into two nonempty groups of total weight five.

Thus the pullback of the pairing along pi kills every generator of C_5 and fixes the target tensor. A nonzero target pairing rules out the proposed equality. The projection need not be canonical or computable; its existence suffices. No period-injectivity assertion is used here.

My exact checks annihilate all 128 ordered basis shuffles (all k=1,2,3,4) for both pairings, and symbolically annihilate both repeated-entry shapes with four independent indeterminate coefficients. Consequently these checks apply over C by polynomial identity, rather than by numerical sampling.

There is a harmless transcription error at `sec_g9_oct5.tex:92`. Exact replacement:

> By bilinearity we may take $s=u$ and $s'=v$, where $u,v$ are words of lengths $k,5-k$.

The present text instead says to take `s=u shuffle v`.

## 2. Lemma topwt and the lower-weight remainder

**Verdict: CORRECT AFTER REPAIR of one false parenthesis; the identity and weight bound are correct.** `sec_g9_oct5.tex:40–59`.

For `J_n(w)=int_0^1 t^n G(w;t)dt`, put `J_n(empty)=1/(n+1)`, of weight zero. The proof establishes

`J_n(w)=A(w,n) - (1/(n+1))*sum_{j=0}^n w_1^(n-j) J_j(tail(w))`.

The lower endpoint vanishes, including the all-zero word. The quotient polynomial has no pole at t=1. If w_1=1, `G(w;t)` has at most logarithmic growth, while `t^(n+1)-1=O(1-t)`, so the displayed limit at :57 is zero. If w_1≠1, its endpoint value converges for this alphabet.

Induction on word length now proves that `J_j(tail(w))` is a finite Q(sqrt(2))-linear combination of endpoint iterated integrals of lengths at most k-1, including the empty integral 1. Indeed, at each step A has the current word length, and the remaining J has a strictly shorter word. There is no infinite iteration in the polynomial exponent. This proves the advertised bound **at most** k-1. After multiplying by X^e Y^f with e+f+k=5, the remainder has weight at most four.

The parenthesis “then G(w;1) diverges” at :46 is false. Counterexample:

`G(1,0;1)=int_0^1 log(t)/(t-1)dt=zeta(2)`.

Nevertheless A((1,0),n)=0 because its prefactor tends to zero; actual divergence is unnecessary. Exact replacement for the parenthesis:

> (the endpoint value may diverge; in all cases $A(w,n)$ is the zero limit in \eqref{g9:eq:A}).

This corrects an explanatory assertion, not the lemma's integration identity.

## 3. Coincident-letter split-entry convention

**Verdict: CORRECT AFTER REPAIR by specifying unit tangents and giving the missing justification.** `sec_g9_oct5.tex:22–27`; implementation `scripts/g9_oct5_symbol.py:103–122`.

The implemented convention is consistent with the standard tangentially regularised symbol, including repeated internal letters and the endpoint collisions in intermediate recursion terms. The supplied weight-two checks alone do not prove this identification. It is not enough to declare an algorithm to be the symbol of an already defined integral.

Exact insertion/replacement for the unsupported sentence at :27:

```latex
We use the unit tangents in the coordinate $t$. The regularised
weight-one factor is
\[
 L(a;b;c)=\begin{cases}
 [c-b]-[a-b],&a\ne b,\ b\ne c,\\
 [c-b],&a=b\ne c,\\
 -[a-b],&a\ne b=c,\\
 0,&a=b=c.
 \end{cases}
\]
Indeed, take the constant term in $\log\varepsilon$ of
$\log((c-b+\varepsilon)/(a-b+\varepsilon))$.
The regularised iterated-integral coproduct gives
\[
 \Delta_{n-1,1}I(a_0;a_1,\ldots,a_n;a_{n+1})
 =\sum_i I(a_0;\ldots,\widehat a_i,\ldots;a_{n+1})
 \otimes I(a_{i-1};a_i;a_{i+1}).
\]
Iterating this formula and substituting $L$ proves the stated
split-entry rule, including coincidences created by deletions.
Signs disappear modulo torsion. Thus this recursion computes
the symbol of the tangentially regularised integral.
```

Cite the regularised coproduct here: [Goncharov, *Galois symmetries of fundamental groupoids and noncommutative geometry*, pp. 1–3, equation (5), Theorem 1.2](https://arxiv.org/pdf/math/0208144). The weight-one formula and coproduct explicitly cover coincident arguments. This supplies the identification missing from the manuscript, beyond declaring an algorithm to be its definition.

I independently implemented this regularised weight-one recursion with an explicit exact table of multiplicative classes. It agrees with the supplied implementation on **all 1,364 words of lengths 1–5** on {0,1,p,q}, including words with regularised endpoints. Moreover **all 5,008 ordered G-word shuffle identities of total weight at most five** pass, including multiplicities. This is an exhaustive finite implementation check; the coproduct argument supplies the mathematical identification with the regularised integral. No inconsistent coincident-letter example was found.

## 4. T_5 versus W_5

**Verdict: CORRECT AFTER REPAIR.** The renaming and subsequent obstruction proof are consistent, but the explicit definition is incomplete. `sec_g9_oct5.tex:38,62–70,75–81,111–145`.

T_5 at :38 is now the whole integral of the leading-log part. W_5 at :64–66 is intended to retain weight five, assigning constant logarithms weight one. The proof at :130–132, the block table, and the later limitations consistently discuss W_5. The distinction has not been lost in those arguments.

However, A(w,n) is defined only for nonempty words in Lemma topwt. The W_5 sum also includes rows with a=b=c=d=0, where the shuffle word is empty. These rows actually occur and have nonzero integrals. For example, the row

`(a,b,c,d,e,f)=(0,0,0,0,0,5)`

in the generated integrand `(3-2sqrt(2))*E_5^top/L` is exactly

`(1/40-sqrt(2)/60)*Y^5 = (3-2sqrt(2))*Y^5/120`,

with t-exponent zero. Its contribution to W_5 is therefore

`(24-16sqrt(2))*Y^5`.

This is a genuine weight-five product, not a lower-weight remainder. At :64 the present formula calls the undefined `A(empty,0)`; at :66 its expanded symbol formula also calls the nonexistent first letter of the empty word. The supplied code handles this case correctly. It pairs to zero, so the two obstruction numbers survive either way, but the full definition and symbol still need fixing.

Exact insertion before the W_5 formula:

> For the empty word we set $A(\emptyset,n)=1/(n+1)$, of weight zero, and $G(\emptyset;1)=1$. Its product with $X^eY^f$ is retained when $e+f=5$. Thus the formula below also includes the rows with $a=b=c=d=0$.

Exact replacement for the first sentence at :62:

> For a nonempty word $w$ of length $k$, the component of weight $k$ of $\int_0^1t^nG(w;t)\dd t$ is $A(w,n)$. After multiplication by $X^eY^f$ with $e+f+k=5$, its component of weight five is $X^eY^f A(w,n)$; the remaining terms have weight at most four.

The current sentence incorrectly calls A the “weight-five part” even when k<5. For example A((p),0)=(1-p)log(1-1/p) is a nonzero weight-one term, not a weight-five term.

Replace the symbol formula at :66 with the uniformly defined formula

\[
 \mathcal S W_5=\sum\gamma'\,
 \mathcal S(X^eY^f)\sqcup\!\sqcup\mathcal S A(w,n),
\]

where

\[
 \mathcal S A(w,n)=\begin{cases}
 1/(n+1),&w=\emptyset,\\
 0,&w\ne\emptyset,\ w_1=1,\\
 \dfrac{1-w_1^{n+1}}{n+1}\mathcal S G(w;1),&w\ne\emptyset,\ w_1\ne1.
 \end{cases}
\]

This also avoids writing zero times a divergent endpoint symbol.

Finally the quoted proposition's definition at :70 still describes the degree-five integrand terms rather than explicitly their weight projection. Exact replacement for its second sentence:

> Let $T_5$ be the contribution of the terms of total logarithmic degree five of the edge integrand to $P^{(5)}(1)/u$, and let $W_5$ be its weight-five part defined by the shuffle expansion and Lemma~\ref{g9:lem:topwt} above, with constant logarithms assigned weight one. All symbol assertions below refer to $\Sym W_5$.

The rest of the proposition can then refer unambiguously to “Its symbol”. These changes repair the formal definition without changing the computed pairings.

## 5. Lemma z5 and the factor 120

**Verdict: CORRECT.** `sec_g9_hex5.tex:129–137`, using `sec_g9_order5.tex:55–56,75–83,90–104`.

The bracket in Lemma precise5 is the edge integral before multiplication by 120. Its coefficients are now correctly identified as

\[
 c_5=24\int_0^1\sum_{a+b=5,\ a,b\ge1}g_{ab}(t)\,dt,
 \qquad
 c_{23}=-12\int_0^1(g_{32}+g_{23})(t)\,dt.
\]

The factors follow from the log-Beta terms

\[
 [z_5]\exp\Theta=\epsilon^4\delta+2\epsilon^3\delta^2+
 2\epsilon^2\delta^3+\epsilon\delta^4,
 \qquad
 [z_2z_3]\exp\Theta=-(\epsilon^3\delta^2+\epsilon^2\delta^3)
\]

at total degree five, multiplied by a!b!B_{n,0}. Thus the coefficients are 24B and -12B in exactly the stated cases. Multiplying c_5,c_23 by 120 gives the local-law coefficients 2880 and -1440, and hence Z_e^(5), Z_e^(2,3). Consequently the coefficients *in the bracket* are Z/120, as the lemma states. Corner-log shifts and lower logarithmic degrees cannot produce the required explicit zeta monomials. The non-local contribution is excluded by Lemma nocst for the specified formal representative.

I checked the polynomial identities independently, without importing the supplied z5-factor script or its missing dependencies. I integrated F^-2 in t first:

\[
 \int_0^1\frac{dt}{(1+\alpha tU+\beta(1-t)V)^2}
 =\frac1{(1+\alpha U)(1+\beta V)}.
\]

The square-root numerator of Psi is independent of t. I recovered its coefficients through degree five from its square, divided by the above denominator, converted U=e^u-1, V=e^v-1 using finite differences, and extracted the required four/two mixed coefficients. In Q[alpha,beta,u_0,m] the two differences `120*c_5/L-Z_e^(5)/L` and `120*c_23/L-Z_e^(2,3)/L` are identically zero. All six degree-five Beta coefficient cases also pass exact assertions.

The generic coefficient identity is algebraic in the edge data. The continuation argument at :123 supplies the claimed domain for the whole expression; no unjustified power-series evaluation at r≥2 is used.

## 6. Graded algebra and determinants

**Verdict: CORRECT.** `sec_g9_hex5.tex:79–87,197`.

The revised object `W=Q direct-sum (direct-sum_{ell>=1} W_ell)` is closed under shuffle multiplication; the degree-limited part alone need not be. Every shuffle ends in 1 when both inputs do, and degrees add. The map w→H_w is a homomorphism. Only spanning, not independence of evaluated functions, is required. The words in the lemma are the words ending in 1, as fixed at :9.

I built the five matrices independently by enumerating shuffle positions and weighted exponent tuples in the listed thirteen generators. Their sizes are 1,2,4,8,16, and my determinants are

`1, 2, -6, -96, 5760`.

At weight two my row order is (01,11). Interchanging these rows gives -2. Thus the stated sequence `1,-2,-6,-96,5760` is a valid signed sequence in a compatible ordering, and the invariant absolute determinants and invertibility are verified exactly. The precise enumeration used by the referenced `g9_hex5_s2_atoms.py` is not supplied here, so I do not claim to have rerun that script's specific ordering. This has no bearing on the spanning argument or on the correctness of the listed magnitudes.

## 7. Status of the five round-one must-fix points

**Verdict: CORRECT AFTER REPAIR; the main repairs landed, but two need the completions above.** The five points are `REFEREE_CODEX_ROUND1.md:25–37`.

| Round-one point | Current location | Status |
|---|---|---|
| 1. Factor 120 | `sec_g9_hex5.tex:134` | **Repaired.** Independent exact identities pass. |
| 2. W_5 defined twice | `sec_g9_oct5.tex:38,62–76` | **Renaming repaired; definition still incomplete.** Add the empty-word case and correct the weight-k / weight-five sentence and quoted proposition as in item 4. |
| 3. Divergent endpoint evaluation | `sec_g9_oct5.tex:46,50–59` | **Substantive issue repaired.** The limit is zero and the remainder bound is justified. Change “diverges” to “may diverge” as in item 2. |
| 4. Coincident letters | `sec_g9_oct5.tex:27` | **Implemented rule now stated and mathematically consistent.** Its identification with the regularised integral still lacks the short proof/unit-tangent specification in item 3. |
| 5. Arbitrary classical arguments | `sec_g9_oct5.tex:85–100,132` | **Fully repaired.** The argument works for arbitrary complex arguments and external complex coefficients, in either symbol order, modulo all proper shuffles. |

The other conditional “must fix if intended unconditionally” wording issue in round-one finding 8 is substantially addressed by the presentation/minimality caveats at `sec_g9_hex5.tex:159,180,184`. The separate new overstatement at :24 is discussed next.

## 8. Other errors in the revised text

**Verdict: WRONG for the independence assertion at `sec_g9_hex5.tex:24`; CORRECT AFTER REPAIR for the intended formal interpretation.** This is an explanatory overstatement, not a counterexample to Theorem 26.

The paragraph says that independence of the zeta constants *and the values of the thirteen functions* is open. At r=2, x=-1, its proposed collection already satisfies

\[
 H_{00001}(-1)=\Li_5(-1)=-\frac{15}{16}\zeta(5).
\]

This follows immediately by separating even and odd terms of the absolutely convergent zeta series, and is even printed in this section at :178. Hence this collection is **not** linearly independent, over Q or over any field containing Q. At r=1 all thirteen HPL values are zero, providing another elementary counterexample. For uniqueness of a numerical monomial coefficient, linear independence of just the individual generators would in any event be insufficient: the relevant monomials or reduced constants must be specified.

Exact replacement for the final sentence of :24:

> These are coefficients of the specified formal representative, not intrinsic coefficients of the numerical value. After specialising r and reducing the harmonic-polylogarithm values, coefficients depend on the chosen presentation of constants. Any uniqueness claim requires linear independence of the resulting reduced constants over a specified coefficient field; no such independence is proved here.

This is consistent with the formal theorem, the r=2 discussion, and the gap note at :139. I have no pre-fix source copy in this folder, so I distinguish an error in the revised text from a proven claim about exactly which edit introduced it.

The other defects exposed in this pass are the empty-word omission and the weight-five/weight-k sentence in item 4, and the overstatement about divergence in item 2. I found no further substantive error in the revised `sec_g9_order5.tex`. Its added sheared-parallelogram common-denominator equality at :176 is arithmetically correct: 714025=25*28561, and both radical coefficients reduce as printed.

## Additional verification and scope

The isolated partial rerun of supplied P1 generated 252 rows with 1,302 nonzero coefficients. I omitted only its subsequent external-table comparison. From those rows my own regularised-symbol recursion, shuffle-position enumerator, top-weight assembly and Lie pairing independently reproduce every displayed block and the total

`(-1485+1035*sqrt(2), 900-540*sqrt(2))`.

A separate rerun of supplied P2 reproduces the same totals in both tensor orders, as well as the pair of C. This distinguishes independent assembly/pairing from a supplied-code rerun; the row generator and its Q(sqrt(2)) arithmetic were reused. The nonzero total and non-proportionality to (84,-38sqrt(2)) are unchanged. Products such as the missing pure Y^5 term have zero pairings.

Runnable checks, using the requested interpreter:

```sh
'<project>/.venv/bin/python' -B referee_g9b_checks.py
'<project>/.venv/bin/python' -B referee_g9b_checks.py coincident
'<project>/.venv/bin/python' -B referee_g9b_checks.py oct-rerun
```

The first command imports no supplied modules. The second compares the independently implemented regularised recursion with the supplied symbol module. The third writes fresh row data only under `referee_g9b_work/`, performs the independent assembly just described, and then reruns supplied P2. Assertions fail on mismatches.

## What I did not check

- I did not rerun the complete hexagon S1–S6 chain, the 307-line table fold, the RUL1 certificates, or their external data comparisons. The supplied factor/determinant wrappers reference source/data absent from this folder; that absence is **not** a finding. Items 5–6 were checked by new computations instead.
- I did not independently regenerate the entire octagon geometry/integrand by a second new engine, compare P1 with the external pickle, or recompute M_4's calibration. P1 was rerun as supplied, stopping before its external comparison.
- I did not verify the full classical-plus-C symbol decomposition, rank 32/32 claim, numerical perimeter derivatives, historical precision claims, PSLQ/LLL searches, or numerical closures for other polygons.
- I did not prove any period conjecture or numerical independence/minimality statement. Complex scalar coefficients in item 1 are external coefficients of formal tensors; that statement is not an assertion about arbitrary representations of a numerical complex value.
- I did not rebuild the PDF or audit external sections g4. The review concerns the supplied g9 mathematics and the targeted fix pass.
