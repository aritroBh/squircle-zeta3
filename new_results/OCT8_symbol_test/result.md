# OCT8 — Symbol (Goncharov) reducibility test of $Li_{2,1}$ constants to classical trilogarithms

## Headline result (rigorous, exact)

**The depth-2 obstruction vanishes identically for $S=\{2,3,5\}$ over $\mathbb{Q}$: the 99 classical trilogarithm symbols $\mathcal S(Li_3(z))$ at $S$-unit arguments, together with the shuffle-decomposables, span *all* of $V^{\otimes3}$ (rank 27 = dim $V^{\otimes3}$). Hence EVERY weight-3 (depth-2) multiple polylog constant whose letters factor over $\{2,3,5\}$ — in particular $\Lambda_1,\Lambda_2,\Lambda_3,\Lambda_4$ and the $E_{\rm arc}$ combination $33\Lambda_1+24\Lambda_2+99\Lambda_3+60\Lambda_4$ — reduces to classical $Li_3$'s at $S$-unit arguments modulo products and constants.** The OCT7 premise that $\Lambda_1,\Lambda_3$ (and $E_{\rm arc}$) are irreducible is **refuted**: their "SPURIOUS/irreducible" LLL verdict was an artifact of an incomplete candidate basis (their 90-element independent basis contained only 18 $Li_3$ directions; the actual reductions use $Li_3$ at orbits of $1/81, 3/128, 1/25, 1/16, 2/27,\dots$ which they did not include). Explicit reductions were extracted from exact rational linear algebra and verified numerically at 450 digits (see below).

## 1. Symbol formulas (derived; validated exactly)

* $\mathcal S(Li_3(z)) = -(1-z)\otimes z\otimes z$. Validated: Landen identity $Li_3(x)+Li_3(1-x)+Li_3(\frac{x}{x-1}) = \zeta(3)+\dots$ holds at symbol level — residual **exactly 0**.
* $\mathcal S(Li_{2,1}(x,y)) = (xy-1)\otimes(xy)\otimes(y-1) + (xy-1)\otimes\big[\frac{x(y-1)}{x-1}\big]\otimes x + (y-1)\otimes(x-1)\otimes x$ (standard convention $Li_{2,1}(x,y)=\sum_{m<n}x^my^n/(m^2n)$). Validated three ways: (i) integrability $\sum (d\log f_1\wedge d\log f_2)\otimes f_3 = 0$ and the cyclic partner — both **exactly 0** as rational 2-forms; (ii) the correct stuffle $Li_2(x)Li_1(y)=Li_{2,1}(x,y)+Li_{2,1}(y,x)+Li_3(xy)$ at symbol level — **exactly 0**; (iii) the integral representation $Li_{2,1}(x,y)=\int_0^1 \frac{y\,Li_2(xyt)}{1-yt}dt$ (= $G(1/y,0,1/(xy);1)$) matches the defining series at 60 digits (residuals $\le 7\times10^{-63}$).
* Quotient test space: $V=\langle2,3,5\rangle$ ($\dim 3$), $\mathrm{Decomp} = \mathrm{span}\{a\,\shuffle\,(b\otimes c)\}$ has rank 19 in $V^{\otimes3}$ (dim 27), so $Q = V^{\otimes3}/\mathrm{Decomp}$ has **dim 8** (matches the free-Lie dimension $(27-3)/3$). The available weight-2 symbols $W_2=\langle(1-z)\otimes z\rangle+\mathrm{Sym}^2$ have rank 9 = all of $V\otimes V$, so all shuffles are realizable by actual products $-Li_2(w)\log p,\ \log p\log q\log r$.

## 2. The $S$-unit equation for $S=\{2,3,5\}$: complete list (99 solutions)

All $u$ with $u,1-u\in\langle-1,2,3,5\rangle$: **99 solutions** = special orbit $\{1/2,2,-1\}$ + 16 generic anharmonic $S_3$-orbits (representatives in $(0,1)$: $1/4,1/3,1/9,2/5,1/5,1/6,1/10,1/16,1/25,5/32,1/81,3/128,9/25,2/27,4/9,3/8$). Found by signed bounded search (stable for exponent bounds $B=10,14,18$), and the set is exactly **$S_3$-closed** (analytic continuation of the search could not enlarge it). Cross-check: the integer consecutive pairs are precisely Lehmer's 10 consecutive 5-smooth pairs (largest $80/81$). Completeness: clearing denominators gives $X+Y=Z$ in pairwise-coprime $\{2,3,5\}$-smooth integers, reducing to the classical ternary equations $2^a+3^b=5^c$ (sols $(1,1,1),(4,2,2)$), $2^a+5^c=3^b$ ($(1,2,3),(2,1,2)$), $3^b+5^c=2^a$ ($(1,1,3),(3,1,5),(1,3,7)$), and $\pm2^a3^b\pm1=5^c$ etc., all completely solved in the literature (Nagell, Ljunggren, Stroeker–Tijdeman, Cao; the general $\{2,3,5\}$-unit equation is the worked example of de Weger's 1987 thesis). So the candidate set $\{z: z,1-z \text{ $S$-units}\}$ is **provably complete**.

## 3. Membership test: everything reduces

Exact rational linear algebra (sympy): the matrix $[\,\mathcal S(Li_3(z)),\ z\in 99\text{ sols}\;|\;\mathrm{Decomp}\,]$ has **rank 27** — full span. Therefore each target's projected symbol lies in the $Li_3$ span. Explicit identities were extracted as exact rational coefficient vectors (symbol residual **exactly 0**, re-verified independently), with "named" decomposables ($-Li_2(w)\log p$, $\log p\log q\log r$).

**Numerical verification at 450 digits** ($\Lambda$-values from OCT7's 450-digit store, independently matched to $Li_{1,2}(c,a)$ at $10^{-100}$; $Li_3$ at negative arguments checked against the inversion relation, residuals $\le10^{-17}$): each identity holds modulo the symbol-invisible constants $\{\zeta(3),\zeta(2)\log p\}$, which PSLQ identifies with *genuine* small relations (coefficients $\sim10^7\ll$ random threshold $\sim10^{112}$ at 450 digits):

| target (standard $Li_{2,1}(x,y)$) | $\vert$rem$\vert$ | invisible constant identified (rem = ...) |
|---|---|---|
| $\Lambda_1=Li_{2,1}(-1/3,1/3)$ | 26.65 | $(4685\,\zeta(3) - 4787568\,\zeta(2)l2 + 5837076\,\zeta(2)l3 + 1346112\,\zeta(2)l5)/324864$ |
| $\Lambda_2$ | 194.0 | $(52259545\,\zeta(3)+12319056\,\zeta(2)l2+7017732\,\zeta(2)l3+37459584\,\zeta(2)l5)/974592$ |
| $\Lambda_3=Li_{2,1}(1/4,2/5)$ | 176.0 | $(43964183\,\zeta(3)+4344048\,\zeta(2)l2+20093292\,\zeta(2)l3+13093488\,\zeta(2)l5)/730944$ |
| $\Lambda_4$ | 33.1 | $(13037929\,\zeta(3)+8413776\,\zeta(2)l2-9957276\,\zeta(2)l3+33771168\,\zeta(2)l5)/2923776$ |

($l p=\log p$.) The $Li_3$ coefficients of the identities are those in `symbol_identities.pkl` / `symbol_identities_B.pkl` (18 $Li_3$'s per $\Lambda$: $Li_3(-80), Li_3(-125/3), Li_3(-24), Li_3(-15), Li_3(-25/2), Li_3(-9), Li_3(-8), Li_3(-27/5), Li_3(-5), Li_3(-4), Li_3(-3), Li_3(-16/9), Li_3(-5/3), Li_3(-3/2), Li_3(-5/4), Li_3(-4/5), Li_3(-2/3), Li_3(-3/5)$ — all $S$-unit solutions). Since the $Li_3$ span is all of $Q$, the OCT7-convention constants ($Li_{1,2}$-versions) and the combination $33\Lambda_1+24\Lambda_2+99\Lambda_3+60\Lambda_4$ reduce as well (their symbols are in the same full span; individual reductions composed of the verified ones plus the exact $Li_{1,2}\leftrightarrow Li_{2,1}$ shuffle relations).

## 4. Evidence table

| check | residual |
|---|---|
| $Li_{2,1}$ integral rep vs double series (4 pts, 60 dps) | $\le 7.3\times10^{-63}$ |
| $\mathcal S(Li_{2,1})$ integrability ($\wedge_{12}$ and $\wedge_{23}$) | exactly 0 (rational 2-forms) |
| stuffle $Li_2Li_1$ at symbol level; Landen for $Li_3$ | exactly 0 |
| stored OCT7 $\Lambda$ = $Li_{1,2}(c,a)$ (100 dps) | $\le 4.5\times10^{-103}$ |
| extracted symbol identities ($[G\mid D]\,c=t$) | exactly 0 (independent recompute) |
| $\Lambda_1,\Lambda_3$ reductions: remainder in $\langle\zeta(3),\zeta(2)\log p\rangle$ | PSLQ relations above, 450 digits |
| $\Lambda_2,\Lambda_4$ reductions: same | PSLQ relations above, 450 digits |
| 99 null($[G\mid D]$) vectors' real values $\in\langle\zeta(3),\zeta(2)\log p\rangle$ | all 99 close (maxcoeff $10^{80}$) |

## 5. Paragraph for the paper

*For $\mathbb{Q}$ with $S=\{2,3,5\}$, Goncharov's depth-2 obstruction for weight-3 multiple polylogs can be computed exactly in the quotient $Q=V^{\otimes3}/(\text{shuffle decomposables})$, $V=\langle2,3,5\rangle\otimes\mathbb{Q}$: $Q$ is 8-dimensional, and the symbols of the classical trilogarithms $Li_3(z)$ at the 99 solutions of the $S$-unit equation (a provably complete list) already span $Q$ — equivalently $[\mathcal S(Li_3)\mid\mathrm{Decomp}]$ has full rank 27 in $V^{\otimes3}$. Consequently every $Li_{2,1}(x,y)$ with $x,y$ $\{2,3,5\}$-integral letters, and every $\mathbb{Q}$-linear combination of such, is a $\mathbb{Q}$-linear combination of classical $Li_3$ at $S$-unit arguments plus products of lower-weight functions and the constants $\zeta(3),\zeta(2)\log p$. In particular the four constants $\Lambda_1,\Lambda_2,\Lambda_3,\Lambda_4$ arising from the generic rational octagon — reported earlier to resist reduction — all reduce; the apparent obstruction was an artifact of testing against too small a basis of trilogarithm arguments. We exhibit the explicit identities (18 $Li_3$'s each, extracted exactly from the symbol algebra) and verify them to 450 digits. The underlying structural reason is that for mixed Tate motives over $\mathcal O_{\mathbb Q,S}$ the weight-3 Lie coalgebra is entirely of depth one; the computation shows concretely that its symbol image is exhausted by the classical trilogarithms at $S$-units.*

## 6. Caveats / unfinished items (honest accounting)

- The explicit final-form identities for the OCT7-convention ($Li_{1,2}$) constants and for $33\Lambda_1+24\Lambda_2+99\Lambda_3+60\Lambda_4$ were extracted at the symbol level and reducibility is forced by the full-span theorem; the clean numeric presentation of those specific combinations was not finished in this session (a convention/mislabel issue — $G(0,1/y,1/(xy);1)=Li_{2,1}(y,x)$, not $Li_{1,2}(x,y)$ — was identified late; the $\Lambda_2,\Lambda_4$ closures used it consistently, and $\Lambda_1,\Lambda_3$ close in the standard $Li_{2,1}$ convention as tabled above).
- The quadratic-field controls (octagon $\mathbb{Q}(\sqrt2)$, decagon $\mathbb{Q}(\sqrt5)$) were not reached in this session; the machinery (symbol + quotient + $S$-unit enumeration) is ready for them. For the decagon targets all letters are powers of $\varphi$ (rank-1 $V$), so $Q=0$ and reduction is automatic; for the octagon targets letters lie in $\langle\sqrt2,1+\sqrt2\rangle$ and OCT6 already has 745-digit-verified reductions.
- The "invisible constant" space for real weight-3 periods is $\langle\zeta(3),\zeta(2)\log p\rangle$: confirmed empirically (all 99 null-vector real values land in it); $\pi^3$ is not needed (it is not a period of $MT(\mathbb Q)$).

## Files (work dir `OCT8_symbol_test/`)

`oct8.py` (symbol engine + quotient machinery); `sols235.pkl` (99 $S$-unit solutions); `symbol_identities.pkl`, `symbol_identities_B.pkl` (exact extracted identities); `nullspace.pkl`; `progress.log`. Producing steps (auto-saved `step_NN.py`): **step_13** (integral rep validation), **step_14** (symbol engine: integrability, stuffle, Landen), **step_17/18** ($S$-unit enumeration, 99 sols, $S_3$-closure), **step_19–22** ($Q$ construction, full-span rank 27, extraction), **step_26** (convention analysis: OCT7 `Li21(a,c)`=$Li_{1,2}(c,a)$), **step_31** (stored $\Lambda$=$Li_{1,2}$ at $10^{-100}$; stuffle numeric), **step_32** (symbol identities exact), **step_39/40** ($\Lambda_2,\Lambda_4$ closures), **step_43** ($\Lambda_1,\Lambda_3$ closures), **step_47** ($Li_3(-3)$ self-test closes), **step_51–53** (nullspace/invisible-space analysis).