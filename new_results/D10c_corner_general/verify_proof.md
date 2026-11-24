# D10c verification, key = `proof` (Part 1 (a)–(f) and Part 2 of `result.md`)

All scripts and outputs are in `<scratch>/verify_proof_D10c/` (`v1_algebra.py` … `v8_simplified_logpart.py`, each with a `.out`). Nothing was imported from the D10c scripts except where noted as a comparison target (`step_19.M1_expanded`, `step_22.reduced_closed_form`). Nothing under the project was modified except this file.

## Verdict

| item | status |
|---|---|
| Part 1 closed form for $M_1(\alpha)$, $0<\alpha<\pi$ | **Correct.** Every algebraic step (a), (b), (c), (e) re-derived independently and exact. |
| Part 1 proof as written | **Incomplete.** (d)/(f) use complex $\operatorname{Li}_2$ evaluated *on its branch cut* along the whole path for 10 of the 19 off-diagonal pairs. Reality and correctness of the $-i\pi$ convention were only checked numerically (open question 5). Two analytic steps are also unstated. All gaps are closed below by a proof that uses only the real dilogarithm. |
| (5) "the $\zeta(2)$ coefficient is exactly $-1$" | **True only as bookkeeping in the D10c representation, and misleading as stated.** It is not an invariant of $M_1$. In the real form it is $+1$. In D10c's own reduced form it is $-1+c$. |
| Part 2: all coefficients in $\mathbb Q+\mathbb Qc$ | **Conclusion true** (indeed in $\mathbb Z[\tfrac12]+\mathbb Z[\tfrac12]c$). The stated mechanism ("reductions only permute coefficients") is **false as worded**, but the corrected Q-linearity argument works. |
| Part 2 side claim "maxcoeff ≤ 2 … in fact only $\{\pm\frac12,\pm1,\pm2\}\cup\{\pm\frac c2,\pm c,\pm2c\}$ occur" | **False** for the Li₂ and ζ(2) coefficients of D10c's own reduced form. It holds only for the unmerged log-product coefficients. |

## Independent checks performed

| script | check | result |
|---|---|---|
| v1 | (0) $\partial_u$ of the corner integrand at $u=0$, straight from the definition, minus $\frac1g\big[A_0(\frac{1-c_ag}{S_0}\log g-\log(1+g))-\log g\big]$ | exact 0 |
| v1 | (a) $S_0(g(t))=P^2/4(t-c_a)^2$; $t-g=P/2(t-c_a)$; $1+g$; $dg/dt$; $P=(t-c_a)^2+4c^2s^2$; solving $(t-g)^2=S_0$ returns exactly $g(t)$; $g(1)=0$; $g(1+2s)=1$; $S_0-(g-c_a)^2=4c^2s^2$; the two further identities of (a) (v9) | all exact 0 |
| v1 | (b) coefficients $X,Y$ of $\log g,\log(1+g)$, built directly from the $g$-integrand, minus $A^*,B^*$; the numerator $N$ | exact 0; also $X$ reduces to $-[(t+1)^3-4c^2(t+1)-8c^4(t-1)]/[(t+1)(t-c_a)(t-r_1)(t-r_2)]$ |
| v1 | (c) 10 residues computed from $X,Y$ (not from the claimed $A^*$), plus both full partial-fraction identities | exact 0 |
| v1 | pole gaps factored: $1-c_a=2(1-c^2)$, $1-r_1=2(1-c)$, $c_a-r_1=2c(c-1)$, $-1-r_2=2c$, … | ordering $r_2<-1<c_a<r_1<1$ for $c\in(0,1)$ |
| v2 | real pair primitives, 3 cases, positive symbols, **no** `force=True` | exact 0 |
| v2 | D10c complex primitive vs real primitive: $I-J=\tfrac12\log^2(\sigma-\rho)-i\pi\log(\sigma-\rho)-\pi^2/3$, independent of $t$; mpmath branch checked ($\operatorname{polylog}(2,y)=\frac{\pi^2}3-\frac12\log^2y-\operatorname{Li}_2(1/y)-i\pi\log y$ at $y=1.5,2,5$) | certificate for the $-i\pi$ convention |
| v3 | $dF/dt-(A^*\log g+B^*\log(1+g))$ with logs split into logs of positive factors | exact 0 |
| v3 | real formula vs **own** $g$-quadrature (dps 60), α = 0.01, 0.2, π/2, 2.0, 3.0, π−0.01 | 56.5–60.9 digits |
| v3 | real formula vs Richardson central difference of $M(u;\alpha)$ itself ($h=10^{-9}$). This tests "definition → integrand" independently. | 37.7–44.2 digits (the $O(h^4)$ limit) |
| v3 | real formula vs D10c `M1_expanded` (Re) | 55.5–60.1 digits; D10c $|\mathrm{Im}|\le 4\times10^{-60}$ |
| v4 | explicit 27-term real formula vs quadrature (dps 80), α = 0.01, 0.3, π/2, 2π/3, 2.5, 3.1, π−10⁻³; also 10⁻⁶ and π−10⁻⁶ | 76.1–81.0 digits; 67.7 and 78.7 at the extremes |
| v6 | **exact** term-by-term identity: explicit formula = pair-sum $F(t_1)-F(1)$. Li₂ arguments matched exactly via $c=\frac{1-m^2}{1+m^2}$, $s=\frac{2m}{1+m^2}$. Log part compared in a letter basis modulo $s^2=(1-c)(1+c)$ and $(1+s-c)(1+s+c)=2s(1+s)$. | exact 0; $\Lambda_1=0$; ζ(2) coefficient $+1$ |
| v8 | simplified log part (the one stated in the theorem below) vs quadrature, 6 angles | 57.1–60.8 digits (dps 60) |
| v5, v5b | D10c reduced form (`step_22`) rerun with their own recognizer at dps 220 | Li₂ coefficients $\pm4c$, $1+c$, $-1+c$, $-1-c$, $2-2c$, $2+2c$, $-2-2c$ occur; ζ(2) coefficient $-1+c$ at 5π/6 and π/12, and $-2$ at α = 1.3 |
| v7 | minimal polynomials of $\cos(\alpha/2)$ at 5π/6, π/12, 3π/4, π/6 | degrees 4, 8, 4, 4, matching `result.md` |
| v7 | sanity limits | $M_1(0^+)=\log^22$ (by direct integration at α = 0), $M_1(\pi-\varepsilon)\sim\varepsilon^2/8$ |

## Gaps (in order of the task list)

**G1 (proof gap, minor; Part 1 setup).** `result.md` never justifies $M_1=\partial_uM|_{u=0}=\int_0^1\partial_u(\cdot)\,dg$. It is true: the dominating bound $|\partial_u h|\le 8(1+|\log g|)g^{-1/4}$ for $|u|\le\frac14$ is given in Step 0 below. This step was also confirmed numerically by the finite difference of $M(u)$ itself (v3).

**G2 (proof gap, minor; item 1).** The change of variables is stated as "path $t:1\to t_1$" without proof. Missing pieces:
- Monotonicity: $dt/dg=1+(g-c_a)/\sqrt{S_0}>0$ because $S_0-(g-c_a)^2=\sin^2\alpha>0$.
- Inverse: $(t-g)^2=S_0$ gives the unique $g=(t^2-1)/2(t-c_a)$, with $t>1>c_a$.
- Sign of the root: $\sqrt{S_0}=+P/2(t-c_a)$ with $P=(t-c_a)^2+\sin^2\alpha>0$.
- Endpoint: $t_1=1+\sqrt{2-2c_a}=1+2s$.

All are verified (v1 a1–a9) and written out in Step 1. Item (1) itself is correct.

**G3 (main proof gap; item 3).** (d) is labelled "PROVED (formal-log)" and the reality of (f) is NUMERICAL. The problem is concrete.
- For every pair with $\rho<\sigma$, $z=(t-\rho)/(\sigma-\rho)>1$ along the **entire** path, since $t>1\ge\sigma$. So $-\operatorname{Li}_2(z)$ is evaluated on its cut $[1,\infty)$.
- The pair identity $\log\frac{\sigma-t}{\sigma-\rho}+\log(\rho-\sigma)=\log(t-\sigma)$ takes logs of negative numbers on both sides. It holds only modulo $2\pi i$ (here: with an $i\pi$ shift fixed by which side of the cut is used).
- `step_14 [1]` proves it with `force=True`, which proves nothing about branches.

Fix: use the real primitives of Step 3. They need no branch at all, since every Li₂ argument lies in $(-\infty,1]$. They also give an exact certificate for the D10c complex form. With mpmath's boundary value $\operatorname{Li}_2(y)=\operatorname{Li}_2(y-i0)$ for $y>1$ (checked), each pair satisfies $I_{\rm D10c}=J+\tfrac12\log^2(\sigma-\rho)-i\pi\log(\sigma-\rho)-\pi^2/3$, a $t$-independent constant. So $I(t_1)-I(1)=J(t_1)-J(1)$ pair by pair. This also holds at $z=1$, since Li₂ is continuous there.

So the D10c formula **is** exactly real and exactly equal to $M_1$, and open question 5 is closed. Two corrections to the wording:
- The $-i\pi$ pieces cancel per pair, not just "globally".
- The $-i\pi$ choice is tied to the lower-side convention: a CAS returning $\operatorname{Li}_2(y+i0)$ would need $+i\pi$. "Take Re" is convention-free.

For the paper, use the real form only.

**G4 (presentation; item 4).** The endpoint cancellation is correct, and stronger than stated. The coefficient of $\log(t-1)$ in the $\rho=1$ row is
$$\Lambda_1=\log\frac{2(1-c_a)}{(1-r_1)(1-r_2)}=\log1=0$$
**identically in $t$**, not only in the limit. So the $\rho=1$ row is exactly $F_1(t)=-\sum_\sigma W(1,\sigma)\operatorname{Li}_2\!\big(-\frac{t-1}{1-\sigma}\big)$ with $F_1(1)=0$, and its upper-endpoint log terms vanish too.

As printed, however, the boxed formula contains $\log(1-\rho)$ at $\rho=1$, i.e. $\log0$, in its last sum and in the $I(1^+;1,\sigma)$ terms. "Regularized to 0" must be replaced by the exact statement above: the $\rho=1$ row is written separately, with no log terms.

Also needed and unstated: continuity of $F$ at $t=1$. For $\rho\ne1$ the only pairs reaching $t=\sigma$ are $(\rho,1)$, where $\operatorname{Li}_2$ is continuous at 1. The FTC for the improper integral is applied in Step 4.

**G5 (misstatement; item 5).** "$\operatorname{Li}_2(1)=\zeta(2)$ occurs with total coefficient $-1$ (proved)" is a correct *count* of the literal $\operatorname{Li}_2(1)$ terms in the D10c representation. It is not a property of $M_1$:
- In that representation, $\operatorname{Re}\operatorname{Li}_2(y)=2\zeta(2)-\frac12\log^2y-\operatorname{Li}_2(1/y)$ $(y>1)$ hides further $\zeta(2)$'s.
- In the real representation below the coefficient is $+1$. Converting each $\operatorname{Re}\operatorname{Li}_2(y)$ with $y>1$ to $\operatorname{Li}_2(1/y)$ adds $-2\zeta(2)W$ at the upper endpoint and $+2\zeta(2)W$ at the lower endpoint of every pair $\rho<\sigma$. These cancel, except for the four pairs $(\rho,1)$, whose lower value $y_0=1$ is already the literal $\operatorname{Li}_2(1)$. The net change is $-2\zeta(2)\sum_{\rho\ne1}a_\rho=+2\zeta(2)$, so $-1\to+1$.
- A single reflection inside the real form changes it again. $\operatorname{Li}_2\!\frac{s}{1+s}+\operatorname{Li}_2\!\frac1{1+s}=\zeta(2)-\log\frac s{1+s}\log\frac1{1+s}$, and both terms occur, which turns $+1$ into $+2$.
- D10c's own canonical reduction (`step_22`, args in $(0,\frac12]$) gives $-1+c$ at 5π/6 and π/12 and $-2$ at α = 1.3. So within that scheme the coefficient even depends on α.

What *is* proved and invariant: $\sum_{\rho\ne1}a_\rho=-1$ is independent of $c$, because $a(r_1)+a(r_2)=0$. A paper should say "in the representation of Theorem X the constant term is $\zeta(2)$" and nothing stronger.

**G6 (Part 2 mechanism; item 6).** The statement "reduction … only *permutes* coefficients … appends ζ(2)- and log-product terms carrying the *same* single coefficient" is false as worded, on several counts:
- Inversion for $x>1$ appends $-2W\zeta(2)$ and $+\frac W2\log^2x$.
- Landen appends $\frac W2\log^2(1-x)$.
- Reflection appends $W\zeta(2)$ and $-W\log x\log(1-x)$, and flips the Li₂ sign.
- Merging equal arguments adds coefficients, e.g. $1-c$, $-(1+c)$, $\pm4c$, $2\pm2c$.
- Expanding $\log$ of products of letters multiplies by integer exponents.

Also, "each pair gets exactly one universal residue" is wrong for the $\sigma=c_a$ column, where $W(\rho,c_a)=-(a_\rho+b_\rho)$ is a sum of two residues.

None of this breaks the conclusion. Every one of these operations is $\mathbb Z[\frac12]$-linear on the coefficient vector, and the unreduced coefficients lie in $\mathbb Z+\mathbb Zc$. So every coefficient stays in $\mathbb Z[\frac12]+\mathbb Z[\frac12]c$. The real reason is linearity in $W$: no two coefficients are ever multiplied. The corrected proof is in the Corollary below. Note too that "$\mathbb Q\oplus\mathbb Qc$" is a direct sum only when $c\notin\mathbb Q$; for example α = 2π/3 gives $c=\frac12$. And the property is about these representations, since individual coefficients are not invariants (cf. G5).

**G7 (Part 2 side claim false).** Rerunning `step_22.reduced_closed_form` with their `rec_in_Qc` at dps 220 (v5b):
- **5π/6:** Li₂ coefficients $-1+c$ (×2), $1+c$ (×2), $2-2c$, $-4c$, $-2-2c$, $4c$, $-2c$, $2$; ζ(2) coefficient $-1+c$.
- **π/12:** includes $-4c$, $2+2c$, $4c$, $-1-c$.

So "maxcoeff ≤ 2" fails (±4c), and "only $\{\pm\frac12,\pm1,\pm2\}\cup\{\pm\frac c2,\pm c,\pm2c\}$ occur" fails for Li₂ and ζ(2). It is true only for the 12 *unmerged* log-product coefficients, which are trivially $\pm W,\pm W/2$. All values are still in $\mathbb Q+\mathbb Qc$, consistent with G6.

**G8 (cosmetic).**
- The argument list in (f) has loose signs. Only $-s(1+s)/c^2$ and $-s^2/c^2$ occur, never $+$. The pairings are $+\frac{s(1+s)}{c(1-c)}$ and $-\frac{s(1+s)}{c(1+c)}$, and $\mp\frac{1+s\mp c}{c}$.
- The count of 31 Li₂ terms is right: 19 upper plus 12 lower non-ζ(2) terms. They merge into 27 distinct arguments.
- (b) should say that the $(t-1)$ factor of $N$ cancels the $(t-1)$ in $D$, so $A^*$ has no pole at 1. That is what makes $a_1=0$.

---

## Corrected proof (paper-ready; real dilogarithms only)

*(Written to drop into the corner-asymptotics subsection of "Investigating ζ(3)"; numbering to be adjusted.)*

Throughout, $0<\alpha<\pi$, $c=\cos\frac\alpha2$, $s=\sin\frac\alpha2$ (so $c,s\in(0,1)$, $s^2=1-c^2$), and $c_a=\cos\alpha=2c^2-1$. By $\operatorname{Li}_2$ we always mean the **real** dilogarithm $\operatorname{Li}_2(x)=-\int_0^x\frac{\log(1-u)}u\,du$ on $(-\infty,1]$. It is real-analytic on $(-\infty,1)$ with $\operatorname{Li}_2'(x)=-\log(1-x)/x$, and continuous at $1$ with $\operatorname{Li}_2(1)=\zeta(2)$. We never evaluate it anywhere else.

We want the first Taylor coefficient $M_1(\alpha)$ of
$$M(u;\alpha)=\int_0^1\Big[(1+g)^{-1-u}\sqrt{1+g^{2-2u}-2c_a\,g^{1-u}}-1\Big]g^{u-1}\,dg .$$

### Step 0: one integral for $M_1$

Let's start by differentiating under the integral sign, which we have to justify. Call the integrand $h(u,g)$ and put $Q=1+g^{2-2u}-2c_ag^{1-u}=(g^{1-u}-c_a)^2+\sin^2\alpha\ge\sin^2\alpha>0$. So $h$ is smooth in $u$.

For $|u|\le\frac14$ and $0<g\le1$ we have $Q\le4$ and $|g^{1-u}-c_a|\le\sqrt Q$. We also have $|\sqrt Q-1|\le|Q-1|\le3g^{1-u}$ and $|(1+g)^{-1-u}-1|\le\frac54g$, which give $|h|\le5g^{-1/4}$. Since
$$\partial_uh=g^{u-1}(1+g)^{-1-u}\Big[-\log(1+g)\sqrt Q-\log g\,\frac{g^{1-u}(g^{1-u}-c_a)}{\sqrt Q}\Big]+\log g\cdot h ,$$
we get $|\partial_uh|\le 8(1+|\log g|)\,g^{-1/4}$, which is integrable on $(0,1]$. By dominated convergence $M$ is $C^1$ near $u=0$, and
$$M_1(\alpha)=\int_0^1\frac1g\Big[A_0\Big(\frac{1-c_ag}{S_0}\log g-\log(1+g)\Big)-\log g\Big]dg,\qquad S_0=1-2c_ag+g^2,\quad A_0=\frac{\sqrt{S_0}}{1+g}. \tag{0}$$
The integrand behaves like $-\log g$ at $g=0$, so it is absolutely integrable.

### Step 1: Euler's substitution

Put $t=g+\sqrt{S_0}$. Because $S_0-(g-c_a)^2=\sin^2\alpha>0$, we have $dt/dg=1+(g-c_a)/\sqrt{S_0}>0$. So $t$ increases strictly from $t(0)=1$ to $t(1)=1+\sqrt{2-2c_a}=1+2s=:t_1$, and the substitution is a smooth bijection $[0,1]\to[1,t_1]$.

Squaring $(t-g)^2=S_0$ gives $g=\frac{t^2-1}{2(t-c_a)}$, with $t\ge1>c_a$. With $P(t)=t^2-2c_at+1=(t-c_a)^2+\sin^2\alpha>0$, $r_1=2c-1$ and $r_2=-1-2c$, a direct computation gives
$$\sqrt{S_0}=t-g=\frac{P}{2(t-c_a)},\quad 1+g=\frac{(t-r_1)(t-r_2)}{2(t-c_a)},\quad \frac{dg}{dt}=\frac{P}{2(t-c_a)^2},\quad 1-c_ag=\frac{-c_at^2+2t-c_a}{2(t-c_a)}.$$
The five numbers $\mathcal P=\{1,-1,c_a,r_1,r_2\}$ satisfy $r_2<-1<c_a<r_1<1$, since $-1-r_2=2c$, $c_a+1=2c^2$, $r_1-c_a=2c(1-c)$ and $1-r_1=2(1-c)$. So on $(1,t_1]$ every $t-\sigma$ ($\sigma\in\mathcal P$) is positive. Writing $\ell_\sigma=\log(t-\sigma)$,
$$\log g=\ell_1+\ell_{-1}-\ell_{c_a}-\log2,\qquad \log(1+g)=\ell_{r_1}+\ell_{r_2}-\ell_{c_a}-\log2 .$$

### Step 2: the rational coefficients

Substituting into (0),
$$M_1(\alpha)=\int_1^{t_1}\big[A(t)\log g+B(t)\log(1+g)\big]dt,$$
$$A(t)=-\frac{(t+1)^3-4c^2(t+1)-8c^4(t-1)}{(t+1)(t-c_a)(t-r_1)(t-r_2)},\qquad B(t)=-\frac{P(t)^2}{(t^2-1)(t-c_a)(t-r_1)(t-r_2)}.$$
Both are proper fractions with simple poles in $\mathcal P$, and
$$A=\sum_{\rho\in\mathcal P}\frac{a_\rho}{t-\rho},\quad B=\sum_{\rho\in\mathcal P}\frac{b_\rho}{t-\rho},\qquad (a_\rho)=(0,-2,1,-c,c),\quad (b_\rho)=(-1,1,-1,2c,-2c),$$
in the order $\rho=1,-1,c_a,r_1,r_2$. The residues are linear in $c$; that is the whole reason for the coefficient rule below. Hence
$$A\log g+B\log(1+g)=\sum_{\rho,\sigma\in\mathcal P}W(\rho,\sigma)\frac{\ell_\sigma}{t-\rho}+\log2\sum_{\rho}\frac{\omega_\rho}{t-\rho},$$
where $W(\rho,1)=W(\rho,-1)=a_\rho$, $W(\rho,c_a)=-(a_\rho+b_\rho)$, $W(\rho,r_1)=W(\rho,r_2)=b_\rho$, and $\omega_\rho=-(a_\rho+b_\rho)=(1,1,0,-c,c)$.

### Step 3: real primitives

**Lemma.** Let $\rho,\sigma$ be real and $t>\max(\rho,\sigma)$. Then $\frac{d}{dt}J_{\rho\sigma}(t)=\frac{\log(t-\sigma)}{t-\rho}$, where
$$J_{\rho\rho}=\tfrac12\log^2(t-\rho);\qquad J_{\rho\sigma}=-\operatorname{Li}_2\!\Big(-\frac{t-\rho}{\rho-\sigma}\Big)+\log(\rho-\sigma)\log(t-\rho)\ \ (\rho>\sigma);$$
$$J_{\rho\sigma}=\operatorname{Li}_2\!\Big(\frac{\sigma-\rho}{t-\rho}\Big)+\tfrac12\log^2(t-\rho)\ \ (\rho<\sigma).$$
In the second case the argument is negative, in the third it lies in $(0,1)$, and when $\rho<\sigma$ the function $J_{\rho\sigma}$ extends continuously to $t=\sigma$.

*Proof.* For $\rho>\sigma$, set $x=t-\rho$ and $d=\rho-\sigma$. Then $-\frac{d}{dx}\operatorname{Li}_2(-x/d)=\frac{\log(1+x/d)}{x}=\frac{\log(x+d)-\log d}{x}$, and adding $\frac{\log d}{x}$ gives $\frac{\log(t-\sigma)}{t-\rho}$.

For $\rho<\sigma$, set $y=\frac{\sigma-\rho}{t-\rho}\in(0,1)$. Then $\frac{d}{dt}\operatorname{Li}_2(y)=\frac{\log(1-y)}{t-\rho}=\frac{\log(t-\sigma)-\log(t-\rho)}{t-\rho}$, and adding $\frac{d}{dt}\frac12\log^2(t-\rho)$ gives the claim. All logarithms are of positive numbers. $\square$

### Step 4: the primitive and the endpoint $t=1$

For $\rho\neq1$ put $F_\rho(t)=\sum_\sigma W(\rho,\sigma)J_{\rho\sigma}(t)+\omega_\rho\log2\,\log(t-\rho)$. This is continuous on $[1,t_1]$: $\rho<1\le t$, and the only $\sigma$ that $t$ reaches is $\sigma=1>\rho$, the continuous case of the Lemma.

For $\rho=1$, every $\sigma\ne1$ lies below $1$, and $W(1,\cdot)=(0,0,1,-1,-1)$, $\omega_1=1$. So
$$\sum_\sigma W(1,\sigma)J_{1\sigma}+\omega_1\log2\,\log(t-1)=-\sum_\sigma W(1,\sigma)\operatorname{Li}_2\!\Big(-\frac{t-1}{1-\sigma}\Big)+\Lambda_1\log(t-1),$$
with
$$\Lambda_1=\log(1-c_a)-\log(1-r_1)-\log(1-r_2)+\log2=\log\frac{2\cdot2s^2}{2(1-c)\cdot2(1+c)}=0 .$$
So the would-be divergence cancels *identically*, and $F_1(t):=-\sum_\sigma W(1,\sigma)\operatorname{Li}_2\big(-\frac{t-1}{1-\sigma}\big)$ is continuous with $F_1(1)=0$.

Then $F=\sum_\rho F_\rho$ is continuous on $[1,t_1]$ and $F'=A\log g+B\log(1+g)$ on $(1,t_1]$. The integrand has only a logarithmic singularity at $t=1$, so
$$M_1(\alpha)=F(t_1)-F(1).$$

### Step 5: evaluation

We need the values
- $t_1-1=2s$, $t_1+1=2(1+s)$, $t_1-c_a=2s(1+s)$, $t_1-r_1=2(1+s-c)$, $t_1-r_2=2(1+s+c)$;
- $1-c_a=2s^2$, $1-r_1=2(1-c)$, $1-r_2=2(1+c)$;
- the gaps listed in Step 1.

With these, Step 4 becomes explicit:

**Theorem.** For $0<\alpha<\pi$, with $c=\cos\frac\alpha2$, $s=\sin\frac\alpha2$, $X_-=\frac{1+s-c}{1-c}$ and $X_+=\frac{1+s+c}{1+c}$,
$$M_1(\alpha)=\zeta(2)+U(\alpha)+V(c)+\mathcal L(\alpha),$$
where
$$\begin{aligned}
U={}&-\operatorname{Li}_2\!\big(-\tfrac1s\big)+\operatorname{Li}_2\!\big(-\tfrac{s}{1-c}\big)+\operatorname{Li}_2\!\big(-\tfrac{s}{1+c}\big)
-\operatorname{Li}_2\!\big(-\tfrac{1+s}{c}\big)-\operatorname{Li}_2\!\big(-\tfrac{s(1+s)}{c^2}\big)+\operatorname{Li}_2\!\big(-\tfrac{s(1+s)}{c(1+c)}\big)\\
&+c\operatorname{Li}_2\!\big(-\tfrac{1+s-c}{c}\big)+c\operatorname{Li}_2\!\big(-\tfrac{1+s-c}{c(1-c)}\big)-2c\operatorname{Li}_2\!\big(-\tfrac{1+s-c}{2c}\big)\\
&-2\operatorname{Li}_2\!\big(\tfrac1{1+s}\big)+\operatorname{Li}_2\!\big(\tfrac{c^2}{1+s}\big)+\operatorname{Li}_2\!\big(\tfrac{c}{1+s}\big)+\operatorname{Li}_2\!\big(\tfrac{s}{1+s}\big)-\operatorname{Li}_2\!\big(\tfrac{c(1-c)}{s(1+s)}\big)\\
&-c\operatorname{Li}_2\!\big(\tfrac{1-c}{1+s-c}\big)+c\operatorname{Li}_2\!\big(\tfrac{1+c}{1+s+c}\big)+c\operatorname{Li}_2\!\big(\tfrac{c}{1+s+c}\big)+c\operatorname{Li}_2\!\big(\tfrac{c(1+c)}{1+s+c}\big)-2c\operatorname{Li}_2\!\big(\tfrac{2c}{1+s+c}\big),
\end{aligned}$$
$$V=(1-c)\operatorname{Li}_2\!\big(-\tfrac1c\big)+\operatorname{Li}_2\!\big(1-\tfrac1{c^2}\big)-(1+c)\operatorname{Li}_2\!\big(-\tfrac{1-c}{c}\big)+2c\operatorname{Li}_2\!\big(-\tfrac{1-c}{2c}\big)-\operatorname{Li}_2(c^2)-(1+c)\operatorname{Li}_2(c)+(1-c)\operatorname{Li}_2\!\big(\tfrac{c}{1+c}\big)+2c\operatorname{Li}_2\!\big(\tfrac{2c}{1+c}\big),$$
$$\mathcal L=\log(1+s)\log\tfrac{c}{1+s}+\log\tfrac{c}{1+c}\log\tfrac{1+s}{s}+c\log4\,\log X_-+\tfrac c2\log^2X_--\tfrac c2\log X_+\,\log\big((1+c)(1+s+c)\big).$$
Every $\operatorname{Li}_2$ argument lies in $(-\infty,1)$ and every logarithm is of a positive number.

$U$ comes from the 19 upper-endpoint pairs, $V$ from the 12 lower-endpoint pairs other than $(\rho,1)$ (merged into 8 arguments that depend on $c$ only), and $\zeta(2)$ from $\operatorname{Li}_2(1)$ at the four pairs $(\rho,1)$, $\rho\ne1$. Its coefficient is $-\sum_{\rho\ne1}a_\rho=-(-2+1-c+c)=1$.

*Checks.*
- The Theorem equals $F(t_1)-F(1)$ exactly, term by term (v6).
- It agrees with direct quadrature of (0) to 76–81 digits at α = 0.01, 0.3, π/2, 2π/3, 2.5, 3.1, π−10⁻³, and with the finite-difference derivative of $M(u;\alpha)$ itself (v3, v4, v8).
- Limits: $M_1(0^+)=\log^22$ and $M_1(\pi-\varepsilon)\sim\varepsilon^2/8$.

**Remark (the earlier complex form).** If one instead uses $I_{\rho\sigma}=-\operatorname{Li}_2\big(\frac{t-\rho}{\sigma-\rho}\big)+\big(\log|\rho-\sigma|-i\pi\,[\rho<\sigma]\big)\log(t-\rho)$, with $\operatorname{Li}_2(y)$ meaning $\operatorname{Li}_2(y-i0)$ for $y>1$, then
$$I_{\rho\sigma}=J_{\rho\sigma}+\tfrac12\log^2(\sigma-\rho)-i\pi\log(\sigma-\rho)-\tfrac{\pi^2}{3}\qquad(\rho<\sigma),$$
and $I_{\rho\sigma}=J_{\rho\sigma}$ for $\rho\ge\sigma$. The constants cancel between the endpoints pair by pair, so that form is exactly real and equal to the Theorem. Its literal $\zeta(2)$ coefficient is $-1$ rather than $+1$, because each $\operatorname{Re}\operatorname{Li}_2(y)$ with $y>1$ carries $2\zeta(2)$.

**Remark (the ζ(2) coefficient is not canonical).** Both $\operatorname{Li}_2\frac{s}{1+s}$ and $\operatorname{Li}_2\frac{1}{1+s}$ occur in $U$, and $\operatorname{Li}_2\frac{s}{1+s}+\operatorname{Li}_2\frac1{1+s}=\zeta(2)-\log\frac{s}{1+s}\log\frac{1}{1+s}$. Rewriting with this changes the constant term from $\zeta(2)$ to $2\zeta(2)$. So we only ever speak of the coefficient of $\zeta(2)$ *in a stated representation*.

### Corollary (the coefficient rule at weight 2)

In the Theorem, the coefficient of every $\operatorname{Li}_2$ term, of every product of logarithms of the letters $2,c,s,1\pm c,1+s,1+s\pm c$, and of $\zeta(2)$ lies in $\mathbb Z[\tfrac12]+\mathbb Z[\tfrac12]\,c\subset\mathbb Q+\mathbb Q\cos\frac\alpha2$. The same holds for every representation obtained from it by the reflection, inversion, Landen or duplication formulas. In particular this covers the reduction of all arguments to $(0,\frac12]$.

*Proof.* Each term of $F(t_1)-F(1)$ is a coefficient from $\{0,\pm1,\pm2,\pm c,\pm2c\}$ (a value of $W$ or $\omega$), possibly halved by a $\frac12\log^2$, times one fixed transcendental function: a $\operatorname{Li}_2$, a $\log^2$, or a product of two logs. No coefficient ever multiplies another. The operations that remain are:
- factoring each logarithm into letters, which uses integer exponents;
- merging equal arguments, which adds coefficients;
- applying the standard functional equations, which contribute $\operatorname{Li}_2$, $\zeta(2)$, $\log^2$ and $\log\cdot\log$ terms with factors $\pm1,\pm\frac12,\pm2$.

All of these are $\mathbb Z[\frac12]$-linear maps on the coefficient vector, so they preserve the $\mathbb Z[\frac12]$-module $\mathbb Z[\frac12]+\mathbb Z[\frac12]c$. $\square$

(The sum is direct only when $c\notin\mathbb Q$; for example α = 2π/3 gives $c=\frac12$. As noted above, individual coefficients depend on the representation, so the rule is a statement about these representations, not about invariants. Under $\tau:c\mapsto-c$ with $s$ fixed, $r_1\leftrightarrow r_2$ while $1,-1,c_a,t_1$ are fixed, and the residue vectors are $\tau$-covariant.)

## Not covered here

Part 3 (weight 3, $M_2$) was out of scope for this key and was not audited.
