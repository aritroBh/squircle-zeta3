wd = 'research/08_theory/PRV1_depth2_proofs'
result = r'''# PRV1: Derivative proofs of the depth-2 reductions — regular octagon and decagon

**Status: PROVED (exact symbolic) + numerically verified.** Every identity below is proved
exactly: function identities are proved by formal differentiation (exact rational-function
algebra in sympy) plus an exact constant pin at an endpoint; point identities are proved by
exact algebra over $\mathbb Q(\sqrt2)$, resp. $\mathbb Q(\sqrt5)$. All numerics (mpmath,
120–250 dps) are independent checks, not part of the proofs. One file runs everything:
`verify_full.py` → prints `ALL CHECKS PASSED` (log: `verify_full.out`).

**Objects.** The only depth-2 constants of the regular octagon (OCT6) and decagon (DEC1) are
$$\Lambda=\operatorname{Li}_{2,1}(x_1,\,x_2/x_1),\qquad
  \mathrm{swap}=\operatorname{Li}_{2,1}(x_2,\,x_1/x_2),\qquad
  \operatorname{Li}_{2,1}(x,y)=\sum_{n>m\ge1}\frac{x^n y^m}{n^2 m},$$
with the **Landen pair** $x_2=x_1/(x_1-1)$ (equivalently $x_1x_2=x_1+x_2$):
octagon $x_1=1-\sqrt2$, $x_2=1-\sqrt2/2$; decagon $x_1=-1/\varphi$, $x_2=1/\varphi^2$.

**Results.**
* **Theorem O.** Both octagon values reduce to depth 1: $\Lambda_{\rm oct}$ and
  $\mathrm{swap}_{\rm oct}$ equal the explicit $\{\zeta(3),\operatorname{Li}_3,\operatorname{Li}_2,\log\}$
  combinations of §1.5, which are **exactly** OCT6's stored reductions (`reductions.pkl`;
  symbolic difference $0$).
* **Theorem D.** Both decagon values reduce to depth 1: $\Lambda_{\rm dec}$ and
  $\mathrm{swap}_{\rm dec}$ equal DEC1's LLL forms
  $$\Lambda_{\rm dec}=\tfrac25\zeta(3)+\tfrac45\,l_\varphi z_2-l_\varphi^3
    +3\operatorname{Li}_3(-1/\varphi)+\operatorname{Li}_3(1/\varphi),$$
  $$\mathrm{swap}_{\rm dec}=\tfrac35\zeta(3)-\tfrac75\,l_\varphi z_2+2l_\varphi^3
    +\operatorname{Li}_3(-1/\varphi)+\operatorname{Li}_3(1/\varphi),$$
  $l_\varphi=\log\varphi$, $z_2=\zeta(2)=\pi^2/6$ (exact symbolic difference $0$ modulo the
  golden $\operatorname{Li}_3$ sum relation, which is itself proved). Previously these two were
  numerical only (LLL, 830 digits, residuals $3.6{\cdot}10^{-832}$/$1.8{\cdot}10^{-830}$).
* **Lemma L (§3).** Both are the members $c=\tfrac12$, $c=1$ of one line-curve family
  through any Landen pair; the derivative always reduces to depth 1 along the curve, and the
  primitive closes in single polylogarithms for $c\in\{\tfrac12,1\}$, covering the regular
  octagon ($m=4$) and decagon ($m=5$). The $c=\tfrac23$ member (hexagon-type Landen pair
  $x_1=-\tfrac12$, $x_2=\tfrac13$) shows **no** classical reduction (PSLQ, §3.3) — and the
  hexagon's depth-2 constant is in fact of Nielsen $S_{1,2}$ type, already reduced by K1e/K1f.

---

## 1. The octagon

### 1.1 (Task a) Partial derivatives of $\operatorname{Li}_{2,1}$ and the $\operatorname{Li}_{1,1}$ closed form

**Lemma 1 (partial derivatives).** On the polydisc $|x|<1$, $|y|<1$ the defining series of
$\operatorname{Li}_{2,1}$ converges absolutely and uniformly on compact subsets, may be
differentiated term by term, and
$$\partial_x\operatorname{Li}_{2,1}(x,y)=\frac{\operatorname{Li}_{1,1}(x,y)}{x},
\qquad
\partial_y\operatorname{Li}_{2,1}(x,y)=\frac{\operatorname{Li}_2(x)-\operatorname{Li}_2(xy)/y}{1-y},$$
where $\operatorname{Li}_{1,1}(x,y)=\sum_{n>m\ge1}x^ny^m/(nm)$.

*Proof.* The first identity is $\sum_{n>m}x^{n-1}y^m/m$ termwise. For the second,
$\partial_y\sum_{n>m}\frac{x^ny^m}{n^2m}=\sum_{n>m}\frac{x^ny^{m-1}}{n^2}
=\sum_n\frac{x^n}{n^2}\frac{1-y^{n-1}}{1-y}
=\frac{1}{1-y}\big(\operatorname{Li}_2(x)-\operatorname{Li}_2(xy)/y\big)$
(geometric sum over $1\le m\le n-1$). $\square$
[`verify_full.py` Part 0.1: both identities checked as **exact finite-sum identities**
(truncation $N=25$; the $y$-identity in the multiplied-out polynomial form
$(1-y)\,\partial_y S_N=\sum_{n<N}\frac{x^n}{n^2}(1-y^{n-1})$, residual $0$); the division by
$1-y$ is then formal. Domain of validity: $|x|<1$, $|y|<1$, $y\ne1$; the formulas extend
continuously to $y\to0$.]

**Lemma 2 ($\operatorname{Li}_{1,1}$ closed form).** For real $-1<x<1$, $-1<y<1$,
$$\operatorname{Li}_{1,1}(x,y)
=\log(1-x)\log(1-y)+\operatorname{Li}_2\!\Big(\frac{-y}{1-y}\Big)
-\operatorname{Li}_2\!\Big(\frac{-y(1-x)}{1-y}\Big).$$

*Proof.* From the series, $\partial_x\operatorname{Li}_{1,1}(x,y)
=\sum_{m\ge1}\frac{y^m}{m}\sum_{n>m}x^{n-1}
=\frac{1}{1-x}\sum_{m\ge1}\frac{(xy)^m}{m}=-\frac{\log(1-xy)}{1-x}$.
Differentiating the right side (with $z=-y(1-x)/(1-y)$, $z'/z=-1/(1-x)$,
$1-z=(1-xy)/(1-y)$):
$$\partial_x{\rm RHS}=-\frac{\log(1-y)}{1-x}
+\log\!\frac{1-xy}{1-y}\cdot\frac{-1}{1-x}=-\frac{\log(1-xy)}{1-x},$$
using the real-log split $\log\frac{1-xy}{1-y}=\log(1-xy)-\log(1-y)$ — valid since all
arguments are positive for $-1<x,y<1$. Both sides vanish at $x=0$ (RHS:
$\log1\cdot\log(1-y)+\operatorname{Li}_2-\operatorname{Li}_2=0$), so they agree.
$\square$
[Part 0.2: derivative difference is **exactly** $0$ in sympy after the branch-correct log
split; the vanishing at $x=0$ is exact.]

### 1.2 (Task b) The one-parameter families and the chain rule

Define on $0<t<1$
$$G(t)=\operatorname{Li}_{2,1}\!\Big(-t,\,-\frac{1-t}{2t}\Big),\qquad
H(t)=\operatorname{Li}_{2,1}\!\Big(\frac{1-t}{2},\,-\frac{2t}{1-t}\Big).$$
In $(x,p)$-coordinates with $p=xy$ these are the straight lines $(x,p)=(-t,\tfrac{1-t}{2})$
and its swap $(x,p)=(\tfrac{1-t}{2},-t)$. On $0<t<1$: for $G$, $|x|=t<1$ and
$|xy|=(1-t)/2<1$; for $H$, $|x|=(1-t)/2<1$ and $|xy|=t<1$; hence both series converge and
$G,H\in C^1(0,1)$. The chain rule with Lemma 1 gives
$F'(t)=\partial_x\operatorname{Li}_{2,1}\,x'(t)+\partial_y\operatorname{Li}_{2,1}\,y'(t)$.

On the curves the $\operatorname{Li}_{1,1}$ arguments simplify rationally
(**all checked exactly**, `verify_full.py` CHAIN-RULE ASSEMBLY):
$$\begin{aligned}
G:&\quad 1-x=1+t,\quad 1-y=\tfrac{1+t}{2t},\quad \tfrac{-y}{1-y}=\tfrac{1-t}{1+t}=:c,\quad
\tfrac{-y(1-x)}{1-y}=1-t,\quad xy=\tfrac{1-t}{2}=:b,\\
H:&\quad 1-x=\tfrac{1+t}{2}=:h,\quad 1-y=\tfrac{1+t}{1-t},\quad \tfrac{-y}{1-y}=\tfrac{2t}{1+t}=:u,\quad
\tfrac{-y(1-x)}{1-y}=t,\quad xy=-t .
\end{aligned}$$
With $R=\log t$, $B=\log(1+t)$, $q=\log\frac{1+t}{2}$, $l_c=\log c$, $U=\log u$ (all real
logs of positive arguments on $0<t<1$):
$$G'(t)=\frac{J_G}{t}+\frac{\operatorname{Li}_2(-t)}{t(1+t)}+\frac{2\operatorname{Li}_2(b)}{1-t^2},
\qquad J_G=(q-R)B+\operatorname{Li}_2(c)-\operatorname{Li}_2(1-t),$$
$$H'(t)=-\frac{J_H}{1-t}-\frac{\operatorname{Li}_2(-t)}{t(1+t)}-\frac{2\operatorname{Li}_2(b)}{1-t^2},
\qquad J_H=-l_cq+\operatorname{Li}_2(u)-\operatorname{Li}_2(t).$$
[These are *assembled* from Lemmas 1–2 inside the verifier and match the stated formulas with
residual $0$, exactly.]

### 1.3 (Task c) The closed forms and the exact derivative identities

$$G_{\rm cl}(t)=\tfrac34\zeta(3)+\operatorname{Li}_3(t)+\operatorname{Li}_3(-t)
+\operatorname{Li}_3(-c)-\operatorname{Li}_3(u)
-q\Big(\operatorname{Li}_2(t)+\tfrac12\operatorname{Li}_2(b)-\tfrac12\operatorname{Li}_2(h)\Big),$$
$$H_{\rm cl}(t)=\operatorname{Li}_3(1-t)+\operatorname{Li}_3(b)+\operatorname{Li}_3(w)
-\operatorname{Li}_3(c)-\operatorname{Li}_3(\tfrac12)
-B\big(\operatorname{Li}_2(c)+\operatorname{Li}_2(b)\big)+\tfrac12B^2U,$$
with $w=t/(1+t)$.

**Classical function identities used** (each proved in Part 1 by differentiation — exact
residual $0$ — plus an exact endpoint pin):
$$\begin{aligned}
{\rm E}:&\ \operatorname{Li}_2(t)+\operatorname{Li}_2(1-t)=\zeta(2)-\log t\,\log(1-t)
   &&\text{pin }t\to0^+:\ \operatorname{Li}_2(1)=\zeta(2),\\
{\rm L1}:&\ \operatorname{Li}_2(-t)+\operatorname{Li}_2(w)=-\tfrac12B^2
   &&\text{pin }t=0,\\
{\rm L2}:&\ \operatorname{Li}_2(-c)+\operatorname{Li}_2(b)=-\tfrac12q^2
   &&\text{pin }t=0:\ -\tfrac{\zeta(2)}2+\big(\tfrac{\zeta(2)}2-\tfrac{\log^22}2\big)+\tfrac{\log^22}2=0,\\
{\rm DUP2}:&\ \operatorname{Li}_2(t^2)=2\operatorname{Li}_2(t)+2\operatorname{Li}_2(-t)
   &&\text{pin }t=0,\\
{\rm OC}:&\ \operatorname{Li}_2(t)-\operatorname{Li}_2(-t)+\operatorname{Li}_2(c)-\operatorname{Li}_2(-c)
   =\tfrac32\zeta(2)-\log t\,\log c
   &&\text{pin }t\to1^-:\ \zeta(2)+\tfrac{\zeta(2)}2=\tfrac32\zeta(2),\\
{\rm DUP3}:&\ \operatorname{Li}_3(t^2)=4\operatorname{Li}_3(t)+4\operatorname{Li}_3(-t)
   &&\text{pin }t=0.
\end{aligned}$$
($\operatorname{Li}_2(-1)=-\zeta(2)/2$ from DUP2 at $t=1$;
$\operatorname{Li}_2(1/2)=\zeta(2)/2-\log^22/2$ from E at $t=1/2$;
$\operatorname{Li}_3(-1)=-3\zeta(3)/4$ from DUP3 at $t=1$. The pin for L2 uses the first two;
$\log t\log c\to0$ as $t\to1^-$.)

**Derivative identities (exact).**
$$\frac{d}{dt}H_{\rm cl}(t)-H'(t)\equiv0\pmod{{\rm E,\,L1,\,L2}},\qquad
\frac{d}{dt}G_{\rm cl}(t)-G'(t)\equiv-\frac{{\rm OC}(t)}{1+t}\pmod{{\rm E,\,L1,\,L2}};$$
since ${\rm OC}\equiv0$, both closed forms have the same derivatives as $G,H$ on $(0,1)$.
[Part 2: formal differentiation of the closed forms in the letter algebra
$\{R,A,B,C,\zeta_2,\zeta_3,\operatorname{Li}_2(\cdot),\operatorname{Li}_3(\cdot)\}$ with the
derivative table from $\frac d{dz}\operatorname{Li}_n=\operatorname{Li}_{n-1}/z$; the
differences simplify **exactly** to $0$, resp. to $-{\rm OC}/(1+t)$, after substitution of
the proved identities E, L1, L2.]

### 1.4 (Task d) The constants

$|G(t)|\le\operatorname{Li}_2(|x|)\,|y|/(1-|y|)\le\zeta(2)\,|y|/(1-|y|)\to0$ as $t\to1^-$
(since $y=-(1-t)/(2t)\to0$); and
$$G_{\rm cl}(1)=\tfrac34\zeta(3)+\zeta(3)+\operatorname{Li}_3(-1)+0-\zeta(3)-0
=\tfrac34\zeta(3)+\zeta(3)-\tfrac34\zeta(3)-\zeta(3)=0.$$
$|H(t)|\to0$ as $t\to0^+$ (same estimate, $y=-2t/(1-t)\to0$); and
$$H_{\rm cl}(0)=\zeta(3)+\operatorname{Li}_3(\tfrac12)+0-\zeta(3)
-\operatorname{Li}_3(\tfrac12)-0+0=0.$$
Two $C^1$ functions on $(0,1)$ with equal derivatives and equal endpoint limits coincide:
$$G\equiv G_{\rm cl},\qquad H\equiv H_{\rm cl}\quad\text{on }(0,1).\qquad\square$$

### 1.5 (Task e) Specialisation at $t=k=\sqrt2-1$

At $t=k$ (all argument identities **exact over $\mathbb Q(\sqrt2)$**, Part 2b):
$b(k)=x_2$, $h(k)=\sqrt2/2$, $c(k)=k$, $u(k)=2-\sqrt2=:d$, $w(k)=x_2$,
$k(1+\sqrt2)=1$, $1-k=\sqrt2\,k$, $x_2=x_1/(x_1-1)$, and the second arguments of $G,H$ are
$-\sqrt2/2=x_2/x_1$, $-\sqrt2=x_1/x_2$ respectively: so $G(k)=\Lambda_{\rm oct}$,
$H(k)=\mathrm{swap}_{\rm oct}$. The silver dilog rules are exact specialisations of E and
Landen at points (log identities verified as exact *multiplicative* identities among the
positive arguments, e.g. $d(1+\sqrt2)=\sqrt2\iff\log d=\tfrac12\log2-\log(1+\sqrt2)$).
With $l_2=\log2$, $l_u=\log(1+\sqrt2)$:
$$\boxed{\;\Lambda_{\rm oct}=\tfrac34\zeta(3)+\operatorname{Li}_3(k)+2\operatorname{Li}_3(x_1)
-\operatorname{Li}_3(d)+\tfrac{l_2}{2}\operatorname{Li}_2(k)
+\tfrac{l_2}{4}\big(\operatorname{Li}_2(x_2)-\operatorname{Li}_2(h)\big)\;}$$
$$\boxed{\;\mathrm{swap}_{\rm oct}=\operatorname{Li}_3(d)+2\operatorname{Li}_3(x_2)
-\operatorname{Li}_3(k)-4\operatorname{Li}_3(h)-4\operatorname{Li}_3(-h)
-\tfrac{l_2}{2}\big(\operatorname{Li}_2(k)+\operatorname{Li}_2(x_2)\big)
+\tfrac{l_2^2}{8}\big(\tfrac{l_2}{2}-l_u\big)\;}$$
($x_1=1-\sqrt2=-k$, $x_2=1-\sqrt2/2$, $d=2-\sqrt2$, $h=\sqrt2/2$; the $\operatorname{Li}_3(1/2)$
in $H_{\rm cl}$ is rewritten by DUP3 at $h$, since $h^2=1/2$.)
**These equal OCT6's stored reductions exactly** (Part 2b: symbolic difference with
`OCT6_octagon_closed/reductions.pkl` is $0$ for both, after the silver dilog rules — an exact
rigid comparison, not a fit). Numeric cross-check (120 dps): series
$\Lambda_{\rm oct}=-0.0275883613456533384687572087639\ldots$ vs $G_{\rm cl}(k)$ residual
$1.7{\cdot}10^{-121}$; $\mathrm{swap}_{\rm oct}=-0.0321696900151107721899470611759\ldots$ vs
$H_{\rm cl}(k)$ residual $1.6{\cdot}10^{-121}$. ∎

---

## 2. The decagon

### 2.1 The families, chain rule, closed forms

Same shape with $c=1$ (§3): on $0<t<1$,
$$G_{\rm dec}(t)=\operatorname{Li}_{2,1}\!\Big(-t,\,-\frac{1-t}{t}\Big),\qquad
H_{\rm dec}(t)=\operatorname{Li}_{2,1}\!\Big(1-t,\,-\frac{t}{1-t}\Big),$$
curves $(x,p)=(-t,\,1-t)$ and its swap $(1-t,\,-t)$; convergence and $C^1$ as before
($|x|,|xy|\in\{t,1-t\}<1$). Curve simplifications (checked exactly):
for $G_{\rm dec}$: $1-y=1/t$, $-y/(1-y)=1-t$, $-y(1-x)/(1-y)=1-t^2$, $xy=1-t$;
for $H_{\rm dec}$: $1-x=t$, $1-y=1/(1-t)$, $-y/(1-y)=t$, $-y(1-x)/(1-y)=t^2$, $xy=-t$.
$$G_{\rm dec}'(t)=\frac{-RB+\operatorname{Li}_2(1-t)-\operatorname{Li}_2(1-t^2)}{t}
+\frac{\operatorname{Li}_2(-t)}{t}+\frac{\operatorname{Li}_2(1-t)}{1-t},$$
$$H_{\rm dec}'(t)=-\frac{-RA+\operatorname{Li}_2(t)-\operatorname{Li}_2(t^2)}{1-t}
-\frac{\operatorname{Li}_2(1-t)}{1-t}-\frac{\operatorname{Li}_2(-t)}{t},\qquad A=\log(1-t).$$

**Closed forms.**
$$G_{{\rm dec},{\rm cl}}(t)=\zeta(3)-\tfrac12(\log t)\operatorname{Li}_2(t^2)
+\operatorname{Li}_3(t^2)-2\operatorname{Li}_3(t)-\operatorname{Li}_3(1-t),$$
$$H_{{\rm dec},{\rm cl}}(t)=\operatorname{Li}_3(1-t^2)-\operatorname{Li}_3\!\tfrac1{1+t}
+\operatorname{Li}_3\!\tfrac{t}{1+t}
-B\Big(\tfrac12\operatorname{Li}_2(-t)+\tfrac32\operatorname{Li}_2(1-t)
+\tfrac14\operatorname{Li}_2(1-t^2)+\tfrac14\zeta(2)\Big).$$

**Derivative identities (exact):**
$\frac{d}{dt}G_{{\rm dec},{\rm cl}}-G_{\rm dec}'\equiv0$ and
$\frac{d}{dt}H_{{\rm dec},{\rm cl}}-H_{\rm dec}'\equiv0$
modulo $\{$DUP2, E, E at $t^2$, L1, E at $v=t/(1+t)\}$ (Part 3: residuals exactly $0$).

### 2.2 Constants

$G_{\rm dec}(t)\to0$ as $t\to1^-$ (same bound, $y=-(1-t)/t\to0$);
$G_{{\rm dec},{\rm cl}}(1)=\zeta(3)-0+\zeta(3)-2\zeta(3)-0=0$.
$H_{\rm dec}(t)\to0$ as $t\to0^+$ ($y=-t/(1-t)\to0$);
$H_{{\rm dec},{\rm cl}}(0)=\zeta(3)-\zeta(3)+0-0=0$.
Hence $G_{\rm dec}\equiv G_{{\rm dec},{\rm cl}}$, $H_{\rm dec}\equiv H_{{\rm dec},{\rm cl}}$
on $(0,1)$. $\square$

### 2.3 Golden lattice and specialisation at $t=1/\varphi$

**Golden dilogs** (E at $1/\varphi^2$, Landen at $1/\varphi^2$, DUP2 at $1/\varphi$ —
three exact function identities at three points; the $3\times3$ linear system is solved
exactly):
$$\operatorname{Li}_2(1/\varphi^2)=\tfrac25z_2-l_\varphi^2,\qquad
\operatorname{Li}_2(1/\varphi)=\tfrac35z_2-l_\varphi^2,\qquad
\operatorname{Li}_2(-1/\varphi)=-\tfrac25z_2+\tfrac12l_\varphi^2.$$
**Golden trilogs:** the 3-term Landen identity
$$\operatorname{Li}_3(t)+\operatorname{Li}_3(1-t)+\operatorname{Li}_3\!\Big(\frac{-t}{1-t}\Big)
=\zeta(3)+\zeta(2)\log(1-t)-\tfrac12\log t\,\log^2(1-t)+\tfrac16\log^3(1-t)$$
(proved in Part 3b by differentiation + pin at $t\to0^+$, exactly, using E and Landen) at
$t=1/\varphi^2$ gives
$\operatorname{Li}_3(1/\varphi^2)+\operatorname{Li}_3(1/\varphi)+\operatorname{Li}_3(-1/\varphi)
=\zeta(3)-l_\varphi z_2+\tfrac56l_\varphi^3$;
with DUP3 at $1/\varphi$ ($\operatorname{Li}_3(1/\varphi^2)=4\operatorname{Li}_3(1/\varphi)
+4\operatorname{Li}_3(-1/\varphi)$):
$$\operatorname{Li}_3(1/\varphi^2)=\tfrac45\big(\zeta(3)-l_\varphi z_2\big)+\tfrac23l_\varphi^3
\ \text{(Lewin)},\qquad
\operatorname{Li}_3(1/\varphi)+\operatorname{Li}_3(-1/\varphi)
=\tfrac15\big(\zeta(3)-l_\varphi z_2\big)+\tfrac16l_\varphi^3.$$
(The difference $\operatorname{Li}_3(1/\varphi)-\operatorname{Li}_3(-1/\varphi)$ is genuinely
undetermined — classically these two are independent constants; DEC1 confirmed no relation at
830 digits.)

**Specialisation.** At $t=1/\varphi$: $G_{\rm dec}(1/\varphi)
=\operatorname{Li}_{2,1}(-1/\varphi,-1/\varphi)=\Lambda_{\rm dec}$ and
$H_{\rm dec}(1/\varphi)=\operatorname{Li}_{2,1}(1/\varphi^2,-\varphi)=\mathrm{swap}_{\rm dec}$
(second arguments: $-(1-t)/t=-1/\varphi$, $-t/(1-t)=-\varphi$, exact over
$\mathbb Q(\sqrt5)$). Therefore
$$\Lambda_{\rm dec}=\zeta(3)+\tfrac12l_\varphi\operatorname{Li}_2(1/\varphi^2)
-2\operatorname{Li}_3(1/\varphi)
=\zeta(3)+\tfrac15l_\varphi z_2-\tfrac12l_\varphi^3-2\operatorname{Li}_3(1/\varphi),$$
$$\mathrm{swap}_{\rm dec}=\operatorname{Li}_3(1/\varphi^2)+\tfrac32l_\varphi^3
-\tfrac45l_\varphi z_2
=\tfrac45\zeta(3)-\tfrac85l_\varphi z_2+\tfrac{13}6l_\varphi^3\quad(\text{via Lewin}).$$
Reducing modulo the proved sum relation gives **exactly DEC1's forms** (Part 3c: symbolic
differences $0$):
$$\Lambda_{\rm dec}=\tfrac25\zeta(3)+\tfrac45l_\varphi z_2-l_\varphi^3
+3\operatorname{Li}_3(-1/\varphi)+\operatorname{Li}_3(1/\varphi),$$
$$\mathrm{swap}_{\rm dec}=\tfrac35\zeta(3)-\tfrac75l_\varphi z_2+2l_\varphi^3
+\operatorname{Li}_3(-1/\varphi)+\operatorname{Li}_3(1/\varphi).\qquad\square$$
Numeric cross-check (170 dps): series $\Lambda_{\rm dec}
=-0.05126163058647632738126326036571123343019\ldots$ vs DEC1 form residual
$4.5{\cdot}10^{-172}$; $\mathrm{swap}_{\rm dec}
=-0.06341782168216461215085970440865279340954\ldots$ vs DEC1 form residual
$1.2{\cdot}10^{-171}$. ∎

---

## 3. The general lemma

**Lemma L (Landen-pair line curve).** Let $x_1\in(-1,0)$, $x_2=x_1/(x_1-1)\in(0,1)$
(equivalently $x_1x_2=x_1+x_2$), and put
$$c=\frac{x_2}{1+x_1}>0.$$
(i) The line $p=c(1+x)$ in $(x,p=xy)$-space passes through $(x_1,x_2)$ and $(-1,0)$; the
parametrisation $x=-t$ gives the $C^1$ family on $0<t<1$
$$G_c(t)=\operatorname{Li}_{2,1}\!\Big(-t,\,-c\,\frac{1-t}{t}\Big),\qquad
H_c(t)=\operatorname{Li}_{2,1}\!\Big(c(1-t),\,-\frac{t}{c(1-t)}\Big)\ \ (x\leftrightarrow p\
\text{swap}),$$
with $G_c(t)\to0$ as $t\to1^-$ and $H_c(t)\to0$ as $t\to0^+$.
(ii) For every $c>0$ the derivatives $G_c',H_c'$ reduce **exactly to depth 1** (logs and
dilogs at rational arguments): with $T_c(t)=\dfrac{c(1-t)}{c+(1-c)t}$ — an involution,
$T_c\circ T_c=\mathrm{id}$ (checked exactly) — the $\operatorname{Li}_{1,1}$ arguments on the
curve are $T_c(t)$ and $(1+t)T_c(t)$, so Lemma 2 expresses $G_c',H_c'$ through
$\operatorname{Li}_2$ at $t$, $T_c(t)$, $(1+t)T_c(t)$ and logs.
(iii) The primitive **closes in single polylogarithms** (one-dimensional
$\operatorname{Li}_3/\operatorname{Li}_2/\log$ combinations) for $c=\tfrac12$ and $c=1$:
* $c=\tfrac12$: $T_{1/2}(t)=\frac{1-t}{1+t}$ is the Landen involution, $(1+t)T_{1/2}=1-t$;
  $G_{1/2},H_{1/2}=G,H$ of §1, closed forms $G_{\rm cl},H_{\rm cl}$, closure modulo
  $\{$E, L1, L2$\}$ + OC. At $t=\sqrt2-1$ the curve hits the octagon's Landen pair
  $(x_1,x_2)=(1-\sqrt2,\,1-\sqrt2/2)$. **(Theorem O.)**
* $c=1$: $T_1(t)=1-t$, $(1+t)T_1=1-t^2$; $G_1,H_1=G_{\rm dec},H_{\rm dec}$ of §2, closed
  forms $G_{{\rm dec},{\rm cl}},H_{{\rm dec},{\rm cl}}$, closure modulo
  $\{$E, L1, DUP2$\}$. At $t=1/\varphi$ the curve hits the decagon's Landen pair
  $(-1/\varphi,\,1/\varphi^2)$. **(Theorem D.)**

**Which regular $2m$-gons this covers.** The two-ratio edges of a regular $2m$-gon produce
exactly the pair $\{\Lambda,\mathrm{swap}\}$ with $x_2=x_1/(x_1-1)$ (DEC1 §2, the same
structure as the octagon). Lemma L applies with $c=x_2/(1+x_1)$:
* $m=4$ (octagon): $c=\tfrac12$ — **proved** (Theorem O);
* $m=5$ (decagon): $c=1$ — **proved** (Theorem D);
* $m=2$ (square): no depth-2 constants at all (elementary, E1);
* $m=3$ (hexagon): the Landen pair would be $(x_1,x_2)=(-\tfrac12,\tfrac13)$, $c=\tfrac23$;
  **no classical closure** — a PSLQ search (150 dps, 19-element weight-3 basis
  $\{\zeta(3),\operatorname{Li}_3(\pm\tfrac12,\pm\tfrac13,\pm\tfrac23,\pm\tfrac25,\pm\tfrac35,
  \tfrac14,\tfrac34),\zeta(2)\log2,\zeta(2)\log3,\log^3,\log^2\log\}$) found no relation for
  $\operatorname{Li}_{2,1}(-\tfrac12,-\tfrac23)$: the only candidate had coefficients
  $\sim5{\cdot}10^7$ and residual $2.2{\cdot}10^{-110}$ **unchanged at 250 dps**, i.e. a
  lattice artifact. The regular hexagon's actual depth-2 content is of Nielsen $S_{1,2}$ type
  and was reduced by the Nielsen–Kölbig identity (K1e/K1f), a different mechanism.
* $m\ge6$: no closing cases of the present family are known; the closure condition is the
  rational-integrability of $\operatorname{Li}_2(T_c(t))$-type terms, special to
  $c\in\{\tfrac12,1\}$ among the cases tested.

So the two proved cases are exactly the two regular polygons whose depth-2 constants were
outstanding: the depth-2 part of the regular-octagon and regular-decagon $\zeta(3)$ problems
is now proved, by one and the same method.

---

## 4. Computer verification

`verify_full.py` (extends CX1's `verify.py`; run: `PY verify_full.py` → `ALL CHECKS PASSED`,
log `verify_full.out`):

| Part | Content | Nature |
|---|---|---|
| 0 | Lemma 1 (truncated-series identities, residual 0); Lemma 2 ($d/dx$ residual 0, pin at $x=0$) | exact sympy |
| 1 | E, L1, L2, DUP2, OC derivatives exactly 0; DUP3 after DUP2; all endpoint pins exact | exact sympy |
| CR | 20 curve argument identities + assembly of $G',H',G_{\rm dec}',H_{\rm dec}'$ from Lemmas 1–2 (residuals 0) | exact sympy |
| 2 | $d/dt\,G_{\rm cl}-G'=-{\rm OC}/(1+t)$, $d/dt\,H_{\rm cl}-H'=0$ (mod E,L1,L2); constants | exact sympy |
| 2b | specialisation at $\sqrt2-1$ over $\mathbb Q(\sqrt2)$; **exact match with OCT6 `reductions.pkl`** | exact sympy |
| 3 | decagon derivative identities (mod DUP2,E,L1); constants | exact sympy |
| 3b | 3-term Landen trilog identity (derivative 0, pin) | exact sympy |
| 3c | golden dilog/trilog lattice; specialisation at $1/\varphi$; **exact match with DEC1 forms** | exact sympy |
| 4 | independent numeric evidence (below) | mpmath, checks only |

**Numeric evidence table** (residual = |series value − closed form|; all at stated precision):

| Identity | dps | worst/actual residual |
|---|---|---|
| E, L1, L2, DUP2, OC, DUP3, Landen-3 at $t=0.3,0.55,0.8$ | 120 | $3.9{\cdot}10^{-121}$ |
| Lemma 2 vs series at $(0.3,0.4),(-0.4,0.5),(0.1,-0.6)$ | 120 | $5.9{\cdot}10^{-122}$ |
| $G=G_{\rm cl}$, $H=H_{\rm cl}$ at $t=0.17\ldots0.82$ (5 pts) | 120 | $1.9{\cdot}10^{-121}$ |
| $\Lambda_{\rm oct}=G_{\rm cl}(k)$, $\mathrm{swap}_{\rm oct}=H_{\rm cl}(k)$ | 120 | $1.7{\cdot}10^{-121}$ |
| $G_{\rm dec}=G_{{\rm dec},{\rm cl}}$, $H_{\rm dec}=H_{{\rm dec},{\rm cl}}$ (5 pts) | 120 | $1.9{\cdot}10^{-121}$ |
| $\Lambda_{\rm dec}=$ DEC1 form | 170 | $4.5{\cdot}10^{-172}$ |
| $\mathrm{swap}_{\rm dec}=$ DEC1 form | 170 | $1.2{\cdot}10^{-171}$ |
| golden lattice (Lewin, sum relation, 3 dilogs) | 170 | $5.2{\cdot}10^{-172}$ |
| $c=\tfrac23$ member: no relation (PSLQ artifact check) | 150→250 | residual frozen at $2.2{\cdot}10^{-110}$ ⇒ spurious |

Series engine: $\operatorname{Li}_{2,1}(x,y)=\sum_{n>m}x^ny^m/(n^2m)$ summed by the
recurrence $h_n=xh_{n-1}+x(xy)^n/n$, truncation $N$ with
$|\text{tail}|\le\rho^{N+1}/(1-\rho)$, $\rho=\max(|x|,|xy|)$ (same engine as CX1's verifier).

**Relation to prior work.** CX1 (`CX1_octagon_depth2_proof/verify.py`) proved the octagon
derivative identities + pins and matched OCT6 numerically; the present work (i) makes every
step exact (Part 0–3c are symbolic; the only numerical content is the independent Part 4
evidence), (ii) re-derives the chain-rule formulas from first principles inside the verifier
(CR), (iii) proves the decagon reductions (previously 830-digit LLL conjectures), and
(iv) isolates the general mechanism (Lemma L) and its scope.

**Files.** `verify_full.py` (+ `verify_full.out`), `result.md`, `progress.log`,
`step_01.py`–`step_NN.py` (derivation trail: exploration of the families, closed-form
guessing by symbolic integration, golden-lattice algebra, PSLQ scope tests).
'''
with open(f'{wd}/result.md','w') as f:
    f.write(result)
print("result.md written:", len(result), "chars")
