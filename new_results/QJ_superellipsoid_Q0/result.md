# QJ: the closed form of $J$ and the integrals $P_1,P_{2a},P_{2c}$ (vertex constant $Q_0$ of thm:S)

**Verdict: PROVED BY HAND + EXACT CAS CHECKS.** Both boxed steps of `paper/ieee/proofs/sec_g3_S.tex` now have
written proofs: (1) the closed form of $J$ (referee item G2), and (2) the three evaluations of
Lemma `g3:lem:P` (referee item E5). Every substitution and every antiderivative is checked exactly
by sympy (derivative minus integrand $=0$, or the transformed integrand equals the claimed one).
Every branch is stated. The dilogarithm values are reduced by **15 explicit instances** of reflection,
inversion, Landen, duplication and Abel's five-term relation (with their validity domains checked
exactly), plus $\Li_2(\tfrac12)$, $\Li_2(-1)$ and Lemma `g3:lem:A` of the supplement. The 54-point
linear-algebra engine of the notes is no longer needed, and neither is $\Li_2(\tfrac13)$.
Assembly: $6(P_1-P_{2a}-P_{2c})$ minus the printed $Q_0$ is **exactly 0** (sympy).

Numerical guards (none of them is used in the proof):

| check | precision | residual |
|---|---|---|
| $J$ closed form vs quadrature, at $q=1/\sqrt2,0.75,0.85,0.95,0.999$ | 80 dps | $\le2.4\cdot10^{-80}$ |
| each change of variables (P2c, P2a: 5 forms of $W$; P1: $v$-form) vs 100-digit reference | 90 dps | $\le6\cdot10^{-92}$ |
| $P_1,P_{2a},P_{2c}$ (exact, unreduced) vs 100-digit quadrature | 90 dps | $\le5\cdot10^{-91}$ |
| each complex boundary value (8 instances) | 90 dps | $\le2.5\cdot10^{-91}$ |
| each of the 15 real identity instances | 110 dps | $\le6.7\cdot10^{-111}$ |
| $P_1,P_{2a},P_{2c}$ (reduced closed forms) vs quadrature | 110 dps | $\le5\cdot10^{-101}$ |
| $Q_0$ closed form vs $-\pi^{-1/2}\int_0^\infty h_0^3$ | 75 digits | $3.1\cdot10^{-80}$ |
| $Q_0$ closed form vs the simplex integral of Lemma `g3:lem:Q0simplex` (2-D quadrature) | 65 digits | $4.3\cdot10^{-70}$ |
| $Q_0$ closed form vs $6\int_0^1J(q(s))\dd s$ (quadrature of the $J$ formula) | 100 digits | $2.2\cdot10^{-101}$ |

$Q_0=-2.832935683494847760462057251515769215937683323767941228331713531285211573488576468799720436503511083\ldots$

## What is still only numerical

Nothing in the chain of proof. What a reader must accept rather than recheck:

1. The classical functional equations, with the domains stated below (Lewin, *Polylogarithms and
   associated functions*, ch. 1; Zagier, *The dilogarithm function*, §1). For complex arguments we use
   only: reflection and inversion of $\Li_2$, Landen, duplication, and two consequences for
   $\chi_2$ (inversion, Landen), whose domains we derive by analytic continuation.
2. The value $\operatorname{Re}\Li_2(e^{i\theta})=\frac{\pi^2}6-\frac{\theta(2\pi-\theta)}4$, $0\le\theta\le2\pi$
   (the Fourier series of the Bernoulli polynomial $B_2$).
3. Lemma `g3:lem:A` of the supplement (already proved there by hand).
4. That sympy's `simplify` returning $0$ is correct. Every such identity is rational-function algebra in
   $\QQ(\sqrt2)(v)$, or in one variable with one square root of a positive quantity, so it can also be
   checked by hand. Where logarithms are combined, we work with the logarithms as independent
   symbols and use only $\log(xy)=\log x+\log y$ for $x,y>0$.

One check is numerical only and is not needed: in `04_P1.py`, Step 1, the equality of the expanded
integrand with $\frac{v^2+1}{2v^2}f(q(v))$ is checked at three points. It follows exactly from the two
factorisations $1+q=\frac{(v+a)(v+b)}{2\sqrt2v}$, $2q=\frac{v^2+1}{\sqrt2v}$ (checked exactly) and
$\log$ of a product of positive numbers.

## Files

| script | what it proves / checks | log |
|---|---|---|
| `common.py` | exact arithmetic in $\QQ(\sqrt2)$, exact logs of $S$-units ($\log x=\frac i2\lambda+j\log3+kL$), numerics | |
| `00_reference.py` | 100-digit quadratures of $P_1,P_{2a},P_{2c}$ (independent of everything else) | `00_reference.log`, `reference_100.json` |
| `01_J.py` | closed form of $J$ (Part A) | `01_J.log` |
| `02_P2c.py` | $P_{2c}$ (Part B), exact closed form | `02_P2c.log`, `P2c_exact.json` |
| `03_P2a.py` | $P_{2a}$ (Part C) up to real dilogarithms | `03_P2a.log`, `P2a_unreduced.json` |
| `04_P1.py` | $P_1$ (Part D) up to real dilogarithms, with a term-by-term table | `04_P1.log`, `P1_unreduced.json` |
| `05_reduce.py` | the 15-instance reduction (Part E), exact closed forms of $P_1,P_{2a},P_{2c}$, exact assembly of $Q_0$ (Part F) | `05_reduce.log`, `05_reduce_out.json` |
| `09_Q0_direct.py` | $Q_0$ from $-\pi^{-1/2}\int k_0^3$ (arg `75 k`) and from the 2-D simplex integral (arg `65 s`) | `09_Q0_k3.log`, `09_Q0_simplex.log`, `Q0_k3.json`, `Q0_simplex_65.json` |
| `run_all.sh` | runs everything in order (about one minute) | |

Python: `<project>/.venv/bin/python` (sympy 1.14, mpmath 1.3).

---

# LaTeX for the supplement

The text below replaces the two `\gapnote` boxes of Section `g3:sec:S:Q0` (the one after the display of
$J$ and the one after Lemma `g3:lem:P`). Notation as in the supplement:
$a=1+\sqrt2$, $b=\sqrt2-1=1/a$, $L=\log a$, $\lambda=\log2$, $\alpha_c=\arctan\frac1{\sqrt2}$, $q=q(s)$.
We also write $\chi_2(z)=\frac12(\Li_2(z)-\Li_2(-z))$. All logarithms and dilogarithms are principal:
$\log$ is analytic on $\CC\setminus(-\infty,0]$ and $\Li_2$ on $\CC\setminus[1,\infty)$, with
$\frac{\dd}{\dd z}\Li_2(z)=-\frac{\log(1-z)}z$ and $\Li_2(\bar z)=\overline{\Li_2(z)}$.
Scripts are in `new_results/QJ_superellipsoid_Q0/`.

```latex
\subsubsection*{Step 2: the closed form of $J$}

\begin{lemma}\label{g3:lem:J}
For $0<q<1$, with $\beta=\sqrt{1+q^2}$,
\[
J(q)=\frac{2\log(1+q)-(1+q)\log(2q)}{1-q^2}-\frac{\operatorname{arsinh}q+\operatorname{arsinh}(1/q)}{\beta}.
\]
\end{lemma}

\begin{proof}
Put $R=\sqrt{q^2+r^2}$ and $u=r+R$. The map $r\mapsto u$ is a $C^1$ increasing bijection of $(0,\infty)$
onto $(q,\infty)$, with $r=\frac{u^2-q^2}{2u}$, $R=\frac{u^2+q^2}{2u}$ and $\dd r=\frac{u^2+q^2}{2u^2}\dd u$.
Since $1+r+R=1+u$, the integrand of $J$ becomes a rational function of $u$:
\[
J(q)=\int_q^\infty g(u)\dd u,\qquad g(u)=-\frac{2(u^2+q^2)}{(1+q)(1+u)(u+q)(u^2+2u-q^2)} .
\]
We have $u^2+2u-q^2=(u+1-\beta)(u+1+\beta)$, so $g$ has the simple poles $-1,-q,-1\pm\beta$. Its partial
fractions are
\begin{multline*}
g(u)=\frac{2}{(q^2-1)(u+1)}-\frac{2q}{(q^2-1)(u+q)}\\
-\frac{2(\beta-1)}{\beta(1+q)(\beta+q-1)}\,\frac1{u+1-\beta}
+\frac{2(\beta+1)}{\beta(1+q)(\beta+1-q)}\,\frac1{u+1+\beta}.
\end{multline*}
The four residues add up to $0$, so $F(u)=\sum_kc_k\log(u-u_k)$ ($c_k$ the residue at $u_k$) satisfies
$F'=g$ and $F(U)=\sum_kc_k\log(1-u_k/U)\to0$ as $U\to\infty$. For $u\ge q$ every $u-u_k$ is positive:
$q+1>0$, $2q>0$, $q+1+\beta>0$, and $q+1-\beta=\frac{2q}{q+1+\beta}>0$ because $(q+1)^2-\beta^2=2q$. So
\begin{multline*}
J(q)=-F(q)=\frac{2\log(1+q)-2q\log(2q)}{1-q^2}\\
+\Bigl(\frac1\beta-\frac1{1+q}\Bigr)\log(q+1-\beta)
-\Bigl(\frac1\beta+\frac1{1+q}\Bigr)\log(q+1+\beta).
\end{multline*}
(The two coefficients are the residues simplified with $\beta^2=1+q^2$.) Now
$\log(q+1-\beta)+\log(q+1+\beta)=\log(2q)$, and $(q+1+\beta)^2=2(q+\beta)(1+\beta)$ gives
$\log\frac{q+1+\beta}{q+1-\beta}=\log\frac{(q+\beta)(1+\beta)}q=\operatorname{arsinh}q+\operatorname{arsinh}(1/q)$.
The coefficient of $\log(2q)$ becomes $-\frac{2q}{1-q^2}-\frac1{1+q}=-\frac{1+q}{1-q^2}$, which gives the claim.
\end{proof}

Every identity in this proof is checked exactly by \path{01_J.py} (the substitution, the partial
fractions, $F'=g$, and the final comparison with the logarithms treated as independent symbols); the
closed form agrees with quadrature to $10^{-79}$ at five values of $q$. Lemma~\ref{g3:lem:step1} and
Lemma~\ref{g3:lem:J} give $Q_0=6(P_1-P_{2a}-P_{2c})$; the set $\{s:q(s)=1\}=\{0,1\}$ has measure zero.

\subsubsection*{Step 3: the three integrals}

The integrands of $P_1,P_{2a},P_{2c}$ are continuous functions of $q\in[\frac1{\sqrt2},1]$ ($P_1$ after
removing the removable singularity at $q=1$), and $q(s)=q(1-s)$. On $(0,\frac12)$ we have
$s=\frac12(1-\sqrt{2q^2-1})$ and $\dd s=-\frac{q\,\dd q}{\sqrt{2q^2-1}}$, so for each of them
\begin{equation}\label{g3:eq:sq}
\int_0^1\Phi(q(s))\dd s=2\int_{1/\sqrt2}^1\Phi(q)\frac{q\,\dd q}{\sqrt{2q^2-1}} .
\end{equation}
We use the following classical identities. For real arguments:
reflection $\Li_2(x)+\Li_2(1-x)=\frac{\pi^2}6-\log x\log(1-x)$ ($0<x<1$);
inversion $\Li_2(x)+\Li_2(1/x)=-\frac{\pi^2}6-\frac12\log^2(-x)$ ($x<0$);
Landen $\Li_2(x)+\Li_2\bigl(\frac x{x-1}\bigr)=-\frac12\log^2(1-x)$ ($x<1$);
duplication $\Li_2(x)+\Li_2(-x)=\frac12\Li_2(x^2)$ ($|x|\le1$);
Abel's relation
\begin{multline*}
\Li_2\Bigl(\frac x{1-y}\Bigr)+\Li_2\Bigl(\frac y{1-x}\Bigr)-\Li_2\Bigl(\frac{xy}{(1-x)(1-y)}\Bigr)\\
-\Li_2(x)-\Li_2(y)=\log(1-x)\log(1-y)
\end{multline*}
on the convex set $x<1$, $y<1$, $x+y<1$ (there all five arguments are $<1$, and the identity holds near
$(0,0)$). For complex arguments: reflection on $\CC\setminus((-\infty,0]\cup[1,\infty))$, inversion
on $\CC\setminus[0,\infty)$, Landen on $\CC\setminus[1,\infty)$ (both sides are analytic there, and these
sets are connected and contain an interval where the real identity holds), duplication for $|z|<1$,
and $\operatorname{Re}\Li_2(e^{i\theta})=\frac{\pi^2}6-\frac{\theta(2\pi-\theta)}4$ for $0\le\theta\le2\pi$.

\paragraph{$P_{2c}$.}
Put $\beta=\sqrt{1+q^2}$ and $c=\cos\alpha_c=\sqrt{2/3}$. Then $q\,\dd q=\beta\,\dd\beta$, $2q^2-1=2\beta^2-3$,
and $\operatorname{arsinh}(1/q)=\log\frac{1+\beta}q=\frac12\log\frac{\beta+1}{\beta-1}$. Substitute
$\beta=\frac{x^2+1}{2cx}$, $x\in[1,\sqrt3]$ (so $x=1\leftrightarrow q=\frac1{\sqrt2}$ and
$x=\sqrt3\leftrightarrow q=1$). Then $\sqrt{2\beta^2-3}=\frac{\sqrt3(x^2-1)}{2x}$,
$\frac{\dd\beta}{\sqrt{2\beta^2-3}}=\frac{\dd x}{\sqrt2\,x}$ and
$\frac{\beta+1}{\beta-1}=\frac{x^2+2cx+1}{x^2-2cx+1}$. By \eqref{g3:eq:sq},
\[
P_{2c}=\frac1{\sqrt2}\int_1^{\sqrt3}\bigl[\log(x^2+2cx+1)-\log(x^2-2cx+1)\bigr]\frac{\dd x}x .
\]
Let $\omega=e^{i\alpha_c}=\frac{\sqrt2+i}{\sqrt3}$. For $x>0$ the points $\pm x\omega$ are not real, so
$\Li_2(\pm x\omega)$ is analytic in $x$, and $\frac{\dd}{\dd x}\Li_2(x\omega)=-\frac{\log(1-x\omega)}x$.
Taking real parts and using $|1\mp x\omega|^2=x^2\mp2cx+1$, the function
$G(x)=\sqrt2\,[\operatorname{Re}\Li_2(x\omega)-\operatorname{Re}\Li_2(-x\omega)]=2\sqrt2\operatorname{Re}\chi_2(x\omega)$
is a $C^1$ antiderivative of the integrand on $[1,\sqrt3]$. Since $\sqrt3\,\omega=\sqrt2+i$,
\[
P_{2c}=2\sqrt2\,\bigl[\operatorname{Re}\chi_2(\sqrt2+i)-\operatorname{Re}\chi_2(e^{i\alpha_c})\bigr].
\]
(a) With $-e^{i\alpha_c}=e^{i(\alpha_c+\pi)}$ and the Fourier value,
$\operatorname{Re}\chi_2(e^{i\alpha_c})=\frac{\pi^2}8-\frac{\pi\alpha_c}4$.
(b) Inversion for $\Li_2$ at $z$ and at $-z$ gives $\chi_2(z)+\chi_2(1/z)=\frac14(\log^2z-\log^2(-z))$ for
$z\notin\RR$. At $z=\sqrt2+i$: $1/z=\frac{\sqrt2-i}3=:t$, $\log z=\frac12\log3+i\alpha_c$,
$\log(-z)=\frac12\log3+i(\alpha_c-\pi)$, so
$\operatorname{Re}\chi_2(\sqrt2+i)+\operatorname{Re}\chi_2(t)=\frac{\pi^2}4-\frac{\pi\alpha_c}2$.
(c) Landen's identity for $\chi_2$, $\chi_2(t)+\chi_2\bigl(\frac{1-t}{1+t}\bigr)=\frac{\pi^2}8+\frac12\log t\,\log\frac{1+t}{1-t}$,
holds for $0<t<1$; both sides are analytic on the connected set $\CC\setminus((-\infty,0]\cup[1,\infty))$,
so it holds there. At $t=\frac{\sqrt2-i}3$ we have $\frac{1-t}{1+t}=b\,e^{i\pi/4}=:y$,
$\log t=-\frac12\log3-i\alpha_c$ and $\log\frac{1+t}{1-t}=\log(a\,e^{-i\pi/4})=L-\frac{i\pi}4$, so
$\operatorname{Re}\chi_2(t)+\operatorname{Re}\chi_2(y)=\frac{\pi^2}8-\frac{L\log3}4-\frac{\pi\alpha_c}8$.
(d) By duplication $\chi_2(y)=\Li_2(y)-\frac14\Li_2(y^2)$, and $y^2=ib^2$, $\operatorname{Re}\Li_2(ib^2)=\frac14\Li_2(-b^4)$
(power series). So $\operatorname{Re}\chi_2(y)=\operatorname{Re}\Li_2(y)-\frac1{16}\Li_2(-b^4)$.
(e) Lemma~\ref{g3:lem:A}: $\operatorname{Re}\Li_2(y)=\frac12\Li_2(b)+\frac{\pi^2}{192}$.
Solving (b)--(e) and inserting (a):
\[
P_{2c}=\sqrt2\Bigl[\Li_2(b)-\tfrac18\Li_2(-b^4)+\tfrac{L\log3}2+\tfrac{\pi^2}{96}-\tfrac{\pi\alpha_c}4\Bigr].
\]

\paragraph{$P_{2a}$.}
Substitute $q=\sinh\varphi$, $\varphi\in[\varphi_0,L]$, $\varphi_0=\operatorname{arsinh}\frac1{\sqrt2}$ (note
$\sinh L=1$). Then $\operatorname{arsinh}q=\varphi$, $\sqrt{1+q^2}=\cosh\varphi$, and by \eqref{g3:eq:sq}
$P_{2a}=2\int_{\varphi_0}^L\varphi\,\frac{\sinh\varphi}{\sqrt{2\sinh^2\varphi-1}}\dd\varphi$.
Let $G(\varphi)=\operatorname{arcosh}(c\cosh\varphi)$, $c=\sqrt{2/3}$. Since $c^2\cosh^2\varphi-1=\frac13(2\sinh^2\varphi-1)$,
$G'(\varphi)=\frac{\sqrt2\sinh\varphi}{\sqrt{2\sinh^2\varphi-1}}$, which is integrable on $[\varphi_0,L]$;
$G$ is continuous there, $G(\varphi_0)=\operatorname{arcosh}1=0$ and $G(L)=\operatorname{arcosh}\frac2{\sqrt3}=\frac12\log3$.
Integrating by parts,
\[
P_{2a}=\sqrt2\Bigl[\frac{L\log3}2-W\Bigr],\qquad W=\int_{\varphi_0}^L\operatorname{arcosh}(c\cosh\varphi)\dd\varphi .
\]
For $\varphi\in[\varphi_0,L]$, $\operatorname{arcosh}(c\cosh\varphi)=\int_{1/\cosh\varphi}^c\frac{\cosh\varphi\,\dd k}{\sqrt{k^2\cosh^2\varphi-1}}$.
The integrand is positive, so by Tonelli we may integrate first in $\varphi\in[\operatorname{arcosh}(1/k),L]$,
for $k\in[\frac1{\sqrt2},c]$. With $\mu=\sqrt{1-k^2}/k$ we have $k^2\cosh^2\varphi-1=k^2(\sinh^2\varphi-\mu^2)$,
the inner integral is $\frac1k\operatorname{arcosh}\frac{\sinh\varphi}\mu\big|_{\sinh\varphi=\mu}^{\sinh\varphi=1}$, and
\[
W=\int_{1/\sqrt2}^c\frac1k\operatorname{arcosh}\frac k{\sqrt{1-k^2}}\dd k .
\]
Now put $\tau=\sqrt{2k^2-1}/k$, i.e. $k=(2-\tau^2)^{-1/2}$, $\tau\in[0,\frac1{\sqrt2}]$. Then
$\frac k{\sqrt{1-k^2}}=\frac1{\sqrt{1-\tau^2}}=\cosh(\operatorname{artanh}\tau)$ and $\frac{\dd k}k=\frac{\tau\,\dd\tau}{2-\tau^2}$, so
$W=\int_0^{1/\sqrt2}\operatorname{artanh}\tau\,\frac{\tau\,\dd\tau}{2-\tau^2}$. Write
$\operatorname{artanh}\tau=\frac12(\log(1+\tau)-\log(1-\tau))$ and $\frac{\tau}{2-\tau^2}=\frac12K(\tau)$ with the odd
function $K(\tau)=\frac1{\sqrt2-\tau}-\frac1{\sqrt2+\tau}$; the substitution $\tau\to-\tau$ in the
$\log(1-\tau)$ part gives
\[
W=\frac14\int_{-1/\sqrt2}^{1/\sqrt2}\log(1+\tau)\Bigl[\frac1{\sqrt2-\tau}-\frac1{\sqrt2+\tau}\Bigr]\dd\tau .
\]
Put $w=1+\tau\in[\frac b{\sqrt2},\frac a{\sqrt2}]$; then $\sqrt2-\tau=a-w$ and $\sqrt2+\tau=w+b$. On this interval
$w>0$, $0<w/a<1$ and $-w/b<0$, and the functions
\begin{align*}
F_-(w)&=-\log w\,\log\Bigl(1-\frac wa\Bigr)-\Li_2\Bigl(\frac wa\Bigr),\\
F_+(w)&=\log w\,\log\Bigl(1+\frac wb\Bigr)+\Li_2\Bigl(-\frac wb\Bigr)
\end{align*}
are real, continuous, and satisfy $F_-'(w)=\frac{\log w}{a-w}$, $F_+'(w)=\frac{\log w}{w+b}$. So
$W=\frac14\bigl([F_-]-[F_+]\bigr)_{b/\sqrt2}^{a/\sqrt2}$. The arguments are $\frac wa\in\{\frac{b^2}{\sqrt2},\frac1{\sqrt2}\}$,
$-\frac wb\in\{-\frac1{\sqrt2},-\frac{a^2}{\sqrt2}\}$, and $1-\frac1{\sqrt2}=\frac b{\sqrt2}$, $1-\frac{b^2}{\sqrt2}=\frac{3b}{\sqrt2}$,
$1+\frac1{\sqrt2}=\frac a{\sqrt2}$, $1+\frac{a^2}{\sqrt2}=\frac{3a}{\sqrt2}$. Collecting the logarithms,
\begin{multline*}
W=\frac{\lambda L}2-\frac{L\log3}2\\
+\frac14\Bigl[\Li_2\Bigl(-\frac1{\sqrt2}\Bigr)-\Li_2\Bigl(\frac1{\sqrt2}\Bigr)
+\Li_2\Bigl(\frac{b^2}{\sqrt2}\Bigr)-\Li_2\Bigl(-\frac{a^2}{\sqrt2}\Bigr)\Bigr].
\end{multline*}
By instances (2), (3) and (15) of Table~\ref{g3:tab:dilog} below, the bracket equals
$-12\Li_2(b)+\Li_2(-b^4)-2L^2+\frac56\pi^2$. Hence
\[
P_{2a}=\sqrt2\Bigl[3\Li_2(b)-\tfrac14\Li_2(-b^4)-\tfrac{\lambda L}2+L\log3+\tfrac{L^2}2-\tfrac{5\pi^2}{24}\Bigr].
\]

\paragraph{$P_1$.}
Let $f(q)=\frac{2\log(1+q)-(1+q)\log(2q)}{1-q^2}$. Substitute $q=\frac{v^2+1}{2\sqrt2\,v}$, $v\in[1,a]$
($q$ increases from $\frac1{\sqrt2}$ to $1$). Then $\sqrt{2q^2-1}=\frac{v^2-1}{2v}$ and
$\frac{q\,\dd q}{\sqrt{2q^2-1}}=\frac{v^2+1}{4v^2}\dd v$, so by \eqref{g3:eq:sq}
$P_1=\int_1^a\frac{v^2+1}{2v^2}f(q(v))\dd v$. Moreover
\[
1+q=\frac{(v+a)(v+b)}{2\sqrt2\,v},\qquad 1-q=-\frac{(v-a)(v-b)}{2\sqrt2\,v},\qquad 2q=\frac{v^2+1}{\sqrt2\,v},
\]
so the integrand is $\sum_\ell C_\ell(v)\log\ell(v)+C_\lambda(v)\lambda$ over the letters
$\ell\in\{v+a,\,v+b,\,v,\,v^2+1\}$, with
\begin{align*}
C_{v+a}=C_{v+b}&=-\frac2{v-a}+\frac2{v+a}+\frac2{v-b}-\frac2{v+b},\\
C_v&=-\frac{\sqrt2}v-\frac2{v+a}+\frac2{v+b},\\
C_{v^2+1}&=\frac{\sqrt2}v+\frac2{v-a}-\frac2{v-b},\\
C_\lambda&=-\frac{\sqrt2}{2v}+\frac2{v-a}-\frac3{v+a}-\frac2{v-b}+\frac3{v+b}.
\end{align*}
The pole at the end point $v=a$ cancels: $-2\log(2a)-2\log(2\sqrt2)+2\log(a^2+1)+2\lambda=0$, since
$a^2+1=2\sqrt2\,a$. So we may replace $\log\ell(v)$ by $\log\ell(v)-\log\ell(a)$ in every term with the
pole $v=a$. The other poles $c\in\{0,-a,b,-b\}$ satisfy $c<1$, so $v-c>0$ on $[1,a]$. We use six
primitives; each is checked by differentiation, under the stated hypotheses, which hold for every term
(\path{04_P1.py} checks them exactly):
\begin{itemize}
\item[T1] $\int\frac{\log(v+p)-\log(a+p)}{v-a}\dd v=-\Li_2\bigl(\frac{a-v}{a+p}\bigr)$, for $v+p>0$, $a+p>0$, $v\le a$ (argument in $[0,1)$);
\item[T2] $\int\frac{\log(v+p)}{v-c}\dd v=\log d\,\log(v-c)-\Li_2\bigl(-\frac{v-c}d\bigr)$, for $d=c+p>0$, $v>c$;
\item[T3] the same integral $=\frac12\log^2(v-c)+\Li_2\bigl(\frac{|d|}{v-c}\bigr)$, for $d<0$, $v-c>|d|$;
\item[T4] the same integral $=\frac12\log^2(v-c)$, for $d=0$;
\item[T5] $\int\frac{\log(v^2+1)}v\dd v=-\frac12\Li_2(-v^2)$;
\item[T6] $\int\frac{\log(v^2+1)-\log(c^2+1)}{v-c}\dd v=-2\operatorname{Re}\Li_2\bigl(\frac{c-v}{c+i}\bigr)$, for real $c$.
\end{itemize}
For T6, $Z(v)=\frac{c-v}{c+i}$ is real only at $v=c$, where $Z=0$; so $Z(v)\notin[1,\infty)$, $\Li_2(Z(v))$ is
analytic in $v$, $\frac{Z'}{Z}=\frac1{v-c}$ is real and $|1-Z|^2=\frac{v^2+1}{c^2+1}$, which gives the
derivative. The complex letter $v^2+1$ meets only the poles $0$, $a$, $b$; T6 is used with $c=a$
(end values $0$ and $\operatorname{Re}\Li_2(Z_1)$) and with $c=b$ (end values $\operatorname{Re}\Li_2(Z_2)$ and $\operatorname{Re}\Li_2(Z_3)$), where
\[
Z_1=\frac{a-1}{a+i}=\frac{1-ib}2,\qquad Z_2=\frac{b-a}{b+i}=\frac{-1+ia}{\sqrt2},\qquad Z_3=\frac{b-1}{b+i}=\frac{-b+i}2 .
\]
These values are real dilogarithms (we use $\tan\frac\pi8=b$, $\tan\frac{3\pi}8=a$):
\begin{itemize}
\item $\operatorname{Re}Z_1=\frac12$ and $1-Z_1=\bar Z_1$, $|Z_1|^2=\frac b{\sqrt2}$, $\arg Z_1=-\frac\pi8$. Reflection gives
$2\operatorname{Re}\Li_2(Z_1)=\frac{\pi^2}6-|\log Z_1|^2$, i.e.
$\operatorname{Re}\Li_2(Z_1)=\frac{29\pi^2}{384}-\frac{L^2}8-\frac{\lambda L}8-\frac{\lambda^2}{32}$.
\item $1-Z_2=a\,e^{-i\pi/4}$, $|Z_2|^2=\sqrt2\,a$, $\arg Z_2=\frac{5\pi}8$. Reflection gives
\begin{multline*}\operatorname{Re}\Li_2(Z_2)=\tfrac{\pi^2}6-\bigl(\tfrac\lambda4+\tfrac L2\bigr)L-\tfrac{5\pi^2}{32}\\-\operatorname{Re}\Li_2(a e^{-i\pi/4}).\end{multline*}
Inversion at $a e^{-i\pi/4}$ ($-ae^{-i\pi/4}=ae^{3i\pi/4}$) and Lemma~\ref{g3:lem:A} give
\[\operatorname{Re}\Li_2(ae^{-i\pi/4})=-\tfrac{L^2}2-\tfrac12\Li_2(b)+\tfrac{7\pi^2}{64},\]
so \[\operatorname{Re}\Li_2(Z_2)=\tfrac12\Li_2(b)-\tfrac{\lambda L}4-\tfrac{19\pi^2}{192}.\]
\item $\frac{Z_3}{Z_3-1}=b\,e^{-i\pi/4}$, $|1-Z_3|^2=\frac a{\sqrt2}$, $\arg(1-Z_3)=-\frac\pi8$. Landen and Lemma~\ref{g3:lem:A} give
\[\operatorname{Re}\Li_2(Z_3)=-\tfrac12\Li_2(b)-\tfrac{L^2}8+\tfrac{\lambda L}8-\tfrac{\lambda^2}{32}+\tfrac{\pi^2}{384}.\]
\end{itemize}
Adding the eighteen terms (listed one by one in \path{04_P1.log}) gives
\begin{multline*}
P_1=6\Li_2(b)-4\Li_2(\tfrac12)+2\Li_2(\tfrac1{\sqrt2})-2\Li_2(-\tfrac1{\sqrt2})+2\Li_2(-\sqrt2)\\
-2\Li_2(-a)-4\Li_2(\tfrac b{\sqrt2})-2\Li_2(\sqrt2\,b)+2\Li_2(-\tfrac b2)+2\Li_2(\tfrac b{2\sqrt2})\\
+\tfrac{\sqrt2}2\Li_2(-1)-\tfrac{\sqrt2}2\Li_2(-a^2)
-\tfrac{\sqrt2}2L^2-\lambda L-\tfrac{\sqrt2}2\lambda L+\tfrac{\lambda^2}4-\tfrac{5\pi^2}{48}.
\end{multline*}
With the values of Table~\ref{g3:tab:dilog} this is
\[
P_1=3\sqrt2\Li_2(b)-\tfrac{\lambda^2}2-\tfrac{\sqrt2}2\lambda L+\bigl(1+\tfrac{3\sqrt2}2\bigr)L^2-\bigl(\tfrac1{48}+\tfrac{3\sqrt2}{16}\bigr)\pi^2 .
\]

\paragraph{The real dilogarithms.}
All arguments above are $S$-units $x$ of $\QQ(\sqrt2)$, $S=\{2,3\}$, with $1-x$ also an $S$-unit; so all
logarithms are integer combinations of $\frac\lambda2$, $\log3$ and $L$. Table~\ref{g3:tab:dilog} lists
the fifteen identity instances we use, in order; each line uses only values found on earlier lines,
$\Li_2(b)$, $\Li_2(-b^4)$, and the classical values $\Li_2(\frac12)=\frac{\pi^2}{12}-\frac{\lambda^2}2$,
$\Li_2(-1)=-\frac{\pi^2}{12}$. We write $\Lambda=\Li_2(b)$ and $M=\Li_2(-b^4)$.

\begin{table}[ht]
\centering\footnotesize\setlength{\tabcolsep}{3pt}
\caption{Reduction of the real dilogarithms ($\Lambda=\Li_2(b)$, $M=\Li_2(-b^4)$).}\label{g3:tab:dilog}
\begin{tabular}{@{}rll@{}}
 & identity, arguments & value obtained \\\hline
1 & inversion, $x=-a^4$ & $\Li_2(-a^4)=-M-8L^2-\frac{\pi^2}6$\\
2 & Landen, $x=-\frac1{\sqrt2}$ ($\frac x{x-1}=b$) & $\Li_2(-\frac1{\sqrt2})=-\Lambda-\frac12(L-\frac\lambda2)^2$\\
3 & duplication, $x=\frac1{\sqrt2}$ & $\Li_2(\frac1{\sqrt2})=\Lambda+\frac{L^2}2-\frac{\lambda L}2-\frac{\lambda^2}8+\frac{\pi^2}{24}$\\
4 & Landen, $x=-a$ ($\frac x{x-1}=\frac1{\sqrt2}$) & $\Li_2(-a)=-\Lambda-L^2-\frac{\pi^2}{24}$\\
5 & inversion, $x=-b$ & $\Li_2(-b)=\Lambda+\frac{L^2}2-\frac{\pi^2}8$\\
6 & duplication, $x=b$ & $\Li_2(b^2)=4\Lambda+L^2-\frac{\pi^2}4$\\
7 & reflection, $x=b^2$ ($1-x=2b$) & $\Li_2(2b)=-4\Lambda-3L^2+2\lambda L+\frac{5\pi^2}{12}$\\
8 & Landen, $x=2b$ ($\frac x{x-1}=-2a$) & $\Li_2(-2a)=4\Lambda+L^2-2\lambda L-\frac{5\pi^2}{12}$\\
9 & Abel, $x=-a^2$, $y=-2a$ & $\Li_2(-a^2)=-6\Lambda-4L^2+\frac{7\pi^2}{24}$\\
10 & reflection, $x=b$ ($1-x=\sqrt2\,b$) & $\Li_2(\sqrt2\,b)=-\Lambda-L^2+\frac{\lambda L}2+\frac{\pi^2}6$\\
11 & Landen, $x=\sqrt2\,b$ ($\frac x{x-1}=-\sqrt2$) & $\Li_2(-\sqrt2)=\Lambda+\frac{L^2}2-\frac{\lambda L}2-\frac{\pi^2}6$\\
12 & Landen, $x=-\frac b2$ ($\frac x{x-1}=b^2$) & $\Li_2(-\frac b2)=-4\Lambda-\frac32L^2+\lambda L-\frac{\lambda^2}2+\frac{\pi^2}4$\\
13 & Abel, $x=-a$, $y=\frac12$ & $\Li_2(\frac b{2\sqrt2})=-6\Lambda-\frac52L^2+\frac32\lambda L-\frac98\lambda^2+\frac{11\pi^2}{24}$\\
14 & Landen, $x=-b$ ($\frac x{x-1}=\frac b{\sqrt2}$) & $\Li_2(\frac b{\sqrt2})=-\Lambda-\frac{L^2}2-\frac{\lambda^2}8+\frac{\pi^2}8$\\
15 & Abel, $x=-a^2$, $y=2b$ & $\Li_2(\frac{b^2}{\sqrt2})-\Li_2(-\frac{a^2}{\sqrt2})=$\\
 & & \quad$M-10\Lambda-L^2-\lambda L+\frac{7\pi^2}8$
\end{tabular}
\end{table}

The Abel instances are admissible: in (9) $x+y=-a^2-2a<1$ and the three other arguments are
$\frac x{1-y}=-1$, $\frac y{1-x}=-\frac1{\sqrt2}$, $\frac{xy}{(1-x)(1-y)}=\frac1{\sqrt2}$, with right side
$\log(1+a^2)\log(1+2a)=2L(L+\frac32\lambda)$; in (13) $x+y=\frac12-a<1$, the arguments are $-2a$,
$\frac b{2\sqrt2}$, $-\frac1{\sqrt2}$, and the right side is $-\lambda(L+\frac\lambda2)$; in (15)
$x+y=-5<1$ (as $1-2b=b^2$ and $1+a^2=2\sqrt2a$), the arguments are $-a^4$, $\frac{b^2}{\sqrt2}$,
$-\frac{a^2}{\sqrt2}$, and the right side is $-2L(L+\frac32\lambda)$. Line (15) is the only place where
$\Li_2(\frac{b^2}{\sqrt2})$ and $\Li_2(-\frac{a^2}{\sqrt2})$ occur, and they occur in $W$ only through
their difference. No value $\Li_2(\frac13)$ is needed.

\paragraph{End of Step 3.}
Inserting the table into $W$ and into $P_1$ gives the closed forms of $P_{2a}$ and $P_1$ above, which
are those of Lemma~\ref{g3:lem:P}; with $P_{2c}$ this proves Lemma~\ref{g3:lem:P}.
```

Then, in "End of the proof of Theorem~\ref{g3:thm:S}", replace "By Lemma~\ref{g3:lem:step1} and the
closed form of $J$" by "By Lemma~\ref{g3:lem:step1} and Lemma~\ref{g3:lem:J}". The addition of the
three closed forms already written there is unchanged (and is checked exactly by `05_reduce.py`).
In the Plan paragraph, "Steps~2 and~3 (evaluation) by computer assistance, described exactly in the box
after Lemma~\ref{g3:lem:P}" can become "Steps~2 and~3 (evaluation) by hand; every substitution,
antiderivative and identity instance is also checked exactly by computer algebra (\path{new_results/QJ_superellipsoid_Q0})".
The paper's sentence after the theorem ("proved, with computer assistance, ...") can drop "with
computer assistance" (the CAS checks are confirmations, not part of the argument).

## Notes for the planner

* Supplement statement of Lemma `g3:lem:P` is unchanged and confirmed exactly.
* Compared with the notes (`13_proof_Q0.py`): our $P_{2a}$ route is shorter ($W=\int_0^{1/\sqrt2}\operatorname{artanh}\tau\,\frac{\tau\,\dd\tau}{2-\tau^2}$
  instead of the notes' form with $\Li_2(\frac23)$; no letter $3$ appears before the last Abel instance),
  and the 54-point / 346-instance elimination is replaced by the 15 instances of the table.
* The table instances (1)-(15) are printed with their exact right sides and 100-digit residuals in
  `05_reduce.log` (the labels there use $r2=\sqrt2$, so $(-3-2r2)=-a^2$, $(-2+2r2)=2b$, etc.).
