
md = r"""# AH3 — Ahmed's Problem 3 at every angle: the length of $x\,u+f(x)^n\,v$

**Setup.** $\gamma_n(x)=x\,u+f(x)^n\,v$, $x\in[0,1]$, $u,v$ unit vectors, turning angle
$\theta$ (interior angle $\alpha=\pi-\theta$), $c=\cos\theta$, $z=\sin^2\frac\theta2$.
The angle between $u$ and $v$ as placed tail-to-tail is $\theta$, so
$|u+sv|=\sqrt{1+2sc+s^2}$; at $\theta=\pi/2$ this is the paper's graph length $\sqrt{1+s^2}$.
$L_n=\int_0^1|u+y'(x)v|\,dx$, $y=f^n$, and $2-L_n=\int_0^1 h_c(y'(x))\,dx$ with the paper's
density $h_c(s)=1+s-\sqrt{1+2cs+s^2}$ (since $\int_0^1(1+y')dx=2$).

## Task 1 — the $\theta$-dependent slope law (proved)

Layer cake on $h_c(s)=\int_0^s\phi_\theta(t)\,dt$ with
$$\phi_\theta(s)=h_c'(s)=1-\frac{s+\cos\theta}{\sqrt{1+2s\cos\theta+s^2}},\qquad s>0,$$
gives, with $m_y(s)=|\{x:y'(x)>s\}|$,
$$2-L_n=\int_0^\infty\phi_\theta(s)\,m_y(s)\,ds=(1-\cos\theta)\,\mathbb E[m_y(S_\theta)].$$
$\phi_\theta$ has mass $1-\cos\theta$ (checked exactly: $h_c(\infty)-h_c(0)=1-c$).
The (proper) slope law $S_\theta$ is identified in four equivalent ways:

* **Tail:** $\mathbb P(S_\theta>s)=\dfrac{\sqrt{1+2cs+s^2}-s-c}{1-c}$ (verified to $2\cdot10^{-59}$).
* **Explicit sampler:** $S_\theta=\dfrac{\cos^2\frac\theta2}{U}-\cos\theta-\sin^2\frac\theta2\,U$,
  $U\sim\mathrm{Unif}(0,1)$; equivalently $S_\theta=\dfrac{\sin(\theta-B)}{\sin B}$ with
  $\tan(B/2)\sim\mathrm{Unif}(0,\tan\frac\theta2)$. Geometrically $\phi_\theta(s)=1-\cos\beta(s)$,
  $\beta(s)$ = angle between the chord $u+sv$ and $v$.
* **Mellin transform = Legendre integral:**
  $\mathbb E[S_\theta^w]=\dfrac{\pi w}{\sin\pi w}\,\dfrac{\Lambda_w(\theta)}{1-\cos\theta}$,
  $-1<w<1$, with $\Lambda_w(\theta)=\int_{\cos\theta}^1P_w(t)\,dt
  =2\sin^2\frac\theta2\,{}_2F_1(-w,1+w;2;\sin^2\frac\theta2)$ (verified to $6\cdot10^{-36}$).
  *Proof:* the paper's first identity is $\int_0^\infty\rho_c(s)s^{w-1}ds
  =\frac{\pi}{\sin\pi w}\Lambda_w(\theta)$ with $\rho_c(s)=\sqrt{1+2cs+s^2}-s-c=\int_s^\infty\phi_\theta$;
  integration by parts turns it into the Mellin transform of $\phi_\theta$.
* **Cumulant generating function:** $\log\mathbb E[S_\theta^w]=\sum_{j\ge1}a_jw^j$ with the
  paper's $a_j$ from $\log R_\theta(\varepsilon)=\log(2z)+\sum a_j\varepsilon^j$: the power-law
  corner coefficients *are* the cumulants of $\log S_\theta$. Geometric mean slope
  $e^{a_1}=(\cos^2\frac\theta2)^{-\cot^2\frac\theta2}/e$ ($=2/e$ at $\theta=\pi/2$; verified
  to $7\cdot10^{-61}$).

## Task 2 — leading term (proved)

The level-set window of the paper's thm:fn proof is $\theta$-independent, so the leading term
carries over with the new slope law:
$$2-L_n=\frac{(1-\cos\theta)\bigl(\log n+\log f'(1)+b(\theta)\bigr)}{f'(1)\,(n-1)}
+O\Bigl(\frac{\log^2n}{n^2}\Bigr),$$
$$c(\theta)=1-\cos\theta=2\sin^2\tfrac\theta2,\qquad
b(\theta)=k_1(\theta)=1+\cot^2\tfrac\theta2\,\log\cos^2\tfrac\theta2.$$
Checks: $\theta=\pi/2$ gives $c=1$, $b=1-\log2$ (paper's thm:fn). For $f(x)=x$ this is the
paper's power-law-corner leading term $2z(\log n-a_1)/(n-1)$ with $b=-a_1$
(and $C(\theta)=-2za_1$ of the paper). At $\theta\to0$ both $c,b\to0$ — no corner, no deficit,
consistent with $2-L_n\equiv0$ for a straight graph.

## Task 3 — cumulants and $\mathcal J_{1,2,3}$ (derived; exact at $f(x)=x$; verified numerically)

Cumulants of $\log(1/S_\theta)$, i.e. $k_j=(-1)^jj!\,a_j$:
$$k_1(\theta)=1+\cot^2\tfrac\theta2\,\log\cos^2\tfrac\theta2,$$
$$k_2(\theta)=1+\frac{\pi^2}{3}-2\operatorname{Li}_2\!\bigl(\sin^2\tfrac\theta2\bigr)
-\cot^4\tfrac\theta2\,\log^2\cos^2\tfrac\theta2,$$
$$k_3(\theta)=2+2\cot^6\tfrac\theta2\,\log^3\cos^2\tfrac\theta2
-12\cot^2\tfrac\theta2\,S_{1,2}\!\bigl(\sin^2\tfrac\theta2\bigr),$$
with $S_{1,2}(z)=\zeta(3)-\operatorname{Li}_3(1-z)+\log(1-z)\operatorname{Li}_2(1-z)
+\tfrac12\log^2(1-z)\log z$. **$\zeta(3)$ enters $k_3$ through $S_{1,2}(\sin^2\frac\theta2)$
with coefficient $-12\cot^2\frac\theta2$, exactly as in the power-law corner's
$a_3=2\cot^2\frac\theta2\,S_{1,2}(z)-\frac13(1+u^3)$** ($k_3=-3!\,a_3$). At
$\theta=\pi/2$: $k_1=1-\log2$, $k_2=1+\pi^2/6$, $k_3=2-\frac32\zeta(3)$ (paper's values;
$k_4=6+\pi^4/15$ also checked).

With $\varepsilon=1/(n-1)$, $\Xi=\log(nf'(1))+k_1(\theta)$, $a=f'(1)$, $A=f''(1)$, $B=f'''(1)$,
$\mu_1=\Xi$, $\mu_2=\Xi^2+k_2$, $\mu_3=\Xi^3+3k_2\Xi+k_3$:
$$2-L_n=\varepsilon\mathcal J_1+\varepsilon^2\mathcal J_2+\varepsilon^3\mathcal J_3
+O(\varepsilon^4\log^4n),$$
$$\mathcal J_1=\frac{2z}{a}\mu_1,\qquad
\mathcal J_2=\frac{z}{a^3}\bigl[A(\mu_2-2\mu_1)-a^2\mu_2+2a(1-z)\bigr],$$
$$\mathcal J_3=\frac{z}{3a^5}\Bigl[a^4\mu_3-(3Aa^2+Ba)(\mu_3-3\mu_2)
+3A^2(\mu_3-4\mu_2+2\mu_1)\Bigr]-\frac{4Az(1-z)}{a^4}-\frac{cz(1-z)}{a^3}.$$
$\zeta(3)$ enters $\mathcal J_3$ only through $k_3$, with coefficient
$\frac{z}{3a^5}(a^4-3Aa^2-Ba+3A^2)\cdot\bigl(-12\cot^2\frac\theta2\bigr)$ on the $\zeta(3)$
letter of $S_{1,2}$. The constant (non-log) terms $2a(1-z)$, $-4Az(1-z)/a^4$,
$-cz(1-z)/a^3$ are the slope-law tail terms: they come from the large-$s$ expansion
$\phi_\theta(s)\sim\sum_{j\ge2}A_js^{-j}$, $A_j=(-1)^j(j-1)\Lambda_{j-1}(\theta)$, i.e. from
the $m\ge1$ terms of Theorem thm:P2. At $\theta=\pi/2$ this reduces to the paper's expansion
($\mathcal J_1=\Xi/a$); at $f(x)=x$ it agrees **symbolically, at general $\theta$, to all
three orders** with the exact Legendre series of thm:P2 (sympy: differences exactly $0$).

## Task 4 — numerics (quadrature at 55–60 digits; $f=x$ cross-checked against the exact
Legendre series of thm:P2, max error $2.2\cdot10^{-48}$ over 12 points)

Values of $2-L_n$ (i.e. $L_n=2-(\cdot)$), expansion through $\varepsilon^3$, residual, and
residual$/(\varepsilon^4\log^4n)$:

### $f(x)=x$
| $\theta$ | $n$ | $2-L_n$ (quadrature) | residual | res$/(\varepsilon^4\log^4n)$ |
|---|---|---|---|---|
| $\pi/3$ | $10^2$ | 0.0233505479024578712763534890130987 | $-1.9833\cdot10^{-7}$ | $-0.04236$ |
| $\pi/3$ | $10^3$ | 0.0035131047714130985185226370634085 | $-7.0475\cdot10^{-11}$ | $-0.03083$ |
| $\pi/3$ | $10^4$ | 0.0004671892753766815868369206232515 | $-1.9229\cdot10^{-14}$ | $-0.02671$ |
| $\pi/3$ | $10^5$ | 0.0000582465480746667584306975678219 | $-4.3506\cdot10^{-18}$ | $-0.02476$ |
| $\pi/2$ | $10^2$ | 0.0483282461489732465745922128270199 | $-4.2407\cdot10^{-7}$ | $-0.09057$ |
| $\pi/2$ | $10^3$ | 0.0071950006295693692906668090698611 | $-1.4881\cdot10^{-10}$ | $-0.06510$ |
| $\pi/2$ | $10^4$ | 0.0009513534551530651244699352654004 | $-4.0306\cdot10^{-14}$ | $-0.05599$ |
| $\pi/2$ | $10^5$ | 0.0001181918973786406655349150085220 | $-9.0677\cdot10^{-18}$ | $-0.05161$ |
| $2\pi/3$ | $10^2$ | 0.0758203405035608833239618400119346 | $-6.8640\cdot10^{-7}$ | $-0.14660$ |
| $2\pi/3$ | $10^3$ | 0.0111369021686870495406280570853439 | $-2.3914\cdot10^{-10}$ | $-0.10461$ |
| $2\pi/3$ | $10^4$ | 0.0014616578211853715580862583277139 | $-6.4285\cdot10^{-14}$ | $-0.08930$ |
| $2\pi/3$ | $10^5$ | 0.0001807532049102043292079408141650 | $-1.4364\cdot10^{-17}$ | $-0.08175$ |

### $f(x)=(e^x-1)/(e-1)$  ($a=A=B=e/(e-1)$)
| $\theta$ | $n$ | $2-L_n$ (quadrature) | residual | res$/(\varepsilon^4\log^4n)$ |
|---|---|---|---|---|
| $\pi/3$ | $10^2$ | 0.0163375368239376000174628000326445 | $-2.1789\cdot10^{-8}$ | $-0.00465$ |
| $\pi/3$ | $10^3$ | 0.0023690934486110753330947916788309 | $-8.6509\cdot10^{-13}$ | $-0.00038$ |
| $\pi/3$ | $10^4$ | 0.0003098829827122085514691657553698 | $6.1074\cdot10^{-16}$ | $0.00085$ |
| $\pi/3$ | $10^5$ | 0.0000382696167394225719903596010328 | $2.3672\cdot10^{-19}$ | $0.00135$ |
| $\pi/2$ | $10^2$ | 0.0337262279127231873049165236750074 | $-4.7653\cdot10^{-8}$ | $-0.01018$ |
| $\pi/2$ | $10^3$ | 0.0048452655222469162219962904583079 | $-1.7969\cdot10^{-12}$ | $-0.00079$ |
| $\pi/2$ | $10^4$ | 0.0006305015448495836791806853012918 | $1.3010\cdot10^{-15}$ | $0.00181$ |
| $\pi/2$ | $10^5$ | 0.0000776131502464841899724249013393 | $4.9726\cdot10^{-19}$ | $0.00283$ |
| $2\pi/3$ | $10^2$ | 0.0527378175407306891043224128589428 | $-7.7204\cdot10^{-8}$ | $-0.01649$ |
| $2\pi/3$ | $10^3$ | 0.0074863686559331208596004596169628 | $-2.5350\cdot10^{-12}$ | $-0.00111$ |
| $2\pi/3$ | $10^4$ | 0.0009676520365125585599366746788522 | $2.1501\cdot10^{-15}$ | $0.00299$ |
| $2\pi/3$ | $10^5$ | 0.0001186103906373137898505924367710 | $7.9901\cdot10^{-19}$ | $0.00455$ |

Residuals scale like $\varepsilon^4\log^4n$ with a small bounded coefficient, as required.

### $f(x)=\sin(\pi x/2)$ — degenerate: $f'(1)=0$ (outside the theorem's hypotheses)
Here $-\log f(1-\tau)\sim\beta\tau^2$ with $\beta=\pi^2/8$: higher contact $m=2$, and the
$\varepsilon$-expansion does **not** apply. The deficit follows the contact law
$$2-L_n\sim(1-\cos\theta)\Bigl(\frac{\log n}{m\,\beta\,n}\Bigr)^{1/m}
=(1-\cos\theta)\frac{2\sqrt{\log n}}{\pi\sqrt n},$$
and more precisely $2-L_n=(1-\cos\theta)\,u^*/\sqrt{\beta n}\,(1+o(1))$ with
$Vu^*e^{-(u^*)^2}=1$, $V=2\sqrt{\beta n}$ (ratio $\to1$; error $O((\log u^*)/u^{*2})$-ish,
numerically). Two independent quadratures agree to $10^{-25}$–$10^{-45}$.

| $\theta$ | $n$ | $2-L_n$ | /leading | /refined $u^*$ |
|---|---|---|---|---|
| $\pi/3$ | $10^2$ | 0.0800741306063264 | 1.1722 | 0.9169 |
| $\pi/3$ | $10^3$ | 0.0310198078089139 | 1.1725 | 0.9684 |
| $\pi/3$ | $10^4$ | 0.0111997147407071 | 1.1594 | 0.9892 |
| $\pi/3$ | $10^5$ | 0.0039056885431476 | 1.1435 | 0.9974 |
| $\pi/2$ | $10^2$ | 0.167451173508204 | 1.2257 | 0.9587 |
| $\pi/2$ | $10^3$ | 0.0637846643271565 | 1.2055 | 0.9957 |
| $\pi/2$ | $10^4$ | 0.0228248979961081 | 1.1814 | 1.0080 |
| $\pi/2$ | $10^5$ | 0.0079221284964806 | 1.1598 | 1.0116 |
| $2\pi/3$ | $10^2$ | 0.265681347474239 | 1.2965 | 1.0141 |
| $2\pi/3$ | $10^3$ | 0.0990702332188297 | 1.2483 | 1.0310 |
| $2\pi/3$ | $10^4$ | 0.0350710883795769 | 1.2102 | 1.0326 |
| $2\pi/3$ | $10^5$ | 0.0121028664923090 | 1.1812 | 1.0303 |

## Task 5 — paragraph and theorem for the paper

See `paper_addition.tex` (ready to paste). Status, honestly:

* Slope law, expectation representation, tail/sampler/Mellin identifications,
  cumulants $k_{1,2,3}$, geometric-mean law: **proved** (the Mellin identity follows from the
  paper's proved thm:P2 first identity by one integration by parts; the sampler is algebra).
* Leading asymptotic (Task 2): **proved** (the paper's thm:fn window argument is
  $\theta$-independent; only the slope average changes).
* $\varepsilon$-expansion with $\mathcal J_{1,2,3}$: **derived** (formal level-set expansion);
  **proved** at $f(x)=x$ at every angle (exact symbolic match with thm:P2); numerically
  verified to residual $O(\varepsilon^4\log^4n)$ for $f=x$ and $f=(e^x-1)/(e-1)$ at three
  angles to $n=10^5$. A full proof needs the paper's thm:fn error analysis pushed one order
  further — same method, not written here.
* Contact law at general $\theta$: **derived + numerical** (the $\theta=\pi/2$ case is the
  paper's proved contact law; the $(1-\cos\theta)$ factor comes from the plateau value
  $h_c\to1-\cos\theta$).

## Evidence table

| check | result | digits / agreement | step files |
|---|---|---|---|
| mass of $\phi_\theta$ = $1-\cos\theta$ | exact | $0$ (algebraic) | step_10, step_35 |
| tail $\mathbb P(S_\theta>s)$ identity | exact | $1.9\cdot10^{-59}$ (25 pairs) | step_11, step_35 |
| Mellin $\mathbb E[S_\theta^w]$ = Legendre | verified | $6.2\cdot10^{-36}$ @ dps 90 (15 pairs) | step_12, step_37, step_38 |
| $\mathbb E\log S_\theta=a_1$ (geometric mean) | verified | $6.7\cdot10^{-61}$ (5 angles) | step_35 |
| $k_{1,2,3}$ closed forms vs numeric cumulants | verified | $5.9\cdot10^{-60}$ (5 angles) | step_14, step_36 |
| right-angle $k_1..k_4$ = paper's values | exact/numeric | $1-\log2$, $1+\pi^2/6$, $2-\frac32\zeta(3)$, $6+\pi^4/15$ | step_14 |
| $\mathcal J_{1,2,3}$ vs exact thm:P2 series ($f=x$, general $\theta$) | **exact** | sympy difference $=0$ | step_17, step_19 |
| quadrature vs exact thm:P2 series ($f=x$) | verified | $\le2.2\cdot10^{-48}$ (12 points, to $n=10^5$) | step_15 |
| $\varepsilon^3$ expansion residual scaling | verified | table above, 24 points | step_22 |
| $\sin$ contact law (2 quadratures) | verified | agree $10^{-25}$–$10^{-45}$; ratios above | step_28, step_29 |

All computations with mpmath at dps 55–90 on darwin; every number above was printed by the
listed `step_NN.py` in this directory (`research/08_theory/AH3_fn_any_angle/`).
"""
open('result.md','w').write(md)
print('result.md written:', len(md), 'chars')
