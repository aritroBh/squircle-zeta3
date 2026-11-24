# AH3: Ahmed's Problem 3 at every angle — $L_n$ of $\gamma_n(x)=x\,u+f(x)^n\,v$

**Conventions.** Turning angle $\theta$ (interior $\alpha=\pi-\theta$), $c=\cos\theta$, $z=\sin^2\frac\theta2$; $u\cdot v=c$ so $|u+sv|=\sqrt{1+2sc+s^2}$ ($\theta=\pi/2\Rightarrow\sqrt{1+s^2}$, the paper's case). Since $\int_0^1(1+y')\,dx=2$ for $y=f^n$,
$$2-L_n=\int_0^1 h_c(y'(x))\,dx,\qquad h_c(s)=1+s-\sqrt{1+2cs+s^2}.$$

## Task 1 — random-slope representation (PROVED)

Layer cake with $\phi_\theta=h_c'$:
$$2-L_n=\int_0^\infty\phi_\theta(s)\,m_y(s)\,ds=(1-\cos\theta)\,\mathbb E[m_y(S_\theta)],\qquad \phi_\theta(s)=1-\frac{s+\cos\theta}{\sqrt{1+2s\cos\theta+s^2}},$$
a density of mass $1-\cos\theta$. The $\theta$-slope law $S_\theta$ is identified four ways:
- **Tail:** $\mathbb P(S_\theta>s)=\bigl(\sqrt{1+2cs+s^2}-s-c\bigr)/(1-c)$.
- **Sampler:** $S_\theta=\dfrac{\cos^2\frac\theta2}{U}-\cos\theta-\sin^2\frac\theta2\,U$, $U\sim\mathrm{Unif}(0,1)$; equivalently $S_\theta=\sin(\theta-B)/\sin B$ with $\tan(B/2)\sim\mathrm{Unif}(0,\tan\frac\theta2)$, $B$ the angle the chord $u+sv$ makes with $v$ (so $\phi_\theta=1-\cos B$).
- **Mellin:** $\mathbb E[S_\theta^w]=\dfrac{\pi w}{\sin\pi w}\dfrac{\Lambda_w(\theta)}{1-\cos\theta}$, $-1<w<1$, $\Lambda_w(\theta)=\int_{\cos\theta}^1P_w(t)\,dt=2z\,{}_2F_1(-w,1+w;2;z)$. *Proof:* the paper's proved identity $\int_0^\infty\rho_c s^{w-1}ds=\frac{\pi}{\sin\pi w}\Lambda_w$ with $\rho_c=\int_s^\infty\phi_\theta$, integrated by parts.
- **CGF:** $\log\mathbb E[S_\theta^w]=\sum_{j\ge1}a_jw^j$ — the paper's corner coefficients $a_j$ *are* the cumulants of $\log S_\theta$; geometric mean $e^{a_1}=(\cos^2\frac\theta2)^{-\cot^2\frac\theta2}/e$ ($=2/e$ at $\pi/2$).

## Task 2 — leading term (PROVED)

$$2-L_n\sim\frac{c(\theta)\,\bigl(\log n+\log f'(1)+b(\theta)\bigr)}{f'(1)\,n},\qquad c(\theta)=1-\cos\theta,\quad b(\theta)=1+\cot^2\tfrac\theta2\log\cos^2\tfrac\theta2.$$
$\theta=\pi/2$: $c=1$, $b=1-\log2$ ✓ (paper's thm:fn). $f(x)=x$: matches thm:P2asym ($b=-a_1$, $C(\theta)=-2za_1$) ✓. $\theta\to0$: $c,b\to0$, no deficit ✓.

## Task 3 — cumulants and $\mathcal J_k$ (DERIVED; exact at $f=x$; verified)

$k_j=(-1)^jj!\,a_j$ for $\log(1/S_\theta)$:
$$k_1=1+\cot^2\tfrac\theta2\log\cos^2\tfrac\theta2,\qquad k_2=1+\tfrac{\pi^2}{3}-2\mathrm{Li}_2(z)-\cot^4\tfrac\theta2\log^2\cos^2\tfrac\theta2,$$
$$k_3=2+2\cot^6\tfrac\theta2\log^3\cos^2\tfrac\theta2-12\cot^2\tfrac\theta2\,S_{1,2}(z).$$
**$\zeta(3)$ enters $k_3$ through $S_{1,2}(\sin^2\frac\theta2)$ with coefficient $-12\cot^2\frac\theta2$**, exactly as in the corner's $a_3$ ($k_3=-3!a_3$). At $\pi/2$: $1-\log2$, $1+\pi^2/6$, $2-\frac32\zeta(3)$ ✓.

With $\varepsilon=1/(n-1)$, $\Xi=\log(nf'(1))+k_1(\theta)$, $a,A,B=f',f'',f'''(1)$, $\mu_1=\Xi$, $\mu_2=\Xi^2+k_2$, $\mu_3=\Xi^3+3k_2\Xi+k_3$:
$$2-L_n=\varepsilon\mathcal J_1+\varepsilon^2\mathcal J_2+\varepsilon^3\mathcal J_3+O(\varepsilon^4\log^4n),\qquad \mathcal J_1=\frac{2z}{a}\mu_1,$$
$$\mathcal J_2=\frac{z}{a^3}\bigl[A(\mu_2-2\mu_1)-a^2\mu_2+2a(1-z)\bigr],$$
$$\mathcal J_3=\frac{z}{3a^5}\bigl[a^4\mu_3-(3Aa^2+Ba)(\mu_3-3\mu_2)+3A^2(\mu_3-4\mu_2+2\mu_1)\bigr]-\frac{4Az(1-z)}{a^4}-\frac{cz(1-z)}{a^3}.$$
$\zeta(3)$ enters $\mathcal J_3$ only via $k_3$ (coefficient $\frac{z}{3a^5}(a^4-3Aa^2-Ba+3A^2)(-12\cot^2\frac\theta2)$ on the $\zeta(3)$ letter). The constant terms are slope-tail ($m\ge1$ Legendre) terms. At $f=x$ this matches the exact thm:P2 series **symbolically at general $\theta$** (sympy difference $=0$).

## Task 4 — numerics (dps 55–90; all values printed by code)

Quadrature validated against the exact thm:P2 series for $f=x$: max error $\le2.2\cdot10^{-48}$ over 12 $(\theta,n)$ points up to $n=10^5$. Expansion through $\varepsilon^3$: residual scales as $\varepsilon^4\log^4n$ with bounded coefficient:

| $f$ | $\theta$ | $2-L_{100}$ | $2-L_{1000}$ | resid@$10^4$ | resid@$10^5$ | scaled@$10^5$ |
|---|---|---|---|---|---|---|
| $x$ | $\pi/3$ | 0.02335054790245787128 | 0.00351310477141309852 | $-1.92{\cdot}10^{-14}$ | $-4.35{\cdot}10^{-18}$ | $-0.0248$ |
| $x$ | $\pi/2$ | 0.04832824614897324657 | 0.00719500062956936929 | $-4.03{\cdot}10^{-14}$ | $-9.07{\cdot}10^{-18}$ | $-0.0516$ |
| $x$ | $2\pi/3$ | 0.07582034050356088332 | 0.01113690216868704954 | $-6.43{\cdot}10^{-14}$ | $-1.44{\cdot}10^{-17}$ | $-0.0818$ |
| $e^x$ | $\pi/3$ | 0.01633753682393760002 | 0.00236909344861107533 | $6.11{\cdot}10^{-16}$ | $2.37{\cdot}10^{-19}$ | $0.0013$ |
| $e^x$ | $\pi/2$ | 0.03372622791272318730 | 0.00484526552224691622 | $1.30{\cdot}10^{-15}$ | $4.97{\cdot}10^{-19}$ | $0.0028$ |
| $e^x$ | $2\pi/3$ | 0.05273781754073068910 | 0.00748636865593312086 | $2.15{\cdot}10^{-15}$ | $7.99{\cdot}10^{-19}$ | $0.0045$ |

Also validated for generic $f$ ($(x+x^2)/2$, $x+\frac25x^2(x-1)$; scaled residuals $\le0.036$).

$f=\sin(\pi x/2)$ has $f'(1)=0$ (outside the theorem): contact regime $-\log f(1-\tau)\sim\beta\tau^2$, $\beta=\pi^2/8$, and
$$2-L_n\sim(1-\cos\theta)\sqrt{\log n/(2\beta n)}\;,\quad\text{refined: }(1-\cos\theta)\,u^*/\sqrt{\beta n},\ \ Vu^*e^{-(u^*)^2}=1,\ V=2\sqrt{\beta n}.$$
Ratios at $n=10^5$ (two quadratures agreeing to $10^{-25}$–$10^{-45}$): leading law $1.14,1.16,1.18\!\downarrow\!1$; refined $u^*$ law $0.997,1.012,1.030$ at $\theta=\pi/3,\pi/2,2\pi/3$.

## Task 5 — status

`paper_addition.tex` contains the paragraph + `thm:fnangle`. Status: slope law/cumulants/geometric-mean/Mellin **proved** (from the paper's proved thm:P2 identity by parts); leading term **proved** (paper's window argument is $\theta$-independent); $\mathcal J_{1,2,3}$ **derived**, **proved at $f=x$** (exact symbolic match), numerically verified at residual $O(\varepsilon^4\log^4n)$ for three further $f$'s; general-$\theta$ contact law **derived + numerical** (the $(1-\cos\theta)$ factor from the plateau $h_c\to1-\cos\theta$).

## Evidence table

| check | status | agreement | step files |
|---|---|---|---|
| mass $\int\phi_\theta=1-\cos\theta$ | proved | exact ($0$) | step_10, step_35 |
| tail identity | proved | $1.9\cdot10^{-59}$ (25 pairs) | step_11, step_35 |
| Mellin=Legendre | verified | $6.2\cdot10^{-36}$ @dps90 | step_12, step_37, step_38 |
| $\mathbb E\log S_\theta=a_1$ | verified | $6.7\cdot10^{-61}$ | step_35 |
| $k_{1,2,3}$ closed forms | verified | $5.9\cdot10^{-60}$ (5 angles) | step_14, step_36 |
| $\mathcal J_{1,2,3}$ vs exact series ($f=x$) | proved | sympy $=0$ | step_17, step_19 |
| quadrature vs exact series | verified | $\le2.2\cdot10^{-48}$ | step_15 |
| $\varepsilon^3$ residuals (4 $f$'s × 3 $\theta$ × 4 $n$) | verified | tables above | step_22, step_40, step_42 |
| $\sin$ contact law | verified | ratios above | step_28, step_29 |

Deliverables in `new_results/AH3_fn_any_angle/`: `result.md` (full report), `paper_addition.tex`, `results_main.json`, `results_sin.json`, `partA.json`, `step_01.py`–`step_45.py` (every number above was printed by them).