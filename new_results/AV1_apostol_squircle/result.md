# AV1: Apostol's $\zeta(2)$ double integral on the squircle $D_p=\{x,y\ge0:x^p+y^p\le1\}$

**Verdict up front:** Yes — the non-linear squircle deformation produces $\zeta(3)$, in several distinct, precisely located places. At the square end, $\zeta(2)-A(1/\varepsilon)$ is a *pure asymptotic power series* in $\varepsilon=1/p$ (no $\varepsilon^k\log^j\varepsilon$ terms) whose $\varepsilon^k$ coefficient is a multiple zeta value of weight $k$; **$\zeta(3)$ first appears at order $\varepsilon^3$ with proved coefficient $\zeta(3)/6$**, and every later odd order brings richer MZV content (up through $\zeta(3)^3$, $\zeta(3)^2\zeta(5)$, $\zeta(2)\zeta(3)^3$ at orders $9,11$, computed here). At the triangle end $\zeta(3)$ first enters at order $\tau^3$. The quarter-disc has the shockingly simple value $A(2)=7\zeta(2)/12$, and $p=1/2$ evaluates through the golden ratio: $A(1/2)=\pi^2/9-4\log^2\varphi$.

All work is in `steps/step_01.py … step_25.py` (one log line per step in `progress.log`). Status codes: **PROVED** = exact derivation (also verified numerically to the stated digits); **NUMERICAL** = PSLQ/LLL-identified and re-verified at a second precision or by an independent method; **OPEN** = no closed form found.

---

## 0. Two exact representations (everything follows from these)

**(A) Dirichlet series.** Expanding $1/(1-xy)=\sum_{n\ge0}(xy)^n$ and using Dirichlet's integral $\iint_{D_p}x^{a-1}y^{b-1}\,dxdy=\frac{\Gamma(a/p)\Gamma(b/p)}{p^2\Gamma((a+b)/p+1)}$:
$$A(p)=\sum_{m\ge1}\frac{\Gamma(1+m/p)^2}{m^2\,\Gamma(1+2m/p)},\qquad\text{i.e.}\qquad A(1/\varepsilon)=\sum_{m\ge1}\frac{f(m\varepsilon)}{m^2},\quad f(z)=\frac{\Gamma(1+z)^2}{\Gamma(1+2z)}.$$
(For $\iint_{D_p}(xy)^s/(1-xy)$: replace $m$ by $m+s$; for $\iint_{D_p}1/(1-x^2y^2)$: sum over $2n+1$.) Termwise $p\to\infty$ gives $\sum 1/m^2=\zeta(2)$ since $\Gamma(1+z)^2/\Gamma(1+2z)\to1$.

**(B) Corner-layer form of the deficit (exact).** From $A(p)=\int_0^1-\log(1-x(1-x^p)^{1/p})\,dx/x$ and $\zeta(2)=\int_0^1-\log(1-x)\,dx/x$, the substitution $x=e^{-\varepsilon u}$ gives, *exactly*,
$$\zeta(2)-A(1/\varepsilon)=\varepsilon\int_0^\infty F(\varepsilon,u)\,du,\qquad F(\varepsilon,u)=\log\frac{1-e^{\varepsilon\Phi(u)}}{1-e^{-\varepsilon u}},\qquad \Phi(u)=-u+\log(1-e^{-u}).$$
Verified against (A) at $p=2,3,7,25,100$ (agreement $10^{-11}$, quadrature-limited; `step_02`, `step_03`).

---

## Task 1 — Values of $A(p)$, two methods, $\ge100$ digits

`step_01`: 1-D quadrature of $\int_0^1-\log(1-x(1-x^p)^{1/p})/x\,dx$ (mpmath, 110 dps) vs the Dirichlet series (A). Agreement:

| $p$ | $A(p)$ (30 digits) | $p$ | $A(p)$ |
|---|---|---|---|
| 0.5 | 0.170363428923373387660560492838 | 2 | 0.959544872328132087942242180544 |
| 0.75 | 0.373900643565931501765258146052 | 3 | 1.151306019012375699866132050923 |
| 1 | 0.548311355616075478824138388882 | 4 | 1.259975230470078248857102777029 |
| 1.25 | 0.687396980291360689177705342441 | 8 | 1.440769068582192232131588616595 |
| 1.5 | 0.797919249584069508661173205204 | 16 | 1.539785272975062116491982967175 |
| 64 | 1.618055934762283721990570092836 | 256 | 1.638177063389875072671921339765 |

$|\text{quad}-\text{series}|\in[5.7\times10^{-108},\,4.5\times10^{-104}]$ over all twelve $p$; $\zeta(2)-A(256)=-6.757\times10^{-3}$, approaching $\zeta(2)$ from below ($A$ is increasing in $p$ since $D_p$ is nested). **Status: NUMERICAL at 100+ digits, two independent methods.**

---

## Task 2 — Square-end expansion in $\varepsilon=1/p$: $\zeta(3)$ enters at order $\varepsilon^3$

### Form of the expansion: pure powers of $\varepsilon$, no $\varepsilon^k\log^j\varepsilon$

Using $\log\frac{e^z-1}{z}=\frac z2+\sum_{j\ge1}\frac{B_{2j}}{2j(2j)!}z^{2j}$ (Bernoulli parity: *even powers only* beyond the linear term) in representation (B):
$$F(\varepsilon,u)=\log\frac{-\Phi(u)}{u}+\frac{\varepsilon}{2}\log(1-e^{-u})+\sum_{j\ge1}\frac{B_{2j}}{2j(2j)!}\,\varepsilon^{2j}\bigl(\Phi(u)^{2j}-u^{2j}\bigr).$$
Every coefficient function is integrable on $(0,\infty)$ (near $u=0$: powers of $\log u$; near $u=\infty$: $\Phi^{2j}-u^{2j}\sim e^{-u}\mathrm{poly}(u)$), so termwise integration gives a complete asymptotic expansion — **no logarithms of $\varepsilon$ anywhere**. The corner cut-off is visible: the layer curve $e^{-u}+e^{-t}=1$ passes through $(\log2,\log2)$, i.e. the corner is clipped at distance $\sim\varepsilon\log2$ as the prompt anticipates, but the profile is analytic in $\varepsilon$ so only pure powers arise. Numerically this is confirmed by pure-power least-squares fits that succeed to $10^{-15}$–$10^{-24}$ and by $d_4=d_6=d_8=d_{10}=0$ (fit values $2\times10^{-14},6\times10^{-10},\dots$).

### The coefficients (deficit $=\zeta(2)-A(1/\varepsilon)=\sum_{k\ge1}d_k\varepsilon^k$; $d_{2j}=0$ for $j=2,3,4,5$)

$$d_k=\frac{B_{k-1}}{(k-1)(k-1)!}\int_0^\infty\!\!\bigl(\Phi(u)^{k-1}-(-u)^{k-1}\bigr)du\ (k\ge2),\qquad d_1=\int_0^\infty\log\frac{u-\log(1-e^{-u})}{u}\,du .$$

| coeff | value | status |
|---|---|---|
| $d_1$ | $=1.73300259020287278417012270424573195250974925287880319940548\ldots$ (85 digits) | **OPEN** (PSLQ-negative, pools incl. $1,\pi,\gamma,\log2,\log2\pi,\zeta_2,\zeta_3,\log^2 2$, Glaisher, at 90 and 130 dps; `step_19`, `step_20`, `step_21`) |
| $d_2$ | $-\zeta(2)/2$ | **PROVED** ($\int_0^\infty\log(1-e^{-u})=-\zeta_2$; verified $1.4\times10^{-11}$ quad, then exactly via $I_{0,1}$ to $10^{-45}$) |
| $d_3$ | $\boxed{\zeta(3)/6}$ | **PROVED** ($\int(\Phi^2-u^2)=\int(-2uL+L^2)=2\zeta_3+2\zeta_3=4\zeta_3$, times $B_2/(2\cdot2!)=1/24$; **$\zeta(3)$'s first appearance, at order $\varepsilon^3=p^{-3}$**) |
| $d_5$ | $(\zeta(2)\zeta(3)-3\zeta(5))/60$ | NUMERICAL (PSLQ residual $4.4\times10^{-73}$; integral formula proved; `step_05`) |
| $d_7$ | $\zeta(7)/14-\zeta(2)\zeta(5)/42-\zeta(2)^2\zeta(3)/140$ | NUMERICAL (PSLQ, assembly-consistent; `step_07`) |
| $d_9$ | $-\frac{17}{72}\zeta(9)+\frac{3}{40}\zeta(2)\zeta(7)+\frac{79}{1814400}\pi^6\zeta(3)+\frac{1}{1600}\pi^4\zeta(5)-\frac{1}{180}\zeta(3)^3$ | NUMERICAL (PSLQ residual $1.4\times10^{-88}$ + independent residual-Richardson $|{\rm diff}|=8.3\times10^{-18}$; `step_08`, `step_19`) |
| $d_{11}$ | $\frac{31}{22}\zeta(11)-\frac{85}{198}\zeta(2)\zeta(9)-\frac{3}{880}\pi^4\zeta(7)-\frac{79}{332640}\pi^6\zeta(5)-\frac{2339}{119750400}\pi^8\zeta(3)+\frac{1}{11}\zeta(3)^2\zeta(5)-\frac{1}{99}\zeta(2)\zeta(3)^3$ | NUMERICAL (assembled by exact rational arithmetic from the proved $I_{a,b}$ relations, residual $1.5\times10^{-74}$ + independent Richardson $2.0\times10^{-13}$; `step_23`, `step_21`) |

Values: $d_5=-0.018891315252213561013283516964766282$, $d_7=0.0081811775464432526051573580006987106$, $d_9=-0.0083883762048471802267965724837618138$, $d_{11}=0.015174766531794468712069761251378249$.

**Mechanism of the weights:** expanding $\Phi^{k-1}-(-u)^{k-1}$ gives $\int_0^\infty u^a\log^b(1-e^{-u})\,du$ with $a+b=k-1$, and
$$\int_0^\infty u^a\bigl(-\log(1-e^{-u})\bigr)^bdu=a!\!\!\sum_{k_1,\dots,k_b\ge1}\frac{1}{k_1\cdots k_b(k_1+\cdots+k_b)^{a+1}},$$
a Mordell–Tornheim/multiple-zeta sum of weight $a+b+1=k$. Hence the $\varepsilon^k$ coefficient is a **weight-$k$ MZV**; odd weights bring in $\zeta(3)$ and its products. (All $I_{a,b}$ needed were also individually identified — e.g. $I_{2,2}=8\zeta_5-4\zeta_2\zeta_3$, $I_{0,6}=720\zeta_7$, $I_{0,10}=10!\,\zeta_{11}$ — with residuals $10^{-84}$–$10^{-97}$; `step_06`, `step_07`, `step_22`.)

**Global check:** the truncated series through $d_{11}$ reproduces the deficit to $6.3\times10^{-23}$ at $p=40$, $3.2\times10^{-25}$ at $p=60$, $7.7\times10^{-27}$ at $p=80$ (`step_24`). The series is asymptotic (radius 0, coefficients grow like $(2\pi)^{-k}$ from $L(z)$'s $2\pi$-radius); optimal truncation at $\varepsilon=1/2$ leaves error $\sim10^{-5}$, as observed.

**Answer for Task 2: $\zeta(3)$ first appears at order $\varepsilon^3$, coefficient $\zeta(3)/6$ (proved); there are no $\varepsilon^k\log^j\varepsilon$ terms.**

---

## Task 3 — Triangle end ($p=1$): $A(1)=\zeta(2)/3$; $\zeta(3)$ at order $\tau^3$

**Closed form (PROVED, two ways; verified $1.2\times10^{-106}$):**
$$A(1)=\iint_{x+y\le1}\frac{dx\,dy}{1-xy}=\int_0^1\frac{-\log(1-x+x^2)}{x}dx=\int_0^1\frac{\log(1+x)-\log(1+x^3)}{x}dx=\frac{\zeta(2)}2-\frac{\zeta(2)}6=\frac{\zeta(2)}3=\frac{\pi^2}{18},$$
using $1-x+x^2=(1+x^3)/(1+x)$; equivalently $A(1)=\sum_{m\ge1}\frac{1}{m^2\binom{2m}{m}}=\zeta(2)/3$ (classical, from $\arcsin^2$).

**Expansion in $\tau=1-1/p$** ($A$ analytic in $\tau$ at 0): $c_k=\frac{(-1)^k}{k!}\sum_{m\ge1}m^{k-2}f^{(k)}(m)$ with $f(z)=\Gamma(1+z)^2/\Gamma(1+2z)$; via $g=\log f$, $g^{(j)}(m)=2\psi^{(j-1)}(m+1)-2^j\psi^{(j-1)}(2m+1)$:

$$c_0=0.548311355616075478824138388882\;(=\zeta(2)/3),\quad c_1=0.623642263728653966476642755554,$$
$$c_2=0.335712456134004377918919527989,\quad c_3=0.111148100693333219514134820899,$$
$$c_4=0.0241049373555511495627640367923,\quad c_5=0.00294602807596330098449447849793,\quad c_6=-0.000258072152253727774978675619578.$$
(Correct Bell-polynomial values — an early version of `step_11` dropped Bell products; fixed in `step_14`/`step_24`, and independently confirmed by a direct $\tau$-fit of $A(1/(1-\tau))$: $c_1,c_2,c_3$ to 12/9/7 digits, and 7-term series reproduces $A(1.25)$ to $5\times10^{-9}$, `step_17`, `step_25`.) **Status: NUMERICAL, 30 digits.** $c_1,c_2,c_3$ are central-binomial harmonic sums; PSLQ against $\{1,\pi,\log2,\log3,\zeta_2,\zeta_3,\pi/\sqrt3,\mathrm{Cl}_2(\pi/3),\pi\mathrm{Cl}_2,\zeta_2\log2\}$ at $10^{-40}$, maxcoeff $10^8$ found nothing (all apparent hits were spurious — different coefficients at different precisions; `step_12`, `step_13`, `step_14`). Only the plain sums have closed forms (proved): $\sum f(m)=\tfrac13+\tfrac{2\pi}{9\sqrt3}$, $\sum f(m)/m=\tfrac{\pi}{3\sqrt3}$, $\sum mf(m)=\tfrac23+\tfrac{2\pi}{9\sqrt3}$.

**Where $\zeta(3)$ appears:** since $g^{(j)}(m)$ involves $\zeta(j)$, $\zeta(3)$ can first enter $c_3$ (order $\tau^3$), and it does, with a provable nonzero component: $g'''(m)=12\zeta(3)+\{\text{rational}\}$ gives
$$c_3\ \supset\ -\tfrac{1}{6}\cdot 12\zeta(3)\textstyle\sum mf(m)=-2\bigl(\tfrac23+\tfrac{2\pi}{9\sqrt3}\bigr)\zeta(3)=-2.5717603360902256648,$$
the remainder being $W=c_3+2(\tfrac23+\tfrac{2\pi}{9\sqrt3})\zeta(3)=2.68290843678355888427196073142$ (a central-binomial sum, no closed form found; no $\zeta(3)$ detected in it by PSLQ). No $\zeta(3)$ detected at orders $\tau^1,\tau^2$ (PSLQ with $\zeta(3)$ in pool). **Answer: $\zeta(3)$ first appears at order $\tau^3$, with proved coefficient component $-(4/3+4\pi/(9\sqrt3))\zeta(3)$, mixed with non-elementary sums.**

---

## Task 4 — Special points

| point | value | status / proof |
|---|---|---|
| $p=1$ | $\zeta(2)/3=\pi^2/18$ | PROVED (above; residual $1.2\times10^{-106}$) |
| $p=2$ (quarter disc) | $\boxed{7\zeta(2)/12=7\pi^2/72}$ | PROVED: even/odd split of series (A): even part $\frac14\sum\frac{1}{k^2\binom{2k}{k}}=\frac{\zeta(2)}{12}$; odd part $\pi\sum_k\frac{\binom{2k}{k}}{16^k}\frac{2k}{(2k-1)^2}=\frac{\pi^2}{12}$, using $\sum_k\frac{\binom{2k}{k}}{16^k(2k-1)}=\frac{2-\sqrt3}{2}$ and $\sum_k\frac{\binom{2k}{k}}{16^k(2k-1)^2}=\frac{\pi}{12}-\frac{2-\sqrt3}{2}$ (elementary integrals: $\int_0^{1/2}\frac{(1-u^2)^{-1/2}-1}{u^2}du=\tan\frac{\pi}{12}$, and integration by parts + $\tan\varphi\cot2\varphi=(1-\tan^2\varphi)/2$). All three verified $10^{-96}$–$10^{-102}$ (`step_15`). Residual $7.2\times10^{-106}$. |
| $p=1/2$ | $\boxed{2\zeta(2)/3-4\log^2\varphi=\pi^2/9-4\log^2\varphi}$, $\varphi=(1+\sqrt5)/2$ | PROVED: $x=u^2,y=v^2$ gives $A(1/2)=4\iint_{\rm tri}\frac{uv}{1-u^2v^2}=2A(1)-2C$ with $C=\iint_{\rm tri}\frac{dudv}{1+uv}=\int_0^1\frac{\log(1+x-x^2)}{x}dx$; factoring $1+x-x^2=(1-x/\varphi)(1+\varphi x)$ gives $C=-\operatorname{Li}_2(1/\varphi)-\operatorname{Li}_2(-\varphi)=2\log^2\varphi$ by the classical golden-ratio dilog identities $\operatorname{Li}_2(1/\varphi)=\pi^2/10-\log^2\varphi$, $\operatorname{Li}_2(-\varphi)=-\pi^2/10-\log^2\varphi$. Verified: $A(1/2)$ residual $4.7\times10^{-97}$; $C=2\log^2\varphi$ residual $1.1\times10^{-101}$ (`step_15`). |
| $p=3/2$ | $0.797919249584069508661173205204$ | OPEN (PSLQ-negative vs $\{1,\pi,\log2,\zeta_2,\zeta_3,\pi/\sqrt3,\Gamma(1/3)^3/\sqrt3\}$) |

---

## Task 5 — Variants

### (a) $A_s(p)=\iint_{D_p}(xy)^s/(1-xy)\,dxdy$  (→ $\psi'(1+s)=\zeta(2,1{+}s)$ as $p\to\infty$)

Same corner-layer analysis with $J_s(u)=\int_0^u t^s/(1-t)dt=C_s+\log(1/\delta)+s\delta-\frac{s(s-1)}4\delta^2+\cdots$ near $u=1$ (PROVED formulas, numerically verified):
$$\zeta(2,1{+}s)-A_s(1/\varepsilon)=d_1\varepsilon-(s+\tfrac12)\zeta(2)\,\varepsilon^2+\frac{6s^2+6s+1}{6}\,\zeta(3)\,\varepsilon^3-\tfrac38\,s(s+1)(2s+1)\,\zeta(4)\,\varepsilon^4+O(\varepsilon^5).$$
Same $d_1$ for all $s$; **$\zeta(3)$ at order $\varepsilon^3$ with polynomial coefficient $(6s^2+6s+1)/6$** (for $s=0$: $\zeta_3/6$ ✓). Richardson extraction at $s=\pm\tfrac12,1,2$ (`step_16`): $d_2$ to 9 digits, $d_3$ to 6–8 digits, and the *predicted zeros* at $s=-1/2$ for both the $\varepsilon^2$ and $\varepsilon^4$ coefficients confirmed ($6.9\times10^{-12}$, $2.3\times10^{-7}$). At $p=1$, integer $s\ge1$: $A_s(1)=\frac{\zeta(2)}3-\sum_{j=0}^{s-1}\frac{1}{(j+1)^2\binom{2j+2}{j+1}}$ (telescoping $xy/(1-xy)=1/(1-xy)-1$): $A_1(1)=\zeta(2)/3-1/2$, $A_2(1)=\zeta(2)/3-13/24$ — PROVED, residuals $10^{-100}$+.

### (b) $B(p)=\iint_{D_p}1/(1-x^2y^2)\,dxdy$  (→ $\sum(2n+1)^{-2}=3\zeta(2)/4=\pi^2/8$)

Here $\log(1+e^{-x})+\frac x2=\log(2\cosh\frac x2)$ is **even**, so the $\varepsilon$-linear pieces of the two logs cancel and *all even powers of $\varepsilon$ vanish identically* (PROVED + verified):
$$\tfrac34\zeta(2)-B(1/\varepsilon)=\frac{d_1}{2}\,\varepsilon-\frac{\zeta(3)}{6}\,\varepsilon^3+\frac{7}{60}\bigl(3\zeta(5)-\zeta(2)\zeta(3)\bigr)\varepsilon^5+O(\varepsilon^7).$$
Verified by degree-8 Richardson (`step_16`, `step_17`): $d_1/2$ to 10 digits; $\varepsilon^2$ coefficient $3.9\times10^{-15}$; $-\zeta(3)/6$ to 9 digits; $\varepsilon^4$ coefficient $9\times10^{-10}$; $\varepsilon^5$ coefficient $0.132239017$ vs predicted $0.132239206765$ (fit-limited, consistent). **$\zeta(3)$ at order $\varepsilon^3$ again, coefficient $-\zeta(3)/6$; odd powers and odd zeta values only — the purest display of the mechanism.**

Special values: $B(2)=\dfrac{\pi}{4}\sum_{n\ge0}\dfrac{\binom{2n}{n}}{16^n(2n+1)}=\dfrac{\pi}{4}\cdot\dfrac{\pi}{3}=\boxed{\pi^2/12}$ PROVED (residual $2.4\times10^{-106}$); $B(1)=\frac12A(1)+\frac12C=\boxed{\zeta(2)/6+\log^2\varphi}$ PROVED (residual $5.8\times10^{-108}$); $B(1/2)=0.166787133798766760630220405622$ NUMERICAL/OPEN.

---

## Task 6 — The mechanism, in plain words

1. **Everything happens in the corner.** For $\varepsilon=1/p\to0$, $D_p$ coincides with the unit square up to $e^{-\mathrm{const}/\varepsilon}$ except in a layer of width $\sim\varepsilon$ around $(1,1)$, where Apostol's kernel $1/(1-xy)$ is singular. The whole deficit $\zeta(2)-A(p)$ lives there.
2. **The non-linear deformation leaves a *profile*.** Under $x=e^{-\varepsilon u}$ the squircle boundary becomes the *fixed* convex curve $e^{-u}+e^{-t}=1$ (through $(\log2,\log2)$ — the "distance $\varepsilon\log2$"), encoded by $\Phi(u)=-u+\log(1-e^{-u})$. A **linear** homotopy would give a profile without the $\log(1-e^{-u})$ part; it is exactly this logarithm that produces $\int_0^\infty-u\log(1-e^{-u})\,du=\zeta(3)$, the classical trilog integral. This is the quantitative content of "non-linear homotopy breaks parity and lets $\zeta(3)$ in."
3. **The expansion is pure power, no $\varepsilon^k\log^j\varepsilon$:** the layer cross-section is analytic in $\varepsilon$ and all coefficient integrals converge; the would-be singular endpoint behaviors are integrable.
4. **Parity survives as zero coefficients, not as exclusion of $\zeta(3)$:** $\log\frac{e^z-1}{z}-\frac z2$ is even (odd Bernoulli numbers vanish), so $d_4=d_6=d_8=d_{10}=0$; and the $\varepsilon^k$ coefficient is a weight-$k$ MZV — $\zeta(3)$, the first odd zeta, enters at the first open slot, $\varepsilon^3$, with coefficient $\zeta(3)/6$. Weights then escalate: $\{\zeta_5,\zeta_2\zeta_3\}$ at $\varepsilon^5$; $\{\zeta_7,\zeta_2\zeta_5,\zeta_2^2\zeta_3\}$ at $\varepsilon^7$; $\{\zeta_9,\zeta_2\zeta_7,\pi^6\zeta_3,\pi^4\zeta_5,\zeta_3^3\}$ at $\varepsilon^9$; $\{\zeta_{11},\dots,\zeta_3^2\zeta_5,\zeta_2\zeta_3^3\}$ at $\varepsilon^{11}$.
5. **Comparison with Aritro's squircle-perimeter story.** Same two pillars. (i) *Corner layer:* as in the perimeter problem, the deformation is felt only in a corner boundary layer; there it produced corner coefficients $M_k$, here it produces the profile integrals $\int u^a\log^b(1-e^{-u})$. (ii) *Same $\Gamma$-source of $\zeta(3)$:* even the area has it — $\log\operatorname{area}(D_p)=2\log\Gamma(1+\varepsilon)-\log\Gamma(1+2\varepsilon)=\sum_{k\ge2}(-1)^k(2-2^k)\zeta(k)\varepsilon^k/k$, whose $\varepsilon^3$ coefficient is $2\zeta(3)$; the area itself is $1-\zeta(2)\varepsilon^2+2\zeta(3)\varepsilon^3+\cdots$. And at the triangle end, $\zeta(3)$ enters at $\tau^3$ via $\psi''$(integers) — again "third derivative of $\log\Gamma$." The contrast with Aritro's Sections VI–VII is sharp: on the Fourier side, parity kills the $1/n^3$ term outright; on this non-linear path, parity only *relocates* — it forces $d_{2j}=0$ ($j\ge2$) and, for the $1/(1-x^2y^2)$ variant, kills every even power — while the full odd-zeta tower flows through the open odd-weight channels.

---

## Evidence table (all numbers printed by the stated step files; two-method or two-precision verification throughout)

| Item | Value / identity | Residual & digits | Status | File |
|---|---|---|---|---|
| $A(p)$, 12 values | table above | quad vs series: $10^{-108}$–$10^{-104}$, 100+ digits | NUMERICAL | step_01 |
| deficit $=\varepsilon\int F$ | exact identity | vs series: $10^{-11}$–$10^{-13}$ (quad-limited) | PROVED | step_02, step_03 |
| $d_1$ | $1.73300259020287278417012270424573195250974925287880319940548\ldots$ | quad 85 digits; Richardson match 24 digits; PSLQ-negative 2 precisions | OPEN | step_19, step_20, step_21 |
| $d_2=-\zeta_2/2$, $d_3=\zeta_3/6$ | proved integrals | series-fit matches $10^{-12}$; exact $I_{a,b}$ checks $10^{-45}$ | PROVED | step_03, step_04, step_05 |
| $d_4=d_6=0$ | fits: $2\times10^{-14}$, $6\times10^{-10}$ | — | PROVED (Bernoulli parity) | step_03, step_04 |
| $d_5=(\zeta_2\zeta_3-3\zeta_5)/60$ | $-0.018891315252213561013283516964766282$ | PSLQ residual $4.4\times10^{-73}$ | NUMERICAL | step_05, step_06 |
| $d_7=\zeta_7/14-\zeta_2\zeta_5/42-\zeta_2^2\zeta_3/140$ | $0.0081811775464432526051573580006987106$ | PSLQ + assembly, individuals resid. $10^{-94}$+ | NUMERICAL | step_07 |
| $d_9$ (formula above) | $-0.0083883762048471802267965724837618138$ | PSLQ $1.4\times10^{-88}$; resid.-Richardson $8.3\times10^{-18}$; two forms $1.3\times10^{-113}$ | NUMERICAL | step_08, step_19, step_24 |
| $d_{11}$ (formula above) | $0.015174766531794468712069761251378249$ | assembly resid. $1.5\times10^{-74}$; Richardson $2.0\times10^{-13}$ | NUMERICAL | step_21, step_22, step_23 |
| series vs deficit, $p=40,60,80$ | through $d_{11}$ | $6.3\times10^{-23},3.2\times10^{-25},7.7\times10^{-27}$ | — | step_24 |
| $A(1)=\zeta_2/3$ | $0.548311355616075478824138388882$ | $1.2\times10^{-106}$ | PROVED | step_01, step_24 |
| $c_1..c_6$ (triangle) | values above | direct $\tau$-fit confirms 12/9/7 digits; 7-term series vs $A(1.25)$: $5.2\times10^{-9}$ | NUMERICAL | step_14, step_17, step_24, step_25 |
| $c_3$'s $\zeta_3$ part | $-2(\frac23+\frac{2\pi}{9\sqrt3})\zeta_3=-2.5717603360902256648$; $W=2.68290843678355888427196073142$ | recombination exact by construction | PROVED component | step_14 |
| $A(2)=7\zeta_2/12$ | $0.959544872328132087942242180544$ | $7.2\times10^{-106}$; binomial lemmas $10^{-102}$ | PROVED | step_15 |
| $A(1/2)=2\zeta_2/3-4\log^2\varphi$ | $0.170363428923373387660560492838$ | $4.7\times10^{-97}$; $C=2\log^2\varphi$: $1.1\times10^{-101}$ | PROVED | step_15 |
| $A(3/2)$ | $0.797919249584069508661173205204$ | PSLQ-negative | OPEN | step_15 |
| $A_s$ expansion ($d_2,d_3,d_4$ in $s$) | formulas above | 9/8/3-4 digits + exact zeros at $s=-1/2$ | PROVED + verified | step_16 |
| $B$ expansion | $\frac{d_1}2\varepsilon-\frac{\zeta_3}6\varepsilon^3+\frac7{60}(3\zeta_5-\zeta_2\zeta_3)\varepsilon^5$ | 10 / 9 / 7 digits; even coefficients $\le10^{-9}$ | PROVED + verified | step_16, step_17 |
| $B(2)=\pi^2/12$, $B(1)=\zeta_2/6+\log^2\varphi$ | $0.822467033424113218236207583323$, $0.505720498385232131908998265673$ | $2.4\times10^{-106}$, $5.8\times10^{-108}$ | PROVED | step_15, step_24 |
| $A_1(1)=\zeta_2/3-1/2$, $A_2(1)=\zeta_2/3-13/24$ | $0.048311355616075478824138388882$, $0.00664468894940881215747172221534$ | $10^{-100}$+ | PROVED | step_15, step_17 |

**Files:** `steps/step_01.py` (Task 1), `step_02–03` (corner-layer identity), `step_04–10` (square-end coefficients $d_1..d_9$, $I_{a,b}$ machinery), `step_11–14` (triangle end; note `step_11`'s $c_{2..6}$ superseded by `step_14/24` after a Bell-polynomial bug), `step_15` (special points), `step_16–17` (variants), `step_18–19` (diagnosis + corrected 85-digit $d_1$; the `step_04` 45-digit $d_1$ was accurate only to $10^{-11}$), `step_20–23` ($d_1$ PSLQ, $d_{11}$), `step_24` (consolidation/evidence), `step_25` (close-out); `progress.log` has one line per step.