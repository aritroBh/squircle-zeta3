# AV2 — Beukers' ζ(3) integrals on squircle domains: report

## 0. Setup and notation

$$D_p=\{x,y\ge0:x^p+y^p\le1\},\qquad B_p=\{x,y,z\ge0:x^p+y^p+z^p\le1\},$$
$$B2(p)=-\tfrac12\iint_{D_p}\frac{\log(xy)}{1-xy}\,dx\,dy,\qquad B3(p)=\iiint_{B_p}\frac{dx\,dy\,dz}{1-(1-xy)z}.$$

Checks (proved): $B2(\infty)=\zeta(3)$ (Beukers), and since $\int_0^1\frac{dz}{1-(1-xy)z}=-\frac{\log(xy)}{1-xy}$, $B3(\infty)=\iint_{[0,1]^2}\frac{-\log(xy)}{1-xy}=2\zeta(3)$. So the claimed normalisation is $B3(\infty)=2\zeta(3)$, not $\zeta(3)$.

**Proved reductions used throughout.**

* (R1) Dirichlet moments: $\iint_{D_p}x^a y^b\,dxdy=\frac1{p^2}\frac{\Gamma(\frac{a+1}p)\Gamma(\frac{b+1}p)}{\Gamma(\frac{a+b+2}p+1)}$. Expanding $1/(1-xy)=\sum(xy)^n$ and differentiating in the exponent gives the exact series
$$B2(p)=\sum_{n\ge1}\frac{\varphi(n\varepsilon)}{n^3},\quad \varepsilon=\tfrac1p,\qquad \varphi(w)=w^3\,\frac{\Gamma(w)^2}{\Gamma(2w+1)}\big(\psi(2w+1)-\psi(w)\big),\ \ \varphi(0)=1.$$
Terms decay like $4^{-n/p}$, so the series is computable to $100+$ digits for all $p\le10^2$ and (as $\varphi(w)\sim\frac{\sqrt\pi}{2}\ln2\,w^{5/2}4^{-w}$) terminates geometrically for any $p$.

* (R2) $p$-polar coordinates $x=r(\cos\theta)^{2/p},\,y=r(\sin\theta)^{2/p}$; doing the $r$-integral and using Landen's identity $\operatorname{Li}_2(c)+\log c\log(1-c)=\frac{\pi^2}6-\operatorname{Li}_2(1-c)$:
$$B2(p)=\frac1p\int_0^{\pi/4}(\cos\theta\sin\theta)^{2/p-1}g(c)\,d\theta,\quad c=(\cos\theta\sin\theta)^{2/p},\quad g(c)=\frac{\operatorname{Li}_2(c)+\log c\,\log(1-c)}{c}.$$
The substitution $t=\tan\theta=e^{-pv/2}$ removes the boundary layer that hides at the endpoints:
$$B2(p)=\frac12\int_0^\infty u^{2/p}g\!\left(u^{2/p}\right)dv,\qquad u=\frac{1}{2\cosh(pv/2)},$$
which tends pointwise to $e^{-v}g(e^{-v})$, giving $B2(\infty)=\frac12\int_0^1 g(c)\,dc=\zeta(3)$. Inverting $c=u^{2/p}$ gives the third form
$$B2(p)=\frac12\int_0^{2^{-2/p}}\frac{g(c)}{\sqrt{1-4c^p}}\,dc.$$

* (R3) For $B3$: integrate $z$ first, then proceed as in (R2):
$$B3(p)=\iint_{D_p}\frac{-\log(1-(1-xy)Z)}{1-xy}\,dxdy,\quad Z=(1-x^p-y^p)^{1/p},$$
$$B3(p)=2\int_0^{2^{-2/p}}\frac{F(c,p)}{\sqrt{1-4c^p}}\,dc,\qquad F(c,p)=\int_0^1 r\,\frac{-\log\big(1-(1-r^2c)(1-r^p)^{1/p}\big)}{1-r^2c}\,dr,$$
with $F(c,\infty)=g(c)/2$ (hence $B3(\infty)=2\zeta(3)$). Numerically we use the peak-subtracted form $B3=2F(c^*)c^*\,\frac{B(1/p,1/2)}p+2\int\frac{F(c)-F(c^*)}{\sqrt{1-4c^p}}dc$ with $c^*=2^{-2/p}$. Also the identity (proved): $B3(p)=2B2(p)+E(p)$, $E(p)=\iint_{D_p}\frac{\log\big(xy/(1-(1-xy)Z)\big)}{1-xy}dxdy$.

## 1. Task 1 — values of $B2(p)$ and $B3(p)$

$B2$ was computed three ways (series R1; R2 v-form; R2 c-form), $B3$ two ways (nested tanh-sinh in the $c$-form at two outer mesh levels, with two different inner grids; earlier cross-checked against the independent square-mapped 2D form and the polar 2D form).

**$B2(p)$, 110-digit computations.** Residuals |series $-$ v-form| (120 dps, step_45): $1.2{\rm e}{-111}$ ($p=1$) to $6.8{\rm e}{-109}$ ($p=96$); |series $-$ c-form| $\sim1{\rm e}{-57}$. 

| p | $B2(p)$ (45 digits) | | p | $B2(p)$ |
|---|---|---|---|---|
| 1 | 0.801371268773062856933158774340966660509990862 | | 16 | 1.19898778233408413459289384106860577200972839 |
| 1.5 | 0.977008724901748547728148406764179328357820245 | | 24 | 1.20067187469952688506555398368144703577126711 |
| 2 | 1.05936833579935324104250490229980343260203366 | | 32 | 1.20127186112860731492248094105458534685376165 |
| 3 | 1.13047216302441973748452688594081533652620718 | | 48 | 1.20170532604138940110845720603655138923048972 |
| 4 | 1.15925940580251098657697122601537056186909813 | | 64 | 1.20185838619021559916506692405955534618233326 |
| 6 | 1.18183729788147160484575802508698300546524082 | | 96 | 1.20196833671646377223906288327484935028108814 |
| 8 | 1.19033032885031205106689324357724655851328124 | | 256 | 1.20204438908545826150839567196842678487610384 |
| 12 | 1.19668336914869755222378156555314197734018898 | | 4096 | 1.20205685414548431279127052047077410856228922 |

**$B3(p)$** (nested-1D tanh-sinh; “digits” = agreement of two outer levels / two inner grids; these are conservative):

| p | $B3(p)$ (45 digits) | digits |
|---|---|---|
| 1 | 0.246617302384923740565168068839592055997186102 | 123 |
| 1.5 | 0.653691227818656823802291538113784838296777177 | 120 |
| 2 | 1.04295375290715279640754254963911132144066641 | 119 |
| 3 | 1.57125186516477580085341788240955003235761282 | 82–95 |
| 4 | 1.85286140923924092547370791577919461117780849 | 62–78 |
| 6 | 2.11146752212026265180823997003908646194291723 | 130 |
| 8 | 2.22171526236031689465391367945661381829721835 | 88 |
| 12 | 2.31244346244700387525547824957110523811172406 | 97 |
| 16 | 2.34847229470482016223514740103697419247670491 | 83 |
| 24 | 2.37689843685618841153873587726478616207814553 | 65 |
| 32 | 2.38783282707398643183188225739441218243449347 | 55 |
| 48 | 2.39627739161807185873916452537538317260730248 | 47 |
| 64 | 2.39946870176394404469322124229337448914970570 | 46 |
| 256 | 2.40375304072207313640466037331344429925296330 | 77 |
| 1024 | 2.40408689560606935868850025202647004088848179 | 51 |
| 4096 | 2.40411185234934951387116467363437252585637404 | 60 |
| 16384 | 2.40411366720276711842143589210663900764285605 | 56 |
| 262144 | 2.40411380564302646005690361738559556012165390 | 52 |
| 1048576 | 2.40411380627278047055850514742926544146997212 | 49 |

So $B3$ is at $\ge80$ digits for $p\in\{1,1.5,2,3,6,8,12,16\}$ (and 62–78 at $p=4$), 45–77 digits at the remaining 11 values, all by $\ge2$ methods/levels. $B3$ is monotone increasing in $p$, approaching $2\zeta(3)=2.4041138063191886$ from below.

## 2. Task 2 — expansions

### 2a. Square end for $B2$: a complete, clean theory (PROVED structure)

Since $B2(p)-\zeta(3)=\sum_{n\ge1}(\varphi(n\varepsilon)-1)n^{-3}$ and $\varphi$ is $C^\infty$, exponentially decaying, with Taylor series $\varphi(w)=1+\sum\varphi_k w^k$ satisfying (exactly, by expanding $\log\Gamma$ and $\psi$ — verified with sympy and to 25 digits numerically)
$$\varphi_1=\varphi_2=0,\quad \varphi_3=-\zeta(3),\quad \varphi_4=\tfrac{9}{4}\zeta(4)=\tfrac{\pi^4}{40},\quad \varphi_5=3\zeta(2)\zeta(3)-9\zeta(5),\quad \varphi_6=\tfrac{79\pi^6}{7560}-4\zeta(3)^2,\ \dots$$
applying Euler–Maclaurin to $h(x)=(\varphi(x\varepsilon)-1)x^{-3}$ gives the asymptotic expansion, valid to all orders:
$$\boxed{\;B2(p)=\zeta(3)+C_2\varepsilon^2+\frac{\zeta(3)}2\varepsilon^3-\frac{3\zeta(4)}{16}\varepsilon^4+0\cdot\varepsilon^5+\frac{\varphi_6}{120}\varepsilon^6+0\cdot\varepsilon^7-\frac{\varphi_8}{252}\varepsilon^8+\cdots\;}$$
i.e. $B2(p)=\zeta(3)+C_2\varepsilon^2+\sum_{k\ge3}\varphi_k\zeta(3-k)\varepsilon^k$, where $\zeta(3-k)$ vanishes for odd $k\ge5$ (so only $\varepsilon^3,\varepsilon^4,\varepsilon^6,\varepsilon^8,\varepsilon^{10},\dots$ appear). **There is no $\varepsilon$ term, no $\varepsilon^2\log\varepsilon$ term, and no logarithm anywhere** — because $\varphi_2=0$ exactly kills the would-be double pole that produces $\varepsilon^2\log\varepsilon$ in this Mellin/EM analysis. The only non-zeta coefficient is
$$C_2=\int_0^\infty\frac{\varphi(w)-1}{w^3}\,dw \;\stackrel{\text{NUMERICAL}}{=}\; -\frac{\zeta(2)}2=-\frac{\pi^2}{12},$$
confirmed two independent ways: the integral above at 81 digits ($|C_2+\zeta(2)/2|=4{\rm e}{-81}$) and the $B2(p)$ fit at 54 digits. A Frullani/Laplace proof of $C_2=-\pi^2/12$ looks within reach ($C_2=\frac12\int_0^\infty \varphi''(w)w^{-1}dw$ since $\varphi'(0)=\varphi''(0)=0$) but was not completed.

Numerical verification (fits on $p=2^8..2^{24}$ at 135 digits, step_04): $c_2$ vs $-\zeta(2)/2$: $1.4{\rm e}{-54}$; $c_3$ vs $\zeta(3)/2$: $1.2{\rm e}{-47}$; $c_4$ vs $-3\zeta(4)/16$: $3.6{\rm e}{-41}$; $c_5=0$: $4{\rm e}{-35}$; $c_6$ vs $(79\pi^6/7560-4\zeta(3)^2)/120$: $2.4{\rm e}{-29}$; $c_7=0$: $6{\rm e}{-24}$; $c_8$ vs $-\varphi_8/252$ with $\varphi_8=-36\zeta(3)\zeta(5)+\pi^2\zeta(3)^2+\frac{2339\pi^8}{604800}$: $7{\rm e}{-19}$; $c_{10}$: $8{\rm e}{-10}$. Note the task's guessed constants $\zeta(5)$ and $\zeta(2)\zeta(3)$ do appear — inside $\varphi_k$ ($k\ge5$) — not $\operatorname{Li}$ at $1/2$.

### 2b. Square end for $B3$: the logarithm that $B2$ dodged (NUMERICAL + mechanism)

$B3$ has a genuinely different expansion. With $E=B3-2B2$ (proved identity) fitted on $p=2^4..2^{20}$ at up to 130 digits (steps_28–30):
$$\boxed{\;B3(p)=2\zeta(3)+2\zeta(2)\,\varepsilon^2\log\varepsilon+A_2\,\varepsilon^2+4\zeta(3)\,\varepsilon^3+O(\varepsilon^4)\;}$$
- coefficient of $\varepsilon^2\log\varepsilon$: $3.28986813370845$ vs $2\zeta(2)=3.28986813369645$ — agreement $1.2{\rm e}{-11}$ (NUMERICAL, very strong). Since $\log\varepsilon<0$, $B3<2\zeta(3)$ ✓.
- $A_2=-5.4189939975$ (two fits agree to $1.2{\rm e}{-9}$): **OPEN constant** (PSLQ against $\{1,\gamma,\zeta(2),\zeta(3),\zeta(2)^2,\zeta(2)\zeta(3),\log2,\zeta(2)\log2,\zeta(3)\log2,\log^22\}$ found nothing).
- $\varepsilon^3$ coefficient $=4.808228$ vs $4\zeta(3)=4.8082276047$ (NUMERICAL, ~5 digits). No $\log$ term detected at $\varepsilon^3$ or $\varepsilon^4$ ($|f_3|,|f_4|<3{\rm e}{-3}$).
- Equivalently $E(p)=[-2\zeta(2)\log p-3.7740599]\varepsilon^2+3\zeta(3)\varepsilon^3+\cdots$ with $e_2=-3.7740599$ OPEN.

Mechanism (supports the $2\zeta(2)$ log term): in $E=2\int_0^{c^*}H(c)\,(1-4c^p)^{-1/2}dc$, $H=F-g/2$, the inner function $F(c,p)$ switches near $c\sim e^{-p}$ from $\int r(-\log(1-(1-r^p)^{1/p}))dr$ to $g(c)/2$; the layer contributes a $\operatorname{Li}_2$-type integral $\int_0^1(\operatorname{Li}_2((1-t)^{1/p})-\zeta(2))\,dt/t$ whose expansion produces exactly $\pm\zeta(2)\log p/p^2$. A fully rigorous uniform-asymptotics derivation was not finished; status NUMERICAL.

**Why does $B3$ get a log when $B2$ doesn't?** $B2=\sum n^{-3}\varphi(n\varepsilon)$ is a pure harmonic sum whose $\varepsilon^2\log\varepsilon$ coefficient is $\varphi_2=0$. $B3$ has no such harmonic-sum form; its correction lives in the boundary layer near the corner $(0,0,1)$ of the superball (where the integrand $1/(1-(1-xy)z)$ is singular), and the layer's thickness is logarithmic in $p$ — hence $\varepsilon^2\log\varepsilon$ survives. Note the dichotomy with Aritro's model (Section VIII–XV of the paper): at the square end the deformation is gentle for $B2$ (power corrections only) but the triple integral $B3$ feels the corner singularity directly and responds with a $\log$.

### 2c. Triangle/simplex end ($p=1$) and the $\tau=1-1/p$ expansion

* **$B2(1)=2\zeta(3)/3$** (NUMERICAL: $3B2(1)=2\zeta(3)$ to 101 digits by series and quadrature independently + PSLQ). Two reductions for a (routine) proof: (i) $B2(1)=\frac14\int_0^1\frac{\zeta(2)-\operatorname{Li}_2(1-s+s^2)}{s(1-s)}ds$, which by Landen becomes Beta/dilog integrals; (ii) the exact series $B2(1)=\sum_{n\ge1}\frac{(H_{2n+2}-H_n)}{(2n+1)(2n+2)\binom{2n}{n}}$, a classical central-binomial-harmonic sum. Status: NUMERICAL (proof reduces to routine but lengthy evaluations; not written out).
* **$B3(1)=0.24661730238492374056516806883959205599718610151827$** (123 digits). **OPEN.** Not in the span of $\{1,\zeta(2),\zeta(3),\zeta'(2),\zeta'(3),\log2,\log3,\zeta(2)\log2,\zeta(3)\log2,\log^22,\log^32,\operatorname{Li}_3(\tfrac12)\}$ or the Clausen family. Decomposition: $B3(1)=P1+P2$ with $P1=\int_0^1\frac{-4\log s\,\arcsin(s/2)}{s\sqrt{4-s^2}}ds=1.0200208006525428$ and $P2=-0.7734034982676190$, both OPEN individually.
* $\tau$-expansion of $B2$: $B2=\frac{2\zeta(3)}3+d_1\tau+d_2\tau^2+\cdots$ with $d_1=-B2'(\varepsilon)|_{\varepsilon=1}=0.524664448852056928181428100965551051111928417616155814312367$ (**OPEN**, exact per-term derivative of the series, geometric convergence; no PSLQ match), $d_2\approx0.0427434684$ (fit, NUMERICAL ~10 digits).
* $\tau$-expansion of $B3$: $d_1=+B3'(p)|_{p=1}=0.7266961$ (exact derivative of the moment series $\sum_m\sum_j(-1)^j\binom{m}{j}J(j,j,m,1)$ to $m=600$ + tail model; NUMERICAL ~6 digits).

### 2d. Special value $p=2$ (the disk)

$B2(2)=1.05936833579935324104250490229980343260203366$. Splitting the Dirichlet series by parity of $n$:
- even-index part $=\zeta(3)/12$ exactly (NUMERICAL, 100 digits);
- odd-index part $=0.959196927202720383925860055507182600038284800017559691920732$ (**OPEN**), which decomposes classically into $\frac\pi8(1+2\log2)+\frac\pi8\big(S_1-S_1'+S_2+2\log2\,S_3\big)$ with the inverse-central-binomial sums $S_2=\sum_{m\ge1}\binom{2m}{m}/(16^m(2m+1)^2)=\operatorname{Cl}_2(\pi/3)-1$ and $S_3=\sum\binom{2m}{m}/(16^m(2m+1))=(\pi-3)/3$ (both NUMERICAL, 106 digits), while $S_1=\sum\binom{2m}{m}H_m/(16^m(2m+1))=0.0502870873772546$ and $S_1'=\sum\binom{2m}{m}H_{2m}/(16^m(2m+1))=0.0743778361833922$ are **OPEN**.

## 3. Task 3 — Legendre-weighted Beukers integrals

$L_n(p)=\iint_{D_p}\frac{-\log(xy)}{1-xy}P_n(x)P_n(y)\,dxdy$ with shifted Legendre $P_n$, computed from exact moment series $M(j,k,p)=-\frac1p\sum_{m\ge0}I(j{+}m,k{+}m,p)\big[\psi(\frac{j+m+1}p)+\psi(\frac{k+m+1}p)-2\psi(\frac{j+k+2m+2}p+1)\big]$ ($m$-tail $4^{-m/p}$; $m\le220p$ used), cross-validated against direct 2D quadrature ($L_2(2),L_3(2)$ agree to 30 digits).

**At $p=\infty$ (Beukers verified, $n\le12$).** Partial fractions give $M(j,j,\infty)=2(\zeta(3)-H_j^{(3)})$ and $M(j,k,\infty)\in\mathbb Q$ for $j\ne k$ (the $\zeta(2)$-pieces cancel), hence immediately $L_n(\infty)=2\beta_n\zeta(3)+(\text{rational})$ with $\beta_n=\sum_k\binom{n}{k}^2\binom{n+k}{k}^2$ the Apéry numbers — the mechanism of Beukers' proof made transparent. Numerically, $d_n^3L_n(\infty)=A_n\zeta(3)+B_n$ with $A_n=2\beta_nd_n^3$ and $B_n\in\mathbb Z$ exactly for all $n\le12$ (e.g. $B_0=0,B_1=-12,B_2=-1404,B_3=-750372,B_4=-137096340$; $n\ge7$ verified to nsum precision). Smallness: $|L_n|^{1/n}\to0.0294=(\sqrt2-1)^4$ (measured $0.0223$ at $n=12$, still converging down).

**At finite $p$: the arithmetic structure does NOT survive.** Testing the ansatz $L_n(p)\stackrel?=2\beta_nB2(p)+(\text{rational})$:
- $L_0(p)=2B2(p)$ exactly at all $p$ ✓ (trivial).
- $r_n=L_n(p)-2\beta_nB2(p)$ is **not rational** at any finite $p$ tested. E.g. $r_1(p)=L_1(p)-10B2(p)$: $-11.9935823215308$ at $p=32$ (moments good to $1{\rm e}{-89}$; the deviation $0.0064$ from $-12$ is real), $-10.6238626362997$ at $p=2$, $-7.93655656786325$ at $p=1$. $d_n^3r_n$ is nowhere near an integer.
- Reason: at finite $p$ the diagonal moments $M(j,j,p)$ are *distinct* $p$-dependent constants (periods of the squircle); the collapse $M(j,j,\infty)=2\zeta(3)-2H_j^{(3)}$ — one $\zeta(3)$ shared by all $j$, with the off-diagonals rational — is a miracle of the square, and it is what lets $\zeta(3)$ factor out with coefficient $2\beta_n$. On the squircle there is no single "ζ(3)"; the whole vector of moments $\{M(j,k,p)\}$ plays that role.

**Decay rates $|L_n(p)|^{1/n}$** (accurate moments, step_43):

| n | p=1 | p=2 | p=4 | p=8 | p=32 | p=∞ |
|---|---|---|---|---|---|---|
| 1 | 0.07716 | 0.03018 | 0.01935 | 0.00438 | 0.01914 | 0.02057 |
| 2 | 0.11246 | 0.17432 | 0.02617 | 0.07903 | 0.02962 | 0.01755 |
| 3 | 0.16850 | 0.21361 | 0.19356 | 0.04933 | 0.09542 | 0.01785 |
| 6 | 0.30148 | 0.39163 | 0.29245 | 0.26979 | 0.19855 | 0.01978 |
| 9 | 0.39624 | 0.47791 | 0.43664 | 0.40247 | 0.36545 | 0.02124 |
| 12 | 0.46632 | 0.51805 | 0.50300 | 0.48769 | 0.38582 | 0.02228 |

The Beukers smallness degrades as $p$ decreases: the rate grows from $(\sqrt2-1)^4\approx0.0294$ at the square toward $\sim0.5$ at the triangle. Mechanism: on $D_p$ with finite $p$, $xy\le2^{-2/p}<1$, so the kernel's $(1,1)$-corner singularity (the engine of Beukers' cancellation) is cut off. So on this path the Legendre-Beukers machine loses both its exact arithmetic form and (much of) its smallness — different from Aritro's squircle path, where $\zeta(3)$ survives at the diamond end.

## 4. Task 4 — plain-words answer

Deforming the domain from the square/cube to the squircle/superball deforms the *weight*, not the kernel. For $B2$, the deformation is remarkably gentle: the exact Dirichlet series $B2(p)=\sum\varphi(n\varepsilon)/n^3$ shows that the entire $p\to\infty$ asymptotics is controlled by the Taylor coefficients of one scalar function $\varphi$; since $\varphi_1=\varphi_2=0$, the corrections are pure powers of $\varepsilon=1/p$ with zeta-value coefficients (one anomalous constant $C_2=-\pi^2/12$), and **no new transcendental enters** — $\zeta(3)$ merely slides down to $2\zeta(3)/3$ at the triangle. For $B3$, the deformation punches the corner singularity at $(0,0,1)$: a logarithmic boundary layer appears, producing $\varepsilon^2\log\varepsilon$ with coefficient $2\zeta(2)$ and at least one new constant $A_2=-5.4189939975$ (OPEN), plus the OPEN simplex value $B3(1)$. On the Legendre side, Beukers' exact $(a_n\zeta(3)+b_n)/d_n^3$ survives only at the square; at finite $p$ it dissolves into the vector of squircle periods $M(j,k,p)$, and the decay rate degrades $\sim16\times$ by the triangle. Net: on this path, $\zeta(3)$ is deformed but not replaced for $B2$; for $B3$ genuinely new constants ($A_2$, $B3(1)$, the $E$-constant $e_2=-3.7740599$) appear.

## 5. Evidence table (constants, statuses)

| # | Statement | Value | Digits / residual | Status |
|---|---|---|---|---|
| 1 | $B2(\infty)=\zeta(3)$, $B3(\infty)=2\zeta(3)$ | — | proved | PROVED |
| 2 | $B2(p)=\sum\varphi(n/p)/n^3$; $\varphi(w)=w^3\Gamma(w)^2\Gamma(2w+1)^{-1}(\psi(2w+1)-\psi(w))$ | — | series vs quadrature residuals $1{\rm e}{-109}$ | PROVED |
| 3 | $B2(p)=\frac12\int_0^\infty u^{2/p}g(u^{2/p})dv=\frac12\int_0^{2^{-2/p}}g(c)/\sqrt{1-4c^p}\,dc$ | — | same | PROVED |
| 4 | $B3(p)=2\int_0^{c^*}F(c,p)/\sqrt{1-4c^p}\,dc$ | — | two outer levels + two inner grids agree (table) | PROVED |
| 5 | $B2(p)=\zeta(3)+C_2\varepsilon^2+\frac{\zeta(3)}2\varepsilon^3-\frac{3\zeta(4)}{16}\varepsilon^4+\frac{\varphi_6}{120}\varepsilon^6-\frac{\varphi_8}{252}\varepsilon^8+\cdots$ | — | coefficients fit to $1{\rm e}{-10}$–$1{\rm e}{-54}$ | PROVED (EM asymptotics) |
| 6 | $\varphi_2=0\Rightarrow$ no $\varepsilon$, no $\varepsilon^2\log\varepsilon$ in $B2$ | — | $c_5=4{\rm e}{-35}$, $c_7=6{\rm e}{-24}$ | PROVED |
| 7 | $C_2=\int_0^\infty(\varphi(w)-1)w^{-3}dw=-\pi^2/12$ | $-0.82246703342411321823620758332301259461$ | $4{\rm e}{-81}$ (integral), $1.4{\rm e}{-54}$ (fit) | NUMERICAL |
| 8 | $B3(p)=2\zeta(3)+2\zeta(2)\varepsilon^2\log\varepsilon+A_2\varepsilon^2+4\zeta(3)\varepsilon^3+\cdots$ | — | log-coefficient agrees $1.2{\rm e}{-11}$ | NUMERICAL |
| 9 | $A_2$ | $-5.4189939975$ | $1.2{\rm e}{-9}$ (two fits) | OPEN |
| 10 | $B3$ $\varepsilon^3$ coefficient $=4\zeta(3)$ | $4.808228$ vs $4.8082276047$ | $3{\rm e}{-5}$ | NUMERICAL |
| 11 | $B2(1)=2\zeta(3)/3$ | $0.80137126877306285693315877434096666051$ | $2{\rm e}{-101}$, PSLQ $[3,0,0,-2]$ | NUMERICAL |
| 12 | $B3(1)$ | $0.24661730238492374056516806883959205599718610$ | $3.2{\rm e}{-123}$ | OPEN |
| 13 | $B2$ $\tau$-derivative $d_1$ | $0.52466444885205692818142810096555105111$ | exact series, 60 digits | OPEN |
| 14 | $B3$ $\tau$-derivative $d_1$ | $0.7266961$ | moment series + tail | NUMERICAL |
| 15 | $B2(2)$ even part $=\zeta(3)/12$ | $0.10017140859663285711664484679262083256$ | $5{\rm e}{-100}$, PSLQ | NUMERICAL |
| 16 | $B2(2)$ odd part | $0.95919692720272038392586005550718260004$ | 105 digits | OPEN |
| 17 | $S_2=\operatorname{Cl}_2(\pi/3)-1$, $S_3=(\pi-3)/3$ | — | $1{\rm e}{-106}$ | NUMERICAL |
| 18 | $d_n^3L_n(\infty)=2\beta_nd_n^3\zeta(3)+B_n$, $B_n\in\mathbb Z$ | $B_1=-12,B_2=-1404,B_3=-750372,\dots$ | exact, $n\le12$ | NUMERICAL (mechanism PROVED) |
| 19 | $|L_n(\infty)|^{1/n}\to(\sqrt2-1)^4$ | $0.02228$ at $n=12$ | trend | NUMERICAL |
| 20 | $r_n=L_n(p)-2\beta_nB2(p)\notin\mathbb Q$ for finite $p$ | $r_1(32)=-11.9935823215308$ | moments $1{\rm e}{-89}$ | NUMERICAL |

## 6. Files produced

- `av2core.py` — all formulas and integrators (Dirichlet series; three 1D forms of $B2$; nested tanh-sinh machinery for $B3$ with corner-stable evaluation, peak subtraction, cached inner grids; moment sums $M(j,k,p)$; custom tensor tanh-sinh).
- `b2_table.json`, `b2_largep.json` — $B2(p)$ at 24 values of $p$ + large-$p$ fit data; `c2.json` — $C_2$.
- `b3_final.json`, `b3_hifi.json`, `b3_table.json` — $B3(p)$ at 19+ values with two-level residuals; `b3run.py`, `b3run2.py` — production drivers.
- `legendre_inf.json`, `legendre_finite.json` — Legendre tables at $p=\infty$ and $p=1,2,4,8,32$.
- Computation scripts (each logs to `progress.log`): `step_02` ($\varphi_k$ discovery, first checks), `step_03` (three-method $B2$ cross-checks, $C_2$), `step_04` ($\varphi_k$ exact via sympy + full eps-expansion fit), `step_10` ($B2$ table), `step_11` ($B3$ table part 1), `step_12`–`step_24` (B3 reduction development: c-form, v-form, split inner, fast cached inner, peak subtraction, convergence studies), `step_26`–`step_30` ($B3$ large-$p$ fits, $E=B3-2B2$), `step_31`–`step_37` (closed-form hunts: $B2(1)$, $B2(2)$ decomposition, $S_2$, $S_3$), `step_38`, `step_39` ($\tau$-expansion of $B2$), `step_40`–`step_43` (Legendre at $\infty$ and finite $p$, quadrature cross-check), `step_44` ($B3(1)=P1+P2$), `step_45` ($B3$ $\tau$-grid; $B2$ residual table), `step_46`, `step_47` (exact $dB3/dp$ at $p=1$). (The harness additionally saved interactive cells as `step_01.py` and `step_48.py`–`step_99.py`.) `summary.md` mirrors this report's results.

Every number above was printed by the cited scripts; nothing is quoted from memory or literature.