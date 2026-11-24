# SRF1: Radii of convergence of the superellipsoid's series at both ends

**Work directory:** `<project>/new_results/SRF1_radii_3d/` (all scripts auto-saved as `step_NN.py`; artifacts `w_exact_n3_K12.json`, `SRF1_n3_K60.json`, `EV_coeffs_K30.json`). Every number below was printed by my code. Note: the ε=−1 pole-residue direct continuation (`step_58.py`) still had a residual cancellation bug at corner nodes and is reported as unfinished; the pole is instead established by coefficient asymptotics, which are decisive.

---

## 1. Results for the paper (one paragraph)

**Octahedron end (p=1).** The series of $S(1+d)$ in $d=p-1$, of $S$ in $\tau=1-1/p=d/(1+d)$, and the pure-weight log-series $\log[(2-p)S/(4\sqrt3)]=\sum w_k\tau^k$ all have **radius of convergence ZERO** (status: NUMERICAL, same as the paper's 2D Theorem [thm:zero] proof-sketch status). The coefficients grow factorially, $a_k\sim C\,k!\,A^k k^{\beta}$ with $A=2/\sqrt{\pi^2+\log^2 2}=0.621668\ldots$; the Borel transform has its nearest singularities at $t_{\rm edge}=(-\log2\pm i\pi)/2$, $|t_{\rm edge}|=\sqrt{\pi^2+\log^2 2}/2=1.608575\ldots$, coming from the *edge* configurations of the face simplex (one Gamma variable $\to0$), where $\sum_i e^{2tm_i}=0$ has its nearest root. This is the 3D analogue of the 2D Borel points $\pm i\pi/2$ (control: n=2 reproduces $A=2/\pi$, radius zero, in both $d$ and $\tau$ — the paper states the radius-zero theorem in $\eta$ and $d$; since $\tau=d/(1+d)$ is biholomorphic at 0, the $\tau$-series has radius zero too). The Gamma-prefactor's real poles (at $\tau=n/(n-1)=3/2$, i.e. $p=-2$) and the moment blow-ups (real $\tau\ge 3/4$, and $B_3$ diverging at $\tau=1$) are *not* what stops the series: the obstruction is the complex square-root branch of $B_3(\tau)=\mathbb E[(\frac13\sum e^{2\tau m_i})^{1/2}]$ integrated against the exponential tail of the log-Gamma variables — factorial coefficients with radius zero, exactly the 2D mechanism.

**Cube end (p=∞).** Write $\varepsilon=1/p$. Then (3D analogue of Lemma [lem:holo])
$$S(p)=24+48\,\varepsilon\,E(\varepsilon)+48\,\varepsilon^2\,V(\varepsilon),\qquad \varepsilon=1/p,$$
with $E(\varepsilon)=\int_0^1 a^{2\varepsilon-1}[(1+a)^{-1-2\varepsilon}\sqrt{1+a^{2-2\varepsilon}}-1]\,da$ (the 12 edges) and $V(\varepsilon)=\iint_{[0,1]^2} a^{2\varepsilon-1}c^{\varepsilon-1}[g(a,c)-g(a,0)]\,da\,dc$, $g=N^{-1-2\varepsilon}\sqrt{1+a^{2-2\varepsilon}+(ac)^{2-2\varepsilon}}$, $N=1+a+ac$ (the 8 vertices), both **holomorphic on $|\varepsilon|<1$** (subtract the small-$a$, small-$c$ linear terms — explicit rational pieces in $\varepsilon$; the bracket is then $O(a^2c+ac^{2-2\varepsilon})$, uniformly integrable on $|\varepsilon|\le r<1$; the only singularities the continuation produces are at $|\varepsilon|\ge1$). Hence $S(p)=\sum s_k p^{-k}$ **converges for every $p>1$ and its radius is EXACTLY ONE** (status: radius $\ge1$ PROVED modulo the $V$-bounds lemma, which follows the paper's 2D lemma technique; radius $\le1$ NUMERICAL — $\varepsilon=1$ is singular, since $S(p)$ is not analytic at $p=1$ by the octahedron result). Two new 3D phenomena: (i) $E$ and $V$ have *poles accumulating at* $\varepsilon=1^+$ (at $\varepsilon=k/(k-1)=2,3/2,4/3,\ldots$ from the $a^{(2-2\varepsilon)k}$ terms, and $\varepsilon=2j/(2j-1)$ from the $c$-side) — the "radius-zero at $p=1$" mechanism made visible; (ii) the asymptotics $s_k\approx48R\,(-1)^k$, $|s_k|\to48|R|\approx90$, are **dominated by a simple pole of $V$ at $\varepsilon=-1$** ($p=-1$), $R=\lim(-1)^kV_k\approx-1.87$, while the edge coefficients $E_k$ decay stretched-exponentially like the 2D $M_k$. So at the cube end the nearest singularity in *modulus* is shared by $\varepsilon=\pm1$, but the coefficient asymptotics are governed by $\varepsilon=-1$.

---

## 2. Details and proofs

### Task 1(a): coefficients and radius estimates (n=3)

Computed three ways: (i) the exact D16 engine (sympy, formal $Z_k$) at $K=12$ — reproduces the paper's thm:oct values $w_2..w_6$ **exactly** (diff 0.0, `step_19.py`); (ii) my independent mpmath clone of the D16 composition-sum formulas at $K=60$, dps=220 (`step_23.py`) — worst deviation vs exact $K=12$: **5.5e-194**, and vs H1's independent Route-A 80-order table: worst relative **3.6e-150**; (iii) H1's 150-digit 80-term tables (validated by (ii)).

Root/ratio/Hankel analysis (`step_24,33,34,35,36.py`). The coefficients oscillate (conjugate Borel pair), so plain ratios fail; the Borel-rescaled Hankel determinants $u_k=b_{k+1}^2-b_kb_{k+2}=A^{2k+2}k^{2\beta}R^2\sin^2\theta(1+O(1/k))$ with $b_k=a_k/k!$ give clean monotone estimators, extrapolated by $Y(1+(\beta\log k+c)/k)$:

| series | $A_\infty$ (windows) | predicted $A$ | $\theta_\infty$ | predicted $\theta$ |
|---|---|---|---|---|
| $\tau$-series of $S/(4\sqrt3)$ | 0.62170, 0.62186, 0.62191 | 0.621668 | 1.7885, 1.7898, 1.7886 | 1.787953 |
| $d$-series | 0.62156, 0.62175, 0.62179 | 0.621668 | 1.7896, 1.7877, 1.7885 | 1.787953 |
| log-series $w_k$ | 0.62166, 0.62180 | 0.621668 | 1.7866 | 1.787953 |

The $d$- and $\tau$-series have the *same* rate (the Möbius map $\tau=d/(1+d)$ acts on factorial coefficients as $a_n\mapsto e^{-1/A}a_n(1+o(1))$ — dominated by top coefficients — so it preserves the Gevrey-1 class and the Borel singularity). Since $(|a_k|/k!)^{1/k}\to A>0$, $|a_k|^{1/k}\sim Ak/e\to\infty$: **radius 0**. A model fit $a_k/k!=k^\beta(C_1\cos k\theta+C_2\sin k\theta)$ with $\theta$ fixed at the predicted value has relative residual **1.1e-7** (`step_29.py`), and H1's Padé–Borel ([40/40]) places the poles at $-0.34737\pm1.57606i$ vs predicted $-0.34657\pm1.57080i$ (converging with order).

### Task 1(b): explanation from the representation

$S_3/S_3(1)=2\varepsilon^2\Gamma(1+a)^3/\Gamma(1+2\varepsilon)\cdot B_3(\tau)$, $\varepsilon=1-\tau$, $a=-2\tau/3$. Real-$\tau$ obstructions: $\Gamma(1+a)$ poles at $\tau=3/2,3,\ldots$ ($p=-2,\ldots$); the moment $\mathbb E[G^{c\tau}]$ first fails at $\tau=3/4$; $B_3(\tau)$ diverges at $\tau=1$ ($\int g^{a-\tau/n}dg=\int g^{-\tau}dg$); $\Lambda_3$ alone has radius $1/2$ (log singularity at $\tau=1/2$, i.e. $p=2$ — forced on the log-series by the factor $(2-p)$, not a singularity of $S$). None of these set the radius: the series diverges because for fixed $m=(m_1,m_2,m_3)$ the function $\varphi_m(\tau)=(\frac13\sum e^{2\tau m_i})^{1/2}$ has branch points at the zeros of $\sum e^{2\tau m_i}$; the cheapest large-spread configurations (one $G_i\to0$, cost $e^{-L}$, spread $L$) put the nearest zero at $t\approx(-\log2\pm i\pi)/(2L)$, and the saddle $L\approx k$ gives $a_k\sim k!\,(2/\sqrt{\pi^2+\log^22})^k$ — hence **radius exactly 0** with Borel radius $\sqrt{\pi^2+\log^2 2}/2=1.608575$ (prediction confirmed above to 4 digits).

### Task 1(c): n=2 control

Same Hankel analysis on the 80-term n=2 tables (`step_37.py`): $A_\infty=0.6367$–$0.6369$ vs $2/\pi=0.636620$, $\theta_\infty=1.569$–$1.5726$ vs $\pi/2=1.570796$; $|a_k|^{1/k}/k\to A/e$: **radius zero in both $\tau$ and $d$**. The paper (thm:zero) states radius zero in $\eta$ and $d$; the $\tau$-statement is equivalent (biholomorphic reparametrization) and is what D16's "the series still diverges — the divergence lives entirely in the rational coefficients" records.

### Task 2: cube end

*Holomorphy (the lemma).* $S(p)=24+48\varepsilon E(\varepsilon)+48\varepsilon^2V(\varepsilon)$ (vlib's subtracted form, derived independently in `07_push/S3_superellipsoid/verify/vlib.py`). $E$ extends to $\operatorname{Re}\varepsilon>-1/2$ directly (bracket $=O(a)$ at $a=0$), and to $|\varepsilon|<1$ by subtracting finitely many terms $c_l(\varepsilon)a^l$ (explicit, poles at $\varepsilon=-l/2$ with vanishing residues for $-l/2\in(-1,0)$); its non-trivial singularities are at $\varepsilon=k/(k-1)$ ($k\ge2$, from $a^{2\varepsilon-1}a^{(2-2\varepsilon)k}$, i.e. $\int_0 a^{\mu}da=1/(\mu+1)$) — all $\ge4/3$, accumulating at $1^+$. Same for $V$ (subtract $aA_1+cB_1$ with $A_1=-(1+2\varepsilon)c$, $B_1=-(1+2\varepsilon)a(1+a)^{-2-2\varepsilon}\sqrt{1+a^{2-2\varepsilon}}$; poles from the $c$-side of $N^{-1-2\varepsilon}$ at $\varepsilon=-1,-2,\ldots$ and the $a$/$c$ power families accumulating at $1^+$). So $E,V$ holomorphic on $|\varepsilon|<1$, radius of the $1/p$ series $\ge1$; $\varepsilon=1$ singular (accumulation point; equivalently $S$ not analytic at $p=1$ by Task 1) $\Rightarrow$ **radius exactly 1**.

*Numerics* (`step_47,51,52.py`; flint `arb_series` at every tanh-sinh node, `ctx.cap` bug found and fixed in `step_50-51.py`): $E_2$ (2D $M_k$), $E_3$, $V$ to order 30, h-convergence $4.7{\rm e}{-43}$ ($M_k$), $2.1{\rm e}{-32}$ ($E$), $2.0{\rm e}{-27}$ ($V$). Validated: $s_1,s_2,s_3$ vs paper to 62–65 digits; $s_4$ vs S3's 43-digit value to 4.6e-39. New: $s_5=90.7414583339$, $s_6=-93.2703787987$, $s_7=92.1955840939$, $s_8=-90.7786861648$. Root tests (approach to 1 is the signature of radius 1):

| $k$ | 6 | 10 | 15 | 20 | 25 | 30 |
|---|---|---|---|---|---|---|
| $\|s_k\|^{1/k}$ | 2.130 | 1.567 | 1.348 | 1.251 | 1.197 | 1.161 |
| $\|V_k\|^{1/k}$ | 1.112 | 1.063 | 1.042 | 1.031 | 1.025 | 1.021 |
| $\|M_k\|^{1/k}$ (2D) | 0.309 | 0.469 | 0.538 | 0.552 | 0.629 | 0.656 |
| $\|E_k\|^{1/k}$ | 0.293 | 0.328 | 0.393 | 0.481 | 0.517 | 0.550 |

$V_k=(-1)^k(-1.87+o(1))$ (values $-1.946,1.921,-1.891,1.870,\ldots$ at $k=4..7$): a simple pole at $\varepsilon=-1$ dominates; $E_k$ decays stretched-exponentially (like $M_k\sim e^{-\sqrt{\pi k}}$) from the $\varepsilon=1^+$ pole accumulation. Hence $s_k=48(E_{k-1}+V_{k-2})\approx 48R(-1)^k$ with $|s_k|\to48|R|\approx90$: **radius 1, dominant boundary point $\varepsilon=-1$** — a genuinely 3D feature absent in 2D (where thm:one has $\varepsilon=1$ as the *only* singularity on the circle).

### Evidence table

| claim | status | evidence (residual / agreement) | digits | script |
|---|---|---|---|---|
| $w_2..w_6$ = paper thm:oct | PROVED | exact, diff 0.0 | exact | step_19.py |
| my engine n=3 K=60 vs exact K=12 | NUMERICAL | worst 5.5e-194 | 194 | step_23.py |
| my engine vs H1 Route A (K≤60) | NUMERICAL | worst rel 3.6e-150 | 150 | step_23.py |
| Borel rate $A=2/\sqrt{\pi^2+\log^22}$, radius 0 (n=3) | NUMERICAL | $A_\infty=0.6217$ vs 0.621668; $\theta=1.788$ vs 1.787953; model resid 1.1e-7 | 4 sig. | step_29,33-36.py |
| n=2 control: $A=2/\pi$, radius 0 | NUMERICAL | $A_\infty$ 0.6367 vs 0.636620; $\theta\to\pi/2$ | 4 sig. | step_37.py |
| $s_1..s_4$ vs paper | NUMERICAL | 4.95e-62, 1.86e-65, 3.24e-64, 4.6e-39 | 39–65 | step_47,52.py |
| $s_5..s_{30}$, root test →1 | NUMERICAL | h-convergence ≤2e-27; roots 1.161→1 | 27+ | step_51,52.py |
| $E,V$ holomorphic on $\|\varepsilon\|<1$; radius = 1 | PROVED ≥1 (lemma) / NUMERICAL ≤1 | pole families $\varepsilon=k/(k-1)\to1^+$, $\varepsilon=-s$ | — | step_55-58.py (analysis); residue continuation unfinished |
| $V_k\sim R(-1)^k$, $R\approx-1.87$; pole of $V$ at $\varepsilon=-1$ | NUMERICAL | $\|V_k\|^{1/k}\to1.02\to1$, signs alternate | 2–3 sig. | step_52.py |

**Open items:** the closed form of the $\varepsilon=-1$ residue $R\approx-1.87$ and a bug-free direct continuation of $V$ to $\varepsilon\to-1$ (step_58 draft); the proof that $\varepsilon=1$ is singular for $S$ rests on Task 1's numerical radius-zero (the paper's own 2D thm:zero is at the same proof-sketch status); rigorous $V$-bounds on $|\varepsilon|\le r<1$ follow the 2D lemma's technique but were not written out in full.