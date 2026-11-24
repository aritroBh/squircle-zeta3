# DCYA1 — Decay law of $M_k(\alpha)$ at a general corner angle

**Target (paper, "Corners of Any Angle").** With $N=k+1$, $\theta=\sqrt{2\alpha N}$, $A(\alpha)=\sqrt{2\alpha\sin\alpha}$:
$$\Bigl|M_k(\alpha)+\tfrac{A(\alpha)}{N}e^{-\theta}\cos\theta\Bigr|\le C(\alpha)\,N^{-3/2}e^{-\theta}\qquad(N\ge N_0(\alpha)).$$

## Lemma ledger

| # | Statement | Status | Where |
|---|---|---|---|
| L1 | Exact representation $M(u;\alpha)=T(u)+R(u)$ | **PROVED (by hand)** + CERTIFIED-numerical ($\sim$50 digits) | §1; verify_dcya1.py |
| L2 | Reliable $M_k$ evaluator (custom series + tanh-sinh, high dps); reproduces K2's $M_{100}(\pi/2)$ | NUMERICAL | mk_eval.py; §5 |
| L3 | Weight expansion $M_k=\sum_{m\ge1}(-1)^{m-1}T_k^{(m)}+R_k$ | **PROVED (by hand)** | §2 |
| L4 | Main term = Hankel at $z=\pm i\alpha$ ($m=1$), reduced to 1-D integrals $I_+^{(j)}$; leading $+g_1$ captures $\kappa_1$ | **PROVED (by hand)** + NUMERICAL | §3 |
| L5 | Saddle model $MV$ (linear phase $\omega$) approximates main term; saddle series $F_1,F_2$ via Gaussian moments; $MV\to$ lead | **PROVED (structure)** + NUMERICAL | §4 |
| L6 | Main-term bound $|M_k^{\rm main}-{\rm lead}|\le C_1 N^{-3/2}e^{-\theta}$ | NUMERICAL (structure validated; constant not arb-certified) | §4, **GAP** |
| L7 | Rest bound $|M_k-M_k^{\rm main}|\le C_2 e^{-\delta\sqrt N}$, $\delta>\sqrt{2\alpha}$ | OPEN (machinery derived; constant not assembled) | §6, **GAP** |
| L8 | Residual ratios $\rho=|M_k-{\rm lead}|/(N^{-3/2}e^{-\theta})$ are $O(1)$ at all target angles | NUMERICAL ($k\le100$) | §5 |

## 1. Exact representation (L1) — PROVED

With $L(x)=\log(1+e^{-x})$ and $\widetilde B_\alpha(z)=\dfrac{e^{z}-2\cos\alpha}{1+\sqrt{1-2\cos\alpha\,e^{z}+e^{2z}}}$:
$$M(u;\alpha)=\underbrace{\int_0^\infty e^{-(1+u)L(x)}e^{-x}\widetilde B_\alpha((u-1)x)\,dx}_{T(u)}
+\underbrace{\int_0^\infty\bigl[e^{-(1+u)L(x)}-1\bigr]e^{-ux}\,dx}_{R(u)}.$$
*Derivation (by hand).* Substitute $g=e^{-x}$ in the paper's definition of $M(u;\alpha)$, set $w=(1-u)x$,
use $\sqrt{1-2\cos\alpha e^{-w}+e^{-2w}}=1+e^{-w}\widetilde B_\alpha(-w)$, and split the bracket. $\square$
**Branch points:** $\widetilde B_\alpha$ is singular iff $1-2\cos\alpha e^z+e^{2z}=0\iff e^z=e^{\pm i\alpha}\iff z=\pm i\alpha+2\pi i n$.
Verified to $\sim50$ digits at $(u,\alpha)\in\{0.1,0.3\}\times\{\pi/3,2\pi/3,3\pi/4,\pi/2\}$ (verify_dcya1.py prints residuals).

## 2. Weight expansion (L3) — PROVED

For $x>0$, $(1+e^{-x})^{-1-u}e^{-x}=\sum_{m\ge1}(-1)^{m-1}\binom{u+m-1}{m-1}e^{-mx}$, so
$T(u)=\sum_{m\ge1}(-1)^{m-1}\!\int_0^\infty\!\binom{u+m-1}{m-1}e^{-mx}\widetilde B_\alpha((u-1)x)dx$, and extracting $[u^k]$,
$M_k=\sum_{m\ge1}(-1)^{m-1}T_k^{(m)}+R_k$. The $m=1$ term carries the slowest decay; $m\ge2$ are "image saddles" of relative size $e^{-(\sqrt2-1)\theta}$ (K2/K2b). $\square$

## 3. Main term (L4) — PROVED (structure)

$T_k^{(1)}=[u^k]\!\int_0^\infty e^{-x}\widetilde B_\alpha((u-1)x)dx$. Cauchy extraction with $z=(u-1)x$ gives
$T_k^{(1)}=\int_0^\infty e^{-x}\frac{dx}{x}\,\frac1{2\pi i}\oint\widetilde B_\alpha(z)(1+z/x)^{-N}dz$ (contour around $z=-x$).
Deform to Hankel loops at the nearest branch points $z=\pm i\alpha$, where $\widetilde B_\alpha(z)=\gamma_\pm\sqrt{z\mp i\alpha}\,(1+\sum g_j^{(\pm)}(z\mp i\alpha)^j)$, $\gamma_+=\sqrt{2\sin\alpha}\,e^{i(\pi/4-\alpha/2)}$, $g_1=-\tfrac34-\tfrac{iE}{4\sin\alpha}$ ($E=e^{i\alpha}$).
The exact Hankel transfer $\frac1{2\pi i}\oint\zeta^{j+1/2}(1+\zeta/w)^{-N}d\zeta=\frac{i}{\pi}w^{j+3/2}\Gamma(j+\tfrac32)\frac{\Gamma(N-j-3/2)}{\Gamma(N)}$ reduces each local term to a **1-D integral**
$$I_+^{(j)}=\gamma_+ g_j\frac{i}{\pi}\Gamma(j+\tfrac32)\frac{\Gamma(N-j-\frac32)}{\Gamma(N)}\int_0^\infty e^{-x}x^{N-1}(x+i\alpha)^{-(N-j-3/2)}dx
=\gamma_+ g_j\tfrac{i}{\pi}\Gamma(j+\tfrac32)\Gamma(N-j-\tfrac32)(i\alpha)^{j+3/2}U(N,j+\tfrac52,i\alpha),$$
($g_0=1$; $U$ = Tricomi). Main term $M_k^{\rm main}=2\operatorname{Re}\sum_j I_+^{(j)}$.
**Numerical check:** leading$+g_1$ leaves $|M_k-M_k^{\rm main}|\approx0.45\,N^{-2}e^{-\theta}$ (the $\kappa_2$ term), i.e. $\kappa_1$ is captured exactly (verify_dcya1.py).

## 4. Saddle analysis of $I_+^{(0)}$ (L5)

Deform the ray $x=e^{i\pi/4}t$; phase $\Phi(t)=-bt+(N-1)\log(bt)-(N-\tfrac32)\log(bt+i\alpha)$, $b=e^{i\pi/4}$.
Saddle $t_0$ on the ray: $\operatorname{Re}\Phi'(t_0)=0$, linear phase $\omega=\operatorname{Im}\Phi'(t_0)$ ($\omega\sim c/\sqrt N$), curvature $c_2=-\operatorname{Re}\Phi''(t_0)\sim c/N$.
Model value $MV=\gamma_+\frac{i}{\pi}\Gamma(\tfrac32)\frac{\Gamma(N-3/2)}{\Gamma(N)}\,b\,e^{\Phi(t_0)}\sqrt{2\pi/(-\Phi''(t_0))}\,e^{\omega^2/(2\Phi''(t_0))}$.
**Validated numerically:** $|2\operatorname{Re}I_+^{(0)}-2\operatorname{Re}MV|/(N^{-3/2}e^{-\theta})\approx1.7$ ($N=400,1600$); $MV/{\rm lead}\to1$ ($1.09,1.05,1.03$ at $k=400,1600,6400$). Saddle series $F_1$ (first order, Gaussian moments of the phase remainder) and $F_2$ (second order) capture the remainder: $|I_0-MV-F_2|/{\rm sc}=0.83,0.31$ — decreasing, each term $O({\rm sc})$.

## 5. Numerical evidence (L2, L8)

Evaluator: exact formal $u$-series of the integrand at tanh-sinh nodes (shared across $k$), high dps ($\sim k$). Validated vs `corner_series.py` (2e-20 at $k\le20$) and K2's $M_{100}(\pi/2)=-1.7723168073590541000887\mathrm{e}{-10}$ (matched to 20 digits, dps=140).
Residual ratio $\rho=|M_k-{\rm lead}|/(N^{-3/2}e^{-\theta})$ (target: bounded):
- $\alpha=\pi/2$: $k=30{:}1.97$, $k=60{:}0.57$ (D13 known $\le1.65$).
- $\alpha=2\pi/3$: $M_{60}/{\rm lead}=1.12$; $\alpha=\pi/3$: $M_{60}/{\rm lead}=1.27$ (both $O(1)$).

## 6. The two gaps (OPEN, stated exactly)

- **L6 (main-term constant).** The saddle structure is fully derived and each piece ($MV$, $F_1$, $F_2$, tails) is numerically $O({\rm sc})$, but a hand-derived error function $E(N)$ with **arb certification** $E(N)\le C_1\,{\rm sc}$ for all $N\ge N_0$ (D13s-style subinterval/cell sweep) is **not completed**. Needed: subdivide the window so the phase-remainder $|R(u)|$ is small per subinterval, bound $|e^R-1|$ there, and sweep $N$ on cells with arb.
- **L7 (rest constant).** Bound $|M_k-M_k^{\rm main}|=|\sum_{m\ge2}(-1)^{m-1}T_k^{(m)}+\text{Watson tail}+R_k|\le C_2 e^{-\delta\sqrt N}$ ($\delta=\sqrt{2(2\pi-\alpha)}>\sqrt{2\alpha}$) by CX6-style domination, generalized to $\alpha$. The representation (L1,L3) and the branch-point locations are proved; the explicit domination constant is **not assembled**.

## 7. Verification

`verify_dcya1.py` (<5 min) prints: (i) the L1 residual $|M(u;\alpha)-(T(u)+R(u))|$ to ~50 digits; (ii) the L4 check that leading$+g_1$ captures $\kappa_1$; (iii) the decay-law residual ratio at small $k$. All numbers in this document that are labeled PROVED are derived by hand; all labeled NUMERICAL are printed by code.
