# DCY2 — the first correction of the decay law of $M_k$: proof

**Theorem (DCY2, main term).** Put $N=k+1$, $\theta=\sqrt{\pi N}$, and
$$\varpi_1=\frac{(24+12\pi+\pi^2)+i\,(24-12\pi+\pi^2)}{48\sqrt\pi}
=0.84121508776758439707534863791-0.04501183768517361657373510376\,i .$$
Then for every integer $N\ge 100$,
$$\Bigl|\,M_k^{\rm main}+\frac{\sqrt\pi}{N}\,e^{-\theta}\,
\operatorname{Re}\bigl[e^{i\theta}\bigl(1+\frac{\varpi_1}{\sqrt N}\bigr)\bigr]\Bigr|
\;\le\; C_2\,N^{-2}e^{-\theta},\qquad C_2 = 1{,}000{,}000 .$$

**Corollary (full $M_k$).** With CX6's proved rest bound
$|M_k-M_k^{\rm main}|\le 10{,}600{,}000\,e^{-2\sqrt N}$ ($N\ge100$):
$$\Bigl|\,M_k+\frac{\sqrt\pi}{N}\,e^{-\theta}\,
\operatorname{Re}\bigl[e^{i\theta}\bigl(1+\frac{\varpi_1}{\sqrt N}\bigr)\bigr]\Bigr|
\;\le\; C_3\,N^{-2}e^{-\theta},\qquad C_3 = C_2+1.85\times10^{10},$$
for all $N\ge100$; for $N\ge N_3=25{,}711$ the rest term is $\le N^{-2}e^{-\theta}$, so
$C_3$ can be taken $=C_2+1$ there; for $100\le N\le 400$ the *sharp* residuals are certified
against the table `mk_deg9.tsv` (max res/bound $=0.00141$, i.e. constant $\approx 1409$).

## Lemma ledger

| # | Lemma | Status |
|---|---|---|
| 1 | Ray representation $M_k^{\rm main}=\frac2\pi\operatorname{Im}\tilde I$, $\tilde I=\int_{\rm ray}H(s)e^{i\pi s/2}(s/(1+s))^N ds/s$, $H(s)=\frac12B(\frac s2,\frac32)$ | **PROVED** (CX6 §3) |
| 2 | $\varpi_1$ closed form from the saddle expansion (two formalisms; exact algebra, diff $=0$; numeric match to the paper's 81 digits) | **PROVED** (`step_02_varpi1.py`) |
| 3 | Exponent factorization $J=\frac2\pi e^{\phi(s_0)}\omega\int_{-\rho}^{\infty}e^{\ell(t)}F(t)dt$; exact modulus $M(t)$; single-max lemma ($P(u)$ cubic, $P(0^+)>0$, $P(1)<0$, $P''<0$); tail rates $r_\Delta=\min(-\Lambda'(\Delta),b)$, $\lambda_L$ | **PROVED** (by hand; certified pointwise) |
| 4 | $|\ell_m|\le Ny^{m+1}+y^m/m$; $c_2=-\operatorname{Re}\ell_2\ge0.9b_Ry$; $\operatorname{Re}\ell_1=-y(1+b_Ry/(1+\sqrt2y+y^2))$ exact, $|{\rm Re}\,\ell_1|\le y(1+b_Ry)$ | **PROVED** (hand; $c_2$ certified as a polynomial inequality on $[0,y_1]$) |
| 5 | Watson lemma for $H$: $g(s_0)=1-\frac{3}{4s_0}+O(y^2)$ ($C_g\le1$), $h_j$ expansions, $M_4'$; $C_\rho\le0.5$ | **PROVED** (Watson with remainder) + **CERTIFIED** (`step_03_watson.py`) |
| 6 | Kept-polynomial moments $\mu_p^\infty$ in closed form; model $=T_1+O(N^{-2}e^{-\theta})$; the identities $\sum(\text{moment terms})_{O(y)}=0$ and $S_1+g'_1+B_2/(2b)+\gamma_1=\varpi_1\sqrt{2/\pi}$ | **PROVED** (sympy exact, diff $=0$) + **CERTIFIED** part-(b) constants |
| 7 | Remainder $|E(t)|$ by $|\cdot|$-moments with $[T,\Delta]$-split $e^{\gamma|t|^3}$ factors and damped incomplete-$\Gamma$ edges (implementation matches mpmath to 12+ digits) | **CERTIFIED** (`dcy2_cert.py`); audit vs truth all OK (`dcy2_audit.py`) |
| 8 | Window-tail bounds $\mathrm{TP}_p$ (exact incomplete $\Gamma$); outer tails via $M(\pm\Delta)$ with the single-dip lemma ($W(u)$ convex, $W(1)<0$, $W''>0$ for $u\ge1$) | **CERTIFIED** + **PROVED** (monotonicity structure) |
| 9 | Region B sweep $[400,10^5]$: $\mathrm{LHS}_{\rm cell}\le C_2$ on every one of 279 cells | **CERTIFIED** (`dcy2_sweep.py`, worst = 390,055) |
| 10 | Region C $[10^5,\infty)$: $\mathrm{LHS}(N)\le\mathrm{LHS}(10^5)=74{,}115\le C_2$ | **PROVED** (y-scaling, §4) + **CERTIFIED** (pointwise) |
| 11 | Finite range $[100,400]$: certified quadrature, $\max$ res/bound $=0.00141$ | **CERTIFIED** (`step_06_finite.py`, 301 quadratures) |
| 12 | Part (b): $|{\rm model}-T_1|\le C_bN^{-2}e^{-\theta}$, $C_b=29.5$; constants $C_F=10.15$ | **CERTIFIED** (`step_07_partb.py` + §5 assembly) |
| 13 | Rest combination: $N_3=25{,}711$, $C_3=C_2+10{,}600{,}000\cdot h(309)$, $h(309)=1748.98$ | **CERTIFIED** (arb) |
| 14 | Validation at $N=100,200,400$ vs `mk_deg9.tsv`: residual scales like $N^{-2}e^{-\theta}$ | **CERTIFIED** (`verify_dcy2.py` §3) |

## 1. Setup and the exact representation (L1)

$H(s)=\frac12B(\frac s2,\frac32)=\frac{\sqrt\pi}{4}\frac{\Gamma(s/2)}{\Gamma(3/2+s/2)}
=\int_0^\infty\sqrt{1-e^{-2u}}\,e^{-su}\,du$ (analytic in $\operatorname{Re}s>0$),
$\omega=e^{i\pi/4}$, ray $s=\omega x$, $x\ge0$. Put
$$J:=\frac2\pi\int_0^\infty H(\omega x)e^{i\pi\omega x/2}\Bigl(\frac{\omega x}{1+\omega x}\Bigr)^N\frac{\omega\,dx}{\omega x},
\qquad M_k^{\rm main}=\operatorname{Im}J.$$
We prove the complex estimate $|J-T_1|\le C_2N^{-2}e^{-\theta}$ with
$T_1=-i\frac{\sqrt\pi}{N}e^{-\theta+i\theta}(1+\varpi_1/\sqrt N)$; taking imaginary parts
gives the theorem since $|\operatorname{Im}w|\le|w|$ and
$\operatorname{Im}(-i\sqrt\pi N^{-1}e^{-\theta+i\theta}(1+\varpi_1/\sqrt N))
=-\frac{\sqrt\pi}{N}e^{-\theta}\operatorname{Re}[e^{i\theta}(1+\varpi_1/\sqrt N)]$.
(All code objects named \texttt{It}/\texttt{I\_tilde}/\texttt{model}/\texttt{target} are $J$-quantities.)

Notation: $y=\sqrt{\pi/(2N)}$, $\rho=1/y$, $s_0=\omega/y=\rho\omega$ (ray point nearest the saddle),
$\theta=\sqrt{\pi N}=\pi/(\sqrt2 y)$, $b=\pi\bar\omega/2$, $b_R=\operatorname{Re}b=\pi/(2\sqrt2)$,
$\kappa^2=\kappa_0^2+4\log(N/400)$ with $\kappa_0=2.48$, $\sigma=(2\cdot0.9b_Ry)^{-1/2}$,
$\Delta=\kappa\sigma$ (window half-width), $v=\sqrt y$.

## 2. $\varpi_1$ derived, not imported (L2)

With $s=s_0+\omega t$, $t=\tau/v$, the exponent
$G(y,\tau)=\phi(s_0+\omega t)-\phi(s_0)$ expands as
$G=-b\tau^2+vG_1\tau+v^2G_2\tau^2+\cdots$ with
$$G_1=-\bigl(1+\tfrac{\pi\bar\omega^2}{2}\bigr)\tau+\tfrac{\pi\bar\omega}{2}\tau^3,\qquad
G_2=\bigl(\tfrac12+\tfrac{3\pi\bar\omega^2}{4}\bigr)\tau^2-\tfrac{\pi\bar\omega}{2}\tau^4,$$
and the $H$-ratio $(1+v\tau)^{-3/2}=1+h_1v\tau+h_2v^2\tau^2+\cdots$ ($h_1=-3/2$, $h_2=15/8$).
Then
$$c_1=\Bigl\langle G_2+\frac{G_1^2}{2}+h_2\tau^2+G_1h_1\tau\Bigr\rangle_{e^{-b\tau^2}}
=\frac{(1+i)(15-\pi^2)}{8\sqrt2\,\pi},\qquad
\varpi_1=\sqrt{\frac\pi2}\Bigl(c_1+S_1-\tfrac34\bar\omega\Bigr),$$
$S_1=\pi(1+i)/(6\sqrt2)$ the $O(y)$ term of $e^{\phi(s_0)}=-iye^{-\theta+i\theta}e^{S(y)}$.
**`step_02_varpi1.py` prints $\varpi_1-\frac{(24+12\pi+\pi^2)+i(24-12\pi+\pi^2)}{48\sqrt\pi}=0$
exactly** (sympy), and the 40-digit values match the paper.

Equivalent $t$-formalism (exact $\ell_1,\ell_2$): the moment sum $\sum p_p\mu_p^\infty$ has
$O(y)$ coefficient **exactly $0$** (sympy), and the prefactor identity
$S_1+g'_1+B_2/(2b)+\gamma_1=\varpi_1\sqrt{2/\pi}$ holds **exactly** (sympy, diff $=0$), with
$B_2=\tfrac12+\tfrac34\pi\bar\omega^2$, $\gamma_1=A_1^2/(4b)$, $A_1=1+\pi\bar\omega^2/2$.

## 3. The analytic tail (L3–L9)

**Factorization (L3).** For $t>-\rho$: $\phi(s_0+\omega t)-\phi(s_0)=\ell(t)$ with
$$\ell(t)=(N-1)\log(1+yt)-N\log\Bigl(1+\frac{yt}{1+y\bar\omega}\Bigr)+\frac{i\pi\omega t}{2},$$
so $\tilde I=\frac2\pi e^{\phi(s_0)}\omega\int_{-\rho}^{\infty}e^{\ell(t)}F(t)\,dt$,
$F(t)=H(s_0+\omega t)/H(s_0)$. Modulus $M(t)=|e^{\ell(t)}|$:
$2\log M(t)=2\Lambda(t)=2(N-1)\log(1+yt)+N\log(1+\sqrt2y+y^2)-N\log Q(t)-2b_Rt$ with
$Q(t)=(1+yt)^2+\sqrt2y(1+yt)+y^2$. With $u=1+yt$,
$\operatorname{sgn}\Lambda'(t)=\operatorname{sgn}P(u)$,
$$P(u)=-2b_Ru^3-(\pi+2)yu^2+\bigl(\sqrt2(N-2)y^2-2b_Ry^3\bigr)u+2(N-1)y^3,$$
$P(0^+)>0$, $P(1)<0$, $P''=-12b_Ru-2(\pi+2)y<0$, hence exactly one root $u_c\in(0,1)$:
$M$ increases on $(-\rho,t_c)$, decreases on $(t_c,\infty)$; $t_c\in(-\tfrac12,0)$.
Right tail: $M(t)\le M(\Delta)e^{-r_\Delta(t-\Delta)}$,
$r_\Delta=\min(-\Lambda'(\Delta),b_R)$; the single-dip structure
($W(u)=Nu^2(u^2+\sqrt2yu)-(N-1)Q(u)^2$ convex, $W(1)<0$ certified by hand for all $N\ge100$)
gives $\max_{t\ge\Delta}\Lambda'(t)=\max(\Lambda'(\Delta),-b_R)$. Left tail:
$M(t)\le M(-\Delta)e^{-\lambda_L(-\Delta-t)}$,
$\lambda_L=(\pi y/2)/(2u_DQ(u_D))$, $u_D=1-y\Delta$ (a valid lower rate since
$\pi y/2\le P(u)$ on $(0,u_D]$).

**Window (L7).** On $|t|\le\Delta$: $\ell(t)=\ell_1t+\ell_2t^2+R_3(t)$ with exact
$$\ell_1=(N-1)y-\frac{Ny}{1+y\bar\omega}+\frac{i\pi\omega}{2},\qquad
\ell_2=-\frac{(N-1)y^2}{2}+\frac{Ny^2}{2(1+y\bar\omega)^2}.$$
$F(t)e^{R_3(t)}=P(t)+E(t)$ with the 14-monomial kept polynomial $P$
($w\le2$ plus $w=4$-odd, exact shifted-Gaussian moments $\mu_p^\infty$), and $E(t)$ bounded by
$|\cdot|$-moments with $[T,\Delta]$-split factors $e^{\gamma|t|^3}$, $T=\min(\Delta,c_2/(2\gamma))$,
and damped incomplete-$\Gamma$ edges.

**$\ell$-boxes (L4).** $|\ell_m|\le Ny^{m+1}+y^m/m$ for $m\ge2$;
$c_2=-\operatorname{Re}\ell_2\ge0.9b_Ry$ (polynomial inequality, certified on $[0,y_1]$);
$\operatorname{Re}\ell_1=-y\bigl(1+\frac{b_Ry}{1+\sqrt2y+y^2}\bigr)$ (exact), so
$|{\rm Re}\,\ell_1|\le y(1+b_Ry)\le c_2\Delta/2$.

**$H$ (L5).** Watson with remainder on $H(s)=\sqrt2\int u^{1/2}v(u)e^{-su}du$,
$v(u)=\sqrt{(1-e^{-2u})/(2u)}=1-\frac u2+\frac{5u^2}{24}-\frac{u^3}{16}+\rho_3(u)$,
$|\rho_3|\le C_\rho u^3$ on $[0,1]$ ($C_\rho\le0.5$ certified) and $\le 2.42u^3$ on $[1,\infty)$:
$g(s_0)=1-\frac{3}{4s_0}+r_g$, $|r_g|\le C_gy^2$ ($C_g\le1$);
$h_j=\omega^jH^{(j)}(s_0)/(j!H(s_0))=(-1)^j(\frac32)_jy^j/j!\cdot(1+O(y))$;
$|F(t)-\sum_0^3h_jt^j|\le\frac{M_4'}{24}|t|^4$,
$M_4'=\sqrt2\,\Gamma(\tfrac{11}{2})\bigl(\tfrac{1/y-\Delta}{\sqrt2}\bigr)^{-11/2}/H_{\rm lo}$.

**Assembly (L7–L9).**
$J-{\rm model}=(2/π)e^{\phi(s_0)}\omega\bigl[\int_{|t|\le\Delta}e^{\ell_1t+\ell_2t^2}E(t)dt
-\sum p_p(\mu_p^\infty-\mu_p(\Delta))\bigr]+\text{tails}$,
bounded as in `dcy2_budget.py`/`dcy2_cert.py`; all pieces audited against the true values
(`dcy2_audit.py`: ratios $0.001$–$0.7$).

## 4. Region C: $N\ge10^5$ (L10)

Every piece of $\mathrm{LHS}(N)=E_{\rm total}N^2e^\theta$ is non-increasing for $N\ge10^5$.
Write each piece as $y^{a}\cdot(\text{factors})$: after folding
$e^{-\kappa^2/4}=\frac{800}{\pi}e^{-\kappa_0^2/4}y^2$ into the tail pieces, every piece has
$a\ge0$ with all correction factors non-increasing (checked piecewise in §3 of result.md):
$M_4'$-piece $y^0$; $\alpha^4M_{12}$ $y^0$; $B_6$-piece $y^0$; $\beta_2M_8$ $y^0$;
$\beta_3M_{10}$ $y^0$; $E_{\rm drop}$ pieces $\ge y^0$; $E_{\rm tailmom}$ pieces $\ge y^0$
with $Q_p(\kappa^2/4)$ growing at most logarithmically and $Q'_p/Q_p<1$ (since $Q'_p=Q_p-x^{\deg Q_p}$);
tails: $e^{-\kappa^2/4}\kappa^a$ decreasing for $\kappa^2=28.2\ge2a$; $\kappa\sqrt y$ decreasing
for $N\ge635$; $\kappa^3\sqrt y$ decreasing for $N\ge1221$; $\lambda_L$ decreasing.
The log-derivative lemma: $d/dN[y^aQ_p(\kappa^2/4)]<0$ iff $a\kappa^2/4>\deg(Q_p)$, which holds
($\kappa^2=28.2\ge16$). Pointwise evidence (`dcy2_regionc.py`):
$\mathrm{LHS}=74115,39589,2436,1040,833,803$ at $N=10^5,10^6,10^9,10^{12},10^{16},10^{20}$.

## 5. Part (b): model $-T_1$ (L6, L12)

${\rm model}=(2/π)e^{\phi(s_0)}H(s_0)\omega\sum p_p\mu_p^\infty
=-iye^{-\theta+i\theta}\cdot{\rm ModelRatio}$,
${\rm ModelRatio}=(2/π)e^{S(y)}H(s_0)\omega\sum p_p\mu_p^\infty$;
$T_1=-iye^{-\theta+i\theta}(2y/\sqrt\pi)(1+\varpi_1/\sqrt N)$.
Then $|{\rm model}-T_1|\cdot N^2e^{\theta}=(\pi^{3/2}/2)(1+\varpi_1/\sqrt N)\cdot|\delta_F|$
with $\delta_F$ the $O(y^2)$ remainder of $e^{S}g\,\mathrm{SF}\,\mathrm{EF}\,\mathrm{MSF}$,
$\mathrm{SF}=(1-(B_2/b)y)^{-1/2}$, $\mathrm{EF}=e^{-\ell_2t_0^2}$, $\mathrm{MSF}=1+\sum p_pM_p$:
$|\delta_F|\le C_F=C_S+C_g+C_{S2}+C_E+C_{MS}+\mathrm{cross}$ with
$C_S=e^{S_{\rm sup}}[\pi^2/(72(1-y)^2)+\pi/(8(1-y))]\le0.596$, $C_g\le1$,
$C_{S2}=\frac38|B_2/b|^2/(1-|B_2/b|y_1)^{5/2}\le2.543$, $C_E\le3.164$, $C_{MS}\le0.3555$,
$\mathrm{cross}=|S_1g'_1|+|S_1B_2/(2b)|+|S_1\gamma_1|+|g'_1B_2/(2b)|+|g'_1\gamma_1|+|B_2\gamma_1/(2b)|
\le2.496$; total $C_F\le10.155$, whence $C_b=(\pi^{3/2}/2)(1+|\varpi_1|/\sqrt{400})\cdot10.1=29.3$.
Pointwise: $D(N)=|{\rm model}-T_1|\cdot N^2e^\theta=1.657,1.620,1.601,1.594,1.578,1.565$
at $N=400,1000,2000,3000,10^4,10^5$.

## 6. Finite range and validation (L11, L14)

Certified quadrature of $M_k^{\rm main}$ at every integer $N\in[100,400]$
(segments $[0,1,4,16,64,200]$ + proved tail $\le\frac{2\sqrt2}{\pi aR^2}e^{-aR}$, CX6 eq. (20)),
301 values: max res/bound $=0.00141$ at $N=381$, i.e. the constant needed there is $\approx1409$.
At $N=100,200,400$: residual/($C_2N^{-2}e^{-\theta}$) $=0.000165,0.000751,0.000932$ —
the residual scales like $N^{-2}e^{-\theta}$, as required; $|M_k-M_k^{\rm main}|$ at the same
points $=1.65\times10^{-13},4.2\times10^{-18},4.7\times10^{-25}$ (certified table balls).

## 7. Rest combination (L13)

$10{,}600{,}000\,e^{-2\sqrt N}\le N^{-2}e^{-\sqrt{\pi N}}$ iff
$2\log N+\log(10{,}600{,}000)\le(2-\sqrt\pi)\sqrt N$; certified $N_3=25{,}711$
(the left side decreases for $\sqrt N\ge4/(2-\sqrt\pi)$, i.e. $N\ge309$).
For $100\le N\le N_3$: $C_3=C_2+10{,}600{,}000\cdot\max_{[100,N_3]}N^2e^{-(2-\sqrt\pi)\sqrt N}
=C_2+10{,}600{,}000\cdot h(309)$ with $h(309)=1748.98$, so $C_3=1.85\times10^{10}$.

## 8. Files and reproduction

`verify_dcy2.py` (< 5 min, currently 87 s) reproduces everything:
finite range §3, region B §4, region C §5, part (b) §6, rest §7.
Modules: `dcy2_cell.py`, `dcy2_sweep.py`, `dcy2_regionc.py`,
`step_02_varpi1.py`, `step_03_watson.py`, `step_06_finite.py`, `step_07_partb.py`.
Prototype/audit: `dcy2_common.py`, `dcy2_budget.py`, `dcy2_audit.py`, `dcy2_cert.py`.
All supplied inputs remain read-only.
