# SRF4 — The two uniformity steps of SRF3, written out

**Work directory:** `new_results/SRF4_octahedron_uniformity/`. Statuses: **PROVED** (every
step written here, by hand), **NUMERICAL** (printed by a script in this folder). Nothing outside this
folder was edited.

**Notation (as in SRF3).** $p>0$, $d=p-1$, $\tau=1-\frac1p$. $S(p)$ is the area of
$|x|^p+|y|^p+|z|^p=1$, $S(1)=4\sqrt3$, $S=8(1-\tau)^2T(\tau)$ (SRF3 Theorem A, proved there and refereed),
$$T(\tau)=\iint_\Delta f(\tau;t)\,dt_1dt_2,\qquad f(\tau;t)=(t_1t_2t_3)^{-\tau}\Bigl(\textstyle\sum_i t_i^{2\tau}\Bigr)^{1/2},$$
$\Delta=\{t_1,t_2\ge0,\ t_1+t_2\le1\}$, $t_3=1-t_1-t_2$; area of $\Delta$ is $\tfrac12$.
$z_0=\frac{-\log2+i\pi}2$, $R=|z_0|=\frac12\sqrt{\pi^2+\log^22}=1.6085752558559\ldots$,
$\theta=\arg z_0=1.7879528081\ldots$, $\rho^*=2R$. $a_k=[\tau^k]S/S(1)$, $T_k=[\tau^k]T$, so
$a_k=\frac2{\sqrt3}(T_k-2T_{k-1}+T_{k-2})$. $b_n=[d^n]S(1+d)/S(1)$.

---

## 0. The two steps, located precisely

**Step 1 (corner decomposition).** SRF3 `proof.md` lines 120–128:

> * **Corners** ($u<\delta$, $s<\delta'$ or $s>1-\delta'$; i.e. two small $t_i$): with $x=-\log s$, $y=-\log u$,
> $U=x+y$, $w=x/U$, the zero condition becomes $e^{-2zw}+e^{-2z(1-w)}=-1$. A continuation computation
> (step_06b) shows the root of smallest modulus has $|z|\ge R$ with equality only at $w\to0,1$ […] Hence on any
> compact $w$-subinterval of $(0,1)$ the minimum exceeds $R$, and away from the edge strips the corner pieces
> are analytic in a disc $|z|<R(1+\delta)$ […]

and the open item, lines 323–328:

> 1. **Corner decomposition (§3.1).** Needed: on compact $w$-sets $[w_0,1-w_0]$, $\min|z(w)|\ge R(1+\delta'')$
> […] $\|Q(\tau,s,u)-(\text{fiber})\|=o(1)$ on $|\tau|=R(1+\delta''/2)/U$ uniformly as $U\to\infty$ (Rouché);
> glue to the edge strips ($\delta'<\delta$) by a smooth partition of unity. Status: routine, unwritten […]

It is used in §3.2 (line 171: "The corner pieces get the same split in $U\ge U_0$: tail zeros are at
$\ge R(1+\delta)/U$") for the **upper bound** $\limsup|a_k/k!|^{1/k}\le1/R$.

**Step 2 (Hankel transfer, uniformity and phase).** SRF3 lines 63–70 (per-$(s,u)$ Hankel loop
$B_f(z)\sim-\frac{iE}{2\sqrt\pi}\frac{(\tau^*)^2}{z^{3/2}}e^{z/\tau^*}$), lines 232–267 (integration of that
asymptotic over $(s,x)$), in particular lines 254–267:

> The contour computation determines the singular strength up to the sign of the square-root determination
> […] the correct one […] is fixed independently and decisively by the cosine fit of §4.4

and the open item, lines 329–333:

> 2. **Hankel transfer phases (§4).** Needed: the per-$(s,u)$ Hankel asymptotic […] with error uniform in the
> tail, and analytic tracking of the $(\tau-\tau^*)^{1/2}$ determination (currently Theorem D's sign rests on
> the $10^{-8}$ cosine fit […]). Status: proof-sketch+

It is used in §3.4 (line 198: "By Theorem D (§4), $B_T(z)\sim\ldots$ with *nonzero* strength") for the **lower
bound** $\limsup\ge1/R$. So without Step 2 the claim of the paper (line 1925) was not proved in either
direction's sharp form: Step 1 for $\le$, Step 2 for $\ge$.

**What this note does.** Step 1 is replaced by an *exact* statement with a two-line proof (Lemma 1: the
"fiber" is not an approximation in the right normalization) plus an elementary zero-free lemma valid for
**all** $t$ (Lemma Z). This makes the interior/edge/corner partition and the Rouché gluing unnecessary.
Step 2 is proved as Lemma H (uniform Hankel transfer with an exact remainder formula) plus Lemma B (the
branch, determined analytically) plus Lemma L (the integration over the parameters). The result is a
complete proof of the claim at line 1925 and of the sharper growth law, in an additive-error form
(Theorem G). §8 lists what is still not proved.

---

## 1. Ordered log-ratio coordinates; the exact fiber — PROVED

By symmetry of $f$ under permutations of $(t_1,t_2,t_3)$, $T=6\iint_{t_1\ge t_2\ge t_3}f$. On that region put
$$\lambda=\log\frac{t_1}{t_2},\quad \ell=\log\frac{t_1}{t_3}\quad(0\le\lambda\le\ell),\qquad w=\frac\lambda\ell\in[0,1],\qquad
t_1=\frac1{1+e^{-\lambda}+e^{-\ell}}\in[\tfrac13,1].$$
So $\ell=\log(t_{\max}/t_{\min})$ is the log-spread of $t$.

**Lemma 1.** (a) $dt_2\,dt_3=t_1t_2t_3\,d\lambda\,d\ell=t_1^3e^{-\lambda-\ell}\,d\lambda\,d\ell$; hence
$$T(\tau)=6\int_0^\infty\!\!\int_0^\ell t_1^3e^{-\lambda-\ell}f\,d\lambda\,d\ell .$$
(b) With $P_w(z):=1+e^{-2wz}+e^{-2z}$ and $h_w(z):=e^{(1+w)z}P_w(z)^{1/2}=(e^{2(1+w)z}+e^{2z}+e^{2wz})^{1/2}$,
$$\textstyle\sum_it_i^{2\tau}=t_1^{2\tau}P_w(\tau\ell),\qquad f(\tau;t)=t_1^{-2\tau}h_w(\tau\ell)\quad\text{exactly.}$$

*Proof.* (a) $t_2=e^{-\lambda}/D$, $t_3=e^{-\ell}/D$, $D=1+e^{-\lambda}+e^{-\ell}$: $\partial_\lambda t_2=-t_2(1-t_2)$,
$\partial_\ell t_2=\partial_\lambda t_3=t_2t_3$, $\partial_\ell t_3=-t_3(1-t_3)$, determinant
$t_2t_3[(1-t_2)(1-t_3)-t_2t_3]=t_1t_2t_3$; and $dt_1dt_2=dt_2dt_3$ on $\Delta$. (b) $t_2=t_1e^{-w\ell}$,
$t_3=t_1e^{-\ell}$, so $\sum t_i^{2\tau}=t_1^{2\tau}(1+e^{-2\tau w\ell}+e^{-2\tau\ell})$ and
$(t_1t_2t_3)^{-\tau}=t_1^{-3\tau}e^{(1+w)\tau\ell}$. $\square$
[NUMERICAL, `constants.py`: $6\iint t_1^3e^{-\lambda-\ell}=0.5$ = area of $\Delta$.]

So the zero set of $\sum t_i^{2\tau}$ is *exactly* $\{\tau:\ \tau\ell\in Z(P_w)\}$: one family of entire functions of
one variable, indexed by $w\in[0,1]$. The edge ($t_3\to0$, $t_1\asymp t_2$) is $w\to0$; the corner (two small,
comparable coordinates) is $w\to1$. SRF3's corner "fiber" $e^{-2zw'}+e^{-2z(1-w')}=-1$ is $P$ in a different
normalization ($\tau=z_S/U$, $U=\lambda+\ell$), up to the approximations $t_1\approx s$, $t_2\approx1$; here no
approximation is made.

---

## 2. Lemma Z (uniform zero-free disc) — this is Step 1 — PROVED

**Lemma Z.** Let $w\in[0,1]$ and $\zeta=a+ib$ with $E_w(\zeta):=1+e^{w\zeta}+e^{\zeta}=0$. Then
(i) $|b|\ge\pi$; (ii) $|\zeta|\ge\rho^*=\sqrt{\pi^2+\log^22}$, with equality only if $w\in\{0,1\}$.
Equivalently ($\zeta=-2z$, $P_w(z)=E_w(-2z)$): every zero of $P_w$ has $|\operatorname{Im}z|\ge\pi/2$ and
$|z|\ge R$, with $|z|=R$ only for $w\in\{0,1\}$.

*Proof.* (i) The imaginary part of $E_w(\zeta)$ is $e^{wa}\sin(wb)+e^a\sin b$. If $0<b<\pi$ then $\sin b>0$ and
$0\le wb<\pi$ gives $\sin(wb)\ge0$, so it is $>0$. If $b=0$, $E_w(\zeta)\ge1$. For $-\pi<b<0$ conjugate. So $|b|\ge\pi$.

(ii) If $|b|\ge\rho^*$ then $|\zeta|\ge\rho^*$, with equality only if $a=0$, $|b|=\rho^*$; but $a=0$ makes
$1,e^{iwb},e^{ib}$ three unit vectors with zero sum, so $b\equiv\pm\frac{2\pi}3\pmod{2\pi}$, impossible for
$|b|=\rho^*=3.2171\ldots$ ($\rho^*\bmod 2\pi\notin\{\frac{2\pi}3,\frac{4\pi}3\}$). So assume
$\pi\le|b|<\rho^*$; by conjugation $b=\pi+\beta$, $0\le\beta<\rho^*-\pi=0.0756\ldots$

*Case $a\ge0$.* From $1+e^\zeta=-e^{w\zeta}$: $|1+e^\zeta|^2=e^{2wa}\ge1$. Since $\cos b=-\cos\beta$,
$|1+e^\zeta|^2=1-2e^a\cos\beta+e^{2a}$, so $e^{2a}\ge2e^a\cos\beta$, i.e. $a\ge\log(2\cos\beta)\ (\ge0)$. Hence
$|\zeta|^2\ge g(\beta):=\log^2(2\cos\beta)+(\pi+\beta)^2$, and
$g'(\beta)=-2\log(2\cos\beta)\tan\beta+2(\pi+\beta)\ge2\pi-2\log2\tan(0.0756)>0$, so $g(\beta)\ge g(0)=\rho^{*2}$,
with equality only at $\beta=0$. At $\beta=0$ ($b=\pi$) the imaginary part is $e^{wa}\sin(w\pi)$, which vanishes only
for $w\in\{0,1\}$.

*Case $a<0$.* $E_{1-w}(-\zeta)=e^{-\zeta}E_w(\zeta)=0$ and $\operatorname{Re}(-\zeta)>0$, $|\operatorname{Im}(-\zeta)|\in[\pi,\rho^*)$; apply the
first case to $(1-w,-\zeta)$. $\square$

[NUMERICAL, `lemma_Z_R.py`: all zeros in $|\zeta|<6.5$ for 1001 values of $w$ (Newton count = argument-principle
count for every $w$): $\min|\zeta|/\rho^*=1$ exactly at $w=0,1$, $=1.000500$ at $w=0.001,0.999$, $1.0259$ at $0.05$,
$1.302$ at $\frac12$, profile symmetric under $w\mapsto1-w$; $\min|\operatorname{Im}\zeta|/\pi=1$; the Step-3 inequality
$a\ge\log(2\cos\beta)$ holds at all 94 zeros in the band; $g$ increasing, $\min g'=2\pi$.]

**Corollary Z1 (zero-free region for every $t$).** For every $t\in\Delta$, $\sum_it_i^{2\tau}\ne0$ for
$|\tau|<R/\ell(t)$ and for $|\operatorname{Im}\tau|<\pi/(2\ell(t))$ (all $\tau$ if $\ell=0$). Hence $f(\cdot;t)$ is analytic there,
with the branch that is positive on the real axis.

*Proof.* Lemma 1(b), Lemma Z, and $t_1^{2\tau}\ne0$. The disc and the strip are simply connected and contain $0$,
so the square root has an analytic branch equal to $+\sqrt3\,t_1^{\tau}\cdots$ on $\mathbb R$. $\square$

**Corollary Z2 (SRF3's corner statement, with an explicit constant).** In SRF3's normalization
$\tau=z_S/U$, $U=\lambda+\ell=(1+w)\ell$, every zero satisfies $|z_S|\ge R(1+w)=R/\max(w',1-w')$, where
$w'=\lambda/U$ or $\ell/U$. On $[w_0',1-w_0']$ this gives $\min|z_S|\ge R/(1-w_0')$, i.e. SRF3's
$\delta''=w_0'/(1-w_0')$ (their numerical $1.081R$ at $w_0'=0.05$ is above the bound $1.053R$). And because Lemma 1 is
exact, the "$\|Q-\text{fiber}\|=o(1)$, Rouché, partition of unity" part of Step 1 is not needed at all. $\square$

---

## 3. Upper bound — PROVED

**Lemma D (Taylor coefficients pass under the integral).** $T\in C^\infty(-\frac1{16},\frac1{16})$ and
$T_k=\iint_\Delta[\tau^k]f(\cdot;t)\,dt$ for every $k$.

*Proof.* Let $\tau_0\in[-\frac1{16},\frac1{16}]$, $r=\min(\frac1{16},\frac{\pi}{4\ell})$. The disc $|\tau-\tau_0|\le r$ lies in
the strip of Corollary Z1, and on it, with $s=|\tau_0|+r\le\frac18$ and $t_{\min}\le t_i\le1$:
$|(t_1t_2t_3)^{-\tau}|\le t_{\min}^{-3s}$ and $|\sum t_i^{2\tau}|\le3t_{\min}^{-2s}$, so $|f|\le\sqrt3\,t_{\min}^{-1/2}$. Cauchy:
$|\partial_\tau^kf(\tau_0;t)|\le k!\,r^{-k}\sqrt3\,t_{\min}^{-1/2}\le k!(16+2\log\tfrac1{t_{\min}})^k\sqrt3\,t_{\min}^{-1/2}$
(using $\ell\le\log\frac1{t_{\min}}$, $4/\pi<2$). Under the uniform probability on $\Delta$,
$\Pr(t_{\min}<\eta)\le\sum_i\Pr(t_i<\eta)=3(1-(1-\eta)^2)\le6\eta$, so this majorant is integrable, uniformly in
$\tau_0$; differentiate under the integral $k$ times. $\square$

**Lemma U.** Fix $\varepsilon\in(0,R)$, $\ell_0>0$, $M:=9^{R/\ell_0}\sqrt3\,e^{2R}$. For every $t\in\Delta$ and $k\ge0$,
$$\bigl|[\tau^k]f(\cdot;t)\bigr|\le M\Bigl(\frac{\max(\ell(t),\ell_0)}{R-\varepsilon}\Bigr)^k .$$
*Proof.* Cauchy on $|\tau|=r=(R-\varepsilon)/\max(\ell,\ell_0)<R/\ell$ (Corollary Z1). There
$|t_1^{-2\tau}|\le3^{2r}\le9^{R/\ell_0}$, and $|z|=|\tau\ell|<R$ gives
$|h_w(z)|^2\le e^{2(1+w)|z|}+e^{2|z|}+e^{2w|z|}\le3e^{4R}$. $\square$
[NUMERICAL, `constants.py`: at 64 points $t$ (spreads up to $e^{-20}$), $k\le40$, $\varepsilon=0.05$, $\ell_0=1$:
largest ratio of the left side to the bound is $1.3\cdot10^{-2}$.]

**Proposition U.** $\limsup_k|T_k/k!|^{1/k}\le1/R$, and the same for $a_k$.
*Proof.* By Lemma D and Lemma U, $|T_k|\le\frac M2(R-\varepsilon)^{-k}\mathbb E[\max(\ell,\ell_0)^k]$, and
$\mathbb E[\ell^k]\le\mathbb E[\log^k\frac1{t_{\min}}]=\int_0^\infty kx^{k-1}\Pr(t_{\min}<e^{-x})dx\le6\,k!$. So
$|T_k|\le\frac M2(R-\varepsilon)^{-k}(6k!+\ell_0^k)$. Let $\varepsilon\downarrow0$. For $a_k$:
$|T_{k-j}|/k!\le|T_{k-j}|/(k-j)!$ for $j=1,2$. $\square$

This is the whole upper bound. No interior/edge/corner split, no $V$-split and no kernel convolution are
needed: SRF3's §3.2 failure at moderate $x$ (the referee's $0.72R/x$) happens because $\ell$, not $x=-\log t_3$, is
the correct scale; with $\ell$ the zero-free radius is exactly $R/\ell$, uniformly.

---

## 4. Lemma R (the two nearest zeros, uniformly for small $w$) — PROVED

Fix $\rho_1=2$, $w_0=0.19$.

**Lemma R.** For every $w\in[0,w_0]$, $P_w$ has exactly two zeros in $|z|<\rho_1$ and none on $|z|=\rho_1$. They are
simple and conjugate, $z^*(w)$ (upper half-plane) and $\overline{z^*(w)}$; $z^*$ is continuous on $[0,w_0]$,
analytic near $0$, with $z^*(0)=z_0$, $z^{*\prime}(0)=z_0/2$, $|z^*(w)|\ge R$ with equality only at $w=0$.

*Proof.* $P_0=2+e^{-2z}$ has zeros $-\frac{\log2}2-\frac{i\pi(2m+1)}2$, of moduli $R$ ($m=0,-1$) and
$\ge\frac12\sqrt{\log^22+9\pi^2}=4.725$, so exactly $z_0,\bar z_0$ in $|z|<2$.
*Lower bound for $|P_0|$ on $|z|=2$.* By symmetry take $z=2e^{i\varphi}$, $\varphi\in[0,\pi]$; $e^{-2z}=re^{i\psi}$,
$r=e^{-4\cos\varphi}$, $\psi=-4\sin\varphi$. If $r\le\frac12$: $|P_0|\ge2-r\ge\frac32$. If $r\ge\frac72$: $|P_0|\ge r-2\ge\frac32$.
Otherwise $\cos\varphi\in(-\frac{\log3.5}4,\frac{\log2}4)$, so $\varphi\in(1.3966,1.8895)$, $\sin\varphi>0.94969$ (concavity),
$|\psi|=4\sin\varphi\in(3.79876,4]\subset[\pi,2\pi]$, where $\cos$ increases, so $\cos\psi>\cos(3.79876)=-0.79173$. Hence
$|P_0|^2=4+r^2+4r\cos\psi\ge4+r^2-4r(0.79173)\ge4-(2\cdot0.79173)^2=1.4927$, so $|P_0|\ge1.2217$.
*Rouché.* $|P_w-P_0|=|e^{-2wz}-1|\le e^{2w|z|}-1\le e^{4w_0}-1=1.13828<1.2217\le|P_0|$ on $|z|=2$. So $P_w$ has two
zeros in $|z|<2$, none on the circle. $P_w(\bar z)=\overline{P_w(z)}$ and $P_w>0$ on $\mathbb R$, so they are a
conjugate non-real pair; a double zero would force four. $z^*(w)=\frac1{2\pi i}\oint_{\partial H}z\frac{P_w'}{P_w}dz$
over the boundary of the upper half-disc ($P_w\ne0$ there) is continuous. At $w=0$,
$\partial_zP_0(z_0)=-2e^{-2z_0}=4\ne0$, $\partial_wP(0,z_0)=-2z_0$, so by the implicit function theorem $z^*$ is analytic
near $0$ with $z^{*\prime}(0)=2z_0/4$. The last claim is Lemma Z. $\square$
[NUMERICAL, `lemma_H.py`, `lemma_Z_R.py`: hand bound $1.22175$; true minimum $1.49508$ (grid with Lipschitz
certificate $\ge1.49336$); zero count $=2$ at 501 values $w\in[0,0.05]$; $|z^*|/R$ increases from $1$ to $1.1089$ on
$[0,0.19]$; finite-difference $z^{*\prime}(0)=-0.173286+0.785399i=z_0/2$.]

---

## 5. Step 2: the Hankel transfer — PROVED

For $t$ with $w=\lambda/\ell\in[0,w_0]$ and $\ell\ge\ell_0=1$ (region III below), put
$$F_t(\zeta)=t_1^{-2\zeta/\ell}h_w(\zeta),\qquad\text{so }f(\tau;t)=F_t(\tau\ell),\quad[\tau^k]f=\ell^k[\zeta^k]F_t .$$
Let $\Omega_w=\{|\zeta|<\rho_1\}$ minus the two radial cuts $\Gamma=\{z^*(1+v):v\ge0\}$, $\bar\Gamma$, and
$v_1=\rho_1/|z^*|-1\in[0.121,0.244]$. By Lemma R, $P_w\ne0$ on $\overline{\Omega_w}\setminus\{z^*,\bar z^*\}$, and $\Omega_w$ is simply
connected, so $F_t$ (the continuation of its germ at $0$) is analytic on $\Omega_w$ and bounded:
$|F_t|\le M_F:=3^{2\rho_1/\ell_0}\sqrt3\,e^{(1+w_0)\rho_1}$ on $|\zeta|\le\rho_1$ (same estimate as Lemma U).

Define $D_w(\zeta)=\dfrac{P_w(\zeta)}{1-\zeta/z^*}=-z^*\!\int_0^1P_w'(z^*+s(\zeta-z^*))\,ds$. It is jointly continuous in
$(w,\zeta)$, with continuous $\zeta$-derivative, and $D_w\ne0$ on the segment $[0,z^*]$ and on $\Gamma\cap\{|\zeta|\le\rho_1\}$
($P_w$ has no other zero there, and $D_w(z^*)=-z^*P_w'(z^*)\ne0$ by simplicity). Let $D_w^{1/2}$ be the continuous
square root along $[0,z^*]\cup\Gamma$ with $D_w^{1/2}(0)=+\sqrt3$, and
$$\Psi^h_w(\zeta)=e^{(1+w)\zeta}D_w(\zeta)^{1/2},\qquad\Psi_t(\zeta)=t_1^{-2\zeta/\ell}\Psi^h_w(\zeta).$$

**Lemma B (the branch).** (a) Near $\Gamma$, $F_t(\zeta)=(1-\zeta/z^*)^{1/2}\Psi_t(\zeta)$ with the principal power, and
$\Psi_t$ is analytic across $\Gamma$. (b) $\Psi^h_0(z_0)=\sqrt2\,z_0^{1/2}$ (principal root). (c) $\Psi^h_w(z^*(w))\to\Psi^h_0(z_0)$ as
$w\to0$.

*Proof.* (a) On the segment $\zeta=sz^*$, $0\le s<1$, the principal $(1-s)^{1/2}>0$, so
$F_t/(1-\zeta/z^*)^{1/2}$ is the continuation from $\sqrt3$ of $t_1^{-2\zeta/\ell}e^{(1+w)\zeta}D_w^{1/2}$; both $F_t$ and
$(1-\zeta/z^*)^{1/2}$ change sign across $\Gamma$ (each is a square root of a function with a simple zero at $z^*$ and
no other zero nearby), so the quotient is continuous across $\Gamma$, hence analytic there (Morera).
(b) $D_0(sz_0)=\dfrac{2+2^se^{-i\pi s}}{1-s}$ for $s\in[0,1)$ (because $e^{-2sz_0}=2^se^{-i\pi s}$). Its imaginary part is
$-2^s\sin(\pi s)/(1-s)<0$ on $(0,1)$; it equals $3$ at $s=0$ and tends to $D_0(z_0)=-4z_0=2\log2-2\pi i$. So the
continuous argument stays in $(-\pi,0]$ and the continuous root is the principal one:
$D_0(z_0)^{1/2}=(-4z_0)^{1/2}=2(-z_0)^{1/2}=-2i\,z_0^{1/2}$ (as $\arg(-z_0)=\theta-\pi\in(-\frac\pi2,0)$). With
$e^{z_0}=e^{-\frac{\log2}2}e^{i\pi/2}=\frac i{\sqrt2}$: $\Psi^h_0(z_0)=\frac i{\sqrt2}\cdot(-2i)z_0^{1/2}=\sqrt2\,z_0^{1/2}$.
(c) $(s,w)\mapsto D_w(sz^*(w))$ is continuous and non-vanishing on the contractible set $[0,1]\times[0,w_0]$, so its
root lifted from $\sqrt3$ is jointly continuous. $\square$
[NUMERICAL, `lemma_H.py`: the tracked root at $z^*$ equals the principal root for $w=0,0.05,0.19$;
$\Psi^h_0(z_0)=1.12338847491682+1.39826637166739i=\sqrt2z_0^{1/2}$ to all 15 printed digits.]

This settles SRF3's open sign question analytically. SRF3's own contour computation (lines 255–258) produced
$-\pi\sqrt{3/8}$, i.e. the **wrong** determination (per edge $\frac{\pi i}4e^{z_0}z_0^{1/2}=-\frac{\pi\sqrt2}8z_0^{1/2}$); Lemma B gives
$+\frac{\pi\sqrt2}8z_0^{1/2}$ per edge orientation-pair, which is the sign their fit had selected.

**Uniform constants.** Let $K_{\rm cut}=\{(w,\zeta):w\in[0,w_0],\ \zeta=z^*(w)(1+v),\ 0\le v\le v_1(w)\}$ (compact).
By joint continuity, the following are finite:
$$M_h=\sup_{K_{\rm cut}}|\Psi^h_w|,\qquad L_h=\sup_{K_{\rm cut},\,v>0}\frac{|\Psi^h_w(z^*(1+v))-\Psi^h_w(z^*)|}{|z^*|v}$$
(the second because $\partial_\zeta\Psi^h_w=\Psi^h_w\,[(1+w)+D_w'/(2D_w)]$ is continuous on $K_{\rm cut}$ and $|D_w|$ has a
positive minimum there). Since $|t_1^{-2\zeta/\ell}|\le81$ and $|\partial_\zeta t_1^{-2\zeta/\ell}|\le162\log3$ for $|\zeta|\le2$,
$t_1\in[\frac13,1]$, $\ell\ge1$: $\sup|\Psi_t|\le M_\Psi:=81M_h$ and the secant quotient of $\Psi_t$ is at most
$L:=81L_h+162\log3\,M_h$, for every $t$ in region III.
[NUMERICAL, `constants.py`, 191 values of $w$ × 4000 points per cut: $M_h=1.7936$, $L_h=1.0794$,
$\min_{K_{\rm cut}}|D_w|=6.434$; $\sup_{|\zeta|=2}|h_w|=13.17$.]

**Lemma H (uniform Hankel transfer, exact remainder).** For $t$ in region III and $k\ge2$,
$$[\zeta^k]F_t=-\frac2\pi\,\mathrm B(\tfrac32,k-\tfrac12)\operatorname{Re}\bigl[z^{*-k}\Psi_t(z^*)\bigr]+\epsilon_k(t),\qquad
|\epsilon_k(t)|\le C_H\bigl(|z^*|^{-k}k^{-5/2}+\rho_1^{-k}\bigr),$$
with $C_H=\frac2\pi\bigl(L\rho_1\Gamma(\frac52)\sup_{k\ge2}\frac{k^{5/2}\Gamma(k-\frac32)}{\Gamma(k+1)}+1.25\cdot\frac\pi2M_\Psi\bigr)+M_F$, independent of $t$.

*Proof.* Exact decomposition: deform Cauchy's circle $|\zeta|=r<R$ to the boundary of $\Omega_w$ (outer arcs of
$|\zeta|=\rho_1$, both banks of each cut; $F_t$ is bounded, so the small circles around $z^*,\bar z^*$ contribute $0$ in the
limit). On $\Gamma$, $\zeta=z^*(1+v)$, the banks of $(1-\zeta/z^*)^{1/2}=(-v)^{1/2}$ are $\mp i\sqrt v$, and $\Psi_t$ is continuous
across (Lemma B(a)); traversing the boundary positively gives
$$[\zeta^k]F_t=-\frac1\pi\sum_{\zeta^*\in\{z^*,\bar z^*\}}\zeta^{*-k}J_k(\zeta^*)+C_k,\quad
J_k(\zeta^*)=\int_0^{v_1}\!\sqrt v\,\Psi(\zeta^*(1+v))(1+v)^{-k-1}dv,\quad|C_k|\le M_F\rho_1^{-k},$$
($\Psi$ for $\bar\Gamma$ is $\overline{\Psi_t(\bar\zeta)}$, by $F_t(\bar\zeta)=\overline{F_t(\zeta)}$, so the two cut terms are conjugate).
The sign is fixed by the case $\Psi\equiv1$, $\zeta^*=1$, $\rho_1\to\infty$:
$-\frac1\pi\mathrm B(\frac32,k-\frac12)=-\frac{\Gamma(k-\frac12)}{2\sqrt\pi\,k!}=[x^k](1-x)^{1/2}$.
Then $J_k=\Psi(\zeta^*)\mathrm B(\frac32,k-\frac12)+e_1+e_2$ with
$|e_1|\le L|\zeta^*|\int_0^\infty v^{3/2}(1+v)^{-k-1}dv=L|\zeta^*|\mathrm B(\frac52,k-\frac32)$ and
$|e_2|\le M_\Psi\int_{v_1}^\infty\sqrt v(1+v)^{-k-1}dv\le M_\Psi(1+v_1)^{1-k}\frac\pi2$; note
$|\zeta^*|^{-k}(1+v_1)^{-k}=\rho_1^{-k}$ and $(1+v_1)\le1.25$. $\square$
[NUMERICAL, `lemma_H.py`: the exact decomposition reproduces the Taylor coefficient (computed independently by
series arithmetic) to $10^{-13}$–$10^{-18}$ relative to $|z^*|^{-k}$, for $(w,\ell)=(0,40),(0.05,3),(0.1,10),(0.19,1)$ and
$k=5,10,20,40,60$; the leading-term error is $O(k^{-5/2}|z^*|^{-k})$ with constant below $1$ for $k\ge20$.]

---

## 6. The growth law and the limsup — PROVED

Split the ordered region (×6) into
I: $\ell<\ell_0$; II: $\ell\ge\ell_0$, $\lambda>w_0\ell$; III: $\ell\ge\ell_0$, $\lambda\le w_0\ell$; and $T_k=T_k^{\rm I}+T_k^{\rm II}+T_k^{\rm III}$ (Lemma D).

* **I.** Lemma U with $\ell<\ell_0$: $|T_k^{\rm I}|\le\frac M2(\ell_0/(R-\varepsilon))^k$.
* **II.** Lemma U and $\int_{w_0\ell}^\infty e^{-\lambda}d\lambda=e^{-w_0\ell}$:
  $|T_k^{\rm II}|\le6M(R-\varepsilon)^{-k}\int_0^\infty\ell^ke^{-(1+w_0)\ell}d\ell=\dfrac{6M\,k!}{(1+w_0)^{k+1}(R-\varepsilon)^k}$. Take
  $\varepsilon=0.1$: $(1+w_0)(R-\varepsilon)=1.795>R$.
  (The corners live here: their measure carries $e^{-(1+w)\ell}$, which beats the zeros at $|z|=R$ for $w$ near $1$.)
* **III.** Lemma H. The errors give $6C_H\iint e^{-\lambda-\ell}\ell^k(|z^*|^{-k}k^{-5/2}+\rho_1^{-k})\le6C_H\,k!\,(R^{-k}k^{-5/2}+\rho_1^{-k})$
  (using $|z^*|\ge R$). The main part is
  $$M_k=-\frac{12}\pi\mathrm B(\tfrac32,k-\tfrac12)\operatorname{Re}\bigl[z_0^{-k}I_k\bigr],\quad
  I_k=\iint_{\rm III}t_1^3e^{-\lambda-\ell}\ell^k\Bigl(\frac{z^*(\lambda/\ell)}{z_0}\Bigr)^{-k}\Psi_t\bigl(z^*(\lambda/\ell)\bigr)\,d\lambda\,d\ell .$$

**Lemma L.** $I_k/k!\to\Psi^h_0(z_0)\int_0^\infty e^{-3\lambda/2}(1+e^{-\lambda})^{-3}d\lambda=\sqrt2\,z_0^{1/2}\cdot\frac\pi{16}$.

*Proof.* $I_k/k!=\int_0^\infty e^{-\lambda}\,\mathbb E[\phi_k(X_k,\lambda)]\,d\lambda$ with $X_k\sim\Gamma(k+1,1)$ and
$\phi_k(\ell,\lambda)=\mathbf 1\{\ell\ge\max(\ell_0,\lambda/w_0)\}\,t_1^3\,(z^*(\lambda/\ell)/z_0)^{-k}\,\Psi_t(z^*(\lambda/\ell))$.
*Domination:* $|\phi_k|\le M_\Psi$ ($t_1\le1$; $|z^*|\ge|z_0|$ by Lemma Z). *Pointwise limit for fixed $\lambda$:* fix
$\delta\in(0,\frac12)$; $\Pr(|X_k-k|>\delta k)\le(k+1)/(\delta k)^2\to0$. On $|X_k-k|\le\delta k$, for $k$ large: the indicator is $1$;
$|t_1^3-(1+e^{-\lambda})^{-3}|\le3e^{-\ell}\to0$; $w=\lambda/\ell\to0$, so $\Psi^h_w(z^*(w))\to\Psi^h_0(z_0)$ (Lemma B(c)) and
$t_1^{-2z^*/\ell}\to1$; with $g(w)=\log(z^*(w)/z_0)$, $g(w)=\frac w2+o(w)$ (Lemma R), so
$k\,g(\lambda/\ell)\in[\frac\lambda{2(1+\delta)},\frac\lambda{2(1-\delta)}]+o(1)$. Hence
$\limsup_k|\mathbb E\phi_k-(1+e^{-\lambda})^{-3}e^{-\lambda/2}\Psi^h_0(z_0)|\le M_\Psi\sup_{|x-\lambda/2|\le\lambda\delta}|e^{-x}-e^{-\lambda/2}|$, and
$\delta\downarrow0$. Dominated convergence in $\lambda$ (majorant $M_\Psi e^{-\lambda}$) and Lemma B(b). The integral: $s=(1+e^{-\lambda})^{-1}$
gives $\int_{1/2}^1\sqrt{s(1-s)}\,ds=\frac12\mathrm B(\frac32,\frac32)=\frac\pi{16}$. $\square$
[NUMERICAL: $\int=0.19634954084936\ldots=\pi/16$ to 30 digits, `constants.py`.]

**Theorem G (growth law).** With $\mathrm B(\frac32,k-\frac12)=\frac{\sqrt\pi}2k^{-3/2}(1+O(\frac1k))$,
$$\frac{T_k}{k!}=-\frac{3\sqrt{2\pi}}8\,k^{-3/2}\operatorname{Re}\bigl[z_0^{1/2-k}\bigr]+o\bigl(k^{-3/2}R^{-k}\bigr),\qquad
\boxed{\ \frac{a_k}{k!}=-\sqrt{\frac{3\pi}8}\,k^{-3/2}\operatorname{Re}\bigl[z_0^{1/2-k}\bigr]+o\bigl(k^{-3/2}R^{-k}\bigr).}$$
*Proof.* Sum I–III: I and II are $O(k!\,c^{-k})$ with $c>R$ or $O(C^k)$; III is $M_k+O(k!k^{-5/2}R^{-k})$, and by Lemma L,
$M_k/k!=-\frac{12}\pi\cdot\frac{\sqrt\pi}2k^{-3/2}\operatorname{Re}[z_0^{-k}(\sqrt2z_0^{1/2}\frac\pi{16}+o(1))]$. Then
$a_k=\frac2{\sqrt3}(T_k-2T_{k-1}+T_{k-2})$ and $T_{k-j}/k!=O(k^{-j}k^{-3/2}R^{-k})$; $\frac2{\sqrt3}\cdot\frac{3\sqrt{2\pi}}8=\sqrt{3\pi/8}$. $\square$

The form is additive on purpose: $|\cos((k-\frac12)\theta)|$ comes arbitrarily close to $0$, and then
$a_k/(\text{law})$ need not be near $1$ (computed: the ratio is $0.022$ at $k=91$ and $5.44$ at $k=33$, where
$|\cos|=0.018$, $0.011$). So "$\sim$" in the paper and in SRF3's Theorem D must be read as this additive statement.

**Corollary C (the claim of line 1925).** $\limsup_k|a_k/k!|^{1/k}=1/R=2/\sqrt{\pi^2+\log^22}$.
*Proof.* $\le$: Proposition U. $\ge$: $|\operatorname{Re}z_0^{1/2-k}|=R^{1/2-k}|\cos((k-\frac12)\theta)|$. For any $x$,
$\max(|\cos x|,|\cos(x+\theta)|)\ge\sin\frac{\pi-\theta}2=0.6263$: with $\alpha,\beta$ the distances of $x$, $x+\theta$ to
$\frac\pi2+\pi\mathbb Z$, $\alpha+\beta\ge\operatorname{dist}(\theta,\pi\mathbb Z)=\pi-\theta$ and $|\cos|=\sin(\text{distance})$. So among any two
consecutive $k$ one has $|a_k/k!|\ge(0.6263\sqrt{3\pi/8}R^{1/2}+o(1))k^{-3/2}R^{-k}$. $\square$
(SRF3 used three consecutive $k$; two suffice. [NUMERICAL: minimum $0.626317$ on an $x$-grid, $0.626325$ over $k\le2\cdot10^5$.])

**Corollary Gd (the $d$-series).**
$$\frac{b_n}{n!}=-\sqrt{\frac{3\pi}4}\,n^{-3/2}\operatorname{Im}\bigl[z_0^{1/2-n}\bigr]+o(n^{-3/2}R^{-n}),\qquad\limsup|b_n/n!|^{1/n}=1/R,$$
and both series have radius $0$.
*Proof.* $b_n=\sum_{k=1}^na_k[d^n](\frac d{1+d})^k=\sum_k a_k(-1)^{n-k}\binom{n-1}{k-1}$, so with $\alpha_k=a_k/k!$,
$\frac{b_n}{n!}=\sum_{j=0}^{n-1}\frac{(-1)^j}{j!}(1-\frac jn)\alpha_{n-j}$. Insert Theorem G. For $j\le n/2$ use
$(n-j)^{-3/2}\le2^{3/2}n^{-3/2}$ and $\sum_jR^j/j!<\infty$ (dominated convergence in $j$, termwise limit
$n^{-3/2}z_0^{1/2-n}(-z_0)^j/j!$); the terms $j>n/2$ total $O(R^{-n}\sum_{j>n/2}R^j/j!)$, super-exponentially small. The main sum gives
$\operatorname{Re}[z_0^{1/2-n}e^{-z_0}]$ and $e^{-z_0}=-i\sqrt2$, so $-\sqrt{3\pi/8}\operatorname{Re}[-i\sqrt2\,X]=-\sqrt{3\pi/4}\operatorname{Im}X$. Limsup as in
Corollary C with $|\sin((n-\frac12)\theta)|$ (same two-consecutive bound). Radius $0$ in $d$ also follows directly: a
convergent $d$-series would make the $\tau$-series convergent by composition with $d=\tau/(1-\tau)$. $\square$
(This confirms the referee's correction of SRF3's first $d$-display: the sign is $-\sqrt{3\pi/4}\operatorname{Im}$, as the paper has it.)

---

## 7. The Borel transform — PROVED, with one gap stated

Let $B(z)=\sum_ka_kz^k/k!$ and $B_T(z)=\sum_kT_kz^k/k!$.

**Proposition S.** (a) $B$ is analytic in $|z|<R$ and continuous on $|z|\le R$; its radius is exactly $R$.
(b) Every $z_1$ with $|z_1|=R$, $z_1\notin\{z_0,\bar z_0\}$, has a neighbourhood to which $B$ continues analytically.
(c) $z_0,\bar z_0$ are singular: $|B'(sz_0)|=c(1-s)^{-1/2}(1+o(1))$ as $s\uparrow1$, with
$c=\frac12\sqrt{3\pi/8}\,\sqrt{\pi/R}=0.7584\ldots\ne0$ (numerically $|B'(sz_0)|\sqrt{1-s}=0.470,0.479,0.507$ at $s=0.8,0.85,0.9$, rising slowly towards $c$, as the $O(1)$ terms predict). Moreover,
$B(z)-\pi\sqrt{3/8}\,[z_0^{1/2}(1-z/z_0)^{1/2}+\bar z_0^{1/2}(1-z/\bar z_0)^{1/2}]$ has Taylor coefficients $o(k^{-3/2}R^{-k})$.

*Proof.* (a) Theorem G: $\sum|a_k/k!|R^k<\infty$; radius $R$ by Corollary C.
(b) Work with $B_T$ (then $B=\frac2{\sqrt3}(B_T-2\!\int_0^z\!B_T+\!\int_0^z\!\!\int_0^uB_T)$, and primitives exist on the simply
connected union of the disc and a small disc at $z_1$). $B_T=B_{\rm I}+B_{\rm II}+B_{\rm III}$. $B_{\rm I}$ is entire and $B_{\rm II}$ is
analytic in $|z|<1.795$ (bounds above). For $|z|<R$, summing Lemma H's exact decomposition against $z^k\ell^k/k!$ (Tonelli:
$\sum_k(\ell|z|/R)^k/k!=e^{\ell|z|/R}$ against $e^{-\ell}$) gives
$$B_{\rm III}(z)=6\!\iint_{\rm III}\!t_1^3e^{-\lambda-\ell}\Bigl[\frac1{2\pi i}\!\oint_{|\zeta|=\rho_1}\!\!F_t(\zeta)e^{\ell z/\zeta}\frac{d\zeta}\zeta
-\frac1\pi\sum_{\zeta^*}\int_0^{v_1}\!\!\sqrt v\,\Psi(\zeta^*(1+v))\,e^{\ell z/(\zeta^*(1+v))}\frac{dv}{1+v}\Bigr]d\lambda\,d\ell .$$
The circle term converges locally uniformly for $|z|<\rho_1=2$. For the cut terms: on the compact set
$K=\{z^*(w)(1+v)\}$ (and its conjugate) we have $\operatorname{Re}(z_1/\zeta)\le|z_1|/|\zeta|\le1$, and $=1$ forces $\zeta=z_1$, $|\zeta|=R$, hence
$v=0$, $|z^*(w)|=R$, $w=0$ (Lemma Z), $\zeta=z_0\ne z_1$. So $\max_K\operatorname{Re}(z_1/\zeta)=1-2\delta<1$, and for $|z-z_1|<\delta R$,
$\operatorname{Re}(z/\zeta)\le1-\delta$: the integrand is $O(e^{-\lambda-\delta\ell})$ uniformly, and the integral is analytic there (Morera/Fubini).
It agrees with $B_{\rm III}$ on the overlap with $|z|<R$.
(c) $B'(z)=\sum k\alpha_kz^{k-1}$. With $\alpha_k=-\kappa k^{-3/2}\operatorname{Re}z_0^{1/2-k}+o(k^{-3/2}R^{-k})$, $\kappa=\sqrt{3\pi/8}$:
$B'(sz_0)=-\frac\kappa2z_0^{-1/2}\sum k^{-1/2}s^{k-1}-\frac\kappa2\bar z_0^{1/2}z_0^{-1}\sum k^{-1/2}s^{k-1}e^{2ik\theta}+o((1-s)^{-1/2})$.
The first sum is $\sqrt\pi(1-s)^{-1/2}(1+o(1))$; the second is bounded uniformly in $s$ (Dirichlet test, $e^{2i\theta}\ne1$). If $B$
were analytic at $z_0$, $B'$ would be bounded near $z_0$. The coefficient statement is Theorem G with
$[z^k]z_0^{1/2}(1-z/z_0)^{1/2}=-\frac{z_0^{1/2-k}}{2\sqrt\pi}k^{-3/2}(1+O(\frac1k))$ and $\pi\sqrt{3/8}/\sqrt\pi=\sqrt{3\pi/8}$. $\square$

**What is not proved (the gap).** That $z_0$ is a *square-root branch point in the analytic sense*, i.e.
$B=A+(1-z/z_0)^{1/2}C$ with $A,C$ analytic in a full neighbourhood of $z_0$. The cut representation above diverges
on the set $\{z^*(w)(1+v)(1+w)\}$, a thin region emanating from $z_0$ outward, so it does not continue $B$ there,
and nothing here excludes further singularities of the continuation accumulating at $z_0$ from outside the
disc. What is proved is the square-root *type* in the sense of (c). SRF3's Theorem B ("only singularities on
$|z|=R$ are $z_0,\bar z_0$, both square-root branch points") is now proved in its first half and in this weak
form of its second half. (The 2D analogue, GRW1 T2b, has a similar acknowledged gap.)

---

## 8. Hypotheses SRF3 did not state; errors found

1. **The right scale is $\ell=\log(t_{\max}/t_{\min})$, not $x=-\log t_3$.** SRF3's edge strip $s\in[\delta,1-\delta]$,
   $u=e^{-x}$ does not control $t_{\min}$ for moderate $x$; that is exactly why its §3.2 contour failed (referee: zeros
   at $0.72R/x$). With $\ell$ the zero-free radius is exactly $R/\ell$ for every $t$ (Corollary Z1). The V-split repair is
   then unnecessary.
2. **Step 1's Rouché/partition step is unnecessary** (Lemma 1 is exact). Its needed inequality is Corollary Z2,
   with the explicit $\delta''=w_0'/(1-w_0')$.
3. **Step 2 needs, and SRF3 did not state:** (i) uniqueness and simplicity of the nearest zero *uniformly* in the
   edge parameter (Lemma R, with explicit $\rho_1=2$, $w_0=0.19$); (ii) a uniform Lipschitz bound for the regular
   factor along the whole cut (constants $M_h,L_h$, by compactness); (iii) a continuity argument for the branch
   (Lemma B). SRF3's raw contour sign was wrong (it gives $-\pi\sqrt{3/8}$; the fit-selected $+\pi\sqrt{3/8}$ is the correct
   one, now derived). Also SRF3's "the $s$-integral contributes $B(1+\frac z{2z_0},1+\frac z{2z_0})$" is a heuristic
   interchange; here the profile integral $\pi/16$ appears in Lemma L with dominated convergence.
4. **"$\sim$" is wrong as a ratio statement.** The growth law holds with additive error $o(k^{-3/2}R^{-k})$ (Theorem G);
   the ratio to the law does not tend to $1$ (it is $0.022$ at $k=91$, $5.44$ at $k=33$). This affects SRF3 Theorem D
   and the paper's sentence at line 1925.
5. **The paper says the sharper law is numerical.** With Theorem G and Corollary Gd it is proved (additive form).
6. **"Square-root branch points"** (paper line 1925, SRF3 Theorem B) is proved only in the sense of Proposition S(c).
7. Not an error but an improvement: the limsup lower bound needs two consecutive $k$, not three.

---

## 9. Numerical evidence (all printed by scripts in this folder)

**Coefficients.** `coeffs.py` computes $a_k$ exactly in terms of $\Gamma$-Taylor coefficients, through the identity
$\Gamma(3-2\tau)T(\tau)=\mathbb E[(g_1g_2g_3)^{-\tau}(\sum g_i^{2\tau})^{1/2}]$, $g_i$ i.i.d. Exp(1) (the Dirichlet–Gamma
factorization; $f$ is homogeneous of degree $-2\tau$), expanding $(3+\sum X_i)^{1/2}$, $X_i=g_i^{2\tau}-1$, and using
$\mathbb E[g^{-\tau}X^m]=\sum_k[x^k](e^{-x}(e^{2x}-1)^m)\,\Gamma^{(k)}(1)\tau^k$ (no cancellation in the inner series).
Computed $K=150$ at 330 digits (100 s); $K=80$ at 160 and 220 digits.

| check | result |
|---|---|
| $T_0=\sqrt3/2$; $[d^0]\ldots[d^4]$ vs thm:oct | exact to $10^{-330}$ |
| $\log\frac{(2-p)S}{4\sqrt3}$, $[\tau^2]\ldots[\tau^6]$ vs thm:oct | $\le1.5\cdot10^{-300}$ |
| vs SRF3 `SRF3_hi_coeffs_K120.json` (different engine), $k\le120$ | 149.4 digits ($\tau$), 149.5 ($d$) — the file's precision |
| $K150$/dps330 vs $K80$/dps220 | 191 digits |

**Growth law** (`growth.py`; scale $=k^{-3/2}R^{-k}$; residual $=(a_k/k!-\text{law})/\text{scale}$, should $\to0$):

| $k$ | $\lvert a_k/k!\rvert^{1/k}$ | residual ($\tau$) | $k\cdot$residual | $\lvert b_k/k!\rvert^{1/k}$ | residual ($d$) |
|---|---|---|---|---|---|
| 10 | 0.41802 | 0.1997 | 2.00 | 0.46829 | 0.0034 |
| 20 | 0.50387 | 0.0271 | 0.54 | 0.48045 | 0.0726 |
| 30 | 0.52510 | −0.0477 | −1.43 | 0.52835 | 0.0514 |
| 60 | 0.56316 | 0.0154 | 0.93 | 0.55955 | −0.0265 |
| 100 | 0.57637 | −0.0199 | −1.99 | 0.58360 | 0.0098 |
| 150 | 0.59242 | 0.0039 | 0.58 | 0.58855 | 0.0091 |

$1/R=0.6216681478591508$; $\max_{130\le k\le150}|a_k/k!|^{1/k}R=0.953$ (the $k^{-3/2}$ factor), and with it removed
$\max_{130\le k\le150}|a_k/(k!k^{-3/2})|^{1/k}R=1.0024$. $\max_{100\le k\le150}k|\text{residual}|=2.16$: numerically the $o(\cdot)$ of
Theorem G is $O(k^{-5/2}R^{-k})$. Least-squares fit $a_k/(k!\,\text{scale})=\text{amp}\cos(k\theta+\varphi)+$ two $1/k$-orders, $k=70..150$:
amp $=1.3766016$ (theory $\sqrt{3\pi R/8}=1.3766111$), $\varphi=2.2476202$ (theory $\pi-\theta/2=2.2476162$); $d$-series amp
$1.9468294$ (theory $\sqrt{3\pi R/4}=1.9468221$), phase $-2.4647742$ (theory $-\theta/2-\pi/2=-2.4647727$).
(SRF3/SRF3R got 8–10 digits with more correction terms; not repeated here.)

**Lemmas on grids.** Lemma Z: 1001 values of $w$, all zeros in $|\zeta|<6.5$, counts verified by the argument principle;
min ratio exactly 1 only at the endpoints. Lemma R: hand bound $1.2217$ vs true $1.4951$; 501-$w$ count check.
Lemma B: tracked branch = principal at $w=0,0.05,0.19$; $\Psi^h_0(z_0)=\sqrt2z_0^{1/2}$. Lemma H: exact decomposition to
$10^{-13}$–$10^{-18}$; constants $M_h=1.794$, $L_h=1.079$, $\min|D_w|=6.43$. Lemma U: ratio $\le0.013$ on 64 points × $k\le40$.
Lemma L integral $=\pi/16$, measure $=1/2$. Two-consecutive lemma: min $0.626317=\sin\frac{\pi-\theta}2$.

Scripts: `coeffs.py K dps`, `crosscheck.py`, `growth.py coeffs_K150_dps330.json`, `lemma_Z_R.py`, `lemma_H.py`
(≈4 min), `constants.py`. JSON files hold the coefficients. The virtualenv used to run them has been deleted;
`python3 -m venv .venv && .venv/bin/pip install mpmath numpy` recreates it.

---

## 10. Status table

| # | step | status |
|---|---|---|
| A | $S=8(1-\tau)^2T$ and equivalent forms (SRF3 §1) | PROVED (SRF3, refereed) |
| 1 | ordered log-ratio coordinates; exact fiber $\sum t_i^{2\tau}=t_1^{2\tau}P_w(\tau\ell)$ | PROVED (here) |
| Z | zeros of $P_w$: $|z|\ge R$ (equality iff $w\in\{0,1\}$), $|\operatorname{Im}z|\ge\pi/2$, all $w\in[0,1]$ — **Step 1** | PROVED (here, by hand) |
| Z2 | SRF3's corner bound with $\delta''=w_0'/(1-w_0')$; Rouché/gluing unnecessary | PROVED (here) |
| D | Taylor coefficients under the integral | PROVED (here) |
| U | uniform Cauchy bound; $\limsup\le1/R$ | PROVED (here) |
| SRF3 §3.2 V-split, §3.1 interior/kernel pieces, §3.3 non-cancellation | superseded (not needed) | — |
| R | exactly two simple zeros in $|z|<2$ for $w\le0.19$; $z^{*\prime}(0)=z_0/2$ | PROVED (here, by hand) |
| B | branch of $\sqrt{D_w}$; $\Psi^h_0(z_0)=\sqrt2z_0^{1/2}$ (fixes SRF3's sign) — **Step 2, phase** | PROVED (here) |
| H | uniform Hankel transfer, exact remainder, constants by compactness — **Step 2, uniformity** | PROVED (here); constant values NUMERICAL |
| L | parameter integration, profile $\pi/16$ | PROVED (here) |
| G | $a_k/k!=-\sqrt{3\pi/8}k^{-3/2}\operatorname{Re}z_0^{1/2-k}+o(k^{-3/2}R^{-k})$ | PROVED (here); rate $O(k^{-5/2}R^{-k})$ NUMERICAL |
| C | $\limsup|a_k/k!|^{1/k}=1/R$; radius 0 in $\tau$ | PROVED (here) |
| Gd | $d$-law with $-\sqrt{3\pi/4}\operatorname{Im}$; radius 0 in $d$ | PROVED (here) |
| S(a,b) | Borel transform analytic in $|z|<R$ and across every other point of $|z|=R$ | PROVED (here) |
| S(c) | singular at $z_0,\bar z_0$; square-root type (coefficients, $B'\sim c(1-s)^{-1/2}$) | PROVED (here) |
| S′ | local analytic form $A+(1-z/z_0)^{1/2}C$ near $z_0$ ("branch point" in the strict sense) | NOT PROVED |
| E | Borel summability along $\mathbb R_+$ (SRF3 §5) | NUMERICAL (untouched) |
| — | fitted amplitude/phase to 5–7 digits; residual table | NUMERICAL |

None of the "PROVED (here)" items has yet been checked by an independent referee.

**Verdict.** Yes: with this note the claim at line 1925 is fully proved — $\limsup_k|a_k/k!|^{1/k}=1/|z_0|$, radius zero
in $\tau$ and in $d$, and a Borel transform that is analytic in $|z|<|z_0|$ and singular at $z_0$ and $\bar z_0$. Every step is
written, and no uniformity step is left as a sketch. More is proved than the paper claims (the sharper law, in additive form). One
phrase over-claims: "square-root branch points" holds only in the sense of Proposition S(c).

**Sentence the paper can honestly state** (replacing the last part of the paragraph at line 1925):

> Then $\limsup_k|a_k/k!|^{1/k}=1/|z_0|=2/\sqrt{\pi^2+\log^22}=0.621668\ldots$, so the series in $\tau$, and the one in $d$,
> has radius zero. More precisely,
> $a_k/k!=-\sqrt{3\pi/8}\,k^{-3/2}\operatorname{Re}z_0^{1/2-k}+o(k^{-3/2}|z_0|^{-k})$, and the coefficients in $d$ satisfy the same with
> $-\sqrt{3\pi/4}\,\operatorname{Im}$ in place of $-\sqrt{3\pi/8}\,\operatorname{Re}$. The Borel transform in $\tau$ is analytic in $|z|<|z_0|$ and across every
> other point of the circle $|z|=|z_0|$, and is singular at $z_0$ and $\bar z_0$, where its derivative grows like
> $|1-z/z_0|^{-1/2}$. These points come from the edges of the octahedron, through $1+2e^{2z}=0$ (proved, in our notes; not
> yet independently refereed).

and the status-table row (line 2120) can read "Proved (in our notes, every step written; referee pending)" instead of
"Proved†, refereed", the referee having passed SRF3 only with the two steps at sketch level.
