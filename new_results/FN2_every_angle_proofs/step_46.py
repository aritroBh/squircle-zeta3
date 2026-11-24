import time
part1 = r"""# FN2: Ahmed's Problem 3 at every angle — proofs of $\mathcal J_2,\mathcal J_3$ for general $f$, and the contact law

**Status: both results PROVED.** The $\theta=\pi/2$ proof of
`research/07_push/P3_fn_graph/claude/result.md` (called **[P3]** below) transports to every
angle, because the only $\theta$-dependent inputs are (i) the layer-cake density, (ii) the
moments of the log-slope, and (iii) the tail series of the density at $\infty$; all three are
settled by Lemmas A, B, C. The algebra that turns the master expansion into the paper's
displayed $\mathcal J_2,\mathcal J_3$ is exact (sympy, difference $0$), as is its reduction
to the proved $\theta=\pi/2$ formulas and, for $f(x)=x$, to the exact series.

## 0. Setting

- $\theta\in(0,\pi)$ fixed, $c=\cos\theta$, $z=\sin^2\frac\theta2=\frac{1-c}2$.
- $u,v$ unit vectors, $u\cdot v=c$; $\gamma_n(x)=x\,u+f(x)^n\,v$, $L_n=\operatorname{length}\gamma_n$.
- $f\in C^1[0,1]$ nondecreasing, $f(0)=0$, $f(1)=1$, $a:=f'(1)>0$ (the non-contact case);
  $f\in C^{K+2}$ on $[1-\eta_0,1]$ for the expansion to order $K$.
- $y=y_n=f^n$; $m_y(s)=|\{x:y'(x)>s\}|$.
- $h_c(s)=1+s-\sqrt{1+2cs+s^2}$, so that $2-L_n=\int_0^1 h_c(y'(x))\,dx$ because
  $|u+sv|^2=1+2cs+s^2$ and $\int_0^1(1+y')\,dx=2$.
- $\phi_\theta=h_c'$, the slope density; $S_\theta$ the random slope with density
  $\phi_\theta/(1-c)$; $\xi_\theta=\log(1/S_\theta)$ with cumulants $k_j(\theta)$.
- $\varepsilon=\mu=\frac1{n-1}$, $a=f'(1)$, $A=f''(1)$, $B=f'''(1)$, $D=f^{(4)}(1)$,
  $\Xi=\log(na)+k_1(\theta)$, $\mu_1=\Xi$, $\mu_2=\Xi^2+k_2$, $\mu_3=\Xi^3+3k_2\Xi+k_3$.

## Lemma A (layer cake at angle $\theta$) — PROVED

For every nondecreasing absolutely continuous $y:[0,1]\to[0,1]$ with $y(1)-y(0)=1$,
$$2-\operatorname{length}(x\,u+y(x)\,v)=\int_0^1 h_c(y')\,dx=\int_0^\infty \phi_\theta(s)\,m_y(s)\,ds=(1-c)\,\mathbb E[m_y(S_\theta)].$$

*Proof.* $\int_0^1(1+y')\,dx=2$, so $2-L=\int_0^1(1+y'-|u+y'v|)\,dx=\int_0^1h_c(y')\,dx$.
$h_c(0)=0$ and $h_c(s)=\int_0^s\phi_\theta(t)\,dt$ with $\phi_\theta\ge0$ (Lemma B(i)), so
$h_c(y'(x))=\int_0^\infty\phi_\theta(s)\mathbf 1\{y'(x)>s\}\,ds$ and Tonelli gives the second
identity. $\int_0^\infty\phi_\theta=h_c(\infty)-h_c(0)=1-c$ (Lemma B(i)), so
$\phi_\theta/(1-c)$ is a probability density. $\square$

*Remark.* For $y=f^n$ the left endpoint plays no role: $m_y(s)=0$ for $s\ge y'(1)=na$, and
$m_y$ is computed by Lemma D. The identity is exact; no asymptotics yet.

## Lemma B (the density $\phi_\theta$) — PROVED

(i) $\displaystyle\phi_\theta(s)=1-\frac{s+c}{\sqrt{1+2cs+s^2}}
=\frac{1-c^2}{\sqrt{1+2cs+s^2}\,\bigl(\sqrt{1+2cs+s^2}+s+c\bigr)}>0$,
it is strictly decreasing on $[0,\infty)$ (since
$\phi_\theta'(s)=-(1-c^2)(1+2cs+s^2)^{-3/2}<0$), $\phi_\theta(0)=1-c$,
$\phi_\theta(s)\to0$, hence $0<\phi_\theta\le 1-c$ and
$\int_0^\infty\phi_\theta\,ds=h_c(\infty)=1-c$.

(ii) Bound: for $s\ge1$, $\sqrt{1+2cs+s^2}=\sqrt{(s+c)^2+1-c^2}\ge s+c\ge(1-|c|)s$, so
$$0<\phi_\theta(s)\le\frac{(1-c^2)}{2(1-|c|)^2\,s^2}=\frac{1+|c|}{2(1-|c|)}\,s^{-2}.$$

(iii) Tail series: for $s>1$, with $t=1/s$,
$$\phi_\theta(s)=1-(1+ct)(1+2ct+t^2)^{-1/2}=\sum_{l\ge2}d_l(c)\,s^{-l},
\qquad d_l(c)=\frac{(-1)^l\,(1-c^2)\,P_{l-1}'(c)}{l},$$
$P_l$ the Legendre polynomial. The series converges for $s>1$ (the roots of $1+2ct+t^2$ are
$e^{\pm i(\pi-\theta)}$, modulus $1$), so for every $L$,
$\phi_\theta(s)=\sum_{l=2}^Ld_l s^{-l}+O_{\theta,L}(s^{-L-1})$ uniformly for $s\ge2$.
In particular $d_2=(1-c^2)/2=2z(1-z)$, $d_3=-c(1-c^2)=-4cz(1-z)$,
$d_4=\frac38(1-c^2)(5c^2-1)$; at $c=0$ only even $l$ survive and
$d_{2m}(0)=-\binom{-1/2}m$, the coefficients $b_l$ of [P3].

*Proof.* (i) rationalize; (ii) as displayed; (iii) expand
$(1+2ct+t^2)^{-1/2}=\sum_{k\ge0}(-1)^kP_k(c)t^k$ (generating function at $-c$, using
$P_k(-c)=(-1)^kP_k(c)$), multiply by $1+ct$; the coefficient of $t^k$ is
$(-1)^k\bigl(P_k(c)-cP_{k-1}(c)\bigr)=(-1)^{k+1}(1-c^2)P_{k-1}'(c)/k$ by the Legendre
identity $(k+1)\bigl(P_{k+1}(x)-xP_k(x)\bigr)=(x^2-1)P_k'(x)$ shifted by one; the $k=0,1$
terms cancel ($P_1(c)-cP_0(c)=0$). The $c=0$ values: $P_{2m-1}'(0)/(2m)=-\binom{-1/2}m$.
All identities verified symbolically in `step_24` (sympy, exact). $\square$

## Lemma C (moments and cumulants of the log-slope) — PROVED

Let $\xi_\theta=\log(1/S_\theta)$.

(i) $\mathbb P(S_\theta\le s)\le s$ for $s\ge0$ (since $\phi_\theta\le1-c$), and
$\mathbb P(S_\theta>s)=O_\theta(s^{-1})$ by Lemma B(ii). Hence
$\mathbb E[S_\theta^w]<\infty$ exactly for $-1<w<1$; all moments of $|\xi_\theta|$ are
finite, and $\mathbb E[e^{t|\xi_\theta|}]<\infty$ for $|t|<1$.

(ii) The Mellin transform is the paper's (proved) Legendre integral [paper, Thm. P2, first
identity, integrated by parts]:
$$\mathbb E[S_\theta^w]=\frac{\pi w}{\sin\pi w}\,\frac{\Lambda_w(\theta)}{1-\cos\theta},
\qquad -1<w<1,\quad \Lambda_w(\theta)=\int_{\cos\theta}^1P_w(t)\,dt=2z\,{}_2F_1(-w,1+w;2;z).$$
Both sides are analytic in $-1<\operatorname{Re}w<1$ (the expectation by (i)), so
$\log\mathbb E[S_\theta^w]=\sum_{j\ge1}a_jw^j$ near $0$, with $a_j$ the paper's corner
coefficients (Thm. P2asym); the cumulants are
$k_j(\theta)=(-1)^j j!\,a_j$; in particular, with $w_2=\cot^2\frac\theta2$ and $\mathrm S_{1,2}$
the Nielsen polylogarithm,
$$k_1=1+w_2\log\cos^2\tfrac\theta2,\quad
k_2=1+\tfrac{\pi^2}3-2\operatorname{Li}_2(z)-w_2^2\log^2\cos^2\tfrac\theta2,\quad
k_3=2+2w_2^3\log^3\cos^2\tfrac\theta2-12w_2\,\mathrm S_{1,2}(z),$$
reducing at $\theta=\pi/2$ to $1-\log2$, $1+\pi^2/6$, $2-\frac32\zeta(3)$. The central
moments are $m_2=k_2$, $m_3=k_3$, $m_4=k_4+3k_2^2$.

*Proof.* (i) from Lemma B(i),(ii); (ii) is quoted from the paper (proved there); the
cumulant identification is [P3]'s argument applied to the analytic CGF: differentiate
$\log\mathbb E[e^{w\log S_\theta}]$ at $w=0$. Numerical confirmation to $10^{-60}$ in
`step_28` (7 angles, stable form of $\phi_\theta$). $\square$

## Lemma D (the window; $\theta$-independent) — PROVED in [P3], restated

Assume $f\in C^2[1-\eta,1]$ with $f'\ge a/2$ there, $\theta_f:=f(1-\eta)<1$,
$\sigma_n:=n\theta_f^{\,n-1}\|f'\|_\infty$. For $n\ge n_0(f)$:
$y_n'$ is strictly increasing on $[1-\eta,1]$; for $\sigma_n\le s<na$ the set
$\{y_n'>s\}=(1-\tau_n(s),1]$ with
$$(n-1)\mathcal A(\tau)-\mathcal B(\tau)=\log\tfrac ns,\qquad
\mathcal A(\tau):=-\log f(1-\tau),\quad \mathcal B(\tau):=\log f'(1-\tau);$$
$m_y(s)=0$ for $s\ge na$; and $\bigl|D_n-\int_{\sigma_n}^{na}\phi_\theta(s)\tau_n(s)\,ds\bigr|
\le(1-c)\sigma_n$, exponentially small. (The last bound uses Lemma B(i):
$\int_0^{\sigma_n}\phi_\theta\le(1-c)\sigma_n$; in [P3] the same line with $g'\le1$.)

## Lemma E (the bivariate inverse; $\theta$-independent) — PROVED in [P3], restated

$\Theta(W,\mu)$ solves $\mathcal A(\Theta)-\mu\tilde{\mathcal B}(\Theta)=W$,
$\tilde{\mathcal B}(\tau)=\log(f'(1-\tau)/a)$, $\Theta(0,\mu)=0$; it is $C^{K+1}$ near
$(0,0)$ (analytic if $f$ is), with Taylor coefficients $\theta_{ij}$ depending only on
$f'(1),\dots,f^{(i+j)}(1)$; $\theta_{1j}=a^{-1}(-A/a^2)^j$; $\tau_n(s)=\Theta(W_s,\mu)$ with
$W_s=\mu\log\frac{na}s$, and uniformly for $s\in[\sigma_n,na]$,
$$\tau_n(s)=\sum_{i+j\le K}\theta_{ij}\,\mu^{i+j}\log^i\tfrac{na}s
+O\Bigl(\mu^{K+1}\bigl(1+\bigl|\log\tfrac{na}s\bigr|\bigr)^{K+1}\Bigr).$$
"""
with open("proof.md","w") as f:
    f.write(part1)
print("part1 written", len(part1))
