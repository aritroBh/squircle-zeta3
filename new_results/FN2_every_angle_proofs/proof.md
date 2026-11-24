# FN2: Ahmed's Problem 3 at every angle — proofs of $\mathcal J_2,\mathcal J_3$ for general $f$, and the contact law

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
All identities verified symbolically in `step_25` (sympy, exact). $\square$

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
`step_29` (7 angles, stable form of $\phi_\theta$). $\square$

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

## Theorem 1 (the expansion at every angle) — PROVED

Let $f$ be as in §0, $f\in C^{K+2}$ near $1$. Then
$$\boxed{\;2-L_n=\sum_{k=1}^K\varepsilon^k\mathcal J_k(\Xi)+O\bigl(\varepsilon^{K+1}\log^{K+1}n\bigr)\;}$$
where the implied constant depends only on $\theta$ and the germ of $f$ at $1$, and
$$\mathcal J_k=\sum_{i=1}^k\theta_{i,k-i}\,M_i(\Xi)+T_k,\qquad
M_i(\Xi)=\int_0^\infty\phi_\theta(s)\log^i\tfrac{na}s\,ds
=(1-c)\,\mathbb E\bigl[(\Xi+\xi_\theta-k_1)^i\bigr]=2z\sum_{r=0}^i\binom ir m_r(\theta)\Xi^{i-r},$$
in particular $M_1=2z\,\mu_1$, $M_2=2z\,\mu_2$, $M_3=2z\,\mu_3$; and $T_k$ is the slope-tail
term, the coefficient of $\mu^k$ in
$$\sum_{i,j}\theta_{ij}\mu^{i+j}(-1)^{i+1}i!\sum_{l\ge2}
\frac{d_l(c)}{(l-1)^{i+1}}\Bigl(\frac{\mu}{a(1+\mu)}\Bigr)^{l-1}.$$
$\mathcal J_k$ is a polynomial of degree $k$ in $\Xi$ depending only on
$f'(1),\dots,f^{(k)}(1)$ and $\theta$; its top coefficient is
$2z\,[W^k]\mathcal A^{-1}(W)$.

*Proof.* Insert Lemma E's expansion of $\tau_n(s)$ into Lemma D's
$D_n=\int_{\sigma_n}^{na}\phi_\theta\,\tau_n\,ds+O((1-c)\sigma_n)$ and integrate term by
term. This is exactly [P3]'s proof of its Theorem 2 with $g'$ replaced by $\phi_\theta$.
The **$\theta$-audit** — every property of $g'$ used there, and its replacement:

| # | step of [P3]'s proof | property of $g'$ used | replacement |
|---|---|---|---|
| 1 | layer cake $D_n=\int g'm_y$ | $g'$ a density | Lemma A (exact, every $\theta$) |
| 2 | $[0,\sigma_n]$ cutoff | $0\le g'\le1$ | $0<\phi_\theta\le1-c$ (Lemma B(i)) |
| 3 | $\int_0^\infty g'\log^i\frac{na}s\,ds=M_i$ | all moments of $\log S$ finite; values | Lemma C(i),(ii): $\mathbb E|\log S_\theta|^r<\infty$; $M_i=2z\sum\binom ir m_r(\theta)\Xi^{i-r}$ |
| 4 | $s\ge na$: tail series $g'(s)=\sum b_ls^{-2l}$ | expansion with remainder | Lemma B(iii): $\phi_\theta=\sum_{l\ge2}d_ls^{-l}$, convergent for $s>1$ |
| 5 | the window $\tau_n=\Theta(\mu\log\frac{na}s,\mu)$ | none (no $\theta$ in $\Theta$) | Lemma E, unchanged |
| 6 | uniformity of Taylor's remainder along $[0,(W_s,\mu)]$ | none | Lemma E, unchanged |
| 7 | cumulants $k_j$ | Ahmed's $\log R$ | $k_j(\theta)=(-1)^jj!a_j$ from the paper's proved Legendre identity (Lemma C(ii)) |

Steps 1–4 and 7 are the only $\theta$-sensitive ones; each is settled above. In detail:

- The Taylor remainder integrates to
  $O\bigl(\mu^{K+1}(1-c)\,\mathbb E(1+|\log\tfrac{na}{S_\theta}|)^{K+1}\bigr)
  =O(\mu^{K+1}\log^{K+1}n)$, because $\log\frac{na}{S_\theta}=\Xi+(\xi_\theta-k_1)$ and all
  moments of $\xi_\theta$ are finite (Lemma C(i)).
- The extension $\int_{\sigma_n}^{na}\to\int_0^{na}$ adds $\le(1-c)\sigma_n(C+|\log\sigma_n|)^K$,
  exponentially small.
- The main term of step 3 is exactly $M_i(\Xi)$ as displayed:
  $\int_0^\infty\phi_\theta\log^i\frac{na}s\,ds=(1-c)\mathbb E[(\log na+\xi_\theta)^i]$,
  and $\log na=\Xi-k_1$.
- For $s\ge na$ ($na\ge2$ for large $n$), Lemma B(iii) gives, with
  $\int_{na}^\infty s^{-l}\log^i\frac{na}s\,ds=(na)^{1-l}(-1)^ii!/(l-1)^{i+1}$ and
  $(na)^{1-l}=\bigl(\frac\mu{a(1+\mu)}\bigr)^{l-1}$, exactly the $T$-part (the sign
  $(-1)^{i+1}$: minus for subtracting the $s>na$ range, times $(-1)^i$ from
  $\log^i\frac{na}s=(-1)^i\log^i\frac s{na}$). Summing the tail series first (it converges
  uniformly on $s\ge na\ge2$) and truncating the $\mu$-series at order $K$ leaves an error
  $O(\mu^{K+1})$, absorbed into the remainder.

Collecting powers of $\mu=\varepsilon$ gives the theorem. The dependence and degree
statements are those of [P3]'s Theorem 2 (the $\theta_{ij}$ involve
$f'(1),\dots,f^{(i+j)}(1)$; $T_k$ only $\theta_{ij}$ with $i+j\le k-1$). $\square$

## Corollary (the displayed $\mathcal J_1,\mathcal J_2,\mathcal J_3$) — PROVED (exact algebra)

Specializing Theorem 1 to $K=3$ gives **exactly** the paper's displayed formulas:
$$\mathcal J_1=\frac{2z}{a}\mu_1,\qquad
\mathcal J_2=\frac z{a^3}\Bigl[A(\mu_2-2\mu_1)-a^2\mu_2+2a(1-z)\Bigr],$$
$$\mathcal J_3=\frac z{3a^5}\Bigl[a^4\mu_3-(3Aa^2+Ba)(\mu_3-3\mu_2)+3A^2(\mu_3-4\mu_2+2\mu_1)\Bigr]
-\frac{4Az(1-z)}{a^4}-\frac{cz(1-z)}{a^3},$$
with $2-L_n=\varepsilon\mathcal J_1+\varepsilon^2\mathcal J_2+\varepsilon^3\mathcal J_3
+O(\varepsilon^4\log^4n)$.

*Proof.* The $\theta_{ij}$ for $i+j\le3$ are
$\theta_{10}=\frac1a$, $\theta_{20}=-\frac{a^2-A}{2a^3}$, $\theta_{11}=-\frac A{a^3}$,
$\theta_{30}=\frac{a^4-3a^2A-Ba+3A^2}{6a^5}$,
$\theta_{21}=\frac{3a^2A+Ba-4A^2}{2a^5}$, $\theta_{12}=\frac{A^2}{a^5}$
(`step_30`, agreeing with [P3]). With $M_i=2z\mu_i$ and the tail terms
$$T_2=\frac{d_2}{a^2}=\frac{2z(1-z)}{a^2},\qquad
T_3=\frac{d_2}{a}\Bigl(-2\theta_{20}+\theta_{11}-\theta_{10}\Bigr)+\frac{d_3}{4a^2}\theta_{10}
=-\frac{4Az(1-z)}{a^4}-\frac{cz(1-z)}{a^3},$$
the sum $\sum\theta_{i,k-i}M_i+T_k$ simplifies to the displayed $\mathcal J_k$: sympy gives
difference $0$ at general $(c,z)$ with $2z=1-c$ (`step_31`). At $\theta=\pi/2$
($c=0$, $z=\frac12$) they reduce to [P3]'s proved $Q_1,Q_2,Q_3$ (`step_32`, difference $0$).
For $f(x)=x$ they agree with the exact every-angle series (`step_32`, and independently
re-derived from the layer-cake integral and the CGF in `step_33`: difference $0$). $\square$

*Remark ($\zeta(3)$).* $\zeta(3)$ enters $\mathcal J_3$ only through
$k_3(\theta)=2+2w_2^3\log^3\cos^2\frac\theta2-12w_2\mathrm S_{1,2}(z)$, i.e. through the
Nielsen value $\mathrm S_{1,2}(\sin^2\frac\theta2)$ with coefficient $-12\cot^2\frac\theta2$,
exactly as in the corner coefficient $a_3$ of Thm. P2asym.

*Remark (order 4).* The same specialization at $K=4$ gives
$\theta_{40},\theta_{31},\theta_{22},\theta_{13}$ agreeing with [P3] exactly (`step_40`)
and $T_4\big|_{c=0}$ agreeing with [P3]'s $T_4$ (`step_41`); numerically the residual after
$\mathcal J_3$ equals $\varepsilon^4\mathcal J_4(1+o(1))$ (`step_42`), which pins down the
shape of the error term $O(\varepsilon^4\log^4n)$.

## Theorem 2 (contact of order $m$ at every angle) — PROVED

Assume $f$ as in §0 except for the value of $f'(1)$, and that near $x=1$,
$$\Phi(\tau):=-\log f(1-\tau)=\beta\tau^m\bigl(1+O(\tau^\delta)\bigr),\qquad
\Phi'(\tau)=m\beta\tau^{m-1}\bigl(1+O(\tau^\delta)\bigr),$$
for some $m>0$, $m\ne1$, $\beta>0$, $\delta>0$ (so $f'(1)=0$ if $m>1$, $f'(1)=\infty$ if
$m<1$). Put $N=(\beta n)^{1/m}$, $\Lambda=\log(mN)$, $c_m=1-\frac1m$,
$y_m=c_m(1-\log c_m)$; for $m>1$ let $U(y)>c_m>V(y)$ be the two roots of $u-c_m\log u=y$
(Lambert $W$), and for $m<1$ let $V:=0$, $y_m:=-\infty$. Put $F:=U^{1/m}$. Then

**(a) (exact reduction)**
$$2-L_n=\frac{1-c}{N}\,
\mathbb E\Bigl[\bigl(U^{1/m}-V^{1/m}\bigr)(\Lambda+\xi_\theta)\,;\,\Lambda+\xi_\theta>y_m\Bigr]
\Bigl(1+O\bigl(n^{-\delta/m}\log^{1+\delta/m}n\bigr)\Bigr),$$
with $\delta_N:=C(\log N)^{1+\delta/m}N^{-\delta}$; the sharper relative error
$O(\delta_N/\Lambda)$ of [P3], step 4, holds as well, and for $\delta=m$ (e.g.
$f=\sin\frac{\pi x}2$) it is $O(\log n/n)$.

**(b) (expansion)**
$$2-L_n=\frac{1-c}{N}\Bigl[F(\Lambda+k_1(\theta))+\frac{k_2(\theta)}2F''(\Lambda+k_1(\theta))
+O(\Lambda^{1/m-3})\Bigr]\Bigl(1+O(\delta_N)\Bigr).$$
Equivalently $2-L_n=(1-c)\,\tau_n^+(e^{-k_1(\theta)})\,
\Bigl[1-\frac{(m-1)k_2(\theta)}{2m^2\Lambda^2}+O\bigl(\frac{\log\Lambda}{\Lambda^3}\bigr)\Bigr]$,
where $\tau_n^+(s)$ is the distance from $1$ of the outer point where $(f^n)'=s$.

**(c) (leading law + refined law)**
$$2-L_n=\frac{(1-c)\,\Lambda^{1/m}}{N}\Bigl[1+\frac{(m-1)\log\Lambda+m\,k_1(\theta)}{m^2\Lambda}
+O\Bigl(\frac{\log^2\Lambda}{\Lambda^2}\Bigr)\Bigr]
\;\sim\;(1-\cos\theta)\Bigl(\frac{\log n}{m\beta n}\Bigr)^{1/m}.$$
Since $\Lambda=\frac{\log n+\log(m^m\beta)}m$, the refined law in terms of $\log n$ is
$$\frac{2-L_n}{(1-\cos\theta)\bigl(\log n/(m\beta n)\bigr)^{1/m}}
=1+\frac{(m-1)\log\log n+m\,k_1(\theta)+\log(m\beta)}{m\log n}
+O\Bigl(\Bigl(\frac{\log\log n}{\log n}\Bigr)^2\Bigr).$$
The approach to the leading law is **from above** (for $m>1$) and only logarithmic: this is
the "ratio $1.12$–$1.14$ at $n=10^7$, falling slowly to $1$" of the paper's Table II.

*Proof.* This is [P3]'s proof of its Theorem 4, with the same three replacements as in
Theorem 1. The **$\theta$-audit** of [P3]'s steps 1–6:

| step of [P3] | content | $\theta$-input | replacement |
|---|---|---|---|
| 1 | scaled slope $y_n'(1-\sigma/N)=NH(\sigma)(1+\epsilon_n)$ | none | unchanged |
| 2 | cut-offs: $s<N^{-10}$ contributes $\le N^{-10}$; window inside $\sigma<\sigma_{\max}$ | density bounded near $0$ | $\phi_\theta\le1-c$ (Lemma B(i)) |
| 3 | sandwich $m_n(s)$ between $\frac1N\ell(\Lambda+\log\frac1s\mp\cdots)$ | none | unchanged |
| 4 | $\mathcal F(\Lambda)=\mathbb E\ell(\Lambda+\xi)$; $\mathcal F'=O(\Lambda^{1/m-1})$; the singularity of $\ell$ at $y_m$ meets the density of $\xi$ where it is $O(e^{-\Lambda})$ | left tail $\mathbb P(\xi>t)\le e^{-t}$ and density of $S$ at $\infty$ $O(s^{-2})$ | Lemma C(i): $\mathbb P(S_\theta\le s)\le s$; Lemma B(ii): $\phi_\theta=O(s^{-2})$ |
| 5 | inner root $V^{1/m}\le Ce^{-y/(m-1)}$ gives $O(\Lambda N^{-\min(1,1/(m-1))})$ | density of $S$ at $\infty$ $O(s^{-2})$ | Lemma B(ii) |
| 6 | Taylor-expand $F$ at $\Lambda+k_1$ on $\{\xi>-\Lambda/2\}$; complement has probability $\le e^{-\Lambda/2}$; $\mathbb E(\xi-k_1)=0$ | first two moments + exp. left tail of $\xi$ | Lemma C(i),(ii) |

The only structural change is the mass factor: the density is $\phi_\theta$ with total mass
$1-c$ (not $1$), so every expectation carries the factor $1-c$; the cumulants are
$k_j(\theta)$ instead of $k_j$. The uniform bound
$N\,D_n\in(1-c)\,[\mathcal F(\Lambda-2\delta_N),\mathcal F(\Lambda+2\delta_N)]+O(N^{-9})$
of step 4 then reads, with $\mathcal F(\Lambda)=\mathbb E\,\ell(\Lambda+\xi_\theta)$,
$$2-L_n=\frac{1-c}{N}\,\mathcal F\bigl(\Lambda+O(\delta_N)\bigr)+O(N^{-10})
=\frac{1-c}{N}\,\mathcal F(\Lambda)\bigl(1+O(\delta_N/\Lambda)\bigr),$$
since $\mathcal F'/\mathcal F=O(\Lambda^{-1})$. This is (a). Step 6 gives (b):
$\mathcal F(\Lambda)=F(\Lambda+k_1(\theta))+\frac{k_2(\theta)}2F''(\Lambda+k_1(\theta))
+O(\Lambda^{1/m-3})$, because $\mathbb E[\xi_\theta-k_1]=0$,
$\mathbb E(\xi_\theta-k_1)^2=k_2(\theta)$, the Taylor series of $F$ has
$F^{(j)}(y)=O(y^{1/m-j})$, and the complement
$\{\xi_\theta\le-\Lambda/2\}$ has probability $\mathbb P(S_\theta\ge e^{\Lambda/2})
=O(e^{-\Lambda/2})$ while $F=O(\Lambda^{1/m})$. For (c): $U=y+c_m\log y+O(\frac{\log y}y)$
(solve $u-c_m\log u=y$ by one iteration), so
$F(y)=y^{1/m}\bigl(1+\frac{(m-1)\log y}{m^2y}+O\bigl(\frac{\log^2y}{y^2}\bigr)\bigr)$, and
$F(\Lambda+k_1)=\Lambda^{1/m}\bigl(1+\frac{(m-1)\log\Lambda+m k_1}{m^2\Lambda}
+O(\frac{\log^2\Lambda}{\Lambda^2})\bigr)$; the $k_2$ term of (b) is
$O(\Lambda^{1/m-2})$, absorbed into the $O(\log^2\Lambda/\Lambda^2)$ remainder. The
$\log n$ form follows from $\Lambda=\frac{\log n+\log(m^m\beta)}m$ and
$\log(m^m\beta)-(m-1)\log m=\log(m\beta)$. $\square$

## Proposition (the $u^\ast$ correction; identifies the constant of AH3's step\_29) — PROVED

For the $m=2$ model integral $I(V,c)=\int_0^\infty h_c\bigl(V u e^{-u^2}\bigr)\,du$
($=N\cdot(2-L_n$) for the pure-power model, $V=2N$) and $u^\ast$ defined by
$V u^\ast e^{-(u^\ast)^2}=1$ (equivalently $u^\ast=F(\Lambda)$, $\Lambda=\log V$):
$$\bigl(I(V,c)-(1-c)u^\ast\bigr)u^\ast\longrightarrow\frac{(1-c)\,k_1(\theta)}2
\qquad(V\to\infty).$$

*Proof.* $I(V,c)=(1-c)\,\mathbb E\,\ell(\Lambda+\xi_\theta)$ exactly (Lemma A with
$m(s)=\ell(\Lambda+\log\frac1s)$, $\ell=F-V^{1/2}$; see step 3 of Theorem 2's audit), so by
Theorem 2(b), $I-(1-c)u^\ast=(1-c)\bigl[k_1F'(\Lambda)+O(F''(\Lambda))\bigr]$, and
$F'(\Lambda)u^\ast=\frac{(u^\ast)^2}{2(u^\ast)^2-1}\to\frac12$ while
$F''(\Lambda)u^\ast=O((u^\ast)^{-2})\to0$. $\square$

Consistently, the constant $\mathrm{Klayer}(c)=\frac12\bigl(\int_0^1\frac{h_c(s)}s\,ds
-\int_1^\infty\frac{\rho_c(s)}s\,ds\bigr)$ of AH3 equals $\frac{(1-c)k_1(\theta)}2$
**exactly** (integration by parts: $\int_0^1\frac{h_c}s=-\int_0^1\phi_\theta\log s$ and
$-\int_1^\infty\frac{\rho_c}s=-\int_1^\infty\phi_\theta\log s$; sum
$=-\int_0^\infty\phi_\theta\log s=(1-c)\,\mathbb E\log\frac1{S_\theta}=(1-c)k_1$). AH3's
predicted constant $\mathrm{Klayer}(c)-\frac{1-c^2}4$ is therefore off by $\frac{1-c^2}4$;
numerically the limit is $(1-c)k_1/2$, not $\mathrm{Klayer}(c)-\frac{1-c^2}4$ (`step_45`).

## Numerical evidence (every number printed by code; scripts listed in the ledger)

**T1 (residual after $\mathcal J_3$, scaled by $\varepsilon^4\log^4n$; `step_39`).**
$f=\frac{e^x-1}{e-1}$ and $f=\frac{x+x^2}2$ at $\theta=\frac\pi3,\frac{2\pi}3$ (control
$\frac\pi2$), $n=10^4,10^5,10^6$; quadrature error estimates $\le10^{-69}$:

| $f$ | $\theta$ | scaled@$10^4$ | scaled@$10^5$ | scaled@$10^6$ |
|---|---|---|---|---|
| exp | $\pi/3$ | $0.00085$ | $0.00134$ | $0.00159$ |
| exp | $2\pi/3$ | $0.00299$ | $0.00454$ | $0.00524$ |
| exp | $\pi/2$ | $0.00181$ | $0.00282$ | $0.00330$ |
| quad | $\pi/3$ | $-0.00662$ | $-0.00526$ | $-0.00447$ |
| quad | $2\pi/3$ | $-0.02189$ | $-0.01716$ | $-0.01446$ |
| quad | $\pi/2$ | $-0.01382$ | $-0.01091$ | $-0.00923$ |

Bounded as required; in fact the residual equals $\varepsilon^4\mathcal J_4(1+o(1))$:
ratios resid$/(\varepsilon^4\mathcal J_4)=0.99982,0.99998,0.999998$ at $10^4,10^5,10^6$
(`step_42`), and the new residual is $O(\varepsilon^5\log^5n)$ with scaled values $\le0.005$.

**T2 (contact, $f=\sin\frac{\pi x}2$, $m=2$, $\beta=\pi^2/8$; `step_43`).**
Ratios $D/$law at $\theta=\frac\pi3,\frac\pi2,\frac{2\pi}3$, $n=10^4,\dots,10^8$:

| law | content | $\pi/3$ @$10^5$ | $\pi/2$ @$10^5$ | $2\pi/3$ @$10^5$ | @ $n=10^8$ (3 angles) |
|---|---|---|---|---|---|
| $L_0$ | leading $(1-c)\sqrt{\log n/(2\beta n)}$ | $1.1435$ | $1.1598$ | $1.1812$ | $1.1064,1.1157,1.1282$ |
| $L_1$ | $(1-c)\,u^\ast/N$ (AH3 refined) | $0.9974$ | $1.0116$ | $1.0303$ | $1.0026,1.0110,1.0223$ |
| $L_3$ | (b): $F(\Lambda+k_1)+\frac{k_2}2F''$ | $0.9961$ | $0.9972$ | $0.9984$ | $0.99970,0.99979,0.99990$ |
| full | (a): Lambert-$W$ model | $0.99999$ | $0.99999$ | $0.99999$ | $0.99999998$ (all three) |

$D/L_0$ at $n=10^5$ reproduces AH3's $1.1435,1.1598,1.1812$ exactly, and $D/L_1$ reproduces
$0.9974,1.0116,1.0303$ exactly (cross-validation of the engines). The full-model error
decays like $\log n/n$ ($8.6{\cdot}10^{-5},1.1{\cdot}10^{-5},1.3{\cdot}10^{-6},1.4{\cdot}
10^{-7},2{\cdot}10^{-8}$ at $n=10^4\dots10^8$), matching (a)'s sharper bound. The refined
leading law (c) is quantitative: $(D/L_0-\mathrm{pred})(\Lambda/\log\Lambda)^2$ stays in
$[0.70,0.85]$ over $n=10^4\dots10^8$ (`step_44`). $L_1$'s residual is the missing
$k_1(\theta)$-shift of (b); $L_3$'s error is $O(\Lambda^{-3})$ relative, as proved.

**T3 (moment/Mellin checks).** Cumulants $k_{1,2,3}(\theta)$ vs closed forms, worst
residual $6.1{\cdot}10^{-60}$ over 7 angles; Mellin vs Legendre, worst $7.5{\cdot}10^{-60}$
at $w=\pm0.95,\pm0.5,0.3$, 4 angles (`step_29`). Tail coefficients
$d_l(c)=(-1)^l(1-c^2)P_{l-1}'(c)/l$ symbolic for $l\le8$, remainder constants $<0.4$
(`step_25`,`step_26`).

## Lemma ledger

| label | statement | status | proof / evidence |
|---|---|---|---|
| Lemma A | layer cake $2-L_n=\int\phi_\theta m_y=(1-c)\mathbb E m_y(S_\theta)$ | PROVED | Tonelli; exact |
| Lemma B(i) | $\phi_\theta$ alt form, decreasing, $0<\phi_\theta\le1-c$, mass $1-c$ | PROVED | rationalization; sympy-exact (`step_25`) |
| Lemma B(ii) | $\phi_\theta(s)\le\frac{1+|c|}{2(1-|c|)}s^{-2}$, $s\ge1$ | PROVED | displayed algebra; numeric (`step_26`) |
| Lemma B(iii) | tail series, $d_l(c)=(-1)^l(1-c^2)P_{l-1}'(c)/l$, remainder | PROVED | Legendre generating function + identity; sympy-exact (`step_25`); remainder constants numeric (`step_26`) |
| Lemma C(i) | tails of $S_\theta$: $\mathbb P(S_\theta\le s)\le s$, $\mathbb P(S_\theta>s)=O(s^{-1})$; $\mathbb E S_\theta^w<\infty$ iff $-1<w<1$; exp moments of $\xi_\theta$ for $|t|<1$ | PROVED | from B(i),(ii) |
| Lemma C(ii) | Mellin = Legendre; CGF $=\sum a_jw^j$; $k_j(\theta)=(-1)^jj!a_j$; closed forms $k_{1,2,3}$ | PROVED | paper's Thm. P2 (proved) + differentiation; numeric $6.1{\cdot}10^{-60}$ (`step_29`) |
| Lemma D | window $\tau_n(s)$, cutoffs | PROVED | [P3] Prop. 2 verbatim + B(i) for the $\sigma_n$-bound |
| Lemma E | bivariate inverse $\Theta$, uniform Taylor remainder | PROVED | [P3] Theorem 2's setup (IFT) |
| Theorem 1 | $\sum_{k\le K}\varepsilon^k\mathcal J_k+O(\varepsilon^{K+1}\log^{K+1}n)$ at every $\theta$ | PROVED | §Theorem 1 + $\theta$-audit table |
| Corollary | displayed $\mathcal J_{1,2,3}$; reduce to [P3] at $\pi/2$; exact series at $f=x$ | PROVED | sympy difference $0$ (`step_31`,`step_32`,`step_33`) |
| Theorem 2 | contact: (a) model, (b) expansion, (c) refined leading law | PROVED | §Theorem 2 + $\theta$-audit table |
| Proposition | $(I-(1-c)u^\ast)u^\ast\to(1-c)k_1/2$; Klayer identity | PROVED | Thm 2(b) + integration by parts; numeric (`step_45`) |
| T1 | scaled residuals bounded; resid $=\varepsilon^4\mathcal J_4(1+o(1))$ | NUMERICAL | `step_39`, `step_42` |
| T2 | contact ratios vs $L_0,L_1,L_3$, full model; $n\le10^8$ | NUMERICAL | `step_43`, `step_44` |
| T3 | cumulants/Mellin | NUMERICAL | `step_29` |

Nothing is OPEN. The $\theta=\pi/2$ inputs ([P3] Props. 1,2,6, Thms. 1,2,4 and the paper's
Thm. P2/P2asym) are themselves PROVED.

## Proposed statement for the paper (replaces "derived")

> **Theorem (fn-expansion at every angle).** *Let $u,v$ be unit vectors with
> $u\cdot v=\cos\theta$, $0<\theta<\pi$, and let $f\in C^1[0,1]$ be nondecreasing,
> $f(0)=0$, $f(1)=1$, $a=f'(1)>0$, with $f\in C^{K+2}$ on some $[1-\eta,1]$. Let $L_n$ be
> the length of $x\,u+f(x)^n\,v$ and let $k_j(\theta)$ be the cumulants of $\log(1/S_\theta)$
> above. Then with $\varepsilon=\frac1{n-1}$ and $\Xi=\log(na)+k_1(\theta)$,*
> $$2-L_n=\sum_{k=1}^K\varepsilon^k\mathcal J_k(\Xi)+O\bigl(\varepsilon^{K+1}\log^{K+1}n\bigr),$$
> *where $\mathcal J_k$ is a polynomial of degree $k$ in $\Xi$ depending only on
> $f'(1),\dots,f^{(k)}(1)$ and $\theta$; explicitly
> $\mathcal J_k=\sum_{i=1}^k\theta_{i,k-i}M_i(\Xi)+T_k$ with $M_i=2z\sum_r\binom irm_r(\theta)
> \Xi^{i-r}$ and $T_k$ the slope-tail term with $d_l(c)=(-1)^l(1-c^2)P_{l-1}'(c)/l$. In
> particular $\mathcal J_1,\mathcal J_2,\mathcal J_3$ are the displayed formulas of this
> section (with $\zeta(3)$ in $\mathcal J_3$ through
> $k_3(\theta)=2+2\cot^6\frac\theta2\log^3\cos^2\frac\theta2-12\cot^2\frac\theta2
> \mathrm S_{1,2}(\sin^2\frac\theta2)$).*
>
> *Proof.* The proof of Theorem fn ($\theta=\frac\pi2$) applies verbatim once three
> $\theta$-dependent inputs are replaced: (i) the layer-cake density is $\phi_\theta=h_c'$ of
> mass $1-\cos\theta$ (Lemma A); (ii) all moments of $\log S_\theta$ are finite and its
> cumulants are $k_j(\theta)$, because $\mathbb P(S_\theta\le s)\le s$ and
> $\phi_\theta=O(s^{-2})$ (Lemmas B, C); (iii) the tail series of the density is
> $\phi_\theta(s)=\sum_{l\ge2}d_l(c)s^{-l}$, convergent for $s>1$ (Lemma B). The window
> $\{y'>s\}$ and the bivariate inverse $\Theta$ do not depend on $\theta$. $\square$

> **Proposition (contact at every angle).** *If $-\log f(1-\tau)=\beta\tau^m(1+O(\tau^\delta))$
> with derivative $m\beta\tau^{m-1}(1+O(\tau^\delta))$, $m>0$, $m\ne1$, then with
> $N=(\beta n)^{1/m}$, $\Lambda=\log(mN)$, and $F=U^{1/m}$ for $U$ the outer root of
> $u-(1-\frac1m)\log u=y$,*
> $$2-L_n=\frac{1-\cos\theta}{N}\Bigl[F(\Lambda+k_1(\theta))+\frac{k_2(\theta)}2
> F''(\Lambda+k_1(\theta))+O(\Lambda^{1/m-3})\Bigr]\bigl(1+O(n^{-\delta/m}\log^{1+\delta/m}n)\bigr),$$
> *and in particular*
> $$2-L_n=(1-\cos\theta)\Bigl(\frac{\log n}{m\beta n}\Bigr)^{1/m}
> \Bigl[1+\frac{(m-1)\log\log n+m\,k_1(\theta)+\log(m\beta)}{m\log n}
> +O\Bigl(\frac{\log\log n}{\log n}\Bigr)^2\Bigr].$$
> *The angle enters the leading law only through the factor $1-\cos\theta$; the approach to
> it is from above (for $m>1$) and logarithmically slow.*

Both statements are proved in this note (Theorem 1, Corollary, Theorem 2, and the two
$\theta$-audit tables); every constant was checked numerically (Tables T1–T3).
