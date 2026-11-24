part3 = r"""
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
with the sharper relative error $O(\delta_N/\Lambda)$ of [P3], step 4, valid as well; for
$\delta=m$ (e.g. $f=\sin\frac{\pi x}2$) this is $O(\log n/n)$.

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
numerically the limit is $(1-c)k_1/2$, not $\mathrm{Klayer}(c)-\frac{1-c^2}4$ (`step_43`).

## Numerical evidence (every number printed by code; scripts listed in the ledger)

**T1 (residual after $\mathcal J_3$, scaled by $\varepsilon^4\log^4n$; `step_37`).**
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
(`step_40`), and the new residual is $O(\varepsilon^5\log^5n)$ with scaled values $\le0.005$.

**T2 (contact, $f=\sin\frac{\pi x}2$, $m=2$, $\beta=\pi^2/8$; `step_42`).**
Ratios $D/$law at $\theta=\frac\pi3,\frac\pi2,\frac{2\pi}3$, $n=10^4,\dots,10^8$:

| law | content | $\pi/3$ @$10^5$ | $\pi/2$ @$10^5$ | $2\pi/3$ @$10^5$ | @ $n=10^8$ (3 angles) |
|---|---|---|---|---|---|
| $L_0$ | leading $(1-c)\sqrt{\log n/(2\beta n)}$ | $1.1435$ | $1.1598$ | $1.1812$ | $1.1064,1.1157,1.1282$ |
| $L_1$ | $(1-c)\,u^\ast/N$ (AH3 refined) | $0.9974$ | $1.0116$ | $1.0303$ | $1.0026,1.0110,1.0223$ |
| $L_3$ | (b): $F(\Lambda+k_1)+\frac{k_2}2F''$ | $0.9961$ | $0.9972$ | $0.9984$ | $0.99970,0.99979,0.99990$ |
| full | (a): Lambert-$W$ model | $0.99999$ | $0.99999$ | $0.99999$ | $1.0000000$ (8 digits) |

$D/L_0$ at $n=10^5$ reproduces AH3's $1.1435,1.1598,1.1812$ exactly, and $D/L_1$ reproduces
$0.9974,1.0116,1.0303$ exactly (cross-validation of the engines). The full-model error
decays like $\log n/n$ ($8.6{\cdot}10^{-5},1.1{\cdot}10^{-5},1.3{\cdot}10^{-6},1.4{\cdot}
10^{-7},2{\cdot}10^{-8}$ at $n=10^4\dots10^8$), matching (a)'s sharper bound. The refined
leading law (c) is quantitative: $(D/L_0-\mathrm{pred})(\Lambda/\log\Lambda)^2$ stays in
$[0.70,0.85]$ over $n=10^4\dots10^8$ (`step_44`). $L_1$'s residual is the missing
$k_1(\theta)$-shift of (b); $L_3$'s error is $O(\Lambda^{-3})$ relative, as proved.

**T3 (moment/Mellin checks).** Cumulants $k_{1,2,3}(\theta)$ vs closed forms, worst
residual $6.1{\cdot}10^{-60}$ over 7 angles; Mellin vs Legendre, worst $7.5{\cdot}10^{-60}$
at $w=\pm0.95,\pm0.5,0.3$, 4 angles (`step_28`). Tail coefficients
$d_l(c)=(-1)^l(1-c^2)P_{l-1}'(c)/l$ symbolic for $l\le8$, remainder constants $<0.4$
(`step_24`,`step_25`).
"""
with open("proof.md","a") as f:
    f.write(part3)
print("part3 written", len(part3))
