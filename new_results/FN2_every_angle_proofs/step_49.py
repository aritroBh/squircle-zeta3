part4 = r"""
## Lemma ledger

| label | statement | status | proof / evidence |
|---|---|---|---|
| Lemma A | layer cake $2-L_n=\int\phi_\theta m_y=(1-c)\mathbb E m_y(S_\theta)$ | PROVED | Tonelli; exact |
| Lemma B(i) | $\phi_\theta$ alt form, decreasing, $0<\phi_\theta\le1-c$, mass $1-c$ | PROVED | rationalization; sympy-exact (`step_24`) |
| Lemma B(ii) | $\phi_\theta(s)\le\frac{1+|c|}{2(1-|c|)}s^{-2}$, $s\ge1$ | PROVED | displayed algebra; numeric (`step_25`) |
| Lemma B(iii) | tail series, $d_l(c)=(-1)^l(1-c^2)P_{l-1}'(c)/l$, remainder | PROVED | Legendre generating function + identity; sympy-exact (`step_24`); remainder constants numeric (`step_25`) |
| Lemma C(i) | tails of $S_\theta$: $\mathbb P(S_\theta\le s)\le s$, $\mathbb P(S_\theta>s)=O(s^{-1})$; $\mathbb E S_\theta^w<\infty$ iff $-1<w<1$; exp moments of $\xi_\theta$ for $|t|<1$ | PROVED | from B(i),(ii) |
| Lemma C(ii) | Mellin = Legendre; CGF $=\sum a_jw^j$; $k_j(\theta)=(-1)^jj!a_j$; closed forms $k_{1,2,3}$ | PROVED | paper's Thm. P2 (proved) + differentiation; numeric $6.1{\cdot}10^{-60}$ (`step_28`) |
| Lemma D | window $\tau_n(s)$, cutoffs | PROVED | [P3] Prop. 2 verbatim + B(i) for the $\sigma_n$-bound |
| Lemma E | bivariate inverse $\Theta$, uniform Taylor remainder | PROVED | [P3] Theorem 2's setup (IFT) |
| Theorem 1 | $\sum_{k\le K}\varepsilon^k\mathcal J_k+O(\varepsilon^{K+1}\log^{K+1}n)$ at every $\theta$ | PROVED | §Theorem 1 + $\theta$-audit table |
| Corollary | displayed $\mathcal J_{1,2,3}$; reduce to [P3] at $\pi/2$; exact series at $f=x$ | PROVED | sympy difference $0$ (`step_30`,`step_31`,`step_32`) |
| Theorem 2 | contact: (a) model, (b) expansion, (c) refined leading law | PROVED | §Theorem 2 + $\theta$-audit table |
| Proposition | $(I-(1-c)u^\ast)u^\ast\to(1-c)k_1/2$; Klayer identity | PROVED | Thm 2(b) + integration by parts; numeric (`step_43`) |
| T1 | scaled residuals bounded; resid $=\varepsilon^4\mathcal J_4(1+o(1))$ | NUMERICAL | `step_37`, `step_40` |
| T2 | contact ratios vs $L_0,L_1,L_3$, full model; $n\le10^8$ | NUMERICAL | `step_42`, `step_44` |
| T3 | cumulants/Mellin | NUMERICAL | `step_28` |

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
"""
with open("proof.md","a") as f:
    f.write(part4)
print("part4 written", len(part4))
import os
print("proof.md total:", os.path.getsize("proof.md"))
