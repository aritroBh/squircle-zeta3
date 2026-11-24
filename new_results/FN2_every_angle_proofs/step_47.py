part2 = r"""
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
(`step_29`, agreeing with [P3]). With $M_i=2z\mu_i$ and the tail terms
$$T_2=\frac{d_2}{a^2}=\frac{2z(1-z)}{a^2},\qquad
T_3=\frac{d_2}{a}\Bigl(2\theta_{20}-\theta_{11}+\theta_{10}\Bigr)+\frac{d_3}{4a^2}\theta_{10}
=-\frac{4Az(1-z)}{a^4}-\frac{cz(1-z)}{a^3},$$
the sum $\sum\theta_{i,k-i}M_i+T_k$ simplifies to the displayed $\mathcal J_k$: sympy gives
difference $0$ at general $(c,z)$ with $2z=1-c$ (`step_30`). At $\theta=\pi/2$
($c=0$, $z=\frac12$) they reduce to [P3]'s proved $Q_1,Q_2,Q_3$ (`step_31`, difference $0$).
For $f(x)=x$ they agree with the exact every-angle series (`step_31`, and independently
re-derived from the layer-cake integral and the CGF in `step_32`: difference $0$). $\square$

*Remark ($\zeta(3)$).* $\zeta(3)$ enters $\mathcal J_3$ only through
$k_3(\theta)=2+2w_2^3\log^3\cos^2\frac\theta2-12w_2\mathrm S_{1,2}(z)$, i.e. through the
Nielsen value $\mathrm S_{1,2}(\sin^2\frac\theta2)$ with coefficient $-12\cot^2\frac\theta2$,
exactly as in the corner coefficient $a_3$ of Thm. P2asym.

*Remark (order 4).* The same specialization at $K=4$ gives
$\theta_{40},\theta_{31},\theta_{22},\theta_{13}$ agreeing with [P3] exactly (`step_38`)
and $T_4\big|_{c=0}$ agreeing with [P3]'s $T_4$ (`step_39`); numerically the residual after
$\mathcal J_3$ equals $\varepsilon^4\mathcal J_4(1+o(1))$ (`step_40`), which pins down the
shape of the error term $O(\varepsilon^4\log^4n)$.
"""
with open("proof.md","a") as f:
    f.write(part2)
print("part2 written", len(part2))
