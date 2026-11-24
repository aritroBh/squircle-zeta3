# D16 — The pure-weight theorem near the diamond

**Folder:** `research/08_theory/D16_pure_weight/`. Every number below was printed by a script run in this session (auto-saved as `step_NN.py`; artifacts on disk). Statuses **PROVED / NUMERICAL / CONJECTURE** marked per claim. Notation: $p\ge1$, $\tau=1-1/p$ ($=0$ at the diamond), $d=p-1=\tau/(1-\tau)$, $Z_k$ = formal symbol for $\zeta(k)$ (weight $k$; $\zeta(1)$ never occurs, $\pi^2=6\zeta(2)$).

---

## (1) The pure-weight theorem — PROVED

**Theorem PW.** Let $S_n(p)$ be the surface $((n-1)$-volume) of $\{\sum_i|x_i|^p=1\}$ in $\mathbb{R}^n$, $S_n(1)=2^n\sqrt n/(n-1)!$ the cross-polytope, and define the elementary prefactor
$$\mathrm{Rat}_n(\tau)=(1-\tau)^{\,n-1}\Big/\prod_{i=1}^{n-1}\bigl(1-\tfrac{(n-1)\tau}{i}\bigr)=\prod_{i=1}^{n-2}\frac{i}{i-(n-1-i)d}.$$
Then
$$\frac{S_n(p)}{S_n(1)\,\mathrm{Rat}_n(\tau)}=\exp\Big(\sum_{k\ge2} w_k(n)\,\tau^k\Big),\qquad w_k(n)\in\mathbb Q[\zeta(2),\dots,\zeta(k)]\ \text{homogeneous of weight }k .$$

**Proof (mechanism).** Start from H1's representation (T1, classical Naor–Romik/Schechtman–Zinn): with $\varepsilon=1/p=1-\tau$, $a=-(n-1)\tau/n$,
$$S_n/S_n(1)=(n-1)!\,\varepsilon^{\,n-1}\frac{\Gamma(1+a)^n}{\Gamma(1+(n-1)\varepsilon)}\cdot B_n(\tau),\quad B_n=\mathbb E\Big[\Big(\tfrac1n\sum_i e^{2\tau m_i}\Big)^{1/2}\Big],\quad m_i=\log G_i-\overline{\log G},\ G_i\sim\mathrm{Gamma}(1+a).$$
- **M1 (everything is $\Gamma(1+x)$ at $x=\tau\cdot\mathbb Q$).** Expanding the $\tfrac12$-power, $B_n=\sum_j\binom{1/2}{j}\mathbb E[X^j]$, $X=\frac1n\sum_i(e^{2\tau m_i}-1)$. For a tuple with equality-pattern $\lambda\vdash s$ ($r$ distinct parts), independence of the $G_i$ gives the *completely explicit moment*
$$\mathbb E\Big[\prod_l(e^{2\tau m_{i_l}}-1)\Big]=\sum_{S}(-1)^{j-s}\frac{\prod_v\Gamma\bigl(1+a+2\tau(\nu_v-s/n)\bigr)}{\Gamma(1+a)^n},\qquad \nu_v=\#\{l\in S:i_l=v\}.$$
Hence, with $c(\lambda)=2\lambda-\frac{2s}{n}-\frac{n-1}{n}$,
$$\mathbb E[X^j]=n^{-j}\sum_{s=0}^j\tbinom js(-1)^{j-s}n^{\,j-s}\sum_{\lambda\vdash s}\frac{(n)_r}{\prod \mathrm{mult}!}\frac{s!}{\prod\lambda_i!}\;M(\lambda,s),\quad \log M=\sum_{k\ge2}\frac{(-1)^kZ_k}{k}\Big[\sum_i c(\lambda_i)^k+(n-r)c(0)^k-n\Big(-\tfrac{n-1}{n}\Big)^k\Big]\tau^k .$$
- **M2 ($\log\Gamma(1+x)=-\gamma x+\sum_{k\ge2}(-1)^k\zeta(k)x^k/k$).** The index of $\zeta(k)$ always multiplies $x^k$, and every $x$ is $\tau\times(\text{rational in }n)$: so $\zeta(k)$ enters at exactly $\tau^k$ with weight $k$. Products/exp/log preserve the grading (they are built from $\mathbb Q$-linear combinations of products).
- **M3 ($\gamma$ cancels; integer shifts are rational).** The $\gamma$-part of each moment is proportional to $\sum_v(\nu_v-s/n)=0$; in the flat-face factor, $n(-\gamma a)=\gamma(n-1)\tau$ cancels $-\log\Gamma(1-(n-1)\tau)$'s $\gamma(n-1)\tau$. The recurrence $\Gamma(n-(n-1)\tau)=\Gamma(1-(n-1)\tau)\prod_{i=1}^{n-1}(i-(n-1)\tau)$ produces only finite harmonic sums (rational); these assemble with $(n-1)\log(1-\tau)$ exactly into $\log\mathrm{Rat}_n$; the remaining $\zeta$-part is $\Lambda_n=-\sum_{k\ge2}\zeta(k)(n-1)^k(1-n^{1-k})\tau^k/k$. $\square$

So $\log(S_n/(S_n(1)\mathrm{Rat}_n))=\Lambda_n+\log B_n$ and **$w_k(n)=\lambda_k(n)+b_k(n)$** with $\lambda_k$ above and $b_k=[\tau^k]\log B_n$ from the explicit formula. *Exact verification:* my independent implementation of the $B_n$ formula reproduces H1's Route A data ($n=2$, all 13 coefficients through $\tau^{12}$, formal $Z$'s, difference exactly 0) — `step_18.py`, `step_21.py`.

**The right prefactor.**
- **2D:** $\mathrm{Rat}_2=\frac{1-\tau}{1-\tau}=\mathbf 1$ exactly. So *no prefactor is needed*: $\log\!\big(P(p)/(4\sqrt2)\big)$ itself is pure weight. This is why the perimeter (not just a rescaled version) is pure weight in $\tau$.
- **General $n$:** $\mathrm{Rat}_n$ (ratio of linear factors in $\tau$; in $d$ a nontrivial power series with harmonic coefficients $\sigma_k$).
- **3D:** $\mathrm{Rat}_3=\frac{1-\tau}{1-2\tau}=\frac1{1-d}=(2-p)$ — exactly H1's prefactor.
- **Rhombus $|x/a|^p+|y/b|^p=1$ (general semi-axes):** from D11c's operator representation $P=4r(1+\eta_D)\int_0^1(t(1-t))^{\eta_D/2}\sqrt{\cosh\eta_D l-q\sinh\eta_D l}\,dt$ ($r=\sqrt{a^2+b^2}$, $q=(a^2-b^2)/(a^2+b^2)$, $\eta_D=-\tau$), the same $K,T$-machinery with $c_k(q)\in\mathbb Q[q^2]$ gives **$\log(P(a,b;p)/P(a,b;1))=\sum_{k\ge2}w_k(q)\tau^k$, $w_k$ homogeneous of weight $k$ in $\mathbb Q[q^2,\zeta(2..k)]$; prefactor $=4\sqrt{a^2+b^2}=P(1)$, $\mathrm{Rat}=1$.** (PROVED by the identical M1–M3 argument; consistent with D11c's proved $\kappa_3=1-q^2=4a^2b^2/(a^2+b^2)^2$.)

**Exactly where the weight-breaking in $d$ comes from.** Two sources, both displayed:
(i) *Variable change* $\tau=d/(1+d)$: $[d^n]$ receives $[\tau^j]$ for all $j\le n$ with mixing matrix $[d^n]\tau^j=(-1)^{n-j}\binom{n-1}{j-1}$; since $[\tau^j]$ has weight $j$, weights mix downward. Example (printed, `step_24.py`): $a_3=(\zeta(3))-2(\zeta(2))=\zeta(3)-\pi^2/3$ — the $\pi^2/3$ is pure weight-mixing from $\tau^2$.
(ii) *Prefactor* $\mathrm{Rat}_n(d)=\prod i/(i-(n-1-i)d)$, a rational power series in $d$ (weight 0) multiplying everything — e.g. 3D $c_3=1+\pi^2/9-16\zeta(3)/9$, where the rational $1$ and the split of $\pi^2/9$ are $\sigma$-pollution. In $\tau$ both pollutions are removed at once: divide by $\mathrm{Rat}_n(\tau)$ and take $\log$.

**Invariance lemma (first appearance is coordinate-free) — PROVED.** If $A(\tau)=\sum a_k\tau^k$ is pure weight and $\tau=d+O(d^2)$ any change of variables, then $[Z_n][d^n]A=[Z_n][\tau^n]A$: $[Z_n][\tau^j]=0$ for $j<n$ (weight $j<n$ cannot contain the weight-$n$ symbol), and $[d^n]\tau^n=1$. So the $d$- and $\tau$-versions of the first-appearance coefficient are *identical*, not merely similar. $\square$

---

## (2) The 2D perimeter in $\tau$ through $\tau^{12}$, exactly — PROVED, 3 independent routes

Computed from the derived formula $P/\sqrt2=4\,e^{\Lambda_2}B_2$ (my explicit $\Gamma$-product formula, `step_18.py`), cross-checked exactly against H1's Route A (`step_21.py`) and against a rebuild of D11c's operator formula $4K(\eta_D)T(\eta_D)$ (`step_23.py`). All three agree, all 13 coefficients, formal $Z$'s, difference 0.

$$\begin{aligned}
P/\sqrt2=\;4&+\zeta(2)\tau^2+\zeta(3)\tau^3-\tfrac{33}{16}\zeta(4)\tau^4-\bigl(\tfrac34\zeta(5)+\tfrac74\zeta(2)\zeta(3)\bigr)\tau^5+\bigl(\tfrac{577}{34560}\pi^6-\tfrac78\zeta(3)^2\bigr)\tau^6\\
&+\bigl(\tfrac{129}{16}\zeta(7)+\tfrac{47}{32}\pi^2\zeta(5)+\tfrac{119}{640}\pi^4\zeta(3)\bigr)\tau^7+\bigl(-\tfrac{2891923}{116121600}\pi^8+\tfrac{169}{192}\pi^2\zeta(3)^2+\tfrac{141}{16}\zeta(3)\zeta(5)\bigr)\tau^8+\cdots
\end{aligned}$$
(Full formal-$Z$ table through $\tau^{12}$: `tau_series_2d.json`; reduced: `tau_reduced_2d.json`, printed in `step_25.py`. Each $\tau^n$ coefficient = exactly the weight-$n$ part of D11's $a_n$.)

**Log form (the cleanest statement):**
$$\log\frac{P}{4\sqrt2}=\sum_{k\ge2}u_k(2)\zeta(k)\tau^k+(\text{weight-}k\text{ products}),\quad u_k(2)=\tfrac14\kappa_k:\ \tfrac14,\tfrac14,\tfrac1{32},-\tfrac3{16},\tfrac{91}{192},\tfrac{129}{64},-\tfrac{5039}{1024},-\tfrac{22805}{768},\tfrac{561331}{5120},\tfrac{805173}{1024},-\tfrac{94259879}{24576}.$$

**Checks.**
- *vs D11:* $\tau=d/(1+d)$ substituted, even zetas reduced — matches D11's exact $a_0..a_{10}$ **exactly** (`step_26.py`). Bonus, new exact values: $a_{11}=-34589.61588758435467783021\ldots$, $a_{12}=-184166.4588525884605384399\ldots$ (full closed forms in `tau_reduced_2d.json`).
- *independent numeric engine* (cumulant route, 150 digits, 64 orders): worst relative deviation vs H1's independently computed 80-order table = **3.8e−150** (`step_28.py`).
- *vs direct quadrature* (`step_29.py`, 70 dps; support-function integral and arclength integral agree to $\sim$65 digits): $P(1.02)=5.6577609642376152668664441127402268527837186805047\ldots$, $P(1.05)=5.6622954107876893413977912365487084466151094250548\ldots$

| $p$ | $\tau$ | truncation $N$ | abs. error | digits |
|---|---|---|---|---|
| 1.02 | 0.019608 | 12 | 1.1e−16 | 16 |
| 1.02 | | 40 | 2.5e−31 | 31 |
| 1.02 | | 56 | 3.2e−35 | **35** |
| 1.02 | | 64 | 1.9e−36 | 36 |
| 1.05 | 0.047619 | 12 | 7.2e−12 | 12 |
| 1.05 | | 24 | 1.3e−15 | 15 |
| 1.05 | | **32 (optimal)** | **2.0e−17** | **17** |
| 1.05 | | 48 | 5.9e−15 (growing) | — |
| 1.05 | | 64 | 3.6e−11 (growing) | — |

At $p=1.02$ plain truncation gives ≥30 digits as required. At $p=1.05$ the least term is $\sim$2e−17 (least term at $N^*\!\approx\!32$, errors increase past it — the zero radius of convergence made visible), so **no** truncation can give 30 digits there; Borel summation is required (D11's theorem). My quick `mp.pade`-based Borel attempt mis-evaluated (results discarded; implementation bug, not a discrepancy) — labeled NUMERICAL-not-certified, to be redone with D11's certified Padé code.

---

## (3) The first-appearance law in $\tau$ — PROVED

**Single-appearance lemma (PROVED, structural):** in $\log\!\big(S_n/(S_n(1)\mathrm{Rat}_n)\big)$, the single zeta $\zeta(k)$ occurs **exactly once, at $\tau^k$**, with rational-in-$n$ coefficient $u_k(n)$ (everything else at $\tau^N$, $N>k$, involving $\zeta(k)$ is products). Proof: each moment's $Z_k$-part is homogeneous of degree $k$ in the $\tau$-scaled variables, and $B_n|_{Z=0}=1$ (checked: $\mathbb E[X^j]|_{Z=0}=n^{-j}\sum_s\binom js(-1)^{j-s}n^{j-s}n^s=(1-1)^j=0$). $\square$

**Closed form (PROVED, finite formula):** $u_k(n)=-\dfrac{(n-1)^k(1-n^{1-k})}{k}+\dfrac{(-1)^k}{k}\,U_k(n)$, where
$$U_k(n)=\sum_{j=1}^{k}\tbinom{1/2}{j}n^{-j}\!\!\sum_{s=0}^j\tbinom js(-1)^{j-s}n^{j-s}\Big[\sum_{\lambda\vdash s}\tfrac{(n)_r}{\prod\mathrm{mult}!}\tfrac{s!}{\prod\lambda_i!}\Big(\sum_i c(\lambda_i)^k+(n-r)c(0)^k\Big)-n^s\,n\Big(-\tfrac{n-1}{n}\Big)^k\Big]$$
— finite for each $k$ (terms with $j>k$ vanish: the $s$-sum is the $j$-th finite difference of a degree-$k$ polynomial).

**Comparison with D11c's $\kappa_n$ (answer to "is the $\tau$-version simpler?").** For 2D: $[Z_n][\tau^n](P/\sqrt2)=4u_n(2)=\kappa_n$ — **identical** to the $d$-version (by the invariance lemma; verified exactly for $n=3,5,7,9,11$ against $\kappa_n=\frac{4(-1)^n}{n}(1-2^{1-n}V_n)$, $V_n=n![t^n]\sqrt{(1+e^{4t})/2}$; e.g. $4u_7(2)=4\cdot\frac{129}{64}=\frac{129}{16}=\kappa_7$). So the single-zeta coefficient is **not** simpler in $\tau$ — it is the same number; what becomes simpler is *everything else* (the full $\tau^n$ coefficient is pure weight, and in the log each single zeta appears once). For general $n$: $u_k(n)$ is a rational function of $n$ with denominator $k\,n^{k-1}$, factoring nicely:
$$u_2=-\frac{(n-1)\bigl((n-1)^2-2\bigr)}{2n},\qquad u_3=\frac{(n^2-1)\bigl(2-(n-1)^3\bigr)}{3n^2},\qquad b_2=\frac{n-1}{n},\ b_3=\frac{2(n^2-1)}{3n^2},\ b_4=\frac{(n-1)(n^2-3n+4)}{n^3},$$
(all PROVED, and verified to match H1's general-$n$ table for $n=2..12$). Neither a polynomial nor a single binomial in general — it is a rational function with the $\mathrm{Rat}_n$-companion $-(n-1)^k(1-n^{1-k})/k$; closed forms for $b_k(n)$, $k\ge5$, remain open (symbolic-$n$ engine built but the $k\ge5$ runs were not completed in time).

---

## (4) Zeta-free generating functions, and divergence vs pure weight

**Set all $\zeta=0$ (PROVED):** $\log\Gamma(1+x)\to-\gamma x$; by M3 every moment becomes 1, $B_n\equiv1$, $\Lambda_n=0$: the truncated series is the **rational closed form** $S_n/S_n(1)=\mathrm{Rat}_n(\tau)$ (2D: $P=4\sqrt2$ exactly).

**Keep only $\zeta(2)$ (PROVED as formal-series truncation; numeric spot-checks not run):** $\log\Gamma(1+x)\to-\gamma x+\zeta(2)x^2/2$, so every moment is Gaussian: $M(\lambda,s)=e^{2\zeta(2)\tau^2 Q}$, $Q=\sum_i(c(\lambda_i)^2/2)+\cdots=\sum\nu_v^2-s^2/n$. For $n=2$ the log-$G$ difference $l$ becomes $N(0,2\zeta(2))$ and
$$P/\sqrt2\;\big|_{\zeta(2)\text{ only}}=\;4\sqrt2\,e^{-\zeta(2)\tau^2/4}\sum_{k\ge0}c_k\frac{(2k)!}{k!}\,\bigl(\zeta(2)\tau^2\bigr)^k,\qquad \sqrt{\cosh z}=\sum c_kz^{2k},$$
an explicit single sum — the **Borel–Gauss transform $\mathbb E[\sqrt{\cosh(\sqrt{2x}\,Z)}]$** of $\sqrt{\cosh}$, with $c_k(2k)!/k!=1,\tfrac12,-\tfrac18,\tfrac{19}{48},-\tfrac{559}{384},\dots$ (numerators $1,1,-1,19,-559$ over $2^kk!$). **No elementary (Gamma-ratio/Beta) closed form exists** — because $\sqrt{\cosh}$ solves the algebraic ODE $4f^2f'^2=f^4-1$, whose Borel transform is not classical. **Keep only $\zeta(m)$:** $l$ has only the $m$-th cumulant $\kappa_m=2(-1)^m(m-1)!\zeta(m)$, so $\mathbb E[\sqrt{\cosh(\tau l)}]=\sum_{k:\,m\mid2k}c_k\tau^{2k}\frac{(2k)!}{(2k/m)!}\bigl(\kappa_m/m!\bigr)^{2k/m}$ — again explicit single sums, not elementary. **Status: PROVED as truncations of the exact moment formulas; the Gaussian-moment identity is exact algebra; numerical verification not run in the available time (honest label).**

**How pure weight coexists with factorial divergence.** Weight counts *which* products appear (algebraic grading, from $\log\Gamma$); divergence measures *how large the rational coefficients are* (analytic). The coefficients grow factorially: (i) the single-zeta coefficients themselves, $\kappa_n\sim-4(-1)^{(n-1)/2}(2/\pi)^{n}n!\,n^{-3/2}$ (D11c's proved asymptotics), e.g. $u_{12}(2)=-94259879/24576$; (ii) the coefficients of the pure powers $\zeta(2)^k$ in $\log B_2$: $\tfrac12,-\tfrac14,\tfrac12,-\tfrac{41}{24},\tfrac{17}{2},-\tfrac{10067}{180},\dots$ (moments of the logistic law grow like $(2k)!$); (iii) numerically, the term magnitudes $\sim\sqrt2\,k!(2\tau/\pi)^kk^{-3/2}$ bottom out at $k^*\approx\pi/(2\tau)$ and grow after (my $p=1.05$ table: minimum at $k\approx32$). The Borel singularity at $\pm i\pi/2$ is a zero of $\cosh$ — an analytic fact invisible to the weight grading.

---

## (5) Plain-words summary (Aritro-style)

"Measure the rounding of the diamond by $\tau=1-1/p$, the fraction of the way from the diamond toward the square, and take a logarithm. Then every term in the expansion of the perimeter — and of the cross-polytope's area in every dimension — is built from zeta values of exactly the right size: the coefficient of $\tau^n$ is a rational combination of products $\zeta(n_1)\cdots\zeta(n_r)$ with $n_1+\cdots+n_r=n$, nothing younger and nothing older. The reason is a single identity, $\log\Gamma(1+x)=-\gamma x+\sum\zeta(k)x^k/k$, in which $\zeta(k)$ is born already paired with the $k$-th power; Euler's $\gamma$, the only weight-one interloper, cancels between the faces. The price of the wrong variable $d=p-1$ is binomial mixing $\tau=d/(1+d)$ plus a rational prefactor — and that is precisely where $\zeta(3)-\pi^2/3$ gets its $\pi^2/3$. Each odd zeta makes its debut exactly once: $\zeta(n)$ first appears at $\tau^n$, with the same coefficient as in $d$, and in two dimensions that coefficient is one quarter of $(d/dt)^n\sqrt{(1+e^{4t})/2}$ at zero. The series still diverges — but the divergence lives entirely in the rational coefficients, which grow factorially; the zeta values themselves stand at attention, each at its own order."

---

## Evidence table

| claim | status | evidence | digits/residual | script |
|---|---|---|---|---|
| $B_n$ explicit $\Gamma$-product formula; reproduces H1 Route A $n=2$ (13 coeffs) | PROVED + exact check | difference exactly 0, formal $Z$'s | exact | `step_18.py`,`step_21.py` |
| Theorem PW + prefactors $\mathrm{Rat}_2=1$, $\mathrm{Rat}_3=(2-p)$ | PROVED | derivation M1–M3; matches H1's $\mathrm{Rat}_n$, $\Lambda_n$ | exact | `step_20.py` (engine), text |
| $P/\sqrt2$ in $\tau$ through $\tau^{12}$ exact | PROVED | 3 routes agree exactly (my $B_2$ = H1 Route A = D11c operator) | exact | `step_18/21/23.py` |
| $\tau\circ(d/(1+d))$ = D11's $a_0..a_{10}$ | PROVED | exact, all 11 | exact | `step_26.py` |
| new exact $a_{11},a_{12}$ | PROVED (as consequence of proved $\tau$-series) | closed forms saved | exact | `step_25.py` |
| numeric engine vs H1 80-order table | NUMERICAL | worst rel. dev. 3.8e−150, $k\le64$ | 150 | `step_28.py` |
| quadrature $P(1.02),P(1.05)$ | NUMERICAL | two independent integrals agree | ~65 | `step_29.py` |
| truncation: 1.02 → 35 digits ($N{=}56$); 1.05 optimal → 17 digits ($N^*{=}32$) | NUMERICAL | table above | 35 / 17 | `step_29.py`,`step_30.py` |
| Padé–Borel at $p=1.05$ | not certified | `mp.pade` evaluation bug; results discarded | — | `step_30.py` |
| single-appearance lemma; $u_k(2)$, $k\le12$ | PROVED + exact check | $4u_k(2)=\kappa_k$ exactly, $k=3..11$ | exact | `step_21.py` + D11c values |
| invariance $[Z_n][d^n]=[Z_n][\tau^n]$ | PROVED | one-line weight argument | exact | text |
| $u_2,u_3,b_2,b_3,b_4$ closed forms in $n$ | PROVED | match H1's general-$n$ data $n\le12$ (diff 0 there) | exact | text + H1 data |
| zeta-free truncations (all-zero ⇒ $\mathrm{Rat}_n$; only-$\zeta(2)$ ⇒ Gaussian single sum) | PROVED algebraically; numeric checks not run | derivation in §4 | — | text |
| factorial growth of rational coefficients | NUMERICAL (+D11/D11c proved asymptotics) | growth tables | — | `step_29/30.py` |
| $b_k(n)$, $k\ge5$, closed form; rhombus symbolic $q$-series rerun | open | not completed in time | — | — |

**Files in the working directory:** `d16_engine.py` (exact engine, written by `step_20.py`), `tau_series_2d.json` (exact $P/\sqrt2$, $\log P$, $B_2$, formal $Z$'s, $\tau\le12$), `tau_reduced_2d.json` (reduced series + exact $a_{11},a_{12}$), `hi_tau_numeric_n2.json` (64 coefficients, 150 digits), `quad_ref.json` (quadrature references), `B2_exact.json`; scripts `step_16.py`–`step_30.py` (engine build `step_18/20`, exact series `step_21`, operator-formula route `step_23`, $d$-checks `step_24/26`, reduced series `step_25`, numeric engine `step_28`, quadrature+truncation `step_29/30`).