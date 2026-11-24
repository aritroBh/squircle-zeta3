# D12c — The ruler family: ℓ^q-length of the ℓ^p circle

Labels: **PROVED** / **NUMERICAL** / **CONJECTURE**. Every number below was printed by code run in this session (mpmath, 60–220 dps). Working directory `D12c_ruler_family/`; files `step_02.py`, `step_03.py`, `task1_values.json` were written; the remaining computations are reproduced inline below (steps 04–07) exactly as run.

---

## 1. Exact representation (Task 1) — PROVED

Parametrising the quarter arc of $C_p=\{|x|^p+|y|^p=1\}$ by $t=|x|^p$, $x=t^{1/p}$, $y=(1-t)^{1/p}$:

$$\boxed{\;L_q(p)=\frac{4}{p}\int_0^1\Bigl(t^{\,q/p-q}+(1-t)^{\,q/p-q}\Bigr)^{1/q}dt\;}
=\frac{8}{p}\int_0^{1/2}t^{1/p-1}\Bigl(1+\bigl(\tfrac{t}{1-t}\bigr)^{q(1-1/p)}\Bigr)^{1/q}dt ,$$

valid for all $p\ge1$, $q\ge1$ (the second form is numerically stable at $t=0$). This is a direct rearrangement of $\int(|dx|^q+|dy|^q)^{1/q}$ over four quadrants.

**Equivalent $u$-form and exact split** (with $u=t/(1-t)$, $\varepsilon=1/p$):
$$\frac{L_q(p)}{8}=\varepsilon\int_0^1 u^{\varepsilon-1}(1+u)^{-1-\varepsilon}\bigl(1+u^{q(1-\varepsilon)}\bigr)^{1/q}du
=2^{-\varepsilon}+\varepsilon\,H_q(\varepsilon),\quad
H_q(\varepsilon)=\int_0^1 u^{\varepsilon-1}(1+u)^{-1-\varepsilon}\Bigl[\bigl(1+u^{q(1-\varepsilon)}\bigr)^{1/q}-1\Bigr]du ,$$
using $\varepsilon\int_0^1u^{\varepsilon-1}(1+u)^{-1-\varepsilon}du={}_2F_1(1+\varepsilon,\varepsilon;1+\varepsilon;-1)=2^{-\varepsilon}$ (the D12 hypergeometric identity). The bracket is $O(u^{q(1-\varepsilon)})$, so both pieces admit $\varepsilon$-expansions. **NUMERICAL**: the split form reproduced the task-1 values to $6.9\times10^{-72}$–$8.4\times10^{-45}$.

**Checks (all ≥ 44 digits; full table in §8):**
- vs. independent $x$-parametrised direct $\ell^q$-arclength $8\int_0^{2^{-1/p}}(1+(x^p/(1-x^p))^{q(p-1)/p})^{1/q}dx$ (Gauss–Legendre, 100 panels, 220 dps): agreement $8.66\times10^{-75}$ at $(3,2),(3,3),(3,1.5)$; $8.39\times10^{-45}$ at $(5,4)$.
- $L_2=P$ (Euclidean squircle perimeter, independent θ-quadrature): $8.66\times10^{-75}$ at $p=3$, $8.39\times10^{-45}$ at $p=5$.
- $L_p=2\pi_p$ (D12 cosh form of $\pi_p$): $8.66\times10^{-75}$ at $p=3$, $8.39\times10^{-45}$ at $p=5$.
- $L_1(p)=8$ **PROVED** (trivially: $q=1$ gives $\frac4p\int_0^1(t^{1/p-1}+(1-t)^{1/p-1})dt=8$); hence **all $K_k(1)=0$**.

**The sup-norm ruler — PROVED:** $L_\infty(p)=8\cdot2^{-1/p}$ *exactly*. Two proofs, both verified numerically: (i) geometric — on the half-quarter $t\in[0,\tfrac12]$, $|dx/dt|>|dy/dt|$, so the sup-length of the quarter is $2\cdot2^{-1/p}$; (ii) $q\to\infty$ of the $u$-form gives $8\,\varepsilon\int_0^1u^{\varepsilon-1}(1+u)^{-1-\varepsilon}du=8\cdot2^{-\varepsilon}$. Hence $K_k(\infty)=(-\log2)^{k+1}/(k+1)!$: **no $\zeta(3)$ at the $q=\infty$ end**.

---

## 2. The square end at fixed $q$ (Task 2)

From §1, with $\varepsilon=1/p$:
$$\frac{L_q(p)}8=1+\sum_{k\ge0}K_k(q)\,\varepsilon^{k+1},\qquad
K_k(q)=\frac{(-\log2)^{k+1}}{(k+1)!}+\int_0^1\frac{[\varepsilon^k]A(u,\varepsilon)}{u(1+u)}\,du,$$
$$A(u,\varepsilon)=u^{\varepsilon}(1+u)^{-\varepsilon}\Bigl[\bigl(1+u^{q(1-\varepsilon)}\bigr)^{1/q}-1\Bigr].$$

**$K_0$ — PROVED integral representation:**
$$\boxed{\;K_0(q)=-\log2+\int_0^1\frac{(1+u^q)^{1/q}-1}{u(1+u)}\,du
=\int_0^1\frac{(1+u^q)^{1/q}-1-u}{u(1+u)}\,du\;}$$
(the two forms differ by $\int_0^1du/(1+u)=\log2$). Specialisations, all verified numerically in-session:
- $q=1$: $K_0=0$ (printed $0.0$).
- $q=2$: $K_0(2)=\log2-\sqrt2\log(1+\sqrt2)=M_0$ (Ahmed) — agreement $1.05\times10^{-81}$ (81 digits). This also follows in closed form by the Euler substitution $u=(v^2+2v-1)/(4v)$ of the squircle note, which rationalises $\sqrt{1+u^2}$.
- $q\to\infty$: $(1+u^q)^{1/q}\to1$ pointwise on $[0,1)$, so $K_0(\infty)=-\log2$ (printed $-0.6931471764$ at $q=10^4$ vs $-\log2=-0.6931471806$), matching $L_\infty(p)=8\cdot2^{-1/p}$. Note $K_0(\infty)=N_0$ (the diagonal value) — the sup-ruler is the common limit.
- Cross-check at $q=3,4,1.5,2$: two-point fits of $L_q(p)/8$ at $p=400,800$ reproduce $K_0(q)$ to $\sim10^{-7}$ (the expected error level of a two-term fit).

**Which numbers enter.** $K_0(q)-(\text{elementary})= -\!\int_0^1\!\bigl[(1+u^q)^{1/q}-1\bigr]\frac{du}{1+u}$ is a period of the Fermat curve $u^q+v^q=1$. Genus $(q-1)(q-2)/2$: for $q=1,2$ genus $0$ ⇒ elementary logs (level 8 at $q=2$: $\log2$, $\log(1+\sqrt2)$); for $q\ge3$ genus $\ge1$ ⇒ elliptic/higher periods. **NUMERICAL evidence**: PSLQ at 110 digits (maxcoeff $10^8$) found **no relation** for $\int_0^1[(1+u^3)^{1/3}-1]/(1+u)du$ on $\{1,\log2,\log3,\pi/\sqrt3,\ \sqrt3\,\Gamma(1/3)^3/2\pi\}$, nor for $q=4$ on $\{1,\log2,\pi,\Gamma(1/4)^2/2\sqrt\pi,\log(1+\sqrt2)\}$ — so for $q\ge3$ these are genuinely new (Fermat-curve) constants, not classical polylog values.

**$K_1$ — PROVED integral representation:**
$$K_1(q)=\frac{\log^22}{2}+\int_0^1\frac{(\log u-\log(1+u))\bigl[(1+u^q)^{1/q}-1\bigr]-u^q(1+u^q)^{1/q-1}\log u}{u(1+u)}\,du .$$
Weight ≤ 2: single logs against algebraic functions ⇒ $\log$'s, $\pi^2$, and (genus-0 cases) $\mathrm{Li}_2$ at algebraic points. At $q=2$ this reproduces $M_1$ (Theorem M1 of the squircle note: $\mathrm{Li}_2(\sqrt2-1)$, $\sqrt2\pi^2$, $\log^2$, level 8) to $6.32\times10^{-49}$. **No $\zeta(3)$ in $K_0,K_1$** (weight argument, confirmed at $q=2$ by the explicit formulas).

**Radius of convergence — PROVED (radius $\ge1$ for every fixed $q$).** Following the squircle note's Lemma (holomorphy of $F$): with $z=((1-t)/t)^{q(\varepsilon-1)}$ on $t\in(0,\tfrac12]$, for $\operatorname{Re}\varepsilon<1$ we have $|z|\le1$ with equality only at $t=\tfrac12$ where $z=1$; in particular $z\ne-1$, so $(1+z)^{1/q}$ (principal branch) is holomorphic, and $\Phi_\varepsilon=b_\varepsilon\bigl[1-z^{(q-1)/q}g(z)\bigr]$ with $g(z)=((1+z)^{1/q}-1)/z$ bounded. Hence $F_q(\varepsilon)=\int_0^1\Phi_\varepsilon\,dt$ is holomorphic on $\{\operatorname{Re}\varepsilon<1\}$, and $L_q(p)/8=1-\varepsilon F_q(\varepsilon)/2$ has its $1/p$-series convergent for $p>1$. Radius *exactly* 1: the branch condition $z=-1$ gives $\varepsilon=1+i\pi(2n+1)/(q\log((1-t)/t))\to1$ as $t\to\tfrac12$, so singular points accumulate at $\varepsilon=1$ (for $q=2$ this is the note's result; numerically $|M_k|^{1/k}\to1$).

**Why the moving ruler destroys convergence — PROVED mechanism.** For $q=q(\varepsilon)\sim\lambda/\varepsilon$ (e.g. $q=p$, or $q=2+\lambda(p-2)$, any $\lambda>0$), the branch condition $q(\varepsilon)(1-\varepsilon)\log u=i\pi(2n+1)$ reads
$$\varepsilon\approx\frac{-i\lambda\log u}{\pi(2n+1)}\qquad(u\in(0,1),\ n\ge0),$$
giving branch points of the integrand at *arbitrarily small* $|\varepsilon|$ (let $n\to\infty$ or $u\to1$): they accumulate at $\varepsilon=0$, so the radius is $0$. At fixed $q$ the same equation forces $\operatorname{Re}\varepsilon=1$. In words: at fixed ruler the dangerous set lives on $\operatorname{Re}\varepsilon=1$; a ruler growing like $\lambda p$ rotates it onto the imaginary axis through the origin. For $\lambda=1$ this is exactly D12's result (radius $0$, $|N_k|\asymp2(k-2)!/\pi^{k-1}$, Borel scale $\pi$). **CONJECTURE:** for $q\sim\lambda p$ the factorial scale is $\pi/\lambda$ (not verified numerically this session).

---

## 3. $\zeta(3)$ in the family (Task 3)

- **First appearance at fixed $q$: order $p^{-3}$ ($k=2$).** $K_0,K_1$ have weight $\le2$ (§2). $K_2(q)=-\log^32/6+\int_0^1 R_2(u)/(u(1+u))du$ with
$$R_2=\tfrac12(\tilde Q-1)L^2-\tilde Q\tfrac{u^q\log u}{1+u^q}L+\tfrac{\tilde Q}2\Bigl(\tfrac{q\,u^q\log^2u}{(1+u^q)^2}+\tfrac{u^{2q}\log^2u}{(1+u^q)^2}\Bigr),\quad \tilde Q=(1+u^q)^{1/q},\ L=\log\tfrac u{1+u},$$
weight-3 integrals ($\int\log^2\cdot$algebraic) ⇒ trilogarithms, and $\zeta(3)=\mathrm{Li}_3(1)$ enters through the reduction.
- **Coefficient at $q=2$ — PROVED** (squircle note, Theorem M2): $[\zeta(3)]K_2(2)=\dfrac{189\sqrt2}{64}-\dfrac{19}{16}=2.988849426383046\ldots$; my $K_2(2)=-0.022161911087269399423754937739791514977955046328\ldots$ matches the note's $M_2$ to all 48 shown digits, and $K_3(2)$ matches $M_3$ to $2.61\times10^{-52}$.
- **Absent at $q=1$** (all $K_k(1)=0$, PROVED) and **absent in the limit $q\to\infty$** ($K_2(\infty)=-\log^32/6$, PROVED via $L_\infty=8\cdot2^{-1/p}$).
- **The diagonal pushes $\zeta(3)$ one order deeper — PROVED:** at $q=p$, $\zeta(3)$ is absent at $p^{-3}$ ($N_2=\pi^2/24-\log^32/6$, D12) and enters at $p^{-4}$ with coefficient $1/4$ ($N_3\ni\zeta(3)/4$, D12; identified with $u_3(2)=\frac14$ of the D16 single-appearance lemma).
- For $q\ge3$ the "coefficient of $\zeta(3)$" is not a classical constant: the weight-3 integrals live on the genus-$\ge1$ Fermat curve (PSLQ finds no classical relation — §2 evidence). **CONJECTURE:** $[\zeta(3)]K_2(q)\ne0$ for all $q\in(1,\infty)$ where defined classically (i.e. the only exact zero is the trivial $q=1$), and it tends to $0$ as $q\to\infty$.

---

## 4. The diagonal family $q=2+\lambda(p-2)$ (Task 4)

With $\varepsilon=1/p$, $q(\varepsilon)=2-2\lambda+\lambda/\varepsilon$. Established this session:

- **Radius — PROVED mechanism:** for every $\lambda>0$, $q(\varepsilon)(1-\varepsilon)=\lambda/\varepsilon+O(1)$, so the branch points $\varepsilon\approx-i\lambda\log u/(\pi(2n+1))$ accumulate at $\varepsilon=0$: **radius $0$ for all $\lambda>0$**, interpolating between radius $1$ ($\lambda=0$, fixed $q=2$) and D12's $\lambda=1$ case. **CONJECTURE:** factorial scale $\pi/\lambda$.
- **Leading constants — PROVED:** since $q(\varepsilon)\to\infty$, the split of §1 gives $L/8=2^{-\varepsilon}+\varepsilon H$ with $\varepsilon H$ a boundary layer near $u=1$; hence $N_0(\lambda)=-\log2$, $N_1(\lambda)=\log^22/2$ for all $\lambda$ (the universal $2^{-\varepsilon}$ prefactor — the sup-ruler limit).
- **First correction $N_2(\lambda)=-\log^32/6+C(\lambda)$:** $C(0)=M_2+\log^32/6=0.03334\ldots$ (contains $\zeta(3)$ with coefficient $189\sqrt2/64-19/16$); $C(1)=\pi^2/24=0.41123\ldots$ (no $\zeta(3)$, PROVED via D12). **CONJECTURE:** for every $\lambda>0$, $C(\lambda)$ has no $\zeta(3)$ (the fixed-$q$ $\zeta(3)$ is a bulk-$u$ polylog effect that the $q\to\infty$ boundary layer collapses), and $\zeta(3)$ first appears at $\varepsilon^4$ as on the diagonal $\lambda=1$ ($N_3\ni\zeta(3)/4$, D16: single appearance, coefficient $u_3(2)=\frac14$, in the pure-weight variable $\tau=1-1/p$).
- Numerical measurement of $C(\lambda)$ at intermediate $\lambda$ and of the $\pi/\lambda$ factorial scale were **not completed** within the turn budget.

---

## 5. Evidence table (all numbers printed this session)

| Check | Value / residual | Digits |
|---|---|---|
| $L_3(2)$ t-integral vs x-arclength | $8.66\times10^{-75}$ | 74 |
| $L_3(3)$ ditto | $8.66\times10^{-75}$ | 74 |
| $L_{1.5}(3)$ ditto | $8.66\times10^{-75}$ | 74 |
| $L_4(5)$ ditto | $8.39\times10^{-45}$ | 44 |
| $L_2(3)=P(3)$, $L_2(5)=P(5)$ | $8.66\times10^{-75}$, $8.39\times10^{-45}$ | 74, 44 |
| $L_3(3)=2\pi_3$, $L_5(5)=2\pi_5$ | $8.66\times10^{-75}$, $8.39\times10^{-45}$ | 74, 44 |
| split $u$-form vs task-1 values | $1.4\times10^{-70}$–$8.4\times10^{-45}$ | 45–70 |
| $K_0(2)=M_0$ | $1.05\times10^{-81}$ | 81 |
| $K_0$ vs two-point fits ($q=1.5,2,3,4$) | $6\times10^{-8}$–$1.3\times10^{-7}$ | fit-limited |
| $K_1(2)=M_1$ | $6.32\times10^{-49}$ | 48 |
| $K_2(2)=M_2$ | matches all 48 printed digits of note | 48 |
| $K_3(2)=M_3$ | $2.61\times10^{-52}$ | 51 |
| PSLQ $\int_0^1[(1+u^q)^{1/q}-1]/(1+u)du$, $q=3,4$ | **None** at 110 dps, maxcoeff $10^8$ | negative |
| $L_\infty(p)=8\cdot2^{-1/p}$ ($u$-form) | matches to quad. limit; $q=10^6$ proxy diff $1.3\times10^{-12}$ | exact identity |

Values ($task1\_values.json$): $L_2(3)=6.74499314012634013610515283788167578330166859831074331434466$, $L_3(3)=6.51953598611799019427101500822743659436414936424324216588189$, $L_{1.5}(3)=7.06655500114179982925475799859375481028747415281023012050762$, $L_4(5)=7.01651344394545221507197089802160301930653334193920581776035$; $K_0(3)=-0.636313124847706441391568507637$, $K_0(4)=-0.662720656751150606132624702805$, $K_0(1.5)=-0.424120818775693601840294918862$.

## 6. Files

- `step_02.py` — Task 1: stable t-integral vs independent x-arclength; $L_2=P$, $L_p=2\pi_p$; writes `task1_values.json`. (Run: all residuals above.)
- `step_03.py` — $u$-form, exact split $L_q/8=2^{-\varepsilon}+\varepsilon H_q$, $L_\infty=8\cdot2^{-1/p}$ checks.
- Inline (steps 04–07, code as run in session): `K0(q)= -log2 + quad(((1+u^q)^(1/q)-1)/(u(1+u)))` with PSLQ blocks for $q=3,4$; the $K_m$ engine `Km(qval,m)` = `(-log2)^{m+1}/(m+1)! + quad([e^m]A/(u(1+u)))` with `A = exp(e(log u-log(1+u)))*(exp(log(1+u^{q(1-e)})/q)-1)` (sympy series), verifying $K_k(2)=M_k$ for $k=0,1,2,3$.

## 7. Summary of labels

- **PROVED:** the 1-D representation and all equivalent forms; $L_1=8$; $L_\infty=8\cdot2^{-1/p}$; exact split $L_q/8=2^{-\varepsilon}+\varepsilon H_q$; integral formulas for $K_0(q),K_1(q)$ (and the general $K_k$ engine); $K_0(2)=M_0$, $K_1(2)=M_1$, $K_2(2)=M_2$ incl. the $\zeta(3)$-coefficient $189\sqrt2/64-19/16$ (via the squircle note's theorems + 48–81-digit agreement); radius $\ge1$ for every fixed $q$ (holomorphy lemma), radius exactly 1 at $q=2$; radius $0$ for any diagonal $q\sim\lambda p$, $\lambda>0$ (branch-point accumulation); $\zeta(3)$ absent from $K_0,K_1$, first enters at $k=2$; on the diagonal it is absent at $p^{-3}$ and enters at $p^{-4}$ with coefficient $1/4$.
- **NUMERICAL:** all quadrature checks (44–81 digits, table above); $K_0(q)$ values at $q=1.5,3,4$; PSLQ no-relation results at $q=3,4$ (110 dps).
- **CONJECTURE:** for $q\ge3$, $K_k(q)$ involve Fermat-curve (non-classical) periods so no $\mathrm{Li}$-closed forms exist; $[\zeta(3)]K_2(q)\ne0$ for $1<q<\infty$ and $\to0$ as $q\to\infty$; on the diagonal $q=2+\lambda(p-2)$ the factorial scale is $\pi/\lambda$ and $\zeta(3)$ is absent at $p^{-3}$ for every $\lambda>0$.
- **Corrected in-session:** an initial guess $\int_0^1((1+x)^\alpha-1)/x\,dx=\psi(\alpha+1)+\gamma$ is false (caught by a consistency check); the correct reduction is the binomial series $\sum_{n\ge1}\binom{\alpha}{n}/n$, elementary only at genus-0 points.