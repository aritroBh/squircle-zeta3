# FLD1: the Field rule (Conjecture conj:field) for every $k$ — PROOF

**Theorem (Field rule).** For every integer $k\ge 0$ and every $0<\alpha<\pi$, with $c=\cos\frac\alpha2$,
$$M_k(\alpha)=\sum_{j} z_j\, r_j(c)\, H_{w_j}(c),$$
a finite sum in which
* $z_j$ are monomials in $\zeta(2),\dots,\zeta(k)$;
* $r_j\in\mathbb Q(c)$, with denominators dividing a product $c^a(1-c)^b(1+c)^d$;
* the $H_{w_j}$ are harmonic polylogarithms with letters in $\{0,1,-1\}$;
* the weights are bounded by $\operatorname{weight}(z_j)+|w_j|\le k+1$.

Moreover, for every $k\le 7$ the denominators divide a power of $c$ alone.

**Consequence (the Conjecture).** $H_w(c)$ with letters $\{0,\pm1\}$ are multiple polylogarithms whose arguments lie in
the multiplicative group $\langle 2,\,c,\,1-c,\,1+c\rangle\subset\mathbb Q(c)^{\times}$ — i.e. in the $S$-units with
$S=\{\text{prime ideals of }\mathbb Q(c)\text{ dividing }2,\ c,\ 1-c,\ 1+c\}$. Hence every $M_k(\alpha)$ is a combination of
polylogarithms at $S$-units of $\mathbb Q(\cos\frac\alpha2)$: this is Conjecture **conj:field**, with this explicit $S$.

*Remark.* The power-of-$c$ denominators and the sharp pole set $\{\tfrac12,0\}\cup\{-\tfrac t2:1\le t\le k-1\}$
are exact computations for $k\le 7$ and conjectured for all $k$; they would follow in general from the properness of
the coefficients $\rho_w(n)$ in Lemma A, which for $k\ge8$ is open (it needs a $Z$-sum independence lemma; see Lemma A).

The proof is constructive: it gives an algorithm that produces the HPL form, and proves it correct, for every $k$.
For $k\le 7$ every step is additionally verified by exact computation (Section "Verification"); an independent referee
(FLD1R) re-verified Lemmas A, B, C, E, F for all $k$ and repaired Lemmas A2, C and D as incorporated below.

---

## Setup

$c=\cos\frac\alpha2$. By **Theorem thm:flat** (proved in the paper), for every $k\ge0$,
$$M_k(\alpha)=\sum_{n\ge1}\mu_{k,n}\,c^{2n},\qquad |c|<1,$$
with $\mu_{k,n}\in\mathbb Q[\zeta(2),\dots,\zeta(k)]$. The series converges for $|c|<1$ and defines an even function of $c$.
We prove the sum is an HPL combination.

**Harmonic polylogarithms (Remiddi–Vermaseren).** $H_1=-\log(1-c)$, $H_{-1}=\log(1+c)$, $H_0=\log c$, and
$H_{a,w}(c)=\int_0^c f_a(t)H_w(t)\,dt$ with $f_0=1/t$, $f_1=1/(1-t)$, $f_{-1}=1/(1+t)$. Only words $w$ over $\{0,1,-1\}$
with nonzero last letter occur (these are analytic at $0$, $H_w(0)=0$). The **weight** $|w|$ is the word length.

**Alternating nested sums (Z-sums).** $Z_{a_1,\dots,a_d}(M)=\sum_{M\ge i_1>\cdots>i_d\ge1}\prod_j \operatorname{sgn}(a_j)^{i_j}/i_j^{|a_j|}$,
$a_j\in\mathbb Z\setminus\{0\}$; $Z_{()}=1$; $|w|=\sum|a_j|$. They satisfy the **stuffle** (quasi-shuffle) product
$Z_uZ_v=\sum Z_w$ with $|w|=|u|+|v|$, and $Z_w(M)=Z_w(M-1)+\operatorname{sgn}(a_1)^M M^{-|a_1|}Z_{w'}(M-1)$ (shift).

---

## Lemma A (exact nested-sum form of $\mu_{k,n}$)

For every $k\ge0$ and $n\ge1$,
$$\mu_{k,n}=-\frac1{n(2n-1)}\sum_{(j,i,\ell)} e^{(k)}_{j,i,\ell}(n)\,\frac{(n)_j}{(2n)_j}\,D_{i,\ell}(n+j,n),$$
where $[u^k]E_n=\sum e^{(k)}_{j,i,\ell}(n)\,y^j a^i b^\ell$ (with $E_n=e^{uA}(1+y(e^{u\ell}-1))^{1-2n}$, $A=-(n-1)a+nb$, $\ell=b-a$),
$e^{(k)}_{j,i,\ell}\in\mathbb Q[n]$ nonzero only for $0\le j\le k$ and $i+\ell\le k$, and
$D_{i,\ell}(p,q)=\big(\partial_p^{\,i}\partial_q^{\,\ell}B(p,q)\big)/B(p,q)$.

Consequently, separating by $\zeta$-monomial $z$ (a monomial in $\zeta(2),\dots,\zeta(k)$),
$$[\mu_{k,n}]_z=\sum_{w}\rho^{(z)}_w(n)\,Z_w(2n),$$
a finite sum, in which
* **(A1)** each word $w$ over $\mathbb Z\setminus\{0\}$ has $|w|\le k-\operatorname{weight}(z)$;
* **(A2)** each $\rho^{(z)}_w\in\mathbb Q(n)$ has poles only in $\{\tfrac12,0\}\cup\{-\tfrac t2:1\le t\le 2k-2\}\subset\tfrac12\mathbb Z$
  (proved below in this coarser form; the sharp set with $t\le k-1$ is exact for $k\le7$), and is **proper**:
  $\deg\text{num}\le\deg\text{den}-2$ (proved computationally for $k\le7$; open for $k\ge8$, where it needs a
  $Z$-sum independence lemma).

**Proof.**
(i) *The integral formula.* With $y=g/(1+g)$ and the $g\to1/g$ symmetry, the corner coefficient becomes the Beta
integral above: $\int_0^1 y^{n+j-1}(1-y)^{n-1}a^i b^\ell\,dy=\partial_p^i\partial_q^\ell B(p,q)\big|_{(n+j,n)}
=B(n+j,n)D_{i,\ell}(n+j,n)$, and $B(n+j,n)/B(n,n)=(n)_j/(2n)_j$. This is the derivation of `flat_exact.py` /
`mu3_symbolic.py`, which we re-implement independently for general $k$ (`step_21`, `step_22`).
(ii) *The $u$-expansion.* $(1+y(e^{u\ell}-1))^{1-2n}=\sum_m\binom{1-2n}{m}y^m(e^{u\ell}-1)^m$ and
$(e^{u\ell}-1)^m=m!\sum_{j\ge m}S(j,m)u^j\ell^j/j!$ (Stirling numbers of the second kind), so
$[u^k]E_n=\sum_{i+j=k}\frac{A^i}{i!}\big(\sum_m\binom{1-2n}{m}y^m m! S(j,m)\big)\frac{\ell^j}{j!}$.
Since $\binom{1-2n}{m}\in\mathbb Q[n]$ (degree $m\le j\le k$) and $A^i\ell^j$ has $a,b$-degree $i+j=k$ with
$\mathbb Q[n]$-coefficients of degree $\le k$, the bounds $j\le k$, $i+\ell\le k$, $e\in\mathbb Q[n]$ follow.
(Free sharpening, noted by the referee and checked at $k=2,3,5$: the support has $i+\ell=k$ *exactly*, since every
summand carries $A^i\ell^j$ of $a,b$-degree $i+j=k$.)
(iii) *Cumulants and polygammas.* $D_{i,\ell}$ is the Bell polynomial in the cumulants
$\kappa_{r,0}=\psi^{(r-1)}(p)-\psi^{(r-1)}(p+q)$, $\kappa_{0,m}=\psi^{(m-1)}(q)-\psi^{(m-1)}(p+q)$,
$\kappa_{r,m}=-\psi^{(r+m-1)}(p+q)$ ($r,m\ge1$). At integers,
$\psi^{(m)}(N)=(-1)^{m+1}m!\big(\zeta(m+1)-H^{(m+1)}_{N-1}\big)$ ($m\ge1$) and $\psi(N)+\gamma=H_{N-1}$.
Hence the **pure** cumulants ($r=0$ or $m=0$) are $\zeta$-free (the $\zeta$ cancels in the difference), and the
**mixed** ones carry the single zeta $\zeta(r+m)$. So $D_{i,\ell}(n+j,n)$ is a polynomial in
$H^{(r)}_{n+j-1},H^{(r)}_{n-1},H^{(r)}_{2n+j-1}$ ($r\le i+\ell$) and $\zeta(2),\dots,\zeta(i+\ell)$, with
$\zeta$-weight $+$ harmonic-weight $\le i+\ell\le k$.
(iv) *Reduction to two harmonic points.* $H^{(r)}_{n+j-1}=H^{(r)}_{n-1}+\sum_{t<j}(n+t)^{-r}$ and
$H^{(r)}_{2n+j-1}=H^{(r)}_{2n-1}+\sum_{t<j}(2n+t)^{-r}$, so $\mu_{k,n}\in\mathbb Q(n)[H^{(r)}_{n-1},H^{(r)}_{2n-1},\zeta(2..k)]$,
$\zeta$-weight $+$ harmonic-weight $\le k$.
(v) *Harmonic $\to$ Z-sum at $2n$.* $H^{(r)}_{n-1}=2^{r-1}\big(Z_r(2n)+Z_{-r}(2n)\big)-n^{-r}$ and
$H^{(r)}_{2n-1}=Z_r(2n)-(2n)^{-r}$ (immediate from the definitions). Products are handled by stuffle, which preserves
$|w|$. Thus per $\zeta$-monomial $z$: $[\mu_{k,n}]_z=\sum_w\rho_w(n)Z_w(2n)$ with $|w|\le k-\operatorname{weight}(z)$
(the zeta $z$ uses $\operatorname{weight}(z)$ of the budget; $H^{(r)}\mapsto Z_{\pm r}$ preserves weight). This is (A1).
(vi) *The rational coefficients.* They are assembled from the prefactor $1/(n(2n-1))$ (poles at $0,\frac12$), the
proper ratio $(n)_j/(2n)_j$ (poles at $n=-t/2$, $t\ge1$; there is *no* pole at $n=0$, since the $t=0$ factor is
$n/(2n)=\tfrac12$), the shift corrections $(n+t)^{-r}$ (poles at $n=-t$, $1\le t\le j-1\le k-1$) and
$(2n+t)^{-r}$ (poles at $n=-t/2$), and $e(n)\in\mathbb Q[n]$. All poles therefore lie in
$\{\frac12,0\}\cup\{-\frac t2:1\le t\le 2k-2\}$: the $(n+t)^{-r}$ contribute the integer points $-t=-\frac{2t}{2}$
with $2t\le 2k-2$, which is why the proved set is the coarser one. (The sharp set $\{\frac12,0\}\cup\{-\frac t2:1\le t\le k-1\}$
would require arguing that the integer-pole terms cancel against the Pochhammer numerators; this holds for all $k\le7$
but is not argued in general. Either version suffices downstream, since all poles are in $\frac12\mathbb Z$, hence
integral in $m=2n$.) This proves the pole part of (A2).
*Properness* ($\deg\text{num}-\deg\text{den}=-2$ for every $\rho_w$): verified *exactly* for $k\le7$ — for $k=7$ all
553 word-terms of the $Z$-sum canonical form have $\deg=-2$ (`FLD1R step_32/34`), sharpening our `step_24` ($k\le5$) and
`step_31` ($k\le6$). The heuristic reason is the Beta$(n,n)$-concentration $\mu_{k,n}=O(n^{-2})$; turning this into a
*per-word* proof for all $k$ needs a linear-independence/separation lemma for the family $\{Z_w(2n)\}$ over
$\mathbb Q(n)$ (available in principle via Schneider's difference-ring theory, but not cited or proved), which we leave
**open for $k\ge8$**. Were a polynomial part ever to occur it is handled by $c\frac{d}{dc}$ in Lemma C and produces
only the still-admissible denominators $1\pm c$. $\blacksquare$

---

## Lemma B (the dictionary: HPL word $\leftrightarrow$ nested sum)

**(B1)** For every HPL word $w$ over $\{0,1,-1\}$ with nonzero last letter,
$$[c^N]\,H_w(c)=\operatorname{cf}_w\,\varepsilon_w^{N}\,N^{-k_w}\,Z_{v_w}(N-1),$$
where, writing the nonzero letters of $w$ as $b_1,\dots,b_s$ (left to right, $b_s$ = last letter), with $m_0$ leading
zeros and $m_t$ zeros between $b_t$ and $b_{t+1}$:
$$\varepsilon_w=b_1,\qquad k_w=1+m_0,\qquad \operatorname{cf}_w=b_1\cdots b_s\in\{\pm1\},\qquad
v_w=\big(b_1b_2(1+m_1),\dots,b_{s-1}b_s(1+m_{s-1})\big).$$
The map $w\mapsto(\varepsilon_w,k_w,v_w)$ is a **bijection**
$$\{w\}\;\longleftrightarrow\;\{\pm1\}\times\mathbb Z_{\ge1}\times\{\text{words over }\mathbb Z\setminus\{0\}\}.$$

**(B2)** Under this bijection $|w|=k_w+|v_w|$: the HPL weight equals the nested-sum weight $k+|v|$.

**Proof.** By induction on the word length using $h_{a,w}(N)=\frac1N\sum_{i<N}a^{N-1-i}h_w(i)$ for $a=\pm1$ and
$h_{0,w}(N)=h_w(N)/N$ for a prepended $0$, with the base $h_{(a)}(N)=a^{N+1}/N$. Unfolding gives the displayed
$(\operatorname{cf}_w,\varepsilon_w,k_w,v_w)$. Injectivity: the triple determines $s=|v|+1$, then
$b_1=\varepsilon$, $b_{t+1}=b_t\operatorname{sgn}(v_t)$, $m_0=k-1$, $m_t=|v_t|-1$, which rebuild $w$. Surjectivity: any
$(\varepsilon,k,v)$ has $m_0=k-1\ge0$ and $m_t=|v_t|-1\ge0$, so the recipe yields a valid word. (B2):
$|w|=s+\sum_{t=0}^{s-1}m_t=s+(k-1)+\sum_{t=1}^{s-1}(|v_t|-1)=k+|v|$. Verified against `zsum.hpl_coeff_terms` for all
242 words of length $\le5$ (`step_20`). $\blacksquare$

---

## Lemma C (the construction: even nested sums $\to$ HPL combination)

Let $E(n)=\sum_w\rho_w(n)Z_w(2n)$ be finite with $\rho_w\in\mathbb Q(n)$ proper, poles at half-integers. Define the even
sequence $a(2n)=E(n)$, $a(2n+1)=0$. Then there is an HPL combination $F(c)=\sum_j r_j(c)H_{w_j}(c)$ with
$r_j\in\mathbb Q(c)$ (denominators powers of $c$) and $[c^N]F=a(N)$ for all $N$.

**Proof (the algorithm; each step is exact rational arithmetic).**
*Step 1 (match the leading sums).* Order the words by decreasing weight. For a leading word $w$ with coefficient
$\rho_w(n)$, partial-fraction $\rho_w$ in $m=2n$ (poles are integers; write the proper part as
$\sum_{e,\kappa}c_{e,\kappa}(2n-e)^{-\kappa}$ plus a polynomial part $P$). For each $(e,\kappa)$ take the HPL word
$w_0=\mathrm{bji}(1,\kappa,w)$ given by the Lemma B bijection inverse, and the term
$\frac{c_{e,\kappa}}{\operatorname{cf}_{w_0}}\,c^eH_{w_0}(c)$. By (B1),
$[c^{2n}]\,c^eH_{w_0}=\operatorname{cf}_{w_0}(2n-e)^{-\kappa}Z_w(2n)+\text{(strict-suffix, lower-weight) corrections}$,
so its leading sum is exactly $c_{e,\kappa}(2n-e)^{-\kappa}Z_w(2n)$. A polynomial part $P(n)=\sum_j p_j n^j$ is realized
by applying $(c\frac{d}{dc})^j$ to the corresponding $w_0$-terms, since $[c^{N}](c\frac d{dc})^jF=N^j[c^N]F$ and
$c\frac{d}{dc}$ preserves the class $\{r(c)H_w:\mathrm{den}\mid c(1-c)(1+c)\}$
($c f_0=1$, $c f_1=c/(1-c)$, $c f_{-1}=c/(1+c)$) and the HPL weight. Summing over $(e,\kappa)$ and the $c\frac d{dc}$
terms, the $Z_w(2n)$ terms match $\rho_w(n)Z_w(2n)$ exactly; subtracting, the remainder has strictly smaller maximal
weight (the corrections are strict suffixes). This terminates at the empty word $w=()$, matched by single-letter HPLs.
(The polynomial-part/$c\frac d{dc}$ branch is sound but is never triggered in practice — it is vacuous for $k\le7$,
where every $\rho_w$ is proper.)
The result $F_0$ has $[c^{2n}]F_0=E(n)$ for all $n\ge1$. Here is the required regularity argument (repairing the
earlier, incorrect parenthetical that "no denominator $2n+t$ vanishes at $n\ge1$" — e.g. $2n-2$ *does* vanish at $n=1$):
both sides are, by construction, $\mathbb Q(n)$-linear combinations of the $Z_w(2n)$, and the identity between them is an
identity of *rational functions of $n$* (each coefficient lies in $\mathbb Q(n)$). A priori individual shift terms can carry
a factor $(2n+t)^{-1}$ that is singular at $n=1$ (e.g. $2n-2$), but these factors cancel in the total — equivalently, the
sequence $n\mapsto [c^{2n}]F_0$ and $n\mapsto E(n)$ are both regular at every integer $n\ge1$, the latter because each
$Z_w(2n)$ and each coefficient $\rho_w(n)$ is (Lemma A2: poles only at half-integers $\le\tfrac12$). An identity of
rational functions that is regular at a point holds there; hence the identity, already proved for generic $n$, extends to
*every* integer $n\ge1$, $n=1$ included. (Numerically this is exactly what the exact checks confirm: $[c^{2}]F_0=\mu_{k,1}$
at every verified $k$.) For $k\le7$ no polynomial part occurs (Lemma A2 properness), so the $c\frac d{dc}$ clause is
vacuous there and only powers of $c$ appear as denominators.
*Step 2 (symmetrize).* $H_w(-c)=(-1)^{nz(w)}H_{\nu(w)}(c)$, $\nu$ = negate each nonzero letter (from (B1):
$[c^N]H_w(-c)=(-1)^N[c^N]H_w(c)=\operatorname{cf}_w(-\varepsilon_w)^NN^{-k_w}Z_{v_w}(N-1)$, and the inverse bijection at
$-\varepsilon_w$ negates all nonzero letters). Set $F_1=\frac12\big(F_0(c)+F_0(-c)\big)$; then $F_1$ is even,
$[c^{2n}]F_1=E(n)$, $[c^{2n+1}]F_1=0$.
*Step 3 (correct the low orders).* $F_1-\sum_n E(n)c^{2n}$ is an even Laurent series with zero positive-even
coefficients, hence a finite even Laurent polynomial $L(c)=\sum_{m=0}^{M}d_m c^{-2m}$ ($M$ bounded by the lowest pole
$e\ge-(k-1)$). Replace $F_1$ by $F:=F_1-L(c)$; now $[c^N]F=0$ for $N\le0$, $=E(N/2)$ for $N$ even $\ge2$, $=0$ for $N$ odd.
*Step 4 (coefficients).* Every operation is rational, so $r_j\in\mathbb Q(c)$; the only denominators are powers of $c$
(from the $c^e$ with $e<0$). $\blacksquare$

**Remark.** The correction $L(c)$ is the integration constant / negative-power part fixed by $M_k(0)=0$ (the flat angle
has $M=0$). For the single-zeta parts it is a bare constant $\mp\tfrac12$, matching the fitted forms
(`hpl_fit_k3_(1,0,0).json` contains `$1/2 * 1$`).

---

## Lemma D (weight bound)

$\operatorname{weight}(z_j)+|w_j|\le k+1$ for every term of the constructed form, for **all $k$**.

**Proof (repaired pole-budget, after the referee).** For a nested-sum term $\rho(n)Z_w(2n)$ and a point $p$, let
$\operatorname{ord}_p(\rho)$ be the pole order of $\rho$ at $p$ (0 if regular), and let $\nu(\rho)=\max_p\operatorname{ord}_p(\rho)$.
We show, at *every* location $p$,
$$|w|+\operatorname{ord}_p(\rho_w)\ \le\ k+1-\operatorname{weight}(z)\qquad(\ast)$$
for every term $\rho_w(n)Z_w(2n)$ of $[\mu_{k,n}]_z$; then $|w|+\nu(\rho_w)\le k+1-\operatorname{weight}(z)$.

The sources of poles in $\rho_w$ (Lemma A) split into two disjoint structural classes and the polygamma shifts:
* **the prefactor** $-1/(n(2n-1))$: simple poles at $\{0,\frac12\}$ only;
* **the Pochhammer ratio** $(n)_j/(2n)_j=\prod_{t=0}^{j-1}\frac{n+t}{2n+t}$: simple poles at $n=-t/2$, $t\ge1$ only —
  crucially, *no pole at $n=0$*, because the $t=0$ factor is $\frac{n}{2n}=\frac12$ (removable). Its pole set
  $\{-t/2:t\ge1\}$ is therefore **disjoint** from the prefactor's $\{0,\frac12\}$;
* **the polygamma shifts** $(n+t)^{-r}$ and $(2n+t)^{-r}$ from Lemma A(iv): pole order $r$ at $n=-t$ (integer) or
  $n=-t/2$ (half-integer). A shift arises exactly when a polygamma $\psi^{(r-1)}$ is *not* kept as a full harmonic
  $H^{(r)}\mapsto Z_{\pm r}$; kept, it contributes $r$ to $|w|$ and $0$ to the pole order; shifted, it contributes $0$
  to $|w|$ and $r$ to the pole order. Hence the total shift order, at any single location, is bounded by the part of the
  harmonic budget *not* used in $|w|$, namely $\le (k-\operatorname{weight}(z))-|w|$ (Lemma A1; the budget is additive
  under stuffle, which preserves $|w|$ and adds pole orders).

Now at each location:
* at $p\in\{0,\frac12\}$: $\operatorname{ord}_p\le 1+\big[(k-\mathrm{wt}(z))-|w|\big]$ (prefactor's single pole + shift budget);
* at $p=-t/2$ ($t\ge1$): $\operatorname{ord}_p\le 1+\big[(k-\mathrm{wt}(z))-|w|\big]$ (Pochhammer's single pole + shift budget);
* at $p=-t$ ($t\ge1$ integer): $\operatorname{ord}_p\le (k-\mathrm{wt}(z))-|w|$ (shift budget only; no structural pole).
In all three cases $|w|+\operatorname{ord}_p(\rho_w)\le |w|+1+(k-\mathrm{wt}(z)-|w|)=k+1-\mathrm{wt}(z)$, proving $(\ast)$.
(The "$+1$" is incurred at most once per location, since prefactor and Pochhammer poles are disjoint.)

Finally, in Lemma C the HPL word for a partial-fraction term of pole order $\kappa$ is $w_0=\mathrm{bji}(1,\kappa,w)$ of
weight $|w_0|=\kappa+|w|\le\nu(\rho_w)+|w|$ (Lemma B2), so
$$\operatorname{weight}(z)+|w_0|\ \le\ \operatorname{weight}(z)+\big(k+1-\operatorname{weight}(z)\big)=k+1.$$
And the $c\frac d{dc}$ branch never raises the word length: $c\frac d{dc}H_w=c f_{w_1}H_{w'}$ equals $H_{w'}$
(weight $|w|-1$) if $w_1=0$, and $(\frac1{1\mp c}-1)H_{w'}$ — a prefix-sum HPL of weight $|w'|+1=|w|$ — if $w_1=\pm1$.
So the bound is unaffected by the (vacuous for $k\le7$) polynomial-part clause. $\blacksquare$

**Status.** PROVED for all $k$ (the repair accounts for the Pochhammer poles, which FLD1's budget had forgotten).
Verified per word *and per pole location* at $k=7$: 0 violations, saturated (`FLD1R step_34`; sharpening our
`step_28`/budget check for $k\le6$).

---

## Lemma E (the letters, and the necessity of $-1$)

Only the letters $\{0,1,-1\}$ occur, because Lemma C only ever applies the bijection $\mathrm{bji}$ (whose outputs use
these letters) and the symmetrization $\nu$ (which permutes them). The letter $-1$ **is necessary**: it is produced
already by the even/odd split $H^{(r)}_{n-1}=2^{r-1}(Z_r+Z_{-r})-n^{-r}$, whose alternating sum $Z_{-r}$ maps to
$-1$-letters; in the constructed forms the rational part ($z=1$) uses $-1$ in $19/24$, $63/74$, $195/218$ of its terms
at $k=2,3,4$ (`step_28`).

## Lemma F (the S-units)

$H_w(c)$ with letters $\{0,\pm1\}$ are iterated integrals with kernels $d\log t$, $d\log(1-t)$, $d\log(1+t)$ at $t=c$,
i.e. multiple polylogarithms $\mathrm{Li}_{n_1,\dots,n_d}(a_1,\dots,a_d)$ with arguments in the group generated by
$c,\ 1-c,\ 1+c$ and $2$ (their symbol letters are exactly $c,1-c,1+c,2$; cf. the eight $\mathrm{Li}_3$ arguments
$c,-c,1-c,\frac1{1+c},\frac{1+c}2,\pm\frac{1-c}{1+c},\frac{2c}{1+c}$ in the paper). These are $S$-units of $\mathbb Q(c)$
with $S=\{\text{primes dividing }2,\,c,\,1-c,\,1+c\}$. Since $c=\cos\frac\alpha2$, this is the field
$\mathbb Q(\cos\frac\alpha2)$ of Conjecture conj:field.

---

## Assembly of the Theorem

$M_k=\sum_n\mu_{k,n}c^{2n}$ (thm:flat). Split each $\mu_{k,n}$ by $\zeta$-monomial $z$; by Lemma A, $[\mu_{k,n}]_z$ is a
proper nested-sum expression. By Lemma C it is the even part of an HPL combination with $\mathbb Q(c)$-coefficients and
power-of-$c$ denominators. Hence $M_k=\sum_z z\cdot(\text{HPL combination})=\sum_j z_j r_j(c)H_{w_j}(c)$ as power series,
hence as analytic functions on $|c|<1$ (both sides converge there), and by continuity on $[0,1)$. The weight bound is
Lemma D; the letters Lemma E; the $S$-unit content Lemma F. $\blacksquare$

---

## Lemma ledger

| Lemma | Statement | Status | Evidence |
|---|---|---|---|
| thm:flat | $M_k=\sum\mu_{k,n}c^{2n}$, $\mu_{k,n}\in\mathbb Q[\zeta(2..k)]$ | PROVED (paper) | cited |
| A (formula, A1) | $\mu_{k,n}$ = nested sums; $|w|\le k-\mathrm{wt}(z)$; $i+\ell=k$ exactly | PROVED (all $k$; verified $k\le7$) | `step_21/22/23`; referee `step_26/28/31/42` |
| A2 (poles) | poles in $\{\frac12,0\}\cup\{-\frac t2:t\le2k-2\}\subset\frac12\mathbb Z$ | PROVED (coarser form, all $k$); sharp set $\{t\le k-1\}$ PROVED FOR $k\le7$ | `step_24/31`; referee `step_32/34` |
| A2 (properness) | $\deg\rho_w=-2$ (denominators powers of $c$) | PROVED FOR $k\le7$; OPEN for $k\ge8$ (needs $Z$-sum independence lemma) | `step_24/31`; referee `step_32/34` (553 word-terms, all deg $-2$) |
| B | bijection $w\leftrightarrow(\varepsilon,k,v)$; $|w|=k+|v|$ | PROVED (all $k$) | `step_20` (242 words); referee re-verified independently (728 words len≤6, `step_24`) |
| C | even nested sums $\to$ HPL (algorithm terminates, is correct; regularity at $n=1$) | PROVED (all $k$); verified by construction FOR $k\le7$ | `step_25/32`; referee `step_35/36` (k=7) |
| D | $\mathrm{wt}(z)+|w|\le k+1$ | **PROVED for all $k$** (repaired pole-budget) | this proof; referee `step_34` (per word & pole, 0 violations, saturated) |
| E | letters $\{0,\pm1\}$; $-1$ needed | PROVED (all $k$) | `step_28`; at $k=7$, $-1$ in $4513/4678=96.5\%$ of rational terms (referee) |
| F | $S$-units $\langle2,c,1-c,1+c\rangle$ | PROVED | symbol letters; paper's $\mathrm{Li}_3$ arguments; referee (MPL args $\pm c,\pm1$) |

**The theorem (existence of the HPL form, letters, $S$-units, weight bound) is PROVED for all $k$.** Lemmas A, B, C, E, F
are written for general $k$; Lemma C's polynomial-part clause makes existence independent of the (unproven for $k\ge8$)
properness, costing at most $(1\pm c)$ denominators; and after the referee's repair, Lemma D's weight bound is itself
proved for all $k$. The only clauses not proved for all $k$ are the **power-of-$c$ denominators** and the **sharp pole
set**, both PROVED FOR $k\le7$ and open for $k\ge8$ (they would follow from properness via a $Z$-sum independence lemma).
The whole construction (symbolic $\mu$, $Z$-sum form, HPL form, all-$n$ identity) is verified exactly for $k\le7$; an
independent referee (FLD1R) re-verified Lemmas A, B, C, E, F for all $k$ and supplied the repairs to A2, C and D
incorporated here.

---

## Verification (see verify_fld1.py and result.md)

* **$k=2,3$ reproduction.** The construction reproduces Theorem thm:M2alpha (B2, 24 terms) and B3 (74/15/5 terms)
  *exactly* (Taylor agreement to $c^{40}$; `step_26`), and the all-$n$ Z-sum identities (even part $=\mu$, odd part $=0$)
  hold (`step_25`).
* **$k=4,5$ construction + validation.** HPL forms built by the algorithm: exact $[c^{2n}]=\mu_{k,n}$ for $n\le100$ (all
  $\zeta$-monomials, odd and non-positive parts $0$; `step_27`), and agreement with the corner integral
  (`corner_series`) at 4 angles to **$\ge40$ digits** ($\ge76$ digits at $2\pi/3,\pi/2,3\pi/4$ at N=600/dps=80; and
  **107 digits** at $\pi/3$ with N=2500/dps=110, confirming the form is exact there; `step_27`, `step_30`).
* **Consistency with G2_M4alpha.** Our $k=4$ equals G2's independently-fitted-and-proved B4 (R4, $\Psi_2$,
  $-\zeta(3)\Psi$, $\zeta(4)F_1$, $\zeta(2)^2F_2$) as functions (`step_29`, `step_29b`).
* **$k=6$ (beyond the required range).** Symbolic $\mu_{6,n}$ == exact table ($n\le60$); Z-sum form proper (deg $-2$),
  weight budget $|w|+\nu\le 7-\mathrm{wt}(z)$ (saturated), poles $\{-\frac52,\dots,0,\frac12\}$; and the full
  construction + all-$n$ identity (even $=\mu$, odd $=0$) for all 11 $\zeta$-monomials (`step_31`, `step_32`).
* **$k=7$ (referee FLD1R).** Symbolic $\mu_{7,n}$ (429 monomials) == exact table ($n\le30$, two independent
  implementations identical); the Z-sum canonical form (553 word-terms) is proper (all $\deg=-2$) with poles
  $\{-3,-\frac52,-2,-\frac32,-1,-\frac12,0,\frac12\}$; weight budget per word *and per pole location* has 0 violations
  (saturated at $8-\mathrm{wt}(z)$); the full construction (6553 HPL terms; rational part 4678 terms) satisfies even $=\mu$
  and odd $=0$ as $\mathbb Q(n)$ identities and exact $[c^{2n}]F=\mu_{7,n}$ for $n\le30$; and the HPL form agrees with the
  direct corner integral $M_7(2\pi/3)$ to residual $\mathbf{3.5\cdot10^{-46}}$ ($\sim45$ digits, dps=50, N=500)
  (referee `step_26`–`step_39`). This is the computational content of "power-of-$c$ denominators PROVED FOR $k\le7$".
