# SQSRR — referee report and completion: the squircle sine, $b_n''(2)$, $R$ and $b_1'''(2)$

Referee of `new_results/SQSR_sine_series/result.md` (I did not write it). Everything below was
recomputed with my own code in this folder; I read SQSR's scripts but did not reuse them. The only external number I use is SQSR's 260-digit value of $R$, and only as a check.
Python: `.venv/bin/python` (mpmath, sympy, python-flint).

## 0. Verdicts

| Claim | SQSR status | Verdict |
|---|---|---|
| C1: law for $b_n''(2)$, all odd $n\ge3$ | PROVED | **CONFIRMED AFTER REPAIR** (the proof is right; the regularity paragraph must be written out, Sec. 2.3, and the period constant must be named $q=\frac p2\sin\frac\pi p$, not "$\pi/\pi_p$", Sec. 2.1). Independent check from the definition at $n=3,5,23,41,99$: agreement $\le3\cdot10^{-60}$. |
| C2a: $\sum_{n\ge3}b_n'C_n=(1-\zeta(3))/64$ | PROVED | **CONFIRMED** (independent exact derivation, `r3_R_exact.py`). |
| C2b: $R$ in the $\mathbb Q$-span of level-2 MZVs of weight $\le5$ | PROVED | **CONFIRMED** (pole-order count, Sec. 3); superseded by C3. |
| C3: closed form of $R$ | NUMERICAL IDENTIFICATION | **NOW PROVED (computer-assisted exact)**, Sec. 4. Exactly the SQSR value. $\mathrm{Li}_4(\frac12)$, $\mathrm{Li}_5(\frac12)$ never enter: no weight-4 double sum occurs in $R$ at all. |
| C4: closed form of $b_1'''(2)$ | follows from C3 | **PROVED (computer-assisted exact)**, given C1, C3 and Parseval (Sec. 5). Independent check straight from the definition, without Parseval: agreement $3.5\cdot10^{-57}$. |
| SQSR's table entry "R to 250 digits: PROVED (numerical, rigorous method)" | — | Over-stated: Euler–Maclaurin without an error bound is *numerical*. Harmless now, since $R$ is exact. |
| Paper sentence "$\frac3{64}\pi^2\zeta(3)$ from the products $b_n'b_n''$" | — | **Wrong as a statement about $b_1'''$** (the products give $\frac9{256}\pi^2\zeta(3)$ in total, plus $-\frac{93}{256}\zeta(5)$ and $-\frac{141}{64}\zeta(3)$). SQSR's replacement has the right numbers but the wrong status ("identified numerically"); corrected text in Sec. 6. |

Status separation for the final result:
- **PROVED (by hand, checkable):** ODE and its expansion; Fourier convolution formulas; $(1-n^2)b_n^{(k)}=\widehat{F_k}(n)$; Parseval at third order; reduction of $R$ to the eight double sums of Eq. (4.2).
- **COMPUTER-ASSISTED EXACT:** the symbolic summation giving C1 in closed form (SQSR `symsum.py`/`s02`, reviewed line by line); the reduction of $R$ to Eq. (4.2) and the evaluation of those double sums by double shuffle + distribution (`r3_R_exact.py`, `dshuffle.py`); the assembly of $b_1'''$ (`r5_b1ppp_exact.py`).
- **NUMERICAL (checks only, nothing depends on them):** definition-based quadrature of $b_n'$, $b_n''$, $b_1'''$; 70-digit checks of every double-sum reduction; direct summation of $R$.

## 1. Files

| file | what |
|---|---|
| `r1_ode_norm.py/.log` | sympy: $q(p)$, $q^{-p}$ to $O(e^3)$, the order-1 and order-2 forcings, $N',N'',N'''$ from the Gamma formula, $N(3)$ by quadrature |
| `r0_L0sq.py/.log` | Fourier coefficients of $\log^2\tan\theta$ (ingredient of $F_2$) |
| `r2_defquad.py`, `r2_defquad_A/B.log`, `law.py`, `r2b_compare.py/.log` | $b_n'(2)$, $b_n''(2)$ **from the definition** (incomplete Beta, fixed-node composite Gauss–Legendre, 7-point FD in $p$, arb balls at 480 bits), two independent settings; comparison with C1 |
| `r2c_b1ppp.py`, `r2c_b1ppp_A/B.log` | same with 9-point stencils: $b_1'(2)$, $b_1''(2)$, $b_1'''(2)$, $b_3'''(2)$ from the definition |
| `dshuffle.py/.log` | exact reduction of all depth-2 level-2 sums of weight 2–5 by stuffle, shuffle (incl. linearly regularised) and distribution |
| `r3_R_exact.py/.log`, `r3_R_exact.pkl` | exact $R$, $\sum b_n'C_n$, $\sum_{n\ge3}b_n'^2$ from the law |
| `r3b_R_in_Z.py/.log` | $R$ as an explicit combination of eight double sums before reduction |
| `r4_numeric.py/.log` | 70+-digit numerical checks of every reduction and of $R$ |
| `r5_b1ppp_exact.py/.log` | exact assembly of $b_1''(2)$ and $b_1'''(2)$; odd-weight bookkeeping |

## 2. Part 1 — C1 line by line

### 2.1 Definition, period, ODE
Paper (Sec. "Other Deformations"): $F_p(y)=\int_0^y(1-t^p)^{-1/p}dt$, $\Phi_p=\frac\pi2F_p/F_p(1)$, $S_p=\Phi_p^{-1}$ on $[0,\frac\pi2]$, extended like $\sin$;
$b_n(p)=\frac4{n\pi}\int_0^1\cos(n\Phi_p(y))\,dy$. Note $\Phi_p(y)=\frac\pi2I_{y^p}(\frac1p,1-\frac1p)$ (regularised incomplete Beta), which is what my quadrature uses.

With $\sin_p=F_p^{-1}$, $S_p(\theta)=\sin_p(2F_p(1)\theta/\pi)$, so $\sin_p'=qS_p'$ with
$q=\pi/(2F_p(1))=\frac p2\sin\frac\pi p$ ($F_p(1)=\frac1pB(\frac1p,1-\frac1p)=\pi/(p\sin\frac\pi p)$).
From $|\sin_p|^p+|\sin_p'|^p=1$: $|S|^p+|qS'|^p=1$, and on $(0,\frac\pi2)$ (where $S,S'>0$) differentiating gives
$S''=-q^{-p}(S/S')^{p-2}S$. **Correct.**
*Notational repair.* SQSR (and AV8) write $q=\pi/\pi_p$ with $\pi_p:=2F_p(1)$ (Lindqvist's constant). The paper's $\pi_p$ (Sec. on the $\ell^p$-length, Adler–Tanton/Keller–Vakil) is a different number ($\pi_1=4$ there, while $2F_1(1)=\infty$). The mathematics uses the right $q$; the symbol must not travel into the paper.
`r1`: $q(2)=1$, $q^{-p}=1-e+(\frac14+\frac{\pi^2}{16})e^2+(\frac18-\frac{3\pi^2}{32})e^3+O(e^4)$, $e=p-2$.

### 2.2 Expansion
Independent sympy expansion (`r1`) of $-q^{-p}\exp(e\log(S/S'))S$ with $S=\sin\theta+es_1+\frac{e^2}2s_2$:
order 1 gives $s_1''+s_1=(1-L_0)\sin\theta$, order 2 gives exactly SQSR's
$F_2=-2L_0s_1+2\tan\theta\,s_1'-L_0^2\sin\theta+2L_0\sin\theta-(\frac12+\frac{\pi^2}8)\sin\theta$ (difference simplifies to 0). **Correct.**

### 2.3 Regularity (the repair)
SQSR states in one sentence that $s_1,s_2$ are $C^1$ with integrable second derivative and that $(1-n^2)b_n^{(k)}=\widehat{F_k}(n)$. This needs: (i) $b_n^{(k)}(2)=\frac4\pi\int_0^{\pi/2}s_k\sin n\theta\,d\theta$ with $s_k=\partial_p^kS_p|_{p=2}$; (ii) integration by parts twice with vanishing boundary terms. Argument to be written in the notes:
- $\Phi_p(y)$ is real-analytic in $(p,y)$ on $J\times(0,1)$, $J=[\frac32,\frac52]$, with $\partial_y\Phi_p=q(1-y^p)^{-1/p}>0$, so $S_p(\theta)$ is real-analytic in $(p,\theta)$ on $J\times(0,\frac\pi2)$; $p$- and $\theta$-derivatives commute there, so $s_k''=\partial_p^k S_p''$ and the expansion of 2.2 is legitimate pointwise.
- Near $\theta=\frac\pi2$ ($u=\frac\pi2-\theta$): $\frac\pi2-\Phi_p=\frac\pi2I_x(1-\frac1p,\frac1p)=x^{1-1/p}g(p,x)$ with $x=1-y^p$ and $g$ analytic, $g(p,0)\ne0$. Inverting, $S_p=\Psi(p,u^{p/(p-1)})$ with $\Psi$ analytic and $\Psi(p,0)=1$. Hence $\partial_p^kS_p=O(u^2|\log u|^k)$, $\partial_\theta\partial_p^kS_p=O(u|\log u|^k)$, $\partial_\theta^2\partial_p^kS_p=O(|\log u|^k)$, uniformly for $p\in J$.
- Near $\theta=0$: $\Phi_p(y)=c_py\,{}_2F_1(\frac1p,\frac1p;1+\frac1p;y^p)$ is analytic in $(p,y,y^p)$, so $S_p=\theta/c_p+O(\theta^{p+1})$ and $\partial_p^kS_p=O(\theta)$, the non-smooth part being $O(\theta^3|\log\theta|^k)$.
- So $|\partial_p^kS_p(\theta)|\le C$ on $[0,\frac\pi2]\times J$: dominated convergence gives (i) for $b_n(p)=\frac4\pi\int S_p\sin n\theta$ and for $N(p)=\frac4\pi\int S_p^2$. Also $s_k\in C^1[0,\frac\pi2]$, $s_k''\in L^1$, $s_k(0)=0$, $s_k'(\frac\pi2)=0$; with $\sin0=\cos\frac{n\pi}2=0$ both boundary terms vanish, giving $(1-n^2)b_n^{(k)}=\frac4\pi\int F_k\sin n\theta$ for $n\ge3$ ($F_k\in L^1$: worst terms $\log^2$ at $0$, $\log$ at $\frac\pi2$).
(Alternatively, for (i): $\partial_p^k\Phi_p(y)$ is continuous on $[0,1]\times J$ by the same local forms, so one may differentiate AV8's $\int_0^1\cos(n\Phi_p)dy$ directly.)

### 2.4 Fourier coefficients of $F_2$
Checked by hand:
- $\log\tan\theta=-2\sum_{k\ \rm odd}\cos(2k\theta)/k$ (classical). $\widehat{L_0\sin}(n)=2/(n+1)$ ($n\equiv1$), $-2/(n-1)$ ($n\equiv3$). Correct.
- $L_0s_1$: $\cos(2k\theta)\sin m\theta=\frac12[\sin(m+2k)\theta+\sin(m-2k)\theta]$ gives the coefficient $-\sum_{k\ \rm odd}(\beta_{n-2k}+\beta_{n+2k})/k$, $\beta_{-m}=-\beta_m$; the double series converges absolutely ($\sum|\beta_m|<\infty$). Correct; the bookkeeping of the three ranges of $n-2k$ in `s02` (positive, $\pm1$ giving $\pm\beta_1/(2t+1)$, negative) is correct in both classes.
- $h=\tan\theta\,s_1'$: $h$ is odd and symmetric about $\frac\pi2$ (odd sine series), $\cos\theta\,h=\sin\theta\,s_1'$ gives $h_{2J-1}+h_{2J+1}=c_{J-1}-c_J$ ($J\ge1$, $c_i=(2i+1)\beta_{2i+1}$); SQSR's $h_{2j+1}=c_j+2\sum_{r\ge1}(-1)^rc_{j+r}$ satisfies it (telescoping check), and the homogeneous solution $(-1)^jC$ is not in $\ell^2$ while $h\in L^2$ (only a $\log$ singularity at $\frac\pi2$). Correct; `s02`'s split into $r$ even (same class) and $r$ odd (other class) is correct.
- $L_0^2=\lambda_0+\sum\lambda_m\cos2m\theta$, $\lambda_m=8O_m/m$ ($m$ even), $0$ ($m$ odd): checked by quadrature to $2\cdot10^{-43}$ (`r0`); then $-\widehat{L_0^2\sin}(n)$ is $-8O_{(n-1)/2}/(n-1)$ or $8O_{(n+1)/2}/(n+1)$ as in `s02`. Correct.
- `symsum.py`: $\zeta(s,M+\frac12)=(2^s-1)\zeta(s)-2^sH^{(s)}_{2M}+H^{(s)}_M$, $-\psi(M+\frac12)=\gamma+2\log2-2H_{2M}+H_M$, integer shifts, range reversal in `sum_fin`, and the requirement that $s=1$ coefficients cancel — all correct. The output (`s02b`) is exactly the law C1 in both classes (I re-derived the class-wise coefficients from the task's formula: they match term by term).

### 2.5 Independent check from the definition (`r2*`)
No SQSR or AV8 code or data. $b_n(p)$ on a fixed node set (composite Gauss–Legendre, 40 or 48 nodes per panel, 256 or 200 panels, 140–150 geometric levels at both ends), $p=2+jh$ ($h=10^{-10}$ or $10^{-11}$), central differences of order 6, arb balls at 480 bits.

| $n$ | $b_n''(2)$ (definition) | quad − law, run A | run B |
|---|---|---|---|
| 3 | 0.238425455580538678951779367627325221224477313 | $1.3\cdot10^{-61}$ | $1.3\cdot10^{-61}$ |
| 5 | 0.0281073109426667352443504186445502529451700038 | $2.9\cdot10^{-60}$ | $-4.3\cdot10^{-62}$ |
| 23 | −0.000506845974021496356318988019945478942718104493 | $2.2\cdot10^{-62}$ | $4.4\cdot10^{-64}$ |
| 41 | 0.000140078744040196768386136107986018005448713808 | $7.0\cdot10^{-63}$ | $2.8\cdot10^{-65}$ |
| 99 | −0.000012846940906336124896851291034957677860137604 | $-8.1\cdot10^{-65}$ | $1.9\cdot10^{-65}$ |

$b_n'(2)$ agrees with $\varepsilon/((n^2-1)(2\lfloor\frac{n-1}4\rfloor+1))$ to $\le2\cdot10^{-60}$ in all cases. (So the law holds to $\ge57$ significant digits at $n=23,41,99$.)

**C1 verdict: CONFIRMED AFTER REPAIR** (write out 2.3; rename $q$).

## 3. Part 1 — C2 and the third-order Parseval identity

**Norm.** $N(p)=\sum_nb_n(p)^2=\frac1\pi\int_0^{2\pi}S_p^2=\frac4\pi\int_0^{\pi/2}S_p^2d\theta=\frac4\pi\int_0^1y^2\Phi_p'(y)dy=\frac{2}{F_p(1)}\cdot\frac1pB(\tfrac3p,1-\tfrac1p)=\frac{2\sin(\pi/p)}\pi\frac{\Gamma(3/p)\Gamma(1-1/p)}{\Gamma(1+2/p)}$,
the paper's formula (`r1`: quadrature at $p=3$ agrees). $N'(2)=\log2-1$, $N''(2)=\log^22-3\log2+\frac{5\pi^2}{24}$, $N'''(2)=\frac1{16}[87\zeta(3)+16\log^32-96\log^22+72\log2+10\pi^2\log2-20\pi^2]$ (sympy, exact; matches the claimed values to $10^{-160}$).

**Parseval at third order (derivation).** By 2.3 we may differentiate $N(p)=\frac4\pi\int_0^{\pi/2}S_p^2$ three times under the integral:
$N'''(2)=\frac4\pi\int(2S\,s_3+6s_1s_2)$ with $S=\sin\theta$. Now $\frac4\pi\int\sin\theta\,s_3=b_1'''(2)$, and $\frac4\pi\int s_1s_2=\sum_nb_n'(2)b_n''(2)$ by Parseval for the $L^2$ functions $s_1,s_2$ (both bounded). Hence
$$b_1'''(2)=\tfrac12N'''(2)-3b_1'(2)b_1''(2)-3\sum_{n\ge3}b_n'(2)b_n''(2).$$
No term-by-term differentiation of $\sum b_n^2$ is used. Second order likewise: $b_1''=\frac12N''-\sum_nb_n'^2$.
**Convergence.** $b_n'=O(n^{-3})$; from C1, $A_n$ and $C_n\pi^2$ are each $O(n^{-2})$ (they cancel to $b_n''=O(n^{-3}\log n)$), so $\sum b_n'A_n$, $\sum b_n'C_n$, $\sum b_n'^2$ all converge absolutely ($O(n^{-5})$).

**C2a.** `r3_R_exact.py` (my pipeline, Sec. 4) gives exactly $\sum_{n\ge3}b_n'C_n=\frac1{64}-\frac{\zeta(3)}{64}$ and $\sum_{n\ge3}b_n'^2=\frac{5\zeta(2)}{16}+\frac{15\zeta(4)}{128}-\frac58$; direct summation agrees (Sec. 4.4). **CONFIRMED.**
**C2b.** With $n=4k+\varepsilon$, $a=4k$, $c=4k+2\varepsilon$: $\beta_n=2\varepsilon/(ac^2)$ has pole orders $(1,2)$ at $k=0,-\frac\varepsilon2$; the coefficients of $H_{a/4}$, $H_{a/2}$, $\sum_{i\le a/4}(2i-1)^{-2}$ and $\rho_n$ have orders $(3,0)$, $(3,2)$, $(2,1)$, $(2,3)$. So $\beta_nA_n$ is a sum of $X_k/D^j$ with $D\in\{k,2k\pm1\}$ and weight $j+\mathrm{wt}(X)\le5$. **CONFIRMED.**

## 4. Part 2 — proof of C3

### 4.1 Reduction to alternating double sums (exact; `r3_R_exact.py`)
Notation: $Z(a,\sigma;b,\tau)=\sum_{m>l\ge1}\sigma^m\tau^l m^{-a}l^{-b}$ (write $Z(a^\pm,b^\pm)$), $L(s,\sigma)=\sum_m\sigma^mm^{-s}$; $A_m=\sum_{l\le m}(-1)^{l+1}/l$, $A^{(2)}_m$ likewise.
1. $R=\sum_{\varepsilon=\pm1}\sum_{k\ge1}\beta_n\big[\rho_n-\frac{2\varepsilon}{a^3}H_k+\frac{2\varepsilon(3n^2+2\varepsilon n+3)}{a^3c^2}H_{2k}+\frac{2\varepsilon n}{a^2c}O_k\big]$, $O_k=\sum_{i\le k}(2i-1)^{-2}$ (C1 with $\tau=\frac{\pi^2}8-O_k$, the $\pi^2$ and $\log2$ parts removed).
2. Partial fractions in $k$; poles only at $k=0,\pm\frac12$.
3. Each term $X_k/D^j$ becomes a parity-weighted sum over $m=2k$ ($D=k$) or $m=2k\pm1$, using the exact identities (for the relevant parity of $m$):
   $H_k=H_m-A_m$ ($m=2k$ or $2k+1$), $H_k=H_m-A_m+\frac2{m+1}$ ($m=2k-1$); $H_{2k}=H_m$, $H_m-\frac1m$, $H_m+\frac1{m+1}$; $O_k=\frac12(H^{(2)}_m+A^{(2)}_m)$, minus $\frac1{m^2}$ for $m=2k+1$.
   Rational leftovers with a pole at $m=-1$ are re-indexed ($m'=m+1$, parity flips).
4. $\sum_m\sigma^mH_m/m^s=Z(s,\sigma;1,+)+L(s+1,\sigma)$, $\sum\sigma^mA_m/m^s=-Z(s,\sigma;1,-)-L(s+1,-\sigma)$, and the same with $H^{(2)},A^{(2)}$.
5. Divergent pieces: all partial sums are truncated at $M=2N+O(1)$; with $T:=H_M$, $\sum_{m\le M}1/m=T$ and, by the truncated stuffle, $Z_M(1,+;b,\tau)=TL(b,\tau)-Z(b,\tau;1,+)-L(b+1,\tau)+o(1)$, $Z_M(1,+;1,+)=\frac12(T^2-\zeta(2))+o(1)$. Since $R$ converges, the $T$-dependence must vanish identically; it does, already before any reduction (`r3b`: $T$ and $T^2$ coefficients $0$). This is a strong internal check of steps 2–5.

Result (`r3b_R_in_Z.log`), exact:
$$R=\tfrac12Z(2^+,1^+)-\tfrac14Z(2^+,1^-)-\tfrac18Z(2^-,1^-)+\tfrac1{16}Z(3^+,2^+)-\tfrac1{16}Z(3^+,2^-)+\tfrac18Z(4^+,1^+)-\tfrac1{16}Z(4^+,1^-)-\tfrac1{16}Z(4^-,1^-)$$
$$\qquad+\tfrac5{16}-\tfrac{13}{16}\zeta(2)+\tfrac9{32}\zeta(3)-\tfrac{15}{64}\zeta(4)+\tfrac{31}{256}\zeta(5).\qquad(4.2)$$
No double sum of weight 4 occurs, so $\mathrm{Li}_4(\frac12)$ cannot enter; all double sums have odd weight (3 or 5).

### 4.2 Evaluation of the eight double sums (proved relations; `dshuffle.py`)
For each weight $w$, all $Z(a^\sigma,b^\tau)$, $a+b=w$, $(a,\sigma)\ne(1,+)$, are unknowns, and I impose:
- (St) stuffle $L(a,\sigma)L(b,\tau)=Z(a^\sigma,b^\tau)+Z(b^\tau,a^\sigma)+L(a+b,\sigma\tau)$ — exact for truncated sums, valid whenever all series converge;
- (Sh) shuffle, from the partial fractions $\frac1{m^an^b}=\sum_{j=1}^a\binom{w-j-1}{b-1}\frac1{m^j(m+n)^{w-j}}+\sum_{j=1}^b\binom{w-j-1}{a-1}\frac1{n^j(m+n)^{w-j}}$:
  $L(a,\sigma)L(b,\tau)=\sum_{j\le a}\binom{w-j-1}{b-1}Z((w-j)^\tau,j^{\sigma\tau})+\sum_{j\le b}\binom{w-j-1}{a-1}Z((w-j)^\sigma,j^{\sigma\tau})$; proved with weights $x^{m+n}$ (absolute convergence) and $x\to1^-$ (Abel);
- (ShR) the same when one factor is $L(1,+)$: the divergence is linear, Abel regularisation ($T=-\log(1-x)$) and harmonic regularisation ($T=H_N$) agree for linearly divergent series, so the relation holds with $L(1,+):=T$ and the truncated-stuffle value of $Z(1^+,b^\tau)$; the $T$-parts match identically (asserted in code);
- (Di) distribution $\sum_{\sigma,\tau}Z(a^\sigma,b^\tau)=2^{2-w}Z(a^+,b^+)$, $a\ge2$ (exact for truncated sums);
- depth one: $L(s,+)=\zeta(s)$, $L(s,-)=-(1-2^{1-s})\zeta(s)$, $L(1,-)=-\log2$; Euler $\zeta(2)^2=\frac52\zeta(4)$.

The linear system has a unique solution at weights 2, 3 and 5 and one free element at weight 4 ($Z(3^-,1^-)$, which carries $\mathrm{Li}_4(\frac12)$); this is the classical picture (odd weight reduces to depth one: Borwein–Borwein–Girgensohn 1995, Flajolet–Salvy 1998 Thm 7.2). The values needed in (4.2):

| sum | value |
|---|---|
| $Z(2^+,1^+)=\zeta(2,1)$ | $\zeta(3)$ (Euler) |
| $Z(2^+,1^-)$ | $\zeta(3)-\frac32\zeta(2)\log2$ |
| $Z(2^-,1^-)$ | $\frac32\zeta(2)\log2-\frac{13}8\zeta(3)$ |
| $Z(3^+,2^+)=\zeta(3,2)$ | $3\zeta(2)\zeta(3)-\frac{11}2\zeta(5)$ (Euler) |
| $Z(3^+,2^-)$ | $\frac14\zeta(2)\zeta(3)-\frac{21}{32}\zeta(5)$ |
| $Z(4^+,1^+)=\zeta(4,1)$ | $2\zeta(5)-\zeta(2)\zeta(3)$ (Euler) |
| $Z(4^+,1^-)$ | $2\zeta(5)-\frac38\zeta(2)\zeta(3)-\frac{15}8\zeta(4)\log2$ |
| $Z(4^-,1^-)$ | $\frac{15}8\zeta(4)\log2+\frac34\zeta(2)\zeta(3)-\frac{91}{32}\zeta(5)$ |

Each value is *derived* here from the relations listed, not copied from a table. The three non-alternating ones are Euler's. The alternating ones are instances of the classical odd-weight reduction of linear alternating Euler sums (Sitaramachandrarao 1987; Borwein–Borwein–Girgensohn 1995; Flajolet–Salvy 1998, Thm 7.2; tables in Borwein–Bradley–Broadhurst 1997). I cite these for orientation only: I had no access to the papers (no network), so I do not claim row-by-row agreement with their printed tables. All 31 determined values at weights 2–5 are re-verified numerically (Sec. 4.4).

### 4.3 Result
Substituting (`r3_R_exact.log`):
$$R=\frac5{16}-\frac{13}{16}\zeta(2)+\frac{47}{64}\zeta(3)+\frac3{16}\zeta(2)\log2-\frac{15}{64}\zeta(4)+\frac{31}{256}\zeta(5)+\frac3{128}\zeta(2)\zeta(3),$$
i.e. $R=\frac5{16}-\frac{13\pi^2}{96}+\frac{47}{64}\zeta(3)+\frac{\pi^2\log2}{32}-\frac{\pi^4}{384}+\frac{31}{256}\zeta(5)+\frac{\pi^2\zeta(3)}{256}$; `R - C3 = 0` exactly.
**C3 is PROVED (computer-assisted exact):** an exact symbolic derivation by script from C1 and from stuffle, shuffle (Abel), linearly regularised shuffle, distribution and $\zeta(2)^2=\frac52\zeta(4)$. Nothing stays numerical. $\mathrm{Li}_4(\frac12)$ and $\mathrm{Li}_5(\frac12)$ do not cancel — they never appear.

### 4.4 Numerical re-verification (`r4_numeric.log`)
- Every one of the 31 determined double sums of weights 2–5 (and the 10 weight-4 relations, with the free element evaluated numerically) agrees with an independent evaluation from $Z(a^\sigma,b^\tau)=\frac{\sigma}{\Gamma(a)}\int_0^1(-\log u)^{a-1}\frac{\mathrm{Li}_b(\sigma\tau u)}{1-\sigma u}du$ at 90 digits: worst residual $7\cdot10^{-91}$.
- The classical $\sum_{n\ge1}(-1)^{n+1}H_n/n^3=\frac{11}4\zeta(4)-\frac74\zeta(3)\log2+\frac{\pi^2}{12}\log^22-\frac{\log^42}{12}-2\mathrm{Li}_4(\frac12)$ (the free weight-4 element; *not used*) checks to $4\cdot10^{-91}$.
- Exact $R$ vs SQSR's 260-digit Euler–Maclaurin value (computed from the law by a different route, no reduction): $8\cdot10^{-264}$.
- Fully independent of the law and of Parseval: the definition-based $b_1'''(2)$ (Sec. 5) checks $R$ to about 56 digits.

## 5. Part 3 — the theorem for $b_1'''(2)$

`r5_b1ppp_exact.log`: from $N''$, $N'''$ (Gamma), $b_1'=\frac{\log2-1}2$, $b_1''=\frac12N''-b_1'^2-\sum_{n\ge3}b_n'^2=\frac38-\log2+\frac{\log^22}4+\frac{5\pi^2}{96}-\frac{\pi^4}{768}$, $\sum_{n\ge3}b_n'b_n''=\log2\sum b_n'^2+\pi^2\sum b_n'C_n+R$:

**Theorem.** For the squircle sine $S_p=\sum b_n(p)\sin n\theta$,
$$b_1'''(2)=-\frac38+\frac{33}{16}\log2-\frac{3\pi^2}{16}-\frac98\log^22+\frac{33}{64}\zeta(3)-\frac{\pi^2\log2}{64}+\frac{\log^32}8+\frac{3\pi^4}{512}-\frac{93}{256}\zeta(5)+\frac9{256}\pi^2\zeta(3)-\frac{\pi^4\log2}{512}$$
$=-0.30262285104154317677477849746166540942307522489364848612816171164\ldots$
**Status: PROVED (computer-assisted).** Ingredients: C1 (proved; regularity as in 2.3), Parseval at third order (Sec. 3), the exact reduction of Sec. 4. `b_1''' - C4 = 0` exactly.
Independent check from the definition (no Parseval, no law; `r2c`): run A $-0.302622851041543176774778497461665409423075224893648486128127$, i.e. $3.5\cdot10^{-57}$ from the theorem (run B, $h=10^{-9}$: $3\cdot10^{-53}$). AV8's quadrature (44 digits) and SQSR's 250-digit series value also agree.

Odd-weight constants by source (in the basis $\zeta(3),\pi^2\zeta(3),\zeta(5)$):

| | $\frac12N'''$ | $-3\pi^2\sum b_n'C_n$ | $-3R$ | total |
|---|---|---|---|---|
| $\zeta(3)$ | $\frac{87}{32}$ | 0 | $-\frac{141}{64}$ | $\frac{33}{64}$ |
| $\pi^2\zeta(3)$ | 0 | $\frac3{64}$ | $-\frac3{256}$ | $\frac9{256}$ |
| $\zeta(5)$ | 0 | 0 | $-\frac{93}{256}$ | $-\frac{93}{256}$ |

Caveat for wording: "the coefficient of $\pi^2\zeta(3)$" is read off this closed form; that it is *the* coefficient presumes the $\mathbb Q$-linear independence of the constants involved, which is believed but not proved (the paper already flags this kind of assumption elsewhere). The split between "$\pi^2$ parts" and $R$ is not intrinsic: in C1, $\pi^2$ enters only through the tail $\tau=\frac{\pi^2}8-O_k$.

## 6. Consequence for the paper

Current text (main.tex, Sec. "Other Deformations"): "… Beyond Gamma, $b_1'''(2)$ gets, besides the $\frac{87}{32}\zeta(3)$ from $N'''$, a term $\frac3{64}\pi^2\zeta(3)$ from the products $b_n'b_n''$, through the $\pi^2$ parts of the $b_n''$ (numerical, from the second-order forms for $n\le17$). One more series in $b_1'''$ we have not evaluated, and it could change this coefficient."
This is wrong as a statement about $b_1'''$: the products give $\frac9{256}\pi^2\zeta(3)$ in total, they also give $\zeta(5)$, and they change the $\zeta(3)$ coefficient from $\frac{87}{32}$ to $\frac{33}{64}$.

**Is SQSR's proposed sentence exactly right?** Its constants are all right (now proved). Its status is not: "identified numerically to 250 digits" and the row "$b_1'''(2)$ closed form: Numerical (250 digits)" must become "proved (computer-assisted)". It is right that the closing "$\zeta(3)$ waits for the third" must also mention $\pi^2\zeta(3)$ and $\zeta(5)$. Its table entry "R to 250 digits: PROVED (numerical, rigorous method)" over-states the case, because there is no error bound.

**Corrected sentence** (replaces "Beyond Gamma … could change this coefficient."):
> Beyond Gamma, the products $b_n'b_n''$ bring in more odd zeta values. We prove the second-order law for every odd $n\ge3$$^\dagger$; with it, $\sum b_n'b_n''$ reduces to alternating double Euler sums of weights three and five, and these evaluate exactly:
> $b_1'''(2)=-\frac38+\frac{33}{16}\log2-\frac{3\pi^2}{16}-\frac98\log^22+\frac{33}{64}\zeta(3)-\frac{\pi^2\log2}{64}+\frac{\log^32}8+\frac{3\pi^4}{512}-\frac{93}{256}\zeta(5)+\frac9{256}\pi^2\zeta(3)-\frac{\pi^4\log2}{512}$
> (proved, computer-assisted). So the $\frac{87}{32}\zeta(3)$ from $N'''$ becomes $\frac{33}{64}\zeta(3)$, and $\pi^2\zeta(3)$ and $\zeta(5)$ come in at the same order; $\mathrm{Li}_4(\frac12)$ and $\mathrm{Li}_5(\frac12)$ do not.

Two neighbouring phrases must follow: "(numerical for $n\le17$; $b_1''$ is proved with Parseval's identity)" → "(proved for every $n$$^\dagger$; $b_1''$ with Parseval's identity)"; and the closing "$\zeta(3)$ waits for the third" → "odd zeta values wait for the third".

**Table II row** (replaces the current one):
> Squircle sine near the circle: $b_n'(2)$ and $b_n''(2)$ laws (no odd zeta values through second order), closed form of $b_1'''(2)$ with $\zeta(3)$, $\pi^2\zeta(3)$, $\zeta(5)$ (Sec.~\ref{sec:otherdef}) & $b_n'(2)$, $b_1''(2)$, $b_n''(2)$ (all odd $n\ge3$): Proved$^\dagger$; $N'''(2)$: Proved; $b_1'''(2)$: Proved (computer-assisted) & $b_n'$: 15 (independent, $n\le7$); $b_n''$: 57 (independent, from the definition, $n=23,41,99$); $b_1'''$: 56 (independent, from the definition); $N'''(2)$: 20 (independent) \\
