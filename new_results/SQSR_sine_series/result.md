# SQSR — the series R and the exact value of $b_1'''(2)$ for the squircle sine

Setup and notation as in AV8: $S_p(\theta)=\sin_p(\pi_p\theta/\pi)=\sum_{n\ \rm odd}b_n(p)\sin n\theta$.
Primes are $p$-derivatives at $p=2$. $\beta_n:=b_n'(2)$. For odd $n\ge3$ put
$\varepsilon=(-1)^{(n-1)/2}$, $a=n-\varepsilon$ (a multiple of 4), $c=n+\varepsilon$, so $\beta_n=2\varepsilon/(ac^2)$.

## Status table

| Statement | Status | Evidence |
|---|---|---|
| ODE form: $b_n''=F_{2,n}/(1-n^2)$, $n\ge3$ (Sec. 1) | PROVED | derivation below; quadrature of $F_2$ reproduces AV8's $b_3''$ |
| Law for $b_n''(2)$, all odd $n\ge3$ (Thm 1) | PROVED | exact symbolic summation (`s02`); agrees with AV8 quadrature to $\le4\cdot10^{-56}$ for $n=3..17$, and with fresh quadrature to $\le5\cdot10^{-46}$ for $n=19,21,31,33,61$ (`s05`) |
| $\pi^2$ part $C_n$, $\log2$ part $=\beta_n$ (AV8 laws) | PROVED | special cases of Thm 1 |
| $\sum_{n\ge3}\beta_nC_n=(1-\zeta(3))/64$ | PROVED | exact Hurwitz reduction (`s06`) |
| $R$ lies in the $\mathbb Q$-span of alternating MZVs of weight $\le5$ | PROVED | Thm 1 makes $R$ a finite combination of linear Euler sums of level 2 (Sec. 2) |
| $R$ to 250 digits | PROVED (numerical, rigorous method) | Euler–Maclaurin on the exact law, two cut/order settings agree to $4\cdot10^{-259}$; matches referee's 40 digits |
| Closed form of $R$ (Eq. 2) | NUMERICAL IDENTIFICATION | PSLQ, basis of 20 alternating MZVs of weight $\le5$, 250 digits, residual $10^{-248}$, max coefficient 768 |
| Closed form of $b_1'''(2)$ (Eq. 3) | NUMERICAL IDENTIFICATION (same) | agrees with AV8's independent quadrature value to all 44 printed digits |
| A proof of Eq. 2 | OPEN | reduces to classical linear alternating Euler sums of weight $\le5$ (Sec. 2); not carried out |

## 1. The law for $b_n''(2)$ (Task 1)

**Second-order ODE (proved).** From $|S|^p+|qS'|^p=1$ with $q=\pi/\pi_p$ one gets, on $0<\theta<\pi/2$,
$S''=-q^{-p}\,|S/S'|^{p-2}\,S$. Write $p=2+e$, $S=\sin\theta+e\,s_1+\frac{e^2}2s_2+\dots$, $L_0=\log\tan\theta$,
and use $q^{-p}=1-e+(\frac{\pi^2}{16}+\frac14)e^2+O(e^3)$. Then
- order 1: $s_1''+s_1=(1-L_0)\sin\theta$;
- order 2: $s_2''+s_2=F_2:=-2L_0s_1+2\tan\theta\,s_1'-L_0^2\sin\theta+2L_0\sin\theta-(\frac12+\frac{\pi^2}8)\sin\theta$.

$s_1,s_2$ are $C^1$ with integrable second derivative (only $u^2\log u$-type terms at $\theta=0,\pi/2$), so
$(1-n^2)\,b_n^{(k)}=\widehat{F_k}(n)$ for $n\ge3$. The first-order equation alone reproves the $b_n'$ law in one line:
$L_0=-2\sum_{k\ \rm odd}\cos(2k\theta)/k$, so $\widehat{L_0\sin}(n)=-2/(n-1)$ ($n\equiv3$) or $2/(n+1)$ ($n\equiv1$).

**Fourier coefficients of $F_2$ (proved).** All four pieces are convolutions of known coefficients:
- $L_0s_1$: $2\sum_{k\ \rm odd}(\beta_{n-2k}+\beta_{n+2k})/k$ with $\beta_{-m}=-\beta_m$ (here $\beta_1=\frac{\log2-1}2$ enters once; this is the whole $\log 2$ part);
- $\tan\theta\,s_1'=:h$: from $\cos\theta\,h=\sin\theta\,s_1'$, $h_{2j+1}=c_j+2\sum_{r\ge1}(-1)^rc_{j+r}$, $c_i=(2i+1)\beta_{2i+1}$ (the unique $\ell^2$ solution; $h$ has only a $\log$ singularity);
- $L_0^2$: $\widehat{L_0^2}(2m)=8O_m/m$ for even $m$, $0$ for odd $m$, $O_m=\sum_{i\le m/2}\frac1{2i-1}$.

Each infinite sum is a rational function over an arithmetic progression. Partial fractions and
$\sum_{i\ge0}(i+x)^{-s}=\zeta(s,x)$ (and $-\psi$ for $s=1$) at $x\in t+\frac12\mathbb Z$ give harmonic numbers.
`symsum.py` does this exactly; `s02_bn2_law.py` applies it in each class $n=4t+1$, $n=4t+3$.

**Theorem 1.** For odd $n\ge3$, with $\varepsilon,a,c$ as above, $h=H_{a/4}$, $g=H_{a/2}$ and
$\tau=\sum_{i>a/4}\frac1{(2i-1)^2}=\frac{\pi^2}8-\sum_{i\le a/4}\frac1{(2i-1)^2}$,
$$b_n''(2)=\beta_n\log2-\frac{2\varepsilon n}{a^2c}\,\tau-\frac{2\varepsilon}{a^3}\,h+\frac{2\varepsilon(3n^2+2\varepsilon n+3)}{a^3c^2}\,g+\rho_n,\qquad(1)$$
$\rho_n=-\dfrac{2n(3n-1)}{a^2c^3}$ for $n\equiv1$, $\rho_n=\dfrac{2(3n^2+3n-2)}{a^2c^3}$ for $n\equiv3\pmod4$.

Consequences: the $\log2$ part is $\beta_n$; the $\pi^2$ part is $C_n\pi^2$ with $C_n=-\frac{\varepsilon n}{4a^2c}$ (the AV8 law);
the rest $A_n=\rho_n-\frac{2\varepsilon}{a^3}h+\frac{2\varepsilon(3n^2+2\varepsilon n+3)}{a^3c^2}g+\frac{2\varepsilon n}{a^2c}\sum_{i\le a/4}\frac1{(2i-1)^2}$ is rational.
Example $n=3$: $\frac3{32}+\frac{3\pi^2}{128}-\frac{\log2}8$, as in AV8. The "sporadic" denominators are those of the harmonic numbers.
Natural reading: $\pi^2$ never enters alone; it enters only through the odd-square tail $\tau=\frac14\psi'(\frac a4+\frac12)$.

## 2. The series $R$ (Task 2)

$R=\sum_{n\ge3}\beta_nA_n$ (referee's definition: $b_n''=A_n+\beta_n\log2+C_n\pi^2$). With $n=4k\pm1$, $k\ge1$, Thm 1 gives
$\beta_nb_n''$ as rational functions of $k$ (poles at $k=0,\pm\frac12$) times $1,H_k,H_{2k},\tau_k$. So $R$ is a finite sum of
level-2 linear Euler sums of weight $\le5$, hence an alternating MZV of weight $\le5$ (span of the 20 numbers below).

**Value (250 digits, rigorous method).** `s03b_R_fast.py`: partial sums to the cut, then Euler–Maclaurin with exact
Taylor data (polygamma and partial fractions). Settings $(T,K)=(300,85)$ and $(400,80)$ at 260 digits agree to $4\cdot10^{-259}$.
$R=-0.00922534340238471811655127394523099237323657779726753588695183892914\ldots$ (full value in `R_260_400.txt`).
This does not use any fit. It also returns $\sum\beta_nb_n''$, and the identity $\sum\beta_nb_n''=\log2\sum\beta_n^2+\pi^2\frac{1-\zeta(3)}{64}+R$ holds to $2\cdot10^{-269}$.

**Closed form (numerical identification).** PSLQ against $\{1,\log2,\pi^2,\log^22,\zeta(3),\pi^2\log2,\log^32,\pi^4,\zeta(3)\log2,\pi^2\log^22,\log^42,\mathrm{Li}_4(\frac12),\zeta(5),\pi^2\zeta(3),\pi^4\log2,\zeta(3)\log^22,\pi^2\log^32,\log^52,\mathrm{Li}_4(\frac12)\log2,\mathrm{Li}_5(\frac12)\}$
(20 elements, a spanning set for alternating MZVs of weight $\le5$), 250 digits, finds
$$R=\frac5{16}-\frac{13\pi^2}{96}+\frac{47}{64}\zeta(3)+\frac{\pi^2\log2}{32}-\frac{\pi^4}{384}+\frac{31}{256}\zeta(5)+\frac{\pi^2\zeta(3)}{256}.\qquad(2)$$
Residual $-9.9\cdot10^{-249}$; largest coefficient 768. Checked against both runs to $\le4\cdot10^{-259}$ (`s06`).
No $\mathrm{Li}_4(\frac12)$, $\mathrm{Li}_5(\frac12)$, Catalan or level-4 constants. (Catalan cannot occur: level 2 is proved.)
Why the referee's PSLQ failed: $\zeta(5)$ was not in its basis.

## 3. $b_1'''(2)$ (Task 3)

Parseval at third order: $b_1'''=\frac12N'''-3\beta_1b_1''-3\sum_{n\ge3}\beta_nb_n''$ and
$\sum_{n\ge3}\beta_nb_n''=\log2\sum\beta_n^2+\pi^2\frac{1-\zeta(3)}{64}+R$. With Eq. (2):
$$b_1'''(2)=-\frac38+\frac{33}{16}\log2-\frac{3\pi^2}{16}-\frac98\log^22+\frac{33}{64}\zeta(3)-\frac{\pi^2\log2}{64}+\frac{\log^32}8+\frac{3\pi^4}{512}-\frac{93}{256}\zeta(5)+\frac{9}{256}\pi^2\zeta(3)-\frac{\pi^4\log2}{512}\qquad(3)$$
$=-0.302622851041543176774778497461665409423075224893648486128162\ldots$
PSLQ on $b_1'''$ alone finds the same relation (residual $5\cdot10^{-249}$). AV8's independent quadrature value agrees to all 44 digits.

Odd-weight constants in $b_1'''(2)$, by source:

| Constant | from $\frac12N'''$ (Gamma) | from $-3\pi^2\sum\beta_nC_n$ | from $-3R$ | total |
|---|---|---|---|---|
| $\zeta(3)$ | $\frac{87}{32}$ | 0 | $-\frac{141}{64}$ | $\frac{33}{64}$ |
| $\pi^2\zeta(3)$ | 0 | $\frac3{64}=\frac{12}{256}$ | $-\frac3{256}$ | $\frac9{256}$ |
| $\zeta(5)$ | 0 | 0 | $-\frac{93}{256}$ | $-\frac{93}{256}$ |

The split between the $C_n$ part and $R$ is not natural: in Eq. (1) $\pi^2$ and the rational $\sum1/(2i-1)^2$ come together as the tail $\tau$.

**Verdict on the paper's sentence.** "a term $\frac3{64}\pi^2\zeta(3)$ from the products $b_n'b_n''$" is wrong as a statement about $b_1'''$.
The products give $\frac9{256}\pi^2\zeta(3)$ in total, and they also give $\zeta(5)$ ($-\frac{93}{256}$) and a $\zeta(3)$ term ($-\frac{141}{64}$) that changes the Gamma value $\frac{87}{32}$ to $\frac{33}{64}$.

## What the paper may say

Beyond Gamma, the products $b_n'b_n''$ add more odd zeta values. With the second-order law, which we prove for all $n$,
$b_1'''(2)=-\frac38+\frac{33}{16}\log2-\frac{3\pi^2}{16}-\frac98\log^22+\frac{33}{64}\zeta(3)-\frac{\pi^2\log2}{64}+\frac{\log^32}8+\frac{3\pi^4}{512}-\frac{93}{256}\zeta(5)+\frac9{256}\pi^2\zeta(3)-\frac{\pi^4\log2}{512}$
(identified numerically to 250 digits against all 20 constants of weight $\le5$). So the $\frac{87}{32}\zeta(3)$ from $N'''$ becomes $\frac{33}{64}\zeta(3)$, and $\pi^2\zeta(3)$ (coefficient $\frac9{256}$) and $\zeta(5)$ (coefficient $-\frac{93}{256}$) appear at the same order.
Table row: "$b_n''(2)$ law: Proved; $b_1'''(2)$ closed form: Numerical (250 digits)".
Also "$\zeta(3)$ waits for the third order" stays true, but "with $\zeta(3)$" should read "with $\zeta(3)$, $\pi^2\zeta(3)$ and $\zeta(5)$".

## Files (all in this folder)

- `symsum.py` — exact summation of rational functions over progressions into harmonic numbers.
- `s01_ode_check.py/.log` — first numeric test of the ODE route (its naive `nsum` is wrong in the 6th digit for $\tan\theta\,s_1'$ because the terms alternate between two classes; kept as a warning, superseded by s02).
- `s02_bn2_law.py/.log`, `b2_law.pkl` — Thm 1 per class, checked against AV8 data, $n=3..17$.
- `s02b_law_in_n.py/.log` — Thm 1 rewritten in $n$.
- `s03b_R_fast.py`, `s03b_50.log`, `s03b_260_a.log`, `s03b_260_b.log`, `R_*.txt` — values of $\sum\beta_nb_n''$, $R$, $b_1'''$.
- `s04_pslq.py`, `s04_pslq_250.log` — identification.
- `s05_quad_check.py/.log` — Thm 1 vs direct quadrature of $F_2$ for $n=19,21,31,33,61$.
- `s06_exact_parts.py/.log` — exact $\sum\beta_nC_n$, $\sum\beta_n^2$; closed forms vs 260-digit values.
